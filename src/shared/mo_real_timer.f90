!>
!! utility for real time measurements
!! on AIX a wrapper for the fast read_real_time function is used
!! on SX a CPU counter is used
!! on Linux/Intel a CPU register is read
!!
!! @par Revision History
!! Initial version for AIX by J. Behrens, GWDG, March 2002
!! extended for NEC SX and Linux/Intel by L. Kornblueh, MPI, April 2002
!! NEC SX OpenMP bugfix by S. Shingu, NEC, March 2004
!! adapted for use in ICON by Luis Kornblueh, MPI, May 2004
!! Load balance diagnostics by J. Behrens, DKRZ, November 2009
!! cleanup by L. Kornblueh, MPI, March 2010
!! Hierarchical statistics output and some cleanups by S. Koerner,
!!   DWD (2012-11-27)
!!
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
MODULE mo_real_timer

  USE mo_kind,            ONLY: dp
  USE mo_exception,       ONLY: finish, message, message_text
  USE mo_util_string,     ONLY: separator, sort_and_compress_list, int2string, real2string
  USE mo_util_table,      ONLY: t_table, initialize_table, finalize_table,    &
    &                           set_table_entry, print_table
  USE mo_impl_constants,  ONLY: MAX_CHAR_LENGTH, TIMER_MODE_WRITE_FILES,      &
    &                           TIMER_MODE_DETAILED, TIMER_MODE_AGGREGATED

#ifdef _OPENMP
  USE omp_lib,            ONLY: omp_get_thread_num, omp_get_max_threads, &
                                omp_in_parallel, omp_get_num_threads
#endif

#ifndef NOMPI
  USE mo_mpi,             ONLY: p_recv, p_send, p_barrier, p_real_dp, &
                                p_pe, get_my_mpi_all_comm_size, &
                                p_comm_size, p_n_work
#endif
  USE mo_parallel_config, ONLY: p_test_run

  USE mo_mpi,             ONLY: num_test_procs, get_my_mpi_work_id, &
    &                           get_mpi_comm_world_ranks, p_pe, p_pe_work, &
    &                           p_min, p_max,  &
    &                           num_work_procs, p_sum, p_allgather, mpi_land, &
    &                           p_allreduce, my_process_is_stdio, p_io, &
    &                           p_comm_work, p_comm_work_test, p_gather
  USE mo_master_control,  ONLY: get_my_process_name
  USE mo_run_config,      ONLY: profiling_output


  IMPLICIT NONE

  PRIVATE

#ifndef NOMPI
  INTEGER, PARAMETER :: report_tag = 12345
#else
  INTEGER, PARAMETER :: p_n_work = 1
#endif

  ! raw timers

  PUBLIC :: t_time_mark, set_time_mark, get_time_val

  ! more informative, thread-safe statistics:

  PUBLIC :: new_timer, del_timer, timer_start, timer_stop
  PUBLIC :: timer_val, timer_last, timer_count, timer_average
  PUBLIC :: timer_reset, timer_reset_all
  PUBLIC :: timer_report

  INTEGER, PARAMETER :: timer_max = 512

  INTEGER            :: timer_top = 0

  INTEGER, PARAMETER :: rt_undef_stat = 0
  INTEGER, PARAMETER :: rt_on_stat    = 1
  INTEGER, PARAMETER :: rt_off_stat   = 2


  ! overhead estimations, inner and outer
  ! set in "overhead_estimate"
  REAL(dp) :: delta_i=0.0_dp, delta_o=0.0_dp


  ! elementary time marks

#ifdef __xlC__
  TYPE t_time_mark
    INTEGER :: t(4)                     ! 'raw' timer mark
  END TYPE t_time_mark
  TYPE(t_time_mark), PARAMETER :: init_mark = t_time_mark( (/ 0, 0, 0, 0 /) )
#else
  TYPE t_time_mark
    REAL(dp) :: t
  END TYPE t_time_mark
  TYPE(t_time_mark), PARAMETER :: init_mark = t_time_mark( 0.0_dp )
#endif


  ! thread-shared part of timer
  TYPE t_srt
    LOGICAL           :: reserved       ! usage flag
    CHARACTER(len=80) :: text           ! description of timer
  END TYPE t_srt

  ! thread-private part of timer
  TYPE, BIND(c) :: t_rt
#ifdef __xlC__
    INTEGER  :: mark1(4)
#else
    REAL(dp) :: mark1                   ! last start time
#endif
    REAL(dp) :: tot                     ! total sum of active time
    REAL(dp) :: min                     ! min. active time
    REAL(dp) :: max                     ! max. ..
    REAL(dp) :: last                    ! last measurement
    INTEGER  :: call_n                  ! number of calls
    INTEGER  :: stat                    ! activity status (undefined/on/off)
    INTEGER  :: active_under            ! ID of superordinate timer; -1: "undefined", 0: "none", >0: id
  END TYPE t_rt

  TYPE t_timer_reductions
    REAL(dp), ALLOCATABLE :: val_min(:), val_max(:), val_tot(:), &
         val_tot_min(:), val_tot_max(:), val_call_n(:)
    INTEGER, ALLOCATABLE :: tot_rank_min(:), tot_rank_max(:), &
         rank_min(:), rank_max(:)
    LOGICAL :: inconsistent_timers
  END TYPE t_timer_reductions

  TYPE(t_srt), SAVE      :: srt(timer_max)
  TYPE(t_rt)             :: rt(timer_max)

  ! stack variables for the simultaneously active timer IDs, growing upwards
  ! thread-private
  ! We maintain the stack thread-private, to avoid OMP critical sections which would
  ! lead to more overhead. The more complicated analysis previous to the printout
  ! is acceptable.
  INTEGER, PARAMETER     :: active_timers_max=32       ! max number of simultaneously active timers
  INTEGER                :: active_timers(active_timers_max), &   ! initialization in "init"
                            active_timers_top
  COMMON /th_real_timer/ rt, active_timers, active_timers_top
  !$omp threadprivate(/th_real_timer/)

  TYPE(t_srt), PARAMETER :: srt_init = t_srt(.FALSE., 'noname')
  TYPE(t_rt) , PARAMETER :: rt_init  = t_rt( &
#ifdef __xlC__
        (/ 0, 0, 0, 0/),  &     ! mark1
#else
        0.0_dp, &               ! mark1
#endif
        0.0_dp, &               ! tot
        HUGE(0.0_dp), &         ! min
        0.0_dp, &               ! max
        0.0_dp, &               ! last
        0, &                    ! call_n
        rt_undef_stat, &        ! stat
        -1)                     ! active_under, value of -1 means "undefined"


  LOGICAL                :: need_init = .TRUE.


  INTERFACE
    SUBROUTINE util_init_real_time()
      IMPLICIT NONE
    END SUBROUTINE util_init_real_time

    SUBROUTINE util_get_real_time_size(sz)
      IMPLICIT NONE
      INTEGER, INTENT(out) :: sz
    END SUBROUTINE util_get_real_time_size

#if defined (__xlC__)
    SUBROUTINE util_read_real_time(t)
      IMPLICIT NONE
      INTEGER, INTENT(out) :: t(*)
    END SUBROUTINE util_read_real_time

    SUBROUTINE util_diff_real_time(t1,t2,dt)
      USE mo_kind, ONLY: dp
      IMPLICIT NONE
      INTEGER, INTENT(in)  :: t1(*), t2(*)
      REAL(dp), INTENT(out) :: dt
    END SUBROUTINE util_diff_real_time
#else
    SUBROUTINE util_read_real_time(t)
      USE mo_kind, ONLY: dp
      IMPLICIT NONE
      REAL(dp), INTENT(out) :: t
    END SUBROUTINE util_read_real_time

    SUBROUTINE util_diff_real_time(t1,t2,dt)
      USE mo_kind, ONLY: dp
      IMPLICIT NONE
      REAL(dp), INTENT(in)  :: t1, t2
      REAL(dp), INTENT(out) :: dt
    END SUBROUTINE util_diff_real_time
#endif

  END INTERFACE

CONTAINS

  SUBROUTINE mo_real_timer_init

    INTEGER   :: sz, ii = 0
    REAL (dp) :: dd = 0.0_dp
    INTEGER   :: io_size, integer_byte_size, integer_io_size, realdp_byte_size

  !------------------------------------------------------------------------------------------------
    CALL util_init_real_time()
    CALL util_get_real_time_size(sz)

#if defined (__xlC__)
    IF (BIT_SIZE(rt(1)%mark1)*SIZE(rt(1)%mark1) < sz*8) &
         CALL real_timer_abort(0,'buffer size for time stamps too small')
#else
    integer_byte_size = BIT_SIZE(ii)/8
    INQUIRE (iolength=io_size) ii
    integer_io_size = io_size
    INQUIRE (iolength=io_size) dd
    realdp_byte_size = io_size/integer_io_size*integer_byte_size
    IF (realdp_byte_size < sz) &
         CALL real_timer_abort(0,'buffer size for time stamps too small')
#endif

    ! initialize call stack variables
    !$omp parallel
    active_timers(1)  = 0
    active_timers_top = 1
    !$omp end parallel

    need_init = .FALSE.         ! must be placed before the call to overhead_estimate

    CALL overhead_estimate

  END SUBROUTINE mo_real_timer_init

  !----

  SUBROUTINE timer_reset_all

    INTEGER :: it

    !$omp parallel private(it)
    DO it = 1, timer_top
      rt(it) = rt_init
    ENDDO
    !$omp end parallel

  END SUBROUTINE timer_reset_all

  SUBROUTINE timer_reset(it)

    INTEGER, INTENT(in) :: it

    IF (it < 1 .OR. it > timer_top) &
         CALL real_timer_abort(it,'timer_reset: timer out of bounds')

    !$omp parallel
    rt(it) = rt_init
    !$omp end parallel

  END SUBROUTINE timer_reset


! outdated because of rt's subsequent data type extension
!  SUBROUTINE timer_reset_field(it_field)
!    INTEGER, INTENT(in) :: it_field(:)
!
!    INTEGER :: iit, it
!
!    DO iit = LBOUND(it_field,1), UBOUND(it_field,1)
!      it = it_field(iit)
!      IF (it < 1 .OR. it > timer_top) &
!           CALL real_timer_abort(it,'timer_reset_field: timer out of bounds')
!
!      rt(it)%tot    = 0.0_dp
!      rt(it)%min    = 0.0_dp
!      rt(it)%max    = 0.0_dp
!      rt(it)%last   = 0.0_dp
!      rt(it)%call_n = 0
!
!    ENDDO
!
!  END SUBROUTINE timer_reset_field

  !---

  RECURSIVE INTEGER FUNCTION new_timer(text, timer_id)

    CHARACTER(len=*), INTENT(in)   , OPTIONAL :: text
    INTEGER,          INTENT(inout), OPTIONAL :: timer_id

    INTEGER :: jt

  !------------------------------------------------------------------------------------------------
#ifdef _OPENMP
    IF ( omp_in_parallel() ) &
         CALL real_timer_abort(0,'new_timer called in parallel region')
#endif

    IF (need_init) CALL mo_real_timer_init

    IF (PRESENT(timer_id)) THEN
      new_timer = timer_id
      IF (timer_id > 0) RETURN
    ENDIF

    timer_top = timer_top+1
    IF (timer_top > timer_max) THEN
       CALL message('new_timer','list of timers:')
       DO jt = 1, timer_max
          WRITE (message_text,'(i4,a)')  jt, srt(jt)%text
          CALL message ('new_timer',message_text)
       ENDDO
       CALL message ('new_timer','timer_max is too small')
       CALL real_timer_abort(jt,'new_timer failed')
    ENDIF

    srt(timer_top)          = srt_init
    srt(timer_top)%reserved = .TRUE.
    IF (PRESENT(text)) srt(timer_top)%text = adjustl(text)

    !$omp parallel
    rt(timer_top) = rt_init
    !$omp end parallel

    new_timer = timer_top

    IF (PRESENT(timer_id)) THEN
      timer_id = new_timer
    ENDIF

  END FUNCTION new_timer

  !---

  SUBROUTINE del_timer(it)
    INTEGER, INTENT(in) :: it

    INTEGER :: jt

#ifdef _OPENMP
    IF ( omp_in_parallel() ) &
         CALL real_timer_abort(0,'del_timer called in parallel region')
#endif

    srt(it)%reserved = .FALSE.
    DO jt = timer_top, 1, -1
       IF (srt(jt)%reserved) EXIT
       timer_top = jt-1
    ENDDO

  END SUBROUTINE del_timer

  !---

  SUBROUTINE timer_start(it)

    INTEGER, INTENT(in) :: it

    INTEGER             :: it_sup               ! the actual superordinate timer

  !------------------------------------------------------------------------------------------------
    IF (it < 1 .OR. it > timer_top) &
         CALL real_timer_abort(it,'timer_start: timer out of bounds')
    IF (rt(it)%stat == rt_on_stat) &
         CALL real_timer_abort(it,'timer_start: timer_stop call missing')
    rt(it)%stat = rt_on_stat


    ! call-hierarchy bookkeeping: set <active_under>
    ! The actual superordinate timer is always active_timers(active_timers_top).
    it_sup = active_timers(active_timers_top)
    IF (rt(it)%active_under /= it_sup  .AND.  rt(it)%active_under /= 0) THEN
      IF (rt(it)%active_under == -1) THEN
        rt(it)%active_under = it_sup            ! first start of this timer at all
      ELSE
        rt(it)%active_under = 0                 ! this timer has been called by different
      ENDIF                                     !   superordinate timers, set to top for the printout
    ENDIF


    ! call hierarchy bookkeeping: set <active_timers>
    active_timers_top = active_timers_top + 1
    IF (active_timers_top > active_timers_max) &
         CALL real_timer_abort(it,'timer_start: number of simultaneously active timers higher than ''active_timers_max''')
    active_timers(active_timers_top) = it


    ! read real-time clock after all bookkeepings
    CALL util_read_real_time(rt(it)%mark1)

  END SUBROUTINE timer_start

  !---

  SUBROUTINE set_time_mark(mark0)
    TYPE(t_time_mark), INTENT(out) :: mark0

    ! simple timer - no statistics

    mark0 = init_mark
    CALL util_read_real_time(mark0%t)

  END SUBROUTINE set_time_mark

  REAL(dp) FUNCTION get_time_val(mark0)
    TYPE(t_time_mark), INTENT(in) :: mark0

    ! simple timer - no statistics

    REAL(dp) :: dt

    TYPE(t_time_mark)::mark

    CALL util_read_real_time(mark%t)
    CALL util_diff_real_time(mark0%t,mark%t,dt)

    get_time_val = dt-delta_i

  END FUNCTION get_time_val

  !---

  REAL(dp) FUNCTION timer_val(it)
    INTEGER, INTENT(in) :: it
#if defined (__xlC__)
    INTEGER :: mark2(4)
#else
    REAL(dp) :: mark2
#endif
    REAL(dp) :: dt

    IF (it < 1 .OR. it > timer_top) &
         CALL real_timer_abort(it,'timer_val: invalid timer id')

    timer_val = rt(it)%tot
    IF (rt(it)%stat == rt_on_stat) THEN
      CALL util_read_real_time(mark2)
      CALL util_diff_real_time(rt(it)%mark1,mark2,dt)
      timer_val = timer_val+dt-delta_i
    ENDIF

  END FUNCTION timer_val


  REAL(dp) FUNCTION timer_average(it)
    INTEGER, INTENT(in) :: it

    IF (it < 1 .OR. it > timer_top) &
         CALL real_timer_abort(it,'timer_average: invalid timer id')
    IF (rt(it)%stat == rt_on_stat) &
         CALL real_timer_abort(it,'timer_average: timer still active')
    IF (rt(it)%call_n == 0) &
         CALL real_timer_abort(it,'timer_average: timer never called')

    timer_average = rt(it)%tot/REAL(rt(it)%call_n,dp)

  END FUNCTION timer_average

  !---

  REAL(dp) FUNCTION timer_last(it)
    INTEGER, INTENT(in) :: it

    IF (it < 1 .OR. it > timer_top) &
         CALL real_timer_abort(it,'timer_last: invalid timer id')

    timer_last = rt(it)%last

  END FUNCTION timer_last

  !---

  INTEGER FUNCTION timer_count(it)
    INTEGER, INTENT(in) :: it

    IF (it < 1 .OR. it > timer_top) &
         CALL real_timer_abort(it,'timer_count: invalid timer id')

    timer_count = rt(it)%call_n

  END FUNCTION timer_count

  !---

  SUBROUTINE timer_stop(it)
    INTEGER, INTENT(in) :: it

#if defined (__xlC__)
    INTEGER  :: mark2(4)
#else
    REAL(dp) :: mark2
#endif
    REAL(dp) :: dt

  !------------------------------------------------------------------------------------------------
    ! do the time measurement first to minimize section c overhead
    CALL util_read_real_time(mark2)
    CALL util_diff_real_time(rt(it)%mark1,mark2,dt)

    IF (it < 1 .OR. it > timer_top) &
         CALL real_timer_abort(it,'timer_stop: invalid timer id')

    IF (rt(it)%stat /= rt_on_stat) THEN
      IF (rt(it)%stat == rt_off_stat) THEN
        CALL real_timer_abort(it,'timer_stop: timer_start call missing')
      ELSE
        CALL real_timer_abort(it,'timer_stop: undefined timer')
      ENDIF
    ENDIF


    dt            = dt-delta_i
    rt(it)%last   = dt
    rt(it)%tot    = rt(it)%tot + dt
    rt(it)%call_n = rt(it)%call_n+1
    IF (dt < rt(it)%min) rt(it)%min = dt
    IF (dt > rt(it)%max) rt(it)%max = dt
    rt(it)%stat   = rt_off_stat


    ! timer hierarchy bookkeeping
    ! timer intervals need to be nested
    IF (active_timers(active_timers_top) /= it) &
         CALL real_timer_abort(it,'timer_stop: a subsidary timer is still active, stop that first')
    active_timers_top = active_timers_top - 1


  END SUBROUTINE timer_stop

  !---

  SUBROUTINE timer_report(itimer, short)
    INTEGER, INTENT(in), OPTIONAL :: itimer      ! show this timer if present
    LOGICAL, INTENT(in), OPTIONAL :: short       ! generates condensed output if set
    ! local variables
    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = 'mo_real_timer:timer_report'
    INTEGER :: it

    IF (PRESENT(itimer)) THEN
      it = itimer
    ELSE
      it = -1 ! report of all timers
    ENDIF

    IF (PRESENT(short)) THEN
      IF (short) THEN
        CALL timer_report_short(it)
        RETURN
      ENDIF
    ENDIF

    SELECT CASE(profiling_output)
    CASE(TIMER_MODE_WRITE_FILES, TIMER_MODE_DETAILED)
      CALL timer_report_full(it)
    CASE(TIMER_MODE_AGGREGATED)
      IF (it > 0) THEN
        CALL timer_report_full(it)
      ELSE
        CALL timer_report_aggregated()
      END IF
    CASE DEFAULT
      CALL finish(routine, "Undefined profiling mode!")
    END SELECT

  END SUBROUTINE timer_report

  !---

  SUBROUTINE timer_report_short(itimer)
    INTEGER, INTENT(in) :: itimer

    INTEGER, PARAMETER :: i_sum = 1, i_min = 2, i_max = 3

    REAL(dp), ALLOCATABLE :: sbuf(:,:), sbuf_raw(:,:), rbuf(:,:,:), res(:,:)
    REAL(dp) :: q, avg, alpha, e, t
    INTEGER  :: itpos(timer_top), comm
    INTEGER  :: ip, iit, it, it1, it2, n, thread_id

    CHARACTER(len=12) :: min_str, avg_str, max_str, sum_str, e_str

    IF (itimer > 0) THEN
      it1 = itimer
      it2 = itimer
      CALL real_timer_abort(0,'timer_report_short: itimer>0 not supported (yet)')
    ELSE
      it1 = 1
      it2 = timer_top
    ENDIF

!$omp parallel private(t)
    n = 1
    thread_id = 0
!$  n = omp_get_num_threads()
!$  thread_id = omp_get_thread_num()
!$omp master
    ALLOCATE(sbuf_raw(timer_top,0:n-1), sbuf(timer_top,3), &
         rbuf(timer_top,3,num_test_procs+num_work_procs), res(timer_top,3))
!$omp end master
!$omp barrier
    DO it = 1, timer_top
      t = rt(it)%tot
      sbuf_raw(it,thread_id) = t
    ENDDO
!$omp barrier
!$omp do
    DO it = 1, timer_top
      sbuf(it,i_sum) = SUM(sbuf_raw(it,:))
      sbuf(it,i_min) = MINVAL(sbuf_raw(it,:))
      sbuf(it,i_max) = MAXVAL(sbuf_raw(it,:))
    ENDDO
!$omp end do nowait
!$omp end parallel

    comm = MERGE(p_comm_work, p_comm_work_test, .NOT. p_test_run)
    CALL p_gather(sbuf, rbuf, p_io, comm=comm)

    n = n*(num_test_procs+num_work_procs)
    q = 1.0_dp/REAL(n,dp)
    res(:,:) = rbuf(:,:,1)
    DO ip = 2, num_test_procs+num_work_procs
      DO it = 1, timer_top
        res(i_sum,it) = res(i_sum,it)+rbuf(i_sum,it,ip)
        res(i_min,it) = MIN(res(i_min,it),rbuf(i_min,it,ip))
        res(i_max,it) = MAX(res(i_max,it),rbuf(i_max,it,ip))
      ENDDO
    ENDDO


    IF (my_process_is_stdio()) THEN

      CALL message ('',separator,all_print=.TRUE.)

      WRITE (message_text,'(A,I6,A)') ' Timer report ( tasknum * threadnum = ',n,')'
      CALL message ('',message_text,all_print=.TRUE.)

      WRITE (message_text,'(A,4A10,1X,A6)') ' label                       :  ', &
           't_min',   't_avg',   't_max',   't_sum', 'lbe[%]'

      CALL message ('',message_text,all_print=.TRUE.)
      CALL message ('',separator,all_print=.TRUE.)

      CALL mrgrnk(res(i_sum,:),itpos)

      DO iit = timer_top, 1, -1
        it = itpos(iit)
        IF (rt(it)%stat == rt_undef_stat) CYCLE
        IF (res(i_sum,it) <= 0.0_dp) CYCLE
        avg = res(i_sum,it)*q
        alpha = ABS(res(i_max,it)-avg)/avg; !abs() to avoid -zero irritation
        e = 1.0_dp/(1.0_dp+alpha)
        sum_str = time_sec_str(res(i_sum,it))
        min_str = time_sec_str(res(i_min,it))
        max_str = time_sec_str(res(i_max,it))
        avg_str = time_sec_str(avg)
        WRITE (e_str,'(f6.2)') 100.0_dp*e
        WRITE (message_text,'(a,4a10,1x,a6)') ' '//srt(it)%text(1:27)//' :  ', &
             min_str, avg_str, max_str, sum_str, e_str
        CALL message ('',message_text,all_print=.TRUE.)
      ENDDO
      CALL message ('',separator,all_print=.TRUE.)

    ENDIF

  END SUBROUTINE timer_report_short


  ! --


  SUBROUTINE timer_report_full(itimer)
    INTEGER, INTENT(in) :: itimer

    INTEGER :: it,status
    INTEGER :: timer_file_id
    LOGICAL :: unit_is_occupied

#ifndef NOMPI
    INTEGER :: ibuf(2)
#endif

  !------------------------------------------------------------------------------------------------
#ifdef _OPENMP
    IF ( omp_in_parallel() ) &
         CALL real_timer_abort(0,'timer_report called in parallel region')
#endif

#ifndef NOMPI
    ! order mpi:
    IF (p_pe > 0) THEN
      CALL p_recv(ibuf(1), p_pe-1, report_tag)
    ENDIF
#endif

    IF (profiling_output == TIMER_MODE_WRITE_FILES) THEN
      DO timer_file_id = 500, 5000
        INQUIRE (UNIT=timer_file_id, OPENED=unit_is_occupied)
        IF ( .NOT. unit_is_occupied ) EXIT
      ENDDO
      IF (unit_is_occupied) &
        CALL finish("timer_report_full", "Cannot find available file unit")

      IF (get_my_process_name() /= "" ) THEN
        WRITE(message_text,'(a,a,a,i4.4)') 'timer.', TRIM(get_my_process_name()), ".", &
          &  get_my_mpi_work_id()
  !       write(0,*) "get_my_process_name() /= ''"
      ELSE
        WRITE(message_text,'(a,i4.4)') 'timer.', get_my_mpi_work_id()
      ENDIF
  !     write(0,*) "timer filename=", TRIM(message_text)

      OPEN (UNIT=timer_file_id, FILE=TRIM(message_text), FORM='FORMATTED',IOSTAT=status)
    ENDIF

    !
    !-- start the table output
    !

    CALL message ('','',all_print=.TRUE.)
    CALL message ('',separator,all_print=.TRUE.)

    IF (num_test_procs+num_work_procs > 1) THEN
      WRITE (message_text,'(a,i0)') 'Timer report of PE ', p_pe
    ELSE
      WRITE (message_text,'(a)'   ) 'Timer report'
    ENDIF
    CALL message ('',message_text,all_print=.TRUE.)
    IF (profiling_output == TIMER_MODE_WRITE_FILES) THEN
      WRITE(timer_file_id, '(a)' ) TRIM(message_text)
    END IF

    ! the right-aligned column heads
    WRITE (message_text, &
        '(a,  t4,a  , t32,a    , t46,a  , t54,a      , t70,a  , t86,a)') &
        'th', 'name', '# calls', 't_min', 't_average', 't_max', 't_total'
    CALL message ('',message_text,all_print=.TRUE.)
    IF (profiling_output == TIMER_MODE_WRITE_FILES) THEN
      WRITE(timer_file_id, '(a)') TRIM(message_text)
    END IF

    CALL message ('',separator,all_print=.TRUE.)

    !
    !-- print the timer data lines
    !
    IF (itimer > 0) THEN
      IF (rt(itimer)%stat /= rt_undef_stat) CALL print_report(itimer, timer_file_id, 0)
    ELSE
      DO it = 1, timer_top
        IF (rt(it)%stat /= rt_undef_stat .AND. rt(it)%active_under == 0) THEN
          CALL print_report_hierarchical(it, timer_file_id, 0)! print top-level timers hierarchical
        ENDIF
      ENDDO
    ENDIF

    IF (profiling_output == TIMER_MODE_WRITE_FILES) CLOSE(timer_file_id)

#ifndef NOMPI
    IF (p_pe < get_my_mpi_all_comm_size()-1) THEN
      CALL p_send(ibuf(1), p_pe+1, report_tag)
    ENDIF
    CALL p_barrier(MERGE(p_comm_work, p_comm_work_test, .NOT. p_test_run))
#endif

    CALL message ('',separator,all_print=.TRUE.)

  END SUBROUTINE timer_report_full


  !
  ! print statistics for timer <it> and all of its sub-timers
  !
  RECURSIVE SUBROUTINE print_report_hierarchical(it, timer_file_id, nd)

    INTEGER, INTENT(IN)          :: it, timer_file_id, &
                                    nd                   !<  nesting depth

    INTEGER :: n, &                                      ! number of sub-timers
               subtimer_list(timer_top), &               ! valid entries: 1..n
               k

  !------------------------------------------------------------------------------------------------

    CALL print_report(it, timer_file_id, nd)

    ! How many sub-timers has <it>, and which?
    n = 0
    DO k=1,timer_top
      IF (rt(k)%active_under == it) THEN
        n = n+1
        subtimer_list(n) = k
      ENDIF
    END DO

    ! print subtimers
    DO k=1,n
      CALL print_report_hierarchical(subtimer_list(k), timer_file_id, nd+1)
    ENDDO

  END SUBROUTINE print_report_hierarchical


  !
  !-- print statistics for one timer "it1", all OMP thread-copies
  !
  SUBROUTINE print_report(it1, timer_file_id, nd)
    INTEGER, INTENT(in) :: it1, timer_file_id, &
                           nd              ! nesting depth (determines the print indention)
    INTEGER :: tid
#if defined(_OPENMP)
    INTEGER :: itid
#endif
    ! order omp:
#if defined(_OPENMP)
!$OMP PARALLEL PRIVATE(itid,tid)
!$OMP DO ORDERED
    DO itid = 1, omp_get_num_threads()
      tid = omp_get_thread_num()
!$OMP ORDERED
#else
      tid = 1
#endif

      CALL print_reportline(it1, timer_file_id, nd)
#if defined(_OPENMP)
!$OMP FLUSH
!$OMP END ORDERED
    ENDDO
!$OMP END DO
!$OMP END PARALLEL
#endif
  END SUBROUTINE print_report


  SUBROUTINE print_reportline(it, timer_file_id, nd)
    INTEGER, INTENT(in) :: it, &                ! timer ID
                           timer_file_id, &
                           nd                   ! nesting depth (determines the print indention)

    REAL(dp)            :: avg, total
    CHARACTER(len=12)   :: min_str, avg_str, max_str, tot_str
    INTEGER             :: tid                  ! actual OMP thread ID

  !------------------------------------------------------------------------------------------------
#if defined(_OPENMP)
    tid = omp_get_thread_num()
#else
    tid = 0
#endif

    total = timer_val(it)

    avg = rt(it)%tot/REAL(MAX(1,rt(it)%call_n),dp)

    IF ( rt(it)%call_n > 0 ) THEN
      min_str = time_str(rt(it)%min)
      avg_str = time_str(avg)
      max_str = time_str(rt(it)%max)
    ELSE
      RETURN
      !min_str=''
      !avg_str=''
      !max_str=''
    ENDIF
    tot_str = time_str(total)

    WRITE (message_text, '(i2.2,1x,a29,i6,4a12,f14.5)')  &
        tid, &
        REPEAT('   ',MAX(nd-1,0))//REPEAT(' L ',MIN(nd,1))//srt(it)%text, &
        rt(it)%call_n, min_str, avg_str, max_str, tot_str, total
    CALL message ('',message_text,all_print=.TRUE.)
    IF (profiling_output == TIMER_MODE_WRITE_FILES) THEN
      WRITE(timer_file_id, '(a)') TRIM(message_text)
    END IF

  END SUBROUTINE print_reportline


  !---


  SUBROUTINE timer_report_aggregated()
    ! local variables
    INTEGER                        :: it, irow, max_threads
    INTEGER                        :: nranks
    INTEGER, ALLOCATABLE :: ranks(:)
    CHARACTER(LEN=MAX_CHAR_LENGTH) :: ranklist_str
    TYPE (t_table)                 :: table
    TYPE(t_timer_reductions)        :: tmr

#ifdef _OPENMP
    IF ( omp_in_parallel() ) &
         CALL real_timer_abort(0,'timer_report called in parallel region')
#endif

    !-- start the table output

    IF (p_pe_work == 0) THEN
      CALL message ('','',all_print=.TRUE.)
      IF (p_n_work > 1) THEN
        !-- build list of global ranks within this work communicator
        CALL get_mpi_comm_world_ranks(p_comm_work, ranks, nranks)
        CALL sort_and_compress_list(ranks, ranklist_str)
        WRITE (message_text,'(a,a)') 'Timer report, ranks ', TRIM(ranklist_str(1:MAX_CHAR_LENGTH-23))
      ELSE
        WRITE (message_text,'(a)') 'Timer report'
      END IF
    
      CALL message ('',message_text,all_print=.TRUE.)

      max_threads = 1
#ifdef _OPENMP
!$    max_threads = omp_get_max_threads()
#endif
      IF (max_threads > 1) THEN
        WRITE (message_text,'(a)') '(master thread only)'
        CALL message ('',message_text,all_print=.TRUE.)
      END IF

      !-- initialize print-out table
      CALL initialize_table(table)
      irow = 1
    END IF

    CALL get_timer_reductions(tmr)

    !-- collect the timer data lines
    DO it = 1, timer_top
      IF (rt(it)%stat /= rt_undef_stat .AND. rt(it)%active_under == 0) THEN
        CALL print_report_hierch_agg(table, irow, it, 0,tmr)
      ENDIF
    ENDDO

    IF (p_pe_work == 0) THEN
      !-- print the table
      IF (tmr%inconsistent_timers) &
           CALL message('', 'timer hierarchy inconsistencies detected, &
           &timing report truncated')
      CALL print_table(table, opt_delimiter='   ', opt_hline=.TRUE.)
      CALL message ('','',all_print=.TRUE.)
      !-- clean up
      CALL finalize_table(table)
    END IF

  END SUBROUTINE timer_report_aggregated

  SUBROUTINE get_timer_reductions(tmr)
    TYPE(t_timer_reductions), INTENT(out) :: tmr

    ALLOCATE(tmr%val_min(timer_top), tmr%val_max(timer_top), &
         tmr%val_tot(timer_top), tmr%val_tot_min(timer_top), &
         tmr%val_tot_max(timer_top), tmr%val_call_n(timer_top), &
         tmr%tot_rank_min(timer_top), tmr%tot_rank_max(timer_top), &
         tmr%rank_min(timer_top), tmr%rank_max(timer_top))
!$omp parallel
!$omp master
    tmr%val_min    = rt(1:timer_top)%min
    tmr%val_max    = rt(1:timer_top)%max
    tmr%val_tot    = rt(1:timer_top)%tot
    tmr%val_call_n = REAL(MAX(1,rt(1:timer_top)%call_n), dp)
!$omp end master
!$omp barrier
!$omp sections
!$omp section
    tmr%val_tot_min = tmr%val_tot
!$omp section
    tmr%val_tot_max = tmr%val_tot
!$omp section
    tmr%rank_min       = p_pe_work
!$omp section
    tmr%rank_max       = p_pe_work
!$omp section
    tmr%tot_rank_min   = p_pe_work
!$omp section
    tmr%tot_rank_max   = p_pe_work
!$omp end sections
!$omp end parallel
    ! Note: The following operations are MPI-collective!
    tmr%val_min = p_min(tmr%val_min, proc_id=tmr%rank_min, &
         comm=p_comm_work, root=0)
    tmr%val_max = p_max(tmr%val_max, proc_id=tmr%rank_max, &
         comm=p_comm_work, root=0)
    tmr%val_tot_min = p_min(tmr%val_tot_min, proc_id=tmr%tot_rank_min, &
         comm=p_comm_work, root=0)
    tmr%val_tot_max = p_max(tmr%val_tot_max, proc_id=tmr%tot_rank_max, &
         comm=p_comm_work, root=0)
    tmr%val_tot = p_sum(tmr%val_tot, comm=p_comm_work, root=0)
    tmr%val_call_n = p_sum(tmr%val_call_n, comm=p_comm_work, root=0)

    tmr%inconsistent_timers = .FALSE.
  END SUBROUTINE get_timer_reductions

  !
  ! print statistics for timer <it> and all of its sub-timers
  !
  ! Note: MPI-aggregated implementation.
  !
  RECURSIVE SUBROUTINE print_report_hierch_agg(table, irow, it, nd,tmr)
    TYPE (t_table),   INTENT(INOUT)       :: table       !< print-out table
    INTEGER,          INTENT(INOUT)       :: irow        !< table row index
    INTEGER,          INTENT(IN)          :: it
    INTEGER,          INTENT(IN)          :: nd          !<  nesting depth
    TYPE(t_timer_reductions), INTENT(inout) :: tmr
    ! local variables
    INTEGER :: i, n, &                                   ! number of sub-timers
      &        subtimer_list(timer_top), &               ! valid entries: 1..n
      &        k, ntasks
    INTEGER, ALLOCATABLE :: tcounts(:), tmrlists(:, :)
    LOGICAL :: consistent_sub_timer_counts, consistent_timer_lists

    CALL print_report_aggregated(table, irow, it, nd, tmr)
    irow = irow + 1

    ! How many sub-timers has <it>, and which?
    n = 0
    DO k=1,timer_top
      IF (rt(k)%active_under == it) THEN
        n = n+1
        subtimer_list(n) = k
      ENDIF
    END DO

#ifndef NOMPI
    ntasks = p_comm_size(p_comm_work)
#else
    ntasks = 1
#endif
    ALLOCATE(tcounts(0:ntasks-1))
    CALL p_allgather(n, tcounts, comm=p_comm_work)
    consistent_sub_timer_counts = ALL(tcounts(1:) == tcounts(0))
    consistent_timer_lists = .FALSE.
    consistent_sub_timer_counts &
         = p_allreduce(consistent_sub_timer_counts, mpi_land, p_comm_work)
    IF (consistent_sub_timer_counts) THEN
      ALLOCATE(tmrlists(n, 0:ntasks-1))
      CALL p_allgather(subtimer_list(1:n), tmrlists, comm=p_comm_work)
      consistent_timer_lists = .TRUE.
      DO i = 1, n
        consistent_timer_lists = consistent_timer_lists .AND. &
             ALL(tmrlists(i, 0:ntasks-1) == tmrlists(i, 0))
      END DO
    END IF
    IF (consistent_sub_timer_counts .AND. consistent_timer_lists) THEN
      ! print subtimers
      DO k=1,n
        CALL print_report_hierch_agg(table, irow, subtimer_list(k), nd+1,tmr)
      ENDDO
    ELSE
      IF (my_process_is_stdio()) THEN
        WRITE (0, '(2a,3(", ",i0))') 'problem: sub-timers inconsistent! ', TRIM(srt(it)%text), it, irow - 1, n
        IF(consistent_sub_timer_counts) THEN
            WRITE(0,*) "cause: inconsistent timer lists"
        ELSE
            WRITE(0,*) "cause: inconsistent subtimer counts: ", tcounts
        END IF
      END IF
      tmr%inconsistent_timers = .TRUE.
    END IF
  END SUBROUTINE print_report_hierch_agg


  !
  !-- print statistics for one timer "it1", merges all ranks in MPI
  !-- work communicator and all OMP thread-copies.
  !
  SUBROUTINE print_report_aggregated(table, irow, it, nd, tmr)
    TYPE (t_table),   INTENT(INOUT)       :: table           !< print-out table
    INTEGER,          INTENT(INOUT)       :: irow            !< table row index
    INTEGER,          INTENT(in)          :: it
    INTEGER,          INTENT(in)          :: nd              !< nesting depth (determines the print indention)
    TYPE(t_timer_reductions), INTENT(in)  :: tmr
    ! local variables
    REAL(dp)            :: val_avg
    CHARACTER(len=12)   :: min_str, avg_str, max_str, tot_min_str, tot_max_str

    ! OpenMP: We take the values of the master thread


    IF ((p_pe_work == 0) .AND. ( tmr%val_call_n(it) > 1.e-6_dp )) THEN
      val_avg = tmr%val_tot(it)/tmr%val_call_n(it)
      min_str = time_str(tmr%val_min(it))
      avg_str = time_str(val_avg)
      max_str = time_str(tmr%val_max(it))

      WRITE (tot_min_str, '(f12.3)') tmr%val_tot_min(it)
      WRITE (tot_max_str, '(f12.3)') tmr%val_tot_max(it)

      CALL set_table_entry(table, irow, "name",     &
        &                  REPEAT('   ',MAX(nd-1,0))//REPEAT(' L ',MIN(nd,1))//srt(it)%text)

      ! avoid integer overflow; very large call counts are useless
      ! anyway ...
      IF (tmr%val_call_n(it) <= REAL(HUGE(INT(1)),dp)) THEN
        CALL set_table_entry(table, irow, "# calls",  TRIM(int2string(INT(tmr%val_call_n(it)))))
      ELSE
        CALL set_table_entry(table, irow, "# calls",  TRIM(real2string(tmr%val_call_n(it),'(ES10.1)')))
      END IF
      CALL set_table_entry(table, irow, "t_min",         TRIM(min_str))
      IF (p_n_work > 1) &
        CALL set_table_entry(table, irow, "min rank", "["//TRIM(int2string(tmr%rank_min(it)))//"]")
      CALL set_table_entry(table, irow, "t_avg",         TRIM(avg_str))
      CALL set_table_entry(table, irow, "t_max",         TRIM(max_str))
      IF (p_n_work > 1) &
        CALL set_table_entry(table, irow, "max rank", "["//TRIM(int2string(tmr%rank_max(it)))//"]")
      CALL set_table_entry(table, irow, "total min (s)",           TRIM(tot_min_str))
      IF (p_n_work > 1) &
        CALL set_table_entry(table, irow, "total min rank", "["//TRIM(int2string(tmr%tot_rank_min(it)))//"]")
      CALL set_table_entry(table, irow, "total max (s)",           TRIM(tot_max_str))
      IF (p_n_work > 1) &
        CALL set_table_entry(table, irow, "total max rank", "["//TRIM(int2string(tmr%tot_rank_max(it)))//"]")
    ENDIF

  END SUBROUTINE print_report_aggregated


  !
  !--- set the two overhead times delta_i and delta_o (SAVEd module variables)
  !--- see notes in H8 around 12 November 2012 (SK)
  !
  SUBROUTINE overhead_estimate

    INTEGER, PARAMETER  :: n = 100
    INTEGER             :: k, tt1, tt2
    REAL(dp)            :: dt_min

  !------------------------------------------------------------------------------------------------
    ! allocate two test timers
    tt1 = new_timer()
    tt2 = new_timer()


    ! estimate delta_i
    delta_i = 0.0_dp
    dt_min  = HUGE(dt_min)
    DO k=1,n
      CALL timer_start(tt1)
      CALL timer_stop(tt1)
      dt_min = MIN(dt_min,timer_last(tt1))
    ENDDO
    delta_i = dt_min


    ! estimate delta_o
    dt_min = HUGE(dt_min)
    DO k=1,n
      CALL timer_start(tt1)
        CALL timer_start(tt2)
        CALL timer_stop(tt2)
      CALL timer_stop(tt1)
      dt_min = MIN(dt_min,timer_last(tt1))
      ! the following line is to prevent code elimination of the inner timer calls
      IF (timer_last(tt2) < -5.0_dp) CALL real_timer_abort(reason='estimate_overhead:internal error')
    ENDDO
    delta_o = dt_min


    ! free the two test timers
    CALL del_timer(tt2)
    CALL del_timer(tt1)

  END SUBROUTINE overhead_estimate


  !---

  CHARACTER(len=10) FUNCTION time_sec_str(ts)
    REAL(dp), INTENT(in) :: ts

    IF (ts < 0.0_dp) THEN
      time_sec_str="    ??????"
    ELSE IF(ts < 1.e1_dp) THEN
      WRITE(time_sec_str,'(f10.4)') ts
    ELSE IF(ts < 1.e2_dp) THEN
      WRITE(time_sec_str,'(f10.3)') ts
    ELSE IF(ts < 1.e3_dp) THEN
      WRITE(time_sec_str,'(f10.2)') ts
    ELSE IF(ts < 1.e4_dp) THEN
      WRITE(time_sec_str,'(f10.1)') ts
    ELSE
      WRITE(time_sec_str,'(f10.0)') ts
    ENDIF

  END FUNCTION time_sec_str

  !---

  CHARACTER(len=12) FUNCTION time_str(ts)
    REAL(dp), INTENT(in) :: ts

    REAL(dp) :: rest

    INTEGER :: d, h, m, s

    CHARACTER(len=2) :: d_str, h_str, m_str, s_str
    CHARACTER(len=12) :: x

    ! ensure that we don't crash with illegal inputs
    IF(ts > REAL(99*3600*24,dp)) THEN
        time_str = '>99d'
        RETURN
    ENDIF
    IF(ts < 0) THEN
        time_str = 'negative'
        RETURN
    END IF

    rest = ts

    d = INT(rest/REAL(3600*24,dp))
    rest = rest-REAL(d*(3600*24),dp)
    WRITE(d_str,'(i2.2)') d

    h = INT(rest/3600.0_dp)
    rest = rest-REAL(h*3600,dp)
    WRITE(h_str,'(i2.2)') h

    m = INT(rest/60.0_dp)
    rest = rest-REAL(m*60,dp)
    WRITE(m_str,'(i2.2)') m

    s = INT(rest)
    rest = rest-REAL(s,dp)
    WRITE(s_str,'(i2.2)') s

    IF (d > 0) THEN
      x = TRIM(d_str)//'d'//TRIM(h_str)//'h'
    ELSEIF (h > 0) THEN
      x = TRIM(h_str)//'h'//TRIM(m_str)//'m'
    ELSEIF (m > 0) THEN
      x = TRIM(m_str)//'m'//TRIM(s_str)//'s'
    ELSEIF (ts >= 1.0_dp) THEN
      WRITE(x,'(f7.4,a)') ts, 's'
    ELSE
      WRITE(x,'(f8.5,a)') ts, 's'
    ENDIF
    time_str = ADJUSTR(x)

  END FUNCTION time_str

  !---

  SUBROUTINE real_timer_abort(it,reason)
    INTEGER,          INTENT(in), OPTIONAL :: it
    CHARACTER(len=*), INTENT(in), OPTIONAL :: reason

    WRITE (message_text,*)  'error in module mo_real_timer:'
    CALL message ('', TRIM(message_text))
    IF (PRESENT(it)) THEN
      WRITE (message_text,*) 'timer handle: ', it
      CALL message ('', TRIM(message_text))
      IF (it < 1 .OR. it > timer_top) THEN
        WRITE (message_text,*) 'timer name: unspecified'
      CALL message ('', TRIM(message_text))
    ELSE
        WRITE (message_text,*) 'timer name: ', TRIM(srt(it)%text)
      CALL message ('', TRIM(message_text))
      ENDIF
    ENDIF
    IF (PRESENT(reason)) THEN
      WRITE (message_text,*) '            ', reason
      CALL message ('', TRIM(message_text))
    ENDIF

    CALL finish('real_timer_abort','Abort')

  END SUBROUTINE real_timer_abort

  ! from Orderpack 2.0: ranking of an array
  ! __________________________________________________________
  ! Originally: MRGRNK = Merge-sort ranking of an array

  SUBROUTINE mrgrnk(xdont, irngt)
    REAL (dp), INTENT (in)  :: xdont(:)
    INTEGER,   INTENT (out) :: irngt(:)

    REAL (dp) :: xvala, xvalb

    INTEGER  :: jwrkt(SIZE(irngt))
    INTEGER :: lmtna, lmtnc, irng1, irng2
    INTEGER :: nval, iind, iwrkd, iwrk, iwrkf, jinda, iinda, iindb

    nval = MIN (SIZE(xdont), SIZE(irngt))
    SELECT CASE (nval)
    CASE (:0)
      RETURN
    CASE (1)
      irngt (1) = 1
      RETURN
    CASE default
      CONTINUE
    END SELECT
    !
    !  Fill-in the index array, creating ordered couples
    !
    DO iind = 2, nval, 2
      IF (xdont(iind-1) <= xdont(iind)) THEN
        irngt (iind-1) = iind - 1
        irngt (iind) = iind
      ELSE
        irngt (iind-1) = iind
        irngt (iind) = iind - 1
      END IF
    END DO
    IF (MODULO(nval, 2) /= 0) THEN
      irngt (nval) = nval
    END IF
    !
    !  We will now have ordered subsets A - B - A - B - ...
    !  and merge A and B couples into     C   -   C   - ...
    !
    lmtna = 2
    lmtnc = 4
    !
    !  First iteration. The length of the ordered subsets goes from 2 to 4
    !
    DO
      IF (nval <= 2) EXIT
      !
      !   Loop on merges of A and B into C
      !
      DO iwrkd = 0, nval - 1, 4
        IF ((iwrkd+4) > nval) THEN
          IF ((iwrkd+2) >= nval) EXIT
          !
          !   1 2 3
          !
          IF (xdont(irngt(iwrkd+2)) <= xdont(irngt(iwrkd+3))) EXIT
          !
          !   1 3 2
          !
          IF (xdont(irngt(iwrkd+1)) <= xdont(irngt(iwrkd+3))) THEN
            irng2 = irngt (iwrkd+2)
            irngt (iwrkd+2) = irngt (iwrkd+3)
            irngt (iwrkd+3) = irng2
            !
            !   3 1 2
            !
          ELSE
            irng1 = irngt (iwrkd+1)
            irngt (iwrkd+1) = irngt (iwrkd+3)
            irngt (iwrkd+3) = irngt (iwrkd+2)
            irngt (iwrkd+2) = irng1
          END IF
          EXIT
        END IF
        !
        !   1 2 3 4
        !
        IF (xdont(irngt(iwrkd+2)) <= xdont(irngt(iwrkd+3))) CYCLE
        !
        !   1 3 x x
        !
        IF (xdont(irngt(iwrkd+1)) <= xdont(irngt(iwrkd+3))) THEN
          irng2 = irngt (iwrkd+2)
          irngt (iwrkd+2) = irngt (iwrkd+3)
          IF (xdont(irng2) <= xdont(irngt(iwrkd+4))) THEN
            !   1 3 2 4
            irngt (iwrkd+3) = irng2
          ELSE
            !   1 3 4 2
            irngt (iwrkd+3) = irngt (iwrkd+4)
            irngt (iwrkd+4) = irng2
          END IF
          !
          !   3 x x x
          !
        ELSE
          irng1 = irngt (iwrkd+1)
          irng2 = irngt (iwrkd+2)
          irngt (iwrkd+1) = irngt (iwrkd+3)
          IF (xdont(irng1) <= xdont(irngt(iwrkd+4))) THEN
            irngt (iwrkd+2) = irng1
            IF (xdont(irng2) <= xdont(irngt(iwrkd+4))) THEN
              !   3 1 2 4
              irngt (iwrkd+3) = irng2
            ELSE
              !   3 1 4 2
              irngt (iwrkd+3) = irngt (iwrkd+4)
              irngt (iwrkd+4) = irng2
            END IF
          ELSE
            !   3 4 1 2
            irngt (iwrkd+2) = irngt (iwrkd+4)
            irngt (iwrkd+3) = irng1
            irngt (iwrkd+4) = irng2
          END IF
        END IF
      END DO
      !
      !  The Cs become As and Bs
      !
      lmtna = 4
      EXIT
    END DO
    !
    !  Iteration loop. Each time, the length of the ordered subsets
    !  is doubled.
    !
    DO
      IF (lmtna >= nval) EXIT
      iwrkf = 0
      lmtnc = 2 * lmtnc
      !
      !   Loop on merges of A and B into C
      !
      DO
        iwrk = iwrkf
        iwrkd = iwrkf + 1
        jinda = iwrkf + lmtna
        iwrkf = iwrkf + lmtnc
        IF (iwrkf >= nval) THEN
          IF (jinda >= nval) EXIT
          iwrkf = nval
        END IF
        iinda = 1
        iindb = jinda + 1
        !
        !   Shortcut for the case when the max of A is smaller
        !   than the min of B. This line may be activated when the
        !   initial set is already close to sorted.
        !
        !          if (xdont(irngt(jinda)) <= xdont(irngt(iindb))) cycle
        !
        !  One steps in the C subset, that we build in the final rank array
        !
        !  Make a copy of the rank array for the merge iteration
        !
        jwrkt (1:lmtna) = irngt (iwrkd:jinda)
        !
        xvala = xdont (jwrkt(iinda))
        xvalb = xdont (irngt(iindb))
        !
        DO
          iwrk = iwrk + 1
          !
          !  We still have unprocessed values in both A and B
          !
          IF (xvala > xvalb) THEN
            irngt (iwrk) = irngt (iindb)
            iindb = iindb + 1
            IF (iindb > iwrkf) THEN
              !  Only A still with unprocessed values
              irngt (iwrk+1:iwrkf) = jwrkt (iinda:lmtna)
              EXIT
            END IF
            xvalb = xdont (irngt(iindb))
          ELSE
            irngt (iwrk) = jwrkt (iinda)
            iinda = iinda + 1
            IF (iinda > lmtna) EXIT  ! Only B still with unprocessed values
            xvala = xdont (jwrkt(iinda))
          END IF
          !
        END DO
      END DO
      !
      !  The Cs become As and Bs
      !
      lmtna = 2 * lmtna
    END DO
    !
  END SUBROUTINE mrgrnk

END MODULE mo_real_timer
