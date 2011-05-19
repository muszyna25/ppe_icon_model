MODULE mo_real_timer

  ! utility for real time measurements
  ! on AIX a wrapper for the fast read_real_time function is used
  ! on SX a CPU counter is used
  ! on Linux/Intel a CPU register is read
  !
  ! Authors:
  !
  ! Initial version by:
  !
  ! J. Behrens, GWDG, March 2002, initial version for AIX
  !
  ! History:
  !
  ! L. Kornblueh, MPI, April 2002, extended for NEC SX and Linux/Intel
  ! S. Shingu. NEC, March 2004, bugfix
  ! Luis Kornblueh, MPI, May 2004, adapted for use in ICON
  ! J. Behrens, DKRZ, November 2009, Load balance diagnostics
  ! L. Kornblueh, MPI, March 2010, cleanup
  !

  USE, INTRINSIC :: ISO_C_BINDING

  USE mo_kind,        ONLY: dp
  USE mo_exception,   ONLY: finish, message, message_text
  USE mo_util_String, ONLY: separator

#ifdef _OPENMP
 USE omp_lib,       ONLY: omp_get_thread_num, omp_get_max_threads, &
                          omp_in_parallel, omp_get_num_threads
!                          omp_get_dynamic, omp_set_dynamic
#endif

#ifndef NOMPI
  USE mo_mpi,       ONLY: p_recv, p_send, p_barrier, p_real_dp, &
                          p_pe, p_io
#else
  USE mo_mpi,       ONLY: p_pe, p_io
#endif
  USE mo_parallel_nml, ONLY: num_test_procs, num_work_procs, &
                             p_test_run,                     &
                             p_comm_work, p_comm_work_test

  IMPLICIT NONE

  PRIVATE

#ifndef NOMPI
  INTEGER, PARAMETER :: report_tag = 12345
#endif

  PUBLIC :: t_time_mark_type, set_time_mark, get_time_val

  ! more informative, thread-safe statistics:

  PUBLIC :: new_timer, del_timer, timer_start, timer_stop
  PUBLIC :: timer_val, timer_last, timer_count, timer_average
  PUBLIC :: timer_reset, timer_reset_all
  PUBLIC :: timer_report

  INTEGER, PARAMETER :: timer_max = 256

  INTEGER :: top_timer = 0

  INTEGER, PARAMETER :: rt_undef_stat = 0
  INTEGER, PARAMETER :: rt_on_stat    = 1
  INTEGER, PARAMETER :: rt_off_stat   = 2

  ! minimal internal time needed to do one measuremenet

  REAL(dp) :: tm_shift = 0.0_dp

  ! average total overhead for one timer_start & timer_stop pair

  REAL(dp) :: tm_overhead = 0.0_dp

  ! shared part of timer

  TYPE t_srt_type
     LOGICAL           :: reserved             ! usage flag
     CHARACTER(len=80) :: text                 ! description of timer
  END TYPE t_srt_type

#ifdef __xlC__
  TYPE t_time_mark_type
    INTEGER :: t(4) ! 'raw' timer mark
  END TYPE t_time_mark_type
  TYPE(t_time_mark_type), PARAMETER :: init_mark = t_time_mark_type( (/ 0, 0, 0, 0 /) )
#else
  TYPE t_time_mark_type
    REAL(dp) :: t
  END TYPE t_time_mark_type
  TYPE(t_time_mark_type), PARAMETER :: init_mark = t_time_mark_type( 0.0_dp )
#endif

  ! thread private part of timer

  TYPE t_rt_type
     SEQUENCE
#ifdef __xlC__
     INTEGER  :: mark1(4)
#else
     REAL(dp) :: mark1               ! last start time
#endif
     REAL(dp) :: tot                 ! total sum of active time
     REAL(dp) :: min                 ! min. active time
     REAL(dp) :: max                 ! max. ..
     REAL(dp) :: last                ! last measurement
     INTEGER  :: call_n              ! number of calls
     INTEGER  :: stat                ! status
  END TYPE t_rt_type

  TYPE(t_srt_type), PARAMETER :: srt_init = t_srt_type(.FALSE., 'noname')
  TYPE(t_srt_type) :: srt(timer_max)

  TYPE(t_rt_type), PARAMETER :: rt_init = t_rt_type( &
#ifdef __xlC__
        (/ 0, 0, 0, 0/),  & ! mark1
#else
        0.0_dp, &           ! mark1
#endif
        0.0_dp, &           ! tot
        HUGE(0.0_dp), &     ! min
        0.0_dp, &           ! max
        0.0_dp, &           ! last
        0, &                ! call_n
        rt_undef_stat)      ! stat

  TYPE(t_rt_type) :: rt(timer_max)

#if defined (__SX__) || defined (ES)
! necessary due to a bug in NECs OMP implementation
LOCAL COMMON /th_real_timer/ rt
#else
#ifdef _OPENMP
  COMMON /th_real_timer/ rt
!$omp threadprivate(/th_real_timer/)
#endif
#endif

  LOGICAL :: need_init = .TRUE.

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
      INTEGER, INTENT(in) :: t(*)
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
      REAL(dp), INTENT(in) :: t
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

  SUBROUTINE init
    INTEGER :: sz

    INTEGER   :: ii = 0
    REAL (dp) :: dd = 0.0_dp
    INTEGER :: io_size, integer_byte_size, integer_io_size, realdp_byte_size

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

    need_init = .FALSE.

    CALL estimate_overhead

  END SUBROUTINE init

  !----

  SUBROUTINE estimate_overhead
    INTEGER, PARAMETER :: n = 100 ! tests need about n microsecs on pwr4
    REAL(dp) :: dt_min

    tm_shift = 0.0_dp
    CALL m1(dt_min)
    tm_shift = dt_min

    CALL m2(dt_min)
    tm_overhead = dt_min

  CONTAINS

    SUBROUTINE m1(dt0)
      REAL(dp), INTENT(out) :: dt0
      TYPE(t_time_mark_type) :: mark
      REAL(dp) :: dt
      INTEGER :: i

      dt0 = 1.0_dp
      DO i = 1, n
        CALL set_time_mark(mark)
        dt = get_time_val(mark)
        IF (dt < dt0) dt0 = dt
      ENDDO

    END SUBROUTINE m1

    SUBROUTINE m2(dt0)
      REAL(dp), INTENT(out) :: dt0
      TYPE(t_time_mark_type) :: mark1, mark2
      REAL(dp) :: dt1, dt2
      INTEGER :: i

      dt0 = 1.0_dp
      DO i = 1, n
        CALL set_time_mark(mark2)
        CALL set_time_mark(mark1)
        dt1 = get_time_val(mark1)
        dt2 = get_time_val(mark2)
        IF (dt2 < dt0) dt0 = dt2
        IF (dt2 < dt1) CALL real_timer_abort(reason='estimate_overhead:internal error')
      ENDDO

    END SUBROUTINE m2

  END SUBROUTINE estimate_overhead

  !----

  SUBROUTINE timer_reset_all

    INTEGER :: it

    DO it = 1, top_timer
!$omp parallel
      rt(it) = rt_init
!       rt(it)%tot    = 0.0_dp
!       rt(it)%min    = 0.0_dp
!       rt(it)%max    = 0.0_dp
!       rt(it)%last   = 0.0_dp
!       rt(it)%call_n = 0
!$omp end parallel
    ENDDO

  END SUBROUTINE timer_reset_all

  SUBROUTINE timer_reset(it)

    INTEGER, INTENT(in) :: it

    IF (it < 1 .OR. it > top_timer) &
         CALL real_timer_abort(it,'timer_reset: timer out of bounds')

!$omp parallel
    rt(it) = rt_init
!$omp end parallel
!     rt(it)%tot    = 0.0_dp
!     rt(it)%min    = 0.0_dp
!     rt(it)%max    = 0.0_dp
!     rt(it)%last   = 0.0_dp
!     rt(it)%call_n = 0

  END SUBROUTINE timer_reset


  SUBROUTINE timer_reset_field(it_field)
    INTEGER, INTENT(in) :: it_field(:)

    INTEGER :: iit, it

    DO iit = LBOUND(it_field,1), UBOUND(it_field,1)
      it = it_field(iit)
      IF (it < 1 .OR. it > top_timer) &
           CALL real_timer_abort(it,'timer_reset_field: timer out of bounds')

      rt(it)%tot    = 0.0_dp
      rt(it)%min    = 0.0_dp
      rt(it)%max    = 0.0_dp
      rt(it)%last   = 0.0_dp
      rt(it)%call_n = 0

    ENDDO

  END SUBROUTINE timer_reset_field

  !---

  INTEGER FUNCTION new_timer(text)
    CHARACTER(len=*), INTENT(in), OPTIONAL :: text

    INTEGER::jt

#ifdef _OPENMP
    IF ( omp_in_parallel() ) &
         CALL real_timer_abort(0,'new_timer called in parallel region')
#endif

    IF (need_init) CALL init
    top_timer = top_timer+1

    IF (top_timer > timer_max) THEN
       CALL message('new_timer','list of timers:')
       DO jt = 1, timer_max
          WRITE (message_text,'(i4,a)')  jt, srt(jt)%text
          CALL message ('new_timer',message_text)
       ENDDO
       CALL message ('new_timer','timer_max is too small')
       CALL real_timer_abort(jt,'new_timer failed')
    ENDIF
    srt(top_timer) = srt_init
    srt(top_timer)%reserved = .TRUE.
    IF (PRESENT(text)) srt(top_timer)%text = adjustl(text)

!$omp parallel
    rt(top_timer) = rt_init
!$omp end parallel

    new_timer = top_timer
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
    DO jt = top_timer, 1, -1
       IF (srt(jt)%reserved) EXIT
       top_timer = jt-1
    ENDDO

  END SUBROUTINE del_timer

  !---

  SUBROUTINE timer_start(it)
    INTEGER, INTENT(in) :: it

    IF (it < 1 .OR. it > top_timer) &
         CALL real_timer_abort(it,'timer_start: timer out of bounds')
    IF (rt(it)%stat == rt_on_stat) &
         CALL real_timer_abort(it,'timer_start: timer_stop call missing')

    call util_read_real_time(rt(it)%mark1)

    rt(it)%stat = rt_on_stat

  END SUBROUTINE timer_start

  !---

  SUBROUTINE set_time_mark(mark0)
    TYPE(t_time_mark_type), INTENT(out) :: mark0

    ! simple timer - no statistics

    mark0 = init_mark
    CALL util_read_real_time(mark0%t)

  END SUBROUTINE set_time_mark

  REAL(dp) FUNCTION get_time_val(mark0)
    TYPE(t_time_mark_type), INTENT(in) :: mark0

    ! simple timer - no statistics

    REAL(dp) :: dt

    TYPE(t_time_mark_type)::mark

    CALL util_read_real_time(mark%t)
    CALL util_diff_real_time(mark0%t,mark%t,dt)

    get_time_val = dt-tm_shift

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

    IF (it < 1 .OR. it > top_timer) &
         CALL real_timer_abort(it,'timer_val: invalid timer id')

    timer_val = rt(it)%tot
    IF (rt(it)%stat == rt_on_stat) THEN
      CALL util_read_real_time(mark2)
      CALL util_diff_real_time(rt(it)%mark1,mark2,dt)
      timer_val = timer_val+dt-tm_shift
    ENDIF

  END FUNCTION timer_val


  REAL(dp) FUNCTION timer_average(it)
    INTEGER, INTENT(in) :: it

    IF (it < 1 .OR. it > top_timer) &
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

    IF (it < 1 .OR. it > top_timer) &
         CALL real_timer_abort(it,'timer_last: invalid timer id')

    timer_last = rt(it)%last

  END FUNCTION timer_last

  !---

  INTEGER FUNCTION timer_count(it)
    INTEGER, INTENT(in) :: it

    IF (it < 1 .OR. it > top_timer) &
         CALL real_timer_abort(it,'timer_count: invalid timer id')

    timer_count = rt(it)%call_n

  END FUNCTION timer_count

  !---

  SUBROUTINE timer_stop(it)
    INTEGER, INTENT(in) :: it
#if defined (__xlC__)
    INTEGER :: mark2(4)
#else
    REAL(dp) :: mark2
#endif
    REAL(dp) :: dt

    IF (it < 1 .OR. it > top_timer) &
         CALL real_timer_abort(it,'timer_stop: invalid timer id')

    IF (rt(it)%stat /= rt_on_stat) THEN
      IF (rt(it)%stat == rt_off_stat) THEN
        CALL real_timer_abort(it,'timer_stop: timer_start call missing')
      ELSE
        CALL real_timer_abort(it,'timer_stop: undefined timer')
      ENDIF
    ENDIF

    CALL util_read_real_time(mark2)
    CALL util_diff_real_time(rt(it)%mark1,mark2,dt)

    dt = dt-tm_shift
    rt(it)%last = dt
    rt(it)%tot = rt(it)%tot + dt
    rt(it)%call_n = rt(it)%call_n+1
    IF (dt < rt(it)%min) rt(it)%min = dt
    IF (dt > rt(it)%max) rt(it)%max = dt
    rt(it)%stat = rt_off_stat

  END SUBROUTINE timer_stop

  !---

  SUBROUTINE timer_report(itimer, short)
    INTEGER, INTENT(in), OPTIONAL :: itimer      ! show this timer if present
    LOGICAL, INTENT(in), OPTIONAL :: short       ! generates condensed output if set

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

    CALL timer_report_full(it)

  END SUBROUTINE timer_report

  !---

  SUBROUTINE timer_report_short(itimer)
    INTEGER, INTENT(in) :: itimer

    INTEGER, PARAMETER :: i_sum = 1, i_min = 2, i_max = 3

    REAL(dp)::sbuf(3,top_timer), rbuf(3,top_timer,num_test_procs+num_work_procs), res(3,top_timer)
    REAL(dp) :: q, avg, alpha, e
#ifdef _OPENMP
    REAL(dp) :: t
#endif
    INTEGER  :: p_error, itpos(top_timer)
    INTEGER  :: ip, iit, it, it1, it2, n

    CHARACTER(len=12) :: min_str, avg_str, max_str, sum_str, e_str

    IF (itimer > 0) THEN
      it1 = itimer
      it2 = itimer
      CALL real_timer_abort(0,'timer_report_short: itimer>0 not supported (yet)')
    ELSE
      it1 = 1
      it2 = top_timer
    ENDIF

#ifdef _OPENMP

!$omp parallel private(t)
    DO it = 1, top_timer
!$omp critical
      t = rt(it)%tot
      sbuf(i_sum,it) = sbuf(i_sum,it)+t
      sbuf(i_min,it) = MIN(sbuf(i_min,it),t)
      sbuf(i_max,it) = MAX(sbuf(i_max,it),t)
!$omp end critical
    ENDDO
!$omp end parallel

    n = omp_get_num_threads()

#else
    n = 1
    DO it = 1, top_timer
      sbuf(:,it) = rt(it)%tot
    ENDDO

#endif

#ifndef NOMPI
    IF(p_test_run) THEN
      CALL MPI_GATHER(sbuf, SIZE(sbuf), p_real_dp, &
           rbuf, SIZE(sbuf), p_real_dp, &
           p_io, p_comm_work_test, p_error)
    ELSE
      CALL MPI_GATHER(sbuf, SIZE(sbuf), p_real_dp, &
           rbuf, SIZE(sbuf), p_real_dp, &
           p_io, p_comm_work, p_error)
    ENDIF
#else
    rbuf(:,:,1) = sbuf(:,:)
#endif


    n = n*(num_test_procs+num_work_procs)
    q = 1.0_dp/REAL(n,dp)
    res(:,:) = rbuf(:,:,1)
    DO ip = 2, num_test_procs+num_work_procs
      DO it = 1, top_timer
        res(i_sum,it) = res(i_sum,it)+rbuf(i_sum,it,ip)
        res(i_min,it) = MIN(res(i_min,it),rbuf(i_min,it,ip))
        res(i_max,it) = MAX(res(i_max,it),rbuf(i_max,it,ip))
      ENDDO
    ENDDO


    IF (p_pe == p_io) THEN

      CALL message ('',separator,all_print=.TRUE.)

      WRITE (message_text,'(A,I6,A)') ' Timer report ( tasknum * threadnum = ',n,')'
      CALL message ('',message_text,all_print=.TRUE.)

      WRITE (message_text,'(A,4A10,1X,A6)') ' label                       :  ', &
           't_min',   't_avg',   't_max',   't_sum', 'lbe[%]'

      CALL message ('',message_text,all_print=.TRUE.)
      CALL message ('',separator,all_print=.TRUE.)

      CALL mrgrnk(res(i_sum,:),itpos)

      DO iit = top_timer, 1, -1
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

  !---

  SUBROUTINE timer_report_full(itimer)
    INTEGER, INTENT(in) :: itimer

    INTEGER :: it
#ifndef NOMPI
    INTEGER :: ibuf(2)
#endif

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

    CALL message ('','',all_print=.TRUE.)
    CALL message ('',separator,all_print=.TRUE.)
    WRITE (message_text,'(a)') 'Timer report: '
    CALL message ('',message_text,all_print=.TRUE.)

    IF (num_test_procs+num_work_procs > 1) THEN
      WRITE (message_text,'(a,i5,a)') &
           'PE: ', p_pe, '                calls  t_min       t_average   t_max       t_total    '
      CALL message ('',message_text,all_print=.TRUE.)
    ELSE
      WRITE (message_text,'(a)') &
           '                         calls  t_min       t_average   t_max       t_total    '
      CALL message ('',message_text,all_print=.TRUE.)
    END IF
    CALL message ('',separator,all_print=.TRUE.)

    IF (itimer > 0) THEN
      IF (rt(itimer)%stat /= rt_undef_stat) CALL print_report(itimer)
    ELSE
      DO it = 1, top_timer
        IF (rt(it)%stat /= rt_undef_stat) THEN
          CALL print_report(it)
        ENDIF
      ENDDO
    ENDIF

#ifndef NOMPI
    IF (p_pe < num_test_procs+num_work_procs-1) THEN
      CALL p_send(ibuf(1), p_pe+1, report_tag)
    ENDIF
    IF(p_test_run) THEN
      CALL p_barrier(p_comm_work_test)
    ELSE
      CALL p_barrier(p_comm_work)
    ENDIF
#endif

    CALL message ('',separator,all_print=.TRUE.)

  END SUBROUTINE timer_report_full

  SUBROUTINE print_report(it1)
    INTEGER, INTENT(in) :: it1
    INTEGER :: tid
#ifdef _OPENMP
    INTEGER :: itid
#endif
    ! order omp:
!$OMP PARALLEL PRIVATE(itid,tid)
!$OMP DO ORDERED
#ifdef _OPENMP
    DO itid = 1, omp_get_num_threads()
      tid = omp_get_thread_num()
#else
      tid = 1
#endif

!$OMP ORDERED
      CALL report(it1)
!$OMP FLUSH
!$OMP END ORDERED
#ifdef _OPENMP
    ENDDO
#endif
!$OMP END DO
!$OMP END PARALLEL
  END SUBROUTINE print_report

  SUBROUTINE report(it)
    INTEGER, INTENT(in) :: it

    REAL(dp) :: avg, total
    CHARACTER(len=12) :: min_str, avg_str, max_str, tot_str
    INTEGER :: tid

#ifdef _OPENMP
    tid = omp_get_thread_num()
#else
    tid = -1
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
    IF (tid == -1) THEN
      WRITE (message_text,'(a32,i6,4a12,f10.3)') '    '//srt(it)%text, &
           rt(it)%call_n, min_str, avg_str, max_str, tot_str, total
      CALL message ('',message_text,all_print=.TRUE.)
    ELSE
      WRITE (message_text,'(i2,a32,i6,4a12,f10.3)') tid,': '//srt(it)%text, &
           rt(it)%call_n, min_str, avg_str, max_str, tot_str, total
      CALL message ('',message_text,all_print=.TRUE.)
    ENDIF

  END SUBROUTINE report

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

    rest = ts

    d = INT(rest/REAL(3600*24,dp))
    rest = rest-REAL(d*(3600*24),dp)
    IF (d > 99) THEN
      x = '>99d'
      time_str = ADJUSTR(x)
      RETURN
    ENDIF
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
      WRITE(x,'(f6.3,a)') ts, 's'
    ELSE
      WRITE(x,'(f6.5,a)') ts, 's'
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
      IF (it < 1 .OR. it > top_timer) THEN
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
