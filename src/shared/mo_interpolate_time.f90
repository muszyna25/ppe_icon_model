!----------------------------
#include "omp_definitions.inc"
!----------------------------
MODULE mo_interpolate_time
  USE mo_parallel_config,   ONLY: nproma
  USE mo_util_mtime,        ONLY: t_datetime_ptr, mtime_divide_timedelta
  USE mo_kind,              ONLY: wp
  USE mo_exception,         ONLY: message, message_text, finish
  USE mo_reader_abstract,   ONLY: t_abstract_reader
  USE mo_impl_constants,    ONLY: MAX_CHAR_LENGTH
  USE mtime,             ONLY: datetime, timedelta, &
    & OPERATOR(*), OPERATOR(+), OPERATOR(<), OPERATOR(>), &
    & datetimetostring, max_datetime_str_len, OPERATOR(-)
  USE mo_time_config,    ONLY: time_config
  USE mo_mpi,            ONLY: my_process_is_mpi_workroot, &
    & process_mpi_root_id, p_comm_work, p_bcast, p_pe_work
  
#ifdef _OPENACC
  USE mo_mpi,            ONLY: i_am_accel_node
#endif


  IMPLICIT NONE

  PUBLIC :: t_time_intp

  TYPE t_time_intp
    TYPE(t_datetime_ptr), ALLOCATABLE :: times(:)
    ! tidx and tidx+1 are held by dataold and datanew
    INTEGER                     :: tidx
    ! these two give access to the data
    REAL(wp),           POINTER :: dataold(:,:,:,:)
    REAL(wp),           POINTER :: datanew(:,:,:,:)
    ! these two hold the data and should not be accessed directly
    REAL(wp),       ALLOCATABLE :: dataa(:,:,:,:)
    REAL(wp),       ALLOCATABLE :: datab(:,:,:,:)

    ! Who am I? Passed to netcdf while reading.
    CHARACTER(len=MAX_CHAR_LENGTH)    :: var_name

    ! Options are constant, linear, and weird linear. Default is linear.
    INTEGER :: interpolation_mode

    CLASS(t_abstract_reader), POINTER :: reader
  CONTAINS
    procedure :: init => time_intp_init
    procedure :: intp => time_intp_intp
  END TYPE t_time_intp

  INTEGER, PARAMETER :: intModeConstant    = 0
  INTEGER, PARAMETER :: intModeLinear      = 1
  INTEGER, PARAMETER :: intModeLinearWeird = 11

  CHARACTER(len=*), PARAMETER :: modname = 'mo_interpolate_time'
CONTAINS

  SUBROUTINE time_intp_init(this, reader, local_time, var_name, int_mode)
    CLASS(t_time_intp),       TARGET, INTENT(  out) :: this
    CLASS(t_abstract_reader), TARGET, INTENT(inout) :: reader
    TYPE(datetime),          POINTER, INTENT(in   ) :: local_time
    CHARACTER(*),                     INTENT(in   ) :: var_name
    INTEGER,                OPTIONAL, INTENT(in   ) :: int_mode

    INTEGER :: ntimes
    INTEGER :: i

    CHARACTER(len=max_datetime_str_len)      :: date_str1, date_str2

    TYPE(timedelta) :: timeSinceDataStart


    CHARACTER(*), PARAMETER :: routine = &
      & modname//"::time_intp_init"

    this%reader   => reader
    this%var_name =  var_name

    IF (PRESENT(int_mode)) THEN
      this%interpolation_mode = int_mode
    ELSE
      this%interpolation_mode = intModeLinear
    ENDIF

    CALL this%reader%get_times(this%times)

    ntimes    = size(this%times)
    this%tidx = 1
    DO i = 1,ntimes
      IF (local_time > this%times(i)%ptr) THEN
        this%tidx = i
      ENDIF
    ENDDO

    timeSinceDataStart =  local_time - this%times(1)%ptr

    IF (time_config%tc_startdate < this%times(1)%ptr) THEN
      CALL datetimetostring(time_config%tc_startdate, date_str1)
      CALL datetimetostring(this%times(1)%ptr, date_str2)
      CALL finish(routine, "Start of this run ("//TRIM(date_str1)//") before start of data ("//TRIM(date_str2)//")")
    ENDIF
    IF (time_config%tc_stopdate > this%times(ntimes)%ptr) THEN
      CALL datetimetostring(time_config%tc_stopdate, date_str1)
      CALL datetimetostring(this%times(ntimes)%ptr, date_str2)
      CALL finish(routine, "End of this run ("//TRIM(date_str1)//") after end of data ("//TRIM(date_str2)//")")
    ENDIF

    ! log message
    CALL datetimetostring(this%times(this%tidx)%ptr, date_str1)
    CALL datetimetostring(this%times(this%tidx+1)%ptr, date_str2)
    WRITE(message_text,'(a,i3,a,i3,a)') " loading data, tidx", this%tidx,   " ("//TRIM(date_str1)//") and", &
      &                                                        this%tidx+1, " ("//TRIM(date_str2)//")"
    CALL message(TRIM(routine),message_text)


    !$ACC ENTER DATA COPYIN( this )
    CALL reader%get_one_timelev(this%tidx,   this%var_name, this%dataa)
    this%dataold => this%dataa
    CALL reader%get_one_timelev(this%tidx+1, this%var_name, this%datab)
    this%datanew => this%datab

  END SUBROUTINE time_intp_init

  SUBROUTINE time_intp_intp(this, local_time, interpolated)
    CLASS(t_time_intp), TARGET, INTENT(inout) :: this
    TYPE(datetime),    POINTER, INTENT(in   ) :: local_time
    REAL(wp),      ALLOCATABLE, INTENT(inout) :: interpolated(:,:,:,:)

    TYPE(timedelta) :: curr_delta, dt
    REAL(wp)                 :: weight
    INTEGER                  :: jc,jk,jb,jw
    INTEGER                  :: nlen, nblks, npromz, nlev
    CHARACTER(len=max_datetime_str_len) :: date_str

    CHARACTER(*), PARAMETER :: routine = &
      & modname//"::time_intp_intp"

    IF (local_time > this%times(this%tidx+1)%ptr) THEN
      this%tidx    = this%tidx + 1

      ! log message
      CALL datetimetostring(this%times(this%tidx+1)%ptr, date_str)
      WRITE(message_text,'(a,i3,a)') "loading new data, with tidx:", this%tidx+1, " ("//TRIM(date_str)//")"
      CALL message(TRIM(routine),message_text)

      ! FORTRAN!?!1! this should look like:
      ! DEALLOCATE(this%dataold)
      ! this%dataold => this%datanew
      ! CALL this%reader%get_one_timelev(this%tidx+1, this%var_name, this%datanew)
      ! but fails, since POINTERs are a strange hybrid of pointers and
      ! allocatables, so no chance to implicitly access the allocatable
      ! behind the pointer.
      IF(ASSOCIATED(this%dataold, this%dataa)) THEN
        this%dataold => this%datab
        CALL this%reader%get_one_timelev(this%tidx+1, this%var_name, this%dataa)
        this%datanew => this%dataa
      ELSE
        this%dataold => this%dataa
        CALL this%reader%get_one_timelev(this%tidx+1, this%var_name, this%datab)
        this%datanew => this%datab
      ENDIF
    ENDIF

    IF (this%interpolation_mode == intModeConstant) THEN
      ! A weight of 0 makes the interpolation return dataold.
      weight = 0.d0
    ELSE IF (this%interpolation_mode == intModeLinear) THEN
      curr_delta = local_time - this%times(this%tidx)%ptr
      dt         = this%times(this%tidx+1)%ptr - this%times(this%tidx)%ptr
      CALL mtime_divide_timedelta(curr_delta, dt, weight)
    ELSE IF (this%interpolation_mode == intModeLinearWeird) THEN
      ! This should
      ! a) be renamed. (But what is this?)
      ! b) done similar to calculate_time_interpolation_weights in
      !    shared/mo_bcs_time_interpolation.f90 . Since I do not know how this
      !    works and whether it can be generalized to data with non-monthly
      !    intervals, I leave as this for now.
      CALL finish(routine, "You are weird")
    ENDIF

    ! ATTENTION: This is a trivial Fortran 2008 feature. Does explode with Intel!
    !ALLOCATE(interpolated, MOLD=this%dataa)

    IF (ALLOCATED(interpolated)) THEN
      IF (.NOT. ALL( SHAPE(interpolated) .EQ. SHAPE(this%dataa) )) THEN
        !$ACC EXIT DATA DELETE( interpolated )
        DEALLOCATE(interpolated)
      END IF
    END IF
    IF (.NOT. ALLOCATED(interpolated)) THEN
      ALLOCATE(interpolated(size(this%dataa,1), size(this%dataa,2), size(this%dataa,3), size(this%dataa,4)))
      !$ACC ENTER DATA CREATE( interpolated )
    END IF

    !$ACC DATA PRESENT( interpolated )
    !$ACC KERNELS DEFAULT(NONE) IF (i_am_accel_node)
    interpolated(:,:,:,:) = 0.d0
    !$ACC END KERNELS

    nblks  = this%reader%get_nblks()
    npromz = this%reader%get_npromz()
    nlev   = size(interpolated,2)

    ! we need this mess, since npromz == nproma is not garantueed
    DO jw = 1,size(interpolated,4)
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,nlen,jc) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = 1,nblks
        nlen = MERGE(nproma, npromz, jb /= nblks)
        ! DA: Need to list this%dataxxx in the PRESENT section for attach
        !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) IF (i_am_accel_node) &
        !$ACC          DEFAULT(NONE) PRESENT( this, this%dataold, this%datanew )
        DO jk = 1,nlev
          DO jc = 1,nlen
            interpolated(jc,jk,jb,jw) = (1-weight) * this%dataold(jc,jk,jb,jw) &
              &                          + weight  * this%datanew(jc,jk,jb,jw)
          ENDDO
        ENDDO
      ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
    ENDDO

    !$ACC END DATA

!    if (my_process_is_mpi_workroot()) THEN
!      print *, "blubba weight", weight, interpolated(1,1,1,1), &
!      this%datanew(1,1,1,1), this%dataold(1,1,1,1)
!      print *, (1-weight) * this%dataold(1,1,1,1) + weight*this%datanew(1,1,1,1)
!    ENDIF
  END SUBROUTINE time_intp_intp

END MODULE mo_interpolate_time
