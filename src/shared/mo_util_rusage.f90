!> @file mo_util_rusage
!!
!! @brief Providing a query tool for the maximum resident size.
!!
!! @author Luis Kornblueh, Max Planck Institute for Meteorology
!!
module mo_util_rusage

  use, intrinsic :: iso_c_binding, only: c_int, c_long

  implicit none

  private

  type, bind(c) :: timeval
    integer(c_long) :: tv_sec                 ! seconds
#ifdef __linux__
    integer(c_long)  :: tv_usec                ! and microseconds
#else
    integer(c_int)  :: tv_usec                ! and microseconds
#endif
  end type timeval

  type, bind(c) :: rusage
    type(timeval)   :: ru_utime               ! user time used
    type(timeval)   :: ru_stime               ! system time used
    integer(c_long) :: ru_maxrss              ! max resident set size
    integer(c_long) :: ru_ixrss               ! integral shared text memory size
    integer(c_long) :: ru_idrss               ! integral unshared data size
    integer(c_long) :: ru_isrss               ! integral unshared stack size
    integer(c_long) :: ru_minflt              ! page reclaims
    integer(c_long) :: ru_majflt              ! page faults
    integer(c_long) :: ru_nswap               ! swaps
    integer(c_long) :: ru_inblock             ! block input operations
    integer(c_long) :: ru_oublock             ! block output operations
    integer(c_long) :: ru_msgsnd              ! messages sent
    integer(c_long) :: ru_msgrcv              ! messages received
    integer(c_long) :: ru_nsignals            ! signals received
    integer(c_long) :: ru_nvcsw               ! voluntary context switches
    integer(c_long) :: ru_nivcsw              ! involuntary context switches
  end type rusage

  integer(c_int), parameter :: RUSAGE_SELF     =  0
  integer(c_int), parameter :: RUSAGE_CHILDREN = -1

  interface
    function getrusage(who, r_usage) result(r) bind(c, name='getrusage')
      import :: c_int, c_long, rusage
      integer(c_int) :: r
      integer(c_int), value :: who
      type(rusage), intent(inout) :: r_usage
    end function getrusage
  end interface

  type rss
    type(rusage) :: used_rss
    integer :: idx
  end type rss

  type rss_list
    character(len=32) :: name
    integer :: idx
    integer :: used_rss
    type(rss), allocatable :: rss_usage(:)
    character(len=256) :: filename
    integer :: fileunit
  end type rss_list

  integer, save :: used_rss_lists       = 0
  type(rss_list), allocatable :: rss_lists(:)

  integer, parameter :: max_lists       = 32
  integer, parameter :: max_list_size   = 409600
  integer, parameter :: line_length     = 1024

  public :: add_rss_list
  public :: add_rss_usage
  public :: print_rss_usage
  public :: close_rss_lists

contains

  function add_rss_list(name,tag) result(idx)
    integer :: idx
    character(len=*), intent(in) :: name
    character(len=10), intent(inout),optional :: tag

    type(rss_list), allocatable :: tmp_rss_lists(:)
    integer :: ist

    if (.not. allocated(rss_lists)) then
      allocate(rss_lists(max_lists))
      used_rss_lists = 0
    endif

    if (used_rss_lists == size(rss_lists)) then
      allocate(tmp_rss_lists(2*used_rss_lists))
      tmp_rss_lists(1:size(rss_lists)) = rss_lists(:)
      call move_alloc(to=rss_lists,from=tmp_rss_lists)
    endif

    used_rss_lists = used_rss_lists + 1

    idx                             = used_rss_lists
    rss_lists(idx)%idx              = used_rss_lists
    rss_lists(idx)%name             = name
    allocate(rss_lists(idx)%rss_usage(max_list_size))
    rss_lists(idx)%rss_usage(:)%idx = 0
    rss_lists(idx)%used_rss         = 0

    if (present(tag)) then
      rss_lists(idx)%filename=trim(name)//'_'//trim(tag)//'.log'
      rss_lists(idx)%fileunit = find_next_free_unit(10,999)
      OPEN (UNIT=rss_lists(idx)%fileunit, FILE=rss_lists(idx)%filename, IOSTAT=ist, Recl=line_length)
      write(rss_lists(idx)%fileunit ,'(1x,a)')'idx        maxrss      majflt      minflt       nvcsw      nivcsw'
    else
      rss_lists(idx)%filename = ''
      rss_lists(idx)%fileunit = -1
    endif

  end function add_rss_list

  subroutine add_rss_usage(list_idx)
    integer, intent(in) :: list_idx

    type(rusage) :: ru
    integer :: max_size, idx, ret
    type(rss), allocatable :: tmp_rss_usage(:)
    character(len=line_length)     :: line

    if ((list_idx < 1) .or. (list_idx > used_rss_lists)) then
      print *, 'ERROR: list does not exist ...'
    endif

    max_size = size(rss_lists(list_idx)%rss_usage)
    idx = rss_lists(list_idx)%used_rss

    if (idx == max_list_size) then
      allocate(tmp_rss_usage(2*max_size))
      tmp_rss_usage(1:max_size) = rss_lists(list_idx)%rss_usage(:)
      call move_alloc(to=rss_lists(list_idx)%rss_usage,from=tmp_rss_usage)
    endif

    rss_lists(list_idx)%used_rss = rss_lists(list_idx)%used_rss + 1

    idx = rss_lists(list_idx)%used_rss
    rss_lists(list_idx)%rss_usage%idx  = idx
    ret = getrusage(RUSAGE_SELF, ru)
    rss_lists(list_idx)%rss_usage(idx)%used_rss = ru

    if (-1 .ne. rss_lists(list_idx)%fileunit) then
      line = ''
      write(line, '(1x,i5,5i12)') &
          &          idx, &
          &          rss_lists(list_idx)%rss_usage(idx)%used_rss%ru_maxrss, &
          &          rss_lists(list_idx)%rss_usage(idx)%used_rss%ru_majflt, &
          &          rss_lists(list_idx)%rss_usage(idx)%used_rss%ru_minflt, &
          &          rss_lists(list_idx)%rss_usage(idx)%used_rss%ru_nvcsw, &
          &          rss_lists(list_idx)%rss_usage(idx)%used_rss%ru_nivcsw
      !write(line,'(i6,a,i12)')idx,' ',rss_lists(list_idx)%rss_usage(idx)%used_rss%ru_maxrss
      write(rss_lists(list_idx)%fileunit,'(a)') trim(line)
      flush(rss_lists(list_idx)%fileunit)
    endif

  end subroutine add_rss_usage

  subroutine print_rss_usage()
    integer                        :: il, idx


    do il = 1, used_rss_lists
      print *, 'List: ', trim(rss_lists(il)%name), ' index: ', rss_lists(il)%idx
      print '(1x,a)', '            maxrss      majflt      minflt       nvcsw      nivcsw'
      do idx = 1, rss_lists(il)%used_rss
        print '(1x,i5,5i12)', &
          &          idx, &
          &          rss_lists(il)%rss_usage(idx)%used_rss%ru_maxrss, &
          &          rss_lists(il)%rss_usage(idx)%used_rss%ru_majflt, &
          &          rss_lists(il)%rss_usage(idx)%used_rss%ru_minflt, &
          &          rss_lists(il)%rss_usage(idx)%used_rss%ru_nvcsw, &
          &          rss_lists(il)%rss_usage(idx)%used_rss%ru_nivcsw
      enddo
    enddo

    print *, ''
#ifdef __linux__
    print *, ' maxrss - Maximum resident set size (kbytes)'
#else
    print *, ' maxrss - Maximum resident set size (bytes)'
#endif
    print *,' majflt - Major (requiring I/O) page faults'
    print *,' minflt - Minor (reclaiming a frame) page faults'
    print *,' nvcsw  - Voluntary context switches'
    print *,' nivcsw - Involuntary context switches'

  end subroutine print_rss_usage

  subroutine close_rss_lists(verbose)
    logical, optional :: verbose

    logical :: my_verbose
    integer :: il

    my_verbose = .false.
    if ( present(verbose) ) my_verbose = verbose

    do il = 1, used_rss_lists
      if (-1 .ne. rss_lists(il)%fileunit) then
        if (my_verbose) write(0,*) 'CLOSE: ', trim(rss_lists(il)%filename), ' index: ', rss_lists(il)%idx
        close(unit=rss_lists(il)%fileunit)
      endif
    enddo
  end subroutine close_rss_lists

  FUNCTION find_next_free_unit(istart,istop) RESULT(iunit)
    INTEGER :: iunit
    INTEGER, INTENT(in) :: istart, istop
    !
    INTEGER :: kstart, kstop
    LOGICAL :: lfound, lopened
    INTEGER :: i
    ! 
    lfound = .FALSE.
    !
    kstart = istart
    kstop  = istop
    IF (kstart < 10) kstart = 10
    IF (kstop <= kstart) kstop = kstart+10
    !
    DO i = kstart, kstop
      INQUIRE(unit=i, opened=lopened)
      IF (.NOT. lopened) THEN
        iunit = i
        lfound = .TRUE.
        EXIT
      END IF
    END DO
    !  
    IF (.NOT. lfound) THEN
      iunit = -1
    END IF
    !   
  END FUNCTION find_next_free_unit
end module mo_util_rusage
