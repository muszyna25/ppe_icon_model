!-------------------------------------------------------------------------------------
!>
!! Set of methods for simple statistics
!! NOTE: in order to get correct results make sure you provide the proper range_subset!
!!
!! @author Leonidas Linardakis, MPI-M
!!
!! @par Revision History
!!   First implementation by Leonidas Linardakis, MPI-M, 2012-01-19
!!
!! @par Copyright
!! 2002-2010 by DWD and MPI-M
!! This software is provided for non-commercial use only.
!! See the LICENSE and the WARRANTY conditions.
!!
!! @par License
!! The use of ICON is hereby granted free of charge for an unlimited time,
!! provided the following rules are accepted and applied:
!! <ol>
!! <li> You may use or modify this code for your own non commercial and non
!!    violent purposes.
!! <li> The code may not be re-distributed without the consent of the authors.
!! <li> The copyright notice and statement of authorship must appear in all
!!    copies.
!! <li> You accept the warranty conditions (see WARRANTY).
!! <li> In case you intend to use the code commercially, we oblige you to sign
!!    an according license agreement with DWD and MPI-M.
!! </ol>
!!
!! @par Warranty
!! This code has been tested up to a certain level. Defects and weaknesses,
!! which may be included in the code, do not establish any warranties by the
!! authors.
!! The authors do not make any warranty, express or implied, or assume any
!! liability or responsibility for the use, acquisition or application of this
!! software.
!-------------------------------------------------------------------------------------
MODULE mo_statistics
  !-------------------------------------------------------------------------
  USE mo_kind,               ONLY: wp
  USE mo_exception,          ONLY: message_text, message, warning, finish

#ifndef __ICON_GRID_GENERATOR__
  USE mo_grid_subset,        ONLY: t_subset_range, get_index_range, t_subset_indexed
  USE mo_mpi,                ONLY: process_mpi_stdio_id, get_my_mpi_work_communicator, p_max, p_min, &
    & my_process_is_mpi_parallel, p_sum
#endif
!   USE mo_io_units,           ONLY: nnml, filename_max
!   USE mo_namelist,           ONLY: position_nml, open_nml, positioned

  IMPLICIT NONE

  PRIVATE

  ! !VERSION CONTROL:
  CHARACTER(LEN=*), PARAMETER :: version = '$Id$'

  !-------------------------------------------------------------------------

  ! NOTE: in order to get correct results make sure you provide the proper range_subset (ie, owned)!
#ifndef __ICON_GRID_GENERATOR__
  PUBLIC :: global_minmaxmean, subset_sum

  INTERFACE global_minmaxmean
    MODULE PROCEDURE globalspace_2D_minmaxmean
    MODULE PROCEDURE globalspace_3D_minmaxmean
  END INTERFACE global_minmaxmean

  INTERFACE subset_sum
    MODULE PROCEDURE globalspace_3D_sum
    MODULE PROCEDURE globalspace_3D_sum_max_level_array
  END INTERFACE subset_sum
#endif

  PUBLIC :: construct_statistic_objects, destruct_statistic_objects
  PUBLIC :: new_statistic, delete_statistic
  PUBLIC :: add_statistic_to, add_data_to
  PUBLIC :: min, max, mean
  PUBLIC :: max_statistic_of, min_statistic_of
  PUBLIC :: mean_statistic_of

  PUBLIC :: ADD_MAX_RATIO

  PUBLIC :: new, delete

  PUBLIC :: t_statistic
  !-------------------------------------------------------------------------
  TYPE :: t_statistic
    INTEGER :: id
  END TYPE t_statistic
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  ! Parameters
  INTEGER, PARAMETER :: ADD_VALUE = 1
  INTEGER, PARAMETER :: ADD_MAX_RATIO = 2

  !--------------------------------------------------------------
  ! TYPE definitions
  !> Basic statistics type
  TYPE t_data_statistics
    LOGICAL  :: is_active
    INTEGER  :: no_of_values
    REAL(wp) :: sum_of_values
    REAL(wp) :: max_of_values
    REAL(wp) :: min_of_values
    ! distribution stitistics
    ! at the moment just a
    INTEGER  :: no_of_bars
    REAL(wp) :: min_bars_value
    REAL(wp) :: max_bars_value
    INTEGER, POINTER  :: no_of_values_in_bar(:)

    INTEGER :: mode ! defines how we insert the values,
                          ! ADD_VALUE
                          ! ADD_MAX_RATIO, input twos values

  END TYPE t_data_statistics
  TYPE t_data_statistics_2d
    LOGICAL  :: is_active
    INTEGER  :: no_of_values
    REAL(wp), POINTER :: sum_of_values(:,:)
    REAL(wp), POINTER :: max_of_values(:,:)
    REAL(wp), POINTER :: min_of_values(:,:)
    INTEGER  :: no_of_bars
    REAL(wp) :: min_bars_value
    REAL(wp) :: max_bars_value
    INTEGER, POINTER  :: no_of_values_in_bar(:)
    INTEGER :: mode
  END TYPE t_data_statistics_2d

  TYPE t_data_statistics_3d
    LOGICAL  :: is_active
    INTEGER  :: no_of_values
    REAL(wp), POINTER :: sum_of_values(:,:,:)
    REAL(wp), POINTER :: max_of_values(:,:,:)
    REAL(wp), POINTER :: min_of_values(:,:,:)
    INTEGER  :: no_of_bars
    REAL(wp) :: min_bars_value
    REAL(wp) :: max_bars_value
    INTEGER, POINTER  :: no_of_values_in_bar(:)
    INTEGER :: mode
  END TYPE t_data_statistics_3d

  TYPE t_data_statistics_collection
    INTEGER, ALLOCATABLE :: stats_1d(:)
    INTEGER, ALLOCATABLE :: stats_2d(:)
    INTEGER, ALLOCATABLE :: stats_3d(:)
  END TYPE

  !--------------------------------------------------------------
  !> The maximum number of statistic objects.
  INTEGER, PARAMETER ::  max_no_of_statistic_objects = 50
  !> The array of the statistic objects.
  TYPE(t_data_statistics),    ALLOCATABLE, TARGET :: statistic_object(:)
  TYPE(t_data_statistics_2d), ALLOCATABLE, TARGET :: statistic_object_2d(:)
  TYPE(t_data_statistics_3d), ALLOCATABLE, TARGET :: statistic_object_3d(:)
  !> The number of allocated statistic objects.
  INTEGER :: no_of_allocated_statistics = 0        ! the size of the statistic objects array
  !> The number of actual active statistic objects.
  INTEGER :: active_statistics
  !> The maximum id of the active statistics
  INTEGER :: max_active_statistics
  !> True if the statistic object is active.

  INTERFACE new
    MODULE PROCEDURE new_statistic_operator
  END INTERFACE
  INTERFACE delete
    MODULE PROCEDURE delete_statistic_operator
  END INTERFACE
  INTERFACE add_data_to
    MODULE PROCEDURE add_data_one_value_real
    MODULE PROCEDURE add_data_one_value_int
  END INTERFACE
  INTERFACE min
    MODULE PROCEDURE min_statistic
  END INTERFACE
  INTERFACE max
    MODULE PROCEDURE max_statistic
  END INTERFACE
  INTERFACE mean
    MODULE PROCEDURE mean_statistic
  END INTERFACE
  INTERFACE new_statistic
    MODULE PROCEDURE new_statistic_no_bars
    MODULE PROCEDURE new_statistic_with_bars
    MODULE PROCEDURE new_statistic_with_pointers
  END INTERFACE
  INTERFACE init_statistics
    MODULE PROCEDURE init_statistic_object
    MODULE PROCEDURE init_statistic_object_with_pointers
  END INTERFACE
  INTERFACE add_statistic_to
    MODULE PROCEDURE add_statistic_one_value
    MODULE PROCEDURE add_statistic_two_values
  END INTERFACE

CONTAINS

#ifndef __ICON_GRID_GENERATOR__
  !-----------------------------------------------------------------------
  !>
  FUNCTION globalspace_2D_minmaxmean(values, range_subset) result(minmaxmean)
    REAL(wp), INTENT(in) :: values(:,:)
    TYPE(t_subset_range), TARGET, OPTIONAL :: range_subset
    REAL(wp) :: minmaxmean(3)

    REAL(wp) :: min_in_block, max_in_block, min_value, max_value, sum_value, global_number_of_values
    INTEGER :: block, startidx, endidx, number_of_values
    INTEGER :: communicator

    IF (PRESENT(range_subset)) THEN
      ! init the min, max values
      block = range_subset%start_block
      CALL get_index_range(range_subset, block, startidx, endidx)
      min_value = values(startidx, block)
      max_value = values(startidx, block)
      sum_value = 0
!ICON_OMP_PARALLEL_DO PRIVATE(block, startidx, endidx) FIRSTPRIVATE(min_in_block, max_in_block, sum_value) &
!ICON_OMP  reduction(min:min_value) reduction(max:max_value) reduction(sum:sum_value)
      DO block = range_subset%start_block, range_subset%end_block
        CALL get_index_range(range_subset, block, startidx, endidx)
        min_in_block = MINVAL(values(startidx:endidx, block))
        max_in_block = MAXVAL(values(startidx:endidx, block))
        min_value    = MIN(min_value, min_in_block)
        max_value    = MAX(max_value, max_in_block)
        sum_value    = sum_value + SUM(values(startidx:endidx, block))
      ENDDO
!ICON_OMP_END_PARALLEL_DO
      ! compute the total number of values
      number_of_values = range_subset%size

    ELSE
      min_value = values(1,1)
      max_value = values(1,1)
      sum_value = 0
!ICON_OMP_PARALLEL_DO PRIVATE(block, startidx, endidx) FIRSTPRIVATE(min_in_block, max_in_block, sum_value) &
!ICON_OMP  reduction(min:min_value) reduction(max:max_value) reduction(sum:sum_value)
      DO block = 1,  SIZE(values, 2)
        min_in_block = MINVAL(values(:, block))
        max_in_block = MAXVAL(values(:, block))
        min_value    = MIN(min_value, min_in_block)
        max_value    = MAX(max_value, max_in_block)
        sum_value    = sum_value + SUM(values(:, block))
      ENDDO
!ICON_OMP_END_PARALLEL_DO
      number_of_values = SIZE(values)
    ENDIF


    ! the global min, max, mean, is avaliable only to stdio process
    IF (my_process_is_mpi_parallel()) THEN
      communicator = get_my_mpi_work_communicator()
      minmaxmean(1) = p_min( min_value,  comm=communicator ) ! only mpi_all_reduce is avaliable
      minmaxmean(2) = p_max( max_value,  comm=communicator, root=process_mpi_stdio_id )

      ! these are avaliable to all processes
      global_number_of_values = p_sum( REAL(number_of_values,wp),  comm=communicator)
      minmaxmean(3) = p_sum( sum_value,  comm=communicator) / global_number_of_values


    ELSE

      minmaxmean(1) = min_value
      minmaxmean(2) = max_value
      minmaxmean(3) = sum_value / REAL(number_of_values, wp)

    ENDIF

  END FUNCTION globalspace_2D_minmaxmean
  !-----------------------------------------------------------------------


  !-----------------------------------------------------------------------
  !>
  FUNCTION globalspace_3D_minmaxmean(values, range_subset, start_level, end_level) result(minmaxmean)
    REAL(wp), INTENT(in) :: values(:,:,:)
    TYPE(t_subset_range), TARGET, OPTIONAL :: range_subset
    INTEGER, OPTIONAL :: start_level, end_level
    REAL(wp) :: minmaxmean(3)

    REAL(wp) :: min_in_block, max_in_block, min_value, max_value, sum_value, global_number_of_values
    INTEGER :: block, level, startidx, endidx, start_vertical, end_vertical, number_of_values
    INTEGER :: communicator
!    INTEGER :: idx
    CHARACTER(LEN=*), PARAMETER :: method_name='mo_statistics:globalspace_3D_minmaxmean'


    IF (PRESENT(start_level)) THEN
      start_vertical = start_level
    ELSE
      start_vertical = 1
    ENDIF
    IF (PRESENT(end_level)) THEN
      end_vertical = start_level
    ELSE
      end_vertical = SIZE(values, 2)
    ENDIF
    IF (start_vertical > end_vertical) &
      CALL finish(method_name, "start_vertical > end_vertical")

    IF (PRESENT(range_subset)) THEN
      ! init the min, max values
      block=range_subset%start_block
      CALL get_index_range(range_subset, block, startidx, endidx)
      min_value = values(startidx, start_vertical, block)
      max_value = values(startidx, start_vertical, block)
      sum_value = 0

!ICON_OMP_PARALLEL_DO PRIVATE(block, startidx, endidx) FIRSTPRIVATE(min_in_block, max_in_block, sum_value) &
!ICON_OMP  reduction(min:min_value) reduction(max:max_value) reduction(sum:sum_value)
      DO block = range_subset%start_block, range_subset%end_block
        CALL get_index_range(range_subset, block, startidx, endidx)
        DO level = start_vertical, end_vertical
!          DO idx = startidx, endidx
!            WRITE(0,*) "checking ", block, level, idx
!            WRITE(0,*) "                  =   ", values(idx, level, block)
!          ENDDO
          min_in_block = MINVAL(values(startidx:endidx, level, block))
          max_in_block = MAXVAL(values(startidx:endidx, level, block))
          min_value    = MIN(min_value, min_in_block)
          max_value    = MAX(max_value, max_in_block)
          sum_value    = sum_value + SUM(values(startidx:endidx, level, block))
        ENDDO
      ENDDO
!ICON_OMP_END_PARALLEL_DO

      IF ((range_subset%end_block - range_subset%start_block) > 1) THEN
        number_of_values = (range_subset%end_block - range_subset%start_block -1) * range_subset%block_size
      ELSE
        number_of_values = 0
      ENDIF
      number_of_values = (number_of_values + range_subset%end_index + &
        & (range_subset%block_size - range_subset%start_index + 1)) * &
        & (end_vertical - start_vertical + 1)

    ELSE

      min_value      = values(1, start_vertical, 1)
      max_value      = values(1, start_vertical, 1)
      sum_value = 0

!ICON_OMP_PARALLEL_DO PRIVATE(block, startidx, endidx) FIRSTPRIVATE(min_in_block, max_in_block, sum_value) &
!ICON_OMP  reduction(min:min_value) reduction(max:max_value) reduction(sum:sum_value)
      DO block = 1,  SIZE(values, 3)
        DO level = start_vertical, end_vertical
          min_in_block = MINVAL(values(:, level, block))
          max_in_block = MAXVAL(values(:, level, block))
          min_value    = MIN(min_value, min_in_block)
          max_value    = MAX(max_value, max_in_block)
          sum_value    = sum_value + SUM(values(:, level, block))
        ENDDO
      ENDDO
!ICON_OMP_END_PARALLEL_DO

      number_of_values = SIZE(values(:,start_vertical:end_vertical,:))
    ENDIF

    ! the global min, max is avaliable only to stdio process
    IF (my_process_is_mpi_parallel()) THEN
      communicator = get_my_mpi_work_communicator()
      minmaxmean(1) = p_min( min_value,  comm=communicator ) ! only mpi_all_reduce is avaliable
      minmaxmean(2) = p_max( max_value,  comm=communicator, root=process_mpi_stdio_id )

      ! these are avaliable to all processes
      global_number_of_values = p_sum( REAL(number_of_values,wp),  comm=communicator)
      minmaxmean(3) = p_sum( sum_value,  comm=communicator) / global_number_of_values

    ELSE

      minmaxmean(1) = min_value
      minmaxmean(2) = max_value
      minmaxmean(3) = sum_value / REAL(number_of_values, wp)

    ENDIF

  END FUNCTION globalspace_3D_minmaxmean
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !>
  FUNCTION globalspace_3D_sum(values, indexed_subset, start_level, end_level, weights, &
    & subset_indexed_weights) result(total_sum)
    REAL(wp) :: values(:,:,:)
    TYPE(t_subset_indexed), TARGET, OPTIONAL :: indexed_subset
    INTEGER,  OPTIONAL  :: start_level, end_level
    REAL(wp), OPTIONAL  :: weights(:,:,:)
    REAL(wp), OPTIONAL  :: subset_indexed_weights(:,:)  ! weights but indexed but the subset index
                                                        ! dim: (vertical_levels, indexed_subset%size)
    REAL(wp) :: total_sum

    REAL(wp) :: sum_value
    INTEGER :: i, block, idx, level,  start_vertical, end_vertical
    INTEGER :: communicator
!    INTEGER :: idx
    CHARACTER(LEN=*), PARAMETER :: method_name='mo_statistics:globalspace_3D_sum'


    IF (PRESENT(start_level)) THEN
      start_vertical = start_level
    ELSE
      start_vertical = 1
    ENDIF
    IF (PRESENT(end_level)) THEN
      end_vertical = start_level
    ELSE
      end_vertical = SIZE(values, 2)
    ENDIF
    IF (start_vertical > end_vertical) &
      CALL finish(method_name, "start_vertical > end_vertical")

    IF (.not. PRESENT(indexed_subset)) &
      & CALL finish(method_name,  "Currently requires indexed_subset parameter. Abort.")

    ! init the min, max values
    sum_value = 0
    IF (PRESENT(weights)) THEN

      DO i=1, indexed_subset%size
        block = indexed_subset%block(i)
        idx = indexed_subset%idx(i)
        DO level = start_vertical, end_vertical
          sum_value  = sum_value + values(idx, level, block) * weights(idx, level, block)
        ENDDO
      ENDDO

    ELSEIF (PRESENT(subset_indexed_weights)) THEN

      DO i=1, indexed_subset%size
        block = indexed_subset%block(i)
        idx = indexed_subset%idx(i)
        DO level = start_vertical, end_vertical
          sum_value  = sum_value + values(idx, level, block) * subset_indexed_weights(level, i)
        ENDDO
      ENDDO

    ELSE

      DO i=1, indexed_subset%size
        block = indexed_subset%block(i)
        idx = indexed_subset%idx(i)
        DO level = start_vertical, end_vertical
          sum_value  = sum_value + values(idx, level, block)
        ENDDO
      ENDDO

    ENDIF

    ! the global min, max is avaliable only to stdio process
    IF (my_process_is_mpi_parallel()) THEN

      communicator = indexed_subset%patch%work_communicator
      ! these are avaliable to all processes
      total_sum = p_sum( sum_value,  comm=communicator)

    ELSE

      total_sum = sum_value

    ENDIF

  END FUNCTION globalspace_3D_sum
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !>
  FUNCTION globalspace_3D_sum_max_level_array(values, indexed_subset, weights, &
    & subset_indexed_weights, max_level_array) result(total_sum)
    REAL(wp) :: values(:,:,:)
    TYPE(t_subset_indexed), TARGET, OPTIONAL :: indexed_subset
    REAL(wp), OPTIONAL  :: weights(:,:,:)
    REAL(wp), OPTIONAL  :: subset_indexed_weights(:,:)  ! weights but indexed but the subset index
                                                        ! dim: (vertical_levels, indexed_subset%size)
    INTEGER  :: max_level_array(:,:)
    REAL(wp) :: total_sum

    REAL(wp) :: sum_value
    INTEGER :: i, block, idx, level,  start_vertical, end_vertical
    INTEGER :: communicator
!    INTEGER :: idx
    CHARACTER(LEN=*), PARAMETER :: method_name='mo_statistics:globalspace_3D_sum'

    IF (.not. PRESENT(indexed_subset)) &
      & CALL finish(method_name,  "Currently requires indexed_subset parameter. Abort.")

    ! init the min, max values
    sum_value = 0
    IF (PRESENT(weights)) THEN

      DO i=1, indexed_subset%size
        block = indexed_subset%block(i)
        idx = indexed_subset%idx(i)
        DO level = 1, max_level_array(idx, block)
          sum_value  = sum_value + values(idx, level, block) * weights(idx, level, block)
        ENDDO
      ENDDO

    ELSEIF (PRESENT(subset_indexed_weights)) THEN

      DO i=1, indexed_subset%size
        block = indexed_subset%block(i)
        idx = indexed_subset%idx(i)
        DO level = 1, max_level_array(idx, block)
          sum_value  = sum_value + values(idx, level, block) * subset_indexed_weights(level, i)
        ENDDO
      ENDDO

    ELSE

      DO i=1, indexed_subset%size
        block = indexed_subset%block(i)
        idx = indexed_subset%idx(i)
        DO level = 1, max_level_array(idx, block)
          sum_value  = sum_value + values(idx, level, block)
        ENDDO
      ENDDO

    ENDIF

    ! the global min, max is avaliable only to stdio process
    IF (my_process_is_mpi_parallel()) THEN

      communicator = indexed_subset%patch%work_communicator
      ! these are avaliable to all processes
      total_sum = p_sum( sum_value,  comm=communicator)

    ELSE

      total_sum = sum_value

    ENDIF

  END FUNCTION globalspace_3D_sum_max_level_array
  !-----------------------------------------------------------------------
#endif

  !-----------------------------------------------------------------------
  !>
  SUBROUTINE new_statistic_operator(statistic)
    TYPE(t_statistic), INTENT(inout) :: statistic

    statistic%id = new_statistic_no_bars()

  END SUBROUTINE new_statistic_operator
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !>
  SUBROUTINE delete_statistic_operator(statistic)
    TYPE(t_statistic), INTENT(inout) :: statistic

    CALL delete_statistic(statistic%id)
    statistic%id = -1

  END SUBROUTINE delete_statistic_operator
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !>
  !! Creates a new statistics object and returns its id.
  !! The statistic arrays are not allocated.
  INTEGER FUNCTION new_statistic_with_bars(no_of_bars, min_value, max_value, mode)
    INTEGER,  INTENT(in) :: no_of_bars
    REAL(wp), INTENT(in) :: min_value, max_value
    INTEGER,  INTENT(in), OPTIONAL :: mode

    IF (PRESENT(mode)) THEN
      new_statistic_with_bars = new_statistic(mode)
    ELSE
      new_statistic_with_bars = new_statistic(mode)
    ENDIF
    CALL construct_statistic_bars(new_statistic_with_bars, no_of_bars, min_value, max_value)

  END FUNCTION new_statistic_with_bars
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !>
  !! Creates a new statistics object and returns its id.
  !! The statistic arrays are not allocated.
  INTEGER FUNCTION new_statistic_no_bars(mode)
    INTEGER,  INTENT(in), OPTIONAL :: mode

    INTEGER :: i

    ! check if the statistic objects have been construced
    IF (no_of_allocated_statistics == 0) THEN
      CALL construct_statistic_objects()
    ENDIF

    IF (max_active_statistics /= active_statistics) THEN
      ! we have a hole (inactive statistic) in the list of max_active_statistics
      DO i = 1, max_active_statistics
        IF (.not. statistic_object(i)%is_active) THEN
          new_statistic_no_bars = i
          EXIT
        ENDIF
      ENDDO
    ELSE
      ! add a new statistic object
      IF (max_active_statistics >= no_of_allocated_statistics) THEN
        CALL finish('new_statistic_no_bars', 'exceeded no_of_allocated_statistics');
      ENDIF
      max_active_statistics = max_active_statistics + 1
      new_statistic_no_bars = max_active_statistics
    ENDIF

    active_statistics = active_statistics + 1
    CALL init_statistic_object(new_statistic_no_bars)

    IF (PRESENT(mode)) THEN
      CALL set_mode(new_statistic_no_bars, mode)
    ENDIF


  END FUNCTION new_statistic_no_bars
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !>
  !! Creates a new statistics object and returns its id.
  !! The statistic arrays are not allocated.
  INTEGER FUNCTION new_statistic_with_pointers(target,source,mode)
    REAL(wp), TARGET               :: target,source
    INTEGER,  INTENT(in), OPTIONAL :: mode
    ! GZ: temporary fix to get the code compiled
    new_statistic_with_pointers = 0
  END FUNCTION new_statistic_with_pointers

  !-----------------------------------------------------------------------
  !>
  SUBROUTINE set_mode(statistic_id, mode)
    INTEGER, INTENT(in) :: statistic_id, mode

    CALL check_active_statistic_id(statistic_id)

    SELECT CASE (mode)

    CASE (ADD_VALUE, ADD_MAX_RATIO)
      statistic_object(statistic_id)%mode = mode

    CASE DEFAULT
      CALL finish('set_mode', 'Unkown mode')

    END SELECT
  END SUBROUTINE set_mode
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !>
  !! Deletes the statistic object.
  !! The statistic arrays are deallocated.
  !! The statistic_id is invalid after the delete call.
  !! It maybe associated with another statistic in the future.
  SUBROUTINE delete_statistic(statistic_id)
    INTEGER, INTENT(in) :: statistic_id

    CALL check_active_statistic_id(statistic_id)
!     CALL deallocate_statistic_object(statistic_id)

    CALL destruct_statistic_bars(statistic_id)
    statistic_object(statistic_id)%is_active = .false.
    active_statistics = active_statistics - 1

  END SUBROUTINE delete_statistic
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !>
  SUBROUTINE add_data_one_value_real(statistic, value)
    TYPE(t_statistic), INTENT(inout) :: statistic
    REAL(wp), INTENT(in) :: value

    CALL add_statistic_one_value(statistic%id, value)
  END SUBROUTINE add_data_one_value_real
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !>
  SUBROUTINE add_data_one_value_int(statistic, value)
    TYPE(t_statistic), INTENT(inout) :: statistic
    INTEGER, INTENT(in) :: value

    CALL add_statistic_one_value(statistic%id, REAL(value,wp))
  END SUBROUTINE add_data_one_value_int
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !>
  SUBROUTINE add_statistic_one_value(statistic_id, value)
    INTEGER, INTENT(in) :: statistic_id
    REAL(wp), INTENT(in) :: value

    CALL check_active_statistic_id(statistic_id)

    statistic_object(statistic_id)%sum_of_values  = &
      & statistic_object(statistic_id)%sum_of_values + value
    statistic_object(statistic_id)%max_of_values  = &
      & MAX(statistic_object(statistic_id)%max_of_values, value)
    statistic_object(statistic_id)%min_of_values  = &
      & MIN(statistic_object(statistic_id)%min_of_values, value)

    statistic_object(statistic_id)%no_of_values   = &
      & statistic_object(statistic_id)%no_of_values + 1

  END SUBROUTINE add_statistic_one_value
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !>
  SUBROUTINE add_statistic_two_values(statistic_id, value1, value2)
    INTEGER, INTENT(in) :: statistic_id
    REAL(wp), INTENT(in) :: value1, value2

!     REAL(wp) :: max_value, min_value

!     CALL check_active_statistic_id(statistic_id)

    SELECT CASE (statistic_object(statistic_id)%mode)

    CASE(ADD_MAX_RATIO)

      CALL add_statistic_one_value(statistic_id, &
        MAX(value1, value2) / MIN(value1, value2))

     CASE DEFAULT
      CALL finish('set_mode', 'Unkown mode')

    END SELECT

  END SUBROUTINE add_statistic_two_values
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !>
  REAL(wp) FUNCTION min_statistic(statistic)
    TYPE(t_statistic), INTENT(in) :: statistic

    CALL check_active_statistic_id(statistic%id)
    min_statistic = statistic_object(statistic%id)%min_of_values

  END FUNCTION min_statistic
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !>
  REAL(wp) FUNCTION max_statistic(statistic)
    TYPE(t_statistic), INTENT(in) :: statistic

    CALL check_active_statistic_id(statistic%id)
    max_statistic = statistic_object(statistic%id)%max_of_values

  END FUNCTION max_statistic
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !>
  REAL(wp) FUNCTION mean_statistic(statistic)
    TYPE(t_statistic), INTENT(in) :: statistic

    CALL check_active_statistic_id(statistic%id)
    mean_statistic = statistic_object(statistic%id)%sum_of_values / &
      & REAL(statistic_object(statistic%id)%no_of_values, wp)

  END FUNCTION mean_statistic
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !>
  REAL(wp) FUNCTION min_statistic_of(statistic_id)
    INTEGER, INTENT(in) :: statistic_id

    CALL check_active_statistic_id(statistic_id)
    min_statistic_of = statistic_object(statistic_id)%min_of_values

  END FUNCTION min_statistic_of
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !>
  REAL(wp) FUNCTION mean_statistic_of(statistic_id)
    INTEGER, INTENT(in) :: statistic_id

    CALL check_active_statistic_id(statistic_id)
    mean_statistic_of = statistic_object(statistic_id)%sum_of_values / &
      & REAL(statistic_object(statistic_id)%no_of_values, wp)

  END FUNCTION mean_statistic_of
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !>
  REAL(wp) FUNCTION max_statistic_of(statistic_id)
    INTEGER, INTENT(in) :: statistic_id

    CALL check_active_statistic_id(statistic_id)
    max_statistic_of = statistic_object(statistic_id)%max_of_values

  END FUNCTION max_statistic_of
  !-----------------------------------------------------------------------


  !-----------------------------------------------------------------------
  !>
  !! Deletes all statistic objects
  !! Note: Should be called only at the stop of a program.
  SUBROUTINE destruct_statistic_objects ()

    INTEGER :: i

    IF (no_of_allocated_statistics == 0) THEN
      CALL warning('destruct_statistic_objects', &
        & 'statistic_objects have not been constructed');
      RETURN
    ENDIF

    DO i=1,max_active_statistics
      IF (statistic_object(i)%is_active) THEN
        CALL delete_statistic(i)
      ENDIF
    ENDDO

    DEALLOCATE(statistic_object)

    no_of_allocated_statistics  = 0
    active_statistics     = 0
    max_active_statistics = 0

  END SUBROUTINE destruct_statistic_objects
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !>
  !! "Class construction" method.
  !! Called once to initiallize the "class".
  SUBROUTINE construct_statistic_objects()

    INTEGER :: return_status,i

    CHARACTER(*), PARAMETER :: method_name = "construct_statistic_objects"

    IF (no_of_allocated_statistics /= 0) THEN
      ! do some consistency checks, return if everything ok
      IF (no_of_allocated_statistics /= max_no_of_statistic_objects) THEN
        CALL finish(method_name, &
          & 'no_of_allocated_statistics/= max_no_of_statistic_objects');
      ENDIF
      IF (.not. ALLOCATED(statistic_object)) THEN
        CALL finish(method_name, '.NOT. ALLOCATED(statistic_object)');
      ENDIF
      IF (SIZE(statistic_object) /= no_of_allocated_statistics) THEN
        CALL finish(method_name, &
          & 'SIZE(statistic_object) /= no_of_allocated_statistics');
      ENDIF

      RETURN
    ENDIF

    IF (ALLOCATED(statistic_object)) THEN
      CALL finish(method_name, 'ALLOCATED(statistic_object)');
    ENDIF

    no_of_allocated_statistics  = max_no_of_statistic_objects
    active_statistics     = 0
    max_active_statistics = 0

    ALLOCATE(statistic_object(no_of_allocated_statistics),stat=return_status)
    IF (return_status /= 0) THEN
      CALL finish(method_name, 'failed to ALLOCATE(statistic_object)');
    ENDIF

    DO i=1,no_of_allocated_statistics
      statistic_object(i)%is_active = .false.
    ENDDO

  END SUBROUTINE construct_statistic_objects
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !>
  !! Initializes a new statistic object
  SUBROUTINE init_statistic_object(statistic_id)
    INTEGER, INTENT(in) :: statistic_id

    statistic_object(statistic_id)%no_of_values   = 0
    statistic_object(statistic_id)%sum_of_values  = 0._wp
    statistic_object(statistic_id)%max_of_values  = &
      & -HUGE(statistic_object(statistic_id)%max_of_values)
    statistic_object(statistic_id)%min_of_values  = &
      & HUGE(statistic_object(statistic_id)%min_of_values)
    statistic_object(statistic_id)%mode     = ADD_VALUE

    statistic_object(statistic_id)%no_of_bars     = 0
    statistic_object(statistic_id)%min_bars_value = 0._wp
    statistic_object(statistic_id)%max_bars_value = 0._wp
    NULLIFY(statistic_object(statistic_id)%no_of_values_in_bar)

    statistic_object(statistic_id)%is_active = .true.

  END SUBROUTINE init_statistic_object
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !>
  !! Initializes a new statistic object with source and target pointers
  SUBROUTINE init_statistic_object_with_pointers(statistic_id,source,summation,mean,maximum,minimum)
    INTEGER, INTENT(in)            :: statistic_id
    REAL(wp), POINTER              :: source
    REAL(wp), POINTER, OPTIONAL    :: summation,mean,maximum,minimum

    statistic_object(statistic_id)%no_of_values = 0
   !IF (PRESENT(summation)) THEN
   !  statistic_object(statistic_id)%sum_of_values => summation
   !END IF
   !IF (PRESENT(maximum)) THEN
   !  statistic_object(statistic_id)%max_of_values => maximum
   !END IF
   !IF (PRESENT(minimum)) THEN
   !  statistic_object(statistic_id)%min_of_values => minimum
   !END IF
    statistic_object(statistic_id)%mode = ADD_VALUE
  END SUBROUTINE init_statistic_object_with_pointers

  !-----------------------------------------------------------------------
  !>
  !! Reset values of the statistic object.
  SUBROUTINE reset_statistic(statistic_id)
    INTEGER, INTENT(in) :: statistic_id
    CALL init_statistic_object(statistic_id)
  END SUBROUTINE reset_statistic

  !-----------------------------------------------------------------------
  !>
  !! Initializes a new statistic object
  SUBROUTINE construct_statistic_bars(statistic_id, no_of_bars, min_value, max_value)
    INTEGER,  INTENT(in) :: statistic_id, no_of_bars
    REAL(wp), INTENT(in) :: min_value, max_value

    INTEGER :: return_status
    CHARACTER(*), PARAMETER :: method_name = "construct_statistic_bars"

    CALL check_active_statistic_id(statistic_id)
    IF (no_of_bars < 2) &
      CALL finish(method_name,'no_of_bars < 2');
    IF (min_value >= max_value) &
      CALL finish(method_name,'min_value >= max_value');

    statistic_object(statistic_id)%no_of_bars     = no_of_bars
    statistic_object(statistic_id)%min_bars_value = min_value
    statistic_object(statistic_id)%max_bars_value = max_value

    ALLOCATE(statistic_object(statistic_id)%no_of_values_in_bar(no_of_bars),stat=return_status)
    IF (return_status /= 0) THEN
      CALL finish(method_name, 'failed to ALLOCATE(statistic_object)');
    ENDIF

  END SUBROUTINE construct_statistic_bars
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !>
  !! Initializes a new statistic object
  SUBROUTINE destruct_statistic_bars(statistic_id)
    INTEGER,  INTENT(in) :: statistic_id

!     INTEGER :: return_status
!     CHARACTER(*), PARAMETER :: method_name = "destruct_statistic_bars"

    CALL check_active_statistic_id(statistic_id)

    IF (ASSOCIATED(statistic_object(statistic_id)%no_of_values_in_bar)) THEN
      DEALLOCATE(statistic_object(statistic_id)%no_of_values_in_bar)
    ENDIF

    statistic_object(statistic_id)%no_of_bars     = 0
    NULLIFY(statistic_object(statistic_id)%no_of_values_in_bar)

  END SUBROUTINE destruct_statistic_bars
  !-----------------------------------------------------------------------


  !-----------------------------------------------------------------------
  !>
  !! Checks if a the statistic object associated with the statistic_id is active.
  !! If the statistic object is not active the program stops with an error.
  SUBROUTINE check_active_statistic_id (statistic_id)
    INTEGER :: statistic_id

    IF (statistic_id > max_active_statistics) THEN
      CALL finish('get_statistic', 'statistic_id > max_active_statistics');
    ENDIF
    IF (.not. statistic_object(statistic_id)%is_active) THEN
      CALL finish('get_statistic', 'statistic is not active');
    ENDIF
  END SUBROUTINE check_active_statistic_id
  !-----------------------------------------------------------------------


END MODULE mo_statistics
