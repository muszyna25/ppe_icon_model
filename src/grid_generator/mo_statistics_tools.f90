!-------------------------------------------------------------------------------------
!>
!! Set of methods for simple statistics
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
MODULE mo_statistics_tools
  !-------------------------------------------------------------------------
  USE mo_kind,               ONLY: wp
  USE mo_exception,          ONLY: message_text, message, warning, finish
!   USE mo_io_units,           ONLY: nnml, filename_max
!   USE mo_namelist,           ONLY: position_nml, open_nml, positioned

  IMPLICIT NONE

  PRIVATE

  ! !VERSION CONTROL:
  CHARACTER(LEN=*), PARAMETER :: version = '$Id$'

  PUBLIC :: construct_statistic_objects, destruct_statistic_objects
  PUBLIC :: new_statistic, delete_statistic
  PUBLIC :: add_statistic_to, max_statistic_of, min_statistic_of
  PUBLIC :: mean_statistic_of

  PUBLIC :: ADD_MAX_RATIO

  !-------------------------------------------------------------------------
  ! Parameters
  INTEGER, PARAMETER :: ADD_VALUE = 1
  INTEGER, PARAMETER :: ADD_MAX_RATIO = 2
  
  !--------------------------------------------------------------
  ! TYPE definitions
  !> Basic statistics type
  TYPE t_statistic
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
    
  END TYPE t_statistic

  !--------------------------------------------------------------
  !> The maximum number of statistic objects.
  INTEGER, PARAMETER ::  max_no_of_statistic_objects = 50
  !> The array of the statistic objects.
  TYPE(t_statistic), ALLOCATABLE, TARGET :: statistic_object(:)
  !> The number of allocated statistic objects.
  INTEGER :: no_of_allocated_statistics = 0        ! the size of the statistic objects array
  !> The number of actual active statistic objects.
  INTEGER :: active_statistics
  !> The maximum id of the active statistics
  INTEGER :: max_active_statistics
  !> True if the statistic object is active.
  
  INTERFACE new_statistic
    MODULE PROCEDURE new_statistic_no_bars
    MODULE PROCEDURE new_statistic_with_bars
  END INTERFACE
  
  INTERFACE add_statistic_to
    MODULE PROCEDURE add_statistic_one_value
    MODULE PROCEDURE add_statistic_two_values
  END INTERFACE

CONTAINS

  !-----------------------------------------------------------------------
  !>
  !! Creates a new statistics object and returns its id.
  !! The statistic arrays are not allocated.
  INTEGER FUNCTION new_statistic_with_bars(no_of_bars, min_value, max_value, mode)
    INTEGER,  INTENT(in) :: no_of_bars
    REAL(wp), INTENT(in) :: min_value, max_value
    INTEGER,  INTENT(in), OPTIONAL :: mode

    INTEGER :: statistic_id

    IF (PRESENT(mode)) THEN
      statistic_id = new_statistic_no_bars(mode)
    ELSE
      statistic_id = new_statistic_no_bars(mode)
    ENDIF
    CALL construct_statistic_bars(statistic_id, no_of_bars, min_value, max_value)
    
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
  SUBROUTINE construct_statistic_objects ()

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


END MODULE mo_statistics_tools

