!>
!! This module contains utilty types
!!
!! @par Revision History
!!  Created by Leonidas Linardakis, MPIM (2012-03-06)
!!
!! @par Copyright
!! 2002-2007 by DWD and MPI-M
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
!!
!!
MODULE mo_util_subset

  USE mo_exception,    ONLY: warning

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: version = '$Id$'

  PUBLIC :: t_subset_range, t_subset_range_index

  PUBLIC :: fill_subset, get_index_range

  !----------------------------------------------------
  !> Defines a subset in a range (in terms of blocks)
  TYPE :: t_subset_range
    
    INTEGER :: start_block
    INTEGER :: start_index
    INTEGER :: end_block
    INTEGER :: end_index
    INTEGER :: block_size

    INTEGER :: no_of_holes ! the number of holes in the subset
    
  END TYPE
  !----------------------------------------------------
  

  !----------------------------------------------------
  !> Defines an index for a subset_range 
  TYPE :: t_subset_range_index
  
    INTEGER, POINTER :: block ! the current block in the subset
    INTEGER, POINTER :: index ! the current index in the subset
    
    INTEGER :: start_index ! the current start index within the current block,
    INTEGER :: end_index   ! the current end index within the current block,

    TYPE(t_subset_range), POINTER :: subset_range
  
  END TYPE
  !----------------------------------------------------

CONTAINS

  !----------------------------------------------------
  !>
  ! Fills the subset_range with the indexes which statisfy
  !     start_level <= level <= end_level
  !
  ! Assumes that invalid places in level have value outside the
  ! start_level - end_level
  !
  ! Assumes that level is of the shape (1:,1:)
  SUBROUTINE fill_subset(subset_range, level, start_level, end_level)
    TYPE(t_subset_range), INTENT(inout) :: subset_range
    INTEGER, INTENT(in) :: level(:,:), start_level, end_level

    INTEGER :: levels_size(2) 
    INTEGER :: block, index_in_block, start_index, end_index
    
    CHARACTER(*), PARAMETER :: method_name = "mo_util_subset:fill_subset"

    levels_size = SHAPE(level)

    subset_range%start_block = -1
    subset_range%start_index = -1
    subset_range%end_block   = -2
    subset_range%end_index   = -2
    subset_range%block_size  = levels_size(1)
    subset_range%no_of_holes = 0    
    
    DO block=1,levels_size(2)
      DO index_in_block=1,subset_range%block_size

        IF (level(index_in_block, block) >= start_level .AND. &
            level(index_in_block, block) <= end_level) THEN
          ! we found an elemant in range
          IF (subset_range%start_block < 0) THEN
            ! this is the first element
            subset_range%start_block = block
            subset_range%start_index = index_in_block
          ENDIF
          ! this is the up-to-now last element
          subset_range%end_block = block
          subset_range%end_index = index_in_block
        ENDIF
        
      ENDDO
    ENDDO

!     IF (subset_range%start_block == -1) THEN
!       CALL warning(method_name, "Empty range subset")
!     ENDIF
    
    IF (subset_range%start_block > -1) THEN
      ! count the holes      
      DO block = subset_range%start_block, subset_range%end_block
      
        start_index = 1
        end_index = subset_range%block_size
        IF (block == subset_range%start_block) start_index = subset_range%start_index
        IF (block == subset_range%end_block)   end_index = subset_range%end_index

        DO index_in_block = start_index, end_index
        
          IF (level(index_in_block, block) < start_level .OR. &
              level(index_in_block, block) > end_level) THEN
            ! element is a hole in the range
            subset_range%no_of_holes = subset_range%no_of_holes + 1
          ENDIF  

        ENDDO
      ENDDO
    ENDIF

!     IF (subset_range%no_of_holes > 0) THEN
!       CALL warning(method_name, "We have holes in the range subset")      
!     ENDIF
        
  END SUBROUTINE fill_subset
  !----------------------------------------------------
    
  
  !----------------------------------------------------
  SUBROUTINE get_index_range(subset_range, current_block, start_index, end_index)
    TYPE(t_subset_range), INTENT(in) :: subset_range
    INTEGER, INTENT(in) :: current_block
    INTEGER, INTENT(out) :: start_index, end_index
    
    start_index = 1
    end_index = subset_range%block_size

    IF (current_block == subset_range%start_block) &
      start_index = subset_range%start_index
    IF (current_block == subset_range%end_block) &
      end_index = subset_range%end_index

  END SUBROUTINE get_index_range
  !----------------------------------------------------

END MODULE mo_util_subset

