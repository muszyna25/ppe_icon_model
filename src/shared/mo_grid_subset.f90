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
MODULE mo_grid_subset

  USE mo_exception,    ONLY: warning, finish
  USE mo_model_domain, ONLY: t_patch, t_subset_range, t_subset_range_index
  
  IMPLICIT NONE
  INCLUDE 'netcdf.inc'

  CHARACTER(len=*), PARAMETER, PRIVATE :: version = '$Id$'

  PUBLIC :: t_subset_range, t_subset_range_index

  PUBLIC :: fill_subset, get_index_range
  PUBLIC :: read_subset, write_subset

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
  SUBROUTINE fill_subset(subset_range, patch, level, start_level, end_level, subset_name)
    TYPE(t_subset_range), INTENT(inout) :: subset_range
    TYPE(t_patch), TARGET, INTENT(in) :: patch  ! nag does not return the values in subset_range
    CHARACTER(len=32), OPTIONAL :: subset_name
                                                ! unless the patch is declared INTENT(in)!
    INTEGER, INTENT(in) :: level(:,:), start_level, end_level

    INTEGER :: levels_size(2) 
    INTEGER :: block, index_in_block, start_index, end_index
    
    CHARACTER(*), PARAMETER :: method_name = "mo_grid_subset:fill_subset"

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

    subset_range%patch => patch
    IF (PRESENT(subset_name)) &
      & subset_range%name = TRIM(subset_name)
    
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
  
  !-------------------------------------------------------------------------
  SUBROUTINE read_subset(ncid, subset_range, patch)
    INTEGER, INTENT(in) :: ncid
    TYPE(t_subset_range), INTENT(inout) :: subset_range
    TYPE(t_patch), TARGET, INTENT(in) :: patch  ! nag does not return the values in subset_range
    
    INTEGER :: netcd_status
    CHARACTER(*), PARAMETER :: method_name = "read_subset_range"
            
    netcd_status = nf_get_att_int(ncid, nf_global,TRIM(subset_range%name)//'.start_block', subset_range%start_block)
    IF (netcd_status /= nf_noerr) &
      & CALL finish(method_name, "Could not read start_block")
    netcd_status = nf_get_att_int(ncid, nf_global,TRIM(subset_range%name)//'.start_index', subset_range%start_index)
    IF (netcd_status /= nf_noerr) &
      & CALL finish(method_name, "Could not read start_index")
    
    netcd_status = nf_get_att_int(ncid, nf_global,TRIM(subset_range%name)//'.end_block', subset_range%end_block)
    IF (netcd_status /= nf_noerr) &
      & CALL finish(method_name, "Could not read end_block")
    netcd_status = nf_get_att_int(ncid, nf_global,TRIM(subset_range%name)//'.end_index', subset_range%end_index)
    IF (netcd_status /= nf_noerr) &
      & CALL finish(method_name, "Could not read end_index")
        
    netcd_status = nf_get_att_int(ncid, nf_global,TRIM(subset_range%name)//'.block_size', subset_range%block_size)
    IF (netcd_status /= nf_noerr) &
      & CALL finish(method_name, "Could not read block_size")
    
    netcd_status = nf_get_att_int(ncid, nf_global,TRIM(subset_range%name)//'.entity_type', subset_range%entity_type)
    IF (netcd_status /= nf_noerr) &
      & CALL finish(method_name, "Could not read entity_type")
    
    netcd_status = nf_get_att_int(ncid, nf_global,TRIM(subset_range%name)//'.no_of_holes', subset_range%no_of_holes)
    IF (netcd_status /= nf_noerr) &
      & CALL finish(method_name, "Could not read no_of_holes")
    
    netcd_status = nf_get_att_int(ncid, nf_global,TRIM(subset_range%name)//'.is_in_domain', subset_range%is_in_domain)
    IF (netcd_status /= nf_noerr) &
      & CALL finish(method_name, "Could not read is_in_domain")
    
    subset_range%patch => patch
    
    
  END SUBROUTINE read_subset
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE write_subset(ncid, subset_range)
    INTEGER, INTENT(in) :: ncid
    TYPE(t_subset_range), INTENT(inout) :: subset_range
    
    INTEGER :: netcd_status
    CHARACTER(*), PARAMETER :: method_name = "write_subset_range"
                    
    CALL nf(nf_put_att_int(ncid, nf_global,TRIM(subset_range%name)//'.start_block', nf_int, 1,     &
      &  subset_range%start_block))
    CALL nf(nf_put_att_int(ncid, nf_global,TRIM(subset_range%name)//'.start_index', nf_int, 1,     &
      &  subset_range%start_index))
    
    CALL nf(nf_put_att_int(ncid, nf_global,TRIM(subset_range%name)//'.end_block',   nf_int, 1,     &
      & subset_range%end_block))
    CALL nf(nf_put_att_int(ncid, nf_global,TRIM(subset_range%name)//'.end_index',   nf_int, 1,     &
      & subset_range%end_index))
        
    CALL nf(nf_put_att_int(ncid, nf_global,TRIM(subset_range%name)//'.block_size',  nf_int, 1,     &
      & subset_range%block_size))
    
    CALL nf(nf_put_att_int(ncid, nf_global,TRIM(subset_range%name)//'.entity_type', nf_int, 1,     &
      & subset_range%entity_type))
    
    CALL nf(nf_put_att_int(ncid, nf_global,TRIM(subset_range%name)//'.no_of_holes', nf_int, 1,     &
      & subset_range%no_of_holes))
    
    CALL nf(nf_put_att_int(ncid, nf_global,TRIM(subset_range%name)//'.is_in_domain',nf_int, 1,     &
      & subset_range%is_in_domain))
            
  END SUBROUTINE write_subset
  !-------------------------------------------------------------------------
  
  !--------------------------------------------------------------------
  SUBROUTINE nf(return_status)
    INTEGER, INTENT(in) :: return_status

    IF (return_status /= nf_noerr) THEN
      CALL finish('mo_io_grid netCDF error', nf_strerror(return_status))
    ENDIF

  END SUBROUTINE nf
  !-------------------------------------------------------------------------

END MODULE mo_grid_subset

