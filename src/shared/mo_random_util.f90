!>
!! This module contains the subroutine that Adds random perturbation
!! to the normal wind for the Non-Hyd core
!!
!!
!! @author Pilar Ripodas, DWD
!!  based in mo_ha_prog_util from Hui Wan for the Hydrostatic core
!!
!!
!! @par Revision History
!! First version for non-hydrostatic core by Pilar Ripodas, DWD (2010-09-14)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_random_util

  USE mo_kind,                ONLY: wp
  USE mo_impl_constants,      ONLY: SUCCESS, MAX_CHAR_LENGTH, ON_CELLS, ON_EDGES, ON_VERTICES
  USE mo_exception,           ONLY: message, finish, message_text
  USE mo_grid_subset,         ONLY: t_subset_range, get_index_range
 ! USE mo_mpi,                 ONLY: get_my_global_mpi_id, global_mpi_barrier

  IMPLICIT NONE
  PRIVATE
 
  PUBLIC :: add_random_noise_global

CONTAINS

  !-----------------------------------------------------------------------
  !>
  !! Add random perturbation using global index
  !! for seeds, thus guaranteeing reproducability
  !!
  !! @par Revision History
  !!  Created by Leonidas Linardakis, MPIM, 2014
  !!  Based in  hydro_state_prog_add_random (by Hui)
  !!    and the Non-Hyd. version by Pilar Ripodas, DWD, 2010-09
  !!
  SUBROUTINE add_random_noise_global( &
    & in_subset,                 & ! in
    & in_var,                    & ! inout
    & start_level, end_level,    & ! in
    & noise_scale,               & ! in
    & global_vertical_seed,      & ! input, OPTIONAL
    & add_seed)                    ! input, OPTIONAL

    TYPE(t_subset_range)    :: in_subset
    REAL(wp), INTENT(INOUT) :: in_var(:,:,:) ! variable to perturb
    INTEGER, INTENT(IN)     :: start_level, end_level
    REAL(wp), INTENT(IN)    :: noise_scale ! magnitude of the perturbation
    INTEGER, INTENT(IN), TARGET, OPTIONAL  :: global_vertical_seed(:)
    INTEGER, INTENT(IN), OPTIONAL    :: add_seed ! add this number to the fisrt seed element to obtain the rest of them


    ! LOCAL VARIABLES
    ! INTEGER :: DateTimeArray(8)    ! Holds the date and time
    INTEGER :: status
    INTEGER :: level, block, idx, start_idx, end_idx

    INTEGER :: seed_size, js, seed_trigger
    INTEGER, ALLOCATABLE :: seed_array(:)

    REAL(wp), ALLOCATABLE :: noise_1D(:), noise_3D(:,:,:)
    INTEGER, POINTER  :: vertical_seed(:)
    INTEGER    :: add_thisSeed

    CHARACTER(len=*), PARAMETER ::  &
      &  method_name = 'mo_random_util:add_random_noise_global'

    !-----
    CALL message(TRIM(method_name),'=========== generating random noise based on global index =============')

    !-----------------------------------------------------------
    ! 1. prepare memory for the seed
    !-----------------------------------------------------------

    CALL RANDOM_SEED(SIZE=seed_size)
    WRITE(message_text,*) 'The size of the intrinsic seed is', seed_size
    CALL message(method_name,TRIM(message_text))

    ALLOCATE( seed_array(seed_size), STAT=status)
    IF(status/=SUCCESS)THEN
      CALL finish(method_name,'allocation of seed_array failed')
    ENDIF
 
    IF (PRESENT(global_vertical_seed)) THEN
      vertical_seed => global_vertical_seed
    ELSE
      SELECT CASE (in_subset%entity_location)
      CASE (ON_CELLS )
        vertical_seed => in_subset%patch%cells%decomp_info%glb_index
      CASE (ON_EDGES )
        vertical_seed => in_subset%patch%edges%decomp_info%glb_index
      CASE (ON_VERTICES )
        vertical_seed => in_subset%patch%verts%decomp_info%glb_index
      CASE default
        CALL finish(method_name, "Unknown in_subset%entity_location")
      END SELECT
    END IF
    add_thisSeed = 7
    IF (PRESENT(add_seed)) add_thisSeed = add_seed


    ! use the global index as a seed in order to ensure p_test works
    ! this implies a random sequence for each column
    ALLOCATE(noise_1D(start_level:end_level), STAT=status)

    DO block = in_subset%start_block, in_subset%end_block
      CALL get_index_range( in_subset, block, start_idx, end_idx)
      DO idx = start_idx, end_idx
        seed_trigger = ABS(vertical_seed( (block-1) * in_subset%block_size + idx) )
        DO js = 0, seed_size-1
           seed_array(js+1) = seed_trigger + js * add_thisSeed
        ENDDO
        CALL RANDOM_SEED( PUT=seed_array )
        CALL RANDOM_NUMBER( noise_1D(start_level:end_level))

!          WRITE(0,*) get_my_global_mpi_id(), ":", vertical_seed( (block-1) * in_subset%block_size + idx), &
!            " seeds:", seed_array(:), " noise=", noise_1D(start_level:end_level)

        DO level = start_level, end_level
          in_var(idx,level,block) = in_var(idx,level,block) + (noise_1D(level) - 0.5_wp) * noise_scale
        ENDDO

      ENDDO ! idx = start_idx, end_idx

    ENDDO ! block = in_subset%start_block, in_subset%end_block

    DEALLOCATE( noise_1D )

 !   CALL global_mpi_barrier()
 !   CALL finish(method_name, " global_vertical_seed")

  END SUBROUTINE add_random_noise_global

END MODULE mo_random_util
