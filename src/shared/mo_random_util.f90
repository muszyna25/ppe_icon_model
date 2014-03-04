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
!!
MODULE mo_random_util

  USE mo_kind,                ONLY: wp
  USE mo_impl_constants,      ONLY: SUCCESS, MAX_CHAR_LENGTH
  USE mo_exception,           ONLY: message, finish
  USE mo_grid_subset,         ONLY: t_subset_range, get_index_range
 ! USE mo_mpi,                 ONLY: get_my_global_mpi_id, global_mpi_barrier

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER :: version = '$Id$'
 
  PUBLIC :: add_random_noise

CONTAINS

  !-------------
  !>
  !! SUBROUTINE add_random_noise
  !! Add random perturbation
  !! Currenlty uses only given seeds for reproducability
  !!
  !! @par Revision History
  !!  Based in  hydro_state_prog_add_random (by Hui)
  !!  Adapted to the Non-Hyd. core by Pilar Ripodas, DWD, 2010-09
  !!  
  !! 
  !!
  SUBROUTINE add_random_noise(in_subset,                 & ! in
                              in_var,                    & ! inout
                              start_level, end_level,    & ! in
                              noise_scale,               & ! in
                              global_vertical_seed,      & ! input
                              add_seed)                    ! input

    TYPE(t_subset_range)    :: in_subset
    REAL(wp), INTENT(INOUT) :: in_var(:,:,:) ! variable to perturb
    INTEGER, INTENT(IN)     :: start_level, end_level
    REAL(wp), INTENT(IN)    :: noise_scale ! magnitude of the perturbation
    INTEGER, INTENT(IN), OPTIONAL  :: global_vertical_seed(:)
    INTEGER, INTENT(IN)    :: add_seed ! add this number to the fisrt seed element to obtain the rest of them


    ! LOCAL VARIABLES
    ! INTEGER :: DateTimeArray(8)    ! Holds the date and time
    INTEGER :: status
    INTEGER :: level, block, idx, start_idx, end_idx

    INTEGER :: seed_size, js, seed_trigger
    INTEGER, ALLOCATABLE :: seed_array(:)

    REAL(wp), ALLOCATABLE :: noise_1D(:), noise_3D(:,:,:)

    CHARACTER(len=*), PARAMETER ::  &
      &  method_name = 'mo_random_util:add_random_noise'

    CHARACTER(len=MAX_CHAR_LENGTH) :: string
    !-----
    CALL message(TRIM(method_name),'=========== generating random number =============')

    !-----------------------------------------------------------
    ! 1. prepare memory for the seed
    !-----------------------------------------------------------

    CALL RANDOM_SEED(SIZE=seed_size)
    WRITE(string,*) 'The size of the intrinsic seed is', seed_size
    CALL message(method_name,TRIM(string))

    ALLOCATE( seed_array(seed_size), STAT=status)
    IF(status/=SUCCESS)THEN
      CALL finish(method_name,'allocation of seed_array failed')
    ENDIF

    !----------------------------------------------------------
    ! 2. get seed
    !----------------------------------------------------------
    ! First inquire the current date and time.
    ! The output of the DATE_AND_TIME method_name contains 8 elements:
    !   1: year (e.g. 2008)
    !   2: month
    !   3: day
    !   4: time difference with UTC in minutes
    !   5: hour
    !   6: minute
    !   7: seconds
    !   8: milliseconds

    ! CALL DATE_AND_TIME( VALUES=DateTimeArray )

!   for debugging the seed should be a constant to in theory give reproducible
!   results.  
!    seed_trigger =   DateTimeArray(2)*100000000 + DateTimeArray(3)*1000000 &
!                 & + DateTimeArray(5)*10000     + DateTimeArray(6)*100     &
!                 & + DateTimeArray(7)
!
 
    IF (PRESENT(global_vertical_seed)) THEN
      ! use the global index as a seed in order to ensure p_test works
      ! this implies a random sequence for each column
      ALLOCATE(noise_1D(start_level:end_level), STAT=status)

      DO block = in_subset%start_block, in_subset%end_block
        CALL get_index_range( in_subset, block, start_idx, end_idx)
        DO idx = start_idx, end_idx
          seed_trigger = ABS(global_vertical_seed( (block-1) * in_subset%block_size + idx) )
          DO js = 0, seed_size-1
             seed_array(js+1) = seed_trigger + js * add_seed
          ENDDO
          CALL RANDOM_SEED( PUT=seed_array )
          CALL RANDOM_NUMBER( noise_1D(start_level:end_level))

!          WRITE(0,*) get_my_global_mpi_id(), ":", global_vertical_seed( (block-1) * in_subset%block_size + idx), &
!            " seeds:", seed_array(:), " noise=", noise_1D(start_level:end_level)

          DO level = start_level, end_level
            in_var(idx,level,block) = in_var(idx,level,block) + (noise_1D(level) - 0.5_wp) * noise_scale
          ENDDO

        ENDDO ! idx = start_idx, end_idx

      ENDDO ! block = in_subset%start_block, in_subset%end_block

      DEALLOCATE( noise_1D )

    ELSE
      CALL finish(method_name, "Not implemented without global_vertical_seed")
    ENDIF ! (PRESENT(global_vertical_seed)) THEN

 !   CALL global_mpi_barrier()
 !   CALL finish(method_name, " global_vertical_seed")

  END SUBROUTINE add_random_noise

END MODULE mo_random_util
