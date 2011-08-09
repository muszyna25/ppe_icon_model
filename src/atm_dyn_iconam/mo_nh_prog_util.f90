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
MODULE mo_nh_prog_util

  USE mo_kind,                ONLY: wp
  USE mo_impl_constants,      ONLY: SUCCESS, MAX_CHAR_LENGTH
  USE mo_model_domain,        ONLY: t_patch
  USE mo_nonhydro_state,      ONLY:  t_nh_prog
  USE mo_exception,           ONLY: message, finish

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

  PUBLIC :: nh_prog_add_random

CONTAINS
  !-------------
  !>
  !! SUBROUTINE nh_prog_add_random
  !! Add random perturbation to the normal wind.
  !!
  !! @par Revision History
  !!  Based in  hydro_state_prog_add_random (by Hui)
  !!  Adapted to the Non-Hyd. core by Pilar Ripodas, DWD, 2010-09
  !!
  !!
  !!
  SUBROUTINE nh_prog_add_random(p_patch, & ! in
                                         p_nh_prog,  & ! inout
                                         pscale, nproma, nlev) ! input

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = 'mo_nh_prog_util:nh_prog_add_random'

    TYPE(t_patch)            :: p_patch
    TYPE(t_nh_prog) :: p_nh_prog

    REAL(wp), INTENT(IN) :: pscale ! magnitude of the perturbation
    INTEGER,  INTENT(IN) :: nlev   ! number of vertical layers

    ! LOCAL VARIABLES

    INTEGER :: jk
    INTEGER :: nproma, npromz, nblks

    INTEGER :: DateTimeArray(8)    ! Holds the date and time

    INTEGER :: seed_size, js, seed_trigger, ist

    INTEGER, ALLOCATABLE :: seed_array(:)

    REAL(wp) :: zrand(nproma,p_patch%nblks_e)

    CHARACTER(len=MAX_CHAR_LENGTH) :: string
    !-----
    CALL message(TRIM(routine),'=========== generating random number =============')

    !-----------------------------------------------------------
    ! 1. prepare memory for the seed
    !-----------------------------------------------------------

    CALL RANDOM_SEED(SIZE=seed_size)
    WRITE(string,*) 'The size of the intrinsic seed is', seed_size
    CALL message('',TRIM(string))

    ALLOCATE( seed_array(seed_size), STAT=ist)
    IF(ist/=SUCCESS)THEN
      CALL finish('random number:','allocation of seed_array failed')
    ENDIF

    !----------------------------------------------------------
    ! 2. get seed
    !----------------------------------------------------------
    ! First inquire the current date and time.
    ! The output of the DATE_AND_TIME routine contains 8 elements:
    !   1: year (e.g. 2008)
    !   2: month
    !   3: day
    !   4: time difference with UTC in minutes
    !   5: hour
    !   6: minute
    !   7: seconds
    !   8: milliseconds

    CALL DATE_AND_TIME( VALUES=DateTimeArray )

    seed_trigger =   DateTimeArray(2)*100000000 + DateTimeArray(3)*1000000 &
                 & + DateTimeArray(5)*10000     + DateTimeArray(6)*100     &
                 & + DateTimeArray(7)

    WRITE(string,*) 'the seed trigger is', seed_trigger
    CALL message('',TRIM(string))

    DO js=1,seed_size
       seed_array(js)=ABS(seed_trigger)+(js-1)
    ENDDO

    CALL message('','Seed generated')

    !-----------------------------------------------------------
    ! 3. generate random numbers and perturb the normal wind
    !-----------------------------------------------------------

    CALL RANDOM_SEED( PUT=seed_array )
    CALL RANDOM_NUMBER( zrand )

    zrand = (zrand - 0.5_wp)*pscale

    nblks = p_patch%nblks_e
    npromz = p_patch%npromz_e

    ! add the same random noise to all vertical layers
    DO jk=1,nlev
       p_nh_prog%vn(:,jk,1:nblks-1)    = p_nh_prog%vn(:,jk,1:nblks-1)    + zrand(:,1:nblks-1)
       p_nh_prog%vn(1:npromz,jk,nblks) = p_nh_prog%vn(1:npromz,jk,nblks) + zrand(1:npromz,nblks)
    ENDDO

    !-----------------------------------------------------------
    ! 4. clean up
    !-----------------------------------------------------------

    DEALLOCATE( seed_array, STAT=ist)
    IF(ist/=SUCCESS) CALL finish('random number:','deallocation of seed_array failed')
    CALL message('','=========================================')

  END SUBROUTINE nh_prog_add_random
  !-----------------------------------------------------------------------------
  

END MODULE mo_nh_prog_util
