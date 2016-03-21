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
MODULE mo_nh_prog_util

  USE mo_kind,                ONLY: wp
  USE mo_impl_constants,      ONLY: SUCCESS, MAX_CHAR_LENGTH
  USE mo_model_domain,        ONLY: t_patch
  USE mo_nonhydro_types,      ONLY: t_nh_prog, t_nh_diag
  USE mo_exception,           ONLY: message, finish
  USE mo_parallel_config,     ONLY: nproma

  IMPLICIT NONE
  PRIVATE


  PUBLIC :: nh_prog_add_random
!!$  PUBLIC :: init_nh_state_prog_isoRest  !DR is currently not used at all

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
                                pvar,    & ! inout
                                ctype, pscale, jk_start, jk_end) ! input

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = 'mo_nh_prog_util:nh_prog_add_random'

    TYPE(t_patch)   :: p_patch

    REAL(wp), INTENT(INOUT) :: pvar(:,:,:) ! variable to perturb
    REAL(wp), INTENT(IN)    :: pscale ! magnitude of the perturbation
    CHARACTER(len=*), INTENT(IN)   :: ctype !var type, edge or cell based
    INTEGER,  INTENT(IN)    :: jk_start, jk_end  ! levels between which perturbation is added

    ! LOCAL VARIABLES

    INTEGER :: jk, jb, nlen
    INTEGER :: npromz, nblks

    INTEGER :: DateTimeArray(8)    ! Holds the date and time

    INTEGER :: seed_size, js, seed_trigger, ist

    INTEGER, ALLOCATABLE :: seed_array(:)

    REAL(wp), ALLOCATABLE :: zrand(:,:,:)

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

!    seed_trigger =   DateTimeArray(2)*100000000 + DateTimeArray(3)*1000000 &
!                 & + DateTimeArray(5)*10000     + DateTimeArray(6)*100     &
!                 & + DateTimeArray(7)
    seed_trigger = 704182209

    WRITE(string,*) 'the seed trigger is', seed_trigger
    CALL message('',TRIM(string))

    DO js=1,seed_size
       seed_array(js)=ABS(seed_trigger)+(js-1)
    ENDDO

    CALL message('','Seed generated')

    !-----------------------------------------------------------
    ! 3. generate random numbers and perturb the variable
    !-----------------------------------------------------------

    IF(TRIM(ctype)=="cell")THEN
      nblks  = p_patch%nblks_c
      npromz = p_patch%npromz_c
    ELSEIF(TRIM(ctype)=="edge")THEN
      nblks  = p_patch%nblks_e
      npromz = p_patch%npromz_e
    ELSE
      CALL finish(TRIM(routine),'Wrong ctype for the input variable!')
    END IF

    ALLOCATE(zrand(SIZE(pvar,1),SIZE(pvar,2),SIZE(pvar,3)), STAT=ist)
    IF(ist/=SUCCESS)THEN
      CALL finish(TRIM(routine),'allocation of zrand failed')
    ENDIF

    CALL RANDOM_SEED( PUT=seed_array )
    CALL RANDOM_NUMBER( zrand )

    zrand = (zrand - 0.5_wp)*pscale

    ! add the same random noise to all vertical layers

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,nlen)
    DO jb = 1, nblks
        IF (jb /= nblks) THEN
           nlen = nproma
        ELSE
           nlen = npromz
        ENDIF       
      DO jk=jk_start,jk_end
         pvar(1:nlen,jk,jb) = pvar(1:nlen,jk,jb) + zrand(1:nlen,jk,jb)
      ENDDO
    END DO
!$OMP END DO
!$OMP END PARALLEL

    !-----------------------------------------------------------
    ! 4. clean up
    !-----------------------------------------------------------

    DEALLOCATE( seed_array, STAT=ist)
    IF(ist/=SUCCESS) CALL finish('random number:','deallocation of seed_array failed')
    CALL message('','=========================================')

  END SUBROUTINE nh_prog_add_random
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  !>
  !! Set the prognostic state vector to a resting isothermal state.
  !!
  SUBROUTINE init_nh_state_prog_isoRest( ptemp, ppres_sfc, pt_prog, pt_diag )

    REAL(wp),INTENT(IN) :: ptemp
    REAL(wp),INTENT(IN) :: ppres_sfc

    TYPE(t_nh_prog),    INTENT(INOUT) :: pt_prog      !!the prognostic variables
    TYPE(t_nh_diag),    INTENT(INOUT) :: pt_diag      !!the diagnostic variables

    pt_prog%vn       = 0._wp
    pt_diag%temp     = ptemp
    pt_diag%pres_sfc = ppres_sfc

  END SUBROUTINE init_nh_state_prog_isoRest


END MODULE mo_nh_prog_util
