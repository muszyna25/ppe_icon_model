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
  USE mo_nonhydro_state,      ONLY: t_nh_prog
  USE mo_exception,           ONLY: message, finish

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

  PUBLIC :: nh_prog_add_random, forcing_straka

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
  
  !-----------------------------------------------------------------------------
  !>
  !! forcing_straka
  !!
  !! Computes the tendency from an (artificial) forcing term.
  !! Currently that contains the Straka test case.
  !!
  !! @par Revision History
  !! Initial release by Almut Gassmann, MPI-M (2009-10.16)
  !!
  SUBROUTINE  forcing_straka(p_nh_prog,p_patch,p_int_state,p_metrics,p_nh_diag)

    TYPE(t_patch),TARGET, INTENT(in) :: p_patch    !< single patch
    TYPE(t_int_state),INTENT(in)  :: p_int_state!< single interpolation state
    TYPE(t_nh_metrics),INTENT(in) :: p_metrics  !< single metrics state
    TYPE(t_nh_prog), INTENT(in)   :: p_nh_prog  !< single nh prognostic state
    TYPE(t_nh_diag), INTENT(inout):: p_nh_diag  !< single nh diagnostic state

    REAL(wp), PARAMETER :: z_ny=75.0_wp !< Straka test diffusion coefficient
    REAL(wp) :: z_lapl_theta(nproma,p_patch%nlev  ,p_patch%nblks_c),  &
                z_lapl_w    (nproma,p_patch%nlevp1,p_patch%nblks_c),  &
                z_lapl_vn   (nproma,p_patch%nlev  ,p_patch%nblks_e)
    INTEGER :: jb, jk, nlen, nblks_e, npromz_e, nblks_c, npromz_c, &
               kp1, km1
    INTEGER :: nlev, nlevp1         !< number of full and half levels
    !--------------------------------------------------------------------------

    nblks_c   = p_patch%nblks_int_c
    npromz_c  = p_patch%npromz_int_c
    nblks_e   = p_patch%nblks_int_e
    npromz_e  = p_patch%npromz_int_e

    ! number of vertical levels
    nlev   = p_patch%nlev
    nlevp1 = p_patch%nlevp1

    ! Horizontal part
    CALL nabla2_scalar(p_nh_prog%theta_v, p_patch, p_int_state, z_lapl_theta)
    CALL nabla2_scalar(p_nh_prog%w, p_patch, p_int_state, z_lapl_w, 2, nlev)
    CALL nabla2_vec(p_nh_prog%vn, p_patch, p_int_state, z_lapl_vn)

    ! Vertical part
!$OMP PARALLEL
!$OMP DO PRIVATE(jb, nlen, jk,kp1,km1)
    DO jb = 1, nblks_e
      IF (jb /= nblks_e) THEN
        nlen = nproma
      ELSE
        nlen = npromz_e
      ENDIF
      DO jk = 1, nlev
        kp1=MIN(jk+1,nlev)
        km1=MAX(jk-1,1)
        p_nh_diag%ddt_vn(1:nlen,jk,jb) = z_ny*(z_lapl_vn(1:nlen,jk,jb)+ &
          (p_nh_prog%vn(1:nlen,km1,jb)-2.0_wp*p_nh_prog%vn(1:nlen,jk,jb)&
           +p_nh_prog%vn(1:nlen,kp1,jb))&
           /(p_metrics%ddqz_z_full_e(1:nlen,jk,jb)**2))
      ENDDO
    ENDDO
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL
!$OMP DO PRIVATE(jb, nlen, jk,kp1,km1)
    DO jb = 1, nblks_c
      IF (jb /= nblks_c) THEN
        nlen = nproma
      ELSE
        nlen = npromz_c
      ENDIF
      DO jk = 1, nlev
        kp1=MIN(jk+1,nlev)
        km1=MAX(jk-1,1)
        p_nh_diag%ddt_exner(1:nlen,jk,jb) = z_ny*&
          (z_lapl_theta(1:nlen,jk,jb)+(p_nh_prog%theta_v(1:nlen,km1,jb)&
          -2.0_wp*p_nh_prog%theta_v(1:nlen,jk,jb)&
                 +p_nh_prog%theta_v(1:nlen,kp1,jb))&
          /(p_metrics%ddqz_z_full(1:nlen,jk,jb)**2))&
          /cvd_o_rd*p_nh_prog%exner(1:nlen,jk,jb)&
          /p_nh_prog%theta_v(1:nlen,jk,jb)
      ENDDO
      DO jk = 2, nlev ! bottom and top tendencies do not exist
        kp1=MIN(jk+1,nlevp1)
        km1=MAX(jk-1,1)
        p_nh_diag%ddt_w(1:nlen,jk,jb) = z_ny*&
          (z_lapl_w(1:nlen,jk,jb)+(p_nh_prog%w(1:nlen,km1,jb)&
          -2.0_wp*p_nh_prog%w(1:nlen,jk,jb)+p_nh_prog%w(1:nlen,kp1,jb))&
          /(p_metrics%ddqz_z_half(1:nlen,jk,jb)**2))
      ENDDO
    ENDDO
!$OMP END DO
!$OMP END PARALLEL

  END SUBROUTINE forcing_straka
  !-----------------------------------------------------------------------------

END MODULE mo_nh_prog_util
