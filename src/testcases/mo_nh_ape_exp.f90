!>
!!  Aqua Planet experiments for the NH-Core
!!
!!
!! @par Revision History
!! - first version by P. Ripodas , DWD, (2010-09)
!!
!! @par Literature
!! - Neale, R. B. and Hoskins B. J.(2001): "A standard test for AGCMs
!!   including their physical parametrizations:I: The Proposal", Atmospheric
!!   Science Letters (2001)
!! - Neale, R. B. and Hoskins B. J.(2001): "A standard test for AGCMs
!!   including their physical parametrizations: II: Results for The Met
!!   Office Model", Atmospheric Science Letters (2001)
!! -
!!
!! @par Copyright
!! 2002-2008 by DWD and MPI-M
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
MODULE mo_nh_ape_exp
!-------------------------------------------------------------------------
!
!    ProTeX FORTRAN source: Style 2
!    modified for ICON project, DWD/MPI-M 2008
!
!-------------------------------------------------------------------------
!
!
!

USE mo_kind,                ONLY: wp
USE mo_physical_constants,  ONLY: rd, rd_o_cpd, p0ref, grav, tmelt,  &
                                & cvd_o_rd, vtmpc1
USE mo_model_domain,        ONLY: t_patch
USE mo_ext_data,            ONLY: t_external_data
USE mo_nonhydro_state,      ONLY: t_nh_prog, t_nh_diag, t_nh_metrics
USE mo_run_nml,             ONLY: nproma, iqv, iqcond,ntracer
USE mo_satad,               ONLY:  sat_pres_water, &  !! saturation vapor pressure w.r.t. water
      &                            sat_pres_ice,   &  !! saturation vapor pressure w.r.t. ice
      &                            spec_humi          !! Specific humidity
USE mo_ape_params,          ONLY: ps_dry_init

IMPLICIT NONE

PRIVATE

CHARACTER(LEN=*), PARAMETER :: version = '$Id$'

REAL(wp), PARAMETER :: zp_ape      = 101325._wp            !< surface pressure
REAL(wp), PARAMETER :: zt_ape      = 300._wp               !< atmospheric temperature
REAL(wp), PARAMETER :: ztmc_ape    = 25.006_wp             !< total moisture content


PUBLIC :: init_nh_state_prog_APE

!--------------------------------------------------------------------

CONTAINS

!-------------------------------------------------------------------------
!
  !>
  !!               Initialization of prognostic state vector.
  !!
  !!
  !! @par Revision History
  !!
  !!
  SUBROUTINE init_nh_state_prog_APE( ptr_patch, ptr_nh_prog, ptr_nh_diag, &
    &                                ptr_ext_data, p_metrics, rh_at_1000hpa, qv_max )

    TYPE(t_patch), TARGET,INTENT(INOUT) :: &  !< patch on which computation is performed
      &  ptr_patch

    TYPE(t_nh_prog), INTENT(INOUT)      :: &  !< prognostic state vector
      &  ptr_nh_prog

    TYPE(t_nh_diag), INTENT(INOUT)      :: &  !< diagnostic state vector
      &  ptr_nh_diag

    TYPE(t_external_data), INTENT(INOUT):: &  !< external data
      &  ptr_ext_data

    TYPE(t_nh_metrics), INTENT(IN)      :: p_metrics !< NH metrics state


    INTEGER               :: nblks_e, nblks_c, npromz_e, npromz_c,  &
                             nlen
    INTEGER               :: nlev             !< number of full levels
    INTEGER               :: jk, jb, jjt, jc  ! loop variables

    REAL(wp)              :: zscale_h             !< initialized variables
    REAL(wp), ALLOCATABLE :: z_sfc(:,:,:)
    REAL(wp)              :: zrhf, zsqv, z_help
    REAL(wp)              :: rh_at_1000hpa, qv_max
    REAL(wp)              :: tot_moist, iter


!--------------------------------------------------------------------
!
    zscale_h = rd*zt_ape/grav

    nblks_c  = ptr_patch%nblks_int_c
    npromz_c = ptr_patch%npromz_int_c
    nblks_e  = ptr_patch%nblks_int_e
    npromz_e = ptr_patch%npromz_int_e

    nlev     = ptr_patch%nlev

    ALLOCATE (z_sfc(nproma, 1, ptr_patch%nblks_c))

    ! topography
    ptr_ext_data%atm%topography_c(:,:) = 0.0_wp

    ! init surface pressure
    ptr_nh_diag%pres_sfc(:,:) = zp_ape

    ! init temperature
    ptr_nh_diag%temp(:,:,:)   = zt_ape


!!$OMP PARALLEL
!!$OMP DO PRIVATE(jb,jk,jc,nlen,jjt,z_sfc,zrhf,z_help,zsqv)
    DO jb = 1, nblks_c
      IF (jb /= nblks_c) THEN
         nlen = nproma
      ELSE
         nlen = npromz_c
      ENDIF

      DO jk = 1, nlev

        ! compute full level pressure
        z_sfc(1:nlen,1,jb) = ptr_ext_data%atm%topography_c(1:nlen,jb)
        ptr_nh_diag%pres(1:nlen,jk,jb) = zp_ape                                      &
          &                  * exp(-(p_metrics%z_mc(1:nlen,jk,jb)-z_sfc(1:nlen,1,jb))/zscale_h)
      END DO
     END DO

    DO jb = 1, nblks_c
      IF (jb /= nblks_c) THEN
         nlen = nproma
      ELSE
         nlen = npromz_c
      ENDIF

      DO jk = 1, nlev

        ! Introduce a linear decreasing RH with pressure like in Hui´s test cases
              DO jjt = 1, iqcond

                  IF(jjt == iqv ) THEN

                      DO jc =1, nlen

                       zrhf =  rh_at_1000hpa - 0.5_wp                      &
                              & +ptr_nh_diag%pres(jc,jk,jb)/200000._wp
                        z_help = sat_pres_water( ptr_nh_diag%temp(jc,jk,jb) )
                        zsqv = spec_humi( z_help, ptr_nh_diag%pres(jc,jk,jb) )

                       IF( ptr_nh_diag%temp(jc,jk,jb) <= tmelt) THEN
                         z_help =  sat_pres_ice( ptr_nh_diag%temp(jc,jk,jb) )
                         zsqv = spec_humi(z_help, ptr_nh_diag%pres(jc,jk,jb) )
                       ENDIF

                       ptr_nh_prog%tracer(jc,jk,jb,jjt) = MIN ( zsqv,  zrhf*zsqv )
                       IF (  ptr_nh_diag%pres(jc,jk,jb) <= 10000._wp) &
                         &  ptr_nh_prog%tracer(jc,jk,jb,jjt)      &
                         &  = MIN ( 5.e-6_wp, ptr_nh_prog%tracer(jc,jk,jb,jjt) )

                       ! Limit QV in the tropics
                       !ptr_nh_prog%tracer(jc,jk,jb,jjt) = &
                       !     &   MIN(qv_max,ptr_nh_prog%tracer(jc,jk,jb,jjt))

                     END DO

                  ELSE !other tracers than water vapor zero at start
                    ptr_nh_prog%tracer(:,jk,jb,jjt) = 0._wp
                  ENDIF ! tracer

              ENDDO ! tracer loop


        ! init exner pressure
        ptr_nh_prog%exner(1:nlen,jk,jb) = (ptr_nh_diag%pres(1:nlen,jk,jb)/p0ref)**rd_o_cpd


        ! init rhotheta_v
        ptr_nh_prog%rhotheta_v(1:nlen,jk,jb) =                                   &
                     &  ptr_nh_prog%exner(jc,jk,jb)**cvd_o_rd*p0ref/rd



      ENDDO !jk
    ENDDO !jb
!!$OMP END DO 

! provisional value for rho
   ptr_nh_prog%rho=ptr_nh_prog%rhotheta_v

  DO iter=1,10
! ! Set the total moisture content to 25.006 Kg/m2
    tot_moist=0.0_wp
!!$OMP DO PRIVATE(jb,jk,nlen,jc,tot_moist,z_help)
    DO jb = 1, nblks_c
      IF (jb /= nblks_c) THEN
         nlen = nproma
      ELSE
         nlen = npromz_c
      ENDIF
      DO jk = 1, nlev
         DO jc = 1, nlen
           z_help = p_metrics%ddqz_z_full(jc,jk,jb) * ptr_nh_prog%rho(jc,jk,jb)
           tot_moist= tot_moist + ptr_nh_prog%tracer(jc,jk,jb,iqv) * z_help
         ENDDO 
      ENDDO !jk
    ENDDO !jb

!!$OMP END DO
     IF (tot_moist .GT. 1.e-25) THEN
       ptr_nh_prog%tracer(:,:,:,iqv) = ptr_nh_prog%tracer(:,:,:,iqv) * ztmc_ape / tot_moist
     END IF
!!$OMP DO PRIVATE(jb,jk,nlen)    
    DO jb = 1, nblks_c
      IF (jb /= nblks_c) THEN
         nlen = nproma
      ELSE
         nlen = npromz_c
      ENDIF

      DO jk = 1, nlev

        ! init virtual potential temperature
        ptr_nh_prog%theta_v(1:nlen,jk,jb) = zt_ape &
              &       *( 1._wp +  vtmpc1            &
              &       *ptr_nh_prog%tracer(1:nlen,jk,jb,iqv) ) &
              &        /ptr_nh_prog%exner(1:nlen,jk,jb)

        ! init density field rho
        ptr_nh_prog%rho(1:nlen,jk,jb) = &
              &         ptr_nh_prog%rhotheta_v(1:nlen,jk,jb) &
              &         /ptr_nh_prog%theta_v(1:nlen,jk,jb)
        !
        ! initialize horizontal velocity field
        !

        ptr_nh_prog%vn(1:nlen,jk,jb) =     0.0_wp

        !
        ! initialize vertical velocity field
        !

        ptr_nh_prog%w(1:nlen,jk,jb) =      0.0_wp


      ENDDO !jk
    ENDDO !jb

!!$OMP END DO 

 END DO
!!$OMP END PARALLEL
    !
    WRITE(0,*) "tot_moist ",tot_moist


  END SUBROUTINE init_nh_state_prog_APE


!--------------------------------------------------------------------
  END MODULE mo_nh_ape_exp
