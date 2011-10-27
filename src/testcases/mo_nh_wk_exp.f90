!>
!!  Subroutine to initialized the Weisman Klemp test case 
!!   for the NH-Core in limited area mode
!!
!!
!! @par Revision History
!! - first version by P. Ripodas , DWD, (2011-08)
!!
!! @par Literature
!! -M. L. Weisman and J. B. Klemp, 1982
!!  The Dependence of Numerically Simulated Convective Storms on 
!!  Vertical Wind Shear and Buoyancy.
!!  Monthly Weather Review, 110,504-520
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
MODULE mo_nh_wk_exp
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
   USE mo_physical_constants,  ONLY: rd_o_cpd, p0ref, grav, tmelt,  &
                                   & cvd_o_rd, re, omega,     cpd ,     &
                                     rvd_m_o => vtmpc1 ,                &
                                     rdv,                        &
                                     rd,                         &
                                     cp_d => cpd
   USE mo_model_domain,        ONLY: t_patch
   USE mo_nonhydro_state,      ONLY: t_nh_prog, t_nh_diag, t_nh_metrics
   USE mo_run_config,          ONLY: iqv,iqc,iqcond, ntracer
   USE mo_impl_constants,      ONLY: inwp, MAX_CHAR_LENGTH
   USE mo_parallel_config,     ONLY: nproma, p_test_run
   USE mo_satad,               ONLY:  sat_pres_water, &  !! saturation vapor pressure w.r.t. water
            &                         sat_pres_ice,   &  !! saturation vapor pressure w.r.t. ice
            &                         spec_humi          !! Specific humidity
   USE mo_exception,           ONLY: message, finish, message_text
   USE mo_advection_config,    ONLY: advection_config
   USE mo_interpolation,       ONLY: t_int_state, cells2edges_scalar, edges2cells_scalar
   USE mo_loopindices,         ONLY: get_indices_e
   USE mo_nh_diagnose_pres_temp,ONLY: diagnose_pres_temp
   USE mo_extpar_config,        ONLY: itopo
   USE mo_sync,                 ONLY: global_sum_array, sync_patch_array, SYNC_C, SYNC_E
   USE mo_nh_init_utils,        ONLY: init_w, hydro_adjust, convert_thdvars, virtual_temp

   IMPLICIT NONE

   PUBLIC  :: init_nh_topo_wk, init_nh_env_wk
  
   PRIVATE
! !DEFINED PARAMETERS for Weisman Klemp:
   REAL(wp), PARAMETER ::hmin_wk          = 0.0_wp         ! [m], base height of the profile
   REAL(wp), PARAMETER ::p_base_wk        = 1e5                ! [Pa], pressure at height hmin_wk

   REAL(wp), PARAMETER ::h_tropo_wk       = 12000._wp      ! [m], height of the tropopause
   REAL(wp), PARAMETER :: theta_0_wk      = 300._wp        ! [K], pot. temp. at z=0
   REAL(wp), PARAMETER ::theta_tropo_wk   = 343._wp        ! [K], pot. temp. at z=h_tropo_wk
   REAL(wp), PARAMETER ::expo_theta_wk    = 1.25_wp
   REAL(wp), PARAMETER ::expo_relhum_wk   = 1.25_wp
   REAL(wp), PARAMETER ::t_tropo_wk       = 213.0_wp


  ! Data for the initial profile of relative humidity:
   REAL(wp), PARAMETER :: rh_min_wk   = 0.25_wp        ! [%] ! rel. hum. above the tropopause level
   REAL(wp), PARAMETER :: rh_max_wk   = 1.0_wp         ! [%]  

   INTEGER,  PARAMETER :: niter=20
   REAL(wp), PUBLIC :: qv_max_wk   != 0.012_wp  - 0.016_wp
  ! Data for the wind profile
   REAL(wp), PUBLIC :: u_infty_wk  != 0 - 45  ms-1
   REAL(wp), PARAMETER :: href_wk  = 3000._wp ! Scaling height (height of 70 % windspeed) 
                                    ! for the Weisman-Klemp wind profile [m]
  
!--------------------------------------------------------------------

   CONTAINS
!-------------------------------------------------------------------------
!
  !>
  !! Initialization of topograpphy for the Weisman Klemp test case
  !! The topography is set to 0, but the possibility  to have 
  !!  a different topography is open
  !!
  !! @par Revision History
  !!
  !!
  SUBROUTINE init_nh_topo_wk( ptr_patch, topo_c, topo_v, nblks_c, npromz_c,      &
                             &  nblks_v, npromz_v )

    TYPE(t_patch), TARGET,INTENT(INOUT) :: &  !< patch on which computation is performed
      &  ptr_patch

    INTEGER,  INTENT (IN) ::  nblks_c, nblks_v, npromz_c, npromz_v
    REAL(wp), INTENT(INOUT) :: topo_c    (nproma,nblks_c)
    REAL(wp), INTENT(INOUT) :: topo_v    (nproma,nblks_v)

    ! local variables

  INTEGER        :: jc, jv, jb, nlen
!--------------------------------------------------------------------

      DO jb = 1, nblks_c
        IF (jb /= nblks_c) THEN
          nlen = nproma
        ELSE
          nlen =  npromz_c
        ENDIF
        DO jc = 1, nlen
          IF ( itopo==0 ) THEN
           topo_c(jc,jb) = 0._wp
          END IF
        ENDDO
      ENDDO

      DO jb = 1, nblks_v
        IF (jb /=  nblks_v) THEN
          nlen = nproma
        ELSE
          nlen =  npromz_v
        ENDIF
        DO jv = 1, nlen
          IF ( itopo==0 ) THEN
           topo_v(jv,jb) =0._wp
          END IF
        ENDDO
      ENDDO

  END SUBROUTINE init_nh_topo_wk
!-------------------------------------------------------------------------
!
  !>
  !! Initialization of prognostic state vector for the Weisman Klemp  test case 
  !! 
  !!
  !! @par Revision History
  !!
  !!
  SUBROUTINE init_nh_env_wk( ptr_patch, ptr_nh_prog, ptr_nh_diag, &
    &                                p_metrics, p_int, l_hydro_adjust )

    TYPE(t_patch), TARGET,INTENT(INOUT) :: &  !< patch on which computation is performed
      &  ptr_patch

    TYPE(t_nh_prog), INTENT(INOUT)      :: &  !< prognostic state vector
      &  ptr_nh_prog

    TYPE(t_nh_diag), INTENT(INOUT)      :: &  !< diagnostic state vector
      &  ptr_nh_diag


    TYPE(t_nh_metrics), INTENT(IN)      :: p_metrics !< NH metrics state
    TYPE(t_int_state), INTENT(IN)       :: p_int
    LOGICAL, INTENT(IN)                 :: l_hydro_adjust !if .TRUE. hydrostatically balanced 
                                                         ! initial condition

    INTEGER        ::  jc, jb, jk, je,   &
                    nlen, nblks_e, npromz_e,  nblks_c, npromz_c
    INTEGER        :: iter
    INTEGER        :: i_startidx, i_endidx, i_startblk
    INTEGER        :: nlev        !< number of full levels
    REAL(wp)       :: t_tropo, z_klev, zrdm, zcpm, z_relhum, zesat,     &
                      zsqv, z_u

    REAL(wp), ALLOCATABLE :: theta(:,:,:)
    REAL(wp), ALLOCATABLE :: relhum(:,:,:)

    LOGICAL :: lcond 

!--------------------------------------------------------------------
!

    lcond=.FALSE.
    ! number of vertical levels
    nlev   = ptr_patch%nlev
    ALLOCATE (theta(nproma,nlev,ptr_patch%nblks_c), &
              relhum(nproma,nlev,ptr_patch%nblks_c) )

   nblks_c   = ptr_patch%nblks_int_c
   npromz_c  = ptr_patch%npromz_int_c
   nblks_e   = ptr_patch%nblks_int_e
   npromz_e  = ptr_patch%npromz_int_e


! As first aproximation

   t_tropo=t_tropo_wk

! first analytic expresions for theta and relative humidity
    DO jb = 1, nblks_c
      IF (jb /= nblks_c) THEN
         nlen = nproma
      ELSE
         nlen = npromz_c
      ENDIF

      DO jk = nlev, 1, -1
            DO jc = 1, nlen
              z_klev = 0.5_wp * ( p_metrics%z_ifc(jc,jk,jb) + &
                       p_metrics%z_ifc(jc,jk+1,jb) )
            IF ( z_klev <= h_tropo_wk) THEN
              theta(jc,jk,jb) = theta_0_wk + &
                  (theta_tropo_wk - theta_0_wk)*((z_klev-hmin_wk)/   &
                   (h_tropo_wk-hmin_wk))**(expo_theta_wk)
              relhum(jc,jk,jb)  = rh_max_wk - (rh_max_wk - rh_min_wk)*      &
                  ((z_klev-hmin_wk)/(h_tropo_wk-hmin_wk))**(expo_relhum_wk)
            ELSE
             theta(jc,jk,jb) = theta_tropo_wk*EXP(grav/(cpd*t_tropo)*(z_klev-h_tropo_wk))
             relhum(jc,jk,jb)  = rh_min_wk
            ENDIF  
       !==========================================================================
       !.. First: Integrate atmospheric pressure assuming dry air. 
       !   This will be the starting point for
       !   a fixpoint iteration to include also moisture and condensation 
       !   in supersaturated voxels.
       !==========================================================================
   
       !   The temperature profile is not known a priori,
       !   because we would need the pressure to compute it from the theta-profile.
       !   As a starting point, the temperature will simply be that of the ICAO 
       !   polytrope atmosphere with a constant temperature lapse rate of 0.0065 K/m 
       !   up to 11 km height and a base temperature of 293.16 K. Then, the
       !   pressure estimate for the dry base will be that of the ICAO standard
       !   atmosphere. The below iteration will compute the correct temperature-
       !   and pressure profile afterwards.             
            ptr_nh_diag%temp(jc,jk,jb)  = 293.16_wp - 0.0065_wp * MIN(z_klev, 11000.0_wp)
            ptr_nh_prog%theta_v(jc,jk,jb)= theta(jc,jk,jb)
            ptr_nh_prog%exner(jc,jk,jb) = ptr_nh_diag%temp(jc,jk,jb)/theta(jc,jk,jb) 
            ptr_nh_prog%rho(jc,jk,jb)   = ptr_nh_prog%exner(jc,jk,jb)**cvd_o_rd    &
                                          *p0ref/rd/theta(jc,jk,jb)
            ptr_nh_prog%tracer(jc,jk,jb,:)= 0.0_wp 
            ENDDO !jc
      ENDDO !jk     
     ENDDO !jb
        ! the first guess for the pressure is calculated now
        CALL diagnose_pres_temp ( p_metrics, ptr_nh_prog,     &
            &                     ptr_nh_prog, ptr_nh_diag,   &
            &                     ptr_patch,                  &
            &                     opt_calc_pres=.TRUE.        )


    !==========================================================================
    ! Now humidity and clouds come into play:
    !
    ! The total density depends on qv, qc and therefore also on the pressure (through qv),
    ! but the pressure itself depends on qv via the hydrostatic approximation --> 
    ! iterative solution of this implicit equation for piter necessary!
    ! For this, each layer is again assumed to be a polytrope layer,
    ! this time with a constant virtual temperature gradient:
    !
    ! For the sake of reproducible results, we do a fixed number of iterations
    ! instead of iterating until convergence to a certain accuracy.
    !==========================================================================

    DO iter = 1, niter

    DO jb = 1, nblks_c
      IF (jb /= nblks_c) THEN
         nlen = nproma
      ELSE
         nlen = npromz_c
      ENDIF

      DO jk = nlev, 1, -1
            DO jc = 1, nlen
          ! Limit relhum to its maximum possible value p / E(T):

             zesat    = sat_pres_water(ptr_nh_diag%temp(jc,jk,jb))
             z_relhum = MIN(ptr_nh_diag%pres(jc,jk,jb) /zesat , relhum(jc,jk,jb)    )

             ptr_nh_prog%tracer(jc,jk,jb,iqv)=     &
               qv_Tprelhum(ptr_nh_diag%pres(jc,jk,jb),ptr_nh_diag%temp(jc,jk,jb), &
                   z_relhum, ptr_nh_prog%tracer(jc,jk,jb,iqc))

          IF (lcond .AND. z_relhum > 1.0_wp) THEN
            ! condensation is allowed and physically can happen
            ! and relhum > 1.0, so convert qv -> qc to limit relhum to 1.0:
            zsqv = spec_humi(zesat, ptr_nh_diag%pres(jc,jk,jb))
            ptr_nh_prog%tracer(jc,jk,jb,iqv) = MIN ( zsqv       ,                 &
                             MIN(ptr_nh_prog%tracer(jc,jk,jb,iqv), qv_max_wk) )
            ptr_nh_prog%tracer(jc,jk,jb,iqc) = MAX ( 0.0_wp ,                     &
                             MIN(ptr_nh_prog%tracer(jc,jk,jb,iqv), qv_max_wk) - zsqv )
          ELSE
            ! condensation is not allowed or cannot happen physically at 
            ! that pressure and temperature, so just impose the limit qv_max_wk:
            ptr_nh_prog%tracer(jc,jk,jb,iqv) =                                    &
                                MIN(ptr_nh_prog%tracer(jc,jk,jb,iqv), qv_max_wk)
            ptr_nh_prog%tracer(jc,jk,jb,iqc) = 0.0_wp
          END IF
            ! "moist" r_d:
            zrdm     = rd_moist(ptr_nh_prog%tracer(jc,jk,jb,iqv),                 &
                                ptr_nh_prog%tracer(jc,jk,jb,iqc))
            ! "moist" cp:
!!$            zcpm = cp_moist(qv(i,j,1),qc(i,j,1),0.0_ireals)
            ! COSMO-approximation of "moist" cp:
            zcpm     = cp_moist_cosmo(ptr_nh_prog%tracer(jc,jk,jb,iqv),           &
                                     ptr_nh_prog%tracer(jc,jk,jb,iqc),0.0_wp)
            !recalculate temperature
            ptr_nh_diag%temp(jc,jk,jb) = theta(jc,jk,jb)*                         &
                                         (ptr_nh_diag%pres(jc,jk,jb)/p0ref)**(zrdm/zcpm)
            ! recalculate exner
            ptr_nh_prog%exner(jc,jk,jb) = ptr_nh_diag%temp(jc,jk,jb)/theta(jc,jk,jb)
            ENDDO !jc
      ENDDO !jk     
     ENDDO !jb

        ! recalculate pressure 
        CALL diagnose_pres_temp ( p_metrics, ptr_nh_prog,     &
            &                     ptr_nh_prog, ptr_nh_diag,   &
            &                     ptr_patch,                  &
            &                     opt_calc_pres=.TRUE.        )

    END DO ! niter
    ! after enough iteretaions we have exner and the tracers, we already had theta

    ! use subroutine virtual_temp to calculate theta_v
    CALL virtual_temp ( ptr_patch, theta, ptr_nh_prog%tracer(:,:,:,iqv),             &
                     & temp_v= ptr_nh_prog%theta_v)
     DO jb = 1, nblks_c
      IF (jb /= nblks_c) THEN
         nlen = nproma
      ELSE
         nlen = npromz_c
      ENDIF

      DO jk = nlev, 1, -1
            DO jc = 1, nlen  
             ptr_nh_prog%rho(jc,jk,jb)  = ptr_nh_prog%exner(jc,jk,jb)**cvd_o_rd*p0ref &
                             /rd/ptr_nh_prog%theta_v(jc,jk,jb)
             ptr_nh_prog%rhotheta_v(jc,jk,jb)  = ptr_nh_prog%rho(jc,jk,jb) *          &
                                        ptr_nh_prog%theta_v(jc,jk,jb)
            ENDDO !jc
      ENDDO !jk     
     ENDDO !jb


! initialized horizontal velocities
    i_startblk = ptr_patch%edges%start_blk(2,1)
    ! horizontal normal components of the velocity
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk,je,z_u, z_klev)
    DO jb = i_startblk, nblks_e

      CALL get_indices_e(ptr_patch, jb, i_startblk, nblks_e, &
                         i_startidx, i_endidx, 2)

        DO jk = 1, nlev
          DO je = i_startidx, i_endidx
            z_klev = p_metrics%z_mc_e(je,jk,jb) 
            z_u = u_infty_wk * TANH(z_klev-hmin_wk/href_wk-hmin_wk)  !v component is zero
            ptr_nh_prog%vn(je,jk,jb) = &
             z_u * ptr_patch%edges%primal_normal(je,jb)%v1
          ENDDO !je
        ENDDO !jk
    ENDDO  !jb
!$OMP END DO
!$OMP END PARALLEL

! initialized vertical velocity

   CALL init_w(ptr_patch, p_int, ptr_nh_prog%vn, p_metrics%z_ifc, ptr_nh_prog%w)
   CALL sync_patch_array(SYNC_C, ptr_patch, ptr_nh_prog%w)


  IF (l_hydro_adjust) THEN

   CALL hydro_adjust ( ptr_patch, p_metrics, ptr_nh_prog%rho,     &
                     & ptr_nh_prog%exner, ptr_nh_prog%theta_v,    &
                     & ptr_nh_prog%rhotheta_v  )

  END IF

       
  END SUBROUTINE init_nh_env_wk
!--------------------------------------------------------------------
! Function taken from COSMO

  ! Specific humidity as function of T, p, and relhum 
  !   (and qcrs = sum of all hydrometeor contents):
  ! NOTE: on input, relhum has to be smaller than p / E(T)!
  REAL(wp) FUNCTION qv_Tprelhum(p, temp, relhum, qcrs)
    IMPLICIT NONE
    REAL(wp), INTENT(IN) :: p, temp, relhum, qcrs
    REAL(wp) :: zesat_w, coeff
    !zesat_w  =b1 * EXP( b2w * (temp-b3)/(temp-b4w) )
    zesat_w  = sat_pres_water(temp)
    coeff = relhum * zesat_w * rdv / p
    qv_Tprelhum = coeff * (1.0 + qcrs) / (1.0 - rvd_m_o*coeff)
  END FUNCTION qv_Tprelhum
!--------------------------------------------------------------------
  ! Gas constant of moist air containing hydrometeors (qcrs is the sum of
  ! the specific hydrometeor contents):
  REAL(wp) FUNCTION rd_moist(qv,qcrs)
    IMPLICIT NONE
    REAL(wp), INTENT(in) :: qv, qcrs

    rd_moist = rd * (1.0 + rvd_m_o*qv - qcrs)

  END FUNCTION rd_moist
!--------------------------------------------------------------------
  ! COSMO-APPROXIMATION: CP IS APPROXIMATED TO BE THAT OF DRY AIR
  REAL(wp) FUNCTION cp_moist_cosmo(qv,ql,qi)
    IMPLICIT NONE
    REAL(wp), INTENT(in) :: qv, ql, qi
    REAL(wp)             :: qd

    ! Dry air content:
    qd = 1.0 - qv - ql -qi

    ! cp:
    cp_moist_cosmo = qd*cp_d + qv*cp_d + ql*cp_d + qi*cp_d

  END FUNCTION cp_moist_cosmo
!--------------------------------------------------------------------
  END MODULE mo_nh_wk_exp
