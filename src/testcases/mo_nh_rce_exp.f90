!>
!!  Subroutine to initialize a Radiative Convective Equilibrium Exp 
!!
!!
!! @par Revision History
!! - first version by Levi Silvers , MPIM, (2013-4-24)
!! - revised for global RCE model: start from dry and isothermal atmosphere
!! @par Literature
!! -
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
!!
MODULE mo_nh_rce_exp

  USE mo_kind,                ONLY: wp
  USE mo_model_domain,        ONLY: t_patch
  USE mo_physical_constants,  ONLY: rd, grav, p0ref,rd_o_cpd, vtmpc1, cpd,     &
  &                                 o_m_rdv
  USE mo_nonhydro_types,      ONLY: t_nh_prog, t_nh_diag, t_nh_metrics, t_nh_ref
  USE mo_parallel_config,     ONLY: nproma
  USE mo_nh_testcases_nml,    ONLY: tpe_psfc, tpe_temp, tpe_mu, temp_case,     &
  &                                 temp_fwhm_frac, torus_domain_length,       &
  &                                 th_perturb, lapse_r
  USE mo_nh_wk_exp,           ONLY: bub_ver_width, bubctr_z
  USE mo_impl_constants,      ONLY: max_char_length
  USE mo_run_config,          ONLY: iqv
  USE mo_satad,               ONLY: spec_humi, sat_pres_water

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: init_nh_state_rce_glb

  !DEFINED PARAMETERS (Stevens 2007 JAS):
!! for the initial state, it is unimportant what exact temperature is used
!! this temperature should be consistent with the temperature of the ocean in the RCE case
!! it can be set via the nh_testcase_nml namelist (see the use statement above
!!$  REAL(wp), PARAMETER :: zp0     = 100000._wp !< surface pressure
!!$  REAL(wp), PARAMETER :: zh0     = 0._wp      !< height (m) above which temperature increases
!!$  REAL(wp), PARAMETER :: dtdz    = 0.006_wp   !< lapse rate
!!$  REAL(wp), PARAMETER :: zt0     = 302.15_wp  ! Surface temperature (K)
!  REAL(wp), PARAMETER :: lambda  = 1500._wp   !moist height from Stevens(2007)

!--------------------------------------------------------------------
   CONTAINS
!-------------------------------------------------------------------------
  !>
  !! Initialization of prognostic state vector for the nh RCE test case with a constant T profile
  !!  
  !!
  SUBROUTINE init_nh_state_rce_glb( ptr_patch, ptr_nh_prog,  ptr_nh_ref, ptr_nh_diag,  &
  &                                 ptr_metrics)

    ! INPUT PARAMETERS:
    TYPE(t_patch),TARGET,  INTENT(IN)   :: &  !< patch on which computation is performed
      &  ptr_patch
    TYPE(t_nh_prog),       INTENT(INOUT):: &  !< prognostic state vector
      &  ptr_nh_prog
    TYPE(t_nh_diag),       INTENT(INOUT):: &  !< diagnostic state vector
      &  ptr_nh_diag
    TYPE(t_nh_metrics),    INTENT(IN)   :: &
      &  ptr_metrics                          !< NH metrics state
    TYPE(t_nh_ref),        INTENT(INOUT):: &  !< reference state vector
      &  ptr_nh_ref

    INTEGER  :: jb,jk,jl  ! loop indices
    INTEGER  :: nblks_c,npromz_c
    INTEGER  :: nlen,nlev
    REAL(wp) :: pres, qv, srf_temp, temp, x, sat_pres, tv, fwhm, sigma, z

    ! values for the blocking
    nblks_c  = ptr_patch%nblks_c
    npromz_c = ptr_patch%npromz_c

    ! number of vertical levels
    nlev   = ptr_patch%nlev
  
    ! init surface pressure
    ptr_nh_diag%pres_sfc(:,:) = tpe_psfc
  
    ! Tracers: all zero by default
    ptr_nh_prog%tracer(:,:,:,:) = 0._wp

    ! Parameters for horizontal Gaussian profile
    fwhm = temp_fwhm_frac * torus_domain_length ! full width half max
    sigma = fwhm /  (4._wp * SQRT(4._wp * LOG(2._wp)))! full width half max to sigma
    DO jb = 1, nblks_c

      IF (jb /= nblks_c) THEN
         nlen = nproma
      ELSE
         nlen = npromz_c
      END IF

      DO jk = 1, nlev
        ! start from an isothermal atmosphere with temperature tpe_temp
        ! set pressure and Exner function first
        IF (temp_case == 'blob') THEN
          DO jl = 1, nlen
            x = ptr_patch%cells%cartesian_center(jl, jb)%x(1)
            z = ptr_metrics%z_mc(jl, jk, jb)
            srf_temp = tpe_temp + gaussian(tpe_mu, sigma, 2._wp*th_perturb, x)
            ! Set a 2D patch of moisture at height 'bub_ver_width' [m] above the
            ! centre of the experiment domain. The max moisture is set to
            ! saturation pressure and decreases exponentially towards 0.
            ! Calculate the temperature profile, (instable in the PBL)
            IF ( z <= 2000._wp ) THEN ! we are within the PBL
              temp = srf_temp - lapse_r * 1.5_wp * z
              sat_pres = sat_pres_water(temp) * (2000._wp - z)/2000._wp
            ELSE
              temp = srf_temp - (lapse_r * 2000._wp * 1.5_wp)
              sat_pres = 0._wp
            END IF
              pres     = tpe_psfc * EXP(-grav/rd/temp * ptr_metrics%z_mc(jl,jk,jb))
              sat_pres = sat_pres + (sat_pres_water(tpe_temp + 2._wp*th_perturb)   &
              &           * gaussian(tpe_mu, sigma, 1._wp, x)                  &
              &           * gaussian(bubctr_z, bub_ver_width, 1._wp, z))
              tv                              = temp / (1._wp - (sat_pres / pres) * o_m_rdv)
              qv                              = spec_humi(sat_pres, pres)
              ptr_nh_prog%exner(jl, jk, jb)   = (pres / p0ref)**rd_o_cpd
              ptr_nh_diag%pres(jl, jk, jb)    = pres
              ptr_nh_prog%rho(jl, jk, jb)     = pres / rd / tv
              ptr_nh_prog%theta_v(jl, jk, jb) = tv / ptr_nh_prog%exner(jl, jk, jb) &
                                               * ( 1._wp + vtmpc1 * qv)
              ptr_nh_prog%tracer(jl, jk, jb, iqv) = qv
          END DO !jl
        ELSE ! Constant temperature case
          ptr_nh_diag%pres(1:nlen, jk, jb)    = tpe_psfc * EXP(-grav/rd/tpe_temp*ptr_metrics%z_mc(1:nlen, jk, jb))
          ptr_nh_prog%rho(1:nlen, jk, jb)     = ptr_nh_diag%pres(1:nlen, jk, jb)/rd/tpe_temp
          ptr_nh_prog%exner(1:nlen, jk, jb)   = (ptr_nh_diag%pres(1:nlen, jk, jb)/p0ref)**rd_o_cpd
          ptr_nh_prog%theta_v(1:nlen, jk, jb) = tpe_temp/ptr_nh_prog%exner(1:nlen, jk, jb)
        END IF
      END DO !jk
    END DO ! jb

  !meridional and zonal wind
  ptr_nh_prog%vn = 0._wp
  ptr_nh_ref%vn_ref = ptr_nh_prog%vn

  !vertical wind
  ptr_nh_prog%w = 0._wp
  ptr_nh_ref%w_ref = ptr_nh_prog%w

  END SUBROUTINE init_nh_state_rce_glb

  FUNCTION gaussian(mu, sigma, delta, x)
    !> @brief
    !> calculate Gaussian profile
    !> @param[in] mu
    !> @param[in] sigma
    !> @param[in] delta
    !> @param[in] x
    !
    REAL(wp), INTENT(in) :: mu    !< expectation of the Gaussian
    REAL(wp), INTENT(in) :: sigma !< std-dev. of Gaussian
    REAL(wp), INTENT(in) :: delta !< max. diviation from 0
    REAL(wp), INTENT(in) :: x     !< x value
    REAL(wp) :: gaussian, xx

    xx = ((x - mu) / sigma)
    gaussian = EXP(-.5_wp * xx * xx) * delta
  END FUNCTION gaussian

END MODULE mo_nh_rce_exp

