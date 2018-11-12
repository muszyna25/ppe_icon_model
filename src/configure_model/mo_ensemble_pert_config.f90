!>
!! @brief Ensemble perturbations of nwp physics
!!
!! configuration setup for ensemble physics perturbations
!!
!! @author Guenther Zaengl, DWD
!!
!!
!! @par Revision History
!! Initial revision by Guenther Zaengl, DWD (2015-04-23)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_ensemble_pert_config

  USE mo_kind,               ONLY: wp
  USE mo_impl_constants,     ONLY: max_dom, min_rlcell_int
  USE mo_impl_constants_grf, ONLY: grf_bdywidth_c
  USE mo_math_constants,     ONLY: pi2
  USE mo_loopindices,        ONLY: get_indices_c
  USE mo_nwp_tuning_config,  ONLY: tune_gkwake, tune_gkdrag, tune_gfluxlaun, tune_zvz0i,    &
    &                        tune_entrorg, tune_capdcfac_et, tune_box_liq, tune_rhebc_land, &
    &                        tune_rhebc_ocean, tune_rcucov, tune_texc, tune_qexc,           &
    &                        tune_minsnowfrac, tune_rhebc_land_trop, tune_rhebc_ocean_trop, &
    &                        tune_rcucov_trop, tune_gfrcrit, tune_capdcfac_tr,              &
    &                        tune_lowcapefac, limit_negpblcape, tune_rdepths, tune_rprcon,  &
    &                        tune_box_liq_asy
  USE mo_turbdiff_config,    ONLY: turbdiff_config
  USE mo_gribout_config,     ONLY: gribout_config
  USE mo_atm_phy_nwp_config, ONLY: atm_phy_nwp_config
  USE mo_lnd_nwp_config,     ONLY: ntiles_total, ntiles_lnd, ntiles_water, c_soil, cwimax_ml
  USE mo_grid_config,        ONLY: n_dom
  USE mo_parallel_config,    ONLY: nproma
  USE mo_model_domain,       ONLY: t_patch
  USE mo_nwp_phy_types,      ONLY: t_nwp_phy_diag
  USE mo_ext_data_types,     ONLY: t_external_data
  USE mo_extpar_config,      ONLY: nclass_lu
  USE mo_exception,          ONLY: message_text, message
  USE mtime,                 ONLY: datetime, getDayOfYearFromDateTime

  IMPLICIT NONE
  PRIVATE


  PUBLIC :: use_ensemble_pert, configure_ensemble_pert, compute_ensemble_pert
  PUBLIC :: range_gkwake, range_gkdrag, range_gfluxlaun, range_zvz0i, range_entrorg, range_capdcfac_et, &
            range_box_liq, range_tkhmin, range_tkmmin, range_rlam_heat, range_rhebc, range_texc,        &
            range_minsnowfrac, range_z0_lcc, range_rootdp, range_rsmin, range_laimax, range_charnock,   &
            range_tkred_sfc, range_gfrcrit, range_c_soil, range_cwimax_ml, range_capdcfac_tr,           &
            range_lowcapefac, range_negpblcape, stdev_sst_pert, sst_pert_corrfac, range_rdepths,        &
            range_rprcon, range_qexc, range_turlen, range_a_hshr, range_rain_n0fac, range_box_liq_asy,  &
            itype_pert_gen, timedep_pert

  !!--------------------------------------------------------------------------
  !! Basic configuration setup for ensemble perturbations
  !!--------------------------------------------------------------------------


    ! namelist variables
  REAL(wp) :: &                    !< low level wake drag constant
    &  range_gkwake

  REAL(wp) :: &                    !< gravity wave drag constant
    &  range_gkdrag

  REAL(wp) :: &                    !< critical Froude number used for computing blocking layer depth
    &  range_gfrcrit

  REAL(wp) :: &                    !< gravity wave flux emission
    &  range_gfluxlaun

  REAL(wp) :: &                    !< Terminal fall velocity of ice 
    &  range_zvz0i

  REAL(wp) :: &                    !< Tuning factor for intercept parameter of raindrop size distribution 
    &  range_rain_n0fac

  REAL(wp) :: &                    !< Entrainment parameter for deep convection valid at dx=20 km 
    &  range_entrorg

  REAL(wp) :: &                    !< Maximum shallow convection depth 
    &  range_rdepths

  REAL(wp) :: &                    !< Entrainment parameter for deep convection valid at dx=20 km 
    &  range_rprcon

  REAL(wp) :: &                    !< Fraction of CAPE diurnal cycle correction applied in the extratropics
    &  range_capdcfac_et            ! (relevant only if icapdcycl = 3)

  REAL(wp) :: &                    !< Fraction of CAPE diurnal cycle correction applied in the tropics
    &  range_capdcfac_tr            ! (relevant only if icapdcycl = 3)

  REAL(wp) :: &                    !< Tuning factor for reducing the diurnal cycle correction in low-cape situations
    &  range_lowcapefac            ! (relevant only if icapdcycl = 3)

  REAL(wp) :: &                    !< Minimum allowed negative PBL cape in diurnal cycle correction
    &  range_negpblcape            ! (relevant only if icapdcycl = 3)

  REAL(wp) :: &                    !< RH thresholds for evaporation below cloud base
    &  range_rhebc

  REAL(wp) :: &                    !< Excess value for temperature used in test parcel ascent
    &  range_texc

  REAL(wp) :: &                    !< Excess value for mixing ratio used in test parcel ascent
    &  range_qexc

  REAL(wp) :: &                    !< Minimum value to which the snow cover fraction is artificially reduced
    &  range_minsnowfrac           !  in case of melting show (in case of idiag_snowfrac = 20/30/40)

  REAL(wp) :: &                    !< Fraction of surface area available for bare soil evaporation
    &  range_c_soil

  REAL(wp) :: &                    !< Capacity of interception storage (multiplicative perturbation)
    &  range_cwimax_ml

  REAL(wp) :: &                    !< Box width for liquid clouds assumed in the cloud cover scheme
    &  range_box_liq               ! (in case of inwp_cldcover = 1)

  REAL(wp) :: &                    !< Asymmetry factor for sub-grid scale liquid cloud diagnostic
    &  range_box_liq_asy           ! (in case of inwp_cldcover = 1)

  REAL(wp) :: &                    !< Minimum vertical diffusion for heat/moisture 
    &  range_tkhmin

  REAL(wp) :: &                    !< Minimum vertical diffusion for momentum 
    &  range_tkmmin

  REAL(wp) :: &                    !< Perturbation scale (multiplicative) of reduction of minimum diffusion 
    &  range_tkred_sfc             !  coefficients near the surface 

  REAL(wp), ALLOCATABLE :: &        !< Array of random numbers used for computing the above-mentioned perturbation
    &  rnd_tkred_sfc(:)

  REAL(wp) :: &                    !< Laminar transport resistance parameter 
    &  range_rlam_heat

  REAL(wp) :: &                    !< Maximum turbulent mixing length scale 
    &  range_turlen

  REAL(wp) :: &                    !< Scaling factor for extended horizontal shear term in turbulence scheme 
    &  range_a_hshr

  REAL(wp) :: &                    !< Upper and lower bound of wind-speed dependent Charnock parameter 
    &  range_charnock

  REAL(wp) :: &                    !< Roughness length attributed to land-cover class 
    &  range_z0_lcc

  REAL(wp) :: &                    !< Root depth related to land-cover class
    &  range_rootdp

  REAL(wp) :: &                    !< Minimum stomata resistance related to land-cover class
    &  range_rsmin

  REAL(wp) :: &                    !< Maximum leaf area index related to land-cover class
    &  range_laimax

  REAL(wp) :: &                    !< Standard deviation of SST perturbations specified in SST analysis (K)
    &  stdev_sst_pert              !  this switch controls the correction term sst_pert_corrfac, compensating the systematic
                                   !  increase of evaporation related to the SST perturbations

  REAL(wp) :: &                    !< Tuning factor for saturation pressure over oceans, compensating the systematic
    &  sst_pert_corrfac            !  increase of evaporation related to the SST perturbations

  INTEGER :: itype_pert_gen        !< type of ensemble perturbation generation
  INTEGER :: timedep_pert          !< time dependence of ensemble perturbations
  LOGICAL :: use_ensemble_pert     !< main switch

  CONTAINS


  !>
  !! Application of the ensemble perturbation to the config/namelist variables 
  !!
  !! This is done based on random numbers determined by the ensemble member ID
  !!
  !! @par Revision History
  !! Initial revision by Guenther Zaengl, DWD (2015-04-23)
  !!
  SUBROUTINE configure_ensemble_pert(ext_data, mtime_date)

    TYPE(t_external_data), INTENT(INOUT) :: ext_data(:)
    TYPE(datetime),        POINTER       :: mtime_date

    INTEGER, ALLOCATABLE :: rnd_seed(:)
    INTEGER  :: rnd_size, i, jg, ipn, iy, im, id
    REAL(wp) :: rnd_num, rnd_fac, alpha0_sv, z0_lcc, rootdp, rsmin, laimax, tkfac, svp_pert


    IF (use_ensemble_pert .AND. gribout_config(1)%perturbationNumber >= 1) THEN

      CALL RANDOM_SEED(rnd_size)
      ALLOCATE(rnd_seed(rnd_size))

      ! Initialize randum number generator with an integer sequence depending on the ensemble member ID
      ipn = gribout_config(1)%perturbationNumber
      DO i = 1, rnd_size
        rnd_seed(i) = (135+i)*ipn - (21+i**2)*(5+MOD(ipn,10))**2 + 3*i**3
      ENDDO

      IF (timedep_pert == 1) THEN
        iy = mtime_date%date%year
        im = mtime_date%date%month
        id = mtime_date%date%day
        DO i = 1, rnd_size
          rnd_seed(i) =  rnd_seed(i) + (iy-1800)**2 + im**3 + (id+i)**2
        ENDDO
      ENDIF

      CALL RANDOM_SEED(PUT=rnd_seed)
      DO i = 1, 10+ipn
        CALL RANDOM_NUMBER(rnd_num)
      ENDDO

      ! Apply perturbations to physics tuning parameters

      ! SSO tuning
      CALL random_gen(rnd_num)
      ! perturbations for gkwake and gfrcrit must be correlated
      ! (for gfrcrit, a higher value means a thinner blocking layer)
      tune_gkwake(1:max_dom)  = tune_gkwake(1:max_dom)  + 2._wp*(rnd_num-0.5_wp)*range_gkwake
      tune_gfrcrit(1:max_dom) = tune_gfrcrit(1:max_dom) + 2._wp*(rnd_num-0.5_wp)*range_gfrcrit

      CALL random_gen(rnd_num)
      tune_gkdrag(1:max_dom) = tune_gkdrag(1:max_dom) + 2._wp*(rnd_num-0.5_wp)*range_gkdrag

      ! GWD tuning
      CALL random_gen(rnd_num)
      tune_gfluxlaun = tune_gfluxlaun + 2._wp*(rnd_num-0.5_wp)*range_gfluxlaun

      ! grid-scale microphysics
      CALL random_gen(rnd_num)
      ! perturbations for cloud ice sedimentation and convective precip conversion rate must be anticorrelated
      ! because they systematically alter the temperature bias in the middle/upper troposphere
      tune_zvz0i   = tune_zvz0i   + 2._wp*(rnd_num-0.5_wp)*range_zvz0i
      tune_rprcon  = tune_rprcon  - 2._wp*(rnd_num-0.5_wp)*range_rprcon

      CALL random_gen(rnd_num)
      rnd_fac = range_rain_n0fac**(2._wp*(rnd_num-0.5_wp))
      atm_phy_nwp_config(1:max_dom)%rain_n0_factor = atm_phy_nwp_config(1:max_dom)%rain_n0_factor * rnd_fac

      ! convection
      CALL random_gen(rnd_num)
      tune_entrorg = tune_entrorg + 2._wp*(rnd_num-0.5_wp)*range_entrorg

      CALL random_gen(rnd_num)
      tune_rdepths = tune_rdepths + 2._wp*(rnd_num-0.5_wp)*range_rdepths

      CALL random_gen(rnd_num)
      ! Scale factor for CAPE diurnal cycle correction must be non-negative; for the current default
      ! of tune_capdcfac_et=0, the perturbation is zero for half of the ensemble members
      tune_capdcfac_et = MAX(0._wp, tune_capdcfac_et + 2._wp*(rnd_num-0.5_wp)*range_capdcfac_et)
      CALL random_gen(rnd_num)
      tune_capdcfac_tr = MAX(0._wp, tune_capdcfac_tr + 2._wp*(rnd_num-0.5_wp)*range_capdcfac_tr)

      CALL random_gen(rnd_num)
      tune_lowcapefac  = MAX(0._wp, tune_lowcapefac + 2._wp*(rnd_num-0.5_wp)*range_lowcapefac)
      CALL random_gen(rnd_num)
      limit_negpblcape = MIN(0._wp, limit_negpblcape + 2._wp*(rnd_num-0.5_wp)*range_negpblcape)

      CALL random_gen(rnd_num)
      ! Perturbations for RH thresholds for the onset of evaporation below cloud base and the
      ! convective area fraction must be anticorrelated, i.e. the convective area fraction is reduced
      ! when the RH thresholds are increased
      tune_rhebc_land  = tune_rhebc_land  + 2._wp*(rnd_num-0.5_wp)*range_rhebc
      tune_rhebc_ocean = tune_rhebc_ocean + 2._wp*(rnd_num-0.5_wp)*range_rhebc
      tune_rcucov      = tune_rcucov / (1._wp + 15._wp*range_rhebc*(rnd_num-0.5_wp))
      ! corresponding parameters for the tropics
      tune_rhebc_land_trop  = tune_rhebc_land_trop  + 2._wp*(rnd_num-0.5_wp)*range_rhebc
      tune_rhebc_ocean_trop = tune_rhebc_ocean_trop + 2._wp*(rnd_num-0.5_wp)*range_rhebc
      tune_rcucov_trop      = tune_rcucov_trop / (1._wp + 15._wp*range_rhebc*(rnd_num-0.5_wp))

      CALL random_gen(rnd_num)
      tune_texc = tune_texc + 2._wp*(rnd_num-0.5_wp)*range_texc
      CALL random_gen(rnd_num)
      tune_qexc = tune_qexc + 2._wp*(rnd_num-0.5_wp)*range_qexc

      ! cloud cover
      CALL random_gen(rnd_num)
      tune_box_liq = tune_box_liq + 2._wp*(rnd_num-0.5_wp)*range_box_liq

      CALL random_gen(rnd_num)
      tune_box_liq_asy = tune_box_liq_asy + 2._wp*(rnd_num-0.5_wp)*range_box_liq_asy

      ! turbulence
      CALL random_gen(rnd_num)
      tkfac = turbdiff_config(1)%tkhmin_strat/turbdiff_config(1)%tkhmin
      turbdiff_config(1:max_dom)%tkhmin = turbdiff_config(1:max_dom)%tkhmin + 2._wp*(rnd_num-0.5_wp)*range_tkhmin
      turbdiff_config(1:max_dom)%tkhmin_strat = turbdiff_config(1:max_dom)%tkhmin_strat + 2._wp*tkfac*(rnd_num-0.5_wp)*range_tkhmin

      CALL random_gen(rnd_num)
      tkfac = turbdiff_config(1)%tkmmin_strat/turbdiff_config(1)%tkmmin
      turbdiff_config(1:max_dom)%tkmmin = turbdiff_config(1:max_dom)%tkmmin + 2._wp*(rnd_num-0.5_wp)*range_tkmmin
      turbdiff_config(1:max_dom)%tkmmin_strat = turbdiff_config(1:max_dom)%tkmmin_strat + 2._wp*tkfac*(rnd_num-0.5_wp)*range_tkmmin

      CALL random_gen(rnd_num)
      rnd_fac = range_rlam_heat**(2._wp*(rnd_num-0.5_wp))
      ! The product rlam_heat*rat_sea must stay unchanged in order to keep heat and moisture fluxes over the oceans constant
      turbdiff_config(1:max_dom)%rlam_heat = turbdiff_config(1:max_dom)%rlam_heat * rnd_fac
      turbdiff_config(1:max_dom)%rat_sea   = turbdiff_config(1:max_dom)%rat_sea   / rnd_fac

      CALL random_gen(rnd_num)
      turbdiff_config(1:max_dom)%tur_len = turbdiff_config(1:max_dom)%tur_len + 2._wp*(rnd_num-0.5_wp)*range_turlen

      CALL random_gen(rnd_num)
      turbdiff_config(1:max_dom)%a_hshr = turbdiff_config(1:max_dom)%a_hshr + 2._wp*(rnd_num-0.5_wp)*range_a_hshr

      CALL random_gen(rnd_num)
      rnd_fac   = range_charnock**(2._wp*(rnd_num-0.5_wp))
      alpha0_sv = turbdiff_config(1)%alpha0
      !
      ! Upper and lower bound of the variation range of the wind-speed dependent Charnock parameter
      ! are varied inversely in order to avoid bias changes
      turbdiff_config(1:max_dom)%alpha0     = turbdiff_config(1:max_dom)%alpha0     * rnd_fac
      turbdiff_config(1:max_dom)%alpha0_max = turbdiff_config(1:max_dom)%alpha0_max / rnd_fac

      CALL random_gen(rnd_num)
      ! Additional additive perturbation to Charnock parameter
      turbdiff_config(1:max_dom)%alpha0_pert = (rnd_num-0.5_wp)*alpha0_sv*(range_charnock-1._wp)

      ! TERRA
      CALL random_gen(rnd_num)
      tune_minsnowfrac = tune_minsnowfrac + 2._wp*(rnd_num-0.5_wp)*range_minsnowfrac

      CALL random_gen(rnd_num)
      c_soil = c_soil + 2._wp*(rnd_num-0.5_wp)*range_c_soil
      c_soil = MAX(0._wp,MIN(2._wp,c_soil))

      CALL random_gen(rnd_num)
      rnd_fac = range_cwimax_ml**(2._wp*(rnd_num-0.5_wp))
      cwimax_ml = cwimax_ml * rnd_fac


      ! control output
      WRITE(message_text,'(3f8.4,e11.4,2f8.4)') tune_gkwake(1), tune_gkdrag(1), tune_gfrcrit(1), &
        tune_gfluxlaun, tune_zvz0i, atm_phy_nwp_config(1)%rain_n0_factor
      CALL message('Perturbed values, gkwake, gkdrag, gfrcrit, gfluxlaun, zvz0i, rain_n0fac', TRIM(message_text))

      WRITE(message_text,'(2e11.4,f8.1)') tune_entrorg, tune_rprcon, tune_rdepths
      CALL message('Perturbed values, entrorg, rprcon, rdepths', TRIM(message_text))

      WRITE(message_text,'(3f8.4,f8.1)') tune_capdcfac_et, tune_capdcfac_tr, tune_lowcapefac, limit_negpblcape
      CALL message('Perturbed values, capdcfac_et, capdcfac_tr, lowcapefac, negpblcape', TRIM(message_text))

      WRITE(message_text,'(4f8.4,f8.5)') tune_rhebc_land, tune_rhebc_ocean, tune_rcucov, tune_texc, tune_qexc
      CALL message('Perturbed values, rhebc_land, rhebc_ocean, rcucov, texc, qexc', TRIM(message_text))

      WRITE(message_text,'(3f8.4,f8.3,3f8.4)') turbdiff_config(1)%tkhmin, turbdiff_config(1)%tkmmin, &
        turbdiff_config(1)%rlam_heat , turbdiff_config(1)%rat_sea, turbdiff_config(1)%alpha0,        &
        turbdiff_config(1)%alpha0_max, turbdiff_config(1)%alpha0_pert
      CALL message('Perturbed values, tkhmin, tkmmin, rlam_heat, rat_sea, alpha0_min/max/pert', TRIM(message_text))

      WRITE(message_text,'(f8.2,3f8.4)') turbdiff_config(1)%tur_len, turbdiff_config(1)%a_hshr, &
        tune_box_liq, tune_box_liq_asy
      CALL message('Perturbed values, tur_len, a_hshr, box_liq, box_liq_asy', TRIM(message_text))

      WRITE(message_text,'(2f8.4,e11.4)') tune_minsnowfrac, c_soil, cwimax_ml
      CALL message('Perturbed values, minsnowfrac, c_soil, cwimax_ml', TRIM(message_text))


      ! Reinitialization of randum number generator in order to make external parameter perturbations
      ! independent of the number of RANDOM_NUMBER calls so far
      DO i = 1, rnd_size
        rnd_seed(i) = (139+i)*ipn - (23+i**2)*(5+MOD(ipn,12))**2 + 4*i**3
      ENDDO
      CALL RANDOM_SEED(PUT=rnd_seed)
      DO i = 1, 10+ipn
        CALL RANDOM_NUMBER(rnd_num)
      ENDDO

      ALLOCATE(rnd_tkred_sfc(nclass_lu(1)))

      CALL message('','')
      CALL message('','Perturbed external parameters: roughness length, root depth, min. stomata resistance,&
                   & max. leaf area index; defaults in brackets')
      ! Perturbations for external parameters specified depending on the land cover class
      DO i = 1, nclass_lu(1) ! we assume here that the same land cover dataset is used for all model domains

        ! roughness length
        CALL random_gen(rnd_num)
        z0_lcc  = ext_data(1)%atm%z0_lcc(i) * (1._wp + 2._wp*(rnd_num-0.5_wp)*range_z0_lcc)

        ! root depth
        CALL random_gen(rnd_num)
        rootdp  = ext_data(1)%atm%rootdmax_lcc(i) * (1._wp + 2._wp*(rnd_num-0.5_wp)*range_rootdp)

        ! minimum stomata resistance
        CALL random_gen(rnd_num)
        rsmin   = ext_data(1)%atm%stomresmin_lcc(i) * (1._wp + 2._wp*(rnd_num-0.5_wp)*range_rsmin)

        ! leaf area index
        CALL random_gen(rnd_num)
        laimax  = ext_data(1)%atm%laimax_lcc(i) * (1._wp + 2._wp*(rnd_num-0.5_wp)*range_laimax)

        WRITE(message_text,'(a,i3,2f8.4,f8.1,f8.3,2x,a,2f8.4,f8.1,f8.3,a)') 'Land-cover class:', i,       &
          z0_lcc, rootdp, rsmin, laimax, '(', ext_data(1)%atm%z0_lcc(i), ext_data(1)%atm%rootdmax_lcc(i), &
          ext_data(1)%atm%stomresmin_lcc(i), ext_data(1)%atm%laimax_lcc(i), ' )'
        CALL message('', TRIM(message_text))

        DO jg = 1, n_dom
          ext_data(jg)%atm%z0_lcc(i)         = z0_lcc
          ext_data(jg)%atm%rootdmax_lcc(i)   = rootdp
          ext_data(jg)%atm%stomresmin_lcc(i) = rsmin
          ext_data(jg)%atm%laimax_lcc(i)     = laimax
        ENDDO

        ! store random number for subsequent (in SR compute_ensemble_pert) computation of perturbation of 
        ! minimum diffusion coefficients near the surface
        CALL RANDOM_NUMBER(rnd_num)
        rnd_tkred_sfc(i) = rnd_num

      ENDDO

      DEALLOCATE(rnd_seed)

      ! Calculate correction term for systematically increased evaporation related to SST perturbations
      ! To simplify the calculation, we make the approximation that a 10 K temperature change corresponds
      ! to a factor 2 in the saturation vapor pressure, which is well justified between 0 and 25 deg C
      !
      svp_pert = 0.1_wp*stdev_sst_pert*LOG(2._wp)
      sst_pert_corrfac = 2._wp / (EXP(svp_pert) + EXP(-svp_pert))

      WRITE(message_text,'(f8.5)') sst_pert_corrfac
      CALL message('e_sat correction factor for SST perturbations', TRIM(message_text))

    ELSE

      sst_pert_corrfac = 1._wp

    ENDIF

  END SUBROUTINE configure_ensemble_pert


  !>
  !! Computation of array-based ensemble perturbation fields
  !!
  !! @par Revision History
  !! Initial revision by Guenther Zaengl, DWD (2016-04-08)
  !!
  SUBROUTINE compute_ensemble_pert(p_patch, ext_data, prm_diag, mtime_date)

    TYPE(t_patch),         INTENT(IN)    :: p_patch(:)
    TYPE(t_external_data), INTENT(IN)    :: ext_data(:)
    TYPE(t_nwp_phy_diag),  INTENT(INOUT) :: prm_diag(:)
    TYPE(datetime),        POINTER       :: mtime_date

    INTEGER  :: jg, jb, jc, jt, ilu, iyr
    INTEGER  :: rl_start, rl_end, i_startblk, i_endblk, i_startidx, i_endidx
    REAL(wp) :: wrnd_num(nproma), log_range_tkred, ssny, secyr, phaseshift

    rl_start = grf_bdywidth_c+1
    rl_end   = min_rlcell_int

    log_range_tkred = LOG(range_tkred_sfc)

    ! ssny = seconds since New Year
    ssny = (getDayOfYearFromDateTime(mtime_date)-1)*86400._wp + mtime_date%time%hour*3600._wp + &
           mtime_date%time%minute*60._wp + mtime_date%time%second

    iyr = mtime_date%date%year
    IF (MOD(iyr,4) == 0 .AND. .NOT. (MOD(iyr,100) == 0 .AND. MOD(iyr,400) /= 0)) THEN
      secyr = 366._wp*86400._wp
    ELSE
      secyr = 365._wp*86400._wp
    ENDIF
    phaseshift = 20._wp*ssny/secyr
    phaseshift = phaseshift - INT(phaseshift)

!$OMP PARALLEL PRIVATE(jg,i_startblk,i_endblk)
    DO jg = 1, n_dom

      i_startblk = p_patch(jg)%cells%start_block(rl_start)
      i_endblk   = p_patch(jg)%cells%end_block(rl_end)

      IF (use_ensemble_pert .AND. gribout_config(1)%perturbationNumber >= 1) THEN

!$OMP DO PRIVATE(jb,jc,jt,i_startidx,i_endidx,wrnd_num,ilu)
        DO jb = i_startblk, i_endblk

          CALL get_indices_c(p_patch(jg), jb, i_startblk, i_endblk, &
            i_startidx, i_endidx, rl_start, rl_end)

          wrnd_num(:) = 0._wp
          DO jt = 1, ntiles_total+ntiles_water
            IF (jt <= ntiles_lnd .OR. jt >= ntiles_total+1) THEN
              DO jc = i_startidx, i_endidx
                ilu = MAX(1,ext_data(jg)%atm%lc_class_t(jc,jb,jt))
                wrnd_num(jc) = wrnd_num(jc) + rnd_tkred_sfc(ilu)*ext_data(jg)%atm%lc_frac_t(jc,jb,jt)
              ENDDO
            ENDIF
          ENDDO
          DO jc = i_startidx, i_endidx
            prm_diag(jg)%tkred_sfc(jc,jb) = EXP(log_range_tkred*SIN(pi2*(wrnd_num(jc)+phaseshift)))
          ENDDO

        ENDDO
!$OMP END DO
      ELSE
!$OMP DO
        DO jb = i_startblk, i_endblk
          prm_diag(jg)%tkred_sfc(:,jb) = 1._wp
        ENDDO
!$OMP END DO
      ENDIF
    ENDDO
!$OMP END PARALLEL

  END SUBROUTINE compute_ensemble_pert

  ! Auxiliary routine to switch between equally distributed and discrete random numbers
  SUBROUTINE random_gen(rnd_val)

    REAL(wp), INTENT(OUT) :: rnd_val

    REAL(wp) :: rnd_aux
    CALL RANDOM_NUMBER(rnd_aux)

    IF (itype_pert_gen == 1) THEN
      rnd_val = rnd_aux
    ELSE IF (itype_pert_gen == 2) THEN
      IF (rnd_aux < 0.25_wp) THEN
        rnd_val = 0._wp
      ELSE IF (rnd_aux <= 0.75_wp) THEN
        rnd_val = 0.5_wp
      ELSE
        rnd_val = 1._wp
      ENDIF
    ENDIF

  END SUBROUTINE random_gen

END MODULE mo_ensemble_pert_config
