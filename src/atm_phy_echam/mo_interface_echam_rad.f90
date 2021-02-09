!>
!! @brief Subroutine interface_echam_rad calls the radiative transfer scheme.
!!
!! @author Hui Wan, MPI-M
!! @author Marco Giorgetta, MPI-M
!!
!! @par Revision History
!!  Original version from ECHAM6 (revision 2028)
!!  Modified for ICOHAM by Hui Wan and Marco Giorgetta (2010)
!!  Modified for ICONAM by Marco Giorgetta (2014)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_interface_echam_rad

  USE mo_kind,                ONLY: wp
  USE mtime,                  ONLY: datetime, OPERATOR(<=), OPERATOR(>)

  USE mo_model_domain        ,ONLY: t_patch
!   USE mo_loopindices         ,ONLY: get_indices_c
!   USE mo_impl_constants      ,ONLY: min_rlcell_int, grf_bdywidth_c

  USE mo_parallel_config     ,ONLY: nproma
  USE mo_run_config,          ONLY: nlev, nlevp1
  USE mo_echam_phy_memory,    ONLY: t_echam_phy_field, prm_field

  USE mo_psrad_radiation,     ONLY: psrad_radiation

  USE mo_timer,               ONLY: ltimer, timer_start, timer_stop, timer_radiation

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: interface_echam_rad

CONTAINS

  !-------------------------------------------------------------------
  SUBROUTINE interface_echam_rad(is_in_sd_ed_interval, &
       &                         is_active,            &
       &                         patch,                &
       &                         datetime_old          )

    LOGICAL                 ,INTENT(in)    :: is_in_sd_ed_interval
    LOGICAL                 ,INTENT(in)    :: is_active
    TYPE(t_patch)   ,TARGET ,INTENT(in)    :: patch
    TYPE(datetime)          ,POINTER       :: datetime_old

    ! Local variables
    !
    TYPE(t_echam_phy_field) ,POINTER    :: field
    INTEGER :: jg
    INTEGER :: itype(nproma,patch%nblks_c)              !< type of convection
    LOGICAL :: loland(nproma,patch%nblks_c)
    LOGICAL :: loglac(nproma,patch%nblks_c)

    IF (ltimer) CALL timer_start(timer_radiation)

    ! grid index
    jg = patch%id

    ! associate pointers
    field => prm_field(jg)

    IF ( is_in_sd_ed_interval ) THEN
        !
        IF ( is_active ) THEN
          !
          ! store ts_rad of this radiatiative transfer timestep in ts_rad_rt,
          ! so that it can be reused in radheat in the other timesteps
          field%ts_rad_rt(:,:) = field%ts_rad(:,:)
          !
          itype (:,:) = NINT(field%rtype(:,:))
          loland(:,:) = field% sftlf (:,:) > 0._wp
          loglac(:,:) = field% sftgif(:,:) > 0._wp
          !
          CALL psrad_radiation(                            &
              & patch                                     ,&
              & klev           = nlev                     ,&!< in  number of full levels = number of layers
              & klevp1         = nlevp1                   ,&!< in  number of half levels = number of layer interfaces
              & ktype          = itype(:,:)               ,&!< in  type of convection
              & loland         = loland(:,:)              ,&!< in  land-sea mask. (logical)
              & loglac         = loglac(:,:)              ,&!< in  glacier mask (logical)
              & this_datetime  = datetime_old             ,&!< in  actual time step
              & pcos_mu0       = field%cosmu0_rt(:,:)     ,&!< in  solar zenith angle
              & daylght_frc    = field%daylght_frc_rt(:,:),&!in daylight fraction
              & alb_vis_dir    = field%albvisdir(:,:)     ,&!< in  surface albedo for visible range, direct
              & alb_nir_dir    = field%albnirdir(:,:)     ,&!< in  surface albedo for near IR range, direct
              & alb_vis_dif    = field%albvisdif(:,:)     ,&!< in  surface albedo for visible range, diffuse
              & alb_nir_dif    = field%albnirdif(:,:)     ,&!< in  surface albedo for near IR range, diffuse
              & emissivity     = field%emissivity(:,:)    ,&!< in  surface longwave emissivity
              & tk_sfc         = field%ts_rad_rt(:,:)     ,&!< in  grid box mean surface temperature
              & zf             = field%zf(:,:,:)          ,&!< in  geometric height at full level      [m]
              & zh             = field%zh(:,:,:)          ,&!< in  geometric height at half level      [m]
              & dz             = field%dz(:,:,:)          ,&!< in  geometric height thickness of layer [m]
              & pp_hl          = field%phalf(:,:,:)       ,&!< in  pressure at half levels [Pa]
              & pp_fl          = field%pfull(:,:,:)       ,&!< in  pressure at full levels [Pa]
              & tk_fl          = field%ta(:,:,:)          ,&!< in  tk_fl  = temperature at full level at t-dt
              & xm_dry         = field%mdry(:,:,:)        ,&!< in  dry air mass in layer [kg/m2]
              & xm_trc         = field%mtrc(:,:,:,:)      ,&!< in  tracer  mass in layer [kg/m2]
              & xm_ozn         = field%o3(:,:,:)          ,&!< inout  ozone  mass mixing ratio [kg/kg]
              !
              & cdnc           = field% acdnc(:,:,:)      ,&!< in   cloud droplet number conc
              & cld_frc        = field% aclc(:,:,:)       ,&!< in   cloud fraction [m2/m2]
              & cld_cvr        = field%aclcov(:,:)        ,&!< out  total cloud cover
              !
              & lw_dnw_clr     = field%rldcs_rt(:,:,:)    ,&!< out  Clear-sky net longwave  at all levels
              & lw_upw_clr     = field%rlucs_rt(:,:,:)    ,&!< out  Clear-sky net longwave  at all levels
              & sw_dnw_clr     = field%rsdcs_rt(:,:,:)    ,&!< out  Clear-sky net shortwave at all levels
              & sw_upw_clr     = field%rsucs_rt(:,:,:)    ,&!< out  Clear-sky net shortwave at all levels
              & lw_dnw         = field%rld_rt  (:,:,:)    ,&!< out  All-sky net longwave  at all levels
              & lw_upw         = field%rlu_rt  (:,:,:)    ,&!< out  All-sky net longwave  at all levels
              & sw_dnw         = field%rsd_rt  (:,:,:)    ,&!< out  All-sky net longwave  at all levels
              & sw_upw         = field%rsu_rt  (:,:,:)    ,&!< out  All-sky net longwave  at all levels
              !
              & vis_dn_dir_sfc = field%rvds_dir_rt(:,:)   ,&!< out  all-sky downward direct visible radiation at surface
              & par_dn_dir_sfc = field%rpds_dir_rt(:,:)   ,&!< out  all-sky downward direct PAR     radiation at surface
              & nir_dn_dir_sfc = field%rnds_dir_rt(:,:)   ,&!< out  all-sky downward direct near-IR radiation at surface
              & vis_dn_dff_sfc = field%rvds_dif_rt(:,:)   ,&!< out  all-sky downward diffuse visible radiation at surface
              & par_dn_dff_sfc = field%rpds_dif_rt(:,:)   ,&!< out  all-sky downward diffuse PAR     radiation at surface
              & nir_dn_dff_sfc = field%rnds_dif_rt(:,:)   ,&!< out  all-sky downward diffuse near-IR radiation at surface
              & vis_up_sfc     = field%rvus_rt    (:,:)   ,&!< out  all-sky upward visible radiation at surface
              & par_up_sfc     = field%rpus_rt    (:,:)   ,&!< out  all-sky upward PAR     radiation at surfac
              & nir_up_sfc     = field%rnus_rt    (:,:)   ,&!< out  all-sky upward near-IR radiation at surface
              & aer_aod_533    = field%aer_aod_533 (:,:,:),&!< out  aerosol optical density at 533 nm
              & aer_ssa_533    = field%aer_ssa_533 (:,:,:),&!< out  single scattering albedo at 533 nm
              & aer_asy_533    = field%aer_asy_533 (:,:,:),&!< out  asymmetrie factor at 533 nm
              & aer_aod_2325   = field%aer_aod_2325(:,:,:),&!< out  aerosol optical density at 2325 nm
              & aer_ssa_2325   = field%aer_ssa_2325(:,:,:),&!< out  single scattering albedo at 2325 nm
              & aer_asy_2325   = field%aer_asy_2325(:,:,:),&!< out  asymmetrie factor at 2325 nm
              & aer_aod_9731   = field%aer_aod_9731(:,:,:) &!< out  aerosol optical density at 9731
              )
          !
          END IF
          !
       ELSE
          !
          ! LW
          field%rldcs_rt(:,:,:)  = 0.0_wp !< out  Clear-sky net longwave  at all levels
          field%rlucs_rt(:,:,:)  = 0.0_wp !< out  Clear-sky net longwave  at all levels
          field%rld_rt  (:,:,:)  = 0.0_wp !< out  All-sky net longwave  at all levels
          field%rlu_rt  (:,:,:)  = 0.0_wp !< out  All-sky net longwave  at all levels
          !
          ! SW all
          field%rsdcs_rt(:,:,:)  = 0.0_wp !< out  Clear-sky net shortwave at all levels
          field%rsucs_rt(:,:,:)  = 0.0_wp !< out  Clear-sky net shortwave at all levels
          field%rsd_rt  (:,:,:)  = 0.0_wp !< out  All-sky net longwave  at all levels
          field%rsu_rt  (:,:,:)  = 0.0_wp !< out  All-sky net longwave  at all levels
          !
          ! SW vis, par and nir
          field%rvds_dir_rt(:,:) = 0.0_wp !< out  all-sky downward direct visible radiation at surface
          field%rpds_dir_rt(:,:) = 0.0_wp !< all-sky downward direct PAR     radiation at surface
          field%rnds_dir_rt(:,:) = 0.0_wp !< all-sky downward direct near-IR radiation at surface
          field%rvds_dif_rt(:,:) = 0.0_wp !< all-sky downward diffuse visible radiation at surface
          field%rpds_dif_rt(:,:) = 0.0_wp !< all-sky downward diffuse PAR     radiation at surface
          field%rnds_dif_rt(:,:) = 0.0_wp !< all-sky downward diffuse near-IR radiation at surface
          field%rvus_rt    (:,:) = 0.0_wp !< all-sky upward visible radiation at surface
          field%rpus_rt    (:,:) = 0.0_wp !< all-sky upward PAR     radiation at surfac
          field%rnus_rt    (:,:) = 0.0_wp !< all-sky upward near-IR radiation at surface
          !
          ! total cloud cover diagnostics
          field%aclcov(:,:)      = 0.0_wp !< out  total cloud cover
          !
          ! aerosol optical properties diagnostics
          field%aer_aod_533 (:,:,:) = 0.0_wp
          field%aer_ssa_533 (:,:,:) = 0.0_wp
          field%aer_asy_533 (:,:,:) = 0.0_wp
          field%aer_aod_2325(:,:,:) = 0.0_wp
          field%aer_ssa_2325(:,:,:) = 0.0_wp
          field%aer_asy_2325(:,:,:) = 0.0_wp
          field%aer_aod_9731(:,:,:) = 0.0_wp
       !
       END IF

     IF (ltimer) CALL timer_stop(timer_radiation)

  END SUBROUTINE interface_echam_rad
  !-------------------------------------------------------------------

END MODULE mo_interface_echam_rad
