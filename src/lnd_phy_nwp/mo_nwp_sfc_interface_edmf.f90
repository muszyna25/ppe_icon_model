!>
!! This module is the interface between the EDMF turbulence to the
!! surface parameterisations:
!! inwp_sfc  == 1 == surface scheme TERRA
!! inwp_turb == 3 == turbulence EDMF
!!
!! @author Martin Koehler, DWD, Offenbach (2012-05-04)
!!
!! @par Revision History
!! Initial Martin Koehler, DWD, Offenbach (2012-05-04)
!!
!! Modifications by Dmitrii Mironov, DWD (2016-08-08)
!! - Changes related to the use of prognostic the sea-ice albedo.
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
#if defined __xlC__
@PROCESS SPILL(2160)
#endif
MODULE mo_nwp_sfc_interface_edmf

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: message !, message_text
  USE mo_model_domain,        ONLY: t_patch
  USE mo_impl_constants,      ONLY: zml_soil
  USE mo_ext_data_types,      ONLY: t_external_data
  USE mo_parallel_config,     ONLY: nproma
  USE mo_run_config,          ONLY: msg_level
  USE mo_atm_phy_nwp_config,  ONLY: atm_phy_nwp_config
  USE mo_nwp_phy_state,       ONLY: phy_params, prm_diag
  USE mo_lnd_nwp_config,      ONLY: nlev_soil, nlev_snow, ibot_w_so, ntiles_total,  &
    &                               ntiles_water, lseaice, llake, lmulti_snow,      &
    &                               ntiles_lnd, lsnowtile, isub_water, isub_seaice, &
    &                               isub_lake, itype_interception
  USE mo_satad,               ONLY: sat_pres_water, sat_pres_ice, spec_humi
  USE mo_soil_ml,             ONLY: terra_multlay
  USE mo_nwp_sfc_utils,       ONLY: diag_snowfrac_tg, update_idx_lists_lnd, update_idx_lists_sea
  USE mo_flake,               ONLY: flake_interface
  USE mo_seaice_nwp,          ONLY: seaice_timestep_nwp
  USE mo_phyparam_soil              ! soil and vegetation parameters for TILES
  USE mo_physical_constants,  ONLY: tmelt


  IMPLICIT NONE

  PRIVATE



  PUBLIC  ::  nwp_surface_edmf


#ifdef __SX__
! parameters for loop unrolling
INTEGER, PARAMETER :: nlsoil= 8
INTEGER, PARAMETER :: nlsnow= 2
#endif


CONTAINS
  !!
  !!-------------------------------------------------------------------------
  !!

  SUBROUTINE nwp_surface_edmf (&
                  ext_data         , & !>in
                  jb               , & ! block
                  jg               , & ! patch
                  i_startidx       , & ! start index for computations in the parallel program
                  i_endidx         , & ! end index for computations in the parallel program
                  tcall_sfc_jg     , & ! time step                                     ( s  )

                  u_ex             , & ! zonal wind speed                              ( m/s )
                  v_ex             , & ! meridional wind speed                         ( m/s )
                  t_ex             , & ! temperature                                   (  k  )
                  qv_ex            , & ! specific water vapor content                  (kg/kg)
                  p0_ex            , & !!!! base state pressure                           (Pa)
                  ps_ex            , & ! surface pressure                              ( pa  )

                  t_snow_ex        , & ! temperature of the snow-surface               (  K  )
                  t_snow_mult_ex   , & ! temperature of the snow-surface               (  K  )
                  t_s_ex           , & ! temperature of the ground surface             (  K  )
                  t_g_ex           , & ! weighted surface temperature                  (  K  )
                  qv_s_ex          , & ! specific humidity at the surface              (kg/kg)
                  w_snow_ex        , & ! water content of snow                         (m H2O)
                  rho_snow_ex      , & ! snow density                                  (kg/m**3)
                  rho_snow_mult_ex , & ! snow density                                  (kg/m**3)
                  h_snow_ex        , & ! snow height                                   (  m  )
                  w_i_ex           , & ! water content of interception water           (m H2O)
                  w_p_ex           , & ! water content of pond interception water      (m H2O)
                  w_s_ex           , & ! water content of interception snow            (m H2O)
                  t_so_ex          , & ! soil temperature (main level)                 (  K  )
                  w_so_ex          , & ! total water conent (ice + liquid water)       (m H20)
                  w_so_ice_ex      , & ! ice content                                   (m H20)
!                 t_2m_ex          , & ! temperature in 2m                             (  K  )
                  u_10m_ex         , & ! zonal wind in 10m                             ( m/s )
                  v_10m_ex         , & ! meridional wind in 10m                        ( m/s )

                  freshsnow_ex     , & ! indicator for age of snow in top of snow layer(  -  )
                  snowfrac_lc_ex   , & ! snow-cover fraction (for each land cover)     (  -  )
                  snowfrac_ex      , & ! snow-cover fraction (including snow tiles)    (  -  )
                  wliq_snow_ex     , & ! liquid water content in the snow              (m H2O)
                  wtot_snow_ex     , & ! total (liquid + solid) water content of snow  (m H2O)
                  dzh_snow_ex      , & ! layer thickness between half levels in snow   (  m  )

                  prr_con_ex       , & ! precipitation rate of rain, convective        (kg/m2*s)
                  prs_con_ex       , & ! precipitation rate of snow, convective        (kg/m2*s)
                  prr_gsp_ex       , & ! precipitation rate of rain, grid-scale        (kg/m2*s)
                  prs_gsp_ex       , & ! precipitation rate of snow, grid-scale        (kg/m2*s)

                  tch_ex           , & ! turbulent transfer coefficient for heat       ( -- )
                  tcm_ex           , & ! turbulent transfer coefficient for momentum   ( -- )
                  tfv_ex           , & ! laminar reduction factor for evaporation      ( -- )

                  sobs_ex          , & ! solar radiation at the ground                 ( W/m2)
                  thbs_ex          , & ! thermal radiation at the ground               ( W/m2)
                  pabs_ex          , & !!!! photosynthetic active radiation            ( W/m2)

                  runoff_s_ex      , & ! surface water runoff; sum over forecast       (kg/m2)
                  runoff_g_ex      , & ! soil water runoff; sum over forecast          (kg/m2)

                  t_g              , & ! surface temperature (grid mean)               ( K )
                  qv_s             , & ! surface specific humidity (grid mean)         (kg/kg)

                  t_ice            , & ! sea ice temperature                           (  K  )
                  h_ice            , & ! sea ice height                                (  m  )
                  t_snow_si        , & ! sea ice snow temperature                      (  K  )
                  h_snow_si        , & ! sea ice snow height                           (  m  )
                  alb_si           , & ! sea-ice albedo                                (  -  )
                  fr_seaice        , & ! sea ice fraction                              (  1  )

                  shfl_soil_ex     , & ! sensible heat flux soil/air interface         (W/m2)
                  lhfl_soil_ex     , & ! latent   heat flux soil/air interface         (W/m2)
                  shfl_snow_ex     , & ! sensible heat flux snow/air interface         (W/m2)
                  lhfl_snow_ex     , & ! latent   heat flux snow/air interface         (W/m2)
                  shfl_s_ex        , & ! sensible heat flux                            (W/m2)
                  lhfl_s_ex        , & ! latent   heat flux                            (W/m2)
                  qhfl_s_ex        , & ! moisture      flux                            (W/m2)
                  lhfl_bs_ex       , &
                  lhfl_pl_ex       , &
                  rstom_ex           )


  IMPLICIT NONE

  INTEGER,  INTENT(IN)  ::  &
                  i_startidx       , & ! start index for computations in the parallel program
                  i_endidx         , & ! end index for computations in the parallel program
                  jb               , & ! block
                  jg                   ! patch
  REAL(wp), INTENT(IN)  ::  &
                  tcall_sfc_jg         ! time interval for surface                     (  s  )
  REAL(wp), DIMENSION(nproma), INTENT(IN) :: &
                  u_ex             , & ! zonal wind speed                              ( m/s )
                  v_ex             , & ! meridional wind speed                         ( m/s )
                  t_ex             , & ! temperature                                   (  k  )
                  qv_ex            , & ! specific water vapor content                  (kg/kg)
                  p0_ex            , & !!!! base state pressure                        ( Pa )
                  ps_ex                ! surface pressure                              ( pa  )
  REAL(wp), DIMENSION(nproma,nlev_snow+1,ntiles_total), INTENT(INOUT) :: &
                  t_snow_mult_ex       ! temperature of the snow-surface               (  K  )
  REAL(wp), DIMENSION(nproma,nlev_snow,ntiles_total), INTENT(INOUT) :: &
                  rho_snow_mult_ex     ! snow density                                  (kg/m**3)
  REAL(wp), DIMENSION(nproma,ntiles_total+ntiles_water), INTENT(INOUT) :: &
                  t_g_ex           , & ! weighted surface temperature                  (  K  )
                  qv_s_ex              ! specific humidity at the surface              (kg/kg)
  REAL(wp), DIMENSION(nproma,ntiles_total), INTENT(INOUT) :: &
                  t_snow_ex        , & ! temperature of the snow-surface (K)
                  t_s_ex           , & ! temperature of the ground surface             (  K  )
                  w_snow_ex        , & ! water content of snow                         (m H2O)
                  rho_snow_ex      , & ! snow density                                  (kg/m**3)
                  h_snow_ex        , & ! snow height
                  w_i_ex           , & ! water content of interception water           (m H2O)
                  w_p_ex           , & ! water content of interception water           (m H2O)
                  w_s_ex               ! water content of interception water           (m H2O)

  REAL(wp), DIMENSION(nproma,nlev_soil+1,ntiles_total), INTENT(INOUT) :: &
                  t_so_ex              ! soil temperature (main level)                 (  K  )
  REAL(wp), DIMENSION(nproma,nlev_soil,ntiles_total), INTENT(INOUT) :: &
                  w_so_ex          , & ! total water conent (ice + liquid water)       (m H20)
                  w_so_ice_ex          ! ice content                                   (m H20)
  REAL(wp), DIMENSION(nproma), INTENT(INOUT) :: &
!                 t_2m_ex          , & ! temperature in 2m                             (  K  )
                  u_10m_ex         , & ! zonal wind in 10m                             ( m/s )
                  v_10m_ex             ! meridional wind in 10m                        ( m/s )
  REAL(wp), DIMENSION(nproma,ntiles_total), INTENT(INOUT) :: &
                  freshsnow_ex     , & ! indicator for age of snow in top of snow layer(  -  )
                  snowfrac_lc_ex   , & ! snow-cover fraction                           (  -  )
                  snowfrac_ex          ! snow-cover fraction                           (  -  )
  REAL(wp), DIMENSION(nproma,nlev_snow,ntiles_total), INTENT(INOUT) :: &
                  wliq_snow_ex     , & ! liquid water content in the snow              (m H2O)
                  wtot_snow_ex     , & ! total (liquid + solid) water content of snow  (m H2O)
                  dzh_snow_ex          ! layer thickness between half levels in snow   (  m  )
  REAL(wp), DIMENSION(nproma), INTENT(IN) ::    &
                  prr_con_ex       , & ! precipitation rate of rain, convective        (kg/m2*s)
                  prs_con_ex       , & ! precipitation rate of snow, convective        (kg/m2*s)
                  prr_gsp_ex       , & ! precipitation rate of rain, grid-scale        (kg/m2*s)
                  prs_gsp_ex           ! precipitation rate of snow, grid-scale        (kg/m2*s)
  REAL(wp), DIMENSION(nproma,ntiles_total+ntiles_water), INTENT(INOUT) :: &
                  tch_ex           , & ! turbulent transfer coefficient for heat       ( -- )
                  tcm_ex           , & ! turbulent transfer coefficient for momentum   ( -- )
                  tfv_ex               ! laminar reduction factor for evaporation      ( -- )
  REAL(wp), DIMENSION(nproma,ntiles_total+ntiles_water), INTENT(IN) :: &
                  sobs_ex          , & ! solar radiation at the ground                 ( W/m2)
                  thbs_ex          , & ! thermal radiation at the ground               ( W/m2)
                  pabs_ex              !!!! photosynthetic active radiation            ( W/m2)
  REAL(wp), DIMENSION(nproma,ntiles_total), INTENT(INOUT) :: &
                  runoff_s_ex      , & ! surface water runoff; sum over forecast       (kg/m2)
                  runoff_g_ex          ! soil water runoff; sum over forecast          (kg/m2)
  REAL(wp), DIMENSION(nproma), INTENT(INOUT) :: &
                  t_g              , &
                  qv_s
  REAL(wp), DIMENSION(nproma), INTENT(INOUT) :: &
                  t_ice            , & ! sea ice temperature                           (  K  )
                  h_ice            , & ! sea ice height                                (  m  )
                  t_snow_si        , & ! sea ice snow temperature                      (  K  )
                  h_snow_si        , & ! sea ice snow height                           (  m  )
                  alb_si               ! sea-ice albedo                                (  -  )
  REAL(wp), DIMENSION(nproma), INTENT(INOUT) :: &
                  fr_seaice            ! sea ice fraction                              (  1  )
  REAL(wp), DIMENSION(nproma,ntiles_total+ntiles_water), INTENT(INOUT) :: &
                  shfl_soil_ex     , & ! sensible heat flux soil/air interface         ( W/m2)
                  lhfl_soil_ex     , & ! latent   heat flux soil/air interface         ( W/m2)
                  shfl_snow_ex     , & ! sensible heat flux snow/air interface         ( W/m2)
                  lhfl_snow_ex     , & ! latent   heat flux snow/air interface         ( W/m2)
                  shfl_s_ex        , & ! sensible heat flux                            ( W/m2)
                  lhfl_s_ex        , & ! latent   heat flux                            ( W/m2)
                  qhfl_s_ex        , & ! moisture      flux                            ( W/m2)
                  lhfl_bs_ex       , &
                  rstom_ex
  REAL(wp), DIMENSION(nproma,nlev_soil,ntiles_total+ntiles_water), INTENT(INOUT) :: &
                  lhfl_pl_ex

  TYPE(t_external_data), INTENT(inout) :: ext_data        !< external data


    ! Local array bounds:
    !
    INTEGER :: isubs, isubs_snow

    ! Local scalars:
    !
    INTEGER :: jc,jk      !loop indices


    REAL(wp) :: ps_t       (nproma)
    REAL(wp) :: prr_con_t  (nproma)
    REAL(wp) :: prs_con_t  (nproma)
    REAL(wp) :: conv_frac  (nproma)
    REAL(wp) :: prr_gsp_t  (nproma)
    REAL(wp) :: prs_gsp_t  (nproma)

    REAL(wp) :: u_t (nproma)
    REAL(wp) :: v_t (nproma)
    REAL(wp) :: t_t (nproma)
    REAL(wp) :: qv_t(nproma)
    REAL(wp) :: p0_t(nproma)

    REAL(wp) :: sso_sigma_t(nproma)
    INTEGER  :: lc_class_t (nproma, ntiles_total)

    REAL(wp) :: t_snow_now_t (nproma, ntiles_total)
    REAL(wp) :: t_snow_new_t (nproma, ntiles_total)

    REAL(wp) :: t_s_now_t  (nproma, ntiles_total)
    REAL(wp) :: t_s_new_t  (nproma, ntiles_total)

    REAL(wp) :: t_g_t      (nproma, ntiles_total)
    REAL(wp) :: qv_s_t     (nproma, ntiles_total)

    REAL(wp) :: w_snow_now_t (nproma, ntiles_total)
    REAL(wp) :: w_snow_new_t (nproma, ntiles_total)

    REAL(wp) :: rho_snow_now_t (nproma, ntiles_total)
    REAL(wp) :: rho_snow_new_t (nproma, ntiles_total)

    REAL(wp) :: h_snow_t   (nproma, ntiles_total)
    REAL(wp) :: meltrate   (nproma)

    REAL(wp) :: w_i_now_t  (nproma, ntiles_total)
    REAL(wp) :: w_i_new_t  (nproma, ntiles_total)

    REAL(wp) :: w_p_now_t  (nproma, ntiles_total)
    REAL(wp) :: w_p_new_t  (nproma, ntiles_total)

    REAL(wp) :: w_s_now_t  (nproma, ntiles_total)
    REAL(wp) :: w_s_new_t  (nproma, ntiles_total)

!   REAL(wp) :: t_2m_t     (nproma, ntiles_total)
    REAL(wp) :: u_10m_t    (nproma, ntiles_total)
    REAL(wp) :: v_10m_t    (nproma, ntiles_total)
    REAL(wp) :: freshsnow_t(nproma, ntiles_total)
    REAL(wp) :: snowfrac_lc_t(nproma, ntiles_total)
    REAL(wp) :: snowfrac_t (nproma, ntiles_total)

    REAL(wp) :: tch_t      (nproma, ntiles_total)
    REAL(wp) :: tcm_t      (nproma, ntiles_total)
    REAL(wp) :: tfv_t      (nproma, ntiles_total)

    REAL(wp) :: sobs_t     (nproma, ntiles_total)
    REAL(wp) :: thbs_t     (nproma, ntiles_total)
    REAL(wp) :: pabs_t     (nproma, ntiles_total)

    REAL(wp) :: runoff_s_t (nproma, ntiles_total)
    REAL(wp) :: runoff_g_t (nproma, ntiles_total)

    INTEGER  :: soiltyp_t (nproma, ntiles_total)
    REAL(wp) :: plcov_t   (nproma, ntiles_total)
    REAL(wp) :: rootdp_t  (nproma, ntiles_total)
    REAL(wp) :: sai_t     (nproma, ntiles_total)
    REAL(wp) :: tai_t     (nproma, ntiles_total)
    REAL(wp) :: eai_t     (nproma, ntiles_total)
    REAL(wp) :: rsmin2d_t (nproma, ntiles_total)

    ! local dummy variable for precipitation rate of graupel, grid-scale
    REAL(wp) :: dummy_prg_gsp(nproma)

    REAL(wp) :: t_snow_mult_now_t(nproma, nlev_snow+1, ntiles_total)
    REAL(wp) :: t_snow_mult_new_t(nproma, nlev_snow+1, ntiles_total)

    REAL(wp) :: rho_snow_mult_now_t(nproma, nlev_snow, ntiles_total)
    REAL(wp) :: rho_snow_mult_new_t(nproma, nlev_snow, ntiles_total)

    REAL(wp) :: wliq_snow_now_t(nproma, nlev_snow, ntiles_total)
    REAL(wp) :: wliq_snow_new_t(nproma, nlev_snow, ntiles_total)

    REAL(wp) :: wtot_snow_now_t(nproma, nlev_snow, ntiles_total)
    REAL(wp) :: wtot_snow_new_t(nproma, nlev_snow, ntiles_total)

    REAL(wp) :: dzh_snow_now_t(nproma, nlev_snow, ntiles_total)
    REAL(wp) :: dzh_snow_new_t(nproma, nlev_snow, ntiles_total)

    REAL(wp) :: t_so_now_t(nproma, nlev_soil+1, ntiles_total)
    REAL(wp) :: t_so_new_t(nproma, nlev_soil+1, ntiles_total)

    REAL(wp) :: w_so_now_t(nproma, nlev_soil, ntiles_total)
    REAL(wp) :: w_so_new_t(nproma, nlev_soil, ntiles_total)

    REAL(wp) :: w_so_ice_now_t(nproma, nlev_soil, ntiles_total)
    REAL(wp) :: w_so_ice_new_t(nproma, nlev_soil, ntiles_total)

    INTEGER  :: i_count, i_count_snow, ic, icount_init, is1, is2, init_list(2*nproma), it1(nproma), it2(nproma)
    REAL(wp) :: tmp1, tmp2, tmp3
    REAL(wp) :: frac_sv(nproma), frac_snow_sv(nproma), fact1(nproma), fact2(nproma)
    REAL(wp) :: rain_gsp_rate(nproma, ntiles_total)
    REAL(wp) :: snow_gsp_rate(nproma, ntiles_total)
    REAL(wp) :: rain_con_rate(nproma, ntiles_total)
    REAL(wp) :: snow_con_rate(nproma, ntiles_total)
    REAL(wp), PARAMETER :: small = 1.E-06_wp

    REAL(wp) :: t_g_s(nproma), qv_s_s(nproma)

    REAL(wp) :: shfl_s_t    (nproma, ntiles_total) ! sensible heat flux sfc
    REAL(wp) :: lhfl_s_t    (nproma, ntiles_total) ! latent heat flux sfc
    REAL(wp) :: qhfl_s_t    (nproma, ntiles_total) ! moisture flux sfc
    REAL(wp) :: shfl_soil_t (nproma, ntiles_total) ! sensible heat flux sfc (snow free)
    REAL(wp) :: lhfl_soil_t (nproma, ntiles_total) ! latent heat flux sfc   (snow free)
    REAL(wp) :: shfl_snow_t (nproma, ntiles_total) ! sensible heat flux sfc (snow covered)
    REAL(wp) :: lhfl_snow_t (nproma, ntiles_total) ! latent heat flux sfc   (snow covered)
    REAL(wp) :: lhfl_bs_t   (nproma, ntiles_total)
    REAL(wp) :: lhfl_pl_t   (nproma, nlev_soil, ntiles_total)
    REAL(wp) :: rstom_t     (nproma, ntiles_total)

!--------------------------------------------------------------


    ! initialize dummy variable (precipitation rate of graupel, grid-scale)
    dummy_prg_gsp(1:nproma) = 0._wp

    ! local variables related to the blocking

    IF (msg_level >= 15) THEN
      CALL message('mo_nwp_sfc_interface: ', 'call land-surface scheme')
    ENDIF

    IF( atm_phy_nwp_config(jg)%inwp_surface == 0) THEN
      ! check dry case
      IF( atm_phy_nwp_config(jg)%inwp_satad == 0) THEN
        DO jc = i_startidx, i_endidx
          qv_s (jc) = 0._wp
        ENDDO
      ELSE
        !
        !> adjust  humidity at water surface because of changed surface pressure
        !
        DO jc = i_startidx, i_endidx
          qv_s (jc) = spec_humi( sat_pres_water(t_g(jc)) , ps_ex(jc) )
        ENDDO
      ENDIF
    ENDIF

    DO isubs = 1,ntiles_total
      DO jc = i_startidx, i_endidx
        shfl_soil_ex(jc,isubs) = 0.0_wp    !initialized fluxes over land only
        lhfl_soil_ex(jc,isubs) = 0.0_wp    !not over water/ice/lake
        shfl_snow_ex(jc,isubs) = 0.0_wp    ! -
        lhfl_snow_ex(jc,isubs) = 0.0_wp    ! -
        shfl_s_ex   (jc,isubs) = 0.0_wp    ! -
        lhfl_s_ex   (jc,isubs) = 0.0_wp    ! -
        qhfl_s_ex   (jc,isubs) = 0.0_wp    ! -
        shfl_soil_t (jc,isubs) = 0.0_wp
        lhfl_soil_t (jc,isubs) = 0.0_wp
        shfl_snow_t (jc,isubs) = 0.0_wp
        lhfl_snow_t (jc,isubs) = 0.0_wp
        shfl_s_t    (jc,isubs) = 0.0_wp
        lhfl_s_t    (jc,isubs) = 0.0_wp
        qhfl_s_t    (jc,isubs) = 0.0_wp
        snowfrac_lc_t(jc,isubs)= 0.0_wp
        snowfrac_t  (jc,isubs) = 0.0_wp
      ENDDO
    ENDDO

    IF (  atm_phy_nwp_config(jg)%inwp_surface == 1 ) THEN

       ! Copy precipitation fields for subsequent downscaling
       DO isubs = 1,ntiles_total
         i_count = ext_data%atm%gp_count_t(jb,isubs)
         IF (i_count == 0) CYCLE ! skip loop if the index list for the given tile is empty
!CDIR NODEP,VOVERTAKE,VOB
         DO ic = 1, i_count
           jc = ext_data%atm%idx_lst_t(ic,jb,isubs)
           rain_gsp_rate(jc,isubs) = prr_gsp_ex(jc)
           snow_gsp_rate(jc,isubs) = prs_gsp_ex(jc)
           rain_con_rate(jc,isubs) = prr_con_ex(jc)
           snow_con_rate(jc,isubs) = prs_con_ex(jc)
         END DO
       END DO

!---------- Preparations for TERRA in the case if snow tiles are considered
       IF(lsnowtile) THEN      ! snow is considered as separate tiles

         DO isubs = 1, ntiles_lnd

           isubs_snow = isubs + ntiles_lnd
           i_count_snow = ext_data%atm%gp_count_t(jb,isubs_snow)

!CDIR NODEP,VOVERTAKE,VOB
           DO ic = 1, i_count_snow
             jc = ext_data%atm%idx_lst_t(ic,jb,isubs_snow)

             ! Snow and rain fall onto snow-covered tile surface only,
             ! if
             ! 1) the corresponding snow tile already exists and
             ! 2) the temperature of snow-free tile is below freezing point.
             ! If the temperature of snow-free tile is above freezing point,
             ! precipitation over it will be processed by this tile itself (no snow is created).
             ! If there is no snow tile so far at all, precipitation falls on the snow-free tile,
             ! and the snow tile will be created after TERRA.
             !
             ! DR after discussion with Ekaterina, changed from '>' to '<='
             IF(t_snow_ex(jc,isubs)  <= tmelt) THEN
               rain_gsp_rate(jc,isubs) = 0._wp
               snow_gsp_rate(jc,isubs) = 0._wp
               rain_con_rate(jc,isubs) = 0._wp
               snow_con_rate(jc,isubs) = 0._wp
             END IF
           END DO
         END DO
       END IF

!---------- Copy input fields for each tile

!----------------------------------
       DO isubs = 1,ntiles_total
!----------------------------------

        i_count = ext_data%atm%gp_count_t(jb,isubs)

        IF (i_count == 0) CYCLE ! skip loop if the index list for the given tile is empty

        DO ic = 1, i_count
          jc = ext_data%atm%idx_lst_t(ic,jb,isubs)

          ps_t(ic)      =  ps_ex        (jc)
          prr_con_t(ic) =  rain_con_rate(jc,isubs)
          prs_con_t(ic) =  snow_con_rate(jc,isubs)
          conv_frac(ic) =  phy_params(jg)%rcucov*     (1._wp - prm_diag(jg)%tropics_mask(jc,jb)) + &
                           phy_params(jg)%rcucov_trop*         prm_diag(jg)%tropics_mask(jc,jb)
          prr_gsp_t(ic) =  rain_gsp_rate(jc,isubs)
          prs_gsp_t(ic) =  snow_gsp_rate(jc,isubs)

          u_t(ic)       =  u_ex (jc)
          v_t(ic)       =  v_ex (jc)
          t_t(ic)       =  t_ex (jc)
          qv_t(ic)      =  qv_ex(jc)
          p0_t(ic)      =  p0_ex(jc)

          sso_sigma_t(ic)       = ext_data%atm%sso_stdh  (jc,jb)
          lc_class_t(ic,isubs)  = ext_data%atm%lc_class_t(jc,jb,isubs)

          t_snow_now_t  (ic,isubs) = t_snow_ex   (jc,isubs)
          t_s_now_t     (ic,isubs) = t_s_ex      (jc,isubs)
          t_g_t         (ic,isubs) = t_g_ex      (jc,isubs)
          qv_s_t        (ic,isubs) = qv_s_ex     (jc,isubs)
          w_snow_now_t  (ic,isubs) = w_snow_ex   (jc,isubs)
          rho_snow_now_t(ic,isubs) = rho_snow_ex (jc,isubs)
          w_i_now_t     (ic,isubs) = w_i_ex      (jc,isubs)
          IF (itype_interception == 2) THEN
            w_p_now_t   (ic,isubs) = w_p_ex      (jc,isubs)
            w_s_now_t   (ic,isubs) = w_s_ex      (jc,isubs)
          ELSE
            w_p_now_t   (ic,isubs) = 0._wp
            w_s_now_t   (ic,isubs) = 0._wp
          END IF
          freshsnow_t   (ic,isubs) = freshsnow_ex(jc,isubs)
          snowfrac_lc_t (ic,isubs) = snowfrac_lc_ex (jc,isubs)
          snowfrac_t    (ic,isubs) = snowfrac_ex (jc,isubs)
          runoff_s_t    (ic,isubs) = runoff_s_ex (jc,isubs)
          runoff_g_t    (ic,isubs) = runoff_g_ex (jc,isubs)

!         t_2m_t        (ic,isubs) = t_2m_ex     (jc)
          u_10m_t       (ic,isubs) = u_10m_ex    (jc)
          v_10m_t       (ic,isubs) = v_10m_ex    (jc)
          tch_t         (ic,isubs) = tch_ex      (jc,isubs)
          tcm_t         (ic,isubs) = tcm_ex      (jc,isubs)
          tfv_t         (ic,isubs) = tfv_ex      (jc,isubs)
          sobs_t        (ic,isubs) = sobs_ex     (jc,isubs)
          thbs_t        (ic,isubs) = thbs_ex     (jc,isubs)
          pabs_t        (ic,isubs) = pabs_ex     (jc,isubs)

          soiltyp_t     (ic,isubs) =  ext_data%atm%soiltyp_t(jc,jb,isubs)
          plcov_t       (ic,isubs) =  ext_data%atm%plcov_t  (jc,jb,isubs)
          rootdp_t      (ic,isubs) =  ext_data%atm%rootdp_t (jc,jb,isubs)
          sai_t         (ic,isubs) =  ext_data%atm%sai_t    (jc,jb,isubs)
          tai_t         (ic,isubs) =  ext_data%atm%tai_t    (jc,jb,isubs)
          eai_t         (ic,isubs) =  ext_data%atm%eai_t    (jc,jb,isubs)
          rsmin2d_t     (ic,isubs) =  ext_data%atm%rsmin2d_t(jc,jb,isubs)

          t_so_now_t(ic,nlev_soil+1,isubs) = t_so_ex(jc,nlev_soil+1,isubs)

          IF(lmulti_snow) THEN
            t_snow_mult_now_t(ic,nlev_snow+1,isubs) = t_snow_mult_ex(jc,nlev_snow+1,isubs)
            h_snow_t(ic,isubs)                      = h_snow_ex(jc,isubs)
          ENDIF
        ENDDO

       MSNOWI: IF(lmulti_snow) THEN

#ifdef __LOOP_EXCHANGE
        DO ic = 1, i_count
          jc = ext_data%atm%idx_lst_t(ic,jb,isubs)
          DO jk=1,nlev_snow
#else
!CDIR UNROLL=nlsnow
        DO jk=1,nlev_snow
          DO ic = 1, i_count
            jc = ext_data%atm%idx_lst_t(ic,jb,isubs)
#endif
            t_snow_mult_now_t  (ic,jk,isubs) = t_snow_mult_ex  (jc,jk,isubs)
            rho_snow_mult_now_t(ic,jk,isubs) = rho_snow_mult_ex(jc,jk,isubs)
            wliq_snow_now_t    (ic,jk,isubs) = wliq_snow_ex    (jc,jk,isubs)
            wtot_snow_now_t    (ic,jk,isubs) = wtot_snow_ex    (jc,jk,isubs)
            dzh_snow_now_t     (ic,jk,isubs) = dzh_snow_ex     (jc,jk,isubs)
          ENDDO
        ENDDO

       END IF MSNOWI

#ifdef __LOOP_EXCHANGE
        DO ic = 1, i_count
          jc = ext_data%atm%idx_lst_t(ic,jb,isubs)
          DO jk=1,nlev_soil
#else
!CDIR UNROLL=nlsoil
        DO jk=1,nlev_soil
          DO ic = 1, i_count
            jc = ext_data%atm%idx_lst_t(ic,jb,isubs)
#endif
            t_so_now_t    (ic,jk,isubs) = t_so_ex    (jc,jk,isubs)
            w_so_now_t    (ic,jk,isubs) = w_so_ex    (jc,jk,isubs)
            w_so_ice_now_t(ic,jk,isubs) = w_so_ice_ex(jc,jk,isubs)
          ENDDO
        ENDDO
!
!---------- END Copy index list fields


IF ( .true. ) THEN

        CALL terra_multlay(                                    &
!MR:<
        &  icant=2                                           , & !IN canopy-type (as set before, but maybe "1" is adequate?)
!MR:>
        &  ie=nproma                                         , & !IN array dimensions
        &  istartpar=1,       iendpar=i_count                , & !IN optional start/end indicies
        &  ke_soil=nlev_soil-1, ke_snow=nlev_snow            , & !IN without lowermost (climat.) soil layer
        &  ke_soil_hy   = ibot_w_so                          , & !IN number of hydrological active soil layers
        &  czmls=zml_soil,    ldiag_tg=.FALSE.               , & !IN processing soil level structure
        &  inwp_turb    = atm_phy_nwp_config(jg)%inwp_turb   , & !IN !!! Dangerous HACK !!!
        &  nclass_gscp  = atm_phy_nwp_config(jg)%nclass_gscp , & !IN number of hydrometeor classes
        &  dt=tcall_sfc_jg                                   , & !IN
        &  soiltyp_subs = soiltyp_t(:,isubs)                 , & ! type of the soil (keys 0-9)         --
        &  plcov        = plcov_t(:,isubs)                   , & ! fraction of plant cover             --
        &  rootdp       = rootdp_t(:,isubs)                  , & ! depth of the roots                ( m  )
        &  sai          = sai_t(:,isubs)                     , & ! surface area index                  --
        &  tai          = tai_t(:,isubs)                     , & ! surface area index                  --
        &  eai          = eai_t(:,isubs)                     , & ! surface area index                  --
        &  rsmin2d      = rsmin2d_t(:,isubs)                 , & ! minimum stomata resistance        ( s/m )
!
        &  u  = u_t(:)                                       , & ! zonal wind speed
        &  v  = v_t(:)                                       , & ! meridional wind speed
        &  t  = t_t(:)                                       , & ! temperature                       (  K  )
        &  qv = qv_t(:)                                      , & ! specific water vapor content      (kg/kg)
        &  p0 = p0_t(:)                                      , & ! base state pressure               ( Pa  )
        &  ps = ps_t(:)                                      , & ! surface pressure                  ( Pa  )
!
        &  t_snow_now    = t_snow_now_t(:,isubs)             , & ! temperature of the snow-surface   (  K  )
        &  t_snow_new    = t_snow_new_t(:,isubs)             , & ! temperature of the snow-surface   (  K  )
!
        &  t_snow_mult_now = t_snow_mult_now_t(:,:,isubs)    , & ! temperature of the snow-surface   (  K  )
        &  t_snow_mult_new = t_snow_mult_new_t(:,:,isubs)    , & ! temperature of the snow-surface   (  K  )
!
        &  t_s_now       = t_s_now_t(:,isubs)                , & ! temperature of the ground surface (  K  )
        &  t_s_new       = t_s_new_t(:,isubs)                , & ! temperature of the ground surface (  K  )
!
        &  t_g           =  t_g_t (:,isubs)                  , & ! weighted surface temperature      (  K  )
        &  qv_s          =  qv_s_t(:,isubs)                  , & ! specific humidity at the surface  (kg/kg)
!
        &  w_snow_now    = w_snow_now_t(:,isubs)             , & ! water content of snow             (m H2O)
        &  w_snow_new    = w_snow_new_t(:,isubs)             , & ! water content of snow             (m H2O)
!
        &  rho_snow_now  = rho_snow_now_t(:,isubs)           , & ! snow density                      (kg/m**3)
        &  rho_snow_new  = rho_snow_new_t(:,isubs)           , & ! snow density                      (kg/m**3)
!
        &  rho_snow_mult_now = rho_snow_mult_now_t(:,:,isubs), & ! snow density                      (kg/m**3)
        &  rho_snow_mult_new = rho_snow_mult_new_t(:,:,isubs), & ! snow density                      (kg/m**3)
!
        &  h_snow        = h_snow_t(:,isubs)                 , & ! snow height
        &  h_snow_gp     = h_snow_t(:,isubs)                 , & ! snow height
        &  meltrate      = meltrate(:)                       , & ! snow melting rate
!
        &  w_i_now       = w_i_now_t(:,isubs)                , & ! water content of interception water (m H2O)
        &  w_i_new       = w_i_new_t(:,isubs)                , & ! water content of interception water (m H2O)
!
        &  w_p_now       = w_p_now_t(:,isubs)                , & ! water content of interception water (m H2O)
        &  w_p_new       = w_p_new_t(:,isubs)                , & ! water content of interception water (m H2O)
!
        &  w_s_now       = w_s_now_t(:,isubs)                , & ! water content of interception water (m H2O)
        &  w_s_new       = w_s_new_t(:,isubs)                , & ! water content of interception water (m H2O)
!
        &  t_so_now      = t_so_now_t(:,:,isubs)             , & ! soil temperature (main level)     (  K  )
        &  t_so_new      = t_so_new_t(:,:,isubs)             , & ! soil temperature (main level)     (  K  )
!
        &  w_so_now      = w_so_now_t(:,:,isubs)             , & ! total water content (ice + liquid water) (m H20)
        &  w_so_new      = w_so_new_t(:,:,isubs)             , & ! total water content (ice + liquid water) (m H20)
!
        &  w_so_ice_now  = w_so_ice_now_t(:,:,isubs)         , & ! ice content                       (m H20)
        &  w_so_ice_new  = w_so_ice_new_t(:,:,isubs)         , & ! ice content                       (m H20)
!
        &  u_10m         = u_10m_t(:,isubs)                  , & ! zonal wind in 10m                  ( m/s )
        &  v_10m         = v_10m_t(:,isubs)                  , & ! meridional wind in 10m            ( m/s )
        &  freshsnow     = freshsnow_t(:,isubs)              , & ! indicator for age of snow in top of snow layer (  -  )
        &  zf_snow       = snowfrac_t(:,isubs)               , & ! snow-cover fraction                            (  -  )
!
        &  wliq_snow_now = wliq_snow_now_t(:,:,isubs)        , & ! liquid water content in the snow  (m H2O)
        &  wliq_snow_new = wliq_snow_new_t(:,:,isubs)        , & ! liquid water content in the snow  (m H2O)
!
        &  wtot_snow_now = wtot_snow_now_t(:,:,isubs)        , & ! total (liquid + solid) water content of snow (m H2O)
        &  wtot_snow_new = wtot_snow_new_t(:,:,isubs)        , & ! total (liquid + solid) water content of snow (m H2O)
!
        &  dzh_snow_now  = dzh_snow_now_t(:,:,isubs)         , & ! layer thickness between half levels in snow  (  m  )
        &  dzh_snow_new  = dzh_snow_new_t(:,:,isubs)         , & ! layer thickness between half levels in snow  (  m  )
!
        &  prr_con       = prr_con_t(:)                      , & ! precipitation rate of rain, convective       (kg/m2*s)
        &  prs_con       = prs_con_t(:)                      , & ! precipitation rate of snow, convective       (kg/m2*s)
        &  conv_frac     = conv_frac(:)                      , & ! convective area fraction
        &  prr_gsp       = prr_gsp_t(:)                      , & ! precipitation rate of rain, grid-scale       (kg/m2*s)
        &  prs_gsp       = prs_gsp_t(:)                      , & ! precipitation rate of snow, grid-scale       (kg/m2*s)
        &  prg_gsp       = dummy_prg_gsp(:)                  , & ! precipitation rate of graupel, grid-scale    (kg/m2*s)
!
        &  tch           = tch_t(:,isubs)                    , & ! turbulent transfer coefficient for heat      ( -- )
        &  tcm           = tcm_t(:,isubs)                    , & ! turbulent transfer coefficient for momentum  ( -- )
        &  tfv           = tfv_t(:,isubs)                    , & ! laminar reduction factor for evaporation     ( -- )
!
        &  sobs          = sobs_t(:,isubs)                   , & ! solar radiation at the ground                (W/m2)
        &  thbs          = thbs_t(:,isubs)                   , & ! thermal radiation at the ground              (W/m2)
        &  pabs          = pabs_t(:,isubs)                   , & ! photosynthetic active radiation              (W/m2)
!
        &  runoff_s      = runoff_s_t(:,isubs)               , & ! surface water runoff; sum over forecast      (kg/m2)
        &  runoff_g      = runoff_g_t(:,isubs)               , & ! soil water runoff; sum over forecast         (kg/m2)
!
        &  zshfl_s       = shfl_soil_t(:,isubs)              , & ! sensible heat flux soil/air interface         (W/m2)
        &  zlhfl_s       = lhfl_soil_t(:,isubs)              , & ! latent   heat flux soil/air interface         (W/m2)
        &  zshfl_snow    = shfl_snow_t(:,isubs)              , & ! sensible heat flux snow/air interface         (W/m2)
        &  zlhfl_snow    = lhfl_snow_t(:,isubs)              , & ! latent   heat flux snow/air interface         (W/m2)
        &  lhfl_bs       = lhfl_bs_t  (:,isubs)              , & ! out: latent heat flux from bare soil evap.    (W/m2)
        &  lhfl_pl       = lhfl_pl_t  (:,:,isubs)            , & ! out: latent heat flux from bare soil evap.    (W/m2)
        &  rstom         = rstom_t    (:,isubs)              , & ! out: stomatal resistance                     ( s/m )
        &  zshfl_sfc     = shfl_s_t   (:,isubs)              , & !OUT sensible heat flux surface interface     (W/m2)
        &  zlhfl_sfc     = lhfl_s_t   (:,isubs)              , & !OUT latent   heat flux surface interface     (W/m2)
        &  zqhfl_sfc     = qhfl_s_t   (:,isubs)                ) !OUT moisture flux surface interface          (kg/m2/s)

!DR aggregation of lhfl_bs_t is missing
!DR aggregation of lhfl_pl_t is missing
!DR aggregation of rstom_t is missing
endif


IF (msg_level >= 15) THEN
  DO ic = 1, i_count
    jc = ext_data%atm%idx_lst_t(ic,jb,isubs)
    if ( abs(shfl_soil_t(ic,isubs)) >  800.0 .or. shfl_snow_t(ic,isubs) >  800.0  .or. &
         abs(lhfl_soil_t(ic,isubs)) > 2000.0 .or. lhfl_snow_t(ic,isubs) > 2000.0 ) then
       write(*,*) 'sfc_interface_edmf: fluxes ', &
         ic, jc, isubs, snowfrac_t(ic,isubs), &
         shfl_soil_t(ic,isubs), shfl_snow_t(ic,isubs), &
         lhfl_soil_t(ic,isubs), lhfl_snow_t(ic,isubs)
    endif
    DO jk=1,nlev_soil+1
      if ( abs( t_so_ex(jc,jk,isubs) - t_so_new_t(ic,jk,isubs) ) > 5.0 ) then
        write(*,*) 'sfc_interface_edmf: t_so ', ic, jc, isubs, jk, &
          t_so_ex(jc,jk,isubs), t_so_new_t(ic,jk,isubs), &
          lc_class_t(ic,isubs), snowfrac_t(ic,isubs), &
          shfl_soil_t(ic,isubs), shfl_snow_t(ic,isubs)
      endif
    ENDDO
    DO jk=1,nlev_soil
      if ( abs( w_so_ex(jc,jk,isubs) - w_so_new_t(ic,jk,isubs) ) > 0.5 ) then
        write(*,*) 'sfc_interface_edmf: w_so ', ic, jc, isubs, jk, &
          w_so_ex(jc,jk,isubs), w_so_new_t(ic,jk,isubs)
      endif
      if ( abs( w_so_ice_ex(jc,jk,isubs) - w_so_ice_new_t(ic,jk,isubs) ) > 0.5 ) then
        write(*,*) 'sfc_interface_edmf: w_so_ice ', ic, jc, isubs, jk, &
          w_so_ice_ex(jc,jk,isubs), w_so_ice_new_t(ic,jk,isubs)
      endif
    ENDDO
  ENDDO
ENDIF



if (.true.) then
        CALL diag_snowfrac_tg(                        &
          &  istart = 1, iend = i_count             , & ! start/end indices
          &  lc_class  = lc_class_t        (:,isubs), & ! land-cover class
          &  i_lc_urban = ext_data%atm%i_lc_urban   , & ! land-cover class index for urban areas
          &  t_snow    = t_snow_new_t      (:,isubs), & ! snow temp
          &  t_soiltop = t_s_new_t         (:,isubs), & ! soil top temp
          &  w_snow    = w_snow_new_t      (:,isubs), & ! snow WE
          &  rho_snow  = rho_snow_new_t    (:,isubs), & ! snow depth
          &  freshsnow = freshsnow_t       (:,isubs), & ! fresh snow fraction
          &  meltrate  = meltrate          (:),       & ! snow melting rate
          &  sso_sigma = sso_sigma_t       (:),       & ! sso stdev
          &  tai       = tai_t             (:,isubs), & ! effective leaf area index
          &  snowfrac  = snowfrac_t        (:,isubs), & ! OUT: snow cover fraction
          &  t_g       = t_g_t             (:,isubs)  ) ! OUT: averaged ground temp
endif


!---------- Copy index list fields back to state fields

!CDIR NODEP,VOVERTAKE,VOB
        DO ic = 1, i_count
          jc = ext_data%atm%idx_lst_t(ic,jb,isubs)
          t_snow_ex   (jc,isubs)  = t_snow_new_t  (ic,isubs)
          t_s_ex      (jc,isubs)  = t_s_new_t     (ic,isubs)
          t_g_ex      (jc,isubs)  = t_g_t         (ic,isubs)
          qv_s_ex     (jc,isubs)  = qv_s_t        (ic,isubs)
          w_snow_ex   (jc,isubs)  = w_snow_new_t  (ic,isubs)
          rho_snow_ex (jc,isubs)  = rho_snow_new_t(ic,isubs)
          h_snow_ex   (jc,isubs)  = h_snow_t      (ic,isubs)
          w_i_ex      (jc,isubs)  = w_i_new_t     (ic,isubs)
          IF (itype_interception == 2) THEN
            w_p_ex    (jc,isubs)  = w_p_new_t     (ic,isubs)
            w_s_ex    (jc,isubs)  = w_s_new_t     (ic,isubs)
          END IF
          freshsnow_ex(jc,isubs)  = freshsnow_t   (ic,isubs)
          ! Remark: the two snow-cover fraction variables differ only if lsnowtile=true (see below)
          snowfrac_lc_ex(jc,isubs)= snowfrac_t    (ic,isubs)
          snowfrac_ex (jc,isubs)  = snowfrac_t    (ic,isubs)
          runoff_s_ex (jc,isubs)  = runoff_s_t    (ic,isubs)
          runoff_g_ex (jc,isubs)  = runoff_g_t    (ic,isubs)

          t_so_ex(jc,nlev_soil+1,isubs) = t_so_new_t(ic,nlev_soil+1,isubs)

          lhfl_bs_ex  (jc,isubs)  = lhfl_bs_t     (ic,isubs)
          rstom_ex    (jc,isubs)  = rstom_t       (ic,isubs)

          shfl_s_ex   (jc,isubs)  = shfl_s_t      (ic,isubs)
          lhfl_s_ex   (jc,isubs)  = lhfl_s_t      (ic,isubs)
          qhfl_s_ex   (jc,isubs)  = qhfl_s_t      (ic,isubs)
          shfl_soil_ex(jc,isubs)  = shfl_soil_t   (ic,isubs)
          lhfl_soil_ex(jc,isubs)  = lhfl_soil_t   (ic,isubs)
          shfl_snow_ex(jc,isubs)  = shfl_snow_t   (ic,isubs)
          lhfl_snow_ex(jc,isubs)  = lhfl_snow_t   (ic,isubs)


          IF(lmulti_snow) THEN
            t_snow_mult_ex     (jc,nlev_snow+1,isubs) = &
              t_snow_mult_new_t(ic,nlev_snow+1,isubs)
          ENDIF
        ENDDO

        IF (lsnowtile .AND. isubs > ntiles_lnd) THEN ! copy snowfrac_t to snow-free tile
!CDIR NODEP,VOVERTAKE,VOB                            ! (needed for index list computation)
          DO ic = 1, i_count
            jc = ext_data%atm%idx_lst_t(ic,jb,isubs)
            snowfrac_lc_ex(jc,isubs-ntiles_lnd) = snowfrac_lc_ex(jc,isubs)
          ENDDO
        ENDIF

        MSNOWO: IF(lmulti_snow) THEN

#ifdef __LOOP_EXCHANGE
        DO ic = 1, i_count
          jc = ext_data%atm%idx_lst_t(ic,jb,isubs)
          DO jk=1,nlev_snow
#else
!CDIR UNROLL=nlsnow
        DO jk=1,nlev_snow
!CDIR NODEP,VOVERTAKE,VOB
          DO ic = 1, i_count
            jc = ext_data%atm%idx_lst_t(ic,jb,isubs)
#endif
            t_snow_mult_ex  (jc,jk,isubs) = t_snow_mult_new_t  (ic,jk,isubs)
            rho_snow_mult_ex(jc,jk,isubs) = rho_snow_mult_new_t(ic,jk,isubs)
            wliq_snow_ex    (jc,jk,isubs) = wliq_snow_new_t    (ic,jk,isubs)
            wtot_snow_ex    (jc,jk,isubs) = wtot_snow_new_t    (ic,jk,isubs)
            dzh_snow_ex     (jc,jk,isubs) = dzh_snow_new_t     (ic,jk,isubs)
          ENDDO
        ENDDO

        END IF MSNOWO


#ifdef __LOOP_EXCHANGE
        DO ic = 1, i_count
          jc = ext_data%atm%idx_lst_t(ic,jb,isubs)
          DO jk=1,nlev_soil
#else
!CDIR UNROLL=nlsoil
        DO jk=1,nlev_soil
!CDIR NODEP,VOVERTAKE,VOB
          DO ic = 1, i_count
            jc = ext_data%atm%idx_lst_t(ic,jb,isubs)
#endif
            t_so_ex    (jc,jk,isubs) = t_so_new_t    (ic,jk,isubs)
            w_so_ex    (jc,jk,isubs) = w_so_new_t    (ic,jk,isubs)
            w_so_ice_ex(jc,jk,isubs) = w_so_ice_new_t(ic,jk,isubs)

            ! diagnostic field
            lhfl_pl_ex (jc,jk,isubs) = lhfl_pl_t     (ic,jk,isubs)
          ENDDO
        ENDDO

       END DO ! isubs - loop over tiles

       IF(lsnowtile) THEN      ! snow is considered as separate tiles
         DO isubs = 1, ntiles_lnd

           isubs_snow = isubs + ntiles_lnd

           ! save previous area fractions for subsequent redistribution computations
           frac_sv(:)      = ext_data%atm%frac_t(:,jb,isubs)
           frac_snow_sv(:) = ext_data%atm%frac_t(:,jb,isubs_snow)

           ! update index lists for snow tiles
           CALL update_idx_lists_lnd (idx_lst_lp       = ext_data%atm%idx_lst_lp_t(:,jb,isubs),         &
                                    lp_count           = ext_data%atm%lp_count_t(jb,isubs),             &
                                    idx_lst            = ext_data%atm%idx_lst_t(:,jb,isubs),            &
                                    gp_count           = ext_data%atm%gp_count_t(jb,isubs),             &
                                    idx_lst_snow       = ext_data%atm%idx_lst_t(:,jb,isubs_snow),       &
                                    gp_count_snow      = ext_data%atm%gp_count_t(jb,isubs_snow),        &
                                    lc_frac            = ext_data%atm%lc_frac_t(:,jb,isubs),            &
                                    partial_frac       = ext_data%atm%frac_t(:,jb,isubs),               &
                                    partial_frac_snow  = ext_data%atm%frac_t(:,jb,isubs_snow),          &
                                    snowtile_flag      = ext_data%atm%snowtile_flag_t(:,jb,isubs),      &
                                    snowtile_flag_snow = ext_data%atm%snowtile_flag_t(:,jb,isubs_snow), &
                                    snowfrac           = snowfrac_lc_ex(:,isubs)                      )

           i_count = ext_data%atm%gp_count_t(jb,isubs)
           i_count_snow = ext_data%atm%gp_count_t(jb,isubs_snow)

           ! Check for newly activated grid points that need to be initialized
           icount_init = 0
           DO ic = 1, i_count
             jc = ext_data%atm%idx_lst_t(ic,jb,isubs)
             IF (ext_data%atm%snowtile_flag_t(jc,jb,isubs) == 2) THEN
               icount_init = icount_init + 1
               init_list(icount_init) = jc
               it1(icount_init) = isubs      ! target of copy operation
               it2(icount_init) = isubs_snow ! source of copy operation
             ENDIF
           ENDDO
           DO ic = 1, i_count_snow
             jc = ext_data%atm%idx_lst_t(ic,jb,isubs_snow)
             IF (ext_data%atm%snowtile_flag_t(jc,jb,isubs_snow) == 2) THEN
               icount_init = icount_init + 1
               init_list(icount_init) = jc
               it1(icount_init) = isubs_snow ! target of copy operation
               it2(icount_init) = isubs      ! source of copy operation
             ENDIF
           ENDDO

           DO ic = 1, icount_init
             jc = init_list(ic)
             is1 = it1(ic)
             is2 = it2(ic)
             t_snow_ex      (jc,is1) = t_snow_ex      (jc,is2)
             t_s_ex         (jc,is1) = t_s_ex         (jc,is2)
             t_g_ex         (jc,is1) = t_g_ex         (jc,is2)
             qv_s_ex        (jc,is1) = qv_s_ex        (jc,is2)
             w_snow_ex      (jc,is1) = w_snow_ex      (jc,is2)
             rho_snow_ex    (jc,is1) = rho_snow_ex    (jc,is2)
             h_snow_ex      (jc,is1) = h_snow_ex      (jc,is2)
             w_i_ex         (jc,is1) = w_i_ex         (jc,is2)
             w_p_ex         (jc,is1) = w_p_ex         (jc,is2)
             w_s_ex         (jc,is1) = w_s_ex         (jc,is2)
             freshsnow_ex   (jc,is1) = freshsnow_ex   (jc,is2)
             snowfrac_lc_ex (jc,is1) = snowfrac_lc_ex (jc,is2)
             snowfrac_ex    (jc,is1) = snowfrac_ex    (jc,is2)
             runoff_s_ex    (jc,is1) = runoff_s_ex    (jc,is2)
             runoff_g_ex    (jc,is1) = runoff_g_ex    (jc,is2)

             t_so_ex    (jc,:,is1) = t_so_ex    (jc,:,is2)
             w_so_ex    (jc,:,is1) = w_so_ex    (jc,:,is2)
             w_so_ice_ex(jc,:,is1) = w_so_ice_ex(jc,:,is2)

             IF (lmulti_snow) THEN
               t_snow_mult_ex  (jc,:,is1) = t_snow_mult_ex  (jc,:,is2)
               rho_snow_mult_ex(jc,:,is1) = rho_snow_mult_ex(jc,:,is2)
               wliq_snow_ex    (jc,:,is1) = wliq_snow_ex    (jc,:,is2)
               wtot_snow_ex    (jc,:,is1) = wtot_snow_ex    (jc,:,is2)
               dzh_snow_ex     (jc,:,is1) = dzh_snow_ex     (jc,:,is2)
             ENDIF
           ENDDO

!CDIR NODEP,VOVERTAKE,VOB
           DO ic = 1, i_count_snow
             jc = ext_data%atm%idx_lst_t(ic,jb,isubs_snow)

             IF (ext_data%atm%snowtile_flag_t(jc,jb,isubs_snow) == 1 .AND. &
                 ext_data%atm%snowtile_flag_t(jc,jb,isubs)      == 1) THEN

               ! compute factors for redistribution of heat and moisture
               fact1(jc) = MIN(1._wp,frac_sv(jc)/     MAX(small,ext_data%atm%frac_t(jc,jb,isubs)     ))
               fact2(jc) = MIN(1._wp,frac_snow_sv(jc)/MAX(small,ext_data%atm%frac_t(jc,jb,isubs_snow)))
             ENDIF

           END DO

           ! redistribution of heat and moisture between snow-covered and snow-free tiles
           ! according to their new fractions, in order to keep heat and moisture balances
           DO jk = 1, nlev_soil
!CDIR NODEP,VOVERTAKE,VOB
             DO ic = 1, i_count_snow
               jc = ext_data%atm%idx_lst_t(ic,jb,isubs_snow)

               IF (ext_data%atm%snowtile_flag_t(jc,jb,isubs_snow) == 1 .AND. &
                   ext_data%atm%snowtile_flag_t(jc,jb,isubs)      == 1) THEN

                 tmp1 = t_so_ex(jc,jk,isubs)
                 tmp2 = w_so_ex(jc,jk,isubs)
                 tmp3 = w_so_ice_ex(jc,jk,isubs)

                 t_so_ex    (jc,jk,isubs) = t_so_ex    (jc,jk,isubs)*fact1(jc) &
                   &                      + t_so_ex    (jc,jk,isubs_snow)*(1._wp - fact1(jc))
                 w_so_ex    (jc,jk,isubs) = w_so_ex    (jc,jk,isubs)*fact1(jc) &
                   &                      + w_so_ex    (jc,jk,isubs_snow)*(1._wp - fact1(jc))
                 w_so_ice_ex(jc,jk,isubs) = w_so_ice_ex(jc,jk,isubs)*fact1(jc) &
                   &                      + w_so_ice_ex(jc,jk,isubs_snow)*(1._wp - fact1(jc))

                 t_so_ex    (jc,jk,isubs_snow) = tmp1*(1._wp - fact2(jc)) &
                   &                           + t_so_ex    (jc,jk,isubs_snow)*fact2(jc)
                 w_so_ex    (jc,jk,isubs_snow) = tmp2*(1._wp - fact2(jc)) &
                   &                           + w_so_ex    (jc,jk,isubs_snow)*fact2(jc)
                 w_so_ice_ex(jc,jk,isubs_snow) = tmp3*(1._wp - fact2(jc)) &
                   &                           + w_so_ice_ex(jc,jk,isubs_snow)*fact2(jc)

                 IF (jk == 1) THEN
                   t_s_ex(jc,isubs)      = t_so_ex(jc,jk,isubs)
                   t_s_ex(jc,isubs_snow) = t_so_ex(jc,jk,isubs_snow)
                 ENDIF
               ENDIF

             END DO
           END DO        ! soil layers
!CDIR NODEP,VOVERTAKE,VOB
           DO ic = 1, i_count_snow
             jc = ext_data%atm%idx_lst_t(ic,jb,isubs_snow)

             ! reset field for actual snow-cover for grid points / land-cover classes for which there
             ! are seperate snow-free and snow-covered tiles
             snowfrac_ex(jc,isubs)  = 0._wp
             w_snow_ex(jc,isubs)    = 0._wp
             t_snow_ex(jc,isubs)    = t_s_ex(jc,isubs)
             t_g_ex(jc,isubs)       = t_s_ex(jc,isubs)

             ! to prevent numerical stability problems, we require at least 1 cm of snow in order to
             ! have a snow-cover fraction of 1 on snow tiles (not critical for the single-layer
             ! snow scheme, but the multi-layer snow model becomes numerically unstable within a few
             ! time steps when associating traces of snow with a snow-cover fraction of 1)
             snowfrac_ex(jc,isubs_snow) = MIN(1._wp,h_snow_ex(jc,isubs_snow)/0.01_wp)

             ! Rediagnose t_g according to the modified snow-cover fraction
             t_g_ex(jc,isubs_snow) =  &
               snowfrac_ex(jc,isubs_snow) * t_snow_ex(jc,isubs_snow) + &
               (1._wp-snowfrac_ex(jc,isubs_snow))*t_s_ex(jc,isubs_snow)

             IF (lmulti_snow) THEN
               t_snow_mult_ex(jc,nlev_snow+1,isubs) = t_s_ex(jc,isubs)
             ENDIF
           END DO

           IF (lmulti_snow) THEN
!CDIR UNROLL=nlsnow
             DO jk=1,nlev_snow
!CDIR NODEP,VOVERTAKE,VOB
               DO ic = 1, i_count_snow
                 jc = ext_data%atm%idx_lst_t(ic,jb,isubs_snow)
                 t_snow_mult_ex(jc,jk,isubs) = t_s_ex(jc,isubs)
                 wliq_snow_ex(jc,jk,isubs)   = 0._wp
                 wtot_snow_ex(jc,jk,isubs)   = 0._wp
                 dzh_snow_ex (jc,jk,isubs)   = 0._wp
               ENDDO
             ENDDO
           ENDIF

         END DO

       ENDIF  !snow tiles

    ELSE IF ( atm_phy_nwp_config(jg)%inwp_surface == 2 ) THEN

          !-------------------------------------------------------------------------
          !> ECHAM version
          !-------------------------------------------------------------------------



    ENDIF !inwp_sfc



    !
    ! Call seaice parameterization
    !
    IF ( (atm_phy_nwp_config(jg)%inwp_surface == 1) .AND. (lseaice) ) THEN

      CALL nwp_seaice(ext_data    , jb          , tcall_sfc_jg,                      &
                   &  t_g_ex      , t_s_ex      , qv_s_ex     , ps_ex    , sobs_ex,  &
                   &  thbs_ex,                                                       &
                   &  shfl_soil_ex, lhfl_soil_ex,                                    &
                   &  prr_con_ex  , prs_con_ex  ,  prr_gsp_ex , prs_gsp_ex  ,        &
                   &  t_ice       , h_ice       , t_snow_si   , h_snow_si,           &
                   &  alb_si      ,                                                  &
                   &  fr_seaice                                                      )
    ENDIF



    !
    ! Call Flake model
    !
!    IF ( (atm_phy_nwp_config(jg)%inwp_surface == 1) .AND. (llake) ) THEN
!    ENDIF



    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Final step: aggregate t_g and qv_s            !!
    !                                               !!
    ! Loop over all points (land AND water points)  !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !aggregation of shfl_s and lhfs_s done in turb_sfc_interface

    IF (ntiles_total == 1) THEN
      DO jc = i_startidx, i_endidx
        t_g (jc) = t_g_ex (jc,1)
        qv_s(jc) = qv_s_ex(jc,1)
      ENDDO
    ELSE ! aggregate fields over tiles
      t_g_s (:) =  0._wp
      qv_s_s(:) =  0._wp
      DO isubs = 1,ntiles_total+ntiles_water
        DO jc = i_startidx, i_endidx
          t_g_s(jc)  = t_g_s(jc)  + ext_data%atm%frac_t(jc,jb,isubs)* &
                       t_g_ex(jc,isubs)**4
          qv_s_s(jc) = qv_s_s(jc) + ext_data%atm%frac_t(jc,jb,isubs)* &
                       qv_s_ex(jc,isubs)
        ENDDO
      ENDDO

      DO jc = i_startidx, i_endidx
        t_g(jc)  = SQRT(SQRT(t_g_s(jc)))
        qv_s(jc) = qv_s_s(jc)
      ENDDO

    ENDIF     ! with or without tiles


  END SUBROUTINE nwp_surface_edmf



  !>
  !! Interface for seaice parameterization
  !!
  !! Interface for seaice parameterization. Calls seaice time integration scheme
  !! seaice_timestep_nwp and updates the dynamic seaice index lists.
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2012-08-31)
  !!
  SUBROUTINE nwp_seaice (ext_data    , jb          , dtime       ,                         &
                      &  t_g_ex      , t_s_ex      , qv_s_ex     , ps_ex       , sobs_ex,  &
                      &  thbs_ex,                                                          &
                      &  shfl_soil_ex, lhfl_soil_ex,                                       &
                      &  prr_con_ex  , prs_con_ex  ,  prr_gsp_ex , prs_gsp_ex  ,           &
                      &  t_ice_ex    , h_ice_ex    , t_snow_si_ex, h_snow_si_ex,           &
                      &  alb_si_ex   ,                                                     &
                      &  fr_seaice                                                         )

    INTEGER,                                               INTENT(IN)    :: &
                  jb                   ! block
    TYPE(t_external_data),                                 INTENT(INOUT) :: &
                  ext_data             ! external data
    REAL(wp),                                              INTENT(IN)    :: &
                  dtime                ! time interval for sea ice
    REAL(wp), DIMENSION(nproma,ntiles_total+ntiles_water), INTENT(INOUT) :: &
                  t_g_ex           , & ! weighted surface temperature                  (  K  )
                  t_s_ex           , & ! weighted surface temperature                  (  K  )
                  qv_s_ex              ! specific humidity at the surface              (kg/kg)
    REAL(wp), DIMENSION(nproma),                           INTENT(IN)    :: &
                  ps_ex                ! surface pressure                              ( pa  )
    REAL(wp), DIMENSION(nproma,ntiles_total+ntiles_water), INTENT(IN)    :: &
                  sobs_ex          , & ! solar radiation at the ground                 ( W/m2)
                  thbs_ex              ! thermal radiation at the ground               ( W/m2)
    REAL(wp), DIMENSION(nproma,ntiles_total+ntiles_water), INTENT(IN)    :: &
                  shfl_soil_ex     , & ! sensible heat flux soil/air interface         ( W/m2)
                  lhfl_soil_ex         ! latent   heat flux soil/air interface         ( W/m2)
    REAL(wp), DIMENSION(nproma),                           INTENT(IN)    :: &
                  prr_con_ex       , & ! precipitation rate of rain, convective        (kg/m2*s)
                  prs_con_ex       , & ! precipitation rate of snow, convective        (kg/m2*s)
                  prr_gsp_ex       , & ! precipitation rate of rain, grid-scale        (kg/m2*s)
                  prs_gsp_ex           ! precipitation rate of snow, grid-scale        (kg/m2*s)

    REAL(wp), DIMENSION(nproma),                           INTENT(INOUT) :: &
                  t_ice_ex         , & ! sea ice temperature                           (  K  )
                  h_ice_ex         , & ! sea ice height                                (  m  )
                  t_snow_si_ex     , & ! sea ice snow temperature                      (  K  )
                  h_snow_si_ex     , & ! sea ice snow height                           (  m  )
                  alb_si_ex            ! sea-ice albedo                                (  -  )
    REAL(wp), DIMENSION(nproma),                           INTENT(INOUT) :: &
                  fr_seaice            ! sea ice fraction                              (  1  )

    ! Local arrays  (local copies)
    !
    REAL(wp) :: shfl_s   (nproma)   ! sensible heat flux at the surface                [W/m^2]
    REAL(wp) :: lhfl_s   (nproma)   ! latent heat flux at the surface                  [W/m^2]
    REAL(wp) :: lwflxsfc (nproma)   ! net long-wave radiation flux at the surface      [W/m^2]
    REAL(wp) :: swflxsfc (nproma)   ! net solar radiation flux at the surface          [W/m^2]
    REAL(wp) :: snow_rate(nproma)   ! snow rate (convecive + grid-scale)               [kg/(m^2 s)]
    REAL(wp) :: rain_rate(nproma)   ! rain rate (convecive + grid-scale)               [kg/(m^2 s)]
    REAL(wp) :: tice_now (nproma)   ! temperature of ice upper surface at previous time  [K]
    REAL(wp) :: hice_now (nproma)   ! ice thickness at previous time level               [m]
    REAL(wp) :: tsnow_now(nproma)   ! temperature of snow upper surface at previous time [K]
    REAL(wp) :: hsnow_now(nproma)   ! snow thickness at previous time level              [m]
    REAL(wp) :: albsi_now(nproma)   ! sea-ice albedo at previous time level              [-]
    REAL(wp) :: tice_new (nproma)   ! temperature of ice upper surface at new time       [K]
    REAL(wp) :: hice_new (nproma)   ! ice thickness at new time level                    [m]
    REAL(wp) :: tsnow_new(nproma)   ! temperature of snow upper surface at new time      [K]
    REAL(wp) :: hsnow_new(nproma)   ! snow thickness at new time level                   [m]
    REAL(wp) :: albsi_new(nproma)   ! sea-ice albedo at new time level                   [-]

    REAL(wp) :: t_s_dummy(nproma)   ! dummy for surface temperature
    REAL(wp) :: t_seasfc_dummy(nproma) ! dummy for sea surface temperature

    ! Local scalars:
    !
    INTEGER :: jc, ic               !loop indices
    INTEGER :: i_count

    CHARACTER(len=*), PARAMETER :: routine = 'mo_nwp_sfc_interface:nwp_seaice'
    !-------------------------------------------------------------------------


    IF (msg_level >= 15) THEN
      CALL message('mo_nwp_sfc_interface: ', 'call nwp_seaice scheme')
    ENDIF


      ! Copy input fields
      !
      i_count = ext_data%atm%spi_count(jb)

      DO ic = 1, i_count
        jc = ext_data%atm%idx_lst_spi(ic,jb)

        shfl_s   (ic) = shfl_soil_ex(jc,isub_seaice)     ! sensible heat flux at sfc          [W/m^2]
        lhfl_s   (ic) = lhfl_soil_ex(jc,isub_seaice)     ! latent heat flux at sfc            [W/m^2]
        lwflxsfc (ic) = thbs_ex     (jc,isub_seaice)     ! net lw radiation flux at sfc       [W/m^2]
        swflxsfc (ic) = sobs_ex     (jc,isub_seaice)     ! net solar radiation flux at sfc    [W/m^2]
        snow_rate(ic) = prs_gsp_ex(jc) + prs_con_ex(jc)  ! snow rate (convecive + grid-scale) [kg/(m^2 s)]
        rain_rate(ic) = prr_gsp_ex(jc) + prr_con_ex(jc)  ! rain rate (convecive + grid-scale) [kg/(m^2 s)]
        tice_now (ic) = t_ice_ex    (jc)
        hice_now (ic) = h_ice_ex    (jc)
        tsnow_now(ic) = t_snow_si_ex(jc)
        hsnow_now(ic) = h_snow_si_ex(jc)
        albsi_now(ic) = alb_si_ex   (jc)              ! sea-ice albedo at previous time level [-]
      ENDDO  ! ic

      ! call seaice time integration scheme
      !
      CALL seaice_timestep_nwp (       &
        &   dtime   = dtime,           &
        &   nsigb   = i_count,         & !in
        &   qsen    = shfl_s   (:),    & !in
        &   qlat    = lhfl_s   (:),    & !in
        &   qlwrnet = lwflxsfc (:),    & !in
        &   qsolnet = swflxsfc (:),    & !in
        &   snow_rate = snow_rate(:),  & !in
        &   rain_rate = rain_rate(:),  & !in
        &   tice_p  = tice_now (:),    & !in
        &   hice_p  = hice_now (:),    & !in
        &   tsnow_p = tsnow_now(:),    & !in    ! DUMMY: not used yet
        &   hsnow_p = hsnow_now(:),    & !in    ! DUMMY: not used yet
        &   albsi_p = albsi_now(:),    & !in   
        &   tice_n  = tice_new (:),    & !out
        &   hice_n  = hice_new (:),    & !out
        &   tsnow_n = tsnow_new(:),    & !out   ! DUMMY: not used yet
        &   hsnow_n = hsnow_new(:),    & !out   ! DUMMY: not used yet
        &   albsi_n = albsi_new(:)     ) !out  
      ! optional arguments dticedt, dhicedt, dtsnowdt, dhsnowdt (tendencies) are neglected


      !  Recover fields from index list
      !
      DO ic = 1, i_count
        jc = ext_data%atm%idx_lst_spi(ic,jb)

!debug
        if (tice_new(ic) < 100.0 .or. t_g_ex(jc,isub_seaice) < 100.0 ) then
          write(*,*) 'seaice1: ', tice_new(ic), tice_now(ic), hice_new(ic), &
            & hice_now(ic), shfl_s(ic), lhfl_s(ic), lwflxsfc(ic), swflxsfc(ic)
        endif
!xxxxx

        t_ice_ex    (jc) = tice_new (ic)
        h_ice_ex    (jc) = hice_new (ic)
        t_snow_si_ex(jc) = tsnow_new(ic)
        h_snow_si_ex(jc) = hsnow_new(ic)
        alb_si_ex   (jc) = albsi_new(ic)              ! sea-ice albedo at new time level [-]

        t_g_ex(jc,isub_seaice)  = tice_new(ic)
        ! surface saturation specific humidity (uses saturation water vapor pressure
        ! over ice)
        qv_s_ex(jc,isub_seaice) = spec_humi( sat_pres_ice(tice_new(ic)), ps_ex(jc) )
      ENDDO  ! ic


      ! Update dynamic sea-ice index list
      !
      CALL update_idx_lists_sea (                                    &
        &   hice_n        = h_ice_ex(:),                             &!in
        &   pres_sfc      = ps_ex(:),                                &!in
        &   idx_lst_spw   = ext_data%atm%idx_lst_spw(:,jb),          &!inout
        &   spw_count     = ext_data%atm%spw_count(jb),              &!inout
        &   idx_lst_spi   = ext_data%atm%idx_lst_spi(:,jb),          &!inout
        &   spi_count     = ext_data%atm%spi_count(jb),              &!inout
        &   frac_t_ice    = ext_data%atm%frac_t(:,jb,isub_seaice),   &!inout
        &   frac_t_water  = ext_data%atm%frac_t(:,jb,isub_water),    &!inout
        &   fr_seaice     = fr_seaice(:),                            &!inout
        &   hice_old      = h_ice_ex(:),                             &!inout
        &   tice_old      = t_ice_ex(:),                             &!inout
        &   albsi_now     = alb_si_ex(:),                            &!inout
        &   albsi_new     = alb_si_ex(:),                            &!inout  !DR quick hack (2nd time level miss)
        &   t_g_t_now     = t_g_ex(:,isub_water),                    &!inout  !DR quick hack
        &   t_g_t_new     = t_g_ex(:,isub_water),                    &!inout
        &   t_s_t_now     = t_s_ex(:,isub_water),                    &!inout  !DR quick hack
        &   t_s_t_new     = t_s_ex(:,isub_water),                    &!inout
        &   qv_s_t        = qv_s_ex(:,isub_water),                   &!inout
        &   t_seasfc      = t_seasfc_dummy(:)                        )!inout  ! DUMMY: not used yet
!!DR Currently, t_seasfc is not re-initialized for newly generated sea points when the EDMF scheme is used.
!!DR This is simply because p_lnd_diag%t_seasfc is not readily vailable. This is not critical, however it 
!!DR should however be fixed at some time.

  END SUBROUTINE nwp_seaice


END MODULE mo_nwp_sfc_interface_edmf

