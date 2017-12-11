!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_psrad_interface
  USE mo_kind,                       ONLY: wp
  USE mo_physical_constants,         ONLY: avo, amd, amw, amco2, amch4, amn2o, &
                                           amo3, amo2, amc11, amc12, zemiss_def
  USE mo_model_domain,               ONLY: t_patch
  USE mo_parallel_config,            ONLY: nproma
  USE mo_impl_constants,             ONLY: min_rlcell_int, grf_bdywidth_c
  USE mo_loopindices,                ONLY: get_indices_c
  USE mo_psrad_radiation_parameters, ONLY: rad_perm
  USE mo_psrad_solar_parameters,     ONLY: psctm, ssi_factor
  USE mo_psrad_general,              ONLY: ncfc, ngas, nbndsw, nbndlw, nmixture
  USE mo_psrad_cloud_optics,         ONLY: cloud_optics
  USE mo_bc_aeropt_kinne,            ONLY: set_bc_aeropt_kinne  
  USE mo_bc_aeropt_stenchikov,       ONLY: add_bc_aeropt_stenchikov 
  USE mo_bc_aeropt_splumes,          ONLY: add_bc_aeropt_splumes
  USE mo_psrad_gas_optics,           ONLY: precomputation
  USE mo_psrad_lrtm_driver,          ONLY: lrtm
  USE mo_psrad_srtm_driver,          ONLY: srtm, srtm_diags
  USE mo_psrad_random,               ONLY: rng_seed_size => seed_size
  USE mo_rad_diag,                   ONLY: rad_aero_diag
  USE mtime,                         ONLY: datetime
  USE mo_timer,                      ONLY: ltimer, timer_start, timer_stop, &
   &                                       timer_lrtm, timer_srtm
#ifdef PSRAD_TIMING
  USE mo_timer,                      ONLY: timer_rrtm_coeffs, timer_cloud_optics, &
                                           timer_psrad_scaling, timer_psrad_aerosol
#endif

  IMPLICIT NONE

  PRIVATE

  REAL(wp), PARAMETER :: pressure_scale = 100._wp, &
                         inverse_pressure_scale = 0.01_wp

  PUBLIC :: psrad_interface, pressure_scale, inverse_pressure_scale
  
CONTAINS
  !>
  !! @brief arranges input and calls rrtm sw and lw routines
  !! 
  !! @par Revision History
  !! Original Source Rewritten and renamed by B. Stevens (2009-08)
  !!
  !! @remarks
  !!   Because the RRTM indexes vertical levels differently than ECHAM a chief
  !!   function of thise routine is to reorder the input in the vertical.  In 
  !!   addition some cloud physical properties are prescribed, which are 
  !!   required to derive cloud optical properties
  !!
  !! @par The gases are passed into RRTM via two multi-constituent arrays: 
  !!   zwkl and wx_r. zwkl has ngas and  wx_r has ncfc species  
  !!   The species are identifed as follows.  
  !!     ZWKL [#/cm2]          WX_R [#/cm2]
  !!    index = 1 => H20     index = 1 => n/a
  !!    index = 2 => CO2     index = 2 => CFC11
  !!    index = 3 =>  O3     index = 3 => CFC12
  !!    index = 4 => N2O     index = 4 => n/a
  !!    index = 5 => n/a
  !!    index = 6 => CH4
  !!    index = 7 => O2
  ! Notice that the above are not identical to the current ones in 
  ! psrad_general
  !-------------------------------------------------------------------
  SUBROUTINE psrad_interface(                                               &
      & patch,                                                              &
      & irad_aero       ,klev            ,ktype                            ,& 
      & loland          ,loglac          ,this_datetime                    ,&
      & pcos_mu0        ,daylght_frc                                       ,&
      & alb_vis_dir     ,alb_nir_dir     ,alb_vis_dif     ,alb_nir_dif     ,&
      & zf              ,zh              ,dz                               ,&
      & pp_sfc          ,pp_fl                                             ,&
      & tk_sfc          ,tk_fl           ,tk_hl                            ,&
      & xm_dry          ,xm_vap          ,xm_liq          ,xm_ice          ,&
      & cdnc            ,xc_frc                                            ,&
      & xm_co2          ,xm_ch4          ,xm_n2o          ,xm_cfc          ,&
      & xm_o3           ,xm_o2                                             ,&
      & lw_upw          ,lw_upw_clr      ,lw_dnw          ,lw_dnw_clr      ,&
      & sw_upw          ,sw_upw_clr      ,sw_dnw          ,sw_dnw_clr      ,&
      & vis_dn_dir_sfc  ,par_dn_dir_sfc  ,nir_dn_dir_sfc                   ,&
      & vis_dn_dff_sfc  ,par_dn_dff_sfc  ,nir_dn_dff_sfc                   ,&
      & vis_up_sfc      ,par_up_sfc      ,nir_up_sfc                       )     
     !-------------------------------------------------------------------

    TYPE(t_patch)   ,TARGET ,INTENT(in)    :: patch

    INTEGER,INTENT(IN)  :: &
         irad_aero,        & !< aerosol control
         klev,             & !< number of levels
!!$         ktrac,         & !< number of tracers
         ktype(:,:)          !< type of convection

    LOGICAL,INTENT(IN) ::              &
         loland(:,:),                & !< land sea mask, land=.true.
         loglac(:,:)                   !< glacier mask, glacier=.true.

    TYPE(datetime), POINTER ::  this_datetime !< actual time step

    REAL(WP),INTENT(IN)  :: &
         pcos_mu0(:,:),     & !< mu0 for solar zenith angle
         daylght_frc(:,:),  & !< daylight fraction; with diurnal cycle 0 or 1, with zonal mean in [0,1]
         alb_vis_dir(:,:),  & !< surface albedo for vis range and dir light
         alb_nir_dir(:,:),  & !< surface albedo for NIR range and dir light
         alb_vis_dif(:,:),  & !< surface albedo for vis range and dif light
         alb_nir_dif(:,:),  & !< surface albedo for NIR range and dif light
         zf(:,:,:),         & !< geometric height at full level in m
         zh(:,:,:),         & !< geometric height at half level in m
         dz(:,:,:),         & !< geometric height thickness in m
         pp_sfc(:,:),       & !< surface pressure in Pa
         pp_fl(:,:,:),      & !< full level pressure in Pa
         tk_sfc(:,:),       & !< surface temperature in K
         tk_fl(:,:,:),      & !< full level temperature in K
         tk_hl(:,:,:),      & !< half level temperature in K
         xm_dry(:,:,:),     & !< dry air     mass in kg/m2
         xm_vap(:,:,:),     & !< water vapor mass in kg/m2
         xm_liq(:,:,:),     & !< cloud water mass in kg/m2
         xm_ice(:,:,:),     & !< cloud ice   mass in kg/m2
         cdnc(:,:,:),       & !< cloud nuclei concentration
         xc_frc(:,:,:),     & !< fractional cloud cover
         xm_co2(:,:,:),     & !< co2 mass in kg/m2
         xm_ch4(:,:,:),     & !< ch4 mass in kg/m2
         xm_n2o(:,:,:),     & !< n2o mass in kg/m2
         xm_cfc(:,:,:,:),   & !< cfc mass in kg/m2
         xm_o3(:,:,:),      & !< o3  mass in kg/m2
         xm_o2(:,:,:)       !< o2  mass in kg/m2
!!$         xm_trc(kbdim,klev,ktrac)        !< tracer mass mixing ratios

    REAL(wp), INTENT(OUT)   :: &
      & lw_dnw_clr(:,:,:),& !< Clear-sky downward longwave  at all levels
      & lw_upw_clr(:,:,:),& !< Clear-sky upward   longwave  at all levels
      & sw_dnw_clr(:,:,:),& !< Clear-sky downward shortwave at all levels
      & sw_upw_clr(:,:,:),& !< Clear-sky upward   shortwave at all levels
      & lw_dnw(:,:,:),    & !< All-sky   downward longwave  at all levels
      & lw_upw(:,:,:),    & !< All-sky   upward   longwave  at all levels
      & sw_dnw(:,:,:),    & !< All-sky   downward shortwave at all levels
      & sw_upw(:,:,:)       !< All-sky   upward   shortwave at all levels

    REAL(wp), INTENT(OUT) :: &
         vis_dn_dir_sfc(:,:), & !< Diffuse downward flux surface visible radiation 
         par_dn_dir_sfc(:,:), & !< Diffuse downward flux surface PAR
         nir_dn_dir_sfc(:,:), & !< Diffuse downward flux surface near-infrared radiation
         vis_dn_dff_sfc(:,:), & !< Direct  downward flux surface visible radiation 
         par_dn_dff_sfc(:,:), & !< Direct  downward flux surface PAR
         nir_dn_dff_sfc(:,:), & !< Direct  downward flux surface near-infrared radiation
         vis_up_sfc    (:,:), & !< Upward  flux surface visible radiation 
         par_up_sfc    (:,:), & !< Upward  flux surface PAR
         nir_up_sfc    (:,:)    !< Upward  flux surface near-infrared radiation
 
    INTEGER  :: jg             
    INTEGER  :: i_nchdom, rl_start, rl_end
    INTEGER  :: i_startblk,i_endblk
    INTEGER  :: jb             !< block index
    INTEGER  :: jcs, jce       !< start/end column index within this block

    jg         = patch%id
    rl_start   = grf_bdywidth_c+1
    rl_end     = min_rlcell_int
 
    i_nchdom   = MAX(1,patch%n_childdom)
    i_startblk = patch%cells%start_blk(rl_start,1)
    i_endblk   = patch%cells%end_blk(rl_end,i_nchdom)

!$OMP PARALLEL DO PRIVATE(jcs,jce)
    DO jb = i_startblk,i_endblk
       
      CALL get_indices_c(patch, jb,i_startblk,i_endblk, jcs,jce, rl_start, rl_end)

    !-------------------------------------------------------------------
      CALL psrad_interface_onBlock(            jg              ,jb                   ,&
        & irad_aero       ,jce             ,nproma          ,klev                    ,& 
        & ktype(:,jb)                                                                ,&
        & loland(:,jb)    ,loglac(:,jb)    ,this_datetime                            ,&
        & pcos_mu0(:,jb)  ,daylght_frc(:,jb)                                         ,&
        & alb_vis_dir(:,jb) ,alb_nir_dir(:,jb),alb_vis_dif(:,jb)                     ,&
        & alb_nir_dif(:,jb)                                                          ,&
        & zf(:,:,jb)      ,zh(:,:,jb)      ,dz(:,:,jb)                               ,&
        & pp_sfc(:,jb)    ,pp_fl(:,:,jb)                                             ,&
        & tk_sfc(:,jb)    ,tk_fl(:,:,jb)   ,tk_hl(:,:,jb)                            ,&
        & xm_dry(:,:,jb)  ,xm_vap(:,:,jb)  ,xm_liq(:,:,jb), xm_ice(:,:,jb)           ,&
        & cdnc(:,:,jb)    ,xc_frc(:,:,jb)                                            ,&
        & xm_co2(:,:,jb)  ,xm_ch4(:,:,jb)  ,xm_n2o (:,:,jb) ,xm_cfc (:,:,:,jb)       ,&
        & xm_o3(:,:,jb)   ,xm_o2(:,:,jb)                                             ,&
        & lw_upw (:,:,jb) ,lw_upw_clr (:,:,jb) ,lw_dnw(:,:,jb), lw_dnw_clr (:,:,jb)  ,&
        & sw_upw(:,:,jb)  ,sw_upw_clr(:,:,jb)  ,sw_dnw(:,:,jb) ,sw_dnw_clr(:,:,jb)   ,&
        & vis_dn_dir_sfc(:,jb) ,par_dn_dir_sfc(:,jb) ,nir_dn_dir_sfc(:,jb)           ,&
        & vis_dn_dff_sfc(:,jb) ,par_dn_dff_sfc(:,jb) ,nir_dn_dff_sfc(:,jb)           ,&
        & vis_up_sfc(:,jb)      ,par_up_sfc(:,jb)      ,nir_up_sfc(:,jb)                       )
   END DO
!$OMP END PARALLEL DO  
    !-------------------------------------------------------------------  
  END SUBROUTINE psrad_interface
 ! -------------------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  !>
  !! @brief arranges input and calls rrtm sw and lw routines
  !! 
  !! @par Revision History
  !! Original Source Rewritten and renamed by B. Stevens (2009-08)
  !!
  !! @remarks
  !!   Because the RRTM indexes vertical levels differently than ECHAM a chief
  !!   function of thise routine is to reorder the input in the vertical.  In 
  !!   addition some cloud physical properties are prescribed, which are 
  !!   required to derive cloud optical properties
  !!
  !! @par The gases are passed into RRTM via two multi-constituent arrays: 
  !!   zwkl and wx_r. zwkl has ngas and  wx_r has ncfc species  
  !!   The species are identifed as follows.  
  !!     ZWKL [#/cm2]          WX_R [#/cm2]
  !!    index = 1 => H20     index = 1 => n/a
  !!    index = 2 => CO2     index = 2 => CFC11
  !!    index = 3 =>  O3     index = 3 => CFC12
  !!    index = 4 => N2O     index = 4 => n/a
  !!    index = 5 => n/a
  !!    index = 6 => CH4
  !!    index = 7 => O2
  ! Notice that The above are not identical to the current ones in 
  ! psrad_general

  SUBROUTINE psrad_interface_onBlock(jg   ,krow                             ,&
       & iaero           ,kproma          ,kbdim           ,klev            ,&
!!$       & ktrac                                                              ,&
       & ktype                                                              ,&
       & laland          ,laglac          ,this_datetime                    ,&
       & pcos_mu0        ,daylght_frc                                       ,&
       & alb_vis_dir     ,alb_nir_dir     ,alb_vis_dif     ,alb_nir_dif     ,&
       & zf              ,zh              ,dz                               ,&
       & pp_sfc          ,pp_fl                                             ,&
       & tk_sfc          ,tk_fl           ,tk_hl                            ,&
       & xm_dry          ,xm_vap          ,xm_liq          ,xm_ice          ,&
       & cdnc            ,cld_frc                                           ,&
       & xm_co2          ,xm_ch4          ,xm_n2o          ,xm_cfc          ,&
       & xm_o3           ,xm_o2                                             ,&
!!$       & xm_trc                                                             ,&
       & flx_uplw        ,flx_uplw_clr    ,flx_dnlw        ,flx_dnlw_clr    ,&
       & flx_upsw        ,flx_upsw_clr    ,flx_dnsw        ,flx_dnsw_clr    ,&
       & vis_dn_dir_sfc  ,par_dn_dir_sfc  ,nir_dn_dir_sfc                   ,&
       & vis_dn_dff_sfc  ,par_dn_dff_sfc  ,nir_dn_dff_sfc                   ,&
       & vis_up_sfc      ,par_up_sfc      ,nir_up_sfc                       )

    INTEGER,INTENT(IN)  :: &
         jg,               & !< domain index
         krow,             & !< first dimension of 2-d arrays
         iaero,            & !< aerosol control
         kproma,           & !< number of longitudes
         kbdim,            & !< first dimension of 2-d arrays
         klev,             & !< number of levels
!!$         ktrac,         & !< number of tracers
         ktype(:)        !< type of convection

    LOGICAL,INTENT(IN) :: &
         laland(:),   & !< land sea mask, land=.true.
         laglac(:)      !< glacier mask, glacier=.true.

    TYPE(datetime), POINTER ::  this_datetime !< actual time step

    REAL(WP),INTENT(IN)  ::    &
         pcos_mu0(:),      & !< mu0 for solar zenith angle
         daylght_frc(:),   & !< daylight fraction; with diurnal cycle 0 or 1, with zonal mean in [0,1]
         alb_vis_dir(:),   & !< surface albedo for vis range and dir light
         alb_nir_dir(:),   & !< surface albedo for NIR range and dir light
         alb_vis_dif(:),   & !< surface albedo for vis range and dif light
         alb_nir_dif(:),   & !< surface albedo for NIR range and dif light
         zf(:,:),          & !< geometric height at full level in m
         zh(:,:),          & !< geometric height at half level in m
         dz(:,:),          & !< geometric height thickness in m
         pp_sfc(:),        & !< surface pressure in Pa
         pp_fl(:,:),       & !< full level pressure in Pa
         tk_sfc(:),        & !< surface temperature in K
         tk_fl(:,:),       & !< full level temperature in K
         tk_hl(:,:),    & !< half level temperature in K
         xm_dry(:,:),   & !< dry air     mass in kg/m2
         xm_vap(:,:),   & !< water vapor mass in kg/m2
         xm_liq(:,:),   & !< cloud water mass in kg/m2
         xm_ice(:,:),   & !< cloud ice   mass in kg/m2
         cdnc(:,:),     & !< cloud nuclei concentration
         cld_frc(:,:),  & !< fractional cloud cover
         xm_co2(:,:),   & !< co2 mass in kg/m2
         xm_ch4(:,:),   & !< ch4 mass in kg/m2
         xm_n2o(:,:),   & !< n2o mass in kg/m2
         xm_cfc(kbdim,klev,2), & !< cfc mass in kg/m2
         xm_o3(:,:),    & !< o3  mass in kg/m2
         xm_o2(:,:)       !< o2  mass in kg/m2
!!$         xm_trc(kbdim,klev,ktrac)        !< tracer mass mixing ratios

    REAL (wp), INTENT (OUT) ::       &
         flx_uplw    (:,:), & !<   upward LW flux profile, all sky
         flx_uplw_clr(:,:), & !<   upward LW flux profile, clear sky
         flx_dnlw    (:,:), & !< downward LW flux profile, all sky
         flx_dnlw_clr(:,:), & !< downward LW flux profile, clear sky
         flx_upsw    (:,:), & !<   upward SW flux profile, all sky
         flx_upsw_clr(:,:), & !<   upward SW flux profile, clear sky
         flx_dnsw    (:,:), & !< downward SW flux profile, all sky
         flx_dnsw_clr(:,:)    !< downward SW flux profile, clear sky

    REAL (wp), INTENT (OUT) :: &
         vis_dn_dir_sfc(:) , & !< Diffuse downward flux surface visible radiation 
         par_dn_dir_sfc(:) , & !< Diffuse downward flux surface PAR
         nir_dn_dir_sfc(:) , & !< Diffuse downward flux surface near-infrared radiation
         vis_dn_dff_sfc(:) , & !< Direct  downward flux surface visible radiation 
         par_dn_dff_sfc(:) , & !< Direct  downward flux surface PAR
         nir_dn_dff_sfc(:) , & !< Direct  downward flux surface near-infrared radiation
         vis_up_sfc    (:) , & !< Upward  flux surface visible radiation 
         par_up_sfc    (:) , & !< Upward  flux surface PAR
         nir_up_sfc    (:)     !< Upward  flux surface near-infrared radiation

    ! -----------------------------------------------------------------------

    INTEGER  :: jk, jl !< loop indicies

    REAL(wp) ::                      &
         zsemiss     (kbdim,nbndlw), & !< LW surface emissivity by band
         bnd_wght(nbndlw),           & 
         per_band_flux(kbdim,nbndsw,3)
    ! --- local scaled variables
    INTEGER :: icldlyr_loc(kbdim,klev) !< index for clear or cloudy
    REAL(wp) ::                   &
         col_dry_loc(kbdim,klev), & !< number of molecules/cm2 of
         cdnc_loc   (kbdim,klev), & !< cloud nuclei concentration
         cld_frc_loc(kbdim,klev), & !< secure cloud fraction
         ziwp_loc   (kbdim,klev), & !< in cloud ice water content       [g/m2]
         ziwc_loc   (kbdim,klev), & !< in cloud ice water concentration [g/m3]
         zlwp_loc   (kbdim,klev), & !< in cloud liquid water content    [g/m2]
         zlwc_loc   (kbdim,klev), & !< in cloud liquid water concentration  [g/m3]
         wkl_loc       (kbdim,klev,ngas),& !< number of molecules/cm2 of
         wx_loc        (kbdim,klev,ncfc),& !< number of molecules/cm2 of
         cld_tau_lw_loc(kbdim,klev,nbndlw), & !< LW optical thickness of clouds
         cld_tau_sw_loc(kbdim,klev,nbndsw), & !< extincion
         cld_cg_sw_loc (kbdim,klev,nbndsw), & !< asymmetry factor
         cld_piz_sw_loc(kbdim,klev,nbndsw), & !< single scattering albedo
         aer_tau_lw_loc(kbdim,klev,nbndlw), & !< LW optical thickness of aerosols
         aer_tau_sw_loc(kbdim,klev,nbndsw), & !< aerosol optical thickness
         aer_cg_sw_loc (kbdim,klev,nbndsw), & !< aerosol asymmetry factor
         aer_piz_sw_loc(kbdim,klev,nbndsw) !< aerosol single scattering albedo

    ! --- vertically reversed _loc variables
    REAL(wp) ::                            &
         pm_fl_vr  (kbdim,klev),           & !< full level pressure [mb] 
         cld_frc_vr(kbdim,klev),           & !< secure cloud fraction
         aer_tau_lw_vr(kbdim,klev,nbndlw), & !< LW optical thickness of aerosols
         aer_tau_sw_vr(kbdim,klev,nbndsw), & !< aerosol optical thickness
         aer_cg_sw_vr (kbdim,klev,nbndsw), & !< aerosol asymmetry factor
         aer_piz_sw_vr(kbdim,klev,nbndsw) !< aerosol single scattering albedo

    REAL(wp) ::                &
         re_drop (kbdim,klev), & !< effective radius of liquid
         re_cryst(kbdim,klev), & !< effective radius of ice
         x_cdnc  (kbdim)         !< Scale factor for Cloud Droplet Number Concentration

    !
    ! Random seeds for sampling. Needs to get somewhere upstream 
    !
    INTEGER :: rnseeds(kbdim,rng_seed_size)

    REAL(wp), TARGET  :: actual_scaleminorn2(kbdim,klev), &
                         actual_gases(kbdim,klev,ngas), &
                         actual_ratio(kbdim,2,klev,nmixture) 

    REAL(wp), POINTER :: gases(:,:,:), &
                         ratio(:,:,:,:), &
                         scaleminorn2(:,:)

    REAL(wp)          :: fac(kbdim,2,2,klev)

    INTEGER           :: laytrop(kbdim), & !< tropopause layer index
                         iabs(kbdim,2,2,klev), &
                         jp(kbdim,klev), &
                         indminor(kbdim,klev), &
                         h2o_index(kbdim,klev,2)

    REAL(wp)          :: h2o_factor(kbdim,klev,2), &
                         h2o_fraction(kbdim,klev,2), &
                         colbrd(kbdim,klev), &
                         colmol(kbdim,klev), &
                         minorfrac(kbdim,klev), &
                         scaleminor(kbdim,klev)

    gases => actual_gases
    ratio => actual_ratio
    scaleminorn2 => actual_scaleminorn2
!DIR$ NOFMA

    ! 1.0 Constituent properties 
    !--------------------------------
!IBM* ASSERT(NODEPS)
#ifdef PSRAD_TIMING
    IF (ltimer) CALL timer_start(timer_psrad_scaling)
#endif
    DO jk = 1, klev
      DO jl = 1, kproma
        !
        ! --- Cloud liquid and ice mass: [kg/m2 in cell] --> [g/m2 in cloud]
        !
        cld_frc_loc(jl,jk) = MAX(EPSILON(1.0_wp),cld_frc(jl,jk))
        ziwp_loc(jl,jk)    = xm_ice(jl,jk)*1000.0_wp/cld_frc_loc(jl,jk)
        zlwp_loc(jl,jk)    = xm_liq(jl,jk)*1000.0_wp/cld_frc_loc(jl,jk)
      END DO
    END DO
    !
    ! --- control for zero, infintesimal or negative cloud fractions
    !
    WHERE (cld_frc_loc(1:kproma,:) > 2.0_wp*EPSILON(1.0_wp))
      icldlyr_loc(1:kproma,:) = 1
    ELSEWHERE
      icldlyr_loc(1:kproma,:) = 0
      ziwp_loc(1:kproma,:) = 0.0_wp
      zlwp_loc(1:kproma,:) = 0.0_wp
    END WHERE
    !
    ! --- main constituent vertical reordering and unit conversion
    !
!IBM* ASSERT(NODEPS)
    DO jk = 1, klev
      DO jl = 1, kproma
        !
        ! --- cloud water and ice concentrations [kg/m3]
        !
        ziwc_loc(jl,jk)    = ziwp_loc(jl,jk)/dz(jl,jk)
        zlwc_loc(jl,jk)    = zlwp_loc(jl,jk)/dz(jl,jk)
        !
        ! --- cloud droplet number concentration  [1/m3] --> [1/cm3] ?
        !
        cdnc_loc(jl,jk)    = cdnc(jl,jk)*1.e-6_wp
        !
        ! --- dry air: [kg/m2] --> [molecules/cm2]
        !
        col_dry_loc(jl,jk) = 0.1_wp * avo * xm_dry(jl,jk)/amd
        !
        ! --- H2O, CO2, O3, N2O, CH4, O2: [kg/m2] --> [molecules/cm2]
        !
        wkl_loc(jl,jk,:)   = 0.0_wp
        wkl_loc(jl,jk,1)   = 0.1_wp * avo * xm_vap(jl,jk)/amw
        wkl_loc(jl,jk,2)   = 0.1_wp * avo * xm_co2(jl,jk)/amco2
        wkl_loc(jl,jk,3)   = 0.1_wp * avo * xm_ch4(jl,jk)/amch4
        wkl_loc(jl,jk,4)   = 0.1_wp * avo * xm_o2 (jl,jk)/amo2
        wkl_loc(jl,jk,5)   = 0.1_wp * avo * xm_o3 (jl,jk)/amo3
        wkl_loc(jl,jk,6)   = 0.1_wp * avo * xm_n2o(jl,jk)/amn2o
        !
        ! --- CFC11, CFC12: [kg/m2] --> [molecules/cm2]
        !
        wx_loc(jl,jk,:)    = 0.0_wp
        wx_loc(jl,jk,2)    = 0.1_wp * avo * xm_cfc(jl,jk,1)/amc11
        wx_loc(jl,jk,3)    = 0.1_wp * avo * xm_cfc(jl,jk,2)/amc12
        !
      END DO
    END DO
    CALL flip_ud(kproma, pp_fl,       pm_fl_vr)
    CALL flip_ud(kproma, cld_frc_loc, cld_frc_vr)

    ! 2.0 Surface Properties
    ! --------------------------------
    zsemiss(1:kproma,:) = zemiss_def 
    !
    ! 3.0 Particulate Optical Properties
    ! --------------------------------
#ifdef PSRAD_TIMING
    IF (ltimer) CALL timer_stop(timer_psrad_scaling)
    IF (ltimer) CALL timer_start(timer_psrad_aerosol)
#endif

! IF (aero == ...) THEN
! iaero=0: No aerosol
! iaero=13: only tropospheric Kinne aerosols
! iaero=14: only Stenchikov's volcanic aerosols
! iaero=15: tropospheric Kinne aerosols + volcanic Stenchikov's aerosols
! set all aerosols to zero first
    aer_tau_lw_vr(:,:,:) = 0.0_wp
    aer_tau_sw_vr(:,:,:) = 0.0_wp
    aer_piz_sw_vr(:,:,:) = 1.0_wp
    aer_cg_sw_vr(:,:,:)  = 0.0_wp
    IF (iaero==13 .OR. iaero==15 .OR. iaero==18) THEN
! iaero=13: only Kinne aerosols are used
! iaero=15: Kinne aerosols plus Stenchikov's volcanic aerosols are used
! iaero=18: Kinne background aerosols (of natural origin, 1850) are set
      CALL set_bc_aeropt_kinne( this_datetime                          ,&
           & kproma           ,kbdim                 ,klev             ,&
           & krow             ,nbndsw                ,nbndlw           ,&
           & zf               ,dz                                      ,&
           & aer_tau_sw_vr    ,aer_piz_sw_vr         ,aer_cg_sw_vr     ,&
           & aer_tau_lw_vr                                              )
    END IF
    IF (iaero==14 .OR. iaero==15 .OR. iaero==18) THEN
! iaero=14: only Stechnikov's volcanic aerosols are used (added to zero)
! iaero=15: Stenchikov's volcanic aerosols are added to Kinne aerosols
! iaero=18: Stenchikov's volcanic aerosols are added to Kinne background
!           aerosols (of natural origin, 1850) 
      CALL add_bc_aeropt_stenchikov( this_datetime   ,jg               ,&
           & kproma           ,kbdim                 ,klev             ,&
           & krow             ,nbndsw                ,nbndlw           ,&
           & dz               ,pp_fl                                   ,&
           & aer_tau_sw_vr    ,aer_piz_sw_vr         ,aer_cg_sw_vr     ,&
           & aer_tau_lw_vr                                              )
    END IF
!!$    IF (iaero==16) THEN
!!$      CALL add_aop_volc_ham( &
!!$           & kproma           ,kbdim                 ,klev             ,&
!!$           & krow             ,nbndlw                ,nbndsw           ,&
!!$           & aer_tau_lw_vr    ,aer_tau_sw_vr         ,aer_piz_sw_vr    ,&
!!$           & aer_cg_sw_vr                                               )
!!$    END IF
!!$    IF (iaero==17) THEN
!!$      CALL add_aop_volc_crow( &
!!$           & kproma           ,kbdim                 ,klev             ,&
!!$           & krow             ,nbndlw                ,nbndsw           ,&
!!$           & aer_tau_lw_vr    ,aer_tau_sw_vr         ,aer_piz_sw_vr    ,&
!!$           & aer_cg_sw_vr                                               )
!!$    END IF
    IF (iaero==18) THEN
! iaero=18: Simple plumes are added to Stenchikov's volcanic aerosols 
!           and Kinne background aerosols (of natural origin, 1850) 
      CALL add_bc_aeropt_splumes(jg                                    ,&
           & kproma           ,kbdim                 ,klev             ,&
           & krow             ,nbndsw                ,this_datetime    ,&
           & zf               ,dz                    ,zh(:,klev+1)     ,&
           & aer_tau_sw_vr    ,aer_piz_sw_vr         ,aer_cg_sw_vr     ,&
           & x_cdnc                                                     )
    END IF

    CALL rad_aero_diag (                                  &
      & jg              ,krow            ,kproma          , &
      & kbdim           ,klev            ,nbndlw          , &
      & nbndsw          ,aer_tau_lw_vr   ,aer_tau_sw_vr   , &
      & aer_piz_sw_vr   ,aer_cg_sw_vr                       )

    DO jl = 1,nbndlw
      CALL flip_ud(kproma, aer_tau_lw_vr(:,:,jl), aer_tau_lw_loc(:,:,jl))
    ENDDO
    DO jl = 1,nbndsw
      CALL flip_ud(kproma, aer_tau_sw_vr(:,:,jl), aer_tau_sw_loc(:,:,jl))
      CALL flip_ud(kproma, aer_cg_sw_vr(:,:,jl),  aer_cg_sw_loc(:,:,jl))
      CALL flip_ud(kproma, aer_piz_sw_vr(:,:,jl), aer_piz_sw_loc(:,:,jl))
    ENDDO

#ifdef PSRAD_TIMING
    IF (ltimer) CALL timer_stop(timer_psrad_aerosol)
    IF (ltimer) CALL timer_start(timer_cloud_optics)
#endif
    CALL cloud_optics(                                                      &
         & laglac         ,laland         ,kproma         ,kbdim           ,& 
         & klev           ,ktype                                           ,&       
         & icldlyr_loc    ,zlwp_loc       ,ziwp_loc       ,zlwc_loc        ,&
         & ziwc_loc       ,cdnc_loc       ,cld_tau_lw_loc ,cld_tau_sw_loc  ,&
         & cld_piz_sw_loc ,cld_cg_sw_loc  ,re_drop        ,re_cryst         )  
#ifdef PSRAD_TIMING
    IF (ltimer) CALL timer_stop(timer_cloud_optics)
#endif

    !
    ! 4.0 Radiative Transfer Routines
    ! --------------------------------
    !
    ! Seeds for random numbers come from least significant digits of 
    ! pressure field 
    !
    rnseeds(1:kproma,1:rng_seed_size) = &
      int((inverse_pressure_scale * pm_fl_vr(1:kproma,1:rng_seed_size) -  &
      int(inverse_pressure_scale * pm_fl_vr(1:kproma,1:rng_seed_size)))* 1E9 + rad_perm)
    ! Calculate information needed by the radiative transfer routine
    ! that is specific to this atmosphere, especially some of the 
    ! coefficients and indices needed to compute the optical depths
    ! by interpolating data from stored reference atmospheres. 
    ! The coefficients are functions of temperature and pressure and 
    ! remain the same for all g-point samples.
    ! If gas concentrations, temperatures, or pressures vary with sample (ig) 
    ! the coefficients need to be calculated inside the loop over samples

#ifdef PSRAD_TIMING
    IF (ltimer) CALL timer_start(timer_rrtm_coeffs)
#endif
    CALL precomputation(kproma        ,kbdim         ,klev        , &
         & .false.     ,&
         & pp_fl   ,tk_fl         ,col_dry_loc  , &
         & wkl_loc     ,laytrop       ,jp            ,iabs        , &
         & gases       ,colbrd        ,colmol        ,fac         , &
         & ratio       ,h2o_factor    ,h2o_fraction  ,h2o_index   , &
         & minorfrac   ,scaleminor    ,scaleminorn2  ,indminor)
#ifdef PSRAD_TIMING
    IF (ltimer) CALL timer_stop(timer_rrtm_coeffs)
#endif
    IF (ltimer) CALL timer_start(timer_lrtm)
    CALL lrtm(kproma,      kbdim,          klev,                          &
         & pp_fl,          pp_sfc,         tk_fl,          tk_hl,         &
         & tk_sfc,         wkl_loc,        wx_loc,         col_dry_loc,   &
         & zsemiss,        cld_frc_loc,    cld_tau_lw_loc, aer_tau_lw_loc,&
         & rnseeds,        gases,          ratio,          scaleminorn2,  &
         & fac,            laytrop,        iabs,           jp,            &
         & indminor,       h2o_factor,     h2o_fraction,   h2o_index,     &
         & colbrd,         minorfrac,      scaleminor,     &
         & flx_uplw,       flx_dnlw,       flx_uplw_clr, flx_dnlw_clr )
    IF (ltimer) CALL timer_stop(timer_lrtm)

    !
    ! Reset random seeds so SW doesn't depend on what's happened in LW but is also independent
    !
    rnseeds(1:kproma,1:rng_seed_size) = &
      int((inverse_pressure_scale * pm_fl_vr(1:kproma,rng_seed_size:1:-1) - &
      int(inverse_pressure_scale * pm_fl_vr(1:kproma,rng_seed_size:1:-1)))* 1E9 + rad_perm)

    ! Potential pitfall - we're passing every argument but some may not be present
    IF (ltimer) CALL timer_start(timer_srtm)
    CALL srtm(kproma,       kbdim,          klev,                          &
         &  alb_vis_dir,    alb_vis_dif,    alb_nir_dir,   alb_nir_dif,    &
         &  pcos_mu0,       daylght_frc,    ssi_factor,    psctm,          &
         &  cld_frc_loc,    cld_tau_sw_loc, cld_cg_sw_loc,                 &
         &  cld_piz_sw_loc, aer_tau_sw_loc, aer_cg_sw_loc, aer_piz_sw_loc, & 
         &  rnseeds,        laytrop,        jp,            iabs,           &
         & gases,           colmol,         fac,           h2o_factor,     &
         & h2o_fraction,    h2o_index,                                     &
         &  flx_dnsw,       flx_upsw,       flx_dnsw_clr,  flx_upsw_clr,   &
         &  bnd_wght,       per_band_flux                                  )

    CALL srtm_diags(kproma, kbdim,          per_band_flux,  &
          & vis_dn_dir_sfc, par_dn_dir_sfc, nir_dn_dir_sfc, &
          & vis_dn_dff_sfc, par_dn_dff_sfc, nir_dn_dff_sfc, &
          & vis_up_sfc,     par_up_sfc,     nir_up_sfc)

    IF (ltimer) CALL timer_stop(timer_srtm)

  END SUBROUTINE psrad_interface_onBlock

  SUBROUTINE flip_ud(n, v, u)
    INTEGER,     INTENT(IN) :: n
    REAL(wp),    INTENT(IN) :: v(:,:)
    REAL(wp), INTENT(INOUT) :: u(:,:)
    INTEGER :: m

    m = SIZE(v,2)
    u(1:n,1:m) = v(1:n,m:1:-1)
  END SUBROUTINE flip_ud

END MODULE mo_psrad_interface
