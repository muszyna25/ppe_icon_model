!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_psrad_interface
  USE mo_kind,                       ONLY: wp
  USE mo_physical_constants,         ONLY: avo, amd, amw, amco2, amch4, amn2o, amo3, amo2, amc11, amc12, zemiss_def
  USE mo_exception,                  ONLY: finish, message
  USE mo_model_domain,               ONLY: t_patch
  USE mo_parallel_config            ,ONLY: nproma
  USE mo_impl_constants             ,ONLY: min_rlcell_int, grf_bdywidth_c
  USE mo_loopindices                ,ONLY: get_indices_c
  USE mo_psrad_radiation_parameters, ONLY: rad_perm
  USE mo_psrad_solar_parameters,     ONLY: psctm, ssi_factor
  USE mo_psrad_general,              ONLY: ncfc, ngas, nbndsw, nbndlw
  USE mo_psrad_cloud_optics,         ONLY: cloud_optics
  USE mo_bc_aeropt_kinne,            ONLY: set_bc_aeropt_kinne  
  USE mo_bc_aeropt_stenchikov,       ONLY: add_bc_aeropt_stenchikov 
  USE mo_bc_aeropt_splumes,          ONLY: add_bc_aeropt_splumes
!!$  USE mo_aero_volc,                  ONLY: add_aop_volc
!!$  USE mo_aero_volc_tab,              ONLY: add_aop_volc_ham, add_aop_volc_crow
!!$  USE mo_lrtm_setup,                 ONLY: lrtm_setup
  USE mo_psrad_lrtm_setup,ONLY: setup_lrtm
  USE mo_psrad_srtm_setup,ONLY: setup_srtm ! here, the new name is chosen to distinguish
                                           ! it from lrtm_setup in the "old" version
  USE mo_psrad_lrtm_driver,          ONLY: lrtm
  USE mo_psrad_srtm_driver,          ONLY: srtm
!!$  USE mo_submodel,                   ONLY: lanysubmodel
!!$  USE mo_submodel_interface,         ONLY: radiation_subm_1, radiation_subm_2  
!!$  USE mo_cosp_simulator,             ONLY: cosp_reffl, cosp_reffi,      &
!!$                                           locosp, cosp_f3d, Lisccp_sim,&
!!$                                           cisccp_cldtau3d, cisccp_cldemi3d
  USE mo_psrad_random,               ONLY: seed_size_random
  USE mo_rad_diag,                   ONLY: rad_aero_diag
  USE mtime,                         ONLY: datetime
  USE mo_timer,                      ONLY: ltimer, timer_start, timer_stop, &
    timer_lrtm, timer_srtm

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: setup_psrad, psrad_interface
  
  INTEGER, PARAMETER    :: rng_seed_size = 4
CONTAINS
  !---------------------------------------------------------------------------
  !>
  !! @brief  Sets up (initializes) radation routines
  !! 
  !! @remarks
  !!   Modify preset variables of module MO_RADIATION which control the 
  !!   configuration of the radiation scheme.
  !
  SUBROUTINE setup_psrad
    !
    CALL setup_lrtm
    CALL setup_srtm
    !
    CALL message('','')
    CALL message('','PSrad setup')
    CALL message('','===========')
    CALL message('','- RRTM LW k-distribution data source version = "rrtm_lw_v3.01"')
    CALL message('','- RRTM SW k-distribution data source version = "rrtm_sw_v2.5"')
    CALL message('','- Water vapor continuum version              = "CKD_2.4"')
    CALL message('','')
    !
    IF(seed_size_random() > rng_seed_size) THEN 
      CALL finish('setup_psrad','Random number seed size (rng_seed_size) is smaller than true size')
    END IF
  END SUBROUTINE setup_psrad

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
  ! TODO/BUG?
  ! NOTE: The above are not identical to the current ones in 
  ! psrad_general
  !-------------------------------------------------------------------
  SUBROUTINE psrad_interface(                                               &
      & patch,                                                              &
      & irad_aero     ,klev                                                ,& 
      & ktype                                                              ,&
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

!        & patch,                                                              &
!        & iaero             ,klev                                            ,&
!        & ktype                                                              ,&
!        & laland          ,laglac          ,this_datetime                    ,&
!        & pmu0            ,daylght_frc                                       ,&
!        & alb_vis_dir     ,alb_nir_dir     ,alb_vis_dif     ,alb_nir_dif     ,&
!        & zf              ,zh              ,dz                               ,&
!        & pp_sfc          ,pp_fl                                             ,&
!        & tk_sfc          ,tk_fl           ,tk_hl                            ,&
!        & xm_dry          ,xm_vap          ,xm_liq          ,xm_ice          ,&
!        & cdnc            ,cld_frc                                           ,&
!        & xm_co2          ,xm_ch4          ,xm_n2o          ,xm_cfc          ,&
!        & xm_o3           ,xm_o2                                             ,&
!        & flx_uplw        ,flx_uplw_clr    ,flx_dnlw        ,flx_dnlw_clr    ,&
!        & flx_upsw        ,flx_upsw_clr    ,flx_dnsw        ,flx_dnsw_clr    ,&
!        & vis_dn_dir_sfc  ,par_dn_dir_sfc  ,nir_dn_dir_sfc                   ,&
!        & vis_dn_dff_sfc  ,par_dn_dff_sfc  ,nir_dn_dff_sfc                   ,&
!        & vis_up_sfc      ,par_up_sfc      ,nir_up_sfc                       )

    TYPE(t_patch)   ,TARGET ,INTENT(in)    :: patch

    INTEGER,INTENT(IN)  ::             &
         irad_aero,                    & !< aerosol control
         klev,                         & !< number of levels
!!$         ktrac,                        & !< number of tracers
         ktype(:,:)                    !< type of convection

    LOGICAL,INTENT(IN) ::              &
         loland(:,:),                & !< land sea mask, land=.true.
         loglac(:,:)                   !< glacier mask, glacier=.true.

    TYPE(datetime), POINTER ::  this_datetime !< actual time step

    REAL(WP),INTENT(IN)  ::          &
         pcos_mu0(:,:),              & !< mu0 for solar zenith angle
         daylght_frc(:,:),           & !< daylight fraction; with diurnal cycle 0 or 1, with zonal mean in [0,1]
         alb_vis_dir(:,:),           & !< surface albedo for vis range and dir light
         alb_nir_dir(:,:),           & !< surface albedo for NIR range and dir light
         alb_vis_dif(:,:),           & !< surface albedo for vis range and dif light
         alb_nir_dif(:,:),           & !< surface albedo for NIR range and dif light
         zf(:,:,:),               & !< geometric height at full level in m
         zh(:,:,:),             & !< geometric height at half level in m
         dz(:,:,:),               & !< geometric height thickness in m
         pp_sfc(:,:),                & !< surface pressure in Pa
         pp_fl(:,:,:),            & !< full level pressure in Pa
         tk_sfc(:,:),                & !< surface temperature in K
         tk_fl(:,:,:),            & !< full level temperature in K
         tk_hl(:,:,:),          & !< half level temperature in K
         xm_dry(:,:,:),           & !< dry air     mass in kg/m2
         xm_vap(:,:,:),           & !< water vapor mass in kg/m2
         xm_liq(:,:,:),           & !< cloud water mass in kg/m2
         xm_ice(:,:,:),           & !< cloud ice   mass in kg/m2
         cdnc(:,:,:),             & !< cloud nuclei concentration
         xc_frc(:,:,:),           & !< fractional cloud cover
         xm_co2(:,:,:),           & !< co2 mass in kg/m2
         xm_ch4(:,:,:),           & !< ch4 mass in kg/m2
         xm_n2o(:,:,:),           & !< n2o mass in kg/m2
         xm_cfc(:,:,:,:),         & !< cfc mass in kg/m2
         xm_o3(:,:,:),            & !< o3  mass in kg/m2
         xm_o2(:,:,:)               !< o2  mass in kg/m2
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

    REAL (wp), INTENT (OUT) ::         &
         vis_dn_dir_sfc(:,:)       , & !< Diffuse downward flux surface visible radiation 
         par_dn_dir_sfc(:,:)       , & !< Diffuse downward flux surface PAR
         nir_dn_dir_sfc(:,:)       , & !< Diffuse downward flux surface near-infrared radiation
         vis_dn_dff_sfc(:,:)       , & !< Direct  downward flux surface visible radiation 
         par_dn_dff_sfc(:,:)       , & !< Direct  downward flux surface PAR
         nir_dn_dff_sfc(:,:)       , & !< Direct  downward flux surface near-infrared radiation
         vis_up_sfc    (:,:)       , & !< Upward  flux surface visible radiation 
         par_up_sfc    (:,:)       , & !< Upward  flux surface PAR
         nir_up_sfc    (:,:)           !< Upward  flux surface near-infrared radiation
 
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
        & irad_aero       ,jcs, jce        ,nproma          ,klev                    ,& 
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
  ! TODO/BUG?
  ! NOTE: The above are not identical to the current ones in 
  ! psrad_general

  SUBROUTINE psrad_interface_onBlock(      jg              ,krow            ,&
       & iaero           ,jcs, kproma     ,kbdim           ,klev            ,&
!!$       & ktrac                                                              ,&
       & ktype                                                              ,&
       & laland          ,laglac          ,this_datetime                    ,&
       & pmu0            ,daylght_frc                                       ,&
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

    INTEGER,INTENT(IN)  ::             &
         jg,                           & !< domain index
         krow,                         & !< block index
         iaero,                        & !< aerosol control
         jcs,                          & !< start and
         kproma,                       & !< end index in block
         kbdim,                        & !< first dimension of 2-d arrays
         klev,                         & !< number of levels
!!$         ktrac,                        & !< number of tracers
         ktype(kbdim)                    !< type of convection

    LOGICAL,INTENT(IN) ::              &
         laland(kbdim),                & !< land sea mask, land=.true.
         laglac(kbdim)                   !< glacier mask, glacier=.true.

    TYPE(datetime), POINTER ::  this_datetime !< actual time step

    REAL(WP),INTENT(IN)  ::            &
         pmu0(kbdim),                  & !< mu0 for solar zenith angle
         daylght_frc(kbdim),           & !< daylight fraction; with diurnal cycle 0 or 1, with zonal mean in [0,1]
         alb_vis_dir(kbdim),           & !< surface albedo for vis range and dir light
         alb_nir_dir(kbdim),           & !< surface albedo for NIR range and dir light
         alb_vis_dif(kbdim),           & !< surface albedo for vis range and dif light
         alb_nir_dif(kbdim),           & !< surface albedo for NIR range and dif light
         zf(kbdim,klev),               & !< geometric height at full level in m
         zh(kbdim,klev+1),             & !< geometric height at half level in m
         dz(kbdim,klev),               & !< geometric height thickness in m
         pp_sfc(kbdim),                & !< surface pressure in Pa
         pp_fl(kbdim,klev),            & !< full level pressure in Pa
         tk_sfc(kbdim),                & !< surface temperature in K
         tk_fl(kbdim,klev),            & !< full level temperature in K
         tk_hl(kbdim,klev+1),          & !< half level temperature in K
         xm_dry(kbdim,klev),           & !< dry air     mass in kg/m2
         xm_vap(kbdim,klev),           & !< water vapor mass in kg/m2
         xm_liq(kbdim,klev),           & !< cloud water mass in kg/m2
         xm_ice(kbdim,klev),           & !< cloud ice   mass in kg/m2
         cdnc(kbdim,klev),             & !< cloud nuclei concentration
         cld_frc(kbdim,klev),          & !< fractional cloud cover
         xm_co2(kbdim,klev),           & !< co2 mass in kg/m2
         xm_ch4(kbdim,klev),           & !< ch4 mass in kg/m2
         xm_n2o(kbdim,klev),           & !< n2o mass in kg/m2
         xm_cfc(kbdim,klev,2),         & !< cfc mass in kg/m2
         xm_o3(kbdim,klev),            & !< o3  mass in kg/m2
         xm_o2(kbdim,klev)               !< o2  mass in kg/m2
!!$         xm_trc(kbdim,klev,ktrac)        !< tracer mass mixing ratios

    REAL (wp), INTENT (OUT) ::         &
         flx_uplw    (kbdim,klev+1),   & !<   upward LW flux profile, all sky
         flx_uplw_clr(kbdim,klev+1),   & !<   upward LW flux profile, clear sky
         flx_dnlw    (kbdim,klev+1),   & !< downward LW flux profile, all sky
         flx_dnlw_clr(kbdim,klev+1),   & !< downward LW flux profile, clear sky
         flx_upsw    (kbdim,klev+1),   & !<   upward SW flux profile, all sky
         flx_upsw_clr(kbdim,klev+1),   & !<   upward SW flux profile, clear sky
         flx_dnsw    (kbdim,klev+1),   & !< downward SW flux profile, all sky
         flx_dnsw_clr(kbdim,klev+1)      !< downward SW flux profile, clear sky

    REAL (wp), INTENT (OUT) ::         &
         vis_dn_dir_sfc(kbdim)       , & !< Diffuse downward flux surface visible radiation 
         par_dn_dir_sfc(kbdim)       , & !< Diffuse downward flux surface PAR
         nir_dn_dir_sfc(kbdim)       , & !< Diffuse downward flux surface near-infrared radiation
         vis_dn_dff_sfc(kbdim)       , & !< Direct  downward flux surface visible radiation 
         par_dn_dff_sfc(kbdim)       , & !< Direct  downward flux surface PAR
         nir_dn_dff_sfc(kbdim)       , & !< Direct  downward flux surface near-infrared radiation
         vis_up_sfc    (kbdim)       , & !< Upward  flux surface visible radiation 
         par_up_sfc    (kbdim)       , & !< Upward  flux surface PAR
         nir_up_sfc    (kbdim)           !< Upward  flux surface near-infrared radiation

    ! -------------------------------------------------------------------------------------
    INTEGER  :: jk, jl, jkb,              & !< loop indicies
         icldlyr(kbdim,klev)                !< index for clear or cloudy

    REAL(wp) ::                           &
         zsemiss     (kbdim,nbndlw),      & !< LW surface emissivity by band
         pm_sfc      (kbdim)                !< surface pressure in mb
    !
    ! --- vertically reversed _vr variables
    !
    REAL(wp) ::                             &
         col_dry_vr(kbdim,klev),           & !< number of molecules/cm2 of
         pm_fl_vr  (kbdim,klev),           & !< full level pressure [mb] 
         tk_fl_vr  (kbdim,klev),           & !< full level temperature [K]
         tk_hl_vr  (kbdim,klev+1),         & !< half level temperature [K]
         cdnc_vr   (kbdim,klev),           & !< cloud nuclei concentration
         cld_frc_vr(kbdim,klev),           & !< secure cloud fraction
         ziwp_vr   (kbdim,klev),           & !< in cloud ice water content        [g/m2]
         ziwc_vr   (kbdim,klev),           & !< in cloud ice water concentration  [g/m3]
         zlwp_vr   (kbdim,klev),           & !< in cloud liquid water content        [g/m2]
         zlwc_vr   (kbdim,klev),           & !< in cloud liquid water concentration  [g/m3]
         re_drop   (kbdim,klev),           & !< effective radius of liquid
         re_cryst  (kbdim,klev),           & !< effective radius of ice
         wkl_vr    (kbdim,klev,ngas),      & !< number of molecules/cm2 of
         wx_vr     (kbdim,ncfc,klev),      & !< number of molecules/cm2 of
         cld_tau_lw_vr(kbdim,klev,nbndlw), & !< LW optical thickness of clouds
         cld_tau_sw_vr(kbdim,klev,nbndsw), & !< extincion
         cld_cg_sw_vr (kbdim,klev,nbndsw), & !< asymmetry factor
         cld_piz_sw_vr(kbdim,klev,nbndsw), & !< single scattering albedo
         aer_tau_lw_vr(kbdim,klev,nbndlw), & !< LW optical thickness of aerosols
         aer_tau_sw_vr(kbdim,klev,nbndsw), & !< aerosol optical thickness
         aer_cg_sw_vr (kbdim,klev,nbndsw), & !< aerosol asymmetry factor
         aer_piz_sw_vr(kbdim,klev,nbndsw), & !< aerosol single scattering albedo
         x_cdnc       (kbdim)                !< Scale factor for Cloud Droplet Number Concentration
    REAL(wp) ::                            &
         flx_uplw_vr(kbdim,klev+1),        & !< upward flux, total sky
         flx_uplw_clr_vr(kbdim,klev+1),    & !< upward flux, clear sky
         flx_dnlw_vr(kbdim,klev+1),        & !< downward flux, total sky
         flx_dnlw_clr_vr(kbdim,klev+1)       !< downward flux, clear sky

    !
    ! Random seeds for sampling. Needs to get somewhere upstream 
    !
    INTEGER :: rnseeds(kbdim,rng_seed_size)

    ! 1.0 Constituent properties 
    !--------------------------------
!IBM* ASSERT(NODEPS)
    DO jk = 1, klev
      jkb = klev+1-jk
      DO jl = jcs, kproma
        !
        ! --- Cloud liquid and ice mass: [kg/m2 in cell] --> [g/m2 in cloud]
        !
        cld_frc_vr(jl,jk) = MAX(EPSILON(1.0_wp),cld_frc(jl,jkb))
        ziwp_vr(jl,jk)    = xm_ice(jl,jkb)*1000.0_wp/cld_frc_vr(jl,jk)
        zlwp_vr(jl,jk)    = xm_liq(jl,jkb)*1000.0_wp/cld_frc_vr(jl,jk)
      END DO
    END DO
    !
    ! --- control for zero, infintesimal or negative cloud fractions
    !
    WHERE (cld_frc_vr(jcs:kproma,:) > 2.0_wp*EPSILON(1.0_wp))
      icldlyr(jcs:kproma,:) = 1
    ELSEWHERE
      icldlyr(jcs:kproma,:) = 0
      ziwp_vr(jcs:kproma,:) = 0.0_wp
      zlwp_vr(jcs:kproma,:) = 0.0_wp
    END WHERE
    !
    ! --- main constituent vertical reordering and unit conversion
    !
!IBM* ASSERT(NODEPS)
    DO jl = jcs, kproma
      tk_hl_vr(jl,klev+1) = tk_hl(jl,1)
      pm_sfc(jl)          = 0.01_wp*pp_sfc(jl)
    END DO

!IBM* ASSERT(NODEPS)
    wkl_vr = 0
    wx_vr = 0
    DO jk = 1, klev
      jkb = klev+1-jk
      DO jl = jcs, kproma
        !
        ! --- thermodynamic arrays
        !
        pm_fl_vr(jl,jk)   = 0.01_wp*pp_fl(jl,jkb)
        tk_hl_vr(jl,jk)   = tk_hl(jl,jkb+1)
        tk_fl_vr(jl,jk)   = tk_fl(jl,jkb)
        !
        ! --- cloud water and ice concentrations [kg/m3]
        !
        ziwc_vr(jl,jk)    = ziwp_vr(jl,jk)/dz(jl,jkb)
        zlwc_vr(jl,jk)    = zlwp_vr(jl,jk)/dz(jl,jkb)
        !
        ! --- cloud droplet number concentration  [1/m3] --> [1/cm3] ?
        !
        cdnc_vr(jl,jk)    = cdnc(jl,jkb)*1.e-6_wp
        !
        ! --- dry air: [kg/m2] --> [molecules/cm2]
        !
        col_dry_vr(jl,jk) = 0.1_wp * avo * xm_dry(jl,jkb)/amd
        !
        ! --- H2O, CO2, O3, N2O, CH4, O2: [kg/m2] --> [molecules/cm2]
        !
        wkl_vr(jl,jk,1)   = 0.1_wp * avo * xm_vap(jl,jkb)/amw
        wkl_vr(jl,jk,2)   = 0.1_wp * avo * xm_co2(jl,jkb)/amco2
        wkl_vr(jl,jk,3)   = 0.1_wp * avo * xm_ch4(jl,jkb)/amch4
        wkl_vr(jl,jk,4)   = 0.1_wp * avo * xm_o2 (jl,jkb)/amo2
        wkl_vr(jl,jk,5)   = 0.1_wp * avo * xm_o3 (jl,jkb)/amo3
        wkl_vr(jl,jk,6)   = 0.1_wp * avo * xm_n2o(jl,jkb)/amn2o
        !
        ! --- CFC11, CFC12: [kg/m2] --> [molecules/cm2]
        !
        wx_vr(jl,2,jk)    = 0.1_wp * avo * xm_cfc(jl,jkb,1)/amc11
        wx_vr(jl,3,jk)    = 0.1_wp * avo * xm_cfc(jl,jkb,2)/amc12
        !
      END DO
    END DO

    !
    ! 2.0 Surface Properties
    ! --------------------------------
    zsemiss(jcs:kproma,:) = zemiss_def 

    !
    ! 3.0 Particulate Optical Properties
    ! --------------------------------

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
      CALL set_bc_aeropt_kinne( this_datetime        ,jcs              ,&
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
           & jcs, kproma      ,kbdim                 ,klev             ,&
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
     CALL add_bc_aeropt_splumes(jg                   ,jcs              ,&
           & kproma           ,kbdim                 ,klev             ,&
           & krow             ,nbndsw                ,this_datetime    ,&
           & zf               ,dz                    ,zh(:,klev+1)     ,&
           & aer_tau_sw_vr    ,aer_piz_sw_vr         ,aer_cg_sw_vr     ,&
           & x_cdnc                                                     )
   END IF

      CALL rad_aero_diag (                jg              , &
      & krow            ,jcs             ,kproma          , &
      & kbdim           ,klev            ,nbndlw          , &
      & nbndsw          ,aer_tau_lw_vr   ,aer_tau_sw_vr   , &
      & aer_piz_sw_vr   ,aer_cg_sw_vr                       )

    CALL cloud_optics(                                  laglac         ,&
         & laland        ,jcs           ,kproma        ,kbdim          ,& 
         & klev          , ktype,&       
         & icldlyr       ,zlwp_vr       ,ziwp_vr       ,zlwc_vr        ,&
         & ziwc_vr       ,cdnc_vr       ,cld_tau_lw_vr ,cld_tau_sw_vr  ,&
         & cld_piz_sw_vr ,cld_cg_sw_vr  ,re_drop       ,re_cryst    )  

!!$    IF (locosp) THEN
!!$      DO jk=1,klev
!!$        jkb = klev+1-jk
!!$        DO jl=1,kproma
!!$          !
!!$          ! To fix - particle sizes are local variables that aren't set! 
!!$          !
!!$          cosp_reffl(jl,jk,krow) = re_drop(jl,jkb)
!!$          cosp_reffi(jl,jk,krow) = re_cryst(jl,jkb)
!!$          cosp_f3d(jl,jk,krow)   = cld_frc(jl,jk) ! Or cld_frc_vr?
!!$        END DO
!!$      END DO
!!$      IF ( Lisccp_sim ) THEN   
!!$        !don't really need vert. re-arrange ...
!!$        DO jk=1,klev
!!$          jkb = klev+1-jk
!!$          DO jl=1,kproma
!!$            cisccp_cldtau3d(jl,jk,krow) = &
!!$                 cld_tau_sw_vr(jl,jkb,9)  ! band 9 is 625 - 778 nm, needed is 670 nm
!!$          END DO
!!$        END DO
!!$        DO jk=1,klev
!!$          jkb = klev+1-jk
!!$          DO jl=1,kproma
!!$            !this is ONLY o.k. as long as wp equals dp, else conversion needed
!!$            cisccp_cldemi3d(jl,jk,krow) = &
!!$                 1._wp - exp(-1._wp*cld_tau_lw_vr(jl,jkb,6)) ! band 6 is 820 - 980 cm-1, 
!!$            !  we need 10.5 µm channel
!!$          END DO
!!$        END DO
!!$      END IF
!!$    END IF
!!$
!!$    !
!!$    ! 3.5 Interface for submodels that provide aerosol and/or cloud radiative properties:
!!$    ! -----------------------------------------------------------------------------------
!!$    IF (lanysubmodel)  CALL radiation_subm_1(                     &
!!$         kproma           ,kbdim            ,klev         ,krow  ,&
!!$         ktrac            ,iaero            ,nbndlw       ,nbndsw,&
!!$         aer_tau_sw_vr    ,aer_piz_sw_vr    ,aer_cg_sw_vr        ,&
!!$         aer_tau_lw_vr    ,......           ,xm_trc               )
    !
    ! 4.0 Radiative Transfer Routines
    ! --------------------------------
    !
    ! Seeds for random numbers come from least significant digits of pressure field 
    !
    rnseeds(jcs:kproma,1:rng_seed_size) = (pm_fl_vr(jcs:kproma,1:rng_seed_size) -  &
         int(pm_fl_vr(jcs:kproma,1:rng_seed_size)))* 1E9 + rad_perm

    IF (ltimer) CALL timer_start(timer_lrtm)
    CALL lrtm(              jcs             ,kproma                           ,&
         & kbdim           ,klev            ,pm_fl_vr        ,pm_sfc          ,&
         & tk_fl_vr        ,tk_hl_vr        ,tk_sfc          ,wkl_vr          ,&
         & wx_vr           ,col_dry_vr      ,zsemiss         ,cld_frc_vr      ,&
         & cld_tau_lw_vr   ,aer_tau_lw_vr   , &
         & rnseeds         ,&
         & flx_uplw_vr     ,flx_dnlw_vr     ,flx_uplw_clr_vr ,flx_dnlw_clr_vr )
    IF (ltimer) CALL timer_stop(timer_lrtm)
    !
    ! Reset random seeds so SW doesn't depend on what's happened in LW but is also independent
    !
    rnseeds(jcs:kproma,1:rng_seed_size) = (pm_fl_vr(jcs:kproma,rng_seed_size:1:-1) - &
         int(pm_fl_vr(jcs:kproma,rng_seed_size:1:-1)))* 1E9 + rad_perm

    ! 
    ! Potential pitfall - we're passing every argument but some may not be present
    !
    IF (ltimer) CALL timer_start(timer_srtm)
    CALL srtm(               jcs             ,kproma                           , & 
         &  kbdim           ,klev            ,pm_fl_vr        ,tk_fl_vr        , &
         &  wkl_vr          ,col_dry_vr                                        , &
         &  alb_vis_dir     ,alb_vis_dif     ,alb_nir_dir     ,alb_nir_dif     , &
         &  pmu0            ,daylght_frc     ,ssi_factor      ,psctm           , &
         &  cld_frc_vr      ,cld_tau_sw_vr   ,cld_cg_sw_vr                     , &
         &  cld_piz_sw_vr   ,aer_tau_sw_vr   ,aer_cg_sw_vr    ,aer_piz_sw_vr   , & 
         &  rnseeds         ,&
         &  flx_dnsw        , flx_upsw        ,flx_dnsw_clr    ,flx_upsw_clr   , &
         &  vis_dn_dir_sfc  ,par_dn_dir_sfc  ,nir_dn_dir_sfc                   , &
         &  vis_dn_dff_sfc  ,par_dn_dff_sfc  ,nir_dn_dff_sfc                   , &
         &  vis_up_sfc      ,par_up_sfc      ,nir_up_sfc                       )
    IF (ltimer) CALL timer_stop(timer_srtm)
    !
    ! 5.0 Post Processing
    ! --------------------------------
    !
    ! Lw fluxes are vertically reversed but SW fluxes are not
    !
    flx_uplw    (jcs:kproma,1:klev+1) = flx_uplw_vr    (jcs:kproma,klev+1:1:-1) 
    flx_uplw_clr(jcs:kproma,1:klev+1) = flx_uplw_clr_vr(jcs:kproma,klev+1:1:-1)
    flx_dnlw    (jcs:kproma,1:klev+1) = flx_dnlw_vr    (jcs:kproma,klev+1:1:-1)
    flx_dnlw_clr(jcs:kproma,1:klev+1) = flx_dnlw_clr_vr(jcs:kproma,klev+1:1:-1)

!!$    !
!!$    ! 6.0 Interface for submodel diagnosics after radiation calculation:
!!$    ! ------------------------------------------------------------------
!!$    IF (lanysubmodel)  CALL radiation_subm_2(kproma, kbdim, krow, klev, ktrac, iaero, xm_trc)
!!$
  END SUBROUTINE psrad_interface_onBlock

END MODULE mo_psrad_interface
