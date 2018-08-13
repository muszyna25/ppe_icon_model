!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
!! @author L. Linardakis, MPI-M, Hamburg
!!
!! @remarks Test the concurrent ps_rad in sequential mode
!!
MODULE mo_psrad_interface_test

  USE mo_kind,                       ONLY: wp
  USE mo_exception,                  ONLY: finish
  USE mo_model_domain,               ONLY: t_patch, p_patch
  USE mo_parallel_config            ,ONLY: nproma
  USE mtime,                         ONLY: datetime
  USE mo_psrad_general,              ONLY: ncfc, ngas, nbndsw, nbndlw  ! constants
  USE mo_psrad_interface,            ONLY: psrad_interface
  USE mo_psrad_interface_namelist,   ONLY: configure_ps_radiation_test, number_of_levels
  USE mo_psrad_interface_memory,     ONLY: t_psrad_interface, construct_psrad_interface_memory, &
    & destruct_psrad_interface_memory, psrad_interface_memory
  USE mo_master_control,             ONLY: get_my_namelist_filename
  USE mo_io_units,                   ONLY: filename_max

  IMPLICIT NONE

  PRIVATE
  LOGICAL :: is_first_call = .true.

  PUBLIC :: psrad_interface_test
  
CONTAINS

  !-----------------------------------------------------------------------------
  !>
  SUBROUTINE psrad_interface_test(                                               &
      & patch,                                                              &
      & irad_aero     ,klev                                                ,& 
      & ktype                                                              ,&
      & psctm, ssi_factor,                                                  &
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

    INTEGER,INTENT(IN)  ::             &
         irad_aero,                    & !< aerosol control
         klev,                         & !< number of levels
!!$         ktrac,                        & !< number of tracers
         ktype(:,:)                    !< type of convection

    REAL(wp),INTENT(IN) :: psctm                         !< orbit and time dependent solar constant for radiation time step
    REAL(wp),INTENT(IN) :: ssi_factor(nbndsw)            !< fraction of TSI in the 14 RRTM SW bands

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
 
    CHARACTER(len=filename_max) :: my_namelist_filename
    CHARACTER(len=filename_max) :: master_namelist_filename="icon_master.namelist"
    TYPE(t_psrad_interface),POINTER :: this_memory  !< shape: (n_dom)
 
    CHARACTER(len=*), PARAMETER :: method_name="psrad_interface_test"
    !---------------------------------------------------------

    IF (is_first_call) THEN
      my_namelist_filename = get_my_namelist_filename()
      CALL configure_ps_radiation_test(my_namelist_filename,master_namelist_filename)
      CALL construct_psrad_interface_memory( p_patch )

      this_memory => psrad_interface_memory(1)
      IF (this_memory%const%irad_aero /= irad_aero) &
        CALL finish(method_name, "this_memory%const%irad_aero /= irad_aero")
      IF (this_memory%const%no_of_levels /= klev) &
        CALL finish(method_name, "this_memory%const%no_of_levels /= klev")

      this_memory%const%zf = zf
      this_memory%const%zh = zh
      this_memory%const%dz = dz

      is_first_call = .false.
    ENDIF

    this_memory => psrad_interface_memory(1)
    
    !---------------------------------------------
    this_memory%parameterized%xm_ch4  = xm_ch4
    this_memory%parameterized%xm_n2o  = xm_n2o
    this_memory%parameterized%xm_cfc  = xm_cfc
    this_memory%parameterized%xm_o2   = xm_o2
    this_memory%parameterized%xm_o3   = xm_o3

    !---------------------------------------------
    this_memory%in%psctm =psctm
    this_memory%in%ssi_factor = ssi_factor
    this_memory%in%this_datetime => this_datetime
    this_memory%in%convection_type = ktype
    this_memory%in%loland = loland
    this_memory%in%loglac = loglac
    this_memory%in%pcos_mu0 = pcos_mu0
    this_memory%in%daylght_frc = daylght_frc
    this_memory%in%alb_vis_dir = alb_vis_dir
    this_memory%in%alb_nir_dir = alb_nir_dir
    this_memory%in%alb_vis_dif = alb_vis_dif
    this_memory%in%alb_nir_dif = alb_nir_dif
    this_memory%in%pp_sfc = pp_sfc
    this_memory%in%pp_fl = pp_fl
    this_memory%in%tk_sfc = tk_sfc
    this_memory%in%tk_fl = tk_fl
    this_memory%in%tk_hl = tk_hl
    this_memory%in%xm_dry = xm_dry
    this_memory%in%xm_vap = xm_vap
    this_memory%in%xm_liq = xm_liq
    this_memory%in%xm_ice = xm_ice
    this_memory%in%cdnc = cdnc
    this_memory%in%xc_frc = xc_frc
    this_memory%in%xm_co2 = xm_co2

    !------------------------------------------------------------------------
    CALL psrad_interface(                                                   &
      & this_memory%const%patch,                                           &
      & this_memory%const%irad_aero     ,this_memory%const%no_of_levels    ,& 
      & this_memory%in%convection_type                                     ,&
      & this_memory%in%psctm, this_memory%in%ssi_factor,                    &
      & this_memory%in%loland          ,this_memory%in%loglac          ,    &
      & this_memory%in%this_datetime                                      ,&
      & this_memory%in%pcos_mu0        ,this_memory%in%daylght_frc          ,&
      & this_memory%in%alb_vis_dir     ,this_memory%in%alb_nir_dir     ,    &
      & this_memory%in%alb_vis_dif     ,this_memory%in%alb_nir_dif     ,&
      & this_memory%const%zf, this_memory%const%zh, this_memory%const%dz   ,&
      & this_memory%in%pp_sfc          ,this_memory%in%pp_fl               ,&
      & this_memory%in%tk_sfc          ,this_memory%in%tk_fl,               &
      & this_memory%in%tk_hl                            ,&
      & this_memory%in%xm_dry          ,this_memory%in%xm_vap          ,   &
      & this_memory%in%xm_liq          ,this_memory%in%xm_ice             ,&
      & this_memory%in%cdnc            ,this_memory%in%xc_frc              ,&
      & this_memory%in%xm_co2          ,this_memory%parameterized%xm_ch4    , &
      & this_memory%parameterized%xm_n2o,this_memory%parameterized%xm_cfc   ,&
      & this_memory%parameterized%xm_o3,this_memory%parameterized%xm_o2     ,&
      & lw_upw          ,lw_upw_clr      ,lw_dnw          ,lw_dnw_clr      ,&
      & sw_upw          ,sw_upw_clr      ,sw_dnw          ,sw_dnw_clr      ,&
      & vis_dn_dir_sfc  ,par_dn_dir_sfc  ,nir_dn_dir_sfc                   ,&
      & vis_dn_dff_sfc  ,par_dn_dff_sfc  ,nir_dn_dff_sfc                   ,&
      & vis_up_sfc      ,par_up_sfc      ,nir_up_sfc                       ) 

  END SUBROUTINE psrad_interface_test
 ! -------------------------------------------------------------------------------------
 

END MODULE mo_psrad_interface_test
