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
MODULE mo_atmo_psrad_interface

  USE iso_c_binding,                 ONLY: c_loc

  USE mo_kind,                       ONLY: wp
  USE mo_exception,                  ONLY: finish, message
  USE mo_model_domain,               ONLY: t_patch, p_patch
  USE mo_parallel_config,            ONLY: nproma
  USE mtime,                         ONLY: datetime
  USE mo_master_control,             ONLY: isRestart
  USE mo_psrad_interface,            ONLY: psrad_interface
  USE mo_psrad_interface_memory,     ONLY: t_psrad_interface, construct_psrad_interface_memory, &
    & destruct_psrad_interface_memory, psrad_interface_memory
  USE mo_mpi,                        ONLY: global_mpi_barrier

  USE mo_psrad_general,              ONLY: ncfc, ngas, nbndsw, nbndlw, nmixture
  USE mo_io_units,                   ONLY: filename_max
  USE mo_master_control,             ONLY: get_my_namelist_filename, process_exists, ps_radiation_process
  USE mo_psrad_interface_namelist,   ONLY: configure_ps_radiation_test
  USE mo_echam_phy_config,           ONLY: echam_phy_tc
  USE mo_psrad_communication,        ONLY: setup_atmo_2_psrad_communication, &
                                           exchange_data_atmo_2_psrad,       &
                                           exchange_data_psrad_2_atmo,       &
                                           free_atmo_psrad_communication 
  USE mo_master_control,             ONLY: process_exists, atmo_process,     &
                                        ps_radiation_process

  USE mo_echam_phy_memory,           ONLY: t_echam_phy_field, prm_field


  USE mtime                         ,ONLY: datetime, newDatetime, deallocateDatetime, getTotalSecondsTimeDelta
  USE mo_timer,                      ONLY: ltimer, timer_start, timer_stop, &
   &  timer_extra1,  timer_extra2, timer_extra21
  IMPLICIT NONE

  PRIVATE

  PUBLIC :: atmo_psrad_interface, psrad_concurrent_interface, dtrad_shift, &
    & setup_atmo_2_psrad, finalize_atmo_radation, finalize_psrad_concurrent


  REAL(wp) :: dtrad_shift = 0.0_wp ! this holds the radation time shift due to concurrency = dt_rad

  LOGICAL :: is_first_timestep = .true.

  LOGICAL :: is_sequential_test = .false.  ! deactivates the concurrent radation, even if activated in the script, allows to compare the results to the is_concurrent_test
  LOGICAL :: is_concurrent_test = .false.  ! .false. ! runs the concurrent radation in sequential mode
  
CONTAINS

  !-----------------------------------------------------------------------------
  !>
  ! called from the atmo model
  SUBROUTINE setup_atmo_2_psrad()

    TYPE(datetime), POINTER    ::  reference_date

    dtrad_shift = 0.0_wp
    
    IF (process_exists(ps_radiation_process)) THEN
      CALL setup_atmo_2_psrad_communication()

      IF (.not. is_sequential_test .and. .not. is_concurrent_test) THEN
        ! For conversion to milliseconds, we need an anchor date,
        ! although we know that "dtime" is too small for this to be
        ! relevant.
        reference_date => newDatetime("1980-06-01T00:00:00.000")
        dtrad_shift = getTotalSecondsTimeDelta(echam_phy_tc(1)%dt_rad,reference_date)
        CALL deallocateDatetime(reference_date)
      END IF

    END IF

  END SUBROUTINE setup_atmo_2_psrad
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  !>
  ! called from the atmo model
  SUBROUTINE finalize_atmo_radation()

    CHARACTER(len=*), PARAMETER :: method_name="finalize_atmo_radation"

    TYPE(t_echam_phy_field) ,POINTER    :: field
    REAL(wp), POINTER   :: &
      & lw_dnw_clr(:,:,:),& !< Clear-sky downward longwave  at all levels
      & lw_upw_clr(:,:,:),& !< Clear-sky upward   longwave  at all levels
      & sw_dnw_clr(:,:,:),& !< Clear-sky downward shortwave at all levels
      & sw_upw_clr(:,:,:),& !< Clear-sky upward   shortwave at all levels
      & lw_dnw(:,:,:),    & !< All-sky   downward longwave  at all levels
      & lw_upw(:,:,:),    & !< All-sky   upward   longwave  at all levels
      & sw_dnw(:,:,:),    & !< All-sky   downward shortwave at all levels
      & sw_upw(:,:,:)       !< All-sky   upward   shortwave at all levels

    REAL (wp), POINTER ::         &
         vis_dn_dir_sfc(:,:)       , & !< Diffuse downward flux surface visible radiation 
         par_dn_dir_sfc(:,:)       , & !< Diffuse downward flux surface PAR
         nir_dn_dir_sfc(:,:)       , & !< Diffuse downward flux surface near-infrared radiation
         vis_dn_dff_sfc(:,:)       , & !< Direct  downward flux surface visible radiation 
         par_dn_dff_sfc(:,:)       , & !< Direct  downward flux surface PAR
         nir_dn_dff_sfc(:,:)       , & !< Direct  downward flux surface near-infrared radiation
         vis_up_sfc    (:,:)       , & !< Upward  flux surface visible radiation 
         par_up_sfc    (:,:)       , & !< Upward  flux surface PAR
         nir_up_sfc    (:,:)           !< Upward  flux surface near-infrared radiation
    INTEGER :: grid_id

    IF (process_exists(ps_radiation_process)) THEN

      IF (.not. is_sequential_test .and. .not. is_concurrent_test) THEN
        CALL message(method_name, "starts communication...")
        ! recieve the final results
        IF (ltimer) CALL timer_start(timer_extra2)


        field => prm_field(1)
        lw_dnw_clr     => field%rldcs_rt(:,:,:)    !< out  Clear-sky net longwave  at all levels
        lw_upw_clr     => field%rlucs_rt(:,:,:)    !< out  Clear-sky net longwave  at all levels
        sw_dnw_clr     => field%rsdcs_rt(:,:,:)    !< out  Clear-sky net shortwave at all levels
        sw_upw_clr     => field%rsucs_rt(:,:,:)    !< out  Clear-sky net shortwave at all levels
        lw_dnw         => field%rld_rt  (:,:,:)    !< out  All-sky net longwave  at all levels
        lw_upw         => field%rlu_rt  (:,:,:)    !< out  All-sky net longwave  at all levels
        sw_dnw         => field%rsd_rt  (:,:,:)    !< out  All-sky net longwave  at all levels
        sw_upw         => field%rsu_rt  (:,:,:)    !< out  All-sky net longwave  at all levels
      !
        vis_dn_dir_sfc => field%rvds_dir_rt(:,:)   !< out  all-sky downward direct visible radiation at surface
        par_dn_dir_sfc => field%rpds_dir_rt(:,:)   !< out  all-sky downward direct PAR     radiation at surface
        nir_dn_dir_sfc => field%rnds_dir_rt(:,:)   !< out  all-sky downward direct near-IR radiation at surface
        vis_dn_dff_sfc => field%rvds_dif_rt(:,:)   !< out  all-sky downward diffuse visible radiation at surface
        par_dn_dff_sfc => field%rpds_dif_rt(:,:)   !< out  all-sky downward diffuse PAR     radiation at surface
        nir_dn_dff_sfc => field%rnds_dif_rt(:,:)   !< out  all-sky downward diffuse near-IR radiation at surface
        vis_up_sfc     => field%rvus_rt    (:,:)   !< out  all-sky upward visible radiation at surface
        par_up_sfc     => field%rpus_rt    (:,:)   !< out  all-sky upward PAR     radiation at surfac
        nir_up_sfc     => field%rnus_rt    (:,:)   !< out  all-sky upward near-IR radiation at surface

        grid_id = 1
        CALL exchange_data_psrad_2_atmo(grid_id,  &
                      c_loc(lw_upw(1,1,1)),       &
                      c_loc(lw_upw_clr(1,1,1)),   &
                      c_loc(lw_dnw(1,1,1)),       &
                      c_loc(lw_dnw_clr(1,1,1)),   &
                      c_loc(sw_upw(1,1,1)),       &
                      c_loc(sw_upw_clr(1,1,1)),   &
                      c_loc(sw_dnw(1,1,1)),       &
                      c_loc(sw_dnw_clr(1,1,1)),   &
                      c_loc(vis_dn_dir_sfc(1,1)), &
                      c_loc(par_dn_dir_sfc(1,1)), &
                      c_loc(nir_dn_dir_sfc(1,1)), &
                      c_loc(vis_dn_dff_sfc(1,1)), &
                      c_loc(par_dn_dff_sfc(1,1)), &
                      c_loc(nir_dn_dff_sfc(1,1)), &
                      c_loc(vis_up_sfc(1,1)),     &
                      c_loc(par_up_sfc(1,1)),     &
                      c_loc(nir_up_sfc(1,1)))
        IF (ltimer) CALL timer_stop(timer_extra2)
        CALL message(method_name, "ends communication")

      END IF

!      CALL free_atmo_psrad_communication()

    END IF

  END SUBROUTINE finalize_atmo_radation
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  !>
  SUBROUTINE atmo_psrad_interface(                                          &
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

    ! OUT
    REAL(wp), INTENT(INOUT)   :: &
      & lw_dnw_clr(:,:,:),& !< Clear-sky downward longwave  at all levels
      & lw_upw_clr(:,:,:),& !< Clear-sky upward   longwave  at all levels
      & sw_dnw_clr(:,:,:),& !< Clear-sky downward shortwave at all levels
      & sw_upw_clr(:,:,:),& !< Clear-sky upward   shortwave at all levels
      & lw_dnw(:,:,:),    & !< All-sky   downward longwave  at all levels
      & lw_upw(:,:,:),    & !< All-sky   upward   longwave  at all levels
      & sw_dnw(:,:,:),    & !< All-sky   downward shortwave at all levels
      & sw_upw(:,:,:)       !< All-sky   upward   shortwave at all levels

    REAL (wp), INTENT (INOUT) ::         &
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
 
    CHARACTER(len=*), PARAMETER :: method_name="atmo_psrad_interface"

!     CALL message(method_name, " starts...")
    !---------------------------------------------------------
    IF (is_sequential_test .or. .not. process_exists(ps_radiation_process)) THEN
      !-------------------------------------------------------------------
      ! this is the original sequential call to psrad
      ! the concurrent radiation will be diactivated if is running
      CALL psrad_interface(                                                   &
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

    ELSE
        
      !-------------------------------------------------------------------
      ! run psrad concurrently
      CALL atmo_psrad_concurrent_interface(                                   &
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

    ENDIF

    is_first_timestep = .false.

!     CALL message(method_name, " ended.")

  END SUBROUTINE atmo_psrad_interface
  !-------------------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  !>
  SUBROUTINE atmo_psrad_concurrent_interface(                               &
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
         klev                            !< number of levels
    INTEGER,INTENT(IN), TARGET  ::             &
         ktype(:,:)                    !< type of convection

    REAL(wp),INTENT(IN) :: psctm                         !< orbit and time dependent solar constant for radiation time step
    REAL(wp),INTENT(IN), TARGET :: ssi_factor(nbndsw)            !< fraction of TSI in the 14 RRTM SW bands

    LOGICAL,INTENT(IN), TARGET ::              &
         loland(:,:),                & !< land sea mask, land=.true.
         loglac(:,:)                   !< glacier mask, glacier=.true.

    TYPE(datetime), POINTER ::  this_datetime !< actual time step

    REAL(WP),INTENT(IN), TARGET  ::  &
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

    ! output
    REAL(wp), INTENT(INOUT), TARGET   :: &
      & lw_dnw_clr(:,:,:),& !< Clear-sky downward longwave  at all levels
      & lw_upw_clr(:,:,:),& !< Clear-sky upward   longwave  at all levels
      & sw_dnw_clr(:,:,:),& !< Clear-sky downward shortwave at all levels
      & sw_upw_clr(:,:,:),& !< Clear-sky upward   shortwave at all levels
      & lw_dnw(:,:,:),    & !< All-sky   downward longwave  at all levels
      & lw_upw(:,:,:),    & !< All-sky   upward   longwave  at all levels
      & sw_dnw(:,:,:),    & !< All-sky   downward shortwave at all levels
      & sw_upw(:,:,:)       !< All-sky   upward   shortwave at all levels

    ! output
    REAL (wp), INTENT (INOUT), TARGET ::         &
         vis_dn_dir_sfc(:,:)       , & !< Diffuse downward flux surface visible radiation 
         par_dn_dir_sfc(:,:)       , & !< Diffuse downward flux surface PAR
         nir_dn_dir_sfc(:,:)       , & !< Diffuse downward flux surface near-infrared radiation
         vis_dn_dff_sfc(:,:)       , & !< Direct  downward flux surface visible radiation 
         par_dn_dff_sfc(:,:)       , & !< Direct  downward flux surface PAR
         nir_dn_dff_sfc(:,:)       , & !< Direct  downward flux surface near-infrared radiation
         vis_up_sfc    (:,:)       , & !< Upward  flux surface visible radiation 
         par_up_sfc    (:,:)       , & !< Upward  flux surface PAR
         nir_up_sfc    (:,:)           !< Upward  flux surface near-infrared radiation


    INTEGER, TARGET ::  irad_aero_target(1), klev_target(1)
    REAL(wp), TARGET ::  psctm_target(1)

    CHARACTER(len=*), PARAMETER :: method_name="atmo_psrad_concurrent_interface"


!     CALL message (method_name, " starts...")
    !---------------------------------------------------------
    ! this is only a test !
!     CALL global_mpi_barrier
!     CALL message (method_name, " barrier done")

! ! #ifndef NAGFOR
    IF (ltimer) CALL timer_start(timer_extra1)

   ! send to radiation
    irad_aero_target(1) = irad_aero
    klev_target(1)      = klev
    psctm_target(1)     = psctm

    CALL exchange_data_atmo_2_psrad(patch%id,  &
                  c_loc(irad_aero_target(1)),  &
                  c_loc(klev_target(1)),   &
                  c_loc(ktype(1,1)),       &
                  c_loc(psctm_target(1)),  &
                  c_loc(ssi_factor(1)),    &
                  c_loc(loland(1,1)),      &
                  c_loc(loglac(1,1)),      &
                  c_loc(pcos_mu0(1,1)),    &
                  c_loc(daylght_frc(1,1)), &
                  c_loc(alb_vis_dir(1,1)), &
                  c_loc(alb_nir_dir(1,1)), &
                  c_loc(alb_vis_dif(1,1)), &
                  c_loc(alb_nir_dif(1,1)), &
                  c_loc(zf(1,1,1)),        &
                  c_loc(zh(1,1,1)),        &
                  c_loc(dz(1,1,1)),        &
                  c_loc(pp_sfc(1,1)),      &
                  c_loc(pp_fl(1,1,1)),     &
                  c_loc(tk_sfc(1,1)),      &
                  c_loc(tk_fl(1,1,1)),     &
                  c_loc(tk_hl(1,1,1)),     &
                  c_loc(xm_dry(1,1,1)),    &
                  c_loc(xm_vap(1,1,1)),    &
                  c_loc(xm_liq(1,1,1)),    &
                  c_loc(xm_ice(1,1,1)),    &
                  c_loc(cdnc(1,1,1)),      &
                  c_loc(xc_frc(1,1,1)),    &
                  c_loc(xm_co2(1,1,1)),    &
                  c_loc(xm_ch4(1,1,1)),    &
                  c_loc(xm_n2o(1,1,1)),    &
                  c_loc(xm_cfc(1,1,1,1)),  &
                  c_loc(xm_o3(1,1,1)),     &
                  c_loc(xm_o2(1,1,1)))

    IF (ltimer) CALL timer_stop(timer_extra1)

    
    IF (.not. is_first_timestep .or. is_concurrent_test) THEN
      IF (ltimer) CALL timer_start(timer_extra2)
      ! receive from radiation
      CALL exchange_data_psrad_2_atmo(patch%id,                   &
                    c_loc(lw_upw(1,1,1)),       &
                    c_loc(lw_upw_clr(1,1,1)),   &
                    c_loc(lw_dnw(1,1,1)),       &
                    c_loc(lw_dnw_clr(1,1,1)),   &
                    c_loc(sw_upw(1,1,1)),       &
                    c_loc(sw_upw_clr(1,1,1)),   &
                    c_loc(sw_dnw(1,1,1)),       &
                    c_loc(sw_dnw_clr(1,1,1)),   &
                    c_loc(vis_dn_dir_sfc(1,1)), &
                    c_loc(par_dn_dir_sfc(1,1)), &
                    c_loc(nir_dn_dir_sfc(1,1)), &
                    c_loc(vis_dn_dff_sfc(1,1)), &
                    c_loc(par_dn_dff_sfc(1,1)), &
                    c_loc(nir_dn_dff_sfc(1,1)), &
                    c_loc(vis_up_sfc(1,1)),     &
                    c_loc(par_up_sfc(1,1)),     &
                    c_loc(nir_up_sfc(1,1)))
      IF (ltimer) CALL timer_stop(timer_extra2)

      
    ELSEIF (.not. isRestart()) THEN
      ! in the very first timestep, call the radiation also sequentially in order to get the first output values 
      ! in the case of restart the values will be read from the restart
      IF (ltimer) CALL timer_start(timer_extra21)
      CALL psrad_interface(                                                   &
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
      IF (ltimer) CALL timer_stop(timer_extra21)

    ENDIF
! ! #else
! !     CALL finish(method_name, "NAG does not compile the c_loc function")
! ! #endif

!     CALL message (method_name, " ended")

  END SUBROUTINE atmo_psrad_concurrent_interface
  !-----------------------------------------------------------------------------



  !-----------------------------------------------------------------------------
  ! this is the concurrent radiation
  !>
  SUBROUTINE psrad_concurrent_interface(mtime_current, patch)
    TYPE(datetime),  POINTER  :: mtime_current     ! current datetime (mtime)
    TYPE(t_patch),   POINTER  :: patch    ! Patch

    TYPE(t_psrad_interface),POINTER :: this_memory
    CHARACTER(len=*), PARAMETER :: method_name="psrad_concurrent_interface"

    IF (is_sequential_test) RETURN ! nothing to do here

!     CALL message (method_name, " starts...")
    !------------------------------------------------------------------------
    ! this is only a test !
!     CALL global_mpi_barrier
!     CALL message (method_name, " barrier done")

    this_memory => psrad_interface_memory(1)

!    these are filled initially
!     this_memory%const%patch      => patch
!     psrad_interface_fields%const%no_of_levels = no_of_levels

!   these are communicatted (only needed once)
!     this_memory%const%zf = 1.0_wp
!     this_memory%const%zh = 1.0_wp
!     this_memory%const%dz = 1.0_wp

    this_memory%in%this_datetime => mtime_current

! ! #ifndef NAGFOR
    IF (ltimer) CALL timer_start(timer_extra1)
    ! recieve from atmo
    CALL exchange_data_atmo_2_psrad(this_memory%const%patch%id,                       &
                                    c_loc(this_memory%const%irad_aero),               &
                                    c_loc(this_memory%const%no_of_levels),            &
                                    c_loc(this_memory%in%convection_type(1,1)),       &
                                    c_loc(this_memory%in%psctm),                      &
                                    c_loc(this_memory%in%ssi_factor(1)),              &
                                    c_loc(this_memory%in%loland(1,1)),                &
                                    c_loc(this_memory%in%loglac(1,1)),                &
                                    c_loc(this_memory%in%pcos_mu0(1,1)),              &
                                    c_loc(this_memory%in%daylght_frc(1,1)),           &
                                    c_loc(this_memory%in%alb_vis_dir(1,1)),           &
                                    c_loc(this_memory%in%alb_nir_dir(1,1)),           &
                                    c_loc(this_memory%in%alb_vis_dif(1,1)),           &
                                    c_loc(this_memory%in%alb_nir_dif(1,1)),           &
                                    c_loc(this_memory%const%zf(1,1,1)),               &
                                    c_loc(this_memory%const%zh(1,1,1)),               &
                                    c_loc(this_memory%const%dz(1,1,1)),               &
                                    c_loc(this_memory%in%pp_sfc(1,1)),                &
                                    c_loc(this_memory%in%pp_fl(1,1,1)),               &
                                    c_loc(this_memory%in%tk_sfc(1,1)),                &
                                    c_loc(this_memory%in%tk_fl(1,1,1)),               &
                                    c_loc(this_memory%in%tk_hl(1,1,1)),               &
                                    c_loc(this_memory%in%xm_dry(1,1,1)),              &
                                    c_loc(this_memory%in%xm_vap(1,1,1)),              &
                                    c_loc(this_memory%in%xm_liq(1,1,1)),              &
                                    c_loc(this_memory%in%xm_ice(1,1,1)),              &
                                    c_loc(this_memory%in%cdnc(1,1,1)),                &
                                    c_loc(this_memory%in%xc_frc(1,1,1)),              &
                                    c_loc(this_memory%in%xm_co2(1,1,1)),              &
                                    c_loc(this_memory%parameterized%xm_ch4(1,1,1)),   &
                                    c_loc(this_memory%parameterized%xm_n2o(1,1,1)),   &
                                    c_loc(this_memory%parameterized%xm_cfc(1,1,1,1)), &
                                    c_loc(this_memory%parameterized%xm_o3(1,1,1)),    &
                                    c_loc(this_memory%parameterized%xm_o2(1,1,1)))

     IF (ltimer) CALL timer_stop(timer_extra1)


    !---------------------------------------------
    IF (is_concurrent_test) THEN
       IF (ltimer) CALL timer_start(timer_extra21)
      ! run the radiation in between receiving the input and sending the output
      CALL psrad_interface(                                                  &
       & this_memory%const%patch,                                            &
       & this_memory%const%irad_aero     ,this_memory%const%no_of_levels    ,& 
       & this_memory%in%convection_type                                     ,&
       & this_memory%in%psctm, this_memory%in%ssi_factor,                    &
       & this_memory%in%loland          ,this_memory%in%loglac          ,    &
       & this_memory%in%this_datetime                                      , &
       & this_memory%in%pcos_mu0        ,this_memory%in%daylght_frc         ,&
       & this_memory%in%alb_vis_dir     ,this_memory%in%alb_nir_dir     ,    &
       & this_memory%in%alb_vis_dif     ,this_memory%in%alb_nir_dif     ,    &
       & this_memory%const%zf, this_memory%const%zh, this_memory%const%dz   ,&
       & this_memory%in%pp_sfc          ,this_memory%in%pp_fl               ,&
       & this_memory%in%tk_sfc          ,this_memory%in%tk_fl,               &
       & this_memory%in%tk_hl                            ,                   &
       & this_memory%in%xm_dry          ,this_memory%in%xm_vap          ,    &
       & this_memory%in%xm_liq          ,this_memory%in%xm_ice             , &
       & this_memory%in%cdnc            ,this_memory%in%xc_frc              ,&
       & this_memory%in%xm_co2          ,this_memory%parameterized%xm_ch4  , &
       & this_memory%parameterized%xm_n2o,this_memory%parameterized%xm_cfc  ,&
       & this_memory%parameterized%xm_o3,this_memory%parameterized%xm_o2    ,&
       & this_memory%out%lw_upw          ,this_memory%diagnostics%lw_upw_clr,&
       & this_memory%out%lw_dnw          ,this_memory%diagnostics%lw_dnw_clr,&
       & this_memory%out%sw_upw          ,this_memory%diagnostics%sw_upw_clr,&
       & this_memory%out%sw_dnw          ,this_memory%diagnostics%sw_dnw_clr,&
       & this_memory%out%vis_dn_dir_sfc  ,this_memory%out%par_dn_dir_sfc  ,  &
       & this_memory%out%nir_dn_dir_sfc                   ,                  &
       & this_memory%out%vis_dn_dff_sfc  ,this_memory%out%par_dn_dff_sfc    ,&
       & this_memory%out%nir_dn_dff_sfc,  this_memory%out%vis_up_sfc        ,&
       & this_memory%out%par_up_sfc      ,this_memory%out%nir_up_sfc         )
       IF (ltimer) CALL timer_stop(timer_extra21)
    ENDIF
    !---------------------------------------------

    ! send to atmo
    IF (.not. is_first_timestep .or. is_concurrent_test) THEN
      IF (ltimer) CALL timer_start(timer_extra2)

      CALL exchange_data_psrad_2_atmo(this_memory%const%patch%id, &
                                    c_loc(this_memory%out%lw_upw(1,1,1)), &
                                    c_loc(this_memory%diagnostics%lw_upw_clr(1,1,1)), &
                                    c_loc(this_memory%out%lw_dnw(1,1,1)), &
                                    c_loc(this_memory%diagnostics%lw_dnw_clr(1,1,1)), &
                                    c_loc(this_memory%out%sw_upw(1,1,1)), &
                                    c_loc(this_memory%diagnostics%sw_upw_clr(1,1,1)), &
                                    c_loc(this_memory%out%sw_dnw(1,1,1)), &
                                    c_loc(this_memory%diagnostics%sw_dnw_clr(1,1,1)), &
                                    c_loc(this_memory%out%vis_dn_dir_sfc(1,1)), &
                                    c_loc(this_memory%out%par_dn_dir_sfc(1,1)), &
                                    c_loc(this_memory%out%nir_dn_dir_sfc(1,1)), &
                                    c_loc(this_memory%out%vis_dn_dff_sfc(1,1)), &
                                    c_loc(this_memory%out%par_dn_dff_sfc(1,1)), &
                                    c_loc(this_memory%out%nir_dn_dff_sfc(1,1)), &
                                    c_loc(this_memory%out%vis_up_sfc(1,1)), &
                                    c_loc(this_memory%out%par_up_sfc(1,1)), &
                                    c_loc(this_memory%out%nir_up_sfc(1,1)))

       IF (ltimer) CALL timer_stop(timer_extra2)
   
    ENDIF
 
    !---------------------------------------------
    IF ( .not. is_concurrent_test) THEN
       IF (ltimer) CALL timer_start(timer_extra21)
     ! run the radiation after receiving the input and sending the output
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
       & this_memory%parameterized%xm_o3,this_memory%parameterized%xm_o2    ,&
       & this_memory%out%lw_upw          ,this_memory%diagnostics%lw_upw_clr      ,  &
       & this_memory%out%lw_dnw          ,this_memory%diagnostics%lw_dnw_clr        ,&
       & this_memory%out%sw_upw          ,this_memory%diagnostics%sw_upw_clr      ,  &
       & this_memory%out%sw_dnw          ,this_memory%diagnostics%sw_dnw_clr      ,&
       & this_memory%out%vis_dn_dir_sfc  ,this_memory%out%par_dn_dir_sfc  ,&
       & this_memory%out%nir_dn_dir_sfc                   ,&
       & this_memory%out%vis_dn_dff_sfc  ,this_memory%out%par_dn_dff_sfc  ,&
       & this_memory%out%nir_dn_dff_sfc,  this_memory%out%vis_up_sfc      ,&
       & this_memory%out%par_up_sfc      ,this_memory%out%nir_up_sfc                       )

      IF (ltimer) CALL timer_stop(timer_extra21)
    ENDIF
    !---------------------------------------------
! ! #else
! !     CALL finish(method_name, "NAG does not compile the c_loc function")
! ! #endif

!     CALL message (method_name, " ended")

    is_first_timestep = .false.

  END SUBROUTINE psrad_concurrent_interface
  !-----------------------------------------------------------------------------
 
  !-----------------------------------------------------------------------------
  !>
  SUBROUTINE finalize_psrad_concurrent(patch)
    TYPE(t_patch),   POINTER  :: patch    ! Patch

    TYPE(t_psrad_interface),POINTER :: this_memory
    CHARACTER(len=*), PARAMETER :: method_name="finalize_psrad_concurrent"

    this_memory => psrad_interface_memory(1)

    !---------------------------------------------
    IF (.not. is_sequential_test .and. .not. is_concurrent_test) THEN
      CALL message(method_name, "starts communication...")
      ! send to atmo
      IF (ltimer) CALL timer_start(timer_extra2)

      CALL exchange_data_psrad_2_atmo(this_memory%const%patch%id, &
                c_loc(this_memory%out%lw_upw(1,1,1)), &
                c_loc(this_memory%diagnostics%lw_upw_clr(1,1,1)), &
                c_loc(this_memory%out%lw_dnw(1,1,1)), &
                c_loc(this_memory%diagnostics%lw_dnw_clr(1,1,1)), &
                c_loc(this_memory%out%sw_upw(1,1,1)), &
                c_loc(this_memory%diagnostics%sw_upw_clr(1,1,1)), &
                c_loc(this_memory%out%sw_dnw(1,1,1)), &
                c_loc(this_memory%diagnostics%sw_dnw_clr(1,1,1)), &
                c_loc(this_memory%out%vis_dn_dir_sfc(1,1)), &
                c_loc(this_memory%out%par_dn_dir_sfc(1,1)), &
                c_loc(this_memory%out%nir_dn_dir_sfc(1,1)), &
                c_loc(this_memory%out%vis_dn_dff_sfc(1,1)), &
                c_loc(this_memory%out%par_dn_dff_sfc(1,1)), &
                c_loc(this_memory%out%nir_dn_dff_sfc(1,1)), &
                c_loc(this_memory%out%vis_up_sfc(1,1)), &
                c_loc(this_memory%out%par_up_sfc(1,1)), &
                c_loc(this_memory%out%nir_up_sfc(1,1)))

       IF (ltimer) CALL timer_stop(timer_extra2)
   
      CALL message(method_name, "ends communication.")
    ENDIF 
    !---------------------------------------------

!     CALL free_atmo_psrad_communication()

  END SUBROUTINE finalize_psrad_concurrent
  !-----------------------------------------------------------------------------
  
 
  !-----------------------------------------------------------------------------
  ! tests in sequantial mode the interface (not used any more, but kept for future testing of changes)
  ! create the psrad model and run through its interface sequentially
  !>
  SUBROUTINE psrad_interface_test(                                          &
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
    TYPE(t_psrad_interface),POINTER :: this_memory  

    LOGICAL, SAVE :: is_first_call = .true.
 
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
      & this_memory%parameterized%xm_o3,this_memory%parameterized%xm_o2    ,&
      & this_memory%out%lw_upw          ,this_memory%diagnostics%lw_upw_clr      ,  &
      & this_memory%out%lw_dnw          ,this_memory%diagnostics%lw_dnw_clr        ,&
      & this_memory%out%sw_upw          ,this_memory%diagnostics%sw_upw_clr      ,  &
      & this_memory%out%sw_dnw          ,this_memory%diagnostics%sw_dnw_clr      ,&
      & this_memory%out%vis_dn_dir_sfc  ,this_memory%out%par_dn_dir_sfc  ,&
      & this_memory%out%nir_dn_dir_sfc                   ,&
      & this_memory%out%vis_dn_dff_sfc  ,this_memory%out%par_dn_dff_sfc  ,&
      & this_memory%out%nir_dn_dff_sfc,  this_memory%out%vis_up_sfc      ,&
      & this_memory%out%par_up_sfc      ,this_memory%out%nir_up_sfc                       ) 

      lw_dnw = this_memory%out%lw_dnw !< All-sky   downward longwave  at all levels
      lw_upw = this_memory%out%lw_upw!< All-sky   upward   longwave  at all levels
      sw_dnw = this_memory%out%sw_dnw!< All-sky   downward shortwave at all levels
      sw_upw = this_memory%out%sw_upw!< All-sky   upward   shortwave at all levels

      ! 2D 
      vis_dn_dir_sfc = this_memory%out%vis_dn_dir_sfc !< Direct downward flux surface visible radiation 
      par_dn_dir_sfc = this_memory%out%par_dn_dir_sfc !< Direct downward flux surface PAR
      nir_dn_dir_sfc = this_memory%out%nir_dn_dir_sfc !< Direct downward flux surface near-infrared radiation
      vis_dn_dff_sfc = this_memory%out%vis_dn_dff_sfc !< Diffuse  downward flux surface visible radiation 
      par_dn_dff_sfc = this_memory%out%par_dn_dff_sfc !< Diffuse downward flux surface PAR
      nir_dn_dff_sfc = this_memory%out%nir_dn_dff_sfc !< Diffuse downward flux surface near-infrared radiation
      vis_up_sfc     = this_memory%out%vis_up_sfc     !< Upward  flux surface visible radiation 
      par_up_sfc     = this_memory%out%par_up_sfc     !< Upward  flux surface PAR
      nir_up_sfc     = this_memory%out%nir_up_sfc     !< Upward  flux surface near-infrared radiation

      ! diagnostics !
      lw_dnw_clr = this_memory%diagnostics%lw_dnw_clr !< Clear-sky downward longwave  at all levels
      lw_upw_clr = this_memory%diagnostics%lw_upw_clr !< Clear-sky upward   longwave  at all levels
      sw_dnw_clr = this_memory%diagnostics%sw_dnw_clr !< Clear-sky downward shortwave at all levels
      sw_upw_clr = this_memory%diagnostics%sw_upw_clr !< Clear-sky upward   shortwave at all levels

  END SUBROUTINE psrad_interface_test
 ! -------------------------------------------------------------------------------------

END MODULE mo_atmo_psrad_interface
