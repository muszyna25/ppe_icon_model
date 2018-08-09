!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_psrad_interface
  USE mo_kind,                       ONLY: wp
  USE mo_physical_constants,         ONLY: avo, amd, amw, amco2, amch4, &
                                           amn2o, amo3, amo2, amc11, amc12, &
                                           zemiss_def
  USE mo_exception,                  ONLY: finish, message, warning, message_text
  USE mo_model_domain,               ONLY: t_patch
  USE mo_parallel_config,            ONLY: nproma
  USE mo_impl_constants,             ONLY: min_rlcell_int, grf_bdywidth_c
  USE mo_loopindices,                ONLY: get_indices_c
  USE mo_psrad_radiation_parameters, ONLY: rad_perm
  USE mo_psrad_general,              ONLY: ncfc, ngas, nbndsw, nbndlw, nmixture, finish_cb, message_cb, warning_cb
  USE mo_psrad_cloud_optics,         ONLY: cloud_optics
  USE mo_bc_aeropt_kinne,            ONLY: set_bc_aeropt_kinne  
  USE mo_bc_aeropt_stenchikov,       ONLY: add_bc_aeropt_stenchikov 
  USE mo_bc_aeropt_splumes,          ONLY: add_bc_aeropt_splumes
  USE mo_psrad_gas_optics,           ONLY: precomputation
  USE mo_psrad_lrtm_driver,          ONLY: lrtm
  USE mo_psrad_srtm_driver,          ONLY: srtm, srtm_diags
  USE mo_psrad_random,               ONLY: seed_size
  USE mo_rad_diag,                   ONLY: rad_aero_diag
  USE mtime,                         ONLY: datetime
  USE mo_timer,                      ONLY: ltimer, timer_start, timer_stop, &
   &                                       timer_lrtm, timer_srtm

  USE mo_radiation_config,           ONLY: lrad_aero_diag

  USE mo_psrad_setup,                ONLY : psrad_basic_setup
  USE mo_namelist,                   ONLY: open_nml, position_nml, close_nml, POSITIONED
  USE mo_io_units,                   ONLY: nnml, nnml_output
  USE mo_echam_cld_config,           ONLY: echam_cld_config
  USE mo_run_config,                 ONLY: nlev
  USE mo_mpi,                        ONLY: my_process_is_stdio
  USE mo_restart_namelist,    ONLY: open_tmpfile, store_and_close_namelist

#ifdef PSRAD_TIMING
  USE mo_timer,                      ONLY: timer_rrtm_coeffs,   &
                                           timer_cloud_optics,  &
                                           timer_psrad_scaling, &
                                           timer_psrad_aerosol
#endif

! These need to be sent once in the init phae from the atmo to the ps_rad
!          zf(kbdim,klev),               & !< geometric height at full level in m
!          zh(kbdim,klev+1),             & !< geometric height at half level in m
!          dz(kbdim,klev),               & !< geometric height thickness in m

  IMPLICIT NONE

  PRIVATE

  REAL(wp), PARAMETER :: pressure_scale = 100._wp,         &
                         inverse_pressure_scale = 0.01_wp, &
                         droplet_scale = 1.0e2

  PUBLIC :: psrad_interface, pressure_scale, inverse_pressure_scale, &
            droplet_scale, setup_psrad_radiation
  
CONTAINS
  !-------------------------------------------------------------------
  SUBROUTINE setup_psrad_radiation(file_name)

    CHARACTER(len=*), INTENT(IN)      :: file_name
    INTEGER :: istat, funit

    NAMELIST /psrad_nml/ rad_perm

    ! 1.0 Read psrad_nml namelist 
    ! --------------------------------
    CALL open_nml(TRIM(file_name))
    CALL position_nml ('psrad_nml', status=istat)
    SELECT CASE (istat)
      CASE (POSITIONED) 
        READ (nnml, psrad_nml)
    END SELECT
    CALL close_nml
    ! store namelist for restart
    IF(my_process_is_stdio()) THEN
      funit = open_tmpfile()
      WRITE(funit,NML=psrad_nml)
      CALL store_and_close_namelist(funit, 'psrad_nml')
    END IF
    ! write the contents of the namelist to an ASCII file
    IF (my_process_is_stdio()) THEN
      WRITE(nnml_output,nml=psrad_nml)
    END IF
    CALL psrad_basic_setup(.false., nlev, pressure_scale, droplet_scale, &
     & echam_cld_config(1)%cinhoml1 ,echam_cld_config(1)%cinhoml2, &
     & echam_cld_config(1)%cinhoml3 ,echam_cld_config(1)%cinhomi)

    finish_cb  => finish
    message_cb => warning
    warning_cb => warning

  END SUBROUTINE setup_psrad_radiation
  !-------------------------------------------------------------------

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
#ifdef __INTEL_COMPILER
!DIR$ OPTIMIZE:1
#endif
     !-------------------------------------------------------------------

    TYPE(t_patch)   ,TARGET ,INTENT(in)    :: patch

    INTEGER,INTENT(IN)  :: &
         irad_aero,        & !< aerosol control
         klev,             & !< number of levels
         ktype(:,:)          !< type of convection

    REAL(wp),INTENT(IN) :: psctm                         !< orbit and time dependent solar constant for radiation time step
    REAL(wp),INTENT(IN) :: ssi_factor(nbndsw)            !< fraction of TSI in the 14 RRTM SW bands

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
    INTEGER  :: klevp1

    jg         = patch%id
    rl_start   = grf_bdywidth_c+1
    rl_end     = min_rlcell_int
 
    i_nchdom   = MAX(1,patch%n_childdom)
    i_startblk = patch%cells%start_blk(rl_start,1)
    i_endblk   = patch%cells%end_blk(rl_end,i_nchdom)

    klevp1     = klev+1

!$OMP PARALLEL DO PRIVATE(jcs,jce)
    DO jb = i_startblk,i_endblk
       
      CALL get_indices_c(patch, jb,i_startblk,i_endblk, jcs,jce, rl_start, rl_end)

      IF (jcs>jce) CYCLE

      IF (jcs==1) THEN
         !
         CALL psrad_interface_onBlock(                                          &
            & jg,                   jb,                                         &
            &                       jce,                                        &
            & nproma,               klev,                                       &
            !
            & irad_aero,            ktype(:,jb),                                &
            & psctm,                ssi_factor,                                 &
            & loland(:,jb),         loglac(:,jb),         this_datetime,        &
            & pcos_mu0(:,jb),       daylght_frc(:,jb),                          &
            & alb_vis_dir(:,jb),    alb_nir_dir(:,jb),                          &
            & alb_vis_dif(:,jb),    alb_nir_dif(:,jb),                          &
            & zf(:,:,jb),           zh(:,:,jb),           dz(:,:,jb),           &
            & pp_sfc(:,jb),         pp_fl(:,:,jb),                              &
            & tk_sfc(:,jb),         tk_fl(:,:,jb),        tk_hl(:,:,jb),        &
            & xm_dry(:,:,jb),       xm_vap(:,:,jb),       xm_liq(:,:,jb),       &
            & xm_ice(:,:,jb),       cdnc(:,:,jb),         xc_frc(:,:,jb),       &
            & xm_co2(:,:,jb),       xm_ch4(:,:,jb),       xm_n2o (:,:,jb),      &
            & xm_cfc(:,:,:,jb),     xm_o3(:,:,jb),        xm_o2(:,:,jb),        &
            !
            & lw_upw(:,:,jb),       lw_upw_clr (:,:,jb),                        &
            & lw_dnw(:,:,jb),       lw_dnw_clr (:,:,jb),                        &
            & sw_upw(:,:,jb),       sw_upw_clr(:,:,jb),                         &
            & sw_dnw(:,:,jb),       sw_dnw_clr(:,:,jb),                         &
            & vis_dn_dir_sfc(:,jb), par_dn_dir_sfc(:,jb), nir_dn_dir_sfc(:,jb), &
            & vis_dn_dff_sfc(:,jb), par_dn_dff_sfc(:,jb), nir_dn_dff_sfc(:,jb), &
            & vis_up_sfc(:,jb),     par_up_sfc(:,jb),     nir_up_sfc(:,jb)      )
         !
      ELSE
         !
         CALL shift_and_call_psrad_interface_onBlock(                           &
            & jg,                   jb,                                         &
            & jcs,                  jce,                                        &
            & nproma,               klev,                 klevp1,               &
            !
            & irad_aero,            ktype(:,jb),                                &
            & psctm,                ssi_factor,                                 &
            & loland(:,jb),         loglac(:,jb),         this_datetime,        &
            & pcos_mu0(:,jb),       daylght_frc(:,jb),                          &
            & alb_vis_dir(:,jb),    alb_nir_dir(:,jb),                          &
            & alb_vis_dif(:,jb),    alb_nir_dif(:,jb),                          &
            & zf(:,:,jb),           zh(:,:,jb),           dz(:,:,jb),           &
            & pp_sfc(:,jb),         pp_fl(:,:,jb),                              &
            & tk_sfc(:,jb),         tk_fl(:,:,jb),        tk_hl(:,:,jb),        &
            & xm_dry(:,:,jb),       xm_vap(:,:,jb),       xm_liq(:,:,jb),       &
            & xm_ice(:,:,jb),       cdnc(:,:,jb),         xc_frc(:,:,jb),       &
            & xm_co2(:,:,jb),       xm_ch4(:,:,jb),       xm_n2o (:,:,jb),      &
            & xm_cfc(:,:,:,jb),     xm_o3(:,:,jb),        xm_o2(:,:,jb),        &
            !
            & lw_upw(:,:,jb),       lw_upw_clr (:,:,jb),                        &
            & lw_dnw(:,:,jb),       lw_dnw_clr (:,:,jb),                        &
            & sw_upw(:,:,jb),       sw_upw_clr(:,:,jb),                         &
            & sw_dnw(:,:,jb),       sw_dnw_clr(:,:,jb),                         &
            & vis_dn_dir_sfc(:,jb), par_dn_dir_sfc(:,jb), nir_dn_dir_sfc(:,jb), &
            & vis_dn_dff_sfc(:,jb), par_dn_dff_sfc(:,jb), nir_dn_dff_sfc(:,jb), &
            & vis_up_sfc(:,jb),     par_up_sfc(:,jb),     nir_up_sfc(:,jb)      )
         !
      END IF

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

  SUBROUTINE psrad_interface_onBlock(                        &
       & jg,             jb,                                 &
       &                 kproma,                             &
       & kbdim,          klev,                               &
       & iaero,          ktype,                              &
       & psctm,          ssi_factor,                         &
       & laland,         laglac,         this_datetime,      &
       & pcos_mu0,       daylght_frc,                        &
       & alb_vis_dir,    alb_nir_dir,                        &
       & alb_vis_dif,    alb_nir_dif,                        &
       & zf,             zh,             dz,                 &
       & pp_sfc,         pp_fl,                              &
       & tk_sfc,         tk_fl,          tk_hl,              &
       & xm_dry,         xm_vap,         xm_liq,             &
       & xm_ice,         cdnc,           cld_frc,            &
       & xm_co2,         xm_ch4,         xm_n2o ,            &
       & xm_cfc ,        xm_o3,          xm_o2,              &
       & flx_uplw,       flx_uplw_clr,                       &
       & flx_dnlw,       flx_dnlw_clr,                       &
       & flx_upsw,       flx_upsw_clr,                       &
       & flx_dnsw,       flx_dnsw_clr,                       &
       & vis_dn_dir_sfc, par_dn_dir_sfc, nir_dn_dir_sfc,     &
       & vis_dn_dff_sfc, par_dn_dff_sfc, nir_dn_dff_sfc,     &
       & vis_up_sfc,     par_up_sfc,     nir_up_sfc          )

#ifdef __INTEL_COMPILER
!DIR$ OPTIMIZE:1
#endif

    INTEGER,INTENT(IN)  :: &
         jg,               & !< domain index
         jb,               & !< first dimension of 2-d arrays
         kproma,           & !< number of longitudes
         kbdim,            & !< first dimension of 2-d arrays
         klev,             & !< number of levels
         iaero,            & !< aerosol control
         ktype(:)            !< type of convection

    REAL(wp),INTENT(IN) :: psctm                         !< orbit and time dependent solar constant for radiation time step
    REAL(wp),INTENT(IN) :: ssi_factor(nbndsw)            !< fraction of TSI in the 14 RRTM SW bands

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
         zf(:,:),       & !< geometric height at full level in m
         zh(:,:),     & !< geometric height at half level in m
         dz(:,:),       & !< geometric height thickness in m
         pp_sfc(:),        & !< surface pressure in Pa
         pp_fl(:,:),    & !< full level pressure in Pa
         tk_sfc(:),        & !< surface temperature in K
         tk_fl(:,:),    & !< full level temperature in K
         tk_hl(:,:),  & !< half level temperature in K
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
         cld_frc_loc(kbdim,klev), & !< secure cloud fraction
         ziwp_loc   (kbdim,klev), & !< in cloud ice water content       [g/m2]
         ziwc_loc   (kbdim,klev), & !< in cloud ice water concentration [g/m3]
         zlwp_loc   (kbdim,klev), & !< in cloud liquid water content    [g/m2]
         zlwc_loc   (kbdim,klev), & !< in cloud liquid water concentration  [g/m3]
         wx_loc        (kbdim,klev,ncfc),& !< number of molecules/cm2 of
         cld_tau_lw_loc(kbdim,klev,nbndlw), & !< LW optical thickness of clouds
         cld_tau_sw_loc(kbdim,klev,nbndsw), & !< extincion
         cld_cg_sw_loc (kbdim,klev,nbndsw), & !< asymmetry factor
         cld_piz_sw_loc(kbdim,klev,nbndsw), & !< single scattering albedo
         aer_tau_lw_loc(kbdim,klev,nbndlw), & !< LW optical thickness of aerosols
         aer_tau_sw_loc(kbdim,klev,nbndsw), & !< aerosol optical thickness
         aer_cg_sw_loc (kbdim,klev,nbndsw), & !< aerosol asymmetry factor
         aer_piz_sw_loc(kbdim,klev,nbndsw) !< aerosol single scattering albedo

    REAL(wp), TARGET ::                   &
         gases         (kbdim,klev,ngas)  !< number of molecules/cm2 of
    ! --- vertically reversed _loc variables
    REAL(wp) ::                             &
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
    INTEGER :: rnseeds(kbdim,seed_size)

    REAL(wp), TARGET :: actual_scaleminorn2(KBDIM,klev)
    REAL(wp), TARGET :: actual_ratio(KBDIM,2,klev,nmixture) 
    REAL(wp), POINTER :: ratio(:,:,:,:), &
      scaleminorn2(:,:)
    REAL(wp) :: fac(KBDIM,2,2,klev)
    INTEGER :: laytrop(KBDIM), & !< tropopause layer index
      iabs(KBDIM,2,2,klev)
    INTEGER, DIMENSION(KBDIM,klev) :: jp, indminor
    REAL(wp), DIMENSION(KBDIM,klev,2) :: h2o_factor,h2o_fraction
    INTEGER, DIMENSION(KBDIM,klev,2) :: h2o_index
    REAL(wp), DIMENSION(KBDIM,klev) :: colbrd, colmol, &
      minorfrac, scaleminor

    ratio => actual_ratio
    scaleminorn2 => actual_scaleminorn2
    !
    ! Number of g-points per time step. Determine here to allow automatic array allocation in 
    !   lrtm, srtm subroutines. 
    !
    !INTEGER :: n_gpts_ts

    ! 1.0 Constituent properties 
    !--------------------------------
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
    icldlyr_loc(1:kproma,:) = 1
    WHERE (cld_frc_loc(1:kproma,:) <= 2.0_wp*EPSILON(1.0_wp))
      icldlyr_loc(1:kproma,:) = 0
      ziwp_loc(1:kproma,:) = 0.0_wp
      zlwp_loc(1:kproma,:) = 0.0_wp
    END WHERE
    !
    ! --- main constituent vertical reordering and unit conversion
    !
    DO jk = 1, klev
      DO jl = 1, kproma
        !
        ! --- cloud water and ice concentrations [kg/m3]
        !
        ziwc_loc(jl,jk)    = ziwp_loc(jl,jk)/dz(jl,jk)
        zlwc_loc(jl,jk)    = zlwp_loc(jl,jk)/dz(jl,jk)
        !
        ! --- dry air: [kg/m2] --> [molecules/cm2]
        !
        col_dry_loc(jl,jk) = 0.1_wp * avo * xm_dry(jl,jk)/amd
        !
        ! --- H2O, CO2, O3, N2O, CH4, O2: [kg/m2] --> [molecules/cm2]
        !
        gases(jl,jk,:)   = 0.0_wp
        gases(jl,jk,1)   = 0.1_wp * avo * xm_vap(jl,jk)/amw
        gases(jl,jk,2)   = 0.1_wp * avo * xm_co2(jl,jk)/amco2
        gases(jl,jk,3)   = 0.1_wp * avo * xm_ch4(jl,jk)/amch4
        gases(jl,jk,4)   = 0.1_wp * avo * xm_o2 (jl,jk)/amo2
        gases(jl,jk,5)   = 0.1_wp * avo * xm_o3 (jl,jk)/amo3
        gases(jl,jk,6)   = 0.1_wp * avo * xm_n2o(jl,jk)/amn2o
        !
        ! --- CFC11, CFC12: [kg/m2] --> [molecules/cm2]
        !
        wx_loc(jl,jk,:)    = 0.0_wp
        wx_loc(jl,jk,2)    = 0.1_wp * avo * xm_cfc(jl,jk,1)/amc11
        wx_loc(jl,jk,3)    = 0.1_wp * avo * xm_cfc(jl,jk,2)/amc12
        !
      END DO
    END DO
    CALL flip_ud(kproma, pp_fl, pm_fl_vr)
    CALL flip_ud(kproma, cld_frc_loc, cld_frc_vr)

    ! 2.0 Surface Properties
    ! --------------------------------
    zsemiss(1:kproma,:) = zemiss_def 
    !
    ! 3.0 Particulate Optical Properties
    ! --------------------------------
#ifdef PSRAD_TIMING
    IF (ltimer) CALL timer_stop(timer_psrad_scaling)
#endif
#ifdef PSRAD_TIMING
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
      CALL set_bc_aeropt_kinne(this_datetime,                           &
           & jg,                                                        &
           & 1, kproma,      kbdim,                 klev,                  &
           & jb,               nbndsw,                nbndlw,         &
           & zf,               dz,                                      &
           & aer_tau_sw_vr,    aer_piz_sw_vr,         aer_cg_sw_vr,     &
           & aer_tau_lw_vr                                              )
    END IF
    IF (iaero==14 .OR. iaero==15 .OR. iaero==18) THEN
! iaero=14: only Stechnikov's volcanic aerosols are used (added to zero)
! iaero=15: Stenchikov's volcanic aerosols are added to Kinne aerosols
! iaero=18: Stenchikov's volcanic aerosols are added to Kinne background
!           aerosols (of natural origin, 1850) 
      CALL add_bc_aeropt_stenchikov(this_datetime,    jg,               &
           & 1, kproma,           kbdim,                 klev,             &
           & jb,             nbndsw,                nbndlw,           &
           & dz,               pp_fl,                                   &
           & aer_tau_sw_vr,    aer_piz_sw_vr,         aer_cg_sw_vr,     &
           & aer_tau_lw_vr                                              )
    END IF
!!$    IF (iaero==16) THEN
!!$      CALL add_aop_volc_ham( &
!!$           & kproma,           kbdim,                 klev,             &
!!$           & jb,             nbndlw,                nbndsw,           &
!!$           & aer_tau_lw_vr,    aer_tau_sw_vr,         aer_piz_sw_vr,    &
!!$           & aer_cg_sw_vr                                               )
!!$    END IF
!!$    IF (iaero==17) THEN
!!$      CALL add_aop_volc_crow( &
!!$           & kproma,           kbdim,                 klev,             &
!!$           & jb,             nbndlw,                nbndsw,           &
!!$           & aer_tau_lw_vr,    aer_tau_sw_vr,         aer_piz_sw_vr,    &
!!$           & aer_cg_sw_vr                                               )
!!$    END IF
    IF (iaero==18) THEN
! iaero=18: Simple plumes are added to Stenchikov's volcanic aerosols 
!           and Kinne background aerosols (of natural origin, 1850) 
      CALL add_bc_aeropt_splumes(jg,                                     &
           & 1, kproma,           kbdim,                 klev,             &
           & jb,             nbndsw,                this_datetime,    &
           & zf,               dz,                    zh(:,klev+1),     &
           & aer_tau_sw_vr,    aer_piz_sw_vr,         aer_cg_sw_vr,     &
           & x_cdnc                                                     )
    END IF

    ! this should be decativated in the concurrent version and make the aer_* global variables for output
    IF (lrad_aero_diag) THEN
      CALL rad_aero_diag (                                  &
        & jg,              jb,         1, kproma,           &
        & kbdim,           klev,            nbndlw,           &
        & nbndsw,          aer_tau_lw_vr,   aer_tau_sw_vr,    &
        & aer_piz_sw_vr,   aer_cg_sw_vr                       )
    ENDIF

    DO jl = 1,nbndlw
      CALL flip_ud(kproma, aer_tau_lw_vr(:,:,jl), aer_tau_lw_loc(:,:,jl))
    ENDDO
    DO jl = 1,nbndsw
      CALL flip_ud(kproma, aer_tau_sw_vr(:,:,jl), aer_tau_sw_loc(:,:,jl))
      CALL flip_ud(kproma, aer_cg_sw_vr(:,:,jl), aer_cg_sw_loc(:,:,jl))
      CALL flip_ud(kproma, aer_piz_sw_vr(:,:,jl), aer_piz_sw_loc(:,:,jl))
    ENDDO


#ifdef PSRAD_TIMING
    IF (ltimer) CALL timer_stop(timer_psrad_aerosol)
#endif

#ifdef PSRAD_TIMING
    IF (ltimer) CALL timer_start(timer_cloud_optics)
#endif
    CALL cloud_optics(                                                     &
         & laglac,         laland,         kproma,         kbdim,          &
         & klev,           ktype,                                          &
         & icldlyr_loc,    zlwp_loc,       ziwp_loc,       zlwc_loc,       &
         & ziwc_loc,       cdnc,           cld_tau_lw_loc, cld_tau_sw_loc, &
         & cld_piz_sw_loc, cld_cg_sw_loc,  re_drop,        re_cryst        )
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
    rnseeds(1:kproma,1:seed_size) = &
      int((inverse_pressure_scale * pm_fl_vr(1:kproma,1:seed_size) -  &
      int(inverse_pressure_scale * pm_fl_vr(1:kproma,1:seed_size)))* 1E9 + rad_perm)
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
    CALL precomputation(kproma,      kbdim,         klev,        &
         & .false.,   pp_fl,         tk_fl,         col_dry_loc, &
         & gases,     laytrop,       jp,            iabs,        &
         & colbrd,    colmol,        fac,                        &
         & ratio,     h2o_factor,    h2o_fraction,  h2o_index,   &
         & minorfrac, scaleminor,    scaleminorn2,  indminor     )
#ifdef PSRAD_TIMING
    IF (ltimer) CALL timer_stop(timer_rrtm_coeffs)
#endif
    IF (ltimer) CALL timer_start(timer_lrtm)
    CALL lrtm(kproma,      kbdim,           klev,                        &
         & pp_fl,          pp_sfc,                                       &
         & tk_fl,          tk_hl,          tk_sfc,                       &
         & gases, wx_loc,  col_dry_loc,    zsemiss,        cld_frc_loc,  &
         & cld_tau_lw_loc, aer_tau_lw_loc, rnseeds,                      &
         & ratio,          scaleminorn2,   fac,                          &
         & laytrop,        iabs,           jp,             indminor,     &
         & h2o_factor,     h2o_fraction,   h2o_index,      colbrd,       &
         & minorfrac,      scaleminor,                                  &
         & flx_uplw,       flx_dnlw,       flx_uplw_clr,   flx_dnlw_clr )
    IF (ltimer) CALL timer_stop(timer_lrtm)

    !
    ! Reset random seeds so SW doesn't depend on what's happened in LW but is also independent
    !
    rnseeds(1:kproma,1:seed_size) = &
      int((inverse_pressure_scale * pm_fl_vr(1:kproma,seed_size:1:-1) - &
      int(inverse_pressure_scale * pm_fl_vr(1:kproma,seed_size:1:-1)))* 1E9 + rad_perm)

    ! Potential pitfall - we're passing every argument but some may not be present
    IF (ltimer) CALL timer_start(timer_srtm)
    CALL srtm(kproma,       kbdim,          klev,                           &
         &  alb_vis_dir,    alb_vis_dif,    alb_nir_dir,   alb_nir_dif,     &
         &  pcos_mu0,       daylght_frc,    ssi_factor,    psctm,           &
         &  cld_frc_loc,    cld_tau_sw_loc, cld_cg_sw_loc,                  &
         &  cld_piz_sw_loc, aer_tau_sw_loc, aer_cg_sw_loc, aer_piz_sw_loc,  & 
         &  rnseeds,                                                        &
         &  laytrop,        jp,             iabs,          gases,           &
         &  colmol,         fac,            h2o_factor,    h2o_fraction,    &
         &  h2o_index,                                                      &
         &  flx_dnsw,       flx_upsw,       flx_dnsw_clr,  flx_upsw_clr,    &
         &  bnd_wght,       per_band_flux)

    CALL srtm_diags(kproma, kbdim, per_band_flux,     &
      vis_dn_dir_sfc, par_dn_dir_sfc, nir_dn_dir_sfc, &
      vis_dn_dff_sfc, par_dn_dff_sfc, nir_dn_dff_sfc, &
      vis_up_sfc,     par_up_sfc,     nir_up_sfc      )

    IF (ltimer) CALL timer_stop(timer_srtm)

  END SUBROUTINE psrad_interface_onBlock


  SUBROUTINE flip_ud(n, v, u)
    INTEGER, INTENT(IN) :: n
    REAL(wp), INTENT(IN) :: v(:,:)
    REAL(wp), INTENT(INOUT) :: u(:,:)
    INTEGER :: m
    m = SIZE(v,2)
    u(1:n,1:m) = v(1:n,m:1:-1)
  END SUBROUTINE flip_ud


  SUBROUTINE shift_and_call_psrad_interface_onBlock(     &
       & jg,             jb,                             &
       & jcs,            jce,                            &
       & kbdim,          klev,           klevp1,         &
       !
       & iaero,          ktype,                          &
       & psctm,          ssi_factor,                     &
       & laland,         laglac,         this_datetime,  &
       & pcos_mu0,       daylght_frc,                    &
       & alb_vis_dir,    alb_nir_dir,                    &
       & alb_vis_dif,    alb_nir_dif,                    &
       & zf,             zh,             dz,             &
       & pp_sfc,         pp_fl,                          &
       & tk_sfc,         tk_fl,          tk_hl,          &
       & xm_dry,         xm_vap,         xm_liq,         &
       & xm_ice,         cdnc,           xc_frc,         &
       & xm_co2,         xm_ch4,         xm_n2o,         &
       & xm_cfc,         xm_o3,          xm_o2,          &
       !
       & lw_upw,         lw_upw_clr,                     &
       & lw_dnw,         lw_dnw_clr,                     &
       & sw_upw,         sw_upw_clr,                     &
       & sw_dnw,         sw_dnw_clr,                     &
       & vis_dn_dir_sfc, par_dn_dir_sfc, nir_dn_dir_sfc, &
       & vis_dn_dff_sfc, par_dn_dff_sfc, nir_dn_dff_sfc, &
       & vis_up_sfc,     par_up_sfc,     nir_up_sfc      )

    INTEGER,INTENT(IN)  :: &
         & jg,             & !< domain index
         & jb,               & !< first dimension of 2-d arrays
         & jcs,              & !< cell index, start
         & jce,              & !< cell index, end
         & kbdim,            & !< first dimension of 2-d arrays
         & klev,             & !< number of full levels
         & klevp1,           & !< number of half levels
         & iaero,            & !< aerosol control
         & ktype(:)            !< type of convection

    REAL(wp),INTENT(IN) :: psctm                         !< orbit and time dependent solar constant for radiation time step
    REAL(wp),INTENT(IN) :: ssi_factor(nbndsw)            !< fraction of TSI in the 14 RRTM SW bands

    LOGICAL,INTENT(IN) :: &
         & laland(:),   & !< land sea mask, land=.true.
         & laglac(:)      !< glacier mask, glacier=.true.

    TYPE(datetime), POINTER ::  this_datetime !< actual time step

    REAL(WP),INTENT(IN)  ::    &
         & pcos_mu0(:),      & !< mu0 for solar zenith angle
         & daylght_frc(:),   & !< daylight fraction; with diurnal cycle 0 or 1, with zonal mean in [0,1]
         & alb_vis_dir(:),   & !< surface albedo for vis range and dir light
         & alb_nir_dir(:),   & !< surface albedo for NIR range and dir light
         & alb_vis_dif(:),   & !< surface albedo for vis range and dif light
         & alb_nir_dif(:),   & !< surface albedo for NIR range and dif light
         & zf(:,:),       & !< geometric height at full level in m
         & zh(:,:),     & !< geometric height at half level in m
         & dz(:,:),       & !< geometric height thickness in m
         & pp_sfc(:),        & !< surface pressure in Pa
         & pp_fl(:,:),    & !< full level pressure in Pa
         & tk_sfc(:),        & !< surface temperature in K
         & tk_fl(:,:),    & !< full level temperature in K
         & tk_hl(:,:),  & !< half level temperature in K
         & xm_dry(:,:),   & !< dry air     mass in kg/m2
         & xm_vap(:,:),   & !< water vapor mass in kg/m2
         & xm_liq(:,:),   & !< cloud water mass in kg/m2
         & xm_ice(:,:),   & !< cloud ice   mass in kg/m2
         & cdnc(:,:),     & !< cloud nuclei concentration
         & xc_frc(:,:),   & !< fractional cloud cover
         & xm_co2(:,:),   & !< co2 mass in kg/m2
         & xm_ch4(:,:),   & !< ch4 mass in kg/m2
         & xm_n2o(:,:),   & !< n2o mass in kg/m2
         & xm_cfc(kbdim,klev,2), & !< cfc mass in kg/m2
         & xm_o3(:,:),    & !< o3  mass in kg/m2
         & xm_o2(:,:)       !< o2  mass in kg/m2

    REAL (wp), INTENT (OUT) ::       &
         & lw_upw    (:,:), & !<   upward LW flux profile, all sky
         & lw_upw_clr(:,:), & !<   upward LW flux profile, clear sky
         & lw_dnw    (:,:), & !< downward LW flux profile, all sky
         & lw_dnw_clr(:,:), & !< downward LW flux profile, clear sky
         & sw_upw    (:,:), & !<   upward SW flux profile, all sky
         & sw_upw_clr(:,:), & !<   upward SW flux profile, clear sky
         & sw_dnw    (:,:), & !< downward SW flux profile, all sky
         & sw_dnw_clr(:,:)    !< downward SW flux profile, clear sky

    REAL (wp), INTENT (OUT) :: &
         & vis_dn_dir_sfc(:) , & !< Diffuse downward flux surface visible radiation 
         & par_dn_dir_sfc(:) , & !< Diffuse downward flux surface PAR
         & nir_dn_dir_sfc(:) , & !< Diffuse downward flux surface near-infrared radiation
         & vis_dn_dff_sfc(:) , & !< Direct  downward flux surface visible radiation 
         & par_dn_dff_sfc(:) , & !< Direct  downward flux surface PAR
         & nir_dn_dff_sfc(:) , & !< Direct  downward flux surface near-infrared radiation
         & vis_up_sfc    (:) , & !< Upward  flux surface visible radiation 
         & par_up_sfc    (:) , & !< Upward  flux surface PAR
         & nir_up_sfc    (:)     !< Upward  flux surface near-infrared radiation



    INTEGER ::                 &
         & s_jce               !< cell index, shifted end

    ! Shifted input arguments
    !
    INTEGER ::                 &
         & s_ktype          (kbdim)            !< type of convection
    !
    LOGICAL ::                 &
         & s_laland         (kbdim),   & !< land sea mask, land=.true.
         & s_laglac         (kbdim)      !< glacier mask, glacier=.true.
    !
    REAL(wp)  ::                    &
         & s_pcos_mu0       (kbdim),      & !< mu0 for solar zenith angle
         & s_daylght_frc    (kbdim),   & !< daylight fraction; with diurnal cycle 0 or 1, with zonal mean in [0,1]
         & s_alb_vis_dir    (kbdim),   & !< surface albedo for vis range and dir light
         & s_alb_nir_dir    (kbdim),   & !< surface albedo for NIR range and dir light
         & s_alb_vis_dif    (kbdim),   & !< surface albedo for vis range and dif light
         & s_alb_nir_dif    (kbdim),   & !< surface albedo for NIR range and dif light
         & s_zf             (kbdim,klev),       & !< geometric height at full level in m
         & s_zh             (kbdim,klevp1),     & !< geometric height at half level in m
         & s_dz             (kbdim,klev),       & !< geometric height thickness in m
         & s_pp_sfc         (kbdim),        & !< surface pressure in Pa
         & s_pp_fl          (kbdim,klev),    & !< full level pressure in Pa
         & s_tk_sfc         (kbdim),        & !< surface temperature in K
         & s_tk_fl          (kbdim,klev),    & !< full level temperature in K
         & s_tk_hl          (kbdim,klevp1),  & !< half level temperature in K
         & s_xm_dry         (kbdim,klev),   & !< dry air     mass in kg/m2
         & s_xm_vap         (kbdim,klev),   & !< water vapor mass in kg/m2
         & s_xm_liq         (kbdim,klev),   & !< cloud water mass in kg/m2
         & s_xm_ice         (kbdim,klev),   & !< cloud ice   mass in kg/m2
         & s_cdnc           (kbdim,klev),     & !< cloud nuclei concentration
         & s_xc_frc         (kbdim,klev),   & !< fractional cloud cover
         & s_xm_co2         (kbdim,klev),   & !< co2 mass in kg/m2
         & s_xm_ch4         (kbdim,klev),   & !< ch4 mass in kg/m2
         & s_xm_n2o         (kbdim,klev),   & !< n2o mass in kg/m2
         & s_xm_cfc         (kbdim,klev,2),  & !< cfc mass in kg/m2
         & s_xm_o3          (kbdim,klev),    & !< o3  mass in kg/m2
         & s_xm_o2          (kbdim,klev)       !< o2  mass in kg/m2

    ! Shifted output arguments
    !
    REAL (wp) ::                          &
         & s_lw_upw       (kbdim,klevp1), & !<   upward LW flux profile, all sky
         & s_lw_upw_clr   (kbdim,klevp1), & !<   upward LW flux profile, clear sky
         & s_lw_dnw       (kbdim,klevp1), & !< downward LW flux profile, all sky
         & s_lw_dnw_clr   (kbdim,klevp1), & !< downward LW flux profile, clear sky
         & s_sw_upw       (kbdim,klevp1), & !<   upward SW flux profile, all sky
         & s_sw_upw_clr   (kbdim,klevp1), & !<   upward SW flux profile, clear sky
         & s_sw_dnw       (kbdim,klevp1), & !< downward SW flux profile, all sky
         & s_sw_dnw_clr   (kbdim,klevp1)    !< downward SW flux profile, clear sky
    !
    REAL (wp) ::                     &
         & s_vis_dn_dir_sfc (kbdim), & !< Diffuse downward flux surface visible radiation 
         & s_par_dn_dir_sfc (kbdim), & !< Diffuse downward flux surface PAR
         & s_nir_dn_dir_sfc (kbdim), & !< Diffuse downward flux surface near-infrared radiation
         & s_vis_dn_dff_sfc (kbdim), & !< Direct  downward flux surface visible radiation 
         & s_par_dn_dff_sfc (kbdim), & !< Direct  downward flux surface PAR
         & s_nir_dn_dff_sfc (kbdim), & !< Direct  downward flux surface near-infrared radiation
         & s_vis_up_sfc     (kbdim), & !< Upward  flux surface visible radiation 
         & s_par_up_sfc     (kbdim), & !< Upward  flux surface PAR
         & s_nir_up_sfc     (kbdim)    !< Upward  flux surface near-infrared radiation

    ! Shift input arguments
    !
    s_jce = jce-jcs+1
    !
    s_ktype       (1:s_jce)     = ktype       (jcs:jce)
    s_laland      (1:s_jce)     = laland      (jcs:jce)
    s_laglac      (1:s_jce)     = laglac      (jcs:jce)
    s_pcos_mu0    (1:s_jce)     = pcos_mu0    (jcs:jce)
    s_daylght_frc (1:s_jce)     = daylght_frc (jcs:jce)
    s_alb_vis_dir (1:s_jce)     = alb_vis_dir (jcs:jce)
    s_alb_nir_dir (1:s_jce)     = alb_nir_dir (jcs:jce)
    s_alb_vis_dif (1:s_jce)     = alb_vis_dif (jcs:jce)
    s_alb_nir_dif (1:s_jce)     = alb_nir_dif (jcs:jce)
    s_zf          (1:s_jce,:)   = zf          (jcs:jce,:)
    s_zh          (1:s_jce,:)   = zh          (jcs:jce,:)
    s_dz          (1:s_jce,:)   = dz          (jcs:jce,:)
    s_pp_sfc      (1:s_jce)     = pp_sfc      (jcs:jce)
    s_pp_fl       (1:s_jce,:)   = pp_fl       (jcs:jce,:)
    s_tk_sfc      (1:s_jce)     = tk_sfc      (jcs:jce)
    s_tk_fl       (1:s_jce,:)   = tk_fl       (jcs:jce,:)
    s_tk_hl       (1:s_jce,:)   = tk_hl       (jcs:jce,:)
    s_xm_dry      (1:s_jce,:)   = xm_dry      (jcs:jce,:)
    s_xm_vap      (1:s_jce,:)   = xm_vap      (jcs:jce,:)
    s_xm_liq      (1:s_jce,:)   = xm_liq      (jcs:jce,:)
    s_xm_ice      (1:s_jce,:)   = xm_ice      (jcs:jce,:)
    s_cdnc        (1:s_jce,:)   = cdnc        (jcs:jce,:)
    s_xc_frc      (1:s_jce,:)   = xc_frc      (jcs:jce,:)
    s_xm_co2      (1:s_jce,:)   = xm_co2      (jcs:jce,:)
    s_xm_ch4      (1:s_jce,:)   = xm_ch4      (jcs:jce,:)
    s_xm_n2o      (1:s_jce,:)   = xm_n2o      (jcs:jce,:)
    s_xm_cfc      (1:s_jce,:,:) = xm_cfc      (jcs:jce,:,:)
    s_xm_o3       (1:s_jce,:)   = xm_o3       (jcs:jce,:)
    s_xm_o2       (1:s_jce,:)   = xm_o2       (jcs:jce,:)

    ! Call radiation with shifted input arguments and receive shifted output arguments
    !
    CALL psrad_interface_onBlock(                                         &
         &   jg,                  jb,                                     &
         &                      s_jce,                                    &
         &   kbdim,               klev,                                   &
         !
         &   iaero,             s_ktype(:),                               &
         &   psctm,               ssi_factor,                             &
         & s_laland(:),         s_laglac(:),         this_datetime,       &
         & s_pcos_mu0(:),       s_daylght_frc(:),                         &
         & s_alb_vis_dir(:),    s_alb_nir_dir(:),                         &
         & s_alb_vis_dif(:),    s_alb_nir_dif(:),                         &
         & s_zf(:,:),           s_zh(:,:),           s_dz(:,:),           &
         & s_pp_sfc(:),         s_pp_fl(:,:),                             &
         & s_tk_sfc(:),         s_tk_fl(:,:),        s_tk_hl(:,:),        &
         & s_xm_dry(:,:),       s_xm_vap(:,:),       s_xm_liq(:,:),       &
         & s_xm_ice(:,:),       s_cdnc(:,:),         s_xc_frc(:,:),       &
         & s_xm_co2(:,:),       s_xm_ch4(:,:),       s_xm_n2o (:,:),      &
         & s_xm_cfc(:,:,:),     s_xm_o3(:,:),        s_xm_o2(:,:),        &
         !
         & s_lw_upw(:,:),       s_lw_upw_clr (:,:),                       &
         & s_lw_dnw(:,:),       s_lw_dnw_clr (:,:),                       &
         & s_sw_upw(:,:),       s_sw_upw_clr(:,:),                        &
         & s_sw_dnw(:,:),       s_sw_dnw_clr(:,:),                        &
         & s_vis_dn_dir_sfc(:), s_par_dn_dir_sfc(:), s_nir_dn_dir_sfc(:), &
         & s_vis_dn_dff_sfc(:), s_par_dn_dff_sfc(:), s_nir_dn_dff_sfc(:), &
         & s_vis_up_sfc(:),     s_par_up_sfc(:),     s_nir_up_sfc(:)      )

    ! Shift output arguments
    !
    lw_upw         (jcs:jce,:) = s_lw_upw         (1:s_jce,:)
    lw_upw_clr     (jcs:jce,:) = s_lw_upw_clr     (1:s_jce,:)
    lw_dnw         (jcs:jce,:) = s_lw_dnw         (1:s_jce,:)
    lw_dnw_clr     (jcs:jce,:) = s_lw_dnw_clr     (1:s_jce,:)
    sw_upw         (jcs:jce,:) = s_sw_upw         (1:s_jce,:)
    sw_upw_clr     (jcs:jce,:) = s_sw_upw_clr     (1:s_jce,:)
    sw_dnw         (jcs:jce,:) = s_sw_dnw         (1:s_jce,:)
    sw_dnw_clr     (jcs:jce,:) = s_sw_dnw_clr     (1:s_jce,:)
    vis_dn_dir_sfc (jcs:jce)   = s_vis_dn_dir_sfc (1:s_jce)
    par_dn_dir_sfc (jcs:jce)   = s_par_dn_dir_sfc (1:s_jce)
    nir_dn_dir_sfc (jcs:jce)   = s_nir_dn_dir_sfc (1:s_jce)
    vis_dn_dff_sfc (jcs:jce)   = s_vis_dn_dff_sfc (1:s_jce)
    par_dn_dff_sfc (jcs:jce)   = s_par_dn_dff_sfc (1:s_jce)
    nir_dn_dff_sfc (jcs:jce)   = s_nir_dn_dff_sfc (1:s_jce)
    vis_up_sfc     (jcs:jce)   = s_vis_up_sfc     (1:s_jce)
    par_up_sfc     (jcs:jce)   = s_par_up_sfc     (1:s_jce)
    nir_up_sfc     (jcs:jce)   = s_nir_up_sfc     (1:s_jce)


  END SUBROUTINE shift_and_call_psrad_interface_onBlock
      
END MODULE mo_psrad_interface
