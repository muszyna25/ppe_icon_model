!>
!! This module is the interface between ICON:nwp_radiation to the radiation scheme ecRad
!!
!! - There are two interfaces within this module: nwp_ecrad_radiation and
!!   nwp_ecrad_radiation_reduced. The former provides the interface to ecrad on the full
!!   ICON dynamics grid. The latter one provides the interface on a grid with a reduced
!!   spatial resolution.
!! - The decision which of the two interface routines is used, is done via the namelist
!!   switch lredgrid_phys. Based on the value of lredgrid_phys, the correct interface
!!   is called by mo_nwp_rad_interface:nwp_radiation.
!! - The interfaces have to fill the different ecRad input types (ecrad_aerosol, 
!!   ecrad_single_level, ecrad_thermodynamics, ecrad_gas, ecrad_cloud) with data from 
!!   ICON variables. Then, the ecRad radiation code is called. At the end, the fluxes
!!   calculated by ecRad and stored in the ecrad_flux structure are copied to ICON variables.
!! - The difference between nwp_ecrad_radiation and nwp_ecrad_radiation_reduced is mostly
!!   an upscaling at the beginning and a downscaling at the end of the interface.
!! - The transfer of data from ICON to ecRad and vice versa is performed within
!!   routines from mo_nwp_ecrad_utilities and mo_nwp_ecrad_prep_aerosol, independent of
!!   the choice to use a reduced radiation grid or not.
!!
!! @author Daniel Rieger, Deutscher Wetterdienst, Offenbach
!!
!! @par Revision History
!! Initial release by Daniel Rieger, Deutscher Wetterdienst, Offenbach (2019-01-31)
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

MODULE mo_nwp_ecrad_interface

  USE mo_kind,                   ONLY: wp
  USE mo_exception,              ONLY: finish, message
  USE mo_math_constants,         ONLY: pi
  USE mo_model_domain,           ONLY: t_patch, p_patch_local_parent
  USE mo_impl_constants,         ONLY: min_rlcell_int, MAX_CHAR_LENGTH, nexlevs_rrg_vnest
  USE mo_impl_constants_grf,     ONLY: grf_bdywidth_c, grf_ovlparea_start_c, grf_fbk_start_c
  USE mo_fortran_tools,          ONLY: init
  USE mo_parallel_config,        ONLY: nproma
  USE mo_loopindices,            ONLY: get_indices_c
  USE mo_grid_config,            ONLY: l_limited_area
  USE mo_ext_data_types,         ONLY: t_external_data
  USE mo_nwp_lnd_types,          ONLY: t_lnd_prog
  USE mo_nonhydro_types,         ONLY: t_nh_diag
  USE mo_nwp_phy_types,          ONLY: t_nwp_phy_diag
  USE mo_physical_constants,     ONLY: rhoh2o
  USE mo_run_config,             ONLY: msg_level, iqv, iqi, iqc
  USE mo_atm_phy_nwp_config,     ONLY: atm_phy_nwp_config
  USE mo_radiation_config,       ONLY: irad_aero
  USE mo_phys_nest_utilities,    ONLY: upscale_rad_input, downscale_rad_output
  USE mtime,                     ONLY: datetime
#ifdef __ECRAD
  USE mo_ecrad,                  ONLY: ecrad,                                    &
                                   &   t_ecrad_conf, t_ecrad_aerosol_type,       &
                                   &   t_ecrad_single_level_type,                &
                                   &   t_ecrad_thermodynamics_type,              &
                                   &   t_ecrad_gas_type, t_ecrad_flux_type,      &
                                   &   t_ecrad_cloud_type
  USE mo_nwp_ecrad_prep_aerosol, ONLY: nwp_ecrad_prep_aerosol
  USE mo_nwp_ecrad_utilities,    ONLY: ecrad_set_single_level,                   &
                                   &   ecrad_set_thermodynamics,                 &
                                   &   ecrad_set_clouds,                         &
                                   &   ecrad_set_gas,                            &
                                   &   ecrad_store_fluxes, add_3D_diffuse_rad
#endif


  IMPLICIT NONE

  PRIVATE
#ifdef __ECRAD
  !> module name string
  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_nwp_ecrad_interface'


  PUBLIC :: nwp_ecrad_radiation
  PUBLIC :: nwp_ecrad_radiation_reduced


CONTAINS


  !---------------------------------------------------------------------------------------
  !>
  !! SUBROUTINE nwp_ecrad_radiation:
  !! Interface to ecRad on full grid. This routine
  !!  ... allocates the ecRad data types
  !!  ... fills the ecRad data types with current atmospheric and external data
  !!  ... saves the output to ICON physics data structure
  !!
  !! @par Revision History
  !! Initial release by Daniel Rieger, Deutscher Wetterdienst, Offenbach (2019-01-31)
  !!
  SUBROUTINE nwp_ecrad_radiation ( current_datetime, pt_patch, ext_data,                    &
    &  zaeq1, zaeq2, zaeq3, zaeq4, zaeq5, pt_diag, prm_diag, lnd_prog, ecrad_conf )

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER::  &
      &  routine = modname//'::nwp_ecrad_radiation'

    TYPE(datetime),          INTENT(in)    :: current_datetime !< Current date and time

    TYPE(t_patch), TARGET,   INTENT(in)    :: pt_patch         !< Current domain info
    TYPE(t_external_data),   INTENT(in)    :: ext_data         !< External data container

    REAL(wp),                INTENT(in)    ::             &
      & zaeq1(nproma,pt_patch%nlev,pt_patch%nblks_c),     & !< Climatological aerosol (Tegen)
      & zaeq2(nproma,pt_patch%nlev,pt_patch%nblks_c),     & !< Climatological aerosol (Tegen)
      & zaeq3(nproma,pt_patch%nlev,pt_patch%nblks_c),     & !< Climatological aerosol (Tegen)
      & zaeq4(nproma,pt_patch%nlev,pt_patch%nblks_c),     & !< Climatological aerosol (Tegen)
      & zaeq5(nproma,pt_patch%nlev,pt_patch%nblks_c)        !< Climatological aerosol (Tegen)

    TYPE(t_nh_diag), TARGET, INTENT(in)    :: pt_diag       !< ICON diagnostic variables
    TYPE(t_nwp_phy_diag),    INTENT(inout) :: prm_diag      !< ICON physics diagnostics
    TYPE(t_lnd_prog),        INTENT(inout) :: lnd_prog      !< ICON prognostic land state

    TYPE(t_ecrad_conf),      INTENT(in)    :: ecrad_conf    !< ecRad configuration object
! Local variables
    TYPE(t_ecrad_aerosol_type)        :: &
      &  ecrad_aerosol                     !< ecRad aerosol information (input)
    TYPE(t_ecrad_single_level_type)   :: &
      &  ecrad_single_level                !< ecRad single level information (input)
    TYPE(t_ecrad_thermodynamics_type) :: &
      &  ecrad_thermodynamics              !< ecRad thermodynamics information (input)
    TYPE(t_ecrad_gas_type)            :: &
      &  ecrad_gas                         !< ecRad gas information (input)
    TYPE(t_ecrad_cloud_type)          :: &
      &  ecrad_cloud                       !< ecRad cloud information (input)
    TYPE(t_ecrad_flux_type)           :: &
      &  ecrad_flux                        !< ecRad flux information (output)
    REAL(wp)                 :: &
      &  fact_reffc               !< Factor in the calculation of cloud droplet effective radius
    INTEGER                  :: &
      &  jc, jb,                & !< Loop indices
      &  jg,                    & !< Domain index
      &  nlev, nlevp1,          & !< Number of vertical levels (full, half)
      &  rl_start, rl_end,      & !< 
      &  i_startblk, i_endblk,  & !< blocks
      &  i_startidx, i_endidx     !< slices
    LOGICAL, ALLOCATABLE     :: &
      &  cosmu0mask(:)            !< Mask if cosmu0 > 0

    nlev      = pt_patch%nlev
    nlevp1    = nlev+1
    jg        = pt_patch%id

    fact_reffc = (3.0e-9_wp/(4.0_wp*pi*rhoh2o))**(1.0_wp/3.0_wp)
    
    ALLOCATE(cosmu0mask(nproma))

    IF (msg_level >= 7) &
      &       CALL message(routine, 'ecrad radiation on full grid')

    CALL ecrad_single_level%allocate(nproma, 2, 1, .true.) !< use_sw_albedo_direct, 2 bands
    ecrad_single_level%solar_irradiance = 1._wp            !< Obtain normalized fluxes which corresponds to the 
                                                           !< transmissivity needed in the following

    IF (ecrad_conf%use_spectral_solar_scaling) THEN
      ALLOCATE(ecrad_single_level%spectral_solar_scaling(ecrad_conf%n_bands_sw))
      ecrad_single_level%spectral_solar_scaling = (/  1.0_wp, 1.0_wp, 1.0_wp, 1.0478_wp, 1.0404_wp, 1.0317_wp, &
         &   1.0231_wp, 1.0054_wp, 0.98413_wp, 0.99863_wp, 0.99907_wp, 0.90589_wp, 0.92213_wp, 1.0_wp /)
    ENDIF

    CALL ecrad_thermodynamics%allocate(nproma, nlev, use_h2o_sat=.false., rrtm_pass_temppres_fl=.true.)

    CALL ecrad_gas%allocate(nproma, nlev)
    
    CALL ecrad_cloud%allocate(nproma, nlev)
    ! Currently hardcoded values for FSD
    CALL ecrad_cloud%create_fractional_std(nproma, nlev, 1._wp)

    IF ( ecrad_conf%use_aerosols ) THEN
      ! Allocate aerosol container
      CALL ecrad_aerosol%allocate_direct(ecrad_conf, nproma, 1, nlev)
    ENDIF

    CALL ecrad_flux%allocate(ecrad_conf, 1, nproma, nlev)

    rl_start = grf_bdywidth_c+1
    rl_end   = min_rlcell_int

    i_startblk = pt_patch%cells%start_block(rl_start)
    i_endblk   = pt_patch%cells%end_block(rl_end)

!$OMP PARALLEL PRIVATE(jb,jc,i_startidx,i_endidx,cosmu0mask)                         &
!$OMP          FIRSTPRIVATE(ecrad_aerosol,ecrad_single_level, ecrad_thermodynamics,  &
!$OMP                       ecrad_gas, ecrad_cloud,ecrad_flux)
!$OMP DO ICON_OMP_GUIDED_SCHEDULE
    DO jb = i_startblk, i_endblk
      CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
        &                i_startidx, i_endidx, rl_start, rl_end)

      ! It may happen that an MPI patch contains only nest boundary points
      ! In this case, no action is needed
      IF (i_startidx > i_endidx) CYCLE

      prm_diag%tsfctrad(i_startidx:i_endidx,jb) = lnd_prog%t_g(i_startidx:i_endidx,jb)

      cosmu0mask(:) = .FALSE.
      DO jc = i_startidx, i_endidx
        IF ( prm_diag%cosmu0(jc,jb) > 0._wp ) THEN
          cosmu0mask(jc) = .TRUE.
        ENDIF
      ENDDO

! Fill single level configuration type
      CALL ecrad_set_single_level(ecrad_single_level, current_datetime, pt_patch%cells%center(:,jb),           &
        &                         prm_diag%cosmu0(:,jb), prm_diag%tsfctrad(:,jb), prm_diag%albvisdif(:,jb),    &
        &                         prm_diag%albnirdif(:,jb), prm_diag%albvisdir(:,jb),                          &
        &                         prm_diag%albnirdir(:,jb), prm_diag%lw_emiss(:,jb),                           &
        &                         i_startidx, i_endidx)
      
! Fill thermodynamics configuration type
      CALL ecrad_set_thermodynamics(ecrad_thermodynamics, pt_diag%temp(:,:,jb), pt_diag%pres(:,:,jb),          &
        &                           pt_diag%pres_ifc(:,:,jb), prm_diag%tsfctrad(:,jb),                         &
        &                           nlev, nlevp1, i_startidx, i_endidx)

! Fill gas configuration type
      CALL ecrad_set_gas(ecrad_gas, ecrad_conf, ext_data%atm%o3(:,:,jb), prm_diag%tot_cld(:,:,jb,iqv), &
        &                pt_diag%pres(:,:,jb), i_startidx, i_endidx, nlev)

! Fill clouds configuration type
      CALL ecrad_set_clouds(ecrad_cloud, ecrad_thermodynamics, prm_diag%tot_cld(:,:,jb,iqc),       &
        &                   prm_diag%tot_cld(:,:,jb,iqi), prm_diag%clc(:,:,jb),                    &
        &                   pt_diag%temp(:,:,jb), pt_diag%pres(:,:,jb), prm_diag%acdnc(:,:,jb),    &
        &                   ext_data%atm%fr_glac_smt(:,jb), ext_data%atm%fr_land_smt(:,jb),        &
        &                   fact_reffc, ecrad_conf%cloud_fraction_threshold, nlev, i_startidx, i_endidx)

! Fill aerosol configuration type
      SELECT CASE (irad_aero)
        CASE(0)
          ! No aerosol, nothing to do
        CASE(2)
          ! Case 2: Constant aerosol
          !         Arguments can be added to fill ecrad_aerosol with actual values. For the time being,
          !         we stay consistent with RRTM where irad_aero=2 does not add any aerosol
          CALL nwp_ecrad_prep_aerosol(ecrad_conf, ecrad_aerosol)
        CASE(5,6)
          ! Fill aerosol configuration type with Tanre or Tegen aerosol
          CALL nwp_ecrad_prep_aerosol(1, nlev, i_startidx, i_endidx,     &
            &                         zaeq1(:,:,jb), zaeq2(:,:,jb),      &
            &                         zaeq3(:,:,jb), zaeq4(:,:,jb),      &
            &                         zaeq5(:,:,jb),                     &
            &                         ecrad_conf, ecrad_aerosol)
        CASE(9)
          ! Use ART aerosol
          CALL nwp_ecrad_prep_aerosol(1, nlev, i_startidx, i_endidx, jb, pt_patch%id, &
            &                         zaeq1(:,:,jb), zaeq2(:,:,jb),                   &
            &                         zaeq3(:,:,jb), zaeq4(:,:,jb),                   &
            &                         zaeq5(:,:,jb),                                  &
            &                         ecrad_conf, ecrad_aerosol)
          CALL finish(TRIM(routine),'irad_aero = 9 not yet fully implemented for ecRad')
        CASE DEFAULT
          CALL finish(TRIM(routine),'irad_aero not valid for ecRad')
      END SELECT

      ecrad_flux%cloud_cover_sw(:) = 0._wp
      ecrad_flux%cloud_cover_lw(:) = 0._wp

!---------------------------------------------------------------------------------------
! Call the radiation scheme ecRad
!---------------------------------------------------------------------------------------
      CALL ecrad(nproma, nlev, i_startidx, i_endidx,      & !< Array and loop bounds (input)
        &        ecrad_conf,                              & !< General ecRad configuration object (input)
        &        ecrad_single_level,                      & !< ecRad single level configuration object (input)
        &        ecrad_thermodynamics,                    & !< ecRad thermodynamics configuration object (input)
        &        ecrad_gas,                               & !< ecRad gas configuration object (input)
        &        ecrad_cloud,                             & !< ecRad cloud configuration object (input)
        &        ecrad_aerosol,                           & !< ecRad aerosol configuration object (input)
        &        ecrad_flux                               ) !< ecRad fluxes in the longwave BUT flux/solar constant in the shortwave (output)

!---------------------------------------------------------------------------------------

! Update ICON variables with fluxes from ecRad
      CALL ecrad_store_fluxes(jg, ecrad_flux, prm_diag%cosmu0(:,jb), prm_diag%trsolall    (:,:,jb),  &
        &                     prm_diag%trsol_up_toa   (:,jb), prm_diag%trsol_up_sfc     (:,jb),  &
        &                     prm_diag%trsol_par_sfc  (:,jb), prm_diag%trsol_dn_sfc_diff(:,jb),  &
        &                     prm_diag%trsolclr_sfc   (:,jb), prm_diag%lwflxall       (:,:,jb),  &
        &                     prm_diag%lwflx_up_sfc_rs(:,jb), prm_diag%lwflxclr_sfc     (:,jb),  &
        &                     prm_diag%lwflx_up     (:,:,jb), prm_diag%lwflx_dn       (:,:,jb),  &
        &                     prm_diag%swflx_up     (:,:,jb), prm_diag%swflx_dn       (:,:,jb),  &
        &                     prm_diag%lwflx_up_clr (:,:,jb), prm_diag%lwflx_dn_clr   (:,:,jb),  &
        &                     prm_diag%swflx_up_clr (:,:,jb), prm_diag%swflx_dn_clr   (:,:,jb),  &  
        &                     cosmu0mask(:), i_startidx, i_endidx, nlevp1)

      ! Add 3D contribution to diffuse radiation
      CALL add_3D_diffuse_rad(ecrad_flux, prm_diag%clc(:,:,jb), pt_diag%pres(:,:,jb), pt_diag%temp(:,:,jb),      &
        &                     prm_diag%cosmu0(:,jb), prm_diag%trsol_dn_sfc_diff(:,jb), i_startidx, i_endidx, nlev)

    ENDDO ! jb
!$OMP END DO NOWAIT
!$OMP END PARALLEL

! CLEANUP
    CALL ecrad_single_level%deallocate
    CALL ecrad_thermodynamics%deallocate
    CALL ecrad_gas%deallocate
    CALL ecrad_cloud%deallocate
    IF ( ecrad_conf%use_aerosols ) CALL ecrad_aerosol%deallocate
    CALL ecrad_flux%deallocate
    DEALLOCATE(cosmu0mask)
  

  END SUBROUTINE nwp_ecrad_radiation
  !---------------------------------------------------------------------------------------


  !---------------------------------------------------------------------------------------
  !>
  !! SUBROUTINE nwp_ecrad_radiation_reduced:
  !! Interface to ecRad on reduced radiation grid. This routine
  !!  ... allocates the ecRad data types
  !!  ... fills the ecRad data types with current atmospheric and external data
  !!  ... saves the output to ICON physics data structure
  !!
  !! @par Revision History
  !! Initial release by Daniel Rieger, Deutscher Wetterdienst, Offenbach (2019-01-31)
  !! Open TODOs: dust_tunefac not considered so far
  !!
  SUBROUTINE nwp_ecrad_radiation_reduced (current_datetime, pt_patch, pt_par_patch, ext_data, &
    &                                     zaeq1,zaeq2,zaeq3,zaeq4,zaeq5,                      &
    &                                     pt_diag,prm_diag,lnd_prog, ecrad_conf )

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER::  &
      &  routine = modname//'::nwp_ecrad_radiation_reduced'

    TYPE(datetime),          INTENT(in)    :: current_datetime !< Current date and time

    TYPE(t_patch), TARGET,   INTENT(in)    :: pt_patch         !< Current domain info
    TYPE(t_patch), TARGET,   INTENT(in)    :: pt_par_patch     !< Parent domain info
    TYPE(t_external_data),   INTENT(in)    :: ext_data         !< External data container

    REAL(wp),                INTENT(in)    ::             &
      & zaeq1(nproma,pt_patch%nlev,pt_patch%nblks_c),     & !< Climatological aerosol (Tegen)
      & zaeq2(nproma,pt_patch%nlev,pt_patch%nblks_c),     & !< Climatological aerosol (Tegen)
      & zaeq3(nproma,pt_patch%nlev,pt_patch%nblks_c),     & !< Climatological aerosol (Tegen)
      & zaeq4(nproma,pt_patch%nlev,pt_patch%nblks_c),     & !< Climatological aerosol (Tegen)
      & zaeq5(nproma,pt_patch%nlev,pt_patch%nblks_c)        !< Climatological aerosol (Tegen)

    TYPE(t_nh_diag), TARGET, INTENT(in)    :: pt_diag       !< ICON diagnostic variables
    TYPE(t_nwp_phy_diag),    INTENT(inout) :: prm_diag      !< ICON physics diagnostics
    TYPE(t_lnd_prog),        INTENT(inout) :: lnd_prog      !< ICON prognostic land state

    TYPE(t_ecrad_conf),      INTENT(in)    :: ecrad_conf    !< ecRad configuration object
! Local variables
    TYPE(t_patch), POINTER            :: &
      &  ptr_pp                            !< Pointer to parent patch of current domain
    TYPE(t_ecrad_aerosol_type)        :: &
      &  ecrad_aerosol                     !< ecRad aerosol information (input)
    TYPE(t_ecrad_single_level_type)   :: &
      &  ecrad_single_level                !< ecRad single level information (input)
    TYPE(t_ecrad_thermodynamics_type) :: &
      &  ecrad_thermodynamics              !< ecRad thermodynamics information (input)
    TYPE(t_ecrad_gas_type)            :: &
      &  ecrad_gas                         !< ecRad gas information (input)
    TYPE(t_ecrad_cloud_type)          :: &
      &  ecrad_cloud                       !< ecRad cloud information (input)
    TYPE(t_ecrad_flux_type)           :: &
      &  ecrad_flux                        !< ecRad flux information (output)
    REAL(wp)                 :: &
      &  fact_reffc               !< Factor in the calculation of cloud droplet effective radius
    INTEGER                  :: &
      &  nblks_par_c,           & !< nblks for reduced grid (parent domain)
      &  nblks_lp_c,            & !< nblks for reduced grid (local parent)
      &  jb, jc, jk,            & !< loop indices
      &  jg,                    & !< domain id
      &  nlev,                  & !< number of full levels
      &  nlev_rg, nlev_rgp1,    & !< number of full and half levels at reduced grid
      &  rl_start, rl_end,      & !< 
      &  i_startblk, i_endblk,  & !< blocks
      &  i_startidx, i_endidx,  & !< slices
      &  np, nl                   !< dimension variables for allocation (3d fluxes)
    ! For radiation on reduced grid
    ! These fields need to be allocatable because they have different dimensions for
    ! the global grid and nested grids, and for runs with/without MPI parallelization
    ! Input fields
    REAL(wp), ALLOCATABLE, TARGET :: &
      &  zrg_cosmu0(:,:),            & !< Cosine of solar zenith angle on reduced grid
      &  zrg_tsfc(:,:),              & !< Surface temperature on reduced grid
      &  zrg_emis_rad(:,:),          & !< Longwave surface emissivity on reduced grid
      &  zrg_albvisdir(:,:),         & !< Surface albedo for visible range (direct) on reduced grid
      &  zrg_albnirdir(:,:),         & !< Surface albedo for near IR range (direct) on reduced grid
      &  zrg_albvisdif(:,:),         & !< Surface albedo for visible range (diffuse) on reduced grid
      &  zrg_albnirdif(:,:),         & !< Surface albedo for near IR range (diffuse) on reduced grid
      &  zrg_pres(:,:,:),            & !< Pressure at full levels
      &  zrg_pres_ifc(:,:,:),        & !< Pressure at half levels
      &  zrg_temp(:,:,:),            & !< Temperature at full levels
      &  zrg_o3(:,:,:),              & !< Ozone mass mixing ratio on reduced grid
      &  zrg_tot_cld(:,:,:,:),       & !< Mass mixing ratio of water vapor, cloud water and cloud ice on reduced grid
      &  zrg_clc(:,:,:),             & !< Cloud cover on reduced grid
      &  zrg_fr_land(:,:),           & !< Land fraction on reduced grid
      &  zrg_fr_glac(:,:),           & !< Glacier fraction on reduced grid
      &  zrg_acdnc(:,:,:),           & !< Cloud droplet numb. conc. (m-3) on reduced grid
      &  zrg_aeq1(:,:,:),            & !< Climatological aerosol on reduced grid
      &  zrg_aeq2(:,:,:),            & !< Climatological aerosol on reduced grid
      &  zrg_aeq3(:,:,:),            & !< Climatological aerosol on reduced grid
      &  zrg_aeq4(:,:,:),            & !< Climatological aerosol on reduced grid
      &  zrg_aeq5(:,:,:),            & !< Climatological aerosol on reduced grid
      &  zrg_trsolall(:,:,:),        & !< solar transmissivity, all sky, net down on reduced grid
      &  zrg_lwflxall(:,:,:),        & !< Terrestrial flux, all sky, net down on reduced grid
      &  zrg_lwflx_up_sfc(:,:),      & !< Longwave upward flux at surface on reduced grid
      &  zrg_trsol_up_toa(:,:),      & !< Upward solar transmissivity at TOA on reduced grid
      &  zrg_trsol_up_sfc(:,:),      & !< Upward solar transmissivity at surface on reduced grid
      &  zrg_trsol_dn_sfc_diff(:,:), & !< Downward diffuse solar transmissivity at surface on reduced grid
      &  zrg_trsol_clr_sfc(:,:),     & !< Clear-sky net transmissvity at surface on reduced grid
      &  zrg_aclcov(:,:),            & !< Cloud cover on reduced grid
      &  zrg_trsol_par_sfc(:,:),     & !< Photosynthetically active radiation
      &  zrg_lwflx_clr_sfc(:,:),     & !< clear-sky net LW flux at surface
      &  zrg_lwflx_up    (:,:,:),    & !< longwave  3D upward   flux          
      &  zrg_lwflx_dn    (:,:,:),    & !< longwave  3D downward flux           
      &  zrg_swflx_up    (:,:,:),    & !< shortwave 3D upward   flux          
      &  zrg_swflx_dn    (:,:,:),    & !< shortwave 3D downward flux          
      &  zrg_lwflx_up_clr(:,:,:),    & !< longwave  3D upward   flux clear-sky
      &  zrg_lwflx_dn_clr(:,:,:),    & !< longwave  3D downward flux clear-sky
      &  zrg_swflx_up_clr(:,:,:),    & !< shortwave 3D upward   flux clear-sky
      &  zrg_swflx_dn_clr(:,:,:)       !< shortwave 3D downward flux clear-sky
    ! Some unused variables to be up- and downscaled (to not change the interface to up- and downscale)
    REAL(wp), ALLOCATABLE, TARGET :: &
      &  aclcov(:,:),                & !< Cloud cover
      &  zrg_albdif(:,:),            &
      &  zrg_rtype(:,:),             &
      &  zlp_pres_ifc(:,:,:),        &
      &  zlp_tot_cld(:,:,:,:)
    LOGICAL, ALLOCATABLE          :: &
      &  cosmu0mask(:)                 !< Mask if cosmu0 > 0

    jg         = pt_patch%id
    nlev       = pt_patch%nlev

    fact_reffc = (3.0e-9_wp/(4.0_wp*pi*rhoh2o))**(1.0_wp/3.0_wp)

    IF (msg_level >= 7) &
      &       CALL message(routine, 'ecrad radiation on reduced grid')


    IF (jg == 1 .AND. .NOT. l_limited_area) THEN
      ptr_pp      => pt_par_patch
      nblks_par_c =  pt_par_patch%nblks_c
      nblks_lp_c  =  p_patch_local_parent(jg)%nblks_c
    ELSE ! Nested domain with MPI parallelization
      ptr_pp      => p_patch_local_parent(jg)
      nblks_par_c =  ptr_pp%nblks_c
      nblks_lp_c  =  ptr_pp%nblks_c
    ENDIF

    ! Add extra layer for atmosphere above model top if requested
    IF (atm_phy_nwp_config(jg)%latm_above_top) THEN
      IF (jg == 1 .OR. pt_patch%nshift == 0) THEN
        nlev_rg = nlev + 1
      ELSE ! add a specified number levels up to the top of the parent domain in case of vertical nesting
        nlev_rg = MIN(nlev+nexlevs_rrg_vnest, pt_par_patch%nlev)
      ENDIF
    ELSE
      nlev_rg = nlev
    ENDIF
    nlev_rgp1 = nlev_rg+1

    CALL ecrad_single_level%allocate(nproma, 2, 1, .true.) !< use_sw_albedo_direct, 2 bands
    ecrad_single_level%solar_irradiance = 1._wp            !< Obtain normalized fluxes which corresponds to the 
                                                           !< transmissivity needed in the following

    IF (ecrad_conf%use_spectral_solar_scaling) THEN
      ALLOCATE(ecrad_single_level%spectral_solar_scaling(ecrad_conf%n_bands_sw))
      ecrad_single_level%spectral_solar_scaling = (/  1.0_wp, 1.0_wp, 1.0_wp, 1.0478_wp, 1.0404_wp, 1.0317_wp, &
         &   1.0231_wp, 1.0054_wp, 0.98413_wp, 0.99863_wp, 0.99907_wp, 0.90589_wp, 0.92213_wp, 1.0_wp /)
    ENDIF

    CALL ecrad_thermodynamics%allocate(nproma, nlev_rg, use_h2o_sat=.false., rrtm_pass_temppres_fl=.true.)

    CALL ecrad_gas%allocate(nproma, nlev_rg)
    
    CALL ecrad_cloud%allocate(nproma, nlev_rg)
    ! Currently hardcoded values for FSD
    CALL ecrad_cloud%create_fractional_std(nproma, nlev_rg, 1._wp)

    IF ( ecrad_conf%use_aerosols ) THEN
      ! Allocate aerosol container
      CALL ecrad_aerosol%allocate_direct(ecrad_conf, nproma, 1, nlev_rg)
    ENDIF

    CALL ecrad_flux%allocate(ecrad_conf, 1, nproma, nlev_rg)

    ALLOCATE(cosmu0mask(nproma))

    ! Allocate for reduced radiation grid
    ALLOCATE(zrg_cosmu0           (nproma,nblks_par_c),     &
      &      zrg_tsfc             (nproma,nblks_par_c),     &
      &      zrg_fr_land          (nproma,nblks_par_c),     &
      &      zrg_fr_glac          (nproma,nblks_par_c),     &
      &      zrg_emis_rad         (nproma,nblks_par_c),     &
      &      zrg_albvisdir        (nproma,nblks_par_c),     &
      &      zrg_albnirdir        (nproma,nblks_par_c),     &
      &      zrg_albvisdif        (nproma,nblks_par_c),     &
      &      zrg_albnirdif        (nproma,nblks_par_c),     &
      &      zrg_aclcov           (nproma,nblks_par_c),     &
      &      zrg_lwflx_up_sfc     (nproma,nblks_par_c),     &
      &      zrg_trsol_up_toa     (nproma,nblks_par_c),     &
      &      zrg_trsol_up_sfc     (nproma,nblks_par_c),     &
      &      zrg_trsol_par_sfc    (nproma,nblks_par_c),     &
      &      zrg_trsol_dn_sfc_diff(nproma,nblks_par_c),     &
      &      zrg_trsol_clr_sfc    (nproma,nblks_par_c),     &
      &      zrg_lwflx_clr_sfc    (nproma,nblks_par_c),     &
      &      aclcov               (nproma,pt_patch%nblks_c))
      
    ! Set dimensions for 3D radiative flux variables
    IF (atm_phy_nwp_config(jg)%l_3d_rad_fluxes) THEN
       np = nproma
       nl = nlev_rgp1
    ELSE
       np = 1
       nl = 1
    END IF

    ALLOCATE(zrg_pres_ifc    (nproma,nlev_rgp1,nblks_par_c),&
      &      zrg_lwflxall    (nproma,nlev_rgp1,nblks_par_c),&
      &      zrg_trsolall    (nproma,nlev_rgp1,nblks_par_c),&
      &      zrg_lwflx_up    (np, nl, nblks_par_c),&
      &      zrg_lwflx_dn    (np, nl, nblks_par_c),&   
      &      zrg_swflx_up    (np, nl, nblks_par_c),&
      &      zrg_swflx_dn    (np, nl, nblks_par_c),&
      &      zrg_lwflx_up_clr(np, nl, nblks_par_c),&
      &      zrg_lwflx_dn_clr(np, nl, nblks_par_c),&
      &      zrg_swflx_up_clr(np, nl, nblks_par_c),&
      &      zrg_swflx_dn_clr(np, nl, nblks_par_c) )

    ALLOCATE(zrg_pres     (nproma,nlev_rg  ,nblks_par_c),   &
      &      zrg_temp     (nproma,nlev_rg  ,nblks_par_c),   &
      &      zrg_o3       (nproma,nlev_rg  ,nblks_par_c),   &
      &      zrg_aeq1     (nproma,nlev_rg  ,nblks_par_c),   &
      &      zrg_aeq2     (nproma,nlev_rg  ,nblks_par_c),   &
      &      zrg_aeq3     (nproma,nlev_rg  ,nblks_par_c),   &
      &      zrg_aeq4     (nproma,nlev_rg  ,nblks_par_c),   &
      &      zrg_aeq5     (nproma,nlev_rg  ,nblks_par_c),   &
      &      zrg_acdnc    (nproma,nlev_rg  ,nblks_par_c),   &
      &      zrg_clc      (nproma,nlev_rg  ,nblks_par_c))

    ALLOCATE(zrg_tot_cld  (nproma,nlev_rg  ,nblks_par_c,3))

    ! Unused variables, allocated to not change the upscale_rad_input interface
    ALLOCATE(zrg_albdif   (nproma,nblks_par_c),             &
      &      zrg_rtype    (nproma,nblks_par_c),             &
      &      zlp_pres_ifc (nproma,nlev_rgp1,nblks_lp_c ),   &
      &      zlp_tot_cld  (nproma,nlev_rg  ,nblks_lp_c,3))



    rl_start = 1 ! SR radiation is not set up to handle boundaries of nested domains
    rl_end   = min_rlcell_int
    i_startblk = pt_patch%cells%start_block(rl_start)
    i_endblk   = pt_patch%cells%end_block(rl_end)

!$OMP PARALLEL PRIVATE(jb,i_startidx,i_endidx)

    ! Initialize output fields
    CALL init(zrg_trsolall(:,:,:), 0._wp)
    CALL init(zrg_lwflxall(:,:,:), 0._wp)
    CALL init(zrg_trsol_up_toa(:,:)     )
    CALL init(zrg_trsol_up_sfc(:,:)     )
    CALL init(zrg_trsol_par_sfc(:,:)    )
    CALL init(zrg_trsol_dn_sfc_diff(:,:))
    CALL init(zrg_trsol_clr_sfc(:,:)    )
    CALL init(zrg_lwflx_up_sfc(:,:)     )
    CALL init(zrg_lwflx_clr_sfc(:,:)    )
    CALL init(zrg_aclcov(:,:)           )
 
    IF (atm_phy_nwp_config(jg)%l_3d_rad_fluxes) THEN
      CALL init(zrg_lwflx_up    (:,:,:), 0._wp)   
      CALL init(zrg_lwflx_dn    (:,:,:), 0._wp)     
      CALL init(zrg_swflx_up    (:,:,:), 0._wp)  
      CALL init(zrg_swflx_dn    (:,:,:), 0._wp) 
      CALL init(zrg_lwflx_up_clr(:,:,:), 0._wp)
      CALL init(zrg_lwflx_dn_clr(:,:,:), 0._wp)
      CALL init(zrg_swflx_up_clr(:,:,:), 0._wp)
      CALL init(zrg_swflx_dn_clr(:,:,:), 0._wp)
    END IF

!$OMP DO ICON_OMP_GUIDED_SCHEDULE
    DO jb = i_startblk, i_endblk
      CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
        &                       i_startidx, i_endidx, rl_start, rl_end)
      prm_diag%tsfctrad(i_startidx:i_endidx,jb) = lnd_prog%t_g(i_startidx:i_endidx,jb)
    ENDDO ! jb
!$OMP END DO NOWAIT

!$OMP END PARALLEL

! Upscale ICON input fields from full grid to reduced radiation grid
    CALL upscale_rad_input(pt_patch%id, pt_par_patch%id,                                 &
      &                    nlev_rg, ext_data%atm%fr_land_smt, ext_data%atm%fr_glac_smt,  &
      &                    prm_diag%lw_emiss, prm_diag%cosmu0,                           &
      &                    prm_diag%albvisdir, prm_diag%albnirdir, prm_diag%albvisdif,   &
      &                    prm_diag%albnirdif, prm_diag%albdif, prm_diag%tsfctrad,       &
      &                    prm_diag%ktype, pt_diag%pres_ifc, pt_diag%pres,               &
      &                    pt_diag%temp,prm_diag%acdnc, prm_diag%tot_cld, prm_diag%clc,  &
      &                    ext_data%atm%o3, zaeq1, zaeq2, zaeq3, zaeq4, zaeq5,           &
      &                    zrg_fr_land, zrg_fr_glac, zrg_emis_rad,                       &
      &                    zrg_cosmu0, zrg_albvisdir, zrg_albnirdir, zrg_albvisdif,      &
      &                    zrg_albnirdif, zrg_albdif, zrg_tsfc, zrg_rtype, zrg_pres_ifc, &
      &                    zrg_pres, zrg_temp, zrg_acdnc, zrg_tot_cld, zrg_clc, zrg_o3,  &
      &                    zrg_aeq1, zrg_aeq2, zrg_aeq3, zrg_aeq4, zrg_aeq5,             &
      &                    zlp_pres_ifc, zlp_tot_cld, prm_diag%buffer_rrg)

! Set indices for reduced grid loop
    IF (jg == 1 .AND. l_limited_area) THEN
      rl_start = grf_fbk_start_c
    ELSE
      rl_start = grf_ovlparea_start_c
    ENDIF
    rl_end     = min_rlcell_int
    i_startblk = ptr_pp%cells%start_block(rl_start)
    i_endblk   = ptr_pp%cells%end_block(rl_end)

!$OMP PARALLEL PRIVATE(jb,jc,i_startidx,i_endidx,cosmu0mask)                         &
!$OMP          FIRSTPRIVATE(ecrad_aerosol,ecrad_single_level, ecrad_thermodynamics,  &
!$OMP                       ecrad_gas, ecrad_cloud, ecrad_flux)
!$OMP DO ICON_OMP_GUIDED_SCHEDULE
    DO jb = i_startblk, i_endblk
      CALL get_indices_c(ptr_pp, jb, i_startblk, i_endblk, &
        &                i_startidx, i_endidx, rl_start, rl_end)

      ! It may happen that an MPI patch contains only nest boundary points
      ! In this case, no action is needed
      IF (i_startidx > i_endidx) CYCLE

      ! GZ, 2020-05-08: provisional workaround for indexing error hidden somewhere in ecRad:
      ! Start indices larger than 1 lead to erroneous access of array elements and break processor invariance.
      ! This problem was already present in ecRad version 1.3 prior to changes/optimizations implemented at DWD
      !
      IF (i_startidx > 1) THEN
        DO jc = 1, i_startidx-1
          zrg_cosmu0(jc,jb)    = zrg_cosmu0(i_startidx,jb)
          zrg_tsfc(jc,jb)      = zrg_tsfc(i_startidx,jb)
          zrg_albvisdif(jc,jb) = zrg_albvisdif(i_startidx,jb)
          zrg_albnirdif(jc,jb) = zrg_albnirdif(i_startidx,jb)
          zrg_albvisdir(jc,jb) = zrg_albvisdir(i_startidx,jb)
          zrg_albnirdir(jc,jb) = zrg_albnirdir(i_startidx,jb)
          zrg_fr_land(jc,jb)   = zrg_fr_land(i_startidx,jb)
          zrg_fr_glac(jc,jb)   = zrg_fr_glac(i_startidx,jb)
          zrg_emis_rad(jc,jb)  = zrg_emis_rad(i_startidx,jb)
          zrg_tsfc(jc,jb)      = zrg_tsfc(i_startidx,jb)
          zrg_pres_ifc(jc,nlev_rgp1,jb) = zrg_pres_ifc(i_startidx,nlev_rgp1,jb)
        ENDDO
        DO jk = 1,nlev_rg
          DO jc = 1, i_startidx-1
            zrg_temp(jc,jk,jb)        = zrg_temp(i_startidx,jk,jb)
            zrg_pres(jc,jk,jb)        = zrg_pres(i_startidx,jk,jb)
            zrg_pres_ifc(jc,jk,jb)    = zrg_pres_ifc(i_startidx,jk,jb)
            zrg_o3(jc,jk,jb)          = zrg_o3(i_startidx,jk,jb)
            zrg_tot_cld(jc,jk,jb,iqv) = zrg_tot_cld(i_startidx,jk,jb,iqv)
            zrg_tot_cld(jc,jk,jb,iqc) = zrg_tot_cld(i_startidx,jk,jb,iqc)
            zrg_tot_cld(jc,jk,jb,iqi) = zrg_tot_cld(i_startidx,jk,jb,iqi)
            zrg_clc(jc,jk,jb)         = zrg_clc(i_startidx,jk,jb)
            zrg_acdnc(jc,jk,jb)       = zrg_acdnc(i_startidx,jk,jb)
            zrg_aeq1(jc,jk,jb)        = zrg_aeq1(i_startidx,jk,jb)
            zrg_aeq2(jc,jk,jb)        = zrg_aeq2(i_startidx,jk,jb)
            zrg_aeq3(jc,jk,jb)        = zrg_aeq3(i_startidx,jk,jb)
            zrg_aeq4(jc,jk,jb)        = zrg_aeq4(i_startidx,jk,jb)
            zrg_aeq5(jc,jk,jb)        = zrg_aeq5(i_startidx,jk,jb)
          ENDDO
        ENDDO
      ENDIF
      i_startidx = 1
      !
      ! end of workaround

      cosmu0mask(:) = .FALSE.
      DO jc = i_startidx, i_endidx
        IF ( zrg_cosmu0(jc,jb) > 0._wp ) THEN
          cosmu0mask(jc) = .TRUE.
        ENDIF
      ENDDO

! Fill single level configuration type
      CALL ecrad_set_single_level(ecrad_single_level, current_datetime, ptr_pp%cells%center(:,jb),            &
        &                         zrg_cosmu0(:,jb), zrg_tsfc(:,jb), zrg_albvisdif(:,jb), zrg_albnirdif(:,jb), &
        &                         zrg_albvisdir(:,jb), zrg_albnirdir(:,jb), zrg_emis_rad(:,jb),               &
        &                         i_startidx, i_endidx)

! Fill thermodynamics configuration type
      CALL ecrad_set_thermodynamics(ecrad_thermodynamics, zrg_temp(:,:,jb), zrg_pres(:,:,jb),     &
        &                           zrg_pres_ifc(:,:,jb), zrg_tsfc(:,jb),                         &
        &                           nlev_rg, nlev_rgp1, i_startidx, i_endidx)

! Fill gas configuration type
      CALL ecrad_set_gas(ecrad_gas, ecrad_conf, zrg_o3(:,:,jb), zrg_tot_cld(:,:,jb,iqv), &
        &                zrg_pres(:,:,jb), i_startidx, i_endidx, nlev_rg)

! Fill clouds configuration type
      CALL ecrad_set_clouds(ecrad_cloud, ecrad_thermodynamics, zrg_tot_cld(:,:,jb,iqc),  &
        &                   zrg_tot_cld(:,:,jb,iqi), zrg_clc(:,:,jb),                    &
        &                   zrg_temp(:,:,jb), zrg_pres(:,:,jb), zrg_acdnc(:,:,jb),       &
        &                   zrg_fr_glac(:,jb), zrg_fr_land(:,jb),                        &
        &                   fact_reffc, ecrad_conf%cloud_fraction_threshold, nlev_rg, i_startidx, i_endidx)

! Fill aerosol configuration type
      SELECT CASE (irad_aero)
        CASE(0)
          ! No aerosol, nothing to do
        CASE(2)
          ! Case 2: Constant aerosol
          !         Arguments can be added to fill ecrad_aerosol with actual values. For the time being,
          !         we stay consistent with RRTM where irad_aero=2 does not add any aerosol
          CALL nwp_ecrad_prep_aerosol(ecrad_conf, ecrad_aerosol)
        CASE(5,6)
          ! Fill aerosol configuration type with Tanre or Tegen aerosol
          CALL nwp_ecrad_prep_aerosol(1, nlev_rg, i_startidx, i_endidx,     &
            &                         zrg_aeq1(:,:,jb), zrg_aeq2(:,:,jb),   &
            &                         zrg_aeq3(:,:,jb), zrg_aeq4(:,:,jb),   &
            &                         zrg_aeq5(:,:,jb),                     &
            &                         ecrad_conf, ecrad_aerosol)
        CASE DEFAULT
          CALL finish(TRIM(routine),'irad_aero not valid for ecRad')
      END SELECT

      ! Reset output values
      ecrad_flux%cloud_cover_sw(:) = 0._wp
      ecrad_flux%cloud_cover_lw(:) = 0._wp

!---------------------------------------------------------------------------------------
! Call the radiation scheme ecRad
!---------------------------------------------------------------------------------------
      CALL ecrad(nproma, nlev_rg, i_startidx, i_endidx,   & !< Array and loop bounds (input)
        &        ecrad_conf,                              & !< General ecRad configuration object (input)
        &        ecrad_single_level,                      & !< ecRad single level configuration object (input)
        &        ecrad_thermodynamics,                    & !< ecRad thermodynamics configuration object (input)
        &        ecrad_gas,                               & !< ecRad gas configuration object (input)
        &        ecrad_cloud,                             & !< ecRad cloud configuration object (input)
        &        ecrad_aerosol,                           & !< ecRad aerosol configuration object (input)
        &        ecrad_flux                               ) !< ecRad fluxes (output)
!---------------------------------------------------------------------------------------

! Update ICON variables with fluxes from ecRad
      CALL ecrad_store_fluxes(jg, ecrad_flux, zrg_cosmu0(:,jb), zrg_trsolall   (:,:,jb),    &
        &                     zrg_trsol_up_toa  (:,jb), zrg_trsol_up_sfc     (:,jb),    &
        &                     zrg_trsol_par_sfc (:,jb), zrg_trsol_dn_sfc_diff(:,jb),    &
        &                     zrg_trsol_clr_sfc (:,jb), zrg_lwflxall       (:,:,jb),    &
        &                     zrg_lwflx_up_sfc  (:,jb), zrg_lwflx_clr_sfc    (:,jb),    &
        &                     zrg_lwflx_up    (:,:,jb), zrg_lwflx_dn       (:,:,jb),    &
        &                     zrg_swflx_up    (:,:,jb), zrg_swflx_dn       (:,:,jb),    &
        &                     zrg_lwflx_up_clr(:,:,jb), zrg_lwflx_dn_clr   (:,:,jb),    &
        &                     zrg_swflx_up_clr(:,:,jb), zrg_swflx_dn_clr   (:,:,jb),    &  
        &                     cosmu0mask(:), i_startidx, i_endidx, nlev_rgp1)

      ! Add 3D contribution to diffuse radiation
      CALL add_3D_diffuse_rad(ecrad_flux, zrg_clc(:,:,jb), zrg_pres(:,:,jb), zrg_temp(:,:,jb),            &
        &                     zrg_cosmu0(:,jb), zrg_trsol_dn_sfc_diff(:,jb), i_startidx, i_endidx, nlev_rg)

    ENDDO !jb
!$OMP END DO NOWAIT
!$OMP END PARALLEL

! Downscale radiative fluxes from reduced radiation grid to full grid
    CALL downscale_rad_output(pt_patch%id, pt_par_patch%id,                                         &
      &  nlev_rg, zrg_aclcov, zrg_lwflxall, zrg_trsolall, zrg_trsol_clr_sfc, zrg_lwflx_clr_sfc,     &
      &  zrg_lwflx_up_sfc, zrg_trsol_up_toa, zrg_trsol_up_sfc, zrg_trsol_par_sfc,                   &
      &  zrg_trsol_dn_sfc_diff, zrg_tsfc, zrg_albdif, zrg_emis_rad, zrg_cosmu0, zrg_tot_cld,        &
      &  zlp_tot_cld, zrg_pres_ifc, zlp_pres_ifc, prm_diag%tsfctrad, prm_diag%albdif, aclcov,       &
      &  prm_diag%lwflxall, prm_diag%trsolall, prm_diag%lwflx_up_sfc_rs, prm_diag%trsol_up_toa,     &
      &  prm_diag%trsol_up_sfc, prm_diag%trsol_par_sfc, prm_diag%trsol_dn_sfc_diff,                 &
      &  prm_diag%trsolclr_sfc, prm_diag%lwflxclr_sfc,                                              & 
      &  zrg_lwflx_up         , zrg_lwflx_dn         , zrg_swflx_up         , zrg_swflx_dn,         &
      &  zrg_lwflx_up_clr     , zrg_lwflx_dn_clr     , zrg_swflx_up_clr     , zrg_swflx_dn_clr,     &
      &  prm_diag%lwflx_up    , prm_diag%lwflx_dn    , prm_diag%swflx_up    , prm_diag%swflx_dn,    &
      &  prm_diag%lwflx_up_clr, prm_diag%lwflx_dn_clr, prm_diag%swflx_up_clr, prm_diag%swflx_dn_clr )

! CLEANUP
    CALL ecrad_single_level%deallocate
    CALL ecrad_thermodynamics%deallocate
    CALL ecrad_gas%deallocate
    CALL ecrad_cloud%deallocate
    IF ( ecrad_conf%use_aerosols ) CALL ecrad_aerosol%deallocate
    CALL ecrad_flux%deallocate

    DEALLOCATE (zrg_cosmu0, zrg_tsfc, zrg_emis_rad, zrg_albvisdir, zrg_albnirdir, zrg_albvisdif,   &
      &         zrg_albnirdif, zrg_pres_ifc, zrg_o3, zrg_aeq1, zrg_aeq2, zrg_aeq3, zrg_clc,        &
      &         zrg_aeq4, zrg_aeq5, zrg_tot_cld, zrg_pres, zrg_temp, zrg_trsolall, zrg_lwflxall,   &
      &         zrg_lwflx_up_sfc, zrg_trsol_up_toa, zrg_trsol_up_sfc, zrg_trsol_dn_sfc_diff,       &
      &         zrg_trsol_clr_sfc, zrg_aclcov, zrg_trsol_par_sfc, zrg_fr_land, zrg_fr_glac, aclcov,&
      &         zrg_acdnc, zrg_albdif, zrg_rtype, zlp_pres_ifc, zlp_tot_cld, zrg_lwflx_clr_sfc,    &
      &         zrg_lwflx_up    , zrg_lwflx_dn    , zrg_swflx_up    , zrg_swflx_dn,                &
      &         zrg_lwflx_up_clr, zrg_lwflx_dn_clr, zrg_swflx_up_clr, zrg_swflx_dn_clr             )

    DEALLOCATE(cosmu0mask)

  END SUBROUTINE nwp_ecrad_radiation_reduced
  !---------------------------------------------------------------------------------------

#endif
END MODULE mo_nwp_ecrad_interface
