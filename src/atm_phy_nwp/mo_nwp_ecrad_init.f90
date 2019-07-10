!>
!! In this module the configuration state for the ecRad radiation code is being set up.
!!
!! - Setup information is stored inside the object ecrad_conf of the derived type t_ecrad_conf
!!   containing all the configuration information needed to run the radiation scheme.
!! - The intention is that this part of the configuration is fixed for a given model run.
!! - ICON namelist settings are translated to ecRad conform settings, if unsupported values
!!   are provided via namelist, the user gets an error. (These values should already be
!!   checked by the nml_crosscheck)
!! - Currently, only the McICA Solver is supported.
!! - Please note that only a subset of the available configuration options of ecRad is
!!   filled by this routine. E.g., options only connected to the SPARTACUS Solver are 
!!   currently not changed from the default as the SPARTACUS Solver was not tested in
!!   ICON so far. For a full list of ecRad settings, please have a look at
!!   externals/ecrad/radiation/radiation_config.F90
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
MODULE mo_nwp_ecrad_init

  USE mo_kind,                 ONLY: wp
  USE mo_exception,            ONLY: finish, message, message_text
  USE mo_impl_constants,       ONLY: MAX_CHAR_LENGTH
  USE mo_radiation_config,     ONLY: icld_overlap, irad_aero, ecrad_data_path,           &
                                 &   llw_cloud_scat, iliquid_scat, iice_scat
#ifdef __ECRAD
  USE mo_ecrad,                ONLY: t_ecrad_conf, ecrad_setup,                          &
                                 &   ISolverHomogeneous, ISolverMcICA, ISolverSpartacus, &
                                 &   ISolverTripleclouds,                                &
                                 &   IGasModelMonochromatic, IGasModelIFSRRTMG,          &
                                 &   IGasModelPSRRTMG, ILiquidModelSOCRATES,             &
                                 &   ILiquidModelSlingo, IIceModelFu, IIceModelBaran2016,&
                                 &   IOverlapMaximumRandom, IOverlapExponentialRandom,   &
                                 &   IOverlapExponential,                                &
                                 &   nweight_par_ecrad, iband_par_ecrad, weight_par_ecrad
#endif


  IMPLICIT NONE

  PRIVATE

  !> module name string
  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_nwp_ecrad_init'

#ifdef __ECRAD
  PUBLIC :: setup_ecrad


CONTAINS


  !---------------------------------------------------------------------------------------
  !>
  !! @par Revision History
  !! Initial release by Daniel Rieger, Deutscher Wetterdienst, Offenbach (2019-01-31)
  !!
  SUBROUTINE setup_ecrad ( ecrad_conf )

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
      &  routine = modname//'::setup_ecrad' 

    TYPE(t_ecrad_conf), INTENT(inout) :: &
      &  ecrad_conf           !< ecRad configuration state

    ! Local variables
    REAL(wp)                 :: &
      &  wavelength_bound(1)      !< Wavelength bound between VIS and NIR albedo (m)
    INTEGER                  :: &
      &  i_band_in(2)             !< The albedo band indices corresponding to each interval

    WRITE (message_text,'(A)') 'Setup of ecRad'
    CALL message('',message_text)

    !---------------------------------------------------------------------------------------
    ! Checks
    !---------------------------------------------------------------------------------------

    ! Compatibility check wp and JPRB. If this check fails, JPRB has to be adapted manually to wp.
    IF (PRECISION(wavelength_bound) /= PRECISION(ecrad_conf%cloud_fraction_threshold)) &
      &  CALL finish(TRIM(routine),'ICON working precision (wp) does not match ecRad precision.')
    IF (EPSILON(wavelength_bound(1)) /= EPSILON(ecrad_conf%cloud_fraction_threshold)) &
      &  CALL finish(TRIM(routine),'Smallest number in working precision (wp) is different from ecRad precision.')

    !---------------------------------------------------------------------------------------
    ! Configuration based on ICON namelist settings
    !---------------------------------------------------------------------------------------

    ! Directory with all input data required by ecRad
    ecrad_conf%directory_name = TRIM(ecrad_data_path)

    ! Overlap scheme
    SELECT CASE (icld_overlap)
      CASE (1)
        ecrad_conf%i_overlap_scheme = IOverlapMaximumRandom
      CASE (2)
        ecrad_conf%i_overlap_scheme = IOverlapExponentialRandom
      CASE (5)
        ecrad_conf%i_overlap_scheme = IOverlapExponential
      CASE DEFAULT
        CALL finish(TRIM(routine),'Only values of 1 (MAX-RAN), 2 (EXP-RAN) and 5 (EXP) are valid for icld_overlap')
    END SELECT

    ! Aerosol climatology
    SELECT CASE (irad_aero)
      CASE (0) ! No aerosol
        ecrad_conf%use_aerosols = .false.
      CASE (2,5,6) ! Constant, Tanre, Tegen
        ecrad_conf%use_aerosols = .true.
      CASE DEFAULT
        CALL finish(TRIM(routine),'irad_aero not valid for ecRad')
    END SELECT

    ! LW scattering due to clouds
    ecrad_conf%do_lw_cloud_scattering = llw_cloud_scat
    
    ! Liquid cloud particle scattering properties
    SELECT CASE (iliquid_scat)
      CASE(0)
        ecrad_conf%i_liq_model = ILiquidModelSOCRATES
      CASE(1)
        ecrad_conf%i_liq_model = ILiquidModelSlingo
      CASE DEFAULT
        CALL finish(TRIM(routine),'i_liquid_scat not valid for ecRad')
    END SELECT
    
    ! Ice cloud particle scattering properties
    SELECT CASE (iice_scat)
      CASE(0)
        ecrad_conf%i_ice_model = IIceModelFu
      CASE(1)
        ecrad_conf%i_ice_model = IIceModelBaran2016
      CASE DEFAULT
        CALL finish(TRIM(routine),'i_ice_scat not valid for ecRad')
    END SELECT

    !---------------------------------------------------------------------------------------
    ! Currently hardcoded configuration
    !---------------------------------------------------------------------------------------
  
    ecrad_conf%do_lw                       = .true.       !< Do we compute longwave radiation?
    !
    ecrad_conf%do_sw                       = .true.       !< Do we compute shortwave radiation?
    !
    ecrad_conf%do_clear                    = .true.       !< Do we compute clear-sky fluxes?
    !
    ecrad_conf%do_sw_direct                = .true.       !< Do we compute solar direct fluxes?
    !
    ecrad_conf%do_3d_effects               = .false.      !< Do we include 3D effects?
    !
    ecrad_conf%do_lw_aerosol_scattering    = .false.      !< LW scattering due to aerosol
    !
    ecrad_conf%use_beta_overlap            = .false.      !< Use Shonk et al. (2010) "beta" overlap parameter
                                                          !< instead of "alpha" (Hogan and Illingworth, 2000)
    !
    ecrad_conf%i_solver_sw                 = ISolverMcICA !< Short-wave solver
    !
    ecrad_conf%i_solver_lw                 = ISolverMcICA !< Long-wave solver
    !
    ecrad_conf%iverbosesetup               = 0            !< Verbosity (0: none,     1: warning,  2: info,
    ecrad_conf%iverbose                    = 0            !<            3: progress, 4: detailed, 5: debug)
    !
    ecrad_conf%do_surface_sw_spectral_flux = .true.       !< Save the surface downwelling shortwave fluxes in each band
                                                          !< Needed for photosynthetic active radiation
    !
    ecrad_conf%do_fu_lw_ice_optics_bug     = .false.      !< In the IFS environment there was a bug in the Fu longwave
                                                          !< ice optics producing better results than the fixed version
    !
    ecrad_conf%do_sw_delta_scaling_with_gases = .false.   !< Do SW delta-Eddington scaling on cloud-aerosol-gas mixture (.true.)
                                                          !< More correct approach of separately scaling the cloud and aerosol 
                                                          !< scattering properties before merging with gases (.false.)
    !
    ecrad_conf%i_gas_model           = IGasModelIFSRRTMG  !< Use RRTM gas model (only available option)
    IF (ecrad_conf%i_gas_model == IGasModelIFSRRTMG) THEN
      ecrad_conf%do_setup_ifsrrtm = .true.
    ELSE
      ecrad_conf%do_setup_ifsrrtm = .false.
    ENDIF
    !
    ecrad_conf%cloud_fraction_threshold    = 1.0e-6_wp    !< Cloud is present in a layer if cloud fraction exceeds this value
    ecrad_conf%cloud_mixing_ratio_threshold= 1.0e-9_wp    !< ..and total cloud water mixing ratio exceeds this value
    !
    ecrad_conf%cloud_inhom_decorr_scaling  = 0.5_wp       !< Ratio of the overlap decorrelation length for cloud inhomogeneities
                                                          !< to the overlap decorrelation length for cloud boundaries.
                                                          !< Observations suggest this has a value of 0.5
    !
    ecrad_conf%min_gas_od_lw               = 1.0e-15_wp   !< Minimum gas optical depth, for stability (long-wave)
    !
    ecrad_conf%min_gas_od_sw               = 0.0_wp       !< Minimum gas optical depth, for stability (short-wave)
    !
    ecrad_conf%max_cloud_od                = 20.0_wp      !< Maximum total optical depth of a cloudy region, for stability

    !---------------------------------------------------------------------------------------
    ! Call to ecRad setup routine. This also consolidates the configuration
    !---------------------------------------------------------------------------------------

    CALL ecrad_setup(ecrad_conf)

    !---------------------------------------------------------------------------------------
    ! Tell ecRad about the wavelength bounds of ICON data
    !---------------------------------------------------------------------------------------

    ! Set up the photosynthetically active radiation wavelength bounds
    CALL ecrad_conf%get_sw_weights(0.4e-6_wp, 0.7e-6_wp, &
      &  nweight_par_ecrad, iband_par_ecrad, weight_par_ecrad, &
      &  'photosynthetically active radiation, PAR')

    ! ICON external parameters have SW albedo for two different wavelength bands, visible and near infrared. The following call to
    ! ecrad_conf%define_sw_albedo_intervals tells ecrad about the two bands and the wavelength bound which is at 700 nm (according
    ! to a comment in mo_nwp_phy_types).
    wavelength_bound(1) = 0.7_wp * 1.e-6_wp !< 700 nm
    i_band_in(1)        = 1
    i_band_in(2)        = 2
    CALL ecrad_conf%define_sw_albedo_intervals(2, wavelength_bound=wavelength_bound, i_band_in=i_band_in) 

    ! Similar to the short wave albedo bands, ecRad needs to know the number of longwave emissivity bands provided from ICON
    ! external data. As the number is 1, no other arguments are needed
    CALL ecrad_conf%define_lw_emiss_intervals(1)
    
  END SUBROUTINE setup_ecrad
  !---------------------------------------------------------------------------------------

#endif
END MODULE mo_nwp_ecrad_init
