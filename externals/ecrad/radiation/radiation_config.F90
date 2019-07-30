! radiation_config.f90 - Derived type to configure the radiation scheme
!
! Copyright (C) 2014-2018 ECMWF
!
! Author:  Robin Hogan
! Email:   r.j.hogan@ecmwf.int
! License: see the COPYING file for details
!
! Modifications
!   2017-07-22  R. Hogan  Added Yi et al. ice optics model
!   2017-10-23  R. Hogan  Renamed single-character variables
!   2018-03-15  R. Hogan  Added logicals controlling surface spectral treatment

module radiation_config

  use parkind1,                      only : jprb

  use radiation_cloud_optics_data,   only : cloud_optics_type
  use radiation_aerosol_optics_data, only : aerosol_optics_type
  use radiation_pdf_sampler,         only : pdf_sampler_type
  use radiation_cloud_cover,         only : OverlapName, &
       & IOverlapMaximumRandom, IOverlapExponentialRandom, IOverlapExponential

  implicit none

  ! Configuration codes: use C-style enumerators to avoid having to
  ! remember the numbers

  ! Solvers: can be specified for longwave and shortwave
  ! independently, except for "Homogeneous", which must be the same
  ! for both
  enum, bind(c) 
     enumerator ISolverHomogeneous, ISolverMcICA, ISolverSpartacus, &
          &     ISolverTripleclouds
  end enum
  character(len=*), parameter :: SolverName(0:3) = (/ 'Homogeneous ', &
       &                                              'McICA       ', &
       &                                              'SPARTACUS   ', &
       &                                              'Tripleclouds' /)

  ! SPARTACUS shortwave solver can treat the reflection of radiation
  ! back up into different regions in various ways: "minimum" permits
  ! lateral transport between regions but not within them, "maximum"
  ! assumes complete horizontal homogenization of radiation within
  ! regions (the original SPARTACUS assumption), and "computed"
  ! estimates the characteristic horizontal migration distance and
  ! uses this to estimate the fraction of radiation reflected back up
  ! into another region.
  enum, bind(c) 
     enumerator IEncroachmentMinimum, &
          &     IEncroachmentMaximum, IEncroachmentComputed
  end enum
  character(len=*), parameter :: EncroachmentName(0:2) = (/ 'Minimum ', &
       &                                                    'Maximum ', &
       &                                                    'Computed' /)

  ! Two-stream models
  ! This is not configurable at run-time

  ! Gas models
  enum, bind(c) 
     enumerator IGasModelMonochromatic, IGasModelIFSRRTMG, IGasModelPSRRTMG
  end enum
  character(len=*), parameter :: GasModelName(0:2) = (/ 'Monochromatic', &
       &                                                'RRTMG-IFS    ', &
       &                                                'RRTMG-PSRAD  ' /)

  ! Hydrometeor scattering models
  enum, bind(c) 
     enumerator ILiquidModelMonochromatic, ILiquidModelHuStamnesPSRAD, &
          &     ILiquidModelSOCRATES, ILiquidModelSlingo
  end enum
  character(len=*), parameter :: LiquidModelName(0:3) = (/ 'Monochromatic', &
       &                                                   'HuStamnes    ', &
       &                                                   'SOCRATES     ', &
       &                                                   'Slingo       ' /)

  enum, bind(c) 
     enumerator IIceModelMonochromatic, IIceModelFuPSRAD, IIceModelFu, &
          &  IIceModelBaran, IIceModelBaran2016, IIceModelBaran2017,   &
          &  IIceModelYi
  end enum
  character(len=*), parameter :: IceModelName(0:6) = (/ 'Monochromatic', &
       &                                                'Fu-PSRAD     ', &
       &                                                'Fu-IFS       ', &
       &                                                'Baran        ', &
       &                                                'Baran2016    ', &
       &                                                'Baran2017    ', &
       &                                                'Yi           ' /)

  ! Cloud PDF distribution shapes
  enum, bind(c)
    enumerator IPdfShapeLognormal, IPdfShapeGamma
  end enum
  character(len=*), parameter :: PdfShapeName(0:1) = (/ 'Lognormal', &
       &                                                'Gamma    ' /)


  integer, parameter :: NMaxAerosolTypes = 256

  ! Length of string buffer for printing config information
  integer, parameter :: NPrintStringLen = 60

  !---------------------------------------------------------------------
  ! Derived type containing all the configuration information needed
  ! to run the radiation scheme.  The intention is that this is fixed
  ! for a given model run.  The parameters are to list first those
  ! quantities that can be set directly by the user, for example using a
  ! namelist, and second those quantities that are computed afterwards
  ! from the user-supplied numbers, especially the details of the gas
  ! optics model.
  type config_type
    ! USER-CONFIGURABLE PARAMETERS

    ! Override default solar spectrum
    logical :: use_spectral_solar_scaling = .false.

    ! Directory in which gas, cloud and aerosol data files are to be
    ! found
    character(len=511) :: directory_name = '.'

    ! Cloud is deemed to be present in a layer if cloud fraction
    ! exceeds this value
    real(jprb) :: cloud_fraction_threshold = 1.0e-6_jprb
    ! ...and total cloud water mixing ratio exceeds this value
    real(jprb) :: cloud_mixing_ratio_threshold = 1.0e-9_jprb

    ! Overlap scheme
    integer :: i_overlap_scheme = IOverlapExponentialRandom

    ! Use the Shonk et al. (2010) "beta" overlap parameter, rather
    ! than the "alpha" overlap parameter of Hogan and Illingworth
    ! (2000)?
    logical :: use_beta_overlap = .false.

    ! Shape of sub-grid cloud water PDF
    integer :: i_cloud_pdf_shape = IPdfShapeGamma

    ! The ratio of the overlap decorrelation length for cloud
    ! inhomogeneities to the overlap decorrelation length for cloud
    ! boundaries.  Observations suggest this has a value of 0.5
    ! (e.g. from the decorrelation lengths of Hogan and Illingworth
    ! 2003 and Hogan and Illingworth 2000).
    real(jprb) :: cloud_inhom_decorr_scaling = 0.5_jprb

    ! Factor controlling how much of the cloud edge length interfaces
    ! directly between the clear-sky region (region a) and the
    ! optically thick cloudy region (region c).  If Lxy is the length
    ! of the interfaces between regions x and y, and Lab and Lbc have
    ! been computed already, then
    !   Lac=clear_to_thick_fraction*min(Lab,Lbc).
    real(jprb) :: clear_to_thick_fraction = 0.0_jprb

    ! Factor allowing lateral transport when the sun is close to
    ! overhead; consider atand(overhead_sun_factor) to be the number
    ! of degrees that the sun angle is perturbed from zenith for the
    ! purposes of computing lateral transport.  A value of up to 0.1
    ! seems to be necessary to account for the fact that some forward
    ! scattered radiation is treated as unscattered by delta-Eddington
    ! scaling; therefore it ought to have the chance to escape.
    real(jprb) :: overhead_sun_factor = 0.0_jprb

    ! Minimum gas optical depth in a single layer at any wavelength,
    ! for stability
    real(jprb) :: min_gas_od_lw = 1.0e-15_jprb
    real(jprb) :: min_gas_od_sw = 0.0_jprb

    ! Maximum gas optical depth in a layer before that g-point will
    ! not be considered for 3D treatment: a limit is required to avoid
    ! expensive computation of matrix exponentials on matrices with
    ! large elements
    real(jprb) :: max_gas_od_3d = 8.0_jprb

    ! Maximum total optical depth of a cloudy region for stability:
    ! optical depth will be capped at this value
    real(jprb) :: max_cloud_od = 18.0_jprb

    ! How much longwave scattering is included?
    logical :: do_lw_cloud_scattering = .true.
    logical :: do_lw_aerosol_scattering = .true.

    ! Number of regions used to describe clouds and clear skies. A
    ! value of 2 means one clear and one cloudy region, so clouds are
    ! horizontally homogeneous, while a value of 3 means two cloudy
    ! regions with different optical depth, thereby representing
    ! inhomogeneity via the Shonk & Hogan (2008) "Tripleclouds"
    ! method.
    integer :: nregions = 3

    ! Code specifying the solver to be used: use the enumerations
    ! defined above
    integer :: i_solver_sw = ISolverMcICA
    integer :: i_solver_lw = ISolverMcICA

    ! Do shortwave delta-Eddington scaling on the cloud-aerosol-gas
    ! mixture (as in the original IFS scheme), rather than the more
    ! correct approach of separately scaling the cloud and aerosol
    ! scattering properties before merging with gases.  Note that
    ! .true. is not compatible with the SPARTACUS solver.
    logical :: do_sw_delta_scaling_with_gases = .false.

    ! Codes describing the gas and cloud scattering models to use, the
    ! latter of which is currently not used
    integer :: i_gas_model = IGasModelIFSRRTMG
    !     integer :: i_cloud_model

    ! Optics if i_gas_model==IGasModelMonochromatic.
    ! The wavelength to use for the Planck function in metres. If this
    ! is positive then the output longwave fluxes will be in units of
    ! W m-2 um-1.  If this is zero or negative (the default) then
    ! sigma*T^4 will be used and the output longwave fluxes will be in
    ! W m-2.
    real(jprb) :: mono_lw_wavelength = -1.0_jprb
    ! Total zenith optical depth of the atmosphere in the longwave and
    ! shortwave, distributed vertically according to the pressure.
    ! Default is zero.
    real(jprb) :: mono_lw_total_od = 0.0_jprb
    real(jprb) :: mono_sw_total_od = 0.0_jprb

    ! Codes describing particle scattering models
    integer :: i_liq_model = ILiquidModelSOCRATES
    integer :: i_ice_model = IIceModelBaran
    logical :: use_psrad_cloud_optics = .false.
    
    ! Number of separate spectral regions in which surface shortwave
    ! albedo and longwave emissivity are defined: for example, there
    ! may be just two albedos for albedo, a UV+Visible one and a
    ! near-IR one). For the moment we assume that albedo and
    ! emissivity are constant across the spectrum, so these are set
    ! equal to 1 in radiation_interface.f90
    !integer :: n_albedo_bands = 1
    !integer :: n_emiss_bands  = 1

    ! Arrays of length n_bands_lw/n_bands_sw providing an index to the
    ! surface albedo and emissivity bands
    integer, allocatable, dimension(:) :: i_albedo_from_band_sw
    integer, allocatable, dimension(:) :: i_emiss_from_band_lw

    ! Do we compute longwave and/or shortwave radiation?
    logical :: do_lw = .true.
    logical :: do_sw = .true.

    ! Do we compute clear-sky fluxes and/or solar direct fluxes?
    logical :: do_clear = .true.
    logical :: do_sw_direct = .true.

    ! Do we include 3D effects?
    logical :: do_3d_effects = .true.
    
    ! To what extent do we include "encroachment" effects in the
    ! SPARTACUS solver? This essentially means that in a situation
    ! like this
    !
    ! 000111
    ! 222222
    !
    ! radiation downwelling from region 1 may be reflected back into
    ! region 0 due to some degree of homogenization of the radiation
    ! in region 2.  Hogan and Shonk (2013) referred to this as
    ! "anomalous horizontal transport" for a 1D model, although for 3D
    ! calculations it is desirable to include at least some of it. In
    ! the longwave it is "on" or "off":
    logical :: do_3d_lw_multilayer_effects = .false.
    ! In shortwave we have more control over the process, with
    ! "minimum", "maximum" and "computed" being the options
    integer :: i_3d_sw_encroachment = IEncroachmentComputed

    ! Do we account for the effective emissivity of the side of
    ! clouds?
    logical :: do_lw_side_emissivity = .true.

    ! The 3D transfer rate "X" is such that if transport out of a
    ! region was the only process occurring then by the base of a
    ! layer only exp(-X) of the original flux would remain in that
    ! region. The transfer rate computed geometrically can be very
    ! high for the clear-sky regions in layers with high cloud
    ! fraction.  For stability reasons it is necessary to provide a
    ! maximum possible 3D transfer rate.
    real(jprb) :: max_3d_transfer_rate = 10.0_jprb

    ! By default, the Meador & Weaver (1980) expressions are used
    ! instead of the matrix exponential whenever 3D effects can be
    ! neglected (e.g. cloud-free layers or clouds with infinitely
    ! large effective cloud size), but setting the following to true
    ! uses the matrix exponential everywhere, enabling the two
    ! methods to be compared. Note that Meador & Weaver will still be
    ! used for very optically thick g points where the matrix
    ! exponential can produce incorrect results.
    logical :: use_expm_everywhere = .false.

    ! Aerosol descriptors: aerosol_type_mapping must be of length
    ! n_aerosol_types, and contains 0 if that type is to be ignored,
    ! positive numbers to map on to the indices of hydrophobic
    ! aerosols in the aerosol optics configuration file, and negative
    ! numbers to map on to (the negative of) the indices of
    ! hydrophilic aerosols in the configuration file.
    logical :: use_aerosols = .false.
    integer :: n_aerosol_types = 0
    integer :: i_aerosol_type_map(NMaxAerosolTypes)

    ! Save the gas and cloud optical properties for each g point in
    ! "radiative_properties.nc"?
    logical :: do_save_radiative_properties = .false.

    ! Save the flux profiles in each band?
    logical :: do_save_spectral_flux = .false.

    ! Save the surface downwelling shortwave fluxes in each band?
    logical :: do_surface_sw_spectral_flux = .true.

    ! Compute the longwave derivatives needed to apply the approximate
    ! radiation updates of Hogan and Bozzo (2015)
    logical :: do_lw_derivatives = .false.

    ! Save the flux profiles in each g-point (overrides
    ! do_save_spectral_flux if TRUE)?
    logical :: do_save_gpoint_flux = .false.

    ! In the IFS environment, setting up RRTM has already been done
    ! so not needed to do it again
    logical :: do_setup_ifsrrtm = .true.

    ! In the IFS environment the old scheme has a bug in the Fu
    ! longwave ice optics whereby the single scattering albedo is one
    ! minus what it should be.  Unfortunately fixing it makes
    ! forecasts worse. Setting the following to true reproduces the
    ! bug.
    logical :: do_fu_lw_ice_optics_bug = .false.

    ! Control verbosity: 0=none (no output to standard output; write
    ! to standard error only if an error occurs), 1=warning, 2=info,
    ! 3=progress, 4=detailed, 5=debug.  Separate settings for the
    ! setup of the scheme and the execution of it.
    integer :: iverbosesetup = 3
    integer :: iverbose = 1

    ! Radiative transfer in complex surface canopies
    ! (streets/vegetation): do we use the full spectrum as in the
    ! atmosphere, or just the reduced spectrum in which the shortwave
    ! albedo and longwave emissivity are provided?
    logical :: use_canopy_full_spectrum_sw = .false.
    logical :: use_canopy_full_spectrum_lw = .false.
    ! Do we treat gas radiative transfer in streets/vegetation?
    logical :: do_canopy_gases_sw = .false.
    logical :: do_canopy_gases_lw = .false.

    ! Optics file names for overriding the ones generated from the
    ! other options. If these remain empty then the generated names
    ! will be used (see the "consolidate_config" routine below). If
    ! the user assigns one of these and it starts with a '/' character
    ! then that will be used instead. If the user assigns one and it
    ! doesn't start with a '/' character then it will be prepended by
    ! the contents of directory_name.
    character(len=511) :: ice_optics_override_file_name = ''
    character(len=511) :: liq_optics_override_file_name = ''
    character(len=511) :: aerosol_optics_override_file_name = ''

    ! Optionally override the look-up table file for the cloud-water
    ! PDF used by the McICA solver
    character(len=511) :: cloud_pdf_override_file_name = ''

    ! COMPUTED PARAMETERS
    ! Users of this library should not edit these parameters directly

    ! Wavenumber range for each band, in cm-1, which will be allocated
    ! to be of length n_bands_sw or n_bands_lw
    real(jprb), allocatable, dimension(:) :: wavenumber1_sw
    real(jprb), allocatable, dimension(:) :: wavenumber2_sw
    real(jprb), allocatable, dimension(:) :: wavenumber1_lw
    real(jprb), allocatable, dimension(:) :: wavenumber2_lw

    ! Arrays of length the number of g-points that convert from
    ! g-point to the band index
    integer, allocatable, dimension(:) :: i_band_from_g_lw
    integer, allocatable, dimension(:) :: i_band_from_g_sw

    ! We allow for the possibility for g-points to be ordered in terms
    ! of likely absorption (weakest to strongest) across the shortwave
    ! or longwave spectrum, in order that in SPARTACUS we select only
    ! the first n g-points that will not have too large an absorption,
    ! and therefore matrix exponentials that are both finite and not
    ! too expensive to compute.  The following two arrays map the
    ! reordered g-points to the original ones.
    integer, allocatable, dimension(:) :: i_g_from_reordered_g_lw
    integer, allocatable, dimension(:) :: i_g_from_reordered_g_sw

    ! The following map the reordered g-points to the bands
    integer, allocatable, dimension(:) :: i_band_from_reordered_g_lw
    integer, allocatable, dimension(:) :: i_band_from_reordered_g_sw

    ! The following map the reordered g-points to the spectral
    ! information being saved: if do_save_gpoint_flux==TRUE then this
    ! will map on to the original g points, but if only
    ! do_save_spectral_flux==TRUE then this will map on to the bands
    integer, pointer, dimension(:) :: i_spec_from_reordered_g_lw
    integer, pointer, dimension(:) :: i_spec_from_reordered_g_sw

    ! Number of spectral intervals used for the canopy radiative
    ! transfer calculation; they are either equal to
    ! n_albedo_bands/n_emiss_bands or n_g_sw/n_g_lw
    integer :: n_canopy_bands_sw = 1
    integer :: n_canopy_bands_lw = 1

    ! Data structure containing cloud scattering data
    type(cloud_optics_type)      :: cloud_optics

    ! Data structure containing aerosol scattering data
    type(aerosol_optics_type)    :: aerosol_optics

    ! Object for sampling from a gamma or lognormal distribution
    type(pdf_sampler_type)       :: pdf_sampler

    ! Optics file names
    character(len=511) :: ice_optics_file_name, &
         &                liq_optics_file_name, &
         &                aerosol_optics_file_name
    
    ! McICA PDF look-up table file name
    character(len=511) :: cloud_pdf_file_name

    ! Number of gpoints and bands in the shortwave and longwave - set
    ! to zero as will be set properly later
    integer :: n_g_sw = 0, n_g_lw = 0
    integer :: n_bands_sw = 0, n_bands_lw = 0

    ! Number of spectral points to save (equal either to the number of
    ! g points or the number of bands
    integer :: n_spec_sw = 0, n_spec_lw = 0

    ! Dimensions to store variables that are only needed if longwave
    ! scattering is included. "n_g_lw_if_scattering" is equal to
    ! "n_g_lw" if aerosols are allowed to scatter in the longwave,
    ! and zero otherwise. "n_bands_lw_if_scattering" is equal to
    ! "n_bands_lw" if clouds are allowed to scatter in the longwave,
    ! and zero otherwise.
    integer :: n_g_lw_if_scattering = 0, n_bands_lw_if_scattering = 0

    ! Treat clouds as horizontally homogeneous within the gribox
    logical :: is_homogeneous = .false.

   contains
     procedure :: read => read_config_from_namelist
     procedure :: consolidate => consolidate_config
     procedure :: set  => set_config
     procedure :: print => print_config
     procedure :: get_sw_weights
     procedure :: define_sw_albedo_intervals
     procedure :: define_lw_emiss_intervals

  end type config_type

!  procedure, private :: print_logical, print_real, print_int

contains


  !---------------------------------------------------------------------
  ! This subroutine reads configuration data from a namelist file, and
  ! anything that is not in the namelists will be set to default
  ! values. If optional output argument "is_success" is present, then
  ! on error (e.g. missing file) it will be set to .false.; if this
  ! argument is missing then on error the program will be aborted. You
  ! may either specify the file_name or the unit of an open file to
  ! read, but not both.
  subroutine read_config_from_namelist(this, file_name, unit, is_success)

    use ecradhook,      only : lhook, dr_hook
    use radiation_io, only : nulerr, nulrad, radiation_abort

    class(config_type), intent(inout)         :: this
    character(*),       intent(in),  optional :: file_name
    integer,            intent(in),  optional :: unit
    logical,            intent(out), optional :: is_success

    integer :: iosopen, iosread ! Status after calling open and read

    ! The following variables are read from the namelists and map
    ! directly onto members of the config_type derived type

    ! To be read from the radiation_config namelist 
    logical :: do_sw, do_lw, do_clear, do_sw_direct
    logical :: do_3d_effects, use_expm_everywhere, use_aerosols
    logical :: do_lw_side_emissivity
    logical :: do_3d_lw_multilayer_effects, do_fu_lw_ice_optics_bug
    logical :: do_lw_aerosol_scattering, do_lw_cloud_scattering
    logical :: do_save_radiative_properties, do_save_spectral_flux
    logical :: do_save_gpoint_flux, do_surface_sw_spectral_flux
    logical :: use_beta_overlap, do_lw_derivatives
    logical :: do_sw_delta_scaling_with_gases
    logical :: use_canopy_full_spectrum_sw, use_canopy_full_spectrum_lw
    logical :: do_canopy_gases_sw, do_canopy_gases_lw
    integer :: n_regions, iverbose, iverbosesetup, n_aerosol_types
    real(jprb):: mono_lw_wavelength, mono_lw_total_od, mono_sw_total_od
    real(jprb):: cloud_inhom_decorr_scaling, cloud_fraction_threshold
    real(jprb):: clear_to_thick_fraction, max_gas_od_3d, max_cloud_od
    real(jprb):: cloud_mixing_ratio_threshold, overhead_sun_factor
    real(jprb):: max_3d_transfer_rate
    character(511) :: directory_name
    character(511) :: cloud_pdf_override_file_name
    character(63)  :: liquid_model_name, ice_model_name, gas_model_name
    character(63)  :: sw_solver_name, lw_solver_name, overlap_scheme_name
    character(63)  :: sw_encroachment_name, cloud_pdf_shape_name
    integer :: i_aerosol_type_map(NMaxAerosolTypes) ! More than 256 is an error

    integer :: iunit ! Unit number of namelist file

    namelist /radiation/ do_sw, do_lw, do_sw_direct, &
         &  do_3d_effects, do_lw_side_emissivity, do_clear, &
         &  do_save_radiative_properties, sw_encroachment_name, &
         &  do_3d_lw_multilayer_effects, do_fu_lw_ice_optics_bug, &
         &  do_save_spectral_flux, do_save_gpoint_flux, &
         &  do_surface_sw_spectral_flux, do_lw_derivatives, &
         &  do_lw_aerosol_scattering, do_lw_cloud_scattering, &
         &  n_regions, directory_name, gas_model_name, cloud_pdf_override_file_name, &
         &  liquid_model_name, ice_model_name, max_3d_transfer_rate, &
         &  use_canopy_full_spectrum_sw, use_canopy_full_spectrum_lw, &
         &  do_canopy_gases_sw, do_canopy_gases_lw, &
         &  do_sw_delta_scaling_with_gases, overlap_scheme_name, &
         &  sw_solver_name, lw_solver_name, use_beta_overlap, &
         &  use_expm_everywhere, iverbose, iverbosesetup, &
         &  cloud_inhom_decorr_scaling, cloud_fraction_threshold, &
         &  clear_to_thick_fraction, max_gas_od_3d, max_cloud_od, &
         &  cloud_mixing_ratio_threshold, overhead_sun_factor, &
         &  n_aerosol_types, i_aerosol_type_map, use_aerosols, &
         &  mono_lw_wavelength, mono_lw_total_od, mono_sw_total_od, &
         &  cloud_pdf_shape_name

    real(jprb) :: hook_handle

    if (lhook) call dr_hook('radiation_config:read',0,hook_handle)

    ! Copy default values from the original structure 
    do_sw = this%do_sw
    do_lw = this%do_lw
    do_sw_direct = this%do_sw_direct
    do_3d_effects = this%do_3d_effects
    do_3d_lw_multilayer_effects = this%do_3d_lw_multilayer_effects
    do_lw_side_emissivity = this%do_lw_side_emissivity
    do_clear = this%do_clear
    do_lw_aerosol_scattering = this%do_lw_aerosol_scattering
    do_lw_cloud_scattering = this%do_lw_cloud_scattering
    do_sw_delta_scaling_with_gases = this%do_sw_delta_scaling_with_gases
    do_fu_lw_ice_optics_bug = this%do_fu_lw_ice_optics_bug
    use_canopy_full_spectrum_sw = this%use_canopy_full_spectrum_sw
    use_canopy_full_spectrum_lw = this%use_canopy_full_spectrum_lw
    do_canopy_gases_sw = this%do_canopy_gases_sw
    do_canopy_gases_lw = this%do_canopy_gases_lw
    n_regions = this%nregions
    directory_name = this%directory_name
    cloud_pdf_override_file_name = this%cloud_pdf_override_file_name
    use_expm_everywhere = this%use_expm_everywhere
    use_aerosols = this%use_aerosols
    do_save_radiative_properties = this%do_save_radiative_properties
    do_save_spectral_flux = this%do_save_spectral_flux
    do_save_gpoint_flux = this%do_save_gpoint_flux
    do_lw_derivatives = this%do_lw_derivatives
    do_surface_sw_spectral_flux = this%do_surface_sw_spectral_flux
    iverbose = this%iverbose
    iverbosesetup = this%iverbosesetup
    cloud_fraction_threshold = this%cloud_fraction_threshold
    cloud_mixing_ratio_threshold = this%cloud_mixing_ratio_threshold
    use_beta_overlap = this%use_beta_overlap
    cloud_inhom_decorr_scaling = this%cloud_inhom_decorr_scaling
    clear_to_thick_fraction = this%clear_to_thick_fraction
    overhead_sun_factor = this%overhead_sun_factor
    max_gas_od_3d = this%max_gas_od_3d
    max_cloud_od = this%max_cloud_od
    max_3d_transfer_rate = this%max_3d_transfer_rate
    gas_model_name = '' !DefaultGasModelName
    liquid_model_name = '' !DefaultLiquidModelName
    ice_model_name = '' !DefaultIceModelName
    sw_solver_name = '' !DefaultSwSolverName
    lw_solver_name = '' !DefaultLwSolverName
    sw_encroachment_name = ''
    overlap_scheme_name = ''
    cloud_pdf_shape_name = ''
    n_aerosol_types = this%n_aerosol_types
    mono_lw_wavelength = this%mono_lw_wavelength
    mono_lw_total_od = this%mono_lw_total_od
    mono_sw_total_od = this%mono_sw_total_od
    i_aerosol_type_map = this%i_aerosol_type_map

    if (present(file_name) .and. present(unit)) then
      write(nulerr,'(a,a,a)') '*** Error: cannot specify both file_name and unit in call to config_type%read'
      call radiation_abort('Radiation configuration error')
    else if (.not. present(file_name) .and. .not. present(unit)) then
      write(nulerr,'(a,a,a)') '*** Error: neither file_name nor unit specified in call to config_type%read'
      call radiation_abort('Radiation configuration error')
    end if

    if (present(file_name)) then
      ! Open the namelist file
      iunit = nulrad
      open(unit=iunit, iostat=iosopen, file=trim(file_name))
    else
      ! Assume that iunit represents and open file
      iosopen = 0
      iunit = unit
    end if

    if (iosopen /= 0) then
      ! An error occurred opening the file
      if (present(is_success)) then
        is_success = .false.
        ! We now continue the subroutine so that the default values
        ! are placed in the config structure
      else
        write(nulerr,'(a,a,a)') '*** Error: namelist file "', &
             &                trim(file_name), '" not found'
        call radiation_abort('Radiation configuration error')
      end if
    else
      read(unit=iunit, iostat=iosread, nml=radiation)
      if (iosread /= 0) then
        ! An error occurred reading the file
        if (present(is_success)) then
          is_success = .false.
          ! We now continue the subroutine so that the default values
          ! are placed in the config structure
        else if (present(file_name)) then
          write(nulerr,'(a,a,a)') '*** Error reading namelist "radiation" from file "', &
               &      trim(file_name), '"'
          close(unit=iunit)
          call radiation_abort('Radiation configuration error')
        else
          write(nulerr,'(a,i0)') '*** Error reading namelist "radiation" from unit ', &
               &      iunit
          call radiation_abort('Radiation configuration error')
        end if
      end if

      if (present(file_name)) then
        close(unit=iunit)
      end if
    end if

    ! Copy namelist data into configuration object
    this%do_lw = do_lw
    this%do_sw = do_sw
    this%do_clear = do_clear
    this%do_sw_direct = do_sw_direct
    this%do_3d_effects = do_3d_effects
    this%do_3d_lw_multilayer_effects = do_3d_lw_multilayer_effects
    this%do_lw_side_emissivity = do_lw_side_emissivity
    this%use_expm_everywhere = use_expm_everywhere
    this%use_aerosols = use_aerosols
    this%do_lw_cloud_scattering = do_lw_cloud_scattering
    this%do_lw_aerosol_scattering = do_lw_aerosol_scattering
    this%nregions = n_regions
    this%do_surface_sw_spectral_flux = do_surface_sw_spectral_flux
    this%do_sw_delta_scaling_with_gases = do_sw_delta_scaling_with_gases
    this%do_fu_lw_ice_optics_bug = do_fu_lw_ice_optics_bug
    this%use_canopy_full_spectrum_sw = use_canopy_full_spectrum_sw
    this%use_canopy_full_spectrum_lw = use_canopy_full_spectrum_lw
    this%do_canopy_gases_sw = do_canopy_gases_sw
    this%do_canopy_gases_lw = do_canopy_gases_lw
    this%mono_lw_wavelength = mono_lw_wavelength
    this%mono_lw_total_od = mono_lw_total_od
    this%mono_sw_total_od = mono_sw_total_od
    this%use_beta_overlap = use_beta_overlap
    this%cloud_inhom_decorr_scaling = cloud_inhom_decorr_scaling
    this%clear_to_thick_fraction = clear_to_thick_fraction
    this%overhead_sun_factor = overhead_sun_factor
    this%max_gas_od_3d = max_gas_od_3d
    this%max_cloud_od = max_cloud_od
    this%max_3d_transfer_rate = max_3d_transfer_rate
    this%directory_name = directory_name
    this%cloud_pdf_override_file_name = cloud_pdf_override_file_name
    this%cloud_fraction_threshold = cloud_fraction_threshold
    this%cloud_mixing_ratio_threshold = cloud_mixing_ratio_threshold
    this%n_aerosol_types = n_aerosol_types
    this%do_save_radiative_properties = do_save_radiative_properties
    this%do_lw_derivatives = do_lw_derivatives
    this%do_save_spectral_flux = do_save_spectral_flux
    this%do_save_gpoint_flux = do_save_gpoint_flux
    if (do_save_gpoint_flux) then
      ! Saving the fluxes every g-point overrides saving as averaged
      ! in a band, but this%do_save_spectral_flux needs to be TRUE as
      ! it is tested inside the solver routines to decide whether to
      ! save anything
      this%do_save_spectral_flux = .true.
    end if

    ! Determine liquid optics model
    call get_enum_code(liquid_model_name, LiquidModelName, &
         &            'liquid_model_name', this%i_liq_model)

    ! Determine ice optics model
    call get_enum_code(ice_model_name, IceModelName, &
         &            'ice_model_name', this%i_ice_model)

    ! Determine gas optics model
    call get_enum_code(gas_model_name, GasModelName, &
         &            'gas_model_name', this%i_gas_model)

    ! Determine solvers
    call get_enum_code(sw_solver_name, SolverName, &
         &            'sw_solver_name', this%i_solver_sw)
    call get_enum_code(lw_solver_name, SolverName, &
         &            'lw_solver_name', this%i_solver_lw)

    call get_enum_code(sw_encroachment_name, EncroachmentName, &
         &            'sw_encroachment_name', this%i_3d_sw_encroachment)
    ! Determine overlap scheme
    call get_enum_code(overlap_scheme_name, OverlapName, &
         &             'overlap_scheme_name', this%i_overlap_scheme)
    
    ! Determine cloud PDF shape 
    call get_enum_code(cloud_pdf_shape_name, PdfShapeName, &
         &             'cloud_pdf_shape_name', this%i_cloud_pdf_shape)

    this%i_aerosol_type_map = 0
    if (this%use_aerosols) then
      this%i_aerosol_type_map(1:n_aerosol_types) &
           &  = i_aerosol_type_map(1:n_aerosol_types)
    end if

    ! Ensure verbose levels within limits
    if (iverbose < 0) then
      iverbose = 0
    end if
    this%iverbose = iverbose

    if (iverbosesetup < 0) then
      iverbosesetup = 0
    end if
    this%iverbosesetup = iverbosesetup

    ! Normal subroutine exit
    if (present(is_success)) then
      is_success = .true.
    end if

    if (lhook) call dr_hook('radiation_config:read',1,hook_handle)

  end subroutine read_config_from_namelist


  !---------------------------------------------------------------------
  ! This routine is called by radiation_interface:setup_radiation and
  ! it converts the user specified options into some more specific
  ! data such as data file names
  subroutine consolidate_config(this)

    use ecradhook,      only : lhook, dr_hook
    use radiation_io, only : nulout, nulerr, radiation_abort

    class(config_type), intent(inout)         :: this

    real(jprb) :: hook_handle

    if (lhook) call dr_hook('radiation_config:consolidate',0,hook_handle)

    ! Check consistency of models
    if ((this%i_liq_model == ILiquidModelHuStamnesPSRAD &
         .and. this%i_ice_model /= IIceModelFuPSRAD) .or. &
         (this%i_liq_model /= ILiquidModelHuStamnesPSRAD &
         .and. this%i_ice_model == IIceModelFuPSRAD)) then
      write(nulerr,'(a,a,a)') '*** Error: liquid model "HuStamnesPSRAD" must be used with ice model "FuPSRAD"'
      call radiation_abort('Radiation configuration error')
    else if (this%i_liq_model == ILiquidModelHuStamnesPSRAD) then
      this%use_psrad_cloud_optics = .true.
    else
      this%use_psrad_cloud_optics = .false.
    end if

    ! SPARTACUS only works with Exp-Ran overlap scheme
    if ((       this%i_solver_sw == ISolverSPARTACUS &
         & .or. this%i_solver_lw == ISolverSPARTACUS &
         & .or. this%i_solver_sw == ISolverTripleclouds &
         & .or. this%i_solver_lw == ISolverTripleclouds) &
         & .and. this%i_overlap_scheme /= IOverlapExponentialRandom) then
      write(nulerr,'(a)') '*** Error: SPARTACUS/Tripleclouds solvers can only do Exponential-Random overlap'
      call radiation_abort('Radiation configuration error')

    end if

    ! Set aerosol optics file name
    if (len_trim(this%aerosol_optics_override_file_name) > 0) then
      if (this%aerosol_optics_override_file_name(1:1) == '/') then
        this%aerosol_optics_file_name = trim(this%aerosol_optics_override_file_name)
      else
        this%aerosol_optics_file_name = trim(this%directory_name) &
             &  // '/' // trim(this%aerosol_optics_override_file_name)
      end if
    else
      this%aerosol_optics_file_name &
           &   = trim(this%directory_name) // "/aerosol_ifs_rrtm_43R3.nc"
    end if

    ! Set liquid optics file name
    if (len_trim(this%liq_optics_override_file_name) > 0) then
      if (this%liq_optics_override_file_name(1:1) == '/') then
        this%liq_optics_file_name = trim(this%liq_optics_override_file_name)
      else
        this%liq_optics_file_name = trim(this%directory_name) &
             &  // '/' // trim(this%liq_optics_override_file_name)
      end if
    else if (this%i_liq_model == ILiquidModelSOCRATES) then
      this%liq_optics_file_name &
           &  = trim(this%directory_name) // "/socrates_droplet_scattering_rrtm.nc"
    else if (this%i_liq_model == ILiquidModelSlingo) then
      this%liq_optics_file_name &
           &  = trim(this%directory_name) // "/slingo_droplet_scattering_rrtm.nc"
    end if

    ! Set ice optics file name
    if (len_trim(this%ice_optics_override_file_name) > 0) then
      if (this%ice_optics_override_file_name(1:1) == '/') then
        this%ice_optics_file_name = trim(this%ice_optics_override_file_name)
      else
        this%ice_optics_file_name = trim(this%directory_name) &
             &  // '/' // trim(this%ice_optics_override_file_name)
      end if
    else if (this%i_ice_model == IIceModelFu) then
      this%ice_optics_file_name &
           &   = trim(this%directory_name) // "/fu_ice_scattering_rrtm.nc"
    else if (this%i_ice_model == IIceModelBaran) then
      this%ice_optics_file_name &
           &   = trim(this%directory_name) // "/baran_ice_scattering_rrtm.nc"
    else if (this%i_ice_model == IIceModelBaran2016) then
      this%ice_optics_file_name &
           &   = trim(this%directory_name) // "/baran2016_ice_scattering_rrtm.nc"
    else if (this%i_ice_model == IIceModelBaran2017) then
      this%ice_optics_file_name &
           &   = trim(this%directory_name) // "/baran2017_ice_scattering_rrtm.nc"
    else if (this%i_ice_model == IIceModelYi) then
      this%ice_optics_file_name &
           &   = trim(this%directory_name) // "/yi_ice_scattering_rrtm.nc"
    end if

    ! Set cloud-water PDF look-up table file name
    if (len_trim(this%cloud_pdf_override_file_name) > 0) then
      if (this%cloud_pdf_override_file_name(1:1) == '/') then
        this%cloud_pdf_file_name = trim(this%cloud_pdf_override_file_name)
      else
        this%cloud_pdf_file_name = trim(this%directory_name) &
             &  // '/' // trim(this%cloud_pdf_override_file_name)
      end if
    elseif (this%i_cloud_pdf_shape == IPdfShapeLognormal) then
      this%cloud_pdf_file_name = trim(this%directory_name) // "/mcica_lognormal.nc"
    else
      this%cloud_pdf_file_name = trim(this%directory_name) // "/mcica_gamma.nc"
    end if

    ! Aerosol data
    if (this%n_aerosol_types < 0 &
         &  .or. this%n_aerosol_types > NMaxAerosolTypes) then
      write(nulerr,'(a,i0)') '*** Error: number of aerosol types must be between 0 and ', &
           &  NMaxAerosolTypes
      call radiation_abort('Radiation configuration error')
    end if

!    if (this%use_aerosols .and. this%n_aerosol_types == 0) then
!      this%use_aerosols = .false.
!      if (this%iverbose >= 1) then
!        write(nulout, '(a)') 'Warning: no aerosol types specified so aerosols switched off'
!      end if
!    end if

    ! In the monochromatic case we need to override the liquid, ice
    ! and aerosol models to ensure compatibility
    if (this%i_gas_model == IGasModelMonochromatic) then
      this%use_psrad_cloud_optics = .false.
      this%i_liq_model = ILiquidModelMonochromatic
      this%i_ice_model = IIceModelMonochromatic
      this%use_aerosols = .false.
    end if

    ! McICA solver currently can't store full profiles of spectral fluxes
    if (this%i_solver_sw == ISolverMcICA) then
      this%do_save_spectral_flux = .false.
    end if

    if (this%i_solver_sw == ISolverSPARTACUS .and. this%do_sw_delta_scaling_with_gases) then
      write(nulerr,'(a)') '*** Error: SW delta-Eddington scaling with gases not possible with SPARTACUS solver'
      call radiation_abort('Radiation configuration error')
    end if

    if ((this%do_lw .and. this%do_sw) .and. &
         & (     (      this%i_solver_sw == ISolverHomogeneous  &
         &        .and. this%i_solver_lw /= ISolverHomogeneous) &
         &  .or. (      this%i_solver_sw /= ISolverHomogeneous  &
         &        .and. this%i_solver_lw == ISolverHomogeneous) &
         & ) ) then
      write(nulerr,'(a)') '*** Error: if one solver is "Homogeneous" then the other must be'
      call radiation_abort('Radiation configuration error')
    end if

    ! Set is_homogeneous if the active solvers are homogeneous, since
    ! this affects how "in-cloud" water contents are computed
    if (        (this%do_sw .and. this%i_solver_sw == ISolverHomogeneous) &
         & .or. (this%do_lw .and. this%i_solver_lw == ISolverHomogeneous)) then
      this%is_homogeneous = .true.
    end if

    if (lhook) call dr_hook('radiation_config:consolidate',1,hook_handle)

  end subroutine consolidate_config


  !---------------------------------------------------------------------
  ! This subroutine sets members of the configuration object via
  ! optional arguments, and any member not specified is left
  ! untouched. Therefore, this should be called after taking data from
  ! the namelist.
  subroutine set_config(config, directory_name, &
       &  do_lw, do_sw, &
       &  do_lw_aerosol_scattering, do_lw_cloud_scattering, &
       &  do_sw_direct)

    class(config_type), intent(inout):: config
    character(len=*), intent(in), optional  :: directory_name
    logical, intent(in), optional           :: do_lw, do_sw
    logical, intent(in), optional           :: do_lw_aerosol_scattering
    logical, intent(in), optional           :: do_lw_cloud_scattering
    logical, intent(in), optional           :: do_sw_direct

    if (present(do_lw)) then
       config%do_lw = do_lw
    end if

    if(present(do_sw)) then
       config%do_sw = do_sw
    end if

    if (present(do_sw_direct)) then
       config%do_sw_direct = do_sw_direct
    end if

    if (present(directory_name)) then
       config%directory_name = trim(directory_name)
    end if

    if (present(do_lw_aerosol_scattering)) then
       config%do_lw_aerosol_scattering = .true.
    end if

    if (present(do_lw_cloud_scattering)) then
       config%do_lw_cloud_scattering = .true.
    end if

  end subroutine set_config


  !---------------------------------------------------------------------
  ! Print configuration information to standard output
  subroutine print_config(this, iverbose)

    use radiation_io, only : nulout

    class(config_type), intent(in) :: this

    integer, optional,  intent(in) :: iverbose
    integer                        :: i_local_verbose

    if (present(iverbose)) then
      i_local_verbose = iverbose
    else
      i_local_verbose = this%iverbose
    end if

    if (i_local_verbose >= 2) then
      !---------------------------------------------------------------------
      write(nulout, '(a)') 'General settings:'
      write(nulout, '(a,a,a)') '  Data files expected in "', &
           &                   trim(this%directory_name), '"'
      call print_logical('  Clear-sky calculations are', 'do_clear', this%do_clear)
      call print_logical('  Saving intermediate radiative properties', &
           &   'do_save_radiative_properties', this%do_save_radiative_properties)
      call print_logical('  Saving spectral flux profiles', &
           &   'do_save_spectral_flux', this%do_save_spectral_flux)
      if (this%do_sw) then
        call print_logical('  Saving surface shortwave spectral fluxes', &
             &   'do_surface_sw_spectral_flux', this%do_surface_sw_spectral_flux)
      end if
      call print_enum('  Gas model is', GasModelName, 'i_gas_model', &
           &          this%i_gas_model)
      call print_logical('  Aerosols are', 'use_aerosols', this%use_aerosols)
      call print_logical('  Longwave derivative calculation is', &
           &   'do_lw_derivatives',this%do_lw_derivatives)


      !---------------------------------------------------------------------
      write(nulout, '(a)') 'Cloud settings:'
      call print_real('  Cloud fraction threshold', &
           &   'cloud_fraction_threshold', this%cloud_fraction_threshold)
      call print_real('  Cloud mixing-ratio threshold', &
           &   'cloud_mixing_ratio_threshold', this%cloud_mixing_ratio_threshold)
      call print_enum('  Liquid optics scheme is', LiquidModelName, &
           &          'i_liq_model',this%i_liq_model)
      call print_enum('  Ice optics scheme is', IceModelName, &
           &          'i_ice_model',this%i_ice_model)
      if (this%i_ice_model == IIceModelFu) then
        call print_logical('  Longwave ice optics bug in Fu scheme is', &
             &   'do_fu_lw_ice_optics_bug',this%do_fu_lw_ice_optics_bug)
      end if
      call print_enum('  Cloud overlap scheme is', OverlapName, &
           &          'i_overlap_scheme',this%i_overlap_scheme)
      call print_logical('  Use "beta" overlap parameter is', &
           &   'use_beta_overlap', this%use_beta_overlap)
      call print_enum('  Cloud PDF shape is', PdfShapeName, &
           &          'i_cloud_pdf_shape',this%i_cloud_pdf_shape)
      call print_real('  Cloud inhom decorrelation scaling', &
           &   'cloud_inhom_decorr_scaling', this%cloud_inhom_decorr_scaling)

      !---------------------------------------------------------------------
      write(nulout, '(a)') 'Solver settings:'
      if (this%do_sw) then
        call print_enum('  Shortwave solver is', SolverName, &
             &          'i_solver_sw', this%i_solver_sw)
        
        if (this%i_gas_model == IGasModelMonochromatic) then
          call print_real('  Shortwave atmospheric optical depth', &
               &   'mono_sw_total_od', this%mono_sw_total_od)  
        end if
        call print_logical('  Shortwave delta scaling after merge with gases', &
             &   'do_sw_delta_scaling_with_gases', &
             &   this%do_sw_delta_scaling_with_gases)
      else
        call print_logical('  Shortwave calculations are','do_sw',this%do_sw)
      end if

      if (this%do_lw) then
        call print_enum('  Longwave solver is', SolverName, 'i_solver_lw', &
             &          this%i_solver_lw)

        if (this%i_gas_model == IGasModelMonochromatic) then
          if (this%mono_lw_wavelength > 0.0_jprb) then
            call print_real('  Longwave effective wavelength (m)', &
                 &   'mono_lw_wavelength', this%mono_lw_wavelength)
          else
            write(nulout,'(a)') '  Longwave fluxes are broadband                              (mono_lw_wavelength<=0)'
          end if
          call print_real('  Longwave atmospheric optical depth', &
               &   'mono_lw_total_od', this%mono_lw_total_od)  
        end if
        call print_logical('  Longwave cloud scattering is', &
             &   'do_lw_cloud_scattering',this%do_lw_cloud_scattering)
        call print_logical('  Longwave aerosol scattering is', &
             &   'do_lw_aerosol_scattering',this%do_lw_aerosol_scattering)
      else
        call print_logical('  Longwave calculations are','do_lw', this%do_lw)
      end if

      if (this%i_solver_sw == ISolverSpartacus &
           &  .or. this%i_solver_lw == ISolverSpartacus) then
        write(nulout, '(a)') '  SPARTACUS options:'
        call print_integer('    Number of regions', 'nregions', this%nregions)
        call print_real('    Max cloud optical depth per layer', &
             &   'max_cloud_od', this%max_cloud_od)
        call print_enum('    Shortwave encroachment is', EncroachmentName, &
             &          'i_3d_sw_encroachment', this%i_3d_sw_encroachment)
        call print_logical('    Multilayer longwave horizontal transport is', &
             'do_3d_lw_multilayer_effects', this%do_3d_lw_multilayer_effects)
        call print_logical('    Use matrix exponential everywhere is', &
             &   'use_expm_everywhere', this%use_expm_everywhere)
        call print_logical('    3D effects are', 'do_3d_effects', &
             &             this%do_3d_effects)

        if (this%do_3d_effects) then
          call print_logical('    Longwave side emissivity parameterization is', &
               &  'do_lw_side_emissivity', this%do_lw_side_emissivity)
          call print_real('    Clear-to-thick edge fraction is', &
               &  'clear_to_thick_fraction', this%clear_to_thick_fraction)
          call print_real('    Overhead sun factor is', &
               &  'overhead_sun_factor', this%overhead_sun_factor)
          call print_real('    Max gas optical depth for 3D effects', &
               &   'max_gas_od_3d', this%max_gas_od_3d)
          call print_real('    Max 3D transfer rate', &
               &   'max_3d_transfer_rate', this%max_3d_transfer_rate)
        end if
      end if
            
    end if
    
  end subroutine print_config



  !---------------------------------------------------------------------
  ! In order to estimate UV and photosynthetically active radiation,
  ! we need weighted sum of fluxes considering wavelength range
  ! required.  This routine returns information for how to correctly
  ! weight output spectral fluxes for a range of input wavelengths.
  ! Note that this is approximate; internally it may be assumed that
  ! the energy is uniformly distributed in wavenumber space, for
  ! example.  If the character string "weighting_name" is present, and
  ! iverbose>=2, then information on the weighting will be provided on
  ! nulout.
  subroutine get_sw_weights(this, wavelength1, wavelength2, &
       &                    nweights, iband, weight, weighting_name)

    use parkind1, only : jprb
    use radiation_io, only : nulout, nulerr, radiation_abort

    class(config_type), intent(in) :: this
    ! Range of wavelengths to get weights for (m)
    real(jprb), intent(in) :: wavelength1, wavelength2
    ! Output number of weights needed
    integer,    intent(out)   :: nweights
    ! Only write to the first nweights of these arrays: they contain
    ! the indices to the non-zero bands, and the weight in each of
    ! those bands
    integer,    intent(out)   :: iband(:)
    real(jprb), intent(out)   :: weight(:)
    character(len=*), optional, intent(in) :: weighting_name

    ! Internally we deal with wavenumber
    real(jprb) :: wavenumber1, wavenumber2 ! cm-1

    integer :: jband ! Loop index for spectral band

    if (this%n_bands_sw <= 0) then
      write(nulerr,'(a)') '*** Error: get_sw_weights called before number of shortwave bands set'
      call radiation_abort()      
    end if

    ! Convert wavelength range (m) to wavenumber (cm-1)
    wavenumber1 = 0.01_jprb / wavelength2
    wavenumber2 = 0.01_jprb / wavelength1

    nweights = 0

    do jband = 1,this%n_bands_sw
      if (wavenumber1 < this%wavenumber2_sw(jband) &
           &  .and. wavenumber2 > this%wavenumber1_sw(jband)) then
        nweights = nweights+1
        iband(nweights) = jband
        weight(nweights) = (min(wavenumber2,this%wavenumber2_sw(jband)) &
             &         - max(wavenumber1,this%wavenumber1_sw(jband))) &
             & / (this%wavenumber2_sw(jband) - this%wavenumber1_sw(jband))
      end if
    end do

    if (nweights == 0) then
      write(nulerr,'(a,e8.4,a,e8.4,a)') '*** Error: wavelength range ', &
           &  wavelength1, ' to ', wavelength2, ' m is outside shortwave band'
      call radiation_abort()
    else if (this%iverbose >= 2 .and. present(weighting_name)) then
      write(nulout,'(a,a,a,f6.0,a,f6.0,a)') 'Spectral weights for ', &
           &  weighting_name, ' (', wavenumber1, ' to ', &
           &  wavenumber2, ' cm-1):'
      do jband = 1, nweights
        write(nulout, '(a,i0,a,f6.0,a,f6.0,a,f8.4)') '  Shortwave band ', &
             &  iband(jband), ' (', this%wavenumber1_sw(iband(jband)), ' to ', &
             &  this%wavenumber2_sw(iband(jband)), ' cm-1): ', weight(jband)
      end do
    end if

  end subroutine get_sw_weights


  !---------------------------------------------------------------------
  ! The input shortwave surface albedo coming in is likely to be in
  ! different spectral intervals to the gas model in the radiation
  ! scheme. We assume that the input albedo is defined within
  ! "ninterval" spectral intervals covering the wavelength range 0 to
  ! infinity, but allow for the possibility that two intervals may be
  ! indexed back to the same albedo band.  
  subroutine define_sw_albedo_intervals(this, ninterval, wavelength_bound, &
       &                                i_band_in)

    use radiation_io, only : nulout, nulerr, radiation_abort

    class(config_type),   intent(inout) :: this
    ! Number of spectral intervals in which albedo is defined
    integer,              intent(in)    :: ninterval
    ! Monotonically increasing wavelength bounds between intervals,
    ! not including the outer bounds (which are assumed to be zero and
    ! infinity)
    real(jprb), optional, intent(in)    :: wavelength_bound(ninterval-1)
    ! The albedo band indices corresponding to each interval
    integer,    optional, intent(in)    :: i_band_in(ninterval)

    integer :: jband
    integer :: iinterval
    real(jprb) :: wavenumber_bound ! cm-1
    real(jprb) :: wavenumber_mid   ! cm-1

    if (.not. allocated(this%i_albedo_from_band_sw)) then
      allocate(this%i_albedo_from_band_sw(this%n_bands_sw))
    end if

    if (ninterval == 1) then 
      ! Only one interval, so wavelength bounds are not needed
      if (present(i_band_in)) then
        this%i_albedo_from_band_sw = i_band_in(1)
      else
        ! If there is only one albedo interval then the index to it
        ! should always be 1
        this%i_albedo_from_band_sw = 1
      end if
    else if (.not. present(i_band_in) .or. .not. present(wavelength_bound)) then
      write(nulerr, '(a,a)') '*** Error: define_sw_albedo_intervals called', &
           ' with ninterval>1, but wavelength_bound or i_band_in not provided'
      call radiation_abort()
    else
      ! Check wavelength is monotonically increasing
      do jband = 2,ninterval-1
        if (wavelength_bound(jband) <= wavelength_bound(jband-1)) then
          write(nulerr, '(a,a)') '*** Error: wavenumber_bound passed to ', &
               &  'define_sw_albedo_intervals must be monotonically increasing'
          call radiation_abort()
        end if
      end do

      ! Loop over shortwave bands
      do jband = 1,this%n_bands_sw
        ! Compute mid-point of band in wavenumber space (cm-1)
        wavenumber_mid = 0.5_jprb * (this%wavenumber1_sw(jband) &
             &                     + this%wavenumber2_sw(jband))
        iinterval = 1
        ! Convert wavelength (m) into wavenumber (cm-1) at the lower
        ! bound of the albedo interval
        wavenumber_bound = 0.01_jprb / wavelength_bound(iinterval)
        ! Find the albedo interval that has the largest overlap with
        ! the shortwave band; this approach assumes that the albedo
        ! intervals are larger than the shortwave spectral bands
        do while (wavenumber_bound >= wavenumber_mid &
             &    .and. iinterval < ninterval)
          iinterval = iinterval + 1
          if (iinterval < ninterval) then
            wavenumber_bound = 0.01_jprb / wavelength_bound(iinterval)
          else
            ! For the last interval there is no lower bound
            wavenumber_bound = 0.0_jprb
          end if
        end do
        ! Save the index of the band corresponding to the albedo
        ! interval and move onto the next shortwave band
        this%i_albedo_from_band_sw(jband) = i_band_in(iinterval)
      end do
    end if

    ! Define how many bands to use for reporting surface downwelling
    ! fluxes for canopy radiation scheme
    if (this%use_canopy_full_spectrum_sw) then
      this%n_canopy_bands_sw = this%n_g_sw
    else if (allocated(this%i_albedo_from_band_sw)) then
      this%n_canopy_bands_sw = maxval(this%i_albedo_from_band_sw)
    else
      this%n_canopy_bands_sw = 1
    end if

    if (this%iverbose >= 2) then
      write(nulout, '(a,i0,a)',advance="no") 'Mapping from ', this%n_bands_sw, &
           &  ' shortwave bands to albedo bands:'
      do jband = 1,this%n_bands_sw
        write(nulout,'(a,i0)',advance="no") ' ',this%i_albedo_from_band_sw(jband)
      end do
      write(nulout, '()')
    end if

  end subroutine define_sw_albedo_intervals

  !---------------------------------------------------------------------
  ! The input longwave surface emissivity coming in is likely to be in
  ! different spectral intervals to the gas model in the radiation
  ! scheme. We assume that the input emissivity is defined within
  ! "ninterval" spectral intervals covering the wavelength range 0 to
  ! infinity, but allow for the possibility that two intervals may be
  ! indexed back to the same emissivity band (which can happen if the
  ! atmospheric window uses one value and all other parts of the
  ! spectrum another).
  subroutine define_lw_emiss_intervals(this, ninterval, wavelength_bound, &
       &                                i_band_in)

    use radiation_io, only : nulout, nulerr, radiation_abort

    class(config_type),   intent(inout) :: this
    ! Number of spectral intervals in which emissivity is defined
    integer,              intent(in)    :: ninterval
    ! Monotonically increasing wavelength bounds between intervals,
    ! not including the outer bounds (which are assumed to be zero and
    ! infinity)
    real(jprb), optional, intent(in)    :: wavelength_bound(ninterval-1)
    ! The emissivity band indices corresponding to each interval
    integer,    optional, intent(in)    :: i_band_in(ninterval)

    integer :: jband
    integer :: iinterval
    real(jprb) :: wavenumber_bound ! cm-1
    real(jprb) :: wavenumber_mid   ! cm-1 
    if (.not. allocated(this%i_emiss_from_band_lw)) then
      allocate(this%i_emiss_from_band_lw(this%n_bands_lw))
    end if

    if (ninterval == 1) then 
      ! Only one interval, so wavelength bounds are not needed
      if (present(i_band_in)) then
        this%i_emiss_from_band_lw = i_band_in(1)
      else
        ! If there is only one emissivity interval then the index to it
        ! should always be 1
        this%i_emiss_from_band_lw = 1
      end if
    else if (.not. present(i_band_in) .or. .not. present(wavelength_bound)) then
      write(nulerr, '(a,a)') '*** Error: define_lw_emiss_intervals called', &
           ' with ninterval>1, but wavelength_bound or i_band_in not provided'
      call radiation_abort()
    else
      ! Check wavelength is monotonically increasing
      do jband = 2,ninterval-1
        if (wavelength_bound(jband) <= wavelength_bound(jband-1)) then
          write(nulerr, '(a,a)') '*** Error: wavenumber_bound passed to ', &
               &  'define_lw_emiss_intervals must be monotonically increasing'
          call radiation_abort()
        end if
      end do

      ! Loop over longwave bands
      do jband = 1,this%n_bands_lw
        ! Compute mid-point of band in wavenumber space (cm-1)
        wavenumber_mid = 0.5_jprb * (this%wavenumber1_lw(jband) &
             &                     + this%wavenumber2_lw(jband))
        iinterval = 1
        ! Convert wavelength (m) into wavenumber (cm-1) at the lower
        ! bound of the emissivity interval
        wavenumber_bound = 0.01_jprb / wavelength_bound(iinterval)
        ! Find the emissivity interval that has the largest overlap with
        ! the longwave band; this approach assumes that the emissivity
        ! intervals are larger than the longwave spectral bands
        do while (wavenumber_bound >= wavenumber_mid &
             &    .and. iinterval < ninterval)
          iinterval = iinterval + 1
          if (iinterval < ninterval) then
            wavenumber_bound = 0.01_jprb / wavelength_bound(iinterval)
          else
            ! For the last interval there is no lower bound
            wavenumber_bound = 0.0_jprb
          end if
        end do
        ! Save the index of the band corresponding to the emissivity
        ! interval and move onto the next longwave band
        this%i_emiss_from_band_lw(jband) = i_band_in(iinterval)
      end do
    end if

    ! Define how many bands to use for reporting surface downwelling
    ! fluxes for canopy radiation scheme
    if (this%use_canopy_full_spectrum_lw) then
      this%n_canopy_bands_lw = this%n_g_lw
    elseif (allocated(this%i_emiss_from_band_lw)) then
      this%n_canopy_bands_lw = maxval(this%i_emiss_from_band_lw)
    else
      this%n_canopy_bands_lw = 1
    end if

    if (this%iverbose >= 2) then
      write(nulout, '(a,i0,a)',advance="no") 'Mapping from ', this%n_bands_lw, &
           &  ' longwave bands to emissivity bands:'
      do jband = 1,this%n_bands_lw
        write(nulout,'(a,i0)',advance="no") ' ',this%i_emiss_from_band_lw(jband)
      end do
      write(nulout, '()')
    end if

  end subroutine define_lw_emiss_intervals


  !---------------------------------------------------------------------
  ! Return the 0-based index for str in enum_str, or abort if it is
  ! not found
  subroutine get_enum_code(str, enum_str, var_name, icode)

    use radiation_io, only : nulerr, radiation_abort

    character(len=*), intent(in)  :: str
    character(len=*), intent(in)  :: enum_str(0:)
    character(len=*), intent(in)  :: var_name
    integer,          intent(out) :: icode

    integer :: jc
    logical :: is_not_found

    ! If string is empty then we don't modify icode but assume it has
    ! a sensible default value
    if (len_trim(str) > 1) then
      is_not_found = .true.

      do jc = 0,size(enum_str)-1
        if (trim(str) == trim(enum_str(jc))) then
          icode = jc
          is_not_found = .false.
          exit
        end if
      end do
      if (is_not_found) then
        write(nulerr,'(a,a,a,a,a)',advance="no") '*** Error: ', trim(var_name), &
             &  ' must be one of: "', enum_str(0), '"'
        do jc = 1,size(enum_str)-1
          write(nulerr,'(a,a,a)',advance="no") ', "', trim(enum_str(jc)), '"'
        end do
        write(nulerr,'(a)') ''
        call radiation_abort('Radiation configuration error')
      end if
    end if

  end subroutine get_enum_code


  !---------------------------------------------------------------------
  ! Print one line of information: logical
  subroutine print_logical(message_str, name, val)
    use radiation_io, only : nulout
    character(len=*),   intent(in) :: message_str
    character(len=*),   intent(in) :: name
    logical,            intent(in) :: val
    character(4)                   :: on_or_off
    character(NPrintStringLen)     :: str
    if (val) then
      on_or_off = ' ON '
    else
      on_or_off = ' OFF'
    end if
    write(str, '(a,a4)') message_str, on_or_off
    write(nulout,'(a,a,a,a,l1,a)') str, ' (', name, '=', val,')'
  end subroutine print_logical


  !---------------------------------------------------------------------
  ! Print one line of information: integer
  subroutine print_integer(message_str, name, val)
    use radiation_io, only : nulout
    character(len=*),   intent(in) :: message_str
    character(len=*),   intent(in) :: name
    integer,            intent(in) :: val
    character(NPrintStringLen)     :: str
    write(str, '(a,a,i0)') message_str, ' = ', val
    write(nulout,'(a,a,a,a)') str, ' (', name, ')'
  end subroutine print_integer


  !---------------------------------------------------------------------
  ! Print one line of information: real
  subroutine print_real(message_str, name, val)
    use parkind1,     only : jprb
    use radiation_io, only : nulout
    character(len=*),   intent(in) :: message_str
    character(len=*),   intent(in) :: name
    real(jprb),         intent(in) :: val
    character(NPrintStringLen)     :: str
    write(str, '(a,a,g8.3)') message_str, ' = ', val
    write(nulout,'(a,a,a,a)') str, ' (', name, ')'
  end subroutine print_real


  !---------------------------------------------------------------------
  ! Print one line of information: enum
  subroutine print_enum(message_str, enum_str, name, val)
    use radiation_io, only : nulout
    character(len=*),   intent(in) :: message_str
    character(len=*),   intent(in) :: enum_str(0:)
    character(len=*),   intent(in) :: name
    integer,            intent(in) :: val

    character(NPrintStringLen)     :: str
    write(str, '(a,a,a)') message_str, ' ', enum_str(val)
    write(nulout,'(a,a,a,a,i0,a)') str, ' (', name, '=', val,')'
  end subroutine print_enum

end module radiation_config
