!>
!! <Short description of module for listings and indices>
!!
!! <Describe the concepts of the procedures and algorithms used in the module.>
!! <Details of procedures are documented below with their definitions.>
!! <Include any applicable external references inline as module::procedure,>
!! <external_procedure(), or by using @see.>
!! <Don't forget references to literature.>
!!
!! @author <name, affiliation>
!! @author <name, affiliation>
!!
!!
!! @par Revision History
!! <Description of activity> by <name, affiliation> (<YYYY-MM-DD>)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_radiation_config

  USE mo_kind,           ONLY: wp
  USE mo_io_units,       ONLY: filename_max
  USE mo_impl_constants, ONLY: MAX_CHAR_LENGTH

  IMPLICIT NONE
  PUBLIC

  !--------------------------------------------------------------------------
  ! Basic configuration setup for radiation
  !--------------------------------------------------------------------------

  !TYPE t_radiation_config
    !
    ! -- Switches for solar irradiation
    !
    LOGICAL :: ldiur        !< .TRUE. : with diurnal cycle
    !                       !< .FALSE.: zonally averaged irradiation
    !
    ! -- Switches for Earth orbit
    !
    INTEGER :: nmonth       !< i=0    : Earth circles on orbit, i.e. with annual cycle
    !                                !< i=1-12 : Earth orbit position fixed for month i
    !
    LOGICAL :: lyr_perp     !< .FALSE.: transient Earth orbit following vsop87
    !                       !  .TRUE. : Earth orbit of year yr_perp of the vsop87 orbit
    !                       !           is perpetuated
    INTEGER :: yr_perp      !< year used for lyr_perp = .TRUE.
    !
    !
    LOGICAL :: lradforcing(2) = (/.FALSE.,.FALSE./) !< diagnostic of instantaneous
    !                                               !< aerosol solar (lradforcing(1)) and
    !                                               !< thermal (lradforcing(2)) radiation forcing 
    ! nmonth currently works for zonal mean ozone and the orbit (year 1987) only
    INTEGER :: isolrad = 3        !< mode of solar constant calculation
    !< default is rrtm solar constant
    !
    INTEGER :: albedo_type ! 1: albedo based on surface-type specific set of constants
                           !    (see )
                           ! 2: Modis albedo
                           ! 3: fixed albedo with value albedo_fixed

    REAL(wp) :: albedo_fixed   ! value of fixed albedo for albedo_type=3

    INTEGER :: direct_albedo   ! 1: SZA dependence according to Ritter and Geleyn (1992)
                               ! 2: limitation to diffuse albedo according to Zaengl 
                               !    applied to all land points
                               !    Ritter-Geleyn for ice 
                               ! 3: Parameterization after Yang et al (2008) for snow-free land points
                               !    limitation after Zaengl for snow-coverer points
                               !    Ritter-Geleyn implementation for ice
                               ! 4: Parameterization after Briegleb and Ramanathan (1992) for snow-free land points
                               !    limitation after Zaengl for snow-coverer points
  
    INTEGER :: direct_albedo_water ! 1: Ritter and Geleyn (1992)
                                   ! 2: Yang et al (2008)
                                   ! 3: Taylor et al (1996) as in IFS

    INTEGER :: albedo_whitecap ! 0: no whitecap albedo
                               ! 1: Seferian et al (2018) whitecaps albedo from breaking ocean waves

    INTEGER :: icld_overlap    ! method for cloud overlap calculation in shortwave part of RRTM
                               ! 1: maximum-random overlap
                               ! 2: generalized overlap (Hogan, Illingworth, 2000)
                               ! 3: maximum overlap
                               ! 4: random overlap

    INTEGER :: islope_rad      ! slope correction for surface radiation
                               ! 0: none
                               ! 1: slope correction for solar radiation without shading effects
                               ! option 2 is reserved for slope-dependent radiation with shading (not yet implemented)

    ! --- Switches for radiative agents
    !     irad_x=0 : radiation uses tracer x = 0
    !     irad_x=1 : radiation uses tracer x from a tracer variable
    !     irad_x>1 : radiation uses tracer x following external specifications of various kinds:
    !                - globally constant  or spatially varying
    !                - constant in time, constant annual cycle, or transient
    !
    INTEGER  :: irad_h2o    !< water vapor, clouds and ice for radiation
    INTEGER  :: irad_co2    !< CO2
    INTEGER  :: irad_ch4    !< CH4
    INTEGER  :: irad_n2o    !< N2O
    INTEGER  :: irad_o3     !< O3
    INTEGER  :: irad_o2     !< O2
    INTEGER  :: irad_cfc11  !< CFC 11
    INTEGER  :: irad_cfc12  !< CFC 12
    INTEGER  :: irad_aero   !< aerosols
    LOGICAL  :: lrad_aero_diag  !< diagnose aerosols
    !
    ! --- Name of the file that contains  dynamic greenhouse values
    !
    CHARACTER(LEN=filename_max)  :: ghg_filename
    !
    ! --- Default gas mixing ratios - 1990 values (CMIP5)
    !
    REAL(wp) :: vmr_co2  , mmr_co2   !< CO2
    REAL(wp) :: vmr_n2o  , mmr_n2o   !< N20
    REAL(wp) :: vmr_o2   , mmr_o2    !< O2
    REAL(wp) :: vmr_ch4  , mmr_ch4   !< CH4
    REAL(wp) :: vmr_cfc11, mmr_cfc11 !< CFC 11
    REAL(wp) :: vmr_cfc12, mmr_cfc12 !< CFC 12
    !
    ! --- Different specifications of the zenith angle
    INTEGER  :: izenith           ! circular orbit, no seasonal cycle but with diurnal cycle 
    REAL(wp) :: cos_zenith_fixed  ! fixed cosine of zenith angle for izenith=6
  
    ! 2.0 Non NAMELIST global variables and parameters
    ! --------------------------------

    REAL(wp) :: rad_csalbw(10) ! slope of solar albedo with respect to soil water content
                               ! as a function of depth of upper soil layer
  
    
    ! vertical profile parameters (vpp) of CH4 and N2O
    REAL(wp), PARAMETER :: vpp_ch4(3) = (/1.25e-01_wp,  683.0_wp, -1.43_wp/)
    REAL(wp), PARAMETER :: vpp_n2o(3) = (/1.20e-02_wp, 1395.0_wp, -1.43_wp/)
    !
    !
    ! --- solar activity for radiation time step (this time step can be
    !     different from the actual integration time step, in general it
    !     is in the future relative to actual time step)
    !
    REAL(wp) :: ssi_radt(14)  !< spectrally resolved solar irradiance (SSI) 
    !                         !< [W/m2] at 1 AU distance from the sun
  
    REAL(wp) :: tsi_radt !< total solar irradiance (TSI) [W/m2]
    !                    !< at 1 AU distance from the sun
    !                    !< = SUM(ssi_radt(:))
    !
    ! --- solar activity for actual time step (for calculation of heating
    !     rates)
    !
    REAL(wp) :: tsi
    !
    ! ecRad specific configuration
    LOGICAL  :: llw_cloud_scat          !< Do long wave cloud scattering?
    INTEGER  :: iliquid_scat            !< Optical properties for liquid cloud scattering
                                        !< 0: SOCRATES
                                        !< 1: Slingo (1989)
    INTEGER  :: iice_scat               !< Optical properties for ice cloud scattering
                                        !< 0: Fu et al. (1996)
                                        !< 1: Baran et al. (2016)
    CHARACTER(len=MAX_CHAR_LENGTH) :: &
      &  ecrad_data_path                !< Folder containing optical properties
    INTEGER  :: nproma_rad              !< subblock size used for the ecrad interface.
                                        !< If set negativ, the absolute value is considered as the number of subblocks.
    !

    !$ACC DECLARE COPYIN(vpp_ch4, vpp_n2o)
  
  !END TYPE t_radiation_config
  !>
  !!
  !TYPE(t_radiation_config) :: radiation_config(max_dom)

END MODULE mo_radiation_config
