!>
!! This module provides parameters controlling the radiation interface.
!!
!! @author Bjorn Stevens, MPI-M, Hamburg (2009-09-19):
!!
!!
!! @par Revision History
!! - New module, extracted from mo_radiation, Martin Schultz, FZJ, Juelich (2010-04-13)
!! - Added parameter for local solar constant, Hauke Schmidt, MPI-M, Hamburg (2010-0?-??)
!! - Added decl_sun_cur (for MOZ photolysis), Martin Schultz, FZJ, Juelich (2010-06-02)
!! - Modified for ICON, Marco Giorgetta, MPI-M, Hamburg (2010-07-24)
!!   - added subroutine read_radiation_nml
!! - Modified for ICON, Hui Wan, MPI-M, Hamburg (2010-11-06)
!!   - added namelist variable dt_rad
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_radiation_nml

    USE mo_radiation_config, ONLY: config_ldiur      => ldiur,       &
                                 & config_nmonth     => nmonth,      &
                                 & config_lyr_perp   => lyr_perp,    &
                                 & config_yr_perp    => yr_perp,     &
                                 & config_isolrad    => isolrad,     &
                                 & config_albedo_type=> albedo_type, &
                                 & config_direct_albedo       => direct_albedo,       &
                                 & config_direct_albedo_water => direct_albedo_water, &
                                 & config_albedo_whitecap     => albedo_whitecap,     &
                                 & config_icld_overlap        => icld_overlap,        &
                                 & config_islope_rad => islope_rad,  &
                                 & config_irad_h2o   => irad_h2o,    &
                                 & config_irad_co2   => irad_co2,    &
                                 & config_irad_ch4   => irad_ch4,    &
                                 & config_irad_n2o   => irad_n2o,    &
                                 & config_irad_o3    => irad_o3,     &
                                 & config_irad_o2    => irad_o2,     &
                                 & config_irad_cfc11 => irad_cfc11,  &
                                 & config_irad_cfc12 => irad_cfc12,  &
                                 & config_irad_aero  => irad_aero,   &
                                 & config_lrad_aero_diag => lrad_aero_diag,  &
                                 & config_ighg       => ighg,        &
                                 & config_ghg_filename   => ghg_filename,    &
                                 & config_vmr_co2    => vmr_co2,     &
                                 & config_vmr_ch4    => vmr_ch4,     &
                                 & config_vmr_n2o    => vmr_n2o,     &
                                 & config_vmr_o2     => vmr_o2,      &
                                 & config_vmr_cfc11  => vmr_cfc11,   &
                                 & config_vmr_cfc12  => vmr_cfc12,   &
                                 & config_izenith    => izenith,     &
                                 & config_mmr_co2    => mmr_co2,     &
                                 & config_mmr_ch4    => mmr_ch4,     &
                                 & config_mmr_n2o    => mmr_n2o,     &
                                 & config_mmr_o2     => mmr_o2,      &
                                 & config_mmr_cfc11  => mmr_cfc11,   &
                                 & config_mmr_cfc12  => mmr_cfc12,   &
                                 & config_llw_cloud_scat => llw_cloud_scat, &
                                 & config_iliquid_scat => iliquid_scat, &
                                 & config_iice_scat => iice_scat,    &
                                 & config_ecrad_data_path => ecrad_data_path

  USE mo_kind,               ONLY: wp
  USE mo_impl_constants,     ONLY: MAX_CHAR_LENGTH
  USE mo_mpi,                ONLY: my_process_is_stdio
  USE mo_namelist,           ONLY: position_nml, positioned, open_nml, close_nml
  USE mo_io_units,           ONLY: nnml, nnml_output
  USE mo_physical_constants, ONLY: amd, amco2, amch4, amn2o, amo2, amc11, amc12
  USE sfc_terra_data,        ONLY: csalb, csalb1, csalb2
  USE mo_master_control,     ONLY: use_restart_namelists
  USE mo_restart_namelist,   ONLY: open_tmpfile, store_and_close_namelist, &
                                 & open_and_restore_namelist, close_tmpfile
  USE mo_nml_annotate,       ONLY: temp_defaults, temp_settings
  USE mo_io_units,           ONLY: filename_max

  IMPLICIT NONE
  PRIVATE
  PUBLIC:: read_radiation_namelist

  !-----------------------------------
  ! namelist variables and parameters
  !-----------------------------------
  !
  ! -- Switches for solar irradiation
  !
  LOGICAL :: ldiur     !< .TRUE. : with diurnal cycle
  !                    !< .FALSE.: zonally averaged irradiation
  !
  ! -- Switches for Earth orbit
  !
  INTEGER :: nmonth    !< i=0    : Earth circles on orbit, i.e. with annual cycle
  !                    !< i=1-12 : Earth orbit position fixed for month i
  !
  LOGICAL :: lyr_perp  !< .FALSE.: transient Earth orbit following vsop87
  !                    !  .TRUE. : Earth orbit of year yr_perp of the vsop87 orbit
  !                    !           is perpetuated
  INTEGER :: yr_perp   !< year used for lyr_perp = .TRUE.

  ! nmonth currently works for zonal mean ozone and the orbit (year 1987) only
  INTEGER :: isolrad   !< mode of solar constant calculation
  !< default is rrtm solar constant
  !
  INTEGER :: albedo_type ! 1: albedo based on surface-type specific set of constants
                         !    (see )
                         ! 2: Modis albedo

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
                                 ! 4: Yang et al (2008) for sea water, Ritter and Geleyn (1992) for lakes

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
  INTEGER  :: irad_h2o
  INTEGER  :: irad_co2
  INTEGER  :: irad_ch4
  INTEGER  :: irad_n2o
  INTEGER  :: irad_o3
  INTEGER  :: irad_o2
  INTEGER  :: irad_cfc11
  INTEGER  :: irad_cfc12
  INTEGER  :: irad_aero
  LOGICAL  :: lrad_aero_diag
  !
  ! --- Select dynamic greenhouse gases scenario (read from file)
  !     ighg = 0 : select default gas volume mixing ratios - 1990 values (CMIP5)
  !     ighg = 1 : transient CMIP5 scenario from file
  !
  INTEGER  :: ighg
  !
  ! --- Name of the file that contains  dynamic greenhouse values
  !
  CHARACTER(LEN=filename_max)  :: ghg_filename
  !
  ! --- Default gas volume mixing ratios - 1990 values (CMIP5)
  !
!DR preliminary restart fix
#ifdef __SX__
  INTEGER, PARAMETER :: qp = SELECTED_REAL_KIND(24, 307)
  REAL(qp) :: vmr_co2
  REAL(qp) :: vmr_ch4
  REAL(qp) :: vmr_n2o
  REAL(qp) :: vmr_o2
  REAL(qp) :: vmr_cfc11
  REAL(qp) :: vmr_cfc12
#else
  REAL(wp) :: vmr_co2
  REAL(wp) :: vmr_ch4
  REAL(wp) :: vmr_n2o
  REAL(wp) :: vmr_o2
  REAL(wp) :: vmr_cfc11
  REAL(wp) :: vmr_cfc12
#endif
  !
  ! --- Time control
  !
  !
  ! --- Different specifications of the zenith angle
  INTEGER  :: izenith
  !
  ! ecRad specific configuration
  LOGICAL  :: llw_cloud_scat
  INTEGER  :: iliquid_scat
  INTEGER  :: iice_scat
  CHARACTER(len=MAX_CHAR_LENGTH) :: ecrad_data_path
  !
  NAMELIST /radiation_nml/ ldiur, nmonth,         &
    &                      lyr_perp, yr_perp,     &
    &                      isolrad,               &
    &                      albedo_type,           &
    &                      direct_albedo,         &
    &                      direct_albedo_water,   &
    &                      albedo_whitecap,       &
    &                      irad_h2o,              &
    &                      irad_co2,   vmr_co2,   &
    &                      irad_ch4,   vmr_ch4,   &
    &                      irad_n2o,   vmr_n2o,   &
    &                      irad_o3,               &
    &                      irad_o2,    vmr_o2,    &
    &                      irad_cfc11, vmr_cfc11, &
    &                      irad_cfc12, vmr_cfc12, &
    &                      irad_aero,             &
    &                      lrad_aero_diag,        &
    &                      ighg,                  &
    &                      ghg_filename,          &
    &                      izenith, icld_overlap, &
    &                      islope_rad,            &
    &                      llw_cloud_scat,        &
    &                      iliquid_scat,          &
    &                      iice_scat,             &
    &                      ecrad_data_path

CONTAINS

  !>
  !! Read Namelist for radiation. 
  !!
  !! This subroutine 
  !! - reads the Namelist for radiation
  !! - sets default values
  !! - potentially overwrites the defaults by values used in a 
  !!   previous integration (if this is a resumed run)
  !! - reads the user's (new) specifications
  !! - stores the Namelist for restart
  !! - fills the configuration state (partly)    
  !!
  !! @par Revision History
  !!  by Daniel Reinert, DWD (2011-06-07)
  !!
  SUBROUTINE read_radiation_namelist( filename )

    CHARACTER(LEN=*), INTENT(IN) :: filename
    INTEGER :: istat, funit
    INTEGER :: iunit

    !0!CHARACTER(len=*), PARAMETER ::  &
    !0!  &  routine = 'mo_radiation_nml:read_radiation_namelist'

    !-----------------------
    ! 1. default settings   
    !-----------------------
    ldiur          = .TRUE.
    nmonth         =  0   
    lyr_perp       = .FALSE.
    yr_perp        = -99999

    isolrad        = 0
    albedo_type    = 1
    direct_albedo  = 4      ! Parameterization after Briegleb and Ramanathan (1992)
    direct_albedo_water = 4 ! Parameterization after Yang et al (2008) for sea water, and Ritter and Geleyn (1992) for lakes
    albedo_whitecap= 0      ! no whitecap albedo from breaking ocean waves
    icld_overlap   = 2      ! generalized random overlap
    islope_rad     = 0      ! no slope correction

    irad_h2o    = 1
    irad_co2    = 2
    irad_ch4    = 3
    irad_n2o    = 3
    irad_o3     = 0
    irad_o2     = 2
    irad_cfc11  = 2
    irad_cfc12  = 2
    irad_aero   = 2
    lrad_aero_diag = .FALSE.

    ighg        = 0
    ghg_filename= 'bc_greenhouse_gases.nc'

    vmr_co2     = 348.0e-06_wp
    vmr_ch4     = 1650.0e-09_wp
    vmr_n2o     =  306.0e-09_wp
    vmr_o2      =    0.20946_wp
    vmr_cfc11   =  214.5e-12_wp
    vmr_cfc12   =  371.1e-12_wp

    izenith     = 4  ! Default: seasonal orbit and diurnal cycle

    llw_cloud_scat  = .FALSE.
    iliquid_scat    = 0
    iice_scat       = 0
    ecrad_data_path = '.'

    !------------------------------------------------------------------
    ! 2. If this is a resumed integration, overwrite the defaults above 
    !    by values used in the previous integration.
    !------------------------------------------------------------------

    IF (use_restart_namelists()) THEN
      funit = open_and_restore_namelist('radiation_nml')
      READ(funit,NML=radiation_nml)
      CALL close_tmpfile(funit)
    END IF


    !--------------------------------------------------------------------
    ! 3. Read user's (new) specifications (Done so far by all MPI processes)
    !--------------------------------------------------------------------
    CALL open_nml(TRIM(filename))
    CALL position_nml ('radiation_nml', status=istat)
    IF (my_process_is_stdio()) THEN
      iunit = temp_defaults()
      WRITE(iunit, radiation_nml)  ! write defaults to temporary text file
    END IF
    SELECT CASE (istat)
    CASE (POSITIONED)
      READ (nnml, radiation_nml)                                      ! overwrite default settings
      IF (my_process_is_stdio()) THEN
        iunit = temp_settings()
        WRITE(iunit, radiation_nml)  ! write settings to temporary text file
      END IF
    END SELECT
    CALL close_nml

    !----------------------------------------------------
    ! 4. Fill the configuration state
    !----------------------------------------------------

    config_ldiur      = ldiur
    config_nmonth     = nmonth
    config_lyr_perp   = lyr_perp
    config_yr_perp    = yr_perp
    config_isolrad    = isolrad
    config_albedo_type= albedo_type
    config_direct_albedo       = direct_albedo
    config_direct_albedo_water = direct_albedo_water
    config_albedo_whitecap     = albedo_whitecap
    config_icld_overlap        = icld_overlap
    config_islope_rad = islope_rad
    config_irad_h2o   = irad_h2o
    config_irad_co2   = irad_co2
    config_irad_ch4   = irad_ch4
    config_irad_n2o   = irad_n2o
    config_irad_o3    = irad_o3
    config_irad_o2    = irad_o2
    config_irad_cfc11 = irad_cfc11
    config_irad_cfc12 = irad_cfc12
    config_irad_aero  = irad_aero
    config_lrad_aero_diag = lrad_aero_diag
    config_ighg       = ighg
    config_ghg_filename   = ghg_filename
    config_vmr_co2    = vmr_co2
    config_vmr_ch4    = vmr_ch4
    config_vmr_n2o    = vmr_n2o
    config_vmr_o2     = vmr_o2
    config_vmr_cfc11  = vmr_cfc11
    config_vmr_cfc12  = vmr_cfc12
    config_mmr_co2    = vmr_co2   * amco2/amd
    config_mmr_ch4    = vmr_ch4   * amch4/amd
    config_mmr_n2o    = vmr_n2o   * amn2o/amd
    config_mmr_o2     = vmr_o2    * amo2 /amd
    config_mmr_cfc11  = vmr_cfc11 * amc11/amd
    config_mmr_cfc12  = vmr_cfc12 * amc12/amd

    config_izenith    = izenith

    config_llw_cloud_scat  = llw_cloud_scat
    config_iliquid_scat    = iliquid_scat
    config_iice_scat       = iice_scat
    config_ecrad_data_path = TRIM(ecrad_data_path)

    IF ( direct_albedo_water == 3 ) THEN
      csalb => csalb2
    ELSE
      csalb => csalb1
    ENDIF

    !-----------------------------------------------------
    ! 5. Store the namelist for restart
    !-----------------------------------------------------
    IF(my_process_is_stdio())  THEN
      funit = open_tmpfile()
      WRITE(funit,NML=radiation_nml)
      CALL store_and_close_namelist(funit, 'radiation_nml') 
    ENDIF
    ! 6. write the contents of the namelist to an ASCII file
    !
    IF(my_process_is_stdio()) WRITE(nnml_output,nml=radiation_nml)

  END SUBROUTINE read_radiation_namelist

END MODULE mo_radiation_nml
