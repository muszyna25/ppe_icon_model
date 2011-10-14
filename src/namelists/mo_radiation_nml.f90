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
!! @par Copyright
!! 2002-2010 by DWD and MPI-M
!! This software is provided for non-commercial use only.
!! See the LICENSE and the WARRANTY conditions.
!!
!! @par License
!! The use of ICON is hereby granted free of charge for an unlimited time,
!! provided the following rules are accepted and applied:
!! <ol>
!! <li> You may use or modify this code for your own non commercial and non
!!    violent purposes.
!! <li> The code may not be re-distributed without the consent of the authors.
!! <li> The copyright notice and statement of authorship must appear in all
!!    copies.
!! <li> You accept the warranty conditions (see WARRANTY).
!! <li> In case you intend to use the code commercially, we oblige you to sign
!!    an according license agreement with DWD and MPI-M.
!! </ol>
!!
!! @par Warranty
!! This code has been tested up to a certain level. Defects and weaknesses,
!! which may be included in the code, do not establish any warranties by the
!! authors.
!! The authors do not make any warranty, express or implied, or assume any
!! liability or responsibility for the use, acquisition or application of this
!! software.
!!
MODULE mo_radiation_nml

    USE mo_radiation_config, ONLY: config_ldiur      => ldiur,      & 
                                 & config_nmonth     => nmonth,     &
                                 & config_lyr_perp   => lyr_perp,   &
                                 & config_yr_perp    => yr_perp,    &
                                 & config_isolrad    => isolrad,    &
                                 & config_irad_h2o   => irad_h2o,   &
                                 & config_irad_co2   => irad_co2,   &
                                 & config_irad_ch4   => irad_ch4,   &
                                 & config_irad_n2o   => irad_n2o,   &
                                 & config_irad_o3    => irad_o3,    &
                                 & config_irad_o2    => irad_o2,    &
                                 & config_irad_cfc11 => irad_cfc11, &
                                 & config_irad_cfc12 => irad_cfc12, &
                                 & config_irad_aero  => irad_aero,  &
                                 & config_vmr_co2    => vmr_co2,    &
                                 & config_vmr_ch4    => vmr_ch4,    &
                                 & config_vmr_n2o    => vmr_n2o,    &
                                 & config_vmr_o2     => vmr_o2,     &
                                 & config_vmr_cfc11  => vmr_cfc11,  &
                                 & config_vmr_cfc12  => vmr_cfc12,  &
                                 & config_dt_rad     => dt_rad,     &
                                 & config_izenith    => izenith,    &
                                 & config_mmr_co2    => mmr_co2,    &
                                 & config_mmr_ch4    => mmr_ch4,    &
                                 & config_mmr_n2o    => mmr_n2o,    &
                                 & config_mmr_o2     => mmr_o2

  USE mo_kind,               ONLY: wp
  USE mo_mpi,                ONLY: my_process_is_stdio
  USE mo_namelist,           ONLY: position_nml, positioned, open_nml, close_nml
  USE mo_io_units,           ONLY: nnml, nnml_output
  USE mo_physical_constants, ONLY: amd, amco2, amch4, amn2o, amo2
  USE mo_master_control,     ONLY: is_restart_run
  USE mo_io_restart_namelist,ONLY: open_tmpfile, store_and_close_namelist, &
                                 & open_and_restore_namelist, close_tmpfile

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
  !
  ! --- Default gas volume mixing ratios - 1990 values (CMIP5)
  !
  REAL(wp) :: vmr_co2
  REAL(wp) :: vmr_ch4
  REAL(wp) :: vmr_n2o
  REAL(wp) :: vmr_o2
  REAL(wp) :: vmr_cfc11
  REAL(wp) :: vmr_cfc12
  !
  ! --- Time control
  !
  REAL(wp) :: dt_rad  !< time interval of full radiation computation given in seconds 
  !
  ! --- Different specifications of the zenith angle
  INTEGER  :: izenith
  !
  NAMELIST /radiation_nml/ ldiur, nmonth,         &
    &                      lyr_perp, yr_perp,     &
    &                      isolrad,               &
    &                      irad_h2o,              &
    &                      irad_co2,   vmr_co2,   &
    &                      irad_ch4,   vmr_ch4,   &
    &                      irad_n2o,   vmr_n2o,   &
    &                      irad_o3,               &
    &                      irad_o2,    vmr_o2,    &
    &                      irad_cfc11, vmr_cfc11, &
    &                      irad_cfc12, vmr_cfc12, &
    &                      irad_aero,             &
    &                      dt_rad, izenith

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

    !0!CHARACTER(len=*), PARAMETER ::  &
    !0!  &  routine = 'mo_radiation_nml:read_radiation_namelist'

    !-----------------------
    ! 1. default settings   
    !-----------------------
    ldiur          = .TRUE.
    nmonth         =  0   
    lyr_perp       = .FALSE.
    yr_perp        = -99999

    isolrad    = 0

    irad_h2o   = 1
    irad_co2   = 2
    irad_ch4   = 3
    irad_n2o   = 3
    irad_o3    = 0
    irad_o2    = 2
    irad_cfc11 = 2
    irad_cfc12 = 2
    irad_aero  = 2

    vmr_co2    =  353.9e-06_wp
    vmr_ch4    = 1693.6e-09_wp
    vmr_n2o    =  309.5e-09_wp
    vmr_o2     =    0.20946_wp
    vmr_cfc11  =  252.8e-12_wp
    vmr_cfc12  =  466.2e-12_wp

    dt_rad = 7200._wp

    izenith = 3  ! Default: get it working KF
!    izenith = 4  ! Default: seasonal orbit and diurnal cycle

    !------------------------------------------------------------------
    ! 2. If this is a resumed integration, overwrite the defaults above 
    !    by values used in the previous integration.
    !------------------------------------------------------------------
    IF (is_restart_run()) THEN
      funit = open_and_restore_namelist('radiation_nml')
      READ(funit,NML=radiation_nml)
      CALL close_tmpfile(funit)
    END IF

    !--------------------------------------------------------------------
    ! 3. Read user's (new) specifications (Done so far by all MPI processes)
    !--------------------------------------------------------------------
    CALL open_nml(TRIM(filename))
    CALL position_nml ('radiation_nml', status=istat)
    SELECT CASE (istat)
    CASE (POSITIONED)
      READ (nnml, radiation_nml)
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
    config_irad_h2o   = irad_h2o
    config_irad_co2   = irad_co2
    config_irad_ch4   = irad_ch4
    config_irad_n2o   = irad_n2o
    config_irad_o3    = irad_o3
    config_irad_o2    = irad_o2
    config_irad_cfc11 = irad_cfc11
    config_irad_cfc12 = irad_cfc12
    config_irad_aero  = irad_aero
    config_vmr_co2    = vmr_co2
    config_vmr_ch4    = vmr_ch4
    config_vmr_n2o    = vmr_n2o
    config_vmr_o2     = vmr_o2
    config_vmr_cfc11  = vmr_cfc11
    config_vmr_cfc12  = vmr_cfc12
    config_mmr_co2    = vmr_co2 * amco2/amd
    config_mmr_ch4    = vmr_ch4 * amch4/amd
    config_mmr_n2o    = vmr_n2o * amn2o/amd
    config_mmr_o2     = vmr_o2  * amo2 /amd
    config_dt_rad     = dt_rad
    config_izenith    = izenith

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
