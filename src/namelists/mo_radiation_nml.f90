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

  USE mo_kind,               ONLY: wp
  USE mo_exception,          ONLY: finish
  USE mo_impl_constants,     ONLY: MAX_CHAR_LENGTH
  USE mo_mpi,                ONLY: p_pe, p_io
  USE mo_run_nml,            ONLY: iforcing, inwp
  USE mo_namelist,           ONLY: position_nml, positioned
  USE mo_io_units,           ONLY: nnml, nnml_output
  USE mo_physical_constants, ONLY: amd, amco2, amch4, amn2o, amo2
  USE mo_master_nml,         ONLY: lrestart
  USE mo_io_restart_namelist,ONLY: open_tmpfile, store_and_close_namelist, & 
                                 & open_and_restore_namelist, close_tmpfile

  IMPLICIT NONE
  PUBLIC

  ! 1.0 NAMELIST global variables and parameters
  ! --------------------------------
  !
  ! -- Switches for solar irradiation
  !
  LOGICAL :: ldiur       = .TRUE.  !< .TRUE. : with diurnal cycle
  !                                !< .FALSE.: zonally averaged irradiation
  !
  ! -- Switches for Earth orbit
  !
  INTEGER :: nmonth      =  0      !< i=0    : Earth circles on orbit, i.e. with annual cycle
  !                                !< i=1-12 : Earth orbit position fixed for month i
  !
  LOGICAL :: lyr_perp    = .FALSE. !< .FALSE.: transient Earth orbit following vsop87
  !                                !  .TRUE. : Earth orbit of year yr_perp of the vsop87 orbit
  !                                !           is perpetuated
  INTEGER :: yr_perp     = -99999  !< year used for lyr_perp = .TRUE.
  !
  !

  LOGICAL :: lradforcing(2) = (/.FALSE.,.FALSE./) !< diagnostic of instantaneous
  !                                               !< aerosol solar (lradforcing(1)) and
  !                                               !< thermal (lradforcing(2)) radiation forcing
  ! nmonth currently works for zonal mean ozone and the orbit (year 1987) only
  INTEGER :: isolrad     =  0      !< mode of solar constant calculation
  !< default is rrtm solar constant
  !
  ! --- Switches for radiative agents
  !     irad_x=0 : radiation uses tracer x = 0
  !     irad_x=1 : radiation uses tracer x from a tracer variable
  !     irad_x>1 : radiation uses tracer x following external specifications of various kinds:
  !                - globally constant  or spatially varying
  !                - constant in time, constant annual cycle, or transient
  !
  INTEGER  :: irad_h2o   = 1  !< water vapor, clouds and ice for radiation
  INTEGER  :: irad_co2   = 2  !< CO2
  INTEGER  :: irad_ch4   = 3  !< CH4
  INTEGER  :: irad_n2o   = 3  !< N2O
  INTEGER  :: irad_o3    = 0  !< O3
  INTEGER  :: irad_o2    = 2  !< O2
  INTEGER  :: irad_cfc11 = 2  !< CFC 11
  INTEGER  :: irad_cfc12 = 2  !< CFC 12
  INTEGER  :: irad_aero  = 2  !< aerosols
  !
  ! --- Default gas volume mixing ratios - 1990 values (CMIP5)
  !
  REAL(wp) :: vmr_co2    =  353.9e-06_wp     !< CO2
  REAL(wp) :: vmr_ch4    = 1693.6e-09_wp     !< CH4
  REAL(wp) :: vmr_n2o    =  309.5e-09_wp     !< N20
  REAL(wp) :: vmr_o2     =    0.20946_wp     !< O2
  REAL(wp) :: vmr_cfc11  =  252.8e-12_wp     !< CFC 11
  REAL(wp) :: vmr_cfc12  =  466.2e-12_wp     !< CFC 12
  !
  ! --- Time control
  !
  REAL(wp) :: dt_rad = 7200._wp        !< time interval of full radiation computation 
                                       !< given in seconds 
  !
  ! --- Different specifications of the zenith angle
  INTEGER  :: izenith = 3  ! circular orbit, no seasonal cycle but with diurnal cycle 
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

  ! 2.0 Non NAMELIST global variables
  ! ---------------------------------
  !
  ! --- solar ativity
  !
  REAL(wp) :: ssi(14)  !< spectrally resolved solar irradiance (SSI) [W/m2]
  !                    !< at 1 AU distance from the sun

  REAL(wp) :: tsi      !< total solar irradiance (TSI) [W/m2]
  !                    !< at 1 AU distance from the sun
  !                    !< = SUM(ssi(:))
  !
!!$  ! --- other parameters
!!$  !
!!$  REAL(wp) :: flx_ratio_cur
!!$  REAL(wp) :: flx_ratio_rad
!!$  REAL(wp) :: decl_sun_cur                  !< solar declination at current time step
  !
  !
  ! 3.0 Variables computed by routines in mo_radiation (export to submodels)
  ! --------------------------------
  !
  REAL(wp) :: mmr_co2, mmr_ch4, mmr_n2o, mmr_o2                ! setup_radiation


  PUBLIC:: read_radiation_namelist

CONTAINS
  !>
  !! "read_radiation_nml" reads the radiation_nml namelist from the namelist file.
  !!
  !! "read_radiation_nml" should be called if "lrad", the main switch for using radiative
  !! forcing, is set to true.
  !!
  !! In case of a restart, the namelist is first read from the global attributes of the
  !! restart file
  !!
  !! @par Revision History
  !! <Description of activity> by <name, affiliation> (<YYYY-MM-DD>)
  !!
  SUBROUTINE read_radiation_nml

    INTEGER :: ist, funit   !< status variable for namelist positioning

    ! For nwp, we want to have seasonal orbit and diurnal cycle as default:
    IF (iforcing==inwp) izenith=4
   
    !----------------------------------------------------------------
    ! If this is a resumed integration, overwrite the defaults above
    ! by values in the previous integration.
    !----------------------------------------------------------------
    IF (lrestart) THEN
      funit = open_and_restore_namelist('radiation_nml')
      READ(funit,NML=radiation_nml)
      CALL close_tmpfile(funit)
    END IF

    !---------------------------------------------------------------------
    ! Read user's (new) specifications (Done so far by all MPI processes)
    !---------------------------------------------------------------------
    CALL position_nml ('radiation_nml', STATUS=ist)
    SELECT CASE (ist)
    CASE (positioned)
      READ (nnml, radiation_nml)
    END SELECT

    ! Checks (none so far)

    !-----------------------------------------------------
    ! Store the namelist for restart
    !-----------------------------------------------------
    funit = open_tmpfile()
    WRITE(funit,NML=radiation_nml)
    CALL store_and_close_namelist(funit, 'radiation_nml')

    ! Set dependent variables
    ! -----------------------
    mmr_co2 = vmr_co2 * amco2/amd
    mmr_ch4 = vmr_ch4 * amch4/amd
    mmr_n2o = vmr_n2o * amn2o/amd
    mmr_o2  = vmr_o2  * amo2 /amd

  END SUBROUTINE read_radiation_nml

!!$  !>
!!$  !! "read_restart_radiation_nml" reads the parameters of the radiation
!!$  !! namelist from the global attributes of a restart file in netcdf format.
!!$  !!
!!$  !! @par Revision History
!!$  !! <Description of activity> by <name, affiliation> (<YYYY-MM-DD>)
!!$  !!
!!$  SUBROUTINE read_restart_radiation_nml
!!$
!!$    CALL finish('mo_radiation_nml/read_restart_radiation_nml', &
!!$      &         'This subroutine is not yet available')
!!$
!!$  END SUBROUTINE read_restart_radiation_nml
!!$
!!$  !>
!!$  !! "write_restart_radiation_nml" writes the parameters of the radiation
!!$  !! namelist as global attributes to a restart file in netcdf format.
!!$  !!
!!$  !! @par Revision History
!!$  !! <Description of activity> by <name, affiliation> (<YYYY-MM-DD>)
!!$  !!
!!$  SUBROUTINE write_restart_radiation_nml
!!$
!!$    CALL finish('mo_radiation_nml/write_restart_radiation_nml', &
!!$      &         'This subroutine is not yet available')
!!$
!!$  END SUBROUTINE write_restart_radiation_nml
!!$
!!$  !>
!!$  !! "write_rawdata_radiation_nml" writes the parameters of the radiation
!!$  !! namelist as global attributes to a raw data file in netcdf format.
!!$  !!
!!$  !! @par Revision History
!!$  !! <Description of activity> by <name, affiliation> (<YYYY-MM-DD>)
!!$  !!
!!$  SUBROUTINE write_rawdata_radiation_nml
!!$
!!$    CALL finish('mo_radiation_nml/write_rawdata_radiation_nml', &
!!$      &         'This subroutine is not yet available')
!!$
!!$  END SUBROUTINE write_rawdata_radiation_nml

  !-------------------------------------------------------------------------
  !
  !
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
  SUBROUTINE read_radiation_namelist
    !
    INTEGER :: istat, funit

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = 'mo_radiation_nml: read_radiation_namelist'

    !-----------------------------------------------------------------------

    !-----------------------!
    ! 1. default settings   !
    !-----------------------!

    ! For nwp, we want to have seasonal orbit and diurnal cycle as default:
    IF (iforcing==inwp) izenith=4


    !------------------------------------------------------------------
    ! 2. If this is a resumed integration, overwrite the defaults above 
    !    by values used in the previous integration.
    !------------------------------------------------------------------
    IF (lrestart) THEN
      funit = open_and_restore_namelist('radiation_nml')
      READ(funit,NML=radiation_nml)
      CALL close_tmpfile(funit)
    END IF


    !--------------------------------------------------------------------
    ! 3. Read user's (new) specifications (Done so far by all MPI processes)
    !--------------------------------------------------------------------
    CALL position_nml ('radiation_nml', status=istat)
    SELECT CASE (istat)
    CASE (POSITIONED)
      READ (nnml, radiation_nml)
    END SELECT


    !----------------------------------------------------
    ! 4. Fill the configuration state
    !----------------------------------------------------


    !-----------------------------------------------------
    ! 5. Store the namelist for restart
    !-----------------------------------------------------
    funit = open_tmpfile()
    WRITE(funit,NML=radiation_nml)                    
    CALL store_and_close_namelist(funit, 'radiation_nml') 


    ! 6. write the contents of the namelist to an ASCII file
    !
    IF(p_pe == p_io) WRITE(nnml_output,nml=radiation_nml)


  END SUBROUTINE read_radiation_namelist

END MODULE mo_radiation_nml
