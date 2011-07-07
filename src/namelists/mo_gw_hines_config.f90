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
!! @par Copyright
!! 2002-2011 by DWD and MPI-M
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
MODULE mo_gw_hines_config

  USE mo_kind,           ONLY: wp
  USE mo_impl_constants, ONLY: MAX_NTRACER, MAX_CHAR_LENGTH, max_dom

  IMPLICIT NONE
  PUBLIC
  CHARACTER(len=*),PARAMETER,PRIVATE :: version = '$Id$'


  !!--------------------------------------------------------------------------
  !! Basic configuration setup for radiation
  !!--------------------------------------------------------------------------


  TYPE t_gw_hines_config
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

  LOGICAL :: lradforcing(2) !< diagnostic of instantaneous
  !                         !< aerosol solar (lradforcing(1)) and
  !                         !< thermal (lradforcing(2)) radiation forcing
  ! nmonth currently works for zonal mean ozone and the orbit (year 1987) only
  INTEGER :: isolrad        !< mode of solar constant calculation
  !< default is rrtm solar constant
  !
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
  !
  ! --- Default gas volume mixing ratios - 1990 values (CMIP5)
  !
  REAL(wp) :: vmr_co2     !< CO2
  REAL(wp) :: vmr_ch4     !< CH4
  REAL(wp) :: vmr_n2o     !< N20
  REAL(wp) :: vmr_o2      !< O2
  REAL(wp) :: vmr_cfc11   !< CFC 11
  REAL(wp) :: vmr_cfc12   !< CFC 12
  !
  ! --- Time control
  !
  REAL(wp) :: dt_rad      !< time interval of full radiation computation 
                          !< given in seconds 
  !
  ! --- Different specifications of the zenith angle
  INTEGER  :: izenith     ! circular orbit, no seasonal cycle but with diurnal cycle 


  END TYPE t_gw_hines_config
  !>
  !!
  TYPE(t_gw_hines_config) :: gw_hines_config(max_dom)

END MODULE mo_gw_hines_config
