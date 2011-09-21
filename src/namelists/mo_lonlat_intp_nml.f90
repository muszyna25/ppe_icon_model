!>
!! Contains the setup of the variables for lonlat interpolation of output.
!!
!! @par Revision History
!!
!! @par Copyright
!! 2002-2006 by DWD and MPI-M
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
!!
MODULE mo_lonlat_intp_nml

  USE mo_kind,               ONLY: wp
  USE mo_impl_constants,     ONLY: max_dom
  USE mo_io_units,           ONLY: nnml, nnml_output
  USE mo_namelist,           ONLY: position_nml, positioned, open_nml, close_nml
  USE mo_mpi,                ONLY: my_process_is_stdio, p_n_work
  USE mo_master_control,     ONLY: is_restart_run
  USE mo_io_restart_namelist,ONLY: open_tmpfile, store_and_close_namelist,   &
                                 & open_and_restore_namelist, close_tmpfile
  USE mo_lonlat_intp_config, ONLY: lonlat_intp_config
  USE mo_exception,          ONLY: finish
  USE mo_math_utilities,     ONLY: t_lon_lat_grid

  IMPLICIT NONE
  PUBLIC :: read_lonlat_intp_namelist
  CHARACTER(len=*), PARAMETER, PRIVATE :: version = '$Id$'

  !-------------------------------------------------------------------------
  ! Namelist variables
  !-------------------------------------------------------------------------

  LOGICAL              :: llonlat_enabled(max_dom) ! Flag. True, if interpolation onto lon-lat grid is enabled
  CHARACTER (len=1024) :: lonlat_var_list        ! string with a list of variables for lon-lat interpolation
  REAL(wp)             :: lon_delta(max_dom)     ! lon-lat grid resolution,                unit:DEGREE
  REAL(wp)             :: lat_delta(max_dom)     ! lon-lat grid resolution,                unit:DEGREE
  REAL(wp)             :: lon_sw_corner(max_dom) ! south western corner of area (lon/lat), unit:DEGREE
  REAL(wp)             :: lat_sw_corner(max_dom) ! south western corner of area (lon/lat), unit:DEGREE
  REAL(wp)             :: lon_poleN(max_dom)     ! position of north pole (lon,lat),       unit:DEGREE
  REAL(wp)             :: lat_poleN(max_dom)     ! position of north pole (lon,lat),       unit:DEGREE
  INTEGER              :: lon_dimen(max_dom)     ! grid dimensions
  INTEGER              :: lat_dimen(max_dom)     ! grid dimensions

  !> Namelist for interpolation of output variables
  NAMELIST/lonlat_intp_nml/ llonlat_enabled, lonlat_var_list,         &
    &                       lon_delta, lon_sw_corner, lon_poleN, lon_dimen,    &
    &                       lat_delta, lat_sw_corner, lat_poleN, lat_dimen

CONTAINS
  !>
  !! Read Namelist for lon-lat interpolation of output.
  !!
  !! This subroutine
  !! - reads the Namelist for lon-lat interpolation
  !! - sets default values
  !! - potentially overwrites the defaults by values used in a
  !!   previous integration (if this is a resumed run)
  !! - reads the user's (new) specifications
  !! - stores the Namelist for restart
  !! - fills the configuration state (partly)
  !!
  !! @par Revision History
  !!  by F. Prill, DWD (2011-09-15)
  !!
  SUBROUTINE read_lonlat_intp_namelist( filename )

    CHARACTER(LEN=*), INTENT(IN)   :: filename
    ! local variables
    CHARACTER(len=*), PARAMETER    :: routine = &
      &  'mo_lonlat_intp_nml:read_lonlat_intp_namelist'
    INTEGER                        :: istat, funit, idom
    REAL(wp)                       :: pi_180
    TYPE (t_lon_lat_grid), POINTER :: grid ! lon-lat grid

    pi_180 = ATAN(1._wp)/45._wp

    !-----------------------
    ! 1. default settings
    !-----------------------

    llonlat_enabled    = .FALSE.
    lonlat_var_list    = "'PS', 'Q7', 'normal_velocity'"
    lon_delta(:)       =    2._wp
    lat_delta(:)       =    2._wp
    lon_sw_corner(:)   = -180._wp
    lat_sw_corner(:)   =  -90._wp
    lon_poleN(:)       =    0._wp
    lat_poleN(:)       =   90._wp
    lon_dimen(:)       = 181
    lat_dimen(:)       =  91

    !------------------------------------------------------------------
    ! 2. If this is a resumed integration, overwrite the defaults above
    !    by values used in the previous integration.
    !------------------------------------------------------------------
    IF (is_restart_run()) THEN
      funit = open_and_restore_namelist('lonlat_intp_nml')
      READ(funit,NML=lonlat_intp_nml)
      CALL close_tmpfile(funit)
    END IF

    !-------------------------------------------------------------------------
    ! 3. Read user's (new) specifications (Done so far by all MPI processes)
    !-------------------------------------------------------------------------
    CALL open_nml(TRIM(filename))
    CALL position_nml ('lonlat_intp_nml', status=istat)
    SELECT CASE (istat)
    CASE (POSITIONED)
      READ (nnml, lonlat_intp_nml)
    END SELECT
    CALL close_nml

    !----------------------------------------------------
    ! 4. Fill the configuration state
    !----------------------------------------------------

    ! fill in values for each model domain:
    DO idom=1,max_dom
      lonlat_intp_config(idom)%l_enabled = llonlat_enabled(idom)

      ! string with a list of variable names due for lon-lat
      ! interpolation
      ! Note: This list is identical on all patches.
      lonlat_intp_config(idom)%zlist = lonlat_var_list

      grid => lonlat_intp_config(idom)%lonlat_grid

      grid%delta(1:2)     = (/ lon_delta(idom),     lat_delta(idom)     /) * pi_180
      grid%sw_corner(1:2) = (/ lon_sw_corner(idom), lat_sw_corner(idom) /) * pi_180
      grid%poleN(1:2)     = (/ lon_poleN(idom),     lat_poleN(idom)     /) * pi_180
      grid%dimen(1:2)     = (/ lon_dimen(idom),     lat_dimen(idom)     /)

    END DO

    !-----------------------------------------------------
    ! 5. Store the namelist for restart
    !-----------------------------------------------------
    IF(my_process_is_stdio())  THEN
      funit = open_tmpfile()
      WRITE(funit,NML=lonlat_intp_nml)
      CALL store_and_close_namelist(funit, 'lonlat_intp_nml')
    ENDIF

    !-----------------------------------------------------
    ! 6. write the contents of the namelist to an ASCII file
    !-----------------------------------------------------
    IF(my_process_is_stdio()) THEN
      WRITE(nnml_output,nml=lonlat_intp_nml)
    END IF

  END SUBROUTINE read_lonlat_intp_namelist

END MODULE mo_lonlat_intp_nml
