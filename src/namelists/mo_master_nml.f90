!>
!! @brief Master namelist.
!!        
!! @par Revision History
!! Created by Rene Redler (2011-03-22)
!!
!! @par Copyright
!! 2010-2011 by DWD and MPI-M
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
MODULE mo_master_nml

  USE mo_impl_constants, ONLY: MAX_CHAR_LENGTH, SUCCESS
  USE mo_exception,      ONLY: warning, message_text, finish
  USE mo_io_units,       ONLY: filename_max, nnml
  USE mo_namelist,       ONLY: open_nml, position_nml, POSITIONED

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: version = '$Id$'

  PUBLIC

  !-------------------------------------------------------------------------
  ! Namelist variables 
  !-------------------------------------------------------------------------
  LOGICAL :: lrestart

  ! Component models

  CHARACTER(len=132) :: atmo_name
  CHARACTER(len=filename_max) :: atmo_namelist_filename
  CHARACTER(len=filename_max) :: atmo_restart_info_filename
  LOGICAL :: l_atmo_active
  INTEGER :: atmo_min_rank
  INTEGER :: atmo_max_rank
  INTEGER :: atmo_inc_rank

  CHARACTER(len=132) :: ocean_name
  CHARACTER(len=filename_max) :: ocean_namelist_filename
  CHARACTER(len=filename_max) :: ocean_restart_info_filename
  LOGICAL :: l_ocean_active
  INTEGER :: ocean_min_rank
  INTEGER :: ocean_max_rank
  INTEGER :: ocean_inc_rank

  NAMELIST /master_nml/ l_atmo_active, atmo_name, atmo_namelist_filename,    &
                      & atmo_min_rank, atmo_max_rank, atmo_inc_rank,         &
                      & atmo_restart_info_filename,                          &
                      & l_ocean_active, ocean_name, ocean_namelist_filename, &
                      & ocean_min_rank, ocean_max_rank, ocean_inc_rank,      &
                      & ocean_restart_info_filename,                         &
                      & lrestart

CONTAINS
  !>
  !! Initialization of variables that contain general information
  !! about the coupled model run. The configuration is read from
  !! namelist 'master_nml'.
  !!
  !! @par Revision History
  !!
  INTEGER FUNCTION read_master_namelist(namelist_filename)
    
    CHARACTER(LEN=*), INTENT(in) :: namelist_filename
    !
    ! Local variables
    !
    INTEGER :: istat
    CHARACTER(len=MAX_CHAR_LENGTH) :: routine = 'mo_master_nml:read_master_namelist'

    !------------------------------------------------------------
    ! 1. Set default values
    !------------------------------------------------------------
    lrestart      = .FALSE.

    atmo_name     = ""
    atmo_namelist_filename     = "NAMELIST_ICON"
    atmo_restart_info_filename = "restart.info"
    l_atmo_active = .FALSE.
    atmo_min_rank = 0
    atmo_max_rank = -1
    atmo_inc_rank = 1

    ocean_name     = ""
    ocean_namelist_filename     = "NAMELIST_ICON"
    ocean_restart_info_filename = "restart.info"
    l_ocean_active = .FALSE.
    ocean_min_rank = 0
    ocean_max_rank = -1
    ocean_inc_rank = 1

    !------------------------------------------------------------------
    ! 2. Read user's specifications (done so far by all MPI processes)
    !------------------------------------------------------------------

    OPEN( nnml, FILE=TRIM(namelist_filename), IOSTAT=istat, &
        & STATUS='old', ACTION='read', DELIM='apostrophe')

    IF (istat/=0) THEN
      CALL warning(namelist_filename,"not found")
      read_master_namelist=-1
      RETURN
    ENDIF
    
    CALL position_nml('master_nml',STATUS=istat)
    
    IF (istat/=POSITIONED) THEN

      CALL finish( TRIM(routine), &
                 & 'Namelist master_nml not found in file '  &
                 & //TRIM(namelist_filename) )

      read_master_namelist=-2
      RETURN      
    ENDIF
        
    READ (nnml, master_nml)
    CLOSE (nnml, IOSTAT=istat)

    read_master_namelist=SUCCESS

  END FUNCTION read_master_namelist

END MODULE mo_master_nml
