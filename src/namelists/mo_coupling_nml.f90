!>
!!        Contains the variables to set up the coupling.
!!
!!        
!! @par Revision History
!!   Created by Rene Redler (2011-03-22)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
!!

MODULE mo_coupling_nml

  !-------------------------------------------------------------------------
  !
  !    ProTeX FORTRAN source: Style 2
  !    modified for ICON project, DWD/MPI-M 2006
  !
  !-------------------------------------------------------------------------

  USE mo_impl_constants,  ONLY: max_char_length
  USE mo_io_units,        ONLY: nnml
  USE mo_namelist,        ONLY: open_nml, close_nml, position_nml, POSITIONED

  USE mo_coupling_config, ONLY: config_coupled_mode

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: read_coupling_namelist

CONTAINS

  !>
  !!  Initialization of variables that contain general information.
  !!
  !!               Initialization of variables that contain general information
  !!               about the coupled model run. The configuration is read from
  !!               namelist 'icon_cpl'.
  !!
  !! @par Revision History
  !!

  SUBROUTINE read_coupling_namelist (namelist_filename)

    CHARACTER(LEN=*), INTENT(in) :: namelist_filename

    !
    ! Local variables
    !

    LOGICAL :: coupled_mode
    INTEGER :: istat

    CHARACTER(len=max_char_length), PARAMETER :: &
         &   routine = 'mo_coupling_nml:read_coupling_namelist'

    NAMELIST /coupling_mode_nml/ coupled_mode

    !--------------------------------------------------------------------
    ! 1. Set default values
    !--------------------------------------------------------------------

    coupled_mode  = .FALSE.

    !--------------------------------------------------------------------
    ! 2. Read user's (new) specifications (done so far by all MPI processes)
    !--------------------------------------------------------------------

#ifdef YAC_coupling

    CALL open_nml (TRIM(namelist_filename))
    
    CALL position_nml('coupling_mode_nml',STATUS=istat)
    IF (istat==POSITIONED) THEN
      READ (nnml, coupling_mode_nml)
    ENDIF

    CALL close_nml

#endif

    config_coupled_mode = coupled_mode

  END SUBROUTINE read_coupling_namelist

END MODULE mo_coupling_nml
