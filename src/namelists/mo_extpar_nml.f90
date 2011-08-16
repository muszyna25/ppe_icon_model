!>
!!        
!! @par Revision History
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
!!
MODULE mo_extpar_nml

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: finish
  USE mo_io_units,            ONLY: nnml, nnml_output
  USE mo_namelist,            ONLY: position_nml, positioned, open_nml, close_nml
  USE mo_mpi,                 ONLY: my_process_is_stdio
  USE mo_master_control,      ONLY: is_restart_run

  USE mo_io_restart_namelist, ONLY: open_tmpfile, store_and_close_namelist         , &
                                  & open_and_restore_namelist, close_tmpfile

  USE mo_extpar_config,       ONLY: config_itopo              => itopo             , &
                                  & config_fac_smooth_topo    => fac_smooth_topo   , &
                                  & config_n_iter_smooth_topo => n_iter_smooth_topo, &
                                  & config_l_emiss            => l_emiss

  IMPLICIT NONE
  PRIVATE
  PUBLIC read_extpar_namelist

  CHARACTER(len=*), PARAMETER, PRIVATE :: version = '$Id$'

  !------------------------------------------------------------------------
  ! Namelist variables
  !------------------------------------------------------------------------
  INTEGER  :: itopo   ! 0: topography specified by analytical functions,
                      ! 1: topography read from netcdf files

  REAL(wp) :: fac_smooth_topo
  INTEGER  :: n_iter_smooth_topo
  LOGICAL  :: l_emiss ! if true: read external emissivity map

  NAMELIST /extpar_nml/ itopo, fac_smooth_topo,n_iter_smooth_topo,l_emiss

CONTAINS
  !>
  !!
  SUBROUTINE read_extpar_namelist( filename )

    CHARACTER(LEN=*), INTENT(IN) :: filename
    INTEGER :: istat, funit
    CHARACTER(LEN=*), PARAMETER :: routine = 'mo_extpar_nml:read_extpar_namelist'

    !------------------------------------------------------------
    ! Default settings
    !------------------------------------------------------------
    itopo              = 0
    fac_smooth_topo    = 0.015625_wp
    n_iter_smooth_topo = 2
    l_emiss            = .TRUE.

    !------------------------------------------------------------------
    ! If this is a resumed integration, overwrite the defaults above 
    ! by values used in the previous integration.
    !------------------------------------------------------------------
    IF (is_restart_run()) THEN
      funit = open_and_restore_namelist('extpar_nml')
      READ(funit,NML=extpar_nml)
      CALL close_tmpfile(funit)
    END IF

    !------------------------------------------------------------------------
    ! Read user's (new) specifications. (Done so far by all MPI processes)
    !------------------------------------------------------------------------
    CALL open_nml(TRIM(filename))
    CALL position_nml ('extpar_nml', status=istat)
    SELECT CASE (istat)
    CASE (POSITIONED)
      READ (nnml, extpar_nml)
    END SELECT
    CALL close_nml

    !----------------------------------------------------
    ! Sanity check
    !----------------------------------------------------
    SELECT CASE (itopo)
    CASE (0,1,2) !OK
    CASE default
      CALL finish(TRIM(routine),'Wrong value for itopo. Must be 0 or 1.')
    END SELECT

    !----------------------------------------------------
    ! Fill the configuration state
    !----------------------------------------------------
    config_itopo              = itopo 
    config_fac_smooth_topo    = fac_smooth_topo 
    config_n_iter_smooth_topo = n_iter_smooth_topo
    config_l_emiss            = l_emiss

    !-----------------------------------------------------
    ! Store the namelist for restart
    !-----------------------------------------------------
    IF(my_process_is_stdio())  THEN
      funit = open_tmpfile()
      WRITE(funit,NML=extpar_nml)
      CALL store_and_close_namelist(funit, 'extpar_nml')
    ENDIF
    !write the contents of the namelist to an ASCII file
    IF(my_process_is_stdio()) WRITE(nnml_output,nml=extpar_nml)

  END SUBROUTINE read_extpar_namelist

END MODULE mo_extpar_nml
