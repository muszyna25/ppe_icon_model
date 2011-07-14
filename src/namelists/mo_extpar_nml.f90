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

  USE mo_kind,           ONLY: wp
  USE mo_exception,      ONLY: finish
  USE mo_io_units,       ONLY: nnml, nnml_output
  USE mo_namelist,       ONLY: position_nml, positioned, open_nml, close_nml
  USE mo_mpi,            ONLY: p_pe, p_io
  USE mo_master_nml,     ONLY: lrestart

  USE mo_io_restart_attributes, ONLY: get_restart_attribute
  USE mo_io_restart_namelist,   ONLY: open_tmpfile, store_and_close_namelist,   &
                                    & open_and_restore_namelist, close_tmpfile

  IMPLICIT NONE
 !PRIVATE
  PUBLIC
  PUBLIC read_extpar_namelist

  CHARACTER(len=*), PARAMETER, PRIVATE :: version = '$Id$'

  !------------------------------------------------------------------------
  ! Namelist variables
  !------------------------------------------------------------------------
  INTEGER  :: nml_itopo   ! 0: topography specified by analytical functions,
                          ! 1: topography read from netcdf files provided by Herrmann Asensio

  REAL(wp):: nml_fac_smooth_topo
  INTEGER :: nml_n_iter_smooth_topo

  NAMELIST/extpar_nml/ nml_itopo, nml_fac_smooth_topo,nml_n_iter_smooth_topo

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
    nml_itopo  = 0
    nml_fac_smooth_topo    = 0.015625_wp
    nml_n_iter_smooth_topo = 35

    !------------------------------------------------------------------
    ! If this is a resumed integration, overwrite the defaults above 
    ! by values used in the previous integration.
    !------------------------------------------------------------------
    IF (lrestart) THEN
      funit = open_and_restore_namelist('extpar_nml')
      READ(funit,NML=extpar_nml)
      CALL close_tmpfile(funit)
    END IF

    !------------------------------------------------------------------------
    ! Read user's (new) specifications. (Done so far by all MPI processors)
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
    SELECT CASE (nml_itopo)
    CASE (0,1) !OK
    CASE default
      CALL finish(TRIM(routine),'Wrong value for nml_itopo. Must be 0 or 1.')
    END SELECT

    !----------------------------------------------------
    ! Fill the configuration state
    !----------------------------------------------------


    !-----------------------------------------------------
    ! Store the namelist for restart
    !-----------------------------------------------------
    funit = open_tmpfile()
    WRITE(funit,NML=extpar_nml)
    CALL store_and_close_namelist(funit, 'extpar_nml')

    !write the contents of the namelist to an ASCII file
    IF(p_pe == p_io) WRITE(nnml_output,nml=extpar_nml)

  END SUBROUTINE read_extpar_namelist

END MODULE mo_extpar_nml
