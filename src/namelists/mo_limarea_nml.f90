!>
!! Contains the setup of variables related to the boundary relaxation of limited
!! area models from an external dataset
!!
!! @author S. Brdar (DWD)
!!
!!
!! @par Revision History
!! Initial release by S. Brdar (2013-06-13)
!!
!! @par Copyright
!! 2002-2013 by DWD and MPI-M
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
MODULE mo_limarea_nml

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: finish
  USE mo_io_units,            ONLY: nnml, nnml_output, filename_max
  USE mo_namelist,            ONLY: position_nml, positioned, open_nml, close_nml
  USE mo_mpi,                 ONLY: my_process_is_stdio
  USE mo_master_control,      ONLY: is_restart_run
  USE mo_impl_constants,      ONLY: max_dom
  USE mo_io_restart_namelist, ONLY: open_tmpfile, store_and_close_namelist     , &
                                  & open_and_restore_namelist, close_tmpfile
  USE mo_limarea_config,      ONLY: latbc_config
  USE mo_util_string,         ONLY: MAX_STRING_LEN
  USE mo_nml_annotate,        ONLY: temp_defaults, temp_settings

  IMPLICIT NONE
  PRIVATE
  PUBLIC read_limarea_namelist

  CHARACTER(len=*), PARAMETER, PRIVATE :: version = '$Id$'

  !------------------------------------------------------------------------
  ! Namelist variables
  !------------------------------------------------------------------------
  INTEGER                         :: itype_latbc    ! type of limited area boundary nudging
  REAL(wp)                        :: dtime_latbc    ! dt between two consequtive external latbc files
  INTEGER                         :: nlev_latbc     ! number of vertical levels in boundary data
  CHARACTER(LEN=filename_max)     :: latbc_filename ! prefix of latbc files
  CHARACTER(LEN=MAX_STRING_LEN)   :: latbc_path     ! directory containing external latbc files

  NAMELIST /limarea_nml/ itype_latbc, dtime_latbc, nlev_latbc, latbc_filename, latbc_path

CONTAINS
  !>
  !!
  SUBROUTINE read_limarea_namelist( filename )

    CHARACTER(LEN=*), INTENT(IN) :: filename
    INTEGER :: istat, funit
    CHARACTER(LEN=*), PARAMETER :: routine = 'mo_limarea_nml:read_limarea_namelist'

    !------------------------------------------------------------
    ! Default settings
    !------------------------------------------------------------
    itype_latbc      = 0
    dtime_latbc      = 43200._wp
    nlev_latbc       = 0
    latbc_filename   = "<path>prepiconR<nroot>B<jlev>_DOM<dom>_<timestamp>.nc"
    latbc_path       = "<path>"

    !------------------------------------------------------------------
    ! If this is a resumed integration, overwrite the defaults above 
    ! by values used in the previous integration.
    !------------------------------------------------------------------
    IF (is_restart_run()) THEN
      funit = open_and_restore_namelist('limarea_nml')
      READ(funit,NML=limarea_nml)
      CALL close_tmpfile(funit)
    END IF

    !------------------------------------------------------------------------
    ! Read user's (new) specifications. (Done so far by all MPI processes)
    !------------------------------------------------------------------------
    CALL open_nml(TRIM(filename))
    CALL position_nml ('limarea_nml', status=istat)
    IF (my_process_is_stdio()) WRITE(temp_defaults(), limarea_nml)  ! write defaults to temporary text file
    SELECT CASE (istat)
    CASE (POSITIONED)
      READ (nnml, limarea_nml)
      IF (my_process_is_stdio()) WRITE(temp_settings(), limarea_nml)  ! write settings to temporary text file
    END SELECT
    CALL close_nml

    !----------------------------------------------------
    ! Fill the configuration state
    !----------------------------------------------------
    latbc_config% itype_latbc     = itype_latbc
    latbc_config% dtime_latbc     = dtime_latbc
    latbc_config% nlev_in         = nlev_latbc
    latbc_config% latbc_filename  = latbc_filename
    latbc_config% latbc_path      = TRIM(latbc_path)//'/'

    !-----------------------------------------------------
    ! Store the namelist for restart
    !-----------------------------------------------------
    IF(my_process_is_stdio()) THEN
      funit = open_tmpfile()
      WRITE(funit,NML=limarea_nml)
      CALL store_and_close_namelist(funit, 'limarea_nml')
    ENDIF

    !-----------------------------------------------------
    !write the contents of the namelist to an ASCII file
    !-----------------------------------------------------
    IF(my_process_is_stdio()) WRITE(nnml_output,nml=limarea_nml)

  END SUBROUTINE read_limarea_namelist

END MODULE mo_limarea_nml
