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
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
!!
MODULE mo_limarea_nml

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: finish
  USE mo_io_units,            ONLY: nnml, nnml_output, filename_max
  USE mo_namelist,            ONLY: position_nml, positioned, open_nml, close_nml
  USE mo_mpi,                 ONLY: my_process_is_stdio
  USE mo_master_control,      ONLY: use_restart_namelists
  USE mo_impl_constants,      ONLY: max_dom, MAX_CHAR_LENGTH
  USE mo_restart_namelist,    ONLY: open_tmpfile, store_and_close_namelist     , &
                                  & open_and_restore_namelist, close_tmpfile
  USE mo_limarea_config,      ONLY: latbc_config
  USE mo_nml_annotate,        ONLY: temp_defaults, temp_settings
  USE mtime,                  ONLY: MAX_DATETIME_STR_LEN

  IMPLICIT NONE
  PRIVATE
  PUBLIC read_limarea_namelist

  !------------------------------------------------------------------------
  ! Namelist variables
  !------------------------------------------------------------------------
  INTEGER                         :: itype_latbc         ! type of limited area boundary nudging
  REAL(wp)                        :: dtime_latbc         ! dt between two consequtive external latbc files
  INTEGER                         :: nlev_latbc          ! number of vertical levels in boundary data
  CHARACTER(LEN=filename_max)     :: latbc_filename      ! prefix of latbc files
  CHARACTER(LEN=MAX_CHAR_LENGTH)  :: latbc_path          ! directory containing external latbc files
  CHARACTER(LEN=FILENAME_MAX)     :: latbc_boundary_grid ! grid file defining the lateral boundary

  ! dictionary which maps internal variable names onto
  ! GRIB2 shortnames or NetCDF var names used for lateral boundary nudging.
  CHARACTER(LEN=filename_max) :: latbc_varnames_map_file  

  NAMELIST /limarea_nml/ itype_latbc, dtime_latbc, nlev_latbc, latbc_filename, &
    &                    latbc_path, latbc_boundary_grid, latbc_varnames_map_file

CONTAINS
  !>
  !!
  SUBROUTINE read_limarea_namelist( filename )

    CHARACTER(LEN=*), INTENT(IN) :: filename
    INTEGER :: istat, funit
    INTEGER :: iunit

    !------------------------------------------------------------
    ! Default settings
    !------------------------------------------------------------
    itype_latbc         = 0
    dtime_latbc         = 10800._wp
    nlev_latbc          = 0
    latbc_filename      = "prepiconR<nroot>B<jlev>_<y><m><d><h>.nc"
    latbc_path          = "./"
    latbc_boundary_grid = ""  ! empty string means: whole domain is read for lateral boundary
    latbc_varnames_map_file = " "

    !------------------------------------------------------------------
    ! If this is a resumed integration, overwrite the defaults above 
    ! by values used in the previous integration.
    !------------------------------------------------------------------
    IF (use_restart_namelists()) THEN
      funit = open_and_restore_namelist('limarea_nml')
      READ(funit,NML=limarea_nml)
      CALL close_tmpfile(funit)
    END IF

    !------------------------------------------------------------------------
    ! Read user's (new) specifications. (Done so far by all MPI processes)
    !------------------------------------------------------------------------
    CALL open_nml(TRIM(filename))
    CALL position_nml ('limarea_nml', status=istat)
    IF (my_process_is_stdio()) THEN
      iunit = temp_defaults()
      WRITE(iunit, limarea_nml)  ! write defaults to temporary text file
    END IF
    SELECT CASE (istat)
    CASE (POSITIONED)
      READ (nnml, limarea_nml)
      IF (my_process_is_stdio()) THEN
        iunit = temp_settings()
        WRITE(iunit, limarea_nml)  ! write settings to temporary text file
      END IF
    END SELECT
    CALL close_nml

    !----------------------------------------------------
    ! Fill the configuration state
    !----------------------------------------------------
    latbc_config% itype_latbc         = itype_latbc
    latbc_config% dtime_latbc         = dtime_latbc
    latbc_config% nlev_in             = nlev_latbc
    latbc_config% latbc_filename      = latbc_filename
    latbc_config% latbc_path          = TRIM(latbc_path)//'/'
    latbc_config% latbc_boundary_grid = latbc_boundary_grid
    latbc_config% lsparse_latbc       = (LEN_TRIM(latbc_boundary_grid) > 0)
    latbc_config% latbc_varnames_map_file = latbc_varnames_map_file

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
