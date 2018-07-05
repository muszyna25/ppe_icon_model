!>
!! Read configuration parameters as Fortran namelist from an external file. 
!!
!! @author Monika Esch, MPI-M
!!
!! @par Revision History
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_ccycle_nml

  USE mo_ccycle_config,      ONLY: name, ccycle_config
  
  USE mo_mpi,                ONLY: my_process_is_stdio
  USE mo_impl_constants,     ONLY: max_char_length
  USE mo_io_units,           ONLY: nnml, nnml_output
  USE mo_namelist,           ONLY: position_nml, positioned, open_nml, close_nml
  USE mo_master_control,     ONLY: use_restart_namelists
  USE mo_restart_namelist,   ONLY: open_tmpfile, store_and_close_namelist, &
                                 & open_and_restore_namelist, close_tmpfile
  USE mo_nml_annotate,       ONLY: temp_defaults, temp_settings

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: read_ccycle_nml

  NAMELIST /ccycle_nml/ ccycle_config

CONTAINS
  !>
  !!
  SUBROUTINE read_ccycle_nml( filename )

    CHARACTER(LEN=*), INTENT(IN)   :: filename
    CHARACTER(LEN=max_char_length) :: nml_name

    INTEGER :: istat
    INTEGER :: funit
    INTEGER :: iunit

    nml_name=name//'_nml'

    !------------------------------------------------------------------
    ! 1. If this is a resumed integration, overwrite the initialization
    !    values with values used in the previous integration.
    !------------------------------------------------------------------
    IF (use_restart_namelists()) THEN
      funit = open_and_restore_namelist(TRIM(nml_name))
      !
      READ(funit,NML=ccycle_nml)
      !
      CALL close_tmpfile(funit)
    END IF
   
    !------------------------------------------------------------------
    ! 2. Read user's (new) specifications (done by all MPI processes)
    !------------------------------------------------------------------
    CALL open_nml(TRIM(filename))
    CALL position_nml (TRIM(nml_name), STATUS=istat)
    IF (my_process_is_stdio()) THEN
      iunit = temp_defaults()
      !
      WRITE(iunit, NML=ccycle_nml)      ! write defaults to temporary text file
      !
    END IF
    SELECT CASE (istat)
    CASE (positioned)
      !
      READ (nnml, NML=ccycle_nml)       ! overwrite default settings
      !
      IF (my_process_is_stdio()) THEN
        iunit = temp_settings()
        !
        WRITE(iunit, NML=ccycle_nml)    ! write settings to temporary text file
        !
      END IF
    END SELECT
    CALL close_nml

    !------------------------------------------------------------------
    ! 3. Store the namelist for restart
    !------------------------------------------------------------------
    IF (my_process_is_stdio())  THEN
       funit = open_tmpfile()
       !
       WRITE(funit,NML=ccycle_nml)
       !
       CALL store_and_close_namelist(funit, TRIM(nml_name))
    END IF

    !------------------------------------------------------------------
    ! 4. Write the namelist to an ASCII file
    !------------------------------------------------------------------
    IF ( my_process_is_stdio() ) THEN
       !
       WRITE(nnml_output,NML=ccycle_nml)
       !
    END IF

  END SUBROUTINE read_ccycle_nml

END MODULE mo_ccycle_nml
