!>
!! Main switches of the MPI physics package, for turning on/off
!! the main parameterized processes.
!!
!! All the switches have default values, but can be changed by model
!! user via namelist "mpi_phy_nml". This module contains a subroutine
!! that reads the namelist and makes necessary modifications of the
!! default or user-specified values. Note that output of the namelist values
!! to the stdout and the ASCII file is not done here, but after
!! all namelists have been read in and the cross-check has been finished.
!!
!! @author Hui Wan, MPI-M
!! @author Marco, Giorgetta, MPI-M
!!
!! @par Revision History
!! -Variables taken from mo_control and mo_param_switches of ECHAM6
!!  by Hui Wan, MPI (2010-07)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_mpi_phy_nml

  USE mo_mpi_phy_config,     ONLY: mpi_phy_config
  
  USE mo_mpi,                ONLY: my_process_is_stdio
  USE mo_io_units,           ONLY: nnml, nnml_output
  USE mo_namelist,           ONLY: position_nml, positioned, open_nml, close_nml
  USE mo_master_control,     ONLY: use_restart_namelists
  USE mo_restart_namelist,   ONLY: open_tmpfile, store_and_close_namelist, &
                                 & open_and_restore_namelist, close_tmpfile
  USE mo_nml_annotate,       ONLY: temp_defaults, temp_settings

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: read_mpi_phy_namelist

  NAMELIST /mpi_phy_nml/ mpi_phy_config

CONTAINS
  !>
  !!
  SUBROUTINE read_mpi_phy_namelist( filename )

    CHARACTER(LEN=*), INTENT(IN) :: filename

    INTEGER :: istat
    INTEGER :: funit
    INTEGER :: iunit


    !------------------------------------------------------------------
    ! 1. If this is a resumed integration, overwrite the initialization
    !    values with values used in the previous integration.
    !------------------------------------------------------------------
    IF (use_restart_namelists()) THEN
      funit = open_and_restore_namelist('mpi_phy_nml')
      READ(funit,NML=mpi_phy_nml)
      CALL close_tmpfile(funit)
    END IF
   
    !------------------------------------------------------------------
    ! 2. Read user's (new) specifications (done by all MPI processes)
    !------------------------------------------------------------------
    CALL open_nml(TRIM(filename))
    CALL position_nml ('mpi_phy_nml', STATUS=istat)
    IF (my_process_is_stdio()) THEN
      iunit = temp_defaults()
      WRITE(iunit, mpi_phy_nml)      ! write defaults to temporary text file
    END IF
    SELECT CASE (istat)
    CASE (positioned)
      READ (nnml, mpi_phy_nml)       ! overwrite default settings
      IF (my_process_is_stdio()) THEN
        iunit = temp_settings()
        WRITE(iunit, mpi_phy_nml)    ! write settings to temporary text file
      END IF
    END SELECT
    CALL close_nml

    !------------------------------------------------------------------
    ! 3. Store the namelist for restart
    !------------------------------------------------------------------
    IF(my_process_is_stdio())  THEN
      funit = open_tmpfile()
      WRITE(funit,NML=mpi_phy_nml)
      CALL store_and_close_namelist(funit, 'mpi_phy_nml')
    ENDIF

    !------------------------------------------------------------------
    ! 4. Write the namelist to an ASCII file
    !------------------------------------------------------------------
    IF ( my_process_is_stdio() ) WRITE(nnml_output,nml=mpi_phy_nml)

  END SUBROUTINE read_mpi_phy_namelist

END MODULE mo_mpi_phy_nml
