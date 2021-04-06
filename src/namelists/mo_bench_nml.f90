!! Namelist for flags that disable parts of the ICON model inside the TIME_LOOP.
!! These are used at MeteoSwiss for bencharking while not all needed parts are ported to GPU

MODULE mo_bench_nml

  USE mo_io_units,            ONLY: nnml, nnml_output
  USE mo_namelist,            ONLY: position_nml, POSITIONED, open_nml, close_nml
  USE mo_master_control,      ONLY: use_restart_namelists
  USE mo_mpi,                 ONLY: my_process_is_stdio
  USE mo_nml_annotate,        ONLY: temp_defaults, temp_settings
  USE mo_bench_config,        ONLY: bench_config
  USE mo_restart_nml_and_att, ONLY: open_tmpfile, store_and_close_namelist,     &
  &                                 open_and_restore_namelist, close_tmpfile
    
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: read_bench_namelist

  LOGICAL :: d_unpb !! disable update_nwp_phy_bcs
  LOGICAL :: d_ndfo !! disable nwp_diag_for_output
  LOGICAL :: d_rld  !! disable recv_latbc_data
  LOGICAL :: d_n    !! disable nudging
  LOGICAL :: d_wnlo !! disable write_name_list_output
   
  NAMELIST /bench_nml/ d_unpb, d_ndfo, d_rld, d_n, d_wnlo

  CONTAINS

  SUBROUTINE read_bench_namelist( filename )

    CHARACTER(LEN=*), INTENT(IN) :: filename

    INTEGER :: istat, funit, iunit
    INTEGER :: jg

    !-----------------------
    ! 1. default settings   
    !-----------------------
    d_unpb = .FALSE.
    d_ndfo = .FALSE.
    d_rld = .FALSE.
    d_n = .FALSE.
    d_wnlo = .FALSE.

    !------------------------------------------------------------------
    ! 2. If this is a resumed integration, overwrite the defaults above 
    !    by values used in the previous integration.
    !------------------------------------------------------------------
    IF (use_restart_namelists()) THEN
      funit = open_and_restore_namelist('bench_nml')
      READ(funit,NML=bench_nml)
      CALL close_tmpfile(funit)
    END IF

    !--------------------------------------------------------------------
    ! 3. Read user's (new) specifications (Done so far by all MPI processes)
    !--------------------------------------------------------------------
    CALL open_nml(TRIM(filename))
    CALL position_nml ('bench_nml', STATUS=istat)
    IF (my_process_is_stdio()) THEN
      iunit = temp_defaults()
      WRITE(iunit, bench_nml)   ! write defaults to temporary text file
    END IF
    SELECT CASE (istat)
    CASE (POSITIONED)
      READ (nnml, bench_nml)                                        ! overwrite default settings
      IF (my_process_is_stdio()) THEN
        iunit = temp_settings()
        WRITE(iunit, bench_nml)    ! write settings to temporary text file
      END IF
    END SELECT
    CALL close_nml

    !----------------------------------------------------
    ! 4. Sanity check
    !----------------------------------------------------
    ! -

    !----------------------------------------------------
    ! 5. Fill the configuration state
    !----------------------------------------------------
      bench_config%d_unpb = d_unpb
      bench_config%d_ndfo = d_ndfo
      bench_config%d_rld  = d_rld
      bench_config%d_n    = d_n
      bench_config%d_wnlo = d_wnlo

    !-----------------------------------------------------
    ! 6. Store the namelist for restart
    !-----------------------------------------------------
    IF(my_process_is_stdio())  THEN
      funit = open_tmpfile()
      WRITE(funit,NML=bench_nml)                    
      CALL store_and_close_namelist(funit, 'bench_nml')             
    ENDIF

    !-----------------------------------------------------
    ! 7. write the contents of the namelist to an ASCII file
    !-----------------------------------------------------
    IF(my_process_is_stdio()) WRITE(nnml_output,nml=bench_nml)

  END SUBROUTINE read_bench_namelist

END MODULE mo_bench_nml

