!>
!! Contains the setup of the variables for meteogram output.
!!
!! @par Revision History
!!  by F. Prill, DWD (2011-09-28)
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
MODULE mo_meteogram_nml

  USE mo_kind,               ONLY: wp
  USE mo_parallel_config,    ONLY: nproma
  USE mo_impl_constants,     ONLY: max_dom
  USE mo_exception,          ONLY: finish
  USE mo_io_units,           ONLY: nnml, nnml_output
  USE mo_namelist,           ONLY: position_nml, positioned, open_nml, close_nml
  USE mo_mpi,                ONLY: my_process_is_stdio
  USE mo_master_control,     ONLY: use_restart_namelists
  USE mo_io_restart_namelist,ONLY: open_tmpfile, store_and_close_namelist,   &
                                 & open_and_restore_namelist, close_tmpfile
  USE mo_meteogram_config,   ONLY: t_station_list, t_meteogram_output_config, &
    &                              meteogram_output_config, MAX_NVARS, &
    &                              MAX_NAME_LENGTH, MAX_NUM_STATIONS, FTYPE_NETCDF
  USE mo_nml_annotate,       ONLY: temp_defaults, temp_settings

  IMPLICIT NONE
  PUBLIC :: read_meteogram_namelist

  !-------------------------------------------------------------------------
  ! Namelist variables
  !-------------------------------------------------------------------------

  TYPE t_list_of_stations
    REAL(wp)              :: lat, lon
    CHARACTER (LEN=48)    :: zname
  END TYPE t_list_of_stations

  LOGICAL                         :: lmeteogram_enabled(max_dom)  !> Flag. True, if meteogram output is enabled
  CHARACTER (len=MAX_NAME_LENGTH) :: zprefix(max_dom)         !> string with file name prefix for output file
  INTEGER                         :: ftype(max_dom)           !< file type (NetCDF, ...)
  LOGICAL                         :: ldistributed(max_dom)    !< Flag. Separate files for each PE
  INTEGER                         :: n0_mtgrm(max_dom)        !> intitial time step for meteogram output
  INTEGER                         :: ninc_mtgrm(max_dom)      !> output interval (in time steps)
  LOGICAL                         :: loutput_tiles            !> activate tile specific output

  ! same for all patches:
  TYPE(t_list_of_stations) :: stationlist_tot(MAX_NUM_STATIONS)   !> list of meteogram stations

  ! Positive-list of variables (optional). Only variables contained in
  ! this list are included in this meteogram. If the default list is
  ! not changed by user input, then all available variables are added
  ! to the meteogram
  CHARACTER(len=MAX_NAME_LENGTH)    :: var_list(MAX_NVARS)

  !> Namelist for meteogram output
  NAMELIST/meteogram_output_nml/ lmeteogram_enabled, zprefix, ldistributed,            &
    &                            loutput_tiles, n0_mtgrm, ninc_mtgrm, stationlist_tot, &
    &                            var_list

CONTAINS
  !>
  !! Read Namelist for meteogram output.
  !!
  !! This subroutine
  !! - reads the Namelist for meteogram output
  !! - sets default values
  !! - potentially overwrites the defaults by values used in a
  !!   previous integration (if this is a resumed run)
  !! - reads the user's (new) specifications
  !! - stores the Namelist for restart
  !! - fills the configuration state (partly)
  !!
  !! @par Revision History
  !!  by F. Prill, DWD (2011-09-28)
  !!
  SUBROUTINE read_meteogram_namelist( filename )

    CHARACTER(LEN=*), INTENT(IN)   :: filename
    ! local variables
    CHARACTER(len=*), PARAMETER ::  &
      &  routine = 'mo_mtgrm_nml::read_meteogram_namelist'

    INTEGER                        :: istat, funit, idom, istation, &
      &                               jb, jc, nblks, npromz, nstations, idx
    INTEGER                        :: iunit

    !-----------------------
    ! 1. default settings
    !-----------------------

    lmeteogram_enabled(:)    =          .FALSE.
    zprefix(:)               =     "METEOGRAM_"
    ftype(:)                 =     FTYPE_NETCDF
    ldistributed(:)          =           .TRUE.
    loutput_tiles            =          .FALSE.
    n0_mtgrm(:)              =               0
    ninc_mtgrm(:)            =               1
    var_list(:)              =             " "
    stationlist_tot(:)%lon   = 0._wp
    stationlist_tot(:)%lat   = 0._wp
    stationlist_tot(:)%zname = ""
    ! To determine the length of the provided list, we set the default
    ! value of the longitude to "-1."
    stationlist_tot(:)%lon   = -1._wp
    stationlist_tot(1)%zname = "Hamburg"
    stationlist_tot(1)%lon   = 53.633_wp
    stationlist_tot(1)%lat   =  9.983_wp

    !------------------------------------------------------------------
    ! 2. If this is a resumed integration, overwrite the defaults above
    !    by values used in the previous integration.
    !------------------------------------------------------------------
    IF (use_restart_namelists()) THEN
      funit = open_and_restore_namelist('meteogram_output_nml')
      READ(funit,NML=meteogram_output_nml)
      CALL close_tmpfile(funit)
    END IF

    !-------------------------------------------------------------------------
    ! 3. Read user's (new) specifications (Done so far by all MPI processes)
    !-------------------------------------------------------------------------
    CALL open_nml(TRIM(filename))
    CALL position_nml ('meteogram_output_nml', status=istat)
    IF (my_process_is_stdio()) THEN
      iunit = temp_defaults()
      WRITE(iunit, meteogram_output_nml)   ! write defaults to temporary text file
    END IF
    SELECT CASE (istat)
    CASE (POSITIONED)
      READ (nnml, meteogram_output_nml)                                       ! overwrite default settings
      IF (my_process_is_stdio()) THEN
        iunit = temp_settings()
        WRITE(iunit, meteogram_output_nml)   ! write settings to temporary text file
      END IF
    END SELECT
    CALL close_nml

    !----------------------------------------------------
    ! 4. Fill the configuration state
    !----------------------------------------------------

    ! determine length of station list:
    nstations = 0
    DO
      nstations = nstations + 1
      IF (nstations > MAX_NUM_STATIONS) EXIT
      IF (stationlist_tot(nstations)%lon == -1._wp) EXIT
    END DO
    nstations = nstations - 1

    ! fill in values for each model domain:
    DO idom=1,max_dom
      meteogram_output_config(idom)%lenabled     = lmeteogram_enabled(idom)
      meteogram_output_config(idom)%zprefix      = zprefix(idom)
      meteogram_output_config(idom)%ftype        = ftype(idom)
#ifdef NOMPI
      meteogram_output_config(idom)%ldistributed = .TRUE.
#else
      meteogram_output_config(idom)%ldistributed = ldistributed(idom)
#endif
      meteogram_output_config(idom)%loutput_tiles= loutput_tiles
      meteogram_output_config(idom)%n0_mtgrm     = n0_mtgrm(idom)
      meteogram_output_config(idom)%ninc_mtgrm   = ninc_mtgrm(idom)
      meteogram_output_config(idom)%nstations    = nstations
      meteogram_output_config(idom)%var_list     = var_list

      nblks   = nstations/nproma + 1
      npromz  = nstations - nproma*(nblks-1)
      meteogram_output_config(idom)%nblks   = nblks
      meteogram_output_config(idom)%npromz  = npromz 

      ALLOCATE(meteogram_output_config(idom)%station_list(nproma,nblks))

      jc = 0
      jb = 1
      DO istation=1,nstations
        jc = jc + 1
        IF (jc>nproma) THEN
          jc = 1
          jb = jb + 1
        END IF
        meteogram_output_config(idom)%station_list(jc,jb)%zname        = &
          &  stationlist_tot(istation)%zname
        meteogram_output_config(idom)%station_list(jc,jb)%location%lat = &
          &  stationlist_tot(istation)%lat
        meteogram_output_config(idom)%station_list(jc,jb)%location%lon = &
          &  stationlist_tot(istation)%lon
      END DO

    END DO

    ! consistency check
    idx = 0
    DO idom=1,max_dom
      IF (meteogram_output_config(idom)%lenabled) THEN 
        idx = idom
        EXIT
      END IF
    END DO
    DO idom=(idx+1),max_dom
      IF (.NOT. meteogram_output_config(idom)%lenabled) CYCLE
      IF (meteogram_output_config(idom)%ldistributed .NEQV. meteogram_output_config(idx)%ldistributed) THEN
        CALL finish( TRIM(routine), "Inconsistent namelist setting for domains (ldistributed)!")
      END IF
    END DO

    !-----------------------------------------------------
    ! 5. Store the namelist for restart
    !-----------------------------------------------------
    IF(my_process_is_stdio())  THEN
      funit = open_tmpfile()
      WRITE(funit,NML=meteogram_output_nml)
      CALL store_and_close_namelist(funit, 'meteogram_output_nml')
    ENDIF

    !-----------------------------------------------------
    ! 6. write the contents of the namelist to an ASCII file
    !-----------------------------------------------------
    IF(my_process_is_stdio()) THEN
      WRITE(nnml_output,nml=meteogram_output_nml)
    END IF

  END SUBROUTINE read_meteogram_namelist

END MODULE mo_meteogram_nml
