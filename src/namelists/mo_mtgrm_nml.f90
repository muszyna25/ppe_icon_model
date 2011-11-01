!>
!! Contains the setup of the variables for meteogram output.
!!
!! @par Revision History
!!  by F. Prill, DWD (2011-09-28)
!!
!! @par Copyright
!! 2002-2006 by DWD and MPI-M
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
MODULE mo_meteogram_nml

  USE mo_kind,               ONLY: wp
  USE mo_parallel_config,    ONLY: nproma
  USE mo_impl_constants,     ONLY: max_dom
  USE mo_io_units,           ONLY: nnml, nnml_output
  USE mo_namelist,           ONLY: position_nml, positioned, open_nml, close_nml
  USE mo_mpi,                ONLY: my_process_is_stdio
  USE mo_master_control,     ONLY: is_restart_run
  USE mo_io_restart_namelist,ONLY: open_tmpfile, store_and_close_namelist,   &
                                 & open_and_restore_namelist, close_tmpfile
  USE mo_meteogram_config,   ONLY: t_station_list, t_meteogram_output_config, &
    &                              meteogram_output_config, &
    &                              MAX_NAME_LENGTH, MAX_NUM_STATIONS, FTYPE_NETCDF

  IMPLICIT NONE
  PUBLIC :: read_meteogram_namelist
  CHARACTER(len=*), PARAMETER, PRIVATE :: version = '$Id$'

  !-------------------------------------------------------------------------
  ! Namelist variables
  !-------------------------------------------------------------------------

  TYPE t_list_of_stations
    REAL(wp)              :: lat, lon
    CHARACTER (LEN=20)    :: zname
  END TYPE t_list_of_stations

  LOGICAL                         :: lmeteogram_enabled(max_dom)  !> Flag. True, if meteogram output is enabled
  CHARACTER (len=MAX_NAME_LENGTH) :: zprefix(max_dom)         !> string with file name prefix for output file
  INTEGER                         :: ftype(max_dom)           !< file type (NetCDF, ...)
  LOGICAL                         :: ldistributed(max_dom)    !< Flag. Separate files for each PE
  INTEGER                         :: n0_mtgrm(max_dom)        !> intitial time step for meteogram output
  INTEGER                         :: ninc_mtgrm(max_dom)      !> output interval (in time steps)

  ! same for all patches:
  TYPE(t_list_of_stations) :: stationlist_tot(MAX_NUM_STATIONS)   !> list of meteogram stations

  !> Namelist for meteogram output
  NAMELIST/meteogram_output_nml/ lmeteogram_enabled, zprefix, ldistributed, &
    &                            n0_mtgrm, ninc_mtgrm, stationlist_tot

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
    INTEGER                        :: istat, funit, idom, istation, &
      &                               jb, jc, nblks, npromz, nstations

    !-----------------------
    ! 1. default settings
    !-----------------------

    lmeteogram_enabled(:)    =          .FALSE.
    zprefix(:)               =     "METEOGRAM_"
    ftype(:)                 =     FTYPE_NETCDF
    ldistributed(:)          =           .TRUE.
    n0_mtgrm(:)              =               1
    ninc_mtgrm(:)            =               1
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
    IF (is_restart_run()) THEN
      funit = open_and_restore_namelist('meteogram_output_nml')
      READ(funit,NML=meteogram_output_nml)
      CALL close_tmpfile(funit)
    END IF

    !-------------------------------------------------------------------------
    ! 3. Read user's (new) specifications (Done so far by all MPI processes)
    !-------------------------------------------------------------------------
    CALL open_nml(TRIM(filename))
    CALL position_nml ('meteogram_output_nml', status=istat)
    SELECT CASE (istat)
    CASE (POSITIONED)
      READ (nnml, meteogram_output_nml)
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
      meteogram_output_config(idom)%n0_mtgrm     = n0_mtgrm(idom)
      meteogram_output_config(idom)%ninc_mtgrm   = ninc_mtgrm(idom)
      meteogram_output_config(idom)%nstations    = nstations

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
