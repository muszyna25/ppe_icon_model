!>
!! Contains the setup of the variables for io.
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
!!
MODULE mo_io_nml
!-------------------------------------------------------------------------
!
!    ProTeX FORTRAN source: Style 2
!    modified for ICON project, DWD/MPI-M 2006
!
!-------------------------------------------------------------------------
  USE mo_kind,               ONLY: wp
  USE mo_impl_constants,     ONLY: max_char_length, max_ntracer, max_dom, &
    &                              PRES_MSL_METHOD_GME, RH_METHOD_WMO
  USE mo_io_units,           ONLY: nnml, nnml_output, filename_max
  USE mo_namelist,           ONLY: position_nml, positioned, open_nml, close_nml
  USE mo_mpi,                ONLY: my_process_is_stdio, p_n_work
  USE mo_master_control,     ONLY: use_restart_namelists
  USE mo_io_restart_namelist,ONLY: open_tmpfile, store_and_close_namelist,   &
                                 & open_and_restore_namelist, close_tmpfile
  USE mo_nml_annotate,       ONLY: temp_defaults, temp_settings
  USE mo_io_config,          ONLY: config_lkeep_in_sync           => lkeep_in_sync          , &
                                 & config_dt_diag                 => dt_diag                , &
                                 & config_dt_checkpoint           => dt_checkpoint          , &
                                 & config_inextra_2d              => inextra_2d             , &
                                 & config_inextra_3d              => inextra_3d             , &
                                 & config_lflux_avg               => lflux_avg              , &
                                 & config_itype_pres_msl          => itype_pres_msl         , &
                                 & config_output_nml_dict         => output_nml_dict        , &
                                 & config_netcdf_dict             => netcdf_dict            , &
                                 & config_lnetcdf_flt64_output    => lnetcdf_flt64_output   , &
                                 & config_itype_rh                => itype_rh               , &
                                 & config_restart_file_type       => restart_file_type      , &
                                 & config_write_initial_state     => write_initial_state    , &
                                 & config_write_last_restart      => write_last_restart     , &
                                 & config_timeSteps_per_outputStep  => timeSteps_per_outputStep ,&
                                 & config_lmask_boundary            => lmask_boundary

  USE mo_exception,        ONLY: finish


  IMPLICIT NONE
  PUBLIC :: read_io_namelist

  
CONTAINS
  !>
  !! Read Namelist for I/O.
  !!
  !! This subroutine
  !! - reads the Namelist for I/O
  !! - sets default values
  !! - potentially overwrites the defaults by values used in a
  !!   previous integration (if this is a resumed run)
  !! - reads the user's (new) specifications
  !! - stores the Namelist for restart
  !! - fills the configuration state (partly)
  !!
  !! @par Revision History
  !!  by Daniel Reinert, DWD (2011-06-07)
  !!
  SUBROUTINE read_io_namelist( filename )

    CHARACTER(LEN=*), INTENT(IN)   :: filename
    INTEGER                        :: istat, funit
    INTEGER                        :: iunit

    !-------------------------------------------------------------------------
    ! Namelist variables
    !-------------------------------------------------------------------------

    LOGICAL :: lkeep_in_sync              ! if .true., sync stream after each timestep
    REAL(wp):: dt_diag                    ! diagnostic output timestep [seconds]
    REAL(wp):: dt_checkpoint              ! timestep [seconds] for triggering new restart file

    INTEGER :: inextra_2d                 ! number of extra output fields for debugging
    INTEGER :: inextra_3d                 ! number of extra output fields for debugging
    LOGICAL :: lflux_avg                  ! if .FALSE. the output fluxes are accumulated
                                          !  from the beginning of the run
                                          ! if .TRUE. the output fluxex are average values
                                          !  from the beginning of the run, except of
                                          !  TOT_PREC that would be accumulated
    INTEGER :: itype_pres_msl             ! Specifies method for computation of mean sea level pressure
                                          ! 1: GME-type extrapolation
                                          ! 2: stepwise analytical integration
                                          ! 3: IFS method
                                          ! 4: IFS method with consistency correction
    INTEGER :: itype_rh                   ! Specifies method for computation of relative humidity
                                          ! 1: WMO: water only (e_s=e_s_water)
                                          ! 2: IFS: mixed phases (e_s=a*e_s_water + b*e_s_ice)


    CHARACTER(LEN=filename_max) :: &
      &        output_nml_dict,    &      !< maps variable names onto the internal ICON names.
      &        netcdf_dict                !< maps internal variable names onto names in output file (NetCDF only).

    LOGICAL :: lnetcdf_flt64_output       !< if .TRUE. floating point valued NetCDF output
                                          !  is written in 64-bit instead of 32-bit accuracy


    INTEGER :: restart_file_type

    LOGICAL :: write_initial_state

    LOGICAL :: write_last_restart

    INTEGER :: timeSteps_per_outputStep

    LOGICAL :: lmask_boundary ! flag: true, if interpolation zone should be masked *in output*
    
    NAMELIST/io_nml/ lkeep_in_sync, dt_diag, dt_checkpoint,  &
      &              inextra_2d, inextra_3d,                 &
      &              lflux_avg, itype_pres_msl, itype_rh,    &
      &              output_nml_dict, netcdf_dict,           &
      &              lnetcdf_flt64_output,                   &
      &              restart_file_type, write_initial_state, &
      &              write_last_restart, timeSteps_per_outputStep, &
      &              lmask_boundary

    !-----------------------
    ! 1. default settings
    !-----------------------
    lkeep_in_sync           = .FALSE.

    dt_diag                 = 86400._wp    !  1 day

    ! Note: The default needs to be empty, since there exists a
    ! concurrent namelist parameter to specify this value:
    dt_checkpoint           = 0._wp  ! unspecified

    inextra_2d              = 0     ! no extra output 2D fields
    inextra_3d              = 0     ! no extra output 3D fields
    lflux_avg               = .TRUE.
    itype_pres_msl          = PRES_MSL_METHOD_GME
    itype_rh                = RH_METHOD_WMO       ! WMO: water only
    output_nml_dict         = ' '
    netcdf_dict             = ' '
    lnetcdf_flt64_output    = .FALSE.

    restart_file_type       = config_restart_file_type
    write_initial_state     = config_write_initial_state
    write_last_restart      = config_write_last_restart
    timeSteps_per_outputStep        = config_timeSteps_per_outputStep

    lmask_boundary          = .FALSE.

    !------------------------------------------------------------------
    ! 2. If this is a resumed integration, overwrite the defaults above
    !    by values used in the previous integration.
    !------------------------------------------------------------------
    IF (use_restart_namelists()) THEN
      funit = open_and_restore_namelist('io_nml')
      READ(funit,NML=io_nml)
      CALL close_tmpfile(funit)
    END IF

    !-------------------------------------------------------------------------
    ! 3. Read user's (new) specifications (Done so far by all MPI processes)
    !-------------------------------------------------------------------------
    CALL open_nml(TRIM(filename))
    CALL position_nml ('io_nml', status=istat)
    IF (my_process_is_stdio()) THEN
      iunit = temp_defaults()
      WRITE(iunit, io_nml)   ! write defaults to temporary text file
    END IF
    SELECT CASE (istat)
    CASE (POSITIONED)
      READ (nnml, io_nml)                                       ! overwrite default settings
      IF (my_process_is_stdio()) THEN
        iunit = temp_settings()
        WRITE(iunit, io_nml)   ! write settings to temporary text file
      END IF
    END SELECT
    CALL close_nml

    !----------------------------------------------------
    ! 4. Fill the configuration state
    !----------------------------------------------------

    config_lkeep_in_sync           = lkeep_in_sync
    config_dt_diag                 = dt_diag
    config_dt_checkpoint           = dt_checkpoint
    config_inextra_2d              = inextra_2d
    config_inextra_3d              = inextra_3d
    config_lflux_avg               = lflux_avg
    config_itype_pres_msl          = itype_pres_msl
    config_itype_rh                = itype_rh
    config_output_nml_dict         = output_nml_dict
    config_netcdf_dict             = netcdf_dict
    config_lnetcdf_flt64_output    = lnetcdf_flt64_output
    config_restart_file_type       = restart_file_type
    config_write_initial_state     = write_initial_state
    config_timeSteps_per_outputStep= timeSteps_per_outputStep
    config_write_last_restart      = write_last_restart
    config_lmask_boundary          = lmask_boundary

    !-----------------------------------------------------
    ! 5. Store the namelist for restart
    !-----------------------------------------------------
    IF(my_process_is_stdio())  THEN
      funit = open_tmpfile()
      WRITE(funit,NML=io_nml)
      CALL store_and_close_namelist(funit, 'io_nml')
    ENDIF
    !-----------------------------------------------------
    ! 6. write the contents of the namelist to an ASCII file
    !-----------------------------------------------------
    IF(my_process_is_stdio()) THEN
      WRITE(nnml_output,nml=io_nml)
    END IF


  END SUBROUTINE read_io_namelist

END MODULE mo_io_nml
