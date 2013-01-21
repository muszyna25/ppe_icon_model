!>
!! Contains the setup of the variables for io.
!!
!! @par Revision History
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
MODULE mo_io_nml
!-------------------------------------------------------------------------
!
!    ProTeX FORTRAN source: Style 2
!    modified for ICON project, DWD/MPI-M 2006
!
!-------------------------------------------------------------------------
!
!
!
!
  USE mo_kind,               ONLY: wp
  USE mo_impl_constants,     ONLY: max_char_length, max_ntracer, max_dom, &
    &                              PRES_MSL_METHOD_GME
  USE mo_io_units,           ONLY: nnml, nnml_output, filename_max
  USE mo_namelist,           ONLY: position_nml, positioned, open_nml, close_nml
  USE mo_mpi,                ONLY: my_process_is_stdio, p_n_work
  USE mo_master_control,     ONLY: is_restart_run
  USE mo_io_restart_namelist,ONLY: open_tmpfile, store_and_close_namelist,   &
                                 & open_and_restore_namelist, close_tmpfile
  USE mo_io_config,          ONLY: config_out_expname             => out_expname            , &
                                 & config_out_filetype            => out_filetype           , &
                                 & config_lkeep_in_sync           => lkeep_in_sync          , &
                                 & config_dt_data                 => dt_data                , &
                                 & config_dt_diag                 => dt_diag                , &
                                 & config_dt_file                 => dt_file                , &
                                 & config_dt_checkpoint           => dt_checkpoint          , &
                                 & config_lwrite_initial          => lwrite_initial         , &
                                 & config_lwrite_oce_timestepping => lwrite_oce_timestepping, &
                                 & config_lwrite_vorticity        => lwrite_vorticity       , &
                                 & config_lwrite_divergence       => lwrite_divergence      , &
                                 & config_lwrite_omega            => lwrite_omega           , &
                                 & config_lwrite_pres             => lwrite_pres            , &
                                 & config_lwrite_z3               => lwrite_z3              , &
                                 & config_lwrite_tracer           => lwrite_tracer          , &
                                 & config_lwrite_tend_phy         => lwrite_tend_phy        , &
                                 & config_lwrite_radiation        => lwrite_radiation       , &
                                 & config_lwrite_precip           => lwrite_precip          , &
                                 & config_lwrite_cloud            => lwrite_cloud           , &
                                 & config_lwrite_tke              => lwrite_tke             , &
                                 & config_lwrite_surface          => lwrite_surface         , &
                                 & config_lwrite_extra            => lwrite_extra           , &
                                 & config_inextra_2d              => inextra_2d             , &
                                 & config_inextra_3d              => inextra_3d             , &
                                 & config_lflux_avg               => lflux_avg              , &
                                 & config_lwrite_dblprec          => lwrite_dblprec         , &
                                 & config_lwrite_decomposition    => lwrite_decomposition   , &
                                 & config_itype_pres_msl          => itype_pres_msl         , &
                                 & config_varnames_map_file       => varnames_map_file      , &
                                 & config_out_varnames_map_file   => out_varnames_map_file

  USE mo_exception,        ONLY: message, message_text, finish
  USE mo_parallel_config,  ONLY: nproma

  IMPLICIT NONE
  PUBLIC :: read_io_namelist
  CHARACTER(len=*), PARAMETER, PRIVATE :: version = '$Id$'

  !-------------------------------------------------------------------------
  ! Namelist variables
  !-------------------------------------------------------------------------

  CHARACTER(len=max_char_length) :: out_expname
  INTEGER :: out_filetype               ! 1 - GRIB1, 2 - netCDF
  LOGICAL :: lkeep_in_sync              ! if .true., sync stream after each timestep
  REAL(wp):: dt_data                    ! output timestep [seconds]
  REAL(wp):: dt_diag                    ! diagnostic output timestep [seconds]
  REAL(wp):: dt_file                    ! timestep [seconds] for triggering new output file
  REAL(wp):: dt_checkpoint              ! timestep [seconds] for triggering new restart file

  LOGICAL :: lwrite_initial             ! if .true., write out initial state
  LOGICAL :: lwrite_oce_timestepping    ! if .true., write out intermediate ocean variables
  LOGICAL :: lwrite_dblprec             ! if .true., write out in double precision
  LOGICAL :: lwrite_decomposition       ! if .true., write field with MPI_RANK
  LOGICAL :: lwrite_vorticity           ! if .true., write out vorticity
  LOGICAL :: lwrite_divergence          ! if .true., write out divergence
  LOGICAL :: lwrite_pres                ! if .true., write out full level pressure
  LOGICAL :: lwrite_tend_phy            ! if .true., write out physics-induced tendencies
  LOGICAL :: lwrite_radiation           ! if .true., write out fields related to radiation
  LOGICAL :: lwrite_precip              ! if .true., write out precip
  LOGICAL :: lwrite_cloud               ! if .true., write out cloud variables
  LOGICAL :: lwrite_z3                  ! if .true., write out geopotential on full levels
  LOGICAL :: lwrite_omega               ! if .true., write out the vertical velocity
                                        ! in pressure coordinate
  LOGICAL :: lwrite_tke                 ! if .true., write out TKE
  LOGICAL :: lwrite_surface             ! if .true., write out surface related fields
  LOGICAL :: lwrite_tracer(max_ntracer) ! for each tracer, if .true. write out
                                        ! tracer on full levels
  LOGICAL :: lwrite_extra               ! if .true., write out extra fields
  INTEGER :: inextra_2d                 ! number of extra output fields for debugging
  INTEGER :: inextra_3d                 ! number of extra output fields for debugging
  LOGICAL :: lflux_avg                  ! if .FALSE. the output fluxes are accumulated 
                                        !  from the beginning of the run
                                        ! if .TRUE. the output fluxex are average values 
                                        !  from the beginning of the run, except of 
                                        !  TOT_PREC that would be accumulated
  INTEGER :: itype_pres_msl             ! Specifies method for computation of mean sea level pressure

  CHARACTER(LEN=filename_max) :: &
    &        varnames_map_file,   &     !< maps variable names onto the internal ICON names.
    &        out_varnames_map_file      !< maps internal variable names onto names in output file (NetCDF only).

  NAMELIST/io_nml/ out_expname, out_filetype, lkeep_in_sync,             &
    &              dt_data, dt_diag, dt_file, dt_checkpoint,             &
    &              lwrite_initial, lwrite_dblprec, lwrite_decomposition, &
    &              lwrite_vorticity, lwrite_divergence, lwrite_omega,    &
    &              lwrite_pres, lwrite_z3, lwrite_tracer,                &
    &              lwrite_tend_phy, lwrite_radiation, lwrite_precip,     &
    &              lwrite_cloud, lwrite_tke, lwrite_surface,             &
    &              lwrite_extra, inextra_2d, inextra_3d,                 &
    &              lflux_avg, lwrite_oce_timestepping, itype_pres_msl,   &
    &              varnames_map_file, out_varnames_map_file
  
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

    !-----------------------
    ! 1. default settings
    !-----------------------
    out_expname             = 'IIIEEEETTTT'
    out_filetype            = 2
    lkeep_in_sync           = .FALSE.

    dt_data                 = 21600.0_wp   !  6 hours
    dt_diag                 = 86400._wp    !  1 day
    dt_file                 = 2592000._wp  ! 30 days
    dt_checkpoint           = 2592000._wp  ! 30 days

    lwrite_initial          = .TRUE.
    lwrite_dblprec          = .FALSE.
    lwrite_decomposition    = .FALSE.
    lwrite_vorticity        = .TRUE.
    lwrite_divergence       = .TRUE.
    lwrite_omega            = .TRUE.
    lwrite_pres             = .TRUE.
    lwrite_z3               = .TRUE.
    lwrite_tracer(:)        = .TRUE.

    lwrite_tend_phy         = .FALSE.
    lwrite_radiation        = .FALSE.
    lwrite_precip           = .FALSE.
    lwrite_cloud            = .FALSE.
    lwrite_tke              = .FALSE.
    lwrite_surface          = .FALSE.
    lwrite_extra            = .FALSE.
    inextra_2d              = 0     ! no extra output 2D fields
    inextra_3d              = 0     ! no extra output 3D fields
    lflux_avg               = .TRUE.
    lwrite_oce_timestepping = .FALSE.
    itype_pres_msl          = PRES_MSL_METHOD_GME
    varnames_map_file       = ' '
    out_varnames_map_file   = ' '

    !------------------------------------------------------------------
    ! 2. If this is a resumed integration, overwrite the defaults above
    !    by values used in the previous integration.
    !------------------------------------------------------------------
    IF (is_restart_run()) THEN
      funit = open_and_restore_namelist('io_nml')
      READ(funit,NML=io_nml)
      CALL close_tmpfile(funit)
    END IF

    !-------------------------------------------------------------------------
    ! 3. Read user's (new) specifications (Done so far by all MPI processes)
    !-------------------------------------------------------------------------
    CALL open_nml(TRIM(filename))
    CALL position_nml ('io_nml', status=istat)
    SELECT CASE (istat)
    CASE (POSITIONED)
      READ (nnml, io_nml)
    END SELECT
    CALL close_nml

    !----------------------------------------------------
    ! 4. Fill the configuration state
    !----------------------------------------------------

    config_out_expname             = out_expname
    config_out_filetype            = out_filetype
    config_lkeep_in_sync           = lkeep_in_sync
    config_dt_data                 = dt_data
    config_dt_diag                 = dt_diag
    config_dt_file                 = dt_file
    config_dt_checkpoint           = dt_checkpoint
    config_lwrite_initial          = lwrite_initial
    config_lwrite_vorticity        = lwrite_vorticity
    config_lwrite_divergence       = lwrite_divergence
    config_lwrite_omega            = lwrite_omega
    config_lwrite_pres             = lwrite_pres
    config_lwrite_z3               = lwrite_z3
    config_lwrite_tracer           = lwrite_tracer
    config_lwrite_tend_phy         = lwrite_tend_phy
    config_lwrite_radiation        = lwrite_radiation
    config_lwrite_precip           = lwrite_precip
    config_lwrite_cloud            = lwrite_cloud
    config_lwrite_tke              = lwrite_tke
    config_lwrite_surface          = lwrite_surface
    config_lwrite_extra            = lwrite_extra
    config_inextra_2d              = inextra_2d
    config_inextra_3d              = inextra_3d
    config_lflux_avg               = lflux_avg
    config_lwrite_dblprec          = lwrite_dblprec
    config_lwrite_decomposition    = lwrite_decomposition
    config_lwrite_oce_timestepping = lwrite_oce_timestepping
    config_itype_pres_msl          = itype_pres_msl
    config_varnames_map_file       = varnames_map_file
    config_out_varnames_map_file   = out_varnames_map_file

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
