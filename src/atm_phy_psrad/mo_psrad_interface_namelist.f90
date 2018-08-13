!>
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
!! @author L. Linardakis, MPI-M, Hamburg
!!
!! @remarks
!!  
MODULE mo_psrad_interface_namelist

  USE mo_impl_constants, ONLY: max_dom
  USE mo_exception,      ONLY: finish, message, warning, message_text
  USE mo_mpi,            ONLY: my_process_is_stdio
  USE mo_run_nml,        ONLY: read_run_namelist
  USE mo_io_nml         ,ONLY: read_io_namelist
  USE mo_dbg_nml        ,ONLY: read_dbg_namelist
  USE mo_grid_nml       ,ONLY: read_grid_namelist
  USE mo_name_list_output_init,ONLY: read_name_list_output_namelists
  USE mo_parallel_nml   ,ONLY: read_parallel_namelist
  USE mo_time_nml       ,ONLY: read_time_namelist
  USE mo_namelist,       ONLY: position_nml, positioned, open_nml, close_nml, &
    & open_nml_output, close_nml_output, POSITIONED
  USE mo_nml_annotate   ,ONLY: log_nml_settings
  USE mo_io_units,       ONLY: nnml, nnml_output
  USE mo_nml_annotate,   ONLY: temp_defaults, temp_settings
 
!   USE mo_radiation_nml,  ONLY: read_radiation_namelist
  USE mo_echam_cld_nml  ,ONLY: process_echam_cld_nml
  USE mo_echam_rad_nml  ,ONLY: process_echam_rad_nml
  USE mo_echam_phy_nml,  ONLY: process_echam_phy_nml

!   USE mo_psrad_radiation_parameters, ONLY : rad_perm ! read by namelist here

  USE mo_echam_phy_config, ONLY: eval_echam_phy_config, eval_echam_phy_tc
  USE mo_time_config,    ONLY: time_config, dt_restart
  USE mo_run_config,     ONLY: nlev, num_lev, configure_run
  USE mo_grid_config,    ONLY: init_grid_configuration, n_dom
  USE mo_time_management,ONLY: compute_timestep_settings,                        &
    &                          compute_restart_settings,                         &
    &                          compute_date_settings
  USE mo_event_manager,   ONLY: initEventManager

  USE mo_psrad_interface, ONLY: setup_psrad_radiation
  USE mo_load_restart,    ONLY: read_restart_header
  USE mo_master_config,   ONLY: isRestart

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: configure_ps_radiation, configure_ps_radiation_test
  PUBLIC :: number_of_levels

  INTEGER :: number_of_levels(max_dom) = -1

  INTEGER :: nsteps = -1
  !--------------------------------------------------------------------------

CONTAINS

  !---------------------------------------------------------------------
  !>
  SUBROUTINE configure_ps_radiation(ps_rad_namelist_filename,shr_namelist_filename)
    CHARACTER(LEN=*), INTENT(in) :: ps_rad_namelist_filename
    CHARACTER(LEN=*), INTENT(in) :: shr_namelist_filename

    IF (isRestart()) THEN
      CALL message('configure_ps_radiation','Read restart file meta data ...')
      CALL read_restart_header("atm") ! we read the same file as the atmosphere
    ENDIF

    CALL read_ps_radiation_all_namelists(ps_rad_namelist_filename,shr_namelist_filename)

    !---------------------------------------------------------------------
    ! 2. Call configure_run to finish filling the run_config state.
    !    This needs to be done very early (but anyway after atm_crosscheck)
    !    because some component of the state, e.g., num_lev, may be
    !    modified in this subroutine which affect the following CALLs.
    !---------------------------------------------------------------------
    CALL configure_run

    CALL init_grid_configuration()    ! so that the number of grids is known

    !--------------------------------------------------------------------
    ! Compute date/time/time step settings
    !--------------------------------------------------------------------
    !
    ! Note that the ordering of the following three calls must not be
    ! changed, since they rely on previous results:
    !
    CALL compute_timestep_settings()
    CALL compute_restart_settings()
    CALL compute_date_settings("atm", dt_restart, nsteps)
    CALL initEventManager(time_config%tc_exp_refdate)

    CALL  eval_echam_phy_config
    CALL  eval_echam_phy_tc

    number_of_levels = num_lev
    
  END SUBROUTINE configure_ps_radiation

  !---------------------------------------------------------------------
  !>
  SUBROUTINE read_ps_radiation_all_namelists(ps_rad_namelist_filename,shr_namelist_filename)

    CHARACTER(LEN=*), INTENT(in) :: ps_rad_namelist_filename
    CHARACTER(LEN=*), INTENT(in) :: shr_namelist_filename

    !-----------------------------------------------------------------
    ! Create a new file in which all the namelist variables and their
    ! actual values used in the model run will be stored.
    !-----------------------------------------------------------------

    IF(my_process_is_stdio()) CALL open_nml_output('NAMELIST_ICON_output_radiation')

    !-----------------------------------------------------------------
    ! Read namelists that are shared by all components of the model.
    ! This means that the same namelists with the same values are
    ! read by all components of a coupled system.
    !-----------------------------------------------------------------

    CALL read_time_namelist           (TRIM(shr_namelist_filename))

    !-----------------------------------------------------------------
    ! Read namelists that are specific to the oce model.
    ! In case of a coupled simulation, the atmosphere model may also
    ! read some of these namelists, but probably from a different
    ! ASCII file containing different values.
    !-----------------------------------------------------------------

    ! General
    !
    CALL read_parallel_namelist       (TRIM(ps_rad_namelist_filename))
    CALL read_run_namelist            (TRIM(ps_rad_namelist_filename))
    CALL read_io_namelist             (TRIM(ps_rad_namelist_filename))
    CALL read_name_list_output_namelists (TRIM(ps_rad_namelist_filename))
    CALL read_dbg_namelist            (TRIM(ps_rad_namelist_filename))

    ! Grid
    CALL read_grid_namelist           (TRIM(ps_rad_namelist_filename))

    !    
    CALL process_echam_phy_nml        (TRIM(ps_rad_namelist_filename)) 
    CALL process_echam_cld_nml        (TRIM(ps_rad_namelist_filename))
    CALL process_echam_rad_nml        (TRIM(ps_rad_namelist_filename))

    ! needs to be done after process_echam_cld_nml
    CALL setup_psrad_radiation(TRIM(ps_rad_namelist_filename))
   !-------------------------------------------------------------------


    IF (my_process_is_stdio()) CALL close_nml_output

    ! write an annotate table of all namelist settings to a text file
    IF (my_process_is_stdio()) CALL log_nml_settings("nml.ps_radiation.log")

  END SUBROUTINE read_ps_radiation_all_namelists
  !-------------------------------------------------------------------------


  !---------------------------------------------------------------------
  !>
  SUBROUTINE configure_ps_radiation_test(ps_rad_namelist_filename,shr_namelist_filename)
    CHARACTER(LEN=*), INTENT(in) :: ps_rad_namelist_filename
    CHARACTER(LEN=*), INTENT(in) :: shr_namelist_filename

!     CALL read_ps_radiation_all_namelists(ps_rad_namelist_filename,shr_namelist_filename)

    !---------------------------------------------------------------------
    ! 2. Call configure_run to finish filling the run_config state.
    !    This needs to be done very early (but anyway after atm_crosscheck)
    !    because some component of the state, e.g., num_lev, may be
    !    modified in this subroutine which affect the following CALLs.
    !---------------------------------------------------------------------
!     CALL configure_run
! 
!     CALL init_grid_configuration()    ! so that the number of grids is known

    !--------------------------------------------------------------------
    ! Compute date/time/time step settings
    !--------------------------------------------------------------------
    !
    ! Note that the ordering of the following three calls must not be
    ! changed, since they rely on previous results:
    !
!     CALL compute_timestep_settings()
!     CALL compute_restart_settings()
!     CALL compute_date_settings("atm", dt_restart, nsteps)

    number_of_levels = num_lev
    
    write(0,*) "================================================="
    write(0,*) " nsteps=", nsteps
    write(0,*) " number_of_levels=", number_of_levels
    write(0,*) "================================================="

  END SUBROUTINE configure_ps_radiation_test

END MODULE mo_psrad_interface_namelist
