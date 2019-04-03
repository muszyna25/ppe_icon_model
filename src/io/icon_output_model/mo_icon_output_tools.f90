!>
!! @page pagecontrolmodelf901 Main program for the ICON icon_output model
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
MODULE mo_icon_output_tools

  USE mo_exception,           ONLY: message, finish
  USE mo_master_config,       ONLY: isRestart
  USE mo_parallel_config,     ONLY: p_test_run, l_test_openmp, num_io_procs, &
       &                            pio_type, num_test_pe, num_prefetch_proc
  USE mo_mpi,                 ONLY: set_mpi_work_communicators, process_mpi_io_size, &
       &                            stop_mpi, my_process_is_io, my_process_is_mpi_test,   &
       &                            set_mpi_work_communicators, process_mpi_io_size
  USE mo_impl_constants,      ONLY: pio_type_async
#ifdef HAVE_CDI_PIO
  USE mo_mpi,                 ONLY: mpi_comm_null, p_comm_work_io
  USE mo_impl_constants,      ONLY: pio_type_cdipio
  USE yaxt,                   ONLY: xt_initialize, xt_initialized
  USE mo_cdi,                 ONLY: namespacegetactive
  USE mo_cdi_pio_interface,   ONLY: nml_io_cdi_pio_namespace, &
    &                                   cdi_base_namespace, &
    &                                   nml_io_cdi_pio_client_comm, &
    &                                   nml_io_cdi_pio_conf_handle
  USE mo_name_list_output_init, ONLY: init_cdipio_cb
  USE mo_name_list_output,    ONLY: write_ready_files_cdipio
  USE mo_impl_constants,      ONLY: pio_type_cdipio
  USE mo_cdi,                 ONLY: namespaceGetActive, namespaceSetActive
  USE mo_cdi_pio_interface,   ONLY: nml_io_cdi_pio_namespace
#endif
  USE mo_timer,               ONLY: init_timer, timer_start, timer_stop, print_timer, &
       &                            timer_model_init
  USE mo_name_list_output_init, ONLY: init_name_list_output, parse_variable_groups, &
    &                                 create_vertical_axes, output_file
  USE mo_derived_variable_handling, ONLY: init_statistics_streams, finish_statistics_streams
  USE mo_name_list_output,    ONLY: close_name_list_output, name_list_io_main_proc
  USE mo_name_list_output_config,  ONLY: use_async_name_list_io
  USE mo_level_selection, ONLY: create_mipz_level_selections
  USE mo_zaxis_type,          ONLY: zaxisTypeList, t_zaxisTypeList

  !  USE mo_advection_config,    ONLY: configure_advection
  USE mo_run_config,          ONLY: output_mode, dtime
  USE mo_gribout_config,      ONLY: configure_gribout

  ! Control parameters: run control, dynamics, i/o
  !


  ! Horizontal grid
  USE mo_model_domain,        ONLY: t_patch_3d, p_patch_local_parent
  !
  USE mo_grid_config,         ONLY: n_dom, use_dummy_cell_closure

  USE mo_build_decomposition, ONLY: build_decomposition
  USE mo_complete_subdivision,ONLY: setup_phys_patches

  USE mo_util_dbg_prnt,       ONLY: init_dbg_index
  USE mo_impl_constants,      ONLY: success

  USE mo_alloc_patches,        ONLY: destruct_patches, destruct_comm_patterns
  USE mo_icon_output_read_namelists, ONLY: read_icon_output_namelists
  USE mo_load_restart,         ONLY: read_restart_header, read_restart_files
  USE mo_restart_nml_and_att,  ONLY: getAttributesForRestarting
  USE mo_icon_comm_interface,  ONLY: construct_icon_communication, destruct_icon_communication
  USE mo_output_event_types,   ONLY: t_sim_step_info
  USE mo_io_config,            ONLY: restartWritingParameters, write_initial_state
  USE mo_key_value_store,     ONLY: t_key_value_store
  !-------------------------------------------------------------
 
  IMPLICIT NONE

  PRIVATE
#ifdef HAVE_CDI_PIO
  INCLUDE 'cdipio.inc'
#endif

    PUBLIC :: init_io_processes
    PUBLIC :: prepare_output

  CONTAINS
  !-------------------------------------------------------------------------
  SUBROUTINE init_io_processes()
    USE mo_time_config,         ONLY: time_config

    TYPE(t_sim_step_info)   :: sim_step_info
    INTEGER                 :: jstep0
    TYPE(t_key_value_store), POINTER :: restartAttributes
    CHARACTER(LEN=*), PARAMETER :: &
      & method_name = 'mo_ocean_model:init_io_processes'
    
    IF (process_mpi_io_size < 1) THEN
      IF (output_mode%l_nml) THEN
        ! -----------------------------------------
        ! non-asynchronous I/O (performed by PE #0)
        ! -----------------------------------------
        CALL message(method_name,'synchronous namelist I/O scheme is enabled.')
      ENDIF
      ! nothing to do
      RETURN
    ENDIF
!..pa ++ .AND. pio_type == pio_type_async
!    IF (process_mpi_io_size > 0 .AND. pio_type == pio_type_async) THEN
!..pa
!
!    ! Decide whether async vlist or name_list IO is to be used,
!    ! only one of both may be enabled!
!
!     IF (output_mode%l_nml) THEN
!      ! -----------------------------------------
!      ! asynchronous I/O
!      ! -----------------------------------------
!      !
!      use_async_name_list_io = .TRUE.
!      CALL message(method_name,'asynchronous namelist I/O scheme is enabled.')
!      ! consistency check
!      IF (my_process_is_io() .AND. (.NOT. my_process_is_mpi_test())) THEN
!
!        ! compute sim_start, sim_end
!        CALL datetimeToString(time_config%tc_exp_startdate, sim_step_info%sim_start)
!        CALL datetimeToString(time_config%tc_exp_stopdate, sim_step_info%sim_end)
!        CALL datetimeToString(time_config%tc_startdate, sim_step_info%run_start)
!        CALL datetimeToString(time_config%tc_stopdate, sim_step_info%restart_time)
!
!        sim_step_info%dtime      = dtime
!        jstep0 = 0
!
!        restartAttributes => getAttributesForRestarting()
!        IF (ASSOCIATED(restartAttributes)) THEN
!
!          ! get start counter for time loop from restart file:
!          jstep0 = restartAttributes%getInteger("jstep")
!        END IF
!        sim_step_info%jstep0    = jstep0
!!         CALL name_list_io_main_proc(sim_step_info, isample=1)
!!pa
!        write(0,*)"Before name_list_io_main_proc"
!        CALL name_list_io_main_proc(sim_step_info)
!        write(0,*)"After name_list_io_main_proc"
!!pa
!      END IF
!     ELSE IF (my_process_is_io() .AND. (.NOT. my_process_is_mpi_test())) THEN
!      ! Shut down MPI
!        CALL stop_mpi
!        STOP
!     ENDIF
!    ENDIF

#ifdef HAVE_CDI_PIO
    
    IF (process_mpi_io_size > 0 .AND. pio_type == pio_type_cdipio) THEN

      CALL message(method_name,'Collective asynchronous namelist I/O scheme is enabled.')
      IF (my_process_is_io() .AND. (.NOT. my_process_is_mpi_test())) THEN

        ! compute sim_start, sim_end
        sim_step_info%sim_start = time_config%tc_exp_startdate
        sim_step_info%sim_end = time_config%tc_exp_stopdate
        sim_step_info%run_start = time_config%tc_startdate
        sim_step_info%restart_time = time_config%tc_stopdate

        sim_step_info%dtime      = dtime
        jstep0 = 0

        CALL getAttributesForRestarting(restartAttributes)
        ! get start counter for time loop from restart file:
        IF (ASSOCIATED(restartAttributes)) CALL restartAttributes%get("jstep", jstep0)
        sim_step_info%jstep0    = jstep0
        CALL init_statistics_streams
      ENDIF

     IF (.NOT. xt_initialized()) CALL xt_initialize(p_comm_work_io)


      cdi_base_namespace = namespaceGetActive()
!pa      write(0,*)"Before 1st cdiPioConfSetCallBackActions"
      CALL cdiPioConfSetCallBackActions(nml_io_cdi_pio_conf_handle, &
        cdipio_callback_postcommsetup, init_cdipio_cb)
!pa      write(0,*)"After 1st cdiPioConfSetCallBackActions"
!pa      write(0,*)"Before 2nd cdiPioConfSetCallBackActions"
      CALL cdiPioConfSetCallBackActions(nml_io_cdi_pio_conf_handle, &
        cdipio_callback_postwritebatch, write_ready_files_cdipio)
!pa      write(0,*)"After 2nd cdiPioConfSetCallBackActions"
!pa      write(0,*)"Before cdiPioInit"
      nml_io_cdi_pio_client_comm = &
        &   cdiPioInit(p_comm_work_io, nml_io_cdi_pio_conf_handle, &
        &              nml_io_cdi_pio_namespace)
!pa      write(0,*)"After cdiPioInit"
      IF (nml_io_cdi_pio_client_comm == mpi_comm_null) THEN
        ! todo: terminate program cleanly here
!pa      write(0,*)"init_io_processes: nml_io_cdi_pio_client_comm=", &
!pa      &          nml_io_cdi_pio_client_comm
!pa      write(0,*)"init_io_processes: p_comm_work_io=", &
!pa      &          p_comm_work_io
!pa      write(0,*)"init_io_processes: nml_io_cdi_pio_namespace", &
!pa      &          nml_io_cdi_pio_namespace
        CALL stop_mpi
        STOP
      END IF
    ELSE IF (my_process_is_io() .AND. (.NOT. my_process_is_mpi_test())) THEN
      ! Shut down MPI
      CALL stop_mpi
      STOP
    ENDIF
#else
      CALL finish(method_name, 'CDI-PIO requested but unavailable')
#endif


  END SUBROUTINE init_io_processes
  !-------------------------------------------------------------------

  !--------------------------------------------------------------------------
  SUBROUTINE prepare_output()
    USE mo_time_config,         ONLY: time_config

    CHARACTER(*), PARAMETER :: method_name = "mo_ocean_model:prepare_output"

    TYPE(t_sim_step_info)               :: sim_step_info
    INTEGER                             :: jstep0
    TYPE(t_key_value_store), POINTER :: restartAttributes

    !------------------------------------------------------------------
    ! Initialize output file if necessary;
    ! Write out initial conditions.
    !------------------------------------------------------------------

    IF (output_mode%l_nml) THEN
!       WRITE(0,*)'process_mpi_io_size:',process_mpi_io_size
!       IF (process_mpi_io_size > 0) use_async_name_list_io = .TRUE.
      CALL parse_variable_groups()
      ! compute sim_start, sim_end
      sim_step_info%sim_start = time_config%tc_exp_startdate
      sim_step_info%sim_end = time_config%tc_exp_stopdate
      sim_step_info%run_start = time_config%tc_startdate
      sim_step_info%restart_time = time_config%tc_stopdate

      sim_step_info%dtime      = dtime
      jstep0 = 0

      CALL getAttributesForRestarting(restartAttributes)
      ! get start counter for time loop from restart file:
      IF (ASSOCIATED(restartAttributes)) CALL restartAttributes%get("jstep", jstep0)

      sim_step_info%jstep0    = jstep0
      CALL init_statistics_streams
      CALL init_name_list_output(sim_step_info, opt_lprintlist=.TRUE.,opt_l_is_ocean=.TRUE.)
      CALL create_mipz_level_selections(output_file)
      CALL create_vertical_axes(output_file)
    ENDIF
  
  END SUBROUTINE prepare_output
  !--------------------------------------------------------------------------

END MODULE mo_icon_output_tools

