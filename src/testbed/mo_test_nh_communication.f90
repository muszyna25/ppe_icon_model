!>
!! @brief Main program for the ICON atmospheric model
!!
!! @author
!!  Leonidas Linardakis (MPI-M)
!!  Hui Wan             (MPI-M)
!!
!! @par Copyright
!! 2002-2011 by DWD and MPI-M
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
MODULE mo_test_nh_communication

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: message, message_text, finish
  USE mo_mpi,                 ONLY: work_mpi_barrier, my_process_is_stdio
  USE mo_timer,               ONLY: init_timer, ltimer, new_timer, timer_start, timer_stop, &
    & print_timer

  USE mo_master_control,      ONLY: get_my_process_name, get_my_model_no
  USE mo_icon_testbed_config, ONLY: testbed_iterations, calculate_iterations

  USE mo_model_domain,        ONLY:  p_patch  
  USE mo_atmo_model,          ONLY: construct_atmo_model, destruct_atmo_model
  
  USE mo_parallel_config,    ONLY: itype_comm, iorder_sendrecv
  USE mo_sync,               ONLY: SYNC_C, SYNC_E, SYNC_V, sync_patch_array, sync_patch_array_mult
!   USE mo_icon_comm_interface,ONLY: construct_icon_communication, destruct_icon_communication
  USE mo_icon_comm_lib
  USE mo_atmo_nonhydrostatic, ONLY: construct_atmo_nonhydrostatic, destruct_atmo_nonhydrostatic

  !-------------------------------------------------------------------------
  ! for the nh run
  USE mo_kind,                ONLY: wp
  USE mo_nonhydro_types,      ONLY: t_nh_state, t_nh_prog, t_nh_diag, t_nh_metrics
  USE mo_nonhydro_state,      ONLY: bufr
  USE mo_nonhydrostatic_config,ONLY: iadv_rcf, lhdiff_rcf, l_nest_rcf, itime_scheme
  USE mo_diffusion_config,     ONLY: diffusion_config
  USE mo_dynamics_config,      ONLY: nnow,nnew, nnow_rcf, nnew_rcf, nsav1, nsav2
  USE mo_io_config,            ONLY: l_outputtime, is_checkpoint_time,&
    &                                istime4output
  USE mo_parallel_config,      ONLY: nproma, itype_comm
  USE mo_run_config,           ONLY: ltestcase, dtime, dtime_adv, nsteps,     &
    &                                ltransport, ntracer, lforcing, iforcing, &
    &                                msg_level
  USE mo_timer,               ONLY: ltimer, timers_level, timer_start, timer_stop,   &
    &                               timer_model_init, timer_nudging,                 &
    &                               timer_bdy_interp, timer_feedback, timer_nesting, &
    &                               timer_integrate_nh, timer_nh_diagnostics
  USE mo_grid_config,         ONLY: global_cell_type
  USE mo_atm_phy_nwp_config,  ONLY: dt_phy, atm_phy_nwp_config
  USE mo_nwp_phy_state,       ONLY: prm_diag, prm_nwp_tend, phy_params
  USE mo_lnd_nwp_config,      ONLY: nlev_soil, nlev_snow
  USE mo_model_domain,        ONLY: t_patch
  USE mo_grid_config,         ONLY: n_dom, lfeedback, ifeedback_type, l_limited_area, &
    &                               n_dom_start, lredgrid_phys
  USE mo_nh_testcases,        ONLY: init_nh_testtopo, init_nh_testcase, nh_test_name, &
    &                               rotate_axis_deg
  USE mo_nh_pa_test,          ONLY: set_nh_w_rho
  USE mo_nh_df_test,          ONLY: get_nh_df_velocity
  USE mo_integrate_density_pa,ONLY: integrate_density_pa
  USE mo_nh_hex_util,         ONLY: forcing_straka, momentum_adv
  USE mo_nh_supervise,        ONLY: supervise_total_integrals_nh
  USE mo_intp_data_strc,      ONLY: t_int_state, t_lon_lat_intp
  USE mo_intp_rbf,            ONLY: rbf_vec_interpol_cell
  USE mo_intp,                ONLY: edges2cells_scalar, verts2edges_scalar, edges2verts_scalar, &
    &                               verts2cells_scalar
  USE mo_grf_intp_data_strc,  ONLY: t_gridref_state
  USE mo_grf_bdyintp,         ONLY: interpol_scal_grf
  USE mo_nh_nest_utilities,   ONLY: compute_tendencies, boundary_interpolation,    &
                                    complete_nesting_setup, prep_bdy_nudging,      &
                                    outer_boundary_nudging, nest_boundary_nudging, &
                                    prep_rho_bdy_nudging, density_boundary_nudging
  USE mo_nh_feedback,         ONLY: feedback, relax_feedback
  USE mo_datetime,            ONLY: t_datetime, print_datetime, add_time
  USE mo_timer,               ONLY: timer_total, timer_start, timer_stop
  USE mo_output,              ONLY: init_output_files, write_output,  &
    &                               create_restart_file
  USE mo_io_restart,          ONLY: write_restart_info_file
  USE mo_exception,           ONLY: message, message_text, finish
  USE mo_impl_constants,      ONLY: SUCCESS, MAX_CHAR_LENGTH, iphysproc,    &
    &                               iphysproc_short, itconv, itccov, itrad, &
    &                               itradheat, itsso, itsatad, itgwd, inwp, &
    &                               itupdate, itturb, itgscp, itsfc, min_rlcell_int, &
                                    min_rledge_int
  USE mo_divergent_modes,     ONLY: divergent_modes_5band
  USE mo_math_divrot,         ONLY: div_avg, div
  USE mo_solve_nonhydro,      ONLY: solve_nh
  USE mo_advection_stepping,  ONLY: step_advection
  USE mo_nh_dtp_interface,    ONLY: prepare_tracer
  USE mo_nh_diffusion,        ONLY: diffusion_tria, diffusion_hex
  USE mo_mpi,                 ONLY: my_process_is_stdio, my_process_is_mpi_parallel, &
    &                               proc_split, push_glob_comm, pop_glob_comm
#ifdef NOMPI
  USE mo_mpi,                 ONLY: my_process_is_mpi_all_seq
#endif
  
  USE mo_sync,                ONLY: sync_patch_array_mult, &
                                    global_max, &
                                    SYNC_C, SYNC_E, sync_patch_array
  USE mo_nh_interface_nwp,    ONLY: nwp_nh_interface
  USE mo_phys_nest_utilities, ONLY: interpol_phys_grf, feedback_phys_diag, interpol_rrg_grf
  USE mo_vertical_grid,       ONLY: set_nh_metrics
  USE mo_nh_diagnose_pres_temp,ONLY: diagnose_pres_temp
  USE mo_nh_held_suarez_interface, ONLY: held_suarez_nh_interface
  USE mo_vertical_coord_table,ONLY: vct
  USE mo_master_control,      ONLY: is_restart_run
  USE mo_io_restart_attributes,ONLY: get_restart_attribute

!   USE mo_nwp_mpiomp_rrtm_interface, ONLY: nwp_start_radiation_ompthread, model_end_ompthread, &
!     & init_ompthread_radiation
  USE mo_parallel_config,     ONLY: parallel_radiation_omp, nh_stepping_ompthreads
  USE mo_name_list_output,    ONLY: write_name_list_output, istime4name_list_output, &
    &                               output_file
  !-------------------------------------------------------------------------

!-------------------------------------------------------------------------
IMPLICIT NONE
PRIVATE

PUBLIC :: test_nh_communication

CONTAINS
  

  !-------------------------------------------------------------------------
  !>
  !!
  SUBROUTINE test_nh_communication(namelist_filename,shr_namelist_filename)

    CHARACTER(LEN=*), INTENT(in) :: namelist_filename
    CHARACTER(LEN=*), INTENT(in) :: shr_namelist_filename

    ! 2D variables    

    
    INTEGER :: patch_no, i

    
    CHARACTER(*), PARAMETER :: method_name = "mo_test_communication:test_nh_communication"

    !---------------------------------------------------------------------
    write(0,*) TRIM(get_my_process_name()), ': Start of ', method_name
    
    !---------------------------------------------------------------------

    ltimer = .false.
    CALL construct_atmo_model(namelist_filename,shr_namelist_filename)
    CALL construct_atmo_nonhydrostatic()
        
    CALL work_mpi_barrier()
    
    !---------------------------------------------------------------------
    DO i=1,testbed_iterations
!       CALL integrate_nh_test_comm(p_nh_state, p_patch, p_int_state, datetime, p_grf_state, &
!         &               1, jstep, dtime, dtime_adv, sim_time, 1,                 &
!         &               l_compute_diagnostic_quants                              )
    ENDDO
    !---------------------------------------------------------------------
    

    !---------------------------------------------------------------------
    ! Carry out the shared clean-up processes
    !---------------------------------------------------------------------
    CALL destruct_atmo_nonhydrostatic()
    CALL destruct_atmo_model()
!     CALL destruct_icon_communication()
    CALL message(TRIM(method_name),'clean-up finished')

    !---------------------------------------------------------------------
    ! print the timers
!    IF (my_process_is_stdio()) THEN
      CALL message("===================", "=======================")
      WRITE(message_text,*) "Communication Iterations=", testbed_iterations
      CALL message(method_name, TRIM(message_text))
      CALL print_timer()
!    ENDIF
    !---------------------------------------------------------------------
     

  END SUBROUTINE test_nh_communication
  !-------------------------------------------------------------------------
  

  !-------------------------------------------------------------------------
#ifdef __NO_NESTING__
   SUBROUTINE integrate_nh_test_comm (p_nh_state, p_patch, p_int_state, datetime,  &
  &        p_grf_state, jg, nstep_global, dt_loc, dtadv_loc, sim_time,   &
  &        num_steps, l_compute_diagnostic_quants)
#else
  RECURSIVE SUBROUTINE integrate_nh_test_comm (p_nh_state, p_patch, p_int_state, datetime,  &
  &        p_grf_state, jg, nstep_global, dt_loc, dtadv_loc, sim_time, &
  &        num_steps, l_compute_diagnostic_quants)
#endif

    CHARACTER(len=*), PARAMETER ::  &
      &  routine = 'mo_nh_stepping:integrate_nh'

    TYPE(t_patch), TARGET, INTENT(in)    :: p_patch(n_dom_start:n_dom)    !< patch
    TYPE(t_int_state),TARGET,INTENT(in)  :: p_int_state(n_dom_start:n_dom)!< interpolation state
    TYPE(t_datetime), INTENT(in)         :: datetime
    TYPE(t_nh_state), TARGET, INTENT(inout) :: p_nh_state(n_dom) !< nonhydrostatic state
    TYPE(t_gridref_state), INTENT(INOUT) :: p_grf_state(n_dom_start:n_dom)!< gridref state

    INTEGER , INTENT(IN)    :: jg           !< current grid level
    INTEGER , INTENT(IN)    :: nstep_global !< counter of global time step
    INTEGER , INTENT(IN)    :: num_steps    !< number of time steps to be executed
    REAL(wp), INTENT(IN)    :: dt_loc       !< time step applicable to local grid level
    REAL(wp), INTENT(IN)    :: dtadv_loc    !< advective time step applicable to 
                                            !< local grid level
    REAL(wp), INTENT(INOUT) :: sim_time(n_dom) !< elapsed simulation time on each
                                               !< grid level
    LOGICAL, INTENT(IN) :: l_compute_diagnostic_quants    !< computation of diagnostic quantities

    ! Local variables

    ! Time levels
    INTEGER :: n_now_grf, n_now, n_new, n_save, n_temp
    INTEGER :: n_now_rcf, n_new_rcf, n_upt_rcf  ! accounts for reduced calling frequencies (rcf)

    INTEGER :: jstep, jgp, jgc, jn
    INTEGER :: nsteps_nest ! number of time steps executed in nested domain

    REAL(wp):: dt_sub, dtadv_sub ! (advective) timestep for next finer grid level
    REAL(wp):: rdt_loc,  rdtadv_loc ! inverse time step for local grid level

    REAL(wp), DIMENSION(:,:,:), POINTER  :: p_vn   => NULL()
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE  :: z_tmp_e

    LOGICAL, PARAMETER :: l_straka=.FALSE.
    LOGICAL :: l_predictor
    LOGICAL :: l_bdy_nudge
    LOGICAL :: linit_vertnest(2)
    LOGICAL :: l_recompute
    LOGICAL :: lclean_mflx   ! for reduced calling freqency: determines whether
                             ! mass-fluxes and trajectory-velocities are reset to zero
                             ! i.e. for starting new integration sweep
    LOGICAL :: lcall_hdiff

    ! Switch to determine manner of OpenMP parallelization in interpol_scal_grf
!     LOGICAL :: lpar_fields=.FALSE.

    ! Switch to determine if nested domains are called at a given time step
    LOGICAL :: l_call_nests = .FALSE.

!$  INTEGER :: num_threads_omp, omp_get_max_threads

   END SUBROUTINE integrate_nh_test_comm 
  !-------------------------------------------------------------------------

END MODULE mo_test_nh_communication

