!>
!! Contains the main stepping routine the 3-dim hydrostatic ocean model.
!!
!! @author Peter Korn, Stephan Lorenz, MPI
!!
!! @par Revision History
!!  Initial version by Stephan Lorenz (MPI-M), (2010-04).
!!   - renaming and adjustment of hydrostatic ocean model V1.0.3 to ocean domain and patch_oce
!!  Modification by Stephan Lorenz, MPI-M, 2010-10
!!   - new module mo_hydro_ocean_run including updated reconstructions
!
!
!! @par Copyright
!! 2002-2007 by DWD and MPI-M
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
MODULE mo_hydro_ocean_run
  !-------------------------------------------------------------------------
  USE mo_kind,                   ONLY: wp
  USE mo_impl_constants,         ONLY: max_char_length
  USE mo_model_domain,           ONLY: t_patch, t_patch_3d, t_subset_range, t_patch_vert
  USE mo_grid_config,            ONLY: n_dom
  USE mo_grid_subset,            ONLY: get_index_range
  USE mo_sync,                   ONLY: sync_patch_array, sync_e!, sync_c, sync_v
  USE mo_ocean_nml,              ONLY: iswm_oce, n_zlev, no_tracer, &
    & diagnostics_level, &
    & eos_type, i_sea_ice, l_staggered_timestep, gibraltar,l_time_marching
  USE mo_dynamics_config,        ONLY: nold, nnew
  USE mo_io_config,              ONLY: n_checkpoints
  USE mo_run_config,             ONLY: nsteps, dtime, ltimer, output_mode
  USE mo_exception,              ONLY: message, message_text, finish
  USE mo_ext_data_types,         ONLY: t_external_data
  !USE mo_io_units,               ONLY: filename_max
  USE mo_datetime,               ONLY: t_datetime, print_datetime, add_time, datetime_to_string
  USE mo_timer,                  ONLY: timer_start, timer_stop, timer_total, timer_solve_ab,  &
    & timer_tracer_ab, timer_vert_veloc, timer_normal_veloc, &
    & timer_upd_phys, timer_upd_flx
  USE mo_oce_ab_timestepping,    ONLY: solve_free_surface_eq_ab, &
    & calc_normal_velocity_ab,  &
    & calc_vert_velocity,       &
    & update_time_indices
 USE mo_oce_types,              ONLY: t_hydro_ocean_state, t_hydro_ocean_acc, t_hydro_ocean_diag, &
    & t_hydro_ocean_prog
 USE mo_oce_state,              ONLY: destruct_hydro_ocean_state,            &
    & ocean_restart_list
  USE mo_ocean_initialization,   ONLY: set_lateral_boundary_values
  USE mo_oce_math_operators,     ONLY: calc_thickness
  USE mo_operator_ocean_coeff_3d,ONLY: t_operator_coeff, update_diffusion_matrices
  USE mo_scalar_product,         ONLY: calc_scalar_product_veloc_3d
  USE mo_oce_tracer,             ONLY: advect_tracer_ab
  USE mo_io_restart,             ONLY: write_restart_info_file, create_restart_file
  USE mo_oce_bulk,               ONLY: update_sfcflx
  USE mo_oce_forcing,            ONLY: destruct_ocean_forcing
  USE mo_sea_ice,                ONLY: destruct_atmos_for_ocean,&
    & destruct_atmos_fluxes,&
    & destruct_sea_ice,  &
    & update_ice_statistic, compute_mean_ice_statistics, reset_ice_statistics
  USE mo_sea_ice_types,          ONLY: t_sfc_flx, t_atmos_fluxes, t_atmos_for_ocean, &
    & t_sea_ice
  USE mo_oce_physics,            ONLY: t_ho_params, &
    & destruct_ho_params, update_ho_params
  USE mo_oce_thermodyn,          ONLY: calc_density_mpiom_func, calc_density_lin_eos_func,&
    & calc_density_jmdwfg06_eos_func, calc_potential_density, &
    & calc_density
  USE mo_name_list_output,       ONLY: write_name_list_output, istime4name_list_output
  USE mo_oce_diagnostics,        ONLY: calc_slow_oce_diagnostics, calc_fast_oce_diagnostics, &
    & construct_oce_diagnostics,&
    & destruct_oce_diagnostics, t_oce_timeseries, &
    & calc_moc, calc_psi
  USE mo_oce_ab_timestepping_mimetic, ONLY: init_ho_lhs_fields_mimetic
  USE mo_linked_list,            ONLY: t_list_element, find_list_element
  USE mo_var_list,               ONLY: print_var_list
  USE mo_io_restart_attributes,  ONLY: get_restart_attribute
  USE mo_mpi,                    ONLY: my_process_is_stdio
  USE mo_time_config,            ONLY: time_config
  USE mo_master_control,         ONLY: is_restart_run
  USE mo_statistics
  USE mo_sea_ice_nml,            ONLY: i_ice_dyn
  USE mo_util_dbg_prnt,          ONLY: dbg_print
  
  IMPLICIT NONE
  
  PRIVATE
  
  !VERSION CONTROL:
  CHARACTER(LEN=*), PARAMETER :: version = '$Id$'
  
  ! public interface
  !
  ! public subroutines
  PUBLIC :: perform_ho_stepping
  PUBLIC :: prepare_ho_stepping
  PUBLIC :: finalise_ho_integration
  PRIVATE:: update_intermediate_tracer_vars
  !
  CHARACTER(LEN=12)  :: str_module = 'HYDRO-ocerun'  ! Output of module for 1 line debug
  INTEGER :: idt_src                      ! Level of detail for 1 line debug
  !
  !-------------------------------------------------------------------------
  
CONTAINS
  SUBROUTINE prepare_ho_stepping(patch_3d, operators_coefficients, ocean_state, p_phys_param, is_restart)
    TYPE(t_patch_3d ), INTENT(in)     :: patch_3d
    TYPE(t_operator_coeff)            :: operators_coefficients
    TYPE(t_hydro_ocean_state), TARGET :: ocean_state
    TYPE (t_ho_params)                :: p_phys_param
    LOGICAL, INTENT(in)               :: is_restart
    
    IF (is_restart) THEN
      ! Prepare ocean_state%p_prog, since it is needed by the dynamical sea ice model
      IF ( i_sea_ice > 0 .AND. i_ice_dyn == 1 )         &
        & CALL calc_scalar_product_veloc_3d( patch_3d,  &
        & ocean_state%p_prog(nnew(1))%vn,             &
        & ocean_state%p_prog(nnew(1))%vn,             &
        & ocean_state%p_diag,                         &
        & operators_coefficients)
    ELSE
    ENDIF

    CALL update_diffusion_matrices( patch_3d,         &
      & p_phys_param,                 &
      & operators_coefficients%matrix_vert_diff_e,&
      & operators_coefficients%matrix_vert_diff_c)
  END SUBROUTINE prepare_ho_stepping
  
  !-------------------------------------------------------------------------
  !>
  !! Main stepping routine for call of hydrostatic ocean model
  !!
  !! @par Revision History
  !! Developed by Peter Korn, MPI-M  (2008-2010).
  !! Initial release by Stephan Lorenz, MPI-M (2010-07)
  !
  !
  SUBROUTINE perform_ho_stepping( patch_3d, ocean_state, p_ext_data,          &
    & datetime, lwrite_restart,            &
    & p_sfc_flx, p_phys_param,             &
    & p_as, p_atm_f, p_ice,operators_coefficients)
    
    TYPE(t_patch_3d ),TARGET, INTENT(inout)          :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET, INTENT(inout) :: ocean_state(n_dom)
    TYPE(t_external_data), TARGET, INTENT(in)        :: p_ext_data(n_dom)
    TYPE(t_datetime), INTENT(inout)                  :: datetime
    LOGICAL, INTENT(in)                              :: lwrite_restart
    TYPE(t_sfc_flx)                                  :: p_sfc_flx
    TYPE (t_ho_params)                               :: p_phys_param
    TYPE(t_atmos_for_ocean),  INTENT(inout)          :: p_as
    TYPE(t_atmos_fluxes ),    INTENT(inout)          :: p_atm_f
    TYPE (t_sea_ice),         INTENT(inout)          :: p_ice
    TYPE(t_operator_coeff),   INTENT(inout)          :: operators_coefficients
    
    ! local variables
    INTEGER :: jstep, jg, jtrc
    INTEGER :: nsteps_since_last_output
    INTEGER :: ocean_statistics
    !LOGICAL                         :: l_outputtime
    CHARACTER(LEN=32)               :: datestring, plaindatestring
    TYPE(t_oce_timeseries), POINTER :: oce_ts
    TYPE(t_patch), POINTER :: patch_2d
    TYPE(t_patch_vert), POINTER :: patch_1d
    INTEGER, POINTER :: dolic(:,:)
    REAL(wp), POINTER :: prism_thickness(:,:,:)
    INTEGER :: jstep0 ! start counter for time loop
    
    !CHARACTER(LEN=filename_max)  :: outputfile, gridfile
    CHARACTER(LEN=max_char_length), PARAMETER :: &
      & routine = 'mo_hydro_ocean_run:perform_ho_stepping'
    !------------------------------------------------------------------
    
    nsteps_since_last_output = 1
    CALL init_ho_lhs_fields_mimetic   ( patch_3d )
    
    !------------------------------------------------------------------
    ! no grid refinement allowed here so far
    !------------------------------------------------------------------
    IF (n_dom > 1 ) THEN
      CALL finish(TRIM(routine), ' N_DOM > 1 is not allowed')
    END IF
    jg = n_dom
    
    patch_2d => patch_3d%p_patch_2d(jg)
    
    CALL datetime_to_string(datestring, datetime)
    
    IF (diagnostics_level == 1) &
      & CALL construct_oce_diagnostics( patch_3d, ocean_state(jg), oce_ts, datestring)
    
    ! IF (ltimer) CALL timer_start(timer_total)
    CALL timer_start(timer_total)
    
    time_config%sim_time(:) = 0.0_wp
    
    !------------------------------------------------------------------
    ocean_statistics = new_statistic()
    
    jstep0 = 0
    IF (is_restart_run() .AND. .NOT. time_config%is_relative_time) THEN
      ! get start counter for time loop from restart file:
      CALL get_restart_attribute("jstep", jstep0)
    END IF
    !------------------------------------------------------------------
    ! write initial
    !------------------------------------------------------------------
    IF (output_mode%l_nml) THEN
      ! in general nml output is writen based on the nnew status of the
      ! prognostics variables. Unfortunately, the initialization has to be written
      ! to the nold state. That's why the following manual copying is nec.
      IF (.NOT. is_restart_run()) THEN
        ocean_state(jg)%p_prog(nnew(1))%tracer = ocean_state(jg)%p_prog(nold(1))%tracer
        ! copy old tracer values to spot value fields for propper initial timestep
        ! output
        IF(no_tracer>=1)THEN
          ocean_state(jg)%p_diag%t = ocean_state(jg)%p_prog(nold(1))%tracer(:,:,:,1)
        ENDIF
        IF(no_tracer>=2)THEN
          ocean_state(jg)%p_diag%s = ocean_state(jg)%p_prog(nold(1))%tracer(:,:,:,2)
        ENDIF
        ocean_state(jg)%p_diag%h = ocean_state(jg)%p_prog(nold(1))%h
        CALL calc_potential_density( patch_3d,                     &
          & ocean_state(jg)%p_prog(nold(1))%tracer,&
          & ocean_state(jg)%p_diag%rhopot )
        CALL calc_density( patch_3d,                        &
          & ocean_state(jg)%p_prog(nold(1))%tracer, &
          & ocean_state(jg)%p_diag%rho )
        
        CALL update_ocean_statistics(ocean_state(1),                              &
          & p_sfc_flx,                            &
          & patch_3d%p_patch_2d(1)%cells%owned,   &
          & patch_3d%p_patch_2d(1)%edges%owned,   &
          & patch_3d%p_patch_2d(1)%verts%owned,   &
          & n_zlev)
        IF (i_sea_ice >= 1) CALL update_ice_statistic(p_ice%acc, p_ice,patch_3d%p_patch_2d(1)%cells%owned)
      ENDIF
      
      CALL write_name_list_output(jstep=jstep0)
      
      IF (.NOT. is_restart_run()) THEN
        CALL reset_ocean_statistics(ocean_state(1)%p_acc,p_sfc_flx)
        IF (i_sea_ice >= 1) CALL reset_ice_statistics(p_ice%acc)
      ENDIF
      
    ENDIF ! output_mode%l_nml
    !------------------------------------------------------------------
    ! call the dynamical core: start the time loop
    !------------------------------------------------------------------
    IF(.NOT.l_time_marching)THEN

      !IF(itestcase_oce==28)THEN
      DO jstep = (jstep0+1), (jstep0+nsteps)
      
        CALL datetime_to_string(datestring, datetime)
        WRITE(message_text,'(a,i10,2a)') '  Begin of timestep =',jstep,'  datetime:  ', datestring
        CALL message (TRIM(routine), message_text)
 
!          IF(jstep==1)THEN
!          ocean_state(jg)%p_diag%vn_time_weighted = ocean_state(jg)%p_prog(nold(1))%vn
!          ocean_state(jg)%p_prog(nnew(1))%vn = ocean_state(jg)%p_prog(nold(1))%vn
!          ocean_state(jg)%p_diag%w        =  0.0_wp!0.0833_wp!0.025_wp
!          ocean_state(jg)%p_diag%w(:,:,:) = -0.0833_wp!0.025_wp
!          ENDIF

        !CALL calc_thickness( patch_3d, ocean_state(jg), p_ext_data(jg))
        !CALL calc_vert_velocity(patch_3d, ocean_state(jg),operators_coefficients)
        CALL advect_tracer_ab( patch_3d, ocean_state(jg),  &
          & p_phys_param,p_sfc_flx,&
          & operators_coefficients,&
          & jstep)
        ! One integration cycle finished on the lowest grid level (coarsest
        ! resolution). Set model time.
        CALL add_time(dtime,0,0,0,datetime)
      
        ! Not nice, but the name list output requires this
        time_config%sim_time(1) = time_config%sim_time(1) + dtime
      
        ! update accumulated vars
        CALL update_ocean_statistics(ocean_state(1),&
        & p_sfc_flx,                                &
        & patch_3d%p_patch_2d(1)%cells%owned,       &
        & patch_3d%p_patch_2d(1)%edges%owned,       &
        & patch_3d%p_patch_2d(1)%verts%owned,       &
        & n_zlev)
          

        IF (istime4name_list_output(jstep)) THEN
          IF (diagnostics_level == 1 ) THEN
            CALL calc_slow_oce_diagnostics( patch_3d,      &
            & ocean_state(jg),      &
            & p_sfc_flx,     &
            & p_ice,         &
            & jstep-jstep0,  &
            & datetime,      &
            & oce_ts)
                    
          ENDIF
        
          CALL compute_mean_ocean_statistics(ocean_state(1)%p_acc,p_sfc_flx,nsteps_since_last_output)
          CALL compute_mean_ice_statistics(p_ice%acc,nsteps_since_last_output)
        
          ! set the output variable pointer to the correct timelevel
          CALL set_output_pointers(nnew(1), ocean_state(jg)%p_diag, ocean_state(jg)%p_prog(nnew(1)))
        
          IF (output_mode%l_nml) THEN
            CALL write_name_list_output(jstep)
          ENDIF
        
          CALL message (TRIM(routine),'Write output at:')
          CALL print_datetime(datetime)
        
          ! reset accumulation vars
          CALL reset_ocean_statistics(ocean_state(1)%p_acc,p_sfc_flx,nsteps_since_last_output)
          IF (i_sea_ice >= 1) CALL reset_ice_statistics(p_ice%acc)
        
        END IF
      
        ! Shift time indices for the next loop
        ! this HAS to ge into the restart files, because the start with the following loop
        CALL update_time_indices(jg)
        ! update intermediate timestepping variables for the tracers
        CALL update_intermediate_tracer_vars(ocean_state(jg))
      
        ! write a restart or checkpoint file
        IF (MOD(jstep,n_checkpoints())==0 .OR. ((jstep==(jstep0+nsteps)) .AND. lwrite_restart)) THEN
          CALL create_restart_file( patch=patch_2d,        &
            & datetime=datetime,                           &
            & jstep=jstep,                                 &
            & model_type="oce",                            &
            & opt_depth=n_zlev,                            &
            & opt_sim_time=time_config%sim_time(1),        &
            & opt_nice_class=1)
          ! Create the master (meta) file in ASCII format which contains
          ! info about which files should be read in for a restart run.
          CALL write_restart_info_file
        END IF
      
        nsteps_since_last_output = nsteps_since_last_output + 1
          
      END DO
    ELSEIF(l_time_marching)THEN
    
      time_loop: DO jstep = (jstep0+1), (jstep0+nsteps)
      
        CALL datetime_to_string(datestring, datetime)
        WRITE(message_text,'(a,i10,2a)') '  Begin of timestep =',jstep,'  datetime:  ', datestring
        CALL message (TRIM(routine), message_text)
      
          !In case of a time-varying forcing:
          IF (ltimer) CALL timer_start(timer_upd_flx)
          CALL update_sfcflx( patch_3d, ocean_state(jg), p_as, p_ice, p_atm_f, p_sfc_flx, &
          & jstep, datetime, operators_coefficients)
        
          IF(.NOT.l_staggered_timestep)THEN
          
            CALL calc_thickness( patch_3d, ocean_state(jg), p_ext_data(jg))
          
            CALL set_lateral_boundary_values( patch_3d, ocean_state(jg)%p_prog(nold(1))%vn)
            CALL sync_patch_array(sync_e, patch_3d%p_patch_2d(jg), ocean_state(jg)%p_prog(nold(1))%vn)
          
            CALL calc_scalar_product_veloc_3d( patch_3d,  &
            & ocean_state(jg)%p_prog(nold(1))%vn,         &
            & ocean_state(jg)%p_prog(nold(1))%vn,         &
            & ocean_state(jg)%p_diag,                     &
            & operators_coefficients)
          
            ! activate for calc_scalar_product_veloc_3D
            !---------DEBUG DIAGNOSTICS-------------------------------------------
            idt_src=2  ! output print level (1-5, fix)
            CALL dbg_print('on entry: h-old'           ,ocean_state(jg)%p_prog(nold(1))%h ,str_module,idt_src, &
            & patch_2d%cells%owned )
            idt_src=1  ! output print level (1-5, fix)
            CALL dbg_print('on entry: h-new'           ,ocean_state(jg)%p_prog(nnew(1))%h ,str_module,idt_src, &
            & patch_2d%cells%owned )
            idt_src=3  ! output print level (1-5, fix)
            CALL dbg_print('HydOce: ScaProdVel kin'    ,ocean_state(jg)%p_diag%kin        ,str_module,idt_src, &
            & patch_2d%cells%owned )
            CALL dbg_print('HydOce: ScaProdVel ptp_vn' ,ocean_state(jg)%p_diag%ptp_vn     ,str_module,idt_src, &
            & patch_2d%edges%owned )
            !---------------------------------------------------------------------          
          IF (ltimer) CALL timer_stop(timer_upd_flx)
          
          
          IF (ltimer) CALL timer_start(timer_upd_phys)        
          SELECT CASE (eos_type)
          CASE(1)
            CALL update_ho_params(patch_3d, ocean_state(jg), p_sfc_flx, p_phys_param,&
            & calc_density_lin_eos_func)
          CASE(2)
            CALL update_ho_params(patch_3d, ocean_state(jg), p_sfc_flx, p_phys_param,&
            & calc_density_mpiom_func)
          CASE(3)
            CALL update_ho_params(patch_3d,ocean_state(jg), p_sfc_flx, p_phys_param,&
            & calc_density_jmdwfg06_eos_func)
          CASE default
          END SELECT
          IF (ltimer) CALL timer_stop(timer_upd_phys)
        
          CALL update_diffusion_matrices( patch_3d,   &
          & p_phys_param,                             &
          & operators_coefficients%matrix_vert_diff_e,&
          & operators_coefficients%matrix_vert_diff_c)
        
          !------------------------------------------------------------------------
          ! solve for new free surface
          IF (ltimer) CALL timer_start(timer_solve_ab)
          CALL solve_free_surface_eq_ab (patch_3d, ocean_state(jg), p_ext_data(jg), &
          & p_sfc_flx, p_phys_param, jstep, operators_coefficients)!, p_int(jg))
          IF (ltimer) CALL timer_stop(timer_solve_ab)
        
          !------------------------------------------------------------------------
          ! Step 4: calculate final normal velocity from predicted horizontal
          ! velocity vn_pred and updated surface height
          IF (ltimer) CALL timer_start(timer_normal_veloc)
          CALL calc_normal_velocity_ab(patch_3d, ocean_state(jg),&
          & operators_coefficients, p_ext_data(jg), p_phys_param)
          IF (ltimer) CALL timer_stop(timer_normal_veloc)
        
          !------------------------------------------------------------------------
          ! Step 5: calculate vertical velocity from continuity equation under
          ! incompressiblity condition in the non-shallow-water case
          IF ( iswm_oce /= 1 ) THEN
            IF (ltimer) CALL timer_start(timer_vert_veloc)
            CALL calc_vert_velocity( patch_3d, ocean_state(jg),operators_coefficients)
            IF (ltimer) CALL timer_stop(timer_vert_veloc)
          ENDIF
        
          !------------------------------------------------------------------------
          ! Step 6: transport tracers and diffuse them
          IF (no_tracer>=1) THEN
            IF (ltimer) CALL timer_start(timer_tracer_ab)
            CALL advect_tracer_ab( patch_3d, ocean_state(jg), p_phys_param,&
            & p_sfc_flx,&
            & operators_coefficients,&
            & jstep)
            IF (ltimer) CALL timer_stop(timer_tracer_ab)
          ENDIF
        
        ENDIF
      
        ! One integration cycle finished on the lowest grid level (coarsest
        ! resolution). Set model time.
        CALL add_time(dtime,0,0,0,datetime)
      
        ! Not nice, but the name list output requires this
        time_config%sim_time(1) = time_config%sim_time(1) + dtime
      
        ! perform accumulation for special variables
        IF (no_tracer>=1) THEN
          CALL calc_potential_density( patch_3d,                            &
          & ocean_state(jg)%p_prog(nold(1))%tracer,                         &
          & ocean_state(jg)%p_diag%rhopot )
          CALL calc_psi (patch_2d,patch_3d, ocean_state(jg)%p_diag%u(:,:,:),&
          & ocean_state(jg)%p_prog(nold(1))%h(:,:),                         &
          & ocean_state(jg)%p_diag%u_vint, datetime)
        ENDIF 
        ! update accumulated vars
        CALL update_ocean_statistics(ocean_state(1),&
        & p_sfc_flx,                                &
        & patch_3d%p_patch_2d(1)%cells%owned,       &
        & patch_3d%p_patch_2d(1)%edges%owned,       &
        & patch_3d%p_patch_2d(1)%verts%owned,       &
        & n_zlev)
        IF (i_sea_ice >= 1) CALL update_ice_statistic(p_ice%acc,p_ice,patch_3d%p_patch_2d(1)%cells%owned)
      
        dolic           => patch_3d%p_patch_1d(1)%dolic_c
        prism_thickness => patch_3d%p_patch_1d(1)%prism_thick_c
        CALL calc_fast_oce_diagnostics( patch_3d%p_patch_2d(1),      &
        & dolic, &
        & prism_thickness, &
        & patch_3d%p_patch_1d(1)%zlev_m, &
        & ocean_state(jg)%p_diag)
write(0,*)'istime4name_list_output',jstep,istime4name_list_output(jstep)                
        IF (istime4name_list_output(jstep).OR.jstep>0) THEN
          IF (diagnostics_level == 1 ) THEN
            CALL calc_slow_oce_diagnostics( patch_3d,      &
            & ocean_state(jg),      &
            & p_sfc_flx,     &
            & p_ice,         &
            & jstep-jstep0,  &
            & datetime,      &
            & oce_ts)
            IF (no_tracer>=2) THEN
              CALL calc_moc (patch_2d,patch_3d, ocean_state(jg)%p_diag%w(:,:,:), datetime)
            ENDIF
          ENDIF
          ! compute mean values for output interval
          !TODO [ram] src/io/shared/mo_output_event_types.f90 for types to use
          !TODO [ram] nsteps_since_last_output =
          !TODO [ram] output_event%event_step(output_event%i_event_step)%i_sim_step - output_event%event_step(output_event%i_event_step-1)%i_sim_step
        
          CALL compute_mean_ocean_statistics(ocean_state(1)%p_acc,p_sfc_flx,nsteps_since_last_output)
          CALL compute_mean_ice_statistics(p_ice%acc,nsteps_since_last_output)
        
          ! set the output variable pointer to the correct timelevel
          CALL set_output_pointers(nnew(1), ocean_state(jg)%p_diag, ocean_state(jg)%p_prog(nnew(1)))
        
          IF (output_mode%l_nml) THEN
            CALL write_name_list_output(jstep)
          ENDIF
        
          CALL message (TRIM(routine),'Write output at:')
          CALL print_datetime(datetime)
        
          ! reset accumulation vars
          CALL reset_ocean_statistics(ocean_state(1)%p_acc,p_sfc_flx,nsteps_since_last_output)
          IF (i_sea_ice >= 1) CALL reset_ice_statistics(p_ice%acc)
        
        END IF
      
        ! Shift time indices for the next loop
        ! this HAS to ge into the restart files, because the start with the following loop
        CALL update_time_indices(jg)
        ! update intermediate timestepping variables for the tracers
        CALL update_intermediate_tracer_vars(ocean_state(jg))
      
        ! write a restart or checkpoint file
        IF (MOD(jstep,n_checkpoints())==0 .OR. ((jstep==(jstep0+nsteps)) .AND. lwrite_restart)) THEN
          CALL create_restart_file( patch = patch_2d,       &
            & datetime=datetime,                            &
            & jstep=jstep,                                  &
            & model_type="oce",                             &
            & opt_depth=n_zlev,                             &
            & opt_sim_time=time_config%sim_time(1),         &
            & opt_nice_class=1)
          ! Create the master (meta) file in ASCII format which contains
          ! info about which files should be read in for a restart run.
          CALL write_restart_info_file
        END IF
      
        nsteps_since_last_output = nsteps_since_last_output + 1
      ENDDO time_loop
      
    ENDIF!(l_no_time_marching)THEN   
    
    IF (diagnostics_level==1) CALL destruct_oce_diagnostics(oce_ts)
    CALL delete_statistic(ocean_statistics)
    
    CALL timer_stop(timer_total)
    
  END SUBROUTINE perform_ho_stepping
  !-------------------------------------------------------------------------
  
  
  !-------------------------------------------------------------------------
  !>
  !! Simple routine for finalising integration of hydrostatic ocean model.
  !!
  !! Simple routine for finalising integration of hydrostatic ocean model.
  !! Calls basic routines ...
  !!
  !!
  !! @par Revision History
  !! Initial release by Stephan Lorenz, MPI-M (2010-07)
  !
  !
  SUBROUTINE finalise_ho_integration(ocean_state, p_phys_param, p_as, p_atm_f, p_ice, p_sfc_flx)
    TYPE(t_hydro_ocean_state), INTENT(inout) :: ocean_state(n_dom)
    TYPE (t_ho_params),        INTENT(inout) :: p_phys_param
    TYPE(t_atmos_for_ocean),   INTENT(inout) :: p_as
    TYPE(t_atmos_fluxes ),     INTENT(inout) :: p_atm_f
    TYPE (t_sea_ice),          INTENT(inout) :: p_ice
    TYPE(t_sfc_flx),           INTENT(inout) :: p_sfc_flx
    
    
    !------------------------------------------------------------------
    ! destruct ocean physics and forcing
    ! destruct ocean state is in control_model
    !------------------------------------------------------------------
    CALL destruct_hydro_ocean_state(ocean_state)
    !CALL destruct_hydro_ocean_base(v_base)
    CALL destruct_ho_params(p_phys_param)
    
    IF(no_tracer>0) CALL destruct_ocean_forcing(p_sfc_flx)
    CALL destruct_sea_ice(p_ice)
    CALL destruct_atmos_for_ocean(p_as)
    CALL destruct_atmos_fluxes(p_atm_f)
    
    
  END SUBROUTINE finalise_ho_integration
  
  SUBROUTINE update_intermediate_tracer_vars(ocean_state)
    TYPE(t_hydro_ocean_state), INTENT(inout) :: ocean_state
    
    ! velocity
    ocean_state%p_aux%g_nm1 = ocean_state%p_aux%g_n
    ocean_state%p_aux%g_n   = 0.0_wp
  END SUBROUTINE update_intermediate_tracer_vars
  
  SUBROUTINE update_ocean_statistics(ocean_state,p_sfc_flx,cells,edges,verts,max_zlev)
    TYPE(t_hydro_ocean_state), INTENT(inout) :: ocean_state
    TYPE(t_sfc_flx),           INTENT(inout) :: p_sfc_flx
    TYPE(t_subset_range),      INTENT(in)    :: cells,edges,verts
    INTEGER, INTENT(in)                      :: max_zlev
    
    INTEGER :: jtrc,i
    
    
    ! update ocean state accumulated values
    CALL add_fields(ocean_state%p_acc%h     , ocean_state%p_prog(nnew(1))%h, cells)
    CALL add_fields(ocean_state%p_acc%u     , ocean_state%p_diag%u         , cells)
    CALL add_fields(ocean_state%p_acc%v     , ocean_state%p_diag%v         , cells)
    CALL add_fields(ocean_state%p_acc%rhopot, ocean_state%p_diag%rhopot    , cells)
    DO jtrc=1,no_tracer
      CALL add_fields(ocean_state%p_acc%tracer(:,:,:,jtrc),           &
        & ocean_state%p_prog(nnew(1))%tracer(:,:,:,jtrc), &
        & cells)
    END DO
    CALL add_fields(ocean_state%p_acc%u_vint        , ocean_state%p_diag%u_vint        , cells)
    CALL add_fields(ocean_state%p_acc%w             , ocean_state%p_diag%w             , cells,max_zlev+1)
    CALL add_fields(ocean_state%p_acc%div_mass_flx_c, ocean_state%p_diag%div_mass_flx_c, cells)
    CALL add_fields(ocean_state%p_acc%rho           , ocean_state%p_diag%rho           , cells)
    CALL add_fields(ocean_state%p_acc%vt            , ocean_state%p_diag%vt            , edges)
    CALL add_fields(ocean_state%p_acc%mass_flx_e    , ocean_state%p_diag%mass_flx_e    , edges)
    CALL add_fields(ocean_state%p_acc%vort          , ocean_state%p_diag%vort          , verts,max_zlev)
    CALL add_fields(ocean_state%p_acc%kin           , ocean_state%p_diag%kin           , cells)
    
    ! update forcing accumulated values
    CALL add_fields(p_sfc_flx%forc_wind_u_acc     , p_sfc_flx%forc_wind_u     , cells)
    CALL add_fields(p_sfc_flx%forc_wind_v_acc     , p_sfc_flx%forc_wind_v     , cells)
    CALL add_fields(p_sfc_flx%forc_swflx_acc      , p_sfc_flx%forc_swflx      , cells)
    CALL add_fields(p_sfc_flx%forc_lwflx_acc      , p_sfc_flx%forc_lwflx      , cells)
    CALL add_fields(p_sfc_flx%forc_ssflx_acc      , p_sfc_flx%forc_ssflx      , cells)
    CALL add_fields(p_sfc_flx%forc_slflx_acc      , p_sfc_flx%forc_slflx      , cells)
    CALL add_fields(p_sfc_flx%forc_precip_acc     , p_sfc_flx%forc_precip     , cells)
    CALL add_fields(p_sfc_flx%forc_evap_acc       , p_sfc_flx%forc_evap       , cells)
    CALL add_fields(p_sfc_flx%forc_runoff_acc     , p_sfc_flx%forc_runoff     , cells)
    CALL add_fields(p_sfc_flx%forc_fw_bc_acc      , p_sfc_flx%forc_fw_bc      , cells)
    CALL add_fields(p_sfc_flx%forc_fwrelax_acc    , p_sfc_flx%forc_fwrelax    , cells)
    CALL add_fields(p_sfc_flx%forc_fw_bc_oce_acc  , p_sfc_flx%forc_fw_bc_oce  , cells)
    CALL add_fields(p_sfc_flx%forc_fw_bc_ice_acc  , p_sfc_flx%forc_fw_bc_ice  , cells)
    CALL add_fields(p_sfc_flx%forc_fw_ice_vol_acc,  p_sfc_flx%forc_fw_ice_vol , cells)
    CALL add_fields(p_sfc_flx%forc_fw_tot_acc     , p_sfc_flx%forc_fw_tot     , cells)
    CALL add_fields(p_sfc_flx%forc_hfrelax_acc    , p_sfc_flx%forc_hfrelax    , cells)
    CALL add_fields(p_sfc_flx%forc_hflx_acc       , p_sfc_flx%forc_hflx       , cells)
    DO jtrc=1,no_tracer
      CALL add_fields(p_sfc_flx%forc_tracer_acc(:,:,jtrc), p_sfc_flx%forc_tracer(:,:,jtrc), cells)
      CALL add_fields(p_sfc_flx%forc_tracer_relax_acc(:,:,jtrc), p_sfc_flx%forc_tracer_relax(:,:,jtrc), cells)
    END DO
    
  END SUBROUTINE update_ocean_statistics
  
  SUBROUTINE compute_mean_ocean_statistics(p_acc,p_sfc_flx,nsteps_since_last_output)
    TYPE(t_hydro_ocean_acc), INTENT(inout) :: p_acc
    TYPE(t_sfc_flx),         INTENT(inout) :: p_sfc_flx
    INTEGER,INTENT(in)                     :: nsteps_since_last_output
    
    
    !TODO [ram] trigger aoutput wrt to multiple output intervals
    !TODO [ram] CALL collect_group(TRIM(grp_name), grp_vars_fg_sfc, ngrp_vars_fg_sfc,    &
    !TODO [ram]       &                loutputvars_only=.FALSE.,lremap_lonlat=.FALSE.)
    
    !TODO [ram]   IF (ALLOCATED(output_file)) THEN
    !TODO [ram]      DO i = 1, nvar_lists
    !TODO [ram]         element => var_lists(i)%p%first_list_element
    !TODO [ram]         DO
    !TODO [ram]           IF(.NOT. ASSOCIATED(element)) EXIT
    !TODO [ram]           if (.not. element%isteptype==TSTEP_ACCUM) cycle
    !TODO [ram]           DO i = 1, SIZE(output_file)
    !TODO [ram]              if (istime4output*output_file(i)%out_event) then
    !TODO [ram]                 if (one_of(element%name, output_file(i)%namelist%ml_varlist(1:output_file(i)%num_vars)) then
    !TODO [ram]               ! do accumulation for variable varlist(j)
    !TODO [ram]                   element%field%rptr = element%field%rptr /
    !TODO [ram]            end if
    !TODO [ram]         end if
    !TODO [ram]      end DO
    !TODO [ram]   end IF
    
    
    p_acc%tracer                    = p_acc%tracer                   /REAL(nsteps_since_last_output,wp)
    p_acc%h                         = p_acc%h                        /REAL(nsteps_since_last_output,wp)
    p_acc%u                         = p_acc%u                        /REAL(nsteps_since_last_output,wp)
    p_acc%v                         = p_acc%v                        /REAL(nsteps_since_last_output,wp)
    p_acc%rhopot                    = p_acc%rhopot                   /REAL(nsteps_since_last_output,wp)
    p_acc%u_vint                    = p_acc%u_vint                   /REAL(nsteps_since_last_output,wp)
    p_acc%w                         = p_acc%w                        /REAL(nsteps_since_last_output,wp)
    p_acc%div_mass_flx_c            = p_acc%div_mass_flx_c           /REAL(nsteps_since_last_output,wp)
    p_acc%rho                       = p_acc%rho                      /REAL(nsteps_since_last_output,wp)
    p_acc%vt                        = p_acc%vt                       /REAL(nsteps_since_last_output,wp)
    p_acc%mass_flx_e                = p_acc%mass_flx_e               /REAL(nsteps_since_last_output,wp)
    p_acc%vort                      = p_acc%vort                     /REAL(nsteps_since_last_output,wp)
    p_acc%kin                       = p_acc%kin                      /REAL(nsteps_since_last_output,wp)
    p_sfc_flx%forc_wind_u_acc       = p_sfc_flx%forc_wind_u_acc      /REAL(nsteps_since_last_output,wp)
    p_sfc_flx%forc_wind_v_acc       = p_sfc_flx%forc_wind_v_acc      /REAL(nsteps_since_last_output,wp)
    p_sfc_flx%forc_swflx_acc        = p_sfc_flx%forc_swflx_acc       /REAL(nsteps_since_last_output,wp)
    p_sfc_flx%forc_lwflx_acc        = p_sfc_flx%forc_lwflx_acc       /REAL(nsteps_since_last_output,wp)
    p_sfc_flx%forc_ssflx_acc        = p_sfc_flx%forc_ssflx_acc       /REAL(nsteps_since_last_output,wp)
    p_sfc_flx%forc_slflx_acc        = p_sfc_flx%forc_slflx_acc       /REAL(nsteps_since_last_output,wp)
    p_sfc_flx%forc_precip_acc       = p_sfc_flx%forc_precip_acc      /REAL(nsteps_since_last_output,wp)
    p_sfc_flx%forc_evap_acc         = p_sfc_flx%forc_evap_acc        /REAL(nsteps_since_last_output,wp)
    p_sfc_flx%forc_runoff_acc       = p_sfc_flx%forc_runoff_acc      /REAL(nsteps_since_last_output,wp)
    p_sfc_flx%forc_fw_bc_acc        = p_sfc_flx%forc_fw_bc_acc       /REAL(nsteps_since_last_output,wp)
    p_sfc_flx%forc_fw_bc_oce_acc    = p_sfc_flx%forc_fw_bc_oce_acc   /REAL(nsteps_since_last_output,wp)
    p_sfc_flx%forc_fw_bc_ice_acc    = p_sfc_flx%forc_fw_bc_ice_acc   /REAL(nsteps_since_last_output,wp)
    p_sfc_flx%forc_fwrelax_acc      = p_sfc_flx%forc_fwrelax_acc     /REAL(nsteps_since_last_output,wp)
    p_sfc_flx%forc_fw_ice_vol_acc   = p_sfc_flx%forc_fw_ice_vol_acc  /REAL(nsteps_since_last_output,wp)
    p_sfc_flx%forc_fw_tot_acc       = p_sfc_flx%forc_fw_tot_acc      /REAL(nsteps_since_last_output,wp)
    p_sfc_flx%forc_hfrelax_acc      = p_sfc_flx%forc_hfrelax_acc     /REAL(nsteps_since_last_output,wp)
    p_sfc_flx%forc_hflx_acc         = p_sfc_flx%forc_hflx_acc        /REAL(nsteps_since_last_output,wp)
    IF(no_tracer>0)THEN
      p_sfc_flx%forc_tracer_acc       = p_sfc_flx%forc_tracer_acc      /REAL(nsteps_since_last_output,wp)
      p_sfc_flx%forc_tracer_relax_acc = p_sfc_flx%forc_tracer_relax_acc/REAL(nsteps_since_last_output,wp)
    ENDIF  
  END SUBROUTINE compute_mean_ocean_statistics
  
  SUBROUTINE reset_ocean_statistics(p_acc,p_sfc_flx,nsteps_since_last_output)
    TYPE(t_hydro_ocean_acc), INTENT(inout) :: p_acc
    TYPE(t_sfc_flx),         INTENT(inout) :: p_sfc_flx
    INTEGER,OPTIONAL,        INTENT(inout) :: nsteps_since_last_output
    
    IF (PRESENT(nsteps_since_last_output)) nsteps_since_last_output        = 0
    
    p_acc%tracer                    = 0.0_wp
    p_acc%h                         = 0.0_wp
    p_acc%u                         = 0.0_wp
    p_acc%v                         = 0.0_wp
    p_acc%rhopot                    = 0.0_wp
    p_acc%u_vint                    = 0.0_wp
    p_acc%w                         = 0.0_wp
    p_acc%div_mass_flx_c            = 0.0_wp
    p_acc%rho                       = 0.0_wp
    p_acc%vt                        = 0.0_wp
    p_acc%mass_flx_e                = 0.0_wp
    p_acc%vort                      = 0.0_wp
    p_acc%kin                       = 0.0_wp
    p_sfc_flx%forc_wind_u_acc       = 0.0_wp
    p_sfc_flx%forc_wind_v_acc       = 0.0_wp
    p_sfc_flx%forc_swflx_acc        = 0.0_wp
    p_sfc_flx%forc_lwflx_acc        = 0.0_wp
    p_sfc_flx%forc_ssflx_acc        = 0.0_wp
    p_sfc_flx%forc_slflx_acc        = 0.0_wp
    p_sfc_flx%forc_precip_acc       = 0.0_wp
    p_sfc_flx%forc_evap_acc         = 0.0_wp
    p_sfc_flx%forc_runoff_acc       = 0.0_wp
    p_sfc_flx%forc_fw_bc_acc        = 0.0_wp
    p_sfc_flx%forc_fw_bc_oce_acc    = 0.0_wp
    p_sfc_flx%forc_fw_bc_ice_acc    = 0.0_wp
    p_sfc_flx%forc_fwrelax_acc      = 0.0_wp
    p_sfc_flx%forc_fw_ice_vol_acc   = 0.0_wp
    p_sfc_flx%forc_fw_tot_acc       = 0.0_wp
    p_sfc_flx%forc_hfrelax_acc      = 0.0_wp
    p_sfc_flx%forc_hflx_acc         = 0.0_wp
    IF(no_tracer>0)THEN
      p_sfc_flx%forc_tracer_acc       = 0.0_wp
      p_sfc_flx%forc_tracer_relax_acc = 0.0_wp
    ENDIF
  END SUBROUTINE reset_ocean_statistics
  
  SUBROUTINE new_ocean_statistics()
  END SUBROUTINE new_ocean_statistics
  
  SUBROUTINE set_output_pointers(timelevel,p_diag,p_prog)
    INTEGER, INTENT(in) :: timelevel
    TYPE(t_hydro_ocean_diag) :: p_diag
    TYPE(t_hydro_ocean_prog) :: p_prog
    
    TYPE(t_list_element), POINTER :: output_var => NULL()
    TYPE(t_list_element), POINTER :: prog_var   => NULL()
    CHARACTER(LEN=max_char_length) :: timelevel_str
    !-------------------------------------------------------------------------
    WRITE(timelevel_str,'(a,i2.2)') '_TL',timelevel
    !write(0,*)'>>>>>>>>>>>>>>>> T timelevel_str:',TRIM(timelevel_str)
    
    !CALL print_var_list(ocean_restart_list)
    !prog_var               => find_list_element(ocean_restart_list,'h'//TRIM(timelevel_str))
    !output_var             => find_list_element(ocean_restart_list,'h')
    !output_var%field%r_ptr => prog_var%field%r_ptr
    !p_diag%h               => prog_var%field%r_ptr(:,:,1,1,1)
    p_diag%h               =  p_prog%h
    
    !output_var             => find_list_element(ocean_restart_list,'vn')
    !prog_var               => find_list_element(ocean_restart_list,'vn'//TRIM(timelevel_str))
    !output_var%field%r_ptr => prog_var%field%r_ptr
    !p_diag%vn              => prog_var%field%r_ptr(:,:,:,1,1)
    p_diag%vn(:,:,:)       =  p_prog%vn
    
    !output_var             => find_list_element(ocean_restart_list,'t')
    !prog_var               => find_list_element(ocean_restart_list,'t'//TRIM(timelevel_str))
    !output_var%field%r_ptr => prog_var%field%r_ptr
    !p_diag%t               => prog_var%field%r_ptr(:,:,:,1,1)
    IF(no_tracer>0)p_diag%t(:,:,:)        =  p_prog%tracer(:,:,:,1)
    
    !output_var             => find_list_element(ocean_restart_list,'s')
    !prog_var               => find_list_element(ocean_restart_list,'s'//TRIM(timelevel_str))
    !output_var%field%r_ptr => prog_var%field%r_ptr
    !p_diag%s               => prog_var%field%r_ptr(:,:,:,1,1)
     IF(no_tracer>1)p_diag%s(:,:,:)        =  p_prog%tracer(:,:,:,2)
  END SUBROUTINE set_output_pointers
  
END MODULE mo_hydro_ocean_run
