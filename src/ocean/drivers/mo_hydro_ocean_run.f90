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
  USE mo_sync,                   ONLY: sync_patch_array, sync_e, sync_c !, sync_v
  USE mo_ocean_nml,              ONLY: iswm_oce, n_zlev, no_tracer, &
    & eos_type, i_sea_ice, l_staggered_timestep, gibraltar
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
  USE mo_oce_state,              ONLY: ocean_restart_list
 ! USE mo_ocean_initialization,   ONLY: set_lateral_boundary_values
  USE mo_oce_math_operators,     ONLY: calculate_thickness
  USE mo_operator_ocean_coeff_3d,ONLY: t_operator_coeff! , update_diffusion_matrices
  USE mo_scalar_product,         ONLY: calc_scalar_product_veloc_3d
  USE mo_oce_tracer,             ONLY: advect_tracer_ab
  USE mo_io_restart,             ONLY: write_restart_info_file, create_restart_file
  USE mo_oce_bulk,               ONLY: update_sfcflx
  USE mo_sea_ice,                ONLY: update_ice_statistic, compute_mean_ice_statistics, reset_ice_statistics
  USE mo_sea_ice_types,          ONLY: t_sfc_flx, t_atmos_fluxes, t_atmos_for_ocean, &
    & t_sea_ice
  USE mo_oce_physics,            ONLY: t_ho_params, update_ho_params
  USE mo_oce_thermodyn,          ONLY: calc_density_mpiom_func, calc_density_lin_eos_func,&
    & calc_density_jmdwfg06_eos_func, calc_potential_density, &
    & calc_density
  USE mo_name_list_output,       ONLY: write_name_list_output, istime4name_list_output
  USE mo_oce_diagnostics,        ONLY: calc_fast_oce_diagnostics, &
    & t_oce_timeseries, calc_moc, calc_psi
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
  USE mo_ocean_statistics
  USE mo_ocean_output
  
  IMPLICIT NONE
  
  PRIVATE
  
  !VERSION CONTROL:
  CHARACTER(LEN=*), PARAMETER :: version = '$Id$'
  
  ! public interface
  !
  ! public subroutines
  PUBLIC :: perform_ho_stepping
  PUBLIC :: prepare_ho_stepping
  PRIVATE:: update_intermediate_tracer_vars
  !
  CHARACTER(LEN=12)  :: str_module = 'HYDRO-ocerun'  ! Output of module for 1 line debug
  INTEGER :: idt_src                      ! Level of detail for 1 line debug
  !
  !-------------------------------------------------------------------------
  
CONTAINS

  !-------------------------------------------------------------------------
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

!    CALL update_diffusion_matrices( patch_3d,         &
!      & p_phys_param,                 &
!      & operators_coefficients%matrix_vert_diff_e,&
!      & operators_coefficients%matrix_vert_diff_c)

     CALL init_ho_lhs_fields_mimetic   ( patch_3d )

  END SUBROUTINE prepare_ho_stepping
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  !>
  !! Main stepping routine for call of hydrostatic ocean model
  !!
  !! @par Revision History
  !! Developed by Peter Korn, MPI-M  (2008-2010).
  !! Initial release by Stephan Lorenz, MPI-M (2010-07)
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
    !LOGICAL                         :: l_outputtime
    CHARACTER(LEN=32)               :: datestring, plaindatestring
    TYPE(t_patch), POINTER :: patch_2d
    TYPE(t_patch_vert), POINTER :: patch_1d
    INTEGER, POINTER :: dolic(:,:)
    REAL(wp), POINTER :: prism_thickness(:,:,:)
    INTEGER :: jstep0 ! start counter for time loop
    
    !CHARACTER(LEN=filename_max)  :: outputfile, gridfile
    CHARACTER(LEN=max_char_length), PARAMETER :: &
      & routine = 'mo_hydro_ocean_run:perform_ho_stepping'
    !------------------------------------------------------------------
    
    patch_2D      => patch_3d%p_patch_2d(1)
    
    !------------------------------------------------------------------
    ! no grid refinement allowed here so far
    !------------------------------------------------------------------
    IF (n_dom > 1 ) THEN
      CALL finish(TRIM(routine), ' N_DOM > 1 is not allowed')
    END IF
    jg = n_dom
    
    patch_2d => patch_3d%p_patch_2d(jg)
    
    CALL datetime_to_string(datestring, datetime)

    ! IF (ltimer) CALL timer_start(timer_total)
    CALL timer_start(timer_total)
    
    time_config%sim_time(:) = 0.0_wp
    
    !------------------------------------------------------------------
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
          & patch_2D%cells%owned,   &
          & patch_2D%edges%owned,   &
          & patch_2D%verts%owned,   &
          & n_zlev)
        IF (i_sea_ice >= 1) CALL update_ice_statistic(p_ice%acc, p_ice,patch_2D%cells%owned)

        CALL write_name_list_output(jstep=jstep0)

        CALL reset_ocean_statistics(ocean_state(1)%p_acc,p_sfc_flx)
        IF (i_sea_ice >= 1) CALL reset_ice_statistics(p_ice%acc)
      ENDIF

    ENDIF ! output_mode%l_nml
    !------------------------------------------------------------------
    ! call the dynamical core: start the time loop
    !------------------------------------------------------------------
    
      time_loop: DO jstep = (jstep0+1), (jstep0+nsteps)
      
        CALL datetime_to_string(datestring, datetime)
        WRITE(message_text,'(a,i10,2a)') '  Begin of timestep =',jstep,'  datetime:  ', datestring
        CALL message (TRIM(routine), message_text)
      
          !In case of a time-varying forcing:
          IF (ltimer) CALL timer_start(timer_upd_flx)
          CALL update_sfcflx( patch_3d, ocean_state(jg), p_as, p_ice, p_atm_f, p_sfc_flx, &
          & jstep, datetime, operators_coefficients)

          ! update_sfcflx has changed p_prog(nold(1))%h
          CALL sync_patch_array(sync_c, patch_2D, ocean_state(jg)%p_prog(nold(1))%h)
        
          IF(.NOT.l_staggered_timestep)THEN
          
            CALL calculate_thickness( patch_3d, ocean_state(jg), p_ext_data(jg), operators_coefficients)
          
          !  CALL set_lateral_boundary_values( patch_3d, ocean_state(jg)%p_prog(nold(1))%vn)
          !  CALL sync_patch_array(sync_e, patch_3d%p_patch_2d(jg), ocean_state(jg)%p_prog(nold(1))%vn)
          
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
        
!          CALL update_diffusion_matrices( patch_3d,   &
!          & p_phys_param,                             &
!          & operators_coefficients%matrix_vert_diff_e,&
!          & operators_coefficients%matrix_vert_diff_c)
        
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
            &                          ocean_state(jg)%p_prog(nold(1))%tracer,                         &
            &                          ocean_state(jg)%p_diag%rhopot )
          CALL calc_psi (patch_2d,patch_3d, ocean_state(jg)%p_diag%u(:,:,:),&
            &                          ocean_state(jg)%p_prog(nold(1))%h(:,:),                         &
            &                          ocean_state(jg)%p_diag%u_vint, datetime)
        ENDIF 
        ! update accumulated vars
        CALL update_ocean_statistics(ocean_state(1),&
          &                          p_sfc_flx,           &
          &                          patch_2D%cells%owned,&
          &                          patch_2D%edges%owned,&
          &                          patch_2D%verts%owned,&
        & n_zlev)
        IF (i_sea_ice >= 1) CALL update_ice_statistic(p_ice%acc,p_ice,patch_2D%cells%owned)
      
        dolic           => patch_3d%p_patch_1d(1)%dolic_c
        prism_thickness => patch_3d%p_patch_1d(1)%prism_thick_c
        CALL calc_fast_oce_diagnostics( patch_2D,      &
          &                             dolic, &
          &                             prism_thickness, &
          &                             patch_3d%p_patch_1d(1)%zlev_m, &
          &                             ocean_state(jg)%p_diag)

        CALL output_ocean( patch_3d, ocean_state, p_ext_data,          &
          &                datetime, lwrite_restart,            &
          &                p_sfc_flx, p_phys_param,             &
          &                p_as, p_atm_f, p_ice,operators_coefficients, &
          &                jstep, jstep0)
      
        ! Shift time indices for the next loop
        ! this HAS to ge into the restart files, because the start with the following loop
        CALL update_time_indices(jg)
        ! update intermediate timestepping variables for the tracers
        CALL update_intermediate_tracer_vars(ocean_state(jg))
      
        ! write a restart or checkpoint file
        IF (MOD(jstep,n_checkpoints())==0 .OR. ((jstep==(jstep0+nsteps)) .AND. lwrite_restart)) THEN
          CALL create_restart_file( patch = patch_2d,       &
            &                       datetime=datetime,      &
            &                       jstep=jstep,            &
            &                       model_type="oce",       &
            &                       opt_depth=n_zlev,       &
            &                       opt_sim_time=time_config%sim_time(1),&
            &                       opt_nice_class=1)
          ! Create the master (meta) file in ASCII format which contains
          ! info about which files should be read in for a restart run.
          CALL write_restart_info_file
        END IF
      
      ENDDO time_loop
      
    CALL timer_stop(timer_total)
    
  END SUBROUTINE perform_ho_stepping
  !-------------------------------------------------------------------------
  

  !-------------------------------------------------------------------------
  SUBROUTINE update_intermediate_tracer_vars(ocean_state)
    TYPE(t_hydro_ocean_state), INTENT(inout) :: ocean_state
    
    ! velocity
    ocean_state%p_aux%g_nm1 = ocean_state%p_aux%g_n
    ocean_state%p_aux%g_n   = 0.0_wp
  END SUBROUTINE update_intermediate_tracer_vars
  !-------------------------------------------------------------------------
  
  
    
END MODULE mo_hydro_ocean_run
