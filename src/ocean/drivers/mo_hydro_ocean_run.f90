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
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
!!
MODULE mo_hydro_ocean_run
  !-------------------------------------------------------------------------
  USE mo_kind,                   ONLY: wp
  USE mo_impl_constants,         ONLY: max_char_length
  USE mo_model_domain,           ONLY: t_patch, t_patch_3d
  USE mo_grid_config,            ONLY: n_dom
  USE mo_ocean_nml,              ONLY: iswm_oce, n_zlev, no_tracer, &
    & i_sea_ice, cfl_check, cfl_threshold, cfl_stop_on_violation,   &
    & cfl_write, surface_module
  USE mo_ocean_nml,              ONLY: iforc_oce, Coupled_FluxFromAtmo
  USE mo_dynamics_config,        ONLY: nold, nnew
  USE mo_io_config,              ONLY: n_checkpoints, write_last_restart
  USE mo_run_config,             ONLY: nsteps, dtime, ltimer, output_mode, debug_check_level
  USE mo_exception,              ONLY: message, message_text, finish
  USE mo_ext_data_types,         ONLY: t_external_data
  !USE mo_io_units,               ONLY: filename_max
  USE mo_datetime,               ONLY: t_datetime, add_time, datetime_to_string
  USE mo_timer,                  ONLY: timer_start, timer_stop, timer_total, timer_solve_ab,  &
    & timer_tracer_ab, timer_vert_veloc, timer_normal_veloc,     &
    & timer_upd_phys, timer_upd_flx, timer_extra20, timers_level, &
    & timer_scalar_prod_veloc, timer_extra21, timer_extra22
  USE mo_ocean_ab_timestepping,    ONLY: solve_free_surface_eq_ab, &
    & calc_normal_velocity_ab,  &
    & calc_vert_velocity,       &
    & update_time_indices
  USE mo_ocean_types,              ONLY: t_hydro_ocean_state, &
    & t_operator_coeff, t_solvercoeff_singleprecision
  USE mo_ocean_math_operators,   ONLY: update_height_depdendent_variables, check_cfl_horizontal, check_cfl_vertical
  USE mo_scalar_product,         ONLY: calc_scalar_product_veloc_3d
  USE mo_ocean_tracer,             ONLY: advect_tracer_ab
  USE mo_io_restart,             ONLY: create_restart_file
  USE mo_ocean_bulk,             ONLY: update_surface_flux
  USE mo_ocean_surface,          ONLY: update_ocean_surface
  USE mo_ocean_surface_types,    ONLY: t_ocean_surface
  USE mo_sea_ice,                ONLY: update_ice_statistic, reset_ice_statistics
  USE mo_sea_ice_types,          ONLY: t_sfc_flx, t_atmos_fluxes, t_atmos_for_ocean, &
    & t_sea_ice
  USE mo_ocean_physics,            ONLY: t_ho_params, update_ho_params
  USE mo_ocean_thermodyn,          ONLY: calc_potential_density, &
    & calculate_density! , ocean_correct_ThermoExpansion
  USE mo_name_list_output,       ONLY: write_name_list_output
  USE mo_ocean_diagnostics,        ONLY: calc_fast_oce_diagnostics, calc_psi, calc_psi_vn
  USE mo_ocean_ab_timestepping_mimetic, ONLY: construct_ho_lhs_fields_mimetic, destruct_ho_lhs_fields_mimetic
  USE mo_io_restart_attributes,  ONLY: get_restart_attribute
  USE mo_time_config,            ONLY: time_config
  USE mo_master_config,          ONLY: isRestart
!  USE mo_sea_ice_nml,            ONLY: i_ice_dyn
  USE mo_util_dbg_prnt,          ONLY: dbg_print, debug_printValue
  USE mo_dbg_nml,                ONLY: idbg_mxmn
  USE mo_statistics
  USE mo_ocean_statistics
  USE mo_derived_variable_handling, ONLY: perform_accumulation
  USE mo_ocean_output
  USE mo_ocean_coupling,         ONLY: couple_ocean_toatmo_fluxes

  IMPLICIT NONE

  PRIVATE

  PUBLIC  :: perform_ho_stepping
  PUBLIC  :: prepare_ho_stepping, end_ho_stepping
  PUBLIC  :: write_initial_ocean_timestep
  
  CHARACTER(LEN=12)  :: str_module = 'HYDRO-ocerun'  ! Output of module for 1 line debug
  INTEGER            :: idt_src    = 1               ! Level of detail for 1 line debug
  !-------------------------------------------------------------------------

CONTAINS

  !-------------------------------------------------------------------------
  !<Optimize:inUse>
  SUBROUTINE prepare_ho_stepping(patch_3d, operators_coefficients, ocean_state, ext_data, is_restart, &
    & solvercoeff_sp)
    TYPE(t_patch_3d ), INTENT(in)     :: patch_3d
    TYPE(t_operator_coeff)            :: operators_coefficients
    TYPE(t_hydro_ocean_state), TARGET :: ocean_state
    TYPE(t_external_data), TARGET, INTENT(in) :: ext_data
! !   TYPE (t_ho_params)                :: p_phys_param
    LOGICAL, INTENT(in)               :: is_restart
    TYPE(t_solvercoeff_singleprecision), INTENT(inout) :: solvercoeff_sp
! 
!     IF (is_restart) THEN
!       ! Prepare ocean_state%p_prog, since it is needed by the sea ice model (e.g. wind stress computation)
!       IF ( i_sea_ice > 0 )         &
!       CALL update_height_depdendent_variables( patch_3d, ocean_state, ext_data, operators_coefficients, solvercoeff_sp)
!       
!       CALL calc_scalar_product_veloc_3d( patch_3d,  &
!         & ocean_state%p_prog(nold(1))%vn,         &
!         & ocean_state%p_diag,                     &
!         & operators_coefficients)
!     ELSE
!     ENDIF
! 
!     !    CALL update_diffusion_matrices( patch_3d,         &
!     !      & p_phys_param,                 &
!     !      & operators_coefficients%matrix_vert_diff_e,&
!     !      & operators_coefficients%matrix_vert_diff_c)
! 
     CALL construct_ho_lhs_fields_mimetic   ( patch_3d )
! 
  END SUBROUTINE prepare_ho_stepping
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !<Optimize:inUse>
  SUBROUTINE end_ho_stepping()

     CALL destruct_ho_lhs_fields_mimetic()
    
  END SUBROUTINE end_ho_stepping
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !>
  !! Main stepping routine for call of hydrostatic ocean model
  !!
  !! @par Revision History
  !! Developed by Peter Korn, MPI-M  (2008-2010).
  !! Initial release by Stephan Lorenz, MPI-M (2010-07)
  !
!<Optimize:inUse>
  SUBROUTINE perform_ho_stepping( patch_3d, ocean_state, p_ext_data,          &
    & datetime,                                    &
    & surface_fluxes, p_sfc, p_phys_param,              &
    & p_as, p_atm_f, sea_ice,operators_coefficients, &
    & solvercoeff_sp)

    TYPE(t_patch_3d ),TARGET, INTENT(inout)          :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET, INTENT(inout) :: ocean_state(n_dom)
    TYPE(t_external_data), TARGET, INTENT(in)        :: p_ext_data(n_dom)
    TYPE(t_datetime), INTENT(inout)                  :: datetime
    TYPE(t_sfc_flx)                                  :: surface_fluxes
    TYPE(t_ocean_surface)                            :: p_sfc
    TYPE (t_ho_params)                               :: p_phys_param
    TYPE(t_atmos_for_ocean),  INTENT(inout)          :: p_as
    TYPE(t_atmos_fluxes ),    INTENT(inout)          :: p_atm_f
    TYPE (t_sea_ice),         INTENT(inout)          :: sea_ice
    TYPE(t_operator_coeff),   INTENT(inout)          :: operators_coefficients
    TYPE(t_solvercoeff_singleprecision), INTENT(inout) :: solvercoeff_sp

    ! local variables
    INTEGER :: jstep, jg, return_status
    !LOGICAL                         :: l_outputtime
    CHARACTER(LEN=32)               :: datestring
    TYPE(t_patch), POINTER :: patch_2d
    INTEGER :: jstep0 ! start counter for time loop
    REAL(wp) :: mean_height, old_mean_height
    REAL(wp) :: verticalMeanFlux(n_zlev+1)
    INTEGER :: level
    !CHARACTER(LEN=filename_max)  :: outputfile, gridfile
    CHARACTER(LEN=max_char_length), PARAMETER :: &
      & routine = 'mo_hydro_ocean_run:perform_ho_stepping'
    !------------------------------------------------------------------

    patch_2d      => patch_3d%p_patch_2d(1)

    !------------------------------------------------------------------
    ! no grid refinement allowed here so far
    !------------------------------------------------------------------
    IF (n_dom > 1 ) THEN
      CALL finish(TRIM(routine), ' N_DOM > 1 is not allowed')
    END IF
    jg = n_dom

    patch_2d => patch_3d%p_patch_2d(jg)

    ! CALL datetime_to_string(datestring, datetime)

    time_config%sim_time(:) = 0.0_wp

    !------------------------------------------------------------------
    jstep0 = 0
    IF (isRestart() .AND. .NOT. time_config%is_relative_time) THEN
      ! get start counter for time loop from restart file:
      CALL get_restart_attribute("jstep", jstep0)
    END IF
    IF (isRestart() .AND. mod(nold(jg),2) /=1 ) THEN
      ! swap the g_n and g_nm1
      CALL update_time_g_n(ocean_state(jg))
    ENDIF
    !------------------------------------------------------------------
    ! call the dynamical core: start the time loop
    !------------------------------------------------------------------
    ! IF (ltimer) CALL timer_start(timer_total)
    CALL timer_start(timer_total)

    time_loop: DO jstep = (jstep0+1), (jstep0+nsteps)
      ! write(0,*) "nold nnew=", nold(1), nnew(1)

      CALL datetime_to_string(datestring, datetime)
      WRITE(message_text,'(a,i10,2a)') '  Begin of timestep =',jstep,'  datetime:  ', datestring
      CALL message (TRIM(routine), message_text)
      
      ! Set model time - now before doing the calculation:
      !  - datetime refers to the currently calculating timestep
      CALL add_time(dtime,0,0,0,datetime)
      ! Not nice, but the name list output requires this - needed?
      time_config%sim_time(1) = time_config%sim_time(1) + dtime
      
      IF (timers_level > 2)  CALL timer_start(timer_extra22)
      CALL update_height_depdendent_variables( patch_3d, ocean_state(jg), p_ext_data(jg), operators_coefficients, solvercoeff_sp)
      IF (timers_level > 2)  CALL timer_stop(timer_extra22)
      
      IF (timers_level > 2) CALL timer_start(timer_scalar_prod_veloc)
      CALL calc_scalar_product_veloc_3d( patch_3d,  &
        & ocean_state(jg)%p_prog(nold(1))%vn,         &
        & ocean_state(jg)%p_diag,                     &
        & operators_coefficients)
      IF (timers_level > 2) CALL timer_stop(timer_scalar_prod_veloc)
      
      !In case of a time-varying forcing:
      ! update_surface_flux or update_ocean_surface has changed p_prog(nold(1))%h, SST and SSS
      IF (ltimer) CALL timer_start(timer_upd_flx)
      IF (surface_module == 1) THEN
        CALL update_surface_flux( patch_3d, ocean_state(jg), p_as, sea_ice, p_atm_f, surface_fluxes, &
          & jstep, datetime, operators_coefficients)
      ELSEIF (surface_module == 2) THEN
        CALL update_ocean_surface( patch_3d, ocean_state(jg), p_as, sea_ice, p_atm_f, surface_fluxes, p_sfc, &
          & jstep, datetime, operators_coefficients)
      ENDIF
      IF (ltimer) CALL timer_stop(timer_upd_flx)

      IF (timers_level > 2)  CALL timer_start(timer_extra22)
      CALL update_height_depdendent_variables( patch_3d, ocean_state(jg), p_ext_data(jg), operators_coefficients, solvercoeff_sp)
      IF (timers_level > 2)  CALL timer_stop(timer_extra22)

!       IF (timers_level > 2) CALL timer_start(timer_scalar_prod_veloc)
!       CALL calc_scalar_product_veloc_3d( patch_3d,  &
!         & ocean_state(jg)%p_prog(nold(1))%vn,         &
!         & ocean_state(jg)%p_diag,                     &
!         & operators_coefficients)
!       IF (timers_level > 2) CALL timer_stop(timer_scalar_prod_veloc)

      !---------DEBUG DIAGNOSTICS-------------------------------------------
      idt_src=3  ! output print level (1-5, fix)
      CALL dbg_print('on entry: h-old'           ,ocean_state(jg)%p_prog(nold(1))%h ,str_module,idt_src, &
        & patch_2d%cells%owned )
      CALL dbg_print('on entry: h-new'           ,ocean_state(jg)%p_prog(nnew(1))%h ,str_module,idt_src, &
        & patch_2d%cells%owned )
      CALL dbg_print('HydOce: ScaProdVel kin'    ,ocean_state(jg)%p_diag%kin        ,str_module,idt_src, &
        & patch_2d%cells%owned )
      CALL dbg_print('HydOce: ScaProdVel ptp_vn' ,ocean_state(jg)%p_diag%ptp_vn     ,str_module,idt_src, &
        & patch_2d%edges%owned )
      !---------------------------------------------------------------------

      CALL update_ho_params(patch_3d, ocean_state(jg), p_as%fu10, sea_ice%concsum, p_phys_param, operators_coefficients)

      !------------------------------------------------------------------------
      IF (debug_check_level > 5) THEN
        CALL horizontal_mean(values=ocean_state(jg)%p_prog(nold(1))%h(:,:), weights=patch_2d%cells%area(:,:), &
          & in_subset=patch_2d%cells%owned, mean=old_mean_height)
      END IF
      !------------------------------------------------------------------------
      ! solve for new free surface
      IF (ltimer) CALL timer_start(timer_solve_ab)
      CALL solve_free_surface_eq_ab (patch_3d, ocean_state(jg), p_ext_data(jg), &
        & surface_fluxes, p_phys_param, jstep, operators_coefficients, solvercoeff_sp, return_status)!, p_int(jg))
      IF (return_status /= 0) THEN
       CALL output_ocean(              &
         & patch_3d=patch_3d,          &
         & ocean_state=ocean_state,    &
         & datetime=datetime,          &
         & surface_fluxes=surface_fluxes, &
         & sea_ice=sea_ice,            &
         & jstep=jstep, jstep0=jstep0, &
         & force_output=.true.)
        CALL finish(TRIM(routine), 'solve_free_surface_eq_ab  returned error')
      ENDIF
      
      IF (ltimer) CALL timer_stop(timer_solve_ab)

      !------------------------------------------------------------------------
      ! Step 4: calculate final normal velocity from predicted horizontal
      ! velocity vn_pred and updated surface height
      IF (ltimer) CALL timer_start(timer_normal_veloc)
      CALL calc_normal_velocity_ab(patch_3d, ocean_state(jg),&
        & operators_coefficients, solvercoeff_sp,  p_ext_data(jg), p_phys_param)
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

      IF (idbg_mxmn >= 2 .OR. debug_check_level > 5) THEN
        CALL horizontal_mean(values=ocean_state(jg)%p_prog(nnew(1))%h(:,:), weights=patch_2d%cells%area(:,:), &
          & in_subset=patch_2d%cells%owned, mean=mean_height)
        CALL debug_printValue(description="Mean Height", value=mean_height, detail_level=2)
      ENDIF
      IF (debug_check_level > 5 .AND. idbg_mxmn >= 2) THEN
        ! check difference from old_mean_height
        CALL debug_printValue(description="Old/New Mean Height", value=old_mean_height, &
          & value1=mean_height, value2=(mean_height-old_mean_height) / old_mean_height, &
          & detail_level=2)
        ! check if vertical and horizontal fluxes add to 0
!         ocean_state(jg)%p_diag%w
        CALL horizontal_mean(values=ocean_state(jg)%p_diag%w, weights=patch_2d%cells%area(:,:), &
          & in_subset=patch_2d%cells%owned, mean=verticalMeanFlux, start_level=2, end_level=n_zlev-1)
        
        DO level=2, n_zlev-1
          CALL debug_printValue(description="Mean vertical flux at", value=REAL(level,wp),  &
            & value1=verticalMeanFlux(level), detail_level=2)
        ENDDO         
      END IF
      !------------------------------------------------------------------------

      !------------------------------------------------------------------------
      ! Step 6: transport tracers and diffuse them
      IF (no_tracer>=1) THEN
        IF (ltimer) CALL timer_start(timer_tracer_ab)
        CALL advect_tracer_ab( patch_3d, ocean_state(jg), p_phys_param,&
          & surface_fluxes,&
          & operators_coefficients,&
          & jstep)
        IF (ltimer) CALL timer_stop(timer_tracer_ab)
      ENDIF

!       ! One integration cycle finished. Set model time.
!       CALL add_time(dtime,0,0,0,datetime)
! 
!       ! Not nice, but the name list output requires this
!       time_config%sim_time(1) = time_config%sim_time(1) + dtime

      ! perform accumulation for special variables
      IF (timers_level > 2)  CALL timer_start(timer_extra20)
      
      IF (no_tracer>=1) THEN
        CALL calc_potential_density( patch_3d,                            &
          & ocean_state(jg)%p_prog(nold(1))%tracer,                       &
          & ocean_state(jg)%p_diag%rhopot )
          
        ! calculate diagnostic barotropic stream function
        CALL calc_psi (patch_3d, ocean_state(jg)%p_diag%u(:,:,:),         &
          & patch_3D%p_patch_1d(1)%prism_thick_c(:,:,:),                  &
          & ocean_state(jg)%p_diag%u_vint, datetime)
        CALL dbg_print('calc_psi: u_vint' ,ocean_state(jg)%p_diag%u_vint, str_module, 3, in_subset=patch_2d%cells%owned)
          
        ! calculate diagnostic barotropic stream function with vn
    !  not yet mature
    !   CALL calc_psi_vn (patch_3d, ocean_state(jg)%p_prog(nold(1))%vn,   &
    !     & patch_3D%p_patch_1d(1)%prism_thick_e(:,:,:),                  &
    !     & operators_coefficients,                                       &
    !     & ocean_state(jg)%p_diag%u_vint, ocean_state(jg)%p_diag%v_vint, datetime)
    !   CALL dbg_print('calc_psi_vn: u_vint' ,ocean_state(jg)%p_diag%u_vint, str_module, 5, in_subset=patch_2d%cells%owned)
    !   CALL dbg_print('calc_psi_vn: v_vint' ,ocean_state(jg)%p_diag%v_vint, str_module, 5, in_subset=patch_2d%cells%owned)
      ENDIF

      ! update accumulated vars
      CALL update_ocean_statistics(ocean_state(1),&
        & surface_fluxes, &
        & patch_2d%cells%owned,&
        & patch_2d%edges%owned,&
        & patch_2d%verts%owned,&
        & n_zlev,p_phys_param=p_phys_param)
        
      IF (i_sea_ice >= 1) CALL update_ice_statistic(sea_ice%acc,sea_ice,patch_2d%cells%owned)

      CALL calc_fast_oce_diagnostics( patch_2d,      &
        & patch_3d%p_patch_1d(1)%dolic_c, &
        & patch_3d%p_patch_1d(1)%prism_thick_c, &
        & patch_3d%p_patch_1d(1)%zlev_m, &
        & ocean_state(jg)%p_diag)

      IF (timers_level > 2)  CALL timer_stop(timer_extra20)

      CALL perform_accumulation(nnew(1),0)

      CALL output_ocean( patch_3d, ocean_state, &
        &                datetime,              &
        &                surface_fluxes,             &
        &                sea_ice,                 &
        &                jstep, jstep0)
      
      ! receive coupling fluxes for ocean at the end of time stepping loop
      IF (iforc_oce == Coupled_FluxFromAtmo) &  !  14
        &  CALL couple_ocean_toatmo_fluxes(patch_3D, ocean_state(jg), sea_ice, p_atm_f, datetime)

      IF (timers_level > 2)  CALL timer_start(timer_extra21)
      ! Shift time indices for the next loop
      ! this HAS to ge into the restart files, because the start with the following loop
      CALL update_time_indices(jg)
      ! update intermediate timestepping variables for the tracers
      CALL update_time_g_n(ocean_state(jg))

      ! write a restart or checkpoint file
      IF (MOD(jstep,n_checkpoints())==0) THEN
        CALL create_restart_file( patch = patch_2d,       &
          & datetime=datetime,      &
          & jstep=jstep,            &
          & model_type="oce",       &
          & opt_sim_time=time_config%sim_time(1),&
          & opt_nice_class=1,       &
          & ocean_zlevels=n_zlev,                                         &
          & ocean_zheight_cellmiddle = patch_3d%p_patch_1d(1)%zlev_m(:),  &
          & ocean_zheight_cellinterfaces = patch_3d%p_patch_1d(1)%zlev_i(:))
      END IF

      ! check cfl criterion
      IF (cfl_check) THEN
        CALL check_cfl_horizontal(ocean_state(jg)%p_prog(nnew(1))%vn, &
          & patch_2d%edges%inv_dual_edge_length, &
          & dtime, &
          & patch_2d%edges%ALL, &
          & cfl_threshold, &
          & ocean_state(jg)%p_diag%cfl_horz, &
          & cfl_stop_on_violation,&
          & cfl_write)
        CALL check_cfl_vertical(ocean_state(jg)%p_diag%w, &
          & patch_3d%p_patch_1d(1)%prism_center_dist_c, &
          & dtime, &
          & patch_2d%cells%ALL,&
          & cfl_threshold, &
          & ocean_state(jg)%p_diag%cfl_vert, &
          & cfl_stop_on_violation,&
          & cfl_write)
      END IF
      
     IF (timers_level > 2)  CALL timer_stop(timer_extra21)

    ENDDO time_loop

    IF (write_last_restart) &
      & CALL create_restart_file( patch = patch_2d,       &
        & datetime=datetime,      &
        & jstep=jstep,            &
        & model_type="oce",       &
        & opt_sim_time=time_config%sim_time(1),&
        & opt_nice_class=1,       &
        & ocean_zlevels=n_zlev,                                         &
        & ocean_zheight_cellmiddle = patch_3d%p_patch_1d(1)%zlev_m(:),  &
        & ocean_zheight_cellinterfaces = patch_3d%p_patch_1d(1)%zlev_i(:))
       
    CALL timer_stop(timer_total)

  END SUBROUTINE perform_ho_stepping
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
!<Optimize:inUse>
  SUBROUTINE write_initial_ocean_timestep(patch_3d,ocean_state,surface_fluxes,sea_ice, operators_coefficients, p_phys_param)
    TYPE(t_patch_3D), INTENT(IN) :: patch_3d
    TYPE(t_hydro_ocean_state), INTENT(INOUT)    :: ocean_state
    TYPE(t_sfc_flx) , INTENT(INOUT)             :: surface_fluxes
    TYPE(t_sea_ice),          INTENT(INOUT)     :: sea_ice
    TYPE(t_operator_coeff),   INTENT(inout)     :: operators_coefficients    
    TYPE(t_ho_params), INTENT(IN), OPTIONAL     :: p_phys_param
    
    TYPE(t_patch), POINTER :: patch_2d

    patch_2d => patch_3d%p_patch_2d(1)

    ! in general nml output is writen based on the nnew status of the
    ! prognostics variables. Unfortunately, the initialization has to be written
    ! to the nold state. That's why the following manual copying is nec.
    ocean_state%p_prog(nnew(1))%h      = ocean_state%p_prog(nold(1))%h
    
    ocean_state%p_prog(nnew(1))%vn     = ocean_state%p_prog(nold(1))%vn    
    
    CALL calc_scalar_product_veloc_3d( patch_3d,  ocean_state%p_prog(nnew(1))%vn,&
    & ocean_state%p_diag, operators_coefficients)
    ! CALL update_height_depdendent_variables( patch_3d, ocean_state, p_ext_data, operators_coefficients, solvercoeff_sp)
    
    ! copy old tracer values to spot value fields for propper initial timestep
    ! output
    IF(no_tracer>=1)THEN
      ocean_state%p_diag%t = ocean_state%p_prog(nold(1))%tracer(:,:,:,1)
      ! in general nml output is writen based on the nnew status of the
      ! prognostics variables. Unfortunately, the initialization has to be written
      ! to the nold state. That's why the following manual copying is nec.
      ocean_state%p_prog(nnew(1))%tracer = ocean_state%p_prog(nold(1))%tracer
    ENDIF
    IF(no_tracer>=2)THEN
      ocean_state%p_diag%s = ocean_state%p_prog(nold(1))%tracer(:,:,:,2)
    ENDIF
    ocean_state%p_diag%h = ocean_state%p_prog(nold(1))%h
    IF(no_tracer>=1)THEN
      CALL calc_potential_density( patch_3d,                     &
        & ocean_state%p_prog(nold(1))%tracer,&
        & ocean_state%p_diag%rhopot )
      CALL calculate_density( patch_3d,                        &
        & ocean_state%p_prog(nold(1))%tracer, &
        & ocean_state%p_diag%rho )
    ENDIF

    CALL update_ocean_statistics( &
      & ocean_state,            &
      & surface_fluxes,              &
      & patch_2d%cells%owned,   &
      & patch_2d%edges%owned,   &
      & patch_2d%verts%owned,   &
      & n_zlev,p_phys_param=p_phys_param)
    IF (i_sea_ice >= 1) CALL update_ice_statistic(sea_ice%acc, sea_ice,patch_2d%cells%owned)

    CALL write_name_list_output(jstep=0)

    CALL reset_ocean_statistics(ocean_state%p_acc,ocean_state%p_diag,surface_fluxes)
    IF (i_sea_ice >= 1) CALL reset_ice_statistics(sea_ice%acc)

  END SUBROUTINE write_initial_ocean_timestep
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
!<Optimize:inUse>
  SUBROUTINE update_time_g_n(ocean_state)
    TYPE(t_hydro_ocean_state), INTENT(inout) :: ocean_state
    REAL(wp), POINTER ::  tmp(:,:,:)

    ! velocity
    ! just exchange the pointers
    ! ocean_state%p_aux%g_nm1 = ocean_state%p_aux%g_n
    ! ocean_state%p_aux%g_n   = 0.0_wp
    tmp => ocean_state%p_aux%g_n
    ocean_state%p_aux%g_n => ocean_state%p_aux%g_nm1
    ocean_state%p_aux%g_nm1 => tmp
    
  END SUBROUTINE update_time_g_n
  !-------------------------------------------------------------------------

END MODULE mo_hydro_ocean_run
