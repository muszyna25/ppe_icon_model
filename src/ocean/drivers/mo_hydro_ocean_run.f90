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
!
!    ProTeX FORTRAN source: Style 2
!    modified for ICON project, DWD/MPI-M 2006
!
!-------------------------------------------------------------------------
!
USE mo_kind,                   ONLY: wp
USE mo_impl_constants,         ONLY: max_char_length
USE mo_model_domain,           ONLY: t_patch, t_patch_3D, t_subset_range
USE mo_grid_config,            ONLY: n_dom, use_dummy_cell_closure
USE mo_grid_subset,            ONLY: get_index_range
USE mo_sync,                   ONLY: sync_patch_array, sync_e!, sync_c, sync_v
USE mo_ocean_nml,              ONLY: iswm_oce, n_zlev, no_tracer, &
  &                                  itestcase_oce, idiag_oce, init_oce_prog, init_oce_relax, &
  &                                  EOS_TYPE, i_sea_ice, l_staggered_timestep, gibraltar
USE mo_dynamics_config,        ONLY: nold, nnew
USE mo_io_config,              ONLY: n_files, n_checkpoints, is_output_time!, istime4newoutputfile
USE mo_run_config,             ONLY: nsteps, dtime, ltimer, output_mode
USE mo_exception,              ONLY: message, message_text, finish
USE mo_ext_data_types,         ONLY: t_external_data
!USE mo_io_units,               ONLY: filename_max
USE mo_datetime,               ONLY: t_datetime, print_datetime, add_time, datetime_to_string
USE mo_timer,                  ONLY: timer_start, timer_stop, timer_total, timer_solve_ab,  &
  &                                  timer_tracer_ab, timer_vert_veloc, timer_normal_veloc, &
  &                                  timer_upd_phys, timer_upd_flx
USE mo_oce_ab_timestepping,    ONLY: solve_free_surface_eq_ab, &
  &                                  calc_normal_velocity_ab,  &
  &                                  calc_vert_velocity,       &
  &                                  update_time_indices
USE mo_oce_init,               ONLY: init_ho_testcases, init_ho_prog, init_ho_coupled,&
  &                                  init_ho_recon_fields, init_ho_relaxation, init_oce_index
USE mo_util_dbg_prnt,          ONLY: init_dbg_index, dbg_print
USE mo_oce_state,              ONLY: t_hydro_ocean_state, t_hydro_ocean_acc, t_hydro_ocean_diag, &
  &                                  t_hydro_ocean_prog, &
  &                                  init_ho_base, init_ho_basins, v_base, &
  &                                  construct_hydro_ocean_base, &! destruct_hydro_ocean_base, &
  &                                  construct_hydro_ocean_state, destruct_hydro_ocean_state, &
  &                                  init_coriolis_oce, init_oce_config, &
  &                                  set_lateral_boundary_values, construct_patch_3D, init_patch_3D, &
  &                                  setup_ocean_namelists, ocean_default_list, ocean_restart_list
USE mo_oce_math_operators,     ONLY: calc_thickness 
USE mo_operator_ocean_coeff_3d,ONLY: t_operator_coeff, allocate_exp_coeff,par_init_operator_coeff,&
  &                                  update_diffusion_matrices
USE mo_scalar_product,         ONLY: calc_scalar_product_veloc_3D
USE mo_oce_tracer,             ONLY: advect_tracer_ab
USE mo_io_restart,             ONLY: write_restart_info_file, create_restart_file
USE mo_oce_bulk,               ONLY: update_sfcflx
USE mo_sea_ice,                ONLY: construct_sfcflx,destruct_sfcflx,&
  &                                  construct_atmos_for_ocean,&
  &                                  destruct_atmos_for_ocean,&
  &                                  construct_atmos_fluxes, destruct_atmos_fluxes,&
  &                                  construct_sea_ice, destruct_sea_ice, ice_init, &
  &                                  update_ice_statistic, compute_mean_ice_statistics, reset_ice_statistics
USE mo_sea_ice_types,          ONLY: t_sfc_flx, t_atmos_fluxes, t_atmos_for_ocean, &
  &                                  t_sea_ice
USE mo_oce_forcing,            ONLY: init_sfcflx
USE mo_oce_physics,            ONLY: t_ho_params, &
  &                                  construct_ho_params, init_ho_params, &
  &                                  destruct_ho_params, update_ho_params
USE mo_oce_thermodyn,          ONLY: calc_density_MPIOM_func, calc_density_lin_EOS_func,&
  &                                  calc_density_JMDWFG06_EOS_func, calc_potential_density, &
  &                                  calc_density
USE mo_name_list_output,       ONLY: write_name_list_output, istime4name_list_output
USE mo_oce_diagnostics,        ONLY: calculate_oce_diagnostics,&
  &                                  construct_oce_diagnostics,&
  &                                  destruct_oce_diagnostics, t_oce_timeseries, &
  &                                  calc_moc, calc_psi
USE mo_oce_ab_timestepping_mimetic, ONLY: init_ho_lhs_fields_mimetic
USE mo_linked_list,            ONLY: t_list_element, find_list_element
USE mo_var_list,               ONLY: print_var_list
USE mo_io_restart_attributes,  ONLY: get_restart_attribute
USE mo_mpi,                    ONLY: my_process_is_stdio
USE mo_time_config,            ONLY: time_config
USE mo_master_control,         ONLY: is_restart_run
USE mo_statistics
USE mo_grid_tools,             ONLY: create_dummy_cell_closure
USE mo_sea_ice_nml,            ONLY: i_ice_dyn
USE mo_ocean_nml,              ONLY: i_sea_ice

IMPLICIT NONE

PRIVATE
INTEGER, PARAMETER :: kice = 1

!VERSION CONTROL:
CHARACTER(LEN=*), PARAMETER :: version = '$Id$'

! public interface
!
! public subroutines
PUBLIC :: perform_ho_stepping
PUBLIC :: prepare_ho_stepping
PUBLIC :: construct_ocean_states
PUBLIC :: finalise_ho_integration
PRIVATE:: update_intermediate_tracer_vars
!
CHARACTER(len=12)  :: str_module = 'HYDRO-ocerun'  ! Output of module for 1 line debug
INTEGER            :: idt_src                      ! Level of detail for 1 line debug
!
!-------------------------------------------------------------------------

CONTAINS
  SUBROUTINE prepare_ho_stepping(p_patch_3D, operators_coefficients, ocean_state, is_restart)
    TYPE(t_patch_3D ), INTENT(IN)     :: p_patch_3D
    TYPE(t_hydro_ocean_state), TARGET :: ocean_state
    LOGICAL, INTENT(IN)               :: is_restart
    TYPE(t_operator_coeff)            :: operators_coefficients

    IF (is_restart) THEN
      ! Prepare p_os%p_prog, since it is needed by the dynamical sea ice model
      IF ( i_sea_ice > 0 .and. i_ice_dyn == 1 )         &
        CALL calc_scalar_product_veloc_3D( p_patch_3D,  &
          & ocean_state%p_prog(nnew(1))%vn,             &
          & ocean_state%p_prog(nnew(1))%vn,             &
          & ocean_state%p_diag,                         &
          & operators_coefficients)
    ELSE
    ENDIF
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
  SUBROUTINE perform_ho_stepping( patch_3D, p_os, p_ext_data,          &
                                & datetime, lwrite_restart,            &
                                & p_sfc_flx, p_phys_param,             &
                                & p_as, p_atm_f, p_ice,p_op_coeff)

  TYPE(t_patch_3D ),TARGET, INTENT(INOUT)          :: patch_3D
  TYPE(t_hydro_ocean_state), TARGET, INTENT(INOUT) :: p_os(n_dom)
  TYPE(t_external_data), TARGET, INTENT(IN)        :: p_ext_data(n_dom)
  TYPE(t_datetime), INTENT(INOUT)                  :: datetime
  LOGICAL, INTENT(IN)                              :: lwrite_restart
  TYPE(t_sfc_flx)                                  :: p_sfc_flx
  TYPE (t_ho_params)                               :: p_phys_param
  TYPE(t_atmos_for_ocean),  INTENT(INOUT)          :: p_as
  TYPE(t_atmos_fluxes ),    INTENT(INOUT)          :: p_atm_f
  TYPE (t_sea_ice),         INTENT(INOUT)          :: p_ice
  TYPE(t_operator_coeff),   INTENT(INOUT)          :: p_op_coeff

  ! local variables
  INTEGER                         :: jstep, jg, jtrc
  INTEGER                         :: nsteps_since_last_output
  INTEGER                         :: ocean_statistics
  !LOGICAL                         :: l_outputtime
  CHARACTER(len=32)               :: datestring, plaindatestring
  TYPE(t_oce_timeseries), POINTER :: oce_ts
  TYPE(t_patch), POINTER          :: patch_2D
  INTEGER                         :: jstep0 ! start counter for time loop

  !CHARACTER(LEN=filename_max)  :: outputfile, gridfile
  CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
    &      routine = 'mo_hydro_ocean_run:perform_ho_stepping'
  !------------------------------------------------------------------

  nsteps_since_last_output = 1

  !------------------------------------------------------------------
  ! no grid refinement allowed here so far
  !------------------------------------------------------------------
  IF (n_dom > 1 ) THEN
    CALL finish(TRIM(routine), ' N_DOM > 1 is not allowed')
  END IF
  jg = n_dom

  patch_2D => patch_3D%p_patch_2D(jg)

  CALL datetime_to_string(datestring, datetime)

  IF (idiag_oce == 1) &
    & CALL construct_oce_diagnostics( patch_3D, p_os(jg), oce_ts, datestring)

  ! IF (ltimer) CALL timer_start(timer_total)
  CALL timer_start(timer_total)

  time_config%sim_time(:) = 0.0_wp

  !------------------------------------------------------------------
  ocean_statistics = new_statistic()

  !------------------------------------------------------------------
  ! write initial
  !------------------------------------------------------------------
  IF (output_mode%l_nml) THEN
    ! in general nml output is writen based on the nnew status of the
    ! prognostics variables. Unfortunately, the initialization has to be written
    ! to the nold state. That's why the following manual copying is nec.
    IF (.NOT. is_restart_run()) THEN
      p_os(jg)%p_prog(nnew(1))%tracer = p_os(jg)%p_prog(nold(1))%tracer
      ! copy old tracer values to spot value fields for propper initial timestep
      ! output
      p_os(jg)%p_diag%t = p_os(jg)%p_prog(nold(1))%tracer(:,:,:,1)
      p_os(jg)%p_diag%s = p_os(jg)%p_prog(nold(1))%tracer(:,:,:,2)
      p_os(jg)%p_diag%h = p_os(jg)%p_prog(nold(1))%h
      CALL calc_potential_density( patch_3D,                     &
        &                          p_os(jg)%p_prog(nold(1))%tracer,&
        &                          p_os(jg)%p_diag%rhopot )
      CALL calc_density( patch_3D,                        &
        &                p_os(jg)%p_prog(nold(1))%tracer, &
        &                p_os(jg)%p_diag%rho )

      p_os(jg)%p_acc%tracer(:,:,:,1)  = p_os(jg)%p_prog(nold(1))%tracer(:,:,:,1)
      p_os(jg)%p_acc%tracer(:,:,:,2)  = p_os(jg)%p_prog(nold(1))%tracer(:,:,:,2)
      p_os(jg)%p_acc%h                = p_os(jg)%p_prog(nold(1))%h
      p_os(jg)%p_acc%rhopot           = p_os(jg)%p_diag%rhopot
      p_os(jg)%p_acc%rho              = p_os(jg)%p_diag%rho
      IF (i_sea_ice >= 1) THEN
        CALL update_ice_statistic(p_ice%acc, p_ice,patch_3D%p_patch_2D(1)%cells%owned)
      ENDIF
    ENDIF
    CALL write_name_list_output(jstep=0)
    IF (.NOT. is_restart_run()) THEN
      ! reset the accs to zero
      p_os(jg)%p_acc%tracer(:,:,:,1)  = 0.0_wp
      p_os(jg)%p_acc%tracer(:,:,:,2)  = 0.0_wp
      p_os(jg)%p_acc%h                = 0.0_wp
      p_os(jg)%p_acc%rhopot          = 0.0_wp
      p_os(jg)%p_acc%rho             = 0.0_wp
      IF (i_sea_ice >= 1) THEN
        CALL reset_ice_statistics(p_ice%acc)
      ENDIF
    ENDIF
  ENDIF
  !------------------------------------------------------------------
  ! call the dynamical core: start the time loop
  !------------------------------------------------------------------

  jstep0 = 0
  IF (is_restart_run()) THEN
    ! get start counter for time loop from restart file:
    CALL get_restart_attribute("jstep", jstep0)
  END IF

  TIME_LOOP: DO jstep = (jstep0+1), (jstep0+nsteps)

    CALL datetime_to_string(datestring, datetime)
    WRITE(message_text,'(a,i10,2a)') '  Begin of timestep =',jstep,'  datetime:  ', datestring
    CALL message (TRIM(routine), message_text)

    IF(itestcase_oce==28)THEN
      CALL calc_thickness( patch_3D, p_os(jg), p_ext_data(jg))
      CALL calc_vert_velocity(patch_3D, p_os(jg),p_op_coeff)
      CALL advect_tracer_ab( patch_3D, p_os(jg),  &
                           & p_phys_param,p_sfc_flx,&
                           & p_op_coeff,&
                           & jstep)
    ELSE
      !In case of a time-varying forcing:
      IF (ltimer) CALL timer_start(timer_upd_flx)
      CALL update_sfcflx( patch_3D, p_os(jg), p_as, p_ice, p_atm_f, p_sfc_flx, &
        &                jstep, datetime, p_op_coeff)

      IF(.NOT.l_staggered_timestep)THEN

        CALL calc_thickness( patch_3D, p_os(jg), p_ext_data(jg))

        CALL set_lateral_boundary_values( patch_3D, p_os(jg)%p_prog(nold(1))%vn)
        CALL sync_patch_array(sync_e, patch_3D%p_patch_2D(jg), p_os(jg)%p_prog(nold(1))%vn)

        CALL calc_scalar_product_veloc_3D( patch_3D, &
          & p_os(jg)%p_prog(nold(1))%vn,         &
          & p_os(jg)%p_prog(nold(1))%vn,         &
          & p_os(jg)%p_diag,                     &
          & p_op_coeff)

        ! activate for calc_scalar_product_veloc_3D
        !---------DEBUG DIAGNOSTICS-------------------------------------------
        idt_src=2  ! output print level (1-5, fix)
        CALL dbg_print('on entry: h-old'           ,p_os(jg)%p_prog(nold(1))%h ,str_module,idt_src, &
          patch_2D%cells%owned )
        idt_src=1  ! output print level (1-5, fix)
        CALL dbg_print('on entry: h-new'           ,p_os(jg)%p_prog(nnew(1))%h ,str_module,idt_src, &
          patch_2D%cells%owned )
        idt_src=3  ! output print level (1-5, fix)
        CALL dbg_print('HydOce: ScaProdVel kin'    ,p_os(jg)%p_diag%kin        ,str_module,idt_src, &
          patch_2D%cells%owned )
        CALL dbg_print('HydOce: ScaProdVel ptp_vn' ,p_os(jg)%p_diag%ptp_vn     ,str_module,idt_src, &
          patch_2D%edges%owned )
        !---------------------------------------------------------------------

      ENDIF
      IF (ltimer) CALL timer_stop(timer_upd_flx)

      IF (ltimer) CALL timer_start(timer_upd_phys)

      SELECT CASE (EOS_TYPE)
        CASE(1)
          CALL update_ho_params(patch_3D, p_os(jg), p_sfc_flx, p_phys_param,&
            &                   calc_density_lin_EOS_func)
        CASE(2)
          CALL update_ho_params(patch_3D, p_os(jg), p_sfc_flx, p_phys_param,&
            &                   calc_density_MPIOM_func)
        CASE(3)
          CALL update_ho_params(patch_3D,p_os(jg), p_sfc_flx, p_phys_param,&
            &                   calc_density_JMDWFG06_EOS_func)
      CASE DEFAULT
      END SELECT
      IF (ltimer) CALL timer_stop(timer_upd_phys)

      CALL update_diffusion_matrices( patch_3D,                   &
                                    & p_os(jg),                     &
                                    & p_phys_param,                 &
                                    & p_op_coeff%matrix_vert_diff_e,&
                                    & p_op_coeff%matrix_vert_diff_c)

      !------------------------------------------------------------------------
      ! solve for new free surface
      IF (ltimer) CALL timer_start(timer_solve_ab)
      CALL solve_free_surface_eq_ab (patch_3D, p_os(jg), p_ext_data(jg), &
        &                            p_sfc_flx, p_phys_param, jstep, p_op_coeff)!, p_int(jg))
      IF (ltimer) CALL timer_stop(timer_solve_ab)

      !------------------------------------------------------------------------
      ! Step 4: calculate final normal velocity from predicted horizontal
      ! velocity vn_pred and updated surface height
      IF (ltimer) CALL timer_start(timer_normal_veloc)
      CALL calc_normal_velocity_ab(patch_3D, p_os(jg),&
                                  &p_op_coeff, p_ext_data(jg), p_phys_param)
      IF (ltimer) CALL timer_stop(timer_normal_veloc)

      !------------------------------------------------------------------------
      ! Step 5: calculate vertical velocity from continuity equation under
      ! incompressiblity condition in the non-shallow-water case
      IF ( iswm_oce /= 1 ) THEN
        IF (ltimer) CALL timer_start(timer_vert_veloc)
        CALL calc_vert_velocity( patch_3D, p_os(jg),p_op_coeff)
        IF (ltimer) CALL timer_stop(timer_vert_veloc)
      ENDIF

      !------------------------------------------------------------------------
      ! Step 6: transport tracers and diffuse them
      IF (no_tracer>=1) THEN
        IF (ltimer) CALL timer_start(timer_tracer_ab)
        CALL advect_tracer_ab( patch_3D, p_os(jg), p_phys_param,&
                             & p_sfc_flx,&
                             & p_op_coeff,&
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
    CALL calc_potential_density( patch_3D,                       &
      &                          p_os(jg)%p_prog(nold(1))%tracer,&
      &                          p_os(jg)%p_diag%rhopot )
    CALL calc_psi (patch_2D,patch_3D, p_os(jg)%p_diag%u(:,:,:), &
      &            p_os(jg)%p_prog(nold(1))%h(:,:),             &
      &            p_os(jg)%p_diag%u_vint, datetime)

    ! update accumulated vars
    CALL update_ocean_statistics(p_os(1),                              &
      &                          p_sfc_flx,                            &
      &                          patch_3D%p_patch_2D(1)%cells%owned,   &
      &                          patch_3D%p_patch_2D(1)%edges%owned,   &
      &                          patch_3D%p_patch_2D(1)%verts%owned,   &
      &                          n_zlev)
    CALL update_ice_statistic(p_ice%acc,p_ice,patch_3D%p_patch_2D(1)%cells%owned)


    IF (istime4name_list_output(jstep)) THEN
      IF (idiag_oce == 1 ) THEN
        CALL calculate_oce_diagnostics( patch_3D,    &
          &                             p_os(jg),      &
          &                             p_sfc_flx,     &
          &                             p_ice,         &
          &                             p_phys_param,  &
          &                             jstep-jstep0,         &
          &                             datetime,      &
          &                             oce_ts)

        CALL calc_moc (patch_2D,patch_3D, p_os(jg)%p_diag%w(:,:,:), datetime)

      ENDIF
      ! compute mean values for output interval
      !TODO [ram] src/io/shared/mo_output_event_types.f90 for types to use
      !TODO [ram] nsteps_since_last_output =
      !TODO [ram] output_event%event_step(output_event%i_event_step)%i_sim_step - output_event%event_step(output_event%i_event_step-1)%i_sim_step

      CALL compute_mean_ocean_statistics(p_os(1)%p_acc,p_sfc_flx,nsteps_since_last_output)
      CALL compute_mean_ice_statistics(p_ice%acc,nsteps_since_last_output)

      ! set the output variable pointer to the correct timelevel
      CALL set_output_pointers(nnew(1), p_os(jg)%p_diag, p_os(jg)%p_prog(nnew(1)))

      IF (output_mode%l_nml) THEN
        CALL write_name_list_output(jstep)
      ENDIF

      CALL message (TRIM(routine),'Write output at:')
      CALL print_datetime(datetime)

      ! reset accumulation vars
      CALL reset_ocean_statistics(p_os(1)%p_acc,p_sfc_flx,nsteps_since_last_output)
      CALL reset_ice_statistics(p_ice%acc)

    END IF

    ! Shift time indices for the next loop
    ! this HAS to ge into the restart files, because the start with the following loop
    CALL update_time_indices(jg)
    ! update intermediate timestepping variables for the tracers
    CALL update_intermediate_tracer_vars(p_os(jg))

    ! write a restart or checkpoint file
    IF (MOD(jstep,n_checkpoints())==0 .OR. ((jstep==(jstep0+nsteps)) .AND. lwrite_restart)) THEN
      CALL create_restart_file( patch_2D, datetime,                             &
                              & jstep, opt_depth=n_zlev,                       &
                              & opt_sim_time=time_config%sim_time(1),          &
                              & opt_nice_class=1)
      ! Create the master (meta) file in ASCII format which contains
      ! info about which files should be read in for a restart run.
      CALL write_restart_info_file
    END IF

    nsteps_since_last_output = nsteps_since_last_output + 1
  ENDDO TIME_LOOP

  IF (idiag_oce==1) CALL destruct_oce_diagnostics(oce_ts)
  CALL delete_statistic(ocean_statistics)

  CALL timer_stop(timer_total)

  END SUBROUTINE perform_ho_stepping
 !-------------------------------------------------------------------------
  !>
  !! Simple routine for preparing hydrostatic ocean model.
  !!
  !! Simple routine for preparing hydrostatic ocean model.
  !! Calls basic routines ...
  !!
  !!
  !! @par Revision History
  !! Initial release by Stephan Lorenz, MPI-M (2010-07)
  SUBROUTINE construct_ocean_states(patch_3D, p_os, p_ext_data, p_sfc_flx, &
                                  & p_phys_param, p_as,&
                                  & p_atm_f, p_ice, p_op_coeff)

    TYPE(t_patch_3D ),TARGET,   INTENT(INOUT)  :: patch_3D
    TYPE(t_hydro_ocean_state),  INTENT(INOUT)  :: p_os(n_dom)
    TYPE(t_external_data),      INTENT(INOUT)  :: p_ext_data(n_dom)
    TYPE(t_sfc_flx),            INTENT(INOUT)  :: p_sfc_flx
    TYPE(t_ho_params),          INTENT(INOUT)  :: p_phys_param
    TYPE(t_atmos_for_ocean ),   INTENT(INOUT)  :: p_as
    TYPE(t_atmos_fluxes ),      INTENT(INOUT)  :: p_atm_f
    TYPE(t_sea_ice),            INTENT(INOUT)  :: p_ice
    TYPE(t_operator_coeff),     INTENT(INOUT)  :: p_op_coeff

    ! local variables
    INTEGER :: jg
    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
      &      routine = 'mo_test_hydro_ocean:construct_ocean_states'

    CALL message (TRIM(routine),'start')
    !------------------------------------------------------------------
    ! no grid refinement allowed here so far
    !------------------------------------------------------------------

    IF (n_dom > 1 ) THEN
      CALL finish(TRIM(routine), ' N_DOM > 1 is not allowed')
    END IF
    jg = n_dom

    !------------------------------------------------------------------
    ! construct ocean state and physics
    !------------------------------------------------------------------
    CALL init_oce_config

    ! initialize ocean indices for debug output (before ocean state, no 3-dim)
    CALL init_dbg_index(patch_3D%p_patch_2D(jg))!(patch_2D(jg))

    ! hydro_ocean_base contains the 3-dimensional structures for the ocean state

    CALL construct_patch_3D(patch_3D)

    CALL construct_hydro_ocean_base(patch_3D%p_patch_2D(jg), v_base)
    CALL init_ho_base     (patch_3D%p_patch_2D(jg), p_ext_data(jg), v_base)
    CALL init_ho_basins   (patch_3D%p_patch_2D(jg),                 v_base)
    CALL init_coriolis_oce(patch_3D%p_patch_2D(jg) )
    CALL init_patch_3D    (patch_3D,                p_ext_data(jg), v_base)
    !CALL init_patch_3D(patch_3D, v_base)

    !------------------------------------------------------------------
    ! construct ocean state and physics
    !------------------------------------------------------------------

    ! patch_2D and p_os have dimension n_dom
    CALL construct_hydro_ocean_state(patch_3D%p_patch_2D, p_os)

    ! initialize ocean indices for debug output (including 3-dim lsm)
    CALL init_oce_index( patch_3D%p_patch_2D,patch_3D, p_os, p_ext_data )

    CALL construct_ho_params(patch_3D%p_patch_2D(jg), p_phys_param)
    CALL init_ho_params(patch_3D, p_phys_param)

    !------------------------------------------------------------------
    ! construct ocean forcing and testcases
    !------------------------------------------------------------------

    CALL construct_sfcflx(patch_3D%p_patch_2D(jg),p_sfc_flx, ocean_default_list)
    CALL      init_sfcflx(patch_3D, p_sfc_flx)

    CALL construct_sea_ice(patch_3D, p_ice, kice)
    CALL construct_atmos_for_ocean(patch_3D%p_patch_2D(jg), p_as)
    CALL construct_atmos_fluxes(patch_3D%p_patch_2D(jg), p_atm_f, kice)

    IF (init_oce_prog == 0) THEN
      CALL init_ho_testcases(patch_3D%p_patch_2D(jg),patch_3D, p_os(jg), p_ext_data(jg), p_op_coeff,p_sfc_flx)
    ELSE IF (init_oce_prog == 1) THEN

      CALL init_ho_prog(patch_3D%p_patch_2D(jg),patch_3D, p_os(jg), p_sfc_flx)
    END IF

    IF (init_oce_relax == 1) THEN
      CALL init_ho_relaxation(patch_3D%p_patch_2D(jg),patch_3D, p_os(jg), p_sfc_flx)
    END IF

    CALL init_ho_coupled(patch_3D%p_patch_2D(jg), p_os(jg))
    IF (i_sea_ice >= 1) &
      &   CALL ice_init(patch_3D, p_os(jg), p_ice)

    CALL allocate_exp_coeff     ( patch_3D%p_patch_2D(jg), p_op_coeff, ocean_default_list)
    CALL par_init_operator_coeff( patch_3D, p_os(jg),p_phys_param, p_op_coeff)
    CALL init_ho_recon_fields   ( patch_3D%p_patch_2D(jg),patch_3D, p_os(jg), p_op_coeff)

    CALL init_ho_lhs_fields_mimetic   ( patch_3D )

    IF (use_dummy_cell_closure) CALL create_dummy_cell_closure(patch_3D)

    CALL message (TRIM(routine),'end')

  END SUBROUTINE construct_ocean_states

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
  SUBROUTINE finalise_ho_integration(p_os, p_phys_param, p_as, p_atm_f, p_ice, p_sfc_flx)
    TYPE(t_hydro_ocean_state), INTENT(INOUT) :: p_os(n_dom)
    TYPE (t_ho_params),        INTENT(INOUT) :: p_phys_param
    TYPE(t_atmos_for_ocean),   INTENT(INOUT) :: p_as
    TYPE(t_atmos_fluxes ),     INTENT(INOUT) :: p_atm_f
    TYPE (t_sea_ice),          INTENT(INOUT) :: p_ice
    TYPE(t_sfc_flx),           INTENT(INOUT) :: p_sfc_flx


    !------------------------------------------------------------------
    ! destruct ocean physics and forcing
    ! destruct ocean state is in control_model
    !------------------------------------------------------------------
     CALL destruct_hydro_ocean_state(p_os)
     !CALL destruct_hydro_ocean_base(v_base)
     CALL destruct_ho_params(p_phys_param)

     IF(no_tracer>0) CALL destruct_sfcflx(p_sfc_flx)
     CALL destruct_sea_ice(p_ice)
     CALL destruct_atmos_for_ocean(p_as)
     CALL destruct_atmos_fluxes(p_atm_f)


  END SUBROUTINE finalise_ho_integration

  SUBROUTINE update_intermediate_tracer_vars(p_os)
    TYPE(t_hydro_ocean_state), INTENT(INOUT) :: p_os

    ! velocity
    p_os%p_aux%g_nm1 = p_os%p_aux%g_n
    p_os%p_aux%g_n   = 0.0_wp
  END SUBROUTINE update_intermediate_tracer_vars

  SUBROUTINE update_ocean_statistics(p_os,p_sfc_flx,cells,edges,verts,max_zlev)
    TYPE(t_hydro_ocean_state), INTENT(INOUT) :: p_os
    TYPE(t_sfc_flx),           INTENT(INOUT) :: p_sfc_flx
    TYPE(t_subset_range),      INTENT(IN)    :: cells,edges,verts
    INTEGER, INTENT(IN)                      :: max_zlev

    INTEGER :: jtrc,i


    ! update ocean state accumulated values
    CALL add_fields(p_os%p_acc%h     , p_os%p_prog(nnew(1))%h, cells)
    CALL add_fields(p_os%p_acc%u     , p_os%p_diag%u         , cells)
    CALL add_fields(p_os%p_acc%v     , p_os%p_diag%v         , cells)
    CALL add_fields(p_os%p_acc%rhopot, p_os%p_diag%rhopot    , cells)
    DO jtrc=1,no_tracer
    CALL add_fields(p_os%p_acc%tracer(:,:,:,jtrc),           &
      &             p_os%p_prog(nnew(1))%tracer(:,:,:,jtrc), &
      &             cells)
    END DO
    CALL add_fields(p_os%p_acc%u_vint        , p_os%p_diag%u_vint        , cells)
    CALL add_fields(p_os%p_acc%w             , p_os%p_diag%w             , cells,max_zlev+1)
    CALL add_fields(p_os%p_acc%div_mass_flx_c, p_os%p_diag%div_mass_flx_c, cells)
    CALL add_fields(p_os%p_acc%rho           , p_os%p_diag%rho           , cells)
    CALL add_fields(p_os%p_acc%vt            , p_os%p_diag%vt            , edges)
    CALL add_fields(p_os%p_acc%mass_flx_e    , p_os%p_diag%mass_flx_e    , edges)
    CALL add_fields(p_os%p_acc%vort          , p_os%p_diag%vort          , verts,max_zlev)
    CALL add_fields(p_os%p_acc%kin           , p_os%p_diag%kin           , cells)

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
    TYPE(t_hydro_ocean_acc), INTENT(INOUT) :: p_acc
    TYPE(t_sfc_flx),         INTENT(INOUT) :: p_sfc_flx
    INTEGER,INTENT(IN)                     :: nsteps_since_last_output


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
    p_sfc_flx%forc_tracer_acc       = p_sfc_flx%forc_tracer_acc      /REAL(nsteps_since_last_output,wp)
    p_sfc_flx%forc_tracer_relax_acc = p_sfc_flx%forc_tracer_relax_acc/REAL(nsteps_since_last_output,wp)
  END SUBROUTINE compute_mean_ocean_statistics

  SUBROUTINE reset_ocean_statistics(p_acc,p_sfc_flx,nsteps_since_last_output)
    TYPE(t_hydro_ocean_acc), INTENT(INOUT) :: p_acc
    TYPE(t_sfc_flx),         INTENT(INOUT) :: p_sfc_flx
    INTEGER,                 INTENT(INOUT) :: nsteps_since_last_output

    nsteps_since_last_output        = 0
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
    p_sfc_flx%forc_tracer_acc       = 0.0_wp
    p_sfc_flx%forc_tracer_relax_acc = 0.0_wp

  END SUBROUTINE reset_ocean_statistics

  SUBROUTINE new_ocean_statistics()
  END SUBROUTINE new_ocean_statistics

  SUBROUTINE set_output_pointers(timelevel,p_diag,p_prog)
    INTEGER, INTENT(IN) :: timelevel
    TYPE(t_hydro_ocean_diag) :: p_diag
    TYPE(t_hydro_ocean_prog) :: p_prog

    TYPE(t_list_element), POINTER  :: output_var => null()
    TYPE(t_list_element), POINTER  :: prog_var   => null()
    CHARACTER(len=max_char_length) :: timelevel_str
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
    p_diag%t(:,:,:)        =  p_prog%tracer(:,:,:,1)

   !output_var             => find_list_element(ocean_restart_list,'s')
   !prog_var               => find_list_element(ocean_restart_list,'s'//TRIM(timelevel_str))
   !output_var%field%r_ptr => prog_var%field%r_ptr
   !p_diag%s               => prog_var%field%r_ptr(:,:,:,1,1)
    p_diag%s(:,:,:)        =  p_prog%tracer(:,:,:,2)
  END SUBROUTINE set_output_pointers

END MODULE mo_hydro_ocean_run
