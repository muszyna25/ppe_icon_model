!>
!! Contains the main stepping method_name the 3-dim hydrostatic ocean model.
!!
!! @author Peter Korn, Stephan Lorenz, MPI
!!
!! @par Revision History
!!  Initial version by Stephan Lorenz (MPI-M), (2010-04).
!!   - renaming and adjustment of hydrostatic ocean model V1.0.3 to ocean domain and patch_oce
!!  Modification by Stephan Lorenz, MPI-M, 2010-10
!!   - new module mo_ocean_testbed_modules including updated reconstructions
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
MODULE mo_ocean_testbed_modules
  !-------------------------------------------------------------------------
  USE mo_kind,                   ONLY: wp
  USE mo_impl_constants,         ONLY: max_char_length
  USE mo_model_domain,           ONLY: t_patch, t_patch_3d, t_subset_range, t_patch_vert
  USE mo_grid_config,            ONLY: n_dom
  USE mo_grid_subset,            ONLY: get_index_range
  USE mo_sync,                   ONLY: sync_patch_array, sync_e, sync_c !, sync_v
  USE mo_ocean_nml,              ONLY: iswm_oce, n_zlev, no_tracer, &
    & diagnostics_level, &
    & eos_type, i_sea_ice, l_staggered_timestep, gibraltar
  USE mo_dynamics_config,        ONLY: nold, nnew
  USE mo_io_config,              ONLY: n_checkpoints
  USE mo_run_config,             ONLY: nsteps, dtime, ltimer, output_mode, test_mode
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
    & t_hydro_ocean_prog, t_ocean_tracer
  USE mo_oce_math_operators,     ONLY: calculate_thickness
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
  USE mo_oce_physics,            ONLY: t_ho_params
  USE mo_oce_thermodyn,          ONLY: calc_density_mpiom_func, calc_density_lin_eos_func,&
    & calc_density_jmdwfg06_eos_func, calc_potential_density, &
    & calc_density
  USE mo_name_list_output,       ONLY: write_name_list_output, istime4name_list_output
  USE mo_oce_diagnostics,        ONLY: calc_slow_oce_diagnostics, calc_fast_oce_diagnostics, &
    & construct_oce_diagnostics,&
    & destruct_oce_diagnostics, t_oce_timeseries, &
    & calc_moc, calc_psi
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
  USE mo_oce_diffusion,          ONLY:  tracer_diffusion_vert_implicit, veloc_diffusion_vert_implicit
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ocean_test_modules
  
  !VERSION CONTROL:
  CHARACTER(LEN=*), PARAMETER :: version = '$Id$'
  CHARACTER(len=12)           :: debug_string = 'testbed     '  ! Output of module for 1 line debug
  
  !-------------------------------------------------------------------------
  
CONTAINS

  !-------------------------------------------------------------------------
  !>
  SUBROUTINE ocean_test_modules( patch_3d, ocean_state, external_data,          &
    & datetime, surface_fluxes, physics_parameters,             &
    & oceans_atmosphere, oceans_atmosphere_fluxes, ocean_ice,operators_coefficients)

    TYPE(t_patch_3d ),TARGET, INTENT(inout)          :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET, INTENT(inout) :: ocean_state(n_dom)
    TYPE(t_external_data), TARGET, INTENT(in)        :: external_data(n_dom)
    TYPE(t_datetime), INTENT(inout)                  :: datetime
    TYPE(t_sfc_flx)                                  :: surface_fluxes
    TYPE (t_ho_params)                               :: physics_parameters
    TYPE(t_atmos_for_ocean),  INTENT(inout)          :: oceans_atmosphere
    TYPE(t_atmos_fluxes ),    INTENT(inout)          :: oceans_atmosphere_fluxes
    TYPE (t_sea_ice),         INTENT(inout)          :: ocean_ice
    TYPE(t_operator_coeff),   INTENT(inout)          :: operators_coefficients

    CHARACTER(LEN=*), PARAMETER ::  method_name = "ocean_test_modules"

    SELECT CASE (test_mode)  !  1 - 99 test ocean modules
      CASE (1)
        CALL ocean_test_advection( patch_3d, ocean_state, external_data,   &
          & datetime, surface_fluxes, physics_parameters,             &
          & oceans_atmosphere, oceans_atmosphere_fluxes, ocean_ice,operators_coefficients)

      CASE (2)
        CALL test_tracer_diffusion_vert_implicit( patch_3d, ocean_state, physics_parameters,  &
           & operators_coefficients)

      CASE (3)
        CALL test_velocity_diffusion_vert_implicit( patch_3d, ocean_state, physics_parameters,  &
           & operators_coefficients)

      CASE DEFAULT
        CALL finish(method_name, "Unknown test_mode")

    END SELECT



  END SUBROUTINE ocean_test_modules
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE test_tracer_diffusion_vert_implicit( patch_3d, ocean_state, physics_parameters,       &
    & operators_coefficients)

    TYPE(t_patch_3d ),TARGET, INTENT(inout)          :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET, INTENT(inout) :: ocean_state(n_dom)
    TYPE (t_ho_params)                               :: physics_parameters
    TYPE(t_operator_coeff),   INTENT(inout)          :: operators_coefficients

    INTEGER :: inner_iter, outer_iter
    TYPE(t_subset_range), POINTER :: cells_in_domain
    TYPE(t_ocean_tracer), POINTER :: ocean_tracer

    ocean_tracer => ocean_state(1)%p_prog(nold(1))%ocean_tracers(1)
    cells_in_domain => patch_3d%p_patch_2D(1)%cells%in_domain
    !---------------------------------------------------------------------

    ocean_tracer%concentration(:,:,:) = 0.0_wp
    ocean_tracer%concentration(:,1,:) = 1.0_wp

    CALL dbg_print('tracer', ocean_tracer%concentration,  debug_string, 1, in_subset=cells_in_domain)
    DO outer_iter=1,1000
      DO inner_iter=1,1000
        CALL tracer_diffusion_vert_implicit( patch_3D, ocean_tracer, ocean_state(1)%p_prog(nold(1))%h, &
          & physics_parameters%A_tracer_v(:,:,:, 1), operators_coefficients)
      ENDDO
      WRITE(message_text,'(i6,a)') outer_iter, 'x1000 iter, tracer'
      CALL dbg_print(message_text, ocean_tracer%concentration,  debug_string, 1, in_subset=cells_in_domain)
    ENDDO

  END SUBROUTINE test_tracer_diffusion_vert_implicit
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE test_velocity_diffusion_vert_implicit( patch_3d, ocean_state, physics_parameters,       &
    & operators_coefficients)

    TYPE(t_patch_3d ),TARGET, INTENT(inout)          :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET, INTENT(inout) :: ocean_state(n_dom)
    TYPE (t_ho_params)                               :: physics_parameters
    TYPE(t_operator_coeff),   INTENT(inout)          :: operators_coefficients

    INTEGER :: inner_iter, outer_iter
    TYPE(t_subset_range), POINTER :: edges_in_domain
    REAL(wp), POINTER :: vn_inout(:,:,:), vn_out(:,:,:)

    edges_in_domain => patch_3d%p_patch_2D(1)%edges%in_domain
    !---------------------------------------------------------------------
    vn_inout  => ocean_state(1)%p_diag%vn_pred
    vn_out    => ocean_state(1)%p_diag%vn_impl_vert_diff
    ! vn_inout(:,1,:) = 1.0
    vn_inout(:,:,:) = 0.0
    vn_inout(:,1,:) = 1.0

    CALL dbg_print('vn', vn_inout,  debug_string, 1, in_subset=edges_in_domain)

    DO outer_iter=1,1000
      DO inner_iter=1,1000
        CALL veloc_diffusion_vert_implicit( patch_3d,  &
          & vn_inout,                                  &
          & physics_parameters%a_veloc_v,              &
          & operators_coefficients,                    &
          & vn_out)
      ENDDO
      WRITE(message_text,'(i6,a)') outer_iter, 'x1000 iter, vn'
      CALL dbg_print(message_text, vn_inout,  debug_string, 1, in_subset=edges_in_domain)
    ENDDO
    return

  END SUBROUTINE test_velocity_diffusion_vert_implicit
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  SUBROUTINE ocean_test_advection( patch_3d, ocean_state, external_data,          &
    & datetime, surface_fluxes, physics_parameters,             &
    & oceans_atmosphere, oceans_atmosphere_fluxes, ocean_ice,operators_coefficients)
    
    TYPE(t_patch_3d ),TARGET, INTENT(inout)          :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET, INTENT(inout) :: ocean_state(n_dom)
    TYPE(t_external_data), TARGET, INTENT(in)        :: external_data(n_dom)
    TYPE(t_datetime), INTENT(inout)                  :: datetime
    TYPE(t_sfc_flx)                                  :: surface_fluxes
    TYPE (t_ho_params)                               :: physics_parameters
    TYPE(t_atmos_for_ocean),  INTENT(inout)          :: oceans_atmosphere
    TYPE(t_atmos_fluxes ),    INTENT(inout)          :: oceans_atmosphere_fluxes
    TYPE (t_sea_ice),         INTENT(inout)          :: ocean_ice
    TYPE(t_operator_coeff),   INTENT(inout)          :: operators_coefficients
    
    ! local variables
    INTEGER :: jstep, jg, jtrc
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
      & method_name = 'mo_ocean_testbed_modules:ocean_test_advection'
    !------------------------------------------------------------------
    
    patch_2D      => patch_3d%p_patch_2d(1)
    CALL datetime_to_string(datestring, datetime)

    ! IF (ltimer) CALL timer_start(timer_total)
    CALL timer_start(timer_total)

    time_config%sim_time(:) = 0.0_wp
    jstep0 = 0
    !------------------------------------------------------------------
    ! IF(.NOT.l_time_marching)THEN

      !IF(itestcase_oce==28)THEN
      DO jstep = (jstep0+1), (jstep0+nsteps)
      
        CALL datetime_to_string(datestring, datetime)
        WRITE(message_text,'(a,i10,2a)') '  Begin of timestep =',jstep,'  datetime:  ', datestring
        CALL message (TRIM(method_name), message_text)
 
!          IF(jstep==1)THEN
!          ocean_state(jg)%p_diag%vn_time_weighted = ocean_state(jg)%p_prog(nold(1))%vn
!          ocean_state(jg)%p_prog(nnew(1))%vn = ocean_state(jg)%p_prog(nold(1))%vn
!          ocean_state(jg)%p_diag%w        =  0.0_wp!0.0833_wp!0.025_wp
!          ocean_state(jg)%p_diag%w(:,:,:) = -0.0833_wp!0.025_wp
!          ENDIF

        !CALL calculate_thickness( patch_3d, ocean_state(jg), external_data(jg))
        !CALL calc_vert_velocity(patch_3d, ocean_state(jg),operators_coefficients)
        CALL advect_tracer_ab( patch_3d, ocean_state(jg),  &
          & physics_parameters,surface_fluxes,&
          & operators_coefficients,&
          & jstep)
        ! One integration cycle finished on the lowest grid level (coarsest
        ! resolution). Set model time.
        CALL add_time(dtime,0,0,0,datetime)
      
        ! Not nice, but the name list output requires this
        time_config%sim_time(1) = time_config%sim_time(1) + dtime
      
        ! update accumulated vars
        CALL update_ocean_statistics(ocean_state(1),&
        & surface_fluxes,                                &
        & patch_2D%cells%owned,       &
        & patch_2D%edges%owned,       &
        & patch_2D%verts%owned,       &
        & n_zlev)
          
        CALL output_ocean( patch_3d, ocean_state, external_data,          &
          & datetime, .false.,                  &
          & surface_fluxes, physics_parameters,             &
          & oceans_atmosphere, oceans_atmosphere_fluxes, ocean_ice,operators_coefficients, &
          & jstep, jstep0)

        ! Shift time indices for the next loop
        ! this HAS to ge into the restart files, because the start with the following loop
        CALL update_time_indices(jg)
        ! update intermediate timestepping variables for the tracers
        ! velocity
        ocean_state(jg)%p_aux%g_nm1 = ocean_state(jg)%p_aux%g_n
        ocean_state(jg)%p_aux%g_n   = 0.0_wp

      END DO
    ! ENDIF!(l_no_time_marching)THEN
    
    CALL timer_stop(timer_total)
    
  END SUBROUTINE ocean_test_advection
  !-------------------------------------------------------------------------
  
  
END MODULE mo_ocean_testbed_modules
