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
  USE mo_kind,                   ONLY: wp, sp
  USE mo_impl_constants,         ONLY: max_char_length
  USE mo_model_domain,           ONLY: t_patch, t_patch_3d, t_subset_range, t_patch_vert
  USE mo_grid_config,            ONLY: n_dom
  USE mo_grid_subset,            ONLY: get_index_range
  USE mo_sync,                   ONLY: sync_patch_array, sync_e, sync_c !, sync_v
  USE mo_ocean_nml,              ONLY: iswm_oce, n_zlev, no_tracer, &
    & diagnostics_level, use_tracer_x_height, &
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
  USE mo_parallel_config,        ONLY: nproma
  USE mo_math_utility_solvers,   ONLY: apply_triangular_matrix
  USE mo_ocean_initial_conditions, ONLY: fill_tracer_x_height

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ocean_test_modules
  
  !VERSION CONTROL:
  CHARACTER(LEN=*), PARAMETER :: version = '$Id$'
  CHARACTER(len=12)           :: debug_string = 'testbed     '  ! Output of module for 1 line debug
  
  !-------------------------------------------------------------------------
  INTEGER :: vertical_diffusion_resIterations = 0
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

    CALL calculate_thickness( patch_3D, ocean_state(1), external_data(1), operators_coefficients)

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
    REAL(wp) :: residual(1:nproma,1:n_zlev,1:patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    REAL(wp) :: sum_trace

    ocean_tracer => ocean_state(1)%p_prog(nold(1))%ocean_tracers(1)
    cells_in_domain => patch_3d%p_patch_2D(1)%cells%in_domain

    !---------------------------------------------------------------------
    ! ocean_tracer%concentration(:,:,:) = 0.0_wp
    ! ocean_tracer%concentration(:,17,:) = 1.0_wp

    CALL dbg_print('tracer', ocean_tracer%concentration,  debug_string, 1, in_subset=cells_in_domain)

    CALL fill_tracer_x_height(patch_3d, ocean_state(1))
!    write(0,*) ocean_tracer%concentration_x_height(:, :, 1)
    sum_trace = SUM(ocean_tracer%concentration_x_height(1, :, 1))
    WRITE(message_text,'(f18.10)') sum_trace
    CALL message("sum=", message_text)

    DO outer_iter=1,1000
      DO inner_iter=1,1000
        CALL tracer_diffusion_vert_implicit( patch_3D, ocean_tracer, ocean_state(1)%p_prog(nold(1))%h, &
          & physics_parameters%A_tracer_v(:,:,:, 1), operators_coefficients) !, residual)
      ENDDO

      WRITE(message_text,'(i6,a)') outer_iter, 'x1000 iter, tracer'
      CALL dbg_print(message_text, ocean_tracer%concentration,  debug_string, 1, in_subset=cells_in_domain)
      CALL fill_tracer_x_height(patch_3d, ocean_state(1))
      sum_trace = SUM(ocean_tracer%concentration_x_height(1, :, 1))
      WRITE(message_text,'(f18.10)') sum_trace
      CALL message("sum=", message_text)
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
    REAL(wp) :: sum_vn

    edges_in_domain => patch_3d%p_patch_2D(1)%edges%in_domain
    !---------------------------------------------------------------------
    vn_inout  => ocean_state(1)%p_diag%vn_pred
    vn_out    => ocean_state(1)%p_diag%vn_impl_vert_diff
     vn_inout(:,:,:) = 10.0
    !vn_inout(:,1,:) = 1.0

    CALL dbg_print('vn', vn_inout,  debug_string, 1, in_subset=edges_in_domain)
    sum_vn = SUM(vn_inout(1, :, 1) * patch_3d%p_patch_1d(1)%prism_thick_flat_sfc_e(1, :, 1))
    WRITE(message_text,'(f18.10)') sum_vn
    CALL message("sum=", message_text)

    DO outer_iter=1,1000
      DO inner_iter=1,1000

        CALL update_diffusion_matrices( patch_3d,     &
          & physics_parameters,                       &
          & operators_coefficients%matrix_vert_diff_e,&
          & operators_coefficients%matrix_vert_diff_c)

        CALL veloc_diffusion_vert_implicit( patch_3d,  &
          & vn_inout,                                  &
          & physics_parameters%a_veloc_v,              &
          & operators_coefficients,                    &
          & vn_out)
      ENDDO
      WRITE(message_text,'(i6,a)') outer_iter, 'x1000 iter, vn'
      CALL dbg_print(message_text, vn_inout,  debug_string, 1, in_subset=edges_in_domain)
      sum_vn = SUM(vn_inout(1, :, 1) * patch_3d%p_patch_1d(1)%prism_thick_flat_sfc_e(1, :, 1))
      WRITE(message_text,'(f18.10)') sum_vn
      CALL message("sum=", message_text)

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
  

  !-------------------------------------------------------------------------
  SUBROUTINE tracer_diffusion_vert_implicit_r1( patch_3D,               &
                                           & ocean_tracer, h_c, A_v,   &
                                           & operators_coefficients, &
                                           & residual_3D  ) !,  &
                                          ! & diff_column)

    TYPE(t_patch_3D ),TARGET, INTENT(IN) :: patch_3D
    TYPE(t_ocean_tracer), TARGET         :: ocean_tracer
    REAL(wp), INTENT(inout)              :: h_c(:,:)  !surface height, relevant for thickness of first cell, in
    REAL(wp), INTENT(inout)              :: A_v(:,:,:)
    TYPE(t_operator_coeff),TARGET     :: operators_coefficients
    REAL(wp), INTENT(inout) :: residual_3D(1:nproma,1:n_zlev,1:patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    ! REAL(wp), INTENT(inout)           :: diff_column(1:nproma,1:n_zlev,1:patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    !
    !Local variables
    INTEGER :: slev
    INTEGER :: jc, jk, jb
    INTEGER :: i_startidx_c, i_endidx_c
    REAL(wp) :: a(1:n_zlev), b(1:n_zlev), c(1:n_zlev)
    REAL(wp) :: fact(1:n_zlev)
    ! REAL(wp) :: inv_prisms_center_distance(1:n_zlev)
    ! REAL(wp) :: inv_prism_thickness(1:n_zlev)
    REAL(wp) :: dt_inv, dt
    REAL(wp), POINTER   :: field_column(:,:,:)
    REAL(wp) :: column_tracer(1:n_zlev), residual(1:n_zlev), old_tracer(1:n_zlev)
    REAL(wp) :: residual_fraction(1:n_zlev)
    REAL(wp) :: prism_thickness(1:n_zlev), inv_prism_thickness(1:n_zlev), inv_prisms_center_distance(1:n_zlev)
    REAL(wp) :: prisms_center_distance(1:n_zlev), unit(1:n_zlev)
    INTEGER  :: z_dolic
    INTEGER  :: resIter
    TYPE(t_subset_range), POINTER :: cells_in_domain
    TYPE(t_patch), POINTER         :: patch_2D
    ! CHARACTER(len=max_char_length), PARAMETER :: &
    !        & routine = ('mo_oce_diffusion:tracer_diffusion_impl')
    !-----------------------------------------------------------------------
    patch_2D         => patch_3D%p_patch_2D(1)
    cells_in_domain => patch_2D%cells%in_domain
    field_column    => ocean_tracer%concentration
    !-----------------------------------------------------------------------
    slev   = 1
    dt_inv = 1.0_wp/dtime
    dt = dtime
    unit(1:n_zlev) = dt

    DO jb = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, jb, i_startidx_c, i_endidx_c)
      DO jc = i_startidx_c, i_endidx_c
        z_dolic = patch_3D%p_patch_1D(1)%dolic_c(jc,jb)

        IF (z_dolic > 0 ) THEN

          ! recalculate coefficients
          prism_thickness(1)     = patch_3D%p_patch_1D(1)%del_zlev_m(1) + h_c(jc,jb)
          inv_prism_thickness(1) = 1.0_wp / prism_thickness(1)
          DO jk=2,z_dolic
            prism_thickness(jk)            = patch_3D%p_patch_1D(1)%del_zlev_m(jk)
            ! inv_prism_thickness(jk)        = patch_3D%p_patch_1D(1)%inv_prism_thick_c(jc,jk,jb)
            inv_prism_thickness(jk)        = 1.0_wp / patch_3D%p_patch_1D(1)%del_zlev_m(jk)
            inv_prisms_center_distance(jk) = 2.0_wp / (prism_thickness(jk-1) + prism_thickness(jk))
            prisms_center_distance(jk)     = (prism_thickness(jk-1) + prism_thickness(jk)) * 0.5_wp
            !inv_prisms_center_distance(jk) = 1.0_wp / patch_3D%p_patch_1D(1)%del_zlev_i(jk)
          ENDDO

          !------------------------------------
             ! top level
            a(1) = 0.0_wp
            c(1) = -A_v(jc,2,jb) * inv_prism_thickness(1) * inv_prisms_center_distance(2)
            b(1) = dt_inv - c(1)
            !Fill triangular matrix
            !b is diagonal a is the lower diagonal, c is the upper
            DO jk = 2, z_dolic-1
              a(jk) = - A_v(jc,jk,jb)   * inv_prism_thickness(jk) * inv_prisms_center_distance(jk)
              c(jk) = - A_v(jc,jk+1,jb) * inv_prism_thickness(jk) * inv_prisms_center_distance(jk+1)
              b(jk) = dt_inv - a(jk) - c(jk)
            END DO
            !bottom
            a(z_dolic) = -A_v(jc,z_dolic,jb) * inv_prism_thickness(z_dolic) * inv_prisms_center_distance(z_dolic)
            c(z_dolic) = 0.0_wp
            b(z_dolic) = dt_inv - a(z_dolic)

            ! without dt_iv
             ! top level
!            a(1) = 0.0_wp
!            c(1) = - A_v(jc,2,jb) * inv_prism_thickness(1) * inv_prisms_center_distance(2)
!            b(1) = 1.0_wp - c(1)
!            !Fill triangular matrix
!            !b is diagonal a is the lower diagonal, c is the upper
!            DO jk = 2, z_dolic-1
!              a(jk) = - A_v(jc,jk,jb) * inv_prism_thickness(jk) * inv_prisms_center_distance(jk)
!              c(jk) = - A_v(jc,jk+1,jb) * inv_prism_thickness(jk) * inv_prisms_center_distance(jk+1)
!              b(jk) = 1.0_wp - a(jk) - c(jk)
!            END DO
!            !bottom
!            a(z_dolic) = -dt * A_v(jc,z_dolic,jb) * inv_prism_thickness(z_dolic) * inv_prisms_center_distance(z_dolic)
!            c(z_dolic) = 0.0_wp
!            b(z_dolic) = 1.0_wp - a(z_dolic)
!            column_tracer(1:z_dolic) = field_column(jc,1:z_dolic,jb)

!
            ! get locally the column tracer / dt
            CALL apply_triangular_matrix(z_dolic,c,b,a,unit(1:z_dolic),residual_fraction)
!            write(0,*) residual
!            stop
 !           column_tracer(1:z_dolic) = column_tracer(1:z_dolic) * residual(1:z_dolic)
!            column_tracer(1:z_dolic) = field_column(jc,1:z_dolic,jb) * dt_inv / residual(1:z_dolic)
            column_tracer(1:z_dolic) = field_column(jc,1:z_dolic,jb) * dt_inv

            column_tracer(1:z_dolic) = field_column(jc,1:z_dolic,jb) * dt_inv

            ! multiply with the diagonal (1, b1, b2, ...)
            DO jk = 2, z_dolic
              column_tracer(jk) = column_tracer(jk) * b(jk-1)
              a(jk) = a(jk) *  b(jk-1)
              b(jk) = b(jk) *  b(jk-1)
              ! c(jk) = c(jk) *  b(jk-1)
              c(jk) = dt_inv * b(jk-1) - a(jk) - b(jk)
            ENDDO
            c(z_dolic) = 0.0_wp
!

          !------------------------------------

          ! solver from lapack
          !
          ! eliminate lower diagonal
          DO jk=slev, z_dolic-1
            fact(jk+1) = a( jk+1 ) / b( jk )
            b( jk+1 ) = b( jk+1 ) - fact(jk+1) * c( jk )
            column_tracer( jk+1 ) = column_tracer( jk+1 ) - fact(jk+1) * column_tracer( jk )
 !           b( jk+1 ) = b( jk+1 ) - c( jk ) * a( jk+1 ) / b( jk )
 !           column_tracer( jk+1 ) = column_tracer( jk+1 ) - column_tracer( jk ) * a( jk+1 ) / b( jk )
          ENDDO
  !        DO jk=slev+1, z_dolic
  !          a(jk) = 0.0_wp
  !        ENDDO
          old_tracer(:) = column_tracer(:)

          !     Back solve with the matrix U from the factorization.
          column_tracer( z_dolic ) = column_tracer( z_dolic ) / b( z_dolic )
          DO jk =  z_dolic-1, 1, -1
            column_tracer( jk ) = ( column_tracer( jk ) - c( jk ) * column_tracer( jk+1 ) ) / b( jk )
          ENDDO

          ! check residual
          a(:) = 0.0_wp
          DO resIter=1, vertical_diffusion_resIterations
            CALL apply_triangular_matrix(z_dolic,c,b,a,column_tracer,residual)
            residual(1:z_dolic) = (old_tracer(1:z_dolic) - residual(1:z_dolic)) * residual_fraction(1:z_dolic)
            !     Back solve with the matrix U from the factorization.
            column_tracer( z_dolic ) =  column_tracer( z_dolic ) + residual( z_dolic ) / b( z_dolic )
            DO jk =  z_dolic-1, 1, -1
              column_tracer( jk ) = column_tracer( jk ) + ( residual( jk ) - c( jk ) * residual( jk+1 ) ) / b( jk )
            ENDDO
          ENDDO

          DO jk = 1, z_dolic
            ocean_tracer%concentration(jc,jk,jb) = column_tracer(jk)
          ENDDO

          ! check new residual
          CALL apply_triangular_matrix(z_dolic,c,b,a,column_tracer,residual)
          DO jk = 1, z_dolic
            residual_3D(jc,jk,jb) = old_tracer(jk) - residual(jk)
          ENDDO



        ENDIF ! z_dolic > 0

      END DO ! jc = i_startidx_c, i_endidx_c
    END DO ! jb = cells_in_domain%start_block, cells_in_domain%end_block

    ! CALL sync_patch_array(SYNC_C, patch_2D, diff_column)

  END SUBROUTINE tracer_diffusion_vert_implicit_r1
  !------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !!Subroutine implements implicit vertical diffusion for scalar fields.
  !>
  !! sbr identical to sbr above but now with homogeneous boundary conditions
  !!
  !! @par Revision History
  !! Developed  by  Peter Korn, MPI-M (2011).
  !! Preconditioning by Leonidas Linardakis, MPI-M (2014)
  !!
  !! The result ocean_tracer%concetration is calculated on domain_cells
  !-------------------------------------------------------------------------
  SUBROUTINE tracer_diffusion_vert_implicit_r2( patch_3D,               &
                                           & ocean_tracer, h_c, A_v,   &
                                           & operators_coefficients ) !,  &
                                          ! & diff_column)

    TYPE(t_patch_3D ),TARGET, INTENT(IN) :: patch_3D
    TYPE(t_ocean_tracer), TARGET         :: ocean_tracer
    REAL(wp), INTENT(inout)              :: h_c(:,:)  !surface height, relevant for thickness of first cell, in
    REAL(wp), INTENT(inout)              :: A_v(:,:,:)
    TYPE(t_operator_coeff),TARGET     :: operators_coefficients
    !
    !
    INTEGER :: slev
    INTEGER :: jc, jk, jb
    INTEGER :: i_startidx_c, i_endidx_c
    REAL(wp) :: a(1:n_zlev), b(1:n_zlev), c(1:n_zlev)
    REAL(wp) :: fact(1:n_zlev)
    ! REAL(wp) :: inv_prisms_center_distance(1:n_zlev)
    ! REAL(wp) :: inv_prism_thickness(1:n_zlev)
    REAL(wp) :: dt_inv
    REAL(wp), POINTER   :: field_column(:,:,:)
    REAL(wp) :: column_tracer(1:n_zlev)
    REAL(wp) :: inv_prism_thickness(1:n_zlev), inv_prisms_center_distance(1:n_zlev)
    INTEGER  :: z_dolic
    TYPE(t_subset_range), POINTER :: cells_in_domain
    TYPE(t_patch), POINTER         :: patch_2D
    REAL(wp) :: prism_thickness(1:n_zlev)
    REAL(wp) :: prisms_center_distance(1:n_zlev), unit(1:n_zlev)
    !-----------------------------------------------------------------------
    patch_2D        => patch_3D%p_patch_2D(1)
    cells_in_domain => patch_2D%cells%in_domain
    field_column    => ocean_tracer%concentration
    !-----------------------------------------------------------------------
    slev   = 1
    dt_inv = 1.0_wp/dtime

    DO jb = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, jb, i_startidx_c, i_endidx_c)
      DO jc = i_startidx_c, i_endidx_c
        z_dolic = patch_3D%p_patch_1D(1)%dolic_c(jc,jb)

        IF (z_dolic <= 0 ) CYCLE

        DO jk=1,z_dolic
          inv_prism_thickness(jk)        = patch_3D%p_patch_1D(1)%inv_prism_thick_c(jc,jk,jb)
          inv_prisms_center_distance(jk) = patch_3d%p_patch_1d(1)%inv_prism_center_dist_c(jc,jk,jb)
        ENDDO

        !------------------------------------
        ! Fill triangular matrix
        ! b is diagonal, a is the lower diagonal, c is the upper
        !   top level
        a(1) = 0.0_wp
        c(1) = -A_v(jc,2,jb) * inv_prism_thickness(1) * inv_prisms_center_distance(2)
        b(1) = dt_inv - c(1)
        DO jk = 2, z_dolic-1
          a(jk) = - A_v(jc,jk,jb)   * inv_prism_thickness(jk) * inv_prisms_center_distance(jk)
          c(jk) = - A_v(jc,jk+1,jb) * inv_prism_thickness(jk) * inv_prisms_center_distance(jk+1)
          b(jk) = dt_inv - a(jk) - c(jk)
        END DO
        ! bottom
        a(z_dolic) = -A_v(jc,z_dolic,jb) * inv_prism_thickness(z_dolic) * inv_prisms_center_distance(z_dolic)
        b(z_dolic) = dt_inv - a(z_dolic)

        ! precondition: multiply with the diagonal (1, b1, b2, ...)
        column_tracer(1) = field_column(jc,1,jb) * dt_inv
        DO jk = 2, z_dolic
          a(jk) = a(jk) *  b(jk-1)
          b(jk) = b(jk) *  b(jk-1)
          c(jk) = dt_inv * b(jk-1) - a(jk) - b(jk)

          column_tracer(jk) = field_column(jc,jk,jb) * dt_inv * b(jk-1)

        ENDDO
        c(z_dolic) = 0.0_wp

        !------------------------------------
        ! solver from lapack
        !
        ! eliminate upper diagonal
        DO jk=z_dolic-1, slev, -1
          fact(jk+1)  = c( jk ) / b( jk+1 )
          b( jk ) = b( jk ) - fact(jk+1) * a( jk +1 )
          column_tracer( jk ) = column_tracer( jk ) - fact(jk+1) * column_tracer( jk+1 )
        ENDDO

        !     Back solve with the matrix U from the factorization.
        column_tracer( 1 ) = column_tracer( 1 ) / b( 1 )
        DO jk =  2,z_dolic
          column_tracer( jk ) = ( column_tracer( jk ) - a( jk ) * column_tracer( jk-1 ) ) / b( jk )
        ENDDO

        DO jk = 1, z_dolic
          ocean_tracer%concentration(jc,jk,jb) = column_tracer(jk)
        ENDDO

      END DO ! jc = i_startidx_c, i_endidx_c
    END DO ! jb = cells_in_domain%start_block, cells_in_domain%end_block

  END SUBROUTINE tracer_diffusion_vert_implicit_r2
  !------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !!Subroutine implements implicit vertical diffusion for scalar fields.
  !>
  !! sbr identical to sbr above but now with homogeneous boundary conditions
  !!
  !! @par Revision History
  !! Developed  by  Peter Korn, MPI-M (2011).
  !! Preconditioning by Leonidas Linardakis, MPI-M (2014)
  !!
  !! The result ocean_tracer%concetration is calculated on domain_cells
  SUBROUTINE tracer_diffusion_vert_implicit_r0( patch_3D,               &
                                           & ocean_tracer, h_c, A_v,   &
                                           & operators_coefficients )

    TYPE(t_patch_3D ),TARGET, INTENT(IN) :: patch_3D
    TYPE(t_ocean_tracer), TARGET         :: ocean_tracer
    REAL(wp), INTENT(inout)              :: h_c(:,:)  !surface height, relevant for thickness of first cell, in
    REAL(wp), INTENT(inout)              :: A_v(:,:,:)
    TYPE(t_operator_coeff),TARGET     :: operators_coefficients
    !
    !
    INTEGER :: slev
    INTEGER :: jc, jk, jb
    INTEGER :: i_startidx_c, i_endidx_c
    REAL(wp) :: a(1:n_zlev), b(1:n_zlev), c(1:n_zlev)
    REAL(wp) :: fact(1:n_zlev)
    ! REAL(wp) :: inv_prisms_center_distance(1:n_zlev)
    ! REAL(wp) :: inv_prism_thickness(1:n_zlev)
    REAL(wp) :: dt_inv
    REAL(wp), POINTER   :: field_column(:,:,:)
    REAL(wp) :: column_tracer(1:n_zlev)
    REAL(wp) :: inv_prism_thickness(1:n_zlev), inv_prisms_center_distance(1:n_zlev)
    INTEGER  :: z_dolic
    TYPE(t_subset_range), POINTER :: cells_in_domain
    TYPE(t_patch), POINTER         :: patch_2D
    REAL(wp) :: prism_thickness(1:n_zlev)
    REAL(wp) :: prisms_center_distance(1:n_zlev), unit(1:n_zlev)
    !-----------------------------------------------------------------------
    patch_2D        => patch_3D%p_patch_2D(1)
    cells_in_domain => patch_2D%cells%in_domain
    field_column    => ocean_tracer%concentration
    !-----------------------------------------------------------------------
    slev   = 1
    dt_inv = 1.0_wp/dtime

    DO jb = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, jb, i_startidx_c, i_endidx_c)
      DO jc = i_startidx_c, i_endidx_c
        z_dolic = patch_3D%p_patch_1D(1)%dolic_c(jc,jb)

        IF (z_dolic <= 0 ) CYCLE

        DO jk=1,z_dolic
          inv_prism_thickness(jk)        = patch_3D%p_patch_1D(1)%inv_prism_thick_c(jc,jk,jb)
          inv_prisms_center_distance(jk) = patch_3d%p_patch_1d(1)%inv_prism_center_dist_c(jc,jk,jb)
        ENDDO

        !------------------------------------
        ! Fill triangular matrix
        ! b is diagonal, a is the lower diagonal, c is the upper
        !   top level
        a(1) = 0.0_wp
        c(1) = -A_v(jc,2,jb) * inv_prism_thickness(1) * inv_prisms_center_distance(2)
        b(1) = dt_inv - c(1)
        DO jk = 2, z_dolic-1
          a(jk) = - A_v(jc,jk,jb)   * inv_prism_thickness(jk) * inv_prisms_center_distance(jk)
          c(jk) = - A_v(jc,jk+1,jb) * inv_prism_thickness(jk) * inv_prisms_center_distance(jk+1)
          b(jk) = dt_inv - a(jk) - c(jk)
        END DO
        ! bottom
        a(z_dolic) = -A_v(jc,z_dolic,jb) * inv_prism_thickness(z_dolic) * inv_prisms_center_distance(z_dolic)
        b(z_dolic) = dt_inv - a(z_dolic)

        ! precondition: multiply with the diagonal (1, b1, b2, ...)
        column_tracer(1) = field_column(jc,1,jb) * dt_inv
        DO jk = 2, z_dolic
          a(jk) = a(jk) *  b(jk-1)
          b(jk) = b(jk) *  b(jk-1)
          c(jk) = dt_inv * b(jk-1) - a(jk) - b(jk)

          column_tracer(jk) = field_column(jc,jk,jb) * dt_inv * b(jk-1)

        ENDDO
        c(z_dolic) = 0.0_wp

        !------------------------------------
        ! solver from lapack
        !
        ! eliminate lower diagonal
        DO jk=slev, z_dolic-1
          fact(jk+1) = a( jk+1 ) / b( jk )
          b( jk+1 ) = b( jk+1 ) - fact(jk+1) * c( jk )
          column_tracer( jk+1 ) = column_tracer( jk+1 ) - fact(jk+1) * column_tracer( jk )
        ENDDO

        !     Back solve with the matrix U from the factorization.
        column_tracer( z_dolic ) = column_tracer( z_dolic ) / b( z_dolic )
        DO jk =  z_dolic-1, 1, -1
          column_tracer( jk ) = ( column_tracer( jk ) - c( jk ) * column_tracer( jk+1 ) ) / b( jk )
        ENDDO

        DO jk = 1, z_dolic
          ocean_tracer%concentration(jc,jk,jb) = column_tracer(jk)
        ENDDO

      END DO ! jc = i_startidx_c, i_endidx_c
    END DO ! jb = cells_in_domain%start_block, cells_in_domain%end_block

  END SUBROUTINE tracer_diffusion_vert_implicit_r0
  !------------------------------------------------------------------------

  
END MODULE mo_ocean_testbed_modules
