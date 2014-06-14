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
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
!!
MODULE mo_ocean_testbed_modules
  !-------------------------------------------------------------------------
  USE mo_kind,                   ONLY: wp
  USE mo_impl_constants,         ONLY: max_char_length
  USE mo_model_domain,           ONLY: t_patch, t_patch_3d
  USE mo_grid_config,            ONLY: n_dom
  USE mo_ocean_nml,              ONLY: n_zlev
  USE mo_dynamics_config,        ONLY: nold
  USE mo_run_config,             ONLY: nsteps, dtime, output_mode, test_mode
  USE mo_exception,              ONLY: message, message_text, finish
  !USE mo_io_units,               ONLY: filename_max
  USE mo_datetime,               ONLY: t_datetime, add_time, datetime_to_string
  USE mo_timer,                  ONLY: timer_start, timer_stop, timer_total
  USE mo_oce_ab_timestepping,    ONLY: &
!    solve_free_surface_eq_ab, &
!   & calc_vert_velocity,       &
    & update_time_indices
  USE mo_oce_types,              ONLY: t_hydro_ocean_state
  USE mo_operator_ocean_coeff_3d,ONLY: t_operator_coeff! , update_diffusion_matrices
  USE mo_oce_tracer,             ONLY: advect_tracer_ab
  USE mo_oce_bulk,               ONLY: update_surface_flux
  USE mo_sea_ice,                ONLY: ice_budgets
  USE mo_sea_ice_types,          ONLY: t_sfc_flx, t_atmos_fluxes, t_atmos_for_ocean, &
    & t_sea_ice
  USE mo_physical_constants,     ONLY: rhoi, rhos, rho_ref
  USE mo_oce_physics,            ONLY: t_ho_params
  USE mo_oce_thermodyn,          ONLY: calc_neutralslope_coeff, calc_neutralslope_coeff_func
  USE mo_time_config,            ONLY: time_config
  USE mo_statistics
  USE mo_util_dbg_prnt,          ONLY: dbg_print
  USE mo_ocean_statistics
  USE mo_ocean_output
  USE mo_parallel_config,        ONLY: nproma
  USE mo_statistics
  USE mo_ocean_testbed_vertical_diffusion
  USE mo_hydro_ocean_run,        ONLY: write_initial_ocean_timestep

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ocean_test_modules
  
  CHARACTER(len=12)           :: debug_string = 'testbed     '  ! Output of module for 1 line debug
  
  !-------------------------------------------------------------------------
CONTAINS

  !-------------------------------------------------------------------------
  !>
  SUBROUTINE ocean_test_modules( patch_3d, ocean_state,  &
    & datetime, surface_fluxes, physics_parameters,             &
    & oceans_atmosphere, oceans_atmosphere_fluxes, ocean_ice,operators_coefficients)

    TYPE(t_patch_3d ),TARGET, INTENT(inout)          :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET, INTENT(inout) :: ocean_state(n_dom)
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
        CALL ocean_test_advection( patch_3d, ocean_state, &
          & datetime, surface_fluxes, physics_parameters,             &
          & ocean_ice,operators_coefficients)

      CASE (2)
        CALL test_tracer_diffusion_vertical_implicit( patch_3d, ocean_state, physics_parameters,  &
           & operators_coefficients)

      CASE (3)
        CALL test_velocity_diffusion_vert_implicit( patch_3d, ocean_state, physics_parameters,  &
           & operators_coefficients)

      CASE (4)
        CALL test_surface_flux( patch_3d, ocean_state,  &
          & datetime, surface_fluxes,             &
          & oceans_atmosphere, oceans_atmosphere_fluxes, ocean_ice, operators_coefficients)

      CASE (5)
        CALL test_neutralcoeff( patch_3d, ocean_state)

      CASE DEFAULT
        CALL finish(method_name, "Unknown test_mode")

    END SELECT



  END SUBROUTINE ocean_test_modules
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !>
  SUBROUTINE ocean_test_advection( patch_3d, ocean_state, &
    & datetime, surface_fluxes, physics_parameters,             &
    & ocean_ice,operators_coefficients)
    
    TYPE(t_patch_3d ),TARGET, INTENT(inout)          :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET, INTENT(inout) :: ocean_state(n_dom)
    TYPE(t_datetime), INTENT(inout)                  :: datetime
    TYPE(t_sfc_flx)                                  :: surface_fluxes
    TYPE (t_ho_params)                               :: physics_parameters
    TYPE (t_sea_ice),         INTENT(inout)          :: ocean_ice
    TYPE(t_operator_coeff),   INTENT(inout)          :: operators_coefficients
    
    ! local variables
    INTEGER :: jstep, jg
    !LOGICAL                         :: l_outputtime
    CHARACTER(LEN=32)               :: datestring
    TYPE(t_patch), POINTER :: patch_2d
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
          
        CALL output_ocean( patch_3d, &
          & ocean_state,             &
          & datetime,                &
          & surface_fluxes,          &
          & ocean_ice,               &
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
  !>
  SUBROUTINE test_surface_flux( patch_3d, p_os, &
    & datetime, surface_fluxes,         &
    & p_as, atmos_fluxes, p_ice, operators_coefficients)
    
    TYPE(t_patch_3d ),TARGET, INTENT(inout)          :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET, INTENT(inout) :: p_os(n_dom)
    TYPE(t_datetime), INTENT(inout)                  :: datetime
    TYPE(t_sfc_flx)                                  :: surface_fluxes
    TYPE(t_atmos_for_ocean),  INTENT(inout)          :: p_as
    TYPE(t_atmos_fluxes ),    INTENT(inout)          :: atmos_fluxes
    TYPE(t_sea_ice),          INTENT(inout)          :: p_ice
    TYPE(t_operator_coeff),   INTENT(inout)          :: operators_coefficients
    
    ! local variables
    REAL(wp), DIMENSION(nproma,patch_3D%p_patch_2D(1)%alloc_cell_blocks) &
      &                                              :: draft, &
      &                                                 saltBefore, saltAfter, saltBudget
    INTEGER :: jstep
    !INTEGER :: ocean_statistics
    !LOGICAL                         :: l_outputtime
    CHARACTER(LEN=32)               :: datestring
    TYPE(t_patch), POINTER :: patch_2d
    INTEGER :: jstep0
    
    CHARACTER(LEN=max_char_length), PARAMETER :: &
      & method_name = 'mo_ocean_testbed_modules:test_sea_ice'
    !------------------------------------------------------------------
    patch_2D      => patch_3d%p_patch_2d(1)

    CALL datetime_to_string(datestring, datetime)

    time_config%sim_time(:) = 0.0_wp
    jstep0 = 0
    !------------------------------------------------------------------

    !  Overwrite init:
    WHERE (p_ice%hi(:,1,:) > 1.0_wp)
      p_ice%hi(:,1,:) = 0.3_wp
      p_ice%hs(:,1,:) = 0.1_wp
      p_ice%conc(:,1,:) = 0.5_wp
      p_os(n_dom)%p_prog(nold(1))%tracer(:,1,:,1) = -1.7_wp
    ENDWHERE

    ! TODO: use prism_thick_flat_sfc_c instead of del_zlev_m
    draft(:,:)           = (rhos * p_ice%hs(:,1,:) + rhoi * p_ice%hi(:,1,:)) / rho_ref
    p_ice%zUnderIce(:,:) = patch_3d%p_patch_1d(1)%prism_thick_flat_sfc_c(:,1,:) +  p_os(n_dom)%p_prog(nold(1))%h(:,:) &
      &                  - draft(:,:) * p_ice%conc(:,1,:)

    CALL dbg_print('sfcflx: draft    ' ,draft          ,debug_string, 4, in_subset=patch_3d%p_patch_2D(1)%cells%owned)
    CALL dbg_print('sfcflx: zUnderIce' ,p_ice%zUnderIce,debug_string, 4, in_subset=patch_3d%p_patch_2D(1)%cells%owned)

    !------------------------------------------------------------------
    ! write initial
    !------------------------------------------------------------------
    IF (output_mode%l_nml) THEN
      CALL write_initial_ocean_timestep(patch_3D,p_os(n_dom),surface_fluxes,p_ice,jstep0)
    ENDIF
    DO jstep = (jstep0+1), (jstep0+nsteps)
    
      p_os(n_dom)%p_prog(nold(1))%h(:,:) = 0.0_wp  !  do not change h
      CALL datetime_to_string(datestring, datetime)
      WRITE(message_text,'(a,i10,2a)') '  Begin of timestep =',jstep,'  datetime:  ', datestring
      CALL message (TRIM(method_name), message_text)

      ! Set model time.
      CALL add_time(dtime,0,0,0,datetime)

      ! diagnostics BEFORE
      saltBefore = ice_budgets(patch_2D, p_ice, p_os(n_dom),'BEFORE')

      ! call component
      CALL update_surface_flux(patch_3D, p_os(n_dom), p_as, p_ice, atmos_fluxes, surface_fluxes, jstep, datetime, &
        &  operators_coefficients)

      ! diagnostics AFTER
      saltAfter =  ice_budgets(patch_2D, p_ice, p_os(n_dom),'AFTER')

      ! compute budget
      saltBudget = saltAfter - saltBefore

      CALL dbg_print('IceBudget: salt diff', &
      &              saltBudget , &
      &              debug_string, 4, in_subset=patch_3d%p_patch_2D(1)%cells%owned)

      time_config%sim_time(1) = time_config%sim_time(1) + dtime

      ! add values to output field
      p_ice%budgets%salt_00(:,:) = saltBudget(:,:)

      CALL dbg_print('IceBudget: salt_00' ,p_ice%budgets%salt_00 ,debug_string, 4, in_subset=patch_3d%p_patch_2D(1)%cells%owned)

      CALL output_ocean( patch_3D,   &
        &                p_os(n_dom),&
        &                datetime,   &
        &                surface_fluxes,  &
        &                p_ice,      &
        &                jstep, jstep0)

    END DO
    
  END SUBROUTINE test_surface_flux
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  SUBROUTINE test_neutralcoeff( patch_3d, p_os)
    CHARACTER(LEN=*), PARAMETER ::  routine = "testbed: neutralcoeff"
    
    TYPE(t_patch_3d ),TARGET, INTENT(inout)          :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET, INTENT(inout) :: p_os(n_dom)

    ! local variables
    REAL(wp):: t, s, p, co(2), aob
    REAL(wp):: alph(1:nproma,1:n_zlev,1:patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    REAL(wp):: beta(1:nproma,1:n_zlev,1:patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    !INTEGER :: jk

    alph(:,:,:) = 0.0_wp
    beta(:,:,:) = 0.0_wp

    CALL calc_neutralslope_coeff( &
      &    patch_3d,              &
      &    p_os(n_dom)%p_prog(nold(1))%tracer(:,:,:,:), &
      &    p_os(n_dom)%p_prog(nold(1))%h(:,:), &
      &    alph, beta)

    !  test values
    t = 10.0_wp
    s = 40.0_wp
    p = 4000.0_wp    !  4000 dbar = 400 bar
    co = calc_neutralslope_coeff_func(t,s,p)
    aob = co(1)/co(2)

    WRITE(message_text,'(3(a,1pg18.8))') '  Parameter: alpha = ',co(1), ' beta = ',co(2), ' alpha/beta = ',aob
    CALL message (TRIM(routine), message_text)

  END SUBROUTINE test_neutralcoeff
  !-------------------------------------------------------------------------
  
END MODULE mo_ocean_testbed_modules
