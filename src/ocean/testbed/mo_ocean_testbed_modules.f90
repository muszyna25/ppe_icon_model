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
  USE mo_model_domain,           ONLY: t_patch, t_patch_3d,t_subset_range
  USE mo_grid_config,            ONLY: n_dom
  USE mo_ocean_nml,              ONLY: n_zlev, GMRedi_configuration, GMRedi_combined, Cartesian_Mixing, &
    & atmos_flux_analytical_type, no_tracer, surface_module, OceanReferenceDensity
  USE mo_sea_ice_nml,            ONLY: init_analytic_conc_param, t_heat_base
  USE mo_dynamics_config,        ONLY: nold, nnew
  USE mo_run_config,             ONLY: nsteps, dtime, output_mode, test_mode !, test_param
  USE mo_exception,              ONLY: message, message_text, finish
  USE mo_ext_data_types,         ONLY: t_external_data
  !USE mo_io_units,               ONLY: filename_max
  USE mo_timer,                  ONLY: timer_start, timer_stop, timer_total
  USE mo_ocean_ab_timestepping,  ONLY: &
!    solve_free_surface_eq_ab, &
!   & calc_vert_velocity,       &
    & update_time_indices
  USE mo_random_util,            ONLY: add_random_noise_global
  USE mo_ocean_types,            ONLY: t_hydro_ocean_state, t_operator_coeff, t_solvercoeff_singleprecision
  USE mo_hamocc_types,           ONLY: t_hamocc_state
  USE mo_restart,                ONLY: t_RestartDescriptor, createRestartDescriptor, deleteRestartDescriptor
  USE mo_restart_attributes,     ONLY: t_RestartAttributeList, getAttributesForRestarting
  USE mo_io_config,              ONLY: n_checkpoints, write_last_restart
  USE mo_operator_ocean_coeff_3d,ONLY: t_operator_coeff! , update_diffusion_matrices
  USE mo_ocean_tracer,           ONLY: advect_ocean_tracers
  USE mo_ocean_bulk,             ONLY: update_surface_flux
  USE mo_ocean_surface,          ONLY: update_ocean_surface
  USE mo_ocean_surface_types,    ONLY: t_ocean_surface
  USE mo_sea_ice,                ONLY: salt_content_in_surface, energy_content_in_surface
  USE mo_sea_ice_types,          ONLY: t_sfc_flx, t_atmos_fluxes, t_atmos_for_ocean, &
    & t_sea_ice
  USE mo_physical_constants,     ONLY: rhoi, rhos, clw, alf, Tf
  USE mo_ocean_physics,          ONLY: t_ho_params
  USE mo_master_config,          ONLY: isRestart
  USE mo_ocean_GM_Redi,          ONLY: calc_neutralslope_coeff, calc_neutralslope_coeff_func_onColumn, &
  &                                    prepare_ocean_physics,calc_ocean_physics
  USE mo_ocean_diagnostics,      ONLY: calc_fast_oce_diagnostics, calc_psi
  USE mo_ocean_thermodyn,        ONLY: calc_potential_density, calculate_density
  USE mo_time_config,            ONLY: time_config
  USE mo_statistics
  USE mo_util_dbg_prnt,          ONLY: dbg_print
  USE mo_ocean_statistics
  USE mo_ocean_output
  USE mo_parallel_config,        ONLY: nproma
  USE mo_statistics
  USE mo_ocean_testbed_vertical_diffusion
  USE mo_ocean_math_operators
  USE mo_ocean_ab_timestepping
  USE mo_grid_subset,            ONLY: t_subset_range, get_index_range 
  USE mo_ocean_diffusion,        ONLY: tracer_diffusion_vertical_implicit,tracer_diffusion_horz
  USE mo_scalar_product,         ONLY: calc_scalar_product_veloc_3d
  USE mo_ocean_tracer_transport_horz, ONLY: diffuse_horz
  USE mo_hydro_ocean_run
  USE mo_ocean_physics

  USE mtime,                     ONLY: datetime, newDatetime, deallocateDatetime, datetimeToString, &
       &                               timedelta, newTimedelta, deallocateTimedelta,                &
       &                               MAX_DATETIME_STR_LEN, newDatetime,                           &
       &                               MAX_MTIME_ERROR_STR_LEN, no_error, mtime_strerror,           &
       &                               OPERATOR(-), OPERATOR(+), OPERATOR(>), OPERATOR(*),          &
       &                               ASSIGNMENT(=), OPERATOR(==), OPERATOR(>=), OPERATOR(/=),     &
       &                               event, eventGroup, newEvent,                                 &
       &                               addEventToEventGroup, isCurrentEventActive
  USE mo_event_manager,          ONLY: initEventManager, addEventGroup, getEventGroup, printEventGroup

  USE mo_hamocc_types,          ONLY: t_hamocc_state
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ocean_test_modules
  
  CHARACTER(len=12)           :: debug_string = 'testbedMod  '  ! Output of module for 1 line debug
  
  !-------------------------------------------------------------------------
CONTAINS

  !-------------------------------------------------------------------------
  !>
  SUBROUTINE ocean_test_modules( patch_3d, ocean_state,  external_data,  &
    & this_datetime, surface_fluxes, ocean_surface, physics_parameters,             &
    & oceans_atmosphere, oceans_atmosphere_fluxes, ocean_ice, operators_coefficients, &
    & solvercoeff_sp)

    TYPE(t_patch_3d ),TARGET, INTENT(inout)          :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET, INTENT(inout) :: ocean_state(n_dom)
    TYPE(t_external_data), TARGET, INTENT(in)        :: external_data(n_dom)
    TYPE(datetime), POINTER                          :: this_datetime
    TYPE(t_sfc_flx)                                  :: surface_fluxes
    TYPE(t_ocean_surface)                            :: ocean_surface
    TYPE(t_ho_params)                                :: physics_parameters
    TYPE(t_atmos_for_ocean),  INTENT(inout)          :: oceans_atmosphere
    TYPE(t_atmos_fluxes ),    INTENT(inout)          :: oceans_atmosphere_fluxes
    TYPE(t_sea_ice),          INTENT(inout)          :: ocean_ice
    TYPE(t_operator_coeff),   INTENT(inout)          :: operators_coefficients
    TYPE(t_solvercoeff_singleprecision), INTENT(inout) :: solvercoeff_sp

    CHARACTER(LEN=*), PARAMETER ::  method_name = "ocean_test_modules"
    TYPE (t_hamocc_state)        :: hamocc_State

    SELECT CASE (test_mode)  !  1 - 99 test ocean modules
      CASE (1)
        CALL ocean_test_advection( patch_3d, ocean_state, &
          & this_datetime, surface_fluxes, physics_parameters,             &
          & ocean_ice,operators_coefficients)

      CASE (2)
        CALL test_tracer_diffusion_vertical_implicit( patch_3d, ocean_state, physics_parameters,  &
           & operators_coefficients)

      CASE (3)
        CALL test_velocity_diffusion_vert_implicit( patch_3d, ocean_state, physics_parameters,  &
           & operators_coefficients)

      CASE (4)
        CALL test_surface_flux( patch_3d, ocean_state,  &
          & this_datetime, surface_fluxes,             &
          & oceans_atmosphere, oceans_atmosphere_fluxes, ocean_ice, operators_coefficients)

      CASE (41)
        CALL test_surface_flux_slo( patch_3d, ocean_state,  &
          & this_datetime, surface_fluxes, ocean_surface,        &
          & oceans_atmosphere, oceans_atmosphere_fluxes, ocean_ice, operators_coefficients)

      CASE (5)
        CALL test_neutralcoeff( patch_3d, ocean_state)


      CASE (10)
        CALL ocean_test_GMRedi( patch_3d, ocean_state, &
          & this_datetime, surface_fluxes, physics_parameters,             &
          & ocean_ice,operators_coefficients)

      CASE (11) ! surface only processing to get quasi output fast
        CALL test_output( patch_3d, ocean_state,  &
          & this_datetime, surface_fluxes,             &
          & physics_parameters,                   &
          & oceans_atmosphere, oceans_atmosphere_fluxes, ocean_ice, hamocc_state,operators_coefficients)

      CASE (12) ! surface only processing to get quasi output fast
        CALL test_events( patch_3d, ocean_state,  &
          & external_data , &
          & this_datetime, &
          & surface_fluxes, &
          & ocean_surface,   &
          & physics_parameters, &
          & oceans_atmosphere, &
          & oceans_atmosphere_fluxes, &
          & ocean_ice, &
          & hamocc_state, &
          & operators_coefficients, &
          & solvercoeff_sp)

      CASE DEFAULT
        CALL finish(method_name, "Unknown test_mode")

    END SELECT



  END SUBROUTINE ocean_test_modules
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !>
  SUBROUTINE ocean_test_GMRedi( patch_3d, ocean_state, &
    & this_datetime, surface_fluxes, physics_parameters,             &
    & ocean_ice,operators_coefficients)
    
    TYPE(t_patch_3d ),TARGET, INTENT(inout)          :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET, INTENT(inout) :: ocean_state(n_dom)
    TYPE(datetime), POINTER                          :: this_datetime
    TYPE(t_sfc_flx)                                  :: surface_fluxes
    TYPE (t_ho_params)                               :: physics_parameters
    TYPE (t_sea_ice),         INTENT(inout)          :: ocean_ice
    TYPE(t_operator_coeff),   INTENT(inout)          :: operators_coefficients
    
    ! local variables
    TYPE (t_hamocc_state)        :: hamocc_State
    INTEGER :: jstep, jg
    !LOGICAL                         :: l_outputtime
    CHARACTER(LEN=32)               :: datestring
    TYPE(t_patch), POINTER :: patch_2d
    INTEGER :: jstep0 ! start counter for time loop
    INTEGER :: tracer_index
    
    INTEGER :: jc,level,jb
    INTEGER :: start_cell_index, end_cell_index
    TYPE(t_subset_range), POINTER :: cells_in_domain, edges_in_domain
    REAL(wp) :: delta_t
    
    REAL(wp) :: z_diff_flux_h(nproma, n_zlev, patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp) :: div_diff_flux_horz(nproma,n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp) :: div_diff_flx_vert(nproma, n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp) :: trac_cart(nproma,n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp) :: div_diff_flux_horz_cart(nproma,n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks)

    TYPE(timedelta), POINTER :: model_time_step => NULL()

    CHARACTER(LEN=max_char_length), PARAMETER :: &
      & method_name = 'mo_ocean_testbed_modules:ocean_test_advection'
    !------------------------------------------------------------------
    tracer_index=2!test is on salinity

    
    
    patch_2D      => patch_3d%p_patch_2d(1)
   
    cells_in_domain => patch_2D%cells%in_domain
    edges_in_domain => patch_2D%edges%in_domain
    delta_t = dtime

    z_diff_flux_h(1:nproma,1:n_zlev,1:patch_3d%p_patch_2d(1)%nblks_e)=0.0_wp
    div_diff_flux_horz(1:nproma,1:n_zlev,1:patch_2d%alloc_cell_blocks)=0.0_wp
    div_diff_flx_vert (1:nproma,1:n_zlev,1:patch_2d%alloc_cell_blocks)=0.0_wp 
    trac_cart(1:nproma,1:n_zlev, 1:patch_3d%p_patch_2d(1)%alloc_cell_blocks)=0.0_wp
    div_diff_flux_horz_cart(1:nproma,1:n_zlev, 1:patch_3d%p_patch_2d(1)%alloc_cell_blocks)=0.0_wp
    !---------------------------------------------------------------------   
    CALL datetimeToString(this_datetime, datestring)

    ! IF (ltimer) CALL timer_start(timer_total)
    CALL timer_start(timer_total)
        
    jstep0 = 0
    jg=1
    !------------------------------------------------------------------
    ! IF(.NOT.l_time_marching)THEN

      !IF(itestcase_oce==28)THEN
      DO jstep = (jstep0+1), (jstep0+nsteps)
      
        CALL datetimeToString(this_datetime, datestring)
        WRITE(message_text,'(a,i10,2a)') '  Begin of timestep =',jstep,'  datetime:  ', datestring
        CALL message (TRIM(method_name), message_text)
 
!          IF(jstep==1)THEN
!          ocean_state(jg)%p_diag%vn_time_weighted = ocean_state(jg)%p_prog(nold(1))%vn
!          ocean_state(jg)%p_prog(nnew(1))%vn = ocean_state(jg)%p_prog(nold(1))%vn
!          ocean_state(jg)%p_diag%w        =  0.0_wp!0.0833_wp!0.025_wp
!          ocean_state(jg)%p_diag%w(:,:,:) = -0.0833_wp!0.025_wp
!          ENDIF  
       CALL calc_scalar_product_veloc_3d( patch_3d,  &
        & ocean_state(n_dom)%p_prog(nold(1))%vn,     &
        & ocean_state(n_dom)%p_diag,                 &
        & operators_coefficients)
        
        IF(GMRedi_configuration/=Cartesian_Mixing)THEN    
          CALL prepare_ocean_physics(patch_3d, &
            & ocean_state(n_dom),    &
            & physics_parameters, &
            & operators_coefficients)
        ENDIF
        
!        DO tracer_index=1,2
        
          IF(GMRedi_configuration/=Cartesian_Mixing)THEN
        
            CALL calc_ocean_physics( patch_3d, &
                                   & ocean_state(n_dom),     &
                                   &  physics_parameters,    &
                                   &  operators_coefficients,&
                                   &  tracer_index)
            CALL div_oce_3d( ocean_state(n_dom)%p_diag%GMRedi_flux_horz(:,:,:,tracer_index),&
                     &   patch_3D, &
                     &   operators_coefficients%div_coeff, &
                     &   div_diff_flux_horz )
            !vertical div of GMRedi-flux
            CALL verticalDiv_scalar_onFullLevels( patch_3d, &
                                            & ocean_state(n_dom)%p_diag%GMRedi_flux_vert(:,:,:,tracer_index), &
                                            & div_diff_flx_vert)
                                   
         ELSE
          CALL tracer_diffusion_horz(patch_3D,&
                                   & ocean_state(n_dom)%p_prog(nold(1))%ocean_tracers(tracer_index)%concentration,&
                                   & ocean_state(n_dom), z_diff_flux_h, physics_parameters%k_tracer_h(:,:,:,tracer_index ))

            CALL div_oce_3d( z_diff_flux_h,&
                     &   patch_3D, &
                     &   operators_coefficients%div_coeff, &
                     &   div_diff_flux_horz_cart )
    
         ENDIF

         !cart
         DO jb = cells_in_domain%start_block, cells_in_domain%end_block
           CALL get_index_range(cells_in_domain, jb, start_cell_index, end_cell_index)
           DO jc = start_cell_index, end_cell_index

              DO level = 1, patch_3d%p_patch_1d(1)%dolic_c(jc,jb)
              ocean_state(n_dom)%p_prog(nnew(1))%ocean_tracers(tracer_index-1)%concentration(jc,level,jb) &
                & = ocean_state(n_dom)%p_prog(nold(1))%ocean_tracers(tracer_index-1)%concentration(jc,level,jb)-  &
                &  (delta_t /  patch_3D%p_patch_1D(1)%prism_thick_c(jc,level,jb))  &
                &    * ( - (div_diff_flux_horz_cart(jc,level,jb)))
! write(123,*)'details',level,  ocean_state(n_dom)%p_prog(nnew(1))%ocean_tracers(tracer_index)%concentration(jc,level,jb),&
! & ocean_state(n_dom)%p_prog(nold(1))%ocean_tracers(tracer_index)%concentration(jc,level,jb),&
! & div_diff_flux_horz(jc,level,jb)

           ENDDO
         END DO
       END DO


       !GM                               
       DO jb = cells_in_domain%start_block, cells_in_domain%end_block
         CALL get_index_range(cells_in_domain, jb, start_cell_index, end_cell_index)
         DO jc = start_cell_index, end_cell_index

           DO level = 1, patch_3d%p_patch_1d(1)%dolic_c(jc,jb)
             ocean_state(n_dom)%p_prog(nnew(1))%ocean_tracers(tracer_index)%concentration(jc,level,jb) &
                & = ocean_state(n_dom)%p_prog(nold(1))%ocean_tracers(tracer_index)%concentration(jc,level,jb)-  &
                &  (delta_t /  patch_3D%p_patch_1D(1)%prism_thick_c(jc,level,jb))  &
                &    * ( - (div_diff_flux_horz(jc,level,jb)-div_diff_flx_vert(jc,level,jb)))
!  write(123,*)'details',level,  ocean_state(n_dom)%p_prog(nnew(1))%ocean_tracers(tracer_index)%concentration(jc,level,jb),&
!  & ocean_state(n_dom)%p_prog(nold(1))%ocean_tracers(tracer_index)%concentration(jc,level,jb),&
!  & div_diff_flux_horz(jc,level,jb),div_diff_flx_vert(jc,level,jb),&
! &ocean_state(n_dom)%p_aux%slopes_squared(jc,level,jb),&
! &ocean_state(n_dom)%p_aux%taper_function_1(jc,level,jb),&
! &ocean_state(n_dom)%p_aux%taper_function_2(jc,level,jb)!,&

           ENDDO
         END DO
       END DO


      !cart    
      CALL tracer_diffusion_vertical_implicit(                         &
      & patch_3d,                                                      &
      & ocean_state(n_dom)%p_prog(nnew(1))%ocean_tracers(tracer_index-1),&
      & physics_parameters%a_tracer_v(:,:,:, tracer_index),            &
      & operators_coefficients)

          
      !GM    
      CALL tracer_diffusion_vertical_implicit(                         &
      & patch_3d,                                                      &
      & ocean_state(n_dom)%p_prog(nnew(1))%ocean_tracers(tracer_index),&
      & physics_parameters%a_tracer_v(:,:,:, tracer_index),            &
      & operators_coefficients)
      
      
      ocean_state(n_dom)%p_diag%rho=ocean_state(n_dom)%p_prog(nnew(1))%ocean_tracers(tracer_index-1)%concentration&
      &-ocean_state(n_dom)%p_prog(nnew(1))%ocean_tracers(tracer_index)%concentration
         
IF(tracer_index==2)THEN
DO level = 1, 5
write(0,*)'tracer:GM:Cart',&
& maxval( ocean_state(n_dom)%p_prog(nnew(1))%ocean_tracers(tracer_index)%concentration(:,level,:)),&
& minval( ocean_state(n_dom)%p_prog(nnew(1))%ocean_tracers(tracer_index)%concentration(:,level,:)),&
& maxval( ocean_state(n_dom)%p_prog(nnew(1))%ocean_tracers(tracer_index-1)%concentration(:,level,:)),&
& minval( ocean_state(n_dom)%p_prog(nnew(1))%ocean_tracers(tracer_index-1)%concentration(:,level,:)),&
& maxval( ocean_state(n_dom)%p_diag%rho(:,level,:)),&
& minval( ocean_state(n_dom)%p_diag%rho(:,level,:)),&
& maxval( div_diff_flux_horz(:,level,:)),&
& minval( div_diff_flux_horz(:,level,:))
END DO
ENDIF
      !END DO       

        ! One integration cycle finished on the lowest grid level (coarsest
        ! resolution). Set model time.

        model_time_step => newTimedelta('+', 0, 0, 0, 0, 0, NINT(dtime), 0)
        this_datetime = this_datetime + model_time_step
        CALL deallocateTimedelta(model_time_step)        
!        CALL add_time(dtime,0,0,0,this_datetime)
     
        ! update accumulated vars
        CALL update_ocean_statistics(ocean_state(1),     &
        & surface_fluxes,                                &
        & patch_2D%cells%owned,       &
        & patch_2D%edges%owned,       &
        & patch_2D%verts%owned,       &
        & n_zlev)
          
        CALL output_ocean( patch_3d, &
          & ocean_state,             &
          & this_datetime,                &
          & surface_fluxes,          &
          & ocean_ice,               &
          & hamocc_state,               &
          & jstep, jstep0)

        ! Shift time indices for the next loop
        ! this HAS to ge into the restart files, because the start with the following loop
        CALL update_time_indices(jg)
        ! update intermediate timestepping variables for the tracers
        ! velocity

      END DO
    ! ENDIF!(l_no_time_marching)THEN
    
    CALL timer_stop(timer_total)
    
  END SUBROUTINE ocean_test_GMRedi
  !-------------------------------------------------------------------------






  !-------------------------------------------------------------------------
  !>
  SUBROUTINE ocean_test_advection( patch_3d, ocean_state, &
    & this_datetime, surface_fluxes, physics_parameters,             &
    & ocean_ice,operators_coefficients)
    
    TYPE(t_patch_3d ),TARGET, INTENT(inout)          :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET, INTENT(inout) :: ocean_state(n_dom)
    TYPE(datetime), POINTER                          :: this_datetime
    TYPE(t_sfc_flx)                                  :: surface_fluxes
    TYPE (t_ho_params)                               :: physics_parameters
    TYPE (t_sea_ice),         INTENT(inout)          :: ocean_ice
    TYPE(t_operator_coeff),   INTENT(inout)          :: operators_coefficients
    
    ! local variables
    TYPE (t_hamocc_state)        :: hamocc_State
    INTEGER :: jstep, jg
    !LOGICAL                         :: l_outputtime
    CHARACTER(LEN=32)               :: datestring
    TYPE(t_patch), POINTER :: patch_2d
    INTEGER :: jstep0 ! start counter for time loop

    TYPE(timedelta), POINTER :: model_time_step => NULL()
    
    !CHARACTER(LEN=filename_max)  :: outputfile, gridfile
    CHARACTER(LEN=max_char_length), PARAMETER :: &
      & method_name = 'mo_ocean_testbed_modules:ocean_test_advection'
    !------------------------------------------------------------------
    
    patch_2D      => patch_3d%p_patch_2d(1)
    CALL datetimeToString(this_datetime, datestring)

    ! IF (ltimer) CALL timer_start(timer_total)
    CALL timer_start(timer_total)

    jstep0 = 0
    !------------------------------------------------------------------
    ! IF(.NOT.l_time_marching)THEN

      !IF(itestcase_oce==28)THEN
      DO jstep = (jstep0+1), (jstep0+nsteps)
      
        CALL datetimeToString(this_datetime, datestring)
        WRITE(message_text,'(a,i10,2a)') '  Begin of timestep =',jstep,'  datetime:  ', datestring
        CALL message (TRIM(method_name), message_text)
 
!          IF(jstep==1)THEN
!          ocean_state(jg)%p_diag%vn_time_weighted = ocean_state(jg)%p_prog(nold(1))%vn
!          ocean_state(jg)%p_prog(nnew(1))%vn = ocean_state(jg)%p_prog(nold(1))%vn
!          ocean_state(jg)%p_diag%w        =  0.0_wp!0.0833_wp!0.025_wp
!          ocean_state(jg)%p_diag%w(:,:,:) = -0.0833_wp!0.025_wp
!          ENDIF

        !CALL calc_vert_velocity(patch_3d, ocean_state(jg),operators_coefficients)
        CALL advect_ocean_tracers( patch_3d, ocean_state(jg),  &
          & physics_parameters,surface_fluxes,&
          & operators_coefficients,&
          & jstep)
        ! One integration cycle finished on the lowest grid level (coarsest
        ! resolution). Set model time.
        model_time_step => newTimedelta('+', 0, 0, 0, 0, 0, NINT(dtime), 0)
        this_datetime = this_datetime + model_time_step
        CALL deallocateTimedelta(model_time_step) 
!        CALL add_time(dtime,0,0,0,this_datetime)
      
        ! update accumulated vars
        CALL update_ocean_statistics(ocean_state(1),&
        & surface_fluxes,                                &
        & patch_2D%cells%owned,       &
        & patch_2D%edges%owned,       &
        & patch_2D%verts%owned,       &
        & n_zlev)
          
        CALL output_ocean( patch_3d, &
          & ocean_state,             &
          & this_datetime,                &
          & surface_fluxes,          &
          & ocean_ice,               &
          & hamocc_state,               &
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
  !> surface only call, regular output, restart, checkpoints
  SUBROUTINE test_output( patch_3d, p_os,           &
    & this_datetime, surface_fluxes, physics_parameters, &
    & p_as, atmos_fluxes, p_ice, hamocc_state,operators_coefficients)
    
    TYPE(t_patch_3d ),TARGET, INTENT(inout)          :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET, INTENT(inout) :: p_os(n_dom)
    TYPE(datetime), POINTER                          :: this_datetime
    TYPE(t_sfc_flx)                                  :: surface_fluxes
    TYPE (t_ho_params)                               :: physics_parameters
    TYPE(t_atmos_for_ocean),  INTENT(inout)          :: p_as
    TYPE(t_atmos_fluxes ),    INTENT(inout)          :: atmos_fluxes
    TYPE(t_sea_ice),          INTENT(inout)          :: p_ice
    TYPE(t_hamocc_state),          INTENT(inout)      ::hamocc_state
    TYPE(t_operator_coeff),   INTENT(inout)          :: operators_coefficients
    
    ! local variables
    INTEGER                       :: jstep
    INTEGER                       :: block, cell, cellStart,cellEnd
    TYPE(t_subset_range), POINTER :: subset
    !INTEGER                      :: ocean_statistics
    !LOGICAL                      :: l_outputtime
    CHARACTER(LEN=32)             :: datestring
    TYPE(t_patch), POINTER        :: patch_2d
    INTEGER                       :: jstep0,levels
    REAL(wp)                      :: delta_z
    LOGICAL                       :: lwrite_restart

    TYPE(timedelta), POINTER :: model_time_step => NULL()
    
    CLASS(t_RestartDescriptor), POINTER :: restartDescriptor
    CHARACTER(LEN=max_char_length), PARAMETER :: method_name = 'mo_ocean_testbed_modules:test_output'

    TYPE(datetime), POINTER                      :: current_date
    !------------------------------------------------------------------
    patch_2D      => patch_3d%p_patch_2d(1)

    levels = 40
    ! IF (ltimer) CALL timer_start(timer_total)
    CALL timer_start(timer_total)

    CALL datetimeToString(this_datetime, datestring)

    jstep0                                                               = 0
    !------------------------------------------------------------------
    ! write initial
    ! this is done 
      IF (output_mode%l_nml) THEN
        CALL write_initial_ocean_timestep(patch_3D,p_os(n_dom),surface_fluxes,p_ice,hamocc_state, operators_coefficients)
      ENDIF

    restartDescriptor => createRestartDescriptor("oce")

    ! timeloop
    DO jstep = (jstep0+1), (jstep0+nsteps)
      CALL datetimeToString(this_datetime, datestring)
      WRITE(message_text,'(a,i10,2a)') '  Begin of timestep =',jstep,'  datetime:  ', datestring
      CALL message (TRIM(method_name), message_text)

      ! Set model time.
        model_time_step => newTimedelta('+', 0, 0, 0, 0, 0, NINT(dtime), 0)
        this_datetime = this_datetime + model_time_step
        CALL deallocateTimedelta(model_time_step) 
!      CALL add_time(dtime,0,0,0,this_datetime)


      !------------------------------------------------------------------------
      ! call surface model
      CALL update_surface_flux(patch_3D, p_os(n_dom), p_as, p_ice, atmos_fluxes, surface_fluxes, jstep, this_datetime, &
        &  operators_coefficients)

      !------------------------------------------------------------------------
      ! output: TODO not working for 3d prognostics

      ! add values to output field

      p_os(n_dom)%p_prog(nnew(1))%tracer = p_os(n_dom)%p_prog(nold(1))%tracer
      p_os(n_dom)%p_prog(nnew(1))%h      = p_os(n_dom)%p_prog(nold(1))%h
      p_os(n_dom)%p_diag%t               = p_os(n_dom)%p_prog(nold(1))%tracer(:,:,:,1)
      p_os(n_dom)%p_diag%s               = p_os(n_dom)%p_prog(nold(1))%tracer(:,:,:,2)
      p_os(n_dom)%p_diag%h               = p_os(n_dom)%p_prog(nold(1))%h
      ! add noise {{{     
      CALL add_random_noise_global(in_subset=patch_2D%cells%all, &
        &  in_var=p_os(n_dom)%p_acc%tracer(:,:,:,1), &
        & start_level=1,end_level=levels,noise_scale=10.0_wp,debug=.FALSE.)
      CALL add_random_noise_global(in_subset=patch_2D%cells%all, &
        &  in_var=p_os(n_dom)%p_acc%tracer(:,:,:,2), &
        & start_level=1,end_level=levels,noise_scale=10.0_wp,debug=.FALSE.)
      CALL add_random_noise_global(in_subset=patch_2D%cells%all, &
        &  in_var=p_os(n_dom)%p_acc%u(:,:,:), &
        & start_level=1,end_level=levels,noise_scale=10.0_wp,debug=.FALSE.)
      CALL add_random_noise_global(in_subset=patch_2D%cells%all, &
        &  in_var=p_os(n_dom)%p_acc%v(:,:,:), &
        & start_level=1,end_level=levels,noise_scale=10.0_wp,debug=.FALSE.)
      CALL add_random_noise_global(in_subset=patch_2D%cells%all, &
        &  in_var=p_os(n_dom)%p_acc%w(:,:,:), &
        & start_level=1,end_level=levels,noise_scale=10.0_wp,debug=.FALSE.)
                                                                                                
      IF (no_tracer>=1) THEN
        CALL calc_potential_density( patch_3d,                            &
          & p_os(n_dom)%p_prog(nold(1))%tracer,                       &
          & p_os(n_dom)%p_diag%rhopot )
          
        ! calculate diagnostic barotropic stream function
        CALL calc_psi (patch_3d, p_os(n_dom)%p_diag%u(:,:,:),         &
          & patch_3D%p_patch_1d(1)%prism_thick_c(:,:,:),                  &
          & p_os(n_dom)%p_diag%u_vint, this_datetime)
        CALL dbg_print('calc_psi: u_vint' ,p_os(n_dom)%p_diag%u_vint, debug_string, 3, in_subset=patch_2d%cells%owned)
      END IF
      ! update accumulated vars
      CALL update_ocean_statistics(p_os(1),&
        & surface_fluxes, &
        & patch_2d%cells%owned,&
        & patch_2d%edges%owned,&
        & patch_2d%verts%owned,&
        & n_zlev,p_phys_param=physics_parameters)
      CALL calc_fast_oce_diagnostics( patch_2d,      &
        & patch_3d%p_patch_1d(1)%dolic_c, &
        & patch_3d%p_patch_1d(1)%prism_thick_c, &
        & patch_3d%p_patch_1d(1)%zlev_m, &
        & p_os(n_dom)%p_diag)
      ! }}}
      CALL output_ocean( patch_3D,   &
        &                p_os(n_dom),&
        &                this_datetime,   &
        &                surface_fluxes,  &
        &                p_ice,      &
        &                hamocc_state,      &
        &                jstep, jstep0)

      CALL update_time_indices(n_dom)
      ! write a restart or checkpoint file
      lwrite_restart = (nsteps == INT(time_config%dt_restart/dtime))
      IF (MOD(jstep,n_checkpoints())==0 .OR. ((jstep==(jstep0+nsteps)) .AND. lwrite_restart)) THEN

          CALL restartDescriptor%updatePatch(patch_2d, &
                                            &opt_nice_class=1, &
                                            &opt_ocean_zlevels=n_zlev, &
                                            &opt_ocean_zheight_cellmiddle = patch_3d%p_patch_1d(1)%zlev_m(:), &
                                            &opt_ocean_zheight_cellinterfaces = patch_3d%p_patch_1d(1)%zlev_i(:))
          CALL restartDescriptor%writeRestart(this_datetime, jstep)

        END IF
    END DO

    CALL deleteRestartDescriptor(restartDescriptor)

    CALL timer_stop(timer_total)
  END SUBROUTINE test_output

  SUBROUTINE test_events( patch_3d, ocean_state, p_ext_data,  &
    & this_datetime, surface_fluxes, p_sfc, p_phys_param, &
    & p_as, p_atm_f, sea_ice, &
    & hamocc_state,operators_coefficients,solvercoeff_sp)
    
    TYPE(t_patch_3d ),TARGET, INTENT(inout)          :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET, INTENT(inout) :: ocean_state(n_dom)
    TYPE(t_external_data), TARGET, INTENT(in)        :: p_ext_data(n_dom)
    TYPE(datetime), POINTER                          :: this_datetime
    TYPE(t_sfc_flx)                                  :: surface_fluxes
    TYPE(t_ocean_surface)                            :: p_sfc
    TYPE(t_ho_params)                                :: p_phys_param
    TYPE(t_atmos_for_ocean),  INTENT(inout)          :: p_as
    TYPE(t_atmos_fluxes ),    INTENT(inout)          :: p_atm_f
    TYPE(t_sea_ice),          INTENT(inout)          :: sea_ice
    TYPE(t_hamocc_state),     INTENT(inout)          :: hamocc_state
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
      & routine = 'mo_ocean_testbed_modules:test_output'

    TYPE(eventGroup), POINTER           :: checkpointEventGroup => NULL()
    TYPE(timedelta), POINTER            :: model_time_step => NULL()
    TYPE(datetime), POINTER             :: mtime_current   => NULL()
    TYPE(datetime), POINTER             :: eventRefDate    => NULL(), eventStartDate  => NULL(), eventEndDate    => NULL()
    TYPE(timedelta), POINTER            :: eventInterval   => NULL()
    TYPE(event), POINTER                :: checkpointEvent => NULL()
    TYPE(event), POINTER                :: restartEvent    => NULL()
    
    INTEGER                             :: checkpointEvents, ierr
    LOGICAL                             :: lwrite_checkpoint, lret

    CHARACTER(LEN=MAX_DATETIME_STR_LEN)   :: dstring
    CHARACTER(len=MAX_MTIME_ERROR_STR_LEN):: errstring

    CLASS(t_RestartDescriptor), POINTER :: restartDescriptor
    TYPE(t_RestartAttributeList), POINTER :: restartAttributes
    TYPE(datetime), POINTER               :: current_date
    !------------------------------------------------------------------
    patch_2D      => patch_3d%p_patch_2d(1)

    IF (n_dom > 1 ) THEN
      CALL finish(TRIM(routine), ' N_DOM > 1 is not allowed')
    END IF
    jg = n_dom

    jstep0 = 0
    restartAttributes => getAttributesForRestarting()
    IF (ASSOCIATED(restartAttributes)) THEN
      ! get start counter for time loop from restart file:
      jstep0 = restartAttributes%getInteger("jstep")
    END IF
    IF (isRestart() .AND. mod(nold(jg),2) /=1 ) THEN
      ! swap the g_n and g_nm1
      CALL update_time_g_n(ocean_state(jg))
    ENDIF

    ! set events, group and the events

    CALL message('','')

    eventRefDate   => time_config%tc_exp_refdate
    eventStartDate => time_config%tc_exp_startdate
    eventEndDate   => time_config%tc_exp_stopdate

    ! use start/end setup from the restart
    IF (isRestart()) THEN
      eventRefDate   => time_config%tc_startdate
      eventStartDate => time_config%tc_startdate
      eventEndDate   => time_config%tc_stopdate
    ENDIF


    ! create an event manager, ie. a collection of different events
    CALL initEventManager(time_config%tc_exp_refdate)

    ! --- create an event group for checkpointing and restart
    checkpointEvents =  addEventGroup('checkpointEventGroup')
    checkpointEventGroup => getEventGroup(checkpointEvents)
    
    ! --- --- create checkpointing event
    eventInterval  => time_config%tc_dt_checkpoint
    checkpointEvent => newEvent('checkpoint', eventRefDate, eventStartDate, eventEndDate, eventInterval, errno=ierr)
    IF (ierr /= no_Error) THEN
       CALL mtime_strerror(ierr, errstring)
       CALL finish('perform_ho_timeloop', errstring)
    ENDIF
    lret = addEventToEventGroup(checkpointEvent, checkpointEventGroup)

    ! --- --- create restart event, ie. checkpoint + model stop
    eventInterval  => time_config%tc_dt_restart
    restartEvent => newEvent('restart', eventRefDate, eventStartDate, eventEndDate, eventInterval, errno=ierr)
    IF (ierr /= no_Error) THEN
       CALL mtime_strerror(ierr, errstring)
       CALL finish('perform_ho_timeloop', errstring)
    ENDIF
    lret = addEventToEventGroup(restartEvent, checkpointEventGroup)

    CALL printEventGroup(checkpointEvents)

    ! set time loop properties
    model_time_step => time_config%tc_dt_model

    mtime_current => this_datetime
    
    CALL message('','')
    CALL datetimeToString(mtime_current, dstring)
    WRITE(message_text,'(a,a)') 'Start date of this run: ', dstring
    CALL message('',message_text)
    CALL datetimeToString(time_config%tc_stopdate, dstring)
    WRITE(message_text,'(a,a)') 'Stop date of this run:  ', dstring
    CALL message('',message_text)
    CALL message('','')


    !------------------------------------------------------------------
    ! call the dynamical core: start the time loop
    !------------------------------------------------------------------
    CALL timer_start(timer_total)

    restartDescriptor => createRestartDescriptor("oce")
    
    jstep = jstep0
    TIME_LOOP: DO 

      ! update model date and time mtime based
      mtime_current = mtime_current + model_time_step
      jstep = jstep + 1
      IF (mtime_current > time_config%tc_stopdate) then
#ifdef _MTIME_DEBUG
        ! consistency check: compare step counter to expected end step
        if (jstep /= (jstep0+nsteps)) then
           call finish(routine, 'Step counter does not match expected end step: '//int2string(jstep,'(i0)')&
               &//' /= '//int2string((jstep0+nsteps),'(i0)'))
        end if
#endif
        ! leave time loop
        EXIT TIME_LOOP
      END IF

      CALL datetimeToString(mtime_current, datestring)
      WRITE(message_text,'(a,i10,2a)') '  Begin of timestep =',jstep,'  datetime:  ', datestring
      CALL message (TRIM(routine), message_text)
            
      if (jstep == jstep0 ) &
        & CALL update_height_depdendent_variables( patch_3d, ocean_state(jg), p_ext_data(jg), operators_coefficients, &
        & solvercoeff_sp)
      
      if ( jstep == jstep0 ) CALL calc_scalar_product_veloc_3d( patch_3d,  &
        & ocean_state(jg)%p_prog(nold(1))%vn,         &
        & ocean_state(jg)%p_diag,                     &
        & operators_coefficients)
      
      !In case of a time-varying forcing:
      ! update_surface_flux or update_ocean_surface has changed p_prog(nold(1))%h, SST and SSS
      if (jstep == jstep0 ) then
        IF (surface_module == 1) THEN
          CALL update_surface_flux( patch_3d, ocean_state(jg), p_as, sea_ice, p_atm_f, surface_fluxes, &
            & jstep, mtime_current, operators_coefficients)
        ELSEIF (surface_module == 2) THEN
          CALL update_ocean_surface( patch_3d, ocean_state(jg), p_as, sea_ice, p_atm_f, surface_fluxes, p_sfc, &
            & jstep, mtime_current, operators_coefficients)
        ENDIF
      endif

    
      if (jstep == jstep0 ) then
        CALL update_height_depdendent_variables( patch_3d, ocean_state(jg), p_ext_data(jg), operators_coefficients, &
          & solvercoeff_sp)
      endif

      !---------------------------------------------------------------------

      if (jstep == jstep0 ) &
        & CALL update_ho_params(patch_3d, ocean_state(jg), p_as%fu10, sea_ice%concsum, p_phys_param, operators_coefficients)

      !------------------------------------------------------------------------
      ! solve for new free surface
      !CALL solve_free_surface_eq_ab (patch_3d, ocean_state(jg), p_ext_data(jg), &
      !  & surface_fluxes, p_phys_param, jstep, operators_coefficients, solvercoeff_sp, return_status)!, p_int(jg))
      !IF (return_status /= 0) THEN
      ! CALL output_ocean(              &
      !   & patch_3d=patch_3d,          &
      !   & ocean_state=ocean_state,    &
      !   & this_datetime=mtime_current, &
      !   & surface_fluxes=surface_fluxes, &
      !   & sea_ice=sea_ice,            &
      !   & hamocc=hamocc_state,        &
      !   & jstep=jstep, jstep0=jstep0, &
      !   & force_output=.true.)
      !  CALL finish(TRIM(routine), 'solve_free_surface_eq_ab  returned error')
      !ENDIF
      

      !------------------------------------------------------------------------
      ! Step 4: calculate final normal velocity from predicted horizontal
      ! velocity vn_pred and updated surface height
      !CALL calc_normal_velocity_ab(patch_3d, ocean_state(jg),&
      ! & operators_coefficients, solvercoeff_sp,  p_ext_data(jg), p_phys_param)

      !------------------------------------------------------------------------
      ! Step 5: calculate vertical velocity from continuity equation under
      ! incompressiblity condition in the non-shallow-water case
      !CALL calc_vert_velocity( patch_3d, ocean_state(jg),operators_coefficients)
      !------------------------------------------------------------------------

      !------------------------------------------------------------------------
      ! Step 6: transport tracers and diffuse them
      !IF (no_tracer>=1) THEN
      !  CALL advect_ocean_tracers( patch_3d, ocean_state(jg), p_phys_param,&
      !    & surface_fluxes,&
      !    & operators_coefficients,&
      !    & jstep)
      !ENDIF

      ! perform accumulation for special variables

      ! update accumulated vars

      CALL output_ocean( patch_3d, ocean_state, &
        &                mtime_current,              &
        &                surface_fluxes,             &
        &                sea_ice,                 &
        &                hamocc_state,            &
        &                jstep, jstep0)
      
  
      ! Shift time indices for the next loop
      ! this HAS to ge into the restart files, because the start with the following loop
      CALL update_time_indices(jg)
  
      ! update intermediate timestepping variables for the tracers
      CALL update_time_g_n(ocean_state(jg))

!       CALL message('','')
      ! trigger creation of a restart file ...
      lwrite_checkpoint = .FALSE.      
      IF ( &
           !   ... CASE A: if normal checkpoint cycle has been reached ...
           &       isCurrentEventActive(checkpointEvent, mtime_current)               &
           !          or restart cycle has been reached, i.e. checkpoint+model stop
           & .OR.  isCurrentEventActive(restartEvent, mtime_current)) THEN
        lwrite_checkpoint = .TRUE.
      ENDIF

      ! if this is the first timestep (it cannot occur), or output is disabled, do not write the restart
      IF ( (time_config%tc_startdate == mtime_current) .OR. output_mode%l_none ) THEN
        lwrite_checkpoint = .FALSE.
      ENDIF

!       CALL message('','')
      
      ! write a restart or checkpoint file
!      IF (MOD(jstep,n_checkpoints())==0) THEN
      IF (lwrite_checkpoint) THEN
          CALL restartDescriptor%updatePatch(patch_2d, &
               &                             opt_nice_class=1, &
               &                             opt_ocean_zlevels=n_zlev, &
               &                             opt_ocean_zheight_cellmiddle = patch_3d%p_patch_1d(1)%zlev_m(:), &
               &                             opt_ocean_zheight_cellinterfaces = patch_3d%p_patch_1d(1)%zlev_i(:))
          CALL restartDescriptor%writeRestart(mtime_current, jstep)
      END IF


    ENDDO time_loop
    
    IF (write_last_restart .and. .not. lwrite_checkpoint) THEN
      CALL restartDescriptor%updatePatch(patch_2d, &
           &                             opt_nice_class=1, &
           &                             opt_ocean_zlevels=n_zlev, &
           &                             opt_ocean_zheight_cellmiddle = patch_3d%p_patch_1d(1)%zlev_m(:), &
           &                             opt_ocean_zheight_cellinterfaces = patch_3d%p_patch_1d(1)%zlev_i(:))
      CALL restartDescriptor%writeRestart(mtime_current, jstep)
    END IF
  
    CALL timer_stop(timer_total)
  END SUBROUTINE test_events

  SUBROUTINE test_surface_flux( patch_3d, p_os, &
    & this_datetime, surface_fluxes,         &
    & p_as, atmos_fluxes, p_ice, operators_coefficients)
    
    TYPE(t_patch_3d ),TARGET, INTENT(inout)          :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET, INTENT(inout) :: p_os(n_dom)
    TYPE(datetime), POINTER                          :: this_datetime
    TYPE(t_sfc_flx)                                  :: surface_fluxes
    TYPE(t_atmos_for_ocean),  INTENT(inout)          :: p_as
    TYPE(t_atmos_fluxes ),    INTENT(inout)          :: atmos_fluxes
    TYPE(t_sea_ice),          INTENT(inout)          :: p_ice
    TYPE(t_operator_coeff),   INTENT(inout)          :: operators_coefficients
    
    ! local variables
    TYPE (t_hamocc_state)        :: hamocc_State
    REAL(wp), DIMENSION(nproma,patch_3D%p_patch_2D(1)%alloc_cell_blocks) &
      &                                              :: draft, &
      &                                                 saltBefore, saltAfter, saltBudget, &
      &                                                 salinityBefore, salinityAfter, salinityBudget, &
      &                                                 zUnderIceBefore, zUnderIceAfter
    INTEGER :: jstep
    INTEGER :: block, cell, cellStart,cellEnd
    TYPE(t_subset_range), POINTER :: subset
    !INTEGER :: ocean_statistics
    !LOGICAL                         :: l_outputtime
    CHARACTER(LEN=32)               :: datestring
    TYPE(t_patch), POINTER :: patch_2d
    INTEGER :: jstep0,computation_type
    REAL(wp) :: delta_z

    TYPE(timedelta), POINTER :: model_time_step => NULL()
    
    CHARACTER(LEN=max_char_length), PARAMETER :: &
      & method_name = 'mo_ocean_testbed_modules:test_sea_ice'
    !------------------------------------------------------------------
    patch_2D      => patch_3d%p_patch_2d(1)

    ! IF (ltimer) CALL timer_start(timer_total)
    CALL timer_start(timer_total)

    CALL datetimeToString(this_datetime, datestring)

    jstep0                                                               = 0
    draft(1:nproma,1:patch_3D%p_patch_2D(1)%alloc_cell_blocks)           = 0.0_wp
    saltBefore(1:nproma,1:patch_3D%p_patch_2D(1)%alloc_cell_blocks)      = 0.0_wp
    saltAfter(1:nproma,1:patch_3D%p_patch_2D(1)%alloc_cell_blocks)       = 0.0_wp
    saltBudget(1:nproma,1:patch_3D%p_patch_2D(1)%alloc_cell_blocks)      = 0.0_wp
    salinityBefore(1:nproma,1:patch_3D%p_patch_2D(1)%alloc_cell_blocks)  = 0.0_wp
    salinityAfter(1:nproma,1:patch_3D%p_patch_2D(1)%alloc_cell_blocks)   = 0.0_wp
    salinityBudget(1:nproma,1:patch_3D%p_patch_2D(1)%alloc_cell_blocks)  = 0.0_wp
    zUnderIceBefore(1:nproma,1:patch_3D%p_patch_2D(1)%alloc_cell_blocks) = 0.0_wp
    zUnderIceAfter(1:nproma,1:patch_3D%p_patch_2D(1)%alloc_cell_blocks)  = 0.0_wp
    !------------------------------------------------------------------

    ! x-check if the saltbudget is correctly computed {{{
    ! p_os(n_dom)%p_prog(nold(1))%tracer(:,:,:,1) = 30.0_wp
    ! p_os(n_dom)%p_prog(nold(1))%tracer(:,1,:,2) = 30.0_wp
    ! p_os(n_dom)%p_prog(nold(1))%h(:,:) =  0.0_wp
    ! p_ice%hi(:,:,:) = 0.0_wp
    ! p_ice%conc(:,:,:) = 0.0_wp
    ! p_ice%concSum(:,:) = 0.0_wp
    ! }}}

    CALL dbg_print('TB sfcFlx: hi' ,p_ice%hi          ,debug_string, 4, in_subset=patch_2D%cells%owned)

    !  Overwrite init:
    ! p_ice%hs(:,1,:) = 0.0_wp

  ! draft(:,:)           = (rhos * p_ice%hs(:,1,:) + rhoi * p_ice%hi(:,1,:)) / OceanReferenceDensity
  ! p_ice%zUnderIce(:,:) = patch_3d%p_patch_1d(1)%prism_thick_flat_sfc_c(:,1,:) +  p_os(n_dom)%p_prog(nold(1))%h(:,:) &
  !   &                  - draft(:,:) * p_ice%conc(:,1,:)

    CALL dbg_print('sfcflx: draft    ' ,draft          ,debug_string, 4, in_subset=patch_2D%cells%owned)
    CALL dbg_print('sfcflx: zUnderIce' ,p_ice%zUnderIce,debug_string, 4, in_subset=patch_2D%cells%owned)

    !------------------------------------------------------------------
    ! write initial
    ! this is done 
!     IF (output_mode%l_nml) THEN
!       CALL write_initial_ocean_timestep(patch_3D,p_os(n_dom),surface_fluxes,p_ice)
!     ENDIF


    ! make sure, that h is zero at start
    !p_os(n_dom)%p_prog(nold(1))%h(:,:) = 0.0_wp  !  do not change h

    ! timeloop
    DO jstep = (jstep0+1), (jstep0+nsteps)
      CALL message('test_surface_flux','IceBudget === BEGIN TIMESTEP ======================================================&
        &==============================================')
    

      CALL datetimeToString(this_datetime, datestring)
      WRITE(message_text,'(a,i10,2a)') '  Begin of timestep =',jstep,'  datetime:  ', datestring
      CALL message (TRIM(method_name), message_text)

      ! Set model time.
      model_time_step => newTimedelta('+', 0, 0, 0, 0, 0, NINT(dtime), 0)
      this_datetime = this_datetime + model_time_step
      CALL deallocateTimedelta(model_time_step) 
!      CALL add_time(dtime,0,0,0,this_datetime)

      ! print out 3d salinity
      CALL dbg_print('IceBudget: saltinity BEFORE' ,&
        &            p_os(n_dom)%p_prog(nold(1))%tracer(:,:,:,2) ,&
        &            debug_string, 4, in_subset=patch_2D%cells%owned)
      CALL dbg_print('IceBudget: zUnderIce  BEFORE' ,&
        &             p_ice%zUnderIce(:,:),debug_string,4,in_subset=patch_2D%cells%owned)

      !------------------------------------------------------------------------
      ! computation_type = test_param
      computation_type = 5  ! #slo# merging ocean_sea-ice-thermodyn r208xx
      ! BEFOR : {{{
      zUnderIceBefore = p_ice%zUnderIce
      !salinity
      salinityBefore  = p_os(n_dom)%p_prog(nold(1))%tracer(:,1,:,2)
      ! salt
      saltBefore      = salt_content_in_surface(patch_2D, &
        &                                       patch_3d%p_patch_1d(1)%prism_thick_flat_sfc_c(:,1,:),&
        &                                       p_ice, p_os(n_dom),surface_fluxes,zUnderIceBefore,&
        &                                       computation_type=computation_type,info='BEFORE')
      ! liquid water height
      !}}}

      !------------------------------------------------------------------------
      ! call surface model
      CALL update_surface_flux(patch_3D, p_os(n_dom), p_as, p_ice, atmos_fluxes, surface_fluxes, jstep, this_datetime, &
        &  operators_coefficients)
      CALL dbg_print('IceBudget: saltinity  AFTER' ,&
        &             p_os(n_dom)%p_prog(nold(1))%tracer(:,1,:,2),debug_string,4,in_subset=patch_2D%cells%owned)

      ! simplified update of old temperature with surface forcing
      subset => patch_2D%cells%owned
      DO block = subset%start_block, subset%end_block
        CALL get_index_range(subset, block, cellStart, cellEnd)
        DO cell = cellStart, cellEnd
          IF (subset%vertical_levels(cell,block) < 1) CYCLE
          delta_z = patch_3D%p_patch_1D(1)%prism_thick_flat_sfc_c(cell,1,block)+p_os(n_dom)%p_prog(nold(1))%h(cell,block)
            p_os(n_dom)%p_prog(nold(1))%tracer(cell,1,block,1) = p_os(n_dom)%p_prog(nold(1))%tracer(cell,1,block,1) + &
              & (dtime/delta_z ) * surface_fluxes%topBoundCond_Temp_vdiff(cell,block)
        END DO
      END DO
      CALL dbg_print('sfcflx: trac_new ' ,p_os(n_dom)%p_prog(nold(1))%tracer(:,1,:,1),debug_string, 4, &
        &  in_subset=patch_3d%p_patch_2D(1)%cells%owned)

      !------------------------------------------------------------------------
      ! AFTER {{{
      zUnderIceAfter = p_ice%zUnderIce
      saltAfter      = salt_content_in_surface(patch_2D, patch_3d%p_patch_1d(1)%prism_thick_flat_sfc_c(:,1,:),&
        &                                      p_ice, p_os(n_dom),surface_fluxes,zUnderIceAfter,&
        &                                      computation_type=computation_type,info='AFTER')
      salinityAfter  = p_os(n_dom)%p_prog(nold(1))%tracer(:,1,:,2)
      !}}}

      ! print out 3d salinity
      CALL dbg_print('IceBudget: saltinity  AFTER' ,&
        &             p_os(n_dom)%p_prog(nold(1))%tracer(:,:,:,2),debug_string,4,in_subset=patch_2D%cells%owned)
      CALL dbg_print('IceBudget: FrshFlux_TotalSalt  AFTER' ,&
        &             surface_fluxes%FrshFlux_TotalSalt,debug_string,4,in_subset=patch_2D%cells%owned)
      CALL dbg_print('IceBudget: zUnderIce  AFTER' ,&
        &             p_ice%zUnderIce(:,:),debug_string,4,in_subset=patch_2D%cells%owned)
      CALL dbg_print('IceBudget: zUnderIce  DIFF ' ,&
        &             zUnderIceAfter(:,:) - zUnderIceBefore(:,:),debug_string,4,in_subset=patch_2D%cells%owned)

      ! compute budget
      saltBudget     = saltAfter - saltBefore        ! this discribes the saltbudget in kg, which has to be zero
      salinityBudget = salinityAfter - salinityBefore ! is not allowed to be changed by the sea ice model

      ! SALT-FRESHWATER-CHECK:
      ! (1) check if the FreshwaterFlux (created by the ice model only) is consistent with the zUnderIce variable:
      !     (11) apply the freshwater flux to the upper most cell
      !     (12) for that: switch off all external fresh water fluxes
      !     (12)           ignore snow, i.e. set hs:=0 (no possible snow-to-ice conversion)
      !     (13) check, if the new height correspondes with zUnderIce
      !     PROBLEM: h is based on liquid representatoin of ice, a fresh water flux
      !     PROBLEM: the total ich volume change is put into the
      !     FrshFlux_VolumeIce and this is applied to the h. This seems to be
      !     wrong, because the height does not change by ice groth and melt
      ! (2) check salt content
      !     (21) apply the saltinityFux to the uppermost cell
      !          with (a) constant thickness
      !               (b) zUnderIce
      !          and compute the total salt content in each grid cell: has to be constant over time

      !------------------------------------------------------------------------
      ! output: TODO not working for 3d prognostics

      ! add values to output field
      subset => patch_2D%cells%owned
      DO block = subset%start_block, subset%end_block
        CALL get_index_range(subset, block, cellStart, cellEnd)
        DO cell = cellStart, cellEnd
          IF (subset%vertical_levels(cell,block) < 1) CYCLE
          p_ice%budgets%salt_00(cell,block) = saltBudget(cell,block)
        ENDDO
      ENDDO

      p_os(n_dom)%p_prog(nnew(1))%tracer = p_os(n_dom)%p_prog(nold(1))%tracer
      p_os(n_dom)%p_prog(nnew(1))%h      = p_os(n_dom)%p_prog(nold(1))%h
      p_os(n_dom)%p_diag%t               = p_os(n_dom)%p_prog(nold(1))%tracer(:,:,:,1)
      p_os(n_dom)%p_diag%s               = p_os(n_dom)%p_prog(nold(1))%tracer(:,:,:,2)
      p_os(n_dom)%p_diag%h               = p_os(n_dom)%p_prog(nold(1))%h
      CALL output_ocean( patch_3D,   &
        &                p_os(n_dom),&
        &                this_datetime,   &
        &                surface_fluxes,  &
        &                p_ice,      &
        &                hamocc_state,      &
        &                jstep, jstep0)

      CALL update_time_indices(n_dom)

      CALL dbg_print('IceBudget: salt     diff',saltBudget ,debug_string, 4, in_subset=patch_2D%cells%owned)
      CALL dbg_print('IceBudget: salt    After',saltAfter ,debug_string, 4, in_subset=patch_2D%cells%owned)
      CALL dbg_print('IceBudget: salt   Before',saltBefore ,debug_string, 4, in_subset=patch_2D%cells%owned)
      CALL dbg_print('IceBudget: salinity diff',salinityBudget,debug_string,4,in_subset=patch_2D%cells%owned)
      CALL dbg_print('IceBudget: salt_00' ,p_ice%budgets%salt_00 ,debug_string, 4, in_subset=patch_2D%cells%owned)
    END DO
    
    CALL timer_stop(timer_total)
    
  END SUBROUTINE test_surface_flux
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !>
  SUBROUTINE test_surface_flux_slo( patch_3d, p_os, &
    & this_datetime, surface_fluxes, p_oce_sfc,        &
    & p_as, atmos_fluxes, p_ice, operators_coefficients)
    
    TYPE(t_patch_3d ),TARGET, INTENT(inout)          :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET, INTENT(inout) :: p_os(n_dom)
    TYPE(datetime),           POINTER                :: this_datetime
    TYPE(t_sfc_flx),          INTENT(inout)          :: surface_fluxes
    TYPE(t_ocean_surface),    INTENT(inout)          :: p_oce_sfc
    TYPE(t_atmos_for_ocean),  INTENT(inout)          :: p_as
    TYPE(t_atmos_fluxes ),    INTENT(inout)          :: atmos_fluxes
    TYPE(t_sea_ice),          INTENT(inout)          :: p_ice
    TYPE(t_operator_coeff),   INTENT(inout)          :: operators_coefficients
    
    ! local variables
    TYPE (t_hamocc_state)        :: hamocc_State
    REAL(wp), DIMENSION(nproma,patch_3D%p_patch_2D(1)%alloc_cell_blocks) &
      &                           :: energyCheck, energyCh2, energySav, energyDiff, energyDits, &
      &                              sstCheck, hCheck, meltdraft, conc_old, sst_old, fwfcheck,  &
      &                              saltBefore, saltAfter, saltBudget
    REAL(wp), POINTER             :: sst(:,:), sss(:,:), flat(:,:)
    REAL(wp)                      :: delta_z, t_base, sst_flux, ice_conc
    INTEGER                       :: jstep, jstep0, jc, jb, jk, start_cell_index, end_cell_index
    INTEGER                       :: budget_type_salt, budget_type_energy
    CHARACTER(LEN=32)             :: datestring
    TYPE(t_patch),        POINTER :: patch_2d
    TYPE(t_subset_range), POINTER :: cells_in_domain

    TYPE(timedelta), POINTER :: model_time_step => NULL()
    
    CHARACTER(LEN=max_char_length), PARAMETER :: &
      & method_name = 'mo_ocean_testbed_modules:test_sea_ice'
    !------------------------------------------------------------------
    patch_2D        => patch_3d%p_patch_2d(1)
    cells_in_domain => patch_2D%cells%in_domain

    flat            => patch_3D%p_patch_1D(1)%prism_thick_flat_sfc_c(:,1,:)

    ! IF (ltimer) CALL timer_start(timer_total)
    CALL timer_start(timer_total)

    !t_base = Tf
    !t_base = -5.0_wp
    t_base = t_heat_base
    energyCheck(:,:) = 0.0_wp
    energySav  (:,:) = 0.0_wp
    energyDiff (:,:) = 0.0_wp
    energyDits (:,:) = 0.0_wp
    conc_old   (:,:) = 0.0_wp
    sst_old    (:,:) = 0.0_wp
    fwfcheck   (:,:) = 0.0_wp

    ! sst after init:
    sstCheck(:,:) = p_os(n_dom)%p_prog(nold(1))%tracer(:,1,:,1)

    ! constant thickness for check: initial h; hi, hs in water equivalent to be added to dz+h
  ! hCheck(:,:)    = flat(:,:) + p_os(n_dom)%p_prog(nold(1))%h(:,:) &
  !   &              + (p_ice%hi(:,1,:)*rhoi + p_ice%hs(:,1,:)*rhos) *p_ice%conc(:,1,:)/OceanReferenceDensity

    ! constant thickness for check: initial h, hi, hs in water equivalent are part of initial water height
  ! hCheck(:,:)    = flat(:,:) + p_os(n_dom)%p_prog(nold(1))%h(:,:)

    ! initialized thickness for check: zUnderIce; hi, hs in water equivalent to be subtracted from dz+h
    hCheck(:,:)    = flat(:,:) + p_os(n_dom)%p_prog(nold(1))%h(:,:) &
      &              - (p_ice%hi(:,1,:)*rhoi + p_ice%hs(:,1,:)*rhos)*p_ice%conc(:,1,:)/OceanReferenceDensity

    ! initial energyCh2 - same as energyCheck 
    ! meltdraft: energy content of ice and snow: ((Tf-t_base)*clw-alf) * draftave
  ! energyCh2 = energy_content_in_surface(patch_2d, flat(:,:), p_os(n_dom)%p_prog(nold(1))%h(:,:), &
  !   &         p_ice, sstCheck(:,:), computation_type=computation_type, info='INITIAL')
  ! draft(:,:)           = (rhos * p_ice%hs(:,1,:) + rhoi * p_ice%hi(:,1,:)) / OceanReferenceDensity
    meltdraft(:,:) = ((Tf-t_base)*clw - alf) * (p_ice%hi(:,1,:)*rhoi + p_ice%hs(:,1,:)*rhos)*p_ice%conc(:,1,:)
    energyCh2(:,:) = (sstCheck(:,:) - t_base) * hCheck(:,:)*OceanReferenceDensity*clw + meltdraft(:,:)

    CALL dbg_print('TB.SfcFlux: heightCH2 INI' ,hCheck         (:,:),debug_string, 2, in_subset=patch_2D%cells%owned)
    CALL dbg_print('TB.SfcFlux: energyCh2 INI' ,energyCh2      (:,:),debug_string, 2, in_subset=patch_2D%cells%owned)
    CALL dbg_print('TB.SfcFlux: meltdraft INI' ,meltdraft      (:,:),debug_string, 2, in_subset=patch_2D%cells%owned)

    !  test for draftave, calculated in ice_init
    hCheck(:,:)    = flat(:,:) + p_os(n_dom)%p_prog(nold(1))%h(:,:) - p_ice%draftave(:,:)
    CALL dbg_print('TB.SfcFlux: heightCH2 tst' ,hCheck         (:,:),debug_string, 4, in_subset=patch_2D%cells%owned)

    CALL datetimeToString(this_datetime, datestring)

    jstep0 = 0
    !------------------------------------------------------------------

    DO jstep = (jstep0+1), (jstep0+nsteps)

      ! update pointer every timestep due to changing nold
      sst             => p_os(n_dom)%p_prog(nold(1))%tracer(:,1,:,1)
      sss             => p_os(n_dom)%p_prog(nold(1))%tracer(:,1,:,2)

      ! p_os(n_dom)%p_prog(nold(1))%h(:,:) = 0.0_wp  !  do not change h

      CALL datetimeToString(this_datetime, datestring)
      WRITE(message_text,'(a,i10,2a)') '  Begin of timestep =',jstep,'  datetime:  ', datestring
      CALL message (TRIM(method_name), message_text)

      ! Set model time.
      model_time_step => newTimedelta('+', 0, 0, 0, 0, 0, NINT(dtime), 0)
      this_datetime = this_datetime + model_time_step
      CALL deallocateTimedelta(model_time_step) 
!      CALL add_time(dtime,0,0,0,this_datetime)

      budget_type_salt   = 5  ! #slo# merging ocean_sea-ice-thermodyn r208xx
      budget_type_energy = 0

      !---  energy  -----------------------------------------------------------
      energyCheck = energy_content_in_surface(patch_2d, flat(:,:), p_os(n_dom)%p_prog(nold(1))%h(:,:), &
        &             p_ice, sst(:,:), computation_type=budget_type_energy, info='BEFORE')

      energysav(:,:) = energyDiff(:,:)

      !---  salt    -----------------------------------------------------------
      saltBefore  = salt_content_in_surface(patch_2D, flat(:,:),p_ice, p_os(n_dom),surface_fluxes, &
        &            p_ice%zUnderIce,computation_type=budget_type_salt,info='BEFORE')

      !---------DEBUG DIAGNOSTICS-------------------------------------------
      CALL dbg_print('TB.SfcFlux: hi        BEF' ,p_ice%hi     (:,:,:),debug_string, 2, in_subset=patch_2D%cells%owned)
      CALL dbg_print('TB.SfcFlux: hs        BEF' ,p_ice%hs     (:,:,:),debug_string, 2, in_subset=patch_2D%cells%owned)
      CALL dbg_print('TB.SfcFlux: conc      BEF' ,p_ice%conc   (:,:,:),debug_string, 2, in_subset=patch_2D%cells%owned)
      CALL dbg_print('TB.SfcFlux: zUnderIce BEF' ,p_ice%zUnderIce(:,:),debug_string, 4, in_subset=patch_2D%cells%owned)
      CALL dbg_print('TB.SfcFlux: sst       BEF' ,sst            (:,:),debug_string, 4, in_subset=patch_2D%cells%owned)
      CALL dbg_print('TB.SfcFlux: energy    BEF' ,energyCheck    (:,:),debug_string, 4, in_subset=patch_2D%cells%owned)
      CALL dbg_print('TB.SfcFlux: sstCheck  BEF' ,sstCheck       (:,:),debug_string, 4, in_subset=patch_2D%cells%owned)
      CALL dbg_print('TB.SfcFlux: energyCh2 BEF' ,energyCh2      (:,:),debug_string, 4, in_subset=patch_2D%cells%owned)
      CALL dbg_print('TB.SfcFlux: energyDiffBEF' ,energyDiff     (:,:),debug_string, 4, in_subset=patch_2D%cells%owned)
      CALL dbg_print('TB.SfcFlux: trac_old  BEF' ,p_os(n_dom)%p_prog(nold(1))%tracer(:,1,:,1),debug_string, 5, &
        &  in_subset=patch_3d%p_patch_2D(1)%cells%owned)
      CALL dbg_print('TB.SfcFlux: sss       BEF' ,sss            (:,:),debug_string, 4, in_subset=patch_2D%cells%owned)
      CALL dbg_print('TB.SfcFlux: saltBefore   ' ,saltBefore     (:,:),debug_string, 2, in_subset=patch_2D%cells%owned)
      !---------------------------------------------------------------------

      conc_old(:,:) = p_ice%conc(:,1,:)
      sst_old (:,:) = sst(:,:)

      !-----------------------------------------------------------------------------------------------------------------
      ! call component
      IF (surface_module == 1) THEN
        CALL update_surface_flux(patch_3D, p_os(n_dom), p_as, p_ice, atmos_fluxes, surface_fluxes, jstep, this_datetime, &
          &  operators_coefficients)
      ELSEIF (surface_module == 2) THEN
        CALL update_ocean_surface(patch_3D, p_os(n_dom), p_as, p_ice, atmos_fluxes, surface_fluxes, p_oce_sfc, &
          &  jstep, this_datetime, operators_coefficients)
      ENDIF
      !-----------------------------------------------------------------------------------------------------------------

      !---  energy  -----------------------------------------------------------
      energyCheck = energy_content_in_surface(patch_2d, flat(:,:), p_os(n_dom)%p_prog(nold(1))%h(:,:), &
        &             p_ice, sst(:,:), computation_type=budget_type_energy, info='UpdFlux')

      !---------DEBUG DIAGNOSTICS-------------------------------------------
      CALL dbg_print('TB.SfcFlux: energy  UpdFl' ,energyCheck    (:,:),debug_string, 2, in_subset=patch_2D%cells%owned)
      CALL dbg_print('TB.SfcFlux: sst AFT UpdFl' ,sst            (:,:),debug_string, 2, in_subset=patch_2D%cells%owned)
      CALL dbg_print('TB.SfcFlux: trc AFT UpdFl' ,p_os(n_dom)%p_prog(nold(1))%tracer(:,1,:,1), &
        &  debug_string, 5, in_subset=patch_2D%cells%owned)
      CALL dbg_print('TB.SfcFlux: Flx AFT UpdFl' ,surface_fluxes%topBoundCond_Temp_vdiff(:,:), &
        &  debug_string, 5, in_subset=patch_2D%cells%owned)
      !---------------------------------------------------------------------

      ! simplified update of old temperature with surface forcing
    ! DO jb = cells_in_domain%start_block, cells_in_domain%end_block
    !   CALL get_index_range(cells_in_domain, jb, start_cell_index, end_cell_index)
    !   DO jc = start_cell_index, end_cell_index
    !     delta_z     = patch_3D%p_patch_1D(1)%prism_thick_flat_sfc_c(jc,1,jb)+p_os(n_dom)%p_prog(nold(1))%h(jc,jb)
    !     ! now correct delta_z: use zunderIce
    !     delta_z     = p_ice%zUnderIce(jc,jb)
    !     DO jk = 1, MIN(patch_3d%p_patch_1d(1)%dolic_c(jc,jb),1)  ! this at most should be 1
    !       p_os(n_dom)%p_prog(nold(1))%tracer(jc,jk,jb,1) = p_os(n_dom)%p_prog(nold(1))%tracer(jc,jk,jb,1) + &
    !         & (dtime/delta_z ) * surface_fluxes%topBoundCond_Temp_vdiff(jc,jb)
    !       sst(jc,jb) = sst(jc,jb) + (dtime/delta_z ) * surface_fluxes%topBoundCond_Temp_vdiff(jc,jb)
    !     END DO
    !   END DO
    ! END DO

      !---  energy  -----------------------------------------------------------
      energyCheck = energy_content_in_surface(patch_2d, flat(:,:), p_os(n_dom)%p_prog(nold(1))%h(:,:), &
        &             p_ice, sst(:,:), computation_type=budget_type_energy, info='UpdSST')

      ! check energy input via atmospheric fluxes into surface layer - all atmos_fluxes enter, no flux leaves
      DO jb = cells_in_domain%start_block, cells_in_domain%end_block
        CALL get_index_range(cells_in_domain, jb, start_cell_index, end_cell_index)
        DO jc = start_cell_index, end_cell_index
          DO jk = 1, MIN(patch_3d%p_patch_1d(1)%dolic_c(jc,jb),1)  ! this at most should be 1
  !         SELECT CASE (atmos_flux_analytical_type)
  !         CASE(101,103)
  !           IF (use_ice_concCh) ice_conc = conc_old(jc,jb)
              !  shortwave flux on ice and ocean - whole flux enters the system, no ice concentration matters
              sst_flux = atmos_fluxes%swnetw(jc,jb)+atmos_fluxes%lwnetw(jc,jb) &
                       + atmos_fluxes%latw  (jc,jb)+atmos_fluxes%sensw (jc,jb)
  !  old adjustments, not corrected
  !           IF ( atmos_fluxes%swnetw(jc,jb) == atmos_fluxes%swnet(jc,1,jb) ) THEN
  !             sst_flux = atmos_fluxes%swnetw(jc,jb)
  !             sst_flux = (p_ice%Qtop(jc,1,jb)+p_ice%Qbot(jc,1,jb))*ice_conc + (1.0_wp-ice_conc)*atmos_fluxes%swnetw(jc,jb)
  !           ENDIF
              !  shortwave flux on ice only (Qtop and Qbot calculated)
  !           IF ( atmos_fluxes%swnetw(jc,jb) == 0.0_wp ) THEN
  !             sst_flux = (p_ice%Qtop(jc,1,jb)+p_ice%Qbot(jc,1,jb))*ice_conc
  !           ENDIF
  !         CASE(102)
  !           ! constant initial ice concentration:
  !           ice_conc=init_analytic_conc_param
  !           ! with change in concentration - check with concentration at begin of timestep
  !     !     IF (use_ice_concCh) ice_conc = p_ice%conc(jc,1,jb)
  !           IF (use_ice_concCh) ice_conc = conc_old(jc,jb)
  !           sst_flux = p_ice%Qtop(jc,1,jb)*ice_conc
  !         END SELECT

            ! precipitation: add heat content minus latent heat of frozen rpreci to meltdraft 
            meltdraft(jc,jb) = meltdraft(jc,jb) + ((Tf-t_base)*clw - alf)*atmos_fluxes%rpreci(jc,jb)*dtime* &
              & conc_old(jc,jb)*OceanReferenceDensity

            ! heat: add energy due to fluxes, using old height
            sstCheck(jc,jb)  = sstCheck(jc,jb) + sst_flux*dtime/(clw*OceanReferenceDensity*hCheck(jc,jb))
            ! add energy due to additional water column: precip over open water + precip through ice (rprecw) with sst_old
            fwfcheck(jc,jb)  = p_as%FrshFlux_Precipitation(jc,jb)*(1.0_wp-conc_old(jc,jb))*dtime &
              &                + atmos_fluxes%rprecw(jc,jb)*dtime*conc_old(jc,jb)

            ! calculate theoretical energy for comparison
            !  - energy to check using SST/energy of previous timestep + energy from flux on SST
            !  - additonal freshwater flux yields old SST
            !  - additonal freshwater flux yields new original SST
            !  - when SST and zunderice are changed in one step, a tiny correction term due to meltwater entering
            !    at new SST and not at Tf should be considered here using delhice aus upper_ocean_TS:
            !    T_meltcorr = Delhice*rhoi/OceanReferenceDensity*conc*(sst-tf)*OceanReferenceDensity*clw
            energyCh2(jc,jb) = (sstCheck(jc,jb) - t_base) * hCheck(jc,jb)*OceanReferenceDensity*clw   &  ! new SST with old height
          !   &              + (sst_old (jc,jb) - t_base) * fwfcheck(jc,jb)*OceanReferenceDensity*clw &  ! old SST with added height
              &              + (sst     (jc,jb) - t_base) * fwfcheck(jc,jb)*OceanReferenceDensity*clw &  ! added height receives real SST
          !   &              + (sst     (jc,jb) - t_base) * Delhice*rhoi/OceanReferenceDensity*clw &  ! added height receives real SST
              &              + meltdraft(jc,jb)

            ! needs update of theoretical height and sst for next timestep:
            sstCheck(jc,jb)  = (sstCheck(jc,jb)*hCheck(jc,jb) + sst(jc,jb)*fwfCheck(jc,jb))/(hCheck(jc,jb)+fwfCheck(jc,jb))
            hCheck(jc,jb)    = hCheck(jc,jb) + fwfCheck(jc,jb)  !  update height

          END DO
        END DO
      END DO

      energyDiff(:,:) = energyCheck(:,:) - energyCh2(:,:)
      energyDits(:,:) = energyDiff(:,:)  - energySav(:,:)

      !---  salt    -----------------------------------------------------------
      saltAfter   = salt_content_in_surface(patch_2D, flat(:,:),p_ice, p_os(n_dom),surface_fluxes, &
        &            p_ice%zUnderIce,computation_type=budget_type_salt,info='AFTER')

      ! compute budget
      saltBudget     = saltAfter - saltBefore        ! this discribes the saltbudget in kg, which has to be zero

      !---------DEBUG DIAGNOSTICS-------------------------------------------
      CALL dbg_print('TB.SfcFlux: saltAfter    ' ,saltAfter      (:,:),debug_string, 2, in_subset=patch_2D%cells%owned)
      CALL dbg_print('TB.SfcFlux: saltBudget   ' ,saltBudget     (:,:),debug_string, 2, in_subset=patch_2D%cells%owned)
      CALL dbg_print('TB.SfcFlux: sss      uSST' ,sss            (:,:),debug_string, 2, in_subset=patch_2D%cells%owned)
      CALL dbg_print('TB.SfcFlux: sst      uSST' ,sst            (:,:),debug_string, 2, in_subset=patch_2D%cells%owned)
      CALL dbg_print('TB.SfcFlux: zUnderIc uSST' ,p_ice%zUnderIce(:,:),debug_string, 4, in_subset=patch_2D%cells%owned)
      CALL dbg_print('TB.SfcFlux: conc_old uSST' ,conc_old       (:,:),debug_string, 4, in_subset=patch_2D%cells%owned)
      CALL dbg_print('TB.SfcFlux: energy   uSST' ,energyCheck    (:,:),debug_string, 2, in_subset=patch_2D%cells%owned)
      CALL dbg_print('TB.SfcFlux: sstCheck uSST' ,sstCheck       (:,:),debug_string, 2, in_subset=patch_2D%cells%owned)
      CALL dbg_print('TB.SfcFlux: energyCh2uSST' ,energyCh2      (:,:),debug_string, 2, in_subset=patch_2D%cells%owned)
      CALL dbg_print('TB.SfcFlux: energyDifuSST' ,energyDiff     (:,:),debug_string, 2, in_subset=patch_2D%cells%owned)
      CALL dbg_print('TB.SfcFlux: heightCH2uSST' ,hCheck         (:,:),debug_string, 4, in_subset=patch_2D%cells%owned)
      CALL dbg_print('TB.SfcFlux: fwfcheck uSST' ,fwfcheck       (:,:),debug_string, 4, in_subset=patch_2D%cells%owned)
      CALL dbg_print('TB.SfcFlux: meltdraftuSST' ,meltdraft      (:,:),debug_string, 4, in_subset=patch_2D%cells%owned)
      CALL dbg_print('TB.SfcFlux: energySavuSST' ,energySav      (:,:),debug_string, 4, in_subset=patch_2D%cells%owned)
      CALL dbg_print('TB.SfcFlux: energyDitsSST' ,energyDits     (:,:),debug_string, 2, in_subset=patch_2D%cells%owned)
      CALL dbg_print('TB.SfcFlux: trac_old uSST' ,p_os(n_dom)%p_prog(nold(1))%tracer(:,1,:,1),debug_string, 5, &
        &  in_subset=patch_3d%p_patch_2D(1)%cells%owned)
      !---------------------------------------------------------------------

      ! hack for writing energy/salt budgets on ice%u/v:
   !  p_ice%u(:,:) = energyDits(:,:)
   !  p_ice%v(:,:) = saltBudget(:,:)
   !  CALL dbg_print('TB.SfcFlux: ice%u     END' ,p_ice%u(:,:),debug_string, 3, in_subset=patch_2D%cells%owned)
   !  CALL dbg_print('TB.SfcFlux: ice%v     END' ,p_ice%v(:,:),debug_string, 3, in_subset=patch_2D%cells%owned)

      p_os(n_dom)%p_prog(nnew(1))%tracer = p_os(n_dom)%p_prog(nold(1))%tracer
      p_os(n_dom)%p_prog(nnew(1))%h      = p_os(n_dom)%p_prog(nold(1))%h
      p_os(n_dom)%p_diag%t               = p_os(n_dom)%p_prog(nold(1))%tracer(:,:,:,1)
      p_os(n_dom)%p_diag%s               = p_os(n_dom)%p_prog(nold(1))%tracer(:,:,:,2)
      p_os(n_dom)%p_diag%h               = p_os(n_dom)%p_prog(nold(1))%h
      
      ! update accumulated vars
      CALL update_ocean_statistics(p_os(n_dom), &
        & surface_fluxes,                       &
        & patch_2D%cells%owned,                 &
        & patch_2D%edges%owned,                 &
        & patch_2D%verts%owned,                 &
        & n_zlev)

      CALL output_ocean( patch_3D,   &
        &                p_os(n_dom),&
        &                this_datetime,   &
        &                surface_fluxes,  &
        &                p_ice,      &
        &                hamocc_state,      &
        &                jstep, jstep0)

      CALL update_time_indices(n_dom)

    END DO

    CALL timer_stop(timer_total)
    
  END SUBROUTINE test_surface_flux_slo
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  SUBROUTINE test_neutralcoeff( patch_3d, p_os)
    CHARACTER(LEN=*), PARAMETER ::  routine = "testbed: neutralcoeff"
    
    TYPE(t_patch_3d ),TARGET, INTENT(inout)          :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET, INTENT(inout) :: p_os(n_dom)

    ! local variables
    REAL(wp):: t(n_zlev), s(n_zlev), p(n_zlev), co(n_zlev,2), aob
    REAL(wp):: alph(1:nproma,1:n_zlev,1:patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    REAL(wp):: beta(1:nproma,1:n_zlev,1:patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    TYPE (t_hamocc_state)        :: hamocc_State
    !INTEGER :: jk

    alph(:,:,:) = 0.0_wp
    beta(:,:,:) = 0.0_wp

    CALL calc_neutralslope_coeff( &
      &    patch_3d,              &
      &    p_os(n_dom)%p_prog(nold(1))%tracer(:,:,:,:), &
!       &    p_os(n_dom)%p_prog(nold(1))%h(:,:), &
      &    alph, beta)

    !  test values
    t = 10.0_wp
    s = 40.0_wp
    p = 4000.0_wp    !  4000 dbar = 400 bar
    co = calc_neutralslope_coeff_func_onColumn(t,s,p,n_zlev)
    aob = co(1,1)/co(1,2)

    WRITE(message_text,'(3(a,1pg18.8))') '  Parameter: alpha = ',co(1,1), ' beta = ',co(1,2), ' alpha/beta = ',aob
    CALL message (TRIM(routine), message_text)

  END SUBROUTINE test_neutralcoeff
  !-------------------------------------------------------------------------
  
END MODULE mo_ocean_testbed_modules
