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
  USE mo_ocean_nml,              ONLY: n_zlev, GMRedi_configuration,GMRedi_combined,  GM_only,Redi_only ,Cartesian_Mixing
  USE mo_dynamics_config,        ONLY: nold, nnew
  USE mo_run_config,             ONLY: nsteps, dtime, output_mode, test_mode, test_param
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
  USE mo_sea_ice,                ONLY: salt_content_in_surface
  USE mo_sea_ice_types,          ONLY: t_sfc_flx, t_atmos_fluxes, t_atmos_for_ocean, &
    & t_sea_ice
  USE mo_physical_constants,     ONLY: rhoi, rhos, rho_ref
  USE mo_oce_physics,            ONLY: t_ho_params
  USE mo_oce_GM_Redi,            ONLY: calc_neutralslope_coeff, calc_neutralslope_coeff_func,&
  &                                    prepare_ocean_physics,calc_ocean_physics
  USE mo_time_config,            ONLY: time_config
  USE mo_statistics
  USE mo_util_dbg_prnt,          ONLY: dbg_print
  USE mo_ocean_statistics
  USE mo_ocean_output
  USE mo_parallel_config,        ONLY: nproma
  USE mo_statistics
  USE mo_ocean_testbed_vertical_diffusion
  USE mo_oce_math_operators,        ONLY: div_oce_3d, verticalDiv_scalar_midlevel 
  USE mo_grid_subset,               ONLY: t_subset_range, get_index_range 
  USE mo_oce_diffusion,             ONLY: tracer_diffusion_vertical_implicit,tracer_diffusion_horz
  USE mo_scalar_product,         ONLY: calc_scalar_product_veloc_3d
  USE mo_oce_tracer_transport_horz, ONLY: diffuse_horz
!   USE mo_hydro_ocean_run,        ONLY: write_initial_ocean_timestep

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


      CASE (10)
        CALL ocean_test_GMRedi( patch_3d, ocean_state, &
          & datetime, surface_fluxes, physics_parameters,             &
          & ocean_ice,operators_coefficients)

      CASE DEFAULT
        CALL finish(method_name, "Unknown test_mode")

    END SELECT



  END SUBROUTINE ocean_test_modules
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !>
  SUBROUTINE ocean_test_GMRedi( patch_3d, ocean_state, &
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
    INTEGER :: tracer_index
    
    INTEGER :: jc,level,jb, je
    INTEGER :: z_dolic
    INTEGER :: start_cell_index, end_cell_index
    TYPE(t_subset_range), POINTER :: cells_in_domain, edges_in_domain
    REAL(wp) :: delta_t
    
    REAL(wp) :: z_diff_flux_h(nproma, n_zlev, patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp) :: div_diff_flux_horz(nproma,n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp) :: div_diff_flx_vert(nproma, n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp) :: trac_cart(nproma,n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp) :: div_diff_flux_horz_cart(nproma,n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
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
    CALL datetime_to_string(datestring, datetime)

    ! IF (ltimer) CALL timer_start(timer_total)
    CALL timer_start(timer_total)
        
    time_config%sim_time(:) = 0.0_wp
    jstep0 = 0
    jg=1
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
            CALL verticalDiv_scalar_midlevel( patch_3d, &
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
        CALL add_time(dtime,0,0,0,datetime)
      
        ! Not nice, but the name list output requires this
        time_config%sim_time(1) = time_config%sim_time(1) + dtime
      
        ! update accumulated vars
        CALL update_ocean_statistics(ocean_state(1),     &
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

      END DO
    ! ENDIF!(l_no_time_marching)THEN
    
    CALL timer_stop(timer_total)
    
  END SUBROUTINE ocean_test_GMRedi
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
    
    CHARACTER(LEN=max_char_length), PARAMETER :: &
      & method_name = 'mo_ocean_testbed_modules:test_sea_ice'
    !------------------------------------------------------------------
    patch_2D      => patch_3d%p_patch_2d(1)

    CALL datetime_to_string(datestring, datetime)

    time_config%sim_time(:)                                              = 0.0_wp
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

  ! draft(:,:)           = (rhos * p_ice%hs(:,1,:) + rhoi * p_ice%hi(:,1,:)) / rho_ref
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
    

      CALL datetime_to_string(datestring, datetime)
      WRITE(message_text,'(a,i10,2a)') '  Begin of timestep =',jstep,'  datetime:  ', datestring
      CALL message (TRIM(method_name), message_text)

      ! Set model time.
      CALL add_time(dtime,0,0,0,datetime)

      ! print out 3d salinity
      CALL dbg_print('IceBudget: saltinity BEFORE' ,&
        &            p_os(n_dom)%p_prog(nold(1))%tracer(:,:,:,2) ,&
        &            debug_string, 4, in_subset=patch_2D%cells%owned)
      CALL dbg_print('IceBudget: zUnderIce  BEFORE' ,&
        &             p_ice%zUnderIce(:,:),debug_string,4,in_subset=patch_2D%cells%owned)

      !------------------------------------------------------------------------
      computation_type = test_param
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
      CALL update_surface_flux(patch_3D, p_os(n_dom), p_as, p_ice, atmos_fluxes, surface_fluxes, jstep, datetime, &
        &  operators_coefficients)
      CALL dbg_print('IceBudget: saltinity  AFTER' ,&
        &             p_os(n_dom)%p_prog(nold(1))%tracer(:,1,:,2),debug_string,4,in_subset=patch_2D%cells%owned)

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
      time_config%sim_time(1) = time_config%sim_time(1) + dtime
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
        &                datetime,   &
        &                surface_fluxes,  &
        &                p_ice,      &
        &                jstep, jstep0)

      CALL update_time_indices(n_dom)

      CALL dbg_print('IceBudget: salt     diff',saltBudget ,debug_string, 4, in_subset=patch_2D%cells%owned)
      CALL dbg_print('IceBudget: salt    After',saltAfter ,debug_string, 4, in_subset=patch_2D%cells%owned)
      CALL dbg_print('IceBudget: salt   Before',saltBefore ,debug_string, 4, in_subset=patch_2D%cells%owned)
      CALL dbg_print('IceBudget: salinity diff',salinityBudget,debug_string,4,in_subset=patch_2D%cells%owned)
      CALL dbg_print('IceBudget: salt_00' ,p_ice%budgets%salt_00 ,debug_string, 4, in_subset=patch_2D%cells%owned)
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
!       &    p_os(n_dom)%p_prog(nold(1))%h(:,:), &
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
