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
#include "omp_definitions.inc"
#include "iconfor_dsl_definitions.inc"

MODULE mo_ocean_testbed_modules
  !-------------------------------------------------------------------------
  USE mo_kind,                   ONLY: wp
  USE mo_impl_constants,         ONLY: max_char_length, sea_boundary, zero_coriolis
  USE mo_model_domain,           ONLY: t_patch, t_patch_3d,t_subset_range
  USE mo_grid_config,            ONLY: n_dom, grid_sphere_radius, grid_angular_velocity
  USE mo_math_constants,         ONLY: pi, pi_2, rad2deg, deg2rad, dbl_eps
  USE mo_math_types,             ONLY: t_cartesian_coordinates
  USE mo_ocean_nml,              ONLY: n_zlev, GMRedi_configuration, GMRedi_combined, Cartesian_Mixing, &
    & atmos_flux_analytical_type, no_tracer, OceanReferenceDensity, l_with_vert_tracer_advection, &
    & tracer_update_mode, use_none, l_edge_based
  USE mo_sea_ice_nml,            ONLY: init_analytic_conc_param, t_heat_base
  USE mo_dynamics_config,        ONLY: nold, nnew
  USE mo_run_config,             ONLY: nsteps, dtime, output_mode, test_mode !, test_param
  USE mo_exception,              ONLY: message, message_text, finish
  USE mo_ext_data_types,         ONLY: t_external_data
  !USE mo_io_units,               ONLY: filename_max
  USE mo_timer,                  ONLY: timer_start, timer_stop, timer_total
  USE mo_ocean_ab_timestepping,  ONLY: update_time_indices , &
     & solve_free_surface_eq_ab, calc_vert_velocity
  USE mo_random_util,            ONLY: add_random_noise_global
  USE mo_ocean_types,            ONLY: t_hydro_ocean_state, t_operator_coeff, t_solvercoeff_singleprecision
  USE mo_hamocc_types,           ONLY: t_hamocc_state
  USE mo_restart,                ONLY: t_RestartDescriptor, createRestartDescriptor, deleteRestartDescriptor
  USE mo_restart_attributes,     ONLY: t_RestartAttributeList, getAttributesForRestarting
  USE mo_io_config,              ONLY: n_checkpoints, write_last_restart
  USE mo_operator_ocean_coeff_3d,ONLY: t_operator_coeff, no_primal_edges
  USE mo_ocean_tracer,           ONLY: advect_ocean_tracers
  USE mo_ocean_surface_refactor, ONLY: update_ocean_surface_refactor
  USE mo_ocean_surface_types,    ONLY: t_ocean_surface, t_atmos_for_ocean
  USE mo_sea_ice,                ONLY: salt_content_in_surface, energy_content_in_surface
  USE mo_sea_ice_types,          ONLY: t_atmos_fluxes, t_sea_ice
  USE mo_ice_diagnostics,        ONLY: energy_in_surface, salt_in_surface
  USE mo_physical_constants,     ONLY: rhoi, rhos, clw, alf, Tf
  USE mo_ocean_physics_types,    ONLY: t_ho_params
  USE mo_master_config,          ONLY: isRestart
  USE mo_ocean_GM_Redi,          ONLY: prepare_ocean_physics,calc_ocean_physics
  USE mo_ocean_diagnostics,      ONLY: calc_fast_oce_diagnostics, calc_psi
  USE mo_ocean_thermodyn,        ONLY: calc_potential_density, calculate_density,&
  &                                    calc_neutralslope_coeff_func_onColumn,calc_neutralslope_coeff_func_onColumn_UNESCO
  USE mo_time_config,            ONLY: time_config
  USE mo_statistics
  USE mo_util_dbg_prnt,          ONLY: dbg_print
  USE mo_ocean_statistics
  USE mo_ocean_output
  USE mo_parallel_config,        ONLY: nproma
  USE mo_statistics
  USE mo_ocean_testbed_vertical_diffusion
  USE mo_ocean_math_operators
  USE mo_grid_subset,            ONLY: t_subset_range, get_index_range 
  USE mo_scalar_product,         ONLY: calc_scalar_product_veloc_3d, &
      & map_edges2cell_3d, map_scalar_center2prismtop, map_vec_prismtop2center_on_block, &
      & map_cell2edges_3D, map_edges2edges_viacell_3d_const_z
  USE mo_ocean_tracer_transport_horz, ONLY: diffuse_horz
  USE mo_hydro_ocean_run
  USE mo_var_list
  USE mo_linked_list
  USE mo_cdi
  use mo_cdi_constants
  use mo_zaxis_type
  use mo_cf_convention
  use mo_grib2

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
  USE mo_ocean_physics,         ONLY: update_ho_params

  USE mo_ocean_tracer_transport_types

  !! Needed to test advection of velocity 
  USE mo_ocean_velocity_advection, ONLY: veloc_adv_horz_mimetic, &
      & veloc_adv_vert_mimetic

  USE mo_ocean_tracer_transport_horz, ONLY: advect_horz, diffuse_horz
  USE mo_ocean_tracer_transport_vert, ONLY: advect_flux_vertical
  USE mo_sync,                        ONLY: sync_c, sync_c1, sync_patch_array, sync_patch_array_mult

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ocean_test_modules
  
  CHARACTER(len=12)           :: debug_string = 'testbedMod  '  ! Output of module for 1 line debug
  
  !-------------------------------------------------------------------------
CONTAINS

  !-------------------------------------------------------------------------
  !>
  SUBROUTINE ocean_test_modules( patch_3d, ocean_state,  external_data,  &
    & this_datetime, ocean_surface, physics_parameters,             &
    & oceans_atmosphere, oceans_atmosphere_fluxes, ocean_ice, operators_coefficients, &
    & solvercoeff_sp)

    TYPE(t_patch_3d ), POINTER, INTENT(in)          :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET, INTENT(inout) :: ocean_state(n_dom)
    TYPE(t_external_data), TARGET, INTENT(in)        :: external_data(n_dom)
    TYPE(datetime), POINTER                          :: this_datetime
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
          & this_datetime, ocean_surface, physics_parameters,             &
          & ocean_ice,operators_coefficients)

      CASE (2)
        CALL test_tracer_diffusion_vertical_implicit( patch_3d, ocean_state, physics_parameters,  &
           & operators_coefficients)

      CASE (3)
        CALL test_velocity_diffusion_vert_implicit( patch_3d, ocean_state, physics_parameters,  &
           & operators_coefficients)

      CASE (4)
        CALL finish(method_name, "This test mode was removed - STOP")

      CASE (41)
        CALL test_surface_flux_slo( patch_3d, ocean_state,  &
          & this_datetime, ocean_surface,        &
          & oceans_atmosphere, oceans_atmosphere_fluxes, ocean_ice, operators_coefficients)

      CASE (42)
        CALL test_sea_ice( patch_3d, ocean_state,  &
          & this_datetime, ocean_surface,        &
          & oceans_atmosphere, oceans_atmosphere_fluxes, ocean_ice, operators_coefficients)

      CASE (5)
        CALL test_neutralcoeff( patch_3d, ocean_state)


      CASE (10)
        CALL finish(method_name, "The GMRedi test has been disabled - STOP")
!         CALL ocean_test_GMRedi( patch_3d, ocean_state, &
!           & this_datetime, ocean_surface, physics_parameters,             &
!           & ocean_ice,operators_coefficients)

      CASE (11) ! surface only processing to get quasi output fast
        CALL test_output( patch_3d, ocean_state,  &
          & this_datetime, physics_parameters,                   &
          & oceans_atmosphere, oceans_atmosphere_fluxes, ocean_surface, ocean_ice, hamocc_state,operators_coefficients)

      CASE (12) ! surface only processing to get quasi output fast
        CALL test_events( patch_3d, ocean_state,  &
          & external_data , &
          & this_datetime, &
          & ocean_surface, &
          & physics_parameters, &
          & oceans_atmosphere, &
          & oceans_atmosphere_fluxes, &
          & ocean_ice, &
          & hamocc_state, &
          & operators_coefficients, &
          & solvercoeff_sp)

      CASE (13) ! check the handling to var_list entries wrt. the their output name
        CALL checkVarlistsForOutput(patch_3d, ocean_state(1), &
          & this_datetime, physics_parameters,                   &
          & oceans_atmosphere, oceans_atmosphere_fluxes, ocean_surface, ocean_ice, hamocc_state,operators_coefficients)

      CASE (14)
        CALL checkVarlistKeys(patch_3d%p_patch_2D(1))
    
    
      CASE (91) ! Test zstar transformation of operators 
        CALL test_zstar( patch_3d, ocean_state,  &
          & external_data , &
          & this_datetime, &
          & ocean_surface, &
          & physics_parameters, &
          & oceans_atmosphere, &
          & oceans_atmosphere_fluxes, &
          & ocean_ice, &
          & hamocc_state, &
          & operators_coefficients, &
          & solvercoeff_sp)

      CASE (92) ! Test advection part of momentum  
        CALL test_mom( patch_3d, ocean_state,  &
          & external_data , &
          & this_datetime, &
          & ocean_surface, &
          & physics_parameters, &
          & oceans_atmosphere, &
          & oceans_atmosphere_fluxes, &
          & ocean_ice, &
          & hamocc_state, &
          & operators_coefficients, &
          & solvercoeff_sp)

      
      CASE (93)
        CALL ocean_test_zstar_advection( patch_3d, ocean_state, &
          & this_datetime, ocean_surface, physics_parameters,             &
          & ocean_ice,operators_coefficients)


      CASE DEFAULT
        CALL finish(method_name, "Unknown test_mode")

    END SELECT



  END SUBROUTINE ocean_test_modules
  !-------------------------------------------------------------------------

!   !-------------------------------------------------------------------------
!   !>
!   SUBROUTINE ocean_test_GMRedi( patch_3d, ocean_state, &
!     & this_datetime, ocean_surface, physics_parameters,             &
!     & ocean_ice,operators_coefficients)
!     
!     TYPE(t_patch_3d ),TARGET, INTENT(inout)          :: patch_3d
!     TYPE(t_hydro_ocean_state), TARGET, INTENT(inout) :: ocean_state(n_dom)
!     TYPE(datetime), POINTER                          :: this_datetime
!     TYPE(t_ocean_surface)                            :: ocean_surface
!     TYPE (t_ho_params)                               :: physics_parameters
!     TYPE (t_sea_ice),         INTENT(inout)          :: ocean_ice
!     TYPE(t_operator_coeff),   INTENT(inout)          :: operators_coefficients
!     
!     ! local variables
!     TYPE (t_hamocc_state)        :: hamocc_State
!     INTEGER :: jstep, jg
!     !LOGICAL                         :: l_outputtime
!     CHARACTER(LEN=32)               :: datestring
!     TYPE(t_patch), POINTER :: patch_2d
!     INTEGER :: jstep0 ! start counter for time loop
!     INTEGER :: tracer_index
!     
!     INTEGER :: jc,level,jb
!     INTEGER :: start_cell_index, end_cell_index
!     TYPE(t_subset_range), POINTER :: cells_in_domain, edges_in_domain
!     REAL(wp) :: delta_t
!     
!     REAL(wp) :: z_diff_flux_h(nproma, n_zlev, patch_3d%p_patch_2d(1)%nblks_e)
!     !REAL(wp) :: div_diff_flux_horz(nproma,n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
!     !REAL(wp) :: div_diff_flx_vert(nproma, n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)
!     REAL(wp) :: density_backup(nproma,n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
!     !REAL(wp) :: div_diff_flux_horz_cart(nproma,n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
!     REAL(wp),POINTER :: fu10   (:,:)
!     REAL(wp), POINTER :: concsum(:,:)
! 
!     TYPE(timedelta), POINTER :: model_time_step => NULL()
! 
!     CHARACTER(LEN=max_char_length), PARAMETER :: &
!       & method_name = 'mo_ocean_testbed_modules:ocean_test_advection'
!     !------------------------------------------------------------------
!     tracer_index=1!test is on salinity
! 
!     
!     
!     patch_2D      => patch_3d%p_patch_2d(1)
!    
!     cells_in_domain => patch_2D%cells%in_domain
!     edges_in_domain => patch_2D%edges%in_domain
!     delta_t = dtime
! 
!     !z_diff_flux_h(1:nproma,1:n_zlev,1:patch_3d%p_patch_2d(1)%nblks_e)=0.0_wp
!     !div_diff_flux_horz(1:nproma,1:n_zlev,1:patch_2d%alloc_cell_blocks)=0.0_wp
!     !!div_diff_flx_vert (1:nproma,1:n_zlev,1:patch_2d%alloc_cell_blocks)=0.0_wp 
!     !trac_cart(1:nproma,1:n_zlev, 1:patch_3d%p_patch_2d(1)%alloc_cell_blocks)=0.0_wp
!     !div_diff_flux_horz_cart(1:nproma,1:n_zlev, 1:patch_3d%p_patch_2d(1)%alloc_cell_blocks)=0.0_wp
!     !---------------------------------------------------------------------   
!     CALL datetimeToString(this_datetime, datestring)
! 
!     ! IF (ltimer) CALL timer_start(timer_total)
!     CALL timer_start(timer_total)
!         
!     jstep0 = 0
!     jg=1
!     density_backup=ocean_state(n_dom)%p_diag%rho
!     ocean_state(n_dom)%p_diag%rho_GM=density_backup
!     !------------------------------------------------------------------
!     ! IF(.NOT.l_time_marching)THEN
! 
!       !IF(itestcase_oce==28)THEN
!       DO jstep = (jstep0+1), (jstep0+nsteps)
!       
! !        density_backup=ocean_state(n_dom)%p_diag%rho
! !        ocean_state(n_dom)%p_diag%rho_GM=density_backup
!       
!         CALL update_ho_params(patch_3d, ocean_state(n_dom), fu10, concsum, physics_parameters, operators_coefficients) 
! !        ocean_state(n_dom)%p_diag%rho_GM=density_backup             
! !        ocean_state(n_dom)%p_diag%rho   =density_backup             
! ! ocean_state(n_dom)%p_prog(nold(1))%ocean_tracers(1)%concentration=ocean_state(n_dom)%p_diag%rho_GM       
!  
!         CALL datetimeToString(this_datetime, datestring)
!         WRITE(message_text,'(a,i10,2a)') '  Begin of timestep =',jstep,'  datetime:  ', datestring
!         CALL message (TRIM(method_name), message_text)
! ! physics_parameters%a_tracer_v=k_pot_temp_v
! !          IF(jstep==1)THEN
! !          ocean_state(jg)%p_diag%vn_time_weighted = ocean_state(jg)%p_prog(nold(1))%vn
! !          ocean_state(jg)%p_prog(nnew(1))%vn = ocean_state(jg)%p_prog(nold(1))%vn
! !          ocean_state(jg)%p_diag%w        =  0.0_wp!0.0833_wp!0.025_wp
! !          ocean_state(jg)%p_diag%w(:,:,:) = -0.0833_wp!0.025_wp
! !          ENDIF  
! !       CALL calc_scalar_product_veloc_3d( patch_3d,  &
! !        & ocean_state(n_dom)%p_prog(nold(1))%vn,     &
! !        & ocean_state(n_dom)%p_diag,                 &
!         !calculate some information that is used for all tracers
! 
!         CALL prepare_tracer_transport( patch_3d, &
!           & ocean_state(jg), operators_coefficients)
! !        & operators_coefficients)
! 
! 
! CALL advect_ocean_tracers(patch_3d, ocean_state(n_dom), physics_parameters, operators_coefficients,jstep)
! !IF(GMRedi_configuration/=Cartesian_Mixing)THEN 
! !      CALL prepare_ocean_physics(patch_3d, &
! !        & ocean_state(n_dom),    &
! !        & physics_parameters, &
! !        & operators_coefficients)
! !ENDIF
! !DO tracer_index=1,no_tracer
! !         CALL advect_diffuse_tracer( patch_3d, &
! !           & ocean_state(n_dom)%p_prog(nold(1))%ocean_tracers(tracer_index),&
! !           & ocean_state(n_dom),            &
! !           & operators_coefficients,      &
! !           & ocean_state(n_dom)%p_aux%bc_top_tracer(:,:,tracer_index),   &
! !           & ocean_state(n_dom)%p_aux%bc_bot_tracer,   &
! !           & physics_parameters,         &
! !           & physics_parameters%k_tracer_h(:,:,:,1),             &
! !           & physics_parameters%a_tracer_v(:,:,:,1),             &
! !           & ocean_state(n_dom)%p_prog(nnew(1))%ocean_tracers(tracer_index),&
! !           & tracer_index)
! !END DO
! 
! !      CALL calculate_density( patch_3d,                         &
! !       & ocean_state(n_dom)%p_prog(nold(1))%tracer(:,:,:,1:no_tracer),&
! !       & ocean_state(n_dom)%p_diag%rho(:,:,:) )
! 
! 
! 
! 
!         
! !        IF(GMRedi_configuration/=Cartesian_Mixing)THEN    
! !          CALL prepare_ocean_physics(patch_3d, &
! !            & ocean_state(n_dom),    &
! !            & physics_parameters, &
! !            & operators_coefficients)
! !        ENDIF
! !        
! !!        DO tracer_index=1,2
! !        
! !          IF(GMRedi_configuration/=Cartesian_Mixing)THEN
! !        
! !            CALL calc_ocean_physics( patch_3d, &
! !                                   & ocean_state(n_dom),     &
! !                                   &  physics_parameters,    &
! !                                   &  operators_coefficients,&
! !                                   &  tracer_index)
! !            CALL div_oce_3d( ocean_state(n_dom)%p_diag%GMRedi_flux_horz(:,:,:,tracer_index),&
! !                     &   patch_3D, &
! !                     &   operators_coefficients%div_coeff, &
! !                     &   div_diff_flux_horz )
! !            !vertical div of GMRedi-flux
! !            CALL verticalDiv_scalar_onFullLevels( patch_3d, &
! !                                            & ocean_state(n_dom)%p_diag%GMRedi_flux_vert(:,:,:,tracer_index), &
! !                                            & div_diff_flx_vert)
! !                                   
! !         ELSE
! !          CALL tracer_diffusion_horz(patch_3D,&
! !                                   & ocean_state(n_dom)%p_prog(nold(1))%ocean_tracers(tracer_index)%concentration,&
! !                                   & ocean_state(n_dom), z_diff_flux_h, physics_parameters%k_tracer_h(:,:,:,tracer_index ))
! !
! !            CALL div_oce_3d( z_diff_flux_h,&
! !                     &   patch_3D, &
! !                     &   operators_coefficients%div_coeff, &
! !                     &   div_diff_flux_horz_cart )
! !    
! !         ENDIF
! !
! !         !cart
! !         DO jb = cells_in_domain%start_block, cells_in_domain%end_block
! !           CALL get_index_range(cells_in_domain, jb, start_cell_index, end_cell_index)
! !           DO jc = start_cell_index, end_cell_index
! !
! !              DO level = 1, patch_3d%p_patch_1d(1)%dolic_c(jc,jb)
! !              ocean_state(n_dom)%p_prog(nnew(1))%ocean_tracers(tracer_index-1)%concentration(jc,level,jb) &
! !                & = ocean_state(n_dom)%p_prog(nold(1))%ocean_tracers(tracer_index-1)%concentration(jc,level,jb)-  &
! !                &  (delta_t /  patch_3D%p_patch_1D(1)%prism_thick_c(jc,level,jb))  &
! !                &    * ( - (div_diff_flux_horz_cart(jc,level,jb)))
! !! write(123,*)'details',level,  ocean_state(n_dom)%p_prog(nnew(1))%ocean_tracers(tracer_index)%concentration(jc,level,jb),&
! !! & ocean_state(n_dom)%p_prog(nold(1))%ocean_tracers(tracer_index)%concentration(jc,level,jb),&
! !! & div_diff_flux_horz(jc,level,jb)
! !
! !           ENDDO
! !         END DO
! !       END DO
! !
! !
! !       !GM                               
! !       DO jb = cells_in_domain%start_block, cells_in_domain%end_block
! !         CALL get_index_range(cells_in_domain, jb, start_cell_index, end_cell_index)
! !         DO jc = start_cell_index, end_cell_index
! !
! !           DO level = 1, patch_3d%p_patch_1d(1)%dolic_c(jc,jb)
! !             ocean_state(n_dom)%p_prog(nnew(1))%ocean_tracers(tracer_index)%concentration(jc,level,jb) &
! !                & = ocean_state(n_dom)%p_prog(nold(1))%ocean_tracers(tracer_index)%concentration(jc,level,jb)-  &
! !                &  (delta_t /  patch_3D%p_patch_1D(1)%prism_thick_c(jc,level,jb))  &
! !                &    * ( - (div_diff_flux_horz(jc,level,jb)-div_diff_flx_vert(jc,level,jb)))
! !!  write(123,*)'details',level,  ocean_state(n_dom)%p_prog(nnew(1))%ocean_tracers(tracer_index)%concentration(jc,level,jb),&
! !!  & ocean_state(n_dom)%p_prog(nold(1))%ocean_tracers(tracer_index)%concentration(jc,level,jb),&
! !!  & div_diff_flux_horz(jc,level,jb),div_diff_flx_vert(jc,level,jb),&
! !! &ocean_state(n_dom)%p_aux%slopes_squared(jc,level,jb),&
! !! &ocean_state(n_dom)%p_aux%taper_function_1(jc,level,jb),&
! !! &ocean_state(n_dom)%p_aux%taper_function_2(jc,level,jb)!,&
! !
! !           ENDDO
! !         END DO
! !       END DO
! !
! !
! !      !cart    
! !      CALL tracer_diffusion_vertical_implicit(                         &
! !      & patch_3d,                                                      &
! !      & ocean_state(n_dom)%p_prog(nnew(1))%ocean_tracers(tracer_index-1),&
! !      & physics_parameters%a_tracer_v(:,:,:, tracer_index),            &
! !      & operators_coefficients)
! !
! !          
! !      !GM    
! !      CALL tracer_diffusion_vertical_implicit(                         &
! !      & patch_3d,                                                      &
! !      & ocean_state(n_dom)%p_prog(nnew(1))%ocean_tracers(tracer_index),&
! !      & physics_parameters%a_tracer_v(:,:,:, tracer_index),            &
! !      & operators_coefficients)
! !      
! !      
! !      ocean_state(n_dom)%p_diag%rho=ocean_state(n_dom)%p_prog(nnew(1))%ocean_tracers(tracer_index-1)%concentration&
! !      &-ocean_state(n_dom)%p_prog(nnew(1))%ocean_tracers(tracer_index)%concentration
! !         
!       !END DO       
! 
!         ! One integration cycle finished on the lowest grid level (coarsest
!         ! resolution). Set model time.
! 
!         model_time_step => newTimedelta('+', 0, 0, 0, 0, 0, NINT(dtime), 0)
!         this_datetime = this_datetime + model_time_step
!         CALL deallocateTimedelta(model_time_step)        
! !        CALL add_time(dtime,0,0,0,this_datetime)
!      
!         ! update accumulated vars
! !       CALL update_ocean_statistics(ocean_state(1),     &
! !       & ocean_surface,                                &
! !       & patch_2D%cells%owned,       &
! !       & patch_2D%edges%owned,       &
! !       & patch_2D%verts%owned,       &
! !       & n_zlev)
!           
!         CALL output_ocean( patch_3d, &
!           & ocean_state,             &
!           & this_datetime,                &
!           & ocean_surface,          &
!           & ocean_ice,               &
!           & jstep, jstep0)
! 
!         ! Shift time indices for the next loop
!         ! this HAS to ge into the restart files, because the start with the following loop
!         CALL update_time_indices(jg)
!         ! update intermediate timestepping variables for the tracers
!         ! velocity
! !IF(tracer_index==1)THEN
! !DO level = 1, n_zlev
! !write(0,*)'tracer:rho',&
! !& maxval( ocean_state(n_dom)%p_prog(nold(1))%ocean_tracers(tracer_index)%concentration(:,level,:)),&
! !& minval( ocean_state(n_dom)%p_prog(nold(1))%ocean_tracers(tracer_index)%concentration(:,level,:)),&
! !!& maxval( ocean_state(n_dom)%p_prog(nnew(1))%ocean_tracers(tracer_index-1)%concentration(:,level,:)),&
! !!& minval( ocean_state(n_dom)%p_prog(nnew(1))%ocean_tracers(tracer_index-1)%concentration(:,level,:)),&
! !& maxval( ocean_state(n_dom)%p_diag%rho(:,level,:)),&
! !& minval( ocean_state(n_dom)%p_diag%rho(:,level,:))!,&
! !!& maxval( div_diff_flux_horz(:,level,:)),&
! !!& minval( div_diff_flux_horz(:,level,:))
! !
! !END DO
! !ENDIF
! 
!       END DO
!     ! ENDIF!(l_no_time_marching)THEN
!     
!     CALL timer_stop(timer_total)
!     
!   END SUBROUTINE ocean_test_GMRedi
!   !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !>
  SUBROUTINE ocean_test_advection( patch_3d, ocean_state, &
    & this_datetime, ocean_surface, physics_parameters,             &
    & ocean_ice,operators_coefficients)
    
    TYPE(t_patch_3d), POINTER, INTENT(in)          :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET, INTENT(inout) :: ocean_state(n_dom)
    TYPE(datetime), POINTER                          :: this_datetime
    TYPE(t_ocean_surface)                            :: ocean_surface
    TYPE (t_ho_params)                               :: physics_parameters
    TYPE (t_sea_ice),         INTENT(inout)          :: ocean_ice
    TYPE(t_operator_coeff),   INTENT(in)          :: operators_coefficients
    
    ! local variables
    TYPE (t_hamocc_state)        :: hamocc_State
    INTEGER :: jstep, jg
    !LOGICAL                         :: l_outputtime
    CHARACTER(LEN=32)               :: datestring
    TYPE(t_patch), POINTER :: patch_2d
    INTEGER :: jstep0 ! start counter for time loop
    INTEGER :: i
    TYPE(timedelta), POINTER :: model_time_step => NULL()
    
    !CHARACTER(LEN=filename_max)  :: outputfile, gridfile
    TYPE(t_tracer_collection) , POINTER              :: old_tracer_collection, new_tracer_collection
    TYPE(t_ocean_transport_state)                    :: transport_state

    CHARACTER(LEN=max_char_length), PARAMETER :: &
      & method_name = 'mo_ocean_testbed_modules:ocean_test_advection'
    !------------------------------------------------------------------
    
    patch_2D      => patch_3d%p_patch_2d(1)
    CALL datetimeToString(this_datetime, datestring)

    ! IF (ltimer) CALL timer_start(timer_total)
    CALL timer_start(timer_total)
    
    IF (n_dom > 1 ) THEN
      CALL finish(TRIM(method_name), ' N_DOM > 1 is not allowed')
    END IF
    jg = n_dom

    jstep0 = 0
      
    DO jstep = (jstep0+1), (jstep0+nsteps)
      
        CALL datetimeToString(this_datetime, datestring)
        WRITE(message_text,'(a,i10,2a)') '  Begin of timestep =',jstep,'  datetime:  ', datestring
        CALL message (TRIM(method_name), message_text)
 
        old_tracer_collection => ocean_state(jg)%p_prog(nold(1))%tracer_collection
        new_tracer_collection => ocean_state(jg)%p_prog(nnew(1))%tracer_collection

        IF (no_tracer>=1) THEN

          !! Update vertical velocity w 
          CALL calc_vert_velocity( patch_3d, ocean_state(jg),operators_coefficients)
          !! Update mass_flx_e for tracer advection
          CALL map_edges2edges_viacell_3d_const_z( patch_3d, ocean_state(jg)%p_prog(nold(1))%vn, &
                         & operators_coefficients, ocean_state(jg)%p_diag%mass_flx_e)
          ! fill transport_state
          transport_state%patch_3d    => patch_3d
          transport_state%h_old       => ocean_state(jg)%p_prog(nold(1))%h
          transport_state%h_new       => ocean_state(jg)%p_prog(nnew(1))%h
          transport_state%vn          => ocean_state(jg)%p_prog(nold(1))%vn
          transport_state%w           => ocean_state(jg)%p_diag%w
          transport_state%mass_flux_e => ocean_state(jg)%p_diag%mass_flx_e

          ! fill diffusion coefficients
          old_tracer_collection%tracer(1)%hor_diffusion_coeff => physics_parameters%TracerDiffusion_coeff(:,:,:,1)
          old_tracer_collection%tracer(1)%ver_diffusion_coeff => physics_parameters%a_tracer_v(:,:,:,1)
          DO i = 2, old_tracer_collection%no_of_tracers
            old_tracer_collection%tracer(i)%hor_diffusion_coeff => physics_parameters%TracerDiffusion_coeff(:,:,:,2)
            old_tracer_collection%tracer(i)%ver_diffusion_coeff => physics_parameters%a_tracer_v(:,:,:,2)
          ENDDO
        ENDIF
        !------------------------------------------------------------------------

        CALL advect_ocean_tracers(old_tracer_collection, new_tracer_collection, transport_state, operators_coefficients)

        ! One integration cycle finished on the lowest grid level (coarsest
        ! resolution). Set model time.
        model_time_step => newTimedelta('+', 0, 0, 0, 0, 0, NINT(dtime), 0)
        this_datetime = this_datetime + model_time_step
        CALL deallocateTimedelta(model_time_step) 
          
        CALL output_ocean( patch_3d, &
          & ocean_state,             &
          & this_datetime,                &
          & ocean_surface,          &
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
    
    CALL timer_stop(timer_total)
    
  END SUBROUTINE ocean_test_advection
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !> surface only call, regular output, restart, checkpoints
  SUBROUTINE test_output(patch_3d, p_os,           &
    & this_datetime, physics_parameters, &
    & p_as, atmos_fluxes, p_oce_sfc, p_ice, hamocc_state,operators_coefficients)
    
    TYPE(t_patch_3d), POINTER, INTENT(in) :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET, INTENT(inout) :: p_os(n_dom)
    TYPE(datetime), POINTER                          :: this_datetime
    TYPE (t_ho_params)                               :: physics_parameters
    TYPE(t_atmos_for_ocean),  INTENT(inout)          :: p_as
    TYPE(t_atmos_fluxes ),    INTENT(inout)          :: atmos_fluxes
    TYPE(t_ocean_surface),    INTENT(inout)          :: p_oce_sfc
    TYPE(t_sea_ice),          INTENT(inout)          :: p_ice
    TYPE(t_hamocc_state),          INTENT(inout)      ::hamocc_state
    TYPE(t_operator_coeff), INTENT(in) :: operators_coefficients
    
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
        CALL write_initial_ocean_timestep(patch_3D,p_os(n_dom),p_oce_sfc,p_ice, operators_coefficients)
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

      CALL update_ocean_surface_refactor(patch_3D, p_os(n_dom), p_as, p_ice, atmos_fluxes, p_oce_sfc, &
             & this_datetime, operators_coefficients)

      !------------------------------------------------------------------------
      ! output: TODO not working for 3d prognostics
      ! add noise {{{     
                                                                                                
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
!TODO     CALL calc_fast_oce_diagnostics( patch_2d,      &
!TODO       & patch_3d%p_patch_1d(1)%dolic_c, &
!TODO       & patch_3d%p_patch_1d(1)%prism_thick_c, &
!TODO       & patch_3d%p_patch_1d(1)%zlev_m, &
!TODO       & p_os(n_dom)%p_diag)
!TODO     ! }}}
      CALL output_ocean( patch_3D,   &
        &                p_os(n_dom),&
        &                this_datetime,   &
        &                p_oce_sfc,  &
        &                p_ice,      &
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
    & this_datetime, p_oce_sfc, physics_parameters, &
    & p_as, p_atm_f, sea_ice, &
    & hamocc_state,operators_coefficients,solvercoeff_sp)
    
    TYPE(t_patch_3d ), POINTER, INTENT(in)          :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET, INTENT(inout) :: ocean_state(n_dom)
    TYPE(t_external_data), TARGET, INTENT(in)        :: p_ext_data(n_dom)
    TYPE(datetime), POINTER                          :: this_datetime
    TYPE(t_ocean_surface)                            :: p_oce_sfc
    TYPE(t_ho_params)                                :: physics_parameters
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

      IF (jstep == jstep0 ) THEN
        CALL update_ocean_surface_refactor(patch_3D, ocean_state(n_dom), p_as, sea_ice, p_atm_f, p_oce_sfc, &
             & this_datetime, operators_coefficients)
      ENDIF

      if (jstep == jstep0 ) then
        CALL update_height_depdendent_variables( patch_3d, ocean_state(jg), p_ext_data(jg), operators_coefficients, &
          & solvercoeff_sp)
      endif

      !---------------------------------------------------------------------

      if (jstep == jstep0 ) &
        & CALL update_ho_params(patch_3d, ocean_state(jg), p_as%fu10, sea_ice%concsum, physics_parameters, operators_coefficients)

      !------------------------------------------------------------------------
      ! solve for new free surface
      !CALL solve_free_surface_eq_ab (patch_3d, ocean_state(jg), p_ext_data(jg), &
      !  & ocean_surface, physics_parameters, jstep, operators_coefficients, solvercoeff_sp, return_status)!, p_int(jg))
      !IF (return_status /= 0) THEN
      ! CALL output_ocean(              &
      !   & patch_3d=patch_3d,          &
      !   & ocean_state=ocean_state,    &
      !   & this_datetime=mtime_current, &
      !   & ocean_surface=ocean_surface, &
      !   & sea_ice=sea_ice,            &
      !   & jstep=jstep, jstep0=jstep0, &
      !   & force_output=.true.)
      !  CALL finish(TRIM(routine), 'solve_free_surface_eq_ab  returned error')
      !ENDIF
      

      !------------------------------------------------------------------------
      ! Step 4: calculate final normal velocity from predicted horizontal
      ! velocity vn_pred and updated surface height
      !CALL calc_normal_velocity_ab(patch_3d, ocean_state(jg),&
      ! & operators_coefficients, solvercoeff_sp,  p_ext_data(jg), physics_parameters)

      !------------------------------------------------------------------------
      ! Step 5: calculate vertical velocity from continuity equation under
      ! incompressiblity condition in the non-shallow-water case
      !CALL calc_vert_velocity( patch_3d, ocean_state(jg),operators_coefficients)
      !------------------------------------------------------------------------

      !------------------------------------------------------------------------
      ! Step 6: transport tracers and diffuse them
      !IF (no_tracer>=1) THEN
      !  CALL advect_ocean_tracers( patch_3d, ocean_state(jg), physics_parameters,&
      !    & ocean_surface,&
      !    & operators_coefficients,&
      !    & jstep)
      !ENDIF

      ! perform accumulation for special variables

      ! update accumulated vars

      CALL output_ocean( patch_3d, ocean_state, &
        &                mtime_current,              &
        &                p_oce_sfc,             &
        &                sea_ice,                 &
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


  !-------------------------------------------------------------------------
  !>
  SUBROUTINE test_surface_flux_slo( patch_3d, p_os, &
    & this_datetime, p_oce_sfc,        &
    & p_as, atmos_fluxes, p_ice, operators_coefficients)
    
    TYPE(t_patch_3d ), POINTER, INTENT(in)          :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET, INTENT(inout) :: p_os(n_dom)
    TYPE(datetime),           POINTER                :: this_datetime
    TYPE(t_ocean_surface),    INTENT(inout)          :: p_oce_sfc
    TYPE(t_atmos_for_ocean),  INTENT(inout)          :: p_as
    TYPE(t_atmos_fluxes ),    INTENT(inout)          :: atmos_fluxes
    TYPE(t_sea_ice),          INTENT(inout)          :: p_ice
    TYPE(t_operator_coeff),   INTENT(in)          :: operators_coefficients
    
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
      & method_name = 'mo_ocean_testbed_modules:test_surface_flux_slo'
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
      saltBefore  = salt_content_in_surface(patch_2D, flat(:,:),p_ice, p_os(n_dom),p_oce_sfc, &
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
      CALL update_ocean_surface_refactor(patch_3D, p_os(n_dom), p_as, p_ice, atmos_fluxes, p_oce_sfc, &
           & this_datetime, operators_coefficients)
      !-----------------------------------------------------------------------------------------------------------------

      !---  energy  -----------------------------------------------------------
      energyCheck = energy_content_in_surface(patch_2d, flat(:,:), p_os(n_dom)%p_prog(nold(1))%h(:,:), &
        &             p_ice, sst(:,:), computation_type=budget_type_energy, info='UpdFlux')

      !---------DEBUG DIAGNOSTICS-------------------------------------------
      CALL dbg_print('TB.SfcFlux: energy  UpdFl' ,energyCheck    (:,:),debug_string, 2, in_subset=patch_2D%cells%owned)
      CALL dbg_print('TB.SfcFlux: sst AFT UpdFl' ,sst            (:,:),debug_string, 2, in_subset=patch_2D%cells%owned)
      CALL dbg_print('TB.SfcFlux: trc AFT UpdFl' ,p_os(n_dom)%p_prog(nold(1))%tracer(:,1,:,1), &
        &  debug_string, 5, in_subset=patch_2D%cells%owned)
      CALL dbg_print('TB.SfcFlux: Flx AFT UpdFl' ,p_oce_sfc%TopBC_Temp_vdiff(:,:), &
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
    !         & (dtime/delta_z ) * p_oce_sfc%TopBC_Temp_vdiff(jc,jb)
    !       sst(jc,jb) = sst(jc,jb) + (dtime/delta_z ) * p_oce_sfc%TopBC_Temp_vdiff(jc,jb)
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
      saltAfter   = salt_content_in_surface(patch_2D, flat(:,:),p_ice, p_os(n_dom),p_oce_sfc, &
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

      
      ! update accumulated vars
!     CALL update_ocean_statistics(p_os(n_dom), &
!       & p_oce_sfc,                       &
!       & patch_2D%cells%owned,                 &
!       & patch_2D%edges%owned,                 &
!       & patch_2D%verts%owned,                 &
!       & n_zlev)

      CALL output_ocean( patch_3D,   &
        &                p_os(n_dom),&
        &                this_datetime,   &
        &                p_oce_sfc,  &
        &                p_ice,      &
        &                jstep, jstep0)

      CALL update_time_indices(n_dom)

    END DO

    CALL timer_stop(timer_total)
    
  END SUBROUTINE test_surface_flux_slo
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  SUBROUTINE test_sea_ice( patch_3d, p_os, &
    & this_datetime, p_oce_sfc,        &
    & p_as, atmos_fluxes, p_ice, operators_coefficients)

    TYPE(t_patch_3d ), POINTER, INTENT(IN) :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET, INTENT(inout) :: p_os(n_dom)
    TYPE(t_ocean_surface),    INTENT(inout)          :: p_oce_sfc
    TYPE(t_atmos_for_ocean),  INTENT(inout)          :: p_as
    TYPE(t_atmos_fluxes ),    INTENT(inout)          :: atmos_fluxes
    TYPE(t_sea_ice),          INTENT(inout)          :: p_ice
    TYPE(t_operator_coeff),   INTENT(IN) :: operators_coefficients
    TYPE(datetime), POINTER                          :: this_datetime

    TYPE (t_hamocc_state)        :: hamocc_State

    ! local variables
    INTEGER                         :: jstep0, jstep, jg
    INTEGER                         :: budget_type_salt, budget_type_energy
    CHARACTER(LEN=32)               :: datestring
    REAL(wp), POINTER               :: sss(:,:), sst(:,:), ssh(:,:)!, area(:,:)
    REAL(wp), DIMENSION(nproma,patch_3D%p_patch_2D(1)%alloc_cell_blocks) &
      &                             :: energyBefore, energyAfter, saltBefore, saltAfter, energyDiff, saltDiff

    TYPE(t_patch), POINTER          :: p_patch
    TYPE(t_subset_range), POINTER   :: owned_cells
    TYPE(t_RestartAttributeList), POINTER :: restartAttributes
    CLASS(t_RestartDescriptor), POINTER :: restartDescriptor
    CHARACTER(LEN = *), PARAMETER   :: routine = 'mo_ocean_testbed_modules:test_sea_ice'

    TYPE(eventGroup), POINTER           :: checkpointEventGroup => NULL()
    TYPE(timedelta), POINTER            :: model_time_step => NULL()

    TYPE(datetime), POINTER             :: mtime_current     => NULL()
    TYPE(datetime), POINTER             :: eventRefDate      => NULL(), &
         &                                 eventStartDate    => NULL(), &
         &                                 eventEndDate      => NULL()
    TYPE(datetime), POINTER             :: checkpointRefDate => NULL(), &
         &                                 restartRefDate    => NULL()

    TYPE(timedelta), POINTER            :: eventInterval   => NULL()
    TYPE(event), POINTER                :: checkpointEvent => NULL()
    TYPE(event), POINTER                :: restartEvent    => NULL()

    INTEGER                             :: checkpointEvents, ierr
    LOGICAL                             :: lwrite_checkpoint, lret

    CHARACTER(LEN=MAX_DATETIME_STR_LEN)    :: dstring
    CHARACTER(len=MAX_MTIME_ERROR_STR_LEN) :: errstring

    LOGICAL :: l_isStartdate, l_isExpStopdate, l_isRestart, l_isCheckpoint, l_doWriteRestart

    !------------------------------------------------------------------
    ! no grid refinement allowed here so far
    !------------------------------------------------------------------
    IF (n_dom > 1 ) THEN
      CALL finish(TRIM(routine), ' N_DOM > 1 is not allowed')
    END IF
    jg = n_dom
    !------------------------------------------------------------------

    p_patch         => patch_3d%p_patch_2d(jg)
    owned_cells     => p_patch%cells%owned
!    flat            => patch_3D%p_patch_1D(1)%prism_thick_flat_sfc_c(:,1,:)

    energyBefore(:,:) = 0.0_wp
    saltBefore  (:,:) = 0.0_wp
    energyAfter (:,:) = 0.0_wp
    saltAfter   (:,:) = 0.0_wp
    energyDiff  (:,:) = 0.0_wp
    saltDiff    (:,:) = 0.0_wp

    budget_type_salt   = 5  ! #slo# merging ocean_sea-ice-thermodyn r208xx
    budget_type_energy = 0

    !------------------------------------------------------------------



    !------------------------------------------------------------------
    jstep0 = 0

    restartAttributes => getAttributesForRestarting()
    IF (ASSOCIATED(restartAttributes)) THEN
      ! get start counter for time loop from restart file:
      jstep0 = restartAttributes%getInteger("jstep")
    END IF
    IF (isRestart() .AND. mod(nold(jg),2) /=1 ) THEN
      ! swap the g_n and g_nm1
      CALL update_time_g_n(p_os(jg))
    ENDIF

    restartDescriptor => createRestartDescriptor("oce")

    !-------------------------- MTIME setup ---------------------------

    ! set events, group and the events
    CALL message('','')

    eventRefDate   => time_config%tc_exp_refdate
    eventStartDate => time_config%tc_exp_startdate
    eventEndDate   => time_config%tc_exp_stopdate

    ! for debugging purposes the referenece (anchor) date for checkpoint
    ! and restart may be switched to be relative to current jobs start
    ! date instead of the experiments start date.

    IF (time_config%is_relative_time) THEN
      checkpointRefDate => time_config%tc_startdate
      restartRefDate    => time_config%tc_startdate
    ELSE
      checkpointRefDate => time_config%tc_exp_startdate
      restartRefDate    => time_config%tc_exp_startdate
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
    ! start the time loop
    !------------------------------------------------------------------
    CALL timer_start(timer_total)

    jstep = jstep0
    TIME_LOOP: DO

      jstep = jstep + 1
      ! update model date and time mtime based
      mtime_current = mtime_current + model_time_step

      CALL datetimeToString(mtime_current, datestring)
      WRITE(message_text,'(a,i10,2a)') '  Begin of timestep =',jstep,'  datetime:  ', datestring
      CALL message (TRIM(routine), message_text)

      ! update pointer every timestep due to changing nold
      sst   => p_os(n_dom)%p_prog(nold(1))%tracer(:,1,:,1)
      sss   => p_os(n_dom)%p_prog(nold(1))%tracer(:,1,:,2)
      ssh   => p_os(n_dom)%p_prog(nold(1))%h(:,:)


      !---------DEBUG DIAGNOSTICS-------------------------------------------
      energyBefore = energy_in_surface(p_patch, p_ice, ssh(:,:), sst(:,:), computation_type=budget_type_energy, &
      &                               dbg_lev=2, info='test_sea_ice: before surface call')

      saltBefore  = salt_in_surface(p_patch, p_ice, sss(:,:), computation_type=budget_type_salt, &
        &                           dbg_lev=2, info='test_sea_ice: before surface call')
      !---------------------------------------------------------------------

      CALL update_ocean_surface_refactor(patch_3D, p_os(n_dom), p_as, p_ice, atmos_fluxes, p_oce_sfc, &
           & this_datetime, operators_coefficients)

      !---------------------------------------------------------------------

      !---------DEBUG DIAGNOSTICS-------------------------------------------
      energyAfter = energy_in_surface(p_patch, p_ice, ssh(:,:), sst(:,:), computation_type=budget_type_energy, &
      &                               dbg_lev=2, info='test_sea_ice: after surface call')

      saltAfter  = salt_in_surface(p_patch, p_ice, sss(:,:), computation_type=budget_type_salt, &
        &                           dbg_lev=2, info='test_sea_ice: after surface call')
      !---------------------------------------------------------------------

!      energyDiff = energyAfter - energyBefore
      saltDiff = saltAfter - saltBefore   ! this discribes the saltbudget in kg, which has to be zero
      CALL dbg_print('TB.SfcFlux: saltDiff', saltDiff(:,:), debug_string, 2, in_subset=owned_cells)

      ! hack for writing energy/salt budgets on p_ice%u/v:
      !  p_ice%u(:,:) = energyDits(:,:)
      !  p_ice%v(:,:) = saltBudget(:,:)
      !---------------------------------------------------------------------

!      ! update accumulated vars
!      CALL update_ocean_statistics(p_os(jg),p_oce_sfc, owned_cells, p_patch%edges%owned,&
!        &                                              p_patch%verts%owned, n_zlev)
!
      CALL output_ocean( patch_3D, p_os(jg), mtime_current, p_oce_sfc,  &
        &                p_ice, jstep, jstep0)

      ! Shift time indices for the next loop
      ! this HAS to ge into the restart files, because the start with the following loop
      CALL update_time_indices(jg)

      ! update intermediate timestepping variables for the tracers
      CALL update_time_g_n(p_os(jg))

      ! check whether time has come for writing restart file
      ! default is to assume we do not write a checkpoint/restart file
      lwrite_checkpoint = .FALSE.
      ! if thwe model is not supposed to write output, do not write checkpoints
      IF (.NOT. output_mode%l_none ) THEN
        ! to clarify the decision tree we use shorter and more expressive names:

        l_isStartdate    = (time_config%tc_startdate == mtime_current)
        l_isExpStopdate  = (time_config%tc_exp_stopdate == mtime_current)
        l_isRestart      = isCurrentEventActive(restartEvent, mtime_current)
        l_isCheckpoint   = isCurrentEventActive(checkpointEvent, mtime_current)
        l_doWriteRestart = time_config%tc_write_restart

        IF ( &
             !  if normal checkpoint or restart cycle has been reached, i.e. checkpoint+model stop
             &         (l_isRestart .OR. l_isCheckpoint)                     &
             &  .AND.                                                        &
             !  and the current date differs from the start date
             &        .NOT. l_isStartdate                                    &
             &  .AND.                                                        &
             !  and end of run has not been reached or restart writing has been disabled
             &        (.NOT. l_isExpStopdate .OR. l_doWriteRestart)          &
             & ) THEN
          lwrite_checkpoint = .TRUE.
        END IF
      END IF

      IF (lwrite_checkpoint) THEN
          CALL restartDescriptor%updatePatch(p_patch, &
                                            &opt_nice_class=1, &
                                            &opt_ocean_zlevels=n_zlev, &
                                            &opt_ocean_zheight_cellmiddle = patch_3d%p_patch_1d(1)%zlev_m(:), &
                                            &opt_ocean_zheight_cellinterfaces = patch_3d%p_patch_1d(1)%zlev_i(:))
          CALL restartDescriptor%writeRestart(mtime_current, jstep)
      END IF

      IF (mtime_current >= time_config%tc_stopdate) THEN
        ! leave time loop
        EXIT TIME_LOOP
      END IF

    ENDDO TIME_LOOP

    CALL timer_stop(timer_total)
  

  END SUBROUTINE test_sea_ice
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

!    CALL calc_neutralslope_coeff( &
!      &    patch_3d,              &
!      &    p_os(n_dom)%p_prog(nold(1))%tracer(:,:,:,:), &
!!       &    p_os(n_dom)%p_prog(nold(1))%h(:,:), &
!      &    alph, beta)

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

  !-------------------------------------------------------------------------
  !>
  SUBROUTINE checkVarlistsForOutput(patch_3d, ocean_state, &
    & this_datetime, physics_parameters, &
    & p_as, atmos_fluxes, p_oce_sfc, p_ice, hamocc_state,operators_coefficients)

    TYPE(t_patch_3d), POINTER, INTENT(in) :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET, INTENT(inout) :: ocean_state
    TYPE(datetime), POINTER                          :: this_datetime
    TYPE (t_ho_params)                               :: physics_parameters
    TYPE(t_atmos_for_ocean),  INTENT(inout)          :: p_as
    TYPE(t_atmos_fluxes ),    INTENT(inout)          :: atmos_fluxes
    TYPE(t_ocean_surface),    INTENT(inout)          :: p_oce_sfc
    TYPE(t_sea_ice),          INTENT(inout)          :: p_ice
    TYPE(t_hamocc_state),          INTENT(inout)      ::hamocc_state
    TYPE(t_operator_coeff), INTENT(in) :: operators_coefficients

    IF (output_mode%l_nml) THEN
      CALL write_initial_ocean_timestep(patch_3d,ocean_state, &
          &  p_oce_sfc,p_ice, operators_coefficients)
    ENDIF
  END SUBROUTINE checkVarlistsForOutput
  
  SUBROUTINE checkVarlistKeys(patch_2d)
    TYPE(t_patch), TARGET, INTENT(in) :: patch_2d
    
    CHARACTER(LEN=max_char_length) :: listname
    TYPE(t_var_list)     :: varnameCheckList
    integer :: alloc_cell_blocks

    REAL(wp), POINTER :: var0(:,:,:)
    REAL(wp), POINTER :: var1(:,:,:)
    REAL(wp), POINTER :: var2(:,:,:)
    REAL(wp), POINTER :: var3(:,:,:)
    
    WRITE(listname,'(a)')  'varnameCheck_list'
    CALL new_var_list(varnameCheckList, listname, patch_id=patch_2d%id)
    CALL default_var_list_settings( varnameCheckList,  &
      & lrestart=.TRUE.,loutput=.TRUE.,&
      & model_type='oce' )

    alloc_cell_blocks = patch_2d%alloc_cell_blocks
    call add_var(varnamechecklist,'h',var0,grid_unstructured_cell, za_depth_below_sea_half, &
        & t_cf_var('h','m','h',datatype_flt64,'ssh'),&
        & grib2_var(255, 255, 255, datatype_pack16, grid_unstructured, grid_cell),&
        & ldims=(/nproma, n_zlev, alloc_cell_blocks/))
    call add_var(varnamechecklist,'H',var0,grid_unstructured_cell, za_depth_below_sea_half, &
        & t_cf_var('h','m','h',datatype_flt64,'ssh'),&
        & grib2_var(255, 255, 255, datatype_pack16, grid_unstructured, grid_cell),&
        & ldims=(/nproma, n_zlev, alloc_cell_blocks/))
    call add_var(varnameCheckList,'t',var1,grid_unstructured_cell, za_depth_below_sea_half, &
        & t_cf_var('t','m','t',DATATYPE_FLT64,'t'),&
        & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL),&
        & ldims=(/nproma, n_zlev, alloc_cell_blocks/))
    call add_var(varnameCheckList,'s',var2,grid_unstructured_cell, za_depth_below_sea_half, &
        & t_cf_var('s','m','s',DATATYPE_FLT64,'s'),&
        & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL),&
        & ldims=(/nproma, n_zlev, alloc_cell_blocks/))
    call add_var(varnameCheckList,'sS',var2,grid_unstructured_cell, za_depth_below_sea_half, &
        & t_cf_var('s','m','s',DATATYPE_FLT64,'s'),&
        & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL),&
        & ldims=(/nproma, n_zlev, alloc_cell_blocks/))
    call add_var(varnameCheckList,'ss',var2,grid_unstructured_cell, za_depth_below_sea_half, &
        & t_cf_var('s','m','s',DATATYPE_FLT64,'s'),&
        & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL),&
        & ldims=(/nproma, n_zlev, alloc_cell_blocks/))


    call print_var_list(varnameCheckList)
    call delete_var_list(varnameCheckList)
  END SUBROUTINE checkVarlistKeys
 
  
  !>
  !! Testing for operators in the zstar transformed co-ordinates 
  !!
  !!
  SUBROUTINE test_zstar( patch_3d, ocean_state, p_ext_data,  &
    & this_datetime, p_oce_sfc, physics_parameters, &
    & p_as, p_atm_f, sea_ice, &
    & hamocc_state,operators_coefficients,solvercoeff_sp)
    
    TYPE(t_patch_3d ), POINTER, INTENT(in)          :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET, INTENT(inout) :: ocean_state(n_dom)
    TYPE(t_external_data), TARGET, INTENT(in)        :: p_ext_data(n_dom)
    TYPE(datetime), POINTER                          :: this_datetime
    TYPE(t_ocean_surface)                            :: p_oce_sfc
    TYPE(t_ho_params)                                :: physics_parameters
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
   
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
    !! New variables to test stretching 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
    REAL(wp) :: vn_test(nproma, n_zlev, patch_3d%p_patch_2d(1)%nblks_e) 
    TYPE(t_subset_range), POINTER :: all_edges, edges_in_domain
    TYPE(t_subset_range), POINTER :: all_cells
    REAL(wp) :: point_lon, point_lat     ! latitude of point
    REAL(wp) :: uu, vv      ! zonal,  meridional velocity
    REAL(wp) :: angle1, angle2, edge_vn, COS_angle1, SIN_angle1
    REAL(wp) :: u0, sphere_radius 
    INTEGER  :: edge_index, edge_block, start_edges_index, end_edges_index
    INTEGER  :: blockNo, cell_index, cell_block, start_cells_index, end_cells_index
    TYPE(t_cartesian_coordinates) :: cell_center, edge_center
    INTEGER  :: start_index, end_index, neigbor 
    INTEGER  :: bot_level 
    INTEGER  :: zcoord_type 
    INTEGER,  DIMENSION(:,:,:), POINTER :: idx, blk
    INTEGER  :: id1, id2, bl1, bl2 
    INTEGER  :: bt_level 
    REAL(wp) :: ht_edge, ht1, ht2, ht3 
    REAL(wp) :: pt_l1, pt_l2 
    REAL(wp) :: z1, z2, H1, H2, eta1, eta2 
    REAL(wp) :: eta(nproma, patch_3d%p_patch_2d(1)%alloc_cell_blocks) 
    REAL(wp) :: mom_grad(nproma, n_zlev, patch_3d%p_patch_2d(1)%nblks_e) 
    REAL(wp) :: wstar(nproma, n_zlev + 1, patch_3d%p_patch_2d(1)%alloc_cell_blocks) 
    REAL(wp) :: zgrad(nproma, n_zlev, patch_3d%p_patch_2d(1)%nblks_e) 
    TYPE(t_cartesian_coordinates) :: p_zgrad(nproma, n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks) 
    REAL(wp) :: p_zgrad_sc(nproma, n_zlev ,patch_3d%p_patch_2d(1)%alloc_cell_blocks) 
    REAL(wp) :: p_zgrad_top(nproma, n_zlev + 1, patch_3d%p_patch_2d(1)%alloc_cell_blocks) 
    REAL(wp) :: phi(nproma, n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks) 
    REAL(wp) :: dphi_dz1, dphi_dz2, dz_dn, dphi_dn, dpv_dn, ke1, ke2, dk_dn
    REAL(wp) :: pv_x1, pv_x2, pv_y1, pv_y2, dpvx_dn, dpvy_dn
    
    TYPE(t_cartesian_coordinates) :: z_adv_u_i(nproma, n_zlev+1)
    TYPE(t_cartesian_coordinates) :: z_adv_u_m(nproma, n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    TYPE(t_cartesian_coordinates) :: pvn(nproma, n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks) 
    REAL(wp) :: vn_z(1:nproma,1:n_zlev,patch_3D%p_patch_2D(1)%nblks_e)
    REAL(wp) :: ptpvn(1:nproma,1:n_zlev,patch_3D%p_patch_2D(1)%nblks_e)
    REAL(wp) :: vert_der_e(1:nproma,1:n_zlev,patch_3D%p_patch_2D(1)%nblks_e)
    REAL(wp) :: div_v(nproma, n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks) 
    REAL(wp) :: div_v_z(nproma, n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks) 
    REAL(wp) :: coeff, dz_dy
 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

    TYPE(t_subset_range), POINTER :: owned_edges, owned_cells
    !------------------------------------------------------------------
    patch_2D      => patch_3d%p_patch_2d(1)
    idx           => patch_3D%p_patch_2D(1)%edges%cell_idx
    blk           => patch_3D%p_patch_2D(1)%edges%cell_blk
    edges_in_domain => patch_3D%p_patch_2D(1)%edges%in_domain

    IF (n_dom > 1 ) THEN
      CALL finish(TRIM(routine), ' N_DOM > 1 is not allowed')
    END IF
    jg = n_dom

    vn_test = 0.0 

    all_edges   => patch_2d%edges%ALL
    all_cells   => patch_2d%cells%ALL
    owned_cells => patch_2D%cells%owned
    
    sphere_radius = grid_sphere_radius
    u0            = (2.0_wp*pi*sphere_radius)/(12.0_wp*24.0_wp*3600.0_wp)

    bot_level = n_zlev + 1

    !! If zcoord type is 1, then we set eta equal to
    !! sea surface height type 201 which is what we compare against
    eta = 0.  
    ! Initialize eta for zstar
    ! #slo#: simple elevation between 30W and 30E (pi/3.)
    DO cell_block = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, cell_block, start_index, end_index)
      DO cell_index = start_index, end_index
!        IF ( patch_3d%lsm_c(cell_index,1,cell_block) <= sea_boundary ) THEN
    
!           eta(cell_index, cell_block) = 10.0_wp * &
!            & SIN(patch_2d%cells%center(cell_index, cell_block)%lon * 6.0_wp) &
!            & * COS(patch_2d%cells%center(cell_index, cell_block)%lat * 3.0_wp)
               
           eta(cell_index, cell_block) = 2.0_wp*patch_2d%cells%center(cell_index, cell_block)%lat 

!        ENDIF
      END DO
    END DO
 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! 0 for z, 1 for zstar
    !! Make sure that the height is not initialized if zcoord_type = 1
    zcoord_type = 1

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! FIXME: For now we assume that the height of the edges is not 
    !! available as a separate array
    !! Initialize velocity
    !! Note that the velocity only has a z dependence. However on a zstar grid
    !! with arbitrary order variation in z, derivatives like dv/dn will not
    !! be exact for constant in x, y and linear in z initialization
    DO edge_block = all_edges%start_block, all_edges%end_block
      CALL get_index_range(all_edges, edge_block, start_edges_index, end_edges_index)
      DO edge_index = start_edges_index, end_edges_index
        id1 = (idx(edge_index, edge_block, 1))
        id2 = (idx(edge_index, edge_block, 2))
        bl1 = (blk(edge_index, edge_block, 1))
        bl2 = (blk(edge_index, edge_block, 2))


        point_lon = patch_2d%edges%center(edge_index,edge_block)%lon  
        point_lat = patch_2d%edges%center(edge_index,edge_block)%lat  
        pt_l1     = patch_2d%cells%center(id1,bl1)%lat  
        pt_l2     = patch_2d%cells%center(id2,bl2)%lat  

        bt_level = patch_3d%p_patch_1d(1)%dolic_e(edge_index,edge_block)      
        DO level = 1, bt_level 
          ht_edge = 0.
          IF (zcoord_type == 0) THEN
              !! Get height of the cell center at mid point from the bottom
              ht1 = patch_3d%p_patch_1d(1)%depth_CellInterface(id1, bt_level + 1, bl1)& 
                  & - patch_3d%p_patch_1d(1)%depth_CellMiddle(id1, level, bl1) 
              ht2 = patch_3d%p_patch_1d(1)%depth_CellInterface(id2, bt_level + 1, bl2)& 
                  & - patch_3d%p_patch_1d(1)%depth_CellMiddle(id2, level, bl2) 
    
              !! Get height of edge as average of the centers of the 2 adjoining cells
              ht_edge = 0.5_wp*( ht1 + ht2 )
          ELSE IF (zcoord_type == 1) THEN
              !! Get height of the cell center at mid point from the bottom
              H1  = patch_3d%p_patch_1d(1)%depth_CellInterface(id1,bt_level+1,bl1)
              ht1 = H1 - patch_3d%p_patch_1d(1)%depth_CellMiddle(id1, level, bl1) 
              H2  = patch_3d%p_patch_1d(1)%depth_CellInterface(id2,bt_level+1,bl2) 
              ht2 = H2 - patch_3d%p_patch_1d(1)%depth_CellMiddle(id2, level, bl2) 
              
              eta1 = eta(id1, bl1) 
              eta2 = eta(id2, bl2) 

              !! Transform to z from zstar
              z1   = ht1*(H1 + eta1)/H1 + eta1
              z2   = ht2*(H2 + eta2)/H2 + eta2
   
              !! Get height of edge as average of the centers of the 2 adjoining cells
              ht_edge = 0.5_wp*( z1 + z2 )

          ENDIF

          uu = 0.0_wp 
          vv = 1.0_wp + 0.05_wp*ht_edge

          edge_vn = uu * patch_2d%edges%primal_normal(edge_index,edge_block)%v1 &
            & + vv * patch_2d%edges%primal_normal(edge_index,edge_block)%v2

          ocean_state(jg)%p_prog(nold(1))%vn(edge_index, level, edge_block) &
              & = edge_vn
          
          vn_z(edge_index, level, edge_block) =  edge_vn*ht_edge
          
        ENDDO

      ENDDO
    ENDDO


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!Do tests on functions and operators
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    !! Show that cell centers and edge centers are not in the same z plane
    DO cell_block = owned_cells%start_block, owned_cells%end_block
      CALL get_index_range(owned_cells, cell_block, start_index, end_index)
      DO cell_index = start_index, end_index

        cell_center%x = patch_2D%cells%cartesian_center(cell_index, cell_block)%x

        !-------------------------------
        DO neigbor=1, patch_2D%cells%num_edges(cell_index,cell_block)!no_primal_edges

          edge_index = patch_2D%cells%edge_idx(cell_index, cell_block, neigbor)
          edge_block = patch_2D%cells%edge_blk(cell_index, cell_block, neigbor)

          IF (edge_block > 0 ) THEN
!            write(*, *)  patch_2D%edges%cartesian_center(edge_index,edge_block)%x(1), &
!              &  patch_2D%edges%cartesian_center(edge_index,edge_block)%x(2), &
!              &  patch_2D%edges%cartesian_center(edge_index,edge_block)%x(3), &
!              & cell_center%x(1), cell_center%x(2), cell_center%x(3) 

          ENDIF !(edge_block > 0 )
        ENDDO !neigbor=1,patch_2D%num_edges
        !-------------------------------
      ENDDO ! cell_index = start_index, end_index
    ENDDO !cell_block = owned_cells%start_block, owned_cells%end_block

    
    !! Show that Pv is not constant for a constant input velocity
    DO blockNo = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, blockNo, start_cells_index, end_cells_index)
        DO cell_index =  start_cells_index, end_cells_index
          DO level = 1, patch_3d%p_patch_1d(1)%dolic_c(cell_index,blockNo)      

!              write(*, *)ocean_state(jg)%p_diag%p_vn(cell_index,level,blockNo)%x
              
          ENDDO
        ENDDO
    ENDDO


    !! Initialize a function that is constant in z plane and 
    !! linear in height
    phi = 0.  
    DO blockNo = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, cell_block, start_index, end_index)
      DO cell_index = start_index, end_index
 
          bt_level = patch_3d%p_patch_1d(1)%dolic_c(cell_index,blockNo)      
          H1       = patch_3d%p_patch_1d(1)%depth_CellInterface(cell_index, bt_level+1, blockNo)

          DO level = 1, patch_3d%p_patch_1d(1)%dolic_c(cell_index,blockNo)      
   
              ht1  = H1 - patch_3d%p_patch_1d(1)%depth_CellMiddle(cell_index, level, blockNo) 
              eta1 = eta(cell_index, blockNo) 
              z1   = ht1*(H1 + eta1)/H1 + eta1

              phi(cell_index, level, blockNo) = 0.5_wp + 0.1_wp*z1

          END DO

      END DO
    END DO

    !! Show that the chain rule d/dn = dz/dn*ds/dz*d/ds works for 
    !! functions with linear z dependence
    DO edge_block = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(all_edges, edge_block, start_edges_index, end_edges_index)
      DO edge_index = start_edges_index, end_edges_index
        bt_level = patch_3d%p_patch_1d(1)%dolic_e(edge_index,edge_block)      
        DO level = 2, bt_level - 1 ! FIXME: starts from second level
          id1 = (idx(edge_index, edge_block, 1))
          id2 = (idx(edge_index, edge_block, 2))
          bl1 = (blk(edge_index, edge_block, 1))
          bl2 = (blk(edge_index, edge_block, 2))
    
          !! Get height of the cell center at mid point from the bottom
          H1  = patch_3d%p_patch_1d(1)%depth_CellInterface(id1,bt_level+1,bl1)
          ht1 = H1 - patch_3d%p_patch_1d(1)%depth_CellMiddle(id1, level, bl1) 
          H2  = patch_3d%p_patch_1d(1)%depth_CellInterface(id2,bt_level+1,bl2) 
          ht2 = H2 - patch_3d%p_patch_1d(1)%depth_CellMiddle(id2, level, bl2) 
          
          eta1 = eta(id1, bl1) 
          eta2 = eta(id2, bl2) 
    
          !! Transform to z from zstar
          z1   = ht1*(H1 + eta1)/H1 + eta1
          z2   = ht2*(H2 + eta2)/H2 + eta2
   
          !! Note that ds/dz is the analytic value
          !! We are using the interface value instead of the middle value for the test 
          dphi_dz1 =  ( H1/(H1 + eta1)  ) *                                                       & 
             & ( -phi(id1, level, bl1) + phi(id1, level - 1, bl1) ) *                             &
             & patch_3D%p_patch_1D(1)%constantPrismCenters_invZdistance(id1, level, bl1)
          dphi_dz2 =  ( H2/(H2 + eta2)  ) *                                                       & 
             & ( -phi(id2, level, bl2) + phi(id2, level - 1, bl2) ) *                             &
             & patch_3D%p_patch_1D(1)%constantPrismCenters_invZdistance(id2, level, bl2)
    
          dz_dn   =  (z1 - z2)*operators_coefficients%grad_coeff(edge_index,level,edge_block)
          
          dphi_dn =  (phi(id1,level,bl1) - phi(id2,level,bl2))*                                 &
              & operators_coefficients%grad_coeff(edge_index,level,edge_block)
    
!          write(*, *) dz_dn, dphi_dz1, dphi_dn 
    
    
        ENDDO
    
      ENDDO
    ENDDO
    
    CALL map_edges2cell_3d(patch_3d, ocean_state(jg)%p_prog(nold(1))%vn, &
        & operators_coefficients, pvn) 
    CALL map_cell2edges_3D( patch_3D, pvn, ptpvn, operators_coefficients)

    CALL div_oce_3d( ocean_state(jg)%p_prog(nold(1))%vn, patch_3D, operators_coefficients%div_coeff, div_v)
    CALL div_oce_3d( vn_z, patch_3D, operators_coefficients%div_coeff, div_v_z)

    DO blockNo = all_cells%start_block, all_cells%end_block
      !vertical derivative at ocean interior Surface is handled below
      ! this does not include h
      CALL verticalDeriv_vec_midlevel_on_block( patch_3d, &
                                              & pvn(:,:,blockNo),  &
                                              & z_adv_u_i(:,:),&
                                              & 1+1,             & ! FIXME: starts from second level
                                              & blockNo, start_index, end_index)

      !! TODO: Need to multiply with ds/dz to get values at cell center
      !! Then need to multiply with edge values, which will require transformation
      CALL get_index_range(all_cells, cell_block, start_index, end_index)

      do cell_index = start_index, end_index
 
          bt_level = patch_3d%p_patch_1d(1)%dolic_c(cell_index,blockNo)      
          H1       = patch_3d%p_patch_1d(1)%depth_CellInterface(cell_index, bt_level+1, blockNo)
          eta1     = eta(cell_index, blockNo) 

          z_adv_u_i(cell_index, :)%x(1) = z_adv_u_i(cell_index, :)%x(1)*( H1/(H1 + eta1)  )
          z_adv_u_i(cell_index, :)%x(2) = z_adv_u_i(cell_index, :)%x(2)*( H1/(H1 + eta1)  )
          z_adv_u_i(cell_index, :)%x(3) = z_adv_u_i(cell_index, :)%x(3)*( H1/(H1 + eta1)  )

          DO level = 1, patch_3d%p_patch_1d(1)%dolic_c(cell_index,blockNo)      
          
              ht1   = H1 - patch_3d%p_patch_1d(1)%depth_CellMiddle(cell_index, level, blockNo) 
              z1    = ht1*(H1 + eta1)/H1 + eta1

              !! Multiply with u.grad(z) = div(uz) - z*div(u)
              coeff  = div_v_z(cell_index, level, blockNo) - z1*div_v(cell_index, level, blockNo)  

              dz_dy  = (2._wp*(1._wp + ht1/H1))/44106._wp

              !! Show that v.grad(z) is exact
!              write(*, *) level, coeff, (1._wp + 0.05_wp*z1)*dz_dy 
!              write(*, *) level, div_v_z(cell_index, level, blockNo), & 
!                  & dz_dy + 0.1*z1*dz_dy
!              write(*, *) level, div_v(cell_index, level, blockNo) , & 
!                  & 0.05_wp*dz_dy

              z_adv_u_i(cell_index, level)%x(1) = z_adv_u_i(cell_index, level)%x(1)*coeff  
              z_adv_u_i(cell_index, level)%x(2) = z_adv_u_i(cell_index, level)%x(2)*coeff
              z_adv_u_i(cell_index, level)%x(3) = z_adv_u_i(cell_index, level)%x(3)*coeff

          END DO

      END DO
     
      !! Project values to prism middle from top
      CALL map_vec_prismtop2center_on_block(patch_3d, z_adv_u_i, z_adv_u_m(:,:,blockNo), &
        & blockNo, start_index, end_index)

    END DO

    !! Get the derivatives at the cell edges    
    CALL map_cell2edges_3D( patch_3D, z_adv_u_m, vert_der_e, operators_coefficients)

    
    !! Get kinetic energy
    CALL calc_scalar_product_veloc_3d( patch_3d,  &
        & ocean_state(jg)%p_prog(nold(1))%vn,    &
        & ocean_state(jg)%p_diag,                     &
        & operators_coefficients)

    ! STEP 1: horizontal advection
    CALL veloc_adv_horz_mimetic( patch_3d,         &
      & ocean_state(jg)%p_prog(nold(1))%vn,    &
      & ocean_state(jg)%p_prog(nold(1))%vn,    &
      & ocean_state(jg)%p_diag,                &
      & ocean_state(jg)%p_diag%veloc_adv_horz, &
      & operators_coefficients)

    !! Horizontal total derivative
    mom_grad = ocean_state(jg)%p_diag%grad + ocean_state(jg)%p_diag%veloc_adv_horz
          
    DO edge_block = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(all_edges, edge_block, start_edges_index, end_edges_index)
      DO edge_index = start_edges_index, end_edges_index
        bt_level = patch_3d%p_patch_1d(1)%dolic_e(edge_index,edge_block)      
        DO level = 2, bt_level - 1 ! FIXME: starts from second level
          id1 = (idx(edge_index, edge_block, 1))
          id2 = (idx(edge_index, edge_block, 2))
          bl1 = (blk(edge_index, edge_block, 1))
          bl2 = (blk(edge_index, edge_block, 2))
    
          !! Get height of the cell center at mid point from the bottom
          H1  = patch_3d%p_patch_1d(1)%depth_CellInterface(id1,bt_level+1,bl1)
          ht1 = H1 - patch_3d%p_patch_1d(1)%depth_CellMiddle(id1, level, bl1) 
          H2  = patch_3d%p_patch_1d(1)%depth_CellInterface(id2,bt_level+1,bl2) 
          ht2 = H2 - patch_3d%p_patch_1d(1)%depth_CellMiddle(id2, level, bl2) 
          
          eta1 = eta(id1, bl1) 
          eta2 = eta(id2, bl2) 
    
          !! Transform to z from zstar
          z1   = ht1*(H1 + eta1)/H1 + eta1
          z2   = ht2*(H2 + eta2)/H2 + eta2
          
          !! Get height of edge as average of the centers of the 2 adjoining cells
          ht_edge = 0.5_wp*( z1 + z2 )

          pv_x1  = pvn(id1, level, bl1)%x(1) 
          pv_x2  = pvn(id2, level, bl2)%x(1)
          pv_y1  = pvn(id1, level, bl1)%x(2) 
          pv_y2  = pvn(id2, level, bl2)%x(2)

          ke1    = 0.5_wp*( pv_x1**2._wp + pv_y1**2._wp )
          ke2    = 0.5_wp*( pv_x2**2._wp + pv_y2**2._wp )
  
          dz_dn    =  (z1 - z2)*operators_coefficients%grad_coeff(edge_index,level,edge_block)
          dpvx_dn  =  (pv_x1 - pv_x2)*operators_coefficients%grad_coeff(edge_index,level,edge_block)
          dpvy_dn  =  (pv_y1 - pv_y2)*operators_coefficients%grad_coeff(edge_index,level,edge_block)
          
          dpv_dn   = dpvx_dn*patch_2d%edges%primal_normal(edge_index,edge_block)%v1 + &
              & dpvy_dn*patch_2d%edges%primal_normal(edge_index,edge_block)%v2 
          
          dk_dn    =  (ke1 - ke2)*operators_coefficients%grad_coeff(edge_index,level,edge_block)
              
!          !! Show that after vn and PTP(vn) are the same 
!          write(*, *) level, ocean_state(jg)%p_prog(nold(1))%vn(edge_index, level, edge_block), & 
!              & ptpvn(edge_index, level, edge_block) 

          !! Show that after reconstruction using P, averages are conserved but not difference 
!          write(*, *) level,  pv_y1 + pv_y2, 2._wp + 0.05_wp*( z1 + z2 )
!          write(*, *) level,  pv_y1 - pv_y2,         0.05_wp*( z1 - z2 )
!          write(*, *) level,  ( pv_y1 - pv_y2 )-0.05_wp*( z1 - z2 )

          !! Show that the derivative is exact only if normal is  
          !! aligned with the velocity
!          write(*, *) level, dpv_dn, 0.05_wp*dz_dn,  &
!              & patch_2d%edges%primal_normal(edge_index,edge_block)%v1, &
!              &  patch_2d%edges%primal_normal(edge_index,edge_block)%v2

!           !! Show that u.grad(u) NE ( u.grad(z) ) du/dz
!           !! NE -> not equal 
!           write(*, *) level, mom_grad(edge_index, level, edge_block), &
!              &  vert_der_e(edge_index, level, edge_block)

        ENDDO
    
      ENDDO
    ENDDO
    

  
  END SUBROUTINE test_zstar

 
  
  !>
  !! Testing for the advection of velocity in the momentum equation 
  !!
  !!
  SUBROUTINE test_mom( patch_3d, ocean_state, p_ext_data,  &
    & this_datetime, p_oce_sfc, physics_parameters, &
    & p_as, p_atm_f, sea_ice, &
    & hamocc_state,operators_coefficients,solvercoeff_sp)
    
    TYPE(t_patch_3d ), POINTER, INTENT(in)          :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET, INTENT(inout) :: ocean_state(n_dom)
    TYPE(t_external_data), TARGET, INTENT(in)        :: p_ext_data(n_dom)
    TYPE(datetime), POINTER                          :: this_datetime
    TYPE(t_ocean_surface)                            :: p_oce_sfc
    TYPE(t_ho_params)                                :: physics_parameters
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
   
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
    !! New variables to test stretching 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
    REAL(wp) :: vn_test(nproma, n_zlev, patch_3d%p_patch_2d(1)%nblks_e) 
    TYPE(t_subset_range), POINTER :: all_edges, edges_in_domain
    TYPE(t_subset_range), POINTER :: all_cells
    REAL(wp) :: point_lon, point_lat     ! latitude of point
    REAL(wp) :: uu, vv      ! zonal,  meridional velocity
    REAL(wp) :: angle1, angle2, edge_vn, COS_angle1, SIN_angle1
    REAL(wp) :: u0, sphere_radius 
    INTEGER  :: edge_index, edge_block, start_edges_index, end_edges_index
    INTEGER  :: blockNo, cell_index, cell_block, start_cells_index, end_cells_index
    TYPE(t_cartesian_coordinates) :: cell_center, edge_center
    INTEGER  :: start_index, end_index, neigbor 
    INTEGER  :: bot_level 
    INTEGER  :: zcoord_type 
    INTEGER,  DIMENSION(:,:,:), POINTER :: idx, blk
    INTEGER  :: id1, id2, bl1, bl2 
    INTEGER  :: bt_level 
    REAL(wp) :: ht_edge, ht1, ht2, ht3 
    REAL(wp) :: z1, z2, H1, H2, eta1, eta2 
    REAL(wp) :: eta(nproma, patch_3d%p_patch_2d(1)%alloc_cell_blocks) 
    REAL(wp) :: mom_grad(nproma, n_zlev, patch_3d%p_patch_2d(1)%nblks_e) 
    REAL(wp) :: wstar(nproma, n_zlev + 1, patch_3d%p_patch_2d(1)%alloc_cell_blocks) 
    REAL(wp) :: zgrad(nproma, n_zlev, patch_3d%p_patch_2d(1)%nblks_e) 
    TYPE(t_cartesian_coordinates) :: p_zgrad(nproma, n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks) 
    REAL(wp) :: p_zgrad_sc(nproma, n_zlev ,patch_3d%p_patch_2d(1)%alloc_cell_blocks) 
    REAL(wp) :: p_zgrad_top(nproma, n_zlev + 1, patch_3d%p_patch_2d(1)%alloc_cell_blocks) 
    REAL(wp) :: phi(nproma, n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks) 
    REAL(wp) :: dphi_dz1, dphi_dz2, dz_dn, dphi_dn
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

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
    
    TYPE(t_subset_range), POINTER :: owned_edges, owned_cells
    !------------------------------------------------------------------
    patch_2D      => patch_3d%p_patch_2d(1)
    idx           => patch_3D%p_patch_2D(1)%edges%cell_idx
    blk           => patch_3D%p_patch_2D(1)%edges%cell_blk
    edges_in_domain => patch_3D%p_patch_2D(1)%edges%in_domain

    IF (n_dom > 1 ) THEN
      CALL finish(TRIM(routine), ' N_DOM > 1 is not allowed')
    END IF
    jg = n_dom

    vn_test = 0.0 

    all_edges   => patch_2d%edges%ALL
    all_cells   => patch_2d%cells%ALL
    owned_cells => patch_2D%cells%owned
    
    sphere_radius = grid_sphere_radius
    u0            = (2.0_wp*pi*sphere_radius)/(12.0_wp*24.0_wp*3600.0_wp)

    bot_level = n_zlev + 1

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! 0 for z, 1 for zstar
    !! Make sure that the height is not initialized if zcoord_type = 1
    zcoord_type = 1 

    IF (zcoord_type == 1) THEN
        !! If zcoord type is 1, then we set eta equal to
        !! sea surface height type 201 which is what we compare against
        eta = 0.  
        ! Initialize eta for zstar
        ! #slo#: simple elevation between 30W and 30E (pi/3.)
        DO cell_block = all_cells%start_block, all_cells%end_block
          CALL get_index_range(all_cells, cell_block, start_index, end_index)
          DO cell_index = start_index, end_index
            IF ( patch_3d%lsm_c(cell_index,1,cell_block) <= sea_boundary ) THEN
    
               eta(cell_index, cell_block) = 10.0_wp * &
                & SIN(patch_2d%cells%center(cell_index, cell_block)%lon * 6.0_wp) &
                & * COS(patch_2d%cells%center(cell_index, cell_block)%lat * 3.0_wp)

            ENDIF
          END DO
        END DO
       
    ENDIF

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! FIXME: For now we assume that the height of the edges is not 
    !! available as a separate array
    !! Initialize velocity
    DO edge_block = all_edges%start_block, all_edges%end_block
      CALL get_index_range(all_edges, edge_block, start_edges_index, end_edges_index)
      DO edge_index = start_edges_index, end_edges_index
        point_lon = patch_2d%edges%center(edge_index,edge_block)%lon  
        point_lat = patch_2d%edges%center(edge_index,edge_block)%lat  

        bt_level = patch_3d%p_patch_1d(1)%dolic_e(edge_index,edge_block)      
        DO level = 1, bt_level 
          id1 = (idx(edge_index, edge_block, 1))
          id2 = (idx(edge_index, edge_block, 2))
          bl1 = (blk(edge_index, edge_block, 1))
          bl2 = (blk(edge_index, edge_block, 2))

          ht_edge = 0.
          IF (zcoord_type == 0) THEN
              !! Get height of the cell center at mid point from the bottom
              ht1 = patch_3d%p_patch_1d(1)%depth_CellInterface(id1, bt_level + 1, bl1)& 
                  & - patch_3d%p_patch_1d(1)%depth_CellMiddle(id1, level, bl1) 
              ht2 = patch_3d%p_patch_1d(1)%depth_CellInterface(id2, bt_level + 1, bl2)& 
                  & - patch_3d%p_patch_1d(1)%depth_CellMiddle(id2, level, bl2) 
    
              !! Get height of edge as average of the centers of the 2 adjoining cells
              ht_edge = 0.5*( ht1 + ht2 )
          ELSE IF (zcoord_type == 1) THEN
              !! Get height of the cell center at mid point from the bottom
              H1  = patch_3d%p_patch_1d(1)%depth_CellInterface(id1,bt_level+1,bl1)
              ht1 = H1 - patch_3d%p_patch_1d(1)%depth_CellMiddle(id1, level, bl1) 
              H2  = patch_3d%p_patch_1d(1)%depth_CellInterface(id2,bt_level+1,bl2) 
              ht2 = H2 - patch_3d%p_patch_1d(1)%depth_CellMiddle(id2, level, bl2) 
              
              eta1 = eta(id1, bl1) 
              eta2 = eta(id2, bl2) 

              !! Transform to z from zstar
              z1   = ht1*(H1 + eta1)/H1 + eta1
              z2   = ht2*(H2 + eta2)/H2 + eta2
    
              !! Get height of edge as average of the centers of the 2 adjoining cells
              ht_edge = 0.5*( z1 + z2 )

          ENDIF

          uu = 1.0_wp + 0.01*ht_edge
          vv = 0.0_wp 

          edge_vn = uu * patch_2d%edges%primal_normal(edge_index,edge_block)%v1 &
            & + vv * patch_2d%edges%primal_normal(edge_index,edge_block)%v2

          vn_test(edge_index, level, edge_block) = edge_vn

          ocean_state(jg)%p_prog(nold(1))%vn(edge_index, level, edge_block) &
              & = edge_vn

        ENDDO

      ENDDO
    ENDDO


    !! Get kinetic energy
    CALL calc_scalar_product_veloc_3d( patch_3d,  &
        & ocean_state(jg)%p_prog(nold(1))%vn,    &
        & ocean_state(jg)%p_diag,                     &
        & operators_coefficients)

    ! STEP 1: horizontal advection
    CALL veloc_adv_horz_mimetic( patch_3d,         &
      & ocean_state(jg)%p_prog(nold(1))%vn,    &
      & ocean_state(jg)%p_prog(nold(1))%vn,    &
      & ocean_state(jg)%p_diag,                &
      & ocean_state(jg)%p_diag%veloc_adv_horz, &
      & operators_coefficients)
    
    CALL dbg_print('test_mom: coriolis  ' ,ocean_state(jg)%p_diag%veloc_adv_horz ,&
        & '', 3, patch_2D%edges%owned )
    CALL dbg_print('test_mom: w',          ocean_state(jg)%p_diag%w,              &
        & '', 3, patch_2D%cells%owned )

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Transform vertical velocity using 
    !! wstar = ( w - u.grad( z ) )*dzstar/dz
    !! Only valid for static grid 
    IF (zcoord_type == 1) THEN
        DO blockNo = all_cells%start_block, all_cells%end_block
            CALL get_index_range(all_cells, blockNo, start_cells_index, end_cells_index)
            DO cell_index =  start_cells_index, end_cells_index
                H1   = patch_3d%p_patch_1d(1)%depth_CellInterface(cell_index, bt_level+1, blockNo)
                eta1 = eta(cell_index, blockNo) 

!                DO level = 1, patch_3d%p_patch_1d(1)%dolic_c(cell_index,blockNo)+1 
!                    ocean_state(jg)%p_diag%w(cell_index, level, blockNo) = &
!                      &  ( ocean_state(jg)%p_diag%w(cell_index, level, blockNo) &
!                      & - p_zgrad_top(cell_index, level, blockNo) )*( H1/(H1+eta1) )
!                ENDDO

            ENDDO
        ENDDO

    ENDIF

    CALL dbg_print('test_mom: w corr' , ocean_state(jg)%p_diag%w ,&
        & '', 3, patch_2D%cells%owned )


    ! STEP 2: compute 3D contributions: vertical velocity advection
    CALL veloc_adv_vert_mimetic(          &
      & patch_3d,                         &
      & ocean_state(jg)%p_diag,           &
      & operators_coefficients,           &
      & ocean_state(jg)%p_diag%veloc_adv_vert )

    !! Add kin. energy gradient, vorticity term and vertical mom. gradient
    mom_grad = ocean_state(jg)%p_diag%grad + ocean_state(jg)%p_diag%veloc_adv_horz&
        & + ocean_state(jg)%p_diag%veloc_adv_vert 

    CALL dbg_print('test_mom: vel. terms' ,mom_grad ,&
        & '', 3, patch_2D%edges%owned )

  END SUBROUTINE test_mom



  !-------------------------------------------------------------------------
  !>
  !! FIXME: Notes and overview
  !! Routines for testing advection with zstar are below
  !! Routines have been copied from other files and modified for zstar
  !! 1. Ideally, coefficients for the modified discretization should be calculated at one place
  !! and passed as arguments to repeat calculations
  !! 2. One routine can be used to calculate both low and high order flux for speedup
  !! 3. Variable to calculate depth needs to be clarified
  !! 4. This can be converted to using generalized vertical co-ordinates by using only 
  !! the coefficient of dz as variables to be modified
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! Flux limiter for horizontal advection
  !!
  !! Zalesak Flux-Limiter (Flux corrected transport)
  !! The corrected flux is a weighted average of the low order flux and the
  !! given high order flux. The high order flux is used to the greatest extent
  !! possible without introducing overshoots and undershoots.
  !! In vicinity of a lateral boundary only the low order flux is used: The criterion 
  !! is that at least one of the edges of the two neighboring cells of
  !! a central edges is a boundary edge.
  !! Note: This limiter is positive definite and almost monotone (but not strictly).
  !!
  !! @par Literature:
  !! - Zalesak, S.T. (1979): Fully Multidimensional Flux-corrected Transport
  !!   Algorithms for Fluids. JCP, 31, 335-362
  !!
  !! Adapted for zstar
  !! FIXME: The limiter assumes no knowledge of eta for next time step which 
  !! would be required if the formulation were to be correct
  !!
  SUBROUTINE limiter_ocean_zalesak_horz_zstar( patch_3d,&
    & vert_velocity,          &
    & tracer,                 &
    & p_mass_flx_e,           &
    & flx_tracer_low,         &    
    & flx_tracer_high,        &
    & div_adv_flux_vert,      &   
    & eta,                    &   
    & operators_coefficients, &
    & flx_tracer_final )       
    
    TYPE(t_patch_3d ),TARGET, INTENT(in):: patch_3d
    REAL(wp),INTENT(inout)              :: vert_velocity(nproma,n_zlev+1,patch_3d%p_patch_2d(1)%alloc_cell_blocks)    
    REAL(wp), INTENT(inout)             :: tracer           (nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp), INTENT(inout)             :: p_mass_flx_e     (nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp), INTENT(inout)             :: flx_tracer_low   (nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)     
    REAL(wp), INTENT(inout)             :: flx_tracer_high  (nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e) 
    REAL(wp), INTENT(inout)             :: flx_tracer_final (nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)     
    REAL(wp), INTENT(inout)             :: div_adv_flux_vert(nproma,n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks)    
    REAL(wp), INTENT(in)                 :: eta(nproma, patch_3d%p_patch_2d(1)%alloc_cell_blocks) !Surface ht 
    TYPE(t_operator_coeff),INTENT(in)   :: operators_coefficients
    
    !Local variables
    REAL(wp) :: z_mflx_anti(patch_3d%p_patch_2d(1)%cells%max_connectivity)
    REAL(wp) :: z_fluxdiv_c     !< flux divergence at cell center
    REAL(wp) :: z_anti          (nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)          !< antidiffusive tracer mass flux (F_H - F_L)    
    REAL(wp) :: z_tracer_new_low(nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)!< new tracer field after transport, if low order fluxes are used
    REAL(wp) :: z_tracer_max    (nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)!< local maximum of current tracer value and low order update
    REAL(wp) :: z_tracer_min    (nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)!< local minimum of current tracer value and low order update
    REAL(wp) :: r_p             (nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)!< fraction which must multiply all in/out fluxes of cell jc to guarantee
    REAL(wp) :: r_m             (nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)!< no overshoot/undershoot
    REAL(wp) :: z_tracer_update_horz(nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)!< new tracer field after transport, if low order fluxes are used    
    REAL(wp) :: r_frac          !< computed minimum fraction which must multiply< the flux at the edge
    REAL(wp) :: z_min, z_max    !< minimum/maximum value in cell and neighboring cells
    REAL(wp) :: z_signum        !< sign of antidiffusive velocity
    REAL(wp) :: p_p, p_m        !< sum of antidiffusive fluxes into and out of cell jc
    REAL(wp) :: prism_thick_old(n_zlev), inv_prism_thick_new(n_zlev)
    REAL(wp) :: delta_z_new, delta_z
    INTEGER, DIMENSION(:,:,:), POINTER ::  cellOfEdge_idx, cellOfEdge_blk
    INTEGER, DIMENSION(:,:,:), POINTER :: neighbor_cell_idx, neighbor_cell_blk
    INTEGER, DIMENSION(:,:,:), POINTER :: edge_of_cell_idx, edge_of_cell_blk
    INTEGER :: start_level, end_level            
    INTEGER :: start_index, end_index
    INTEGER :: edge_index, level, blockNo, jc,  cell_connect, sum_lsm_quad_edge, ctr
    TYPE(t_subset_range), POINTER :: edges_in_domain,  cells_in_domain
    TYPE(t_patch), POINTER :: patch_2d    
    
    INTEGER  :: bt_lev 
    REAL(wp) :: H_l, eta_l
    REAL(wp) :: coeff_l(nproma, patch_3d%p_patch_2d(1)%alloc_cell_blocks)      
    TYPE(t_subset_range), POINTER :: all_cells 

    !-------------------------------------------------------------------------
    patch_2d        => patch_3d%p_patch_2d(1)
    edges_in_domain => patch_2d%edges%in_domain
    cells_in_domain => patch_2d%cells%in_domain
    all_cells       => patch_2d%cells%ALL
    !-------------------------------------------------------------------------
    start_level = 1
    end_level   = n_zlev
    cellOfEdge_idx  => patch_2d%edges%cell_idx
    cellOfEdge_blk  => patch_2d%edges%cell_blk
    edge_of_cell_idx  => patch_2d%cells%edge_idx
    edge_of_cell_blk  => patch_2d%cells%edge_blk
    neighbor_cell_idx => patch_2d%cells%neighbor_idx
    neighbor_cell_blk => patch_2d%cells%neighbor_blk
    
#ifdef NAGFOR
    z_tracer_max(:,:,:) = 0.0_wp
    z_tracer_min(:,:,:) = 0.0_wp
    r_m(:,:,:)          = 0.0_wp
    r_p(:,:,:)          = 0.0_wp
#endif
 
    !-----------------------------------------------------------------------
    
    DO blockNo = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, blockNo, start_index, end_index)
      DO jc = start_index, end_index
        bt_lev = patch_3d%p_patch_1d(1)%dolic_c(jc,blockNo)      
 
        H_l   = patch_3d%p_patch_1d(1)%depth_CellInterface(jc, bt_lev + 1, blockNo)

        eta_l = eta(jc, blockNo) 

        coeff_l(jc, blockNo) = (H_l + eta_l)/H_l
        
      END DO
    END DO
    
    !-----------------------------------------------------------------------

  
!ICON_OMP_PARALLEL

!ICON_OMP_DO PRIVATE(start_index, end_index, edge_index, level) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, blockNo, start_index, end_index)
      
      z_anti(:,:,blockNo)     = 0.0_wp
      DO edge_index = start_index, end_index
        DO level = start_level, MIN(patch_3d%p_patch_1d(1)%dolic_e(edge_index,blockNo), end_level)
          
          ! calculate antidiffusive flux for each edge
          z_anti(edge_index,level,blockNo) = flx_tracer_high(edge_index,level,blockNo)&
                                          &- flx_tracer_low(edge_index,level,blockNo)
        END DO  ! end loop over edges
      END DO  ! end loop over levels
    END DO  ! end loop over blocks
!ICON_OMP_END_DO

    
!ICON_OMP_DO PRIVATE(start_index, end_index, jc, level, delta_z, delta_z_new, &
!ICON_OMP z_fluxdiv_c) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, blockNo, start_index, end_index)
      
      z_tracer_new_low(:,:,blockNo)    = 0.0_wp
      z_tracer_update_horz(:,:,blockNo)= 0.0_wp
      z_tracer_max(:,:,blockNo)        = 0.0_wp
      z_tracer_min(:,:,blockNo)        = 0.0_wp

      DO jc = start_index, end_index
        IF (patch_3d%p_patch_1d(1)%dolic_c(jc,blockNo) < 1) CYCLE
        
        ! 3. Compute the complete (with horizontal and vertical divergence) updated low order solution z_tracer_new_low
        DO level = start_level  , MIN(patch_3d%p_patch_1d(1)%dolic_c(jc,blockNo), end_level)       
          !  compute divergence of low order fluxes
          z_fluxdiv_c = 0
          DO cell_connect = 1, patch_2d%cells%num_edges(jc,blockNo)
            z_fluxdiv_c =  z_fluxdiv_c + &
              & flx_tracer_low(edge_of_cell_idx(jc,blockNo,cell_connect),level,edge_of_cell_blk(jc,blockNo,cell_connect)) * &
              & operators_coefficients%div_coeff(jc,level,blockNo,cell_connect)
          ENDDO

          delta_z     = coeff_l(jc, blockNo)*patch_3d%p_patch_1D(1)%prism_thick_flat_sfc_c(jc,level,blockNo)
          delta_z_new = coeff_l(jc, blockNo)*patch_3d%p_patch_1D(1)%prism_thick_flat_sfc_c(jc,level,blockNo)
          !
              
          z_tracer_new_low(jc,level,blockNo) = (tracer(jc,level,blockNo) * delta_z                     &
            & - dtime * (z_fluxdiv_c+div_adv_flux_vert(jc,level,blockNo)))/delta_z_new
            
        ENDDO
      ENDDO
      
      ! precalculate local maximum/minimum of current tracer value and low order
      ! updated value
      z_tracer_max(:,:,blockNo) =            &
        & MAX(          tracer(:,:,blockNo), &
        &     z_tracer_new_low(:,:,blockNo))
      z_tracer_min(:,:,blockNo) =            &
        & MIN(          tracer(:,:,blockNo), &
        &     z_tracer_new_low(:,:,blockNo))

    ENDDO
!ICON_OMP_END_DO

!ICON_OMP_MASTER
    CALL sync_patch_array_mult(sync_c1, patch_2d, 2, z_tracer_max, z_tracer_min)
!ICON_OMP_END_MASTER
!ICON_OMP_BARRIER
    ! 4. Limit the antidiffusive fluxes z_mflx_anti, such that the updated tracer
    !    field is free of any new extrema.    
!ICON_OMP_DO PRIVATE(start_index, end_index, jc, level, inv_prism_thick_new, &
!ICON_OMP z_mflx_anti, z_max, z_min, cell_connect, p_p, p_m) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block

      ! this is only needed for the parallel test setups
      ! it will try  tocheck the uninitialized (land) parts
      r_m(:,:,blockNo) = 0.0_wp
      r_p(:,:,blockNo) = 0.0_wp
        
      CALL get_index_range(cells_in_domain, blockNo, start_index, end_index)
      DO jc = start_index, end_index
        
        ! get prism thickness
        DO level = start_level  , MIN(patch_3d%p_patch_1d(1)%dolic_c(jc,blockNo), end_level)
          inv_prism_thick_new(level) = patch_3D%p_patch_1d(1)%inv_prism_thick_c(jc,level,blockNo) &
            & /coeff_l(jc, blockNo)
        ENDDO
        
        DO level = start_level, MIN(patch_3d%p_patch_1d(1)%dolic_c(jc,blockNo), end_level)         
          ! 2. Define "antidiffusive" fluxes A(jc,level,blockNo,edge_index) for each cell. It is the difference
          !    between the high order fluxes (given by the FFSL-scheme) and the low order
          !    ones. Multiply with geometry factor to have units [kg/kg] and the correct sign.
          !    - positive for outgoing fluxes
          !    - negative for incoming fluxes
          !    this sign convention is related to the definition of the divergence operator.          
          z_mflx_anti(:) = 0.0_wp
          z_max = z_tracer_max(jc,level,blockNo)
          z_min = z_tracer_min(jc,level,blockNo)
          p_p = 0.0_wp
          p_m = 0_wp
          DO cell_connect = 1, patch_2d%cells%num_edges(jc,blockNo)
            IF (patch_3d%p_patch_1d(1)% &
              & dolic_c(neighbor_cell_idx(jc,blockNo,cell_connect), neighbor_cell_blk(jc,blockNo,cell_connect)) >= level) THEN
              
              z_max = MAX(z_max, &
                & z_tracer_max(neighbor_cell_idx(jc,blockNo,cell_connect),level,neighbor_cell_blk(jc,blockNo,cell_connect)))
              z_min = MIN(z_min, &
                & z_tracer_min(neighbor_cell_idx(jc,blockNo,cell_connect),level,neighbor_cell_blk(jc,blockNo,cell_connect)))

              z_mflx_anti(cell_connect) =                                                        &
                & dtime * operators_coefficients%div_coeff(jc,level,blockNo,cell_connect) * inv_prism_thick_new(level)  &
                & * z_anti(edge_of_cell_idx(jc,blockNo,cell_connect),level,edge_of_cell_blk(jc,blockNo,cell_connect))

              ! Sum of all incoming antidiffusive fluxes into cell jc
              ! outgoing fluxes carry a positive sign, incoming a negative
              p_p = p_p - MIN(0._wp, z_mflx_anti(cell_connect))
              ! Sum of all outgoing antidiffusive fluxes out of cell jc
              p_m = p_m + MAX(0._wp, z_mflx_anti(cell_connect))
            ENDIF
          ENDDO                
          ! fraction which must multiply all fluxes out of cell jc to guarantee no
          ! undershoot
          ! Nominator: maximum allowable decrease of tracer
          r_m(jc,level,blockNo) = (z_tracer_new_low(jc,level,blockNo) - z_min ) / (p_m + dbl_eps)!&
          !
          ! fraction which must multiply all fluxes into cell jc to guarantee no
          ! overshoot
          ! Nominator: maximum allowable increase of tracer
          r_p(jc,level,blockNo) = (z_max - z_tracer_new_low(jc,level,blockNo)) / (p_p + dbl_eps)!&
          !
          !update old tracer with low-order flux
        ENDDO
      ENDDO
    ENDDO
!ICON_OMP_END_DO

    
!ICON_OMP_MASTER
    ! Synchronize r_m and r_p
    CALL sync_patch_array_mult(sync_c1, patch_2d, 2, r_m, r_p)
!ICON_OMP_END_MASTER
!ICON_OMP_BARRIER   

    ! 5. Now loop over all edges and determine the minimum fraction which must
    !    multiply the antidiffusive flux at the edge.
    !    At the end, compute new, limited fluxes which are then passed to the main
!ICON_OMP_DO PRIVATE(start_index, end_index, edge_index, level, z_signum, r_frac) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, blockNo, start_index, end_index)
      flx_tracer_final(:,:,blockNo) = 0.0_wp
      DO edge_index = start_index, end_index
      
        DO level = start_level, MIN(patch_3d%p_patch_1d(1)%dolic_e(edge_index,blockNo), end_level)
        
          IF( operators_coefficients%edges_SeaBoundaryLevel(edge_index,level,blockNo) > -2)THEN! edge < 2nd order boundary
          
            flx_tracer_final(edge_index,level,blockNo) = flx_tracer_low(edge_index,level,blockNo)
            
          ELSE!IF(sum_lsm_quad_edge==all_water_edges)THEN
          
            !z_anti>0 returns  1: here z_anti is outgoing, i.e. flux_high>flux_low
            !z_anti<0 returns -1: here z_anti is ingoing, i.e. flux_high<flux_low
            z_signum = SIGN(1._wp, z_anti(edge_index,level,blockNo))
                    
          ! This does the same as an IF (z_signum > 0) THEN ... ELSE ... ENDIF,
          ! but is computationally more efficient
          r_frac = 0.5_wp * (       &
            & (1._wp + z_signum) * & !<- active for z_signum=1
            & MIN(r_m(cellOfEdge_idx(edge_index,blockNo,1),level,cellOfEdge_blk(edge_index,blockNo,1)),  &
            &     r_p(cellOfEdge_idx(edge_index,blockNo,2),level,cellOfEdge_blk(edge_index,blockNo,2)))  &
            &+(1._wp - z_signum) * & !<- active for z_signum=-1
            & MIN(r_m(cellOfEdge_idx(edge_index,blockNo,2),level,cellOfEdge_blk(edge_index,blockNo,2)),  &
            &     r_p(cellOfEdge_idx(edge_index,blockNo,1),level,cellOfEdge_blk(edge_index,blockNo,1)))  )
          
          ! Limited flux
          flx_tracer_final(edge_index,level,blockNo) = flx_tracer_low(edge_index,level,blockNo)&
           & + MIN(1.0_wp,r_frac) *z_anti(edge_index,level,blockNo)      

            ENDIF
        END DO
       ENDDO
    ENDDO
!ICON_OMP_END_DO NOWAIT
!ICON_OMP_END_PARALLEL


  END SUBROUTINE limiter_ocean_zalesak_horz_zstar
  !-------------------------------------------------------------------------

  

  !-----------------------------------------------------------------------------
  ! the map_edges2edges_viacell_3d_mlev_constZs optimized for triangles
  ! modified for zstar 
  SUBROUTINE map_edges2edges_zstar( patch_3d, vn_e, eta, scalar, operators_coefficients, out_vn_e)
    
    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    REAL(wp), INTENT(in)                 :: vn_e(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp), INTENT(in)                 :: eta(nproma, patch_3d%p_patch_2d(1)%alloc_cell_blocks) !Surface ht 
    TYPE(t_operator_coeff), INTENT(in)   :: operators_coefficients
    REAL(wp), INTENT(inout)              :: out_vn_e(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp), INTENT(in)                 :: scalar(nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    !Local variables
    INTEGER :: startLevel, endLevel
    TYPE(t_subset_range), POINTER :: edges_in_domain
    TYPE(t_subset_range), POINTER :: all_edges 
    TYPE(t_subset_range), POINTER :: all_cells 
    TYPE(t_patch), POINTER :: patch_2d
    REAL(wp), POINTER :: coeffs(:,:,:,:)

    ! omp private
    ! defined in a struct for omp clearness
    TYPE omp_local_private
    
      INTEGER :: start_edge_index, end_edge_index
      INTEGER :: cell_1_index, cell_2_index, cell_1_block, cell_2_block
      INTEGER :: edge_11_index, edge_12_index, edge_13_index ! edges of cell_1
      INTEGER :: edge_11_block, edge_12_block, edge_13_block
      INTEGER :: edge_21_index, edge_22_index, edge_23_index ! edges of cell_2
      INTEGER :: edge_21_block, edge_22_block, edge_23_block
      
    END TYPE omp_local_private

    REAL(wp) :: coeff_e(nproma, patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp) :: H1, H2, H_e, eta_e 
    INTEGER  :: bt_lev 
    REAL(wp) :: H_l, eta_l
    REAL(wp) :: coeff_l(nproma, patch_3d%p_patch_2d(1)%alloc_cell_blocks)      
    INTEGER  :: start_index, end_index

    INTEGER, DIMENSION(:,:,:), POINTER :: idx, blk
    INTEGER  :: id1, id2, bl1, bl2 

!     TYPE(omp_local_private) :: omp_this 
    INTEGER :: jc, je, blockNo, level
    INTEGER :: start_edge_index, end_edge_index
    INTEGER :: cell_1_index, cell_2_index, cell_1_block, cell_2_block
    INTEGER :: edge_11_index, edge_12_index, edge_13_index ! edges of cell_1
    INTEGER :: edge_11_block, edge_12_block, edge_13_block
    INTEGER :: edge_21_index, edge_22_index, edge_23_index ! edges of cell_2
    INTEGER :: edge_21_block, edge_22_block, edge_23_block
   !-----------------------------------------------------------------------
    IF (no_primal_edges /= 3) &
      & CALL finish ('map_edges2edges_viacell triangle version', 'no_primal_edges /= 3')
    
    !-----------------------------------------------------------------------
    patch_2d   => patch_3d%p_patch_2d(1)
    all_edges  => patch_2d%edges%ALL
    all_cells  => patch_2d%cells%ALL
    edges_in_domain => patch_2d%edges%in_domain
    startLevel = 1
    endLevel = n_zlev
    coeffs => operators_coefficients%edge2edge_viacell_coeff
    idx      => patch_3D%p_patch_2D(1)%edges%cell_idx
    blk      => patch_3D%p_patch_2D(1)%edges%cell_blk
 
    !-----------------------------------------------------------------------
    DO blockNo = all_edges%start_block, all_edges%end_block
      CALL get_index_range(all_edges, blockNo, start_edge_index, end_edge_index)
      DO je =  start_edge_index, end_edge_index
        !! FIXME: Assuming levels are the same at the bottom
        bt_lev = patch_3d%p_patch_1d(1)%dolic_e(je,blockNo)      
 
        id1 = idx(je, blockNo, 1)
        id2 = idx(je, blockNo, 2)
        bl1 = blk(je, blockNo, 1)
        bl2 = blk(je, blockNo, 2)
 
        H1  = patch_3d%p_patch_1d(1)%depth_CellInterface(id1, bt_lev + 1, bl1)
        H2  = patch_3d%p_patch_1d(1)%depth_CellInterface(id2, bt_lev + 1, bl2)

        H_e = 0.5_wp*(H1 + H2)

        eta_e = 0.5_wp*( eta(id1, bl1) + eta(id2, bl2) )

        coeff_e(je, blockNo) = (H_e + eta_e)/H_e
        
      END DO
    END DO

    !-----------------------------------------------------------------------
    
    DO blockNo = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, blockNo, start_index, end_index)
      DO jc = start_index, end_index
        bt_lev = patch_3d%p_patch_1d(1)%dolic_c(jc,blockNo)      
 
        H_l   = patch_3d%p_patch_1d(1)%depth_CellInterface(jc, bt_lev + 1, blockNo)

        eta_l = eta(jc, blockNo) 

        coeff_l(jc, blockNo) = (H_l + eta_l)/H_l
        
      END DO
    END DO
    
    !-----------------------------------------------------------------------


    !-----------------------------------------------------------------------
    
!ICON_OMP_PARALLEL_DO PRIVATE( je, level, start_edge_index, end_edge_index, &
!ICON_OMP cell_1_index, cell_2_index, cell_1_block, cell_2_block, &
!ICON_OMP edge_11_index, edge_12_index, edge_13_index, &
!ICON_OMP edge_11_block, edge_12_block, edge_13_block, &
!ICON_OMP edge_21_index, edge_22_index, edge_23_index, &
!ICON_OMP edge_21_block, edge_22_block, edge_23_block) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, blockNo, start_edge_index, end_edge_index)
      
      DO je =  start_edge_index, end_edge_index
        
        cell_1_index = patch_2d%edges%cell_idx(je,blockNo,1)
        cell_1_block = patch_2d%edges%cell_blk(je,blockNo,1)
        cell_2_index = patch_2d%edges%cell_idx(je,blockNo,2)
        cell_2_block = patch_2d%edges%cell_blk(je,blockNo,2)
        
        edge_11_index = patch_2d%cells%edge_idx(cell_1_index, cell_1_block, 1)
        edge_12_index = patch_2d%cells%edge_idx(cell_1_index, cell_1_block, 2)
        edge_13_index = patch_2d%cells%edge_idx(cell_1_index, cell_1_block, 3)
        edge_11_block = patch_2d%cells%edge_blk(cell_1_index, cell_1_block, 1)
        edge_12_block = patch_2d%cells%edge_blk(cell_1_index, cell_1_block, 2)
        edge_13_block = patch_2d%cells%edge_blk(cell_1_index, cell_1_block, 3)
        
        edge_21_index = patch_2d%cells%edge_idx(cell_2_index, cell_2_block, 1)
        edge_22_index = patch_2d%cells%edge_idx(cell_2_index, cell_2_block, 2)
        edge_23_index = patch_2d%cells%edge_idx(cell_2_index, cell_2_block, 3)
        edge_21_block = patch_2d%cells%edge_blk(cell_2_index, cell_2_block, 1)
        edge_22_block = patch_2d%cells%edge_blk(cell_2_index, cell_2_block, 2)
        edge_23_block = patch_2d%cells%edge_blk(cell_2_index, cell_2_block, 3)
                
        ! levels
        DO level = startLevel, patch_3d%p_patch_1d(1)%dolic_e(je,blockNo)
          
          out_vn_e(je, level, blockNo) =  &
            & (  vn_e(edge_11_index, level, edge_11_block) * coeffs(je, level, blockNo, 1)      &
            &    * patch_3d%p_patch_1d(1)%prism_thick_e(edge_11_index, level, edge_11_block)    &
            & * coeff_e(edge_11_index, edge_11_block)                                           &
            &  + vn_e(edge_12_index, level, edge_12_block) * coeffs(je, level, blockNo, 2)      &
            &    * patch_3d%p_patch_1d(1)%prism_thick_e(edge_12_index, level, edge_12_block)    &
            & * coeff_e(edge_12_index, edge_12_block)                                           &
            &  + vn_e(edge_13_index, level, edge_13_block) * coeffs(je, level, blockNo, 3)      &
            &    * patch_3d%p_patch_1d(1)%prism_thick_e(edge_13_index, level, edge_13_block)    &
            & * coeff_e(edge_13_index, edge_13_block)                                           &
            & ) * scalar(cell_1_index, level, cell_1_block)                                     &
            & + &
            & (  vn_e(edge_21_index, level, edge_21_block) * coeffs(je, level, blockNo, 4)           &
            &   * patch_3d%p_patch_1d(1)%prism_thick_e(edge_21_index, level, edge_21_block)     &
            & * coeff_e(edge_21_index, edge_21_block)                                           &
            &  + vn_e(edge_22_index, level, edge_22_block) * coeffs(je, level, blockNo, 5)      &
            &  * patch_3d%p_patch_1d(1)%prism_thick_e(edge_22_index, level, edge_22_block)      &
            & * coeff_e(edge_22_index, edge_22_block)                                           &
            &  + vn_e(edge_23_index, level, edge_23_block) * coeffs(je, level, blockNo, 6)      &
            & * patch_3d%p_patch_1d(1)%prism_thick_e(edge_23_index, level, edge_23_block)       &
            & * coeff_e(edge_23_index, edge_23_block)                                           &
            & ) * scalar(cell_2_index, level, cell_2_block)                                      
      
        END DO !levels
    
      END DO
      
    END DO ! blockNo = edges_in_domain%start_block, edges_in_domain%end_block
!ICON_OMP_END_PARALLEL_DO
     
  END SUBROUTINE map_edges2edges_zstar
  !-----------------------------------------------------------------------------
  




  
  SUBROUTINE upwind_zstar_hflux_oce( patch_3d, cell_value, eta, edge_vn, edge_upwind_flux, opt_start_level, opt_end_level )
    
    TYPE(t_patch_3d ),TARGET, INTENT(in)   :: patch_3d
    REAL(wp), INTENT(in)              :: cell_value   (nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)      !< advected cell centered variable
    REAL(wp), INTENT(in)              :: eta(nproma, patch_3d%p_patch_2d(1)%alloc_cell_blocks) !Surface ht 
    REAL(wp), INTENT(in)              :: edge_vn    (nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)       !< normal velocity on edges
    REAL(wp), INTENT(inout)           :: edge_upwind_flux(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)   !< variable in which the upwind flux is stored
    INTEGER, INTENT(in), OPTIONAL :: opt_start_level    ! optional vertical start level
    INTEGER, INTENT(in), OPTIONAL :: opt_end_level    ! optional vertical end level
    ! local variables
    INTEGER, DIMENSION(:,:,:), POINTER :: iilc,iibc  ! pointer to line and block indices
    INTEGER, DIMENSION(:,:,:), POINTER :: idx, blk
    INTEGER  :: start_level, end_level
    INTEGER  :: start_index, end_index
    INTEGER  :: edge_index, level, blockNo         !< index of edge, vert level, block
    INTEGER  :: id1, id2, bl1, bl2 
    INTEGER  :: bt_level 
    REAL(wp) :: H1, H2, coeff1, coeff2, eta1, eta2 
    TYPE(t_subset_range), POINTER :: edges_in_domain
    TYPE(t_patch), POINTER :: patch_2d
    !-----------------------------------------------------------------------
    patch_2d        => patch_3d%p_patch_2d(1)
    edges_in_domain => patch_2d%edges%in_domain
    idx             => patch_3D%p_patch_2D(1)%edges%cell_idx
    blk             => patch_3D%p_patch_2D(1)%edges%cell_blk
    !-----------------------------------------------------------------------
    IF ( PRESENT(opt_start_level) ) THEN
      start_level = opt_start_level
    ELSE
      start_level = 1
    END IF
    IF ( PRESENT(opt_end_level) ) THEN
      end_level = opt_end_level
    ELSE
      end_level = n_zlev
    END IF
    !
    ! advection is done with 1st order upwind scheme,
    ! i.e. a piecewise constant approx. of the cell centered values
    ! is used.
    !
!ICON_OMP_PARALLEL PRIVATE(iilc, iibc)
    ! line and block indices of two neighboring cells
    iilc => patch_2d%edges%cell_idx
    iibc => patch_2d%edges%cell_blk
    
    ! loop through all patch edges (and blocks)
!ICON_OMP_DO PRIVATE(start_index, end_index, edge_index, level) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, blockNo, start_index, end_index)
      edge_upwind_flux(:,:,blockNo) = 0.0_wp
      DO edge_index = start_index, end_index
        bt_level = patch_3d%p_patch_1d(1)%dolic_e(edge_index,blockNo)      

        id1 = idx(edge_index, blockNo, 1)
        id2 = idx(edge_index, blockNo, 2)
        bl1 = blk(edge_index, blockNo, 1)
        bl2 = blk(edge_index, blockNo, 2)
 
        !! Get height of the cell center at mid point from the bottom
        H1  = patch_3d%p_patch_1d(1)%depth_CellInterface(id1,bt_level+1,bl1)
        H2  = patch_3d%p_patch_1d(1)%depth_CellInterface(id2,bt_level+1,bl2) 
        
        eta1 = eta(id1, bl1) 
        eta2 = eta(id2, bl2) 

        !! Transform to z from zstar
        coeff1  = (H1 + eta1)/H1 
        coeff2  = (H2 + eta2)/H2 

        DO level = start_level, MIN(patch_3d%p_patch_1d(1)%dolic_e(edge_index,blockNo), end_level)
          !
          ! compute the first order upwind flux; notice
          ! that multiplication by edge length is avoided to
          ! compute final conservative update using the discrete
          ! div operator
          edge_upwind_flux(edge_index,level,blockNo) =  &
             0.5_wp * (        edge_vn(edge_index,level,blockNo)  *           &
               & ( coeff1*cell_value(iilc(edge_index,blockNo,1),level,iibc(edge_index,blockNo,1)) + &
               &   coeff2*cell_value(iilc(edge_index,blockNo,2),level,iibc(edge_index,blockNo,2)) ) &
               &   - ABS( edge_vn(edge_index,level,blockNo) ) *               &
               & ( coeff2*cell_value(iilc(edge_index,blockNo,2),level,iibc(edge_index,blockNo,2)) - &
               &   coeff1*cell_value(iilc(edge_index,blockNo,1),level,iibc(edge_index,blockNo,1)) ) )
          
        END DO  ! end loop over edges
      END DO  ! end loop over levels
    END DO  ! end loop over blocks
!ICON_OMP_END_DO NOWAIT
!ICON_OMP_END_PARALLEL
    
  END SUBROUTINE upwind_zstar_hflux_oce
  !-----------------------------------------------------------------------
   


  !-------------------------------------------------------------------------
  !>
  !! Computation of new vertical velocity using continuity equation
  !! Calculate diagnostic vertical velocity from horizontal velocity using the
  !! incommpressibility condition in the continuity equation.
  !! vertical velocity is integrated from bottom to topLevel
  !! vertical velocity is negative for positive divergence
  !! of horizontal velocity
  !!
  SUBROUTINE calc_vert_velocity_zstar( patch_3d, ocean_state, op_coeffs, eta)
    TYPE(t_patch_3d), TARGET :: patch_3d       ! patch on which computation is performed
    TYPE(t_hydro_ocean_state) :: ocean_state
    TYPE(t_operator_coeff), INTENT(in) :: op_coeffs
    REAL(wp), INTENT(in)               :: eta(nproma, patch_3d%p_patch_2d(1)%alloc_cell_blocks) !Surface ht 
    ! Local variables
    INTEGER :: jc, jk, blockNo, je, z_dolic, start_index, end_index
    REAL(wp) :: z_c(nproma,patch_3d%p_patch_2d(1)%alloc_cell_blocks), z_abort
    TYPE(t_subset_range), POINTER :: cells_in_domain, edges_in_domain, all_cells, cells_owned
    TYPE(t_patch), POINTER :: patch_2D
    REAL(wp), POINTER :: vertical_velocity(:,:,:)
    REAL(wp) :: H_l, eta_l, coeff_l 
    INTEGER  :: bt_level 
    REAL(wp) :: z_adv_flux_h (nproma, n_zlev, patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp) :: div_z_c(nproma,n_zlev)
    REAL(wp) :: div_z_depth_int_c(nproma)
    REAL(wp) :: temp(nproma,n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp), POINTER :: cell_thickness(:,:,:)

    CHARACTER(len=*), PARAMETER :: method_name='mo_ocean_ab_timestepping_mimetic:alc_vert_velocity_mim_bottomup'
    !-----------------------------------------------------------------------
    patch_2D         => patch_3d%p_patch_2d(1)
    cells_in_domain  => patch_2D%cells%in_domain
    cells_owned      => patch_2D%cells%owned
    all_cells        => patch_2D%cells%all
    edges_in_domain  => patch_2D%edges%in_domain
    vertical_velocity=> ocean_state%p_diag%w
    cell_thickness   => patch_3D%p_patch_1d(1)%prism_thick_c
    ! due to nag -nan compiler-option:
    !------------------------------------------------------------------
    ! Step 1) Calculate divergence of horizontal velocity at all levels
    !------------------------------------------------------------------
    !-------------------------------------------------------------------------------
 
    CALL map_edges2edges_viacell_3d_const_z( patch_3d, ocean_state%p_prog(nold(1))%vn, &
      & op_coeffs, ocean_state%p_diag%mass_flx_e)

    !! Use trick of using constant coefficient to get mass flux 
    temp = 1.0_wp
  
    CALL map_edges2edges_zstar( patch_3d, ocean_state%p_prog(nold(1))%vn, eta, &
      & temp, op_coeffs, z_adv_flux_h)

    vertical_velocity = 0.0_wp

!ICON_OMP_PARALLEL_DO PRIVATE(start_index,end_index, jc, jk) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, blockNo, start_index, end_index)
      
      CALL div_oce_3D_onTriangles_onBlock(z_adv_flux_h, patch_3D, op_coeffs%div_coeff, &
        & div_z_c(:,:), blockNo=blockNo, start_index=start_index, &
        & end_index=end_index, start_level=1, end_level=n_zlev)

      DO jc = start_index, end_index
        !! Get summation over depth of divergence for RHS of sfc equation 
        div_z_depth_int_c(jc) = SUM(div_z_c(jc, 1:patch_3d%p_patch_1d(1)%dolic_c(jc,blockNo)))
       
        bt_level = patch_3d%p_patch_1d(1)%dolic_c(jc, blockNo)      

        !! Get height of the cell center at mid point from the bottom
        H_l    = patch_3d%p_patch_1d(1)%depth_CellInterface(jc, bt_level+1, blockNo)
        
        eta_l  = eta(jc, blockNo) 

        !! Transform to z from zstar
        coeff_l = (H_l + eta_l)/H_l 

        !! Trick here is that dw/dz = -div(v) or w2-w1=-dz.div(v)=-div(dz.v)
        !! In z*, it becomes d(J.w)/dz = -div(J.v) or w2-w1=-dz*.div(J.v)/J=-div(dz*.J.v)/J
        !! J being the vertical co-ordinate Jacobian
        !! d_t eta = -div_z_depth_int_c(jc)/H_l
        DO jk = patch_3d%p_patch_1d(1)%dolic_c(jc,blockNo), 1, -1
        
          vertical_velocity(jc,jk,blockNo) = vertical_velocity(jc,jk+1,blockNo) - &
            & ( div_z_c(jc,jk) - &
            & cell_thickness(jc, jk, blockNo)*div_z_depth_int_c(jc)/H_l)/coeff_l
          !! Switch to the below lines to use w^* for d_t eta = 0
          !! This will give incorrect results for topmost layer where a zero
          !! flux boundary condn is used in vertical advection
!          vertical_velocity(jc,jk,blockNo) = vertical_velocity(jc,jk+1,blockNo) - &
!            &  div_z_c(jc,jk)/coeff_l


        END DO
      END DO
    END DO ! blockNo
!ICON_OMP_END_PARALLEL_DO
    
    CALL sync_patch_array(sync_c,patch_2D,vertical_velocity)

  END SUBROUTINE calc_vert_velocity_zstar
  !-------------------------------------------------------------------------
  

  
  !-------------------------------------------------------------------------
  !> Setup a test case that advects tracers for testing with zstar
  !  Should start by using the low order horizontal advection
  SUBROUTINE ocean_test_zstar_advection( patch_3d, ocean_state, &
    & this_datetime, ocean_surface, physics_parameters,             &
    & ocean_ice,operators_coefficients)
    
    TYPE(t_patch_3d), POINTER, INTENT(in)          :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET, INTENT(inout) :: ocean_state(n_dom)
    TYPE(datetime), POINTER                          :: this_datetime
    TYPE(t_ocean_surface)                            :: ocean_surface
    TYPE (t_ho_params)                               :: physics_parameters
    TYPE (t_sea_ice),         INTENT(inout)          :: ocean_ice
    TYPE(t_operator_coeff),   INTENT(in)          :: operators_coefficients
    
    ! local variables
    TYPE (t_hamocc_state)        :: hamocc_State
    INTEGER :: jstep, jg
    !LOGICAL                         :: l_outputtime
    CHARACTER(LEN=32)               :: datestring
    TYPE(t_patch), POINTER :: patch_2d
    INTEGER :: jstep0 ! start counter for time loop
    INTEGER :: i
    INTEGER :: tracer_index 
    TYPE(timedelta), POINTER :: model_time_step => NULL()
    
    !CHARACTER(LEN=filename_max)  :: outputfile, gridfile
    TYPE(t_tracer_collection) , POINTER              :: old_tracer_collection, new_tracer_collection
    TYPE(t_ocean_transport_state)                    :: transport_state
    
    INTEGER  :: jb, jc, je, level 
    REAL(wp) :: eta(nproma, patch_3d%p_patch_2d(1)%alloc_cell_blocks) 
    REAL(wp) :: delta_t, delta_z,delta_z_new, delta_z1,delta_z_new1
    REAL(wp) :: div_adv_flux_horz(nproma,n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp) :: div_diff_flux_horz(nproma,n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp) :: flux_horz(nproma,n_zlev, patch_3d%p_patch_2D(1)%nblks_e)
    REAL(wp) :: div_adv_flux_vert(nproma,n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp) :: top_bc(nproma)
    INTEGER  :: start_index, end_index
    INTEGER  :: start_cell_index, end_cell_index
    REAL(wp) :: z_adv_flux_h (nproma, n_zlev, patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp) :: z_adv_low (nproma, n_zlev, patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp) :: z_adv_high(nproma, n_zlev, patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp) :: z2(nproma, n_zlev, patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp) :: H_l, eta_l, coeff_l 
    INTEGER  :: bt_level 
    
    TYPE(t_ocean_tracer), POINTER :: new_tracer
    TYPE(t_ocean_tracer), POINTER :: old_tracer

    TYPE(t_subset_range), POINTER :: cells_in_domain, edges_in_domain
    
    REAL(wp) :: temp(nproma,n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks)

    CHARACTER(LEN=max_char_length), PARAMETER :: &
      & method_name = 'mo_ocean_testbed_modules:ocean_test_zstar_advection'
    !------------------------------------------------------------------
    
    patch_2D        => patch_3d%p_patch_2d(1)
    cells_in_domain => patch_2D%cells%in_domain
    edges_in_domain => patch_2d%edges%in_domain
    CALL datetimeToString(this_datetime, datestring)

    ! IF (ltimer) CALL timer_start(timer_total)
    CALL timer_start(timer_total)
    
    IF (n_dom > 1 ) THEN
      CALL finish(TRIM(method_name), ' N_DOM > 1 is not allowed')
    END IF
    jg = n_dom

    !! sea surface height type 201 set explicitly since we don't want
    !! the grid to change
    eta = 0.  
    ! Initialize eta for zstar
    ! #slo#: simple elevation between 30W and 30E (pi/3.)
    DO jb = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, jb, start_cell_index, end_cell_index)
      DO jc = start_cell_index, end_cell_index
        IF ( patch_3d%lsm_c(jc, 1, jb) <= sea_boundary ) THEN
           eta(jc, jb) = 10.0_wp * &
            & SIN(patch_2d%cells%center(jc, jb)%lon * 6.0_wp) &
            & * COS(patch_2d%cells%center(jc, jb)%lat * 3.0_wp)
        ENDIF
      END DO
    END DO

!    !---------------------------------------------------------------------
!    !-FIXME: test divergence of constant fn
!    !---------------------------------------------------------------------
!  
!    ! calc_vert_vel uses vn_time_weighter instead of vn
!    ocean_state(jg)%p_diag%vn_time_weighted = ocean_state(jg)%p_prog(nold(1))%vn
!    
!    !! Update mass_flux and w 
!    CALL calc_vert_velocity_zstar( patch_3d, ocean_state(jg),operators_coefficients, eta)
!
!    ! fill transport_state
!    transport_state%patch_3d    => patch_3d
!    transport_state%h_old       => ocean_state(jg)%p_prog(nold(1))%h
!    transport_state%h_new       => ocean_state(jg)%p_prog(nnew(1))%h
!    transport_state%vn          => ocean_state(jg)%p_prog(nold(1))%vn
!    transport_state%mass_flux_e => ocean_state(jg)%p_diag%mass_flx_e
!    transport_state%w           => ocean_state(jg)%p_diag%w
!
!    temp = 1.0_wp
!
!    CALL advect_flux_vertical( patch_3d,&
!        & temp, &
!        & transport_state,                           &
!        & operators_coefficients,                     &
!        & div_adv_flux_vert)
!
!    DO jb = cells_in_domain%start_block, cells_in_domain%end_block
!      CALL get_index_range(cells_in_domain, jb, start_cell_index, end_cell_index)
!      DO jc = start_cell_index, end_cell_index
!        bt_level = patch_3d%p_patch_1d(1)%dolic_c(jc, jb)      
!
!        !! Get height of the cell center at mid point from the bottom
!        H_l    = patch_3d%p_patch_1d(1)%depth_CellInterface(jc, bt_level+1, jb)
!        
!        eta_l  = eta(jc, jb) 
!
!        !! Transform to z from zstar
!        coeff_l = (H_l + eta_l)/H_l 
!
!        div_adv_flux_vert(jc, :, jb) = coeff_l*div_adv_flux_vert(jc, :, jb)
!      END DO
!    END DO
!
!    !---------------------------------------------------------------------
!    !-Horizontal  advection
!    !---------------------------------------------------------------------
!    CALL upwind_zstar_hflux_oce( patch_3d,  &
!      & temp, &
!      & eta, &
!      & transport_state%mass_flux_e,         &
!      & z_adv_flux_h)                         
! 
!    !Calculate divergence of advective fluxes
!    CALL div_oce_3d( z_adv_flux_h, patch_3D, operators_coefficients%div_coeff, &
!      & div_adv_flux_horz, subset_range=cells_in_domain )
!
!    CALL map_edges2edges_zstar( patch_3d, transport_state%vn, eta, &
!      & temp, operators_coefficients, z2)
!
!    DO jb = cells_in_domain%start_block, cells_in_domain%end_block
!      CALL get_index_range(cells_in_domain, jb, start_cell_index, end_cell_index)
!      DO jc = start_cell_index, end_cell_index
!        DO level = 1, patch_3d%p_patch_1d(1)%dolic_c(jc,jb)
!
!          if ( (jb == 4) .AND. (jc == 10) ) THEN
!!            write(*, *) level, div_adv_flux_horz(jc,level,jb) , temp(jc,level,jb)
!!              & transport_state%w(jb, level + 1, jc) - transport_state%w(jb, level, jc), &
!!              & patch_3d%p_patch_1d(1)%prism_thick_c(jb, level, jc)
!          ENDIF
!    
!        ENDDO
!      ENDDO
!    ENDDO
! 
!    DO jb = edges_in_domain%start_block, edges_in_domain%end_block
!      CALL get_index_range(edges_in_domain, jb, start_index, end_index)
!      DO je = start_index, end_index
!        DO level = 1, patch_3d%p_patch_1d(1)%dolic_e(je,jb)
!           if ( (jb == 4) .AND. (je == 10) ) THEN
!            write(*, *) level, z_adv_flux_h(je,level,jb) , z2(je,level,jb)
!          ENDIF
!    
!        ENDDO
!      ENDDO
!    ENDDO
!
!
!    !---------------------------------------------------------------------
!    !-FIXME: test end 
!    !---------------------------------------------------------------------
 
    jstep0 = 0
    DO jstep = (jstep0+1), (jstep0+nsteps)

        CALL datetimeToString(this_datetime, datestring)
        WRITE(message_text,'(a,i10,2a)') '  Begin of timestep =',jstep,'  datetime:  ', datestring
        CALL message (TRIM(method_name), message_text)
 
        old_tracer_collection => ocean_state(jg)%p_prog(nold(1))%tracer_collection
        new_tracer_collection => ocean_state(jg)%p_prog(nnew(1))%tracer_collection

        IF (no_tracer>=1) THEN

          ! calc_vert_vel uses vn_time_weighter instead of vn
          ocean_state(jg)%p_diag%vn_time_weighted = ocean_state(jg)%p_prog(nold(1))%vn

          !! Update w and mass_flx_e for tracer advection
          CALL calc_vert_velocity_zstar( patch_3d, ocean_state(jg),operators_coefficients, eta)

          ! fill transport_state
          transport_state%patch_3d    => patch_3d
          transport_state%h_old       => ocean_state(jg)%p_prog(nold(1))%h
          transport_state%h_new       => ocean_state(jg)%p_prog(nnew(1))%h
          transport_state%vn          => ocean_state(jg)%p_prog(nold(1))%vn
          transport_state%w           => ocean_state(jg)%p_diag%w
          transport_state%mass_flux_e => ocean_state(jg)%p_diag%mass_flx_e

          ! fill diffusion coefficients
          old_tracer_collection%tracer(1)%hor_diffusion_coeff => physics_parameters%TracerDiffusion_coeff(:,:,:,1)
          old_tracer_collection%tracer(1)%ver_diffusion_coeff => physics_parameters%a_tracer_v(:,:,:,1)
          DO i = 2, old_tracer_collection%no_of_tracers
            old_tracer_collection%tracer(i)%hor_diffusion_coeff => physics_parameters%TracerDiffusion_coeff(:,:,:,2)
            old_tracer_collection%tracer(i)%ver_diffusion_coeff => physics_parameters%a_tracer_v(:,:,:,2)
          ENDDO
        ENDIF
        !------------------------------------------------------------------------

        DO tracer_index = 1, old_tracer_collection%no_of_tracers
          
          old_tracer => old_tracer_collection%tracer(tracer_index)
          new_tracer => new_tracer_collection%tracer(tracer_index)
          IF ( old_tracer%is_advected) THEN
            
            !-------------------------------------------------------------------------------
            patch_2D        => patch_3d%p_patch_2d(1)
            cells_in_domain => patch_2D%cells%in_domain
            edges_in_domain => patch_2D%edges%in_domain
            delta_t = dtime
            !---------------------------------------------------------------------
         
            ! these are probably not necessary
            div_adv_flux_vert = 0.0_wp
            div_adv_flux_horz = 0.0_wp
            div_diff_flux_horz = 0.0_wp
            
            !---------------------------------------------------------------------
            !-Vertical advection
            !---------------------------------------------------------------------
            IF ( l_with_vert_tracer_advection ) THEN
        
              CALL advect_flux_vertical( patch_3d,&
                & old_tracer%concentration, &
                & transport_state,                           &
                & operators_coefficients,                     &
                & div_adv_flux_vert)
        
            ENDIF  ! l_with_vert_tracer_advection
 
            !---------------------------------------------------------------------
            !-Horizontal  advection
            !---------------------------------------------------------------------
            CALL upwind_zstar_hflux_oce( patch_3d,  &
              & old_tracer%concentration, &
              & eta, &
              & transport_state%mass_flux_e,         &
              & z_adv_flux_h)
            z_adv_low = z_adv_flux_h
 
            CALL map_edges2edges_zstar( patch_3d, transport_state%vn, eta, &
              & old_tracer%concentration, operators_coefficients, z_adv_flux_h)
            z_adv_high = z_adv_flux_h

            CALL limiter_ocean_zalesak_horz_zstar( patch_3d,   &
              & transport_state%w,           &
              & old_tracer%concentration,              &
              & transport_state%mass_flux_e,           &
              & z_adv_low,                             &
              & z_adv_high,                            &
              & div_adv_flux_vert,                     &            
              & eta,                                   &            
              & operators_coefficients,                &
              & z_adv_flux_h)                          

            !Calculate divergence of advective fluxes
            CALL div_oce_3d( z_adv_flux_h, patch_3D, operators_coefficients%div_coeff, &
              & div_adv_flux_horz, subset_range=cells_in_domain )

        !ICON_OMP_PARALLEL_DO PRIVATE(start_cell_index, end_cell_index, jc, level, &
        !ICON_OMP delta_z, delta_z_new, top_bc) ICON_OMP_DEFAULT_SCHEDULE
            DO jb = cells_in_domain%start_block, cells_in_domain%end_block
              CALL get_index_range(cells_in_domain, jb, start_cell_index, end_cell_index)
              DO jc = start_cell_index, end_cell_index
                bt_level = patch_3d%p_patch_1d(1)%dolic_c(jc, jb)      

                !! Get height of the cell center at mid point from the bottom
                H_l    = patch_3d%p_patch_1d(1)%depth_CellInterface(jc, bt_level+1, jb)
        
                eta_l  = eta(jc, jb) 

                !! Transform to z from zstar
                coeff_l = (H_l + eta_l)/H_l 

                !! d_z*(coeff*w*C) = coeff*d_z(w*C) since coeff is constant for each column
                div_adv_flux_vert(jc, :, jb) = coeff_l*div_adv_flux_vert(jc, :, jb)

                DO level = 1, patch_3d%p_patch_1d(1)%dolic_c(jc,jb)
        
                  new_tracer%concentration(jc,level,jb) =                          &
                    &  old_tracer%concentration(jc,level,jb) -                     &
                    &  (delta_t /  ( coeff_l*patch_3d%p_patch_1D(1)%prism_thick_c(jc,level,jb) ) )    &
                    & * (div_adv_flux_horz(jc,level,jb)  + div_adv_flux_vert(jc,level,jb))

!                  write(*, *) level, ocean_state(jg)%p_diag%w (jc,level,jb), &
!                    & div_adv_flux_horz(jc,level,jb) , div_adv_flux_vert(jc,level,jb)

                ENDDO
        
              END DO
            END DO
        !ICON_OMP_END_PARALLEL_DO
        
        
            CALL sync_patch_array(sync_c, patch_2D, new_tracer%concentration)
        
          ENDIF
        END DO

        ! One integration cycle finished on the lowest grid level (coarsest
        ! resolution). Set model time.
        model_time_step => newTimedelta('+', 0, 0, 0, 0, 0, NINT(dtime), 0)
        this_datetime = this_datetime + model_time_step
        CALL deallocateTimedelta(model_time_step) 
          
        CALL output_ocean( patch_3d, &
          & ocean_state,             &
          & this_datetime,                &
          & ocean_surface,          &
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
    
    CALL timer_stop(timer_total)
    
  END SUBROUTINE ocean_test_zstar_advection
  !-------------------------------------------------------------------------



END MODULE mo_ocean_testbed_modules
