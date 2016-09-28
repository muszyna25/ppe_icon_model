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
!----------------------------
#include "icon_definitions.inc"
!----------------------------
MODULE mo_hydro_ocean_run
  !-------------------------------------------------------------------------
  USE mo_kind,                   ONLY: wp
  USE mo_impl_constants,         ONLY: max_char_length
  USE mo_model_domain,           ONLY: t_patch, t_patch_3d
  USE mo_grid_config,            ONLY: n_dom
  USE mo_ocean_nml,              ONLY: iswm_oce, n_zlev, no_tracer, lhamocc,&
    & i_sea_ice, cfl_check, cfl_threshold, cfl_stop_on_violation,   &
    & cfl_write, surface_module
  USE mo_ocean_nml,              ONLY: iforc_oce, Coupled_FluxFromAtmo
  USE mo_dynamics_config,        ONLY: nold, nnew
  USE mo_io_config,              ONLY: n_checkpoints, write_last_restart
  USE mo_run_config,             ONLY: nsteps, dtime, ltimer, output_mode, debug_check_level
  USE mo_exception,              ONLY: message, message_text, finish
  USE mo_ext_data_types,         ONLY: t_external_data
  !USE mo_io_units,               ONLY: filename_max
  !  USE mo_datetime,               ONLY: t_datetime, add_time, datetime_to_string
  USE mo_timer,                  ONLY: timer_start, timer_stop, timer_total, timer_solve_ab,  &
    & timer_tracer_ab, timer_vert_veloc, timer_normal_veloc,     &
    & timer_upd_phys, timer_upd_flx, timer_extra20, timers_level, &
    & timer_scalar_prod_veloc, timer_extra21, timer_extra22, timer_bgc_ini, &
    & timer_bgc_inv, timer_bgc_tot
  USE mo_ocean_ab_timestepping,    ONLY: solve_free_surface_eq_ab, &
    & calc_normal_velocity_ab,  &
    & calc_vert_velocity,       &
    & update_time_indices
  USE mo_ocean_types,              ONLY: t_hydro_ocean_state, &
    & t_operator_coeff, t_solvercoeff_singleprecision
  USE mo_ocean_math_operators,   ONLY: update_height_depdendent_variables, check_cfl_horizontal, check_cfl_vertical
  USE mo_scalar_product,         ONLY: calc_scalar_product_veloc_3d
  USE mo_ocean_tracer,             ONLY: advect_ocean_tracers
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
  USE mo_ocean_diagnostics,        ONLY: calc_fast_oce_diagnostics, calc_psi
  USE mo_ocean_ab_timestepping_mimetic, ONLY: construct_ho_lhs_fields_mimetic, destruct_ho_lhs_fields_mimetic
  USE mo_io_restart_attributes,  ONLY: get_restart_attribute
  USE mo_time_config,            ONLY: time_config
  USE mo_master_config,          ONLY: isRestart
!  USE mo_sea_ice_nml,            ONLY: i_ice_dyn
  USE mo_util_dbg_prnt,          ONLY: dbg_print, debug_printValue
  USE mo_dbg_nml,                ONLY: idbg_mxmn
  USE mo_statistics
  USE mo_ocean_statistics
  USE mo_hamocc_statistics,     ONLY: update_hamocc_statistics, reset_hamocc_statistics
  USE mo_hamocc_types,          ONLY: t_hamocc_state
  USE mo_derived_variable_handling, ONLY: perform_accumulation, reset_accumulation
  USE mo_ocean_output
  USE mo_ocean_coupling,         ONLY: couple_ocean_toatmo_fluxes

  USE mtime,                     ONLY: datetime, datetimeToString, deallocateDatetime,              &
       &                               timedelta, newTimedelta, deallocateTimedelta,                &
       &                               MAX_DATETIME_STR_LEN, newDatetime,                           &
       &                               MAX_MTIME_ERROR_STR_LEN, no_error, mtime_strerror,           &
       &                               OPERATOR(-), OPERATOR(+), OPERATOR(>), OPERATOR(*),          &
       &                               ASSIGNMENT(=), OPERATOR(==), OPERATOR(>=), OPERATOR(/=),     &
       &                               event, eventGroup, newEvent,                                 &
       &                               addEventToEventGroup, isCurrentEventActive
  USE mo_event_manager,          ONLY: initEventManager, addEventGroup, getEventGroup, printEventGroup
  
  USE mo_bgc_bcond,          ONLY: ext_data_bgc, update_bgc_bcond
  USE mo_hamocc_diagnostics,    ONLY: get_inventories
  USE mo_hamocc_nml,         ONLY: io_stdo_bgc


  IMPLICIT NONE

  PRIVATE

  PUBLIC  :: perform_ho_stepping
  PUBLIC  :: prepare_ho_stepping, end_ho_stepping
  PUBLIC  :: write_initial_ocean_timestep
  PUBLIC  :: update_time_g_n, update_time_indices
  
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
    if(lhamocc)then
      if(ltimer)call timer_start(timer_bgc_ini)
      CALL ini_bgc_icon(patch_3d,ocean_state,is_restart)
      if(ltimer)call timer_stop(timer_bgc_ini)
    endif
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
     CALL update_height_depdendent_variables( patch_3d, ocean_state, ext_data, operators_coefficients, solvercoeff_sp)
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
    & this_datetime,                                    &
    & surface_fluxes, p_sfc, p_phys_param,              &
    & p_as, p_atm_f, sea_ice, hamocc_state, operators_coefficients, &
    & solvercoeff_sp)

    TYPE(t_patch_3d ),TARGET, INTENT(inout)          :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET, INTENT(inout) :: ocean_state(n_dom)
    TYPE(t_external_data), TARGET, INTENT(in)        :: p_ext_data(n_dom)
    TYPE(datetime), POINTER                          :: this_datetime
    TYPE(t_sfc_flx)                                  :: surface_fluxes
    TYPE(t_ocean_surface)                            :: p_sfc
    TYPE (t_ho_params)                               :: p_phys_param
    TYPE(t_atmos_for_ocean),  INTENT(inout)          :: p_as
    TYPE(t_atmos_fluxes ),    INTENT(inout)          :: p_atm_f
    TYPE (t_sea_ice),         INTENT(inout)          :: sea_ice
    TYPE(t_hamocc_state), INTENT(INOUT)                :: hamocc_state
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

    TYPE(eventGroup), POINTER           :: checkpointEventGroup => NULL()
    TYPE(timedelta), POINTER            :: model_time_step => NULL()
    TYPE(datetime), POINTER             :: mtime_current   => NULL()
    TYPE(datetime), POINTER             :: eventRefDate    => NULL(), eventStartDate  => NULL(), eventEndDate    => NULL()
    TYPE(timedelta), POINTER            :: eventInterval   => NULL()
    TYPE(event), POINTER                :: checkpointEvent => NULL()
    TYPE(event), POINTER                :: restartEvent    => NULL()
    
    INTEGER                             :: checkpointEvents, ierr
    LOGICAL                             :: lwrite_checkpoint, lret

    CHARACTER(LEN=MAX_DATETIME_STR_LEN)    :: dstring
    CHARACTER(len=MAX_MTIME_ERROR_STR_LEN) :: errstring

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

    ! CALL datetime_to_string(datestring, this_datetime)

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

    ! set events, group and the events

    CALL message('','')

    eventRefDate   => time_config%tc_exp_refdate
    eventStartDate => time_config%tc_exp_startdate
    eventEndDate   => time_config%tc_exp_stopdate

    ! use start/end setup from the restart
    IF (isRestart() .AND. time_config%is_relative_time) THEN
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

    jstep = jstep0
    TIME_LOOP: DO 
    IF (lhamocc) THEN
      IF (ltimer) CALL timer_start(timer_bgc_inv)
      CALL message ('start of time loop', 'HAMOCC inventories', io_stdo_bgc)
      IF (ltimer) CALL get_inventories(hamocc_state,ocean_state(1),patch_3d,nold(1))
      CALL timer_stop(timer_bgc_inv)
    ENDIF

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
            
      start_detail_timer(timer_extra22,6)
      CALL update_height_depdendent_variables( patch_3d, ocean_state(jg), p_ext_data(jg), operators_coefficients, solvercoeff_sp)
      stop_detail_timer(timer_extra22,6)
      
      start_timer(timer_scalar_prod_veloc,2)
      CALL calc_scalar_product_veloc_3d( patch_3d,  &
        & ocean_state(jg)%p_prog(nold(1))%vn,         &
        & ocean_state(jg)%p_diag,                     &
        & operators_coefficients)
      stop_timer(timer_scalar_prod_veloc,2)
      
      !In case of a time-varying forcing:
      ! update_surface_flux or update_ocean_surface has changed p_prog(nold(1))%h, SST and SSS
      start_timer(timer_upd_flx,3)
      IF (surface_module == 1) THEN
        CALL update_surface_flux( patch_3d, ocean_state(jg), p_as, sea_ice, p_atm_f, surface_fluxes, &
          & jstep, mtime_current, operators_coefficients)
      ELSEIF (surface_module == 2) THEN
        CALL update_ocean_surface( patch_3d, ocean_state(jg), p_as, sea_ice, p_atm_f, surface_fluxes, p_sfc, &
          & jstep, mtime_current, operators_coefficients)
      ENDIF

    


      IF(lhamocc)CALL update_bgc_bcond( patch_3d, ext_data_bgc, jstep, this_datetime)
      stop_timer(timer_upd_flx,3)

      start_detail_timer(timer_extra22,4)
      CALL update_height_depdendent_variables( patch_3d, ocean_state(jg), p_ext_data(jg), operators_coefficients, solvercoeff_sp)
      stop_detail_timer(timer_extra22,4)

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
      start_timer(timer_solve_ab,1)
      CALL solve_free_surface_eq_ab (patch_3d, ocean_state(jg), p_ext_data(jg), &
        & surface_fluxes, p_phys_param, jstep, operators_coefficients, solvercoeff_sp, return_status)!, p_int(jg))
      IF (return_status /= 0) THEN
       CALL output_ocean(              &
         & patch_3d=patch_3d,          &
         & ocean_state=ocean_state,    &
         & this_datetime=mtime_current, &
         & surface_fluxes=surface_fluxes, &
         & sea_ice=sea_ice,            &
         & hamocc=hamocc_state,        &
         & jstep=jstep, jstep0=jstep0, &
         & force_output=.true.)
        CALL finish(TRIM(routine), 'solve_free_surface_eq_ab  returned error')
      ENDIF
      
      stop_timer(timer_solve_ab,1)

      !------------------------------------------------------------------------
      ! Step 4: calculate final normal velocity from predicted horizontal
      ! velocity vn_pred and updated surface height
      start_timer(timer_normal_veloc,4)
      CALL calc_normal_velocity_ab(patch_3d, ocean_state(jg),&
        & operators_coefficients, solvercoeff_sp,  p_ext_data(jg), p_phys_param)
      stop_timer(timer_normal_veloc,4)

      !------------------------------------------------------------------------
      ! Step 5: calculate vertical velocity from continuity equation under
      ! incompressiblity condition in the non-shallow-water case
      IF ( iswm_oce /= 1 ) THEN
        start_timer(timer_vert_veloc,4)
        CALL calc_vert_velocity( patch_3d, ocean_state(jg),operators_coefficients)
        stop_timer(timer_vert_veloc,4)
      ENDIF
      !------------------------------------------------------------------------

      IF (idbg_mxmn >= 2 .OR. debug_check_level > 5) THEN
        CALL horizontal_mean(values=ocean_state(jg)%p_prog(nnew(1))%h(:,:), weights=patch_2d%cells%area(:,:), &
          & in_subset=patch_2d%cells%owned, mean=mean_height)
        CALL debug_printValue(description="Mean Height", value=mean_height, detail_level=2)
      ENDIF
      IF (debug_check_level > 5 .AND. idbg_mxmn >= 2) THEN
        ! check difference from old_mean_height
!         CALL debug_printValue(description="Old/New Mean Height", value=old_mean_height, &
!           & value1=mean_height, value2=(mean_height-old_mean_height) / old_mean_height, &
!           & detail_level=2)
        CALL debug_printValue(description="Old/New Mean Height", &
          & value=old_mean_height, value1=mean_height, detail_level=2)
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


         ! Step : call HAMOCC
      if(lhamocc)then
        if(ltimer) call timer_start(timer_bgc_tot)
        CALL bgc_icon(patch_3d,ocean_state(jg),p_as,sea_ice)
        if(ltimer) call timer_stop(timer_bgc_tot)
      endif
      !------------------------------------------------------------------------
      ! Step 6: transport tracers and diffuse them
      IF (no_tracer>=1) THEN
        start_timer(timer_tracer_ab,1)
        CALL advect_ocean_tracers( patch_3d, ocean_state(jg), p_phys_param,&
          & surface_fluxes,&
          & operators_coefficients,&
          & jstep)
        stop_timer(timer_tracer_ab,1)
      ENDIF

      ! perform accumulation for special variables
      start_detail_timer(timer_extra20,5)     
      IF (no_tracer>=1) THEN
        CALL calc_potential_density( patch_3d,                            &
          & ocean_state(jg)%p_prog(nold(1))%tracer,                       &
          & ocean_state(jg)%p_diag%rhopot )
          
        ! calculate diagnostic barotropic stream function
        CALL calc_psi (patch_3d, ocean_state(jg)%p_diag%u(:,:,:),         &
          & patch_3D%p_patch_1d(1)%prism_thick_c(:,:,:),                  &
          & ocean_state(jg)%p_diag%u_vint, mtime_current)
        CALL dbg_print('calc_psi: u_vint' ,ocean_state(jg)%p_diag%u_vint, str_module, 3, in_subset=patch_2d%cells%owned)
          
        ! calculate diagnostic barotropic stream function with vn
    !  not yet mature
    !   CALL calc_psi_vn (patch_3d, ocean_state(jg)%p_prog(nold(1))%vn,   &
    !     & patch_3D%p_patch_1d(1)%prism_thick_e(:,:,:),                  &
    !     & operators_coefficients,                                       &
    !     & ocean_state(jg)%p_diag%u_vint, ocean_state(jg)%p_diag%v_vint, mtime_current)
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

      IF(lhamocc)CALL update_hamocc_statistics(hamocc_state,&
        & patch_2d%cells%owned,&
        & patch_2d%edges%owned,&
        & patch_2d%verts%owned,&
        & n_zlev)

      CALL calc_fast_oce_diagnostics( patch_2d,      &
        & patch_3d%p_patch_1d(1)%dolic_c, &
        & patch_3d%p_patch_1d(1)%prism_thick_c, &
        & patch_3d%p_patch_1d(1)%zlev_m, &
        & ocean_state(jg)%p_diag)

      stop_detail_timer(timer_extra20,5)

      CALL perform_accumulation(nnew(1),0)

      CALL output_ocean( patch_3d, ocean_state, &
        &                mtime_current,              &
        &                surface_fluxes,             &
        &                sea_ice,                 &
        &                hamocc_state,            &
        &                jstep, jstep0)
      
    CALL reset_accumulation
      ! send and receive coupling fluxes for ocean at the end of time stepping loop
      IF (iforc_oce == Coupled_FluxFromAtmo) &  !  14
        &  CALL couple_ocean_toatmo_fluxes(patch_3D, ocean_state(jg), sea_ice, p_atm_f, p_as, mtime_current)
!       &  CALL couple_ocean_toatmo_fluxes(patch_3D, ocean_state(jg), sea_ice, p_atm_f, p_as%fu10, mtime_current)
!       &  CALL couple_ocean_toatmo_fluxes(patch_3D, ocean_state(jg), sea_ice, p_atm_f, mtime_current)

  
      start_detail_timer(timer_extra21,5)
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
        CALL create_restart_file( patch = patch_2d,       &
             & current_date=mtime_current, &
             & jstep=jstep,            &
             & model_type="oce",       &
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
      
      stop_detail_timer(timer_extra21,5)

    ENDDO time_loop
    
    IF(lhamocc) THEN
     if(ltimer)CALL timer_start(timer_bgc_inv)
     CALL message ('end of time loop', 'HAMOCC inventories', io_stdo_bgc)
     CALL get_inventories(hamocc_state,ocean_state(1),patch_3d,nold(1))
     if(ltimer)CALL timer_stop(timer_bgc_inv)
    ENDIF

    IF (write_last_restart .and. .not. lwrite_checkpoint) &
        & CALL create_restart_file( patch = patch_2d,       &
        & current_date=mtime_current, &
        & jstep=jstep-1,          &
        & model_type="oce",       &
        & opt_nice_class=1,       &
        & ocean_zlevels=n_zlev,                                         &
        & ocean_zheight_cellmiddle = patch_3d%p_patch_1d(1)%zlev_m(:),  &
        & ocean_zheight_cellinterfaces = patch_3d%p_patch_1d(1)%zlev_i(:))
  
    CALL timer_stop(timer_total)
  
  END SUBROUTINE perform_ho_stepping
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
!<Optimize:inUse>
  SUBROUTINE write_initial_ocean_timestep(patch_3d,ocean_state,surface_fluxes,sea_ice, &
& hamocc_state, operators_coefficients, p_phys_param)

    TYPE(t_patch_3D), INTENT(IN) :: patch_3d
    TYPE(t_hydro_ocean_state), INTENT(INOUT)    :: ocean_state
    TYPE(t_sfc_flx) , INTENT(INOUT)             :: surface_fluxes
    TYPE(t_sea_ice),          INTENT(INOUT)     :: sea_ice
    TYPE(t_hamocc_state), INTENT(INOUT)         :: hamocc_state
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
      IF (lhamocc) CALL update_hamocc_statistics( &
      & hamocc_state,            &
      & patch_2d%cells%owned,   &
      & patch_2d%edges%owned,   &
      & patch_2d%verts%owned,   &
      & n_zlev)

    CALL perform_accumulation(nnew(1),0)

    CALL write_name_list_output(jstep=0)

    CALL reset_ocean_statistics(ocean_state%p_acc,ocean_state%p_diag,surface_fluxes)
    CALL reset_accumulation
    IF (i_sea_ice >= 1) CALL reset_ice_statistics(sea_ice%acc)
    IF (lhamocc) CALL reset_hamocc_statistics(hamocc_state%p_acc)

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
