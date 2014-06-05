!>
!! Contains the main stepping routine the 3-dim hydrostatic ocean model.
!!
!! @author Peter Korn, Stephan Lorenz, MPI
!!
!! @par Revision History
!!  Initial version by Stephan Lorenz (MPI-M), (2010-04).
!!   - renaming and adjustment of hydrostatic ocean model V1.0.3 to ocean domain and patch_oce
!!  Modification by Stephan Lorenz, MPI-M, 2010-10
!!   - new module mo_ocean_output including updated reconstructions
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
MODULE mo_ocean_output
  !-------------------------------------------------------------------------
  USE mo_kind,                   ONLY: wp
  USE mo_impl_constants,         ONLY: max_char_length
  USE mo_model_domain,           ONLY: t_patch, t_patch_3d, t_subset_range, t_patch_vert
  USE mo_grid_config,            ONLY: n_dom
  USE mo_grid_subset,            ONLY: get_index_range
  USE mo_sync,                   ONLY: sync_patch_array, sync_e, sync_c !, sync_v
  USE mo_ocean_nml,              ONLY: iswm_oce, n_zlev, no_tracer, &
    & diagnostics_level, &
    & eos_type, i_sea_ice, gibraltar
  USE mo_dynamics_config,        ONLY: nold, nnew
  USE mo_io_config,              ONLY: n_checkpoints
  USE mo_run_config,             ONLY: nsteps, dtime, ltimer, output_mode
  USE mo_exception,              ONLY: message, message_text, finish
  USE mo_ext_data_types,         ONLY: t_external_data
  !USE mo_io_units,               ONLY: filename_max
  USE mo_datetime,               ONLY: t_datetime, print_datetime, add_time, datetime_to_string
  USE mo_oce_types,              ONLY: t_hydro_ocean_state, t_hydro_ocean_acc, t_hydro_ocean_diag, &
    & t_hydro_ocean_prog
  USE mo_oce_state,              ONLY: ocean_restart_list
  USE mo_operator_ocean_coeff_3d,ONLY: t_operator_coeff
  USE mo_oce_tracer,             ONLY: advect_tracer_ab
  USE mo_sea_ice,                ONLY: compute_mean_ice_statistics, reset_ice_statistics
  USE mo_sea_ice_types,          ONLY: t_sfc_flx, t_atmos_fluxes, t_atmos_for_ocean, &
    & t_sea_ice
  USE mo_oce_physics,            ONLY: t_ho_params
  USE mo_name_list_output,       ONLY: write_name_list_output, istime4name_list_output
  USE mo_oce_diagnostics,        ONLY: calc_slow_oce_diagnostics, calc_fast_oce_diagnostics, &
    & construct_oce_diagnostics,&
    & destruct_oce_diagnostics, t_oce_timeseries, &
    & calc_moc, calc_psi
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
  
  IMPLICIT NONE
  
  PRIVATE
  PUBLIC :: output_ocean
  
  INTEGER :: nsteps_since_last_output = 0
  
  !-------------------------------------------------------------------------
  
CONTAINS
  
  !-------------------------------------------------------------------------
  !>
  !! @par Revision History
  !! Developed by Peter Korn, MPI-M  (2008-2010).
  !! Initial release by Stephan Lorenz, MPI-M (2010-07)
  !
!<Optimize:inUse>
  SUBROUTINE output_ocean( patch_3d, ocean_state,     &
    & datetime,    &
    & p_sfc_flx,             &
    & p_ice, &
    & jstep, jstep0)

    TYPE(t_patch_3d ),TARGET, INTENT(inout)          :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET, INTENT(inout) :: ocean_state(n_dom)
    TYPE(t_datetime), INTENT(inout)                  :: datetime
    TYPE(t_sfc_flx)                                  :: p_sfc_flx
    TYPE (t_sea_ice),         INTENT(inout)          :: p_ice
    INTEGER,   INTENT(in)                            :: jstep, jstep0

    

    ! local variables
    INTEGER :: jg, jtrc
    INTEGER :: ocean_statistics
    !LOGICAL                         :: l_outputtime
    CHARACTER(LEN=32)               :: datestring, plaindatestring
!    TYPE(t_oce_timeseries), POINTER :: oce_ts
    TYPE(t_patch), POINTER :: patch_2d
    TYPE(t_patch_vert), POINTER :: patch_1d
    INTEGER, POINTER :: dolic(:,:)
    REAL(wp), POINTER :: prism_thickness(:,:,:)
    
    !CHARACTER(LEN=filename_max)  :: outputfile, gridfile
    CHARACTER(LEN=max_char_length), PARAMETER :: &
      & routine = 'mo_ocean_output:output_ocean'
    !------------------------------------------------------------------
    
    patch_2D      => patch_3d%p_patch_2d(1)
    jg = 1
    nsteps_since_last_output = nsteps_since_last_output + 1
    !------------------------------------------------------------------
     IF (.not. istime4name_list_output(jstep)) RETURN


    !------------------------------------------------------------------
        IF (istime4name_list_output(jstep))THEN!.OR.jstep>0) THEN
          IF (diagnostics_level == 1 ) THEN
            CALL calc_slow_oce_diagnostics( patch_3d       , &
              &                             ocean_state(jg), &
              &                             p_sfc_flx      , &
              &                             p_ice          , &
              &                             jstep-jstep0   , &
              &                             datetime) !    , &
            ! &                             oce_ts)
            IF (no_tracer>=2) THEN
              CALL calc_moc (patch_2d,patch_3d, ocean_state(jg)%p_diag%w(:,:,:), datetime)
            ENDIF
          ENDIF
          ! compute mean values for output interval
          !TODO [ram] src/io/shared/mo_output_event_types.f90 for types to use
          !TODO [ram] nsteps_since_last_output =
          !TODO [ram] output_event%event_step(output_event%i_event_step)%i_sim_step - output_event%event_step(output_event%i_event_step-1)%i_sim_step
        
          CALL compute_mean_ocean_statistics(ocean_state(1)%p_acc,p_sfc_flx,nsteps_since_last_output)
          CALL compute_mean_ice_statistics(p_ice%acc,nsteps_since_last_output)
        
          ! set the output variable pointer to the correct timelevel
          CALL set_output_pointers(nnew(1), ocean_state(jg)%p_diag, ocean_state(jg)%p_prog(nnew(1)))
        
          IF (output_mode%l_nml) THEN
            CALL write_name_list_output(jstep)
          ENDIF
        
          CALL message (TRIM(routine),'Write output at:')
          CALL print_datetime(datetime)
        
          ! reset accumulation vars
          CALL reset_ocean_statistics(ocean_state(1)%p_acc,p_sfc_flx,nsteps_since_last_output)
          IF (i_sea_ice >= 1) CALL reset_ice_statistics(p_ice%acc)
        
        END IF

  END SUBROUTINE output_ocean
  !-------------------------------------------------------------------------
    
  !-------------------------------------------------------------------------
!<Optimize:inUse>
  SUBROUTINE set_output_pointers(timelevel,p_diag,p_prog)
    INTEGER, INTENT(in) :: timelevel
    TYPE(t_hydro_ocean_diag) :: p_diag
    TYPE(t_hydro_ocean_prog) :: p_prog
    
    TYPE(t_list_element), POINTER :: output_var => NULL()
    TYPE(t_list_element), POINTER :: prog_var   => NULL()
    CHARACTER(LEN=max_char_length) :: timelevel_str
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
    IF(no_tracer>0)p_diag%t(:,:,:)        =  p_prog%tracer(:,:,:,1)
    
    !output_var             => find_list_element(ocean_restart_list,'s')
    !prog_var               => find_list_element(ocean_restart_list,'s'//TRIM(timelevel_str))
    !output_var%field%r_ptr => prog_var%field%r_ptr
    !p_diag%s               => prog_var%field%r_ptr(:,:,:,1,1)
     IF(no_tracer>1)p_diag%s(:,:,:)        =  p_prog%tracer(:,:,:,2)
  END SUBROUTINE set_output_pointers
  
END MODULE mo_ocean_output
