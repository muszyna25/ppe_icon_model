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
USE mo_grid_config,            ONLY: n_dom
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
USE mo_oce_state,              ONLY: t_hydro_ocean_state, t_hydro_ocean_acc, &
  &                                  init_ho_base, init_ho_basins, v_base, &
  &                                  construct_hydro_ocean_base, &! destruct_hydro_ocean_base, &
  &                                  construct_hydro_ocean_state, destruct_hydro_ocean_state, &
  &                                  init_coriolis_oce, init_oce_config, &
  &                                  set_lateral_boundary_values, construct_patch_3D, init_patch_3D, &
  &                                  setup_ocean_namelists, check_ocean_subsets, &
  &                                  ocean_default_list, ocean_restart_list
USE mo_oce_math_operators,     ONLY: calc_thickness! , height_related_quantities
USE mo_operator_ocean_coeff_3d,ONLY: t_operator_coeff, allocate_exp_coeff,par_init_operator_coeff,&
  &                                  update_diffusion_matrices
USE mo_scalar_product,         ONLY: calc_scalar_product_veloc_3D
USE mo_oce_tracer,             ONLY: advect_tracer_ab
USE mo_io_restart,             ONLY: write_restart_info_file
USE mo_oce_bulk,               ONLY: update_sfcflx
USE mo_sea_ice,                ONLY: construct_sfcflx,destruct_sfcflx,&
  &                                  construct_atmos_for_ocean,&
  &                                  destruct_atmos_for_ocean,&
  &                                  construct_atmos_fluxes, destruct_atmos_fluxes,&
  &                                  construct_sea_ice, destruct_sea_ice, ice_init
USE mo_sea_ice_types,          ONLY: t_sfc_flx, t_atmos_fluxes, t_atmos_for_ocean, &
  &                                  t_sea_ice
USE mo_oce_forcing,            ONLY: init_sfcflx
USE mo_oce_physics,            ONLY: t_ho_params, &
  &                                  construct_ho_params, init_ho_params, &
  &                                  destruct_ho_params, update_ho_params
USE mo_oce_thermodyn,          ONLY: calc_density_MPIOM_func, calc_density_lin_EOS_func,&
  &                                  calc_density_JMDWFG06_EOS_func, calc_potential_density
USE mo_output,                 ONLY: init_output_files, &
  &                                  create_restart_file, write_output_oce! , write_output
USE mo_fortran_tools,          ONLY: assign_if_present
USE mo_name_list_output,       ONLY: write_name_list_output, istime4name_list_output
USE mo_oce_diagnostics,        ONLY: calculate_oce_diagnostics,&
  &                                  construct_oce_diagnostics,&
  &                                  destruct_oce_diagnostics, t_oce_timeseries, &
  &                                  calc_moc, calc_psi
USE mo_oce_ab_timestepping_mimetic, ONLY: init_ho_lhs_fields_mimetic
USE mo_linked_list,            ONLY: t_list_element, find_list_element
USE mo_var_list,               ONLY: print_var_list
  USE mo_mpi,                               ONLY: my_process_is_stdio
  USE mo_time_config,          ONLY: time_config
  USE mo_master_control,       ONLY: is_restart_run
  USE mo_statistics




IMPLICIT NONE

PRIVATE
INTEGER, PARAMETER :: kice = 1

!VERSION CONTROL:
CHARACTER(LEN=*), PARAMETER :: version = '$Id$'

! public interface
!
! public subroutines
PUBLIC :: perform_ho_stepping
PUBLIC :: prepare_ho_integration
PUBLIC :: finalise_ho_integration
PRIVATE:: update_intermediate_tracer_vars
!
INTERFACE add_fields
  MODULE PROCEDURE add_fields_3d
  MODULE PROCEDURE add_fields_2d
END INTERFACE add_fields
!
!-------------------------------------------------------------------------

CONTAINS

  !-------------------------------------------------------------------------
  !>
  !! Main stepping routine for call of hydrostatic ocean model
  !!
  !! @par Revision History
  !! Developed by Peter Korn, MPI-M  (2008-2010).
  !! Initial release by Stephan Lorenz, MPI-M (2010-07)
  !
  !
  SUBROUTINE perform_ho_stepping( p_patch_3D, p_os, p_ext_data,          &
                                & datetime, n_io, jfile, lwrite_restart, &
                                & p_sfc_flx, p_phys_param,               &
                                & p_as, p_atm_f, p_ice,p_op_coeff,       &
                                & l_have_output)

  TYPE(t_patch_3D ),TARGET, INTENT(INOUT)          :: p_patch_3D
  TYPE(t_hydro_ocean_state), TARGET, INTENT(INOUT) :: p_os(n_dom)
  TYPE(t_external_data), TARGET, INTENT(IN)        :: p_ext_data(n_dom)
  TYPE(t_datetime), INTENT(INOUT)                  :: datetime
  INTEGER, INTENT(IN)                              :: n_io
  INTEGER, INTENT(INOUT)                           :: jfile
  LOGICAL, INTENT(IN)                              :: lwrite_restart
  TYPE(t_sfc_flx)                                  :: p_sfc_flx
  TYPE (t_ho_params)                               :: p_phys_param
  TYPE(t_atmos_for_ocean),  INTENT(INOUT)          :: p_as
  TYPE(t_atmos_fluxes ),    INTENT(INOUT)          :: p_atm_f
  TYPE (t_sea_ice),         INTENT(INOUT)          :: p_ice
  TYPE(t_operator_coeff),   INTENT(INOUT)          :: p_op_coeff
  LOGICAL,                  INTENT(INOUT)          :: l_have_output

  ! local variables
  INTEGER                         :: jstep, jg, jtrc
  INTEGER                         :: nsteps_since_last_output
  INTEGER                         :: ocean_statistics
  !LOGICAL                         :: l_outputtime
  CHARACTER(len=32)               :: datestring, plaindatestring
  TYPE(t_oce_timeseries), POINTER :: oce_ts
  TYPE(t_patch), POINTER          :: p_patch

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

  p_patch => p_patch_3D%p_patch_2D(jg)

  CALL datetime_to_string(datestring, datetime)

  IF (idiag_oce == 1) &
    & CALL construct_oce_diagnostics( p_patch_3D, p_os(jg), oce_ts, datestring)

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
    IF (.NOT. is_restart_run()) p_os(jg)%p_prog(nnew(1))%tracer = p_os(jg)%p_prog(nold(1))%tracer
    CALL write_name_list_output( datetime, time_config%sim_time(1), last_step=.FALSE., initial_step=.not.is_restart_run())
  ENDIF
  !------------------------------------------------------------------
  ! call the dynamical core: start the time loop
  !------------------------------------------------------------------
  TIME_LOOP: DO jstep = 1, nsteps

    CALL datetime_to_string(datestring, datetime)
    WRITE(message_text,'(a,i6,2a)') '  Begin of timestep =',jstep,'  datetime:  ', datestring
    CALL message (TRIM(routine), message_text)

    IF(itestcase_oce==28)THEN
      CALL calc_thickness( p_patch_3D, p_os(jg), p_ext_data(jg))
      !CALL height_related_quantities( p_patch_3D, p_os(jg), p_ext_data(jg))
      CALL calc_vert_velocity(p_patch_3D, p_os(jg),p_op_coeff)
      CALL advect_tracer_ab( p_patch_3D, p_os(jg),  &
                           & p_phys_param,p_sfc_flx,&
                           & p_op_coeff,&
                           & jstep)
    ELSE
      !In case of a time-varying forcing:
      IF (ltimer) CALL timer_start(timer_upd_flx)
      CALL update_sfcflx( p_patch_3D, p_os(jg), p_as, p_ice, p_atm_f, p_sfc_flx, &
        &                jstep, datetime)

      IF(.NOT.l_STAGGERED_TIMESTEP)THEN

        CALL calc_thickness( p_patch_3D, p_os(jg), p_ext_data(jg))
        !CALL height_related_quantities( p_patch_3D, p_os(jg), p_ext_data(jg))

        CALL set_lateral_boundary_values( p_patch_3D, p_os(jg)%p_prog(nold(1))%vn)
        CALL sync_patch_array(sync_e,  p_patch_3D%p_patch_2D(jg),  p_os(jg)%p_prog(nold(1))%vn)

        CALL calc_scalar_product_veloc_3D( p_patch_3D, &
          & p_os(jg)%p_prog(nold(1))%vn,         &
          & p_os(jg)%p_prog(nold(1))%vn,         &
          & p_os(jg)%p_diag,                     &
          & p_op_coeff)

      ENDIF
      IF (ltimer) CALL timer_stop(timer_upd_flx)

      IF (ltimer) CALL timer_start(timer_upd_phys)

      SELECT CASE (EOS_TYPE)
        CASE(1)
          CALL update_ho_params(p_patch_3D, p_os(jg), p_sfc_flx, p_phys_param,&
            &                   calc_density_lin_EOS_func)
        CASE(2)
          CALL update_ho_params(p_patch_3D, p_os(jg), p_sfc_flx, p_phys_param,&
            &                   calc_density_MPIOM_func)
        CASE(3)
          CALL update_ho_params(p_patch_3D,p_os(jg), p_sfc_flx, p_phys_param,&
            &                   calc_density_JMDWFG06_EOS_func)
      CASE DEFAULT
      END SELECT
      IF (ltimer) CALL timer_stop(timer_upd_phys)

      CALL update_diffusion_matrices( p_patch_3D,                   &
                                    & p_os(jg),                     &
                                    & p_phys_param,                 &
                                    & p_op_coeff%matrix_vert_diff_e,&
                                    & p_op_coeff%matrix_vert_diff_c)

      ! solve for new free surface
      IF (ltimer) CALL timer_start(timer_solve_ab)
      CALL solve_free_surface_eq_ab (p_patch_3D, p_os(jg), p_ext_data(jg), &
        &                            p_sfc_flx, p_phys_param, jstep, p_op_coeff)!, p_int(jg))
      IF (ltimer) CALL timer_stop(timer_solve_ab)

      ! Step 4: calculate final normal velocity from predicted horizontal velocity vn_pred
      !         and updated surface height
      IF (ltimer) CALL timer_start(timer_normal_veloc)
      CALL calc_normal_velocity_ab(p_patch_3D, p_os(jg),&
                                  &p_op_coeff, p_ext_data(jg), p_phys_param)
      IF (ltimer) CALL timer_stop(timer_normal_veloc)

      ! Step 5: calculate vertical velocity from continuity equation under incompressiblity condition
      ! in the non-shallow-water case
      IF ( iswm_oce /= 1 ) THEN
        IF (ltimer) CALL timer_start(timer_vert_veloc)
        CALL calc_vert_velocity( p_patch_3D, p_os(jg),p_op_coeff)
        IF (ltimer) CALL timer_stop(timer_vert_veloc)
      ENDIF

      ! Step 6: transport tracers and diffuse them
      IF (no_tracer>=1) THEN
        IF (ltimer) CALL timer_start(timer_tracer_ab)
        CALL advect_tracer_ab( p_patch_3D, p_os(jg), p_phys_param,&
                             & p_sfc_flx,&
                             & p_op_coeff,&
                             & jstep)
        IF (ltimer) CALL timer_stop(timer_tracer_ab)
      ENDIF

    ENDIF  ! testcase 28

    ! One integration cycle finished on the lowest grid level (coarsest
    ! resolution). Set model time.
    CALL add_time(dtime,0,0,0,datetime)

    ! Not nice, but the name list output requires this
    time_config%sim_time(1) = time_config%sim_time(1) + dtime

    ! perform accumulation for special variables
    CALL calc_potential_density( p_patch_3D,                     &
      &                          p_os(jg)%p_prog(nold(1))%tracer,&
      &                          p_os(jg)%p_diag%rhopot )

    ! update accumulated vars
    CALL update_ocean_statistics(p_os(1),p_sfc_flx,p_patch_3D%p_patch_2D(1)%cells%owned)

    IF (is_output_time(jstep) .OR. istime4name_list_output(time_config%sim_time(1))) THEN
      IF (idiag_oce == 1 ) THEN
        CALL calculate_oce_diagnostics( p_patch_3D,    &
          &                             p_os(jg),      &
          &                             p_sfc_flx,     &
          &                             p_ice,         &
          &                             p_phys_param,  &
          &                             jstep,         &
          &                             datetime,      &
          &                             oce_ts)

        CALL calc_moc (p_patch,p_patch_3D, p_os(jg)%p_diag%w(:,:,:), datetime)
        CALL calc_psi (p_patch,p_patch_3D, p_os(jg)%p_diag%u(:,:,:), &
          &                        p_os(jg)%p_prog(nold(1))%h(:,:), &
          &                        p_os(jg)%p_diag%u_vint, datetime)

      ENDIF
      ! compute mean values for output interval
      CALL compute_mean_ocean_statistics(p_os(1)%p_acc,p_sfc_flx,nsteps_since_last_output)

      ! set the output variable pointer to the correct timelevel
!TODO (ram)      CALL set_output_pointers

      IF (output_mode%l_nml) THEN
        CALL write_name_list_output( datetime, time_config%sim_time(1), jstep==nsteps)
      ENDIF
      IF (output_mode%l_vlist) THEN
          CALL write_output_oce( datetime, time_config%sim_time(1),p_patch_3D, p_os)
      ENDIF

      CALL message (TRIM(routine),'Write output at:')
      CALL print_datetime(datetime)
      l_have_output = .TRUE.

      ! reset accumulation vars
      CALL reset_ocean_statistics(p_os(1)%p_acc,p_sfc_flx,nsteps_since_last_output)

    END IF

    ! If it's time, close the current output file and trigger a new one
    IF (jstep/=1.AND.(MOD(jstep,n_files())==0).AND.jstep/=nsteps  .AND. output_mode%l_vlist ) THEN
      jfile = jfile +1
      CALL init_output_files(jfile,lclose=l_have_output,p_patch_2D=p_patch_3D%p_patch_2D)
    ENDIF

!   ! close the current output file and trigger a new one
!   IF (istime4newoutputfile(jstep)) THEN
!     jfile = jfile +1
!     CALL init_output_files(jfile,lclose=l_have_output,p_patch_2D=p_patch_3D%p_patch_2D)
!   END IF

    ! Shift time indices for the next loop
    ! this HAS to ge into the restart files, because the start with the following loop
    CALL update_time_indices(jg)
    ! update intermediate timestepping variables for the tracers
    CALL update_intermediate_tracer_vars(p_os(jg))

    ! write a restart or checkpoint file
    IF (MOD(jstep,n_checkpoints())==0 .OR. (jstep==nsteps .AND. lwrite_restart)) THEN
      CALL create_restart_file( p_patch, datetime,                 &
                              & jfile, l_have_output,opt_depth=n_zlev,&
                              & opt_sim_time=time_config%sim_time(1), &
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
  !
  !
  SUBROUTINE prepare_ho_integration(p_patch_3D, p_os, p_ext_data, p_sfc_flx, &
                                  & p_phys_param, p_as,&
                                  & p_atm_f, p_ice, p_op_coeff)

    TYPE(t_patch_3D ),TARGET,   INTENT(INOUT)  :: p_patch_3D
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
      &      routine = 'mo_test_hydro_ocean:prepare_ho_integration'

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
    CALL init_dbg_index(p_patch_3D%p_patch_2D(jg))!(p_patch(jg))

    ! hydro_ocean_base contains the 3-dimensional structures for the ocean state

    CALL construct_patch_3D(p_patch_3D)

    CALL construct_hydro_ocean_base(p_patch_3D%p_patch_2D(jg), v_base)
    CALL init_ho_base     (p_patch_3D%p_patch_2D(jg), p_ext_data(jg), v_base)
    CALL init_ho_basins   (p_patch_3D%p_patch_2D(jg),                 v_base)
    CALL init_coriolis_oce(p_patch_3D%p_patch_2D(jg) )
    CALL init_patch_3D    (p_patch_3D,                p_ext_data(jg), v_base)
    !CALL init_patch_3D(p_patch_3D, v_base)
    CALL check_ocean_subsets(p_patch_3D)

    !------------------------------------------------------------------
    ! construct ocean state and physics
    !------------------------------------------------------------------

    ! p_patch and p_os have dimension n_dom
    CALL construct_hydro_ocean_state(p_patch_3D%p_patch_2D, p_os)

    ! initialize ocean indices for debug output (including 3-dim lsm)
    CALL init_oce_index( p_patch_3D%p_patch_2D,p_patch_3D, p_os, p_ext_data )

    CALL construct_ho_params(p_patch_3D%p_patch_2D(jg), p_phys_param)
    CALL init_ho_params(p_patch_3D, p_phys_param)

    !------------------------------------------------------------------
    ! construct ocean forcing and testcases
    !------------------------------------------------------------------

    CALL construct_sfcflx(p_patch_3D%p_patch_2D(jg),p_sfc_flx, ocean_default_list)
    CALL      init_sfcflx(p_patch_3D, p_sfc_flx)

    CALL construct_sea_ice(p_patch_3D%p_patch_2D(jg), p_ice, kice)
    CALL construct_atmos_for_ocean(p_patch_3D%p_patch_2D(jg), p_as)
    CALL construct_atmos_fluxes(p_patch_3D%p_patch_2D(jg), p_atm_f, kice)

    IF (init_oce_prog == 0) THEN
      CALL init_ho_testcases(p_patch_3D%p_patch_2D(jg),p_patch_3D, p_os(jg), p_ext_data(jg), p_op_coeff,p_sfc_flx)
    ELSE IF (init_oce_prog == 1) THEN

      CALL init_ho_prog(p_patch_3D%p_patch_2D(jg),p_patch_3D, p_os(jg), p_sfc_flx)
    END IF

    IF (init_oce_relax == 1) THEN
      CALL init_ho_relaxation(p_patch_3D%p_patch_2D(jg),p_patch_3D, p_os(jg), p_sfc_flx)
    END IF

    CALL init_ho_coupled(p_patch_3D%p_patch_2D(jg), p_os(jg))
    IF (i_sea_ice >= 1) &
      &   CALL ice_init(p_patch_3D%p_patch_2D(jg), p_os(jg), p_ice)

    CALL allocate_exp_coeff     ( p_patch_3D%p_patch_2D(jg), p_op_coeff)
    CALL par_init_operator_coeff( p_patch_3D, p_os(jg),p_phys_param, p_op_coeff)
    CALL init_ho_recon_fields   ( p_patch_3D%p_patch_2D(jg),p_patch_3D, p_os(jg), p_op_coeff)

    CALL init_ho_lhs_fields_mimetic   ( p_patch_3D )

    CALL message (TRIM(routine),'end')

  END SUBROUTINE prepare_ho_integration

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

  SUBROUTINE update_ocean_statistics(p_os,p_sfc_flx,subset)
    TYPE(t_hydro_ocean_state), INTENT(INOUT) :: p_os
    TYPE(t_sfc_flx),           INTENT(INOUT) :: p_sfc_flx
    TYPE(t_subset_range),INTENT(IN) :: subset

    INTEGER :: jtrc,i


    ! update ocean state accumulated values
    CALL add_fields(p_os%p_acc%h, p_os%p_prog(nnew(1))%h, subset)
    CALL add_fields(p_os%p_acc%u, p_os%p_diag%u, subset)
    CALL add_fields(p_os%p_acc%v, p_os%p_diag%v, subset)
    CALL add_fields(p_os%p_acc%rhopot,p_os%p_diag%rhopot,subset)
    DO jtrc=1,no_tracer
    CALL add_fields(p_os%p_acc%tracer(:,:,:,jtrc), &
      &                  p_os%p_prog(nnew(1))%tracer(:,:,:,jtrc), &
      &                  subset)
    END DO
    CALL add_fields(p_os%p_acc%u_vint        , p_os%p_diag%u_vint        , subset)
    CALL add_fields(p_os%p_acc%w             , p_os%p_diag%w             , subset,n_zlev+1)
    CALL add_fields(p_os%p_acc%div_mass_flx_c, p_os%p_diag%div_mass_flx_c, subset)

    ! update forcing accumulated values
    CALL add_fields(p_sfc_flx%forc_wind_u_acc  , p_sfc_flx%forc_wind_u  , subset)
    CALL add_fields(p_sfc_flx%forc_wind_v_acc  , p_sfc_flx%forc_wind_v  , subset)
    CALL add_fields(p_sfc_flx%forc_swflx_acc   , p_sfc_flx%forc_swflx   , subset)
    CALL add_fields(p_sfc_flx%forc_lwflx_acc   , p_sfc_flx%forc_lwflx   , subset)
    CALL add_fields(p_sfc_flx%forc_ssflx_acc   , p_sfc_flx%forc_ssflx   , subset)
    CALL add_fields(p_sfc_flx%forc_slflx_acc   , p_sfc_flx%forc_slflx   , subset)
    CALL add_fields(p_sfc_flx%forc_precip_acc  , p_sfc_flx%forc_precip  , subset)
    CALL add_fields(p_sfc_flx%forc_evap_acc    , p_sfc_flx%forc_evap    , subset)
    CALL add_fields(p_sfc_flx%forc_runoff_acc  , p_sfc_flx%forc_runoff  , subset)
    CALL add_fields(p_sfc_flx%forc_fwbc_acc    , p_sfc_flx%forc_fwbc    , subset)
    CALL add_fields(p_sfc_flx%forc_fwrelax_acc , p_sfc_flx%forc_fwrelax , subset)
    CALL add_fields(p_sfc_flx%forc_fwfx_acc    , p_sfc_flx%forc_fwfx    , subset)
    CALL add_fields(p_sfc_flx%forc_hfrelax_acc , p_sfc_flx%forc_hfrelax , subset)
    CALL add_fields(p_sfc_flx%forc_hflx_acc    , p_sfc_flx%forc_hflx    , subset)
    DO jtrc=1,no_tracer
      CALL add_fields(p_sfc_flx%forc_tracer_acc(:,:,jtrc), p_sfc_flx%forc_tracer(:,:,jtrc), subset)
      CALL add_fields(p_sfc_flx%forc_tracer_relax_acc(:,:,jtrc), p_sfc_flx%forc_tracer_relax(:,:,jtrc), subset)
    END DO
  END SUBROUTINE update_ocean_statistics

  SUBROUTINE compute_mean_ocean_statistics(p_acc,p_sfc_flx,nsteps_since_last_output)
    TYPE(t_hydro_ocean_acc), INTENT(INOUT) :: p_acc
    TYPE(t_sfc_flx),         INTENT(INOUT) :: p_sfc_flx
    INTEGER,INTENT(IN)                     :: nsteps_since_last_output

    p_acc%tracer                    = p_acc%tracer                   /REAL(nsteps_since_last_output,wp)
    p_acc%h                         = p_acc%h                        /REAL(nsteps_since_last_output,wp)
    p_acc%u                         = p_acc%u                        /REAL(nsteps_since_last_output,wp)
    p_acc%v                         = p_acc%v                        /REAL(nsteps_since_last_output,wp)
    p_acc%rhopot                    = p_acc%rhopot                   /REAL(nsteps_since_last_output,wp)
    p_acc%u_vint                    = p_acc%u_vint                   /REAL(nsteps_since_last_output,wp)
    p_acc%w                         = p_acc%w                        /REAL(nsteps_since_last_output,wp)
    p_acc%div_mass_flx_c            = p_acc%div_mass_flx_c           /REAL(nsteps_since_last_output,wp)
    p_sfc_flx%forc_wind_u_acc       = p_sfc_flx%forc_wind_u_acc      /REAL(nsteps_since_last_output,wp)
    p_sfc_flx%forc_wind_v_acc       = p_sfc_flx%forc_wind_v_acc      /REAL(nsteps_since_last_output,wp)
    p_sfc_flx%forc_swflx_acc        = p_sfc_flx%forc_swflx_acc       /REAL(nsteps_since_last_output,wp)
    p_sfc_flx%forc_lwflx_acc        = p_sfc_flx%forc_lwflx_acc       /REAL(nsteps_since_last_output,wp)
    p_sfc_flx%forc_ssflx_acc        = p_sfc_flx%forc_ssflx_acc       /REAL(nsteps_since_last_output,wp)
    p_sfc_flx%forc_slflx_acc        = p_sfc_flx%forc_slflx_acc       /REAL(nsteps_since_last_output,wp)
    p_sfc_flx%forc_precip_acc       = p_sfc_flx%forc_precip_acc      /REAL(nsteps_since_last_output,wp)
    p_sfc_flx%forc_evap_acc         = p_sfc_flx%forc_evap_acc        /REAL(nsteps_since_last_output,wp)
    p_sfc_flx%forc_runoff_acc       = p_sfc_flx%forc_runoff_acc      /REAL(nsteps_since_last_output,wp)
    p_sfc_flx%forc_fwbc_acc         = p_sfc_flx%forc_fwbc_acc        /REAL(nsteps_since_last_output,wp)
    p_sfc_flx%forc_fwrelax_acc      = p_sfc_flx%forc_fwrelax_acc     /REAL(nsteps_since_last_output,wp)
    p_sfc_flx%forc_fwfx_acc         = p_sfc_flx%forc_fwfx_acc        /REAL(nsteps_since_last_output,wp)
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
    p_sfc_flx%forc_wind_u_acc       = 0.0_wp
    p_sfc_flx%forc_wind_v_acc       = 0.0_wp
    p_sfc_flx%forc_swflx_acc        = 0.0_wp
    p_sfc_flx%forc_lwflx_acc        = 0.0_wp
    p_sfc_flx%forc_ssflx_acc        = 0.0_wp
    p_sfc_flx%forc_slflx_acc        = 0.0_wp
    p_sfc_flx%forc_precip_acc       = 0.0_wp
    p_sfc_flx%forc_evap_acc         = 0.0_wp
    p_sfc_flx%forc_runoff_acc       = 0.0_wp
    p_sfc_flx%forc_fwbc_acc         = 0.0_wp
    p_sfc_flx%forc_fwrelax_acc      = 0.0_wp
    p_sfc_flx%forc_fwfx_acc         = 0.0_wp
    p_sfc_flx%forc_hfrelax_acc      = 0.0_wp
    p_sfc_flx%forc_hflx_acc         = 0.0_wp
    p_sfc_flx%forc_tracer_acc       = 0.0_wp
    p_sfc_flx%forc_tracer_relax_acc = 0.0_wp

  END SUBROUTINE reset_ocean_statistics

  SUBROUTINE add_fields_3d(f_a,f_b,subset,levels)
    REAL(wp),INTENT(INOUT)          :: f_a(:,:,:)
    REAL(wp),INTENT(IN)             :: f_b(:,:,:)
    TYPE(t_subset_range),INTENT(IN) :: subset
    INTEGER,INTENT(IN),OPTIONAL     :: levels

    INTEGER :: mylevels
    INTEGER :: jb,jc,jk,jc_start,jc_end

    mylevels = n_zlev
    CALL assign_if_present(mylevels,levels)

    DO jb = subset%start_block, subset%end_block
      CALL get_index_range(subset, jb, jc_start, jc_end)
      DO jk=1,mylevels
        DO jc = jc_start, jc_end
          f_a(jc,jk,jb) = f_a(jc,jk,jb) + f_b(jc,jk,jb)
        END DO
      END DO
    END DO
  END SUBROUTINE add_fields_3d

  SUBROUTINE add_fields_2d(f_a,f_b,subset)
    REAL(wp),INTENT(INOUT)          :: f_a(:,:)
    REAL(wp),INTENT(IN)             :: f_b(:,:)
    TYPE(t_subset_range),INTENT(IN) :: subset

    INTEGER :: jb,jc,jc_start,jc_end

    DO jb = subset%start_block, subset%end_block
      CALL get_index_range(subset, jb, jc_start, jc_end)
      DO jc = jc_start, jc_end
        f_a(jc,jb) = f_a(jc,jb) + f_b(jc,jb)
      END DO
    END DO
  END SUBROUTINE add_fields_2d

  SUBROUTINE new_ocean_statistics()
  END SUBROUTINE new_ocean_statistics

  SUBROUTINE set_output_pointers

    TYPE(t_list_element), POINTER :: output_var, prog_var
    CHARACTER(len=max_char_length) :: timelevel
   !-------------------------------------------------------------------------
    WRITE(timelevel,'(a,i2.2)') '_TL',nnew(1)

    !CALL print_var_list(ocean_restart_list)
    output_var             =  find_list_element(ocean_restart_list,'h')
    prog_var               =  find_list_element(ocean_restart_list,'h'//TRIM(timelevel))
    IF (ASSOCIATED (output_var) .AND. ASSOCIATED(prog_var)) THEN
      output_var%field%r_ptr => prog_var%field%r_ptr
    ELSE
      IF (my_process_is_stdio()) THEN
        CALL finish('set_output_pointers', 'output_var h or h'//TRIM(timelevel)//'not ASSOCIATED')
      ENDIF
    ENDIF

    output_var             =  find_list_element(ocean_restart_list,'vn')
    prog_var               =  find_list_element(ocean_restart_list,'vn'//TRIM(timelevel))
    output_var%field%r_ptr => prog_var%field%r_ptr

    output_var             =  find_list_element(ocean_restart_list,'t')
    prog_var               =  find_list_element(ocean_restart_list,'t'//TRIM(timelevel))
    output_var%field%r_ptr => prog_var%field%r_ptr

    output_var             =  find_list_element(ocean_restart_list,'s')
    prog_var               =  find_list_element(ocean_restart_list,'s'//TRIM(timelevel))
    output_var%field%r_ptr => prog_var%field%r_ptr
  END SUBROUTINE set_output_pointers

END MODULE mo_hydro_ocean_run
