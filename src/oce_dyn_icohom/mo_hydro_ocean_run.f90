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
USE mo_model_domain,           ONLY: t_patch
USE mo_grid_config,            ONLY: n_dom
USE mo_sync,                   ONLY: sync_e, sync_c, sync_v, sync_patch_array
USE mo_ocean_nml,              ONLY: iswm_oce, n_zlev, no_tracer, &
  &                                  itestcase_oce, idiag_oce, init_oce_prog, init_oce_relax, &
  &                                  EOS_type, i_sea_ice, l_staggered_timestep
USE mo_dynamics_config,        ONLY: nold, nnew
USE mo_io_config,              ONLY: out_expname, istime4output, istime4newoutputfile,&
  &                                  is_checkpoint_time, n_checkpoints
USE mo_run_config,             ONLY: nsteps, dtime, ltimer
USE mo_exception,              ONLY: message, message_text, finish
USE mo_ext_data_types,         ONLY: t_external_data
USE mo_io_units,               ONLY: filename_max
USE mo_datetime,               ONLY: t_datetime, print_datetime, add_time, datetime_to_string
USE mo_timer,                  ONLY: timer_start, timer_stop, timer_total, timer_solve_ab,  &
  &                                  timer_tracer_ab, timer_vert_veloc, timer_normal_veloc, &
  &                                  timer_upd_phys, timer_upd_flx  !,timer_oce_init
USE mo_oce_ab_timestepping,    ONLY: solve_free_surface_eq_ab, &
  &                                  calc_normal_velocity_ab,  &
  &                                  calc_vert_velocity,       &
  &                                  update_time_indices
USE mo_oce_init,               ONLY: init_ho_testcases, init_ho_prog, init_ho_coupled,&
  &                                  init_ho_recon_fields, init_ho_relaxation, init_oce_index
USE mo_util_dbg_prnt,          ONLY: init_dbg_index
USE mo_oce_state,              ONLY: t_hydro_ocean_state, t_hydro_ocean_base, &
  &                                  init_ho_base, init_ho_basins, v_base, &
  &                                  construct_hydro_ocean_base, destruct_hydro_ocean_base, &
  &                                  construct_hydro_ocean_state, destruct_hydro_ocean_state, &
  &                                  init_coriolis_oce, init_oce_config, &
  &                                  set_lateral_boundary_values
USE mo_oce_math_operators,     ONLY: height_related_quantities
USE mo_operator_ocean_coeff_3d,ONLY: t_operator_coeff, allocate_exp_coeff, par_init_operator_coeff
USE mo_scalar_product,         ONLY: calc_scalar_product_veloc_3D
USE mo_oce_tracer,             ONLY: advect_tracer_ab
USE mo_io_restart,             ONLY: write_restart_info_file
USE mo_intp_data_strc,         ONLY: t_int_state
USE mo_oce_bulk,               ONLY: update_sfcflx
USE mo_sea_ice,                ONLY: construct_sfcflx,destruct_sfcflx,&
  &                                  construct_atmos_for_ocean,&
  &                                  destruct_atmos_for_ocean,&
  &                                  construct_atmos_fluxes, destruct_atmos_fluxes,&
  &                                  construct_sea_ice, destruct_sea_ice, &
  &                                  ice_init, ice_slow
! #
USE mo_sea_ice_types,          ONLY: t_sfc_flx, t_atmos_fluxes, t_atmos_for_ocean, &
  &                                  t_sea_ice
USE mo_oce_forcing,            ONLY: init_sfcflx
USE mo_oce_physics,            ONLY: t_ho_params, &
  &                                  construct_ho_params, init_ho_params, &
  &                                  destruct_ho_params, update_ho_params
USE mo_oce_thermodyn,          ONLY: calc_density_MPIOM_func, calc_density_lin_EOS_func,&
  &                                  calc_density_JMDWFG06_EOS_func, calc_density
USE mo_output,                 ONLY: init_output_files, write_output, &
  &                                  create_restart_file
USE mo_oce_diagnostics,        ONLY: calculate_oce_diagnostics,&
  &                                  construct_oce_diagnostics,&
  &                                  destruct_oce_diagnostics, t_oce_timeseries, &
  &                                  calc_moc, calc_psi
USE mo_mpi,                    ONLY: my_process_is_mpi_all_parallel
IMPLICIT NONE

PRIVATE
INTEGER, PARAMETER :: kice = 1

INCLUDE 'cdi.inc'

!VERSION CONTROL:
CHARACTER(LEN=*), PARAMETER :: version = '$Id$'

! public interface
!
! public subroutines
PUBLIC :: perform_ho_stepping
PUBLIC :: prepare_ho_integration
PUBLIC :: finalise_ho_integration
!
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
  SUBROUTINE perform_ho_stepping( ppatch, pstate_oce, p_ext_data,               &
                                & datetime, n_io, jfile, lwrite_restart, p_int, &
                                & p_sfc_flx, p_phys_param,                      &
                                & p_as, p_atm_f, p_ice,ptr_op_coeff,            &
                                & l_have_output)

  TYPE(t_patch),             TARGET, INTENT(INOUT) :: ppatch(n_dom)
  TYPE(t_hydro_ocean_state), TARGET, INTENT(INOUT) :: pstate_oce(n_dom)
  TYPE(t_external_data), TARGET, INTENT(IN)        :: p_ext_data(n_dom)
  TYPE(t_datetime), INTENT(INOUT)                  :: datetime
  INTEGER, INTENT(IN)                              :: n_io
  INTEGER, INTENT(INOUT)                           :: jfile
  LOGICAL, INTENT(IN)                              :: lwrite_restart
  TYPE(t_int_state),TARGET,INTENT(inout)           :: p_int(n_dom)
  TYPE(t_sfc_flx)                                  :: p_sfc_flx
  TYPE (t_ho_params)                               :: p_phys_param
  TYPE(t_atmos_for_ocean),  INTENT(INOUT)          :: p_as
  TYPE(t_atmos_fluxes ),    INTENT(INOUT)          :: p_atm_f
  TYPE (t_sea_ice),         INTENT(INOUT)          :: p_ice
  TYPE(t_operator_coeff),   INTENT(INOUT)          :: ptr_op_coeff
  LOGICAL,                  INTENT(INOUT)          :: l_have_output



  ! local variables
  INTEGER                         :: jstep, jg
  LOGICAL                         :: l_outputtime
  CHARACTER(len=32)               :: datestring
  CHARACTER(len=36)               :: moc_fname
  TYPE(t_oce_timeseries), POINTER :: oce_ts
  !TYPE(t_operator_coeff)          :: ptr_op_coeff

  !CHARACTER(LEN=filename_max)  :: outputfile, gridfile
  CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
    &      routine = 'mo_hydro_ocean_run:perform_ho_stepping'
  !------------------------------------------------------------------

  !------------------------------------------------------------------
  ! no grid refinement allowed here so far
  !------------------------------------------------------------------

  IF (n_dom > 1 ) THEN
    CALL finish(TRIM(routine), ' N_DOM > 1 is not allowed')
  END IF
  jg = n_dom

!   CALL allocate_exp_coeff( ppatch(jg), ptr_op_coeff)
!   CALL par_init_operator_coeff( ppatch(jg), ptr_op_coeff, p_int(jg))
  ! As an alternative a serial 2D version of the coefficient computation can be used
!   IF (my_process_is_mpi_all_parallel()) THEN
!     CALL par_init_operator_coeff( ppatch(jg), ptr_op_coeff, p_int(jg))
!   ELSE
!     CALL init_operator_coeff( ppatch(jg), ptr_op_coeff)
!   ENDIF

!   CALL init_ho_recon_fields( ppatch(jg), pstate_oce(jg), ptr_op_coeff)

  IF (idiag_oce == 1) &
    & CALL construct_oce_diagnostics( ppatch(jg), pstate_oce(jg), p_ext_data(jg), oce_ts)

  IF (ltimer) CALL timer_start(timer_total)

  ! open file for MOC - extraordinary at this time
  CALL datetime_to_string(datestring, datetime)
  moc_fname='MOC.'//TRIM(datestring)
  !IF (my_process_is_stdio()) THEN
  OPEN (77,file=moc_fname,form='unformatted')
  WRITE(message_text,'(2a)') ' MOC-file opened successfully, filename=',TRIM(moc_fname)
  CALL message (TRIM(routine), message_text)
  !END IF

  ! call of MOC before time loop
  !CALL calc_moc (ppatch(jg), pstate_oce(jg)%p_diag%w(:,:,:), datetime)


  !------------------------------------------------------------------
  ! call the dynamical core: start the time loop
  !------------------------------------------------------------------
  TIME_LOOP: DO jstep = 1, nsteps

    CALL datetime_to_string(datestring, datetime)
    WRITE(message_text,'(a,i6,2a)') '  Begin of timestep =',jstep,'  datetime:  ', datestring
    CALL message (TRIM(routine), message_text)

    IF(itestcase_oce==28)THEN
      CALL height_related_quantities(ppatch(jg), pstate_oce(jg), p_ext_data(jg))
      CALL calc_vert_velocity( ppatch(jg), pstate_oce(jg),ptr_op_coeff)
      CALL advect_tracer_ab(ppatch(jg), pstate_oce(jg), &
                           &p_phys_param,p_sfc_flx,&
                           & ptr_op_coeff,&
                           & jstep)
    ELSE


      !In case of a time-varying forcing:
      IF (ltimer) CALL timer_start(timer_upd_flx)
      CALL update_sfcflx(ppatch(jg), pstate_oce(jg), p_as, p_ice, p_atm_f, p_sfc_flx, &
        &                jstep, datetime)

      IF(.NOT.l_STAGGERED_TIMESTEP)THEN

        CALL height_related_quantities(ppatch(jg), pstate_oce(jg), p_ext_data(jg))

        !This is required in top boundary condition for
        !vertical velocity: the time derivative of the surface height
        !is used there and needs special treatment in the first timestep.
        !see sbr top_bound_cond_vert_veloc in mo_ho_boundcond
        pstate_oce(jg)%p_prog(nnew(1))%h = pstate_oce(jg)%p_prog(nold(1))%h

        Call set_lateral_boundary_values(ppatch(jg), pstate_oce(jg)%p_prog(nold(1))%vn)
        CALL sync_patch_array(sync_e, ppatch(jg),  pstate_oce(jg)%p_prog(nold(1))%vn)

!          CALL calc_scalar_product_veloc( ppatch(jg), &
!            & pstate_oce(jg)%p_prog(nold(1))%vn,&
!            & pstate_oce(jg)%p_prog(nold(1))%vn,&
!            & pstate_oce(jg)%p_diag%h_e,        &
!            & pstate_oce(jg)%p_diag)

         CALL calc_scalar_product_veloc_3D( ppatch(jg), &
           & pstate_oce(jg)%p_prog(nold(1))%vn,         &
           & pstate_oce(jg)%p_prog(nold(1))%vn,         &
           & pstate_oce(jg)%p_diag%h_e,                 &
           & pstate_oce(jg)%p_diag,                     &
           & ptr_op_coeff)

      ENDIF
      IF (ltimer) CALL timer_stop(timer_upd_flx)

      IF (ltimer) CALL timer_start(timer_upd_phys)
      ! Calc density globally
      !CALL calc_density(ppatch(jg),&
      !  &               pstate_oce(jg)%p_prog(nold(1))%tracer,&
      !  &               pstate_oce(jg)%p_diag%rho)

      SELECT CASE (EOS_TYPE)
      CASE(1)
        CALL update_ho_params(ppatch(jg), pstate_oce(jg), p_sfc_flx, p_phys_param,&
          &                   calc_density_lin_EOS_func)

      CASE(2)
        CALL update_ho_params(ppatch(jg), pstate_oce(jg), p_sfc_flx, p_phys_param,&
          &                   calc_density_MPIOM_func)
      CASE(3)
        CALL update_ho_params(ppatch(jg), pstate_oce(jg), p_sfc_flx, p_phys_param,&
          &                   calc_density_JMDWFG06_EOS_func)
      CASE DEFAULT
      END SELECT
      IF (ltimer) CALL timer_stop(timer_upd_phys)

      ! solve for new free surface
      IF (ltimer) CALL timer_start(timer_solve_ab)
      CALL solve_free_surface_eq_ab (ppatch(jg), pstate_oce(jg), p_ext_data(jg), &
        &                            p_sfc_flx, p_phys_param, jstep, ptr_op_coeff, p_int(jg))
      IF (ltimer) CALL timer_stop(timer_solve_ab)

      ! Step 4: calculate final normal velocity from predicted horizontal velocity vn_pred
      !         and updated surface height
      IF (ltimer) CALL timer_start(timer_normal_veloc)
      CALL calc_normal_velocity_ab(ppatch(jg), pstate_oce(jg),&
                                  &ptr_op_coeff, p_ext_data(jg), p_phys_param)
      IF (ltimer) CALL timer_stop(timer_normal_veloc)

      ! Step 5: calculate vertical velocity from continuity equation under incompressiblity condition
      ! in the non-shallow-water case
      IF ( iswm_oce /= 1 ) THEN
        IF (ltimer) CALL timer_start(timer_vert_veloc)
        CALL calc_vert_velocity( ppatch(jg), pstate_oce(jg),ptr_op_coeff)
        IF (ltimer) CALL timer_stop(timer_vert_veloc)
      ENDIF

      ! Step 6: transport tracers and diffuse them
      IF (no_tracer>=1) THEN
        IF (ltimer) CALL timer_start(timer_tracer_ab)
        CALL advect_tracer_ab(ppatch(jg), pstate_oce(jg), p_phys_param,&
                             &p_sfc_flx,&
                             &ptr_op_coeff,&
                             & jstep)
        IF (ltimer) CALL timer_stop(timer_tracer_ab)
      ENDIF

    ENDIF  ! testcase 28

   ! Actually diagnostics for 3D not implemented, PK March 2011
    IF (idiag_oce == 1 ) THEN
      CALL calculate_oce_diagnostics( ppatch(jg),    &
                                    & pstate_oce(jg),&
                                    & p_sfc_flx,     &
                                    & p_phys_param,  &
                                    & jstep,         &
                                    & oce_ts)
    ENDIF

    ! One integration cycle finished on the lowest grid level (coarsest
    ! resolution). Set model time.
    CALL add_time(dtime,0,0,0,datetime)

    l_outputtime = (MOD(jstep,n_io) == 0)
    IF ( l_outputtime ) THEN

      CALL calc_moc (ppatch(jg), pstate_oce(jg)%p_diag%w(:,:,:), datetime)
      CALL calc_psi (ppatch(jg), pstate_oce(jg)%p_diag%u(:,:,:), &
        &                        pstate_oce(jg)%p_prog(nold(1))%h(:,:), &
        &                        pstate_oce(jg)%p_diag%u_vint, datetime)

      CALL write_output( datetime )
      CALL message (TRIM(routine),'Write output at:')
      CALL print_datetime(datetime)
      l_have_output = .TRUE.

    END IF

    ! close the current output file and trigger a new one
    IF (istime4newoutputfile(jstep)) THEN
      jfile = jfile +1
      CALL init_output_files(jfile,lclose=l_have_output)
    END IF

    ! Shift time indices for the next loop
    ! this HAS to ge into the restart files, because the start with the following loop
    CALL update_time_indices(jg)
    ! update intermediate timestepping variables for the tracers
    CALL update_intermediate_tracer_vars(pstate_oce(jg))

    ! write a restart or checkpoint file
    IF (MOD(jstep,n_checkpoints())==0 .OR. (jstep==nsteps .AND. lwrite_restart)) THEN
      CALL create_restart_file( ppatch(jg), datetime,  &
                              & jfile, l_have_output,opt_depth=n_zlev)
      ! Create the master (meta) file in ASCII format which contains
      ! info about which files should be read in for a restart run.
      CALL write_restart_info_file
    END IF


  ENDDO TIME_LOOP

  IF (idiag_oce==1) CALL destruct_oce_diagnostics(oce_ts)

  IF (ltimer) CALL timer_stop(timer_total)

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
  SUBROUTINE prepare_ho_integration(ppatch, pstate_oce, p_ext_data, p_sfc_flx, &
                                  & p_phys_param, p_as,&
                                  & p_atm_f, p_ice, ptr_op_coeff, p_int)

    TYPE(t_patch),                INTENT(INOUT)  :: ppatch(n_dom)
    TYPE(t_hydro_ocean_state),    INTENT(INOUT)  :: pstate_oce(n_dom)
    TYPE(t_external_data),        INTENT(INOUT)  :: p_ext_data(n_dom)
    TYPE(t_sfc_flx),              INTENT(INOUT)  :: p_sfc_flx
    TYPE (t_ho_params),           INTENT(INOUT)  :: p_phys_param
    TYPE(t_atmos_for_ocean ),     INTENT(INOUT)  :: p_as
    TYPE(t_atmos_fluxes ),        INTENT(INOUT)  :: p_atm_f
    TYPE (t_sea_ice),             INTENT(INOUT)  :: p_ice
    TYPE(t_operator_coeff),       INTENT(INOUT)  :: ptr_op_coeff
    TYPE(t_int_state),TARGET,     INTENT(inout)  :: p_int(n_dom)

    ! local variables
    !TYPE(t_hydro_ocean_base)                  :: p_base
    INTEGER :: jg
    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
      &      routine = 'mo_test_hydro_ocean:prepare_ho_integration'

    CALL message (TRIM(routine),'start')
    !------------------------------------------------------------------
    ! no grid refinement allowed here so far
    !------------------------------------------------------------------

 !  IF (ltimer) CALL timer_start(timer_oce_init)

    IF (n_dom > 1 ) THEN
      CALL finish(TRIM(routine), ' N_DOM > 1 is not allowed')
    END IF
    jg = n_dom

    !------------------------------------------------------------------
    ! construct ocean state and physics
    !------------------------------------------------------------------
    CALL init_oce_config

    ! initialize ocean indices for debug output (before ocean state, no 3-dim)
    CALL init_dbg_index(ppatch(1))

    ! hydro_ocean_base contains the 3-dimensional structures for the ocean state
    CALL construct_hydro_ocean_base(ppatch(jg), v_base)
    CALL init_ho_base     (ppatch(jg), p_ext_data(jg), v_base)
    CALL init_ho_basins   (ppatch(jg),                 v_base)
    CALL init_coriolis_oce( ppatch(jg) )

    !------------------------------------------------------------------
    ! construct ocean state and physics
    !------------------------------------------------------------------

    ! ppatch and pstate_oce have dimension n_dom
    CALL construct_hydro_ocean_state(ppatch, pstate_oce)

    ! initialize ocean indices for debug output (including 3-dim lsm)
    CALL init_oce_index ( ppatch, pstate_oce, p_ext_data )

    CALL construct_ho_params(ppatch(jg), p_phys_param)
    CALL init_ho_params(ppatch(jg), p_phys_param)

    !------------------------------------------------------------------
    ! construct ocean forcing and testcases
    !------------------------------------------------------------------

    CALL construct_sfcflx(ppatch(jg), p_sfc_flx)
    CALL      init_sfcflx(ppatch(jg), p_sfc_flx)

    CALL construct_sea_ice(ppatch(jg), p_ice, kice)
    CALL construct_atmos_for_ocean(ppatch(jg), p_as)
    CALL construct_atmos_fluxes(ppatch(jg), p_atm_f, kice)

    IF (init_oce_prog == 0) THEN
      CALL init_ho_testcases(ppatch(jg), pstate_oce(jg), p_ext_data(jg), ptr_op_coeff,p_sfc_flx)
    ELSE IF (init_oce_prog == 1) THEN
      CALL init_ho_prog(ppatch(jg), pstate_oce(jg), p_sfc_flx)
    END IF

    IF (init_oce_relax == 1) THEN
      CALL init_ho_relaxation(ppatch(jg), pstate_oce(jg), p_sfc_flx)
    END IF

    CALL init_ho_coupled(ppatch(jg), pstate_oce(jg))
    IF (i_sea_ice >= 1) &
      &   CALL ice_init(ppatch(jg), pstate_oce(jg), p_ice)


    CALL allocate_exp_coeff( ppatch(jg), ptr_op_coeff)
    CALL par_init_operator_coeff( ppatch(jg), ptr_op_coeff, p_int(jg))
    CALL init_ho_recon_fields( ppatch(jg), pstate_oce(jg), ptr_op_coeff)

  ! IF (ltimer) CALL timer_stop(timer_oce_init)
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
     CALL destruct_ho_params(p_phys_param)

     CALL destruct_sfcflx(p_sfc_flx)
     CALL destruct_sea_ice(p_ice)
     CALL destruct_atmos_for_ocean(p_as)
     CALL destruct_atmos_fluxes(p_atm_f)


  END SUBROUTINE finalise_ho_integration

  SUBROUTINE update_intermediate_tracer_vars(p_os)
    TYPE(t_hydro_ocean_state), INTENT(INOUT) :: p_os

    INTEGER :: it

    ! tracer updates
    DO it = 1,no_tracer
      !horiz
      p_os%p_aux%g_nm1_c_h(:,:,:,it)  = p_os%p_aux%g_n_c_h(:,:,:,it)
      p_os%p_aux%g_n_c_h(:,:,:,it)    = 0.0_wp
      p_os%p_aux%g_nimd_c_h(:,:,:,it) = 0.0_wp

      !vert
      p_os%p_aux%g_nm1_c_v(:,:,:,it)  = p_os%p_aux%g_n_c_v(:,:,:,it)
      p_os%p_aux%g_n_c_v(:,:,:,it)    = 0.0_wp
      p_os%p_aux%g_nimd_c_v(:,:,:,it) = 0.0_wp
    END DO
    ! vertical velocity
    p_os%p_aux%g_nm1 = p_os%p_aux%g_n
    p_os%p_aux%g_n   = 0.0_wp
  END SUBROUTINE update_intermediate_tracer_vars

END MODULE mo_hydro_ocean_run

