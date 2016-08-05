!>
!! Contains subroutines to prepare and organize time steps in the hydrostatic.
!!
!! Contains subroutines to prepare and organize time steps in the hydrostatic
!! model. Formerly, this was contained in the main program hydro_atmos, but as
!! the main program should only be a driver for different kinds of equations,
!! the present module is provided.
!!
!! @par Revision History
!! initial release by Almut Gassmann, MPI-M, (2009-03-04)
!!
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
#include "omp_definitions.inc"
!----------------------------
MODULE mo_ha_stepping

  USE mo_kind,                ONLY: wp
  USE mo_icoham_dyn_types,    ONLY: t_hydro_atm
  USE mo_icoham_dyn_memory,   ONLY: construct_icoham_dyn_state
  USE mo_intp_data_strc,      ONLY: t_int_state
  USE mo_datetime,            ONLY: t_datetime, print_datetime, add_time
  USE mo_exception,           ONLY: message, message_text
  USE mo_model_domain,        ONLY: t_patch
  USE mo_ext_data_state,      ONLY: ext_data
  USE mo_grid_config,         ONLY: n_dom
  USE mo_dynamics_config,     ONLY: lshallow_water, ltwotime, nnow, nold
  USE mo_ha_dyn_config,       ONLY: ha_dyn_config, configure_ha_dyn
  USE mo_io_config,           ONLY: l_diagtime, l_outputtime, is_checkpoint_time, &
    &                               n_chkpt, n_diag
  USE mo_run_config,          ONLY: nsteps, dtime, ntracer,  &
                                  & ldynamics, ltransport, msg_level,   &
                                  & ltestcase, output_mode
  USE mo_master_config,       ONLY: isRestart
  USE mo_ha_testcases,        ONLY: init_testcase
  USE mo_si_correction,       ONLY: init_si_params
  USE mo_ha_rungekutta,       ONLY: init_RungeKutta
  USE mo_ha_diagnostics,      ONLY: supervise_total_integrals, &
                                  & init_total_integrals
  USE mo_ha_prog_util,        ONLY: copy_prog_state
  USE mo_ha_diag_util,        ONLY: update_diag_state, update_dyn_output
  USE mo_expensive_functions, ONLY: convert_t2theta
  USE mo_hierarchy_management,ONLY: process_grid, interpolate_diagnostics
  USE mo_grf_intp_data_strc,  ONLY: t_gridref_state
  USE mo_impl_constants,      ONLY: LEAPFROG_EXPL, LEAPFROG_SI, &
                                    RK4, SSPRK54, MAX_CHAR_LENGTH
  USE mo_timer,               ONLY: ltimer, timer_start, timer_stop, timer_total, timer_intrp_diagn
  USE mo_sync,                ONLY: global_max
  USE mo_vertical_coord_table,ONLY: vct
  
  USE mo_icon_comm_lib,       ONLY: icon_comm_sync_all
  USE mo_parallel_config,     ONLY: use_icon_comm, use_async_restart_output
  USE mo_name_list_output,    ONLY: write_name_list_output, istime4name_list_output
  USE mo_restart,             ONLY: t_RestartDescriptor, createRestartDescriptor, deleteRestartDescriptor
  USE mo_restart_attributes,  ONLY: t_RestartAttributeList, getAttributesForRestarting
  USE mo_time_config,         ONLY: time_config

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: prepare_ha_dyn, initcond_ha_dyn
  PUBLIC :: perform_ha_stepping

CONTAINS

  !>
  !!
  !! Hui Wan (MPI-M, 2011-05)
  !!
  SUBROUTINE prepare_ha_dyn( p_patch )

    TYPE(t_patch),TARGET,INTENT(IN) :: p_patch(n_dom)
    INTEGER :: nadd, ntl

    !-----------------------------------
    ! Set up scheme-specific constants
    !-----------------------------------
    CALL configure_ha_dyn

    SELECT CASE (ha_dyn_config%itime_scheme)
    CASE (LEAPFROG_SI)
      CALL init_si_params( ha_dyn_config%lsi_3d,            &
                           ha_dyn_config%si_offctr,         &
                           ha_dyn_config%si_cmin,           &
                           lshallow_water                   )
    CASE (RK4,SSPRK54)
      CALL init_RungeKutta(ha_dyn_config%itime_scheme)
    END SELECT

    !---------------------------------------------------------------
    ! Open file and allocate memory for diagnosing global integrals
    !---------------------------------------------------------------
    CALL init_total_integrals

    !-------------------------------------------------------------
    ! Allocate memory for the state vector of the dynamical core
    !-------------------------------------------------------------
    ! In case of grid refinement, an extra "time level" is added
    ! for computing boundary tendencies and feedback increments

    IF (n_dom > 1) THEN ; nadd = 1
    ELSE                ; nadd = 0
    ENDIF

    SELECT CASE (ha_dyn_config%itime_scheme)
    CASE (LEAPFROG_EXPL, LEAPFROG_SI, RK4, SSPRK54)
    ! 3 time level scheme, or 2 time level scheme with multi stages
      ntl = 3 + nadd
    CASE DEFAULT
      ntl = 2 + nadd
    END SELECT

    CALL construct_icoham_dyn_state( ntl, ntracer, p_patch )

  END SUBROUTINE prepare_ha_dyn
  !-------------
  !>
  !!
  SUBROUTINE initcond_ha_dyn( p_patch, p_int_state, p_grf_state, p_hydro_state )

    TYPE(t_patch),        TARGET,INTENT(INout) :: p_patch(n_dom)
    TYPE(t_int_state),    TARGET,INTENT(IN)    :: p_int_state(n_dom)
    TYPE(t_gridref_state),TARGET,INTENT(IN)    :: p_grf_state(n_dom)
    TYPE(t_hydro_atm),    TARGET,INTENT(INOUT) :: p_hydro_state(n_dom)

    INTEGER :: jg, jn, jgc

    !--------------------------------------------------
    ! Assign initial condition to prognostic variables
    !--------------------------------------------------
    IF (ltestcase) CALL init_testcase( p_patch, p_hydro_state, p_int_state, ext_data)

    ! Convert temperature to potential temperature if the latter is
    ! chosen as a prognostic variable.

    IF (ha_dyn_config%ltheta_dyn) THEN
      DO jg = 1, n_dom

        CALL convert_t2theta(p_patch(jg),                      &
                             p_hydro_state(jg)%prog(nnow(jg)), &
                             p_hydro_state(jg)%diag           )

        IF (.NOT.ltwotime)                 &
        CALL convert_t2theta(p_patch(jg),                      &
                             p_hydro_state(jg)%prog(nold(jg)), &
                             p_hydro_state(jg)%diag           )
      ENDDO
    ENDIF

   !--------------------------------------------------------------
   ! Copy initial condition to output buffer; update diagnostic
   ! variables for output
   !--------------------------------------------------------------
   DO jg=1,n_dom
     CALL copy_prog_state( p_hydro_state(jg)%prog(nnow(jg)), &! in
                           p_hydro_state(jg)%prog_out,       &! out
                           ha_dyn_config%ltheta_dyn,         &! in
                           ltransport                        )! in

     CALL update_diag_state( p_hydro_state(jg)%prog_out,    &! in
                             p_patch(jg), p_int_state(jg),  &! in
                             ext_data(jg),                  &! in
                             p_hydro_state(jg)%diag_out )    ! out

     CALL update_dyn_output( p_patch(jg), p_int_state(jg),  &! in
                             p_hydro_state(jg)%prog_out,    &! in
                             p_hydro_state(jg)%diag_out )    ! inout

     IF (use_icon_comm) CALL icon_comm_sync_all()

     ! Fill boundary cells of nested domains
     DO jn = 1, p_patch(jg)%n_childdom
       jgc = p_patch(jg)%child_id(jn)
       CALL interpolate_diagnostics(p_patch,                      &
                                    p_hydro_state(jg )%diag_out,  &
                                    p_hydro_state(jgc)%diag_out,  &
                                    p_int_state,p_grf_state,jg,jgc)
     ENDDO
   ENDDO

  END SUBROUTINE initcond_ha_dyn
  !--------
  !>
  !! Organizes the hydrostatic model time stepping.
  !!
  !!
  !! @par Revision History
  !! Initial release by Almut Gassmann (2009-03-04)
  !!
  SUBROUTINE perform_ha_stepping( p_patch, p_int_state,               &
                                & p_grf_state,                        &
                                & p_hydro_state, datetime             )

  CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = 'mo_ha_stepping:perform_ha_stepping'

  TYPE(t_patch),         TARGET, INTENT(IN)    :: p_patch(n_dom)
  TYPE(t_int_state),     TARGET, INTENT(IN)    :: p_int_state(n_dom)
  TYPE(t_gridref_state), TARGET, INTENT(INOUT) :: p_grf_state(n_dom)

  TYPE(t_datetime), INTENT(INOUT)              :: datetime
  TYPE(t_hydro_atm), TARGET, INTENT(INOUT)     :: p_hydro_state(n_dom)

  REAL(wp), DIMENSION(:,:,:), POINTER          :: p_vn  => NULL()
  REAL(wp)                                     :: sim_time(n_dom)
  INTEGER                                      :: jg, jn, jgc, jstep
  LOGICAL                                      :: l_nml_output
  LOGICAL                                      :: l_3tl_init(n_dom)
  INTEGER                                      :: jstep0 ! start counter for time loop
  TYPE(t_RestartAttributeList), POINTER        :: restartAttributes
  CLASS(t_RestartDescriptor), POINTER          :: restartDescriptor

#ifdef _OPENMP
  INTEGER  :: jb
  REAL(wp) :: vnmax, vn_aux(p_patch(1)%nblks_e)
#else
  REAL(wp) :: vnmax(1)
#endif


!-----------------------------------------------------------------------
  sim_time(:) = 0._wp

  IF (ltimer) CALL timer_start(timer_total)

  restartDescriptor => createRestartDescriptor()

  jstep0 = 0
  restartAttributes => getAttributesForRestarting()
  IF (ASSOCIATED(restartAttributes) .AND. .NOT. time_config%is_relative_time) THEN
    ! get start counter for time loop from restart file:
    jstep0 = restartAttributes%getInteger("jstep")
  END IF

  TIME_LOOP: DO jstep = (jstep0+1), (jstep0+nsteps)

    !--------------------------------------------------------------------------
    ! Send to stdout information about the current integration cycle
    ! (maximum normal wind).
    !--------------------------------------------------------------------------

    IF (msg_level >= 5) THEN ! compute maximum velocity in global model domain
      p_vn => p_hydro_state(1)%prog(nnow(1))%vn
      vnmax = -HUGE(vnmax)
#if (!defined __xlC__) && (defined _OPENMP)
!$OMP PARALLEL DO PRIVATE(jb) ICON_OMP_DEFAULT_SCHEDULE
! GZ: The reduction operation that used to be here tended to be slow, this one seems to be faster
      DO jb = 1, p_patch(1)%nblks_e
        vn_aux(jb) = MAXVAL(ABS(p_vn(:,:,jb)))
      ENDDO
!$OMP END PARALLEL DO
      vnmax = MAXVAL(vn_aux)
#else
      vnmax = MAXVAL(ABS(p_vn))
#endif
      vnmax = global_max(vnmax) ! Get max over all PEs
      WRITE(message_text,'(a,e14.6)') 'MAXABS VN ', vnmax
      CALL message(TRIM(routine),message_text)
    ENDIF

    !--------------------------------------------------------------------------
    ! Time integration from time step n to n+1
    !--------------------------------------------------------------------------
    ! Whether to apply a special initial treatment for
    ! 3-time-level schemes. In case of a restart run, the treatment is
    ! not necessary, thus the variable is set to .TRUE.

    l_3tl_init(1:n_dom) = (.NOT.ltwotime).AND.(.NOT.isRestart()).AND.((jstep+jstep0)==1)

    ! Call recursive subroutine 'process_grid_level', which executes
    ! one timestep for the global domain and calls itself in the presence
    ! of nested domains with recursively halved time steps

    CALL process_grid( p_patch, p_hydro_state, p_int_state, p_grf_state,    &
      &                ext_data, 1, jstep, l_3tl_init, dtime, sim_time,     &
      &                1, datetime )

    !--------------------------------------------------------------------------
    ! One integration cycle finished on the lowest grid level (coarsest
    ! resolution). Set model time.
    !--------------------------------------------------------------------------
    CALL add_time(dtime,0,0,0,datetime)
!!$    ! Not nice, but the name list output requires this
!!$    sim_time(1) = MODULO(sim_time(1) + dtime, 86400.0_wp)
    sim_time(1) = sim_time(1) + dtime   ! RS: is this correct? process_grid already advances sim_time by dtime !

    !--------------------------------------------------------------------------
    ! Set output flags
    !--------------------------------------------------------------------------
    ! The integration cycle "jstep" computes the atmospheric state variables
    ! at the new time step (n+1) using the already known step n (and possibly
    ! n-1, n-2, ... as well, depending on the user's choice of time stepping
    ! scheme). In some schemes (e.g. leapfrog), the step n values will be
    ! modified within this integration cycle. Therefore at the end of the
    ! cycle, we write out variables of time step n rather than n+1.

    l_nml_output   = output_mode%l_nml     .AND. &
      & ((jstep==(nsteps+jstep0)) .OR. istime4name_list_output(jstep))
    l_outputtime = l_nml_output

    IF ( MOD(jstep-1,n_diag)==0 .OR. jstep==(nsteps+jstep0) ) THEN
      l_diagtime = .TRUE. ! Diagnostic output is done at the end of the
    ELSE                  ! time step
      l_diagtime = .FALSE.
    ENDIF
    IF(.NOT.ldynamics.AND.lshallow_water)l_diagtime=.FALSE.

    !--------------------------------------------------------------------------
    ! Write output (prognostic and diagnostic variables)
    !--------------------------------------------------------------------------
    IF (l_outputtime) THEN

      !====================
      ! Prepare for output
      !====================
      DO jg = 1, n_dom
        CALL copy_prog_state( p_hydro_state(jg)%prog(nnow(jg)),  &! in
          &                   p_hydro_state(jg)%prog_out,     &! out
          &                   ha_dyn_config%ltheta_dyn,       &! in
          &                   ltransport                     ) ! in

        CALL update_diag_state( p_hydro_state(jg)%prog_out,   &! in
          &                     p_patch(jg), p_int_state(jg), &! in
          &                     ext_data(jg),                 &! in
          &                     p_hydro_state(jg)%diag_out )   ! out

        CALL update_dyn_output( p_patch(jg), p_int_state(jg), &! in
          &                     p_hydro_state(jg)%prog_out,   &! in
          &                     p_hydro_state(jg)%diag_out )   ! inout
      ENDDO

      IF (ltimer) CALL timer_start(timer_intrp_diagn)

      ! Interpolate diagnostic variables to nest boundaries
      IF (n_dom > 1) THEN
        DO jg = 1, n_dom-1
          DO jn = 1, p_patch(jg)%n_childdom
            jgc = p_patch(jg)%child_id(jn)
            CALL interpolate_diagnostics(p_patch,                      &
                                         p_hydro_state(jg )%diag_out,  &
                                         p_hydro_state(jgc)%diag_out,  &
                                         p_int_state,p_grf_state,jg,jgc)
          ENDDO
        ENDDO
      ENDIF

      IF (l_nml_output) THEN
        CALL message(TRIM(routine),'Output (name_list) at:')
        CALL print_datetime(datetime)
        CALL write_name_list_output(jstep)
      ENDIF

      IF (ltimer) CALL timer_stop(timer_intrp_diagn)

    ENDIF !l_outputtime

    !--------------------------------------------------------------------------
    ! Diagnose global integrals
    !--------------------------------------------------------------------------

    IF (l_diagtime) &
    CALL supervise_total_integrals( jstep-jstep0, p_patch, p_hydro_state, &
                                    ext_data(1:n_dom), nnow(1:n_dom), jstep == (nsteps+jstep0) )
    !--------------------------------------------------------------------------
    ! Write restart file
    !--------------------------------------------------------------------------
    IF (is_checkpoint_time(jstep,n_chkpt,(nsteps+jstep0))) THEN
        DO jg = 1, n_dom
            CALL restartDescriptor%updatePatch(p_patch(jg), opt_pvct = vct)
        ENDDO
        CALL restartDescriptor%writeRestart(datetime, jstep, "atm")
    END IF

  ENDDO TIME_LOOP

  CALL deleteRestartDescriptor(restartDescriptor)

  IF (ltimer) CALL timer_stop(timer_total)

  END SUBROUTINE perform_ha_stepping

END MODULE mo_ha_stepping
