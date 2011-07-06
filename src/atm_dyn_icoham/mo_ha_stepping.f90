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
!! @par Copyright
!! 2002-2009 by DWD and MPI-M
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
MODULE mo_ha_stepping

  USE mo_kind,                ONLY: wp
  USE mo_icoham_dyn_types,    ONLY: t_hydro_atm
  USE mo_icoham_dyn_memory,   ONLY: construct_icoham_dyn_state
  USE mo_interpolation,       ONLY: t_int_state
  USE mo_datetime,            ONLY: t_datetime, print_datetime, add_time
  USE mo_exception,           ONLY: message, message_text
  USE mo_model_domain,        ONLY: t_patch
  USE mo_ext_data,            ONLY: ext_data
  USE mo_model_domain_import, ONLY: n_dom
  USE mo_atm_dyn_nml,         ONLY: ltwotime, itime_scheme, nold, nnow
  USE mo_io_nml,              ONLY: l_outputtime, lprepare_output, l_diagtime,  &
                                  & l_checkpoint_time
  USE mo_run_nml,             ONLY: lshallow_water, nsteps, dtime,ltheta_dyn,   &
                                  & ldynamics, ltransport, msg_level, ltimer,   &
                                  & ltestcase, lrestart
  USE mo_hydro_testcases,     ONLY: init_testcase
  USE mo_si_correction,       ONLY: init_si_params
  USE mo_ha_rungekutta,       ONLY: init_RungeKutta
  USE mo_ha_diagnostics,      ONLY: supervise_total_integrals, &
                                  & init_total_integrals
  USE mo_ha_prog_util,        ONLY: copy_prog_state
  USE mo_ha_diag_util,        ONLY: update_diag_state, update_dyn_output
  USE mo_expensive_functions, ONLY: convert_t2theta
  USE mo_output,              ONLY: init_output_files, write_output, &
                                  & create_restart_file
  USE mo_hierarchy_management,ONLY: process_grid, interpolate_diagnostics
  USE mo_grf_interpolation,   ONLY: t_gridref_state
  USE mo_impl_constants,      ONLY: LEAPFROG_EXPL, LEAPFROG_SI, &
                                    RK4, SSPRK54, MAX_CHAR_LENGTH
  USE mo_timer,               ONLY: timer_total, timer_start, timer_stop
  USE mo_sync,                ONLY: global_max
  USE mo_vertical_coord_table,ONLY: vct
  USE mo_io_restart,          ONLY: write_restart_info_file

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

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
    SELECT CASE (itime_scheme)
    CASE (LEAPFROG_SI) ; CALL init_si_params
    CASE (RK4,SSPRK54) ; CALL init_RungeKutta(itime_scheme)
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

    SELECT CASE (itime_scheme)
    CASE (LEAPFROG_EXPL, LEAPFROG_SI, RK4, SSPRK54)
    ! 3 time level scheme, or 2 time level scheme with multi stages
      ntl = 3 + nadd
    CASE DEFAULT
      ntl = 2 + nadd
    END SELECT

    CALL construct_icoham_dyn_state( ntl, p_patch )

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
    IF (ltestcase) CALL init_testcase( p_patch, p_hydro_state, p_int_state)

    ! Convert temperature to potential temperature if the latter is
    ! chosen as a prognostic variable.

    IF (ltheta_dyn) THEN
      DO jg = 1, n_dom

        CALL convert_t2theta(p_patch(jg),                      &
                             p_hydro_state(jg)%prog(nnow(jg)), &
                             p_hydro_state(jg)%diag           )

        IF (.NOT.ltwotime)                                     &
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
                           ltheta_dyn, ltransport        )    ! in

     CALL update_diag_state( p_hydro_state(jg)%prog_out,    &! in
                             p_patch(jg), p_int_state(jg),  &! in
                             ext_data(jg),                  &! in
                             p_hydro_state(jg)%diag_out )    ! out

     CALL update_dyn_output( p_patch(jg), p_int_state(jg),  &! in
                             p_hydro_state(jg)%prog_out,    &! in
                             p_hydro_state(jg)%diag_out )    ! inout

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
  SUBROUTINE perform_ha_stepping( p_patch, p_int_state, p_grf_state,  &
                                & p_hydro_state, datetime,            &
                                & n_io, n_file, n_checkpoint, n_diag, &
                                & jfile, l_have_output                )

  CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = 'mo_ha_stepping:perform_ha_stepping'

  TYPE(t_patch), TARGET, INTENT(IN)         :: p_patch(n_dom)
  TYPE(t_int_state), TARGET, INTENT(IN)     :: p_int_state(n_dom)
  TYPE(t_gridref_state), TARGET, INTENT(INOUT) :: p_grf_state(n_dom)
  INTEGER, INTENT(IN) :: n_io, n_file, n_checkpoint, n_diag
  INTEGER, INTENT(INOUT) :: jfile
  LOGICAL, INTENT(INOUT) :: l_have_output

  TYPE(t_datetime), INTENT(INOUT)         :: datetime
  TYPE(t_hydro_atm), TARGET, INTENT(INOUT):: p_hydro_state(n_dom)

  REAL(wp), DIMENSION(:,:,:), POINTER     :: p_vn  => NULL()
  REAL(wp) :: sim_time(n_dom)
  INTEGER  :: jg, jn, jgc, jstep
  LOGICAL  :: l_3tl_init(n_dom)

#ifdef _OPENMP
  INTEGER  :: jb
  REAL(wp) :: vnmax, vn_aux(p_patch(1)%nblks_e)
#else
  REAL(wp) :: vnmax(1)
#endif


!-----------------------------------------------------------------------
  sim_time(:) = 0._wp

  IF (ltimer) CALL timer_start(timer_total)

  TIME_LOOP: DO jstep = 1, nsteps

    !--------------------------------------------------------------------------
    ! Send to stdout information about the current integration cycle
    ! (maximum normal wind).
    !--------------------------------------------------------------------------

    IF (msg_level >= 5) THEN ! compute maximum velocity in global model domain
      p_vn => p_hydro_state(1)%prog(nnow(1))%vn
      vnmax = -HUGE(vnmax)
#if (!defined __xlC__) && (defined _OPENMP)
!$OMP PARALLEL DO PRIVATE(jb)
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
    ! Set output flags
    !--------------------------------------------------------------------------
    ! The integration cycle "jstep" computes the atmospheric state variables
    ! at the new time step (n+1) using the already known step n (and possibly
    ! n-1, n-2, ... as well, depending on the user's choice of time stepping
    ! scheme). In some schemes (e.g. leapfrog), the step n values will be
    ! modified within this integration cycle. Therefore at the end of the
    ! cycle, we write out variables of time step n rather than n+1.

    IF ( MOD(jstep-1,n_io)==0 .AND. jstep/=1 ) THEN
      l_outputtime       = .TRUE. ! Output at the end of the time step
      lprepare_output(:) = .TRUE. ! Prepare output values on all grid levels
    ELSE
      l_outputtime       = .FALSE.
      lprepare_output(:) = .FALSE.
    ENDIF

    IF ( MOD(jstep,n_checkpoint)==0 .AND. jstep/=nsteps ) THEN
      l_checkpoint_time = .TRUE.
    ELSE
      l_checkpoint_time = .FALSE.
    ENDIF

    IF ( MOD(jstep-1,n_diag)==0 .OR. jstep==nsteps ) THEN
      l_diagtime = .TRUE. ! Diagnostic output is done at the end of the
    ELSE                  ! time step
      l_diagtime = .FALSE.
    ENDIF
    IF(.NOT.ldynamics.AND.lshallow_water)l_diagtime=.FALSE.

    !--------------------------------------------------------------------------
    ! Time integration from time step n to n+1
    !--------------------------------------------------------------------------
    ! Whether to apply a special initial treatment for 
    ! 3-time-level schemes. In case of a restart run, the treatment is 
    ! not necessary, thus the variable is set to .TRUE. 

    l_3tl_init(1:n_dom) = (.NOT.ltwotime).AND.(.NOT.lrestart).AND.(jstep==1)

    ! Call recursive subroutine 'process_grid_level', which executes
    ! one timestep for the global domain and calls itself in the presence
    ! of nested domains with recursively halved time steps

    CALL process_grid( p_patch, p_hydro_state, p_int_state, p_grf_state,    &
      &                1, jstep, l_3tl_init, dtime, sim_time, 1, datetime )

    !--------------------------------------------------------------------------
    ! Write output (prognostic and diagnostic variables) 
    !--------------------------------------------------------------------------
    IF (l_outputtime) THEN

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

      CALL write_output( datetime )
      CALL message(TRIM(routine),'Output at:')
      CALL print_datetime(datetime)
      l_have_output = .TRUE.

    ENDIF !l_outputtime

    ! If it's time, close the current output file and trigger a new one

    IF (jstep/=1.AND.(MOD(jstep-1,n_file)==0).AND.jstep/=nsteps) THEN
      jfile = jfile +1
      CALL init_output_files(jfile,lclose=l_have_output)
    ENDIF

    !--------------------------------------------------------------------------
    ! Diagnose global integrals
    !--------------------------------------------------------------------------
    IF (l_diagtime) &
    CALL supervise_total_integrals( jstep, p_patch, p_hydro_state, nnow )

    !--------------------------------------------------------------------------
    ! One integration cycle finished on the lowest grid level (coarsest
    ! resolution). Set model time.
    !--------------------------------------------------------------------------
    CALL add_time(dtime,0,0,0,datetime)

    !--------------------------------------------------------------------------
    ! Write restart file
    !--------------------------------------------------------------------------
    IF (l_checkpoint_time) THEN
      DO jg = 1, n_dom
        CALL create_restart_file( p_patch(jg), datetime, vct,        &
                                & jfile, l_have_output              )
      END DO

      ! Create the master (meta) file in ASCII format which contains
      ! info about which files should be read in for a restart run.

      CALL write_restart_info_file
    END IF

  ENDDO TIME_LOOP

  IF (ltimer) CALL timer_stop(timer_total)

  END SUBROUTINE perform_ha_stepping

END MODULE mo_ha_stepping

