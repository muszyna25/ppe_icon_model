#ifdef __xlC__
!@PROCESS nosmp
!@PROCESS NOOPTimize
#endif
!>
!!  This module contains the routines needed for managing flow control.
!!
!!  This module contains the routines needed for managing flow control
!!  with mesh refinement. The main routine, process_level, is a recursive
!!  subroutine that is called from sw_atmos for the global mesh and calls
!!  itself recursively for the refined meshes. It contains the whole time
!!  stepping management that was previously located in sw_atmos.
!!  Further subroutines serve for interpolating the time tendencies to the
!!  lateral boundaries of the refined meshes, for smoothing these fields,
!!  and for feedback.
!!
!! @par Revision History
!!  Developed and tested by Guenther Zaengl, DWD (2008-07)
!!  Modified by Marco Giorgetta, MPI-M (2009-02-26)
!!  - renamed ltracer to ltransport
!!  Modification by Guenther Zaengl, DWD (2009-06-22)
!!  - preparation for generalized grid refinement (affects all subroutines)
!!  Modification by Daniel Reinert, DWD (2010-02-22)
!!  - call of modified transport interface including call of prepare_tracer
!!    for each time integration scheme. Included new local variables for
!!    horizontal and vertical mass fluxes, velocities and pressure values
!!    at half and full levels.
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

MODULE mo_hierarchy_management

  USE mo_kind,                ONLY: wp, dp
  USE mo_exception,           ONLY: message_text, message, finish
  USE mo_model_domain,        ONLY: t_patch
  USE mo_ext_data_types,      ONLY: t_external_data
  USE mo_grid_config,         ONLY: n_dom, lfeedback, l_limited_area
  USE mo_intp_data_strc,      ONLY: t_int_state
  USE mo_grf_intp_data_strc,  ONLY: t_gridref_state
  USE mo_gridref_config,      ONLY: grf_intmethod_c, grf_intmethod_ct
  USE mo_grf_bdyintp,         ONLY: interpol_scal_grf
  USE mo_dynamics_config,     ONLY: ltwotime, lshallow_water,            &
                                    nold, nnow, nnew, nsav1, nsav2
  USE mo_ha_dyn_config,       ONLY: ha_dyn_config 
  USE mo_diffusion_config,    ONLY: diffusion_config
!   USE mo_io_config,           ONLY: lprepare_output
  USE mo_parallel_config,     ONLY: nproma, p_test_run
  USE mo_run_config,          ONLY: ldynamics, ltransport, &
    &                               nlev, nlevp1, ntracer, iforcing, lforcing
  USE mo_icoham_dyn_types,    ONLY: t_hydro_atm
  USE mo_ha_prog_util,        ONLY: update_prog_state, diag_tend !, copy_prog_state
  USE mo_ha_diag_util,        ONLY: update_pres !, update_diag_state, update_dyn_output
  USE mo_advection_stepping,  ONLY: step_advection
  USE mo_ha_leapfrog,         ONLY: step_leapfrog_expl, asselin, &
    &                               leapfrog_update_prog
  USE mo_ha_rungekutta,       ONLY: step_rungekutta
  USE mo_ha_2tl_si,           ONLY: step_2tl_si
  USE mo_hdiff,               ONLY: hdiff_expl
  USE mo_hdiff_hyb_lin,       ONLY: hdiff_hyb_lin
  USE mo_si_correction,       ONLY: si_correction
  USE mo_pa_test,             ONLY: set_vertical_velocity
  USE mo_sv_test,             ONLY: get_sv_tracer
  USE mo_df_test,             ONLY: get_df_velocity, get_departure_points,   &
    &                               prep_departure_points_err
  USE mo_ha_testcases,        ONLY: ctest_name,rotate_axis_deg
  USE mo_impl_constants,      ONLY: success, MAX_CHAR_LENGTH
  USE mo_expensive_functions, ONLY: convert_t2theta_lin, convert_theta2t_lin
!!$  USE mo_interface_icoham_echam, ONLY: interface_icoham_echam
  USE mo_loopindices,         ONLY: get_indices_c, get_indices_e
  USE mo_impl_constants_grf,  ONLY: grf_bdywidth_c, grf_bdywidth_e
  USE mo_impl_constants,      ONLY: tracer_only, two_tl_si,         &
    &                               leapfrog_expl, leapfrog_si,     &
    &                               rk4, ssprk54,                   &
    &                               min_rlcell_int, min_rledge_int, &
    &                               iheldsuarez, iecham, ildf_echam,& 
    &                               ildf_dry
  USE mo_ha_dtp_interface,    ONLY: prepare_tracer, prepare_tracer_rk, &
    &                               prepare_tracer_leapfrog!!$, prepare_echam_phy
  USE mo_held_suarez_interface, ONLY: held_suarez_interface
  USE mo_hierarchy_management_intp
  USE mo_mpi,                 ONLY: push_glob_comm, pop_glob_comm, proc_split
  USE mo_ldf_test,            ONLY: ldf_temp
  USE mtime,                  ONLY: datetime, timeDelta, MAX_DATETIME_STR_LEN, &
    &                               newTimeDelta, newDatetime, OPERATOR(*),    &
    &                               OPERATOR(+), deallocateDatetime,           &
    &                               deallocateTimedelta, datetimeToString,     &
    &                               getTimedeltaFromDatetime,                  &
    &                               getTotalMillisecondsTimedelta
  USE mo_time_config,         ONLY: time_config

  IMPLICIT NONE
  PRIVATE

  ! Local variables
  REAL(wp),POINTER :: z_mflx_me(:,:,:)      !< mass flux at full level edges
  !                                             !< at time step n+\alpha
  REAL(wp),ALLOCATABLE :: z_mflx_ic(:,:,:)      !< mass flux at half level centers
  !                                             !< at time step n+\alpha
  REAL(wp),ALLOCATABLE :: z_vn_traj(:,:,:)      !< horizontal velocity at edges for
  !                                             !< calculation of backward trajectories
  !                                             !< at time step n+0.5
  REAL(wp),ALLOCATABLE :: z_omega_traj (:,:,:)  !< vertical velocity at half level center
  !                                             !< for calculation of backward trajectories
  !                                             !< at time step n+0.5
  REAL(wp),ALLOCATABLE :: z_delp_mc_now(:,:,:)  !< pressure thickness at full levels at time
  !                                             !< step n
  REAL(wp),ALLOCATABLE :: z_pres_mc_now(:,:,:)  !< full level pressure at time step n
  REAL(wp),ALLOCATABLE :: z_pres_ic_now(:,:,:)  !< half level pressure at time step n

  PUBLIC :: process_grid, interpolate_diagnostics

CONTAINS


  !>
  !!
  !! Manages time stepping of a given grid (model domain) including
  !! interpolation to the lateral boundary of the child domain (if present)
  !! and feedback to the parent domain.
  !!
  !! @par Revision History
  !! Developed  by Guenther Zaengl, DWD, 2008-04-15
  !! Revised 2008-10-24 for use in hydrostatic model by ??
  !! Revised for physics-dynamics coupling by Hui Wan, MPI-M (2009-02-15)
  !!
  RECURSIVE SUBROUTINE process_grid( p_patch, p_hydro_state,         &
                                   & p_int_state, p_grf_state,       &
                                   & ext_data, jg, nstep_global,     &
                                   & l_3tl_init, dt_loc, mtime_step, &
                                   & nsteps,                         &
                                   & mtime_current )

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = 'mo_hierarchy_management:process_grid'

    TYPE(t_patch),TARGET, INTENT(IN)           :: p_patch(n_dom)
    TYPE(t_hydro_atm),  TARGET,INTENT(INOUT)   :: p_hydro_state(n_dom)
    TYPE(t_int_state),TARGET,INTENT(IN)        :: p_int_state(n_dom)
    TYPE(t_gridref_state),INTENT(INOUT)        :: p_grf_state(n_dom)
    TYPE(t_external_data), INTENT(INOUT)       :: ext_data(n_dom)

    INTEGER, INTENT(IN)    :: jg           ! current grid level
    INTEGER, INTENT(IN)    :: nstep_global ! number of global time step
    INTEGER, INTENT(IN)    :: nsteps       ! number of time steps to be executed

    REAL(wp),        INTENT(IN) :: dt_loc          ! time step applicable to local grid level
    TYPE(timeDelta), POINTER    :: mtime_step      ! time step (mtime object)

    LOGICAL, INTENT(INOUT) :: l_3tl_init(n_dom)

    TYPE(datetime),   POINTER    :: mtime_current     ! current datetime (mtime)

    TYPE(datetime),   POINTER    :: grid_datetime

    INTEGER  :: ns, jb, jstep, jn, jgc
    INTEGER  :: n_old, n_now, n_new, n_sav1, n_sav2, n_temp
    INTEGER  :: nlen
    REAL(wp) :: dt_sub, zdtime, rdt_loc

    TYPE(timeDelta), POINTER :: mtime_step_sub

    REAL(wp),DIMENSION ( nproma, nlev, p_patch(jg)%nblks_c ) :: temp_save
    REAL(wp), ALLOCATABLE :: z_pres_sfc(:,:,:)

    REAL(wp),DIMENSION(:,:),    POINTER :: p_psfc     => NULL()
    REAL(wp),DIMENSION(:,:),    POINTER :: p_psfc_sv  => NULL()
    REAL(wp),DIMENSION(:,:,:),  POINTER :: p_vn       => NULL()
    REAL(wp),DIMENSION(:,:,:),  POINTER :: p_vn_sv    => NULL()
    REAL(wp),DIMENSION(:,:,:),  POINTER :: p_temp     => NULL()
    REAL(wp),DIMENSION(:,:,:),  POINTER :: p_temp_sv  => NULL()
    REAL(wp),DIMENSION(:,:,:),  POINTER :: p_theta    => NULL()
    REAL(wp),DIMENSION(:,:,:),  POINTER :: p_theta_sv => NULL()
    REAL(wp),DIMENSION(:,:,:,:),POINTER :: p_trac     => NULL()
    REAL(wp),DIMENSION(:,:,:,:),POINTER :: p_trac_sv  => NULL()

    REAL(dp) :: tu, ts
    REAL(dp) :: t_tot, t_0, t_1
    INTEGER :: iret
    INTEGER :: nblks_e, nblks_c
    INTEGER :: ist
    INTEGER, EXTERNAL :: util_cputime
    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: dstring
    TYPE(timeDelta), POINTER            :: time_diff
    REAL(wp)                            :: sim_time     !< elapsed simulation time on this grid level

    ! calculate elapsed simulation time in seconds
    time_diff  => newTimedelta("PT0S")
    time_diff  =  getTimeDeltaFromDateTime(mtime_current, time_config%tc_startdate)
    sim_time   =  getTotalMillisecondsTimedelta(time_diff, mtime_current)*1.e-3_wp
    CALL deallocateTimedelta(time_diff)

    !---------------

    nblks_c = p_patch(jg)%nblks_c
    nblks_e = p_patch(jg)%nblks_e

    ! set up datetime variable valid within process_grid.
    grid_datetime => newDatetime(mtime_current)

    DO ns = 1, nsteps

      ! information on global time step, grid index and local time step
      WRITE (message_text,'(a,i6,a,i2,a,i2)')                  &
        & '.        step :',nstep_global,'  grid :',jg,' substep :',ns
      CALL message(TRIM(routine),message_text)
      CALL datetimeToString(mtime_current, dstring)
      CALL message(TRIM(routine),TRIM(dstring))

      iret = util_cputime(tu, ts)
      t_0 = tu+ts
      t_tot = 0.0_dp

      ALLOCATE( z_mflx_me    ( nproma,nlev  ,nblks_e ), &
        &       z_mflx_ic    ( nproma,nlevp1,nblks_c ), &
        &       z_vn_traj    ( nproma,nlev  ,nblks_e ), &
        &       z_omega_traj ( nproma,nlevp1,nblks_c ), &
        &       z_delp_mc_now( nproma,nlev  ,nblks_c ), &
        &       z_pres_mc_now( nproma,nlev  ,nblks_c ), &
        &       z_pres_ic_now( nproma,nlevp1,nblks_c ), &
        &       STAT=ist )
        
      IF (p_test_run) THEN
        z_mflx_me    ( : , : , : ) = 0.0_wp
        z_mflx_ic    ( : , : , : ) =  0.0_wp
        z_vn_traj    ( : , : , : ) =  0.0_wp
        z_omega_traj ( : , : , : ) =  0.0_wp
        z_delp_mc_now( : , : , : ) =  0.0_wp
        z_pres_mc_now( : , : , : ) =  0.0_wp
        z_pres_ic_now( : , : , : ) =  0.0_wp
      ENDIF
      
      IF (ist /= success) THEN
        CALL finish ( 'mo_hierarchy_management: process_grid',                &
          &           'allocation for z_mflx_me, z_mflx_ic, z_vn_traj, '    //&
          &           'z_omega_traj, z_delp_mc_now, z_pres_mc_now,     '    //&
          &           'z_pres_ic_now failed' )
      ENDIF

      IF ( jg < n_dom ) THEN
        ! Save prognostic variables at current timestep to compute
        ! nest boundary tendencies (not needed at highest nest level)
        nsav1(jg) = nnow(jg)
      ENDIF



      IF ( ns==1 .AND. jg>1 ) THEN !the 1st substep in a refined domain

        ! Save prognostic variables at current timestep to compute
        ! feedback increments (not needed in global domain)
        n_now = nnow(jg)
        n_sav2 = nsav2(jg)
        ! The pointers are needed because the subsequent parallelization
        ! causes trouble otherwise
        p_psfc    => p_hydro_state(jg)%prog(n_now)%pres_sfc
        p_psfc_sv => p_hydro_state(jg)%prog(n_sav2)%pres_sfc
        p_vn      => p_hydro_state(jg)%prog(n_now)%vn
        p_vn_sv   => p_hydro_state(jg)%prog(n_sav2)%vn
        p_temp    => p_hydro_state(jg)%prog(n_now)%temp
        p_temp_sv => p_hydro_state(jg)%prog(n_sav2)%temp
        IF (ha_dyn_config%ltheta_dyn) THEN
          p_theta   => p_hydro_state(jg)%prog(n_now)%theta
          p_theta_sv=> p_hydro_state(jg)%prog(n_sav2)%theta
        ENDIF
        IF(ntracer > 0) THEN
          p_trac    => p_hydro_state(jg)%prog(n_now)%tracer
          p_trac_sv => p_hydro_state(jg)%prog(n_sav2)%tracer
        ENDIF

!$OMP PARALLEL
!$OMP WORKSHARE
        p_psfc_sv = p_psfc
        p_vn_sv   = p_vn
        p_temp_sv = p_temp
!$OMP END WORKSHARE

        IF (ha_dyn_config%ltheta_dyn) THEN
!$OMP WORKSHARE
          p_theta_sv = p_theta
!$OMP END WORKSHARE
        ENDIF

        IF(ntracer > 0) THEN
!$OMP WORKSHARE
          p_trac_sv = p_trac
!$OMP END WORKSHARE
        ENDIF
!$OMP END PARALLEL

      ENDIF !the 1st substep in a refined domain


      ! LL : These values should be already initialized.
      !      In any case, they cannot be initilized here
      !       when running in restart mode
!       IF ((nstep_global == 1) .AND. (jg == 1)) THEN
! 
!         ! For the 1st step of on the coarsest grid level, set all tendencies
!         ! to zero
! 
!         p_psfc => p_hydro_state(jg)%tend_dyn%pres_sfc
!         p_vn   => p_hydro_state(jg)%tend_dyn%vn
!         p_temp => p_hydro_state(jg)%tend_dyn%temp
!         IF (ntracer > 0) THEN
!           p_trac => p_hydro_state(jg)%tend_dyn%tracer
!         ENDIF
! 
! !$OMP PARALLEL
! !$OMP WORKSHARE
!         p_psfc = 0._wp
!         p_vn   = 0._wp
!         p_temp = 0._wp
! !$OMP END WORKSHARE
!         IF (ntracer > 0) THEN
! !$OMP WORKSHARE
!           p_trac = 0._wp
! !$OMP END WORKSHARE
!         ENDIF
! !$OMP END PARALLEL
! 
!       ENDIF !(nstep_global == 1) .AND. (jg == 1)

      IF ( l_limited_area .AND. (jg == 1) ) THEN
        ! Perform interpolation of lateral boundary tendencies from external data
        ! Currently, this is only a dummy routine that sets the boundary tendencies
        ! to zero
        CALL boundary_tendencies ( p_patch(jg), p_hydro_state(jg) )
      ENDIF

      ! jstep is used for step_advection; it needs to alternate between
      ! even and odd numbers for subsequent time steps
      IF (jg == 1) THEN
        jstep = nstep_global
      ELSE
        jstep = ns
      ENDIF

      IF ( (nstep_global==1) .AND. l_3tl_init(jg) ) THEN
        !==========================================================================
        ! Special treatment for 3 time level schemes
        !==========================================================================

        CALL message(TRIM(routine),' special treatment for 3 time level schemes')

        WRITE(message_text,'(a,i10)') 'TIME STEP n: ', nstep_global
        CALL message(TRIM(routine),message_text)

        CALL leapfrog_startup( ha_dyn_config%ileapfrog_startup,            &
          &                    p_patch, p_int_state, ext_data, jg, dt_loc, &
          &                    p_hydro_state, n_old, n_now, n_new, n_sav1  )

        l_3tl_init(jg) = .FALSE.

        zdtime = dt_loc !needed in case "update_prog_state" is called in the "fast physics" part

      ELSE !no special treatment for later time steps or 2 time level schemes
        !==========================================================================
        ! Now start the normal integration steps
        !==========================================================================

        n_old  = nold(jg)
        n_now  = nnow(jg)
        n_new  = nnew(jg)
        n_sav1 = nsav1(jg)

        zdtime = dt_loc

        SELECT CASE (ha_dyn_config%itime_scheme)
          !------------------
          ! Pure advection
          !------------------
        CASE (tracer_only)

          SELECT CASE ( TRIM(ctest_name) )

          CASE ('PA') ! solid body rotation
            IF (.NOT.lshallow_water) THEN
              ! set time-variant vertical velocity
              CALL set_vertical_velocity( p_patch(jg), p_hydro_state(jg)%diag,  &
                &                         jstep, sim_time )
            ENDIF

          CASE ('SV') ! static vortex
            ! since the second tracer field is only allocated in order to store
            ! the analytic reference solution, accessing 'step_advection' with
            ! ntracer=2 would lead to an unnecessary overhead. This is avoided
            ! by setting ntracer=1 before and ntracer=2 after the call of
            ! 'step_advection'
            ntracer = 1

          CASE ('DF1', 'DF2', 'DF3', 'DF4') ! deformational flow

            ! get velocity field
            CALL get_df_velocity( p_patch(jg), p_hydro_state(jg)%prog(n_new),      &
              &                   ctest_name, rotate_axis_deg, sim_time+zdtime )

            ! compute barycenter of departure region
            CALL get_departure_points( p_patch(jg), p_hydro_state(jg)%prog(n_new), &
              &                        ctest_name, rotate_axis_deg, zdtime,        &
              &                        sim_time+zdtime )
          END SELECT

          ! Diagnose some pressure- and velocity-related quantities for
          ! the tracer transport scheme
          CALL prepare_tracer( p_patch(jg), p_int_state(jg),   &! in
            &                  p_hydro_state(jg)%prog(n_now),  &! in
            &                  p_hydro_state(jg)%prog(n_new),  &! in
            &                  ha_dyn_config%itime_scheme,     &! in
            &                  ha_dyn_config%si_2tls,          &! in
            &                  p_hydro_state(jg)%diag,         &! inout
            &                  z_mflx_me, z_vn_traj,           &! out
            &                  z_mflx_ic, z_omega_traj,        &! out
            &                  z_delp_mc_now,                  &! out
            &                  z_pres_mc_now, z_pres_ic_now    )! out


          ! tracer advection
          CALL step_advection( p_patch(jg), p_int_state(jg), zdtime, jstep, & !in
            &                  p_hydro_state(jg)%prog(n_now)%tracer,        & !in
            &                  z_mflx_me, z_vn_traj,                        & !in
            &                  p_hydro_state(jg)%diag%weta,                 & !in
            &                  p_hydro_state(jg)%diag%weta,                 & !in
            &                  z_delp_mc_now,                               & !in
            &                  p_hydro_state(jg)%diag%delp_c,               & !in
            &                  z_delp_mc_now,                               & !in
            &                  p_hydro_state(jg)%tend_dyn%tracer,           & !inout
            &                  p_hydro_state(jg)%prog(n_new)%tracer,        & !inout
            &                  p_hydro_state(jg)%diag%hfl_tracer,           & !out
            &                  p_hydro_state(jg)%diag%vfl_tracer )            !out


          SELECT CASE ( TRIM(ctest_name) )

          CASE ('SV') ! static vortex

            ntracer = 2
            ! save analytic solution (tracer field) for error assessment
            CALL get_sv_tracer( p_patch(jg),                                   &
              &                 p_hydro_state(jg)%prog(n_new)%tracer(:,:,:,2), &
              &                 sim_time )

          CASE ('DF1', 'DF2', 'DF3', 'DF4') ! deformational flow

            ! compute error of trajectory calculation and store in additional 
            ! (dummy) tracer field (sum over cell edges)
            CALL prep_departure_points_err( p_patch(jg),                                  &
              &                      p_hydro_state(jg)%prog(n_new)%tracer(:,:,:,ntracer) )

          END SELECT

          !---------------------------------------
          ! Semi-implicit two time level scheme
          !---------------------------------------
        CASE (TWO_TL_SI)

          !! Note that for 2 time level schemes, time step n and n+1 
          !! are referred to as "n_now" and "n_new", respectively.
          !!
          IF ( iforcing==iheldsuarez ) THEN

             CALL update_pres( p_hydro_state(jg)%prog(n_now),  &! in
               &               p_patch(jg), p_int_state(jg),   &! in
               &               p_hydro_state(jg)%diag        )  ! inout

             CALL held_suarez_interface( p_patch(jg), p_int_state(jg),   &! in
               &                         p_hydro_state(jg)%prog(n_now),  &! in
               &                         p_hydro_state(jg)%diag,         &! in
               &                         p_hydro_state(jg)%tend_phy    )  ! inout

             CALL update_prog_state( zdtime, p_patch(jg),             &! in
               &                     p_hydro_state(jg)%tend_phy,      &! in
               &                     .FALSE.,ha_dyn_config%ltheta_dyn,&! in
               &                     .FALSE.,                         &! in
               &                     p_hydro_state(jg)%prog(n_now)   ) ! inout

           ENDIF !( iforcing==iheldsuarez )          

          
          ! Dynamical core
          CALL step_2tl_si( ha_dyn_config%si_expl_scheme,    &! in
            &               ha_dyn_config%si_2tls,           &! in
            &               ha_dyn_config%si_rtol,           &! in
            &               ha_dyn_config%ltheta_dyn,        &! in
            &               zdtime,                          &! in
            &               p_patch(jg), p_int_state(jg),    &! in
            &               p_hydro_state(jg)%prog(n_now),   &! in
            &               ext_data(jg),                    &! in
            &               p_hydro_state(jg)%prog(n_new),   &! inout
            &               p_hydro_state(jg)%diag,          &! out
            &               p_hydro_state(jg)%tend_dyn     )  ! inout

          ! tracer advection
          IF ( ltransport ) THEN

            ! Diagnose some pressure- and velocity-related quantities for
            ! the tracer transport scheme
            CALL prepare_tracer( p_patch(jg), p_int_state(jg),    &! in
              &           p_hydro_state(jg)%prog(n_now),          &! in
              &           p_hydro_state(jg)%prog(n_new),          &! in
              &           ha_dyn_config%itime_scheme,             &! in
              &           ha_dyn_config%si_2tls,                  &! in
              &           p_hydro_state(jg)%diag,                 &! inout
              &           z_mflx_me, z_vn_traj,                   &! out
              &           z_mflx_ic, z_omega_traj,                &! out
              &           z_delp_mc_now,                          &! out
              &           z_pres_mc_now, z_pres_ic_now            )! out


            CALL step_advection( p_patch(jg), p_int_state(jg), zdtime, jstep,&! in
              &                 p_hydro_state(jg)%prog(n_now)%tracer,        &! in
              &                 z_mflx_me, z_vn_traj,                        &! in
              &                 z_mflx_ic,                                   &! in
              &                 z_mflx_ic,                                   &! in
              &                 z_delp_mc_now,                               &! in
              &                 p_hydro_state(jg)%diag%delp_c,               &! in
              &                 z_delp_mc_now,                               &! in
              &                 p_hydro_state(jg)%tend_dyn%tracer,           &! inout*
              &                 p_hydro_state(jg)%prog(n_new)%tracer,        &! inout
              &                 p_hydro_state(jg)%diag%hfl_tracer,           &! out
              &                 p_hydro_state(jg)%diag%vfl_tracer            )! out

             !Note: the argument * is used only in case of grid refinement.
          ENDIF

!!$          ! ECHAM physics. Meteorological fields of time step n are provided to the 
!!$          ! the physics package, which also takes over the dynamics/transport-induced 
!!$          ! tendencies (diagnosed in subroutine "diag_tend" and stored in the 
!!$          ! the variable %tend_phy), use them in the parameterization schemes, 
!!$          ! and add new contributions to them.
!!$
!!$          IF (iforcing==iecham) THEN
!!$
!!$            CALL diag_tend( p_patch(jg),                           &! in
!!$              &             ha_dyn_config%ltheta_dyn,              &! in
!!$              &             ltransport, zdtime,                    &! in
!!$              &             p_hydro_state(jg)%prog(n_now),         &! in
!!$              &             p_hydro_state(jg)%prog(n_new),         &! in
!!$              &             p_hydro_state(jg)%tend_phy           )  ! inout
!!$
!!$            CALL prepare_echam_phy( ha_dyn_config%ltheta_dyn,      &! in
!!$              &                     p_patch(jg), p_int_state(jg),  &! in
!!$              &                     ext_data(jg),                  &! in
!!$              &                     p_hydro_state(jg)%prog(n_now), &! in
!!$              &                     p_hydro_state(jg)%diag       )  ! inout
!!$
!!$            CALL interface_icoham_echam( zdtime, zdtime,                &! in
!!$              &                          grid_datetime,                 &! in
!!$              &                          p_patch(jg),                   &! in
!!$              &                          p_int_state(jg),               &! in
!!$              &                          p_hydro_state(jg)%prog(n_now), &! inout
!!$              &                          p_hydro_state(jg)%diag,        &! in
!!$              &                          p_hydro_state(jg)%prog(n_new), &! in
!!$              &                          p_hydro_state(jg)%tend_phy )    ! inout
!!$
!!$
!!$            ! At this point the %tend_phy state contains the accumulated tendencies.
!!$            ! Get the prognostic state at time step n+1 ("n_new").
!!$
!!$            CALL update_prog_state( zdtime, p_patch(jg),             &! in
!!$              &                     p_hydro_state(jg)%tend_phy,      &! in
!!$              &                     .FALSE.,ha_dyn_config%ltheta_dyn,&! in
!!$              &                     ltransport,                      &! in
!!$              &                     p_hydro_state(jg)%prog(n_now),   &! in
!!$              &                     p_hydro_state(jg)%prog(n_new)    )! inout
!!$
!!$          ENDIF ! iforcing==iecham

          ! Horizontal diffusion is called later in this subroutine.

          !------------------------------------------
          ! Leapfrog schemes
          !------------------------------------------
        CASE (leapfrog_expl, leapfrog_si)

          ! Because of the Asselin filter, the only filtered time step
          ! at this point is n-1, which should be used for
          ! the slow processes in the physics parameterisations.
          ! For the Held-Suarez test we couple physics and dynamics in the
          ! "time-split" manner (sequential splitting), the prognostic state
          ! n_old is updated before calling the dynamics.

          IF ( iforcing==iheldsuarez ) THEN

            CALL update_pres( p_hydro_state(jg)%prog(n_old),  &! in
              &               p_patch(jg), p_int_state(jg),   &! in
              &               p_hydro_state(jg)%diag        )  ! inout

            CALL held_suarez_interface( p_patch(jg), p_int_state(jg), &! in
              &               p_hydro_state(jg)%prog(n_old),          &! in
              &               p_hydro_state(jg)%diag,                 &! in
              &               p_hydro_state(jg)%tend_phy    )          ! inout

            CALL update_prog_state( 2._wp*zdtime, p_patch(jg),       &! in
              &                     p_hydro_state(jg)%tend_phy,      &! in
              &                     .FALSE.,ha_dyn_config%ltheta_dyn,&! in
              &                     .FALSE.,                         &! in
              &                     p_hydro_state(jg)%prog(n_old)   ) ! inout
            
          ENDIF !( iforcing==iheldsuarez )

          ! First part of the dynamical core: the explicit leapfrog scheme.
          ! Compute tendencies of time step n, and the preliminary n+1 values.

          CALL step_leapfrog_expl( zdtime, zdtime,                   &! in
            &                      ha_dyn_config%ltheta_dyn,         &! in
            &                      p_patch(jg), p_int_state(jg),     &! in
            &                      p_hydro_state(jg)%prog(n_old),    &! in
            &                      ext_data(jg),                     &! in
            &                      p_hydro_state(jg)%prog(n_now),    &! in
            &                      p_hydro_state(jg)%diag,           &! out
            &                      p_hydro_state(jg)%prog(n_new),    &! inout
            &                      p_hydro_state(jg)%tend_dyn     )   ! inout

          !Local Diabatic Forcing Testcase
          IF( iforcing==ildf_dry .OR. iforcing==ildf_echam ) THEN

            CALL update_pres( p_hydro_state(jg)%prog(n_old),    &! in
              &               p_patch(jg), p_int_state(jg),   &! in
              &               p_hydro_state(jg)%diag        )  ! inout

            CALL ldf_temp   ( p_patch(jg),                            &! in
              &               p_hydro_state(jg)%prog(n_old),          &! in
              &               p_hydro_state(jg)%diag,                 &! in
              &               p_hydro_state(jg)%tend_dyn    )          ! inout

            CALL leapfrog_update_prog( p_hydro_state(jg)%prog(n_new),    &! inout
              &                        p_hydro_state(jg)%prog(n_now),    &! in, for refinement
              &                        p_hydro_state(jg)%prog(n_old),    &! in
              &                        p_hydro_state(jg)%tend_dyn,       &! in
              &                        2._wp*zdtime, zdtime, ltransport, &! in
              &                        ha_dyn_config%ltheta_dyn,         &! in
              &                        p_patch(jg)                       )! in

          END IF ! (iforcing==ildf_dry) .OR. (iforcing==ildf_echam) )

          ! If running the model with ECHAM6 physics, perform tracer transport
          ! and call the physics package. The tracer tendencies and the 
          ! leapfrog tendencies will be passed on to the ECHAM physics package 
          ! and used therein. The purpose of this implementation is to get 
          ! a ICOHAM version similar to ECHAM.

!!$          IF ( iforcing==iecham .OR. iforcing==ildf_echam ) THEN
!!$
!!$            ! Tracer transport.
!!$            ! Like in ECHAM, the transport scheme computes the tracer specific
!!$            ! density at time step n+1 by using tracer and air density of n-1,
!!$            ! as well as wind fields (vn, weta) at step n. Note that at this 
!!$            ! point the surface pressure at step n+1 has not yet been 
!!$            ! corrected by the semi-implicit scheme, thus this implementation
!!$            ! guarantees neither mass conservation nor consistency with continuity.
!!$
!!$            IF (ltransport) THEN
!!$              CALL prepare_tracer_leapfrog( p_patch(jg),                     &! in
!!$                &                           p_hydro_state(jg)%prog(n_old),   &! in
!!$                &                           p_hydro_state(jg)%prog(n_now),   &! in
!!$                &                           p_hydro_state(jg)%prog(n_new),   &! in
!!$                &                           p_hydro_state(jg)%diag,          &! inout
!!$                &                           z_vn_traj,                       &! out ("now")
!!$                &                           z_delp_mc_now,                   &! out (in fact "old")
!!$                &                           z_pres_mc_now, z_pres_ic_now     )! out (in fact "old")
!!$
!!$
!!$              CALL step_advection( p_patch(jg), p_int_state(jg),               &!in
!!$                &                  2._wp*zdtime, jstep,                        &!in
!!$                &                  p_hydro_state(jg)%prog(n_old)%tracer,       &!in
!!$                &                  p_hydro_state(jg)%diag%mass_flux_e,         &!in ("now")
!!$                &                  z_vn_traj,                                  &!in ("now")
!!$                &                  p_hydro_state(jg)%diag%weta,                &!in ("now")
!!$                &                  p_hydro_state(jg)%diag%weta,                &!in ("now")
!!$                &                  z_delp_mc_now,                              &!in (in fact"old")
!!$                &                  p_hydro_state(jg)%diag%delp_c,              &!in ("new")
!!$                &                  z_delp_mc_now,                              &!in (in fact "old")
!!$                &                  p_hydro_state(jg)%tend_dyn%tracer,          &!in (ONLY for ref.)
!!$                &                  p_hydro_state(jg)%prog(n_new)%tracer,       &!inout
!!$                &                  p_hydro_state(jg)%diag%hfl_tracer,          &!out
!!$                &                  p_hydro_state(jg)%diag%vfl_tracer,          &!out
!!$                &                  opt_ddt_tracer_adv =                        &!optional dummy
!!$                &                  p_hydro_state(jg)%tend_dyn%tracer   )        !out
!!$            ELSE
!!$              p_hydro_state(jg)%tend_dyn%tracer = 0._wp   ! Q&D. Not appropriate for refinement
!!$            ENDIF
!!$
!!$            ! ECHAM physics. The meteorological fields of step n-1 are (Asselin) filtered values.
!!$            ! The physics package takes over the dynamics/transport-induced tendencies,
!!$            ! use them in the parameterization schemes, and add new contributions to them.
!!$
!!$            CALL prepare_echam_phy( ha_dyn_config%ltheta_dyn,      &! in
!!$              &                     p_patch(jg), p_int_state(jg),  &! in
!!$              &                     ext_data(jg),                  &! in
!!$              &                     p_hydro_state(jg)%prog(n_old), &! in
!!$              &                     p_hydro_state(jg)%diag         )! inout
!!$
!!$            CALL interface_icoham_echam( zdtime, 2._wp*zdtime,          &! in
!!$              &                          grid_datetime,                 &! in
!!$              &                          p_patch(jg),                   &! in
!!$              &                          p_int_state(jg),               &! in
!!$              &                          p_hydro_state(jg)%prog(n_old), &! inout
!!$              &                          p_hydro_state(jg)%diag,        &! in
!!$              &                          p_hydro_state(jg)%prog(n_new), &! in
!!$              &                          p_hydro_state(jg)%tend_dyn     )! inout
!!$
!!$            ! Explicit dynamics, tracer transport and physics computation have all
!!$            ! been finished. The %tend_dyn state contains the accumulated tendencies.
!!$            ! Now get the prognostic variables of time step n+1 using the leapfrog formula.
!!$            ! The n+1 values will later be corrected by the semi-implicit scheme.
!!$
!!$            CALL leapfrog_update_prog( p_hydro_state(jg)%prog(n_new),    &! inout
!!$              &                        p_hydro_state(jg)%prog(n_now),    &! in, for refinement
!!$              &                        p_hydro_state(jg)%prog(n_old),    &! in
!!$              &                        p_hydro_state(jg)%tend_dyn,       &! in
!!$              &                        2._wp*zdtime, zdtime, ltransport, &! in
!!$              &                        ha_dyn_config%ltheta_dyn,         &! in
!!$              &                        p_patch(jg)                       )! in
!!$
!!$          ENDIF ! (iforcing==iecham)

          ! Second part of the dynamical core: semi-implicit correction

          IF (ha_dyn_config%itime_scheme == leapfrog_si) THEN ! Add semi-implicit correction

            IF (ha_dyn_config%ltheta_dyn) THEN
              CALL convert_theta2t_lin( p_patch(jg),                   &
                &                       p_hydro_state(jg)%prog(n_now), &
                &                       p_hydro_state(jg)%prog(n_new), &
                &                       p_hydro_state(jg)%diag        )

              p_temp => p_hydro_state(jg)%prog(n_new)%temp
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,nlen) ICON_OMP_DEFAULT_SCHEDULE
              DO jb = 1,p_patch(jg)%nblks_c

                IF (jb /= p_patch(jg)%nblks_c) THEN
                  nlen = nproma
                ELSE
                  nlen = p_patch(jg)%npromz_c
                ENDIF

                temp_save(1:nlen,1:nlev,jb) = p_temp(1:nlen,1:nlev,jb)
              ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
            ENDIF !ltheta_dyn

            CALL si_correction( ha_dyn_config%lsi_3d,                 &! in
              &                 ha_dyn_config%si_coeff,               &! in 
              &                 ha_dyn_config%si_rtol,                &! in
              &                 lshallow_water,                       &! in
              &                 zdtime, p_patch(jg), p_int_state(jg), &! in
              &                 p_hydro_state(jg)%prog(n_old),        &! in
              &                 p_hydro_state(jg)%prog(n_now),        &! in
              &                 p_hydro_state(jg)%prog(n_new)         )! inout

            IF (ha_dyn_config%ltheta_dyn) THEN
              CALL convert_t2theta_lin( p_patch(jg),                         &
                &                       p_hydro_state(jg)%prog(n_new),       &
                &                       p_hydro_state(jg)%diag, temp_save )
            ENDIF

          ENDIF !with(out) semi-implicit correction

          ! Tracer transport (only when not using ECHAM physics).
          ! This implementation is mass conserving, but not consistent with continuity
          ! because the dry dynamics uses a 3-time-level scheme.

          IF ( (ltransport).AND.( iforcing/=iecham .AND. iforcing/=ildf_echam ) ) THEN

            ! Diagnose some pressure- and velocity-related quantities for
            ! the tracer transport scheme
            CALL prepare_tracer( p_patch(jg), p_int_state(jg),    &! in
              &                  p_hydro_state(jg)%prog(n_now),   &! in
              &                  p_hydro_state(jg)%prog(n_new),   &! in
              &                  ha_dyn_config%itime_scheme,      &! in
              &                  ha_dyn_config%si_2tls,           &! in
              &                  p_hydro_state(jg)%diag,          &! inout
              &                  z_mflx_me, z_vn_traj,            &! out
              &                  z_mflx_ic, z_omega_traj,         &! out
              &                  z_delp_mc_now,                   &! out
              &                  z_pres_mc_now, z_pres_ic_now     )! out


            CALL step_advection( p_patch(jg), p_int_state(jg), zdtime, jstep, &! in
              &                  p_hydro_state(jg)%prog(n_now)%tracer,        &! in
              &                  z_mflx_me, z_vn_traj,                        &! in
              &                  z_mflx_ic,                                   &! in
              &                  z_mflx_ic,                                   &! in
              &                  z_delp_mc_now,                               &! in
              &                  p_hydro_state(jg)%diag%delp_c,               &! in
              &                  z_delp_mc_now,                               &! in
              &                  p_hydro_state(jg)%tend_dyn%tracer,           &! inout
              &                  p_hydro_state(jg)%prog(n_new)%tracer,        &! inout
              &                  p_hydro_state(jg)%diag%hfl_tracer,           &! out
              &                  p_hydro_state(jg)%diag%vfl_tracer   )         ! out
          ENDIF

          ! Horizontal diffusion and fast processes in the physics package
          ! are called later in this subroutine before entering finer grid;
          ! Asselin filter is called after feedback

          !---------------------
          ! Runge-Kutta method
          !---------------------
        CASE (rk4,ssprk54)

          !! Note that for 2 time level schemes, time step n and n+1 
          !! are referred to as "n_now" and "n_new", respectively.

          IF ( iforcing==iheldsuarez ) THEN

            CALL update_pres( p_hydro_state(jg)%prog(n_now),  &! in
              &               p_patch(jg), p_int_state(jg),   &! in
              &               p_hydro_state(jg)%diag        )  ! inout

            CALL held_suarez_interface( p_patch(jg), p_int_state(jg),   &! in
              &                         p_hydro_state(jg)%prog(n_now),  &! in
              &                         p_hydro_state(jg)%diag,         &! in
              &                         p_hydro_state(jg)%tend_phy    )  ! inout

            ! At this point we've got the tendencies induced by the HS forcing.
            ! Here we couple physics and dynamics in the "time-split" manner 
            ! (sequential splitting). Thus the prognostic state
            ! n_now is updated here before calling the dynamics.

            CALL update_prog_state( zdtime, p_patch(jg),             &! in
              &                     p_hydro_state(jg)%tend_phy,      &! in
              &                     .FALSE.,ha_dyn_config%ltheta_dyn,&! in
              &                     .FALSE.,                         &! in
              &                     p_hydro_state(jg)%prog(n_now)   ) ! inout

          ENDIF ! iforcing==iheldsuarez

          ! Dynamical core
          CALL step_rungekutta( zdtime, ha_dyn_config%ltheta_dyn,     & ! in
            &                   p_patch(jg), p_int_state(jg),         & ! in
            &                   p_hydro_state(jg)%prog(n_now),        & ! in
            &                   ext_data(jg),                         & ! in
            &                   p_hydro_state(jg)%prog(n_old),        & ! tmp
            &                   p_hydro_state(jg)%prog(n_new),        & ! inout
            &                   p_hydro_state(jg)%diag,               & ! out
            &                   p_hydro_state(jg)%tend_dyn)             ! inout

          ! tracer advection
          IF ( ltransport ) THEN

            CALL prepare_tracer_rk( p_patch(jg), p_int_state(jg),      &! in
              &                     p_hydro_state(jg)%prog(n_now),     &! in
              &                     p_hydro_state(jg)%prog(n_new),     &! in
              &                     p_hydro_state(jg)%diag,            &! inout
              &                     z_vn_traj, z_delp_mc_now,          &! out
              &                     z_pres_mc_now, z_pres_ic_now       )! out


            CALL step_advection( p_patch(jg), p_int_state(jg), zdtime, jstep, &! in
              &                  p_hydro_state(jg)%prog(n_now)%tracer,        &! in
              &                  p_hydro_state(jg)%diag%mass_flux_e,          &! in
              &                  z_vn_traj,                                   &! in
              &                  p_hydro_state(jg)%diag%weta,                 &! in
              &                  p_hydro_state(jg)%diag%weta,                 &! in
              &                  z_delp_mc_now,                               &! in
              &                  p_hydro_state(jg)%diag%delp_c,               &! in
              &                  z_delp_mc_now,                               &! in
              &                  p_hydro_state(jg)%tend_dyn%tracer,           &! inout*
              &                  p_hydro_state(jg)%prog(n_new)%tracer,        &! inout
              &                  p_hydro_state(jg)%diag%hfl_tracer,           &! out
              &                  p_hydro_state(jg)%diag%vfl_tracer            )! out

          ENDIF

!!$          ! ECHAM physics. The physics package not only uses meteorological 
!!$          ! fields of step n, but also takes over the dynamics/transport-
!!$          ! induced tendencies (diagnosed in subroutine "diag_tend" and stored 
!!$          ! in the variable %tend_phy), use them in the parameterization schemes,
!!$          ! and add new contributions to them.
!!$          IF ( iforcing==iecham ) THEN
!!$
!!$            CALL diag_tend( p_patch(jg),                           &! in
!!$              &             ha_dyn_config%ltheta_dyn,              &! in
!!$              &             ltransport, zdtime,                    &! in
!!$              &             p_hydro_state(jg)%prog(n_now),         &! in
!!$              &             p_hydro_state(jg)%prog(n_new),         &! in
!!$              &             p_hydro_state(jg)%tend_phy           )  ! inout
!!$ 
!!$            CALL prepare_echam_phy( ha_dyn_config%ltheta_dyn,      &! in
!!$              &                     p_patch(jg), p_int_state(jg),  &! in
!!$              &                     ext_data(jg),                  &! in
!!$              &                     p_hydro_state(jg)%prog(n_now), &! in
!!$              &                     p_hydro_state(jg)%diag       )  ! inout
!!$
!!$            CALL interface_icoham_echam( zdtime, zdtime,                &! in
!!$              &                          grid_datetime,                 &! in
!!$              &                          p_patch(jg),                   &! in
!!$              &                          p_int_state(jg),               &! in
!!$              &                          p_hydro_state(jg)%prog(n_now), &! inout
!!$              &                          p_hydro_state(jg)%diag,        &! in
!!$              &                          p_hydro_state(jg)%prog(n_new), &! in
!!$              &                          p_hydro_state(jg)%tend_phy )    ! inout
!!$
!!$            CALL update_prog_state( zdtime, p_patch(jg),             &! in
!!$              &                     p_hydro_state(jg)%tend_phy,      &! in
!!$              &                     .FALSE.,ha_dyn_config%ltheta_dyn,&! in
!!$              &                     ltransport,                      &! in
!!$              &                     p_hydro_state(jg)%prog(n_now),   &! in
!!$              &                     p_hydro_state(jg)%prog(n_new)   ) ! inout
!!$
!!$          ENDIF ! iforcing==iecham

          ! Horizontal diffusion is called later in this subroutine 
          ! before entering finer grid.

        CASE DEFAULT
          CALL message(TRIM(routine),'wrong itime_scheme')
        END SELECT !( itime_scheme )

        !======================
        ! Horizontal diffusion
        !======================
        IF (ldynamics) THEN
          SELECT CASE(diffusion_config(jg)%hdiff_order)
          CASE(-1)
            ! No diffusion

          CASE(24,42)
            CALL hdiff_hyb_lin( jg, p_patch(jg), p_int_state(jg), &! in
                              & p_hydro_state(jg)%diag,           &! in
                              & p_hydro_state(jg)%prog(n_new)     )! inout

          CASE DEFAULT
            CALL hdiff_expl( p_patch(jg), p_hydro_state(jg)%prog(n_new), &
                           & p_hydro_state(jg)%diag, p_int_state(jg),    &
                           & zdtime, jg )
          END SELECT
        ENDIF

      ENDIF !if apply special treatment to the 1st step of a 3 time level scheme

      !=======================================
      ! Fast processes in the physics package
      !=======================================
      ! In case operator splitting is needed within the physics package,
      ! we calculate here the fast processes using the meteorological and
      ! tracer fields at time step n+1. The same is done for all time stepping
      ! schemes.

      IF ( lforcing ) THEN

        SELECT CASE (iforcing)
!!$        CASE(iecham,ildf_echam)
!!$          CALL prepare_physics( p_patch(jg), p_int_state(jg),  &! in
!!$            &                   p_hydro_state(jg)%prog(n_new), &! in
!!$            &                   p_hydro_state(jg)%diag       )  ! inout
!!$
!!$          CALL interface_icoham_echam( p_patch(jg), p_int_state(jg),  &! in
!!$            &                          p_hydro_state(jg)%prog(n_new), &! in
!!$            &                          p_hydro_state(jg)%diag,        &! in
!!$            &                          p_hydro_state(jg)%tend_phy )    ! inout
!!$
!!$           CALL update_prog_state( zdtime, p_patch(jg),              &! in
!!$             &                     p_hydro_state(jg)%tend_phy,       &! in
!!$             &                     .FALSE., ha_dyn_config%ltheta_dyn,&! in
!!$             &                     ltransport,                       &! in
!!$             &                     p_hydro_state(jg)%prog(n_new)    ) ! inout


        END SELECT !(iforcing)
      ENDIF ! lforcing


      ! deallocate local arrays used by tracer transport
      DEALLOCATE( z_mflx_me, z_mflx_ic, z_vn_traj, z_omega_traj,        &
        &         z_delp_mc_now, z_pres_mc_now, z_pres_ic_now, STAT=ist )
      IF (ist /= success) THEN
       CALL finish ( 'mo_hierarchy_management: process_grid',                      &
         &           'deallocation for z_mflx_me, z_mflx_ic, z_vn_traj, '       // &
         &           'z_omega_traj, z_delp_mc_now, z_pres_mc_now, z_pres_ic_now'// &
         &           ' failed' )
      ENDIF

      !==========================================================================
      ! If there are nested domains...
      !==========================================================================
      IF (p_patch(jg)%n_childdom > 0) THEN

        dt_sub         =  dt_loc/2._wp
        rdt_loc        =  1._wp/dt_loc
        mtime_step_sub => newTimedelta(mtime_step)
        mtime_step_sub =  mtime_step_sub * 0.5_wp

        ! Compute time tendencies for interpolation to refined mesh boundaries
        CALL compute_tendencies( ha_dyn_config%ltheta_dyn, p_patch(jg), &
                                 p_hydro_state(jg),n_new,n_sav1,rdt_loc )

        ! Loop over nested domains
        DO jn = 1, p_patch(jg)%n_childdom

          jgc = p_patch(jg)%child_id(jn)

          ! Interpolate tendencies to lateral boundaries of refined mesh (jgc)
          CALL interpolate_tendencies(p_patch,p_hydro_state,p_int_state,p_grf_state,jg,jgc)

          ! In case of gradient-based interpolation, full-field interpolation
          ! is needed for cell-based variables at the beginning of each time step
          ! (i.e. n_sav1 -> n_sav1) to ensure numerical stability
          IF (grf_intmethod_c == 2) THEN

            IF (ha_dyn_config%ltheta_dyn) THEN
              CALL interpol_scal_grf ( p_patch(jg), p_patch(jgc),                        &
                &    p_grf_state(jg)%p_dom(jn), 1, p_hydro_state(jg)%prog(n_sav1)%theta, &
                &    p_hydro_state(jgc)%prog(nnow(jgc))%theta)
            ELSE IF (.NOT. lshallow_water) THEN
              CALL interpol_scal_grf ( p_patch(jg), p_patch(jgc),                       &
                &    p_grf_state(jg)%p_dom(jn), 1, p_hydro_state(jg)%prog(n_sav1)%temp, &
                &    p_hydro_state(jgc)%prog(nnow(jgc))%temp)
            ENDIF

            ALLOCATE (z_pres_sfc(nproma,1,p_patch(jgc)%nblks_c))
            CALL interpol_scal_grf (p_patch(jg), p_patch(jgc), p_grf_state(jg)%p_dom(jn), 1,                  &
              &  RESHAPE(p_hydro_state(jg)%prog(n_sav1)%pres_sfc,(/nproma,1,p_patch(jg)%nblks_c/)), z_pres_sfc)

            p_hydro_state(jgc)%prog(nnow(jgc))%pres_sfc(:,:) = z_pres_sfc(:,1,:)
            DEALLOCATE (z_pres_sfc)
          ENDIF

          IF(ltransport .AND. grf_intmethod_ct == 2) THEN
            CALL interpol_scal_grf (p_patch(jg), p_patch(jgc), p_grf_state(jg)%p_dom(jn), ntracer,                        &
              &  f4din1=p_hydro_state(jg)%prog(n_sav1)%tracer,                                   &
              &  f4dout1=p_hydro_state(jgc)%prog(nnow(jgc))%tracer)
          ENDIF

        ENDDO

        iret = util_cputime(tu, ts)
        t_1 = tu+ts
        t_tot = t_tot + t_1 - t_0

        DO jn = 1, p_patch(jg)%n_childdom

          jgc = p_patch(jg)%child_id(jn)

          IF(p_patch(jgc)%n_patch_cells > 0) THEN
            IF(proc_split) CALL push_glob_comm(p_patch(jgc)%comm, p_patch(jgc)%proc0)

            ! Recursive call to process_grid_level for child grid level
            CALL process_grid( p_patch, p_hydro_state,             &
              &                p_int_state, p_grf_state,           &
              &                ext_data, jgc, nstep_global,        &
              &                l_3tl_init, dt_sub, mtime_step_sub, &
              &                2, grid_datetime)
            IF(proc_split) CALL pop_glob_comm()
          ENDIF

        ENDDO

        iret = util_cputime(tu, ts)
        t_0 = tu+ts

        DO jn = 1, p_patch(jg)%n_childdom

          jgc = p_patch(jg)%child_id(jn)

          ! Call feedback to copy averaged prognostic variables from refined mesh
          ! back to the coarse mesh (i.e. from jgc to jg)
          IF (lfeedback(jgc)) THEN
            CALL feedback(ha_dyn_config%ltheta_dyn, p_patch, p_hydro_state, &
                          p_int_state, p_grf_state, jgc, jg)
          ENDIF

        ENDDO

        CALL deallocateTimedelta(mtime_step_sub)

      ENDIF
      
      !==========================================================================
      ! Interaction with parent and child domains done. Finish this time step
      ! and prepare for the next
      !==========================================================================
      ! Asselin filter for the leapfrog scheme
      !========================================
      IF ( (ha_dyn_config%itime_scheme==leapfrog_expl).OR. &
           (ha_dyn_config%itime_scheme==leapfrog_si) ) THEN

        CALL asselin( ha_dyn_config%asselin_coeff,     &
                      ha_dyn_config%ltheta_dyn,        &
                      p_hydro_state(jg)%prog(n_old),   &
                      p_hydro_state(jg)%prog(n_new),   &
                      p_hydro_state(jg)%prog(n_now) )
      ENDIF

      !====================
      ! Prepare for output
      !====================
      ! LL: this is moved into the mo_ha_stepping along with the output command
!       IF (lprepare_output(jg)) THEN
! 
!         CALL copy_prog_state( p_hydro_state(jg)%prog(n_now),  &! in
!           &                   p_hydro_state(jg)%prog_out,     &! out
!           &                   ha_dyn_config%ltheta_dyn,       &! in
!           &                   ltransport                     ) ! in
! 
!         CALL update_diag_state( p_hydro_state(jg)%prog_out,   &! in
!           &                     p_patch(jg), p_int_state(jg), &! in
!           &                     ext_data(jg),                 &! in
!           &                     p_hydro_state(jg)%diag_out )   ! out
! 
!         CALL update_dyn_output( p_patch(jg), p_int_state(jg), &! in
!           &                     p_hydro_state(jg)%prog_out,   &! in
!           &                     p_hydro_state(jg)%diag_out )   ! inout
! 
!         ! For nested domains, output preparation is not needed during the
!         ! remaining substeps of the same big step
! 
!         lprepare_output(jg) = .FALSE.
! 
!       ENDIF

      !===================
      ! Swap time indices
      !===================
      IF (ltwotime) THEN

        n_temp    = nnow(jg)
        nnow(jg)  = nnew(jg)
        nnew(jg)  = n_temp

      ELSE
        n_temp    = nold(jg)
        nold(jg)  = nnow(jg)
        nnow(jg)  = nnew(jg)
        nnew(jg)  = n_temp

      ENDIF

      !===================================
      ! Check CPU time used by this patch
      !===================================
      iret = util_cputime(tu, ts)
      t_1 = tu+ts
      t_tot = t_tot + t_1 - t_0

      WRITE(message_text,'(''Total time patch '',i2,'': '',f10.3)') jg, t_tot
      CALL message(TRIM(routine),message_text)

      ! update here "grid_datetime" for next time step in this loop
      grid_datetime = grid_datetime + mtime_step

    ENDDO !Loop over all substeps on the current grid level

    CALL deallocateDatetime(grid_datetime)
    
    !---------------
  END SUBROUTINE process_grid


  !-------------------------------------------------------------------------
  !>
  !! Computes the time tendencies of the prognostic variables needed for
  !! interpolation to the lateral boundaries of the nested domains
  !!
  !! @par Revision History
  !! Developed  by Guenther Zaengl, DWD, 2009-06-22
  !!
  SUBROUTINE compute_tendencies (ltheta_dyn, p_patch,p_hydro_state,n_new,n_sav1,rdt)

    LOGICAL,INTENT(IN) :: ltheta_dyn

    TYPE(t_patch),     TARGET, INTENT(IN)    ::  p_patch
    TYPE(t_hydro_atm), TARGET, INTENT(INOUT) ::  p_hydro_state

    ! Time levels from which tendencies are computed
    INTEGER,  INTENT(IN) ::  n_new,n_sav1
    ! Inverse value of time step needed for computing the tendencies
    REAL(wp), INTENT(IN) ::  rdt

    ! local variables

    REAL(wp), DIMENSION(:,:,:), POINTER :: p_temp_new  => NULL()
    REAL(wp), DIMENSION(:,:,:), POINTER :: p_temp_save => NULL()

    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx,       &
      &        jb, jc, je, jk, jt, i_nchdom
    !-----------------------------------------------------------------------

    i_nchdom = MAX(1,p_patch%n_childdom)

    IF (ltheta_dyn) THEN
      p_temp_new  => p_hydro_state%prog(n_new)%theta
      p_temp_save => p_hydro_state%prog(n_sav1)%theta
    ELSE
      p_temp_new  => p_hydro_state%prog(n_new)%temp
      p_temp_save => p_hydro_state%prog(n_sav1)%temp
    ENDIF

!$OMP PARALLEL PRIVATE(i_startblk,i_endblk)

    ! cell-based variables
    i_startblk = p_patch%cells%start_blk(grf_bdywidth_c+1,1)
    i_endblk   = p_patch%cells%end_blk(min_rlcell_int-2,i_nchdom)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk,                          &
        &                i_startidx, i_endidx, grf_bdywidth_c+1, min_rlcell_int-2)

      DO jc = i_startidx, i_endidx
        p_hydro_state%tend_dyn%pres_sfc(jc,jb) =                     &
          &      ( p_hydro_state%prog(n_new )%pres_sfc(jc,jb)        &
          &       -p_hydro_state%prog(n_sav1)%pres_sfc(jc,jb) )  *rdt
      ENDDO

      DO jk = 1, nlev
        DO jc = i_startidx, i_endidx
          p_hydro_state%tend_dyn%temp(jc,jk,jb) =                    &
            &    ( p_temp_new(jc,jk,jb) - p_temp_save(jc,jk,jb) )*rdt
        ENDDO
      ENDDO

    ENDDO
!$OMP END DO

    IF (ltransport) THEN

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk,jt) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk, i_endblk

        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk,                          &
          &                i_startidx, i_endidx, grf_bdywidth_c+1, min_rlcell_int-2)

        DO jt = 1,ntracer
          DO jk = 1, nlev
            DO jc = i_startidx, i_endidx
              p_hydro_state%tend_dyn%tracer(jc,jk,jb,jt) =                &
                &   ( p_hydro_state%prog(n_new )%tracer(jc,jk,jb,jt)      &
                &    -p_hydro_state%prog(n_sav1)%tracer(jc,jk,jb,jt) )*rdt
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!$OMP END DO

    ENDIF

    ! edge-based variables
    i_startblk = p_patch%edges%start_blk(grf_bdywidth_e+1,1)
    i_endblk   = p_patch%edges%end_blk(min_rledge_int-3,i_nchdom)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,je,jk) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk,                          &
        &                i_startidx, i_endidx, grf_bdywidth_e+1, min_rledge_int-3)

      DO jk = 1, nlev
        DO je = i_startidx, i_endidx
          p_hydro_state%tend_dyn%vn(je,jk,jb) =               &
            &  ( p_hydro_state%prog(n_new )%vn(je,jk,jb)      &
            &   -p_hydro_state%prog(n_sav1)%vn(je,jk,jb) )*rdt
        ENDDO
      ENDDO
    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE compute_tendencies



  !>
  !! First time step in the leapfrog scheme.
  !!
  !! The time integration is performed using either a simple Euler forward
  !! scheme (ileapfrog_startup = 1) or a series of sub-steps ( = 2).
  !!
  SUBROUTINE leapfrog_startup( ileapfrog_startup,     &! in
    &                          p_patch, p_int_state,  &! in
    &                          ext_data,              &! in
    &                          jg, dt_loc,            &! in
    &                          p_hydro_state,         &! inout
    &                          n_old, n_now, n_new,   &! inout
    &                          n_sav1               )  ! inout

    ! Arguments

    INTEGER, INTENT(IN) :: ileapfrog_startup
    TYPE(t_patch),TARGET, INTENT(IN)    ::  p_patch(n_dom)
    TYPE(t_int_state),TARGET,INTENT(IN) ::  p_int_state(n_dom)
    TYPE(t_external_data), INTENT(INOUT)::  ext_data(n_dom)

    INTEGER, INTENT(IN) :: jg      !< current grid level
    REAL(wp),INTENT(IN) :: dt_loc  !< time step applicable to local grid level

    TYPE(t_hydro_atm),TARGET,INTENT(INOUT) ::  p_hydro_state(n_dom)

    INTEGER,INTENT(inout)  :: n_old, n_now, n_new, n_sav1

    ! Local variables

    INTEGER  :: n_temp           !< time index
    INTEGER  :: jstep            !< step counter for tracer transport
    INTEGER  :: jb, nlen         !< related to nproma blocking
    REAL(wp) :: zdtime           !< length of time step

    REAL(wp) :: temp_save ( nproma, nlev, p_patch(jg)%nblks_c )

    REAL(wp),POINTER :: p_temp(:,:,:) => NULL()


    SELECT CASE(ileapfrog_startup)
      !==========================================================================
    CASE(1) ! Step 1 is a simple forward step
      !==========================================================================
      ! old,now       new
      !    |-----------|-----------|
      !    0           1           2
      !    * ========> *

      nsav1(jg) = nold(jg)

      n_old  = nold(jg)
      n_now  = nnow(jg)
      n_new  = nnew(jg)
      n_sav1 = nsav1(jg)
      jstep  = 0

      zdtime = 0.5_wp*dt_loc  ! it will be doubled in step_leapfrog_expl

      CALL step_leapfrog_expl( zdtime, 2._wp*zdtime,          &! in
        &                      ha_dyn_config%ltheta_dyn,      &! in
        &                      p_patch(jg), p_int_state(jg),  &! in
        &                      p_hydro_state(jg)%prog(n_old), &! in
        &                      ext_data(jg),                  &! in
        &                      p_hydro_state(jg)%prog(n_now), &! in
        &                      p_hydro_state(jg)%diag,        &! out
        &                      p_hydro_state(jg)%prog(n_new), &! inout
        &                      p_hydro_state(jg)%tend_dyn    ) ! inout

      ! tracer advection
      IF ( ltransport ) THEN

        ! Diagnose some pressure- and velocity-related quantities for
        ! the tracer transport scheme
        CALL prepare_tracer( p_patch(jg), p_int_state(jg),    &! in
          &                  p_hydro_state(jg)%prog(n_now),   &! in
          &                  p_hydro_state(jg)%prog(n_new),   &! in
          &                  ha_dyn_config%itime_scheme,      &! in
          &                  ha_dyn_config%si_2tls,           &! in
          &                  p_hydro_state(jg)%diag,          &! inout
          &                  z_mflx_me, z_vn_traj,            &! out
          &                  z_mflx_ic, z_omega_traj,         &! out
          &                  z_delp_mc_now,                   &! out
          &                  z_pres_mc_now, z_pres_ic_now     )! out



        CALL step_advection( p_patch(jg), p_int_state(jg), zdtime, jstep, & !in
          &                  p_hydro_state(jg)%prog(n_now)%tracer,        & !in
          &                  z_mflx_me, z_vn_traj,                        & !in
          &                  z_mflx_ic,                                   & !in
          &                  z_mflx_ic,                                   & !in
          &                  z_delp_mc_now,                               & !in
          &                  p_hydro_state(jg)%diag%delp_c,               & !in
          &                  z_delp_mc_now,                               & !in
          &                  p_hydro_state(jg)%tend_dyn%tracer,           & !inout
          &                  p_hydro_state(jg)%prog(n_new)%tracer,        & !inout
          &                  p_hydro_state(jg)%diag%hfl_tracer,           & !out
          &                  p_hydro_state(jg)%diag%vfl_tracer )            !out

      END IF

      SELECT CASE(diffusion_config(jg)%hdiff_order)
      CASE(-1)    ! No diffusion
      CASE(24,42) ! Hybrid linear
        CALL hdiff_hyb_lin( jg, p_patch(jg), p_int_state(jg), &! in
                          & p_hydro_state(jg)%diag,           &! in
                          & p_hydro_state(jg)%prog(n_new)     )! inout
      CASE DEFAULT
        CALL hdiff_expl( p_patch(jg), p_hydro_state(jg)%prog(n_new), &
                       & p_hydro_state(jg)%diag, p_int_state(jg),    &
                       & zdtime, jg )
      END SELECT

      !==========================================================================
    CASE(2) ! A scheme using a series of shorter sub-steps.
      !==========================================================================
      ! Step 0.5* is a simple forward step
      !
      !   old,now  new
      !       |-----^-----|-----------|
      !       0    0.5*   1           2
      !       * ==> *

      nsav1(jg) = nold(jg)

      n_old  = nold(jg)
      n_now  = nnow(jg)
      n_new  = nnew(jg)
      n_sav1 = nsav1(jg)
      jstep  = 0

      zdtime = 0.25_wp*dt_loc  ! it will be doubled in step_leapfrog_expl

      CALL step_leapfrog_expl( zdtime, 2._wp*zdtime,          &! in
        &                      ha_dyn_config%ltheta_dyn,      &! in
        &                      p_patch(jg), p_int_state(jg),  &! in
        &                      p_hydro_state(jg)%prog(n_old), &! in
        &                      ext_data(jg),                  &! in
        &                      p_hydro_state(jg)%prog(n_now), &! in
        &                      p_hydro_state(jg)%diag,        &! out
        &                      p_hydro_state(jg)%prog(n_new), &! inout
        &                      p_hydro_state(jg)%tend_dyn    ) ! inout

      ! tracer advection
      IF ( ltransport ) THEN

        ! Diagnose some pressure- and velocity-related quantities for
        ! the tracer transport scheme
        CALL prepare_tracer( p_patch(jg), p_int_state(jg),    &! in
          &                  p_hydro_state(jg)%prog(n_now),   &! in
          &                  p_hydro_state(jg)%prog(n_new),   &! in
          &                  ha_dyn_config%itime_scheme,      &! in
          &                  ha_dyn_config%si_2tls,           &! in
          &                  p_hydro_state(jg)%diag,          &! inout
          &                  z_mflx_me, z_vn_traj,            &! out
          &                  z_mflx_ic, z_omega_traj,         &! out
          &                  z_delp_mc_now,                   &! out
          &                  z_pres_mc_now, z_pres_ic_now     )! out


        CALL step_advection( p_patch(jg), p_int_state(jg), zdtime, jstep,& !in
          &                 p_hydro_state(jg)%prog(n_now)%tracer,        & !in
          &                 z_mflx_me, z_vn_traj,                        & !in
          &                 z_mflx_ic,                                   & !in
          &                 z_mflx_ic,                                   & !in
          &                 z_delp_mc_now,                               & !in
          &                 p_hydro_state(jg)%diag%delp_c,               & !in
          &                 z_delp_mc_now,                               & !in
          &                 p_hydro_state(jg)%tend_dyn%tracer,           & !inout
          &                 p_hydro_state(jg)%prog(n_new)%tracer,        & !inout
          &                 p_hydro_state(jg)%diag%hfl_tracer,           & !out
          &                 p_hydro_state(jg)%diag%vfl_tracer )            !out

      END IF

      SELECT CASE(diffusion_config(jg)%hdiff_order)
      CASE(-1)    ! No diffusion
      CASE(24,42) ! Hybrid linear
        CALL hdiff_hyb_lin( jg, p_patch(jg), p_int_state(jg), &! in
                          & p_hydro_state(jg)%diag,           &! in
                          & p_hydro_state(jg)%prog(n_new)     )! inout
      CASE DEFAULT
        CALL hdiff_expl( p_patch(jg), p_hydro_state(jg)%prog(n_new), &
                       & p_hydro_state(jg)%diag, p_int_state(jg),    &
                       & zdtime, jg )
      END SELECT


      ! Step 0.5 is a backward step
      !
      !nold  = nold           ! step 0        old  now,new
      n_temp     = nnow(jg)  !                 |-----^-----|-----------|
      nnow(jg)  = nnew(jg)   ! step 0.5*       0    0.5    1           2
      nnew(jg)  = n_temp     ! step 0.5        * ==> *

      n_now = nnow(jg)
      n_new = nnew(jg)

      zdtime = 0.25_wp*dt_loc  ! it will be doubled in step_leapfrog_expl

      CALL step_leapfrog_expl( zdtime, 0._wp,                     & !in
        &                      ha_dyn_config%ltheta_dyn,          &! in
        &                      p_patch(jg), p_int_state(jg),      & !in
        &                      p_hydro_state(jg)%prog(n_old),     & !in
        &                      ext_data(jg),                      & !in
        &                      p_hydro_state(jg)%prog(n_now),     & !in
        &                      p_hydro_state(jg)%diag,            & !out
        &                      p_hydro_state(jg)%prog(n_new),     & !inout
        &                      p_hydro_state(jg)%tend_dyn       )   !inout

      ! tracer advection
      IF ( ltransport ) THEN

        ! Diagnose some pressure- and velocity-related quantities for
        ! the tracer transport scheme
        CALL prepare_tracer( p_patch(jg), p_int_state(jg),    &! in
          &                  p_hydro_state(jg)%prog(n_now),   &! in
          &                  p_hydro_state(jg)%prog(n_new),   &! in
          &                  ha_dyn_config%itime_scheme,      &! in
          &                  ha_dyn_config%si_2tls,           &! in
          &                  p_hydro_state(jg)%diag,          &! inout
          &                  z_mflx_me, z_vn_traj,            &! out
          &                  z_mflx_ic, z_omega_traj,         &! out
          &                  z_delp_mc_now,                   &! out
          &                  z_pres_mc_now, z_pres_ic_now     )! out


        CALL step_advection( p_patch(jg), p_int_state(jg), zdtime, jstep,& !in
          &                 p_hydro_state(jg)%prog(n_now)%tracer,        & !in
          &                 z_mflx_me, z_vn_traj,                        & !in
          &                 z_mflx_ic,                                   & !in
          &                 z_mflx_ic,                                   & !in
          &                 z_delp_mc_now,                               & !in
          &                 p_hydro_state(jg)%diag%delp_c,               & !in
          &                 z_delp_mc_now,                               & !in
          &                 p_hydro_state(jg)%tend_dyn%tracer,           & !inout
          &                 p_hydro_state(jg)%prog(n_new)%tracer,        & !inout
          &                 p_hydro_state(jg)%diag%hfl_tracer,           & !out
          &                 p_hydro_state(jg)%diag%vfl_tracer )            !out

      END IF

      SELECT CASE(diffusion_config(jg)%hdiff_order)
      CASE(-1)    ! No diffusion
      CASE(24,42) ! Hybrid linear
        CALL hdiff_hyb_lin( jg, p_patch(jg), p_int_state(jg), &! in
                          & p_hydro_state(jg)%diag,           &! in
                          & p_hydro_state(jg)%prog(n_new)     )! inout
      CASE DEFAULT
        CALL hdiff_expl( p_patch(jg), p_hydro_state(jg)%prog(n_new), &
                       & p_hydro_state(jg)%diag, p_int_state(jg),    &
                       & zdtime, jg )
      END SELECT

      ! Step 1 is a leapfrog step with half time-interval
      !
      !nold  = nold          !            old   now   new
      n_temp     = nnow(jg)  !             |-----^-----|-----------|
      nnow(jg)  = nnew(jg)   ! step 0.5    0    0.5    1           2
      nnew(jg)  = n_temp     ! step 1      * ========> *

      n_now = nnow(jg)
      n_new = nnew(jg)
      jstep = 1

      zdtime = 0.5_wp*dt_loc

      CALL step_leapfrog_expl( zdtime, zdtime,                 & !in
        &                      ha_dyn_config%ltheta_dyn,       &! in
        &                      p_patch(jg), p_int_state(jg),   & !in
        &                      p_hydro_state(jg)%prog(n_old),  & !in
        &                      ext_data(jg),                   & !in
        &                      p_hydro_state(jg)%prog(n_now),  & !in
        &                      p_hydro_state(jg)%diag,         & !out
        &                      p_hydro_state(jg)%prog(n_new),  & !inout
        &                      p_hydro_state(jg)%tend_dyn     )  !inout

      IF (ha_dyn_config%itime_scheme==leapfrog_si) THEN

        IF (ha_dyn_config%ltheta_dyn) THEN
          CALL convert_theta2t_lin( p_patch(jg), p_hydro_state(jg)%prog(n_now),&
            &                       p_hydro_state(jg)%prog(n_new),             &
            &                       p_hydro_state(jg)%diag )

          p_temp => p_hydro_state(jg)%prog(n_new)%temp

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,nlen) ICON_OMP_DEFAULT_SCHEDULE
          DO jb = 1,p_patch(jg)%nblks_c

            IF (jb /= p_patch(jg)%nblks_c) THEN
              nlen = nproma
            ELSE
              nlen = p_patch(jg)%npromz_c
            ENDIF

            temp_save(1:nlen,1:nlev,jb) = p_temp(1:nlen,1:nlev,jb)

          ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
        ENDIF

        CALL si_correction( ha_dyn_config%lsi_3d,                 &! in
          &                 ha_dyn_config%si_coeff,               &! in 
          &                 ha_dyn_config%si_rtol,                &! in
          &                 lshallow_water,                       &! in
          &                 zdtime, p_patch(jg), p_int_state(jg), &! in
          &                 p_hydro_state(jg)%prog(n_old),        &! in
          &                 p_hydro_state(jg)%prog(n_now),        &! in
          &                 p_hydro_state(jg)%prog(n_new))         ! inout

        IF (ha_dyn_config%ltheta_dyn) THEN
          CALL convert_t2theta_lin( p_patch(jg),p_hydro_state(jg)%prog(n_new),&
            &                       p_hydro_state(jg)%diag, temp_save )
        ENDIF

      ENDIF ! with semi-implicit correction

      ! tracer advection
      IF ( ltransport ) THEN

        ! Diagnose some pressure- and velocity-related quantities for
        ! the tracer transport scheme
        CALL prepare_tracer( p_patch(jg), p_int_state(jg),  &! in
          &           p_hydro_state(jg)%prog(n_now),        &! in
          &           p_hydro_state(jg)%prog(n_new),        &! in
          &           ha_dyn_config%itime_scheme,           &! in
          &           ha_dyn_config%si_2tls,                &! in
          &           p_hydro_state(jg)%diag,               &! inout
          &           z_mflx_me, z_vn_traj,                 &! out
          &           z_mflx_ic, z_omega_traj,              &! out
          &           z_delp_mc_now,                        &! out
          &           z_pres_mc_now, z_pres_ic_now          )! out


        CALL step_advection( p_patch(jg), p_int_state(jg), zdtime, jstep,& !in
          &                 p_hydro_state(jg)%prog(n_now)%tracer,        & !in
          &                 z_mflx_me, z_vn_traj,                        & !in
          &                 z_mflx_ic,                                   & !in
          &                 z_mflx_ic,                                   & !in
          &                 z_delp_mc_now,                               & !in
          &                 p_hydro_state(jg)%diag%delp_c,               & !in
          &                 z_delp_mc_now,                               & !in
          &                 p_hydro_state(jg)%tend_dyn%tracer,           & !inout
          &                 p_hydro_state(jg)%prog(n_new)%tracer,        & !inout
          &                 p_hydro_state(jg)%diag%hfl_tracer,           & !out
          &                 p_hydro_state(jg)%diag%vfl_tracer )            !out

      END IF

      SELECT CASE(diffusion_config(jg)%hdiff_order)
      CASE(-1)    ! No diffusion
      CASE(24,42) ! Hybrid linear
        CALL hdiff_hyb_lin( jg, p_patch(jg), p_int_state(jg), &! in
                          & p_hydro_state(jg)%diag,           &! in
                          & p_hydro_state(jg)%prog(n_new)     )! inout
      CASE DEFAULT
        CALL hdiff_expl( p_patch(jg), p_hydro_state(jg)%prog(n_new), &
                       & p_hydro_state(jg)%diag, p_int_state(jg),    &
                       & zdtime, jg )
      END SELECT

      ! Asselin filter is now called after feedback

    CASE DEFAULT
      CALL finish('mo_hierarchy_management', &
                  'unknown choice of ha_dyn_config%ileapfrog_startup')
    END SELECT

  END SUBROUTINE leapfrog_startup


END MODULE mo_hierarchy_management
