#ifdef __xlC__
!@PROCESS NOHOT
!@PROCESS NOOPTimize
#endif
!>
!! Initializes and controls the time stepping in the nonhydrostatic model.
!!
!!
!! @par Revision History
!! Initial release by Almut Gassmann, MPI-M (2009-02-06)
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
MODULE mo_nh_stepping
!-------------------------------------------------------------------------
!
!    ProTeX FORTRAN source: Style 2
!    modified for ICON project, DWD/MPI-M 2006
!
!-------------------------------------------------------------------------
!
!

  USE mo_kind,                ONLY: wp
  USE mo_nonhydro_state,      ONLY: t_nh_state, t_nh_prog, t_nh_diag, t_nh_metrics, &
                                    construct_nh_state, bufr
  USE mo_nonhydrostatic_config,ONLY: iadv_rcf, l_nest_rcf, itime_scheme
  USE mo_diffusion_config,    ONLY: diffusion_config
  USE mo_dynamics_config,     ONLY: nnow,nnew, nnow_rcf, nnew_rcf, nsav1,nsav2
  USE mo_io_config,           ONLY: l_outputtime, l_diagtime, l_checkpoint_time, &
    &                               lwrite_pzlev
  USE mo_parallel_config,     ONLY: nproma, itype_comm
  USE mo_run_config,          ONLY: ltestcase, dtime, nsteps,  &
    &                               ltransport, ntracer, lforcing, iforcing, &
    &                               msg_level, ltimer
  USE mo_grid_config,         ONLY: global_cell_type
    
  USE mo_atm_phy_nwp_config,  ONLY: tcall_phy
  USE mo_nwp_phy_init,        ONLY: init_nwp_phy
  USE mo_nwp_phy_state,       ONLY: prm_diag, prm_nwp_tend, mean_charlen
  USE mo_nwp_lnd_state,       ONLY: p_lnd_state
  USE mo_ext_data,            ONLY: ext_data
  USE mo_model_domain,        ONLY: t_patch
  USE mo_model_domain_import, ONLY: n_dom, lfeedback, l_limited_area, &
    &                               n_dom_start, lredgrid_phys
  USE mo_nh_testcases,        ONLY: init_nh_testtopo, init_nh_testcase, nh_test_name, &
    &                               rotate_axis_deg
  USE mo_nh_pa_test,          ONLY: set_nh_w_rho
  USE mo_nh_df_test,          ONLY: get_nh_df_velocity, get_nh_df_mflx_rho
  USE mo_nh_hex_util,         ONLY: forcing_straka, momentum_adv
  USE mo_nh_supervise,        ONLY: supervise_total_integrals_nh
  USE mo_interpolation,       ONLY: t_int_state, rbf_vec_interpol_cell, &
                                    edges2cells_scalar,&
                                    verts2edges_scalar,&
                                    edges2verts_scalar
  USE mo_grf_interpolation,   ONLY: t_gridref_state
  USE mo_grf_bdyintp,         ONLY: interpol_scal_grf
  USE mo_nh_nest_utilities,   ONLY: compute_tendencies, boundary_interpolation, &
                                    complete_nesting_setup, prep_bdy_nudging,   &
                                    outer_boundary_nudging, nest_boundary_nudging
  USE mo_nh_feedback,         ONLY: feedback
  USE mo_datetime,            ONLY: t_datetime, print_datetime, add_time
  USE mo_timer,               ONLY: timer_total, timer_start, timer_stop
  USE mo_output,              ONLY: init_output_files, write_output,  &
    &                               create_restart_file
  USE mo_io_restart,          ONLY: write_restart_info_file
  USE mo_exception,           ONLY: message, message_text, finish
  USE mo_impl_constants,      ONLY: SUCCESS, MAX_CHAR_LENGTH,iphysproc,itconv,   &
    &                               itccov, itrad, itradheat, itsso, itsatad,    &
    &                               inwp
  USE mo_divergent_modes,     ONLY: divergent_modes_5band
  USE mo_math_operators,      ONLY: div_avg, div
  USE mo_solve_nonhydro,      ONLY: solve_nh
  USE mo_solve_nh_async,      ONLY: solve_nh_ahc
  USE mo_advection_stepping,  ONLY: step_advection
  USE mo_nh_dtp_interface,    ONLY: prepare_tracer
  USE mo_nh_vert_interp,      ONLY: intp_to_p_and_z_levels
  USE mo_nh_diffusion,        ONLY: diffusion_tria, diffusion_hex
  USE mo_mpi,                 ONLY: my_process_is_stdio, my_process_is_mpi_parallel
#ifdef NOMPI
  USE mo_mpi,                 ONLY: my_process_is_mpi_all_seq
#endif
  
  USE mo_sync,                ONLY: sync_patch_array_mult, &
                                    push_glob_comm, pop_glob_comm, global_max, &
                                    SYNC_C, SYNC_E, sync_patch_array
  USE mo_communication,       ONLY: start_delayed_exchange, do_delayed_exchange
  USE mo_subdivision,         ONLY: proc_split
  USE mo_nh_interface_nwp,    ONLY: nwp_nh_interface
  USE mo_phys_nest_utilities, ONLY: interpol_phys_grf, feedback_phys_diag
  USE mo_vertical_grid,       ONLY: set_nh_metrics
  USE mo_nh_diagnose_pres_temp,ONLY: diagnose_pres_temp
  USE mo_nh_held_suarez_interface, ONLY: held_suarez_nh_interface
  USE mo_vertical_coord_table,ONLY: vct
  USE mo_master_control,      ONLY: is_restart_run
  USE mo_io_restart_attributes,ONLY: get_restart_attribute

  USE mo_nwp_mpiomp_rrtm_interface, ONLY: nwp_start_radiation_ompthread, model_end_ompthread, &
    & init_ompthread_radiation
  USE mo_parallel_config,     ONLY: parallel_radiation_omp, nh_stepping_ompthreads

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

  ! for preparation of transport with optional reduced calling frequency
  TYPE :: t_prepare_adv
    ! mass flux at full level edges (currently at N+1\2)
    REAL(wp), ALLOCATABLE :: mass_flx_me(:,:,:)

    !mass flux at half level centers (currently at N+1\2)
    REAL(wp), ALLOCATABLE :: mass_flx_ic(:,:,:)

    ! horizontal velocity at edges for computation of backward trajectories
    ! (currently at N+1/2)
    REAL(wp), ALLOCATABLE :: vn_traj(:,:,:)

    ! vertical velocity at half level centers for computation of
    ! backward trajectories (currently at N+1/2)
    REAL(wp), ALLOCATABLE :: w_traj(:,:,:)

    ! density times layer thickness at cell center at time step N
    REAL(wp), ALLOCATABLE :: rhodz_mc_now(:,:,:)

    !< density times layer thickness at cell center at time step N+1
    REAL(wp), ALLOCATABLE :: rhodz_mc_new(:,:,:)

    !< density at half levels (currently at N+1/2)
    REAL(wp), ALLOCATABLE :: rho_ic(:,:,:)

    !< vertical tracer flux at domain top (time average; n+1/2)
    REAL(wp), ALLOCATABLE :: topflx_tra(:,:,:)

  END TYPE t_prepare_adv


  ! counter for 'reduced calling frequency' and Marchuk splitting
  TYPE :: t_step_adv
    ! Counts total number of dynamics time steps for each patch (necessary for
    ! generalization of rcf to arbitrary (even) number of iadv_rcf)
    INTEGER :: ntsteps

    ! Determines sequence of operations for Marchuk-splitting (for transport)
    INTEGER :: marchuk_order
  END TYPE t_step_adv


  TYPE(t_prepare_adv), ALLOCATABLE :: prep_adv(:)  ! n_dom

  TYPE(t_step_adv),    ALLOCATABLE :: jstep_adv(:) ! n_dom


  ! additional flow control variables that need to be dimensioned with the
  ! number of model domains
  LOGICAL, ALLOCATABLE :: lstep_adv(:)   ! determines whether tracer continuity equations
                                         ! should be integrated (.true.) or not (.false.)

  LOGICAL, ALLOCATABLE :: lcall_phy(:,:) ! contains information which physics package
                                         ! must be called at the current timestep
                                         ! and on the current domain.

  REAL(wp), ALLOCATABLE :: t_elapsed_phy(:,:)  ! time (in s) since the last call of
                                               ! the corresponding physics package

  LOGICAL, ALLOCATABLE :: linit_slowphy(:) ! determines whether slow physics has already been initialized

  LOGICAL, ALLOCATABLE :: linit_dyn(:)  ! determines whether dynamics has already been initialized


  REAL(wp), ALLOCATABLE :: sim_time(:)  ! elapsed simulation time


  PUBLIC :: prepare_nh_integration, perform_nh_stepping

  CONTAINS

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!
!
  !>
  !! Initialisation of the nonhydrostatic state and initial conditions.
  !!
  !!
  !! @par Revision History
  !! Initial release by Almut Gassmann, (2009-03-06)
  !!
  SUBROUTINE prepare_nh_integration (p_patch, p_nh_state, p_int, p_grf)
!
  TYPE(t_patch), TARGET, INTENT(INOUT) :: p_patch(n_dom)
  TYPE(t_int_state), INTENT(IN)        :: p_int(n_dom)
  TYPE(t_gridref_state), INTENT(INOUT) :: p_grf(n_dom)

  TYPE(t_nh_state), TARGET, INTENT(INOUT):: p_nh_state(n_dom)

  INTEGER :: ntl

!-----------------------------------------------------------------------

  ! for the split explict scheme, ntl is always 2
  ntl = 2
  CALL construct_nh_state(p_patch, p_nh_state, ntl)

  IF (ltestcase) THEN
    CALL init_nh_testtopo(p_patch)    ! set analytic topography
  ENDIF

  CALL set_nh_metrics(p_patch, p_nh_state, p_int)

  IF (n_dom > 1) THEN
    CALL complete_nesting_setup(p_patch, p_nh_state, p_grf)
  ENDIF

  IF (ltestcase) THEN
    CALL init_nh_testcase(p_patch, p_nh_state, p_int, ntl)
  ENDIF

  END SUBROUTINE prepare_nh_integration


  !-------------------------------------------------------------------------
  !>
  !! Organizes nonhydrostatic time stepping
  !! Currently we assume to have only one grid level.
  !!
  !! @par Revision History
  !! Initial release by Almut Gassmann, (2009-04-15)
  !!
  SUBROUTINE perform_nh_stepping (p_patch, p_int_state, p_grf_state, p_nh_state, &
                                  datetime, n_io, n_file, n_checkpoint, n_diag,  &
                                  l_have_output )
!
  TYPE(t_patch), TARGET, INTENT(IN)            :: p_patch(n_dom_start:n_dom)
  TYPE(t_int_state), TARGET, INTENT(IN)        :: p_int_state(n_dom_start:n_dom)
  TYPE(t_gridref_state), TARGET, INTENT(INOUT) :: p_grf_state(n_dom_start:n_dom)
  INTEGER, INTENT(IN)                          :: n_io, n_file, n_checkpoint, n_diag
  LOGICAL, INTENT(INOUT) :: l_have_output

  TYPE(t_nh_state), TARGET, INTENT(INOUT):: p_nh_state(n_dom)
  TYPE(t_datetime), INTENT(INOUT)      :: datetime

  INTEGER                              :: jg

!$  INTEGER omp_get_num_threads
!$  INTEGER omp_get_max_threads
!$  INTEGER omp_get_max_active_levels
!-----------------------------------------------------------------------

  CALL allocate_nh_stepping (p_patch)

  IF (iforcing == inwp .AND. is_restart_run()) THEN
    DO jg=1, n_dom
      CALL init_nwp_phy( dtime                 ,&
           & p_patch(jg)                       ,&
           & p_nh_state(jg)%metrics            ,&
           & p_nh_state(jg)%prog(nnow(jg))     ,&
           & p_nh_state(jg)%prog(nnew(jg))     ,&
           & p_nh_state(jg)%diag               ,&
           & prm_diag(jg)                      ,&
           & prm_nwp_tend(jg)                  ,&
           & p_lnd_state(jg)%prog_lnd(nnow(jg)),&
           & p_lnd_state(jg)%prog_lnd(nnew(jg)),&
           & p_lnd_state(jg)%diag_lnd          ,&
           & ext_data(jg)                      ,&
           & mean_charlen(jg)                   )
    ENDDO
  ELSE IF (iforcing == inwp) THEN ! for cold start, use atmospheric fields at time level nnow only
    DO jg=1, n_dom
      CALL init_nwp_phy( dtime                 ,&
           & p_patch(jg)                       ,&
           & p_nh_state(jg)%metrics            ,&
           & p_nh_state(jg)%prog(nnow(jg))     ,&
           & p_nh_state(jg)%prog(nnow(jg))     ,&
           & p_nh_state(jg)%diag               ,&
           & prm_diag(jg)                      ,&
           & prm_nwp_tend(jg)                  ,&
           & p_lnd_state(jg)%prog_lnd(nnow(jg)),&
           & p_lnd_state(jg)%prog_lnd(nnew(jg)),&
           & p_lnd_state(jg)%diag_lnd          ,&
           & ext_data(jg)                      ,&
           & mean_charlen(jg)                   )
    ENDDO
  ENDIF

  IF (parallel_radiation_omp) THEN

    !---------------------------------------
    CALL init_ompthread_radiation()
    
!$    CALL omp_set_nested(.true.)
!$    CALL omp_set_num_threads(2)
!$    write(0,*) 'omp_get_max_active_levels=',omp_get_max_active_levels
!$    write(0,*) 'omp_get_max_threads=',omp_get_max_threads()
!$OMP PARALLEL SECTIONS
!$OMP SECTION
!$  CALL omp_set_num_threads(nh_stepping_ompthreads)
!$    write(0,*) 'This is the nh_timeloop, max threads=',omp_get_max_threads()
!$    write(0,*) 'omp_get_num_threads=',omp_get_num_threads()

    CALL perform_nh_timeloop (p_patch, p_int_state, p_grf_state, p_nh_state, &
                              datetime, n_io, n_file, n_checkpoint, n_diag,  &
                              l_have_output )
    CALL model_end_ompthread()

!$OMP SECTION
!$  write(0,*) 'This is the nwp_parallel_radiation_thread, max threads=',&
!$    omp_get_max_threads()
  CALL nwp_start_radiation_ompthread()
!$OMP END PARALLEL SECTIONS

  ELSE
    !---------------------------------------

    CALL perform_nh_timeloop (p_patch, p_int_state, p_grf_state, p_nh_state, &
                              datetime, n_io, n_file, n_checkpoint, n_diag,  &
                              l_have_output )
  ENDIF


  CALL deallocate_nh_stepping ()


  END SUBROUTINE perform_nh_stepping
  !-------------------------------------------------------------------------

 !-------------------------------------------------------------------------
  !>
  !! Organizes nonhydrostatic time stepping
  !! Currently we assume to have only one grid level.
  !!
  !! @par Revision History
  !! Initial release by Almut Gassmann, (2009-04-15)
  !!
  SUBROUTINE perform_nh_timeloop (p_patch, p_int_state, p_grf_state, p_nh_state, &
                                  datetime, n_io, n_file, n_checkpoint, n_diag,  &
                                  l_have_output )
!
    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = 'mo_nh_stepping:perform_nh_timeloop'

  TYPE(t_patch), TARGET, INTENT(IN)            :: p_patch(n_dom_start:n_dom)
  TYPE(t_int_state), TARGET, INTENT(IN)        :: p_int_state(n_dom_start:n_dom)
  TYPE(t_gridref_state), TARGET, INTENT(INOUT) :: p_grf_state(n_dom_start:n_dom)
  INTEGER, INTENT(IN)                          :: n_io, n_file, n_checkpoint, n_diag
  LOGICAL, INTENT(INOUT) :: l_have_output

  TYPE(t_nh_state), TARGET, INTENT(INOUT):: p_nh_state(n_dom)
  TYPE(t_datetime), INTENT(INOUT)      :: datetime

  INTEGER                              :: jfile, jstep, jb, nlen, jg
  REAL(wp)                             :: vnmax, wmax
  REAL(wp)                             :: vn_aux(p_patch(1)%nblks_int_e)
  REAL(wp)                             :: w_aux(p_patch(1)%nblks_int_c)
  REAL(wp), DIMENSION(:,:,:), POINTER  :: p_vn, p_w

!$  INTEGER omp_get_num_threads
!-----------------------------------------------------------------------

  jfile = 1

  IF (ltimer) CALL timer_start(timer_total)

  TIME_LOOP: DO jstep = 1, nsteps

    CALL add_time(dtime,0,0,0,datetime)

    WRITE(message_text,'(a,i10)') 'TIME STEP n: ', jstep
    CALL message(TRIM(routine),message_text)

    IF (msg_level >= 5) THEN ! print maximum velocities in global domain

      p_vn => p_nh_state(1)%prog(nnow(1))%vn
      p_w  => p_nh_state(1)%prog(nnow(1))%w

!$OMP PARALLEL
!$OMP DO PRIVATE(jb, nlen)
      DO jb = 1, p_patch(1)%nblks_int_e
        IF (jb /= p_patch(1)%nblks_int_e) THEN
          nlen = nproma
        ELSE
          nlen = p_patch(1)%npromz_int_e
        ENDIF
        vn_aux(jb) = MAXVAL(ABS(p_vn(1:nlen,:,jb)))
      ENDDO
!$OMP END DO
!$OMP DO PRIVATE(jb, nlen)
      DO jb = 1, p_patch(1)%nblks_int_c
        IF (jb /= p_patch(1)%nblks_int_c) THEN
          nlen = nproma
        ELSE
          nlen = p_patch(1)%npromz_int_c
        ENDIF
        w_aux(jb) = MAXVAL(ABS(p_w(1:nlen,:,jb)))
      ENDDO
!$OMP END DO
!$OMP END PARALLEL

      vnmax = MAXVAL(vn_aux)
      wmax  = MAXVAL(w_aux)

      vnmax = global_max(vnmax) ! Get max over all PEs
      wmax  = global_max(wmax) ! Get max over all PEs

      WRITE(message_text,'(a,2e14.6)') 'MAXABS VN, W ', vnmax, wmax
      CALL message(TRIM(routine),message_text)

    ENDIF ! msg_level >= 5

    ! Store first old exner pressure
    ! (to prepare some kind of divergence damping, or to account for
    ! physically based 'implicit weights' in forward backward time stepping)
    IF (jstep == 1 .AND. .NOT. is_restart_run()) THEN
!$OMP PARALLEL PRIVATE(jg)
!   write(0,*) 'Entering perform_nh_timeloop, threads=',omp_get_num_threads()
      DO jg = 1, n_dom
!$OMP WORKSHARE
        p_nh_state(jg)%diag%exner_old(:,:,:)=&
        & p_nh_state(jg)%prog(nnow(1))%exner(:,:,:)
!$OMP END WORKSHARE
      ENDDO
!$OMP END PARALLEL
    ENDIF


    !--------------------------------------------------------------------------
    ! Set output flags
    !--------------------------------------------------------------------------
    IF (MOD(jstep,n_io) == 0 .OR. jstep==nsteps) THEN
      l_outputtime = .TRUE. ! Output is written at the end of the time step,
    ELSE                    ! thus diagnostic quantities need to be computed
      l_outputtime = .FALSE.
    ENDIF

    IF (jstep == 1 .OR. MOD(jstep,n_diag) == 0 .OR. jstep==nsteps) THEN
      l_diagtime = .TRUE. ! Diagnostic output is written at the end of the time step,
                          ! thus diagnostic quantities need to be computed
    ENDIF

    IF ( MOD(jstep,n_checkpoint)==0 ) THEN
      l_checkpoint_time = .TRUE.
    ELSE
      l_checkpoint_time = .FALSE.
    ENDIF

    !
    ! dynamics stepping
    !
    CALL integrate_nh(p_nh_state, p_patch, p_int_state, p_grf_state, &
                      1, jstep, dtime, sim_time, 1)

    ! output of results
    ! note: nnew has been replaced by nnow here because the update
    IF (l_outputtime) THEN

      ! Interpolate selected fields to p- and/or z-levels
      IF (lwrite_pzlev) THEN
        CALL intp_to_p_and_z_levels(p_patch(1:), prm_diag, p_nh_state)
      ENDIF

      CALL write_output( datetime, sim_time(1) )
      CALL message('','Output at:')
      CALL print_datetime(datetime)
      l_have_output = .TRUE.

    ENDIF



    ! Diagnostics computation is not yet properly MPI-parallelized
#ifdef NOMPI
    IF(global_cell_type == 3) THEN
      IF (l_diagtime .AND. my_process_is_mpi_all_seq() .AND. &
        & (lstep_adv(1) .OR. jstep==nsteps))  THEN
        IF (jstep == iadv_rcf) THEN
          CALL supervise_total_integrals_nh( 1, p_patch(1:), p_nh_state,      &
                                           & nnow(1:n_dom), nnow_rcf(1:n_dom))
        ELSE
          CALL supervise_total_integrals_nh(jstep, p_patch(1:), p_nh_state,   &
                                           & nnow(1:n_dom), nnow_rcf(1:n_dom))
        ENDIF
        l_diagtime = .FALSE.
      ENDIF
    ENDIF
#endif
    IF(global_cell_type == 6 .AND. l_diagtime) THEN
      CALL supervise_total_integrals_nh(jstep, p_patch(1:), p_nh_state,   &
                                       & nnow(1:n_dom),  nnow_rcf(1:n_dom))
    ENDIF


    ! close the current output file and trigger a new one
    IF (MOD(jstep,n_file) == 0 .and. jstep/=nsteps) THEN

      jfile = jfile +1
      call init_output_files(jfile,lclose=l_have_output)

    ENDIF


    !--------------------------------------------------------------------------
    ! Write restart file
    !--------------------------------------------------------------------------
    IF (l_checkpoint_time) THEN
      DO jg = 1, n_dom
        CALL create_restart_file( p_patch(jg), datetime,        &
                                & vct, jfile, l_have_output,    &
                                & t_elapsed_phy,                &
                                & lcall_phy, sim_time(jg),      &
                                & jstep_adv(jg)%ntsteps,        &
                                & jstep_adv(jg)%marchuk_order   )
      END DO

      ! Create the master (meta) file in ASCII format which contains
      ! info about which files should be read in for a restart run.
      CALL write_restart_info_file
    END IF


  ENDDO TIME_LOOP

  IF (ltimer) CALL timer_stop(timer_total)


  END SUBROUTINE perform_nh_timeloop
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! @par Revision History
  !!
  SUBROUTINE deallocate_nh_stepping ()
!

  INTEGER                              ::  jg, ist

!-----------------------------------------------------------------------
  !
  ! deallocate auxiliary fields for tracer transport and rcf
  !
  DO jg = 1, n_dom
    DEALLOCATE( prep_adv(jg)%mass_flx_me, prep_adv(jg)%mass_flx_ic,    &
      &         prep_adv(jg)%vn_traj, prep_adv(jg)%w_traj,             &
      &         prep_adv(jg)%rhodz_mc_now, prep_adv(jg)%rhodz_mc_new,  &
      &         prep_adv(jg)%rho_ic, prep_adv(jg)%topflx_tra, STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ( 'mo_nh_stepping: perform_nh_stepping',            &
        &    'deallocation for mass_flx_me, mass_flx_ic, vn_traj,' // &
        &    'w_traj, rhodz_mc_now, rhodz_mc_new, rho_ic, '        // &
        &    'topflx_tra failed' )
    ENDIF
    DEALLOCATE(bufr(jg)%send_c1, &
      &      bufr(jg)%recv_c1,bufr(jg)%send_c3,bufr(jg)%recv_c3, &
      &      bufr(jg)%send_e1,bufr(jg)%recv_e1,bufr(jg)%send_e2, &
      &      bufr(jg)%recv_e2,bufr(jg)%send_e3,bufr(jg)%recv_e3, &
      &      bufr(jg)%send_v2,bufr(jg)%recv_v2, STAT=ist )
    IF (ist /= SUCCESS) &
      CALL finish ( 'mo_nh_stepping: perform_nh_stepping',  &
      &    'deallocation of MPI exchange buffers failed' )
  ENDDO

  DEALLOCATE( bufr, STAT=ist )
  IF (ist /= SUCCESS) &
    CALL finish ( 'mo_nh_stepping: perform_nh_stepping',            &
      &    'deallocation of bufr failed' )

  DEALLOCATE( prep_adv, STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ( 'mo_nh_stepping: perform_nh_stepping',              &
      &    'deallocation for prep_adv failed' )
  ENDIF

  DEALLOCATE( jstep_adv, STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ( 'mo_nh_stepping: perform_nh_stepping',              &
      &    'deallocation for jstep_adv failed' )
  ENDIF

  !
  ! deallocate flow control variables
  !
  DEALLOCATE( lstep_adv, lcall_phy, linit_slowphy, linit_dyn, t_elapsed_phy, &
    &         sim_time, STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ( 'mo_nh_stepping: perform_nh_stepping',          &
      &    'deallocation for lstep_adv, lcall_phy,' //            &
      &    't_elapsed_phy failed' )
  ENDIF

  END SUBROUTINE deallocate_nh_stepping
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! @par Revision History
  !!
  SUBROUTINE allocate_nh_stepping (p_patch)
!
  TYPE(t_patch), TARGET, INTENT(IN)            :: p_patch(n_dom_start:n_dom)

  INTEGER                              :: jg, jp !, nlen
  INTEGER                              :: ist
  CHARACTER(len=MAX_CHAR_LENGTH)       :: attname   ! attribute name

!-----------------------------------------------------------------------
  ! Allocate global buffers for MPI communication
  ALLOCATE(bufr(n_dom), STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ( 'mo_nh_stepping: perform_nh_stepping',           &
    &      'allocation of buffer for MPI communication failed' )
  ENDIF
  !
  ! allocate axiliary fields for transport
  !
  ALLOCATE(prep_adv(n_dom), STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ( 'mo_nh_stepping: perform_nh_stepping',           &
    &      'allocation for prep_adv failed' )
  ENDIF

  ALLOCATE(jstep_adv(n_dom), STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ( 'mo_nh_stepping: perform_nh_stepping',           &
    &      'allocation for jstep_adv failed' )
  ENDIF


  ! allocate flow control variables for transport and slow physics calls
  ALLOCATE(lstep_adv(n_dom),lcall_phy(n_dom,iphysproc),linit_slowphy(n_dom), &
    &      linit_dyn(n_dom),t_elapsed_phy(n_dom,iphysproc), sim_time(n_dom), &
    &      STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ( 'mo_nh_stepping: perform_nh_stepping',           &
    &      'allocation for flow control variables failed' )
  ENDIF
  !
  ! initialize
  IF (is_restart_run()) THEN
    linit_slowphy(:)  = .FALSE.
    !
    ! Get sim_time, t_elapsed_phy and lcall_phy from restart file
    DO jg = 1,n_dom
      WRITE(attname,'(a,i2.2)') 'jstep_adv_ntsteps_DOM',jg
      CALL get_restart_attribute(TRIM(attname), jstep_adv(jg)%ntsteps)
      WRITE(attname,'(a,i2.2)') 'jstep_adv_marchuk_order_DOM',jg
      CALL get_restart_attribute(TRIM(attname), jstep_adv(jg)%marchuk_order)
      WRITE(attname,'(a,i2.2)') 'sim_time_DOM',jg
      CALL get_restart_attribute(TRIM(attname), sim_time(jg))
      DO jp = 1,iphysproc
        WRITE(attname,'(a,i2.2,a,i2.2)') 't_elapsed_phy_DOM',jg,'_PHY',jp
        CALL get_restart_attribute(TRIM(attname), t_elapsed_phy(jg,jp))

        WRITE(attname,'(a,i2.2,a,i2.2)') 'lcall_phy_DOM',jg,'_PHY',jp
        CALL get_restart_attribute(TRIM(attname), lcall_phy(jg,jp))
      ENDDO
    ENDDO
    linit_dyn(:)      = .FALSE.
  ELSE
    jstep_adv(:)%ntsteps       = 0
    jstep_adv(:)%marchuk_order = 0
    sim_time(:)                = 0._wp
    linit_slowphy(:)           = .TRUE.
    t_elapsed_phy(:,:)         = 0._wp
    linit_dyn(:)               = .TRUE.
  ENDIF


  DO jg=1, n_dom
    ALLOCATE(                                                                      &
      &  prep_adv(jg)%mass_flx_me (nproma,p_patch(jg)%nlev  ,p_patch(jg)%nblks_e), &
      &  prep_adv(jg)%mass_flx_ic (nproma,p_patch(jg)%nlevp1,p_patch(jg)%nblks_c), &
      &  prep_adv(jg)%vn_traj     (nproma,p_patch(jg)%nlev,  p_patch(jg)%nblks_e), &
      &  prep_adv(jg)%w_traj      (nproma,p_patch(jg)%nlevp1,p_patch(jg)%nblks_c), &
      &  prep_adv(jg)%rhodz_mc_now(nproma,p_patch(jg)%nlev  ,p_patch(jg)%nblks_c), &
      &  prep_adv(jg)%rhodz_mc_new(nproma,p_patch(jg)%nlev  ,p_patch(jg)%nblks_c), &
      &  prep_adv(jg)%rho_ic      (nproma,p_patch(jg)%nlevp1,p_patch(jg)%nblks_c), &
      &  prep_adv(jg)%topflx_tra  (nproma,p_patch(jg)%nblks_c,ntracer),            &
      &       STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ( 'mo_nh_stepping: perform_nh_stepping',           &
      &      'allocation for mass_flx_me, mass_flx_ic, vn_traj, ' // &
      &      'w_traj, rhodz_mc_now, rhodz_mc_new, rho_ic, '       // &
      &      'topflx_tra failed' )
    ENDIF
    !
    ! initialize (as long as restart output is synchroinzed with advection, 
    ! these variables do not need to go into the restart file)
!$OMP PARALLEL
!$OMP WORKSHARE
    prep_adv(jg)%mass_flx_me (:,:,:) = 0._wp
    prep_adv(jg)%mass_flx_ic (:,:,:) = 0._wp
    prep_adv(jg)%vn_traj     (:,:,:) = 0._wp
    prep_adv(jg)%w_traj      (:,:,:) = 0._wp
    prep_adv(jg)%rhodz_mc_now(:,:,:) = 0._wp
    prep_adv(jg)%rhodz_mc_new(:,:,:) = 0._wp
    prep_adv(jg)%rho_ic      (:,:,:) = 0._wp
    prep_adv(jg)%topflx_tra  (:,:,:) = 0._wp
!$OMP END WORKSHARE
!$OMP END PARALLEL

    ! Allocate global buffers for MPI communication
    IF (itype_comm >= 2 .AND. my_process_is_mpi_parallel()) THEN
      ALLOCATE(bufr(jg)%send_c1 (p_patch(jg)%nlevp1,   p_patch(jg)%comm_pat_c%n_send), &
        &      bufr(jg)%recv_c1 (p_patch(jg)%nlevp1,   p_patch(jg)%comm_pat_c%n_recv), &
        &      bufr(jg)%send_c3 (3*p_patch(jg)%nlev+1, p_patch(jg)%comm_pat_c%n_send), &
        &      bufr(jg)%recv_c3 (3*p_patch(jg)%nlev+1, p_patch(jg)%comm_pat_c%n_recv), &
        &      bufr(jg)%send_e1 (p_patch(jg)%nlev,     p_patch(jg)%comm_pat_e%n_send), &
        &      bufr(jg)%recv_e1 (p_patch(jg)%nlev,     p_patch(jg)%comm_pat_e%n_recv), &
        &      bufr(jg)%send_e2 (2*p_patch(jg)%nlev,   p_patch(jg)%comm_pat_e%n_send), &
        &      bufr(jg)%recv_e2 (2*p_patch(jg)%nlev,   p_patch(jg)%comm_pat_e%n_recv), &
        &      bufr(jg)%send_e3 (3*p_patch(jg)%nlev,   p_patch(jg)%comm_pat_e%n_send), &
        &      bufr(jg)%recv_e3 (3*p_patch(jg)%nlev,   p_patch(jg)%comm_pat_e%n_recv), &
        &      bufr(jg)%send_v2 (2*p_patch(jg)%nlev,   p_patch(jg)%comm_pat_v%n_send), &
        &      bufr(jg)%recv_v2 (2*p_patch(jg)%nlev,   p_patch(jg)%comm_pat_v%n_recv)  )
    ELSE
      ALLOCATE(bufr(jg)%send_c1 (1,1), &
        &      bufr(jg)%recv_c1 (1,1), &
        &      bufr(jg)%send_c3 (1,1), &
        &      bufr(jg)%recv_c3 (1,1), &
        &      bufr(jg)%send_e1 (1,1), &
        &      bufr(jg)%recv_e1 (1,1), &
        &      bufr(jg)%send_e2 (1,1), &
        &      bufr(jg)%recv_e2 (1,1), &
        &      bufr(jg)%send_e3 (1,1), &
        &      bufr(jg)%recv_e3 (1,1), &
        &      bufr(jg)%send_v2 (1,1), &
        &      bufr(jg)%recv_v2 (1,1)  )
    ENDIF

  ENDDO

  END SUBROUTINE allocate_nh_stepping

!-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  !>
  !! integrate_nh
  !!
  !! Performs dynamics time stepping:  Rotational modes (helicity bracket) and
  !! divergent modes (Poisson bracket) are splitted using Strang splitting.
  !!
  !!
  !! @par Revision History
  !! Initial release by Almut Gassmann, MPI-M (2009-08-25)
  !! Adaptation for grid refinement by Guenther Zaengl, DWD (2010-02-09)
  !! Modification by Daniel Reinert, DWD (2010-04-15)
  !!  - Implementation of tracer transport
  !! Modification by Daniel Reinert, DWD (2010-07-23)
  !!  - optional reduced calling frequency for transport and physics
  !!
#ifdef __NO_NESTING__
   SUBROUTINE integrate_nh (p_nh_state, p_patch, p_int_state,        &
  &        p_grf_state, jg, nstep_global, dt_loc, sim_time, num_steps)
#else
  RECURSIVE SUBROUTINE integrate_nh (p_nh_state, p_patch, p_int_state,  &
  &        p_grf_state, jg, nstep_global, dt_loc, sim_time, num_steps)
#endif

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = 'mo_nh_stepping:integrate_nh'

    TYPE(t_patch), TARGET, INTENT(in)    :: p_patch(n_dom_start:n_dom)    !< patch
    TYPE(t_int_state),TARGET,INTENT(in)  :: p_int_state(n_dom_start:n_dom)!< interpolation state
    TYPE(t_nh_state), TARGET, INTENT(inout) :: p_nh_state(n_dom) !< nonhydrostatic state
    TYPE(t_gridref_state), INTENT(INOUT) :: p_grf_state(n_dom_start:n_dom)!< gridref state

    INTEGER, INTENT(IN)     :: jg           !< current grid level
    INTEGER, INTENT(IN)     :: nstep_global !< counter of global time step
    INTEGER, INTENT(IN)     :: num_steps    !< number of time steps to be executed
    REAL(wp), INTENT(IN)    :: dt_loc       !< time step applicable to local grid level
    REAL(wp), INTENT(INOUT) :: sim_time(n_dom) !< elapsed simulation time on each
                                               !< grid level

    ! Local variables

    ! Time levels
    INTEGER :: n_now_grf, n_now, n_new, n_save, n_temp
    INTEGER :: n_now_rcf, n_new_rcf, n_upt_rcf  ! accounts for reduced calling frequencies (rcf)

    INTEGER :: jstep, jgp, jgc, jn
    INTEGER :: nsteps_nest ! number of time steps executed in nested domain

    REAL(wp):: dt_sub, rdt_loc
    REAL(wp):: dt_rcf, rdt_rcf   ! time step for advection and fast physics

    REAL(wp), DIMENSION(:,:,:), POINTER  :: p_vn   => NULL()
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE  :: z_tmp_e

    LOGICAL, PARAMETER :: l_straka=.FALSE.
    LOGICAL :: l_predictor
    LOGICAL :: l_bdy_nudge
    LOGICAL :: linit_vertnest(2)
    LOGICAL :: lclean_mflx   ! for reduced calling freqency: determines whether
                             ! mass-fluxes and trajectory-velocities are reset to zero
                             ! i.e. for starting new integration sweep

    ! Switch to determine manner of OpenMP parallelization in interpol_scal_grf
!     LOGICAL :: lpar_fields=.FALSE.

    ! Switch to determine if nested domains are called at a given time step
    LOGICAL :: l_call_nests = .FALSE.

!$  INTEGER :: num_threads_omp, omp_get_max_threads

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------
    ! settings for calling frequency for slow physics
    !--------------------------------------------------------------------------

    IF (jg == 1 .AND. linit_dyn(jg)) THEN

      IF (.NOT. l_nest_rcf) nsav1(1:n_dom) = nnow(1:n_dom)
    ENDIF

    !--------------------------------------------------------------------------

!$  num_threads_omp = omp_get_max_threads()

    ! Determine parent domain ID
    IF ( jg > 1) THEN
      jgp = p_patch(jg)%parent_id
    ELSE IF (n_dom_start == 0) THEN
      jgp = 0
    ELSE
      jgp = 1
    ENDIF

    ! If the limited-area mode is used, save initial state in the coarse domain
    ! The save time level is later on used for boundary relaxation in the case of
    ! fixed boundary conditions.
    ! If time-dependent data from a driving model are provided (not yet implemented),
    ! they should be written to the save time level, so that the relaxation routine
    ! automatically does the right thing

    IF (jg == 1 .AND. l_limited_area .AND. linit_dyn(jg)) THEN
      n_now = nnow(jg)
      n_save = nsav2(jg)

      WRITE(message_text,'(a)') 'save initial fields for outer boundary nudging'
       CALL message(TRIM(routine), TRIM(message_text))

      p_nh_state(jg)%prog(n_save)%vn      = p_nh_state(jg)%prog(n_now)%vn
      p_nh_state(jg)%prog(n_save)%w       = p_nh_state(jg)%prog(n_now)%w
      p_nh_state(jg)%prog(n_save)%rho     = p_nh_state(jg)%prog(n_now)%rho
      p_nh_state(jg)%prog(n_save)%theta_v = p_nh_state(jg)%prog(n_now)%theta_v

      ! tracer(nsav2) is currently not allocated in mo_nonhydro_state;
      ! thus, there is no tracer relaxation
      ! IF (ltransport) &
      !   p_nh_state(jg)%prog(n_save)%tracer = p_nh_state(jg)%prog(n_now)%tracer

    ENDIF

    ! This executes one time step for the global domain and two steps for nested domains
    DO jstep = 1, num_steps

      IF ( jstep == 1 .AND. jg > 1 ) THEN
        ! Save prognostic variables at current timestep to compute
        ! feedback increments (not needed in global domain)
        n_now = nnow(jg)
        n_save = nsav2(jg)
!$OMP PARALLEL
!$OMP WORKSHARE
        p_nh_state(jg)%prog(n_save)%vn      = p_nh_state(jg)%prog(n_now)%vn
        p_nh_state(jg)%prog(n_save)%w       = p_nh_state(jg)%prog(n_now)%w
        p_nh_state(jg)%prog(n_save)%rho     = p_nh_state(jg)%prog(n_now)%rho
        p_nh_state(jg)%prog(n_save)%theta_v = p_nh_state(jg)%prog(n_now)%theta_v
!$OMP END WORKSHARE
!$OMP END PARALLEL
      ENDIF

      ! Set 'reduced calling frequency'-time step (for transport and some physics parts)
      dt_rcf = REAL(iadv_rcf,wp)*dt_loc

      ! update several switches which decide upon
      ! - calling transport running with rcf (lstep_adv)
      ! - re-initializing temporary transport fields for rcf (lclean_mflx)
      ! - switching order of operators in case of Marchuk-splitting
      jstep_adv(jg)%ntsteps = jstep_adv(jg)%ntsteps + 1

      IF (iadv_rcf == 1) THEN
        lstep_adv(jg)  = .TRUE.   ! no reduced calling frequency
        lclean_mflx    = .TRUE.   ! no mass flux/velocity summation/averaging
        jstep_adv(jg)%marchuk_order = jstep_adv(jg)%marchuk_order + 1
      ELSE IF (MOD(jstep_adv(jg)%ntsteps,iadv_rcf) == 0 ) THEN
        lstep_adv(jg)  = .TRUE.   ! do call transport and physics
        lclean_mflx    = .FALSE.  ! do NOT re-initialize mass fluxes and velocities
        jstep_adv(jg)%marchuk_order = jstep_adv(jg)%marchuk_order + 1
      ELSE IF (MOD(jstep_adv(jg)%ntsteps,iadv_rcf) == 1 ) THEN
        lstep_adv(jg)  = .FALSE.  ! do not call transport and physics
        lclean_mflx    = .TRUE.   ! re-initialize mass fluxes and velocities
      ELSE
        lstep_adv(jg)  = .FALSE.  ! do not call transport and physics
        lclean_mflx    = .FALSE.  ! do NOT re-initialize mass fluxes and velocities
      ENDIF


      IF ( l_nest_rcf .AND. n_dom > 1) THEN
        IF (jg == 1 .AND. MOD(nstep_global,iadv_rcf) == 1 .OR. &
            jg > 1 .AND. MOD(jstep,iadv_rcf) == 1 ) THEN

          ! Save prognostic variables at current timestep to compute
          ! interpolation tendencies
          n_now  = nnow(jg)
          n_save = nsav1(jg)
!$OMP PARALLEL
!$OMP WORKSHARE
          p_nh_state(jg)%prog(n_save)%vn      = p_nh_state(jg)%prog(n_now)%vn
          p_nh_state(jg)%prog(n_save)%w       = p_nh_state(jg)%prog(n_now)%w
          p_nh_state(jg)%prog(n_save)%rho     = p_nh_state(jg)%prog(n_now)%rho
          p_nh_state(jg)%prog(n_save)%theta_v = p_nh_state(jg)%prog(n_now)%theta_v
!$OMP END WORKSHARE
!$OMP END PARALLEL
        ENDIF
      ENDIF

      ! Set local variables for time levels
      n_now  = nnow(jg)
      n_new  = nnew(jg)

      ! Set local variable for rcf-time levels
      n_now_rcf = nnow_rcf(jg)
      n_new_rcf = nnew_rcf(jg)
      ! the next time level is essential for physics packages, which are not
      ! synchronized with transport (i.e. which are not called for each advection
      ! step or for each ith advection step). Those unsynchronized physics-routines
      ! need to read from and update to timelevel n_upt_rcf and NOT n_new_rcf !!
      IF (lstep_adv(jg) ) THEN
        n_upt_rcf = nnew_rcf(jg)
      ELSE
        n_upt_rcf = nnow_rcf(jg)
      ENDIF



      IF ( l_limited_area .AND. jg == 1 ) THEN
        ! Perform interpolation of lateral boundary data from a driving model
        ! This routine still has to be written...
        ! Boundary data should be written to time level nnow, boundary tendencies
        ! and nudging increments (or, alternatively, the full fields) should be
        ! written to the grf_tend fields

!        CALL boundary_data ( p_patch(jg), p_nh_state(jg), ... )

        ! Apply nudging at the lateral boundaries if the limited-area-mode is used

        CALL outer_boundary_nudging (p_patch(jg), p_nh_state(jg), p_int_state(jg), &
          &                          nnow(jg),nnow_rcf(jg),nsav2(jg),lstep_adv(jg))
      ENDIF
      ! Note: boundary nudging in nested domains (if feedback is turned off) is
      ! applied in solve_nh for velocity components and in SR nest_boundary_nudging
      ! (called at the advective time step) for thermodynamical and tracer variables


      ! PR: The update of sim_time is moved to the beginning of the time step. 
      !  Then it has the updated value when it is passed to the physics interface subroutine
      !  Daniel, you have to adjust the advection experiments to it.
      !  I should discuss again with Thorsten if it is OK for the radiation
      !
      ! counter for simulation time in seconds
      sim_time(jg) = sim_time(jg) + dt_loc

      IF (itime_scheme == 1) THEN
        !------------------
        ! Pure advection
        !------------------

        SELECT CASE ( TRIM(nh_test_name) )

        CASE ('PA') ! solid body rotation

          ! set time-variant vertical velocity
          CALL set_nh_w_rho( p_patch(jg),p_nh_state(jg)%metrics,       &! in
            & jstep_adv(jg)%marchuk_order, dt_loc, sim_time(jg)-dt_loc,&! in
            &               p_nh_state(jg)%prog(n_new)%w,              &! inout
            &               p_nh_state(jg)%diag%pres,                  &! inout
            &               p_nh_state(jg)%diag%rho_ic                 )! inout

        CASE ('DF1', 'DF2', 'DF3', 'DF4') ! deformational flow

          ! get velocity field
          CALL get_nh_df_velocity( p_patch(jg), p_nh_state(jg)%prog(n_new), &
            &                     nh_test_name, rotate_axis_deg,            &
            &                     sim_time(jg)-dt_loc+dt_rcf )

          ! get mass flux and new \rho. The latter one is only computed,
          ! if the density equation is re-integrated.
          CALL get_nh_df_mflx_rho(p_patch(jg), p_int_state(jg), & !in
            &                     p_nh_state(jg)%prog(n_now),   & !in
            &                     p_nh_state(jg)%prog(n_new),   & !in
            &                     p_nh_state(jg)%metrics,       & !in
            &                     p_nh_state(jg)%diag, dt_rcf   ) !inout,in

        END SELECT


        ! Diagnose some velocity-related quantities for the tracer
        ! transport scheme
        CALL prepare_tracer( p_patch(jg), p_nh_state(jg)%prog(n_now),     &! in
          &         p_nh_state(jg)%prog(n_new),                           &! in
          &         p_nh_state(jg)%metrics, p_int_state(jg),              &! in
          &         iadv_rcf, lstep_adv(jg), lclean_mflx,                 &! in
          &         p_nh_state(jg)%diag,                                  &! inout
          &         prep_adv(jg)%vn_traj, prep_adv(jg)%mass_flx_me,       &! inout
          &         prep_adv(jg)%w_traj, prep_adv(jg)%mass_flx_ic,        &! inout
          &         prep_adv(jg)%rhodz_mc_now, prep_adv(jg)%rhodz_mc_new, &! inout
          &         prep_adv(jg)%rho_ic, prep_adv(jg)%topflx_tra          )! inout,out


        IF (lstep_adv(jg)) THEN
          CALL step_advection( p_patch(jg), p_int_state(jg), dt_rcf,       & !in
            &        jstep_adv(jg)%marchuk_order,                          & !in
            &        p_nh_state(jg)%prog(n_now_rcf)%tracer,                & !in
            &        prep_adv(jg)%mass_flx_me, prep_adv(jg)%vn_traj,       & !in
            &        prep_adv(jg)%mass_flx_ic, prep_adv(jg)%w_traj,        & !in
            &        p_nh_state(jg)%metrics%ddqz_z_full,                   & !in
            &        prep_adv(jg)%rhodz_mc_new, prep_adv(jg)%rhodz_mc_now, & !in
            &        p_nh_state(jg)%metrics%z_mc,                          & !in
            &        p_nh_state(jg)%metrics%z_ifc,                         & !in
            &        p_nh_state(jg)%diag%grf_tend_tracer,                  & !inout
            &        p_nh_state(jg)%prog(n_new_rcf)%tracer,                & !inout
            &        p_nh_state(jg)%diag%hfl_tracer,                       & !out
            &        p_nh_state(jg)%diag%vfl_tracer,                       & !out
            &        opt_rho_ic=prep_adv(jg)%rho_ic,                       & !in
            &        opt_topflx_tra=prep_adv(jg)%topflx_tra,               & !in
            &        opt_q_int=p_nh_state(jg)%diag%q_int,                  & !out
            &        opt_ddt_tracer_adv=p_nh_state(jg)%diag%ddt_tracer_adv ) !out
        ENDIF


      ELSE  ! itime_scheme /= 1

        ! artificial forcing (Straka)
        IF (l_straka) THEN
          CALL forcing_straka(p_nh_state(jg)%prog(n_now), p_patch(jg), p_int_state(jg), &
                              p_nh_state(jg)%metrics, p_nh_state(jg)%diag)
        ENDIF

        ! artificial forcing (Held-Suarez test forcing)
        IF ( lforcing .AND. iforcing == 1) THEN
          CALL held_suarez_nh_interface (p_nh_state(jg)%prog(n_now), p_patch(jg), &
                                         p_int_state(jg),p_nh_state(jg)%metrics,  &
                                         p_nh_state(jg)%diag)
        ENDIF


        IF ( linit_slowphy(jg) .AND. iforcing == inwp ) THEN

          CALL time_ctrl_physics ( tcall_phy, lstep_adv, dt_loc, jg,  &! in
            &                      nstep_global, .TRUE.,              &! in
            &                      t_elapsed_phy,                     &! inout
            &                      lcall_phy )                         ! out

          IF (msg_level >= 10) THEN
            WRITE(message_text,'(a,i2,11l3)') 'initial call of slow physics',jg , lcall_phy(jg,:11)
            CALL message(TRIM(routine), TRIM(message_text))
          ENDIF

          ! NOTE (DR): To me it is not clear yet, which timestep should be
          ! used for the first call of the slow_physics part dt_rcf, dt_loc,
          ! tcall_phy(jg,:) ...?
          CALL nwp_nh_interface(lcall_phy(jg,:),                   & !in
            &                  lredgrid_phys(jg),                  & !in
            &                  nstep_global,                       & !in
            &                  tcall_phy(jg,:),                    & !in
            &                  sim_time(jg),                       & !in 
            &                  p_patch(jg)  ,                      & !in
            &                  p_int_state(jg),                    & !in
            &                  p_nh_state(jg)%metrics ,            & !in
            &                  p_patch(jgp),                       & !in
            &                  p_int_state(jgp),                   & !in
            &                  p_grf_state(jgp),                   & !in
            &                  ext_data(jg)           ,            & !in
            &                  mean_charlen(jg)       ,            & !in
            &                  p_nh_state(jg)%prog(n_now) ,        & !inout
            &                  p_nh_state(jg)%prog(n_now_rcf) ,    & !inout
            &                  p_nh_state(jg)%prog(n_now_rcf) ,    & !inout
            &                  p_nh_state(jg)%diag,                & !inout
            &                  prm_diag  (jg),                     & !inout
            &                  prm_nwp_tend(jg)                ,   &
            &                  p_lnd_state(jg)%diag_lnd,           &
            &                  p_lnd_state(jg)%prog_lnd(n_now_rcf),& !inout
            &                  p_lnd_state(jg)%prog_lnd(n_now_rcf) ) !inout

          linit_slowphy(jg) = .FALSE. ! no further initialization calls needed
        ENDIF


       ! Determine which physics packages must be called/not called at the current
       ! time step
        IF ( iforcing == inwp ) THEN
          CALL time_ctrl_physics ( tcall_phy, lstep_adv, dt_loc, jg,  &! in
            &                      nstep_global, .FALSE.,             &! in
            &                      t_elapsed_phy,                     &! inout
            &                      lcall_phy )                         ! out

          IF (msg_level >= 10) THEN
            WRITE(message_text,'(a,i2,11l3)') 'call phys. packages DOM:',jg , lcall_phy(jg,:11)
            CALL message(TRIM(routine), TRIM(message_text))
            IF(ltransport) THEN
              WRITE(message_text,'(a,i4,l4)') 'call advection',jg , lstep_adv(jg)
              CALL message('integrate_nh', TRIM(message_text))
            ELSE IF(.NOT. ltransport) THEN
              WRITE(message_text,'(a,l4)') 'no advection, ltransport=', ltransport
              CALL message('integrate_nh', TRIM(message_text))
            ENDIF
          ENDIF
        ENDIF


        IF (p_patch(jg)%cell_type == 3) THEN

          IF (jg > 1 .AND. .NOT. lfeedback(jg)) THEN
            l_bdy_nudge = .TRUE. ! apply boundary nudging if feedback is turned off
          ELSE
            l_bdy_nudge = .FALSE.
          ENDIF

          IF (.NOT. l_nest_rcf .OR. lclean_mflx) THEN
            linit_vertnest(1) = .TRUE.
          ELSE
            linit_vertnest(1) = .FALSE.
          ENDIF

          IF (.NOT. l_nest_rcf .OR. lstep_adv(jg)) THEN
            linit_vertnest(2) = .TRUE.
          ELSE
            linit_vertnest(2) = .FALSE.
          ENDIF

          ! For real-data runs, perform an extra diffusion call before the first time
          ! step because no other filtering of the interpolated velocity field is done
          IF (.NOT.ltestcase .AND. linit_dyn(jg) .AND. diffusion_config(jg)%lhdiff_vn) THEN
            CALL diffusion_tria(p_nh_state(jg)%prog(n_now), p_nh_state(jg)%diag,            &
              p_nh_state(jg)%metrics, p_patch(jg), p_int_state(jg), bufr(jg), dt_loc, .TRUE.)
          ENDIF

          IF (itype_comm <= 2) THEN

            CALL solve_nh(p_nh_state(jg), p_patch(jg), p_int_state(jg), bufr(jg),       &
                        n_now, n_new, linit_dyn(jg), linit_vertnest, l_bdy_nudge, dt_loc)

          IF (diffusion_config(jg)%lhdiff_vn) &
            CALL diffusion_tria(p_nh_state(jg)%prog(n_new), p_nh_state(jg)%diag,             &
              p_nh_state(jg)%metrics, p_patch(jg), p_int_state(jg), bufr(jg), dt_loc, .FALSE.)
          ELSE
            ! call version for asynchronous halo communication, 
            ! combining solve and Smagorinsky diffusion
            IF( diffusion_config(jg)%hdiff_order /= 5) &
              CALL finish ( 'mo_nh_stepping:perform_nh_stepping',  &
              'asynchronous halo communication requires hdiff_order=5' )

            CALL solve_nh_ahc(p_nh_state(jg), p_patch(jg), p_int_state(jg), bufr(jg),   &
                        n_now, n_new, linit_dyn(jg), linit_vertnest, l_bdy_nudge, dt_loc)
          ENDIF

        ELSE ! hexagonal case

          l_predictor=.TRUE.
          ! 1. nonlinear advection terms (know) for the predictor step
          !-----------------------------------------------------------

          CALL momentum_adv(p_nh_state(jg), p_patch(jg), p_int_state(jg), n_now, n_new, &
                            l_predictor)
          ! 2. predictor step -> knew values
          !---------------------------------
          CALL divergent_modes_5band(p_nh_state(jg), p_patch(jg), p_int_state(jg), n_now, &
                                     n_new, l_predictor)

          linit_dyn(jg) = .FALSE.
          l_predictor = .FALSE.
          ! 3. nonlinear advection terms (knew) for the main step
          !------------------------------------------------------
          CALL momentum_adv(p_nh_state(jg), p_patch(jg), p_int_state(jg), n_now, &
                             n_new, l_predictor)
          ! 4. main prediction step
          !------------------------
          CALL divergent_modes_5band(p_nh_state(jg), p_patch(jg), p_int_state(jg), n_now, &
                                     n_new, l_predictor)
        ENDIF



        ! 5. tracer advection
        !-----------------------
        IF ( ltransport ) THEN

          ! Diagnose some velocity-related quantities for the tracer
          ! transport scheme
          CALL prepare_tracer( p_patch(jg), p_nh_state(jg)%prog(n_now),     &! in
            &         p_nh_state(jg)%prog(n_new),                           &! in
            &         p_nh_state(jg)%metrics, p_int_state(jg),              &! in
            &         iadv_rcf, lstep_adv(jg), lclean_mflx,                 &! in
            &         p_nh_state(jg)%diag,                                  &! inout
            &         prep_adv(jg)%vn_traj, prep_adv(jg)%mass_flx_me,       &! inout
            &         prep_adv(jg)%w_traj,  prep_adv(jg)%mass_flx_ic,       &! inout
            &         prep_adv(jg)%rhodz_mc_now, prep_adv(jg)%rhodz_mc_new, &! inout
            &         prep_adv(jg)%rho_ic, prep_adv(jg)%topflx_tra          )! inout,out


          IF (lstep_adv(jg)) THEN

            CALL step_advection( p_patch(jg), p_int_state(jg), dt_rcf,         & !in
              &          jstep_adv(jg)%marchuk_order,                          & !in
              &          p_nh_state(jg)%prog(n_now_rcf)%tracer,                & !in
              &          prep_adv(jg)%mass_flx_me, prep_adv(jg)%vn_traj,       & !in
              &          prep_adv(jg)%mass_flx_ic, prep_adv(jg)%w_traj,        & !in
              &          p_nh_state(jg)%metrics%ddqz_z_full,                   & !in
              &          prep_adv(jg)%rhodz_mc_new, prep_adv(jg)%rhodz_mc_now, & !in
              &          p_nh_state(jg)%metrics%z_mc,                          & !in
              &          p_nh_state(jg)%metrics%z_ifc,                         & !in
              &          p_nh_state(jg)%diag%grf_tend_tracer,                  & !inout
              &          p_nh_state(jg)%prog(n_new_rcf)%tracer,                & !inout
              &          p_nh_state(jg)%diag%hfl_tracer,                       & !out
              &          p_nh_state(jg)%diag%vfl_tracer,                       & !out
              &          opt_rho_ic=prep_adv(jg)%rho_ic,                       & !in
              &          opt_topflx_tra=prep_adv(jg)%topflx_tra,               & !in
              &          opt_q_int=p_nh_state(jg)%diag%q_int,                  & !out
              &          opt_ddt_tracer_adv=p_nh_state(jg)%diag%ddt_tracer_adv ) !out


!            IF (  iforcing==inwp .AND. inwp_turb == 1) THEN
!              !> KF preliminary relabeling of TKE as long as there is no advection for it
!              p_nh_state(jg)%prog(n_new_rcf)%tke =  p_nh_state(jg)%prog(n_now)%tke
!            ENDIF

           ENDIF  !lstep_adv

        ENDIF

        ! Apply boundary nudging in case of one-way nesting
        IF (lstep_adv(jg) .AND. .NOT. lfeedback(jg)) THEN
          CALL nest_boundary_nudging(p_patch(jg), p_nh_state(jg), p_int_state(jg), &
            &                        nnew(jg),nnew_rcf(jg),REAL(iadv_rcf,wp))
        ENDIF

        IF (  iforcing==inwp .AND. lstep_adv(jg) ) THEN

          !> moist tracer update is now synchronized with advection and satad

          CALL nwp_nh_interface(lcall_phy(jg,:),                   & !in
            &                  lredgrid_phys(jg),                  & !in
            &                  nstep_global,                       & !in
            &                  t_elapsed_phy(jg,:),                & !in
            &                  sim_time(jg),                       & !in  
            &                  p_patch(jg)  ,                      & !in
            &                  p_int_state(jg),                    & !in
            &                  p_nh_state(jg)%metrics ,            & !in
            &                  p_patch(jgp),                       & !in
            &                  p_int_state(jgp),                   & !in
            &                  p_grf_state(jgp),                   & !in
            &                  ext_data(jg)           ,            & !in
            &                  mean_charlen(jg)       ,            & !in
            &                  p_nh_state(jg)%prog(n_new) ,        & !inout
            &                  p_nh_state(jg)%prog(n_now_rcf),     & !in for tke
            &                  p_nh_state(jg)%prog(n_new_rcf) ,    & !inout
            &                  p_nh_state(jg)%diag ,               & !inout
            &                  prm_diag  (jg),                     & !inout
            &                  prm_nwp_tend(jg),                   &
            &                  p_lnd_state(jg)%diag_lnd,           &
            &                  p_lnd_state(jg)%prog_lnd(n_now_rcf),& !inout
            &                  p_lnd_state(jg)%prog_lnd(n_new_rcf) ) !inout

        ENDIF !iforcing


        ! This calls 4th order diffusion for the hexagonal case
        IF (diffusion_config(jg)%lhdiff_vn .AND. &
          &  p_patch(jg)%cell_type == 6) THEN
          CALL diffusion_hex(p_nh_state(jg)%prog(n_new), p_nh_state(jg)%metrics,&
          &                  p_patch(jg), p_int_state(jg), dt_loc)
        ENDIF

      ENDIF  ! itime_scheme

      ! If there are nested domains...
      IF (l_nest_rcf .AND. lstep_adv(jg))  THEN
        l_call_nests = .TRUE.
        rdt_loc = 1._wp/(dt_loc*REAL(iadv_rcf,wp))
        n_now_grf    = nsav1(jg)
        nsteps_nest  = 2*iadv_rcf
      ELSE IF (.NOT. l_nest_rcf) THEN
        l_call_nests = .TRUE.
        rdt_loc = 1._wp/dt_loc
        n_now_grf    = n_now
        nsteps_nest  = 2
      ELSE
        l_call_nests = .FALSE.
      ENDIF

#ifndef __NO_NESTING__
      IF (l_call_nests .AND. p_patch(jg)%n_childdom > 0) THEN

        dt_sub = dt_loc/2._wp
        rdt_rcf = 1._wp/dt_rcf

        ! Compute time tendencies for interpolation to refined mesh boundaries
        CALL compute_tendencies (p_patch(jg),p_nh_state(jg),n_new,n_now_grf,n_new_rcf, &
          &                      n_now_rcf,rdt_loc,rdt_rcf,lstep_adv(jg))

        ! Please note:
        ! The use of start_delayed_exchange/do_delayed_exchange is not restricted
        ! to processor splitting only. If it turns out that the code below
        ! works faster with these routines also in the case that processor
        ! splitting is not in effect, it can be used always.
        ! start_delayed_exchange/do_delayed_exchang can also be removed completely
        ! if it should have any negative effects.

        ! GZ: it turned out that the delayed exchange is detrimental on the NEC
#ifndef __SX__
        IF(proc_split) CALL start_delayed_exchange ! Data exchanges will be buffered
#endif

        ! Loop over nested domains
        DO jn = 1, p_patch(jg)%n_childdom

          jgc = p_patch(jg)%child_id(jn)

          ! Interpolate tendencies to lateral boundaries of refined mesh (jgc)
          CALL boundary_interpolation(p_patch,p_nh_state,p_int_state,p_grf_state, &
            &     jg,jgc,n_now_grf,nnow(jgc),n_now_rcf,nnow_rcf(jgc),lstep_adv(jg))

        ENDDO

#ifndef __SX__
        IF(proc_split) CALL do_delayed_exchange ! actually execute exchanges
#endif

        ! prep_bdy_nudging can not be called using delayed requests!
        DO jn = 1, p_patch(jg)%n_childdom

          jgc = p_patch(jg)%child_id(jn)
          ! If feedback is turned off for child domain, compute parent-child
          ! differences for boundary nudging
          ! *** prep_bdy_nudging adapted for reduced calling frequency of tracers ***
          IF (.NOT. lfeedback(jgc)) THEN
            CALL prep_bdy_nudging( p_patch,p_nh_state,p_int_state,p_grf_state, &
              &                   jg,jgc,lstep_adv(jg) )
          ENDIF
        ENDDO

        DO jn = 1, p_patch(jg)%n_childdom

          jgc = p_patch(jg)%child_id(jn)

          IF(p_patch(jgc)%n_patch_cells > 0) THEN
            IF(proc_split) CALL push_glob_comm(p_patch(jgc)%comm, p_patch(jgc)%proc0)
            ! Recursive call to process_grid_level for child grid level
            CALL integrate_nh( p_nh_state, p_patch, p_int_state,           &
              p_grf_state, jgc, nstep_global, dt_sub, sim_time, nsteps_nest)
            IF(proc_split) CALL pop_glob_comm()
          ENDIF

        ENDDO

        DO jn = 1, p_patch(jg)%n_childdom

          ! Call feedback to copy averaged prognostic variables from refined mesh back
          ! to the coarse mesh (i.e. from jgc to jg)

          jgc = p_patch(jg)%child_id(jn)
          IF (lfeedback(jgc)) THEN
            CALL feedback(p_patch, p_nh_state, p_int_state, p_grf_state, jgc, &
                          jg, lstep_adv(jg))
            ! Note: the last argument of "feedback" ensures that tracer feedback is
            ! only done for those time steps in which transport and microphysics are called
          ENDIF
        ENDDO

      ENDIF
#endif

      IF (l_outputtime) THEN ! compute diagnostic quantities
        p_vn  => p_nh_state(jg)%prog(n_new)%vn
        SELECT CASE (p_patch(jg)%cell_type)
        CASE (3)
          CALL rbf_vec_interpol_cell(p_vn,p_patch(jg),p_int_state(jg),&
                                     p_nh_state(jg)%diag%u,p_nh_state(jg)%diag%v)
          CALL div_avg(p_vn, p_patch(jg), p_int_state(jg), p_int_state(jg)%c_bln_avg, &
                 p_nh_state(jg)%diag%div)

          ! Fill boundaries of nested domains
          IF (p_patch(jg)%n_childdom > 0) THEN
            CALL sync_patch_array_mult(SYNC_C, p_patch(jg), 3, p_nh_state(jg)%diag%u,      &
              p_nh_state(jg)%diag%v, p_nh_state(jg)%diag%div)
          ENDIF
          DO jn = 1, p_patch(jg)%n_childdom
            jgc = p_patch(jg)%child_id(jn)

            CALL interpol_scal_grf (p_patch(jg), p_patch(jgc), p_int_state(jg),       &
                 p_grf_state(jg)%p_dom(jn), jn, 2, p_nh_state(jg)%diag%u,             &
                 p_nh_state(jgc)%diag%u, p_nh_state(jg)%diag%v, p_nh_state(jgc)%diag%v)

            CALL interpol_scal_grf (p_patch(jg), p_patch(jgc), p_int_state(jg),       &
                 p_grf_state(jg)%p_dom(jn), jn, 1, p_nh_state(jg)%diag%div,           &
                 p_nh_state(jgc)%diag%div)

            IF ( iforcing == inwp ) THEN
              CALL interpol_phys_grf(p_patch(jg), p_patch(jgc),                  &
                                     p_int_state(jg), p_grf_state(jg)%p_dom(jn), &
                                     jg, jgc, jn )

              IF (lfeedback(jgc)) CALL feedback_phys_diag(p_patch, p_grf_state, jgc, jg)

              ! Fill lateral boundaries of TKE field; note: time-level-switching has
              ! already been done at child level, but not yet at parent level
              CALL sync_patch_array(SYNC_C, p_patch(jg), p_nh_state(jg)%prog(nnew_rcf(jg))%tke)

              CALL interpol_scal_grf (p_patch(jg), p_patch(jgc), p_int_state(jg),        &
                 p_grf_state(jg)%p_dom(jn), jn, 1, p_nh_state(jg)%prog(nnew_rcf(jg))%tke,&
                 p_nh_state(jgc)%prog(nnow_rcf(jgc))%tke)

            ENDIF

          ENDDO
        CASE (6)
          CALL edges2cells_scalar(p_vn,p_patch(jg),p_int_state(jg)%hex_east ,&
                                  p_nh_state(jg)%diag%u)
          CALL edges2cells_scalar(p_vn,p_patch(jg),p_int_state(jg)%hex_north,&
                                  p_nh_state(jg)%diag%v)
          CALL div(p_vn, p_patch(jg), p_int_state(jg), p_nh_state(jg)%diag%div)

          ALLOCATE(z_tmp_e(nproma,p_patch(jg)%nlev,p_patch(jg)%nblks_e))
          CALL verts2edges_scalar(p_nh_state(jg)%diag%omega_z,p_patch(jg), &
          &                       p_int_state(jg)%tria_aw_rhom,z_tmp_e)
          CALL sync_patch_array(SYNC_E,p_patch(jg),z_tmp_e)
          CALL edges2verts_scalar(z_tmp_e,p_patch(jg),p_int_state(jg)%e_1o3_v,&
                                  p_nh_state(jg)%diag%omega_z)
          DEALLOCATE(z_tmp_e)

        END SELECT
      ENDIF
      IF ((l_outputtime.AND.(p_patch(jg)%cell_type==3)).OR.&
        & (p_patch(jg)%cell_type==6)) THEN
        CALL diagnose_pres_temp (p_nh_state(jg)%metrics, p_nh_state(jg)%prog(nnew(jg)), &
          &                      p_nh_state(jg)%prog(nnew_rcf(jg)),                     &
          &                      p_nh_state(jg)%diag,p_patch(jg),                       &
          &                      opt_calc_temp=.TRUE.,                                  &
          &                      opt_calc_pres=.TRUE.                                 )
      ENDIF


      ! Finally, switch between time levels now and new for next time step
      n_temp   = nnow(jg)
      nnow(jg) = nnew(jg)
      IF (.NOT. l_nest_rcf) nsav1(jg) = nnow(jg)
      nnew(jg) = n_temp

      ! Special treatment for processes (i.e. advection) which can be treated with
      ! reduced calling frequency. Switch between time levels now and new immediately
      ! AFTER the last transport timestep.
      IF (lstep_adv(jg)) THEN
        n_temp       = nnow_rcf(jg)
        nnow_rcf(jg) = nnew_rcf(jg)
        nnew_rcf(jg) = n_temp
      ENDIF

    ENDDO

  END SUBROUTINE integrate_nh
  !-----------------------------------------------------------------------------
  

  !-------------------------------------------------------------------------
  !>
  !! Physics time control
  !!
  !! Time control for slow and fast physics. This function provides a 2D array
  !! of type LOGICAL. For each physical process there is one column of length
  !! n_dom. A physical process (iphys) on domain (jg) will be called, if
  !! lcall_phy(jg,iphys)=.TRUE.. Whether it is .TRUE. or .FALSE. depends
  !! on the current time and the prescribed calling period listed in
  !! tcall_phy.
  !!
  !! @par Revision History
  !! Developed by Daniel Reinert, DWD (2010-09-23)
  !!
  SUBROUTINE time_ctrl_physics ( tcall_phy, lstep_adv, dt_loc, jg, nstep_global, &
    &                            linit, t_elapsed_phy, lcall_phy )

    REAL(wp), INTENT(IN)    ::   &      !< Field of calling-time interval (seconds) for
      &  tcall_phy(:,:)                 !< each domain and physical process

    LOGICAL, INTENT(IN)     ::   &      !< determines whether this is a timestep with
      &  lstep_adv(:)                   !< (.TRUE.) or without (.FALSE.) scalar transport

    LOGICAL, INTENT(IN)     ::   &      !< special initialization of lcall_phy and
      &  linit                          !< t_elapsed_phy for the first call of
                                        !< nwp_nh_interface before the first dynamcs step

    REAL(wp), INTENT(INOUT) ::   &      !< elapsed time after the last call of physics
      &  t_elapsed_phy(:,:)             !< packages

    LOGICAL, INTENT(OUT)    ::   &
      &  lcall_phy(:,:)

    REAL(wp), INTENT(IN) :: dt_loc      !< dynamics time step

    INTEGER, INTENT(IN) :: jg           !< domain number

    INTEGER, INTENT(IN) :: nstep_global !< counter of global time step

    INTEGER :: ip                       !< loop index


  !-------------------------------------------------------------------------

    ! special treatment for the first physics call prior to the very first
    ! dynamics step
    IF (linit) THEN

      ! Initialize lcall_phy as false. Only those set explicitly below are true
      lcall_phy(jg,:)  = .FALSE.

      ! slow physics
      IF ( tcall_phy(jg,itconv) > 0._wp ) lcall_phy(jg,itconv) = .TRUE.

      IF ( tcall_phy(jg,itccov) > 0._wp ) lcall_phy(jg,itccov) = .TRUE.

      IF ( tcall_phy(jg,itrad) > 0._wp )  lcall_phy(jg,itrad) = .TRUE.

      IF ( tcall_phy(jg,itradheat) > 0._wp )  lcall_phy(jg,itradheat) = .TRUE.

      IF ( tcall_phy(jg,itsso) > 0._wp )  lcall_phy(jg,itsso) = .TRUE.
 
    ELSE
      !
      ! all physics packages (except saturation adjustment) are forced to run
      ! at a multiple of the advective time step
      !
      DO ip = 1, iphysproc

        ! If a physics package has been called at previous timestep,
        ! reset the time counter.
        IF ( lcall_phy(jg,ip) ) THEN
           t_elapsed_phy(jg,ip)   = 0._wp
        ENDIF

        ! update time counter
        t_elapsed_phy(jg,ip) = t_elapsed_phy(jg,ip) + dt_loc  ! dynamics timestep !!


        IF (.NOT. lstep_adv(jg) .AND. ip /= itsatad ) THEN
          lcall_phy(jg,ip) = .FALSE.
        ELSE
          IF( t_elapsed_phy(jg,ip) >= tcall_phy(jg,ip) .AND. &
            & tcall_phy(jg,ip) > 0._wp) THEN

            ! Please note: tcall_phy(jg,ip) = 0._wp indicates that the process
            ! with number ip will not be called at all (see mo_nwp_phy_nml)
            lcall_phy(jg,ip)       = .TRUE.
          ELSE
            lcall_phy(jg,ip)       = .FALSE.
          ENDIF

        ENDIF

      ENDDO
    ENDIF

  END SUBROUTINE time_ctrl_physics

END MODULE mo_nh_stepping

