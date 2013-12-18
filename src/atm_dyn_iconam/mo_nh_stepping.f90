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

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_nh_stepping
!-------------------------------------------------------------------------
!
!    ProTeX FORTRAN source: Style 2
!    modified for ICON project, DWD/MPI-M 2006
!
!-------------------------------------------------------------------------
!
!

  USE mo_kind,                ONLY: wp, i8
  USE mo_nonhydro_state,      ONLY: p_nh_state
  USE mo_nonhydrostatic_config,ONLY: iadv_rcf, lhdiff_rcf, l_nest_rcf, itime_scheme, &
    & nest_substeps, divdamp_order, divdamp_fac, divdamp_fac_o2
  USE mo_diffusion_config,     ONLY: diffusion_config
  USE mo_dynamics_config,      ONLY: nnow,nnew, nnow_rcf, nnew_rcf, nsav1, nsav2
  USE mo_io_config,            ONLY: l_outputtime, l_diagtime, is_checkpoint_time
  USE mo_parallel_config,      ONLY: nproma, itype_comm, iorder_sendrecv, use_async_restart_output
  USE mo_run_config,           ONLY: ltestcase, dtime, dtime_adv, nsteps,     &
    &                                ldynamics, ltransport, ntracer, lforcing, iforcing, &
    &                                msg_level, test_mode, output_mode
  USE mo_radiation_config,     ONLY: albedo_type
  USE mo_timer,               ONLY: ltimer, timers_level, timer_start, timer_stop,   &
    &                               timer_total, timer_model_init, timer_nudging,    &
    &                               timer_bdy_interp, timer_feedback, timer_nesting, &
    &                               timer_integrate_nh, timer_nh_diagnostics
  USE mo_atm_phy_nwp_config,  ONLY: dt_phy, atm_phy_nwp_config
  USE mo_nwp_phy_init,        ONLY: init_nwp_phy
  USE mo_nwp_phy_state,       ONLY: prm_diag, prm_nwp_tend, phy_params, prm_nwp_diag_list, &
                                    prm_nwp_tend_list
  USE mo_lnd_nwp_config,      ONLY: nlev_soil, nlev_snow, sstice_mode, lseaice
  USE mo_nwp_lnd_state,       ONLY: p_lnd_state
  USE mo_ext_data_state,      ONLY: ext_data, interpol_monthly_mean
  USE mo_lnd_jsbach_config,   ONLY: lnd_jsbach_config
  USE mo_extpar_config,       ONLY: itopo
  USE mo_limarea_config,      ONLY: latbc_config
  USE mo_model_domain,        ONLY: p_patch
  USE mo_time_config,         ONLY: time_config
  USE mo_grid_config,         ONLY: n_dom, lfeedback, ifeedback_type, l_limited_area, &
    &                               n_dom_start, lredgrid_phys, start_time, end_time
  USE mo_nh_testcases,        ONLY: init_nh_testtopo, init_nh_testcase 
  USE mo_nh_testcases_nml,    ONLY: nh_test_name, rotate_axis_deg, lcoupled_rho
  USE mo_nh_pa_test,          ONLY: set_nh_w_rho
  USE mo_nh_df_test,          ONLY: get_nh_df_velocity
  USE mo_nh_supervise,        ONLY: supervise_total_integrals_nh, print_maxwinds
  USE mo_intp_data_strc,      ONLY: t_int_state, t_lon_lat_intp, p_int_state
  USE mo_intp_rbf,            ONLY: rbf_vec_interpol_cell
  USE mo_intp,                ONLY: edges2cells_scalar, verts2edges_scalar, edges2verts_scalar, &
    &                               verts2cells_scalar
  USE mo_grf_intp_data_strc,  ONLY: p_grf_state
  USE mo_gridref_config,      ONLY: l_density_nudging, grf_intmethod_e
  USE mo_grf_bdyintp,         ONLY: interpol_scal_grf
  USE mo_nh_nest_utilities,   ONLY: compute_tendencies, boundary_interpolation,    &
                                    complete_nesting_setup, prep_bdy_nudging,      &
                                    outer_boundary_nudging, nest_boundary_nudging, &
                                    prep_rho_bdy_nudging, density_boundary_nudging,&
                                    prep_outer_bdy_nudging
  USE mo_nh_feedback,         ONLY: feedback, relax_feedback
  USE mo_datetime,            ONLY: t_datetime, print_datetime, add_time, check_newday, &
                                    rdaylen, date_to_time
  USE mo_io_restart,          ONLY: write_restart_info_file, create_restart_file
  USE mo_exception,           ONLY: message, message_text, finish
  USE mo_impl_constants,      ONLY: SUCCESS, MAX_CHAR_LENGTH, iphysproc,    &
    &                               iphysproc_short, itconv, itccov, itrad, &
    &                               itradheat, itsso, itsatad, itgwd, inwp, iecham, &
    &                               itupdate, itturb, itgscp, itsfc, min_rlcell_int, &
                                    min_rledge_int, MODE_DWDANA, MODIS, icosmo
  USE mo_math_divrot,         ONLY: div, div_avg, rot_vertex
  USE mo_solve_nonhydro,      ONLY: solve_nh
  USE mo_update_dyn,          ONLY: add_slowphys
  USE mo_advection_stepping,  ONLY: step_advection
  USE mo_integrate_density_pa,ONLY: integrate_density_pa
  USE mo_nh_dtp_interface,    ONLY: prepare_tracer
  USE mo_nh_diffusion,        ONLY: diffusion
  USE mo_mpi,                 ONLY: my_process_is_stdio, my_process_is_mpi_parallel, &
    &                               proc_split, push_glob_comm, pop_glob_comm,       &
    &                               get_my_mpi_all_id
#ifdef NOMPI
  USE mo_mpi,                 ONLY: my_process_is_mpi_all_seq
#endif
  
  USE mo_sync,                ONLY: sync_patch_array_mult, &
                                    SYNC_C, SYNC_E, sync_patch_array
  USE mo_nh_interface_nwp,    ONLY: nwp_nh_interface
  USE mo_phys_nest_utilities, ONLY: interpol_phys_grf, feedback_phys_diag, interpol_rrg_grf
  USE mo_vertical_grid,       ONLY: set_nh_metrics
  USE mo_nh_diagnose_pres_temp,ONLY: diagnose_pres_temp
  USE mo_nh_held_suarez_interface, ONLY: held_suarez_nh_interface
  USE mo_master_control,      ONLY: is_restart_run
  USE mo_io_restart_attributes,ONLY: get_restart_attribute

!   USE mo_nwp_mpiomp_rrtm_interface, ONLY: nwp_start_radiation_ompthread, model_end_ompthread, &
!     & init_ompthread_radiation
!   USE mo_parallel_config,     ONLY: parallel_radiation_omp, nh_stepping_ompthreads
  USE mo_meteogram_config,    ONLY: meteogram_output_config
  USE mo_meteogram_output,    ONLY: meteogram_sample_vars, meteogram_is_sample_step
  USE mo_name_list_output,    ONLY: write_name_list_output, istime4name_list_output
  USE mo_pp_scheduler,        ONLY: new_simulation_status, pp_scheduler_process
  USE mo_pp_tasks,            ONLY: t_simulation_status
  USE mo_art_emission_interface, ONLY: art_emission_interface
  USE mo_art_sedi_interface,  ONLY: art_sedi_interface
  USE mo_art_tools_interface, ONLY: art_tools_interface
  USE mo_art_config,          ONLY: art_config
  USE mo_nwp_sfc_utils,       ONLY: aggregate_landvars, update_sstice, update_ndvi
  USE mo_nh_init_nest_utils,  ONLY: initialize_nest, topo_blending_and_fbk
  USE mo_nh_init_utils,       ONLY: hydro_adjust_downward
  USE mo_td_ext_data,         ONLY: set_actual_td_ext_data,  &
                                  & read_td_ext_data_file
  USE mo_initicon_config,     ONLY: init_mode
  USE mo_ls_forcing_nml,      ONLY: is_ls_forcing
  USE mo_ls_forcing,          ONLY: init_ls_forcing
  USE mo_nh_latbc,            ONLY: prepare_latbc_data , read_latbc_data, &
    &                               deallocate_latbc_data, p_latbc_data,   &
    &                               read_latbc_tlev, last_latbc_tlev, &
    &                               update_lin_interc
  USE mo_interface_les,       ONLY: les_phy_interface
  USE mo_io_restart_async,    ONLY: prepare_async_restart, write_async_restart, &
    &                               close_async_restart, set_data_async_restart
  USE mo_nh_prepadv_types,    ONLY: prep_adv, jstep_adv

#ifdef MESSY
  USE messy_main_channel_bi,   ONLY: messy_channel_write_output &
    &                              , IOMODE_RST
#ifdef MESSYTIMER
  USE messy_main_timer_bi,     ONLY: messy_timer_reset_time 
#endif
#endif

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER :: version = '$Id$'



  ! additional flow control variables that need to be dimensioned with the
  ! number of model domains
  LOGICAL, ALLOCATABLE :: lstep_adv(:)   ! determines whether tracer continuity equations
                                         ! should be integrated (.true.) or not (.false.)

  LOGICAL, ALLOCATABLE :: lcall_phy(:,:) ! contains information which physics package
                                         ! must be called at the current timestep
                                         ! and on the current domain.

  REAL(wp), ALLOCATABLE :: t_elapsed_phy(:,:)  ! time (in s) since the last call of
                                               ! the corresponding physics package
                                               ! (fast physics packages are treated as one)

  LOGICAL, ALLOCATABLE :: linit_slowphy(:) ! determines whether slow physics has already been initialized

  LOGICAL, ALLOCATABLE :: linit_dyn(:)  ! determines whether dynamics has already been initialized


  ! additional time control variables which are not dimensioned with the number 
  ! of model domains
  LOGICAL :: map_phyproc(iphysproc,iphysproc_short) !< mapping matrix
  INTEGER :: iproclist(iphysproc)  !< x-axis of mapping matrix

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
  SUBROUTINE prepare_nh_integration
!
  INTEGER :: ntl

!-----------------------------------------------------------------------

  ! for the split explict scheme, ntl is always 2
  ntl = 2

  IF (ltestcase) THEN
    CALL init_nh_testtopo(p_patch(1:), ext_data)   ! set analytic topography
  ENDIF

  IF (n_dom > 1) CALL topo_blending_and_fbk(1)

  CALL set_nh_metrics(p_patch(1:), p_nh_state, p_int_state(1:), ext_data)

  IF (n_dom > 1) THEN
    CALL complete_nesting_setup()
  ENDIF

  IF (ltestcase) THEN
    CALL init_nh_testcase(p_patch(1:), p_nh_state, p_int_state(1:), p_lnd_state(1:), &
      & ext_data, ntl)
     
    IF(is_ls_forcing) &
       CALL init_ls_forcing(p_nh_state(1)%metrics)
  ENDIF

  CALL setup_time_ctrl_physics( )

  END SUBROUTINE prepare_nh_integration
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !>
  !! Organizes nonhydrostatic time stepping
  !! Currently we assume to have only one grid level.
  !!
  !! @par Revision History
  !! Initial release by Almut Gassmann, (2009-04-15)
  !!
  SUBROUTINE perform_nh_stepping (datetime, n_checkpoint, n_diag )
!
  INTEGER, INTENT(IN)                          :: n_checkpoint, n_diag


  TYPE(t_datetime), INTENT(INOUT)      :: datetime
  TYPE(t_simulation_status)            :: simulation_status

!!$  CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
!!$    &  routine = 'mo_nh_stepping:perform_nh_stepping'

  INTEGER                              :: jg

!$  INTEGER omp_get_num_threads
!$  INTEGER omp_get_max_threads
!$  INTEGER omp_get_max_active_levels
!-----------------------------------------------------------------------

  IF (timers_level > 3) CALL timer_start(timer_model_init)

  CALL allocate_nh_stepping ()

  ! Compute diagnostic dynamics fields for initial output and physics initialization
  IF (.NOT.is_restart_run()) CALL diag_for_output_dyn (linit=.TRUE.)

  IF (sstice_mode > 1 .AND. iforcing == inwp) THEN
    ! t_seasfc and fr_seaice have to be set again from the ext_td_data files
    !  the values from the analysis have to be overwritten
    CALL set_actual_td_ext_data (.TRUE.,datetime,datetime,sstice_mode,  &
                                &  p_patch(1:), ext_data, p_lnd_state)
  END IF

  IF (iforcing == inwp .AND. is_restart_run()) THEN
    DO jg=1, n_dom
      IF (.NOT. p_patch(jg)%ldom_active) CYCLE
      CALL init_nwp_phy( dtime                     ,&
           & p_patch(jg)                           ,&
           & p_nh_state(jg)%metrics                ,&
           & p_nh_state(jg)%prog(nnow(jg))         ,&
           & p_nh_state(jg)%diag                   ,&
           & prm_diag(jg)                          ,&
           & prm_nwp_tend(jg)                      ,&
           & p_lnd_state(jg)%prog_lnd(nnow_rcf(jg)),&
           & p_lnd_state(jg)%prog_lnd(nnew_rcf(jg)),&
           & p_lnd_state(jg)%prog_wtr(nnow_rcf(jg)),&
           & p_lnd_state(jg)%prog_wtr(nnew_rcf(jg)),&
           & p_lnd_state(jg)%diag_lnd              ,&
           & ext_data(jg)                          ,&
           & phy_params(jg)                         )
    ENDDO
  ELSE IF (iforcing == inwp) THEN ! for cold start, use atmospheric fields at time level nnow only
    DO jg=1, n_dom
      IF (.NOT. p_patch(jg)%ldom_active) CYCLE
      CALL init_nwp_phy( dtime                     ,&
           & p_patch(jg)                           ,&
           & p_nh_state(jg)%metrics                ,&
           & p_nh_state(jg)%prog(nnow(jg))         ,&
           & p_nh_state(jg)%diag                   ,&
           & prm_diag(jg)                          ,&
           & prm_nwp_tend(jg)                      ,&
           & p_lnd_state(jg)%prog_lnd(nnow_rcf(jg)),&
           & p_lnd_state(jg)%prog_lnd(nnew_rcf(jg)),&
           & p_lnd_state(jg)%prog_wtr(nnow_rcf(jg)),&
           & p_lnd_state(jg)%prog_wtr(nnew_rcf(jg)),&
           & p_lnd_state(jg)%diag_lnd              ,&
           & ext_data(jg)                          ,&
           & phy_params(jg)                         )
    ENDDO
    ! Compute diagnostic physics fields
    CALL diag_for_output_phys
    ! Initial call of (slow) physics schemes, including computation of transfer coefficients
    CALL init_slowphysics (datetime, 1, dtime, dtime_adv, time_config%sim_time)
  ENDIF

  !------------------------------------------------------------------
  !  get and write out some of the initial values
  !------------------------------------------------------------------
  IF (.NOT.is_restart_run()) THEN

    !--------------------------------------------------------------------------
    ! loop over the list of internal post-processing tasks, e.g.
    ! interpolate selected fields to p- and/or z-levels
    simulation_status = new_simulation_status(l_first_step   = .TRUE.,                  &
      &                                       l_output_step  = .TRUE.,                  &
      &                                       l_dom_active   = p_patch(1:)%ldom_active, &
      &                                       i_timelevel    = nnow)
    CALL pp_scheduler_process(simulation_status)

    IF (output_mode%l_nml) THEN
      CALL write_name_list_output(jstep=0)
    END IF

#ifdef MESSY
    ! MESSy initial output
!    CALL messy_write_output
#endif

  END IF ! not is_restart_run()


  IF (timers_level > 3) CALL timer_stop(timer_model_init)

!   IF (parallel_radiation_omp) THEN
! 
!     !---------------------------------------
!     CALL init_ompthread_radiation()
!     
! !$    CALL omp_set_nested(.true.)
! !$    CALL omp_set_num_threads(2)
! !$    write(0,*) 'omp_get_max_active_levels=',omp_get_max_active_levels
! !$    write(0,*) 'omp_get_max_threads=',omp_get_max_threads()
! !$OMP PARALLEL SECTIONS
! !$OMP SECTION
! !$  CALL omp_set_num_threads(nh_stepping_ompthreads)
! !$    write(0,*) 'This is the nh_timeloop, max threads=',omp_get_max_threads()
! !$    write(0,*) 'omp_get_num_threads=',omp_get_num_threads()
! 
!     CALL perform_nh_timeloop (datetime, jfile, n_checkpoint, n_diag, l_have_output )
!     CALL model_end_ompthread()
! 
! !$OMP SECTION
! !$  write(0,*) 'This is the nwp_parallel_radiation_thread, max threads=',&
! !$    omp_get_max_threads()
!   CALL nwp_start_radiation_ompthread()
! !$OMP END PARALLEL SECTIONS
! 
!   ELSE
    !---------------------------------------

    CALL perform_nh_timeloop (datetime, n_checkpoint, n_diag )
!   ENDIF

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
  SUBROUTINE perform_nh_timeloop (datetime, n_checkpoint, n_diag )
!
    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = 'mo_nh_stepping:perform_nh_timeloop'

  INTEGER, INTENT(IN)                          :: n_checkpoint, n_diag

  TYPE(t_datetime), INTENT(INOUT)      :: datetime

  INTEGER                              :: jstep, jg, kstep
  INTEGER                              :: ierr
  LOGICAL                              :: l_compute_diagnostic_quants,  &
    &                                     l_nml_output, &
    &                                     l_supervise_total_integrals,  &
    &                                     lwrite_checkpoint, ldom_active(n_dom)
  TYPE(t_simulation_status)            :: simulation_status

  TYPE(t_datetime)                     :: datetime_old

  INTEGER           :: nsoil(n_dom), nsnow
  REAL(wp)          :: elapsed_time_global
  INTEGER           :: jstep0 ! start counter for time loop

!$  INTEGER omp_get_num_threads
!-----------------------------------------------------------------------

  IF (ltimer) CALL timer_start(timer_total)

  lwrite_checkpoint = .FALSE.

  ! Prepare number of soil/snow layers for TERRA/JSBACH to be used for restart file creation below.
  ! AD: Initialize with 0 to avoid errors with certain compilers
  nsoil(:) = 0
  nsnow    = 0
  IF (iforcing == inwp) THEN
    DO jg=1,n_dom
      nsoil(jg) = nlev_soil
      nsnow     = nlev_snow
    END DO
  ELSE IF (iforcing == iecham) THEN
    DO jg=1,n_dom
      nsoil(jg) = lnd_jsbach_config(jg)%nsoil
      nsnow     = 0
    END DO
  END IF

  ! If the testbed mode is selected, reset iorder_sendrecv to 0 in order to suppress
  ! MPI communication from now on. 
  IF (test_mode > 0) iorder_sendrecv = 0

  datetime_old = datetime

  IF (use_async_restart_output) THEN
    CALL prepare_async_restart(opt_t_elapsed_phy_size = SIZE(t_elapsed_phy, 2), &
      &                        opt_lcall_phy_size     = SIZE(lcall_phy, 2))
  ENDIF

  jstep0 = 0
  IF (is_restart_run() .AND. .NOT. time_config%is_relative_time) THEN
    ! get start counter for time loop from restart file:
    CALL get_restart_attribute("jstep", jstep0)
  END IF

  TIME_LOOP: DO jstep = (jstep0+1), (jstep0+nsteps)

    CALL add_time(dtime,0,0,0,datetime)

    ! read boundary data if necessary
    IF (l_limited_area .AND. (latbc_config%itype_latbc > 0)) &
      CALL read_latbc_data(p_patch(1), p_nh_state(1), p_int_state(1), ext_data(1), datetime)

    IF (msg_level >= 2 .OR. (jstep == (jstep0+1)) .OR. MOD(jstep,100) == 0) THEN
      WRITE(message_text,'(a,i10)') 'TIME STEP n: ', jstep
      CALL message(TRIM(routine),message_text)
    ENDIF

    IF ( check_newday(datetime_old,datetime) ) THEN

      WRITE(message_text,'(a,i10,a,i10)') 'New day  day_old: ', datetime_old%day, &
                &                 'day: ', datetime%day
      CALL message(TRIM(routine),message_text)

      !Update ndvi normalized differential vegetation index
      IF (itopo == 1 .AND. iforcing == inwp .AND.                  &
        & ALL(atm_phy_nwp_config(1:n_dom)%inwp_surface >= 1)) THEN
        DO jg=1, n_dom
          CALL interpol_monthly_mean(p_patch(jg), datetime,          &! in
            &                        ext_data(jg)%atm_td%ndvi_mrat,  &! in
            &                        ext_data(jg)%atm%ndviratio      )! out
        ENDDO

        ! after updating ndvi_mrat, probably plcov_t and tai_t have to be updated also.
        ! So it is better not to update ndvi_mrat till this is clarified 
        CALL update_ndvi(p_patch(1:), ext_data)
      END IF

      !Check if the the SST and Sea ice fraction have to be updated (sstice_mode 2,3,4)
      IF (sstice_mode > 1 .AND. iforcing == inwp  ) THEN

        CALL set_actual_td_ext_data (.FALSE., datetime,datetime_old,sstice_mode,  &
                                  &  p_patch(1:), ext_data, p_lnd_state)

        CALL update_sstice( p_patch(1:),           &
                        & ext_data, p_lnd_state, p_nh_state )

      END IF  !sstice_mode>1


      ! Check if MODIS albedo needs to be updated
      IF ( albedo_type == MODIS) THEN
        ! Note that here only an update of the external parameter fields is 
        ! performed. The actual update happens in mo_albedo.
        DO jg = 1, n_dom
          CALL interpol_monthly_mean(p_patch(jg), datetime,            &! in
            &                        ext_data(jg)%atm_td%alb_dif,      &! in
            &                        ext_data(jg)%atm%alb_dif          )! out

          CALL interpol_monthly_mean(p_patch(jg), datetime,            &! in
            &                        ext_data(jg)%atm_td%albuv_dif,    &! in
            &                        ext_data(jg)%atm%albuv_dif        )! out

          CALL interpol_monthly_mean(p_patch(jg), datetime,            &! in
            &                        ext_data(jg)%atm_td%albni_dif,    &! in
            &                        ext_data(jg)%atm%albni_dif        )! out
        ENDDO
      ENDIF

      datetime_old = datetime

    END IF !newday 
! end SST and sea ice fraction update

    ! Print control output for maximum horizontal and vertical wind speed
    IF (msg_level >= 5 .AND. MOD(jstep,iadv_rcf) == 1 .OR. msg_level >= 8) THEN 
      CALL print_maxwinds(p_patch(1), p_nh_state(1)%prog(nnow(1))%vn, p_nh_state(1)%prog(nnow(1))%w)
    ENDIF

    ! Store first old exner pressure
    ! (to prepare some kind of divergence damping, or to account for
    ! physically based 'implicit weights' in forward backward time stepping)
    IF (jstep == (jstep0+1) .AND. .NOT. is_restart_run()) THEN
!$OMP PARALLEL PRIVATE(jg)
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

    l_nml_output   = output_mode%l_nml   .AND. &
      &              ((jstep==(nsteps+jstep0)) .OR. &
      &              ((MOD(jstep_adv(1)%ntsteps+1,iadv_rcf)==0  .AND.  &
      &                  istime4name_list_output(jstep))))
    ! "l_outputtime", "l_diagtime": global flags used by other subroutines
    l_outputtime   = l_nml_output
    l_diagtime     = (.NOT. output_mode%l_none) .AND. &
      & ((jstep == (jstep0+1)) .OR. (MOD(jstep,n_diag) == 0) .OR. (jstep==(nsteps+jstep0)))

    ! Computation of diagnostic quantities may also be necessary for
    ! meteogram sampling:
    l_compute_diagnostic_quants = l_nml_output
    DO jg = 1, n_dom
      l_compute_diagnostic_quants = l_compute_diagnostic_quants .OR. &
        &          meteogram_is_sample_step(meteogram_output_config(jg), jstep,&
        &          jstep_adv(1)%ntsteps, iadv_rcf)
    END DO
    l_compute_diagnostic_quants = l_compute_diagnostic_quants .AND. &
      &                           .NOT. output_mode%l_none
    
    ! This serves to ensure that postprocessing is executed both immediately after
    ! activating and immediately after terminating a model domain
    ldom_active(1:n_dom) = p_patch(1:n_dom)%ldom_active

    ! This serves for enhancing the sound wave damping during the first 2 hours of integration
    ! If mixed second-order - fourth-order divergence damping is chosen (divdamp_order=24),
    ! the coefficient for second-order divergence damping is also updated
    elapsed_time_global = REAL(jstep,wp)*dtime
    IF (elapsed_time_global <= 7200._wp+REAL(iadv_rcf,wp)*dtime .AND. MOD(jstep,iadv_rcf) == 1  &
       .AND. .NOT. is_restart_run() .AND. .NOT. ltestcase) THEN
      CALL update_w_offctr(elapsed_time_global)
    ENDIF

    !--------------------------------------------------------------------------
    !
    ! dynamics stepping
    !
    CALL integrate_nh(datetime, 1, jstep, dtime, dtime_adv, 1)

    ldom_active(1:n_dom) = ldom_active(1:n_dom) .OR. p_patch(1:n_dom)%ldom_active

    ! Compute diagnostics for output if necessary
    IF (l_compute_diagnostic_quants) THEN
      CALL diag_for_output_dyn ( linit=.FALSE. )
      IF (iforcing == inwp) CALL diag_for_output_phys

      ! Unit conversion for output from mass mixing ratios to densities
      !
      DO jg = 1, n_dom
        CALL art_tools_interface('unit_conversion',p_nh_state(jg),jg)
      END DO
    ENDIF

    
    !--------------------------------------------------------------------------
    ! loop over the list of internal post-processing tasks, e.g.
    ! interpolate selected fields to p- and/or z-levels
    simulation_status = new_simulation_status(l_output_step  = l_nml_output,             &
      &                                       l_last_step    = (jstep==(nsteps+jstep0)), &
      &                                       l_dom_active   = ldom_active,     &
      &                                       i_timelevel    = nnow)
    CALL pp_scheduler_process(simulation_status)

#ifdef MESSY
    CALL messy_write_output
#endif

    ! output of results
    ! note: nnew has been replaced by nnow here because the update
    IF (l_nml_output) THEN
      CALL write_name_list_output(jstep)
    ENDIF

    ! sample meteogram output
    DO jg = 1, n_dom
      IF (.NOT. output_mode%l_none .AND. &    ! meteogram output is not initialized for output=none
        & meteogram_is_sample_step(meteogram_output_config(jg), jstep, jstep_adv(1)%ntsteps, iadv_rcf)) THEN
        CALL meteogram_sample_vars(jg, jstep, datetime, ierr)
        IF (ierr /= SUCCESS) THEN
          CALL finish (routine, 'Error in meteogram sampling! Sampling buffer too small?')
        ENDIF
      END IF
    END DO

    ! Diagnostics computation is not yet properly MPI-parallelized

    ! calls to "supervise_total_integrals_nh":
    ! - in the first time step (or the first time step after restart), or
    ! - if (MOD(jstep,n_diag) == 0), or
    ! - in the very last time step (jstep==nsteps)

    l_supervise_total_integrals =((lstep_adv(1) .AND. (jstep <= iadv_rcf)) .OR. &
      &                           (MOD(jstep,n_diag) == 0)                 .OR. &
      &                           (jstep==(nsteps+jstep0)))                .AND.&
      &                           (output_mode%l_totint)
    kstep = jstep-jstep0
    IF (jstep <= iadv_rcf)  kstep=1     !DR: necessary to work properly in combination with restart

    IF (l_supervise_total_integrals) THEN
#ifdef NOMPI
      IF (my_process_is_mpi_all_seq()) &
#endif
        CALL supervise_total_integrals_nh( kstep, p_patch(1:), p_nh_state,  &
        &                                  nnow(1:n_dom), nnow_rcf(1:n_dom), jstep == (nsteps+jstep0))
    ENDIF

    !--------------------------------------------------------------------------
    ! Write restart file
    !--------------------------------------------------------------------------
    IF (is_checkpoint_time(jstep,n_checkpoint) .AND. .NOT. output_mode%l_none) THEN
      lwrite_checkpoint = .TRUE.
    ENDIF

    ! Enforce that checkpointing files are written at the end of a physics time step
    ! (no reproducibility otherwise)
    IF (lwrite_checkpoint .AND. MOD(jstep_adv(1)%ntsteps,iadv_rcf)==0) THEN
      IF (use_async_restart_output) THEN
        DO jg = 1, n_dom
          CALL set_data_async_restart(p_patch(jg)%id, p_patch(jg)%ldom_active, &
            & opt_t_elapsed_phy          = t_elapsed_phy(jg,:),        &
            & opt_lcall_phy              = lcall_phy(jg,:),            &
            & opt_sim_time               = time_config%sim_time(jg),   &
            & opt_jstep_adv_ntstep       = jstep_adv(jg)%ntsteps,      &
            & opt_jstep_adv_marchuk_order= jstep_adv(jg)%marchuk_order,&
            & opt_depth_lnd              = nsoil(jg),                  &
            & opt_nlev_snow              = nsnow)
        ENDDO
        CALL write_async_restart(datetime, jstep)
      ELSE
        DO jg = 1, n_dom
          IF (.NOT. p_patch(jg)%ldom_active) CYCLE
          CALL create_restart_file( patch= p_patch(jg),datetime= datetime,                   &
                                  & jstep                      = jstep,                      &
                                  & opt_t_elapsed_phy          = t_elapsed_phy,              &
                                  & opt_lcall_phy              = lcall_phy,                  &
                                  & opt_sim_time               = time_config%sim_time(jg),   &
                                  & opt_jstep_adv_ntsteps      = jstep_adv(jg)%ntsteps,      &
                                  & opt_jstep_adv_marchuk_order= jstep_adv(jg)%marchuk_order,&
                                  & opt_depth_lnd              = nsoil(jg),                  &
                                  & opt_nlev_snow              = nsnow )
        END DO

        ! Create the master (meta) file in ASCII format which contains
        ! info about which files should be read in for a restart run.
        CALL write_restart_info_file
#ifdef MESSY
        CALL messy_channel_write_output(IOMODE_RST)
!        CALL messy_ncregrid_write_restart
#endif
      END IF

      lwrite_checkpoint = .FALSE.
    END IF

#ifdef MESSYTIMER
    ! timer sync
    CALL messy_timer_reset_time
#endif
  ENDDO TIME_LOOP

  IF (use_async_restart_output) CALL close_async_restart

  IF (ltimer) CALL timer_stop(timer_total)

  END SUBROUTINE perform_nh_timeloop
  !-------------------------------------------------------------------------

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
  RECURSIVE SUBROUTINE integrate_nh (datetime, jg, nstep_global,   &
    &                                dt_loc, dtadv_loc, num_steps )

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = 'mo_nh_stepping:integrate_nh'

    TYPE(t_datetime), INTENT(INOUT)         :: datetime

    INTEGER , INTENT(IN)    :: jg           !< current grid level
    INTEGER , INTENT(IN)    :: nstep_global !< counter of global time step
    INTEGER , INTENT(IN)    :: num_steps    !< number of time steps to be executed
    REAL(wp), INTENT(IN)    :: dt_loc       !< time step applicable to local grid level
    REAL(wp), INTENT(IN)    :: dtadv_loc    !< advective time step applicable to 
                                            !< local grid level

    ! Local variables

    ! Time levels
    INTEGER :: n_now_grf, n_now, n_new, n_save, n_temp
    INTEGER :: n_now_rcf, n_new_rcf         ! accounts for reduced calling frequencies (rcf)

    INTEGER :: jstep, jgp, jgc, jn
    INTEGER :: nsteps_nest ! number of time steps executed in nested domain

    REAL(wp):: dt_sub, dtadv_sub ! (advective) timestep for next finer grid level
    REAL(wp):: rdt_loc,  rdtadv_loc, rdtmflx_loc ! inverse time step for local grid level

    LOGICAL :: l_bdy_nudge
    INTEGER :: idyn_timestep
    LOGICAL :: l_recompute, lsave_mflx, lprep_adv, lfull_comp
    LOGICAL :: lclean_mflx   ! for reduced calling freqency: determines whether
                             ! mass-fluxes and trajectory-velocities are reset to zero
                             ! i.e. for starting new integration sweep
    LOGICAL :: lcall_hdiff
    LOGICAL :: lnest_active

    ! Switch to determine manner of OpenMP parallelization in interpol_scal_grf
!     LOGICAL :: lpar_fields=.FALSE.

    ! Switch to determine if nested domains are called at a given time step
    LOGICAL :: l_call_nests = .FALSE.

    !--------------------------------------------------------------------------
    ! This timer must not be called in nested domain because the model crashes otherwise
    IF (jg == 1 .AND. ltimer) CALL timer_start(timer_integrate_nh)
    
    !--------------------------------------------------------------------------
    ! settings for calling frequency for slow physics
    !--------------------------------------------------------------------------
    IF (jg == 1 .AND. linit_dyn(jg)) THEN

      IF (.NOT. l_nest_rcf) nsav1(1:n_dom) = nnow(1:n_dom)
    ENDIF
    !--------------------------------------------------------------------------
      
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
    ! If time-dependent data from a driving model are provided,
    ! they should be written to the save time level, so that the relaxation routine
    ! automatically does the right thing

    IF (jg == 1 .AND. l_limited_area .AND. linit_dyn(jg)) THEN

      n_save = nsav2(jg)
      n_now = nnow(jg)
!$OMP PARALLEL
!$OMP WORKSHARE
      p_nh_state(jg)%prog(n_save)%vn      = p_nh_state(jg)%prog(n_now)%vn
      p_nh_state(jg)%prog(n_save)%w       = p_nh_state(jg)%prog(n_now)%w
      p_nh_state(jg)%prog(n_save)%rho     = p_nh_state(jg)%prog(n_now)%rho
      p_nh_state(jg)%prog(n_save)%theta_v = p_nh_state(jg)%prog(n_now)%theta_v
!$OMP END WORKSHARE
!$OMP END PARALLEL
        
    ENDIF

    ! This executes one time step for the global domain and two steps for nested domains
    DO jstep = 1, num_steps

      ! Print control output for maximum horizontal and vertical wind speed
      ! (remark: the output for the global domain is called from perform_nh_timeloop)
      IF (jg > 1 .AND. msg_level >= 12) THEN 
        CALL print_maxwinds(p_patch(jg), p_nh_state(jg)%prog(nnow(jg))%vn, &
          p_nh_state(jg)%prog(nnow(jg))%w)
      ENDIF

      IF (ifeedback_type == 1 .AND. (jstep == 1) .AND. jg > 1 ) THEN
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

#ifdef MESSY
      CALL messy_global_start
#endif      
      !
      ! counter for simulation time in seconds
      time_config%sim_time(jg) = time_config%sim_time(jg) + dt_loc

      IF (itime_scheme == 1) THEN
        !------------------
        ! Pure advection
        !------------------

        lfull_comp = .TRUE.  ! full set of computations in prepare_tracer required

        SELECT CASE ( TRIM(nh_test_name) )

        CASE ('PA') ! solid body rotation

          ! set time-variant vertical velocity
          CALL set_nh_w_rho( p_patch(jg),p_nh_state(jg)%metrics,                    &! in
            & jstep_adv(jg)%marchuk_order, dt_loc, time_config%sim_time(jg)-dt_loc, &! in
            &               p_nh_state(jg)%prog(n_new)%w,                           &! inout
            &               p_nh_state(jg)%diag%pres,                               &! inout
            &               p_nh_state(jg)%diag%rho_ic                              )! inout

        CASE ('DF1', 'DF2', 'DF3', 'DF4') ! deformational flow

          ! get velocity field
          CALL get_nh_df_velocity( p_patch(jg), p_nh_state(jg)%prog(n_new), &
            &                     nh_test_name, rotate_axis_deg,            &
            &                     time_config%sim_time(jg)-dt_loc+dtadv_loc )


          ! get mass flux and new \rho. The latter one is only computed,
          ! if the density equation is re-integrated.
          CALL integrate_density_pa(p_patch(jg), p_int_state(jg),  & !in
            &                     p_nh_state(jg)%prog(n_now),      & !in
            &                     p_nh_state(jg)%prog(n_new),      & !in
            &                     p_nh_state(jg)%metrics,          & !in
            &                     p_nh_state(jg)%diag, dtadv_loc,  & !inout,in
            &                     jstep_adv(jg)%marchuk_order,     & !in
            &                     lcoupled_rho                     )
        END SELECT


        ! Diagnose some velocity-related quantities for the tracer
        ! transport scheme
        CALL prepare_tracer( p_patch(jg), p_nh_state(jg)%prog(n_now),     &! in
          &         p_nh_state(jg)%prog(n_new),                           &! in
          &         p_nh_state(jg)%metrics, p_int_state(jg),              &! in
          &         iadv_rcf, lstep_adv(jg), lclean_mflx, lfull_comp,     &! in
          &         p_nh_state(jg)%diag,                                  &! inout
          &         prep_adv(jg)%vn_traj, prep_adv(jg)%mass_flx_me,       &! inout
          &         prep_adv(jg)%w_traj, prep_adv(jg)%mass_flx_ic,        &! inout
          &         prep_adv(jg)%rhodz_mc_now, prep_adv(jg)%rhodz_mc_new, &! inout
          &         prep_adv(jg)%topflx_tra                               )! out


        IF (lstep_adv(jg)) THEN
          CALL step_advection( p_patch(jg), p_int_state(jg), dtadv_loc,    & !in
            &        jstep_adv(jg)%marchuk_order,                          & !in
            &        p_nh_state(jg)%prog(n_now_rcf)%tracer,                & !in
            &        prep_adv(jg)%mass_flx_me, prep_adv(jg)%vn_traj,       & !in
            &        prep_adv(jg)%mass_flx_ic, prep_adv(jg)%w_traj,        & !in
            &        p_nh_state(jg)%metrics%ddqz_z_full,                   & !in
            &        prep_adv(jg)%rhodz_mc_new, prep_adv(jg)%rhodz_mc_now, & !in
            &        p_nh_state(jg)%diag%grf_tend_tracer,                  & !inout
            &        p_nh_state(jg)%prog(n_new_rcf)%tracer,                & !inout
            &        p_nh_state(jg)%diag%hfl_tracer,                       & !out
            &        p_nh_state(jg)%diag%vfl_tracer,                       & !out
            &        opt_topflx_tra=prep_adv(jg)%topflx_tra,               & !in
            &        opt_q_int=p_nh_state(jg)%diag%q_int,                  & !out
            &        opt_ddt_tracer_adv=p_nh_state(jg)%diag%ddt_tracer_adv ) !out
        ENDIF


      ELSE  ! itime_scheme /= 1


        ! artificial forcing (Held-Suarez test forcing)
        IF ( lforcing .AND. iforcing == 1) THEN
          CALL held_suarez_nh_interface (p_nh_state(jg)%prog(n_now), p_patch(jg), &
                                         p_int_state(jg),p_nh_state(jg)%metrics,  &
                                         p_nh_state(jg)%diag)
        ENDIF

        ! Determine which physics packages must be called/not called at the current
        ! time step
        IF ( iforcing == inwp ) THEN
          CALL time_ctrl_physics ( dt_phy, lstep_adv, dt_loc, jg,  &! in
            &                      .FALSE.,                        &! in
            &                      t_elapsed_phy,                  &! inout
            &                      lcall_phy )                      ! out

          IF (msg_level >= 13) THEN
            WRITE(message_text,'(a,i2,a,5l2,a,6l2)') 'call phys. proc DOM:', &
              &  jg ,'   SP:', lcall_phy(jg,1:5), '   FP:',lcall_phy(jg,6:11)

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


        IF (jg > 1 .AND. .NOT. lfeedback(jg) .OR. jg == 1 .AND. l_limited_area) THEN
          ! apply boundary nudging if feedback is turned off and in limited-area mode
          l_bdy_nudge = .TRUE. 
        ELSE
          l_bdy_nudge = .FALSE.
        ENDIF

        ! index of dynamics time step within a large time step (ranges from 1 to iadv_rcf)
        idyn_timestep = MOD(jstep_adv(jg)%ntsteps,iadv_rcf)
        IF (idyn_timestep == 0) idyn_timestep = iadv_rcf

        IF (iforcing == inwp .AND. lclean_mflx) THEN
          l_recompute = .TRUE. ! always recompute velocity tendencies for predictor
        ELSE                   ! step after a physics call
          l_recompute = .FALSE.
        ENDIF

        IF (diffusion_config(jg)%lhdiff_vn .AND.                     &
          (.NOT. lhdiff_rcf .OR. lhdiff_rcf .AND. lstep_adv(jg))) THEN
          lcall_hdiff = .TRUE.
        ELSE
          lcall_hdiff = .FALSE.
        ENDIF

        IF (p_patch(jg)%n_childdom > 0 .AND. (jg == 1 .AND. MOD(nstep_global,iadv_rcf) == 1 &
            .OR. jg > 1 .AND. MOD(jstep,iadv_rcf) == 1) ) THEN
          lsave_mflx = .TRUE.
        ELSE
          lsave_mflx = .FALSE.
        ENDIF

        IF ( ltransport .OR. p_patch(jg)%n_childdom > 0 .AND. grf_intmethod_e >= 5) THEN
          lprep_adv = .TRUE.
        ELSE
          lprep_adv = .FALSE.
        ENDIF

        lfull_comp = .FALSE.  ! do not perform full set of computations in prepare_tracer

        ! For real-data runs, perform an extra diffusion call before the first time
        ! step because no other filtering of the interpolated velocity field is done
        IF (.NOT.ltestcase .AND. linit_dyn(jg) .AND. diffusion_config(jg)%lhdiff_vn .AND. &
            init_mode /= MODE_DWDANA) THEN
          CALL diffusion(p_nh_state(jg)%prog(n_now), p_nh_state(jg)%diag,       &
            p_nh_state(jg)%metrics, p_patch(jg), p_int_state(jg), dt_loc, .TRUE.)
        ENDIF

        IF (itype_comm == 1) THEN

          IF (ldynamics) THEN
            CALL solve_nh(p_nh_state(jg), p_patch(jg), p_int_state(jg), prep_adv(jg),        &
              n_now, n_new, linit_dyn(jg), l_recompute, lsave_mflx, lprep_adv, lstep_adv(jg),&
              lclean_mflx, idyn_timestep, jstep-1, l_bdy_nudge, dt_loc)
            
            IF (lcall_hdiff) &
              CALL diffusion(p_nh_state(jg)%prog(n_new), p_nh_state(jg)%diag,        &
                p_nh_state(jg)%metrics, p_patch(jg), p_int_state(jg), dt_loc, .FALSE.)
          ELSE
            CALL add_slowphys(p_nh_state(jg), p_patch(jg), n_now, n_new, dt_loc, &
                              lstep_adv(jg), n_now_rcf, n_new_rcf)
          ENDIF   
        ELSE
          CALL finish ( 'mo_nh_stepping', 'itype_comm /= 1 currently not implemented')
        ENDIF



        ! 5. tracer advection
        !-----------------------
        IF ( ltransport) THEN

          ! Diagnose some velocity-related quantities for the tracer
          ! transport scheme
          CALL prepare_tracer( p_patch(jg), p_nh_state(jg)%prog(n_now),     &! in
            &         p_nh_state(jg)%prog(n_new),                           &! in
            &         p_nh_state(jg)%metrics, p_int_state(jg),              &! in
            &         iadv_rcf, lstep_adv(jg), lclean_mflx, lfull_comp,     &! in
            &         p_nh_state(jg)%diag,                                  &! inout
            &         prep_adv(jg)%vn_traj, prep_adv(jg)%mass_flx_me,       &! inout
            &         prep_adv(jg)%w_traj,  prep_adv(jg)%mass_flx_ic,       &! inout
            &         prep_adv(jg)%rhodz_mc_now, prep_adv(jg)%rhodz_mc_new, &! inout
            &         prep_adv(jg)%topflx_tra                               )! out

        ENDIF

        IF ( ltransport ) THEN
          IF (lstep_adv(jg)) THEN

            IF (art_config(jg)%lart) THEN
              CALL art_emission_interface( p_patch(jg),          &!in
                &      dtadv_loc,                                &!in
                &      datetime,                                 &!in   
                &      p_nh_state(jg)%prog_list(n_now_rcf),      &!in
                &      p_nh_state(jg)%prog(n_new)%rho,           &!in 
                &      p_nh_state(jg)%prog(n_now_rcf)%tracer)     !inout
            ENDIF   


            CALL step_advection( p_patch(jg), p_int_state(jg), dtadv_loc,      & !in
              &          jstep_adv(jg)%marchuk_order,                          & !in
              &          p_nh_state(jg)%prog(n_now_rcf)%tracer,                & !in
              &          prep_adv(jg)%mass_flx_me, prep_adv(jg)%vn_traj,       & !in
              &          prep_adv(jg)%mass_flx_ic, prep_adv(jg)%w_traj,        & !in
              &          p_nh_state(jg)%metrics%ddqz_z_full,                   & !in
              &          prep_adv(jg)%rhodz_mc_new, prep_adv(jg)%rhodz_mc_now, & !in
              &          p_nh_state(jg)%diag%grf_tend_tracer,                  & !inout
              &          p_nh_state(jg)%prog(n_new_rcf)%tracer,                & !inout
              &          p_nh_state(jg)%diag%hfl_tracer,                       & !out
              &          p_nh_state(jg)%diag%vfl_tracer,                       & !out
              &          opt_topflx_tra=prep_adv(jg)%topflx_tra,               & !in
              &          opt_q_int=p_nh_state(jg)%diag%q_int,                  & !out
              &          opt_ddt_tracer_adv=p_nh_state(jg)%diag%ddt_tracer_adv ) !out

            IF (art_config(jg)%lart) THEN
!              CALL art_sedi_interface( p_patch(jg),             &!in
!                 &      dtadv_loc,                              &!in
!                 &      p_nh_state(jg)%prog_list(n_new_rcf),    &!in
!                 &      p_nh_state(jg)%metrics,                 &!in
!                 &      p_nh_state(jg)%prog(n_new)%rho,         &!in
!                 &      p_nh_state(jg)%diag,                    &!in
!                 &      p_nh_state(jg)%prog(n_new_rcf)%tracer,  &!inout
!                 &      p_nh_state(jg)%metrics%ddqz_z_full,     &!in
!                 &      prep_adv(jg)%rhodz_mc_new,              &!in
!                 &      opt_topflx_tra=prep_adv(jg)%topflx_tra)  !in
            ENDIF
                 
!            IF (  iforcing==inwp .AND. inwp_turb == icosmo) THEN
!              !> KF preliminary relabeling of TKE as long as there is no advection for it
!              p_nh_state(jg)%prog(n_new_rcf)%tke =  p_nh_state(jg)%prog(n_now)%tke
!            ENDIF


           ENDIF  !lstep_adv

        ENDIF


        IF (  iforcing==inwp .AND. lstep_adv(jg) ) THEN
       
          !Call interface for LES physics
          IF(atm_phy_nwp_config(jg)%is_les_phy)THEN     

            CALL les_phy_interface(lcall_phy(jg,:), .FALSE.,         & !in
              &                  lredgrid_phys(jg),                  & !in
              &                  dt_loc,                             & !in
              &                  dtadv_loc,                          & !in
              &                  t_elapsed_phy(jg,:),                & !in
              &                  time_config%sim_time(jg),           & !in
              &                  datetime,                           & !in
              &                  p_patch(jg)  ,                      & !in
              &                  p_int_state(jg),                    & !in
              &                  p_nh_state(jg)%metrics ,            & !in
              &                  p_patch(jgp),                       & !in
              &                  ext_data(jg)           ,            & !in
              &                  p_nh_state(jg)%prog(n_new) ,        & !inout
              &                  p_nh_state(jg)%prog(n_now_rcf),     & !in for tke
              &                  p_nh_state(jg)%prog(n_new_rcf) ,    & !inout
              &                  p_nh_state(jg)%diag ,               & !inout
              &                  prm_diag  (jg),                     & !inout
              &                  prm_nwp_tend(jg),                   &
              &                  p_lnd_state(jg)%diag_lnd,           &
              &                  p_lnd_state(jg)%prog_lnd(n_now_rcf),& !inout
              &                  p_lnd_state(jg)%prog_lnd(n_new_rcf),& !inout
              &                  p_lnd_state(jg)%prog_wtr(n_now_rcf),& !inout
              &                  p_lnd_state(jg)%prog_wtr(n_new_rcf),& !inout
              &                  p_nh_state(jg)%prog_list(n_new_rcf) ) !in

          ELSE !operational nwp physics
        
            !> moist tracer update is now synchronized with advection and satad
            CALL nwp_nh_interface(lcall_phy(jg,:), .FALSE.,          & !in
              &                  lredgrid_phys(jg),                  & !in
              &                  dt_loc,                             & !in
              &                  dtadv_loc,                          & !in
              &                  t_elapsed_phy(jg,:),                & !in
              &                  time_config%sim_time(jg),           & !in
              &                  datetime,                           & !in
              &                  p_patch(jg)  ,                      & !in
              &                  p_int_state(jg),                    & !in
              &                  p_nh_state(jg)%metrics ,            & !in
              &                  p_patch(jgp),                       & !in
              &                  ext_data(jg)           ,            & !in
              &                  p_nh_state(jg)%prog(n_new) ,        & !inout
              &                  p_nh_state(jg)%prog(n_now_rcf),     & !in for tke
              &                  p_nh_state(jg)%prog(n_new_rcf) ,    & !inout
              &                  p_nh_state(jg)%diag ,               & !inout
              &                  prm_diag  (jg),                     & !inout
              &                  prm_nwp_tend(jg),                   &
              &                  p_lnd_state(jg)%diag_lnd,           &
              &                  p_lnd_state(jg)%prog_lnd(n_now_rcf),& !inout
              &                  p_lnd_state(jg)%prog_lnd(n_new_rcf),& !inout
              &                  p_lnd_state(jg)%prog_wtr(n_now_rcf),& !inout
              &                  p_lnd_state(jg)%prog_wtr(n_new_rcf),& !inout
              &                  p_nh_state(jg)%prog_list(n_new_rcf),& !in
              &                  prep_adv(jg)%rhodz_mc_now,          & !in
              &                  prep_adv(jg)%rhodz_mc_new           ) !in

          END IF

          ! Boundary interpolation of land state variables entering into radiation computation
          ! if a reduced grid is used in the child domain(s)
          IF (ltimer)            CALL timer_start(timer_nesting)
          IF (timers_level >= 2) CALL timer_start(timer_bdy_interp)
          DO jn = 1, p_patch(jg)%n_childdom

            jgc = p_patch(jg)%child_id(jn)

            ! Remark: ideally, we should check for lcall_phy(jgc,itrad) here, but
            ! it is not yet known at the time of the call whether one of the 
            ! two physics time steps in the upcoming nest call will be a radiation
            ! step. Calling interpol_rrg_grf from the nested domain would be better
            ! from a flow control perspective but does not work with processor splitting
            !
            IF (lredgrid_phys(jgc) .AND. atm_phy_nwp_config(jgc)%inwp_surface >= 1 &
              .AND. lcall_phy(jg,itrad) ) THEN

              CALL interpol_rrg_grf(jg, jgc, jn, nnew_rcf(jg))
            ENDIF
          ENDDO
          IF (timers_level >= 2) CALL timer_stop(timer_bdy_interp)
          IF (ltimer)            CALL timer_stop(timer_nesting)

        ENDIF !iforcing

        ! Apply boundary nudging in case of one-way nesting
        IF (jg > 1 .AND. lstep_adv(jg)) THEN
          IF (ltimer)            CALL timer_start(timer_nesting)
          IF (timers_level >= 2) CALL timer_start(timer_nudging)

          IF (lfeedback(jg) .AND. l_density_nudging) THEN
            CALL density_boundary_nudging(jg,nnew(jg),REAL(iadv_rcf,wp))
          ELSE IF (.NOT. lfeedback(jg)) THEN
            CALL nest_boundary_nudging(jg,nnew(jg),nnew_rcf(jg),REAL(iadv_rcf,wp))
          ENDIF

          IF (timers_level >= 2) CALL timer_stop(timer_nudging)
          IF (ltimer)            CALL timer_stop(timer_nesting)
        ENDIF

      ENDIF  ! itime_scheme

      ! Update nudging tendency fields for limited-area mode
      IF (jg == 1 .AND. l_limited_area .AND. lstep_adv(jg)) THEN

        IF (latbc_config%itype_latbc > 0) THEN ! use time-dependent boundary data

          ! update the coefficients for the linear interpolation
          CALL update_lin_interc(datetime)

          CALL prep_outer_bdy_nudging(p_patch(jg),p_nh_state(jg)%prog(n_new),p_nh_state(jg)%prog(n_new_rcf), &
            p_nh_state(jg)%metrics,p_nh_state(jg)%diag,p_latbc_old=p_latbc_data(last_latbc_tlev)%atm,        &
            p_latbc_new=p_latbc_data(read_latbc_tlev)%atm)

        ELSE ! constant lateral boundary data

          CALL prep_outer_bdy_nudging(p_patch(jg),p_nh_state(jg)%prog(n_new),p_nh_state(jg)%prog(n_new_rcf), &
            p_nh_state(jg)%metrics,p_nh_state(jg)%diag,p_latbc_const=p_nh_state(jg)%prog(nsav2(jg)))

        ENDIF

        ! Apply nudging at the lateral boundaries
        CALL outer_boundary_nudging (jg, n_new, n_new_rcf, REAL(iadv_rcf,wp))

      ENDIF

      ! If there are nested domains...
      IF (l_nest_rcf .AND. lstep_adv(jg))  THEN
        l_call_nests = .TRUE.
        rdt_loc = 1._wp/(dt_loc*REAL(iadv_rcf,wp))  ! = 1._wp/dtadv_loc
        n_now_grf    = nsav1(jg)
        nsteps_nest  = 2*iadv_rcf
      ELSE IF (.NOT. l_nest_rcf) THEN
        l_call_nests = .TRUE.
        rdt_loc = 1._wp/dt_loc
        n_now_grf    = n_now
        nsteps_nest  = nest_substeps
      ELSE
        l_call_nests = .FALSE.
      ENDIF

      ! Check if at least one of the nested domains is active
      IF (l_call_nests .AND. p_patch(jg)%n_childdom > 0) THEN
        lnest_active = .FALSE.
        DO jn = 1, p_patch(jg)%n_childdom
          jgc = p_patch(jg)%child_id(jn)
          IF (p_patch(jgc)%ldom_active) lnest_active = .TRUE.
        ENDDO
        IF (.NOT. lnest_active) l_call_nests = .FALSE.
      ENDIF

      IF (l_call_nests .AND. p_patch(jg)%n_childdom > 0) THEN

        dt_sub     = dt_loc/2._wp    ! dyn. time step on next refinement level
        dtadv_sub  = dtadv_loc/2._wp ! adv. time step on next refinement level
        rdtadv_loc = 1._wp/dtadv_loc
        rdtmflx_loc = 1._wp/(dt_loc*REAL(MAX(1,iadv_rcf-1),wp))

        IF (ltimer)            CALL timer_start(timer_nesting)
        IF (timers_level >= 2) CALL timer_start(timer_bdy_interp)

        ! Compute time tendencies for interpolation to refined mesh boundaries
        CALL compute_tendencies (jg,n_new,n_now_grf,n_new_rcf,n_now_rcf,     &
          &                      rdt_loc,rdtadv_loc,rdtmflx_loc,lstep_adv(jg))

        ! Loop over nested domains
        DO jn = 1, p_patch(jg)%n_childdom

          jgc = p_patch(jg)%child_id(jn)

          ! Interpolate tendencies to lateral boundaries of refined mesh (jgc)
          IF (p_patch(jgc)%ldom_active) THEN
            CALL boundary_interpolation(jg, jgc,                         &
              &  n_now_grf,nnow(jgc),n_now_rcf,nnow_rcf(jgc),lstep_adv(jg),&
              &  prep_adv(jg)%mass_flx_me,prep_adv(jgc)%mass_flx_me)
          ENDIF

        ENDDO
        IF (timers_level >= 2) CALL timer_stop(timer_bdy_interp)

        IF (timers_level >= 2) CALL timer_start(timer_nudging)
        ! prep_bdy_nudging can not be called using delayed requests!
        DO jn = 1, p_patch(jg)%n_childdom

          jgc = p_patch(jg)%child_id(jn)
          IF (.NOT. p_patch(jgc)%ldom_active) CYCLE
          ! If feedback is turned off for child domain, compute parent-child
          ! differences for boundary nudging
          ! *** prep_bdy_nudging adapted for reduced calling frequency of tracers ***
          IF (lfeedback(jgc) .AND. l_density_nudging) THEN
            CALL prep_rho_bdy_nudging(jg,jgc)
          ELSE IF (.NOT. lfeedback(jgc)) THEN
            CALL prep_bdy_nudging(jg,jgc,lstep_adv(jg))
          ENDIF
        ENDDO
        IF (timers_level >= 2) CALL timer_stop(timer_nudging)
        IF (ltimer)            CALL timer_stop(timer_nesting)

        DO jn = 1, p_patch(jg)%n_childdom

          jgc = p_patch(jg)%child_id(jn)
          IF (.NOT. p_patch(jgc)%ldom_active) CYCLE

          IF(p_patch(jgc)%n_patch_cells > 0) THEN
            IF(proc_split) CALL push_glob_comm(p_patch(jgc)%comm, p_patch(jgc)%proc0)
            ! Recursive call to process_grid_level for child grid level
            CALL integrate_nh( datetime, jgc, nstep_global, dt_sub, &
              dtadv_sub, nsteps_nest )
            IF(proc_split) CALL pop_glob_comm()
          ENDIF

        ENDDO

        IF (ltimer)            CALL timer_start(timer_nesting)
        IF (timers_level >= 2) CALL timer_start(timer_feedback)
        DO jn = 1, p_patch(jg)%n_childdom

          ! Call feedback to copy averaged prognostic variables from refined mesh back
          ! to the coarse mesh (i.e. from jgc to jg)
          jgc = p_patch(jg)%child_id(jn)
          IF (.NOT. p_patch(jgc)%ldom_active) CYCLE

          IF (lfeedback(jgc)) THEN
            IF (ifeedback_type == 1) THEN
              CALL feedback(p_patch, p_nh_state, p_int_state, p_grf_state, p_lnd_state, jgc, &
                            jg, lstep_adv(jg))
            ELSE
              CALL relax_feedback(  p_patch(n_dom_start:n_dom),          &
                & p_nh_state(1:n_dom), p_int_state(n_dom_start:n_dom),   &
                & p_grf_state(n_dom_start:n_dom), jgc, jg, lstep_adv(jg))
            ENDIF
            ! Note: the last argument of "feedback" ensures that tracer feedback is
            ! only done for those time steps in which transport and microphysics are called
          ENDIF
        ENDDO
        IF (timers_level >= 2) CALL timer_stop(timer_feedback)
        IF (ltimer)            CALL timer_stop(timer_nesting)

      ENDIF

      IF (test_mode <= 0) THEN ! ... normal execution of time stepping
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
      ENDIF

      ! Check if nested domains have to activated or deactivated
      IF (lstep_adv(jg) .AND. p_patch(jg)%n_childdom > 0) THEN

        ! Loop over nested domains
        DO jn = 1, p_patch(jg)%n_childdom
          jgc = p_patch(jg)%child_id(jn)

          IF (p_patch(jgc)%ldom_active .AND. time_config%sim_time(jg) >= end_time(jgc)) THEN
            p_patch(jgc)%ldom_active = .FALSE.
            WRITE(message_text,'(a,i2,a,f12.2)') 'domain ',jgc,' stopped at time ',time_config%sim_time(jg)
            CALL message('integrate_nh', TRIM(message_text))
          ENDIF

          IF (.NOT. p_patch(jgc)%ldom_active .AND. time_config%sim_time(jg) >= start_time(jgc) .AND. &
              time_config%sim_time(jg) < end_time(jgc)) THEN
            p_patch(jgc)%ldom_active = .TRUE.

            IF (  atm_phy_nwp_config(jgc)%inwp_surface == 1 ) THEN
              CALL aggregate_landvars(p_patch(jg), ext_data(jg),                &
                p_lnd_state(jg)%prog_lnd(nnow_rcf(jg)), p_lnd_state(jg)%diag_lnd)
            ENDIF

            CALL initialize_nest(jg, jgc)

            ! Apply hydrostatic adjustment, using downward integration
            CALL hydro_adjust_downward(p_patch(jgc), p_nh_state(jgc)%metrics,                     &
              p_nh_state(jgc)%prog(nnow(jgc))%rho, p_nh_state(jgc)%prog(nnow(jgc))%exner,         &
              p_nh_state(jgc)%prog(nnow(jgc))%theta_v )

            CALL init_nwp_phy( dtime                    ,&
              & p_patch(jgc)                            ,&
              & p_nh_state(jgc)%metrics                 ,&
              & p_nh_state(jgc)%prog(nnow(jgc))         ,&
              & p_nh_state(jgc)%diag                    ,&
              & prm_diag(jgc)                           ,&
              & prm_nwp_tend(jgc)                       ,&
              & p_lnd_state(jgc)%prog_lnd(nnow_rcf(jgc)),&
              & p_lnd_state(jgc)%prog_lnd(nnew_rcf(jgc)),&
              & p_lnd_state(jgc)%prog_wtr(nnow_rcf(jgc)),&
              & p_lnd_state(jgc)%prog_wtr(nnew_rcf(jgc)),&
              & p_lnd_state(jgc)%diag_lnd               ,&
              & ext_data(jgc)                           ,&
              & phy_params(jgc)                          )

            time_config%sim_time(jgc) = time_config%sim_time(jg)

            WRITE(message_text,'(a,i2,a,f12.2)') 'domain ',jgc,' started at time ',time_config%sim_time(jg)
            CALL message('integrate_nh', TRIM(message_text))

          ENDIF
        ENDDO
      ENDIF

    ENDDO
    
    IF (jg == 1 .AND. ltimer) CALL timer_stop(timer_integrate_nh)

#ifdef MESSY
    CALL messy_global_end
#endif

  END SUBROUTINE integrate_nh

  !-------------------------------------------------------------------------
  !>
  !! Driver routine for initial call of physics routines.
  !! Apart from the full set of slow physics parameterizations, also turbulent transfer is 
  !! called, in order to have proper transfer coefficients available at the initial time step.
  !!
  !! This had to be moved ahead of the initial output for the physics fields to be more complete
  !!
  !! @par Revision History
  !! Developed by Guenther Zaengl, DWD (2013-01-04)
  !!
  RECURSIVE SUBROUTINE init_slowphysics (datetime, jg, dt_loc, dtadv_loc, sim_time)

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = 'mo_nh_stepping:init_slowphysics'

    TYPE(t_datetime), INTENT(in)         :: datetime

    INTEGER , INTENT(IN)    :: jg           !< current grid level
    REAL(wp), INTENT(IN)    :: dt_loc       !< time step applicable to local grid level
    REAL(wp), INTENT(IN)    :: dtadv_loc    !< advective time step applicable to 
                                            !< local grid level
    REAL(wp), INTENT(INOUT) :: sim_time(n_dom) !< elapsed simulation time on each
                                               !< grid level

    ! Local variables

    ! Time levels
    INTEGER :: n_now,n_now_rcf

    INTEGER :: jgp, jgc, jn

    REAL(wp):: dt_sub, dtadv_sub ! (advective) timestep for next finer grid level

    ! Determine parent domain ID
    IF ( jg > 1) THEN
      jgp = p_patch(jg)%parent_id
    ELSE IF (n_dom_start == 0) THEN
      jgp = 0
    ELSE
      jgp = 1
    ENDIF

    ! Set local variables for time levels
    n_now  = nnow(jg)

    ! Set local variable for rcf-time levels
    n_now_rcf = nnow_rcf(jg)

    CALL time_ctrl_physics ( dt_phy, lstep_adv, dt_loc, jg,  &! in
      &                      .TRUE.,                         &! in
      &                      t_elapsed_phy,                  &! inout
      &                      lcall_phy )                      ! out

    IF (msg_level >= 7) THEN
      WRITE(message_text,'(a,i2)') 'initial call of (slow) physics, domain ', jg
      CALL message(TRIM(routine), TRIM(message_text))
    ENDIF

    IF(atm_phy_nwp_config(jg)%is_les_phy)THEN!LES physics

      CALL les_phy_interface(lcall_phy(jg,:), .TRUE.,          & !in
        &                  lredgrid_phys(jg),                  & !in
        &                  dt_loc,                             & !in
        &                  dtadv_loc,                          & !in
        &                  dt_phy(jg,:),                       & !in
        &                  time_config%sim_time(jg),           & !in
        &                  datetime,                           & !in
        &                  p_patch(jg)  ,                      & !in
        &                  p_int_state(jg),                    & !in
        &                  p_nh_state(jg)%metrics ,            & !in
        &                  p_patch(jgp),                       & !in
        &                  ext_data(jg)           ,            & !in
        &                  p_nh_state(jg)%prog(n_now) ,        & !inout
        &                  p_nh_state(jg)%prog(n_now_rcf) ,    & !inout
        &                  p_nh_state(jg)%prog(n_now_rcf) ,    & !inout
        &                  p_nh_state(jg)%diag,                & !inout
        &                  prm_diag  (jg),                     & !inout
        &                  prm_nwp_tend(jg)                ,   &
        &                  p_lnd_state(jg)%diag_lnd,           &
        &                  p_lnd_state(jg)%prog_lnd(n_now_rcf),& !inout
        &                  p_lnd_state(jg)%prog_lnd(n_now_rcf),& !inout
        &                  p_lnd_state(jg)%prog_wtr(n_now_rcf),& !inout
        &                  p_lnd_state(jg)%prog_wtr(n_now_rcf),& !inout
        &                  p_nh_state(jg)%prog_list(n_now_rcf) ) !in
  
    ELSE !operational nwp physics
  
      CALL nwp_nh_interface(lcall_phy(jg,:), .TRUE.,           & !in
        &                  lredgrid_phys(jg),                  & !in
        &                  dt_loc,                             & !in
        &                  dtadv_loc,                          & !in
        &                  dt_phy(jg,:),                       & !in
        &                  time_config%sim_time(jg),           & !in
        &                  datetime,                           & !in
        &                  p_patch(jg)  ,                      & !in
        &                  p_int_state(jg),                    & !in
        &                  p_nh_state(jg)%metrics ,            & !in
        &                  p_patch(jgp),                       & !in
        &                  ext_data(jg)           ,            & !in
        &                  p_nh_state(jg)%prog(n_now) ,        & !inout
        &                  p_nh_state(jg)%prog(n_now_rcf) ,    & !inout
        &                  p_nh_state(jg)%prog(n_now_rcf) ,    & !inout
        &                  p_nh_state(jg)%diag,                & !inout
        &                  prm_diag  (jg),                     & !inout
        &                  prm_nwp_tend(jg)                ,   &
        &                  p_lnd_state(jg)%diag_lnd,           &
        &                  p_lnd_state(jg)%prog_lnd(n_now_rcf),& !inout
        &                  p_lnd_state(jg)%prog_lnd(n_now_rcf),& !inout
        &                  p_lnd_state(jg)%prog_wtr(n_now_rcf),& !inout
        &                  p_lnd_state(jg)%prog_wtr(n_now_rcf),& !inout
        &                  p_nh_state(jg)%prog_list(n_now_rcf),& !in
        &                  prep_adv(jg)%rhodz_mc_now,          & !in
        &                  prep_adv(jg)%rhodz_mc_new           ) !in 

    END IF

    ! Boundary interpolation of land state variables entering into radiation computation
    ! if a reduced grid is used in the child domain(s)
    DO jn = 1, p_patch(jg)%n_childdom

      jgc = p_patch(jg)%child_id(jn)
      IF (.NOT. p_patch(jgc)%ldom_active) CYCLE

      IF ( lredgrid_phys(jgc) ) THEN
        CALL interpol_rrg_grf(jg, jgc, jn, nnow_rcf(jg))
      ENDIF
    ENDDO

    IF (p_patch(jg)%n_childdom > 0) THEN

      dt_sub     = dt_loc/2._wp    ! dyn. time step on next refinement level
      dtadv_sub  = dtadv_loc/2._wp ! adv. time step on next refinement level

      DO jn = 1, p_patch(jg)%n_childdom

        jgc = p_patch(jg)%child_id(jn)
        IF (.NOT. p_patch(jgc)%ldom_active) CYCLE

        IF(p_patch(jgc)%n_patch_cells > 0) THEN
          IF(proc_split) CALL push_glob_comm(p_patch(jgc)%comm, p_patch(jgc)%proc0)
          CALL init_slowphysics( datetime, jgc, dt_sub, dtadv_sub, sim_time)
          IF(proc_split) CALL pop_glob_comm()
        ENDIF

      ENDDO

    ENDIF

  END SUBROUTINE init_slowphysics

  !-------------------------------------------------------------------------
  !>
  !! Diagnostic computations for output - dynamics fields
  !!
  !! This routine encapsulates calls to diagnostic computations required at output
  !! times only
  !!
  !! @par Revision History
  !! Developed by Guenther Zaengl, DWD (2012-05-09)
  !!
  SUBROUTINE diag_for_output_dyn (linit)

    LOGICAL, INTENT(IN) :: linit ! switch for computing additional diagnostics for initial output

!!$    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
!!$      &  routine = 'mo_nh_stepping:diag_for_output_dyn'

    ! Local variables
    INTEGER :: jg, jgc, jn ! loop indices

    REAL(wp), DIMENSION(:,:,:), POINTER  :: p_vn   => NULL()

    IF (ltimer) CALL timer_start(timer_nh_diagnostics)

    DO jg = 1, n_dom

      IF(p_patch(jg)%n_patch_cells == 0) CYCLE
      IF (.NOT. p_patch(jg)%ldom_active) CYCLE

      p_vn  => p_nh_state(jg)%prog(nnow(jg))%vn

        
      CALL rbf_vec_interpol_cell(p_vn,p_patch(jg),p_int_state(jg),&
                                 p_nh_state(jg)%diag%u,p_nh_state(jg)%diag%v)

      !CALL div(p_vn, p_patch(jg), p_int_state(jg), p_nh_state(jg)%diag%div)
      CALL div_avg(p_vn, p_patch(jg), p_int_state(jg),p_int_state(jg)%c_bln_avg,&
                                                          p_nh_state(jg)%diag%div)

      IF (linit) THEN
        CALL rot_vertex (p_vn, p_patch(jg), p_int_state(jg), p_nh_state(jg)%diag%omega_z)
      ENDIF

      ! Diagnose relative vorticity on cells
      CALL verts2cells_scalar(p_nh_state(jg)%diag%omega_z, p_patch(jg), &
        p_int_state(jg)%verts_aw_cells, p_nh_state(jg)%diag%vor)


      CALL diagnose_pres_temp (p_nh_state(jg)%metrics, p_nh_state(jg)%prog(nnow(jg)), &
        &                      p_nh_state(jg)%prog(nnow_rcf(jg)),                     &
        &                      p_nh_state(jg)%diag,p_patch(jg),                       &
        &                      opt_calc_temp=.TRUE.,                                  &
        &                      opt_calc_pres=.TRUE.                                   )

    ENDDO ! jg-loop

    ! Fill boundaries of nested domains
    DO jg = n_dom, 1, -1

      IF(p_patch(jg)%n_patch_cells == 0 .OR. p_patch(jg)%n_childdom == 0) CYCLE
      IF (.NOT. p_patch(jg)%ldom_active) CYCLE

      CALL sync_patch_array_mult(SYNC_C, p_patch(jg), 3, p_nh_state(jg)%diag%u,      &
        p_nh_state(jg)%diag%v, p_nh_state(jg)%diag%div)


      DO jn = 1, p_patch(jg)%n_childdom
        jgc = p_patch(jg)%child_id(jn)
        IF (.NOT. p_patch(jgc)%ldom_active) CYCLE

        CALL interpol_scal_grf (p_patch(jg), p_patch(jgc), p_grf_state(jg)%p_dom(jn), 3, &
             p_nh_state(jg)%diag%u, p_nh_state(jgc)%diag%u, p_nh_state(jg)%diag%v,       &
             p_nh_state(jgc)%diag%v, p_nh_state(jg)%diag%div, p_nh_state(jgc)%diag%div   )

      ENDDO

    ENDDO ! jg-loop

    IF (ltimer) CALL timer_stop(timer_nh_diagnostics)

  END SUBROUTINE diag_for_output_dyn

  !-------------------------------------------------------------------------
  !>
  !! Diagnostic computations for output - physics fields
  !!
  !! This routine encapsulates calls to diagnostic computations required at output
  !! times only
  !!
  !! @par Revision History
  !! Developed by Guenther Zaengl, DWD (2013-01-04)
  !!
  SUBROUTINE diag_for_output_phys

!!$    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
!!$      &  routine = 'mo_nh_stepping:diag_for_output_phys'

    ! Local variables
    INTEGER :: jg, jgc, jn ! loop indices

    IF (ltimer) CALL timer_start(timer_nh_diagnostics)

    DO jg = 1, n_dom

      IF(p_patch(jg)%n_patch_cells == 0) CYCLE
      IF (.NOT. p_patch(jg)%ldom_active) CYCLE

      IF (  atm_phy_nwp_config(jg)%inwp_surface == 1 ) THEN
        CALL aggregate_landvars( p_patch(jg), ext_data(jg),                 &
             p_lnd_state(jg)%prog_lnd(nnow_rcf(jg)), p_lnd_state(jg)%diag_lnd)
      ENDIF

    ENDDO ! jg-loop

    ! Fill boundaries of nested domains
    DO jg = n_dom, 1, -1

      IF(p_patch(jg)%n_patch_cells == 0 .OR. p_patch(jg)%n_childdom == 0) CYCLE
      IF (.NOT. p_patch(jg)%ldom_active) CYCLE

      CALL sync_patch_array(SYNC_C, p_patch(jg), p_nh_state(jg)%prog(nnow_rcf(jg))%tke)

      DO jn = 1, p_patch(jg)%n_childdom
        jgc = p_patch(jg)%child_id(jn)
        IF (.NOT. p_patch(jgc)%ldom_active) CYCLE

        CALL interpol_phys_grf(jg, jgc, jn) 

        IF (lfeedback(jgc)) CALL feedback_phys_diag(jgc, jg)

        CALL interpol_scal_grf (p_patch(jg), p_patch(jgc), p_grf_state(jg)%p_dom(jn), 1, &
           p_nh_state(jg)%prog(nnow_rcf(jg))%tke, p_nh_state(jgc)%prog(nnow_rcf(jgc))%tke)

      ENDDO

    ENDDO ! jg-loop

    IF (ltimer) CALL timer_stop(timer_nh_diagnostics)

  END SUBROUTINE diag_for_output_phys

  !-------------------------------------------------------------------------
  !>
  !! Update of vertical wind offcentering
  !!
  !! This routine handles the increased sound-wave damping (by increasing the vertical
  !! wind offcentering) during the initial spinup phase
  !!
  !! @par Revision History
  !! Developed by Guenther Zaengl, DWD (2013-06-04)
  !!
  SUBROUTINE update_w_offctr(elapsed_time)

    REAL(wp), INTENT(IN) :: elapsed_time
    INTEGER :: jg
    REAL(wp) :: min_vwind_impl_wgt, time1, time2

    time1 = 1800._wp  ! enhanced damping during the first half hour of integration
    time2 = 7200._wp  ! linear decrease of enhanced damping until time2

    IF (elapsed_time <= time1) THEN ! apply slightly super-implicit weights
      min_vwind_impl_wgt = 1.1_wp
      IF (divdamp_order == 24) divdamp_fac_o2 = 8._wp*divdamp_fac
    ELSE IF (elapsed_time <= time2) THEN ! linearly decrease minimum weights to 0.5
      min_vwind_impl_wgt = 0.5_wp + 0.6_wp*(time2-elapsed_time)/(time2-time1)
      IF (divdamp_order == 24) divdamp_fac_o2 = 8._wp*divdamp_fac*(time2-elapsed_time)/(time2-time1)
    ELSE
      min_vwind_impl_wgt = 0.5_wp
      IF (divdamp_order == 24) divdamp_fac_o2 = 0._wp
    ENDIF

    DO jg = 1, n_dom
      p_nh_state(jg)%metrics%vwind_impl_wgt(:,:) = MAX(min_vwind_impl_wgt, &
        p_nh_state(jg)%metrics%vwind_impl_wgt_sv(:,:))
      p_nh_state(jg)%metrics%vwind_expl_wgt(:,:) = 1._wp-p_nh_state(jg)%metrics%vwind_impl_wgt(:,:)
    ENDDO

  END SUBROUTINE update_w_offctr

  !-------------------------------------------------------------------------
  !>
  !! Physics time control
  !!
  !! Time control for slow and fast physics. This function provides a 2D array
  !! of type LOGICAL. For each physical process there is one column of length
  !! n_dom. A physical process (iphys) on domain (jg) will be called, if
  !! lcall_phy(jg,iphys)=.TRUE.. Whether it is .TRUE. or .FALSE. depends
  !! on the current time and the prescribed calling period listed in
  !! dt_phy.
  !!
  !! @par Revision History
  !! Developed by Daniel Reinert, DWD (2010-09-23)
  !! Modification by Daniel Reinert (2011-12-14)
  !! - replaced zoo of fast physics time steps by a single time step, named 
  !!   dt_fastphy.
  !!
  SUBROUTINE time_ctrl_physics ( dt_phy, lstep_adv, dt_loc, jg,  &
    &                            linit, t_elapsed_phy, lcall_phy )

    REAL(wp), INTENT(IN)    ::   &      !< Field of calling-time interval (seconds) for
      &  dt_phy(:,:)                    !< each domain and physical process

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

    INTEGER :: ip, ips                  !< loop index


  !-------------------------------------------------------------------------

    ! special treatment for the first physics call prior to the very first
    ! dynamics step
    IF (linit) THEN

      ! Initialize lcall_phy with .false. Only slow physics will be set to .true. initially;
      ! turbulent transfer is called by specifying the initialization mode of the NWP interface
      lcall_phy(jg,:)  = .FALSE.

      ! slow physics
      IF ( atm_phy_nwp_config(jg)%lproc_on(itconv) ) lcall_phy(jg,itconv) = .TRUE.

      IF ( atm_phy_nwp_config(jg)%lproc_on(itccov) ) lcall_phy(jg,itccov) = .TRUE.

      IF ( atm_phy_nwp_config(jg)%lproc_on(itrad)  ) lcall_phy(jg,itrad)  = .TRUE.

      IF ( atm_phy_nwp_config(jg)%lproc_on(itradheat) ) lcall_phy(jg,itradheat) = .TRUE.

      IF ( atm_phy_nwp_config(jg)%lproc_on(itsso)  ) lcall_phy(jg,itsso)  = .TRUE.

      IF ( atm_phy_nwp_config(jg)%lproc_on(itgwd)  ) lcall_phy(jg,itgwd)  = .TRUE.
 
    ELSE
      !
      ! all physical processes are forced to run at a multiple of 
      ! the advective time step. Note that fast physics are treated 
      ! as a combined process in the case of t_elapsed_phy and dt_phy, 
      ! but treated individually in the case of lcall_phy.
      !
      DO ips = 1, iphysproc_short

        ! If a physics package has been called at previous timestep,
        ! reset the time counter.
        !
        IF ( ANY(lcall_phy(jg,PACK(iproclist,map_phyproc(1:iphysproc,ips)))) ) THEN
           t_elapsed_phy(jg,ips)  = 0._wp
        ENDIF

        ! update time counter
        !
        t_elapsed_phy(jg,ips) = t_elapsed_phy(jg,ips) + dt_loc  ! dynamics timestep !!


        IF ( .NOT. lstep_adv(jg) ) THEN
          lcall_phy(jg,PACK(iproclist,map_phyproc(1:iphysproc,ips))) = .FALSE.
        ELSE

          ! The purpose of the 0.99999999 factor is to prevent pathological cases of 
          ! truncation error accumulation; will become obsolete when changing to 
          ! integer arithmetics for time control
          IF( t_elapsed_phy(jg,ips) >= 0.99999999_wp*dt_phy(jg,ips) ) THEN
            lcall_phy(jg,PACK(iproclist,map_phyproc(1:iphysproc,ips)))  = .TRUE.
          ELSE
            lcall_phy(jg,PACK(iproclist,map_phyproc(1:iphysproc,ips)))  = .FALSE.
          ENDIF

        ENDIF

      ENDDO  ! ips

      ! In addition, it must be checked, whether the individual processes 
      ! are switched on at all (lproc_on =.TRUE.). If not, lcall_phy is 
      ! reset to false.
      !
      DO ip = 1, iphysproc   ! not that we have to loop over ALL processes
        IF (.NOT. atm_phy_nwp_config(jg)%lproc_on(ip) ) THEN
          lcall_phy(jg,ip)  = .FALSE.
        ENDIF
      ENDDO  ! ip

    ENDIF

  END SUBROUTINE time_ctrl_physics


  !-------------------------------------------------------------------------
  !>
  !! Setup of physics time control
  !!
  !! Setup of time control for slow and fast physics. The mapping matrix is 
  !! initialized, which provides mapping rules required when mapping 
  !! between variables of size iphysproc and iphysproc_short. Typical examples 
  !! are lcall_phy and t_elapsed_phy.
  !!
  !! @par Revision History
  !! Developed by Daniel Reinert, DWD (2011-12-14)
  !!
  SUBROUTINE setup_time_ctrl_physics ( )

  !-------------------------------------------------------------------------

  ! list of physical processes (x-axis of mapping matrix)
    iproclist = (/ itconv,itccov,itrad,itsso,itgwd,itupdate,itsatad,itturb,&
      &            itgscp,itsfc,itradheat /)

    map_phyproc(1:iphysproc,1:iphysproc_short) = .FALSE. ! initialization

    map_phyproc(1,1)    = .TRUE.  ! simple one to one mapping
    map_phyproc(2,2)    = .TRUE.  ! simple one to one mapping
    map_phyproc(3,3)    = .TRUE.  ! simple one to one mapping
    map_phyproc(4,4)    = .TRUE.  ! simple one to one mapping
    map_phyproc(5,5)    = .TRUE.  ! simple one to one mapping
    map_phyproc(6:11,6) = .TRUE.  ! mapping of fast physics processes to single one

  END SUBROUTINE setup_time_ctrl_physics
  !-------------------------------------------------------------------------
  

  !-------------------------------------------------------------------------
  !>
  !! @par Revision History
  !!
  SUBROUTINE deallocate_nh_stepping

  INTEGER                              ::  jg, ist

  !-----------------------------------------------------------------------
  !
  ! deallocate auxiliary fields for tracer transport and rcf
  !
  DO jg = 1, n_dom
    DEALLOCATE( prep_adv(jg)%mass_flx_me, prep_adv(jg)%mass_flx_ic,    &
      &         prep_adv(jg)%vn_traj, prep_adv(jg)%w_traj,             &
      &         prep_adv(jg)%rhodz_mc_now, prep_adv(jg)%rhodz_mc_new,  &
      &         prep_adv(jg)%topflx_tra, STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ( 'mo_nh_stepping: perform_nh_stepping',            &
        &    'deallocation for mass_flx_me, mass_flx_ic, vn_traj,' // &
        &    'w_traj, rhodz_mc_now, rhodz_mc_new, topflx_tra failed' )
    ENDIF
  ENDDO

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
  DEALLOCATE( lstep_adv, lcall_phy, linit_slowphy, linit_dyn, &
    &         t_elapsed_phy, STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ( 'mo_nh_stepping: perform_nh_stepping',          &
      &    'deallocation for lstep_adv, lcall_phy,' //            &
      &    't_elapsed_phy failed' )
  ENDIF

  IF (l_limited_area .AND. (latbc_config%itype_latbc > 0)) THEN
    CALL deallocate_latbc_data(p_patch(1))
  ENDIF

  END SUBROUTINE deallocate_nh_stepping
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! @par Revision History
  !!
  SUBROUTINE allocate_nh_stepping
!
  INTEGER                              :: jg, jp !, nlen
  INTEGER                              :: ist
  CHARACTER(len=MAX_CHAR_LENGTH)       :: attname   ! attribute name

!-----------------------------------------------------------------------

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
    &      linit_dyn(n_dom),t_elapsed_phy(n_dom,iphysproc_short),            &
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
      CALL get_restart_attribute(TRIM(attname), time_config%sim_time(jg))
      DO jp = 1,iphysproc_short
        WRITE(attname,'(a,i2.2,a,i2.2)') 't_elapsed_phy_DOM',jg,'_PHY',jp
        CALL get_restart_attribute(TRIM(attname), t_elapsed_phy(jg,jp))
      ENDDO
      DO jp = 1,iphysproc
        WRITE(attname,'(a,i2.2,a,i2.2)') 'lcall_phy_DOM',jg,'_PHY',jp
        CALL get_restart_attribute(TRIM(attname), lcall_phy(jg,jp))
      ENDDO
    ENDDO
    linit_dyn(:)      = .FALSE.
  ELSE
    jstep_adv(:)%ntsteps       = 0
    jstep_adv(:)%marchuk_order = 0
    time_config%sim_time(:)    = 0._wp
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
      &  prep_adv(jg)%topflx_tra  (nproma,p_patch(jg)%nblks_c,MAX(1,ntracer)),     &
      &       STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ( 'mo_nh_stepping: perform_nh_stepping',           &
      &      'allocation for mass_flx_me, mass_flx_ic, vn_traj, ' // &
      &      'w_traj, rhodz_mc_now, rhodz_mc_new, topflx_tra failed' )
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
    prep_adv(jg)%topflx_tra  (:,:,:) = 0._wp
!$OMP END WORKSHARE
!$OMP END PARALLEL

  ENDDO


  IF (l_limited_area .AND. (latbc_config%itype_latbc > 0)) &
    &   CALL prepare_latbc_data(p_patch(1), p_int_state(1), p_nh_state(1), ext_data(1))

END SUBROUTINE allocate_nh_stepping
  !-----------------------------------------------------------------------------

END MODULE mo_nh_stepping

