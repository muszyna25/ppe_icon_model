!>
!! @author <name, affiliation>
!!
!! @par Revision History
!! <Description of activity> by <name, affiliation> (<YYYY-MM-DD>)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_timer

#if (defined(__GFORTRAN__) || defined(_CRAYFTN) || defined(__PGIF90__))
  USE iso_fortran_env, ONLY: compiler_version, compiler_options
#endif
  
#ifdef __SCT__
  USE sct, ONLY: new_timer     => sct_new_timer,             &
       &         timer_start   => sct_start,                 &
       &         timer_stop    => sct_stop,                  &
       &         cleanup_timer => sct_reset_timer,           &
       &         delete_timer  => sct_del_timer,             &
       &         sct_init, sct_report, sct_add_report_attribute
#else
  USE mo_real_timer, ONLY: new_timer,                        &
       &                   timer_start,                      &
       &                   timer_stop,                       &
       &                   timer_report,                     &
       &                   cleanup_timer => timer_reset,     &
       &                   delete_timer => del_timer

#endif

  USE mo_run_config, ONLY: ltimer, timers_level,  activate_sync_timers, iforcing

  USE mo_impl_constants, ONLY: iecham
  
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ltimer, timers_level, activate_sync_timers
  PUBLIC :: new_timer, timer_start, timer_stop                !< procedures imported from mo_real_timer
  PUBLIC :: print_timer, cleanup_timer, delete_timer          !< procedures imported and renamed
  PUBLIC :: init_timer                                        !< procedure of this module

  PUBLIC :: timer_total                         !< IDs of timers
  PUBLIC :: timer_exch_data, timer_exch_data_rv, timer_exch_data_async, timer_exch_data_wait
  PUBLIC :: timer_global_sum, timer_omp_global_sum, timer_ordglb_sum, timer_omp_ordglb_sum
  PUBLIC :: timer_icon_comm_sync
  PUBLIC :: timer_icon_comm_fillrecv, timer_icon_comm_wait, timer_icon_comm_isend,       &
       &    timer_icon_comm_ircv, timer_icon_comm_fillsend, timer_icon_comm_fillandsend, &
       &    timer_icon_comm_barrier_2, timer_icon_comm_send
  PUBLIC :: timer_barrier

  PUBLIC :: timer_integrate_nh
  PUBLIC :: timer_solve_nh, timer_solve_nh_veltend, timer_solve_nh_cellcomp, timer_solve_nh_edgecomp, &
       &    timer_solve_nh_vnupd, timer_solve_nh_vimpl, timer_solve_nh_exch
  PUBLIC :: timer_physics
                        !< IDs of timers
  PUBLIC :: timer_radiaton_recv, timer_radiaton_comp, timer_radiaton_send, &
       &    timer_preradiaton, timer_synsat

  PUBLIC :: timer_div, timer_grad, timer_gmres, timer_lhs, timer_lhs_sp
  PUBLIC :: timer_corio, timer_intp
  PUBLIC :: timer_transport
  PUBLIC :: timer_back_traj
  PUBLIC :: timer_nh_hdiffusion

  PUBLIC :: timer_update_prog_phy
  PUBLIC :: timer_diagnose_pres_temp
  PUBLIC :: timer_nh_diagnostics

  ! atmosphere - ocean coupling
  PUBLIC :: timer_coupling, timer_coupling_init
  PUBLIC :: timer_coupling_1stget, timer_coupling_get
  PUBLIC :: timer_coupling_put

  ! iconam - echam coupling
  PUBLIC :: timer_iconam_echam
  PUBLIC :: timer_dyn2phy, timer_d2p_prep, timer_d2p_sync, timer_d2p_couple
  PUBLIC :: timer_echam_bcs, timer_echam_phy
  PUBLIC :: timer_phy2dyn, timer_p2d_prep, timer_p2d_sync, timer_p2d_couple
  !
  ! echam physics
  PUBLIC :: timer_cov
  PUBLIC :: timer_rad , timer_rht
  PUBLIC :: timer_vdf , timer_vdf_dn , timer_vdf_sf , timer_vdf_up
  PUBLIC :: timer_gwd !!$, timer_sso
  PUBLIC :: timer_cnv , timer_cld
  PUBLIC :: timer_car , timer_mox
  PUBLIC :: timer_wmo
  !
  ! echam radiation
  PUBLIC :: timer_rrtm_prep, timer_rrtm_post
  PUBLIC :: timer_lrtm, timer_srtm
#ifdef PSRAD_TIMING
  PUBLIC :: timer_gas_optics_lw, timer_gas_optics_sw, timer_cloud_optics, &
    timer_sample_cloud_lw, timer_sample_cloud_sw, &
    timer_rrtm_coeffs, timer_psrad_scaling, timer_psrad_aerosol
#endif

  PUBLIC :: timer_satad_v_3D
  PUBLIC :: timer_phys_exner
  PUBLIC :: timer_phys_u_v
  PUBLIC :: timer_nwp_turbulence, timer_nwp_surface
  PUBLIC :: timer_nwp_microphysics
  PUBLIC :: timer_phys_sync_patch
  PUBLIC :: timer_fast_phys
  PUBLIC :: timer_nwp_convection
  PUBLIC :: timer_nwp_radiation
  PUBLIC :: timer_pre_radiation_nwp
  PUBLIC :: timer_phys_acc, timer_phys_acc_1,timer_phys_acc_2, timer_phys_dpsdt
  PUBLIC :: timer_phys_sync_tracers
  PUBLIC :: timer_phys_sync_tempv
  PUBLIC :: timer_phys_acc_par
  PUBLIC :: timer_phys_sync_ddt_u
  PUBLIC :: timer_phys_sync_vn

  PUBLIC :: timer_held_suarez_intr

!   PUBLIC :: timer_sync_wait
!   PUBLIC :: timer_sync_delay,timer_sync_outbuffer
!   PUBLIC :: timer_sync_psend_1, timer_sync_isend_2, timer_sync_recv_2,timer_sync_isend_3

  PUBLIC :: timer_sso
  PUBLIC :: timer_cover_koe
  PUBLIC :: timer_radiation
  PUBLIC :: timer_radheat
  PUBLIC :: timer_omp_radiation
  PUBLIC :: timer_lonlat_setup
  PUBLIC :: timer_write_output
  PUBLIC :: timer_model_init, timer_init_latbc
  PUBLIC :: timer_domain_decomp, timer_compute_coeffs, timer_ext_data, timer_init_icon, timer_read_restart
  PUBLIC :: timer_solve_ab, timer_tracer_ab, timer_vert_veloc, timer_normal_veloc
  PUBLIC :: timer_upd_phys, timer_upd_flx, timer_calc_moc
  PUBLIC :: timer_ab_expl, timer_ab_rhs4sfc
  PUBLIC :: timer_adv_horz, timer_dif_horz, timer_hflx_lim
  PUBLIC :: timer_adv_vert, timer_dif_vert, timer_ppm_slim, timer_adpo_vert
  PUBLIC :: timer_dbg_prnt
  PUBLIC :: timer_si_correction
  PUBLIC :: timer_cube_root
  PUBLIC :: timer_RK_tend, timer_RK_update, timer_step_RK

  PUBLIC :: timer_intrp_diagn
  PUBLIC :: timer_step_2tl_si
  PUBLIC :: timer_prep_echam_phy
  PUBLIC :: timer_prep_phy
  PUBLIC :: timer_prep_tracer_leapfrog
  PUBLIC :: timer_prep_tracer
  PUBLIC :: timer_prep_tracer_RK
  PUBLIC :: timer_hdiff_expl
  PUBLIC :: timer_dyn_theta, timer_dyn_temp

  PUBLIC :: timer_con_l_theta2t, timer_con_l_t2theta, timer_con_theta2t, timer_con_t2theta

  PUBLIC :: timer_nesting
  PUBLIC :: timer_nudging
  PUBLIC :: timer_bdy_interp
  PUBLIC :: timer_feedback

  ! ocean
  PUBLIC :: timer_gmres_p_sum
  PUBLIC :: timer_scalar_prod_veloc

  ! Timer IDs for sea ice
  PUBLIC :: timer_ice_fast, timer_ice_slow, timer_ice_slow2, timer_ice_momentum, &
       &    timer_ice_interp, timer_ice_advection

  ! Timer IDs for HAMOCC

  PUBLIC :: timer_bgc_up_bgc, timer_bgc_swr, timer_bgc_wea, timer_bgc_depo, &
       &    timer_bgc_chemcon, timer_bgc_ocprod, timer_bgc_sett, timer_bgc_cya,&
       &    timer_bgc_gx, timer_bgc_calc, timer_bgc_powach, timer_bgc_up_ic, &
       &    timer_bgc_tend,timer_bgc_ini, timer_bgc_inv, timer_bgc_tot 

  ! restart timers
  PUBLIC :: timer_load_restart
  PUBLIC :: timer_load_restart_io
  PUBLIC :: timer_load_restart_comm_setup
  PUBLIC :: timer_load_restart_communication
  PUBLIC :: timer_load_restart_get_var_id
  PUBLIC :: timer_write_restart
  PUBLIC :: timer_write_restart_io
  PUBLIC :: timer_write_restart_communication
  PUBLIC :: timer_write_restart_wait
  PUBLIC :: timer_write_restart_setup
  PUBLIC :: timer_restart_collector_setup
  PUBLIC :: timer_restart_indices_setup

  ! Timer for data assimilation
  PUBLIC :: timer_datass, timer_lhn

  PUBLIC :: timer_extra1,  timer_extra2,  timer_extra3,  timer_extra4,  timer_extra5,  &
            timer_extra6,  timer_extra7,  timer_extra8,  timer_extra9,  timer_extra10, &
            timer_extra11, timer_extra12, timer_extra13, timer_extra14, timer_extra15, &
            timer_extra16, timer_extra17, timer_extra18, timer_extra19, timer_extra20, &
            timer_extra21, timer_extra22, timer_extra23, timer_extra24, timer_extra25, &
            timer_extra26, timer_extra27, timer_extra28, timer_extra29, timer_extra30, &
            timer_extra31, timer_extra32, timer_extra33, timer_extra34, timer_extra35, &
            timer_extra36, timer_extra37, timer_extra38, timer_extra39, timer_extra40

  ! ART timers of the interfaces
  PUBLIC :: timer_art, timer_art_emissInt, timer_art_reacInt,                                  &
            timer_art_cldInt, timer_art_diagInt, timer_art_sedInt, timer_art_toolInt,          &
            timer_art_tracInt, timer_art_turbdiffInt, timer_art_washoutInt, timer_art_initInt, &
            timer_art_radInt, timer_art_photo, timer_art_losschem

  ! low level timing routine
  PUBLIC :: tic, toc
  PUBLIC :: timer_ls_forcing 

  !-------------------
  ! Module variables
  !-------------------

  ! ID of timer for total model integration time
  INTEGER :: timer_total
  INTEGER :: timer_exch_data, timer_exch_data_rv, timer_exch_data_async, timer_exch_data_wait
  INTEGER :: timer_global_sum, timer_omp_global_sum, timer_ordglb_sum, timer_omp_ordglb_sum
  INTEGER :: timer_icon_comm_sync
  INTEGER :: timer_icon_comm_fillrecv, timer_icon_comm_wait, timer_icon_comm_isend, &
       &     timer_icon_comm_ircv, timer_icon_comm_fillsend,timer_icon_comm_fillandsend,   &
       &     timer_icon_comm_barrier_2, timer_icon_comm_send
  INTEGER :: timer_barrier
  INTEGER :: timer_gmres_p_sum

  INTEGER :: timer_nh_hdiffusion

  INTEGER :: timer_integrate_nh
  INTEGER :: timer_solve_nh, timer_solve_nh_veltend, timer_solve_nh_cellcomp, timer_solve_nh_edgecomp, &
       &     timer_solve_nh_vnupd, timer_solve_nh_vimpl, timer_solve_nh_exch
  INTEGER :: timer_physics
  INTEGER :: timer_update_prog_phy

  INTEGER :: timer_nh_diagnostics
  INTEGER :: timer_diagnose_pres_temp
  INTEGER :: timer_satad_v_3D
  INTEGER :: timer_phys_exner
  INTEGER :: timer_phys_u_v
  INTEGER :: timer_nwp_turbulence, timer_nwp_surface
  INTEGER :: timer_nwp_microphysics
  INTEGER :: timer_phys_sync_patch
  INTEGER :: timer_fast_phys
  INTEGER :: timer_nwp_convection
  INTEGER :: timer_nwp_radiation
  INTEGER :: timer_synsat
  INTEGER :: timer_radiaton_recv, timer_radiaton_comp, timer_radiaton_send, &
       &     timer_preradiaton
  INTEGER :: timer_pre_radiation_nwp
  INTEGER :: timer_phys_acc, timer_phys_acc_1,timer_phys_acc_2, timer_phys_dpsdt
  INTEGER :: timer_phys_sync_tracers
  INTEGER :: timer_phys_sync_tempv
  INTEGER :: timer_phys_acc_par
  INTEGER :: timer_phys_sync_ddt_u
  INTEGER :: timer_phys_sync_vn
  INTEGER :: timer_dyn_theta, timer_dyn_temp
!   INTEGER :: timer_sync_wait
!   INTEGER :: timer_sync_delay,timer_sync_outbuffer
!   INTEGER :: timer_sync_psend_1, timer_sync_isend_2, timer_sync_recv_2,timer_sync_isend_3

  INTEGER :: timer_sso
  INTEGER :: timer_cover_koe
  INTEGER :: timer_radiation
  INTEGER :: timer_radheat

  ! Timer ID's for horizontal operators
  INTEGER :: timer_div
  INTEGER :: timer_grad
  INTEGER :: timer_gmres
  INTEGER :: timer_lhs, timer_lhs_sp
  INTEGER :: timer_corio
  INTEGER :: timer_intp

  INTEGER :: timer_transport    ! tracer transport
  INTEGER :: timer_back_traj

  ! Timer ID's for forcings and testcases
  INTEGER :: timer_held_suarez_intr

  ! Timer ID's for ocean-atmosphere coupling
  INTEGER :: timer_coupling, timer_coupling_init
  INTEGER :: timer_coupling_get, timer_coupling_1stget, timer_coupling_put

  ! Timer ID's for physics-dynamics coupling

  ! iconam - echam coupling
  INTEGER :: timer_iconam_echam
  INTEGER :: timer_dyn2phy, timer_d2p_prep, timer_d2p_sync, timer_d2p_couple
  INTEGER :: timer_echam_bcs, timer_echam_phy
  INTEGER :: timer_phy2dyn, timer_p2d_prep, timer_p2d_sync, timer_p2d_couple
  !
  ! echam physics
  INTEGER :: timer_cov
  INTEGER :: timer_rad , timer_rht
  INTEGER :: timer_vdf , timer_vdf_dn , timer_vdf_sf , timer_vdf_up
  INTEGER :: timer_gwd !!$, timer_sso
  INTEGER :: timer_cnv , timer_cld
  INTEGER :: timer_car , timer_mox
  INTEGER :: timer_wmo
  !
  ! echam radiation
  INTEGER :: timer_rrtm_prep, timer_rrtm_post
  INTEGER :: timer_lrtm, timer_srtm
#ifdef PSRAD_TIMING
  INTEGER :: timer_gas_optics_lw, timer_gas_optics_sw, timer_cloud_optics, &
    timer_sample_cloud_lw, timer_sample_cloud_sw, &
    timer_rrtm_coeffs, timer_psrad_scaling, timer_psrad_aerosol
#endif

  INTEGER :: timer_omp_radiation
  INTEGER :: timer_write_output
  INTEGER :: timer_model_init, timer_init_latbc
  INTEGER :: timer_domain_decomp, timer_compute_coeffs, timer_ext_data, timer_init_icon, timer_read_restart
  INTEGER :: timer_solve_ab, timer_tracer_ab, timer_vert_veloc, timer_normal_veloc
  INTEGER :: timer_upd_phys, timer_upd_flx, timer_calc_moc
  INTEGER :: timer_ab_expl, timer_ab_rhs4sfc
  INTEGER :: timer_adv_horz, timer_dif_horz, timer_hflx_lim
  INTEGER :: timer_adv_vert, timer_dif_vert, timer_ppm_slim, timer_adpo_vert
  INTEGER :: timer_dbg_prnt
  INTEGER :: timer_si_correction
  INTEGER :: timer_cube_root
  INTEGER :: timer_RK_tend, timer_RK_update, timer_step_RK

  INTEGER :: timer_intrp_diagn
  INTEGER :: timer_step_2tl_si
  INTEGER :: timer_prep_echam_phy
  INTEGER :: timer_prep_phy
  INTEGER :: timer_prep_tracer_leapfrog
  INTEGER :: timer_prep_tracer
  INTEGER :: timer_prep_tracer_RK
  INTEGER :: timer_hdiff_expl

  ! Timer ID for optional lon-lat interpolation
  INTEGER :: timer_lonlat_setup

  ! Timer IDs for boundary interpolation, feedback & nudging
  INTEGER :: timer_nesting
  INTEGER :: timer_nudging
  INTEGER :: timer_bdy_interp
  INTEGER :: timer_feedback

  INTEGER :: timer_con_l_theta2t, timer_con_l_t2theta, timer_con_theta2t, timer_con_t2theta

  ! ocean
  INTEGER :: timer_scalar_prod_veloc
  ! Timer IDs for sea ice
  INTEGER :: timer_ice_fast, timer_ice_slow, timer_ice_slow2, timer_ice_momentum,       &
       &     timer_ice_interp, timer_ice_advection

  ! Timer IDs HAMOCC
  INTEGER :: timer_bgc_up_bgc, timer_bgc_swr, timer_bgc_wea, timer_bgc_depo, &
       &     timer_bgc_chemcon, timer_bgc_ocprod, timer_bgc_sett, timer_bgc_cya,&
       &     timer_bgc_gx, timer_bgc_calc, timer_bgc_powach, timer_bgc_up_ic, &
       &     timer_bgc_tend, timer_bgc_ini, timer_bgc_inv, timer_bgc_tot

  ! restart timers
  INTEGER :: timer_load_restart
  INTEGER :: timer_load_restart_io
  INTEGER :: timer_load_restart_comm_setup
  INTEGER :: timer_load_restart_communication
  INTEGER :: timer_load_restart_get_var_id
  INTEGER :: timer_write_restart
  INTEGER :: timer_write_restart_io
  INTEGER :: timer_write_restart_communication
  INTEGER :: timer_write_restart_wait
  INTEGER :: timer_write_restart_setup
  INTEGER :: timer_restart_collector_setup
  INTEGER :: timer_restart_indices_setup

  ! Data assimilation
  INTEGER :: timer_datass, timer_lhn

  ! The purpose of these "extra" timers is to have otherwise unused timers available for
  ! special-purpose measurements. Please do not remove them and do not use them permanently.
  INTEGER :: timer_extra1,  timer_extra2,  timer_extra3,  timer_extra4,  timer_extra5,  &
       &     timer_extra6,  timer_extra7,  timer_extra8,  timer_extra9,  timer_extra10, &
       &     timer_extra11, timer_extra12, timer_extra13, timer_extra14, timer_extra15, &
       &     timer_extra16, timer_extra17, timer_extra18, timer_extra19, timer_extra20, &
       &     timer_extra21, timer_extra22, timer_extra23, timer_extra24, timer_extra25, &
       &     timer_extra26, timer_extra27, timer_extra28, timer_extra29, timer_extra30, &
       &     timer_extra31, timer_extra32, timer_extra33, timer_extra34, timer_extra35, &
       &     timer_extra36, timer_extra37, timer_extra38, timer_extra39, timer_extra40

  INTEGER :: timer_ls_forcing

  ! ART timers around the ART interfaces
  INTEGER :: timer_art, timer_art_emissInt, timer_art_reacInt,                                  &
             timer_art_cldInt, timer_art_diagInt, timer_art_sedInt, timer_art_toolInt,          &
             timer_art_tracInt, timer_art_turbdiffInt, timer_art_washoutInt, timer_art_initInt, &
             timer_art_radInt, timer_art_photo, timer_art_losschem

CONTAINS

  SUBROUTINE print_timer
#ifdef __SCT__

    USE mo_util_sysinfo,    ONLY: util_user_name, util_os_system, util_node_name
    USE mo_util_vcs,        ONLY: util_repository_url, util_branch_name, util_revision_key
    USE mtime,              ONLY: timedelta, newTimedelta, deallocateTimedelta, &
         &                        OPERATOR(-), timedeltaToString, max_timedelta_str_len
    USE mo_time_config,     ONLY: time_config 
    USE mo_parallel_config, ONLY: get_nproma
    USE mo_run_config,      ONLY: nlev
    USE mo_grid_config,     ONLY: nroot, start_lev
    
    INTEGER :: istat
    
    CHARACTER(len=256) :: repository  = ''
    CHARACTER(len=256) :: branch      = ''
    CHARACTER(len=256) :: revision    = ''
    CHARACTER(len=256) :: executable  = ''
    CHARACTER(len=256) :: user_name   = ''
    CHARACTER(len=256) :: os_name     = ''
    CHARACTER(len=256) :: host_name   = ''
    CHARACTER(len=256) :: expname     = ''
    CHARACTER(len=256) :: jobid       = ''
    CHARACTER(len=256) :: jobname     = ''
    CHARACTER(len=256) :: submit_date = ''        
    CHARACTER(len=256) :: tmp_string  = ''
    
    INTEGER :: nlen, nlena, nlenb, nlenc, nlend

    CHARACTER(len=max_timedelta_str_len) :: tdstring
    TYPE(timedelta), POINTER :: length_of_run

    CHARACTER(len=6) :: gridstring
    
    nlen = 256
    call util_repository_url(repository, nlen)
    nlen = 256
    call util_branch_name(branch, nlen)
    nlen = 256
    call util_revision_key(revision, nlen)
    
    tmp_string = ''
    CALL util_os_system (tmp_string, nlena)
    os_name = tmp_string(1:nlena)
    
    tmp_string = ''
    CALL util_user_name (tmp_string, nlenb)
    user_name = tmp_string(1:nlenb)
    
    tmp_string = ''
    CALL util_node_name (tmp_string, nlenc)
    host_name = tmp_string(1:nlenc)

    CALL get_command_argument(0, executable, nlend)
    
    CALL get_environment_variable('EXPNAME', expname, status=istat) 
    CALL get_environment_variable('SCT_JOB_ID', jobid, status=istat)
    CALL get_environment_variable('SCT_JOB_NAME', jobname, status=istat)
    CALL get_environment_variable('SCT_SUBMIT_DATE', submit_date, status=istat)         
    
    ! sct end date             missing, to be done in sct

    ! model simulation time    tc_stopdate-tc_startdate
    length_of_run => newTimedelta("PT0S")
    length_of_run = time_config%tc_stopdate-time_config%tc_startdate
    CALL timedeltaToString(length_of_run, tdstring)

    write(gridstring,'(a,i2.2,a,i2.2)') 'R', nroot, 'B', start_lev
    
    CALL sct_add_report_attribute('model',                 'icon')
    CALL sct_add_report_attribute('executable',            executable)
    CALL sct_add_report_attribute('repository',            repository)
    CALL sct_add_report_attribute('branch',                branch)
    CALL sct_add_report_attribute('revision',              revision)
    CALL sct_add_report_attribute('user name',             user_name)
    CALL sct_add_report_attribute('operating system name', os_name)
    CALL sct_add_report_attribute('experiment name',       expname)
    CALL sct_add_report_attribute('job id',                jobid)
    CALL sct_add_report_attribute('job name',              jobname)
    CALL sct_add_report_attribute('submit date',           submit_date)    
    CALL sct_add_report_attribute('run length',            tdstring)
    CALL sct_add_report_attribute('vertical levels',       nlev)
    CALL sct_add_report_attribute('horizontal grid',       gridstring)
    CALL sct_add_report_attribute('nproma',                get_nproma())    
#if (defined(__GFORTRAN__) || defined(_CRAYFTN) || defined(__PGIF90__))
    CALL sct_add_report_attribute('compiler version',      compiler_version())
    CALL sct_add_report_attribute('compiler options',      compiler_options())
#else
    CALL sct_add_report_attribute('compiler version',      'unknown')
    CALL sct_add_report_attribute('compiler options',      'unknown')    
#endif
    
    CALL sct_report()

    CALL deallocateTimedelta(length_of_run)
#else
    CALL timer_report()
#endif
  END SUBROUTINE print_timer
    
  SUBROUTINE init_timer

#ifdef __SCT__
    CALL sct_init(timer_max=512)
#endif

    ! major timers
    timer_total        = new_timer("total")

    IF (.NOT. ltimer)  return

    timer_barrier  = new_timer("mpi_barrier")
    timer_exch_data = new_timer("exch_data")
    timer_exch_data_rv = new_timer("exch_data_rv")
    timer_exch_data_async = new_timer("exch_data_async")
    timer_exch_data_wait = new_timer("exch_data.wait")
    timer_global_sum = new_timer("global_sum")
    timer_omp_global_sum = new_timer("omp_global_sum")
    timer_ordglb_sum = new_timer("ordglb_sum")
!     timer_omp_ordglb_sum = new_timer("omp_ordglb_sum")
    timer_icon_comm_sync         = new_timer("icon_comm_sync")
    timer_icon_comm_fillsend     = new_timer("comm_fillsend")
    timer_icon_comm_fillandsend  = new_timer("comm_fillandsend")
    timer_icon_comm_fillrecv     = new_timer("comm_fillrecv")
    timer_icon_comm_ircv         = new_timer("comm_ircv")
    timer_icon_comm_isend        = new_timer("comm_isend")
    timer_icon_comm_send         = new_timer("comm_send")
    timer_icon_comm_wait         = new_timer("comm_wait")
    timer_icon_comm_barrier_2    = new_timer("comm_barrier_2")
    timer_gmres_p_sum            = new_timer("gmres_p_sum")

    timer_write_output  = new_timer("wrt_output")

    timer_integrate_nh      = new_timer  ("integrate_nh")
    timer_solve_nh          = new_timer  ("nh_solve")
    timer_solve_nh_veltend  = new_timer  ("nh_solve.veltend")
    timer_solve_nh_cellcomp = new_timer  ("nh_solve.cellcomp")
    timer_solve_nh_edgecomp = new_timer  ("nh_solve.edgecomp")
    timer_solve_nh_vnupd    = new_timer  ("nh_solve.vnupd")
    timer_solve_nh_vimpl    = new_timer  ("nh_solve.vimpl")
    timer_solve_nh_exch     = new_timer  ("nh_solve.exch")

    timer_step_2tl_si = new_timer("2tl_si_solve")
    timer_step_RK     = new_timer("RK_solve")
    timer_nh_hdiffusion= new_timer("nh_hdiff")

    timer_physics   = new_timer("physics")
    timer_transport = new_timer("transport")
    timer_back_traj = new_timer("back_traj")
    timer_dyn_theta = new_timer("dyn_theta")
    timer_dyn_temp  = new_timer("dyn_temp")

    timer_held_suarez_intr = new_timer("held_suarez_intr")

    ! dynamics timers
    timer_RK_tend = new_timer("RK_tend")
    timer_RK_update = new_timer("RK_update")
    timer_si_correction = new_timer("si_correction")

    timer_intrp_diagn = new_timer   ("intrp_diagn")
    timer_prep_tracer = new_timer   ("prep_tracer")
    timer_prep_tracer_RK = new_timer("prep_tracer_RK")
    timer_hdiff_expl = new_timer    ("hdiff_expl")
    timer_prep_tracer_leapfrog = new_timer("prep_trc_leapfrog")
    timer_div       = new_timer("div")
    timer_grad      = new_timer("grad")
    timer_corio     = new_timer("corio")
    timer_intp      = new_timer("intp")

    ! atmosphere-ocean coupling
    timer_coupling        = new_timer("coupling")
    timer_coupling_init   = new_timer("coupling_init")
    timer_coupling_1stget = new_timer("coupling_1stget")
    timer_coupling_get    = new_timer("coupling_get")
    timer_coupling_put    = new_timer("coupling_put")

    IF (iforcing == iecham) THEN
       !
       ! iconam - echam coupling
       timer_iconam_echam= new_timer("iconam_echam")
       timer_dyn2phy     = new_timer("dyn2phy")
       timer_d2p_prep    = new_timer("d2p_prep")
       timer_d2p_sync    = new_timer("d2p_sync")
       timer_d2p_couple  = new_timer("d2p_couple")
       timer_echam_bcs   = new_timer("echam_bcs")
       timer_echam_phy   = new_timer("echam_phy")
       timer_phy2dyn     = new_timer("phy2dyn")
       timer_p2d_prep    = new_timer("p2d_prep")
       timer_p2d_sync    = new_timer("p2d_sync")
       timer_p2d_couple  = new_timer("p2d_couple")
       !
       ! echam physics
       timer_cov    = new_timer("interface_echam_cov")
       timer_rad    = new_timer("interface_echam_rad")
       timer_rht    = new_timer("interface_echam_rht")
       timer_vdf    = new_timer("interface_echam_vdf")
       timer_vdf_dn = new_timer("vdiff_down")
       timer_vdf_sf = new_timer("update_surface")
       timer_vdf_up = new_timer("vdiff_up")
       timer_gwd    = new_timer("interface_echam_gwd")
       timer_sso    = new_timer("interface_echam_sso")
       timer_cnv    = new_timer("interface_echam_cnv")
       timer_cld    = new_timer("interface_echam_cld")
       timer_car    = new_timer("interface_echam_car")
       timer_mox    = new_timer("interface_echam_mox")
       timer_wmo    = new_timer("interface_echam_wmo")
       !
    END IF
    !
    ! radiation
    timer_rrtm_prep = new_timer("rrtm_prep")
    timer_rrtm_post = new_timer("rrtm_post")
    timer_lrtm      = new_timer("lrtm")
    timer_srtm      = new_timer("srtm")
#ifdef PSRAD_TIMING
    timer_gas_optics_lw = new_timer("gas_optics_lw")
    timer_gas_optics_sw = new_timer("gas_optics_sw")
    timer_rrtm_coeffs = new_timer("rrtm_coeffs")
    timer_cloud_optics = new_timer("cloud_optics")
    timer_sample_cloud_lw = new_timer("sample_cloud")
    timer_sample_cloud_sw = new_timer("sample_cloud")
    timer_psrad_scaling = new_timer("psrad_scaling")
    timer_psrad_aerosol = new_timer("psrad_aerosol")
#endif

    ! nwp physics timers
    timer_omp_radiation = new_timer("omp_radiation")
    timer_nwp_radiation = new_timer("nwp_radiation")
    timer_radiaton_recv = new_timer("radiaton_recv")
    timer_radiaton_comp  = new_timer("radiaton_comp")
    timer_radiaton_send  = new_timer("radiaton_send")
    timer_preradiaton = new_timer("preradiaton")
    timer_phys_acc = new_timer("phys_acc_sync")
    timer_phys_exner = new_timer("phys_exner")
    timer_phys_acc_1 = new_timer("phys_acc_1")
    timer_phys_acc_2 = new_timer("phys_acc_2")
    timer_phys_dpsdt = new_timer("phys_dpsdt")
    timer_phys_sync_tracers = new_timer("phys_sync_tracer")
    timer_phys_sync_tempv    = new_timer("phys_sync_tempv")
    timer_phys_acc_par  = new_timer("phys_acc_par")
    timer_phys_sync_ddt_u  = new_timer("phys_sync_ddt_u")
    timer_phys_sync_vn  = new_timer("phys_sync_vn")
    timer_prep_echam_phy = new_timer("prep_echam_phy")
    timer_prep_phy = new_timer("prep_phy")

    timer_update_prog_phy = new_timer("update_prog_phy")
    timer_nh_diagnostics = new_timer("nh_diagnostics")
    timer_diagnose_pres_temp = new_timer("diagnose_pres_temp")
    timer_satad_v_3D = new_timer("satad")
    timer_phys_u_v = new_timer("phys_u_v")
    timer_nwp_turbulence = new_timer("nwp_turbulence")
    timer_nwp_surface = new_timer("nwp_surface")
    timer_nwp_microphysics = new_timer("nwp_microphysics")
    timer_phys_sync_patch = new_timer("phys_sync_patch")
    timer_fast_phys = new_timer("rediag_prog_vars")
    timer_nwp_convection = new_timer("nwp_convection")
    timer_pre_radiation_nwp = new_timer("pre_radiation_nwp")
    IF (iforcing/=iecham) timer_sso = new_timer("sso")
    timer_cover_koe = new_timer("cloud_cover")
    timer_radiation = new_timer("radiation")
    timer_radheat   = new_timer("radheat")
    timer_synsat    = new_timer("synsat")

    timer_model_init    = new_timer("model_init")
    timer_domain_decomp = new_timer("compute_domain_decomp")
    timer_compute_coeffs = new_timer("compute_intp_coeffs")
    timer_ext_data      = new_timer("init_ext_data")
    timer_init_icon     = new_timer("init_icon") 
    timer_init_latbc    = new_timer("init_latbc") 
    timer_read_restart  = new_timer("read_restart_files")
    timer_solve_ab      = new_timer("solve_ab")
    timer_upd_phys      = new_timer("upd_phys_param")
    timer_upd_flx       = new_timer("upd_flx")
    timer_calc_moc      = new_timer("calc_moc")
    timer_ab_expl       = new_timer("ab_expl")
    timer_ab_rhs4sfc    = new_timer("ab_rhs4sfc")
    timer_tracer_ab     = new_timer("tracer_ab")
    timer_adv_horz      = new_timer("adv_horiz")
    timer_dif_horz      = new_timer("dif_horiz")
    timer_hflx_lim      = new_timer("hflx_lim")
    timer_adv_vert      = new_timer("adv_vert")
    timer_dif_vert      = new_timer("dif_vert")
    timer_ppm_slim      = new_timer("ppm_slim")
    timer_adpo_vert     = new_timer("adpo_vert")
    timer_vert_veloc    = new_timer("vert_veloc")
    timer_normal_veloc  = new_timer("normal_veloc")
    timer_dbg_prnt      = new_timer("dbg_prnt")

    timer_cube_root = new_timer("cube_root")
    timer_lonlat_setup = new_timer("lonlat_setup")

    timer_con_l_theta2t = new_timer("con_l_theta2t")
    timer_con_l_t2theta = new_timer("con_l_t2theta")
    timer_con_theta2t = new_timer("con_theta2t")
    timer_con_t2theta = new_timer("con_t2theta")

    ! timers for boundary interpolation, feedback & nudging
    timer_nesting    = new_timer("nesting")
    timer_nudging    = new_timer("nesting.nudging")
    timer_bdy_interp = new_timer("nesting.bdy_interp")
    timer_feedback   = new_timer("nesting.feedback")

    !ocean timers
    timer_gmres     = new_timer("gmres")
    timer_lhs       = new_timer("lhs")
    timer_lhs_sp    = new_timer("lhs_sp")
    timer_scalar_prod_veloc =new_timer("veloc_prod")
    
    ! Timer IDs for sea ice
    timer_ice_fast      = new_timer("ice_fast")
    timer_ice_slow      = new_timer("ice_slow")
    timer_ice_slow2     = new_timer("ice_slow2")
    timer_ice_momentum  = new_timer("ice_momentum")
    timer_ice_advection = new_timer("ice_advection")
    timer_ice_interp    = new_timer("ice_interp")
  
    ! Timer IDs for HAMOCC
    timer_bgc_up_bgc  = new_timer("hamocc_update_bgc")
    timer_bgc_swr     = new_timer("hamocc_swr")
    timer_bgc_wea     = new_timer("hamocc_weathering")
    timer_bgc_depo    = new_timer("hamocc_dust_depo")
    timer_bgc_chemcon = new_timer("hamocc_chemcon")
    timer_bgc_ocprod  = new_timer("hamocc_ocprod")
    timer_bgc_sett    = new_timer("hamocc_settling") 
    timer_bgc_cya     = new_timer("hamocc_cyanos")
    timer_bgc_gx      = new_timer("hamocc_gas_ex") 
    timer_bgc_calc    = new_timer("hamocc_calc_dissol")
    timer_bgc_powach  = new_timer("hamocc_powach") 
    timer_bgc_up_ic   = new_timer("hamocc_update_icon")
    timer_bgc_tend    = new_timer("hamocc_up_tendencies") 
    timer_bgc_ini     = new_timer("hamocc_ini") 
    timer_bgc_inv     = new_timer("hamocc_inventories") 
    timer_bgc_tot     = new_timer("hamocc_total") 

    ! timers for restart writing/loading
    timer_load_restart = new_timer("load_restart")
    timer_load_restart_io = new_timer("load_restart_io")
    timer_load_restart_comm_setup = new_timer("load_restart_comm_setup")
    timer_load_restart_communication = new_timer("load_restart_communication")
    timer_load_restart_get_var_id = new_timer("load_restart_get_var_id")
    timer_write_restart = new_timer("write_restart")
    timer_write_restart_io = new_timer("write_restart_io")
    timer_write_restart_communication = new_timer("write_restart_communication")
    timer_write_restart_wait = new_timer("write_restart_wait")
    timer_write_restart_setup = new_timer("write_restart_setup")
    timer_restart_indices_setup = new_timer("write_restart_indices")
    timer_restart_collector_setup = new_timer("write_restart_collectors")


    ! Timer for data assimilation
    timer_datass  = new_timer("datass")
    timer_lhn     = new_timer("lhn")

  ! extra timers for on-demand (non-permanent) timings
    timer_extra1  = new_timer("extra1")
    timer_extra2  = new_timer("extra2")
    timer_extra3  = new_timer("extra3")
    timer_extra4  = new_timer("extra4")
    timer_extra5  = new_timer("extra5")
    timer_extra6  = new_timer("extra6")
    timer_extra7  = new_timer("extra7")
    timer_extra8  = new_timer("extra8")
    timer_extra9  = new_timer("extra9")
    timer_extra10 = new_timer("extra10")
    timer_extra11 = new_timer("extra11")
    timer_extra12 = new_timer("extra12")
    timer_extra13 = new_timer("extra13")
    timer_extra14 = new_timer("extra14")
    timer_extra15 = new_timer("extra15")
    timer_extra16 = new_timer("extra16")
    timer_extra17 = new_timer("extra17")
    timer_extra18 = new_timer("extra18")
    timer_extra19 = new_timer("extra19")
    timer_extra20 = new_timer("extra20")
    timer_extra21 = new_timer("extra21")
    timer_extra22 = new_timer("extra22")
    timer_extra23 = new_timer("extra23")
    timer_extra24 = new_timer("extra24")
    timer_extra25 = new_timer("extra25")
    timer_extra26 = new_timer("extra26")
    timer_extra27 = new_timer("extra27")
    timer_extra28 = new_timer("extra28")
    timer_extra29 = new_timer("extra29")
    timer_extra30 = new_timer("extra30")
    timer_extra31 = new_timer("extra31")
    timer_extra32 = new_timer("extra32")
    timer_extra33 = new_timer("extra33")
    timer_extra34 = new_timer("extra34")
    timer_extra35 = new_timer("extra35")
    timer_extra36 = new_timer("extra36")
    timer_extra37 = new_timer("extra37")
    timer_extra38 = new_timer("extra38")
    timer_extra39 = new_timer("extra39")
    timer_extra40 = new_timer("extra40")

    timer_ls_forcing = new_timer("ls_forcing")
   
    ! ART timers around the ART interfaces
    timer_art = new_timer("ART")
    timer_art_emissInt = new_timer("art_emissInt") 
    timer_art_reacInt = new_timer("art_reacInt") 
    timer_art_photo = new_timer("art_photo") 
    timer_art_losschem = new_timer("art_losschem")
    timer_art_cldInt = new_timer("art_cldInt") 
    timer_art_diagInt = new_timer("art_diagInt") 
    timer_art_initInt = new_timer("art_initInt")
    timer_art_radInt = new_timer("art_radInt")
    timer_art_sedInt = new_timer("art_sedInt") 
    timer_art_toolInt = new_timer("art_toolInt")
    timer_art_tracInt = new_timer("art_tracInt") 
    timer_art_turbdiffInt = new_timer("art_turbdiffInt") 
    timer_art_washoutInt = new_timer("art_washoutInt")
  END SUBROUTINE init_timer


  !> Low-level timing routine: start timing.
  !
  !  @note Currently implemented for multi-threaded runs only!
  SUBROUTINE tic(time_s)
!$  USE OMP_LIB

    REAL, INTENT(OUT) :: time_s
    ! local variables:
    LOGICAL :: lopenmp

    lopenmp = .FALSE.
!$  lopenmp = .TRUE.

    IF (.NOT. lopenmp) THEN
      ! do nothing
      time_s = 0.
    ELSE
!$    time_s = REAL(omp_get_wtime())
    END IF
  END SUBROUTINE tic

  !> Low-level timing routine: stop timing, return elapsed time in
  !> seconds.
  !
  !  @note Currently implemented for multi-threaded runs only!
  FUNCTION toc(time_s)
!$  USE OMP_LIB

    REAL :: toc
    REAL, INTENT(IN) :: time_s
    ! local variables:
    LOGICAL :: lopenmp

    lopenmp = .FALSE.
!$  lopenmp = .TRUE.

    IF (.NOT. lopenmp) THEN
      ! do nothing
      toc = 0.
    ELSE
!$    toc = REAL(omp_get_wtime()) - time_s
    END IF
  END FUNCTION toc

END MODULE mo_timer
