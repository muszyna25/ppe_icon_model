!>
!! @author <name, affiliation>
!!
!! @par Revision History
!! <Description of activity> by <name, affiliation> (<YYYY-MM-DD>)
!!
!! @par Copyright
!! 2002-2010 by DWD and MPI-M
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
MODULE mo_timer

  USE mo_real_timer, ONLY: new_timer,                        &
       &                   timer_start,                      &
       &                   timer_stop,                       &
       &                   print_timer   => timer_report,    &
       &                   cleanup_timer => timer_reset_all, &
       &                   delete_timer => del_timer

   USE mo_run_config, ONLY: ltimer, timers_level,  activate_sync_timers

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ltimer, timers_level, activate_sync_timers
  PUBLIC :: new_timer, timer_start, timer_stop  !< procedures imported from mo_real_timer
  PUBLIC :: print_timer, cleanup_timer, delete_timer          !< procedures imported and renamed
  PUBLIC :: init_timer                          !< procedure of this module

  PUBLIC :: timer_total                         !< IDs of timers
  PUBLIC :: timer_exch_data, timer_exch_data_rv, timer_exch_data_async
  PUBLIC :: timer_global_sum, timer_omp_global_sum, timer_ordglb_sum, timer_omp_ordglb_sum
  PUBLIC :: timer_icon_comm_sync  
  PUBLIC :: timer_barrier  

  PUBLIC :: timer_integrate_nh
  PUBLIC :: timer_solve_nh
  PUBLIC :: timer_physics
                        !< IDs of timers
  PUBLIC :: timer_radiation

  PUBLIC :: timer_lrtm_1, timer_lrtm_2

  PUBLIC :: timer_div, timer_grad, timer_gmres  !< IDs of timers
  PUBLIC :: timer_corio, timer_intp             !< IDS of timers
  PUBLIC :: timer_transport
  PUBLIC :: timer_cover, timer_cloud
  PUBLIC :: timer_cucall
  PUBLIC :: timer_vdiff
  PUBLIC :: timer_gw_hines
  PUBLIC :: timer_echam_phy
  PUBLIC :: timer_dyn2phy, timer_phy2dyn
  PUBLIC :: timer_echam_sync_temp,timer_echam_sync_tracers
  PUBLIC :: timer_nh_hdiffusion
  
  PUBLIC :: timer_update_prog_phy
  PUBLIC :: timer_diagnose_pres_temp
  PUBLIC :: timer_nh_diagnostics
  
  PUBLIC :: timer_satad_v_3D
  PUBLIC :: timer_phys_exner
  PUBLIC :: timer_phys_u_v
  PUBLIC :: timer_nwp_turbulence
  PUBLIC :: timer_nwp_microphysics
  PUBLIC :: timer_phys_sync_patch
  PUBLIC :: timer_fast_phys
  PUBLIC :: timer_nwp_convection
  PUBLIC :: timer_nwp_radiation
  PUBLIC :: timer_pre_radiation_nwp
  PUBLIC :: timer_radheat
  PUBLIC :: timer_phys_acc, timer_phys_acc_1,timer_phys_acc_2
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
  PUBLIC :: timer_omp_radiation
  PUBLIC :: timer_lonlat_setup
  PUBLIC :: timer_write_restart_file
  PUBLIC :: timer_write_output
  PUBLIC :: timer_model_init
  PUBLIC :: timer_solve_ab, timer_tracer_ab, timer_vert_veloc, timer_normal_veloc, timer_oce_init
  PUBLIC :: timer_upd_phys, timer_upd_flx
  PUBLIC :: timer_ab_expl, timer_ab_rhs4sfc
  PUBLIC :: timer_adv_horz, timer_dif_horz, timer_hflx_lim
  PUBLIC :: timer_adv_vert, timer_dif_vert, timer_ppm_slim
  PUBLIC :: timer_print_mxmn
  PUBLIC :: timer_si_correction
  PUBLIC :: timer_cube_root
  PUBLIC :: timer_coupling
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

  !-------------------
  ! Module variables
  !-------------------

  ! ID of timer for total model integration time
  INTEGER :: timer_total
  INTEGER :: timer_exch_data, timer_exch_data_rv, timer_exch_data_async
  INTEGER :: timer_global_sum, timer_omp_global_sum, timer_ordglb_sum, timer_omp_ordglb_sum
  INTEGER :: timer_icon_comm_sync
  INTEGER :: timer_barrier
  INTEGER :: timer_nh_hdiffusion

  INTEGER :: timer_integrate_nh
  INTEGER :: timer_solve_nh
  INTEGER :: timer_physics
  INTEGER :: timer_update_prog_phy

  INTEGER :: timer_nh_diagnostics
  INTEGER :: timer_diagnose_pres_temp  
  INTEGER :: timer_satad_v_3D
  INTEGER :: timer_phys_exner
  INTEGER :: timer_phys_u_v
  INTEGER :: timer_nwp_turbulence
  INTEGER :: timer_nwp_microphysics
  INTEGER :: timer_phys_sync_patch
  INTEGER :: timer_fast_phys
  INTEGER :: timer_nwp_convection
  INTEGER :: timer_nwp_radiation
  INTEGER :: timer_pre_radiation_nwp
  INTEGER :: timer_phys_acc, timer_phys_acc_1,timer_phys_acc_2
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
  
  ! Timer ID's for horizontal operators
  INTEGER :: timer_div
  INTEGER :: timer_grad
  INTEGER :: timer_gmres
  INTEGER :: timer_corio
  INTEGER :: timer_intp

  INTEGER :: timer_transport    ! tracer transport

  ! Timer ID's for ECHAM6 physics
  INTEGER :: timer_cover
  INTEGER :: timer_cloud
  INTEGER :: timer_radiation
  INTEGER :: timer_lrtm_1, timer_lrtm_2
  INTEGER :: timer_radheat
  INTEGER :: timer_cucall
  INTEGER :: timer_vdiff
  INTEGER :: timer_gw_hines
  INTEGER :: timer_echam_phy

  ! Timer ID's for forcings and testcases
  INTEGER :: timer_held_suarez_intr
  
  ! Timer ID's for physics-dynamics coupling
  
  INTEGER :: timer_dyn2phy
  INTEGER :: timer_phy2dyn
  INTEGER :: timer_echam_sync_temp, timer_echam_sync_tracers

  INTEGER :: timer_omp_radiation
  INTEGER :: timer_write_restart_file
  INTEGER :: timer_write_output
  INTEGER :: timer_model_init
  INTEGER :: timer_solve_ab, timer_tracer_ab, timer_vert_veloc, timer_normal_veloc, timer_oce_init
  INTEGER :: timer_upd_phys, timer_upd_flx
  INTEGER :: timer_ab_expl, timer_ab_rhs4sfc
  INTEGER :: timer_adv_horz, timer_dif_horz, timer_hflx_lim
  INTEGER :: timer_adv_vert, timer_dif_vert, timer_ppm_slim
  INTEGER :: timer_print_mxmn
  INTEGER :: timer_si_correction
  INTEGER :: timer_cube_root
  INTEGER :: timer_coupling
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

CONTAINS

  SUBROUTINE init_timer

    ! major timers
    timer_total        = new_timer("total")

    timer_exch_data = new_timer("exch_data")
    timer_exch_data_rv = new_timer("exch_data_rv")
    timer_exch_data_async = new_timer("exch_data_async")
    timer_global_sum = new_timer("global_sum")
    timer_omp_global_sum = new_timer("omp_global_sum")
    timer_ordglb_sum = new_timer("ordglb_sum")
    timer_omp_ordglb_sum = new_timer("omp_ordglb_sum")
    timer_icon_comm_sync = new_timer("icon_comm_sync")
    timer_barrier  = new_timer("barrier")
      
    timer_coupling      = new_timer("coupling")
    timer_write_output  = new_timer("wrt_output")
    timer_write_restart_file = new_timer("wrt_restart")
 
    timer_integrate_nh= new_timer  ("ntegrate_nh")
    timer_solve_nh    = new_timer  ("nh_solve")
    timer_step_2tl_si = new_timer("2tl_si_solve")
    timer_step_RK     = new_timer("RK_solve")
    timer_nh_hdiffusion= new_timer("nh_hdiff")
   
    timer_physics   = new_timer("physics")
    timer_echam_phy = new_timer("echam_phy")

    timer_transport = new_timer("transport")
    timer_dyn_theta = new_timer("dyn_theta")
    timer_dyn_temp  = new_timer("dyn_temp")
    
    timer_held_suarez_intr = new_timer("held_suarez_intr")
    
    timer_gw_hines  = new_timer("gw_hines")

    ! dynamics timers
    timer_gmres     = new_timer("gmres")
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

    ! physics timers
    timer_radiation = new_timer("radiation")    
    timer_lrtm_1    = new_timer("rad_lrtm_1")
    timer_lrtm_2    = new_timer("rad_lrtm_2")
    timer_omp_radiation = new_timer("omp_radiation")
    timer_nwp_radiation = new_timer("nwp_radiation")
    timer_radheat = new_timer("radheat")
    timer_cover     = new_timer("cover")
    timer_cloud     = new_timer("cloud")
    timer_cucall    = new_timer("cucall")
    timer_vdiff     = new_timer("vdiff")
    timer_dyn2phy   = new_timer("dyn2phy")
    timer_phy2dyn   = new_timer("phy2dyn")
    timer_echam_sync_temp= new_timer("echam_sync_temp")
    timer_echam_sync_tracers= new_timer("echam_sync_tracers")
    timer_phys_acc = new_timer("phys_acc")
    timer_phys_exner = new_timer("phys_exner")
    timer_phys_acc_1 = new_timer("phys_acc_1")
    timer_phys_acc_2 = new_timer("phys_acc_2")
    timer_phys_sync_tracers = new_timer("phys_sync_tracer")
    timer_phys_sync_tempv    = new_timer("phys_sync_tempv")
    timer_phys_acc_par  = new_timer("phys_acc_par")
    timer_phys_sync_ddt_u  = new_timer("phys_sync_ddt_u")
    timer_phys_sync_vn  = new_timer("phys_sync_vn")
    timer_prep_echam_phy = new_timer("prep_echam_phy")
    timer_prep_phy = new_timer("prep_phy")
    
    timer_update_prog_phy = new_timer("update_prog_phy")
 
    

!     timer_sync_delay = new_timer("sync_delay")
!     timer_sync_outbuffer = new_timer("sync_outbuffer")
!     timer_sync_psend_1 = new_timer("sync_psend_1")
!     timer_sync_isend_2 = new_timer("sync_isend_2")
!     timer_sync_recv_2 = new_timer("sync_recv_2")
!     timer_sync_isend_3 = new_timer("sync_isend_3")
!     timer_sync_wait = new_timer("sync_wait")

    timer_nh_diagnostics = new_timer("nh_diagnostics")
    
    timer_diagnose_pres_temp = new_timer("diagnose_pres_temp")
    timer_satad_v_3D = new_timer("satad_v_3D")
    timer_phys_u_v = new_timer("phys_u_v")
    timer_nwp_turbulence = new_timer("nwp_turbulence")
    timer_nwp_microphysics = new_timer("nwp_microphysics")
    timer_phys_sync_patch = new_timer("phys_sync_patch")
    timer_fast_phys = new_timer("fast_phys")
    timer_nwp_convection = new_timer("nwp_convection")
    timer_pre_radiation_nwp = new_timer("pre_radiation_nwp")
    timer_sso = new_timer("sso")
    timer_cover_koe = new_timer("cover_koe")
        

    timer_model_init    = new_timer("model_init")
    timer_oce_init      = new_timer("oce_init")
    timer_solve_ab      = new_timer("solve_ab")
    timer_upd_phys      = new_timer("upd_phys")
    timer_upd_flx       = new_timer("upd_flx")
    timer_ab_expl       = new_timer("ab_expl")
    timer_ab_rhs4sfc    = new_timer("ab_rhs4sfc")
    timer_tracer_ab     = new_timer("tracer_ab")
    timer_adv_horz      = new_timer("adv_horz")
    timer_dif_horz      = new_timer("dif_horz")
    timer_hflx_lim      = new_timer("hflx_lim")
    timer_adv_vert      = new_timer("adv_vert")
    timer_dif_vert      = new_timer("dif_vert")
    timer_ppm_slim      = new_timer("ppm_slim")
    timer_vert_veloc    = new_timer("vert_veloc")
    timer_normal_veloc  = new_timer("normal_veloc")
    timer_print_mxmn    = new_timer("print_mxmn")
  
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

  END SUBROUTINE init_timer

END MODULE mo_timer








