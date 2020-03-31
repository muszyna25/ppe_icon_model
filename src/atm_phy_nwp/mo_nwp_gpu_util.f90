MODULE mo_nwp_gpu_util

  USE mo_ext_data_types,          ONLY: t_external_data
  USE mo_nwp_phy_types,           ONLY: t_nwp_phy_diag
  USE mo_model_domain,            ONLY: t_patch
  USE mo_dynamics_config,         ONLY: nnow, nnew, nnow_rcf, nnew_rcf
  USE mo_turbdiff_config,         ONLY: turbdiff_config
  USE mo_intp_data_strc,          ONLY: t_int_state
  USE mo_nh_prepadv_types,        ONLY: t_prepare_adv
  USE mo_grf_intp_data_strc,      ONLY: t_gridref_state, t_gridref_single_state
  USE mo_nwp_parameters,          ONLY: t_phy_params
  USE mo_nonhydrostatic_config,   ONLY: kstart_moist, kstart_tracer
  USE mo_grid_config,             ONLY: n_dom
  USE mo_nwp_phy_state,           ONLY: phy_params

#ifdef _OPENACC
  USE mo_var_list_gpu,            ONLY: gpu_h2d_var_list, gpu_d2h_var_list
#endif

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: gpu_d2h_nh_nwp, gpu_h2d_nh_nwp, devcpy_nwp, hostcpy_nwp

  CONTAINS

  SUBROUTINE gpu_d2h_nh_nwp(pt_patch, prm_diag, ext_data, p_int, prep_adv, p_grf, p_grf_single)

    TYPE(t_patch), TARGET, INTENT(in) :: pt_patch
    TYPE(t_nwp_phy_diag), INTENT(inout) :: prm_diag
    TYPE(t_external_data), OPTIONAL, INTENT(inout):: ext_data
    TYPE(t_int_state), OPTIONAL, INTENT(inout) :: p_int
    TYPE(t_prepare_adv), OPTIONAL, INTENT(inout) :: prep_adv
    TYPE(t_gridref_state), OPTIONAL, INTENT(inout) :: p_grf
    TYPE(t_gridref_single_state), OPTIONAL, INTENT(inout) :: p_grf_single

    INTEGER :: jg

    !$ACC UPDATE HOST(prm_diag%qrs_flux)

    !$ACC UPDATE HOST(ext_data%atm%gp_count_t, ext_data%atm%lp_count_t, ext_data%atm%list_seaice%ncount, &
    !$ACC             ext_data%atm%list_seaice%idx, ext_data%atm%list_lake%ncount, ext_data%atm%list_lake%idx, &
    !$ACC             ext_data%atm%list_seawtr%ncount, ext_data%atm%list_seawtr%idx, ext_data%atm%emis_rad, &
    !$ACC             ext_data%atm%z0_lcc, ext_data%atm%z0_lcc_min, ext_data%atm%plcovmax_lcc, &
    !$ACC             ext_data%atm%laimax_lcc, ext_data%atm%rootdmax_lcc, ext_data%atm%stomresmin_lcc, &
    !$ACC             ext_data%atm%snowalb_lcc, ext_data%atm%snowtile_lcc ) &
    !$acc        IF(PRESENT(ext_data))

    !$ACC UPDATE HOST(p_int%lsq_high, p_int%lsq_lin, &
    !$ACC             p_int%c_bln_avg, p_int%c_lin_e, p_int%cells_aw_verts, &
    !$ACC             p_int%e_bln_c_s, p_int%e_flx_avg, p_int%geofac_div, &
    !$ACC             p_int%geofac_grdiv, p_int%geofac_grg, p_int%geofac_n2s, &
    !$ACC             p_int%geofac_rot, p_int%lsq_high%lsq_blk_c, &
    !$ACC             p_int%lsq_high%lsq_dim_stencil, p_int%lsq_high%lsq_idx_c, &
    !$ACC             p_int%lsq_high%lsq_moments, p_int%lsq_high%lsq_moments_hat, &
    !$ACC             p_int%lsq_high%lsq_pseudoinv, p_int%lsq_high%lsq_qtmat_c, &
    !$ACC             p_int%lsq_high%lsq_rmat_utri_c, p_int%lsq_high%lsq_weights_c, &
    !$ACC             p_int%lsq_lin%lsq_blk_c, &
    !$ACC             p_int%lsq_lin%lsq_dim_stencil, p_int%lsq_lin%lsq_idx_c, &
    !$ACC             p_int%lsq_lin%lsq_moments, p_int%lsq_lin%lsq_moments_hat, &
    !$ACC             p_int%lsq_lin%lsq_pseudoinv, p_int%lsq_lin%lsq_qtmat_c, &
    !$ACC             p_int%lsq_lin%lsq_rmat_utri_c, p_int%lsq_lin%lsq_weights_c, &
    !$ACC             p_int%nudgecoeff_e, p_int%pos_on_tplane_e, &
    !$ACC             p_int%rbf_c2grad_blk, p_int%rbf_c2grad_idx, p_int%rbf_c2grad_coeff, &
    !$ACC             p_int%rbf_vec_blk_c, p_int%rbf_vec_idx_c, p_int%rbf_vec_coeff_c, &
    !$ACC             p_int%rbf_vec_blk_e, p_int%rbf_vec_idx_e, p_int%rbf_vec_coeff_e, &
    !$ACC             p_int%rbf_vec_blk_v, p_int%rbf_vec_idx_v, p_int%rbf_vec_coeff_v, &
    !$ACC             p_int%verts_aw_cells) &
    !$ACC        IF(PRESENT(p_int))        

    !$ACC UPDATE HOST(prep_adv%vn_traj,prep_adv%mass_flx_me,prep_adv%mass_flx_ic,prep_adv%topflx_tra) &
    !$ACC        IF(PRESENT(prep_adv))

    !$ACC UPDATE HOST(p_grf%fbk_wgt_aw, p_grf%fbk_wgt_bln, p_grf%fbk_wgt_e, p_grf%fbk_dom_area, &
    !$ACC             p_grf%mask_ovlp_c, p_grf%mask_ovlp_ch, p_grf%mask_ovlp_e, p_grf%mask_ovlp_v, &
    !$ACC             p_grf%idxlist_bdyintp_src_c, p_grf%idxlist_bdyintp_src_e, p_grf%blklist_bdyintp_src_c, &
    !$ACC             p_grf%blklist_bdyintp_src_e,p_grf%p_dom)  &
    !$ACC        IF(PRESENT(p_grf))    

    !$ACC UPDATE HOST(p_grf_single%grf_dist_pc2cc, p_grf_single%grf_dist_pe2ce, p_grf_single%idxlist_bdyintp_c, &
    !$ACC             p_grf_single%idxlist_bdyintp_e, p_grf_single%idxlist_ubcintp_c, p_grf_single%idxlist_ubcintp_e, p_grf_single%blklist_bdyintp_c, &
    !$ACC             p_grf_single%blklist_bdyintp_e, p_grf_single%blklist_ubcintp_c, p_grf_single%blklist_ubcintp_e, p_grf_single%idxlist_rbfintp_v, &
    !$ACC             p_grf_single%blklist_rbfintp_v, p_grf_single%edge_vert_idx, p_grf_single%coeff_bdyintp_c, p_grf_single%coeff_ubcintp_c, &
    !$ACC             p_grf_single%dist_pc2cc_bdy, p_grf_single%dist_pc2cc_ubc, p_grf_single%prim_norm, p_grf_single%coeff_bdyintp_e12, &
    !$ACC             p_grf_single%coeff_bdyintp_e34, p_grf_single%dist_pe2ce, p_grf_single%coeff_ubcintp_e12, p_grf_single%coeff_ubcintp_e34, &
    !$ACC             p_grf_single%coeff_rbf_v ) &
    !$ACC        IF(PRESENT(p_grf_single))        

#ifdef _OPENACC
    jg = pt_patch%id

    ! Update NWP fields
    CALL gpu_d2h_var_list('prm_diag_of_domain_', domain=jg)
    CALL gpu_d2h_var_list('prm_tend_of_domain_', domain=jg)
    CALL gpu_d2h_var_list('lnd_prog_of_domain_', domain=jg, substr='_and_timelev_', timelev=nnow(jg))
    CALL gpu_d2h_var_list('lnd_prog_of_domain_', domain=jg, substr='_and_timelev_', timelev=nnew(jg))
    CALL gpu_d2h_var_list('lnd_diag_of_domain_', domain=jg)
    CALL gpu_d2h_var_list('wtr_prog_of_domain_', domain=jg, substr='_and_timelev_', timelev=nnow(jg))
    CALL gpu_d2h_var_list('wtr_prog_of_domain_', domain=jg, substr='_and_timelev_', timelev=nnew(jg))
    CALL gpu_d2h_var_list('ext_data_atm_D', domain=jg)

    ! Update dynamics fields
    CALL gpu_d2h_var_list('nh_state_metrics_of_domain_', domain=jg)
    CALL gpu_d2h_var_list('nh_state_diag_of_domain_', domain=jg)
    CALL gpu_d2h_var_list('nh_state_prog_of_domain_', domain=jg, substr='_and_timelev_', timelev=nnew(jg)) !p_prog
    CALL gpu_d2h_var_list('nh_state_prog_of_domain_', domain=jg, substr='_and_timelev_', timelev=nnow_rcf(jg)) !p_prog_now_rcf
    CALL gpu_d2h_var_list('nh_state_prog_of_domain_', domain=jg, substr='_and_timelev_', timelev=nnew_rcf(jg)) !p_prog_rcf
#endif

  END SUBROUTINE gpu_d2h_nh_nwp

  !-------------------------------------------------------------------------
  !-------------------------------------------------------------------------

  SUBROUTINE gpu_h2d_nh_nwp(pt_patch, prm_diag, ext_data, p_int, prep_adv, p_grf, p_grf_single)

    TYPE(t_patch), TARGET, INTENT(in) :: pt_patch
    TYPE(t_nwp_phy_diag), INTENT(inout) :: prm_diag
    TYPE(t_external_data), OPTIONAL, INTENT(inout):: ext_data
    TYPE(t_int_state), OPTIONAL, INTENT(inout) :: p_int
    TYPE(t_prepare_adv), OPTIONAL, INTENT(inout) :: prep_adv
    TYPE(t_gridref_state), OPTIONAL, INTENT(inout) :: p_grf
    TYPE(t_gridref_single_state), OPTIONAL, INTENT(inout) :: p_grf_single


    INTEGER :: jg

    !$ACC UPDATE DEVICE(prm_diag%qrs_flux)

    !$ACC UPDATE DEVICE(ext_data%atm%gp_count_t, ext_data%atm%lp_count_t, ext_data%atm%list_seaice%ncount, &
    !$ACC               ext_data%atm%list_seaice%idx, ext_data%atm%list_lake%ncount, ext_data%atm%list_lake%idx, &
    !$ACC               ext_data%atm%list_seawtr%ncount, ext_data%atm%list_seawtr%idx, ext_data%atm%emis_rad, &
    !$ACC               ext_data%atm%z0_lcc, ext_data%atm%z0_lcc_min, ext_data%atm%plcovmax_lcc, &
    !$ACC               ext_data%atm%laimax_lcc, ext_data%atm%rootdmax_lcc, ext_data%atm%stomresmin_lcc, &
    !$ACC               ext_data%atm%snowalb_lcc, ext_data%atm%snowtile_lcc ) &
    !$acc        IF(PRESENT(ext_data))

    !$ACC UPDATE DEVICE(p_int%lsq_high, p_int%lsq_lin, &
    !$ACC               p_int%c_bln_avg, p_int%c_lin_e, p_int%cells_aw_verts, &
    !$ACC               p_int%e_bln_c_s, p_int%e_flx_avg, p_int%geofac_div, &
    !$ACC               p_int%geofac_grdiv, p_int%geofac_grg, p_int%geofac_n2s, &
    !$ACC               p_int%geofac_rot, p_int%lsq_high%lsq_blk_c, &
    !$ACC               p_int%lsq_high%lsq_dim_stencil, p_int%lsq_high%lsq_idx_c, &
    !$ACC               p_int%lsq_high%lsq_moments, p_int%lsq_high%lsq_moments_hat, &
    !$ACC               p_int%lsq_high%lsq_pseudoinv, p_int%lsq_high%lsq_qtmat_c, &
    !$ACC               p_int%lsq_high%lsq_rmat_utri_c, p_int%lsq_high%lsq_weights_c, &
    !$ACC               p_int%lsq_lin%lsq_blk_c, &
    !$ACC               p_int%lsq_lin%lsq_dim_stencil, p_int%lsq_lin%lsq_idx_c, &
    !$ACC               p_int%lsq_lin%lsq_moments, p_int%lsq_lin%lsq_moments_hat, &
    !$ACC               p_int%lsq_lin%lsq_pseudoinv, p_int%lsq_lin%lsq_qtmat_c, &
    !$ACC               p_int%lsq_lin%lsq_rmat_utri_c, p_int%lsq_lin%lsq_weights_c, &
    !$ACC               p_int%nudgecoeff_e, p_int%pos_on_tplane_e, &
    !$ACC               p_int%rbf_c2grad_blk, p_int%rbf_c2grad_idx, p_int%rbf_c2grad_coeff, &
    !$ACC               p_int%rbf_vec_blk_c, p_int%rbf_vec_idx_c, p_int%rbf_vec_coeff_c, &
    !$ACC               p_int%rbf_vec_blk_e, p_int%rbf_vec_idx_e, p_int%rbf_vec_coeff_e, &
    !$ACC               p_int%rbf_vec_blk_v, p_int%rbf_vec_idx_v, p_int%rbf_vec_coeff_v, &
    !$ACC               p_int%verts_aw_cells) &
    !$ACC        IF(PRESENT(p_int))        

    !$ACC UPDATE DEVICE(prep_adv%vn_traj,prep_adv%mass_flx_me,prep_adv%mass_flx_ic,prep_adv%topflx_tra) &
    !$ACC        IF(PRESENT(prep_adv))

    !$ACC UPDATE DEVICE(p_grf%fbk_wgt_aw, p_grf%fbk_wgt_bln, p_grf%fbk_wgt_e, p_grf%fbk_dom_area, &
    !$ACC               p_grf%mask_ovlp_c, p_grf%mask_ovlp_ch, p_grf%mask_ovlp_e, p_grf%mask_ovlp_v, &
    !$ACC               p_grf%idxlist_bdyintp_src_c, p_grf%idxlist_bdyintp_src_e, p_grf%blklist_bdyintp_src_c, &
    !$ACC               p_grf%blklist_bdyintp_src_e,p_grf%p_dom)  &
    !$ACC        IF(PRESENT(p_grf))    

    !$ACC UPDATE DEVICE(p_grf_single%grf_dist_pc2cc, p_grf_single%grf_dist_pe2ce, p_grf_single%idxlist_bdyintp_c, &
    !$ACC               p_grf_single%idxlist_bdyintp_e, p_grf_single%idxlist_ubcintp_c, p_grf_single%idxlist_ubcintp_e, p_grf_single%blklist_bdyintp_c, &
    !$ACC               p_grf_single%blklist_bdyintp_e, p_grf_single%blklist_ubcintp_c, p_grf_single%blklist_ubcintp_e, p_grf_single%idxlist_rbfintp_v, &
    !$ACC               p_grf_single%blklist_rbfintp_v, p_grf_single%edge_vert_idx, p_grf_single%coeff_bdyintp_c, p_grf_single%coeff_ubcintp_c, &
    !$ACC               p_grf_single%dist_pc2cc_bdy, p_grf_single%dist_pc2cc_ubc, p_grf_single%prim_norm, p_grf_single%coeff_bdyintp_e12, &
    !$ACC               p_grf_single%coeff_bdyintp_e34, p_grf_single%dist_pe2ce, p_grf_single%coeff_ubcintp_e12, p_grf_single%coeff_ubcintp_e34, &
    !$ACC               p_grf_single%coeff_rbf_v ) &
    !$ACC        IF(PRESENT(p_grf_single))        
    
#ifdef _OPENACC
    jg = pt_patch%id

    ! Update NWP fields
    CALL gpu_h2d_var_list('prm_diag_of_domain_', domain=jg)
    CALL gpu_h2d_var_list('prm_tend_of_domain_', domain=jg)
    CALL gpu_h2d_var_list('lnd_prog_of_domain_', domain=jg, substr='_and_timelev_', timelev=nnow(jg))
    CALL gpu_h2d_var_list('lnd_prog_of_domain_', domain=jg, substr='_and_timelev_', timelev=nnew(jg))
    CALL gpu_h2d_var_list('lnd_diag_of_domain_', domain=jg)
    CALL gpu_h2d_var_list('wtr_prog_of_domain_', domain=jg, substr='_and_timelev_', timelev=nnow(jg))
    CALL gpu_h2d_var_list('wtr_prog_of_domain_', domain=jg, substr='_and_timelev_', timelev=nnew(jg))
    CALL gpu_h2d_var_list('lnd_diag_of_domain_', domain=jg)
    CALL gpu_h2d_var_list('ext_data_atm_D', domain=jg)

    ! Update dynamics fields
    CALL gpu_h2d_var_list('nh_state_metrics_of_domain_', domain=jg)
    CALL gpu_h2d_var_list('nh_state_diag_of_domain_', domain=jg)
    CALL gpu_h2d_var_list('nh_state_prog_of_domain_', domain=jg, substr='_and_timelev_', timelev=nnew(jg)) !p_prog
    CALL gpu_h2d_var_list('nh_state_prog_of_domain_', domain=jg, substr='_and_timelev_', timelev=nnow_rcf(jg)) !p_prog_now_rcf
    CALL gpu_h2d_var_list('nh_state_prog_of_domain_', domain=jg, substr='_and_timelev_', timelev=nnew_rcf(jg)) !p_prog_new_rcf
#endif

  END SUBROUTINE gpu_h2d_nh_nwp

  SUBROUTINE devcpy_nwp()

    !$ACC UPDATE DEVICE(phy_params)

    !$ACC ENTER DATA COPYIN(kstart_moist, kstart_tracer)

  END SUBROUTINE devcpy_nwp

  SUBROUTINE hostcpy_nwp()

    !$ACC EXIT DATA DELETE(kstart_moist, kstart_tracer)

  END SUBROUTINE hostcpy_nwp

END MODULE mo_nwp_gpu_util
