!>

!! Type definition for the GPU implementation of the dynamical core of ICONAM.
!!
!! @author William Sawyer (CSCS)
!!
!! @par Revision History
!! Initial release by William Sawyer (2015)
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
MODULE mo_nonhydro_gpu_types

#if defined( _OPENACC )

  USE mo_kind,                 ONLY: wp, vp
  USE mo_mpi,                  ONLY: i_am_accel_node
  USE mo_fortran_tools,        ONLY: t_ptr_2d3d
  USE mo_math_types,           ONLY: t_geographical_coordinates
  USE mo_model_domain,         ONLY: t_patch, t_tangent_vectors
  USE mo_nonhydro_types,       ONLY: t_nh_state, t_nh_diag, t_nh_prog
  USE mo_nh_prepadv_types,     ONLY: t_prepare_adv
  USE mo_intp_data_strc,       ONLY: t_int_state

  IMPLICIT NONE

  PRIVATE 

  PUBLIC :: h2d_icon, d2h_icon

CONTAINS

  SUBROUTINE h2d_icon( p_int_states, p_patches, p_nh_states, prep_advs )

    TYPE ( t_int_state ),  INTENT(INOUT) :: p_int_states(:)
    TYPE ( t_patch ),      INTENT(INOUT) :: p_patches(:)
    TYPE ( t_nh_state ),   INTENT(INOUT) :: p_nh_states(:)
    TYPE ( t_prepare_adv), INTENT(INOUT) :: prep_advs(:)
    INTEGER :: jg
!
! Copy all data need on GPU from host to device
!

#ifndef _CRAYFTN
!$ACC ENTER DATA COPYIN( p_int_states, p_patches, prep_advs, p_nh_states ), IF ( i_am_accel_node  )
#endif

    CALL h2d_int_state( p_int_states )

    CALL h2d_patch( p_patches )

    CALL h2d_prep_adv( prep_advs )

    CALL h2d_nh_state( p_nh_states )

  CONTAINS

    SUBROUTINE h2d_int_state( p_int )

      TYPE ( t_int_state ), TARGET,  INTENT(INOUT) :: p_int(:)
#ifdef _CRAYFTN
      INTEGER,  POINTER, DIMENSION(:,:)    :: h_lsq_dim_stencil, l_lsq_dim_stencil
      INTEGER,  POINTER, DIMENSION(:,:,:)  :: h_lsq_idx_c, h_lsq_blk_c, l_lsq_idx_c, l_lsq_blk_c,           &
                                              rbf_vec_idx_c, rbf_vec_blk_c, rbf_vec_idx_e, rbf_vec_blk_e,   &
                                              rbf_vec_idx_v, rbf_vec_blk_v, rbf_c2grad_idx, rbf_c2grad_blk
      REAL(wp), POINTER, DIMENSION(:,:)    :: nudgecoeff_e
      REAL(wp), POINTER, DIMENSION(:,:,:)  :: c_bln_avg, c_lin_e, cells_aw_verts, e_bln_c_s, e_flx_avg,     &
                                              geofac_div, geofac_grdiv, geofac_n2s, geofac_qdiv, geofac_rot,&
                                              h_lsq_moments, h_lsq_rmat_utri_c, h_lsq_weights_c,            &
                                              l_lsq_moments, l_lsq_rmat_utri_c, l_lsq_weights_c,            &
                                              rbf_vec_coeff_e
      REAL(wp), POINTER, DIMENSION(:,:,:,:):: geofac_grg, h_lsq_moments_hat, h_lsq_pseudoinv, h_lsq_qtmat_c,&
                                              l_lsq_moments_hat, l_lsq_pseudoinv, l_lsq_qtmat_c,            &
                                              pos_on_tplane_e, rbf_vec_coeff_c, rbf_vec_coeff_v,            &
                                              rbf_c2grad_coeff
#endif
      INTEGER  :: j

      DO j=1, SIZE(p_int)

#ifdef _CRAYFTN

        c_bln_avg                            => p_int(j)%c_bln_avg
        c_lin_e                              => p_int(j)%c_lin_e
        cells_aw_verts                       => p_int(j)%cells_aw_verts
        e_bln_c_s                            => p_int(j)%e_bln_c_s
        e_flx_avg                            => p_int(j)%e_flx_avg
        geofac_div                           => p_int(j)%geofac_div
        geofac_grdiv                         => p_int(j)%geofac_grdiv
        geofac_grg                           => p_int(j)%geofac_grg
        geofac_n2s                           => p_int(j)%geofac_n2s
        geofac_qdiv                          => p_int(j)%geofac_qdiv
        geofac_rot                           => p_int(j)%geofac_rot
        h_lsq_blk_c                          => p_int(j)%lsq_high%lsq_blk_c
        h_lsq_dim_stencil                    => p_int(j)%lsq_high%lsq_dim_stencil
        h_lsq_idx_c                          => p_int(j)%lsq_high%lsq_idx_c
        h_lsq_moments                        => p_int(j)%lsq_high%lsq_moments
        h_lsq_moments_hat                    => p_int(j)%lsq_high%lsq_moments_hat
        h_lsq_pseudoinv                      => p_int(j)%lsq_high%lsq_pseudoinv
        h_lsq_qtmat_c                        => p_int(j)%lsq_high%lsq_qtmat_c
        h_lsq_rmat_utri_c                    => p_int(j)%lsq_high%lsq_rmat_utri_c
        h_lsq_weights_c                      => p_int(j)%lsq_high%lsq_weights_c
        l_lsq_blk_c                          => p_int(j)%lsq_lin%lsq_blk_c
        l_lsq_dim_stencil                    => p_int(j)%lsq_lin%lsq_dim_stencil
        l_lsq_idx_c                          => p_int(j)%lsq_lin%lsq_idx_c
        l_lsq_moments                        => p_int(j)%lsq_lin%lsq_moments
        l_lsq_moments_hat                    => p_int(j)%lsq_lin%lsq_moments_hat
        l_lsq_pseudoinv                      => p_int(j)%lsq_lin%lsq_pseudoinv
        l_lsq_qtmat_c                        => p_int(j)%lsq_lin%lsq_qtmat_c
        l_lsq_rmat_utri_c                    => p_int(j)%lsq_lin%lsq_rmat_utri_c
        l_lsq_weights_c                      => p_int(j)%lsq_lin%lsq_weights_c
        nudgecoeff_e                         => p_int(j)%nudgecoeff_e
        pos_on_tplane_e                      => p_int(j)%pos_on_tplane_e
        rbf_c2grad_blk                       => p_int(j)%rbf_c2grad_blk
        rbf_c2grad_idx                       => p_int(j)%rbf_c2grad_idx
        rbf_c2grad_coeff                     => p_int(j)%rbf_c2grad_coeff
        rbf_vec_blk_c                        => p_int(j)%rbf_vec_blk_c
        rbf_vec_idx_c                        => p_int(j)%rbf_vec_idx_c
        rbf_vec_coeff_c                      => p_int(j)%rbf_vec_coeff_c
        rbf_vec_blk_e                        => p_int(j)%rbf_vec_blk_e
        rbf_vec_idx_e                        => p_int(j)%rbf_vec_idx_e
        rbf_vec_coeff_e                      => p_int(j)%rbf_vec_coeff_e
        rbf_vec_blk_v                        => p_int(j)%rbf_vec_blk_v
        rbf_vec_idx_v                        => p_int(j)%rbf_vec_idx_v
        rbf_vec_coeff_v                      => p_int(j)%rbf_vec_coeff_v

!$ACC ENTER DATA &
!$ACC       COPYIN( c_bln_avg, c_lin_e, cells_aw_verts, e_bln_c_s, e_flx_avg,                  &
!$ACC               geofac_div, geofac_grdiv, geofac_grg, geofac_n2s, geofac_qdiv, geofac_rot, &
!$ACC               h_lsq_blk_c, h_lsq_dim_stencil, h_lsq_idx_c, h_lsq_moments,                &
!$ACC               h_lsq_moments_hat, h_lsq_pseudoinv, h_lsq_qtmat_c, h_lsq_rmat_utri_c,      &
!$ACC               h_lsq_weights_c, l_lsq_blk_c, l_lsq_dim_stencil, l_lsq_idx_c,              &
!$ACC               l_lsq_moments, l_lsq_moments_hat, l_lsq_pseudoinv, l_lsq_qtmat_c,          &
!$ACC               l_lsq_rmat_utri_c, l_lsq_weights_c, nudgecoeff_e, pos_on_tplane_e,         &
!$ACC               rbf_c2grad_blk, rbf_c2grad_idx, rbf_c2grad_coeff,                          &
!$ACC               rbf_vec_blk_c, rbf_vec_idx_c, rbf_vec_coeff_c,                             &
!$ACC               rbf_vec_blk_e, rbf_vec_idx_e, rbf_vec_coeff_e,                             &
!$ACC               rbf_vec_blk_v, rbf_vec_idx_v, rbf_vec_coeff_v ),                           &
!$ACC       IF ( i_am_accel_node )

#else

!$ACC ENTER DATA &
!$ACC       COPYIN( p_int(j)%lsq_high, p_int(j)%lsq_lin,                                         &
!$ACC               p_int(j)%c_bln_avg, p_int(j)%c_lin_e, p_int(j)%cells_aw_verts,               &
!$ACC               p_int(j)%e_bln_c_s, p_int(j)%e_flx_avg, p_int(j)%geofac_div,                 &
!$ACC               p_int(j)%geofac_grdiv, p_int(j)%geofac_grg, p_int(j)%geofac_n2s,             &
!$ACC               p_int(j)%geofac_rot, p_int(j)%lsq_high%lsq_blk_c,                            &
!$ACC               p_int(j)%lsq_high%lsq_dim_stencil, p_int(j)%lsq_high%lsq_idx_c,              &
!$ACC               p_int(j)%lsq_high%lsq_moments, p_int(j)%lsq_high%lsq_moments_hat,            &
!$ACC               p_int(j)%lsq_high%lsq_pseudoinv, p_int(j)%lsq_high%lsq_qtmat_c,              &
!$ACC               p_int(j)%lsq_high%lsq_rmat_utri_c, p_int(j)%lsq_high%lsq_weights_c,          &
!$ACC               p_int(j)%lsq_lin%lsq_blk_c,                                                  &
!$ACC               p_int(j)%lsq_lin%lsq_dim_stencil, p_int(j)%lsq_lin%lsq_idx_c,                &
!$ACC               p_int(j)%lsq_lin%lsq_moments, p_int(j)%lsq_lin%lsq_moments_hat,              &
!$ACC               p_int(j)%lsq_lin%lsq_pseudoinv, p_int(j)%lsq_lin%lsq_qtmat_c,                &
!$ACC               p_int(j)%lsq_lin%lsq_rmat_utri_c, p_int(j)%lsq_lin%lsq_weights_c,            &
!$ACC               p_int(j)%nudgecoeff_e, p_int(j)%pos_on_tplane_e,                             &
!$ACC               p_int(j)%rbf_c2grad_blk, p_int(j)%rbf_c2grad_idx, p_int(j)%rbf_c2grad_coeff, &
!$ACC               p_int(j)%rbf_vec_blk_c, p_int(j)%rbf_vec_idx_c, p_int(j)%rbf_vec_coeff_c,    &
!$ACC               p_int(j)%rbf_vec_blk_e, p_int(j)%rbf_vec_idx_e, p_int(j)%rbf_vec_coeff_e,    &
!$ACC               p_int(j)%rbf_vec_blk_v, p_int(j)%rbf_vec_idx_v, p_int(j)%rbf_vec_coeff_v ),  &
!$ACC       IF ( i_am_accel_node )        

#endif

      ENDDO

    END SUBROUTINE h2d_int_state

    SUBROUTINE h2d_patch( p_patch )

      TYPE ( t_patch ), TARGET, INTENT(INOUT) :: p_patch(:)

      INTEGER :: j
#ifdef _CRAYFTN
      INTEGER, POINTER, DIMENSION(:,:)   :: refin_ctrl_c, refin_ctrl_e, refin_ctrl_v
      INTEGER, POINTER, DIMENSION(:,:,:) :: edge_idx, edge_blk, cell_idx, cell_blk, vertex_idx, vertex_blk, &
                                            quad_idx, quad_blk, neighbor_idx, neighbor_blk,                 &
                                            verts_cell_idx, verts_cell_blk, verts_edge_idx, verts_edge_blk
      TYPE(t_tangent_vectors), POINTER   :: primal_normal_cell(:,:,:), dual_normal_cell(:,:,:)
      TYPE(t_tangent_vectors), POINTER   :: primal_normal_vert(:,:,:), dual_normal_vert(:,:,:)
      REAL(wp), POINTER, DIMENSION(:,:)  :: inv_primal_edge_length, inv_dual_edge_length, &
                                            inv_vert_vert_length, tangent_orientation,    &
                                            area, area_edge, f_e
      REAL(wp), POINTER, DIMENSION(:,:,:):: edge_cell_length
      TYPE(t_geographical_coordinates), POINTER   :: center(:,:)
      LOGICAL, POINTER :: owner_mask(:,:)
#endif
!
! Copy the static data structures in p_patch to the device -- this is a small subset of all the components
! The communication patterns are copied over in mo_communication_orig.
!

      DO j=1,SIZE(p_patch)

#ifdef _CRAYFTN

        owner_mask                     =>   p_patch(j)%cells%decomp_info%owner_mask

        area                           =>   p_patch(j)%cells%area
        edge_idx                       =>   p_patch(j)%cells%edge_idx
        edge_blk                       =>   p_patch(j)%cells%edge_blk
        neighbor_idx                   =>   p_patch(j)%cells%neighbor_idx
        neighbor_blk                   =>   p_patch(j)%cells%neighbor_blk
        center                         =>   p_patch(j)%cells%center

        area_edge                      =>   p_patch(j)%edges%area_edge
        cell_idx                       =>   p_patch(j)%edges%cell_idx
        cell_blk                       =>   p_patch(j)%edges%cell_blk
        edge_cell_length               =>   p_patch(j)%edges%edge_cell_length
        f_e                            =>   p_patch(j)%edges%f_e
        quad_idx                       =>   p_patch(j)%edges%quad_idx
        quad_blk                       =>   p_patch(j)%edges%quad_blk
        vertex_idx                     =>   p_patch(j)%edges%vertex_idx
        vertex_blk                     =>   p_patch(j)%edges%vertex_blk

        primal_normal_cell             =>   p_patch(j)%edges%primal_normal_cell
        dual_normal_cell               =>   p_patch(j)%edges%dual_normal_cell
        primal_normal_vert             =>   p_patch(j)%edges%primal_normal_vert
        dual_normal_vert               =>   p_patch(j)%edges%dual_normal_vert

        inv_primal_edge_length         =>   p_patch(j)%edges%inv_primal_edge_length
        inv_dual_edge_length           =>   p_patch(j)%edges%inv_dual_edge_length
        inv_vert_vert_length           =>   p_patch(j)%edges%inv_vert_vert_length
        tangent_orientation            =>   p_patch(j)%edges%tangent_orientation 

        refin_ctrl_c                   =>   p_patch(j)%cells%refin_ctrl
        refin_ctrl_e                   =>   p_patch(j)%edges%refin_ctrl
        refin_ctrl_v                   =>   p_patch(j)%verts%refin_ctrl

        verts_cell_idx                 =>   p_patch(j)%verts%cell_idx
        verts_cell_blk                 =>   p_patch(j)%verts%cell_blk
        verts_edge_idx                 =>   p_patch(j)%verts%edge_idx
        verts_edge_blk                 =>   p_patch(j)%verts%edge_blk

!$ACC ENTER DATA &
!$ACC      COPYIN( area, edge_idx, edge_blk, neighbor_idx, neighbor_blk, center,   &
!$ACC              area_edge, cell_idx, cell_blk, f_e, quad_idx, quad_blk,         &  
!$ACC              vertex_idx, vertex_blk, primal_normal_cell, dual_normal_cell,   &
!$ACC              primal_normal_vert, dual_normal_vert, inv_vert_vert_length,     &
!$ACC              inv_primal_edge_length, inv_dual_edge_length, edge_cell_length, &
!$ACC              owner_mask, refin_ctrl_c, refin_ctrl_e, refin_ctrl_v,           &
!$ACC              tangent_orientation, verts_cell_idx, verts_cell_blk,            &
!$ACC              verts_edge_idx, verts_edge_blk ), &
!$ACC      IF ( i_am_accel_node  )

#else

!$ACC ENTER DATA &
!$ACC      COPYIN( p_patch(j)%cells, p_patch(j)%cells%decomp_info, p_patch(j)%cells%decomp_info%owner_mask, &
!$ACC              p_patch(j)%cells%area, p_patch(j)%cells%edge_idx, p_patch(j)%cells%edge_blk,             &
!$ACC              p_patch(j)%cells%neighbor_idx, p_patch(j)%cells%neighbor_blk,                            &
!$ACC              p_patch(j)%cells%center, p_patch(j)%cells%refin_ctrl,                                    &
!$ACC              p_patch(j)%edges, p_patch(j)%edges%area_edge, p_patch(j)%edges%cell_idx,                 &
!$ACC              p_patch(j)%edges%cell_blk, p_patch(j)%edges%edge_cell_length, p_patch(j)%edges%f_e,      &
!$ACC              p_patch(j)%edges%quad_idx, p_patch(j)%edges%quad_blk, p_patch(j)%edges%vertex_idx,       &
!$ACC              p_patch(j)%edges%vertex_blk, p_patch(j)%edges%primal_normal_cell,                        &
!$ACC              p_patch(j)%edges%dual_normal_cell, p_patch(j)%edges%primal_normal_vert,                  &
!$ACC              p_patch(j)%edges%dual_normal_vert, p_patch(j)%edges%inv_vert_vert_length,                &
!$ACC              p_patch(j)%edges%inv_dual_edge_length, p_patch(j)%edges%inv_primal_edge_length,          &
!$ACC              p_patch(j)%edges%tangent_orientation, p_patch(j)%edges%refin_ctrl,                       &
!$ACC              p_patch(j)%verts, p_patch(j)%verts%cell_idx, p_patch(j)%verts%cell_blk,                  &
!$ACC              p_patch(j)%verts%edge_idx, p_patch(j)%verts%edge_blk, p_patch(j)%verts%refin_ctrl ),     &
!$ACC      IF ( i_am_accel_node  )

#endif

      ENDDO

    END SUBROUTINE h2d_patch


    SUBROUTINE h2d_prep_adv( prep_adv )

      TYPE ( t_prepare_adv ), TARGET, INTENT(INOUT) :: prep_adv(:)

      INTEGER :: j
#ifdef _CRAYFTN
      REAL(wp), POINTER, DIMENSION(:,:,:)  :: vn_traj, mass_flx_me, mass_flx_ic, topflx_tra
#endif

      DO j=1, SIZE(prep_adv)

#ifdef _CRAYFTN

        vn_traj     => prep_adv(j)%vn_traj
        mass_flx_me => prep_adv(j)%mass_flx_me
        mass_flx_ic => prep_adv(j)%mass_flx_ic
        topflx_tra  => prep_adv(j)%topflx_tra
!$ACC ENTER DATA COPYIN( vn_traj, mass_flx_me, mass_flx_ic, topflx_tra ), IF ( i_am_accel_node  )

#else

!$ACC ENTER DATA COPYIN(prep_adv(j)%vn_traj,prep_adv(j)%mass_flx_me,prep_adv(j)%mass_flx_ic,prep_adv(j)%topflx_tra ), &
!$ACC       IF ( i_am_accel_node  )

#endif

      ENDDO    

    END SUBROUTINE h2d_prep_adv


    SUBROUTINE h2d_nh_state( p_nh )

      TYPE ( t_nh_state ), TARGET, INTENT(INOUT) :: p_nh(:)
      INTEGER :: istep, j

#ifdef _CRAYFTN
      INTEGER,  POINTER, DIMENSION(:)      :: bdy_halo_c_blk, bdy_halo_c_idx, bdy_mflx_e_blk, bdy_mflx_e_idx,     &
                                              nudge_e_blk, nudge_e_idx, pg_vertidx, pg_edgeidx, pg_edgeblk
      INTEGER,  POINTER, DIMENSION(:,:)    :: zd_blklist, zd_edgeblk, zd_edgeidx, zd_indlist, zd_vertidx
      INTEGER,  POINTER, DIMENSION(:,:,:,:):: vertidx_gradp
      REAL(wp), POINTER, DIMENSION(:)      :: enhfac_diffu, rayleigh_vn, rayleigh_w, scalfac_dd3d
      REAL(wp), POINTER, DIMENSION(:,:)    :: pres_sfc, vwind_expl_wgt, vwind_impl_wgt,                  &
                                              hmask_dd3d, dvn_ie_ubc, dw_ubc, dtheta_v_ic_ubc,           &
                                              zd_intcoef, zd_geofac, zd_e2cell
      REAL(wp), POINTER, DIMENSION(:,:,:)  :: ddqz_z_full, inv_ddqz_z_full, vn_ref, w_ref, rho_ic,       &
                                              theta_v_ic, dwdx, dwdy, mass_fl_e, div_ic, hdef_ic,        &
                                              exner_pr, mflx_ic_ubc, theta_v_ic_ubc, grf_bdy_mflx,       &
                                              dpres_mc, pres_ifc, pres, temp, tempv, temp_ifc
      REAL(wp), POINTER, DIMENSION(:,:,:,:):: tracer
      REAL(vp), POINTER, DIMENSION(:)      :: pg_exdist
      REAL(vp), POINTER, DIMENSION(:,:,:)  :: d_exner_dz_ref_ic, d2dexdz2_fac1_mc, d2dexdz2_fac2_mc,              &
                                              ddqz_z_half, ddxn_z_full, ddxt_z_full, ddqz_z_full_e,               &
                                              exner_exfac, exner_ref_mc,                                          &
                                              rho_ref_mc, rho_ref_me, theta_ref_ic, theta_ref_mc, theta_ref_me,   &
                                              vt, v, u, vor, vn_ie, w_concorr_c,                                  &
                                              wgtfac_c, wgtfac_e, wgtfacq_c, wgtfacq_e, wgtfacq1_c, wgtfacq1_e,   &
                                              exner_dyn_incr, ddt_exner_phy, ddt_vn_phy, rho_incr, exner_incr,    &
                                              coeff1_dwdz, coeff2_dwdz, coeff_gradekin
                             
      REAL(vp), POINTER, DIMENSION(:,:,:,:):: coeff_gradp, zdiff_gradp, ddt_vn_adv, ddt_w_adv, ddt_tracer_adv

      REAL(wp), POINTER, DIMENSION(:,:,:)  :: exner, rho, theta_v, vn, w

      LOGICAL,  POINTER                    :: mask_prog_halo_c(:,:)
#endif

      DO j = 1, SIZE(p_nh)

#ifdef _CRAYFTN

        bdy_halo_c_blk              => p_nh(j)%metrics%bdy_halo_c_blk
        bdy_halo_c_idx              => p_nh(j)%metrics%bdy_halo_c_idx
        bdy_mflx_e_blk              => p_nh(j)%metrics%bdy_mflx_e_blk
        bdy_mflx_e_idx              => p_nh(j)%metrics%bdy_mflx_e_idx
        coeff1_dwdz                 => p_nh(j)%metrics%coeff1_dwdz
        coeff2_dwdz                 => p_nh(j)%metrics%coeff2_dwdz
        coeff_gradekin              => p_nh(j)%metrics%coeff_gradekin
        coeff_gradp                 => p_nh(j)%metrics%coeff_gradp
        d_exner_dz_ref_ic           => p_nh(j)%metrics%d_exner_dz_ref_ic
        d2dexdz2_fac1_mc            => p_nh(j)%metrics%d2dexdz2_fac1_mc
        d2dexdz2_fac2_mc            => p_nh(j)%metrics%d2dexdz2_fac2_mc
        ddqz_z_full                 => p_nh(j)%metrics%ddqz_z_full
        ddqz_z_full_e               => p_nh(j)%metrics%ddqz_z_full_e
        ddqz_z_half                 => p_nh(j)%metrics%ddqz_z_half
        ddxn_z_full                 => p_nh(j)%metrics%ddxn_z_full
        ddxt_z_full                 => p_nh(j)%metrics%ddxt_z_full
        enhfac_diffu                => p_nh(j)%metrics%enhfac_diffu
        exner_exfac                 => p_nh(j)%metrics%exner_exfac
        exner_exfac                 => p_nh(j)%metrics%exner_exfac
        exner_ref_mc                => p_nh(j)%metrics%exner_ref_mc
        hmask_dd3d                  => p_nh(j)%metrics%hmask_dd3d
        inv_ddqz_z_full             => p_nh(j)%metrics%inv_ddqz_z_full
        mask_prog_halo_c            => p_nh(j)%metrics%mask_prog_halo_c
        nudge_e_blk                 => p_nh(j)%metrics%nudge_e_blk
        nudge_e_idx                 => p_nh(j)%metrics%nudge_e_idx
        pg_exdist                   => p_nh(j)%metrics%pg_exdist
        pg_vertidx                  => p_nh(j)%metrics%pg_vertidx
        pg_edgeidx                  => p_nh(j)%metrics%pg_edgeidx
        pg_edgeblk                  => p_nh(j)%metrics%pg_edgeblk
        rayleigh_vn                 => p_nh(j)%metrics%rayleigh_vn
        rayleigh_w                  => p_nh(j)%metrics%rayleigh_w
        rho_ref_mc                  => p_nh(j)%metrics%rho_ref_mc
        rho_ref_me                  => p_nh(j)%metrics%rho_ref_me
        scalfac_dd3d                => p_nh(j)%metrics%scalfac_dd3d
        theta_ref_ic                => p_nh(j)%metrics%theta_ref_ic
        theta_ref_mc                => p_nh(j)%metrics%theta_ref_mc
        theta_ref_me                => p_nh(j)%metrics%theta_ref_me
        vertidx_gradp               => p_nh(j)%metrics%vertidx_gradp
        vwind_expl_wgt              => p_nh(j)%metrics%vwind_expl_wgt
        vwind_impl_wgt              => p_nh(j)%metrics%vwind_impl_wgt
        wgtfac_c                    => p_nh(j)%metrics%wgtfac_c
        wgtfac_e                    => p_nh(j)%metrics%wgtfac_e
        wgtfacq_c                   => p_nh(j)%metrics%wgtfacq_c
        wgtfacq_e                   => p_nh(j)%metrics%wgtfacq_e
        wgtfacq1_c                  => p_nh(j)%metrics%wgtfacq1_c
        wgtfacq1_e                  => p_nh(j)%metrics%wgtfacq1_e
        zdiff_gradp                 => p_nh(j)%metrics%zdiff_gradp
        zd_indlist                  => p_nh(j)%metrics%zd_indlist
        zd_blklist                  => p_nh(j)%metrics%zd_blklist
        zd_vertidx                  => p_nh(j)%metrics%zd_vertidx
        zd_edgeidx                  => p_nh(j)%metrics%zd_edgeidx
        zd_edgeblk                  => p_nh(j)%metrics%zd_edgeblk
        zd_intcoef                  => p_nh(j)%metrics%zd_intcoef
        zd_e2cell                   => p_nh(j)%metrics%zd_e2cell
        zd_geofac                   => p_nh(j)%metrics%zd_geofac

        vn_ref                      => p_nh(j)%ref%vn_ref
        w_ref                       => p_nh(j)%ref%w_ref

        ddt_exner_phy               => p_nh(j)%diag%ddt_exner_phy
        ddt_vn_adv                  => p_nh(j)%diag%ddt_vn_adv
        ddt_vn_phy                  => p_nh(j)%diag%ddt_vn_phy
        ddt_w_adv                   => p_nh(j)%diag%ddt_w_adv
        dpres_mc                    => p_nh(j)%diag%dpres_mc
        div_ic                      => p_nh(j)%diag%div_ic
        dwdx                        => p_nh(j)%diag%dwdx
        dwdy                        => p_nh(j)%diag%dwdy
        dtheta_v_ic_ubc             => p_nh(j)%diag%dtheta_v_ic_ubc
        dvn_ie_ubc                  => p_nh(j)%diag%dvn_ie_ubc
        dw_ubc                      => p_nh(j)%diag%dw_ubc
        exner_dyn_incr              => p_nh(j)%diag%exner_dyn_incr
        exner_incr                  => p_nh(j)%diag%exner_incr
        exner_pr                    => p_nh(j)%diag%exner_pr
        grf_bdy_mflx                => p_nh(j)%diag%grf_bdy_mflx
        hdef_ic                     => p_nh(j)%diag%hdef_ic
        mass_fl_e                   => p_nh(j)%diag%mass_fl_e
        mflx_ic_ubc                 => p_nh(j)%diag%mflx_ic_ubc
        pres_ifc                    => p_nh(j)%diag%pres_ifc
        pres_sfc                    => p_nh(j)%diag%pres_sfc
        pres                        => p_nh(j)%diag%pres
        rho_ic                      => p_nh(j)%diag%rho_ic
        rho_incr                    => p_nh(j)%diag%rho_incr
        temp                        => p_nh(j)%diag%temp
        tempv                       => p_nh(j)%diag%tempv
        temp_ifc                    => p_nh(j)%diag%temp_ifc
        theta_v_ic                  => p_nh(j)%diag%theta_v_ic
        vn_ie                       => p_nh(j)%diag%vn_ie
        vt                          => p_nh(j)%diag%vt
        v                           => p_nh(j)%diag%v
        u                           => p_nh(j)%diag%u
        vor                         => p_nh(j)%diag%vor
        w_concorr_c                 => p_nh(j)%diag%w_concorr_c
        ddt_tracer_adv              => p_nh(j)%diag%ddt_tracer_adv

!$ACC ENTER DATA &
!$ACC      COPYIN( bdy_halo_c_blk, bdy_halo_c_idx, bdy_mflx_e_blk, bdy_mflx_e_idx, &
!$ACC            coeff1_dwdz, coeff2_dwdz, coeff_gradekin, coeff_gradp,            &
!$ACC            d_exner_dz_ref_ic, d2dexdz2_fac1_mc, d2dexdz2_fac2_mc,            &
!$ACC            ddqz_z_full, ddqz_z_full_e, ddqz_z_half, ddxn_z_full, ddxt_z_full,&
!$ACC            enhfac_diffu, exner_exfac, exner_ref_mc, hmask_dd3d,              &
!$ACC            inv_ddqz_z_full, mask_prog_halo_c, nudge_e_blk, nudge_e_idx,      &
!$ACC            pg_exdist, pg_vertidx, pg_edgeidx, pg_edgeblk,                    &
!$ACC            rayleigh_vn, rayleigh_w, rho_ref_mc, rho_ref_me,                  &
!$ACC            scalfac_dd3d, theta_ref_ic, theta_ref_mc, theta_ref_me,           &
!$ACC            vertidx_gradp, vwind_expl_wgt, vwind_impl_wgt,                    &
!$ACC            wgtfac_c, wgtfac_e, wgtfacq_c, wgtfacq_e, wgtfacq1_c, wgtfacq1_e, &
!$ACC            zdiff_gradp, zd_blklist, zd_e2cell, zd_edgeblk, zd_edgeidx,       &
!$ACC            zd_geofac, zd_indlist, zd_intcoef, zd_vertidx        ),           &
!$ACC      COPYIN( vn_ref, w_ref ),                                                &
!$ACC      COPYIN( ddt_exner_phy, ddt_vn_adv, ddt_vn_phy, ddt_w_adv, div_ic,       &
!$ACC              dpres_mc, dtheta_v_ic_ubc, dvn_ie_ubc, dw_ubc, dwdx, dwdy,      &
!$ACC              exner_dyn_incr, exner_incr, exner_pr, grf_bdy_mflx,             &
!$ACC              hdef_ic, mass_fl_e, mflx_ic_ubc, pres_ifc, pres, pres_sfc,      &
!$ACC              rho_ic, rho_incr, temp, tempv, temp_ifc, theta_v_ic, vt, vn_ie, &
!$ACC              v, u, vor, w_concorr_c, ddt_tracer_adv),                        &
!$ACC      IF ( i_am_accel_node )

        DO istep = 1, SIZE( p_nh(j)%prog )

          exner   =>  p_nh(j)%prog(istep)%exner
          rho     =>  p_nh(j)%prog(istep)%rho
          theta_v =>  p_nh(j)%prog(istep)%theta_v
          tracer  =>  p_nh(j)%prog(istep)%tracer
          vn      =>  p_nh(j)%prog(istep)%vn
          w       =>  p_nh(j)%prog(istep)%w
!$ACC ENTER DATA COPYIN( exner, rho, theta_v, tracer, vn, w ), IF ( i_am_accel_node  )

        ENDDO

#else

!$ACC ENTER DATA &
!$ACC      COPYIN( p_nh(j)%metrics, &
!$ACC              p_nh(j)%metrics%bdy_halo_c_blk, p_nh(j)%metrics%bdy_halo_c_idx,                   &
!$ACC              p_nh(j)%metrics%bdy_mflx_e_blk, p_nh(j)%metrics%bdy_mflx_e_idx,                   &
!$ACC              p_nh(j)%metrics%coeff1_dwdz, p_nh(j)%metrics%coeff2_dwdz,                         &
!$ACC              p_nh(j)%metrics%coeff_gradekin, p_nh(j)%metrics%coeff_gradp,                      &
!$ACC              p_nh(j)%metrics%d_exner_dz_ref_ic, p_nh(j)%metrics%d2dexdz2_fac1_mc,              &
!$ACC              p_nh(j)%metrics%d2dexdz2_fac2_mc, p_nh(j)%metrics%ddqz_z_full,                    &
!$ACC              p_nh(j)%metrics%ddqz_z_full_e, p_nh(j)%metrics%ddqz_z_half,                       &
!$ACC              p_nh(j)%metrics%ddxn_z_full, p_nh(j)%metrics%ddxt_z_full,                         &
!$ACC              p_nh(j)%metrics%enhfac_diffu, p_nh(j)%metrics%exner_exfac,                        &
!$ACC              p_nh(j)%metrics%exner_ref_mc,                                                     &
!$ACC              p_nh(j)%metrics%hmask_dd3d, p_nh(j)%metrics%inv_ddqz_z_full,                      &
!$ACC              p_nh(j)%metrics%mask_prog_halo_c, p_nh(j)%metrics%nudge_e_blk,                    &
!$ACC              p_nh(j)%metrics%nudge_e_idx, p_nh(j)%metrics%pg_exdist,                           &
!$ACC              p_nh(j)%metrics%pg_vertidx, p_nh(j)%metrics%pg_edgeidx,                           &
!$ACC              p_nh(j)%metrics%pg_edgeblk, p_nh(j)%metrics%rayleigh_vn,                          &
!$ACC              p_nh(j)%metrics%rayleigh_w, p_nh(j)%metrics%rho_ref_mc,                           &
!$ACC              p_nh(j)%metrics%rho_ref_me, p_nh(j)%metrics%scalfac_dd3d,                         &
!$ACC              p_nh(j)%metrics%theta_ref_ic, p_nh(j)%metrics%theta_ref_mc,                       &
!$ACC              p_nh(j)%metrics%theta_ref_me, p_nh(j)%metrics%vertidx_gradp,                      &
!$ACC              p_nh(j)%metrics%vwind_expl_wgt, p_nh(j)%metrics%vwind_impl_wgt,                   &
!$ACC              p_nh(j)%metrics%wgtfac_c, p_nh(j)%metrics%wgtfac_e,                               &
!$ACC              p_nh(j)%metrics%wgtfacq_c, p_nh(j)%metrics%wgtfacq_e,                             &
!$ACC              p_nh(j)%metrics%wgtfacq1_c, p_nh(j)%metrics%wgtfacq1_e,                           &
!$ACC              p_nh(j)%metrics%zdiff_gradp, p_nh(j)%metrics%zd_blklist,                          &
!$ACC              p_nh(j)%metrics%zd_e2cell, p_nh(j)%metrics%zd_edgeblk,                            &
!$ACC              p_nh(j)%metrics%zd_edgeidx, p_nh(j)%metrics%zd_geofac, p_nh(j)%metrics%zd_indlist,&
!$ACC              p_nh(j)%metrics%zd_intcoef,p_nh(j)%metrics%zd_vertidx  ),                         &
!$ACC      COPYIN( p_nh(j)%ref, p_nh(j)%ref%vn_ref, p_nh(j)%ref%w_ref ),                             &
!$ACC      COPYIN( p_nh(j)%diag, p_nh(j)%diag%ddt_exner_phy, p_nh(j)%diag%ddt_vn_adv,                &
!$ACC              p_nh(j)%diag%ddt_vn_phy, p_nh(j)%diag%ddt_w_adv,                                  &
!$ACC              p_nh(j)%diag%div_ic, p_nh(j)%diag%dpres_mc,                                       &
!$ACC              p_nh(j)%diag%dtheta_v_ic_ubc, p_nh(j)%diag%dvn_ie_ubc,                            &
!$ACC              p_nh(j)%diag%dw_ubc, p_nh(j)%diag%dwdx, p_nh(j)%diag%dwdy,                        &
!$ACC              p_nh(j)%diag%exner_dyn_incr, p_nh(j)%diag%exner_incr, p_nh(j)%diag%exner_pr,      &
!$ACC              p_nh(j)%diag%grf_bdy_mflx, p_nh(j)%diag%hdef_ic,                                  &
!$ACC              p_nh(j)%diag%mass_fl_e, p_nh(j)%diag%mflx_ic_ubc, p_nh(j)%diag%pres_ifc,          &
!$ACC              p_nh(j)%diag%pres_sfc, p_nh(j)%diag%pres, p_nh(j)%diag%rho_ic,                    &
!$ACC              p_nh(j)%diag%rho_incr, p_nh(j)%diag%temp, p_nh(j)%diag%tempv,                     &
!$ACC              p_nh(j)%diag%temp_ifc, p_nh(j)%diag%theta_v_ic, p_nh(j)%diag%vt,                  &
!$ACC              p_nh(j)%diag%v, p_nh(j)%diag%u, p_nh(j)%diag%vor, p_nh(j)%diag%vn_ie,             &
!$ACC              p_nh(j)%diag%w_concorr_c, p_nh(j)%diag%ddt_tracer_adv ),                          &
!$ACC      IF ( i_am_accel_node )

!$ACC ENTER DATA COPYIN( p_nh(j)%prog ), IF ( i_am_accel_node  )

        DO istep = 1, SIZE( p_nh(j)%prog )

!$ACC ENTER DATA COPYIN( p_nh(j)%prog(istep)%exner, p_nh(j)%prog(istep)%rho,                         &
!$ACC                    p_nh(j)%prog(istep)%tracer,  p_nh(j)%prog(istep)%theta_v,                   &
!$ACC                    p_nh(j)%prog(istep)%vn, p_nh(j)%prog(istep)%w ),                            &
!$ACC      IF ( i_am_accel_node  )

        ENDDO

#endif

      ENDDO

    END SUBROUTINE h2d_nh_state

  END SUBROUTINE h2d_icon

  SUBROUTINE d2h_icon( p_int_states, p_patches, p_nh_states, prep_advs )

    TYPE ( t_int_state ),  INTENT(INOUT) :: p_int_states(:)
    TYPE ( t_patch ),      INTENT(INOUT) :: p_patches(:)
    TYPE ( t_nh_state ),   INTENT(INOUT) :: p_nh_states(:)
    TYPE ( t_prepare_adv), INTENT(INOUT) :: prep_advs(:)

    REAL(wp), POINTER, DIMENSION(:,:,:)  :: vn_traj, mass_flx_me, mass_flx_ic

!
! Delete all data on GPU
!

    CALL d2h_prep_adv( prep_advs )
    CALL d2h_nh_state( p_nh_states )
    CALL d2h_patch( p_patches )
    CALL d2h_int_state( p_int_states )

#ifndef _CRAYFTN
!$ACC EXIT DATA DELETE( prep_advs, p_nh_states, p_patches, p_int_states ), IF ( i_am_accel_node  )
#endif

  CONTAINS

    SUBROUTINE d2h_int_state( p_int )

      TYPE ( t_int_state ), TARGET,  INTENT(INOUT) :: p_int(:)
#ifdef _CRAYFTN
      INTEGER,  POINTER, DIMENSION(:,:)    :: h_lsq_dim_stencil, l_lsq_dim_stencil
      INTEGER,  POINTER, DIMENSION(:,:,:)  :: h_lsq_idx_c, h_lsq_blk_c, l_lsq_idx_c, l_lsq_blk_c,           &
                                              rbf_vec_idx_c, rbf_vec_blk_c, rbf_vec_idx_e, rbf_vec_blk_e,   &
                                              rbf_vec_idx_v, rbf_vec_blk_v, rbf_c2grad_idx, rbf_c2grad_blk
      REAL(wp), POINTER, DIMENSION(:,:)    :: nudgecoeff_e
      REAL(wp), POINTER, DIMENSION(:,:,:)  :: c_bln_avg, c_lin_e, cells_aw_verts, e_bln_c_s, e_flx_avg,     &
                                              geofac_div, geofac_grdiv, geofac_n2s, geofac_qdiv, geofac_rot,&
                                              h_lsq_moments, h_lsq_rmat_utri_c, h_lsq_weights_c,            &
                                              l_lsq_moments, l_lsq_rmat_utri_c, l_lsq_weights_c,            &
                                              rbf_vec_coeff_e
      REAL(wp), POINTER, DIMENSION(:,:,:,:):: geofac_grg, h_lsq_moments_hat, h_lsq_pseudoinv, h_lsq_qtmat_c,&
                                              l_lsq_moments_hat, l_lsq_pseudoinv, l_lsq_qtmat_c,            &
                                              pos_on_tplane_e, rbf_vec_coeff_c, rbf_vec_coeff_v,            &
                                              rbf_c2grad_coeff
#endif
      INTEGER  :: j

      DO j=1, SIZE(p_int)

#ifdef _CRAYFTN

        c_bln_avg                            => p_int(j)%c_bln_avg
        c_lin_e                              => p_int(j)%c_lin_e
        cells_aw_verts                       => p_int(j)%cells_aw_verts
        e_bln_c_s                            => p_int(j)%e_bln_c_s
        e_flx_avg                            => p_int(j)%e_flx_avg
        geofac_div                           => p_int(j)%geofac_div
        geofac_grdiv                         => p_int(j)%geofac_grdiv
        geofac_grg                           => p_int(j)%geofac_grg
        geofac_n2s                           => p_int(j)%geofac_n2s
        geofac_qdiv                          => p_int(j)%geofac_qdiv
        geofac_rot                           => p_int(j)%geofac_rot
        h_lsq_blk_c                          => p_int(j)%lsq_high%lsq_blk_c
        h_lsq_dim_stencil                    => p_int(j)%lsq_high%lsq_dim_stencil
        h_lsq_idx_c                          => p_int(j)%lsq_high%lsq_idx_c
        h_lsq_moments                        => p_int(j)%lsq_high%lsq_moments
        h_lsq_moments_hat                    => p_int(j)%lsq_high%lsq_moments_hat
        h_lsq_pseudoinv                      => p_int(j)%lsq_high%lsq_pseudoinv
        h_lsq_qtmat_c                        => p_int(j)%lsq_high%lsq_qtmat_c
        h_lsq_rmat_utri_c                    => p_int(j)%lsq_high%lsq_rmat_utri_c
        h_lsq_weights_c                      => p_int(j)%lsq_high%lsq_weights_c
        l_lsq_blk_c                          => p_int(j)%lsq_lin%lsq_blk_c
        l_lsq_dim_stencil                    => p_int(j)%lsq_lin%lsq_dim_stencil
        l_lsq_idx_c                          => p_int(j)%lsq_lin%lsq_idx_c
        l_lsq_moments                        => p_int(j)%lsq_lin%lsq_moments
        l_lsq_moments_hat                    => p_int(j)%lsq_lin%lsq_moments_hat
        l_lsq_pseudoinv                      => p_int(j)%lsq_lin%lsq_pseudoinv
        l_lsq_qtmat_c                        => p_int(j)%lsq_lin%lsq_qtmat_c
        l_lsq_rmat_utri_c                    => p_int(j)%lsq_lin%lsq_rmat_utri_c
        l_lsq_weights_c                      => p_int(j)%lsq_lin%lsq_weights_c
        nudgecoeff_e                         => p_int(j)%nudgecoeff_e
        pos_on_tplane_e                      => p_int(j)%pos_on_tplane_e
        rbf_c2grad_blk                       => p_int(j)%rbf_c2grad_blk
        rbf_c2grad_idx                       => p_int(j)%rbf_c2grad_idx
        rbf_c2grad_coeff                     => p_int(j)%rbf_c2grad_coeff
        rbf_vec_blk_c                        => p_int(j)%rbf_vec_blk_c
        rbf_vec_idx_c                        => p_int(j)%rbf_vec_idx_c
        rbf_vec_coeff_c                      => p_int(j)%rbf_vec_coeff_c
        rbf_vec_blk_e                        => p_int(j)%rbf_vec_blk_e
        rbf_vec_idx_e                        => p_int(j)%rbf_vec_idx_e
        rbf_vec_coeff_e                      => p_int(j)%rbf_vec_coeff_e
        rbf_vec_blk_v                        => p_int(j)%rbf_vec_blk_v
        rbf_vec_idx_v                        => p_int(j)%rbf_vec_idx_v
        rbf_vec_coeff_v                      => p_int(j)%rbf_vec_coeff_v

!$ACC EXIT DATA &
!$ACC      DELETE(  c_bln_avg, c_lin_e, cells_aw_verts, e_bln_c_s, e_flx_avg,                  &
!$ACC               geofac_div, geofac_grdiv, geofac_grg, geofac_n2s, geofac_qdiv, geofac_rot, &
!$ACC               h_lsq_blk_c, h_lsq_dim_stencil, h_lsq_idx_c, h_lsq_moments,                &
!$ACC               h_lsq_moments_hat, h_lsq_pseudoinv, h_lsq_qtmat_c, h_lsq_rmat_utri_c,      &
!$ACC               h_lsq_weights_c, l_lsq_blk_c, l_lsq_dim_stencil, l_lsq_idx_c,              &
!$ACC               l_lsq_moments, l_lsq_moments_hat, l_lsq_pseudoinv, l_lsq_qtmat_c,          &
!$ACC               l_lsq_rmat_utri_c, l_lsq_weights_c, nudgecoeff_e, pos_on_tplane_e,         &
!$ACC               rbf_c2grad_blk, rbf_c2grad_idx, rbf_c2grad_coeff,                          &
!$ACC               rbf_vec_blk_c, rbf_vec_idx_c, rbf_vec_coeff_c,                             &
!$ACC               rbf_vec_blk_e, rbf_vec_idx_e, rbf_vec_coeff_e,                             &
!$ACC               rbf_vec_blk_v, rbf_vec_idx_v, rbf_vec_coeff_v ),                           &
!$ACC       IF ( i_am_accel_node )

#else

!$ACC EXIT DATA &
!$ACC      DELETE(  p_int(j)%c_bln_avg, p_int(j)%c_lin_e, p_int(j)%cells_aw_verts,               &
!$ACC               p_int(j)%e_bln_c_s, p_int(j)%e_flx_avg, p_int(j)%geofac_div,                 &
!$ACC               p_int(j)%geofac_grdiv, p_int(j)%geofac_grg, p_int(j)%geofac_n2s,             &
!$ACC               p_int(j)%geofac_rot, p_int(j)%lsq_high%lsq_blk_c,                            &
!$ACC               p_int(j)%lsq_high%lsq_dim_stencil, p_int(j)%lsq_high%lsq_idx_c,              &
!$ACC               p_int(j)%lsq_high%lsq_moments, p_int(j)%lsq_high%lsq_moments_hat,            &
!$ACC               p_int(j)%lsq_high%lsq_pseudoinv, p_int(j)%lsq_high%lsq_qtmat_c,              &
!$ACC               p_int(j)%lsq_high%lsq_rmat_utri_c, p_int(j)%lsq_high%lsq_weights_c,          &
!$ACC               p_int(j)%lsq_lin%lsq_blk_c,                                                  &
!$ACC               p_int(j)%lsq_lin%lsq_dim_stencil, p_int(j)%lsq_lin%lsq_idx_c,                &
!$ACC               p_int(j)%lsq_lin%lsq_moments, p_int(j)%lsq_lin%lsq_moments_hat,              &
!$ACC               p_int(j)%lsq_lin%lsq_pseudoinv, p_int(j)%lsq_lin%lsq_qtmat_c,                &
!$ACC               p_int(j)%lsq_lin%lsq_rmat_utri_c, p_int(j)%lsq_lin%lsq_weights_c,            &
!$ACC               p_int(j)%nudgecoeff_e, p_int(j)%pos_on_tplane_e,                             &
!$ACC               p_int(j)%rbf_c2grad_blk, p_int(j)%rbf_c2grad_idx, p_int(j)%rbf_c2grad_coeff, &
!$ACC               p_int(j)%rbf_vec_blk_c, p_int(j)%rbf_vec_idx_c, p_int(j)%rbf_vec_coeff_c,    &
!$ACC               p_int(j)%rbf_vec_blk_e, p_int(j)%rbf_vec_idx_e, p_int(j)%rbf_vec_coeff_e,    &
!$ACC               p_int(j)%rbf_vec_blk_v, p_int(j)%rbf_vec_idx_v, p_int(j)%rbf_vec_coeff_v,    &
!$ACC               p_int(j)%lsq_high, p_int(j)%lsq_lin )                                        &
!$ACC       IF ( i_am_accel_node )        

#endif

      ENDDO

    END SUBROUTINE d2h_int_state

    SUBROUTINE d2h_patch( p_patch )

      TYPE ( t_patch ), TARGET, INTENT(INOUT) :: p_patch(:)
      INTEGER :: j
#ifdef _CRAYFTN
      INTEGER, POINTER, DIMENSION(:,:)   :: refin_ctrl_c, refin_ctrl_e, refin_ctrl_v
      INTEGER, POINTER, DIMENSION(:,:,:) :: edge_idx, edge_blk, cell_idx, cell_blk, vertex_idx, vertex_blk, &
                                            quad_idx, quad_blk, neighbor_idx, neighbor_blk,                 &
                                            verts_cell_idx, verts_cell_blk, verts_edge_idx, verts_edge_blk
      TYPE(t_tangent_vectors), POINTER   :: primal_normal_cell(:,:,:), dual_normal_cell(:,:,:)
      TYPE(t_tangent_vectors), POINTER   :: primal_normal_vert(:,:,:), dual_normal_vert(:,:,:)
      REAL(wp), POINTER, DIMENSION(:,:)  :: inv_primal_edge_length, inv_dual_edge_length, &
                                            inv_vert_vert_length, tangent_orientation,    &
                                            area, area_edge, f_e
      REAL(wp), POINTER, DIMENSION(:,:,:):: edge_cell_length
      TYPE(t_geographical_coordinates), POINTER   :: center(:,:)
      LOGICAL, POINTER :: owner_mask(:,:)
#endif

!
! Copy the static data structures in p_patch to the device -- this is a small subset of all the components
! The communication patterns are copied over in mo_communication_orig.
!

      DO j=1,SIZE(p_patch)

#ifdef _CRAYFTN

        owner_mask                     =>   p_patch(j)%cells%decomp_info%owner_mask

        area                           =>   p_patch(j)%cells%area
        edge_idx                       =>   p_patch(j)%cells%edge_idx
        edge_blk                       =>   p_patch(j)%cells%edge_blk
        neighbor_idx                   =>   p_patch(j)%cells%neighbor_idx
        neighbor_blk                   =>   p_patch(j)%cells%neighbor_blk
        center                         =>   p_patch(j)%cells%center
        refin_ctrl_c                   =>   p_patch(j)%cells%refin_ctrl

        area_edge                      =>   p_patch(j)%edges%area_edge
        cell_idx                       =>   p_patch(j)%edges%cell_idx
        cell_blk                       =>   p_patch(j)%edges%cell_blk
        edge_cell_length               =>   p_patch(j)%edges%edge_cell_length
        f_e                            =>   p_patch(j)%edges%f_e
        quad_idx                       =>   p_patch(j)%edges%quad_idx
        quad_blk                       =>   p_patch(j)%edges%quad_blk
        vertex_idx                     =>   p_patch(j)%edges%vertex_idx
        vertex_blk                     =>   p_patch(j)%edges%vertex_blk

        primal_normal_cell             =>   p_patch(j)%edges%primal_normal_cell
        dual_normal_cell               =>   p_patch(j)%edges%dual_normal_cell
        primal_normal_vert             =>   p_patch(j)%edges%primal_normal_vert
        dual_normal_vert               =>   p_patch(j)%edges%dual_normal_vert

        inv_primal_edge_length         =>   p_patch(j)%edges%inv_primal_edge_length
        inv_dual_edge_length           =>   p_patch(j)%edges%inv_dual_edge_length
        inv_vert_vert_length           =>   p_patch(j)%edges%inv_vert_vert_length
        tangent_orientation            =>   p_patch(j)%edges%tangent_orientation 

        refin_ctrl_e                   =>   p_patch(j)%edges%refin_ctrl

        verts_cell_idx                 =>   p_patch(j)%verts%cell_idx
        verts_cell_blk                 =>   p_patch(j)%verts%cell_blk
        verts_edge_idx                 =>   p_patch(j)%verts%edge_idx
        verts_edge_blk                 =>   p_patch(j)%verts%edge_blk
        refin_ctrl_v                   =>   p_patch(j)%verts%refin_ctrl

!$ACC EXIT DATA &
!$ACC      DELETE( area, edge_idx, edge_blk, neighbor_idx, neighbor_blk, center,   &
!$ACC              area_edge, cell_idx, cell_blk, f_e, quad_idx, quad_blk,         &  
!$ACC              vertex_idx, vertex_blk, primal_normal_cell, dual_normal_cell,   &
!$ACC              primal_normal_vert, dual_normal_vert, inv_vert_vert_length,     &
!$ACC              inv_primal_edge_length, inv_dual_edge_length, edge_cell_length, &
!$ACC              owner_mask, refin_ctrl_c, refin_ctrl_e, refin_ctrl_v,           &
!$ACC              tangent_orientation, verts_cell_idx, verts_cell_blk,            &
!$ACC              verts_edge_idx, verts_edge_blk ),                               &
!$ACC      IF ( i_am_accel_node  )

#else

!$ACC EXIT DATA &
!$ACC      DELETE( p_patch(j)%cells%decomp_info, p_patch(j)%cells%decomp_info%owner_mask,               &
!$ACC              p_patch(j)%cells%area, p_patch(j)%cells%edge_idx, p_patch(j)%cells%edge_blk,         &
!$ACC              p_patch(j)%cells%neighbor_idx, p_patch(j)%cells%neighbor_blk,                        &
!$ACC              p_patch(j)%cells%center, p_patch(j)%cells%refin_ctrl, p_patch(j)%cells,              &
!$ACC              p_patch(j)%edges%area_edge, p_patch(j)%edges%cell_idx,                               &
!$ACC              p_patch(j)%edges%cell_blk, p_patch(j)%edges%edge_cell_length, p_patch(j)%edges%f_e,  &
!$ACC              p_patch(j)%edges%quad_idx, p_patch(j)%edges%quad_blk, p_patch(j)%edges%vertex_idx,   &
!$ACC              p_patch(j)%edges%vertex_blk, p_patch(j)%edges%primal_normal_cell,                    &
!$ACC              p_patch(j)%edges%dual_normal_cell, p_patch(j)%edges%primal_normal_vert,              &
!$ACC              p_patch(j)%edges%dual_normal_vert, p_patch(j)%edges%inv_vert_vert_length,            &
!$ACC              p_patch(j)%edges%inv_dual_edge_length, p_patch(j)%edges%inv_primal_edge_length,      &
!$ACC              p_patch(j)%edges%tangent_orientation, p_patch(j)%edges%refin_ctrl, p_patch(j)%edges, &
!$ACC              p_patch(j)%verts%edge_idx, p_patch(j)%verts%edge_blk, p_patch(j)%verts%refin_ctrl,   &
!$ACC              p_patch(j)%verts  ), &
!$ACC      IF ( i_am_accel_node  )

#endif

      ENDDO

    END SUBROUTINE d2h_patch

    SUBROUTINE d2h_prep_adv( prep_adv )

      TYPE ( t_prepare_adv ), TARGET, INTENT(INOUT) :: prep_adv(:)
      INTEGER :: j
#ifdef _CRAYFTN
      REAL(wp), POINTER, DIMENSION(:,:,:)  :: vn_traj, mass_flx_me, mass_flx_ic, topflx_tra
#endif

      DO j = 1, SIZE( prep_adv )

#ifdef _CRAYFTN

        vn_traj     => prep_adv(j)%vn_traj
        mass_flx_me => prep_adv(j)%mass_flx_me
        mass_flx_ic => prep_adv(j)%mass_flx_ic
        topflx_tra  => prep_adv(j)%topflx_tra
!$ACC EXIT DATA DELETE( vn_traj, mass_flx_me, mass_flx_ic, topflx_tra ), IF ( i_am_accel_node  )

#else
!$ACC EXIT DATA DELETE(prep_adv(j)%vn_traj,prep_adv(j)%mass_flx_me,prep_adv(j)%mass_flx_ic,prep_adv(j)%topflx_tra ), &
!$ACC      IF ( i_am_accel_node  )
#endif

      ENDDO

    END SUBROUTINE d2h_prep_adv

    SUBROUTINE d2h_nh_state( p_nh )

      TYPE ( t_nh_state ), TARGET, INTENT(INOUT) :: p_nh(:)
      INTEGER :: istep, j

#ifdef _CRAYFTN
      INTEGER,  POINTER, DIMENSION(:)      :: bdy_halo_c_blk, bdy_halo_c_idx, bdy_mflx_e_blk, bdy_mflx_e_idx,     &
                                              nudge_e_blk, nudge_e_idx, pg_vertidx, pg_edgeidx, pg_edgeblk
      INTEGER,  POINTER, DIMENSION(:,:)    :: zd_blklist, zd_edgeblk, zd_edgeidx, zd_indlist, zd_vertidx
      INTEGER,  POINTER, DIMENSION(:,:,:,:):: vertidx_gradp
      REAL(wp), POINTER, DIMENSION(:)      :: enhfac_diffu, rayleigh_vn, rayleigh_w, scalfac_dd3d
      REAL(wp), POINTER, DIMENSION(:,:)    :: pres_sfc, vwind_expl_wgt, vwind_impl_wgt,                     &
                                              hmask_dd3d, dvn_ie_ubc, dw_ubc, dtheta_v_ic_ubc,              &
                                              zd_intcoef, zd_geofac, zd_e2cell
      REAL(wp), POINTER, DIMENSION(:,:,:)  :: inv_ddqz_z_full, vn_ref, w_ref, rho_ic, theta_v_ic,           &
                                              dwdx, dwdy, mass_fl_e, div_ic, hdef_ic, ddqz_z_full,          &
                                              exner_pr, mflx_ic_ubc, theta_v_ic_ubc, grf_bdy_mflx,          &
                                              dpres_mc, pres_ifc, pres, temp, tempv, temp_ifc
      REAL(wp), POINTER, DIMENSION(:,:,:,:):: tracer
      REAL(vp), POINTER, DIMENSION(:)      :: pg_exdist
      REAL(vp), POINTER, DIMENSION(:,:,:)  :: d_exner_dz_ref_ic, d2dexdz2_fac1_mc, d2dexdz2_fac2_mc,        &
                                              ddqz_z_half, ddxn_z_full, ddxt_z_full, ddqz_z_full_e,         &
                                              exner_exfac, exner_ref_mc,                                    &
                                              rho_ref_mc, rho_ref_me, theta_ref_ic, theta_ref_mc,           &
                                              theta_ref_me, vt, v, u, vor, vn_ie, w_concorr_c,              &
                                              wgtfac_c, wgtfac_e, wgtfacq_c,                                &
                                              wgtfacq_e, wgtfacq1_c, wgtfacq1_e,                            &
                                              exner_dyn_incr, ddt_exner_phy, ddt_vn_phy,                    &
                                              rho_incr, exner_incr, coeff1_dwdz, coeff2_dwdz, coeff_gradekin
                             
      REAL(vp), POINTER, DIMENSION(:,:,:,:):: coeff_gradp, zdiff_gradp, ddt_vn_adv, ddt_w_adv, ddt_tracer_adv


      REAL(wp), POINTER, DIMENSION(:,:,:)  :: exner, rho, theta_v, vn, w

      LOGICAL,  POINTER                    :: mask_prog_halo_c(:,:)
#endif

      DO j = 1, SIZE(p_nh)

#ifdef _CRAYFTN

        bdy_halo_c_blk              => p_nh(j)%metrics%bdy_halo_c_blk
        bdy_halo_c_idx              => p_nh(j)%metrics%bdy_halo_c_idx
        bdy_mflx_e_blk              => p_nh(j)%metrics%bdy_mflx_e_blk
        bdy_mflx_e_idx              => p_nh(j)%metrics%bdy_mflx_e_idx
        coeff1_dwdz                 => p_nh(j)%metrics%coeff1_dwdz
        coeff2_dwdz                 => p_nh(j)%metrics%coeff2_dwdz
        coeff_gradekin              => p_nh(j)%metrics%coeff_gradekin
        coeff_gradp                 => p_nh(j)%metrics%coeff_gradp
        d_exner_dz_ref_ic           => p_nh(j)%metrics%d_exner_dz_ref_ic
        d2dexdz2_fac1_mc            => p_nh(j)%metrics%d2dexdz2_fac1_mc
        d2dexdz2_fac2_mc            => p_nh(j)%metrics%d2dexdz2_fac2_mc
        ddqz_z_full                 => p_nh(j)%metrics%ddqz_z_full
        ddqz_z_full_e               => p_nh(j)%metrics%ddqz_z_full_e
        ddqz_z_half                 => p_nh(j)%metrics%ddqz_z_half
        ddxn_z_full                 => p_nh(j)%metrics%ddxn_z_full
        ddxt_z_full                 => p_nh(j)%metrics%ddxt_z_full
        enhfac_diffu                => p_nh(j)%metrics%enhfac_diffu
        exner_exfac                 => p_nh(j)%metrics%exner_exfac
        exner_ref_mc                => p_nh(j)%metrics%exner_ref_mc
        hmask_dd3d                  => p_nh(j)%metrics%hmask_dd3d
        inv_ddqz_z_full             => p_nh(j)%metrics%inv_ddqz_z_full
        mask_prog_halo_c            => p_nh(j)%metrics%mask_prog_halo_c
        nudge_e_blk                 => p_nh(j)%metrics%nudge_e_blk
        nudge_e_idx                 => p_nh(j)%metrics%nudge_e_idx
        pg_exdist                   => p_nh(j)%metrics%pg_exdist
        pg_vertidx                  => p_nh(j)%metrics%pg_vertidx
        pg_edgeidx                  => p_nh(j)%metrics%pg_edgeidx
        pg_edgeblk                  => p_nh(j)%metrics%pg_edgeblk
        rayleigh_vn                 => p_nh(j)%metrics%rayleigh_vn
        rayleigh_w                  => p_nh(j)%metrics%rayleigh_w
        rho_ref_mc                  => p_nh(j)%metrics%rho_ref_mc
        rho_ref_me                  => p_nh(j)%metrics%rho_ref_me
        scalfac_dd3d                => p_nh(j)%metrics%scalfac_dd3d
        theta_ref_ic                => p_nh(j)%metrics%theta_ref_ic
        theta_ref_mc                => p_nh(j)%metrics%theta_ref_mc
        theta_ref_me                => p_nh(j)%metrics%theta_ref_me
        vertidx_gradp               => p_nh(j)%metrics%vertidx_gradp
        vwind_expl_wgt              => p_nh(j)%metrics%vwind_expl_wgt
        vwind_impl_wgt              => p_nh(j)%metrics%vwind_impl_wgt
        wgtfac_c                    => p_nh(j)%metrics%wgtfac_c
        wgtfac_e                    => p_nh(j)%metrics%wgtfac_e
        wgtfacq_c                   => p_nh(j)%metrics%wgtfacq_c
        wgtfacq_e                   => p_nh(j)%metrics%wgtfacq_e
        wgtfacq1_c                  => p_nh(j)%metrics%wgtfacq1_c
        wgtfacq1_e                  => p_nh(j)%metrics%wgtfacq1_e
        zdiff_gradp                 => p_nh(j)%metrics%zdiff_gradp
        zd_indlist                  => p_nh(j)%metrics%zd_indlist
        zd_blklist                  => p_nh(j)%metrics%zd_blklist
        zd_vertidx                  => p_nh(j)%metrics%zd_vertidx
        zd_edgeidx                  => p_nh(j)%metrics%zd_edgeidx
        zd_edgeblk                  => p_nh(j)%metrics%zd_edgeblk
        zd_intcoef                  => p_nh(j)%metrics%zd_intcoef
        zd_e2cell                   => p_nh(j)%metrics%zd_e2cell
        zd_geofac                   => p_nh(j)%metrics%zd_geofac

        vn_ref                      => p_nh(j)%ref%vn_ref
        w_ref                       => p_nh(j)%ref%w_ref

        dpres_mc                    => p_nh(j)%diag%dpres_mc
        ddt_exner_phy               => p_nh(j)%diag%ddt_exner_phy
        ddt_vn_adv                  => p_nh(j)%diag%ddt_vn_adv
        ddt_vn_phy                  => p_nh(j)%diag%ddt_vn_phy
        ddt_w_adv                   => p_nh(j)%diag%ddt_w_adv
        div_ic                      => p_nh(j)%diag%div_ic
        dwdx                        => p_nh(j)%diag%dwdx
        dwdy                        => p_nh(j)%diag%dwdy
        dtheta_v_ic_ubc             => p_nh(j)%diag%dtheta_v_ic_ubc
        dvn_ie_ubc                  => p_nh(j)%diag%dvn_ie_ubc
        dw_ubc                      => p_nh(j)%diag%dw_ubc
        exner_dyn_incr              => p_nh(j)%diag%exner_dyn_incr
        exner_incr                  => p_nh(j)%diag%exner_incr
        exner_pr                    => p_nh(j)%diag%exner_pr
        grf_bdy_mflx                => p_nh(j)%diag%grf_bdy_mflx
        hdef_ic                     => p_nh(j)%diag%hdef_ic
        mass_fl_e                   => p_nh(j)%diag%mass_fl_e
        mflx_ic_ubc                 => p_nh(j)%diag%mflx_ic_ubc
        pres_ifc                    => p_nh(j)%diag%pres_ifc
        pres_sfc                    => p_nh(j)%diag%pres_sfc
        pres                        => p_nh(j)%diag%pres
        rho_ic                      => p_nh(j)%diag%rho_ic
        rho_incr                    => p_nh(j)%diag%rho_incr
        temp                        => p_nh(j)%diag%temp
        tempv                       => p_nh(j)%diag%tempv
        temp_ifc                    => p_nh(j)%diag%temp_ifc
        theta_v_ic                  => p_nh(j)%diag%theta_v_ic
        vn_ie                       => p_nh(j)%diag%vn_ie
        vt                          => p_nh(j)%diag%vt
        v                           => p_nh(j)%diag%v
        u                           => p_nh(j)%diag%u
        vor                         => p_nh(j)%diag%vor
        w_concorr_c                 => p_nh(j)%diag%w_concorr_c
        ddt_tracer_adv              => p_nh(j)%diag%ddt_tracer_adv

!$ACC EXIT DATA &
!$ACC      DELETE( bdy_halo_c_blk, bdy_halo_c_idx, bdy_mflx_e_blk, bdy_mflx_e_idx, &
!$ACC            coeff1_dwdz, coeff2_dwdz, coeff_gradekin, coeff_gradp,            &
!$ACC            d_exner_dz_ref_ic, d2dexdz2_fac1_mc, d2dexdz2_fac2_mc,            &
!$ACC            ddqz_z_full, ddqz_z_full_e, ddqz_z_half, ddxn_z_full, ddxt_z_full,&
!$ACC            enhfac_diffu, exner_exfac, exner_ref_mc, hmask_dd3d,              &
!$ACC            inv_ddqz_z_full, mask_prog_halo_c, nudge_e_blk, nudge_e_idx,      &
!$ACC            pg_exdist, pg_vertidx, pg_edgeidx, pg_edgeblk,                    &
!$ACC            rayleigh_vn, rayleigh_w, rho_ref_mc, rho_ref_me,                  &
!$ACC            scalfac_dd3d, theta_ref_ic, theta_ref_mc, theta_ref_me,           &
!$ACC            vertidx_gradp, vwind_expl_wgt, vwind_impl_wgt,                    &
!$ACC            wgtfac_c, wgtfac_e, wgtfacq_c, wgtfacq_e, wgtfacq1_c, wgtfacq1_e, &
!$ACC            zdiff_gradp, zd_blklist, zd_e2cell, zd_edgeblk, zd_edgeidx,       &
!$ACC            zd_geofac, zd_indlist, zd_intcoef, zd_vertidx        ),           &
!$ACC      DELETE( vn_ref, w_ref ),                                                &
!$ACC      DELETE( ddt_exner_phy, ddt_vn_adv, ddt_vn_phy, ddt_w_adv, div_ic,       &
!$ACC              dpres_mc, dtheta_v_ic_ubc, dvn_ie_ubc, dw_ubc, dwdx, dwdy,      &
!$ACC              exner_dyn_incr, exner_incr, exner_pr, grf_bdy_mflx,             &
!$ACC              hdef_ic, mass_fl_e, mflx_ic_ubc, pres_ifc, pres, pres_sfc,      &
!$ACC              rho_ic, rho_incr, temp, tempv, temp_ifc, theta_v_ic, vt, vn_ie, &
!$ACC              u, v, vor, w_concorr_c, ddt_tracer_adv ),                       &
!$ACC      IF ( i_am_accel_node )

        DO istep = 1, SIZE( p_nh(j)%prog )

          exner   =>  p_nh(j)%prog(istep)%exner
          rho     =>  p_nh(j)%prog(istep)%rho
          theta_v =>  p_nh(j)%prog(istep)%theta_v
          vn      =>  p_nh(j)%prog(istep)%vn
          w       =>  p_nh(j)%prog(istep)%w
!$ACC EXIT DATA DELETE( exner, rho, theta_v, vn, w ), IF ( i_am_accel_node  )

        ENDDO

#else

!$ACC EXIT DATA &
!$ACC      DELETE( p_nh(j)%metrics%bdy_halo_c_blk, p_nh(j)%metrics%bdy_halo_c_idx,                   &
!$ACC              p_nh(j)%metrics%bdy_mflx_e_blk, p_nh(j)%metrics%bdy_mflx_e_idx,                   &
!$ACC              p_nh(j)%metrics%coeff1_dwdz, p_nh(j)%metrics%coeff2_dwdz,                         &
!$ACC              p_nh(j)%metrics%coeff_gradekin, p_nh(j)%metrics%coeff_gradp,                      &
!$ACC              p_nh(j)%metrics%d_exner_dz_ref_ic, p_nh(j)%metrics%d2dexdz2_fac1_mc,              &
!$ACC              p_nh(j)%metrics%d2dexdz2_fac2_mc, p_nh(j)%metrics%ddqz_z_full,                    &
!$ACC              p_nh(j)%metrics%ddqz_z_full_e, p_nh(j)%metrics%ddqz_z_half,                       &
!$ACC              p_nh(j)%metrics%ddxn_z_full, p_nh(j)%metrics%ddxt_z_full,                         &
!$ACC              p_nh(j)%metrics%exner_ref_mc,                                                     &
!$ACC              p_nh(j)%metrics%hmask_dd3d, p_nh(j)%metrics%inv_ddqz_z_full,                      &
!$ACC              p_nh(j)%metrics%mask_prog_halo_c, p_nh(j)%metrics%nudge_e_blk,                    &
!$ACC              p_nh(j)%metrics%nudge_e_idx, p_nh(j)%metrics%pg_exdist,                           &
!$ACC              p_nh(j)%metrics%pg_vertidx, p_nh(j)%metrics%pg_edgeidx,                           &
!$ACC              p_nh(j)%metrics%pg_edgeblk, p_nh(j)%metrics%rayleigh_vn,                          &
!$ACC              p_nh(j)%metrics%rayleigh_w, p_nh(j)%metrics%rho_ref_mc,                           &
!$ACC              p_nh(j)%metrics%rho_ref_me, p_nh(j)%metrics%scalfac_dd3d,                         &
!$ACC              p_nh(j)%metrics%theta_ref_ic, p_nh(j)%metrics%theta_ref_mc,                       &
!$ACC              p_nh(j)%metrics%theta_ref_me, p_nh(j)%metrics%vertidx_gradp,                      &
!$ACC              p_nh(j)%metrics%vwind_expl_wgt, p_nh(j)%metrics%vwind_impl_wgt,                   &
!$ACC              p_nh(j)%metrics%wgtfac_c, p_nh(j)%metrics%wgtfac_e,                               &
!$ACC              p_nh(j)%metrics%wgtfacq_c, p_nh(j)%metrics%wgtfacq_e,                             &
!$ACC              p_nh(j)%metrics%wgtfacq1_c, p_nh(j)%metrics%wgtfacq1_e,                           &
!$ACC              p_nh(j)%metrics%zdiff_gradp, p_nh(j)%metrics%zd_blklist,                          &
!$ACC              p_nh(j)%metrics%zd_e2cell, p_nh(j)%metrics%zd_edgeblk,                            &
!$ACC              p_nh(j)%metrics%zd_edgeidx, p_nh(j)%metrics%zd_geofac, p_nh(j)%metrics%zd_indlist,&
!$ACC              p_nh(j)%metrics%zd_intcoef,p_nh(j)%metrics%zd_vertidx, p_nh(j)%metrics ),         &
!$ACC      DELETE( p_nh(j)%ref%vn_ref, p_nh(j)%ref%w_ref, p_nh(j)%ref ),                             &
!$ACC      DELETE( p_nh(j)%diag%ddt_exner_phy, p_nh(j)%diag%ddt_vn_adv,                              &
!$ACC              p_nh(j)%diag%ddt_vn_phy, p_nh(j)%diag%ddt_w_adv,                                  &
!$ACC              p_nh(j)%diag%div_ic, p_nh(j)%diag%dpres_mc,                                       &
!$ACC              p_nh(j)%diag%dtheta_v_ic_ubc, p_nh(j)%diag%dvn_ie_ubc,                            &
!$ACC              p_nh(j)%diag%dw_ubc, p_nh(j)%diag%dwdx, p_nh(j)%diag%dwdy,                        &
!$ACC              p_nh(j)%diag%exner_dyn_incr, p_nh(j)%diag%exner_incr, p_nh(j)%diag%exner_pr,      &
!$ACC              p_nh(j)%diag%grf_bdy_mflx, p_nh(j)%diag%hdef_ic,                                  &
!$ACC              p_nh(j)%diag%mass_fl_e, p_nh(j)%diag%mflx_ic_ubc, p_nh(j)%diag%pres_ifc,          &
!$ACC              p_nh(j)%diag%pres_sfc, p_nh(j)%diag%pres, p_nh(j)%diag%rho_ic,                    &
!$ACC              p_nh(j)%diag%rho_incr, p_nh(j)%diag%temp, p_nh(j)%diag%tempv,                     &
!$ACC              p_nh(j)%diag%temp_ifc, p_nh(j)%diag%theta_v_ic, p_nh(j)%diag%vt,                  &
!$ACC              p_nh(j)%diag%vn_ie, p_nh(j)%diag%v, p_nh(j)%diag%u, p_nh(j)%diag%vor,             &
!$ACC              p_nh(j)%diag%w_concorr_c, p_nh(j)%diag%ddt_tracer_adv ),                          &
!$ACC      IF ( i_am_accel_node )

        DO istep = 1, SIZE( p_nh(j)%prog )

!$ACC EXIT DATA DELETE( p_nh(j)%prog(istep)%exner, p_nh(j)%prog(istep)%rho, p_nh(j)%prog(istep)%theta_v, &
!$ACC                   p_nh(j)%prog%tracer, p_nh(j)%prog(istep)%vn, p_nh(j)%prog(istep)%w ),            &
!$ACC      IF ( i_am_accel_node  )

        ENDDO

!$ACC EXIT DATA DELETE( p_nh(j)%prog ), IF ( i_am_accel_node  )

#endif

      ENDDO

      END SUBROUTINE d2h_nh_state

    END SUBROUTINE d2h_icon

#endif

END MODULE mo_nonhydro_gpu_types
