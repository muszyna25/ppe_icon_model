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
  USE mo_impl_constants,       ONLY: MAX_CHAR_LENGTH, inwp, iecham
  USE mo_mpi,                  ONLY: i_am_accel_node
  USE mo_fortran_tools,        ONLY: t_ptr_2d3d
  USE mo_math_types,           ONLY: t_geographical_coordinates
  USE mo_model_domain,         ONLY: t_patch, t_tangent_vectors
  USE mo_nonhydro_types,       ONLY: t_nh_state, t_nh_diag, t_nh_prog
  USE mo_nh_prepadv_types,     ONLY: t_prepare_adv
  USE mo_advection_config,     ONLY: t_advection_config
  USE mo_intp_data_strc,       ONLY: t_int_state
  USE mo_grf_intp_data_strc,   ONLY: t_gridref_single_state, t_gridref_state
  USE mo_var_list_gpu,         ONLY: gpu_h2d_var_list, gpu_d2h_var_list
  USE mo_run_config,           ONLY: ltestcase

  IMPLICIT NONE

  PRIVATE 

  PUBLIC :: h2d_icon, d2h_icon, devcpy_grf_state

CONTAINS

  SUBROUTINE h2d_icon( p_int_states, p_patches, p_nh_states, prep_advs, advection_config, iforcing )

    TYPE ( t_int_state ),       INTENT(INOUT) :: p_int_states(:)
    TYPE ( t_patch ),           INTENT(INOUT) :: p_patches(:)
    TYPE ( t_nh_state ),        INTENT(INOUT) :: p_nh_states(:)
    TYPE ( t_prepare_adv),      INTENT(INOUT) :: prep_advs(:)
    TYPE ( t_advection_config), INTENT(INOUT) :: advection_config(:)
    INTEGER, INTENT(IN)                       :: iforcing 
    INTEGER :: jg
!
! Copy all data need on GPU from host to device
!

!$ACC ENTER DATA COPYIN( p_int_states, p_patches, prep_advs, advection_config ), IF ( i_am_accel_node  )

    CALL transfer_int_state( p_int_states, .TRUE. )

    CALL transfer_patch( p_patches, .TRUE. )

    CALL transfer_prep_adv( prep_advs, .TRUE. )

    CALL transfer_nh_state( p_nh_states, .TRUE. )

    CALL transfer_advection_config( advection_config, .TRUE. )

    IF( iforcing == iecham ) THEN
      CALL transfer_echam( p_patches, .TRUE. )
    END IF

  END SUBROUTINE h2d_icon

  SUBROUTINE d2h_icon( p_int_states, p_patches, p_nh_states, prep_advs, advection_config, iforcing )

    TYPE ( t_int_state ),  INTENT(INOUT)      :: p_int_states(:)
    TYPE ( t_patch ),      INTENT(INOUT)      :: p_patches(:)
    TYPE ( t_nh_state ),   INTENT(INOUT)      :: p_nh_states(:)
    TYPE ( t_prepare_adv), INTENT(INOUT)      :: prep_advs(:)
    TYPE ( t_advection_config), INTENT(INOUT) :: advection_config(:)
    INTEGER, INTENT(IN)                       :: iforcing 

!
! Delete all data on GPU
!
    CALL transfer_nh_state( p_nh_states, .FALSE. )
    CALL transfer_prep_adv( prep_advs, .FALSE. )
    CALL transfer_patch( p_patches, .FALSE. )
    CALL transfer_int_state( p_int_states, .FALSE. )
    CALL transfer_advection_config( advection_config, .FALSE. )

    IF( iforcing == iecham ) THEN
      CALL transfer_echam( p_patches, .FALSE. )
    END IF

!$ACC EXIT DATA DELETE( advection_config, prep_advs, p_patches, p_int_states ), IF ( i_am_accel_node  )

  END SUBROUTINE d2h_icon

  SUBROUTINE transfer_int_state( p_int, host_to_device )

    LOGICAL, INTENT(IN)                        :: host_to_device     !   .TRUE. : h2d   .FALSE. : d2h
    TYPE ( t_int_state ), TARGET,  INTENT(INOUT) :: p_int(:)

    INTEGER  :: j

    DO j=1, SIZE(p_int)

      IF ( host_to_device ) THEN

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
!$ACC               p_int(j)%rbf_vec_blk_v, p_int(j)%rbf_vec_idx_v, p_int(j)%rbf_vec_coeff_v,    &
!$ACC               p_int(j)%verts_aw_cells ),                                                   &
!$ACC       IF ( i_am_accel_node )        

      ELSE

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
!$ACC               p_int(j)%verts_aw_cells, p_int(j)%lsq_high, p_int(j)%lsq_lin )               &
!$ACC       IF ( i_am_accel_node )        

      ENDIF

    ENDDO

  END SUBROUTINE transfer_int_state


  SUBROUTINE transfer_patch( p_patch, host_to_device )

    LOGICAL, INTENT(IN)                        :: host_to_device     !   .TRUE. : h2d   .FALSE. : d2h
    TYPE ( t_patch ), TARGET, INTENT(INOUT)    :: p_patch(:)

    INTEGER :: j

!
! Copy the static data structures in p_patch to the device -- this is a small subset of all the components
! The communication patterns are copied over in mo_communication_orig.
!

      DO j=1,SIZE(p_patch)

        IF ( host_to_device ) THEN

!$ACC ENTER DATA &
!$ACC      COPYIN( p_patch(j)%cells, p_patch(j)%cells%decomp_info, p_patch(j)%cells%decomp_info%owner_mask, &
!$ACC              p_patch(j)%cells%area, p_patch(j)%cells%edge_idx, p_patch(j)%cells%edge_blk,             &
!$ACC              p_patch(j)%cells%neighbor_idx, p_patch(j)%cells%neighbor_blk,                            &
!$ACC              p_patch(j)%cells%center, p_patch(j)%cells%refin_ctrl,                                    &
!$ACC              p_patch(j)%cells%vertex_blk, p_patch(j)%cells%vertex_idx,                                &
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
     
        ELSE

!$ACC EXIT DATA &
!$ACC      DELETE( p_patch(j)%cells%decomp_info, p_patch(j)%cells%decomp_info%owner_mask,               &
!$ACC              p_patch(j)%cells%area, p_patch(j)%cells%edge_idx, p_patch(j)%cells%edge_blk,         &
!$ACC              p_patch(j)%cells%neighbor_idx, p_patch(j)%cells%neighbor_blk,                        &
!$ACC              p_patch(j)%cells%center, p_patch(j)%cells%refin_ctrl, p_patch(j)%cells,              &
!$ACC              p_patch(j)%cells%vertex_blk, p_patch(j)%cells%vertex_idx,                            &
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

      ENDIF   

    ENDDO

  END SUBROUTINE transfer_patch



  SUBROUTINE transfer_prep_adv( prep_adv, host_to_device )

    LOGICAL, INTENT(IN)                        :: host_to_device     !   .TRUE. : h2d   .FALSE. : d2h
    TYPE ( t_prepare_adv ), TARGET, INTENT(INOUT) :: prep_adv(:)

    INTEGER :: j

    DO j=1, SIZE(prep_adv)

      IF ( host_to_device ) THEN

!$ACC ENTER DATA COPYIN(prep_adv(j)%vn_traj,prep_adv(j)%mass_flx_me,prep_adv(j)%mass_flx_ic,prep_adv(j)%topflx_tra ), &
!$ACC       IF ( i_am_accel_node  )

      ELSE

!$ACC EXIT DATA DELETE(prep_adv(j)%vn_traj,prep_adv(j)%mass_flx_me,prep_adv(j)%mass_flx_ic,prep_adv(j)%topflx_tra ), &
!$ACC      IF ( i_am_accel_node  )

      ENDIF

    ENDDO    

  END SUBROUTINE transfer_prep_adv

  SUBROUTINE transfer_advection_config( advection_config, host_to_device )

    LOGICAL, INTENT(IN)                        :: host_to_device     !   .TRUE. : h2d   .FALSE. : d2h
    TYPE ( t_advection_config ), TARGET, INTENT(INOUT) :: advection_config(:)

    INTEGER :: j

    DO j=1, SIZE(advection_config)

      IF ( host_to_device ) THEN

!$ACC ENTER DATA COPYIN( advection_config(j)%trHydroMass%list ) &
!$ACC            IF ( i_am_accel_node  )

      ELSE

!$ACC EXIT DATA DELETE( advection_config(j)%trHydroMass%list )  &
!$ACC           IF ( i_am_accel_node  )

      ENDIF

    ENDDO

  END SUBROUTINE transfer_advection_config

  SUBROUTINE transfer_nh_state( p_nh, host_to_device )

    LOGICAL, INTENT(IN)                        :: host_to_device     !   .TRUE. : h2d   .FALSE. : d2h
    TYPE ( t_nh_state ), TARGET, INTENT(INOUT) :: p_nh(:)
    INTEGER :: istep, jg
    CHARACTER(len=MAX_CHAR_LENGTH) :: metrics, prog, ref, diag

!
! At this point, p_nh and all its underlying subtypes have been created on the device
!
    DO jg = 1, SIZE(p_nh)

      WRITE(metrics,'(a)') 'nh_state_metrics_of_domain_'
      WRITE(diag,'(a)') 'nh_state_diag_of_domain_'
      WRITE(ref,'(a)') 'nh_state_ref_of_domain_'
      WRITE(prog,'(a)') 'nh_state_prog_of_domain_'

      IF ( host_to_device ) THEN

        CALL gpu_h2d_var_list(TRIM(metrics), domain=jg )
        IF (ltestcase) CALL gpu_h2d_var_list(TRIM(ref), domain=jg )
        CALL gpu_h2d_var_list(TRIM(diag), domain=jg )
        DO istep = 1, SIZE(p_nh(jg)%prog)
          CALL gpu_h2d_var_list(TRIM(prog), domain=jg, substr='_and_timelev_', timelev=istep )
        ENDDO

      ELSE

! 
! WS:  currently it appears to be unnecessary to update any of these values back to the host 
!      after the end of the time loop.  Dycore variables are updated in ACC_VALIDATE mode individually
!      BUT: there should be a way to delete all the variables with DEL_VAR
!
#if 0
        CALL gpu_d2h_var_list(TRIM(metrics), domain=jg )
        CALL gpu_d2h_var_list(TRIM(diag), domain=jg )
        CALL gpu_d2h_var_list(TRIM(ref), domain=jg )

        DO istep = 1, SIZE(p_nh(jg)%prog)
          WRITE(prog,'(a,i2.2,a,i2.2)') 'nh_state_prog_of_domain_',jg, '_and_timelev_'
          CALL gpu_d2h_var_list(TRIM(prog), timelev=istep )
        ENDDO
#endif

      ENDIF

    ENDDO

  END SUBROUTINE transfer_nh_state

  SUBROUTINE transfer_echam( p_patches, host_to_device )
    TYPE ( t_patch ),      INTENT(INOUT) :: p_patches(:)
    LOGICAL, INTENT(IN)                  :: host_to_device     !   .TRUE. : h2d   .FALSE. : d2h
    INTEGER :: jg

    DO jg = 1, SIZE(p_patches)
      IF( host_to_device ) THEN
        CALL gpu_h2d_var_list('prm_field_D', domain=jg)
        CALL gpu_h2d_var_list('prm_tend_D', domain=jg)
      ELSE
        CALL gpu_d2h_var_list('prm_field_D', domain=jg)
        CALL gpu_d2h_var_list('prm_tend_D', domain=jg)
      END IF
    END DO

  END SUBROUTINE transfer_echam

  SUBROUTINE devcpy_grf_state( p_grf, l_h2d )

      TYPE ( t_gridref_state ), TARGET,  INTENT(INOUT) :: p_grf(:)
      LOGICAL, INTENT(IN) :: l_h2d    ! true host-to-device, false device-to-host

      INTEGER  :: j,k

#ifndef _CRAYFTN
!$ACC ENTER DATA COPYIN( p_grf ), IF ( i_am_accel_node  )
#endif

      DO j=1, SIZE(p_grf)

        IF (l_h2d) THEN

!$ACC ENTER DATA &
!$ACC       COPYIN( p_grf(j)%fbk_wgt_aw, p_grf(j)%fbk_wgt_bln, p_grf(j)%fbk_wgt_e, p_grf(j)%fbk_dom_area,  &
!$ACC               p_grf(j)%mask_ovlp_c, p_grf(j)%mask_ovlp_ch, p_grf(j)%mask_ovlp_e, p_grf(j)%mask_ovlp_v,        &
!$ACC               p_grf(j)%idxlist_bdyintp_src_c, p_grf(j)%idxlist_bdyintp_src_e, p_grf(j)%blklist_bdyintp_src_c, &
!$ACC               p_grf(j)%blklist_bdyintp_src_e,p_grf(j)%p_dom ),  &
!$ACC       IF ( i_am_accel_node )        

        DO k = 1, SIZE(p_grf(j)%p_dom)
          CALL devcpy_grf_single_state( p_grf(j)%p_dom, l_h2d )
        ENDDO

        ELSE

        DO k = 1, SIZE(p_grf(j)%p_dom)
          CALL devcpy_grf_single_state( p_grf(j)%p_dom, l_h2d )
        ENDDO

!$ACC EXIT DATA &
!$ACC      DELETE(  p_grf(j)%fbk_wgt_aw, p_grf(j)%fbk_wgt_bln, p_grf(j)%fbk_wgt_e, p_grf(j)%fbk_dom_area,           &
!$ACC               p_grf(j)%mask_ovlp_c, p_grf(j)%mask_ovlp_ch, p_grf(j)%mask_ovlp_e, p_grf(j)%mask_ovlp_v,        &
!$ACC               p_grf(j)%idxlist_bdyintp_src_c, p_grf(j)%idxlist_bdyintp_src_e, p_grf(j)%blklist_bdyintp_src_c, &
!$ACC               p_grf(j)%blklist_bdyintp_src_e,p_grf(j)%p_dom )                                        &
!$ACC       IF ( i_am_accel_node )        

#ifndef _CRAYFTN
!$ACC EXIT DATA DELETE( p_grf ), IF ( i_am_accel_node  )
#endif

        ENDIF

      ENDDO

    END SUBROUTINE devcpy_grf_state

    SUBROUTINE devcpy_grf_single_state( p_grf, l_h2d )

      TYPE ( t_gridref_single_state ), TARGET,  INTENT(INOUT) :: p_grf(:)
      LOGICAL, INTENT(IN) :: l_h2d    ! true host-to-device, false device-to-host

      INTEGER  :: j


      DO j=1, SIZE(p_grf)

        IF (l_h2d) THEN

!$ACC ENTER DATA &
!$ACC       COPYIN( p_grf(j)%grf_dist_pc2cc, p_grf(j)%grf_dist_pe2ce, p_grf(j)%idxlist_bdyintp_c,                  &
!$ACC       p_grf(j)%idxlist_bdyintp_e, p_grf(j)%idxlist_ubcintp_c, p_grf(j)%idxlist_ubcintp_e, p_grf(j)%blklist_bdyintp_c, &
!$ACC       p_grf(j)%blklist_bdyintp_e, p_grf(j)%blklist_ubcintp_c, p_grf(j)%blklist_ubcintp_e, p_grf(j)%idxlist_rbfintp_v, &
!$ACC       p_grf(j)%blklist_rbfintp_v, p_grf(j)%edge_vert_idx, p_grf(j)%coeff_bdyintp_c, p_grf(j)%coeff_ubcintp_c,         &
!$ACC       p_grf(j)%dist_pc2cc_bdy, p_grf(j)%dist_pc2cc_ubc, p_grf(j)%prim_norm, p_grf(j)%coeff_bdyintp_e12,               &
!$ACC       p_grf(j)%coeff_bdyintp_e34, p_grf(j)%dist_pe2ce, p_grf(j)%coeff_ubcintp_e12, p_grf(j)%coeff_ubcintp_e34,        &
!$ACC       p_grf(j)%coeff_rbf_v ),    IF ( i_am_accel_node )        

        ELSE

!$ACC EXIT DATA &
!$ACC      DELETE(  p_grf(j)%grf_dist_pc2cc, p_grf(j)%grf_dist_pe2ce, p_grf(j)%idxlist_bdyintp_c,                           &
!$ACC       p_grf(j)%idxlist_bdyintp_e, p_grf(j)%idxlist_ubcintp_c, p_grf(j)%idxlist_ubcintp_e, p_grf(j)%blklist_bdyintp_c, &
!$ACC       p_grf(j)%blklist_bdyintp_e, p_grf(j)%blklist_ubcintp_c, p_grf(j)%blklist_ubcintp_e, p_grf(j)%idxlist_rbfintp_v, &
!$ACC       p_grf(j)%blklist_rbfintp_v, p_grf(j)%edge_vert_idx, p_grf(j)%coeff_bdyintp_c, p_grf(j)%coeff_ubcintp_c,         &
!$ACC       p_grf(j)%dist_pc2cc_bdy, p_grf(j)%dist_pc2cc_ubc, p_grf(j)%prim_norm, p_grf(j)%coeff_bdyintp_e12,               &
!$ACC       p_grf(j)%coeff_bdyintp_e34, p_grf(j)%dist_pe2ce, p_grf(j)%coeff_ubcintp_e12, p_grf(j)%coeff_ubcintp_e34,        &
!$ACC       p_grf(j)%coeff_rbf_v )                                        &
!$ACC       IF ( i_am_accel_node )        

        ENDIF

      ENDDO

    END SUBROUTINE devcpy_grf_single_state

#endif



END MODULE mo_nonhydro_gpu_types
