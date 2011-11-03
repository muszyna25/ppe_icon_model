
#ifdef __xlC__
! @PROCESS nosmp
! @PROCESS NOOPTimize
@PROCESS smp=noopt
@PROCESS noopt
#endif
#ifdef __PGI
!pgi$g opt=1
#endif
!>
!! Contains the implementation of interpolation and reconstruction.
!!
!! Contains the implementation of interpolation and reconstruction
!! routines used by the shallow water model, including the RBF
!! reconstruction routines.
!!
!! @par Revision History
!! Developed  by Luca Bonaventura and Will Sawyer (2002-4).
!! Modified to ProTeX-style by  Luca Bonaventura and Thomas Heinze (2004).
!! Adapted to new data structure by Thomas Heinze,
!! Peter Korn and Luca Bonaventura (2005).
!! Modification by Thomas Heinze (2006-02-21):
!! - renamed m_modules to mo_modules
!! Modification by Thomas Heinze (2006-07-05):
!! - modified cell2edge_lin_int_coeff
!! - created cc_dot_product
!! Modification by Peter Korn and Luca Bonaventura(2006-07-28):
!! - moved several auxiliary functions to mo_math_utilities
!! - introduced recoded rbf interpolation for vector fields
!! - added lraviart switch to force RT interpolation to be used
!! Modification by Thomas Heinze  and Luca Bonaventura(2006-10-05):
!! - merged with 'Milano' version by P. Korn
!! Modification by Pilar Ripodas (2006-11):
!! - new subroutine rbf_vec_interpol_car with the cartesian
!!   coordinates as output
!! Modification by Peter Korn, MPI-M, (2006-11-23):
!! - replacements in TYPE patch: ic by l2g_c, ie by l2g_e, iv by l2g_v,
!!   iic by g2l_c, iie by g2l_e, iiv by g2l_v
!! - replaced edge_index by edge_idx
!! - replaced vertex_index by vertex_idx
!! - replaced cell_index by cell_idx
!! - replaced neighbor_index by neighbor_idx
!! Modification by Pilar Ripodas (2006-12):
!! - dt_tan_vec and dt_tan_rt_vec are wrong. They are renamed to
!!   dt_tan_vec_old and dt_tan_rt_vec_old and should not be used
!! - New subroutines dt_tan_vec_h and dt_tan_vec_kin and
!!   dt_tan_vec_gen are produced and
!!   moved to mo_sw_state.f90
!!  Modification by Peter Korn, MPI-M (2007-02)
!!  Modification by Hui Wan, MPI-M (2007-02-22)
!!  - changes in the USE section because
!!    the coordinate types had been move from mo_model_domain
!!    to mo_math_utilities;
!!  Modification by Almut Gassmann, MPI-M (2007-04)
!!  - removed reference to unused halo_verts
!!  - summing over all halos of the various parallel patches (Quick and Dirty!)
!!  Modification by Almut Gassmann, MPI-M (2007-04)
!!  - abandon grid for the sake of patch
!!  Modification by Thomas Heinze, DWD (2007-07-26)
!!  - including all the improvements of Tobias Ruppert's diploma thesis
!!  - several changes according to the programming guide
!!  Modification by Pilar Ripodas, DWD (2007-07):
!!  - substruct the outgoing component of the reconstructed
!!    vector in subroutine "rbf_vec_interpol_car"
!!  Modification by Thomas Heinze, DWD (2007-08-02)
!!  - replaced rbf_kern_dim by rbf_kern_dim_c
!!  - replaced rbf_vec_dim by rbf_vec_dim_c
!!  - replaced rbf_mat_dim by rbf_mat_dim_c
!!  - replaced rbf_vec_scale by rbf_vec_scale_c
!!  - replaced rbf_vec_pdeg_c by rbf_vec_rbf_vec_pdeg_c_c
!!  Modification by Hui Wan, MPI-M (2007-08-02; 2007-11-30)
!!  - added interpolation coefficients c_aw_e and e_aw_c
!!    and the initialization subroutine aw_int_coeff.
!!  - added subroutine edges2cells_scalar
!!  Modification by Jochen Foerstner, DWD (2008-05-05)
!!  - four new subroutines
!!      rbf_vec_index_vertex
!!      rbf_vec_compute_coeff_vertex
!!      rbf_vec_interpol_car_vertex
!!      prepare_simpson
!!    to reconstruct a Cartesian vector at the vertices using
!!    RBF interpolation and to prepare quadrature via the
!!    Simpson's rule.
!!  Modification by Marco Restelli, MPI (2008-07-17)
!!  - included the subroutines
!!      cells2vertex_scalar, cells2vertex_coeff, ravtom_normgrad2,
!!      ls_normgrad2, ls_normgrad2_ii, edges2points_vector
!!    to compute polynomial fitting with sufficient accuracy as
!!    required in SW-alpha model.
!!  Modification by Jochen Foerstner, DWD (2008-09-12)
!!  - moved SUBROUTINE ravtom_normgrad2 to mo_math_operators
!!    because of conflicting use statements.
!!  Modification by Almut Gassmann, MPI-M (2008-10-09)
!!  - added features for helicity bracket reconstruction
!!  Modification by Guenther Zaengl, DWD (2008-10-23)
!!  - added interpolation routines needed for mesh refinement
!!  Modification by Almut Gassmann, MPI-M (2009-01-29)
!!  - conforming scalar interpolation routines and adjusting coefficients
!!  Modification by Guenther Zaengl, DWD (2009-02-11)
!!  - all routines needed for grid refinement are moved into the new
!!    module mo_grf_interpolation
!!  Modification by Guenther Zaengl, DWD (2009-02-13)
!!  - RBFs are changed to direct reconstruction of velocity components on
!!    the sphere, avoiding the detour over the 3D Cartesian space
!!  Modification by Almut Gassmann, DWD (2009-03-17)
!!  - remove lraviart
!!  Modification by Almut Gassmann, MPI-M (2009-04-23)
!!  - remove all Raviart Thomas stuff, add edge to verts averaging
!!  Modification by Daniel Reinert, DWD (2009-07-20)
!!  - added subroutine grad_lsq_compute_coeff_cell to prepare
!!    (2D) gradient reconstruction at circumcenter via the least squares
!!    method.
!!  Modification by Almut Gassmann, MPI-M (2009-10-05)
!!  - set RBF vec dimensions to predefined values (edges:4,vertices:6,cells:9);
!!    All other switches and belongings are deleted. The reason is that
!!    the Hollingsworth instability requires 4 edges, cell reconstruction
!!    is only needed for output and vertices are only used in the bracket
!!    version, where the dimension at the vertices should be 6
!!  Modification by Daniel Reinert, DWD (2009-12-10)
!!  - replaced grad_lsq_compute_coeff_cell by lsq_compute_coeff_cell
!!    which initializes either a second order or a third order least squares
!!    reconstruction.
!!  Modification by Almut Gassmann, MPI-M (2010-01-12)
!!  - generalize p_int%primal_normal_ec and p_int%edge_cell_length to hexagons
!!
!! @par Copyright
!! 2002-2007 by DWD and MPI-M
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
MODULE mo_intp_coeffs
!-------------------------------------------------------------------------
!
!    ProTeX FORTRAN source: Style 2
!    modified for ICON project, DWD/MPI-M 2006
!
!-------------------------------------------------------------------------
!
!
!
USE mo_kind,                ONLY: wp
USE mo_math_constants,      ONLY: pi2, pi_2,deg2rad
USE mo_physical_constants,  ONLY: re,omega
USE mo_exception,           ONLY: message, finish
USE mo_impl_constants,      ONLY: min_rlcell, min_rledge, min_rlvert, MAX_CHAR_LENGTH,&
  &  beta_plane_coriolis,full_coriolis,min_rledge_int,min_rlcell_int,min_rlvert_int
USE mo_impl_constants_grf,  ONLY: grf_nudge_start_c, grf_nudge_start_e
USE mo_model_domain,        ONLY: t_patch, t_grid_edges, t_grid_vertices, t_grid_cells
USE mo_model_domain_import, ONLY: lplane, lfeedback
USE mo_math_utilities,      ONLY: gc2cc, cc2gc, gnomonic_proj,               &
                                & gvec2cvec, cvec2gvec,                      &
                                & t_cartesian_coordinates,                   &
                                & rotate_latlon, arc_length,                 &
                                & t_geographical_coordinates
USE mo_dynamics_config,     ONLY: divavg_cntrwgt
USE mo_parallel_config,  ONLY: nproma
USE mo_loopindices,         ONLY: get_indices_c, get_indices_e, get_indices_v
Use mo_sync,                ONLY: SYNC_C, SYNC_E, SYNC_V, sync_patch_array, sync_idx, global_max

USE mo_intp_data_strc
USE mo_intp_coeffs_lsq_bln
!USE mo_interpol_nml
USE mo_interpol_config
USE mo_ocean_nml,           ONLY: n_zlev, dzlev_m, no_tracer, t_ref, s_ref,          &
  &                               CORIOLIS_TYPE, basin_center_lat, basin_height_deg

IMPLICIT NONE

PRIVATE

PUBLIC ::  lsq_stencil_create, lsq_compute_coeff_cell, scalar_int_coeff,      &
          & bln_int_coeff_e2c, compute_heli_bra_coeff_idx, init_cellavg_wgt,  &
          & init_geo_factors, complete_patchinfo, init_tplane_e,              &
          & init_geo_factors_oce, init_scalar_product_oce,                    &
          & init_nudgecoeffs, tri_quadrature_pts

CONTAINS


!-------------------------------------------------------------------------
!
!
!>
!! Coefficients and indices for vorticity flux term
!!
!! Computes the coefficents and index informations for the computation of
!! the 2 dimensional projection of the vortex bracket for the hexagonal grid.
!! With that, the omega x v term is computed.
!!
!! @par Revision History
!! Initial release by Almut Gassmann (2008-10-08)
!! Extension to the hexagonal grid by A. Gassmann (2008-11-05)
!! Extension to the complete philosophy on triangles by A. Gassmann (2008-12-18)
!! Modification by Almut Gassmann (2009-02-05)
!! - include new geometry (use of new edge_cell_length, edge_vert_length)
!! Modification by Almut Gassmann (2009-12-20)
!! - changed structure of the code to be more readable and efficient
!! Modification by Almut Gassmann (2010-02-03)
!! - Optimizing the code and avoiding possible errors in vector projections.
!!   This improved the result remarkbly.
!! Modification by Almut Gassmann (2010-03-16)
!! - Implementation of the method of Thuburn/Ringler/Skamarock/Klemp (JCP 228),
!!   which is different from my approach. This is done for comparison.
!!
SUBROUTINE compute_heli_bra_coeff_idx (ptr_patch, ptr_int)

  IMPLICIT NONE

  TYPE(t_patch), TARGET, INTENT(inout) :: ptr_patch    ! patch

  ! interpolation state
  TYPE(t_int_state), TARGET, INTENT(inout):: ptr_int

  TYPE(t_cartesian_coordinates) :: z_cart_no, z_cart_ea
  INTEGER    :: jb, je, jc, nlen, nblks_e, npromz_e, nblks_c, npromz_c, &
  &             n_edges, imyself, ineigh, incr, ile1_v(2), ibe1_v(2),&
  &             je1, ile1, ibe1, je2, ile2, ibe2, je3, ile3, ibe3, je4, je5, &
  &             ile5, ibe5, jc1, ilc1, ibc1, jc3, ilc3, ibc3, jc4, jv2, jv3, &
  &             ilv3, ibv3, jv4,ilc,ibc, npromz_v, nblks_v
  REAL(wp)   :: z_tan_proj, z_nor_proj, z_metric, z_metric_new, z_norm, &
                z_metric_a, z_metric_b, z_metric_c, z_metric_d
  REAL(wp),ALLOCATABLE :: z_frac_area(:,:,:), &
  &                       z_quad_north(:,:,:), z_quad_east(:,:,:)
  TYPE(t_grid_edges),    POINTER :: p_ed
  TYPE(t_grid_cells),    POINTER :: p_ce
  TYPE(t_grid_vertices), POINTER :: p_ve
  TYPE(t_patch),         POINTER :: p_pa
  TYPE(t_int_state),     POINTER :: p_in
  LOGICAL :: l_cycle, l_found
  INTEGER :: ile_c(6), ibe_c(6), n_verts, jv, ilv, ibv, jm, ilv1,ibv1, &
  &          ineigh_c, ineigh_v, ilv_c(6), ibv_c(6), incrc,            &
  &          jvindex, jeindex, ile, ibe
  REAL(wp):: z_coeff(0:6), zorient_ce, zorient_ve

!-----------------------------------------------------------------------

  !i_cori_method = 1 ! Almut's solution for reconstruction, for PV as in TRSK
  !i_cori_method = 2 ! Thuburn, Skamarock, Klemp, Ringler method
  !i_cori_method = 3,4 ! Almut's solution for reconstruction and for PV

  !USA/GB solution
  !---------------
  IF (i_cori_method == 2) THEN

    p_pa => ptr_patch
    p_ed => ptr_patch%edges
    p_ce => ptr_patch%cells
    p_ve => ptr_patch%verts
    p_in => ptr_int

    ALLOCATE(z_frac_area(6,nproma,p_pa%nblks_c))
    z_frac_area = 0._wp

    ! Fractional areas R_{i,v}
    nblks_c  = p_pa%nblks_int_c
    npromz_c = p_pa%npromz_int_c
    DO jb = 1, nblks_c
      IF (jb /= nblks_c) THEN
        nlen = nproma
      ELSE
        nlen = npromz_c
      ENDIF
      DO jc = 1, nlen

        IF(.NOT.p_ce%owner_mask(jc,jb)) CYCLE

        n_edges = p_ce%num_edges(jc,jb)
        DO je1 = 1, n_edges
          ! edges of cells
          ile_c(je1) = p_ce%edge_idx(jc,jb,je1)
          ibe_c(je1) = p_ce%edge_blk(jc,jb,je1)
        ENDDO
        n_verts = n_edges
        DO jv = 1, n_verts
          ! Corner
          z_frac_area(jv,jc,jb) = 0.0_wp
          ilv = p_ce%vertex_idx(jc,jb,jv)
          ibv = p_ce%vertex_blk(jc,jb,jv)
          ! edges of corner
          DO je1 = 1, 3
            ile1 = p_ve%edge_idx(ilv,ibv,je1)
            ibe1 = p_ve%edge_blk(ilv,ibv,je1)
            DO je2 = 1, 6
              IF(ile1 == ile_c(je2) .AND. ibe1 == ibe_c(je2)) THEN
                IF(jc == p_ed%cell_idx(ile1,ibe1,1) .AND. &
                   jb == p_ed%cell_blk(ile1,ibe1,1)) THEN
                  ineigh_c = 1
                ELSE
                  ineigh_c = 2
                ENDIF
                IF(ilv == p_ed%vertex_idx(ile1,ibe1,1) .AND. &
                   ibv == p_ed%vertex_blk(ile1,ibe1,1)) THEN
                  ineigh_v = 1
                ELSE
                  ineigh_v = 2
                ENDIF
                z_frac_area(jv,jc,jb) = z_frac_area(jv,jc,jb) &
                & + p_ed%edge_cell_length(ile1,ibe1,ineigh_c) &
                & * p_ed%edge_vert_length(ile1,ibe1,ineigh_v) * 0.5_wp &
                & / p_ce%area(jc,jb)
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    DO jv = 1, UBOUND(z_frac_area,1)
      CALL sync_patch_array(SYNC_C,ptr_patch,z_frac_area(jv,:,:))
    ENDDO

    ! heli coeffs
    nblks_e  = p_pa%nblks_int_e
    npromz_e = p_pa%npromz_int_e
    DO jb = 1, nblks_e
      IF (jb /= nblks_e) THEN
        nlen = nproma
      ELSE
        nlen = npromz_e
      ENDIF
      DO je = 1, nlen

        IF(.NOT.p_ed%owner_mask(je,jb)) CYCLE

        incr = 0
        DO jc = 1, 2

          ilc = p_ed%cell_idx(je,jb,jc)
          ibc = p_ed%cell_blk(je,jb,jc)

          n_verts = p_ce%num_edges(ilc,ibc)
          n_edges = n_verts
          ! index of corners
          DO jv = 1, n_verts
            ilv_c(jv) = p_ce%vertex_idx(ilc,ibc,jv)
            ibv_c(jv) = p_ce%vertex_blk(ilc,ibc,jv)
          ENDDO
          ! index of edges
          DO je1 = 1, n_edges
            ile_c(je1) = p_ce%edge_idx(ilc,ibc,je1)
            ibe_c(je1) = p_ce%edge_blk(ilc,ibc,je1)
            ! orientation of the base edge
            IF((ile_c(je1) == je) .and. (ibe_c(je1)==jb) ) THEN
              zorient_ce= p_ce%edge_orientation(ilc,ibc,je1)
            ENDIF
          ENDDO

          ile1 = je
          ibe1 = jb

          incrc = 0
          z_coeff(incrc) = -0.5_wp

          ! do the circle around
          DO jm = 1, n_edges-1

            IF (jm == 1) THEN
              ! select an arbitrary corner in the first step
              ilv1=p_ed%vertex_idx(ile1,ibe1,1)
              ibv1=p_ed%vertex_blk(ile1,ibe1,1)
            ELSE
              ! select the next corner
              DO jv = 1,2
                IF(.NOT.(p_ed%vertex_idx(ile1,ibe1,jv) ==ilv1 .AND. &
                         p_ed%vertex_blk(ile1,ibe1,jv) ==ibv1) ) THEN
                  jvindex = jv
                ENDIF
              ENDDO
              ilv1=p_ed%vertex_idx(ile1,ibe1,jvindex)
              ibv1=p_ed%vertex_blk(ile1,ibe1,jvindex)
            ENDIF
            ! select the corner index with repect to cell
            DO jv = 1, n_verts
              IF (ilv1 == ilv_c(jv) .AND. ibv1 == ibv_c(jv)) THEN
                jvindex = jv
                EXIT
              ENDIF
            ENDDO
            ! select the next edge index with respect to cell and
            ! determine the edge orientation with respect to vertex
            l_found=.FALSE.
            DO je1 = 1, 3
              ile = p_ve%edge_idx(ilv1,ibv1,je1)
              ibe = p_ve%edge_blk(ilv1,ibv1,je1)
              IF(ile == ile1 .and. ibe == ibe1) CYCLE ! that was the previous
              DO je2 = 1, n_edges
                IF(ile == ile_c(je2) .AND. ibe==ibe_c(je2))THEN
                  jeindex = je2
                  zorient_ve=p_ve%edge_orientation(ilv1,ibv1,je1)
                  l_found = .TRUE.
                  EXIT
                ENDIF
              ENDDO
              IF (l_found) EXIT
            ENDDO

            incr = incr+1
            incrc = incrc+1

            z_coeff(incrc) = z_coeff(incrc-1)+z_frac_area(jvindex,ilc,ibc)

            p_in%heli_vn_idx(incr,je,jb) = ile_c(jeindex)
            p_in%heli_vn_blk(incr,je,jb) = ibe_c(jeindex)
            p_in%heli_coeff (incr,je,jb) = - 0.5_wp*z_coeff(incrc)*zorient_ce&
            &     *p_ed%primal_edge_length(ile_c(jeindex),ibe_c(jeindex)) &
            &     /p_ed%dual_edge_length(je,jb)*zorient_ve

            ile1 = ile_c(jeindex)
            ibe1 = ibe_c(jeindex)

          ENDDO
          ! If there is a pentagon
          IF (jc==2.and.incr==9) THEN
            incr = incr+1
            p_in%heli_vn_idx(incr,je,jb) = je
            p_in%heli_vn_blk(incr,je,jb) = jb
            p_in%heli_coeff (incr,je,jb) = 0.0_wp
          ENDIF
        ENDDO
      ENDDO
    ENDDO

    DO incr = 1, UBOUND(p_in%heli_coeff,1)
      CALL sync_patch_array(SYNC_E,ptr_patch,p_in%heli_coeff(incr,:,:))
      CALL sync_idx(SYNC_E,SYNC_E,ptr_patch,p_in%heli_vn_idx(incr,:,:),p_in%heli_vn_blk(incr,:,:))
    ENDDO

  ENDIF

  p_pa => ptr_patch
  p_ed => ptr_patch%edges
  p_ce => ptr_patch%cells
  p_ve => ptr_patch%verts
  p_in => ptr_int

  ALLOCATE(z_quad_east (6,nproma,p_pa%nblks_e))
  ALLOCATE(z_quad_north(6,nproma,p_pa%nblks_e))

  ! Vector reconstruction on hexagons
  !==================================
  nblks_c  = p_pa%nblks_int_c
  npromz_c = p_pa%npromz_int_c
  DO jb = 1, nblks_c
    IF (jb /= nblks_c) THEN
      nlen = nproma
    ELSE
      nlen = npromz_c
    ENDIF
    DO jc = 1, nlen

      IF(.NOT.p_ce%owner_mask(jc,jb)) CYCLE

      IF (.NOT. lplane) THEN
        CALL gvec2cvec(0.0_wp,1.0_wp,p_ce%center(jc,jb)%lon,p_ce%center(jc,jb)%lat, &
        &              z_cart_no%x(1),z_cart_no%x(2),z_cart_no%x(3))
        z_norm = SQRT( DOT_PRODUCT(z_cart_no%x(1:3),z_cart_no%x(1:3)))
        z_cart_no%x(1:3) = z_cart_no%x(1:3) /z_norm
        CALL gvec2cvec(1.0_wp,0.0_wp,p_ce%center(jc,jb)%lon,p_ce%center(jc,jb)%lat, &
        &              z_cart_ea%x(1),z_cart_ea%x(2),z_cart_ea%x(3))
        z_norm = SQRT( DOT_PRODUCT(z_cart_ea%x(1:3),z_cart_ea%x(1:3)))
        z_cart_ea%x(1:3) = z_cart_ea%x(1:3) /z_norm
      ENDIF
      n_edges = p_ce%num_edges(jc,jb)
      DO je1 = 1,n_edges
        ile1 = p_ce%edge_idx(jc,jb,je1)
        ibe1 = p_ce%edge_blk(jc,jb,je1)
        ! Projection
        !-----------
        ! Inner product of T(jc,jb)*N(ile1,ibe1)
        IF (.NOT. lplane) THEN
          z_tan_proj = DOT_PRODUCT(p_ed%primal_cart_normal(ile1,ibe1)%x,z_cart_no%x)
          z_nor_proj = DOT_PRODUCT(p_ed%primal_cart_normal(ile1,ibe1)%x,z_cart_ea%x)
        ELSE
          z_tan_proj = p_ed%primal_normal(ile1,ibe1)%v2
          z_nor_proj = p_ed%primal_normal(ile1,ibe1)%v1
        ENDIF
        z_metric = 0.5_wp*p_ed%dual_edge_length(ile1,ibe1)&
        &                *p_ed%primal_edge_length(ile1,ibe1)     &
        &                /p_ce%area(jc,jb)
        ! north projection
        p_in%hex_north(jc,je1,jb) = z_tan_proj*z_metric
        ! east projection
        p_in%hex_east(jc,je1,jb)  = z_nor_proj*z_metric
      ENDDO
    ENDDO
  ENDDO

  CALL sync_patch_array(SYNC_C,ptr_patch,p_in%hex_north)
  CALL sync_patch_array(SYNC_C,ptr_patch,p_in%hex_east)

  ! Vector reconstruction on triangles
  !===================================
  nblks_v  = p_pa%nblks_int_v
  npromz_v = p_pa%npromz_int_v
  DO jb = 1, nblks_v
    IF (jb /= nblks_v) THEN
      nlen = nproma
    ELSE
      nlen = npromz_v
    ENDIF
    DO jv = 1, nlen

      IF(.NOT.p_ve%owner_mask(jv,jb)) CYCLE

      IF (.NOT. lplane) THEN
        CALL gvec2cvec(0.0_wp,1.0_wp,p_ve%vertex(jv,jb)%lon,p_ve%vertex(jv,jb)%lat, &
        &              z_cart_no%x(1),z_cart_no%x(2),z_cart_no%x(3))
        z_norm = SQRT( DOT_PRODUCT(z_cart_no%x(1:3),z_cart_no%x(1:3)))
        z_cart_no%x(1:3) = z_cart_no%x(1:3) /z_norm
        CALL gvec2cvec(1.0_wp,0.0_wp,p_ve%vertex(jv,jb)%lon,p_ve%vertex(jv,jb)%lat, &
        &              z_cart_ea%x(1),z_cart_ea%x(2),z_cart_ea%x(3))
        z_norm = SQRT( DOT_PRODUCT(z_cart_ea%x(1:3),z_cart_ea%x(1:3)))
        z_cart_ea%x(1:3) = z_cart_ea%x(1:3) /z_norm
      ENDIF
      DO je1 = 1, 3
        ile1 = p_ve%edge_idx(jv,jb,je1)
        ibe1 = p_ve%edge_blk(jv,jb,je1)
        IF((p_ed%vertex_idx(ile1,ibe1,1)==jv) .AND. &
        &  (p_ed%vertex_blk(ile1,ibe1,1)==jb)) THEN
          ineigh = 1
        ELSE
          ineigh = 2
        ENDIF
        ! Projection
        !-----------
        IF (.NOT. lplane) THEN
          z_tan_proj = DOT_PRODUCT(p_ed%primal_cart_normal(ile1,ibe1)%x,z_cart_no%x)
          z_nor_proj = DOT_PRODUCT(p_ed%primal_cart_normal(ile1,ibe1)%x,z_cart_ea%x)
        ELSE
          z_tan_proj = p_ed%primal_normal(ile1,ibe1)%v2
          z_nor_proj = p_ed%primal_normal(ile1,ibe1)%v1
        ENDIF
        z_metric = p_ed%dual_edge_length(ile1,ibe1)       &
        &         *p_ed%edge_vert_length(ile1,ibe1,ineigh)&
        &         /p_ve%dual_area(jv,jb)
        ! north projection
        p_in%tria_north(je1,jv,jb) = z_tan_proj*z_metric
        ! east projection
        p_in%tria_east(je1,jv,jb)  = z_nor_proj*z_metric
      ENDDO
    ENDDO
  ENDDO

  DO je1 = 1, 3
    CALL sync_patch_array(SYNC_V,ptr_patch,p_in%tria_north(je1,:,:))
    CALL sync_patch_array(SYNC_V,ptr_patch,p_in%tria_east(je1,:,:))
  ENDDO

  IF (i_cori_method /= 2 ) THEN

    ! Vector reconstruction on rhombi
    !================================
    nblks_e  = p_pa%nblks_int_e
    npromz_e = p_pa%npromz_int_e
    DO jb = 1, nblks_e
      IF (jb /= nblks_e) THEN
        nlen = nproma
      ELSE
        nlen = npromz_e
      ENDIF
      DO je = 1, nlen

        IF(.NOT.p_ed%owner_mask(je,jb)) CYCLE

        IF (.NOT. lplane) THEN
          CALL gvec2cvec(0.0_wp,1.0_wp,p_ed%center(je,jb)%lon,p_ed%center(je,jb)%lat, &
          &              z_cart_no%x(1),z_cart_no%x(2),z_cart_no%x(3))
          z_norm = SQRT( DOT_PRODUCT(z_cart_no%x(1:3),z_cart_no%x(1:3)))
          z_cart_no%x(1:3) = z_cart_no%x(1:3) /z_norm
          CALL gvec2cvec(1.0_wp,0.0_wp,p_ed%center(je,jb)%lon,p_ed%center(je,jb)%lat, &
          &              z_cart_ea%x(1),z_cart_ea%x(2),z_cart_ea%x(3))
          z_norm = SQRT( DOT_PRODUCT(z_cart_ea%x(1:3),z_cart_ea%x(1:3)))
          z_cart_ea%x(1:3) = z_cart_ea%x(1:3) /z_norm
        ENDIF
        DO je1 = 1,6
          IF (je1 <=4 )THEN
            ile1 = p_ed%quad_idx(je,jb,je1)
            ibe1 = p_ed%quad_blk(je,jb,je1)
            IF(((p_ed%vertex_idx(ile1,ibe1,1)==p_ed%vertex_idx(je,jb,1)) .AND. &
            &   (p_ed%vertex_blk(ile1,ibe1,1)==p_ed%vertex_blk(je,jb,1))) .OR. &
            &  ((p_ed%vertex_idx(ile1,ibe1,1)==p_ed%vertex_idx(je,jb,2)) .AND. &
            &   (p_ed%vertex_blk(ile1,ibe1,1)==p_ed%vertex_blk(je,jb,2))))THEN
              ineigh = 1
            ELSE
              ineigh = 2
            ENDIF
          ELSE
            ile1 = je
            ibe1 = jb
            ineigh = 7-je1
          ENDIF
          ! Projection
          !-----------
          ! Inner product of T(je,jb)*N(je1,je1)
          IF (.NOT. lplane) THEN
            z_tan_proj = DOT_PRODUCT(p_ed%primal_cart_normal(ile1,ibe1)%x,z_cart_no%x)
            z_nor_proj = DOT_PRODUCT(p_ed%primal_cart_normal(ile1,ibe1)%x,z_cart_ea%x)
          ELSE
            z_tan_proj = p_ed%primal_normal(ile1,ibe1)%v2
            z_nor_proj = p_ed%primal_normal(ile1,ibe1)%v1
          ENDIF
          z_metric = p_ed%edge_vert_length(ile1,ibe1,ineigh)  &
          &         *p_ed%dual_edge_length(ile1,ibe1)         &
          &         /p_ed%quad_area(je,jb)
          ! north projection
          z_quad_north(je1,je,jb) = z_tan_proj*z_metric
          ! east projection
          z_quad_east(je1,je,jb)  = z_nor_proj*z_metric
        ENDDO
        IF(i_cori_method >= 3) THEN
          p_in%quad_north(1:4,je,jb) = z_quad_north(1:4,je,jb)
          p_in%quad_east (1:4,je,jb) = z_quad_east (1:4,je,jb)
          p_in%quad_north(5,je,jb) = z_quad_north(5,je,jb)+ z_quad_north(6,je,jb)
          p_in%quad_east (5,je,jb) = z_quad_east (5,je,jb)+ z_quad_east (6,je,jb)
        ENDIF
      ENDDO
    ENDDO

    IF(i_cori_method >= 3) THEN
      DO je1 = 1, 5
        CALL sync_patch_array(SYNC_E,ptr_patch,p_in%quad_north(je1,:,:))
        CALL sync_patch_array(SYNC_E,ptr_patch,p_in%quad_east(je1,:,:))
      ENDDO
    ENDIF
    DO je1 = 1, 6
      CALL sync_patch_array(SYNC_E,ptr_patch,z_quad_north(je1,:,:))
      CALL sync_patch_array(SYNC_E,ptr_patch,z_quad_east (je1,:,:))
    ENDDO

    ! Computation of the coefficients for the vorticity flux term
    !============================================================
    nblks_e  = p_pa%nblks_int_e
    npromz_e = p_pa%npromz_int_e

    DO jb = 1, nblks_e
      IF (jb /= nblks_e) THEN
        nlen = nproma
      ELSE
        nlen = npromz_e
      ENDIF
      DO je = 1, nlen

        IF(.NOT.p_ed%owner_mask(je,jb)) CYCLE

        incr    = 0

        !============================================================
        ! Look for those points which do not appear also in a rhombus
        !============================================================
        DO jc1 = 1, 2

          ilc1 = p_ed%cell_idx(je,jb,jc1)
          ibc1 = p_ed%cell_blk(je,jb,jc1)
          n_edges = p_ce%num_edges(ilc1,ibc1)

          DO je2 = 1, n_edges ! loop over the edges
            ile2 = p_ce%edge_idx(ilc1,ibc1,je2)
            ibe2 = p_ce%edge_blk(ilc1,ibc1,je2)
            IF (ile2 == je .AND. ibe2==jb) THEN
              imyself = je2
            ENDIF
          ENDDO

          !=========
          ! method 3
          !=========
          IF (i_cori_method >= 3) THEN
            p_in%heli_coeff(2*jc1-1,je,jb) =  p_in%hex_east (ilc1,imyself,ibc1) &
            &                *0.5_wp*p_ce%area(ilc1,ibc1)/p_ed%area_edge(je,jb)
            p_in%heli_coeff(2*jc1  ,je,jb) = -p_in%hex_north(ilc1,imyself,ibc1) &
            &                *0.5_wp*p_ce%area(ilc1,ibc1)/p_ed%area_edge(je,jb)
          ENDIF

          !=========
          ! method 1
          !=========
          IF (i_cori_method == 1) THEN

            DO je2 = 1, n_edges ! loop over the edges

              ile2 = p_ce%edge_idx(ilc1,ibc1,je2)
              ibe2 = p_ce%edge_blk(ilc1,ibc1,je2)

              IF (ile2 == je .AND. ibe2==jb) CYCLE ! same edges
              l_cycle=.FALSE.
              DO je3 = 1, 4
                ile3 = p_ed%quad_idx(je,jb,je3)
                ibe3 = p_ed%quad_blk(je,jb,je3)
                DO je4 = 1, 4
                  IF (p_ed%quad_idx(ile3,ibe3,je4)==ile2.AND.&
                      p_ed%quad_blk(ile3,ibe3,je4)==ibe2)THEN
                    l_cycle = .TRUE.
                  ENDIF
                ENDDO
              ENDDO
              IF (l_cycle) CYCLE

              incr = incr + 1  ! count the edges

              ! indices for normal velocity components
              p_in%heli_vn_idx(incr,je,jb) = ile2
              p_in%heli_vn_blk(incr,je,jb) = ibe2

              z_metric_new = 0.5_wp*p_ce%area(ilc1,ibc1)/p_ed%area_edge(je,jb)   &
              & *(p_in%hex_east (ilc1,imyself,ibc1)*p_in%hex_north(ilc1,je2,ibc1) &
              &  -p_in%hex_north(ilc1,imyself,ibc1)*p_in%hex_east (ilc1,je2,ibc1))

                ! here, the factor 1/2 for averaging the PVs is already incorporated
                p_in%heli_coeff(incr,je,jb)=z_metric_new*0.5_wp

            ENDDO ! je2
            IF (incr == 1 .AND. jc1 ==2) THEN
              !This would appear if one of the cells
              !is a pentagon. Do a dummy contribution here.
              incr = incr + 1
              ! indices for normal velocity components
              p_in%heli_vn_idx(incr,je,jb) = je
              p_in%heli_vn_blk(incr,je,jb) = jb
              ! coefficients for helicity term
              p_in%heli_coeff(incr,je,jb)  = 0.0_wp
            ENDIF

          ENDIF ! cori_method==1

        ENDDO  ! jc1

        !=========
        ! method 3
        !=========
        IF (i_cori_method>=3) THEN
          p_in%heli_coeff(5,je,jb) =  p_in%quad_east (5,je,jb)/6.0_wp &
          &              *p_ed%quad_area(je,jb)/p_ed%area_edge(je,jb)
          p_in%heli_coeff(6,je,jb) = -p_in%quad_north(5,je,jb)/6.0_wp &
          &              *p_ed%quad_area(je,jb)/p_ed%area_edge(je,jb)
        ENDIF


        ! Contribution from rhombus neighbors
        !====================================
        DO je1 = 1, 4  ! loop over neighboring quads

          ! Indices
          ile1 = p_ed%quad_idx(je,jb,je1)
          ibe1 = p_ed%quad_blk(je,jb,je1)

          !=========
          ! method 3
          !=========
          IF (i_cori_method >= 3) THEN
            DO je2 = 1, 4 ! loop over the edges
              ile2 = p_ed%quad_idx(ile1,ibe1,je2)
              ibe2 = p_ed%quad_blk(ile1,ibe1,je2)
              IF (ile2 == je .AND. ibe2==jb) THEN
                imyself = je2
              ENDIF
            ENDDO
            p_in%heli_coeff(2*je1+5,je,jb) =  p_in%quad_east (imyself,ile1,ibe1)/6.0_wp &
            &                          *p_ed%quad_area(ile1,ibe1)/p_ed%area_edge(je,jb)
            p_in%heli_coeff(2*je1+6,je,jb) = -p_in%quad_north(imyself,ile1,ibe1)/6.0_wp &
            &                          *p_ed%quad_area(ile1,ibe1)/p_ed%area_edge(je,jb)

          ENDIF

          !=========
          ! method 1
          !=========
          IF (i_cori_method== 1) THEN

          ! vertices of rhombus e1
          DO jv2 = 1, 2
            ile1_v(jv2) = p_ed%vertex_idx(ile1,ibe1,jv2)
            ibe1_v(jv2) = p_ed%vertex_blk(ile1,ibe1,jv2)
          ENDDO

          ! Now counting all the outer edges of the rhombus

          DO je2 = 1, 4

            ile2 = p_ed%quad_idx(ile1,ibe1,je2)
            ibe2 = p_ed%quad_blk(ile1,ibe1,je2)
            IF (je == ile2 .AND. jb == ibe2) CYCLE ! same edge

            incr = incr + 1

            !======================================================
            ! Rhombus e1, this is always applicable (trice per e1).
            !======================================================

            DO je3 = 1, 4
              ile3 = p_ed%quad_idx(ile1,ibe1,je3)
              ibe3 = p_ed%quad_blk(ile1,ibe1,je3)
              IF (ile3 == je .AND. ibe3 == jb) THEN
                imyself = je3
              ENDIF
            ENDDO
            z_metric_a = 1.0_wp/6.0_wp &
            & *p_ed%quad_area(ile1,ibe1)/p_ed%area_edge(je,jb) &
            & *(z_quad_east (imyself,ile1,ibe1)*z_quad_north(je2,ile1,ibe1) &
            &  -z_quad_north(imyself,ile1,ibe1)*z_quad_east (je2,ile1,ibe1))
            z_metric_new = 0.5_wp*z_metric_a

            !=========================================
            ! Hexagonal part, this occurs twice per e1
            !=========================================
            DO jc3 = 1, 2
              ilc3 = p_ed%cell_idx(ile2,ibe2,jc3)
              ibc3 = p_ed%cell_blk(ile2,ibe2,jc3)
              DO jc4 = 1, 2
                IF( ilc3 == p_ed%cell_idx(je,jb,jc4) .AND. &
                &   ibc3 == p_ed%cell_blk(je,jb,jc4) ) THEN

                  DO je5 = 1, p_ce%num_edges(ilc3,ibc3)
                    ile5 = p_ce%edge_idx(ilc3,ibc3,je5)
                    ibe5 = p_ce%edge_blk(ilc3,ibc3,je5)
                    IF (ile5 == je   .AND. ibe5 == jb  ) THEN
                      imyself = je5
                    ELSE IF (ile5 == ile2 .AND. ibe5 == ibe2) THEN
                      ineigh  = je5
                    ENDIF
                  ENDDO
                  z_metric_b = 0.5_wp*p_ce%area(ilc3,ibc3)/p_ed%area_edge(je,jb)          &
                  & *(p_in%hex_east (ilc3,imyself,ibc3)*p_in%hex_north(ilc3,ineigh,ibc3) &
                  &  -p_in%hex_north(ilc3,imyself,ibc3)*p_in%hex_east (ilc3,ineigh,ibc3))
                  z_metric_new = z_metric_new + 0.5_wp*z_metric_b
                ENDIF
              ENDDO
            ENDDO

            ! The short distance e2 point goes into rhombi e2 and je
            DO jv3 = 1, 2
              ilv3 = p_ed%vertex_idx(je,jb,jv3)
              ibv3 = p_ed%vertex_blk(je,jb,jv3)
              DO jv4 = 1, 2
                IF(ilv3==p_ed%vertex_idx(ile2,ibe2,jv4).AND.&
                   ibv3==p_ed%vertex_blk(ile2,ibe2,jv4)) THEN

                  !========================================
                  ! Rhombus e2 (e2 is the center edge here)
                  !========================================

                  DO je5 = 1, 4
                    ile5 = p_ed%quad_idx(ile2,ibe2,je5)
                    ibe5 = p_ed%quad_blk(ile2,ibe2,je5)
                    IF (ile5 == je .AND. ibe5 == jb) THEN
                      imyself = je5
                    ENDIF
                  ENDDO
                  z_metric_c = 1.0_wp/6.0_wp &
                  & *p_ed%quad_area(ile2,ibe2)/p_ed%area_edge(je,jb) &
                  & *(z_quad_east(imyself,ile2,ibe2) &
                  &  *(z_quad_north(5,ile2,ibe2)+z_quad_north(6,ile2,ibe2))&
                  &  -z_quad_north(imyself,ile2,ibe2) &
                  &  *(z_quad_east (5,ile2,ibe2)+z_quad_east (6,ile2,ibe2)))
                  z_metric_new = z_metric_new+0.5_wp*z_metric_c

                  !======================
                  ! Rhombus je
                  !======================

                  DO je5 = 1, 4
                    ile5 = p_ed%quad_idx(je,jb,je5)
                    ibe5 = p_ed%quad_blk(je,jb,je5)
                    IF (ile5 == ile2 .AND. ibe5 == ibe2) THEN
                      ineigh = je5
                    ENDIF
                  ENDDO
                  z_metric_d = 1.0_wp/6.0_wp &
                  & *p_ed%quad_area(je,jb)/p_ed%area_edge(je,jb) &
                  & *((z_quad_east (5,je,jb)+z_quad_east (6,je,jb)) &
                  &   *z_quad_north(ineigh,je,jb)  &
                  &  -(z_quad_north(5,je,jb)+z_quad_north(6,je,jb)) &
                  &   *z_quad_east (ineigh,je,jb))
                  z_metric_new = z_metric_new+0.5_wp*z_metric_d

                ENDIF
              ENDDO
            ENDDO

            ! indices for normal velocity components
            p_in%heli_vn_idx(incr,je,jb) = ile2
            p_in%heli_vn_blk(incr,je,jb) = ibe2
            ! coefficient
            p_in%heli_coeff(incr,je,jb) = z_metric_new


          ENDDO ! je1
          ENDIF ! cori_method == 1
        ENDDO ! jq
      ENDDO ! je
    ENDDO ! jb

    DO incr = 1, UBOUND(p_in%heli_coeff,1)
      CALL sync_patch_array(SYNC_E,ptr_patch,p_in%heli_coeff(incr,:,:))
      IF(i_cori_method < 3) &
     & CALL sync_idx(SYNC_E,SYNC_E,ptr_patch,p_in%heli_vn_idx(incr,:,:),p_in%heli_vn_blk(incr,:,:))
    ENDDO

  ENDIF

  DEALLOCATE(z_quad_east)
  DEALLOCATE(z_quad_north)

END SUBROUTINE compute_heli_bra_coeff_idx
!-------------------------------------------------------------------------

!-------------------------------------------------------------------------
!
!
!>
!! Computes the weighting coefficients for cell averaging with.
!!
!! Computes the weighting coefficients for cell averaging with
!! variable interpolation factors. Results are stored in ptr_patch%cells%avg_wgt
!!
!! @par Revision History
!!  developed by Guenther Zaengl, 2008-12-05
!! @par
!!  modification by Guenther Zaengl, 2009-09-02
!!  revised weights to achieve mass conservation
!!
SUBROUTINE init_cellavg_wgt( ptr_patch, ptr_int )
!
!
!  patch on which computation is performed
!
TYPE(t_patch), TARGET, INTENT(inout) :: ptr_patch

! Interpolation state
TYPE(t_int_state), INTENT(inout):: ptr_int
!

INTEGER :: jc, je, jb, iter, niter
INTEGER :: rl_start, rl_end
INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx, i_nchdom

INTEGER :: ilc1, ibc1, ilc2, ibc2, ilc3, ibc3, inb1, inb2, inb3, ie4, ie5
INTEGER :: ile1, ibe1, ile2, ibe2, ile3, ibe3, ile4, ibe4
INTEGER, DIMENSION(nproma) :: iie1, iie2, iie3, iie4
REAL(wp), DIMENSION (nproma,3) :: z_nx1, z_nx2, z_nx3, z_nx4, z_nx5
REAL(wp) :: checksum(nproma,ptr_patch%nblks_e)

REAL(wp) :: xtemp,ytemp,wgt(3),xloc,yloc,x(3),y(3), &
            pollat,pollon,relax_coeff,wgt_loc,maxwgt_loc,minwgt_loc

REAL(wp), DIMENSION(nproma,ptr_patch%nblks_c)  :: wgt_loc_sum, resid
INTEGER, DIMENSION(nproma,ptr_patch%nblks_c,3) :: inv_neighbor_id
REAL(wp), DIMENSION(nproma,ptr_patch%nblks_c,3) :: z_inv_neighbor_id

#ifdef DEBUG_COEFF
REAL(wp) :: sum1
#endif
!-----------------------------------------------------------------------

! Number of iterations for computation of bilinear weights
niter = 1000

! Relaxation coefficient for adaptation of local weight (empirically determined)
relax_coeff = 0.46_wp

! Initial weighting factor of the local grid point
wgt_loc = divavg_cntrwgt

! Maximum/minimum  weighting factors of the local grid point
maxwgt_loc = wgt_loc + 0.003_wp
minwgt_loc = wgt_loc - 0.003_wp

! Initialization of the residuum  field
resid(:,:) = 0._wp

! values for the blocking
rl_start = 2
rl_end = min_rlcell

i_nchdom   = MAX(1,ptr_patch%n_childdom)
i_startblk = ptr_patch%cells%start_blk(rl_start,1)
i_endblk   = ptr_patch%cells%end_blk(rl_end,i_nchdom)
!
! Compute inverse neighbor ID's
! The inverse neigbor ID of a neighbor cell (ilc1,ibc1) is the neighbor ID
! the local cell (jc,jb) has from the point of view of the neighbor cell

inv_neighbor_id = 0

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx,ilc1,ibc1,ilc2,ibc2,ilc3,ibc3)
DO jb = i_startblk, i_endblk

  CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
                     i_startidx, i_endidx, rl_start, rl_end)

  DO jc = i_startidx, i_endidx

    IF(.NOT.ptr_patch%cells%owner_mask(jc,jb)) CYCLE

    ! line and block indices of the neighbouring cells

    ilc1 = ptr_patch%cells%neighbor_idx(jc,jb,1)
    ibc1 = ptr_patch%cells%neighbor_blk(jc,jb,1)
    ilc2 = ptr_patch%cells%neighbor_idx(jc,jb,2)
    ibc2 = ptr_patch%cells%neighbor_blk(jc,jb,2)
    ilc3 = ptr_patch%cells%neighbor_idx(jc,jb,3)
    ibc3 = ptr_patch%cells%neighbor_blk(jc,jb,3)

    IF ( (ilc1>0) .AND. (ibc1>0) ) THEN
      IF ((ptr_patch%cells%neighbor_idx(ilc1,ibc1,1) == jc) .AND. &
          (ptr_patch%cells%neighbor_blk(ilc1,ibc1,1) == jb)) THEN
        inv_neighbor_id(jc,jb,1) = 1
      ELSE IF ((ptr_patch%cells%neighbor_idx(ilc1,ibc1,2) == jc) .AND. &
               (ptr_patch%cells%neighbor_blk(ilc1,ibc1,2) == jb)) THEN
        inv_neighbor_id(jc,jb,1) = 2
      ELSE IF ((ptr_patch%cells%neighbor_idx(ilc1,ibc1,3) == jc) .AND. &
               (ptr_patch%cells%neighbor_blk(ilc1,ibc1,3) == jb)) THEN
        inv_neighbor_id(jc,jb,1) = 3
      ENDIF
    ENDIF
    IF ( (ilc2>0) .AND. (ibc2>0) ) THEN
      IF ((ptr_patch%cells%neighbor_idx(ilc2,ibc2,1) == jc) .AND. &
          (ptr_patch%cells%neighbor_blk(ilc2,ibc2,1) == jb)) THEN
        inv_neighbor_id(jc,jb,2)  = 1
      ELSE IF ((ptr_patch%cells%neighbor_idx(ilc2,ibc2,2) == jc) .AND. &
               (ptr_patch%cells%neighbor_blk(ilc2,ibc2,2) == jb)) THEN
        inv_neighbor_id(jc,jb,2)  = 2
      ELSE IF ((ptr_patch%cells%neighbor_idx(ilc2,ibc2,3) == jc) .AND. &
               (ptr_patch%cells%neighbor_blk(ilc2,ibc2,3) == jb)) THEN
        inv_neighbor_id(jc,jb,2)  = 3
      ENDIF
    ENDIF
    IF ( (ilc3>0) .AND. (ibc3>0) ) THEN
      IF ((ptr_patch%cells%neighbor_idx(ilc3,ibc3,1) == jc) .AND. &
          (ptr_patch%cells%neighbor_blk(ilc3,ibc3,1) == jb)) THEN
        inv_neighbor_id(jc,jb,3)  = 1
      ELSE IF ((ptr_patch%cells%neighbor_idx(ilc3,ibc3,2) == jc) .AND. &
               (ptr_patch%cells%neighbor_blk(ilc3,ibc3,2) == jb)) THEN
        inv_neighbor_id(jc,jb,3)  = 2
      ELSE IF ((ptr_patch%cells%neighbor_idx(ilc3,ibc3,3) == jc) .AND. &
               (ptr_patch%cells%neighbor_blk(ilc3,ibc3,3) == jb)) THEN
        inv_neighbor_id(jc,jb,3)  = 3
      ENDIF
    ENDIF

  ENDDO !cell loop

END DO !block loop
!$OMP END DO

! Compute coefficients for bilinear interpolation

!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx,yloc,xloc,pollat,pollon,&
!$OMP            ilc1,ibc1,ilc2,ibc2,ilc3,ibc3,xtemp,ytemp,wgt,x,y)
  DO jb = i_startblk, i_endblk

    CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, rl_start, rl_end)

    DO jc = i_startidx, i_endidx

      IF(.NOT.ptr_patch%cells%owner_mask(jc,jb)) CYCLE

      yloc = ptr_patch%cells%center(jc,jb)%lat
      xloc = ptr_patch%cells%center(jc,jb)%lon

      ! Rotate local point into the equator for better accuracy of bilinear weights
      IF (yloc >= 0._wp) THEN
        pollat = yloc - pi2/4._wp
      ELSE
        pollat = yloc + pi2/4._wp
      ENDIF
      pollon = xloc

      CALL rotate_latlon( yloc, xloc, pollat, pollon )

      ! line and block indices of the neighbouring cells

      ilc1 = ptr_patch%cells%neighbor_idx(jc,jb,1)
      ibc1 = ptr_patch%cells%neighbor_blk(jc,jb,1)
      ilc2 = ptr_patch%cells%neighbor_idx(jc,jb,2)
      ibc2 = ptr_patch%cells%neighbor_blk(jc,jb,2)
      ilc3 = ptr_patch%cells%neighbor_idx(jc,jb,3)
      ibc3 = ptr_patch%cells%neighbor_blk(jc,jb,3)

      ! x and y are the zonal and meridional distances from the local
      ! cell point (ignoring the earth's radius, which drops out anyway)

      xtemp = ptr_patch%cells%center(ilc1,ibc1)%lon
      ytemp = ptr_patch%cells%center(ilc1,ibc1)%lat
      CALL rotate_latlon( ytemp, xtemp, pollat, pollon )

      y(1)  = ytemp-yloc
      x(1)  = xtemp-xloc
      ! This is needed when the date line is crossed
      IF (x(1) >  3.5_wp) x(1) = x(1) - pi2
      IF (x(1) < -3.5_wp) x(1) = x(1) + pi2

      xtemp = ptr_patch%cells%center(ilc2,ibc2)%lon
      ytemp = ptr_patch%cells%center(ilc2,ibc2)%lat
      CALL rotate_latlon( ytemp, xtemp, pollat, pollon )

      y(2)  = ytemp-yloc
      x(2)  = xtemp-xloc
      ! This is needed when the date line is crossed
      IF (x(2) >  3.5_wp) x(2) = x(2) - pi2
      IF (x(2) < -3.5_wp) x(2) = x(2) + pi2

      xtemp = ptr_patch%cells%center(ilc3,ibc3)%lon
      ytemp = ptr_patch%cells%center(ilc3,ibc3)%lat
      CALL rotate_latlon( ytemp, xtemp, pollat, pollon )

      y(3)  = ytemp-yloc
      x(3)  = xtemp-xloc
      ! This is needed when the date line is crossed
      IF (x(3) >  3.5_wp) x(3) = x(3) - pi2
      IF (x(3) < -3.5_wp) x(3) = x(3) + pi2

      ! The weighting factors are based on the requirement that sum(w(i)*x(i)) = 0
      ! and sum(w(i)*y(i)) = 0, which ensures that linear horizontal gradients
      ! are not aliased into a checkerboard pattern between upward- and downward
      ! directed cells. The third condition is sum(w(i)) = 1., and the weight
      ! of the local point is 0.5 (see above). Analytical elimination yields...

      wgt(3) = 1._wp/( (y(3)-y(1)) - (x(3)-x(1))*(y(2)-y(1))/(x(2)-x(1)) ) * &
                  (1._wp-wgt_loc)*( -y(1) + x(1)*(y(2)-y(1))/(x(2)-x(1)) )
      wgt(2) = (-(1._wp-wgt_loc)*x(1) - wgt(3)*(x(3)-x(1)))/(x(2)-x(1))
      wgt(1) = 1._wp - wgt_loc - wgt(2) - wgt(3)

      ! Store results in ptr_patch%cells%avg_wgt
      ptr_int%c_bln_avg(jc,1,jb) = wgt_loc
      ptr_int%c_bln_avg(jc,2,jb) = wgt(1)
      ptr_int%c_bln_avg(jc,3,jb) = wgt(2)
      ptr_int%c_bln_avg(jc,4,jb) = wgt(3)

    ENDDO !cell loop

  END DO !block loop
!$OMP END DO
!$OMP END PARALLEL

z_inv_neighbor_id = inv_neighbor_id
CALL sync_patch_array(SYNC_C,ptr_patch,z_inv_neighbor_id(:,:,1))
CALL sync_patch_array(SYNC_C,ptr_patch,z_inv_neighbor_id(:,:,2))
CALL sync_patch_array(SYNC_C,ptr_patch,z_inv_neighbor_id(:,:,3))
inv_neighbor_id = z_inv_neighbor_id

CALL sync_patch_array(SYNC_C,ptr_patch,ptr_int%c_bln_avg)


! The coefficients for bilinear interpolation are now iteratively modified
! in order to obtain mass conservation.
! The criterion for conservation is that the three-point divergence
! calculated for any given grid point is used with a total factor of 1

DO iter = 1, niter

  ! Compute sum of weighting coefficients with which
  ! each local divergence value is used
  ! Note: the summation needs to be split into 4 loops in order to
  ! allow for vectorization and parallelization

  wgt_loc_sum = 0

  rl_start = 2
  rl_end = min_rlcell
  i_startblk = ptr_patch%cells%start_blk(rl_start,1)
  i_endblk   = ptr_patch%cells%end_blk(rl_end,i_nchdom)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx,ilc1,ibc1,ilc2,ibc2,ilc3,ibc3,inb1,inb2,inb3)
  DO jb = i_startblk, i_endblk

    CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, rl_start, rl_end)

    DO jc = i_startidx, i_endidx

      IF(.NOT. ptr_patch%cells%owner_mask(jc,jb)) CYCLE

      ilc1 = ptr_patch%cells%neighbor_idx(jc,jb,1)
      ibc1 = ptr_patch%cells%neighbor_blk(jc,jb,1)
      ilc2 = ptr_patch%cells%neighbor_idx(jc,jb,2)
      ibc2 = ptr_patch%cells%neighbor_blk(jc,jb,2)
      ilc3 = ptr_patch%cells%neighbor_idx(jc,jb,3)
      ibc3 = ptr_patch%cells%neighbor_blk(jc,jb,3)
      inb1 = inv_neighbor_id(jc,jb,1) + 1
      inb2 = inv_neighbor_id(jc,jb,2) + 1
      inb3 = inv_neighbor_id(jc,jb,3) + 1

      wgt_loc_sum(jc,jb) = &
        ptr_int%c_bln_avg(jc,1,jb)*ptr_patch%cells%area(jc,jb)          + &
        ptr_int%c_bln_avg(ilc1,inb1,ibc1)*ptr_patch%cells%area(ilc1,ibc1)  + &
        ptr_int%c_bln_avg(ilc2,inb2,ibc2)*ptr_patch%cells%area(ilc2,ibc2)  + &
        ptr_int%c_bln_avg(ilc3,inb3,ibc3)*ptr_patch%cells%area(ilc3,ibc3)

    ENDDO !cell loop

  END DO !block loop
!$OMP END DO
!$OMP END PARALLEL

  rl_start = 3
  i_startblk = ptr_patch%cells%start_blk(rl_start,1)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx)
  DO jb = i_startblk, i_endblk

    CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, rl_start, rl_end)

    DO jc = i_startidx, i_endidx

      IF(.NOT. ptr_patch%cells%owner_mask(jc,jb)) CYCLE

      ! For mass conservation, wgt_loc_sum/area should be 1 for each cell
      ! The deviation therefrom is termed residuum here.

      resid(jc,jb) = wgt_loc_sum(jc,jb)/ptr_patch%cells%area(jc,jb)-1._wp

    ENDDO !cell loop

  END DO !block loop
!$OMP END DO
!$OMP END PARALLEL

CALL sync_patch_array(SYNC_C,ptr_patch,resid)

IF (iter < niter) THEN ! Apply iterative correction to weighting coefficients
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx,ilc1,ibc1,ilc2,ibc2,ilc3,ibc3)
  DO jb = i_startblk, i_endblk

    CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, rl_start, rl_end)

    DO jc = i_startidx, i_endidx

      IF(.NOT. ptr_patch%cells%owner_mask(jc,jb)) CYCLE

      ! line and block indices of the neighbouring cells

      ilc1 = ptr_patch%cells%neighbor_idx(jc,jb,1)
      ibc1 = ptr_patch%cells%neighbor_blk(jc,jb,1)
      ilc2 = ptr_patch%cells%neighbor_idx(jc,jb,2)
      ibc2 = ptr_patch%cells%neighbor_blk(jc,jb,2)
      ilc3 = ptr_patch%cells%neighbor_idx(jc,jb,3)
      ibc3 = ptr_patch%cells%neighbor_blk(jc,jb,3)

      ! Modify weighting coefficients

      ptr_int%c_bln_avg(jc,1,jb) = ptr_int%c_bln_avg(jc,1,jb) - relax_coeff*resid(jc,jb)
      ptr_int%c_bln_avg(jc,2,jb) = ptr_int%c_bln_avg(jc,2,jb) - relax_coeff*resid(ilc1,ibc1)
      ptr_int%c_bln_avg(jc,3,jb) = ptr_int%c_bln_avg(jc,3,jb) - relax_coeff*resid(ilc2,ibc2)
      ptr_int%c_bln_avg(jc,4,jb) = ptr_int%c_bln_avg(jc,4,jb) - relax_coeff*resid(ilc3,ibc3)

      wgt_loc_sum(jc,jb) = SUM(ptr_int%c_bln_avg(jc,1:4,jb)) - 1._wp

      ptr_int%c_bln_avg(jc,1,jb) = ptr_int%c_bln_avg(jc,1,jb) - 0.25_wp*wgt_loc_sum(jc,jb)
      ptr_int%c_bln_avg(jc,2,jb) = ptr_int%c_bln_avg(jc,2,jb) - 0.25_wp*wgt_loc_sum(jc,jb)
      ptr_int%c_bln_avg(jc,3,jb) = ptr_int%c_bln_avg(jc,3,jb) - 0.25_wp*wgt_loc_sum(jc,jb)
      ptr_int%c_bln_avg(jc,4,jb) = ptr_int%c_bln_avg(jc,4,jb) - 0.25_wp*wgt_loc_sum(jc,jb)

      ! To be safe: Avoid runaway of central weight
      ptr_int%c_bln_avg(jc,1,jb) = MAX(ptr_int%c_bln_avg(jc,1,jb),minwgt_loc)
      ptr_int%c_bln_avg(jc,1,jb) = MIN(ptr_int%c_bln_avg(jc,1,jb),maxwgt_loc)

    ENDDO !cell loop

  END DO !block loop
!$OMP END DO
!$OMP END PARALLEL

  CALL sync_patch_array(SYNC_C,ptr_patch,ptr_int%c_bln_avg)

ELSE ! In the last iteration, enforce the mass conservation condition
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx)
  DO jb = i_startblk, i_endblk

    CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, rl_start, rl_end)

    DO jc = i_startidx, i_endidx

      IF(.NOT. ptr_patch%cells%owner_mask(jc,jb)) CYCLE

      ! Modify weighting coefficients

      ptr_int%c_bln_avg(jc,1,jb) = ptr_int%c_bln_avg(jc,1,jb) - resid(jc,jb)

    ENDDO !cell loop

  END DO !block loop
!$OMP END DO
!$OMP END PARALLEL

  CALL sync_patch_array(SYNC_C,ptr_patch,ptr_int%c_bln_avg)

! Compute coefficients needed to reconstruct averaged mass fluxes
! for approximately mass-consistent transport with divergence-averaging
! They can alternatively be used to average the velocity going into the divergence
! computation (without div averaging), yielding exact mass consistency but somewhat
! larger discretization errors for divergence

  rl_start = 4
  rl_end   = min_rledge
  i_startblk = ptr_patch%edges%start_blk(rl_start,1)
  i_endblk   = ptr_patch%edges%end_blk(rl_end,i_nchdom)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,je,i_startidx,i_endidx,ilc1,ibc1,ilc2,ibc2,inb1,inb2,inb3,ie4,ie5)
  DO jb = i_startblk, i_endblk

    CALL get_indices_e(ptr_patch, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, rl_start, rl_end)

    DO je = i_startidx, i_endidx

      IF(.NOT. ptr_patch%edges%owner_mask(je,jb)) CYCLE

      ilc1 = ptr_patch%edges%cell_idx(je,jb,1)
      ibc1 = ptr_patch%edges%cell_blk(je,jb,1)
      ilc2 = ptr_patch%edges%cell_idx(je,jb,2)
      ibc2 = ptr_patch%edges%cell_blk(je,jb,2)

      IF (je == ptr_patch%cells%edge_idx(ilc1,ibc1,1) .AND. &
          jb == ptr_patch%cells%edge_blk(ilc1,ibc1,1)) THEN

        inb1 = inv_neighbor_id(ilc1,ibc1,1)
        ie4  = MOD(inb1,  3)+1
        ie5  = MOD(inb1+1,3)+1

        ptr_int%e_flx_avg(je,2,jb) = ptr_int%c_bln_avg(ilc2,inb1+1,ibc2) &
          *ptr_int%geofac_div(ilc1,2,ibc1)/ptr_int%geofac_div(ilc2,inb1,ibc2)

        ptr_int%e_flx_avg(je,3,jb) = ptr_int%c_bln_avg(ilc2,inb1+1,ibc2) &
          *ptr_int%geofac_div(ilc1,3,ibc1)/ptr_int%geofac_div(ilc2,inb1,ibc2)

        ptr_int%e_flx_avg(je,4,jb) = ptr_int%c_bln_avg(ilc1,2,ibc1) &
          *ptr_int%geofac_div(ilc2,ie4,ibc2)/ptr_int%geofac_div(ilc1,1,ibc1)

        ptr_int%e_flx_avg(je,5,jb) = ptr_int%c_bln_avg(ilc1,2,ibc1) &
          *ptr_int%geofac_div(ilc2,ie5,ibc2)/ptr_int%geofac_div(ilc1,1,ibc1)


      ELSE IF (je == ptr_patch%cells%edge_idx(ilc1,ibc1,2) .AND. &
               jb == ptr_patch%cells%edge_blk(ilc1,ibc1,2)) THEN

        inb2 = inv_neighbor_id(ilc1,ibc1,2)
        ie4  = MOD(inb2  ,3)+1
        ie5  = MOD(inb2+1,3)+1

        ptr_int%e_flx_avg(je,2,jb) = ptr_int%c_bln_avg(ilc2,inb2+1,ibc2) &
          *ptr_int%geofac_div(ilc1,3,ibc1)/ptr_int%geofac_div(ilc2,inb2,ibc2)

        ptr_int%e_flx_avg(je,3,jb) = ptr_int%c_bln_avg(ilc2,inb2+1,ibc2) &
          *ptr_int%geofac_div(ilc1,1,ibc1)/ptr_int%geofac_div(ilc2,inb2,ibc2)

        ptr_int%e_flx_avg(je,4,jb) = ptr_int%c_bln_avg(ilc1,3,ibc1) &
          *ptr_int%geofac_div(ilc2,ie4,ibc2)/ptr_int%geofac_div(ilc1,2,ibc1)

        ptr_int%e_flx_avg(je,5,jb) = ptr_int%c_bln_avg(ilc1,3,ibc1) &
          *ptr_int%geofac_div(ilc2,ie5,ibc2)/ptr_int%geofac_div(ilc1,2,ibc1)


      ELSE IF (je == ptr_patch%cells%edge_idx(ilc1,ibc1,3) .AND. &
               jb == ptr_patch%cells%edge_blk(ilc1,ibc1,3)) THEN

        inb3 = inv_neighbor_id(ilc1,ibc1,3)
        ie4  = MOD(inb3  ,3)+1
        ie5  = MOD(inb3+1,3)+1

        ptr_int%e_flx_avg(je,2,jb) = ptr_int%c_bln_avg(ilc2,inb3+1,ibc2) &
          *ptr_int%geofac_div(ilc1,1,ibc1)/ptr_int%geofac_div(ilc2,inb3,ibc2)

        ptr_int%e_flx_avg(je,3,jb) = ptr_int%c_bln_avg(ilc2,inb3+1,ibc2) &
          *ptr_int%geofac_div(ilc1,2,ibc1)/ptr_int%geofac_div(ilc2,inb3,ibc2)

        ptr_int%e_flx_avg(je,4,jb) = ptr_int%c_bln_avg(ilc1,4,ibc1) &
          *ptr_int%geofac_div(ilc2,ie4,ibc2)/ptr_int%geofac_div(ilc1,3,ibc1)

        ptr_int%e_flx_avg(je,5,jb) = ptr_int%c_bln_avg(ilc1,4,ibc1) &
          *ptr_int%geofac_div(ilc2,ie5,ibc2)/ptr_int%geofac_div(ilc1,3,ibc1)

      ENDIF

    ENDDO !edge loop

  END DO !block loop
!$OMP END DO
!$OMP END PARALLEL

  CALL sync_patch_array(SYNC_E,ptr_patch,ptr_int%e_flx_avg)

  rl_start = 5
  i_startblk = ptr_patch%edges%start_blk(rl_start,1)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,je,i_startidx,i_endidx,ilc1,ibc1,ilc2,ibc2,inb1,inb2,inb3,ie4,ie5, &
!$OMP            ile1,ibe1,ile2,ibe2,ile3,ibe3,ile4,ibe4,iie1,iie2,iie3,iie4)
  DO jb = i_startblk, i_endblk

    CALL get_indices_e(ptr_patch, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, rl_start, rl_end)

    DO je = i_startidx, i_endidx

      IF(.NOT. ptr_patch%edges%owner_mask(je,jb)) CYCLE

      ilc1 = ptr_patch%edges%cell_idx(je,jb,1)
      ibc1 = ptr_patch%edges%cell_blk(je,jb,1)
      ilc2 = ptr_patch%edges%cell_idx(je,jb,2)
      ibc2 = ptr_patch%edges%cell_blk(je,jb,2)

      ile1 = ptr_patch%edges%quad_idx(je,jb,1)
      ibe1 = ptr_patch%edges%quad_blk(je,jb,1)
      ile2 = ptr_patch%edges%quad_idx(je,jb,2)
      ibe2 = ptr_patch%edges%quad_blk(je,jb,2)
      ile3 = ptr_patch%edges%quad_idx(je,jb,3)
      ibe3 = ptr_patch%edges%quad_blk(je,jb,3)
      ile4 = ptr_patch%edges%quad_idx(je,jb,4)
      ibe4 = ptr_patch%edges%quad_blk(je,jb,4)

      IF (ptr_patch%edges%cell_idx(ile1,ibe1,1) == ilc1 .AND. &
          ptr_patch%edges%cell_blk(ile1,ibe1,1) == ibc1 ) THEN
        iie1(je) = 3
      ELSE IF (ptr_patch%edges%cell_idx(ile1,ibe1,2) == ilc1 .AND. &
               ptr_patch%edges%cell_blk(ile1,ibe1,2) == ibc1 ) THEN
        iie1(je) = 5
      ENDIF
      IF (ptr_patch%edges%cell_idx(ile2,ibe2,1) == ilc1 .AND. &
          ptr_patch%edges%cell_blk(ile2,ibe2,1) == ibc1 ) THEN
        iie2(je) = 2
      ELSE IF (ptr_patch%edges%cell_idx(ile2,ibe2,2) == ilc1 .AND. &
               ptr_patch%edges%cell_blk(ile2,ibe2,2) == ibc1 ) THEN
        iie2(je) = 4
      ENDIF
      IF (ptr_patch%edges%cell_idx(ile3,ibe3,1) == ilc2 .AND. &
          ptr_patch%edges%cell_blk(ile3,ibe3,1) == ibc2 ) THEN
        iie3(je) = 3
      ELSE IF (ptr_patch%edges%cell_idx(ile3,ibe3,2) == ilc2 .AND. &
               ptr_patch%edges%cell_blk(ile3,ibe3,2) == ibc2 ) THEN
        iie3(je) = 5
      ENDIF
      IF (ptr_patch%edges%cell_idx(ile4,ibe4,1) == ilc2 .AND. &
          ptr_patch%edges%cell_blk(ile4,ibe4,1) == ibc2 ) THEN
        iie4(je) = 2
      ELSE IF (ptr_patch%edges%cell_idx(ile4,ibe4,2) == ilc2 .AND. &
               ptr_patch%edges%cell_blk(ile4,ibe4,2) == ibc2 ) THEN
        iie4(je) = 4
      ENDIF

      IF (je == ptr_patch%cells%edge_idx(ilc1,ibc1,1) .AND. &
          jb == ptr_patch%cells%edge_blk(ilc1,ibc1,1)) THEN

        inb1 = inv_neighbor_id(ilc1,ibc1,1)
        ie4  = MOD(inb1  ,3)+1
        ie5  = MOD(inb1+1,3)+1

        ptr_int%e_flx_avg(je,1,jb) = 0.5_wp *(                                        &
          ( ptr_int%geofac_div(ilc1,1,ibc1)*ptr_int%c_bln_avg(ilc1,1,ibc1)            &
          + ptr_int%geofac_div(ilc2,inb1,ibc2)*ptr_int%c_bln_avg(ilc1,2,ibc1)         &
          - ptr_int%e_flx_avg(ile1,iie1(je),ibe1)*ptr_int%geofac_div(ilc1,2,ibc1)     &
          - ptr_int%e_flx_avg(ile2,iie2(je),ibe2)*ptr_int%geofac_div(ilc1,3,ibc1) )   &
          / ptr_int%geofac_div(ilc1,1,ibc1)   +                                       &
          ( ptr_int%geofac_div(ilc2,inb1,ibc2)*ptr_int%c_bln_avg(ilc2,1,ibc2)         &
          + ptr_int%geofac_div(ilc1,1,ibc1)*ptr_int%c_bln_avg(ilc2,inb1+1,ibc2)       &
          - ptr_int%e_flx_avg(ile3,iie3(je),ibe3)*ptr_int%geofac_div(ilc2,ie4,ibc2)   &
          - ptr_int%e_flx_avg(ile4,iie4(je),ibe4)*ptr_int%geofac_div(ilc2,ie5,ibc2) ) &
          / ptr_int%geofac_div(ilc2,inb1,ibc2) )


      ELSE IF (je == ptr_patch%cells%edge_idx(ilc1,ibc1,2) .AND. &
               jb == ptr_patch%cells%edge_blk(ilc1,ibc1,2)) THEN

        inb2 = inv_neighbor_id(ilc1,ibc1,2)
        ie4  = MOD(inb2  ,3)+1
        ie5  = MOD(inb2+1,3)+1

        ptr_int%e_flx_avg(je,1,jb) = 0.5_wp *(                                        &
          ( ptr_int%geofac_div(ilc1,2,ibc1)*ptr_int%c_bln_avg(ilc1,1,ibc1)            &
          + ptr_int%geofac_div(ilc2,inb2,ibc2)*ptr_int%c_bln_avg(ilc1,3,ibc1)         &
          - ptr_int%e_flx_avg(ile1,iie1(je),ibe1)*ptr_int%geofac_div(ilc1,3,ibc1)     &
          - ptr_int%e_flx_avg(ile2,iie2(je),ibe2)*ptr_int%geofac_div(ilc1,1,ibc1) )   &
          / ptr_int%geofac_div(ilc1,2,ibc1)   +                                       &
          ( ptr_int%geofac_div(ilc2,inb2,ibc2)*ptr_int%c_bln_avg(ilc2,1,ibc2)         &
          + ptr_int%geofac_div(ilc1,2,ibc1)*ptr_int%c_bln_avg(ilc2,inb2+1,ibc2)       &
          - ptr_int%e_flx_avg(ile3,iie3(je),ibe3)*ptr_int%geofac_div(ilc2,ie4,ibc2)   &
          - ptr_int%e_flx_avg(ile4,iie4(je),ibe4)*ptr_int%geofac_div(ilc2,ie5,ibc2) ) &
          / ptr_int%geofac_div(ilc2,inb2,ibc2) )


      ELSE IF (je == ptr_patch%cells%edge_idx(ilc1,ibc1,3) .AND. &
               jb == ptr_patch%cells%edge_blk(ilc1,ibc1,3)) THEN

        inb3 = inv_neighbor_id(ilc1,ibc1,3)
        ie4  = MOD(inb3  ,3)+1
        ie5  = MOD(inb3+1,3)+1

        ptr_int%e_flx_avg(je,1,jb) = 0.5_wp *(                                        &
          ( ptr_int%geofac_div(ilc1,3,ibc1)*ptr_int%c_bln_avg(ilc1,1,ibc1)            &
          + ptr_int%geofac_div(ilc2,inb3,ibc2)*ptr_int%c_bln_avg(ilc1,4,ibc1)         &
          - ptr_int%e_flx_avg(ile1,iie1(je),ibe1)*ptr_int%geofac_div(ilc1,1,ibc1)     &
          - ptr_int%e_flx_avg(ile2,iie2(je),ibe2)*ptr_int%geofac_div(ilc1,2,ibc1) )   &
          / ptr_int%geofac_div(ilc1,3,ibc1)   +                                       &
          ( ptr_int%geofac_div(ilc2,inb3,ibc2)*ptr_int%c_bln_avg(ilc2,1,ibc2)         &
          + ptr_int%geofac_div(ilc1,3,ibc1)*ptr_int%c_bln_avg(ilc2,inb3+1,ibc2)       &
          - ptr_int%e_flx_avg(ile3,iie3(je),ibe3)*ptr_int%geofac_div(ilc2,ie4,ibc2)   &
          - ptr_int%e_flx_avg(ile4,iie4(je),ibe4)*ptr_int%geofac_div(ilc2,ie5,ibc2) ) &
          / ptr_int%geofac_div(ilc2,inb3,ibc2) )

      ENDIF

    ENDDO !edge loop

  END DO !block loop
!$OMP END DO

  ! Finally, the weighting coefficients are scaled in order to
  ! yield the right result for a constant wind field

!$OMP DO PRIVATE(jb,je,i_startidx,i_endidx,ile1,ibe1,ile2,ibe2,ile3,ibe3,ile4,ibe4, &
!$OMP            z_nx1,z_nx2,z_nx3,z_nx4,z_nx5)
  DO jb = i_startblk, i_endblk

    CALL get_indices_e(ptr_patch, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, rl_start, rl_end)

    DO je = i_startidx, i_endidx

      IF(.NOT. ptr_patch%edges%owner_mask(je,jb)) CYCLE

      ile1 = ptr_patch%edges%quad_idx(je,jb,1)
      ibe1 = ptr_patch%edges%quad_blk(je,jb,1)
      ile2 = ptr_patch%edges%quad_idx(je,jb,2)
      ibe2 = ptr_patch%edges%quad_blk(je,jb,2)
      ile3 = ptr_patch%edges%quad_idx(je,jb,3)
      ibe3 = ptr_patch%edges%quad_blk(je,jb,3)
      ile4 = ptr_patch%edges%quad_idx(je,jb,4)
      ibe4 = ptr_patch%edges%quad_blk(je,jb,4)

      z_nx1(je,1:3) = ptr_patch%edges%primal_cart_normal(je,jb)%x(1:3)
      z_nx2(je,1:3) = ptr_patch%edges%primal_cart_normal(ile1,ibe1)%x(1:3)
      z_nx3(je,1:3) = ptr_patch%edges%primal_cart_normal(ile2,ibe2)%x(1:3)
      z_nx4(je,1:3) = ptr_patch%edges%primal_cart_normal(ile3,ibe3)%x(1:3)
      z_nx5(je,1:3) = ptr_patch%edges%primal_cart_normal(ile4,ibe4)%x(1:3)

      ! The sum of the coefficients - multiplied by the projection factors -
      ! is enforced to be 1 so that a constant vector field is processed correctly

      checksum(je,jb) = ptr_int%e_flx_avg(je,1,jb)                            &
        + DOT_PRODUCT(z_nx1(je,1:3),z_nx2(je,1:3))*ptr_int%e_flx_avg(je,2,jb) &
        + DOT_PRODUCT(z_nx1(je,1:3),z_nx3(je,1:3))*ptr_int%e_flx_avg(je,3,jb) &
        + DOT_PRODUCT(z_nx1(je,1:3),z_nx4(je,1:3))*ptr_int%e_flx_avg(je,4,jb) &
        + DOT_PRODUCT(z_nx1(je,1:3),z_nx5(je,1:3))*ptr_int%e_flx_avg(je,5,jb)

      ptr_int%e_flx_avg(je,1,jb) = ptr_int%e_flx_avg(je,1,jb)/checksum(je,jb)
      ptr_int%e_flx_avg(je,2,jb) = ptr_int%e_flx_avg(je,2,jb)/checksum(je,jb)
      ptr_int%e_flx_avg(je,3,jb) = ptr_int%e_flx_avg(je,3,jb)/checksum(je,jb)
      ptr_int%e_flx_avg(je,4,jb) = ptr_int%e_flx_avg(je,4,jb)/checksum(je,jb)
      ptr_int%e_flx_avg(je,5,jb) = ptr_int%e_flx_avg(je,5,jb)/checksum(je,jb)

    ENDDO !edge loop

  END DO !block loop
!$OMP END DO
!$OMP END PARALLEL

  CALL sync_patch_array(SYNC_C,ptr_patch,ptr_int%c_bln_avg)
  CALL sync_patch_array(SYNC_E,ptr_patch,ptr_int%e_flx_avg)

ENDIF ! end of last iteration
ENDDO ! iteration loop

! Optional debug output for bilinear averaging coefficients
#ifdef DEBUG_COEFF

rl_start = 2
rl_end = min_rlcell
i_startblk = ptr_patch%cells%start_blk(rl_start,1)
i_endblk   = ptr_patch%cells%end_blk(rl_end,i_nchdom)

sum1 = 0._wp
wgt_loc_sum = 1._wp

DO jb = i_startblk, i_endblk

  CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
                     i_startidx, i_endidx, rl_start, rl_end)

  DO jc = i_startidx, i_endidx

    sum1 = sum1 + resid(jc,jb)**2
    wgt_loc_sum(jc,jb) = SUM(ptr_int%c_bln_avg(jc,1:4,jb))

    WRITE(710+ptr_patch%id,'(2i5,5f12.6,e13.5)') jb,jc,ptr_int%c_bln_avg(jc,1:4,jb),&
       wgt_loc_sum(jc,jb),resid(jc,jb)

  END DO
END DO
WRITE(710+ptr_patch%id,'(4e13.5)') MAXVAL(resid),SQRT(sum1/ptr_patch%n_patch_cells),&
                                   MAXVAL(wgt_loc_sum)-1._wp,MINVAL(wgt_loc_sum)-1._wp
CLOSE (710+ptr_patch%id)

! Debug output for mass flux averaging weights

rl_start = 5
rl_end = min_rledge
i_startblk = ptr_patch%edges%start_blk(rl_start,1)
i_endblk   = ptr_patch%edges%end_blk(rl_end,i_nchdom)

DO jb = i_startblk, i_endblk

  CALL get_indices_e(ptr_patch, jb, i_startblk, i_endblk, &
                     i_startidx, i_endidx, rl_start, rl_end)

  DO je = i_startidx, i_endidx

    WRITE(720+ptr_patch%id,'(2i5,6f12.6)') jb,je,ptr_int%e_flx_avg(je,1:5,jb),&
       checksum(je,jb)

  END DO
END DO
CLOSE (720+ptr_patch%id)

#endif

END SUBROUTINE init_cellavg_wgt


!-------------------------------------------------------------------------
!
!
!>
!! Computes the coefficients for lateral boundary nudging needed for
!! one-way nesting and the limited-area mode
!! The nudging coefficients are defined via three namelist variables:
!! nudge_max_coeff: Maximum relaxation coefficient in the cell row bordering to
!! the boundary interpolation zone
!! nudge_efold_width: e-folding width of exponential decay of coefficients
!! (in units of grid cell rows)
!! nudge_zone_width: Total width of nudging zone (in units of grid cell rows)
!!
!! @par Revision History
!!  developed by Guenther Zaengl, 2010-06-21
!!
SUBROUTINE init_nudgecoeffs( ptr_patch, ptr_int )
!
!
!  patch on which computation is performed
!
TYPE(t_patch), TARGET, INTENT(inout) :: ptr_patch

! Interpolation state
TYPE(t_int_state), INTENT(inout):: ptr_int
!

INTEGER :: jc, je, jb
INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx, i_nchdom

INTEGER :: max_rlval

!-----------------------------------------------------------------------

! Check if required refin_ctrl information is available
max_rlval = MAXVAL(ptr_patch%cells%refin_ctrl(:,:))
max_rlval = NINT(global_max(REAL(max_rlval,wp)))

IF (max_rlval < nudge_zone_width+grf_nudge_start_c-1) THEN
    CALL finish('init_nudgecoeffs',&
                'bdy_indexing_depth in prepare_gridref must be at least nudge_zone_width+4')
ENDIF

i_nchdom   = MAX(1,ptr_patch%n_childdom)

!$OMP PARALLEL PRIVATE(i_startblk,i_endblk)

! a) Nudging coefficients for cells
i_startblk = ptr_patch%cells%start_blk(grf_nudge_start_c,1)
i_endblk   = ptr_patch%cells%end_blk(0,i_nchdom)

!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx)
DO jb = i_startblk, i_endblk

  CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
                     i_startidx, i_endidx, grf_nudge_start_c, 0)

  DO jc = i_startidx, i_endidx

    IF (ptr_patch%cells%refin_ctrl(jc,jb) > 0 .AND. &
        ptr_patch%cells%refin_ctrl(jc,jb) <= nudge_zone_width+grf_nudge_start_c-1) THEN
      ptr_int%nudgecoeff_c(jc,jb) = &
        nudge_max_coeff*EXP(-REAL(ptr_patch%cells%refin_ctrl(jc,jb)-grf_nudge_start_c,wp) / &
                            nudge_efold_width)
    ENDIF

  ENDDO !cell loop

END DO !block loop
!$OMP END DO

! b) Nudging coefficients for edges
i_startblk = ptr_patch%edges%start_blk(grf_nudge_start_e,1)
i_endblk   = ptr_patch%edges%end_blk(0,i_nchdom)

IF (ptr_patch%id > 1 .AND. lfeedback(ptr_patch%id)) THEN
  ! Use nudging coefficients optimized for velocity boundary diffusion
!$OMP DO PRIVATE(jb,je,i_startidx,i_endidx)
  DO jb = i_startblk, i_endblk

    CALL get_indices_e(ptr_patch, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, grf_nudge_start_e, 0)

    DO je = i_startidx, i_endidx

      IF (ptr_patch%edges%refin_ctrl(je,jb) > 0 .AND. &
          ptr_patch%edges%refin_ctrl(je,jb) <= grf_nudge_start_e+9) THEN
        ptr_int%nudgecoeff_e(je,jb) = nudge_max_coeff* &
          EXP(-REAL(ptr_patch%edges%refin_ctrl(je,jb)-grf_nudge_start_e,wp) / 4._wp)
      ENDIF

    ENDDO !edge loop

  END DO !block loop
!$OMP END DO
ELSE
  ! Use nudging coefficients from namelist
!$OMP DO PRIVATE(jb,je,i_startidx,i_endidx)
  DO jb = i_startblk, i_endblk

    CALL get_indices_e(ptr_patch, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, grf_nudge_start_e, 0)

    DO je = i_startidx, i_endidx

      IF (ptr_patch%edges%refin_ctrl(je,jb) > 0 .AND. &
          ptr_patch%edges%refin_ctrl(je,jb) <= 2*nudge_zone_width+grf_nudge_start_e-3) THEN
        ptr_int%nudgecoeff_e(je,jb) = &
          nudge_max_coeff*EXP(-REAL(ptr_patch%edges%refin_ctrl(je,jb)-grf_nudge_start_e,wp) / &
                              (2._wp*nudge_efold_width))
      ENDIF

    ENDDO !edge loop

  END DO !block loop
!$OMP END DO
ENDIF

!$OMP END PARALLEL

CALL sync_patch_array(SYNC_C,ptr_patch,ptr_int%nudgecoeff_c)
CALL sync_patch_array(SYNC_E,ptr_patch,ptr_int%nudgecoeff_e)

END SUBROUTINE init_nudgecoeffs

!-------------------------------------------------------------------------
!
!
!>
!! Precomputes the geometrical factors used in the divergence, rotation.
!!
!! Precomputes the geometrical factors used in the divergence, rotation
!! and nabla_2_scalar operators in order to improve computational efficiency
!!
!! @par Revision History
!!  developed by Guenther Zaengl, 2009-03-17
!!  Modification by Almut Gassmann, 2009-12-19
!!  - Vorticity is computed on quads in case of the hexagonal grid
!!  Modification by Almut Gassmann, 2010-02-05
!!  - Added feature for poor men's 3rd order advection, where a directional
!!    laplace is needed at the edges.
!!  Modification by Almut Gassmann, 2010-06-11
!!  - Further works for providing directional gradients of velocities
!!
SUBROUTINE init_geo_factors( ptr_patch, ptr_int )
!
IMPLICIT NONE
!
!  patch on which computation is performed
!
TYPE(t_patch), TARGET, INTENT(inout) :: ptr_patch

! Interpolation state
TYPE(t_int_state), INTENT(inout):: ptr_int
!

INTEGER :: jc, jb, je, jv, je1, jn,jm, nincr
INTEGER :: rl_start, rl_end
INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx, i_nchdom

INTEGER :: ile, ibe, ilc1, ibc1, ilc2, ibc2, ifac, ic, ilnc, ibnc, &
           ilv1, ilv2, ibv1, ibv2, jr1, ilr, ibr
TYPE(t_cartesian_coordinates)::z_pn_k,z_pn_j,z_pt_k,z_cart_no,z_cart_ea
REAL(wp) :: z_proj, z_norm, z_lon, z_lat


!-----------------------------------------------------------------------

i_nchdom   = MAX(1,ptr_patch%n_childdom)


!$OMP PARALLEL PRIVATE(rl_start,rl_end,i_startblk,i_endblk,ifac)
! a) Geometrical factor for divergence
rl_start = 1
rl_end = min_rlcell

! values for the blocking
i_startblk = ptr_patch%cells%start_blk(rl_start,1)
i_endblk   = ptr_patch%cells%end_blk(rl_end,i_nchdom)
!
! loop through all patch cells (and blocks)
!
!$OMP DO PRIVATE(jb,je,jc,i_startidx,i_endidx,ile,ibe)
DO jb = i_startblk, i_endblk

  CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
                     i_startidx, i_endidx, rl_start, rl_end)

  DO je = 1, ptr_patch%cell_type
    DO jc = i_startidx, i_endidx

      IF (je > ptr_patch%cells%num_edges(jc,jb)) CYCLE ! relevant for hexagons

      ile = ptr_patch%cells%edge_idx(jc,jb,je)
      ibe = ptr_patch%cells%edge_blk(jc,jb,je)

      ptr_int%geofac_div(jc,je,jb) = &
        ptr_patch%edges%primal_edge_length(ile,ibe) * &
        ptr_patch%cells%edge_orientation(jc,jb,je)  / &
        ptr_patch%cells%area(jc,jb)

    ENDDO !cell loop
  ENDDO

END DO !block loop
!$OMP END DO

! b) Geometrical factor for rotation
rl_start = 2
rl_end = min_rlvert

! Vorticity should have the right sign
  SELECT CASE (ptr_patch%cell_type)
  CASE (3)
    ifac = 1
  CASE (6)
    ifac = -1
  END SELECT
  ! values for the blocking
  i_startblk = ptr_patch%verts%start_blk(rl_start,1)
  i_endblk   = ptr_patch%verts%end_blk(rl_end,i_nchdom)
  !
  ! loop through all patch cells (and blocks)
  !
!$OMP DO PRIVATE(jb,je,jv,i_startidx,i_endidx,ile,ibe)
  DO jb = i_startblk, i_endblk

    CALL get_indices_v(ptr_patch, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, rl_start, rl_end)

    DO je = 1, 9-ptr_patch%cell_type
      DO jv = i_startidx, i_endidx

        IF(.NOT. ptr_patch%verts%owner_mask(jv,jb)) CYCLE

        IF (je > ptr_patch%verts%num_edges(jv,jb)) CYCLE

        ile = ptr_patch%verts%edge_idx(jv,jb,je)
        ibe = ptr_patch%verts%edge_blk(jv,jb,je)

        ptr_int%geofac_rot(jv,je,jb) =                &
          ptr_patch%edges%dual_edge_length(ile,ibe) * &
          ptr_patch%verts%edge_orientation(jv,jb,je)/ &
          ptr_patch%verts%dual_area(jv,jb) * REAL(ifac,wp)

      ENDDO !vertex loop
    ENDDO

  END DO !block loop
!$OMP END DO

! c) Geometrical factor for nabla2_scalar
rl_start = 2
rl_end = min_rlcell

! values for the blocking
i_startblk = ptr_patch%cells%start_blk(rl_start,1)
i_endblk   = ptr_patch%cells%end_blk(rl_end,i_nchdom)
!
! loop through all patch cells (and blocks)
!
!$OMP DO PRIVATE(jb,je,jc,ic,i_startidx,i_endidx,ile,ibe,ilc1,ibc1,&
!$OMP    ilc2,ibc2,ilnc,ibnc)
DO jb = i_startblk, i_endblk

  CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
                     i_startidx, i_endidx, rl_start, rl_end)

  DO je = 1, ptr_patch%cell_type
    DO jc = i_startidx, i_endidx

      ile = ptr_patch%cells%edge_idx(jc,jb,je)
      ibe = ptr_patch%cells%edge_blk(jc,jb,je)

      ilc1 = ptr_patch%edges%cell_idx(ile,ibe,1)
      ibc1 = ptr_patch%edges%cell_blk(ile,ibe,1)
      ilc2 = ptr_patch%edges%cell_idx(ile,ibe,2)
      ibc2 = ptr_patch%edges%cell_blk(ile,ibe,2)

      IF (jc == ilc1 .AND. jb == ibc1) THEN
        IF (ptr_patch%cell_type == 3) THEN
          ptr_int%geofac_n2s(jc,1,jb) = ptr_int%geofac_n2s(jc,1,jb) - &
            ptr_int%geofac_div(jc,je,jb) /                            &
            ptr_patch%edges%dual_edge_length(ile,ibe)
        ELSE IF (ptr_patch%cell_type == 6) THEN
          ptr_int%geofac_n2s(jc,1,jb) = ptr_int%geofac_n2s(jc,1,jb) - &
            ptr_int%geofac_div(jc,je,jb) /                            &
            ptr_patch%edges%dual_edge_length(ile,ibe)*                &
            ptr_patch%edges%system_orientation(ile,ibe)
        ENDIF
      ELSE IF (jc == ilc2 .AND. jb == ibc2) THEN
        IF (ptr_patch%cell_type == 3) THEN
          ptr_int%geofac_n2s(jc,1,jb) = ptr_int%geofac_n2s(jc,1,jb) + &
            ptr_int%geofac_div(jc,je,jb) /                            &
            ptr_patch%edges%dual_edge_length(ile,ibe)
        ELSE IF (ptr_patch%cell_type == 6) THEN
          ptr_int%geofac_n2s(jc,1,jb) = ptr_int%geofac_n2s(jc,1,jb) + &
            ptr_int%geofac_div(jc,je,jb) /                            &
            ptr_patch%edges%dual_edge_length(ile,ibe)*                &
            ptr_patch%edges%system_orientation(ile,ibe)
        ENDIF
      ENDIF
      DO ic = 1, ptr_patch%cell_type
        ilnc = ptr_patch%cells%neighbor_idx(jc,jb,ic)
        ibnc = ptr_patch%cells%neighbor_blk(jc,jb,ic)
        IF (ilnc == ilc1 .AND. ibnc == ibc1) THEN
          IF (ptr_patch%cell_type == 3) THEN
            ptr_int%geofac_n2s(jc,ic+1,jb) = ptr_int%geofac_n2s(jc,ic+1,jb) - &
              ptr_int%geofac_div(jc,je,jb) /                                  &
              ptr_patch%edges%dual_edge_length(ile,ibe)
          ELSE IF (ptr_patch%cell_type == 6) THEN
            ptr_int%geofac_n2s(jc,ic+1,jb) = ptr_int%geofac_n2s(jc,ic+1,jb) - &
              ptr_int%geofac_div(jc,je,jb) /                                  &
              ptr_patch%edges%dual_edge_length(ile,ibe)*                      &
              ptr_patch%edges%system_orientation(ile,ibe)
          ENDIF
        ELSE IF (ilnc == ilc2 .AND. ibnc == ibc2) THEN
          IF (ptr_patch%cell_type == 3) THEN
            ptr_int%geofac_n2s(jc,ic+1,jb) = ptr_int%geofac_n2s(jc,ic+1,jb) + &
              ptr_int%geofac_div(jc,je,jb) /                                  &
              ptr_patch%edges%dual_edge_length(ile,ibe)
          ELSE IF (ptr_patch%cell_type == 6) THEN
            ptr_int%geofac_n2s(jc,ic+1,jb) = ptr_int%geofac_n2s(jc,ic+1,jb) + &
              ptr_int%geofac_div(jc,je,jb) /                                  &
              ptr_patch%edges%dual_edge_length(ile,ibe)*                      &
              ptr_patch%edges%system_orientation(ile,ibe)
          ENDIF
        ENDIF
      ENDDO

      ! To ensure that dummy edges have a factor of 0:
      IF (je > ptr_patch%cells%num_edges(jc,jb)) THEN
        ptr_int%geofac_n2s(jc,je+1,jb) = 0._wp
      ENDIF

    ENDDO !cell loop
  ENDDO

END DO !block loop
!$OMP END DO

! d) Geometrical factor for quad-cell divergence (triangles only)

IF (ptr_patch%cell_type == 3) THEN

  rl_start = 2
  rl_end = min_rledge

  ! values for the blocking
  i_startblk = ptr_patch%edges%start_blk(rl_start,1)
  i_endblk   = ptr_patch%edges%end_blk(rl_end,i_nchdom)

!$OMP DO PRIVATE(jb,je,je1,i_startidx,i_endidx,ile,ibe)
  DO jb = i_startblk, i_endblk

    CALL get_indices_e(ptr_patch, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, rl_start, rl_end)

    DO je1 = 1, 4
      DO je = i_startidx, i_endidx

      IF(.NOT. ptr_patch%edges%owner_mask(je,jb)) CYCLE

      ile = ptr_patch%edges%quad_idx(je,jb,je1)
      ibe = ptr_patch%edges%quad_blk(je,jb,je1)

      ptr_int%geofac_qdiv(je,je1,jb) = &
        ptr_patch%edges%primal_edge_length(ile,ibe) * &
        ptr_patch%edges%quad_orientation(je,jb,je1)  / &
        ptr_patch%edges%quad_area(je,jb)

      ENDDO !edge loop
    ENDDO

  END DO !block loop
!$OMP END DO

ENDIF

    ! e) coefficients for directional gradient of a normal vector quantity
    ! at the same edge (gives directional laplacian if gradient psi is assumed as input)

    IF (ptr_patch%cell_type == 6) THEN

      ! Now compute coefficients
      !-------------------------
      rl_start = 2
      rl_end = min_rledge
      ! values for the blocking
      i_startblk = ptr_patch%edges%start_blk(rl_start,1)
      i_endblk   = ptr_patch%edges%end_blk(rl_end,i_nchdom)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,je,z_pn_k,z_pt_k,&
!$OMP ilc1,ibc1,ilc2,ibc2,ilv1,ibv1,ilv2,ibv2,je1,ile,ibe,z_proj,z_pn_j,&
!$OMP jm,jn,nincr,ilr,ibr,jr1,z_lon,z_lat,z_norm,z_cart_no,z_cart_ea)
      DO jb = i_startblk, i_endblk
        CALL get_indices_e(ptr_patch, jb, i_startblk, i_endblk, &
        &                  i_startidx, i_endidx, rl_start, rl_end)

        DO je = i_startidx, i_endidx

          IF(.NOT. ptr_patch%edges%owner_mask(je,jb)) CYCLE

          ! transform primal normal to cartesian vector z_pn_k
          z_pn_k%x=ptr_patch%edges%primal_cart_normal(je,jb)%x

          ! transform primal tangential (=dual normal) to cartesian vector z_pt_k
          z_pt_k%x=ptr_patch%edges%dual_cart_normal(je,jb)%x

          ! determine neighbor cells
          ilc1 = ptr_patch%edges%cell_idx(je,jb,1)
          ibc1 = ptr_patch%edges%cell_blk(je,jb,1)
          ilc2 = ptr_patch%edges%cell_idx(je,jb,2)
          ibc2 = ptr_patch%edges%cell_blk(je,jb,2)

          ! neighbor cell 1
          ! cartesian vector of cell center
          IF (.NOT.lplane) THEN
            z_lon = ptr_patch%cells%center(ilc1,ibc1)%lon
            z_lat = ptr_patch%cells%center(ilc1,ibc1)%lat
          ELSE
            z_lon = -pi_2 !-90
            z_lat =  pi_2 !+90
          ENDIF
          CALL gvec2cvec(0.0_wp,1.0_wp,z_lon,z_lat,z_cart_no%x(1),z_cart_no%x(2),z_cart_no%x(3))
          z_norm = SQRT( DOT_PRODUCT(z_cart_no%x(1:3),z_cart_no%x(1:3)))
          z_cart_no%x(1:3) = z_cart_no%x(1:3) /z_norm
          CALL gvec2cvec(1.0_wp,0.0_wp,z_lon,z_lat,z_cart_ea%x(1),z_cart_ea%x(2),z_cart_ea%x(3))
          z_norm = SQRT( DOT_PRODUCT(z_cart_ea%x(1:3),z_cart_ea%x(1:3)))
          z_cart_ea%x(1:3) = z_cart_ea%x(1:3) /z_norm
          ! projection of local cell north(east) onto the edge normal
          ptr_int%cno_en(je,1,jb) = &
          &  DOT_PRODUCT(z_cart_no%x(:),ptr_patch%edges%primal_cart_normal(je,jb)%x(:))
          ptr_int%cea_en(je,1,jb) = &
          &  DOT_PRODUCT(z_cart_ea%x(:),ptr_patch%edges%primal_cart_normal(je,jb)%x(:))

          ! neighbor cell 2
          ! cartesian vector of cell center
          IF (.NOT.lplane) THEN
            z_lon = ptr_patch%cells%center(ilc2,ibc2)%lon
            z_lat = ptr_patch%cells%center(ilc2,ibc2)%lat
          ELSE
            z_lon = -pi_2 !-90
            z_lat =  pi_2 !+90
          ENDIF
          CALL gvec2cvec(0.0_wp,1.0_wp,z_lon,z_lat,z_cart_no%x(1),z_cart_no%x(2),z_cart_no%x(3))
          z_norm = SQRT( DOT_PRODUCT(z_cart_no%x(1:3),z_cart_no%x(1:3)))
          z_cart_no%x(1:3) = z_cart_no%x(1:3) /z_norm
          CALL gvec2cvec(1.0_wp,0.0_wp,z_lon,z_lat,z_cart_ea%x(1),z_cart_ea%x(2),z_cart_ea%x(3))
          z_norm = SQRT( DOT_PRODUCT(z_cart_ea%x(1:3),z_cart_ea%x(1:3)))
          z_cart_ea%x(1:3) = z_cart_ea%x(1:3) /z_norm
          ! projection of local cell north(east) onto the edge normal
          ptr_int%cno_en(je,2,jb) = &
          &  DOT_PRODUCT(z_cart_no%x(:),ptr_patch%edges%primal_cart_normal(je,jb)%x(:))
          ptr_int%cea_en(je,2,jb) = &
          &  DOT_PRODUCT(z_cart_ea%x(:),ptr_patch%edges%primal_cart_normal(je,jb)%x(:))

          ! neighbor cell 1
          DO je1 = 1, ptr_patch%cells%num_edges(ilc1,ibc1)

            ! get indices of edges of cell
            ile = ptr_patch%cells%edge_idx(ilc1,ibc1,je1)
            ibe = ptr_patch%cells%edge_blk(ilc1,ibc1,je1)

            ! transform primal normal to cartesian vector z_pn_j
            z_pn_j%x=ptr_patch%edges%primal_cart_normal(ile,ibe)%x

            ! UX coeff (directional Laplace)
            ! get projection of both vectors
            z_proj = DOT_PRODUCT(z_pn_j%x(:),z_pn_k%x(:))
            ! coefficient
            ptr_int%dir_gradhux_c1(je1,je,jb) = z_proj*z_proj &
            & *ptr_int%geofac_div(ilc1,je1,ibc1)

            ! strain deformation coeff
            ! get projection of tangential k and normal j vectors
            z_proj = DOT_PRODUCT(z_pn_j%x(:),z_pt_k%x(:))
            ! coefficient
            ptr_int%strain_def_c1(je1,je,jb) = &
            & ptr_int%dir_gradhux_c1(je1,je,jb) &
            & - z_proj*z_proj*ptr_int%geofac_div(ilc1,je1,ibc1)

            ! indices
            ptr_int%dir_gradh_i1(je1,je,jb) = ile
            ptr_int%dir_gradh_b1(je1,je,jb) = ibe

          ENDDO
          IF (ptr_patch%cells%num_edges(ilc1,ibc1) == 5) THEN
            ptr_int%dir_gradh_i1(6,je,jb) = je
            ptr_int%dir_gradh_b1(6,je,jb) = jb
            ptr_int%dir_gradhux_c1(6,je,jb) = 0._wp
            ptr_int%strain_def_c1(6,je,jb) = 0._wp
          ENDIF

          ! neighbor cell 2
          DO je1 = 1, ptr_patch%cells%num_edges(ilc2,ibc2)

            ! get indices of edges of cell
            ile = ptr_patch%cells%edge_idx(ilc2,ibc2,je1)
            ibe = ptr_patch%cells%edge_blk(ilc2,ibc2,je1)

            ! transform primal normal to cartesian vector z_pn_j
            z_pn_j%x=ptr_patch%edges%primal_cart_normal(ile,ibe)%x

            ! UX coeff (directional Laplace)
            ! get projection of both vectors
            z_proj = DOT_PRODUCT(z_pn_j%x(:),z_pn_k%x(:))
            ! coefficient
            ptr_int%dir_gradhux_c2(je1,je,jb) = z_proj*z_proj &
            & *ptr_int%geofac_div(ilc2,je1,ibc2)

            ! strain deformation coeff
            ! get projection of tangential k and normal j vectors
            z_proj = DOT_PRODUCT(z_pn_j%x(:),z_pt_k%x(:))
            ! coefficient
            ptr_int%strain_def_c2(je1,je,jb) = &
            & ptr_int%dir_gradhux_c2(je1,je,jb) &
            & - z_proj*z_proj*ptr_int%geofac_div(ilc2,je1,ibc2)

            ! indices
            ptr_int%dir_gradh_i2(je1,je,jb) = ile
            ptr_int%dir_gradh_b2(je1,je,jb) = ibe

          ENDDO
          IF (ptr_patch%cells%num_edges(ilc2,ibc2) == 5) THEN
            ptr_int%dir_gradh_i2(6,je,jb) = je
            ptr_int%dir_gradh_b2(6,je,jb) = jb
            ptr_int%dir_gradhux_c2(6,je,jb) = 0._wp
            ptr_int%strain_def_c2(6,je,jb) = 0._wp
          ENDIF

          ! determine neighbor verts
          ilv1 = ptr_patch%edges%vertex_idx(je,jb,1)
          ibv1 = ptr_patch%edges%vertex_blk(je,jb,1)
          ilv2 = ptr_patch%edges%vertex_idx(je,jb,2)
          ibv2 = ptr_patch%edges%vertex_blk(je,jb,2)

          ! neighbor triangle 1 (ilv1,ibv1)
          nincr = 0
          DO jr1 = 1,3

            ! get rhombus indices
            ilr = ptr_patch%verts%edge_idx(ilv1,ibv1,jr1)
            ibr = ptr_patch%verts%edge_blk(ilv1,ibv1,jr1)

            DO je1 = 1, 4
              ! get indices of edges of rhombus
              ile = ptr_patch%edges%quad_idx(ilr,ibr,je1)
              ibe = ptr_patch%edges%quad_blk(ilr,ibr,je1)

              ! transform dual normal to cartesian vector z_pn_j
              z_pn_j%x=ptr_patch%edges%primal_cart_normal(ile,ibe)%x

              ! get projection of both vectors
              z_proj = DOT_PRODUCT(z_pn_j%x(:),z_pn_k%x(:))

              jm = 0
              DO jn=1,nincr
                IF(ile == ptr_int%dir_gradt_i1(jn,je,jb) .AND. &
                   ibe == ptr_int%dir_gradt_b1(jn,je,jb)) THEN
                  ! This was already counted for another quad
                  jm = jn
                  EXIT
                ENDIF
              ENDDO
              IF (jm == 0 ) THEN
                ! This has not yet been counted
                nincr = nincr+1
                jm = nincr
                ptr_int%dir_gradt_i1(jm,je,jb) = ile
                ptr_int%dir_gradt_b1(jm,je,jb) = ibe
              ENDIF
              ! All coefficients have to be initialized with zero in advance!

              ! XY coeff
              ptr_int%dir_gradtxy_v1(jm,je,jb) = ptr_int%dir_gradtxy_v1(jm,je,jb) &
              & - z_proj*z_proj/3.0_wp &
              & /ptr_patch%edges%quad_area(ilr,ibr) &
              & *ptr_patch%edges%dual_edge_length(ile,ibe) &
              & *ptr_patch%edges%quad_orientation(ilr,ibr,je1)*(-1.0_wp)
              ! note here the -1 factor because this is similar to the rotation computation

              ! get projection of both vectors
              z_proj = DOT_PRODUCT(z_pn_j%x(:),z_pt_k%x(:))
              ! YX coeff
              ptr_int%dir_gradtyx_v1(jm,je,jb) = ptr_int%dir_gradtyx_v1(jm,je,jb) &
              & + z_proj*z_proj/3.0_wp &
              & /ptr_patch%edges%quad_area(ilr,ibr) &
              & *ptr_patch%edges%dual_edge_length(ile,ibe) &
              & *ptr_patch%edges%quad_orientation(ilr,ibr,je1)*(-1.0_wp)
              ! note here the -1 factor because this is similar to the rotation computation

            ENDDO
          ENDDO
          DO jm = 1, 9
            ptr_int%shear_def_v1(jm,je,jb) = &
            & ptr_int%dir_gradtxy_v1(jm,je,jb)+ptr_int%dir_gradtyx_v1(jm,je,jb)
          ENDDO

          ! neighbor triangle 2 (ilv2,ibv2)
          nincr = 0
          DO jr1 = 1,3

            ! get rhombus indices
            ilr = ptr_patch%verts%edge_idx(ilv2,ibv2,jr1)
            ibr = ptr_patch%verts%edge_blk(ilv2,ibv2,jr1)

            DO je1 = 1, 4
              ! get indices of edges of rhombus
              ile = ptr_patch%edges%quad_idx(ilr,ibr,je1)
              ibe = ptr_patch%edges%quad_blk(ilr,ibr,je1)

              ! transform dual normal to cartesian vector z_pn_j
              z_pn_j%x=ptr_patch%edges%primal_cart_normal(ile,ibe)%x

              ! get projection of both vectors
              z_proj = DOT_PRODUCT(z_pn_j%x(:),z_pn_k%x(:))

              jm = 0
              DO jn=1,nincr
                IF(ile == ptr_int%dir_gradt_i2(jn,je,jb) .AND. &
                   ibe == ptr_int%dir_gradt_b2(jn,je,jb)) THEN
                  ! This was already counted for another quad
                  jm = jn
                  EXIT
                ENDIF
              ENDDO
              IF (jm == 0 ) THEN
                ! This has not yet been counted
                nincr = nincr+1
                jm = nincr
                ptr_int%dir_gradt_i2(jm,je,jb) = ile
                ptr_int%dir_gradt_b2(jm,je,jb) = ibe
              ENDIF
              ! All coefficients have to be initialized with zero in advance!

              ! XY coeff
              ptr_int%dir_gradtxy_v2(jm,je,jb) = ptr_int%dir_gradtxy_v2(jm,je,jb) &
              & -z_proj*z_proj/3.0_wp &
              & /ptr_patch%edges%quad_area(ilr,ibr) &
              & *ptr_patch%edges%dual_edge_length(ile,ibe) &
              & *ptr_patch%edges%quad_orientation(ilr,ibr,je1)*(-1.0_wp)

              ! get projection of both vectors
              z_proj = DOT_PRODUCT(z_pn_j%x(:),z_pt_k%x(:))
              ptr_int%dir_gradtyx_v2(jm,je,jb) = ptr_int%dir_gradtyx_v2(jm,je,jb) &
              & +z_proj*z_proj/3.0_wp &
              & /ptr_patch%edges%quad_area(ilr,ibr) &
              & *ptr_patch%edges%dual_edge_length(ile,ibe) &
              & *ptr_patch%edges%quad_orientation(ilr,ibr,je1)*(-1.0_wp)

            ENDDO !je1
          ENDDO !jr1
          DO jm = 1, 9
            ptr_int%shear_def_v2(jm,je,jb) = &
            & ptr_int%dir_gradtxy_v2(jm,je,jb)+ptr_int%dir_gradtyx_v2(jm,je,jb)
          ENDDO

        ENDDO ! je

      ENDDO !jb
!$OMP END DO

    ENDIF

! f) Geometrical factor for Green-Gauss gradient
rl_start = 2
rl_end = min_rlcell

! values for the blocking
i_startblk = ptr_patch%cells%start_blk(rl_start,1)
i_endblk   = ptr_patch%cells%end_blk(rl_end,i_nchdom)
!
! loop through all patch cells (and blocks)
!
!$OMP DO PRIVATE(jb,je,jc,ic,i_startidx,i_endidx,ile,ibe,ilc1,ibc1,&
!$OMP    ilc2,ibc2,ilnc,ibnc)
DO jb = i_startblk, i_endblk

  CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
                     i_startidx, i_endidx, rl_start, rl_end)

  DO je = 1, ptr_patch%cell_type
    DO jc = i_startidx, i_endidx

      IF(.NOT. ptr_patch%cells%owner_mask(jc,jb)) CYCLE

      ile = ptr_patch%cells%edge_idx(jc,jb,je)
      ibe = ptr_patch%cells%edge_blk(jc,jb,je)

      ilc1 = ptr_patch%edges%cell_idx(ile,ibe,1)
      ibc1 = ptr_patch%edges%cell_blk(ile,ibe,1)
      ilc2 = ptr_patch%edges%cell_idx(ile,ibe,2)
      ibc2 = ptr_patch%edges%cell_blk(ile,ibe,2)

      IF (jc == ilc1 .AND. jb == ibc1) THEN
        ptr_int%geofac_grg(jc,1,jb,1) = ptr_int%geofac_grg(jc,1,jb,1) +      &
          ptr_int%primal_normal_ec(jc,jb,je,1)*ptr_int%geofac_div(jc,je,jb)* &
          ptr_int%c_lin_e(ile,1,ibe)
        ptr_int%geofac_grg(jc,1,jb,2) = ptr_int%geofac_grg(jc,1,jb,2) +      &
          ptr_int%primal_normal_ec(jc,jb,je,2)*ptr_int%geofac_div(jc,je,jb)* &
          ptr_int%c_lin_e(ile,1,ibe)
      ELSE IF (jc == ilc2 .AND. jb == ibc2) THEN
        ptr_int%geofac_grg(jc,1,jb,1) = ptr_int%geofac_grg(jc,1,jb,1) +      &
          ptr_int%primal_normal_ec(jc,jb,je,1)*ptr_int%geofac_div(jc,je,jb)* &
          ptr_int%c_lin_e(ile,2,ibe)
        ptr_int%geofac_grg(jc,1,jb,2) = ptr_int%geofac_grg(jc,1,jb,2) +      &
          ptr_int%primal_normal_ec(jc,jb,je,2)*ptr_int%geofac_div(jc,je,jb)* &
          ptr_int%c_lin_e(ile,2,ibe)
      ENDIF
      DO ic = 1, ptr_patch%cell_type
        ilnc = ptr_patch%cells%neighbor_idx(jc,jb,ic)
        ibnc = ptr_patch%cells%neighbor_blk(jc,jb,ic)
        IF (ilnc == ilc1 .AND. ibnc == ibc1) THEN
          ptr_int%geofac_grg(jc,ic+1,jb,1) = ptr_int%geofac_grg(jc,ic+1,jb,1)+ &
            ptr_int%primal_normal_ec(jc,jb,je,1)*ptr_int%geofac_div(jc,je,jb)* &
            ptr_int%c_lin_e(ile,1,ibe)
          ptr_int%geofac_grg(jc,ic+1,jb,2) = ptr_int%geofac_grg(jc,ic+1,jb,2)+ &
            ptr_int%primal_normal_ec(jc,jb,je,2)*ptr_int%geofac_div(jc,je,jb)* &
            ptr_int%c_lin_e(ile,1,ibe)
        ELSE IF (ilnc == ilc2 .AND. ibnc == ibc2) THEN
          ptr_int%geofac_grg(jc,ic+1,jb,1) = ptr_int%geofac_grg(jc,ic+1,jb,1)+ &
            ptr_int%primal_normal_ec(jc,jb,je,1)*ptr_int%geofac_div(jc,je,jb)* &
            ptr_int%c_lin_e(ile,2,ibe)
          ptr_int%geofac_grg(jc,ic+1,jb,2) = ptr_int%geofac_grg(jc,ic+1,jb,2)+ &
            ptr_int%primal_normal_ec(jc,jb,je,2)*ptr_int%geofac_div(jc,je,jb)* &
            ptr_int%c_lin_e(ile,2,ibe)
        ENDIF
      ENDDO

      ! To ensure that dummy edges have a factor of 0:
      IF (je > ptr_patch%cells%num_edges(jc,jb)) THEN
        ptr_int%geofac_grg(jc,je+1,jb,1:2) = 0._wp
      ENDIF

    ENDDO !cell loop
  ENDDO

END DO !block loop
!$OMP END DO

!$OMP END PARALLEL

CALL sync_patch_array(SYNC_C,ptr_patch,ptr_int%geofac_div)
CALL sync_patch_array(SYNC_V,ptr_patch,ptr_int%geofac_rot)
CALL sync_patch_array(SYNC_C,ptr_patch,ptr_int%geofac_n2s)

IF (ptr_patch%cell_type == 3) THEN
  CALL sync_patch_array(SYNC_E,ptr_patch,ptr_int%geofac_qdiv)
ENDIF

IF (ptr_patch%cell_type == 6) THEN
  CALL sync_patch_array(SYNC_E,ptr_patch,ptr_int%cno_en)
  CALL sync_patch_array(SYNC_E,ptr_patch,ptr_int%cea_en)

  DO jm = 1, 6
    CALL sync_idx(SYNC_E,SYNC_E,ptr_patch,ptr_int%dir_gradh_i1(jm,:,:), &
                                        & ptr_int%dir_gradh_b1(jm,:,:))
    CALL sync_idx(SYNC_E,SYNC_E,ptr_patch,ptr_int%dir_gradh_i2(jm,:,:), &
                                        & ptr_int%dir_gradh_b2(jm,:,:))
    CALL sync_patch_array(SYNC_E,ptr_patch,ptr_int%dir_gradhux_c1(jm, :, :))
    CALL sync_patch_array(SYNC_E,ptr_patch,ptr_int%dir_gradhux_c2(jm, :, :))
    CALL sync_patch_array(SYNC_E,ptr_patch,ptr_int%strain_def_c1(jm, :, :))
    CALL sync_patch_array(SYNC_E,ptr_patch,ptr_int%strain_def_c2(jm, :, :))
  ENDDO

  DO jm = 1, 9
    CALL sync_idx(SYNC_E,SYNC_E,ptr_patch,ptr_int%dir_gradt_i1(jm,:,:), &
                                        & ptr_int%dir_gradt_b1(jm,:,:))
    CALL sync_idx(SYNC_E,SYNC_E,ptr_patch,ptr_int%dir_gradt_i2(jm,:,:), &
                                        & ptr_int%dir_gradt_b2(jm,:,:))
    CALL sync_patch_array(SYNC_E,ptr_patch,ptr_int%dir_gradtxy_v1(jm, :, :))
    CALL sync_patch_array(SYNC_E,ptr_patch,ptr_int%dir_gradtxy_v2(jm, :, :))
    CALL sync_patch_array(SYNC_E,ptr_patch,ptr_int%dir_gradtyx_v1(jm, :, :))
    CALL sync_patch_array(SYNC_E,ptr_patch,ptr_int%dir_gradtyx_v2(jm, :, :))
    CALL sync_patch_array(SYNC_E,ptr_patch,ptr_int%shear_def_v1(jm, :, :))
    CALL sync_patch_array(SYNC_E,ptr_patch,ptr_int%shear_def_v2(jm, :, :))
  ENDDO
ENDIF

CALL sync_patch_array(SYNC_C,ptr_patch,ptr_int%geofac_grg(:,:,:,1))
CALL sync_patch_array(SYNC_C,ptr_patch,ptr_int%geofac_grg(:,:,:,2))



END SUBROUTINE init_geo_factors



!-------------------------------------------------------------------------
!
!
!>
!! Computes the local orientation of the edge primal normal and dual normal.
!!
!! Computes the local orientation of the edge primal normal and dual normal
!! at the location of the cell centers and vertices.
!! Moreover, the Cartesian orientation vectors of the edge primal normals
!! are stored for use in the RBF initialization routines, and inverse
!! primal and dual edge lengths are computed
!!
!! @par Revision History
!!  developed by Guenther Zaengl, 2009-03-31
!!
SUBROUTINE complete_patchinfo( ptr_patch, ptr_int )
!

!
!  patch on which computation is performed
!
TYPE(t_patch), TARGET, INTENT(inout) :: ptr_patch

! Interpolation state
TYPE(t_int_state),     INTENT(inout) :: ptr_int
!

INTEGER :: jb, je, jc
INTEGER :: rl_start, rl_end
INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx, i_nchdom

INTEGER :: ilc1, ibc1, ilv1, ibv1, ilc2, ibc2, ilv2, ibv2, &
           ilv3, ibv3, ilv4, ibv4, ile1, ibe1

REAL(wp) :: z_nu, z_nv, z_lon, z_lat, z_nx1(3), z_nx2(3), z_norm

TYPE(t_cartesian_coordinates) :: cc_edge, cc_cell, cc_ev3, cc_ev4

!-----------------------------------------------------------------------

i_nchdom   = MAX(1,ptr_patch%n_childdom)

!$OMP PARALLEL  PRIVATE(rl_start,rl_end,i_startblk,i_endblk)
rl_start = 1
rl_end = min_rlcell

! values for the blocking
i_startblk = ptr_patch%cells%start_blk(rl_start,1)
i_endblk   = ptr_patch%cells%end_blk(rl_end,i_nchdom)
!
! First step: compute Cartesian coordinates of cell centers on full domain
! this is needed for least squares gradient reconstruction;
!
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,cc_cell)
DO jb = i_startblk, i_endblk

  CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
                     i_startidx, i_endidx, rl_start, rl_end)

  DO jc =  i_startidx, i_endidx

    IF(.NOT.ptr_patch%cells%owner_mask(jc,jb)) CYCLE

    ! compute Cartesian coordinates (needed for RBF initialization)

    cc_cell = gc2cc(ptr_patch%cells%center(jc,jb))
    ptr_int%cart_cell_coord(jc,jb,1:3)= cc_cell%x(1:3)

  ENDDO

END DO !block loop
!$OMP END DO


rl_start = 1
rl_end = min_rledge

! values for the blocking
i_startblk = ptr_patch%edges%start_blk(rl_start,1)
i_endblk   = ptr_patch%edges%end_blk(rl_end,i_nchdom)
!
! First step: compute Cartesian coordinates and Cartesian vectors on full domain
! this is needed to vectorize RBF initialization; the existing field carrying
! the Cartesian orientation vectors (primal_cart_normal) did not work for that
! because it is a derived data type
! In addition, the fields for the inverse primal and dual edge lengths are
! initialized here.
!
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,je,cc_edge)
DO jb = i_startblk, i_endblk

  CALL get_indices_e(ptr_patch, jb, i_startblk, i_endblk, &
                     i_startidx, i_endidx, rl_start, rl_end)

  DO je =  i_startidx, i_endidx

    IF(.NOT.ptr_patch%edges%owner_mask(je,jb)) CYCLE

    ! compute Cartesian coordinates (needed for RBF initialization)

    cc_edge = gc2cc(ptr_patch%edges%center(je,jb))
    ptr_int%cart_edge_coord(je,jb,1:3)= cc_edge%x(1:3)

    ! finally, compute inverse primal edge length
    ! (dual follows below in the rl_start=2 section)

    ptr_patch%edges%inv_primal_edge_length(je,jb) = &
      1._wp/ptr_patch%edges%primal_edge_length(je,jb)

  ENDDO

END DO !block loop
!$OMP END DO

rl_start = 2
rl_end = min_rledge

! Second step: computed projected orientation vectors and related information
i_startblk = ptr_patch%edges%start_blk(rl_start,1)
i_endblk   = ptr_patch%edges%end_blk(rl_end,i_nchdom)

! Initialization of lateral boundary points 
IF (ptr_patch%id > 1) THEN
!$OMP WORKSHARE
  ptr_patch%edges%inv_dual_edge_length(:,1:i_startblk)    = 0._wp
  ptr_patch%edges%vertex_idx(:,1:i_startblk,3)            = 0
  ptr_patch%edges%vertex_idx(:,1:i_startblk,4)            = 0
  ptr_patch%edges%vertex_blk(:,1:i_startblk,3)            = 0
  ptr_patch%edges%vertex_blk(:,1:i_startblk,4)            = 0
  ptr_patch%edges%inv_vert_vert_length(:,1:i_startblk)    = 0._wp
  ptr_patch%edges%primal_normal_cell(:,1:i_startblk,:)%v1 = 0._wp
  ptr_patch%edges%dual_normal_cell  (:,1:i_startblk,:)%v1 = 0._wp
  ptr_patch%edges%primal_normal_vert(:,1:i_startblk,:)%v1 = 0._wp
  ptr_patch%edges%dual_normal_vert  (:,1:i_startblk,:)%v1 = 0._wp
  ptr_patch%edges%primal_normal_cell(:,1:i_startblk,:)%v2 = 0._wp
  ptr_patch%edges%dual_normal_cell  (:,1:i_startblk,:)%v2 = 0._wp
  ptr_patch%edges%primal_normal_vert(:,1:i_startblk,:)%v2 = 0._wp
  ptr_patch%edges%dual_normal_vert  (:,1:i_startblk,:)%v2 = 0._wp
!$OMP END WORKSHARE
ENDIF
!
! loop through all patch edges
!
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,je,ilc1,ibc1,ilv1,ibv1,ilc2,ibc2,ilv2, &
!$OMP            ibv2,ilv3,ibv3,ilv4,ibv4,z_nu,z_nv,z_lon,z_lat,z_nx1,z_nx2,   &
!$OMP            cc_ev3,cc_ev4,z_norm)
DO jb = i_startblk, i_endblk

  CALL get_indices_e(ptr_patch, jb, i_startblk, i_endblk, &
                     i_startidx, i_endidx, rl_start, rl_end)

  DO je =  i_startidx, i_endidx

    IF(.NOT.ptr_patch%edges%owner_mask(je,jb)) CYCLE

    ! compute inverse dual edge length (undefined for refin_ctrl=1)

    ptr_patch%edges%inv_dual_edge_length(je,jb) = &
      1._wp/ptr_patch%edges%dual_edge_length(je,jb)

    ! compute edge-vertex indices (and blocks) 3 and 4, which
    ! are the outer vertices of cells 1 and 2, respectively,
    ! and the inverse length bewtween vertices 3 and 4

    ilc1 = ptr_patch%edges%cell_idx(je,jb,1)
    ibc1 = ptr_patch%edges%cell_blk(je,jb,1)
    ilc2 = ptr_patch%edges%cell_idx(je,jb,2)
    ibc2 = ptr_patch%edges%cell_blk(je,jb,2)

    ilv1 = ptr_patch%edges%vertex_idx(je,jb,1)
    ibv1 = ptr_patch%edges%vertex_blk(je,jb,1)
    ilv2 = ptr_patch%edges%vertex_idx(je,jb,2)
    ibv2 = ptr_patch%edges%vertex_blk(je,jb,2)

    IF ((ptr_patch%cells%vertex_idx(ilc1,ibc1,1) /= &
         ptr_patch%edges%vertex_idx(je,jb,1) .OR.  &
         ptr_patch%cells%vertex_blk(ilc1,ibc1,1) /= &
         ptr_patch%edges%vertex_blk(je,jb,1)) .AND.  &
        (ptr_patch%cells%vertex_idx(ilc1,ibc1,1) /= &
         ptr_patch%edges%vertex_idx(je,jb,2) .OR.  &
         ptr_patch%cells%vertex_blk(ilc1,ibc1,1) /= &
         ptr_patch%edges%vertex_blk(je,jb,2)) )        THEN

      ptr_patch%edges%vertex_idx(je,jb,3) = ptr_patch%cells%vertex_idx(ilc1,ibc1,1)
      ptr_patch%edges%vertex_blk(je,jb,3) = ptr_patch%cells%vertex_blk(ilc1,ibc1,1)

    ELSE IF ((ptr_patch%cells%vertex_idx(ilc1,ibc1,2) /= &
              ptr_patch%edges%vertex_idx(je,jb,1) .OR.  &
              ptr_patch%cells%vertex_blk(ilc1,ibc1,2) /= &
              ptr_patch%edges%vertex_blk(je,jb,1)) .AND.  &
             (ptr_patch%cells%vertex_idx(ilc1,ibc1,2) /= &
              ptr_patch%edges%vertex_idx(je,jb,2) .OR.  &
              ptr_patch%cells%vertex_blk(ilc1,ibc1,2) /= &
              ptr_patch%edges%vertex_blk(je,jb,2)) )        THEN

      ptr_patch%edges%vertex_idx(je,jb,3) = ptr_patch%cells%vertex_idx(ilc1,ibc1,2)
      ptr_patch%edges%vertex_blk(je,jb,3) = ptr_patch%cells%vertex_blk(ilc1,ibc1,2)

    ELSE IF ((ptr_patch%cells%vertex_idx(ilc1,ibc1,3) /= &
              ptr_patch%edges%vertex_idx(je,jb,1) .OR.  &
              ptr_patch%cells%vertex_blk(ilc1,ibc1,3) /= &
              ptr_patch%edges%vertex_blk(je,jb,1)) .AND.  &
             (ptr_patch%cells%vertex_idx(ilc1,ibc1,3) /= &
              ptr_patch%edges%vertex_idx(je,jb,2) .OR.  &
              ptr_patch%cells%vertex_blk(ilc1,ibc1,3) /= &
              ptr_patch%edges%vertex_blk(je,jb,2)) )        THEN

      ptr_patch%edges%vertex_idx(je,jb,3) = ptr_patch%cells%vertex_idx(ilc1,ibc1,3)
      ptr_patch%edges%vertex_blk(je,jb,3) = ptr_patch%cells%vertex_blk(ilc1,ibc1,3)

    ENDIF

    IF ((ptr_patch%cells%vertex_idx(ilc2,ibc2,1) /= &
         ptr_patch%edges%vertex_idx(je,jb,1) .OR.  &
         ptr_patch%cells%vertex_blk(ilc2,ibc2,1) /= &
         ptr_patch%edges%vertex_blk(je,jb,1)) .AND.  &
        (ptr_patch%cells%vertex_idx(ilc2,ibc2,1) /= &
         ptr_patch%edges%vertex_idx(je,jb,2) .OR.  &
         ptr_patch%cells%vertex_blk(ilc2,ibc2,1) /= &
         ptr_patch%edges%vertex_blk(je,jb,2)) )        THEN

      ptr_patch%edges%vertex_idx(je,jb,4) = ptr_patch%cells%vertex_idx(ilc2,ibc2,1)
      ptr_patch%edges%vertex_blk(je,jb,4) = ptr_patch%cells%vertex_blk(ilc2,ibc2,1)

    ELSE IF ((ptr_patch%cells%vertex_idx(ilc2,ibc2,2) /= &
              ptr_patch%edges%vertex_idx(je,jb,1) .OR.  &
              ptr_patch%cells%vertex_blk(ilc2,ibc2,2) /= &
              ptr_patch%edges%vertex_blk(je,jb,1)) .AND.  &
             (ptr_patch%cells%vertex_idx(ilc2,ibc2,2) /= &
              ptr_patch%edges%vertex_idx(je,jb,2) .OR.  &
              ptr_patch%cells%vertex_blk(ilc2,ibc2,2) /= &
              ptr_patch%edges%vertex_blk(je,jb,2)) )        THEN

      ptr_patch%edges%vertex_idx(je,jb,4) = ptr_patch%cells%vertex_idx(ilc2,ibc2,2)
      ptr_patch%edges%vertex_blk(je,jb,4) = ptr_patch%cells%vertex_blk(ilc2,ibc2,2)

    ELSE IF ((ptr_patch%cells%vertex_idx(ilc2,ibc2,3) /= &
              ptr_patch%edges%vertex_idx(je,jb,1) .OR.  &
              ptr_patch%cells%vertex_blk(ilc2,ibc2,3) /= &
              ptr_patch%edges%vertex_blk(je,jb,1)) .AND.  &
             (ptr_patch%cells%vertex_idx(ilc2,ibc2,3) /= &
              ptr_patch%edges%vertex_idx(je,jb,2) .OR.  &
              ptr_patch%cells%vertex_blk(ilc2,ibc2,3) /= &
              ptr_patch%edges%vertex_blk(je,jb,2)) )        THEN

      ptr_patch%edges%vertex_idx(je,jb,4) = ptr_patch%cells%vertex_idx(ilc2,ibc2,3)
      ptr_patch%edges%vertex_blk(je,jb,4) = ptr_patch%cells%vertex_blk(ilc2,ibc2,3)

    ENDIF

    ilv3 = ptr_patch%edges%vertex_idx(je,jb,3)
    ibv3 = ptr_patch%edges%vertex_blk(je,jb,3)
    ilv4 = ptr_patch%edges%vertex_idx(je,jb,4)
    ibv4 = ptr_patch%edges%vertex_blk(je,jb,4)

    cc_ev3 = gc2cc(ptr_patch%verts%vertex(ilv3,ibv3))
    cc_ev4 = gc2cc(ptr_patch%verts%vertex(ilv4,ibv4))

    ! inverse length bewtween vertices 3 and 4
    IF (ptr_patch%cell_type == 3 ) THEN
      ptr_patch%edges%inv_vert_vert_length(je,jb) = 1._wp/(re*arc_length(cc_ev3,cc_ev4))
    ENDIF

    ! next step: compute projected orientation vectors for cells and vertices
    ! bordering to each edge (incl. vertices 3 and 4 intorduced above)

    ! transform orientation vectors at local edge center to Cartesian space
    z_lon = ptr_patch%edges%center(je,jb)%lon
    z_lat = ptr_patch%edges%center(je,jb)%lat

    ! transform primal normal to cartesian vector z_nx1
    z_nx1(:)=ptr_patch%edges%primal_cart_normal(je,jb)%x(:)

    ! transform dual normal to cartesian vector z_nx2
    z_nx2(:)=ptr_patch%edges%dual_cart_normal(je,jb)%x(:)

    ! get location of cell 1

    z_lon = ptr_patch%cells%center(ilc1,ibc1)%lon
    z_lat = ptr_patch%cells%center(ilc1,ibc1)%lat

    ! compute local primal and dual normals at cell 1

    CALL cvec2gvec(z_nx1(1),z_nx1(2),z_nx1(3),z_lon,z_lat,z_nu,z_nv)
    z_norm = SQRT(z_nu*z_nu+z_nv*z_nv)

    ptr_patch%edges%primal_normal_cell(je,jb,1)%v1 = z_nu/z_norm
    ptr_patch%edges%primal_normal_cell(je,jb,1)%v2 = z_nv/z_norm

    CALL cvec2gvec(z_nx2(1),z_nx2(2),z_nx2(3),z_lon,z_lat,z_nu,z_nv)
    z_norm = SQRT(z_nu*z_nu+z_nv*z_nv)

    ptr_patch%edges%dual_normal_cell(je,jb,1)%v1 = z_nu/z_norm
    ptr_patch%edges%dual_normal_cell(je,jb,1)%v2 = z_nv/z_norm

    ! get location of cell 2

    z_lon = ptr_patch%cells%center(ilc2,ibc2)%lon
    z_lat = ptr_patch%cells%center(ilc2,ibc2)%lat

    ! compute local primal and dual normals at cell 2

    CALL cvec2gvec(z_nx1(1),z_nx1(2),z_nx1(3),z_lon,z_lat,z_nu,z_nv)
    z_norm = SQRT(z_nu*z_nu+z_nv*z_nv)

    ptr_patch%edges%primal_normal_cell(je,jb,2)%v1 = z_nu/z_norm
    ptr_patch%edges%primal_normal_cell(je,jb,2)%v2 = z_nv/z_norm

    CALL cvec2gvec(z_nx2(1),z_nx2(2),z_nx2(3),z_lon,z_lat,z_nu,z_nv)
    z_norm = SQRT(z_nu*z_nu+z_nv*z_nv)

    ptr_patch%edges%dual_normal_cell(je,jb,2)%v1 = z_nu/z_norm
    ptr_patch%edges%dual_normal_cell(je,jb,2)%v2 = z_nv/z_norm

    ! get location of vertex 1

    z_lon = ptr_patch%verts%vertex(ilv1,ibv1)%lon
    z_lat = ptr_patch%verts%vertex(ilv1,ibv1)%lat

    ! compute local primal and dual normals at vertex 1

    CALL cvec2gvec(z_nx1(1),z_nx1(2),z_nx1(3),z_lon,z_lat,z_nu,z_nv)
    z_norm = SQRT(z_nu*z_nu+z_nv*z_nv)

    ptr_patch%edges%primal_normal_vert(je,jb,1)%v1 = z_nu/z_norm
    ptr_patch%edges%primal_normal_vert(je,jb,1)%v2 = z_nv/z_norm

    CALL cvec2gvec(z_nx2(1),z_nx2(2),z_nx2(3),z_lon,z_lat,z_nu,z_nv)
    z_norm = SQRT(z_nu*z_nu+z_nv*z_nv)

    ptr_patch%edges%dual_normal_vert(je,jb,1)%v1 = z_nu/z_norm
    ptr_patch%edges%dual_normal_vert(je,jb,1)%v2 = z_nv/z_norm

    ! get location of vertex 2

    z_lon = ptr_patch%verts%vertex(ilv2,ibv2)%lon
    z_lat = ptr_patch%verts%vertex(ilv2,ibv2)%lat

    ! compute local primal and dual normals at vertex 2

    CALL cvec2gvec(z_nx1(1),z_nx1(2),z_nx1(3),z_lon,z_lat,z_nu,z_nv)
    z_norm = SQRT(z_nu*z_nu+z_nv*z_nv)

    ptr_patch%edges%primal_normal_vert(je,jb,2)%v1 = z_nu/z_norm
    ptr_patch%edges%primal_normal_vert(je,jb,2)%v2 = z_nv/z_norm

    CALL cvec2gvec(z_nx2(1),z_nx2(2),z_nx2(3),z_lon,z_lat,z_nu,z_nv)
    z_norm = SQRT(z_nu*z_nu+z_nv*z_nv)

    ptr_patch%edges%dual_normal_vert(je,jb,2)%v1 = z_nu/z_norm
    ptr_patch%edges%dual_normal_vert(je,jb,2)%v2 = z_nv/z_norm

    ! get location of vertex 3

    z_lon = ptr_patch%verts%vertex(ilv3,ibv3)%lon
    z_lat = ptr_patch%verts%vertex(ilv3,ibv3)%lat

    ! compute local primal and dual normals at vertex 3

    CALL cvec2gvec(z_nx1(1),z_nx1(2),z_nx1(3),z_lon,z_lat,z_nu,z_nv)
    z_norm = SQRT(z_nu*z_nu+z_nv*z_nv)

    ptr_patch%edges%primal_normal_vert(je,jb,3)%v1 = z_nu/z_norm
    ptr_patch%edges%primal_normal_vert(je,jb,3)%v2 = z_nv/z_norm

    CALL cvec2gvec(z_nx2(1),z_nx2(2),z_nx2(3),z_lon,z_lat,z_nu,z_nv)
    z_norm = SQRT(z_nu*z_nu+z_nv*z_nv)

    ptr_patch%edges%dual_normal_vert(je,jb,3)%v1 = z_nu/z_norm
    ptr_patch%edges%dual_normal_vert(je,jb,3)%v2 = z_nv/z_norm

    ! get location of vertex 4

    z_lon = ptr_patch%verts%vertex(ilv4,ibv4)%lon
    z_lat = ptr_patch%verts%vertex(ilv4,ibv4)%lat

    ! compute local primal and dual normals at vertex 2

    CALL cvec2gvec(z_nx1(1),z_nx1(2),z_nx1(3),z_lon,z_lat,z_nu,z_nv)
    z_norm = SQRT(z_nu*z_nu+z_nv*z_nv)

    ptr_patch%edges%primal_normal_vert(je,jb,4)%v1 = z_nu/z_norm
    ptr_patch%edges%primal_normal_vert(je,jb,4)%v2 = z_nv/z_norm

    CALL cvec2gvec(z_nx2(1),z_nx2(2),z_nx2(3),z_lon,z_lat,z_nu,z_nv)
    z_norm = SQRT(z_nu*z_nu+z_nv*z_nv)

    ptr_patch%edges%dual_normal_vert(je,jb,4)%v1 = z_nu/z_norm
    ptr_patch%edges%dual_normal_vert(je,jb,4)%v2 = z_nv/z_norm

  ENDDO

END DO !block loop
!$OMP END DO

!$OMP END PARALLEL

  ! primal_normal_cell must be sync'd before next loop,
  ! so do a sync for all above calculated quantities

  CALL sync_patch_array(SYNC_C,ptr_patch,ptr_int%cart_cell_coord(:,:,1))
  CALL sync_patch_array(SYNC_C,ptr_patch,ptr_int%cart_cell_coord(:,:,2))
  CALL sync_patch_array(SYNC_C,ptr_patch,ptr_int%cart_cell_coord(:,:,3))
  CALL sync_patch_array(SYNC_E,ptr_patch,ptr_int%cart_edge_coord(:,:,1))
  CALL sync_patch_array(SYNC_E,ptr_patch,ptr_int%cart_edge_coord(:,:,2))
  CALL sync_patch_array(SYNC_E,ptr_patch,ptr_int%cart_edge_coord(:,:,3))
  CALL sync_patch_array(SYNC_E,ptr_patch,ptr_patch%edges%inv_primal_edge_length)

  CALL sync_idx(SYNC_E,SYNC_V,ptr_patch,ptr_patch%edges%vertex_idx(:,:,3), &
                                      & ptr_patch%edges%vertex_blk(:,:,3))
  CALL sync_idx(SYNC_E,SYNC_V,ptr_patch,ptr_patch%edges%vertex_idx(:,:,4), &
                                      & ptr_patch%edges%vertex_blk(:,:,4))

  CALL sync_patch_array(SYNC_E,ptr_patch,ptr_patch%edges%inv_dual_edge_length)
  CALL sync_patch_array(SYNC_E,ptr_patch,ptr_patch%edges%inv_vert_vert_length)

  CALL sync_patch_array(SYNC_E,ptr_patch,ptr_patch%edges%primal_normal_cell(:,:,1)%v1)
  CALL sync_patch_array(SYNC_E,ptr_patch,ptr_patch%edges%primal_normal_cell(:,:,2)%v1)
  CALL sync_patch_array(SYNC_E,ptr_patch,ptr_patch%edges%primal_normal_vert(:,:,1)%v1)
  CALL sync_patch_array(SYNC_E,ptr_patch,ptr_patch%edges%primal_normal_vert(:,:,2)%v1)
  CALL sync_patch_array(SYNC_E,ptr_patch,ptr_patch%edges%primal_normal_vert(:,:,3)%v1)
  CALL sync_patch_array(SYNC_E,ptr_patch,ptr_patch%edges%primal_normal_vert(:,:,4)%v1)

  CALL sync_patch_array(SYNC_E,ptr_patch,ptr_patch%edges%dual_normal_cell(:,:,1)%v1)
  CALL sync_patch_array(SYNC_E,ptr_patch,ptr_patch%edges%dual_normal_cell(:,:,2)%v1)
  CALL sync_patch_array(SYNC_E,ptr_patch,ptr_patch%edges%dual_normal_vert(:,:,1)%v1)
  CALL sync_patch_array(SYNC_E,ptr_patch,ptr_patch%edges%dual_normal_vert(:,:,2)%v1)
  CALL sync_patch_array(SYNC_E,ptr_patch,ptr_patch%edges%dual_normal_vert(:,:,3)%v1)
  CALL sync_patch_array(SYNC_E,ptr_patch,ptr_patch%edges%dual_normal_vert(:,:,4)%v1)

  CALL sync_patch_array(SYNC_E,ptr_patch,ptr_patch%edges%primal_normal_cell(:,:,1)%v2)
  CALL sync_patch_array(SYNC_E,ptr_patch,ptr_patch%edges%primal_normal_cell(:,:,2)%v2)
  CALL sync_patch_array(SYNC_E,ptr_patch,ptr_patch%edges%primal_normal_vert(:,:,1)%v2)
  CALL sync_patch_array(SYNC_E,ptr_patch,ptr_patch%edges%primal_normal_vert(:,:,2)%v2)
  CALL sync_patch_array(SYNC_E,ptr_patch,ptr_patch%edges%primal_normal_vert(:,:,3)%v2)
  CALL sync_patch_array(SYNC_E,ptr_patch,ptr_patch%edges%primal_normal_vert(:,:,4)%v2)

  CALL sync_patch_array(SYNC_E,ptr_patch,ptr_patch%edges%dual_normal_cell(:,:,1)%v2)
  CALL sync_patch_array(SYNC_E,ptr_patch,ptr_patch%edges%dual_normal_cell(:,:,2)%v2)
  CALL sync_patch_array(SYNC_E,ptr_patch,ptr_patch%edges%dual_normal_vert(:,:,1)%v2)
  CALL sync_patch_array(SYNC_E,ptr_patch,ptr_patch%edges%dual_normal_vert(:,:,2)%v2)
  CALL sync_patch_array(SYNC_E,ptr_patch,ptr_patch%edges%dual_normal_vert(:,:,3)%v2)
  CALL sync_patch_array(SYNC_E,ptr_patch,ptr_patch%edges%dual_normal_vert(:,:,4)%v2)


!$OMP PARALLEL  PRIVATE(rl_start,rl_end,i_startblk,i_endblk)

rl_start = 2
rl_end = min_rlcell

! Final step: store primal_normal_cell also with respect to cell points
! in order to reduce indirect addressing during runtime
i_startblk = ptr_patch%cells%start_blk(rl_start,1)
i_endblk   = ptr_patch%cells%end_blk(rl_end,i_nchdom)
!
! loop through all patch cells
!
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,ile1,ibe1)
DO jb = i_startblk, i_endblk

  CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
                     i_startidx, i_endidx, rl_start, rl_end)

  DO jc =  i_startidx, i_endidx

    IF(.NOT.ptr_patch%cells%owner_mask(jc,jb)) CYCLE

    DO je = 1, ptr_patch%cells%num_edges(jc,jb)

      ile1 = ptr_patch%cells%edge_idx(jc,jb,je)
      ibe1 = ptr_patch%cells%edge_blk(jc,jb,je)


      IF ((ptr_patch%edges%cell_idx(ile1,ibe1,1) == jc) .AND. &
          (ptr_patch%edges%cell_blk(ile1,ibe1,1) == jb)) THEN

        ptr_int%primal_normal_ec(jc,jb,je,1) = &
          ptr_patch%edges%primal_normal_cell(ile1,ibe1,1)%v1
        ptr_int%primal_normal_ec(jc,jb,je,2) = &
          ptr_patch%edges%primal_normal_cell(ile1,ibe1,1)%v2
        ptr_int%edge_cell_length(jc,jb,je) = &
          ptr_patch%edges%edge_cell_length(ile1,ibe1,1)

      ELSE IF ((ptr_patch%edges%cell_idx(ile1,ibe1,2) == jc) .AND. &
               (ptr_patch%edges%cell_blk(ile1,ibe1,2) == jb)) THEN

        ptr_int%primal_normal_ec(jc,jb,je,1) = &
          ptr_patch%edges%primal_normal_cell(ile1,ibe1,2)%v1
        ptr_int%primal_normal_ec(jc,jb,je,2) = &
          ptr_patch%edges%primal_normal_cell(ile1,ibe1,2)%v2
        ptr_int%edge_cell_length(jc,jb,je) = &
          ptr_patch%edges%edge_cell_length(ile1,ibe1,2)

      ENDIF

    ENDDO

  ENDDO

END DO !block loop
!$OMP END DO

!$OMP END PARALLEL

  DO je = 1, ptr_patch%cell_type
    CALL sync_patch_array(SYNC_C,ptr_patch,ptr_int%primal_normal_ec(:,:,je,1))
    CALL sync_patch_array(SYNC_C,ptr_patch,ptr_int%primal_normal_ec(:,:,je,2))
    CALL sync_patch_array(SYNC_C,ptr_patch,ptr_int%edge_cell_length(:,:,je))
  ENDDO

END SUBROUTINE complete_patchinfo

  !----------------------------------------------------------------------------
  !>
  !! Initializes a tangential plane at each edge midpoint. Necessary  for efficient
  !! calculation of backward trajectories and the corresponding 'upstream area'.
  !!
  !! For FFSL-schemes like Miura it is necessary to calculate backward trajectories
  !! in order to determine an approximation to the upstream area which is advected
  !! across each cell edge during the time step $\delta t$. In our case, this
  !! calculation is perfomed on a plane which is tangent to the edge midpoint.
  !! The coordinate axes point to the local normal and tangential direction.
  !!
  !! The position of additional points on this tangential plane (like the
  !! circumcenters of the neighbour cells, the edge vertices and the edge midpoints of
  !! the corresponding quadrilateral) is precomputed using the gnomonic projection
  !! including a subsequent rotation.
  !!
  !! For a trajectory computation of second order accuracy the computation of scalar
  !! products between primal/dual normals at quadrilateral edges and the inner edge of
  !! the quadrilateral has been added.
  !!
  !! Order of storage for ptr_int%pos_on_tplane_e(nproma,nblks_e,8,2)
  !! pos_on_tplane_e(:,:,1:2,:) :: neighboring cell centers
  !! pos_on_tplane_e(:,:,3:6,:) :: edge midpoints of the quadrilateral
  !! pos_on_tplane_e(:,:,7:8,:) :: edge vertices
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2010-03-12)
  !! Modification by Daniel Reinert, DWD (2010-05-12)
  !! - added projection of edge vertices onto tangential plane
  !!
  SUBROUTINE init_tplane_e (ptr_patch, ptr_int)

    TYPE(t_patch), INTENT(in) :: ptr_patch  !< patch

    TYPE(t_int_state), INTENT(inout) :: ptr_int  !< interpolation state

    REAL(wp) ::                  &    !< geographical coordinates of edge midpoint
      &  xyloc_edge(2)

    REAL(wp) ::                  &    !< geographical coordinates of neighbouring cell
      &  xyloc_n1(2), xyloc_n2(2)     !< centers

    REAL(wp) ::                  &    !< coordinates of neighbouring cell centers on plane
      &  xyloc_plane_n1(2), xyloc_plane_n2(2)

    REAL(wp) ::                  &    !< geographical coordinates of edge midpoints for the
      &  xyloc_quad(4,2)              !< corresponding quadrilateral cell

    REAL(wp) ::                  &    !< coordinates of edge midpoints for the quadrilateral
      &  xyloc_plane_quad(4,2)        !< cell on plane

    REAL(wp) ::                  &    !< geographical coordinates of edge vertices
      &  xyloc_ve(2,2)

    REAL(wp) ::                  &    !< coordinates of edge vertices on plane
      &  xyloc_plane_ve(2,2)

    REAL(wp) ::                  &    !< primal/dual normal in cartesian coordinates
      &  z_nx(3), z_ny(3)

    REAL(wp) :: z_nx_quad(3),    &    !< primal/dual normal at quadrilateral
      &         z_ny_quad(3)          !< edges in cartesian coordinates

    INTEGER :: ilc1, ilc2, ibc1, ibc2 !< line and block indices of neighbour
                                      !< cell centers
    INTEGER :: ilq, ibq               !< line and block indices of quadrilateral edges
    INTEGER :: ilv, ibv               !< line and block indices of edge vertices
    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER :: i_rcstartlev
    INTEGER :: jb, je                 !< loop indices for block and edges
    INTEGER :: ne, nv                 !< loop index for quadrilateral edges and
                                      !< edge vertices
  !-------------------------------------------------------------------------

    CALL message('mo_interpolation:init_tplane_e', '')

    i_rcstartlev = 2

      ! start and end block
    i_startblk = ptr_patch%edges%start_blk(i_rcstartlev,1)
    i_endblk   = ptr_patch%nblks_int_e


    !
    ! compute position of adjacent cell circumcenters, edge vertices and edge
    ! midpoints of the quadrilateral cell on the tangent plane. The gnomonic
    ! projection is used. Note that we first project the points onto a local
    ! geographical (\lambda-\Phi) system. Then we rotate the coordinates of
    ! these points into a new system with coordinate directions normal and
    ! tangential to the edge.

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,je,ne,nv,ilc1,ibc1,ilc2,ibc2,ilq,ibq,ilv,ibv,i_startidx, &
!$OMP            i_endidx,xyloc_edge,xyloc_n1,xyloc_n2,xyloc_plane_n1,       &
!$OMP            xyloc_plane_n2,xyloc_quad,xyloc_plane_quad,xyloc_ve,        &
!$OMP            xyloc_plane_ve)
    DO jb = i_startblk, i_endblk

      CALL get_indices_e(ptr_patch, jb, i_startblk, i_endblk, &
                        i_startidx, i_endidx, i_rcstartlev)

      DO je = i_startidx, i_endidx

        IF(.NOT.ptr_patch%edges%owner_mask(je,jb)) CYCLE

        !
        ! 1. neighboring cell centers
        !

        ! get geographical coordinates of edge midpoint
        xyloc_edge(1) = ptr_patch%edges%center(je,jb)%lon
        xyloc_edge(2) = ptr_patch%edges%center(je,jb)%lat

        ! get line and block indices of neighbour cells
        ilc1 = ptr_patch%edges%cell_idx(je,jb,1)
        ibc1 = ptr_patch%edges%cell_blk(je,jb,1)
        ilc2 = ptr_patch%edges%cell_idx(je,jb,2)
        ibc2 = ptr_patch%edges%cell_blk(je,jb,2)

        ! get geographical coordinates of first cell center
        xyloc_n1(1)   = ptr_patch%cells%center(ilc1,ibc1)%lon
        xyloc_n1(2)   = ptr_patch%cells%center(ilc1,ibc1)%lat

        ! projection first cell center into local \lambda-\Phi-system
        CALL gnomonic_proj( xyloc_edge(1), xyloc_edge(2), xyloc_n1(1), xyloc_n1(2), &! in
          &                 xyloc_plane_n1(1), xyloc_plane_n1(2) )                   ! out


        ! get geographical coordinates of second cell center
        xyloc_n2(1)   = ptr_patch%cells%center(ilc2,ibc2)%lon
        xyloc_n2(2)   = ptr_patch%cells%center(ilc2,ibc2)%lat

        ! projection second cell center into local \lambda-\Phi-system
        CALL gnomonic_proj( xyloc_edge(1), xyloc_edge(2), xyloc_n2(1), xyloc_n2(2), &! in
          &                 xyloc_plane_n2(1), xyloc_plane_n2(2) )                   ! out



        !
        ! 2. Edge midpoints of the quadrilateral
        !
        DO ne = 1,4

          ! get line and block indices of edge midpoints
          ilq = ptr_patch%edges%quad_idx(je,jb,ne)
          ibq = ptr_patch%edges%quad_blk(je,jb,ne)

          ! get geographical coordinates of edge midpoints
          xyloc_quad(ne,1)   = ptr_patch%edges%center(ilq,ibq)%lon
          xyloc_quad(ne,2)   = ptr_patch%edges%center(ilq,ibq)%lat

          ! projection of edge midpoint into local \lambda-\Phi-system
          CALL gnomonic_proj( xyloc_edge(1), xyloc_edge(2), xyloc_quad(ne,1), &! in
            &                 xyloc_quad(ne,2),                               &! in
            &                 xyloc_plane_quad(ne,1), xyloc_plane_quad(ne,2) ) ! out

        END DO


        !
        ! 3. Edge vertices
        !
        DO nv = 1,2

          ! get line and block indices of edge vertices
          ilv = ptr_patch%edges%vertex_idx(je,jb,nv)
          ibv = ptr_patch%edges%vertex_blk(je,jb,nv)

          ! get geographical coordinates of edge vertices
          xyloc_ve(nv,1)   = ptr_patch%verts%vertex(ilv,ibv)%lon
          xyloc_ve(nv,2)   = ptr_patch%verts%vertex(ilv,ibv)%lat

          ! projection of edge vertices into local \lambda-\Phi-system
          CALL gnomonic_proj( xyloc_edge(1), xyloc_edge(2), xyloc_ve(nv,1), &! in
            &                 xyloc_ve(nv,2),                               &! in
            &                 xyloc_plane_ve(nv,1), xyloc_plane_ve(nv,2)   ) ! out

        END DO



        !
        ! 4. rotate these vectors into a new local cartesian system. In this rotated
        !    system the coordinate axes point into the local normal and tangential
        !    direction at each edge.
        !

        ! centers
        !
        ptr_int%pos_on_tplane_e(je,jb,1,1) = re * (                     &
          &     xyloc_plane_n1(1)  * ptr_patch%edges%primal_normal(je,jb)%v1  &
          &   + xyloc_plane_n1(2)  * ptr_patch%edges%primal_normal(je,jb)%v2 )

        ptr_int%pos_on_tplane_e(je,jb,1,2) = re * (                     &
          &     xyloc_plane_n1(1)  * ptr_patch%edges%dual_normal(je,jb)%v1    &
          &   + xyloc_plane_n1(2)  * ptr_patch%edges%dual_normal(je,jb)%v2 )

        ptr_int%pos_on_tplane_e(je,jb,2,1) = re * (                     &
          &     xyloc_plane_n2(1)  * ptr_patch%edges%primal_normal(je,jb)%v1  &
          &   + xyloc_plane_n2(2)  * ptr_patch%edges%primal_normal(je,jb)%v2 )

        ptr_int%pos_on_tplane_e(je,jb,2,2) = re * (                     &
          &     xyloc_plane_n2(1)  * ptr_patch%edges%dual_normal(je,jb)%v1    &
          &   + xyloc_plane_n2(2)  * ptr_patch%edges%dual_normal(je,jb)%v2 )


        ! edges
        !
        DO ne = 1,4
          ptr_int%pos_on_tplane_e(je,jb,2+ne,1) = re * (                       &
            &     xyloc_plane_quad(ne,1)  * ptr_patch%edges%primal_normal(je,jb)%v1  &
            &   + xyloc_plane_quad(ne,2)  * ptr_patch%edges%primal_normal(je,jb)%v2 )

          ptr_int%pos_on_tplane_e(je,jb,2+ne,2) = re * (                       &
            &     xyloc_plane_quad(ne,1)  * ptr_patch%edges%dual_normal(je,jb)%v1    &
            &   + xyloc_plane_quad(ne,2)  * ptr_patch%edges%dual_normal(je,jb)%v2 )
        END DO


        ! vertices
        !
        DO nv = 1,2
          ptr_int%pos_on_tplane_e(je,jb,6+nv,1) = re * (                    &
            &     xyloc_plane_ve(nv,1)  * ptr_patch%edges%primal_normal(je,jb)%v1 &
            &   + xyloc_plane_ve(nv,2)  * ptr_patch%edges%primal_normal(je,jb)%v2 )

          ptr_int%pos_on_tplane_e(je,jb,6+nv,2) = re * (                    &
            &     xyloc_plane_ve(nv,1)  * ptr_patch%edges%dual_normal(je,jb)%v1   &
            &   + xyloc_plane_ve(nv,2)  * ptr_patch%edges%dual_normal(je,jb)%v2 )
        END DO

      ENDDO ! edges
    ENDDO  ! blocks
!$OMP END DO



    !
    ! For each of the 4 rhomboidal edges transform normal and tangential
    ! unit vectors into cartesian system. Then compute dot product
    ! between these unit vectors and the unit vectors of the inner edge.
    ! normalization not necessary fo cartesian vectors since these are
    ! exactly =1.
    !
!$OMP DO PRIVATE(je,jb,ne,ilq,ibq,i_startidx,i_endidx,z_nx,z_ny,z_nx_quad,z_ny_quad)
    DO jb = i_startblk, i_endblk

      CALL get_indices_e(ptr_patch, jb, i_startblk, i_endblk, &
                        i_startidx, i_endidx, i_rcstartlev)

      DO je = i_startidx, i_endidx

        IF(.NOT.ptr_patch%edges%owner_mask(je,jb)) CYCLE

        !
        ! For the current edge transform normal and tangential unit vectors
        ! into cartesian system
        !

        ! transform primal normal to cartesian vector z_nx
        z_nx(:) = ptr_patch%edges%primal_cart_normal(je,jb)%x(:)

        ! transform dual normal to cartesian vector z_ny
        z_ny(:) = ptr_patch%edges%dual_cart_normal(je,jb)%x(:)

        ! for each of the 4 rhomboidal edges transform normal and tangential
        ! unit vectors into cartesian system. Then compute dot product
        ! between these unit vectors and the unit vectors of the inner edge.
        ! normalization not necessary fo cartesian vectors since these are
        ! exactly =1.
        DO ne=1,4
          ilq = ptr_patch%edges%quad_idx(je,jb,ne)
          ibq = ptr_patch%edges%quad_blk(je,jb,ne)

          z_nx_quad(:)=ptr_patch%edges%primal_cart_normal(ilq,ibq)%x(:)
          z_ny_quad(:)=ptr_patch%edges%dual_cart_normal(ilq,ibq)%x(:)

          ! Compute Dot Products
          ptr_int%tplane_e_dotprod(je,jb,ne,1)= DOT_PRODUCT(z_nx_quad(1:3),z_nx(1:3))
          ptr_int%tplane_e_dotprod(je,jb,ne,2)= DOT_PRODUCT(z_ny_quad(1:3),z_nx(1:3))
          ptr_int%tplane_e_dotprod(je,jb,ne,3)= DOT_PRODUCT(z_nx_quad(1:3),z_ny(1:3))
          ptr_int%tplane_e_dotprod(je,jb,ne,4)= DOT_PRODUCT(z_ny_quad(1:3),z_ny(1:3))

        ENDDO
      ENDDO
    ENDDO
!$OMP END DO
!$OMP END PARALLEL

    DO ne=1,8
      call sync_patch_array(SYNC_E, ptr_patch, ptr_int%pos_on_tplane_e(:,:,ne,1))
      call sync_patch_array(SYNC_E, ptr_patch, ptr_int%pos_on_tplane_e(:,:,ne,2))
    ENDDO
    DO ne=1,4
      call sync_patch_array(SYNC_E, ptr_patch, ptr_int%tplane_e_dotprod(:,:,ne,1))
      call sync_patch_array(SYNC_E, ptr_patch, ptr_int%tplane_e_dotprod(:,:,ne,2))
      call sync_patch_array(SYNC_E, ptr_patch, ptr_int%tplane_e_dotprod(:,:,ne,3))
      call sync_patch_array(SYNC_E, ptr_patch, ptr_int%tplane_e_dotprod(:,:,ne,4))
    ENDDO


  END SUBROUTINE init_tplane_e


  !----------------------------------------------------------------------------
  !>
  !! Primal cell quadrature points and weights
  !!
  !! Computes quadrature points and weights for triangular grid cells.
  !! Quadrature points and weights are provided for accurately integrating 
  !! linear, quadratic and cubic functions. This is necessary for initializing 
  !! idealized testcases.
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2010-11-16)
  !!
  !! @par Literature
  !! Numerical Methods in Engineering with Python, Jaan Kiusalaas (2005), 
  !! 233-247
  SUBROUTINE tri_quadrature_pts (ptr_patch, ptr_int)

    TYPE(t_patch), INTENT(IN) :: ptr_patch  !< patch

    TYPE(t_int_state), INTENT(INOUT) :: ptr_int  !< interpolation state

    REAL(wp) ::  alpha_l(3),        & !< area coordinates for quadrature up to 
      &  alpha_q(3,3), alpha_c(3,4)   !< fourth order
                                      !< (n_area_coords,n_pts))

    TYPE(t_cartesian_coordinates)    :: z_vert_cc(3) ! cell vertices in cartesian
                                                     ! coordinates
    TYPE(t_cartesian_coordinates)    :: z_quad_cc  ! triangle quadrature point in cartesian
                                                   ! coordinates
    TYPE(t_geographical_coordinates) :: z_quad_gg  ! triangle quadrature point in geographical
                                                   ! coordinates

    INTEGER :: ilv, ibv               !< line and block indices of cell vertices
    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER :: i_rcstartlev
    INTEGER :: nq, nv                 !< loop index for quadrature points and
                                      !< cell vertices
    INTEGER :: jc, jb                 !< loop index for cells

  !-------------------------------------------------------------------------

    CALL message('mo_interpolation:tri_quadrature_pts', '')

    ! set area coordinates
    !
    ! linear
    alpha_l(1) = 1._wp/3._wp
    alpha_l(2) = 1._wp/3._wp
    alpha_l(3) = 1._wp/3._wp

    !
    ! quadratic
    ! 
    alpha_q(1,1) = 0.5_wp 
    alpha_q(2,1) = 0._wp
    alpha_q(3,1) = 0.5_wp

    alpha_q(1,2) = 0.5_wp
    alpha_q(2,2) = 0.5_wp
    alpha_q(3,2) = 0._wp

    alpha_q(1,3) = 0._wp
    alpha_q(2,3) = 0.5_wp
    alpha_q(3,3) = 0.5_wp
 
    !
    ! cubic
    !
    alpha_c(1,1) = 1._wp/3._wp 
    alpha_c(2,1) = 1._wp/3._wp
    alpha_c(3,1) = 1._wp/3._wp

    alpha_c(1,2) = 1._wp/5._wp
    alpha_c(2,2) = 1._wp/5._wp
    alpha_c(3,2) = 3._wp/5._wp

    alpha_c(1,3) = 3._wp/5._wp
    alpha_c(2,3) = 1._wp/5._wp
    alpha_c(3,3) = 1._wp/5._wp

    alpha_c(1,4) = 1._wp/5._wp
    alpha_c(2,4) = 3._wp/5._wp
    alpha_c(3,4) = 1._wp/5._wp

    ! note that the linear weighting factor is 1 (not stored)
    ptr_int%gquad%weights_tri_q(1:3) = 1._wp/3._wp
    ptr_int%gquad%weights_tri_c(1)   = -27._wp/48._wp
    ptr_int%gquad%weights_tri_c(2:4) = 25._wp/48._wp


    i_rcstartlev = 2

    ! start and end block
    i_startblk = ptr_patch%cells%start_blk(i_rcstartlev,1)
    i_endblk   = ptr_patch%nblks_int_c

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,nv,nq,i_startidx,i_endidx,ilv,ibv,z_vert_cc,z_quad_cc, &
!$OMP           z_quad_gg)
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
                        i_startidx, i_endidx, i_rcstartlev)

      DO jc = i_startidx, i_endidx

        IF(.NOT.ptr_patch%cells%owner_mask(jc,jb)) CYCLE

       ! loop over triangle vertices
!CDIR EXPAND=3
       DO nv=1,3
         ! get line and block indices of cell vertices
         ilv= ptr_patch%cells%vertex_idx(jc,jb,nv)
         ibv= ptr_patch%cells%vertex_blk(jc,jb,nv)

         ! Transform geographical coordinates to cartesian coordinates for vertices
          z_vert_cc(nv)=gc2cc(ptr_patch%verts%vertex(ilv,ibv))
       ENDDO

       !
       ! Linear
       !
       ! Compute quadrature point in cartesian coordinates (= triangle centroid)
       ! i.e. map area coordinates into cartesian coordinates
       z_quad_cc%x(1)= alpha_l(1)*z_vert_cc(1)%x(1)  &
         &           + alpha_l(2)*z_vert_cc(2)%x(1)  &
         &           + alpha_l(3)*z_vert_cc(3)%x(1)

       z_quad_cc%x(2)= alpha_l(1)*z_vert_cc(1)%x(2)  &
         &           + alpha_l(2)*z_vert_cc(2)%x(2)  &
         &           + alpha_l(3)*z_vert_cc(3)%x(2)

       z_quad_cc%x(3)= alpha_l(1)*z_vert_cc(1)%x(3)  &
         &           + alpha_l(2)*z_vert_cc(2)%x(3)  &
         &           + alpha_l(3)*z_vert_cc(3)%x(3)


       ! Transform back to geographical coordinates
       z_quad_gg = cc2gc(z_quad_cc)

       ! store
       ptr_int%gquad%qpts_tri_l(jc,jb)%lat = z_quad_gg%lat
       ptr_int%gquad%qpts_tri_l(jc,jb)%lon = z_quad_gg%lon


       !
       ! quadratic
       !
       ! Loop over quadrature points
!CDIR EXPAND=3
       DO nq=1,3
         ! map area coordinates into cartesian coordinates
         z_quad_cc%x(1)= alpha_q(1,nq)*z_vert_cc(1)%x(1)  &
           &           + alpha_q(2,nq)*z_vert_cc(2)%x(1)  &
           &           + alpha_q(3,nq)*z_vert_cc(3)%x(1)

         z_quad_cc%x(2)= alpha_q(1,nq)*z_vert_cc(1)%x(2)  &
           &           + alpha_q(2,nq)*z_vert_cc(2)%x(2)  &
           &           + alpha_q(3,nq)*z_vert_cc(3)%x(2)

         z_quad_cc%x(3)= alpha_q(1,nq)*z_vert_cc(1)%x(3)  &
           &           + alpha_q(2,nq)*z_vert_cc(2)%x(3)  &
           &           + alpha_q(3,nq)*z_vert_cc(3)%x(3)


         ! Transform back to geographical coordinates
         z_quad_gg = cc2gc(z_quad_cc)

         ! store
         ptr_int%gquad%qpts_tri_q(jc,jb,nq)%lat = z_quad_gg%lat
         ptr_int%gquad%qpts_tri_q(jc,jb,nq)%lon = z_quad_gg%lon
       ENDDO


       !
       ! cubic
       !
       ! Loop over quadrature points
!CDIR EXPAND=4
       DO nq=1,4
         ! map area coordinates into cartesian coordinates
         z_quad_cc%x(1)= alpha_c(1,nq)*z_vert_cc(1)%x(1)  &
           &           + alpha_c(2,nq)*z_vert_cc(2)%x(1)  &
           &           + alpha_c(3,nq)*z_vert_cc(3)%x(1)

         z_quad_cc%x(2)= alpha_c(1,nq)*z_vert_cc(1)%x(2)  &
           &           + alpha_c(2,nq)*z_vert_cc(2)%x(2)  &
           &           + alpha_c(3,nq)*z_vert_cc(3)%x(2)

         z_quad_cc%x(3)= alpha_c(1,nq)*z_vert_cc(1)%x(3)  &
           &           + alpha_c(2,nq)*z_vert_cc(2)%x(3)  &
           &           + alpha_c(3,nq)*z_vert_cc(3)%x(3)


         ! Transform back to geographical coordinates
         z_quad_gg = cc2gc(z_quad_cc)

         ! store
         ptr_int%gquad%qpts_tri_c(jc,jb,nq)%lat = z_quad_gg%lat
         ptr_int%gquad%qpts_tri_c(jc,jb,nq)%lon = z_quad_gg%lon
       ENDDO

      END DO
    ENDDO
!$OMP END DO
!$OMP END PARALLEL

    CALL sync_patch_array(SYNC_C,ptr_patch,ptr_int%gquad%qpts_tri_l(:,:)%lat)
    CALL sync_patch_array(SYNC_C,ptr_patch,ptr_int%gquad%qpts_tri_l(:,:)%lon)
    DO nq=1,3
      CALL sync_patch_array(SYNC_C,ptr_patch,ptr_int%gquad%qpts_tri_q(:,:,nq)%lat)
      CALL sync_patch_array(SYNC_C,ptr_patch,ptr_int%gquad%qpts_tri_q(:,:,nq)%lon)
    ENDDO
    DO nq=1,4
      CALL sync_patch_array(SYNC_C,ptr_patch,ptr_int%gquad%qpts_tri_c(:,:,nq)%lat)
      CALL sync_patch_array(SYNC_C,ptr_patch,ptr_int%gquad%qpts_tri_c(:,:,nq)%lon)
    ENDDO

  END SUBROUTINE tri_quadrature_pts

  !-------------------------------------------------------------------------
  !
  !>
  !! Computes the coefficients that determine the scalar product on the primal grid. This
  !! scalar product depends on the grid geometry only and  is used to formulate the primitive
  !! equations in weak form. The coefficients are applied in module "mo_scalar_product".
  !! The following components of the data type "ocean_patch" are filled:
  !!   edge2cell_coeff  : coefficients for edge to cell mapping
  !!   edge2cell_coeff_t: coefficients for transposed of edge to cell mappings
  !!   edge2vert_coeff  : coefficients for edge to vertex mapping
  !!   edge2vert_coeff_t: coefficients for transposed of edge to vertex mappings
  !!   fixed_vol_norm   : summed volume weight of moved cell
  !!   variable_vol_norm: volume weight at the edges of moved cell
  !!
  !! @par Revision History
  !!  developed by Peter Korn, MPI-M  2010-09
  !!  Modification by Stephan Lorenz, 2010-11
  !!
  SUBROUTINE init_scalar_product_oce( ptr_patch, ptr_intp)

    !  patch on which computation is performed
    !
    TYPE(t_patch)    , TARGET, INTENT(INOUT) :: ptr_patch
    TYPE(t_int_state),         INTENT(INOUT) :: ptr_intp

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
    & routine = ('mo_intp_coeffs:init_scalar_product_oce')

    INTEGER, PARAMETER :: no_cell_edges = 3
    INTEGER, PARAMETER :: no_vert_edges = 6
    INTEGER :: jb, je, jv, ie, ie_1, ie_2, icc
    INTEGER :: il_e,ib_e,k  
    INTEGER :: il_c1, ib_c1, il_c2, ib_c2
    INTEGER :: il_v1, il_v2, ib_v1, ib_v2

    INTEGER :: iil_c1(no_cell_edges), iil_c2(no_cell_edges)
    INTEGER :: iib_c1(no_cell_edges), iib_c2(no_cell_edges)

    INTEGER :: jil_c1, jib_c1,jil_c2, jib_c2
    INTEGER :: rl_start_e, rl_end_e
    INTEGER :: i_startblk_e, i_endblk_e, i_startidx_e, i_endidx_e
    INTEGER :: rl_start, rl_end
    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx

    REAL(wp) :: z_tmp
    REAL(wp) :: norm_c1_c2, norm_v1_v2, norm
    REAL(wp) :: dual_edge_length(no_vert_edges)
    REAL(wp) :: vert_edge_dist(no_vert_edges,2)
    REAL(wp) :: vert_dual_mid_dist(no_vert_edges,2)
    REAL(wp) :: vert_edge_distance, vert_dual_mid_distance

    TYPE(t_geographical_coordinates) :: gc_mid_dual_edge(no_vert_edges)
    TYPE(t_geographical_coordinates) :: gc1,gc2

    TYPE(t_cartesian_coordinates)    :: cc_dual_edge(no_vert_edges), cc_edge(no_cell_edges)
    TYPE(t_cartesian_coordinates)    :: xx1,xx2
    TYPE(t_cartesian_coordinates)    :: vert1_midedge_cc(nproma,ptr_patch%nblks_v,no_vert_edges)
    TYPE(t_cartesian_coordinates)    :: vert2_midedge_cc(nproma,ptr_patch%nblks_v,no_vert_edges)
    TYPE(t_cartesian_coordinates)    :: cell2cell_cc
    TYPE(t_cartesian_coordinates)    :: cc_e0, cc_c1,cc_c2,cc_v0
    TYPE(t_cartesian_coordinates)    :: cv_c1_e0, cv_c2_e0, cv_c1_c2
    TYPE(t_cartesian_coordinates)    :: cc_mid_dual_edge(no_vert_edges)
    TYPE(t_cartesian_coordinates)    :: recon_vec_cc
    TYPE(t_cartesian_coordinates)    :: z_vec_c1(no_cell_edges),z_vec_c2(no_cell_edges)
    TYPE(t_cartesian_coordinates)    :: recon_vec_cc_v1(no_vert_edges) 
    TYPE(t_cartesian_coordinates)    :: recon_vec_cc_v2(no_vert_edges)

    REAL(wp) :: z_edge_length(no_cell_edges)
    REAL(wp) :: z_cell_edge_dist_c1(no_cell_edges,2),z_cell_edge_dist_c2(no_cell_edges,2)
    REAL(wp) :: z_y

    REAL(wp) :: z_sync_c(nproma,ptr_patch%nblks_c)
    REAL(wp) :: z_sync_e(nproma,ptr_patch%nblks_e)
    REAL(wp) :: z_sync_v(nproma,ptr_patch%nblks_v)

    LOGICAL, PARAMETER :: MID_POINT_DUAL_EDGE = .TRUE. !Please do not change this unless
                                                       !you are sure, you know what you do.
    LOGICAL, PARAMETER :: LARC_LENGTH = .FALSE.
    !-----------------------------------------------------------------------
    CALL message (TRIM(routine), 'start')

    rl_start     = 1
    rl_end       = min_rlcell_int
    i_startblk   = ptr_patch%cells%start_blk(rl_start,1)
    i_endblk     = ptr_patch%cells%end_blk(rl_end,1)

    rl_start_e   = 1
    rl_end_e     = min_rledge_int
    i_startblk_e = ptr_patch%edges%start_blk(rl_start_e,1)
    i_endblk_e   = ptr_patch%edges%end_blk(rl_end_e,1)

    !-----------------------------------------------------------------------
    !STEP 1: edge2cell and cell2edge coefficients
    EDGE_BLK_LOOP_PRIMAL: DO jb = i_startblk_e, i_endblk_e

      CALL get_indices_e(ptr_patch, jb,&
                       & i_startblk_e, i_endblk_e,&
                       & i_startidx_e, i_endidx_e,&
                       & rl_start_e, rl_end_e)

      EDGE_IDX_LOOP_PRIMAL: DO je =  i_startidx_e, i_endidx_e

        !Get indices of two adjacent triangles
        il_c1 = ptr_patch%edges%cell_idx(je,jb,1)
        ib_c1 = ptr_patch%edges%cell_blk(je,jb,1)
        il_c2 = ptr_patch%edges%cell_idx(je,jb,2)
        ib_c2 = ptr_patch%edges%cell_blk(je,jb,2)

        cc_c1 = gc2cc(ptr_patch%cells%center(il_c1, ib_c1))
        cc_c2 = gc2cc(ptr_patch%cells%center(il_c2, ib_c2))

        z_cell_edge_dist_c1 = 0.0_wp
        z_cell_edge_dist_c2 = 0.0_wp

        !normals in cell 1
        DO ie = 1, no_cell_edges

          !actual edges of cell c1
          iil_c1(ie) = ptr_patch%cells%edge_idx(il_c1,ib_c1,ie)
          iib_c1(ie) = ptr_patch%cells%edge_blk(il_c1,ib_c1,ie)

          cc_edge(ie)   = gc2cc(ptr_patch%edges%center(iil_c1(ie),iib_c1(ie)))

          !calculate edge length
          !get vertex indices adjacent to actual edge
          il_v1 = ptr_patch%edges%vertex_idx(iil_c1(ie),iib_c1(ie),1)
          ib_v1 = ptr_patch%edges%vertex_blk(iil_c1(ie),iib_c1(ie),1)
          il_v2 = ptr_patch%edges%vertex_idx(iil_c1(ie),iib_c1(ie),2)
          ib_v2 = ptr_patch%edges%vertex_blk(iil_c1(ie),iib_c1(ie),2)

          !get vertex positions
          xx1 = gc2cc(ptr_patch%verts%vertex(il_v1,ib_v1))
          xx2 = gc2cc(ptr_patch%verts%vertex(il_v2,ib_v2))
 
          IF(LARC_LENGTH)THEN
            norm=SQRT(SUM(xx1%x*xx1%x))
            xx1%x= xx1%x/norm
            norm=SQRT(SUM(xx2%x*xx2%x))
            xx2%x= xx2%x/norm
            z_edge_length(ie) = arc_length(xx2,xx1)
            !z_edge_length(ie) = ptr_patch%edges%primal_edge_length(iil_c1(ie),iib_c1(ie))/re
          ELSE
            z_edge_length(ie) = SQRT(SUM((xx2%x-xx1%x)*(xx2%x-xx1%x)))
           !write(*,*)'length 3:', z_edge_length(ie),arc_length(xx1,xx2) 
          ENDIF

          !calculate cell-edge distance as half of cell-cell distance
          !get cell indices adjacent to actual edge
          jil_c1 = ptr_patch%edges%cell_idx(iil_c1(ie),iib_c1(ie),1)
          jib_c1 = ptr_patch%edges%cell_blk(iil_c1(ie),iib_c1(ie),1)
          jil_c2 = ptr_patch%edges%cell_idx(iil_c1(ie),iib_c1(ie),2)
          jib_c2 = ptr_patch%edges%cell_blk(iil_c1(ie),iib_c1(ie),2)

          !get cell positions
          xx1 = gc2cc(ptr_patch%cells%center(jil_c1,jib_c1))
          xx2 = gc2cc(ptr_patch%cells%center(jil_c2,jib_c2))

          IF(jil_c1==il_c1.AND.jib_c1==ib_c1)THEN
           k=1
          ELSEIF(jil_c2==il_c1.AND.jib_c2==ib_c1)THEN
           k=2
          ENDIF

          IF(LARC_LENGTH)THEN
            norm=SQRT(SUM(xx1%x*xx1%x))
            xx1%x= xx1%x/norm
            norm=SQRT(SUM(xx2%x*xx2%x))
            xx2%x= xx2%x/norm
            norm          = SQRT(SUM(cc_edge(ie)%x*cc_edge(ie)%x))
            cc_edge(ie)%x = cc_edge(ie)%x/norm
            z_cell_edge_dist_c1(ie,1) = arc_length(cc_edge(ie),xx1)
            z_cell_edge_dist_c1(ie,2) = arc_length(cc_edge(ie),xx2)
            !z_cell_edge_dist_c1(ie,1) = ptr_patch%edges%edge_cell_length(iil_c1(ie),iib_c1(ie),k)/re 
            !z_cell_edge_dist_c1(ie,2) = ptr_patch%edges%edge_cell_length(iil_c1(ie),iib_c1(ie),k)/re 
          ELSE
            z_cell_edge_dist_c1(ie,1) = SQRT(SUM((cc_edge(ie)%x-xx1%x)*(cc_edge(ie)%x-xx1%x)))
            z_cell_edge_dist_c1(ie,2) = SQRT(SUM((cc_edge(ie)%x-xx2%x)*(cc_edge(ie)%x-xx2%x)))  
            !write(*,*)'length 4',z_cell_edge_dist_c1(ie,1), z_cell_edge_dist_c1(ie,1),&
            !&ptr_patch%edges%edge_cell_length(iil_c1(ie),iib_c1(ie),1)/re 
          ENDIF
          ptr_intp%dist_cell2edge(iil_c1(ie),iib_c1(ie),1) = z_cell_edge_dist_c1(ie,1)
          ptr_intp%dist_cell2edge(iil_c1(ie),iib_c1(ie),2) = z_cell_edge_dist_c1(ie,2)

          z_vec_c1(ie)%x = cc_edge(ie)%x - cc_c1%x     !ptr_patch%edges%primal_cart_normal(iil_c1(ie),iib_c1(ie))
          norm           = SQRT(SUM( z_vec_c1(ie)%x* z_vec_c1(ie)%x))

          ptr_intp%edge2cell_coeff_cc(il_c1,ib_c1,ie)%x = &
            & z_vec_c1(ie)%x*ptr_patch%cells%edge_orientation(il_c1,ib_c1,ie)*z_edge_length(ie)

          ptr_intp%fixed_vol_norm(il_c1,ib_c1)       = ptr_intp%fixed_vol_norm(il_c1,ib_c1) + &
            &                                        0.5_wp*norm*z_edge_length(ie)
          ptr_intp%variable_vol_norm(il_c1,ib_c1,ie) = 0.5_wp*norm*z_edge_length(ie)

          !write(*,*)'edge length   :',z_edge_length(ie),ptr_patch%edges%primal_edge_length(iil_c1(ie),iib_c1(ie))/re
          !write(*,*)'cell-edge dist:', z_cell_edge_dist_c1(ie,k),ptr_patch%edges%edge_cell_length(iil_c1(ie),iib_c1(ie),k)/re
        END DO

        !normals in cell 2
        DO ie = 1, no_cell_edges

          !actual edges of cell c2
          iil_c2(ie) = ptr_patch%cells%edge_idx(il_c2,ib_c2,ie)
          iib_c2(ie) = ptr_patch%cells%edge_blk(il_c2,ib_c2,ie)
          !write(0,*)'iil_c2(ie):',iil_c2(ie),' iib_c2(ie):',iib_c2(ie)


          cc_edge(ie) = gc2cc(ptr_patch%edges%center(iil_c2(ie),iib_c2(ie)))

          !calculate edge length
          !get vertex indices adjacent to actual edge
          il_v1 = ptr_patch%edges%vertex_idx(iil_c2(ie),iib_c2(ie),1)
          ib_v1 = ptr_patch%edges%vertex_blk(iil_c2(ie),iib_c2(ie),1)
          il_v2 = ptr_patch%edges%vertex_idx(iil_c2(ie),iib_c2(ie),2)
          ib_v2 = ptr_patch%edges%vertex_blk(iil_c2(ie),iib_c2(ie),2)

          !get vertex positions
          xx1 = gc2cc(ptr_patch%verts%vertex(il_v1,ib_v1))
          xx2 = gc2cc(ptr_patch%verts%vertex(il_v2,ib_v2))

          IF(LARC_LENGTH)THEN
            norm              = SQRT(SUM(xx1%x*xx1%x))
            xx1%x             = xx1%x/norm

            norm              = SQRT(SUM(xx2%x*xx2%x))
            xx2%x             = xx2%x/norm

            z_edge_length(ie) = arc_length(xx2,xx1)
            !z_edge_length(ie) = ptr_patch%edges%primal_edge_length(iil_c2(ie),iib_c2(ie))/re
             !write(*,*)'arc length',arc_length(xx2,xx1),z_edge_length(ie),SQRT(SUM((xx2%x-xx1%x)*(xx2%x-xx1%x)))
          ELSE
            z_edge_length(ie) = SQRT(SUM((xx2%x-xx1%x)*(xx2%x-xx1%x)))
          ENDIF
          !calculate cell-edge distance as half of cell-cell distance
          !get cell indices adjacent to actual edge
          jil_c1 = ptr_patch%edges%cell_idx(iil_c2(ie),iib_c2(ie),1)
          jib_c1 = ptr_patch%edges%cell_blk(iil_c2(ie),iib_c2(ie),1)
          jil_c2 = ptr_patch%edges%cell_idx(iil_c2(ie),iib_c2(ie),2)
          jib_c2 = ptr_patch%edges%cell_blk(iil_c2(ie),iib_c2(ie),2)

          !write(0,*)'jil_c1:',jil_c1,' jib_c1:',jib_c1,' jil_c2:',jil_c2,' jib_c2:',jib_c2
          if (jil_c2 < 0) THEN
            write(0,*)'ptr_patch%edges%cell_idx:',ptr_patch%edges%cell_idx
          ENDIF
          !get cell positions
          xx1 = gc2cc(ptr_patch%cells%center(jil_c1,jib_c1))
          xx2 = gc2cc(ptr_patch%cells%center(jil_c2,jib_c2))

          IF(jil_c1==il_c2.AND.jib_c1==ib_c2)THEN
            k=1
          ELSEIF(jil_c2==il_c2.AND.jib_c2==ib_c2)THEN
            k=2
          ENDIF  

          IF(LARC_LENGTH)THEN
            norm=SQRT(SUM(xx1%x*xx1%x))
            xx1%x= xx1%x/norm
            norm=SQRT(SUM(xx2%x*xx2%x))
            xx2%x= xx2%x/norm
            norm=SQRT(SUM(cc_edge(ie)%x*cc_edge(ie)%x))
            cc_edge(ie)%x =  cc_edge(ie)%x/norm
            z_cell_edge_dist_c2(ie,1) = arc_length(cc_edge(ie),xx1)
            z_cell_edge_dist_c2(ie,2) = arc_length(cc_edge(ie),xx2)
            !z_cell_edge_dist_c2(ie,1) = ptr_patch%edges%edge_cell_length(iil_c2(ie),iib_c2(ie),1)/re
            !z_cell_edge_dist_c2(ie,2) = ptr_patch%edges%edge_cell_length(iil_c2(ie),iib_c2(ie),2)/re
            !write(*,*)'arc length',0.5_wp*arc_length(xx2,xx1),z_cell_edge_dist_c2(ie,k),0.5_wp*SQRT(SUM((xx2%x-xx1%x)*(xx2%x-xx1%x)))
          ELSE 
            z_cell_edge_dist_c2(ie,1) = SQRT(SUM((cc_edge(ie)%x-xx1%x)*(cc_edge(ie)%x-xx1%x)))
            z_cell_edge_dist_c2(ie,2) = SQRT(SUM((cc_edge(ie)%x-xx2%x)*(cc_edge(ie)%x-xx2%x)))
          ENDIF
          ptr_intp%dist_cell2edge(iil_c2(ie),iib_c2(ie),1) = z_cell_edge_dist_c1(ie,1)
          ptr_intp%dist_cell2edge(iil_c2(ie),iib_c2(ie),2) = z_cell_edge_dist_c1(ie,2)

          z_vec_c2(ie)%x = cc_edge(ie)%x - cc_c2%x  !ptr_patch%edges%primal_cart_normal(iil_c2(ie),iib_c2(ie))
          norm           = SQRT(SUM( z_vec_c2(ie)%x* z_vec_c2(ie)%x))

          ptr_intp%edge2cell_coeff_cc(il_c2,ib_c2,ie)%x&
            & = z_vec_c2(ie)%x*ptr_patch%cells%edge_orientation(il_c2,ib_c2,ie)*z_edge_length(ie)

          ptr_intp%fixed_vol_norm(il_c2,ib_c2)       = ptr_intp%fixed_vol_norm(il_c2,ib_c2) + &
            &                                        0.5_wp*norm*z_edge_length(ie)
          ptr_intp%variable_vol_norm(il_c2,ib_c2,ie) = 0.5_wp*norm*z_edge_length(ie)

        END DO
      END DO EDGE_IDX_LOOP_PRIMAL
    END DO EDGE_BLK_LOOP_PRIMAL
    !In he edge loop above each triangle is visisted three times. Since the "fixed_vol_norm" is
    !accumulated we correct its value here:
    ptr_intp%fixed_vol_norm = ptr_intp%fixed_vol_norm/3.0_wp

    rl_start   = 1
    rl_end     = min_rledge_int
    i_startblk = ptr_patch%edges%start_blk(rl_start,1)
    i_endblk   = ptr_patch%edges%end_blk(rl_end,1)

    EDGE_BLK_LOOP_SECONDARY: DO jb = i_startblk, i_endblk
      CALL get_indices_e(ptr_patch, jb, i_startblk, i_endblk,&
                       & i_startidx, i_endidx, rl_start, rl_end)
      EDGE_IDX_LOOP_SECONDARY: DO je =  i_startidx, i_endidx

        !Get indices of two adjacent triangles
        il_c1 = ptr_patch%edges%cell_idx(je,jb,1)
        ib_c1 = ptr_patch%edges%cell_blk(je,jb,1)
        il_c2 = ptr_patch%edges%cell_idx(je,jb,2)
        ib_c2 = ptr_patch%edges%cell_blk(je,jb,2)

        !cartesian coordinates of edge and neighbor cells on 1-sphere
        cc_e0 = gc2cc(ptr_patch%edges%center(je,jb))
        cc_c1 = gc2cc(ptr_patch%cells%center(il_c1,ib_c1))
        cc_c2 = gc2cc(ptr_patch%cells%center(il_c2,ib_c2))

        !cartesian vectors from:
        !cell 2 to cell 1, cell 1 to edge je and cell 2 to edge je
        cv_c1_c2%x = cc_c1%x - cc_c2%x
        cv_c1_e0%x = cc_e0%x - cc_c1%x
        cv_c2_e0%x = cc_e0%x - cc_c2%x

        IF(LARC_LENGTH)THEN
          norm=SQRT(SUM(cc_e0%x*cc_e0%x))
          cc_e0%x= cc_e0%x/norm
          norm=SQRT(SUM(cc_c1%x*cc_c1%x))
          cc_c1%x= cc_c1%x/norm
          norm=SQRT(SUM(cc_c2%x*cc_c2%x))
          cc_c2%x= cc_c2%x/norm
          norm_c1_c2 = arc_length(cc_c1, cc_c2)
          !norm_c1_c2 = ptr_patch%edges%dual_edge_length(je,jb)/re 
                       !ptr_patch%edges%edge_cell_length(je,jb,1)/re&
                      !&+ptr_patch%edges%edge_cell_length(je,jb,2)/re
        ELSE
          norm_c1_c2 = SQRT(SUM(cv_c1_e0%x*cv_c1_e0%x))+SQRT(SUM(cv_c2_e0%x*cv_c2_e0%x)) !SQRT(SUM(cv_c1_c2%x*cv_c1_c2%x))!
        ENDIF

        !Determine which edge of both of the two adjacent cells corresponds to the
        !actual edge "je". This information is used below for the edge-orientation.
        DO ie = 1, no_cell_edges
          IF (ptr_patch%cells%edge_idx(il_c1,ib_c1,ie) == je.AND.&
            & ptr_patch%cells%edge_blk(il_c1,ib_c1,ie) == jb) THEN
            ie_1 = ie
          END IF
          IF (ptr_patch%cells%edge_idx(il_c2,ib_c2,ie) == je.AND.&
            & ptr_patch%cells%edge_blk(il_c2,ib_c2,ie) == jb) THEN
            ie_2 = ie
          END IF
        END DO

        ptr_intp%edge2cell_coeff_cc_t(je,jb,1)%x&
        & = cv_c1_e0%x * ptr_patch%cells%edge_orientation(il_c1,ib_c1,ie_1)/norm_c1_c2

        ptr_intp%edge2cell_coeff_cc_t(je,jb,2)%x&
        & = cv_c2_e0%x * ptr_patch%cells%edge_orientation(il_c2,ib_c2,ie_2)/norm_c1_c2

      END DO EDGE_IDX_LOOP_SECONDARY
    END DO EDGE_BLK_LOOP_SECONDARY

    !------------------------------------------------------------------------------
    !STEP 2: edge2vert coefficients for dual grid
    !------------------------------------------------------------------------------

    rl_start = 1
    rl_end   = min_rlvert
   rl_end   = min_rlvert_int  ! inner part of decomposition only - no halo (!!)

    i_startblk = ptr_patch%verts%start_blk(rl_start,1)
    i_endblk   = ptr_patch%verts%end_blk(rl_end,1)

    VERT_BLK_LOOP: DO jb = i_startblk, i_endblk
      CALL get_indices_v(ptr_patch, jb, i_startblk, i_endblk,&
                       & i_startidx, i_endidx, rl_start, rl_end)
      VERT_IDX_LOOP: DO jv =  i_startidx, i_endidx
        ! current number of edges around vertex (5 or 6)
        cc_v0        = gc2cc(ptr_patch%verts%vertex(jv,jb))

        DO ie = 1, no_vert_edges  ! #slo# it_vertedges ??

          il_e = ptr_patch%verts%edge_idx(jv,jb,ie)
          ib_e = ptr_patch%verts%edge_blk(jv,jb,ie)

          ! #slo# - I assume all geographical coordinates are already synchronized

          cc_dual_edge(ie) = gc2cc(ptr_patch%edges%center(il_e,ib_e))
          !Parts of this code parrallels the implementation in the grid-generator
          !module "mo_geometry".
          !
          !1) determine normal vector from adjacent cell to adjacent cell
          !   in cartesian coordinate for moved dual cell
          !Get indices of two adjacent triangles
          il_c1 = ptr_patch%edges%cell_idx(il_e,ib_e,1)
          ib_c1 = ptr_patch%edges%cell_blk(il_e,ib_e,1)
          il_c2 = ptr_patch%edges%cell_idx(il_e,ib_e,2)
          ib_c2 = ptr_patch%edges%cell_blk(il_e,ib_e,2)

          xx1 = gc2cc(ptr_patch%cells%center(il_c1,ib_c1))
          norm=SQRT(SUM(xx1%x*xx1%x))
          xx1%x= xx1%x/norm

          xx2 = gc2cc(ptr_patch%cells%center(il_c2,ib_c2))
          norm=SQRT(SUM(xx2%x*xx2%x))
          xx2%x= xx2%x/norm

          cell2cell_cc%x       = xx2%x - xx1%x
          IF(LARC_LENGTH)THEN
            norm_c1_c2 = arc_length(xx1,xx2)
          ELSE
            norm_c1_c2 = SQRT(SUM(cell2cell_cc%x*cell2cell_cc%x))
          ENDIF
          dual_edge_length(ie) = norm_c1_c2
!          cell2cell_cc%x       = cell2cell_cc%x/norm_c1_c2

          IF(MID_POINT_DUAL_EDGE)THEN
            cc_mid_dual_edge(ie)%x = 0.5_wp*(xx2%x+xx1%x)
            gc_mid_dual_edge(ie)   = cc2gc(cc_mid_dual_edge(ie))

            IF(CORIOLIS_TYPE==full_coriolis)THEN
              ptr_patch%edges%f_e(il_e, ib_e) = 2._wp*omega*SIN(gc_mid_dual_edge(ie)%lat)
            ELSEIF(CORIOLIS_TYPE==BETA_PLANE_CORIOLIS)THEN
              gc1%lat = basin_center_lat* deg2rad - 0.5_wp*basin_height_deg*deg2rad
              gc1%lon = 0.0_wp
              xx1     = gc2cc(gc1)

              gc2%lat = gc_mid_dual_edge(ie)%lat!*deg2rad
              gc2%lon = 0.0_wp
              xx2     = gc2cc(gc2)
              z_y     = re*arc_length(xx2,xx1)
 
              !z_y = ptr_patch%edges%center(je,jb)%lat - z_lat_basin_center
              ptr_patch%edges%f_e(il_e, ib_e) = 2.0_wp*omega*( sin(basin_center_lat * deg2rad) &
                &                             + (cos(basin_center_lat * deg2rad)/re)*z_y)
            ENDIF
          ELSE
            cc_mid_dual_edge(ie)%x = cc_dual_edge(ie)%x
            gc_mid_dual_edge(ie)   = cc2gc(cc_mid_dual_edge(ie)) 
          ENDIF

          !2) determine vector from adjacent vertex to adjacent vertex
          !   in cartesian coordinate for moved dual cell
          !Get indices of two adjacent vertices
          il_v1 = ptr_patch%edges%vertex_idx(il_e,ib_e,1)
          ib_v1 = ptr_patch%edges%vertex_blk(il_e,ib_e,1)
          il_v2 = ptr_patch%edges%vertex_idx(il_e,ib_e,2)
          ib_v2 = ptr_patch%edges%vertex_blk(il_e,ib_e,2)

          xx1 = gc2cc(ptr_patch%verts%vertex(il_v1,ib_v1))
          norm=SQRT(SUM(xx1%x*xx1%x))
          xx1%x= xx1%x/norm

          xx2 = gc2cc(ptr_patch%verts%vertex(il_v2,ib_v2))
          norm=SQRT(SUM(xx2%x*xx2%x))
          xx2%x= xx2%x/norm

          vert1_midedge_cc(jv, jb, ie)%x = cc_mid_dual_edge(ie)%x - xx1%x
          vert2_midedge_cc(jv, jb, ie)%x = cc_mid_dual_edge(ie)%x - xx2%x
          !vert2vert_cc(jv, jb, ie)%x = cc_mid_dual_edge(ie)%x - xx1%x !xx2%x - xx1%x
!           norm                       = SQRT(SUM(vert2vert_cc(jv,jb,ie)%x*vert2vert_cc(jv,jb,ie)%x))
!           vert2vert_cc(jv,jb,ie)%x   = vert2vert_cc(jv, jb, ie)%x/norm
          norm = SQRT(SUM(vert1_midedge_cc(jv,jb,ie)%x*vert1_midedge_cc(jv,jb,ie)%x))
          vert1_midedge_cc(jv, jb, ie)%x = vert1_midedge_cc(jv, jb, ie)%x/norm

          norm = SQRT(SUM(vert2_midedge_cc(jv,jb,ie)%x*vert2_midedge_cc(jv,jb,ie)%x))
          vert2_midedge_cc(jv, jb, ie)%x = vert2_midedge_cc(jv, jb, ie)%x/norm


          !calculate vertex edge distance 
          IF(LARC_LENGTH)THEN
            vert_edge_dist(ie,1) = arc_length (cc_dual_edge(ie), xx1) 
            vert_edge_dist(ie,2) = arc_length (cc_dual_edge(ie), xx2) 
            vert_dual_mid_dist(ie,1)= arc_length (cc_mid_dual_edge(ie), xx1) 
            vert_dual_mid_dist(ie,2)= arc_length (cc_mid_dual_edge(ie), xx2)
          ELSE
            vert_edge_dist(ie,1)&
            & = SQRT(SUM((cc_dual_edge(ie)%x - xx1%x)*(cc_dual_edge(ie)%x - xx1%x)))
            vert_edge_dist(ie,2)&
            & = SQRT(SUM((cc_dual_edge(ie)%x - xx2%x)*(cc_dual_edge(ie)%x - xx2%x)))
            vert_dual_mid_dist(ie,1)&
            & = SQRT(SUM((cc_mid_dual_edge(ie)%x - xx1%x)*(cc_mid_dual_edge(ie)%x - xx1%x)))
            vert_dual_mid_dist(ie,2)&
            & = SQRT(SUM((cc_mid_dual_edge(ie)%x - xx2%x)*(cc_mid_dual_edge(ie)%x - xx2%x)))
          ENDIF

          !calculate normal vector that is perpendicular to vertex-vertex- and edge position vector
          !If one uses the edge position vector this results in the moved primal normal. Later
          !edge position vector has to be replaced by the midpoint of the dual edge.
          recon_vec_cc_v1(ie)   = vector_product(vert1_midedge_cc(jv, jb, ie),&
            &                                    cc_mid_dual_edge(ie))
          norm                  = SQRT(SUM(recon_vec_cc_v1(ie)%x*recon_vec_cc_v1(ie)%x))
          recon_vec_cc_v1(ie)%x = recon_vec_cc_v1(ie)%x/norm

          recon_vec_cc_v2(ie)   = vector_product(vert2_midedge_cc(jv,jb,ie),&
            &                                    cc_mid_dual_edge(ie))
          norm                  = SQRT(SUM(recon_vec_cc_v2(ie)%x*recon_vec_cc_v2(ie)%x))
          recon_vec_cc_v2(ie)%x = recon_vec_cc_v2(ie)%x/norm

          !Fix orientation
          z_tmp = DOT_PRODUCT(recon_vec_cc_v1(ie)%x,&
            &                 ptr_patch%edges%primal_cart_normal(il_e,ib_e)%x)
          IF (z_tmp <0._wp) recon_vec_cc_v1(ie)%x = -1._wp * recon_vec_cc_v1(ie)%x

          z_tmp = DOT_PRODUCT(recon_vec_cc_v2(ie)%x,&
            &                 ptr_patch%edges%primal_cart_normal(il_e,ib_e)%x)
          IF (z_tmp <0._wp) recon_vec_cc_v2(ie)%x = -1._wp * recon_vec_cc_v2(ie)%x


          IF ( (ptr_patch%edges%vertex_idx(il_e,ib_e,1) == jv) .and. &
               (ptr_patch%edges%vertex_blk(il_e,ib_e,1) == jb)         ) THEN

            vert_edge_distance     = vert_edge_dist(ie,1)
            vert_dual_mid_distance = vert_dual_mid_dist(ie,1)
            recon_vec_cc           = recon_vec_cc_v1(ie)

            ptr_intp%edge2vert_vector_cc(jv,jb,ie)=&
            &vector_product(vert1_midedge_cc(jv, jb, ie), cc_mid_dual_edge(ie))
            z_tmp = DOT_PRODUCT(ptr_intp%edge2vert_vector_cc(jv,jb,ie)%x,&
            &ptr_patch%edges%primal_cart_normal(il_e,ib_e)%x)
            IF (z_tmp <0._wp) ptr_intp%edge2vert_vector_cc(jv,jb,ie)%x&
            & = -1._wp * ptr_intp%edge2vert_vector_cc(jv,jb,ie)%x


          ELSE IF ( (ptr_patch%edges%vertex_idx(il_e,ib_e,2) == jv) .and. &
                    (ptr_patch%edges%vertex_blk(il_e,ib_e,2) == jb) ) THEN

            vert_edge_distance     = vert_edge_dist(ie,2)
            vert_dual_mid_distance = vert_dual_mid_dist(ie,2)
            recon_vec_cc           = recon_vec_cc_v2(ie)

            ptr_intp%edge2vert_vector_cc(jv,jb,ie)=&
            &vector_product(vert2_midedge_cc(jv, jb, ie), cc_mid_dual_edge(ie))
            z_tmp = DOT_PRODUCT(ptr_intp%edge2vert_vector_cc(jv,jb,ie)%x,&
            & ptr_patch%edges%primal_cart_normal(il_e,ib_e)%x)
            IF (z_tmp <0._wp) ptr_intp%edge2vert_vector_cc(jv,jb,ie)%x =&
            & -1._wp * ptr_intp%edge2vert_vector_cc(jv,jb,ie)%x
          ELSE
            CALL message (TRIM(routine), 'WARNING - vert_edge_distance not found')
            write(*,'(a,7i5)') 'jv, jb, edge, ptr_patch%edges%vertex_idx/blk(il_e,ib_e,1-2)=', &
              &         jv, jb, ie, &
              &         ptr_patch%edges%vertex_idx(il_e,ib_e,1), &
              &         ptr_patch%edges%vertex_blk(il_e,ib_e,1), &
              &         ptr_patch%edges%vertex_idx(il_e,ib_e,2), &
              &         ptr_patch%edges%vertex_blk(il_e,ib_e,2)

          END IF

          ptr_intp%variable_dual_vol_norm(jv,jb,ie) = &
            &                           0.5_wp*dual_edge_length(ie)*vert_dual_mid_distance
                                               !vert_edge_distance*dual_edge_length(ie)!

          ptr_intp%edge2vert_coeff_cc(jv,jb,ie)%x   = &
            &                      recon_vec_cc%x*dual_edge_length(ie)*vert_dual_mid_distance

          norm_v1_v2 = SQRT(SUM(vert1_midedge_cc(jv, jb, ie)%x*vert1_midedge_cc(jv, jb, ie)%x))&
                    &+ SQRT(SUM(vert2_midedge_cc(jv, jb, ie)%x*vert2_midedge_cc(jv, jb, ie)%x))
          
          ptr_intp%edge2vert_coeff_cc_t(il_e,ib_e,1)%x = vert1_midedge_cc(jv, jb, ie)%x * &
            &    ( ptr_patch%edges%system_orientation(il_e,ib_e)/norm_v1_v2 )
          
          ptr_intp%edge2vert_coeff_cc_t(il_e,ib_e,2)%x = vert2_midedge_cc(jv, jb, ie)%x * &
            &    ( ptr_patch%edges%system_orientation(il_e,ib_e)/norm_v1_v2 )

        END DO

      ENDDO VERT_IDX_LOOP
    END DO VERT_BLK_LOOP

!TODOram notnec    !--------------------------------------------------------------------------
!TODOram notnec    ! SYNCHRONIZE ALL ELEMENTS OF V_BASE:
!TODOram notnec    ! synchronize elements on cells
!TODOram notnec    DO ie = 1, no_cell_edges
!TODOram notnec      DO icc = 1, 3
!TODOram notnec        z_sync_c(:,:) =  ptr_intp%edge2cell_coeff_cc(:,:,ie)%x(icc)
!TODOram notnec        CALL sync_patch_array(SYNC_C, ptr_patch, z_sync_c(:,:))
!TODOram notnec        ptr_intp%edge2cell_coeff_cc(:,:,ie)%x(icc) = z_sync_c(:,:)
!TODOram notnec      END DO
!TODOram notnec      z_sync_c(:,:) = ptr_intp%variable_vol_norm(:,:,ie)
!TODOram notnec      CALL sync_patch_array(SYNC_C, ptr_patch, z_sync_c(:,:))
!TODOram notnec      ptr_intp%variable_vol_norm(:,:,ie) = z_sync_c(:,:)
!TODOram notnec    END DO
!TODOram notnec    CALL sync_patch_array(SYNC_C, ptr_patch,ptr_intp%fixed_vol_norm)
!TODOram notnec
!TODOram notnec    ! synchronize elements on edges
!TODOram notnec    DO ie = 1, 2
!TODOram notnec      DO icc = 1, 3
!TODOram notnec        z_sync_e(:,:) =  ptr_intp%edge2vert_coeff_cc_t(:,:,ie)%x(icc)
!TODOram notnec        CALL sync_patch_array(SYNC_E, ptr_patch, z_sync_e(:,:))
!TODOram notnec        ptr_intp%edge2vert_coeff_cc_t(:,:,ie)%x(icc) = z_sync_e(:,:)
!TODOram notnec
!TODOram notnec        z_sync_e(:,:) =  ptr_intp%edge2cell_coeff_cc_t(:,:,ie)%x(icc)
!TODOram notnec        CALL sync_patch_array(SYNC_E, ptr_patch, z_sync_e(:,:))
!TODOram notnec        ptr_intp%edge2cell_coeff_cc_t(:,:,ie)%x(icc) = z_sync_e(:,:)
!TODOram notnec      END DO
!TODOram notnec    END DO
!TODOram notnec
!TODOram notnec    ! synchronize cartesian coordinates on vertices:
!TODOram notnec    DO ie = 1, no_vert_edges
!TODOram notnec      DO icc = 1, 3
!TODOram notnec        z_sync_v(:,:) =  ptr_intp%edge2vert_vector_cc(:,:,ie)%x(icc)
!TODOram notnec        CALL sync_patch_array(SYNC_V, ptr_patch, z_sync_v(:,:))
!TODOram notnec        ptr_intp%edge2vert_vector_cc(:,:,ie)%x(icc) = z_sync_v(:,:)
!TODOram notnec
!TODOram notnec        z_sync_v(:,:) = ptr_intp%edge2vert_coeff_cc(:,:,ie)%x(icc)
!TODOram notnec        CALL sync_patch_array(SYNC_V, ptr_patch, z_sync_v(:,:))
!TODOram notnec        ptr_intp%edge2vert_coeff_cc(:,:,ie)%x(icc) = z_sync_v(:,:)
!TODOram notnec      END DO
!TODOram notnec      z_sync_v(:,:) = ptr_intp%variable_dual_vol_norm(:,:,ie)
!TODOram notnec      CALL sync_patch_array(SYNC_V, ptr_patch, z_sync_v(:,:))
!TODOram notnec      ptr_intp%variable_dual_vol_norm(:,:,ie) = z_sync_v(:,:)
!TODOram notnec    END DO

    CALL message (TRIM(routine), 'end')

  END SUBROUTINE init_scalar_product_oce
!-------------------------------------------------------------------------
  !
  !
  !>
  !! Precomputes the geometrical factors used in the divergence, rotation.
  !!
  !! @par Revision History
  !!  developed by Guenther Zaengl, 2009-03-17
  !!  Modification by Almut Gassmann, 2009-12-19
  !!  - Vorticity is computed on quads in case of the hexagonal grid
  !!  Modification by Almut Gassmann, 2010-02-05
  !!  - Added feature for poor men's 3rd order advection, where a directional
  !!    laplace is needed at the edges.
  !!  Modification by Stephan Lorenz, 2010-06-02
  !!  - Storage moved from int_state into patch_oce since it is static
  !!    geometric information used in the ocean model
  !!  Modification by Peter Korn, 2010-11
  !!  - Calculation of cell area changed to achieve compatibility with
  !!    sw-model (cell area and consequently divergence different)
  !!  Modification by Stephan Lorenz, 2011-07
  !!   - 3-dim structures moved from patch_oce to hydro_ocean_base for parallelization
  !!
  SUBROUTINE init_geo_factors_oce( ptr_patch, ptr_intp )
    !
    IMPLICIT NONE
    !
    !  patch on which computation is performed
    !
    TYPE(t_patch), TARGET, INTENT(INOUT) :: ptr_patch
    TYPE(t_int_state),     INTENT(INOUT) :: ptr_intp
    !

    INTEGER :: jc, jb, je, jv, je1, ie
    INTEGER :: rl_start, rl_end
    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx, i_nchdom

    INTEGER :: ile, ibe, ilc1, ibc1, ilc2, ibc2, ifac, ic, ilnc, ibnc
    INTEGER :: ile1, ibe1,ile2,ibe2,ile3,ibe3
    INTEGER, PARAMETER :: i_cell_type = 3
    !TYPE(cartesian_coordinates)::z_pn_k,z_pn_j
    !REAL(wp) :: z_lon, z_lat, z_nu, z_nv, z_proj
    REAL(wp) :: cell_area

    REAL(wp) :: z_sync_c(nproma,ptr_patch%nblks_c)
    REAL(wp) :: z_sync_e(nproma,ptr_patch%nblks_e)
    REAL(wp) :: z_sync_v(nproma,ptr_patch%nblks_v)

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
    & routine = ('mo_oce_state:init_geo_factors_oce')

    !-----------------------------------------------------------------------

    CALL message (TRIM(routine), 'start')

    i_nchdom   = MAX(1,ptr_patch%n_childdom)


!$OMP PARALLEL PRIVATE(rl_start,rl_end,i_startblk,i_endblk,ifac)
    ! a) Geometrical factor for divergence
    rl_start = 1
    rl_end = min_rlcell

    ! values for the blocking
    i_startblk = ptr_patch%cells%start_blk(rl_start,1)
    i_endblk   = ptr_patch%cells%end_blk(rl_end,i_nchdom)
    !
    ! loop through all patch cells (and blocks)
    !
!$OMP DO PRIVATE(jb,je,jc,i_startidx,i_endidx,ile,ibe)
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk,      &
        &                i_startidx, i_endidx, rl_start, rl_end)

      DO jc = i_startidx, i_endidx

        ile1 = ptr_patch%cells%edge_idx(jc,jb,1)
        ibe1 = ptr_patch%cells%edge_blk(jc,jb,1)
        ile2 = ptr_patch%cells%edge_idx(jc,jb,2)
        ibe2 = ptr_patch%cells%edge_blk(jc,jb,2)
        ile3 = ptr_patch%cells%edge_idx(jc,jb,3)
        ibe3 = ptr_patch%cells%edge_blk(jc,jb,3)

        cell_area =  0.25_wp&
  & *( ptr_patch%edges%primal_edge_length(ile1,ibe1)*ptr_patch%edges%dual_edge_length(ile1,ibe1)&
  &   +ptr_patch%edges%primal_edge_length(ile2,ibe2)*ptr_patch%edges%dual_edge_length(ile2,ibe2)&
  &   +ptr_patch%edges%primal_edge_length(ile3,ibe3)*ptr_patch%edges%dual_edge_length(ile3,ibe3))


       DO je = 1, i_cell_type

          IF (je > ptr_patch%cells%num_edges(jc,jb)) CYCLE ! relevant for hexagons

           ile = ptr_patch%cells%edge_idx(jc,jb,je)
           ibe = ptr_patch%cells%edge_blk(jc,jb,je)

           ptr_intp%geofac_div(jc,je,jb) =                &
         &   ptr_patch%edges%primal_edge_length(ile,ibe) * &
         &   ptr_patch%cells%edge_orientation(jc,jb,je)  / &
         &   ptr_patch%cells%area(jc,jb)

           ptr_intp%geofac_div(jc,je,jb) =                &
         &   ptr_patch%edges%primal_edge_length(ile,ibe) * &
         &   ptr_patch%cells%edge_orientation(jc,jb,je)  / &
         &   ptr_patch%cells%area(jc,jb)!cell_area

        ENDDO !edge loop

      ENDDO !idx loop

    END DO !block loop
!$OMP END DO

    ! b) Geometrical factor for rotation
    rl_start = 1  ! #slo# changed to 1 - 2010-12-07
    rl_end = min_rlvert_int

    ! Vorticity should have the right sign
      ifac = 0
      SELECT CASE (i_cell_type)
      CASE (3)
        ifac = 1
      CASE (6)
        ifac = -1
      END SELECT
      ! values for the blocking
      i_startblk = ptr_patch%verts%start_blk(rl_start,1)
      i_endblk   = ptr_patch%verts%end_blk(rl_end,i_nchdom)
      !
      ! loop through all patch cells (and blocks)
      !
!$OMP DO PRIVATE(jb,je,jv,i_startidx,i_endidx,ile,ibe)
      DO jb = i_startblk, i_endblk

        CALL get_indices_v(ptr_patch, jb, i_startblk, i_endblk, &
          &                i_startidx, i_endidx, rl_start, rl_end)

        DO je = 1, 9-i_cell_type
          DO jv = i_startidx, i_endidx

            IF (je > ptr_patch%verts%num_edges(jv,jb)) CYCLE   ! relevant for hexagons

            ile = ptr_patch%verts%edge_idx(jv,jb,je)
            ibe = ptr_patch%verts%edge_blk(jv,jb,je)

            ptr_intp%geofac_rot(jv,je,jb) =              &
         &    ptr_patch%edges%dual_edge_length(ile,ibe) * &
         &    ptr_patch%verts%edge_orientation(jv,jb,je)/ &
         &    ptr_patch%verts%dual_area(jv,jb) * REAL(ifac,wp)

          ENDDO !vertex loop
        ENDDO

      END DO !block loop
!$OMP END DO

      ! c) Geometrical factor for nabla2_scalar
      rl_start = 1  ! #slo# changed to 1 - 2010-12-07
      rl_end = min_rlcell_int

      ! values for the blocking
      i_startblk = ptr_patch%cells%start_blk(rl_start,1)
      i_endblk   = ptr_patch%cells%end_blk(rl_end,i_nchdom)
      !
      ! loop through all patch cells (and blocks)
      !
!$OMP DO PRIVATE(jb,je,jc,ic,i_startidx,i_endidx,ile,ibe,ilc1,ibc1,&
!$OMP    ilc2,ibc2,ilnc,ibnc)
      DO jb = i_startblk, i_endblk

        CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
                           i_startidx, i_endidx, rl_start, rl_end)

        DO je = 1, i_cell_type
          DO jc = i_startidx, i_endidx

            ile = ptr_patch%cells%edge_idx(jc,jb,je)
            ibe = ptr_patch%cells%edge_blk(jc,jb,je)

            ilc1 = ptr_patch%edges%cell_idx(ile,ibe,1)
            ibc1 = ptr_patch%edges%cell_blk(ile,ibe,1)
            ilc2 = ptr_patch%edges%cell_idx(ile,ibe,2)
            ibc2 = ptr_patch%edges%cell_blk(ile,ibe,2)

            IF (jc == ilc1 .AND. jb == ibc1) THEN
              IF (i_cell_type == 3) THEN
                ptr_intp%geofac_n2s(jc,1,jb)     =  &
                &  ptr_intp%geofac_n2s(jc,1,jb)  -  &
                &  ptr_intp%geofac_div(jc,je,jb) /  &
                &  ptr_patch%edges%dual_edge_length(ile,ibe)
            ELSE IF (i_cell_type == 6) THEN
              ptr_intp%geofac_n2s(jc,1,jb)       =  &
                &  ptr_intp%geofac_n2s(jc,1,jb)  -  &
                &  ptr_intp%geofac_div(jc,je,jb) /  &
                &  ptr_patch%edges%dual_edge_length(ile,ibe)*  &
                &  ptr_patch%edges%system_orientation(ile,ibe)
            ENDIF
          ELSE IF (jc == ilc2 .AND. jb == ibc2) THEN
            IF (i_cell_type == 3) THEN
              ptr_intp%geofac_n2s(jc,1,jb)       =  &
                &  ptr_intp%geofac_n2s(jc,1,jb)  +  &
                &  ptr_intp%geofac_div(jc,je,jb) /  &
                &  ptr_patch%edges%dual_edge_length(ile,ibe)
            ELSE IF (i_cell_type == 6) THEN
              ptr_intp%geofac_n2s(jc,1,jb)       =  &
                &  ptr_intp%geofac_n2s(jc,1,jb)  +  &
                &  ptr_intp%geofac_div(jc,je,jb) /  &
                &  ptr_patch%edges%dual_edge_length(ile,ibe)*  &
                &  ptr_patch%edges%system_orientation(ile,ibe)
            ENDIF
          ENDIF
          DO ic = 1, i_cell_type
            ilnc = ptr_patch%cells%neighbor_idx(jc,jb,ic)
            ibnc = ptr_patch%cells%neighbor_blk(jc,jb,ic)
            IF (ilnc == ilc1 .AND. ibnc == ibc1) THEN
              IF (i_cell_type == 3) THEN
                ptr_intp%geofac_n2s(jc,ic+1,jb)     = &
                  &  ptr_intp%geofac_n2s(jc,ic+1,jb)- &
                  &  ptr_intp%geofac_div(jc,je,jb)  / &
                  &  ptr_patch%edges%dual_edge_length(ile,ibe)
              ELSE IF (i_cell_type == 6) THEN
                ptr_intp%geofac_n2s(jc,ic+1,jb)     = &
                  &  ptr_intp%geofac_n2s(jc,ic+1,jb)- &
                  &  ptr_intp%geofac_div(jc,je,jb)  / &
                  &  ptr_patch%edges%dual_edge_length(ile,ibe) * &
                  &  ptr_patch%edges%system_orientation(ile,ibe)
              ENDIF
            ELSE IF (ilnc == ilc2 .AND. ibnc == ibc2) THEN
              IF (i_cell_type == 3) THEN
                ptr_intp%geofac_n2s(jc,ic+1,jb)     = &
                  &  ptr_intp%geofac_n2s(jc,ic+1,jb)+ &
                  &  ptr_intp%geofac_div(jc,je,jb)  / &
                  &  ptr_patch%edges%dual_edge_length(ile,ibe)
              ELSE IF (i_cell_type == 6) THEN
                ptr_intp%geofac_n2s(jc,ic+1,jb)     = &
                  &  ptr_intp%geofac_n2s(jc,ic+1,jb)+ &
                  &  ptr_intp%geofac_div(jc,je,jb)  / &
                  &  ptr_patch%edges%dual_edge_length(ile,ibe) * &
                  &  ptr_patch%edges%system_orientation(ile,ibe)
              ENDIF
            ENDIF
          ENDDO

          ! To ensure that dummy edges have a factor of 0:
          IF (je > ptr_patch%cells%num_edges(jc,jb)) THEN
            ptr_intp%geofac_n2s(jc,je+1,jb) = 0._wp
          ENDIF

        ENDDO !cell loop
      ENDDO

    END DO !block loop
!$OMP END DO

    ! d) Geometrical factor for quad-cell divergence (triangles only)
    IF (i_cell_type == 3) THEN

      rl_start = 1  ! #slo# changed to 1 - 2010-12-07
      rl_end = min_rledge_int

      ! values for the blocking
      i_startblk = ptr_patch%edges%start_blk(rl_start,1)
      i_endblk   = ptr_patch%edges%end_blk(rl_end,i_nchdom)

!$OMP DO PRIVATE(jb,je,je1,i_startidx,i_endidx,ile,ibe)
      DO jb = i_startblk, i_endblk

        CALL get_indices_e(ptr_patch, jb, i_startblk, i_endblk, &
                           i_startidx, i_endidx, rl_start, rl_end)

        DO je1 = 1, 4
          DO je = i_startidx, i_endidx

          ile = ptr_patch%edges%quad_idx(je,jb,je1)
          ibe = ptr_patch%edges%quad_blk(je,jb,je1)

          ptr_intp%geofac_qdiv(je,je1,jb) =               &
            ptr_patch%edges%primal_edge_length(ile,ibe) *  &
            ptr_patch%edges%quad_orientation(je,jb,je1) /  &
            ptr_patch%edges%quad_area(je,jb)

          ENDDO !edge loop
        ENDDO

      END DO !block loop
!$OMP END DO

    ENDIF

    ! f) compute inverse dual edge length (used in math_operators for the ocean)

    rl_start = 1  ! #slo# changed to 1 - 2010-12-07
    rl_end = min_rledge_int

    ! Second step: computed projected orientation vectors and related information
    i_startblk = ptr_patch%edges%start_blk(rl_start,1)
    i_endblk   = ptr_patch%edges%end_blk(rl_end,i_nchdom)
    !
    ! loop through all patch edges
    !
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,je)
    DO jb = i_startblk, i_endblk

      CALL get_indices_e(ptr_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, rl_start, rl_end)

      DO je =  i_startidx, i_endidx

        ! compute inverse dual edge length (undefined for refin_ctrl=1)

        ptr_patch%edges%inv_dual_edge_length(je,jb) = &
          1._wp/ptr_patch%edges%dual_edge_length(je,jb)

      ENDDO

    END DO !block loop
!$OMP END DO

!$OMP END PARALLEL

    ! synchronize all elements of ptr_intp:

    DO ie = 1, i_cell_type

        z_sync_c(:,:) = ptr_intp%geofac_div(:,ie,:)
        CALL sync_patch_array(SYNC_C, ptr_patch, z_sync_c(:,:))
        ptr_intp%geofac_div(:,ie,:) = z_sync_c(:,:)

        z_sync_c(:,:) = ptr_intp%geofac_n2s(:,ie,:)
        CALL sync_patch_array(SYNC_C, ptr_patch, z_sync_c(:,:))
        ptr_intp%geofac_n2s(:,ie,:) = z_sync_c(:,:)

    END DO

    DO ie = 1, 4

        z_sync_e(:,:) = ptr_intp%geofac_qdiv(:,ie,:)
        CALL sync_patch_array(SYNC_e, ptr_patch, z_sync_e(:,:))
        ptr_intp%geofac_qdiv(:,ie,:) = z_sync_e(:,:)

    END DO

    CALL sync_patch_array(SYNC_E, ptr_patch, ptr_patch%edges%inv_dual_edge_length(:,:))

    DO ie = 1, 9-i_cell_type

        z_sync_v(:,:) = ptr_intp%geofac_rot(:,ie,:)
        CALL sync_patch_array(SYNC_v, ptr_patch, z_sync_v(:,:))
        ptr_intp%geofac_rot(:,ie,:) = z_sync_v(:,:)

    END DO

    CALL message (TRIM(routine), 'end')

  END SUBROUTINE init_geo_factors_oce
!-------------------------------------------------------------------------  



END MODULE mo_intp_coeffs
