
#ifdef __xlC__
! @PROCESS nosmp
! @PROCESS NOOPTimize
! @PROCESS smp=noopt
@process noopt
#endif
!#ifdef __PGI
! !pgi$g opt=1
!#endif
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
!!  Modification by Almut Gassmann, MPI-M (2012-04-19)
!!  - added routine init_tplane_c, which projects vertices and mass points onto
!!    a plane tangent to cell centers.
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

MODULE mo_intp_coeffs
  !-------------------------------------------------------------------------
  !
  !    ProTeX FORTRAN source: Style 2
  !    modified for ICON project, DWD/MPI-M 2006
  !
  !-------------------------------------------------------------------------
  !
  USE mo_kind,                ONLY: wp
  USE mo_math_constants,      ONLY: pi2, pi_2
  USE mo_exception,           ONLY: message, finish
  USE mo_impl_constants,      ONLY: min_rlcell, min_rledge, min_rlvert
  USE mo_impl_constants_grf,  ONLY: grf_nudge_start_c, grf_nudge_start_e
  USE mo_model_domain,        ONLY: t_patch, t_grid_edges, t_grid_vertices, t_grid_cells
  USE mo_grid_config,         ONLY: lplane, lfeedback, grid_sphere_radius
  USE mo_math_types,          ONLY: t_cartesian_coordinates, t_geographical_coordinates
  USE mo_math_utilities,      ONLY: gc2cc, cc2gc, gnomonic_proj,               &
    & gvec2cvec, cvec2gvec,                      &
    & rotate_latlon, arc_length,                 &
    & plane_torus_closest_coordinates
  USE mo_dynamics_config,     ONLY: divavg_cntrwgt
  USE mo_parallel_config,     ONLY: nproma
  USE mo_loopindices,         ONLY: get_indices_c, get_indices_e, get_indices_v
  USE mo_sync,                ONLY: sync_c, sync_e, sync_v, sync_patch_array, sync_idx, global_max
  USE mo_intp_data_strc,      ONLY: t_int_state
  USE mo_interpol_config,     ONLY: i_cori_method, nudge_zone_width, nudge_max_coeff, &
    & nudge_efold_width

  USE mo_grid_subset,         ONLY: get_index_range
  USE mo_grid_geometry_info,  ONLY: planar_torus_geometry, sphere_geometry
  USE mo_grid_subset,         ONLY: get_index_range
  USE mo_fortran_tools,       ONLY: init

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: compute_heli_bra_coeff_idx
  PUBLIC :: init_cellavg_wgt
  PUBLIC :: init_geo_factors
  PUBLIC :: complete_patchinfo
  PUBLIC :: init_tplane_e
  PUBLIC :: init_tplane_c
  PUBLIC :: init_nudgecoeffs
  PUBLIC :: tri_quadrature_pts

  ! flags for computing ocean coefficients
  LOGICAL, PARAMETER :: mid_point_dual_edge = .TRUE. !Please do not change this unless
  !you are sure, you know what you do.
  LOGICAL, PARAMETER :: larc_length = .FALSE.

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
    INTEGER :: jb, je, jc, nlen, nblks_e, npromz_e, nblks_c, npromz_c, &
      & n_edges, imyself, ineigh, incr, ile1_v(2), ibe1_v(2),&
      & je1, ile1, ibe1, je2, ile2, ibe2, je3, ile3, ibe3, je4, je5, &
      & ile5, ibe5, jc1, ilc1, ibc1, jc3, ilc3, ibc3, jc4, jv2, jv3, &
      & ilv3, ibv3, jv4,ilc,ibc, npromz_v, nblks_v
    REAL(wp)   :: z_tan_proj, z_nor_proj, z_metric, z_metric_new, z_norm, &
      & z_metric_a, z_metric_b, z_metric_c, z_metric_d
    REAL(wp),ALLOCATABLE :: z_frac_area(:,:,:), &
      & z_quad_north(:,:,:), z_quad_east(:,:,:)
    TYPE(t_grid_edges),    POINTER :: p_ed
    TYPE(t_grid_cells),    POINTER :: p_ce
    TYPE(t_grid_vertices), POINTER :: p_ve
    TYPE(t_patch),         POINTER :: p_pa
    TYPE(t_int_state),     POINTER :: p_in
    LOGICAL :: l_cycle, l_found
    INTEGER :: ile_c(6), ibe_c(6), n_verts, jv, ilv, ibv, jm, ilv1,ibv1, &
      & ineigh_c, ineigh_v, ilv_c(6), ibv_c(6), incrc,            &
      & jvindex, jeindex, ile, ibe
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
      nblks_c  = p_pa%nblks_c
      npromz_c = p_pa%npromz_c
      DO jb = 1, nblks_c
        IF (jb /= nblks_c) THEN
          nlen = nproma
        ELSE
          nlen = npromz_c
        ENDIF
        DO jc = 1, nlen

          IF(.NOT.p_ce%decomp_info%owner_mask(jc,jb)) CYCLE

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
                    & jb == p_ed%cell_blk(ile1,ibe1,1)) THEN
                    ineigh_c = 1
                  ELSE
                    ineigh_c = 2
                  ENDIF
                  IF(ilv == p_ed%vertex_idx(ile1,ibe1,1) .AND. &
                    & ibv == p_ed%vertex_blk(ile1,ibe1,1)) THEN
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
        CALL sync_patch_array(sync_c,ptr_patch,z_frac_area(jv,:,:))
      ENDDO

      ! heli coeffs
      nblks_e  = p_pa%nblks_e
      npromz_e = p_pa%npromz_e
      DO jb = 1, nblks_e
        IF (jb /= nblks_e) THEN
          nlen = nproma
        ELSE
          nlen = npromz_e
        ENDIF
        DO je = 1, nlen

          IF(.NOT.p_ed%decomp_info%owner_mask(je,jb)) CYCLE

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
              IF((ile_c(je1) == je) .AND. (ibe_c(je1)==jb) ) THEN
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
                    & p_ed%vertex_blk(ile1,ibe1,jv) ==ibv1) ) THEN
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
                IF(ile == ile1 .AND. ibe == ibe1) CYCLE ! that was the previous
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
                & *p_ed%primal_edge_length(ile_c(jeindex),ibe_c(jeindex)) &
                & /p_ed%dual_edge_length(je,jb)*zorient_ve

              ile1 = ile_c(jeindex)
              ibe1 = ibe_c(jeindex)

            ENDDO
            ! If there is a pentagon
            IF (jc==2.AND.incr==9) THEN
              incr = incr+1
              p_in%heli_vn_idx(incr,je,jb) = je
              p_in%heli_vn_blk(incr,je,jb) = jb
              p_in%heli_coeff (incr,je,jb) = 0.0_wp
            ENDIF
          ENDDO
        ENDDO
      ENDDO

      DO incr = 1, UBOUND(p_in%heli_coeff,1)
        CALL sync_patch_array(sync_e,ptr_patch,p_in%heli_coeff(incr,:,:))
        CALL sync_idx(sync_e,sync_e,ptr_patch,p_in%heli_vn_idx(incr,:,:),p_in%heli_vn_blk(incr,:,:))
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
    nblks_c  = p_pa%nblks_c
    npromz_c = p_pa%npromz_c
    DO jb = 1, nblks_c
      IF (jb /= nblks_c) THEN
        nlen = nproma
      ELSE
        nlen = npromz_c
      ENDIF
      DO jc = 1, nlen

        IF(.NOT.p_ce%decomp_info%owner_mask(jc,jb)) CYCLE

        IF (.NOT. lplane) THEN
          CALL gvec2cvec(0.0_wp,1.0_wp,p_ce%center(jc,jb)%lon,p_ce%center(jc,jb)%lat, &
            & z_cart_no%x(1),z_cart_no%x(2),z_cart_no%x(3))
          z_norm = SQRT( DOT_PRODUCT(z_cart_no%x(1:3),z_cart_no%x(1:3)))
          z_cart_no%x(1:3) = z_cart_no%x(1:3) /z_norm
          CALL gvec2cvec(1.0_wp,0.0_wp,p_ce%center(jc,jb)%lon,p_ce%center(jc,jb)%lat, &
            & z_cart_ea%x(1),z_cart_ea%x(2),z_cart_ea%x(3))
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
            & *p_ed%primal_edge_length(ile1,ibe1)     &
            & /p_ce%area(jc,jb)
          ! north projection
          p_in%hex_north(jc,je1,jb) = z_tan_proj*z_metric
          ! east projection
          p_in%hex_east(jc,je1,jb)  = z_nor_proj*z_metric
        ENDDO
      ENDDO
    ENDDO

    CALL sync_patch_array(sync_c,ptr_patch,p_in%hex_north)
    CALL sync_patch_array(sync_c,ptr_patch,p_in%hex_east)

    ! Vector reconstruction on triangles
    !===================================
    nblks_v  = p_pa%nblks_v
    npromz_v = p_pa%npromz_v
    DO jb = 1, nblks_v
      IF (jb /= nblks_v) THEN
        nlen = nproma
      ELSE
        nlen = npromz_v
      ENDIF
      DO jv = 1, nlen

        IF(.NOT.p_ve%decomp_info%owner_mask(jv,jb)) CYCLE

        IF (.NOT. lplane) THEN
          CALL gvec2cvec(0.0_wp,1.0_wp,p_ve%vertex(jv,jb)%lon,p_ve%vertex(jv,jb)%lat, &
            & z_cart_no%x(1),z_cart_no%x(2),z_cart_no%x(3))
          z_norm = SQRT( DOT_PRODUCT(z_cart_no%x(1:3),z_cart_no%x(1:3)))
          z_cart_no%x(1:3) = z_cart_no%x(1:3) /z_norm
          CALL gvec2cvec(1.0_wp,0.0_wp,p_ve%vertex(jv,jb)%lon,p_ve%vertex(jv,jb)%lat, &
            & z_cart_ea%x(1),z_cart_ea%x(2),z_cart_ea%x(3))
          z_norm = SQRT( DOT_PRODUCT(z_cart_ea%x(1:3),z_cart_ea%x(1:3)))
          z_cart_ea%x(1:3) = z_cart_ea%x(1:3) /z_norm
        ENDIF
        DO je1 = 1, 3
          ile1 = p_ve%edge_idx(jv,jb,je1)
          ibe1 = p_ve%edge_blk(jv,jb,je1)
          IF((p_ed%vertex_idx(ile1,ibe1,1)==jv) .AND. &
            & (p_ed%vertex_blk(ile1,ibe1,1)==jb)) THEN
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
            & *p_ed%edge_vert_length(ile1,ibe1,ineigh)&
            & /p_ve%dual_area(jv,jb)
          ! north projection
          p_in%tria_north(je1,jv,jb) = z_tan_proj*z_metric
          ! east projection
          p_in%tria_east(je1,jv,jb)  = z_nor_proj*z_metric
        ENDDO
      ENDDO
    ENDDO

    DO je1 = 1, 3
      CALL sync_patch_array(sync_v,ptr_patch,p_in%tria_north(je1,:,:))
      CALL sync_patch_array(sync_v,ptr_patch,p_in%tria_east(je1,:,:))
    ENDDO

    IF (i_cori_method /= 2 ) THEN

      ! Vector reconstruction on rhombi
      !================================
      nblks_e  = p_pa%nblks_e
      npromz_e = p_pa%npromz_e
      DO jb = 1, nblks_e
        IF (jb /= nblks_e) THEN
          nlen = nproma
        ELSE
          nlen = npromz_e
        ENDIF
        DO je = 1, nlen

          IF(.NOT.p_ed%decomp_info%owner_mask(je,jb)) CYCLE

          IF (.NOT. lplane) THEN
            CALL gvec2cvec(0.0_wp,1.0_wp,p_ed%center(je,jb)%lon,p_ed%center(je,jb)%lat, &
              & z_cart_no%x(1),z_cart_no%x(2),z_cart_no%x(3))
            z_norm = SQRT( DOT_PRODUCT(z_cart_no%x(1:3),z_cart_no%x(1:3)))
            z_cart_no%x(1:3) = z_cart_no%x(1:3) /z_norm
            CALL gvec2cvec(1.0_wp,0.0_wp,p_ed%center(je,jb)%lon,p_ed%center(je,jb)%lat, &
              & z_cart_ea%x(1),z_cart_ea%x(2),z_cart_ea%x(3))
            z_norm = SQRT( DOT_PRODUCT(z_cart_ea%x(1:3),z_cart_ea%x(1:3)))
            z_cart_ea%x(1:3) = z_cart_ea%x(1:3) /z_norm
          ENDIF
          DO je1 = 1,6
            IF (je1 <=4 )THEN
              ile1 = p_ed%quad_idx(je,jb,je1)
              ibe1 = p_ed%quad_blk(je,jb,je1)
              IF(((p_ed%vertex_idx(ile1,ibe1,1)==p_ed%vertex_idx(je,jb,1)) .AND. &
                & (p_ed%vertex_blk(ile1,ibe1,1)==p_ed%vertex_blk(je,jb,1))) .OR. &
                & ((p_ed%vertex_idx(ile1,ibe1,1)==p_ed%vertex_idx(je,jb,2)) .AND. &
                & (p_ed%vertex_blk(ile1,ibe1,1)==p_ed%vertex_blk(je,jb,2))))THEN
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
              & *p_ed%dual_edge_length(ile1,ibe1)         &
              & /p_ed%quad_area(je,jb)
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
          CALL sync_patch_array(sync_e,ptr_patch,p_in%quad_north(je1,:,:))
          CALL sync_patch_array(sync_e,ptr_patch,p_in%quad_east(je1,:,:))
        ENDDO
      ENDIF
      DO je1 = 1, 6
        CALL sync_patch_array(sync_e,ptr_patch,z_quad_north(je1,:,:))
        CALL sync_patch_array(sync_e,ptr_patch,z_quad_east (je1,:,:))
      ENDDO

      ! Computation of the coefficients for the vorticity flux term
      !============================================================
      nblks_e  = p_pa%nblks_e
      npromz_e = p_pa%npromz_e

      DO jb = 1, nblks_e
        IF (jb /= nblks_e) THEN
          nlen = nproma
        ELSE
          nlen = npromz_e
        ENDIF
        DO je = 1, nlen

          IF(.NOT.p_ed%decomp_info%owner_mask(je,jb)) CYCLE

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
                & *0.5_wp*p_ce%area(ilc1,ibc1)/p_ed%area_edge(je,jb)
              p_in%heli_coeff(2*jc1  ,je,jb) = -p_in%hex_north(ilc1,imyself,ibc1) &
                & *0.5_wp*p_ce%area(ilc1,ibc1)/p_ed%area_edge(je,jb)
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
                      & p_ed%quad_blk(ile3,ibe3,je4)==ibe2)THEN
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
                  & -p_in%hex_north(ilc1,imyself,ibc1)*p_in%hex_east (ilc1,je2,ibc1))

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
              & *p_ed%quad_area(je,jb)/p_ed%area_edge(je,jb)
            p_in%heli_coeff(6,je,jb) = -p_in%quad_north(5,je,jb)/6.0_wp &
              & *p_ed%quad_area(je,jb)/p_ed%area_edge(je,jb)
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
                & *p_ed%quad_area(ile1,ibe1)/p_ed%area_edge(je,jb)
              p_in%heli_coeff(2*je1+6,je,jb) = -p_in%quad_north(imyself,ile1,ibe1)/6.0_wp &
                & *p_ed%quad_area(ile1,ibe1)/p_ed%area_edge(je,jb)

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
                  & -z_quad_north(imyself,ile1,ibe1)*z_quad_east (je2,ile1,ibe1))
                z_metric_new = 0.5_wp*z_metric_a

                !=========================================
                ! Hexagonal part, this occurs twice per e1
                !=========================================
                DO jc3 = 1, 2
                  ilc3 = p_ed%cell_idx(ile2,ibe2,jc3)
                  ibc3 = p_ed%cell_blk(ile2,ibe2,jc3)
                  DO jc4 = 1, 2
                    IF( ilc3 == p_ed%cell_idx(je,jb,jc4) .AND. &
                      & ibc3 == p_ed%cell_blk(je,jb,jc4) ) THEN

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
                        & -p_in%hex_north(ilc3,imyself,ibc3)*p_in%hex_east (ilc3,ineigh,ibc3))
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
                      & ibv3==p_ed%vertex_blk(ile2,ibe2,jv4)) THEN

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
                        & *(z_quad_north(5,ile2,ibe2)+z_quad_north(6,ile2,ibe2))&
                        & -z_quad_north(imyself,ile2,ibe2) &
                        & *(z_quad_east (5,ile2,ibe2)+z_quad_east (6,ile2,ibe2)))
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
                        & *z_quad_north(ineigh,je,jb)  &
                        & -(z_quad_north(5,je,jb)+z_quad_north(6,je,jb)) &
                        & *z_quad_east (ineigh,je,jb))
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
        CALL sync_patch_array(sync_e,ptr_patch,p_in%heli_coeff(incr,:,:))
        IF(i_cori_method < 3) &
          & CALL sync_idx(sync_e,sync_e,ptr_patch,p_in%heli_vn_idx(incr,:,:),p_in%heli_vn_blk(incr,:,:))
      ENDDO

    ENDIF

    DEALLOCATE(z_quad_east)
    DEALLOCATE(z_quad_north)

  END SUBROUTINE compute_heli_bra_coeff_idx
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
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
    TYPE(t_patch), TARGET, INTENT(inout) :: ptr_patch
    TYPE(t_int_state), INTENT(inout):: ptr_int

    INTEGER                     :: max_iter !max no. of iterations for forcing
                                            !mass conservation
    CHARACTER(LEN=*), PARAMETER :: method_name = 'mo_intp_coeffs:init_cellavg_wgt'

    SELECT CASE(ptr_patch%geometry_info%geometry_type)

    CASE (planar_torus_geometry)
      max_iter = 1
      CALL calculate_uniform_bilinear_cellavg_wgt( ptr_patch, ptr_int )
      CALL force_mass_conservation_to_cellavg_wgt( ptr_patch, ptr_int, max_iter)

    CASE (sphere_geometry)
      max_iter = 1000
      CALL calculate_bilinear_cellavg_wgt( ptr_patch, ptr_int )
      CALL force_mass_conservation_to_cellavg_wgt( ptr_patch, ptr_int, max_iter )

    CASE default
      CALL finish(method_name, "Undefined geometry type")

    END SELECT

  END SUBROUTINE init_cellavg_wgt
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! Computes the weighting coefficients for cell averaging with
  !! variable interpolation factors. Results are stored in ptr_patch%cells%avg_wgt
  !!
  !! The weighting factors are based on the requirement that sum(w(i)*x(i)) = 0
  !! and sum(w(i)*y(i)) = 0, which ensures that linear horizontal gradients
  !! are not aliased into a checkerboard pattern between upward- and downward
  !! directed cells. The third condition is sum(w(i)) = 1., and the weight
  !! of the local point is 0.5 (see above).
  !!
  !! @par Revision History
  !!  developed by Anurag Dipankar, MPI-M(2012-12-12)
  !!  adapted by Leonidas Linardakis, MPI-M(2012-12-12)
  SUBROUTINE calculate_uniform_bilinear_cellavg_wgt( patch, interpolation_state )
    !  patch on which computation is performed
    TYPE(t_patch), TARGET, INTENT(inout) :: patch
    ! Interpolation state
    TYPE(t_int_state), INTENT(inout):: interpolation_state

    INTEGER :: cell_block, cell_index, start_index, end_index
    REAL(wp) :: local_weight, neigbor_weight

    !-----------------------------------------------------------------------
    ! Initial weighting factor of the local grid point
    local_weight = divavg_cntrwgt
    neigbor_weight = (1.0_wp - local_weight) / 3.0_wp

!$OMP PARALLEL
!$OMP DO PRIVATE(cell_block, cell_index, start_index, end_index) ICON_OMP_DEFAULT_SCHEDULE
    DO cell_block = patch%cells%all%start_block, patch%cells%all%end_block
      CALL get_index_range(patch%cells%all, cell_block, start_index, end_index)
      DO cell_index = start_index, end_index

        ! Simple for plane torus
        interpolation_state%c_bln_avg(cell_index,1,  cell_block) = local_weight
        interpolation_state%c_bln_avg(cell_index,2:4,cell_block) = neigbor_weight

      ENDDO !cell loop
  END DO !block loop
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  ! Note: no sync is required since the weights are calculated for all cells
  !       including halos
  END SUBROUTINE calculate_uniform_bilinear_cellavg_wgt
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! Computes the weighting coefficients for cell averaging with
  !! variable interpolation factors. Results are stored in ptr_patch%cells%avg_wgt
  !!
  !! The weighting factors are based on the requirement that sum(w(i)*x(i)) = 0
  !! and sum(w(i)*y(i)) = 0, which ensures that linear horizontal gradients
  !! are not aliased into a checkerboard pattern between upward- and downward
  !! directed cells. The third condition is sum(w(i)) = 1., and the weight
  !! of the local point is 0.5 (see above).
  !!
  !! @par Revision History
  !!  developed by Guenther Zaengl, 2008-12-05
  !!
  SUBROUTINE calculate_bilinear_cellavg_wgt( ptr_patch, ptr_int )
    !  patch on which computation is performed
    TYPE(t_patch), TARGET, INTENT(inout) :: ptr_patch
    ! Interpolation state
    TYPE(t_int_state), INTENT(inout):: ptr_int

    INTEGER :: jc, jb
    INTEGER :: rl_start, rl_end
    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx, i_nchdom

    INTEGER :: ilc1, ibc1, ilc2, ibc2, ilc3, ibc3

    REAL(wp) :: xtemp,ytemp,wgt(3),xloc,yloc,x(3),y(3), &
      & pollat,pollon,wgt_loc

    REAL(wp) :: cell_area   ! area of triangular cell made up by 3 neighboring cell centers
    REAL(wp) :: mfac        ! 1 or -1, depending on the sign of yloc
                            ! takes care of correct sign of B-matrix
    !-----------------------------------------------------------------------

    ! Initial weighting factor of the local grid point
    wgt_loc = divavg_cntrwgt

    ! values for the blocking
    rl_start = 2
    rl_end = min_rlcell

    i_nchdom   = MAX(1,ptr_patch%n_childdom)
    i_startblk = ptr_patch%cells%start_blk(rl_start,1)
    i_endblk   = ptr_patch%cells%end_blk(rl_end,i_nchdom)

    ! Compute coefficients for bilinear interpolation
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx,yloc,xloc,pollat,pollon,&
!$OMP            ilc1,ibc1,ilc2,ibc2,ilc3,ibc3,xtemp,ytemp,wgt,x,y,cell_area,mfac) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
        & i_startidx, i_endidx, rl_start, rl_end)

      DO jc = i_startidx, i_endidx

        IF(.NOT.ptr_patch%cells%decomp_info%owner_mask(jc,jb)) CYCLE

        yloc = ptr_patch%cells%center(jc,jb)%lat
        xloc = ptr_patch%cells%center(jc,jb)%lon

        ! Rotate local point into the equator for better accuracy of bilinear weights
        IF (yloc >= 0._wp) THEN
          pollat = yloc - pi2/4._wp
          mfac   = -1._wp
        ELSE
          pollat = yloc + pi2/4._wp
          mfac   = 1._wp
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

        IF (ABS(x(2)-x(1)) > 1.e-11_wp .AND. ABS(y(3)-y(1)) > 1.e-11_wp ) THEN
          wgt(3) = 1._wp/( (y(3)-y(1)) - (x(3)-x(1))*(y(2)-y(1))/(x(2)-x(1)) ) * &
            & (1._wp-wgt_loc)*( -y(1) + x(1)*(y(2)-y(1))/(x(2)-x(1)) )
          wgt(2) = (-(1._wp-wgt_loc)*x(1) - wgt(3)*(x(3)-x(1)))/(x(2)-x(1))
          wgt(1) = 1._wp - wgt_loc - wgt(2) - wgt(3)
        ELSE
          wgt(2) = 1._wp/( (y(2)-y(1)) - (x(2)-x(1))*(y(3)-y(1))/(x(3)-x(1)) ) * &
            & (1._wp-wgt_loc)*( -y(1) + x(1)*(y(3)-y(1))/(x(3)-x(1)) )
          wgt(3) = (-(1._wp-wgt_loc)*x(1) - wgt(2)*(x(2)-x(1)))/(x(3)-x(1))
          wgt(1) = 1._wp - wgt_loc - wgt(2) - wgt(3)
        ENDIF

        ! Store results in ptr_patch%cells%avg_wgt
        ptr_int%c_bln_avg(jc,1,jb) = wgt_loc
        ptr_int%c_bln_avg(jc,2,jb) = wgt(1)
        ptr_int%c_bln_avg(jc,3,jb) = wgt(2)
        ptr_int%c_bln_avg(jc,4,jb) = wgt(3)


        ! B-matrix for cell based gradient (based on linear triangular finite element)
        !
        ! !!Attention!!: pollat IF-statement flips sign of (xi-xj), (yi-yj) terms
        ! That's why we need to multiply by -1 for yloc >= 0._wp
        !
        cell_area = 0.5_wp*grid_sphere_radius**2  &
          &       *((x(2)*y(3)-x(3)*y(2)) - (x(1)*y(3)-x(3)*y(1)) + (x(1)*y(2)-x(2)*y(1)))
        ! compute b-matrix
        ptr_int%gradc_bmat(jc,1,1,jb) = mfac*grid_sphere_radius*(y(2) - y(3))/(2._wp*cell_area)
        ptr_int%gradc_bmat(jc,1,2,jb) = mfac*grid_sphere_radius*(y(3) - y(1))/(2._wp*cell_area)
        ptr_int%gradc_bmat(jc,1,3,jb) = mfac*grid_sphere_radius*(y(1) - y(2))/(2._wp*cell_area)

        ptr_int%gradc_bmat(jc,2,1,jb) = mfac*grid_sphere_radius*(x(3) - x(2))/(2._wp*cell_area)
        ptr_int%gradc_bmat(jc,2,2,jb) = mfac*grid_sphere_radius*(x(1) - x(3))/(2._wp*cell_area)
        ptr_int%gradc_bmat(jc,2,3,jb) = mfac*grid_sphere_radius*(x(2) - x(1))/(2._wp*cell_area)

      ENDDO !cell loop

    END DO !block loop
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    CALL sync_patch_array(sync_c,ptr_patch,ptr_int%c_bln_avg)
    CALL sync_patch_array(sync_c,ptr_patch,ptr_int%gradc_bmat(:,1,:,:))
    CALL sync_patch_array(sync_c,ptr_patch,ptr_int%gradc_bmat(:,2,:,:))

  END SUBROUTINE calculate_bilinear_cellavg_wgt
  !-----------------------------------------------------------------------


  !-----------------------------------------------------------------------
  !>
  ! The coefficients for bilinear interpolation are iteratively modified
  ! in order to obtain mass conservation.
  ! The criterion for conservation is that the three-point divergence
  ! calculated for any given grid point is used with a total factor of 1
  SUBROUTINE force_mass_conservation_to_cellavg_wgt( ptr_patch, ptr_int, niter )
    !  patch on which computation is performed
    TYPE(t_patch), TARGET, INTENT(inout) :: ptr_patch
    ! Interpolation state
    TYPE(t_int_state), INTENT(inout)     :: ptr_int
    ! max number of iterations
    INTEGER, INTENT(in)                  ::  niter

    INTEGER :: jc, je, jb
    INTEGER :: rl_start, rl_end
    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx, i_nchdom

    INTEGER :: ilc1, ibc1, ilc2, ibc2, ilc3, ibc3, inb1, inb2, inb3, ie4, ie5
    INTEGER :: ile1, ibe1, ile2, ibe2, ile3, ibe3, ile4, ibe4
    INTEGER, DIMENSION(nproma) :: iie1, iie2, iie3, iie4
    REAL(wp), DIMENSION (nproma,3) :: z_nx1, z_nx2, z_nx3, z_nx4, z_nx5
    REAL(wp) :: checksum(nproma,ptr_patch%nblks_e)

    REAL(wp) :: relax_coeff
    INTEGER ::  iter

    REAL(wp) :: maxwgt_loc,minwgt_loc

    REAL(wp), DIMENSION(nproma,ptr_patch%nblks_c)  :: wgt_loc_sum, resid

    INTEGER, DIMENSION(nproma,ptr_patch%nblks_c,3) :: inv_neighbor_id
    REAL(wp), DIMENSION(nproma,ptr_patch%nblks_c,3) :: z_inv_neighbor_id

#ifdef DEBUG_COEFF
    REAL(wp) :: sum1
#endif
    CHARACTER(LEN=*), PARAMETER :: method_name = 'mo_intp_coeffs:force_mass_conservation_to_cellavg_wgt'

    ! write(0,*) "force_mass_conservation_to_cellavg_wgt, processing ", ptr_patch%grid_filename
    !-----------------------------------------------------------------------
    ! The coefficients for bilinear interpolation are now iteratively modified
    ! in order to obtain mass conservation.
    ! The criterion for conservation is that the three-point divergence
    ! calculated for any given grid point is used with a total factor of 1
    ! Number of iterations for computation of bilinear weights
    i_nchdom   = MAX(1,ptr_patch%n_childdom)

    !niter = 1000 !now this value is passed as an arguement

    ! Relaxation coefficient for adaptation of local weight (empirically determined)
    relax_coeff = 0.46_wp

    ! Maximum/minimum  weighting factors of the local grid point
    maxwgt_loc = divavg_cntrwgt + 0.003_wp
    minwgt_loc = divavg_cntrwgt - 0.003_wp

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
!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx,ilc1,ibc1,ilc2,ibc2,ilc3,&
!$OMP ibc3) ICON_OMP_DEFAULT_SCHEDULE
!    DO jb = i_startblk, i_endblk
!
!      CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
!        & i_startidx, i_endidx, rl_start, rl_end)
!
!      DO jc = i_startidx, i_endidx
     DO jb = ptr_patch%cells%in_domain%start_block, ptr_patch%cells%in_domain%end_block
       CALL get_index_range(ptr_patch%cells%in_domain, jb, i_startidx, i_endidx)
       DO jc = i_startidx, i_endidx

        IF(.NOT.ptr_patch%cells%decomp_info%owner_mask(jc,jb)) CYCLE

        ! line and block indices of the neighbouring cells

        ilc1 = ptr_patch%cells%neighbor_idx(jc,jb,1)
        ibc1 = ptr_patch%cells%neighbor_blk(jc,jb,1)
        ilc2 = ptr_patch%cells%neighbor_idx(jc,jb,2)
        ibc2 = ptr_patch%cells%neighbor_blk(jc,jb,2)
        ilc3 = ptr_patch%cells%neighbor_idx(jc,jb,3)
        ibc3 = ptr_patch%cells%neighbor_blk(jc,jb,3)

        IF ( (ilc1>0) .AND. (ibc1>0) ) THEN
          IF ((ptr_patch%cells%neighbor_idx(ilc1,ibc1,1) == jc) .AND. &
            & (ptr_patch%cells%neighbor_blk(ilc1,ibc1,1) == jb)) THEN
            inv_neighbor_id(jc,jb,1) = 1
          ELSE IF ((ptr_patch%cells%neighbor_idx(ilc1,ibc1,2) == jc) .AND. &
            & (ptr_patch%cells%neighbor_blk(ilc1,ibc1,2) == jb)) THEN
            inv_neighbor_id(jc,jb,1) = 2
          ELSE IF ((ptr_patch%cells%neighbor_idx(ilc1,ibc1,3) == jc) .AND. &
            & (ptr_patch%cells%neighbor_blk(ilc1,ibc1,3) == jb)) THEN
            inv_neighbor_id(jc,jb,1) = 3
          ELSE
            CALL finish(method_name, "Undefined inv_neighbor_id 1")
          ENDIF
        ENDIF
        IF ( (ilc2>0) .AND. (ibc2>0) ) THEN
          IF ((ptr_patch%cells%neighbor_idx(ilc2,ibc2,1) == jc) .AND. &
            & (ptr_patch%cells%neighbor_blk(ilc2,ibc2,1) == jb)) THEN
            inv_neighbor_id(jc,jb,2)  = 1
          ELSE IF ((ptr_patch%cells%neighbor_idx(ilc2,ibc2,2) == jc) .AND. &
            & (ptr_patch%cells%neighbor_blk(ilc2,ibc2,2) == jb)) THEN
            inv_neighbor_id(jc,jb,2)  = 2
          ELSE IF ((ptr_patch%cells%neighbor_idx(ilc2,ibc2,3) == jc) .AND. &
            & (ptr_patch%cells%neighbor_blk(ilc2,ibc2,3) == jb)) THEN
            inv_neighbor_id(jc,jb,2)  = 3
          ELSE
            CALL finish(method_name, "Undefined inv_neighbor_id 2")
          ENDIF
        ENDIF
        IF ( (ilc3>0) .AND. (ibc3>0) ) THEN
          IF ((ptr_patch%cells%neighbor_idx(ilc3,ibc3,1) == jc) .AND. &
            & (ptr_patch%cells%neighbor_blk(ilc3,ibc3,1) == jb)) THEN
            inv_neighbor_id(jc,jb,3)  = 1
          ELSE IF ((ptr_patch%cells%neighbor_idx(ilc3,ibc3,2) == jc) .AND. &
            & (ptr_patch%cells%neighbor_blk(ilc3,ibc3,2) == jb)) THEN
            inv_neighbor_id(jc,jb,3)  = 2
          ELSE IF ((ptr_patch%cells%neighbor_idx(ilc3,ibc3,3) == jc) .AND. &
            & (ptr_patch%cells%neighbor_blk(ilc3,ibc3,3) == jb)) THEN
            inv_neighbor_id(jc,jb,3)  = 3
          ELSE
            CALL finish(method_name, "Undefined inv_neighbor_id 3")
          ENDIF
        ENDIF

      ENDDO !cell loop

    END DO !block loop
!$OMP END DO NOWAIT
!$OMP END PARALLEL
    z_inv_neighbor_id = REAL(inv_neighbor_id,wp)
    CALL sync_patch_array(sync_c,ptr_patch,z_inv_neighbor_id(:,:,1))
    CALL sync_patch_array(sync_c,ptr_patch,z_inv_neighbor_id(:,:,2))
    CALL sync_patch_array(sync_c,ptr_patch,z_inv_neighbor_id(:,:,3))
    inv_neighbor_id = NINT(z_inv_neighbor_id)

    DO iter = 1, niter

      ! Compute sum of weighting coefficients with which
      ! each local divergence value is used
      ! Note: the summation needs to be split into 4 loops in order to
      ! allow for vectorization and parallelization
      wgt_loc_sum = 0._wp

      rl_start = 2
      rl_end = min_rlcell
      i_startblk = ptr_patch%cells%start_blk(rl_start,1)
      i_endblk   = ptr_patch%cells%end_blk(rl_end,i_nchdom)


!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx,ilc1,ibc1,ilc2,ibc2,ilc3,ibc3,inb1,&
!$OMP inb2,inb3) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk, i_endblk

        CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
          & i_startidx, i_endidx, rl_start, rl_end)

        DO jc = i_startidx, i_endidx

          IF(.NOT. ptr_patch%cells%decomp_info%owner_mask(jc,jb)) CYCLE

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
            & ptr_int%c_bln_avg(jc,1,jb)*ptr_patch%cells%area(jc,jb)          + &
            & ptr_int%c_bln_avg(ilc1,inb1,ibc1)*ptr_patch%cells%area(ilc1,ibc1)  + &
            & ptr_int%c_bln_avg(ilc2,inb2,ibc2)*ptr_patch%cells%area(ilc2,ibc2)  + &
            & ptr_int%c_bln_avg(ilc3,inb3,ibc3)*ptr_patch%cells%area(ilc3,ibc3)

        ENDDO !cell loop

      END DO !block loop
!$OMP END DO NOWAIT
!$OMP END PARALLEL

      rl_start = 3
      i_startblk = ptr_patch%cells%start_blk(rl_start,1)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk, i_endblk

        CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
          & i_startidx, i_endidx, rl_start, rl_end)

        DO jc = i_startidx, i_endidx

          IF(.NOT. ptr_patch%cells%decomp_info%owner_mask(jc,jb)) CYCLE

          ! For mass conservation, wgt_loc_sum/area should be 1 for each cell
          ! The deviation therefrom is termed residuum here.

          resid(jc,jb) = wgt_loc_sum(jc,jb)/ptr_patch%cells%area(jc,jb)-1._wp

        ENDDO !cell loop

      END DO !block loop
!$OMP END DO NOWAIT
!$OMP END PARALLEL

      CALL sync_patch_array(sync_c,ptr_patch,resid)

      IF (iter < niter) THEN ! Apply iterative correction to weighting coefficients
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx,ilc1,ibc1,ilc2,ibc2,&
!$OMP ilc3,ibc3) ICON_OMP_DEFAULT_SCHEDULE
        DO jb = i_startblk, i_endblk

          CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
            & i_startidx, i_endidx, rl_start, rl_end)

          DO jc = i_startidx, i_endidx

            IF(.NOT. ptr_patch%cells%decomp_info%owner_mask(jc,jb)) CYCLE

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
!$OMP END DO NOWAIT
!$OMP END PARALLEL

        CALL sync_patch_array(sync_c,ptr_patch,ptr_int%c_bln_avg)

      ELSE ! In the last iteration, enforce the mass conservation condition
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
        DO jb = i_startblk, i_endblk

          CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
            & i_startidx, i_endidx, rl_start, rl_end)

          DO jc = i_startidx, i_endidx

            IF(.NOT. ptr_patch%cells%decomp_info%owner_mask(jc,jb)) CYCLE

            ! Modify weighting coefficients

            ptr_int%c_bln_avg(jc,1,jb) = ptr_int%c_bln_avg(jc,1,jb) - resid(jc,jb)

          ENDDO !cell loop

        END DO !block loop
!$OMP END DO NOWAIT
!$OMP END PARALLEL

        CALL sync_patch_array(sync_c,ptr_patch,ptr_int%c_bln_avg)

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
!$OMP DO PRIVATE(jb,je,i_startidx,i_endidx,ilc1,ibc1,ilc2,ibc2,inb1,&
!$OMP inb2,inb3,ie4,ie5) ICON_OMP_DEFAULT_SCHEDULE
        DO jb = i_startblk, i_endblk

          CALL get_indices_e(ptr_patch, jb, i_startblk, i_endblk, &
            & i_startidx, i_endidx, rl_start, rl_end)

          DO je = i_startidx, i_endidx

            IF(.NOT. ptr_patch%edges%decomp_info%owner_mask(je,jb)) CYCLE

            ilc1 = ptr_patch%edges%cell_idx(je,jb,1)
            ibc1 = ptr_patch%edges%cell_blk(je,jb,1)
            ilc2 = ptr_patch%edges%cell_idx(je,jb,2)
            ibc2 = ptr_patch%edges%cell_blk(je,jb,2)

            IF (je == ptr_patch%cells%edge_idx(ilc1,ibc1,1) .AND. &
              & jb == ptr_patch%cells%edge_blk(ilc1,ibc1,1)) THEN

              inb1 = inv_neighbor_id(ilc1,ibc1,1)
              ie4  = MOD(inb1,  3)+1
              ie5  = MOD(inb1+1,3)+1

              ptr_int%e_flx_avg(je,2,jb) = ptr_int%c_bln_avg(ilc2,inb1+1,ibc2) &
                & *ptr_int%geofac_div(ilc1,2,ibc1)/ptr_int%geofac_div(ilc2,inb1,ibc2)

              ptr_int%e_flx_avg(je,3,jb) = ptr_int%c_bln_avg(ilc2,inb1+1,ibc2) &
                & *ptr_int%geofac_div(ilc1,3,ibc1)/ptr_int%geofac_div(ilc2,inb1,ibc2)

              ptr_int%e_flx_avg(je,4,jb) = ptr_int%c_bln_avg(ilc1,2,ibc1) &
                & *ptr_int%geofac_div(ilc2,ie4,ibc2)/ptr_int%geofac_div(ilc1,1,ibc1)

              ptr_int%e_flx_avg(je,5,jb) = ptr_int%c_bln_avg(ilc1,2,ibc1) &
                & *ptr_int%geofac_div(ilc2,ie5,ibc2)/ptr_int%geofac_div(ilc1,1,ibc1)


            ELSE IF (je == ptr_patch%cells%edge_idx(ilc1,ibc1,2) .AND. &
              & jb == ptr_patch%cells%edge_blk(ilc1,ibc1,2)) THEN

              inb2 = inv_neighbor_id(ilc1,ibc1,2)
              ie4  = MOD(inb2  ,3)+1
              ie5  = MOD(inb2+1,3)+1

              ptr_int%e_flx_avg(je,2,jb) = ptr_int%c_bln_avg(ilc2,inb2+1,ibc2) &
                & *ptr_int%geofac_div(ilc1,3,ibc1)/ptr_int%geofac_div(ilc2,inb2,ibc2)

              ptr_int%e_flx_avg(je,3,jb) = ptr_int%c_bln_avg(ilc2,inb2+1,ibc2) &
                & *ptr_int%geofac_div(ilc1,1,ibc1)/ptr_int%geofac_div(ilc2,inb2,ibc2)

              ptr_int%e_flx_avg(je,4,jb) = ptr_int%c_bln_avg(ilc1,3,ibc1) &
                & *ptr_int%geofac_div(ilc2,ie4,ibc2)/ptr_int%geofac_div(ilc1,2,ibc1)

              ptr_int%e_flx_avg(je,5,jb) = ptr_int%c_bln_avg(ilc1,3,ibc1) &
                & *ptr_int%geofac_div(ilc2,ie5,ibc2)/ptr_int%geofac_div(ilc1,2,ibc1)


            ELSE IF (je == ptr_patch%cells%edge_idx(ilc1,ibc1,3) .AND. &
              & jb == ptr_patch%cells%edge_blk(ilc1,ibc1,3)) THEN

              inb3 = inv_neighbor_id(ilc1,ibc1,3)
              ie4  = MOD(inb3  ,3)+1
              ie5  = MOD(inb3+1,3)+1

              ptr_int%e_flx_avg(je,2,jb) = ptr_int%c_bln_avg(ilc2,inb3+1,ibc2) &
                & *ptr_int%geofac_div(ilc1,1,ibc1)/ptr_int%geofac_div(ilc2,inb3,ibc2)

              ptr_int%e_flx_avg(je,3,jb) = ptr_int%c_bln_avg(ilc2,inb3+1,ibc2) &
                & *ptr_int%geofac_div(ilc1,2,ibc1)/ptr_int%geofac_div(ilc2,inb3,ibc2)

              ptr_int%e_flx_avg(je,4,jb) = ptr_int%c_bln_avg(ilc1,4,ibc1) &
                & *ptr_int%geofac_div(ilc2,ie4,ibc2)/ptr_int%geofac_div(ilc1,3,ibc1)

              ptr_int%e_flx_avg(je,5,jb) = ptr_int%c_bln_avg(ilc1,4,ibc1) &
                & *ptr_int%geofac_div(ilc2,ie5,ibc2)/ptr_int%geofac_div(ilc1,3,ibc1)

            ENDIF

          ENDDO !edge loop

        END DO !block loop
!$OMP END DO NOWAIT
!$OMP END PARALLEL

        CALL sync_patch_array(sync_e,ptr_patch,ptr_int%e_flx_avg)

        rl_start = 5
        i_startblk = ptr_patch%edges%start_blk(rl_start,1)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,je,i_startidx,i_endidx,ilc1,ibc1,ilc2,ibc2,inb1,inb2,inb3,ie4,ie5, &
!$OMP            ile1,ibe1,ile2,ibe2,ile3,ibe3,ile4,ibe4,iie1,iie2,iie3,iie4) ICON_OMP_DEFAULT_SCHEDULE
        DO jb = i_startblk, i_endblk

          CALL get_indices_e(ptr_patch, jb, i_startblk, i_endblk, &
            & i_startidx, i_endidx, rl_start, rl_end)

          DO je = i_startidx, i_endidx

            IF(.NOT. ptr_patch%edges%decomp_info%owner_mask(je,jb)) CYCLE

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
              & ptr_patch%edges%cell_blk(ile1,ibe1,1) == ibc1 ) THEN
              iie1(je) = 3
            ELSE IF (ptr_patch%edges%cell_idx(ile1,ibe1,2) == ilc1 .AND. &
              & ptr_patch%edges%cell_blk(ile1,ibe1,2) == ibc1 ) THEN
              iie1(je) = 5
            ENDIF
            IF (ptr_patch%edges%cell_idx(ile2,ibe2,1) == ilc1 .AND. &
              & ptr_patch%edges%cell_blk(ile2,ibe2,1) == ibc1 ) THEN
              iie2(je) = 2
            ELSE IF (ptr_patch%edges%cell_idx(ile2,ibe2,2) == ilc1 .AND. &
              & ptr_patch%edges%cell_blk(ile2,ibe2,2) == ibc1 ) THEN
              iie2(je) = 4
            ENDIF
            IF (ptr_patch%edges%cell_idx(ile3,ibe3,1) == ilc2 .AND. &
              & ptr_patch%edges%cell_blk(ile3,ibe3,1) == ibc2 ) THEN
              iie3(je) = 3
            ELSE IF (ptr_patch%edges%cell_idx(ile3,ibe3,2) == ilc2 .AND. &
              & ptr_patch%edges%cell_blk(ile3,ibe3,2) == ibc2 ) THEN
              iie3(je) = 5
            ENDIF
            IF (ptr_patch%edges%cell_idx(ile4,ibe4,1) == ilc2 .AND. &
              & ptr_patch%edges%cell_blk(ile4,ibe4,1) == ibc2 ) THEN
              iie4(je) = 2
            ELSE IF (ptr_patch%edges%cell_idx(ile4,ibe4,2) == ilc2 .AND. &
              & ptr_patch%edges%cell_blk(ile4,ibe4,2) == ibc2 ) THEN
              iie4(je) = 4
            ENDIF

            IF (je == ptr_patch%cells%edge_idx(ilc1,ibc1,1) .AND. &
              & jb == ptr_patch%cells%edge_blk(ilc1,ibc1,1)) THEN

              inb1 = inv_neighbor_id(ilc1,ibc1,1)
              ie4  = MOD(inb1  ,3)+1
              ie5  = MOD(inb1+1,3)+1

              ptr_int%e_flx_avg(je,1,jb) = 0.5_wp *(                                        &
                & ( ptr_int%geofac_div(ilc1,1,ibc1)*ptr_int%c_bln_avg(ilc1,1,ibc1)            &
                & + ptr_int%geofac_div(ilc2,inb1,ibc2)*ptr_int%c_bln_avg(ilc1,2,ibc1)         &
                & - ptr_int%e_flx_avg(ile1,iie1(je),ibe1)*ptr_int%geofac_div(ilc1,2,ibc1)     &
                & - ptr_int%e_flx_avg(ile2,iie2(je),ibe2)*ptr_int%geofac_div(ilc1,3,ibc1) )   &
                & / ptr_int%geofac_div(ilc1,1,ibc1)   +                                       &
                & ( ptr_int%geofac_div(ilc2,inb1,ibc2)*ptr_int%c_bln_avg(ilc2,1,ibc2)         &
                & + ptr_int%geofac_div(ilc1,1,ibc1)*ptr_int%c_bln_avg(ilc2,inb1+1,ibc2)       &
                & - ptr_int%e_flx_avg(ile3,iie3(je),ibe3)*ptr_int%geofac_div(ilc2,ie4,ibc2)   &
                & - ptr_int%e_flx_avg(ile4,iie4(je),ibe4)*ptr_int%geofac_div(ilc2,ie5,ibc2) ) &
                & / ptr_int%geofac_div(ilc2,inb1,ibc2) )


            ELSE IF (je == ptr_patch%cells%edge_idx(ilc1,ibc1,2) .AND. &
              & jb == ptr_patch%cells%edge_blk(ilc1,ibc1,2)) THEN

              inb2 = inv_neighbor_id(ilc1,ibc1,2)
              ie4  = MOD(inb2  ,3)+1
              ie5  = MOD(inb2+1,3)+1

              ptr_int%e_flx_avg(je,1,jb) = 0.5_wp *(                                        &
                & ( ptr_int%geofac_div(ilc1,2,ibc1)*ptr_int%c_bln_avg(ilc1,1,ibc1)            &
                & + ptr_int%geofac_div(ilc2,inb2,ibc2)*ptr_int%c_bln_avg(ilc1,3,ibc1)         &
                & - ptr_int%e_flx_avg(ile1,iie1(je),ibe1)*ptr_int%geofac_div(ilc1,3,ibc1)     &
                & - ptr_int%e_flx_avg(ile2,iie2(je),ibe2)*ptr_int%geofac_div(ilc1,1,ibc1) )   &
                & / ptr_int%geofac_div(ilc1,2,ibc1)   +                                       &
                & ( ptr_int%geofac_div(ilc2,inb2,ibc2)*ptr_int%c_bln_avg(ilc2,1,ibc2)         &
                & + ptr_int%geofac_div(ilc1,2,ibc1)*ptr_int%c_bln_avg(ilc2,inb2+1,ibc2)       &
                & - ptr_int%e_flx_avg(ile3,iie3(je),ibe3)*ptr_int%geofac_div(ilc2,ie4,ibc2)   &
                & - ptr_int%e_flx_avg(ile4,iie4(je),ibe4)*ptr_int%geofac_div(ilc2,ie5,ibc2) ) &
                & / ptr_int%geofac_div(ilc2,inb2,ibc2) )


            ELSE IF (je == ptr_patch%cells%edge_idx(ilc1,ibc1,3) .AND. &
              & jb == ptr_patch%cells%edge_blk(ilc1,ibc1,3)) THEN

              inb3 = inv_neighbor_id(ilc1,ibc1,3)
              ie4  = MOD(inb3  ,3)+1
              ie5  = MOD(inb3+1,3)+1

              ptr_int%e_flx_avg(je,1,jb) = 0.5_wp *(                                        &
                & ( ptr_int%geofac_div(ilc1,3,ibc1)*ptr_int%c_bln_avg(ilc1,1,ibc1)            &
                & + ptr_int%geofac_div(ilc2,inb3,ibc2)*ptr_int%c_bln_avg(ilc1,4,ibc1)         &
                & - ptr_int%e_flx_avg(ile1,iie1(je),ibe1)*ptr_int%geofac_div(ilc1,1,ibc1)     &
                & - ptr_int%e_flx_avg(ile2,iie2(je),ibe2)*ptr_int%geofac_div(ilc1,2,ibc1) )   &
                & / ptr_int%geofac_div(ilc1,3,ibc1)   +                                       &
                & ( ptr_int%geofac_div(ilc2,inb3,ibc2)*ptr_int%c_bln_avg(ilc2,1,ibc2)         &
                & + ptr_int%geofac_div(ilc1,3,ibc1)*ptr_int%c_bln_avg(ilc2,inb3+1,ibc2)       &
                & - ptr_int%e_flx_avg(ile3,iie3(je),ibe3)*ptr_int%geofac_div(ilc2,ie4,ibc2)   &
                & - ptr_int%e_flx_avg(ile4,iie4(je),ibe4)*ptr_int%geofac_div(ilc2,ie5,ibc2) ) &
                & / ptr_int%geofac_div(ilc2,inb3,ibc2) )

            ENDIF

          ENDDO !edge loop

        END DO !block loop
!$OMP END DO

        ! Finally, the weighting coefficients are scaled in order to
        ! yield the right result for a constant wind field

!$OMP DO PRIVATE(jb,je,i_startidx,i_endidx,ile1,ibe1,ile2,ibe2,ile3,ibe3,ile4,ibe4, &
!$OMP            z_nx1,z_nx2,z_nx3,z_nx4,z_nx5) ICON_OMP_DEFAULT_SCHEDULE
        DO jb = i_startblk, i_endblk

          CALL get_indices_e(ptr_patch, jb, i_startblk, i_endblk, &
            & i_startidx, i_endidx, rl_start, rl_end)

          DO je = i_startidx, i_endidx

            IF(.NOT. ptr_patch%edges%decomp_info%owner_mask(je,jb)) CYCLE

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
              & + DOT_PRODUCT(z_nx1(je,1:3),z_nx2(je,1:3))*ptr_int%e_flx_avg(je,2,jb) &
              & + DOT_PRODUCT(z_nx1(je,1:3),z_nx3(je,1:3))*ptr_int%e_flx_avg(je,3,jb) &
              & + DOT_PRODUCT(z_nx1(je,1:3),z_nx4(je,1:3))*ptr_int%e_flx_avg(je,4,jb) &
              & + DOT_PRODUCT(z_nx1(je,1:3),z_nx5(je,1:3))*ptr_int%e_flx_avg(je,5,jb)

            ptr_int%e_flx_avg(je,1,jb) = ptr_int%e_flx_avg(je,1,jb)/checksum(je,jb)
            ptr_int%e_flx_avg(je,2,jb) = ptr_int%e_flx_avg(je,2,jb)/checksum(je,jb)
            ptr_int%e_flx_avg(je,3,jb) = ptr_int%e_flx_avg(je,3,jb)/checksum(je,jb)
            ptr_int%e_flx_avg(je,4,jb) = ptr_int%e_flx_avg(je,4,jb)/checksum(je,jb)
            ptr_int%e_flx_avg(je,5,jb) = ptr_int%e_flx_avg(je,5,jb)/checksum(je,jb)

          ENDDO !edge loop

        END DO !block loop
!$OMP END DO NOWAIT
!$OMP END PARALLEL

        CALL sync_patch_array(sync_c,ptr_patch,ptr_int%c_bln_avg)
        CALL sync_patch_array(sync_e,ptr_patch,ptr_int%e_flx_avg)

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
        & i_startidx, i_endidx, rl_start, rl_end)

      DO jc = i_startidx, i_endidx

        sum1 = sum1 + resid(jc,jb)**2
        wgt_loc_sum(jc,jb) = SUM(ptr_int%c_bln_avg(jc,1:4,jb))

        WRITE(710+ptr_patch%id,'(2i5,5f12.6,e13.5)') jb,jc,ptr_int%c_bln_avg(jc,1:4,jb),&
          & wgt_loc_sum(jc,jb),resid(jc,jb)

      END DO
    END DO
    WRITE(710+ptr_patch%id,'(4e13.5)') MAXVAL(resid),SQRT(sum1/ptr_patch%n_patch_cells),&
      & MAXVAL(wgt_loc_sum)-1._wp,MINVAL(wgt_loc_sum)-1._wp
    CLOSE (710+ptr_patch%id)

    ! Debug output for mass flux averaging weights

    rl_start = 5
    rl_end = min_rledge
    i_startblk = ptr_patch%edges%start_blk(rl_start,1)
    i_endblk   = ptr_patch%edges%end_blk(rl_end,i_nchdom)

    DO jb = i_startblk, i_endblk

      CALL get_indices_e(ptr_patch, jb, i_startblk, i_endblk, &
        & i_startidx, i_endidx, rl_start, rl_end)

      DO je = i_startidx, i_endidx

        WRITE(720+ptr_patch%id,'(2i5,6f12.6)') jb,je,ptr_int%e_flx_avg(je,1:5,jb),&
          & checksum(je,jb)

      END DO
    END DO
    CLOSE (720+ptr_patch%id)

#endif

  END SUBROUTINE force_mass_conservation_to_cellavg_wgt
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
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

    ! write(0,*) nudge_zone_width, grf_nudge_start_c, max_rlval

    IF (max_rlval < nudge_zone_width+grf_nudge_start_c-1) THEN
      CALL finish('init_nudgecoeffs',&
        & 'bdy_indexing_depth in prepare_gridref must be at least nudge_zone_width+4')
    ENDIF

    i_nchdom   = MAX(1,ptr_patch%n_childdom)

!$OMP PARALLEL PRIVATE(i_startblk,i_endblk)

    ! a) Nudging coefficients for cells
    i_startblk = ptr_patch%cells%start_blk(grf_nudge_start_c,1)
    i_endblk   = ptr_patch%cells%end_blk(0,i_nchdom)

!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
        & i_startidx, i_endidx, grf_nudge_start_c, 0)

      DO jc = i_startidx, i_endidx

        IF (ptr_patch%cells%refin_ctrl(jc,jb) > 0 .AND. &
          & ptr_patch%cells%refin_ctrl(jc,jb) <= nudge_zone_width+grf_nudge_start_c-1) THEN
          ptr_int%nudgecoeff_c(jc,jb) = &
            & nudge_max_coeff*EXP(-REAL(ptr_patch%cells%refin_ctrl(jc,jb)-grf_nudge_start_c,wp) / &
            & nudge_efold_width)
        ENDIF

      ENDDO !cell loop

    END DO !block loop
!$OMP END DO

    ! b) Nudging coefficients for edges
    i_startblk = ptr_patch%edges%start_blk(grf_nudge_start_e,1)
    i_endblk   = ptr_patch%edges%end_blk(0,i_nchdom)

    IF (ptr_patch%id > 1 .AND. lfeedback(ptr_patch%id)) THEN
      ! Use nudging coefficients optimized for velocity boundary diffusion
!$OMP DO PRIVATE(jb,je,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk, i_endblk

        CALL get_indices_e(ptr_patch, jb, i_startblk, i_endblk, &
          & i_startidx, i_endidx, grf_nudge_start_e, 0)

        DO je = i_startidx, i_endidx

          IF (ptr_patch%edges%refin_ctrl(je,jb) > 0 .AND. &
            & ptr_patch%edges%refin_ctrl(je,jb) <= grf_nudge_start_e+9) THEN
            ptr_int%nudgecoeff_e(je,jb) = nudge_max_coeff* &
              & EXP(-REAL(ptr_patch%edges%refin_ctrl(je,jb)-grf_nudge_start_e,wp) / 4._wp)
          ENDIF

        ENDDO !edge loop

      END DO !block loop
!$OMP END DO NOWAIT
    ELSE
      ! Use nudging coefficients from namelist
!$OMP DO PRIVATE(jb,je,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk, i_endblk

        CALL get_indices_e(ptr_patch, jb, i_startblk, i_endblk, &
          & i_startidx, i_endidx, grf_nudge_start_e, 0)

        DO je = i_startidx, i_endidx

          IF (ptr_patch%edges%refin_ctrl(je,jb) > 0 .AND. &
            & ptr_patch%edges%refin_ctrl(je,jb) <= 2*nudge_zone_width+grf_nudge_start_e-3) THEN
            ptr_int%nudgecoeff_e(je,jb) = &
              & nudge_max_coeff*EXP(-REAL(ptr_patch%edges%refin_ctrl(je,jb)-grf_nudge_start_e,wp) / &
              & (2._wp*nudge_efold_width))
          ENDIF

        ENDDO !edge loop

      END DO !block loop
!$OMP END DO NOWAIT
    ENDIF
!$OMP END PARALLEL

    CALL sync_patch_array(sync_c,ptr_patch,ptr_int%nudgecoeff_c)
    CALL sync_patch_array(sync_e,ptr_patch,ptr_int%nudgecoeff_e)

  END SUBROUTINE init_nudgecoeffs
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
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
      & ilv1, ilv2, ibv1, ibv2, jr1, ilr, ibr,                  &
      & ile1, ibe1, ile2, ibe2, ile3, ibe3
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
!$OMP DO PRIVATE(jb,je,jc,i_startidx,i_endidx,ile,ibe) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
        & i_startidx, i_endidx, rl_start, rl_end)

      DO je = 1, ptr_patch%geometry_info%cell_type
        DO jc = i_startidx, i_endidx

          IF (je > ptr_patch%cells%num_edges(jc,jb)) CYCLE ! relevant for hexagons
!             write(0,*) jc, jb, ":", ptr_patch%cells%num_edges(jc,jb)
!             STOP
!           ENDIF

          ile = ptr_patch%cells%edge_idx(jc,jb,je)
          ibe = ptr_patch%cells%edge_blk(jc,jb,je)

          ptr_int%geofac_div(jc,je,jb) = &
            & ptr_patch%edges%primal_edge_length(ile,ibe) * &
            & ptr_patch%cells%edge_orientation(jc,jb,je)  / &
            & ptr_patch%cells%area(jc,jb)

        ENDDO !cell loop
      ENDDO

    END DO !block loop
!$OMP END DO

    ! b) Geometrical factor for rotation
    rl_start = 2
    rl_end = min_rlvert

    ! Vorticity should have the right sign
    SELECT CASE (ptr_patch%geometry_info%cell_type)
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
!$OMP DO PRIVATE(jb,je,jv,i_startidx,i_endidx,ile,ibe) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_v(ptr_patch, jb, i_startblk, i_endblk, &
        & i_startidx, i_endidx, rl_start, rl_end)

      DO je = 1, 9-ptr_patch%geometry_info%cell_type
        DO jv = i_startidx, i_endidx

          IF(.NOT. ptr_patch%verts%decomp_info%owner_mask(jv,jb)) CYCLE

          IF (je > ptr_patch%verts%num_edges(jv,jb)) CYCLE

          ile = ptr_patch%verts%edge_idx(jv,jb,je)
          ibe = ptr_patch%verts%edge_blk(jv,jb,je)

          ptr_int%geofac_rot(jv,je,jb) =                &
            & ptr_patch%edges%dual_edge_length(ile,ibe) * &
            & ptr_patch%verts%edge_orientation(jv,jb,je)/ &
            & ptr_patch%verts%dual_area(jv,jb) * REAL(ifac,wp)

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
!$OMP    ilc2,ibc2,ilnc,ibnc) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
        & i_startidx, i_endidx, rl_start, rl_end)

      DO je = 1, ptr_patch%geometry_info%cell_type
        DO jc = i_startidx, i_endidx

          ile = ptr_patch%cells%edge_idx(jc,jb,je)
          ibe = ptr_patch%cells%edge_blk(jc,jb,je)

          ilc1 = ptr_patch%edges%cell_idx(ile,ibe,1)
          ibc1 = ptr_patch%edges%cell_blk(ile,ibe,1)
          ilc2 = ptr_patch%edges%cell_idx(ile,ibe,2)
          ibc2 = ptr_patch%edges%cell_blk(ile,ibe,2)

          IF (jc == ilc1 .AND. jb == ibc1) THEN
            IF (ptr_patch%geometry_info%cell_type == 3) THEN
              ptr_int%geofac_n2s(jc,1,jb) = ptr_int%geofac_n2s(jc,1,jb) - &
                & ptr_int%geofac_div(jc,je,jb) /                            &
                & ptr_patch%edges%dual_edge_length(ile,ibe)
            ELSE IF (ptr_patch%geometry_info%cell_type == 6) THEN
              ptr_int%geofac_n2s(jc,1,jb) = ptr_int%geofac_n2s(jc,1,jb) - &
                & ptr_int%geofac_div(jc,je,jb) /                            &
                & ptr_patch%edges%dual_edge_length(ile,ibe)*                &
                & ptr_patch%edges%tangent_orientation(ile,ibe)
            ENDIF
          ELSE IF (jc == ilc2 .AND. jb == ibc2) THEN
            IF (ptr_patch%geometry_info%cell_type == 3) THEN
              ptr_int%geofac_n2s(jc,1,jb) = ptr_int%geofac_n2s(jc,1,jb) + &
                & ptr_int%geofac_div(jc,je,jb) /                            &
                & ptr_patch%edges%dual_edge_length(ile,ibe)
            ELSE IF (ptr_patch%geometry_info%cell_type == 6) THEN
              ptr_int%geofac_n2s(jc,1,jb) = ptr_int%geofac_n2s(jc,1,jb) + &
                & ptr_int%geofac_div(jc,je,jb) /                            &
                & ptr_patch%edges%dual_edge_length(ile,ibe)*                &
                & ptr_patch%edges%tangent_orientation(ile,ibe)
            ENDIF
          ENDIF
          DO ic = 1, ptr_patch%geometry_info%cell_type
            ilnc = ptr_patch%cells%neighbor_idx(jc,jb,ic)
            ibnc = ptr_patch%cells%neighbor_blk(jc,jb,ic)
            IF (ilnc == ilc1 .AND. ibnc == ibc1) THEN
              IF (ptr_patch%geometry_info%cell_type == 3) THEN
                ptr_int%geofac_n2s(jc,ic+1,jb) = ptr_int%geofac_n2s(jc,ic+1,jb) - &
                  & ptr_int%geofac_div(jc,je,jb) /                                  &
                  & ptr_patch%edges%dual_edge_length(ile,ibe)
              ELSE IF (ptr_patch%geometry_info%cell_type == 6) THEN
                ptr_int%geofac_n2s(jc,ic+1,jb) = ptr_int%geofac_n2s(jc,ic+1,jb) - &
                  & ptr_int%geofac_div(jc,je,jb) /                                  &
                  & ptr_patch%edges%dual_edge_length(ile,ibe)*                      &
                  & ptr_patch%edges%tangent_orientation(ile,ibe)
              ENDIF
            ELSE IF (ilnc == ilc2 .AND. ibnc == ibc2) THEN
              IF (ptr_patch%geometry_info%cell_type == 3) THEN
                ptr_int%geofac_n2s(jc,ic+1,jb) = ptr_int%geofac_n2s(jc,ic+1,jb) + &
                  & ptr_int%geofac_div(jc,je,jb) /                                  &
                  & ptr_patch%edges%dual_edge_length(ile,ibe)
              ELSE IF (ptr_patch%geometry_info%cell_type == 6) THEN
                ptr_int%geofac_n2s(jc,ic+1,jb) = ptr_int%geofac_n2s(jc,ic+1,jb) + &
                  & ptr_int%geofac_div(jc,je,jb) /                                  &
                  & ptr_patch%edges%dual_edge_length(ile,ibe)*                      &
                  & ptr_patch%edges%tangent_orientation(ile,ibe)
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

    IF (ptr_patch%geometry_info%cell_type == 3) THEN

      rl_start = 2
      rl_end = min_rledge

      ! values for the blocking
      i_startblk = ptr_patch%edges%start_blk(rl_start,1)
      i_endblk   = ptr_patch%edges%end_blk(rl_end,i_nchdom)

!$OMP DO PRIVATE(jb,je,je1,i_startidx,i_endidx,ile,ibe) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk, i_endblk

        CALL get_indices_e(ptr_patch, jb, i_startblk, i_endblk, &
          & i_startidx, i_endidx, rl_start, rl_end)

        DO je1 = 1, 4
          DO je = i_startidx, i_endidx

            IF(.NOT. ptr_patch%edges%decomp_info%owner_mask(je,jb)) CYCLE

            ile = ptr_patch%edges%quad_idx(je,jb,je1)
            ibe = ptr_patch%edges%quad_blk(je,jb,je1)

            ptr_int%geofac_qdiv(je,je1,jb) = &
              & ptr_patch%edges%primal_edge_length(ile,ibe) * &
              & ptr_patch%edges%quad_orientation(je,jb,je1)  / &
              & ptr_patch%edges%quad_area(je,jb)

          ENDDO !edge loop
        ENDDO

      END DO !block loop
!$OMP END DO

    ENDIF

    ! e) Geometrical factor for gradient of divergence (triangles only)

    ! sync does not work on patch 0 for some unknown reason. But we don't need this
    ! field on the radiation grid anyway, so let's just skip it
    IF (ptr_patch%geometry_info%cell_type == 3 .AND. ptr_patch%id >= 1) THEN

      rl_start = 2
      rl_end = min_rledge

      ! values for the blocking
      i_startblk = ptr_patch%edges%start_blk(rl_start,1)
      i_endblk   = ptr_patch%edges%end_blk(rl_end,i_nchdom)

!$OMP DO PRIVATE(jb,je,je1,i_startidx,i_endidx,ile,ibe,ile1,ibe1,ile2,ibe2,ile3,ibe3,&
!$OMP            ilc1,ilc2,ibc1,ibc2) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk, i_endblk

        CALL get_indices_e(ptr_patch, jb, i_startblk, i_endblk, &
          & i_startidx, i_endidx, rl_start, rl_end)

        DO je = i_startidx, i_endidx

          IF(.NOT. ptr_patch%edges%decomp_info%owner_mask(je,jb)) CYCLE

          ilc1 = ptr_patch%edges%cell_idx(je,jb,1)
          ibc1 = ptr_patch%edges%cell_blk(je,jb,1)
          ilc2 = ptr_patch%edges%cell_idx(je,jb,2)
          ibc2 = ptr_patch%edges%cell_blk(je,jb,2)

          ! First consider edges of neighbor cell 2
          ile1 = ptr_patch%cells%edge_idx(ilc2,ibc2,1)
          ibe1 = ptr_patch%cells%edge_blk(ilc2,ibc2,1)
          ile2 = ptr_patch%cells%edge_idx(ilc2,ibc2,2)
          ibe2 = ptr_patch%cells%edge_blk(ilc2,ibc2,2)
          ile3 = ptr_patch%cells%edge_idx(ilc2,ibc2,3)
          ibe3 = ptr_patch%cells%edge_blk(ilc2,ibc2,3)

          IF (je == ile1 .AND. jb == ibe1) THEN
            ptr_int%geofac_grdiv(je,1,jb) = ptr_int%geofac_div(ilc2,1,ibc2)
          ELSE IF (je == ile2 .AND. jb == ibe2) THEN
            ptr_int%geofac_grdiv(je,1,jb) = ptr_int%geofac_div(ilc2,2,ibc2)
          ELSE IF (je == ile3 .AND. jb == ibe3) THEN
            ptr_int%geofac_grdiv(je,1,jb) = ptr_int%geofac_div(ilc2,3,ibc2)
          ENDIF

          ! Now consider edges of neighbor cell 1 and compute gradient
          ile1 = ptr_patch%cells%edge_idx(ilc1,ibc1,1)
          ibe1 = ptr_patch%cells%edge_blk(ilc1,ibc1,1)
          ile2 = ptr_patch%cells%edge_idx(ilc1,ibc1,2)
          ibe2 = ptr_patch%cells%edge_blk(ilc1,ibc1,2)
          ile3 = ptr_patch%cells%edge_idx(ilc1,ibc1,3)
          ibe3 = ptr_patch%cells%edge_blk(ilc1,ibc1,3)

          IF (je == ile1 .AND. jb == ibe1) THEN
            ptr_int%geofac_grdiv(je,1,jb) = (ptr_int%geofac_grdiv(je,1,jb) - &
              & ptr_int%geofac_div(ilc1,1,ibc1))*ptr_patch%edges%inv_dual_edge_length(je,jb)
          ELSE IF (je == ile2 .AND. jb == ibe2) THEN
            ptr_int%geofac_grdiv(je,1,jb) = (ptr_int%geofac_grdiv(je,1,jb) - &
              & ptr_int%geofac_div(ilc1,2,ibc1))*ptr_patch%edges%inv_dual_edge_length(je,jb)
          ELSE IF (je == ile3 .AND. jb == ibe3) THEN
            ptr_int%geofac_grdiv(je,1,jb) = (ptr_int%geofac_grdiv(je,1,jb) - &
              & ptr_int%geofac_div(ilc1,3,ibc1))*ptr_patch%edges%inv_dual_edge_length(je,jb)
          ENDIF
        ENDDO

        ! The quad edge indices are computed such that edges 1 and 2 border
        ! to neighbor cell 1 and edges 3 and 4 border to neighbor cell 2.
        ! Thus, splitting the following loop saves case discriminations
        DO je1 = 1, 2
          DO je = i_startidx, i_endidx

            IF(.NOT. ptr_patch%edges%decomp_info%owner_mask(je,jb)) CYCLE

            ile = ptr_patch%edges%quad_idx(je,jb,je1)
            ibe = ptr_patch%edges%quad_blk(je,jb,je1)

            ilc1 = ptr_patch%edges%cell_idx(je,jb,1)
            ibc1 = ptr_patch%edges%cell_blk(je,jb,1)

            ile1 = ptr_patch%cells%edge_idx(ilc1,ibc1,1)
            ibe1 = ptr_patch%cells%edge_blk(ilc1,ibc1,1)
            ile2 = ptr_patch%cells%edge_idx(ilc1,ibc1,2)
            ibe2 = ptr_patch%cells%edge_blk(ilc1,ibc1,2)
            ile3 = ptr_patch%cells%edge_idx(ilc1,ibc1,3)
            ibe3 = ptr_patch%cells%edge_blk(ilc1,ibc1,3)

            IF (ile == ile1 .AND. ibe == ibe1) THEN
              ptr_int%geofac_grdiv(je,je1+1,jb) = - &
                & ptr_int%geofac_div(ilc1,1,ibc1)*ptr_patch%edges%inv_dual_edge_length(je,jb)
            ELSE IF (ile == ile2 .AND. ibe == ibe2) THEN
              ptr_int%geofac_grdiv(je,je1+1,jb) = - &
                & ptr_int%geofac_div(ilc1,2,ibc1)*ptr_patch%edges%inv_dual_edge_length(je,jb)
            ELSE IF (ile == ile3 .AND. ibe == ibe3) THEN
              ptr_int%geofac_grdiv(je,je1+1,jb) = - &
                & ptr_int%geofac_div(ilc1,3,ibc1)*ptr_patch%edges%inv_dual_edge_length(je,jb)
            ENDIF

          ENDDO !edge loop
        ENDDO

        DO je1 = 3, 4
          DO je = i_startidx, i_endidx

            IF(.NOT. ptr_patch%edges%decomp_info%owner_mask(je,jb)) CYCLE

            ile = ptr_patch%edges%quad_idx(je,jb,je1)
            ibe = ptr_patch%edges%quad_blk(je,jb,je1)

            ilc2 = ptr_patch%edges%cell_idx(je,jb,2)
            ibc2 = ptr_patch%edges%cell_blk(je,jb,2)

            ile1 = ptr_patch%cells%edge_idx(ilc2,ibc2,1)
            ibe1 = ptr_patch%cells%edge_blk(ilc2,ibc2,1)
            ile2 = ptr_patch%cells%edge_idx(ilc2,ibc2,2)
            ibe2 = ptr_patch%cells%edge_blk(ilc2,ibc2,2)
            ile3 = ptr_patch%cells%edge_idx(ilc2,ibc2,3)
            ibe3 = ptr_patch%cells%edge_blk(ilc2,ibc2,3)

            IF (ile == ile1 .AND. ibe == ibe1) THEN
              ptr_int%geofac_grdiv(je,je1+1,jb) = &
                & ptr_int%geofac_div(ilc2,1,ibc2)*ptr_patch%edges%inv_dual_edge_length(je,jb)
            ELSE IF (ile == ile2 .AND. ibe == ibe2) THEN
              ptr_int%geofac_grdiv(je,je1+1,jb) = &
                & ptr_int%geofac_div(ilc2,2,ibc2)*ptr_patch%edges%inv_dual_edge_length(je,jb)
            ELSE IF (ile == ile3 .AND. ibe == ibe3) THEN
              ptr_int%geofac_grdiv(je,je1+1,jb) = &
                & ptr_int%geofac_div(ilc2,3,ibc2)*ptr_patch%edges%inv_dual_edge_length(je,jb)
            ENDIF

          ENDDO !edge loop
        ENDDO

      END DO !block loop
!$OMP END DO

    ENDIF

    ! f) coefficients for directional gradient of a normal vector quantity
    ! at the same edge (gives directional laplacian if gradient psi is assumed as input)

    IF (ptr_patch%geometry_info%cell_type == 6) THEN

      ! Now compute coefficients
      !-------------------------
      rl_start = 2
      rl_end = min_rledge
      ! values for the blocking
      i_startblk = ptr_patch%edges%start_blk(rl_start,1)
      i_endblk   = ptr_patch%edges%end_blk(rl_end,i_nchdom)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,je,z_pn_k,z_pt_k,&
!$OMP ilc1,ibc1,ilc2,ibc2,ilv1,ibv1,ilv2,ibv2,je1,ile,ibe,z_proj,z_pn_j,&
!$OMP jm,jn,nincr,ilr,ibr,jr1,z_lon,z_lat,z_norm,z_cart_no,z_cart_ea) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk, i_endblk
        CALL get_indices_e(ptr_patch, jb, i_startblk, i_endblk, &
          & i_startidx, i_endidx, rl_start, rl_end)

        DO je = i_startidx, i_endidx

          IF(.NOT. ptr_patch%edges%decomp_info%owner_mask(je,jb)) CYCLE

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
            & DOT_PRODUCT(z_cart_no%x(:),ptr_patch%edges%primal_cart_normal(je,jb)%x(:))
          ptr_int%cea_en(je,1,jb) = &
            & DOT_PRODUCT(z_cart_ea%x(:),ptr_patch%edges%primal_cart_normal(je,jb)%x(:))

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
            & DOT_PRODUCT(z_cart_no%x(:),ptr_patch%edges%primal_cart_normal(je,jb)%x(:))
          ptr_int%cea_en(je,2,jb) = &
            & DOT_PRODUCT(z_cart_ea%x(:),ptr_patch%edges%primal_cart_normal(je,jb)%x(:))

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
                  & ibe == ptr_int%dir_gradt_b1(jn,je,jb)) THEN
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
                  & ibe == ptr_int%dir_gradt_b2(jn,je,jb)) THEN
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

    ! g) Geometrical factor for Green-Gauss gradient
    rl_start = 2
    rl_end = min_rlcell

    ! values for the blocking
    i_startblk = ptr_patch%cells%start_blk(rl_start,1)
    i_endblk   = ptr_patch%cells%end_blk(rl_end,i_nchdom)
    !
    ! loop through all patch cells (and blocks)
    !
!$OMP DO PRIVATE(jb,je,jc,ic,i_startidx,i_endidx,ile,ibe,ilc1,ibc1,&
!$OMP    ilc2,ibc2,ilnc,ibnc) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
        & i_startidx, i_endidx, rl_start, rl_end)

      DO je = 1, ptr_patch%geometry_info%cell_type
        DO jc = i_startidx, i_endidx

          IF(.NOT. ptr_patch%cells%decomp_info%owner_mask(jc,jb)) CYCLE

          ile = ptr_patch%cells%edge_idx(jc,jb,je)
          ibe = ptr_patch%cells%edge_blk(jc,jb,je)

          ilc1 = ptr_patch%edges%cell_idx(ile,ibe,1)
          ibc1 = ptr_patch%edges%cell_blk(ile,ibe,1)
          ilc2 = ptr_patch%edges%cell_idx(ile,ibe,2)
          ibc2 = ptr_patch%edges%cell_blk(ile,ibe,2)

          IF (jc == ilc1 .AND. jb == ibc1) THEN
            ptr_int%geofac_grg(jc,1,jb,1) = ptr_int%geofac_grg(jc,1,jb,1) +      &
              & ptr_int%primal_normal_ec(jc,jb,je,1)*ptr_int%geofac_div(jc,je,jb)* &
              & ptr_int%c_lin_e(ile,1,ibe)
            ptr_int%geofac_grg(jc,1,jb,2) = ptr_int%geofac_grg(jc,1,jb,2) +      &
              & ptr_int%primal_normal_ec(jc,jb,je,2)*ptr_int%geofac_div(jc,je,jb)* &
              & ptr_int%c_lin_e(ile,1,ibe)
          ELSE IF (jc == ilc2 .AND. jb == ibc2) THEN
            ptr_int%geofac_grg(jc,1,jb,1) = ptr_int%geofac_grg(jc,1,jb,1) +      &
              & ptr_int%primal_normal_ec(jc,jb,je,1)*ptr_int%geofac_div(jc,je,jb)* &
              & ptr_int%c_lin_e(ile,2,ibe)
            ptr_int%geofac_grg(jc,1,jb,2) = ptr_int%geofac_grg(jc,1,jb,2) +      &
              & ptr_int%primal_normal_ec(jc,jb,je,2)*ptr_int%geofac_div(jc,je,jb)* &
              & ptr_int%c_lin_e(ile,2,ibe)
          ENDIF
          DO ic = 1, ptr_patch%geometry_info%cell_type
            ilnc = ptr_patch%cells%neighbor_idx(jc,jb,ic)
            ibnc = ptr_patch%cells%neighbor_blk(jc,jb,ic)
            IF (ilnc == ilc1 .AND. ibnc == ibc1) THEN
              ptr_int%geofac_grg(jc,ic+1,jb,1) = ptr_int%geofac_grg(jc,ic+1,jb,1)+ &
                & ptr_int%primal_normal_ec(jc,jb,je,1)*ptr_int%geofac_div(jc,je,jb)* &
                & ptr_int%c_lin_e(ile,1,ibe)
              ptr_int%geofac_grg(jc,ic+1,jb,2) = ptr_int%geofac_grg(jc,ic+1,jb,2)+ &
                & ptr_int%primal_normal_ec(jc,jb,je,2)*ptr_int%geofac_div(jc,je,jb)* &
                & ptr_int%c_lin_e(ile,1,ibe)
            ELSE IF (ilnc == ilc2 .AND. ibnc == ibc2) THEN
              ptr_int%geofac_grg(jc,ic+1,jb,1) = ptr_int%geofac_grg(jc,ic+1,jb,1)+ &
                & ptr_int%primal_normal_ec(jc,jb,je,1)*ptr_int%geofac_div(jc,je,jb)* &
                & ptr_int%c_lin_e(ile,2,ibe)
              ptr_int%geofac_grg(jc,ic+1,jb,2) = ptr_int%geofac_grg(jc,ic+1,jb,2)+ &
                & ptr_int%primal_normal_ec(jc,jb,je,2)*ptr_int%geofac_div(jc,je,jb)* &
                & ptr_int%c_lin_e(ile,2,ibe)
            ENDIF
          ENDDO

          ! To ensure that dummy edges have a factor of 0:
          IF (je > ptr_patch%cells%num_edges(jc,jb)) THEN
            ptr_int%geofac_grg(jc,je+1,jb,1:2) = 0._wp
          ENDIF

        ENDDO !cell loop
      ENDDO

    END DO !block loop
!$OMP END DO NOWAIT

!$OMP END PARALLEL

    CALL sync_patch_array(sync_c,ptr_patch,ptr_int%geofac_div)
    CALL sync_patch_array(sync_v,ptr_patch,ptr_int%geofac_rot)
    CALL sync_patch_array(sync_c,ptr_patch,ptr_int%geofac_n2s)

    IF (ptr_patch%geometry_info%cell_type == 3) THEN
      CALL sync_patch_array(sync_e,ptr_patch,ptr_int%geofac_qdiv)
      IF (ptr_patch%id >= 1) CALL sync_patch_array(sync_e,ptr_patch,ptr_int%geofac_grdiv)
    ENDIF

    IF (ptr_patch%geometry_info%cell_type == 6) THEN
      CALL sync_patch_array(sync_e,ptr_patch,ptr_int%cno_en)
      CALL sync_patch_array(sync_e,ptr_patch,ptr_int%cea_en)

      DO jm = 1, 6
        CALL sync_idx(sync_e,sync_e,ptr_patch,ptr_int%dir_gradh_i1(jm,:,:), &
          & ptr_int%dir_gradh_b1(jm,:,:))
        CALL sync_idx(sync_e,sync_e,ptr_patch,ptr_int%dir_gradh_i2(jm,:,:), &
          & ptr_int%dir_gradh_b2(jm,:,:))
        CALL sync_patch_array(sync_e,ptr_patch,ptr_int%dir_gradhux_c1(jm, :, :))
        CALL sync_patch_array(sync_e,ptr_patch,ptr_int%dir_gradhux_c2(jm, :, :))
        CALL sync_patch_array(sync_e,ptr_patch,ptr_int%strain_def_c1(jm, :, :))
        CALL sync_patch_array(sync_e,ptr_patch,ptr_int%strain_def_c2(jm, :, :))
      ENDDO

      DO jm = 1, 9
        CALL sync_idx(sync_e,sync_e,ptr_patch,ptr_int%dir_gradt_i1(jm,:,:), &
          & ptr_int%dir_gradt_b1(jm,:,:))
        CALL sync_idx(sync_e,sync_e,ptr_patch,ptr_int%dir_gradt_i2(jm,:,:), &
          & ptr_int%dir_gradt_b2(jm,:,:))
        CALL sync_patch_array(sync_e,ptr_patch,ptr_int%dir_gradtxy_v1(jm, :, :))
        CALL sync_patch_array(sync_e,ptr_patch,ptr_int%dir_gradtxy_v2(jm, :, :))
        CALL sync_patch_array(sync_e,ptr_patch,ptr_int%dir_gradtyx_v1(jm, :, :))
        CALL sync_patch_array(sync_e,ptr_patch,ptr_int%dir_gradtyx_v2(jm, :, :))
        CALL sync_patch_array(sync_e,ptr_patch,ptr_int%shear_def_v1(jm, :, :))
        CALL sync_patch_array(sync_e,ptr_patch,ptr_int%shear_def_v2(jm, :, :))
      ENDDO
    ENDIF

    CALL sync_patch_array(sync_c,ptr_patch,ptr_int%geofac_grg(:,:,:,1))
    CALL sync_patch_array(sync_c,ptr_patch,ptr_int%geofac_grg(:,:,:,2))



  END SUBROUTINE init_geo_factors
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !>
  !! Computes the local orientation of the edge primal normal and dual normal
  !! at the location of the cell centers and vertices.
  !! Moreover, the Cartesian orientation vectors of the edge primal normals
  !! are stored for use in the RBF initialization routines, and inverse
  !! primal and dual edge lengths are computed
  !!
  !! @par Revision History
  !!  developed by Guenther Zaengl, 2009-03-31
  !! Modification by Anurag Dipankar, MPIM, 2012-12-28
  !! -Removed the calculation of ptr_int%cart_edge/cell_coord as this
  !!  information is now stored in patch. Also adapted some calculations
  !!  for geometry based info
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
      & ilv3, ibv3, ilv4, ibv4, ile1, ibe1

    REAL(wp) :: z_nu, z_nv, z_lon, z_lat, z_nx1(3), z_nx2(3), z_norm

    TYPE(t_cartesian_coordinates) :: cc_cell, cc_ev3, cc_ev4, cc_v1, cc_v2, cc_v3, &
      & cc_dis1, cc_dis2, cc_dis3

    CHARACTER(LEN=*), PARAMETER :: method_name = 'mo_intp_coeffs:complete_patchinfo'

    !-----------------------------------------------------------------------

    i_nchdom   = MAX(1,ptr_patch%n_childdom)

!$OMP PARALLEL  PRIVATE(rl_start,rl_end,i_startblk,i_endblk)
    rl_start = 1
    rl_end = min_rledge

    ! values for the blocking
    i_startblk = ptr_patch%edges%start_blk(rl_start,1)
    i_endblk   = ptr_patch%edges%end_blk(rl_end,i_nchdom)
    !
    ! The fields for the inverse primal and dual edge lengths are
    ! initialized here.
    !
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,je) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_e(ptr_patch, jb, i_startblk, i_endblk, &
        & i_startidx, i_endidx, rl_start, rl_end)

      DO je =  i_startidx, i_endidx

        IF(.NOT.ptr_patch%edges%decomp_info%owner_mask(je,jb)) CYCLE

        ! compute inverse primal edge length
        ! (dual follows below in the rl_start=2 section)

        ptr_patch%edges%inv_primal_edge_length(je,jb) = &
          & 1._wp/ptr_patch%edges%primal_edge_length(je,jb)

      ENDDO

    END DO !block loop
!$OMP END DO

    rl_start = 2
    rl_end = min_rledge

    ! Second step: computed projected orientation vectors and related information
    i_startblk = ptr_patch%edges%start_blk(rl_start,1)
    ! i_endblk   = ptr_patch%edges%end_blk(rl_end,i_nchdom)
    i_endblk   = ptr_patch%nblks_e

    ! Initialization of lateral boundary points
    IF (ptr_patch%id > 1) THEN
      CALL init(ptr_patch%edges%inv_dual_edge_length(:,1:i_startblk))
      CALL init(ptr_patch%edges%vertex_idx(:,1:i_startblk,3))
      CALL init(ptr_patch%edges%vertex_idx(:,1:i_startblk,4))
      CALL init(ptr_patch%edges%vertex_blk(:,1:i_startblk,3))
      CALL init(ptr_patch%edges%vertex_blk(:,1:i_startblk,4))
      CALL init(ptr_patch%edges%inv_vert_vert_length(:,1:i_startblk))
      CALL init(ptr_patch%edges%primal_normal_cell(:,1:i_startblk,:)%v1)
      CALL init(ptr_patch%edges%dual_normal_cell  (:,1:i_startblk,:)%v1)
      CALL init(ptr_patch%edges%primal_normal_vert(:,1:i_startblk,:)%v1)
      CALL init(ptr_patch%edges%dual_normal_vert  (:,1:i_startblk,:)%v1)
      CALL init(ptr_patch%edges%primal_normal_cell(:,1:i_startblk,:)%v2)
      CALL init(ptr_patch%edges%dual_normal_cell  (:,1:i_startblk,:)%v2)
      CALL init(ptr_patch%edges%primal_normal_vert(:,1:i_startblk,:)%v2)
      CALL init(ptr_patch%edges%dual_normal_vert  (:,1:i_startblk,:)%v2)
!$OMP BARRIER
    ENDIF
    !
    ! loop through all patch edges
    !
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,je,ilc1,ibc1,ilv1,ibv1,ilc2,ibc2,ilv2, &
!$OMP            ibv2,ilv3,ibv3,ilv4,ibv4,z_nu,z_nv,z_lon,z_lat,z_nx1,z_nx2,   &
!$OMP            cc_ev3,cc_ev4,z_norm) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_e(ptr_patch, jb, i_startblk, i_endblk, &
        & i_startidx, i_endidx, rl_start, rl_end)

      DO je =  i_startidx, i_endidx

        IF(.NOT.ptr_patch%edges%decomp_info%owner_mask(je,jb)) CYCLE

        ! compute inverse dual edge length (undefined for refin_ctrl=1)

        ptr_patch%edges%inv_dual_edge_length(je,jb) = &
          & 1._wp/ptr_patch%edges%dual_edge_length(je,jb)

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
          & ptr_patch%edges%vertex_idx(je,jb,1) .OR.  &
          & ptr_patch%cells%vertex_blk(ilc1,ibc1,1) /= &
          & ptr_patch%edges%vertex_blk(je,jb,1)) .AND.  &
          & (ptr_patch%cells%vertex_idx(ilc1,ibc1,1) /= &
          & ptr_patch%edges%vertex_idx(je,jb,2) .OR.  &
          & ptr_patch%cells%vertex_blk(ilc1,ibc1,1) /= &
          & ptr_patch%edges%vertex_blk(je,jb,2)) )        THEN

          ptr_patch%edges%vertex_idx(je,jb,3) = ptr_patch%cells%vertex_idx(ilc1,ibc1,1)
          ptr_patch%edges%vertex_blk(je,jb,3) = ptr_patch%cells%vertex_blk(ilc1,ibc1,1)

        ELSE IF ((ptr_patch%cells%vertex_idx(ilc1,ibc1,2) /= &
          & ptr_patch%edges%vertex_idx(je,jb,1) .OR.  &
          & ptr_patch%cells%vertex_blk(ilc1,ibc1,2) /= &
          & ptr_patch%edges%vertex_blk(je,jb,1)) .AND.  &
          & (ptr_patch%cells%vertex_idx(ilc1,ibc1,2) /= &
          & ptr_patch%edges%vertex_idx(je,jb,2) .OR.  &
          & ptr_patch%cells%vertex_blk(ilc1,ibc1,2) /= &
          & ptr_patch%edges%vertex_blk(je,jb,2)) )        THEN

          ptr_patch%edges%vertex_idx(je,jb,3) = ptr_patch%cells%vertex_idx(ilc1,ibc1,2)
          ptr_patch%edges%vertex_blk(je,jb,3) = ptr_patch%cells%vertex_blk(ilc1,ibc1,2)

        ELSE IF ((ptr_patch%cells%vertex_idx(ilc1,ibc1,3) /= &
          & ptr_patch%edges%vertex_idx(je,jb,1) .OR.  &
          & ptr_patch%cells%vertex_blk(ilc1,ibc1,3) /= &
          & ptr_patch%edges%vertex_blk(je,jb,1)) .AND.  &
          & (ptr_patch%cells%vertex_idx(ilc1,ibc1,3) /= &
          & ptr_patch%edges%vertex_idx(je,jb,2) .OR.  &
          & ptr_patch%cells%vertex_blk(ilc1,ibc1,3) /= &
          & ptr_patch%edges%vertex_blk(je,jb,2)) )        THEN

          ptr_patch%edges%vertex_idx(je,jb,3) = ptr_patch%cells%vertex_idx(ilc1,ibc1,3)
          ptr_patch%edges%vertex_blk(je,jb,3) = ptr_patch%cells%vertex_blk(ilc1,ibc1,3)

        ELSE
          CALL finish(method_name, "Unresolved edges%vertex(3)")
        ENDIF

        IF ((ptr_patch%cells%vertex_idx(ilc2,ibc2,1) /= &
          & ptr_patch%edges%vertex_idx(je,jb,1) .OR.  &
          & ptr_patch%cells%vertex_blk(ilc2,ibc2,1) /= &
          & ptr_patch%edges%vertex_blk(je,jb,1)) .AND.  &
          & (ptr_patch%cells%vertex_idx(ilc2,ibc2,1) /= &
          & ptr_patch%edges%vertex_idx(je,jb,2) .OR.  &
          & ptr_patch%cells%vertex_blk(ilc2,ibc2,1) /= &
          & ptr_patch%edges%vertex_blk(je,jb,2)) )        THEN

          ptr_patch%edges%vertex_idx(je,jb,4) = ptr_patch%cells%vertex_idx(ilc2,ibc2,1)
          ptr_patch%edges%vertex_blk(je,jb,4) = ptr_patch%cells%vertex_blk(ilc2,ibc2,1)

        ELSE IF ((ptr_patch%cells%vertex_idx(ilc2,ibc2,2) /= &
          & ptr_patch%edges%vertex_idx(je,jb,1) .OR.  &
          & ptr_patch%cells%vertex_blk(ilc2,ibc2,2) /= &
          & ptr_patch%edges%vertex_blk(je,jb,1)) .AND.  &
          & (ptr_patch%cells%vertex_idx(ilc2,ibc2,2) /= &
          & ptr_patch%edges%vertex_idx(je,jb,2) .OR.  &
          & ptr_patch%cells%vertex_blk(ilc2,ibc2,2) /= &
          & ptr_patch%edges%vertex_blk(je,jb,2)) )        THEN

          ptr_patch%edges%vertex_idx(je,jb,4) = ptr_patch%cells%vertex_idx(ilc2,ibc2,2)
          ptr_patch%edges%vertex_blk(je,jb,4) = ptr_patch%cells%vertex_blk(ilc2,ibc2,2)

        ELSE IF ((ptr_patch%cells%vertex_idx(ilc2,ibc2,3) /= &
          & ptr_patch%edges%vertex_idx(je,jb,1) .OR.  &
          & ptr_patch%cells%vertex_blk(ilc2,ibc2,3) /= &
          & ptr_patch%edges%vertex_blk(je,jb,1)) .AND.  &
          & (ptr_patch%cells%vertex_idx(ilc2,ibc2,3) /= &
          & ptr_patch%edges%vertex_idx(je,jb,2) .OR.  &
          & ptr_patch%cells%vertex_blk(ilc2,ibc2,3) /= &
          & ptr_patch%edges%vertex_blk(je,jb,2)) )        THEN

          ptr_patch%edges%vertex_idx(je,jb,4) = ptr_patch%cells%vertex_idx(ilc2,ibc2,3)
          ptr_patch%edges%vertex_blk(je,jb,4) = ptr_patch%cells%vertex_blk(ilc2,ibc2,3)

        ELSE
          CALL finish(method_name, "Unresolved edges%vertex(4)")
        ENDIF

        ilv3 = ptr_patch%edges%vertex_idx(je,jb,3)
        ibv3 = ptr_patch%edges%vertex_blk(je,jb,3)
        ilv4 = ptr_patch%edges%vertex_idx(je,jb,4)
        ibv4 = ptr_patch%edges%vertex_blk(je,jb,4)

        cc_ev3 = ptr_patch%verts%cartesian(ilv3,ibv3)
        cc_ev4 = ptr_patch%verts%cartesian(ilv4,ibv4)

        ! inverse length bewtween vertices 3 and 4
        IF (ptr_patch%geometry_info%cell_type == 3 ) THEN
          ptr_patch%edges%inv_vert_vert_length(je,jb) = 1._wp/&
            & (grid_sphere_radius*arc_length(cc_ev3,cc_ev4,ptr_patch%geometry_info))
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

        CALL cvec2gvec(z_nx1(1),z_nx1(2),z_nx1(3),z_lon,z_lat,z_nu,z_nv,ptr_patch%geometry_info)
        z_norm = SQRT(z_nu*z_nu+z_nv*z_nv)

        ptr_patch%edges%primal_normal_cell(je,jb,1)%v1 = z_nu/z_norm
        ptr_patch%edges%primal_normal_cell(je,jb,1)%v2 = z_nv/z_norm

        CALL cvec2gvec(z_nx2(1),z_nx2(2),z_nx2(3),z_lon,z_lat,z_nu,z_nv,ptr_patch%geometry_info)
        z_norm = SQRT(z_nu*z_nu+z_nv*z_nv)

        ptr_patch%edges%dual_normal_cell(je,jb,1)%v1 = z_nu/z_norm
        ptr_patch%edges%dual_normal_cell(je,jb,1)%v2 = z_nv/z_norm

        ! get location of cell 2

        z_lon = ptr_patch%cells%center(ilc2,ibc2)%lon
        z_lat = ptr_patch%cells%center(ilc2,ibc2)%lat

        ! compute local primal and dual normals at cell 2

        CALL cvec2gvec(z_nx1(1),z_nx1(2),z_nx1(3),z_lon,z_lat,z_nu,z_nv,ptr_patch%geometry_info)
        z_norm = SQRT(z_nu*z_nu+z_nv*z_nv)

        ptr_patch%edges%primal_normal_cell(je,jb,2)%v1 = z_nu/z_norm
        ptr_patch%edges%primal_normal_cell(je,jb,2)%v2 = z_nv/z_norm

        CALL cvec2gvec(z_nx2(1),z_nx2(2),z_nx2(3),z_lon,z_lat,z_nu,z_nv,ptr_patch%geometry_info)
        z_norm = SQRT(z_nu*z_nu+z_nv*z_nv)

        ptr_patch%edges%dual_normal_cell(je,jb,2)%v1 = z_nu/z_norm
        ptr_patch%edges%dual_normal_cell(je,jb,2)%v2 = z_nv/z_norm

        ! get location of vertex 1

        z_lon = ptr_patch%verts%vertex(ilv1,ibv1)%lon
        z_lat = ptr_patch%verts%vertex(ilv1,ibv1)%lat

        ! compute local primal and dual normals at vertex 1

        CALL cvec2gvec(z_nx1(1),z_nx1(2),z_nx1(3),z_lon,z_lat,z_nu,z_nv,ptr_patch%geometry_info)
        z_norm = SQRT(z_nu*z_nu+z_nv*z_nv)

        ptr_patch%edges%primal_normal_vert(je,jb,1)%v1 = z_nu/z_norm
        ptr_patch%edges%primal_normal_vert(je,jb,1)%v2 = z_nv/z_norm

        CALL cvec2gvec(z_nx2(1),z_nx2(2),z_nx2(3),z_lon,z_lat,z_nu,z_nv,ptr_patch%geometry_info)
        z_norm = SQRT(z_nu*z_nu+z_nv*z_nv)

        ptr_patch%edges%dual_normal_vert(je,jb,1)%v1 = z_nu/z_norm
        ptr_patch%edges%dual_normal_vert(je,jb,1)%v2 = z_nv/z_norm

        ! get location of vertex 2

        z_lon = ptr_patch%verts%vertex(ilv2,ibv2)%lon
        z_lat = ptr_patch%verts%vertex(ilv2,ibv2)%lat

        ! compute local primal and dual normals at vertex 2

        CALL cvec2gvec(z_nx1(1),z_nx1(2),z_nx1(3),z_lon,z_lat,z_nu,z_nv,ptr_patch%geometry_info)
        z_norm = SQRT(z_nu*z_nu+z_nv*z_nv)

        ptr_patch%edges%primal_normal_vert(je,jb,2)%v1 = z_nu/z_norm
        ptr_patch%edges%primal_normal_vert(je,jb,2)%v2 = z_nv/z_norm

        CALL cvec2gvec(z_nx2(1),z_nx2(2),z_nx2(3),z_lon,z_lat,z_nu,z_nv,ptr_patch%geometry_info)
        z_norm = SQRT(z_nu*z_nu+z_nv*z_nv)

        ptr_patch%edges%dual_normal_vert(je,jb,2)%v1 = z_nu/z_norm
        ptr_patch%edges%dual_normal_vert(je,jb,2)%v2 = z_nv/z_norm

        ! get location of vertex 3

        z_lon = ptr_patch%verts%vertex(ilv3,ibv3)%lon
        z_lat = ptr_patch%verts%vertex(ilv3,ibv3)%lat

        ! compute local primal and dual normals at vertex 3

        CALL cvec2gvec(z_nx1(1),z_nx1(2),z_nx1(3),z_lon,z_lat,z_nu,z_nv,ptr_patch%geometry_info)
        z_norm = SQRT(z_nu*z_nu+z_nv*z_nv)

        ptr_patch%edges%primal_normal_vert(je,jb,3)%v1 = z_nu/z_norm
        ptr_patch%edges%primal_normal_vert(je,jb,3)%v2 = z_nv/z_norm

        CALL cvec2gvec(z_nx2(1),z_nx2(2),z_nx2(3),z_lon,z_lat,z_nu,z_nv,ptr_patch%geometry_info)
        z_norm = SQRT(z_nu*z_nu+z_nv*z_nv)

        ptr_patch%edges%dual_normal_vert(je,jb,3)%v1 = z_nu/z_norm
        ptr_patch%edges%dual_normal_vert(je,jb,3)%v2 = z_nv/z_norm

        ! get location of vertex 4

        z_lon = ptr_patch%verts%vertex(ilv4,ibv4)%lon
        z_lat = ptr_patch%verts%vertex(ilv4,ibv4)%lat

        ! compute local primal and dual normals at vertex 2

        CALL cvec2gvec(z_nx1(1),z_nx1(2),z_nx1(3),z_lon,z_lat,z_nu,z_nv,ptr_patch%geometry_info)
        z_norm = SQRT(z_nu*z_nu+z_nv*z_nv)

        ptr_patch%edges%primal_normal_vert(je,jb,4)%v1 = z_nu/z_norm
        ptr_patch%edges%primal_normal_vert(je,jb,4)%v2 = z_nv/z_norm

        CALL cvec2gvec(z_nx2(1),z_nx2(2),z_nx2(3),z_lon,z_lat,z_nu,z_nv,ptr_patch%geometry_info)
        z_norm = SQRT(z_nu*z_nu+z_nv*z_nv)

        ptr_patch%edges%dual_normal_vert(je,jb,4)%v1 = z_nu/z_norm
        ptr_patch%edges%dual_normal_vert(je,jb,4)%v2 = z_nv/z_norm

      ENDDO

    END DO !block loop
!$OMP END DO NOWAIT

!$OMP END PARALLEL

    ! primal_normal_cell must be sync'd before next loop,
    ! so do a sync for all above calculated quantities

    CALL sync_patch_array(sync_e,ptr_patch,ptr_patch%edges%inv_primal_edge_length)

    CALL sync_idx(sync_e,sync_v,ptr_patch,ptr_patch%edges%vertex_idx(:,:,3), &
      & ptr_patch%edges%vertex_blk(:,:,3))
    CALL sync_idx(sync_e,sync_v,ptr_patch,ptr_patch%edges%vertex_idx(:,:,4), &
      & ptr_patch%edges%vertex_blk(:,:,4))

    CALL sync_patch_array(sync_e,ptr_patch,ptr_patch%edges%inv_dual_edge_length)
    CALL sync_patch_array(sync_e,ptr_patch,ptr_patch%edges%inv_vert_vert_length)

    CALL sync_patch_array(sync_e,ptr_patch,ptr_patch%edges%primal_normal_cell(:,:,1)%v1)
    CALL sync_patch_array(sync_e,ptr_patch,ptr_patch%edges%primal_normal_cell(:,:,2)%v1)
    CALL sync_patch_array(sync_e,ptr_patch,ptr_patch%edges%primal_normal_vert(:,:,1)%v1)
    CALL sync_patch_array(sync_e,ptr_patch,ptr_patch%edges%primal_normal_vert(:,:,2)%v1)
    CALL sync_patch_array(sync_e,ptr_patch,ptr_patch%edges%primal_normal_vert(:,:,3)%v1)
    CALL sync_patch_array(sync_e,ptr_patch,ptr_patch%edges%primal_normal_vert(:,:,4)%v1)

    CALL sync_patch_array(sync_e,ptr_patch,ptr_patch%edges%dual_normal_cell(:,:,1)%v1)
    CALL sync_patch_array(sync_e,ptr_patch,ptr_patch%edges%dual_normal_cell(:,:,2)%v1)
    CALL sync_patch_array(sync_e,ptr_patch,ptr_patch%edges%dual_normal_vert(:,:,1)%v1)
    CALL sync_patch_array(sync_e,ptr_patch,ptr_patch%edges%dual_normal_vert(:,:,2)%v1)
    CALL sync_patch_array(sync_e,ptr_patch,ptr_patch%edges%dual_normal_vert(:,:,3)%v1)
    CALL sync_patch_array(sync_e,ptr_patch,ptr_patch%edges%dual_normal_vert(:,:,4)%v1)

    CALL sync_patch_array(sync_e,ptr_patch,ptr_patch%edges%primal_normal_cell(:,:,1)%v2)
    CALL sync_patch_array(sync_e,ptr_patch,ptr_patch%edges%primal_normal_cell(:,:,2)%v2)
    CALL sync_patch_array(sync_e,ptr_patch,ptr_patch%edges%primal_normal_vert(:,:,1)%v2)
    CALL sync_patch_array(sync_e,ptr_patch,ptr_patch%edges%primal_normal_vert(:,:,2)%v2)
    CALL sync_patch_array(sync_e,ptr_patch,ptr_patch%edges%primal_normal_vert(:,:,3)%v2)
    CALL sync_patch_array(sync_e,ptr_patch,ptr_patch%edges%primal_normal_vert(:,:,4)%v2)

    CALL sync_patch_array(sync_e,ptr_patch,ptr_patch%edges%dual_normal_cell(:,:,1)%v2)
    CALL sync_patch_array(sync_e,ptr_patch,ptr_patch%edges%dual_normal_cell(:,:,2)%v2)
    CALL sync_patch_array(sync_e,ptr_patch,ptr_patch%edges%dual_normal_vert(:,:,1)%v2)
    CALL sync_patch_array(sync_e,ptr_patch,ptr_patch%edges%dual_normal_vert(:,:,2)%v2)
    CALL sync_patch_array(sync_e,ptr_patch,ptr_patch%edges%dual_normal_vert(:,:,3)%v2)
    CALL sync_patch_array(sync_e,ptr_patch,ptr_patch%edges%dual_normal_vert(:,:,4)%v2)


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
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,ile1,ibe1) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
        & i_startidx, i_endidx, rl_start, rl_end)

      DO jc =  i_startidx, i_endidx

        IF(.NOT.ptr_patch%cells%decomp_info%owner_mask(jc,jb)) CYCLE

        DO je = 1, ptr_patch%cells%num_edges(jc,jb)

          ile1 = ptr_patch%cells%edge_idx(jc,jb,je)
          ibe1 = ptr_patch%cells%edge_blk(jc,jb,je)


          IF ((ptr_patch%edges%cell_idx(ile1,ibe1,1) == jc) .AND. &
            & (ptr_patch%edges%cell_blk(ile1,ibe1,1) == jb)) THEN

            ptr_int%primal_normal_ec(jc,jb,je,1) = &
              & ptr_patch%edges%primal_normal_cell(ile1,ibe1,1)%v1
            ptr_int%primal_normal_ec(jc,jb,je,2) = &
              & ptr_patch%edges%primal_normal_cell(ile1,ibe1,1)%v2
            ptr_int%edge_cell_length(jc,jb,je) = &
              & ptr_patch%edges%edge_cell_length(ile1,ibe1,1)

          ELSE IF ((ptr_patch%edges%cell_idx(ile1,ibe1,2) == jc) .AND. &
            & (ptr_patch%edges%cell_blk(ile1,ibe1,2) == jb)) THEN

            ptr_int%primal_normal_ec(jc,jb,je,1) = &
              & ptr_patch%edges%primal_normal_cell(ile1,ibe1,2)%v1
            ptr_int%primal_normal_ec(jc,jb,je,2) = &
              & ptr_patch%edges%primal_normal_cell(ile1,ibe1,2)%v2
            ptr_int%edge_cell_length(jc,jb,je) = &
              & ptr_patch%edges%edge_cell_length(ile1,ibe1,2)

          ENDIF

        ENDDO

      ENDDO

    END DO !block loop
!$OMP END DO

    rl_start = 1
    rl_end = min_rlcell

    ! Compute cell-vertex distances for gradient limiter
    i_startblk = ptr_patch%cells%start_blk(rl_start,1)
    i_endblk   = ptr_patch%cells%end_blk(rl_end,i_nchdom)
    !
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,cc_cell,cc_v1,cc_v2,cc_v3,cc_dis1,cc_dis2,cc_dis3, &
!$OMP            z_lon,z_lat,z_nx1,z_nx2,z_norm,ilv1,ibv1,ilv2,ibv2,ilv3,ibv3) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
        & i_startidx, i_endidx, rl_start, rl_end)

      DO jc =  i_startidx, i_endidx

        cc_cell = ptr_patch%cells%cartesian_center(jc,jb)
        z_lon   = ptr_patch%cells%center(jc,jb)%lon
        z_lat   = ptr_patch%cells%center(jc,jb)%lat

        CALL gvec2cvec(1._wp,0._wp,z_lon,z_lat,z_nx1(1),z_nx1(2),z_nx1(3),ptr_patch%geometry_info)
        z_norm = SQRT( DOT_PRODUCT(z_nx1(1:3),z_nx1(1:3)) )
        z_nx1(1:3)  = 1._wp/z_norm * z_nx1(1:3)

        CALL gvec2cvec(0._wp,1._wp,z_lon,z_lat,z_nx2(1),z_nx2(2),z_nx2(3),ptr_patch%geometry_info)
        z_norm = SQRT( DOT_PRODUCT(z_nx2(1:3),z_nx2(1:3)) )
        z_nx2(1:3)  = 1._wp/z_norm * z_nx2(1:3)

        ilv1 = ptr_patch%cells%vertex_idx(jc,jb,1)
        ibv1 = ptr_patch%cells%vertex_blk(jc,jb,1)
        ilv2 = ptr_patch%cells%vertex_idx(jc,jb,2)
        ibv2 = ptr_patch%cells%vertex_blk(jc,jb,2)
        ilv3 = ptr_patch%cells%vertex_idx(jc,jb,3)
        ibv3 = ptr_patch%cells%vertex_blk(jc,jb,3)

        cc_v1  =  ptr_patch%verts%cartesian(ilv1,ibv1)
        cc_v2  =  ptr_patch%verts%cartesian(ilv2,ibv2)
        cc_v3  =  ptr_patch%verts%cartesian(ilv3,ibv3)

        !Find closest coordinate for flat torus case: no need as
        !these points belong to 1 triangle so they can not be cyclic
        !IF(ptr_patch%geometry_info%geometry_type==planar_torus_geometry)THEN
        !  cc_v1  = plane_torus_closest_coordinates(cc_cell%x(1:3),cc_v1%x(1:3),ptr_patch%geometry_info)
        !  cc_v2  = plane_torus_closest_coordinates(cc_cell%x(1:3),cc_v2%x(1:3),ptr_patch%geometry_info)
        !  cc_v3  = plane_torus_closest_coordinates(cc_cell%x(1:3),cc_v3%x(1:3),ptr_patch%geometry_info)
        !END IF

        cc_dis1%x(1:3) = cc_v1%x(1:3) - cc_cell%x(1:3)
        z_norm = SQRT( DOT_PRODUCT(cc_dis1%x(1:3),cc_dis1%x(1:3)) )
        cc_dis1%x(1:3) = cc_dis1%x(1:3)/z_norm

        cc_dis2%x(1:3) = cc_v2%x(1:3) - cc_cell%x(1:3)
        z_norm = SQRT( DOT_PRODUCT(cc_dis2%x(1:3),cc_dis2%x(1:3)) )
        cc_dis2%x(1:3) = cc_dis2%x(1:3)/z_norm

        cc_dis3%x(1:3) = cc_v3%x(1:3) - cc_cell%x(1:3)
        z_norm = SQRT( DOT_PRODUCT(cc_dis3%x(1:3),cc_dis3%x(1:3)) )
        cc_dis3%x(1:3) = cc_dis3%x(1:3)/z_norm

        ptr_int%cell_vert_dist(jc,1,1,jb) = grid_sphere_radius* &
          & arc_length(cc_cell,cc_v1,ptr_patch%geometry_info)*DOT_PRODUCT(z_nx1(1:3),cc_dis1%x(1:3))
        ptr_int%cell_vert_dist(jc,2,1,jb) = grid_sphere_radius* &
          & arc_length(cc_cell,cc_v2,ptr_patch%geometry_info)*DOT_PRODUCT(z_nx1(1:3),cc_dis2%x(1:3))
        ptr_int%cell_vert_dist(jc,3,1,jb) = grid_sphere_radius* &
          & arc_length(cc_cell,cc_v3,ptr_patch%geometry_info)*DOT_PRODUCT(z_nx1(1:3),cc_dis3%x(1:3))

        ptr_int%cell_vert_dist(jc,1,2,jb) = grid_sphere_radius* &
          & arc_length(cc_cell,cc_v1,ptr_patch%geometry_info)*DOT_PRODUCT(z_nx2(1:3),cc_dis1%x(1:3))
        ptr_int%cell_vert_dist(jc,2,2,jb) = grid_sphere_radius* &
          & arc_length(cc_cell,cc_v2,ptr_patch%geometry_info)*DOT_PRODUCT(z_nx2(1:3),cc_dis2%x(1:3))
        ptr_int%cell_vert_dist(jc,3,2,jb) = grid_sphere_radius* &
          & arc_length(cc_cell,cc_v3,ptr_patch%geometry_info)*DOT_PRODUCT(z_nx2(1:3),cc_dis3%x(1:3))

      ENDDO

    END DO !block loop
!$OMP END DO NOWAIT

!$OMP END PARALLEL

    DO je = 1, ptr_patch%geometry_info%cell_type
      CALL sync_patch_array(sync_c,ptr_patch,ptr_int%primal_normal_ec(:,:,je,1))
      CALL sync_patch_array(sync_c,ptr_patch,ptr_int%primal_normal_ec(:,:,je,2))
      CALL sync_patch_array(sync_c,ptr_patch,ptr_int%edge_cell_length(:,:,je))
    ENDDO

  END SUBROUTINE complete_patchinfo
  !----------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! Calls routines to calculate coefficients "ptr_int%pos_on_tplane_e" for backward tracjectory
  !! calculations depending on grid geometry.
  !!
  SUBROUTINE init_tplane_e( ptr_patch, ptr_int )
    TYPE(t_patch),     INTENT(inout) :: ptr_patch
    TYPE(t_int_state), INTENT(inout) :: ptr_int

    CHARACTER(LEN=*), PARAMETER :: method_name = 'mo_intp_coeffs:init_tplane_e'

    !
    SELECT CASE(ptr_patch%geometry_info%geometry_type)

    CASE (planar_torus_geometry)
      CALL calculate_planar_distance_at_edge( ptr_patch, ptr_int )
      CALL calculate_dotproduct_at_edge(ptr_patch, ptr_int)

    CASE (sphere_geometry)
      CALL calculate_tangent_plane_at_edge( ptr_patch, ptr_int )
      CALL calculate_dotproduct_at_edge(ptr_patch, ptr_int)

    CASE DEFAULT
      CALL finish(method_name, "Undefined geometry type")

    END SELECT

  END SUBROUTINE init_tplane_e
  !-------------------------------------------------------------------------

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
  SUBROUTINE calculate_tangent_plane_at_edge (ptr_patch, ptr_int)

    TYPE(t_patch), INTENT(inout) :: ptr_patch  !< patch

    TYPE(t_int_state), INTENT(inout) :: ptr_int  !< interpolation state

    REAL(wp) ::                  &    !< geographical coordinates of edge midpoint
      & xyloc_edge(2)

    REAL(wp) ::                  &    !< geographical coordinates of neighbouring cell
      & xyloc_n1(2), xyloc_n2(2)     !< centers

    REAL(wp) ::                  &    !< coordinates of neighbouring cell centers on plane
      & xyloc_plane_n1(2), xyloc_plane_n2(2)

    REAL(wp) ::                  &    !< geographical coordinates of edge midpoints for the
      & xyloc_quad(4,2)              !< corresponding quadrilateral cell

    REAL(wp) ::                  &    !< coordinates of edge midpoints for the quadrilateral
      & xyloc_plane_quad(4,2)        !< cell on plane

    REAL(wp) ::                  &    !< geographical coordinates of edge vertices
      & xyloc_ve(2,2)

    REAL(wp) ::                  &    !< coordinates of edge vertices on plane
      & xyloc_plane_ve(2,2)

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

    CALL message('mo_intp_coeffs:calculate_tangent_plane_at_edge', '')

    i_rcstartlev = 2

    ! start and end block
    i_startblk = ptr_patch%edges%start_blk(i_rcstartlev,1)
    i_endblk   = ptr_patch%nblks_e


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
!$OMP            xyloc_plane_ve) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_e(ptr_patch, jb, i_startblk, i_endblk, &
        & i_startidx, i_endidx, i_rcstartlev)

      DO je = i_startidx, i_endidx

        IF(.NOT.ptr_patch%edges%decomp_info%owner_mask(je,jb)) CYCLE

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
          & xyloc_plane_n1(1), xyloc_plane_n1(2) )                   ! out


        ! get geographical coordinates of second cell center
        xyloc_n2(1)   = ptr_patch%cells%center(ilc2,ibc2)%lon
        xyloc_n2(2)   = ptr_patch%cells%center(ilc2,ibc2)%lat

        ! projection second cell center into local \lambda-\Phi-system
        CALL gnomonic_proj( xyloc_edge(1), xyloc_edge(2), xyloc_n2(1), xyloc_n2(2), &! in
          & xyloc_plane_n2(1), xyloc_plane_n2(2) )                   ! out


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
            & xyloc_quad(ne,2),                               &! in
            & xyloc_plane_quad(ne,1), xyloc_plane_quad(ne,2) ) ! out

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
            & xyloc_ve(nv,2),                               &! in
            & xyloc_plane_ve(nv,1), xyloc_plane_ve(nv,2)   ) ! out

        END DO


        !
        ! 4. rotate these vectors into a new local cartesian system. In this rotated
        !    system the coordinate axes point into the local normal and tangential
        !    direction at each edge.
        !

        ! centers
        !
        ptr_int%pos_on_tplane_e(je,jb,1,1) = grid_sphere_radius * (      &
          & xyloc_plane_n1(1)  * ptr_patch%edges%primal_normal(je,jb)%v1  &
          & + xyloc_plane_n1(2)  * ptr_patch%edges%primal_normal(je,jb)%v2 )

        ptr_int%pos_on_tplane_e(je,jb,1,2) = grid_sphere_radius * (      &
          & xyloc_plane_n1(1)  * ptr_patch%edges%dual_normal(je,jb)%v1    &
          & + xyloc_plane_n1(2)  * ptr_patch%edges%dual_normal(je,jb)%v2 )

        ptr_int%pos_on_tplane_e(je,jb,2,1) = grid_sphere_radius * (      &
          & xyloc_plane_n2(1)  * ptr_patch%edges%primal_normal(je,jb)%v1  &
          & + xyloc_plane_n2(2)  * ptr_patch%edges%primal_normal(je,jb)%v2 )

        ptr_int%pos_on_tplane_e(je,jb,2,2) = grid_sphere_radius * (      &
          & xyloc_plane_n2(1)  * ptr_patch%edges%dual_normal(je,jb)%v1    &
          & + xyloc_plane_n2(2)  * ptr_patch%edges%dual_normal(je,jb)%v2 )


        ! edges
        !
        DO ne = 1,4
          ptr_int%pos_on_tplane_e(je,jb,2+ne,1) = grid_sphere_radius * (        &
            & xyloc_plane_quad(ne,1)  * ptr_patch%edges%primal_normal(je,jb)%v1  &
            & + xyloc_plane_quad(ne,2)  * ptr_patch%edges%primal_normal(je,jb)%v2 )

          ptr_int%pos_on_tplane_e(je,jb,2+ne,2) = grid_sphere_radius * (        &
            & xyloc_plane_quad(ne,1)  * ptr_patch%edges%dual_normal(je,jb)%v1    &
            & + xyloc_plane_quad(ne,2)  * ptr_patch%edges%dual_normal(je,jb)%v2 )
        END DO


        ! vertices
        !
        DO nv = 1,2
          ptr_int%pos_on_tplane_e(je,jb,6+nv,1) = grid_sphere_radius * (     &
            & xyloc_plane_ve(nv,1)  * ptr_patch%edges%primal_normal(je,jb)%v1 &
            & + xyloc_plane_ve(nv,2)  * ptr_patch%edges%primal_normal(je,jb)%v2 )

          ptr_int%pos_on_tplane_e(je,jb,6+nv,2) = grid_sphere_radius * (     &
            & xyloc_plane_ve(nv,1)  * ptr_patch%edges%dual_normal(je,jb)%v1   &
            & + xyloc_plane_ve(nv,2)  * ptr_patch%edges%dual_normal(je,jb)%v2 )
        END DO

      ENDDO ! edges
    ENDDO  ! blocks
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    DO ne=1,8
      CALL sync_patch_array(sync_e, ptr_patch, ptr_int%pos_on_tplane_e(:,:,ne,1))
      CALL sync_patch_array(sync_e, ptr_patch, ptr_int%pos_on_tplane_e(:,:,ne,2))
    ENDDO

  END SUBROUTINE calculate_tangent_plane_at_edge
  !----------------------------------------------------------------------------

  !----------------------------------------------------------------------------
  ! @AD: This part was initially in the init_tplane_e but now separated into two
  !
  ! For each of the 4 rhomboidal edges transform normal and tangential
  ! unit vectors into cartesian system. Then compute dot product
  ! between these unit vectors and the unit vectors of the inner edge.
  ! normalization not necessary fo cartesian vectors since these are
  ! exactly =1.
  !
  SUBROUTINE calculate_dotproduct_at_edge (ptr_patch, ptr_int)

    TYPE(t_patch), INTENT(inout) :: ptr_patch  !< patch

    TYPE(t_int_state), INTENT(inout) :: ptr_int  !< interpolation state

    REAL(wp) ::                  &    !< primal/dual normal in cartesian coordinates
      & z_nx(3), z_ny(3)

    REAL(wp) :: z_nx_quad(3),    &    !< primal/dual normal at quadrilateral
      & z_ny_quad(3)          !< edges in cartesian coordinates

    !< cell centers
    INTEGER :: ilq, ibq               !< line and block indices of quadrilateral edges
    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER :: i_rcstartlev
    INTEGER :: jb, je                 !< loop indices for block and edges
    INTEGER :: ne                     !< loop index for quadrilateral edges
    !-------------------------------------------------------------------------

    CALL message('mo_intp_coeffs:calculate_dotproduct_at_edge', '')

    i_rcstartlev = 2

    ! start and end block
    i_startblk = ptr_patch%edges%start_blk(i_rcstartlev,1)
    i_endblk   = ptr_patch%nblks_e

!$OMP PARALLEL
!$OMP DO PRIVATE(je,jb,ne,ilq,ibq,i_startidx,i_endidx,z_nx,z_ny,&
!$OMP z_nx_quad,z_ny_quad) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_e(ptr_patch, jb, i_startblk, i_endblk, &
        & i_startidx, i_endidx, i_rcstartlev)

      DO je = i_startidx, i_endidx

        IF(.NOT.ptr_patch%edges%decomp_info%owner_mask(je,jb)) CYCLE

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
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    DO ne=1,4
      CALL sync_patch_array(sync_e, ptr_patch, ptr_int%tplane_e_dotprod(:,:,ne,1))
      CALL sync_patch_array(sync_e, ptr_patch, ptr_int%tplane_e_dotprod(:,:,ne,2))
      CALL sync_patch_array(sync_e, ptr_patch, ptr_int%tplane_e_dotprod(:,:,ne,3))
      CALL sync_patch_array(sync_e, ptr_patch, ptr_int%tplane_e_dotprod(:,:,ne,4))
    ENDDO


  END SUBROUTINE calculate_dotproduct_at_edge
  !----------------------------------------------------------------------------

  !----------------------------------------------------------------------------
  !>
  !! @par Revision History
  !! Taken from "calculate_tangent_plane_at_edge" for flat torus case
  !!
  SUBROUTINE calculate_planar_distance_at_edge (ptr_patch, ptr_int)

    TYPE(t_patch), INTENT(inout) :: ptr_patch  !< patch

    TYPE(t_int_state), INTENT(inout) :: ptr_int  !< interpolation state

    !CC of points on the plane torus grid
    TYPE(t_cartesian_coordinates) :: cc_n1, cc_n2, cc_quad(4), cc_ve(2), cc_edge

    !relative location of those points w.r.t cc_edge
    TYPE(t_cartesian_coordinates) :: cc_plane_n1, cc_plane_n2, cc_plane_quad(4), &
                                     cc_plane_ve(2)

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

    CALL message('mo_intp_coeffs:calculate_planar_distance_at_edge', '')

    i_rcstartlev = 2

    ! start and end block
    i_startblk = ptr_patch%edges%start_blk(i_rcstartlev,1)
    i_endblk   = ptr_patch%nblks_e

    !<< AD
    ! Modification for is_plane_torus for HDCP2: the planar distance between any two
    ! points is same as the Gnomonic projected distance
    ! AD >>

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,je,ne,nv,ilc1,ibc1,ilc2,ibc2,ilq,ibq,ilv,ibv,i_startidx, &
!$OMP            i_endidx,cc_edge,cc_n1,cc_n2,cc_plane_n1,       &
!$OMP            cc_plane_n2,cc_quad,cc_plane_quad,cc_ve,        &
!$OMP            cc_plane_ve) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_e(ptr_patch, jb, i_startblk, i_endblk, &
                        i_startidx, i_endidx, i_rcstartlev)

      DO je = i_startidx, i_endidx

        IF(.NOT.ptr_patch%edges%decomp_info%owner_mask(je,jb)) CYCLE

        ! 1. CC of neighboring cell centers; z index is 0
        cc_edge = ptr_patch%edges%cartesian_center(je,jb)

        ! get line and block indices of neighbour cells
        ilc1 = ptr_patch%edges%cell_idx(je,jb,1)
        ibc1 = ptr_patch%edges%cell_blk(je,jb,1)
        ilc2 = ptr_patch%edges%cell_idx(je,jb,2)
        ibc2 = ptr_patch%edges%cell_blk(je,jb,2)

        ! get CC of the cell centers
        cc_n1 = ptr_patch%cells%cartesian_center(ilc1,ibc1)
        cc_n2 = ptr_patch%cells%cartesian_center(ilc2,ibc2)

        !now calculate the separation vector between the edge and cell centers
        cc_n1 =  plane_torus_closest_coordinates(cc_edge%x,cc_n1%x,ptr_patch%geometry_info)
        cc_plane_n1%x(:) = cc_n1%x(:) - cc_edge%x(:)

        cc_n2 =  plane_torus_closest_coordinates(cc_edge%x,cc_n2%x,ptr_patch%geometry_info)
        cc_plane_n2%x(:) = cc_n2%x(:) - cc_edge%x(:)

        ! 2. Edge midpoints of the quadrilateral
        DO ne = 1,4

          ! get line and block indices of edge midpoints
          ilq = ptr_patch%edges%quad_idx(je,jb,ne)
          ibq = ptr_patch%edges%quad_blk(je,jb,ne)

          ! get CC coordinates of edge midpoints
          cc_quad(ne) = ptr_patch%edges%cartesian_center(ilq,ibq)

          !now calculate the separation vector between the edge and cell centers
          cc_quad(ne) =  plane_torus_closest_coordinates(cc_edge%x,cc_quad(ne)%x,ptr_patch%geometry_info)
          cc_plane_quad(ne)%x(:) = cc_quad(ne)%x(:) - cc_edge%x(:)

        END DO
        !
        ! 3. Edge vertices
        !
        DO nv = 1,2

          ! get line and block indices of edge vertices
          ilv = ptr_patch%edges%vertex_idx(je,jb,nv)
          ibv = ptr_patch%edges%vertex_blk(je,jb,nv)

          ! get CC coordinates of edge vertices
          cc_ve(nv) = ptr_patch%verts%cartesian(ilv,ibv)

          !now calculate the separation vector between the edge and cell centers
          cc_ve(nv) =  plane_torus_closest_coordinates(cc_edge%x,cc_ve(nv)%x,ptr_patch%geometry_info)
          cc_plane_ve(nv)%x(:) = cc_ve(nv)%x(:) - cc_edge%x(:)

        END DO

        !
        ! 4. rotate these vectors into a new local cartesian system. In this rotated
        !    system the coordinate axes point into the local normal and tangential
        !    direction at each edge. All these vectors are 2D so no point using the
        !    z coordinate
        !

        ! centers
        !
        ptr_int%pos_on_tplane_e(je,jb,1,1) = SUM(      &
          &  cc_plane_n1%x(1:2)  * ptr_patch%edges%primal_cart_normal(je,jb)%x(1:2) )

        ptr_int%pos_on_tplane_e(je,jb,1,2) = SUM(      &
          &  cc_plane_n1%x(1:2)  * ptr_patch%edges%dual_cart_normal(je,jb)%x(1:2) )

        ptr_int%pos_on_tplane_e(je,jb,2,1) = SUM(      &
          &  cc_plane_n2%x(1:2)  * ptr_patch%edges%primal_cart_normal(je,jb)%x(1:2) )

        ptr_int%pos_on_tplane_e(je,jb,2,2) = SUM(      &
          &  cc_plane_n2%x(1:2)  * ptr_patch%edges%dual_cart_normal(je,jb)%x(1:2) )

        ! edges
        !
        DO ne = 1,4
          ptr_int%pos_on_tplane_e(je,jb,2+ne,1) = SUM(      &
            &  cc_plane_quad(ne)%x(1:2)  * ptr_patch%edges%primal_cart_normal(je,jb)%x(1:2) )

          ptr_int%pos_on_tplane_e(je,jb,2+ne,2) = SUM(      &
            &  cc_plane_quad(ne)%x(1:2)  * ptr_patch%edges%dual_cart_normal(je,jb)%x(1:2) )
        END DO

        ! vertices
        !
        DO nv = 1,2
          ptr_int%pos_on_tplane_e(je,jb,6+nv,1) = SUM(      &
            &  cc_plane_ve(nv)%x(1:2)  * ptr_patch%edges%primal_cart_normal(je,jb)%x(1:2) )

          ptr_int%pos_on_tplane_e(je,jb,6+nv,2) = SUM(      &
            &  cc_plane_ve(nv)%x(1:2)  * ptr_patch%edges%dual_cart_normal(je,jb)%x(1:2) )
        END DO

      ENDDO ! edges
    ENDDO  ! blocks
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    DO ne=1,8
      call sync_patch_array(SYNC_E, ptr_patch, ptr_int%pos_on_tplane_e(:,:,ne,1))
      call sync_patch_array(SYNC_E, ptr_patch, ptr_int%pos_on_tplane_e(:,:,ne,2))
    ENDDO


  END SUBROUTINE calculate_planar_distance_at_edge
  !----------------------------------------------------------------------------

  !----------------------------------------------------------------------------
  !>
  !! Calls routines to calculate coefficients "ptr_int%pos_on_tplane_c_edge".
  !! Bifurcation for calculations depending on grid geometry.
  !!
  SUBROUTINE init_tplane_c( ptr_patch, ptr_int )
    TYPE(t_patch),     INTENT(inout) :: ptr_patch
    TYPE(t_int_state), INTENT(inout) :: ptr_int

    CHARACTER(LEN=*), PARAMETER :: method_name = 'mo_intp_coeffs:init_tplane_c'

    !
    SELECT CASE(ptr_patch%geometry_info%geometry_type)

    CASE (planar_torus_geometry)
      CALL init_tplane_c_torus ( ptr_patch, ptr_int )
    CASE (sphere_geometry)
      CALL init_tplane_c_sphere ( ptr_patch, ptr_int )
    CASE DEFAULT
      CALL finish(method_name, "Undefined geometry type")
    END SELECT

  END SUBROUTINE init_tplane_c


  !----------------------------------------------------------------------------
  !>
  !!AD> This is old "init_tplane_c" routine now renamed!
  !!
  !! Initializes a tangential plane at each cell circumcenter. Necessary for efficient
  !! computation of flux areas and overlap regions between flux areas and the
  !! model grid.
  !!
  !! The position of cell vertices is precomputed using the gnomonic projection.
  !! These vertices are then stored in an edge-based structure. Moreover, in order
  !! to avoid inconsistencies, an additional "projection error" is added when
  !! projecting the vertices. Instead of performing the projection on a plane
  !! tangent to the cell center (cell-based system), the projection is first
  !! performed on a plane tangent to the cell edge (edge-based system). Then, the
  !! projected points are transformed back into a cell based system, assuming
  !! co-planarity of the two systems.
  !!
  !! Order of storage for pos_on_tplane_c_edge:
  !! pos_on_tplane_c_edge(nproma,nblks_e,ncells=2,npts=5)%lon/lat
  !! - cell ordering according to edge%cell_idx/blk
  !! - npts 1-3: cell vertices
  !!   - ordering of first 2 vertices according to edge%vertex_idx/blk
  !!   - ordering of coordinates: 1=x, 2=y
  !! - npts 4-5: coordinates of neighboring cell centers (share vertex 1 and 2, respectively)
  !!   - only those 2 neighbors that do not share the given edge.
  !!   - no "projection error" added to cell center
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2012-03-12)
  !! Modification by Daniel Reinert, DWD, (2012-04-12):
  !! - added projection of cell centers
  !!
  SUBROUTINE init_tplane_c_sphere (ptr_patch, ptr_int)

    TYPE(t_patch),     INTENT(inout) :: ptr_patch  !< patch

    TYPE(t_int_state), INTENT(inout) :: ptr_int    !< interpolation state

    REAL(wp) ::   &              !< geographical coords of edge midpoint
      & xyloc_edge(2)

    REAL(wp) ::   &              !< geographical coords of edge-vertices
      & xyloc_v(4,2)

    REAL(wp) ::   &              !< coords of vertices when projected onto
      & xyloc_plane_v(4,2)      !< edge-based plane
    REAL(wp) ::   &              !< same, but for rotated system (normal-tangential)
      & xyloc_plane_nt_v(4,2)

    REAL(wp) ::   &              !< geographical coords of neighboring cell centers
      & xyloc_n1(2), xyloc_n2(2)

    REAL(wp) ::   &              !< coords of cell centers when projected onto
      & xyloc_plane_n1(2), xyloc_plane_n2(2) !< edge-based plane
    REAL(wp) ::   &              !< same but for rotated system (normal-tangential)
      & xyloc_plane_nt_n1(2), xyloc_plane_nt_n2(2)

    REAL(wp) ::   &              !< coords of edge-vertices in translated system
      & xyloc_trans1_v(4,2), xyloc_trans2_v(4,2)

    REAL(wp) ::   &              !< primal and dual normals for neighboring cells
      & pn_cell1(2), pn_cell2(2), dn_cell1(2), dn_cell2(2)

    REAL(wp) ::   &              !< geographical coords of butterfly neighbors
      & xyloc_bf1(2,2),       & !< for edge-neighbor 1 and 2
      & xyloc_bf2(2,2)

    REAL(wp) ::   &              !< coords of butterfly neighbors after projection
      & xyloc_plane_bf1(2,2), & !< onto cell-based plane
      & xyloc_plane_bf2(2,2)

    INTEGER :: ilv(4), ibv(4)
    INTEGER :: ilc1, ilc2, ibc1, ibc2
    INTEGER :: ilc_bf1(2), ilc_bf2(2), ibc_bf1(2), ibc_bf2(2)
    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER :: i_rcstartlev
    INTEGER :: nv, nc            !< loop indices for vertices and cells
    INTEGER :: jb, je            !< loop indices for block and edges

    !-------------------------------------------------------------------------

    CALL message('mo_interpolation:init_tplane_c_sphere','')

    i_rcstartlev = 2

    ! start and end block
    i_startblk = ptr_patch%edges%start_blk(i_rcstartlev,1)
    i_endblk   = ptr_patch%nblks_e



!$OMP PARALLEL
!$OMP DO PRIVATE(jb,je,nv,i_startidx,i_endidx,xyloc_edge,ilv,ibv,          &
!$OMP            xyloc_v,xyloc_plane_v,ilc1,ilc2,ibc1,ibc2,xyloc_n1,       &
!$OMP            xyloc_n2,xyloc_plane_n1,xyloc_plane_n2,xyloc_plane_nt_n1, &
!$OMP            xyloc_plane_nt_n2,xyloc_plane_nt_v,xyloc_trans1_v,        &
!$OMP            xyloc_trans2_v,pn_cell1,pn_cell2,dn_cell1,dn_cell2)
    DO jb = i_startblk, i_endblk

      CALL get_indices_e(ptr_patch, jb, i_startblk, i_endblk, &
        & i_startidx, i_endidx, i_rcstartlev)

      DO je = i_startidx, i_endidx

        IF(.NOT.ptr_patch%edges%decomp_info%owner_mask(je,jb)) CYCLE

        !
        ! 1. project edge-vertices and the centers of the neighboring cells
        !    onto a plane tangent to the edge midpoint.
        !

        ! get geographical coordinates of edge midpoint
        xyloc_edge(1) = ptr_patch%edges%center(je,jb)%lon
        xyloc_edge(2) = ptr_patch%edges%center(je,jb)%lat

        ! get line and block indices of edge vertices (including the
        ! non-edge-aligned vertices of the neighboring cells
        ilv(1:4)=ptr_patch%edges%vertex_idx(je,jb,1:4)
        ibv(1:4)=ptr_patch%edges%vertex_blk(je,jb,1:4)

        ! get geographical coordinates of edge vertices
        DO nv=1,4
          xyloc_v(nv,1)=ptr_patch%verts%vertex(ilv(nv),ibv(nv))%lon
          xyloc_v(nv,2)=ptr_patch%verts%vertex(ilv(nv),ibv(nv))%lat
        ENDDO


        ! project vertices into edge-based local system
        CALL gnomonic_proj( xyloc_edge(1), xyloc_edge(2), xyloc_v(1,1), xyloc_v(1,2), &! in
          & xyloc_plane_v(1,1), xyloc_plane_v(1,2) )                   ! out

        CALL gnomonic_proj( xyloc_edge(1), xyloc_edge(2), xyloc_v(2,1), xyloc_v(2,2), &! in
          & xyloc_plane_v(2,1), xyloc_plane_v(2,2) )                   ! out

        CALL gnomonic_proj( xyloc_edge(1), xyloc_edge(2), xyloc_v(3,1), xyloc_v(3,2), &! in
          & xyloc_plane_v(3,1), xyloc_plane_v(3,2) )                   ! out

        CALL gnomonic_proj( xyloc_edge(1), xyloc_edge(2), xyloc_v(4,1), xyloc_v(4,2), &! in
          & xyloc_plane_v(4,1), xyloc_plane_v(4,2) )                   ! out


        ! get line and block indices of neighbour cells
        ilc1 = ptr_patch%edges%cell_idx(je,jb,1)
        ibc1 = ptr_patch%edges%cell_blk(je,jb,1)
        ilc2 = ptr_patch%edges%cell_idx(je,jb,2)
        ibc2 = ptr_patch%edges%cell_blk(je,jb,2)

        ! get geographical coordinates of first cell center
        xyloc_n1(1)   = ptr_patch%cells%center(ilc1,ibc1)%lon
        xyloc_n1(2)   = ptr_patch%cells%center(ilc1,ibc1)%lat

        ! get geographical coordinates of second cell center
        xyloc_n2(1)   = ptr_patch%cells%center(ilc2,ibc2)%lon
        xyloc_n2(2)   = ptr_patch%cells%center(ilc2,ibc2)%lat

        ! project first cell center into edge-based local system
        CALL gnomonic_proj( xyloc_edge(1), xyloc_edge(2), xyloc_n1(1), xyloc_n1(2), &! in
          & xyloc_plane_n1(1), xyloc_plane_n1(2) )                   ! out

        ! project second cell center into edge-based local system local
        CALL gnomonic_proj( xyloc_edge(1), xyloc_edge(2), xyloc_n2(1), xyloc_n2(2), &! in
          & xyloc_plane_n2(1), xyloc_plane_n2(2) )                   ! out


        !
        ! 2. rotate these vectors into a new local cartesian system. In this rotated
        !    system the coordinate axes point into the local normal and tangential
        !    direction at each edge.
        !

        ! centers
        !
        xyloc_plane_nt_n1(1) =                                                &
          & xyloc_plane_n1(1)  * ptr_patch%edges%primal_normal(je,jb)%v1  &
          & + xyloc_plane_n1(2)  * ptr_patch%edges%primal_normal(je,jb)%v2

        xyloc_plane_nt_n1(2) =                                                &
          & xyloc_plane_n1(1)  * ptr_patch%edges%dual_normal(je,jb)%v1    &
          & + xyloc_plane_n1(2)  * ptr_patch%edges%dual_normal(je,jb)%v2

        xyloc_plane_nt_n2(1) =                                                &
          & xyloc_plane_n2(1)  * ptr_patch%edges%primal_normal(je,jb)%v1  &
          & + xyloc_plane_n2(2)  * ptr_patch%edges%primal_normal(je,jb)%v2

        xyloc_plane_nt_n2(2) =                                                &
          & xyloc_plane_n2(1)  * ptr_patch%edges%dual_normal(je,jb)%v1    &
          & + xyloc_plane_n2(2)  * ptr_patch%edges%dual_normal(je,jb)%v2

        ! vertices
        !
        DO nv = 1,4
          xyloc_plane_nt_v(nv,1) =                                              &
            & xyloc_plane_v(nv,1) * ptr_patch%edges%primal_normal(je,jb)%v1 &
            & + xyloc_plane_v(nv,2) * ptr_patch%edges%primal_normal(je,jb)%v2

          xyloc_plane_nt_v(nv,2) =                                              &
            & xyloc_plane_v(nv,1) * ptr_patch%edges%dual_normal(je,jb)%v1   &
            & + xyloc_plane_v(nv,2) * ptr_patch%edges%dual_normal(je,jb)%v2
        END DO



        ! 3. Calculate position of vertices in a translated coordinate system.
        !    This is done twice. The origin is located once at the circumcenter
        !    of the neighboring cell 1 and once cell 2. The distance vectors point
        !    from the cell center to the vertices.
        xyloc_trans1_v(1:4,1) = xyloc_plane_nt_v(1:4,1) - xyloc_plane_nt_n1(1)
        xyloc_trans1_v(1:4,2) = xyloc_plane_nt_v(1:4,2) - xyloc_plane_nt_n1(2)

        xyloc_trans2_v(1:4,1) = xyloc_plane_nt_v(1:4,1) - xyloc_plane_nt_n2(1)
        xyloc_trans2_v(1:4,2) = xyloc_plane_nt_v(1:4,2) - xyloc_plane_nt_n2(2)



        ! 4. Rotate points into coordinate system pointing into local north
        !    and local east direction. Store in edge-based data structure. This
        !    is done twice (for both neighboring cells).
        ! e_n= pn_cell1(1)*e_\lambda + pn_cell1(2)*e_\phi
        ! e_t= dn_cell1(1)*e_\lambda + dn_cell1(2)*e_\phi
        pn_cell1(1) = ptr_patch%edges%primal_normal_cell(je,jb,1)%v1
        pn_cell1(2) = ptr_patch%edges%primal_normal_cell(je,jb,1)%v2
        dn_cell1(1) = ptr_patch%edges%dual_normal_cell(je,jb,1)%v1
        dn_cell1(2) = ptr_patch%edges%dual_normal_cell(je,jb,1)%v2

        pn_cell2(1) = ptr_patch%edges%primal_normal_cell(je,jb,2)%v1
        pn_cell2(2) = ptr_patch%edges%primal_normal_cell(je,jb,2)%v2
        dn_cell2(1) = ptr_patch%edges%dual_normal_cell(je,jb,2)%v1
        dn_cell2(2) = ptr_patch%edges%dual_normal_cell(je,jb,2)%v2

        ! components in longitudinal direction (cell 1)
        ptr_int%pos_on_tplane_c_edge(je,jb,1,1:3)%lon =  grid_sphere_radius          &
          & *( xyloc_trans1_v(1:3,1) * pn_cell1(1)       &
          & +  xyloc_trans1_v(1:3,2) * dn_cell1(1) )

        ! components in latitudinal direction (cell 1)
        ptr_int%pos_on_tplane_c_edge(je,jb,1,1:3)%lat =  grid_sphere_radius          &
          & *( xyloc_trans1_v(1:3,1) * pn_cell1(2)       &
          & +  xyloc_trans1_v(1:3,2) * dn_cell1(2) )

        ! components in longitudinal direction (cell 2)
        ptr_int%pos_on_tplane_c_edge(je,jb,2,1:3)%lon =  grid_sphere_radius          &
          & *( xyloc_trans2_v((/1,2,4/),1) * pn_cell2(1) &
          & +  xyloc_trans2_v((/1,2,4/),2) * dn_cell2(1) )

        ! components in latitudinal direction (cell 2)
        ptr_int%pos_on_tplane_c_edge(je,jb,2,1:3)%lat =  grid_sphere_radius          &
          & *( xyloc_trans2_v((/1,2,4/),1) * pn_cell2(2) &
          & +  xyloc_trans2_v((/1,2,4/),2) * dn_cell2(2) )

      ENDDO  ! je
    ENDDO  ! jb
!$OMP END DO
!$OMP END PARALLEL



    !
    ! Projection of butterfly cell centers
    !

    i_rcstartlev = 3

    ! start and end block
    i_startblk = ptr_patch%edges%start_blk(i_rcstartlev,1)
    i_endblk   = ptr_patch%nblks_e

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,je,i_startidx,i_endidx,ilc1,ilc2,ibc1,ibc2,xyloc_n1,&
!$OMP            xyloc_n2,ilc_bf1,ibc_bf1,ilc_bf2,ibc_bf2,xyloc_bf1,    &
!$OMP            xyloc_bf2,xyloc_plane_bf1,xyloc_plane_bf2)
    DO jb = i_startblk, i_endblk

      CALL get_indices_e(ptr_patch, jb, i_startblk, i_endblk, &
        & i_startidx, i_endidx, i_rcstartlev)

      DO je = i_startidx, i_endidx

        IF(.NOT.ptr_patch%edges%decomp_info%owner_mask(je,jb)) CYCLE


        ! get line and block indices of edge-neighbours (cells)
        ilc1 = ptr_patch%edges%cell_idx(je,jb,1)
        ibc1 = ptr_patch%edges%cell_blk(je,jb,1)
        ilc2 = ptr_patch%edges%cell_idx(je,jb,2)
        ibc2 = ptr_patch%edges%cell_blk(je,jb,2)

        ! get geographical coordinates of edge-neighbor 1
        xyloc_n1(1)   = ptr_patch%cells%center(ilc1,ibc1)%lon
        xyloc_n1(2)   = ptr_patch%cells%center(ilc1,ibc1)%lat

        ! get geographical coordinates of edge-neighbor 2
        xyloc_n2(1)   = ptr_patch%cells%center(ilc2,ibc2)%lon
        xyloc_n2(2)   = ptr_patch%cells%center(ilc2,ibc2)%lat



        ! 5a. project cell centers adjacent to the cells sharing edge je.
        !    (from butterfly_idx)


        ! get line and block indices of the neighbors of edge-neighbor 1
        ilc_bf1(1:2) = ptr_patch%edges%butterfly_idx(je,jb,1,1:2)
        ibc_bf1(1:2) = ptr_patch%edges%butterfly_blk(je,jb,1,1:2)

        ! get line and block indices of the neighbors of edge-neighbor 2
        ilc_bf2(1:2) = ptr_patch%edges%butterfly_idx(je,jb,2,1:2)
        ibc_bf2(1:2) = ptr_patch%edges%butterfly_blk(je,jb,2,1:2)


        ! get geographical coordinates of cell centers (neighbor 1)
        xyloc_bf1(1,1) = ptr_patch%cells%center(ilc_bf1(1),ibc_bf1(1))%lon
        xyloc_bf1(1,2) = ptr_patch%cells%center(ilc_bf1(1),ibc_bf1(1))%lat
        xyloc_bf1(2,1) = ptr_patch%cells%center(ilc_bf1(2),ibc_bf1(2))%lon
        xyloc_bf1(2,2) = ptr_patch%cells%center(ilc_bf1(2),ibc_bf1(2))%lat

        ! get geographical coordinates of cell centers (neighbor 2)
        xyloc_bf2(1,1) = ptr_patch%cells%center(ilc_bf2(1),ibc_bf2(1))%lon
        xyloc_bf2(1,2) = ptr_patch%cells%center(ilc_bf2(1),ibc_bf2(1))%lat
        xyloc_bf2(2,1) = ptr_patch%cells%center(ilc_bf2(2),ibc_bf2(2))%lon
        xyloc_bf2(2,2) = ptr_patch%cells%center(ilc_bf2(2),ibc_bf2(2))%lat


        ! project cell centers into cell-based local system (edge-neighbor 1)
        CALL gnomonic_proj( xyloc_n1(1), xyloc_n1(2), xyloc_bf1(1,1), xyloc_bf1(1,2), &! in
          & xyloc_plane_bf1(1,1), xyloc_plane_bf1(1,2) )               ! out

        CALL gnomonic_proj( xyloc_n1(1), xyloc_n1(2), xyloc_bf1(2,1), xyloc_bf1(2,2), &! in
          & xyloc_plane_bf1(2,1), xyloc_plane_bf1(2,2) )               ! out



        ! project cell centers into cell-based local system (edge-neighbor 2)
        CALL gnomonic_proj( xyloc_n2(1), xyloc_n2(2), xyloc_bf2(1,1), xyloc_bf2(1,2), &! in
          & xyloc_plane_bf2(1,1), xyloc_plane_bf2(1,2) )               ! out

        CALL gnomonic_proj( xyloc_n2(1), xyloc_n2(2), xyloc_bf2(2,1), xyloc_bf2(2,2), &! in
          & xyloc_plane_bf2(2,1), xyloc_plane_bf2(2,2) )               ! out


        ! components in longitudinal direction (cell 1)
        ptr_int%pos_on_tplane_c_edge(je,jb,1,4:5)%lon =  grid_sphere_radius * &
          & xyloc_plane_bf1(1:2,1)

        ! components in latitudinal direction (cell 1)
        ptr_int%pos_on_tplane_c_edge(je,jb,1,4:5)%lat =  grid_sphere_radius * &
          & xyloc_plane_bf1(1:2,2)


        ! components in longitudinal direction (cell 2)
        ptr_int%pos_on_tplane_c_edge(je,jb,2,4:5)%lon =  grid_sphere_radius * &
          & xyloc_plane_bf2(1:2,1)

        ! components in latitudinal direction (cell 2)
        ptr_int%pos_on_tplane_c_edge(je,jb,2,4:5)%lat =  grid_sphere_radius * &
          & xyloc_plane_bf2(1:2,2)

      ENDDO  ! je
    ENDDO  ! jb
!$OMP END DO
!$OMP END PARALLEL


    DO nv=1,5
      DO nc=1,2
        CALL sync_patch_array(sync_e, ptr_patch, ptr_int%pos_on_tplane_c_edge(:,:,nc,nv)%lon)
        CALL sync_patch_array(sync_e, ptr_patch, ptr_int%pos_on_tplane_c_edge(:,:,nc,nv)%lat)
      ENDDO
    ENDDO


  END SUBROUTINE init_tplane_c_sphere

  !----------------------------------------------------------------------------

  !----------------------------------------------------------------------------
  !>
  !!AD> This is torus version of "init_tplane_c_sphere"- everything remains same
  !!    just that for flat torus geometry there is no need of gnomonic projection.
  !!    Radial vector between two point is just coordinate difference
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2012-03-12)
  !! Modified by Anurag Dipankar, MPIM (2013-04) for torus geometry

  SUBROUTINE init_tplane_c_torus (ptr_patch, ptr_int)

    TYPE(t_patch),     INTENT(inout) :: ptr_patch  !< patch

    TYPE(t_int_state), INTENT(inout) :: ptr_int    !< interpolation state

    !CC of points in the stencil
    TYPE(t_cartesian_coordinates) :: cc_edge, cc_vert(4), cc_cell(2), &
                                     cc_bf1(2), cc_bf2(2)

    !Distance vectors
    TYPE(t_cartesian_coordinates) :: dist_edge_cell(2), dist_edge_vert(4), &
                                     dist_bf1_cell1(2), dist_bf2_cell2(2)

    REAL(wp) ::   &              !< same, but for rotated system (normal-tangential)
      & xyloc_plane_nt_v(4,2)

    REAL(wp) ::   &              !< same but for rotated system (normal-tangential)
      & xyloc_plane_nt_n1(2), xyloc_plane_nt_n2(2)

    REAL(wp) ::   &              !< coords of edge-vertices in translated system
      & xyloc_trans1_v(4,2), xyloc_trans2_v(4,2)

    REAL(wp) ::   &              !< primal and dual normals for neighboring cells
      & pn_cell1(2), pn_cell2(2), dn_cell1(2), dn_cell2(2)

    INTEGER :: ilv(4), ibv(4)
    INTEGER :: ilc1, ilc2, ibc1, ibc2
    INTEGER :: ilc_bf1(2), ilc_bf2(2), ibc_bf1(2), ibc_bf2(2)
    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER :: i_rcstartlev
    INTEGER :: nv, nc            !< loop indices for vertices and cells
    INTEGER :: jb, je            !< loop indices for block and edges

    !-------------------------------------------------------------------------

    CALL message('mo_interpolation:', 'init_tplane_c_torus')

    i_rcstartlev = 2

    ! start and end block
    i_startblk = ptr_patch%edges%start_blk(i_rcstartlev,1)
    i_endblk   = ptr_patch%nblks_e

    !<< AD
    ! Modification for is_plane_torus for HDCP2: the planar distance between any two
    ! points is same as the Gnomonic projected distance
    ! AD >>

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,je,nv,i_startidx,i_endidx,cc_edge,ilv,ibv,          &
!$OMP            cc_vert,dist_edge_vert,ilc1,ilc2,ibc1,ibc2,cc_cell,    &
!$OMP            dist_edge_cell,xyloc_plane_nt_n1,xyloc_plane_nt_n2,    &
!$OMP            xyloc_plane_nt_v,xyloc_trans1_v,        &
!$OMP            xyloc_trans2_v,pn_cell1,pn_cell2,dn_cell1,dn_cell2)
    DO jb = i_startblk, i_endblk

      CALL get_indices_e(ptr_patch, jb, i_startblk, i_endblk, &
        & i_startidx, i_endidx, i_rcstartlev)

      DO je = i_startidx, i_endidx

        IF(.NOT.ptr_patch%edges%decomp_info%owner_mask(je,jb)) CYCLE

        !
        ! 1. Get distance vector of edge-vertices and the centers of the neighboring cells
        !    w.r.t to the edge midpoint.
        !

        ! get CC of edge midpoint
        cc_edge = ptr_patch%edges%cartesian_center(je,jb)

        ! get line and block indices of edge vertices (including the
        ! non-edge-aligned vertices of the neighboring cells
        ilv(1:4)=ptr_patch%edges%vertex_idx(je,jb,1:4)
        ibv(1:4)=ptr_patch%edges%vertex_blk(je,jb,1:4)

        ! get CC of edge vertices
        DO nv=1,4
          cc_vert(nv) = ptr_patch%verts%cartesian(ilv(nv),ibv(nv))

          !Get the actual location of the vertex w.r.t the edge cc
          cc_vert(nv) = plane_torus_closest_coordinates(cc_edge%x,cc_vert(nv)%x, &
                                                         ptr_patch%geometry_info)
          !Get the distance vector between cc_edge and cc_vets
          dist_edge_vert(nv)%x(:) = cc_vert(nv)%x(:) - cc_edge%x(:)
        ENDDO

        ! get line and block indices of neighbour cells
        ilc1 = ptr_patch%edges%cell_idx(je,jb,1)
        ibc1 = ptr_patch%edges%cell_blk(je,jb,1)
        ilc2 = ptr_patch%edges%cell_idx(je,jb,2)
        ibc2 = ptr_patch%edges%cell_blk(je,jb,2)

        ! get CC of the two cell centers
        cc_cell(1)   = ptr_patch%cells%cartesian_center(ilc1,ibc1)
        cc_cell(2)   = ptr_patch%cells%cartesian_center(ilc2,ibc2)

        !Get the actual location of the cell w.r.t the edge cc
        cc_cell(1) = plane_torus_closest_coordinates(cc_edge%x,cc_cell(1)%x, &
                                                     ptr_patch%geometry_info)
        cc_cell(2) = plane_torus_closest_coordinates(cc_edge%x,cc_cell(2)%x, &
                                                     ptr_patch%geometry_info)

        !Get the distance vector between cc_edge and cc_cell
        dist_edge_cell(1)%x(:) = cc_cell(1)%x(:) - cc_edge%x(:)
        dist_edge_cell(2)%x(:) = cc_cell(2)%x(:) - cc_edge%x(:)


        !
        ! 2. rotate these vectors into a new local cartesian system. In this rotated
        !    system the coordinate axes point into the local normal and tangential
        !    direction at each edge. All these vectors are 2D so no point using the
        !    z coordinate
        !

        ! centers
        !
        xyloc_plane_nt_n1(1) =  SUM( dist_edge_cell(1)%x(1:2) *  &
             ptr_patch%edges%primal_cart_normal(je,jb)%x(1:2) )

        xyloc_plane_nt_n1(2) =  SUM( dist_edge_cell(1)%x(1:2) *  &
             ptr_patch%edges%dual_cart_normal(je,jb)%x(1:2) )

        xyloc_plane_nt_n2(1) =  SUM( dist_edge_cell(2)%x(1:2) *  &
             ptr_patch%edges%primal_cart_normal(je,jb)%x(1:2) )

        xyloc_plane_nt_n2(2) =  SUM( dist_edge_cell(2)%x(1:2) *  &
             ptr_patch%edges%dual_cart_normal(je,jb)%x(1:2) )

        ! vertices
        !
        DO nv = 1,4
           xyloc_plane_nt_v(nv,1) =  SUM( dist_edge_vert(nv)%x(1:2) *  &
                ptr_patch%edges%primal_cart_normal(je,jb)%x(1:2) )

           xyloc_plane_nt_v(nv,2) =  SUM( dist_edge_vert(nv)%x(1:2) *  &
                ptr_patch%edges%dual_cart_normal(je,jb)%x(1:2) )
        END DO

        !
        !AD: From this point on the sphere and the torus version are same
        !    except for the multiplication with grid_sphere_radius

        ! 3. Calculate position of vertices in a translated coordinate system.
        !    This is done twice. The origin is located once at the circumcenter
        !    of the neighboring cell 1 and once cell 2. The distance vectors point
        !    from the cell center to the vertices.
        xyloc_trans1_v(1:4,1) = xyloc_plane_nt_v(1:4,1) - xyloc_plane_nt_n1(1)
        xyloc_trans1_v(1:4,2) = xyloc_plane_nt_v(1:4,2) - xyloc_plane_nt_n1(2)

        xyloc_trans2_v(1:4,1) = xyloc_plane_nt_v(1:4,1) - xyloc_plane_nt_n2(1)
        xyloc_trans2_v(1:4,2) = xyloc_plane_nt_v(1:4,2) - xyloc_plane_nt_n2(2)



        ! 4. Rotate points into coordinate system pointing into local north
        !    and local east direction. Store in edge-based data structure. This
        !    is done twice (for both neighboring cells).
        ! e_n= pn_cell1(1)*e_\lambda + pn_cell1(2)*e_\phi
        ! e_t= dn_cell1(1)*e_\lambda + dn_cell1(2)*e_\phi
        pn_cell1(1) = ptr_patch%edges%primal_normal_cell(je,jb,1)%v1
        pn_cell1(2) = ptr_patch%edges%primal_normal_cell(je,jb,1)%v2
        dn_cell1(1) = ptr_patch%edges%dual_normal_cell(je,jb,1)%v1
        dn_cell1(2) = ptr_patch%edges%dual_normal_cell(je,jb,1)%v2

        pn_cell2(1) = ptr_patch%edges%primal_normal_cell(je,jb,2)%v1
        pn_cell2(2) = ptr_patch%edges%primal_normal_cell(je,jb,2)%v2
        dn_cell2(1) = ptr_patch%edges%dual_normal_cell(je,jb,2)%v1
        dn_cell2(2) = ptr_patch%edges%dual_normal_cell(je,jb,2)%v2

        ! components in longitudinal direction (cell 1)
        ptr_int%pos_on_tplane_c_edge(je,jb,1,1:3)%lon =  &
          &  ( xyloc_trans1_v(1:3,1) * pn_cell1(1)       &
          & +  xyloc_trans1_v(1:3,2) * dn_cell1(1) )

        ! components in latitudinal direction (cell 1)
        ptr_int%pos_on_tplane_c_edge(je,jb,1,1:3)%lat =   &
          &  ( xyloc_trans1_v(1:3,1) * pn_cell1(2)        &
          & +  xyloc_trans1_v(1:3,2) * dn_cell1(2) )

        ! components in longitudinal direction (cell 2)
        ptr_int%pos_on_tplane_c_edge(je,jb,2,1:3)%lon =  &
          &  ( xyloc_trans2_v((/1,2,4/),1) * pn_cell2(1) &
          & +  xyloc_trans2_v((/1,2,4/),2) * dn_cell2(1) )

        ! components in latitudinal direction (cell 2)
        ptr_int%pos_on_tplane_c_edge(je,jb,2,1:3)%lat =  &
          &  ( xyloc_trans2_v((/1,2,4/),1) * pn_cell2(2) &
          & +  xyloc_trans2_v((/1,2,4/),2) * dn_cell2(2) )

      ENDDO  ! je
    ENDDO  ! jb
!$OMP END DO
!$OMP END PARALLEL



    !
    ! Projection of butterfly cell centers
    !

    i_rcstartlev = 3

    ! start and end block
    i_startblk = ptr_patch%edges%start_blk(i_rcstartlev,1)
    i_endblk   = ptr_patch%nblks_e

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,je,i_startidx,i_endidx,ilc1,ilc2,ibc1,ibc2,cc_cell,&
!$OMP            ilc_bf1,ibc_bf1,ilc_bf2,ibc_bf2,cc_bf1,    &
!$OMP            cc_bf2,dist_bf1_cell1,dist_bf2_cell2)
    DO jb = i_startblk, i_endblk

      CALL get_indices_e(ptr_patch, jb, i_startblk, i_endblk, &
        & i_startidx, i_endidx, i_rcstartlev)

      DO je = i_startidx, i_endidx

        IF(.NOT.ptr_patch%edges%decomp_info%owner_mask(je,jb)) CYCLE


        ! get line and block indices of edge-neighbours (cells)
        ilc1 = ptr_patch%edges%cell_idx(je,jb,1)
        ibc1 = ptr_patch%edges%cell_blk(je,jb,1)
        ilc2 = ptr_patch%edges%cell_idx(je,jb,2)
        ibc2 = ptr_patch%edges%cell_blk(je,jb,2)

        ! get CC of edge-neighbor 1
        cc_cell(1) = ptr_patch%cells%cartesian_center(ilc1,ibc1)

        ! get CC of edge-neighbor 2
        cc_cell(2) = ptr_patch%cells%cartesian_center(ilc2,ibc2)


        ! 5a. project cell centers adjacent to the cells sharing edge je.
        !    (from butterfly_idx)


        ! get line and block indices of the neighbors of edge-neighbor 1
        ilc_bf1(1:2) = ptr_patch%edges%butterfly_idx(je,jb,1,1:2)
        ibc_bf1(1:2) = ptr_patch%edges%butterfly_blk(je,jb,1,1:2)

        ! get line and block indices of the neighbors of edge-neighbor 2
        ilc_bf2(1:2) = ptr_patch%edges%butterfly_idx(je,jb,2,1:2)
        ibc_bf2(1:2) = ptr_patch%edges%butterfly_blk(je,jb,2,1:2)


        ! get CC cell centers (neighbor 1)
        cc_bf1(1) = ptr_patch%cells%cartesian_center(ilc_bf1(1),ibc_bf1(1))
        cc_bf1(2) = ptr_patch%cells%cartesian_center(ilc_bf1(2),ibc_bf1(2))

        ! get CC cell centers (neighbor 1)
        cc_bf2(1) = ptr_patch%cells%cartesian_center(ilc_bf2(1),ibc_bf2(1))
        cc_bf2(2) = ptr_patch%cells%cartesian_center(ilc_bf2(2),ibc_bf2(2))


        ! project cell centers into cell-based local system (edge-neighbor 1)
        cc_bf1(1) = plane_torus_closest_coordinates(cc_cell(1)%x, cc_bf1(1)%x, &
                                                   ptr_patch%geometry_info)
        dist_bf1_cell1(1)%x(:) = cc_bf1(1)%x(:) - cc_cell(1)%x(:)

        cc_bf1(2) = plane_torus_closest_coordinates(cc_cell(1)%x, cc_bf1(2)%x, &
                                                   ptr_patch%geometry_info)
        dist_bf1_cell1(2)%x(:) = cc_bf1(2)%x(:) - cc_cell(1)%x(:)


        ! project cell centers into cell-based local system (edge-neighbor 2)

        cc_bf2(1) = plane_torus_closest_coordinates(cc_cell(2)%x, cc_bf2(1)%x, &
                                                   ptr_patch%geometry_info)
        dist_bf2_cell2(1)%x(:) = cc_bf2(1)%x(:) - cc_cell(2)%x(:)

        cc_bf2(2) = plane_torus_closest_coordinates(cc_cell(2)%x, cc_bf2(2)%x, &
                                                   ptr_patch%geometry_info)
        dist_bf2_cell2(2)%x(:) = cc_bf2(2)%x(:) - cc_cell(2)%x(:)


        ! components in longitudinal direction (or X for torus) (cell 1)
        ptr_int%pos_on_tplane_c_edge(je,jb,1,4:5)%lon = dist_bf1_cell1(1:2)%x(1)

        ! components in latitudinal direction (or Y for torus) (cell 1)
        ptr_int%pos_on_tplane_c_edge(je,jb,1,4:5)%lat = dist_bf1_cell1(1:2)%x(2)


        ! components in longitudinal (or X for torus) direction (cell 2)
        ptr_int%pos_on_tplane_c_edge(je,jb,2,4:5)%lon =  dist_bf2_cell2(1:2)%x(1)

        ! components in latitudinal (or Y for torus) direction (cell 2)
        ptr_int%pos_on_tplane_c_edge(je,jb,2,4:5)%lat =  dist_bf2_cell2(1:2)%x(2)

      ENDDO  ! je
    ENDDO  ! jb
!$OMP END DO
!$OMP END PARALLEL

    DO nv=1,5
      DO nc=1,2
        CALL sync_patch_array(sync_e, ptr_patch, ptr_int%pos_on_tplane_c_edge(:,:,nc,nv)%lon)
        CALL sync_patch_array(sync_e, ptr_patch, ptr_int%pos_on_tplane_c_edge(:,:,nc,nv)%lat)
      ENDDO
    ENDDO


  END SUBROUTINE init_tplane_c_torus

  !----------------------------------------------------------------------------


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

    TYPE(t_patch), INTENT(inout) :: ptr_patch  !< patch

    TYPE(t_int_state), INTENT(inout) :: ptr_int  !< interpolation state

    REAL(wp) ::  alpha_l(3),        & !< area coordinates for quadrature up to
      & alpha_q(3,3), alpha_c(3,4)   !< fourth order
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
    i_endblk   = ptr_patch%nblks_c

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,nv,nq,i_startidx,i_endidx,ilv,ibv,z_vert_cc,z_quad_cc, &
!$OMP           z_quad_gg) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
        & i_startidx, i_endidx, i_rcstartlev)

      DO jc = i_startidx, i_endidx

        IF(.NOT.ptr_patch%cells%decomp_info%owner_mask(jc,jb)) CYCLE

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
          & + alpha_l(2)*z_vert_cc(2)%x(1)  &
          & + alpha_l(3)*z_vert_cc(3)%x(1)

        z_quad_cc%x(2)= alpha_l(1)*z_vert_cc(1)%x(2)  &
          & + alpha_l(2)*z_vert_cc(2)%x(2)  &
          & + alpha_l(3)*z_vert_cc(3)%x(2)

        z_quad_cc%x(3)= alpha_l(1)*z_vert_cc(1)%x(3)  &
          & + alpha_l(2)*z_vert_cc(2)%x(3)  &
          & + alpha_l(3)*z_vert_cc(3)%x(3)


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
            & + alpha_q(2,nq)*z_vert_cc(2)%x(1)  &
            & + alpha_q(3,nq)*z_vert_cc(3)%x(1)

          z_quad_cc%x(2)= alpha_q(1,nq)*z_vert_cc(1)%x(2)  &
            & + alpha_q(2,nq)*z_vert_cc(2)%x(2)  &
            & + alpha_q(3,nq)*z_vert_cc(3)%x(2)

          z_quad_cc%x(3)= alpha_q(1,nq)*z_vert_cc(1)%x(3)  &
            & + alpha_q(2,nq)*z_vert_cc(2)%x(3)  &
            & + alpha_q(3,nq)*z_vert_cc(3)%x(3)


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
            & + alpha_c(2,nq)*z_vert_cc(2)%x(1)  &
            & + alpha_c(3,nq)*z_vert_cc(3)%x(1)

          z_quad_cc%x(2)= alpha_c(1,nq)*z_vert_cc(1)%x(2)  &
            & + alpha_c(2,nq)*z_vert_cc(2)%x(2)  &
            & + alpha_c(3,nq)*z_vert_cc(3)%x(2)

          z_quad_cc%x(3)= alpha_c(1,nq)*z_vert_cc(1)%x(3)  &
            & + alpha_c(2,nq)*z_vert_cc(2)%x(3)  &
            & + alpha_c(3,nq)*z_vert_cc(3)%x(3)


          ! Transform back to geographical coordinates
          z_quad_gg = cc2gc(z_quad_cc)

          ! store
          ptr_int%gquad%qpts_tri_c(jc,jb,nq)%lat = z_quad_gg%lat
          ptr_int%gquad%qpts_tri_c(jc,jb,nq)%lon = z_quad_gg%lon
        ENDDO

      END DO
    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    CALL sync_patch_array(sync_c,ptr_patch,ptr_int%gquad%qpts_tri_l(:,:)%lat)
    CALL sync_patch_array(sync_c,ptr_patch,ptr_int%gquad%qpts_tri_l(:,:)%lon)
    DO nq=1,3
      CALL sync_patch_array(sync_c,ptr_patch,ptr_int%gquad%qpts_tri_q(:,:,nq)%lat)
      CALL sync_patch_array(sync_c,ptr_patch,ptr_int%gquad%qpts_tri_q(:,:,nq)%lon)
    ENDDO
    DO nq=1,4
      CALL sync_patch_array(sync_c,ptr_patch,ptr_int%gquad%qpts_tri_c(:,:,nq)%lat)
      CALL sync_patch_array(sync_c,ptr_patch,ptr_int%gquad%qpts_tri_c(:,:,nq)%lon)
    ENDDO

  END SUBROUTINE tri_quadrature_pts
  !-------------------------------------------------------------------------


END MODULE mo_intp_coeffs
