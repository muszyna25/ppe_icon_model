
! #ifdef __xlC__
! @PROCESS HOT
! #endif
#ifdef __PGI
!pgi$g opt=1
#endif
!>
!! Contains the interpolation routines needed for grid refinement.
!!
!! These had originally been included in mo_grf_interpolation but then were
!! packed into a separate module to clean up the code
!!
!! @par Revision History
!! Created by Guenther Zaengl, DWD (2009-02-09)
!! Modification by Guenther Zaengl, DWD (2009-06-22)
!! - preparation for generalized grid refinement (affects all subroutines)
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

MODULE mo_grf_intp_coeffs
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
USE mo_exception,           ONLY: finish
USE mo_impl_constants,      ONLY: SUCCESS, min_rlcell_int, min_rledge_int
USE mo_model_domain,        ONLY: t_patch, t_grid_edges, t_grid_cells, t_grid_vertices, &
 &                                p_patch_local_parent

USE mo_grid_config,         ONLY: n_dom, n_dom_start, grid_sphere_radius

USE mo_math_utilities,      ONLY: gc2cc, gvec2cvec, arc_length, &
                                  t_cartesian_coordinates, arc_length_v
USE mo_math_utility_solvers, ONLY: solve_chol_v, choldec_v

USE mo_impl_constants_grf,  ONLY: grf_bdyintp_start_c, grf_bdyintp_start_e,  &
                                  grf_fbk_start_c

USE mo_parallel_config,     ONLY: nproma
USE mo_loopindices,         ONLY: get_indices_c, get_indices_e
USE mo_mpi,                 ONLY: my_process_is_mpi_parallel
USE mo_communication,       ONLY: exchange_data
USE mo_sync,                ONLY: global_sum_array

USE mo_grf_intp_data_strc
USE mo_gridref_config

IMPLICIT NONE

PRIVATE

TYPE(t_patch), POINTER :: p_patch(:)

PUBLIC :: grf_intp_coeffs_setpatch
PUBLIC :: compute_pc2cc_distances, compute_pe2ce_distances, gridref_info,         &
        & init_fbk_wgt, grf_index, rbf_compute_coeff_grf_e, idw_compute_coeff_grf_e
CONTAINS

#include "intp_functions.inc"


SUBROUTINE grf_intp_coeffs_setpatch(in_patch)

  TYPE(t_patch), TARGET, INTENT(IN)  :: in_patch(n_dom_start:)

  p_patch => in_patch

END SUBROUTINE grf_intp_coeffs_setpatch


!-------------------------------------------------------------------------
!
!>
!! This routine computes the distances between the parent cell and its children,
!! which are needed for interpolation to the child cells using gradients at
!! the parent cell center
!!
!! @par Revision History
!! Developed  by Guenther. Zaengl, DWD, 2009-12-16
!!
SUBROUTINE compute_pc2cc_distances(p_grf)
!
! patch and gridref state

TYPE(t_gridref_state), TARGET, INTENT(INOUT) :: p_grf  (n_dom_start:)

! local variables

TYPE(t_patch),      POINTER :: p_pp => NULL()
TYPE(t_patch),      POINTER :: p_pc => NULL()
TYPE(t_grid_cells), POINTER :: p_cp => NULL()
TYPE(t_grid_cells), POINTER :: p_cc => NULL()

TYPE(t_gridref_single_state), POINTER :: p_grfs


INTEGER :: jb, jc, jg, jcd, jgc, i_startblk, i_endblk, &
           i_startidx, i_endidx, ici1, icb1, ici2, icb2, ici3, icb3, ici4, icb4

REAL(wp), DIMENSION (3) :: z_nx1, z_nx2
REAL(wp) :: z_lon, z_lat, z_norm
TYPE(t_cartesian_coordinates) :: cc_center, cc_ch1, cc_ch2, cc_ch3, cc_ch4, &
                               cc_dis1, cc_dis2, cc_dis3, cc_dis4


LEV_LOOP: DO jg = n_dom_start, n_dom-1

 CD_LOOP: DO jcd = 1, p_patch(jg)%n_childdom

  jgc    =  p_patch(jg)%child_id(jcd)
  p_pc   => p_patch(jgc)

  p_pp   => p_patch_local_parent(jgc)
  p_grfs => p_grf_state_local_parent(jgc)%p_dom(jcd)

  p_cp   => p_pp%cells
  p_cc   => p_pc%cells

  ! Start and end blocks for which coefficients ar needed
  i_startblk = p_cp%start_blk(grf_bdyintp_start_c,jcd)
  i_endblk   = p_cp%end_blk(min_rlcell_int,jcd)


  DO jb = i_startblk, i_endblk

    CALL get_indices_c(p_pp, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, grf_bdyintp_start_c, min_rlcell_int)

    DO jc = i_startidx, i_endidx

      cc_center = gc2cc(p_cp%center(jc,jb))
      z_lon     = p_cp%center(jc,jb)%lon
      z_lat     = p_cp%center(jc,jb)%lat

      CALL gvec2cvec(1._wp,0._wp,z_lon,z_lat,z_nx1(1),z_nx1(2),z_nx1(3))
      z_norm = SQRT( DOT_PRODUCT(z_nx1(1:3),z_nx1(1:3)) )
      z_nx1(1:3)  = 1._wp/z_norm * z_nx1(1:3)

      CALL gvec2cvec(0._wp,1._wp,z_lon,z_lat,z_nx2(1),z_nx2(2),z_nx2(3))
      z_norm = SQRT( DOT_PRODUCT(z_nx2(1:3),z_nx2(1:3)) )
      z_nx2(1:3)  = 1._wp/z_norm * z_nx2(1:3)

      ici1 = p_cp%child_idx(jc,jb,1)
      icb1 = p_cp%child_blk(jc,jb,1)
      ici2 = p_cp%child_idx(jc,jb,2)
      icb2 = p_cp%child_blk(jc,jb,2)
      ici3 = p_cp%child_idx(jc,jb,3)
      icb3 = p_cp%child_blk(jc,jb,3)
      ici4 = p_cp%child_idx(jc,jb,4)
      icb4 = p_cp%child_blk(jc,jb,4)

      cc_ch1 = gc2cc(p_cc%center(ici1,icb1))
      cc_ch2 = gc2cc(p_cc%center(ici2,icb2))
      cc_ch3 = gc2cc(p_cc%center(ici3,icb3))
      cc_ch4 = gc2cc(p_cc%center(ici4,icb4))

      cc_dis1%x(1:3) = cc_ch1%x(1:3) - cc_center%x(1:3)
      z_norm = SQRT( DOT_PRODUCT(cc_dis1%x(1:3),cc_dis1%x(1:3)) )
      cc_dis1%x(1:3) = cc_dis1%x(1:3)/z_norm

      cc_dis2%x(1:3) = cc_ch2%x(1:3) - cc_center%x(1:3)
      z_norm = SQRT( DOT_PRODUCT(cc_dis2%x(1:3),cc_dis2%x(1:3)) )
      cc_dis2%x(1:3) = cc_dis2%x(1:3)/z_norm

      cc_dis3%x(1:3) = cc_ch3%x(1:3) - cc_center%x(1:3)
      z_norm = SQRT( DOT_PRODUCT(cc_dis3%x(1:3),cc_dis3%x(1:3)) )
      cc_dis3%x(1:3) = cc_dis3%x(1:3)/z_norm

      cc_dis4%x(1:3) = cc_ch4%x(1:3) - cc_center%x(1:3)
      z_norm = SQRT( DOT_PRODUCT(cc_dis4%x(1:3),cc_dis4%x(1:3)) )
      cc_dis4%x(1:3) = cc_dis4%x(1:3)/z_norm

      p_grfs%grf_dist_pc2cc(jc,1,1,jb) = grid_sphere_radius* &
        arc_length(cc_center,cc_ch1)*DOT_PRODUCT(z_nx1(1:3),cc_dis1%x(1:3))
      p_grfs%grf_dist_pc2cc(jc,2,1,jb) = grid_sphere_radius* &
        arc_length(cc_center,cc_ch2)*DOT_PRODUCT(z_nx1(1:3),cc_dis2%x(1:3))
      p_grfs%grf_dist_pc2cc(jc,3,1,jb) = grid_sphere_radius* &
        arc_length(cc_center,cc_ch3)*DOT_PRODUCT(z_nx1(1:3),cc_dis3%x(1:3))
      p_grfs%grf_dist_pc2cc(jc,4,1,jb) = grid_sphere_radius* &
        arc_length(cc_center,cc_ch4)*DOT_PRODUCT(z_nx1(1:3),cc_dis4%x(1:3))

      p_grfs%grf_dist_pc2cc(jc,1,2,jb) = grid_sphere_radius* &
        arc_length(cc_center,cc_ch1)*DOT_PRODUCT(z_nx2(1:3),cc_dis1%x(1:3))
      p_grfs%grf_dist_pc2cc(jc,2,2,jb) = grid_sphere_radius* &
        arc_length(cc_center,cc_ch2)*DOT_PRODUCT(z_nx2(1:3),cc_dis2%x(1:3))
      p_grfs%grf_dist_pc2cc(jc,3,2,jb) = grid_sphere_radius* &
        arc_length(cc_center,cc_ch3)*DOT_PRODUCT(z_nx2(1:3),cc_dis3%x(1:3))
      p_grfs%grf_dist_pc2cc(jc,4,2,jb) = grid_sphere_radius* &
        arc_length(cc_center,cc_ch4)*DOT_PRODUCT(z_nx2(1:3),cc_dis4%x(1:3))

    ENDDO
  ENDDO

 ENDDO CD_LOOP
ENDDO LEV_LOOP

END SUBROUTINE compute_pc2cc_distances


!-------------------------------------------------------------------------
!
!>
!! This routine computes the distances between the parent edge and its children,
!! which are needed for interpolation to the child edges using gradients at
!! the parent edge
!!
!! @par Revision History
!! Developed  by Guenther. Zaengl, DWD, 2010-03-11
!!
SUBROUTINE compute_pe2ce_distances( p_grf)
!
! patch and gridref state

TYPE(t_gridref_state), TARGET, INTENT(INOUT) :: p_grf  (n_dom_start:)

! local variables

TYPE(t_patch),      POINTER :: p_pp => NULL()
TYPE(t_patch),      POINTER :: p_pc => NULL()
TYPE(t_grid_edges), POINTER :: p_ep => NULL()
TYPE(t_grid_edges), POINTER :: p_ec => NULL()

TYPE(t_gridref_single_state), POINTER :: p_grfs


INTEGER :: jb, je, jg, jcd, jgc, i_startblk, i_endblk, &
           i_startidx, i_endidx, ici1, icb1, ici2, icb2

REAL(wp), DIMENSION (3) :: z_nx
REAL(wp) :: z_norm
TYPE(t_cartesian_coordinates) :: cc_center, cc_ch1, cc_ch2, cc_dis1, cc_dis2


LEV_LOOP: DO jg = n_dom_start, n_dom-1

 CD_LOOP: DO jcd = 1, p_patch(jg)%n_childdom

  jgc    =  p_patch(jg)%child_id(jcd)
  p_pc   => p_patch(jgc)

  p_pp   => p_patch_local_parent(jgc)
  p_grfs => p_grf_state_local_parent(jgc)%p_dom(jcd)

  p_ep   => p_pp%edges
  p_ec   => p_pc%edges

  ! Start and end blocks for which coefficients are needed
  i_startblk = p_ep%start_blk(grf_bdyintp_start_e,jcd)
  i_endblk   = p_ep%end_blk(min_rledge_int,jcd)


  DO jb = i_startblk, i_endblk

    CALL get_indices_e(p_pp, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, grf_bdyintp_start_e, min_rledge_int)

    DO je = i_startidx, i_endidx

      cc_center = gc2cc(p_ep%center(je,jb))

      z_nx(:) = p_ep%dual_cart_normal(je,jb)%x(:)

      ! child edges 1 and 2
      ici1 = p_ep%child_idx(je,jb,1)
      icb1 = p_ep%child_blk(je,jb,1)
      ici2 = p_ep%child_idx(je,jb,2)
      icb2 = p_ep%child_blk(je,jb,2)

      cc_ch1 = gc2cc(p_ec%center(ici1,icb1))
      cc_ch2 = gc2cc(p_ec%center(ici2,icb2))

      cc_dis1%x(1:3) = cc_ch1%x(1:3) - cc_center%x(1:3)
      z_norm = SQRT( DOT_PRODUCT(cc_dis1%x(1:3),cc_dis1%x(1:3)) )
      cc_dis1%x(1:3) = cc_dis1%x(1:3)/z_norm

      cc_dis2%x(1:3) = cc_ch2%x(1:3) - cc_center%x(1:3)
      z_norm = SQRT( DOT_PRODUCT(cc_dis2%x(1:3),cc_dis2%x(1:3)) )
      cc_dis2%x(1:3) = cc_dis2%x(1:3)/z_norm

      ! The distance in tangential direction is needed normalized with the primal edge length
      p_grfs%grf_dist_pe2ce(je,1,jb) = grid_sphere_radius*p_ep%tangent_orientation(je,jb)* &
        arc_length(cc_center,cc_ch1)*DOT_PRODUCT(z_nx(1:3),cc_dis1%x(1:3))/&
        (p_ep%primal_edge_length(je,jb)*ABS(DOT_PRODUCT(z_nx(1:3),cc_dis1%x(1:3))))
      p_grfs%grf_dist_pe2ce(je,2,jb) = grid_sphere_radius*p_ep%tangent_orientation(je,jb)* &
        arc_length(cc_center,cc_ch2)*DOT_PRODUCT(z_nx(1:3),cc_dis2%x(1:3))/&
        (p_ep%primal_edge_length(je,jb)*ABS(DOT_PRODUCT(z_nx(1:3),cc_dis2%x(1:3))))

    ENDDO
  ENDDO

 ENDDO CD_LOOP
ENDDO LEV_LOOP

END SUBROUTINE compute_pe2ce_distances

!-------------------------------------------------------------------------
!
!
!
!>
!! This routine computes fbk_dom_area
!!
!! @par Revision History
!! Developed  by Guenther. Zaengl, DWD, 2009-03-19
!! Reduced to the computation of fbk_dom_area, Rainer Johanni, 2011-12-20
!!
SUBROUTINE gridref_info( p_grf)


TYPE(t_gridref_state), TARGET, INTENT(INOUT) :: p_grf  (n_dom_start:)

! local variables

TYPE(t_grid_cells), POINTER :: p_gcp => NULL()
TYPE(t_patch),      POINTER :: p_pp => NULL()

INTEGER :: jb, jc, jg, ji, i_startblk, i_endblk,         &
           i_startidx, i_endidx, i_nchdom, icid


REAL(wp), ALLOCATABLE :: z_area(:,:)

!-----------------------------------------------------------------------
!

DO jg = n_dom_start, n_dom-1

  p_gcp  => p_patch(jg)%cells
  p_pp   => p_patch(jg)

  i_nchdom = p_pp%n_childdom
  IF (i_nchdom == 0) CYCLE

! c) Compute area of feedback domains
  DO ji = 1, i_nchdom

    icid = p_pp%child_id(ji)
    p_grf(jg)%fbk_dom_area(ji) = 0._wp

    ALLOCATE(z_area(nproma,p_pp%nblks_c))
    z_area = 0._wp

    i_startblk = p_gcp%start_blk(1,1)
    i_endblk   = p_gcp%end_blk(min_rlcell_int,i_nchdom)

    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_pp, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, 1, min_rlcell_int)

      DO jc = i_startidx, i_endidx

        IF (p_gcp%refin_ctrl(jc,jb) <= grf_fbk_start_c .AND. p_gcp%child_id(jc,jb) == icid .AND. &
           p_gcp%decomp_info%owner_mask(jc,jb)) z_area(jc,jb) = p_gcp%area(jc,jb)

      ENDDO
    ENDDO

    p_grf(jg)%fbk_dom_area(ji) = global_sum_array(z_area)
    DEALLOCATE(z_area)

  ENDDO

ENDDO


END SUBROUTINE gridref_info

!-------------------------------------------------------------------------
!
!
!
!>
!! This routine computes the feedback weighting coefficients needed for.
!!
!! This routine computes the feedback weighting coefficients needed for
!! cell-based variables. Arithmetic or area-weighted averaging turned out to be
!! inadequate because it aliases horizontal gradients on the fine mesh into a
!! checkerboard noise pattern between upward- and downward-directed triangles
!! on the coarse mesh
!!
!! @par Revision History
!! Developed  by Guenther. Zaengl, DWD, 2008-09-12
!! Moved from hierarchy_management to grf_interpolation (2009-03-12)
!! Add feedback weight computation for edge-based variables (2009-03-19)
!!
SUBROUTINE init_fbk_wgt( p_grf)


TYPE(t_gridref_state), TARGET, INTENT(INOUT) :: p_grf  (n_dom_start:)

! local variables

TYPE(t_grid_cells), POINTER :: p_gcp => NULL()
TYPE(t_grid_cells), POINTER :: p_gcc => NULL()
TYPE(t_grid_edges), POINTER :: p_gep => NULL()
TYPE(t_grid_edges), POINTER :: p_gec => NULL()
TYPE(t_patch),      POINTER :: p_pp => NULL()
TYPE(t_patch),      POINTER :: p_pc => NULL()

TYPE(t_gridref_state), POINTER :: p_grfp

INTEGER  :: jb, jc, jg, je, j, ki(2), kb(2), i_startblk, i_endblk, &
            i_startidx, i_endidx, j1, j2, js1, js2, i_chidx
REAL(wp) :: sum1,wgt(4),x(4),y(4)

INTEGER  ::  ici1, icb1, ici2, icb2, ici3, icb3, ici4, icb4, ierror

INTEGER, ALLOCATABLE :: ierrcount(:)

REAL(wp), DIMENSION (3) :: z_nx1, z_nx2
REAL(wp) :: z_lon, z_lat, z_norm
TYPE(t_cartesian_coordinates) :: cc_center, cc_ch1, cc_ch2, cc_ch3, cc_ch4, &
                               cc_dis1, cc_dis2, cc_dis3, cc_dis4

!-----------------------------------------------------------------------
!
! Part 1: Feedback weights for cell-based variables

DO jg = n_dom_start+1, n_dom

  p_pc => p_patch(jg)

  p_pp   => p_patch_local_parent(jg)
  p_grfp => p_grf_state_local_parent(jg)

  i_chidx = p_pc%parent_child_index

  p_gcp  => p_pp%cells
  p_gcc  => p_pc%cells
  p_gep  => p_pp%edges
  p_gec  => p_pc%edges

  i_chidx = p_pc%parent_child_index

  ALLOCATE (ierrcount(p_pp%nblks_c))
  ierrcount(:) = 0

  ! If the nested domain is barely larger than the boundary interpolation zone,
  ! the setting of the start and end indices may fail in the presence of
  ! multiple nests per nesting level. As such small nests do not make sense anyway,
  ! we just stop here in such pathological cases.
  ierror = 0
  i_startblk = p_gcp%start_blk(grf_bdyintp_start_c,i_chidx)
  i_endblk   = p_gcp%end_blk(min_rlcell_int,i_chidx)
  DO jb = i_startblk, i_endblk

    CALL get_indices_c(p_pp, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, grf_bdyintp_start_c, min_rlcell_int)

    DO jc = i_startidx, i_endidx
      IF (p_gcp%child_id(jc,jb) /= jg) ierror = ierror + 1
    ENDDO
  ENDDO
  IF (ierror > 0) THEN
    write(0,*) 'Error in domain ID: ',jg
    CALL finish ('mo_grf_interpolation:init_fbk_wgt',  &
      &          'size of nested domain is too small')
  ENDIF

!$OMP PARALLEL PRIVATE(wgt,i_startblk,i_endblk)
  i_startblk = p_gcp%start_blk(grf_bdyintp_start_c,i_chidx)
  i_endblk   = p_gcp%end_blk(min_rlcell_int,i_chidx)

  ! 1a) Area-weighted feedback weights

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,j,ici1,icb1) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_pp, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, grf_bdyintp_start_c, min_rlcell_int)

      DO jc = i_startidx, i_endidx

        DO j = 1, 4
          ici1 = p_gcp%child_idx(jc,jb,j)
          icb1 = p_gcp%child_blk(jc,jb,j)

          p_grfp%fbk_wgt_aw(jc,jb,1:4) = p_gcc%area(ici1,icb1)/p_gcp%area(jc,jb)
        ENDDO

      ENDDO
    ENDDO
!$OMP END DO

  ! 1b) Bilinear feedback weights

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,z_lon,z_lat,z_norm,z_nx1,z_nx2,x,y, &
!$OMP            ici1,icb1,ici2,icb2,ici3,icb3,ici4,icb4,cc_center,            &
!$OMP            cc_ch1,cc_ch2,cc_ch3,cc_ch4,cc_dis1,cc_dis2,cc_dis3,          &
!$OMP  cc_dis4) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_pp, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, grf_bdyintp_start_c, min_rlcell_int)

      DO jc = i_startidx, i_endidx

        cc_center = gc2cc(p_gcp%center(jc,jb))
        z_lon     = p_gcp%center(jc,jb)%lon
        z_lat     = p_gcp%center(jc,jb)%lat

        CALL gvec2cvec(1._wp,0._wp,z_lon,z_lat,z_nx1(1),z_nx1(2),z_nx1(3))
        z_norm = SQRT( DOT_PRODUCT(z_nx1(1:3),z_nx1(1:3)) )
        z_nx1(1:3)  = 1._wp/z_norm * z_nx1(1:3)

        CALL gvec2cvec(0._wp,1._wp,z_lon,z_lat,z_nx2(1),z_nx2(2),z_nx2(3))
        z_norm = SQRT( DOT_PRODUCT(z_nx2(1:3),z_nx2(1:3)) )
        z_nx2(1:3)  = 1._wp/z_norm * z_nx2(1:3)

        ici1 = p_gcp%child_idx(jc,jb,1)
        icb1 = p_gcp%child_blk(jc,jb,1)
        ici2 = p_gcp%child_idx(jc,jb,2)
        icb2 = p_gcp%child_blk(jc,jb,2)
        ici3 = p_gcp%child_idx(jc,jb,3)
        icb3 = p_gcp%child_blk(jc,jb,3)
        ici4 = p_gcp%child_idx(jc,jb,4)
        icb4 = p_gcp%child_blk(jc,jb,4)

        cc_ch1 = gc2cc(p_gcc%center(ici1,icb1))
        cc_ch2 = gc2cc(p_gcc%center(ici2,icb2))
        cc_ch3 = gc2cc(p_gcc%center(ici3,icb3))
        cc_ch4 = gc2cc(p_gcc%center(ici4,icb4))

        cc_dis1%x(1:3) = cc_ch1%x(1:3) - cc_center%x(1:3)
        z_norm = SQRT( DOT_PRODUCT(cc_dis1%x(1:3),cc_dis1%x(1:3)) )
        cc_dis1%x(1:3) = cc_dis1%x(1:3)/z_norm

        cc_dis2%x(1:3) = cc_ch2%x(1:3) - cc_center%x(1:3)
        z_norm = SQRT( DOT_PRODUCT(cc_dis2%x(1:3),cc_dis2%x(1:3)) )
        cc_dis2%x(1:3) = cc_dis2%x(1:3)/z_norm

        cc_dis3%x(1:3) = cc_ch3%x(1:3) - cc_center%x(1:3)
        z_norm = SQRT( DOT_PRODUCT(cc_dis3%x(1:3),cc_dis3%x(1:3)) )
        cc_dis3%x(1:3) = cc_dis3%x(1:3)/z_norm

        cc_dis4%x(1:3) = cc_ch4%x(1:3) - cc_center%x(1:3)
        z_norm = SQRT( DOT_PRODUCT(cc_dis4%x(1:3),cc_dis4%x(1:3)) )
        cc_dis4%x(1:3) = cc_dis4%x(1:3)/z_norm

        x(1) = arc_length(cc_center,cc_ch1)*DOT_PRODUCT(z_nx1(1:3),cc_dis1%x(1:3))
        x(2) = arc_length(cc_center,cc_ch2)*DOT_PRODUCT(z_nx1(1:3),cc_dis2%x(1:3))
        x(3) = arc_length(cc_center,cc_ch3)*DOT_PRODUCT(z_nx1(1:3),cc_dis3%x(1:3))
        x(4) = arc_length(cc_center,cc_ch4)*DOT_PRODUCT(z_nx1(1:3),cc_dis4%x(1:3))

        y(1) = arc_length(cc_center,cc_ch1)*DOT_PRODUCT(z_nx2(1:3),cc_dis1%x(1:3))
        y(2) = arc_length(cc_center,cc_ch2)*DOT_PRODUCT(z_nx2(1:3),cc_dis2%x(1:3))
        y(3) = arc_length(cc_center,cc_ch3)*DOT_PRODUCT(z_nx2(1:3),cc_dis3%x(1:3))
        y(4) = arc_length(cc_center,cc_ch4)*DOT_PRODUCT(z_nx2(1:3),cc_dis4%x(1:3))

        ! Use area fraction for weight of central child cell
        wgt(3) = p_gcc%area(ici3,icb3) /                   &
          (p_gcc%area(ici1,icb1) + p_gcc%area(ici2,icb2) + &
           p_gcc%area(ici3,icb3) + p_gcc%area(ici4,icb4))

        ! The weighting factors are based on the requirement that sum(w(i)*x(i)) = 0
        ! and sum(w(i)*y(i)) = 0, which ensures that linear horizontal gradients
        ! are not aliased into a checkerboard pattern between upward- and downward
        ! directed cells. The third condition is sum(w(i)) = 1., and one coefficient
        ! can be freely chosen (see above). Analytical elimination yields...
        !
        IF (ABS(x(2)-x(1)) > 1.e-11_wp) THEN
          wgt(4) = 1.0_wp/( (y(4)-y(1)) - (y(2)-y(1))*(x(4)-x(1))/(x(2)-x(1)) )* &
                   (-y(1) + x(1)*(y(2)-y(1))/(x(2)-x(1)) - &
                   wgt(3) *( (y(3)-y(1)) - (y(2)-y(1))*(x(3)-x(1))/(x(2)-x(1))) )
          wgt(2) = (-x(1) - wgt(3)*(x(3)-x(1)) - wgt(4)*(x(4)-x(1)))/(x(2)-x(1))
          wgt(1) = 1.0_wp - SUM(wgt(2:4))
        ELSE IF (ABS(x(4)-x(1)) > 1.e-11_wp) THEN
          wgt(2) = 1.0_wp/( (y(2)-y(1)) - (y(4)-y(1))*(x(2)-x(1))/(x(4)-x(1)) )* &
                   (-y(1) + x(1)*(y(4)-y(1))/(x(4)-x(1)) - &
                   wgt(3) *( (y(3)-y(1)) - (y(4)-y(1))*(x(3)-x(1))/(x(4)-x(1))) )
          wgt(4) = (-x(1) - wgt(3)*(x(3)-x(1)) - wgt(2)*(x(2)-x(1)))/(x(4)-x(1))
          wgt(1) = 1.0_wp - SUM(wgt(2:4))
        ELSE
          ierrcount(jb) = ierrcount(jb) + 1
        ENDIF

        IF (MINVAL(wgt(1:4)) < 0.0_wp) ierrcount(jb) = ierrcount(jb) + 1

        ! Save the weighting factors in fbk_wgt_bln
        p_grfp%fbk_wgt_bln(jc,jb,1:4) = wgt(1:4)
      ENDDO
    ENDDO
!$OMP END DO


! Part 2: Feedback weights for edge-based variables

  IF (grf_velfbk == 1) THEN ! Inverse-distance weighted interpolation to parent edge

    i_startblk = p_gep%start_blk(grf_bdyintp_start_e,i_chidx)
    i_endblk   = p_gep%end_blk(min_rledge_int,i_chidx)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,je,j,sum1,ki,kb) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_e(p_pp, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, grf_bdyintp_start_e, min_rledge_int)

      DO je = i_startidx, i_endidx
        sum1 = 0._wp
        DO j = 1, 2
          ki(j) = p_gep%child_idx(je,jb,j)
          kb(j) = p_gep%child_blk(je,jb,j)
          sum1 = sum1 + p_gec%primal_edge_length(ki(j),kb(j))
        ENDDO
        DO j = 1, 2
          ki(j) = p_gep%child_idx(je,jb,j)
          kb(j) = p_gep%child_blk(je,jb,j)
          p_grfp%fbk_wgt_e(je,jb,j) = 1._wp - p_gec%primal_edge_length(ki(j),kb(j))/sum1
        ENDDO
      ENDDO
    ENDDO
!$OMP END DO

  ELSE IF (grf_velfbk == 2) THEN ! Second-order interpolation of normal velocities
                                 ! using RBF reconstruction to child vertices

    i_startblk = p_gep%start_blk(grf_bdyintp_start_e,i_chidx)
    i_endblk   = p_gep%end_blk(min_rledge_int,i_chidx)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,je,j,j1,j2,ki,kb,js1,js2) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_e(p_pp, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, grf_bdyintp_start_e, min_rledge_int)

      js1 = 0
      js2 = 0

      DO je = i_startidx, i_endidx

        DO j = 1, 2
          ki(j) = p_gep%child_idx(je,jb,j)
          kb(j) = p_gep%child_blk(je,jb,j)
        ENDDO
        DO j1 = 1, 2
          DO j2 = 1, 2
            IF(p_gec%vertex_idx(ki(1),kb(1),j1)==p_gec%vertex_idx(ki(2),kb(2),j2)) THEN
              js1 = j1
              js2 = j2
            ENDIF
          ENDDO
        ENDDO
        p_grfp%fbk_wgt_e(je,jb,1:2) = 32._wp/90._wp
        IF (js1 == 1) THEN
          p_grfp%fbk_wgt_e(je,jb,3) = 6._wp/90._wp
          p_grfp%fbk_wgt_e(je,jb,4) = 7._wp/90._wp
        ELSE
          p_grfp%fbk_wgt_e(je,jb,3) = 7._wp/90._wp
          p_grfp%fbk_wgt_e(je,jb,4) = 6._wp/90._wp
        ENDIF
        IF (js2 == 1) THEN
          p_grfp%fbk_wgt_e(je,jb,5) = 6._wp/90._wp
          p_grfp%fbk_wgt_e(je,jb,6) = 7._wp/90._wp
        ELSE
          p_grfp%fbk_wgt_e(je,jb,5) = 7._wp/90._wp
          p_grfp%fbk_wgt_e(je,jb,6) = 6._wp/90._wp
        ENDIF
      ENDDO
    ENDDO
!$OMP END DO NOWAIT

  ENDIF
!$OMP END PARALLEL

  IF (ANY(ierrcount(:) > 0)) THEN
    CALL finish ('init_fbk_wgt',  &
      &          'negative feedback coefficients occurred - change grid optimization')
  ENDIF

  IF(my_process_is_mpi_parallel()) THEN
    DO j = 1, 4
      CALL exchange_data(p_pp%comm_pat_c, p_grfp%fbk_wgt_aw(:,:,j))
    ENDDO
    DO j = 1, 4
      CALL exchange_data(p_pp%comm_pat_c, p_grfp%fbk_wgt_bln(:,:,j))
    ENDDO
    DO j = 1, 6
      CALL exchange_data(p_pp%comm_pat_e, p_grfp%fbk_wgt_e(:,:,j))
    ENDDO
  ENDIF

  DEALLOCATE (ierrcount)
ENDDO

END SUBROUTINE init_fbk_wgt

!-------------------------------------------------------------------------
!
!
!>
!! This routine initializes the indices used to define the stencils.
!!
!! This routine initializes the indices used to define the stencils
!! for coarse-to-fine-mesh interpolation. The stencils can be combined
!! with RBF and IDW interpolation.
!! Because both scalars and wind vectors need to be interpolated to the
!! nest boundaries, different groups of indices are computed.
!!
!! @par Revision History
!! Developed and tested  by G. Zaengl (June 2008)
!! Rewritten for vectorization by G. Zaengl (May 2010)
!!
SUBROUTINE grf_index( ptr_grf_state)
!

TYPE(t_gridref_state), TARGET, INTENT(inout) ::  &
  &  ptr_grf_state(n_dom_start:) ! State type for grid refinement coefficients

TYPE(t_patch),      POINTER :: p_pp => NULL()
TYPE(t_patch),      POINTER :: p_pc => NULL()

TYPE(t_grid_edges),POINTER  :: ptr_ep, ptr_ec ! pointer to parent and child edges
TYPE(t_grid_cells),POINTER  :: ptr_cp, ptr_cc ! pointer to parent and child cells
TYPE(t_grid_vertices),POINTER  :: ptr_vp      ! pointer to parent vertices

TYPE(t_gridref_single_state),POINTER :: ptr_grf ! pointer to gridref_state for a
                                              ! single child domain


! Loop indices and other index variables
INTEGER  :: jg, jb, je, jcd, jgc,                      &
            iipv1, ibpv1, iipv2, ibpv2,                &
            iice, ibce, iicc, ibcc, iipc, ibpc, iie, ibe

INTEGER, DIMENSION(nproma) :: ii2, ii3, ii4, ii5, ierror

INTEGER  :: i_startblk       ! start block
INTEGER  :: i_endblk         ! end block
INTEGER  :: i_startidx       ! start index
INTEGER  :: i_endidx         ! end index

! Auxiliaries for computing the cartesian orientation of the edge points
REAL(wp) :: z_nx1(3),z_nx2(3), thresh_dotpr

TYPE(t_patch), POINTER :: ptr_patch(:)! indices of all data points
    CHARACTER(LEN=*), PARAMETER :: method_name = 'mo_grf_intp_coeffs:grf_index'
!--------------------------------------------------------------------

ptr_patch => p_patch

LEV_LOOP: DO jg = n_dom_start, n_dom-1

 CD_LOOP: DO jcd = 1, ptr_patch(jg)%n_childdom

    jgc  = ptr_patch(jg)%child_id(jcd)
    p_pc => p_patch(jgc)

    p_pp    => p_patch_local_parent(jgc)
    ptr_grf => p_grf_state_local_parent(jgc)%p_dom(jcd)

    ptr_ep   => p_pp%edges
    ptr_ec   => p_pc%edges
    ptr_cp   => p_pp%cells
    ptr_cc   => p_pc%cells
    ptr_vp   => p_pp%verts

    ! This is to check the correct orientation of stencil points 4 and 5
    ! for child edges 3 and 4
    thresh_dotpr = 0.85_wp

    ! Start and end blocks for which vector interpolation is needed
    i_startblk = p_pp%edges%start_blk(grf_bdyintp_start_e,jcd)
    i_endblk   = p_pp%edges%end_blk(min_rledge_int,jcd)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,je,i_startidx,i_endidx,iipv1,ibpv1,iipv2,ibpv2,ierror, &
!$OMP    iice,ibce,iicc,ibcc,iipc,ibpc,ii2,ii3,ii4,ii5,iie,ibe,            &
!$OMP    z_nx1,z_nx2) ICON_OMP_DEFAULT_SCHEDULE

    DO jb = i_startblk, i_endblk

      CALL get_indices_e(p_pp, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, grf_bdyintp_start_e, min_rledge_int)

      ! Part 1: child edges aligned with the parent edge
      !
      ! The search algorithm is as follows:
      ! First, we determine the parent vertex point adjacent to the target edge;
      ! then, we take the six (five for pentagon points) parent edges
      ! bordering to this vertex

      DO je = i_startidx, i_endidx

        ! vertices of parent edge
        iipv1 = ptr_ep%vertex_idx(je,jb,1)
        ibpv1 = ptr_ep%vertex_blk(je,jb,1)
        iipv2 = ptr_ep%vertex_idx(je,jb,2)
        ibpv2 = ptr_ep%vertex_blk(je,jb,2)

        IF (ptr_grf%grf_dist_pe2ce(je,1,jb) > 0._wp) THEN

          ! in this case, parent vertex 2(1) belongs to child edge 1(2)

          ! Note: at pentagon points, edge index/block 6 is filled
          ! with the values of edge 5 in model_domain_import. Thus, no case
          ! discrimination is needed here.

          ptr_grf%grf_vec_stencil_1a(je,jb) = ptr_vp%num_edges(iipv2,ibpv2)
!CDIR EXPAND=6
          ptr_grf%grf_vec_ind_1a(je,1:6,jb) = ptr_vp%edge_idx(iipv2,ibpv2,1:6)
!CDIR EXPAND=6
          ptr_grf%grf_vec_blk_1a(je,1:6,jb) = ptr_vp%edge_blk(iipv2,ibpv2,1:6)

          ptr_grf%grf_vec_stencil_1b(je,jb) = ptr_vp%num_edges(iipv1,ibpv1)
!CDIR EXPAND=6
          ptr_grf%grf_vec_ind_1b(je,1:6,jb) = ptr_vp%edge_idx(iipv1,ibpv1,1:6)
!CDIR EXPAND=6
          ptr_grf%grf_vec_blk_1b(je,1:6,jb) = ptr_vp%edge_blk(iipv1,ibpv1,1:6)

        ELSE ! ptr_grf%grf_dist_pe2ce(je,1,jb) < 0

          ! in this case, parent vertex 1(2) belongs to child edge 1(2)

          ptr_grf%grf_vec_stencil_1a(je,jb) = ptr_vp%num_edges(iipv1,ibpv1)

!CDIR EXPAND=6
          ptr_grf%grf_vec_ind_1a(je,1:6,jb) = ptr_vp%edge_idx(iipv1,ibpv1,1:6)
!CDIR EXPAND=6
          ptr_grf%grf_vec_blk_1a(je,1:6,jb) = ptr_vp%edge_blk(iipv1,ibpv1,1:6)

          ptr_grf%grf_vec_stencil_1b(je,jb) = ptr_vp%num_edges(iipv2,ibpv2)

!CDIR EXPAND=6
          ptr_grf%grf_vec_ind_1b(je,1:6,jb) = ptr_vp%edge_idx(iipv2,ibpv2,1:6)
!CDIR EXPAND=6
          ptr_grf%grf_vec_blk_1b(je,1:6,jb) = ptr_vp%edge_blk(iipv2,ibpv2,1:6)

        ENDIF

      ENDDO  ! computation for child edges 1 and 2 finished


      ! Part 2: interior child edges (i.e. child edges 3 and 4)
      !
      ! The interpolation stencil consists of 5 points, namely the 3 edges
      ! of the parent cell in which the child cell is located, plus the 2
      ! edges of the neighbor cells of the parent cell that have (approximately)
      ! the same orientation as the child edge

      ierror = 0

      DO je = i_startidx, i_endidx

        ptr_grf%grf_vec_stencil_2a(je,jb) = 5

        ! The first point of the stencil is the parent edge
        ptr_grf%grf_vec_ind_2a(je,1,jb) = je
        ptr_grf%grf_vec_blk_2a(je,1,jb) = jb

        ! Determine index of the parent cell in which child edge 3 is located

        iice = ABS(ptr_ep%child_idx(je,jb,3))
        ibce = ptr_ep%child_blk(je,jb,3)
        ! In the inner domain, it does not matter if we take cell 1 or 2
        ! At the boundary row, we have to take the one which is included in the child boundary
        iicc = ptr_ec%cell_idx(iice,ibce,1)
        ibcc = ptr_ec%cell_blk(iice,ibce,1)
        IF ( .NOT. (ptr_cc%edge_idx(iicc,ibcc,1) == iice .AND. ptr_cc%edge_blk(iicc,ibcc,1) == ibce .OR. &
                    ptr_cc%edge_idx(iicc,ibcc,2) == iice .AND. ptr_cc%edge_blk(iicc,ibcc,2) == ibce .OR. &
                    ptr_cc%edge_idx(iicc,ibcc,3) == iice .AND. ptr_cc%edge_blk(iicc,ibcc,3) == ibce )) THEN
          iicc = ptr_ec%cell_idx(iice,ibce,2)
          ibcc = ptr_ec%cell_blk(iice,ibce,2)
        ENDIF
        iipc = ptr_cc%parent_loc_idx(iicc,ibcc)
        ibpc = ptr_cc%parent_loc_blk(iicc,ibcc)

        ! Determine stencil points 2 and 3 (remaining edges of the parent cell)
        IF (ptr_cp%edge_idx(iipc,ibpc,1) == je .AND. ptr_cp%edge_blk(iipc,ibpc,1) == jb ) THEN
          ii2(je) = 2
          ii3(je) = 3
        ELSE IF (ptr_cp%edge_idx(iipc,ibpc,2) == je .AND. ptr_cp%edge_blk(iipc,ibpc,2) == jb ) THEN
          ii2(je) = 3
          ii3(je) = 1
        ELSE IF (ptr_cp%edge_idx(iipc,ibpc,3) == je .AND. ptr_cp%edge_blk(iipc,ibpc,3) == jb ) THEN
          ii2(je) = 1
          ii3(je) = 2
        ELSE
          CALL finish(method_name, "Undefined stencil points 2 and 3 (remaining edges of the parent cell)")
        ENDIF

        ptr_grf%grf_vec_ind_2a(je,2,jb) = ptr_cp%edge_idx(iipc,ibpc,ii2(je))
        ptr_grf%grf_vec_blk_2a(je,2,jb) = ptr_cp%edge_blk(iipc,ibpc,ii2(je))
        ptr_grf%grf_vec_ind_2a(je,3,jb) = ptr_cp%edge_idx(iipc,ibpc,ii3(je))
        ptr_grf%grf_vec_blk_2a(je,3,jb) = ptr_cp%edge_blk(iipc,ibpc,ii3(je))

        ! Finally, determine stencil points 4 and 5
        iie = ptr_grf%grf_vec_ind_2a(je,2,jb)
        ibe = ptr_grf%grf_vec_blk_2a(je,2,jb)
        IF (ptr_ep%quad_idx(iie,ibe,1) == je .AND. ptr_ep%quad_blk(iie,ibe,1) == jb) THEN
          ii4(je) = 3
        ELSE IF (ptr_ep%quad_idx(iie,ibe,2) == je .AND. ptr_ep%quad_blk(iie,ibe,2) == jb) THEN
          ii4(je) = 4
        ELSE IF (ptr_ep%quad_idx(iie,ibe,3) == je .AND. ptr_ep%quad_blk(iie,ibe,3) == jb) THEN
          ii4(je) = 1
        ELSE IF (ptr_ep%quad_idx(iie,ibe,4) == je .AND. ptr_ep%quad_blk(iie,ibe,4) == jb) THEN
          ii4(je) = 2
        ELSE
          CALL finish(method_name, "Undefined stencil point 4 (remaining edges of the parent cell)")
        ENDIF

        ptr_grf%grf_vec_ind_2a(je,4,jb) = ptr_ep%quad_idx(iie,ibe,ii4(je))
        ptr_grf%grf_vec_blk_2a(je,4,jb) = ptr_ep%quad_blk(iie,ibe,ii4(je))

        iie = ptr_grf%grf_vec_ind_2a(je,3,jb)
        ibe = ptr_grf%grf_vec_blk_2a(je,3,jb)
        IF (ptr_ep%quad_idx(iie,ibe,1) == je .AND. ptr_ep%quad_blk(iie,ibe,1) == jb) THEN
          ii5(je) = 3
        ELSE IF (ptr_ep%quad_idx(iie,ibe,2) == je .AND. ptr_ep%quad_blk(iie,ibe,2) == jb) THEN
          ii5(je) = 4
        ELSE IF (ptr_ep%quad_idx(iie,ibe,3) == je .AND. ptr_ep%quad_blk(iie,ibe,3) == jb) THEN
          ii5(je) = 1
        ELSE IF (ptr_ep%quad_idx(iie,ibe,4) == je .AND. ptr_ep%quad_blk(iie,ibe,4) == jb) THEN
          ii5(je) = 2
        ELSE
          CALL finish(method_name, "Undefined stencil point 5 (remaining edges of the parent cell)")
        ENDIF

        ptr_grf%grf_vec_ind_2a(je,5,jb) = ptr_ep%quad_idx(iie,ibe,ii5(je))
        ptr_grf%grf_vec_blk_2a(je,5,jb) = ptr_ep%quad_blk(iie,ibe,ii5(je))

        ! ** Check orientation of stencil points 4 and 5 **

        ! Compute cartesian orientation of child edge
        z_nx1(:) = ptr_ec%primal_cart_normal(iice,ibce)%x(:)

        ! Compute cartesian orientation of stencil point 4
        iie   = ptr_grf%grf_vec_ind_2a(je,4,jb)
        ibe   = ptr_grf%grf_vec_blk_2a(je,4,jb)

        z_nx2(:) = ptr_ep%primal_cart_normal(iie,ibe)%x(:)

        IF (ABS(DOT_PRODUCT(z_nx1,z_nx2)) < thresh_dotpr) ierror(je) = ierror(je) + 1

        ! Compute cartesian orientation of stencil point 5
        iie   = ptr_grf%grf_vec_ind_2a(je,5,jb)
        ibe   = ptr_grf%grf_vec_blk_2a(je,5,jb)

        z_nx2(:) = ptr_ep%primal_cart_normal(iie,ibe)%x(:)

        IF (ABS(DOT_PRODUCT(z_nx1,z_nx2)) < thresh_dotpr) ierror(je) = ierror(je) + 1


        ! Now process child edge 4
        IF (ptr_ep%refin_ctrl(je,jb) == -1) THEN ! in this case, child edge 4 does not exist
          ptr_grf%grf_vec_stencil_2b(je,jb) = 0
          CYCLE
        ENDIF
        ptr_grf%grf_vec_stencil_2b(je,jb) = 5

        ! The first point of the stencil is the parent edge
        ptr_grf%grf_vec_ind_2b(je,1,jb) = je
        ptr_grf%grf_vec_blk_2b(je,1,jb) = jb

        ! Determine index of the parent cell in which child edge 4 is located

        iice = ABS(ptr_ep%child_idx(je,jb,4))
        ibce = ptr_ep%child_blk(je,jb,4)
        ! In the inner domain, it does not matter if we take cell 1 or 2
        ! At the boundary row, we have to take the one which is included in the child boundary
        iicc = ptr_ec%cell_idx(iice,ibce,1)
        ibcc = ptr_ec%cell_blk(iice,ibce,1)
        IF ( .NOT. (ptr_cc%edge_idx(iicc,ibcc,1) == iice .AND. ptr_cc%edge_blk(iicc,ibcc,1) == ibce .OR. &
                    ptr_cc%edge_idx(iicc,ibcc,2) == iice .AND. ptr_cc%edge_blk(iicc,ibcc,2) == ibce .OR. &
                    ptr_cc%edge_idx(iicc,ibcc,3) == iice .AND. ptr_cc%edge_blk(iicc,ibcc,3) == ibce )) THEN
          iicc = ptr_ec%cell_idx(iice,ibce,2)
          ibcc = ptr_ec%cell_blk(iice,ibce,2)
        ENDIF
        iipc = ptr_cc%parent_loc_idx(iicc,ibcc)
        ibpc = ptr_cc%parent_loc_blk(iicc,ibcc)

        ! Determine stencil points 2 and 3 (remaining edges of the parent cell)
        IF (ptr_cp%edge_idx(iipc,ibpc,1) == je .AND. ptr_cp%edge_blk(iipc,ibpc,1) == jb ) THEN
          ii2(je) = 2
          ii3(je) = 3
        ELSE IF (ptr_cp%edge_idx(iipc,ibpc,2) == je .AND. ptr_cp%edge_blk(iipc,ibpc,2) == jb ) THEN
          ii2(je) = 3
          ii3(je) = 1
        !ELSE
        ELSE IF (ptr_cp%edge_idx(iipc,ibpc,3) == je .AND. ptr_cp%edge_blk(iipc,ibpc,3) == jb ) THEN
          ii2(je) = 1
          ii3(je) = 2
        ELSE
          CALL finish(method_name, "Undefined stencil point 2 and 3 (2nd remaining edges of the parent cell)")
        ENDIF

        ptr_grf%grf_vec_ind_2b(je,2,jb) = ptr_cp%edge_idx(iipc,ibpc,ii2(je))
        ptr_grf%grf_vec_blk_2b(je,2,jb) = ptr_cp%edge_blk(iipc,ibpc,ii2(je))
        ptr_grf%grf_vec_ind_2b(je,3,jb) = ptr_cp%edge_idx(iipc,ibpc,ii3(je))
        ptr_grf%grf_vec_blk_2b(je,3,jb) = ptr_cp%edge_blk(iipc,ibpc,ii3(je))

        ! Finally, determine stencil points 4 and 5
        iie = ptr_grf%grf_vec_ind_2b(je,2,jb)
        ibe = ptr_grf%grf_vec_blk_2b(je,2,jb)
        IF (ptr_ep%quad_idx(iie,ibe,1) == je .AND. ptr_ep%quad_blk(iie,ibe,1) == jb) THEN
          ii4(je) = 3
        ELSE IF (ptr_ep%quad_idx(iie,ibe,2) == je .AND. ptr_ep%quad_blk(iie,ibe,2) == jb) THEN
          ii4(je) = 4
        ELSE IF (ptr_ep%quad_idx(iie,ibe,3) == je .AND. ptr_ep%quad_blk(iie,ibe,3) == jb) THEN
          ii4(je) = 1
        ELSE IF (ptr_ep%quad_idx(iie,ibe,4) == je .AND. ptr_ep%quad_blk(iie,ibe,4) == jb) THEN
          ii4(je) = 2
        ELSE
          CALL finish(method_name, "Undefined stencil point 4 (2nd remaining edges of the parent cell)")
        ENDIF

        ptr_grf%grf_vec_ind_2b(je,4,jb) = ptr_ep%quad_idx(iie,ibe,ii4(je))
        ptr_grf%grf_vec_blk_2b(je,4,jb) = ptr_ep%quad_blk(iie,ibe,ii4(je))

        iie = ptr_grf%grf_vec_ind_2b(je,3,jb)
        ibe = ptr_grf%grf_vec_blk_2b(je,3,jb)
        IF (ptr_ep%quad_idx(iie,ibe,1) == je .AND. ptr_ep%quad_blk(iie,ibe,1) == jb) THEN
          ii5(je) = 3
        ELSE IF (ptr_ep%quad_idx(iie,ibe,2) == je .AND. ptr_ep%quad_blk(iie,ibe,2) == jb) THEN
          ii5(je) = 4
        ELSE IF (ptr_ep%quad_idx(iie,ibe,3) == je .AND. ptr_ep%quad_blk(iie,ibe,3) == jb) THEN
          ii5(je) = 1
        ELSE IF (ptr_ep%quad_idx(iie,ibe,4) == je .AND. ptr_ep%quad_blk(iie,ibe,4) == jb) THEN
          ii5(je) = 2
        ELSE
          CALL finish(method_name, "Undefined stencil point 5 (2nd remaining edges of the parent cell)")
        ENDIF

        ptr_grf%grf_vec_ind_2b(je,5,jb) = ptr_ep%quad_idx(iie,ibe,ii5(je))
        ptr_grf%grf_vec_blk_2b(je,5,jb) = ptr_ep%quad_blk(iie,ibe,ii5(je))

        ! ** Check orientation of stencil points 4 and 5 **

        ! Compute cartesian orientation of child edge
        z_nx1(:) = ptr_ec%primal_cart_normal(iice,ibce)%x(:)

        ! Compute cartesian orientation of stencil point 4
        iie   = ptr_grf%grf_vec_ind_2b(je,4,jb)
        ibe   = ptr_grf%grf_vec_blk_2b(je,4,jb)

        ! write(0,*) iie, ibe
        z_nx2(:) = ptr_ep%primal_cart_normal(iie,ibe)%x(:)

        IF (ABS(DOT_PRODUCT(z_nx1,z_nx2)) < thresh_dotpr) ierror(je) = ierror(je) + 1

        ! Compute cartesian orientation of stencil point 5
        iie   = ptr_grf%grf_vec_ind_2b(je,5,jb)
        ibe   = ptr_grf%grf_vec_blk_2b(je,5,jb)

        z_nx2(:) = ptr_ep%primal_cart_normal(iie,ibe)%x(:)

        IF (ABS(DOT_PRODUCT(z_nx1,z_nx2)) < thresh_dotpr) ierror(je) = ierror(je) + 1

      ENDDO

      IF (MAXVAL(ierror) > 0) THEN
        write(0,*) 'Number of errors: ',SUM(ierror(1:nproma))
        CALL finish ('mo_grf_interpolation:grf_index',  &
          &          'orientation of edge points incorrect')
      ENDIF

    ENDDO ! block loop
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  ENDDO CD_LOOP
ENDDO LEV_LOOP

END SUBROUTINE grf_index

!-------------------------------------------------------------------------
!
!
!>
!! This routine computes the coefficients needed for vector RBF interpolation.
!!
!! This routine computes the coefficients needed for vector RBF interpolation to
!! child edges, which is required along the lateral boundaries
!! of nested model domains.
!! This computation involves the inversion of the interpolation matrix, which
!! is performed by a Cholesky decomposition.
!! The Cholesky decomposition is currently implemented by a home made routine
!! which can be substituted by a call to a numerical library, if available.
!!
!! @par Revision History
!! Developed and tested by Guenther Zaengl (May 2008)
!!
SUBROUTINE rbf_compute_coeff_grf_e (ptr_grf_state)
!

TYPE(t_gridref_state), TARGET, INTENT(inout) :: ptr_grf_state(n_dom_start:)

TYPE(t_patch),      POINTER :: p_pp => NULL()
TYPE(t_patch),      POINTER :: p_pc => NULL()

TYPE(t_grid_edges),POINTER  :: ptr_ep, ptr_ec ! pointer to parent and child edges

TYPE(t_gridref_single_state),POINTER :: ptr_grf ! pointer to gridref_state for a
                                              ! single child domain

TYPE(t_patch), POINTER      :: ptr_patch(:)
                                                             ! coordinates of ...
TYPE(t_cartesian_coordinates) :: cc_e1, cc_e2, cc_childedge    ! edge midpoints
REAL(wp)           :: cc_cer(nproma,3), cc_e2r(3)

REAL(wp) :: z_nx1(nproma,3),z_nx2(nproma,3) ! 3d  normal velocity
                                            ! vectors at edge midpoints

REAL(wp)           :: z_dist                ! distance between data points

REAL(wp)           :: z_nxprod              ! scalar product of normal
                                            ! velocity vectors

REAL(wp), ALLOCATABLE :: z_rbfmat(:,:,:)    ! RBF interpolation matrix

REAL(wp), ALLOCATABLE :: z_diag(:,:)        ! diagonal of cholesky
                                            ! decomposition matrix

REAL(wp), ALLOCATABLE :: z_rbfval(:,:)      ! RBF function value


! Index variables
INTEGER :: jg, jb, je, jgc, jcd, jce, &
           iie1, ibe1, iie2, ibe2, je1, je2, iiec, ibec

! Other local variables
INTEGER :: ist, ishift, max_points, istencil(nproma)


INTEGER :: i_startblk                ! start block
INTEGER :: i_endblk                  ! end index
INTEGER :: i_startidx                ! start index
INTEGER :: i_endidx                  ! end index

REAL(wp) ::  checksum_vt             ! to check if sum of interpolation coefficients is correct

! Pointers to index and coefficient fields
INTEGER, DIMENSION(:,:,:), POINTER :: iidx, iblk
INTEGER, DIMENSION(:,:),   POINTER :: ptr_stc
REAL(wp), DIMENSION(:,:,:),POINTER :: ptr_coeff

!--------------------------------------------------------------------

ptr_patch => p_patch

LEV_LOOP: DO jg = n_dom_start, n_dom-1

 CD_LOOP: DO jcd = 1, ptr_patch(jg)%n_childdom

  jgc  = ptr_patch(jg)%child_id(jcd)
  p_pc => p_patch(jgc)

  p_pp    => p_patch_local_parent(jgc)
  ptr_grf => p_grf_state_local_parent(jgc)%p_dom(jcd)

  ptr_ep   => p_pp%edges
  ptr_ec   => p_pc%edges

  CE_LOOP: DO jce = 1, 4  ! loop over child edges

    IF (jce == 1) THEN
      ptr_coeff  => ptr_grf%grf_vec_coeff_1a
      ptr_stc    => ptr_grf%grf_vec_stencil_1a
      iidx       => ptr_grf%grf_vec_ind_1a
      iblk       => ptr_grf%grf_vec_blk_1a
      ishift     =  0
      max_points =  6
    ELSE IF (jce == 2) THEN
      ptr_coeff  => ptr_grf%grf_vec_coeff_1b
      ptr_stc    => ptr_grf%grf_vec_stencil_1b
      iidx       => ptr_grf%grf_vec_ind_1b
      iblk       => ptr_grf%grf_vec_blk_1b
      ishift     =  0
      max_points =  6
    ELSE IF (jce == 3) THEN
      ptr_coeff  => ptr_grf%grf_vec_coeff_2a
      ptr_stc    => ptr_grf%grf_vec_stencil_2a
      iidx       => ptr_grf%grf_vec_ind_2a
      iblk       => ptr_grf%grf_vec_blk_2a
      ishift     =  0
      max_points =  5
    ELSE IF (jce == 4) THEN
      ptr_coeff  => ptr_grf%grf_vec_coeff_2b
      ptr_stc    => ptr_grf%grf_vec_stencil_2b
      iidx       => ptr_grf%grf_vec_ind_2b
      iblk       => ptr_grf%grf_vec_blk_2b
      ishift     =  1
      max_points =  5
    ENDIF

!$OMP PARALLEL PRIVATE(z_rbfmat,z_diag,z_rbfval,ist,i_startblk,i_endblk)

    ALLOCATE( z_rbfmat(nproma,max_points,max_points),  &
              z_diag(nproma,max_points),               &
              z_rbfval(nproma,max_points), STAT=ist )
    IF (ist /= SUCCESS) THEN
          CALL finish ('mo_grf_interpolation:rbf_compute_coeff_grf',      &
        &             'allocation for working arrays failed')
    ENDIF

    ! Start and end blocks for which vector interpolation is needed
    i_startblk = p_pp%edges%start_blk(grf_bdyintp_start_e-ishift,jcd)
    i_endblk   = p_pp%edges%end_blk(min_rledge_int,jcd)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,je1,je2,je,istencil,iie1,ibe1,iie2,ibe2,       &
!$OMP            cc_e1,cc_e2,z_nx1,z_nx2,z_dist,z_nxprod, &
!$OMP            iiec,ibec,cc_childedge,cc_cer,cc_e2r,checksum_vt) ICON_OMP_DEFAULT_SCHEDULE
    DO jb =  i_startblk, i_endblk

      CALL get_indices_e(p_pp, jb, i_startblk, i_endblk, i_startidx, i_endidx,&
                         grf_bdyintp_start_e-ishift, min_rledge_int)

      !
      ! for each edge, build the vector RBF interpolation matrix
      !
      DO je1 = 1, max_points

        DO je2 = 1, je1

          DO je = i_startidx, i_endidx

            istencil(je) = ptr_stc(je,jb)
            !
            IF ( (je1 > istencil(je)) .OR. (je2 > istencil(je)) ) CYCLE
            !
            ! line and block indices for each edge je1 and je2 of RBF stencil
            !
            iie1 = iidx(je,je1,jb)
            ibe1 = iblk(je,je1,jb)
            iie2 = iidx(je,je2,jb)
            ibe2 = iblk(je,je2,jb)
            !
            ! get Cartesian coordinates
            !
            cc_e1 = gc2cc(ptr_ep%center(iie1,ibe1))
            cc_e2 = gc2cc(ptr_ep%center(iie2,ibe2))

            ! normal component and coordinates of edge midpoints for each edge of RBF stencil

            z_nx1(je,:) = ptr_ep%primal_cart_normal(iie1,ibe1)%x(:)

            z_nx2(je,:) = ptr_ep%primal_cart_normal(iie2,ibe2)%x(:)

            z_nxprod = DOT_PRODUCT(z_nx1(je,:),z_nx2(je,:))
            z_dist   = arc_length(cc_e1,cc_e2)

            ! set up interpolation matrix
            !
            IF      (rbf_vec_kern_grf_e == 1) THEN
              z_rbfmat(je,je1,je2) = z_nxprod * gaussi(z_dist,rbf_scale_grf_e(MAX(jg,1)))
            ELSE IF (rbf_vec_kern_grf_e == 2) THEN
              z_rbfmat(je,je1,je2) = z_nxprod * inv_multiq2(z_dist,rbf_scale_grf_e(MAX(jg,1)))
            ELSE IF (rbf_vec_kern_grf_e == 3) THEN
              z_rbfmat(je,je1,je2) = z_nxprod * inv_multiq(z_dist,rbf_scale_grf_e(MAX(jg,1)))
            ENDIF

            IF (je1 > je2) z_rbfmat(je,je2,je1) = z_rbfmat(je,je1,je2)

          END DO
        END DO
      END DO


      ! apply Cholesky decomposition to matrix
      !
!CDIR NOIEXPAND
#ifdef __SX__
      CALL choldec_v(i_startidx,i_endidx,istencil,max_points,z_rbfmat,z_diag)
#else
      CALL choldec_v(i_startidx,i_endidx,istencil,           z_rbfmat,z_diag)
#endif

      DO je = i_startidx, i_endidx

        !
        ! Solve immediately for coefficients
        !
        ! convert coordinates of edge midpoint to cartesian vector
        !

        iiec = ABS(ptr_ep%child_idx(je,jb,jce)) ! current child edge
        ibec = ptr_ep%child_blk(je,jb,jce)
        cc_childedge = gc2cc(ptr_ec%center(iiec,ibec))
        cc_cer(je,1:3) = cc_childedge%x(1:3)

        z_nx1(je,:) = ptr_ec%primal_cart_normal(iiec,ibec)%x(:)

      END DO

      !
      ! set up right hand side for interpolation system
      !
      DO je2 = 1, max_points

        DO je = i_startidx, i_endidx

          IF (je2 > istencil(je)) CYCLE
          !
          ! get indices and coordinates of edge midpoints and compute distance
          ! to the edge where the vector is reconstructed
          !
          iie2 = iidx(je,je2,jb)
          ibe2 = iblk(je,je2,jb)

          cc_e2 = gc2cc(ptr_ep%center(iie2,ibe2))
          cc_e2r(1:3) = cc_e2%x(1:3)
          z_dist = arc_length_v(cc_cer(je,:), cc_e2r)
          !
          ! get Cartesian orientation vector
          z_nx2(je,:) = ptr_ep%primal_cart_normal(iie2,ibe2)%x(:)

          z_nxprod = DOT_PRODUCT(z_nx1(je,:),z_nx2(je,:))

          IF      (rbf_vec_kern_grf_e == 1) THEN
            z_rbfval(je,je2) = z_nxprod * gaussi(z_dist,rbf_scale_grf_e(MAX(jg,1)))
          ELSE IF (rbf_vec_kern_grf_e == 2) THEN
            z_rbfval(je,je2) = z_nxprod * inv_multiq2(z_dist,rbf_scale_grf_e(MAX(jg,1)))
          ELSE IF (rbf_vec_kern_grf_e == 3) THEN
            z_rbfval(je,je2) = z_nxprod * inv_multiq(z_dist,rbf_scale_grf_e(MAX(jg,1)))
          ENDIF

        END DO
      END DO

      !
      ! compute vector coefficients
      !
!CDIR NOIEXPAND
#ifdef __SX__
      CALL solve_chol_v(i_startidx, i_endidx, istencil, max_points, z_rbfmat, &
                        z_diag, z_rbfval, ptr_coeff(:,:,jb))
#else
      CALL solve_chol_v(i_startidx, i_endidx, istencil,             z_rbfmat, &
                        z_diag, z_rbfval, ptr_coeff(:,:,jb))
#endif

      DO je = i_startidx, i_endidx

        ! Ensure that sum of interpolation coefficients is correct

        checksum_vt = 0._wp

        DO je1 = 1, istencil(je)
          iie1   = iidx(je,je1,jb)
          ibe1   = iblk(je,je1,jb)

          z_nx2(je,:) = ptr_ep%primal_cart_normal(iie1,ibe1)%x(:)

          checksum_vt = checksum_vt+ptr_coeff(je1,je,jb)*DOT_PRODUCT(z_nx1(je,:),z_nx2(je,:))
        ENDDO

        DO je1 = 1, istencil(je)
          ptr_coeff(je1,je,jb) = ptr_coeff(je1,je,jb) / checksum_vt
        ENDDO

      END DO
    END DO
!$OMP END DO NOWAIT

    DEALLOCATE( z_rbfmat, z_diag, z_rbfval, STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_grf_interpolation:rbf_compute_coeff_grf',      &
      &             'deallocation for working arrays failed')
    ENDIF

!$OMP END PARALLEL
  ENDDO CE_LOOP

! Optional debug output for RBF coefficients
#ifdef DEBUG_COEFF
  i_startblk = ptr_patch(jg)%edges%start_blk(grf_bdyintp_start_e,jcd)
  i_endblk   = ptr_patch(jg)%edges%end_blk(min_rledge_int,jcd)

  DO jb =  i_startblk, i_endblk

    CALL get_indices_e(ptr_patch(jg), jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, grf_bdyintp_start_e, min_rledge_int)

    DO je = i_startidx, i_endidx

      istencil(je) = ptr_grf%grf_vec_stencil_1a(je,jb)
      write(570+jg,'(2i5,25f12.6)') jb,je,ptr_grf%grf_vec_coeff_1a(1:istencil(je),je,jb)
      istencil(je) = ptr_grf%grf_vec_stencil_1b(je,jb)
      write(570+jg,'(2i5,25f12.6)') jb,je,ptr_grf%grf_vec_coeff_1b(1:istencil(je),je,jb)
      istencil(je) = ptr_grf%grf_vec_stencil_2a(je,jb)
      write(570+jg,'(2i5,25f12.6)') jb,je,ptr_grf%grf_vec_coeff_2a(1:istencil(je),je,jb)
      IF (ptr_ep%refin_ctrl(je,jb) == -1) CYCLE
      istencil(je) = ptr_grf%grf_vec_stencil_2b(je,jb)
      write(570+jg,'(2i5,25f12.6)') jb,je,ptr_grf%grf_vec_coeff_2b(1:istencil(je),je,jb)

    END DO
  END DO
#endif

 END DO CD_LOOP
#ifdef DEBUG_COEFF
  CLOSE (570+jg)
#endif
END DO LEV_LOOP

END SUBROUTINE rbf_compute_coeff_grf_e

!-------------------------------------------------------------------------
!
!
!>
!! This routine computes the coefficients needed for inverse-distance-weighting.
!!
!! This routine computes the coefficients needed for inverse-distance-weighting
!! (IDW) interpolation to child cells and edges, which is required along the
!! lateral boundaries of nested model domains. Unlike RBF interpolation, the
!! IDW interpolation assumes a locally plain surface tangential to the sphere at
!! the target point of the interpolation (i.e. there is no transformation to
!! 3D cartesian space).
!!
!! @par Revision History
!! Developed and tested  by Guenther Zaengl (May 2008)
!!
SUBROUTINE idw_compute_coeff_grf_e ( ptr_grf_state)
!

TYPE(t_gridref_state), TARGET, INTENT(inout) :: ptr_grf_state(n_dom_start:)

TYPE(t_patch),      POINTER :: p_pp => NULL()
TYPE(t_patch),      POINTER :: p_pc => NULL()

TYPE(t_grid_edges),POINTER  :: ptr_ep, ptr_ec ! pointer to parent and child edges

TYPE(t_gridref_single_state),POINTER :: ptr_grf ! pointer to gridref_state for a
                                              ! single child domain

                                                             ! coordinates of ...
TYPE(t_cartesian_coordinates) :: cc_e1, cc_e, cc_childedge    ! edge midpoints

REAL(wp)           :: z_nx1(3),z_nx2(3)     ! 3d  normal velocity
                                            ! vectors at edge midpoints

REAL(wp)           :: z_dist, z_mindist     ! distance/minimum distance between data points

REAL(wp)           :: z_nxprod              ! scalar product of normal
                                            ! velocity vectors
REAL(wp)           :: z_wgtsum              ! sum of weighting factors (for normalization)

REAL(wp) :: z_idwwgt(6)      ! IDW weighting factors

! Index variables
INTEGER :: jg, jb, je, jcd, jgc, &
           iie1, ibe1, je1, iiec, ibec
INTEGER :: istencil              ! number of edges for the stencil
INTEGER :: ist                   ! status variable
INTEGER :: i_startblk                ! start block
INTEGER :: i_endblk                  ! end index
INTEGER :: i_startidx                ! start index
INTEGER :: i_endidx                  ! end index

TYPE(t_patch), POINTER      :: ptr_patch(:)
!--------------------------------------------------------------------

ptr_patch => p_patch

LEV_LOOP: DO jg = n_dom_start, n_dom-1

 CD_LOOP: DO jcd = 1, ptr_patch(jg)%n_childdom

  jgc  = ptr_patch(jg)%child_id(jcd)
  p_pc => p_patch(jgc)

  p_pp    => p_patch_local_parent(jgc)
  ptr_grf => p_grf_state_local_parent(jgc)%p_dom(jcd)

  ptr_ep   => p_pp%edges
  ptr_ec   => p_pc%edges

  ! Start and end blocks for which vector interpolation is needed
  i_startblk = p_pp%edges%start_blk(grf_bdyintp_start_e,jcd)
  i_endblk   = p_pp%edges%end_blk(min_rledge_int,jcd)

!$OMP PARALLEL PRIVATE (z_idwwgt,istencil,ist,jb, i_startidx, i_endidx)
  DO jb =  i_startblk, i_endblk

    CALL get_indices_e(p_pp, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, grf_bdyintp_start_e, min_rledge_int)

! Note: OMP parallelization is done over je because for reasonable (=efficient)
! choices of nproma, the boundary interpolation zone extends over no more
! than 2 blocks.
! Note: OMP directives do not work in combination with subroutine inlining on the NEC
!$OMP DO PRIVATE (je,je1,iie1,ibe1,cc_e,cc_e1,                       &
!$OMP    z_mindist,z_nx1,z_wgtsum,      &
!$OMP    z_nx2,z_nxprod,z_dist,iiec,ibec,cc_childedge) ICON_OMP_DEFAULT_SCHEDULE
    DO je = i_startidx, i_endidx

      ! 1a) child edge 1
      !
      ! get current number of points for stencil and
      ! allocate temporary arrays
      istencil = ptr_grf%grf_vec_stencil_1a(je,jb)

     ! compute minimum distance between local point and other stencil points
      cc_e = gc2cc(ptr_ep%center(je,jb))
      z_mindist = 1.e38_wp
      DO je1 = 1, istencil
        iie1 = ptr_grf%grf_vec_ind_1a(je,je1,jb)
        ibe1 = ptr_grf%grf_vec_blk_1a(je,je1,jb)
        IF ((je /= iie1) .OR. (jb /= ibe1)) THEN
          cc_e1 = gc2cc(ptr_ep%center(iie1,ibe1))
          z_mindist = MIN(z_mindist,arc_length(cc_e1,cc_e))
        ENDIF
      END DO

      !  compute the vector IDW interpolation coefficients

      iiec = ptr_ep%child_idx(je,jb,1) ! child edge 1
      ibec = ptr_ep%child_blk(je,jb,1)
      cc_childedge = gc2cc(ptr_ec%center(iiec,ibec))

      z_nx1(:) = ptr_ec%primal_cart_normal(iiec,ibec)%x(:)

      z_wgtsum = 0.0_wp

      DO je1 = 1, istencil
        !
        ! indices for each edge je1 of IDW stencil
        !
        iie1 = ptr_grf%grf_vec_ind_1a(je,je1,jb)
        ibe1 = ptr_grf%grf_vec_blk_1a(je,je1,jb)
        !
        ! transform to cartesian vector
        !
        cc_e1 = gc2cc(ptr_ep%center(iie1,ibe1))
        z_dist = arc_length(cc_e1,cc_childedge)/z_mindist

        !
        ! normal component and coordinates of edge midpoints for each edge je1
        ! of IDW stencil
        ! transform to cartesian vectors
        ! compute unit normal to edge je1
        z_nx2(:) = ptr_ep%primal_cart_normal(iie1,ibe1)%x(:)

        z_nxprod = DOT_PRODUCT(z_nx1,z_nx2)

        z_idwwgt(je1) = z_nxprod * invdwgt(z_dist,grf_idw_exp_e12)
        z_wgtsum = z_wgtsum + z_idwwgt(je1)*z_nxprod

      END DO

      ! Store normalized weighting factors in output array
      DO je1 = 1, istencil
        ptr_grf%grf_vec_coeff_1a(je1,je,jb) = z_idwwgt(je1)/z_wgtsum
      END DO

      ! 1a) child edge 2
      !
      ! get current number of points for stencil and
      ! allocate temporary arrays
      istencil = ptr_grf%grf_vec_stencil_1b(je,jb)

      ! compute minimum distance between local point and other stencil points
      cc_e = gc2cc(ptr_ep%center(je,jb))
      z_mindist = 1.e38_wp
      DO je1 = 1, istencil
        iie1 = ptr_grf%grf_vec_ind_1b(je,je1,jb)
        ibe1 = ptr_grf%grf_vec_blk_1b(je,je1,jb)
        IF ((je /= iie1) .OR. (jb /= ibe1)) THEN
          cc_e1 = gc2cc(ptr_ep%center(iie1,ibe1))
          z_mindist = MIN(z_mindist,arc_length(cc_e1,cc_e))
        ENDIF
      END DO

      !  compute the vector IDW interpolation coefficients

      iiec = ptr_ep%child_idx(je,jb,2) ! child edge 2
      ibec = ptr_ep%child_blk(je,jb,2)
      cc_childedge = gc2cc(ptr_ec%center(iiec,ibec))

      z_nx1(:) = ptr_ec%primal_cart_normal(iiec,ibec)%x(:)

      z_wgtsum = 0.0_wp

      DO je1 = 1, istencil
        !
        ! indices for each edge je1 of IDW stencil
        !
        iie1 = ptr_grf%grf_vec_ind_1b(je,je1,jb)
        ibe1 = ptr_grf%grf_vec_blk_1b(je,je1,jb)
        !
        ! transform to cartesian vector
        !
        cc_e1 = gc2cc(ptr_ep%center(iie1,ibe1))
        z_dist = arc_length(cc_e1,cc_childedge)/z_mindist

        !
        ! normal component and coordinates of edge midpoints for each edge je1
        ! of IDW stencil
        ! transform to cartesian vectors
        ! compute unit normal to edge je1
        z_nx2(:) = ptr_ep%primal_cart_normal(iie1,ibe1)%x(:)

        z_nxprod = DOT_PRODUCT(z_nx1,z_nx2)

        z_idwwgt(je1) = z_nxprod * invdwgt(z_dist,grf_idw_exp_e12)
        z_wgtsum = z_wgtsum + z_idwwgt(je1)*z_nxprod

      END DO

      ! Store normalized weighting factors in output array
      DO je1 = 1, istencil
        ptr_grf%grf_vec_coeff_1b(je1,je,jb) = z_idwwgt(je1)/z_wgtsum
      END DO

      ! 1c) child edge 3
      !
      ! get current number of points for stencil and
      ! allocate temporary arrays
      istencil = ptr_grf%grf_vec_stencil_2a(je,jb)

      ! compute minimum distance between local point and other stencil points
      cc_e = gc2cc(ptr_ep%center(je,jb))
      z_mindist = 1.e38_wp
      DO je1 = 1, istencil
        iie1 = ptr_grf%grf_vec_ind_2a(je,je1,jb)
        ibe1 = ptr_grf%grf_vec_blk_2a(je,je1,jb)
        IF ((je /= iie1) .OR. (jb /= ibe1)) THEN
          cc_e1 = gc2cc(ptr_ep%center(iie1,ibe1))
          z_mindist = MIN(z_mindist,arc_length(cc_e1,cc_e))
        ENDIF
      END DO

      !  compute the vector IDW interpolation coefficients

      iiec = ABS(ptr_ep%child_idx(je,jb,3)) ! child edge 3
      ibec = ptr_ep%child_blk(je,jb,3)
      cc_childedge = gc2cc(ptr_ec%center(iiec,ibec))

      z_nx1(:) = ptr_ec%primal_cart_normal(iiec,ibec)%x(:)

      z_wgtsum = 0.0_wp

      DO je1 = 1, istencil
        !
        ! indices for each edge je1 of IDW stencil
        !
        iie1 = ptr_grf%grf_vec_ind_2a(je,je1,jb)
        ibe1 = ptr_grf%grf_vec_blk_2a(je,je1,jb)
        !
        ! transform to cartesian vector
        !
        cc_e1 = gc2cc(ptr_ep%center(iie1,ibe1))
        z_dist = arc_length(cc_e1,cc_childedge)/z_mindist

        !
        ! normal component and coordinates of edge midpoints for each edge je1
        ! of IDW stencil
        ! transform to cartesian vectors
        ! compute unit normal to edge je1
        z_nx2(:) = ptr_ep%primal_cart_normal(iie1,ibe1)%x(:)

        z_nxprod = DOT_PRODUCT(z_nx1,z_nx2)

        z_idwwgt(je1) = z_nxprod * invdwgt(z_dist,grf_idw_exp_e34)
        z_wgtsum = z_wgtsum + z_idwwgt(je1)*z_nxprod

      END DO

      ! Store normalized weighting factors in output array
      DO je1 = 1, istencil
        ptr_grf%grf_vec_coeff_2a(je1,je,jb) = z_idwwgt(je1)/z_wgtsum
      END DO

      ! 1d) child edge 4
      IF (ptr_ep%refin_ctrl(je,jb) == -1) CYCLE
      !
      ! get current number of points for stencil and
      ! allocate temporary arrays
      istencil = ptr_grf%grf_vec_stencil_2b(je,jb)

      ! compute minimum distance between local point and other stencil points
      cc_e = gc2cc(ptr_ep%center(je,jb))
      z_mindist = 1.e38_wp
      DO je1 = 1, istencil
        iie1 = ptr_grf%grf_vec_ind_2b(je,je1,jb)
        ibe1 = ptr_grf%grf_vec_blk_2b(je,je1,jb)
        IF ((je /= iie1) .OR. (jb /= ibe1)) THEN
          cc_e1 = gc2cc(ptr_ep%center(iie1,ibe1))
          z_mindist = MIN(z_mindist,arc_length(cc_e1,cc_e))
        ENDIF
      END DO

      !  compute the vector IDW interpolation coefficients

      iiec = ABS(ptr_ep%child_idx(je,jb,4)) ! child edge 4
      ibec = ptr_ep%child_blk(je,jb,4)
      cc_childedge = gc2cc(ptr_ec%center(iiec,ibec))

      z_nx1(:) = ptr_ec%primal_cart_normal(iiec,ibec)%x(:)

      z_wgtsum = 0.0_wp

      DO je1 = 1, istencil
        !
        ! indices for each edge je1 of IDW stencil
        !
        iie1 = ptr_grf%grf_vec_ind_2b(je,je1,jb)
        ibe1 = ptr_grf%grf_vec_blk_2b(je,je1,jb)
        !
        ! transform to cartesian vector
        !
        cc_e1 = gc2cc(ptr_ep%center(iie1,ibe1))
        z_dist = arc_length(cc_e1,cc_childedge)/z_mindist

        !
        ! normal component and coordinates of edge midpoints for each edge je1
        ! of IDW stencil
        ! transform to cartesian vectors
        ! compute unit normal to edge je1
        z_nx2(:) = ptr_ep%primal_cart_normal(iie1,ibe1)%x(:)

        z_nxprod = DOT_PRODUCT(z_nx1,z_nx2)

        z_idwwgt(je1) = z_nxprod * invdwgt(z_dist,grf_idw_exp_e34)
        z_wgtsum = z_wgtsum + z_idwwgt(je1)*z_nxprod

      END DO

      ! Store normalized weighting factors in output array
      DO je1 = 1, istencil
        ptr_grf%grf_vec_coeff_2b(je1,je,jb) = z_idwwgt(je1)/z_wgtsum
      END DO

    ENDDO
!$OMP END DO NOWAIT
  ENDDO ! blocks

!$OMP END PARALLEL

! Optional debug output for IDW coefficients
#ifdef DEBUG_COEFF
  DO jb =  i_startblk, i_endblk

    CALL get_indices_e(ptr_patch(jg), jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, grf_bdyintp_start_e, min_rledge_int)

    DO je = i_startidx, i_endidx

      istencil = ptr_grf%grf_vec_stencil_1a(je,jb)
      write(610+jg,'(2i5,25f12.6)') jb,je,ptr_grf%grf_vec_coeff_1a(1:istencil,je,jb)
      istencil = ptr_grf%grf_vec_stencil_1b(je,jb)
      write(610+jg,'(2i5,25f12.6)') jb,je,ptr_grf%grf_vec_coeff_1b(1:istencil,je,jb)
      istencil = ptr_grf%grf_vec_stencil_2a(je,jb)
      write(610+jg,'(2i5,25f12.6)') jb,je,ptr_grf%grf_vec_coeff_2a(1:istencil,je,jb)
      IF (ptr_ep%refin_ctrl(je,jb) == -1) CYCLE
      istencil = ptr_grf%grf_vec_stencil_2b(je,jb)
      write(610+jg,'(2i5,25f12.6)') jb,je,ptr_grf%grf_vec_coeff_2b(1:istencil,je,jb)

    END DO
  END DO
#endif

 END DO CD_LOOP
#ifdef DEBUG_COEFF
  CLOSE (610+jg)
#endif
END DO LEV_LOOP

END SUBROUTINE idw_compute_coeff_grf_e


!-------------------------------------------------------------------------
END MODULE mo_grf_intp_coeffs
