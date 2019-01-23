!>
!! Contains averaging and interpolation routines (ICON <--> FEM)
!!
!! @par Revision History
!! Based on: *list modules which were used for this development*
!! Developed  by Einar (2014), Vladimir (2015)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_ice_fem_interpolation
  !-------------------------------------------------------------------------
  !
  USE mo_kind,                ONLY: wp
  USE mo_parallel_config,     ONLY: nproma

  USE mo_impl_constants,      ONLY: min_rlvert
  USE mo_run_config,          ONLY: ltimer
  USE mo_loopindices,         ONLY: get_indices_v
  USE mo_timer,               ONLY: timer_start, timer_stop, timer_intp

  USE mo_model_domain,        ONLY: t_patch, t_patch_3D
  USE mo_math_types,          ONLY: t_cartesian_coordinates
  USE mo_math_utilities,      ONLY: cc_norm, gvec2cvec, cvec2gvec
  USE mo_grid_subset,         ONLY: t_subset_range, get_index_range
  USE mo_impl_constants,      ONLY: sea_boundary

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: gvec2cvec_c_2d
  PUBLIC :: cvec2gvec_c_2d
  PUBLIC :: rotate_cvec_v
  PUBLIC :: cvec2gvec_v_fem
  PUBLIC :: gvec2cvec_v_fem
  PUBLIC :: map_edges2verts
  PUBLIC :: map_verts2edges
  PUBLIC :: cells2verts_scalar_seaice

CONTAINS

  !-------------------------------------------------------------------------
  !
  !> Convert to cartesian coordinates lat-lon velocity vector on cells centers
  !!
  !! @par Revision History
  !! Developed by Vladimir Lapin, MPI-M (2015-10-13)
  !<Optimize:inUse>
  SUBROUTINE gvec2cvec_c_2d(patch_3d, gvec_u, gvec_v, cvec)

    TYPE(t_patch_3d),TARGET, INTENT(in)       :: patch_3d
    REAL(wp), INTENT(in)                      :: gvec_u(:,:), gvec_v(:,:)
    TYPE(t_cartesian_coordinates),INTENT(out) :: cvec(nproma,patch_3d%p_patch_2D(1)%alloc_cell_blocks)

   ! Local variables
    ! Patch and ranges
    TYPE(t_patch), POINTER :: p_patch
    TYPE(t_subset_range), POINTER :: all_cells

    ! Indexing
    INTEGER  :: i_startidx_c, i_endidx_c, jc, jb!, jk
    !-----------------------------------------------------------------------

    p_patch   => patch_3d%p_patch_2d(1)
    all_cells => p_patch%cells%all

!ICON_OMP_PARALLEL_DO PRIVATE(i_startidx_c,i_endidx_c, jc) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
      DO jc = i_startidx_c, i_endidx_c
        IF(patch_3d%lsm_c(jc,1,jb) <= sea_boundary)THEN
          CALL gvec2cvec(  gvec_u(jc,jb), gvec_v(jc,jb), &
                         & p_patch%cells%center(jc,jb)%lon,     &
                         & p_patch%cells%center(jc,jb)%lat,     &
                         & cvec(jc,jb)%x(1),cvec(jc,jb)%x(2),cvec(jc,jb)%x(3))
        ELSE
          cvec(jc,jb)%x    = 0.0_wp
        ENDIF
      END DO
    END DO
!ICON_OMP_END_PARALLEL_DO

  END SUBROUTINE gvec2cvec_c_2d
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !
  !> Inverse of gvec2cvec_c_2d. Convert cc vector on cell centers to lat-lon
  !!
  !! @par Revision History
  !! Developed by Vladimir Lapin, MPI-M (2015-10-13)
  !<Optimize:inUse>
  SUBROUTINE cvec2gvec_c_2d(patch_3d, cvec, gvec_u, gvec_v)

    TYPE(t_patch_3d),TARGET, INTENT(in)       :: patch_3d
    TYPE(t_cartesian_coordinates),INTENT(in)  :: cvec(nproma,patch_3d%p_patch_2D(1)%alloc_cell_blocks)
    REAL(wp), INTENT(out)                     :: gvec_u(nproma,patch_3d%p_patch_2D(1)%alloc_cell_blocks), &
                                               & gvec_v(nproma,patch_3d%p_patch_2D(1)%alloc_cell_blocks)

   ! Local variables
    ! Patch and ranges
    TYPE(t_patch), POINTER :: p_patch
    TYPE(t_subset_range), POINTER :: all_cells

    ! Indexing
    INTEGER  :: i_startidx_c, i_endidx_c, jc, jb!, jk
    !-----------------------------------------------------------------------

    p_patch   => patch_3d%p_patch_2d(1)
    all_cells => p_patch%cells%all

!ICON_OMP_PARALLEL_DO PRIVATE(i_startidx_c,i_endidx_c, jc) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
      DO jc = i_startidx_c, i_endidx_c
        IF(patch_3d%lsm_c(jc,1,jb) <= sea_boundary)THEN
          CALL cvec2gvec(  cvec(jc,jb)%x(1),cvec(jc,jb)%x(2),cvec(jc,jb)%x(3), &
                         & p_patch%cells%center(jc,jb)%lon,     &
                         & p_patch%cells%center(jc,jb)%lat,     &
                         & gvec_u(jc,jb), gvec_v(jc,jb))
        ELSE
          gvec_u(jc,jb) = 0._wp
          gvec_v(jc,jb) = 0._wp
        ENDIF
      END DO
    END DO
!ICON_OMP_END_PARALLEL_DO

  END SUBROUTINE cvec2gvec_c_2d
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !
  !> Rotate cartesian velocity vector on verts to the rotated FEM grid
  !!
  !! @par Revision History
  !! Developed by Vladimir Lapin, MPI-M (2015-10-25)
  !<Optimize:inUse>
  SUBROUTINE rotate_cvec_v(p_patch, cvec_in, rot_mat_3D, cvec_out)

    TYPE(t_patch), TARGET, INTENT(in)         :: p_patch
    TYPE(t_cartesian_coordinates),INTENT(in)  :: cvec_in (nproma,p_patch%nblks_v)
    REAL(wp)                                  :: rot_mat_3D(3,3)
    TYPE(t_cartesian_coordinates),INTENT(out) :: cvec_out(nproma,p_patch%nblks_v)

   ! Local variables
    TYPE(t_subset_range), POINTER :: all_verts

    ! Indexing
    INTEGER  :: i_startidx_v, i_endidx_v, jv, jb
    !-----------------------------------------------------------------------

    all_verts => p_patch%verts%all

!ICON_OMP_PARALLEL_DO PRIVATE(i_startidx_v,i_endidx_v, jv) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = all_verts%start_block, all_verts%end_block
      CALL get_index_range(all_verts, jb, i_startidx_v, i_endidx_v)
      DO jv = i_startidx_v, i_endidx_v
        ! Intrinsic function matmul not applied, due to poor performance.
        cvec_out(jv,jb)%x(1) = DOT_PRODUCT(rot_mat_3D(1,:),cvec_in(jv,jb)%x(:))
        cvec_out(jv,jb)%x(2) = DOT_PRODUCT(rot_mat_3D(2,:),cvec_in(jv,jb)%x(:))
        cvec_out(jv,jb)%x(3) = DOT_PRODUCT(rot_mat_3D(3,:),cvec_in(jv,jb)%x(:))
      END DO
    END DO
!ICON_OMP_END_PARALLEL_DO

  END SUBROUTINE rotate_cvec_v
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !
  !> Convert cartesian velocity vector to lat-lon vector (on the FEM grid)
  !!
  !! @par Revision History
  !! Developed by Vladimir Lapin, MPI-M (2015-10-25)
  !
  SUBROUTINE cvec2gvec_v_fem(p_patch, cvec, gvec_u, gvec_v)

    USE mo_ice_fem_mesh,           ONLY: coord_nod2D

    TYPE(t_patch), TARGET, INTENT(in)       :: p_patch
    TYPE(t_cartesian_coordinates),INTENT(in):: cvec(nproma,p_patch%nblks_v)
    REAL(wp), INTENT(out)                   :: gvec_u(:), gvec_v(:)

   ! Local variables
    TYPE(t_subset_range), POINTER :: all_verts

    ! Indexing
    INTEGER  :: i_startidx_v, i_endidx_v, jv, jb, jk
    !-----------------------------------------------------------------------

    all_verts => p_patch%verts%all

    jk=0
    DO jb = all_verts%start_block, all_verts%end_block
      CALL get_index_range(all_verts, jb, i_startidx_v, i_endidx_v)
      DO jv = i_startidx_v, i_endidx_v
        jk=jk+1
          CALL cvec2gvec(cvec(jv,jb)%x(1), cvec(jv,jb)%x(2), cvec(jv,jb)%x(3), &
                       & coord_nod2D(1,jk), coord_nod2D(2,jk), & ! lon, lat
                       & gvec_u(jk), gvec_v(jk))
      END DO
    END DO

  END SUBROUTINE cvec2gvec_v_fem
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !
  !> Inverse of cvec2gvec_v_fem. Convert lat-lon vector to cartesian (on the FEM grid)
  !!
  !! @par Revision History
  !! Developed by Vladimir Lapin, MPI-M (2015-10-25)
  !
  SUBROUTINE gvec2cvec_v_fem(p_patch, gvec_u, gvec_v, cvec)

    USE mo_ice_fem_mesh,           ONLY: coord_nod2D

    TYPE(t_patch), TARGET, INTENT(in)        :: p_patch
    REAL(wp), INTENT(in)                     :: gvec_u(:), gvec_v(:)
    TYPE(t_cartesian_coordinates),INTENT(out):: cvec(nproma,p_patch%nblks_v)

    ! Local variables
    TYPE(t_subset_range), POINTER :: all_verts

    ! Indexing
    INTEGER  :: i_startidx_v, i_endidx_v, jv, jb, jk
    !-----------------------------------------------------------------------

    all_verts => p_patch%verts%all

    jk=0
    DO jb = all_verts%start_block, all_verts%end_block
      CALL get_index_range(all_verts, jb, i_startidx_v, i_endidx_v)
      DO jv = i_startidx_v, i_endidx_v
        jk=jk+1
          CALL gvec2cvec(  gvec_u(jk), gvec_v(jk),                   &
                         & coord_nod2D(1,jk), coord_nod2D(2,jk), & ! lon, lat
                         & cvec(jv,jb)%x(1), cvec(jv,jb)%x(2), cvec(jv,jb)%x(3) )
      END DO
    END DO

  END SUBROUTINE gvec2cvec_v_fem
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !
  !> Map vectors from edges to vertices
  !! Based on map_edges2vert_3d in ocean/math/mo_ocean_math_operators.f90
  !!      and edges2verts_scalar in shr_horizontal/mo_icon_interpolation_scalar.f90
  !!
  !! @par Revision History
  !! Developed by Einar Olason, MPI-M (2013-08-05)
  !! Modified by Vladimir Lapin, MPI-M (2015-10-13)
  !<Optimize:inUse>
  SUBROUTINE map_edges2verts(p_patch, vn, edge2vert_coeff_cc, p_vn_dual)

    TYPE(t_patch), TARGET, INTENT(in)         :: p_patch
    REAL(wp), INTENT(in)                      :: vn(:,:)
    TYPE(t_cartesian_coordinates),INTENT(in)  :: edge2vert_coeff_cc(:,:,:,:)
    TYPE(t_cartesian_coordinates),INTENT(out) :: p_vn_dual(:,:)

    ! Local variables
    TYPE(t_subset_range), POINTER :: verts_in_domain

    ! Indexing
    INTEGER :: jv, jb,jev
    INTEGER :: ile, ibe
    INTEGER :: i_startidx_v, i_endidx_v

    !-----------------------------------------------------------------------

    verts_in_domain => p_patch%verts%in_domain
    ! Set to zero for nag compiler
    p_vn_dual(:,:)%x(1) = 0._wp
    p_vn_dual(:,:)%x(2) = 0._wp
    p_vn_dual(:,:)%x(3) = 0._wp

!ICON_OMP_PARALLEL_DO PRIVATE(i_startidx_v, i_endidx_v, jv, jev, ile, ibe) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = verts_in_domain%start_block, verts_in_domain%end_block
      CALL get_index_range(verts_in_domain, jb, i_startidx_v, i_endidx_v)
        DO jv = i_startidx_v, i_endidx_v

          p_vn_dual(jv,jb)%x = 0.0_wp
          DO jev = 1, p_patch%verts%num_edges(jv,jb)

            ! get line and block indices of edge jev around vertex jv
            ile = p_patch%verts%edge_idx(jv,jb,jev)
            ibe = p_patch%verts%edge_blk(jv,jb,jev)

        ! Sea-land boundary is taken into account by coeffcients.
            p_vn_dual(jv,jb)%x = p_vn_dual(jv,jb)%x + edge2vert_coeff_cc(jv,1,jb,jev)%x * vn(ile,ibe)

        END DO
      END DO
    END DO
!ICON_OMP_END_PARALLEL_DO

  END SUBROUTINE map_edges2verts
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !
  !> Map vectors from vertices to edges
  !! Based on ideas from rot_vertex_ocean_3d in ocean/math/mo_ocean_math_operators.f90
  !!
  !! @par Revision History
  !! Developed by Vladimir Lapin, MPI-M (2015-08-13)
  !<Optimize:inUse>
  SUBROUTINE map_verts2edges(p_patch, p_vn_dual, edge2vert_coeff_cc_t, vn)

    TYPE(t_patch), TARGET, INTENT(in)        :: p_patch
    TYPE(t_cartesian_coordinates),INTENT(in) :: p_vn_dual(:,:) !(nproma,p_patch%nblks_v)
    TYPE(t_cartesian_coordinates),INTENT(in) :: edge2vert_coeff_cc_t(:,:,:,:)

    REAL(wp), INTENT(inout)                  :: vn(:,:) !(nproma,p_patch%nblks_e)

    ! Local variables
    TYPE(t_subset_range), POINTER :: edges_in_domain

    ! Indexing
    INTEGER :: edge_index, edge_block
    INTEGER :: il_v1, il_v2,ib_v1, ib_v2
    INTEGER :: start_index_e, end_index_e

    TYPE(t_cartesian_coordinates)   :: p_vn_dual_e

    !-----------------------------------------------------------------------
    edges_in_domain => p_patch%edges%in_domain

    ! loop through all edges and add contribution from neighboring vertices

!ICON_OMP_PARALLEL_DO PRIVATE(start_index_e,end_index_e,edge_index, &
!ICON_OMP  il_v1,ib_v1,il_v2,ib_v2,p_vn_dual_e) ICON_OMP_DEFAULT_SCHEDULE
    DO edge_block = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, edge_block, start_index_e, end_index_e)
        DO edge_index = start_index_e, end_index_e

            !!--------------------------------------------------------------------
            !! Sea-land boundary is taken into account by coeffcients.
            !! It is also assumed here that p_vn_dual is already zero where it should be zero.
            !!--------------------------------------------------------------------
              ! two neighboring verts for each edge
              il_v1 = p_patch%edges%vertex_idx(edge_index,edge_block,1)
              ib_v1 = p_patch%edges%vertex_blk(edge_index,edge_block,1)
              il_v2 = p_patch%edges%vertex_idx(edge_index,edge_block,2)
              ib_v2 = p_patch%edges%vertex_blk(edge_index,edge_block,2)

              ! full cartesian velocity at the edge center
              p_vn_dual_e%x = cc_norm(edge2vert_coeff_cc_t(edge_index,1,edge_block,1))*p_vn_dual(il_v1,ib_v1)%x &
                  & + cc_norm(edge2vert_coeff_cc_t(edge_index,1,edge_block,2))*p_vn_dual(il_v2,ib_v2)%x
              ! project to get the normal component only
              vn(edge_index,edge_block) = DOT_PRODUCT(p_vn_dual_e%x, p_patch%edges%primal_cart_normal(edge_index,edge_block)%x)
        END DO
    END DO
!ICON_OMP_END_PARALLEL_DO

  END SUBROUTINE map_verts2edges
  !-------------------------------------------------------------------------

!------------------------------------------------------------------------
!>
!! Computes  average of scalar fields from centers of cells to vertices.
!! Based on edges2verts_scalar in shr_horizontal/mo_icon_interpolation_scalar.f90
!!
!! Sea-ice module used to call this funciton instead of the standard cells2verts_scalar
!! because in mo_ocean_nml_crosscheck <use_duplicated_connectivity> was set to .FALSE.
!! The modified version includes an ad-hoc (slow!) fix that checks for zero indices in ptr_patch%verts%cell_idx.
!! By default empty indices (e.g. 6th vertex in pentagons) are replaced by last non-zero values when
!! use_duplicated_connectivity = .TRUE. See CALL move_dummies_to_end_idxblk in subroutine complete_patches
!!
!! Usage is depriciated.
!!
!! @par Revision History
!! Added by Vladimir Lapin, MPI-M (2015-08-13)

!!
SUBROUTINE cells2verts_scalar_seaice( p_cell_in, ptr_patch, c_int, p_vert_out,  &
  &                            opt_slev, opt_elev, opt_rlstart, opt_rlend )
!

TYPE(t_patch), TARGET, INTENT(in) :: ptr_patch

! cell based scalar input field
REAL(wp), INTENT(in) :: p_cell_in(:,:,:)   ! dim: (nproma,nlev,nblks_c)

! coefficients for interpolation
REAL(wp), INTENT(in) :: c_int(:,:,:)       ! dim: (nproma,9-cell_type,nblks_v)

INTEGER, INTENT(in), OPTIONAL :: opt_slev  ! optional vertical start level

INTEGER, INTENT(in), OPTIONAL :: opt_elev  ! optional vertical end level

! start and end values of refin_ctrl flag
INTEGER, INTENT(in), OPTIONAL ::  opt_rlstart, opt_rlend

! vertex based scalar output field
REAL(wp), INTENT(inout) :: p_vert_out(:,:,:) ! dim: (nproma,nlev,nblks_v)

INTEGER :: slev, elev     ! vertical start and end level
INTEGER :: jv, jk, jb, ji
INTEGER :: rl_start, rl_end
INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx, i_nchdom
INTEGER :: cell_index, cell_block

INTEGER,  DIMENSION(:,:,:),   POINTER :: iidx, iblk

!-----------------------------------------------------------------------

! check optional arguments
IF ( PRESENT(opt_slev) ) THEN
  slev = opt_slev
ELSE
  slev = 1
END IF
IF ( PRESENT(opt_elev) ) THEN
  elev = opt_elev
ELSE
  elev = UBOUND(p_cell_in,2)
END IF

IF ( PRESENT(opt_rlstart) ) THEN
  rl_start = opt_rlstart
ELSE
  rl_start = 2
END IF
IF ( PRESENT(opt_rlend) ) THEN
  rl_end = opt_rlend
ELSE
  rl_end = min_rlvert
END IF

iidx => ptr_patch%verts%cell_idx
iblk => ptr_patch%verts%cell_blk

! values for the blocking
i_nchdom   = MAX(1,ptr_patch%n_childdom)
i_startblk = ptr_patch%verts%start_blk(rl_start,1)
i_endblk   = ptr_patch%verts%end_blk(rl_end,i_nchdom)


IF (ltimer) CALL timer_start(timer_intp)

IF (ptr_patch%geometry_info%cell_type == 6) THEN
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jv,jk) ICON_OMP_DEFAULT_SCHEDULE
  DO jb = i_startblk, i_endblk

    CALL get_indices_v(ptr_patch, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, rl_start, rl_end)


#ifdef __LOOP_EXCHANGE
    DO jv = i_startidx, i_endidx
      DO jk = slev, elev
#else
!CDIR UNROLL=6
    DO jk = slev, elev
      DO jv = i_startidx, i_endidx
#endif

         p_vert_out(jv,jk,jb) =                       &
           c_int(jv,1,jb) * p_cell_in(iidx(jv,jb,1),jk,iblk(jv,jb,1)) + &
           c_int(jv,2,jb) * p_cell_in(iidx(jv,jb,2),jk,iblk(jv,jb,2)) + &
           c_int(jv,3,jb) * p_cell_in(iidx(jv,jb,3),jk,iblk(jv,jb,3))

      ENDDO
    ENDDO

  ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
ELSE IF (ptr_patch%geometry_info%cell_type == 3) THEN
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jv,jk) ICON_OMP_DEFAULT_SCHEDULE
  DO jb = i_startblk, i_endblk

    CALL get_indices_v(ptr_patch, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, rl_start, rl_end)


#ifdef __LOOP_EXCHANGE
    DO jv = i_startidx, i_endidx
      DO jk = slev, elev
#else
!CDIR UNROLL=6
    DO jk = slev, elev
      DO jv = i_startidx, i_endidx
#endif

        p_vert_out(jv,jk,jb) = 0.0_wp

        DO ji = 1, 6
          cell_index = iidx(jv,jb,ji)
          cell_block = iblk(jv,jb,ji)
          IF (cell_index > 0)                                      &
            & p_vert_out(jv,jk,jb) = p_vert_out(jv,jk,jb) +                                        &
            &      c_int(jv,ji,jb) * p_cell_in(cell_index,jk,cell_block)
        ENDDO

!       Optimized version of the code that relies on non-zero values in iidx => ptr_patch%verts%cell_idx
!         p_vert_out(jv,jk,jb) =                       &
!           c_int(jv,1,jb) * p_cell_in(iidx(jv,jb,1),jk,iblk(jv,jb,1)) + &
!           c_int(jv,2,jb) * p_cell_in(iidx(jv,jb,2),jk,iblk(jv,jb,2)) + &
!           c_int(jv,3,jb) * p_cell_in(iidx(jv,jb,3),jk,iblk(jv,jb,3)) + &
!           c_int(jv,4,jb) * p_cell_in(iidx(jv,jb,4),jk,iblk(jv,jb,4)) + &
!           c_int(jv,5,jb) * p_cell_in(iidx(jv,jb,5),jk,iblk(jv,jb,5)) + &
!           c_int(jv,6,jb) * p_cell_in(iidx(jv,jb,6),jk,iblk(jv,jb,6))

      ENDDO
    ENDDO

  ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
ENDIF

IF (ltimer) CALL timer_stop(timer_intp)


END SUBROUTINE cells2verts_scalar_seaice
!------------------------------------------------------------------------

END MODULE mo_ice_fem_interpolation
