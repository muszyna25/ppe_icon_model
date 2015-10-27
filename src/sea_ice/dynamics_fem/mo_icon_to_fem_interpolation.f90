!>
!! Contains the implementation of interpolation from ICON grid to FEM grid for the sea-ice dynamics
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
!!

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_icon_to_fem_interpolation
  !-------------------------------------------------------------------------
  !
  USE mo_kind,                ONLY: wp
  USE mo_parallel_config,     ONLY: nproma
  USE mo_exception,           ONLY: finish,message

  USE mo_impl_constants,      ONLY: min_rlcell, min_rledge, min_rlvert
!  USE mo_grid_config,         ONLY: l_limited_area
  USE mo_run_config,          ONLY: ltimer
  USE mo_loopindices,         ONLY: get_indices_v!, get_indices_c, get_indices_e
  USE mo_timer,               ONLY: timer_start, timer_stop, timer_intp

  USE mo_model_domain,        ONLY: t_patch, t_patch_3D
  USE mo_math_utilities,      ONLY: t_cartesian_coordinates, cc_norm
  USE mo_grid_subset,         ONLY: t_subset_range, get_index_range
  USE mo_impl_constants,      ONLY: sea_boundary

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: map_edges2verts
  PUBLIC :: map_verts2edges_einar
  PUBLIC :: map_verts2edges
  PUBLIC :: cells2verts_scalar_seaice

CONTAINS

!-----------------------------------------------------------------------
!
!  ! averaging and interpolation routines and
!  ! routines needed to compute the coefficients therein
!
!-----------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !
  !> Map vectors from edges to vertices
  !! Based on map_edges2vert_3d in oce_dyn_icohom/mo_oce/math_operators.f90
  !!      and edges2verts_scalar in shr_horizontal/mo_icon_interpolation_scalar.f90
  !!
  !! @par Revision History
  !! Developed by Einar Olason, MPI-M (2013-08-05)
  !! Modified by Vladimir Lapin, MPI-M (2015-10-13)
  !
  SUBROUTINE map_edges2verts(p_patch, vn, edge2vert_coeff_cc, p_vn_dual)

    TYPE(t_patch), TARGET, INTENT(in)         :: p_patch
    REAL(wp), INTENT(in)                      :: vn(:,:)
    TYPE(t_cartesian_coordinates),INTENT(in)  :: edge2vert_coeff_cc(:,:,:,:)
    TYPE(t_cartesian_coordinates),INTENT(out) :: p_vn_dual(nproma,p_patch%nblks_v)

    ! Local variables

    ! Sub-set
    TYPE(t_subset_range), POINTER :: verts_in_domain

    ! Indexing
    INTEGER :: vertexIndex, blockNo, vertexConnect
    INTEGER :: edgeOfVertex_index, edgeOfVertex_block
    INTEGER :: start_index_v, end_index_v
    !-----------------------------------------------------------------------

    verts_in_domain => p_patch%verts%in_domain

    ! Set to zero outside the loop. Otherwise arithmetic error with nag compiler
    p_vn_dual(:,:)%x(1) = 0.0_wp
    p_vn_dual(:,:)%x(2) = 0.0_wp
    p_vn_dual(:,:)%x(3) = 0.0_wp

!ICON_OMP_PARALLEL_DO PRIVATE(blockNo,start_index_v,end_index_v, vertexIndex, vertexConnect, &
!ICON_OMP edgeOfVertex_index, edgeOfVertex_block) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = verts_in_domain%start_block, verts_in_domain%end_block
      CALL get_index_range(verts_in_domain, blockNo, start_index_v, end_index_v)
!      p_vn_dual(:,blockNo)%x(1) = 0.0_wp
!      p_vn_dual(:,blockNo)%x(2) = 0.0_wp
!      p_vn_dual(:,blockNo)%x(3) = 0.0_wp
      DO vertexIndex = start_index_v, end_index_v

        DO vertexConnect = 1, p_patch%verts%num_edges(vertexIndex,blockNo)
         ! get line and block indices of edge vertexConnect around vertex vertexIndex
          edgeOfVertex_index = p_patch%verts%edge_idx(vertexIndex,blockNo,vertexConnect)
          edgeOfVertex_block = p_patch%verts%edge_blk(vertexIndex,blockNo,vertexConnect)

            p_vn_dual(vertexIndex,blockNo)%x = p_vn_dual(vertexIndex,blockNo)%x        &
              & + edge2vert_coeff_cc(vertexIndex,1,blockNo,vertexConnect)%x & ! level = 1
              & * vn(edgeOfVertex_index,edgeOfVertex_block)
        END DO
      END DO
    END DO
!ICON_OMP_END_PARALLEL_DO

  END SUBROUTINE map_edges2verts

!--------------------------------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !
  !> Map vectors from vertices to edges
  !! Based on ideas from rot_vertex_ocean_3d in oce_dyn_icohom/mo_ocean_math_operators.f90
  !!
  !! @par Revision History
  !! Developed by Vladimir Lapin, MPI-M (2015-08-13)
  !
  SUBROUTINE map_verts2edges(p_patch_3D, p_vn_dual, edge2vert_coeff_cc_t, vn)

    TYPE(t_patch_3D), TARGET, INTENT(in)     :: p_patch_3D
    TYPE(t_cartesian_coordinates),INTENT(in) :: p_vn_dual(:,:) !(nproma,patch_3D%p_patch_2D(1)%nblks_v)
    TYPE(t_cartesian_coordinates),INTENT(in) :: edge2vert_coeff_cc_t(:,:,:,:)

    REAL(wp), INTENT(inout)                  :: vn(:,:) !(nproma,patch_3D%p_patch_2D(1)%nblks_e)

    ! Local variables

    ! Patch and sub-set
    TYPE(t_subset_range), POINTER :: edges_in_domain
    TYPE(t_patch), POINTER        :: p_patch

    ! Indexing
    INTEGER :: edge_index, edge_block
    INTEGER :: il_v1, il_v2,ib_v1, ib_v2
    INTEGER :: start_index_e, end_index_e

    TYPE(t_cartesian_coordinates)   :: p_vn_dual_e

    !-----------------------------------------------------------------------

    p_patch         => p_patch_3D%p_patch_2D(1)
    edges_in_domain => p_patch%edges%in_domain

    ! loop through all edges and add contribution from neighboring vertices

!ICON_OMP_PARALLEL_DO PRIVATE(edge_block,start_index_e,end_index_e,edge_index, &
!ICON_OMP  il_v1,ib_v1,il_v2,ib_v2,p_vn_dual_e,vn) ICON_OMP_DEFAULT_SCHEDULE
    DO edge_block = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, edge_block, start_index_e, end_index_e)
        DO edge_index = start_index_e, end_index_e

            !!--------------------------------------------------------------------
            !! Sea-land boundary is taken into account by coeffcients.
            !! It is also assumed here that p_vn_dual is already zero where it should be zero.
            !!--------------------------------------------------------------------
!            IF(p_patch_3D%lsm_e(edge_index,1,edge_block) <= sea_boundary) THEN
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
!            ELSE
!              vn(edge_index,edge_block) = 0.0_wp
!            ENDIF
        END DO
    END DO
!ICON_OMP_END_PARALLEL_DO

  END SUBROUTINE map_verts2edges

!------------------------------------------------------------------------
!>
!!  Computes  average of scalar fields from centers of cells to vertices.
!!
!!
!! @par Revision History
!! Developed  by Almut Gassmann, MPI-M (2009-01-28)
!!
!!
!>
!! Computes  average of scalar fields from centers of cells to vertices.
!! Based on edges2verts_scalar in shr_horizontal/mo_icon_interpolation_scalar.f90
!!
!! Includes a check for pentagons (very slow!). Ad-hoc fix for sea ice dynamics.
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

    ! Sea-ice module calls this funciton to interpolate from ICON grid to FEM P1P1 grid
    ! An indexing error used to occur because of zeros in iidx(jv,jb,:) at the pentagons centers and boundary vertices

    ! In the atmosphere module this is done by calling this function somewhere at the initialization step
    !    CALL move_dummies_to_end_idxblk( &
    !        ptr_patch%verts%cell_idx(:,:,1:max_verts_connectivity), &
    !        ptr_patch%n_patch_verts, max_verts_connectivity, &
    !        use_duplicated_connectivity)

!    if ( (iidx(jv,jb,6)==0) .or. (iblk(jv,jb,6)==0) ) then
!        print *, 'lat, lon =', ptr_patch%verts%vertex(jv,jb)
!    !    print *, '6 interp coeff =', c_int(jv,1:6,jb)
!        print *, 'iidx(jv,jb,:)', iidx(jv,jb,:)
!        print *, 'iblk(jv,jb,:)', iblk(jv,jb,:)
!        print *, 'max_connectivity', ptr_patch%verts%max_connectivity
!        print *, 'ptr_patch%verts%num_edges(jv,jb)', ptr_patch%verts%num_edges(jv,jb)
!    endif

        p_vert_out(jv,jk,jb) = 0.0_wp

        DO ji = 1, 6
          cell_index = iidx(jv,jb,ji)
          cell_block = iblk(jv,jb,ji)
          IF (cell_index > 0)                                      &
            & p_vert_out(jv,jk,jb) = p_vert_out(jv,jk,jb) +                                        &
            &      c_int(jv,ji,jb) * p_cell_in(cell_index,jk,cell_block)
        ENDDO

!       Old version of the code that relies on non-zero values in iidx => ptr_patch%verts%cell_idx
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

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!Orinal interpolation routines developed by by Einar Olason. Depreciated.
!------------------------------------------------------------------------
!------------------------------------------------------------------------
  !-------------------------------------------------------------------------
  !
  !> Map vectors from vertices to edges
  !! Based on ideas from rot_vertex_ocean_3d in oce_dyn_icohom/mo_ocean_math_operators.f90
  !!
  !! @par Revision History
  !! Developed by Einar Olason, MPI-M (2013-08-05)
  !!
  !! Uses incorrect mapping coefficients: edge2cell_coeff_cc_t
  !! Does not include OMP parallelization
  !! Replaced by map_verts2edges, Vladimir Lapin, MPI-M (2015-08-05)
  !
  SUBROUTINE map_verts2edges_einar(p_patch_3D, p_vn_dual, edge2cell_coeff_cc_t, vn)

    TYPE(t_patch_3D), TARGET, INTENT(in)      :: p_patch_3D
    TYPE(t_cartesian_coordinates),INTENT(IN) :: p_vn_dual(:,:)
    TYPE(t_cartesian_coordinates),INTENT(in):: edge2cell_coeff_cc_t(:,:,:,:)
    REAL(wp), INTENT(inOUT)           :: vn(:,:)

    ! Local variables

    ! Patch and sub-set
    TYPE(t_subset_range), POINTER :: verts_in_domain
    TYPE(t_patch), POINTER        :: p_patch

    ! Indexing
    INTEGER :: jb, jv, jev
    INTEGER :: ile, ibe
    INTEGER :: il_v1, il_v2,ib_v1, ib_v2
    INTEGER :: i_startidx_v, i_endidx_v

    !-----------------------------------------------------------------------

    p_patch         => p_patch_3D%p_patch_2D(1)
    verts_in_domain => p_patch%verts%in_domain


    ! Copied from oce_dyn_icohom/mo_ocean_math_operators.f90:rot_vertex_ocean_3d
    ! Replaced the coefficient and sign to get normal velocity instead of tangental
    DO jb = verts_in_domain%start_block, verts_in_domain%end_block
      CALL get_index_range(verts_in_domain, jb, i_startidx_v, i_endidx_v)
        DO jv = i_startidx_v, i_endidx_v

          DO jev = 1, p_patch%verts%num_edges(jv,jb)

            ! get line and block indices of edge jev around vertex jv
            ile = p_patch%verts%edge_idx(jv,jb,jev)
            ibe = p_patch%verts%edge_blk(jv,jb,jev)

            IF(p_patch_3D%lsm_e(ile,1,ibe) <= sea_boundary)THEN
              !calculate normal velocity
              il_v1 = p_patch%edges%vertex_idx(ile,ibe,1)
              ib_v1 = p_patch%edges%vertex_blk(ile,ibe,1)
              il_v2 = p_patch%edges%vertex_idx(ile,ibe,2)
              ib_v2 = p_patch%edges%vertex_blk(ile,ibe,2)

              vn(ile,ibe) = &
                &   DOT_PRODUCT(p_vn_dual(il_v1,ib_v1)%x,edge2cell_coeff_cc_t(ile,1,ibe,1)%x) &
                & + DOT_PRODUCT(p_vn_dual(il_v2,ib_v2)%x,edge2cell_coeff_cc_t(ile,1,ibe,2)%x)
            ELSE
              vn(ile,ibe) = 0._wp
            ENDIF
          END DO
      END DO
    END DO

  END SUBROUTINE map_verts2edges_einar



END MODULE mo_icon_to_fem_interpolation
