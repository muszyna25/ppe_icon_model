!>
!!   Contains the implementation of the mathematical operators for the ocean.
!!
!!   Contains the implementation of the mathematical operators for the ocean.
!!
!! @par Revision History
!!  Developed  by Peter Korn and Stephan Lorenz 2010-04
!!  Modified by Stephan Lorenz                  2011-02
!!    correct implementation of ocean boundaries
!!
!! @par To Do
!! Boundary exchange, nblks in presence of halos and dummy edge
!!
!! @par Copyright
!! 2002-2006 by DWD and MPI-M
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
MODULE mo_oce_math_operators
  !-------------------------------------------------------------------------
  USE mo_kind,               ONLY: wp
  USE mo_parallel_config,    ONLY: nproma
  USE mo_run_config,         ONLY: ltimer
  USE mo_math_constants
  USE mo_physical_constants
  USE mo_impl_constants,     ONLY: boundary, sea, sea_boundary !,land, land_boundary, sea, max_char_length, &
  USE mo_model_domain,       ONLY: t_patch
  USE mo_ext_data_types,     ONLY: t_external_data
  USE mo_ocean_nml,          ONLY: n_zlev, iswm_oce
  USE mo_dynamics_config,    ONLY: nold
  USE mo_util_dbg_prnt,      ONLY: dbg_print
  USE mo_exception,          ONLY: finish, message
#ifndef __SX__
  USE mo_timer,              ONLY: timer_start, timer_stop, timer_div, timer_grad
#endif
  USE mo_oce_state,          ONLY: t_hydro_ocean_state, v_base
  USE mo_intp_data_strc,     ONLY: p_int_state
  USE mo_math_utilities,     ONLY: t_cartesian_coordinates, vector_product !, gc2cc
  USE mo_operator_ocean_coeff_3d, ONLY: t_operator_coeff
  USE mo_grid_subset,         ONLY: t_subset_range, get_index_range
  USE mo_sync,                ONLY: SYNC_C, SYNC_E, SYNC_V, sync_patch_array

  IMPLICIT NONE

  PRIVATE

  CHARACTER(LEN=*), PARAMETER :: version = '$Id$'
  CHARACTER(len=12)           :: str_module    = 'oceMathOps  '  ! Output of module for 1 line debug
  INTEGER                     :: idt_src       = 1               ! Level of detail for 1 line debug


  PUBLIC :: grad_fd_norm_oce_3d
  PUBLIC :: div_oce_3d
  PUBLIC :: rot_vertex_ocean_3d

  PUBLIC :: grad_fd_norm_oce_2d, grad_fd_norm_oce_2d_3d

  PUBLIC :: rot_vertex_ocean
  PUBLIC :: rot_vertex_ocean_rbf
  PUBLIC :: height_related_quantities
  PUBLIC :: map_edges2vert_3D

  INTERFACE div_oce_3d
    MODULE PROCEDURE div_oce_3d_mlevels
    MODULE PROCEDURE div_oce_3d_1level
  END INTERFACE

CONTAINS

  !-------------------------------------------------------------------------
  !>
  !! mpi parallelized by LL
  SUBROUTINE map_edges2vert_3d(p_patch, vn, h_e, edge2vert_coeff_cc, p_vn_dual)
    
    TYPE(t_patch), TARGET, INTENT(in)      :: p_patch
    REAL(wp), INTENT(in)           :: vn(:,:,:)
    TYPE(t_cartesian_coordinates),INTENT(in):: edge2vert_coeff_cc(:,:,:,:)
    REAL(wp), INTENT(in)           :: h_e(:,:)
    TYPE(t_cartesian_coordinates)  :: p_vn_dual(nproma,n_zlev,p_patch%nblks_v)


    !Local variables
    !

    !REAL(wp) :: zarea_fraction
    !REAL(wp) :: z_area_scaled
    INTEGER :: slev, elev     ! vertical start and end level
    INTEGER :: jv, jk, jb,jev
    INTEGER :: ile, ibe
    INTEGER :: i_startidx_v, i_endidx_v
    !INTEGER :: icell_idx_1, icell_blk_1
    !INTEGER :: icell_idx_2, icell_blk_2

    INTEGER :: i_v_ctr(nproma,n_zlev,p_patch%nblks_v)
    INTEGER,PARAMETER :: ino_dual_edges = 6
    TYPE(t_subset_range), POINTER :: verts_in_domain
    !-----------------------------------------------------------------------
    verts_in_domain => p_patch%verts%in_domain

    i_v_ctr(:,:,:) = 0
    slev         = 1
    elev         = n_zlev

    DO jb = verts_in_domain%start_block, verts_in_domain%end_block
      CALL get_index_range(verts_in_domain, jb, i_startidx_v, i_endidx_v)
      DO jk = slev, elev
        DO jv = i_startidx_v, i_endidx_v

          p_vn_dual(jv,jk,jb)%x = 0.0_wp
          DO jev = 1, p_patch%verts%num_edges(jv,jb)

            ! get line and block indices of edge jev around vertex jv
            ile = p_patch%verts%edge_idx(jv,jb,jev)
            ibe = p_patch%verts%edge_blk(jv,jb,jev)

              p_vn_dual(jv,jk,jb)%x = p_vn_dual(jv,jk,jb)%x        &
              & +edge2vert_coeff_cc(jv,jk,jb,jev)%x &
              & *vn(ile,jk,ibe)!/(p_patch%verts%dual_area(jv,jb)/(re*re))
          END DO
        END DO ! jv = i_startidx_v, i_endidx_v
      END DO ! jk = slev, elev
    END DO ! jb = verts_in_domain%start_block, verts_in_domain%end_block

    ! sync the result
    CALL sync_patch_array(SYNC_V, p_patch, p_vn_dual(:,:,:)%x(1))
    CALL sync_patch_array(SYNC_V, p_patch, p_vn_dual(:,:,:)%x(2))
    CALL sync_patch_array(SYNC_V, p_patch, p_vn_dual(:,:,:)%x(3))

  END SUBROUTINE map_edges2vert_3d
  !-------------------------------------------------------------------------------------
          !
!   !-------------------------------------------------------------------------
!   !>
!   !!  no-mpi parallelized
!   SUBROUTINE map_edges2vert(p_patch, vn, h_e, p_vn_dual, subset_range)
!     
!     TYPE(t_patch), TARGET, INTENT(in)      :: p_patch
!     REAL(wp), INTENT(in)           :: vn(:,:,:)
!     REAL(wp), INTENT(in)           :: h_e(:,:)
!     TYPE(t_cartesian_coordinates)  :: p_vn_dual(nproma,n_zlev,p_patch%nblks_v)
!     TYPE(t_subset_range), OPTIONAL :: subset_range
! 
!     !Local variables
!     !
!     !REAL(wp) :: z_weight(nproma,n_zlev,p_patch%nblks_v)
!     REAL(wp) :: zarea_fraction
!     REAL(wp) :: z_area_scaled
! 
!     INTEGER :: slev, elev     ! vertical start and end level
!     INTEGER :: jv, jk, jb,jev
!     INTEGER :: ile, ibe
!     INTEGER :: i_startblk_v, i_endblk_v, i_startidx_v, i_endidx_v
!     !INTEGER :: i_startblk_e, i_endblk_e, i_startidx_e, i_endidx_e
!     INTEGER,PARAMETER :: rl_start_v = 2
!     INTEGER,PARAMETER :: rl_end_v   = min_rlvert
! 
!     INTEGER :: icell_idx_1, icell_blk_1
!     INTEGER :: icell_idx_2, icell_blk_2
!     !INTEGER :: il_v1, il_v2,ib_v1, ib_v2
!     INTEGER :: i_v_ctr(nproma,n_zlev,p_patch%nblks_v)
!     TYPE(t_cartesian_coordinates) :: cell1_cc, cell2_cc, vertex_cc
!     INTEGER,PARAMETER :: ino_dual_edges = 6
! 
!     TYPE(t_subset_range), POINTER :: verts_in_domain
!     !-----------------------------------------------------------------------
!     verts_in_domain => p_patch%verts%in_domain
! 
!     i_v_ctr(:,:,:) = 0
!     slev         = 1
!     elev         = n_zlev
! 
!     i_startblk_v = p_patch%verts%start_blk(rl_start_v,1)
!     i_endblk_v   = p_patch%verts%end_blk(rl_end_v,1)
! 
!     DO jb = verts_in_domain%start_block, verts_in_domain%end_block
!       CALL get_index_range(verts_in_domain, jb, i_startidx_v, i_endidx_v)
!       DO jk = slev, elev
!         DO jv = i_startidx_v, i_endidx_v
! 
!           zarea_fraction      = 0.0_wp
!           z_area_scaled       = 0.0_wp
!           !i_bdr_ctr           = 0
!           !z_weight(jv,jk,jb) = 0.0_wp
!           p_vn_dual(jv,jk,jb)%x = 0.0_wp
! 
!           vertex_cc = gc2cc(p_patch%verts%vertex(jv,jb))
!           DO jev = 1, p_patch%verts%num_edges(jv,jb)
! 
!             ! get line and block indices of edge jev around vertex jv
!             ile = p_patch%verts%edge_idx(jv,jb,jev)
!             ibe = p_patch%verts%edge_blk(jv,jb,jev)
!             !Check, if edge is sea or boundary edge and take care of dummy edge
!             ! edge with indices ile, ibe is sea edge
!             IF ( v_base%lsm_oce_e(ile,jk,ibe) == sea) THEN
! 
!               p_vn_dual(jv,jk,jb)%x = p_vn_dual(jv,jk,jb)%x        &
!                 & +p_int_state(1)%edge2vert_coeff_cc(jv,jb,jev)%x &
!                 & *vn(ile,jk,ibe)!*z_thick
! 
!               !z_weight might be an alternative to dual_area and can include
!               !varying height in top layer. Differences have to be explored.
!               !z_weight(jv,jk,jb) = z_weight(jv,jk,jb) &
!               !&+ p_int_state(1)%variable_dual_vol_norm(jv,jb,jev)!*z_thick
! 
!               !increase wet edge ctr
!               i_v_ctr(jv,jk,jb)=i_v_ctr(jv,jk,jb)+1
!             END IF
!           END DO
!           !
!           !divide by hex/pentagon area, if all dual cells are in the ocean interior
!           !divide by apropriate fraction if boundaries are involved
! 
!           IF ( i_v_ctr(jv,jk,jb) == p_patch%verts%num_edges(jv,jb) ) THEN
! 
!             z_area_scaled         = p_patch%verts%dual_area(jv,jb)/(re*re)
!             p_vn_dual(jv,jk,jb)%x = p_vn_dual(jv,jk,jb)%x/z_area_scaled!z_weight(jv,jk,jb)
! 
! 
!           ELSEIF(i_v_ctr(jv,jk,jb)/=0)THEN!boundary edges are involved
! 
!             !Modified area calculation
!             DO jev = 1, p_patch%verts%num_edges(jv,jb)
!               ! get line and block indices of edge jev around vertex jv
!               ile = p_patch%verts%edge_idx(jv,jb,jev)
!               ibe = p_patch%verts%edge_blk(jv,jb,jev)
!               !get neighbor cells
!               icell_idx_1 = p_patch%edges%cell_idx(ile,ibe,1)
!               icell_idx_2 = p_patch%edges%cell_idx(ile,ibe,2)
!               icell_blk_1 = p_patch%edges%cell_blk(ile,ibe,1)
!               icell_blk_2 = p_patch%edges%cell_blk(ile,ibe,2)
!               cell1_cc    = gc2cc(p_patch%cells%center(icell_idx_1,icell_blk_1))
!               cell2_cc    = gc2cc(p_patch%cells%center(icell_idx_2,icell_blk_2))
!               !Check, if edge is sea or boundary edge and take care of dummy edge
!               ! edge with indices ile, ibe is sea edge
!               !Add up for wet dual area.
!               IF ( v_base%lsm_oce_e(ile,jk,ibe) <= sea_boundary ) THEN
!                 zarea_fraction = zarea_fraction  &
!                   & + triangle_area(cell1_cc, vertex_cc, cell2_cc)
!                 ! edge with indices ile, ibe is boundary edge
!               ELSE IF ( v_base%lsm_oce_e(ile,jk,ibe) == boundary ) THEN
!                 zarea_fraction = zarea_fraction  &
!                   & + 0.5_wp*triangle_area(cell1_cc, vertex_cc, cell2_cc)
!               END IF
!             END DO
! 
!             ! no division by zero
!             IF (zarea_fraction /= 0.0_wp) THEN
!               !z_area_scaled   = zarea_fraction
!               z_area_scaled       = p_patch%verts%dual_area(jv,jb)/(re*re)
!               p_vn_dual(jv,jk,jb)%x  = p_vn_dual(jv,jk,jb)%x/z_area_scaled!z_weight(jv,jk,jb)!
!             ENDIF
!           ENDIF
! 
!         END DO
!       END DO
!     END DO
! 
!     IF (PRESENT(subset_range)) THEN
!       IF (.NOT.  subset_range%is_in_domain) THEN
!         CALL sync_patch_array(SYNC_E, p_patch, p_vn_dual(:,:,:)%x(1))
!         CALL sync_patch_array(SYNC_E, p_patch, p_vn_dual(:,:,:)%x(2))
!         CALL sync_patch_array(SYNC_E, p_patch, p_vn_dual(:,:,:)%x(3))
!       ENDIF
!     ENDIF
! 
!   END SUBROUTINE map_edges2vert
!   !-------------------------------------------------------------------------



  !-------------------------------------------------------------------------
  !>
  !!  Computes directional derivative of a cell centered variable in presence of lateral boundaries
  !!  as in the ocean model setting.
  !!
  !!  Computes directional  derivative of a cell centered variable
  !!  with respect to direction normal to triangle edge.
  !! input: lives on centres of triangles
  !! output:  lives on edges (velocity points)
  !!
  !! @par Revision History
  !! Developed  by  Luca Bonaventura, MPI-M (2002-5).
  !! Adapted to new data structure by Peter Korn
  !! and Luca Bonaventura, MPI-M (2005).
  !! Modifications by P. Korn, MPI-M(2007-2)
  !! -Switch fom array arguments to pointers
  !! Modification by Almut Gassmann, MPI-M (2007-04-20)
  !! - abandon grid for the sake of patch
  !!Boundary handling for triangles by P. Korn (2009)
  !!  mpi note: the result is not synced. Should be done in the calling method if required
  SUBROUTINE grad_fd_norm_oce_3d( psi_c, ptr_patch, grad_coeff, grad_norm_psi_e)

    !  patch on which computation is performed
    !
    TYPE(t_patch), TARGET, INTENT(in) :: ptr_patch
    REAL(wp), INTENT(in)              :: grad_coeff(:,:,:)
    !  cell based variable of which normal derivative is computed
    REAL(wp), INTENT(in)              ::  psi_c(:,:,:)       ! dim: (nproma,n_zlev,nblks_c)
    !  edge based variable in which normal derivative is stored
    REAL(wp), INTENT(inout) ::  grad_norm_psi_e(:,:,:)  ! dim: (nproma,n_zlev,nblks_e)

    !
    !Local variables
    !
    INTEGER :: slev, elev     ! vertical start and end level
    INTEGER :: je, jk, jb
    !INTEGER :: rl_start, rl_end
    !INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER :: i_startidx, i_endidx
    INTEGER,  DIMENSION(:,:,:),   POINTER :: iidx, iblk
    TYPE(t_subset_range), POINTER :: edges_in_domain

    !-----------------------------------------------------------------------
    edges_in_domain => ptr_patch%edges%in_domain

    slev = 1
    elev = n_zlev

    iidx => ptr_patch%edges%cell_idx
    iblk => ptr_patch%edges%cell_blk
    !
    !  loop through all patch edges (and blocks)
    !
#ifndef __SX__
    IF (ltimer) CALL timer_start(timer_grad)
#endif
!$OMP PARALLEL
    ! The special treatment of 2D fields is essential for efficiency on the NEC

#ifdef __SX__
    IF (slev > 1) THEN
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,je,jk)
    DO jb = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, jb, i_startidx, i_endidx)
#ifdef _URD2
!CDIR UNROLL=_URD2
#endif
        DO jk = slev, elev
          DO je = i_startidx, i_endidx
            IF ( v_base%lsm_oce_e(je,jk,jb) <= sea_boundary) THEN
              grad_norm_psi_e(je,jk,jb) = grad_coeff(je,jk,jb)* &
              & ( psi_c(iidx(je,jb,2),jk,iblk(je,jb,2))-        &
              & psi_c(iidx(je,jb,1),jk,iblk(je,jb,1)) )
            ENDIF
          ENDDO
        END DO
      END DO
!$OMP END DO
    ELSE
#endif

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,je,jk)
    DO jb = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, jb, i_startidx, i_endidx)
#ifdef _URD
!CDIR UNROLL=_URD
#endif
        DO jk = slev, elev
          DO je = i_startidx, i_endidx
            IF ( v_base%lsm_oce_e(je,jk,jb) <= sea_boundary) THEN
              grad_norm_psi_e(je,jk,jb) = grad_coeff(je,jk,jb)* &
                & ( psi_c(iidx(je,jb,2),jk,iblk(je,jb,2)) -     &
                & psi_c(iidx(je,jb,1),jk,iblk(je,jb,1)) )
            ENDIF 
          ENDDO
        END DO
      END DO
!$OMP END DO
#ifdef __SX__
    ENDIF
#endif


!$OMP END PARALLEL
#ifndef __SX__
    IF (ltimer) CALL timer_stop(timer_grad)
#endif
  END SUBROUTINE grad_fd_norm_oce_3d
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! Computes discrete divergence of a vector field in presence of lateral boundaries as in ocean setting.
  !!
  !! Computes discrete divergence of a vector field
  !! given by its components in the directions normal to triangle edges.
  !! The midpoint rule is used for quadrature.
  !! input:  lives on edges (velocity points)
  !! output: lives on centers of triangles
  !!
  !! @par Revision History
  !! Developed  by  Luca Bonaventura, MPI-M (2002-5).
  !! Changes according to programming guide by Thomas Heinze, DWD (2006-08-18).
  !! Modification by Thomas Heinze, DWD (2006-09-11):
  !! - loop only over the inner cells of a patch, not any more over halo cells
  !! Modifications by P. Korn, MPI-M(2007-2)
  !! - Switch fom array arguments to pointers
  !! Modification by Almut Gassmann, MPI-M (2007-04-20)
  !! - abandon grid for the sake of patch
  !! Modification by Guenther Zaengl, DWD (2009-03-17)
  !! - vector optimization
  !! Modification by Peter Korn, MPI-M    (2009)
  !! - Boundary treatment for the ocean
  !! Modification by Stephan Lorenz, MPI-M (2010-08-05)
  !! - New boundary definition with inner and boundary points on land/sea
  !!
  !!  mpi parallelized LL (no sync required)
  SUBROUTINE div_oce_3d_mlevels( vec_e, ptr_patch, div_coeff, div_vec_c, opt_slev, opt_elev, &
    & subset_range)
    !
    !
    !  patch on which computation is performed
    !
    TYPE(t_patch), TARGET, INTENT(in) :: ptr_patch
    !
    ! edge based variable of which divergence
    ! is computed
    !
    REAL(wp), INTENT(in)          :: vec_e(:,:,:) ! dim: (nproma,n_zlev,nblks_e)
    REAL(wp), INTENT(in)          :: div_coeff(:,:,:,:)
    REAL(wp), INTENT(inout)       :: div_vec_c(:,:,:) ! dim: (nproma,n_zlev,nblks_c)
    INTEGER, INTENT(in), OPTIONAL :: opt_slev       ! optional vertical start level
    INTEGER, INTENT(in), OPTIONAL :: opt_elev       ! optional vertical end level
    TYPE(t_subset_range), TARGET, INTENT(in), OPTIONAL :: subset_range

    INTEGER :: slev, elev     ! vertical start and end level
    INTEGER :: jc, jk, jb
    INTEGER ::i_startidx, i_endidx
    INTEGER,  DIMENSION(:,:,:),   POINTER :: iidx, iblk
    TYPE(t_subset_range), POINTER :: all_cells
    !-----------------------------------------------------------------------
    IF (PRESENT(subset_range)) THEN
      all_cells => subset_range
    ELSE
      all_cells => ptr_patch%cells%all
    ENDIF

    IF ( PRESENT(opt_slev) ) THEN
      slev = opt_slev
    ELSE
      slev = 1
    END IF
    IF ( PRESENT(opt_elev) ) THEN
      elev = opt_elev
    ELSE
      elev = n_zlev
    END IF

#ifndef __SX__
    IF (ltimer) CALL timer_start(timer_div)
#endif
!$OMP PARALLEL

    iidx => ptr_patch%cells%edge_idx
    iblk => ptr_patch%cells%edge_blk

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk)
    DO jb = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, jb, i_startidx, i_endidx)
#ifdef __SX__
!CDIR UNROLL=6
#endif
      DO jk = slev, elev
        DO jc = i_startidx, i_endidx
          ! compute the discrete divergence for cell jc by finite volume
          ! approximation (see Bonaventura and Ringler MWR 2005);
          ! multiplication of the normal vector component vec_e at the edges
          ! by the appropriate cell based edge_orientation is required to
          ! obtain the correct value for the application of Gauss theorem
          ! (which requires the scalar product of the vector field with the
          ! OUTWARD pointing unit vector with respect to cell jc; since the
          ! positive direction for the vector components is not necessarily
          ! the outward pointing one with respect to cell jc, a correction
          ! coefficient (equal to +-1) is necessary, given by
          ! ptr_patch%grid%cells%edge_orientation)
          !
          ! Distinghuish: case of a land cell (put div to zero), and
          ! cases where one of the edges are boundary or land
          ! (put corresponding velocity to zero).
          ! sea, sea_boundary, boundary (edges only), land_boundary, land =
          !  -2,      -1,         0,                  1,             2
          !This information is stored inside the divergence coefficients.
          div_vec_c(jc,jk,jb) =  &
            & vec_e(iidx(jc,jb,1),jk,iblk(jc,jb,1)) * div_coeff(jc,jk,jb,1) + &
            & vec_e(iidx(jc,jb,2),jk,iblk(jc,jb,2)) * div_coeff(jc,jk,jb,2) + &
            & vec_e(iidx(jc,jb,3),jk,iblk(jc,jb,3)) * div_coeff(jc,jk,jb,3)
        END DO
      END DO
    END DO
!$OMP END DO

!$OMP END PARALLEL
#ifndef __SX__
    IF (ltimer) CALL timer_stop(timer_div)
#endif
  END SUBROUTINE div_oce_3d_mlevels

  !-------------------------------------------------------------------------
  !>
  !! Computes discrete divergence of a vector field in presence of lateral boundaries as in ocean setting.
  !!
  !! Computes discrete divergence of a vector field
  !! given by its components in the directions normal to triangle edges.
  !! The midpoint rule is used for quadrature.
  !! input:  lives on edges (velocity points)
  !! output: lives on centers of triangles
  !!
  !! @par Revision History
  !! Developed  by  Luca Bonaventura, MPI-M (2002-5).
  !! Changes according to programming guide by Thomas Heinze, DWD (2006-08-18).
  !! Modification by Thomas Heinze, DWD (2006-09-11):
  !! - loop only over the inner cells of a patch, not any more over halo cells
  !! Modifications by P. Korn, MPI-M(2007-2)
  !! - Switch fom array arguments to pointers
  !! Modification by Almut Gassmann, MPI-M (2007-04-20)
  !! - abandon grid for the sake of patch
  !! Modification by Guenther Zaengl, DWD (2009-03-17)
  !! - vector optimization
  !! Modification by Peter Korn, MPI-M    (2009)
  !! - Boundary treatment for the ocean
  !! Modification by Stephan Lorenz, MPI-M (2010-08-05)
  !! - New boundary definition with inner and boundary points on land/sea
  !!
  !!  mpi parallelized LL (no sync required)
  SUBROUTINE div_oce_3d_1level( vec_e, ptr_patch, div_coeff, div_vec_c,  &
    & level, subset_range)
    !
    !
    !  patch on which computation is performed
    !
    TYPE(t_patch), TARGET, INTENT(in) :: ptr_patch
    !
    ! edge based variable of which divergence
    ! is computed
    !
    REAL(wp), INTENT(inout)       :: vec_e(:,:) ! dim: (nproma,n_zlev,nblks_e)
    REAL(wp), INTENT(in)          :: div_coeff(:,:,:,:)
    REAL(wp), INTENT(inout)       :: div_vec_c(:,:) ! dim: (nproma,n_zlev,nblks_c)
    INTEGER,  INTENT(in)          :: level
    TYPE(t_subset_range), TARGET, INTENT(in), OPTIONAL :: subset_range

    INTEGER :: jc, jb
    INTEGER :: i_startidx, i_endidx
    !INTEGER :: nlen, npromz_c, nblks_c
    INTEGER,  DIMENSION(:,:,:),   POINTER :: iidx, iblk
    TYPE(t_subset_range), POINTER :: all_cells
    !-----------------------------------------------------------------------
    IF (PRESENT(subset_range)) THEN
      all_cells => subset_range
    ELSE
      all_cells => ptr_patch%cells%all
    ENDIF

#ifndef __SX__
    IF (ltimer) CALL timer_start(timer_div)
#endif
!$OMP PARALLEL

    iidx => ptr_patch%cells%edge_idx
    iblk => ptr_patch%cells%edge_blk

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc)
    DO jb = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, jb, i_startidx, i_endidx)
#ifdef __SX__
!CDIR UNROLL=6
#endif
        DO jc = i_startidx, i_endidx

          ! compute the discrete divergence for cell jc by finite volume
          ! approximation (see Bonaventura and Ringler MWR 2005);
          ! multiplication of the normal vector component vec_e at the edges
          ! by the appropriate cell based edge_orientation is required to
          ! obtain the correct value for the application of Gauss theorem
          ! (which requires the scalar product of the vector field with the
          ! OUTWARD pointing unit vector with respect to cell jc; since the
          ! positive direction for the vector components is not necessarily
          ! the outward pointing one with respect to cell jc, a correction
          ! coefficient (equal to +-1) is necessary, given by
          ! ptr_patch%grid%cells%edge_orientation)
          !
          ! Distinghuish: case of a land cell (put div to zero), and
          ! cases where one of the edges are boundary or land
          ! (put corresponding velocity to zero).
          ! sea, sea_boundary, boundary (edges only), land_boundary, land =
          !  -2,      -1,         0,                  1,             2
          !This information is stored inside the divergence coefficients.
          div_vec_c(jc,jb) =  &
            & vec_e(iidx(jc,jb,1),iblk(jc,jb,1)) * div_coeff(jc,level,jb,1) + &
            & vec_e(iidx(jc,jb,2),iblk(jc,jb,2)) * div_coeff(jc,level,jb,2) + &
            & vec_e(iidx(jc,jb,3),iblk(jc,jb,3)) * div_coeff(jc,level,jb,3)
        END DO
    END DO
!$OMP END DO

!$OMP END PARALLEL
#ifndef __SX__
    IF (ltimer) CALL timer_stop(timer_div)
#endif
  END SUBROUTINE div_oce_3d_1level
!   !-------------------------------------------------------------------------
  !
  !>
  !!  Computes directional derivative of a cell centered variable in presence of lateral boundaries
  !!  as in the ocean model setting.
  !!
  !!  Computes directional  derivative of a cell centered variable
  !!  with respect to direction normal to triangle edge.
  !! input: lives on centres of triangles
  !! output:  lives on edges (velocity points)
  !!
  !! @par Revision History
  !! Developed  by  Luca Bonaventura, MPI-M (2002-5).
  !! Adapted to new data structure by Peter Korn
  !! and Luca Bonaventura, MPI-M (2005).
  !! Modifications by P. Korn, MPI-M(2007-2)
  !! -Switch fom array arguments to pointers
  !! Modification by Almut Gassmann, MPI-M (2007-04-20)
  !! - abandon grid for the sake of patch
  !! Boundary handling for triangles by P. Korn (2009)
  !!  mpi note: the result is not synced. Should be done in the calling method if required
  !!
  SUBROUTINE grad_fd_norm_oce_2d_3d( psi_c, ptr_patch, grad_coeff, grad_norm_psi_e)
    !
    TYPE(t_patch), TARGET, INTENT(in) :: ptr_patch
    REAL(wp), INTENT(in)    :: psi_c(:,:)             ! dim: (nproma,nblks_c)
    REAL(wp), INTENT(in)    :: grad_coeff(:,:,:)
    REAL(wp), INTENT(inout) ::  grad_norm_psi_e(:,:)  ! dim: (nproma,nblks_e)

    !!
    !!local variables
    INTEGER :: slev, elev     ! vertical start and end level
    INTEGER :: je, jb
    INTEGER :: i_startidx, i_endidx
    INTEGER,  DIMENSION(:,:,:),   POINTER :: iidx, iblk
    TYPE(t_subset_range), POINTER :: edges_in_domain
    !-----------------------------------------------------------------------
    edges_in_domain => ptr_patch%edges%in_domain
    slev = 1
    elev = 1

    iidx => ptr_patch%edges%cell_idx
    iblk => ptr_patch%edges%cell_blk
        !
    !  loop through all patch edges (and blocks)
    !
#ifndef __SX__
    IF (ltimer) CALL timer_start(timer_grad)
#endif
!$OMP PARALLEL
    ! The special treatment of 2D fields is essential for efficiency on the NEC
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,je)
    DO jb = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, jb, i_startidx, i_endidx)
#ifdef _URD
!CDIR UNROLL=_URD
#endif

      DO je = i_startidx, i_endidx
        ! compute the normal derivative
        ! by the finite difference approximation
        ! (see Bonaventura and Ringler MWR 2005)
        grad_norm_psi_e(je,jb) =  &
          & (psi_c(iidx(je,jb,2),iblk(je,jb,2))-psi_c(iidx(je,jb,1),iblk(je,jb,1)))& 
          & * grad_coeff(je,slev,jb)
      END DO
    END DO
!$OMP END DO
!$OMP END PARALLEL
#ifndef __SX__
    IF (ltimer) CALL timer_stop(timer_grad)
#endif
  END SUBROUTINE grad_fd_norm_oce_2d_3d
  !-------------------------------------------------------------------------
  !
  !>
  !!  Computes directional derivative of a cell centered variable in presence of lateral boundaries
  !!  as in the ocean model setting.
  !!
  !!  Computes directional  derivative of a cell centered variable
  !!  with respect to direction normal to triangle edge.
  !! input: lives on centres of triangles
  !! output:  lives on edges (velocity points)
  !!
  !! @par Revision History
  !! Developed  by  Luca Bonaventura, MPI-M (2002-5).
  !! Adapted to new data structure by Peter Korn
  !! and Luca Bonaventura, MPI-M (2005).
  !! Modifications by P. Korn, MPI-M(2007-2)
  !! -Switch fom array arguments to pointers
  !! Modification by Almut Gassmann, MPI-M (2007-04-20)
  !! - abandon grid for the sake of patch
  !! Boundary handling for triangles by P. Korn (2009)
  !!   mpi note: this will no compute the halo values.
  !!             if necessary they have to be synced by the calling method
  !!
  SUBROUTINE grad_fd_norm_oce_2d( psi_c, ptr_patch, grad_norm_psi_e)
    !
    !
    !  patch on which computation is performed
    !
    TYPE(t_patch), TARGET, INTENT(in) :: ptr_patch
    !
    !  cell based variable of which normal derivative is computed
    !
    REAL(wp), INTENT(in)    :: psi_c(:,:)             ! dim: (nproma,nblks_c)
    !REAL(wp), INTENT(in)    :: grad_coeff(:,:,:)
    !
    !  edge based variable in which normal derivative is stored
    !
    !REAL(wp), INTENT(out) ::  &
    REAL(wp), INTENT(inout) ::  grad_norm_psi_e(:,:)  ! dim: (nproma,nblks_e)

    INTEGER :: slev, elev     ! vertical start and end level
    INTEGER :: je, jb
    !INTEGER :: rl_start, rl_end
    INTEGER :: i_startidx, i_endidx
    INTEGER :: nblks_e, npromz_e
    INTEGER,  DIMENSION(:,:,:),   POINTER :: iidx, iblk
    !
    TYPE(t_subset_range), POINTER :: edges_in_domain
    !-----------------------------------------------------------------------
    edges_in_domain => ptr_patch%edges%in_domain

    slev = 1
    elev = 1

    iidx => ptr_patch%edges%cell_idx
    iblk => ptr_patch%edges%cell_blk

    !
    !  loop through all patch edges (and blocks)
    !
#ifndef __SX__
    IF (ltimer) CALL timer_start(timer_grad)
#endif
!$OMP PARALLEL
    ! The special treatment of 2D fields is essential for efficiency on the NEC

    SELECT CASE (ptr_patch%cell_type)

    CASE (3) ! (cell_type == 3)
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,je)
      DO jb = edges_in_domain%start_block, edges_in_domain%end_block
        CALL get_index_range(edges_in_domain, jb, i_startidx, i_endidx)
#ifdef _URD
!CDIR UNROLL=_URD
#endif
        DO je = i_startidx, i_endidx
          ! compute the normal derivative
          ! by the finite difference approximation
          ! (see Bonaventura and Ringler MWR 2005)
          !
          !IF ( v_base%lsm_oce_e(je,1,jb) <= sea_boundary ) THEN
          grad_norm_psi_e(je,jb) =  &
            & ( psi_c(iidx(je,jb,2),iblk(je,jb,2)) - psi_c(iidx(je,jb,1),iblk(je,jb,1)) ) &
            & * ptr_patch%edges%inv_dual_edge_length(je,jb) &
            & *v_base%wet_e(je,1,jb)
          !ELSE
          !  grad_norm_psi_e(je,jb) =  0.0_wp
          !ENDIF
        END DO
      END DO
!$OMP END DO

    CASE (6) ! (cell_type == 6)
      ! no grid refinement in hexagonal model
      nblks_e   = ptr_patch%nblks_int_e
      npromz_e  = ptr_patch%npromz_int_e

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,je)
      DO jb = edges_in_domain%start_block, edges_in_domain%end_block
        CALL get_index_range(edges_in_domain, jb, i_startidx, i_endidx)
        DO je = i_startidx, i_endidx
          !
          ! compute the normal derivative
          ! by the finite difference approximation
          ! (see Bonaventura and Ringler MWR 2005)
          !
          grad_norm_psi_e(je,jb) =  &
            & ( psi_c(iidx(je,jb,2),iblk(je,jb,2)) - psi_c(iidx(je,jb,1),iblk(je,jb,1)) ) &
            & * ptr_patch%edges%inv_dual_edge_length(je,jb) &
            & * ptr_patch%edges%system_orientation(je,jb)
        ENDDO
      END DO
!$OMP END DO
    END SELECT
!$OMP END PARALLEL
#ifndef __SX__
    IF (ltimer) CALL timer_stop(timer_grad)
#endif
  END SUBROUTINE grad_fd_norm_oce_2d
  !-------------------------------------------------------------------------


    !
!   !-------------------------------------------------------------------------
  !! Computes the discrete rotation at vertices in presence of boundaries as in the ocean setting.
  !! Computes in presence of boundaries the discrete rotation at vertices
  !! of triangle cells (centers of dual grid cells) from a vector field
  !! given by its components in the directions normal to triangle edges and
  !! takes the presence of boundaries into account.
  !!
  !! This sbr calculates the vorticity for mimetic discretization. A second one for the RBF-branch
  !! can be found below. The two sbr use the same approach to calculate the curl, but they differ
  !! in the calculation of the tangential velocity, which is only need at lateral boundaries. Mimetic
  !! does the tangential velocity calculate from velocity vector at vertices (p_vn_dual), while RBF uses
  !! a specific routine for that purpose.
  !!   mpi note: the results is not synced. should be done by the calling method if necessary
  SUBROUTINE rot_vertex_ocean_3d( p_patch, vn, p_vn_dual, p_op_coeff, rot_vec_v)
    !>
    !!
    TYPE(t_patch), TARGET, INTENT(in)         :: p_patch
    REAL(wp), INTENT(in)                      :: vn(:,:,:)
    TYPE(t_cartesian_coordinates), INTENT(in) :: p_vn_dual(nproma,n_zlev,p_patch%nblks_v)
    TYPE(t_operator_coeff),TARGET, INTENT(in) :: p_op_coeff
    REAL(wp), INTENT(inout)                   :: rot_vec_v(:,:,:)

    !Local variables
    !
    REAL(wp) :: z_vort_int
    REAL(wp) :: z_vort_boundary(nproma,n_zlev,p_patch%nblks_v)
    INTEGER :: slev, elev
    INTEGER :: jv, jk, jb, jev
    INTEGER :: ile, ibe
    INTEGER :: il_v1, il_v2,ib_v1, ib_v2
    INTEGER :: i_startidx_v, i_endidx_v

    REAL(wp) :: z_vt(nproma,n_zlev,p_patch%nblks_e)
    !INTEGER,PARAMETER :: ino_dual_edges = 6
    INTEGER, POINTER :: ibnd_edge_idx(:,:,:,:), ibnd_edge_blk(:,:,:,:), i_edge_idx(:,:,:,:)
    REAL(wp), POINTER :: z_orientation(:,:,:,:)

    TYPE(t_subset_range), POINTER :: verts_in_domain

    !-----------------------------------------------------------------------
    verts_in_domain => p_patch%verts%in_domain
    slev         = 1
    elev         = n_zlev


    z_vort_boundary  = 0.0_wp
    rot_vec_v(:,:,:) = 0.0_wp
    z_vt(:,:,:)      = 0.0_wp

    !set pointer that carry edge information
    ibnd_edge_idx    => p_op_coeff%bnd_edge_idx
    ibnd_edge_blk    => p_op_coeff%bnd_edge_blk
    i_edge_idx       => p_op_coeff%edge_idx
    z_orientation    => p_op_coeff%orientation


    !In this loop tangential velocity is calulated at boundaries
    DO jb = verts_in_domain%start_block, verts_in_domain%end_block
      CALL get_index_range(verts_in_domain, jb, i_startidx_v, i_endidx_v)
      DO jk = slev, elev
        !!$OMP PARALLEL DO SCHEDULE(runtime) DEFAULT(PRIVATE)  &
        !!$OMP   SHARED(u_vec_e,v_vec_e,ptr_patch,rot_vec_v,jb) FIRSTPRIVATE(jk)
        DO jv = i_startidx_v, i_endidx_v

          DO jev = 1, p_patch%verts%num_edges(jv,jb)

            ! get line and block indices of edge jev around vertex jv
            ile = p_patch%verts%edge_idx(jv,jb,jev)
            ibe = p_patch%verts%edge_blk(jv,jb,jev)

            IF ( v_base%lsm_oce_e(ile,jk,ibe) == boundary ) THEN
              !calculate tangential velocity
              il_v1 = p_patch%edges%vertex_idx(ile,ibe,1)
              ib_v1 = p_patch%edges%vertex_blk(ile,ibe,1)
              il_v2 = p_patch%edges%vertex_idx(ile,ibe,2)
              ib_v2 = p_patch%edges%vertex_blk(ile,ibe,2)

              z_vt(ile,jk,ibe)= &
                & - DOT_PRODUCT(p_vn_dual(il_v1,jk,ib_v1)%x,&
                & p_op_coeff%edge2vert_coeff_cc_t(ile,jk,ibe,1)%x)&
                !& p_int_state(1)%edge2vert_coeff_cc_t(ile,ibe,2)%x)&
                & + DOT_PRODUCT(p_vn_dual(il_v2,jk,ib_v2)%x,&
                & p_op_coeff%edge2vert_coeff_cc_t(ile,jk,ibe,2)%x)
                !& p_int_state(1)%edge2vert_coeff_cc_t(ile,ibe,1)%x)
            ENDIF
          END DO
        END DO
        !!$OMP END PARALLEL DO
      END DO
    END DO



    !In this loop vorticity at vertices is calculated
    DO jb = verts_in_domain%start_block, verts_in_domain%end_block
      CALL get_index_range(verts_in_domain, jb, i_startidx_v, i_endidx_v)
      DO jk = slev, elev
        !!$OMP PARALLEL DO SCHEDULE(runtime) DEFAULT(PRIVATE)  &
        !!$OMP   SHARED(u_vec_e,v_vec_e,ptr_patch,rot_vec_v,jb) FIRSTPRIVATE(jk)
        DO jv = i_startidx_v, i_endidx_v

          z_vort_int = 0.0_wp

          DO jev = 1, p_patch%verts%num_edges(jv,jb)

            ! get line and block indices of edge jev around vertex jv
            ile = p_patch%verts%edge_idx(jv,jb,jev)
            ibe = p_patch%verts%edge_blk(jv,jb,jev)
            !add contribution of normal velocity at edge (ile,ibe) to rotation
            !IF ( v_base%lsm_oce_e(ile,jk,ibe) == sea) THEN
            ! sea, sea_boundary, boundary (edges only), land_boundary, land =
            !  -2,      -1,         0,                  1,             2
            !Distinction between sea-lean-boundary is taken into account by coeffcients.
            !It is assumed here that vn is already zero at boundary edges.

            z_vort_int = z_vort_int + vn(ile,jk,ibe)*p_op_coeff%rot_coeff(jv,jk,jb,jev)
          END DO

          !Finalize vorticity calculation by closing the dual loop along boundary edges
          !IF(i_v_bnd_edge_ctr(jv,jk,jb)==2)THEN
          IF(p_op_coeff%bnd_edges_per_vertex(jv,jk,jb)==2)THEN

            z_vort_boundary(jv,jk,jb)=&
              & z_vt(ibnd_edge_idx(jv,jk,jb,1),jk,ibnd_edge_blk(jv,jk,jb,1))&
              & *p_op_coeff%rot_coeff(jv,jk,jb,i_edge_idx(jv,jk,jb,1))&
              & +&
              & z_vt(ibnd_edge_idx(jv,jk,jb,2),jk,ibnd_edge_blk(jv,jk,jb,2))&
              & *p_op_coeff%rot_coeff(jv,jk,jb,i_edge_idx(jv,jk,jb,2))

            !ELSEIF(i_v_bnd_edge_ctr(jv,jk,jb)==4)THEN
          ELSEIF(p_op_coeff%bnd_edges_per_vertex(jv,jk,jb)==4)THEN
            !In case of 4 boundary edges within a dual loop, we have 2 land triangles
            !around the vertex. these two land triangles have one vertex in common and are
            !seperated by two wet triangles

            z_vort_boundary(jv,jk,jb)=&
              & z_vt(ibnd_edge_idx(jv,jk,jb,1),jk,ibnd_edge_blk(jv,jk,jb,1))&
              & *p_op_coeff%rot_coeff(jv,jk,jb,i_edge_idx(jv,jk,jb,1))&
              & +&
              & z_vt(ibnd_edge_idx(jv,jk,jb,2),jk,ibnd_edge_blk(jv,jk,jb,2))&
              & *p_op_coeff%rot_coeff(jv,jk,jb,i_edge_idx(jv,jk,jb,2))&
              & +&
              & z_vt(ibnd_edge_idx(jv,jk,jb,3),jk,ibnd_edge_blk(jv,jk,jb,3))&
              & *p_op_coeff%rot_coeff(jv,jk,jb,i_edge_idx(jv,jk,jb,3))&
              & +&
              & z_vt(ibnd_edge_idx(jv,jk,jb,4),jk,ibnd_edge_blk(jv,jk,jb,4))&
              & *p_op_coeff%rot_coeff(jv,jk,jb,i_edge_idx(jv,jk,jb,4))
          ENDIF

          !Final vorticity calculation
!TODO ram
    !       rot_vec_v(jv,jk,jb) = (z_vort_int + z_vort_boundary(jv,jk,jb)) / &
    !         & p_patch%verts%dual_area(jv,jb)
          rot_vec_v(jv,jk,jb) = z_vort_int + z_vort_boundary(jv,jk,jb)

        END DO
        !!$OMP END PARALLEL DO
      END DO
    END DO

    ! DO jb = i_startblk_v, i_endblk_v
    !   CALL get_indices_v(p_patch, jb, i_startblk_v, i_endblk_v, &
    !                      i_startidx_v, i_endidx_v, rl_start_v, rl_end_v)
    !   DO jk = slev, elev
    !     DO jv = i_startidx_v, i_endidx_v
    ! IF(rot_vec_v(jv,jk,jb)/=0.0_wp)THEN
    ! write(123456,*)'rot 3D:',jk,jv,jb,rot_vec_v(jv,jk,jb),&
    ! &z_vort_boundary(jv,jk,jb), p_op_coeff%bnd_edges_per_vertex(jv,jk,jb)
    ! ENDIF
    !     END DO
    !   END DO
    ! END DO
  END SUBROUTINE rot_vertex_ocean_3d
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !
  !! Computes the discrete rotation at vertices in presence of boundaries as in the ocean setting.
  !! Computes in presence of boundaries the discrete rotation at vertices
  !! of triangle cells (centers of dual grid cells) from a vector field
  !! given by its components in the directions normal to triangle edges and
  !! takes the presence of boundaries into account.
  !!
  !! This sbr calculates the vorticity for mimetic discretization. A second one for the RBF-branch
  !! can be found below. The two sbr use the same approach to calculate the curl, but they differ
  !! in the calculation of the tangential velocity, which is only need at lateral boundaries. Mimetic
  !! does the tangential velocity calculate from velocity vector at vertices (p_vn_dual), while RBF uses
  !! a specific routine for that purpose.
  !!  mpi note: the result is not synced. Should be done in the calling method if required
  !!
  SUBROUTINE rot_vertex_ocean( p_patch, vn, p_vn_dual, rot_vec_v)
    !>
    !!
    TYPE(t_patch), TARGET, INTENT(in) :: p_patch
    REAL(wp), INTENT(inout)           :: vn(:,:,:)
    TYPE(t_cartesian_coordinates)  :: p_vn_dual(nproma,n_zlev,p_patch%nblks_v)
    REAL(wp), INTENT(inout)        :: rot_vec_v(:,:,:)

    !Local variables
    !
    REAL(wp) :: z_vort_tmp, z_vort_tmp_boundary(nproma,n_zlev,p_patch%nblks_v)
    !REAL(wp) :: z_weight(nproma,n_zlev,p_patch%nblks_v)
    REAL(wp) :: zarea_fraction
    !REAL(wp) :: z_area_scaled

    INTEGER :: slev, elev     ! vertical start and end level
    INTEGER :: jv, jk, jb, jev
    INTEGER :: ile, ibe
    INTEGER :: i_startidx_v, i_endidx_v
    !INTEGER :: i_bdr_ctr
    !INTEGER :: icell_idx_1, icell_blk_1
    !INTEGER :: icell_idx_2, icell_blk_2
    INTEGER :: il_v1, il_v2,ib_v1, ib_v2
    INTEGER :: i_v_ctr(nproma,n_zlev,p_patch%nblks_v)
    INTEGER :: i_v_bnd_edge_ctr(nproma,n_zlev,p_patch%nblks_v)
    INTEGER :: ibnd_edge_idx(4), ibnd_edge_blk(4)  !maximal 4 boundary edges in a dual loop.
    INTEGER :: i_edge_idx(4)
    REAL(wp) :: z_orientation(4)
    REAL(wp) :: z_vt(nproma,n_zlev,p_patch%nblks_e)
    !TYPE(t_cartesian_coordinates) :: vertex_cc!cell1_cc, cell2_cc
    !INTEGER,PARAMETER :: ino_dual_edges = 6

    TYPE(t_subset_range), POINTER :: verts_in_domain
    !-----------------------------------------------------------------------
    verts_in_domain => p_patch%verts%in_domain

    slev         = 1
    elev         = n_zlev

    ! #slo# due to nag -nan compiler-option
    i_v_ctr(:,:,:)          = 0
    i_v_bnd_edge_ctr(:,:,:) = 0
    ibnd_edge_idx(1:4)      = 0
    ibnd_edge_blk(1:4)      = 0
    rot_vec_v(:,:,:)        = 0.0_wp
    z_vt(:,:,:)             = 0.0_wp
    z_orientation(1:4)      = 0.0_wp
    z_vort_tmp_boundary     = 0.0_wp

    !In this loop vorticity at vertices is calculated
    DO jb = verts_in_domain%start_block, verts_in_domain%end_block
      CALL get_index_range(verts_in_domain, jb, i_startidx_v, i_endidx_v)
      DO jk = slev, elev
        !!$OMP PARALLEL DO SCHEDULE(runtime) DEFAULT(PRIVATE)  &
        !!$OMP   SHARED(u_vec_e,v_vec_e,ptr_patch,rot_vec_v,jb) FIRSTPRIVATE(jk)
        DO jv = i_startidx_v, i_endidx_v

          z_vort_tmp          = 0.0_wp
          zarea_fraction      = 0.0_wp

          !vertex_cc = gc2cc(p_patch%verts%vertex(jv,jb))
          DO jev = 1, p_patch%verts%num_edges(jv,jb)

            ! get line and block indices of edge jev around vertex jv
            ile = p_patch%verts%edge_idx(jv,jb,jev)
            ibe = p_patch%verts%edge_blk(jv,jb,jev)
            !Check, if edge is sea or boundary edge and take care of dummy edge
            ! edge with indices ile, ibe is sea edge

            IF ( v_base%lsm_oce_e(ile,jk,ibe) == sea) THEN
              !Distinguish the following cases
              ! edge ie_k is
              !a) ocean edge: compute as usual,
              !b) land edge: do not consider it
              !c) boundary edge take:
              !  no-slip boundary condition:  normal and tangential velocity at boundary are zero
              ! sea, sea_boundary, boundary (edges only), land_boundary, land =
              !  -2,      -1,         0,                  1,             2
              !add contribution of normal velocity at edge (ile,ibe) to rotation
              z_vort_tmp = z_vort_tmp + vn(ile,jk,ibe)                   &
                & * p_patch%edges%dual_edge_length(ile,ibe)  &
                & * p_patch%verts%edge_orientation(jv,jb,jev)


              !increase wet edge ctr
              i_v_ctr(jv,jk,jb)=i_v_ctr(jv,jk,jb)+1

              ! edge with indices ile, ibe is boundary edge
            ELSE IF ( v_base%lsm_oce_e(ile,jk,ibe) == boundary ) THEN

              !calculate tangential velocity
              il_v1 = p_patch%edges%vertex_idx(ile,ibe,1)
              ib_v1 = p_patch%edges%vertex_blk(ile,ibe,1)
              il_v2 = p_patch%edges%vertex_idx(ile,ibe,2)
              ib_v2 = p_patch%edges%vertex_blk(ile,ibe,2)

              z_vt(ile,jk,ibe)= &
                & - DOT_PRODUCT(p_vn_dual(il_v1,jk,ib_v1)%x,&
                & p_int_state(1)%edge2vert_coeff_cc_t(ile,ibe,2)%x)&
                & + DOT_PRODUCT(p_vn_dual(il_v2,jk,ib_v2)%x,&
                & p_int_state(1)%edge2vert_coeff_cc_t(ile,ibe,1)%x)

              !increase boundary edge counter
              i_v_bnd_edge_ctr(jv,jk,jb)=i_v_bnd_edge_ctr(jv,jk,jb)+1

              !Store actual boundary edge indices
              IF(i_v_bnd_edge_ctr(jv,jk,jb)==1)THEN
                ibnd_edge_idx(1) = ile
                ibnd_edge_blk(1) = ibe
                z_orientation(1) = p_patch%verts%edge_orientation(jv,jb,jev)
                i_edge_idx(1)    = jev
              ELSEIF(i_v_bnd_edge_ctr(jv,jk,jb)==2)THEN
                ibnd_edge_idx(2) = ile
                ibnd_edge_blk(2) = ibe
                z_orientation(2) = p_patch%verts%edge_orientation(jv,jb,jev)
                i_edge_idx(2)    = jev
              ELSEIF(i_v_bnd_edge_ctr(jv,jk,jb)==3)THEN
                ibnd_edge_idx(3) = ile
                ibnd_edge_blk(3) = ibe
                z_orientation(3) = p_patch%verts%edge_orientation(jv,jb,jev)
                i_edge_idx(3)    = jev
              ELSEIF(i_v_bnd_edge_ctr(jv,jk,jb)==4)THEN
                ibnd_edge_idx(4) = ile
                ibnd_edge_blk(4) = ibe
                z_orientation(4) = p_patch%verts%edge_orientation(jv,jb,jev)
                i_edge_idx(4)    = jev
              ELSE
                !maximal 4 boundary edges per dual loop are allowed: somethings wrong withe the grid
                ! write(*,*)'grid error',jv,jk,jb,i_v_bnd_edge_ctr(jv,jk,jb)
                CALL message (TRIM('sbr nonlinear Coriolis'), &
                  & 'more than 4 boundary edges per dual loop: something is wrong with the grid')
                CALL finish ('TRIM(sbr nonlinear Coriolis)','Grid-boundary error !!')
              ENDIF
            END IF
          END DO

          IF ( i_v_ctr(jv,jk,jb) == p_patch%verts%num_edges(jv,jb) ) THEN

            rot_vec_v(jv,jk,jb) = z_vort_tmp/p_patch%verts%dual_area(jv,jb)! (re*re*z_weight(jv,jk,jb))!

            !Finalize vorticity calculation by closing the dual loop along boundary edges
          ELSEIF(i_v_bnd_edge_ctr(jv,jk,jb)==2)THEN!(i_v_ctr(jv,jk,jb)==2)THEN

            z_vort_tmp_boundary(jv,jk,jb) =&
            !& p_patch%edges%system_orientation(ibnd_edge_idx(1),ibnd_edge_blk(1))*&
              & z_orientation(1)*&
              & z_vt(ibnd_edge_idx(1),jk,ibnd_edge_blk(1)) &
              & *p_patch%edges%primal_edge_length(ibnd_edge_idx(1),ibnd_edge_blk(1))&
              & +&
            !& p_patch%edges%system_orientation(ibnd_edge_idx(2),ibnd_edge_blk(2))*&
              & z_orientation(2)*&
              & z_vt(ibnd_edge_idx(2),jk,ibnd_edge_blk(2)) &
              & *p_patch%edges%primal_edge_length(ibnd_edge_idx(2),ibnd_edge_blk(2))

            rot_vec_v(jv,jk,jb) = (z_vort_tmp+0.5_wp*z_vort_tmp_boundary(jv,jk,jb))&
              & /p_patch%verts%dual_area(jv,jb)!z_area_scaled

          ELSEIF(i_v_bnd_edge_ctr(jv,jk,jb)==4)THEN!(i_v_ctr(jv,jk,jb)==4)THEN

            !In case of 4 boundary edges within a dual loop, we have 2 land triangles
            !around the vertex. these two land triangles have one vertex in common and are
            !seperated by two wet triangles.
            z_vort_tmp_boundary(jv,jk,jb) =&
              & z_orientation(1)*&
              & z_vt(ibnd_edge_idx(1),jk,ibnd_edge_blk(1)) &
              & *p_patch%edges%primal_edge_length(ibnd_edge_idx(1),ibnd_edge_blk(1))&
              & +&
              & z_orientation(2)*&
              & z_vt(ibnd_edge_idx(2),jk,ibnd_edge_blk(2)) &
              & *p_patch%edges%primal_edge_length(ibnd_edge_idx(2),ibnd_edge_blk(2))&
              & +&
              & z_orientation(3)*&
              & z_vt(ibnd_edge_idx(3),jk,ibnd_edge_blk(3)) &
              & *p_patch%edges%primal_edge_length(ibnd_edge_idx(3),ibnd_edge_blk(3))&
              & +&
              & z_orientation(4)*&
              & z_vt(ibnd_edge_idx(4),jk,ibnd_edge_blk(4)) &
              & *p_patch%edges%primal_edge_length(ibnd_edge_idx(4),ibnd_edge_blk(4))

            rot_vec_v(jv,jk,jb) = (z_vort_tmp+0.5_wp*z_vort_tmp_boundary(jv,jk,jb))&
              & /p_patch%verts%dual_area(jv,jb)!z_area_scaled

          ENDIF
          !ENDIF
        END DO
        !!$OMP END PARALLEL DO
      END DO
    END DO

    ! DO jb = i_startblk_v, i_endblk_v
    !   CALL get_indices_v(p_patch, jb, i_startblk_v, i_endblk_v, &
    !                      i_startidx_v, i_endidx_v, rl_start_v, rl_end_v)
    !   DO jk = slev, elev
    !     DO jv = i_startidx_v, i_endidx_v
    ! IF(rot_vec_v(jv,jk,jb)/=0.0_wp)THEN
    ! write(12345,*)'rot 2D:',jk,jv,jb,rot_vec_v(jv,jk,jb),&
    ! &0.5_wp*z_vort_tmp_boundary(jv,jk,jb)!, i_v_bnd_edge_ctr(jv,jk,jb)!,i_v_ctr(jv,jk,jb)
    ! ENDIF
    !     END DO
    !   END DO
    ! END DO
  END SUBROUTINE rot_vertex_ocean
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !! Computes the discrete rotation at vertices in presence of boundaries as in the ocean setting.
  !! Computes in presence of boundaries the discrete rotation at vertices
  !! of triangle cells (centers of dual grid cells) from a vector field
  !! given by its components in the directions normal to triangle edges and
  !! takes the presence of boundaries into account.

  !! This sbr calculates the vorticity for RBF discretization. A second one for the Mimetic-branch
  !! can be found above. The two sbr use the same approach to calculate the curl, but they differ
  !! in the calculation of the tangential velocity, which is only need at lateral boundaries. Mimetic
  !! does the tangential velocity calculate from velocity vector at vertices (p_vn_dual), while RBF uses
  !! a specific routine for that purpose. The RBF-sbr has to be called before rot_vertex_ocean_rbf,
  !! and the tangential velocity is passed as an argument.
  !!  mpi note: the result is not synced. Should be done in the calling method if required
  !!
  SUBROUTINE rot_vertex_ocean_rbf( p_patch, vn, vt, rot_vec_v)
    !>
    !!
    TYPE(t_patch), TARGET, INTENT(in) :: p_patch
    REAL(wp), INTENT(in)           :: vn(:,:,:)
    REAL(wp), INTENT(in)           :: vt(:,:,:)
    REAL(wp), INTENT(inout)        :: rot_vec_v(:,:,:)

    !Local variables
    !
    REAL(wp) :: z_vort_tmp, z_vort_tmp_boundary
    !REAL(wp) :: z_weight(nproma,n_zlev,p_patch%nblks_v)
    REAL(wp) :: zarea_fraction

    INTEGER :: slev, elev     ! vertical start and end level
    INTEGER :: jv, jk, jb, jev
    INTEGER :: ile, ibe
    !INTEGER :: rl_start_e, rl_end_e
    INTEGER :: i_startidx_v, i_endidx_v

    !INTEGER :: i_bdr_ctr
    !INTEGER :: icell_idx_1, icell_blk_1
    !INTEGER :: icell_idx_2, icell_blk_2
    !INTEGER :: il_v1, il_v2,ib_v1, ib_v2
    INTEGER :: i_v_ctr(nproma,n_zlev,p_patch%nblks_v)
    INTEGER :: i_v_bnd_edge_ctr(nproma,n_zlev,p_patch%nblks_v)
    INTEGER :: ibnd_edge_idx(4), ibnd_edge_blk(4)  !maximal 4 boundary edges in a dual loop.
    INTEGER :: i_edge_idx(4)
    REAL(wp) :: z_orientation(4)
    TYPE(t_cartesian_coordinates) :: vertex_cc! cell1_cc, cell2_cc
    !INTEGER,PARAMETER :: ino_dual_edges = 6

    !INTEGER,PARAMETER :: rl_start_v = 2
    !INTEGER,PARAMETER :: rl_end_v   = min_rlvert
    TYPE(t_subset_range), POINTER :: verts_in_domain
    !-----------------------------------------------------------------------
    verts_in_domain => p_patch%verts%in_domain

    slev         = 1
    elev         = n_zlev

    ! #slo# due to nag -nan compiler-option
    i_v_ctr(:,:,:)          = 0
    i_v_bnd_edge_ctr(:,:,:) = 0
    ibnd_edge_idx(1:4)      = 0
    ibnd_edge_blk(1:4)      = 0
    rot_vec_v(:,:,:)        = 0.0_wp
    z_orientation(1:4)      = 0.0_wp

    !In this loop vorticity at vertices is calculated
    DO jb = verts_in_domain%start_block, verts_in_domain%end_block
      CALL get_index_range(verts_in_domain, jb, i_startidx_v, i_endidx_v)
      DO jk = slev, elev

        DO jv = i_startidx_v, i_endidx_v

          z_vort_tmp          = 0.0_wp
          zarea_fraction      = 0.0_wp
          !i_bdr_ctr           = 0
          !z_weight(jv,jk,jb) = 0.0_wp

          vertex_cc = p_patch%verts%cartesian(jv,jb)
          DO jev = 1, p_patch%verts%num_edges(jv,jb)

            ! get line and block indices of edge jev around vertex jv
            ile = p_patch%verts%edge_idx(jv,jb,jev)
            ibe = p_patch%verts%edge_blk(jv,jb,jev)
            !Check, if edge is sea or boundary edge and take care of dummy edge
            ! edge with indices ile, ibe is sea edge
            IF ( v_base%lsm_oce_e(ile,jk,ibe) == sea) THEN
              !Distinguish the following cases
              ! edge ie_k is
              !a) ocean edge: compute as usual,
              !b) land edge: do not consider it
              !c) boundary edge take:
              !  no-slip boundary condition:  normal and tangential velocity at boundary are zero
              ! sea, sea_boundary, boundary (edges only), land_boundary, land =
              !  -2,      -1,         0,                  1,             2
              !add contribution of normal velocity at edge (ile,ibe) to rotation
              z_vort_tmp = z_vort_tmp + vn(ile,jk,ibe)                   &
                & * p_patch%edges%dual_edge_length(ile,ibe)  &
                & * p_patch%verts%edge_orientation(jv,jb,jev)

              !z_weight might be an alternative to dual_area and can include
              !varying height in top layer. Differences have to be explored.
              !z_weight(jv,jk,jb) = z_weight(jv,jk,jb) &
              !&+ p_int_state(1)%variable_dual_vol_norm(jv,jb,jev)!*z_thick


              !increase wet edge ctr
              i_v_ctr(jv,jk,jb)=i_v_ctr(jv,jk,jb)+1

              ! edge with indices ile, ibe is boundary edge
            ELSE IF ( v_base%lsm_oce_e(ile,jk,ibe) == boundary ) THEN

              !increase boundary edge counter
              i_v_bnd_edge_ctr(jv,jk,jb) = i_v_bnd_edge_ctr(jv,jk,jb)+1

              !Store actual boundary edge indices
              IF(i_v_bnd_edge_ctr(jv,jk,jb)==1)THEN
                ibnd_edge_idx(1) = ile
                ibnd_edge_blk(1) = ibe
                z_orientation(1) = p_patch%verts%edge_orientation(jv,jb,jev)
                i_edge_idx(1)    = jev
              ELSEIF(i_v_bnd_edge_ctr(jv,jk,jb)==2)THEN
                ibnd_edge_idx(2) = ile
                ibnd_edge_blk(2) = ibe
                z_orientation(2) = p_patch%verts%edge_orientation(jv,jb,jev)
                i_edge_idx(2)    = jev
              ELSEIF(i_v_bnd_edge_ctr(jv,jk,jb)==3)THEN
                ibnd_edge_idx(3) = ile
                ibnd_edge_blk(3) = ibe
                z_orientation(3) = p_patch%verts%edge_orientation(jv,jb,jev)
                i_edge_idx(3)    = jev
              ELSEIF(i_v_bnd_edge_ctr(jv,jk,jb)==4)THEN
                ibnd_edge_idx(4) = ile
                ibnd_edge_blk(4) = ibe
                z_orientation(4) = p_patch%verts%edge_orientation(jv,jb,jev)
                i_edge_idx(4)    = jev
              ELSE
                !only 2 boundary edges per dual loop are allowed: somethings wrong withe the grid
                ! write(*,*)'grid error',jv,jk,jb,i_v_bnd_edge_ctr(jv,jk,jb)
                CALL message (TRIM('sbr nonlinear Coriolis'), &
                  & 'more than 2 boundary edges per dual loop: something is wrong with the grid')
                CALL finish ('TRIM(sbr nonlinear Coriolis)','Grid-boundary error !!')
              ENDIF
            END IF
          END DO

          !write(*,*)'no: sea edges+bnd edges',i_v_ctr(jv,jk,jb),i_v_bnd_edge_ctr(jv,jk,jb)
          !
          !divide by hex/pentagon area, if all dual cells are in the ocean interior
          !divide by apropriate fraction if boundaries are involved

          IF ( i_v_ctr(jv,jk,jb) == p_patch%verts%num_edges(jv,jb) ) THEN

            rot_vec_v(jv,jk,jb) = z_vort_tmp /p_patch%verts%dual_area(jv,jb)! (re*re*z_weight(jv,jk,jb))!


          ELSEIF(i_v_ctr(jv,jk,jb)/=0)THEN!boundary edges are involved

            !Modified area calculation
            ! !         DO jev = 1, p_patch%verts%num_edges(jv,jb)
            ! !           ! get line and block indices of edge jev around vertex jv
            ! !           ile = p_patch%verts%edge_idx(jv,jb,jev)
            ! !           ibe = p_patch%verts%edge_blk(jv,jb,jev)
            ! !           !get neighbor cells
            ! !           icell_idx_1 = p_patch%edges%cell_idx(ile,ibe,1)
            ! !           icell_idx_2 = p_patch%edges%cell_idx(ile,ibe,2)
            ! !           icell_blk_1 = p_patch%edges%cell_blk(ile,ibe,1)
            ! !           icell_blk_2 = p_patch%edges%cell_blk(ile,ibe,2)
            ! !           cell1_cc = gc2cc(p_patch%cells%center(icell_idx_1,icell_blk_1))
            ! !           cell2_cc = gc2cc(p_patch%cells%center(icell_idx_2,icell_blk_2))
            ! !           !Check, if edge is sea or boundary edge and take care of dummy edge
            ! !           ! edge with indices ile, ibe is sea edge
            ! !           !Add up for wet dual area.
            ! !           IF ( v_base%lsm_oce_e(ile,jk,ibe) <= sea_boundary ) THEN
            ! !            zarea_fraction = zarea_fraction  &
            ! !               &     + triangle_area(cell1_cc, vertex_cc, cell2_cc)
            ! !           ! edge with indices ile, ibe is boundary edge
            ! !           ELSE IF ( v_base%lsm_oce_e(ile,jk,ibe) == boundary ) THEN
            ! !            zarea_fraction = zarea_fraction  &
            ! !               &  + 0.5_wp*triangle_area(cell1_cc, vertex_cc, cell2_cc)
            ! !           END IF
            ! !         END DO
            ! !         z_area_scaled   = zarea_fraction*re*re
            !z_area_scaled       = p_patch%verts%dual_area(jv,jb)/(re*re)

            !Finalize vorticity calculation by closing the dual loop along boundary edges
            IF(i_v_ctr(jv,jk,jb)==2)THEN

              z_vort_tmp_boundary =&
              !& p_patch%edges%system_orientation(ibnd_edge_idx(1),ibnd_edge_blk(1))*&
                & z_orientation(1)*&
                & vt(ibnd_edge_idx(1),jk,ibnd_edge_blk(1)) &
                & *p_patch%edges%primal_edge_length(ibnd_edge_idx(1),ibnd_edge_blk(1))&
                & +&
              !& p_patch%edges%system_orientation(ibnd_edge_idx(2),ibnd_edge_blk(2))*&
                & z_orientation(2)*&
                & vt(ibnd_edge_idx(2),jk,ibnd_edge_blk(2)) &
                & *p_patch%edges%primal_edge_length(ibnd_edge_idx(2),ibnd_edge_blk(2))

              rot_vec_v(jv,jk,jb) = (z_vort_tmp+0.5_wp*z_vort_tmp_boundary)&
                & /p_patch%verts%dual_area(jv,jb)!z_area_scaled

              !  write(*,*)'vorticity:',jv,jk,jb,z_vort_tmp,z_vort_tmp/zarea_fraction,&
              !  &z_vort_tmp_boundary,0.5_wp*z_vort_tmp_boundary/zarea_fraction, vort_v(jv,jk,jb)
            ELSEIF(i_v_ctr(jv,jk,jb)==4)THEN

              !In case of 4 boundary edges within a dual loop, we have 2 land triangles
              !around the vertex. these two land triangles have one vertex in common and are
              !seperated by two wet triangles.
              z_vort_tmp_boundary =&
                & z_orientation(1)*&
                & vt(ibnd_edge_idx(1),jk,ibnd_edge_blk(1)) &
                & *p_patch%edges%primal_edge_length(ibnd_edge_idx(1),ibnd_edge_blk(1))&
                & +&
                & z_orientation(2)*&
                & vt(ibnd_edge_idx(2),jk,ibnd_edge_blk(2)) &
                & *p_patch%edges%primal_edge_length(ibnd_edge_idx(2),ibnd_edge_blk(2))&
                & +&
                & z_orientation(3)*&
                & vt(ibnd_edge_idx(3),jk,ibnd_edge_blk(3)) &
                & *p_patch%edges%primal_edge_length(ibnd_edge_idx(3),ibnd_edge_blk(3))&
                & +&
                & z_orientation(4)*&
                & vt(ibnd_edge_idx(4),jk,ibnd_edge_blk(4)) &
                & *p_patch%edges%primal_edge_length(ibnd_edge_idx(4),ibnd_edge_blk(4))

              rot_vec_v(jv,jk,jb) = (z_vort_tmp+0.5_wp*z_vort_tmp_boundary)&
                & /p_patch%verts%dual_area(jv,jb)!z_area_scaled
            ENDIF
          ENDIF
        END DO

      END DO
    END DO
  END SUBROUTINE rot_vertex_ocean_rbf
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  ELEMENTAL FUNCTION triangle_area (x0, x1, x2) result(area)

    TYPE(t_cartesian_coordinates), INTENT(in) :: x0, x1, x2

    REAL(wp) :: area
    REAL(wp) :: z_s12, z_s23, z_s31, z_ca1, z_ca2, z_ca3, z_a1, z_a2, z_a3

    TYPE(t_cartesian_coordinates) :: u12, u23, u31

    ! This variant to calculate the area of a spherical triangle
    ! is more precise.

    !  Compute cross products Uij = Vi x Vj.
    u12 = vector_product (x0, x1)
    u23 = vector_product (x1, x2)
    u31 = vector_product (x2, x0)

    !  Normalize Uij to unit vectors.
    z_s12 = DOT_PRODUCT ( u12%x(1:3), u12%x(1:3) )
    z_s23 = DOT_PRODUCT ( u23%x(1:3), u23%x(1:3) )
    z_s31 = DOT_PRODUCT ( u31%x(1:3), u31%x(1:3) )

    !  Test for a degenerate triangle associated with collinear vertices.
    IF ( z_s12 == 0.0_wp .OR. z_s23 == 0.0_wp  .OR. z_s31 == 0.0_wp ) THEN
      area = 0.0_wp
      RETURN
    END IF

    z_s12 = SQRT(z_s12)
    z_s23 = SQRT(z_s23)
    z_s31 = SQRT(z_s31)

    u12%x(1:3) = u12%x(1:3)/z_s12
    u23%x(1:3) = u23%x(1:3)/z_s23
    u31%x(1:3) = u31%x(1:3)/z_s31

    !  Compute interior angles Ai as the dihedral angles between planes:
    !  CA1 = cos(A1) = -<U12,U31>
    !  CA2 = cos(A2) = -<U23,U12>
    !  CA3 = cos(A3) = -<U31,U23>
    z_ca1 = -u12%x(1)*u31%x(1)-u12%x(2)*u31%x(2)-u12%x(3)*u31%x(3)
    z_ca2 = -u23%x(1)*u12%x(1)-u23%x(2)*u12%x(2)-u23%x(3)*u12%x(3)
    z_ca3 = -u31%x(1)*u23%x(1)-u31%x(2)*u23%x(2)-u31%x(3)*u23%x(3)

    IF (z_ca1 < -1.0_wp) z_ca1 = -1.0_wp
    IF (z_ca1 >  1.0_wp) z_ca1 =  1.0_wp
    IF (z_ca2 < -1.0_wp) z_ca2 = -1.0_wp
    IF (z_ca2 >  1.0_wp) z_ca2 =  1.0_wp
    IF (z_ca3 < -1.0_wp) z_ca3 = -1.0_wp
    IF (z_ca3 >  1.0_wp) z_ca3 =  1.0_wp

    z_a1 = ACOS(z_ca1)
    z_a2 = ACOS(z_ca2)
    z_a3 = ACOS(z_ca3)

    !  Compute areas = z_a1 + z_a2 + z_a3 - pi.
    area = z_a1+z_a2+z_a3-pi

    IF ( area < 0.0_wp ) area = 0.0_wp

  END FUNCTION triangle_area
  !-------------------------------------------------------------------------

  !---------------------------------------------------------------------------------
  !>
  !!
  !!  Calculation of total fluid thickness at cell centers and surface elevation at
  !!  cell edges from prognostic surface height at cell centers. We use height at
  !!  old timelevel "n"
  !!
  !! @par Revision History
  !! Developed  by  Peter Korn, MPI-M (2010).
  !!  mpi parallelized LL
  !!
  SUBROUTINE height_related_quantities( p_patch, p_os, p_ext_data)
    !
    ! Patch on which computation is performed
    TYPE(t_patch), TARGET, INTENT(in) :: p_patch
    !
    ! Type containing ocean state
    TYPE(t_hydro_ocean_state), TARGET :: p_os
    !
    ! Type containing external data
    TYPE(t_external_data), TARGET, INTENT(in) :: p_ext_data

    !  local variables
    INTEGER            :: i_startidx_c, i_endidx_c
    INTEGER            :: i_startidx_e, i_endidx_e
    INTEGER            :: jc, jb, je
    INTEGER            :: il_c1, ib_c1, il_c2, ib_c2

    REAL(wp)           :: z_dist_e_c1, z_dist_e_c2
    INTEGER            :: averaging       = 1
    INTEGER, PARAMETER :: distance_weight = 1
    INTEGER, PARAMETER :: upwind          = 2

    TYPE(t_subset_range), POINTER :: all_cells, edges_in_domain

    !TYPE(t_cartesian_coordinates)    :: cv_c1_e0, cv_c2_e0
    !TYPE(t_cartesian_coordinates)    :: cc_c1, cc_c2, cc_e0
    !0!CHARACTER(len=max_char_length), PARAMETER :: &
    !0!  & routine = ('mo_oce_AB_timestepping_mimetic:height_related_quantities')

    !-------------------------------------------------------------------------------
    !CALL message (TRIM(routine), 'start')

    ! sync before run
    CALL sync_patch_array(sync_c, p_patch, p_os%p_prog(nold(1))%h)

    ! #ifndef __SX__
    ! IF (ltimer) CALL timer_start(timer_height)
    ! #endif
    all_cells => p_patch%cells%all
    edges_in_domain => p_patch%edges%in_domain

    !Step 1: calculate cell-located variables for 2D and 3D case
    !For 3D and for SWE thick_c contains thickness of fluid column
    IF ( iswm_oce == 1 ) THEN  !  SWM

      DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
#ifdef __SX__
!CDIR UNROLL=6
#endif
        !calculate for each fluid colum the total depth, i.e.
        !from bottom boundary to surface height, i.e. using individual bathymetry for SWM
        DO jc = i_startidx_c, i_endidx_c
          IF(v_base%lsm_oce_c(jc,1,jb) <= sea_boundary ) THEN
            p_os%p_diag%thick_c(jc,jb) = p_os%p_prog(nold(1))%h(jc,jb)&
              &  + v_base%zlev_i(v_base%dolic_c(jc,jb)+1)!-p_ext_data%oce%bathymetry_c(jc,jb)
          ELSE
            p_os%p_diag%thick_c(jc,jb) = 0.0_wp
          ENDIF
        END DO
      END DO

    ELSEIF(iswm_oce /= 1 )THEN

      DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
#ifdef __SX__
!CDIR UNROLL=6
#endif
        !calculate for each fluid colum the total depth, i.e.
        !from bottom boundary to surface height
        !the bottom boundary is zlev_i(dolic+1) since zlev_i(1)=0 (air-sea-boundary at h=0
        DO jc = i_startidx_c, i_endidx_c
          IF ( v_base%lsm_oce_c(jc,1,jb) <= sea_boundary ) THEN
            p_os%p_diag%thick_c(jc,jb) = p_os%p_prog(nold(1))%h(jc,jb)&
              & + v_base%zlev_i(v_base%dolic_c(jc,jb)+1)
          ELSE
            p_os%p_diag%thick_c(jc,jb) = 0.0_wp
          ENDIF
        END DO
      END DO
    ENDIF!IF ( iswm_oce == 1 )


    !Step 2: calculate edge-located variables for 2D and 3D case from respective cell variables
    !For SWE and for 3D: thick_e = thickness of fluid column at edges
    !         h_e     = surface elevation at edges, without depth of first layer
    IF ( iswm_oce == 1 ) THEN  !  SWM
      SELECT CASE(averaging)

      CASE(distance_weight)

      DO jb = edges_in_domain%start_block, edges_in_domain%end_block
        CALL get_index_range(edges_in_domain, jb, i_startidx_e, i_endidx_e)
#ifdef __SX__
!CDIR UNROLL=6
#endif
          DO je = i_startidx_e, i_endidx_e

            il_c1 = p_patch%edges%cell_idx(je,jb,1)
            ib_c1 = p_patch%edges%cell_blk(je,jb,1)
            il_c2 = p_patch%edges%cell_idx(je,jb,2)
            ib_c2 = p_patch%edges%cell_blk(je,jb,2)

            z_dist_e_c1 = 0.5_wp!p_int_state(1)%dist_cell2edge(je,jb,1) !0.5_wp
            z_dist_e_c2 = 0.5_wp!p_int_state(1)%dist_cell2edge(je,jb,2) !0.5_wp

            !z_dist_e_c1=p_patch%edges%edge_cell_length(je,jb,1)
            !z_dist_e_c2=p_patch%edges%edge_cell_length(je,jb,2)

            IF ( v_base%lsm_oce_e(je,1,jb) <= sea_boundary ) THEN
              p_os%p_diag%thick_e(je,jb) = ( z_dist_e_c1*p_os%p_diag%thick_c(il_c1,ib_c1)&
                & +   z_dist_e_c2*p_os%p_diag%thick_c(il_c2,ib_c2) )&
                & /(z_dist_e_c1+z_dist_e_c2)

              p_os%p_diag%h_e(je,jb) = ( z_dist_e_c1*p_os%p_prog(nold(1))%h(il_c1,ib_c1)&
                & +   z_dist_e_c2*p_os%p_prog(nold(1))%h(il_c2,ib_c2) )&
                & /(z_dist_e_c1+z_dist_e_c2)

              !write(*,*)'height_e',je,jb, p_os%p_diag%h_e(je,jb), p_os%p_prog(nold(1))%h(il_c1,ib_c1),p_os%p_prog(nold(1))%h(il_c2,ib_c2)
            ELSE
              p_os%p_diag%h_e(je,jb) = 0.0_wp
            ENDIF
          END DO
        END DO
        CALL sync_patch_array(SYNC_E, p_patch, p_os%p_diag%thick_e)
        CALL sync_patch_array(SYNC_E, p_patch, p_os%p_diag%h_e)

      CASE(upwind)

        DO jb = edges_in_domain%start_block, edges_in_domain%end_block
          CALL get_index_range(edges_in_domain, jb, i_startidx_e, i_endidx_e)
#ifdef __SX__
!CDIR UNROLL=6
#endif
          DO je = i_startidx_e, i_endidx_e

            il_c1 = p_patch%edges%cell_idx(je,jb,1)
            ib_c1 = p_patch%edges%cell_blk(je,jb,1)
            il_c2 = p_patch%edges%cell_idx(je,jb,2)
            ib_c2 = p_patch%edges%cell_blk(je,jb,2)


            IF( v_base%lsm_oce_e(je,1,jb) <= sea_boundary ) THEN

              IF(p_os%p_prog(nold(1))%vn(je,1,jb)>0.0_wp)THEN
                p_os%p_diag%thick_e(je,jb) = p_os%p_diag%thick_c(il_c1,ib_c1)
                p_os%p_diag%h_e(je,jb)     = p_os%p_prog(nold(1))%h(il_c1,ib_c1)
              ELSEIF(p_os%p_prog(nold(1))%vn(je,1,jb)<=0.0_wp)THEN
                p_os%p_diag%thick_e(je,jb) = p_os%p_diag%thick_c(il_c2,ib_c2)
                p_os%p_diag%h_e(je,jb)     = p_os%p_prog(nold(1))%h(il_c2,ib_c2)
              ENDIF
              !write(*,*)'height_e',je,jb, p_os%p_diag%h_e(je,jb), p_os%p_prog(nold(1))%h(il_c1,ib_c1),p_os%p_prog(nold(1))%h(il_c2,ib_c2)
            ELSE
              p_os%p_diag%h_e(je,jb) = 0.0_wp
            ENDIF
          END DO
        END DO
        CALL sync_patch_array(SYNC_E, p_patch, p_os%p_diag%thick_e)
        CALL sync_patch_array(SYNC_E, p_patch, p_os%p_diag%h_e)


      END SELECT



    ELSEIF(iswm_oce /= 1 )THEN


      SELECT CASE(averaging)

      CASE(distance_weight)

        DO jb = edges_in_domain%start_block, edges_in_domain%end_block
          CALL get_index_range(edges_in_domain, jb, i_startidx_e, i_endidx_e)
#ifdef __SX__
!CDIR UNROLL=6
#endif
          DO je = i_startidx_e, i_endidx_e

            il_c1 = p_patch%edges%cell_idx(je,jb,1)
            ib_c1 = p_patch%edges%cell_blk(je,jb,1)
            il_c2 = p_patch%edges%cell_idx(je,jb,2)
            ib_c2 = p_patch%edges%cell_blk(je,jb,2)

            z_dist_e_c1 = 0.5_wp!p_int_state(1)%dist_cell2edge(je,jb,1)
            z_dist_e_c2 = 0.5_wp!p_int_state(1)%dist_cell2edge(je,jb,2)
            !z_dist_e_c1 = p_patch%edges%edge_cell_length(je,jb,1)
            !z_dist_e_c2 = p_patch%edges%edge_cell_length(je,jb,2)

            IF( v_base%lsm_oce_e(je,1,jb) <= sea_boundary ) THEN
!TODO ram
              p_os%p_diag%thick_e(je,jb) = ( z_dist_e_c1*p_os%p_diag%thick_c(il_c1,ib_c1)&
                & +   z_dist_e_c2*p_os%p_diag%thick_c(il_c2,ib_c2) )&
                & /(z_dist_e_c1+z_dist_e_c2)

              !  P.K.: Actually h_e is just surface elevation at edges without depth of first layer.
!               !It might make sense to include depth of first layer.
!                p_os%p_diag%h_e(je,jb) = ( z_dist_e_c1*p_os%p_prog(nold(1))%h(il_c1,ib_c1)&
!                  & +   z_dist_e_c2*p_os%p_prog(nold(1))%h(il_c2,ib_c2) )&
!                  & /(z_dist_e_c1+z_dist_e_c2)

            p_os%p_diag%thick_e(je,jb)=&
            &min(p_os%p_diag%thick_c(il_c1,ib_c1),p_os%p_diag%thick_c(il_c2,ib_c2))

            p_os%p_diag%h_e(je,jb)=&
            &min(p_os%p_prog(nold(1))%h(il_c1,ib_c1),p_os%p_prog(nold(1))%h(il_c2,ib_c2))


              !write(*,*)'height_e',je,jb, p_os%p_diag%h_e(je,jb), p_os%p_prog(nold(1))%h(il_c1,ib_c1),p_os%p_prog(nold(1))%h(il_c2,ib_c2)
            ELSE
              p_os%p_diag%h_e(je,jb) = 0.0_wp
            ENDIF
          END DO
        END DO
        CALL sync_patch_array(SYNC_E, p_patch, p_os%p_diag%thick_e)
        CALL sync_patch_array(SYNC_E, p_patch, p_os%p_diag%h_e)

      CASE(upwind)

        DO jb = edges_in_domain%start_block, edges_in_domain%end_block
          CALL get_index_range(edges_in_domain, jb, i_startidx_e, i_endidx_e)
#ifdef __SX__
!CDIR UNROLL=6
#endif
          DO je = i_startidx_e, i_endidx_e

            il_c1 = p_patch%edges%cell_idx(je,jb,1)
            ib_c1 = p_patch%edges%cell_blk(je,jb,1)
            il_c2 = p_patch%edges%cell_idx(je,jb,2)
            ib_c2 = p_patch%edges%cell_blk(je,jb,2)


            IF( v_base%lsm_oce_e(je,1,jb) <= sea_boundary ) THEN

              IF(p_os%p_prog(nold(1))%vn(je,1,jb)>0.0_wp)THEN
                p_os%p_diag%thick_e(je,jb) = p_os%p_diag%thick_c(il_c1,ib_c1)
              ELSEIF(p_os%p_prog(nold(1))%vn(je,1,jb)<=0.0_wp)THEN

                !  P.K.: Actually h_e is just surface elevation at edges without depth of first layer.
                !It might make sense to include depth of first layer.
                p_os%p_diag%h_e(je,jb) = p_os%p_prog(nold(1))%h(il_c2,ib_c2)
              ENDIF
              !write(*,*)'height_e',je,jb, p_os%p_diag%h_e(je,jb), p_os%p_prog(nold(1))%h(il_c1,ib_c1),p_os%p_prog(nold(1))%h(il_c2,ib_c2)
            ELSE
              p_os%p_diag%h_e(je,jb) = 0.0_wp
            ENDIF
          END DO
        END DO
        CALL sync_patch_array(SYNC_E, p_patch, p_os%p_diag%thick_e)
        CALL sync_patch_array(SYNC_E, p_patch, p_os%p_diag%h_e)


      END SELECT
    ENDIF
    ! write(*,*)'max/min thick_c:',maxval(p_os%p_diag%thick_c),minval(p_os%p_diag%thick_c)
    ! write(*,*)'max/min h_c:',maxval(p_os%p_prog(nold(1))%h),minval(p_os%p_prog(nold(1))%h)
    ! write(*,*)'max/min bath_c:',maxval(p_ext_data%oce%bathymetry_c),  &
    !   &        minval(p_ext_data%oce%bathymetry_c)
    ! write(*,*)'max/min thick_e:',maxval(p_os%p_diag%thick_e),minval(p_os%p_diag%thick_e)
    ! write(*,*)'max/min h_e:',maxval(p_os%p_diag%h_e),minval(p_os%p_diag%h_e)
    ! !CALL message (TRIM(routine), 'end')

    ! #ifndef __SX__
    ! IF (ltimer) CALL timer_stop(timer_height)
    ! #endif

    !sync results, already done
!     CALL sync_patch_array(sync_c, p_patch, p_os%p_prog(nold(1))%h)
!     CALL sync_patch_array(sync_e, p_patch, p_os%p_diag%h_e)
!     CALL sync_patch_array(sync_c, p_patch, p_os%p_diag%thick_c)
!     CALL sync_patch_array(sync_e, p_patch, p_os%p_diag%thick_e)

    !---------Debug Diagnostics-------------------------------------------
    idt_src=2  ! output print level (1-5, fix)
    CALL dbg_print('heightRelQuant: h_e'            ,p_os%p_diag%h_e        ,str_module,idt_src)
    idt_src=3  ! output print level (1-5, fix)
    CALL dbg_print('heightRelQuant: h_c'            ,p_os%p_prog(nold(1))%h ,str_module,idt_src)
    CALL dbg_print('heightRelQuant: thick_c'        ,p_os%p_diag%thick_c    ,str_module,idt_src)
    CALL dbg_print('heightRelQuant: thick_e'        ,p_os%p_diag%thick_e    ,str_module,idt_src)
    !---------------------------------------------------------------------


  END SUBROUTINE height_related_quantities
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  ! ! !
  ! ! !! Computes the discrete rotation at vertices in presence of boundaries as in the ocean setting.
  ! ! !! Computes in presence of boundaries the discrete rotation at vertices
  ! ! !! of triangle cells (centers of dual grid cells) from a vector field
  ! ! !! given by its components in the directions normal to triangle edges and
  ! ! !! takes the presence of boundaries into account.
  ! ! !!
  ! ! !! This sbr calculates the vorticity for mimetic discretization. A second one for the RBF-branch
  ! ! !! can be found below. The two sbr use the same approach to calculate the curl, but they differ
  ! ! !! in the calculation of the tangential velocity, which is only need at lateral boundaries. Mimetic
  ! ! !! does the tangential velocity calculate from velocity vector at vertices (p_vn_dual), while RBF uses
  ! ! !! a specific routine for that purpose.
  ! ! !!
  ! ! SUBROUTINE rot_vertex_ocean_3D_backup( p_patch, vn, p_vn_dual, rot_coeff, rot_vec_v)
  ! ! !>
  ! ! !!
  ! !   TYPE(t_patch), INTENT(IN)                 :: p_patch
  ! !    REAL(wp), INTENT(INOUT)                      :: vn(:,:,:)
  ! !   TYPE(t_cartesian_coordinates), INTENT(IN) :: p_vn_dual(nproma,n_zlev,p_patch%nblks_v)
  ! !   REAL(wp), INTENT(in)                      :: rot_coeff(:,:,:,:)
  ! !   REAL(wp), INTENT(inout)                   :: rot_vec_v(:,:,:)
  ! !
  ! ! !Local variables
  ! ! !
  ! ! REAL(wp) :: z_vort_tmp
  ! ! REAL(wp) :: z_vort_tmp_boundary(nproma,n_zlev,p_patch%nblks_v)
  ! ! REAL(wp) :: z_vort_tmp_boundary2(nproma,n_zlev,p_patch%nblks_v)
  ! ! !REAL(wp) :: z_weight(nproma,n_zlev,p_patch%nblks_v)
  ! ! REAL(wp) :: zarea_fraction
  ! ! REAL(wp) :: z_area_scaled
  ! !
  ! ! INTEGER :: slev, elev     ! vertical start and end level
  ! ! INTEGER :: jv, jk, jb, jev,je
  ! ! INTEGER :: ile, ibe, il, ib, ill
  ! ! INTEGER :: rl_start_e, rl_end_e
  ! ! INTEGER :: i_startblk_v, i_endblk_v, i_startidx_v, i_endidx_v
  ! ! INTEGER :: i_startblk_e, i_endblk_e, i_startidx_e, i_endidx_e
  ! !
  ! ! !INTEGER :: i_bdr_ctr
  ! ! !INTEGER :: icell_idx_1, icell_blk_1
  ! ! !INTEGER :: icell_idx_2, icell_blk_2
  ! ! INTEGER :: il_v1, il_v2,ib_v1, ib_v2
  ! ! INTEGER  :: i_v_ctr(nproma,n_zlev,p_patch%nblks_v)
  ! ! INTEGER  :: i_v_bnd_edge_ctr(nproma,n_zlev,p_patch%nblks_v)
  ! ! INTEGER  :: ibnd_edge_idx(4), ibnd_edge_blk(4)  !maximal 4 boundary edges in a dual loop.
  ! ! INTEGER  :: i_edge_idx(4)
  ! ! REAL(wp) :: z_orientation(4),temp1,temp2
  ! ! REAL(wp) :: z_vt(nproma,n_zlev,p_patch%nblks_e)
  ! ! TYPE(t_cartesian_coordinates) :: vertex_cc! cell1_cc, cell2_cc
  ! ! INTEGER,PARAMETER :: ino_dual_edges = 6
  ! !
  ! ! INTEGER,PARAMETER :: rl_start_v = 2
  ! ! INTEGER,PARAMETER :: rl_end_v   = min_rlvert
  ! ! !-----------------------------------------------------------------------
  ! ! slev         = 1
  ! ! elev         = n_zlev
  ! ! rl_start_e   = 1
  ! ! rl_end_e     = min_rledge
  ! !
  ! ! i_startblk_e = p_patch%edges%start_blk(rl_start_e,1)
  ! ! i_endblk_e   = p_patch%edges%end_blk(rl_end_e,1)
  ! !
  ! ! i_startblk_v = p_patch%verts%start_blk(rl_start_v,1)
  ! ! i_endblk_v   = p_patch%verts%end_blk(rl_end_v,1)
  ! !
  ! ! ! #slo# due to nag -nan compiler-option
  ! ! i_v_ctr(:,:,:)          = 0
  ! ! i_v_bnd_edge_ctr(:,:,:) = 0
  ! ! ibnd_edge_idx(1:4)      = 0
  ! ! ibnd_edge_blk(1:4)      = 0
  ! ! i_edge_idx(1:4)         = 0
  ! ! z_orientation(1:4)      = 0.0_wp
  ! ! z_vort_tmp_boundary     = 0.0_wp
  ! ! z_vort_tmp_boundary2    = 0.0_wp
  ! ! rot_vec_v(:,:,:)        = 0.0_wp
  ! ! z_vt(:,:,:)             = 0.0_wp
  ! !
  ! !
  ! ! !In this loop vorticity at vertices is calculated
  ! ! DO jb = i_startblk_v, i_endblk_v
  ! !   CALL get_indices_v(p_patch, jb, i_startblk_v, i_endblk_v, &
  ! !                      i_startidx_v, i_endidx_v, rl_start_v, rl_end_v)
  ! !   DO jk = slev, elev
  ! ! !!$OMP PARALLEL DO SCHEDULE(runtime) DEFAULT(PRIVATE)  &
  ! ! !!$OMP   SHARED(u_vec_e,v_vec_e,ptr_patch,rot_vec_v,jb) FIRSTPRIVATE(jk)
  ! !     DO jv = i_startidx_v, i_endidx_v
  ! !
  ! !       z_vort_tmp          = 0.0_wp
  ! !       !zarea_fraction      = 0.0_wp
  ! !       !i_bdr_ctr           = 0
  ! !       !z_weight(jv,jk,jb) = 0.0_wp
  ! !
  ! !       ibnd_edge_idx(1:4)      = 0
  ! !       ibnd_edge_blk(1:4)      = 0
  ! !       i_edge_idx(1:4)         = 0
  ! !       z_orientation(1:4)      = 0.0_wp
  ! !
  ! !       DO jev = 1, p_patch%verts%num_edges(jv,jb)
  ! !
  ! !         ! get line and block indices of edge jev around vertex jv
  ! !         ile = p_patch%verts%edge_idx(jv,jb,jev)
  ! !         ibe = p_patch%verts%edge_blk(jv,jb,jev)
  ! !
  ! !
  ! !         IF ( v_base%lsm_oce_e(ile,jk,ibe) == sea) THEN
  ! !           !Distinguish the following cases
  ! !           ! edge ie_k is
  ! !           !a) ocean edge: compute as usual,
  ! !           !b) land edge: do not consider it
  ! !           !c) boundary edge take:
  ! !           !  no-slip boundary condition:  normal and tangential velocity at boundary are zero
  ! !           ! sea, sea_boundary, boundary (edges only), land_boundary, land =
  ! !           !  -2,      -1,         0,                  1,             2
  ! !           !add contribution of normal velocity at edge (ile,ibe) to rotation
  ! !           !z_vort_tmp = z_vort_tmp + vn(ile,jk,ibe)                   &
  ! !           !  & * p_patch%edges%dual_edge_length(ile,ibe)  &
  ! !           !  & * p_patch%verts%edge_orientation(jv,jb,jev)
  ! !           z_vort_tmp = z_vort_tmp + vn(ile,jk,ibe)*rot_coeff(jv,jk,jb,jev)
  ! !
  ! !           !increase wet edge ctri_v_ctr(jv,jk,jb)=i_v_ctr(jv,jk,jb)+1
  ! !
  ! !          ! edge with indices ile, ibe is boundary edge
  ! !          ELSE IF ( v_base%lsm_oce_e(ile,jk,ibe) == boundary ) THEN
  ! !
  ! !           !calculate tangential velocity
  ! !           il_v1 = p_patch%edges%vertex_idx(ile,ibe,1)
  ! !           ib_v1 = p_patch%edges%vertex_blk(ile,ibe,1)
  ! !           il_v2 = p_patch%edges%vertex_idx(ile,ibe,2)
  ! !           ib_v2 = p_patch%edges%vertex_blk(ile,ibe,2)
  ! !
  ! !           z_vt(ile,jk,ibe)= &
  ! !           &- DOT_PRODUCT(p_vn_dual(il_v1,jk,ib_v1)%x,&
  ! !           &p_int_state(1)%edge2vert_coeff_cc_t(ile,ibe,2)%x)&
  ! !           &+ DOT_PRODUCT(p_vn_dual(il_v2,jk,ib_v2)%x,&
  ! !           &p_int_state(1)%edge2vert_coeff_cc_t(ile,ibe,1)%x)
  ! !
  ! !           !increase boundary edge counter
  ! !            i_v_bnd_edge_ctr(jv,jk,jb)=i_v_bnd_edge_ctr(jv,jk,jb)+1
  ! !
  ! !            !Store actual boundary edge indices
  ! !            IF(i_v_bnd_edge_ctr(jv,jk,jb)==1)THEN
  ! !              ibnd_edge_idx(1) = ile
  ! !              ibnd_edge_blk(1) = ibe
  ! !              z_orientation(1) = p_patch%verts%edge_orientation(jv,jb,jev)
  ! !              i_edge_idx(1)    = jev
  ! !            ELSEIF(i_v_bnd_edge_ctr(jv,jk,jb)==2)THEN
  ! !              ibnd_edge_idx(2) = ile
  ! !              ibnd_edge_blk(2) = ibe
  ! !              z_orientation(2) = p_patch%verts%edge_orientation(jv,jb,jev)
  ! !              i_edge_idx(2)    = jev
  ! !            ELSEIF(i_v_bnd_edge_ctr(jv,jk,jb)==3)THEN
  ! !              ibnd_edge_idx(3) = ile
  ! !              ibnd_edge_blk(3) = ibe
  ! !              z_orientation(3) = p_patch%verts%edge_orientation(jv,jb,jev)
  ! !              i_edge_idx(3)    = jev
  ! !            ELSEIF(i_v_bnd_edge_ctr(jv,jk,jb)==4)THEN
  ! !              ibnd_edge_idx(4) = ile
  ! !              ibnd_edge_blk(4) = ibe
  ! !              z_orientation(4) = p_patch%verts%edge_orientation(jv,jb,jev)
  ! !              i_edge_idx(4)    = jev
  ! !            ELSE
  ! !            !maximal 4 boundary edges per dual loop are allowed: somethings wrong withe the grid
  ! !            ! write(*,*)'grid error',jv,jk,jb,i_v_bnd_edge_ctr(jv,jk,jb)
  ! !            CALL message (TRIM('sbr nonlinear Coriolis'), &
  ! !            &'more than 4 boundary edges per dual loop: something is wrong with the grid')
  ! !            CALL finish ('TRIM(sbr nonlinear Coriolis)','Grid-boundary error !!')
  ! !            ENDIF
  ! !          END IF
  ! !       END DO
  ! !
  ! !      !Finalize vorticity calculation by closing the dual loop along boundary edges
  ! !      IF(i_v_bnd_edge_ctr(jv,jk,jb)==2)THEN!(i_v_ctr(jv,jk,jb)==2)THEN
  ! ! !            z_vort_tmp_boundary2(jv,jk,jb) =&
  ! ! !            !& p_patch%edges%system_orientation(ibnd_edge_idx(1),ibnd_edge_blk(1))*&
  ! ! !            &0.5_wp*z_orientation(1)*&
  ! ! !            &z_vt(ibnd_edge_idx(1),jk,ibnd_edge_blk(1)) &
  ! ! !            &*p_patch%edges%primal_edge_length(ibnd_edge_idx(1),ibnd_edge_blk(1))&
  ! ! !            &+&
  ! ! !            !& p_patch%edges%system_orientation(ibnd_edge_idx(2),ibnd_edge_blk(2))*&
  ! ! !            &0.5_wp*z_orientation(2)*&
  ! ! !            &z_vt(ibnd_edge_idx(2),jk,ibnd_edge_blk(2)) &
  ! ! !            &*p_patch%edges%primal_edge_length(ibnd_edge_idx(2),ibnd_edge_blk(2))
  ! !
  ! !        z_vort_tmp_boundary(jv,jk,jb)=&
  ! !        &z_vt(ibnd_edge_idx(1),jk,ibnd_edge_blk(1))*rot_coeff(jv,jk,jb,i_edge_idx(1))&
  ! !        &+&
  ! !        &z_vt(ibnd_edge_idx(2),jk,ibnd_edge_blk(2))*rot_coeff(jv,jk,jb,i_edge_idx(2))
  ! !
  ! !      ELSEIF(i_v_bnd_edge_ctr(jv,jk,jb)==4)THEN
  ! !            !In case of 4 boundary edges within a dual loop, we have 2 land triangles
  ! !            !around the vertex. these two land triangles have one vertex in common and are
  ! !            !seperated by two wet triangles.
  ! ! !            z_vort_tmp_boundary2(jv,jk,jb) =&
  ! ! !            &0.5_wp*z_orientation(1)*&
  ! ! !            &z_vt(ibnd_edge_idx(1),jk,ibnd_edge_blk(1)) &
  ! ! !            &*p_patch%edges%primal_edge_length(ibnd_edge_idx(1),ibnd_edge_blk(1))&
  ! ! !            &+&
  ! ! !            &0.5_wp*z_orientation(2)*&
  ! ! !            &z_vt(ibnd_edge_idx(2),jk,ibnd_edge_blk(2)) &
  ! ! !            &*p_patch%edges%primal_edge_length(ibnd_edge_idx(2),ibnd_edge_blk(2))&
  ! ! !            &+&
  ! ! !            &0.5_wp*z_orientation(3)*&
  ! ! !            &z_vt(ibnd_edge_idx(3),jk,ibnd_edge_blk(3)) &
  ! ! !            &*p_patch%edges%primal_edge_length(ibnd_edge_idx(3),ibnd_edge_blk(3))&
  ! ! !            &+&
  ! ! !            &0.5_wp*z_orientation(4)*&
  ! ! !            &z_vt(ibnd_edge_idx(4),jk,ibnd_edge_blk(4)) &
  ! ! !            &*p_patch%edges%primal_edge_length(ibnd_edge_idx(4),ibnd_edge_blk(4))
  ! !
  ! !        z_vort_tmp_boundary(jv,jk,jb)=&
  ! !        &z_vt(ibnd_edge_idx(1),jk,ibnd_edge_blk(1))*rot_coeff(jv,jk,jb,i_edge_idx(1))&
  ! !        &+&
  ! !        &z_vt(ibnd_edge_idx(2),jk,ibnd_edge_blk(2))*rot_coeff(jv,jk,jb,i_edge_idx(2))&
  ! !        &+&
  ! !        &z_vt(ibnd_edge_idx(3),jk,ibnd_edge_blk(3))*rot_coeff(jv,jk,jb,i_edge_idx(3))&
  ! !        &+&
  ! !        &z_vt(ibnd_edge_idx(4),jk,ibnd_edge_blk(4))*rot_coeff(jv,jk,jb,i_edge_idx(4))
  ! !
  ! !      ENDIF
  ! !
  ! !      !Fianl vorticity calculation
  ! !      rot_vec_v(jv,jk,jb) = z_vort_tmp+z_vort_tmp_boundary(jv,jk,jb)
  ! !
  ! !     END DO
  ! ! !!$OMP END PARALLEL DO
  ! !   END DO
  ! ! END DO
  ! ! DO jb = i_startblk_v, i_endblk_v
  ! !   CALL get_indices_v(p_patch, jb, i_startblk_v, i_endblk_v, &
  ! !                      i_startidx_v, i_endidx_v, rl_start_v, rl_end_v)
  ! !   DO jk = slev, elev
  ! !     DO jv = i_startidx_v, i_endidx_v
  ! ! IF(rot_vec_v(jv,jk,jb)/=0.0_wp)THEN
  ! ! write(123456,*)'rot 3D:',jk,jv,jb,rot_vec_v(jv,jk,jb),&
  ! ! &z_vort_tmp_boundary(jv,jk,jb), i_v_bnd_edge_ctr(jv,jk,jb)!,i_v_ctr(jv,jk,jb)
  ! ! ENDIF
  ! !     END DO
  ! !   END DO
  ! ! END DO
  ! !
  ! !
  ! ! END SUBROUTINE rot_vertex_ocean_3D_backup
END MODULE mo_oce_math_operators

