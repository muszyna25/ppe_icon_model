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
  USE mo_exception,          ONLY: finish,message
  USE mo_run_config,         ONLY: ltimer, dtime
  USE mo_math_constants
  USE mo_physical_constants
  USE mo_impl_constants,     ONLY: boundary, sea_boundary !,sea,land, land_boundary, sea, max_char_length, &
  USE mo_model_domain,       ONLY: t_patch, t_patch_3D
  USE mo_ext_data_types,     ONLY: t_external_data
  USE mo_ocean_nml,          ONLY: n_zlev, iswm_oce
  USE mo_dynamics_config,    ONLY: nold
  USE mo_util_dbg_prnt,      ONLY: dbg_print
#ifndef __SX__
  USE mo_timer,              ONLY: timer_start, timer_stop, timer_div, timer_grad
#endif
  USE mo_oce_types,          ONLY: t_hydro_ocean_state
  USE mo_math_utilities,     ONLY: t_cartesian_coordinates, vector_product !, gc2cc
  USE mo_operator_ocean_coeff_3d, ONLY: t_operator_coeff
  USE mo_grid_subset,         ONLY: t_subset_range, get_index_range
  USE mo_sync,                ONLY: SYNC_C, SYNC_E, SYNC_V, sync_patch_array
  USE mo_grid_config,         ONLY: n_dom

  IMPLICIT NONE

  PRIVATE

  CHARACTER(LEN=*), PARAMETER :: version = '$Id$'
  CHARACTER(len=12)           :: str_module    = 'oceMathOps  '  ! Output of module for 1 line debug
  INTEGER                     :: idt_src       = 1               ! Level of detail for 1 line debug


  PUBLIC :: grad_fd_norm_oce_3d
  PUBLIC :: div_oce_3d
  PUBLIC :: rot_vertex_ocean_3d
  PUBLIC :: grad_fd_norm_oce_2d_3d
  !PUBLIC :: rot_vertex_ocean
  !PUBLIC :: rot_vertex_ocean_rbf
  PUBLIC :: calculate_thickness
  PUBLIC :: map_edges2vert_3D
  PUBLIC :: check_cfl_horizontal, check_cfl_vertical


  INTERFACE div_oce_3d
    MODULE PROCEDURE div_oce_3d_mlevels
    MODULE PROCEDURE div_oce_3d_1level
  END INTERFACE

CONTAINS


  !-------------------------------------------------------------------------
  !>
  SUBROUTINE map_edges2vert_3d(patch_2D, vn, edge2vert_coeff_cc, vn_dual)
    
    TYPE(t_patch), TARGET, INTENT(in)       :: patch_2D
    REAL(wp), INTENT(in)                    :: vn(:,:,:)
    TYPE(t_cartesian_coordinates),INTENT(in):: edge2vert_coeff_cc(:,:,:,:)
    TYPE(t_cartesian_coordinates)           :: vn_dual(nproma,n_zlev,patch_2D%nblks_v)

    !Local variables
    !
    INTEGER :: slev, elev     ! vertical start and end level
    INTEGER :: jv, jk, jb,jev
    INTEGER :: ile, ibe
    INTEGER :: i_startidx_v, i_endidx_v
    TYPE(t_subset_range), POINTER :: verts_in_domain
    !-----------------------------------------------------------------------
    verts_in_domain => patch_2D%verts%in_domain

    !i_v_ctr(:,:,:) = 0
    slev         = 1
    elev         = n_zlev

    DO jb = verts_in_domain%start_block, verts_in_domain%end_block
      CALL get_index_range(verts_in_domain, jb, i_startidx_v, i_endidx_v)
      DO jk = slev, elev
        DO jv = i_startidx_v, i_endidx_v

          vn_dual(jv,jk,jb)%x = 0.0_wp
          DO jev = 1, patch_2D%verts%num_edges(jv,jb)

            ! get line and block indices of edge jev around vertex jv
            ile = patch_2D%verts%edge_idx(jv,jb,jev)
            ibe = patch_2D%verts%edge_blk(jv,jb,jev)

            vn_dual(jv,jk,jb)%x = vn_dual(jv,jk,jb)%x        &
            & +edge2vert_coeff_cc(jv,jk,jb,jev)%x &
            & *vn(ile,jk,ibe)
!ENDIF
          END DO
        END DO ! jv = i_startidx_v, i_endidx_v
      END DO ! jk = slev, elev
    END DO ! jb = verts_in_domain%start_block, verts_in_domain%end_block

    ! sync the result
    CALL sync_patch_array(SYNC_V, patch_2D, vn_dual(:,:,:)%x(1))
    CALL sync_patch_array(SYNC_V, patch_2D, vn_dual(:,:,:)%x(2))
    CALL sync_patch_array(SYNC_V, patch_2D, vn_dual(:,:,:)%x(3))

  END SUBROUTINE map_edges2vert_3d
  !-------------------------------------------------------------------------------------

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
  !! - abandon grid for the sake of patch_2D
  !!Boundary handling for triangles by P. Korn (2009)
  !!
  !!  mpi note: the result is on edges_in_domain.
  SUBROUTINE grad_fd_norm_oce_3d( psi_c, patch_3D, grad_coeff, grad_norm_psi_e)

    TYPE(t_patch_3D ),TARGET, INTENT(IN)   :: patch_3D
    REAL(wp), INTENT(IN)                   :: grad_coeff(:,:,:)!(nproma,n_zlev,patch_3D%p_patch_2D(1)%nblks_e)
    REAL(wp), INTENT(IN)                   :: psi_c          (nproma,n_zlev,patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    REAL(wp), INTENT(INOUT)                :: grad_norm_psi_e(nproma,n_zlev,patch_3D%p_patch_2D(1)%nblks_e)

    !
    !Local variables
    !
    INTEGER :: slev, elev     ! vertical start and end level
    INTEGER :: je, jk, jb
    INTEGER :: i_startidx, i_endidx
    INTEGER,  DIMENSION(:,:,:), POINTER :: iidx, iblk
    TYPE(t_subset_range), POINTER       :: edges_in_domain
    TYPE(t_patch), POINTER              :: patch_2D
    !-----------------------------------------------------------------------
    patch_2D           => patch_3D%p_patch_2D(1)
    edges_in_domain => patch_2D%edges%in_domain

    slev = 1
    elev = n_zlev

    iidx => patch_2D%edges%cell_idx
    iblk => patch_2D%edges%cell_blk
    !
    !  loop through all patch_2D edges (and blocks)
    !
#ifndef __SX__
    IF (ltimer) CALL timer_start(timer_grad)
#endif
! !$OMP PARALLEL
    ! The special treatment of 2D fields is essential for efficiency on the NEC

#ifdef __SX__
    IF (slev > 1) THEN
! !$OMP DO PRIVATE(jb,i_startidx,i_endidx,je,jk)
    DO jb = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, jb, i_startidx, i_endidx)
#ifdef _URD2
!CDIR UNROLL=_URD2
#endif
        DO jk = slev, elev
          DO je = i_startidx, i_endidx
            !IF (patch_3D%lsm_e(je,jk,jb) <= sea_boundary) THEN
              grad_norm_psi_e(je,jk,jb) = grad_coeff(je,jk,jb)* &
              & ( psi_c(iidx(je,jb,2),jk,iblk(je,jb,2))-        &
              & psi_c(iidx(je,jb,1),jk,iblk(je,jb,1)) )
            !ENDIF
          ENDDO
        END DO
      END DO
! !$OMP END DO
    ELSE
#endif

! !$OMP DO PRIVATE(jb,i_startidx,i_endidx,je,jk)
    DO jb = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, jb, i_startidx, i_endidx)
#ifdef _URD
!CDIR UNROLL=_URD
#endif
        DO jk = slev, elev
          DO je = i_startidx, i_endidx
            !IF (patch_3D%lsm_e(je,jk,jb) <= sea_boundary) THEN
              grad_norm_psi_e(je,jk,jb) = grad_coeff(je,jk,jb)* &
                & ( psi_c(iidx(je,jb,2),jk,iblk(je,jb,2)) -     &
                & psi_c(iidx(je,jb,1),jk,iblk(je,jb,1)) )
            !ENDIF 
          ENDDO
        END DO
      END DO
! !$OMP END DO
#ifdef __SX__
    ENDIF
#endif


! !$OMP END PARALLEL
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
  !! - loop only over the inner cells of a patch_2D, not any more over halo cells
  !! Modifications by P. Korn, MPI-M(2007-2)
  !! - Switch fom array arguments to pointers
  !! Modification by Almut Gassmann, MPI-M (2007-04-20)
  !! - abandon grid for the sake of patch_2D
  !! Modification by Guenther Zaengl, DWD (2009-03-17)
  !! - vector optimization
  !! Modification by Peter Korn, MPI-M    (2009)
  !! - Boundary treatment for the ocean
  !! Modification by Stephan Lorenz, MPI-M (2010-08-05)
  !! - New boundary definition with inner and boundary points on land/sea
  !!
  !!  mpi parallelized LL (no sync required)
  SUBROUTINE div_oce_3d_mlevels( vec_e, patch_2D,div_coeff, div_vec_c, opt_slev, opt_elev, &
    & subset_range)
    !
    !
    !  patch_2D on which computation is performed
    !
    TYPE(t_patch), TARGET, INTENT(in) :: patch_2D
    !
    ! edge based variable of which divergence
    ! is computed
    !
    REAL(wp), INTENT(in)          :: vec_e(:,:,:) ! dim: (nproma,n_zlev,nblks_e)
    REAL(wp), INTENT(in)          :: div_coeff(:,:,:,:)
    REAL(wp), INTENT(inout)       :: div_vec_c(:,:,:) ! dim: (nproma,n_zlev,alloc_cell_blocks)
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
      all_cells => patch_2D%cells%all
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
! !$OMP PARALLEL

    iidx => patch_2D%cells%edge_idx
    iblk => patch_2D%cells%edge_blk

! !$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk)
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
          ! patch_2D%grid%cells%edge_orientation)
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
! !$OMP END DO

! !$OMP END PARALLEL
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
  !! - loop only over the inner cells of a patch_2D, not any more over halo cells
  !! Modifications by P. Korn, MPI-M(2007-2)
  !! - Switch fom array arguments to pointers
  !! Modification by Almut Gassmann, MPI-M (2007-04-20)
  !! - abandon grid for the sake of patch_2D
  !! Modification by Guenther Zaengl, DWD (2009-03-17)
  !! - vector optimization
  !! Modification by Peter Korn, MPI-M    (2009)
  !! - Boundary treatment for the ocean
  !! Modification by Stephan Lorenz, MPI-M (2010-08-05)
  !! - New boundary definition with inner and boundary points on land/sea
  !!
  SUBROUTINE div_oce_3d_1level( vec_e, patch_2D, div_coeff, div_vec_c,  &
    & level, subset_range)
    !
    !
    !  patch_2D on which computation is performed
    !
    TYPE(t_patch), TARGET, INTENT(in) :: patch_2D
    !
    ! edge based variable of which divergence
    ! is computed
    !
    REAL(wp), INTENT(inout)       :: vec_e(:,:) ! dim: (nproma,n_zlev,nblks_e)
    REAL(wp), INTENT(in)          :: div_coeff(:,:,:,:)
    REAL(wp), INTENT(inout)       :: div_vec_c(:,:) ! dim: (nproma,n_zlev,alloc_cell_blocks)
    INTEGER,  INTENT(in)          :: level
    TYPE(t_subset_range), TARGET, INTENT(in), OPTIONAL :: subset_range

    INTEGER :: jc, jb
    INTEGER :: i_startidx, i_endidx
    INTEGER,  DIMENSION(:,:,:),   POINTER :: iidx, iblk
    TYPE(t_subset_range), POINTER :: all_cells
    !-----------------------------------------------------------------------
    IF (PRESENT(subset_range)) THEN
      all_cells => subset_range
    ELSE
      all_cells => patch_2D%cells%all
    ENDIF

#ifndef __SX__
    IF (ltimer) CALL timer_start(timer_div)
#endif
! !$OMP PARALLEL

    iidx => patch_2D%cells%edge_idx
    iblk => patch_2D%cells%edge_blk

! !$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc)
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
          ! patch_2D%grid%cells%edge_orientation)
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
! !$OMP END DO

! !$OMP END PARALLEL
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
  !! - abandon grid for the sake of patch_2D
  !! Boundary handling for triangles by P. Korn (2009)
  !!  mpi note: the result is not synced. Should be done in the calling method if required
  !!
  SUBROUTINE grad_fd_norm_oce_2d_3d( psi_c, patch_2D, grad_coeff, grad_norm_psi_e)
    !
    TYPE(t_patch), TARGET, INTENT(in) :: patch_2D
    REAL(wp), INTENT(in)    :: psi_c(:,:)             ! dim: (nproma,alloc_cell_blocks)
    REAL(wp), INTENT(in)    :: grad_coeff(:,:)
    REAL(wp), INTENT(inout) ::  grad_norm_psi_e(:,:)  ! dim: (nproma,nblks_e)

    !!
    !!local variables
    INTEGER :: slev, elev     ! vertical start and end level
    INTEGER :: je, jb
    INTEGER :: i_startidx, i_endidx
    INTEGER,  DIMENSION(:,:,:),   POINTER :: iidx, iblk
    TYPE(t_subset_range), POINTER :: edges_in_domain
    !-----------------------------------------------------------------------
    edges_in_domain => patch_2D%edges%in_domain
    slev = 1
    elev = 1

    iidx => patch_2D%edges%cell_idx
    iblk => patch_2D%edges%cell_blk
        !
    !  loop through all patch_2D edges (and blocks)
    !
#ifndef __SX__
    IF (ltimer) CALL timer_start(timer_grad)
#endif
! !$OMP PARALLEL
    ! The special treatment of 2D fields is essential for efficiency on the NEC
! !$OMP DO PRIVATE(jb,i_startidx,i_endidx,je)
    DO jb = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, jb, i_startidx, i_endidx)
#ifdef _URD
!CDIR UNROLL=_URD
#endif

      DO je = i_startidx, i_endidx
        ! compute the normal derivative
        ! by the finite difference approximation
        ! (see Bonaventura and Ringler MWR 2005)
        IF (iidx(je,jb,1) < 1 .or. iidx(je,jb,2) < 1) THEN
!          IF (grad_coeff(je,jb) /= 0.0_wp) THEN
!            CALL finish("grad_fd_norm_oce_2d_3d","grad_coeff(je,jb) /= 0.0_wp")
!          ELSE
            grad_norm_psi_e(je,jb) = 0.0_wp
!          ENDIF
        ELSE
          grad_norm_psi_e(je,jb) =  &
            & (psi_c(iidx(je,jb,2),iblk(je,jb,2))-psi_c(iidx(je,jb,1),iblk(je,jb,1)))&
            & * grad_coeff(je,jb)
        ENDIF

      END DO
    END DO
! !$OMP END DO
! !$OMP END PARALLEL
#ifndef __SX__
    IF (ltimer) CALL timer_stop(timer_grad)
#endif
  END SUBROUTINE grad_fd_norm_oce_2d_3d
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
  !! does the tangential velocity calculate from velocity vector at vertices (vn_dual), while RBF uses
  !! a specific routine for that purpose.
  !!
  !! mpi note: the results is not synced. should be done by the calling method if necessary
  !!     vn, vn_dual must have been synced on level 2 (in_domain + 1)
  SUBROUTINE rot_vertex_ocean_3d( patch_3D, vn, vn_dual, p_op_coeff, rot_vec_v)
    !>
    !!
    TYPE(t_patch_3D ),TARGET, INTENT(IN)      :: patch_3D
    REAL(wp), INTENT(in)                      :: vn(:,:,:)
    TYPE(t_cartesian_coordinates), INTENT(in) :: vn_dual(nproma,n_zlev,patch_3D%p_patch_2D(1)%nblks_v)
    TYPE(t_operator_coeff),TARGET, INTENT(in) :: p_op_coeff
    REAL(wp), INTENT(inout)                   :: rot_vec_v(nproma,n_zlev,patch_3D%p_patch_2D(1)%nblks_v)

    !Local variables
    !
    REAL(wp) :: z_vort_int
    REAL(wp) :: z_vort_boundary(nproma,n_zlev,patch_3D%p_patch_2D(1)%nblks_v)
    INTEGER :: slev, elev
    INTEGER :: jv, jk, jb, jev
    INTEGER :: ile, ibe
    INTEGER :: il_v1, il_v2,ib_v1, ib_v2
    INTEGER :: i_startidx_v, i_endidx_v

    REAL(wp) :: z_vt(nproma,n_zlev,patch_3D%p_patch_2D(1)%nblks_e)
    INTEGER, POINTER :: ibnd_edge_idx(:,:,:,:), ibnd_edge_blk(:,:,:,:), i_edge_idx(:,:,:,:)
    !REAL(wp), POINTER :: z_orientation(:,:,:,:)

    TYPE(t_subset_range), POINTER :: verts_in_domain
    TYPE(t_patch), POINTER        :: patch_2D
    !-----------------------------------------------------------------------
    patch_2D           => patch_3D%p_patch_2D(1)
    verts_in_domain => patch_2D%verts%in_domain
    slev         = 1
    elev         = n_zlev


    z_vort_boundary  = 0.0_wp
    rot_vec_v(:,:,:) = 0.0_wp
    z_vt(:,:,:)      = 0.0_wp

    !set pointer that carry edge information
    ibnd_edge_idx    => p_op_coeff%bnd_edge_idx
    ibnd_edge_blk    => p_op_coeff%bnd_edge_blk
    i_edge_idx       => p_op_coeff%edge_idx
    !z_orientation    => p_op_coeff%orientation


    !In this loop tangential velocity is calulated at boundaries
    DO jb = verts_in_domain%start_block, verts_in_domain%end_block
      CALL get_index_range(verts_in_domain, jb, i_startidx_v, i_endidx_v)
      DO jk = slev, elev
        !!$OMP PARALLEL DO SCHEDULE(runtime) DEFAULT(PRIVATE)  &
        !!$OMP   SHARED(u_vec_e,v_vec_e,patch_2D,rot_vec_v,jb) FIRSTPRIVATE(jk)
        DO jv = i_startidx_v, i_endidx_v

          DO jev = 1, patch_2D%verts%num_edges(jv,jb)

            ! get line and block indices of edge jev around vertex jv
            ile = patch_2D%verts%edge_idx(jv,jb,jev)
            ibe = patch_2D%verts%edge_blk(jv,jb,jev)

            !IF ( v_base%lsm_e(ile,jk,ibe) == boundary ) THEN
            IF(patch_3D%lsm_e(ile,jk,ibe) == boundary)THEN
              !calculate tangential velocity
              il_v1 = patch_2D%edges%vertex_idx(ile,ibe,1)
              ib_v1 = patch_2D%edges%vertex_blk(ile,ibe,1)
              il_v2 = patch_2D%edges%vertex_idx(ile,ibe,2)
              ib_v2 = patch_2D%edges%vertex_blk(ile,ibe,2)

              z_vt(ile,jk,ibe)= &
                & - DOT_PRODUCT(vn_dual(il_v1,jk,ib_v1)%x,&
                & p_op_coeff%edge2vert_coeff_cc_t(ile,jk,ibe,1)%x)&
                & + DOT_PRODUCT(vn_dual(il_v2,jk,ib_v2)%x,&
                & p_op_coeff%edge2vert_coeff_cc_t(ile,jk,ibe,2)%x)
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
        !!$OMP   SHARED(u_vec_e,v_vec_e,patch_2D,rot_vec_v,jb) FIRSTPRIVATE(jk)
        DO jv = i_startidx_v, i_endidx_v

          z_vort_int = 0.0_wp

          DO jev = 1, patch_2D%verts%num_edges(jv,jb)

            ! get line and block indices of edge jev around vertex jv
            ile = patch_2D%verts%edge_idx(jv,jb,jev)
            ibe = patch_2D%verts%edge_blk(jv,jb,jev)
            !add contribution of normal velocity at edge (ile,ibe) to rotation
            !IF ( v_base%lsm_e(ile,jk,ibe) == sea) THEN
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
    !         & patch_2D%verts%dual_area(jv,jb)
          rot_vec_v(jv,jk,jb) = z_vort_int + z_vort_boundary(jv,jk,jb)
          
        END DO
        !!$OMP END PARALLEL DO
      END DO
    END DO

    ! DO jb = i_startblk_v, i_endblk_v
    !   CALL get_indices_v(patch_2D, jb, i_startblk_v, i_endblk_v, &
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
  !---------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------
  !>
  !!
  !!  Calculation of total fluid thickness at cell centers and surface elevation at
  !!  cell edges from prognostic surface height at cell centers. We use height at
  !!  old timelevel "n"
  !!
  !! @par Revision History
  !! Developed  by  Peter Korn, MPI-M (2010).
  !!
  SUBROUTINE calculate_thickness( patch_3D, ocean_state, p_ext_data, operators_coefficients)
  !SUBROUTINE calculate_thickness( p_patch_3D, ocean_state, p_ext_data, ice_hi)
    !
    ! patch_2D on which computation is performed
    TYPE(t_patch_3D ),TARGET, INTENT(IN)   :: patch_3D
    !
    ! Type containing ocean state
    TYPE(t_hydro_ocean_state), TARGET :: ocean_state
    !
    ! Type containing external data
    TYPE(t_external_data), TARGET, INTENT(in) :: p_ext_data
    !REAL(wp), INTENT(IN)                      :: ice_hi(nproma,1,p_patch_3D%p_patch_2D(1)%nblks_c)
    TYPE(t_operator_coeff), INTENT(in)      :: operators_coefficients

    !  local variables
    INTEGER            :: i_startidx_c, i_endidx_c
    INTEGER            :: edge_start_idx, edge_end_idx
    INTEGER            :: jc, jb, je
    INTEGER            :: il_c1, ib_c1, il_c2, ib_c2
    REAL(wp)           :: z_dist_e_c1, z_dist_e_c2
    TYPE(t_subset_range), POINTER :: all_cells, all_edges, edges_in_domain
    TYPE(t_patch), POINTER        :: patch_2D

    INTEGER :: cell_1_index, cell_2_index, cell_1_block, cell_2_block
    REAL(wp) :: top_vn_1, top_vn_2, integrated_vn
    REAL(wp), POINTER :: top_coeffs(:,:,:), integrated_coeffs(:,:,:), sum_to_2D_coeffs(:,:,:)
    REAL(wp), POINTER :: cell_thickeness(:,:,:), edge_thickeness(:,:,:)
    REAL(wp), POINTER :: inv_cell_thickeness(:,:,:), inv_edge_thickeness(:,:,:)
    REAL(wp), POINTER :: inv_prisms_center_distance(:,:,:), inv_faces_center_distance(:,:,:)
    REAL(wp)  :: cell_thickeness_1, cell_thickeness_2
    !-------------------------------------------------------------------------------
    !CALL message (TRIM(routine), 'start')
    patch_2D            => patch_3D%p_patch_2D(1)
    all_cells           => patch_2D%cells%all
    all_edges           => patch_2D%edges%all
    edges_in_domain     => patch_2D%edges%in_domain
    cell_thickeness     => patch_3D%p_patch_1d(n_dom)%prism_thick_c
    inv_cell_thickeness => patch_3D%p_patch_1d(n_dom)%inv_prism_thick_c
    inv_prisms_center_distance => patch_3d%p_patch_1d(1)%inv_prism_center_dist_c
    edge_thickeness     => patch_3D%p_patch_1d(n_dom)%prism_thick_e
    inv_edge_thickeness => patch_3D%p_patch_1d(n_dom)%inv_prism_thick_e
    inv_faces_center_distance => patch_3d%p_patch_1d(1)%inv_prism_center_dist_e

 
     ! already done after update fluxes
!    CALL sync_patch_array(sync_c, patch_2D, ocean_state%p_prog(nold(1))%h)

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
          IF(patch_3D%lsm_c(jc,1,jb) <= sea_boundary)THEN

            ocean_state%p_diag%thick_c(jc,jb) = ocean_state%p_prog(nold(1))%h(jc,jb)&
              &  - p_ext_data%oce%bathymetry_c(jc,jb)
     !        &  - ice_hi(jc,1,jb)
          ELSE
            ocean_state%p_diag%thick_c(jc,jb) = 0.0_wp
          ENDIF
        END DO
      END DO!write(*,*)'bathymetry',maxval(p_ext_data%oce%bathymetry_c),minval(p_ext_data%oce%bathymetry_c)

    ELSEIF(iswm_oce /= 1 )THEN

      DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
        !calculate for each fluid colum the total depth, i.e.
        !from bottom boundary to surface height
        !the bottom boundary is zlev_i(dolic+1) since zlev_i(1)=0 (air-sea-boundary at h=0
        DO jc = i_startidx_c, i_endidx_c
        IF ( patch_3D%p_patch_1D(1)%dolic_c(jc,jb) > 0 ) THEN
            !This corresponds to full cells
            !ocean_state%p_diag%thick_c(jc,jb) = ocean_state%p_prog(nold(1))%h(jc,jb)&
            !    & + p_patch_3D%p_patch_1D(1)%zlev_i(p_patch_3D%p_patch_1D(1)%dolic_c(jc,jb)+1)
            !  prepare for correct partial cells
            ocean_state%p_diag%thick_c(jc,jb) = ocean_state%p_prog(nold(1))%h(jc,jb) + patch_3D%column_thick_c(jc,jb)
     !      &  - ice_hi(jc,1,jb)
          ENDIF
        END DO
      END DO
    ENDIF!IF ( iswm_oce == 1 )


    !Step 2: calculate edge-located variables for 2D and 3D case from respective cell variables
    !For SWE and for 3D: thick_e = thickness of fluid column at edges
    !         h_e     = surface elevation at edges, without depth of first layer
    IF ( iswm_oce == 1 ) THEN  !  SWM

      DO jb = edges_in_domain%start_block, edges_in_domain%end_block
        CALL get_index_range(edges_in_domain, jb, edge_start_idx, edge_end_idx)
          DO je = edge_start_idx, edge_end_idx

            il_c1 = patch_2D%edges%cell_idx(je,jb,1)
            ib_c1 = patch_2D%edges%cell_blk(je,jb,1)
            il_c2 = patch_2D%edges%cell_idx(je,jb,2)
            ib_c2 = patch_2D%edges%cell_blk(je,jb,2)

            z_dist_e_c1 = 0.5_wp!z_dist_e_c1=p_patch%edges%edge_cell_length(je,jb,1)
            z_dist_e_c2 = 0.5_wp!z_dist_e_c2=p_patch%edges%edge_cell_length(je,jb,2)

            IF(patch_3D%lsm_e(je,1,jb) <= sea_boundary)THEN
              ocean_state%p_diag%thick_e(je,jb) = ( z_dist_e_c1*ocean_state%p_diag%thick_c(il_c1,ib_c1)&
                & +   z_dist_e_c2*ocean_state%p_diag%thick_c(il_c2,ib_c2) )&
                & /(z_dist_e_c1+z_dist_e_c2)

              ocean_state%p_diag%h_e(je,jb) = ( z_dist_e_c1*ocean_state%p_prog(nold(1))%h(il_c1,ib_c1)&
                & +   z_dist_e_c2*ocean_state%p_prog(nold(1))%h(il_c2,ib_c2) )&
                & /(z_dist_e_c1+z_dist_e_c2)
            ELSE
              ocean_state%p_diag%h_e(je,jb) = 0.0_wp
            ENDIF
          END DO
        END DO
        CALL sync_patch_array(SYNC_E, patch_2D, ocean_state%p_diag%thick_e)
        CALL sync_patch_array(SYNC_E, patch_2D, ocean_state%p_diag%h_e)

    ELSEIF(iswm_oce /= 1 )THEN

        DO jb = edges_in_domain%start_block, edges_in_domain%end_block
          CALL get_index_range(edges_in_domain, jb, edge_start_idx, edge_end_idx)
          DO je = edge_start_idx, edge_end_idx

            il_c1 = patch_2D%edges%cell_idx(je,jb,1)
            ib_c1 = patch_2D%edges%cell_blk(je,jb,1)
            il_c2 = patch_2D%edges%cell_idx(je,jb,2)
            ib_c2 = patch_2D%edges%cell_blk(je,jb,2)

            z_dist_e_c1 = 0.5_wp!z_dist_e_c1=p_patch%edges%edge_cell_length(je,jb,1)
            z_dist_e_c2 = 0.5_wp!z_dist_e_c2=p_patch%edges%edge_cell_length(je,jb,2)

            IF ( patch_3D%p_patch_1D(1)%dolic_e(je,jb) > 0 ) THEN

                ocean_state%p_diag%h_e(je,jb) = ( z_dist_e_c1*ocean_state%p_prog(nold(1))%h(il_c1,ib_c1)&
                 & +   z_dist_e_c2*ocean_state%p_prog(nold(1))%h(il_c2,ib_c2) )&
                 & /(z_dist_e_c1+z_dist_e_c2)

            !  sea ice thickness on edges ??
            !   ice_hi_e = ( z_dist_e_c1*ice_hi(il_c1,1,ib_c1)&
            !    &       +   z_dist_e_c2*ice_hi(il_c2,1,ib_c2) )&
            !    & /(z_dist_e_c1+z_dist_e_c2)

                ocean_state%p_diag%thick_e(je,jb) = ocean_state%p_diag%h_e(je,jb) &
                &                          + patch_3D%column_thick_e(je,jb)

            !  or:
            !   ocean_state%p_diag%thick_e(je,jb) = ( z_dist_e_c1*ocean_state%p_diag%thick_c(il_c1,ib_c1)&
            !   & +   z_dist_e_c2*ocean_state%p_diag%thick_c(il_c2,ib_c2) )&
            !   & /(z_dist_e_c1+z_dist_e_c2)

            ENDIF
          END DO
        END DO

        ! CALL sync_patch_array_mult(SYNC_E, patch_2D, 2, ocean_state%p_diag%thick_e, ocean_state%p_diag%h_e)
        CALL sync_patch_array(SYNC_E, patch_2D, ocean_state%p_diag%thick_e)
        CALL sync_patch_array(SYNC_E, patch_2D, ocean_state%p_diag%h_e)
    ENDIF 

    !Update prism thickness. The prism-thickness below the surface is
    !not updated it is initialized in construct_hydro_ocean_diag
    !with z-coordinate-thickness.
    !1) Thickness at cells
    DO jb = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
      DO jc = i_startidx_c, i_endidx_c
        IF(patch_3d%lsm_c(jc,1,jb) <= sea_boundary)THEN
          cell_thickeness(jc,1,jb) = &
            & patch_3d%p_patch_1d(1)%prism_thick_flat_sfc_c(jc,1,jb) + ocean_state%p_prog(nold(1))%h(jc,jb)

          patch_3d%p_patch_1d(1)%prism_center_dist_c(jc,2,jb) = 0.5_wp * &
              & (cell_thickeness(jc,1,jb) + cell_thickeness(jc,2,jb))

          patch_3d%p_patch_1d(1)%prism_volume(jc,1,jb) = cell_thickeness(jc,1,jb) * &
            &  patch_2D%cells%area(jc,jb)

          inv_cell_thickeness(jc,1,jb) = 1.0_wp / cell_thickeness(jc,1,jb)

          inv_prisms_center_distance(jc,2,jb) = &
            1.0_wp / patch_3d%p_patch_1d(1)%prism_center_dist_c(jc,2,jb)

        ELSE
          !Surface thickness over land remains zero
          !ocean_state%p_diag%prism_thick_c(jc,1,jb) = 0.0_wp
          cell_thickeness(jc,1,jb)= 0.0_wp
        ENDIF
      END DO
    END DO

    !2) Thickness at edges
    DO jb = all_edges%start_block, all_edges%end_block
      CALL get_index_range(all_edges, jb, edge_start_idx, edge_end_idx)
      DO je = edge_start_idx, edge_end_idx
        IF(patch_3d%lsm_e(je,1,jb) <= sea_boundary)THEN

          edge_thickeness(je,1,jb)&
            & = patch_3d%p_patch_1d(1)%prism_thick_flat_sfc_e(je,1,jb) + ocean_state%p_diag%h_e(je,jb)

          inv_edge_thickeness(je,1,jb)= 1.0_wp / edge_thickeness(je,1,jb)

          inv_faces_center_distance(je,2,jb) = 2.0_wp / &
              & (edge_thickeness(je,1,jb) + edge_thickeness(je,2,jb))

        ELSE
          !Surfacethickness over land remains zero
          edge_thickeness(je,1,jb)= 0.0_wp
        ENDIF
      END DO
    END DO
    
    !---------------------------------------------------------------------
    ! update the coefficients for the edge2edge_viacell_1D fast operator
    top_coeffs        => operators_coefficients%edge2edge_viacell_coeff_top
    integrated_coeffs => operators_coefficients%edge2edge_viacell_coeff_integrated
    sum_to_2D_coeffs  => operators_coefficients%edge2edge_viacell_coeff_all
    DO jb = all_edges%start_block, all_edges%end_block
      CALL get_index_range(all_edges, jb, edge_start_idx, edge_end_idx)
      DO je = edge_start_idx, edge_end_idx

        IF ( patch_3D%p_patch_1D(1)%dolic_e(je,jb) > 0 ) THEN

          ! get the two cells of the edge
          cell_1_index = patch_2D%edges%cell_idx(je,jb,1)
          cell_1_block = patch_2D%edges%cell_blk(je,jb,1)
          cell_2_index = patch_2D%edges%cell_idx(je,jb,2)
          cell_2_block = patch_2D%edges%cell_blk(je,jb,2)
          cell_thickeness_1 = cell_thickeness(cell_1_index, 1, cell_1_block)
          cell_thickeness_2 = cell_thickeness(cell_2_index, 1, cell_2_block)

          sum_to_2D_coeffs(1, je, jb) = top_coeffs(1, je, jb) * cell_thickeness_1 + integrated_coeffs(1, je, jb)
          sum_to_2D_coeffs(2, je, jb) = top_coeffs(2, je, jb) * cell_thickeness_1 + integrated_coeffs(2, je, jb)
          sum_to_2D_coeffs(3, je, jb) = top_coeffs(3, je, jb) * cell_thickeness_1 + integrated_coeffs(3, je, jb)
          sum_to_2D_coeffs(4, je, jb) = top_coeffs(4, je, jb) * cell_thickeness_2 + integrated_coeffs(4, je, jb)
          sum_to_2D_coeffs(5, je, jb) = top_coeffs(5, je, jb) * cell_thickeness_2 + integrated_coeffs(5, je, jb)
          sum_to_2D_coeffs(6, je, jb) = top_coeffs(6, je, jb) * cell_thickeness_2 + integrated_coeffs(6, je, jb)

        ENDIF
      END DO
    END DO ! jb = edges_in_domain%start_block, edges_in_domain%end_block
    
    !---------Debug Diagnostics-------------------------------------------
    idt_src=4  ! output print level (1-5, fix)
    CALL dbg_print('heightRelQuant: h_e'    ,ocean_state%p_diag%h_e        ,str_module,idt_src, &
      & in_subset=patch_2D%edges%owned)
    idt_src=3
    CALL dbg_print('heightRelQuant: h_c'    ,ocean_state%p_prog(nold(1))%h ,str_module,idt_src, &
      & in_subset=patch_2D%cells%owned)
    CALL dbg_print('heightRelQuant: thick_c',ocean_state%p_diag%thick_c    ,str_module,idt_src, &
      & in_subset=patch_2D%cells%owned)
    CALL dbg_print('heightRelQuant: thick_e',ocean_state%p_diag%thick_e    ,str_module,idt_src, &
      & in_subset=patch_2D%edges%owned)
    !---------------------------------------------------------------------
  END SUBROUTINE calculate_thickness
  !-------------------------------------------------------------------------

  SUBROUTINE check_cfl_horizontal(normal_velocity,inv_dual_edge_length,timestep,edges,threshold, &
      &                          cfl_diag, stop_on_violation, output)
    REAL(wp),POINTER      :: normal_velocity(:,:,:)
    TYPE(t_subset_range)  :: edges
    REAL(wp), INTENT(IN)  :: inv_dual_edge_length(:,:), threshold, timestep
    REAL(wp), POINTER :: cfl_diag(:,:,:)
    LOGICAL, INTENT(IN)   :: stop_on_violation, output
    REAL(wp), POINTER     :: cfl(:,:,:)

    INTEGER :: je,jk,jb,i_startidx,i_endidx

    ALLOCATE(cfl(LBOUND(normal_velocity,1):UBOUND(normal_velocity,1),&
      &          LBOUND(normal_velocity,2):UBOUND(normal_velocity,2),&
      &          LBOUND(normal_velocity,3):UBOUND(normal_velocity,3)))
    cfl = 0.0_wp

    DO jb = edges%start_block, edges%end_block
      CALL get_index_range(edges,jb,i_startidx,i_endidx)
      DO je = i_startidx,i_endidx
        DO jk = 1, n_zlev
          cfl(je,jk,jb) = ABS(dtime*normal_velocity(je,jk,jb)*inv_dual_edge_length(je,jb))
        END DO
      END DO
    END DO
    IF (output) THEN
      IF (ASSOCIATED(cfl_diag)) THEN
        cfl_diag(:,:,:) = cfl(:,:,:)
      ELSE
        CALL finish('check_cfl_vertical','cfl_diag pointer for output NOT ASSOCIATED')
      ENDIF
    ENDIF

    CALL dbg_print('check horiz. CFL',cfl ,str_module,3,in_subset=edges)

    CALL check_cfl_threshold(MAXVAL(cfl),threshold,'horz',stop_on_violation)

    DEALLOCATE(cfl)
  END SUBROUTINE check_cfl_horizontal

  SUBROUTINE check_cfl_vertical(vertical_velocity, thicknesses, timestep, cells, threshold, &
      &                        cfl_diag, stop_on_violation, output)
    REAL(wp),POINTER     :: vertical_velocity(:,:,:), thicknesses(:,:,:)
    REAL(wp), INTENT(IN) :: timestep, threshold
    TYPE(t_subset_range) :: cells
    REAL(wp), POINTER    :: cfl_diag(:,:,:)
    LOGICAL, INTENT(IN)  :: stop_on_violation,output
    REAL(wp), POINTER    :: cfl(:,:,:)

    INTEGER  :: jc, jk, jb, i_startidx_c, i_endidx_c

    ALLOCATE(cfl(LBOUND(vertical_velocity,1):UBOUND(vertical_velocity,1),&
      &          LBOUND(vertical_velocity,2):UBOUND(vertical_velocity,2),&
      &          LBOUND(vertical_velocity,3):UBOUND(vertical_velocity,3)))
    cfl = 0.0_wp

    DO jb = cells%start_block, cells%end_block
      CALL get_index_range(cells, jb, i_startidx_c, i_endidx_c)
      DO jc = i_startidx_c, i_endidx_c
        DO jk=1, cells%vertical_levels(jc,jb)
          cfl(jc,jk,jb)=ABS(dtime*vertical_velocity(jc,jk,jb)/thicknesses(jc,jk,jb))
        END DO
      END DO
    END DO

    IF (output) THEN
      IF (ASSOCIATED(cfl_diag)) THEN
        cfl_diag(:,:,:) = cfl(:,:,:)
      ELSE
        CALL finish('check_cfl_vertical','cfl_diag pointer for output NOT ASSOCIATED')
      ENDIF
    ENDIF

    CALL dbg_print('check vert.  CFL',cfl ,str_module,3,in_subset=cells)

    CALL check_cfl_threshold(MAXVAL(cfl),threshold,'vert',stop_on_violation)

    DEALLOCATE(cfl)
  END SUBROUTINE check_cfl_vertical

  SUBROUTINE check_cfl_threshold(maxcfl,threshold,orientation, stop_on_violation)
    REAL(wp),INTENT(IN)          :: maxcfl, threshold
    CHARACTER(LEN=4), INTENT(IN) :: orientation
    LOGICAL, INTENT(IN)          :: stop_on_violation

    IF (threshold < maxcfl) THEN
      IF (stop_on_violation) THEN
        ! location lookup
        ! location print
        ! throw error
        CALL finish('check_cfl','Found violation of CFL ('//TRIM(orientation)//') criterion')
      ELSE
        CALL message('check_cfl','Found violation of CFL ('//TRIM(orientation)//') criterion')
      END IF
    END IF
  END SUBROUTINE check_cfl_threshold
END MODULE mo_oce_math_operators
