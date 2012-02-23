
#ifdef __xlC__
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
MODULE mo_intp_rbf_coeffs
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
USE mo_exception,           ONLY: message, message_text, finish
USE mo_impl_constants,      ONLY: SUCCESS, min_rlcell_int
USE mo_model_domain,        ONLY: t_patch, t_tangent_vectors
USE mo_grid_config,         ONLY: l_limited_area
USE mo_dynamics_config,     ONLY: iequations
USE mo_math_utilities,      ONLY: gc2cc, gvec2cvec, solve_chol_v, choldec_v, &
  &                               arc_length_v, t_cartesian_coordinates,     &
  &                               t_geographical_coordinates, t_lon_lat_grid
USE mo_parallel_config,     ONLY: nproma
USE mo_loopindices,         ONLY: get_indices_c, get_indices_e, get_indices_v
USE mo_intp_data_strc,      ONLY: t_int_state, t_lon_lat_intp
USE mo_interpol_config,     ONLY: rbf_vec_dim_c, rbf_vec_dim_e, rbf_vec_dim_v,     &
  &                               rbf_c2grad_dim, rbf_vec_kern_c, rbf_vec_kern_e,  &
  &                               rbf_vec_kern_v, rbf_vec_scale_c, rbf_vec_scale_e,&
  &                               rbf_vec_scale_v
USE mo_gnat_gridsearch,     ONLY: gnat_init_grid, gnat_destroy, gnat_tree,&
  &                               gnat_query_containing_triangles,        &
  &                               gnat_merge_distributed_queries, gk
USE mo_math_utilities,      ONLY: rotate_latlon_grid
USE mo_physical_constants,  ONLY: re
USE mo_lonlat_intp_config,  ONLY: lonlat_intp_config
USE mo_mpi,                 ONLY: p_gather_field
USE mo_communication,       ONLY: idx_1d, blk_no, idx_no
USE mo_sync,                ONLY: SYNC_C, SYNC_E, SYNC_V, sync_patch_array, sync_idx

IMPLICIT NONE

!> level of output verbosity
INTEGER, PARAMETER  :: dbg_level = 0

PRIVATE

PUBLIC :: rbf_vec_index_cell, rbf_c2grad_index, rbf_vec_compute_coeff_cell,     &
          & rbf_compute_coeff_c2grad, rbf_vec_index_vertex, rbf_vec_index_edge, &
          & rbf_vec_compute_coeff_vertex, rbf_vec_compute_coeff_edge,           &
          & rbf_setup_interpol_lonlat, rbf_setup_interpol_lonlat_grid

CONTAINS

#include "intp_functions.f90.inc"

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
! Part 1: routines for RBF vector reconstruction and related functions
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!
!>
!! This routine initializes the indexes used to define the stencil.
!!
!! This routine initializes the indexes used to define the stencil
!! of the vector RBF interpolator. The stencil is cell based and includes
!! a variable number of edges (rbf_vec_dim) around each cell, in order
!! to reconstruct a vector at the cell center.
!!
!! @par Revision History
!! Developed and tested  by L. Bonaventura  (2004)
!! Adapted to new data structure by L. Bonaventura and P. Korn, 2006
!! Include 15-point stencil by T. Ruppert, DWD, (2007-02-08)
!! Revision by G. Zaengl, DWD (2009-02-16):
!! - change 15-point stencil to 14 points at pentagon points
!!
SUBROUTINE rbf_vec_index_cell( ptr_patch, ptr_int )
!
!
!  patch on which computation is performed
!
TYPE(t_patch), INTENT(in) :: ptr_patch

TYPE(t_int_state), INTENT(inout) :: ptr_int

INTEGER  :: nblks_c
INTEGER  :: jc, jb           ! loop indices
INTEGER  :: ilc, ibc         ! line and block index of neighbour cell

INTEGER  :: i_startblk       ! start block
INTEGER  :: i_startidx       ! start index
INTEGER  :: i_endidx         ! end index
REAL(wp) :: z_stencil(UBOUND(ptr_int%rbf_vec_stencil_c,1),UBOUND(ptr_int%rbf_vec_stencil_c,2))

!--------------------------------------------------------------------

  ! values for the blocking
  nblks_c  = ptr_patch%nblks_int_c

  !
  ! The stencil consists of 9 edges, the cell edges are taken
  ! and the edges of the nearest neighbours are taken
  !
  ! This stencil is not available when refin_ctrl=1
  i_startblk = ptr_patch%cells%start_blk(2,1)

  DO jb = i_startblk, nblks_c

    CALL get_indices_c(ptr_patch, jb, i_startblk, nblks_c, &
                       i_startidx, i_endidx, 2)

    DO jc = i_startidx, i_endidx

      IF(.NOT. ptr_patch%cells%owner_mask(jc,jb)) CYCLE
      !
      ! get the global line and block indices of the edges of each neighbor
      ! cell and store them in the rbf_vec_idx_c and rbf_vec_blk_c components
      !
      ! Sorry for the clumsy code, but it was not possible to get it vectorized
      ! on the NEC with nested do loops and UNROLL directives
      !
      ! First neighbor cell
      ilc = ptr_patch%cells%neighbor_idx(jc,jb,1)
      ibc = ptr_patch%cells%neighbor_blk(jc,jb,1)
      ptr_int%rbf_vec_idx_c(1,jc,jb) = ptr_patch%cells%edge_idx(ilc,ibc,1)
      ptr_int%rbf_vec_blk_c(1,jc,jb) = ptr_patch%cells%edge_blk(ilc,ibc,1)
      ptr_int%rbf_vec_idx_c(2,jc,jb) = ptr_patch%cells%edge_idx(ilc,ibc,2)
      ptr_int%rbf_vec_blk_c(2,jc,jb) = ptr_patch%cells%edge_blk(ilc,ibc,2)
      ptr_int%rbf_vec_idx_c(3,jc,jb) = ptr_patch%cells%edge_idx(ilc,ibc,3)
      ptr_int%rbf_vec_blk_c(3,jc,jb) = ptr_patch%cells%edge_blk(ilc,ibc,3)

      ! Second neighbor cell
      ilc = ptr_patch%cells%neighbor_idx(jc,jb,2)
      ibc = ptr_patch%cells%neighbor_blk(jc,jb,2)
      ptr_int%rbf_vec_idx_c(4,jc,jb) = ptr_patch%cells%edge_idx(ilc,ibc,1)
      ptr_int%rbf_vec_blk_c(4,jc,jb) = ptr_patch%cells%edge_blk(ilc,ibc,1)
      ptr_int%rbf_vec_idx_c(5,jc,jb) = ptr_patch%cells%edge_idx(ilc,ibc,2)
      ptr_int%rbf_vec_blk_c(5,jc,jb) = ptr_patch%cells%edge_blk(ilc,ibc,2)
      ptr_int%rbf_vec_idx_c(6,jc,jb) = ptr_patch%cells%edge_idx(ilc,ibc,3)
      ptr_int%rbf_vec_blk_c(6,jc,jb) = ptr_patch%cells%edge_blk(ilc,ibc,3)

      ! Third neighbor cell
      ilc = ptr_patch%cells%neighbor_idx(jc,jb,3)
      ibc = ptr_patch%cells%neighbor_blk(jc,jb,3)
      ptr_int%rbf_vec_idx_c(7,jc,jb) = ptr_patch%cells%edge_idx(ilc,ibc,1)
      ptr_int%rbf_vec_blk_c(7,jc,jb) = ptr_patch%cells%edge_blk(ilc,ibc,1)
      ptr_int%rbf_vec_idx_c(8,jc,jb) = ptr_patch%cells%edge_idx(ilc,ibc,2)
      ptr_int%rbf_vec_blk_c(8,jc,jb) = ptr_patch%cells%edge_blk(ilc,ibc,2)
      ptr_int%rbf_vec_idx_c(9,jc,jb) = ptr_patch%cells%edge_idx(ilc,ibc,3)
      ptr_int%rbf_vec_blk_c(9,jc,jb) = ptr_patch%cells%edge_blk(ilc,ibc,3)

      ! take care of cells at patch boundaries, then the value of 
      ! ptr_int%rbf_vec_stencil_c might be smaller than rbf_vec_dim_c:
      ptr_int%rbf_vec_stencil_c(jc,jb) = &
        & COUNT(ptr_int%rbf_vec_idx_c(:,jc,jb) /= 0)

    END DO

  END DO

  DO jb = 1, rbf_vec_dim_c
    CALL sync_idx(SYNC_C, SYNC_E, ptr_patch, ptr_int%rbf_vec_idx_c(jb,:,:), &
                                           & ptr_int%rbf_vec_blk_c(jb,:,:))
  ENDDO

  z_stencil(:,:) = ptr_int%rbf_vec_stencil_c(:,:)
  CALL sync_patch_array(SYNC_C,ptr_patch,z_stencil)
  ptr_int%rbf_vec_stencil_c(:,:) = z_stencil(:,:)

END SUBROUTINE rbf_vec_index_cell

!>
!! This routine creates the index list needed for RBF reconstruction
!! of gradients of scalar values at cell centers, using the values at the
!! surrounding cell centers.
!!
!! @par Revision History
!! Developed by G. Zaengl, DWD (2009-12-15)
!!
SUBROUTINE rbf_c2grad_index( ptr_patch, ptr_int )
!
!
!  patch on which computation is performed
!
TYPE(t_patch), TARGET, INTENT(in) :: ptr_patch

TYPE(t_int_state), INTENT(inout) :: ptr_int

INTEGER  :: nblks_c
INTEGER  :: jc, jb           ! loop indices
INTEGER  :: ilc, ibc         ! line and block index of neighbour cell
INTEGER  :: ic1, ic2         ! index counters

INTEGER  :: i_startblk       ! start block
INTEGER  :: i_startidx       ! start index
INTEGER  :: i_endidx         ! end index

! Pointers to neighbor indices/blocks
INTEGER, DIMENSION(:,:,:), POINTER :: inidx, inblk

!--------------------------------------------------------------------

  ! values for the blocking
  nblks_c  = ptr_patch%nblks_int_c

  inidx => ptr_patch%cells%neighbor_idx
  inblk => ptr_patch%cells%neighbor_blk
  !
  ! The stencil has a total of 10 cells, namely the local cell,
  ! its neighbors, and the neighbors of the neighbors (excluding the local cell)
  !
  ! This stencil is not available when refin_ctrl=1
  i_startblk = ptr_patch%cells%start_blk(2,1)

  DO jb = i_startblk, nblks_c

    CALL get_indices_c(ptr_patch, jb, i_startblk, nblks_c, &
                       i_startidx, i_endidx, 2)

    DO jc = i_startidx, i_endidx

      IF(.NOT. ptr_patch%cells%owner_mask(jc,jb)) CYCLE

      ! Local point
      ptr_int%rbf_c2grad_idx(1,jc,jb) = jc
      ptr_int%rbf_c2grad_blk(1,jc,jb) = jb

      ! Direct neighbors
      ptr_int%rbf_c2grad_idx(2,jc,jb) = inidx(jc,jb,1)
      ptr_int%rbf_c2grad_blk(2,jc,jb) = inblk(jc,jb,1)
      ptr_int%rbf_c2grad_idx(5,jc,jb) = inidx(jc,jb,2)
      ptr_int%rbf_c2grad_blk(5,jc,jb) = inblk(jc,jb,2)
      ptr_int%rbf_c2grad_idx(8,jc,jb) = inidx(jc,jb,3)
      ptr_int%rbf_c2grad_blk(8,jc,jb) = inblk(jc,jb,3)

      ! Neighbors of neighbor 1 (excluding local cell)
      ilc = inidx(jc,jb,1)
      ibc = inblk(jc,jb,1)
      IF (inidx(ilc,ibc,1) == jc .AND. inblk(ilc,ibc,1) == jb) THEN
        ic1 = 2
        ic2 = 3
      ELSE IF (inidx(ilc,ibc,2) == jc .AND. inblk(ilc,ibc,2) == jb) THEN
        ic1 = 1
        ic2 = 3
      ELSE
        ic1 = 1
        ic2 = 2
      ENDIF
      ptr_int%rbf_c2grad_idx(3,jc,jb) = inidx(ilc,ibc,ic1)
      ptr_int%rbf_c2grad_blk(3,jc,jb) = inblk(ilc,ibc,ic1)
      ptr_int%rbf_c2grad_idx(4,jc,jb) = inidx(ilc,ibc,ic2)
      ptr_int%rbf_c2grad_blk(4,jc,jb) = inblk(ilc,ibc,ic2)

      ! Neighbors of neighbor 2 (excluding local cell)
      ilc = inidx(jc,jb,2)
      ibc = inblk(jc,jb,2)
      IF (inidx(ilc,ibc,1) == jc .AND. inblk(ilc,ibc,1) == jb) THEN
        ic1 = 2
        ic2 = 3
      ELSE IF (inidx(ilc,ibc,2) == jc .AND. inblk(ilc,ibc,2) == jb) THEN
        ic1 = 1
        ic2 = 3
      ELSE
        ic1 = 1
        ic2 = 2
      ENDIF
      ptr_int%rbf_c2grad_idx(6,jc,jb) = inidx(ilc,ibc,ic1)
      ptr_int%rbf_c2grad_blk(6,jc,jb) = inblk(ilc,ibc,ic1)
      ptr_int%rbf_c2grad_idx(7,jc,jb) = inidx(ilc,ibc,ic2)
      ptr_int%rbf_c2grad_blk(7,jc,jb) = inblk(ilc,ibc,ic2)

      ! Neighbors of neighbor 3 (excluding local cell)
      ilc = inidx(jc,jb,3)
      ibc = inblk(jc,jb,3)
      IF (inidx(ilc,ibc,1) == jc .AND. inblk(ilc,ibc,1) == jb) THEN
        ic1 = 2
        ic2 = 3
      ELSE IF (inidx(ilc,ibc,2) == jc .AND. inblk(ilc,ibc,2) == jb) THEN
        ic1 = 1
        ic2 = 3
      ELSE
        ic1 = 1
        ic2 = 2
      ENDIF
      ptr_int%rbf_c2grad_idx(9,jc,jb) = inidx(ilc,ibc,ic1)
      ptr_int%rbf_c2grad_blk(9,jc,jb) = inblk(ilc,ibc,ic1)
      ptr_int%rbf_c2grad_idx(10,jc,jb) = inidx(ilc,ibc,ic2)
      ptr_int%rbf_c2grad_blk(10,jc,jb) = inblk(ilc,ibc,ic2)
    END DO
  END DO

  DO jb = 1, rbf_c2grad_dim
    CALL sync_idx(SYNC_C, SYNC_C, ptr_patch, ptr_int%rbf_c2grad_idx(jb,:,:), &
                                           & ptr_int%rbf_c2grad_blk(jb,:,:))
  ENDDO

END SUBROUTINE rbf_c2grad_index

!>
!! This routine initializes the indexes used to define the stencil.
!!
!! This routine initializes the indexes used to define the stencil
!! of the vector RBF interpolator e.g. for the Simpson rule. The stencil
!! is vertex based and includes a variable number of edges
!! (rbf_vec_dim_v) around each vertex, in order
!! to reconstruct a vector at the vertices of a triangle.
!!
!! @par Revision History
!! Developed and tested  by J. Foerstner (April 2008)
!!
SUBROUTINE rbf_vec_index_vertex( ptr_patch, ptr_int )

!
TYPE(t_patch), INTENT(in) :: ptr_patch

TYPE(t_int_state), INTENT(inout) :: ptr_int

INTEGER,  ALLOCATABLE :: ile(:), ibe(:) ! edge indices for the stencil

INTEGER :: nblks_v
INTEGER :: istencil                 ! number of edges for the stencil
INTEGER :: jv, jb                   ! loop indices
INTEGER :: jje                      ! additional loop indices
INTEGER :: iie, ipent
INTEGER  :: i_startblk       ! start block
INTEGER  :: i_startidx       ! start index
INTEGER  :: i_endidx         ! end index
LOGICAL :: ll_pent                  ! if .TRUE. vertex is center of a pentagon
INTEGER :: ist                      ! status variable
INTEGER :: jg

REAL(wp) :: z_stencil(UBOUND(ptr_int%rbf_vec_stencil_v,1),UBOUND(ptr_int%rbf_vec_stencil_v,2))

!--------------------------------------------------------------------

  ! values for the blocking
  nblks_v  = ptr_patch%nblks_int_v

  ! init pentagon counter for security check
  ipent = 0

  !
  ! The stencil consists of 6 edges, the edges sharing the vertex are taken.
  ! (for the central vertex of a pentagon we use only 5 edges)
  !
  ! This stencil is not available when refin_ctrl=1
  i_startblk = ptr_patch%verts%start_blk(2,1)

!$OMP PARALLEL PRIVATE(ile,ibe,ist)
  ALLOCATE( ile(rbf_vec_dim_v), ibe(rbf_vec_dim_v), STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:rbf_vec_index_vertex',  &
      &             'allocation for ile, ibe failed')
  ENDIF
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jv,jje,istencil,ll_pent,iie)
  DO jb = i_startblk, nblks_v

    CALL get_indices_v(ptr_patch, jb, i_startblk, nblks_v, &
                       i_startidx, i_endidx, 2)

    DO jv = i_startidx, i_endidx

      IF(.NOT. ptr_patch%verts%owner_mask(jv,jb)) CYCLE

      istencil = 6
      ll_pent = .FALSE.

      ! init counter
      iie = 0

      ! check if vertex is central vertex of a pentagon
      IF ( ptr_patch%verts%num_edges(jv,jb) == 5 ) THEN
        istencil = 5
        ll_pent = .TRUE.
!$OMP ATOMIC
        ipent = ipent+1
      END IF

      DO jje = 1, ptr_patch%verts%num_edges(jv,jb)
        !
        ! get global indices of the edges around the vertex
        !
        iie = iie+1
        ile(iie) = ptr_patch%verts%edge_idx(jv,jb,jje)
        ibe(iie) = ptr_patch%verts%edge_blk(jv,jb,jje)
      END DO

      IF (iie /= istencil) THEN
        PRINT *, "ERROR ==>  iie = ", iie,  " /= istencil = ", istencil
        CALL finish ('mo_interpolation:rbf_vec_index_vertex',  &
          &             'wrong number of stencil points')
      ENDIF

      ! save number of edges for the stencil in rbf_vec_stencil_v
      ptr_int%rbf_vec_stencil_v(jv,jb) = istencil

      DO jje = 1, istencil
        !
        ! store line and block indices of edges for the stencil of vertex jv
        ! in rbf_vec_idx_v, rbf_vec_blk_v
        !
        ptr_int%rbf_vec_idx_v(jje,jv,jb) = ile(jje)
        ptr_int%rbf_vec_blk_v(jje,jv,jb) = ibe(jje)
      END DO

      ! for pentagons fill index arrays with dummy values
      IF ( ll_pent ) THEN
        DO jje = istencil+1, rbf_vec_dim_v
          ptr_int%rbf_vec_idx_v(jje,jv,jb) = ile(jje-istencil)
          ptr_int%rbf_vec_blk_v(jje,jv,jb) = ibe(jje-istencil)
        END DO
      END IF

    END DO ! end vertex loop

  END DO ! end block loop
!$OMP END DO

  DEALLOCATE( ile, ibe, STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:rbf_vec_index_vertex',  &
      &             'deallocation for ile, ibe failed')
  ENDIF
!$OMP END PARALLEL

  DO jb = 1, rbf_vec_dim_v
    CALL sync_idx(SYNC_V, SYNC_E, ptr_patch, ptr_int%rbf_vec_idx_v(jb,:,:), &
                                           & ptr_int%rbf_vec_blk_v(jb,:,:))
  ENDDO

  z_stencil(:,:) = ptr_int%rbf_vec_stencil_v(:,:)
  CALL sync_patch_array(SYNC_V,ptr_patch,z_stencil)
  ptr_int%rbf_vec_stencil_v(:,:) = z_stencil(:,:)

  jg = ptr_patch%id
  IF ((.NOT. l_limited_area) .AND. (jg == 1) .AND. (ipent /= 12)) THEN
!   PRINT *, "ERROR ==>  ipent = ", ipent,  " /= 12 -- no finish (ocean)"
    IF ( iequations == -1 )  THEN    !  hydrostatic ocean: less than  12 pentagons allowed
      WRITE(message_text,'(a,i3,a,i3,a)')  &
        &  'mo_interpolation:rbf_vec_index_vertex: no of pentagons =',ipent, &
        &  ' iequations =',iequations,' ==> ok.'
      CALL message('', TRIM(message_text))
    ELSE                      !  other cases: may stop the run
!Please note: This is the normal case for parallel runs!
      CALL message('mo_interpolation:rbf_vec_index_vertex',  &
        &             'wrong number of detected pentagons')
!     CALL finish ('mo_interpolation:rbf_vec_index_vertex',  &
!       &             'wrong number of detected pentagons')
    ENDIF
  ENDIF

END SUBROUTINE rbf_vec_index_vertex

!-------------------------------------------------------------------------
!
!
!>
!! This routine initializes the indexes used to define the stencil.
!!
!! This routine initializes the indexes used to define the stencil
!! of the vector RBF interpolator. The stencil is edge based and
!! includes a variable number of edges (rbf_vec_dim_e) around
!! each edge, in order to reconstruct a vector at the edge midpoints.
!!
!! @par Revision History
!! Developed and tested  by J. Foerstner (2008-07-15)
!! Modification by G. Zaengl (2009-02-13):
!! Exclude local point from all stencils for direct reconstruction
!! of tangential velocity component
!! Modification by G. Zaengl (2009-02-17):
!! Add 12-point stencil for higher-order reconstruction of tangential
!! velocity component
!!
SUBROUTINE rbf_vec_index_edge(ptr_patch, ptr_int)
!
TYPE(t_patch), INTENT(in) :: ptr_patch

TYPE(t_int_state), INTENT(inout) :: ptr_int

INTEGER :: nblks_e
INTEGER :: istencil                 ! number of edges for the stencil
INTEGER :: je, jb                   ! loop indices

INTEGER  :: i_startblk       ! start block
INTEGER  :: i_startidx       ! start index
INTEGER  :: i_endidx         ! end index

REAL(wp) :: z_stencil(UBOUND(ptr_int%rbf_vec_stencil_e,1),UBOUND(ptr_int%rbf_vec_stencil_e,2))

!--------------------------------------------------------------------

  ! values for the blocking
  nblks_e  = ptr_patch%nblks_int_e

  !
  ! The stencil consists of 4 edges, the edges of two neighboring
  ! triangles, excluding the local point, are taken.
  !
  ! This stencil is not available when refin_ctrl=1
  i_startblk = ptr_patch%edges%start_blk(2,1)

  istencil = 4

  DO jb = i_startblk, nblks_e

    CALL get_indices_e(ptr_patch, jb, i_startblk, nblks_e, &
                       i_startidx, i_endidx, 2)

    DO je = i_startidx, i_endidx

      ! There is not much work to do because the required stencil points
      ! are already stored in quad_idx/quad_blk
      ptr_int%rbf_vec_idx_e(1,je,jb) = ptr_patch%edges%quad_idx(je,jb,1)
      ptr_int%rbf_vec_blk_e(1,je,jb) = ptr_patch%edges%quad_blk(je,jb,1)
      ptr_int%rbf_vec_idx_e(2,je,jb) = ptr_patch%edges%quad_idx(je,jb,2)
      ptr_int%rbf_vec_blk_e(2,je,jb) = ptr_patch%edges%quad_blk(je,jb,2)
      ptr_int%rbf_vec_idx_e(3,je,jb) = ptr_patch%edges%quad_idx(je,jb,3)
      ptr_int%rbf_vec_blk_e(3,je,jb) = ptr_patch%edges%quad_blk(je,jb,3)
      ptr_int%rbf_vec_idx_e(4,je,jb) = ptr_patch%edges%quad_idx(je,jb,4)
      ptr_int%rbf_vec_blk_e(4,je,jb) = ptr_patch%edges%quad_blk(je,jb,4)

      ! save number of edges for the stencil in rbf_vec_stencil_e
      ptr_int%rbf_vec_stencil_e(je,jb) = istencil

    END DO ! end edge loop

  END DO ! end block loop

  DO jb = 1, rbf_vec_dim_e
    CALL sync_idx(SYNC_E, SYNC_E, ptr_patch, ptr_int%rbf_vec_idx_e(jb,:,:), &
                                           & ptr_int%rbf_vec_blk_e(jb,:,:))
  ENDDO

  ! Not really necessary, only for the case that rbf_vec_stencil_e should be changed:
  z_stencil(:,:) = ptr_int%rbf_vec_stencil_e(:,:)
  CALL sync_patch_array(SYNC_E,ptr_patch,z_stencil)
  ptr_int%rbf_vec_stencil_e(:,:) = z_stencil(:,:)

END SUBROUTINE rbf_vec_index_edge

!-------------------------------------------------------------------------
!
!
!>
!! This routine computes the coefficients needed for vector RBF interpolation,.
!!
!! This routine computes the coefficients needed for vector RBF interpolation,
!! which are then stored in the arrays <i>rbf_vec_coeff</i>.
!! This computation involves the inversion of the interpolation
!! matrix, which is performed by a Cholesky decomposition.
!! The Cholesky decomposition is currently implemented by a home made routine
!! which can be substituted by a call to a numerical library, if available.
!!
!! @par Revision History
!! Developed and tested by Will Sawyer, ETHZ, adjusted and tested  for
!! vector rbf by Luca Bonaventura (2005)
!! Adapted to new data structure by L.Bonaventura and P.Korn (2006).
!! Modifications by Tobias Ruppert and Thomas Heinze, DWD (2006-11-14):
!! - LU decomposition replaced by Cholesky decomposition
!! - only lower triangular matrix is calculated, upper is copied from lower
!! Modifications by Tobias Ruppert and Thomas Heinze, DWD (2006-12-05):
!! - distances calculated by arc_length
!! - grid points used in cartesian coordinate system
!! Modifications by Tobias Ruppert, DWD (2007-02-08):
!! - included thin plate splines (rbf_vec_kern_c == 4)
!! Modifications by Almut Gassmann, MPI-M (2007-04-30)
!! - abandon grid for the sake of patch
!! Modification by Guenther Zaengl, DWD (2009-02-13)
!! - change to direct reconstruction of vector components
!! Modification by Guenther Zaengl, DWD (2009-04-20)
!! - vector optimization and removal of unused options (polynomial
!!   component in RBF kernel, multiquadric and thin-plate-spline kernel)
!!
SUBROUTINE rbf_vec_compute_coeff_cell( ptr_patch, ptr_int )
!

!
TYPE(t_patch), INTENT(in) :: ptr_patch

TYPE(t_int_state), INTENT(inout) :: ptr_int

REAL(wp) :: cc_e1(3), cc_e2(3), cc_c(nproma,3)  ! coordinates of edge midpoints

TYPE(t_cartesian_coordinates) :: cc_center    ! coordinates of cell centers

REAL(wp) :: z_lon, z_lat          ! longitude and latitude

! 3d  normal velocity vectors at edge midpoints
REAL(wp), DIMENSION (nproma,3) :: z_nx1, z_nx2, z_nx3

REAL(wp) :: z_norm                ! norm of velocity vectors

REAL(wp) :: z_dist                ! distance between data points

REAL(wp) :: z_nxprod              ! scalar product of normal
                                  ! velocity vectors

REAL(wp),ALLOCATABLE :: z_rbfmat(:,:,:)    ! RBF interpolation matrix

REAL(wp),ALLOCATABLE :: z_diag(:,:)        ! diagonal of cholesky
                                           ! decomposition matrix

REAL(wp),ALLOCATABLE :: z_rbfval(:,:)      ! RBF function value

REAL(wp),ALLOCATABLE :: z_rhs1(:,:), &     ! right hand side of linear
                      & z_rhs2(:,:)        ! interpolation system in 2d

INTEGER :: nblks_c
INTEGER :: jc, jb          ! integer over cells, blocks and levels
INTEGER :: je1, je2        ! integer over edges
INTEGER :: ile1, ibe1, ile2, ibe2  ! edge indices
INTEGER :: ist             ! return value of array allocation
INTEGER :: i_startblk      ! start block
INTEGER :: i_startidx      ! start index
INTEGER :: i_endidx        ! end index
INTEGER :: i_rcstartlev    ! refinement control start level
INTEGER :: istencil(nproma) ! actual number of stencil points
INTEGER :: jg


REAL(wp) ::  checksum_u,checksum_v ! to check if sum of interpolation coefficients is correct

!--------------------------------------------------------------------

  CALL message('mo_interpolation:rbf_vec_compute_coeff_cell', '')

  jg = ptr_patch%id

  i_rcstartlev = 2

  ! values for the blocking
  nblks_c  = ptr_patch%nblks_int_c

  ! The start block depends on the width of the stencil
  i_startblk = ptr_patch%cells%start_blk(i_rcstartlev,1)

!$OMP PARALLEL PRIVATE (z_rbfmat,z_diag,z_rbfval,z_rhs1,z_rhs2, ist)
  ALLOCATE( z_rbfmat(nproma,rbf_vec_dim_c,rbf_vec_dim_c),  &
            z_diag(nproma,rbf_vec_dim_c),                  &
            z_rbfval(nproma,rbf_vec_dim_c),                &
            z_rhs1(nproma,rbf_vec_dim_c),                  &
            z_rhs2(nproma,rbf_vec_dim_c), STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:rbf_vec_compute_coeff_cell',  &
      &             'allocation for working arrays failed')
  ENDIF

!$OMP DO PRIVATE (jb,jc,i_startidx,i_endidx,je1,je2,istencil,      &
!$OMP             ist,ile1,ibe1,cc_e1,z_lon,z_lat,z_norm,      &
!$OMP             z_nx1,ile2,ibe2,cc_e2,cc_c,z_nx2,z_nxprod,z_dist,      &
!$OMP             cc_center,z_nx3,checksum_u,checksum_v)
  DO jb = i_startblk, nblks_c

    CALL get_indices_c(ptr_patch, jb, i_startblk, nblks_c, &
                       i_startidx, i_endidx, i_rcstartlev)

    !
    ! for each cell, build the vector RBF interpolation matrix
    !
    DO je1 = 1, rbf_vec_dim_c

      DO je2 = 1, je1

        DO jc = i_startidx, i_endidx

          IF(.NOT. ptr_patch%cells%owner_mask(jc,jb)) THEN
            ! Avoid the matrix decomposition for boundary cells since the matrix might get singular
            istencil(jc) = 0
            CYCLE
          ENDIF

          ! Get actual number of stencil points
          istencil(jc) = ptr_int%rbf_vec_stencil_c(jc,jb)
          !
          IF ( (je1 > istencil(jc)) .OR. (je2 > istencil(jc)) ) CYCLE
          !
          ! line and block indices for each edge je1 and je2 of RBF stencil
          !
          ile1 = ptr_int%rbf_vec_idx_c(je1,jc,jb)
          ibe1 = ptr_int%rbf_vec_blk_c(je1,jc,jb)
          ile2 = ptr_int%rbf_vec_idx_c(je2,jc,jb)
          ibe2 = ptr_int%rbf_vec_blk_c(je2,jc,jb)
          !
          ! get Cartesian coordinates and orientation vectors
          !
          cc_e1(:) = ptr_int%cart_edge_coord(ile1,ibe1,:)
          cc_e2(:) = ptr_int%cart_edge_coord(ile2,ibe2,:)
          !
          z_nx1(jc,:) = ptr_patch%edges%primal_cart_normal(ile1,ibe1)%x(:)
          z_nx2(jc,:) = ptr_patch%edges%primal_cart_normal(ile2,ibe2)%x(:)
          !
          ! compute dot product of normal vectors and distance between edge midpoints
          !
          z_nxprod = DOT_PRODUCT(z_nx1(jc,:),z_nx2(jc,:))
          z_dist   = arc_length_v(cc_e1,cc_e2)
          !
          ! set up interpolation matrix
          !
          IF      (rbf_vec_kern_c == 1) THEN
            z_rbfmat(jc,je1,je2) = z_nxprod * gaussi(z_dist,rbf_vec_scale_c(MAX(jg,1)))
          ELSE IF (rbf_vec_kern_c == 3) THEN
            z_rbfmat(jc,je1,je2) = z_nxprod * inv_multiq(z_dist,rbf_vec_scale_c(MAX(jg,1)))
          ENDIF

          IF (je1 > je2) z_rbfmat(jc,je2,je1) = z_rbfmat(jc,je1,je2)
        END DO

      END DO

    END DO

    ! apply Cholesky decomposition to matrix
    !
!CDIR NOIEXPAND
#ifdef __SX__
    CALL choldec_v(i_startidx,i_endidx,istencil,rbf_vec_dim_c,z_rbfmat,z_diag)
#else
    CALL choldec_v(i_startidx,i_endidx,istencil,              z_rbfmat,z_diag)
#endif

    DO jc = i_startidx, i_endidx

      IF(.NOT. ptr_patch%cells%owner_mask(jc,jb)) CYCLE
      !
      ! Solve immediately for coefficients
      !
      ! convert coordinates of cell center to cartesian vector
      !
      cc_center = gc2cc(ptr_patch%cells%center(jc,jb))
      cc_c(jc,1:3) = cc_center%x(1:3)

      z_lon = ptr_patch%cells%center(jc,jb)%lon
      z_lat  = ptr_patch%cells%center(jc,jb)%lat

      ! Zonal wind component
      CALL gvec2cvec(1._wp,0._wp,z_lon,z_lat,z_nx1(jc,1),z_nx1(jc,2),z_nx1(jc,3))

      z_norm = SQRT( DOT_PRODUCT(z_nx1(jc,:),z_nx1(jc,:)) )
      z_nx1(jc,:)  = 1._wp/z_norm * z_nx1(jc,:)

      ! Meridional wind component
      CALL gvec2cvec(0._wp,1._wp,z_lon,z_lat,z_nx2(jc,1),z_nx2(jc,2),z_nx2(jc,3))

      z_norm = SQRT( DOT_PRODUCT(z_nx2(jc,:),z_nx2(jc,:)) )
      z_nx2(jc,:)  = 1._wp/z_norm * z_nx2(jc,:)

    END DO

    !
    ! set up right hand side for interpolation system
    !
    DO je2 = 1, rbf_vec_dim_c

      DO jc = i_startidx, i_endidx

        IF(.NOT. ptr_patch%cells%owner_mask(jc,jb)) CYCLE

        IF (je2 > istencil(jc)) CYCLE
        !
        ! get indices and coordinates of edge midpoints and compute distance
        ! to cell center
        !
        ile2   = ptr_int%rbf_vec_idx_c(je2,jc,jb)
        ibe2   = ptr_int%rbf_vec_blk_c(je2,jc,jb)
        !
        cc_e2(:)  = ptr_int%cart_edge_coord(ile2,ibe2,:)

        z_dist = arc_length_v(cc_c(jc,:), cc_e2)

        !
        ! get Cartesian orientation vector
        z_nx3(jc,:) = ptr_patch%edges%primal_cart_normal(ile2,ibe2)%x(:)

        IF (rbf_vec_kern_c == 1) THEN
          z_rbfval(jc,je2) = gaussi(z_dist,rbf_vec_scale_c(MAX(jg,1)))
        ELSE IF (rbf_vec_kern_c == 3) THEN
          z_rbfval(jc,je2) = inv_multiq(z_dist,rbf_vec_scale_c(MAX(jg,1)))
        ENDIF
        !
        ! compute projection on target vector orientation
        !
        z_rhs1(jc,je2) = z_rbfval(jc,je2) * DOT_PRODUCT(z_nx1(jc,:),z_nx3(jc,:))
        z_rhs2(jc,je2) = z_rbfval(jc,je2) * DOT_PRODUCT(z_nx2(jc,:),z_nx3(jc,:))

      END DO

    END DO

    !
    ! compute vector coefficients
    !
!CDIR NOIEXPAND
#ifdef __SX__
    CALL solve_chol_v(i_startidx, i_endidx, istencil, rbf_vec_dim_c, z_rbfmat,  &
      &               z_diag, z_rhs1, ptr_int%rbf_vec_coeff_c(:,1,:,jb))
#else
    CALL solve_chol_v(i_startidx, i_endidx, istencil,                z_rbfmat,  &
      &               z_diag, z_rhs1, ptr_int%rbf_vec_coeff_c(:,1,:,jb))
#endif
!CDIR NOIEXPAND
#ifdef __SX__
    CALL solve_chol_v(i_startidx, i_endidx, istencil, rbf_vec_dim_c, z_rbfmat,  &
      &               z_diag, z_rhs2, ptr_int%rbf_vec_coeff_c(:,2,:,jb))
#else
    CALL solve_chol_v(i_startidx, i_endidx, istencil,                z_rbfmat,  &
      &               z_diag, z_rhs2, ptr_int%rbf_vec_coeff_c(:,2,:,jb))
#endif

    DO jc = i_startidx, i_endidx

      IF(.NOT. ptr_patch%cells%owner_mask(jc,jb)) CYCLE

      ! Ensure that sum of interpolation coefficients is correct

      checksum_u = 0._wp
      checksum_v = 0._wp

      DO je1 = 1, istencil(jc)
        ile1   = ptr_int%rbf_vec_idx_c(je1,jc,jb)
        ibe1   = ptr_int%rbf_vec_blk_c(je1,jc,jb)

        ! get Cartesian orientation vector
        z_nx3(jc,:) = ptr_patch%edges%primal_cart_normal(ile1,ibe1)%x(:)

        checksum_u = checksum_u + ptr_int%rbf_vec_coeff_c(je1,1,jc,jb)* &
          DOT_PRODUCT(z_nx1(jc,:),z_nx3(jc,:))
        checksum_v = checksum_v + ptr_int%rbf_vec_coeff_c(je1,2,jc,jb)* &
          DOT_PRODUCT(z_nx2(jc,:),z_nx3(jc,:))
      ENDDO

      DO je1 = 1, istencil(jc)
        ptr_int%rbf_vec_coeff_c(je1,1,jc,jb) = &
          ptr_int%rbf_vec_coeff_c(je1,1,jc,jb) / checksum_u
        ptr_int%rbf_vec_coeff_c(je1,2,jc,jb) = &
          ptr_int%rbf_vec_coeff_c(je1,2,jc,jb) / checksum_v
      ENDDO

    END DO

  END DO
!$OMP END DO

  DEALLOCATE( z_rbfmat, z_diag, z_rbfval, z_rhs1, z_rhs2,  &
    STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:rbf_vec_compute_coeff_cell',  &
      &             'deallocation for working arrays failed')
  ENDIF
!$OMP END PARALLEL

  DO jb = 1, rbf_vec_dim_c
    call sync_patch_array(SYNC_C, ptr_patch, ptr_int%rbf_vec_coeff_c(jb,1,:,:))
    call sync_patch_array(SYNC_C, ptr_patch, ptr_int%rbf_vec_coeff_c(jb,2,:,:))
  ENDDO


! Optional debug output for RBF coefficients
#ifdef DEBUG_COEFF
  DO jb = i_startblk, nblks_c

    CALL get_indices_c(ptr_patch, jb, i_startblk, nblks_c, &
                       i_startidx, i_endidx, i_rcstartlev)

      DO jc = i_startidx, i_endidx
      istencil(jc) = ptr_int%rbf_vec_stencil_c(jc,jb)
      WRITE(510+jg,'(2i5,25f12.6)') jb,jc,ptr_int%rbf_vec_coeff_c(1:istencil(jc),1,jc,jb)
      WRITE(510+jg,'(2i5,25f12.6)') jb,jc,ptr_int%rbf_vec_coeff_c(1:istencil(jc),2,jc,jb)
    END DO
  END DO
  CLOSE (510+jg)
#endif

END SUBROUTINE rbf_vec_compute_coeff_cell


!>
!! This routine computes the coefficients needed for reconstructing the gradient
!! of a cell-based variable at the cell center. The operations performed here
!! combine taking the centered-difference gradient (like in grad_fd_norm) at
!! the 9 edges entering into the usual 9-point RBF vector reconstruction at
!! cell centers, followed by applying the RBF reconstruction.
!!
!! @par Revision History
!! Developed and tested by Guenther Zaengl, DWD (2009-12-15)
!!
SUBROUTINE rbf_compute_coeff_c2grad (ptr_patch, ptr_int)

!
!  patch on which computation is performed
!
TYPE(t_patch), INTENT(in) :: ptr_patch

! interpolation state
!
TYPE(t_int_state), INTENT(inout) :: ptr_int

INTEGER :: jc, jb, je, jcc
INTEGER :: i_rcstartlev
INTEGER :: i_startblk, nblks_c, i_startidx, i_endidx

INTEGER :: ile, ibe, ilc1, ibc1, ilc2, ibc2, ilcc, ibcc

REAL(wp), DIMENSION(nproma,rbf_c2grad_dim,2) :: aux_coeff

!--------------------------------------------------------------------

  i_rcstartlev = 2

  ! Values for the blocking
  i_startblk = ptr_patch%cells%start_blk(i_rcstartlev,1)
  nblks_c    = ptr_patch%nblks_int_c

  ! loop through all patch cells (and blocks)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,je,jcc,jc,i_startidx,i_endidx,ile,ibe,ilc1,ibc1,&
!$OMP    ilc2,ibc2,ilcc,ibcc,aux_coeff)
  DO jb = i_startblk, nblks_c

    CALL get_indices_c(ptr_patch, jb, i_startblk, nblks_c, &
                       i_startidx, i_endidx, i_rcstartlev)

    aux_coeff(:,:,:) = 0._wp

    DO je = 1, rbf_vec_dim_c
      DO jcc = 1, rbf_c2grad_dim
!CDIR NODEP
        DO jc = i_startidx, i_endidx

          IF(.NOT. ptr_patch%cells%owner_mask(jc,jb)) CYCLE

          ile = ptr_int%rbf_vec_idx_c(je,jc,jb)
          ibe = ptr_int%rbf_vec_blk_c(je,jc,jb)

          ilc1 = ptr_patch%edges%cell_idx(ile,ibe,1)
          ibc1 = ptr_patch%edges%cell_blk(ile,ibe,1)
          ilc2 = ptr_patch%edges%cell_idx(ile,ibe,2)
          ibc2 = ptr_patch%edges%cell_blk(ile,ibe,2)

          ilcc = ptr_int%rbf_c2grad_idx(jcc,jc,jb)
          ibcc = ptr_int%rbf_c2grad_blk(jcc,jc,jb)

          IF (ilcc == ilc1 .AND. ibcc == ibc1) THEN

            aux_coeff(jc,jcc,1) = aux_coeff(jc,jcc,1) - &
              ptr_int%rbf_vec_coeff_c(je,1,jc,jb) /    &
              ptr_patch%edges%dual_edge_length(ile,ibe)

            aux_coeff(jc,jcc,2) = aux_coeff(jc,jcc,2) - &
              ptr_int%rbf_vec_coeff_c(je,2,jc,jb) /    &
              ptr_patch%edges%dual_edge_length(ile,ibe)

          ELSE IF (ilcc == ilc2 .AND. ibcc == ibc2) THEN

            aux_coeff(jc,jcc,1) = aux_coeff(jc,jcc,1) + &
              ptr_int%rbf_vec_coeff_c(je,1,jc,jb) /    &
              ptr_patch%edges%dual_edge_length(ile,ibe)

            aux_coeff(jc,jcc,2) = aux_coeff(jc,jcc,2) + &
              ptr_int%rbf_vec_coeff_c(je,2,jc,jb) /    &
              ptr_patch%edges%dual_edge_length(ile,ibe)

          ENDIF

        ENDDO !cell loop
      ENDDO
    ENDDO

!CDIR NOLOOPCHG
    DO jcc = 1, rbf_c2grad_dim
      DO jc = i_startidx, i_endidx
        ptr_int%rbf_c2grad_coeff(jcc,1,jc,jb) = aux_coeff(jc,jcc,1)
        ptr_int%rbf_c2grad_coeff(jcc,2,jc,jb) = aux_coeff(jc,jcc,2)
      ENDDO
    ENDDO

  END DO !block loop
!$OMP END DO
!$OMP END PARALLEL

  DO jcc = 1, rbf_c2grad_dim
    call sync_patch_array(SYNC_C, ptr_patch, ptr_int%rbf_c2grad_coeff(jcc,1,:,:))
    call sync_patch_array(SYNC_C, ptr_patch, ptr_int%rbf_c2grad_coeff(jcc,2,:,:))
  ENDDO


END SUBROUTINE rbf_compute_coeff_c2grad

!>
!! This routine computes the coefficients needed for vector RBF interpolation to.
!!
!! This routine computes the coefficients needed for vector RBF interpolation to
!! the vertices, which are then stored in the arrays <i>rbf_vec_coeff_v</i>.
!! This computation involves the inversion of the interpolation matrix, which
!! is performed by a Cholesky decomposition.
!! The Cholesky decomposition is currently implemented by a home made routine
!! which can be substituted by a call to a numerical library, if available.
!!
!! @par Revision History
!! Developed and tested  by Jochen Foerstner (April 2008)
!! @par
!! Modification by Guenther Zaengl, DWD (2009-02-13)
!! - change to direct reconstruction of vector components
!! Modification by Guenther Zaengl, DWD (2009-04-20)
!! - vector optimization
!! Modification by Almut Gassmann, MPI-M (2009-11-06)
!! - renaming subroutine and specifying that the result is for normals
!! Modification by Guenther Zaengl, DWD (2009-11-23)
!! - combine routines for normals and tangentials to avoid code duplication
!! Modification by Almut Gassmann, MPI-M (2010-05-04)
!! - we need only reconstruction from given normal components
!!
SUBROUTINE rbf_vec_compute_coeff_vertex( ptr_patch, ptr_int )
!

!
TYPE(t_patch), TARGET, INTENT(in) :: ptr_patch

TYPE(t_int_state), TARGET, INTENT(inout) :: ptr_int

REAL(wp) :: cc_e1(3), cc_e2(3), cc_v(nproma,3) ! coordinates of edge midpoints

TYPE(t_cartesian_coordinates) :: cc_vertex    ! coordinates of vertex

REAL(wp)           :: z_lon, z_lat          ! longitude and latitude

! 3d  normal velocity vectors at edge midpoints
REAL(wp), DIMENSION (nproma,3) :: z_nx1, z_nx2, z_nx3

REAL(wp)           :: z_norm                ! norm of velocity vectors

REAL(wp)           :: z_dist                ! distance between data points

REAL(wp)           :: z_nxprod              ! scalar product of normal
                                            ! velocity vectors

REAL(wp),ALLOCATABLE :: z_rbfmat(:,:,:)    ! RBF interpolation matrix

REAL(wp),ALLOCATABLE :: z_diag(:,:)        ! diagonal of cholesky
                                           ! decomposition matrix

REAL(wp),ALLOCATABLE :: z_rbfval(:,:)      ! RBF function value

REAL(wp),ALLOCATABLE :: z_rhs1(:,:), &     ! right hand side of linear
                      & z_rhs2(:,:)        ! interpolation system in 2d

INTEGER :: nblks_v
INTEGER :: istencil(nproma)          ! number of edges for the stencil
INTEGER :: jb                        ! index of current block
INTEGER :: jv                        ! index of current vertex
INTEGER :: je1, je2                  ! edge indices (loop)
INTEGER :: ile1, ibe1, ile2, ibe2    ! edge indices
INTEGER :: ist                       ! status variable
INTEGER :: i_startblk                ! start block
INTEGER :: i_startidx                ! start index
INTEGER :: i_endidx                  ! end index
INTEGER :: i_rcstartlev              ! refinement control start level
INTEGER :: i_outunit                 ! base unit for optional debug output
INTEGER :: jg

REAL(wp) ::  checksum_u,checksum_v   ! to check if sum of interpolation coefficients is correct

TYPE(t_cartesian_coordinates), DIMENSION(:,:),   POINTER :: ptr_orient ! pointer to orientation vectors
REAL(wp), DIMENSION(:,:,:,:), POINTER :: ptr_coeff  ! pointer to output coefficients

!--------------------------------------------------------------------

  CALL message('mo_interpolation:rbf_vec_compute_coeff_vertex', '')

  jg = ptr_patch%id

  i_rcstartlev = 2

  ! values for the blocking
  nblks_v  = ptr_patch%nblks_int_v

  ! The start block depends on the width of the stencil
  i_startblk = ptr_patch%verts%start_blk(i_rcstartlev,1)

  ! orientation vectors of input components
  ptr_orient => ptr_patch%edges%primal_cart_normal

  ! RBF coefficients (output)
  ptr_coeff  => ptr_int%rbf_vec_coeff_v

  ! base unit for optional debug output
  i_outunit  = 520


!$OMP PARALLEL PRIVATE (z_rbfmat,z_diag,z_rbfval,z_rhs1,z_rhs2, ist)
  ALLOCATE( z_rbfmat(nproma,rbf_vec_dim_v,rbf_vec_dim_v),  &
            z_diag(nproma,rbf_vec_dim_v),                  &
            z_rbfval(nproma,rbf_vec_dim_v),                &
            z_rhs1(nproma,rbf_vec_dim_v),                  &
            z_rhs2(nproma,rbf_vec_dim_v), STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:rbf_compute_coeff_vertex',  &
      &             'allocation for working arrays failed')
  ENDIF

!$OMP DO PRIVATE (jb,jv,i_startidx,i_endidx,je1,je2,istencil,ist,ile1,ibe1, &
!$OMP             cc_e1,z_lon,z_lat,z_norm,z_nx1,ile2,ibe2,cc_e2,cc_v,      &
!$OMP             z_nx2,z_nxprod,z_dist,cc_vertex,z_nx3,checksum_u,checksum_v)
  DO jb = i_startblk, nblks_v

    CALL get_indices_v(ptr_patch, jb, i_startblk, nblks_v, &
                       i_startidx, i_endidx, i_rcstartlev)

    !
    ! for each vertex, build the vector RBF interpolation matrix
    !
    DO je1 = 1, rbf_vec_dim_v

      DO je2 = 1, je1

        DO jv = i_startidx, i_endidx

          IF(.NOT. ptr_patch%verts%owner_mask(jv,jb)) THEN
            ! Avoid the matrix decomposition for boundary cells since the matrix might get singular
            istencil(jv) = 0
            CYCLE
          ENDIF

          ! Get actual number of stencil points
          istencil(jv) = ptr_int%rbf_vec_stencil_v(jv,jb)
          !
          IF ( (je1 > istencil(jv)) .OR. (je2 > istencil(jv)) ) CYCLE
          !
          ! line and block indices for each edge je1 and je2 of RBF stencil
          !
          ile1 = ptr_int%rbf_vec_idx_v(je1,jv,jb)
          ibe1 = ptr_int%rbf_vec_blk_v(je1,jv,jb)
          ile2 = ptr_int%rbf_vec_idx_v(je2,jv,jb)
          ibe2 = ptr_int%rbf_vec_blk_v(je2,jv,jb)
          !
          ! get Cartesian coordinates and orientation vectors
          !
          cc_e1(:) = ptr_int%cart_edge_coord(ile1,ibe1,:)
          cc_e2(:) = ptr_int%cart_edge_coord(ile2,ibe2,:)
          !
          z_nx1(jv,:) = ptr_orient(ile1,ibe1)%x(:)
          z_nx2(jv,:) = ptr_orient(ile2,ibe2)%x(:)
          !
          ! compute dot product of normal vectors and distance between edge midpoints
          !
          z_nxprod = DOT_PRODUCT(z_nx1(jv,:),z_nx2(jv,:))
          z_dist   = arc_length_v(cc_e1,cc_e2)
          !
          ! set up interpolation matrix
          !
          ! up to now only Gaussian or Inverse multiquadric kernel
          ! (without polynomial correction) are implemented
          IF      (rbf_vec_kern_v == 1) THEN
            z_rbfmat(jv,je1,je2) = z_nxprod * gaussi(z_dist,rbf_vec_scale_v(MAX(jg,1)))
          ELSE IF (rbf_vec_kern_v == 3) THEN
            z_rbfmat(jv,je1,je2) = z_nxprod * inv_multiq(z_dist,rbf_vec_scale_v(MAX(jg,1)))
          ENDIF

          IF (je1 > je2) z_rbfmat(jv,je2,je1) = z_rbfmat(jv,je1,je2)

        END DO

      END DO

    END DO

    ! apply Cholesky decomposition to matrix
    !
!CDIR NOIEXPAND
#ifdef __SX__
    CALL choldec_v(i_startidx,i_endidx,istencil,rbf_vec_dim_v,z_rbfmat,z_diag)
#else
    CALL choldec_v(i_startidx,i_endidx,istencil,              z_rbfmat,z_diag)
#endif

    DO jv = i_startidx, i_endidx

      IF(.NOT. ptr_patch%verts%owner_mask(jv,jb)) CYCLE

      !
      ! Solve immediately for coefficients
      !
      ! convert coordinates of vertex to cartesian vector
      !
      cc_vertex = gc2cc(ptr_patch%verts%vertex(jv,jb))
      cc_v(jv,1:3) = cc_vertex%x(1:3)

      z_lon = ptr_patch%verts%vertex(jv,jb)%lon
      z_lat = ptr_patch%verts%vertex(jv,jb)%lat

      ! Zonal wind component
      CALL gvec2cvec(1._wp,0._wp,z_lon,z_lat,z_nx1(jv,1),z_nx1(jv,2),z_nx1(jv,3))

      z_norm = SQRT( DOT_PRODUCT(z_nx1(jv,:),z_nx1(jv,:)) )
      z_nx1(jv,:)  = 1._wp/z_norm * z_nx1(jv,:)

      ! Meridional wind component
      CALL gvec2cvec(0._wp,1._wp,z_lon,z_lat,z_nx2(jv,1),z_nx2(jv,2),z_nx2(jv,3))

      z_norm = SQRT( DOT_PRODUCT(z_nx2(jv,:),z_nx2(jv,:)) )
      z_nx2(jv,:)  = 1._wp/z_norm * z_nx2(jv,:)

    END DO

    !
    ! set up right hand side for interpolation system
    !
    DO je2 = 1, rbf_vec_dim_v

      DO jv = i_startidx, i_endidx

        IF(.NOT. ptr_patch%verts%owner_mask(jv,jb)) CYCLE

        IF (je2 > istencil(jv)) CYCLE

        !
        ! get indices and coordinates of edge midpoints and compute distance
        ! to vertex
        !
        ile2   = ptr_int%rbf_vec_idx_v(je2,jv,jb)
        ibe2   = ptr_int%rbf_vec_blk_v(je2,jv,jb)
        !
        cc_e2(:)  = ptr_int%cart_edge_coord(ile2,ibe2,:)
        z_dist = arc_length_v(cc_v(jv,:), cc_e2)
        !
        ! get Cartesian orientation vector
        z_nx3(jv,:) = ptr_orient(ile2,ibe2)%x(:)

        ! up to now only Gaussian or Inverse multiquadric kernel
        ! (without polynomial correction) are implemented
        IF      (rbf_vec_kern_v == 1) THEN
          z_rbfval(jv,je2) = gaussi(z_dist,rbf_vec_scale_v(MAX(jg,1)))
        ELSE IF (rbf_vec_kern_v == 3) THEN
          z_rbfval(jv,je2) = inv_multiq(z_dist,rbf_vec_scale_v(MAX(jg,1)))
        ENDIF
        !
        ! compute projection on target vector orientation
        !
        z_rhs1(jv,je2) = z_rbfval(jv,je2) * DOT_PRODUCT(z_nx1(jv,:),z_nx3(jv,:))
        z_rhs2(jv,je2) = z_rbfval(jv,je2) * DOT_PRODUCT(z_nx2(jv,:),z_nx3(jv,:))

      END DO

    END DO

    !
    ! compute vector coefficients
    !
!CDIR NOIEXPAND
#ifdef __SX__
    CALL solve_chol_v(i_startidx, i_endidx, istencil, rbf_vec_dim_v, z_rbfmat,  &
      &               z_diag, z_rhs1, ptr_coeff(:,1,:,jb))
#else
    CALL solve_chol_v(i_startidx, i_endidx, istencil,                z_rbfmat,  &
      &               z_diag, z_rhs1, ptr_coeff(:,1,:,jb))
#endif
!CDIR NOIEXPAND
#ifdef __SX__
    CALL solve_chol_v(i_startidx, i_endidx, istencil, rbf_vec_dim_v, z_rbfmat,  &
      &               z_diag, z_rhs2, ptr_coeff(:,2,:,jb))
#else
    CALL solve_chol_v(i_startidx, i_endidx, istencil,                z_rbfmat,  &
      &               z_diag, z_rhs2, ptr_coeff(:,2,:,jb))
#endif

    DO jv = i_startidx, i_endidx

      IF(.NOT. ptr_patch%verts%owner_mask(jv,jb)) CYCLE

      ! Ensure that sum of interpolation coefficients is correct

      checksum_u = 0._wp
      checksum_v = 0._wp

      DO je1 = 1, istencil(jv)
        ile1   = ptr_int%rbf_vec_idx_v(je1,jv,jb)
        ibe1   = ptr_int%rbf_vec_blk_v(je1,jv,jb)

        z_nx3(jv,:) = ptr_orient(ile1,ibe1)%x(:)

        checksum_u = checksum_u + ptr_coeff(je1,1,jv,jb)* &
          DOT_PRODUCT(z_nx1(jv,:),z_nx3(jv,:))
        checksum_v = checksum_v + ptr_coeff(je1,2,jv,jb)* &
          DOT_PRODUCT(z_nx2(jv,:),z_nx3(jv,:))
      ENDDO

      DO je1 = 1, istencil(jv)
        ptr_coeff(je1,1,jv,jb) = ptr_coeff(je1,1,jv,jb) / checksum_u
        ptr_coeff(je1,2,jv,jb) = ptr_coeff(je1,2,jv,jb) / checksum_v
      ENDDO

    END DO

  END DO
!$OMP END DO

  DEALLOCATE( z_rbfmat, z_diag, z_rbfval, z_rhs1, z_rhs2,  &
    STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:rbf_vec_compute_coeff_vertex',  &
      &             'deallocation for working arrays failed')
  ENDIF
!$OMP END PARALLEL

  DO jb = 1, rbf_vec_dim_v
    call sync_patch_array(SYNC_V, ptr_patch, ptr_int%rbf_vec_coeff_v(jb,1,:,:))
    call sync_patch_array(SYNC_V, ptr_patch, ptr_int%rbf_vec_coeff_v(jb,2,:,:))
  ENDDO


! Optional debug output for RBF coefficients
#ifdef DEBUG_COEFF
  DO jb = i_startblk, nblks_v

    CALL get_indices_v(ptr_patch, jb, i_startblk, nblks_v, &
                       i_startidx, i_endidx, i_rcstartlev)

    DO jv = i_startidx, i_endidx

      istencil = ptr_int%rbf_vec_stencil_v(jv,jb)
      WRITE(i_outunit+jg,'(2i5,25f12.6)') jb,jv,ptr_coeff(1:istencil(jv),1,jv,jb)
      WRITE(i_outunit+jg,'(2i5,25f12.6)') jb,jv,ptr_coeff(1:istencil(jv),2,jv,jb)
    END DO
  END DO
  CLOSE (i_outunit+jg)
#endif

END SUBROUTINE rbf_vec_compute_coeff_vertex

!-------------------------------------------------------------------------
!
!
!>
!! This routine computes the coefficients needed for vector RBF interpolation to.
!!
!! This routine computes the coefficients needed for vector RBF interpolation to
!! the edges, which are then stored in the array <i>rbf_vec_coeff_e</i>.
!! This computation involves the inversion of the interpolation matrix, which
!! is performed by a Cholesky decomposition.
!! The Cholesky decomposition is currently implemented by a home made routine
!! which can be substituted by a call to a numerical library, if available.
!!
!! @par Revision History
!! Developed and tested  by Jochen Foerstner (2008-07-15)
!! @par
!! Modification by Guenther Zaengl, DWD (2009-02-13)
!! - change to direct reconstruction of tangential vector component
!! Modification by Guenther Zaengl, DWD (2009-04-20)
!! - vector optimization
!! Modification by Guenther Zaengl, DWD (2009-11-23)
!! - combine routines for normals and tangentials to avoid code duplication
!! Modification by Almut Gassmann, MPI-M (2010-05-04)
!! - we need only reconstruction from given normal components
!!
SUBROUTINE rbf_vec_compute_coeff_edge( ptr_patch, ptr_int )
!

!
TYPE(t_patch), TARGET, INTENT(in) :: ptr_patch

TYPE(t_int_state), TARGET, INTENT(inout) :: ptr_int

REAL(wp) :: cc_e1(3), cc_e2(3), cc_e(nproma,3) ! coordinates of edge midpoints

TYPE(t_cartesian_coordinates) :: cc_edge      ! coordinates of edge

REAL(wp)           :: z_lon, z_lat          ! longitude and latitude
REAL(wp)           :: z_nu, z_nv            ! zonal and meridional component
                                            ! of normal velocity vectors
                                            ! at edge midpoints

! 3d  normal velocity vectors at edge midpoints
REAL(wp), DIMENSION (nproma,3) :: z_nx1, z_nx2

REAL(wp)           :: z_norm                ! norm of velocity vectors

REAL(wp)           :: z_dist                ! distance between data points

REAL(wp)           :: z_nxprod              ! scalar product of normal
                                            ! velocity vectors

REAL(wp), ALLOCATABLE :: z_rbfmat(:,:,:)    ! RBF interpolation matrix

REAL(wp), ALLOCATABLE :: z_diag(:,:)        ! diagonal of cholesky
                                            ! decomposition matrix

REAL(wp), ALLOCATABLE :: z_rbfval(:,:)      ! RBF function value


INTEGER :: nblks_e
INTEGER :: istencil(nproma)          ! number of edges for the stencil
INTEGER :: jb                        ! index of current block
INTEGER :: je                        ! index of current edge
INTEGER :: je1, je2                  ! edge indices (loop)
INTEGER :: ile1, ibe1, ile2, ibe2    ! edge indices
INTEGER :: ist                       ! status variable
INTEGER :: i_startblk                ! start block
INTEGER :: i_startidx                ! start index
INTEGER :: i_endidx                  ! end index
INTEGER :: i_rcstartlev              ! refinement control start level
INTEGER :: i_outunit                 ! base unit for optional debug output
INTEGER :: jg

REAL(wp) ::  checksum_vt             ! to check if sum of interpolation coefficients is correct

TYPE(t_cartesian_coordinates), DIMENSION(:,:), POINTER :: ptr_orient ! pointer to orientation vectors
REAL(wp), DIMENSION(:,:,:), POINTER :: ptr_coeff  ! pointer to output coefficients

! pointer to orientation of output vectors
TYPE(t_tangent_vectors), DIMENSION(:,:), POINTER :: ptr_orient_out

!--------------------------------------------------------------------

  CALL message('mo_interpolation:rbf_vec_compute_coeff_edge', '')

  jg = ptr_patch%id


  i_rcstartlev = 2

  ! values for the blocking
  nblks_e  = ptr_patch%nblks_int_e

  ! The start block depends on the width of the stencil
  i_startblk = ptr_patch%edges%start_blk(i_rcstartlev,1)

  ! orientation vectors of input components
  ptr_orient => ptr_patch%edges%primal_cart_normal

  ! orientation vectors of output components
  ptr_orient_out => ptr_patch%edges%dual_normal

  ! RBF coefficients (output)
  ptr_coeff  => ptr_int%rbf_vec_coeff_e

  ! base unit for optional debug output
  i_outunit  = 530

!$OMP PARALLEL PRIVATE (z_rbfmat,z_diag,z_rbfval, ist)
  ALLOCATE( z_rbfmat(nproma,rbf_vec_dim_e,rbf_vec_dim_e),  &
            z_diag(nproma,rbf_vec_dim_e),                  &
            z_rbfval(nproma,rbf_vec_dim_e), STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:rbf_vec_compute_coeff_edge_t',  &
      &             'allocation for working arrays failed')
  ENDIF

!$OMP DO PRIVATE (jb,je,i_startidx,i_endidx,je1,je2,istencil,        &
!$OMP    ist,ile1,ibe1,cc_e1,z_nu,z_nv,z_lon,z_lat,z_norm,z_nx1,     &
!$OMP    ile2,ibe2,cc_e2,cc_e,z_nx2,z_nxprod,z_dist,cc_edge,checksum_vt)
  DO jb = i_startblk, nblks_e

    CALL get_indices_e(ptr_patch, jb, i_startblk, nblks_e, &
                       i_startidx, i_endidx, i_rcstartlev)

    !
    ! for each edge, build the vector RBF interpolation matrix
    !
    DO je1 = 1, rbf_vec_dim_e

      DO je2 = 1, je1

        DO je = i_startidx, i_endidx

          IF(.NOT. ptr_patch%edges%owner_mask(je,jb)) THEN
            ! Avoid the matrix decomposition for boundary cells since the matrix might get singular
            istencil(je) = 0
            CYCLE
          ENDIF

          istencil(je) = ptr_int%rbf_vec_stencil_e(je,jb)
          !
          IF ( (je1 > istencil(je)) .OR. (je2 > istencil(je)) ) CYCLE
          !
          ! line and block indices for each edge je1 and je2 of RBF stencil
          !
          ile1 = ptr_int%rbf_vec_idx_e(je1,je,jb)
          ibe1 = ptr_int%rbf_vec_blk_e(je1,je,jb)
          ile2 = ptr_int%rbf_vec_idx_e(je2,je,jb)
          ibe2 = ptr_int%rbf_vec_blk_e(je2,je,jb)
          !
          ! get Cartesian coordinates and orientation vectors
          !
          cc_e1(:) = ptr_int%cart_edge_coord(ile1,ibe1,:)
          cc_e2(:) = ptr_int%cart_edge_coord(ile2,ibe2,:)
          !
          z_nx1(je,:) = ptr_orient(ile1,ibe1)%x(:)
          z_nx2(je,:) = ptr_orient(ile2,ibe2)%x(:)
          !
          ! compute dot product of normal vectors and distance between edge midpoints
          !
          z_nxprod = DOT_PRODUCT(z_nx1(je,:),z_nx2(je,:))
          z_dist   = arc_length_v(cc_e1,cc_e2)

          ! set up interpolation matrix
          !
          ! up to now only Gaussian or Inverse multiquadric kernel
          ! (without polynomial correction) are implemented
          IF      (rbf_vec_kern_e == 1) THEN
            z_rbfmat(je,je1,je2) = z_nxprod * gaussi(z_dist,rbf_vec_scale_e(MAX(jg,1)))
          ELSE IF (rbf_vec_kern_e == 3) THEN
            z_rbfmat(je,je1,je2) = z_nxprod * inv_multiq(z_dist,rbf_vec_scale_e(MAX(jg,1)))
          ENDIF

          IF (je1 > je2) z_rbfmat(je,je2,je1) = z_rbfmat(je,je1,je2)

        END DO

      END DO

    END DO


    ! apply Cholesky decomposition to matrix
    !
!CDIR NOIEXPAND
#ifdef __SX__
    CALL choldec_v(i_startidx,i_endidx,istencil,rbf_vec_dim_e,z_rbfmat,z_diag)
#else
    CALL choldec_v(i_startidx,i_endidx,istencil,              z_rbfmat,z_diag)
#endif

    DO je = i_startidx, i_endidx

      IF(.NOT. ptr_patch%edges%owner_mask(je,jb)) CYCLE

      !
      ! Solve immediately for coefficients
      !
      ! convert coordinates of edge midpoint to cartesian vector
      !
      cc_edge = gc2cc(ptr_patch%edges%center(je,jb))
      cc_e(je,1:3) = cc_edge%x(1:3)

      z_nu  = ptr_orient_out(je,jb)%v1
      z_nv  = ptr_orient_out(je,jb)%v2
      z_lon = ptr_patch%edges%center(je,jb)%lon
      z_lat = ptr_patch%edges%center(je,jb)%lat

      CALL gvec2cvec(z_nu,z_nv,z_lon,z_lat,z_nx1(je,1),z_nx1(je,2),z_nx1(je,3))

      z_norm = SQRT( DOT_PRODUCT(z_nx1(je,:),z_nx1(je,:)) )
      z_nx1(je,:)  = 1._wp/z_norm * z_nx1(je,:)

    END DO

    !
    ! set up right hand side for interpolation system
    !
    DO je2 = 1, rbf_vec_dim_e

      DO je = i_startidx, i_endidx

        IF(.NOT. ptr_patch%edges%owner_mask(je,jb)) CYCLE

        IF (je2 > istencil(je)) CYCLE
        !
        ! get indices and coordinates of edge midpoints and compute distance
        ! to the edge where the vector is reconstructed
        !
        ile2   = ptr_int%rbf_vec_idx_e(je2,je,jb)
        ibe2   = ptr_int%rbf_vec_blk_e(je2,je,jb)

        cc_e2(:)  = ptr_int%cart_edge_coord(ile2,ibe2,:)
        z_dist = arc_length_v(cc_e(je,:), cc_e2)
        !
        ! get Cartesian orientation vector
        z_nx2(je,:) = ptr_orient(ile2,ibe2)%x(:)

        z_nxprod = DOT_PRODUCT(z_nx1(je,:),z_nx2(je,:))

        ! up to now only Gaussian or Inverse multiquadric kernel
        ! (without polynomial correction) are implemented
        IF      (rbf_vec_kern_e == 1) THEN
          z_rbfval(je,je2) = gaussi(z_dist,rbf_vec_scale_e(MAX(jg,1))) * z_nxprod
        ELSE IF (rbf_vec_kern_e == 3) THEN
          z_rbfval(je,je2) = inv_multiq(z_dist,rbf_vec_scale_e(MAX(jg,1))) * z_nxprod
        ENDIF

      END DO

    END DO

    !
    ! compute vector coefficients
    !
!CDIR NOIEXPAND
#ifdef __SX__
    CALL solve_chol_v(i_startidx, i_endidx, istencil, rbf_vec_dim_e, z_rbfmat, &
                    z_diag, z_rbfval, ptr_coeff(:,:,jb))
#else
    CALL solve_chol_v(i_startidx, i_endidx, istencil,                z_rbfmat, &
                    z_diag, z_rbfval, ptr_coeff(:,:,jb))
#endif

    DO je = i_startidx, i_endidx

      IF(.NOT. ptr_patch%edges%owner_mask(je,jb)) CYCLE

      ! Ensure that sum of interpolation coefficients is correct

      checksum_vt = 0._wp

      DO je1 = 1, istencil(je)
        ile1   = ptr_int%rbf_vec_idx_e(je1,je,jb)
        ibe1   = ptr_int%rbf_vec_blk_e(je1,je,jb)

        z_nx2(je,:) = ptr_orient(ile1,ibe1)%x(:)

        checksum_vt = checksum_vt+ptr_coeff(je1,je,jb)*&
          DOT_PRODUCT(z_nx1(je,:),z_nx2(je,:))
      ENDDO

      DO je1 = 1, istencil(je)
        ptr_coeff(je1,je,jb) = ptr_coeff(je1,je,jb) / checksum_vt
      ENDDO

    END DO

  END DO
!$OMP END DO

  DEALLOCATE( z_rbfmat, z_diag, z_rbfval, STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:rbf_vec_compute_coeff_edge',  &
      &             'deallocation for working arrays failed')
  ENDIF
!$OMP END PARALLEL

  DO jb = 1, rbf_vec_dim_e
    call sync_patch_array(SYNC_E, ptr_patch, ptr_int%rbf_vec_coeff_e(jb,:,:))
  ENDDO


! Optional debug output for RBF coefficients
#ifdef DEBUG_COEFF
  DO jb = i_startblk, nblks_e

    CALL get_indices_e(ptr_patch, jb, i_startblk, nblks_e, &
                       i_startidx, i_endidx, i_rcstartlev)

    DO je = i_startidx, i_endidx

      istencil(je) = ptr_int%rbf_vec_stencil_e(je,jb)
      WRITE(i_outunit+jg,'(2i5,25f12.6)') jb,je,ptr_coeff(1:istencil(je),je,jb)
    END DO
  END DO
  CLOSE (i_outunit+jg)
#endif

END SUBROUTINE rbf_vec_compute_coeff_edge


!-------------------------------------------------------------------------
!> Computes the RBF coefficients for lon-lat grid points.
!-------------------------------------------------------------------------
!
! This routine is based on mo_intp_rbf_coeffs::rbf_vec_compute_coeff_cell()
! 
! @par Revision History
!      Initial implementation  by  F.Prill, DWD (2011-08)
!
SUBROUTINE rbf_vec_compute_coeff_lonlat( ptr_patch, ptr_int, ptr_int_lonlat, lon_lat_points, &
  &                                      nblks_lonlat, npromz_lonlat )

  ! Input parameters
  TYPE(t_patch),                 INTENT(IN)    :: ptr_patch
  TYPE (t_int_state),    TARGET, INTENT(INOUT) :: ptr_int
  TYPE (t_lon_lat_intp),         INTENT(INOUT) :: ptr_int_lonlat
  REAL(gk),                      INTENT(IN)    :: lon_lat_points(:,:,:)
  INTEGER,                       INTENT(IN)    :: nblks_lonlat, npromz_lonlat ! blocking info
  ! Local parameters
  CHARACTER(*), PARAMETER :: routine = TRIM("mo_interpolation:rbf_vec_compute_coeff_lonlat")
  REAL(wp)                         :: cc_e1(3), cc_e2(3), cc_c(nproma,3)  ! coordinates of edge midpoints
  TYPE(t_cartesian_coordinates)    :: cc_center                  ! coordinates of cell centers
  REAL(wp), DIMENSION (nproma,3)   :: z_nx1, z_nx2, z_nx3        ! 3d  normal velocity vectors at edge midpoints
  REAL(wp)                         :: z_lon, z_lat,            & ! longitude and latitude
    &                                 z_norm,                  & ! norm of velocity vectors
    &                                 z_dist,                  & ! distance between data points
    &                                 z_nxprod                   ! scalar prod of normal veloc.
  REAL(wp),ALLOCATABLE             :: z_rbfmat(:,:,:),         & ! RBF interpolation matrix
    &                                 z_diag(:,:),             & ! diagonal of Cholesky decomp.
    &                                 z_rbfval(:,:),           & ! RBF function value
    &                                 z_rhs1(:,:),             & ! right hand side of linear
    &                                 z_rhs2(:,:)                ! interpolation system in 2d
  INTEGER                          :: &
    &                                 jc, jb,                  & ! integer over lon-lat grid points
    &                                 je1, je2,                & ! integer over edges
    &                                 ile1, ibe1, ile2, ibe2,  & ! edge indices
    &                                 ist,                     & ! return value of array allocation
    &                                 i_startidx, i_endidx,    & ! start/end index
    &                                 i_rcstartlev,            & ! refinement control start level
    &                                 istencil(nproma),        & ! actual number of stencil points
    &                                 jg
  REAL(wp)                         :: checksum_u,checksum_v      ! to check if sum of interpolation coefficients is correct
  TYPE(t_geographical_coordinates) :: grid_point

  !--------------------------------------------------------------------

  CALL message(routine, '')
  IF (ptr_patch%n_patch_cells == 0) RETURN;

  jg = ptr_patch%id

  i_rcstartlev = 2

!$OMP PARALLEL PRIVATE (z_rbfmat,z_diag,z_rbfval,z_rhs1,z_rhs2, ist,   &
!$OMP                   grid_point),                                   &
!$OMP          SHARED  (nproma, rbf_vec_dim_c, nblks_lonlat,           &
!$OMP                   npromz_lonlat, ptr_int_lonlat, ptr_patch,      &
!$OMP                   rbf_vec_kern_c, jg, rbf_vec_scale_c )

  ALLOCATE( z_rbfmat(nproma,rbf_vec_dim_c,rbf_vec_dim_c),  &
            z_diag(nproma,rbf_vec_dim_c),                  &
            z_rbfval(nproma,rbf_vec_dim_c),                &
            z_rhs1(nproma,rbf_vec_dim_c),                  &
            z_rhs2(nproma,rbf_vec_dim_c), STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish (routine, 'allocation for working arrays failed')
  ENDIF

!$OMP DO PRIVATE (jb,jc,i_startidx,i_endidx,je1,je2,istencil,       &
!$OMP             ist,ile1,ibe1,cc_e1,z_lon,z_lat,z_norm,           &
!$OMP             z_nx1,ile2,ibe2,cc_e2,cc_c,z_nx2,z_nxprod,z_dist, &
!$OMP             cc_center,z_nx3,checksum_u,checksum_v)
  BLOCKS: DO jb = 1, nblks_lonlat

    i_startidx = 1
    i_endidx   = nproma
    if (jb == nblks_lonlat) i_endidx = npromz_lonlat

    ! for each cell, build the vector RBF interpolation matrix
    DO je1 = 1, rbf_vec_dim_c
      DO je2 = 1, je1
        DO jc = i_startidx, i_endidx

          ! Get actual number of stencil points
          istencil(jc) = ptr_int_lonlat%rbf_vec_stencil(jc,jb)

          ! paranoia:
          IF ( (je1 > istencil(jc)) .OR. (je2 > istencil(jc)) ) CYCLE

          ! line and block indices for each edge je1 and je2 of RBF stencil
          ile1 = ptr_int_lonlat%rbf_vec_idx(je1,jc,jb)
          ibe1 = ptr_int_lonlat%rbf_vec_blk(je1,jc,jb)
          ile2 = ptr_int_lonlat%rbf_vec_idx(je2,jc,jb)
          ibe2 = ptr_int_lonlat%rbf_vec_blk(je2,jc,jb)

          ! get Cartesian coordinates and orientation vectors
          cc_e1(:) = ptr_int%cart_edge_coord(ile1,ibe1,:)
          cc_e2(:) = ptr_int%cart_edge_coord(ile2,ibe2,:)
          !
          z_nx1(jc,:) = ptr_patch%edges%primal_cart_normal(ile1,ibe1)%x(:)
          z_nx2(jc,:) = ptr_patch%edges%primal_cart_normal(ile2,ibe2)%x(:)

          ! compute dot product of normal vectors and distance between edge midpoints
          z_nxprod = DOT_PRODUCT(z_nx1(jc,:),z_nx2(jc,:))
          z_dist   = arc_length_v(cc_e1,cc_e2)

          ! set up interpolation matrix
          IF      (rbf_vec_kern_c == 1) THEN
            z_rbfmat(jc,je1,je2) = z_nxprod * gaussi(z_dist,rbf_vec_scale_c(MAX(jg,1)))
          ELSE IF (rbf_vec_kern_c == 3) THEN
            z_rbfmat(jc,je1,je2) = z_nxprod * inv_multiq(z_dist,rbf_vec_scale_c(MAX(jg,1)))
          ENDIF

          IF (je1 > je2) z_rbfmat(jc,je2,je1) = z_rbfmat(jc,je1,je2)

        END DO
      END DO
    END DO

    ! apply Cholesky decomposition to matrix
    !
!CDIR NOIEXPAND
#ifdef __SX__
    CALL choldec_v(i_startidx,i_endidx,istencil,rbf_vec_dim_c,z_rbfmat,z_diag)
#else
    CALL choldec_v(i_startidx,i_endidx,istencil,              z_rbfmat,z_diag)
#endif

    DO jc = i_startidx, i_endidx

      !
      ! Solve immediately for coefficients
      !
      ! convert coordinates of cell center to cartesian vector
      !
      grid_point%lon = REAL( lon_lat_points(jc, jb,1), wp)
      grid_point%lat = REAL( lon_lat_points(jc, jb,2), wp)
      cc_center = gc2cc(grid_point)
      cc_c(jc,1:3) = cc_center%x(1:3)

      z_lon = grid_point%lon
      z_lat = grid_point%lat

      ! Zonal wind component
      CALL gvec2cvec(1._wp,0._wp, z_lon,z_lat, &
        &            z_nx1(jc,1),z_nx1(jc,2),z_nx1(jc,3))

      z_norm = SQRT( DOT_PRODUCT(z_nx1(jc,:),z_nx1(jc,:)) )
      z_nx1(jc,:)  = 1._wp/z_norm * z_nx1(jc,:)

      ! Meridional wind component
      CALL gvec2cvec(0._wp,1._wp, z_lon,z_lat, &
        &            z_nx2(jc,1),z_nx2(jc,2),z_nx2(jc,3))

      z_norm = SQRT( DOT_PRODUCT(z_nx2(jc,:),z_nx2(jc,:)) )
      z_nx2(jc,:)  = 1._wp/z_norm * z_nx2(jc,:)

    END DO

    !
    ! set up right hand side for interpolation system
    !
    DO je2 = 1, rbf_vec_dim_c
      DO jc = i_startidx, i_endidx

        IF (je2 > istencil(jc)) CYCLE
      
        ! get indices and coordinates of edge midpoints and compute distance
        ! to cell center
        ile2   = ptr_int_lonlat%rbf_vec_idx(je2,jc,jb)
        ibe2   = ptr_int_lonlat%rbf_vec_blk(je2,jc,jb)

        cc_e2(:)  = ptr_int%cart_edge_coord(ile2,ibe2,:)

        z_dist = arc_length_v(cc_c(jc,:), cc_e2)

        ! get Cartesian orientation vector
        z_nx3(jc,:) = ptr_patch%edges%primal_cart_normal(ile2,ibe2)%x(:)

        IF (rbf_vec_kern_c == 1) THEN
          z_rbfval(jc,je2) = gaussi(z_dist,rbf_vec_scale_c(MAX(jg,1)))
        ELSE IF (rbf_vec_kern_c == 3) THEN
          z_rbfval(jc,je2) = inv_multiq(z_dist,rbf_vec_scale_c(MAX(jg,1)))
        ENDIF

        ! compute projection on target vector orientation
        z_rhs1(jc,je2) = z_rbfval(jc,je2) * &
          &                     DOT_PRODUCT(z_nx1(jc,:),z_nx3(jc,:))
        z_rhs2(jc,je2) = z_rbfval(jc,je2) * &
          &                     DOT_PRODUCT(z_nx2(jc,:),z_nx3(jc,:))

      END DO
    END DO

    ! compute vector coefficients
!CDIR NOIEXPAND
#ifdef __SX__
    CALL solve_chol_v(i_startidx, i_endidx, istencil, rbf_vec_dim_c, z_rbfmat,  &
      &               z_diag, z_rhs1, ptr_int_lonlat%rbf_vec_coeff(:,1,:,jb))
#else
    CALL solve_chol_v(i_startidx, i_endidx, istencil,                z_rbfmat,  &
      &               z_diag, z_rhs1, ptr_int_lonlat%rbf_vec_coeff(:,1,:,jb))
#endif
!CDIR NOIEXPAND
#ifdef __SX__
    CALL solve_chol_v(i_startidx, i_endidx, istencil, rbf_vec_dim_c, z_rbfmat,  &
      &               z_diag, z_rhs2, ptr_int_lonlat%rbf_vec_coeff(:,2,:,jb))
#else
    CALL solve_chol_v(i_startidx, i_endidx, istencil,                z_rbfmat,  &
      &               z_diag, z_rhs2, ptr_int_lonlat%rbf_vec_coeff(:,2,:,jb))
#endif

    DO jc = i_startidx, i_endidx

      ! Ensure that sum of interpolation coefficients is correct
      checksum_u = 0._wp
      checksum_v = 0._wp

      DO je1 = 1, istencil(jc)
        ile1   = ptr_int_lonlat%rbf_vec_idx(je1,jc,jb)
        ibe1   = ptr_int_lonlat%rbf_vec_blk(je1,jc,jb)

        ! get Cartesian orientation vector
        z_nx3(jc,:) = ptr_patch%edges%primal_cart_normal(ile1,ibe1)%x(:)

        checksum_u = checksum_u + ptr_int_lonlat%rbf_vec_coeff(je1,1,jc,jb)* &
          DOT_PRODUCT(z_nx1(jc,:),z_nx3(jc,:))
        checksum_v = checksum_v + ptr_int_lonlat%rbf_vec_coeff(je1,2,jc,jb)* &
          DOT_PRODUCT(z_nx2(jc,:),z_nx3(jc,:))
      ENDDO

      DO je1 = 1, istencil(jc)
        ptr_int_lonlat%rbf_vec_coeff(je1,1,jc,jb) = &
          ptr_int_lonlat%rbf_vec_coeff(je1,1,jc,jb) / checksum_u
      END DO
      DO je1 = 1, istencil(jc)
        ptr_int_lonlat%rbf_vec_coeff(je1,2,jc,jb) = &
          ptr_int_lonlat%rbf_vec_coeff(je1,2,jc,jb) / checksum_v
      END DO

    END DO ! jc

  END DO BLOCKS
!$OMP END DO

  DEALLOCATE( z_rbfmat, z_diag, z_rbfval, z_rhs1, z_rhs2,  STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish (routine, 'deallocation for working arrays failed')
  ENDIF
!$OMP END PARALLEL

END SUBROUTINE rbf_vec_compute_coeff_lonlat


!> Computed combined coefficient/finite difference matrix for lon-lat interpolation.
!! This routine computes the coefficients needed for reconstructing the gradient
!! of a cell-based variable at the cell center. The operations performed here
!! combine taking the centered-difference gradient (like in grad_fd_norm) at
!! the 9 edges entering into the usual 9-point RBF vector reconstruction at
!! cell centers, followed by applying the RBF reconstruction.
!!
!! @par Revision History
!! based on
!!   mo_intp_rbf_coeffs::rbf_compute_coeff_c2grad
!! developed and tested by Guenther Zaengl, DWD (2009-12-15)
!! Restructuring F. Prill, DWD (2011-08-18)
!!
SUBROUTINE rbf_compute_coeff_c2grad_lonlat (ptr_patch,                   &
  &                                         nblks_lonlat, npromz_lonlat, &
  &                                         ptr_int_lonlat)

  TYPE(t_patch),     INTENT(in)    :: ptr_patch                   ! patch on which computation is performed
  INTEGER,           INTENT(IN)    :: nblks_lonlat, npromz_lonlat ! lon-lat grid blocking info
  TYPE (t_lon_lat_intp), TARGET, INTENT(inout) :: ptr_int_lonlat  ! interpolation state

  ! local variables
  INTEGER  :: je, jcc,                                   &
    &         jc, jb,                                    &  ! integer over lon-lat points
    &         i_startidx, i_endidx,                      &
    &         ile, ibe, ilc, ibc, ilcc, ibcc
  REAL(wp) :: inv_length, coeff(2),                      &
    &         aux_coeff(nproma,rbf_c2grad_dim,2)
  REAL(wp), DIMENSION(:,:,:,:), POINTER :: ptr_coeff

  !--------------------------------------------------------------------

  ptr_coeff => ptr_int_lonlat%rbf_vec_coeff

  ! loop through all patch cells (and blocks)

!$OMP PARALLEL 
!$OMP DO PRIVATE(jb,je,jcc,jc,i_startidx,i_endidx,ile,ibe,ilc,ibc,&
!$OMP            aux_coeff, inv_length, coeff, ilcc, ibcc)
  DO jb = 1,nblks_lonlat

    i_startidx = 1
    i_endidx   = nproma
    IF (jb == nblks_lonlat) i_endidx = npromz_lonlat

    aux_coeff(:,:,:) = 0._wp

    DO je = 1, rbf_vec_dim_c
      DO jcc = 1, rbf_c2grad_dim
!CDIR NODEP
        DO jc = i_startidx, i_endidx

          ile = ptr_int_lonlat%rbf_vec_idx(je,jc,jb)
          ibe = ptr_int_lonlat%rbf_vec_blk(je,jc,jb)

          ilcc = ptr_int_lonlat%rbf_c2grad_idx(jcc,jc,jb)
          ibcc = ptr_int_lonlat%rbf_c2grad_blk(jcc,jc,jb)

          inv_length = ptr_patch%edges%inv_dual_edge_length(ile,ibe)
          coeff(1:2) = ptr_coeff(je,1:2,jc,jb) * inv_length

          ! every edge (ile,ibe) of the stencil for cell (jc,jb)
          ! contributes to two entries of the combined coeff/finite
          ! difference matrix block "aux_coeff" (with the sign "-/+1")

          ilc = ptr_patch%edges%cell_idx(ile,ibe,1)
          ibc = ptr_patch%edges%cell_blk(ile,ibe,1)

          IF (ilcc == ilc .AND. ibcc == ibc) THEN
            aux_coeff(jc,jcc,1:2) = aux_coeff(jc,jcc,1:2) - coeff(1:2)
          ELSE
            ilc = ptr_patch%edges%cell_idx(ile,ibe,2)
            ibc = ptr_patch%edges%cell_blk(ile,ibe,2)
            IF (ilcc == ilc .AND. ibcc == ibc) THEN
              aux_coeff(jc,jcc,1:2) = aux_coeff(jc,jcc,1:2) + coeff(1:2)
            END IF
          END IF
        END DO ! jc
      END DO ! jcc
    ENDDO !je

    FORALL (jcc = 1:rbf_c2grad_dim)
      ptr_int_lonlat%rbf_c2grad_coeff(jcc,1,i_startidx:i_endidx,jb) = &
        & aux_coeff(i_startidx:i_endidx,jcc,1)
      ptr_int_lonlat%rbf_c2grad_coeff(jcc,2,i_startidx:i_endidx,jb) = &
        & aux_coeff(i_startidx:i_endidx,jcc,2)
    END FORALL

  END DO !block loop
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE rbf_compute_coeff_c2grad_lonlat


!-------------------------------------------------------------------------
!> Setup routine for RBF reconstruction at lon-lat grid points.
! 
! @par Revision History
!      Initial implementation  by  F.Prill, DWD (2011-08)
!      Changed to a wrapper by Rainer Johanni   (2011-11)
!
!
SUBROUTINE rbf_setup_interpol_lonlat(k_jg, ptr_patch, ptr_int_lonlat, ptr_int)
  
  INTEGER,               INTENT(IN)            :: k_jg           ! patch index
  ! data structure containing grid info:
  TYPE(t_patch), TARGET, INTENT(IN)            :: ptr_patch
  ! Indices of source points and interpolation coefficients
  TYPE (t_lon_lat_intp), TARGET, INTENT(INOUT) :: ptr_int_lonlat
  TYPE (t_int_state),    TARGET, INTENT(INOUT) :: ptr_int

  CALL rbf_setup_interpol_lonlat_grid(lonlat_intp_config(k_jg)%lonlat_grid, &
                                      ptr_patch, ptr_int_lonlat, ptr_int)

END SUBROUTINE rbf_setup_interpol_lonlat

!-------------------------------------------------------------------------
!> Setup routine for RBF reconstruction at lon-lat grid points for an arbitrary grid.
! 
! @par Revision History
!      Initial implementation  by  F.Prill, DWD (2011-08)
!      Changed for abritrary grids by Rainer Johanni (2011-11)
!
! Note: So far, the RBF coefficients and indices are by all compute PEs
!       though they are only required by the writing PEs.
! Note also that asynchronous IO is not yet supported.
!
SUBROUTINE rbf_setup_interpol_lonlat_grid(grid, ptr_patch, ptr_int_lonlat, ptr_int)
  
  TYPE (t_lon_lat_grid), INTENT(INOUT)         :: grid
  ! data structure containing grid info:
  TYPE(t_patch), TARGET, INTENT(IN)            :: ptr_patch
  ! Indices of source points and interpolation coefficients
  TYPE (t_lon_lat_intp), TARGET, INTENT(INOUT) :: ptr_int_lonlat
  TYPE (t_int_state),    TARGET, INTENT(INOUT) :: ptr_int
  ! Local Parameters:
  CHARACTER(*), PARAMETER :: routine = TRIM("mo_intp_rbf_coeffs:rbf_setup_interpol_lonlat")
  REAL(wp), ALLOCATABLE            :: rotated_pts(:,:,:)
  REAL(gk), ALLOCATABLE            :: in_points(:,:,:)
  REAL(gk), ALLOCATABLE            :: min_dist(:,:)       ! minimal distance
  INTEGER                          :: jb, jc, i_startidx, i_endidx,         &
    &                                 nblks_lonlat, npromz_lonlat, errstat, &
    &                                 jb_lonlat, jc_lonlat,                 &
    &                                 block_size, glb_index, i,             &
    &                                 rl_start, rl_end, i_nchdom, i_startblk
  LOGICAL                          :: l_distributed, l_cutoff_local_domains
  LOGICAL, ALLOCATABLE             :: l_cutoff(:,:)
  INTEGER                          :: array_shape_2d(2), array_shape_3d(3), &
    &                                 array_shape_4d(4)
  TYPE(t_geographical_coordinates) :: cell_center, lonlat_pt
  TYPE(t_cartesian_coordinates)    :: p1, p2
  REAL(wp)                         :: point(2), z_norm, z_nx1(3), z_nx2(3), &
    &                                 max_dist
  INTEGER                          :: ithis_local_pts

  !-----------------------------------------------------------------------

  ! TODO[FP] : It shouldn't be harmful to enable distributed
  ! coefficient computation even in sequential mode.
  l_distributed = .TRUE.
  ! Flag: .TRUE., if we want to erase values outside local domains
  l_cutoff_local_domains = .TRUE.

  IF (dbg_level > 1) THEN
    WRITE(message_text,*) "SETUP : rbf_interpol_lonlat"
    CALL message(routine, message_text)
  END IF
  IF (ptr_patch%cell_type == 6) THEN
    CALL finish(routine, "Lon-lat interpolation not yet implemented for cell_type == 6!")
  END IF

  nblks_lonlat  = grid%nblks
  npromz_lonlat = grid%npromz

  !-- compute grid points of rotated lon/lat grid

  ALLOCATE(in_points(nproma, nblks_lonlat, 2),            &
    &      rotated_pts(grid%dimen(1), grid%dimen(2), 2),  &
    &      min_dist(nproma, nblks_lonlat),                &
    &      stat=errstat)
  IF (errstat /= SUCCESS) THEN
    CALL finish (routine, 'allocation for working arrays failed')
  ENDIF

  IF (dbg_level > 1) &
    CALL message(routine, "rotate lon-lat grid points")
  CALL rotate_latlon_grid( grid, rotated_pts )

  in_points(:,:,1) = RESHAPE(REAL(rotated_pts(:,:,1), gk), &
    &                        (/ nproma, nblks_lonlat /), (/ 0._gk /) )
  in_points(:,:,2) = RESHAPE(REAL(rotated_pts(:,:,2), gk), &
    &                        (/ nproma, nblks_lonlat /), (/ 0._gk /) )

  !-- proximity query, build a list of cell indices that contain the
  !   lon/lat grid points

  ! build GNAT data structure
  IF (dbg_level > 1) &
    CALL message(routine, "build GNAT data structure")
  CALL gnat_init_grid(ptr_patch)

  ! Perform query. Note that for distributed patches we receive a
  ! local list of "in_points" actually located on this portion of the
  ! domain.
  IF (dbg_level > 1) &
    CALL message(routine, "proximity query")
  CALL gnat_query_containing_triangles(ptr_patch, gnat_tree, in_points(:,:,:), &
    &                                  nproma, nblks_lonlat, npromz_lonlat,    &
    &                                  ptr_int_lonlat%tri_idx(:,:,:), min_dist(:,:))
  IF (l_distributed) THEN
    CALL gnat_merge_distributed_queries(ptr_patch, grid%total_dim, nproma, grid%nblks, min_dist, &
      &                                 ptr_int_lonlat%tri_idx(:,:,:), in_points(:,:,:),         &
      &                                 ptr_int_lonlat%nlocal_pts(:), ptr_int_lonlat%owner(:),   &
      &                                 ithis_local_pts)
    ! set local values for "nblks" and "npromz"
    nblks_lonlat  = (ithis_local_pts-1)/nproma + 1
    npromz_lonlat = ithis_local_pts - (nblks_lonlat-1)*nproma
  END IF
  ! clean up
  CALL gnat_destroy()

  !-- edge indices of stencil
  ! are available through previous calls to SR rbf_vec_index_cell()
  ! and rbf_c2grad_index()
  
  !-- copy neighbor indices from standard RBF interpolation
  IF (dbg_level > 1) &
    CALL message(routine, "copy neighbor indices from standard RBF interpolation")

!$OMP PARALLEL SHARED(nblks_lonlat, nproma, npromz_lonlat, &
!$OMP                 ptr_int_lonlat, ptr_int),            &
!$OMP          DEFAULT (NONE)
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc), SCHEDULE(runtime)
  DO jb_lonlat = 1,nblks_lonlat

    i_startidx = 1
    i_endidx   = nproma
    IF (jb_lonlat == nblks_lonlat) i_endidx = npromz_lonlat

    DO jc_lonlat = i_startidx, i_endidx

      jc = ptr_int_lonlat%tri_idx(1,jc_lonlat, jb_lonlat)
      jb = ptr_int_lonlat%tri_idx(2,jc_lonlat, jb_lonlat)

      ptr_int_lonlat%rbf_vec_stencil(jc_lonlat,jb_lonlat) = ptr_int%rbf_vec_stencil_c(jc,jb)      
      ptr_int_lonlat%rbf_vec_idx(:,jc_lonlat,jb_lonlat)   = ptr_int%rbf_vec_idx_c(:,jc,jb)
      ptr_int_lonlat%rbf_vec_blk(:,jc_lonlat,jb_lonlat)   = ptr_int%rbf_vec_blk_c(:,jc,jb)
      
      ptr_int_lonlat%rbf_c2grad_idx(:,jc_lonlat,jb_lonlat) = ptr_int%rbf_c2grad_idx(:,jc,jb)
      ptr_int_lonlat%rbf_c2grad_blk(:,jc_lonlat,jb_lonlat) = ptr_int%rbf_c2grad_blk(:,jc,jb)

    ENDDO
  ENDDO
!$OMP END DO
!$OMP END PARALLEL

  !-- compute interpolation coefficients
  IF (dbg_level > 1) &
    CALL message(routine, "compute lon-lat interpolation coefficients")
  CALL rbf_vec_compute_coeff_lonlat( ptr_patch, ptr_int, ptr_int_lonlat,  &
    &                                in_points, nblks_lonlat, npromz_lonlat )

  ! compute distances (x0i - xc)
  DO jb=1,nblks_lonlat
    i_startidx = 1
    i_endidx   = nproma
    IF (jb == nblks_lonlat) i_endidx = npromz_lonlat

    DO jc=i_startidx,i_endidx
      point(:)    = REAL(in_points(jc,jb,:), wp)
      lonlat_pt%lon = point(1)
      lonlat_pt%lat = point(2)
      cell_center = ptr_patch%cells%center(ptr_int_lonlat%tri_idx(1,jc,jb),      &
        &                                  ptr_int_lonlat%tri_idx(2,jc,jb))

      ! convert to Cartesian coordinate system:
      p1 = gc2cc(lonlat_pt)
      p2 = gc2cc(cell_center)
      p1%x(:) = p1%x(:) - p2%x(:)

      ! Zonal component
      CALL gvec2cvec(1._wp,0._wp, point(1), point(2), &
        &            z_nx1(1),z_nx1(2),z_nx1(3))
      z_norm = SQRT( DOT_PRODUCT(z_nx1(:),z_nx1(:)) )
      z_nx1(:)  = 1._wp/z_norm * z_nx1(:)
      ! Meridional component
      CALL gvec2cvec(0._wp, 1._wp, point(1), point(2), &
        &            z_nx2(1),z_nx2(2),z_nx2(3))
      z_norm = SQRT( DOT_PRODUCT(z_nx2(:),z_nx2(:)) )
      z_nx2(:)  = 1._wp/z_norm * z_nx2(:)
      ! projection, scale with earth radius
      ptr_int_lonlat%rdist(1, jc, jb) = re*DOT_PRODUCT(p1%x(:), z_nx1)
      ptr_int_lonlat%rdist(2, jc, jb) = re*DOT_PRODUCT(p1%x(:), z_nx2)
    END DO
  END DO

  CALL rbf_compute_coeff_c2grad_lonlat (ptr_patch, nblks_lonlat, &
    &                                   npromz_lonlat, ptr_int_lonlat)

  IF (l_distributed) THEN

    ! translate local indices to global indices:
    DO jb=1,nblks_lonlat
      i_startidx = 1
      i_endidx   = nproma
      IF (jb == nblks_lonlat) i_endidx = npromz_lonlat

      DO jc=i_startidx,i_endidx
        glb_index = ptr_patch%cells%glb_index(idx_1d(ptr_int_lonlat%tri_idx(1,jc,jb), &
          &                                          ptr_int_lonlat%tri_idx(2,jc,jb)))
        ptr_int_lonlat%tri_idx(:, jc,jb) = (/ idx_no(glb_index), blk_no(glb_index) /)

        DO i=1,rbf_vec_dim_c
          glb_index = ptr_patch%edges%glb_index(idx_1d(ptr_int_lonlat%rbf_vec_idx(i,jc,jb), &
            &                                          ptr_int_lonlat%rbf_vec_blk(i,jc,jb)))
          ptr_int_lonlat%rbf_vec_idx(i,jc,jb) = idx_no(glb_index)
          ptr_int_lonlat%rbf_vec_blk(i,jc,jb) = blk_no(glb_index)
        END DO

        DO i=1,rbf_c2grad_dim
          glb_index = ptr_patch%cells%glb_index(idx_1d(ptr_int_lonlat%rbf_c2grad_idx(i,jc,jb), &
            &                                          ptr_int_lonlat%rbf_c2grad_blk(i,jc,jb)))
          ptr_int_lonlat%rbf_c2grad_idx(i,jc,jb) = idx_no(glb_index)
          ptr_int_lonlat%rbf_c2grad_blk(i,jc,jb) = blk_no(glb_index)
        END DO
      END DO
    END DO

    ! Gather local components - as a result, ALL PEs will own the complete set of coefficients

    CALL p_gather_field(grid%total_dim, ptr_int_lonlat%nlocal_pts, ptr_int_lonlat%owner, &
      &                 ptr_int_lonlat%rbf_vec_coeff)

    CALL p_gather_field(grid%total_dim, ptr_int_lonlat%nlocal_pts, ptr_int_lonlat%owner, &
      &                 ptr_int_lonlat%rbf_c2grad_coeff)

    CALL p_gather_field(grid%total_dim, ptr_int_lonlat%nlocal_pts, ptr_int_lonlat%owner, &
      &                 ptr_int_lonlat%rbf_vec_idx)
    CALL p_gather_field(grid%total_dim, ptr_int_lonlat%nlocal_pts, ptr_int_lonlat%owner, &
      &                 ptr_int_lonlat%rbf_vec_blk)

    CALL p_gather_field(grid%total_dim, ptr_int_lonlat%nlocal_pts, ptr_int_lonlat%owner, &
      &                 ptr_int_lonlat%rbf_vec_stencil)

    CALL p_gather_field(grid%total_dim, ptr_int_lonlat%nlocal_pts, ptr_int_lonlat%owner, &
      &                 ptr_int_lonlat%rbf_c2grad_idx)
    CALL p_gather_field(grid%total_dim, ptr_int_lonlat%nlocal_pts, ptr_int_lonlat%owner, &
      &                 ptr_int_lonlat%rbf_c2grad_blk)

    CALL p_gather_field(grid%total_dim, ptr_int_lonlat%nlocal_pts, ptr_int_lonlat%owner, &
      &                 ptr_int_lonlat%rdist)

    CALL p_gather_field(grid%total_dim, ptr_int_lonlat%nlocal_pts, ptr_int_lonlat%owner, &
      &                 ptr_int_lonlat%tri_idx)

    ! In some cases, some of the lon-lat grid points are associated to none of
    ! the PEs. We get index entries "0" then, which must be corrected by 
    ! setting "1" as dummy index:
    ptr_int_lonlat%rbf_vec_idx    = MAX(1,   ptr_int_lonlat%rbf_vec_idx)
    ptr_int_lonlat%rbf_vec_blk    = MAX(1,   ptr_int_lonlat%rbf_vec_blk)
    ptr_int_lonlat%rbf_c2grad_idx = MAX(1,ptr_int_lonlat%rbf_c2grad_idx)
    ptr_int_lonlat%rbf_c2grad_blk = MAX(1,ptr_int_lonlat%rbf_c2grad_blk)
    ptr_int_lonlat%tri_idx        = MAX(1,       ptr_int_lonlat%tri_idx)
  END IF

  IF ( l_cutoff_local_domains ) THEN
    ! make a sensible guess for maximum distance (just taking a
    ! "randomly chosen" edge length)
    rl_start = 2
    rl_end   = min_rlcell_int
    i_nchdom   = MAX(1,ptr_patch%n_childdom)
    i_startblk = ptr_patch%cells%start_blk(rl_start,1)
    CALL get_indices_e(ptr_patch, i_startblk,  &
      &                i_startblk, i_startblk, &
      &                i_startidx, i_endidx, &
      &                rl_start, rl_end)
    max_dist = 3._gk * REAL(ptr_patch%edges%primal_edge_length(i_startidx,i_startblk)/re, gk)

    ! build a logical array, .true. where lon-lat point is outside of
    ! the current patch:
    ALLOCATE(l_cutoff(nproma, grid%nblks), stat=errstat)
    IF (errstat /= SUCCESS) &
      CALL finish (routine, 'allocation for working arrays failed')
    l_cutoff(:,:) = (MAX(ABS(ptr_int_lonlat%rdist(1,:,:)), &
      &                  ABS(ptr_int_lonlat%rdist(2,:,:))) >= max_dist)

    WHERE (l_cutoff(:,:))
      ptr_int_lonlat%rdist(1,:,:) = 0.
      ptr_int_lonlat%rdist(2,:,:) = 0.
    END WHERE

    ! clean up
    DEALLOCATE(l_cutoff, stat=errstat)
    IF (errstat /= SUCCESS) &
      CALL finish (routine, 'DEALLOCATE of working arrays failed')
  END IF

  DEALLOCATE(in_points, rotated_pts, min_dist, stat=errstat)
  IF (errstat /= SUCCESS) THEN
    CALL finish (routine, 'DEALLOCATE of working arrays failed')
  ENDIF

END SUBROUTINE rbf_setup_interpol_lonlat_grid


END MODULE mo_intp_rbf_coeffs
