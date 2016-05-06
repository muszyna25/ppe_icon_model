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
!!  Modification by Constantin Junk, MPI-M (2011-05-05)
!!  - moved interpol_ctl namelist variables to namelists/mo_interpol_ctl
!!
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
MODULE mo_intp_data_strc
!-------------------------------------------------------------------------
!
!    ProTeX FORTRAN source: Style 2
!    modified for ICON project, DWD/MPI-M 2006
!
!-------------------------------------------------------------------------

USE mo_kind,                ONLY: wp
USE mo_exception,           ONLY: finish
USE mo_impl_constants,      ONLY: min_rlcell, min_rledge, min_rlvert, max_dom
USE mo_math_utilities,      ONLY: t_cartesian_coordinates,   &
  &                               t_geographical_coordinates
USE mo_lonlat_grid,         ONLY: t_lon_lat_grid
USE mo_communication,       ONLY: t_comm_gather_pattern

IMPLICIT NONE


! NOTE: The variables will be use along the mo_interpolation sub-modules
!       They are declared to be public
PUBLIC

CHARACTER(LEN=*), PARAMETER :: modname = 'mo_intp_data_strc'

REAL(wp) :: sick_a, sick_o      ! if i_cori_method >= 2: To avoid the SICK instability
                                ! (Symmetric Instability of Computational Kind or
                                !  Hollingsworth instability), an average of the kinetic
                                ! energy and thus of the mass flux must be defined.
                                ! sick_a is to be given in the namelist as the weight of
                                ! the fully averaged kinetic energy, sick_o=1-sick_a.

TYPE t_lsq_set
  LOGICAL :: l_consv            ! flag to determine whether the least squares
                                ! reconstruction should be conservative

  INTEGER :: dim_c,        &    ! parameter determining the size of the lsq stencil
    &        dim_unk            ! parameter determining the dimension of the solution
                                ! vector (== number of unknowns) of the lsq system

  INTEGER :: wgt_exp            ! least squares weighting exponent
END TYPE t_lsq_set


TYPE t_lsq
  ! fields related to weighted least squares polynomial reconstruction
  !---------------------------------------------------------------------
  INTEGER, ALLOCATABLE  :: lsq_dim_stencil(:,:)  ! stencil size as a function of jc and jb
                                                 ! necessary in order to account for pentagon
                                                 ! points.
  INTEGER, ALLOCATABLE  :: lsq_idx_c(:,:,:)      ! index array defining the stencil for
                                                 ! lsq reconstruction (nproma,nblks_c,lsq_dim_c)
  INTEGER, ALLOCATABLE  :: lsq_blk_c(:,:,:)      ! dito for the blocks
  REAL(wp), ALLOCATABLE :: lsq_weights_c(:,:,:)  ! weights for lsq reconstruction
                                                 ! (nproma,lsq_dim_c,nblks_c)
  REAL(wp), ALLOCATABLE :: lsq_qtmat_c(:,:,:,:)  ! transposed Q of QR-factorization for
                                                 ! lsq reconstruction
                                                 ! (nproma,lsq_dim_unk,lsq_dim_c,nblks_c)
  REAL(wp), ALLOCATABLE :: lsq_rmat_rdiag_c(:,:,:)! reciprocal diagonal elements of R-matrix
                                                  ! resulting from QR-decomposition
                                                  ! (nproma,lsq_dim_unk,nblks_c)
  REAL(wp), ALLOCATABLE :: lsq_rmat_utri_c(:,:,:)! upper triangular elements without diagonal
                                                 ! elements of R-matrix (starting from the bottom
                                                 ! right)
                                                 ! (nproma,(lsq_dim_unk^2-lsq_dim_unk)/2,nblks_c)
  REAL(wp), ALLOCATABLE :: lsq_pseudoinv(:,:,:,:)! pseudo (or Moore-Penrose) inverse of lsq
                                                 ! design matrix A
                                                 ! (nproma,lsq_dim_unk,lsq_dim_c,nblks_c)
  REAL(wp), ALLOCATABLE :: lsq_moments(:,:,:)      ! Moments (x^ny^m)_{i} for control volume
                                                   ! (nproma,nblks_c,lsq_dim_unk)
  REAL(wp), ALLOCATABLE :: lsq_moments_hat(:,:,:,:)! Moments (\hat{x^ny^m}_{ij}) for control volume
                                                   ! (nproma,nblks_c,lsq_dim_c,lsq_dim_unk)
END TYPE t_lsq


TYPE t_gauss_quad
  !
  ! quadrature points for intergration over triangular element
  !
  TYPE(t_geographical_coordinates), ALLOCATABLE :: & !< linear (nproma, nblks_c)
    &  qpts_tri_l(:,:)
  TYPE(t_geographical_coordinates), ALLOCATABLE :: & !< quadratic (nproma, nblks_c,3)
    &  qpts_tri_q(:,:,:)
  TYPE(t_geographical_coordinates), ALLOCATABLE :: & !< cubic (nproma, nblks_c,4)
    &  qpts_tri_c(:,:,:)
  REAL(wp), ALLOCATABLE :: weights_tri_q(:)      !< quadratic weight (3)
  REAL(wp), ALLOCATABLE :: weights_tri_c(:)      !< cubic weights (4)
END TYPE t_gauss_quad


TYPE t_int_state

  ! a) weights which are inconsistent with the Hamiltonian viewpoint
  !-----------------------------------------------------------------

  REAL(wp), ALLOCATABLE :: c_lin_e(:,:,:)   ! coefficient for interpolation
                                            ! from adjacent cells onto edge
                                            ! (nproma,2,nblks_e)

  REAL(wp), ALLOCATABLE :: e_bln_c_s(:,:,:) ! coefficient for bilinear
                                            ! interpolation from edges to cells
                                            ! for scalar quantities

  REAL(wp), ALLOCATABLE :: e_bln_c_u(:,:,:) ! coefficient for bilinear interpolation
                                            ! from edges to cells for vector components
                                            ! (input: v_t, v_n, output: u)

  REAL(wp), ALLOCATABLE :: e_bln_c_v(:,:,:) ! coefficient for bilinear interpolation
                                            ! from edges to cells for vector components
                                            ! (input: v_t, v_n, output: v)

  REAL(wp), ALLOCATABLE :: c_bln_avg(:,:,:) ! coefficients for bilinear divergence
                                            ! averaging (nproma,4,nblks_c)

  REAL(wp), ALLOCATABLE :: e_flx_avg(:,:,:) ! coefficients for related velocity or mass flux
                                            ! averaging (nproma,5,nblks_e)

  REAL(wp), ALLOCATABLE :: v_1o2_e(:,:,:)   ! coefficient for interpolation
                                            ! from vertices onto edges by 1/2
                                            ! weighting (nproma,2,nblks_e),

  REAL(wp), ALLOCATABLE :: cells_plwa_verts(:,:,:)! coefficient for interpolation
                                            ! from cells to verts by pseudo-Laplacian
                                            ! weighted averaging (PLWA)
  REAL(wp), ALLOCATABLE :: gradc_bmat(:,:,:,:) ! Bmatrix for cell centered shape function based
                                            ! gradient (nproma,2,3,nblks_c)


  ! b) weights which are consistent with the Hamiltonian viewpoint
  !---------------------------------------------------------------
  ! The following weights are needed for the mass and theta brackets

  REAL(wp), ALLOCATABLE :: e_inn_c(:,:,:)   ! coefficient for inner product
                                            ! of 2 vector components
                                            ! from edges to cells

  REAL(wp), ALLOCATABLE :: e_inn_v(:,:,:)   ! coefficient for inner product
                                            ! of 2 vector components
                                            ! from edges to verts

  REAL(wp), ALLOCATABLE :: e_aw_c(:,:,:)    ! coefficient for scalar interp
                                            ! from edges to cells

  REAL(wp), ALLOCATABLE :: r_aw_c(:,:,:)    ! coefficient for scalar interp
                                            ! from rhombi to cells

  REAL(wp), ALLOCATABLE :: e_aw_v(:,:,:)    ! coefficient for scalar interp
                                            ! from edges to vertices

  REAL(wp), ALLOCATABLE :: e_1o3_v(:,:,:)   ! coefficient for hexagonal grid
                                            ! averaging from edges to vertices

  REAL(wp), ALLOCATABLE :: tria_aw_rhom(:,:,:)! coefficient for interpolation
                                            ! from triangles to rhombi by area
                                            ! weighting

  REAL(wp), ALLOCATABLE :: verts_aw_cells(:,:,:)! coefficient for interpolation
                                            ! from vertices to cells by
                                            ! area weighting

  REAL(wp), ALLOCATABLE :: cells_aw_verts(:,:,:)! coefficient for interpolation
                                            ! from cells to verts by
                                            ! area weighting

  ! c) RBF related fields
  !----------------------
  INTEGER, ALLOCATABLE  :: rbf_vec_idx_c(:,:,:)  ! index array defining the
                                            ! stencil of surrounding edges for
                                            ! vector rbf interpolation at each
                                            ! cell center
                                            ! (rbf_vec_dim_c,nproma,nblks_c)
  INTEGER, ALLOCATABLE  :: rbf_vec_blk_c(:,:,:)  ! ... dito for the blocks

  INTEGER, ALLOCATABLE  :: rbf_vec_stencil_c(:,:) ! array defining number of
                                            ! surrounding edges in the stencil
                                            ! for vector rbf interpolation at
                                            ! each cell center
                                            ! (nproma,nblks_c)
  REAL(wp), ALLOCATABLE :: rbf_vec_coeff_c(:,:,:,:) ! array containing the
                                            ! coefficients used for
                                            ! vector rbf interpolation
                                            ! at each cell center
                                            ! (rbf_vec_dim_c,2,nproma,nblks_c)

  INTEGER, ALLOCATABLE  :: rbf_c2grad_idx(:,:,:)  ! index array defining the
                                            ! stencil of surrounding cells for
                                            ! 2D gradient reconstruction at each
                                            ! cell center
                                            ! (rbf_c2grad_dim,nproma,nblks_c)
  INTEGER, ALLOCATABLE  :: rbf_c2grad_blk(:,:,:)  ! ... dito for the blocks

  REAL(wp), ALLOCATABLE :: rbf_c2grad_coeff(:,:,:,:) ! array containing the
                                            ! coefficients used for
                                            ! 2D gradient reconstruction
                                            ! at each cell center
                                            ! (rbf_c2grad_dim,2,nproma,nblks_c)

  INTEGER, ALLOCATABLE  :: rbf_vec_idx_v(:,:,:) ! index array defining the
                                            ! stencil of surrounding edges for
                                            ! vector rbf interpolation at each
                                            ! triangle vertex
                                            ! (rbf_vec_dim_v,nproma,nblks_v)
  INTEGER, ALLOCATABLE  :: rbf_vec_blk_v(:,:,:) ! ... dito for the blocks

  INTEGER, ALLOCATABLE  :: rbf_vec_stencil_v(:,:) ! array defining number of
                                            ! surrounding edges in the stencil
                                            ! for vector rbf interpolation at
                                            ! each triangle vertex
                                            ! (nproma,nblks_v)

  REAL(wp), ALLOCATABLE :: rbf_vec_coeff_v(:,:,:,:) ! array containing the
                                            ! coefficients used for vector rbf
                                            ! interpolation at each tringle
                                            ! vertex (input is normal component)
                                            ! (rbf_vec_dim_v,2,nproma,nblks_v)

  INTEGER, ALLOCATABLE  :: rbf_vec_idx_e(:,:,:) ! index array defining the
                                            ! stencil of surrounding edges for
                                            ! vector rbf interpolation at each
                                            ! triangle edge
                                            ! (rbf_vec_dim_e,nproma,nblks_e)
  INTEGER, ALLOCATABLE  :: rbf_vec_blk_e(:,:,:) ! ... dito for the blocks

  INTEGER, ALLOCATABLE  :: rbf_vec_stencil_e(:,:) ! array defining number of
                                            ! surrounding edges in the stencil
                                            ! for vector rbf interpolation at
                                            ! each triangle edge
                                            ! (nproma,nblks_e)

  REAL(wp), ALLOCATABLE :: rbf_vec_coeff_e(:,:,:) ! array containing the
                                            ! coefficients used for rbf inter-
                                            ! polation of the tangential velo-
                                            ! city component (from the
                                            ! surrounding normals) at each
                                            ! triangle edge
                                            ! (rbf_vec_dim_e,nproma,nblks_e)

  ! d) fields needed for reconstructions in hexagonal model
  !--------------------------------------------------------
  REAL(wp), ALLOCATABLE :: heli_coeff(:,:,:)  ! coefficients for lamb_rot computation
  INTEGER, ALLOCATABLE  :: heli_vn_idx(:,:,:) ! len indices for vn in lamb_rot
  INTEGER, ALLOCATABLE  :: heli_vn_blk(:,:,:) ! blk indices for vn in lamb_rot
  REAL(wp), ALLOCATABLE :: hex_north(:,:,:)   ! coeffs for north vector in hexagon center
  REAL(wp), ALLOCATABLE :: hex_east(:,:,:)    ! coeffs for east vector in hexagon center
  REAL(wp), ALLOCATABLE :: tria_north(:,:,:)  ! coeffs for north vector in triangle center
  REAL(wp), ALLOCATABLE :: tria_east(:,:,:)   ! coeffs for east vector in triangle center
  REAL(wp), ALLOCATABLE :: quad_north(:,:,:)  ! coeffs for north vector in rhombus center
  REAL(wp), ALLOCATABLE :: quad_east(:,:,:)   ! coeffs for east vector in rhomus center
  REAL(wp), ALLOCATABLE :: cno_en(:,:,:)      ! coeffs for north vector in edge center
  REAL(wp), ALLOCATABLE :: cea_en(:,:,:)      ! coeffs for east vector in edge center

  ! e) precomputed geometrical factors for mathematical operators (for efficiency)
  !------------------------------------------------------------------------------
  REAL(wp), ALLOCATABLE :: geofac_div(:,:,:)    ! factor for divergence (nproma,cell_type,nblks_c)
  REAL(wp), ALLOCATABLE :: geofac_qdiv(:,:,:)   ! factor for quad-cell divergence (nproma,4,nblks_e)
  REAL(wp), ALLOCATABLE :: geofac_grdiv(:,:,:)  ! factor for gradient of divergence (nproma,5,nblks_e)
  REAL(wp), ALLOCATABLE :: geofac_rot(:,:,:)    ! factor for divergence (nproma,9-cell_type,nblks_v)
  REAL(wp), ALLOCATABLE :: geofac_n2s(:,:,:)    ! factor for nabla2-scalar (nproma,cell_type+1,nblks_c)
  REAL(wp), ALLOCATABLE :: geofac_grg(:,:,:,:)  ! factor for Green-Gauss gradient (nproma,4,nblks_c,2)

  ! f) precomputed Cartesian orientation and location vectors of edge midpoints
  !    and location of cell centers(for efficiency) : it is now computed in grid genrator stored
  !    in p_patch
  !------------------------------------------------------------------------------

  ! g) patch elements restored from edges to cells to reduce frequency of indirect addressing
  !------------------------------------------------------------------------------
  REAL(wp), ALLOCATABLE :: primal_normal_ec(:,:,:,:) ! p_patch%edges%primal_normal_cell stored on
                                                     ! the cell data type (nproma,nblks_c,3,2)
  REAL(wp), ALLOCATABLE :: edge_cell_length(:,:,:)   ! p_patch%edges%edge_cell_length stored on
                                                     ! the cell data type (nproma,nblks_c,3)

  ! h) distance from cells to vertices on local cartesian grid with origin at the cell center
  !    (used for gradient limiter)
  !------------------------------------------------------------------------------
  REAL(wp), ALLOCATABLE :: cell_vert_dist(:,:,:,:)   ! (nproma,3,2,nblks_c)

  ! i) fields related to calculation of backward trajectories on local plane
  !    tangential to the edge midpoint
  !------------------------------------------------------------------------------
  REAL(wp), ALLOCATABLE :: pos_on_tplane_e(:,:,:,:)  ! positions of various points on local plane
                                                     ! tangential to the edge midpoint.
                                                     ! currently: (nproma,nblks_e,8,2)
  REAL(wp), ALLOCATABLE :: tplane_e_dotprod(:,:,:,:) ! Dot product between unit vectors at inner
                                                     ! edge midpoint and edge midpoints of
                                                     ! corresponding quadrilateral cell.
                                                     ! (nproma,nblks_e,4,4)
  TYPE(t_geographical_coordinates), ALLOCATABLE ::  &! positions of vertices and butterfly
    &  pos_on_tplane_c_edge(:,:,:,:)                 ! neighbors on local plane tangential to the
                                                     ! cell circumcenter.
                                                     ! stored in an edge-based data structure
                                                     ! (nproma,nblks_e,2,5)


  ! i) fields related to weighted least squares polynomial reconstruction
  !------------------------------------------------------------------------------
  TYPE(t_lsq) :: lsq_lin,  &  ! coefficients for linear lsq-reconstruction
    &            lsq_high     ! coefficients for higher order lsq-reconstruction


  ! j) fields related to third order advection and Smagorinski diffusion (on hexagons):
  ! directional laplacian ( as du/dx > ux, dv/dy > vy, du/dx=dv/dy > xy )
  !----------------------------------------------------------------------
  INTEGER , ALLOCATABLE :: &
    & dir_gradh_i1(  :,:,:) ! index array for edges of neighbor hexagon cell 1 of considered
                            ! edge (6  , nproma, nblks_e)
  INTEGER , ALLOCATABLE :: &
    & dir_gradh_b1(  :,:,:) ! block array for edges of neighbor hexagon cell 1 of considered
                            ! edge (6  , nproma, nblks_e)
  REAL(wp), ALLOCATABLE :: &
    & dir_gradhux_c1(:,:,:)

  REAL(wp), ALLOCATABLE :: &
    & strain_def_c1( :,:,:) ! coeff array for edges of neighbor hexagon cell 1 of considered
                            ! edge (6  , nproma, nblks_e)
  INTEGER , ALLOCATABLE :: &
    & dir_gradh_i2(  :,:,:) ! index array for edges of neighbor hexagon cell 2 of considered
                            ! edge (6  , nproma, nblks_e)
  INTEGER , ALLOCATABLE :: &
    & dir_gradh_b2(  :,:,:) ! block array for edges of neighbor hexagon cell 2 of considered
                            ! edge (6  , nproma, nblks_e)
  REAL(wp), ALLOCATABLE :: &
    & dir_gradhux_c2(:,:,:)
  REAL(wp), ALLOCATABLE :: &
    & strain_def_c2( :,:,:) ! coeff array for edges of neighbor hexagon cell 2 of considered&
                            ! edge (6  , nproma, nblks_e)

  INTEGER , ALLOCATABLE :: &
    & dir_gradt_i1(  :,:,:) ! index array for edges of neighbor triangle cell 1 of considered
                            ! edge (12, nproma, nblks_e)
  INTEGER , ALLOCATABLE :: &
    & dir_gradt_b1(  :,:,:) ! block array for edges of neighbor triangle cell 1 of considered
                            ! edge (12, nproma, nblks_e)
  REAL(wp), ALLOCATABLE :: &
    & dir_gradtxy_v1(:,:,:), dir_gradtyx_v1(:,:,:)
  REAL(wp), ALLOCATABLE :: &
    & shear_def_v1(  :,:,:) ! coeff array for edges of neighbor triangle cell 1 of considered
                            ! edge (12, nproma, nblks_e)
  INTEGER , ALLOCATABLE :: &
    & dir_gradt_i2(  :,:,:) ! index array for edges of neighbor triangle cell 2 of considered
                            ! edge (12, nproma, nblks_e)
  INTEGER , ALLOCATABLE :: &
    & dir_gradt_b2(  :,:,:) ! block array for edges of neighbor triangle cell 2 of considered
                            ! edge (12, nproma, nblks_e)
  REAL(wp), ALLOCATABLE :: dir_gradtxy_v2(:,:,:), dir_gradtyx_v2(:,:,:)
  REAL(wp), ALLOCATABLE :: &
    & shear_def_v2(  :,:,:) ! coeff array for edges of neighbor triangle cell 2 of considered
                            ! edge (12, nproma, nblks_e)

  ! k) Nudging coefficients used for 1-way nesting and limited-area mode (defined here
  !    rather than in grf_state because the limited-area mode may be used without nesting)
  !------------------------------------------------------------------------------
  REAL(wp), ALLOCATABLE :: nudgecoeff_c(:,:)  !< Nudging coefficient for cells
  REAL(wp), ALLOCATABLE :: nudgecoeff_e(:,:)  !< Nudging coefficient for cells


  ! l) Quadrature points and weights for integration over triangular element
  !--------------------------------------------------------------------------
  TYPE(t_gauss_quad) ::gquad

  ! m) Ocean variables
  ! This will be transfered into an separate variable AFTER we got it running in parallel
  ! The following two arrays are required for the reconstruction process that
  ! is used within the ocean model. Once the new version is implemented this
  ! could eventually be shifted to the edge/vertex datatypes. It is currently
  ! placed here to reduce interference with the atmospheric code (P.K.).
  !
  ! Vector pointing from cell circumcenter to edge midpoint. In the associated
  ! cell2edge_weight-array the cell2edge_vec is multiplied by some other geometric
  ! quantities (edge-length, cell area). The weight is used in the reconstruction
  ! the vector is used in the transposed reconstruction.
  ! index=1,nproma, index2=1,nblks_c, index3=1,3
  ! other choice would be index2=1,nblks_e, index3=1,2
  ! Eventually switch to other second indexing if this is more appropriate

  ! Vector pointing from vertex (dual center) to midpoint of dual edge
  ! (/= midpoint of primal edge).
  ! In the associated vertex2dualedge_mid_weight-array the vertex2dualedge_mid_vec
  ! is multiplied by some other geometric quantities (dual edge-length, dual cell
  ! area). The weight is used in the reconstruction the vector is used in the
  ! transposed reconstruction.
  ! index=1,nproma, index2=1,nblks_v, index3=1,6
  ! other choice index2=1,nblks_e, index3=1,2
  ! Eventually switch to other second indexing if this is more appropriate
  ! new constructs for mimetic core:
  TYPE(t_cartesian_coordinates), ALLOCATABLE :: edge2cell_coeff_cc(:,:,:)
  TYPE(t_cartesian_coordinates), ALLOCATABLE :: edge2cell_coeff_cc_t(:,:,:)


  !REAL(wp),                      ALLOCATABLE :: edge2vert_coeff(:,:,:,:)
  !TYPE(t_cartesian_coordinates), ALLOCATABLE :: edge2vert_coeff_t(:,:,:)

  TYPE(t_cartesian_coordinates), ALLOCATABLE :: edge2vert_coeff_cc(:,:,:)
  TYPE(t_cartesian_coordinates), ALLOCATABLE :: edge2vert_coeff_cc_t(:,:,:)
  TYPE(t_cartesian_coordinates), ALLOCATABLE :: edge2vert_vector_cc(:,:,:)

  REAL(wp), ALLOCATABLE :: fixed_vol_norm(:,:)
  REAL(wp), ALLOCATABLE :: variable_vol_norm(:,:,:)
  REAL(wp), ALLOCATABLE :: variable_dual_vol_norm(:,:,:)


  ! Location of midpoint of dual edge
  !!$    TYPE(t_geographical_coordinates), ALLOCATABLE :: mid_dual_edge(:,:)
  ! Cartesian distance from vertex1 to vertex2 via dual edge midpoint
  REAL(wp), ALLOCATABLE :: dist_cell2edge(:,:,:)
END TYPE t_int_state


!> data structure containing coefficients for (optional) interpolation
!> onto lon-lat grid.
TYPE t_lon_lat_intp

  ! --- Radial Basis Function (RBF) interpolation

  REAL(wp)              :: rbf_scale                  ! RBF shape parameter

  REAL(wp), ALLOCATABLE :: rbf_vec_coeff(:,:,:,:)     ! array containing the
                                                      ! coefficients used for
                                                      ! vector rbf interpolation
                                                      ! at lon-lat grid points
                                                      ! (rbf_vec_dim_c,2,nproma,nblks_lonlat)

  REAL(wp), ALLOCATABLE :: rbf_c2l_coeff(:,:,:)       ! array containing the coefficients used
                                                      ! for direct interpolation from cell centers
                                                      ! to lon-lat points. Uses the same stencil as
                                                      ! the rbf_c2grad algorithm

  INTEGER, ALLOCATABLE  :: rbf_vec_idx(:,:,:)         ! index array defining the
                                                      ! stencil of surrounding edges for
                                                      ! vector rbf interpolation at each
                                                      ! cell center
                                                      ! (rbf_vec_dim_c,nproma,nblks_lonlat)

  INTEGER, ALLOCATABLE  :: rbf_vec_blk(:,:,:)         ! ... dito for the blocks

  INTEGER, ALLOCATABLE  :: rbf_vec_stencil(:,:)       ! array defining number of
                                                      ! surrounding edges in the stencil
                                                      ! for vector rbf interpolation at
                                                      ! each cell center
                                                      ! (nproma,nblks_lonlat)

  ! --- nearest-neighbor interpolation
  !
  !     (we need interpolation weights also here to blank out points)

  REAL(wp), ALLOCATABLE  :: nnb_coeff(:,:,:)          ! array containing nearest-neighbor interpolation 
                                                      ! weights (1,nproma,nblks_lonlat)


  ! --- barycentric interpolation

  INTEGER,  ALLOCATABLE  :: baryctr_idx(:,:,:)        ! index array defining the
                                                      ! stencil for barycentric interpolation
                                                      ! (3,nproma,nblks_lonlat)

  INTEGER,  ALLOCATABLE  :: baryctr_blk(:,:,:)        ! ... dito for the blocks

  REAL(wp), ALLOCATABLE  :: baryctr_coeff(:,:,:)      ! array containing barycentric interpolation 
                                                      ! weights
                                                      ! (3,nproma,nblks_lonlat)


  ! --- direct RBF interpolation from cell centers to lon-lat points:

  INTEGER, ALLOCATABLE  :: rbf_c2l_idx(:,:,:)         ! (rbf_dim_c2l, nproma, nblks_c)
  INTEGER, ALLOCATABLE  :: rbf_c2l_blk(:,:,:)         ! (rbf_dim_c2l, nproma, nblks_c)
  INTEGER, ALLOCATABLE  :: rbf_c2lr_idx(:,:,:)        ! (rbf_dim_c2l, nproma, nblks_lonlat)
  INTEGER, ALLOCATABLE  :: rbf_c2lr_blk(:,:,:)        ! (rbf_dim_c2l, nproma, nblks_lonlat)
  INTEGER, ALLOCATABLE  :: rbf_c2l_stencil(:,:)       ! (nproma, nblks_c)
  INTEGER, ALLOCATABLE  :: rbf_c2lr_stencil(:,:)      ! (nproma, nblks_lonlat)

  ! --- other data fields

  ! list of triangles containing lon-lat grid points (first dim: index and block)
  INTEGER, ALLOCATABLE  :: tri_idx(:,:,:) ! 2, nproma, nblks_lonlat
  ! coordinates of the lon-lat points (nproma,nb nblks_lonlat)
  TYPE(t_geographical_coordinates), ALLOCATABLE :: ll_coord(:,:)

  ! data fields for distributed computations (available on all PEs)
  INTEGER               :: nthis_local_pts  ! number of points local to this PE
  INTEGER, ALLOCATABLE  :: global_idx(:)    ! for each lon-lat point on this PE: global idx

END TYPE t_lon_lat_intp


!> Collects all interpolation coefficients together with the
!> corresponding lon-lat grid and a communication pattern for
!> gathering.
!
TYPE t_lon_lat_data
  TYPE(t_lon_lat_grid)        :: grid
  LOGICAL                     :: l_dom        (max_dom), &
    &                            l_initialized(max_dom)
  TYPE(t_lon_lat_intp)        :: intp         (max_dom)
  TYPE(t_comm_gather_pattern) :: p_pat(max_dom)
END TYPE t_lon_lat_data


TYPE(t_int_state),TARGET,ALLOCATABLE :: p_int_state(:), p_int_state_local_parent(:)

!> maximum no. of lon-lat grids used in this model
!
!  @note This value must be at least as large as the number of used
!        output_nml name lists!
INTEGER, PARAMETER    :: MAX_LONLAT_GRIDS = 30

!> Global list of lon-lat grids and interpolation coefficients.
!
!  All lon-lat grids needed for output are stored in this
!  registry. Several subroutines collaborate by this data structure:
!  coefficient computation, interpolation algorithm, I/O, dump/restore
!
! @note This list has a static size but not all data members are
!       necessarily allocated.
!
! @note On the compute PEs information is *local*, i.e. only for each
!       lon-lat grid only those parts "owned" by the process are
!       stored.
TYPE (t_lon_lat_data), TARGET, SAVE :: lonlat_grid_list(MAX_LONLAT_GRIDS)

!> Actual no. of lon-lat grids currently used in this model
INTEGER               :: n_lonlat_grids

INTERFACE OPERATOR (==)
  MODULE PROCEDURE lonlat_grid_compare
END INTERFACE OPERATOR (==)

CONTAINS

  !---------------------------------------------------------------
  !> @return index with "free" lon-lat grid data
  !
  FUNCTION get_free_lonlat_grid()
    INTEGER :: get_free_lonlat_grid
    CHARACTER(*), PARAMETER :: routine = modname//"::get_free_lonlat_grid"

    IF (n_lonlat_grids == MAX_LONLAT_GRIDS) THEN
      CALL finish(routine, "Maximum number of lon-lat grids exceeded!")
    END IF
    n_lonlat_grids = n_lonlat_grids + 1
    get_free_lonlat_grid = n_lonlat_grids
  END FUNCTION get_free_lonlat_grid


  !---------------------------------------------------------------
  ! @return index in lon-lat grid list for a given grid, -1 if
  !         not in list.
  !
  FUNCTION get_lonlat_grid_ID(grid)
    INTEGER :: get_lonlat_grid_ID
    TYPE(t_lon_lat_grid), INTENT(IN) :: grid
    ! local variables
    CHARACTER(*), PARAMETER :: routine = modname//"::get_lonlat_grid_ID"
    INTEGER :: i

    ! default value: "not found"
    get_lonlat_grid_ID = -1

    DO i=1,n_lonlat_grids
      IF (lonlat_grid_list(i)%grid == grid) THEN
        get_lonlat_grid_ID = i
        EXIT
      END IF
    END DO
  END FUNCTION get_lonlat_grid_ID


  !---------------------------------------------------------------
  ! Test two floating point arrays for equality.
  !
  FUNCTION float_cmp(arr1, arr2, zero_tol)
    LOGICAL :: float_cmp
    REAL(wp), INTENT(IN) :: arr1(:), arr2(:)
    REAL(wp), INTENT(IN) :: zero_tol
    ! local variables
    INTEGER :: i

    float_cmp = .TRUE.
    IF (SIZE(arr1) /= SIZE(arr2)) THEN
      float_cmp = .FALSE.
      RETURN
    END IF
    DO i=1,SIZE(arr1)
      IF (ABS(arr1(i) - arr2(i)) > zero_tol) THEN
        float_cmp = .FALSE.
        RETURN
      END IF
    END DO
  END FUNCTION float_cmp


  !---------------------------------------------------------------
  ! Test two lon-lat grid specifications for equality.
  !
  FUNCTION lonlat_grid_compare(grid1, grid2)
    LOGICAL :: lonlat_grid_compare
    TYPE(t_lon_lat_grid), INTENT(IN) :: grid1, grid2
    ! local variables
    REAL(wp), PARAMETER :: ZERO_TOL = 1.e-15_wp

    lonlat_grid_compare = .TRUE.
    IF  ( (grid1%total_dim    /= grid2%total_dim)                            .OR.  &
      &   (grid1%nblks        /= grid2%nblks    )                            .OR.  &
      &   (grid1%npromz       /= grid2%npromz   )                            .OR.  &
      &   (grid1%lon_dim      /= grid2%lon_dim  )                            .OR.  &
      &   (grid1%lat_dim      /= grid2%lat_dim  )                            .OR.  &
      &   (grid1%reg_def_mode /= grid2%reg_def_mode  )                       .OR.  &
      &   (.NOT. float_cmp(grid1%reg_lon_def, grid2%reg_lon_def , ZERO_TOL)) .OR.  &
      &   (.NOT. float_cmp(grid1%reg_lat_def, grid2%reg_lat_def , ZERO_TOL)) .OR.  &
      &   (.NOT. float_cmp(grid1%north_pole,  grid2%north_pole  , ZERO_TOL)) .OR.  &
      &   (.NOT. float_cmp(grid1%delta,       grid2%delta       , ZERO_TOL)) .OR.  &
      &   (.NOT. float_cmp(grid1%start_corner,grid2%start_corner, ZERO_TOL)) )  THEN
      lonlat_grid_compare = .FALSE.
    END IF

  END FUNCTION lonlat_grid_compare

END MODULE mo_intp_data_strc
