
! #ifdef __xlC__
! @PROCESS HOT
! #endif
#ifdef __PGI
!pgi$g opt=0
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
!!  Modification by Constantin Junk, MPI-M (2011-05-05)
!!  - moved setup of interpol_ctl variables to namelists/mo_interpol_nml
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
MODULE mo_intp_state
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
USE mo_exception,           ONLY: message, finish
USE mo_impl_constants,      ONLY: SUCCESS, MAX_CHAR_LENGTH, ihs_ocean
USE mo_model_domain,        ONLY: t_patch
USE mo_grid_config,         ONLY: n_dom, n_dom_start, lplane, l_limited_area
USE mo_parallel_config,     ONLY: nproma
USE mo_run_config,          ONLY: ltransport
USE mo_dynamics_config,     ONLY: iequations
USE mo_interpol_config,     ONLY: i_cori_method, rbf_vec_dim_c, rbf_c2grad_dim, &
  &                               rbf_vec_dim_v, rbf_vec_dim_e, lsq_lin_set,    &
  &                               lsq_high_set
USE mo_intp_data_strc,      ONLY: t_int_state
USE mo_intp_rbf_coeffs,     ONLY: rbf_vec_index_cell, rbf_vec_index_edge,                &
  &                               rbf_vec_index_vertex, rbf_vec_compute_coeff_cell,      &
  &                               rbf_vec_compute_coeff_edge,                            &
  &                               rbf_vec_compute_coeff_vertex, rbf_c2grad_index,        &
  &                               rbf_compute_coeff_c2grad
USE mo_intp_coeffs,         ONLY: compute_heli_bra_coeff_idx, init_cellavg_wgt,        &
  &                               init_geo_factors, complete_patchinfo, init_tplane_e, &
  &                               init_tplane_c,   tri_quadrature_pts,                 &
  &                               init_nudgecoeffs, tri_quadrature_pts
  !                               init_geo_factors_oce, par_init_scalar_product_oce
USE mo_intp_coeffs_lsq_bln, ONLY: lsq_stencil_create, lsq_compute_coeff_cell,          &
  &                               scalar_int_coeff, bln_int_coeff_e2c
USE mo_sync,                ONLY: SYNC_C, SYNC_E, SYNC_V
USE mo_communication,       ONLY: t_comm_pattern, blk_no, idx_no, idx_1d, &
  &                               delete_comm_pattern, exchange_data
  USE mo_communication_factory, ONLY: setup_comm_pattern
! USE mo_ocean_nml,           ONLY: idisc_scheme
USE mo_decomposition_tools, ONLY: t_grid_domain_decomp_info, get_valid_local_index
USE mo_dist_dir,            ONLY: dist_dir_get_owners



IMPLICIT NONE

PRIVATE

!!CJPUBLIC :: setup_interpol,
PUBLIC :: construct_2d_interpol_state, destruct_2d_interpol_state
PUBLIC :: transfer_interpol_state
PUBLIC :: allocate_int_state, deallocate_int_state

INTERFACE xfer_var
  MODULE PROCEDURE xfer_var_r2
  MODULE PROCEDURE xfer_var_r3
  MODULE PROCEDURE xfer_var_r4
  MODULE PROCEDURE xfer_var_i2
END INTERFACE

INTERFACE xfer_idx
  MODULE PROCEDURE xfer_idx_2
  MODULE PROCEDURE xfer_idx_3
END INTERFACE

CLASS(t_comm_pattern), POINTER :: comm_pat_glb_to_loc_c, &
  &                               comm_pat_glb_to_loc_e, &
  &                               comm_pat_glb_to_loc_v


CONTAINS

!-------------------------------------------------------------------------
!
!
!> Allocation of components of interpolation state.
!!
!! @par Revision History
!! Split off from construct_2d_interpol_state, Rainer Johanni (2010-10-26)
!!
SUBROUTINE allocate_int_state( ptr_patch, ptr_int)
!
  TYPE(t_patch), INTENT(IN) :: ptr_patch

  TYPE(t_int_state), INTENT(inout) :: ptr_int

  INTEGER :: nblks_c, nblks_e, nblks_v, nincr
  INTEGER :: ist,ie
  INTEGER :: idummy

!-----------------------------------------------------------------------

  !
  ! determine size of arrays, i.e.
  ! values for the blocking
  !
  nblks_c  = ptr_patch%nblks_c
  nblks_e  = ptr_patch%nblks_e
  nblks_v  = ptr_patch%nblks_v
  !
  !
  ! allocate interpolation state
  !
  ! c_lin_e
  !
  ALLOCATE (ptr_int%c_lin_e(nproma,2,nblks_e), STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:construct_int_state',  &
    &            'allocation for c_lin_e failed')
  ENDIF

  IF (ptr_patch%geometry_info%cell_type == 3) THEN
    !
    ! e_bln_c_s
    !
    ALLOCATE (ptr_int%e_bln_c_s(nproma,3,nblks_c), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state', &
      &            'allocation for e_bln_c_s failed')
    ENDIF
    !
    ! e_bln_c_u
    !
    ALLOCATE (ptr_int%e_bln_c_u(nproma,6,nblks_c), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state', &
      &            'allocation for e_bln_c_u failed')
    ENDIF
    !
    ! e_bln_c_v
    !
    ALLOCATE (ptr_int%e_bln_c_v(nproma,6,nblks_c), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state', &
      &            'allocation for e_bln_c_v failed')
    ENDIF
    !
    ! c_bln_avg
    !
    ALLOCATE (ptr_int%c_bln_avg(nproma,4,nblks_c), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state', &
      &            'allocation for c_bln_avg failed')
    ENDIF
    !
    ! gradc_bmat
    !
    ALLOCATE (ptr_int%gradc_bmat(nproma,2,3,nblks_c), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state', &
      &            'allocation for gradc_bmat failed')
    ENDIF
    !
    ! e_flx_avg
    !
    ALLOCATE (ptr_int%e_flx_avg(nproma,5,nblks_e), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state', &
      &            'allocation for e_flx_avg failed')
    ENDIF
    !
    !e_aw_v
    !
    ALLOCATE (ptr_int%e_aw_v(nproma,6,nblks_v), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state', &
      &            'allocation for e_aw_v failed')
    ENDIF

  ENDIF
  !
  ! v_1o2_e
  !
  ALLOCATE (ptr_int%v_1o2_e(nproma,2,nblks_e), STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:construct_int_state', &
    &            'allocation for v_1o2_e failed')
  ENDIF
  !
  ! e_inn_c
  !
  ALLOCATE (ptr_int%e_inn_c(nproma,ptr_patch%geometry_info%cell_type,nblks_c), STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:construct_int_state', &
    &            'allocation for e_inn_c failed')
  ENDIF
  !
  IF (ptr_patch%geometry_info%cell_type == 6) THEN
    !
    ! e_inn_v
    !
    ALLOCATE (ptr_int%e_inn_v(nproma,3,nblks_v), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state', &
      &            'allocation for e_inn_v failed')
    ENDIF
    !
    ! e_aw_c
    !
    ALLOCATE (ptr_int%e_aw_c(nproma,6,nblks_c), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state', &
      &            'allocation for e_aw_c failed')
    ENDIF
    !
    ! r_aw_c
    !
    ALLOCATE (ptr_int%r_aw_c(nproma,6,nblks_c), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state', &
      &            'allocation for r_aw_c failed')
    ENDIF
    !
    ! e_aw_v
    !
    ALLOCATE (ptr_int%e_aw_v(nproma,3,nblks_v), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state', &
      &            'allocation for e_aw_v failed')
    ENDIF
    !
    ! e_1o3_v
    !
    ALLOCATE (ptr_int%e_1o3_v(nproma,3,nblks_v), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state', &
      &            'allocation for e_1o3_v failed')
    ENDIF
    !
    ! tria_aw_rhom
    !
    ALLOCATE (ptr_int%tria_aw_rhom(nproma,2,nblks_e), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state', &
      &            'allocation for tria_aw_rhom failed')
    ENDIF
    !
  ENDIF
  !
  ! verts_aw_cells
  !
  ALLOCATE (ptr_int%verts_aw_cells(nproma,ptr_patch%geometry_info%cell_type,nblks_c), STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:construct_int_state',                     &
      &             'allocation for verts_aw_cells failed')
  ENDIF
  !
  ! cells_aw_verts
  !
  ALLOCATE (ptr_int%cells_aw_verts(nproma,9-ptr_patch%geometry_info%cell_type,nblks_v), STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:construct_int_state',                     &
      &             'allocation for cells_aw_verts failed')
  ENDIF
  !
  ! cells_plwa_verts
  !
  ALLOCATE (ptr_int%cells_plwa_verts(nproma,9-ptr_patch%geometry_info%cell_type,nblks_v), STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:construct_int_state',                     &
      &             'allocation for cells_plwa_verts failed')
  ENDIF
  !
  IF( ptr_patch%geometry_info%cell_type == 6 ) THEN
     !
     ! tria_north
     !
     ALLOCATE (ptr_int%tria_north(3,nproma,nblks_v), STAT=ist )
     IF (ist /= SUCCESS) THEN
       CALL finish ('mo_interpolation:construct_int_state',               &
         &             'allocation for tria_north failed')
     ENDIF
     !
     ! tria_east
     !
     ALLOCATE (ptr_int%tria_east(3,nproma,nblks_v), STAT=ist )
     IF (ist /= SUCCESS) THEN
       CALL finish ('mo_interpolation:construct_int_state',               &
         &             'allocation for tria_east failed')
     ENDIF
     !
     ! hex_north
     !
     ALLOCATE (ptr_int%hex_north(nproma,6,nblks_c), STAT=ist )
     IF (ist /= SUCCESS) THEN
       CALL finish ('mo_interpolation:construct_int_state',               &
         &             'allocation for hex_north failed')
     ENDIF
     !
     ! hex_east
     !
     ALLOCATE (ptr_int%hex_east(nproma,6,nblks_c), STAT=ist )
     IF (ist /= SUCCESS) THEN
       CALL finish ('mo_interpolation:construct_int_state',               &
         &             'allocation for hex_east failed')
     ENDIF
     !
     IF(i_cori_method>=3) THEN
       !
       ! quad_north
       !
       ALLOCATE (ptr_int%quad_north(5,nproma,nblks_e), STAT=ist )
       IF (ist /= SUCCESS) THEN
         CALL finish ('mo_interpolation:construct_int_state',               &
           &             'allocation for quad_north failed')
       ENDIF
       !
       ! quad_east
       !
       ALLOCATE (ptr_int%quad_east(5,nproma,nblks_e), STAT=ist )
       IF (ist /= SUCCESS) THEN
         CALL finish ('mo_interpolation:construct_int_state',               &
           &             'allocation for quad_east failed')
       ENDIF
     ENDIF
     !
     ! cno_en
     !
     ALLOCATE (ptr_int%cno_en(nproma,2,nblks_e), STAT=ist )
     IF (ist /= SUCCESS) THEN
       CALL finish ('mo_interpolation:construct_int_state',               &
         &             'allocation for cno_en failed')
     ENDIF
     !
     ! cea_en
     !
     ALLOCATE (ptr_int%cea_en(nproma,2,nblks_e), STAT=ist )
     IF (ist /= SUCCESS) THEN
       CALL finish ('mo_interpolation:construct_int_state',               &
         &             'allocation for cea_en failed')
     ENDIF
     !
     SELECT CASE (i_cori_method)
     CASE (1,3,4)
       nincr = 14
     CASE (2)
       nincr = 10
     END SELECT
     !
     ! heli_coeff
     !
     ALLOCATE (ptr_int%heli_coeff(nincr, nproma, nblks_e), STAT=ist )
     IF (ist /= SUCCESS) THEN
       CALL finish ('mo_interpolation:construct_int_state',               &
         &             'allocation for heli_coeff failed')
     ENDIF

     IF(i_cori_method<3) THEN
       !
       ! heli_vn_idx
       !
       ALLOCATE (ptr_int%heli_vn_idx(nincr,nproma, nblks_e), STAT=ist )
       IF (ist /= SUCCESS) THEN
         CALL finish ('mo_interpolation:construct_int_state',               &
           &             'allocation for heli_vn_idx failed')
       ENDIF
       !
       ! heli_vn_blk
       !
       ALLOCATE (ptr_int%heli_vn_blk(nincr,nproma, nblks_e), STAT=ist )
       IF (ist /= SUCCESS) THEN
         CALL finish ('mo_interpolation:construct_int_state',               &
           &             'allocation for heli_vn_blk failed')
       ENDIF
     ENDIF

  ENDIF

  IF (ptr_patch%geometry_info%cell_type == 3) THEN
    !
    ! rbf_vec_idx_c, rbf_vec_blk_c
    !
    ALLOCATE (ptr_int%rbf_vec_idx_c(rbf_vec_dim_c, nproma, nblks_c), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',&
      &            'allocation for rbf_vec_idx_c failed')
    ENDIF
    ALLOCATE (ptr_int%rbf_vec_blk_c(rbf_vec_dim_c, nproma, nblks_c), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',&
      &            'allocation for rbf_vec_blk_c failed')
    ENDIF
    !
    ! rbf_vec_stencil_c
    !
    ALLOCATE (ptr_int%rbf_vec_stencil_c(nproma, nblks_c), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',&
      &            'allocation for rbf_vec_stencil_c failed')
    ENDIF
    !
    ! rbf_vec_coeff_c
    !
    ALLOCATE (ptr_int%rbf_vec_coeff_c(rbf_vec_dim_c, 2, nproma, nblks_c), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',&
      &            'allocation for rbf_vec_coeff_c failed')
    ENDIF
    !
    ! rbf_c2grad_idx, rbf_c2grad_blk
    !
    ALLOCATE (ptr_int%rbf_c2grad_idx(rbf_c2grad_dim, nproma, nblks_c), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',&
      &            'allocation for rbf_c2grad_idx failed')
    ENDIF
    ALLOCATE (ptr_int%rbf_c2grad_blk(rbf_c2grad_dim, nproma, nblks_c), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',&
      &            'allocation for rbf_c2grad_blk failed')
    ENDIF
    !
    ! rbf_c2grad_coeff
    !
    ALLOCATE (ptr_int%rbf_c2grad_coeff(rbf_c2grad_dim, 2, nproma, nblks_c), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',&
      &            'allocation for rbf_c2grad_coeff failed')
    ENDIF
    !
    ! rbf_vec_idx_v, rbf_vec_blk_v
    !
    ALLOCATE (ptr_int%rbf_vec_idx_v(rbf_vec_dim_v, nproma, nblks_v), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',&
      &            'allocation for rbf_vec_idx_v failed')
    ENDIF
    ALLOCATE (ptr_int%rbf_vec_blk_v(rbf_vec_dim_v, nproma, nblks_v), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',&
        &          'allocation for rbf_vec_blk_v failed')
    ENDIF
    !
    ! rbf_vec_stencil_v
    !
    ALLOCATE (ptr_int%rbf_vec_stencil_v(nproma, nblks_v), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',&
      &            'allocation for rbf_vec_stencil_v failed')
    ENDIF
    !
    ! rbf_vec_coeff_v
    !
    ALLOCATE (ptr_int%rbf_vec_coeff_v(rbf_vec_dim_v, 2, nproma, nblks_v), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',&
      &            'allocation for rbf_vec_coeff_v failed')
    ENDIF
    !
    ! rbf_vec_idx_e, rbf_vec_blk_e
    !
    ALLOCATE (ptr_int%rbf_vec_idx_e(rbf_vec_dim_e, nproma, nblks_e), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',&
      &            'allocation for rbf_vec_idx_e failed')
    ENDIF
    ALLOCATE (ptr_int%rbf_vec_blk_e(rbf_vec_dim_e, nproma, nblks_e), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',&
      &            'allocation for rbf_vec_blk_e failed')
    ENDIF
    !
    ! rbf_vec_stencil_e
    !
    ALLOCATE (ptr_int%rbf_vec_stencil_e(nproma, nblks_e), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',&
      &            'allocation for rbf_vec_stencil_e failed')
    ENDIF
    !
    ! rbf_vec_coeff_e
    !
    ALLOCATE (ptr_int%rbf_vec_coeff_e(rbf_vec_dim_e, nproma, nblks_e), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',&
      &            'allocation for rbf_vec_coeff_e failed')
    ENDIF

  ENDIF

  IF( ltransport .OR. iequations == 3) THEN
    !
    ! pos_on_tplane_e
    !
    ALLOCATE (ptr_int%pos_on_tplane_e(nproma, nblks_e, 8, 2), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',&
      &            'allocation for pos_on_tplane_e failed')
    ENDIF
    !
    ! tplane_e_dotprod
    !
    ALLOCATE (ptr_int%tplane_e_dotprod(nproma, nblks_e, 4, 4), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',&
      &            'allocation for tplane_e_dotprod failed')
    ENDIF

    IF (ptr_patch%geometry_info%cell_type == 3) THEN
      !
      ! pos_on_tplane_c_edge
      !
      ALLOCATE (ptr_int%pos_on_tplane_c_edge(nproma, nblks_e, 2, 5), STAT=ist )
      IF (ist /= SUCCESS) THEN
        CALL finish ('mo_interpolation:construct_int_state',&
        &            'allocation for pos_on_tplane_c_edge failed')
      ENDIF
    ENDIF

    !
    ! Least squares reconstruction
    !
    ! *** linear ***
    !
    !
    ! lsq_dim_stencil
    !
    ALLOCATE (ptr_int%lsq_lin%lsq_dim_stencil(nproma, nblks_c), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',&
        &          'allocation for lsq_dim_stencil failed')
    ENDIF
    !
    ! lsq_idx_c
    !
    ALLOCATE (ptr_int%lsq_lin%lsq_idx_c(nproma, nblks_c, lsq_lin_set%dim_c),          &
      &       STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',&
        &          'allocation for lsq_idx_c failed')
    ENDIF
    !
    ! lsq_blk_c
    !
    ALLOCATE (ptr_int%lsq_lin%lsq_blk_c(nproma, nblks_c, lsq_lin_set%dim_c),          &
      &       STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',&
        &          'allocation for lsq_blk_c failed')
    ENDIF
    !
    ! lsq_weights_c
    !
    ALLOCATE (ptr_int%lsq_lin%lsq_weights_c(nproma, lsq_lin_set%dim_c, nblks_c),      &
      &       STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',&
        &          'allocation for lsq_weights_c failed')
    ENDIF
    !
    ! lsq_qtmat_c
    !
    ALLOCATE (ptr_int%lsq_lin%lsq_qtmat_c(nproma, lsq_lin_set%dim_unk, lsq_lin_set%dim_c, &
      &       nblks_c), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',&
        &          'allocation for lsq_qtmat_c failed')
    ENDIF
    !
    ! lsq_rmat_rdiag_c
    !
    ALLOCATE (ptr_int%lsq_lin%lsq_rmat_rdiag_c(nproma, lsq_lin_set%dim_unk, nblks_c), &
      &       STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',&
        &          'allocation for lsq_rmat_rdiag_c failed')
    ENDIF
    !
    ! lsq_rmat_utri_c
    !
    idummy=(lsq_lin_set%dim_unk*lsq_lin_set%dim_unk - lsq_lin_set%dim_unk)/2
    ALLOCATE (ptr_int%lsq_lin%lsq_rmat_utri_c(nproma, idummy, nblks_c), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',&
        &          'allocation for lsq_rmat_utri_c failed')
    ENDIF
    !
    ! lsq_pseudoinv
    !
    ALLOCATE (ptr_int%lsq_lin%lsq_pseudoinv(nproma, lsq_lin_set%dim_unk,              &
      &       lsq_lin_set%dim_c, nblks_c), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',                            &
        &             'allocation for lsq_pseudoinv failed')
    ENDIF
    !
    ! lsq_moments
    !
    ALLOCATE (ptr_int%lsq_lin%lsq_moments(nproma, nblks_c, lsq_lin_set%dim_unk),      &
      &       STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',                            &
        &             'allocation for lsq_moments failed')
    ENDIF
    !
    ! lsq_moments_hat
    !
    ALLOCATE (ptr_int%lsq_lin%lsq_moments_hat(nproma, nblks_c, lsq_lin_set%dim_c,     &
      &       lsq_lin_set%dim_unk), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',                            &
        &             'allocation for lsq_moments_hat failed')
    ENDIF

    ! *** higher order ***
    !
    !
    ! lsq_dim_stencil
    !
    ALLOCATE (ptr_int%lsq_high%lsq_dim_stencil(nproma, nblks_c), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',&
        &          'allocation for lsq_dim_stencil failed')
    ENDIF
    !
    ! lsq_idx_c
    !
    ALLOCATE (ptr_int%lsq_high%lsq_idx_c(nproma, nblks_c, lsq_high_set%dim_c),        &
      &       STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',&
        &          'allocation for lsq_idx_c failed')
    ENDIF
    !
    ! lsq_blk_c
    !
    ALLOCATE (ptr_int%lsq_high%lsq_blk_c(nproma, nblks_c, lsq_high_set%dim_c),        &
      &       STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',&
        &          'allocation for lsq_blk_c failed')
    ENDIF
    !
    ! lsq_weights_c
    !
    ALLOCATE (ptr_int%lsq_high%lsq_weights_c(nproma, lsq_high_set%dim_c, nblks_c),    &
      &       STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',&
        &          'allocation for lsq_weights_c failed')
    ENDIF
    !
    ! lsq_qtmat_c
    !
    ALLOCATE (ptr_int%lsq_high%lsq_qtmat_c(nproma, lsq_high_set%dim_unk, lsq_high_set%dim_c, &
      &       nblks_c), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',&
        &          'allocation for lsq_qtmat_c failed')
    ENDIF
    !
    ! lsq_rmat_rdiag_c
    !
    ALLOCATE (ptr_int%lsq_high%lsq_rmat_rdiag_c(nproma, lsq_high_set%dim_unk, nblks_c), &
      &       STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',&
        &          'allocation for lsq_rmat_rdiag_c failed')
    ENDIF
    !
    ! lsq_rmat_utri_c
    !
    idummy=(lsq_high_set%dim_unk*lsq_high_set%dim_unk - lsq_high_set%dim_unk)/2
    ALLOCATE (ptr_int%lsq_high%lsq_rmat_utri_c(nproma, idummy, nblks_c), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',&
        &          'allocation for lsq_rmat_utri_c failed')
    ENDIF
    !
    ! lsq_pseudoinv
    !
    ALLOCATE (ptr_int%lsq_high%lsq_pseudoinv(nproma, lsq_high_set%dim_unk,            &
      &       lsq_high_set%dim_c, nblks_c), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',                            &
        &             'allocation for lsq_pseudoinv failed')
    ENDIF
    !
    ! lsq_moments
    !
    ALLOCATE (ptr_int%lsq_high%lsq_moments(nproma, nblks_c, lsq_high_set%dim_unk),    &
      &       STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',                            &
        &             'allocation for lsq_moments failed')
    ENDIF
    !
    ! lsq_moments_hat
    !
    ALLOCATE (ptr_int%lsq_high%lsq_moments_hat(nproma, nblks_c, lsq_high_set%dim_c,   &
      &       lsq_high_set%dim_unk), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',                            &
        &             'allocation for lsq_moments_hat failed')
    ENDIF

  ELSE

    !
    ! Least squares reconstruction
    !

    ! *** lin ***
    !
    !
    ! lsq_dim_stencil
    !
    ALLOCATE (ptr_int%lsq_lin%lsq_dim_stencil(0, 0), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',&
        &          'allocation for lsq_dim_stencil failed')
    ENDIF
    !
    ! lsq_idx_c
    !
    ALLOCATE (ptr_int%lsq_lin%lsq_idx_c(0, 0, 0), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',&
        &          'allocation for lsq_idx_c failed')
    ENDIF
    !
    ! lsq_blk_c
    !
    ALLOCATE (ptr_int%lsq_lin%lsq_blk_c(0, 0, 0), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',&
        &          'allocation for lsq_blk_c failed')
    ENDIF
    !
    ! lsq_weights_c
    !
    ALLOCATE (ptr_int%lsq_lin%lsq_weights_c(0, 0, 0), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',&
        &          'allocation for lsq_weights_c failed')
    ENDIF
    !
    ! lsq_qtmat_c
    !
    ALLOCATE (ptr_int%lsq_lin%lsq_qtmat_c(0, 0, 0, 0), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',&
        &          'allocation for lsq_qtmat_c failed')
    ENDIF
    !
    ! lsq_rmat_rdiag_c
    !
    ALLOCATE (ptr_int%lsq_lin%lsq_rmat_rdiag_c(0, 0, 0), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',&
        &          'allocation for lsq_rmat_rdiag_c failed')
    ENDIF
    !
    ! lsq_rmat_utri_c
    !
    ALLOCATE (ptr_int%lsq_lin%lsq_rmat_utri_c(0, 0, 0), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',&
        &          'allocation for lsq_rmat_utri_c failed')
    ENDIF
    !
    ! lsq_pseudoinv
    !
    ALLOCATE (ptr_int%lsq_lin%lsq_pseudoinv(0, 0, 0, 0), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',                       &
        &             'allocation for lsq_pseudoinv failed')
    ENDIF
    !
    ! lsq_moments
    !
    ALLOCATE (ptr_int%lsq_lin%lsq_moments(0, 0, 0), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',                       &
        &             'allocation for lsq_moments failed')
    ENDIF
    !
    ! lsq_moments_hat
    !
    ALLOCATE (ptr_int%lsq_lin%lsq_moments_hat(0, 0, 0, 0), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',                       &
        &             'allocation for lsq_moments_hat failed')
    ENDIF

    ! *** higher order ***
    !
    !
    ! lsq_dim_stencil
    !
    ALLOCATE (ptr_int%lsq_high%lsq_dim_stencil(0, 0), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',&
        &          'allocation for lsq_dim_stencil failed')
    ENDIF
    !
    ! lsq_idx_c
    !
    ALLOCATE (ptr_int%lsq_high%lsq_idx_c(0, 0, 0), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',&
        &          'allocation for lsq_idx_c failed')
    ENDIF
    !
    ! lsq_blk_c
    !
    ALLOCATE (ptr_int%lsq_high%lsq_blk_c(0, 0, 0), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',&
        &          'allocation for lsq_blk_c failed')
    ENDIF
    !
    ! lsq_weights_c
    !
    ALLOCATE (ptr_int%lsq_high%lsq_weights_c(0, 0, 0), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',&
        &          'allocation for lsq_weights_c failed')
    ENDIF
    !
    ! lsq_qtmat_c
    !
    ALLOCATE (ptr_int%lsq_high%lsq_qtmat_c(0, 0, 0, 0), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',&
        &          'allocation for lsq_qtmat_c failed')
    ENDIF
    !
    ! lsq_rmat_rdiag_c
    !
    ALLOCATE (ptr_int%lsq_high%lsq_rmat_rdiag_c(0, 0, 0), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',&
        &          'allocation for lsq_rmat_rdiag_c failed')
    ENDIF
    !
    ! lsq_rmat_utri_c
    !
    ALLOCATE (ptr_int%lsq_high%lsq_rmat_utri_c(0, 0, 0), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',&
        &          'allocation for lsq_rmat_utri_c failed')
    ENDIF
    !
    ! lsq_pseudoinv
    !
    ALLOCATE (ptr_int%lsq_high%lsq_pseudoinv(0, 0, 0, 0), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',                     &
        &             'allocation for lsq_pseudoinv failed')
    ENDIF
    !
    ! lsq_moments
    !
    ALLOCATE (ptr_int%lsq_high%lsq_moments(0, 0, 0), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',                     &
        &             'allocation for lsq_moments failed')
    ENDIF
    !
    ! lsq_moments_hat
    !
    ALLOCATE (ptr_int%lsq_high%lsq_moments_hat(0, 0, 0, 0), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',                     &
        &             'allocation for lsq_moments_hat failed')
    ENDIF

  END IF

  IF (ptr_patch%geometry_info%cell_type == 3) THEN
    ALLOCATE (ptr_int%geofac_qdiv(nproma, 4, nblks_e), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',                     &
        &             'allocation for geofac_qdiv failed')
    ENDIF
    ALLOCATE (ptr_int%geofac_grdiv(nproma, 5, nblks_e), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',                     &
        &             'allocation for geofac_grdiv failed')
    ENDIF
    ALLOCATE (ptr_int%nudgecoeff_c(nproma, nblks_c), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',                     &
        &             'allocation for nudgecoeff_c failed')
    ENDIF
    ALLOCATE (ptr_int%nudgecoeff_e(nproma, nblks_e), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',                     &
        &             'allocation for nudgecoeff_e failed')
    ENDIF

    !
    ! Quadrature points and weights for integration over triangular element
    !
    ALLOCATE (ptr_int%gquad%qpts_tri_l(nproma, nblks_c), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',                     &
        &             'allocation for qpts_tri_l failed')
    ENDIF
    ALLOCATE (ptr_int%gquad%qpts_tri_q(nproma, nblks_c,3), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',                     &
        &             'allocation for qpts_tri_q failed')
    ENDIF
    ALLOCATE (ptr_int%gquad%qpts_tri_c(nproma, nblks_c,4), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',                     &
        &             'allocation for qpts_tri_c failed')
    ENDIF
    ALLOCATE (ptr_int%gquad%weights_tri_q(3), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',                     &
        &             'allocation for weights_tri_q failed')
    ENDIF
    ALLOCATE (ptr_int%gquad%weights_tri_c(4), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',                     &
        &             'allocation for weights_tri_c failed')
    ENDIF
  ENDIF

  ALLOCATE (ptr_int%geofac_div(nproma, ptr_patch%geometry_info%cell_type, nblks_c), STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:construct_int_state',                       &
      &             'allocation for geofac_div failed')
  ENDIF

  ALLOCATE (ptr_int%geofac_rot(nproma, 9-ptr_patch%geometry_info%cell_type, nblks_v), STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:construct_int_state', &
    &             'allocation for geofac_rot failed')
  ENDIF

  ALLOCATE (ptr_int%geofac_n2s(nproma, ptr_patch%geometry_info%cell_type+1, nblks_c), STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:construct_int_state',                       &
      &             'allocation for geofac_n2s failed')
  ENDIF

  ALLOCATE (ptr_int%geofac_grg(nproma, ptr_patch%geometry_info%cell_type+1, nblks_c, 2), STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:construct_int_state',                       &
      &             'allocation for geofac_grg failed')
  ENDIF

  ALLOCATE (ptr_int%primal_normal_ec(nproma, nblks_c,ptr_patch%geometry_info%cell_type, 2), STAT=ist)
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:construct_int_state',                       &
      &             'allocation for primal_normal_ec failed')
  ENDIF

  ALLOCATE (ptr_int%edge_cell_length(nproma, nblks_c, ptr_patch%geometry_info%cell_type), STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:construct_int_state',                       &
      &             'allocation for edge_cell_length failed')
  ENDIF

  ALLOCATE (ptr_int%cell_vert_dist(nproma, ptr_patch%geometry_info%cell_type, 2, nblks_c), STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:construct_int_state',                       &
      &             'allocation for cell_vert_dist failed')
  ENDIF

  IF (ptr_patch%geometry_info%cell_type == 6) THEN

    ALLOCATE (ptr_int%dir_gradh_i1(6, nproma, nblks_e), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',                       &
        &             'allocation for dir_gradh_i1 failed')
    ENDIF
    ALLOCATE (ptr_int%dir_gradh_i2(6, nproma, nblks_e), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',                       &
        &             'allocation for dir_gradh_i2 failed')
    ENDIF
    ALLOCATE (ptr_int%dir_gradh_b1(6, nproma, nblks_e), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',                       &
        &             'allocation for dir_gradh_b1 failed')
    ENDIF
    ALLOCATE (ptr_int%dir_gradh_b2(6, nproma, nblks_e), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',                       &
        &             'allocation for dir_gradh_b2 failed')
    ENDIF
    ALLOCATE (ptr_int%dir_gradhux_c1(6, nproma, nblks_e), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',                       &
        &             'allocation for dir_gradhux_c1 failed')
    ENDIF
    ALLOCATE (ptr_int%dir_gradhux_c2(6, nproma, nblks_e), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',                       &
        &             'allocation for dir_gradhux_c2 failed')
    ENDIF
    ALLOCATE (ptr_int%strain_def_c1(6, nproma, nblks_e), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',                       &
        &             'allocation for strain_def_c1 failed')
    ENDIF
    ALLOCATE (ptr_int%strain_def_c2(6, nproma, nblks_e), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',                       &
        &             'allocation for strain_def_c2 failed')
    ENDIF

    ALLOCATE (ptr_int%dir_gradt_i1(9, nproma, nblks_e), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',                       &
        &             'allocation for dir_gradt_i1 failed')
    ENDIF
    ALLOCATE (ptr_int%dir_gradt_i2(9, nproma, nblks_e), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',                       &
        &             'allocation for dir_gradt_i2 failed')
    ENDIF
    ALLOCATE (ptr_int%dir_gradt_b1(9, nproma, nblks_e), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',                       &
        &             'allocation for dir_gradt_b1 failed')
    ENDIF
    ALLOCATE (ptr_int%dir_gradt_b2(9, nproma, nblks_e), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',                       &
        &             'allocation for dir_gradt_b2 failed')
    ENDIF
    ALLOCATE (ptr_int%dir_gradtxy_v1(9, nproma, nblks_e), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',                       &
        &             'allocation for dir_gradtxy_v1 failed')
    ENDIF
    ALLOCATE (ptr_int%dir_gradtxy_v2(9, nproma, nblks_e), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',                       &
        &             'allocation for dir_gradtxy_v2 failed')
    ENDIF
    ALLOCATE (ptr_int%dir_gradtyx_v1(9, nproma, nblks_e), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',                       &
        &             'allocation for dir_gradtyx_v1 failed')
    ENDIF
    ALLOCATE (ptr_int%dir_gradtyx_v2(9, nproma, nblks_e), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',                       &
        &             'allocation for dir_gradtyx_v2 failed')
    ENDIF
    ALLOCATE (ptr_int%shear_def_v1(9, nproma, nblks_e), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',                       &
        &             'allocation for shear_def_v1 failed')
    ENDIF
    ALLOCATE (ptr_int%shear_def_v2(9, nproma, nblks_e), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',                       &
        &             'allocation for shear_def_v2 failed')
    ENDIF
  ENDIF

  IF ( iequations == ihs_ocean) THEN
    !
    ! arrays that are required for #slo OLD# reconstruction
    !
    ALLOCATE(ptr_int%dist_cell2edge(nproma,nblks_e,2),STAT=ist)
    IF (ist /= SUCCESS) THEN
      CALL finish ('allocating dist_cell2edge failed')
    ENDIF

    !
    ! arrays that are required for setting up the scalar product
    !
    !coefficients for edge to cell mapping, one half of the scalar product.
    !Dimension: nproma,nblks_c encode number of cells, 1:3 corresponds to number
    !of edges per cell, 1:2 is for u and v component of cell vector
    !     ALLOCATE(ptr_int%edge2cell_coeff(nproma,nblks_c,1:3, 1:2),STAT=ist)
    !     IF (ist /= SUCCESS) THEN
    !       CALL finish ('allocating edge2cell_coeff failed')
    !     ENDIF
    ALLOCATE(ptr_int%edge2cell_coeff_cc(nproma,nblks_c,1:3),STAT=ist)
    IF (ist /= SUCCESS) THEN
      CALL finish ('allocating edge2cell_coeff_cc failed')
    ENDIF

    !coefficients for transposed of edge to cell mapping, second half of the scalar product.
    !Dimension: nproma,nblks_e encode number of edges, 1:2 is for cell neighbors of an edge
    ALLOCATE(ptr_int%edge2cell_coeff_cc_t(nproma,nblks_e,1:2),STAT=ist)
    IF (ist /= SUCCESS) THEN
      CALL finish ('allocating transposed edge2cell_coeff failed')
    ENDIF

    !
    !coefficients for edge to vertex mapping.
    !
    !Dimension: nproma,nblks_v encode number of vertices,
    !1:6 is number of edges of a vertex,
    !1:2 is for u and v component of vertex vector
    ALLOCATE(ptr_int%edge2vert_coeff_cc(nproma,nblks_v,1:6),STAT=ist)
    IF (ist /= SUCCESS) THEN
      CALL finish ('allocating edge2vert_coeff failed')
    ENDIF

    ALLOCATE(ptr_int%edge2vert_coeff_cc_t(nproma,nblks_e,1:2),STAT=ist)
    IF (ist /= SUCCESS) THEN
      CALL finish ('allocating edge2vert_coeff failed')
    ENDIF
    ALLOCATE(ptr_int%edge2vert_vector_cc(nproma,nblks_v,1:6),STAT=ist)
    IF (ist /= SUCCESS) THEN
      CALL finish ('allocating edge2vert_vector failed')
    ENDIF

    !
    !normalizing factors for edge to cell mapping.
    !
    !Either by fixed volume or by variable one taking the surface elevation
    !into account. The later one depends on time and space.
    ALLOCATE(ptr_int%fixed_vol_norm(nproma,nblks_c),STAT=ist)
    IF (ist /= SUCCESS) THEN
      CALL finish ('allocating fixed_vol_norm failed')
    ENDIF
    ALLOCATE(ptr_int%variable_vol_norm(nproma,nblks_c,1:3),STAT=ist)
    IF (ist /= SUCCESS) THEN
      CALL finish ('allocating variable_vol_norm failed')
    ENDIF

    ALLOCATE(ptr_int%variable_dual_vol_norm(nproma,nblks_v,1:6),STAT=ist)
    IF (ist /= SUCCESS) THEN
      CALL finish ('allocating variable_dual_vol_norm failed')
    ENDIF
    DO ie = 1,3
      ptr_int%edge2cell_coeff_cc%x(ie)   = 0._wp
      ptr_int%edge2cell_coeff_cc_t%x(ie) = 0._wp
      ptr_int%edge2vert_coeff_cc%x(ie)   = 0._wp
      ptr_int%edge2vert_coeff_cc_t%x(ie) = 0._wp
      ptr_int%edge2vert_vector_cc%x(ie)  = 0._wp
    END DO

    ptr_int%fixed_vol_norm         = 0._wp
    ptr_int%variable_vol_norm      = 0._wp
    ptr_int%variable_dual_vol_norm = 0._wp

    ptr_int%dist_cell2edge = 0._wp
  ENDIF

  !
  ! initialize all components
  !
  IF (ptr_patch%geometry_info%cell_type == 3 ) THEN
    ptr_int%e_bln_c_s     = 0._wp
    ptr_int%e_bln_c_u     = 0._wp
    ptr_int%e_bln_c_v     = 0._wp
    ptr_int%c_bln_avg     = 0._wp
    ptr_int%gradc_bmat    = 0._wp
    ptr_int%e_flx_avg     = 0._wp
    ptr_int%e_aw_v        = 0._wp
  ENDIF

  ptr_int%v_1o2_e          = 0.5_wp
  ptr_int%c_lin_e          = 0._wp
  ptr_int%e_inn_c          = 0._wp
  ptr_int%verts_aw_cells   = 0._wp
  ptr_int%cells_aw_verts   = 0._wp
  ptr_int%cells_plwa_verts = 0._wp

  IF (ptr_patch%geometry_info%cell_type == 6 ) THEN
    ptr_int%e_inn_v     = 0._wp
    ptr_int%tria_aw_rhom= 0._wp
    ptr_int%e_aw_c      = 0._wp
    ptr_int%r_aw_c      = 0._wp
    ptr_int%e_aw_v      = 0._wp
    ptr_int%e_1o3_v     = 1.0_wp/3.0_wp
    ptr_int%hex_north   = 0.0_wp
    ptr_int%hex_east    = 0.0_wp
    ptr_int%tria_north  = 0.0_wp
    ptr_int%tria_east   = 0.0_wp
    IF(i_cori_method>=3)THEN
      ptr_int%quad_north  = 0.0_wp
      ptr_int%quad_east   = 0.0_wp
    ENDIF
    ptr_int%cno_en      = 0.0_wp
    ptr_int%cea_en      = 0.0_wp
  ENDIF

  IF( ptr_patch%geometry_info%cell_type == 3) THEN
    ptr_int%rbf_vec_idx_c     = 0
    ptr_int%rbf_vec_blk_c     = 0
    ptr_int%rbf_vec_stencil_c = 0
    ptr_int%rbf_vec_coeff_c   = 0._wp

    ptr_int%rbf_c2grad_idx    = 0
    ptr_int%rbf_c2grad_blk    = 0
    ptr_int%rbf_c2grad_coeff  = 0._wp

    ptr_int%rbf_vec_idx_v     = 0
    ptr_int%rbf_vec_blk_v     = 0
    ptr_int%rbf_vec_stencil_v = 0
    ptr_int%rbf_vec_coeff_v   = 0._wp

    ptr_int%rbf_vec_idx_e     = 0
    ptr_int%rbf_vec_blk_e     = 0
    ptr_int%rbf_vec_stencil_e = 0
    ptr_int%rbf_vec_coeff_e   = 0._wp

  ENDIF

  IF( ptr_patch%geometry_info%cell_type == 6 ) THEN
    ptr_int%heli_coeff        = 0._wp
    IF (i_cori_method < 3) THEN
      ptr_int%heli_vn_idx       = 0
      ptr_int%heli_vn_blk       = 0
    ENDIF
  ENDIF

  IF( ltransport .OR. iequations == 3) THEN

    ptr_int%pos_on_tplane_e           = 0._wp
    ptr_int%tplane_e_dotprod          = 0._wp

    IF (ptr_patch%geometry_info%cell_type == 3) THEN
      ptr_int%pos_on_tplane_c_edge(:,:,:,:)%lon = 0._wp
      ptr_int%pos_on_tplane_c_edge(:,:,:,:)%lat = 0._wp
    ENDIF

    ptr_int%lsq_lin%lsq_dim_stencil   = 0
    ptr_int%lsq_lin%lsq_idx_c         = 0
    ptr_int%lsq_lin%lsq_blk_c         = 0
    ptr_int%lsq_lin%lsq_weights_c     = 0._wp
    ptr_int%lsq_lin%lsq_qtmat_c       = 0._wp
    ptr_int%lsq_lin%lsq_rmat_rdiag_c  = 0._wp
    ptr_int%lsq_lin%lsq_rmat_utri_c   = 0._wp
    ptr_int%lsq_lin%lsq_pseudoinv     = 0._wp
    ptr_int%lsq_lin%lsq_moments       = 0._wp
    ptr_int%lsq_lin%lsq_moments_hat   = 0._wp

    ptr_int%lsq_high%lsq_dim_stencil  = 0
    ptr_int%lsq_high%lsq_idx_c        = 0
    ptr_int%lsq_high%lsq_blk_c        = 0
    ptr_int%lsq_high%lsq_weights_c    = 0._wp
    ptr_int%lsq_high%lsq_qtmat_c      = 0._wp
    ptr_int%lsq_high%lsq_rmat_rdiag_c = 0._wp
    ptr_int%lsq_high%lsq_rmat_utri_c  = 0._wp
    ptr_int%lsq_high%lsq_pseudoinv    = 0._wp
    ptr_int%lsq_high%lsq_moments      = 0._wp
    ptr_int%lsq_high%lsq_moments_hat  = 0._wp
  END IF

  IF (ptr_patch%geometry_info%cell_type ==3) THEN
    ptr_int%geofac_qdiv = 0._wp
    ptr_int%geofac_grdiv = 0._wp
    ptr_int%nudgecoeff_c = 0._wp
    ptr_int%nudgecoeff_e = 0._wp

    ptr_int%gquad%qpts_tri_l(:,:)%lat    = 0._wp
    ptr_int%gquad%qpts_tri_l(:,:)%lon    = 0._wp
    ptr_int%gquad%qpts_tri_q(:,:,:)%lat  = 0._wp
    ptr_int%gquad%qpts_tri_q(:,:,:)%lon  = 0._wp
    ptr_int%gquad%qpts_tri_c(:,:,:)%lat  = 0._wp
    ptr_int%gquad%qpts_tri_c(:,:,:)%lon  = 0._wp
    ptr_int%gquad%weights_tri_q(:)       = 0._wp
    ptr_int%gquad%weights_tri_c(:)       = 0._wp
  ENDIF
  ptr_int%geofac_div = 0._wp
  ptr_int%geofac_rot = 0._wp
  ptr_int%geofac_n2s = 0._wp
  ptr_int%geofac_grg = 0._wp
  ptr_int%primal_normal_ec = 0._wp
  ptr_int%edge_cell_length = 0._wp
  ptr_int%cell_vert_dist = 0._wp

  IF(ptr_patch%geometry_info%cell_type==6) THEN
    ptr_int%dir_gradh_i1 = 0
    ptr_int%dir_gradh_i2 = 0
    ptr_int%dir_gradh_b1 = 0
    ptr_int%dir_gradh_b2 = 0
    ptr_int%dir_gradhux_c1 = 0._wp
    ptr_int%dir_gradhux_c2 = 0._wp
    ptr_int%strain_def_c1 = 0._wp
    ptr_int%strain_def_c2 = 0._wp
    ptr_int%dir_gradt_i1 = 0
    ptr_int%dir_gradt_i2 = 0
    ptr_int%dir_gradt_b1 = 0
    ptr_int%dir_gradt_b2 = 0
    ptr_int%dir_gradtxy_v1 = 0._wp
    ptr_int%dir_gradtxy_v2 = 0._wp
    ptr_int%dir_gradtyx_v1 = 0._wp
    ptr_int%dir_gradtyx_v2 = 0._wp
    ptr_int%shear_def_v1 = 0._wp
    ptr_int%shear_def_v2 = 0._wp
  ENDIF

  CALL message ('mo_intp_state:allocate_int_state','memory allocation finished')

END SUBROUTINE allocate_int_state


!-------------------------------------------------------------------------
!
!
!>
!!               Allocation of components of interpolation state.
!!
!!               Initialization of components.
!!
!! @par Revision History
!! Developed  by  Peter Korn, MPI-M (2005).
!! Modified by Luca Bonaventura, MPI-M (2005).
!! Modified by Luca Bonaventura, Polimi, and Peter Korn, MPI-M (2006-07-28)
!! to construct vector rbf structures
!! Modified by Rainer Johanni (2010-10-26): split out allocation
!!
SUBROUTINE construct_2d_interpol_state(ptr_patch, ptr_int_state)
!
TYPE(t_patch), INTENT(INOUT) :: ptr_patch(n_dom_start:)

TYPE(t_int_state),     INTENT(INOUT) :: ptr_int_state(n_dom_start:)

INTEGER :: jg

CHARACTER(len=MAX_CHAR_LENGTH) :: text

  !-----------------------------------------------------------------------

CALL message('mo_intp_state:construct_2d_interpol_state','start to construct int_state')

DO jg = n_dom_start, n_dom

  WRITE(text,'(a,i0)') 'constructing int_state for patch ',jg
  CALL message('mo_intp_state:construct_2d_interpol_state',text)

  CALL allocate_int_state( ptr_patch(jg), ptr_int_state(jg))

  !
  ! initializion of coefficients for averaging of scalars and kinetic energy
  !
  CALL scalar_int_coeff(ptr_patch(jg), ptr_int_state(jg))

  ! Initialization of coefficients for bilinear cell averaging, divergence, rotation
  ! and nabla_2_scalar; transformation of edge orientation vectors to the locations
  ! of cells and vertices; computation of coefficients for bilinear edge-to-cell
  ! interpolation

  CALL complete_patchinfo( ptr_patch(jg), ptr_int_state(jg))
  CALL init_geo_factors(ptr_patch(jg), ptr_int_state(jg))
  IF (ptr_patch(jg)%geometry_info%cell_type==3)THEN
    CALL init_cellavg_wgt(ptr_patch(jg), ptr_int_state(jg))
    CALL bln_int_coeff_e2c( ptr_patch(jg), ptr_int_state(jg) )
    IF(jg>0) THEN
      IF (l_limited_area .AND. jg == 1 .OR. jg > 1) THEN
        CALL init_nudgecoeffs( ptr_patch(jg), ptr_int_state(jg) )
      ENDIF
    ENDIF
  ENDIF

  !
  ! initialization of indices and coefficients for vector rbf interpolation
  !
  IF (ptr_patch(jg)%geometry_info%cell_type == 3) THEN

    ! ... at cell centers
    CALL rbf_vec_index_cell (ptr_patch(jg), ptr_int_state(jg))
    !
    CALL rbf_vec_compute_coeff_cell (ptr_patch(jg), ptr_int_state(jg))
    !
    ! ... at triangle vertices
    CALL rbf_vec_index_vertex (ptr_patch(jg), ptr_int_state(jg))
    !
    CALL rbf_vec_compute_coeff_vertex (ptr_patch(jg), ptr_int_state(jg))
    !
    ! ... at edge midpoints
    CALL rbf_vec_index_edge (ptr_patch(jg), ptr_int_state(jg))
    !
    CALL rbf_vec_compute_coeff_edge (ptr_patch(jg), ptr_int_state(jg))
    !
    ! Compute coefficients needed for gradient reconstruction at cell midpoints
    CALL rbf_c2grad_index (ptr_patch(jg), ptr_int_state(jg))
    CALL rbf_compute_coeff_c2grad (ptr_patch(jg), ptr_int_state(jg))

    ! initialization of quadrature points and weights
    !
    CALL tri_quadrature_pts (ptr_patch(jg), ptr_int_state(jg))
  ENDIF
  !
  ! initialization of coeffs for hexagon
  !
  IF (ptr_patch(jg)%geometry_info%cell_type==6) THEN
    CALL compute_heli_bra_coeff_idx(ptr_patch(jg), ptr_int_state(jg))
  ENDIF
  !
  ! - Initialization of tangential plane (at edge midpoints) for calculation
  !   of backward trajectories.
  ! - Initialization of tangential plane (at cell centers) - for triangular
  !   grid only
  ! - stencil generation
  ! - initialization of coefficients for least squares gradient
  ! reconstruction at cell centers
  !
  IF ( (ltransport .OR. iequations == 3) .AND. (.NOT. lplane)) THEN

    CALL init_tplane_e(ptr_patch(jg), ptr_int_state(jg))

    IF (ptr_patch(jg)%geometry_info%cell_type==3) THEN
      !
      CALL init_tplane_c(ptr_patch(jg), ptr_int_state(jg))

      CALL lsq_stencil_create( ptr_patch(jg), ptr_int_state(jg)%lsq_lin,      &
        &                      lsq_lin_set%dim_c )
      CALL lsq_compute_coeff_cell( ptr_patch(jg), ptr_int_state(jg)%lsq_lin,  &
        &                      lsq_lin_set%l_consv, lsq_lin_set%dim_c,        &
        &                      lsq_lin_set%dim_unk, lsq_lin_set%wgt_exp )
    ENDIF

    CALL lsq_stencil_create( ptr_patch(jg), ptr_int_state(jg)%lsq_high,     &
      &                   lsq_high_set%dim_c )
    CALL lsq_compute_coeff_cell( ptr_patch(jg), ptr_int_state(jg)%lsq_high, &
      &                       lsq_high_set%l_consv, lsq_high_set%dim_c,     &
      &                       lsq_high_set%dim_unk, lsq_high_set%wgt_exp )
  ENDIF

!  IF ( iequations == ihs_ocean) THEN
!    IF (idisc_scheme==1) THEN
!      CALL par_init_scalar_product_oce(ptr_patch(jg), ptr_int_state(jg))
!    ENDIF
!    CALL init_geo_factors_oce(ptr_patch(jg), ptr_int_state(jg))
!  ENDIF
ENDDO

CALL message('mo_intp_state:construct_2d_interpol_state', &
  & 'construction of interpolation state finished')

END SUBROUTINE construct_2d_interpol_state

!-------------------------------------------------------------------------
!
!> xfer_var family: transfer variables from parent to local parent
!!
!! @par Revision History
!! Developed  by  Rainer Johanni (2011-10-26)
!-------------------------------------------------------------------------

SUBROUTINE xfer_var_r2(typ, pos_nproma, pos_nblks, p_p, p_lp, arri, arro)

  INTEGER, INTENT(IN) :: typ, pos_nproma, pos_nblks
  TYPE(t_patch), INTENT(IN) :: p_p, p_lp
  REAL(wp), INTENT(IN)    :: arri(:,:)
  REAL(wp), INTENT(INOUT) :: arro(:,:)
  ! local variables

  IF(typ == SYNC_C) THEN
    CALL exchange_data(comm_pat_glb_to_loc_c, RECV=arro, SEND=arri)
  ELSEIF(typ == SYNC_E) THEN
    CALL exchange_data(comm_pat_glb_to_loc_e, RECV=arro, SEND=arri)
  ELSEIF(typ == SYNC_V) THEN
    CALL exchange_data(comm_pat_glb_to_loc_v, RECV=arro, SEND=arri)
  ELSE
    CALL finish ('mo_interpolation:xfer_var','Illegal type for sync')
  ENDIF

END SUBROUTINE xfer_var_r2

!-------------------------------------------------------------------------

SUBROUTINE xfer_var_r3(typ, pos_nproma, pos_nblks, p_p, p_lp, arri, arro)

  INTEGER, INTENT(IN) :: typ, pos_nproma, pos_nblks
  TYPE(t_patch), INTENT(IN) :: p_p, p_lp
  REAL(wp), INTENT(IN)    :: arri(:,:,:)
  REAL(wp), INTENT(INOUT) :: arro(:,:,:)

  INTEGER :: j

  IF(pos_nproma==1 .and. pos_nblks==2) THEN
    DO j = 1, UBOUND(arri,3)
      CALL xfer_var_r2(typ,1,2,p_p,p_lp,arri(:,:,j),arro(:,:,j))
    ENDDO
  ELSEIF(pos_nproma==1 .and. pos_nblks==3) THEN
    DO j = 1, UBOUND(arri,2)
      CALL xfer_var_r2(typ,1,2,p_p,p_lp,arri(:,j,:),arro(:,j,:))
    ENDDO
  ELSEIF(pos_nproma==2 .and. pos_nblks==3) THEN
    DO j = 1, UBOUND(arri,1)
      CALL xfer_var_r2(typ,1,2,p_p,p_lp,arri(j,:,:),arro(j,:,:))
    ENDDO
  ELSE
    CALL finish ('mo_interpolation:xfer_var','Illegal value for pos_nproma/pos_nblks')
  ENDIF

END SUBROUTINE xfer_var_r3

!-------------------------------------------------------------------------

SUBROUTINE xfer_var_r4(typ, pos_nproma, pos_nblks, p_p, p_lp, arri, arro)

  INTEGER, INTENT(IN) :: typ, pos_nproma, pos_nblks
  TYPE(t_patch), INTENT(IN) :: p_p, p_lp
  REAL(wp), INTENT(IN)    :: arri(:,:,:,:)
  REAL(wp), INTENT(INOUT) :: arro(:,:,:,:)

  INTEGER :: i,j

  ! Variable has 4 dimensions
  IF(pos_nproma == 1 .AND. pos_nblks == 2)  THEN
    DO j = 1, UBOUND(arri,4)
    DO i = 1, UBOUND(arri,3)
      CALL xfer_var_r2(typ,1,2,p_p,p_lp,arri(:,:,i,j),arro(:,:,i,j))
    ENDDO
    ENDDO
  ELSEIF(pos_nproma == 1 .AND. pos_nblks == 3)  THEN
    DO j = 1, UBOUND(arri,4)
    DO i = 1, UBOUND(arri,2)
      CALL xfer_var_r2(typ,1,2,p_p,p_lp,arri(:,i,:,j),arro(:,i,:,j))
    ENDDO
    ENDDO
  ELSEIF(pos_nproma == 1 .AND. pos_nblks == 4)  THEN
    DO j = 1, UBOUND(arri,3)
    DO i = 1, UBOUND(arri,2)
      CALL xfer_var_r2(typ,1,2,p_p,p_lp,arri(:,i,j,:),arro(:,i,j,:))
    ENDDO
    ENDDO
  ELSEIF(pos_nproma == 3 .AND. pos_nblks == 4)  THEN
    DO j = 1, UBOUND(arri,2)
    DO i = 1, UBOUND(arri,1)
      CALL xfer_var_r2(typ,1,2,p_p,p_lp,arri(i,j,:,:),arro(i,j,:,:))
    ENDDO
    ENDDO
  ELSE
    ! Other pos_nproma/pos_nblks combinations are possible but currently not existing!
    CALL finish ('mo_interpolation:xfer_var','unsupported value for pos_nproma/pos_nblks')
  ENDIF


END SUBROUTINE xfer_var_r4
!-------------------------------------------------------------------------

SUBROUTINE xfer_var_i2(typ, pos_nproma, pos_nblks, p_p, p_lp, arri, arro)

  INTEGER, INTENT(IN) :: typ, pos_nproma, pos_nblks
  TYPE(t_patch), INTENT(IN) :: p_p, p_lp
  INTEGER,INTENT(IN) :: arri(:,:)
  INTEGER,INTENT(INOUT) :: arro(:,:)

  REAL(wp) :: r_arri(UBOUND(arri,1),UBOUND(arri,2))
  REAL(wp) :: r_arro(UBOUND(arro,1),UBOUND(arro,2))

  r_arri(:,:) = REAL(arri,wp)
  r_arro(:,:) = 0._wp ! Safety only
  CALL xfer_var_r2(typ,pos_nproma,pos_nblks,p_p,p_lp,r_arri,r_arro)
  arro(:,:) = INT(r_arro)

END SUBROUTINE xfer_var_i2

!-------------------------------------------------------------------------
!
!> xfer_idx family: transfer index variables from parent to local parent
!!
!! @par Revision History
!! Developed  by  Rainer Johanni (2011-10-26)
!-------------------------------------------------------------------------

SUBROUTINE xfer_idx_2(type_arr, type_idx, pos_nproma, pos_nblks, p_p, p_lp, idxi, blki, idxo, blko)

  INTEGER, INTENT(IN) :: type_arr, type_idx, pos_nproma, pos_nblks
  TYPE(t_patch), INTENT(IN), TARGET :: p_p, p_lp
  INTEGER, INTENT(IN)    :: idxi(:,:), blki(:,:)
  INTEGER, INTENT(INOUT) :: idxo(:,:), blko(:,:)

  INTEGER :: jb, jl, i_l, i_g, n_idx_l, n_idx_g
  REAL(wp) :: z_idxi(UBOUND(idxi,1),UBOUND(idxi,2))
  REAL(wp) :: z_idxo(UBOUND(idxo,1),UBOUND(idxo,2))
  INTEGER, POINTER :: glb_index(:)
  TYPE(t_grid_domain_decomp_info), POINTER :: local_decomp_info

  IF(type_idx == SYNC_C) THEN
    glb_index => p_p%cells%decomp_info%glb_index
    local_decomp_info => p_lp%cells%decomp_info
    n_idx_l = p_p%n_patch_cells
    n_idx_g = p_p%n_patch_cells_g
  ELSEIF(type_idx == SYNC_E) THEN
    glb_index => p_p%edges%decomp_info%glb_index
    local_decomp_info => p_lp%edges%decomp_info
    n_idx_l = p_p%n_patch_edges
    n_idx_g = p_p%n_patch_edges_g
  ELSEIF(type_idx == SYNC_V) THEN
    glb_index => p_p%verts%decomp_info%glb_index
    local_decomp_info => p_lp%verts%decomp_info
    n_idx_l = p_p%n_patch_verts
    n_idx_g = p_p%n_patch_verts_g
  ELSE
    CALL finish('xfer_idx','Unsupported type_idx')
  ENDIF

  z_idxi(:,:) = 0._wp
  z_idxo(:,:) = 0._wp

  DO jb = 1, UBOUND(idxi,2)
    DO jl = 1, nproma

      i_l = idx_1d(idxi(jl,jb),blki(jl,jb))

      IF(i_l <= 0 .or. i_l > n_idx_l) THEN
        z_idxi(jl,jb) = 0._wp
      ELSE
        z_idxi(jl,jb) = glb_index(i_l)
      ENDIF

    END DO
  END DO

  CALL xfer_var_r2(type_arr,pos_nproma,pos_nblks,p_p,p_lp,z_idxi,z_idxo)

  DO jb = 1, UBOUND(idxo,2)
    DO jl = 1, nproma

      i_g = INT(z_idxo(jl,jb))

      IF(i_g <= 0 .or. i_g > n_idx_g) THEN
        idxo(jl,jb) = 0
        blko(jl,jb) = 0
      ELSE
        i_l = get_valid_local_index(local_decomp_info%glb2loc_index, i_g)
        idxo(jl,jb) = idx_no(i_l)
        blko(jl,jb) = blk_no(i_l)
      ENDIF

    END DO
  END DO

END SUBROUTINE xfer_idx_2

!-------------------------------------------------------------------------

SUBROUTINE xfer_idx_3(type_arr, type_idx, pos_nproma, pos_nblks, p_p, p_lp, idxi, blki, idxo, blko)

  INTEGER, INTENT(IN) :: type_arr, type_idx, pos_nproma, pos_nblks
  TYPE(t_patch), INTENT(IN) :: p_p, p_lp
  INTEGER, INTENT(IN)    :: idxi(:,:,:), blki(:,:,:)
  INTEGER, INTENT(INOUT) :: idxo(:,:,:), blko(:,:,:)

  INTEGER :: j

  IF(pos_nproma==1 .and. pos_nblks==2) THEN
    DO j = 1, UBOUND(idxi,3)
      CALL xfer_idx_2(type_arr,type_idx,1,2,p_p,p_lp,idxi(:,:,j),blki(:,:,j), &
                                                   & idxo(:,:,j),blko(:,:,j))
    ENDDO
  ELSEIF(pos_nproma==1 .and. pos_nblks==3) THEN
    DO j = 1, UBOUND(idxi,2)
      CALL xfer_idx_2(type_arr,type_idx,1,2,p_p,p_lp,idxi(:,j,:),blki(:,j,:), &
                                                   & idxo(:,j,:),blko(:,j,:))
    ENDDO
  ELSEIF(pos_nproma==2 .and. pos_nblks==3) THEN
    DO j = 1, UBOUND(idxi,1)
      CALL xfer_idx_2(type_arr,type_idx,1,2,p_p,p_lp,idxi(j,:,:),blki(j,:,:), &
                                                   & idxo(j,:,:),blko(j,:,:))
    ENDDO
  ELSE
    CALL finish ('mo_interpolation:xfer_idx','Illegal value for pos_nproma/pos_nblks')
  ENDIF

END SUBROUTINE xfer_idx_3

!-------------------------------------------------------------------------
!!
!>
!!  Transfers interpolation state from parent to local parent
!!
!! @par Revision History
!! Developed  by  Rainer Johanni (2011-10-26)
!-------------------------------------------------------------------------
!
SUBROUTINE transfer_interpol_state(p_p, p_lp, pi, po)
!
  TYPE(t_patch), INTENT(IN)    :: p_p   ! parent
  TYPE(t_patch), INTENT(INOUT) :: p_lp  ! local parent

  TYPE(t_int_state), INTENT(IN)    :: pi ! Interpolation state on parent
  TYPE(t_int_state), INTENT(INOUT) :: po ! Interpolation state on local parent

  INTEGER, ALLOCATABLE :: owner(:)
  INTEGER :: j

  ! Allocate interpolation state for local parent

  CALL allocate_int_state(p_lp, po)

  ! Set up communication patterns for transferring the data to local parents.
  ! Since these communication patterns are not used elsewhere, they are
  ! stored locally and deleted at the end of the routine

  ALLOCATE(owner(MAX(p_lp%n_patch_cells, p_lp%n_patch_verts, &
    &                p_lp%n_patch_edges)))

  owner(1:p_lp%n_patch_cells) = &
    dist_dir_get_owners(p_p%cells%decomp_info%owner_dist_dir, &
      &                 p_lp%cells%decomp_info%glb_index(1:p_lp%n_patch_cells))
  CALL setup_comm_pattern(p_lp%n_patch_cells, owner(1:p_lp%n_patch_cells), &
    &                     p_lp%cells%decomp_info%glb_index,  &
    &                     p_p%cells%decomp_info%glb2loc_index, &
    &                     p_p%n_patch_cells, &
    &                     p_p%cells%decomp_info%owner_local, &
    &                     p_p%cells%decomp_info%glb_index, &
    &                     comm_pat_glb_to_loc_c)

  owner(1:p_lp%n_patch_edges) = &
    dist_dir_get_owners(p_p%edges%decomp_info%owner_dist_dir, &
      &                 p_lp%edges%decomp_info%glb_index(1:p_lp%n_patch_edges))
  CALL setup_comm_pattern(p_lp%n_patch_edges, owner(1:p_lp%n_patch_edges), &
    &                     p_lp%edges%decomp_info%glb_index,  &
    &                     p_p%edges%decomp_info%glb2loc_index, &
    &                     p_p%n_patch_edges, &
    &                     p_p%edges%decomp_info%owner_local, &
    &                     p_p%edges%decomp_info%glb_index, &
    &                     comm_pat_glb_to_loc_e)

  owner(1:p_lp%n_patch_verts) = &
    dist_dir_get_owners(p_p%verts%decomp_info%owner_dist_dir, &
      &                 p_lp%verts%decomp_info%glb_index(1:p_lp%n_patch_verts))
  CALL setup_comm_pattern(p_lp%n_patch_verts, owner(1:p_lp%n_patch_verts), &
    &                     p_lp%verts%decomp_info%glb_index,  &
    &                     p_p%verts%decomp_info%glb2loc_index, &
    &                     p_p%n_patch_verts, &
    &                     p_p%verts%decomp_info%owner_local, &
    &                     p_p%verts%decomp_info%glb_index, &
    &                     comm_pat_glb_to_loc_v)

  DEALLOCATE(owner)

  ! Some edge related values of the patch are only set in
  ! construct_2d_interpol_state (complete_patchinfo) and have to
  ! be set here also

  CALL xfer_var(SYNC_E,1,2,p_p,p_lp,p_p%edges%area_edge,p_lp%edges%area_edge)
  CALL xfer_var(SYNC_E,1,2,p_p,p_lp,p_p %edges%inv_primal_edge_length, &
                                  & p_lp%edges%inv_primal_edge_length)
  CALL xfer_var(SYNC_E,1,2,p_p,p_lp,p_p %edges%inv_dual_edge_length, &
                                  & p_lp%edges%inv_dual_edge_length)
  IF (p_p%geometry_info%cell_type == 3) THEN
  CALL xfer_var(SYNC_E,1,2,p_p,p_lp,p_p %edges%inv_vert_vert_length, &
                                  & p_lp%edges%inv_vert_vert_length)
  ENDIF
  CALL xfer_var(SYNC_E,1,2,p_p,p_lp,p_p %edges%primal_normal_cell(:,:,:)%v1, &
                                  & p_lp%edges%primal_normal_cell(:,:,:)%v1)
  CALL xfer_var(SYNC_E,1,2,p_p,p_lp,p_p %edges%primal_normal_cell(:,:,:)%v2, &
                                  & p_lp%edges%primal_normal_cell(:,:,:)%v2)
  CALL xfer_var(SYNC_E,1,2,p_p,p_lp,p_p %edges%primal_normal_vert(:,:,:)%v1, &
                                  & p_lp%edges%primal_normal_vert(:,:,:)%v1)
  CALL xfer_var(SYNC_E,1,2,p_p,p_lp,p_p %edges%primal_normal_vert(:,:,:)%v2, &
                                  & p_lp%edges%primal_normal_vert(:,:,:)%v2)
  CALL xfer_var(SYNC_E,1,2,p_p,p_lp,p_p %edges%dual_normal_cell(:,:,:)%v1, &
                                  & p_lp%edges%dual_normal_cell(:,:,:)%v1)
  CALL xfer_var(SYNC_E,1,2,p_p,p_lp,p_p %edges%dual_normal_cell(:,:,:)%v2, &
                                  & p_lp%edges%dual_normal_cell(:,:,:)%v2)
  CALL xfer_var(SYNC_E,1,2,p_p,p_lp,p_p %edges%dual_normal_vert(:,:,:)%v1, &
                                  & p_lp%edges%dual_normal_vert(:,:,:)%v1)
  CALL xfer_var(SYNC_E,1,2,p_p,p_lp,p_p %edges%dual_normal_vert(:,:,:)%v2, &
                                  & p_lp%edges%dual_normal_vert(:,:,:)%v2)

  CALL xfer_idx(SYNC_E,SYNC_V,1,2,p_p,p_lp,p_p %edges%vertex_idx,p_p %edges%vertex_blk, &
                                         & p_lp%edges%vertex_idx,p_lp%edges%vertex_blk)

  ! The same for quad_area, quad_orientation, quad_idx which is calculated
  ! after the complete patch has been read

  CALL xfer_var(SYNC_E,1,2,p_p,p_lp,p_p %edges%quad_area, &
                                  & p_lp%edges%quad_area)
  CALL xfer_var(SYNC_E,1,2,p_p,p_lp,p_p %edges%quad_orientation, &
                                  & p_lp%edges%quad_orientation)

  CALL xfer_idx(SYNC_E,SYNC_E,1,2,p_p,p_lp,p_p %edges%quad_idx,p_p %edges%quad_blk, &
                                         & p_lp%edges%quad_idx,p_lp%edges%quad_blk)


  ! Transfer interpolation state

  CALL xfer_var(SYNC_E,1,3,p_p,p_lp,pi%c_lin_e,po%c_lin_e)
  IF (p_p%geometry_info%cell_type == 3) THEN
  CALL xfer_var(SYNC_C,1,3,p_p,p_lp,pi%e_bln_c_s,po%e_bln_c_s)
  CALL xfer_var(SYNC_C,1,3,p_p,p_lp,pi%e_bln_c_u,po%e_bln_c_u)
  CALL xfer_var(SYNC_C,1,3,p_p,p_lp,pi%e_bln_c_v,po%e_bln_c_v)
  CALL xfer_var(SYNC_C,1,3,p_p,p_lp,pi%c_bln_avg,po%c_bln_avg)
  CALL xfer_var(SYNC_C,1,4,p_p,p_lp,pi%gradc_bmat,po%gradc_bmat)
  CALL xfer_var(SYNC_E,1,3,p_p,p_lp,pi%e_flx_avg,po%e_flx_avg)
  CALL xfer_var(SYNC_E,1,3,p_p,p_lp,pi%v_1o2_e,po%v_1o2_e)
  ENDIF
  CALL xfer_var(SYNC_C,1,3,p_p,p_lp,pi%e_inn_c,po%e_inn_c)
!-------------------------------------------------------------------------------
! Please note: for cell_type == 6 there exists no grid refinement
!              and thus no local parents!
!-------------------------------------------------------------------------------
!  IF (p_p%geometry_info%cell_type == 6) THEN
!  CALL xfer_var(SYNC_V,1,3,p_p,p_lp,pi%e_inn_v,po%e_inn_v)
!  CALL xfer_var(SYNC_C,1,3,p_p,p_lp,pi%e_aw_c,po%e_aw_c)
!  CALL xfer_var(SYNC_C,1,3,p_p,p_lp,pi%r_aw_c,po%r_aw_c)
!  CALL xfer_var(SYNC_V,1,3,p_p,p_lp,pi%e_aw_v,po%e_aw_v)
!  CALL xfer_var(SYNC_V,1,3,p_p,p_lp,pi%e_1o3_v,po%e_1o3_v)
!  CALL xfer_var(SYNC_E,1,3,p_p,p_lp,pi%tria_aw_rhom,po%tria_aw_rhom)
!  ENDIF
  CALL xfer_var(SYNC_C,1,3,p_p,p_lp,pi%verts_aw_cells,po%verts_aw_cells)
  CALL xfer_var(SYNC_V,1,3,p_p,p_lp,pi%cells_aw_verts,po%cells_aw_verts)
  CALL xfer_var(SYNC_V,1,3,p_p,p_lp,pi%cells_plwa_verts,po%cells_plwa_verts)
!  IF( p_p%geometry_info%cell_type == 6 ) THEN
!  CALL xfer_var(SYNC_V,2,3,p_p,p_lp,pi%tria_north,po%tria_north)
!  CALL xfer_var(SYNC_V,2,3,p_p,p_lp,pi%tria_east,po%tria_east)
!  CALL xfer_var(SYNC_C,1,3,p_p,p_lp,pi%hex_north,po%hex_north)
!  CALL xfer_var(SYNC_C,1,3,p_p,p_lp,pi%hex_east,po%hex_east)
!  IF (i_cori_method>=3) THEN
!  CALL xfer_var(SYNC_E,2,3,p_p,p_lp,pi%quad_north,po%quad_north)
!  CALL xfer_var(SYNC_E,2,3,p_p,p_lp,pi%quad_east,po%quad_east)
!  ENDIF
!  CALL xfer_var(SYNC_E,1,3,p_p,p_lp,pi%cno_en,po%cno_en)
!  CALL xfer_var(SYNC_E,1,3,p_p,p_lp,pi%cea_en,po%cea_en)
!  CALL xfer_var(SYNC_E,2,3,p_p,p_lp,pi%heli_coeff,po%heli_coeff)
!  IF (i_cori_method < 3) THEN
!  CALL xfer_idx(SYNC_E,SYNC_E,2,3,p_p,p_lp,pi%heli_vn_idx,pi%heli_vn_blk, &
!                                         & po%heli_vn_idx,po%heli_vn_blk)
!  ENDIF
!  ENDIF
  IF (p_p%geometry_info%cell_type == 3) THEN
  CALL xfer_idx(SYNC_C,SYNC_E,2,3,p_p,p_lp,pi%rbf_vec_idx_c,pi%rbf_vec_blk_c, &
                                         & po%rbf_vec_idx_c,po%rbf_vec_blk_c)
  CALL xfer_var(SYNC_C,1,2,p_p,p_lp,pi%rbf_vec_stencil_c,po%rbf_vec_stencil_c)
  CALL xfer_var(SYNC_C,3,4,p_p,p_lp,pi%rbf_vec_coeff_c,po%rbf_vec_coeff_c)
  CALL xfer_idx(SYNC_C,SYNC_C,2,3,p_p,p_lp,pi%rbf_c2grad_idx,pi%rbf_c2grad_blk, &
                                         & po%rbf_c2grad_idx,po%rbf_c2grad_blk)
  CALL xfer_var(SYNC_C,3,4,p_p,p_lp,pi%rbf_c2grad_coeff,po%rbf_c2grad_coeff)
  CALL xfer_idx(SYNC_V,SYNC_E,2,3,p_p,p_lp,pi%rbf_vec_idx_v,pi%rbf_vec_blk_v, &
                                         & po%rbf_vec_idx_v,po%rbf_vec_blk_v)
  CALL xfer_var(SYNC_V,1,2,p_p,p_lp,pi%rbf_vec_stencil_v,po%rbf_vec_stencil_v)
  CALL xfer_var(SYNC_V,3,4,p_p,p_lp,pi%rbf_vec_coeff_v,po%rbf_vec_coeff_v)
  CALL xfer_idx(SYNC_E,SYNC_E,2,3,p_p,p_lp,pi%rbf_vec_idx_e,pi%rbf_vec_blk_e, &
                                         & po%rbf_vec_idx_e,po%rbf_vec_blk_e)
  CALL xfer_var(SYNC_E,1,2,p_p,p_lp,pi%rbf_vec_stencil_e,po%rbf_vec_stencil_e)
  CALL xfer_var(SYNC_E,2,3,p_p,p_lp,pi%rbf_vec_coeff_e,po%rbf_vec_coeff_e)
  ENDIF

  IF (p_p%geometry_info%cell_type == 3) THEN
  CALL xfer_var(SYNC_E,1,3,p_p,p_lp,pi%geofac_qdiv,po%geofac_qdiv)
  CALL xfer_var(SYNC_E,1,3,p_p,p_lp,pi%geofac_grdiv,po%geofac_grdiv)
  CALL xfer_var(SYNC_C,1,2,p_p,p_lp,pi%nudgecoeff_c,po%nudgecoeff_c)
  CALL xfer_var(SYNC_E,1,2,p_p,p_lp,pi%nudgecoeff_e,po%nudgecoeff_e)

  CALL xfer_var(SYNC_C,1,2,p_p,p_lp,pi%gquad%qpts_tri_l(:,:)%lon,po%gquad%qpts_tri_l(:,:)%lon)
  CALL xfer_var(SYNC_C,1,2,p_p,p_lp,pi%gquad%qpts_tri_l(:,:)%lat,po%gquad%qpts_tri_l(:,:)%lat)
  CALL xfer_var(SYNC_C,1,2,p_p,p_lp,pi%gquad%qpts_tri_q(:,:,:)%lon,po%gquad%qpts_tri_q(:,:,:)%lon)
  CALL xfer_var(SYNC_C,1,2,p_p,p_lp,pi%gquad%qpts_tri_q(:,:,:)%lat,po%gquad%qpts_tri_q(:,:,:)%lat)
  CALL xfer_var(SYNC_C,1,2,p_p,p_lp,pi%gquad%qpts_tri_c(:,:,:)%lon,po%gquad%qpts_tri_c(:,:,:)%lon)
  CALL xfer_var(SYNC_C,1,2,p_p,p_lp,pi%gquad%qpts_tri_c(:,:,:)%lat,po%gquad%qpts_tri_c(:,:,:)%lat)

  po%gquad%weights_tri_q(:) = pi%gquad%weights_tri_q(:)
  po%gquad%weights_tri_c(:) = pi%gquad%weights_tri_c(:)

  ENDIF
  CALL xfer_var(SYNC_C,1,3,p_p,p_lp,pi%geofac_div,po%geofac_div)
  CALL xfer_var(SYNC_V,1,3,p_p,p_lp,pi%geofac_rot,po%geofac_rot)
  CALL xfer_var(SYNC_C,1,3,p_p,p_lp,pi%geofac_n2s,po%geofac_n2s)
  CALL xfer_var(SYNC_C,1,3,p_p,p_lp,pi%geofac_grg,po%geofac_grg)
  CALL xfer_var(SYNC_C,1,2,p_p,p_lp,pi%primal_normal_ec,po%primal_normal_ec)
  CALL xfer_var(SYNC_C,1,2,p_p,p_lp,pi%edge_cell_length,po%edge_cell_length)
  CALL xfer_var(SYNC_C,1,4,p_p,p_lp,pi%cell_vert_dist,po%cell_vert_dist)
!  IF (p_p%geometry_info%cell_type == 6) THEN
!  CALL xfer_idx(SYNC_E,SYNC_E,2,3,p_p,p_lp,pi%dir_gradh_i1,po%dir_gradh_i1, &
!                                         & po%dir_gradh_i1,po%dir_gradh_i1)
!  CALL xfer_idx(SYNC_E,SYNC_E,2,3,p_p,p_lp,pi%dir_gradh_i2,po%dir_gradh_i2, &
!                                         & po%dir_gradh_i2,po%dir_gradh_i2)
!  CALL xfer_var(SYNC_E,2,3,p_p,p_lp,pi%dir_gradhux_c1,po%dir_gradhux_c1)
!  CALL xfer_var(SYNC_E,2,3,p_p,p_lp,pi%dir_gradhux_c2,po%dir_gradhux_c2)
!  CALL xfer_var(SYNC_E,2,3,p_p,p_lp,pi%strain_def_c1,po%strain_def_c1)
!  CALL xfer_var(SYNC_E,2,3,p_p,p_lp,pi%strain_def_c2,po%strain_def_c2)
!  CALL xfer_idx(SYNC_E,SYNC_E,2,3,p_p,p_lp,pi%dir_gradt_i1,po%dir_gradt_i1, &
!                                         & po%dir_gradt_i1,po%dir_gradt_i1)
!  CALL xfer_idx(SYNC_E,SYNC_E,2,3,p_p,p_lp,pi%dir_gradt_i2,po%dir_gradt_i2, &
!                                         & po%dir_gradt_i2,po%dir_gradt_i2)
!  CALL xfer_var(SYNC_E,2,3,p_p,p_lp,pi%dir_gradtxy_v1,po%dir_gradtxy_v1)
!  CALL xfer_var(SYNC_E,2,3,p_p,p_lp,pi%dir_gradtxy_v2,po%dir_gradtxy_v2)
!  CALL xfer_var(SYNC_E,2,3,p_p,p_lp,pi%dir_gradtyx_v1,po%dir_gradtyx_v1)
!  CALL xfer_var(SYNC_E,2,3,p_p,p_lp,pi%dir_gradtyx_v2,po%dir_gradtyx_v2)
!  CALL xfer_var(SYNC_E,2,3,p_p,p_lp,pi%shear_def_v1,po%shear_def_v1)
!  CALL xfer_var(SYNC_E,2,3,p_p,p_lp,pi%shear_def_v2,po%shear_def_v2)
!  ENDIF

  ! clean up

!CDIR NOIEXPAND
  CALL delete_comm_pattern(comm_pat_glb_to_loc_c)
!CDIR NOIEXPAND
  CALL delete_comm_pattern(comm_pat_glb_to_loc_e)
!CDIR NOIEXPAND
  CALL delete_comm_pattern(comm_pat_glb_to_loc_v)

END SUBROUTINE transfer_interpol_state
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!
!
!>
!! Deallocation of components of a single 2d interpolation state.
!!
!!
!! @par Revision History
!! Split off from destruct_2d_interpol_state, Rainer Johanni (2010-10-26)
!!
SUBROUTINE deallocate_int_state( ptr_int)
!
TYPE(t_int_state), INTENT(inout) :: ptr_int

INTEGER :: ist

!-----------------------------------------------------------------------

  ! deallocate interpolation state
  !
  ! c_lin_e
  !
  DEALLOCATE (ptr_int%c_lin_e, STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:destruct_int_state',                      &
      &             'deallocation for c_lin_e failed')
  ENDIF
  !
  !
  ! e_bln_c_s
  !
  DEALLOCATE (ptr_int%e_bln_c_s, STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:destruct_int_state',                      &
      &             'deallocation for e_bln_c_s failed')
  ENDIF
  !
  ! e_bln_c_u
  !
  DEALLOCATE (ptr_int%e_bln_c_u, STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:destruct_int_state',                      &
      &             'deallocation for e_bln_c_u failed')
  ENDIF
  !
  ! e_bln_c_v
  !
  DEALLOCATE (ptr_int%e_bln_c_v, STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:destruct_int_state',                      &
      &             'deallocation for e_bln_c_v failed')
  ENDIF
  !
  ! c_bln_avg
  !
  DEALLOCATE (ptr_int%c_bln_avg, STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:destruct_int_state',                      &
      &             'deallocation for c_bln_avg failed')
  ENDIF
  !
  ! gradc_bmat
  !
  DEALLOCATE (ptr_int%gradc_bmat, STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:destruct_int_state',                      &
      &             'deallocation for gradc_bmat failed')
  ENDIF
  !
  ! e_flx_avg
  !
  DEALLOCATE (ptr_int%e_flx_avg, STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:destruct_int_state',                      &
      &             'deallocation for e_flx_avg failed')
  ENDIF
  !
  ! v_1o2_e
  !
  DEALLOCATE (ptr_int%v_1o2_e, STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:destruct_int_state',                      &
      &             'deallocation for v_1o2_e failed')
  ENDIF
  !
  ! e_inn_c
  !
  DEALLOCATE (ptr_int%e_inn_c, STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:destruct_int_state',                      &
      &             'deallocation for e_inn_c failed')
  ENDIF
  !
  !
  ! verts_aw_cells
  !
  DEALLOCATE (ptr_int%verts_aw_cells, STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:destruct_int_state',                      &
      &             'deallocation for verts_aw_cells failed')
  ENDIF
  !
  ! cells_aw_verts
  !
  DEALLOCATE (ptr_int%cells_aw_verts, STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:destruct_int_state',                      &
      &             'deallocation for cells_aw_verts failed')
  ENDIF
  !
  ! cells_plwa_verts
  !
  DEALLOCATE (ptr_int%cells_plwa_verts, STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:destruct_int_state',                      &
      &             'deallocation for cells_plwa_verts failed')
  ENDIF

    !
    ! rbf_vec_idx_c, rbf_vec_blk_c
    !
    DEALLOCATE (ptr_int%rbf_vec_idx_c, STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:destruct_int_state',                      &
        &             'deallocation for rbf_vec_idx_c failed')
    ENDIF
    DEALLOCATE (ptr_int%rbf_vec_blk_c, STAT=ist )
    IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:destruct_int_state',                      &
      &             'deallocation for rbf_vec_blk_c failed')
    ENDIF
    !
    ! rbf_vec_stencil_c
    !
    DEALLOCATE (ptr_int%rbf_vec_stencil_c, STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:destruct_int_state',                      &
        &             'deallocation for rbf_vec_stencil_c failed')
    ENDIF
    !
    ! rbf_vec_coeff_c
    !
    DEALLOCATE (ptr_int%rbf_vec_coeff_c, STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:destruct_int_state',                      &
        &             'deallocation for rbf_vec_coeff_c failed')
    ENDIF
    !
    ! rbf_c2grad_idx, rbf_c2grad_blk
    !
    DEALLOCATE (ptr_int%rbf_c2grad_idx, STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:destruct_int_state',                      &
        &             'deallocation for rbf_c2grad_idx failed')
    ENDIF
    DEALLOCATE (ptr_int%rbf_c2grad_blk, STAT=ist )
    IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:destruct_int_state',                      &
      &             'deallocation for rbf_c2grad_blk failed')
    ENDIF
    !
    ! rbf_c2grad_coeff_c
    !
    DEALLOCATE (ptr_int%rbf_c2grad_coeff, STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:destruct_int_state',                      &
        &             'deallocation for rbf_c2grad_coeff failed')
    ENDIF
    !
    ! rbf_vec_idx_v, rbf_vec_blk_v
    !
    DEALLOCATE (ptr_int%rbf_vec_idx_v, STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:destruct_int_state',                      &
        &             'deallocation for rbf_vec_idx_v failed')
    ENDIF
    DEALLOCATE (ptr_int%rbf_vec_blk_v, STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:destruct_int_state',                      &
        &             'deallocation for rbf_vec_idx_v failed')
    ENDIF
    !
    ! rbf_vec_stencil_v
    !
    DEALLOCATE (ptr_int%rbf_vec_stencil_v, STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:destruct_int_state',                      &
        &             'deallocation for rbf_vec_stencil_v failed')
    ENDIF
    !
    ! rbf_vec_coeff_v
    !
    DEALLOCATE (ptr_int%rbf_vec_coeff_v, STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:destruct_int_state',                      &
        &             'deallocation for rbf_vec_coeff_v failed')
    ENDIF
    !
    ! rbf_vec_idx_e, rbf_vec_blk_e
    !
    DEALLOCATE (ptr_int%rbf_vec_idx_e, STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:destruct_int_state',                      &
        &             'deallocation for rbf_vec_idx_e failed')
    ENDIF
    DEALLOCATE (ptr_int%rbf_vec_blk_e, STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:destruct_int_state',                      &
        &             'deallocation for rbf_vec_blk_e failed')
    ENDIF
    !
    ! rbf_vec_stencil_e
    !
    DEALLOCATE (ptr_int%rbf_vec_stencil_e, STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:destruct_int_state',                      &
        &             'deallocation for rbf_vec_stencil_e failed')
    ENDIF
    !
    ! rbf_vec_coeff_e
    !
    DEALLOCATE (ptr_int%rbf_vec_coeff_e, STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:destruct_int_state',                      &
      &             'deallocation for rbf_vec_coeff_e failed')
    ENDIF


  IF( ltransport .OR. iequations == 3) THEN
    !
    ! pos_on_tplane_e
    !
    DEALLOCATE (ptr_int%pos_on_tplane_e, STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:destruct_int_state',                      &
        &             'deallocation for pos_on_tplane_e failed')
    ENDIF
    !
    ! tplane_e_dotprod
    !
    DEALLOCATE (ptr_int%tplane_e_dotprod, STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:destruct_int_state',                      &
        &             'deallocation for tplane_e_dotprod failed')
    ENDIF

    !
    ! pos_on_tplane_c_edge
    !
    DEALLOCATE (ptr_int%pos_on_tplane_c_edge, STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:destruct_int_state',                      &
        &             'deallocation for pos_on_tplane_c_edge failed')
    ENDIF

    !
    ! Least squares reconstruction
    !

    !
    ! *** linear ***
    !
    ! lsq_dim_stencil
    !
    DEALLOCATE (ptr_int%lsq_lin%lsq_dim_stencil, STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',                   &
        &          'deallocation for lsq_dim_stencil failed')
    ENDIF
    !
    ! lsq_idx_c
    !
    DEALLOCATE (ptr_int%lsq_lin%lsq_idx_c, STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',                   &
        &          'deallocation for lsq_idx_c failed')
    ENDIF
    !
    ! lsq_blk_c
    !
    DEALLOCATE (ptr_int%lsq_lin%lsq_blk_c, STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',                   &
        &          'deallocation for lsq_blk_c failed')
    ENDIF
    !
    ! lsq_weights_c
    !
    DEALLOCATE (ptr_int%lsq_lin%lsq_weights_c, STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:destruct_int_state',                    &
        &          'deallocation for lsq_weights_c failed')
    ENDIF
    !
    ! lsq_qtmat_c
    !
    DEALLOCATE (ptr_int%lsq_lin%lsq_qtmat_c, STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:destruct_int_state',                    &
        &          'deallocation for lsq_qtmat_c failed')
    ENDIF
    !
    ! lsq_rmat_rdiag_c
    !
    DEALLOCATE (ptr_int%lsq_lin%lsq_rmat_rdiag_c, STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:destruct_int_state',                    &
        &          'deallocation for lsq_rmat_rdiag_c failed')
    ENDIF
    !
    ! lsq_rmat_utri_c
    !
    DEALLOCATE (ptr_int%lsq_lin%lsq_rmat_utri_c, STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:destruct_int_state',                    &
        &          'deallocation for lsq_rmat_utri_c failed')
    ENDIF
    !
    ! lsq_pseudoinv
    !
    DEALLOCATE (ptr_int%lsq_lin%lsq_pseudoinv, STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:destruct_int_state',                    &
        &          'deallocation for lsq_pseudoinv failed')
    ENDIF
    !
    ! lsq_moments
    !
    DEALLOCATE (ptr_int%lsq_lin%lsq_moments, STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:destruct_int_state',                    &
        &          'deallocation for lsq_moments failed')
    ENDIF
    !
    ! lsq_moments_hat
    !
    DEALLOCATE (ptr_int%lsq_lin%lsq_moments_hat, STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:destruct_int_state',                    &
        &          'deallocation for lsq_moments_hat failed')
    ENDIF

    !
    ! *** higher order ***
    !
    !
    ! lsq_dim_stencil
    !
    DEALLOCATE (ptr_int%lsq_high%lsq_dim_stencil, STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',                   &
        &          'deallocation for lsq_dim_stencil failed')
    ENDIF
    !
    ! lsq_idx_c
    !
    DEALLOCATE (ptr_int%lsq_high%lsq_idx_c, STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',                   &
        &          'deallocation for lsq_idx_c failed')
    ENDIF
    !
    ! lsq_blk_c
    !
    DEALLOCATE (ptr_int%lsq_high%lsq_blk_c, STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',                   &
        &          'deallocation for lsq_blk_c failed')
    ENDIF
    !
    ! lsq_weights_c
    !
    DEALLOCATE (ptr_int%lsq_high%lsq_weights_c, STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:destruct_int_state',                    &
        &          'deallocation for lsq_weights_c failed')
    ENDIF
    !
    ! lsq_qtmat_c
    !
    DEALLOCATE (ptr_int%lsq_high%lsq_qtmat_c, STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:destruct_int_state',                    &
        &          'deallocation for lsq_qtmat_c failed')
    ENDIF
    !
    ! lsq_rmat_rdiag_c
    !
    DEALLOCATE (ptr_int%lsq_high%lsq_rmat_rdiag_c, STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:destruct_int_state',                    &
        &          'deallocation for lsq_rmat_rdiag_c failed')
    ENDIF
    !
    ! lsq_rmat_utri_c
    !
    DEALLOCATE (ptr_int%lsq_high%lsq_rmat_utri_c, STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:destruct_int_state',                    &
        &          'deallocation for lsq_rmat_utri_c failed')
    ENDIF
    !
    ! lsq_pseudoinv
    !
    DEALLOCATE (ptr_int%lsq_high%lsq_pseudoinv, STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:destruct_int_state',                    &
        &          'deallocation for lsq_pseudoinv failed')
    ENDIF
    !
    ! lsq_moments
    !
    DEALLOCATE (ptr_int%lsq_high%lsq_moments, STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:destruct_int_state',                    &
        &          'deallocation for lsq_moments failed')
    ENDIF
    !
    ! lsq_moments_hat
    !
    DEALLOCATE (ptr_int%lsq_high%lsq_moments_hat, STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:destruct_int_state',                    &
        &          'deallocation for lsq_moments_hat failed')
    ENDIF
  END IF

  DEALLOCATE (ptr_int%geofac_qdiv, STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:destruct_int_state',                      &
      &             'deallocation for geofac_qdiv failed')
  ENDIF
  DEALLOCATE (ptr_int%geofac_grdiv, STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:destruct_int_state',                      &
      &             'deallocation for geofac_grdiv failed')
  ENDIF
  DEALLOCATE (ptr_int%gquad%qpts_tri_l, STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:destruct_int_state',                      &
      &             'deallocation for qpts_tri_l failed')
  ENDIF
  DEALLOCATE (ptr_int%gquad%qpts_tri_q, STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:destruct_int_state',                      &
      &             'deallocation for qpts_tri_q failed')
  ENDIF
  DEALLOCATE (ptr_int%gquad%qpts_tri_c, STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:destruct_int_state',                      &
      &             'deallocation for qpts_tri_q failed')
  ENDIF
  DEALLOCATE (ptr_int%gquad%weights_tri_q, STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:destruct_int_state',                      &
      &             'deallocation for weights_tri_q failed')
  ENDIF
  DEALLOCATE (ptr_int%gquad%weights_tri_c, STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:destruct_int_state',                      &
      &             'deallocation for weights_tri_c failed')
  ENDIF

  DEALLOCATE (ptr_int%geofac_div, STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:destruct_int_state',                      &
      &             'deallocation for geofac_div failed')
  ENDIF

  DEALLOCATE (ptr_int%geofac_rot, STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:destruct_int_state',                      &
      &             'deallocation for geofac_rot failed')
  ENDIF

  DEALLOCATE (ptr_int%geofac_n2s, STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:destruct_int_state',                      &
      &             'deallocation for geofac_n2s failed')
  ENDIF

  DEALLOCATE (ptr_int%geofac_grg, STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:destruct_int_state',                      &
      &             'deallocation for geofac_grg failed')
  ENDIF

  DEALLOCATE (ptr_int%primal_normal_ec, STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:destruct_int_state',                       &
      &             'deallocation for primal_normal_ec failed')
  ENDIF

  DEALLOCATE (ptr_int%edge_cell_length, STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:destruct_int_state',                       &
      &             'deallocation for edge_cell_length failed')
  ENDIF

  DEALLOCATE (ptr_int%cell_vert_dist, STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:destruct_int_state',                       &
      &             'deallocation for cell_vert_dist failed')
  ENDIF

END SUBROUTINE deallocate_int_state


!-------------------------------------------------------------------------
!
!
!>
!!               Deallocation of components of 2d interpolation state.
!!
!!
!! @par Revision History
!! Developed  by  Peter Korn, MPI-M (2005).
!! Modified by Luca Bonaventura, MPI-M (2005).
!!
SUBROUTINE destruct_2d_interpol_state( ptr_int_state)
  !
  TYPE(t_int_state), INTENT(inout) :: ptr_int_state(n_dom_start:)
  ! local variables:
  CHARACTER(*), PARAMETER :: routine = TRIM("mo_interpolation:destruct_int_state")
  INTEGER                 :: jg

  !-----------------------------------------------------------------------

  CALL message(routine, 'start to destruct int state')

  DO jg = n_dom_start, n_dom
    CALL deallocate_int_state(ptr_int_state(jg))
  ENDDO

  CALL message (routine, 'destruction of interpolation state finished')

END SUBROUTINE destruct_2d_interpol_state


END MODULE mo_intp_state
