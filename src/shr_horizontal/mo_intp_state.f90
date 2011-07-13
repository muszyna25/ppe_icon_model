
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
USE mo_io_units,            ONLY: nnml, nnml_output
USE mo_exception,           ONLY: message, finish
USE mo_impl_constants,      ONLY: SUCCESS
USE mo_model_domain,        ONLY: t_patch
USE mo_model_domain_import, ONLY: n_dom, n_dom_start, lplane, l_limited_area, lfeedback
USE mo_namelist,            ONLY: position_nml, POSITIONED
USE mo_impl_constants,      ONLY: MAX_CHAR_LENGTH
USE mo_parallel_configuration,  ONLY: nproma
USE mo_grid_configuration,  ONLY: global_cell_type
USE mo_run_config,          ONLY: ltransport
USE mo_dynamics_config,     ONLY: dynamics_config
USE mo_mpi,                 ONLY: p_pe, p_io

USE mo_interpol_nml
USE mo_intp_data_strc
USE mo_intp_rbf_coeffs
USE mo_intp_coeffs



IMPLICIT NONE

PRIVATE

!!CJPUBLIC :: setup_interpol, 
PUBLIC :: construct_2d_interpol_state, destruct_2d_interpol_state
PUBLIC :: allocate_int_state, deallocate_int_state

CHARACTER(len=*), PARAMETER, PRIVATE :: version = '$Id$'


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
INTEGER :: ist
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
  
  IF (ptr_patch%cell_type == 3) THEN
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
    ! e_flx_avg
    !
    ALLOCATE (ptr_int%e_flx_avg(nproma,5,nblks_e), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state', &
      &            'allocation for e_flx_avg failed')
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
  ALLOCATE (ptr_int%e_inn_c(nproma,ptr_patch%cell_type,nblks_c), STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:construct_int_state', &
    &            'allocation for e_inn_c failed')
  ENDIF
  !
  IF (ptr_patch%cell_type == 6) THEN
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
  ALLOCATE (ptr_int%verts_aw_cells(nproma,ptr_patch%cell_type,nblks_c), STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:construct_int_state',                     &
      &             'allocation for verts_aw_cells failed')
  ENDIF
  !
  ! cells_aw_verts
  !
  ALLOCATE (ptr_int%cells_aw_verts(nproma,9-ptr_patch%cell_type,nblks_v), STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:construct_int_state',                     &
      &             'allocation for cells_aw_verts failed')
  ENDIF
  !
  IF( ptr_patch%cell_type == 6 ) THEN
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

  IF (ptr_patch%cell_type == 3) THEN
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

  IF( ltransport .OR. &
      dynamics_config(ptr_patch%id)%iequations == 3) THEN
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

  IF (ptr_patch%cell_type == 3) THEN
    ALLOCATE (ptr_int%geofac_qdiv(nproma, 4, nblks_e), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',                     &
        &             'allocation for geofac_qdiv failed')
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

  ALLOCATE (ptr_int%geofac_div(nproma, ptr_patch%cell_type, nblks_c), STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:construct_int_state',                       &
      &             'allocation for geofac_div failed')
  ENDIF

  ALLOCATE (ptr_int%geofac_rot(nproma, 9-ptr_patch%cell_type, nblks_v), STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:construct_int_state', &
    &             'allocation for geofac_rot failed')
  ENDIF

  ALLOCATE (ptr_int%geofac_n2s(nproma, ptr_patch%cell_type+1, nblks_c), STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:construct_int_state',                       &
      &             'allocation for geofac_n2s failed')
  ENDIF

  ALLOCATE (ptr_int%geofac_grg(nproma, ptr_patch%cell_type+1, nblks_c, 2), STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:construct_int_state',                       &
      &             'allocation for geofac_grg failed')
  ENDIF

  ALLOCATE (ptr_int%cart_edge_coord(nproma, nblks_e, 3), STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:construct_int_state',                       &
      &             'allocation for cart_edge_coord failed')
  ENDIF

  ALLOCATE (ptr_int%cart_cell_coord(nproma, nblks_c, 3), STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:construct_int_state',                       &
      &             'allocation for cart_cell_coord failed')
  ENDIF

  ALLOCATE (ptr_int%primal_normal_ec(nproma, nblks_c,ptr_patch%cell_type, 2), STAT=ist)
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:construct_int_state',                       &
      &             'allocation for primal_normal_ec failed')
  ENDIF

  ALLOCATE (ptr_int%edge_cell_length(nproma, nblks_c, ptr_patch%cell_type), STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:construct_int_state',                       &
      &             'allocation for edge_cell_length failed')
  ENDIF


  IF (ptr_patch%cell_type == 6) THEN

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

  !
  ! initialize all components
  !
  IF (ptr_patch%cell_type == 3 ) THEN
    ptr_int%e_bln_c_s     = 0._wp
    ptr_int%e_bln_c_u     = 0._wp
    ptr_int%e_bln_c_v     = 0._wp
    ptr_int%c_bln_avg     = 0._wp
    ptr_int%e_flx_avg     = 0._wp
  ENDIF

  ptr_int%v_1o2_e       = 0.5_wp
  ptr_int%c_lin_e       = 0._wp
  ptr_int%e_inn_c       = 0._wp
  ptr_int%verts_aw_cells= 0._wp
  ptr_int%cells_aw_verts= 0._wp

  IF (ptr_patch%cell_type == 6 ) THEN
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

  IF( ptr_patch%cell_type == 3) THEN
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

  IF( ptr_patch%cell_type == 6 ) THEN
    ptr_int%heli_coeff        = 0._wp
    IF (i_cori_method < 3) THEN
      ptr_int%heli_vn_idx       = 0
      ptr_int%heli_vn_blk       = 0
    ENDIF
  ENDIF

  IF( ltransport .OR. &
      dynamics_config(ptr_patch%id)%iequations == 3) THEN

    ptr_int%pos_on_tplane_e   = 0._wp
    ptr_int%tplane_e_dotprod  = 0._wp

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

  IF (ptr_patch%cell_type ==3) THEN
    ptr_int%geofac_qdiv = 0._wp
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
  ptr_int%cart_edge_coord = 0._wp
  ptr_int%cart_cell_coord = 0._wp
  ptr_int%primal_normal_ec = 0._wp
  ptr_int%edge_cell_length = 0._wp

  IF(ptr_patch%cell_type==6) THEN
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

TYPE(t_int_state), INTENT(INOUT) :: ptr_int_state(n_dom_start:)

INTEGER :: jg

CHARACTER(len=MAX_CHAR_LENGTH) :: text

  !-----------------------------------------------------------------------

CALL message('mo_intp_state:construct_2d_interpol_state','start to construct int_state')

DO jg = n_dom_start, n_dom

  WRITE(text,'(a,i0)') 'constructing int_state for patch ',jg
  CALL message('mo_intp_state:construct_2d_interpol_state',text)

  CALL allocate_int_state(ptr_patch(jg), ptr_int_state(jg))
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
  IF (ptr_patch(jg)%cell_type==3)THEN
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
  IF (ptr_patch(jg)%cell_type == 3) THEN

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

    !
    ! initialization of quadrature points and weights
    !
    CALL tri_quadrature_pts (ptr_patch(jg), ptr_int_state(jg))
  ENDIF
  !
  ! initialization of coeffs for hexagon
  !
  IF (ptr_patch(jg)%cell_type==6) THEN
    CALL compute_heli_bra_coeff_idx(ptr_patch(jg), ptr_int_state(jg))
  ENDIF
  !
  ! - Initialization of tangential plane (at edge midpoints) for calculation
  !   of backward trajectories.
  ! - stencil generation
  ! - initialization of coefficients for least squares gradient
  ! reconstruction at cell centers
  !
  IF ( (ltransport .OR. dynamics_config(jg)%iequations == 3) .AND. &
        (.NOT. lplane)) THEN

    CALL init_tplane_e(ptr_patch(jg), ptr_int_state(jg))

    IF (ptr_patch(jg)%cell_type==3) THEN
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

ENDDO

CALL message('mo_intp_state:construct_2d_interpol_state', &
  & 'construction of interpolation state finished')

END SUBROUTINE construct_2d_interpol_state

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
SUBROUTINE deallocate_int_state( iequations, ptr_int )
!
INTEGER,INTENT(IN) :: iequations
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
  IF (global_cell_type == 3 ) THEN
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
    ! e_flx_avg
    !
    DEALLOCATE (ptr_int%e_flx_avg, STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:destruct_int_state',                      &
        &             'deallocation for e_flx_avg failed')
    ENDIF
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
  IF (global_cell_type == 6 ) THEN
    !
    ! e_inn_v
    !
    DEALLOCATE (ptr_int%e_inn_v, STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:destruct_int_state',                      &
        &             'deallocation for e_inn_v failed')
    ENDIF
    !
    ! e_aw_c
    !
    DEALLOCATE (ptr_int%e_aw_c, STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:destruct_int_state',                      &
        &             'deallocation for e_aw_c failed')
    ENDIF
    !
    ! r_aw_c
    !
    DEALLOCATE (ptr_int%r_aw_c, STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:destruct_int_state',                      &
        &             'deallocation for r_aw_c failed')
    ENDIF
    !
    ! e_aw_v
    !
    DEALLOCATE (ptr_int%e_aw_v, STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:destruct_int_state',                      &
        &             'deallocation for e_aw_v failed')
    ENDIF
    !
    ! tria_aw_rhom
    !
    DEALLOCATE (ptr_int%tria_aw_rhom, STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:destruct_int_state',                      &
        &             'deallocation for tria_aw_rhom failed')
    ENDIF
    !
    ! e_1o3_v
    !
    DEALLOCATE (ptr_int%e_1o3_v, STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:destruct_int_state',                    &
        &             'deallocation for e_1o3_v failed')
    ENDIF
    !
    IF (i_cori_method>=3) THEN
      !
      ! quad_north
      !
      DEALLOCATE (ptr_int%quad_north, STAT=ist )
      IF (ist /= SUCCESS) THEN
        CALL finish ('mo_interpolation:destruct_int_state',                   &
          &          'deallocation for quad_north failed')
      ENDIF
      !
      ! quad_east
      !
      DEALLOCATE (ptr_int%quad_east, STAT=ist )
      IF (ist /= SUCCESS) THEN
        CALL finish ('mo_interpolation:destruct_int_state',                   &
          &          'deallocation for quad_east failed')
      ENDIF
    ENDIF
    !
    ! hex_north
    !
    DEALLOCATE (ptr_int%hex_north, STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:destruct_int_state',                   &
        &          'deallocation for hex_north failed')
    ENDIF
    !
    ! hex_east
    !
    DEALLOCATE (ptr_int%hex_east, STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:destruct_int_state',                   &
        &          'deallocation for hex_east failed')
    ENDIF
    !
    ! tria_north
    !
    DEALLOCATE (ptr_int%tria_north, STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:destruct_int_state',                   &
        &          'deallocation for tria_north failed')
    ENDIF
    !
    ! tria_east
    !
    DEALLOCATE (ptr_int%tria_east, STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:destruct_int_state',                   &
        &          'deallocation for tria_east failed')
    ENDIF
    !
    ! cno_en
    !
    DEALLOCATE (ptr_int%cno_en, STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:destruct_int_state',                   &
        &          'deallocation for cno_en failed')
    ENDIF
    !
    ! cea_en
    !
    DEALLOCATE (ptr_int%cea_en, STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:destruct_int_state',                   &
        &          'deallocation for cea_en failed')
    ENDIF
  ENDIF
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

  IF (global_cell_type == 3) THEN
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

  IF(global_cell_type == 6) THEN
    DEALLOCATE (ptr_int%heli_coeff, STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:destruct_int_state',                   &
        &          'deallocation for heli_coeff failed')
    ENDIF
    IF (i_cori_method<3) THEN
      DEALLOCATE (ptr_int%heli_vn_idx, STAT=ist )
      IF (ist /= SUCCESS) THEN
        CALL finish ('mo_interpolation:destruct_int_state',                   &
          &          'deallocation for heli_vn_idx failed')
      ENDIF
      DEALLOCATE (ptr_int%heli_vn_blk, STAT=ist )
      IF (ist /= SUCCESS) THEN
        CALL finish ('mo_interpolation:destruct_int_state',                   &
          &          'deallocation for heli_vn_blk failed')
      ENDIF
    ENDIF
  ENDIF

  IF (global_cell_type == 3) THEN
    DEALLOCATE (ptr_int%geofac_qdiv, STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:destruct_int_state',                      &
        &             'deallocation for geofac_qdiv failed')
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

  DEALLOCATE (ptr_int%cart_edge_coord, STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:destruct_int_state',                       &
      &             'deallocation for cart_edge_coord failed')
  ENDIF

  DEALLOCATE (ptr_int%cart_cell_coord, STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:destruct_int_state',                       &
      &             'deallocation for cart_cell_coord failed')
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


  IF (global_cell_type == 6) THEN

    DEALLOCATE (ptr_int%dir_gradh_i1, STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',                    &
        &             'deallocation for dir_gradh_i1 failed')
    ENDIF
    DEALLOCATE (ptr_int%dir_gradh_i2, STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',                    &
        &             'deallocation for dir_gradh_i2 failed')
    ENDIF
    DEALLOCATE (ptr_int%dir_gradh_b1, STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',                    &
        &             'deallocation for dir_gradh_b1 failed')
    ENDIF
    DEALLOCATE (ptr_int%dir_gradh_b2, STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',                    &
        &             'deallocation for dir_gradh_b2 failed')
    ENDIF
    DEALLOCATE (ptr_int%dir_gradhux_c1, STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',                    &
        &             'deallocation for dir_gradhux_c1 failed')
    ENDIF
    DEALLOCATE (ptr_int%dir_gradhux_c2, STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',                    &
        &             'deallocation for dir_gradhux_c2 failed')
    ENDIF
    DEALLOCATE (ptr_int%strain_def_c1, STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',                    &
        &             'deallocation for strain_def_c1 failed')
    ENDIF
    DEALLOCATE (ptr_int%strain_def_c2, STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',                    &
        &             'deallocation for strain_def_c2 failed')
    ENDIF

    DEALLOCATE (ptr_int%dir_gradt_i1, STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',                    &
        &             'deallocation for dir_gradt_i1 failed')
    ENDIF
    DEALLOCATE (ptr_int%dir_gradt_i2, STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',                    &
        &             'deallocation for dir_gradt_i2 failed')
    ENDIF
    DEALLOCATE (ptr_int%dir_gradt_b1, STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',                    &
        &             'deallocation for dir_gradt_b1 failed')
    ENDIF
    DEALLOCATE (ptr_int%dir_gradt_b2, STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',                    &
        &             'deallocation for dir_gradt_b2 failed')
    ENDIF
    DEALLOCATE (ptr_int%dir_gradtxy_v1, STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',                    &
        &             'deallocation for dir_gradtxy_v1 failed')
    ENDIF
    DEALLOCATE (ptr_int%dir_gradtxy_v2, STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',                    &
        &             'deallocation for dir_gradtxy_v2 failed')
    ENDIF
    DEALLOCATE (ptr_int%dir_gradtyx_v1, STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',                    &
        &             'deallocation for dir_gradtyx_v1 failed')
    ENDIF
    DEALLOCATE (ptr_int%dir_gradtyx_v2, STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',                    &
        &             'deallocation for dir_gradtyx_v2 failed')
    ENDIF
    DEALLOCATE (ptr_int%shear_def_v1, STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',                    &
        &             'deallocation for shear_def_v1 failed')
    ENDIF
    DEALLOCATE (ptr_int%shear_def_v2, STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',                    &
        &             'deallocation for shear_def_v2 failed')
    ENDIF

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
SUBROUTINE destruct_2d_interpol_state( ptr_int_state )
!
TYPE(t_int_state), INTENT(inout) :: ptr_int_state(n_dom_start:)

INTEGER :: jg

!-----------------------------------------------------------------------

CALL message('mo_interpolation:destruct_int_state',                          &
  & 'start to destruct int state')

DO jg = n_dom_start, n_dom
  CALL deallocate_int_state(dynamics_config(jg)%iequations, ptr_int_state(jg))
ENDDO

CALL message ('mo_interpolation:destruct_int_state',                         &
  &              'destruction of interpolation state finished')

END SUBROUTINE destruct_2d_interpol_state


END MODULE mo_intp_state
