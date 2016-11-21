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
!!  Modification by William Sawyer, CSCS (2015-01-27)
!!  - OpenACC implementation
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

MODULE mo_intp_rbf
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
USE mo_impl_constants,      ONLY: min_rlcell_int, min_rledge_int, min_rlvert_int
USE mo_model_domain,        ONLY: t_patch
USE mo_loopindices,         ONLY: get_indices_c, get_indices_e, get_indices_v
USE mo_intp_data_strc,      ONLY: t_int_state
USE mo_fortran_tools,       ONLY: init
#ifdef _OPENACC
USE mo_mpi,                 ONLY: i_am_accel_node
#endif

IMPLICIT NONE


PRIVATE

PUBLIC :: rbf_vec_interpol_cell, rbf_interpol_c2grad,     &
        & rbf_vec_interpol_vertex, rbf_vec_interpol_edge

#if defined( _OPENACC )
#define ACC_DEBUG NOACC
#if defined(__INTP_RBF_NOACC)
  LOGICAL, PARAMETER ::  acc_on = .FALSE.
#else
  LOGICAL, PARAMETER ::  acc_on = .TRUE.
#endif
#endif


CONTAINS

!-------------------------------------------------------------------------
!
!-------------------------------------------------------------------------
!
!
!>
!! Performs vector RBF reconstruction at cell center.
!!
!! Theory described in Narcowich and Ward (Math Comp. 1994) and
!! Bonaventura and Baudisch (Mox Report n. 75).
!! It takes edge based variables as input and combines them
!! into three dimensional cartesian vectors at each cell center.
!!
!! @par Revision History
!! Developed  by L.Bonaventura  (2004).
!! Enhanced efficiency with  direct solve for coefficients of data,
!! introduced by Will Sawyer (2004).
!! Adapted to new data structure by L.Bonaventura and P.Korn (2006).
!! This a new version of the subroutine with the cartesian
!!   components of the reconstructed
!!   vector as output (Pilar Ripodas November 2006)
!! Modifications by Tobias Ruppert, DWD (2007-02-08):
!! - replaced rbf_vec_dim by rbf_mat_dim
!! Modifications by Almut Gassmann, MPI-M (2007-04-30)
!! - abandon grid for the sake of patch
!! Modification by Pilar Ripodas, DWD (2007-07):
!! - substruct the outgoing component of the reconstructed
!!   vector
!! Modification by Guenther Zaengl, DWD (2009-02-13)
!! - change to direct reconstruction of vector components
!!
SUBROUTINE rbf_vec_interpol_cell( p_vn_in, ptr_patch, ptr_int, p_u_out,  &
  &                               p_v_out, opt_slev, opt_elev, opt_rlstart, opt_rlend )

! !INPUT PARAMETERS
!
!  patch on which computation is performed
!
TYPE(t_patch), TARGET, INTENT(in) ::  &
  &  ptr_patch

! input normal components of (velocity) vectors at edge midpoints
REAL(wp), INTENT(IN) ::  &
  &  p_vn_in(:,:,:) ! dim: (nproma,nlev,nblks_e)

! Indices of source points and interpolation coefficients
TYPE(t_int_state), TARGET, INTENT(IN)  ::  ptr_int

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_slev    ! optional vertical start level

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_elev    ! optional vertical end level

! start and end values of refin_ctrl flag
INTEGER, INTENT(in), OPTIONAL :: opt_rlstart, opt_rlend

! reconstructed x-component (u) of velocity vector
REAL(wp),INTENT(INOUT) ::  &
  &  p_u_out(:,:,:) ! dim: (nproma,nlev,nblks_c)

! reconstructed y-component (v) of velocity vector
REAL(wp),INTENT(INOUT) ::  &
  &  p_v_out(:,:,:) ! dim: (nproma,nlev,nblks_c)

! !LOCAL VARIABLES
INTEGER :: slev, elev                ! vertical start and end level
INTEGER :: jc, jk, jb                ! integer over cells, levels, and blocks,

INTEGER :: i_startblk      ! start block
INTEGER :: i_endblk        ! end block
INTEGER :: i_startidx      ! start index
INTEGER :: i_endidx        ! end index
INTEGER :: rl_start, rl_end, i_nchdom


INTEGER,  DIMENSION(:,:,:),   POINTER :: iidx, iblk
REAL(wp), DIMENSION(:,:,:,:), POINTER :: ptr_coeff

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
  elev = UBOUND(p_vn_in,2)
END IF

IF ( PRESENT(opt_rlstart) ) THEN
  rl_start = opt_rlstart
ELSE
  rl_start = 2
END IF
IF ( PRESENT(opt_rlend) ) THEN
  rl_end = opt_rlend
ELSE
  rl_end = min_rlcell_int-1
END IF

iidx => ptr_int%rbf_vec_idx_c
iblk => ptr_int%rbf_vec_blk_c

ptr_coeff => ptr_int%rbf_vec_coeff_c

! values for the blocking
i_nchdom   = MAX(1,ptr_patch%n_childdom)
i_startblk = ptr_patch%cells%start_blk(rl_start,1)
i_endblk   = ptr_patch%cells%end_blk(rl_end,i_nchdom)

#ifdef _OPENACC
!$ACC DATA PCOPYIN( p_vn_in ), PCOPY( p_u_out, p_v_out ), IF( i_am_accel_node .AND. acc_on )
!ACC_DEBUG UPDATE DEVICE( p_vn_in, p_u_out, p_v_out ), IF( i_am_accel_node .AND. acc_on )
!$ACC PARALLEL &
!$ACC PRESENT( ptr_patch, ptr_int, p_vn_in, p_u_out, p_v_out ), &
!$ACC IF( i_am_accel_node .AND. acc_on )

!$ACC LOOP GANG
#else
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk,jc), ICON_OMP_RUNTIME_SCHEDULE
#endif

DO jb = i_startblk, i_endblk

  CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
                     i_startidx, i_endidx, rl_start, rl_end)

!$ACC LOOP VECTOR COLLAPSE(2)
#ifdef __LOOP_EXCHANGE
  DO jc = i_startidx, i_endidx
    DO jk = slev, elev
#else
!CDIR UNROLL=2
  DO jk = slev, elev
    DO jc = i_startidx, i_endidx
#endif

      p_u_out(jc,jk,jb) =  &
        ptr_coeff(1,1,jc,jb)*p_vn_in(iidx(1,jc,jb),jk,iblk(1,jc,jb)) + &
        ptr_coeff(2,1,jc,jb)*p_vn_in(iidx(2,jc,jb),jk,iblk(2,jc,jb)) + &
        ptr_coeff(3,1,jc,jb)*p_vn_in(iidx(3,jc,jb),jk,iblk(3,jc,jb)) + &
        ptr_coeff(4,1,jc,jb)*p_vn_in(iidx(4,jc,jb),jk,iblk(4,jc,jb)) + &
        ptr_coeff(5,1,jc,jb)*p_vn_in(iidx(5,jc,jb),jk,iblk(5,jc,jb)) + &
        ptr_coeff(6,1,jc,jb)*p_vn_in(iidx(6,jc,jb),jk,iblk(6,jc,jb)) + &
        ptr_coeff(7,1,jc,jb)*p_vn_in(iidx(7,jc,jb),jk,iblk(7,jc,jb)) + &
        ptr_coeff(8,1,jc,jb)*p_vn_in(iidx(8,jc,jb),jk,iblk(8,jc,jb)) + &
        ptr_coeff(9,1,jc,jb)*p_vn_in(iidx(9,jc,jb),jk,iblk(9,jc,jb))
      p_v_out(jc,jk,jb) =  &
        ptr_coeff(1,2,jc,jb)*p_vn_in(iidx(1,jc,jb),jk,iblk(1,jc,jb)) + &
        ptr_coeff(2,2,jc,jb)*p_vn_in(iidx(2,jc,jb),jk,iblk(2,jc,jb)) + &
        ptr_coeff(3,2,jc,jb)*p_vn_in(iidx(3,jc,jb),jk,iblk(3,jc,jb)) + &
        ptr_coeff(4,2,jc,jb)*p_vn_in(iidx(4,jc,jb),jk,iblk(4,jc,jb)) + &
        ptr_coeff(5,2,jc,jb)*p_vn_in(iidx(5,jc,jb),jk,iblk(5,jc,jb)) + &
        ptr_coeff(6,2,jc,jb)*p_vn_in(iidx(6,jc,jb),jk,iblk(6,jc,jb)) + &
        ptr_coeff(7,2,jc,jb)*p_vn_in(iidx(7,jc,jb),jk,iblk(7,jc,jb)) + &
        ptr_coeff(8,2,jc,jb)*p_vn_in(iidx(8,jc,jb),jk,iblk(8,jc,jb)) + &
        ptr_coeff(9,2,jc,jb)*p_vn_in(iidx(9,jc,jb),jk,iblk(9,jc,jb))

    ENDDO
  ENDDO

ENDDO
#ifdef _OPENACC
!$ACC END PARALLEL
!ACC_DEBUG UPDATE HOST(p_u_out,p_v_out), IF( i_am_accel_node .AND. acc_on )
!$ACC END DATA
#else
!$OMP END DO NOWAIT
!$OMP END PARALLEL
#endif

END SUBROUTINE rbf_vec_interpol_cell
!====================================================================================


!====================================================================================
SUBROUTINE rbf_interpol_c2grad( p_cell_in, ptr_patch, ptr_int, grad_x,  &
  &                             grad_y, opt_slev, opt_elev, opt_rlstart, opt_rlend )

! !INPUT PARAMETERS
!
!  patch on which computation is performed
!
TYPE(t_patch), TARGET, INTENT(in) ::  &
  &  ptr_patch

! input cell-based variable for which gradient at cell center is computed
REAL(wp), INTENT(IN) ::  &
  &  p_cell_in(:,:,:) ! dim: (nproma,nlev,nblks_c)

! Indices of source points and interpolation coefficients
TYPE(t_int_state), TARGET, INTENT(IN)  ::  ptr_int

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_slev    ! optional vertical start level

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_elev    ! optional vertical end level

! start and end values of refin_ctrl flag
INTEGER, INTENT(in), OPTIONAL :: opt_rlstart, opt_rlend

! reconstructed zonal (x) component of gradient vector
REAL(wp),INTENT(INOUT) ::  &
  &  grad_x(:,:,:) ! dim: (nproma,nlev,nblks_c)

! reconstructed zonal (x) component of gradient vector
REAL(wp),INTENT(INOUT) ::  &
  &  grad_y(:,:,:) ! dim: (nproma,nlev,nblks_c)

! !LOCAL VARIABLES
INTEGER :: slev, elev                ! vertical start and end level
INTEGER :: jc, jk, jb                ! integer over cells, levels, and blocks,

INTEGER :: i_startblk      ! start block
INTEGER :: i_endblk        ! end block
INTEGER :: i_startidx      ! start index
INTEGER :: i_endidx        ! end index
INTEGER :: rl_start, rl_end, i_nchdom


INTEGER,  DIMENSION(:,:,:),   POINTER :: iidx, iblk
REAL(wp), DIMENSION(:,:,:,:), POINTER :: ptr_coeff

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
  rl_end = min_rlcell_int
END IF


iidx => ptr_int%rbf_c2grad_idx
iblk => ptr_int%rbf_c2grad_blk

ptr_coeff => ptr_int%rbf_c2grad_coeff

! values for the blocking
i_nchdom   = MAX(1,ptr_patch%n_childdom)
i_startblk = ptr_patch%cells%start_blk(rl_start,1)
i_endblk   = ptr_patch%cells%end_blk(rl_end,i_nchdom)

#ifndef _OPENACC
!$OMP PARALLEL
#endif

IF (ptr_patch%id > 1) THEN
#ifdef _OPENACC
!$ACC KERNELS IF ( i_am_accel_node .AND. acc_on )
  grad_x(:,:,1:i_startblk) = 0._wp
  grad_y(:,:,1:i_startblk) = 0._wp
!$ACC END KERNELS
#else
  CALL init(grad_x(:,:,1:i_startblk))
  CALL init(grad_y(:,:,1:i_startblk))
!$OMP BARRIER
#endif
ENDIF

#ifdef _OPENACC
!$ACC DATA PCOPYIN( p_cell_in ), PCOPY( grad_x, grad_y ), IF( i_am_accel_node .AND. acc_on )
!ACC_DEBUG UPDATE DEVICE( p_cell_in, grad_x, grad_y ), IF( i_am_accel_node .AND. acc_on )
!$ACC PARALLEL &
!$ACC PRESENT( ptr_patch, ptr_int, p_cell_in, grad_x, grad_y ), &
!$ACC IF( i_am_accel_node .AND. acc_on )

!$ACC LOOP GANG
#else
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk,jc), ICON_OMP_RUNTIME_SCHEDULE
#endif

DO jb = i_startblk, i_endblk

  CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
                     i_startidx, i_endidx, rl_start, rl_end)

!$ACC LOOP VECTOR COLLAPSE(2)
#ifdef __LOOP_EXCHANGE
  DO jc = i_startidx, i_endidx
    DO jk = slev, elev
#else
!CDIR UNROLL=3
  DO jk = slev, elev
    DO jc = i_startidx, i_endidx
#endif

      grad_x(jc,jk,jb) =  &
        ptr_coeff(1,1,jc,jb)*p_cell_in(jc,jk,jb) + &
        ptr_coeff(2,1,jc,jb)*p_cell_in(iidx(2,jc,jb),jk,iblk(2,jc,jb)) + &
        ptr_coeff(3,1,jc,jb)*p_cell_in(iidx(3,jc,jb),jk,iblk(3,jc,jb)) + &
        ptr_coeff(4,1,jc,jb)*p_cell_in(iidx(4,jc,jb),jk,iblk(4,jc,jb)) + &
        ptr_coeff(5,1,jc,jb)*p_cell_in(iidx(5,jc,jb),jk,iblk(5,jc,jb)) + &
        ptr_coeff(6,1,jc,jb)*p_cell_in(iidx(6,jc,jb),jk,iblk(6,jc,jb)) + &
        ptr_coeff(7,1,jc,jb)*p_cell_in(iidx(7,jc,jb),jk,iblk(7,jc,jb)) + &
        ptr_coeff(8,1,jc,jb)*p_cell_in(iidx(8,jc,jb),jk,iblk(8,jc,jb)) + &
        ptr_coeff(9,1,jc,jb)*p_cell_in(iidx(9,jc,jb),jk,iblk(9,jc,jb)) + &
        ptr_coeff(10,1,jc,jb)*p_cell_in(iidx(10,jc,jb),jk,iblk(10,jc,jb))
      grad_y(jc,jk,jb) =  &
        ptr_coeff(1,2,jc,jb)*p_cell_in(jc,jk,jb) + &
        ptr_coeff(2,2,jc,jb)*p_cell_in(iidx(2,jc,jb),jk,iblk(2,jc,jb)) + &
        ptr_coeff(3,2,jc,jb)*p_cell_in(iidx(3,jc,jb),jk,iblk(3,jc,jb)) + &
        ptr_coeff(4,2,jc,jb)*p_cell_in(iidx(4,jc,jb),jk,iblk(4,jc,jb)) + &
        ptr_coeff(5,2,jc,jb)*p_cell_in(iidx(5,jc,jb),jk,iblk(5,jc,jb)) + &
        ptr_coeff(6,2,jc,jb)*p_cell_in(iidx(6,jc,jb),jk,iblk(6,jc,jb)) + &
        ptr_coeff(7,2,jc,jb)*p_cell_in(iidx(7,jc,jb),jk,iblk(7,jc,jb)) + &
        ptr_coeff(8,2,jc,jb)*p_cell_in(iidx(8,jc,jb),jk,iblk(8,jc,jb)) + &
        ptr_coeff(9,2,jc,jb)*p_cell_in(iidx(9,jc,jb),jk,iblk(9,jc,jb)) + &
        ptr_coeff(10,2,jc,jb)*p_cell_in(iidx(10,jc,jb),jk,iblk(10,jc,jb))

    ENDDO
  ENDDO

ENDDO
#ifdef _OPENACC
!$ACC END PARALLEL
!ACC_DEBUG UPDATE HOST( grad_x, grad_y ), IF( i_am_accel_node .AND. acc_on )
!$ACC END DATA
#else
!$OMP END DO NOWAIT
!$OMP END PARALLEL
#endif

END SUBROUTINE rbf_interpol_c2grad

!-------------------------------------------------------------------------
!
!
!>
!! Performs vector RBF reconstruction at triangle vertices.
!!
!! Theory described in Narcowich and Ward (Math Comp. 1994) and
!! Bonaventura and Baudisch (Mox Report n. 75).
!! It takes edge based variables as input and combines them
!! into three dimensional cartesian vectors at each vertex.
!!
!! @par Revision History
!! Developed  by Jochen Foerstner, DWD (2008-05-02)
!! Modification by Guenther Zaengl, DWD (2009-02-13)
!! - change to direct reconstruction of vector components
!! Modification by Almut Gassmann, MPI-M (2009-11-06)
!! - include switch which distinguishes normal or tangential components as input
!!
SUBROUTINE rbf_vec_interpol_vertex( p_e_in, ptr_patch, ptr_int, &
                                    p_u_out, p_v_out,           &
                                    opt_slev, opt_elev, opt_rlstart, opt_rlend )
!
TYPE(t_patch), TARGET, INTENT(in) ::  &
  &  ptr_patch

! input components of velocity or horizontal vorticity vectors at edge midpoints
REAL(wp), INTENT(IN) ::  &
  &  p_e_in(:,:,:) ! dim: (nproma,nlev,nblks_e)

! Indices of source points and interpolation coefficients
TYPE(t_int_state), TARGET, INTENT(IN)  ::  ptr_int

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_slev    ! optional vertical start level

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_elev    ! optional vertical end level

! start and end values of refin_ctrl flag
INTEGER, INTENT(in), OPTIONAL :: opt_rlstart, opt_rlend

! reconstructed x-component (u) of velocity vector
REAL(wp),INTENT(INOUT) ::  &
  &  p_u_out(:,:,:) ! dim: (nproma,nlev,nblks_v)

! reconstructed y-component (v) of velocity vector
REAL(wp),INTENT(INOUT) ::  &
  &  p_v_out(:,:,:) ! dim: (nproma,nlev,nblks_v)

! !LOCAL VARIABLES

INTEGER :: slev, elev                ! vertical start and end level
INTEGER :: jv, jk, jb                ! integer over vertices, levels, and blocks,

INTEGER :: i_startblk      ! start block
INTEGER :: i_endblk        ! end block
INTEGER :: i_startidx      ! start index
INTEGER :: i_endidx        ! end index
INTEGER :: rl_start, rl_end, i_nchdom

INTEGER,  DIMENSION(:,:,:),   POINTER :: iidx, iblk
REAL(wp), DIMENSION(:,:,:,:), POINTER :: ptr_coeff

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
  elev = UBOUND(p_e_in,2)
END IF

IF ( PRESENT(opt_rlstart) ) THEN
  rl_start = opt_rlstart
ELSE
  rl_start = 2
END IF
IF ( PRESENT(opt_rlend) ) THEN
  rl_end = opt_rlend
ELSE
  rl_end = min_rlvert_int-1
END IF

iidx => ptr_int%rbf_vec_idx_v
iblk => ptr_int%rbf_vec_blk_v

ptr_coeff => ptr_int%rbf_vec_coeff_v

! values for the blocking
i_nchdom   = MAX(1,ptr_patch%n_childdom)
i_startblk = ptr_patch%verts%start_blk(rl_start,1)
i_endblk   = ptr_patch%verts%end_blk(rl_end,i_nchdom)

#ifdef _OPENACC
!$ACC DATA PCOPYIN( p_e_in ), PCOPY( p_u_out, p_v_out ), IF( i_am_accel_node .AND. acc_on )
!ACC_DEBUG UPDATE DEVICE( p_e_in, p_u_out, p_v_out ), IF( i_am_accel_node .AND. acc_on )
!$ACC PARALLEL &
!$ACC PRESENT( ptr_patch, ptr_int, p_e_in, p_u_out, p_v_out ), &
!$ACC IF( i_am_accel_node .AND. acc_on )

!$ACC LOOP GANG
#else
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk,jv), ICON_OMP_RUNTIME_SCHEDULE
#endif
DO jb = i_startblk, i_endblk

  CALL get_indices_v(ptr_patch, jb, i_startblk, i_endblk, &
                     i_startidx, i_endidx, rl_start, rl_end)

!$ACC LOOP VECTOR COLLAPSE(2)
#ifdef __LOOP_EXCHANGE
  DO jv = i_startidx, i_endidx
    DO jk = slev, elev
#else
!CDIR UNROLL=6
  DO jk = slev, elev
    DO jv = i_startidx, i_endidx
#endif

      p_u_out(jv,jk,jb) =  &
        ptr_coeff(1,1,jv,jb)*p_e_in(iidx(1,jv,jb),jk,iblk(1,jv,jb)) + &
        ptr_coeff(2,1,jv,jb)*p_e_in(iidx(2,jv,jb),jk,iblk(2,jv,jb)) + &
        ptr_coeff(3,1,jv,jb)*p_e_in(iidx(3,jv,jb),jk,iblk(3,jv,jb)) + &
        ptr_coeff(4,1,jv,jb)*p_e_in(iidx(4,jv,jb),jk,iblk(4,jv,jb)) + &
        ptr_coeff(5,1,jv,jb)*p_e_in(iidx(5,jv,jb),jk,iblk(5,jv,jb)) + &
        ptr_coeff(6,1,jv,jb)*p_e_in(iidx(6,jv,jb),jk,iblk(6,jv,jb))
      p_v_out(jv,jk,jb) =  &
        ptr_coeff(1,2,jv,jb)*p_e_in(iidx(1,jv,jb),jk,iblk(1,jv,jb)) + &
        ptr_coeff(2,2,jv,jb)*p_e_in(iidx(2,jv,jb),jk,iblk(2,jv,jb)) + &
        ptr_coeff(3,2,jv,jb)*p_e_in(iidx(3,jv,jb),jk,iblk(3,jv,jb)) + &
        ptr_coeff(4,2,jv,jb)*p_e_in(iidx(4,jv,jb),jk,iblk(4,jv,jb)) + &
        ptr_coeff(5,2,jv,jb)*p_e_in(iidx(5,jv,jb),jk,iblk(5,jv,jb)) + &
        ptr_coeff(6,2,jv,jb)*p_e_in(iidx(6,jv,jb),jk,iblk(6,jv,jb))

      ENDDO
    ENDDO

ENDDO

#ifdef _OPENACC
!$ACC END PARALLEL
!ACC_DEBUG UPDATE HOST( p_u_out, p_v_out ), IF( i_am_accel_node .AND. acc_on )
!$ACC END DATA
#else
!$OMP END DO NOWAIT
!$OMP END PARALLEL
#endif

END SUBROUTINE rbf_vec_interpol_vertex

!-------------------------------------------------------------------------
!
!
!>
!! Performs vector RBF reconstruction at edge midpoints.
!!
!! Theory described in Narcowich and Ward (Math Comp. 1994) and
!! Bonaventura and Baudisch (Mox Report n. 75).
!! It takes edge based variables as input and combines them
!! into three dimensional cartesian vectors at each edge.
!!
!! @par Revision History
!! Developed  by Jochen Foerstner, DWD (2008-07-15)
!! Modification by Guenther Zaengl, DWD (2009-02-13)
!! - change to direct reconstruction of tangential vector component
!!
SUBROUTINE rbf_vec_interpol_edge( p_vn_in, ptr_patch, ptr_int, p_vt_out,  &
  &                               opt_slev, opt_elev, opt_rlstart, opt_rlend )
!
TYPE(t_patch), TARGET, INTENT(in) ::  &
  &  ptr_patch

! input normal components of velocity vectors at edge midpoints
REAL(wp), INTENT(IN) ::  &
  &  p_vn_in(:,:,:) ! dim: (nproma,nlev,nblks_e)

! Indices of source points and interpolation coefficients
TYPE(t_int_state), TARGET, INTENT(IN)  ::  ptr_int

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_slev    ! optional vertical start level

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_elev    ! optional vertical end level

! start and end values of refin_ctrl flag
INTEGER, INTENT(in), OPTIONAL :: opt_rlstart, opt_rlend

! reconstructed tangential velocity component
REAL(wp),INTENT(INOUT) ::  &
  &  p_vt_out(:,:,:) ! dim: (nproma,nlev,nblks_e)


INTEGER :: slev, elev                ! vertical start and end level
INTEGER :: je, jk, jb                ! integer over edges, levels, and blocks,

INTEGER :: i_startblk      ! start block
INTEGER :: i_endblk        ! end block
INTEGER :: i_startidx      ! start index
INTEGER :: i_endidx        ! end index
INTEGER :: rl_start, rl_end, i_nchdom

INTEGER,  DIMENSION(:,:,:), POINTER :: iidx, iblk
REAL(wp), DIMENSION(:,:,:), POINTER :: ptr_coeff

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
  elev = UBOUND(p_vn_in,2)
END IF

IF ( PRESENT(opt_rlstart) ) THEN
  rl_start = opt_rlstart
ELSE
  rl_start = 2
END IF
IF ( PRESENT(opt_rlend) ) THEN
  rl_end = opt_rlend
ELSE
  rl_end = min_rledge_int-2
END IF

iidx => ptr_int%rbf_vec_idx_e
iblk => ptr_int%rbf_vec_blk_e

ptr_coeff => ptr_int%rbf_vec_coeff_e

! values for the blocking
i_nchdom   = MAX(1,ptr_patch%n_childdom)
i_startblk = ptr_patch%edges%start_blk(rl_start,1)
i_endblk   = ptr_patch%edges%end_blk(rl_end,i_nchdom)

#ifdef _OPENACC
!$ACC DATA PCOPYIN( p_vn_in ), PCOPY( p_vt_out ), IF( i_am_accel_node .AND. acc_on )
!ACC_DEBUG UPDATE DEVICE( p_vn_in, p_vt_out ), IF( i_am_accel_node .AND. acc_on )
!$ACC PARALLEL &
!$ACC PRESENT( ptr_patch, ptr_int, p_vn_in, p_vt_out ), &
!$ACC IF( i_am_accel_node .AND. acc_on )

!$ACC LOOP GANG
#else
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk,je) ICON_OMP_DEFAULT_SCHEDULE
#endif
  DO jb = i_startblk, i_endblk

    CALL get_indices_e(ptr_patch, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, rl_start, rl_end)

!$ACC LOOP VECTOR COLLAPSE(2)
#ifdef __LOOP_EXCHANGE
    DO je = i_startidx, i_endidx
      DO jk = slev, elev
#else
!CDIR UNROLL=3
    DO jk = slev, elev
      DO je = i_startidx, i_endidx
#endif

        p_vt_out(je,jk,jb) =  &
          ptr_coeff(1,je,jb)*p_vn_in(iidx(1,je,jb),jk,iblk(1,je,jb)) + &
          ptr_coeff(2,je,jb)*p_vn_in(iidx(2,je,jb),jk,iblk(2,je,jb)) + &
          ptr_coeff(3,je,jb)*p_vn_in(iidx(3,je,jb),jk,iblk(3,je,jb)) + &
          ptr_coeff(4,je,jb)*p_vn_in(iidx(4,je,jb),jk,iblk(4,je,jb))

      ENDDO
    ENDDO
  ENDDO
#ifdef _OPENACC
!$ACC END PARALLEL
!ACC_DEBUG UPDATE HOST( p_vt_out ), IF( i_am_accel_node .AND. acc_on )
!$ACC END DATA
#else
!$OMP END DO NOWAIT
!$OMP END PARALLEL
#endif

END SUBROUTINE rbf_vec_interpol_edge


END MODULE mo_intp_rbf
