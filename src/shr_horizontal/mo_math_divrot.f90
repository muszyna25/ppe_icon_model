!>
!!   Contains the implementation of the div,rot,recon mathematical operators.
!!
!!   Contains the implementation of the mathematical operators
!!   employed by the shallow water prototype.
!!
!! @par Revision History
!!  Developed  by Luca Bonaventura and Will Sawyer (2002-4).
!!  Modified to ProTeX-style by  Luca Bonaventura and Thomas Heinze (2004).
!!  Adapted to new data structure by Thomas Heinze,
!!  Peter Korn and Luca Bonaventura (2005).
!!  Modification by Thomas Heinze (2006-02-21):
!!  - renamed m_modules to mo_modules
!!  Subroutine for divergence multiplied by area added by P.Korn (2006).
!!  Modification by Peter Korn, MPI-M, (2006-11-23):
!!  - replacements in TYPE patch: ic by l2g_c, ie by l2g_e, iv by l2g_v,
!!    iic by g2l_c, iie by g2l_e, iiv by g2l_v
!!  - replaced edge_index by edge_idx
!!  - replaced vertex_index by vertex_idx
!!  - replaced cell_index by cell_idx
!!  - replaced neighbor_index by neighbor_idx
!!  - replaced child_index by child_idx
!!  Modified by P Ripodas (2007-02):
!!  - include the system orientation factor in the vorticity term of nabla2_vec
!!  - solved errors in nabla4_vec and nabla4_scalar
!!  Modification by Peter Korn, MPI-M, (2006/2007):
!!  -operator overloading of curl operator and nabla2vec to handle atmosphere and ocean version
!!  -change of input/output arguments of subroutines: arrays of fixed size are
!!   changed to pointers, to avoid occurence of not-initialized numbers
!!  Modified by Almut Gassmann, MPI-M, (2007-04)
!!  - removed references to unused halo_verts
!!  - summing over all halos corresponding to different parallel patches
!!  Modified by Hui Wan, MPI-M, (2007-11)
!!  - added subroutine cell_avg
!!  Modification by Jochen Foerstner, DWD, (2008-05-05)
!!  - div and div_times_area are now generic subroutines
!!  - the divergence can now be computed either
!!    using the midpoint rule
!!    (div_midpoint, div_midpoint_times_area) or
!!    using the Simpson's rule
!!    (div_simpson, div_simpson_times_area)
!!  Modification by Jochen Foerstner, DWD, (2008-07-16)
!!  - introduction of several new operators (to be) used in combination with
!!    the tracer advection:
!!    grad_green_gauss_cell, grad_green_gauss_edge and div_quad_twoadjcells
!!    (for the new div operator there is again a version using the midpoint
!!    and a version using the Simpson's rule).
!!    The first operator is used to calculate a cell centered value of the
!!    gradient for the piecewise linear reconstruction. The second and third
!!    will be used in combination with the MPDATA scheme. Both deal with the
!!    quadrilateral control volumes formed by two adjacent triangles.
!!  Modification by Marco Restelli, MPI (2008-07-17)
!!  - included subroutine dtan.
!!  Modification by Jochen Foerstner, DWD (2008-09-12)
!!  - moved SUBROUTINE ravtom_normgrad2 from mo_interpolation to this module
!!    because of conflicting use statements.
!!  Modification by Jochen Foerstner, DWD (2008-09-16)
!!  - removed SUBROUTINE ravtom_normgrad2 (not used)
!!  Modification by Daniel Reinert, DWD (2009-07-20)
!!  - added subroutine grad_lsq_cell for gradient reconstruction via the
!!    least-squares method and grad_green_gauss_gc_cell for Green-Gauss
!!    gradient in geographical coordinates
!!  Modification by Daniel Reinert, DWD (2009-12-14)
!!  - renamed grad_lsq_cell -> recon_lsq_cell_l
!!  Modification by Leonidas Linardakis, MPI-M (2010-21-01)
!!  - splitted mo_math_operators into submodules
!!  Modification by Daniel Reinert, DWD (2010-04-12)
!!  - added subroutine recon_lsq_cell_q for third order accurate least-squares
!!    reconstruction of an arbitrary field. Based on cell centered values.
!!  Modification by Daniel Reinert, DWD (2010-10-14)
!!  - added subroutine recon_lsq_cell_c for fitting a cubic polynomial in a least
!!    squares sense. Based on cell centered values.
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
MODULE mo_math_divrot
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
USE mo_impl_constants,      ONLY: min_rlcell, min_rledge, min_rlvert
USE mo_interpolation,       ONLY: t_int_state, t_lsq, lsq_high_set
USE mo_model_domain,        ONLY: t_patch
USE mo_model_domain_import, ONLY: l_limited_area
USE mo_parallel_config,  ONLY: nproma
USE mo_exception,           ONLY: finish
USE mo_loopindices,         ONLY: get_indices_c, get_indices_e, get_indices_v
! USE mo_timer,              ONLY: timer_start, timer_stop, timer_div

IMPLICIT NONE

PRIVATE

CHARACTER(len=*), PARAMETER :: version = '$Id$'


PUBLIC :: recon_lsq_cell_l, recon_lsq_cell_l_svd
PUBLIC :: recon_lsq_cell_q, recon_lsq_cell_q_svd
PUBLIC :: recon_lsq_cell_cpoor, recon_lsq_cell_cpoor_svd
PUBLIC :: recon_lsq_cell_c, recon_lsq_cell_c_svd
PUBLIC :: div, div_avg
PUBLIC :: div_quad_twoadjcells
PUBLIC :: rot_vertex
PUBLIC :: rot_vertex_atmos


INTERFACE rot_vertex

  MODULE PROCEDURE rot_vertex_atmos

END INTERFACE


INTERFACE div

  MODULE PROCEDURE div3d
  MODULE PROCEDURE div4d

END INTERFACE


CONTAINS


!-------------------------------------------------------------------------
!
!
!>
!! Computes coefficients (i.e. derivatives) for cell centered linear
!! reconstruction.
!!
!! DESCRIPTION:
!! recon: reconstruction of subgrid distribution
!! lsq  : least-squares method
!! cell : solution coefficients defined at cell center
!! l    : linear reconstruction
!!
!! The least squares approach is used. Solves Rx = Q^T d.
!! R: upper triangular matrix (2 x 2)
!! Q: orthogonal matrix (3 x 2)
!! d: input vector (3 x 1)
!! x: solution vector (2 x 1)
!! works only on triangular grid yet
!!
!! @par Revision History
!! Developed and tested by Daniel Reinert, DWD (2009-07-20)
!! Modification by Daniel Reinert, DWD (2010-04-09)
!! - included explicit computation of the constant c0. c0 differs
!!   from the tracer value at the cell barycenter in the case
!!   of a conservative reconstruction.
!!
SUBROUTINE recon_lsq_cell_l( p_cc, ptr_patch, ptr_int_lsq, p_coeff, &
  &                           opt_slev, opt_elev, opt_rlstart,      &
  &                           opt_rlend )

  TYPE(t_patch), TARGET, INTENT(IN) :: &  !< patch on which computation 
    &  ptr_patch                          !<is performed

  TYPE(t_lsq), TARGET, INTENT(IN) :: &  !< data structure for interpolation
    &  ptr_int_lsq

  REAL(wp), INTENT(IN)          ::  &   !<  cell centered variable
    &  p_cc(:,:,:)

  INTEGER, INTENT(IN), OPTIONAL ::  &   !< optional vertical start level
    &  opt_slev

  INTEGER, INTENT(IN), OPTIONAL ::  &   !< optional vertical end level
    &  opt_elev

  INTEGER, INTENT(IN), OPTIONAL ::  &   !< start and end values of refin_ctrl flag
    &  opt_rlstart, opt_rlend

  REAL(wp), INTENT(INOUT) ::  &  !< cell based coefficients (geographical components)
    &  p_coeff(:,:,:,:)          !< (constant and gradients in latitudinal and
                                 !< longitudinal direction)

  REAL(wp)  ::   &               !< weights * difference of scalars i j
    &  z_d(nproma,ptr_patch%nlev,3)
  REAL(wp)  ::   &               !< matrix product of transposed Q matrix and d
    &  z_qt_times_d(2)

  INTEGER, POINTER ::   &             !< Pointer to line and block indices of
    &  iidx(:,:,:), iblk(:,:,:)       !< required stencil
  INTEGER :: slev, elev               !< vertical start and end level
  INTEGER :: jc, jk, jb               !< index of cell, vertical level and block
  INTEGER :: rl_start, rl_end
  INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx, i_nchdom


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
    elev = ptr_patch%nlev
  END IF
  IF ( PRESENT(opt_rlstart) ) THEN
    rl_start = opt_rlstart
  ELSE
    rl_start = 2
  END IF
  IF ( PRESENT(opt_rlend) ) THEN
    rl_end = opt_rlend
  ELSE
    rl_end = min_rlcell
  END IF


  ! values for the blocking
  i_nchdom   = MAX(1,ptr_patch%n_childdom)
  i_startblk = ptr_patch%cells%start_blk(rl_start,1)
  i_endblk   = ptr_patch%cells%end_blk(rl_end,i_nchdom)

  ! pointers to line and block indices of stencil
  iidx => ptr_patch%cells%neighbor_idx
  iblk => ptr_patch%cells%neighbor_blk



  !
  ! 1. reconstruction of cell based gradient (geographical components)
  !
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,jk,i_startidx,i_endidx,z_d,z_qt_times_d), SCHEDULE(runtime)
  DO jb = i_startblk, i_endblk

    CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, rl_start, rl_end)

#ifdef __LOOP_EXCHANGE
    DO jc = i_startidx, i_endidx
      DO jk = slev, elev
#else
!CDIR UNROLL=4
    DO jk = slev, elev
      DO jc = i_startidx, i_endidx
#endif

        ! note that the multiplication with lsq_weights_c(jc,js,jb) at
        ! runtime is now avoided. Instead, the multiplication with
        ! lsq_weights_c(jc,js,jb) has been shifted into the transposed
        ! Q-matrix.
!PK: here the landboundary is taken into account
        z_d(jc,jk,1) = p_cc(iidx(jc,jb,1),jk,iblk(jc,jb,1)) - p_cc(jc,jk,jb)
        z_d(jc,jk,2) = p_cc(iidx(jc,jb,2),jk,iblk(jc,jb,2)) - p_cc(jc,jk,jb)
        z_d(jc,jk,3) = p_cc(iidx(jc,jb,3),jk,iblk(jc,jb,3)) - p_cc(jc,jk,jb)

      END DO ! end loop over cells
    END DO ! end loop over vertical levels

    DO jk = slev, elev
      DO jc = i_startidx, i_endidx

        ! matrix multiplication Q^T d (partitioned into 2 dot products)
        z_qt_times_d(1) = ptr_int_lsq%lsq_qtmat_c(jc,1,1,jb) * z_d(jc,jk,1)  &
          &             + ptr_int_lsq%lsq_qtmat_c(jc,1,2,jb) * z_d(jc,jk,2)  &
          &             + ptr_int_lsq%lsq_qtmat_c(jc,1,3,jb) * z_d(jc,jk,3)
        z_qt_times_d(2) = ptr_int_lsq%lsq_qtmat_c(jc,2,1,jb) * z_d(jc,jk,1)  &
          &             + ptr_int_lsq%lsq_qtmat_c(jc,2,2,jb) * z_d(jc,jk,2)  &
          &             + ptr_int_lsq%lsq_qtmat_c(jc,2,3,jb) * z_d(jc,jk,3)


        ! Solve linear system by backward substitution
        ! Gradient in zonal and meridional direction
        !
        ! meridional
        p_coeff(jc,jk,jb,3) = ptr_int_lsq%lsq_rmat_rdiag_c(jc,2,jb) * z_qt_times_d(2)

        ! zonal
        p_coeff(jc,jk,jb,2) = ptr_int_lsq%lsq_rmat_rdiag_c(jc,1,jb)                  &
          & * (z_qt_times_d(1) - ptr_int_lsq%lsq_rmat_utri_c(jc,1,jb)                &
          & * p_coeff(jc,jk,jb,3))

        ! constant
        p_coeff(jc,jk,jb,1) = p_cc(jc,jk,jb)                                         &
          &                 - p_coeff(jc,jk,jb,2) * ptr_int_lsq%lsq_moments(jc,jb,1) &
          &                 - p_coeff(jc,jk,jb,3) * ptr_int_lsq%lsq_moments(jc,jb,2)


      END DO ! end loop over cells
    END DO ! end loop over vertical levels

  END DO ! end loop over blocks
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE recon_lsq_cell_l



!-------------------------------------------------------------------------
!
!
!>
!! Computes coefficients (i.e. derivatives) for cell centered linear
!! reconstruction.
!!
!! DESCRIPTION:
!! recon: reconstruction of subgrid distribution
!! lsq  : least-squares method
!! cell : solution coefficients defined at cell center
!! l    : linear reconstruction
!!
!! The least squares approach is used. Solves Ax = b via Singular 
!! Value Decomposition (SVD)
!! x = PINV(A) * b
!!
!! Matrices have the following size and shape:
!! PINV(A): Pseudo or Moore-Penrose inverse of A (via SVD) (2 x 3)
!! b: input vector (3 x 1)
!! x: solution vector (2 x 1)
!! only works on triangular grid yet
!!
!! @par Revision History
!! Developed and tested by Daniel Reinert, DWD (2011-05-26)
!!
SUBROUTINE recon_lsq_cell_l_svd( p_cc, ptr_patch, ptr_int_lsq, p_coeff, &
  &                           opt_slev, opt_elev, opt_rlstart,      &
  &                           opt_rlend )

  TYPE(t_patch), TARGET, INTENT(IN) :: &  !< patch on which computation 
    &  ptr_patch                          !< is performed

  TYPE(t_lsq), TARGET, INTENT(IN) :: &  !< data structure for interpolation
    &  ptr_int_lsq

  REAL(wp), INTENT(IN)          ::  &   !<  cell centered variable
    &  p_cc(:,:,:)

  INTEGER, INTENT(IN), OPTIONAL ::  &   !< optional vertical start level
    &  opt_slev

  INTEGER, INTENT(IN), OPTIONAL ::  &   !< optional vertical end level
    &  opt_elev

  INTEGER, INTENT(IN), OPTIONAL ::  &   !< start and end values of refin_ctrl flag
    &  opt_rlstart, opt_rlend

  REAL(wp), INTENT(INOUT) ::  &  !< cell based coefficients (geographical components)
    &  p_coeff(:,:,:,:)          !< (constant and gradients in latitudinal and
                                 !< longitudinal direction)

  REAL(wp)  ::   &               !< weights * difference of scalars i j
    &  z_b(nproma,ptr_patch%nlev,3)

  INTEGER, POINTER ::   &            !< Pointer to line and block indices of
    &  iidx(:,:,:), iblk(:,:,:)      !< required stencil
  INTEGER :: slev, elev              !< vertical start and end level
  INTEGER :: jc, jk, jb              !< index of cell, vertical level and block
  INTEGER :: rl_start, rl_end
  INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx, i_nchdom


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
    elev = ptr_patch%nlev
  END IF
  IF ( PRESENT(opt_rlstart) ) THEN
    rl_start = opt_rlstart
  ELSE
    rl_start = 2
  END IF
  IF ( PRESENT(opt_rlend) ) THEN
    rl_end = opt_rlend
  ELSE
    rl_end = min_rlcell
  END IF


  ! values for the blocking
  i_nchdom   = MAX(1,ptr_patch%n_childdom)
  i_startblk = ptr_patch%cells%start_blk(rl_start,1)
  i_endblk   = ptr_patch%cells%end_blk(rl_end,i_nchdom)

  ! pointers to line and block indices of stencil
  iidx => ptr_patch%cells%neighbor_idx
  iblk => ptr_patch%cells%neighbor_blk



  !
  ! 1. reconstruction of cell based gradient (geographical components)
  !
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,jk,i_startidx,i_endidx,z_b), SCHEDULE(runtime)
  DO jb = i_startblk, i_endblk

    CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, rl_start, rl_end)

#ifdef __LOOP_EXCHANGE
    DO jc = i_startidx, i_endidx
      DO jk = slev, elev
#else
!CDIR UNROLL=4
    DO jk = slev, elev
      DO jc = i_startidx, i_endidx
#endif

        ! note that the multiplication with lsq_weights_c(jc,js,jb) at
        ! runtime is now avoided. Instead, the weights have been shifted 
        ! into the pseudoinverse.
        z_b(jc,jk,1) = p_cc(iidx(jc,jb,1),jk,iblk(jc,jb,1)) - p_cc(jc,jk,jb)
        z_b(jc,jk,2) = p_cc(iidx(jc,jb,2),jk,iblk(jc,jb,2)) - p_cc(jc,jk,jb)
        z_b(jc,jk,3) = p_cc(iidx(jc,jb,3),jk,iblk(jc,jb,3)) - p_cc(jc,jk,jb)

      END DO ! end loop over cells
    END DO ! end loop over vertical levels

    !
    ! 2. compute cell based coefficients for linear reconstruction
    !    calculate matrix vector product PINV(A) * b
    !
    DO jk = slev, elev
      DO jc = i_startidx, i_endidx

        ! meridional
        p_coeff(jc,jk,jb,3) = ptr_int_lsq%lsq_pseudoinv(jc,2,1,jb) * z_b(jc,jk,1)  &
          &                 + ptr_int_lsq%lsq_pseudoinv(jc,2,2,jb) * z_b(jc,jk,2)  &
          &                 + ptr_int_lsq%lsq_pseudoinv(jc,2,3,jb) * z_b(jc,jk,3)

        ! zonal
        p_coeff(jc,jk,jb,2) = ptr_int_lsq%lsq_pseudoinv(jc,1,1,jb) * z_b(jc,jk,1)  &
          &                 + ptr_int_lsq%lsq_pseudoinv(jc,1,2,jb) * z_b(jc,jk,2)  &
          &                 + ptr_int_lsq%lsq_pseudoinv(jc,1,3,jb) * z_b(jc,jk,3)

        ! At the end, the coefficient c0 is derived from the linear constraint
        !
        ! constant
        p_coeff(jc,jk,jb,1) = p_cc(jc,jk,jb)                                         &
          &                 - p_coeff(jc,jk,jb,2) * ptr_int_lsq%lsq_moments(jc,jb,1) &
          &                 - p_coeff(jc,jk,jb,3) * ptr_int_lsq%lsq_moments(jc,jb,2)


      END DO ! end loop over cells
    END DO ! end loop over vertical levels

  END DO ! end loop over blocks
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE recon_lsq_cell_l_svd



!-------------------------------------------------------------------------
!
!
!>
!! Computes coefficients (i.e. derivatives) for cell centered quadratic
!! reconstruction.
!!
!! DESCRIPTION:
!! recon: reconstruction of subgrid distribution
!! lsq  : least-squares method
!! cell : solution coefficients defined at cell center
!! q    : quadratic reconstruction
!!
!! Computes the coefficients (derivatives) for a quadratic reconstruction,
!! using the the least-squares method. The coefficients are provided at
!! cell centers in a local 2D cartesian system (tangential plane).
!! Solves linear system Rx = Q^T d.
!! The matrices have the following size and shape:
!! R  : upper triangular matrix (5 x 5)
!! Q  : orthogonal matrix (9 x 5)
!! Q^T: transposed of Q (5 x 9)
!! d  : input vector (LHS) (9 x 1)
!! x  : solution vector (unknowns) (5 x 1)
!!
!! Coefficients
!! p_coeff(jc,jk,jb,1) : C0
!! p_coeff(jc,jk,jb,2) : C1 (dPhi_dx)
!! p_coeff(jc,jk,jb,3) : C2 (dPhi_dy)
!! p_coeff(jc,jk,jb,4) : C3 (0.5*ddPhi_ddx)
!! p_coeff(jc,jk,jb,5) : C4 (0.5*ddPhi_ddy)
!! p_coeff(jc,jk,jb,6) : C5 (ddPhi_dxdy)
!!
!! works only on triangular grid yet
!!
!! @par Revision History
!!  Developed and tested by Daniel Reinert, DWD (2009-11-13)
!! Modification by Daniel Reinert, DWD (2010-06-04)
!! - some speedup due to rearrangement and compiler directives
!!
!! !LITERATURE
!! Ollivier-Gooch et al (2002): A High-Order-Accurate Unstructured Mesh
!! Finite-Volume Scheme for the Advection-Diffusion Equation, J. Comput. Phys.,
!! 181, 729-752
!!
SUBROUTINE recon_lsq_cell_q( p_cc, ptr_patch, ptr_int_lsq, p_coeff, &
  &                               opt_slev, opt_elev, opt_rlstart,  &
  &                               opt_rlend )

  TYPE(t_patch), INTENT(IN) ::   & !< patch on which computation
    &  ptr_patch                   !< is performed
                                        
  TYPE(t_lsq), TARGET, INTENT(IN) :: &  !< data structure for interpolation
    &  ptr_int_lsq

  REAL(wp), INTENT(IN) ::           & !< cell centered variable
    &  p_cc(:,:,:)

  INTEGER, INTENT(IN), OPTIONAL ::  & !< optional vertical start level
    &  opt_slev

  INTEGER, INTENT(IN), OPTIONAL ::  & !< optional vertical end level
    &  opt_elev

  INTEGER, INTENT(IN), OPTIONAL ::  & !< start and end values of refin_ctrl flag
    &  opt_rlstart, opt_rlend

  REAL(wp), INTENT(INOUT)  ::  &  !< cell based coefficients (geographical components)
    &  p_coeff(:,:,:,:)           !< physically this vector contains gradients, second
                                  !< derivatives, one mixed derivative and a constant
                                  !< coefficient for zonal and meridional direction

  REAL(wp)  ::           &        !< difference of scalars i j
    &  z_d(lsq_high_set%dim_c,nproma,ptr_patch%nlev)
  REAL(wp)  ::           &        !< matrix-vector product of transposed
    &  z_qt_times_d(5)            !< Q matrix and d

  REAL(wp), POINTER ::   &        !< Pointer to reciprocal diagonal R-matrix-elements
    &  ptr_rrdiag(:,:,:)
  REAL(wp), POINTER ::   &        !< Pointer to upper triangular R-matrix-elements
    &  ptr_rutri(:,:,:)

  INTEGER, POINTER  ::   &        !< Pointer to line and block indices of
    &  iidx(:,:,:), iblk(:,:,:)   !< required stencil
  INTEGER :: slev, elev           !< vertical start and end level
  INTEGER :: jc, jk, jb           !< index of cell, vertical level and block
  INTEGER :: js                   !< loop index for cells in stencil
  INTEGER :: rl_start, rl_end
  INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx, i_nchdom

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
    elev = ptr_patch%nlev
  END IF
  IF ( PRESENT(opt_rlstart) ) THEN
    rl_start = opt_rlstart
  ELSE
    rl_start = 2
  END IF
  IF ( PRESENT(opt_rlend) ) THEN
    rl_end = opt_rlend
  ELSE
    rl_end = min_rlcell
  END IF


  ! values for the blocking
  i_nchdom   = MAX(1,ptr_patch%n_childdom)
  i_startblk = ptr_patch%cells%start_blk(rl_start,1)
  i_endblk   = ptr_patch%cells%end_blk(rl_end,i_nchdom)

  ! pointers to line and block indices of required stencil
  iidx => ptr_int_lsq%lsq_idx_c
  iblk => ptr_int_lsq%lsq_blk_c

  ! pointer to reciprocal diagonal R-elements
  ptr_rrdiag => ptr_int_lsq%lsq_rmat_rdiag_c(:,:,:)

  ! pointer to upper triangular R-elements
  ptr_rutri => ptr_int_lsq%lsq_rmat_utri_c(:,:,:)


  IF (ptr_patch%cell_type == 3) THEN


!$OMP PARALLEL

  IF (ptr_patch%id > 1) THEN
!$OMP WORKSHARE
    p_coeff(:,:,1:i_startblk,1:6) = 0._wp
!$OMP END WORKSHARE
  ENDIF

!$OMP DO PRIVATE(jb,jc,jk,js,i_startidx,i_endidx,z_d,z_qt_times_d), SCHEDULE(runtime)
  DO jb = i_startblk, i_endblk

    CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, rl_start, rl_end)

    !
    ! 1. compute right hand side of linear system
    !
#ifdef __LOOP_EXCHANGE
    DO jc = i_startidx, i_endidx
      DO jk = slev, elev
#else
!CDIR UNROLL=5
    DO jk = slev, elev
      DO jc = i_startidx, i_endidx
#endif

!CDIR EXPAND=9
        DO js = 1, 9
          ! original version
          !z_d(js) = ptr_int_lsq%lsq_weights_c(jc,js,jb)                   &
          !  &     * (p_cc(iidx(jc,jb,js),jk,iblk(jc,jb,js)) - p_cc(jc,jk,jb))

          ! note, that the multiplication with lsq_weights_c(jc,js,jb) at
          ! runtime is now avoided. Instead, the multiplication with
          ! lsq_weights_c(jc,js,jb) has been shifted into the transposed
          ! Q-matrix.

          z_d(js,jc,jk) = p_cc(iidx(jc,jb,js),jk,iblk(jc,jb,js)) - p_cc(jc,jk,jb)
        ENDDO
      ENDDO
    ENDDO



    !
    ! 2. compute cell based coefficients for quadratic reconstruction
    !
    DO jk = slev, elev

      DO jc = i_startidx, i_endidx

        ! calculate matrix vector product Q^T d (transposed of Q times LHS)
        ! (intrinsic function matmul not applied, due to massive
        ! performance penalty on the NEC. Instead the intrinsic dot product
        ! function is applied
!CDIR BEGIN EXPAND=9
        z_qt_times_d(1) = DOT_PRODUCT(ptr_int_lsq%lsq_qtmat_c(jc,1,1:9,jb),z_d(1:9,jc,jk))
        z_qt_times_d(2) = DOT_PRODUCT(ptr_int_lsq%lsq_qtmat_c(jc,2,1:9,jb),z_d(1:9,jc,jk))
        z_qt_times_d(3) = DOT_PRODUCT(ptr_int_lsq%lsq_qtmat_c(jc,3,1:9,jb),z_d(1:9,jc,jk))
        z_qt_times_d(4) = DOT_PRODUCT(ptr_int_lsq%lsq_qtmat_c(jc,4,1:9,jb),z_d(1:9,jc,jk))
        z_qt_times_d(5) = DOT_PRODUCT(ptr_int_lsq%lsq_qtmat_c(jc,5,1:9,jb),z_d(1:9,jc,jk))
!CDIR END

        !
        ! Solve linear system Rx=Q^T d by back substitution
        !
        p_coeff(jc,jk,jb,6) = ptr_rrdiag(jc,5,jb) * z_qt_times_d(5)
        p_coeff(jc,jk,jb,5) = ptr_rrdiag(jc,4,jb)                                         &
          &                 * ( z_qt_times_d(4) - ptr_rutri(jc,1,jb)*p_coeff(jc,jk,jb,6) )
        p_coeff(jc,jk,jb,4) = ptr_rrdiag(jc,3,jb)                                         &
          &                 * ( z_qt_times_d(3) - ptr_rutri(jc,2,jb)*p_coeff(jc,jk,jb,5)  &
          &                 - ptr_rutri(jc,3,jb) * p_coeff(jc,jk,jb,6) )
        p_coeff(jc,jk,jb,3) = ptr_rrdiag(jc,2,jb)                                         &
          &                 * ( z_qt_times_d(2) - ptr_rutri(jc,4,jb)*p_coeff(jc,jk,jb,4)  &
          &                 - ptr_rutri(jc,5,jb) * p_coeff(jc,jk,jb,5)                    &
          &                 - ptr_rutri(jc,6,jb) * p_coeff(jc,jk,jb,6) )
        p_coeff(jc,jk,jb,2) = ptr_rrdiag(jc,1,jb)                                         &
          &                 * ( z_qt_times_d(1) - ptr_rutri(jc,7,jb)*p_coeff(jc,jk,jb,3)  &
          &                 - ptr_rutri(jc,8,jb) * p_coeff(jc,jk,jb,4)                    &
          &                 - ptr_rutri(jc,9,jb) * p_coeff(jc,jk,jb,5)                    &
          &                 - ptr_rutri(jc,10,jb)* p_coeff(jc,jk,jb,6) )

        p_coeff(jc,jk,jb,1) = p_cc(jc,jk,jb)                                              &
          &                  - p_coeff(jc,jk,jb,2) * ptr_int_lsq%lsq_moments(jc,jb,1)     &
          &                  - p_coeff(jc,jk,jb,3) * ptr_int_lsq%lsq_moments(jc,jb,2)     &
          &                  - p_coeff(jc,jk,jb,4) * ptr_int_lsq%lsq_moments(jc,jb,3)     &
          &                  - p_coeff(jc,jk,jb,5) * ptr_int_lsq%lsq_moments(jc,jb,4)     &
          &                  - p_coeff(jc,jk,jb,6) * ptr_int_lsq%lsq_moments(jc,jb,5)


      END DO ! end loop over cells

    END DO ! end loop over vertical levels

  END DO ! end loop over blocks
!$OMP END DO
!$OMP END PARALLEL

ELSEIF (ptr_patch%cell_type == 6) THEN

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,jk,js,i_startidx,i_endidx,z_d,z_qt_times_d), SCHEDULE(runtime)
  DO jb = i_startblk, i_endblk

    CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, rl_start, rl_end)

    !
    ! 1. compute right hand side of linear system
    !
!CDIR UNROLL=5
    DO jk = slev, elev

      DO jc = i_startidx, i_endidx

!CDIR EXPAND=6
        DO js = 1, 6
          ! original version
          !z_d(js) = ptr_int_lsq%lsq_weights_c(jc,js,jb)                   &
          !  &     * (p_cc(iidx(jc,jb,js),jk,iblk(jc,jb,js)) - p_cc(jc,jk,jb))

          ! note, that the multiplication with lsq_weights_c(jc,js,jb) at
          ! runtime is now avoided. Instead, the multiplication with
          ! lsq_weights_c(jc,js,jb) has been shifted into the transposed
          ! Q-matrix.

          z_d(js,jc,jk) = p_cc(iidx(jc,jb,js),jk,iblk(jc,jb,js)) - p_cc(jc,jk,jb)
        ENDDO
      ENDDO
    ENDDO

    !
    ! 2. compute cell based coefficients for quadratic reconstruction
    !
    DO jk = slev, elev

      DO jc = i_startidx, i_endidx

        ! calculate matrix vector product Q^T d (transposed of Q times LHS)
        ! (intrinsic function matmul not applied, due to massive
        ! performance penalty on the NEC. Instead the intrinsic dot product
        ! function is applied
!CDIR BEGIN EXPAND=6
!!        z_qt_times_d(1) = DOT_PRODUCT(ptr_int_lsq%lsq_qtmat_c(jc,1,1:6,jb),z_d(1:6,jc,jk))
!!        z_qt_times_d(2) = DOT_PRODUCT(ptr_int_lsq%lsq_qtmat_c(jc,2,1:6,jb),z_d(1:6,jc,jk))
        z_qt_times_d(3) = DOT_PRODUCT(ptr_int_lsq%lsq_qtmat_c(jc,3,1:6,jb),z_d(1:6,jc,jk))
        z_qt_times_d(4) = DOT_PRODUCT(ptr_int_lsq%lsq_qtmat_c(jc,4,1:6,jb),z_d(1:6,jc,jk))
        z_qt_times_d(5) = DOT_PRODUCT(ptr_int_lsq%lsq_qtmat_c(jc,5,1:6,jb),z_d(1:6,jc,jk))
!CDIR END

        !
        ! Solve linear system Rx=Q^T d by back substitution
        !
        p_coeff(jc,jk,jb,6) = ptr_rrdiag(jc,5,jb) * z_qt_times_d(5)
        p_coeff(jc,jk,jb,5) = ptr_rrdiag(jc,4,jb)                                         &
          &                 * ( z_qt_times_d(4) - ptr_rutri(jc,1,jb)*p_coeff(jc,jk,jb,6) )
        p_coeff(jc,jk,jb,4) = ptr_rrdiag(jc,3,jb)                                         &
          &                 * ( z_qt_times_d(3) - ptr_rutri(jc,2,jb)*p_coeff(jc,jk,jb,5)  &
          &                 - ptr_rutri(jc,3,jb) * p_coeff(jc,jk,jb,6) )
        !p_coeff(jc,jk,jb,3) = ptr_rrdiag(jc,2,jb)                                         &
        !  &                 * ( z_qt_times_d(2) - ptr_rutri(jc,4,jb)*p_coeff(jc,jk,jb,4)  &
        !  &                 - ptr_rutri(jc,5,jb) * p_coeff(jc,jk,jb,5)                    &
        !  &                 - ptr_rutri(jc,6,jb) * p_coeff(jc,jk,jb,6) )
        !p_coeff(jc,jk,jb,2) = ptr_rrdiag(jc,1,jb)                                         &
        !  &                 * ( z_qt_times_d(1) - ptr_rutri(jc,7,jb)*p_coeff(jc,jk,jb,3)  &
        !  &                 - ptr_rutri(jc,8,jb) * p_coeff(jc,jk,jb,4)                    &
        !  &                 - ptr_rutri(jc,9,jb) * p_coeff(jc,jk,jb,5)                    &
        !  &                 - ptr_rutri(jc,10,jb)* p_coeff(jc,jk,jb,6) )

        !p_coeff(jc,jk,jb,1) = p_cc(jc,jk,jb)                                              &
        !  &                  - p_coeff(jc,jk,jb,2) * ptr_int_lsq%lsq_moments(jc,jb,1)     &
        !  &                  - p_coeff(jc,jk,jb,3) * ptr_int_lsq%lsq_moments(jc,jb,2)     &
        !  &                  - p_coeff(jc,jk,jb,4) * ptr_int_lsq%lsq_moments(jc,jb,3)     &
        !  &                  - p_coeff(jc,jk,jb,5) * ptr_int_lsq%lsq_moments(jc,jb,4)     &
        !  &                  - p_coeff(jc,jk,jb,6) * ptr_int_lsq%lsq_moments(jc,jb,5)

      END DO ! end loop over cells

    END DO ! end loop over vertical levels

  END DO ! end loop over blocks
!$OMP END DO
!$OMP END PARALLEL

ENDIF

END SUBROUTINE recon_lsq_cell_q



!-------------------------------------------------------------------------
!
!
!>
!! Computes coefficients (i.e. derivatives) for cell centered quadratic
!! reconstruction.
!!
!! DESCRIPTION:
!! recon: reconstruction of subgrid distribution
!! lsq  : least-squares method
!! cell : solution coefficients defined at cell center
!! q    : quadratic reconstruction
!!
!! Computes unknown coefficients (derivatives) of a quadratic polynomial,
!! using the least-squares method. The coefficients are provided at cell 
!! centers in a local 2D cartesian system (tangential plane).
!!
!! Mathematically we solve Ax = b via Singular Value Decomposition (SVD)
!! x = PINV(A) * b
!!
!! Matrices have the following size and shape (triangular grid) :
!! PINV(A): Pseudo or Moore-Penrose inverse of A (via SVD) (5 x 9)
!! b  : input vector (LHS) (9 x 1)
!! x  : solution vector (unknowns) (5 x 1)
!!
!! Coefficients:
!! p_coeff(jc,jk,jb,1) : C0
!! p_coeff(jc,jk,jb,2) : C1 (dPhi_dx)
!! p_coeff(jc,jk,jb,3) : C2 (dPhi_dy)
!! p_coeff(jc,jk,jb,4) : C3 (0.5*ddPhi_ddx)
!! p_coeff(jc,jk,jb,5) : C4 (0.5*ddPhi_ddy)
!! p_coeff(jc,jk,jb,6) : C5 (ddPhi_dxdy)
!!
!!
!! @par Revision History
!!  Developed and tested by Daniel Reinert, DWD (2011-05-31)
!!
!! !LITERATURE
!! Ollivier-Gooch et al (2002): A High-Order-Accurate Unstructured Mesh
!! Finite-Volume Scheme for the Advection-Diffusion Equation, J. Comput. Phys.,
!! 181, 729-752
!!
SUBROUTINE recon_lsq_cell_q_svd( p_cc, ptr_patch, ptr_int_lsq, p_coeff, &
  &                               opt_slev, opt_elev, opt_rlstart,      &
  &                               opt_rlend )

  TYPE(t_patch), INTENT(IN) ::   & !< patch on which computation
    &  ptr_patch                   !< is performed
                                        
  TYPE(t_lsq), TARGET, INTENT(IN) :: &  !< data structure for interpolation
    &  ptr_int_lsq

  REAL(wp), INTENT(IN) ::           & !< cell centered variable
    &  p_cc(:,:,:)

  INTEGER, INTENT(IN), OPTIONAL ::  & !< optional vertical start level
    &  opt_slev

  INTEGER, INTENT(IN), OPTIONAL ::  & !< optional vertical end level
    &  opt_elev

  INTEGER, INTENT(IN), OPTIONAL ::  & !< start and end values of refin_ctrl flag
    &  opt_rlstart, opt_rlend

  REAL(wp), INTENT(INOUT)  ::  &  !< cell based coefficients (geographical components)
    &  p_coeff(:,:,:,:)           !< physically this vector contains gradients, second
                                  !< derivatives, one mixed derivative and a constant
                                  !< coefficient for zonal and meridional direction

  REAL(wp)  ::           &        !< difference of scalars i j
    &  z_b(lsq_high_set%dim_c,nproma,ptr_patch%nlev)

  INTEGER, POINTER  ::   &        !< Pointer to line and block indices of
    &  iidx(:,:,:), iblk(:,:,:)   !< required stencil
  INTEGER :: slev, elev           !< vertical start and end level
  INTEGER :: jc, jk, jb           !< index of cell, vertical level and block
  INTEGER :: js                   !< loop index for cells in stencil
  INTEGER :: rl_start, rl_end
  INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx, i_nchdom

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
    elev = ptr_patch%nlev
  END IF
  IF ( PRESENT(opt_rlstart) ) THEN
    rl_start = opt_rlstart
  ELSE
    rl_start = 2
  END IF
  IF ( PRESENT(opt_rlend) ) THEN
    rl_end = opt_rlend
  ELSE
    rl_end = min_rlcell
  END IF


  ! values for the blocking
  i_nchdom   = MAX(1,ptr_patch%n_childdom)
  i_startblk = ptr_patch%cells%start_blk(rl_start,1)
  i_endblk   = ptr_patch%cells%end_blk(rl_end,i_nchdom)

  ! pointers to line and block indices of required stencil
  iidx => ptr_int_lsq%lsq_idx_c
  iblk => ptr_int_lsq%lsq_blk_c


  IF (ptr_patch%cell_type == 3) THEN


!$OMP PARALLEL

  IF (ptr_patch%id > 1) THEN
!$OMP WORKSHARE
    p_coeff(:,:,1:i_startblk,1:6) = 0._wp
!$OMP END WORKSHARE
  ENDIF

!$OMP DO PRIVATE(jb,jc,jk,js,i_startidx,i_endidx,z_b), SCHEDULE(runtime)
  DO jb = i_startblk, i_endblk

    CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, rl_start, rl_end)

    !
    ! 1. compute right hand side of linear system
    !
#ifdef __LOOP_EXCHANGE
    DO jc = i_startidx, i_endidx
      DO jk = slev, elev
#else
!CDIR UNROLL=5
    DO jk = slev, elev
      DO jc = i_startidx, i_endidx
#endif

!CDIR EXPAND=9
        DO js = 1, 9
          ! original version
          !z_b(js) = ptr_int_lsq%lsq_weights_c(jc,js,jb)                   &
          !  &     * (p_cc(iidx(jc,jb,js),jk,iblk(jc,jb,js)) - p_cc(jc,jk,jb))

          ! note that multiplication with lsq_weights_c(jc,js,jb) during 
          ! runtime is now avoided. Instead, the weights have been shifted 
          ! into the pseudoinverse.
          z_b(js,jc,jk) = p_cc(iidx(jc,jb,js),jk,iblk(jc,jb,js)) - p_cc(jc,jk,jb)
        ENDDO
      ENDDO
    ENDDO



    !
    ! 2. compute cell based coefficients for quadratic reconstruction
    !    calculate matrix vector product PINV(A) * b
    !
    DO jk = slev, elev

      DO jc = i_startidx, i_endidx

        ! (intrinsic function matmul not applied, due to massive
        ! performance penalty on the NEC. Instead the intrinsic dot product
        ! function is used.
!CDIR BEGIN EXPAND=9
        p_coeff(jc,jk,jb,6) = DOT_PRODUCT(ptr_int_lsq%lsq_pseudoinv(jc,5,1:9,jb), &
          &                               z_b(1:9,jc,jk))
        p_coeff(jc,jk,jb,5) = DOT_PRODUCT(ptr_int_lsq%lsq_pseudoinv(jc,4,1:9,jb), &
          &                               z_b(1:9,jc,jk))
        p_coeff(jc,jk,jb,4) = DOT_PRODUCT(ptr_int_lsq%lsq_pseudoinv(jc,3,1:9,jb), &
          &                               z_b(1:9,jc,jk))
        p_coeff(jc,jk,jb,3) = DOT_PRODUCT(ptr_int_lsq%lsq_pseudoinv(jc,2,1:9,jb), &
          &                               z_b(1:9,jc,jk))
        p_coeff(jc,jk,jb,2) = DOT_PRODUCT(ptr_int_lsq%lsq_pseudoinv(jc,1,1:9,jb), &
          &                               z_b(1:9,jc,jk))
!CDIR END


        ! At the end, the coefficient c0 is derived from the linear constraint
        !
        p_coeff(jc,jk,jb,1) = p_cc(jc,jk,jb) - DOT_PRODUCT(p_coeff(jc,jk,jb,2:6), &
          &                   ptr_int_lsq%lsq_moments(jc,jb,1:5))


      END DO ! end loop over cells

    END DO ! end loop over vertical levels

  END DO ! end loop over blocks
!$OMP END DO
!$OMP END PARALLEL

ELSEIF (ptr_patch%cell_type == 6) THEN

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,jk,js,i_startidx,i_endidx,z_b), SCHEDULE(runtime)
  DO jb = i_startblk, i_endblk

    CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, rl_start, rl_end)

    !
    ! 1. compute right hand side of linear system
    !
!CDIR UNROLL=5
    DO jk = slev, elev

      DO jc = i_startidx, i_endidx

!CDIR EXPAND=6
        DO js = 1, 6
          ! original version
          !z_b(js) = ptr_int_lsq%lsq_weights_c(jc,js,jb)                   &
          !  &     * (p_cc(iidx(jc,jb,js),jk,iblk(jc,jb,js)) - p_cc(jc,jk,jb))

          ! note that the multiplication with lsq_weights_c(jc,js,jb) at
          ! runtime is now avoided. Instead, the multiplication has been
          ! shifted into the pseudoinverse.

          z_b(js,jc,jk) = p_cc(iidx(jc,jb,js),jk,iblk(jc,jb,js)) - p_cc(jc,jk,jb)
        ENDDO
      ENDDO
    ENDDO

    !
    ! 2. compute cell based coefficients for quadratic reconstruction
    !    calculate matrix vector product PINV(A) * b
    !
    DO jk = slev, elev

      DO jc = i_startidx, i_endidx

        ! (intrinsic function matmul not applied, due to massive
        ! performance penalty on the NEC. Instead the intrinsic dot product
        ! function is applied
!CDIR BEGIN EXPAND=6
        p_coeff(jc,jk,jb,6) = DOT_PRODUCT(ptr_int_lsq%lsq_pseudoinv(jc,5,1:6,jb), &
          &                               z_b(1:6,jc,jk))
        p_coeff(jc,jk,jb,5) = DOT_PRODUCT(ptr_int_lsq%lsq_pseudoinv(jc,4,1:6,jb), &
          &                               z_b(1:6,jc,jk))
        p_coeff(jc,jk,jb,4) = DOT_PRODUCT(ptr_int_lsq%lsq_pseudoinv(jc,3,1:6,jb), &
          &                               z_b(1:6,jc,jk))

        !p_coeff(jc,jk,jb,3) = DOT_PRODUCT(ptr_int_lsq%lsq_pseudoinv(jc,2,1:6,jb), &
        !  &                               z_b(1:6,jc,jk))
        !p_coeff(jc,jk,jb,2) = DOT_PRODUCT(ptr_int_lsq%lsq_pseudoinv(jc,1,1:6,jb), &
        !  &                               z_b(1:6,jc,jk))
!CDIR END

        !p_coeff(jc,jk,jb,1) = p_cc(jc,jk,jb) - DOT_PRODUCT(p_coeff(jc,jk,jb,2:6), &
        !  &                   ptr_int_lsq%lsq_moments(jc,jb,1:5))


      END DO ! end loop over cells

    END DO ! end loop over vertical levels

  END DO ! end loop over blocks
!$OMP END DO
!$OMP END PARALLEL

ENDIF

END SUBROUTINE recon_lsq_cell_q_svd


!-------------------------------------------------------------------------
!
!
!>
!! Computes coefficients (i.e. derivatives) for cell centered poor man's
!! cubic reconstruction.
!!
!! DESCRIPTION:
!! recon: reconstruction of subgrid distribution
!! lsq  : least-squares method
!! cell : solution coefficients defined at cell center
!! c    : cubic reconstruction
!!
!! Computes the coefficients (derivatives) for a cubic reconstruction
!! without cross derivatives, using the least-squares method. The
!! coefficients are provided at cell centers in a local 2D cartesian
!! system (tangential plane).
!! Solves linear system Rx = Q^T d.
!! The matrices have the following size and shape:
!! R  : upper triangular matrix (7 x 7)
!! Q  : orthogonal matrix (9 x 7)
!! Q^T: transposed of Q (7 x 9)
!! d  : input vector (LHS) (9 x 1)
!! x  : solution vector (unknowns) (7 x 1)
!!
!! Coefficients
!! p_coeff(jc,jk,jb,1) : C0
!! p_coeff(jc,jk,jb,2) : C1 (dPhi_dx)
!! p_coeff(jc,jk,jb,3) : C2 (dPhi_dy)
!! p_coeff(jc,jk,jb,4) : C3 (1/2*ddPhi_ddx)
!! p_coeff(jc,jk,jb,5) : C4 (1/2*ddPhi_ddy)
!! p_coeff(jc,jk,jb,6) : C5 (ddPhi_dxdy)
!! p_coeff(jc,jk,jb,7) : C6 (1/6*dddPhi_dddx)
!! p_coeff(jc,jk,jb,8) : C7 (1/6*dddPhi_dddy)
!!
!! works only on triangular grid yet
!!
!! @par Revision History
!!  Developed and tested by Daniel Reinert, DWD (2010-10-14)
!!
!! !LITERATURE
!! Ollivier-Gooch et al (2002): A High-Order-Accurate Unstructured Mesh
!! Finite-Volume Scheme for the Advection-Diffusion Equation, J. Comput. Phys.,
!! 181, 729-752
!!
SUBROUTINE recon_lsq_cell_cpoor( p_cc, ptr_patch, ptr_int_lsq, p_coeff, &
  &                               opt_slev, opt_elev,               &
  &                               opt_rlstart, opt_rlend )

  TYPE(t_patch), INTENT(IN) :: & !< patch on which computation
    &  ptr_patch                 !< is performed
                                        
  TYPE(t_lsq), TARGET, INTENT(IN) :: &  !< data structure for interpolation
    &  ptr_int_lsq

  REAL(wp), INTENT(IN) ::           & !< cell centered variable
    &  p_cc(:,:,:)

  INTEGER, INTENT(IN), OPTIONAL ::  & !< optional vertical start level
    &  opt_slev

  INTEGER, INTENT(IN), OPTIONAL ::  & !< optional vertical end level
    &  opt_elev

  INTEGER, INTENT(IN), OPTIONAL ::  & !< start and end values of refin_ctrl flag
    &  opt_rlstart, opt_rlend

  REAL(wp), INTENT(INOUT)  ::  &  !< cell based coefficients (geographical components)
    &  p_coeff(:,:,:,:)           !< physically this vector contains gradients, second
                                  !< and third derivatives, one mixed derivative and a
                                  !< constant coefficient for zonal and meridional
                                  !< direction

  REAL(wp)  ::           &        !< difference of scalars i j
    &  z_d(lsq_high_set%dim_c,nproma,ptr_patch%nlev)
  REAL(wp)  ::           &        !< matrix-vector product of transposed
    &  z_qt_times_d(7)            !< Q matrix and d

  REAL(wp), POINTER ::   &        !< Pointer to reciprocal diagonal R-matrix-elements
    &  ptr_rrdiag(:,:,:)
  REAL(wp), POINTER ::   &        !< Pointer to upper triangular R-matrix-elements
    &  ptr_rutri(:,:,:)

  INTEGER, POINTER  ::   &        !< Pointer to line and block indices of
    &  iidx(:,:,:), iblk(:,:,:)   !< required stencil
  INTEGER :: slev, elev           !< vertical start and end level
  INTEGER :: jc, jk, jb           !< index of cell, vertical level and block
  INTEGER :: js                   !< loop index for cells in stencil
  INTEGER :: rl_start, rl_end
  INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx, i_nchdom

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
    elev = ptr_patch%nlev
  END IF
  IF ( PRESENT(opt_rlstart) ) THEN
    rl_start = opt_rlstart
  ELSE
    rl_start = 2
  END IF
  IF ( PRESENT(opt_rlend) ) THEN
    rl_end = opt_rlend
  ELSE
    rl_end = min_rlcell
  END IF


  ! values for the blocking
  i_nchdom   = MAX(1,ptr_patch%n_childdom)
  i_startblk = ptr_patch%cells%start_blk(rl_start,1)
  i_endblk   = ptr_patch%cells%end_blk(rl_end,i_nchdom)

  ! pointers to line and block indices of required stencil
  iidx => ptr_int_lsq%lsq_idx_c
  iblk => ptr_int_lsq%lsq_blk_c

  ! pointer to reciprocal diagonal R-elements
  ptr_rrdiag => ptr_int_lsq%lsq_rmat_rdiag_c(:,:,:)

  ! pointer to upper triangular R-elements
  ptr_rutri => ptr_int_lsq%lsq_rmat_utri_c(:,:,:)



!$OMP PARALLEL

  IF (ptr_patch%id > 1) THEN
!$OMP WORKSHARE
    p_coeff(:,:,1:i_startblk,1:8) = 0._wp
!$OMP END WORKSHARE
  ENDIF

!$OMP DO PRIVATE(jb,jc,jk,js,i_startidx,i_endidx,z_d,z_qt_times_d), SCHEDULE(runtime)
  DO jb = i_startblk, i_endblk

    CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, rl_start, rl_end)

    !
    ! 1. compute right hand side of linear system
    !
#ifdef __LOOP_EXCHANGE
    DO jc = i_startidx, i_endidx
      DO jk = slev, elev
#else
!CDIR UNROLL=5
    DO jk = slev, elev
      DO jc = i_startidx, i_endidx
#endif

!CDIR EXPAND=9
        DO js = 1, 9
          ! original version
          !z_d(js) = ptr_int_lsq%lsq_weights_c(jc,js,jb)                   &
          !  &     * (p_cc(iidx(jc,jb,js),jk,iblk(jc,jb,js)) - p_cc(jc,jk,jb))

          ! note, that the multiplication with lsq_weights_c(jc,js,jb) at
          ! runtime is now avoided. Instead, the multiplication with
          ! lsq_weights_c(jc,js,jb) has been shifted into the transposed
          ! Q-matrix.

          z_d(js,jc,jk) = p_cc(iidx(jc,jb,js),jk,iblk(jc,jb,js)) - p_cc(jc,jk,jb)
        ENDDO
      ENDDO
    ENDDO



    !
    ! 2. compute cell based coefficients for quadratic reconstruction
    !
    DO jk = slev, elev

      DO jc = i_startidx, i_endidx

        ! calculate matrix vector product Q^T d (transposed of Q times LHS)
        ! (intrinsic function matmul not applied, due to massive
        ! performance penalty on the NEC. Instead the intrinsic dot product
        ! function is applied
!CDIR BEGIN EXPAND=9
        z_qt_times_d(1) = DOT_PRODUCT(ptr_int_lsq%lsq_qtmat_c(jc,1,1:9,jb),z_d(1:9,jc,jk))
        z_qt_times_d(2) = DOT_PRODUCT(ptr_int_lsq%lsq_qtmat_c(jc,2,1:9,jb),z_d(1:9,jc,jk))
        z_qt_times_d(3) = DOT_PRODUCT(ptr_int_lsq%lsq_qtmat_c(jc,3,1:9,jb),z_d(1:9,jc,jk))
        z_qt_times_d(4) = DOT_PRODUCT(ptr_int_lsq%lsq_qtmat_c(jc,4,1:9,jb),z_d(1:9,jc,jk))
        z_qt_times_d(5) = DOT_PRODUCT(ptr_int_lsq%lsq_qtmat_c(jc,5,1:9,jb),z_d(1:9,jc,jk))
        z_qt_times_d(6) = DOT_PRODUCT(ptr_int_lsq%lsq_qtmat_c(jc,6,1:9,jb),z_d(1:9,jc,jk))
        z_qt_times_d(7) = DOT_PRODUCT(ptr_int_lsq%lsq_qtmat_c(jc,7,1:9,jb),z_d(1:9,jc,jk))
!CDIR END

        !
        ! Solve linear system Rx=Q^T d by back substitution
        !
        p_coeff(jc,jk,jb,8) = ptr_rrdiag(jc,7,jb) * z_qt_times_d(7)
        p_coeff(jc,jk,jb,7) = ptr_rrdiag(jc,6,jb)                                         &
          &                 * ( z_qt_times_d(6) - ptr_rutri(jc,1,jb)*p_coeff(jc,jk,jb,8) )
        p_coeff(jc,jk,jb,6) = ptr_rrdiag(jc,5,jb)                                         &
          &                 * ( z_qt_times_d(5) - ptr_rutri(jc,2,jb)*p_coeff(jc,jk,jb,7)  &
          &                 - ptr_rutri(jc,3,jb) * p_coeff(jc,jk,jb,8) )
        p_coeff(jc,jk,jb,5) = ptr_rrdiag(jc,4,jb)                                         &
          &                 * ( z_qt_times_d(4) - ptr_rutri(jc,4,jb)*p_coeff(jc,jk,jb,6)  &
          &                 - ptr_rutri(jc,5,jb) * p_coeff(jc,jk,jb,7)                    &
          &                 - ptr_rutri(jc,6,jb) * p_coeff(jc,jk,jb,8) )
        p_coeff(jc,jk,jb,4) = ptr_rrdiag(jc,3,jb)                                         &
          &                 * ( z_qt_times_d(3) - ptr_rutri(jc,7,jb)*p_coeff(jc,jk,jb,5)  &
          &                 - ptr_rutri(jc,8,jb) * p_coeff(jc,jk,jb,6)                    &
          &                 - ptr_rutri(jc,9,jb) * p_coeff(jc,jk,jb,7)                    &
          &                 - ptr_rutri(jc,10,jb)* p_coeff(jc,jk,jb,8) )
        p_coeff(jc,jk,jb,3) = ptr_rrdiag(jc,2,jb)                                         &
          &                 * ( z_qt_times_d(2) - ptr_rutri(jc,11,jb)*p_coeff(jc,jk,jb,4) &
          &                 - ptr_rutri(jc,12,jb) * p_coeff(jc,jk,jb,5)                   &
          &                 - ptr_rutri(jc,13,jb) * p_coeff(jc,jk,jb,6)                   &
          &                 - ptr_rutri(jc,14,jb) * p_coeff(jc,jk,jb,7)                   &
          &                 - ptr_rutri(jc,15,jb) * p_coeff(jc,jk,jb,8) )
        p_coeff(jc,jk,jb,2) = ptr_rrdiag(jc,1,jb)                                         &
          &                 * ( z_qt_times_d(1) - ptr_rutri(jc,16,jb)*p_coeff(jc,jk,jb,3) &
          &                 - ptr_rutri(jc,17,jb) * p_coeff(jc,jk,jb,4)                   &
          &                 - ptr_rutri(jc,18,jb) * p_coeff(jc,jk,jb,5)                   &
          &                 - ptr_rutri(jc,19,jb) * p_coeff(jc,jk,jb,6)                   &
          &                 - ptr_rutri(jc,20,jb) * p_coeff(jc,jk,jb,7)                   &
          &                 - ptr_rutri(jc,21,jb) * p_coeff(jc,jk,jb,8) )


        p_coeff(jc,jk,jb,1) = p_cc(jc,jk,jb)                                              &
          &                  - p_coeff(jc,jk,jb,2) * ptr_int_lsq%lsq_moments(jc,jb,1)     &
          &                  - p_coeff(jc,jk,jb,3) * ptr_int_lsq%lsq_moments(jc,jb,2)     &
          &                  - p_coeff(jc,jk,jb,4) * ptr_int_lsq%lsq_moments(jc,jb,3)     &
          &                  - p_coeff(jc,jk,jb,5) * ptr_int_lsq%lsq_moments(jc,jb,4)     &
          &                  - p_coeff(jc,jk,jb,6) * ptr_int_lsq%lsq_moments(jc,jb,5)     &
          &                  - p_coeff(jc,jk,jb,7) * ptr_int_lsq%lsq_moments(jc,jb,6)     &
          &                  - p_coeff(jc,jk,jb,8) * ptr_int_lsq%lsq_moments(jc,jb,7)


      END DO ! end loop over cells

    END DO ! end loop over vertical levels

  END DO ! end loop over blocks
!$OMP END DO
!$OMP END PARALLEL


END SUBROUTINE recon_lsq_cell_cpoor



!--------------------------------------------------
!
!
!>
!! Computes coefficients (i.e. derivatives) for cell centered poor man's
!! cubic reconstruction.
!!
!! DESCRIPTION:
!! recon: reconstruction of subgrid distribution
!! lsq  : least-squares method
!! cell : solution coefficients defined at cell center
!! c    : cubic reconstruction
!!
!! Computes unknown coefficients (derivatives) of a cubic polynomial,
!! without cross derivatives, using the least-squares method. The 
!! coefficients are provided at cell centers in a local 2D cartesian 
!! system (tangential plane).
!!
!! Mathematically we solve Ax = b via Singular Value Decomposition (SVD)
!! x = PINV(A) * b
!!
!! Matrices have the following size and shape (triangular grid) :
!! PINV(A): Pseudo or Moore-Penrose inverse of A (via SVD) (7 x 9)
!! b  : input vector (LHS) (9 x 1)
!! x  : solution vector (unknowns) (7 x 1)
!!
!! Coefficients
!! p_coeff(jc,jk,jb,1) : C0
!! p_coeff(jc,jk,jb,2) : C1 (dPhi_dx)
!! p_coeff(jc,jk,jb,3) : C2 (dPhi_dy)
!! p_coeff(jc,jk,jb,4) : C3 (1/2*ddPhi_ddx)
!! p_coeff(jc,jk,jb,5) : C4 (1/2*ddPhi_ddy)
!! p_coeff(jc,jk,jb,6) : C5 (ddPhi_dxdy)
!! p_coeff(jc,jk,jb,7) : C6 (1/6*dddPhi_dddx)
!! p_coeff(jc,jk,jb,8) : C7 (1/6*dddPhi_dddy)
!!
!! works only on triangular grid yet
!!
!! @par Revision History
!!  Developed and tested by Daniel Reinert, DWD (2011-05-31)
!!
!! !LITERATURE
!! Ollivier-Gooch et al (2002): A High-Order-Accurate Unstructured Mesh
!! Finite-Volume Scheme for the Advection-Diffusion Equation, J. Comput. Phys.,
!! 181, 729-752
!!
SUBROUTINE recon_lsq_cell_cpoor_svd( p_cc, ptr_patch, ptr_int_lsq, p_coeff, &
  &                               opt_slev, opt_elev, opt_rlstart, opt_rlend )

  TYPE(t_patch), INTENT(IN) :: & !< patch on which computation
    &  ptr_patch                 !< is performed
                                        
  TYPE(t_lsq), TARGET, INTENT(IN) :: &  !< data structure for interpolation
    &  ptr_int_lsq

  REAL(wp), INTENT(IN) ::           & !< cell centered variable
    &  p_cc(:,:,:)

  INTEGER, INTENT(IN), OPTIONAL ::  & !< optional vertical start level
    &  opt_slev

  INTEGER, INTENT(IN), OPTIONAL ::  & !< optional vertical end level
    &  opt_elev

  INTEGER, INTENT(IN), OPTIONAL ::  & !< start and end values of refin_ctrl flag
    &  opt_rlstart, opt_rlend

  REAL(wp), INTENT(INOUT)  ::  &  !< cell based coefficients (geographical components)
    &  p_coeff(:,:,:,:)           !< physically this vector contains gradients, second
                                  !< and third derivatives, one mixed derivative and a
                                  !< constant coefficient for zonal and meridional
                                  !< direction

  REAL(wp)  ::           &        !< difference of scalars i j
    &  z_b(lsq_high_set%dim_c,nproma,ptr_patch%nlev)

  INTEGER, POINTER  ::   &        !< Pointer to line and block indices of
    &  iidx(:,:,:), iblk(:,:,:)   !< required stencil
  INTEGER :: slev, elev           !< vertical start and end level
  INTEGER :: jc, jk, jb           !< index of cell, vertical level and block
  INTEGER :: js                   !< loop index for cells in stencil
  INTEGER :: rl_start, rl_end
  INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx, i_nchdom

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
    elev = ptr_patch%nlev
  END IF
  IF ( PRESENT(opt_rlstart) ) THEN
    rl_start = opt_rlstart
  ELSE
    rl_start = 2
  END IF
  IF ( PRESENT(opt_rlend) ) THEN
    rl_end = opt_rlend
  ELSE
    rl_end = min_rlcell
  END IF


  ! values for the blocking
  i_nchdom   = MAX(1,ptr_patch%n_childdom)
  i_startblk = ptr_patch%cells%start_blk(rl_start,1)
  i_endblk   = ptr_patch%cells%end_blk(rl_end,i_nchdom)

  ! pointers to line and block indices of required stencil
  iidx => ptr_int_lsq%lsq_idx_c
  iblk => ptr_int_lsq%lsq_blk_c


!$OMP PARALLEL

  IF (ptr_patch%id > 1) THEN
!$OMP WORKSHARE
    p_coeff(:,:,1:i_startblk,1:8) = 0._wp
!$OMP END WORKSHARE
  ENDIF

!$OMP DO PRIVATE(jb,jc,jk,js,i_startidx,i_endidx,z_b), SCHEDULE(runtime)
  DO jb = i_startblk, i_endblk

    CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, rl_start, rl_end)

    !
    ! 1. compute right hand side of linear system
    !
#ifdef __LOOP_EXCHANGE
    DO jc = i_startidx, i_endidx
      DO jk = slev, elev
#else
!CDIR UNROLL=5
    DO jk = slev, elev
      DO jc = i_startidx, i_endidx
#endif

!CDIR EXPAND=9
        DO js = 1, 9
          ! original version
          !z_b(js) = ptr_int_lsq%lsq_weights_c(jc,js,jb)                   &
          !  &     * (p_cc(iidx(jc,jb,js),jk,iblk(jc,jb,js)) - p_cc(jc,jk,jb))

          ! note that multiplication with lsq_weights_c(jc,js,jb) during 
          ! runtime is now avoided. Instead, the weights have been shifted 
          ! into the pseudoinverse.

          z_b(js,jc,jk) = p_cc(iidx(jc,jb,js),jk,iblk(jc,jb,js)) - p_cc(jc,jk,jb)
        ENDDO
      ENDDO
    ENDDO



    !
    ! 2. compute cell based coefficients for poor man's cubic reconstruction
    !    calculate matrix vector product PINV(A) * b
    !
    DO jk = slev, elev

      DO jc = i_startidx, i_endidx

        ! (intrinsic function matmul not applied, due to massive
        ! performance penalty on the NEC. Instead the intrinsic dot product
        ! function is used.
!CDIR BEGIN EXPAND=9
        p_coeff(jc,jk,jb,8) = DOT_PRODUCT(ptr_int_lsq%lsq_pseudoinv(jc,7,1:9,jb), &
          &                               z_b(1:9,jc,jk))
        p_coeff(jc,jk,jb,7) = DOT_PRODUCT(ptr_int_lsq%lsq_pseudoinv(jc,6,1:9,jb), &
          &                               z_b(1:9,jc,jk))
        p_coeff(jc,jk,jb,6) = DOT_PRODUCT(ptr_int_lsq%lsq_pseudoinv(jc,5,1:9,jb), &
          &                               z_b(1:9,jc,jk))
        p_coeff(jc,jk,jb,5) = DOT_PRODUCT(ptr_int_lsq%lsq_pseudoinv(jc,4,1:9,jb), &
          &                               z_b(1:9,jc,jk))
        p_coeff(jc,jk,jb,4) = DOT_PRODUCT(ptr_int_lsq%lsq_pseudoinv(jc,3,1:9,jb), &
          &                               z_b(1:9,jc,jk))
        p_coeff(jc,jk,jb,3) = DOT_PRODUCT(ptr_int_lsq%lsq_pseudoinv(jc,2,1:9,jb), &
          &                               z_b(1:9,jc,jk))
        p_coeff(jc,jk,jb,2) = DOT_PRODUCT(ptr_int_lsq%lsq_pseudoinv(jc,1,1:9,jb), &
          &                               z_b(1:9,jc,jk))
!CDIR END

        ! At the end, the coefficient c0 is derived from the linear constraint
        !
!CDIR EXPAND=7
        p_coeff(jc,jk,jb,1) = p_cc(jc,jk,jb) - DOT_PRODUCT(p_coeff(jc,jk,jb,2:8), &
          &                   ptr_int_lsq%lsq_moments(jc,jb,1:7))

      END DO ! end loop over cells

    END DO ! end loop over vertical levels

  END DO ! end loop over blocks
!$OMP END DO
!$OMP END PARALLEL


END SUBROUTINE recon_lsq_cell_cpoor_svd


!-------------------------------------------------------------------------
!
!
!>
!! Computes coefficients (i.e. derivatives) for cell centered cubic
!! reconstruction.
!!
!! DESCRIPTION:
!! recon: reconstruction of subgrid distribution
!! lsq  : least-squares method
!! cell : solution coefficients defined at cell center
!! c    : cubic reconstruction
!!
!! Computes the coefficients (derivatives) for a cubic reconstruction,
!! using the the least-squares method. The coefficients are provided at
!! cell centers in a local 2D cartesian system (tangential plane).
!! Solves linear system Rx = Q^T d.
!! The matrices have the following size and shape:
!! R  : upper triangular matrix (9 x 9)
!! Q  : orthogonal matrix (9 x 9)
!! Q^T: transposed of Q (9 x 9)
!! d  : input vector (LHS) (9 x 1)
!! x  : solution vector (unknowns) (9 x 1)
!!
!! Coefficients
!! p_coeff(jc,jk,jb, 1) : C0
!! p_coeff(jc,jk,jb, 2) : C1 (dPhi_dx)
!! p_coeff(jc,jk,jb, 3) : C2 (dPhi_dy)
!! p_coeff(jc,jk,jb, 4) : C3 (1/2*ddPhi_ddx)
!! p_coeff(jc,jk,jb, 5) : C4 (1/2*ddPhi_ddy)
!! p_coeff(jc,jk,jb, 6) : C5 (ddPhi_dxdy)
!! p_coeff(jc,jk,jb, 7) : C6 (1/6*dddPhi_dddx)
!! p_coeff(jc,jk,jb, 8) : C7 (1/6*dddPhi_dddy)
!! p_coeff(jc,jk,jb, 9) : C8 (1/2*dddPhi_ddxdy)
!! p_coeff(jc,jk,jb,10) : C9 (1/2*dddPhi_dxddy)
!!
!! works only on triangular grid yet
!!
!! @par Revision History
!!  Developed and tested by Daniel Reinert, DWD (2010-10-15)
!!
!! !LITERATURE
!! Ollivier-Gooch et al (2002): A High-Order-Accurate Unstructured Mesh
!! Finite-Volume Scheme for the Advection-Diffusion Equation, J. Comput. Phys.,
!! 181, 729-752
!!
SUBROUTINE recon_lsq_cell_c( p_cc, ptr_patch, ptr_int_lsq, p_coeff, &
  &                               opt_slev, opt_elev, opt_rlstart,  &
  &                               opt_rlend )

  TYPE(t_patch), INTENT(IN) :: & !< patch on which computation
    &  ptr_patch                 !< is performed
                                        
  TYPE(t_lsq), TARGET, INTENT(IN) :: &  !< data structure for interpolation
    &  ptr_int_lsq

  REAL(wp), INTENT(IN) ::           & !< cell centered variable
    &  p_cc(:,:,:)

  INTEGER, INTENT(IN), OPTIONAL ::  & !< optional vertical start level
    &  opt_slev

  INTEGER, INTENT(IN), OPTIONAL ::  & !< optional vertical end level
    &  opt_elev

  INTEGER, INTENT(IN), OPTIONAL ::  & !< start and end values of refin_ctrl flag
    &  opt_rlstart, opt_rlend

  REAL(wp), INTENT(INOUT)  ::  &  !< cell based coefficients (geographical components)
    &  p_coeff(:,:,:,:)           !< physically this vector contains gradients, second
                                  !< and third derivatives, one mixed derivative and a
                                  !< constant coefficient for zonal and meridional
                                  !< direction

  REAL(wp)  ::           &        !< difference of scalars i j
    &  z_d(lsq_high_set%dim_c,nproma,ptr_patch%nlev)
  REAL(wp)  ::           &        !< matrix-vector product of transposed
    &  z_qt_times_d(9)            !< Q matrix and d

  REAL(wp), POINTER ::   &        !< Pointer to reciprocal diagonal R-matrix-elements
    &  ptr_rrdiag(:,:,:)
  REAL(wp), POINTER ::   &        !< Pointer to upper triangular R-matrix-elements
    &  ptr_rutri(:,:,:)

  INTEGER, POINTER  ::   &        !< Pointer to line and block indices of
    &  iidx(:,:,:), iblk(:,:,:)   !< required stencil
  INTEGER :: slev, elev           !< vertical start and end level
  INTEGER :: jc, jk, jb           !< index of cell, vertical level and block
  INTEGER :: js                   !< loop index for cells in stencil
  INTEGER :: rl_start, rl_end
  INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx, i_nchdom

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
    elev = ptr_patch%nlev
  END IF
  IF ( PRESENT(opt_rlstart) ) THEN
    rl_start = opt_rlstart
  ELSE
    rl_start = 2
  END IF
  IF ( PRESENT(opt_rlend) ) THEN
    rl_end = opt_rlend
  ELSE
    rl_end = min_rlcell
  END IF


  ! values for the blocking
  i_nchdom   = MAX(1,ptr_patch%n_childdom)
  i_startblk = ptr_patch%cells%start_blk(rl_start,1)
  i_endblk   = ptr_patch%cells%end_blk(rl_end,i_nchdom)

  ! pointers to line and block indices of required stencil
  iidx => ptr_int_lsq%lsq_idx_c
  iblk => ptr_int_lsq%lsq_blk_c

  ! pointer to reciprocal diagonal R-elements
  ptr_rrdiag => ptr_int_lsq%lsq_rmat_rdiag_c(:,:,:)

  ! pointer to upper triangular R-elements
  ptr_rutri => ptr_int_lsq%lsq_rmat_utri_c(:,:,:)



!$OMP PARALLEL

  IF (ptr_patch%id > 1) THEN
!$OMP WORKSHARE
    p_coeff(:,:,1:i_startblk,1:10) = 0._wp
!$OMP END WORKSHARE
  ENDIF

!$OMP DO PRIVATE(jb,jc,jk,js,i_startidx,i_endidx,z_d,z_qt_times_d), SCHEDULE(runtime)
  DO jb = i_startblk, i_endblk

    CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, rl_start, rl_end)

    !
    ! 1. compute right hand side of linear system
    !
#ifdef __LOOP_EXCHANGE
    DO jc = i_startidx, i_endidx
      DO jk = slev, elev
#else
!CDIR UNROLL=5
    DO jk = slev, elev
      DO jc = i_startidx, i_endidx
#endif

!CDIR EXPAND=9
        DO js = 1, 9
          ! original version
          !z_d(js) = ptr_int_lsq%lsq_weights_c(jc,js,jb)                   &
          !  &     * (p_cc(iidx(jc,jb,js),jk,iblk(jc,jb,js)) - p_cc(jc,jk,jb))

          ! note, that the multiplication with lsq_weights_c(jc,js,jb) at
          ! runtime is now avoided. Instead, the multiplication with
          ! lsq_weights_c(jc,js,jb) has been shifted into the transposed
          ! Q-matrix.

          z_d(js,jc,jk) = p_cc(iidx(jc,jb,js),jk,iblk(jc,jb,js)) - p_cc(jc,jk,jb)
        ENDDO
      ENDDO
    ENDDO



    !
    ! 2. compute cell based coefficients for quadratic reconstruction
    !
    DO jk = slev, elev

      DO jc = i_startidx, i_endidx

        ! calculate matrix vector product Q^T d (transposed of Q times LHS)
        ! (intrinsic function matmul not applied, due to massive
        ! performance penalty on the NEC. Instead the intrinsic dot product
        ! function is applied
!CDIR BEGIN EXPAND=9
        z_qt_times_d(1) = DOT_PRODUCT(ptr_int_lsq%lsq_qtmat_c(jc,1,1:9,jb),z_d(1:9,jc,jk))
        z_qt_times_d(2) = DOT_PRODUCT(ptr_int_lsq%lsq_qtmat_c(jc,2,1:9,jb),z_d(1:9,jc,jk))
        z_qt_times_d(3) = DOT_PRODUCT(ptr_int_lsq%lsq_qtmat_c(jc,3,1:9,jb),z_d(1:9,jc,jk))
        z_qt_times_d(4) = DOT_PRODUCT(ptr_int_lsq%lsq_qtmat_c(jc,4,1:9,jb),z_d(1:9,jc,jk))
        z_qt_times_d(5) = DOT_PRODUCT(ptr_int_lsq%lsq_qtmat_c(jc,5,1:9,jb),z_d(1:9,jc,jk))
        z_qt_times_d(6) = DOT_PRODUCT(ptr_int_lsq%lsq_qtmat_c(jc,6,1:9,jb),z_d(1:9,jc,jk))
        z_qt_times_d(7) = DOT_PRODUCT(ptr_int_lsq%lsq_qtmat_c(jc,7,1:9,jb),z_d(1:9,jc,jk))
        z_qt_times_d(8) = DOT_PRODUCT(ptr_int_lsq%lsq_qtmat_c(jc,8,1:9,jb),z_d(1:9,jc,jk))
        z_qt_times_d(9) = DOT_PRODUCT(ptr_int_lsq%lsq_qtmat_c(jc,9,1:9,jb),z_d(1:9,jc,jk))
!CDIR END


        !
        ! Solve linear system Rx=Q^T d by back substitution
        !
        p_coeff(jc,jk,jb,10) = ptr_rrdiag(jc,9,jb) * z_qt_times_d(9)
        p_coeff(jc,jk,jb,9)  = ptr_rrdiag(jc,8,jb)                                         &
          &                  * ( z_qt_times_d(8) - ptr_rutri(jc,1,jb)*p_coeff(jc,jk,jb,10) )
        p_coeff(jc,jk,jb,8)  = ptr_rrdiag(jc,7,jb)                                         &
          &                  * ( z_qt_times_d(7) - ptr_rutri(jc,2,jb)*p_coeff(jc,jk,jb,9)  &
          &                  - ptr_rutri(jc,3,jb) * p_coeff(jc,jk,jb,10) )
        p_coeff(jc,jk,jb,7)  = ptr_rrdiag(jc,6,jb)                                         &
          &                  * ( z_qt_times_d(6) - ptr_rutri(jc,4,jb)*p_coeff(jc,jk,jb,8)  &
          &                  - ptr_rutri(jc,5,jb) * p_coeff(jc,jk,jb,9)                    &
          &                  - ptr_rutri(jc,6,jb) * p_coeff(jc,jk,jb,10) )
        p_coeff(jc,jk,jb,6)  = ptr_rrdiag(jc,5,jb)                                         &
          &                  * ( z_qt_times_d(5) - ptr_rutri(jc,7,jb)*p_coeff(jc,jk,jb,7)  &
          &                  - ptr_rutri(jc,8,jb) * p_coeff(jc,jk,jb,8)                    &
          &                  - ptr_rutri(jc,9,jb) * p_coeff(jc,jk,jb,9)                    &
          &                  - ptr_rutri(jc,10,jb)* p_coeff(jc,jk,jb,10) )
        p_coeff(jc,jk,jb,5)  = ptr_rrdiag(jc,4,jb)                                         &
          &                  * ( z_qt_times_d(4) - ptr_rutri(jc,11,jb)*p_coeff(jc,jk,jb,6) &
          &                  - ptr_rutri(jc,12,jb) * p_coeff(jc,jk,jb,7)                   &
          &                  - ptr_rutri(jc,13,jb) * p_coeff(jc,jk,jb,8)                   &
          &                  - ptr_rutri(jc,14,jb) * p_coeff(jc,jk,jb,9)                   &
          &                  - ptr_rutri(jc,15,jb) * p_coeff(jc,jk,jb,10) )
        p_coeff(jc,jk,jb,4)  = ptr_rrdiag(jc,3,jb)                                         &
          &                  * ( z_qt_times_d(3) - ptr_rutri(jc,16,jb)*p_coeff(jc,jk,jb,5) &
          &                  - ptr_rutri(jc,17,jb) * p_coeff(jc,jk,jb,6)                   &
          &                  - ptr_rutri(jc,18,jb) * p_coeff(jc,jk,jb,7)                   &
          &                  - ptr_rutri(jc,19,jb) * p_coeff(jc,jk,jb,8)                   &
          &                  - ptr_rutri(jc,20,jb) * p_coeff(jc,jk,jb,9)                   &
          &                  - ptr_rutri(jc,21,jb) * p_coeff(jc,jk,jb,10) )
        p_coeff(jc,jk,jb,3)  = ptr_rrdiag(jc,2,jb)                                         &
          &                  * ( z_qt_times_d(2) - ptr_rutri(jc,22,jb)*p_coeff(jc,jk,jb,4) &
          &                  - ptr_rutri(jc,23,jb) * p_coeff(jc,jk,jb,5)                   &
          &                  - ptr_rutri(jc,24,jb) * p_coeff(jc,jk,jb,6)                   &
          &                  - ptr_rutri(jc,25,jb) * p_coeff(jc,jk,jb,7)                   &
          &                  - ptr_rutri(jc,26,jb) * p_coeff(jc,jk,jb,8)                   &
          &                  - ptr_rutri(jc,27,jb) * p_coeff(jc,jk,jb,9)                   &
          &                  - ptr_rutri(jc,28,jb) * p_coeff(jc,jk,jb,10) )
        p_coeff(jc,jk,jb,2)  = ptr_rrdiag(jc,1,jb)                                         &
          &                  * ( z_qt_times_d(1) - ptr_rutri(jc,29,jb)*p_coeff(jc,jk,jb,3) &
          &                  - ptr_rutri(jc,30,jb) * p_coeff(jc,jk,jb,4)                   &
          &                  - ptr_rutri(jc,31,jb) * p_coeff(jc,jk,jb,5)                   &
          &                  - ptr_rutri(jc,32,jb) * p_coeff(jc,jk,jb,6)                   &
          &                  - ptr_rutri(jc,33,jb) * p_coeff(jc,jk,jb,7)                   &
          &                  - ptr_rutri(jc,34,jb) * p_coeff(jc,jk,jb,8)                   &
          &                  - ptr_rutri(jc,35,jb) * p_coeff(jc,jk,jb,9)                   &
          &                  - ptr_rutri(jc,36,jb) * p_coeff(jc,jk,jb,10) )


        p_coeff(jc,jk,jb,1)  = p_cc(jc,jk,jb)                                              &
          &                  - p_coeff(jc,jk,jb,2)  * ptr_int_lsq%lsq_moments(jc,jb,1)     &
          &                  - p_coeff(jc,jk,jb,3)  * ptr_int_lsq%lsq_moments(jc,jb,2)     &
          &                  - p_coeff(jc,jk,jb,4)  * ptr_int_lsq%lsq_moments(jc,jb,3)     &
          &                  - p_coeff(jc,jk,jb,5)  * ptr_int_lsq%lsq_moments(jc,jb,4)     &
          &                  - p_coeff(jc,jk,jb,6)  * ptr_int_lsq%lsq_moments(jc,jb,5)     &
          &                  - p_coeff(jc,jk,jb,7)  * ptr_int_lsq%lsq_moments(jc,jb,6)     &
          &                  - p_coeff(jc,jk,jb,8)  * ptr_int_lsq%lsq_moments(jc,jb,7)     &
          &                  - p_coeff(jc,jk,jb,9)  * ptr_int_lsq%lsq_moments(jc,jb,8)     &
          &                  - p_coeff(jc,jk,jb,10) * ptr_int_lsq%lsq_moments(jc,jb,9)

      END DO ! end loop over cells

    END DO ! end loop over vertical levels

  END DO ! end loop over blocks
!$OMP END DO
!$OMP END PARALLEL


END SUBROUTINE recon_lsq_cell_c




!--------------------------------------------------
!
!
!>
!! Computes coefficients (i.e. derivatives) for cell centered cubic
!! reconstruction.
!!
!! DESCRIPTION:
!! recon: reconstruction of subgrid distribution
!! lsq  : least-squares method
!! cell : solution coefficients defined at cell center
!! c    : cubic reconstruction
!!
!! Computes unknown coefficients (derivatives) of a cubic polynomial,
!! using the least-squares method. The coefficients are provided at 
!! cell centers in a local 2D cartesian system (tangential plane).
!!
!! Mathematically we solve Ax = b via Singular Value Decomposition (SVD)
!! x = PINV(A) * b
!!
!! Matrices have the following size and shape (triangular grid) :
!! PINV(A): Pseudo or Moore-Penrose inverse of A (via SVD) (9 x 9)
!! b  : input vector (LHS) (9 x 1)
!! x  : solution vector (unknowns) (9 x 1)
!!
!! Coefficients
!! p_coeff(jc,jk,jb, 1) : C0
!! p_coeff(jc,jk,jb, 2) : C1 (dPhi_dx)
!! p_coeff(jc,jk,jb, 3) : C2 (dPhi_dy)
!! p_coeff(jc,jk,jb, 4) : C3 (1/2*ddPhi_ddx)
!! p_coeff(jc,jk,jb, 5) : C4 (1/2*ddPhi_ddy)
!! p_coeff(jc,jk,jb, 6) : C5 (ddPhi_dxdy)
!! p_coeff(jc,jk,jb, 7) : C6 (1/6*dddPhi_dddx)
!! p_coeff(jc,jk,jb, 8) : C7 (1/6*dddPhi_dddy)
!! p_coeff(jc,jk,jb, 9) : C8 (1/2*dddPhi_ddxdy)
!! p_coeff(jc,jk,jb,10) : C9 (1/2*dddPhi_dxddy)
!!
!! works only on triangular grid yet
!!
!! @par Revision History
!!  Developed and tested by Daniel Reinert, DWD (2011-05-31)
!!
!! !LITERATURE
!! Ollivier-Gooch et al (2002): A High-Order-Accurate Unstructured Mesh
!! Finite-Volume Scheme for the Advection-Diffusion Equation, J. Comput. Phys.,
!! 181, 729-752
!!
SUBROUTINE recon_lsq_cell_c_svd( p_cc, ptr_patch, ptr_int_lsq, p_coeff, &
  &                               opt_slev, opt_elev, opt_rlstart,      &
  &                               opt_rlend )

  TYPE(t_patch), INTENT(IN) :: & !< patch on which computation
    &  ptr_patch                 !< is performed
                                        
  TYPE(t_lsq), TARGET, INTENT(IN) :: &  !< data structure for interpolation
    &  ptr_int_lsq

  REAL(wp), INTENT(IN) ::           & !< cell centered variable
    &  p_cc(:,:,:)

  INTEGER, INTENT(IN), OPTIONAL ::  & !< optional vertical start level
    &  opt_slev

  INTEGER, INTENT(IN), OPTIONAL ::  & !< optional vertical end level
    &  opt_elev

  INTEGER, INTENT(IN), OPTIONAL ::  & !< start and end values of refin_ctrl flag
    &  opt_rlstart, opt_rlend

  REAL(wp), INTENT(INOUT)  ::  &  !< cell based coefficients (geographical components)
    &  p_coeff(:,:,:,:)           !< physically this vector contains gradients, second
                                  !< and third derivatives, one mixed derivative and a
                                  !< constant coefficient for zonal and meridional
                                  !< direction

  REAL(wp)  ::           &        !< difference of scalars i j
    &  z_b(lsq_high_set%dim_c,nproma,ptr_patch%nlev)

  INTEGER, POINTER  ::   &        !< Pointer to line and block indices of
    &  iidx(:,:,:), iblk(:,:,:)   !< required stencil
  INTEGER :: slev, elev           !< vertical start and end level
  INTEGER :: jc, jk, jb           !< index of cell, vertical level and block
  INTEGER :: js                   !< loop index for cells in stencil
  INTEGER :: rl_start, rl_end
  INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx, i_nchdom

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
    elev = ptr_patch%nlev
  END IF
  IF ( PRESENT(opt_rlstart) ) THEN
    rl_start = opt_rlstart
  ELSE
    rl_start = 2
  END IF
  IF ( PRESENT(opt_rlend) ) THEN
    rl_end = opt_rlend
  ELSE
    rl_end = min_rlcell
  END IF


  ! values for the blocking
  i_nchdom   = MAX(1,ptr_patch%n_childdom)
  i_startblk = ptr_patch%cells%start_blk(rl_start,1)
  i_endblk   = ptr_patch%cells%end_blk(rl_end,i_nchdom)

  ! pointers to line and block indices of required stencil
  iidx => ptr_int_lsq%lsq_idx_c
  iblk => ptr_int_lsq%lsq_blk_c



!$OMP PARALLEL

  IF (ptr_patch%id > 1) THEN
!$OMP WORKSHARE
    p_coeff(:,:,1:i_startblk,1:10) = 0._wp
!$OMP END WORKSHARE
  ENDIF

!$OMP DO PRIVATE(jb,jc,jk,js,i_startidx,i_endidx,z_b), SCHEDULE(runtime)
  DO jb = i_startblk, i_endblk

    CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, rl_start, rl_end)

    !
    ! 1. compute right hand side of linear system
    !
#ifdef __LOOP_EXCHANGE
    DO jc = i_startidx, i_endidx
      DO jk = slev, elev
#else
!CDIR UNROLL=5
    DO jk = slev, elev
      DO jc = i_startidx, i_endidx
#endif

!CDIR EXPAND=9
        DO js = 1, 9
          ! original version
          !z_d(js) = ptr_int_lsq%lsq_weights_c(jc,js,jb)                   &
          !  &     * (p_cc(iidx(jc,jb,js),jk,iblk(jc,jb,js)) - p_cc(jc,jk,jb))

          ! note that multiplication with lsq_weights_c(jc,js,jb) during 
          ! runtime is now avoided. Instead, the weights have been shifted 
          ! into the pseudoinverse.

          z_b(js,jc,jk) = p_cc(iidx(jc,jb,js),jk,iblk(jc,jb,js)) - p_cc(jc,jk,jb)
        ENDDO
      ENDDO
    ENDDO



    !
    ! 2. compute cell based coefficients for cubic reconstruction
    !    calculate matrix vector product PINV(A) * b
    !
    DO jk = slev, elev

      DO jc = i_startidx, i_endidx

        ! (intrinsic function matmul not applied, due to massive
        ! performance penalty on the NEC. Instead the intrinsic dot product
        ! function is applied
!CDIR BEGIN EXPAND=9

        p_coeff(jc,jk,jb,10) = DOT_PRODUCT(ptr_int_lsq%lsq_pseudoinv(jc,9,1:9,jb), &
          &                               z_b(1:9,jc,jk))
        p_coeff(jc,jk,jb, 9) = DOT_PRODUCT(ptr_int_lsq%lsq_pseudoinv(jc,8,1:9,jb), &
          &                               z_b(1:9,jc,jk))
        p_coeff(jc,jk,jb, 8) = DOT_PRODUCT(ptr_int_lsq%lsq_pseudoinv(jc,7,1:9,jb), &
          &                               z_b(1:9,jc,jk))
        p_coeff(jc,jk,jb, 7) = DOT_PRODUCT(ptr_int_lsq%lsq_pseudoinv(jc,6,1:9,jb), &
          &                               z_b(1:9,jc,jk))
        p_coeff(jc,jk,jb, 6) = DOT_PRODUCT(ptr_int_lsq%lsq_pseudoinv(jc,5,1:9,jb), &
          &                               z_b(1:9,jc,jk))
        p_coeff(jc,jk,jb, 5) = DOT_PRODUCT(ptr_int_lsq%lsq_pseudoinv(jc,4,1:9,jb), &
          &                               z_b(1:9,jc,jk))
        p_coeff(jc,jk,jb, 4) = DOT_PRODUCT(ptr_int_lsq%lsq_pseudoinv(jc,3,1:9,jb), &
          &                               z_b(1:9,jc,jk))
        p_coeff(jc,jk,jb, 3) = DOT_PRODUCT(ptr_int_lsq%lsq_pseudoinv(jc,2,1:9,jb), &
          &                               z_b(1:9,jc,jk))
        p_coeff(jc,jk,jb, 2) = DOT_PRODUCT(ptr_int_lsq%lsq_pseudoinv(jc,1,1:9,jb), &
          &                               z_b(1:9,jc,jk))

!CDIR END


!CDIR BEGIN EXPAND=9
        p_coeff(jc,jk,jb,1)  = p_cc(jc,jk,jb) - DOT_PRODUCT(p_coeff(jc,jk,jb,2:10), &
          &                    ptr_int_lsq%lsq_moments(jc,jb,1:9))


      END DO ! end loop over cells

    END DO ! end loop over vertical levels

  END DO ! end loop over blocks
!$OMP END DO
!$OMP END PARALLEL


END SUBROUTINE recon_lsq_cell_c_svd

!-------------------------------------------------------------------------
!
!
!>
!! Computes discrete divergence of a vector field.
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
!! -Switch fom array arguments to pointers
!! Modification by Almut Gassmann, MPI-M (2007-04-20)
!! - abandon grid for the sake of patch
!! Modification by Guenther Zaengl, DWD (2009-03-17)
!! - vector optimization
!! Modification by Guenther Zaengl, DWD (2010-08-20) :
!! - Option for processing two fields at once for efficiency optimization
!!
SUBROUTINE div3d( vec_e, ptr_patch, ptr_int, div_vec_c, &
  &              opt_slev, opt_elev, opt_in2, opt_out2, opt_rlstart, opt_rlend )
!

!
!  patch on which computation is performed
!
TYPE(t_patch), TARGET, INTENT(in) :: ptr_patch

! Interpolation state
TYPE(t_int_state), INTENT(in)     :: ptr_int
!
! edge based variable of which divergence
! is computed
!
REAL(wp), INTENT(in) ::  &
  &  vec_e(:,:,:) ! dim: (nproma,nlev,nblks_e)

! optional second input field for more efficient processing in NH core
REAL(wp), OPTIONAL, INTENT(in) ::  &
  &  opt_in2(:,:,:) ! dim: (nproma,nlev,nblks_e)

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_slev    ! optional vertical start level

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_elev    ! optional vertical end level

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_rlstart, opt_rlend   ! start and end values of refin_ctrl flag

!
! cell based variable in which divergence is stored
!
!REAL(wp), INTENT(out) ::  &
REAL(wp), INTENT(inout) ::  &
  &  div_vec_c(:,:,:) ! dim: (nproma,nlev,nblks_c)

! optional second output field
REAL(wp), OPTIONAL, INTENT(inout) ::  &
  &  opt_out2(:,:,:) ! dim: (nproma,nlev,nblks_c)

INTEGER :: slev, elev     ! vertical start and end level
INTEGER :: jc, jk, jb
INTEGER :: rl_start, rl_end
INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx, i_nchdom
INTEGER :: nlen, npromz_c, nblks_c

INTEGER,  DIMENSION(:,:,:),   POINTER :: iidx, iblk
LOGICAL :: l2fields

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
  elev = ptr_patch%nlev
END IF

IF ( PRESENT(opt_rlstart) ) THEN
  rl_start = opt_rlstart
ELSE
  rl_start = 1
END IF
IF ( PRESENT(opt_rlend) ) THEN
  rl_end = opt_rlend
ELSE
  rl_end = min_rlcell
END IF

IF ( PRESENT(opt_in2) .AND. PRESENT(opt_out2)) THEN
  l2fields = .TRUE.
ELSE
  l2fields = .FALSE.
ENDIF

iidx => ptr_patch%cells%edge_idx
iblk => ptr_patch%cells%edge_blk

! values for the blocking
i_nchdom   = MAX(1,ptr_patch%n_childdom)
i_startblk = ptr_patch%cells%start_blk(rl_start,1)
i_endblk   = ptr_patch%cells%end_blk(rl_end,i_nchdom)

! loop through all patch cells (and blocks)
!

!IF(ltimer) CALL timer_start(timer_div)

!$OMP PARALLEL


SELECT CASE (ptr_patch%cell_type)

CASE (3) ! (cell_type == 3)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk)
  DO jb = i_startblk, i_endblk

    CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, rl_start, rl_end)

    ! original comment for divergence computation;
    ! everything that follows in this explanation has been combined into geofac_div

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

    IF (l2fields) THEN

#ifdef __LOOP_EXCHANGE
      DO jc = i_startidx, i_endidx
        DO jk = slev, elev
#else
!CDIR UNROLL=5
      DO jk = slev, elev
        DO jc = i_startidx, i_endidx
#endif

          div_vec_c(jc,jk,jb) =  &
            vec_e(iidx(jc,jb,1),jk,iblk(jc,jb,1)) * ptr_int%geofac_div(jc,1,jb) + &
            vec_e(iidx(jc,jb,2),jk,iblk(jc,jb,2)) * ptr_int%geofac_div(jc,2,jb) + &
            vec_e(iidx(jc,jb,3),jk,iblk(jc,jb,3)) * ptr_int%geofac_div(jc,3,jb)

          opt_out2(jc,jk,jb) =  &
            opt_in2(iidx(jc,jb,1),jk,iblk(jc,jb,1)) * ptr_int%geofac_div(jc,1,jb) + &
            opt_in2(iidx(jc,jb,2),jk,iblk(jc,jb,2)) * ptr_int%geofac_div(jc,2,jb) + &
            opt_in2(iidx(jc,jb,3),jk,iblk(jc,jb,3)) * ptr_int%geofac_div(jc,3,jb)

        END DO
      END DO

    ELSE

#ifdef __LOOP_EXCHANGE
      DO jc = i_startidx, i_endidx
        DO jk = slev, elev
#else
!CDIR UNROLL=6
      DO jk = slev, elev
        DO jc = i_startidx, i_endidx
#endif

          div_vec_c(jc,jk,jb) =  &
            vec_e(iidx(jc,jb,1),jk,iblk(jc,jb,1)) * ptr_int%geofac_div(jc,1,jb) + &
            vec_e(iidx(jc,jb,2),jk,iblk(jc,jb,2)) * ptr_int%geofac_div(jc,2,jb) + &
            vec_e(iidx(jc,jb,3),jk,iblk(jc,jb,3)) * ptr_int%geofac_div(jc,3,jb)

        END DO
      END DO

    ENDIF
  END DO
!$OMP END DO

CASE (6) ! (cell_type == 6)

  ! no grid refinement in hexagonal model
  nblks_c   = ptr_patch%nblks_int_c
  npromz_c  = ptr_patch%npromz_int_c

!$OMP DO PRIVATE(jb,nlen,jc,jk)
  DO jb = 1, nblks_c

    IF (jb /= nblks_c) THEN
      nlen = nproma
    ELSE
      nlen = npromz_c
    ENDIF

#ifdef __LOOP_EXCHANGE
    DO jc = 1, nlen
      DO jk = slev, elev
#else
!CDIR UNROLL=6
    DO jk = slev, elev
      DO jc = 1, nlen
#endif

          div_vec_c(jc,jk,jb) =   &
          &   vec_e(iidx(jc,jb,1),jk,iblk(jc,jb,1)) * ptr_int%geofac_div(jc,1,jb) &
          & + vec_e(iidx(jc,jb,2),jk,iblk(jc,jb,2)) * ptr_int%geofac_div(jc,2,jb) &
          & + vec_e(iidx(jc,jb,3),jk,iblk(jc,jb,3)) * ptr_int%geofac_div(jc,3,jb) &
          & + vec_e(iidx(jc,jb,4),jk,iblk(jc,jb,4)) * ptr_int%geofac_div(jc,4,jb) &
          & + vec_e(iidx(jc,jb,5),jk,iblk(jc,jb,5)) * ptr_int%geofac_div(jc,5,jb) &
          & + vec_e(iidx(jc,jb,6),jk,iblk(jc,jb,6)) * ptr_int%geofac_div(jc,6,jb)

      END DO
    END DO
  END DO
!$OMP END DO

END SELECT
!$OMP END PARALLEL

!IF(ltimer) CALL timer_stop(timer_div)

END SUBROUTINE div3d

!-------------------------------------------------------------------------
!
!
!>
!! Special version of div that processes 4D fields in one step
!!
!! See standard routine (div3d) for further description
!!
!! @par Revision History
!! Developed  by Guenther Zaengl, DWD (2011-01-12)
!!
SUBROUTINE div4d( ptr_patch, ptr_int, f4din, f4dout, dim4d, &
  &              opt_slev, opt_elev, opt_rlstart, opt_rlend )
!

!
!  patch on which computation is performed
!
TYPE(t_patch), TARGET, INTENT(in) :: ptr_patch

! Interpolation state
TYPE(t_int_state), INTENT(in)     :: ptr_int
!
! edge based 4D input field of which divergence is computed
!
REAL(wp), INTENT(in) ::  &
  &  f4din(:,:,:,:) ! dim: (nproma,nlev,nblks_e,dim4d)

INTEGER, INTENT(in) :: dim4d ! Last dimension of the input/output fields

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_slev(dim4d)    ! optional vertical start level

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_elev(dim4d)    ! optional vertical end level

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_rlstart, opt_rlend   ! start and end values of refin_ctrl flag

!
! 4D cell based variable in which divergence is stored
!
REAL(wp), INTENT(inout) ::  &
  &  f4dout(:,:,:,:) ! dim: (nproma,nlev,nblks_c,dim4d)

INTEGER :: slev(dim4d), elev(dim4d)     ! vertical start and end level
INTEGER :: jc, jk, jb, ji
INTEGER :: rl_start, rl_end
INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx, i_nchdom
INTEGER :: nlen, npromz_c, nblks_c

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
  elev = ptr_patch%nlev
END IF

IF ( PRESENT(opt_rlstart) ) THEN
  rl_start = opt_rlstart
ELSE
  rl_start = 1
END IF
IF ( PRESENT(opt_rlend) ) THEN
  rl_end = opt_rlend
ELSE
  rl_end = min_rlcell
END IF

iidx => ptr_patch%cells%edge_idx
iblk => ptr_patch%cells%edge_blk

! values for the blocking
i_nchdom   = MAX(1,ptr_patch%n_childdom)
i_startblk = ptr_patch%cells%start_blk(rl_start,1)
i_endblk   = ptr_patch%cells%end_blk(rl_end,i_nchdom)

! loop through all patch cells (and blocks)
!

!IF(ltimer) CALL timer_start(timer_div)

!$OMP PARALLEL


SELECT CASE (ptr_patch%cell_type)

CASE (3) ! (cell_type == 3)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk,ji)
  DO jb = i_startblk, i_endblk

    CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, rl_start, rl_end)

#ifdef __LOOP_EXCHANGE
    DO jc = i_startidx, i_endidx
      DO ji = 1, dim4d
        DO jk = slev(ji), elev(ji)
#else
    DO ji = 1, dim4d
!CDIR UNROLL=6
      DO jk = slev(ji), elev(ji)
        DO jc = i_startidx, i_endidx
#endif

          f4dout(jc,jk,jb,ji) =  &
            f4din(iidx(jc,jb,1),jk,iblk(jc,jb,1),ji) * ptr_int%geofac_div(jc,1,jb) + &
            f4din(iidx(jc,jb,2),jk,iblk(jc,jb,2),ji) * ptr_int%geofac_div(jc,2,jb) + &
            f4din(iidx(jc,jb,3),jk,iblk(jc,jb,3),ji) * ptr_int%geofac_div(jc,3,jb)

        ENDDO
      ENDDO
    ENDDO

  ENDDO
!$OMP END DO

CASE (6) ! (cell_type == 6)

  ! no grid refinement in hexagonal model
  nblks_c   = ptr_patch%nblks_int_c
  npromz_c  = ptr_patch%npromz_int_c

!$OMP DO PRIVATE(jb,nlen,jc,jk,ji)
  DO jb = 1, nblks_c

    IF (jb /= nblks_c) THEN
      nlen = nproma
    ELSE
      nlen = npromz_c
    ENDIF

#ifdef __LOOP_EXCHANGE
    DO jc = 1, nlen
      DO ji = 1, dim4d
        DO jk = slev(ji), elev(ji)
#else
    DO ji = 1, dim4d
!CDIR UNROLL=6
      DO jk = slev(ji), elev(ji)
        DO jc = 1, nlen
#endif

          f4dout(jc,jk,jb,ji) =   &
          &   f4din(iidx(jc,jb,1),jk,iblk(jc,jb,1),ji) * ptr_int%geofac_div(jc,1,jb) &
          & + f4din(iidx(jc,jb,2),jk,iblk(jc,jb,2),ji) * ptr_int%geofac_div(jc,2,jb) &
          & + f4din(iidx(jc,jb,3),jk,iblk(jc,jb,3),ji) * ptr_int%geofac_div(jc,3,jb) &
          & + f4din(iidx(jc,jb,4),jk,iblk(jc,jb,4),ji) * ptr_int%geofac_div(jc,4,jb) &
          & + f4din(iidx(jc,jb,5),jk,iblk(jc,jb,5),ji) * ptr_int%geofac_div(jc,5,jb) &
          & + f4din(iidx(jc,jb,6),jk,iblk(jc,jb,6),ji) * ptr_int%geofac_div(jc,6,jb)

        ENDDO
      ENDDO
    ENDDO

  ENDDO
!$OMP END DO

END SELECT
!$OMP END PARALLEL

!IF(ltimer) CALL timer_stop(timer_div)

END SUBROUTINE div4d


!-------------------------------------------------------------------------
!
!
!>
!! Computes discrete divergence of a vector field.
!!
!! Computes discrete divergence of a vector field
!! given by its components in the directions normal to triangle edges,
!! followed by bilinear averaging to remove checkerboard noise
!! (Combines div_midpoint and cell_avg_varwgt to increase computing efficiency)
!!
!! @par Revision History
!! Developed by Guenther Zaengl, DWD (2009-03-30)
!! Modification by Guenther Zaengl, DWD (2010-04-20) :
!! - Option for processing two fields at once for efficiency optimization
!!
SUBROUTINE div_avg( vec_e, ptr_patch, ptr_int, avg_coeff, div_vec_c, &
  &                 opt_in2, opt_out2, opt_rlstart, opt_rlend )
!
!
!  patch on which computation is performed
!
TYPE(t_patch), TARGET, INTENT(in) :: ptr_patch

! Interpolation state
TYPE(t_int_state), INTENT(in)     :: ptr_int

!  averaging coefficients
REAL(wp), INTENT(in) :: avg_coeff(:,:,:) ! dim: (nproma,nlev,nblks_c)
!
! edge based variable of which divergence
! is computed
!
REAL(wp), INTENT(in) ::  &
  &  vec_e(:,:,:) ! dim: (nproma,nlev,nblks_e)

! optional second input field for more efficient processing in NH core
REAL(wp), OPTIONAL, INTENT(in) ::  &
  &  opt_in2(:,:,:) ! dim: (nproma,nlev,nblks_e)

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_rlstart, opt_rlend   ! start and end values of refin_ctrl flag

!
! cell based variable in which divergence is stored
!
REAL(wp), INTENT(inout) ::  &
  &  div_vec_c(:,:,:) ! dim: (nproma,nlev,nblks_c)

! optional second output field
REAL(wp), OPTIONAL, INTENT(inout) ::  &
  &  opt_out2(:,:,:) ! dim: (nproma,nlev,nblks_c)

INTEGER :: jc, jk, jb
INTEGER :: rl_start, rl_end, rl_start_l2, rl_end_l1
INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx, i_nchdom
INTEGER :: nlev              !< number of full levels

REAL(wp), DIMENSION (nproma,ptr_patch%nlev,ptr_patch%nblks_c) :: aux_c, aux_c2

INTEGER,  DIMENSION(:,:,:),   POINTER :: inidx, inblk, ieidx, ieblk
LOGICAL :: l2fields

!-----------------------------------------------------------------------

! check optional arguments
IF ( PRESENT(opt_rlstart) ) THEN
  rl_start = opt_rlstart
ELSE
  rl_start = 1
END IF
IF ( PRESENT(opt_rlend) ) THEN
  rl_end = opt_rlend
ELSE
  rl_end = min_rlcell
END IF
IF ( PRESENT(opt_in2) .AND. PRESENT(opt_out2)) THEN
  l2fields = .TRUE.
ELSE
  l2fields = .FALSE.
ENDIF

rl_start_l2 = rl_start + 1

IF ( PRESENT(opt_rlend) .AND. rl_end < 0 .AND. rl_end > min_rlcell ) THEN
  rl_end_l1 = rl_end - 1
ELSE
  rl_end_l1 = rl_end
END IF

inidx => ptr_patch%cells%neighbor_idx
inblk => ptr_patch%cells%neighbor_blk
ieidx => ptr_patch%cells%edge_idx
ieblk => ptr_patch%cells%edge_blk

! number of vertical levels
nlev = ptr_patch%nlev

! values for the blocking
i_nchdom   = MAX(1,ptr_patch%n_childdom)

! First compute divergence
!
!$OMP PARALLEL PRIVATE(i_startblk,i_endblk)
i_startblk = ptr_patch%cells%start_blk(rl_start,1)
i_endblk   = ptr_patch%cells%end_blk(rl_end_l1,i_nchdom)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk), SCHEDULE(runtime)
DO jb = i_startblk, i_endblk

  CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
                     i_startidx, i_endidx, rl_start, rl_end_l1)

  IF (l2fields) THEN

#ifdef __LOOP_EXCHANGE
    DO jc = i_startidx, i_endidx
      DO jk = 1, nlev
#else
!CDIR UNROLL=5
    DO jk = 1, nlev
      DO jc = i_startidx, i_endidx
#endif

        aux_c(jc,jk,jb) =  &
          vec_e(ieidx(jc,jb,1),jk,ieblk(jc,jb,1)) * ptr_int%geofac_div(jc,1,jb) + &
          vec_e(ieidx(jc,jb,2),jk,ieblk(jc,jb,2)) * ptr_int%geofac_div(jc,2,jb) + &
          vec_e(ieidx(jc,jb,3),jk,ieblk(jc,jb,3)) * ptr_int%geofac_div(jc,3,jb)

        aux_c2(jc,jk,jb) =  &
          opt_in2(ieidx(jc,jb,1),jk,ieblk(jc,jb,1)) * ptr_int%geofac_div(jc,1,jb) + &
          opt_in2(ieidx(jc,jb,2),jk,ieblk(jc,jb,2)) * ptr_int%geofac_div(jc,2,jb) + &
          opt_in2(ieidx(jc,jb,3),jk,ieblk(jc,jb,3)) * ptr_int%geofac_div(jc,3,jb)

      END DO
    END DO
  ELSE

#ifdef __LOOP_EXCHANGE
    DO jc = i_startidx, i_endidx
      DO jk = 1, nlev
#else
!CDIR UNROLL=6
    DO jk = 1, nlev
      DO jc = i_startidx, i_endidx
#endif
        aux_c(jc,jk,jb) =  &
          vec_e(ieidx(jc,jb,1),jk,ieblk(jc,jb,1)) * ptr_int%geofac_div(jc,1,jb) + &
          vec_e(ieidx(jc,jb,2),jk,ieblk(jc,jb,2)) * ptr_int%geofac_div(jc,2,jb) + &
          vec_e(ieidx(jc,jb,3),jk,ieblk(jc,jb,3)) * ptr_int%geofac_div(jc,3,jb)

      END DO
    END DO
  ENDIF

END DO
!$OMP END DO

IF (l_limited_area .OR. ptr_patch%id > 1) THEN
  ! Fill div_vec_c along the lateral boundaries of nests

  i_startblk = ptr_patch%cells%start_blk(rl_start,1)
  i_endblk   = ptr_patch%cells%end_blk(rl_start_l2,1)
!
  IF (l2fields) THEN
!$OMP WORKSHARE
     div_vec_c(:,:,i_startblk:i_endblk) =  aux_c (:,:,i_startblk:i_endblk)
     opt_out2 (:,:,i_startblk:i_endblk) =  aux_c2(:,:,i_startblk:i_endblk)
!$OMP END WORKSHARE
  ELSE
!$OMP WORKSHARE
     div_vec_c(:,:,i_startblk:i_endblk) =  aux_c(:,:,i_startblk:i_endblk)
!$OMP END WORKSHARE
  ENDIF
ENDIF

!
! Now do averaging with weights given by avg_coeff

! values for the blocking
i_startblk = ptr_patch%cells%start_blk(rl_start_l2,1)
i_endblk   = ptr_patch%cells%end_blk(rl_end,i_nchdom)
!
! loop through all patch cells (and blocks)
!

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk), SCHEDULE(runtime)
DO jb = i_startblk, i_endblk

  CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
                     i_startidx, i_endidx, rl_start_l2, rl_end)

  IF (l2fields) THEN

#ifdef __LOOP_EXCHANGE
    DO jc = i_startidx, i_endidx
      DO jk = 1, nlev
#else
#ifdef _URD
!CDIR UNROLL=_URD
#endif
    DO jk = 1, nlev
      DO jc = i_startidx, i_endidx
#endif

        !  calculate the weighted average
        div_vec_c(jc,jk,jb) =  &
          &    aux_c(jc,jk,jb)                         * avg_coeff(jc,1,jb) &
          &  + aux_c(inidx(jc,jb,1),jk,inblk(jc,jb,1)) * avg_coeff(jc,2,jb) &
          &  + aux_c(inidx(jc,jb,2),jk,inblk(jc,jb,2)) * avg_coeff(jc,3,jb) &
          &  + aux_c(inidx(jc,jb,3),jk,inblk(jc,jb,3)) * avg_coeff(jc,4,jb)

        opt_out2(jc,jk,jb) =  &
          &    aux_c2(jc,jk,jb)                         * avg_coeff(jc,1,jb) &
          &  + aux_c2(inidx(jc,jb,1),jk,inblk(jc,jb,1)) * avg_coeff(jc,2,jb) &
          &  + aux_c2(inidx(jc,jb,2),jk,inblk(jc,jb,2)) * avg_coeff(jc,3,jb) &
          &  + aux_c2(inidx(jc,jb,3),jk,inblk(jc,jb,3)) * avg_coeff(jc,4,jb)

      END DO !cell loop
    END DO !vertical levels loop
  ELSE

#ifdef __LOOP_EXCHANGE
    DO jc = i_startidx, i_endidx
      DO jk = 1, nlev
#else
#ifdef _URD
!CDIR UNROLL=_URD
#endif
    DO jk = 1, nlev
      DO jc = i_startidx, i_endidx
#endif

        !  calculate the weighted average
        div_vec_c(jc,jk,jb) =  &
          &    aux_c(jc,jk,jb)                         * avg_coeff(jc,1,jb) &
          &  + aux_c(inidx(jc,jb,1),jk,inblk(jc,jb,1)) * avg_coeff(jc,2,jb) &
          &  + aux_c(inidx(jc,jb,2),jk,inblk(jc,jb,2)) * avg_coeff(jc,3,jb) &
          &  + aux_c(inidx(jc,jb,3),jk,inblk(jc,jb,3)) * avg_coeff(jc,4,jb)

      END DO !cell loop
    END DO !vertical levels loop
  ENDIF

END DO !block loop
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE div_avg

!-------------------------------------------------------------------------
!
!
!>
!! Computes discrete divergence of a vector field for the quadrilateral.
!!
!! Computes discrete divergence of a vector field for the quadrilateral
!! control volume formed by two adjacent triangles
!! The vector field is given by its components in the directions normal
!! to the edges.
!! The midpoint rule is used for quadrature.
!! input:  lives on edges (velocity points)
!! output: lives on edges
!!
!! @par Revision History
!! Developed  by  Jochen Foerstner, DWD (2008-07-16).
!!
SUBROUTINE div_quad_twoadjcells( vec_e, ptr_patch, ptr_int, div_vec_e,  &
  &                              opt_slev, opt_elev,                    &
                                 opt_rlstart, opt_rlend )
!
!
! patch on which computation is performed
TYPE(t_patch), TARGET, INTENT(in) :: ptr_patch

! Interpolation state
TYPE(t_int_state), INTENT(in)     :: ptr_int
!
! edge based variable of which divergence is computed
REAL(wp), INTENT(in) ::  &
  &  vec_e(:,:,:) ! dim: (nproma,nlev,nblks_e)

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_slev    ! optional vertical start level

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_elev    ! optional vertical end level

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_rlstart, opt_rlend   ! start and end values of refin_ctrl flag

!
! edge based variable in which divergence is stored
!REAL(wp), INTENT(out) ::  &
REAL(wp), INTENT(inout) ::  &
  &  div_vec_e(:,:,:) ! dim: (nproma,nlev,nblks_e)

!

INTEGER :: slev, elev     ! vertical start and end level
INTEGER :: je, jk, jb
INTEGER :: rl_start, rl_end
INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx, i_nchdom

INTEGER,  DIMENSION(:,:,:),   POINTER :: iqidx, iqblk

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
  elev = ptr_patch%nlev
END IF


IF ( PRESENT(opt_rlstart) ) THEN
  IF (opt_rlstart == 1) THEN
    CALL finish ('mo_math_operators:div_quad_twoadjcells_midpoint',  &
          &      'opt_rlstart must not be equal to 1')
  ENDIF
  rl_start = opt_rlstart
ELSE
  rl_start = 2
END IF
IF ( PRESENT(opt_rlend) ) THEN
  rl_end = opt_rlend
ELSE
  rl_end = min_rledge
END IF

! Set pointers to required index lists
iqidx => ptr_patch%edges%quad_idx
iqblk => ptr_patch%edges%quad_blk

! values for the blocking
i_nchdom   = MAX(1,ptr_patch%n_childdom)
i_startblk = ptr_patch%edges%start_blk(rl_start,1)
i_endblk   = ptr_patch%edges%end_blk(rl_end,i_nchdom)

!
! loop through all patch edges (and blocks)
!
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,je,jk,i_startidx,i_endidx)
DO jb = i_startblk, i_endblk

  CALL get_indices_e(ptr_patch, jb, i_startblk, i_endblk, &
                     i_startidx, i_endidx, rl_start, rl_end)

#ifdef __LOOP_EXCHANGE
  DO je = i_startidx, i_endidx
    DO jk = slev, elev
#else
!CDIR UNROLL=3
  DO jk = slev, elev
    DO je = i_startidx, i_endidx
#endif

      div_vec_e(je,jk,jb) =  &
        &    vec_e(iqidx(je,jb,1),jk,iqblk(je,jb,1))*ptr_int%geofac_qdiv(je,1,jb) &
        &  + vec_e(iqidx(je,jb,2),jk,iqblk(je,jb,2))*ptr_int%geofac_qdiv(je,2,jb) &
        &  + vec_e(iqidx(je,jb,3),jk,iqblk(je,jb,3))*ptr_int%geofac_qdiv(je,3,jb) &
        &  + vec_e(iqidx(je,jb,4),jk,iqblk(je,jb,4))*ptr_int%geofac_qdiv(je,4,jb)

    END DO
  END DO

END DO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE div_quad_twoadjcells

!-------------------------------------------------------------------------
!
!>
!! Computes discrete rotation.
!!
!! Computes discrete rotation at
!! (i) vertices of triangle cells (centers of dual grid cells) --or--
!! (ii) edges of hexagonal cells (centers of dual grid cells)
!! from a vector field given by its components in the directions normal
!! to triangle edges.
!! input:  lives on edges (velocity points)
!! output: lives on dual of cells (vertices for triangular grid, edges for
!!         hexagonal grid)
!!
!! @par Revision History
!! Developed and tested  by L.Bonaventura  (2002-4).
!! Adapted to new grid and patch structure by P. Korn
!! and L. Bonaventura (2005).
!! Modifications by P. Korn, MPI-M(2007-2)
!! - Switch fom array arguments to pointers
!! Modifications by Almut Gassmann, MPI-M (2007-04-20)
!! - abandon grid for the sake of patch
!! Modification by Almut Gassmann, MPI-M (2009-12-17)
!! - vorticity of hexagonal grid lives on rhombi
!!
SUBROUTINE rot_vertex_atmos( vec_e, ptr_patch, ptr_int, rot_vec, &
  &                          opt_slev, opt_elev, opt_rlstart, opt_rlend )
!
!  patch on which computation is performed
!
TYPE(t_patch), TARGET, INTENT(in) :: ptr_patch

! Interpolation state
TYPE(t_int_state), INTENT(in)     :: ptr_int
!
!  edge based variable of which rotation is computed
!
REAL(wp), INTENT(in) ::  &
  &  vec_e(:,:,:) ! dim: (nproma,nlev,nblks_e)

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_slev    ! optional vertical start level

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_elev    ! optional vertical end level

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_rlstart, opt_rlend   ! start and end values of refin_ctrl flag

!
!  vertex based variable in which rotation is stored
!
REAL(wp), INTENT(inout) ::  &
  &  rot_vec(:,:,:) ! dim: (nproma,nlev,nblks_v or nblks_e)

INTEGER :: slev, elev     ! vertical start and end level
INTEGER :: jv, jk, jb

INTEGER :: rl_start, rl_end
INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx, i_nchdom
INTEGER :: nlen, npromz_v, nblks_v

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
  elev = ptr_patch%nlev
END IF

IF ( PRESENT(opt_rlstart) ) THEN
  IF (opt_rlstart == 1) THEN
    CALL finish ('mo_math_operators:rot_vertex_atmos',  &
          &      'opt_rlstart must not be equal to 1')
  ENDIF
  rl_start = opt_rlstart
ELSE
  rl_start = 2
END IF
IF ( PRESENT(opt_rlend) ) THEN
  rl_end = opt_rlend
ELSE
  rl_end = min_rlvert
END IF

!
!  loop through over all patch vertices (and blocks)
!
! The special treatment of 2D fields is essential for efficiency on the NEC

SELECT CASE (ptr_patch%cell_type)

CASE (3)

  iidx => ptr_patch%verts%edge_idx
  iblk => ptr_patch%verts%edge_blk

  ! values for the blocking
  i_nchdom   = MAX(1,ptr_patch%n_childdom)
  i_startblk = ptr_patch%verts%start_blk(rl_start,1)
  i_endblk   = ptr_patch%verts%end_blk(rl_end,i_nchdom)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jv,jk), SCHEDULE(runtime)
  DO jb = i_startblk, i_endblk

    CALL get_indices_v(ptr_patch, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, rl_start, rl_end)

    !
    ! compute the discrete rotation for vertex jv by
    ! finite volume approximation
    ! (see Bonaventura and Ringler MWR 2005);
    ! multiplication of the vector component vec_e by
    ! the appropriate dual cell based verts%edge_orientation
    ! is required to obtain the correct value for the
    ! application of Stokes theorem (which requires the scalar
    ! product of the vector field with the tangent unit vectors
    ! going around dual cell jv COUNTERCLOKWISE;
    ! since the positive direction for the vec_e components is
    ! not necessarily the one yelding counterclockwise rotation
    ! around dual cell jv, a correction coefficient (equal to +-1)
    ! is necessary, given by g%verts%edge_orientation
    !

#ifdef __LOOP_EXCHANGE
    DO jv = i_startidx, i_endidx
      DO jk = slev, elev
#else
!CDIR UNROLL=6
    DO jk = slev, elev
      DO jv = i_startidx, i_endidx
#endif
        !
        ! calculate rotation, i.e.
        ! add individual edge contributions to rotation
        !

        rot_vec(jv,jk,jb) =   &
          vec_e(iidx(jv,jb,1),jk,iblk(jv,jb,1)) * ptr_int%geofac_rot(jv,1,jb) + &
          vec_e(iidx(jv,jb,2),jk,iblk(jv,jb,2)) * ptr_int%geofac_rot(jv,2,jb) + &
          vec_e(iidx(jv,jb,3),jk,iblk(jv,jb,3)) * ptr_int%geofac_rot(jv,3,jb) + &
          vec_e(iidx(jv,jb,4),jk,iblk(jv,jb,4)) * ptr_int%geofac_rot(jv,4,jb) + &
          vec_e(iidx(jv,jb,5),jk,iblk(jv,jb,5)) * ptr_int%geofac_rot(jv,5,jb) + &
          vec_e(iidx(jv,jb,6),jk,iblk(jv,jb,6)) * ptr_int%geofac_rot(jv,6,jb)

      END DO

    END DO

  ENDDO
!$OMP END DO
!$OMP END PARALLEL

CASE (6) ! (cell_type == 6)

  iidx => ptr_patch%verts%edge_idx
  iblk => ptr_patch%verts%edge_blk

  ! values for the blocking
  ! no grid refinement in hexagonal model
  nblks_v   = ptr_patch%nblks_int_v
  npromz_v  = ptr_patch%npromz_int_v

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,nlen,jv,jk)
  DO jb = 1, nblks_v
    IF (jb /= nblks_v) THEN
      nlen = nproma
    ELSE
      nlen = npromz_v
    ENDIF
    !
    ! Compute the discrete rotation for a triangle by
    ! application of Stokes theorem (which requires the scalar
    ! product of the vector field with the tangent unit vectors
    ! going around dual cell jv counterclockwise)
    !
#ifdef __LOOP_EXCHANGE
    DO jv = 1, nlen
      DO jk = slev, elev
#else
    DO jk = slev, elev
      DO jv = 1, nlen
#endif
        !
        ! add individual edge contributions to rotation
        !
        rot_vec(jv,jk,jb) =   &
          vec_e(iidx(jv,jb,1),jk,iblk(jv,jb,1)) * ptr_int%geofac_rot(jv,1,jb) + &
          vec_e(iidx(jv,jb,2),jk,iblk(jv,jb,2)) * ptr_int%geofac_rot(jv,2,jb) + &
          vec_e(iidx(jv,jb,3),jk,iblk(jv,jb,3)) * ptr_int%geofac_rot(jv,3,jb)

      END DO

    END DO

  ENDDO
!$OMP END DO
!$OMP END PARALLEL
END SELECT

END SUBROUTINE rot_vertex_atmos

END MODULE mo_math_divrot
