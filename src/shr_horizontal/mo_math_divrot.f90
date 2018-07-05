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
!!  - split mo_math_operators into submodules
!!  Modification by Daniel Reinert, DWD (2010-04-12)
!!  - added subroutine recon_lsq_cell_q for third order accurate least-squares
!!    reconstruction of an arbitrary field. Based on cell centered values.
!!  Modification by Daniel Reinert, DWD (2010-10-14)
!!  - added subroutine recon_lsq_cell_c for fitting a cubic polynomial in a least
!!    squares sense. Based on cell centered values.
!!  Modification by Jens-Olaf Beismann, NEC (2016-03-30), committed by R. Redler MPI-M
!!  - added an alternative for NEC SX-ACE to replace the calls of DOT_PRODUCT 
!!    The DOT_PRODUCE has once been introduced to overcome performance penalties
!!    on older NEC systems. 
!!  Modification by Will Sawyer, CSCS (2016-07-14)
!!  - added OpenACC implementation
!!
!! @par To Do
!! Boundary exchange, nblks in presence of halos and dummy edge
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
USE mo_kind,                ONLY: wp, vp
USE mo_impl_constants,      ONLY: min_rlcell, min_rledge, min_rlvert
USE mo_intp_data_strc,      ONLY: t_int_state, t_lsq
USE mo_interpol_config,     ONLY: lsq_high_set
USE mo_model_domain,        ONLY: t_patch
USE mo_grid_config,         ONLY: l_limited_area
USE mo_parallel_config,     ONLY: nproma
USE mo_exception,           ONLY: finish
USE mo_loopindices,         ONLY: get_indices_c, get_indices_e, get_indices_v
USE mo_fortran_tools,       ONLY: init, copy
#ifdef _OPENACC
USE mo_mpi,                 ONLY: i_am_accel_node
#endif

! USE mo_timer,              ONLY: timer_start, timer_stop, timer_div

IMPLICIT NONE

PRIVATE


PUBLIC :: recon_lsq_cell_l, recon_lsq_cell_l_svd, recon_lsq_cell_l_consv_svd
PUBLIC :: recon_lsq_cell_q, recon_lsq_cell_q_svd
PUBLIC :: recon_lsq_cell_cpoor, recon_lsq_cell_cpoor_svd
PUBLIC :: recon_lsq_cell_c, recon_lsq_cell_c_svd
PUBLIC :: div, div_avg
PUBLIC :: div_quad_twoadjcells
PUBLIC :: rot_vertex, rot_vertex_ri
PUBLIC :: rot_vertex_atmos

INTERFACE rot_vertex

  MODULE PROCEDURE rot_vertex_atmos

END INTERFACE


INTERFACE div

  MODULE PROCEDURE div3d
  MODULE PROCEDURE div3d_2field
  MODULE PROCEDURE div4d

END INTERFACE

#if defined( _OPENACC )
#define ACC_DEBUG NOACC
#if defined(__MATH_DIVROT_NOACC)
  LOGICAL, PARAMETER ::  acc_on = .FALSE.
#else
  LOGICAL, PARAMETER ::  acc_on = .TRUE.
#endif
#endif

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
  &                           opt_rlend, opt_lconsv )

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

  LOGICAL, INTENT(IN), OPTIONAL ::  &   !< if true, conservative reconstruction is used
    &  opt_lconsv

  REAL(wp), INTENT(INOUT) ::  &  !< cell based coefficients (geographical components)
    &  p_coeff(:,:,:,:)          !< (constant and gradients in latitudinal and
                                 !< longitudinal direction)

  REAL(wp)  ::   &               !< weights * difference of scalars i j
    &  z_d(3,nproma,ptr_patch%nlev)
  REAL(wp)  ::   &               !< matrix product of transposed Q matrix and d
    &  z_qt_times_d(2)

  INTEGER, POINTER ::   &             !< Pointer to line and block indices of
    &  iidx(:,:,:), iblk(:,:,:)       !< required stencil
  INTEGER :: slev, elev               !< vertical start and end level
  INTEGER :: jc, jk, jb               !< index of cell, vertical level and block
  INTEGER :: rl_start, rl_end
  INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx, i_nchdom
  LOGICAL :: l_consv

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
    elev = UBOUND(p_cc,2)
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
  IF ( PRESENT(opt_lconsv) ) THEN
    l_consv = opt_lconsv
  ELSE
    l_consv = .FALSE.
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
#ifdef _OPENACC
!$ACC DATA PCOPYIN( ptr_patch, ptr_int_lsq, p_cc ), PCOPY( p_coeff ), CREATE( z_d, z_qt_times_d ), IF( i_am_accel_node .AND. acc_on )
!ACC_DEBUG UPDATE DEVICE ( p_cc, p_coeff ), IF( i_am_accel_node .AND. acc_on )
!$ACC PARALLEL &
!$ACC PRESENT( ptr_patch, ptr_int_lsq, p_cc, p_coeff ), &
!$ACC PRIVATE( z_d, z_qt_times_d ), &
!$ACC IF( i_am_accel_node .AND. acc_on )

!$ACC LOOP GANG PRIVATE(i_startidx, i_endidx)
#else
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,jk,i_startidx,i_endidx,z_d,z_qt_times_d), ICON_OMP_RUNTIME_SCHEDULE
#endif
  DO jb = i_startblk, i_endblk

    CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, rl_start, rl_end)

!$ACC LOOP VECTOR COLLAPSE(2)
#ifdef __LOOP_EXCHANGE
    DO jc = i_startidx, i_endidx
      DO jk = slev, elev
#else
    DO jk = slev, elev
      DO jc = i_startidx, i_endidx
#endif

        ! note that the multiplication with lsq_weights_c(jc,js,jb) at
        ! runtime is now avoided. Instead, the multiplication with
        ! lsq_weights_c(jc,js,jb) has been shifted into the transposed
        ! Q-matrix.
        z_d(1,jc,jk) = p_cc(iidx(jc,jb,1),jk,iblk(jc,jb,1)) - p_cc(jc,jk,jb)
        z_d(2,jc,jk) = p_cc(iidx(jc,jb,2),jk,iblk(jc,jb,2)) - p_cc(jc,jk,jb)
        z_d(3,jc,jk) = p_cc(iidx(jc,jb,3),jk,iblk(jc,jb,3)) - p_cc(jc,jk,jb)

      END DO ! end loop over cells
    END DO ! end loop over vertical levels

!$ACC LOOP VECTOR COLLAPSE(2)
    DO jk = slev, elev
      DO jc = i_startidx, i_endidx

        ! matrix multiplication Q^T d (partitioned into 2 dot products)
        z_qt_times_d(1) = ptr_int_lsq%lsq_qtmat_c(jc,1,1,jb) * z_d(1,jc,jk)  &
          &             + ptr_int_lsq%lsq_qtmat_c(jc,1,2,jb) * z_d(2,jc,jk)  &
          &             + ptr_int_lsq%lsq_qtmat_c(jc,1,3,jb) * z_d(3,jc,jk)
        z_qt_times_d(2) = ptr_int_lsq%lsq_qtmat_c(jc,2,1,jb) * z_d(1,jc,jk)  &
          &             + ptr_int_lsq%lsq_qtmat_c(jc,2,2,jb) * z_d(2,jc,jk)  &
          &             + ptr_int_lsq%lsq_qtmat_c(jc,2,3,jb) * z_d(3,jc,jk)


        ! Solve linear system by backward substitution
        ! Gradient in zonal and meridional direction
        !
        ! meridional
        p_coeff(3,jc,jk,jb) = ptr_int_lsq%lsq_rmat_rdiag_c(jc,2,jb) * z_qt_times_d(2)

        ! zonal
        p_coeff(2,jc,jk,jb) = ptr_int_lsq%lsq_rmat_rdiag_c(jc,1,jb)                  &
          & * (z_qt_times_d(1) - ptr_int_lsq%lsq_rmat_utri_c(jc,1,jb)                &
          & * p_coeff(3,jc,jk,jb))

        ! constant
        p_coeff(1,jc,jk,jb) = p_cc(jc,jk,jb)

      END DO ! end loop over cells
    END DO ! end loop over vertical levels

    IF (l_consv) THEN
!$ACC LOOP VECTOR COLLAPSE(2)
      DO jk = slev, elev
        DO jc = i_startidx, i_endidx
          ! constant
          p_coeff(1,jc,jk,jb) = p_coeff(1,jc,jk,jb)                                    &
            &                 - p_coeff(2,jc,jk,jb) * ptr_int_lsq%lsq_moments(jc,jb,1) &
            &                 - p_coeff(3,jc,jk,jb) * ptr_int_lsq%lsq_moments(jc,jb,2)


        END DO ! end loop over cells
      END DO ! end loop over vertical levels
    ENDIF

  END DO ! end loop over blocks
#ifdef _OPENACC
!$ACC END PARALLEL
!ACC_DEBUG UPDATE HOST(p_coeff), IF( i_am_accel_node .AND. acc_on )
!$ACC END DATA
#else
!$OMP END DO NOWAIT
!$OMP END PARALLEL
#endif

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
  &                           opt_slev, opt_elev, opt_rlstart, opt_rlend )

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

  REAL(vp), INTENT(INOUT) ::  &  !< cell based coefficients (geographical components)
    &  p_coeff(:,:,:,:)          !< (constant and gradients in latitudinal and
                                 !< longitudinal direction)

  REAL(wp)  ::   &               !< weights * difference of scalars i j
    &  z_b(3,nproma,ptr_patch%nlev)

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
    elev = UBOUND(p_cc,2)
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
#ifdef _OPENACC
!$ACC DATA PCOPYIN( ptr_patch, ptr_int_lsq, p_cc ), PCOPY( p_coeff ), CREATE( z_b ), IF( i_am_accel_node .AND. acc_on )
!ACC_DEBUG UPDATE DEVICE ( p_cc, p_coeff ), IF( i_am_accel_node .AND. acc_on )
!$ACC PARALLEL &
!$ACC PRESENT( ptr_patch, ptr_int_lsq, p_cc, p_coeff ), &
!$ACC PRIVATE( z_b ), &
!$ACC IF( i_am_accel_node .AND. acc_on )
!$ACC LOOP GANG PRIVATE(i_startidx, i_endidx)
#else
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,jk,i_startidx,i_endidx,z_b), ICON_OMP_RUNTIME_SCHEDULE
#endif
  DO jb = i_startblk, i_endblk

    CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, rl_start, rl_end)

!$ACC LOOP VECTOR COLLAPSE(2)
#ifdef __LOOP_EXCHANGE
    DO jc = i_startidx, i_endidx
      DO jk = slev, elev
#else
    DO jk = slev, elev
      DO jc = i_startidx, i_endidx
#endif

        ! note that the multiplication with lsq_weights_c(jc,js,jb) at
        ! runtime is now avoided. Instead, the weights have been shifted 
        ! into the pseudoinverse.
        z_b(1,jc,jk) = p_cc(iidx(jc,jb,1),jk,iblk(jc,jb,1)) - p_cc(jc,jk,jb)
        z_b(2,jc,jk) = p_cc(iidx(jc,jb,2),jk,iblk(jc,jb,2)) - p_cc(jc,jk,jb)
        z_b(3,jc,jk) = p_cc(iidx(jc,jb,3),jk,iblk(jc,jb,3)) - p_cc(jc,jk,jb)

      END DO ! end loop over cells
    END DO ! end loop over vertical levels

    !
    ! 2. compute cell based coefficients for linear reconstruction
    !    calculate matrix vector product PINV(A) * b
    !
!$ACC LOOP VECTOR COLLAPSE(2)
    DO jk = slev, elev
      DO jc = i_startidx, i_endidx

        ! meridional
        p_coeff(2,jc,jk,jb) = ptr_int_lsq%lsq_pseudoinv(jc,2,1,jb) * z_b(1,jc,jk)  &
          &                 + ptr_int_lsq%lsq_pseudoinv(jc,2,2,jb) * z_b(2,jc,jk)  &
          &                 + ptr_int_lsq%lsq_pseudoinv(jc,2,3,jb) * z_b(3,jc,jk)

        ! zonal
        p_coeff(1,jc,jk,jb) = ptr_int_lsq%lsq_pseudoinv(jc,1,1,jb) * z_b(1,jc,jk)  &
          &                 + ptr_int_lsq%lsq_pseudoinv(jc,1,2,jb) * z_b(2,jc,jk)  &
          &                 + ptr_int_lsq%lsq_pseudoinv(jc,1,3,jb) * z_b(3,jc,jk)

      END DO ! end loop over cells
    END DO ! end loop over vertical levels


  END DO ! end loop over blocks
#ifdef _OPENACC
!$ACC END PARALLEL
!ACC_DEBUG UPDATE HOST(p_coeff), IF( i_am_accel_node .AND. acc_on )
!$ACC END DATA
#else
!$OMP END DO NOWAIT
!$OMP END PARALLEL
#endif

END SUBROUTINE recon_lsq_cell_l_svd


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
SUBROUTINE recon_lsq_cell_l_consv_svd( p_cc, ptr_patch, ptr_int_lsq, p_coeff, &
  &                                    opt_slev, opt_elev, opt_rlstart,      &
  &                                    opt_rlend, opt_lconsv )

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

  LOGICAL, INTENT(IN), OPTIONAL ::  &   !< if true, conservative reconstruction is used
    &  opt_lconsv

  REAL(wp), INTENT(INOUT) ::  &  !< cell based coefficients (geographical components)
    &  p_coeff(:,:,:,:)          !< (constant and gradients in latitudinal and
                                 !< longitudinal direction)

  REAL(wp)  ::   &               !< weights * difference of scalars i j
    &  z_b(3,nproma,ptr_patch%nlev)

  INTEGER, POINTER ::   &            !< Pointer to line and block indices of
    &  iidx(:,:,:), iblk(:,:,:)      !< required stencil
  INTEGER :: slev, elev              !< vertical start and end level
  INTEGER :: jc, jk, jb              !< index of cell, vertical level and block
  INTEGER :: rl_start, rl_end
  INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx, i_nchdom
  LOGICAL :: l_consv

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
    elev = UBOUND(p_cc,2)
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
  IF ( PRESENT(opt_lconsv) ) THEN
    l_consv = opt_lconsv
  ELSE
    l_consv = .FALSE.
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
#ifdef _OPENACC
!$ACC DATA PCOPYIN( ptr_patch, ptr_int_lsq, p_cc ), PCOPY( p_coeff ), IF( i_am_accel_node .AND. acc_on )
!ACC_DEBUG UPDATE DEVICE ( p_cc, p_coeff ), IF( i_am_accel_node .AND. acc_on )
!$ACC PARALLEL &
!$ACC PRESENT( ptr_patch, ptr_int_lsq, p_cc, p_coeff ), &
!$ACC PRIVATE( z_b ), &
!$ACC IF( i_am_accel_node .AND. acc_on )

!$ACC LOOP GANG PRIVATE(i_startidx, i_endidx)
#else
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,jk,i_startidx,i_endidx,z_b), ICON_OMP_RUNTIME_SCHEDULE
#endif
  DO jb = i_startblk, i_endblk

    CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, rl_start, rl_end)

!$ACC LOOP VECTOR COLLAPSE(2)
#ifdef __LOOP_EXCHANGE
    DO jc = i_startidx, i_endidx
      DO jk = slev, elev
#else
    DO jk = slev, elev
      DO jc = i_startidx, i_endidx
#endif

        ! note that the multiplication with lsq_weights_c(jc,js,jb) at
        ! runtime is now avoided. Instead, the weights have been shifted 
        ! into the pseudoinverse.
        z_b(1,jc,jk) = p_cc(iidx(jc,jb,1),jk,iblk(jc,jb,1)) - p_cc(jc,jk,jb)
        z_b(2,jc,jk) = p_cc(iidx(jc,jb,2),jk,iblk(jc,jb,2)) - p_cc(jc,jk,jb)
        z_b(3,jc,jk) = p_cc(iidx(jc,jb,3),jk,iblk(jc,jb,3)) - p_cc(jc,jk,jb)

      END DO ! end loop over cells
    END DO ! end loop over vertical levels

    !
    ! 2. compute cell based coefficients for linear reconstruction
    !    calculate matrix vector product PINV(A) * b
    ! 
!$ACC LOOP VECTOR COLLAPSE(2)
    DO jk = slev, elev
      DO jc = i_startidx, i_endidx

        ! meridional
        p_coeff(3,jc,jk,jb) = ptr_int_lsq%lsq_pseudoinv(jc,2,1,jb) * z_b(1,jc,jk)  &
          &                 + ptr_int_lsq%lsq_pseudoinv(jc,2,2,jb) * z_b(2,jc,jk)  &
          &                 + ptr_int_lsq%lsq_pseudoinv(jc,2,3,jb) * z_b(3,jc,jk)

        ! zonal
        p_coeff(2,jc,jk,jb) = ptr_int_lsq%lsq_pseudoinv(jc,1,1,jb) * z_b(1,jc,jk)  &
          &                 + ptr_int_lsq%lsq_pseudoinv(jc,1,2,jb) * z_b(2,jc,jk)  &
          &                 + ptr_int_lsq%lsq_pseudoinv(jc,1,3,jb) * z_b(3,jc,jk)

        ! constant
        p_coeff(1,jc,jk,jb) = p_cc(jc,jk,jb)

      END DO ! end loop over cells
    END DO ! end loop over vertical levels

    IF (l_consv) THEN
      DO jk = slev, elev
        DO jc = i_startidx, i_endidx

          ! In the case of a conservative reconstruction, 
          ! the coefficient c0 is derived from the linear constraint
          !
          p_coeff(1,jc,jk,jb) = p_coeff(1,jc,jk,jb)                                    &
            &                 - p_coeff(2,jc,jk,jb) * ptr_int_lsq%lsq_moments(jc,jb,1) &
            &                 - p_coeff(3,jc,jk,jb) * ptr_int_lsq%lsq_moments(jc,jb,2)

        END DO ! end loop over cells
      END DO ! end loop over vertical levels
    ENDIF

  END DO ! end loop over blocks
#ifdef _OPENACC
!$ACC END PARALLEL
!ACC_DEBUG UPDATE HOST(p_coeff), IF( i_am_accel_node .AND. acc_on )
!$ACC END DATA
#else
!$OMP END DO NOWAIT
!$OMP END PARALLEL
#endif

END SUBROUTINE recon_lsq_cell_l_consv_svd



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
    elev = UBOUND(p_cc,2)
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


#ifdef _OPENACC
!$ACC DATA PCOPYIN( ptr_patch, ptr_int_lsq, p_cc ), PCOPY( p_coeff ), &
!$ACC      CREATE(  z_d, z_qt_times_d ), IF( i_am_accel_node .AND. acc_on )
!ACC_DEBUG UPDATE DEVICE ( p_cc, ptr_rrdiag, ptr_rutri, p_coeff ), IF( i_am_accel_node .AND. acc_on )
#else
!$OMP PARALLEL
#endif

  IF (ptr_patch%id > 1 .OR. l_limited_area) THEN
    CALL init(p_coeff(:,:,1:6,1:i_startblk))
#ifndef _OPENACC
!$OMP BARRIER
#endif
  ENDIF

#ifdef _OPENACC
!$ACC PARALLEL &
!$ACC PRESENT( ptr_patch, ptr_int_lsq, p_cc, p_coeff ), &
!$ACC PRIVATE( z_d, z_qt_times_d ), &
!$ACC IF( i_am_accel_node .AND. acc_on )

!$ACC LOOP GANG PRIVATE(i_startidx, i_endidx)
#else
!$OMP DO PRIVATE(jb,jc,jk,i_startidx,i_endidx,z_d,z_qt_times_d), ICON_OMP_RUNTIME_SCHEDULE
#endif
  DO jb = i_startblk, i_endblk

    CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, rl_start, rl_end)

    !
    ! 1. compute right hand side of linear system
    !
!$ACC LOOP VECTOR COLLAPSE(2)
#ifdef __LOOP_EXCHANGE
    DO jc = i_startidx, i_endidx
      DO jk = slev, elev
#else
    DO jk = slev, elev
      DO jc = i_startidx, i_endidx
#endif

        z_d(1,jc,jk) = p_cc(iidx(jc,jb,1),jk,iblk(jc,jb,1)) - p_cc(jc,jk,jb)
        z_d(2,jc,jk) = p_cc(iidx(jc,jb,2),jk,iblk(jc,jb,2)) - p_cc(jc,jk,jb)
        z_d(3,jc,jk) = p_cc(iidx(jc,jb,3),jk,iblk(jc,jb,3)) - p_cc(jc,jk,jb)
        z_d(4,jc,jk) = p_cc(iidx(jc,jb,4),jk,iblk(jc,jb,4)) - p_cc(jc,jk,jb)
        z_d(5,jc,jk) = p_cc(iidx(jc,jb,5),jk,iblk(jc,jb,5)) - p_cc(jc,jk,jb)
        z_d(6,jc,jk) = p_cc(iidx(jc,jb,6),jk,iblk(jc,jb,6)) - p_cc(jc,jk,jb)
        z_d(7,jc,jk) = p_cc(iidx(jc,jb,7),jk,iblk(jc,jb,7)) - p_cc(jc,jk,jb)
        z_d(8,jc,jk) = p_cc(iidx(jc,jb,8),jk,iblk(jc,jb,8)) - p_cc(jc,jk,jb)
        z_d(9,jc,jk) = p_cc(iidx(jc,jb,9),jk,iblk(jc,jb,9)) - p_cc(jc,jk,jb)

      ENDDO
    ENDDO



    !
    ! 2. compute cell based coefficients for quadratic reconstruction
    !
!$ACC LOOP VECTOR COLLAPSE(2)
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
        p_coeff(6,jc,jk,jb) = ptr_rrdiag(jc,5,jb) * z_qt_times_d(5)
        p_coeff(5,jc,jk,jb) = ptr_rrdiag(jc,4,jb)                                         &
          &                 * ( z_qt_times_d(4) - ptr_rutri(jc,1,jb)*p_coeff(6,jc,jk,jb) )
        p_coeff(4,jc,jk,jb) = ptr_rrdiag(jc,3,jb)                                         &
          &                 * ( z_qt_times_d(3) - ptr_rutri(jc,2,jb)*p_coeff(5,jc,jk,jb)  &
          &                 - ptr_rutri(jc,3,jb) * p_coeff(6,jc,jk,jb) )
        p_coeff(3,jc,jk,jb) = ptr_rrdiag(jc,2,jb)                                         &
          &                 * ( z_qt_times_d(2) - ptr_rutri(jc,4,jb)*p_coeff(4,jc,jk,jb)  &
          &                 - ptr_rutri(jc,5,jb) * p_coeff(5,jc,jk,jb)                    &
          &                 - ptr_rutri(jc,6,jb) * p_coeff(6,jc,jk,jb) )
        p_coeff(2,jc,jk,jb) = ptr_rrdiag(jc,1,jb)                                         &
          &                 * ( z_qt_times_d(1) - ptr_rutri(jc,7,jb)*p_coeff(3,jc,jk,jb)  &
          &                 - ptr_rutri(jc,8,jb) * p_coeff(4,jc,jk,jb)                    &
          &                 - ptr_rutri(jc,9,jb) * p_coeff(5,jc,jk,jb)                    &
          &                 - ptr_rutri(jc,10,jb)* p_coeff(6,jc,jk,jb) )

        p_coeff(1,jc,jk,jb) = p_cc(jc,jk,jb)                                              &
          &                  - p_coeff(2,jc,jk,jb) * ptr_int_lsq%lsq_moments(jc,jb,1)     &
          &                  - p_coeff(3,jc,jk,jb) * ptr_int_lsq%lsq_moments(jc,jb,2)     &
          &                  - p_coeff(4,jc,jk,jb) * ptr_int_lsq%lsq_moments(jc,jb,3)     &
          &                  - p_coeff(5,jc,jk,jb) * ptr_int_lsq%lsq_moments(jc,jb,4)     &
          &                  - p_coeff(6,jc,jk,jb) * ptr_int_lsq%lsq_moments(jc,jb,5)


      END DO ! end loop over cells

    END DO ! end loop over vertical levels

  END DO ! end loop over blocks
#ifdef _OPENACC
!$ACC END PARALLEL
!ACC_DEBUG UPDATE HOST(p_coeff), IF( i_am_accel_node .AND. acc_on )
!$ACC END DATA
#else
!$OMP END DO NOWAIT
!$OMP END PARALLEL
#endif


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
    elev = UBOUND(p_cc,2)
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



#ifdef _OPENACC
!$ACC DATA PCOPYIN( ptr_patch, ptr_int_lsq, p_cc ), PCOPY( p_coeff ), IF( i_am_accel_node .AND. acc_on )
!ACC_DEBUG UPDATE DEVICE( p_cc, p_coeff ), IF( i_am_accel_node .AND. acc_on ) 
#else
!$OMP PARALLEL
#endif

  IF (ptr_patch%id > 1 .OR. l_limited_area) THEN
    CALL init(p_coeff(:,:,1:6,1:i_startblk))
#ifndef _OPENACC
!$OMP BARRIER
#endif
  ENDIF

#ifdef _OPENACC
!$ACC PARALLEL &
!$ACC PRESENT( ptr_patch, ptr_int_lsq, p_cc, p_coeff ), &
!$ACC PRIVATE( z_b ), &
!$ACC IF( i_am_accel_node .AND. acc_on )

!$ACC LOOP GANG PRIVATE(i_startidx, i_endidx)
#else
!$OMP DO PRIVATE(jb,jc,jk,i_startidx,i_endidx,z_b), ICON_OMP_RUNTIME_SCHEDULE
#endif
  DO jb = i_startblk, i_endblk

    CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, rl_start, rl_end)

    !
    ! 1. compute right hand side of linear system
    !

!$ACC LOOP VECTOR COLLAPSE(2)
#ifdef __LOOP_EXCHANGE
    DO jc = i_startidx, i_endidx
      DO jk = slev, elev
#else
    DO jk = slev, elev
      DO jc = i_startidx, i_endidx
#endif

        z_b(1,jc,jk) = p_cc(iidx(jc,jb,1),jk,iblk(jc,jb,1)) - p_cc(jc,jk,jb)
        z_b(2,jc,jk) = p_cc(iidx(jc,jb,2),jk,iblk(jc,jb,2)) - p_cc(jc,jk,jb)
        z_b(3,jc,jk) = p_cc(iidx(jc,jb,3),jk,iblk(jc,jb,3)) - p_cc(jc,jk,jb)
        z_b(4,jc,jk) = p_cc(iidx(jc,jb,4),jk,iblk(jc,jb,4)) - p_cc(jc,jk,jb)
        z_b(5,jc,jk) = p_cc(iidx(jc,jb,5),jk,iblk(jc,jb,5)) - p_cc(jc,jk,jb)
        z_b(6,jc,jk) = p_cc(iidx(jc,jb,6),jk,iblk(jc,jb,6)) - p_cc(jc,jk,jb)
        z_b(7,jc,jk) = p_cc(iidx(jc,jb,7),jk,iblk(jc,jb,7)) - p_cc(jc,jk,jb)
        z_b(8,jc,jk) = p_cc(iidx(jc,jb,8),jk,iblk(jc,jb,8)) - p_cc(jc,jk,jb)
        z_b(9,jc,jk) = p_cc(iidx(jc,jb,9),jk,iblk(jc,jb,9)) - p_cc(jc,jk,jb)

      ENDDO
    ENDDO



    !
    ! 2. compute cell based coefficients for quadratic reconstruction
    !    calculate matrix vector product PINV(A) * b
    !
!$ACC LOOP VECTOR COLLAPSE(2)
    DO jk = slev, elev

      DO jc = i_startidx, i_endidx

        ! (intrinsic function matmul not applied, due to massive
        ! performance penalty on the NEC. Instead the intrinsic dot product
        ! function is used.
!CDIR BEGIN EXPAND=9
        p_coeff(6,jc,jk,jb) = DOT_PRODUCT(ptr_int_lsq%lsq_pseudoinv(jc,5,1:9,jb), &
          &                               z_b(1:9,jc,jk))
        p_coeff(5,jc,jk,jb) = DOT_PRODUCT(ptr_int_lsq%lsq_pseudoinv(jc,4,1:9,jb), &
          &                               z_b(1:9,jc,jk))
        p_coeff(4,jc,jk,jb) = DOT_PRODUCT(ptr_int_lsq%lsq_pseudoinv(jc,3,1:9,jb), &
          &                               z_b(1:9,jc,jk))
        p_coeff(3,jc,jk,jb) = DOT_PRODUCT(ptr_int_lsq%lsq_pseudoinv(jc,2,1:9,jb), &
          &                               z_b(1:9,jc,jk))
        p_coeff(2,jc,jk,jb) = DOT_PRODUCT(ptr_int_lsq%lsq_pseudoinv(jc,1,1:9,jb), &
          &                               z_b(1:9,jc,jk))
!CDIR END


        ! At the end, the coefficient c0 is derived from the linear constraint
        !
        p_coeff(1,jc,jk,jb) = p_cc(jc,jk,jb) - DOT_PRODUCT(p_coeff(2:6,jc,jk,jb), &
          &                   ptr_int_lsq%lsq_moments(jc,jb,1:5))


      END DO ! end loop over cells

    END DO ! end loop over vertical levels

  END DO ! end loop over blocks
#ifdef _OPENACC
!$ACC END PARALLEL
!ACC_DEBUG UPDATE HOST(p_coeff) IF( i_am_accel_node .AND. acc_on )
!$ACC END DATA
#else
!$OMP END DO NOWAIT
!$OMP END PARALLEL
#endif

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
    elev = UBOUND(p_cc,2)
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



#ifdef _OPENACC
!$ACC DATA PCOPYIN(  ptr_patch, ptr_int_lsq, p_cc ), PCOPY( p_coeff ), CREATE( z_d, z_qt_times_d ), IF( i_am_accel_node .AND. acc_on )
!ACC_DEBUG UPDATE DEVICE( p_cc, p_coeff ), IF( i_am_accel_node .AND. acc_on )
#else
!$OMP PARALLEL
#endif

  IF (ptr_patch%id > 1 .OR. l_limited_area) THEN
    CALL init(p_coeff(:,:,1:8,1:i_startblk))
#ifndef _OPENACC
!$OMP BARRIER
#endif
  ENDIF

#ifdef _OPENACC
!$ACC PARALLEL &
!$ACC PRESENT( ptr_patch, ptr_int_lsq, p_cc, p_coeff ), &
!$ACC PRIVATE( z_d, z_qt_times_d ), &
!$ACC IF( i_am_accel_node .AND. acc_on )

!$ACC LOOP GANG PRIVATE(i_startidx, i_endidx)
#else
!$OMP DO PRIVATE(jb,jc,jk,i_startidx,i_endidx,z_d,z_qt_times_d), ICON_OMP_RUNTIME_SCHEDULE
#endif
  DO jb = i_startblk, i_endblk

    CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, rl_start, rl_end)

    !
    ! 1. compute right hand side of linear system
    !

!$ACC LOOP VECTOR COLLAPSE(2)
#ifdef __LOOP_EXCHANGE
    DO jc = i_startidx, i_endidx
      DO jk = slev, elev
#else
    DO jk = slev, elev
      DO jc = i_startidx, i_endidx
#endif

        z_d(1,jc,jk) = p_cc(iidx(jc,jb,1),jk,iblk(jc,jb,1)) - p_cc(jc,jk,jb)
        z_d(2,jc,jk) = p_cc(iidx(jc,jb,2),jk,iblk(jc,jb,2)) - p_cc(jc,jk,jb)
        z_d(3,jc,jk) = p_cc(iidx(jc,jb,3),jk,iblk(jc,jb,3)) - p_cc(jc,jk,jb)
        z_d(4,jc,jk) = p_cc(iidx(jc,jb,4),jk,iblk(jc,jb,4)) - p_cc(jc,jk,jb)
        z_d(5,jc,jk) = p_cc(iidx(jc,jb,5),jk,iblk(jc,jb,5)) - p_cc(jc,jk,jb)
        z_d(6,jc,jk) = p_cc(iidx(jc,jb,6),jk,iblk(jc,jb,6)) - p_cc(jc,jk,jb)
        z_d(7,jc,jk) = p_cc(iidx(jc,jb,7),jk,iblk(jc,jb,7)) - p_cc(jc,jk,jb)
        z_d(8,jc,jk) = p_cc(iidx(jc,jb,8),jk,iblk(jc,jb,8)) - p_cc(jc,jk,jb)
        z_d(9,jc,jk) = p_cc(iidx(jc,jb,9),jk,iblk(jc,jb,9)) - p_cc(jc,jk,jb)

      ENDDO
    ENDDO


    !
    ! 2. compute cell based coefficients for quadratic reconstruction
    !
!$ACC LOOP VECTOR COLLAPSE(2)
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
        p_coeff(8,jc,jk,jb) = ptr_rrdiag(jc,7,jb) * z_qt_times_d(7)
        p_coeff(7,jc,jk,jb) = ptr_rrdiag(jc,6,jb)                                         &
          &                 * ( z_qt_times_d(6) - ptr_rutri(jc,1,jb)*p_coeff(8,jc,jk,jb) )
        p_coeff(6,jc,jk,jb) = ptr_rrdiag(jc,5,jb)                                         &
          &                 * ( z_qt_times_d(5) - ptr_rutri(jc,2,jb)*p_coeff(7,jc,jk,jb)  &
          &                 - ptr_rutri(jc,3,jb) * p_coeff(8,jc,jk,jb) )
        p_coeff(5,jc,jk,jb) = ptr_rrdiag(jc,4,jb)                                         &
          &                 * ( z_qt_times_d(4) - ptr_rutri(jc,4,jb)*p_coeff(6,jc,jk,jb)  &
          &                 - ptr_rutri(jc,5,jb) * p_coeff(7,jc,jk,jb)                    &
          &                 - ptr_rutri(jc,6,jb) * p_coeff(8,jc,jk,jb) )
        p_coeff(4,jc,jk,jb) = ptr_rrdiag(jc,3,jb)                                         &
          &                 * ( z_qt_times_d(3) - ptr_rutri(jc,7,jb)*p_coeff(5,jc,jk,jb)  &
          &                 - ptr_rutri(jc,8,jb) * p_coeff(6,jc,jk,jb)                    &
          &                 - ptr_rutri(jc,9,jb) * p_coeff(7,jc,jk,jb)                    &
          &                 - ptr_rutri(jc,10,jb)* p_coeff(8,jc,jk,jb) )
        p_coeff(3,jc,jk,jb) = ptr_rrdiag(jc,2,jb)                                         &
          &                 * ( z_qt_times_d(2) - ptr_rutri(jc,11,jb)*p_coeff(4,jc,jk,jb) &
          &                 - ptr_rutri(jc,12,jb) * p_coeff(5,jc,jk,jb)                   &
          &                 - ptr_rutri(jc,13,jb) * p_coeff(6,jc,jk,jb)                   &
          &                 - ptr_rutri(jc,14,jb) * p_coeff(7,jc,jk,jb)                   &
          &                 - ptr_rutri(jc,15,jb) * p_coeff(8,jc,jk,jb) )
        p_coeff(2,jc,jk,jb) = ptr_rrdiag(jc,1,jb)                                         &
          &                 * ( z_qt_times_d(1) - ptr_rutri(jc,16,jb)*p_coeff(3,jc,jk,jb) &
          &                 - ptr_rutri(jc,17,jb) * p_coeff(4,jc,jk,jb)                   &
          &                 - ptr_rutri(jc,18,jb) * p_coeff(5,jc,jk,jb)                   &
          &                 - ptr_rutri(jc,19,jb) * p_coeff(6,jc,jk,jb)                   &
          &                 - ptr_rutri(jc,20,jb) * p_coeff(7,jc,jk,jb)                   &
          &                 - ptr_rutri(jc,21,jb) * p_coeff(8,jc,jk,jb) )


        p_coeff(1,jc,jk,jb) = p_cc(jc,jk,jb)                                              &
          &                  - p_coeff(2,jc,jk,jb) * ptr_int_lsq%lsq_moments(jc,jb,1)     &
          &                  - p_coeff(3,jc,jk,jb) * ptr_int_lsq%lsq_moments(jc,jb,2)     &
          &                  - p_coeff(4,jc,jk,jb) * ptr_int_lsq%lsq_moments(jc,jb,3)     &
          &                  - p_coeff(5,jc,jk,jb) * ptr_int_lsq%lsq_moments(jc,jb,4)     &
          &                  - p_coeff(6,jc,jk,jb) * ptr_int_lsq%lsq_moments(jc,jb,5)     &
          &                  - p_coeff(7,jc,jk,jb) * ptr_int_lsq%lsq_moments(jc,jb,6)     &
          &                  - p_coeff(8,jc,jk,jb) * ptr_int_lsq%lsq_moments(jc,jb,7)


      END DO ! end loop over cells

    END DO ! end loop over vertical levels

  END DO ! end loop over blocks
#ifdef _OPENACC
!$ACC END PARALLEL
!ACC_DEBUG UPDATE HOST(p_coeff), IF( i_am_accel_node .AND. acc_on )
!$ACC END DATA
#else
!$OMP END DO NOWAIT
!$OMP END PARALLEL
#endif


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
    elev = UBOUND(p_cc,2)
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

#ifdef _OPENACC
!$ACC DATA PCOPYIN( ptr_patch, ptr_int_lsq, p_cc ), PCOPY( p_coeff ), CREATE( z_b ), IF( i_am_accel_node .AND. acc_on )
!ACC_DEBUG UPDATE DEVICE( p_cc, p_coeff ), IF( i_am_accel_node .AND. acc_on )
#else
!$OMP PARALLEL
#endif

  IF (ptr_patch%id > 1 .OR. l_limited_area) THEN
    CALL init(p_coeff(:,:,:,1:i_startblk))
#ifndef _OPENACC
!$OMP BARRIER
#endif
  ENDIF


#ifdef _OPENACC
!$ACC PARALLEL &
!$ACC PRESENT( ptr_patch, ptr_int_lsq, p_cc, p_coeff ), &
!$ACC PRIVATE( z_b ), &
!$ACC IF( i_am_accel_node .AND. acc_on )

!$ACC LOOP GANG PRIVATE(i_startidx, i_endidx)
#else
!$OMP DO PRIVATE(jb,jc,jk,i_startidx,i_endidx,z_b), ICON_OMP_RUNTIME_SCHEDULE
#endif
  DO jb = i_startblk, i_endblk

    CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, rl_start, rl_end)

    !
    ! 1. compute right hand side of linear system
    !
!$ACC LOOP VECTOR COLLAPSE(2)
#ifdef __LOOP_EXCHANGE
    DO jc = i_startidx, i_endidx
      DO jk = slev, elev
#else
    DO jk = slev, elev
      DO jc = i_startidx, i_endidx
#endif

        z_b(1,jc,jk) = p_cc(iidx(jc,jb,1),jk,iblk(jc,jb,1)) - p_cc(jc,jk,jb)
        z_b(2,jc,jk) = p_cc(iidx(jc,jb,2),jk,iblk(jc,jb,2)) - p_cc(jc,jk,jb)
        z_b(3,jc,jk) = p_cc(iidx(jc,jb,3),jk,iblk(jc,jb,3)) - p_cc(jc,jk,jb)
        z_b(4,jc,jk) = p_cc(iidx(jc,jb,4),jk,iblk(jc,jb,4)) - p_cc(jc,jk,jb)
        z_b(5,jc,jk) = p_cc(iidx(jc,jb,5),jk,iblk(jc,jb,5)) - p_cc(jc,jk,jb)
        z_b(6,jc,jk) = p_cc(iidx(jc,jb,6),jk,iblk(jc,jb,6)) - p_cc(jc,jk,jb)
        z_b(7,jc,jk) = p_cc(iidx(jc,jb,7),jk,iblk(jc,jb,7)) - p_cc(jc,jk,jb)
        z_b(8,jc,jk) = p_cc(iidx(jc,jb,8),jk,iblk(jc,jb,8)) - p_cc(jc,jk,jb)
        z_b(9,jc,jk) = p_cc(iidx(jc,jb,9),jk,iblk(jc,jb,9)) - p_cc(jc,jk,jb)

      ENDDO
    ENDDO


    !
    ! 2. compute cell based coefficients for poor man's cubic reconstruction
    !    calculate matrix vector product PINV(A) * b
    !
!$ACC LOOP VECTOR COLLAPSE(2)
    DO jk = slev, elev

      DO jc = i_startidx, i_endidx

        ! (intrinsic function matmul not applied, due to massive
        ! performance penalty on the NEC. Instead the intrinsic dot product
        ! function is used.
!CDIR BEGIN EXPAND=9
        p_coeff(8,jc,jk,jb) = DOT_PRODUCT(ptr_int_lsq%lsq_pseudoinv(jc,7,1:9,jb), &
          &                               z_b(1:9,jc,jk))
        p_coeff(7,jc,jk,jb) = DOT_PRODUCT(ptr_int_lsq%lsq_pseudoinv(jc,6,1:9,jb), &
          &                               z_b(1:9,jc,jk))
        p_coeff(6,jc,jk,jb) = DOT_PRODUCT(ptr_int_lsq%lsq_pseudoinv(jc,5,1:9,jb), &
          &                               z_b(1:9,jc,jk))
        p_coeff(5,jc,jk,jb) = DOT_PRODUCT(ptr_int_lsq%lsq_pseudoinv(jc,4,1:9,jb), &
          &                               z_b(1:9,jc,jk))
        p_coeff(4,jc,jk,jb) = DOT_PRODUCT(ptr_int_lsq%lsq_pseudoinv(jc,3,1:9,jb), &
          &                               z_b(1:9,jc,jk))
        p_coeff(3,jc,jk,jb) = DOT_PRODUCT(ptr_int_lsq%lsq_pseudoinv(jc,2,1:9,jb), &
          &                               z_b(1:9,jc,jk))
        p_coeff(2,jc,jk,jb) = DOT_PRODUCT(ptr_int_lsq%lsq_pseudoinv(jc,1,1:9,jb), &
          &                               z_b(1:9,jc,jk))
!CDIR END

        ! At the end, the coefficient c0 is derived from the linear constraint
        !
!CDIR EXPAND=7
        p_coeff(1,jc,jk,jb) = p_cc(jc,jk,jb) - DOT_PRODUCT(p_coeff(2:8,jc,jk,jb), &
          &                   ptr_int_lsq%lsq_moments(jc,jb,1:7))

      END DO ! end loop over cells

    END DO ! end loop over vertical levels

  END DO ! end loop over blocks
#ifdef _OPENACC
!$ACC END PARALLEL
!ACC_DEBUG UPDATE HOST(p_coeff) IF( i_am_accel_node .AND. acc_on )
!$ACC END DATA
#else
!$OMP END DO NOWAIT
!$OMP END PARALLEL
#endif

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
    elev = UBOUND(p_cc,2)
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


#ifdef _OPENACC
!$ACC DATA PCOPYIN( ptr_patch, ptr_int_lsq, p_cc ), PCOPY( p_coeff ), CREATE( z_d,z_qt_times_d ), IF( i_am_accel_node .AND. acc_on )
!ACC_DEBUG UPDATE DEVICE( p_cc, p_coeff ), IF( i_am_accel_node .AND. acc_on )
#else
!$OMP PARALLEL
#endif

  IF (ptr_patch%id > 1 .OR. l_limited_area) THEN
    CALL init(p_coeff(:,:,1:10,1:i_startblk))
#ifndef _OPENACC
!$OMP BARRIER
#endif
  ENDIF

#ifdef _OPENACC
!$ACC PARALLEL &
!$ACC PRESENT( ptr_patch, p_cc, ptr_int_lsq, p_coeff ), &
!$ACC PRIVATE( z_d, z_qt_times_d ), &
!$ACC IF( i_am_accel_node .AND. acc_on )

!$ACC LOOP GANG PRIVATE(i_startidx, i_endidx)
#else
!$OMP DO PRIVATE(jb,jc,jk,i_startidx,i_endidx,z_d,z_qt_times_d), ICON_OMP_RUNTIME_SCHEDULE
#endif
  DO jb = i_startblk, i_endblk

    CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, rl_start, rl_end)

    !
    ! 1. compute right hand side of linear system
    !
!$ACC LOOP VECTOR COLLAPSE(2)
#ifdef __LOOP_EXCHANGE
    DO jc = i_startidx, i_endidx
      DO jk = slev, elev
#else
    DO jk = slev, elev
      DO jc = i_startidx, i_endidx
#endif

        z_d(1,jc,jk) = p_cc(iidx(jc,jb,1),jk,iblk(jc,jb,1)) - p_cc(jc,jk,jb)
        z_d(2,jc,jk) = p_cc(iidx(jc,jb,2),jk,iblk(jc,jb,2)) - p_cc(jc,jk,jb)
        z_d(3,jc,jk) = p_cc(iidx(jc,jb,3),jk,iblk(jc,jb,3)) - p_cc(jc,jk,jb)
        z_d(4,jc,jk) = p_cc(iidx(jc,jb,4),jk,iblk(jc,jb,4)) - p_cc(jc,jk,jb)
        z_d(5,jc,jk) = p_cc(iidx(jc,jb,5),jk,iblk(jc,jb,5)) - p_cc(jc,jk,jb)
        z_d(6,jc,jk) = p_cc(iidx(jc,jb,6),jk,iblk(jc,jb,6)) - p_cc(jc,jk,jb)
        z_d(7,jc,jk) = p_cc(iidx(jc,jb,7),jk,iblk(jc,jb,7)) - p_cc(jc,jk,jb)
        z_d(8,jc,jk) = p_cc(iidx(jc,jb,8),jk,iblk(jc,jb,8)) - p_cc(jc,jk,jb)
        z_d(9,jc,jk) = p_cc(iidx(jc,jb,9),jk,iblk(jc,jb,9)) - p_cc(jc,jk,jb)

      ENDDO
    ENDDO


    !
    ! 2. compute cell based coefficients for quadratic reconstruction
    !
!$ACC LOOP VECTOR COLLAPSE(2)
    DO jk = slev, elev

      DO jc = i_startidx, i_endidx

        ! calculate matrix vector product Q^T d (transposed of Q times LHS)
        ! (intrinsic function matmul not applied, due to massive
        ! performance penalty on the NEC. Instead the intrinsic dot product
        ! function is applied
!CDIR BEGIN EXPAND=9
!TODO:  these should be nine scalars, since they should reside in registers
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
        p_coeff(10,jc,jk,jb) = ptr_rrdiag(jc,9,jb) * z_qt_times_d(9)
        p_coeff(9,jc,jk,jb)  = ptr_rrdiag(jc,8,jb)                                         &
          &                  * ( z_qt_times_d(8) - ptr_rutri(jc,1,jb)*p_coeff(10,jc,jk,jb) )
        p_coeff(8,jc,jk,jb)  = ptr_rrdiag(jc,7,jb)                                         &
          &                  * ( z_qt_times_d(7) - ptr_rutri(jc,2,jb)*p_coeff(9,jc,jk,jb)  &
          &                  - ptr_rutri(jc,3,jb) * p_coeff(10,jc,jk,jb) )
        p_coeff(7,jc,jk,jb)  = ptr_rrdiag(jc,6,jb)                                         &
          &                  * ( z_qt_times_d(6) - ptr_rutri(jc,4,jb)*p_coeff(8,jc,jk,jb)  &
          &                  - ptr_rutri(jc,5,jb) * p_coeff(9,jc,jk,jb)                    &
          &                  - ptr_rutri(jc,6,jb) * p_coeff(10,jc,jk,jb) )
        p_coeff(6,jc,jk,jb)  = ptr_rrdiag(jc,5,jb)                                         &
          &                  * ( z_qt_times_d(5) - ptr_rutri(jc,7,jb)*p_coeff(7,jc,jk,jb)  &
          &                  - ptr_rutri(jc,8,jb) * p_coeff(8,jc,jk,jb)                    &
          &                  - ptr_rutri(jc,9,jb) * p_coeff(9,jc,jk,jb)                    &
          &                  - ptr_rutri(jc,10,jb)* p_coeff(10,jc,jk,jb) )
        p_coeff(5,jc,jk,jb)  = ptr_rrdiag(jc,4,jb)                                         &
          &                  * ( z_qt_times_d(4) - ptr_rutri(jc,11,jb)*p_coeff(6,jc,jk,jb) &
          &                  - ptr_rutri(jc,12,jb) * p_coeff(7,jc,jk,jb)                   &
          &                  - ptr_rutri(jc,13,jb) * p_coeff(8,jc,jk,jb)                   &
          &                  - ptr_rutri(jc,14,jb) * p_coeff(9,jc,jk,jb)                   &
          &                  - ptr_rutri(jc,15,jb) * p_coeff(10,jc,jk,jb) )
        p_coeff(4,jc,jk,jb)  = ptr_rrdiag(jc,3,jb)                                         &
          &                  * ( z_qt_times_d(3) - ptr_rutri(jc,16,jb)*p_coeff(5,jc,jk,jb) &
          &                  - ptr_rutri(jc,17,jb) * p_coeff(6,jc,jk,jb)                   &
          &                  - ptr_rutri(jc,18,jb) * p_coeff(7,jc,jk,jb)                   &
          &                  - ptr_rutri(jc,19,jb) * p_coeff(8,jc,jk,jb)                   &
          &                  - ptr_rutri(jc,20,jb) * p_coeff(9,jc,jk,jb)                   &
          &                  - ptr_rutri(jc,21,jb) * p_coeff(10,jc,jk,jb) )
        p_coeff(3,jc,jk,jb)  = ptr_rrdiag(jc,2,jb)                                         &
          &                  * ( z_qt_times_d(2) - ptr_rutri(jc,22,jb)*p_coeff(4,jc,jk,jb) &
          &                  - ptr_rutri(jc,23,jb) * p_coeff(5,jc,jk,jb)                   &
          &                  - ptr_rutri(jc,24,jb) * p_coeff(6,jc,jk,jb)                   &
          &                  - ptr_rutri(jc,25,jb) * p_coeff(7,jc,jk,jb)                   &
          &                  - ptr_rutri(jc,26,jb) * p_coeff(8,jc,jk,jb)                   &
          &                  - ptr_rutri(jc,27,jb) * p_coeff(9,jc,jk,jb)                   &
          &                  - ptr_rutri(jc,28,jb) * p_coeff(10,jc,jk,jb) )
        p_coeff(2,jc,jk,jb)  = ptr_rrdiag(jc,1,jb)                                         &
          &                  * ( z_qt_times_d(1) - ptr_rutri(jc,29,jb)*p_coeff(3,jc,jk,jb) &
          &                  - ptr_rutri(jc,30,jb) * p_coeff(4,jc,jk,jb)                   &
          &                  - ptr_rutri(jc,31,jb) * p_coeff(5,jc,jk,jb)                   &
          &                  - ptr_rutri(jc,32,jb) * p_coeff(6,jc,jk,jb)                   &
          &                  - ptr_rutri(jc,33,jb) * p_coeff(7,jc,jk,jb)                   &
          &                  - ptr_rutri(jc,34,jb) * p_coeff(8,jc,jk,jb)                   &
          &                  - ptr_rutri(jc,35,jb) * p_coeff(9,jc,jk,jb)                   &
          &                  - ptr_rutri(jc,36,jb) * p_coeff(10,jc,jk,jb) )


        p_coeff(1,jc,jk,jb)  = p_cc(jc,jk,jb)                                              &
          &                  - p_coeff(2,jc,jk,jb)  * ptr_int_lsq%lsq_moments(jc,jb,1)     &
          &                  - p_coeff(3,jc,jk,jb)  * ptr_int_lsq%lsq_moments(jc,jb,2)     &
          &                  - p_coeff(4,jc,jk,jb)  * ptr_int_lsq%lsq_moments(jc,jb,3)     &
          &                  - p_coeff(5,jc,jk,jb)  * ptr_int_lsq%lsq_moments(jc,jb,4)     &
          &                  - p_coeff(6,jc,jk,jb)  * ptr_int_lsq%lsq_moments(jc,jb,5)     &
          &                  - p_coeff(7,jc,jk,jb)  * ptr_int_lsq%lsq_moments(jc,jb,6)     &
          &                  - p_coeff(8,jc,jk,jb)  * ptr_int_lsq%lsq_moments(jc,jb,7)     &
          &                  - p_coeff(9,jc,jk,jb)  * ptr_int_lsq%lsq_moments(jc,jb,8)     &
          &                  - p_coeff(10,jc,jk,jb) * ptr_int_lsq%lsq_moments(jc,jb,9)

      END DO ! end loop over cells

    END DO ! end loop over vertical levels

  END DO ! end loop over blocks
#ifdef _OPENACC
!$ACC END PARALLEL
!ACC_DEBUG UPDATE HOST(p_coeff) IF( i_am_accel_node .AND. acc_on )
!$ACC END DATA
#else
!$OMP END DO NOWAIT
!$OMP END PARALLEL
#endif


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

#ifndef __SX__
  REAL(wp)  ::           &        !< difference of scalars i j
    &  z_b(lsq_high_set%dim_c,nproma,ptr_patch%nlev)
#else
  REAL(wp)  ::           &
    &  z_b1, z_b2, z_b3, z_b4, z_b5, z_b6, z_b7, z_b8, z_b9, &
    &  zdp(nproma,ptr_patch%nlev)
  INTEGER :: jj
#endif

  INTEGER, POINTER  ::   &        !< Pointer to line and block indices of
    &  iidx(:,:,:), iblk(:,:,:)   !< required stencil
  INTEGER :: slev, elev           !< vertical start and end level
  INTEGER :: jc, jk, jb           !< index of cell, vertical level and block
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
    elev = UBOUND(p_cc,2)
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

#ifdef _OPENACC
!$ACC DATA PCOPYIN( ptr_patch, ptr_int_lsq, p_cc ), PCOPY( p_coeff ), CREATE( z_b), IF( i_am_accel_node .AND. acc_on )
!ACC_DEBUG UPDATE DEVICE ( p_cc, p_coeff  ), IF( i_am_accel_node .AND. acc_on )
#else
!$OMP PARALLEL
#endif

  IF (ptr_patch%id > 1 .OR. l_limited_area) THEN
    CALL init(p_coeff(:,:,:,1:i_startblk))
#ifndef _OPENACC
!$OMP BARRIER
#endif
  ENDIF

#ifdef _OPENACC
!$ACC PARALLEL &
!$ACC PRESENT( ptr_patch, ptr_int_lsq, p_cc, p_coeff ), &
!$ACC PRIVATE( z_b ), &
!$ACC IF( i_am_accel_node .AND. acc_on )

!$ACC LOOP GANG PRIVATE(i_startidx, i_endidx)
#else
#ifndef __SX__
!$OMP DO PRIVATE(jb,jc,jk,i_startidx,i_endidx,z_b), ICON_OMP_RUNTIME_SCHEDULE
#else
!$OMP DO PRIVATE(jb,jc,jk,i_startidx,i_endidx,z_b1,z_b2,z_b3,z_b4,z_b5,z_b6,z_b7,z_b8,z_b9, &
!$OMP            zdp), ICON_OMP_RUNTIME_SCHEDULE
#endif
#endif
  DO jb = i_startblk, i_endblk

    CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, rl_start, rl_end)

    !
    ! 1. compute right hand side of linear system
    !

#ifndef __SX__
!$ACC LOOP VECTOR COLLAPSE(2)
#ifdef __LOOP_EXCHANGE
    DO jc = i_startidx, i_endidx
      DO jk = slev, elev
#else
    DO jk = slev, elev
      DO jc = i_startidx, i_endidx
#endif

        z_b(1,jc,jk) = p_cc(iidx(jc,jb,1),jk,iblk(jc,jb,1)) - p_cc(jc,jk,jb)
        z_b(2,jc,jk) = p_cc(iidx(jc,jb,2),jk,iblk(jc,jb,2)) - p_cc(jc,jk,jb)
        z_b(3,jc,jk) = p_cc(iidx(jc,jb,3),jk,iblk(jc,jb,3)) - p_cc(jc,jk,jb)
        z_b(4,jc,jk) = p_cc(iidx(jc,jb,4),jk,iblk(jc,jb,4)) - p_cc(jc,jk,jb)
        z_b(5,jc,jk) = p_cc(iidx(jc,jb,5),jk,iblk(jc,jb,5)) - p_cc(jc,jk,jb)
        z_b(6,jc,jk) = p_cc(iidx(jc,jb,6),jk,iblk(jc,jb,6)) - p_cc(jc,jk,jb)
        z_b(7,jc,jk) = p_cc(iidx(jc,jb,7),jk,iblk(jc,jb,7)) - p_cc(jc,jk,jb)
        z_b(8,jc,jk) = p_cc(iidx(jc,jb,8),jk,iblk(jc,jb,8)) - p_cc(jc,jk,jb)
        z_b(9,jc,jk) = p_cc(iidx(jc,jb,9),jk,iblk(jc,jb,9)) - p_cc(jc,jk,jb)

      ENDDO
    ENDDO



    !
    ! 2. compute cell based coefficients for cubic reconstruction
    !    calculate matrix vector product PINV(A) * b
    !
!$ACC LOOP VECTOR COLLAPSE(2)
    DO jk = slev, elev
      DO jc = i_startidx, i_endidx

        ! (intrinsic function matmul not applied, due to massive
        ! performance penalty on the NEC. Instead the intrinsic dot product
        ! function is applied
!CDIR BEGIN EXPAND=9

        p_coeff(10,jc,jk,jb) = DOT_PRODUCT(ptr_int_lsq%lsq_pseudoinv(jc,9,1:9,jb), &
          &                               z_b(1:9,jc,jk))
        p_coeff(9, jc,jk,jb) = DOT_PRODUCT(ptr_int_lsq%lsq_pseudoinv(jc,8,1:9,jb), &
          &                               z_b(1:9,jc,jk))
        p_coeff(8, jc,jk,jb) = DOT_PRODUCT(ptr_int_lsq%lsq_pseudoinv(jc,7,1:9,jb), &
          &                               z_b(1:9,jc,jk))
        p_coeff(7, jc,jk,jb) = DOT_PRODUCT(ptr_int_lsq%lsq_pseudoinv(jc,6,1:9,jb), &
          &                               z_b(1:9,jc,jk))
        p_coeff(6, jc,jk,jb) = DOT_PRODUCT(ptr_int_lsq%lsq_pseudoinv(jc,5,1:9,jb), &
          &                               z_b(1:9,jc,jk))
        p_coeff(5, jc,jk,jb) = DOT_PRODUCT(ptr_int_lsq%lsq_pseudoinv(jc,4,1:9,jb), &
          &                               z_b(1:9,jc,jk))
        p_coeff(4, jc,jk,jb) = DOT_PRODUCT(ptr_int_lsq%lsq_pseudoinv(jc,3,1:9,jb), &
          &                               z_b(1:9,jc,jk))
        p_coeff(3, jc,jk,jb) = DOT_PRODUCT(ptr_int_lsq%lsq_pseudoinv(jc,2,1:9,jb), &
          &                               z_b(1:9,jc,jk))
        p_coeff(2, jc,jk,jb) = DOT_PRODUCT(ptr_int_lsq%lsq_pseudoinv(jc,1,1:9,jb), &
          &                               z_b(1:9,jc,jk))

!CDIR END


!CDIR BEGIN EXPAND=9
        p_coeff(1,jc,jk,jb)  = p_cc(jc,jk,jb) - DOT_PRODUCT(p_coeff(2:10,jc,jk,jb), &
          &                    ptr_int_lsq%lsq_moments(jc,jb,1:9))


      END DO ! end loop over cells

    END DO ! end loop over vertical levels

#else

!$ACC LOOP VECTOR COLLAPSE(3)
    DO jj = 2, 10
      DO jk = slev, elev
!CDIR NODEP
        DO jc = i_startidx, i_endidx

          z_b1 = p_cc(iidx(jc,jb,1),jk,iblk(jc,jb,1)) - p_cc(jc,jk,jb)
          z_b2 = p_cc(iidx(jc,jb,2),jk,iblk(jc,jb,2)) - p_cc(jc,jk,jb)
          z_b3 = p_cc(iidx(jc,jb,3),jk,iblk(jc,jb,3)) - p_cc(jc,jk,jb)
          z_b4 = p_cc(iidx(jc,jb,4),jk,iblk(jc,jb,4)) - p_cc(jc,jk,jb)
          z_b5 = p_cc(iidx(jc,jb,5),jk,iblk(jc,jb,5)) - p_cc(jc,jk,jb)
          z_b6 = p_cc(iidx(jc,jb,6),jk,iblk(jc,jb,6)) - p_cc(jc,jk,jb)
          z_b7 = p_cc(iidx(jc,jb,7),jk,iblk(jc,jb,7)) - p_cc(jc,jk,jb)
          z_b8 = p_cc(iidx(jc,jb,8),jk,iblk(jc,jb,8)) - p_cc(jc,jk,jb)
          z_b9 = p_cc(iidx(jc,jb,9),jk,iblk(jc,jb,9)) - p_cc(jc,jk,jb)

          !
          ! 2. compute cell based coefficients for cubic reconstruction
          !    calculate matrix vector product PINV(A) * b
          !
          ! (intrinsic function matmul not applied, due to massive
          ! performance penalty on the NEC. Instead the intrinsic dot product
          ! function is applied

          p_coeff(jj,jc,jk,jb) = ptr_int_lsq%lsq_pseudoinv(jc,jj-1,1,jb) * z_b1 &
                               + ptr_int_lsq%lsq_pseudoinv(jc,jj-1,2,jb) * z_b2 &
                               + ptr_int_lsq%lsq_pseudoinv(jc,jj-1,3,jb) * z_b3 &
                               + ptr_int_lsq%lsq_pseudoinv(jc,jj-1,4,jb) * z_b4 &
                               + ptr_int_lsq%lsq_pseudoinv(jc,jj-1,5,jb) * z_b5 &
                               + ptr_int_lsq%lsq_pseudoinv(jc,jj-1,6,jb) * z_b6 &
                               + ptr_int_lsq%lsq_pseudoinv(jc,jj-1,7,jb) * z_b7 &
                               + ptr_int_lsq%lsq_pseudoinv(jc,jj-1,8,jb) * z_b8 &
                               + ptr_int_lsq%lsq_pseudoinv(jc,jj-1,9,jb) * z_b9

        END DO ! end loop over cells
      END DO ! end loop over vertical levels
    END DO ! jj

    zdp=0._wp

!$ACC LOOP VECTOR COLLAPSE(3)
    DO jj = 2, 10
      DO jk = slev, elev
        DO jc = i_startidx, i_endidx
          zdp(jc,jk) = zdp(jc,jk) + &
                       p_coeff(jj,jc,jk,jb) * ptr_int_lsq%lsq_moments(jc,jb,jj-1)
        END DO ! end loop over cells
      END DO ! end loop over vertical levels
    END DO ! jj

!$ACC LOOP VECTOR COLLAPSE(2)
    DO jk = slev, elev
      DO jc = i_startidx, i_endidx
        p_coeff(1,jc,jk,jb)  = p_cc(jc,jk,jb) - zdp(jc,jk)
      END DO ! end loop over cells
    END DO ! end loop over vertical levels
#endif

  END DO ! end loop over blocks
#ifdef _OPENACC
!$ACC END PARALLEL
!ACC_DEBUG UPDATE HOST(p_coeff) IF( i_am_accel_node .AND. acc_on )
!$ACC END DATA
#else
!$OMP END DO NOWAIT
!$OMP END PARALLEL
#endif


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
!! Modification by Will Sawyer, CSCS (2014-07-18) :
!! - Split out 2 field version into separate subroutine
!!
SUBROUTINE div3d( vec_e, ptr_patch, ptr_int, div_vec_c, &
  &               opt_slev, opt_elev, opt_rlstart, opt_rlend )
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

INTEGER :: slev, elev     ! vertical start and end level
INTEGER :: jc, jk, jb
INTEGER :: rl_start, rl_end
INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx, i_nchdom

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
  elev = UBOUND(vec_e,2)
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


#ifdef _OPENACC
!$ACC DATA PCOPYIN( ptr_patch, ptr_int, vec_e ), PCOPY( div_vec_c ), IF( i_am_accel_node .AND. acc_on )
!ACC_DEBUG UPDATE DEVICE ( vec_e, div_vec_c ) IF( i_am_accel_node .AND. acc_on )
!$ACC PARALLEL &
!$ACC PRESENT( ptr_patch, ptr_int, vec_e, iidx, iblk, div_vec_c ), &
!$ACC IF( i_am_accel_node .AND. acc_on )

!$ACC LOOP GANG PRIVATE(i_startidx, i_endidx)
#else
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk) ICON_OMP_DEFAULT_SCHEDULE
#endif
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

!$ACC LOOP VECTOR COLLAPSE(2)
#ifdef __LOOP_EXCHANGE
      DO jc = i_startidx, i_endidx
        DO jk = slev, elev
#else
      DO jk = slev, elev
        DO jc = i_startidx, i_endidx
#endif

          div_vec_c(jc,jk,jb) =  &
            vec_e(iidx(jc,jb,1),jk,iblk(jc,jb,1)) * ptr_int%geofac_div(jc,1,jb) + &
            vec_e(iidx(jc,jb,2),jk,iblk(jc,jb,2)) * ptr_int%geofac_div(jc,2,jb) + &
            vec_e(iidx(jc,jb,3),jk,iblk(jc,jb,3)) * ptr_int%geofac_div(jc,3,jb)

        END DO
      END DO

  END DO
#ifdef _OPENACC
!$ACC END PARALLEL
!ACC_DEBUG UPDATE HOST(div_vec_c) IF( i_am_accel_node .AND. acc_on )
!$ACC END DATA
#else
!$OMP END DO NOWAIT
!$OMP END PARALLEL
#endif

!IF(ltimer) CALL timer_stop(timer_div)

END SUBROUTINE div3d

SUBROUTINE div3d_2field( vec_e, ptr_patch, ptr_int, div_vec_c, &
  &                      opt_slev, opt_elev, in2, out2, opt_rlstart, opt_rlend )
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

! second input field for more efficient processing in NH core
REAL(wp), INTENT(in) ::  &
  &  in2(:,:,:) ! dim: (nproma,nlev,nblks_e)

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

! second output field
REAL(wp), OPTIONAL, INTENT(inout) ::  &
  &  out2(:,:,:) ! dim: (nproma,nlev,nblks_c)

INTEGER :: slev, elev     ! vertical start and end level
INTEGER :: jc, jk, jb
INTEGER :: rl_start, rl_end
INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx, i_nchdom

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
  elev = UBOUND(vec_e,2)
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


#ifdef _OPENACC
!$ACC DATA PCOPYIN( ptr_patch, ptr_int, vec_e, in2 ), PCOPY( div_vec_c, out2 ), IF( i_am_accel_node .AND. acc_on )
!ACC_DEBUG UPDATE DEVICE ( vec_e, in2, div_vec_c, out2 ), IF( i_am_accel_node .AND. acc_on )
!$ACC PARALLEL &
!$ACC PRESENT( ptr_patch, ptr_int, vec_e, in2, div_vec_c, out2 ), &
!$ACC IF( i_am_accel_node .AND. acc_on )

!$ACC LOOP GANG PRIVATE(i_startidx, i_endidx)
#else
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk) ICON_OMP_DEFAULT_SCHEDULE
#endif
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

!$ACC LOOP VECTOR COLLAPSE(2)
#ifdef __LOOP_EXCHANGE
      DO jc = i_startidx, i_endidx
        DO jk = slev, elev
#else
      DO jk = slev, elev
        DO jc = i_startidx, i_endidx
#endif

          div_vec_c(jc,jk,jb) =  &
            vec_e(iidx(jc,jb,1),jk,iblk(jc,jb,1)) * ptr_int%geofac_div(jc,1,jb) + &
            vec_e(iidx(jc,jb,2),jk,iblk(jc,jb,2)) * ptr_int%geofac_div(jc,2,jb) + &
            vec_e(iidx(jc,jb,3),jk,iblk(jc,jb,3)) * ptr_int%geofac_div(jc,3,jb)

          out2(jc,jk,jb) =  &
            in2(iidx(jc,jb,1),jk,iblk(jc,jb,1)) * ptr_int%geofac_div(jc,1,jb) + &
            in2(iidx(jc,jb,2),jk,iblk(jc,jb,2)) * ptr_int%geofac_div(jc,2,jb) + &
            in2(iidx(jc,jb,3),jk,iblk(jc,jb,3)) * ptr_int%geofac_div(jc,3,jb)

        END DO
      END DO

  END DO
#ifdef _OPENACC
!$ACC END PARALLEL
!ACC_DEBUG UPDATE HOST(div_vec_c, out2) IF( i_am_accel_node .AND. acc_on )
!$ACC END DATA
#else
!$OMP END DO NOWAIT
!$OMP END PARALLEL
#endif

!IF(ltimer) CALL timer_stop(timer_div)

END SUBROUTINE div3d_2field

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
REAL(vp), INTENT(inout) ::  &
  &  f4dout(:,:,:,:) ! dim: (nproma,nlev,nblks_c,dim4d)

INTEGER :: slev(dim4d), elev(dim4d)     ! vertical start and end level
INTEGER :: jc, jk, jb, ji
INTEGER :: rl_start, rl_end
INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx, i_nchdom

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
  elev = UBOUND(f4din,2)
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


#ifdef _OPENACC
!$ACC DATA PCOPYIN( ptr_patch, ptr_int, f4din ), PCOPY( f4dout ),  IF( i_am_accel_node .AND. acc_on )
!ACC_DEBUG UPDATE DEVICE ( f4din, f4dout ), IF( i_am_accel_node .AND. acc_on )
!$ACC PARALLEL &
!$ACC PRESENT( ptr_patch, ptr_int, f4din, f4dout ), &
!$ACC IF( i_am_accel_node .AND. acc_on )

!$ACC LOOP GANG PRIVATE(i_startidx, i_endidx)
#else
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk,ji) ICON_OMP_DEFAULT_SCHEDULE
#endif
  DO jb = i_startblk, i_endblk

    CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, rl_start, rl_end)

!$ACC LOOP VECTOR
#ifdef __LOOP_EXCHANGE
    DO jc = i_startidx, i_endidx
      DO ji = 1, dim4d
        DO jk = slev(ji), elev(ji)
#else
    DO ji = 1, dim4d
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

#ifdef _OPENACC
!$ACC END PARALLEL
!ACC_DEBUG UPDATE HOST(f4dout) IF( i_am_accel_node .AND. acc_on )
!$ACC END DATA
#else
!$OMP END DO NOWAIT
!$OMP END PARALLEL
#endif

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
SUBROUTINE div_avg( vec_e, ptr_patch, ptr_int, avg_coeff, div_vec_c,    &
  &                 opt_in2, opt_out2, opt_slev, opt_elev, opt_rlstart, &
  &                 opt_rlend )
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
  &  opt_slev    ! optional vertical start level

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_elev    ! optional vertical end level

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
INTEGER :: slev, elev     ! vertical start and end level
INTEGER :: rl_start, rl_end, rl_start_l2, rl_end_l1
INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx, i_nchdom

REAL(wp), DIMENSION (nproma,ptr_patch%nlev,ptr_patch%nblks_c) :: aux_c, aux_c2

INTEGER,  DIMENSION(:,:,:),   POINTER :: inidx, inblk, ieidx, ieblk
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
  elev = UBOUND(vec_e,2)
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


! values for the blocking
i_nchdom   = MAX(1,ptr_patch%n_childdom)

! First compute divergence
!
#ifdef _OPENACC
!$ACC DATA PCOPYIN( ptr_patch, ptr_int, vec_e, avg_coeff ), PCOPY( div_vec_c ), CREATE( aux_c ), IF( i_am_accel_node .AND. acc_on )
!$ACC DATA PCOPYIN( opt_in2 ), PCOPY( opt_out2 ), CREATE( aux_c2 ), IF( i_am_accel_node .AND. acc_on .AND. l2fields )
!ACC_DEBUG UPDATE DEVICE ( vec_e, div_vec_c ), IF( i_am_accel_node .AND. acc_on )
!ACC_DEBUG UPDATE DEVICE ( opt_in2, opt_out2 ), IF( i_am_accel_node .AND. acc_on .AND. l2fields )
#else
!$OMP PARALLEL PRIVATE(i_startblk,i_endblk)
#endif
i_startblk = ptr_patch%cells%start_blk(rl_start,1)
i_endblk   = ptr_patch%cells%end_blk(rl_end_l1,i_nchdom)

IF (l2fields) THEN

#ifdef _OPENACC
!$ACC PARALLEL &
!$ACC PRESENT( ptr_patch, ptr_int, vec_e, opt_in2, aux_c, aux_c2 ), &
!$ACC IF( i_am_accel_node .AND. acc_on )

!$ACC LOOP GANG PRIVATE(i_startidx, i_endidx)
#else
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk), ICON_OMP_RUNTIME_SCHEDULE
#endif
  DO jb = i_startblk, i_endblk

    CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
                     i_startidx, i_endidx, rl_start, rl_end_l1)

!$ACC LOOP VECTOR COLLAPSE(2)
#ifdef __LOOP_EXCHANGE
    DO jc = i_startidx, i_endidx
      DO jk = slev, elev
#else
    DO jk = slev, elev
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
  END DO
#ifdef _OPENACC
!$ACC END PARALLEL
#else
!$OMP END DO
#endif

ELSE

#ifdef _OPENACC
!$ACC PARALLEL &
!$ACC PRESENT( ptr_patch, ptr_int, vec_e, aux_c ), &
!$ACC IF( i_am_accel_node .AND. acc_on )

!$ACC LOOP GANG PRIVATE(i_startidx, i_endidx)
#else
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk), ICON_OMP_RUNTIME_SCHEDULE
#endif
  DO jb = i_startblk, i_endblk

    CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
                     i_startidx, i_endidx, rl_start, rl_end_l1)

!$ACC LOOP VECTOR COLLAPSE(2)
#ifdef __LOOP_EXCHANGE
    DO jc = i_startidx, i_endidx
      DO jk = slev, elev
#else
    DO jk = slev, elev
      DO jc = i_startidx, i_endidx
#endif
        aux_c(jc,jk,jb) =  &
          vec_e(ieidx(jc,jb,1),jk,ieblk(jc,jb,1)) * ptr_int%geofac_div(jc,1,jb) + &
          vec_e(ieidx(jc,jb,2),jk,ieblk(jc,jb,2)) * ptr_int%geofac_div(jc,2,jb) + &
          vec_e(ieidx(jc,jb,3),jk,ieblk(jc,jb,3)) * ptr_int%geofac_div(jc,3,jb)

      END DO
    END DO

  END DO
#ifdef _OPENACC
!$ACC END PARALLEL
#else
!$OMP END DO
#endif

ENDIF

IF (l_limited_area .OR. ptr_patch%id > 1) THEN
  ! Fill div_vec_c along the lateral boundaries of nests

  i_startblk = ptr_patch%cells%start_blk(rl_start,1)
  i_endblk   = ptr_patch%cells%end_blk(rl_start_l2,1)
!
  CALL copy(aux_c (:,:,i_startblk:i_endblk), &
       div_vec_c(:,:,i_startblk:i_endblk))
  IF (l2fields) &
       CALL copy(aux_c2(:,:,i_startblk:i_endblk), &
       &         opt_out2 (:,:,i_startblk:i_endblk))
#ifndef _OPENACC
!$OMP BARRIER
#endif
ENDIF

!
! Now do averaging with weights given by avg_coeff

! values for the blocking
i_startblk = ptr_patch%cells%start_blk(rl_start_l2,1)
i_endblk   = ptr_patch%cells%end_blk(rl_end,i_nchdom)
!
! loop through all patch cells (and blocks)
!

IF (l2fields) THEN

#ifdef _OPENACC
!$ACC PARALLEL &
!$ACC PRESENT( ptr_patch, aux_c, aux_c2, avg_coeff, div_vec_c, opt_out2 ), &
!$ACC IF( i_am_accel_node .AND. acc_on )

!$ACC LOOP GANG PRIVATE(i_startidx, i_endidx)
#else
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk), ICON_OMP_RUNTIME_SCHEDULE
#endif
  DO jb = i_startblk, i_endblk

    CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
                     i_startidx, i_endidx, rl_start_l2, rl_end)

!$ACC LOOP VECTOR COLLAPSE(2)
#ifdef __LOOP_EXCHANGE
    DO jc = i_startidx, i_endidx
      DO jk = slev, elev
#else
    DO jk = slev, elev
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
  END DO !block loop

#ifdef _OPENACC
!$ACC END PARALLEL
#else
!$OMP END DO NOWAIT
#endif

ELSE

#ifdef _OPENACC
!$ACC PARALLEL &
!$ACC PRESENT( ptr_patch, aux_c, avg_coeff, div_vec_c), &
!$ACC IF( i_am_accel_node .AND. acc_on )

!$ACC LOOP GANG PRIVATE(i_startidx, i_endidx)
#else
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk), ICON_OMP_RUNTIME_SCHEDULE
#endif
  DO jb = i_startblk, i_endblk

    CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
                     i_startidx, i_endidx, rl_start_l2, rl_end)

!$ACC LOOP VECTOR COLLAPSE(2)
#ifdef __LOOP_EXCHANGE
    DO jc = i_startidx, i_endidx
      DO jk = slev, elev
#else
    DO jk = slev, elev
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

  END DO !block loop

#ifdef _OPENACC
!$ACC END PARALLEL
#else
!$OMP END DO NOWAIT
#endif

ENDIF

#ifdef _OPENACC
!ACC_DEBUG UPDATE HOST(div_vec_c) IF( i_am_accel_node .AND. acc_on )
!ACC_DEBUG UPDATE HOST(opt_out2) IF( i_am_accel_node .AND. acc_on .AND. l2fields )
!$ACC END DATA
!$ACC END DATA
#else
!$OMP END PARALLEL
#endif

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
  elev = UBOUND(vec_e,2)
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
#ifdef _OPENACC
!$ACC DATA PCOPYIN( ptr_patch, ptr_int, vec_e ), PCOPY( div_vec_e ), IF( i_am_accel_node .AND. acc_on )
!ACC_DEBUG UPDATE DEVICE ( vec_e, div_vec_e ) IF( i_am_accel_node .AND. acc_on )
!$ACC PARALLEL &
!$ACC PRESENT( ptr_patch, ptr_int, vec_e, div_vec_e ), &
!$ACC IF( i_am_accel_node .AND. acc_on )

!$ACC LOOP GANG PRIVATE(i_startidx, i_endidx)
#else
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,je,jk,i_startidx,i_endidx) ICON_OMP_RUNTIME_SCHEDULE
#endif
DO jb = i_startblk, i_endblk

  CALL get_indices_e(ptr_patch, jb, i_startblk, i_endblk, &
                     i_startidx, i_endidx, rl_start, rl_end)

!$ACC LOOP VECTOR COLLAPSE(2)
#ifdef __LOOP_EXCHANGE
  DO je = i_startidx, i_endidx
    DO jk = slev, elev
#else
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
#ifdef _OPENACC
!$ACC END PARALLEL
!ACC_DEBUG UPDATE HOST(div_vec_e) IF( i_am_accel_node .AND. acc_on )
!$ACC END DATA
#else
!$OMP END DO NOWAIT
!$OMP END PARALLEL
#endif

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
  elev = UBOUND(vec_e,2)
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

  iidx => ptr_patch%verts%edge_idx
  iblk => ptr_patch%verts%edge_blk

  ! values for the blocking
  i_nchdom   = MAX(1,ptr_patch%n_childdom)
  i_startblk = ptr_patch%verts%start_blk(rl_start,1)
  i_endblk   = ptr_patch%verts%end_blk(rl_end,i_nchdom)

#ifdef _OPENACC
!$ACC DATA PCOPYIN( ptr_patch, ptr_int, vec_e ), PCOPY( rot_vec ), IF( i_am_accel_node .AND. acc_on )
!ACC_DEBUG UPDATE DEVICE ( vec_e, rot_vec ) IF( i_am_accel_node .AND. acc_on )
!$ACC PARALLEL &
!$ACC PRESENT( ptr_patch, ptr_int, vec_e, rot_vec ), &
!$ACC IF( i_am_accel_node .AND. acc_on )
!$ACC LOOP GANG PRIVATE(i_startidx, i_endidx)
#else
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jv,jk), ICON_OMP_RUNTIME_SCHEDULE
#endif
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

!$ACC LOOP VECTOR COLLAPSE(2)
#ifdef __LOOP_EXCHANGE
    DO jv = i_startidx, i_endidx
      DO jk = slev, elev
#else
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
#ifdef _OPENACC
!$ACC END PARALLEL
!ACC_DEBUG UPDATE HOST(rot_vec) IF( i_am_accel_node .AND. acc_on )
!$ACC END DATA
#else
!$OMP END DO NOWAIT
!$OMP END PARALLEL
#endif


END SUBROUTINE rot_vertex_atmos

!>
!! Same as above routine, but expects reversed index order (vertical first)
!! of the output field if __LOOP_EXCHANGE is specified. In addition, the 
!! output field (vorticity) has single precision if __MIXED_PRECISION is specified
!!
!!
SUBROUTINE rot_vertex_ri( vec_e, ptr_patch, ptr_int, rot_vec, &
  &                       opt_slev, opt_elev, opt_rlend )
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
  &  opt_rlend   ! end value of refin_ctrl flag

!
!  vertex based variable in which rotation is stored
!
REAL(vp), INTENT(inout) ::  &
  &  rot_vec(:,:,:) ! dim: (nproma,nlev,nblks_v) or (nlev,nproma,nblks_v)

INTEGER :: slev, elev     ! vertical start and end level
INTEGER :: jv, jk, jb

INTEGER :: rl_start, rl_end
INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx

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
  elev = UBOUND(vec_e,2)
END IF

rl_start = 2

IF ( PRESENT(opt_rlend) ) THEN
  rl_end = opt_rlend
ELSE
  rl_end = min_rlvert
END IF


!
!  loop through over all patch vertices (and blocks)
!

  iidx => ptr_patch%verts%edge_idx
  iblk => ptr_patch%verts%edge_blk

  ! values for the blocking
  i_startblk = ptr_patch%verts%start_block(rl_start)
  i_endblk   = ptr_patch%verts%end_block(rl_end)


#ifdef _OPENACC
!$ACC DATA PCOPYIN( ptr_patch, ptr_int, vec_e ), PCOPY( rot_vec ), IF( i_am_accel_node .AND. acc_on )
!ACC_DEBUG UPDATE DEVICE ( vec_e, rot_vec ), IF( i_am_accel_node .AND. acc_on )
!$ACC PARALLEL &
!$ACC PRESENT( ptr_patch, ptr_int, vec_e, rot_vec ), &
!$ACC IF( i_am_accel_node .AND. acc_on )
!$ACC LOOP GANG PRIVATE(i_startidx, i_endidx)
#else
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jv,jk), ICON_OMP_RUNTIME_SCHEDULE
#endif
  DO jb = i_startblk, i_endblk

    CALL get_indices_v(ptr_patch, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, rl_start, rl_end)

    ! calculate rotation, i.e.
    ! add individual edge contributions to rotation
    !
!$ACC LOOP VECTOR COLLAPSE(2)
#ifdef __LOOP_EXCHANGE
    DO jv = i_startidx, i_endidx
      DO jk = slev, elev
        rot_vec(jk,jv,jb) =   &
#else
    DO jk = slev, elev
      DO jv = i_startidx, i_endidx
        rot_vec(jv,jk,jb) =   &
#endif
          vec_e(iidx(jv,jb,1),jk,iblk(jv,jb,1)) * ptr_int%geofac_rot(jv,1,jb) + &
          vec_e(iidx(jv,jb,2),jk,iblk(jv,jb,2)) * ptr_int%geofac_rot(jv,2,jb) + &
          vec_e(iidx(jv,jb,3),jk,iblk(jv,jb,3)) * ptr_int%geofac_rot(jv,3,jb) + &
          vec_e(iidx(jv,jb,4),jk,iblk(jv,jb,4)) * ptr_int%geofac_rot(jv,4,jb) + &
          vec_e(iidx(jv,jb,5),jk,iblk(jv,jb,5)) * ptr_int%geofac_rot(jv,5,jb) + &
          vec_e(iidx(jv,jb,6),jk,iblk(jv,jb,6)) * ptr_int%geofac_rot(jv,6,jb)

      END DO
    END DO

  ENDDO
#ifdef _OPENACC
!$ACC END PARALLEL
!ACC_DEBUG UPDATE HOST(rot_vec), IF( i_am_accel_node .AND. acc_on )
!$ACC END DATA
#else
!$OMP END DO NOWAIT
!$OMP END PARALLEL
#endif

END SUBROUTINE rot_vertex_ri

END MODULE mo_math_divrot
