!>
!!   Contains the implementation of the mathematical grad operators.
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
!!  Modification by William Sawyer, CSCS (2014-09-23)
!!  - OpenACC implementation
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

MODULE mo_math_gradients
!-------------------------------------------------------------------------
!
!    ProTeX FORTRAN source: Style 2
!    modified for ICON project, DWD/MPI-M 2006
!
!-------------------------------------------------------------------------
!
!
!
USE mo_kind,               ONLY: wp, vp
USE mo_impl_constants,     ONLY: min_rlcell, min_rledge
USE mo_intp_data_strc,     ONLY: t_int_state
USE mo_intp,               ONLY: cells2edges_scalar
USE mo_model_domain,       ONLY: t_patch
USE mo_parallel_config,    ONLY: nproma
USE mo_run_config,         ONLY: timers_level
USE mo_exception,          ONLY: finish
USE mo_timer,              ONLY: timer_start, timer_stop, timer_grad
USE mo_loopindices,        ONLY: get_indices_c, get_indices_e
USE mo_fortran_tools,      ONLY: init
#ifdef _OPENACC
USE mo_mpi,                 ONLY: i_am_accel_node
#endif

IMPLICIT NONE

PRIVATE


PUBLIC :: grad_fd_norm, grad_fd_tang
PUBLIC :: grad_green_gauss_cell
PUBLIC :: grad_fe_cell

INTERFACE grad_green_gauss_cell
  MODULE PROCEDURE grad_green_gauss_cell_adv
  MODULE PROCEDURE grad_green_gauss_cell_dycore
END INTERFACE

INTERFACE grad_fe_cell
  MODULE PROCEDURE grad_fe_cell_adv
  MODULE PROCEDURE grad_fe_cell_adv_2d
  MODULE PROCEDURE grad_fe_cell_dycore
END INTERFACE

#if defined( _OPENACC )
#if defined(__MATH_GRADIENTS_NOACC)
  LOGICAL, PARAMETER ::  acc_on = .FALSE.
#else
  LOGICAL, PARAMETER ::  acc_on = .TRUE.
#endif
  LOGICAL, PARAMETER ::  acc_validate = .FALSE.     !  THIS SHOULD BE .FALSE. AFTER VALIDATION PHASE!
#endif

CONTAINS

!-------------------------------------------------------------------------
!

!-------------------------------------------------------------------------
!
!>
!!  Computes directional  derivative of a cell centered variable.
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
!!
SUBROUTINE grad_fd_norm( psi_c, ptr_patch, grad_norm_psi_e, &
  &                      opt_slev, opt_elev, opt_rlstart, opt_rlend )

!

!
!  patch on which computation is performed
!
TYPE(t_patch), TARGET, INTENT(in) :: ptr_patch
!
!  cell based variable of which normal derivative is computed
!
REAL(wp), INTENT(in) ::  &
  &  psi_c(:,:,:)       ! dim: (nproma,nlev,nblks_c)

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_slev    ! optional vertical start level

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_elev    ! optional vertical end level

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_rlstart, opt_rlend   ! start and end values of refin_ctrl flag

!
!  edge based variable in which normal derivative is stored
!
!REAL(wp), INTENT(out) ::  &
REAL(wp), INTENT(inout) ::  &
  &  grad_norm_psi_e(:,:,:)  ! dim: (nproma,nlev,nblks_e)

INTEGER :: slev, elev     ! vertical start and end level
INTEGER :: je, jk, jb
INTEGER :: rl_start, rl_end
INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx, i_nchdom

INTEGER,  DIMENSION(:,:,:),   POINTER :: iidx, iblk

!
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
  elev = UBOUND(psi_c,2)
END IF

IF ( PRESENT(opt_rlstart) ) THEN
  ! rl_start=1 means edges located along a lateral boundary of a nested
  ! domain. For those, gradient computation is not possible
  IF (opt_rlstart == 1) THEN
    CALL finish ('mo_math_operators:grad_fd_norm',  &
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

iidx => ptr_patch%edges%cell_idx
iblk => ptr_patch%edges%cell_blk

i_nchdom   = MAX(1,ptr_patch%n_childdom)

i_startblk = ptr_patch%edges%start_blk(rl_start,1)
i_endblk   = ptr_patch%edges%end_blk(rl_end,i_nchdom)
!
!  loop through all patch edges (and blocks)
!

IF (timers_level > 5) CALL timer_start(timer_grad)

!$ACC DATA PCOPYIN( psi_c ) PCOPYOUT( grad_norm_psi_e )                    &
!$ACC      PRESENT( ptr_patch, iidx, iblk ) IF( i_am_accel_node .AND. acc_on )
!$ACC UPDATE DEVICE( psi_c ), IF( i_am_accel_node .AND. acc_on .AND. acc_validate )

!$OMP PARALLEL

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,je,jk) ICON_OMP_DEFAULT_SCHEDULE
  DO jb = i_startblk, i_endblk

  CALL get_indices_e(ptr_patch, jb, i_startblk, i_endblk, &
                     i_startidx, i_endidx, rl_start, rl_end)

    !$ACC PARALLEL IF( i_am_accel_node .AND. acc_on )
#ifdef __LOOP_EXCHANGE
    !$ACC LOOP GANG
    DO je = i_startidx, i_endidx
      !$ACC LOOP VECTOR
      DO jk = slev, elev
#else
    !$ACC LOOP GANG
    DO jk = slev, elev
      !$ACC LOOP VECTOR
      DO je = i_startidx, i_endidx
#endif
      !
      ! compute the normal derivative
      ! by the finite difference approximation
      ! (see Bonaventura and Ringler MWR 2005)
      !
       grad_norm_psi_e(je,jk,jb) =  &
          &  ( psi_c(iidx(je,jb,2),jk,iblk(je,jb,2)) - &
          &    psi_c(iidx(je,jb,1),jk,iblk(je,jb,1)) )  &
          &  * ptr_patch%edges%inv_dual_edge_length(je,jb)

      ENDDO

    END DO
    !$ACC END PARALLEL

  END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

!$ACC UPDATE HOST( grad_norm_psi_e ), IF( i_am_accel_node .AND. acc_on .AND. acc_validate )
!$ACC END DATA

IF (timers_level > 5) CALL timer_stop(timer_grad)


END SUBROUTINE grad_fd_norm

!-------------------------------------------------------------------------
!
! RESTRUCT: @Marco: please adjust calls to this routine to your needs.
!>
!! Computes directional derivative of a vertex centered variable with.
!!
!! Computes directional derivative of a vertex centered variable with
!! respect to direction tanget to triangle edge. Notice that the
!! tangential direction is defined by
!!   iorient*(vertex2 - vertex1)
!! input: lives on vertices of triangles
!! output: lives on edges (velocity points)
!!
!! @par Revision History
!! Developed  Marco Restelli (2007-11-23).
!!
SUBROUTINE grad_fd_tang( psi_v, ptr_patch, grad_tang_psi_e,  &
  &                      opt_slev, opt_elev, opt_rlstart, opt_rlend )
!

!
!  patch on which computation is performed
!
TYPE(t_patch), TARGET, INTENT(in) :: ptr_patch
!
! vertex based variable of which tangential derivative is computed
!
REAL(wp), INTENT(in) ::  &
  &  psi_v(:,:,:)

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_slev    ! optional vertical start level

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_elev    ! optional vertical end level

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_rlstart, opt_rlend   ! start and end values of refin_ctrl flag

!
!  edge based variable in which tangential derivative is stored
!
!REAL(wp), INTENT(out) ::  &
REAL(wp), INTENT(inout) ::  &
  &  grad_tang_psi_e(:,:,:)

REAL(wp) :: iorient

INTEGER :: slev, elev     ! vertical start and end level
INTEGER :: je, jk, jb
INTEGER :: rl_start, rl_end
INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx, i_nchdom
INTEGER, DIMENSION(nproma) ::  &
  &  ilv1, ibv1, ilv2, ibv2
!
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
  elev = UBOUND(psi_v,2)
END IF

IF ( PRESENT(opt_rlstart) ) THEN
  ! The possible domain extent depends on the reconstruction algorithm
  ! used to calculate psi_v; the following values are valid for a 6-point
  ! stencil. In the hexagon case (where prognostic variables are located
  ! at vertices), rl_start may be set to 1.
  IF ((opt_rlstart >= 1) .AND. (opt_rlstart <= 2)) THEN
    CALL finish ('mo_math_operators:grad_fd_tang',  &
          &      'opt_rlstart must not be equal to 1 or 2')
  ENDIF
  rl_start = opt_rlstart
ELSE
  rl_start = 3
END IF
IF ( PRESENT(opt_rlend) ) THEN
  rl_end = opt_rlend
ELSE
  rl_end = min_rledge
END IF

! values for the blocking
i_nchdom   = MAX(1,ptr_patch%n_childdom)
i_startblk = ptr_patch%edges%start_blk(rl_start,1)
i_endblk   = ptr_patch%edges%end_blk(rl_end,i_nchdom)

!$ACC DATA PCOPYIN( psi_v ) PCOPYOUT( grad_tang_psi_e ) PRESENT( ptr_patch )   &
!$ACC      CREATE( ilv1, ibv1, ilv2, ibv2 ) IF( i_am_accel_node .AND. acc_on )
!$ACC UPDATE DEVICE( psi_v ), IF( i_am_accel_node .AND. acc_on .AND. acc_validate )

!
! TODO: OpenMP
!

!
!  loop through all patch edges (and blocks)
!
DO jb = i_startblk, i_endblk

  CALL get_indices_e(ptr_patch, jb, i_startblk, i_endblk, &
                     i_startidx, i_endidx, rl_start, rl_end)

  !$ACC PARALLEL IF( i_am_accel_node .AND. acc_on )
  !$ACC LOOP GANG VECTOR
  DO je = i_startidx, i_endidx
    !
    !  get the line and block indices of the vertices of edge je
    !
    ilv1(je) = ptr_patch%edges%vertex_idx(je,jb,1)
    ibv1(je) = ptr_patch%edges%vertex_blk(je,jb,1)
    ilv2(je) = ptr_patch%edges%vertex_idx(je,jb,2)
    ibv2(je) = ptr_patch%edges%vertex_blk(je,jb,2)
  END DO
  !$ACC END PARALLEL

  !$ACC PARALLEL IF( i_am_accel_node .AND. acc_on )
  !$ACC LOOP GANG
  DO jk = slev, elev

    !$ACC LOOP VECTOR PRIVATE( iorient )
    DO je = i_startidx, i_endidx
      !
      ! compute the tangential derivative
      ! by the finite difference approximation
      iorient = ptr_patch%edges%tangent_orientation(je,jb)
      grad_tang_psi_e(je,jk,jb) = iorient  &
        &  * ( psi_v(ilv2(je),jk,ibv2(je)) - psi_v(ilv1(je),jk,ibv1(je)) )  &
        &    / ptr_patch%edges%primal_edge_length(je,jb)
    END DO

  END DO
  !$ACC END PARALLEL

END DO
!
! TODO: OpenMP
!

!$ACC UPDATE HOST( grad_tang_psi_e ), IF( i_am_accel_node .AND. acc_on .AND. acc_validate )
!$ACC END DATA

END SUBROUTINE grad_fd_tang

!-------------------------------------------------------------------------
!
!
!>
!! Computes the cell centered gradient in geographical coordinates.
!!
!! The gradient is computed by taking the derivative of the shape functions
!! for a three-node triangular element (Finite Element thinking).
!!
!! @par Revision History
!!  Initial revision by Daniel Reinert, DWD (2013-11-07)
!!
!! LITERATURE:
!! Fish. J and T. Belytschko, 2007: A first course in finite elements,
!!                                  John Wiley and Sons
!!
!!
SUBROUTINE grad_fe_cell_adv( p_cc, ptr_patch, ptr_int, p_grad, &
  &                      opt_slev, opt_elev, opt_rlstart,  &
  &                      opt_rlend                         )
!
!
!  patch on which computation is performed
!
TYPE(t_patch), TARGET, INTENT(in)     :: ptr_patch
!
!  data structure for interpolation
!
TYPE(t_int_state), TARGET, INTENT(in) :: ptr_int

!
!  cell centered variable
!
REAL(wp), INTENT(in) ::  &
  &  p_cc(:,:,:)

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_slev    ! optional vertical start level

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_elev    ! optional vertical end level

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_rlstart, opt_rlend   ! start and end values of refin_ctrl flag
!
! cell based Green-Gauss reconstructed geographical gradient vector
!
REAL(vp), INTENT(inout) ::  &
  &  p_grad(:,:,:,:)      ! dim:(2,nproma,nlev,nblks_c)

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


iidx => ptr_patch%cells%neighbor_idx
iblk => ptr_patch%cells%neighbor_blk

i_nchdom = MAX(1,ptr_patch%n_childdom)


!
! 2. reconstruction of cell based geographical gradient
!

!$ACC DATA PCOPYIN( p_cc ) PCOPYOUT( p_grad )                                      &
!$ACC      PRESENT( ptr_int, iidx, iblk ) IF( i_am_accel_node .AND. acc_on )
!$ACC UPDATE DEVICE( p_cc ), IF( i_am_accel_node .AND. acc_on .AND. acc_validate )
! Add $ser directives here

!$OMP PARALLEL PRIVATE(i_startblk,i_endblk)

  i_startblk = ptr_patch%cells%start_blk(rl_start,1)
  i_endblk   = ptr_patch%cells%end_blk(rl_end,i_nchdom)

  IF (ptr_patch%id > 1) THEN
  ! Fill nest boundaries with zero to avoid trouble with MPI synchronization

#ifdef _OPENACC
!$ACC KERNELS PRESENT( p_grad ), IF( i_am_accel_node .AND. acc_on )
    p_grad(:,:,:,1:i_startblk) = 0._wp
!$ACC END KERNELS
#else
    CALL init(p_grad(:,:,:,1:i_startblk))
!$OMP BARRIER
#endif
  ENDIF

!$OMP DO PRIVATE(jb,jc,jk,i_startidx,i_endidx), ICON_OMP_RUNTIME_SCHEDULE
  DO jb = i_startblk, i_endblk

    CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, rl_start, rl_end)

    !$ACC PARALLEL IF( i_am_accel_node .AND. acc_on )
#ifdef __LOOP_EXCHANGE
    !$ACC LOOP GANG
    DO jc = i_startidx, i_endidx
!DIR$ IVDEP
      !$ACC LOOP VECTOR
      DO jk = slev, elev
#else
    !$ACC LOOP GANG
    DO jk = slev, elev
      !$ACC LOOP VECTOR
      DO jc = i_startidx, i_endidx
#endif

        ! We do not make use of the intrinsic function DOT_PRODUCT on purpose,
        ! since it is extremely slow on the SX9, when combined with indirect
        ! addressing.

        ! multiply cell-based input values with precomputed grid geometry factor

        ! zonal(u)-component of Green-Gauss gradient
        p_grad(1,jc,jk,jb) = &
          &    ptr_int%gradc_bmat(jc,1,1,jb)*p_cc(iidx(jc,jb,1),jk,iblk(jc,jb,1))  &
          &  + ptr_int%gradc_bmat(jc,1,2,jb)*p_cc(iidx(jc,jb,2),jk,iblk(jc,jb,2))  &
          &  + ptr_int%gradc_bmat(jc,1,3,jb)*p_cc(iidx(jc,jb,3),jk,iblk(jc,jb,3))

        ! meridional(v)-component of Green-Gauss gradient
        p_grad(2,jc,jk,jb) =  &
          &    ptr_int%gradc_bmat(jc,2,1,jb)*p_cc(iidx(jc,jb,1),jk,iblk(jc,jb,1))  &
          &  + ptr_int%gradc_bmat(jc,2,2,jb)*p_cc(iidx(jc,jb,2),jk,iblk(jc,jb,2))  &
          &  + ptr_int%gradc_bmat(jc,2,3,jb)*p_cc(iidx(jc,jb,3),jk,iblk(jc,jb,3))

      END DO ! end loop over cells
    END DO ! end loop over vertical levels
    !$ACC END PARALLEL

  END DO ! end loop over blocks
!$OMP END DO NOWAIT
!$OMP END PARALLEL

!$ACC UPDATE HOST( p_grad ), IF( i_am_accel_node .AND. acc_on .AND. acc_validate )
! Add $ser directives here
!$ACC END DATA


END SUBROUTINE grad_fe_cell_adv




!-------------------------------------------------------------------------
!
!
!>
!! Computes the cell centered gradient in geographical coordinates.
!!
!! The gradient is computed by taking the derivative of the shape functions
!! for a three-node triangular element (Finite Element thinking).
!! 2D version, i.e. for a single vertical level
!!
!! @par Revision History
!!  Initial revision by Daniel Reinert, DWD (2013-11-07)
!!
!! LITERATURE:
!! Fish. J and T. Belytschko, 2007: A first course in finite elements,
!!                                  John Wiley and Sons
!!
!!
SUBROUTINE grad_fe_cell_adv_2d( p_cc, ptr_patch, ptr_int, p_grad, &
  &                             opt_rlstart, opt_rlend            )
!
!
!  patch on which computation is performed
!
TYPE(t_patch), TARGET, INTENT(in)     :: ptr_patch
!
!  data structure for interpolation
!
TYPE(t_int_state), TARGET, INTENT(in) :: ptr_int

!
!  cell centered variable
!
REAL(wp), INTENT(in) ::  &
  &  p_cc(:,:)

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_rlstart, opt_rlend   ! start and end values of refin_ctrl flag
!
! cell based Green-Gauss reconstructed geographical gradient vector
!
REAL(wp), INTENT(inout) ::  &
  &  p_grad(:,:,:)      ! dim:(2,nproma,nblks_c)

INTEGER :: jc, jb
INTEGER :: rl_start, rl_end
INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx, i_nchdom

INTEGER,  DIMENSION(:,:,:),   POINTER :: iidx, iblk

!-----------------------------------------------------------------------

! check optional arguments
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


iidx => ptr_patch%cells%neighbor_idx
iblk => ptr_patch%cells%neighbor_blk

i_nchdom = MAX(1,ptr_patch%n_childdom)


!
! 2. reconstruction of cell based geographical gradient
!
!$OMP PARALLEL PRIVATE(i_startblk,i_endblk)

  i_startblk = ptr_patch%cells%start_blk(rl_start,1)
  i_endblk   = ptr_patch%cells%end_blk(rl_end,i_nchdom)

  IF (ptr_patch%id > 1) THEN
  ! Fill nest boundaries with zero to avoid trouble with MPI synchronization
    CALL init(p_grad(:,:,1:i_startblk))
!$OMP BARRIER
  ENDIF

!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx), ICON_OMP_RUNTIME_SCHEDULE
  DO jb = i_startblk, i_endblk

    CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, rl_start, rl_end)


    DO jc = i_startidx, i_endidx

      ! We do not make use of the intrinsic function DOT_PRODUCT on purpose,
      ! since it is extremely slow on the SX9, when combined with indirect
      ! addressing.

      ! multiply cell-based input values with precomputed grid geometry factor

      ! zonal(u)-component of gradient
      p_grad(1,jc,jb) = &
        &    ptr_int%gradc_bmat(jc,1,1,jb)*p_cc(iidx(jc,jb,1),iblk(jc,jb,1))  &
        &  + ptr_int%gradc_bmat(jc,1,2,jb)*p_cc(iidx(jc,jb,2),iblk(jc,jb,2))  &
        &  + ptr_int%gradc_bmat(jc,1,3,jb)*p_cc(iidx(jc,jb,3),iblk(jc,jb,3))

      ! meridional(v)-component of gradient
      p_grad(2,jc,jb) =  &
        &    ptr_int%gradc_bmat(jc,2,1,jb)*p_cc(iidx(jc,jb,1),iblk(jc,jb,1))  &
        &  + ptr_int%gradc_bmat(jc,2,2,jb)*p_cc(iidx(jc,jb,2),iblk(jc,jb,2))  &
        &  + ptr_int%gradc_bmat(jc,2,3,jb)*p_cc(iidx(jc,jb,3),iblk(jc,jb,3))

    END DO ! end loop over cells

  END DO ! end loop over blocks
!$OMP END DO NOWAIT
!$OMP END PARALLEL


END SUBROUTINE grad_fe_cell_adv_2d



!-------------------------------------------------------------------------
!
!
!>
!! Computes the cell centered gradient in geographical coordinates.
!!
!! The gradient is computed by taking the derivative of the shape functions
!! for a three-node triangular element (Finite Element thinking).
!! Special dycore version, which handles two fields at a time.
!!
!! @par Revision History
!!  Initial revision by Daniel Reinert, DWD (2013-11-07)
!!
!! LITERATURE:
!! Fish. J and T. Belytschko, 2007: A first course in finite elements,
!!                                  John Wiley and Sons
!!
!!
SUBROUTINE grad_fe_cell_dycore( p_ccpr, ptr_patch, ptr_int, p_grad, &
  &                      opt_slev, opt_elev, opt_rlstart,  &
  &                      opt_rlend                         )
!
!
!  patch on which computation is performed
!
TYPE(t_patch), TARGET, INTENT(in)     :: ptr_patch
!
!  data structure for interpolation
!
TYPE(t_int_state), TARGET, INTENT(in) :: ptr_int

REAL(vp), INTENT(in) ::  & ! perturbation fields passed from dycore (nproma,2,nlev,nblks_c)
  &  p_ccpr(:,:,:,:)

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_slev    ! optional vertical start level

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_elev    ! optional vertical end level

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_rlstart, opt_rlend   ! start and end values of refin_ctrl flag
!
! cell based Green-Gauss reconstructed geographical gradient vector
!
REAL(vp), INTENT(inout) ::  &
  &  p_grad(:,:,:,:)      ! dim:(4,nproma,nlev,nblks_c)

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
  elev = UBOUND(p_ccpr,3)
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


iidx => ptr_patch%cells%neighbor_idx
iblk => ptr_patch%cells%neighbor_blk

i_nchdom = MAX(1,ptr_patch%n_childdom)


!
! 2. reconstruction of cell based geographical gradient
!

!$ACC DATA PCOPYIN( p_ccpr ) PCOPYOUT( p_grad )                                      &
!$ACC      PRESENT( ptr_int, iidx, iblk ) IF( i_am_accel_node .AND. acc_on )
!$ACC UPDATE DEVICE( p_ccpr ), IF( i_am_accel_node .AND. acc_on .AND. acc_validate )

!$OMP PARALLEL PRIVATE(i_startblk,i_endblk)

  i_startblk = ptr_patch%cells%start_blk(rl_start,1)
  i_endblk   = ptr_patch%cells%end_blk(rl_end,i_nchdom)


!$OMP DO PRIVATE(jb,jc,jk,i_startidx,i_endidx), ICON_OMP_RUNTIME_SCHEDULE
  DO jb = i_startblk, i_endblk

    CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, rl_start, rl_end)

    !$ACC PARALLEL IF( i_am_accel_node .AND. acc_on )
#ifdef __LOOP_EXCHANGE
    !$ACC LOOP GANG
    DO jc = i_startidx, i_endidx
!DIR$ IVDEP
      !$ACC LOOP VECTOR
      DO jk = slev, elev
#else
    !$ACC LOOP GANG
    DO jk = slev, elev
      !$ACC LOOP VECTOR
      DO jc = i_startidx, i_endidx
#endif

        ! We do not make use of the intrinsic function DOT_PRODUCT on purpose,
        ! since it is extremely slow on the SX9, when combined with indirect
        ! addressing.

        ! multiply cell-based input values with shape function derivatives

        ! zonal(u)-component of gradient, field 1
        p_grad(1,jc,jk,jb) = &
          &    ptr_int%gradc_bmat(jc,1,1,jb)*p_ccpr(1,iidx(jc,jb,1),jk,iblk(jc,jb,1))  &
          &  + ptr_int%gradc_bmat(jc,1,2,jb)*p_ccpr(1,iidx(jc,jb,2),jk,iblk(jc,jb,2))  &
          &  + ptr_int%gradc_bmat(jc,1,3,jb)*p_ccpr(1,iidx(jc,jb,3),jk,iblk(jc,jb,3))

        ! meridional(v)-component of gradient, field 1
        p_grad(2,jc,jk,jb) =  &
          &    ptr_int%gradc_bmat(jc,2,1,jb)*p_ccpr(1,iidx(jc,jb,1),jk,iblk(jc,jb,1))  &
          &  + ptr_int%gradc_bmat(jc,2,2,jb)*p_ccpr(1,iidx(jc,jb,2),jk,iblk(jc,jb,2))  &
          &  + ptr_int%gradc_bmat(jc,2,3,jb)*p_ccpr(1,iidx(jc,jb,3),jk,iblk(jc,jb,3))

        ! zonal(u)-component of gradient, field 2
        p_grad(3,jc,jk,jb) = &
          &    ptr_int%gradc_bmat(jc,1,1,jb)*p_ccpr(2,iidx(jc,jb,1),jk,iblk(jc,jb,1))  &
          &  + ptr_int%gradc_bmat(jc,1,2,jb)*p_ccpr(2,iidx(jc,jb,2),jk,iblk(jc,jb,2))  &
          &  + ptr_int%gradc_bmat(jc,1,3,jb)*p_ccpr(2,iidx(jc,jb,3),jk,iblk(jc,jb,3))

        ! meridional(v)-component of gradient, field 2
        p_grad(4,jc,jk,jb) =  &
          &    ptr_int%gradc_bmat(jc,2,1,jb)*p_ccpr(2,iidx(jc,jb,1),jk,iblk(jc,jb,1))  &
          &  + ptr_int%gradc_bmat(jc,2,2,jb)*p_ccpr(2,iidx(jc,jb,2),jk,iblk(jc,jb,2))  &
          &  + ptr_int%gradc_bmat(jc,2,3,jb)*p_ccpr(2,iidx(jc,jb,3),jk,iblk(jc,jb,3))

      END DO ! end loop over cells
    END DO ! end loop over vertical levels
    !$ACC END PARALLEL

  END DO ! end loop over blocks
!$OMP END DO NOWAIT
!$OMP END PARALLEL

!$ACC UPDATE HOST( p_grad ), IF( i_am_accel_node .AND. acc_on .AND. acc_validate )
!$ACC END DATA

END SUBROUTINE grad_fe_cell_dycore


!-------------------------------------------------------------------------
!
!
!>
!! Computes the cell centered gradient in geographical coordinates.
!!
!! The Green-Gauss approach is used. See for example:
!! http://www.cfd-online.com/Wiki/Gradient_computation
!!
!! @par Revision History
!!  Developed by Daniel Reinert (2009-07-21)
!!  Modification by Almut Gassmann (2010-01-12)
!!  - generalize for hexagons
!!  Modification by Guenther Zaengl (2010-03-09)
!!  - optimization by using precomputed coefficients
!!
SUBROUTINE grad_green_gauss_cell_adv( p_cc, ptr_patch, ptr_int, p_grad, &
  &                                   opt_slev, opt_elev, opt_p_face,   &
  &                                   opt_rlstart, opt_rlend            )
!
!
!  patch on which computation is performed
!
TYPE(t_patch), TARGET, INTENT(in)     :: ptr_patch
!
!  data structure for interpolation
!
TYPE(t_int_state), TARGET, INTENT(in) :: ptr_int

!
!  cell centered variable
!
REAL(wp), INTENT(in) ::  &
  &  p_cc(:,:,:)

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_slev    ! optional vertical start level

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_elev    ! optional vertical end level

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_rlstart, opt_rlend   ! start and end values of refin_ctrl flag
!
! cell based Green-Gauss reconstructed geographical gradient vector
!
REAL(vp), INTENT(inout) ::  &
  &  p_grad(:,:,:,:)      ! dim:(2,nproma,nlev,nblks_c)

! optional: calculated face values of cell centered quantity
REAL(wp), INTENT(inout), OPTIONAL ::  &
  &  opt_p_face(:,:,:)

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


iidx => ptr_patch%cells%neighbor_idx
iblk => ptr_patch%cells%neighbor_blk

i_nchdom = MAX(1,ptr_patch%n_childdom)


! save face values in optional output field
! (the cell-to-edge interpolation is no longer needed otherwise because
!  of using precomputed geometrical factors)
IF ( PRESENT(opt_p_face) ) THEN
  CALL cells2edges_scalar( p_cc, ptr_patch, ptr_int%c_lin_e, opt_p_face,  &
    &                      slev, elev)
ENDIF


!
! 2. reconstruction of cell based geographical gradient
!
!$ACC DATA PCOPYIN( p_cc ) PCOPYOUT( p_grad )                                  &
!$ACC      PRESENT( ptr_int, iidx, iblk ) IF( i_am_accel_node .AND. acc_on )
!$ACC UPDATE DEVICE( p_cc ), IF( i_am_accel_node .AND. acc_on .AND. acc_validate )

!$OMP PARALLEL PRIVATE(i_startblk,i_endblk)

  i_startblk = ptr_patch%cells%start_blk(rl_start,1)
  i_endblk   = ptr_patch%cells%end_blk(rl_end,i_nchdom)

  IF (ptr_patch%id > 1) THEN
  ! Fill nest boundaries with zero to avoid trouble with MPI synchronization
#ifdef _OPENACC
!$ACC KERNELS IF( i_am_accel_node .AND. acc_on )
    p_grad(:,:,:,1:i_startblk) = 0._wp
!$ACC END KERNELS
#else
    CALL init(p_grad(:,:,:,1:i_startblk))
!$OMP BARRIER
#endif
  ENDIF

!$OMP DO PRIVATE(jb,jc,jk,i_startidx,i_endidx), ICON_OMP_RUNTIME_SCHEDULE
  DO jb = i_startblk, i_endblk

    CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, rl_start, rl_end)

    !$ACC PARALLEL IF( i_am_accel_node .AND. acc_on )
#ifdef __LOOP_EXCHANGE
    !$ACC LOOP GANG
    DO jc = i_startidx, i_endidx
!DIR$ IVDEP
      !$ACC LOOP VECTOR
      DO jk = slev, elev
#else
    !$ACC LOOP GANG
    DO jk = slev, elev
      !$ACC LOOP VECTOR
      DO jc = i_startidx, i_endidx
#endif

        ! multiply cell-based input values with precomputed grid geometry factor

        ! zonal(u)-component of Green-Gauss gradient
        p_grad(1,jc,jk,jb) = ptr_int%geofac_grg(jc,1,jb,1)*p_cc(jc,jk,jb)    + &
          ptr_int%geofac_grg(jc,2,jb,1)*p_cc(iidx(jc,jb,1),jk,iblk(jc,jb,1)) + &
          ptr_int%geofac_grg(jc,3,jb,1)*p_cc(iidx(jc,jb,2),jk,iblk(jc,jb,2)) + &
          ptr_int%geofac_grg(jc,4,jb,1)*p_cc(iidx(jc,jb,3),jk,iblk(jc,jb,3))

        ! meridional(v)-component of Green-Gauss gradient
        p_grad(2,jc,jk,jb) = ptr_int%geofac_grg(jc,1,jb,2)*p_cc(jc,jk,jb)    + &
          ptr_int%geofac_grg(jc,2,jb,2)*p_cc(iidx(jc,jb,1),jk,iblk(jc,jb,1)) + &
          ptr_int%geofac_grg(jc,3,jb,2)*p_cc(iidx(jc,jb,2),jk,iblk(jc,jb,2)) + &
          ptr_int%geofac_grg(jc,4,jb,2)*p_cc(iidx(jc,jb,3),jk,iblk(jc,jb,3))

      END DO ! end loop over cells
    END DO ! end loop over vertical levels
    !$ACC END PARALLEL

  END DO ! end loop over blocks
!$OMP END DO NOWAIT
!$OMP END PARALLEL

!$ACC UPDATE HOST( p_grad ), IF( i_am_accel_node .AND. acc_on .AND. acc_validate )
!$ACC END DATA

  END SUBROUTINE grad_green_gauss_cell_adv

  SUBROUTINE grad_green_gauss_cell_dycore(p_ccpr, ptr_patch, ptr_int, p_grad,       &
    &                                     opt_slev, opt_elev, opt_rlstart, opt_rlend)
  !
  !
  !  patch on which computation is performed
  !
  TYPE(t_patch), TARGET, INTENT(in)     :: ptr_patch
  !
  !  data structure for interpolation
  !
  TYPE(t_int_state), TARGET, INTENT(in) :: ptr_int

  !  cell centered I/O variables
  !
  REAL(vp), INTENT(in) :: p_ccpr(:,:,:,:) ! perturbation fields passed from dycore (2,nproma,nlev,nblks_c)

  INTEGER, INTENT(in), OPTIONAL :: opt_slev    ! optional vertical start level

  INTEGER, INTENT(in), OPTIONAL :: opt_elev    ! optional vertical end level

  INTEGER, INTENT(in), OPTIONAL :: opt_rlstart, opt_rlend   ! start and end values of refin_ctrl flag

  !
  ! cell based Green-Gauss reconstructed geographical gradient vector
  !
  REAL(vp), INTENT(inout) :: p_grad(:,:,:,:)      ! dim:(4,nproma,nlev,nblks_c)

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
    elev = UBOUND(p_ccpr,3)
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


  iidx => ptr_patch%cells%neighbor_idx
  iblk => ptr_patch%cells%neighbor_blk

  i_nchdom = MAX(1,ptr_patch%n_childdom)

  !
  ! 2. reconstruction of cell based geographical gradient
  !

!$ACC DATA PCOPYIN( p_ccpr ) PCOPYOUT( p_grad )                                     &
!$ACC      PRESENT( ptr_int, iidx, iblk ) IF( i_am_accel_node .AND. acc_on )
!$ACC UPDATE DEVICE( p_ccpr ), IF( i_am_accel_node .AND. acc_on .AND. acc_validate )

!$OMP PARALLEL PRIVATE(i_startblk,i_endblk)

    i_startblk = ptr_patch%cells%start_blk(rl_start,1)
    i_endblk   = ptr_patch%cells%end_blk(rl_end,i_nchdom)

!$OMP DO PRIVATE(jb,jc,jk,i_startidx,i_endidx), ICON_OMP_RUNTIME_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, rl_start, rl_end)

      !$ACC PARALLEL IF( i_am_accel_node .AND. acc_on )
#ifdef __LOOP_EXCHANGE
      !$ACC LOOP GANG
      DO jc = i_startidx, i_endidx
!DIR$ IVDEP
        !$ACC LOOP VECTOR
        DO jk = slev, elev
#else
      !$ACC LOOP GANG
      DO jk = slev, elev
        !$ACC LOOP VECTOR
        DO jc = i_startidx, i_endidx
#endif

          ! zonal(u)-component of Green-Gauss gradient, field 1
          p_grad(1,jc,jk,jb) = ptr_int%geofac_grg(jc,1,jb,1)*p_ccpr(1,jc,jk,jb)+     &
            ptr_int%geofac_grg(jc,2,jb,1)*p_ccpr(1,iidx(jc,jb,1),jk,iblk(jc,jb,1)) + &
            ptr_int%geofac_grg(jc,3,jb,1)*p_ccpr(1,iidx(jc,jb,2),jk,iblk(jc,jb,2)) + &
            ptr_int%geofac_grg(jc,4,jb,1)*p_ccpr(1,iidx(jc,jb,3),jk,iblk(jc,jb,3))

          ! meridional(v)-component of Green-Gauss gradient, field 1
          p_grad(2,jc,jk,jb) = ptr_int%geofac_grg(jc,1,jb,2)*p_ccpr(1,jc,jk,jb)    + &
            ptr_int%geofac_grg(jc,2,jb,2)*p_ccpr(1,iidx(jc,jb,1),jk,iblk(jc,jb,1)) + &
            ptr_int%geofac_grg(jc,3,jb,2)*p_ccpr(1,iidx(jc,jb,2),jk,iblk(jc,jb,2)) + &
            ptr_int%geofac_grg(jc,4,jb,2)*p_ccpr(1,iidx(jc,jb,3),jk,iblk(jc,jb,3))

          ! zonal(u)-component of Green-Gauss gradient, field 2
          p_grad(3,jc,jk,jb) = ptr_int%geofac_grg(jc,1,jb,1)*p_ccpr(2,jc,jk,jb)    + &
            ptr_int%geofac_grg(jc,2,jb,1)*p_ccpr(2,iidx(jc,jb,1),jk,iblk(jc,jb,1)) + &
            ptr_int%geofac_grg(jc,3,jb,1)*p_ccpr(2,iidx(jc,jb,2),jk,iblk(jc,jb,2)) + &
            ptr_int%geofac_grg(jc,4,jb,1)*p_ccpr(2,iidx(jc,jb,3),jk,iblk(jc,jb,3))

          ! meridional(v)-component of Green-Gauss gradient, field 2
          p_grad(4,jc,jk,jb) = ptr_int%geofac_grg(jc,1,jb,2)*p_ccpr(2,jc,jk,jb)    + &
            ptr_int%geofac_grg(jc,2,jb,2)*p_ccpr(2,iidx(jc,jb,1),jk,iblk(jc,jb,1)) + &
            ptr_int%geofac_grg(jc,3,jb,2)*p_ccpr(2,iidx(jc,jb,2),jk,iblk(jc,jb,2)) + &
            ptr_int%geofac_grg(jc,4,jb,2)*p_ccpr(2,iidx(jc,jb,3),jk,iblk(jc,jb,3))

        END DO ! end loop over cells
      END DO ! end loop over vertical levels
      !$ACC END PARALLEL

    END DO ! end loop over blocks

!$OMP END DO NOWAIT
!$OMP END PARALLEL

!$ACC UPDATE HOST( p_grad ), IF( i_am_accel_node .AND. acc_on .AND. acc_validate )
!$ACC END DATA

  END SUBROUTINE grad_green_gauss_cell_dycore


END MODULE mo_math_gradients
