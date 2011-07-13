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
!!  Modified by Hui Wan, MPI-M, (2008-04-04)
!!  - control variable loce renamed locean
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
USE mo_kind,               ONLY: wp
USE mo_impl_constants,     ONLY: min_rlcell, min_rledge
USE mo_interpolation,      ONLY: t_int_state, cells2edges_scalar
USE mo_model_domain,       ONLY: t_patch
USE mo_nonhydrostatic_nml, ONLY: upstr_beta
USE mo_parallel_configuration,  ONLY: nproma
USE mo_run_config,         ONLY: ltimer
USE mo_exception,          ONLY: finish
USE mo_timer,              ONLY: timer_start, timer_stop, timer_grad
USE mo_loopindices,        ONLY: get_indices_c, get_indices_e

IMPLICIT NONE

PRIVATE

CHARACTER(len=*), PARAMETER :: version = '$Id$'


PUBLIC :: grad_fd_norm, grad_fd_tang
PUBLIC :: grad_green_gauss_cell
PUBLIC :: grad_dir_edge


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
INTEGER :: nlen, nblks_e, npromz_e

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
  elev = ptr_patch%nlev
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

IF (ltimer) CALL timer_start(timer_grad)

!$OMP PARALLEL
! The special treatment of 2D fields is essential for efficiency on the NEC

SELECT CASE (ptr_patch%cell_type)

CASE (3) ! (cell_type == 3)


#ifdef __SX__
IF (slev > 1) THEN
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,je,jk)
  DO jb = i_startblk, i_endblk

  CALL get_indices_e(ptr_patch, jb, i_startblk, i_endblk, &
                     i_startidx, i_endidx, rl_start, rl_end)

#ifdef _URD2
!CDIR UNROLL=_URD2
#endif
#ifdef __LOOP_EXCHANGE
    DO je = i_startidx, i_endidx
      DO jk = slev, elev
#else
    DO jk = slev, elev
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

  END DO
!$OMP END DO
ELSE
#endif

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,je,jk)
  DO jb = i_startblk, i_endblk

  CALL get_indices_e(ptr_patch, jb, i_startblk, i_endblk, &
                     i_startidx, i_endidx, rl_start, rl_end)

#ifdef _URD
!CDIR UNROLL=_URD
#endif

#ifdef __LOOP_EXCHANGE
    DO je = i_startidx, i_endidx
      DO jk = slev, elev
#else
    DO jk = slev, elev
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

  END DO
!$OMP END DO
#ifdef __SX__
ENDIF
#endif
CASE (6) ! (cell_type == 6)

  ! no grid refinement in hexagonal model
  nblks_e   = ptr_patch%nblks_int_e
  npromz_e  = ptr_patch%npromz_int_e

!$OMP DO PRIVATE(jb,nlen,je,jk)
  DO jb = 1, nblks_e

    IF (jb /= nblks_e) THEN
      nlen = nproma
    ELSE
      nlen = npromz_e
    ENDIF

#ifdef __LOOP_EXCHANGE
    DO je = 1, nlen
      DO jk = slev, elev
#else
    DO jk = slev, elev
      DO je = 1, nlen
#endif
      !
      ! compute the normal derivative
      ! by the finite difference approximation
      ! (see Bonaventura and Ringler MWR 2005)
      !
        grad_norm_psi_e(je,jk,jb) =  &
          &  ( psi_c(iidx(je,jb,2),jk,iblk(je,jb,2)) - &
          &    psi_c(iidx(je,jb,1),jk,iblk(je,jb,1)) )  &
          &  * ptr_patch%edges%inv_dual_edge_length(je,jb) &
          &  * ptr_patch%edges%system_orientation(je,jb)

    ENDDO

  END DO

END DO

!$OMP END DO
END SELECT
!$OMP END PARALLEL

IF(ltimer) CALL timer_stop(timer_grad)


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
  elev =ptr_patch%nlev
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

!
!  loop through all patch edges (and blocks)
!
DO jb = i_startblk, i_endblk

  CALL get_indices_e(ptr_patch, jb, i_startblk, i_endblk, &
                     i_startidx, i_endidx, rl_start, rl_end)

  DO je = i_startidx, i_endidx
    !
    !  get the line and block indices of the vertices of edge je
    !
    ilv1(je) = ptr_patch%edges%vertex_idx(je,jb,1)
    ibv1(je) = ptr_patch%edges%vertex_blk(je,jb,1)
    ilv2(je) = ptr_patch%edges%vertex_idx(je,jb,2)
    ibv2(je) = ptr_patch%edges%vertex_blk(je,jb,2)
  END DO

  DO jk = slev, elev

    DO je = i_startidx, i_endidx
      !
      ! compute the tangential derivative
      ! by the finite difference approximation
      iorient = ptr_patch%edges%system_orientation(je,jb)
      grad_tang_psi_e(je,jk,jb) = iorient  &
        &  * ( psi_v(ilv2(je),jk,ibv2(je)) - psi_v(ilv1(je),jk,ibv1(je)) )  &
        &    / ptr_patch%edges%primal_edge_length(je,jb)
    END DO

  END DO

END DO

END SUBROUTINE grad_fd_tang

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
SUBROUTINE grad_green_gauss_cell( p_cc, ptr_patch, ptr_int, p_grad,  &
  &                               opt_slev, opt_elev,   &
  &                               opt_p_face, opt_rlstart, &
  &                               opt_rlend )
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
REAL(wp), INTENT(inout) ::  &
  &  p_grad(:,:,:,:)

! optional: calculated face values of cell centered quantity
REAL(wp), INTENT(inout), OPTIONAL ::  &
  &  opt_p_face(:,:,:)

INTEGER :: slev, elev     ! vertical start and end level
INTEGER :: jc, jk, jb
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
  rl_start = 2
END IF
IF ( PRESENT(opt_rlend) ) THEN
  rl_end = opt_rlend
ELSE
  rl_end = min_rlcell
END IF

iidx => ptr_patch%cells%neighbor_idx
iblk => ptr_patch%cells%neighbor_blk

! values for the blocking
i_nchdom   = MAX(1,ptr_patch%n_childdom)
i_startblk = ptr_patch%cells%start_blk(rl_start,1)
i_endblk   = ptr_patch%cells%end_blk(rl_end,i_nchdom)


! save face values in optional output field
! (the cell-to-edge interpolation is no longer needed otherwise because
!  of using precomputed geometrical factors)
IF ( PRESENT(opt_p_face) ) THEN
  CALL cells2edges_scalar( p_cc, ptr_patch, ptr_int%c_lin_e, opt_p_face,  &
    &                      slev, elev)
ENDIF

SELECT CASE (ptr_patch%cell_type)

  CASE (3) ! (cell_type == 3)
!
! 2. reconstruction of cell based geographical gradient
!
!$OMP PARALLEL

  ! Fill nest boundaries with zero to avoid trouble with MPI synchronization
  IF (ptr_patch%id > 1) THEN
!$OMP WORKSHARE
    p_grad(:,:,1:i_startblk,1) = 0._wp
    p_grad(:,:,1:i_startblk,2) = 0._wp
!$OMP END WORKSHARE
  ENDIF

!$OMP DO PRIVATE(jb,jc,jk,i_startidx,i_endidx), SCHEDULE(runtime)
  DO jb = i_startblk, i_endblk

    CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, rl_start, rl_end)

#ifdef __LOOP_EXCHANGE
    DO jc = i_startidx, i_endidx
      DO jk = slev, elev
#else
!CDIR UNROLL=3
    DO jk = slev, elev
      DO jc = i_startidx, i_endidx
#endif

        ! multiply cell-based input values with precomputed grid geometry factor

        ! zonal(u)-component of Green-Gauss gradient
        p_grad(jc,jk,jb,1) = ptr_int%geofac_grg(jc,1,jb,1)*p_cc(jc,jk,jb)    + &
          ptr_int%geofac_grg(jc,2,jb,1)*p_cc(iidx(jc,jb,1),jk,iblk(jc,jb,1)) + &
          ptr_int%geofac_grg(jc,3,jb,1)*p_cc(iidx(jc,jb,2),jk,iblk(jc,jb,2)) + &
          ptr_int%geofac_grg(jc,4,jb,1)*p_cc(iidx(jc,jb,3),jk,iblk(jc,jb,3))

        ! meridional(v)-component of Green-Gauss gradient
        p_grad(jc,jk,jb,2) = ptr_int%geofac_grg(jc,1,jb,2)*p_cc(jc,jk,jb)    + &
          ptr_int%geofac_grg(jc,2,jb,2)*p_cc(iidx(jc,jb,1),jk,iblk(jc,jb,1)) + &
          ptr_int%geofac_grg(jc,3,jb,2)*p_cc(iidx(jc,jb,2),jk,iblk(jc,jb,2)) + &
          ptr_int%geofac_grg(jc,4,jb,2)*p_cc(iidx(jc,jb,3),jk,iblk(jc,jb,3))

      END DO ! end loop over cells

    END DO ! end loop over vertical levels

  END DO ! end loop over blocks
!$OMP END DO
!$OMP END PARALLEL

CASE(6) ! (cell_type == 6)
!
! 2. reconstruction of cell based geographical gradient
!
  nblks_c = ptr_patch%nblks_int_c
  npromz_c = ptr_patch%npromz_int_c
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,jk,nlen)
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
    DO jk = slev, elev
      DO jc = 1, nlen
#endif

        ! multiply cell-based input values with precomputed grid geometry factor

        ! zonal(u)-component of Green-Gauss gradient
        p_grad(jc,jk,jb,1) = ptr_int%geofac_grg(jc,1,jb,1)*p_cc(jc,jk,jb)    + &
          ptr_int%geofac_grg(jc,2,jb,1)*p_cc(iidx(jc,jb,1),jk,iblk(jc,jb,1)) + &
          ptr_int%geofac_grg(jc,3,jb,1)*p_cc(iidx(jc,jb,2),jk,iblk(jc,jb,2)) + &
          ptr_int%geofac_grg(jc,4,jb,1)*p_cc(iidx(jc,jb,3),jk,iblk(jc,jb,3)) + &
          ptr_int%geofac_grg(jc,5,jb,1)*p_cc(iidx(jc,jb,4),jk,iblk(jc,jb,4)) + &
          ptr_int%geofac_grg(jc,6,jb,1)*p_cc(iidx(jc,jb,5),jk,iblk(jc,jb,5)) + &
          ptr_int%geofac_grg(jc,7,jb,1)*p_cc(iidx(jc,jb,6),jk,iblk(jc,jb,6))

        ! meridional(v)-component of Green-Gauss gradient
        p_grad(jc,jk,jb,2) = ptr_int%geofac_grg(jc,1,jb,2)*p_cc(jc,jk,jb)    + &
          ptr_int%geofac_grg(jc,2,jb,2)*p_cc(iidx(jc,jb,1),jk,iblk(jc,jb,1)) + &
          ptr_int%geofac_grg(jc,3,jb,2)*p_cc(iidx(jc,jb,2),jk,iblk(jc,jb,2)) + &
          ptr_int%geofac_grg(jc,4,jb,2)*p_cc(iidx(jc,jb,3),jk,iblk(jc,jb,3)) + &
          ptr_int%geofac_grg(jc,5,jb,2)*p_cc(iidx(jc,jb,4),jk,iblk(jc,jb,4)) + &
          ptr_int%geofac_grg(jc,6,jb,2)*p_cc(iidx(jc,jb,5),jk,iblk(jc,jb,5)) + &
          ptr_int%geofac_grg(jc,7,jb,2)*p_cc(iidx(jc,jb,6),jk,iblk(jc,jb,6))

      END DO ! end loop over cells

    END DO ! end loop over vertical levels

  END DO ! end loop over blocks
!$OMP END DO
!$OMP END PARALLEL
END SELECT
END SUBROUTINE grad_green_gauss_cell


  !----------------------------------------------------------------------
  !>
  !! directional gradient of a vector(=gradient(scalar)) in the edge direction
  !!
  !! For poor men's third order advection of theta_v or other scalars one
  !! needs the laplacian in the edge direction. Here, the gradient in the edge
  !! direction of a given gradient field is computed. According to the given
  !! velocity field, the gradient is computed using the upstream located cell.
  !! The coefficients are precomputed in mo_intp_coeffs:init_geo_factors.
  !!
  !! CURRENTLY ONLY FOR HEXAGONAL MODEL
  !!
  !! @par Revision History
  !! Initial release by Almut Gassmann, MPI_M (2010-02-08)
  !!
  SUBROUTINE grad_dir_edge (p_vn, p_grads, pt_patch, pt_int, p_gradv)

    TYPE(t_patch), INTENT(in) :: pt_patch           !< patch
    TYPE(t_int_state), TARGET, INTENT(in) :: pt_int !< interpolation state
    REAL(wp), INTENT(in)    :: p_vn(:,:,:)  &
    & ; !< normal velocity field, needed for determination of upstream direction
    REAL(wp), INTENT(in)    :: p_grads(:,:,:) !< gradient of scalar(=vector)
    REAL(wp), INTENT(out)   :: p_gradv(:,:,:) !< gradient of vector(=grads)

    INTEGER :: nblks_e, npromz_e, nlen, jb, jk, je
    INTEGER :: nlev              !< number of full levels
    INTEGER, POINTER :: ii1(:,:,:), ib1(:,:,:), ii2(:,:,:), ib2(:,:,:)
#ifdef __LOOP_EXCHANGE
    REAL(wp) :: z_sign(pt_patch%nlev)
#else
    REAL(wp) :: z_sign(nproma)
#endif
    !-----------------------------------------------------------------

    ii1 => pt_int%dir_gradh_i1
    ib1 => pt_int%dir_gradh_b1
    ii2 => pt_int%dir_gradh_i2
    ib2 => pt_int%dir_gradh_b2

    nblks_e  = pt_patch%nblks_int_e
    npromz_e = pt_patch%npromz_int_e

    ! number of vertical levels
    nlev = pt_patch%nlev

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,nlen,jk,je,z_sign)
    DO jb = 1, nblks_e
      IF (jb /= nblks_e) THEN
        nlen = nproma
      ELSE
        nlen = npromz_e
      ENDIF
#ifdef __LOOP_EXCHANGE
      DO je = 1, nlen
        DO jk = 1, nlev
          z_sign(jk) = upstr_beta*SIGN(0.5_wp,&
          &               p_vn(je,jk,jb)*pt_patch%edges%system_orientation(je,jb))
          p_gradv(je,jk,jb) = (0.5_wp+z_sign(jk))  &
          &     *(pt_int%dir_gradhux_c1(1,je,jb)*p_grads(ii1(1,je,jb),jk,ib1(1,je,jb)) &
          &      +pt_int%dir_gradhux_c1(2,je,jb)*p_grads(ii1(2,je,jb),jk,ib1(2,je,jb)) &
          &      +pt_int%dir_gradhux_c1(3,je,jb)*p_grads(ii1(3,je,jb),jk,ib1(3,je,jb)) &
          &      +pt_int%dir_gradhux_c1(4,je,jb)*p_grads(ii1(4,je,jb),jk,ib1(4,je,jb)) &
          &      +pt_int%dir_gradhux_c1(5,je,jb)*p_grads(ii1(5,je,jb),jk,ib1(5,je,jb)) &
          &      +pt_int%dir_gradhux_c1(6,je,jb)*p_grads(ii1(6,je,jb),jk,ib1(6,je,jb)))&
          &                 + (0.5_wp-z_sign(jk))  &
          &     *(pt_int%dir_gradhux_c2(1,je,jb)*p_grads(ii2(1,je,jb),jk,ib2(1,je,jb)) &
          &      +pt_int%dir_gradhux_c2(2,je,jb)*p_grads(ii2(2,je,jb),jk,ib2(2,je,jb)) &
          &      +pt_int%dir_gradhux_c2(3,je,jb)*p_grads(ii2(3,je,jb),jk,ib2(3,je,jb)) &
          &      +pt_int%dir_gradhux_c2(4,je,jb)*p_grads(ii2(4,je,jb),jk,ib2(4,je,jb)) &
          &      +pt_int%dir_gradhux_c2(5,je,jb)*p_grads(ii2(5,je,jb),jk,ib2(5,je,jb)) &
          &      +pt_int%dir_gradhux_c2(6,je,jb)*p_grads(ii2(6,je,jb),jk,ib2(6,je,jb)))

        ENDDO
      ENDDO
#else
      DO jk = 1, nlev
        DO je = 1, nlen
          z_sign(je) = upstr_beta*SIGN(0.5_wp,&
          &               p_vn(je,jk,jb)*pt_patch%edges%system_orientation(je,jb))
          p_gradv(je,jk,jb) = (0.5_wp+z_sign(je))  &
          &     *(pt_int%dir_gradhux_c1(1,je,jb)*p_grads(ii1(1,je,jb),jk,ib1(1,je,jb)) &
          &      +pt_int%dir_gradhux_c1(2,je,jb)*p_grads(ii1(2,je,jb),jk,ib1(2,je,jb)) &
          &      +pt_int%dir_gradhux_c1(3,je,jb)*p_grads(ii1(3,je,jb),jk,ib1(3,je,jb)) &
          &      +pt_int%dir_gradhux_c1(4,je,jb)*p_grads(ii1(4,je,jb),jk,ib1(4,je,jb)) &
          &      +pt_int%dir_gradhux_c1(5,je,jb)*p_grads(ii1(5,je,jb),jk,ib1(5,je,jb)) &
          &      +pt_int%dir_gradhux_c1(6,je,jb)*p_grads(ii1(6,je,jb),jk,ib1(6,je,jb)))&
          &                 + (0.5_wp-z_sign(je))  &
          &     *(pt_int%dir_gradhux_c2(1,je,jb)*p_grads(ii2(1,je,jb),jk,ib2(1,je,jb)) &
          &      +pt_int%dir_gradhux_c2(2,je,jb)*p_grads(ii2(2,je,jb),jk,ib2(2,je,jb)) &
          &      +pt_int%dir_gradhux_c2(3,je,jb)*p_grads(ii2(3,je,jb),jk,ib2(3,je,jb)) &
          &      +pt_int%dir_gradhux_c2(4,je,jb)*p_grads(ii2(4,je,jb),jk,ib2(4,je,jb)) &
          &      +pt_int%dir_gradhux_c2(5,je,jb)*p_grads(ii2(5,je,jb),jk,ib2(5,je,jb)) &
          &      +pt_int%dir_gradhux_c2(6,je,jb)*p_grads(ii2(6,je,jb),jk,ib2(6,je,jb)))

        ENDDO
      ENDDO
#endif
    ENDDO
!$OMP END DO
!$OMP END PARALLEL

  END SUBROUTINE grad_dir_edge

END MODULE mo_math_gradients
