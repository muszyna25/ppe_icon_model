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
!
!    ProTeX FORTRAN source: Style 2
!    modified for ICON project, DWD/MPI-M 2006
!
!-------------------------------------------------------------------------
!
!
!
USE mo_kind,               ONLY: wp
USE mo_run_nml,            ONLY: nproma, ltimer
USE mo_math_constants
USE mo_physical_constants
USE mo_impl_constants,     ONLY: land_boundary, boundary,  sea_boundary, &!land, sea,
  &                              min_rlcell, min_rledge, min_rlvert,max_char_length
USE mo_model_domain,       ONLY: t_patch
USE mo_ocean_nml,          ONLY: lviscous, n_zlev, iswm_oce
USE mo_run_nml,            ONLY: i_cell_type
USE mo_dynamics_nml,       ONLY: nold
USE mo_exception,          ONLY: finish
#ifndef __SX__
USE mo_timer,              ONLY: timer_start, timer_stop, timer_div, timer_grad
#endif
USE mo_loopindices,        ONLY: get_indices_c, get_indices_e, get_indices_v
USE mo_oce_state,          ONLY: t_hydro_ocean_state, t_hydro_ocean_diag
USE mo_math_utilities,     ONLY: gc2cc,gvec2cvec, t_cartesian_coordinates
IMPLICIT NONE

PRIVATE

CHARACTER(len=*), PARAMETER :: version = '$Id$'


PUBLIC :: grad_fd_norm_oce
PUBLIC :: grad_fd_norm_oce_2D
PUBLIC :: div_oce
PUBLIC :: rot_vertex_ocean
PUBLIC :: nabla2_vec_ocean
PUBLIC :: height_related_quantities

CONTAINS

!-------------------------------------------------------------------------
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
!! - abandon grid for the sake of patch
!!Boundary handling for triangles by P. Korn (2009)
!!
SUBROUTINE grad_fd_norm_oce( psi_c, ptr_patch, grad_norm_psi_e, &
  &                      opt_slev, opt_elev, opt_rlstart, opt_rlend )

!  patch on which computation is performed
!
TYPE(t_patch), TARGET, INTENT(in) :: ptr_patch
!
!  cell based variable of which normal derivative is computed
!
REAL(wp), INTENT(in) ::  &
  &  psi_c(:,:,:)       ! dim: (nproma,n_zlev,nblks_c)

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
  &  grad_norm_psi_e(:,:,:)  ! dim: (nproma,n_zlev,nblks_e)

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
  elev = n_zlev
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
#ifndef __SX__
IF (ltimer) CALL timer_start(timer_grad)
#endif
!$OMP PARALLEL
! The special treatment of 2D fields is essential for efficiency on the NEC

SELECT CASE (i_cell_type)

CASE (3) ! (i_cell_type == 3)


#ifdef __SX__
IF (slev > 1) THEN
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,je,jk)
  DO jb = i_startblk, i_endblk
    CALL get_indices_e(ptr_patch, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, rl_start, rl_end)
#ifdef _URD2
!CDIR UNROLL=_URD2
#endif
    DO jk = slev, elev
      DO je = i_startidx, i_endidx
      !
      ! compute the normal derivative
      ! by the finite difference approximation
      ! (see Bonaventura and Ringler MWR 2005)
      !
!  #slo# 2011-02-28 - correction (no change, see below)
         IF ( ptr_patch%patch_oce%lsm_oce_e(je,jk,jb) <= sea_boundary ) THEN
!          IF (     ptr_patch%patch_oce%lsm_oce_c(iidx(je,jb,2),jk,iblk(je,jb,2)) == sea&
!            &.AND. ptr_patch%patch_oce%lsm_oce_c(iidx(je,jb,1),jk,iblk(je,jb,1))==sea) THEN
          grad_norm_psi_e(je,jk,jb) =  &
            &  ( psi_c(iidx(je,jb,2),jk,iblk(je,jb,2)) - &
            &    psi_c(iidx(je,jb,1),jk,iblk(je,jb,1)) )  &
            &  * ptr_patch%edges%inv_dual_edge_length(je,jb)
        ELSE
          grad_norm_psi_e(je,jk,jb) =  0.0_wp
        ENDIF
!           IF ( (ptr_patch%patch_oce%lsm_oce_c(iidx(je,jb,2),jk,iblk(je,jb,2))&
!           & == ptr_patch%patch_oce%lsm_oce_c(iidx(je,jb,1),jk,iblk(je,jb,1)))&
!           &.AND.ptr_patch%patch_oce%lsm_oce_c(iidx(je,jb,1),jk,iblk(je,jb,1))<0) THEN
! 
!             IF(   ptr_patch%patch_oce%lsm_oce_c(iidx(je,jb,1),jk,iblk(je,jb,1))&
!                & /= ptr_patch%patch_oce%lsm_oce_e(je,jk,jb))THEN
!               write(*,*)'I WARNING: INCONSISTENT LSM', jk, je,jb,iidx(je,jb,2),iblk(je,jb,2),&
!                                                            & iidx(je,jb,1),iblk(je,jb,1)
!               write(*,*)'lsm values:cell1, cell2, edge:',&
!               &ptr_patch%patch_oce%lsm_oce_c(iidx(je,jb,1),jk,iblk(je,jb,1)),&
!               &ptr_patch%patch_oce%lsm_oce_c(iidx(je,jb,2),jk,iblk(je,jb,2)),&
!               &ptr_patch%patch_oce%lsm_oce_e(je,jk,jb)
!             ENDIF
!           ENDIF
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
    DO jk = slev, elev
      DO je = i_startidx, i_endidx
      !
      ! compute the normal derivative
      ! by the finite difference approximation
      ! (see Bonaventura and Ringler MWR 2005)
      !
!  #slo# 2011-02-28 - check of sea is identical, since sea_boundary doesn't exist on edges
          IF ( ptr_patch%patch_oce%lsm_oce_e(je,jk,jb) <= sea_boundary ) THEN
   !  IF (     ptr_patch%patch_oce%lsm_oce_c(iidx(je,jb,2),jk,iblk(je,jb,2)) <= sea_boundary&
   !    &.AND. ptr_patch%patch_oce%lsm_oce_c(iidx(je,jb,1),jk,iblk(je,jb,1))<=sea_boundary) THEN
          grad_norm_psi_e(je,jk,jb) =                     &
            &  ( psi_c(iidx(je,jb,2),jk,iblk(je,jb,2)) -  &
            &    psi_c(iidx(je,jb,1),jk,iblk(je,jb,1)) )  &
            &  * ptr_patch%edges%inv_dual_edge_length(je,jb)

! IF (jk==3 .and. grad_norm_psi_e(je,jk,jb)>0.0_wp) THEN
!   write(*,*) 'gradient:', je,jb,grad_norm_psi_e(je,jk,jb), &
!     &        ptr_patch%patch_oce%lsm_oce_c(iidx(je,jb,2),jk,iblk(je,jb,2)), &
!     &        ptr_patch%patch_oce%lsm_oce_c(iidx(je,jb,1),jk,iblk(je,jb,1)), &
!     &         psi_c(iidx(je,jb,2),jk,iblk(je,jb,2)), psi_c(iidx(je,jb,1),jk,iblk(je,jb,1))
! ENDIF
! IF(jk==3)THEN!.and.iblk(je,jb,2)==900)THEN
! write(900,*)'gradient:', je,jk,jb,grad_norm_psi_e(je,jk,jb),&
! &psi_c(iidx(je,jb,2),jk,iblk(je,jb,2)), psi_c(iidx(je,jb,1),jk,iblk(je,jb,1))
! ENDIF
        ELSE
          grad_norm_psi_e(je,jk,jb) =  0.0_wp
        ENDIF
!           IF ( (ptr_patch%patch_oce%lsm_oce_c(iidx(je,jb,2),jk,iblk(je,jb,2))&
!           & == ptr_patch%patch_oce%lsm_oce_c(iidx(je,jb,1),jk,iblk(je,jb,1)))&
!           &.AND.ptr_patch%patch_oce%lsm_oce_c(iidx(je,jb,1),jk,iblk(je,jb,1))<0) THEN
! 
!             IF(   ptr_patch%patch_oce%lsm_oce_c(iidx(je,jb,1),jk,iblk(je,jb,1))&
!                & /= ptr_patch%patch_oce%lsm_oce_e(je,jk,jb))THEN
!               write(*,*)'II WARNING: INCONSISTENT LSM', jk, je,jb,iidx(je,jb,2),iblk(je,jb,2),&
!                                                            & iidx(je,jb,1),iblk(je,jb,1)
!               write(*,*)'lsm values:cell1, cell2, edge:',&
!               &ptr_patch%patch_oce%lsm_oce_c(iidx(je,jb,1),jk,iblk(je,jb,1)),&
!               &ptr_patch%patch_oce%lsm_oce_c(iidx(je,jb,2),jk,iblk(je,jb,2)),&
!               &ptr_patch%patch_oce%lsm_oce_e(je,jk,jb)
!             ENDIF
!           ENDIF
      ENDDO
    END DO
  END DO
!write(*,*)'inside grad: lev 3:',maxval(grad_norm_psi_e(:,3,:)),minval(grad_norm_psi_e(:,3,:)) 
!$OMP END DO
#ifdef __SX__
ENDIF
#endif


CASE (6) ! (i_cell_type == 6)

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

    DO jk = slev, elev
      DO je = 1, nlen
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
#ifndef __SX__
IF (ltimer) CALL timer_stop(timer_grad)
#endif
END SUBROUTINE grad_fd_norm_oce
!-------------------------------------------------------------------------
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
!! - abandon grid for the sake of patch
!! Boundary handling for triangles by P. Korn (2009)
!!
SUBROUTINE grad_fd_norm_oce_2D( psi_c, ptr_patch, grad_norm_psi_e)

!

!
!  patch on which computation is performed
!
TYPE(t_patch), TARGET, INTENT(in) :: ptr_patch
!
!  cell based variable of which normal derivative is computed
!
REAL(wp), INTENT(in)    :: psi_c(:,:)             ! dim: (nproma,nblks_c)

!
!  edge based variable in which normal derivative is stored
!
!REAL(wp), INTENT(out) ::  &
REAL(wp), INTENT(inout) ::  grad_norm_psi_e(:,:)  ! dim: (nproma,nblks_e)

INTEGER :: slev, elev     ! vertical start and end level
INTEGER :: je, jb
INTEGER :: rl_start, rl_end
INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx
INTEGER :: nlen, nblks_e, npromz_e

INTEGER,  DIMENSION(:,:,:),   POINTER :: iidx, iblk

!
!-----------------------------------------------------------------------


slev = 1
elev = 1
rl_start = 2
rl_end = min_rledge

iidx => ptr_patch%edges%cell_idx
iblk => ptr_patch%edges%cell_blk


i_startblk = ptr_patch%edges%start_blk(rl_start,1)
i_endblk   = ptr_patch%edges%end_blk(rl_end,1)
!
!  loop through all patch edges (and blocks)
!
#ifndef __SX__
IF (ltimer) CALL timer_start(timer_grad)
#endif
!$OMP PARALLEL
! The special treatment of 2D fields is essential for efficiency on the NEC

SELECT CASE (i_cell_type)

CASE (3) ! (i_cell_type == 3)


!$OMP DO PRIVATE(jb,i_startidx,i_endidx,je)
  DO jb = i_startblk, i_endblk

  CALL get_indices_e(ptr_patch, jb, i_startblk, i_endblk, &
                     i_startidx, i_endidx, rl_start, rl_end)

#ifdef _URD
!CDIR UNROLL=_URD
#endif
      DO je = i_startidx, i_endidx
      !
      ! compute the normal derivative
      ! by the finite difference approximation
      ! (see Bonaventura and Ringler MWR 2005)
      !
         IF ( ptr_patch%patch_oce%lsm_oce_e(je,1,jb) <= sea_boundary ) THEN
          grad_norm_psi_e(je,jb) =  &
            &  ( psi_c(iidx(je,jb,2),iblk(je,jb,2)) - &
            &    psi_c(iidx(je,jb,1),iblk(je,jb,1)) )  &
            &  * ptr_patch%edges%inv_dual_edge_length(je,jb)
         ELSE
           grad_norm_psi_e(je,jb) =  0.0_wp
         ENDIF
      ENDDO

  END DO
!$OMP END DO

CASE (6) ! (i_cell_type == 6)

  ! no grid refinement in hexagonal model
  nblks_e   = ptr_patch%nblks_int_e
  npromz_e  = ptr_patch%npromz_int_e

!$OMP DO PRIVATE(jb,nlen,je)
  DO jb = 1, nblks_e

    IF (jb /= nblks_e) THEN
      nlen = nproma
    ELSE
      nlen = npromz_e
    ENDIF

      DO je = 1, nlen
      !
      ! compute the normal derivative
      ! by the finite difference approximation
      ! (see Bonaventura and Ringler MWR 2005)
      !
        grad_norm_psi_e(je,jb) =  &
          &  ( psi_c(iidx(je,jb,2),iblk(je,jb,2)) - &
          &    psi_c(iidx(je,jb,1),iblk(je,jb,1)) )  &
          &  * ptr_patch%edges%inv_dual_edge_length(je,jb) &
          &  * ptr_patch%edges%system_orientation(je,jb)

    ENDDO
  END DO

!$OMP END DO
END SELECT
!$OMP END PARALLEL
#ifndef __SX__
IF (ltimer) CALL timer_stop(timer_grad)
#endif
END SUBROUTINE grad_fd_norm_oce_2D
!-------------------------------------------------------------------------
!
!
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
!! - loop only over the inner cells of a patch, not any more over halo cells
!! Modifications by P. Korn, MPI-M(2007-2)
!! - Switch fom array arguments to pointers
!! Modification by Almut Gassmann, MPI-M (2007-04-20)
!! - abandon grid for the sake of patch
!! Modification by Guenther Zaengl, DWD (2009-03-17)
!! - vector optimization
!! Modification by Peter Korn, MPI-M    (2009)
!! - Boundary treatment for the ocean
!! Modification by Stephan Lorenz, MPI-M (2010-08-05)
!! - New boundary definition with inner and boundary points on land/sea
!!
SUBROUTINE div_oce( vec_e, ptr_patch, div_vec_c, &
  &             opt_slev, opt_elev, opt_rlstart, opt_rlend )
!

!
!  patch on which computation is performed
!
TYPE(t_patch), TARGET, INTENT(in) :: ptr_patch

!
! edge based variable of which divergence
! is computed
!
REAL(wp), INTENT(INOUT) ::  &
  &  vec_e(:,:,:) ! dim: (nproma,n_zlev,nblks_e)

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_slev    ! optional vertical start level

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_elev    ! optional vertical end level

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_rlstart, opt_rlend   ! start and end values of refin_ctrl flag

!
! cell based variable in which divergence is stored
!
REAL(wp), INTENT(INOUT) :: div_vec_c(:,:,:) ! dim: (nproma,n_zlev,nblks_c)

INTEGER :: slev, elev     ! vertical start and end level
INTEGER :: jc, jk, jb, je
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
  elev = n_zlev
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
i_endblk   = ptr_patch%cells%end_blk(rl_end,1)

! loop through all patch cells (and blocks)
!
#ifndef __SX__
IF (ltimer) CALL timer_start(timer_div)
#endif
!$OMP PARALLEL

! The special treatment of 2D fields is essential for efficiency on the NEC

SELECT CASE (i_cell_type)

CASE (3) ! (i_cell_type == 3)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk)
  DO jb = i_startblk, i_endblk

    CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, rl_start, rl_end)

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
        ! ptr_patch%grid%cells%edge_orientation)
        !
        ! Distinghuish: case of a land cell (put div to zero), and
        ! cases where one of the edges are boundary or land
        ! (put corresponding velocity to zero).
        ! sea, sea_boundary, boundary (edges only), land_boundary, land =
        !  -2,      -1,         0,                  1,             2

         IF ( ptr_patch%patch_oce%lsm_oce_c(jc,jk,jb) > sea_boundary ) THEN
           div_vec_c(jc,jk,jb) = 0.0_wp
         ELSE
!            IF ( ptr_patch%patch_oce%lsm_oce_e(iidx(jc,jb,1),jk,iblk(jc,jb,1)) >= boundary ) &
!              &  vec_e(iidx(jc,jb,1),jk,iblk(jc,jb,1)) = 0.0_wp
!            IF ( ptr_patch%patch_oce%lsm_oce_e(iidx(jc,jb,2),jk,iblk(jc,jb,2)) >= boundary ) &
!              &  vec_e(iidx(jc,jb,1),jk,iblk(jc,jb,2)) = 0.0_wp
!            IF ( ptr_patch%patch_oce%lsm_oce_e(iidx(jc,jb,3),jk,iblk(jc,jb,3)) >= boundary ) &
!              &  vec_e(iidx(jc,jb,1),jk,iblk(jc,jb,3)) = 0.0_wp

          div_vec_c(jc,jk,jb) =  &
            vec_e(iidx(jc,jb,1),jk,iblk(jc,jb,1)) * ptr_patch%patch_oce%geofac_div(jc,1,jb) + &
            vec_e(iidx(jc,jb,2),jk,iblk(jc,jb,2)) * ptr_patch%patch_oce%geofac_div(jc,2,jb) + &
            vec_e(iidx(jc,jb,3),jk,iblk(jc,jb,3)) * ptr_patch%patch_oce%geofac_div(jc,3,jb)

        ENDIF

      END DO
    END DO
  END DO
!$OMP END DO

!    do jb=1,1
!    do jc=1,10
!      write(*,*)'vec_e, div_c, geofac:', jc,jb, &
!        &   vec_e(jc,1,jb), div_vec_c(jc,1,jb), (ptr_patch%patch_oce%geofac_div(jc,jk,jb),jk=1,2)
!    enddo
!    enddo
!    do jb=1,1
!    do jc=1,10
!      write(*,*)'edge_len:', jc,jb, &
!        &   (ptr_patch%edges%primal_edge_length(iidx(jc,jb,jk),iblk(jc,jb,jk)),jk=1,3), &
!        &    ptr_patch%cells%area(jc,jb)
!    enddo
!    enddo

CASE (6) ! (i_cell_type == 6)

  ! no grid refinement in hexagonal model
  nblks_c   = ptr_patch%nblks_int_c
  npromz_c  = ptr_patch%npromz_int_c

!$OMP DO PRIVATE(jb,nlen,je,jc,jk)
  DO jb = 1, nblks_c

    IF (jb /= nblks_c) THEN
      nlen = nproma
    ELSE
      nlen = npromz_c
    ENDIF

    div_vec_c(1:nlen,slev:elev,jb) = 0.0_wp

    DO je = 1, i_cell_type

      DO jk = slev, elev
        DO jc = 1, nlen
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
        !
          IF (je > ptr_patch%cells%num_edges(jc,jb)) CYCLE

          div_vec_c(jc,jk,jb) =  div_vec_c(jc,jk,jb) +   &
            &    vec_e(iidx(jc,jb,je),jk,iblk(jc,jb,je)) &
            &  * ptr_patch%patch_oce%geofac_div(jc,je,jb)

        END DO
      END DO
    END DO

  END DO
END SELECT
!$OMP END PARALLEL
#ifndef __SX__
IF (ltimer) CALL timer_stop(timer_div)
#endif
END SUBROUTINE div_oce

!-------------------------------------------------------------------------
!

!>
!! Computes the discrete rotation at vertices in presence of boundaries as in the ocean setting.
!!
!! Computes in presence of boundaries the discrete rotation at vertices
!! of triangle cells (centers of dual grid cells) from a vector field
!! given by its components in the directions normal to triangle edges and
!! takes the presence of boundaries into account.
!! Boundary condition:
!! inviscid case: normal velocity component is zero
!! viscous case: normal and tangential velocity are zero
!!
!! @par Revision History
!! Developed and tested  by L.Bonaventura  (2002-4).
!! Adapted to new grid and patch structure by P. Korn (2005).
!! Some further changes by L. Bonaventura August 2005.
!! Implementation of both boundary conditions by P. Korn (2006)
!! Modifications by P. Korn, MPI-M(2007-2)
!! -Switch fom array arguments to pointers
!! Modifications by Almut Gassmann, MPI-M (2007-04-20)
!! - abandon grid for the sake of patch
!! - DUMMY EDGE COMES NOW AS THE LAST EDGE AFTER THE BOUNDARY EDGES
!! Modifications by Peter Korn, MPI-M (2010-04)
!! Modifications by Stephan Lorenz, MPI-M (2010-06)
!!
SUBROUTINE rot_vertex_ocean( u_vec_e, v_vec_e, ptr_patch, rot_vec_v,  &
  &                          opt_slev, opt_elev, opt_rlstart, opt_rlend )

! input:  lives on edges (velocity points)
! output: lives on vertices of triangle
!

INTEGER, PARAMETER :: ino_dual_edges = 6

!
!  patch on which computation is performed
!
TYPE(t_patch), TARGET, INTENT(in) :: ptr_patch
!
!  edge based variable of which rotation is computed

!   velocity normal to edges
REAL(wp), INTENT(in) ::  &
  &  u_vec_e(:,:,:) ! dim: (nproma,n_zlev,nblks_e)
!
!velocity tangential to edges (only required for inviscid case)
REAL(wp), INTENT(in) ::  &
  &  v_vec_e(:,:,:) ! dim: (nproma,n_zlev,nblks_e)

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_slev    ! optional vertical start level

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_elev    ! optional vertical end level

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_rlstart, opt_rlend   ! start and end values of refin_ctrl flag

!
!  vertex based variable in which rotation is stored
REAL(wp), INTENT(out) :: rot_vec_v(:,:,:) ! dim: (nproma,n_zlev,nblks_v)

REAL(wp) :: ztmp
REAL(wp) :: zarea_fraction

INTEGER :: slev, elev     ! vertical start and end level
INTEGER :: jv, jk, jb, jev
INTEGER :: ile, ibe, il, ib, ill
INTEGER :: ik, ikk
INTEGER :: rl_start, rl_end
INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx, i_nchdom

INTEGER :: iwet_edge_idx(ino_dual_edges)   !line indices of non-boundary edges
INTEGER :: iwet_cell_idx(2*ino_dual_edges) !line indices of cells adjacent to non-boundary edges
                                           !double counting allowed
INTEGER :: iwet_edge_blk(ino_dual_edges)   !block indices of non-boundary edges
INTEGER :: iwet_cell_blk(2*ino_dual_edges) !block indices of cells adjacent to non-boundary edges
                                           !double counting allowed
INTEGER :: iwet_cell_ctr
INTEGER :: icell_idx_1, icell_blk_1
INTEGER :: icell_idx_2, icell_blk_2

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
  elev = n_zlev
END IF

IF ( PRESENT(opt_rlstart) ) THEN
  IF (opt_rlstart == 1) THEN
    CALL finish ('mo_math_operators:rot_vertex_ocean',  &
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

! #slo# due to nag -nan compiler-option
rot_vec_v(:,:,:) = 0.0_wp ! dim: (nproma,n_zlev,nblks_v)


! values for the blocking
i_nchdom   = MAX(1,ptr_patch%n_childdom)
i_startblk = ptr_patch%verts%start_blk(rl_start,1)
i_endblk   = ptr_patch%verts%end_blk(rl_end,i_nchdom)

!
!  loop through over all patch vertices (and blocks)
!
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
  ! At the end local curl is divided by area of dual cell,
  ! if boundaries are present, not the complete area of the
  ! dual cell counts but only the 'wet' part. This wet part
  ! is accumulated in variable 'zarea_fraction'.

  DO jk = slev, elev

!$OMP PARALLEL DO SCHEDULE(runtime) DEFAULT(PRIVATE)  &
!$OMP   SHARED(u_vec_e,v_vec_e,ptr_patch,rot_vec_v,jb) FIRSTPRIVATE(jk)
    DO jv = i_startidx, i_endidx

      ztmp = 0.0_wp

      !init indices
      iwet_cell_idx(:) = 0
      iwet_edge_idx(:) = 0
      iwet_cell_blk(:) = 0
      iwet_edge_blk(:) = 0

      iwet_cell_ctr   = 0

      zarea_fraction  = 0.0_wp

      DO jev = 1, ptr_patch%verts%num_edges(jv,jb)

        !
        ! get line and block indices of edge jev around vertex jv
        !
        ile = ptr_patch%verts%edge_idx(jv,jb,jev)
        ibe = ptr_patch%verts%edge_blk(jv,jb,jev)

        !Check, if edge is sea or boundary edge and take care of dummy edge
        ! edge with indices ile, ibe is sea edge
        IF ( ptr_patch%patch_oce%lsm_oce_e(ile,jk,ibe) <= sea_boundary ) THEN

          !Distinguish the following cases
          ! edge ie_k is
          !a) ocean edge: compute as usual,
          !b) land edge: do not consider it
          !c) boundary edge take:
          !   1) Inviscid fluid: normal velocity at boundary is zero
          !   2) Viscous fluid:  normal and tangential velocity at boundary are zero
          ! sea, sea_boundary, boundary (edges only), land_boundary, land =
          !  -2,      -1,         0,                  1,             2

          !add contribution of normal velocity at edge (ile,ibe) to rotation
          ztmp = ztmp + u_vec_e(ile,jk,ibe)  &
            & * ptr_patch%edges%dual_edge_length(ile,ibe)  &
            & * ptr_patch%verts%edge_orientation(jv,jb,jev)

          !increase ctr, store edge (ile,ibe) as wet edge idx, blk
          !and store simply both of the adjacent cells as wet cells
          !the double counting of cells is handled below and only
          !if boundary edges are present in this dual cell
          iwet_cell_ctr          = iwet_cell_ctr + 1
          iwet_edge_idx(jev)     = ile
          iwet_edge_blk(jev)     = ibe
          iwet_cell_idx(2*jev)   = ptr_patch%edges%cell_idx(ile,ibe,1)
          iwet_cell_idx(2*jev-1) = ptr_patch%edges%cell_idx(ile,ibe,2)
          iwet_cell_blk(2*jev)   = ptr_patch%edges%cell_blk(ile,ibe,1)
          iwet_cell_blk(2*jev-1) = ptr_patch%edges%cell_blk(ile,ibe,2)

        ! edge with indices ile, ibe is boundary edge
        ELSE IF ( ptr_patch%patch_oce%lsm_oce_e(ile,jk,ibe) == boundary ) THEN

          !dual edge length = distance between adjacent triangle centers
          !if edge belongs to boundary only half of this length is taken into account
          !
          !add half contribution of normal velocity at edge "ie" to rotation
          ztmp = ztmp + u_vec_e(ile,jk,ibe)  &
            &  * 0.5_wp * ptr_patch%edges%dual_edge_length(ile,ibe)  &
            &           * ptr_patch%verts%edge_orientation(jv,jb,jev)

          !increase ctr, store ie_k as wet edge idx
          !and store simply both of the adjacent cells as wet cells
          !the double counting of cells is handled below and only
          !if boundary edges are present in this dual cell
          iwet_cell_ctr      = iwet_cell_ctr + 1
          iwet_edge_idx(jev) = ile
          iwet_edge_blk(jev) = ibe
          icell_idx_1        = ptr_patch%edges%cell_idx(ile,ibe,1)
          icell_idx_2        = ptr_patch%edges%cell_idx(ile,ibe,2)
          icell_blk_1        = ptr_patch%edges%cell_blk(ile,ibe,1)
          icell_blk_2        = ptr_patch%edges%cell_blk(ile,ibe,2)

          !check which one is the ocean cell and store it
          IF ( ptr_patch%patch_oce%lsm_oce_c(icell_idx_1,jk,icell_blk_1) <= sea_boundary ) THEN
            iwet_cell_idx(2*jev)  = icell_idx_1
            iwet_cell_idx(2*jev-1)= 0
            iwet_cell_blk(2*jev)  = icell_blk_1
            iwet_cell_blk(2*jev-1)= 0
          ELSE
            iwet_cell_idx(2*jev)  = icell_idx_2
            iwet_cell_idx(2*jev-1)= 0
            iwet_cell_blk(2*jev)  = icell_blk_2
            iwet_cell_blk(2*jev-1)= 0
          END IF

          !In the viscous case with no-slip boundary condition we are done,
          !because at boundary edges the normal velocity was taken into account and
          !the tangential velocity is zero.
          !For the inviscid case the tangential velocity along the boundary has to be considered.
          !This sbr is called from the sbr 'mo_sw_state%compute_sw_diag' after the tangential
          !velocity has been reconstructed, i.e. v_vec_e=vt can be used here without further computations.
          IF ( .NOT.lviscous ) THEN

            !idx of wet cell to which edge 'ie' belongs
            il = iwet_cell_idx(2*jev)
            ib = iwet_cell_blk(2*jev)

            !determine which of the three edges of cell 'il,ib' coincides with edge 'ile,ibe'
            IF ( ptr_patch%cells%edge_idx(il,ib,1) == ile .AND.  &
              &  ptr_patch%cells%edge_blk(il,ib,1) == ibe ) THEN
              ill = 1
            ELSE IF ( ptr_patch%cells%edge_idx(il,ib,2) == ile .AND.  &
              &       ptr_patch%cells%edge_blk(il,ib,2) == ibe ) THEN
              ill = 2
            ELSE IF ( ptr_patch%cells%edge_idx(il,ib,3) == ile .AND.  &
              &       ptr_patch%cells%edge_blk(il,ib,3) == ibe ) THEN
              ill = 3
            ELSE
              ill = 0
              CALL finish ('rot_vertex_ocean:',&
                'Serious error in vorticity computation: none of three edges')
            END IF

            ztmp = ztmp + v_vec_e(ile,jk,ibe)  &
              & * ptr_patch%edges%primal_edge_length(ile,ibe)&
              & * ptr_patch%cells%edge_orientation(il,ib,ill)

          END IF

        END IF

      END DO
      !
      !divide by hex/pentagon area, if all dual cells are in the ocean interior
      !divide by apropriate fraction if boundaries are involved
      IF ( iwet_cell_ctr == ptr_patch%verts%num_edges(jv,jb) ) THEN

        rot_vec_v(jv,jk,jb) = ztmp / ptr_patch%verts%dual_area(jv,jb)

      ELSE

        !eliminate double counted wet cells,
        !i.e. set index to zero, if it occurs twice
        DO ik = 1, 2*ino_dual_edges
          DO ikk = ik+1, 2*ino_dual_edges

            IF ( iwet_cell_idx(ikk) == iwet_cell_idx(ik) .AND.  &
              &  iwet_cell_blk(ikk) == iwet_cell_blk(ik) .AND.  &
              &  iwet_cell_idx(ik) /= 0 .AND. iwet_cell_blk(ik) /= 0 ) THEN

              iwet_cell_idx(ikk) = 0
              iwet_cell_blk(ikk) = 0

            ENDIF

          END DO
        END DO
        !loop over wet cells again and accumulate area of wet cell
        DO ik = 1, 2*ino_dual_edges
          IF ( iwet_cell_idx(ik) /= 0 .AND. iwet_cell_blk(ik) /= 0 ) THEN
            zarea_fraction = zarea_fraction  &
              &            + ptr_patch%cells%area(iwet_cell_idx(ik),iwet_cell_blk(ik))
          ENDIF
        END DO

        ! no division by zero in case of zero-test (#slo# 2010-06-09)
        IF (zarea_fraction == 0.0_wp) THEN
          rot_vec_v(jv,jk,jb) = 0.0_wp
        ELSE
          rot_vec_v(jv,jk,jb) = ztmp / zarea_fraction
        ENDIF

      ENDIF

    END DO
!$OMP END PARALLEL DO

  END DO

END DO

END SUBROUTINE rot_vertex_ocean

!-------------------------------------------------------------------------
!

!>
!!  Computes  laplacian of a vector field.
!!
!! input:  lives on edges (velocity points)
!! output: lives on edges
!!
!! @par Revision History
!! Developed and tested  by L.Bonaventura  (2002-4).
!! Adapted to new grid and patch structure by P. Korn (2005).
!! Modified by Th.Heinze (2006-06-20):
!! - changed u_out(ie1,jn) to u_out(j,jn) according to hint of P.Korn
!! Modifications by P. Korn, MPI-M(2007-2)
!! -Switch fom array arguments to pointers
!! Modified by P Ripodas (2007-02):
!! - include the system orientation factor in the vorticity term
!! Modified by Almut Gassmann, MPI-M (2007-04-20)
!! - abandon grid for the sake of patch
!!
SUBROUTINE nabla2_vec_ocean( u_vec_e, v_vec_e, ptr_patch, nabla2_vec_e,  &
  &                          opt_slev, opt_elev, opt_rlstart, opt_rlend )

!

!
!  patch on which computation is performed
!
TYPE(t_patch), TARGET, INTENT(in) :: ptr_patch
!
!  edge based variable of which laplacian is computed
!  (normal and tangential velocity component)
!
REAL(wp), INTENT(inout) ::  &
  &  u_vec_e(:,:,:), v_vec_e(:,:,:) ! dim: (nproma,n_zlev,nblks_e)

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_slev    ! optional vertical start level

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_elev    ! optional vertical end level

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_rlstart, opt_rlend   ! start and end values of refin_ctrl flag

!
!  edge based variable in which laplacian is stored
!
!REAL(wp), INTENT(out)   ::  &
REAL(wp), INTENT(inout) ::  &
  &  nabla2_vec_e(:,:,:)  ! dim: (nproma,n_zlev,nblks_e)

INTEGER :: slev, elev     ! vertical start and end level
INTEGER :: nblks_c
INTEGER :: je, jk, jb
INTEGER :: rl_start, rl_end
INTEGER :: rl_start_c, rl_end_c, rl_start_v, rl_end_v
INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx, i_nchdom

INTEGER, DIMENSION(nproma) ::  &
  &  ilc1, ibc1, ilc2, ibc2,   &
  &  ilv1, ibv1, ilv2, ibv2

REAL(wp) :: rorient       ! orientation of the system (v2-v1),(c2-c1)
REAL(wp) ::  &
  &  z_div_c(nproma,n_zlev,ptr_patch%nblks_c), &
  &  z_rot_v(nproma,n_zlev,ptr_patch%nblks_v)

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
  elev = n_zlev
END IF

IF ( PRESENT(opt_rlstart) ) THEN
  IF ((opt_rlstart >= 0) .AND. (opt_rlstart <= 2)) THEN
    CALL finish ('mo_math_operators:nabla2_vec_atmos',  &
          &      'opt_rlstart must not be between 0 and 2')
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

rl_start_c = rl_start/2

IF (rl_start > 0) THEN
  rl_start_v = (rl_start+1)/2
ELSE
  rl_start_v = (rl_start-1)/2
ENDIF

IF (rl_end > 0) THEN
  rl_end_c = (rl_end+1)/2
  rl_end_v = rl_end/2+1
ELSE
  rl_end_c = (rl_end-1)/2
  rl_end_v = rl_end/2-1
ENDIF

rl_end_c = MAX(min_rlcell,rl_end_c)
rl_end_v = MAX(min_rlvert,rl_end_v)

! values for the blocking
nblks_c  = ptr_patch%nblks_c

i_nchdom   = MAX(1,ptr_patch%n_childdom)
i_startblk = ptr_patch%edges%start_blk(rl_start,1)
i_endblk   = ptr_patch%edges%end_blk(rl_end,i_nchdom)

! compute divergence of vector field
! the divergence on land is set to zero by div_oce
CALL div_oce( u_vec_e, ptr_patch, z_div_c, slev, elev, &
          opt_rlstart=rl_start_c, opt_rlend=rl_end_c )


! compute rotation of vector field for the ocean
CALL rot_vertex_ocean( u_vec_e, v_vec_e, ptr_patch, z_rot_v, slev, elev, &
  &                    opt_rlstart=rl_start_v, opt_rlend=rl_end_v )

!
!  loop through all patch edges (and blocks)
!
DO jb = i_startblk, i_endblk

  CALL get_indices_e(ptr_patch, jb, i_startblk, i_endblk, &
                     i_startidx, i_endidx, rl_start, rl_end)

  DO je = i_startidx, i_endidx

    !
    !  get the line and block indices of the cells neigbouring edge je
    !
    ilc1(je) = ptr_patch%edges%cell_idx(je,jb,1)
    ibc1(je) = ptr_patch%edges%cell_blk(je,jb,1)
    ilc2(je) = ptr_patch%edges%cell_idx(je,jb,2)
    ibc2(je) = ptr_patch%edges%cell_blk(je,jb,2)

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

      rorient = ptr_patch%edges%system_orientation(je,jb)

      nabla2_vec_e(je,jk,jb) =  &
        & rorient *  &
        &   ( z_rot_v(ilc2(je),jk,ibv2(je)) - z_rot_v(ilc1(je),jk,ibv1(je)) )  &
        &   * ptr_patch%edges%inv_primal_edge_length(je,jb)  &
        & + ( z_div_c(ilv2(je),jk,ibc2(je)) - z_div_c(ilv1(je),jk,ibc1(je)) )  &
        &   * ptr_patch%edges%inv_dual_edge_length(je,jb)

    END DO

  END DO

  !set output on land to zero
  DO jk = slev, elev
    WHERE (ptr_patch%patch_oce%lsm_oce_e(1:nproma,jk,jb) >= land_boundary)
      nabla2_vec_e(1:nproma,jk,jb) = 0.0_wp
    END WHERE
  ENDDO

END DO

END SUBROUTINE nabla2_vec_ocean
!---------------------------------------------------------------------------------
!
!  
!>
!! 
!!  Calculation of total fluid thickness at cell centers and surface elevation at 
!!  cell edges from prognostic surface height at cell centers. We use height at
!!  old timelevel "n"
!! 
!! @par Revision History
!! Developed  by  Peter Korn, MPI-M (2010).
!! 
SUBROUTINE height_related_quantities( p_patch, p_os)
!
! Patch on which computation is performed
TYPE(t_patch), TARGET, INTENT(in) :: p_patch
!
! Type containing ocean state
TYPE(t_hydro_ocean_state), TARGET :: p_os
!
!  local variables
!

INTEGER :: i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c
INTEGER :: rl_start_c, rl_end_c
INTEGER :: i_startblk_e, i_endblk_e, i_startidx_e, i_endidx_e
INTEGER :: rl_start_e, rl_end_e
INTEGER :: jc, jb, je
INTEGER :: il_c1, ib_c1, il_c2, ib_c2

REAL(wp) :: z_dist_e_c1, z_dist_e_c2
INTEGER ::            AVERAGING       =1
INTEGER, PARAMETER :: DISTANCE_WEIGHT =1
INTEGER, PARAMETER :: UPWIND          =2
TYPE(t_cartesian_coordinates)    :: cv_c1_e0, cv_c2_e0
TYPE(t_cartesian_coordinates)    :: cc_c1, cc_c2, cc_e0
CHARACTER(len=max_char_length), PARAMETER :: &
       & routine = ('mo_oce_AB_timestepping_mimetic:height_related_quantities')
!-------------------------------------------------------------------------------
!CALL message (TRIM(routine), 'start')

rl_start_c = 1
rl_end_c = min_rlcell
i_startblk_c = p_patch%cells%start_blk(rl_start_c,1)
i_endblk_c   = p_patch%cells%end_blk(rl_end_c,1)

rl_start_e = 1
rl_end_e = min_rledge
i_startblk_e = p_patch%edges%start_blk(rl_start_e,1)
i_endblk_e   = p_patch%edges%end_blk(rl_end_e,1)

!Step 1: calculate cell-located variables for 2D and 3D case
!For 3D and for SWE thick_c contains thickness of fluid column
IF ( iswm_oce == 1 ) THEN  !  SWM

  DO jb = i_startblk_c, i_endblk_c
    CALL get_indices_c(p_patch, jb, i_startblk_c, i_endblk_c,&
                     & i_startidx_c, i_endidx_c, rl_start_c, rl_end_c)
#ifdef __SX__
!CDIR UNROLL=6
#endif
    !calculate for each fluid colum the total depth, i.e. 
    !from bottom boundary to surface height, i.e. using individual bathymetry for SWM
    DO jc = i_startidx_c, i_endidx_c
      IF(p_patch%patch_oce%lsm_oce_c(jc,1,jb) <= sea_boundary ) THEN
        p_os%p_diag%thick_c(jc,jb) = p_os%p_prog(nold(1))%h(jc,jb)&
                                  &- p_patch%patch_oce%bathymetry_c(jc,jb)
      ELSE
        p_os%p_diag%thick_c(jc,jb) = 0.0_wp
      ENDIF 
    END DO
  END DO

ELSEIF(iswm_oce /= 1 )THEN

  DO jb = i_startblk_c, i_endblk_c
   CALL get_indices_c(p_patch, jb, i_startblk_c, i_endblk_c,&
                    & i_startidx_c, i_endidx_c, rl_start_c, rl_end_c)
#ifdef __SX__
!CDIR UNROLL=6
#endif
    !calculate for each fluid colum the total depth, i.e. 
    !from bottom boundary to surface height 
    !the bottom boundary is zlev_i(dolic+1) since zlev_i(1)=0 (air-sea-boundary at h=0
    DO jc = i_startidx_c, i_endidx_c
      IF ( p_patch%patch_oce%lsm_oce_c(jc,1,jb) <= sea_boundary ) THEN 
        p_os%p_diag%thick_c(jc,jb) = p_os%p_prog(nold(1))%h(jc,jb)&
        &                          + p_patch%patch_oce%zlev_i(p_patch%patch_oce%dolic_c(jc,jb)+1)
      ELSE
        p_os%p_diag%thick_c(jc,jb) = 0.0_wp
      ENDIF
    END DO
  END DO
ENDIF!IF ( iswm_oce == 1 )


!Step 2: calculate edge-located variables for 2D and 3D case from respective cell variables
!For SWE and for 3D: thick_e = thickness of fluid column at edges
!         h_e     = surface elevation at edges, without depth of first layer 
IF ( iswm_oce == 1 ) THEN  !  SWM
SELECT CASE(AVERAGING)

CASE(DISTANCE_WEIGHT)

  DO jb = i_startblk_e, i_endblk_e
    CALL get_indices_e(p_patch, jb, i_startblk_e, i_endblk_e,&
                     & i_startidx_e, i_endidx_e, rl_start_e, rl_end_e)
#ifdef __SX__
!CDIR UNROLL=6
#endif
    DO je = i_startidx_e, i_endidx_e

      il_c1 = p_patch%edges%cell_idx(je,jb,1) 
      ib_c1 = p_patch%edges%cell_blk(je,jb,1)  
      il_c2 = p_patch%edges%cell_idx(je,jb,2) 
      ib_c2 = p_patch%edges%cell_blk(je,jb,2) 

      z_dist_e_c1 = 0.5_wp!p_patch%patch_oce%dist_cell2edge(je,jb,1) !0.5_wp 
      z_dist_e_c2 = 0.5_wp!p_patch%patch_oce%dist_cell2edge(je,jb,2) !0.5_wp 

      !z_dist_e_c1=p_patch%edges%edge_cell_length(je,jb,1)
      !z_dist_e_c2=p_patch%edges%edge_cell_length(je,jb,2)

      IF ( p_patch%patch_oce%lsm_oce_e(je,1,jb) <= sea_boundary ) THEN   
        p_os%p_diag%thick_e(je,jb) = ( z_dist_e_c1*p_os%p_diag%thick_c(il_c1,ib_c1)&
                              &+   z_dist_e_c2*p_os%p_diag%thick_c(il_c2,ib_c2) )&
                              &/(z_dist_e_c1+z_dist_e_c2)

        p_os%p_diag%h_e(je,jb) = ( z_dist_e_c1*p_os%p_prog(nold(1))%h(il_c1,ib_c1)&
                              &+   z_dist_e_c2*p_os%p_prog(nold(1))%h(il_c2,ib_c2) )&
                              &/(z_dist_e_c1+z_dist_e_c2)
        !write(*,*)'height_e',je,jb, p_os%p_diag%h_e(je,jb), p_os%p_prog(nold(1))%h(il_c1,ib_c1),p_os%p_prog(nold(1))%h(il_c2,ib_c2) 
      ELSE
        p_os%p_diag%h_e(je,jb) = 0.0_wp
      ENDIF 
    END DO
  END DO

CASE(upwind)

  DO jb = i_startblk_e, i_endblk_e
    CALL get_indices_e(p_patch, jb, i_startblk_e, i_endblk_e,&
                     & i_startidx_e, i_endidx_e, rl_start_e, rl_end_e)
#ifdef __SX__
!CDIR UNROLL=6
#endif
    DO je = i_startidx_e, i_endidx_e

      il_c1 = p_patch%edges%cell_idx(je,jb,1) 
      ib_c1 = p_patch%edges%cell_blk(je,jb,1)  
      il_c2 = p_patch%edges%cell_idx(je,jb,2) 
      ib_c2 = p_patch%edges%cell_blk(je,jb,2) 


      IF( p_patch%patch_oce%lsm_oce_e(je,1,jb) <= sea_boundary ) THEN

        IF(p_os%p_prog(nold(1))%vn(je,1,jb)>0)THEN
          p_os%p_diag%thick_e(je,jb) = p_os%p_diag%thick_c(il_c1,ib_c1)
          p_os%p_diag%h_e(je,jb)     = p_os%p_prog(nold(1))%h(il_c1,ib_c1)
        ELSEIF(p_os%p_prog(nold(1))%vn(je,1,jb)<=0)THEN
          p_os%p_diag%thick_e(je,jb) = p_os%p_diag%thick_c(il_c2,ib_c2)
          p_os%p_diag%h_e(je,jb)     = p_os%p_prog(nold(1))%h(il_c2,ib_c2)
        ENDIF
        !write(*,*)'height_e',je,jb, p_os%p_diag%h_e(je,jb), p_os%p_prog(nold(1))%h(il_c1,ib_c1),p_os%p_prog(nold(1))%h(il_c2,ib_c2) 
      ELSE
        p_os%p_diag%h_e(je,jb) = 0.0_wp
      ENDIF
    END DO
  END DO
END SELECT



ELSEIF(iswm_oce /= 1 )THEN


SELECT CASE(AVERAGING)

CASE(DISTANCE_WEIGHT)
  DO jb = i_startblk_e, i_endblk_e
    CALL get_indices_e(p_patch, jb, i_startblk_e, i_endblk_e,&
                     & i_startidx_e, i_endidx_e, rl_start_e, rl_end_e)
#ifdef __SX__
!CDIR UNROLL=6
#endif
    DO je = i_startidx_e, i_endidx_e

      il_c1 = p_patch%edges%cell_idx(je,jb,1) 
      ib_c1 = p_patch%edges%cell_blk(je,jb,1)  
      il_c2 = p_patch%edges%cell_idx(je,jb,2) 
      ib_c2 = p_patch%edges%cell_blk(je,jb,2) 

      z_dist_e_c1 = 0.5_wp* p_patch%patch_oce%dist_cell2edge(je,jb,1)
      z_dist_e_c2 = 0.5_wp* p_patch%patch_oce%dist_cell2edge(je,jb,2)
      !z_dist_e_c1 = p_patch%edges%edge_cell_length(je,jb,1)
      !z_dist_e_c2 = p_patch%edges%edge_cell_length(je,jb,2)

      IF( p_patch%patch_oce%lsm_oce_e(je,1,jb) <= sea_boundary ) THEN

        p_os%p_diag%thick_e(je,jb) = ( z_dist_e_c1*p_os%p_diag%thick_c(il_c1,ib_c1)&
                              &+   z_dist_e_c2*p_os%p_diag%thick_c(il_c2,ib_c2) )&
                              &/(z_dist_e_c1+z_dist_e_c2)

        !  P.K.: Actually h_e is just surface elevation at edges without depth of first layer.
        !It might make sense to include depth of first layer.
        p_os%p_diag%h_e(je,jb) = ( z_dist_e_c1*p_os%p_prog(nold(1))%h(il_c1,ib_c1)&
                              &+   z_dist_e_c2*p_os%p_prog(nold(1))%h(il_c2,ib_c2) )&
                              &/(z_dist_e_c1+z_dist_e_c2)
        !write(*,*)'height_e',je,jb, p_os%p_diag%h_e(je,jb), p_os%p_prog(nold(1))%h(il_c1,ib_c1),p_os%p_prog(nold(1))%h(il_c2,ib_c2) 
      ELSE
        p_os%p_diag%h_e(je,jb) = 0.0_wp
      ENDIF
    END DO
  END DO

CASE(upwind)

  DO jb = i_startblk_e, i_endblk_e
    CALL get_indices_e(p_patch, jb, i_startblk_e, i_endblk_e,&
                     & i_startidx_e, i_endidx_e, rl_start_e, rl_end_e)
#ifdef __SX__
!CDIR UNROLL=6
#endif
    DO je = i_startidx_e, i_endidx_e

      il_c1 = p_patch%edges%cell_idx(je,jb,1) 
      ib_c1 = p_patch%edges%cell_blk(je,jb,1)  
      il_c2 = p_patch%edges%cell_idx(je,jb,2) 
      ib_c2 = p_patch%edges%cell_blk(je,jb,2) 


      IF( p_patch%patch_oce%lsm_oce_e(je,1,jb) <= sea_boundary ) THEN

        IF(p_os%p_prog(nold(1))%vn(je,1,jb)>0)THEN
          p_os%p_diag%thick_e(je,jb) = p_os%p_diag%thick_c(il_c1,ib_c1)
        ELSEIF(p_os%p_prog(nold(1))%vn(je,1,jb)<=0)THEN

          !  P.K.: Actually h_e is just surface elevation at edges without depth of first layer.
          !It might make sense to include depth of first layer.
          p_os%p_diag%h_e(je,jb) = p_os%p_prog(nold(1))%h(il_c2,ib_c2)
        ENDIF
        !write(*,*)'height_e',je,jb, p_os%p_diag%h_e(je,jb), p_os%p_prog(nold(1))%h(il_c1,ib_c1),p_os%p_prog(nold(1))%h(il_c2,ib_c2) 
      ELSE
        p_os%p_diag%h_e(je,jb) = 0.0_wp
      ENDIF
    END DO
  END DO
END SELECT
ENDIF
! write(*,*)'max/min thick_c:',maxval(p_os%p_diag%thick_c),minval(p_os%p_diag%thick_c) 
! write(*,*)'max/min h_c:',maxval(p_os%p_prog(nold(1))%h),minval(p_os%p_prog(nold(1))%h) 
! write(*,*)'max/min bath_c:',maxval(p_patch%patch_oce%bathymetry_c),minval(p_patch%patch_oce%bathymetry_c) 
! write(*,*)'max/min thick_e:',maxval(p_os%p_diag%thick_e),minval(p_os%p_diag%thick_e) 
! write(*,*)'max/min h_e:',maxval(p_os%p_diag%h_e),minval(p_os%p_diag%h_e) 
! !CALL message (TRIM(routine), 'end')        
END SUBROUTINE height_related_quantities

END MODULE mo_oce_math_operators

