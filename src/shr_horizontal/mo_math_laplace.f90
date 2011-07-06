!>
!!   Contains the implementation of the nabla mathematical operators.
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
MODULE mo_math_laplace
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
USE mo_interpolation,       ONLY: t_int_state, edges2verts_scalar, verts2edges_scalar
USE mo_model_domain,        ONLY: t_patch
USE mo_model_domain_import, ONLY: l_limited_area
USE mo_parallel_configuration,  ONLY: nproma
USE mo_exception,           ONLY: finish
USE mo_loopindices,         ONLY: get_indices_c, get_indices_e
USE mo_sync,                ONLY: SYNC_C, SYNC_E, SYNC_V, sync_patch_array

USE mo_math_gradients
USE mo_math_divrot


IMPLICIT NONE

PRIVATE

CHARACTER(len=*), PARAMETER :: version = '$Id$'


PUBLIC :: nabla2_vec
PUBLIC :: nabla2_scalar, nabla2_scalar_avg
PUBLIC :: nabla4_vec
PUBLIC :: nabla4_scalar
PUBLIC :: nabla6_vec
PUBLIC :: directional_laplace

INTERFACE nabla2_vec

  MODULE PROCEDURE nabla2_vec_atmos

END INTERFACE

CONTAINS


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
!! Modified by Almut Gassmann (2007-04-20)
!! - abandon grid for the sake of patch
!!
SUBROUTINE nabla2_vec_atmos( vec_e, ptr_patch, ptr_int, nabla2_vec_e, &
  &                          opt_slev, opt_elev, opt_rlstart, opt_rlend )

!
!  patch on which computation is performed
!
TYPE(t_patch), TARGET, INTENT(in) :: ptr_patch

! Interpolation state
TYPE(t_int_state), INTENT(in)     :: ptr_int
!
!  edge based variable of which laplacian is computed
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
!  edge based variable in which laplacian is stored
!
!REAL(wp), INTENT(out) ::  &
REAL(wp), INTENT(inout) ::  &
  &  nabla2_vec_e(:,:,:) ! dim: (nproma,nlev,nblks_e)

INTEGER :: slev, elev     ! vertical start and end level
INTEGER :: je, jk, jb
INTEGER :: rl_start, rl_end
INTEGER :: rl_start_c, rl_end_c, rl_start_v, rl_end_v
INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx, i_nchdom
INTEGER :: nlen, nblks_e, npromz_e

REAL(wp) ::  &
  &  z_div_c(nproma,ptr_patch%nlev,ptr_patch%nblks_c),  &
  &  z_rot_v(nproma,ptr_patch%nlev,ptr_patch%nblks_v),  &
  &  z_rot_e(nproma,ptr_patch%nlev,ptr_patch%nblks_e)

INTEGER,  DIMENSION(:,:,:),   POINTER :: icidx, icblk, ividx, ivblk

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

icidx => ptr_patch%edges%cell_idx
icblk => ptr_patch%edges%cell_blk
ividx => ptr_patch%edges%vertex_idx
ivblk => ptr_patch%edges%vertex_blk

i_nchdom   = MAX(1,ptr_patch%n_childdom)
i_startblk = ptr_patch%edges%start_blk(rl_start,1)
i_endblk   = ptr_patch%edges%end_blk(rl_end,i_nchdom)

! Initialization of unused elements of nabla2_vec_e
DO jb = 1, i_startblk
  nabla2_vec_e(:,:,jb) = 0._wp
ENDDO
DO jb = i_endblk, ptr_patch%nblks_e
  nabla2_vec_e(:,:,jb) = 0._wp
ENDDO

! compute divergence of vector field
CALL div( vec_e, ptr_patch, ptr_int, z_div_c, slev, elev, &
          opt_rlstart=rl_start_c, opt_rlend=rl_end_c )

!
!  loop through over all patch edges (and blocks)
!

! The special treatment of 2D fields is essential for efficiency on the NEC

SELECT CASE (ptr_patch%cell_type)

CASE (3) ! (cell_type == 3)

  ! compute rotation of vector field
  CALL rot_vertex( vec_e, ptr_patch, ptr_int, z_rot_v, slev, elev, &
                   opt_rlstart=rl_start_v, opt_rlend=rl_end_v)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,je,jk)
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

        nabla2_vec_e(je,jk,jb) =  &
          &   ptr_patch%edges%system_orientation(je,jb) *  &
          &   ( z_rot_v(ividx(je,jb,2),jk,ivblk(je,jb,2))  &
          &   - z_rot_v(ividx(je,jb,1),jk,ivblk(je,jb,1)) )  &
          &   * ptr_patch%edges%inv_primal_edge_length(je,jb)  &
          & + ( z_div_c(icidx(je,jb,2),jk,icblk(je,jb,2))    &
          &   - z_div_c(icidx(je,jb,1),jk,icblk(je,jb,1)) )  &
          &   * ptr_patch%edges%inv_dual_edge_length(je,jb)

      END DO
    END DO

  END DO
!$OMP END DO
!$OMP END PARALLEL

CASE (6) ! (cell_type == 6)

  ! compute rotation of vector field
  CALL rot_vertex( vec_e, ptr_patch, ptr_int, z_rot_v, slev, elev)
  CALL sync_patch_array(SYNC_V, ptr_patch, z_rot_v)

  ! compute rhombus vorticities
  CALL verts2edges_scalar(z_rot_v, ptr_patch, ptr_int%tria_aw_rhom, &
                          z_rot_e, slev, elev)
  CALL sync_patch_array(SYNC_E, ptr_patch, z_rot_e)

  ! average the edge values to the vertex values
  CALL edges2verts_scalar(z_rot_e, ptr_patch, ptr_int%e_1o3_v, z_rot_v, &
                          slev, elev)

  ! no grid refinement in hexagonal model
  nblks_e   = ptr_patch%nblks_int_e
  npromz_e  = ptr_patch%npromz_int_e

!$OMP PARALLEL
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

        nabla2_vec_e(je,jk,jb) =  &
          & - ( z_rot_v(ividx(je,jb,2),jk,ivblk(je,jb,2))     &
          & -   z_rot_v(ividx(je,jb,1),jk,ivblk(je,jb,1)) )   &
          &   * ptr_patch%edges%inv_primal_edge_length(je,jb) &
          & + ptr_patch%edges%system_orientation(je,jb) *     &
          &   ( z_div_c(icidx(je,jb,2),jk,icblk(je,jb,2))     &
          & -   z_div_c(icidx(je,jb,1),jk,icblk(je,jb,1)) )     &
          &   * ptr_patch%edges%inv_dual_edge_length(je,jb)

      END DO

    END DO

  END DO
!$OMP END DO
!$OMP END PARALLEL
END SELECT

END SUBROUTINE nabla2_vec_atmos


!-------------------------------------------------------------------------
!

!>
!! Computes biharmonic laplacian @f$\nabla ^4@f$ of a vector field without boundaries as used in atmospheric model.
!!
!! Computes biharmonic laplacian @f$\nabla ^4@f$ of a vector field without boundaries as used in atmospheric model.
!! input:  lives on edges (velocity points)
!! output: lives on edges
!!
!! @par Revision History
!! Developed and tested  by L.Bonaventura  (2002-4).
!! Adapted to new grid and patch structure by P. Korn (2005).
!! Modifications by P. Korn, MPI-M(2007-2)
!! -Switch fom array arguments to pointers
!! Modification by Almut Gassmann, MPI-M (2007-04-20)
!! - abandon grid for the sake of patch
!! Modified by Marco Giorgetta, MPI-M (2009-02-26)
!! - replaced nlev_ocean by nlev
!!
SUBROUTINE nabla4_vec( vec_e, ptr_patch, ptr_int, nabla4_vec_e, &
  &                    opt_nabla2, opt_slev, opt_elev, opt_rlstart, opt_rlend )

!
!  patch on which computation is performed
!
TYPE(t_patch), TARGET, INTENT(in) :: ptr_patch

! Interpolation state
TYPE(t_int_state), INTENT(in)     :: ptr_int
!
!  edge based variable of which laplacian is computed
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
!  edge based variable in which laplacian is stored
!
!REAL(wp), INTENT(out) ::  &
REAL(wp), INTENT(inout) ::  &
  &  nabla4_vec_e(:,:,:) ! dim: (nproma,nlev,nblks_e)

! Optional argument for passing nabla2 to the calling program
! (to avoid double computation for Smagorinsky diffusion and nest boundary diffusion)
REAL(wp), INTENT(inout), TARGET, OPTIONAL  ::  &
  &  opt_nabla2(:,:,:) ! dim: (nproma,nlev,nblks_e)

INTEGER :: slev, elev     ! vertical start and end level
INTEGER :: rl_start, rl_end
INTEGER :: rl_start_s1, rl_end_s1

REAL(wp), ALLOCATABLE, TARGET :: z_nabla2_vec_e(:,:,:) ! dim: (nproma,nlev,ptr_patch%nblks_e)
REAL(wp), POINTER :: p_nabla2(:,:,:)

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
  IF ((opt_rlstart >= 0) .AND. (opt_rlstart <= 4)) THEN
    CALL finish ('mo_math_operators:nabla4_vec',  &
          &      'opt_rlstart must not be between 0 and 4')
  ENDIF
  rl_start = opt_rlstart
ELSE
  rl_start = 5
END IF
IF ( PRESENT(opt_rlend) ) THEN
  rl_end = opt_rlend
ELSE
  rl_end = min_rledge
END IF

IF (rl_start > 0) THEN
  rl_start_s1 = rl_start - 2
ELSE
  rl_start_s1 = rl_start + 2
ENDIF

IF (rl_end > 0) THEN
  rl_end_s1 = rl_end + 2
ELSE
  rl_end_s1 = rl_end - 2
ENDIF

rl_end_s1 = MAX(min_rledge,rl_end_s1)

IF (PRESENT(opt_nabla2) ) THEN
  p_nabla2 => opt_nabla2
ELSE
  ALLOCATE (z_nabla2_vec_e(nproma,ptr_patch%nlev,ptr_patch%nblks_e))
  p_nabla2 => z_nabla2_vec_e
ENDIF

!
! apply second order Laplacian twice
!
CALL nabla2_vec( vec_e, ptr_patch, ptr_int, p_nabla2,  &
  &              slev, elev, opt_rlstart=rl_start_s1, opt_rlend=rl_end_s1 )

CALL sync_patch_array(SYNC_E, ptr_patch, p_nabla2)

CALL nabla2_vec( p_nabla2, ptr_patch, ptr_int, nabla4_vec_e,  &
  &              slev, elev, opt_rlstart=rl_start, opt_rlend=rl_end )

IF (.NOT. PRESENT(opt_nabla2) ) THEN
  DEALLOCATE (z_nabla2_vec_e)
ENDIF

END SUBROUTINE nabla4_vec
!-----------------------------------------------------------------------

!-------------------------------------------------------------------------
!

!>
!!  Computes triharmonic laplacian @f$\nabla ^6 @f$ of a vector field without boundaries as used in atmospheric model.
!!
!!  Computes triharmonic laplacian @f$\nabla ^6 @f$ of a vector field without boundaries as used in atmospheric model.
!! input:  lives on edges (velocity points)
!! output: lives on edges
!!
!! @par Revision History
!! Developed and tested  by L.Bonaventura  (2002-4).
!! Adapted to new grid and patch structure by P. Korn (2005).
!! Modifications by P. Korn, MPI-M(2007-2)
!! -Switch fom array arguments to pointers
!! Modifications by Almut Gassmann, MPI-M (2007-04-20)
!! - abandon grid for the sake of patch
!!
SUBROUTINE nabla6_vec( vec_e, ptr_patch, ptr_int, nabla6_vec_e, &
  &                    opt_slev, opt_elev, opt_rlstart, opt_rlend )

!
!  patch on which computation is performed
!
TYPE(t_patch), TARGET, INTENT(in) :: ptr_patch

! Interpolation state
TYPE(t_int_state), INTENT(in)     :: ptr_int
!
!  edge based variable of which laplacian is computed
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
!  edge based variable in which laplacian is stored
!
!REAL(wp), INTENT(out) ::  &
REAL(wp), INTENT(inout) ::  &
  &  nabla6_vec_e(:,:,:) ! dim: (nproma,nlev,nblks_e)

INTEGER :: slev, elev     ! vertical start and end level
INTEGER :: rl_start, rl_end
INTEGER :: rl_start_s1, rl_end_s1, rl_start_s2, rl_end_s2

REAL(wp) :: &
  &  z_nabla2_vec_e(nproma,ptr_patch%nlev,ptr_patch%nblks_e),  &
  &  z_nabla4_vec_e(nproma,ptr_patch%nlev,ptr_patch%nblks_e)

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
  IF ((opt_rlstart >= 0) .AND. (opt_rlstart <= 6)) THEN
    CALL finish ('mo_math_operators:nabla6_vec',  &
          &      'opt_rlstart must not be between 0 and 6')
  ENDIF
  rl_start = opt_rlstart
ELSE
  rl_start = 7
END IF
IF ( PRESENT(opt_rlend) ) THEN
  rl_end = opt_rlend
ELSE
  rl_end = min_rledge
END IF

IF (rl_start > 0) THEN
  rl_start_s1 = rl_start - 2
  rl_start_s2 = rl_start - 4
ELSE
  rl_start_s1 = rl_start + 2
  rl_start_s2 = rl_start + 4
ENDIF

IF (rl_end > 0) THEN
  rl_end_s1 = rl_end + 2
  rl_end_s2 = rl_end + 4
ELSE
  rl_end_s1 = rl_end - 2
  rl_end_s2 = rl_end - 4
ENDIF

rl_end_s1 = MAX(min_rledge,rl_end_s1)
rl_end_s2 = MAX(min_rledge,rl_end_s2)

!
! apply second order Laplacian three times
!
CALL nabla2_vec( vec_e, ptr_patch, ptr_int, z_nabla2_vec_e,  &
  &              slev, elev, opt_rlstart=rl_start_s2, opt_rlend=rl_end_s2 )

CALL sync_patch_array(SYNC_E, ptr_patch, z_nabla2_vec_e)

CALL nabla2_vec( z_nabla2_vec_e, ptr_patch, ptr_int, z_nabla4_vec_e,  &
  &              slev, elev, opt_rlstart=rl_start_s1, opt_rlend=rl_end_s1 )

CALL sync_patch_array(SYNC_E, ptr_patch, z_nabla4_vec_e)

CALL nabla2_vec( z_nabla4_vec_e, ptr_patch, ptr_int, nabla6_vec_e,  &
  &              slev, elev, opt_rlstart=rl_start, opt_rlend=rl_end )

END SUBROUTINE nabla6_vec
!-----------------------------------------------------------------------

!-------------------------------------------------------------------------
!

!>
!!  Computes laplacian @f$\nabla ^2 @f$ of a scalar field.
!!
!! input:  lives on cells (mass points)
!! output: lives on cells
!!
!! @par Revision History
!! Developed and tested  by L.Bonaventura  (2002-4).
!! Adapted to new grid and patch structure by P.Korn (2005).
!! Derived from nabla4_scalar by Th.Heinze (2006-07-24).
!! Modifications by P. Korn, MPI-M(2007-2)
!! -Switch fom array arguments to pointers
!! Modifications by Almut Gassmann, MPI-M (2007-04-20)
!! -abandon grid for the sake of patch
!!
SUBROUTINE nabla2_scalar( psi_c, ptr_patch, ptr_int, nabla2_psi_c, &
  &                       opt_slev, opt_elev, opt_rlstart, opt_rlend )

!
!  patch on which computation is performed
!
TYPE(t_patch), TARGET, INTENT(in) :: ptr_patch

! Interpolation state
TYPE(t_int_state), INTENT(in)     :: ptr_int
!
!  cells based variable of which biharmonic laplacian is computed
!
REAL(wp), INTENT(in) ::  &
  &  psi_c(:,:,:) ! dim: (nproma,nlev,nblks_c)

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_slev    ! optional vertical start level

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_elev    ! optional vertical end level

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_rlstart, opt_rlend   ! start and end values of refin_ctrl flag

!
!  cell based variable in which biharmonic laplacian is stored
!
!REAL(wp), INTENT(out) ::  &
REAL(wp), INTENT(inout) ::  &
  &  nabla2_psi_c(:,:,:) ! dim: (nproma,nlev,nblks_c)

INTEGER :: slev, elev     ! vertical start and end level
INTEGER :: rl_start, rl_end
INTEGER :: jb, jc, jk, i_startblk, i_endblk, i_startidx, i_endidx, i_nchdom
INTEGER :: nlev                !< number of full levels

REAL(wp) ::  &
  &  z_grad_fd_norm_e(nproma,ptr_patch%nlev,ptr_patch%nblks_e)

INTEGER,  DIMENSION(:,:,:),   POINTER :: iidx, iblk

!-----------------------------------------------------------------------

! number of vertical levels
nlev = ptr_patch%nlev

! check optional arguments
IF ( PRESENT(opt_slev) ) THEN
  slev = opt_slev
ELSE
  slev = 1
END IF
IF ( PRESENT(opt_elev) ) THEN
  elev = opt_elev
ELSE
  elev = nlev
END IF

IF ( PRESENT(opt_rlstart) ) THEN
  IF ((opt_rlstart >= 0) .AND. (opt_rlstart <= 1)) THEN
    CALL finish ('mo_math_operators:nabla2_scalar',  &
          &      'opt_rlstart must not be between 0 and 1')
  ENDIF
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


! The special treatment of 2D fields is essential for efficiency on the NEC

SELECT CASE (ptr_patch%cell_type)

CASE (3) ! (cell_type == 3)

IF (slev == elev) THEN
  jk = slev
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc)
  DO jb = i_startblk, i_endblk

    IF (jb == i_startblk) THEN
      i_startidx = ptr_patch%cells%start_idx(rl_start,1)
      i_endidx   = nproma
      IF (jb == i_endblk) i_endidx = ptr_patch%cells%end_idx(rl_end,i_nchdom)
    ELSE IF (jb == i_endblk) THEN
      i_startidx = 1
      i_endidx   = ptr_patch%cells%end_idx(rl_end,i_nchdom)
    ELSE
      i_startidx = 1
      i_endidx   = nproma
    ENDIF

    DO jc = i_startidx, i_endidx

      !
      !  calculate div(grad) in one step
      !
      nabla2_psi_c(jc,jk,jb) =  &
        &    psi_c(jc,jk,jb)                       * ptr_int%geofac_n2s(jc,1,jb) &
        &  + psi_c(iidx(jc,jb,1),jk,iblk(jc,jb,1)) * ptr_int%geofac_n2s(jc,2,jb) &
        &  + psi_c(iidx(jc,jb,2),jk,iblk(jc,jb,2)) * ptr_int%geofac_n2s(jc,3,jb) &
        &  + psi_c(iidx(jc,jb,3),jk,iblk(jc,jb,3)) * ptr_int%geofac_n2s(jc,4,jb)

    END DO
  END DO
!$OMP END DO
!$OMP END PARALLEL

ELSE

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk)
  DO jb = i_startblk, i_endblk

    CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, rl_start, rl_end)

#ifdef __LOOP_EXCHANGE
    DO jc = i_startidx, i_endidx
      DO jk = slev, elev
#else
#ifdef _URD
!CDIR UNROLL=_URD
#endif
    DO jk = slev, elev
      DO jc = i_startidx, i_endidx
#endif
        !
        !  calculate div(grad) in one step
        !
        nabla2_psi_c(jc,jk,jb) =  &
          &    psi_c(jc,jk,jb)                       * ptr_int%geofac_n2s(jc,1,jb) &
          &  + psi_c(iidx(jc,jb,1),jk,iblk(jc,jb,1)) * ptr_int%geofac_n2s(jc,2,jb) &
          &  + psi_c(iidx(jc,jb,2),jk,iblk(jc,jb,2)) * ptr_int%geofac_n2s(jc,3,jb) &
          &  + psi_c(iidx(jc,jb,3),jk,iblk(jc,jb,3)) * ptr_int%geofac_n2s(jc,4,jb)
      END DO !cell loop

    END DO !vertical levels loop
  END DO
!$OMP END DO
!$OMP END PARALLEL
ENDIF
CASE (6) ! (cell_type == 6) THEN ! Use unoptimized version for the time being

  ! compute finite difference gradient in normal direction
  CALL grad_fd_norm( psi_c, ptr_patch, z_grad_fd_norm_e, slev, elev)

  ! compute divergence of resulting vector field
  CALL div( z_grad_fd_norm_e, ptr_patch, ptr_int, nabla2_psi_c, slev, elev)

END SELECT

END SUBROUTINE nabla2_scalar

!-------------------------------------------------------------------------
!

!>
!!  Computes Laplacian @f$\nabla ^2 @f$ of a scalar field, followed by weighted averaging.
!!
!!  Computes Laplacian @f$\nabla ^2 @f$ of a scalar field, followed by weighted averaging
!!  with the neighboring cells to increase computing efficiency.
!!  NOTE: This optimized routine works for triangular grids only.
!! input:  lives on cells (mass points)
!! output: lives on cells
!!
!! @par Revision History
!! Developed by Guenther Zaengl, DWD, 2009-05-19
!!
SUBROUTINE nabla2_scalar_avg( psi_c, ptr_patch, ptr_int, avg_coeff, nabla2_psi_c, &
  &                           opt_slev, opt_elev )
!
!  patch on which computation is performed
!
TYPE(t_patch), TARGET, INTENT(in) :: ptr_patch

! Interpolation state
TYPE(t_int_state), INTENT(in)     :: ptr_int

!  averaging coefficients
REAL(wp), INTENT(in) :: avg_coeff(:,:,:) ! dim: (nproma,nlev,nblks_c)

!
!  cells based variable of which biharmonic laplacian is computed
!
REAL(wp), INTENT(in) ::  &
  &  psi_c(:,:,:) ! dim: (nproma,nlev,nblks_c)

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_slev    ! optional vertical start level

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_elev    ! optional vertical end level


!
!  cell based variable in which biharmonic laplacian is stored
!
!REAL(wp), INTENT(out) ::  &
REAL(wp), INTENT(inout) ::  &
  &  nabla2_psi_c(:,:,:) ! dim: (nproma,nlev,nblks_c)

INTEGER :: slev, elev     ! vertical start and end level
INTEGER :: rl_start, rl_end, rl_start_l2
INTEGER :: jb, jc, jk, i_startblk, i_endblk, i_startidx, i_endidx, i_nchdom
INTEGER :: nlev              !< number of full levels

REAL(wp), DIMENSION (nproma,ptr_patch%nlev,ptr_patch%nblks_c) :: aux_c

INTEGER,  DIMENSION(:,:,:),   POINTER :: iidx, iblk

!-----------------------------------------------------------------------

! number of vertical levels
nlev = ptr_patch%nlev

! check optional arguments
IF ( PRESENT(opt_slev) ) THEN
  slev = opt_slev
ELSE
  slev = 1
END IF
IF ( PRESENT(opt_elev) ) THEN
  elev = opt_elev
ELSE
  elev = nlev
END IF

iidx => ptr_patch%cells%neighbor_idx
iblk => ptr_patch%cells%neighbor_blk

rl_start = 2
rl_start_l2 = rl_start + 1
rl_end = min_rlcell
i_nchdom   = MAX(1,ptr_patch%n_childdom)

! The special treatment of 2D fields is essential for efficiency on the NEC

SELECT CASE (ptr_patch%cell_type)

CASE (3) ! (cell_type == 3)

IF (slev == elev) THEN
  jk = slev
!$OMP PARALLEL PRIVATE(i_startblk,i_endblk)

! values for the blocking
i_startblk = ptr_patch%cells%start_blk(rl_start,1)
i_endblk   = ptr_patch%cells%end_blk(rl_end,i_nchdom)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc)
  DO jb = i_startblk, i_endblk

    IF (jb == i_startblk) THEN
      i_startidx = ptr_patch%cells%start_idx(rl_start,1)
      i_endidx   = nproma
      IF (jb == i_endblk) i_endidx = ptr_patch%cells%end_idx(rl_end,i_nchdom)
    ELSE IF (jb == i_endblk) THEN
      i_startidx = 1
      i_endidx   = ptr_patch%cells%end_idx(rl_end,i_nchdom)
    ELSE
      i_startidx = 1
      i_endidx   = nproma
    ENDIF

    DO jc = i_startidx, i_endidx

      !
      !  calculate div(grad) in one step
      !
      aux_c(jc,jk,jb) =  &
        &    psi_c(jc,jk,jb)                       * ptr_int%geofac_n2s(jc,1,jb) &
        &  + psi_c(iidx(jc,jb,1),jk,iblk(jc,jb,1)) * ptr_int%geofac_n2s(jc,2,jb) &
        &  + psi_c(iidx(jc,jb,2),jk,iblk(jc,jb,2)) * ptr_int%geofac_n2s(jc,3,jb) &
        &  + psi_c(iidx(jc,jb,3),jk,iblk(jc,jb,3)) * ptr_int%geofac_n2s(jc,4,jb)

    END DO
  END DO
!$OMP END DO

  IF (l_limited_area .OR. ptr_patch%id > 1) THEN
    ! Fill nabla2_psi_c along the lateral boundaries of nests

    i_startblk = ptr_patch%cells%start_blk(rl_start,1)
    i_endblk   = ptr_patch%cells%end_blk(rl_start_l2,1)

!$OMP WORKSHARE
    nabla2_psi_c(:,jk,i_startblk:i_endblk) =  aux_c(:,jk,i_startblk:i_endblk)
!$OMP END WORKSHARE
  ENDIF

!
! Now do averaging with weights given by avg_coeff

  ! values for the blocking
  i_startblk = ptr_patch%cells%start_blk(rl_start_l2,1)
  i_endblk   = ptr_patch%cells%end_blk(rl_end,i_nchdom)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc)
  DO jb = i_startblk, i_endblk

    IF (jb == i_startblk) THEN
      i_startidx = ptr_patch%cells%start_idx(rl_start_l2,1)
      i_endidx   = nproma
      IF (jb == i_endblk) i_endidx = ptr_patch%cells%end_idx(rl_end,i_nchdom)
    ELSE IF (jb == i_endblk) THEN
      i_startidx = 1
      i_endidx   = ptr_patch%cells%end_idx(rl_end,i_nchdom)
    ELSE
      i_startidx = 1
      i_endidx   = nproma
    ENDIF

    DO jc = i_startidx, i_endidx
      !
      !  calculate the weighted average
      !
      nabla2_psi_c(jc,jk,jb) =  &
        &    aux_c(jc,jk,jb)                       * avg_coeff(jc,1,jb) &
        &  + aux_c(iidx(jc,jb,1),jk,iblk(jc,jb,1)) * avg_coeff(jc,2,jb) &
        &  + aux_c(iidx(jc,jb,2),jk,iblk(jc,jb,2)) * avg_coeff(jc,3,jb) &
        &  + aux_c(iidx(jc,jb,3),jk,iblk(jc,jb,3)) * avg_coeff(jc,4,jb)

    END DO !cell loop

  END DO !block loop
!$OMP END DO NOWAIT
!$OMP END PARALLEL

ELSE

!$OMP PARALLEL PRIVATE(i_startblk,i_endblk)

! values for the blocking
i_startblk = ptr_patch%cells%start_blk(rl_start,1)
i_endblk   = ptr_patch%cells%end_blk(rl_end,i_nchdom)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk)
  DO jb = i_startblk, i_endblk

    CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, rl_start, rl_end)

#ifdef __LOOP_EXCHANGE
    DO jc = i_startidx, i_endidx
      DO jk = slev, elev
#else
#ifdef _URD
!CDIR UNROLL=_URD
#endif
    DO jk = slev, elev
      DO jc = i_startidx, i_endidx
#endif
        !
        !  calculate div(grad) in one step
        !
        aux_c(jc,jk,jb) =  &
          &    psi_c(jc,jk,jb)                       * ptr_int%geofac_n2s(jc,1,jb) &
          &  + psi_c(iidx(jc,jb,1),jk,iblk(jc,jb,1)) * ptr_int%geofac_n2s(jc,2,jb) &
          &  + psi_c(iidx(jc,jb,2),jk,iblk(jc,jb,2)) * ptr_int%geofac_n2s(jc,3,jb) &
          &  + psi_c(iidx(jc,jb,3),jk,iblk(jc,jb,3)) * ptr_int%geofac_n2s(jc,4,jb)
      END DO !cell loop

    END DO !vertical levels loop
  END DO
!$OMP END DO

  IF (l_limited_area .OR. ptr_patch%id > 1) THEN
    ! Fill nabla2_psi_c along the lateral boundaries of nests

    i_startblk = ptr_patch%cells%start_blk(rl_start,1)
    i_endblk   = ptr_patch%cells%end_blk(rl_start_l2,1)

!$OMP WORKSHARE
    nabla2_psi_c(:,:,i_startblk:i_endblk) =  aux_c(:,:,i_startblk:i_endblk)
!$OMP END WORKSHARE
  ENDIF

!
! Now do averaging with weights given by avg_coeff

  ! values for the blocking
  i_startblk = ptr_patch%cells%start_blk(rl_start_l2,1)
  i_endblk   = ptr_patch%cells%end_blk(rl_end,i_nchdom)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk)
  DO jb = i_startblk, i_endblk

    CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, rl_start_l2, rl_end)

#ifdef __LOOP_EXCHANGE
    DO jc = i_startidx, i_endidx
      DO jk = slev, elev
#else
#ifdef _URD
!CDIR UNROLL=_URD
#endif
    DO jk = slev, elev
      DO jc = i_startidx, i_endidx
#endif
        !
        !  calculate the weighted average
        !
        nabla2_psi_c(jc,jk,jb) =  &
          &    aux_c(jc,jk,jb)                       * avg_coeff(jc,1,jb) &
          &  + aux_c(iidx(jc,jb,1),jk,iblk(jc,jb,1)) * avg_coeff(jc,2,jb) &
          &  + aux_c(iidx(jc,jb,2),jk,iblk(jc,jb,2)) * avg_coeff(jc,3,jb) &
          &  + aux_c(iidx(jc,jb,3),jk,iblk(jc,jb,3)) * avg_coeff(jc,4,jb)

      END DO !cell loop

    END DO !vertical levels loop

  END DO !block loop
!$OMP END DO NOWAIT

!$OMP END PARALLEL
ENDIF
END SELECT

END SUBROUTINE nabla2_scalar_avg

!-------------------------------------------------------------------------
!
!

!>
!!  Computes biharmonic laplacian @f$\nabla ^4 @f$ of a scalar field.
!!
!! input:  lives on edges (velocity points)
!! output: lives on edges
!!
!! @par Revision History
!! Developed and tested  by L.Bonaventura  (2002-4).
!! Adapted to new grid and patch structure by P. Korn (2005).
!! Modifications by P. Korn, MPI-M(2007-2)
!! -Switch fom array arguments to pointers
!! Modifications by Almut Gassmann, MPI-M (2007-04-20)
!! - abandon grid for the sake of patch
!! - corrected type for temp1 (now edges and no longer cells)
!!
SUBROUTINE nabla4_scalar( psi_c, ptr_patch, ptr_int, nabla4_psi_c, &
  &                       opt_nabla2, opt_slev, opt_elev, opt_rlstart, opt_rlend )

!
!  patch on which computation is performed
!
TYPE(t_patch), TARGET, INTENT(in) :: ptr_patch

! Interpolation state
TYPE(t_int_state), INTENT(in)     :: ptr_int
!
!  cells based variable of which biharmonic laplacian is computed
!
REAL(wp), INTENT(in) ::  &
  &  psi_c(:,:,:) ! dim: (nproma,nlev,nblks_c)

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_slev    ! optional vertical start level

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_elev    ! optional vertical end level

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_rlstart, opt_rlend   ! start and end values of refin_ctrl flag

!
!  cell based variable in which biharmonic laplacian is stored
!
!REAL(wp), INTENT(out) ::  &
REAL(wp), INTENT(inout) ::  &
  &  nabla4_psi_c(:,:,:) ! dim: (nproma,nlev,nblks_c)

! Optional argument for passing nabla2 to the calling program
! (to avoid double computation for Smagorinsky diffusion and nest boundary diffusion)
REAL(wp), INTENT(inout), TARGET, OPTIONAL  ::  &
  &  opt_nabla2(:,:,:) ! dim: (nproma,nlev,nblks_e)

INTEGER :: slev, elev     ! vertical start and end level
INTEGER :: rl_start, rl_end
INTEGER :: rl_start_s1, rl_end_s1


REAL(wp), ALLOCATABLE, TARGET :: z_nab2_c(:,:,:) ! dim: (nproma,nlev,ptr_patch%nblks_c)
REAL(wp), POINTER :: p_nabla2(:,:,:)

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
  IF ((opt_rlstart >= 0) .AND. (opt_rlstart <= 2)) THEN
    CALL finish ('mo_math_operators:nabla4_scalar',  &
          &      'opt_rlstart must not be between 0 and 2')
  ENDIF
  rl_start = opt_rlstart
ELSE
  rl_start = 3
END IF
IF ( PRESENT(opt_rlend) ) THEN
  rl_end = opt_rlend
ELSE
  rl_end = min_rlcell
END IF

IF (rl_start > 0) THEN
  rl_start_s1 = rl_start - 1
ELSE
  rl_start_s1 = rl_start + 1
ENDIF
IF (rl_end > 0) THEN
  rl_end_s1 = rl_end + 1
ELSE
  rl_end_s1 = rl_end - 1
ENDIF

rl_end_s1 = MAX(min_rlcell,rl_end_s1)

IF (PRESENT(opt_nabla2) ) THEN
  p_nabla2 => opt_nabla2
ELSE
  ALLOCATE (z_nab2_c(nproma,ptr_patch%nlev,ptr_patch%nblks_c))
  p_nabla2 => z_nab2_c
ENDIF

! apply second order Laplacian twice

CALL nabla2_scalar( psi_c, ptr_patch, ptr_int, p_nabla2, &
                    slev, elev, opt_rlstart=rl_start_s1, opt_rlend=rl_end_s1 )

IF (cell_type == 6) CALL sync_patch_array(SYNC_C, ptr_patch, p_nabla2)

CALL nabla2_scalar( p_nabla2, ptr_patch, ptr_int, nabla4_psi_c, &
                    slev, elev, opt_rlstart=rl_start, opt_rlend=rl_end )

END SUBROUTINE nabla4_scalar


  !----------------------------------------------------------------------
  !>
  !! directional laplace of a scalar in upstream direction at an edge
  !!
  !! For third order advection of theta_v or other scalars one needs the laplacian in
  !! the edge direction. Here, the laplacian at cell centers is reconstructed
  !! from a least square fit. According to the given velocity field, the evaluation is
  !! done at the upstream located cell.
  !!
  !! CURRENTLY ONLY FOR HEXAGONAL MODEL
  !!
  !! @par Revision History
  !! Initial release by Almut Gassmann, MPI_M (2010-09-06)
  !!
  SUBROUTINE directional_laplace (p_vn, p_scalar, pt_patch, pt_int, p_upstr_beta, &
    &                             p_lapl_e, opt_slev, opt_elev)

    TYPE(t_patch), TARGET, INTENT(in) :: pt_patch           !< patch
    TYPE(t_int_state), TARGET,INTENT(in) :: pt_int !< interpolation state
    REAL(wp), INTENT(in)    :: p_vn(:,:,:)  &
    & ; !< normal velocity field, needed for determination of upstream direction
    REAL(wp), INTENT(in)    :: p_scalar(:,:,:) !< scalar
    REAL(wp), INTENT(out)   :: p_lapl_e(:,:,:) !< upstream Laplace projected onto edge normal
    REAL(wp), INTENT(in)    :: p_upstr_beta

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional vertical start level
      &  opt_slev

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional vertical end level
      &  opt_elev

    REAL(wp) :: coeff(nproma,pt_patch%nlev,pt_patch%nblks_c,6) &
    & ; !< reconstruction coefficients
    INTEGER, POINTER :: iidx(:,:,:), iblk(:,:,:) !< stencil

    INTEGER  :: slev, elev         !< vertical start and end level
    INTEGER :: nblks_e, npromz_e, nlen, jb, jk, je

#ifdef __LOOP_EXCHANGE
    REAL(wp) :: z_sign(pt_patch%nlev)
#else
    REAL(wp) :: z_sign(nproma)
#endif
    !-----------------------------------------------------------------

    ! check optional arguments
    IF ( PRESENT(opt_slev) ) THEN
      slev = opt_slev
    ELSE
      slev = 1
    END IF
    IF ( PRESENT(opt_elev) ) THEN
      elev = opt_elev
    ELSE
      elev = pt_patch%nlev
    END IF

    CALL recon_lsq_cell_q(p_scalar, pt_patch, pt_int%lsq_high, coeff, &
      &                   opt_slev=slev, opt_elev=elev)
    CALL sync_patch_array(SYNC_C,pt_patch,coeff(:,:,:,4))
    CALL sync_patch_array(SYNC_C,pt_patch,coeff(:,:,:,5))
    CALL sync_patch_array(SYNC_C,pt_patch,coeff(:,:,:,6))
    ! Note: The coefficiens have to be multiplied with 2, which is implicitly
    ! done below.

    iidx => pt_patch%edges%cell_idx
    iblk => pt_patch%edges%cell_blk

    nblks_e  = pt_patch%nblks_int_e
    npromz_e = pt_patch%npromz_int_e

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
        DO jk = slev, elev
          z_sign(jk) = p_upstr_beta*SIGN(1.0_wp,&
          &            p_vn(je,jk,jb)*pt_patch%edges%system_orientation(je,jb))
          p_lapl_e(je,jk,jb) = (1.0_wp+z_sign(jk))  &
          &     *(pt_int%cea_en(je,1,jb)**2*coeff(iidx(je,jb,1),jk,iblk(je,jb,1),4) &
          &      +pt_int%cno_en(je,1,jb)**2*coeff(iidx(je,jb,1),jk,iblk(je,jb,1),5) &
          &      +pt_int%cno_en(je,1,jb)*pt_int%cea_en(je,1,jb)&
          &                                *coeff(iidx(je,jb,1),jk,iblk(je,jb,1),6))&
          &                  + (1.0_wp-z_sign(jk))  &
          &     *(pt_int%cea_en(je,2,jb)**2*coeff(iidx(je,jb,2),jk,iblk(je,jb,2),4) &
          &      +pt_int%cno_en(je,2,jb)**2*coeff(iidx(je,jb,2),jk,iblk(je,jb,2),5) &
          &      +pt_int%cno_en(je,2,jb)*pt_int%cea_en(je,2,jb)&
          &                                *coeff(iidx(je,jb,2),jk,iblk(je,jb,2),6))

        ENDDO
      ENDDO
#else
      DO jk = slev, elev
        DO je = 1, nlen
          z_sign(je) = p_upstr_beta*SIGN(1.0_wp,&
          &            p_vn(je,jk,jb)*pt_patch%edges%system_orientation(je,jb))
          p_lapl_e(je,jk,jb) = (1.0_wp+z_sign(je))  &
          &     *(pt_int%cea_en(je,1,jb)**2*coeff(iidx(je,jb,1),jk,iblk(je,jb,1),4) &
          &      +pt_int%cno_en(je,1,jb)**2*coeff(iidx(je,jb,1),jk,iblk(je,jb,1),5) &
          &      +pt_int%cno_en(je,1,jb)*pt_int%cea_en(je,1,jb)&
          &                                *coeff(iidx(je,jb,1),jk,iblk(je,jb,1),6))&
          &                  + (1.0_wp-z_sign(je))  &
          &     *(pt_int%cea_en(je,2,jb)**2*coeff(iidx(je,jb,2),jk,iblk(je,jb,2),4) &
          &      +pt_int%cno_en(je,2,jb)**2*coeff(iidx(je,jb,2),jk,iblk(je,jb,2),5) &
          &      +pt_int%cno_en(je,2,jb)*pt_int%cea_en(je,2,jb)&
          &                                *coeff(iidx(je,jb,2),jk,iblk(je,jb,2),6))

        ENDDO
      ENDDO
#endif
    ENDDO
!$OMP END DO
!$OMP END PARALLEL

  END SUBROUTINE directional_laplace

END MODULE mo_math_laplace
