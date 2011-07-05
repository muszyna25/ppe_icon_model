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
MODULE mo_intp
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
USE mo_exception,           ONLY: finish
USE mo_impl_constants,      ONLY: min_rlcell, min_rledge, min_rlvert
USE mo_model_domain,        ONLY: t_patch
USE mo_parallel_configuration,  ONLY: nproma
USE mo_run_nml,             ONLY: i_cell_type, ltimer
USE mo_loopindices,         ONLY: get_indices_c, get_indices_e, get_indices_v

USE mo_intp_data_strc
USE mo_timer,               ONLY: timer_start, timer_stop, timer_intp

IMPLICIT NONE

PRIVATE

PUBLIC :: verts2edges_scalar, cells2edges_scalar, edges2verts_scalar,   &
          & edges2cells_scalar, edges2cells_vector, cells2verts_scalar, &
          & verts2cells_scalar, cell_avg, edges2edges_scalar
CONTAINS

!-----------------------------------------------------------------------
!
!  ! averaging and interpolation routines and
!  ! routines needed to compute the coefficients therein
!
!-----------------------------------------------------------------------
!
!>
!!  Performs  average of scalar fields from vertices to velocity points.
!!
!!  The coefficients are given by c_int.
!!
!! @par Revision History
!! Developed  by L.Bonaventura  (2004).
!! Adapted to new grid and patch structure by P. Korn (2005)
!! Modifications by Almut Gassmann (2007-04-30)
!! -abandon grid for the sake of patch
!! Modifications by Almut Gassmann (2009-01-28)
!! - Rename (naming structure adapted to other interpolation routines) and
!!   rewrite this routine for accepting given interpolation weights.
!!
SUBROUTINE verts2edges_scalar( p_vertex_in, ptr_patch, c_int, p_edge_out, &
  &                            opt_slev, opt_elev, opt_rlstart, opt_rlend )
!
!

TYPE(t_patch), TARGET, INTENT(in) :: ptr_patch

! vertex based scalar input field
REAL(wp), INTENT(in) ::  p_vertex_in(:,:,:)  ! dim: (nproma,nlev,nblks_v)
! interpolation field
REAL(wp), INTENT(in) ::  c_int(:,:,:)        ! dim: (nproma,2,nblks_e)

INTEGER, INTENT(in), OPTIONAL ::  opt_slev   ! optional vertical start level

INTEGER, INTENT(in), OPTIONAL ::  opt_elev   ! optional vertical end level

! start and end values of refin_ctrl flag
INTEGER, INTENT(in), OPTIONAL :: opt_rlstart, opt_rlend

! edge based scalar output field
REAL(wp), INTENT(inout) :: p_edge_out(:,:,:) ! dim: (nproma,nlev,nblks_e)

INTEGER :: slev, elev     ! vertical start and end level
INTEGER :: i_startblk     ! start block
INTEGER :: i_endblk       ! end block
INTEGER :: i_startidx     ! start index
INTEGER :: i_endidx       ! end index
INTEGER :: rl_start, rl_end, i_nchdom
INTEGER :: je, jk, jb

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
  rl_end = min_rledge
END IF

iidx => ptr_patch%edges%vertex_idx
iblk => ptr_patch%edges%vertex_blk

i_nchdom   = MAX(1,ptr_patch%n_childdom)
i_startblk = ptr_patch%edges%start_blk(rl_start,1)
i_endblk   = ptr_patch%edges%end_blk(rl_end,i_nchdom)

IF (ltimer) CALL timer_start(timer_intp)

! loop over edges and blocks
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,je,jk)
DO jb = i_startblk, i_endblk

  CALL get_indices_e(ptr_patch, jb, i_startblk, i_endblk, &
                     i_startidx, i_endidx, rl_start, rl_end)

#ifdef __LOOP_EXCHANGE
  DO je = i_startidx, i_endidx
    DO jk = slev, elev
#else
#ifdef _URD
!CDIR UNROLL=_URD
#endif
  DO jk = slev, elev
    DO je = i_startidx, i_endidx
#endif

      p_edge_out(je,jk,jb) =  &
        &   c_int(je,1,jb)*p_vertex_in(iidx(je,jb,1),jk,iblk(je,jb,1))  &
        & + c_int(je,2,jb)*p_vertex_in(iidx(je,jb,2),jk,iblk(je,jb,2))

    END DO

  END DO

END DO
!$OMP END DO
!$OMP END PARALLEL

IF (ltimer) CALL timer_stop(timer_intp)


END SUBROUTINE verts2edges_scalar

!------------------------------------------------------------------------
!
!
!>
!!  Computes  average of scalar fields from centers of triangular faces to.
!!
!!  Computes  average of scalar fields from centers of triangular faces to
!!  velocity points.
!!
!! @par Revision History
!! Developed  by L.Bonaventura  (2004).
!! Adapted to new grid and patch structure by P. Korn and  L.Bonaventura  (2005)
!! Modification by Thomas Heinze, DWD (2006-08-18):
!! - changed according to programming guide
!! - replaced index iic1 and ic2 by j at the end of code (after discussion with
!!   Hui Wan and Peter Korn)
!! Modification by Almut Gassmann, MPI-M (2007-04-30)
!! - abandon grid for the sake of patch
!! Modification by Thomas Heinze, DWD (2007-08-06):
!! - changed c_int to POINTER (was array)
!!
SUBROUTINE cells2edges_scalar( p_cell_in, ptr_patch, c_int, p_edge_out,  &
  &                            opt_slev, opt_elev, opt_rlstart, opt_rlend )
!

TYPE(t_patch), TARGET, INTENT(in) :: ptr_patch

! cell based scalar input field
REAL(wp), INTENT(in) :: p_cell_in(:,:,:)   ! dim: (nproma,nlev,nblks_c)

! coefficients for linear interpolation
REAL(wp), INTENT(in) :: c_int(:,:,:)       ! dim: (nproma,2,nblks_e)

INTEGER, INTENT(in), OPTIONAL :: opt_slev  ! optional vertical start level

INTEGER, INTENT(in), OPTIONAL :: opt_elev  ! optional vertical end level

! start and end values of refin_ctrl flag
INTEGER, INTENT(in), OPTIONAL :: opt_rlstart, opt_rlend

! edge based scalar output field
REAL(wp), INTENT(inout) :: p_edge_out(:,:,:) ! dim: (nproma,nlev,nblks_e)

INTEGER :: slev, elev     ! vertical start and end level
INTEGER :: je, jk, jb
INTEGER :: i_startblk                ! start block
INTEGER :: i_endblk                  ! end block
INTEGER :: i_startidx                ! start index
INTEGER :: i_endidx                  ! end index
INTEGER :: rl_start, rl_end, i_nchdom

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
  ! The calculation cannot be done for boundary edges
  IF (opt_rlstart == 1) THEN
    CALL finish ('mo_interpolation:cells2edges_scalar',  &
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


IF (ltimer) CALL timer_start(timer_intp)

!$OMP PARALLEL
IF (slev > 1) THEN
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,je,jk)
  DO jb = i_startblk, i_endblk

    CALL get_indices_e(ptr_patch, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, rl_start, rl_end)

#ifdef __LOOP_EXCHANGE
    DO je = i_startidx, i_endidx
      DO jk = slev, elev
#else
#ifdef _URD2
!CDIR UNROLL=_URD2
#endif
    DO jk = slev, elev
      DO je = i_startidx, i_endidx
#endif

        p_edge_out(je,jk,jb) =  &
          &    c_int(je,1,jb) * p_cell_in(iidx(je,jb,1),jk,iblk(je,jb,1))  &
          &  + c_int(je,2,jb) * p_cell_in(iidx(je,jb,2),jk,iblk(je,jb,2))

      END DO
    END DO

  END DO
!$OMP END DO
ELSE
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,je,jk)
  DO jb = i_startblk, i_endblk

    CALL get_indices_e(ptr_patch, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, rl_start, rl_end)

#ifdef __LOOP_EXCHANGE
    DO je = i_startidx, i_endidx
      DO jk = slev, elev
#else
#ifdef _URD
!CDIR UNROLL=_URD
#endif
    DO jk = slev, elev
      DO je = i_startidx, i_endidx
#endif

        p_edge_out(je,jk,jb) =  &
          &    c_int(je,1,jb) * p_cell_in(iidx(je,jb,1),jk,iblk(je,jb,1))  &
          &  + c_int(je,2,jb) * p_cell_in(iidx(je,jb,2),jk,iblk(je,jb,2))

      END DO
    END DO

  END DO
!$OMP END DO
ENDIF
!$OMP END PARALLEL

IF (ltimer) CALL timer_stop(timer_intp)


END SUBROUTINE cells2edges_scalar


!------------------------------------------------------------------------
!
!>
!!  Computes average of scalar fields from velocity points to.
!!
!!  Computes average of scalar fields from velocity points to
!!  centers of dual faces.
!!
!! @par Revision History
!!  Initial version by Almut Gassmann, MPI-M (2009-04-23)
!!  - exchange calling parameters to have the same shape as for other
!!    interpolation routines
!!
SUBROUTINE edges2verts_scalar( p_edge_in, ptr_patch, v_int, p_vert_out,  &
  &                            opt_slev, opt_elev, opt_rlstart, opt_rlend )
!

TYPE(t_patch), TARGET, INTENT(in) :: ptr_patch

! edge based scalar input field
REAL(wp), INTENT(in) ::  p_edge_in(:,:,:)  ! dim: (nproma,nlev,nblks_e)

! coefficients for (area weighted) interpolation
REAL(wp), INTENT(in) ::  v_int(:,:,:)      ! dim: (nproma,i_cell_type,nblks_v)

INTEGER, INTENT(in), OPTIONAL ::  opt_slev ! optional vertical start level

INTEGER, INTENT(in), OPTIONAL ::  opt_elev ! optional vertical end level

! start and end values of refin_ctrl flag
INTEGER, INTENT(in), OPTIONAL ::  opt_rlstart, opt_rlend

! vertex based scalar output field
REAL(wp), INTENT(inout) :: p_vert_out(:,:,:)  ! dim: (nproma,nlev,nblks_v)

INTEGER :: slev, elev     ! vertical start and end level
INTEGER :: jv, jk, jb
INTEGER :: nblks_v, npromz_v, nlen
INTEGER :: rl_start, rl_end
INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx, i_nchdom

INTEGER,  DIMENSION(:,:,:),   POINTER :: iidx, iblk

!-------------------------------------------------------------------------

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
  rl_end = min_rlvert
END IF

iidx => ptr_patch%verts%edge_idx
iblk => ptr_patch%verts%edge_blk

! values for the blocking
i_nchdom   = MAX(1,ptr_patch%n_childdom)
i_startblk = ptr_patch%verts%start_blk(rl_start,1)
i_endblk   = ptr_patch%verts%end_blk(rl_end,i_nchdom)


IF (ltimer) CALL timer_start(timer_intp)

!loop over blocks and verts
IF (i_cell_type == 6) THEN

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

#ifdef __LOOP_EXCHANGE
    DO jv = 1, nlen
      DO jk = slev, elev
#else
!CDIR UNROLL=6
    DO jk = slev, elev
      DO jv = 1, nlen
#endif

        p_vert_out(jv,jk,jb) =  &
          v_int(jv,1,jb) * p_edge_in(iidx(jv,jb,1),jk,iblk(jv,jb,1)) + &
          v_int(jv,2,jb) * p_edge_in(iidx(jv,jb,2),jk,iblk(jv,jb,2)) + &
          v_int(jv,3,jb) * p_edge_in(iidx(jv,jb,3),jk,iblk(jv,jb,3))

      ENDDO
    ENDDO
  ENDDO  !loop over blocks
!$OMP END DO
!$OMP END PARALLEL

ELSE IF (i_cell_type == 3) THEN

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jv,jk)
  DO jb = i_startblk, i_endblk

    CALL get_indices_v(ptr_patch, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, rl_start, rl_end)

#ifdef __LOOP_EXCHANGE
    DO jv = i_startidx, i_endidx
      DO jk = slev, elev
#else
!CDIR UNROLL=6
    DO jk = slev, elev
      DO jv = i_startidx, i_endidx
#endif

        p_vert_out(jv,jk,jb) =  &
          v_int(jv,1,jb) * p_edge_in(iidx(jv,jb,1),jk,iblk(jv,jb,1)) + &
          v_int(jv,2,jb) * p_edge_in(iidx(jv,jb,2),jk,iblk(jv,jb,2)) + &
          v_int(jv,3,jb) * p_edge_in(iidx(jv,jb,3),jk,iblk(jv,jb,3)) + &
          v_int(jv,4,jb) * p_edge_in(iidx(jv,jb,4),jk,iblk(jv,jb,4)) + &
          v_int(jv,5,jb) * p_edge_in(iidx(jv,jb,5),jk,iblk(jv,jb,5)) + &
          v_int(jv,6,jb) * p_edge_in(iidx(jv,jb,6),jk,iblk(jv,jb,6))

      ENDDO
    ENDDO

  ENDDO  !loop over blocks
!$OMP END DO
!$OMP END PARALLEL
ENDIF

IF (ltimer) CALL timer_stop(timer_intp)


END SUBROUTINE edges2verts_scalar

!------------------------------------------------------------------------
!
!>
!!  Computes interpolation from edges to cells
!!
!!  Computes interpolation of scalar fields from velocity points to
!!  cell centers via given interpolation weights
!!
!! @par Revision History
!!  Original version by Hui Wan (MPI-M, 2006-08-17)
!!  Modification by Almut Gassmann, MPI-M (2009-01-28)
!!  - exchange calling parameters to have the same shape as for other
!!    interpolation routines
!!
SUBROUTINE edges2cells_scalar( p_edge_in, ptr_patch, c_int, p_cell_out,  &
  &                            opt_slev, opt_elev, opt_rlstart, opt_rlend )
!

TYPE(t_patch), TARGET, INTENT(in) :: ptr_patch

! edge based scalar input field
REAL(wp), INTENT(in) ::  p_edge_in(:,:,:)  ! dim: (nproma,nlev,nblks_e)

! coefficients for (area weighted) interpolation
REAL(wp), INTENT(in) ::  c_int(:,:,:)      ! dim: (nproma,i_cell_type,nblks_c)

INTEGER, INTENT(in), OPTIONAL ::  opt_slev ! optional vertical start level

INTEGER, INTENT(in), OPTIONAL ::  opt_elev ! optional vertical end level

! start and end values of refin_ctrl flag
INTEGER, INTENT(in), OPTIONAL ::  opt_rlstart, opt_rlend

! cell based scalar output field
REAL(wp), INTENT(inout) :: p_cell_out(:,:,:)  ! dim: (nproma,nlev,nblks_c)

INTEGER :: slev, elev     ! vertical start and end level
INTEGER :: jc, jk, jb
INTEGER :: rl_start, rl_end
INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx, i_nchdom

INTEGER,  DIMENSION(:,:,:),   POINTER :: iidx, iblk

!-------------------------------------------------------------------------

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


IF (ltimer) CALL timer_start(timer_intp)

!loop over blocks and cells
IF (i_cell_type == 3) THEN
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk)
  DO jb = i_startblk, i_endblk

    CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, rl_start, rl_end)

#ifdef __LOOP_EXCHANGE
    DO jc = i_startidx, i_endidx
      DO jk = slev, elev
#else
!CDIR UNROLL=6
    DO jk = slev, elev
      DO jc = i_startidx, i_endidx
#endif

        p_cell_out(jc,jk,jb) =  &
          c_int(jc,1,jb) * p_edge_in(iidx(jc,jb,1),jk,iblk(jc,jb,1)) + &
          c_int(jc,2,jb) * p_edge_in(iidx(jc,jb,2),jk,iblk(jc,jb,2)) + &
          c_int(jc,3,jb) * p_edge_in(iidx(jc,jb,3),jk,iblk(jc,jb,3))

      ENDDO
    ENDDO

  ENDDO  !loop over blocks
!$OMP END DO
!$OMP END PARALLEL

ELSE IF (i_cell_type == 6) THEN

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk)
  DO jb = i_startblk, i_endblk

    CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, rl_start, rl_end)

#ifdef __LOOP_EXCHANGE
    DO jc = i_startidx, i_endidx
      DO jk = slev, elev
#else
!CDIR UNROLL=6
    DO jk = slev, elev
      DO jc = i_startidx, i_endidx
#endif

        p_cell_out(jc,jk,jb) =  &
          c_int(jc,1,jb) * p_edge_in(iidx(jc,jb,1),jk,iblk(jc,jb,1)) + &
          c_int(jc,2,jb) * p_edge_in(iidx(jc,jb,2),jk,iblk(jc,jb,2)) + &
          c_int(jc,3,jb) * p_edge_in(iidx(jc,jb,3),jk,iblk(jc,jb,3)) + &
          c_int(jc,4,jb) * p_edge_in(iidx(jc,jb,4),jk,iblk(jc,jb,4)) + &
          c_int(jc,5,jb) * p_edge_in(iidx(jc,jb,5),jk,iblk(jc,jb,5)) + &
          c_int(jc,6,jb) * p_edge_in(iidx(jc,jb,6),jk,iblk(jc,jb,6))

      ENDDO
    ENDDO

  ENDDO  !loop over blocks
!$OMP END DO
!$OMP END PARALLEL
ENDIF

IF (ltimer) CALL timer_stop(timer_intp)


END SUBROUTINE edges2cells_scalar

!------------------------------------------------------------------------
!
!>
!!  Bilinear interpolation of normal and tangential velocity components.
!!
!!  Bilinear interpolation of normal and tangential velocity components
!!  at the edges to u and v at the cells
!!  Works only for triangles (bilinear interpolation weights are not implemented
!!  for hexagons)
!!
!! @par Revision History
!!  Developed by Guenther Zaengl, DWD (2009-05-08)
!!
SUBROUTINE edges2cells_vector( p_vn_in, p_vt_in, ptr_patch, p_int, &
  &                            p_uv_out, opt_slev, opt_elev )
!


TYPE(t_patch), TARGET, INTENT(in) :: ptr_patch

! normal velocity component at edges
REAL(wp), INTENT(in) ::  p_vn_in(:,:,:)  ! dim: (nproma,nlev,nblks_e)

! (reconstructed) tangential velocity component at edges
REAL(wp), INTENT(in) ::  p_vt_in(:,:,:)  ! dim: (nproma,nlev,nblks_e)

! Interpolation state
TYPE(t_int_state), INTENT(IN) :: p_int

! cell based output field: u and v
REAL(wp), INTENT(inout) :: p_uv_out(:,:,:,:) ! dim: (nproma,nlev,nblks_c,2)

INTEGER, INTENT(in), OPTIONAL ::  opt_slev ! optional vertical start level

INTEGER, INTENT(in), OPTIONAL ::  opt_elev ! optional vertical end level

INTEGER :: slev, elev     ! vertical start and end level
INTEGER :: jc, jk, jb
INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx, i_nchdom

INTEGER,  DIMENSION(:,:,:),   POINTER :: iidx, iblk

!-------------------------------------------------------------------------

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

iidx => ptr_patch%cells%edge_idx
iblk => ptr_patch%cells%edge_blk

! values for the blocking
i_nchdom   = MAX(1,ptr_patch%n_childdom)
i_startblk = ptr_patch%cells%start_blk(2,1)
i_endblk   = ptr_patch%cells%end_blk(min_rlcell,i_nchdom)


IF (ltimer) CALL timer_start(timer_intp)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk,jc)
DO jb = i_startblk, i_endblk

  CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
                     i_startidx, i_endidx, 2, min_rlcell)


#ifdef __LOOP_EXCHANGE
  DO jc = i_startidx, i_endidx
    DO jk = slev, elev
#else
!CDIR UNROLL=6
  DO jk = slev, elev
    DO jc = i_startidx, i_endidx
#endif

      p_uv_out(jc,jk,jb,1) =  &
        p_int%e_bln_c_u(jc,1,jb)*p_vn_in(iidx(jc,jb,1),jk,iblk(jc,jb,1)) + &
        p_int%e_bln_c_u(jc,2,jb)*p_vt_in(iidx(jc,jb,1),jk,iblk(jc,jb,1)) + &
        p_int%e_bln_c_u(jc,3,jb)*p_vn_in(iidx(jc,jb,2),jk,iblk(jc,jb,2)) + &
        p_int%e_bln_c_u(jc,4,jb)*p_vt_in(iidx(jc,jb,2),jk,iblk(jc,jb,2)) + &
        p_int%e_bln_c_u(jc,5,jb)*p_vn_in(iidx(jc,jb,3),jk,iblk(jc,jb,3)) + &
        p_int%e_bln_c_u(jc,6,jb)*p_vt_in(iidx(jc,jb,3),jk,iblk(jc,jb,3))

      p_uv_out(jc,jk,jb,2) =  &
        p_int%e_bln_c_v(jc,1,jb)*p_vn_in(iidx(jc,jb,1),jk,iblk(jc,jb,1)) + &
        p_int%e_bln_c_v(jc,2,jb)*p_vt_in(iidx(jc,jb,1),jk,iblk(jc,jb,1)) + &
        p_int%e_bln_c_v(jc,3,jb)*p_vn_in(iidx(jc,jb,2),jk,iblk(jc,jb,2)) + &
        p_int%e_bln_c_v(jc,4,jb)*p_vt_in(iidx(jc,jb,2),jk,iblk(jc,jb,2)) + &
        p_int%e_bln_c_v(jc,5,jb)*p_vn_in(iidx(jc,jb,3),jk,iblk(jc,jb,3)) + &
        p_int%e_bln_c_v(jc,6,jb)*p_vt_in(iidx(jc,jb,3),jk,iblk(jc,jb,3))

     ENDDO
   ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL

IF (ltimer) CALL timer_stop(timer_intp)


END SUBROUTINE edges2cells_vector

!------------------------------------------------------------------------
!
!
!>
!!  Computes  average of scalar fields from centers of cells to vertices.
!!
!!
!! @par Revision History
!! Developed  by Almut Gassmann, MPI-M (2009-01-28)
!!
SUBROUTINE cells2verts_scalar( p_cell_in, ptr_patch, c_int, p_vert_out,  &
  &                            opt_slev, opt_elev, opt_rlstart, opt_rlend )
!

TYPE(t_patch), TARGET, INTENT(in) :: ptr_patch

! cell based scalar input field
REAL(wp), INTENT(in) :: p_cell_in(:,:,:)   ! dim: (nproma,nlev,nblks_c)

! coefficients for interpolation
REAL(wp), INTENT(in) :: c_int(:,:,:)       ! dim: (nproma,9-i_cell_type,nblks_v)

INTEGER, INTENT(in), OPTIONAL :: opt_slev  ! optional vertical start level

INTEGER, INTENT(in), OPTIONAL :: opt_elev  ! optional vertical end level

! start and end values of refin_ctrl flag
INTEGER, INTENT(in), OPTIONAL ::  opt_rlstart, opt_rlend

! vertex based scalar output field
REAL(wp), INTENT(inout) :: p_vert_out(:,:,:) ! dim: (nproma,nlev,nblks_v)

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
  rl_end = min_rlvert
END IF

iidx => ptr_patch%verts%cell_idx
iblk => ptr_patch%verts%cell_blk

! values for the blocking
i_nchdom   = MAX(1,ptr_patch%n_childdom)
i_startblk = ptr_patch%verts%start_blk(rl_start,1)
i_endblk   = ptr_patch%verts%end_blk(rl_end,i_nchdom)


IF (ltimer) CALL timer_start(timer_intp)

IF (i_cell_type == 6) THEN
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jv,jk)
  DO jb = i_startblk, i_endblk

    CALL get_indices_v(ptr_patch, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, rl_start, rl_end)


#ifdef __LOOP_EXCHANGE
    DO jv = i_startidx, i_endidx
      DO jk = slev, elev
#else
!CDIR UNROLL=6
    DO jk = slev, elev
      DO jv = i_startidx, i_endidx
#endif

         p_vert_out(jv,jk,jb) =                       &
           c_int(jv,1,jb) * p_cell_in(iidx(jv,jb,1),jk,iblk(jv,jb,1)) + &
           c_int(jv,2,jb) * p_cell_in(iidx(jv,jb,2),jk,iblk(jv,jb,2)) + &
           c_int(jv,3,jb) * p_cell_in(iidx(jv,jb,3),jk,iblk(jv,jb,3))

      ENDDO
    ENDDO

  ENDDO
!$OMP END DO
!$OMP END PARALLEL
ELSE IF (i_cell_type == 3) THEN
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jv,jk)
  DO jb = i_startblk, i_endblk

    CALL get_indices_v(ptr_patch, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, rl_start, rl_end)


#ifdef __LOOP_EXCHANGE
    DO jv = i_startidx, i_endidx
      DO jk = slev, elev
#else
!CDIR UNROLL=6
    DO jk = slev, elev
      DO jv = i_startidx, i_endidx
#endif

         p_vert_out(jv,jk,jb) =                       &
           c_int(jv,1,jb) * p_cell_in(iidx(jv,jb,1),jk,iblk(jv,jb,1)) + &
           c_int(jv,2,jb) * p_cell_in(iidx(jv,jb,2),jk,iblk(jv,jb,2)) + &
           c_int(jv,3,jb) * p_cell_in(iidx(jv,jb,3),jk,iblk(jv,jb,3)) + &
           c_int(jv,4,jb) * p_cell_in(iidx(jv,jb,4),jk,iblk(jv,jb,4)) + &
           c_int(jv,5,jb) * p_cell_in(iidx(jv,jb,5),jk,iblk(jv,jb,5)) + &
           c_int(jv,6,jb) * p_cell_in(iidx(jv,jb,6),jk,iblk(jv,jb,6))

      ENDDO
    ENDDO

  ENDDO
!$OMP END DO
!$OMP END PARALLEL
ENDIF

IF (ltimer) CALL timer_stop(timer_intp)


END SUBROUTINE cells2verts_scalar
!------------------------------------------------------------------------
!
!
!>
!!  Computes  average of scalar fields from vertices to centers of cells.
!!
!!
!! @par Revision History
!! Developed  by Almut Gassmann, MPI-M (2009-01-28)
!!
SUBROUTINE verts2cells_scalar( p_vert_in, ptr_patch, c_int, p_cell_out,  &
  &                            opt_slev, opt_elev )
!

TYPE(t_patch), TARGET, INTENT(in) :: ptr_patch

! cell based scalar input field
REAL(wp), INTENT(in) :: p_vert_in(:,:,:)   ! dim: (nproma,nlev,nblks_v)

! coefficients for interpolation
REAL(wp), INTENT(in) :: c_int(:,:,:)       ! dim: (nproma,i_cell_type,nblks_c)

INTEGER, INTENT(in), OPTIONAL :: opt_slev  ! optional vertical start level

INTEGER, INTENT(in), OPTIONAL :: opt_elev  ! optional vertical end level

! vertex based scalar output field
REAL(wp), INTENT(inout) :: p_cell_out(:,:,:) ! dim: (nproma,nlev,nblks_c)

INTEGER :: slev, elev     ! vertical start and end level
INTEGER :: jk, jb, jc, nlen, nblks_c, npromz_c

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

iidx => ptr_patch%cells%vertex_idx
iblk => ptr_patch%cells%vertex_blk

! values for the blocking
nblks_c  = ptr_patch%nblks_int_c
npromz_c = ptr_patch%npromz_int_c


IF (ltimer) CALL timer_start(timer_intp)

IF (i_cell_type == 3) THEN
!$OMP PARALLEL
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

        p_cell_out(jc,jk,jb) =                                        &
          c_int(jc,1,jb)* p_vert_in(iidx(jc,jb,1),jk,iblk(jc,jb,1)) + &
          c_int(jc,2,jb)* p_vert_in(iidx(jc,jb,2),jk,iblk(jc,jb,2)) + &
          c_int(jc,3,jb)* p_vert_in(iidx(jc,jb,3),jk,iblk(jc,jb,3))

      ENDDO
    ENDDO

  END DO
!$OMP END DO
!$OMP END PARALLEL
ELSE IF (i_cell_type == 6) THEN
!$OMP PARALLEL
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

        p_cell_out(jc,jk,jb) =                                        &
          c_int(jc,1,jb)* p_vert_in(iidx(jc,jb,1),jk,iblk(jc,jb,1)) + &
          c_int(jc,2,jb)* p_vert_in(iidx(jc,jb,2),jk,iblk(jc,jb,2)) + &
          c_int(jc,3,jb)* p_vert_in(iidx(jc,jb,3),jk,iblk(jc,jb,3)) + &
          c_int(jc,4,jb)* p_vert_in(iidx(jc,jb,4),jk,iblk(jc,jb,4)) + &
          c_int(jc,5,jb)* p_vert_in(iidx(jc,jb,5),jk,iblk(jc,jb,5)) + &
          c_int(jc,6,jb)* p_vert_in(iidx(jc,jb,6),jk,iblk(jc,jb,6))

      ENDDO
    ENDDO

  ENDDO
!$OMP END DO
!$OMP END PARALLEL

ENDIF

IF (ltimer) CALL timer_stop(timer_intp)


END SUBROUTINE verts2cells_scalar

!------------------------------------------------------------------------
!
!
!>
!!  Computes interpolation from edges to rhombus centers
!!
!! Interpolation on quads. Coefficients are given on input.
!!
!! @par Revision History
!! Developed  by Almut Gassmann, MPI-M (2010-02-08)
!! Modified by Almut Gassmann, MPI-M (2010-10-18)
!! - center is always counted (remove option l_skip_center)
!!
SUBROUTINE edges2edges_scalar(p_edge_in,ptr_patch,c_int,p_edge_out)
!

TYPE(t_patch), TARGET, INTENT(in) :: ptr_patch

! cell based scalar input field
REAL(wp), INTENT(in) :: p_edge_in(:,:,:)   ! dim: (nproma,nlev,nblks_e)

! coefficients for interpolation
REAL(wp), INTENT(in) :: c_int(:,:,:)       ! dim: (4 or 5,nproma,nblks_e)

! vertex based scalar output field
REAL(wp), INTENT(inout) :: p_edge_out(:,:,:) ! dim: (nproma,nlev,nblks_e)

INTEGER :: nblks_e, npromz_e, nlen, jb, jk, je
INTEGER :: nlev              !< number of full levels

INTEGER,  DIMENSION(:,:,:),   POINTER :: iidx, iblk

!-----------------------------------------------------------------------

iidx => ptr_patch%edges%quad_idx
iblk => ptr_patch%edges%quad_blk

! values for the blocking
nblks_e  = ptr_patch%nblks_int_e
npromz_e = ptr_patch%npromz_int_e

! number of vertical levels
nlev = ptr_patch%nlev

IF (ltimer) CALL timer_start(timer_intp)

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
      DO jk = 1, nlev
#else
    DO jk = 1, nlev
      DO je = 1, nlen
#endif
        p_edge_out(je,jk,jb) =                                        &
        & + c_int(1,je,jb)* p_edge_in(iidx(je,jb,1),jk,iblk(je,jb,1)) &
        & + c_int(2,je,jb)* p_edge_in(iidx(je,jb,2),jk,iblk(je,jb,2)) &
        & + c_int(3,je,jb)* p_edge_in(iidx(je,jb,3),jk,iblk(je,jb,3)) &
        & + c_int(4,je,jb)* p_edge_in(iidx(je,jb,4),jk,iblk(je,jb,4)) &
        & + c_int(5,je,jb)* p_edge_in(je,jk,jb)
      ENDDO
    ENDDO
  END DO
!$OMP END DO
!$OMP END PARALLEL

IF (ltimer) CALL timer_stop(timer_intp)

END SUBROUTINE edges2edges_scalar

!-------------------------------------------------------------------------
!
!
!>
!! Computes the average of a cell-based variable.
!!
!! Computes the average of a cell-based variable
!! over its original location and the neighboring triangles.
!! Version with variable weighting coefficients, computed such that
!! linear horizontal gradients are not aliased into a checkerboard noise
!! input:  lives on centers of triangles
!! output: lives on centers of triangles
!!
!! @par Revision History
!!  developed by Guenther Zaengl, 2008-12-05
!!
SUBROUTINE cell_avg( psi_c, ptr_patch, avg_coeff, avg_psi_c,    &
  &                  opt_slev, opt_elev, opt_rlstart, opt_rlend )
!

!
!  patch on which computation is performed
!
TYPE(t_patch), TARGET, INTENT(in) :: ptr_patch
!
!  averaging coefficients
!
REAL(wp), INTENT(in) :: avg_coeff(:,:,:) ! dim: (nproma,nlev,nblks_c)

!
!  cell based variable before averaging
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
!   cell based variable after averaging
!
!REAL(wp), INTENT(out) ::  &
REAL(wp), INTENT(inout) ::  &
  &  avg_psi_c(:,:,:)  ! dim: (nproma,nlev,nblks_c)


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
  elev = ptr_patch%nlev
END IF

IF ( PRESENT(opt_rlstart) ) THEN
  IF (opt_rlstart == 1) THEN
    CALL finish ('mo_interpolation:cell_avg',  &
          &      'opt_rlstart must not be equal to 1')
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
!
! loop through all patch cells (and blocks)
!

  IF (ltimer) CALL timer_start(timer_intp)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx,jk)
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

        !  calculate the weighted average
        !
        avg_psi_c(jc,jk,jb) =  &
          &    psi_c(jc,jk,jb)                       * avg_coeff(jc,1,jb) &
          &  + psi_c(iidx(jc,jb,1),jk,iblk(jc,jb,1)) * avg_coeff(jc,2,jb) &
          &  + psi_c(iidx(jc,jb,2),jk,iblk(jc,jb,2)) * avg_coeff(jc,3,jb) &
          &  + psi_c(iidx(jc,jb,3),jk,iblk(jc,jb,3)) * avg_coeff(jc,4,jb)

      END DO !cell loop

    END DO !vertical levels loop

  END DO !block loop
!$OMP END DO
!$OMP END PARALLEL

  IF (ltimer) CALL timer_stop(timer_intp)


END SUBROUTINE cell_avg

!-------------------------------------------------------------------------


END MODULE mo_intp
