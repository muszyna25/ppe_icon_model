!>
!! This module provides the nonlinear advection routines that account for.
!!
!! This module provides the nonlinear advection routines that account for
!! @f$\vec{v}\cdot\nabla\vec{v} = \nabla \vec{v}^2 + \vec{\omega_a}\times\vec{v}@f$.
!! Two subroutines are specified,
!! <ul>
!!   <li>  kin_vel_rot computes the kinetic energy, the tangential wind and
!!         the relative vorticity eventually needed for the lamb term. For
!!         diagnosis the geographical horizontal wind comp., divergence and rotation
!!         at cell center are computed
!!   <li>  lamb_rot computes the tendency of the wind due to the second term in
!!         the above equation
!! </ul>
!! This subroutine sees the actual grid and will thus provide options for every
!! cell type (triangle, quadrilateral, hexagon)
!!
!! @par Revision History
!! Initial release by Almut Gassmann, MPI-M (2008-10-07)
!! Modification by Almut Gassmann, MPI-M (2009-01-29)
!! - conforming all scalar interpolation calls
!! Major revision by Almut Gassmann, MPI-M (2009-10-01)
!! - New bracket version for hexagons uses slightly different vorticities, so
!!   that the SICK instability can be prevented also in that case. Kinetic
!!   energy in the hexagonal case is ok.
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

MODULE mo_nonlinear_adv
!-------------------------------------------------------------------------
!
!    ProTeX FORTRAN source: Style 2
!    modified for ICON project, DWD/MPI-M 2008
!
!-------------------------------------------------------------------------
!
!
!

USE mo_kind,               ONLY: wp
USE mo_exception,          ONLY: finish
USE mo_model_domain,       ONLY: t_patch
USE mo_dynamics_config,    ONLY: idiv_method,lshallow_water
USE mo_io_config,          ONLY: l_outputtime, l_diagtime
USE mo_parallel_config,    ONLY: nproma
USE mo_intp_data_strc,     ONLY: t_int_state
USE mo_interpol_config,    ONLY: sick_a, sick_o
USE mo_intp_rbf,           ONLY: rbf_vec_interpol_edge
USE mo_intp,               ONLY: cells2verts_scalar,           &
  &                              edges2cells_scalar,           &
  &                              edges2verts_scalar,           &
  &                              edges2edges_scalar,           &
  &                              verts2edges_scalar,           &
  &                              verts2cells_scalar,           &
  &                              cell_avg
USE mo_interpol_config,    ONLY: l_corner_vort, i_cori_method
USE mo_icoham_dyn_types,   ONLY: t_hydro_atm_diag
USE mo_impl_constants,     ONLY: min_rledge
USE mo_math_divrot,        ONLY: div, rot_vertex
USE mo_loopindices,        ONLY: get_indices_e
USE mo_sync,               ONLY: SYNC_E, SYNC_V, sync_patch_array

IMPLICIT NONE

PRIVATE

PUBLIC :: kin_vel_rot, lamb_rot

CONTAINS

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!
!
!>
!! Computes the specific kinetic energy, the tangential velocity, the vorticity,.
!!
!! Computes the specific kinetic energy, the tangential velocity, the vorticity,
!! the geographical velocity vector, the divergence and the rotation at centers.
!!
!! @par Revision History
!! Initial release by Almut Gassmann, MPI-M (2008-10-09)
!! Modified by Almut Gassman, MPI-M (2009-02-02)
!! - rename this subroutine and add the computation of the vorticity here,
!! - add diagnosis for divergence and vorticity at CELL points
!! Modification by Almut Gassmann, MPI-M (2009-03-17)
!! - remove lraviart and reorder code
!!
SUBROUTINE kin_vel_rot (p_vn, pt_patch, pt_int, pt_diag)
!

REAL(wp), INTENT(IN)         :: p_vn(:,:,:)! normal velocity component at edges
TYPE(t_patch),TARGET,INTENT(IN):: pt_patch   ! actual patch
TYPE(t_int_state), TARGET,INTENT(IN):: pt_int! interpolation coeffs etc.

TYPE(t_hydro_atm_diag), INTENT(INOUT) :: pt_diag       ! diagnostic state

REAL(wp), DIMENSION (nproma,pt_patch%nlev,pt_patch%nblks_c) :: z_aux
REAL(wp), DIMENSION (nproma,pt_patch%nlev,pt_patch%nblks_e) :: z_kin_e
REAL(wp), DIMENSION (nproma,pt_patch%nlev,pt_patch%nblks_v) :: z_kin_v
INTEGER  :: nlen, nblks_c, nblks_e, npromz_c, npromz_e, jb, je, jk
INTEGER  :: i_startblk, i_startidx, i_endidx, i_rcstartlev
INTEGER  :: nlev              !< number of full levels

!-----------------------------------------------------------------------

! number of vertical levels
nlev = pt_patch%nlev

IF(pt_patch%geometry_info%cell_type == 3) THEN
  ! get relative vorticity at vertex
  CALL rot_vertex (p_vn, pt_patch, pt_int, pt_diag%rel_vort)
  ! This needs to be synced for lamb_rot
  CALL sync_patch_array(SYNC_V, pt_patch, pt_diag%rel_vort)
ELSE
  ! get relative vorticity at vertex
  CALL rot_vertex (p_vn, pt_patch, pt_int, pt_diag%rel_vort)
  CALL sync_patch_array(SYNC_V, pt_patch, pt_diag%rel_vort)
  CALL verts2edges_scalar(pt_diag%rel_vort, pt_patch, pt_int%tria_aw_rhom, pt_diag%rel_vort_e)
  ! This needs to be synced for lamb_rot
  CALL sync_patch_array(SYNC_E, pt_patch, pt_diag%rel_vort_e)
  ! for output we need rel_vort at vertex, but averaged
  CALL edges2verts_scalar(pt_diag%rel_vort_e, pt_patch, pt_int%e_1o3_v, pt_diag%rel_vort)
  CALL sync_patch_array(SYNC_V, pt_patch, pt_diag%rel_vort)
ENDIF

IF (l_outputtime) THEN
  ! get divergence
  CALL div (p_vn, pt_patch, pt_int, pt_diag%div)
  ! divergence is averaged according to the method used in the model
  IF(pt_patch%geometry_info%cell_type == 3) THEN
    SELECT CASE (idiv_method)
    CASE(2)
      z_aux(:,:,:) = pt_diag%div(:,:,:)
      CALL cell_avg( z_aux, pt_patch, pt_int%c_bln_avg, pt_diag%div)
    END SELECT
  ENDIF
  !----------------------------------------------------------------------
ENDIF

nblks_c   = pt_patch%nblks_c
npromz_c  = pt_patch%npromz_c
nblks_e   = pt_patch%nblks_e
npromz_e  = pt_patch%npromz_e


SELECT CASE (pt_patch%geometry_info%cell_type)
         !
CASE (3) ! for triangular grid (cell_type == 3)
         !
  !
  ! compute tangential velocity and kinetic energy
  !
  ! calculate the tangential components of the wind at edges
  CALL rbf_vec_interpol_edge( p_vn,            & ! normal wind comp.
    &                         pt_patch,        & ! patch
    &                         pt_int,          & ! interpolation state
    &                         pt_diag%vt       ) ! reconstr. tangential wind


  i_rcstartlev = 2  ! for rbf_vec_dim_edge==4
  i_startblk = pt_patch%edges%start_blk(i_rcstartlev,1)
!$OMP PARALLEL
!$OMP DO PRIVATE(jb, i_startidx,i_endidx,jk,je) ICON_OMP_DEFAULT_SCHEDULE
  DO jb = i_startblk, nblks_e

    CALL get_indices_e(pt_patch, jb, i_startblk, nblks_e, &
                       i_startidx, i_endidx, i_rcstartlev)

    DO jk = 1, nlev
      DO je = i_startidx, i_endidx
        ! calculate kinetic energy at edges from normal and tangential comp.
        z_kin_e(je,jk,jb) =0.5_wp*(pt_diag%vt(je,jk,jb)*pt_diag%vt(je,jk,jb)+&
                                   p_vn(je,jk,jb)*p_vn(je,jk,jb) )
      ENDDO
    ENDDO
  ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
  ! Bilinear interpolation of kinetic energy from the edges to the cells
  CALL edges2cells_scalar( z_kin_e, pt_patch, pt_int%e_bln_c_s,  &
    &                      pt_diag%e_kin, opt_rlstart=2 )

         !
CASE (6) ! for hexagonal/pentagonal grids (cell_type == 6)
         !
  !
  ! tangential velocity is not needed here
  !
  !
  ! kinetic energy is only formed of normal components
  !
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,nlen,jk,je) ICON_OMP_DEFAULT_SCHEDULE
  DO jb = 1, nblks_e
    IF (jb /= nblks_e) THEN
      nlen = nproma
    ELSE
      nlen = npromz_e
    ENDIF
    DO jk = 1, nlev
      DO je =1, nlen
        z_kin_e(je,jk,jb) = 0.5_wp * p_vn(je,jk,jb) * p_vn(je,jk,jb)
      ENDDO
    ENDDO
  ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  CALL edges2cells_scalar( z_kin_e,        & ! kin energy comp at edges
                             pt_patch,       & ! patch
                             pt_int%e_inn_c, & ! interpolation coeffs
                             pt_diag%e_kin   ) ! kine at centers

  IF(i_cori_method>= 2) THEN

    CALL edges2verts_scalar( z_kin_e,        & ! kin energy comp at edges
                             pt_patch,       & ! patch
                             pt_int%e_inn_v, & ! interpolation coeffs
                             z_kin_v         ) ! kine at verts
    CALL sync_patch_array(SYNC_V,pt_patch,z_kin_v)
    CALL verts2cells_scalar( z_kin_v,        & ! kin energy comp at verts
                             pt_patch,       & ! patch
                             pt_int%verts_aw_cells,&! interpolation coeffs
                             z_aux   ) ! kine at centers
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,nlen,jk) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = 1, nblks_c
      IF (jb /= nblks_c) THEN
        nlen = nproma
      ELSE
        nlen = npromz_c
      ENDIF
      DO jk = 1, nlev
        pt_diag%e_kin(1:nlen,jk,jb) = &
        &  sick_a*z_aux(1:nlen,jk,jb) + sick_o*pt_diag%e_kin(1:nlen,jk,jb)
      ENDDO
    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  ENDIF

  IF (l_outputtime) THEN
    !
    ! Reconstruction of geogr winds
    !
    CALL edges2cells_scalar(p_vn,               &! IN: normal wind
      &                     pt_patch,           &! patch
      &                     pt_int%hex_east,    &! interpolation state
      &                     pt_diag%u           )! reconstr. u wind
    CALL edges2cells_scalar(p_vn,               &! IN: normal wind
      &                     pt_patch,           &! patch
      &                     pt_int%hex_north,   &! interpolation state
      &                     pt_diag%v           )! reconstr. v wind
  ENDIF
           !
END SELECT ! cell_type
           !

END SUBROUTINE kin_vel_rot
!-------------------------------------------------------------------------
!
!
!>
!! This subroutine computes the vorticity flux term
!!
!! This subroutine computes the term @f$-frac{\omega_a}{\varrho}\vec{k}\times
!! \vec{v}@f$, abbreviated as the Lamb rotation. The result is the updated
!! tendency for the horizontal normal wind components.
!!
!! @par Revision History
!! Initial release by Almut Gassmann (2008-10-08)
!! Modification by Almut Gassmann, MPI-M (2009-01-29)
!! - further adjustments of interpolations to theory
!!
SUBROUTINE lamb_rot (pt_patch, pt_int, pt_diag, &
                     p_ddt_vn, opt_rlstart, opt_rlend)
!
IMPLICIT NONE

TYPE(t_patch), TARGET, INTENT(IN)   :: pt_patch       ! patch
TYPE(t_int_state),TARGET, INTENT(IN):: pt_int         ! interpolation state
TYPE(t_hydro_atm_diag),INTENT(INOUT)    :: pt_diag        ! diag state

REAL(wp), INTENT(INOUT) :: p_ddt_vn(:,:,:)     ! tendency for vn

INTEGER, INTENT(IN), OPTIONAL ::  &
  &  opt_rlstart, opt_rlend   ! start and end values of refin_ctrl flag

INTEGER   :: jk, jb, je, ji
INTEGER   :: nlen, npromz_e, nblks_e, nincr, nblks_c, npromz_c
INTEGER   :: i_startblk, i_endblk, i_startidx, i_endidx, i_nchdom
INTEGER   :: nlev              !< number of full levels
INTEGER   :: rl_start, rl_end
REAL(wp)  :: z_tmp_v    (nproma,pt_patch%nlev,pt_patch%nblks_v)
REAL(wp)  :: z_potvort_r(nproma,pt_patch%nlev,pt_patch%nblks_e)
REAL(wp)  :: z_potvort_c(nproma,pt_patch%nlev,pt_patch%nblks_c)
REAL(wp)  :: z_potvort_v(nproma,pt_patch%nlev,pt_patch%nblks_v)
REAL(wp)  :: z_u_c      (nproma,pt_patch%nlev,pt_patch%nblks_c)
REAL(wp)  :: z_v_c      (nproma,pt_patch%nlev,pt_patch%nblks_c)
REAL(wp)  :: z_u_e      (nproma,pt_patch%nlev,pt_patch%nblks_e)
REAL(wp)  :: z_v_e      (nproma,pt_patch%nlev,pt_patch%nblks_e)
REAL(wp)  :: z_vort_e   (nproma,pt_patch%nlev,pt_patch%nblks_e)
INTEGER, DIMENSION(:,:,:), POINTER :: ieidx, ieblk, icidx, icblk, ividx, ivblk
!-----------------------------------------------------------------------

! number of vertical levels
nlev = pt_patch%nlev

i_nchdom   = MAX(1,pt_patch%n_childdom)
                          !
SELECT CASE (pt_patch%geometry_info%cell_type)
         !
CASE (3) ! triangles (cell_type == 3)
         !

  ! Note: the rl_start and rl_end values here are only for edges
  IF ( PRESENT(opt_rlstart) ) THEN
    IF ((opt_rlstart >= 1) .AND. (opt_rlstart <= 2)) THEN
      CALL finish ('mo_nonlinear_adv:lamb_rot',  &
            &      'opt_rlstart must not be between 1 and 2')
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

  IF (l_diagtime .AND. lshallow_water) THEN

    ! delp at vertices (hexagons)
    CALL cells2verts_scalar(pt_diag%delp_c, pt_patch, &
                            pt_int%cells_aw_verts, pt_diag%delp_v)

  ENDIF

  ! simple averaged vorticity at exges
  CALL verts2edges_scalar(pt_diag%rel_vort, pt_patch, &
       pt_int%v_1o2_e, z_vort_e, opt_rlstart=3)

  i_startblk = pt_patch%edges%start_blk(rl_start,1)
  i_endblk   = pt_patch%edges%end_blk(rl_end,i_nchdom)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk,je) ICON_OMP_DEFAULT_SCHEDULE
  DO jb = i_startblk, i_endblk

    CALL get_indices_e(pt_patch, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, rl_start, rl_end)

    DO jk = 1, nlev
      DO je = i_startidx, i_endidx
        p_ddt_vn(je,jk,jb) = p_ddt_vn(je,jk,jb)   &
          -(z_vort_e(je,jk,jb)+pt_patch%edges%f_e(je,jb))&
          *pt_diag%vt(je,jk,jb)
      ENDDO
    ENDDO
    !!HW: minus sign here because the tangential direction in the code is
    !!    in fact the opposite from what is written in BR05.

  ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
          !
CASE (6)  ! hexagons/pentagons (cell_type == 6)
          !
  !
  ! Absolute Potential Vorticity at rhombi
  !
  ! no grid refinement in hexagonal model
  nblks_e   = pt_patch%nblks_e
  npromz_e  = pt_patch%npromz_e
  ! loop over blocks and verts
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,nlen,jk) ICON_OMP_DEFAULT_SCHEDULE
  DO jb = 1, nblks_e
    IF (jb /= nblks_e) THEN
      nlen = nproma
    ELSE
      nlen = npromz_e
    ENDIF
    DO jk = 1, nlev
      z_potvort_r(1:nlen,jk,jb) = (pt_diag%rel_vort_e(1:nlen,jk,jb)&
                                  +pt_patch%edges%f_e(1:nlen,jb))&
                                  /pt_diag%delp_e(1:nlen,jk,jb)
    ENDDO
  ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  !
  ! tendencies from vorticity flux term
  !
  SELECT CASE (i_cori_method)

  CASE (1,2)
    ieidx => pt_int%heli_vn_idx
    ieblk => pt_int%heli_vn_blk
    SELECT CASE (i_cori_method)
    CASE (1)
      nincr = 14
    CASE (2)
      nincr = 10
    END SELECT
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,nlen,jk,ji,je) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = 1, nblks_e
      IF (jb /= nblks_e) THEN
        nlen = nproma
      ELSE
        nlen = npromz_e
      ENDIF
#ifdef __LOOP_EXCHANGE
      DO je = 1, nlen
        DO ji = 1, nincr
          DO jk = 1, nlev
#else
      DO jk = 1, nlev
        DO ji = 1, nincr
          DO je = 1, nlen
#endif
            p_ddt_vn(je,jk,jb) = p_ddt_vn(je,jk,jb)                   &
            & +pt_int%heli_coeff        (ji,je,jb)                    &
            & *pt_diag%mass_flux_e(ieidx(ji,je,jb),jk,ieblk(ji,je,jb))&
            & *(z_potvort_r       (ieidx(ji,je,jb),jk,ieblk(ji,je,jb))&
            &  +z_potvort_r       (je             ,jk,jb             ))
          ENDDO
        ENDDO
      ENDDO
    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  CASE (3)

    nblks_c   = pt_patch%nblks_c
    npromz_c  = pt_patch%npromz_c

    ieidx => pt_patch%edges%quad_idx
    ieblk => pt_patch%edges%quad_blk
    icidx => pt_patch%edges%cell_idx
    icblk => pt_patch%edges%cell_blk

    ! first, reconstruct mass flux vectors at centers of rhombi and hexagons
    CALL edges2cells_scalar(pt_diag%mass_flux_e,pt_patch,pt_int%hex_east  ,z_u_c)
    CALL edges2cells_scalar(pt_diag%mass_flux_e,pt_patch,pt_int%hex_north ,z_v_c)
    CALL edges2edges_scalar(pt_diag%mass_flux_e,pt_patch,pt_int%quad_east ,z_u_e)
    CALL edges2edges_scalar(pt_diag%mass_flux_e,pt_patch,pt_int%quad_north,z_v_e)

    ! second, average absolute potential vorticity from rhombi to centers
    IF (l_corner_vort) THEN
      CALL edges2verts_scalar(z_potvort_r,pt_patch,pt_int%e_1o3_v       ,z_tmp_v)
      CALL sync_patch_array(SYNC_V,pt_patch,z_tmp_v)
      CALL verts2cells_scalar(z_tmp_v    ,pt_patch,pt_int%verts_aw_cells,z_potvort_c)
    ELSE
      CALL edges2cells_scalar(z_potvort_r,pt_patch,pt_int%e_aw_c,z_potvort_c)
    ENDIF

    ! third, multiply the absolute vorticities with the velocities,
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,nlen,jk) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = 1, nblks_e
      IF (jb /= nblks_e) THEN
        nlen = nproma
      ELSE
        nlen = npromz_e
      ENDIF
      DO jk = 1,nlev
        z_u_e(1:nlen,jk,jb) = z_u_e(1:nlen,jk,jb)*z_potvort_r(1:nlen,jk,jb)
        z_v_e(1:nlen,jk,jb) = z_v_e(1:nlen,jk,jb)*z_potvort_r(1:nlen,jk,jb)
      ENDDO
    ENDDO
!$OMP END DO NOWAIT
!$OMP DO PRIVATE(jb,nlen,jk) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = 1, nblks_c
      IF (jb /= nblks_c) THEN
        nlen = nproma
      ELSE
        nlen = npromz_c
      ENDIF
      DO jk = 1,nlev
        z_u_c(1:nlen,jk,jb) = z_u_c(1:nlen,jk,jb)*z_potvort_c(1:nlen,jk,jb)
        z_v_c(1:nlen,jk,jb) = z_v_c(1:nlen,jk,jb)*z_potvort_c(1:nlen,jk,jb)
      ENDDO
    ENDDO
!$OMP END DO
    
! fourth, compute vorticity flux term
!$OMP DO PRIVATE(jb,nlen,jk,je) ICON_OMP_DEFAULT_SCHEDULE
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
          p_ddt_vn(je,jk,jb) =  p_ddt_vn(je,jk,jb)             &
          & + pt_int%heli_coeff( 1,je,jb)*z_v_c(icidx(je,jb,1),jk,icblk(je,jb,1))&
          & + pt_int%heli_coeff( 2,je,jb)*z_u_c(icidx(je,jb,1),jk,icblk(je,jb,1))&
          & + pt_int%heli_coeff( 3,je,jb)*z_v_c(icidx(je,jb,2),jk,icblk(je,jb,2))&
          & + pt_int%heli_coeff( 4,je,jb)*z_u_c(icidx(je,jb,2),jk,icblk(je,jb,2))&
          & + pt_int%heli_coeff( 5,je,jb)*z_v_e(je            ,jk,jb            )&
          & + pt_int%heli_coeff( 6,je,jb)*z_u_e(je            ,jk,jb            )&
          & + pt_int%heli_coeff( 7,je,jb)*z_v_e(ieidx(je,jb,1),jk,ieblk(je,jb,1))&
          & + pt_int%heli_coeff( 8,je,jb)*z_u_e(ieidx(je,jb,1),jk,ieblk(je,jb,1))&
          & + pt_int%heli_coeff( 9,je,jb)*z_v_e(ieidx(je,jb,2),jk,ieblk(je,jb,2))&
          & + pt_int%heli_coeff(10,je,jb)*z_u_e(ieidx(je,jb,2),jk,ieblk(je,jb,2))&
          & + pt_int%heli_coeff(11,je,jb)*z_v_e(ieidx(je,jb,3),jk,ieblk(je,jb,3))&
          & + pt_int%heli_coeff(12,je,jb)*z_u_e(ieidx(je,jb,3),jk,ieblk(je,jb,3))&
          & + pt_int%heli_coeff(13,je,jb)*z_v_e(ieidx(je,jb,4),jk,ieblk(je,jb,4))&
          & + pt_int%heli_coeff(14,je,jb)*z_u_e(ieidx(je,jb,4),jk,ieblk(je,jb,4))
        ENDDO
      ENDDO
    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  CASE (4)

    ieidx => pt_patch%edges%quad_idx
    ieblk => pt_patch%edges%quad_blk
    icidx => pt_patch%edges%cell_idx
    icblk => pt_patch%edges%cell_blk
    ividx => pt_patch%edges%vertex_idx
    ivblk => pt_patch%edges%vertex_blk

    ! first, reconstruct mass flux vectors at centers of rhombi and hexagons
    CALL edges2cells_scalar(pt_diag%mass_flux_e,pt_patch,pt_int%hex_east  ,z_u_c)
    CALL edges2cells_scalar(pt_diag%mass_flux_e,pt_patch,pt_int%hex_north ,z_v_c)
    CALL edges2edges_scalar(pt_diag%mass_flux_e,pt_patch,pt_int%quad_east ,z_u_e)
    CALL edges2edges_scalar(pt_diag%mass_flux_e,pt_patch,pt_int%quad_north,z_v_e)

    ! second, average absolute potential vorticity from rhombi to centers
    CALL edges2verts_scalar(z_potvort_r,pt_patch,pt_int%e_1o3_v,z_potvort_v)
    CALL sync_patch_array(SYNC_V,pt_patch,z_potvort_v)
    IF (l_corner_vort) THEN
      CALL verts2cells_scalar(z_potvort_v,pt_patch,pt_int%verts_aw_cells,z_potvort_c)
    ELSE
      CALL edges2cells_scalar(z_potvort_r,pt_patch,pt_int%e_aw_c,z_potvort_c)
    ENDIF

!$OMP PARALLEL
    ! third, compute vorticity flux term
!$OMP DO PRIVATE(jb,nlen,jk,je) ICON_OMP_DEFAULT_SCHEDULE
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
          p_ddt_vn(je,jk,jb) =  p_ddt_vn(je,jk,jb)             &
          & + (pt_int%heli_coeff( 1,je,jb)*z_v_c(icidx(je,jb,1),jk,icblk(je,jb,1)) &
          &   +pt_int%heli_coeff( 2,je,jb)*z_u_c(icidx(je,jb,1),jk,icblk(je,jb,1)))&
          &                         *z_potvort_c(icidx(je,jb,1),jk,icblk(je,jb,1)) &
          & + (pt_int%heli_coeff( 3,je,jb)*z_v_c(icidx(je,jb,2),jk,icblk(je,jb,2)) &
          &   +pt_int%heli_coeff( 4,je,jb)*z_u_c(icidx(je,jb,2),jk,icblk(je,jb,2)))&
          &                         *z_potvort_c(icidx(je,jb,2),jk,icblk(je,jb,2)) &
          & + (pt_int%heli_coeff( 5,je,jb)*z_v_e(je            ,jk,jb            ) &
          & +  pt_int%heli_coeff( 6,je,jb)*z_u_e(je            ,jk,jb            ))&
          &                 *0.5_wp*(z_potvort_v(ividx(je,jb,1),jk,ivblk(je,jb,1)) &
          &                         +z_potvort_v(ividx(je,jb,2),jk,ivblk(je,jb,2)))&
          & + (pt_int%heli_coeff( 7,je,jb)*z_v_e(ieidx(je,jb,1),jk,ieblk(je,jb,1)) &
          &   +pt_int%heli_coeff( 8,je,jb)*z_u_e(ieidx(je,jb,1),jk,ieblk(je,jb,1)) &
          &   +pt_int%heli_coeff( 9,je,jb)*z_v_e(ieidx(je,jb,2),jk,ieblk(je,jb,2)) &
          &   +pt_int%heli_coeff(10,je,jb)*z_u_e(ieidx(je,jb,2),jk,ieblk(je,jb,2)))&
          &                         *z_potvort_v(ividx(je,jb,1),jk,ivblk(je,jb,1)) &
          & + (pt_int%heli_coeff(11,je,jb)*z_v_e(ieidx(je,jb,3),jk,ieblk(je,jb,3)) &
          &   +pt_int%heli_coeff(12,je,jb)*z_u_e(ieidx(je,jb,3),jk,ieblk(je,jb,3)) &
          &   +pt_int%heli_coeff(13,je,jb)*z_v_e(ieidx(je,jb,4),jk,ieblk(je,jb,4)) &
          &   +pt_int%heli_coeff(14,je,jb)*z_u_e(ieidx(je,jb,4),jk,ieblk(je,jb,4)))&
          &                         *z_potvort_v(ividx(je,jb,2),jk,ivblk(je,jb,2))   
        ENDDO
      ENDDO
    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL


  END SELECT

END SELECT ! cell_type

END SUBROUTINE lamb_rot

END MODULE mo_nonlinear_adv
