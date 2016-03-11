!>
!! Contains the interpolation routines needed for grid refinement.
!!
!! These had originally been included in mo_grf_interpolation but then were
!! packed into a separate module to clean up the code
!!
!! @par Revision History
!! Created by Guenther Zaengl, DWD (2009-02-09)
!! Modification by Guenther Zaengl, DWD (2009-06-22)
!! - preparation for generalized grid refinement (affects all subroutines)
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

MODULE mo_grf_nudgintp
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
USE mo_intp_data_strc,      ONLY: t_int_state
USE mo_impl_constants_grf,  ONLY: grf_nudgintp_start_c, grf_nudgintp_start_e, &
                                  grf_nudge_start_e
USE mo_parallel_config,     ONLY: nproma
USE mo_loopindices,         ONLY: get_indices_c, get_indices_e, get_indices_v

USE mo_grf_intp_data_strc


IMPLICIT NONE

PRIVATE

PUBLIC :: interpol_vec_nudging, interpol_scal_nudging

CONTAINS

!
!>
!! Parent-to-child interpolation of normal velocity components needed for boundary nudging
!! Interpolation methods are as in interpol2_vec_grf
!!
!!
!! @par Revision History
!! Developed  by Guenther Zaengl, DWD (2010-06-17)
!!
SUBROUTINE interpol_vec_nudging (ptr_pp, ptr_pc, ptr_int, ptr_grf, ptr_grfc,   &
  &                              i_chidx, nshift, istart_blk, p_vn_in, p_vn_out)

TYPE(t_patch), TARGET, INTENT(in) :: ptr_pp
TYPE(t_patch), TARGET, INTENT(in) :: ptr_pc

! Block index at which p_vn_in starts
INTEGER, INTENT(IN) :: istart_blk

! input normal components of vectors at edge midpoints
REAL(wp), INTENT(IN) :: p_vn_in(:,:,istart_blk:) ! dim: (nproma,nlev,nblks_e)

! Indices of source points and interpolation coefficients
TYPE(t_gridref_single_state),   TARGET, INTENT(IN)    :: ptr_grf
TYPE(t_gridref_state),          TARGET, INTENT(IN)    :: ptr_grfc
TYPE(t_int_state),              TARGET, INTENT(IN)    :: ptr_int

! child domain index as seen from parent domain
INTEGER, INTENT(IN) :: i_chidx

! number of levels by which vertical shifting between input and output fields is needed
INTEGER, INTENT(IN) :: nshift

! reconstructed edge-normal wind component
REAL(wp),INTENT(INOUT) :: p_vn_out(:,:,:) ! dim: (nproma,nlev_c,nblks_e)


INTEGER :: jb, jk, je, jv            ! loop indices
INTEGER :: js                        ! shift parameter
INTEGER :: i_startblk                ! start block
INTEGER :: i_endblk                  ! end index
INTEGER :: i_startidx                ! start index
INTEGER :: i_endidx                  ! end index
INTEGER :: i_nchdom                  ! number of child domains

INTEGER :: nlev_c       !< number of vertical full levels (child domain)

REAL(wp) :: dvn_tang(nproma)

! Auxiliary fields
REAL(wp) :: vn_aux(nproma,ptr_pc%nlev,ptr_pp%edges%start_block(grf_nudgintp_start_e):&
                   MAX(ptr_pp%edges%start_block(grf_nudgintp_start_e),               &
                       ptr_pp%edges%end_block(min_rledge_int)),4)

REAL(wp), DIMENSION(nproma,ptr_pc%nlev,ptr_pp%verts%start_block(grf_nudgintp_start_c):&
                    MAX(ptr_pp%verts%start_block(grf_nudgintp_start_c),               &
                        ptr_pp%verts%end_block(min_rlvert_int))) :: u_vert, v_vert

! Pointers to index fields
INTEGER,  DIMENSION(:,:,:), POINTER :: iidx_2a, iblk_2a, iidx_2b, iblk_2b, ividx, ivblk, &
                                       irvidx, irvblk
INTEGER,  DIMENSION(:,:),   POINTER :: ipeidx, ipeblk, ipcidx

REAL(wp), DIMENSION(:,:,:,:), POINTER :: ptr_rvcoeff
REAL(wp), DIMENSION(:,:,:), POINTER   :: ptr_dist
!-----------------------------------------------------------------------

! Set pointers to index and coefficient fields

iidx_2a => ptr_grf%grf_vec_ind_2a
iblk_2a => ptr_grf%grf_vec_blk_2a
iidx_2b => ptr_grf%grf_vec_ind_2b
iblk_2b => ptr_grf%grf_vec_blk_2b

! vertices of edges
ividx => ptr_pp%edges%vertex_idx
ivblk => ptr_pp%edges%vertex_blk

! for RBF interpolation to vertices
irvidx => ptr_int%rbf_vec_idx_v
irvblk => ptr_int%rbf_vec_blk_v
ptr_rvcoeff => ptr_int%rbf_vec_coeff_v

! normalized distances from parent edges to child edges
ptr_dist  => ptr_grf%grf_dist_pe2ce

! parent edge indices as seen from the child grid level
ipeidx => ptr_pc%edges%parent_loc_idx
ipeblk => ptr_pc%edges%parent_loc_blk
ipcidx => ptr_pc%edges%pc_idx

! number of vertical full levels (child domain)
nlev_c = ptr_pc%nlev

! Shift parameter
js = nshift

! Number of child domains of nested domain
i_nchdom = MAX(1,ptr_pc%n_childdom)

!$OMP PARALLEL PRIVATE (i_startblk,i_endblk)

! Start and end blocks for which vector reconstruction at vertices is needed
! Note: the use of grf_nudgintp_start_c (=-3) is intended
i_startblk = ptr_pp%verts%start_block(grf_nudgintp_start_c)
i_endblk   = ptr_pp%verts%end_block(min_rlvert_int)

!$OMP DO PRIVATE(jb,jk,jv,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
DO jb = i_startblk, i_endblk

  CALL get_indices_v(ptr_pp, jb, i_startblk, i_endblk, &
   i_startidx, i_endidx, grf_nudgintp_start_c, min_rlvert_int)

#ifdef __LOOP_EXCHANGE
  DO jv = i_startidx, i_endidx
    DO jk = 1, nlev_c
#else
!CDIR UNROLL=6
  DO jk = 1, nlev_c
    DO jv = i_startidx, i_endidx
#endif

      u_vert(jv,jk,jb) =  &
        ptr_rvcoeff(1,1,jv,jb)*p_vn_in(irvidx(1,jv,jb),jk+js,irvblk(1,jv,jb)) + &
        ptr_rvcoeff(2,1,jv,jb)*p_vn_in(irvidx(2,jv,jb),jk+js,irvblk(2,jv,jb)) + &
        ptr_rvcoeff(3,1,jv,jb)*p_vn_in(irvidx(3,jv,jb),jk+js,irvblk(3,jv,jb)) + &
        ptr_rvcoeff(4,1,jv,jb)*p_vn_in(irvidx(4,jv,jb),jk+js,irvblk(4,jv,jb)) + &
        ptr_rvcoeff(5,1,jv,jb)*p_vn_in(irvidx(5,jv,jb),jk+js,irvblk(5,jv,jb)) + &
        ptr_rvcoeff(6,1,jv,jb)*p_vn_in(irvidx(6,jv,jb),jk+js,irvblk(6,jv,jb))
      v_vert(jv,jk,jb) =  &
        ptr_rvcoeff(1,2,jv,jb)*p_vn_in(irvidx(1,jv,jb),jk+js,irvblk(1,jv,jb)) + &
        ptr_rvcoeff(2,2,jv,jb)*p_vn_in(irvidx(2,jv,jb),jk+js,irvblk(2,jv,jb)) + &
        ptr_rvcoeff(3,2,jv,jb)*p_vn_in(irvidx(3,jv,jb),jk+js,irvblk(3,jv,jb)) + &
        ptr_rvcoeff(4,2,jv,jb)*p_vn_in(irvidx(4,jv,jb),jk+js,irvblk(4,jv,jb)) + &
        ptr_rvcoeff(5,2,jv,jb)*p_vn_in(irvidx(5,jv,jb),jk+js,irvblk(5,jv,jb)) + &
        ptr_rvcoeff(6,2,jv,jb)*p_vn_in(irvidx(6,jv,jb),jk+js,irvblk(6,jv,jb))

    ENDDO
  ENDDO
ENDDO
!$OMP END DO

! Start and end blocks for which interpolation is needed
i_startblk = ptr_pp%edges%start_block(grf_nudgintp_start_e)
i_endblk   = ptr_pp%edges%end_block(min_rledge_int)

!$OMP DO PRIVATE(jb,jk,je,i_startidx,i_endidx,dvn_tang) ICON_OMP_DEFAULT_SCHEDULE
DO jb =  i_startblk, i_endblk

  CALL get_indices_e(ptr_pp, jb, i_startblk, i_endblk, &
   i_startidx, i_endidx, grf_nudgintp_start_e, min_rledge_int)

! child edges 1 and 2
#ifdef __LOOP_EXCHANGE
  DO je = i_startidx, i_endidx
    DO jk = 1, nlev_c
#else
!CDIR UNROLL=6
  DO jk = 1, nlev_c
    DO je = i_startidx, i_endidx
#endif

      dvn_tang(je) = u_vert(ividx(je,jb,2),jk,ivblk(je,jb,2)) *    &
                     ptr_pp%edges%primal_normal_vert(je,jb,2)%v1 + &
                     v_vert(ividx(je,jb,2),jk,ivblk(je,jb,2)) *    &
                     ptr_pp%edges%primal_normal_vert(je,jb,2)%v2 - &
                    (u_vert(ividx(je,jb,1),jk,ivblk(je,jb,1)) *    &
                     ptr_pp%edges%primal_normal_vert(je,jb,1)%v1 + &
                     v_vert(ividx(je,jb,1),jk,ivblk(je,jb,1)) *    &
                     ptr_pp%edges%primal_normal_vert(je,jb,1)%v2 )

      vn_aux(je,jk,jb,1) = p_vn_in(je,jk+js,jb) + dvn_tang(je)*ptr_dist(je,1,jb)
      vn_aux(je,jk,jb,2) = p_vn_in(je,jk+js,jb) + dvn_tang(je)*ptr_dist(je,2,jb)
    ENDDO
  ENDDO

! child edge 3
#ifdef __LOOP_EXCHANGE
  DO je = i_startidx, i_endidx
    DO jk = 1, nlev_c
#else
!CDIR UNROLL=6
  DO jk = 1, nlev_c
    DO je = i_startidx, i_endidx
#endif
      vn_aux(je,jk,jb,3) = ptr_grf%grf_vec_coeff_2a(1,je,jb) * &
        p_vn_in(iidx_2a(je,1,jb),jk+js,iblk_2a(je,1,jb)) +     &
        ptr_grf%grf_vec_coeff_2a(2,je,jb) *                    &
        p_vn_in(iidx_2a(je,2,jb),jk+js,iblk_2a(je,2,jb)) +     &
        ptr_grf%grf_vec_coeff_2a(3,je,jb) *                    &
        p_vn_in(iidx_2a(je,3,jb),jk+js,iblk_2a(je,3,jb)) +     &
        ptr_grf%grf_vec_coeff_2a(4,je,jb) *                    &
        p_vn_in(iidx_2a(je,4,jb),jk+js,iblk_2a(je,4,jb)) +     &
        ptr_grf%grf_vec_coeff_2a(5,je,jb) *                    &
        p_vn_in(iidx_2a(je,5,jb),jk+js,iblk_2a(je,5,jb))
    ENDDO
  ENDDO

! child edge 4
#ifdef __LOOP_EXCHANGE
  DO je = i_startidx, i_endidx
    DO jk = 1, nlev_c
#else
!CDIR UNROLL=6
  DO jk = 1, nlev_c
    DO je = i_startidx, i_endidx
#endif
      vn_aux(je,jk,jb,4) = ptr_grf%grf_vec_coeff_2b(1,je,jb) * &
        p_vn_in(iidx_2b(je,1,jb),jk+js,iblk_2b(je,1,jb)) +     &
        ptr_grf%grf_vec_coeff_2b(2,je,jb) *                    &
        p_vn_in(iidx_2b(je,2,jb),jk+js,iblk_2b(je,2,jb)) +     &
        ptr_grf%grf_vec_coeff_2b(3,je,jb) *                    &
        p_vn_in(iidx_2b(je,3,jb),jk+js,iblk_2b(je,3,jb)) +     &
        ptr_grf%grf_vec_coeff_2b(4,je,jb) *                    &
        p_vn_in(iidx_2b(je,4,jb),jk+js,iblk_2b(je,4,jb)) +     &
        ptr_grf%grf_vec_coeff_2b(5,je,jb) *                    &
        p_vn_in(iidx_2b(je,5,jb),jk+js,iblk_2b(je,5,jb))
    ENDDO
  ENDDO

ENDDO ! blocks
!$OMP END DO

! Store results in p_vn_out

! Note: flow control has to done from the child level here to avoid conflicts with
! interpolation of boundary tendencies (parent-child differences for nudging and
! boundary tendencies are stored in the same fields)

! For MPI parallelization, parent-edge and parent-child indices have to be mapped
! to the feedback-parent grid level

! Start and end blocks for which interpolation is needed
i_startblk = ptr_pc%edges%start_block(grf_nudge_start_e)
i_endblk   = ptr_pc%edges%end_block(min_rledge_int)

!$OMP DO PRIVATE(jb,jk,je,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
DO jb =  i_startblk, i_endblk

  CALL get_indices_e(ptr_pc, jb, i_startblk, i_endblk, &
    i_startidx, i_endidx, grf_nudge_start_e, min_rledge_int)

#ifdef __LOOP_EXCHANGE
  DO je = i_startidx, i_endidx
      DO jk = 1, nlev_c
        p_vn_out(je,jk,jb) = vn_aux(ipeidx(je,jb),jk,ipeblk(je,jb),ipcidx(je,jb))
      ENDDO
  ENDDO
#else
  DO jk = 1, nlev_c
    DO je = i_startidx, i_endidx
        p_vn_out(je,jk,jb) = vn_aux(ipeidx(je,jb),jk,ipeblk(je,jb),ipcidx(je,jb))
    ENDDO
  ENDDO
#endif
ENDDO
!$OMP END DO NOWAIT

!$OMP END PARALLEL


END SUBROUTINE interpol_vec_nudging


!-------------------------------------------------------------------------
!
!>
!! Parent-to-child interpolation of scalar fields needed for boundary nudging
!! Interpolation methods are as in interpol_scal_grf
!!
!! @par Revision History
!! Developed  by Guenther Zaengl, DWD (2010-06-16)
!!
SUBROUTINE interpol_scal_nudging (ptr_pp, ptr_int, ptr_grf, i_chidx, nshift,     &
                                  nfields, istart_blk, f3din1, f3dout1, f3din2,  &
                                  f3dout2, f3din3, f3dout3, f3din4, f3dout4,     &
                                  f3din5, f3dout5, f4din, f4dout,                &
                                  llimit_nneg, rlimval, overshoot_fac,           &
                                  opt_l_enabled)
!
TYPE(t_patch), TARGET, INTENT(in) :: ptr_pp

! Block index at which the input fields start
INTEGER, INTENT(IN) :: istart_blk

! Indices of source points and interpolation coefficients
TYPE(t_gridref_single_state),   TARGET, INTENT(IN)    ::  ptr_grf
TYPE(t_int_state),              TARGET, INTENT(IN)    ::  ptr_int

! child domain index as seen from parent domain
INTEGER, INTENT(IN) :: i_chidx

! number of levels by which vertical shifting between input and output fields is needed
INTEGER, INTENT(IN) :: nshift

! number of fields provided on input (needed for aux fields and pointer allocation)
INTEGER, INTENT(IN) :: nfields

! input scalar fields at cell points (up to five 3D fields or one 4D field)
REAL(wp), INTENT(IN), OPTIONAL, TARGET ::  & ! dim: (nproma,nlev,nblks_c)
  f3din1(:,:,istart_blk:), f3din2(:,:,istart_blk:), f3din3(:,:,istart_blk:), &
  f3din4(:,:,istart_blk:), f3din5(:,:,istart_blk:),                          &
  f4din(:,:,:,:) ! Note: for f4din, specifying the lower bound of blocks did not work;
                 ! in this case, the full field needs to be allocated in the calling routine

! reconstructed scalar output fields (up to five 3D fields or one 4D field)
REAL(wp), INTENT(INOUT), OPTIONAL, TARGET ::  & ! dim: (nproma,nlev_c,nblks_c)
  f3dout1(:,:,:), f3dout2(:,:,:), f3dout3(:,:,:), f3dout4(:,:,:), f3dout5(:,:,:), &
  f4dout(:,:,:,:)

! logical switch to activate gradient limiter
LOGICAL, INTENT(IN), OPTIONAL :: llimit_nneg(nfields)

! value to which the variable value is limited
REAL(wp), INTENT(IN), OPTIONAL :: rlimval(nfields)

! factor up to to which overshooting is allowed
REAL(wp), INTENT(IN), OPTIONAL :: overshoot_fac

! LOGICAL field: skip level if .FALSE.
LOGICAL, OPTIONAL, INTENT(IN) :: opt_l_enabled(:)

INTEGER :: jb, jk, jc, jn, n         ! loop indices
INTEGER :: js                        ! shift parameter
INTEGER :: i_startblk                ! start block
INTEGER :: i_endblk                  ! end index
INTEGER :: i_startidx                ! start index
INTEGER :: i_endidx                  ! end index
INTEGER :: elev                      ! end index of vertical loop

LOGICAL :: l4d                       ! 4D field is provided as input

! local variables corresponding to optional input variables
LOGICAL :: l_limit_nneg(nfields)
REAL(wp):: r_limval(nfields)

! Auxiliary fields
REAL(wp), DIMENSION(nproma,MAX(32,ptr_pp%nlevp1)) :: grad_x, grad_y, maxval_neighb, minval_neighb
REAL(wp) :: h_aux(nproma,MAX(32,ptr_pp%nlevp1),                       &
                  ptr_pp%cells%start_block(grf_nudgintp_start_c):     &
                  MAX(ptr_pp%cells%start_block(grf_nudgintp_start_c), &
                      ptr_pp%cells%end_block(min_rlcell_int)),4,nfields)

REAL(wp) :: limfac1, limfac2, limfac, min_expval(nproma), max_expval(nproma), ovsht_fac, r_ovsht_fac, &
            relaxed_minval, relaxed_maxval

REAL(wp), PARAMETER :: epsi = 1.e-75_wp
! Pointers to index and coefficient fields
INTEGER,  DIMENSION(:,:,:),   POINTER :: iidx, iblk, ichcidx, ichcblk
REAL(wp), DIMENSION(:,:,:,:), POINTER :: ptr_coeff, ptr_dist

! Allocatable pointer to input and output fields
TYPE t_fieldptr
  REAL(wp), POINTER :: fld(:,:,:)
END TYPE t_fieldptr
TYPE(t_fieldptr) :: p_in(nfields), p_out(nfields)

    LOGICAL, ALLOCATABLE :: l_enabled(:)
    LOGICAL              :: all_enabled

!-----------------------------------------------------------------------

IF (PRESENT(f4din)) THEN
  DO n = 1, nfields
    p_in(n)%fld  => f4din(:,:,:,n)
    p_out(n)%fld => f4dout(:,:,:,n)
  ENDDO
  l4d = .TRUE.
ELSE
  IF (PRESENT(f3din1)) THEN
    p_in(1)%fld  => f3din1
    p_out(1)%fld => f3dout1
  ENDIF
  IF (PRESENT(f3din2)) THEN
    p_in(2)%fld  => f3din2
    p_out(2)%fld => f3dout2
  ENDIF
  IF (PRESENT(f3din3)) THEN
    p_in(3)%fld  => f3din3
    p_out(3)%fld => f3dout3
  ENDIF
  IF (PRESENT(f3din4)) THEN
    p_in(4)%fld  => f3din4
    p_out(4)%fld => f3dout4
  ENDIF
  IF (PRESENT(f3din5)) THEN
    p_in(5)%fld  => f3din5
    p_out(5)%fld => f3dout5
  ENDIF
  l4d = .FALSE.
ENDIF

IF (PRESENT(llimit_nneg)) THEN
  l_limit_nneg(:) = llimit_nneg(:)
ELSE
  l_limit_nneg(:) = .FALSE.
ENDIF

IF (PRESENT(rlimval)) THEN
  r_limval(:) = rlimval(:)
ELSE
  r_limval(:) = 0._wp
ENDIF

! factor of allowed overshooting
IF (PRESENT(overshoot_fac)) THEN
  ovsht_fac = overshoot_fac
ELSE
  ovsht_fac = 1.05_wp
ENDIF

    elev = 1
    DO jn = 1, nfields
      elev = MAX(elev, UBOUND(p_out(jn)%fld,2))
    END DO
    ALLOCATE(l_enabled(elev))
    IF (PRESENT(opt_l_enabled)) THEN
      l_enabled = opt_l_enabled(1:elev)
    ELSE
      l_enabled = .TRUE.
    END IF
    IF (ALL(l_enabled(:))) THEN
      all_enabled = .TRUE.
    ELSE
      all_enabled = .FALSE.
    ENDIF

r_ovsht_fac = 1._wp/ovsht_fac

! Start and end blocks for which scalar interpolation is needed
i_startblk = ptr_pp%cells%start_block(grf_nudgintp_start_c)
i_endblk   = ptr_pp%cells%end_block(min_rlcell_int)

! Pointers to index lists and interpolation coefficients for gradient computation
iidx => ptr_int%rbf_c2grad_idx
iblk => ptr_int%rbf_c2grad_blk

ptr_coeff => ptr_int%rbf_c2grad_coeff
ptr_dist  => ptr_grf%grf_dist_pc2cc

! child cell indices and blocks for non-MPI parent-to-child communication
ichcidx => ptr_pp%cells%child_idx
ichcblk => ptr_pp%cells%child_blk

! Shift parameter
js = nshift

!$OMP PARALLEL PRIVATE(jn,elev)

DO jn = 1, nfields
  elev = UBOUND(p_out(jn)%fld,2)

!$OMP DO PRIVATE (jb,jk,jc,i_startidx,i_endidx,grad_x,grad_y,min_expval, &
!$OMP   max_expval,limfac1,limfac2,limfac,maxval_neighb,minval_neighb,   &
!$OMP   relaxed_minval,relaxed_maxval) ICON_OMP_DEFAULT_SCHEDULE
  DO jb = i_startblk, i_endblk

    CALL get_indices_c(ptr_pp, jb, i_startblk, i_endblk, &
         i_startidx, i_endidx, grf_nudgintp_start_c, min_rlcell_int)

    IF (all_enabled) THEN ! Use vectorizable form with loop reordering
#ifdef __LOOP_EXCHANGE
      DO jc = i_startidx, i_endidx
        DO jk = 1, elev
#else
      DO jk = 1, elev
        DO jc = i_startidx, i_endidx
#endif
          grad_x(jc,jk) =  &
            ptr_coeff(1,1,jc,jb)*p_in(jn)%fld(jc,jk+js,jb) + &
            ptr_coeff(2,1,jc,jb)*p_in(jn)%fld(iidx(2,jc,jb),jk+js,iblk(2,jc,jb)) + &
            ptr_coeff(3,1,jc,jb)*p_in(jn)%fld(iidx(3,jc,jb),jk+js,iblk(3,jc,jb)) + &
            ptr_coeff(4,1,jc,jb)*p_in(jn)%fld(iidx(4,jc,jb),jk+js,iblk(4,jc,jb)) + &
            ptr_coeff(5,1,jc,jb)*p_in(jn)%fld(iidx(5,jc,jb),jk+js,iblk(5,jc,jb)) + &
            ptr_coeff(6,1,jc,jb)*p_in(jn)%fld(iidx(6,jc,jb),jk+js,iblk(6,jc,jb)) + &
            ptr_coeff(7,1,jc,jb)*p_in(jn)%fld(iidx(7,jc,jb),jk+js,iblk(7,jc,jb)) + &
            ptr_coeff(8,1,jc,jb)*p_in(jn)%fld(iidx(8,jc,jb),jk+js,iblk(8,jc,jb)) + &
            ptr_coeff(9,1,jc,jb)*p_in(jn)%fld(iidx(9,jc,jb),jk+js,iblk(9,jc,jb)) + &
            ptr_coeff(10,1,jc,jb)*p_in(jn)%fld(iidx(10,jc,jb),jk+js,iblk(10,jc,jb))
          grad_y(jc,jk) =  &
            ptr_coeff(1,2,jc,jb)*p_in(jn)%fld(jc,jk+js,jb) + &
            ptr_coeff(2,2,jc,jb)*p_in(jn)%fld(iidx(2,jc,jb),jk+js,iblk(2,jc,jb)) + &
            ptr_coeff(3,2,jc,jb)*p_in(jn)%fld(iidx(3,jc,jb),jk+js,iblk(3,jc,jb)) + &
            ptr_coeff(4,2,jc,jb)*p_in(jn)%fld(iidx(4,jc,jb),jk+js,iblk(4,jc,jb)) + &
            ptr_coeff(5,2,jc,jb)*p_in(jn)%fld(iidx(5,jc,jb),jk+js,iblk(5,jc,jb)) + &
            ptr_coeff(6,2,jc,jb)*p_in(jn)%fld(iidx(6,jc,jb),jk+js,iblk(6,jc,jb)) + &
            ptr_coeff(7,2,jc,jb)*p_in(jn)%fld(iidx(7,jc,jb),jk+js,iblk(7,jc,jb)) + &
            ptr_coeff(8,2,jc,jb)*p_in(jn)%fld(iidx(8,jc,jb),jk+js,iblk(8,jc,jb)) + &
            ptr_coeff(9,2,jc,jb)*p_in(jn)%fld(iidx(9,jc,jb),jk+js,iblk(9,jc,jb)) + &
            ptr_coeff(10,2,jc,jb)*p_in(jn)%fld(iidx(10,jc,jb),jk+js,iblk(10,jc,jb))
          maxval_neighb(jc,jk) =                                 &
            MAX(p_in(jn)%fld(jc,jk+js,jb),                       &
                p_in(jn)%fld(iidx(2,jc,jb),jk+js,iblk(2,jc,jb)), &
                p_in(jn)%fld(iidx(3,jc,jb),jk+js,iblk(3,jc,jb)), &
                p_in(jn)%fld(iidx(4,jc,jb),jk+js,iblk(4,jc,jb)), &
                p_in(jn)%fld(iidx(5,jc,jb),jk+js,iblk(5,jc,jb)), &
                p_in(jn)%fld(iidx(6,jc,jb),jk+js,iblk(6,jc,jb)), &
                p_in(jn)%fld(iidx(7,jc,jb),jk+js,iblk(7,jc,jb)), &
                p_in(jn)%fld(iidx(8,jc,jb),jk+js,iblk(8,jc,jb)), &
                p_in(jn)%fld(iidx(9,jc,jb),jk+js,iblk(9,jc,jb)), &
                p_in(jn)%fld(iidx(10,jc,jb),jk+js,iblk(10,jc,jb)))
          minval_neighb(jc,jk) =                                 &
            MIN(p_in(jn)%fld(jc,jk+js,jb),                       &
                p_in(jn)%fld(iidx(2,jc,jb),jk+js,iblk(2,jc,jb)), &
                p_in(jn)%fld(iidx(3,jc,jb),jk+js,iblk(3,jc,jb)), &
                p_in(jn)%fld(iidx(4,jc,jb),jk+js,iblk(4,jc,jb)), &
                p_in(jn)%fld(iidx(5,jc,jb),jk+js,iblk(5,jc,jb)), &
                p_in(jn)%fld(iidx(6,jc,jb),jk+js,iblk(6,jc,jb)), &
                p_in(jn)%fld(iidx(7,jc,jb),jk+js,iblk(7,jc,jb)), &
                p_in(jn)%fld(iidx(8,jc,jb),jk+js,iblk(8,jc,jb)), &
                p_in(jn)%fld(iidx(9,jc,jb),jk+js,iblk(9,jc,jb)), &
                p_in(jn)%fld(iidx(10,jc,jb),jk+js,iblk(10,jc,jb)))
        ENDDO
      ENDDO
    ELSE
      DO jk = 1, elev
        IF (.NOT. l_enabled(jk)) CYCLE
        DO jc = i_startidx, i_endidx
          grad_x(jc,jk) =  &
            ptr_coeff(1,1,jc,jb)*p_in(jn)%fld(jc,jk+js,jb) + &
            ptr_coeff(2,1,jc,jb)*p_in(jn)%fld(iidx(2,jc,jb),jk+js,iblk(2,jc,jb)) + &
            ptr_coeff(3,1,jc,jb)*p_in(jn)%fld(iidx(3,jc,jb),jk+js,iblk(3,jc,jb)) + &
            ptr_coeff(4,1,jc,jb)*p_in(jn)%fld(iidx(4,jc,jb),jk+js,iblk(4,jc,jb)) + &
            ptr_coeff(5,1,jc,jb)*p_in(jn)%fld(iidx(5,jc,jb),jk+js,iblk(5,jc,jb)) + &
            ptr_coeff(6,1,jc,jb)*p_in(jn)%fld(iidx(6,jc,jb),jk+js,iblk(6,jc,jb)) + &
            ptr_coeff(7,1,jc,jb)*p_in(jn)%fld(iidx(7,jc,jb),jk+js,iblk(7,jc,jb)) + &
            ptr_coeff(8,1,jc,jb)*p_in(jn)%fld(iidx(8,jc,jb),jk+js,iblk(8,jc,jb)) + &
            ptr_coeff(9,1,jc,jb)*p_in(jn)%fld(iidx(9,jc,jb),jk+js,iblk(9,jc,jb)) + &
            ptr_coeff(10,1,jc,jb)*p_in(jn)%fld(iidx(10,jc,jb),jk+js,iblk(10,jc,jb))
          grad_y(jc,jk) =  &
            ptr_coeff(1,2,jc,jb)*p_in(jn)%fld(jc,jk+js,jb) + &
            ptr_coeff(2,2,jc,jb)*p_in(jn)%fld(iidx(2,jc,jb),jk+js,iblk(2,jc,jb)) + &
            ptr_coeff(3,2,jc,jb)*p_in(jn)%fld(iidx(3,jc,jb),jk+js,iblk(3,jc,jb)) + &
            ptr_coeff(4,2,jc,jb)*p_in(jn)%fld(iidx(4,jc,jb),jk+js,iblk(4,jc,jb)) + &
            ptr_coeff(5,2,jc,jb)*p_in(jn)%fld(iidx(5,jc,jb),jk+js,iblk(5,jc,jb)) + &
            ptr_coeff(6,2,jc,jb)*p_in(jn)%fld(iidx(6,jc,jb),jk+js,iblk(6,jc,jb)) + &
            ptr_coeff(7,2,jc,jb)*p_in(jn)%fld(iidx(7,jc,jb),jk+js,iblk(7,jc,jb)) + &
            ptr_coeff(8,2,jc,jb)*p_in(jn)%fld(iidx(8,jc,jb),jk+js,iblk(8,jc,jb)) + &
            ptr_coeff(9,2,jc,jb)*p_in(jn)%fld(iidx(9,jc,jb),jk+js,iblk(9,jc,jb)) + &
            ptr_coeff(10,2,jc,jb)*p_in(jn)%fld(iidx(10,jc,jb),jk+js,iblk(10,jc,jb))
          maxval_neighb(jc,jk) =                                 &
            MAX(p_in(jn)%fld(jc,jk+js,jb),                       &
                p_in(jn)%fld(iidx(2,jc,jb),jk+js,iblk(2,jc,jb)), &
                p_in(jn)%fld(iidx(3,jc,jb),jk+js,iblk(3,jc,jb)), &
                p_in(jn)%fld(iidx(4,jc,jb),jk+js,iblk(4,jc,jb)), &
                p_in(jn)%fld(iidx(5,jc,jb),jk+js,iblk(5,jc,jb)), &
                p_in(jn)%fld(iidx(6,jc,jb),jk+js,iblk(6,jc,jb)), &
                p_in(jn)%fld(iidx(7,jc,jb),jk+js,iblk(7,jc,jb)), &
                p_in(jn)%fld(iidx(8,jc,jb),jk+js,iblk(8,jc,jb)), &
                p_in(jn)%fld(iidx(9,jc,jb),jk+js,iblk(9,jc,jb)), &
                p_in(jn)%fld(iidx(10,jc,jb),jk+js,iblk(10,jc,jb)))
          minval_neighb(jc,jk) =                                 &
            MIN(p_in(jn)%fld(jc,jk+js,jb),                       &
                p_in(jn)%fld(iidx(2,jc,jb),jk+js,iblk(2,jc,jb)), &
                p_in(jn)%fld(iidx(3,jc,jb),jk+js,iblk(3,jc,jb)), &
                p_in(jn)%fld(iidx(4,jc,jb),jk+js,iblk(4,jc,jb)), &
                p_in(jn)%fld(iidx(5,jc,jb),jk+js,iblk(5,jc,jb)), &
                p_in(jn)%fld(iidx(6,jc,jb),jk+js,iblk(6,jc,jb)), &
                p_in(jn)%fld(iidx(7,jc,jb),jk+js,iblk(7,jc,jb)), &
                p_in(jn)%fld(iidx(8,jc,jb),jk+js,iblk(8,jc,jb)), &
                p_in(jn)%fld(iidx(9,jc,jb),jk+js,iblk(9,jc,jb)), &
                p_in(jn)%fld(iidx(10,jc,jb),jk+js,iblk(10,jc,jb)))
        ENDDO
      ENDDO
    ENDIF

    DO jk = 1, elev
      IF (.NOT. l_enabled(jk)) CYCLE
      DO jc = i_startidx, i_endidx
        min_expval(jc) = MIN(grad_x(jc,jk)*ptr_dist(jc,1,1,jb) + &
                           grad_y(jc,jk)*ptr_dist(jc,1,2,jb),  &
                           grad_x(jc,jk)*ptr_dist(jc,2,1,jb) + &
                           grad_y(jc,jk)*ptr_dist(jc,2,2,jb),  &
                           grad_x(jc,jk)*ptr_dist(jc,3,1,jb) + &
                           grad_y(jc,jk)*ptr_dist(jc,3,2,jb),  &
                           grad_x(jc,jk)*ptr_dist(jc,4,1,jb) + &
                           grad_y(jc,jk)*ptr_dist(jc,4,2,jb),  &
                           -1.e-80_wp )
        max_expval(jc) = MAX(grad_x(jc,jk)*ptr_dist(jc,1,1,jb) + &
                           grad_y(jc,jk)*ptr_dist(jc,1,2,jb),  &
                           grad_x(jc,jk)*ptr_dist(jc,2,1,jb) + &
                           grad_y(jc,jk)*ptr_dist(jc,2,2,jb),  &
                           grad_x(jc,jk)*ptr_dist(jc,3,1,jb) + &
                           grad_y(jc,jk)*ptr_dist(jc,3,2,jb),  &
                           grad_x(jc,jk)*ptr_dist(jc,4,1,jb) + &
                           grad_y(jc,jk)*ptr_dist(jc,4,2,jb),  &
                           1.e-80_wp )
      ENDDO

      DO jc = i_startidx, i_endidx

        ! Allow a limited amount of over-/undershooting in the downscaled fields
        relaxed_minval = MERGE(r_ovsht_fac, ovsht_fac, &
             minval_neighb(jc,jk) > 0._wp) * minval_neighb(jc,jk)
        relaxed_maxval = MERGE(ovsht_fac, r_ovsht_fac, &
             maxval_neighb(jc,jk) > 0._wp) * maxval_neighb(jc,jk)

        limfac1 = MERGE(1._wp, &
             ABS((relaxed_minval-p_in(jn)%fld(jc,jk+js,jb))/min_expval(jc)), &
             p_in(jn)%fld(jc,jk+js,jb) + min_expval(jc) >= relaxed_minval-epsi)
        limfac2 = MERGE(1._wp, &
             ABS((relaxed_maxval-p_in(jn)%fld(jc,jk+js,jb))/max_expval(jc)), &
             p_in(jn)%fld(jc,jk+js,jb) + max_expval(jc) <= relaxed_maxval+epsi)
        limfac = MIN(limfac1,limfac2)

        grad_x(jc,jk) = grad_x(jc,jk)*limfac
        grad_y(jc,jk) = grad_y(jc,jk)*limfac

      ENDDO
    ENDDO

    DO jk = 1, elev
      IF (.NOT. l_enabled(jk)) CYCLE
      DO jc = i_startidx, i_endidx

        h_aux(jc,jk,jb,1,jn) = p_in(jn)%fld(jc,jk+js,jb) + &
          grad_x(jc,jk)*ptr_dist(jc,1,1,jb)              + &
          grad_y(jc,jk)*ptr_dist(jc,1,2,jb)
        h_aux(jc,jk,jb,2,jn) = p_in(jn)%fld(jc,jk+js,jb) + &
          grad_x(jc,jk)*ptr_dist(jc,2,1,jb)              + &
          grad_y(jc,jk)*ptr_dist(jc,2,2,jb)
        h_aux(jc,jk,jb,3,jn) = p_in(jn)%fld(jc,jk+js,jb) + &
          grad_x(jc,jk)*ptr_dist(jc,3,1,jb)              + &
          grad_y(jc,jk)*ptr_dist(jc,3,2,jb)
        h_aux(jc,jk,jb,4,jn) = p_in(jn)%fld(jc,jk+js,jb) + &
          grad_x(jc,jk)*ptr_dist(jc,4,1,jb)              + &
          grad_y(jc,jk)*ptr_dist(jc,4,2,jb)

      ENDDO
    ENDDO
  ENDDO ! blocks

!$OMP END DO
ENDDO ! fields

! Store results in p_out

! For MPI parallelization, the child cell indices have to be mapped to the
! feedback-parent grid

DO jn = 1, nfields

  elev = UBOUND(p_out(jn)%fld,2)

!$OMP DO PRIVATE (jb,jk,jc,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
  DO jb =  i_startblk, i_endblk

    CALL get_indices_c(ptr_pp, jb, i_startblk, i_endblk, &
         i_startidx, i_endidx, grf_nudgintp_start_c, min_rlcell_int)

    IF (l_limit_nneg(jn)) THEN

#ifdef __LOOP_EXCHANGE
      DO jc = i_startidx, i_endidx
        DO jk = 1, elev
          IF (.NOT. l_enabled(jk)) CYCLE
#else
      DO jk = 1, elev
        IF (.NOT. l_enabled(jk)) CYCLE
!CDIR NODEP,VOVERTAKE,VOB
        DO jc = i_startidx, i_endidx
#endif
          p_out(jn)%fld(ichcidx(jc,jb,1),jk,ichcblk(jc,jb,1)) = &
            MAX(h_aux(jc,jk,jb,1,jn),r_limval(jn))
          p_out(jn)%fld(ichcidx(jc,jb,2),jk,ichcblk(jc,jb,2)) = &
            MAX(h_aux(jc,jk,jb,2,jn),r_limval(jn))
          p_out(jn)%fld(ichcidx(jc,jb,3),jk,ichcblk(jc,jb,3)) = &
            MAX(h_aux(jc,jk,jb,3,jn),r_limval(jn))
          p_out(jn)%fld(ichcidx(jc,jb,4),jk,ichcblk(jc,jb,4)) = &
            MAX(h_aux(jc,jk,jb,4,jn),r_limval(jn))

        ENDDO
      ENDDO

    ELSE

#ifdef __LOOP_EXCHANGE
      DO jc = i_startidx, i_endidx
        DO jk = 1, elev
          IF (.NOT. l_enabled(jk)) CYCLE
#else
      DO jk = 1, elev
        IF (.NOT. l_enabled(jk)) CYCLE
!CDIR NODEP,VOVERTAKE,VOB
        DO jc = i_startidx, i_endidx
#endif

          p_out(jn)%fld(ichcidx(jc,jb,1),jk,ichcblk(jc,jb,1)) = h_aux(jc,jk,jb,1,jn)
          p_out(jn)%fld(ichcidx(jc,jb,2),jk,ichcblk(jc,jb,2)) = h_aux(jc,jk,jb,2,jn)
          p_out(jn)%fld(ichcidx(jc,jb,3),jk,ichcblk(jc,jb,3)) = h_aux(jc,jk,jb,3,jn)
          p_out(jn)%fld(ichcidx(jc,jb,4),jk,ichcblk(jc,jb,4)) = h_aux(jc,jk,jb,4,jn)

        ENDDO
      ENDDO

    ENDIF
  ENDDO
!$OMP END DO NOWAIT
ENDDO

!$OMP END PARALLEL

END SUBROUTINE interpol_scal_nudging

!-------------------------------------------------------------------------
END MODULE mo_grf_nudgintp
