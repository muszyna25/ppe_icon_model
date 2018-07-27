!>
!!  This module contains the routines needed for managing flow control.
!!
!!  This module contains the routines needed for managing flow control
!!  with mesh refinement. The main routine, process_level, is a recursive
!!  subroutine that is called from sw_atmos for the global mesh and calls
!!  itself recursively for the refined meshes. It contains the whole time
!!  stepping management that was previously located in sw_atmos.
!!  Further subroutines serve for interpolating the time tendencies to the
!!  lateral boundaries of the refined meshes, for smoothing these fields,
!!  and for feedback.
!!
!! @par Revision History
!!  Developed and tested by Guenther Zaengl, DWD (2008-07)
!!  Modified by Marco Giorgetta, MPI-M (2009-02-26)
!!  - renamed ltracer to ltransport
!!  Modification by Guenther Zaengl, DWD (2009-06-22)
!!  - preparation for generalized grid refinement (affects all subroutines)
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
MODULE mo_hierarchy_management_intp
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
USE mo_exception,           ONLY: message_text, message
USE mo_model_domain,        ONLY: t_patch, t_grid_cells, t_grid_edges,    &
                                  p_patch_local_parent
USE mo_grid_config,         ONLY: n_dom
USE mo_intp_data_strc,      ONLY: t_int_state
USE mo_intp_rbf,            ONLY: rbf_vec_interpol_vertex
USE mo_grf_intp_data_strc,  ONLY: t_gridref_state
USE mo_gridref_config,      ONLY: grf_velfbk, grf_intmethod_c, grf_intmethod_e,   &
                                  grf_intmethod_ct, grf_scalfbk, grf_tracfbk
USE mo_grf_intp_data_strc , ONLY: p_grf_state_local_parent
USE mo_grf_bdyintp,         ONLY: interpol_scal_grf, interpol_vec_grf, interpol2_vec_grf
USE mo_dynamics_config,     ONLY: nnew, nnow, nsav1, nsav2, lshallow_water
USE mo_parallel_config,  ONLY: nproma, p_test_run
USE mo_run_config,          ONLY: msg_level, ltransport, nlev,    &
                                  ntracer
USE mo_icoham_dyn_types,    ONLY: t_hydro_atm, t_hydro_atm_prog,  &
                                  t_hydro_atm_diag
USE mo_impl_constants,      ONLY: min_rlcell_int, min_rledge, min_rledge_int, &
    &                             MAX_CHAR_LENGTH
USE mo_loopindices,         ONLY: get_indices_c, get_indices_e
USE mo_impl_constants_grf,  ONLY: grf_bdyintp_start_c, &
                                  grf_bdyintp_end_c, &
                                  grf_fbk_start_c, grf_fbk_start_e,        &
                                  grf_bdywidth_c, grf_bdywidth_e
USE mo_communication,       ONLY: exchange_data
USE mo_sync,                ONLY: SYNC_C, SYNC_E, sync_patch_array, check_patch_array, &
                                & global_sum_array2


IMPLICIT NONE

PRIVATE

PUBLIC :: interpolate_diagnostics, interpolate_tendencies, boundary_tendencies, feedback

CONTAINS



!-------------------------------------------------------------------------
!
!
!
!>
!! Interpolates time tendencies of prognostic variables to the lateral boundary.
!!
!! Interpolates time tendencies of prognostic variables to the lateral boundary
!! of a refined mesh
!!
!! @par Revision History
!! Developed  by Guenther Zaengl, DWD, 2008-07-10
!!
SUBROUTINE interpolate_tendencies (p_patch,p_hydro_state,p_int_state,p_grf_state,jg,jgc)

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = 'mo_hierarchy_management_intp:interpolate_tendencies'

TYPE(t_patch),       TARGET, INTENT(IN)    ::  p_patch(n_dom)
TYPE(t_hydro_atm), TARGET, INTENT(INOUT) ::  p_hydro_state(n_dom)
TYPE(t_int_state),     TARGET, INTENT(IN)  ::  p_int_state(n_dom)
TYPE(t_gridref_state), TARGET, INTENT(IN)  ::  p_grf_state(n_dom)

INTEGER, INTENT(IN)     :: jg, jgc      ! domain ID of parent and child grid

! local variables

TYPE(t_hydro_atm_prog), POINTER    :: p_parent_tend => NULL()
TYPE(t_hydro_atm_prog), POINTER    :: p_child_tend => NULL()
TYPE(t_patch), POINTER             :: p_pp => NULL()
TYPE(t_patch), POINTER             :: p_pc => NULL()
TYPE(t_int_state), POINTER         :: p_int => NULL()
TYPE(t_gridref_state), POINTER     :: p_grf => NULL()
TYPE(t_grid_cells), POINTER        :: p_gcp => NULL()
TYPE(t_grid_cells), POINTER        :: p_gcc => NULL()
TYPE(t_grid_edges), POINTER        :: p_gep => NULL()
TYPE(t_grid_edges), POINTER        :: p_gec => NULL()


INTEGER :: jt        ! loop indices

INTEGER :: i_chidx

! Pointers to index fields
INTEGER,  DIMENSION(:,:,:),   POINTER :: iidx, iblk

REAL(wp) :: z_pres_sfc(nproma,1,p_patch(jgc)%nblks_c)

!-----------------------------------------------------------------------

IF (msg_level >= 10) THEN
  WRITE(message_text,'(a,i2,a,i2)') '========= Interpolate:',jg,' =>',jgc
  CALL message(TRIM(routine),message_text)
ENDIF

p_parent_tend => p_hydro_state(jg)%tend_dyn
p_child_tend  => p_hydro_state(jgc)%tend_dyn
p_int         => p_int_state(jg)
p_grf         => p_grf_state(jg)
p_pp          => p_patch(jg)
p_pc          => p_patch(jgc)
p_gcp         => p_patch(jg)%cells
p_gcc         => p_patch(jgc)%cells
p_gep         => p_patch(jg)%edges
p_gec         => p_patch(jgc)%edges

iidx          => p_gcp%child_idx
iblk          => p_gcp%child_blk

i_chidx = p_patch(jgc)%parent_child_index

! Interpolation of cell-based variables

IF (grf_intmethod_c == 1) THEN ! tendency copying for pressure and temperature

    CALL exchange_data(p_pc%comm_pat_interpolation_c,      &
                       RECV=p_child_tend%pres_sfc, &
                       SEND=p_parent_tend%pres_sfc)
    CALL exchange_data(p_pc%comm_pat_interpolation_c,      &
                       RECV=p_child_tend%temp,     &
                       SEND=p_parent_tend%temp)


! grf_intmethod_c = 2, use gradient at cell center for interpolation
ELSE IF (grf_intmethod_c == 2) THEN

  CALL interpol_scal_grf (p_pp, p_pc, p_grf%p_dom(i_chidx), 1,&
                          RESHAPE(p_parent_tend%pres_sfc,(/nproma,1,p_pp%nblks_c/)), z_pres_sfc)

  p_child_tend%pres_sfc(:,:) = z_pres_sfc(:,1,:)

  IF (.NOT. lshallow_water) &
    CALL interpol_scal_grf (p_pp, p_pc, p_grf%p_dom(i_chidx), 1, &
                            p_parent_tend%temp,  p_child_tend%temp)

ENDIF

IF (ltransport .AND. grf_intmethod_ct == 1) THEN

  DO jt = 1, ntracer
    CALL exchange_data(p_pc%comm_pat_interpolation_c, &
                       RECV=p_child_tend%tracer(:,:,:,jt), &
                       SEND=p_parent_tend%tracer(:,:,:,jt))
  ENDDO

IF (p_test_run) CALL check_patch_array(0,p_pc,p_child_tend%tracer,'IT:tracer')

ELSE IF (ltransport .AND. grf_intmethod_ct == 2) THEN

  CALL interpol_scal_grf (p_pp, p_pc, p_grf%p_dom(i_chidx), ntracer,             &
                          f4din1=p_parent_tend%tracer, f4dout1=p_child_tend%tracer)

IF (p_test_run) CALL check_patch_array(0,p_pc,p_child_tend%tracer,'IT:tracer')
ENDIF

! Interpolation of edge-based variables  (velocity components)
IF ((grf_intmethod_e == 1) .OR. (grf_intmethod_e == 2)) THEN

  CALL interpol_vec_grf (p_pp, p_pc, p_grf%p_dom(i_chidx), &
                         p_parent_tend%vn, p_child_tend%vn)

ELSE IF ((grf_intmethod_e == 3) .OR. (grf_intmethod_e == 4)) THEN

  CALL interpol2_vec_grf (p_pp, p_pc, p_grf%p_dom(i_chidx), 1, &
                          p_parent_tend%vn, p_child_tend%vn)

ENDIF

END SUBROUTINE interpolate_tendencies



!-------------------------------------------------------------------------
!
!
!
!>
!! Interpolates diagnostic variables to the lateral boundary of a refined mesh.
!!
!! Interpolates diagnostic variables to the lateral boundary of a refined mesh
!! before writing output.
!!
!! @par Revision History
!! Developed  by Guenther Zaengl, DWD, 2009-02-09
!!
SUBROUTINE interpolate_diagnostics (p_patch, p_diag, p_diagc, &
                                    p_int_state,p_grf_state,jg,jgc)



TYPE(t_patch),       TARGET, INTENT(IN)    :: p_patch(n_dom)
TYPE(t_hydro_atm_diag),  TARGET, INTENT(INOUT) :: p_diag
TYPE(t_hydro_atm_diag),  TARGET, INTENT(INOUT) :: p_diagc
TYPE(t_int_state),   TARGET, INTENT(IN)    :: p_int_state(n_dom)
TYPE(t_gridref_state), TARGET, INTENT(IN)  :: p_grf_state(n_dom)

INTEGER, INTENT(IN)     :: jg,jgc

! local variables

TYPE(t_patch), POINTER             :: p_pp => NULL()
TYPE(t_patch), POINTER             :: p_pc => NULL()
TYPE(t_int_state), POINTER         :: p_int => NULL()
TYPE(t_gridref_state), POINTER     :: p_grf => NULL()

INTEGER :: i_chidx

!-----------------------------------------------------------------------

p_int         => p_int_state(jg)
p_grf         => p_grf_state(jg)
p_pp          => p_patch(jg)
p_pc          => p_patch(jgc)

i_chidx = p_patch(jgc)%parent_child_index

CALL sync_patch_array(SYNC_C, p_pp, p_diag%u)
CALL sync_patch_array(SYNC_C, p_pp, p_diag%v)
CALL sync_patch_array(SYNC_C, p_pp, p_diag%wpres_mc)

CALL interpol_scal_grf (p_pp, p_pc, p_grf%p_dom(i_chidx), 3,      &
                        p_diag%u, p_diagc%u, p_diag%v, p_diagc%v, &
                        p_diag%wpres_mc, p_diagc%wpres_mc)

CALL sync_patch_array(SYNC_C, p_pp, p_diag%div)

CALL interpol_scal_grf (p_pp, p_pc, p_grf%p_dom(i_chidx), 1, &
                        p_diag%div, p_diagc%div)

END SUBROUTINE interpolate_diagnostics




!-------------------------------------------------------------------------
!
!
!
!>
!! Dummy routine needed to test the technical implementation of the limited-area.
!!
!! Dummy routine needed to test the technical implementation of the limited-area
!! mode. Here, the boundary tendencies of the prognostic variables are simply
!! set to zero, implying fixed lateral boundary conditions.
!! This routine will have to be replaced by a real interpolation of lateral
!! boundary tendencies later on.
!!
!! @par Revision History
!! Developed  by Guenther Zaengl, DWD, 2009-06-22
!!
SUBROUTINE boundary_tendencies (p_patch,p_hydro_state)



TYPE(t_patch),       TARGET, INTENT(IN)    ::  p_patch
TYPE(t_hydro_atm), TARGET, INTENT(INOUT) ::  p_hydro_state

! local variables

TYPE(t_hydro_atm_prog), POINTER    :: p_tend => NULL()

INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx,       &
           jb, jc, je, jk, jt
!-----------------------------------------------------------------------

p_tend => p_hydro_state%tend_dyn

! Cell-based variables

! Start and end blocks for which interpolation is needed
i_startblk = p_patch%cells%start_blk(1,1)
i_endblk   = p_patch%cells%end_blk(grf_bdywidth_c,1)

DO jb =  i_startblk, i_endblk

  CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                     i_startidx, i_endidx, 1, grf_bdywidth_c)

  DO jc = i_startidx, i_endidx
    p_tend%pres_sfc(jc,jb) = 0._wp
  ENDDO

  DO jk = 1, nlev
    DO jc = i_startidx, i_endidx
      p_tend%temp(jc,jk,jb) = 0._wp
    ENDDO
  ENDDO

    DO jt = 1, ntracer
      DO jk = 1, nlev
        DO jc = i_startidx, i_endidx
          p_tend%tracer(jc,jk,jb,jt) = 0._wp
        ENDDO
      ENDDO
    ENDDO

ENDDO

! Edge-based variables

! Start and end blocks for which vector interpolation is needed
i_startblk = p_patch%edges%start_blk(1,1)
i_endblk   = p_patch%edges%end_blk(grf_bdywidth_e,1)

DO jb =  i_startblk, i_endblk

  CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                     i_startidx, i_endidx, 1, grf_bdywidth_e)

  DO jk = 1, nlev
    DO je = i_startidx, i_endidx
      p_tend%vn(je,jk,jb) = 0._wp
    ENDDO
  ENDDO
ENDDO

END SUBROUTINE boundary_tendencies



!-------------------------------------------------------------------------
!
!
!
!>
!! This routine computes the feedback of the prognostic variables from the fine mesh.
!!
!! This routine computes the feedback of the prognostic variables from the fine mesh
!! to to the corresponding grid point on the coarse mesh
!! jg in this case denotes the fine mesh level; output goes to parent_id(jg)
!!
!! @par Revision History
!! Developed  by Guenther Zaengl, DWD, 2008-04-15
!! @par
!! Modification by Guenther Zaengl, DWD, 2008-09-12:
!! Change feedback for cell-based variables from area-weighted averaging
!! to using fbk_wgt (see above routine)
!!
SUBROUTINE feedback(ltheta_dyn, p_patch, p_hydro_state, p_int_state, p_grf_state, jg, jgp)

CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = 'mo_hierarchy_management_intp:feedback'

LOGICAL, INTENT(IN) :: ltheta_dyn

TYPE(t_patch),        TARGET, INTENT(INOUT) ::  p_patch(n_dom)
TYPE(t_hydro_atm),    TARGET, INTENT(INOUT) ::  p_hydro_state(n_dom)
TYPE(t_int_state),    TARGET, INTENT(IN)    ::  p_int_state(n_dom)
TYPE(t_gridref_state),TARGET, INTENT(IN)    ::  p_grf_state(n_dom)

INTEGER, INTENT(IN) :: jg   ! child grid level
INTEGER, INTENT(IN) :: jgp  ! parent grid level

! local variables

TYPE(t_hydro_atm_prog),     POINTER :: p_parent_prog => NULL()
TYPE(t_hydro_atm_prog),     POINTER :: p_parent_save => NULL()
TYPE(t_hydro_atm_prog), POINTER :: p_parent_tend => NULL()
TYPE(t_hydro_atm_prog),     POINTER :: p_child_prog => NULL()
TYPE(t_hydro_atm_prog),     POINTER :: p_child_save => NULL()
TYPE(t_hydro_atm_prog), POINTER :: p_child_tend => NULL()
TYPE(t_grid_cells), POINTER     :: p_gcp => NULL()
TYPE(t_grid_cells), POINTER     :: p_gcc => NULL()
TYPE(t_grid_edges), POINTER     :: p_gep => NULL()
TYPE(t_grid_edges), POINTER     :: p_gec => NULL()
TYPE(t_gridref_state), POINTER  :: p_grf => NULL()
TYPE(t_int_state), POINTER      :: p_intc => NULL()
TYPE(t_patch),      POINTER     :: p_pp => NULL()
TYPE(t_patch),      POINTER     :: p_pc => NULL()

REAL(wp), DIMENSION(:,:,:), POINTER :: p_temp_prog => NULL()
REAL(wp), DIMENSION(:,:,:), POINTER :: p_temp_save => NULL()

! Indices
INTEGER :: jb, jc, jk, jt, je, i_nchdom, i_chidx,               &
           i_startblk, i_endblk, i_startidx, i_endidx


REAL(wp), DIMENSION(nproma,nlev,p_patch(jg)%nblks_v) :: z_u, z_v
REAL(wp) :: vn_aux(nproma,nlev,p_patch(jg)%nblks_e,2) ! RBF-reconstructed velocity
REAL(wp), ALLOCATABLE :: feedback_pres_tend(:,:)
REAL(wp), ALLOCATABLE :: feedback_temp_tend(:,:,:)
REAL(wp), ALLOCATABLE :: feedback_vn_tend(:,:,:)
REAL(wp), ALLOCATABLE :: feedback_tracer_tend(:,:,:,:)
REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: fbk_tend, parent_tend
REAL(wp) :: tendency_corr



INTEGER, DIMENSION(:,:,:), POINTER :: iidx, iblk, iidxv, iblkv
REAL(wp), DIMENSION(:,:,:), POINTER :: p_fbkwgt, p_fbkwgt_tr
REAL(wp), DIMENSION(:,:),   POINTER :: p_fbarea

!-----------------------------------------------------------------------

IF (msg_level >= 10) THEN
  WRITE(message_text,'(a,i2,a,i2)') '========= Feedback:',jg,' =>',jgp
  CALL message(TRIM(routine),message_text)
ENDIF

p_parent_prog => p_hydro_state(jgp)%prog(nnew(jgp))
p_parent_save => p_hydro_state(jgp)%prog(nsav1(jgp))
p_parent_tend => p_hydro_state(jgp)%tend_dyn
p_child_prog  => p_hydro_state(jg)%prog(nnow(jg))
p_child_save  => p_hydro_state(jg)%prog(nsav2(jg))
p_child_tend  => p_hydro_state(jg)%tend_dyn
p_intc        => p_int_state(jg)
p_gcc         => p_patch(jg)%cells
p_gec         => p_patch(jg)%edges
p_pc          => p_patch(jg)

p_grf => p_grf_state_local_parent(jg)
p_gcp => p_patch_local_parent(jg)%cells
p_gep => p_patch_local_parent(jg)%edges
p_pp  => p_patch_local_parent(jg)

i_nchdom = MAX(1,p_pc%n_childdom)
i_chidx  = p_pc%parent_child_index

! parent_tend is always calculated on the global parent
! and thus has to be allocated within global parent limits

i_startblk = p_patch(jgp)%cells%start_blk(grf_fbk_start_c,i_chidx)
i_endblk   = p_patch(jgp)%cells%end_blk(min_rlcell_int,i_chidx)

ALLOCATE(parent_tend(nproma, i_startblk:i_endblk))

! Allocation of storage fields
! In parallel runs the lower bound of the feedback_* arrays must be 1 for use in exchange data,
! the lower bound of the fbk_* arrays must be i_startblk for the use in global sum.

i_startblk = p_gcp%start_blk(grf_fbk_start_c,i_chidx)
i_endblk   = p_gcp%end_blk(min_rlcell_int,i_chidx)

ALLOCATE(fbk_tend(nproma, i_startblk:i_endblk))

i_startblk = 1

ALLOCATE(feedback_pres_tend(nproma, i_startblk:i_endblk))
ALLOCATE(feedback_temp_tend(nproma, nlev, i_startblk:i_endblk))
IF(ltransport) &
  ALLOCATE(feedback_tracer_tend(nproma, nlev, i_startblk:i_endblk, ntracer))

i_startblk = p_gep%start_blk(grf_fbk_start_e,i_chidx)
i_endblk   = p_gep%end_blk(min_rledge,i_chidx)

i_startblk = 1

ALLOCATE(feedback_vn_tend(nproma, nlev, i_startblk:i_endblk))

! Set pointers to index and coefficient fields for cell-based variables
iidx => p_gcp%child_idx
iblk => p_gcp%child_blk

IF (grf_scalfbk == 1) THEN
  p_fbkwgt    => p_grf%fbk_wgt_aw
ELSE
  p_fbkwgt    => p_grf%fbk_wgt_bln
ENDIF
IF (grf_tracfbk == 1) THEN
  p_fbkwgt_tr => p_grf%fbk_wgt_aw
ELSE
  p_fbkwgt_tr => p_grf%fbk_wgt_bln
ENDIF

p_fbarea    => p_gcp%area

! Preparation of feedback: compute child tendencies

IF (ltheta_dyn) THEN
  p_temp_prog   => p_child_prog%theta
  p_temp_save   => p_child_save%theta
ELSE
  p_temp_prog   => p_child_prog%temp
  p_temp_save   => p_child_save%temp
ENDIF

!$OMP PARALLEL PRIVATE(i_startblk,i_endblk)
!
! Part 1: cell-based variables
i_startblk = p_gcc%start_blk(3,1)
i_endblk   = p_gcc%end_blk(min_rlcell_int,i_nchdom)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk,jt) ICON_OMP_DEFAULT_SCHEDULE
DO jb = i_startblk, i_endblk

  CALL get_indices_c(p_pc, jb, i_startblk, i_endblk, &
                     i_startidx, i_endidx, 3, min_rlcell_int)

  DO jc = i_startidx, i_endidx
    p_child_tend%pres_sfc(jc,jb) = &
      p_child_prog%pres_sfc(jc,jb) - p_child_save%pres_sfc(jc,jb)
  ENDDO

  DO jk = 1, nlev
    DO jc = i_startidx, i_endidx
      p_child_tend%temp(jc,jk,jb) = &
        p_temp_prog(jc,jk,jb) - p_temp_save(jc,jk,jb)
    ENDDO
  ENDDO

  IF(ltransport) THEN
    DO jt = 1, ntracer
      DO jk = 1, nlev
        DO jc = i_startidx, i_endidx
          p_child_tend%tracer(jc,jk,jb,jt) = &
            p_child_prog%tracer(jc,jk,jb,jt) - &
            p_child_save%tracer(jc,jk,jb,jt)
        ENDDO
      ENDDO
    ENDDO
  ENDIF

ENDDO
!$OMP END DO

! Part 2: Velocity components
i_startblk = p_gec%start_blk(4,1)
i_endblk   = p_gec%end_blk(min_rledge_int-1,i_nchdom)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,je,jk) ICON_OMP_DEFAULT_SCHEDULE
DO jb = i_startblk, i_endblk

  CALL get_indices_e(p_pc, jb, i_startblk, i_endblk, &
                     i_startidx, i_endidx, 4, min_rledge_int-1)

  DO jk = 1, nlev
    DO je = i_startidx, i_endidx
      p_child_tend%vn(je,jk,jb) = &
        p_child_prog%vn(je,jk,jb) - p_child_save%vn(je,jk,jb)
    ENDDO
  ENDDO

ENDDO
!$OMP END DO

! Compute feedback tendency for pressure

! Start/End block in the [local] parent domain
i_startblk = p_gcp%start_blk(grf_fbk_start_c,i_chidx)
i_endblk   = p_gcp%end_blk(min_rlcell_int,i_chidx)


!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc) ICON_OMP_DEFAULT_SCHEDULE
DO jb = i_startblk, i_endblk

  CALL get_indices_c(p_pp, jb, i_startblk, i_endblk, &
                     i_startidx, i_endidx, grf_fbk_start_c, min_rlcell_int)

  fbk_tend(:,jb) = 0._wp

  DO jc = i_startidx, i_endidx

      feedback_pres_tend(jc,jb) =                                             &
        p_child_tend%pres_sfc(iidx(jc,jb,1),iblk(jc,jb,1))*p_fbkwgt(jc,jb,1) + &
        p_child_tend%pres_sfc(iidx(jc,jb,2),iblk(jc,jb,2))*p_fbkwgt(jc,jb,2) + &
        p_child_tend%pres_sfc(iidx(jc,jb,3),iblk(jc,jb,3))*p_fbkwgt(jc,jb,3) + &
        p_child_tend%pres_sfc(iidx(jc,jb,4),iblk(jc,jb,4))*p_fbkwgt(jc,jb,4)

      fbk_tend(jc,jb) = feedback_pres_tend(jc,jb)*p_fbarea(jc,jb)

  ENDDO
  ! Please note that the local parent has no halos, so no need to care about!
ENDDO
!$OMP END DO

! Calculate the area-weighted sum of (p_parent_prog%pres_sfc-p_parent_save%pres_sfc)
! as well as the feedback area - this has always to be done in the parent patch

! Start/End block in the parent domain
i_startblk = p_patch(jgp)%cells%start_blk(grf_fbk_start_c,i_chidx)
i_endblk   = p_patch(jgp)%cells%end_blk(min_rlcell_int,i_chidx)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc) ICON_OMP_DEFAULT_SCHEDULE
DO jb = i_startblk, i_endblk

  parent_tend(:,jb) = 0._wp

  CALL get_indices_c(p_patch(jgp), jb, i_startblk, i_endblk, &
                     i_startidx, i_endidx, grf_fbk_start_c, min_rlcell_int)

  DO jc = i_startidx, i_endidx
    parent_tend(jc,jb) =  &
      (p_parent_prog%pres_sfc(jc,jb)-p_parent_save%pres_sfc(jc,jb))*p_patch(jgp)%cells%area(jc,jb)
  ENDDO
  ! Sum must be taken over inner domain only, the halos must be masked out!
  WHERE(.NOT.p_patch(jgp)%cells%decomp_info%owner_mask(:,jb)) parent_tend(:,jb) = 0._wp
ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

! compute correction for global mass conservation
tendency_corr = (global_sum_array2(parent_tend) -   &
                 global_sum_array2(fbk_tend))/p_grf%fbk_dom_area(i_chidx)

IF (ltheta_dyn) THEN
  p_temp_prog   => p_parent_prog%theta
  p_temp_save   => p_parent_save%theta
ELSE
  p_temp_prog   => p_parent_prog%temp
  p_temp_save   => p_parent_save%temp
ENDIF


!$OMP PARALLEL PRIVATE(i_startblk,i_endblk,jt)
i_startblk = p_gcp%start_blk(grf_fbk_start_c,i_chidx)
i_endblk   = p_gcp%end_blk(min_rlcell_int,i_chidx)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk) ICON_OMP_DEFAULT_SCHEDULE
DO jb = i_startblk, i_endblk

  CALL get_indices_c(p_pp, jb, i_startblk, i_endblk, &
                     i_startidx, i_endidx, grf_fbk_start_c, min_rlcell_int)

  DO jc = i_startidx, i_endidx
    feedback_pres_tend(jc,jb) = feedback_pres_tend(jc,jb) + tendency_corr
  ENDDO

#ifdef _URD
!CDIR UNROLL=_URD
#endif
  DO jk = 1, nlev
    DO jc = i_startidx, i_endidx

      feedback_temp_tend(jc,jk,jb) =                                             &
        p_child_tend%temp(iidx(jc,jb,1),jk,iblk(jc,jb,1))*p_fbkwgt(jc,jb,1) + &
        p_child_tend%temp(iidx(jc,jb,2),jk,iblk(jc,jb,2))*p_fbkwgt(jc,jb,2) + &
        p_child_tend%temp(iidx(jc,jb,3),jk,iblk(jc,jb,3))*p_fbkwgt(jc,jb,3) + &
        p_child_tend%temp(iidx(jc,jb,4),jk,iblk(jc,jb,4))*p_fbkwgt(jc,jb,4)
    ENDDO
  ENDDO
ENDDO
!$OMP END DO

! Tracer feedback
IF (ltransport) THEN

DO jt = 1, ntracer

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk) ICON_OMP_DEFAULT_SCHEDULE
  DO jb = i_startblk, i_endblk

    CALL get_indices_c(p_pp, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, grf_fbk_start_c, min_rlcell_int)

#ifdef _URD
!CDIR UNROLL=_URD
#endif
    DO jk = 1, nlev
      DO jc = i_startidx, i_endidx

        feedback_tracer_tend(jc,jk,jb,jt) =                                              &
          p_child_tend%tracer(iidx(jc,jb,1),jk,iblk(jc,jb,1),jt)*p_fbkwgt_tr(jc,jb,1) + &
          p_child_tend%tracer(iidx(jc,jb,2),jk,iblk(jc,jb,2),jt)*p_fbkwgt_tr(jc,jb,2) + &
          p_child_tend%tracer(iidx(jc,jb,3),jk,iblk(jc,jb,3),jt)*p_fbkwgt_tr(jc,jb,3) + &
          p_child_tend%tracer(iidx(jc,jb,4),jk,iblk(jc,jb,4),jt)*p_fbkwgt_tr(jc,jb,4)
      ENDDO
    ENDDO
  ENDDO
!$OMP END DO
ENDDO

ENDIF
!$OMP END PARALLEL


IF (grf_velfbk == 2) THEN ! Interpolate velocity tendencies in child domain to vertices
  CALL rbf_vec_interpol_vertex( p_child_tend%vn, p_pc,  &
                                p_intc, z_u, z_v)
ENDIF

! Set pointers to index and coefficient fields
! (this needs to be done outside a parallel section!)
iidx => p_gep%child_idx
iblk => p_gep%child_blk

p_fbkwgt => p_grf%fbk_wgt_e

iidxv => p_gec%vertex_idx
iblkv => p_gec%vertex_blk

!$OMP PARALLEL PRIVATE(i_startblk,i_endblk)

! Velocity feedback
IF (grf_velfbk == 1) THEN ! Averaging weighted with child edge lenghts

  i_startblk = p_gep%start_blk(grf_fbk_start_e,i_chidx)
  i_endblk   = p_gep%end_blk(min_rledge_int,i_chidx)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,je,jk) ICON_OMP_DEFAULT_SCHEDULE
  DO jb = i_startblk, i_endblk

    CALL get_indices_e(p_pp, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, grf_fbk_start_e, min_rledge_int)

#ifdef _URD
!CDIR UNROLL=_URD
#endif
    DO jk = 1, nlev
      DO je = i_startidx, i_endidx

        feedback_vn_tend(je,jk,jb) =                                                 &
          p_child_tend%vn(iidx(je,jb,1),jk,iblk(je,jb,1))*p_fbkwgt(je,jb,1) + &
          p_child_tend%vn(iidx(je,jb,2),jk,iblk(je,jb,2))*p_fbkwgt(je,jb,2)
      ENDDO
    ENDDO
  ENDDO
!$OMP END DO

ELSE IF (grf_velfbk == 2) THEN ! Second-order interpolation of normal velocities
                               ! using RBF reconstruction to child vertices

  ! Projection of reconstructed velocities to the edge normals
  i_startblk = p_gec%start_blk(4,1)
  i_endblk   = p_gec%end_blk(min_rledge,i_nchdom)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,je,jk) ICON_OMP_DEFAULT_SCHEDULE
  DO jb = i_startblk, i_endblk

    CALL get_indices_e(p_pc, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, 4, min_rledge)

#ifdef _URD
!CDIR UNROLL=_URD
#endif
    DO jk = 1, nlev
      DO je = i_startidx, i_endidx

        vn_aux(je,jk,jb,1) = z_u(iidxv(je,jb,1),jk,iblkv(je,jb,1)) &
             &             * p_gec%primal_normal_vert(je,jb,1)%v1  &
             &             + z_v(iidxv(je,jb,1),jk,iblkv(je,jb,1)) &
             &             * p_gec%primal_normal_vert(je,jb,1)%v2
        vn_aux(je,jk,jb,2) = z_u(iidxv(je,jb,2),jk,iblkv(je,jb,2)) &
             &             * p_gec%primal_normal_vert(je,jb,2)%v1  &
             &             + z_v(iidxv(je,jb,2),jk,iblkv(je,jb,2)) &
             &             * p_gec%primal_normal_vert(je,jb,2)%v2
      ENDDO
    ENDDO
  ENDDO
!$OMP END DO

  i_startblk = p_gep%start_blk(grf_fbk_start_e,i_chidx)
  i_endblk   = p_gep%end_blk(min_rledge_int,i_chidx)


!$OMP DO PRIVATE(jb,i_startidx,i_endidx,je,jk) ICON_OMP_DEFAULT_SCHEDULE
  DO jb = i_startblk, i_endblk

    CALL get_indices_e(p_pp, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, grf_fbk_start_e, min_rledge_int)

#ifdef _URD
!CDIR UNROLL=_URD
#endif
    DO jk = 1, nlev
      DO je = i_startidx, i_endidx

        feedback_vn_tend(je,jk,jb) =                                             &
        ( p_fbkwgt(je,jb,1)*(p_child_tend%vn(iidx(je,jb,1),jk,iblk(je,jb,1))  + &
                             p_child_tend%vn(iidx(je,jb,2),jk,iblk(je,jb,2))) + &
          p_fbkwgt(je,jb,3)*vn_aux(iidx(je,jb,1),jk,iblk(je,jb,1),1) +                  &
          p_fbkwgt(je,jb,4)*vn_aux(iidx(je,jb,1),jk,iblk(je,jb,1),2) +                  &
          p_fbkwgt(je,jb,5)*vn_aux(iidx(je,jb,2),jk,iblk(je,jb,2),1) +                  &
          p_fbkwgt(je,jb,6)*vn_aux(iidx(je,jb,2),jk,iblk(je,jb,2),2) )
      ENDDO
    ENDDO
  ENDDO
!$OMP END DO

ENDIF
!$OMP END PARALLEL

  CALL exchange_data(p_pp%comm_pat_loc_to_glb_c_fbk, &
                     RECV=p_parent_prog%pres_sfc, &
                     SEND=feedback_pres_tend, &
                     ADD =p_parent_save%pres_sfc)

  CALL exchange_data(p_pp%comm_pat_loc_to_glb_c_fbk, &
                   RECV=p_temp_prog, &
                   SEND=feedback_temp_tend, &
                   ADD= p_temp_save)

  CALL exchange_data(p_pp%comm_pat_loc_to_glb_e_fbk, &
                     RECV=p_parent_prog%vn, &
                     SEND=feedback_vn_tend, &
                     ADD =p_parent_save%vn)

IF (ltransport) THEN
DO jt = 1, ntracer
  CALL exchange_data(p_pp%comm_pat_loc_to_glb_c_fbk, &
                   RECV=p_parent_prog%tracer(:,:,:,jt), &
                   SEND=feedback_tracer_tend(:,:,:,jt), &
                   ADD =p_parent_save%tracer(:,:,:,jt))
ENDDO

ENDIF

CALL sync_patch_array(SYNC_C,p_patch(jgp),p_parent_prog%pres_sfc)
CALL sync_patch_array(SYNC_C,p_patch(jgp),p_temp_prog)
CALL sync_patch_array(SYNC_E,p_patch(jgp),p_parent_prog%vn)


IF (ltransport) THEN
  DO jt = 1, ntracer
    CALL sync_patch_array(SYNC_C,p_patch(jgp),p_parent_prog%tracer(:,:,:,jt))
  ENDDO
ENDIF

DEALLOCATE(feedback_pres_tend)
DEALLOCATE(feedback_temp_tend)
DEALLOCATE(feedback_vn_tend)
IF (ltransport) DEALLOCATE(feedback_tracer_tend)


END SUBROUTINE feedback


END MODULE mo_hierarchy_management_intp
