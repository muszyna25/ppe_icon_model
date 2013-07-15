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
!! $Id: n/a$
!!

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_grf_ubcintp
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
USE mo_impl_constants,      ONLY: min_rlcell_int, min_rledge_int
USE mo_model_domain,        ONLY: t_patch
USE mo_intp_data_strc,      ONLY: t_int_state
USE mo_impl_constants_grf,  ONLY: grf_nudgintp_start_c, grf_nudgintp_start_e
USE mo_parallel_config,     ONLY: nproma,p_test_run 
USE mo_loopindices,         ONLY: get_indices_c, get_indices_e
USE mo_mpi,                 ONLY: my_process_is_mpi_parallel, my_process_is_mpi_seq

USE mo_communication,       ONLY: exchange_data_grf, exchange_data_grf_old

USE mo_grf_intp_data_strc


IMPLICIT NONE

PRIVATE

CHARACTER(len=*), PARAMETER :: version = '$Id$'

PUBLIC :: interpol_vec_ubc, interpol_scal_ubc ! , interpol_vec_ubc_old, interpol_scal_ubc_old

CONTAINS

!
!>
!! Performs gradient-based interpolation from parent edges to child edges 1 and 2
!! and RBF/IDW-based interpolation to child edges 3 and 4 for the upper boundary
!! condition needed for vertical nesting
!!
!!
!! @par Revision History
!! Developed  by Guenther Zaengl, DWD (2011-03-14)
!!
SUBROUTINE interpol_vec_ubc_old (ptr_pp, ptr_pc, ptr_grf, i_chidx, &
                             p_vn_in, p_vn_out)

TYPE(t_patch), TARGET, INTENT(in) :: ptr_pp
TYPE(t_patch), TARGET, INTENT(in) :: ptr_pc

! Indices of source points and interpolation coefficients
TYPE(t_gridref_single_state),   TARGET, INTENT(IN)    ::  ptr_grf

! input: delta vn at interface level
REAL(wp), INTENT(IN) :: p_vn_in(:,:)  ! dim: (nproma,nblks_e)

! child domain index as seen from parent domain
INTEGER, INTENT(IN) :: i_chidx

! reconstructed shift of edge-normal wind component at top
REAL(wp),INTENT(INOUT) :: p_vn_out(nproma,1,ptr_pc%nblks_e)


INTEGER :: jb, je                    ! loop indices
INTEGER :: i_startblk                ! start block
INTEGER :: i_endblk                  ! end index
INTEGER :: i_startidx                ! start index
INTEGER :: i_endidx                  ! end index
INTEGER :: i_nchdom                  ! number of child domains

INTEGER :: nsendtot, nrecvtot        ! for MPI communication call


REAL(wp) :: vn_aux(nproma,1,ptr_pp%edges%start_blk(grf_nudgintp_start_e,i_chidx): &
                   MAX(ptr_pp%edges%start_blk(grf_nudgintp_start_e,i_chidx),      &
                       ptr_pp%edges%end_blk(min_rledge_int,i_chidx)),4)

! Pointers to index fields
INTEGER,  DIMENSION(:,:,:),   POINTER :: iidx_1a, iblk_1a, iidx_1b, iblk_1b, &
                                         iidx_2a, iblk_2a, iidx_2b, iblk_2b
INTEGER,  DIMENSION(:,:,:),   POINTER :: icheidx, icheblk

!-----------------------------------------------------------------------
IF (p_test_run) vn_aux=0._wp

! Set pointers to index and coefficient fields

iidx_1a => ptr_grf%grf_vec_ind_1a
iblk_1a => ptr_grf%grf_vec_blk_1a
iidx_1b => ptr_grf%grf_vec_ind_1b
iblk_1b => ptr_grf%grf_vec_blk_1b

iidx_2a => ptr_grf%grf_vec_ind_2a
iblk_2a => ptr_grf%grf_vec_blk_2a
iidx_2b => ptr_grf%grf_vec_ind_2b
iblk_2b => ptr_grf%grf_vec_blk_2b

icheidx => ptr_pp%edges%child_idx
icheblk => ptr_pp%edges%child_blk

! Number of child domains of nested domain
i_nchdom = MAX(1,ptr_pc%n_childdom)

!$OMP PARALLEL PRIVATE (i_startblk,i_endblk)

! Start and end blocks for which interpolation is needed
i_startblk = ptr_pp%edges%start_blk(grf_nudgintp_start_e,i_chidx)
i_endblk   = ptr_pp%edges%end_blk(min_rledge_int,i_chidx)

!$OMP DO PRIVATE(jb,je,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
DO jb =  i_startblk, i_endblk

  CALL get_indices_e(ptr_pp, jb, i_startblk, i_endblk, &
   i_startidx, i_endidx, grf_nudgintp_start_e, min_rledge_int, i_chidx)

! child edge 1
  DO je = i_startidx, i_endidx
    vn_aux(je,1,jb,1) = ptr_grf%grf_vec_coeff_1a(1,je,jb) * &
      p_vn_in(iidx_1a(je,1,jb),iblk_1a(je,1,jb)) +          &
      ptr_grf%grf_vec_coeff_1a(2,je,jb) *                   &
      p_vn_in(iidx_1a(je,2,jb),iblk_1a(je,2,jb)) +          &
      ptr_grf%grf_vec_coeff_1a(3,je,jb) *                   &
      p_vn_in(iidx_1a(je,3,jb),iblk_1a(je,3,jb)) +          &
      ptr_grf%grf_vec_coeff_1a(4,je,jb) *                   &
      p_vn_in(iidx_1a(je,4,jb),iblk_1a(je,4,jb)) +          &
      ptr_grf%grf_vec_coeff_1a(5,je,jb) *                   &
      p_vn_in(iidx_1a(je,5,jb),iblk_1a(je,5,jb)) +          &
      ptr_grf%grf_vec_coeff_1a(6,je,jb) *                   &
      p_vn_in(iidx_1a(je,6,jb),iblk_1a(je,6,jb))

! child edge 2
    vn_aux(je,1,jb,2) = ptr_grf%grf_vec_coeff_1b(1,je,jb) * &
      p_vn_in(iidx_1b(je,1,jb),iblk_1b(je,1,jb)) +          &
      ptr_grf%grf_vec_coeff_1b(2,je,jb) *                   &
      p_vn_in(iidx_1b(je,2,jb),iblk_1b(je,2,jb)) +          &
      ptr_grf%grf_vec_coeff_1b(3,je,jb) *                   &
      p_vn_in(iidx_1b(je,3,jb),iblk_1b(je,3,jb)) +          &
      ptr_grf%grf_vec_coeff_1b(4,je,jb) *                   &
      p_vn_in(iidx_1b(je,4,jb),iblk_1b(je,4,jb)) +          &
      ptr_grf%grf_vec_coeff_1b(5,je,jb) *                   &
      p_vn_in(iidx_1b(je,5,jb),iblk_1b(je,5,jb)) +          &
      ptr_grf%grf_vec_coeff_1b(6,je,jb) *                   &
      p_vn_in(iidx_1b(je,6,jb),iblk_1b(je,6,jb))

! child edge 3
    vn_aux(je,1,jb,3) = ptr_grf%grf_vec_coeff_2a(1,je,jb) * &
      p_vn_in(iidx_2a(je,1,jb),iblk_2a(je,1,jb)) +          &
      ptr_grf%grf_vec_coeff_2a(2,je,jb) *                   &
      p_vn_in(iidx_2a(je,2,jb),iblk_2a(je,2,jb)) +          &
      ptr_grf%grf_vec_coeff_2a(3,je,jb) *                   &
      p_vn_in(iidx_2a(je,3,jb),iblk_2a(je,3,jb)) +          &
      ptr_grf%grf_vec_coeff_2a(4,je,jb) *                   &
      p_vn_in(iidx_2a(je,4,jb),iblk_2a(je,4,jb)) +          &
      ptr_grf%grf_vec_coeff_2a(5,je,jb) *                   &
      p_vn_in(iidx_2a(je,5,jb),iblk_2a(je,5,jb))

! child edge 4
    vn_aux(je,1,jb,4) = ptr_grf%grf_vec_coeff_2b(1,je,jb) * &
      p_vn_in(iidx_2b(je,1,jb),iblk_2b(je,1,jb)) +          &
      ptr_grf%grf_vec_coeff_2b(2,je,jb) *                   &
      p_vn_in(iidx_2b(je,2,jb),iblk_2b(je,2,jb)) +          &
      ptr_grf%grf_vec_coeff_2b(3,je,jb) *                   &
      p_vn_in(iidx_2b(je,3,jb),iblk_2b(je,3,jb)) +          &
      ptr_grf%grf_vec_coeff_2b(4,je,jb) *                   &
      p_vn_in(iidx_2b(je,4,jb),iblk_2b(je,4,jb)) +          &
      ptr_grf%grf_vec_coeff_2b(5,je,jb) *                   &
      p_vn_in(iidx_2b(je,5,jb),iblk_2b(je,5,jb))
  ENDDO

ENDDO
!$OMP END DO

! Store results in p_vn_out

IF (my_process_is_mpi_seq()) THEN

!$OMP DO PRIVATE(jb,je,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
  DO jb =  i_startblk, i_endblk

  CALL get_indices_e(ptr_pp, jb, i_startblk, i_endblk, &
    i_startidx, i_endidx, grf_nudgintp_start_e, min_rledge_int, i_chidx)

!CDIR NODEP,VOVERTAKE,VOB
    DO je = i_startidx, i_endidx
      p_vn_out(icheidx(je,jb,1),1,icheblk(je,jb,1)) = vn_aux(je,1,jb,1)
      p_vn_out(icheidx(je,jb,2),1,icheblk(je,jb,2)) = vn_aux(je,1,jb,2)
      p_vn_out(icheidx(je,jb,3),1,icheblk(je,jb,3)) = vn_aux(je,1,jb,3)
      p_vn_out(icheidx(je,jb,4),1,icheblk(je,jb,4)) = vn_aux(je,1,jb,4)

    ENDDO

  ENDDO
!$OMP END DO NOWAIT

ENDIF ! not MPI-parallel

!$OMP END PARALLEL

IF (my_process_is_mpi_parallel()) THEN

  nsendtot = SUM(ptr_pc%comm_pat_interpol_vec_ubc(1:4)%n_send)
  nrecvtot = SUM(ptr_pc%comm_pat_interpol_vec_ubc(1:4)%n_recv)
  CALL exchange_data_grf_old(ptr_pc%comm_pat_interpol_vec_ubc,1,1,nsendtot,nrecvtot, &
                         RECV1=p_vn_out,SEND1=vn_aux,SEND_LBOUND3=LBOUND(vn_aux,3))

ENDIF

END SUBROUTINE interpol_vec_ubc_old


!-------------------------------------------------------------------------
!
!>
!! Performs interpolation of scalar upper boundary condition fields from parent 
!! cells to child cells using the 2D gradient at the cell center.
!!
!! @par Revision History
!! Developed  by Guenther Zaengl, DWD (2011-03-16)
!!
SUBROUTINE interpol_scal_ubc_old (ptr_pp, ptr_pc, ptr_int, ptr_grf, i_chidx, nfields, &
                              f3din, f3dout, llimit_nneg )
!
TYPE(t_patch), TARGET, INTENT(in) :: ptr_pp
TYPE(t_patch), TARGET, INTENT(in) :: ptr_pc


! Indices of source points and interpolation coefficients
TYPE(t_gridref_single_state),   TARGET, INTENT(IN)    ::  ptr_grf
TYPE(t_int_state),              TARGET, INTENT(IN)    ::  ptr_int

! child domain index as seen from parent domain
INTEGER, INTENT(IN) :: i_chidx

! number of 2D fields contained in input/output fields
INTEGER, INTENT(IN) :: nfields

! logical switch: if present and true, limit horizontal gradient so as to avoid negative values
LOGICAL, INTENT(IN), OPTIONAL :: llimit_nneg

! input scalar fields at cell points
REAL(wp), INTENT(IN) :: f3din(:,:,:) ! dim: (nproma,nfields,nblks_c)

! interpolated scalar output fields
REAL(wp), INTENT(INOUT) :: f3dout(:,:,:)

INTEGER :: jb, jc, jn                ! loop indices
INTEGER :: i_startblk                ! start block
INTEGER :: i_endblk                  ! end index
INTEGER :: i_startidx                ! start index
INTEGER :: i_endidx                  ! end index

INTEGER :: nsendtot, nrecvtot        ! for MPI communication call

LOGICAL :: l_limit_nneg              ! local variable corresponding to llimit_nneg

! Auxiliary fields
REAL(wp), DIMENSION(nproma,nfields) :: grad_x, grad_y, maxval_neighb, minval_neighb
REAL(wp) :: h_aux(nproma,nfields,ptr_pp%cells%start_blk(grf_nudgintp_start_c,i_chidx):  &
                  MAX(ptr_pp%cells%start_blk(grf_nudgintp_start_c,i_chidx),             &
                      ptr_pp%cells%end_blk(min_rlcell_int,i_chidx)),4)

REAL(wp) :: limfac1, limfac2, limfac, min_expval, max_expval, epsi, ovsht_fac, r_ovsht_fac, &
            relaxed_minval, relaxed_maxval

! Pointers to index and coefficient fields
INTEGER,  DIMENSION(:,:,:),   POINTER :: iidx, iblk, ichcidx, ichcblk
REAL(wp), DIMENSION(:,:,:,:), POINTER :: ptr_coeff, ptr_dist


!-----------------------------------------------------------------------


! Check if gradient limiting is required
IF (PRESENT(llimit_nneg)) THEN
  l_limit_nneg = llimit_nneg
ELSE
  l_limit_nneg = .FALSE.
ENDIF

epsi = 1.e-75_wp
ovsht_fac = 1.0_wp ! factor of allowed overshooting
r_ovsht_fac = 1._wp/ovsht_fac

! Start and end blocks for which scalar interpolation is needed
i_startblk = ptr_pp%cells%start_blk(grf_nudgintp_start_c,i_chidx)
i_endblk   = ptr_pp%cells%end_blk(min_rlcell_int,i_chidx)

! Pointers to index lists and interpolation coefficients for gradient computation
iidx => ptr_int%rbf_c2grad_idx
iblk => ptr_int%rbf_c2grad_blk

ptr_coeff => ptr_int%rbf_c2grad_coeff
ptr_dist  => ptr_grf%grf_dist_pc2cc

! child cell indices and blocks for non-MPI parent-to-child communication
ichcidx => ptr_pp%cells%child_idx
ichcblk => ptr_pp%cells%child_blk


!$OMP PARALLEL PRIVATE(jb,jn,i_startidx,i_endidx)


!$OMP DO PRIVATE (jc,grad_x,grad_y,min_expval,max_expval,limfac1,  &
!$OMP   limfac2,limfac,maxval_neighb,minval_neighb,relaxed_minval, &
!$OMP   relaxed_maxval) ICON_OMP_DEFAULT_SCHEDULE
  DO jb = i_startblk, i_endblk

    CALL get_indices_c(ptr_pp, jb, i_startblk, i_endblk, &
         i_startidx, i_endidx, grf_nudgintp_start_c, min_rlcell_int, i_chidx)

#ifdef __LOOP_EXCHANGE
    DO jc = i_startidx, i_endidx
      DO jn = 1, nfields
#else
!CDIR UNROLL=2
    DO jn = 1, nfields
      DO jc = i_startidx, i_endidx
#endif

        grad_x(jc,jn) =  &
          ptr_coeff(1,1,jc,jb)*f3din(jc,jn,jb) +                       &
          ptr_coeff(2,1,jc,jb)*f3din(iidx(2,jc,jb),jn,iblk(2,jc,jb)) + &
          ptr_coeff(3,1,jc,jb)*f3din(iidx(3,jc,jb),jn,iblk(3,jc,jb)) + &
          ptr_coeff(4,1,jc,jb)*f3din(iidx(4,jc,jb),jn,iblk(4,jc,jb)) + &
          ptr_coeff(5,1,jc,jb)*f3din(iidx(5,jc,jb),jn,iblk(5,jc,jb)) + &
          ptr_coeff(6,1,jc,jb)*f3din(iidx(6,jc,jb),jn,iblk(6,jc,jb)) + &
          ptr_coeff(7,1,jc,jb)*f3din(iidx(7,jc,jb),jn,iblk(7,jc,jb)) + &
          ptr_coeff(8,1,jc,jb)*f3din(iidx(8,jc,jb),jn,iblk(8,jc,jb)) + &
          ptr_coeff(9,1,jc,jb)*f3din(iidx(9,jc,jb),jn,iblk(9,jc,jb)) + &
          ptr_coeff(10,1,jc,jb)*f3din(iidx(10,jc,jb),jn,iblk(10,jc,jb))
        grad_y(jc,jn) =  &
          ptr_coeff(1,2,jc,jb)*f3din(jc,jn,jb) +                       &
          ptr_coeff(2,2,jc,jb)*f3din(iidx(2,jc,jb),jn,iblk(2,jc,jb)) + &
          ptr_coeff(3,2,jc,jb)*f3din(iidx(3,jc,jb),jn,iblk(3,jc,jb)) + &
          ptr_coeff(4,2,jc,jb)*f3din(iidx(4,jc,jb),jn,iblk(4,jc,jb)) + &
          ptr_coeff(5,2,jc,jb)*f3din(iidx(5,jc,jb),jn,iblk(5,jc,jb)) + &
          ptr_coeff(6,2,jc,jb)*f3din(iidx(6,jc,jb),jn,iblk(6,jc,jb)) + &
          ptr_coeff(7,2,jc,jb)*f3din(iidx(7,jc,jb),jn,iblk(7,jc,jb)) + &
          ptr_coeff(8,2,jc,jb)*f3din(iidx(8,jc,jb),jn,iblk(8,jc,jb)) + &
          ptr_coeff(9,2,jc,jb)*f3din(iidx(9,jc,jb),jn,iblk(9,jc,jb)) + &
          ptr_coeff(10,2,jc,jb)*f3din(iidx(10,jc,jb),jn,iblk(10,jc,jb))
        maxval_neighb(jc,jn) =                       &
          MAX(f3din(jc,jn,jb),                       &
              f3din(iidx(2,jc,jb),jn,iblk(2,jc,jb)), &
              f3din(iidx(3,jc,jb),jn,iblk(3,jc,jb)), &
              f3din(iidx(4,jc,jb),jn,iblk(4,jc,jb)), &
              f3din(iidx(5,jc,jb),jn,iblk(5,jc,jb)), &
              f3din(iidx(6,jc,jb),jn,iblk(6,jc,jb)), &
              f3din(iidx(7,jc,jb),jn,iblk(7,jc,jb)), &
              f3din(iidx(8,jc,jb),jn,iblk(8,jc,jb)), &
              f3din(iidx(9,jc,jb),jn,iblk(9,jc,jb)), &
              f3din(iidx(10,jc,jb),jn,iblk(10,jc,jb)))
        minval_neighb(jc,jn) =                       &
          MIN(f3din(jc,jn,jb),                       &
              f3din(iidx(2,jc,jb),jn,iblk(2,jc,jb)), &
              f3din(iidx(3,jc,jb),jn,iblk(3,jc,jb)), &
              f3din(iidx(4,jc,jb),jn,iblk(4,jc,jb)), &
              f3din(iidx(5,jc,jb),jn,iblk(5,jc,jb)), &
              f3din(iidx(6,jc,jb),jn,iblk(6,jc,jb)), &
              f3din(iidx(7,jc,jb),jn,iblk(7,jc,jb)), &
              f3din(iidx(8,jc,jb),jn,iblk(8,jc,jb)), &
              f3din(iidx(9,jc,jb),jn,iblk(9,jc,jb)), &
              f3din(iidx(10,jc,jb),jn,iblk(10,jc,jb)))

      ENDDO
    ENDDO

    DO jn = 1, nfields
      DO jc = i_startidx, i_endidx
        min_expval = MIN(grad_x(jc,jn)*ptr_dist(jc,1,1,jb) + &
                         grad_y(jc,jn)*ptr_dist(jc,1,2,jb),  &
                         grad_x(jc,jn)*ptr_dist(jc,2,1,jb) + &
                         grad_y(jc,jn)*ptr_dist(jc,2,2,jb),  &
                         grad_x(jc,jn)*ptr_dist(jc,3,1,jb) + &
                         grad_y(jc,jn)*ptr_dist(jc,3,2,jb),  &
                         grad_x(jc,jn)*ptr_dist(jc,4,1,jb) + &
                         grad_y(jc,jn)*ptr_dist(jc,4,2,jb),  &
                         -1.e-80_wp )
        max_expval = MAX(grad_x(jc,jn)*ptr_dist(jc,1,1,jb) + &
                         grad_y(jc,jn)*ptr_dist(jc,1,2,jb),  &
                         grad_x(jc,jn)*ptr_dist(jc,2,1,jb) + &
                         grad_y(jc,jn)*ptr_dist(jc,2,2,jb),  &
                         grad_x(jc,jn)*ptr_dist(jc,3,1,jb) + &
                         grad_y(jc,jn)*ptr_dist(jc,3,2,jb),  &
                         grad_x(jc,jn)*ptr_dist(jc,4,1,jb) + &
                         grad_y(jc,jn)*ptr_dist(jc,4,2,jb),  &
                         1.e-80_wp )

        limfac1 = 1._wp
        limfac2 = 1._wp
        ! Allow a limited amount of over-/undershooting in the downscaled fields
        IF (minval_neighb(jc,jn) > 0._wp) THEN
          relaxed_minval = r_ovsht_fac*minval_neighb(jc,jn)
        ELSE
          relaxed_minval = ovsht_fac*minval_neighb(jc,jn)
        ENDIF
        IF (maxval_neighb(jc,jn) > 0._wp) THEN
          relaxed_maxval = ovsht_fac*maxval_neighb(jc,jn)
        ELSE
          relaxed_maxval = r_ovsht_fac*maxval_neighb(jc,jn)
        ENDIF

        IF (f3din(jc,jn,jb) + min_expval < relaxed_minval-epsi) THEN
          limfac1 = ABS((relaxed_minval-f3din(jc,jn,jb))/min_expval)
        ENDIF
        IF (f3din(jc,jn,jb) + max_expval > relaxed_maxval+epsi) THEN
          limfac2 = ABS((relaxed_maxval-f3din(jc,jn,jb))/max_expval)
        ENDIF
        limfac = MIN(limfac1,limfac2)

        grad_x(jc,jn) = grad_x(jc,jn)*limfac
        grad_y(jc,jn) = grad_y(jc,jn)*limfac

      ENDDO
    ENDDO

    IF (l_limit_nneg) THEN
      DO jn = 1, nfields
        DO jc = i_startidx, i_endidx

          h_aux(jc,jn,jb,1) = MAX(0._wp, f3din(jc,jn,jb) + &
            grad_x(jc,jn)*ptr_dist(jc,1,1,jb)            + &
            grad_y(jc,jn)*ptr_dist(jc,1,2,jb))
          h_aux(jc,jn,jb,2) = MAX(0._wp, f3din(jc,jn,jb) + &
            grad_x(jc,jn)*ptr_dist(jc,2,1,jb)            + &
            grad_y(jc,jn)*ptr_dist(jc,2,2,jb))
          h_aux(jc,jn,jb,3) = MAX(0._wp, f3din(jc,jn,jb) + &
            grad_x(jc,jn)*ptr_dist(jc,3,1,jb)            + &
            grad_y(jc,jn)*ptr_dist(jc,3,2,jb))
          h_aux(jc,jn,jb,4) = MAX(0._wp, f3din(jc,jn,jb) + &
            grad_x(jc,jn)*ptr_dist(jc,4,1,jb)            + &
            grad_y(jc,jn)*ptr_dist(jc,4,2,jb))

        ENDDO
      ENDDO
    ELSE
      DO jn = 1, nfields
        DO jc = i_startidx, i_endidx

          h_aux(jc,jn,jb,1) = f3din(jc,jn,jb) + &
            grad_x(jc,jn)*ptr_dist(jc,1,1,jb) + &
            grad_y(jc,jn)*ptr_dist(jc,1,2,jb)
          h_aux(jc,jn,jb,2) = f3din(jc,jn,jb) + &
            grad_x(jc,jn)*ptr_dist(jc,2,1,jb) + &
            grad_y(jc,jn)*ptr_dist(jc,2,2,jb)
          h_aux(jc,jn,jb,3) = f3din(jc,jn,jb) + &
            grad_x(jc,jn)*ptr_dist(jc,3,1,jb) + &
            grad_y(jc,jn)*ptr_dist(jc,3,2,jb)
          h_aux(jc,jn,jb,4) = f3din(jc,jn,jb) + &
            grad_x(jc,jn)*ptr_dist(jc,4,1,jb) + &
            grad_y(jc,jn)*ptr_dist(jc,4,2,jb)

        ENDDO
      ENDDO
    ENDIF

  ENDDO ! blocks
!$OMP END DO


  ! Store results in p_out

  IF (my_process_is_mpi_seq()) THEN

!$OMP DO PRIVATE (jc,jn) ICON_OMP_DEFAULT_SCHEDULE
    DO jb =  i_startblk, i_endblk

      CALL get_indices_c(ptr_pp, jb, i_startblk, i_endblk, &
           i_startidx, i_endidx, grf_nudgintp_start_c, min_rlcell_int, i_chidx)

#ifdef __LOOP_EXCHANGE
      DO jc = i_startidx, i_endidx
        DO jn = 1, nfields
#else
      DO jn = 1, nfields
!CDIR NODEP,VOVERTAKE,VOB
        DO jc = i_startidx, i_endidx
#endif

          f3dout(ichcidx(jc,jb,1),jn,ichcblk(jc,jb,1)) = h_aux(jc,jn,jb,1)
          f3dout(ichcidx(jc,jb,2),jn,ichcblk(jc,jb,2)) = h_aux(jc,jn,jb,2)
          f3dout(ichcidx(jc,jb,3),jn,ichcblk(jc,jb,3)) = h_aux(jc,jn,jb,3)
          f3dout(ichcidx(jc,jb,4),jn,ichcblk(jc,jb,4)) = h_aux(jc,jn,jb,4)

        ENDDO
      ENDDO

    ENDDO
!$OMP END DO NOWAIT

  ENDIF
!$OMP END PARALLEL

IF (my_process_is_mpi_parallel()) THEN

  nsendtot = SUM(ptr_pc%comm_pat_interpol_scal_ubc(1:4)%n_send)
  nrecvtot = SUM(ptr_pc%comm_pat_interpol_scal_ubc(1:4)%n_recv)


  CALL exchange_data_grf_old(ptr_pc%comm_pat_interpol_scal_ubc,1,nfields,nsendtot, &
                         nrecvtot,RECV1=f3dout,SEND1=h_aux,SEND_LBOUND3=LBOUND(h_aux,3))

ENDIF

END SUBROUTINE interpol_scal_ubc_old


!
!>
!! Performs RBF-based interpolation from parent edges to child edges
!! for the upper boundary condition needed for vertical nesting
!!
!!
!! @par Revision History
!! Developed  by Guenther Zaengl, DWD (2011-03-14)
!!
SUBROUTINE interpol_vec_ubc(p_pp, p_pc, p_grf, p_vn_in, p_vn_out)
  !
  TYPE(t_patch), TARGET, INTENT(in) :: p_pp
  TYPE(t_patch), TARGET, INTENT(in) :: p_pc

  ! input: delta vn at interface level
  REAL(wp), INTENT(IN) :: p_vn_in(:,:) ! dim: (nproma,nblks_e)

  ! Indices of source points and interpolation coefficients
  TYPE(t_gridref_single_state),   TARGET, INTENT(IN) ::  p_grf

  ! reconstructed edge-normal wind component
  REAL(wp),INTENT(INOUT) :: p_vn_out(nproma,1,p_pc%nblks_e)


  INTEGER :: jb, je                    ! loop indices
  INTEGER :: nproma_ubcintp, nblks_ubcintp, npromz_ubcintp, nlen, nshift
  INTEGER :: ib, ie

  INTEGER :: nsendtot, nrecvtot        ! for MPI communication call

  REAL(wp) :: vn_aux(1,p_grf%npoints_ubcintp_e,4)

  ! Pointers to index fields/lists
  INTEGER,  DIMENSION(:,:,:), POINTER :: icheidx, icheblk
  INTEGER,  DIMENSION(:,:),   POINTER :: iidx, iblk

  INTEGER :: nlev_c       !< number of vertical full levels (child domain)
  !-----------------------------------------------------------------------

  ! Set pointers to child edge indices (needed for non-MPI case only)
  icheidx => p_pp%edges%child_idx
  icheblk => p_pp%edges%child_blk

  ! Set pointers to index lists
  iidx    => p_grf%idxlist_ubcintp_e
  iblk    => p_grf%blklist_ubcintp_e

  ! Compute values for dynamic nproma blocking
  nproma_ubcintp = MIN(nproma,256)
  nblks_ubcintp  = INT(p_grf%npoints_ubcintp_e/nproma_ubcintp)
  npromz_ubcintp = MOD(p_grf%npoints_ubcintp_e,nproma_ubcintp)
  IF (npromz_ubcintp > 0) THEN
    nblks_ubcintp = nblks_ubcintp + 1
  ELSE
    npromz_ubcintp = nproma_ubcintp
  ENDIF

!$OMP PARALLEL
!$OMP DO PRIVATE (jb,je,nlen,nshift) ICON_OMP_DEFAULT_SCHEDULE
  DO jb = 1, nblks_ubcintp
    IF (jb == nblks_ubcintp) THEN
      nlen = npromz_ubcintp
    ELSE
      nlen = nproma_ubcintp
    ENDIF
    nshift = (jb-1)*nproma_ubcintp

    DO je = nshift+1, nshift+nlen

      ! child edge 1
      vn_aux(1,je,1) = p_grf%coeff_ubcintp_e12(1,je) * &
        p_vn_in(iidx(1,je),iblk(1,je)) +               &
        p_grf%coeff_ubcintp_e12(2,je) *                &
        p_vn_in(iidx(2,je),iblk(2,je)) +               &
        p_grf%coeff_ubcintp_e12(3,je)*                 &
        p_vn_in(iidx(6,je),iblk(6,je)) +               &
        p_grf%coeff_ubcintp_e12(4,je)*                 &
        p_vn_in(iidx(10,je),iblk(10,je)) +             &
        p_grf%coeff_ubcintp_e12(5,je) *                &
        p_vn_in(iidx(11,je),iblk(11,je)) +             &
        p_grf%coeff_ubcintp_e12(6,je)*                 &
        p_vn_in(iidx(12,je),iblk(12,je)) 

      ! child edge 2
      vn_aux(1,je,2) = p_grf%coeff_ubcintp_e12(7,je) * &
        p_vn_in(iidx(1,je),iblk(1,je)) +               &
        p_grf%coeff_ubcintp_e12(8,je) *                &
        p_vn_in(iidx(3,je),iblk(3,je)) +               &
        p_grf%coeff_ubcintp_e12(9,je)*                 &
        p_vn_in(iidx(7,je),iblk(7,je)) +               &
        p_grf%coeff_ubcintp_e12(10,je)*                &
        p_vn_in(iidx(13,je),iblk(13,je)) +             &
        p_grf%coeff_ubcintp_e12(11,je) *               &
        p_vn_in(iidx(14,je),iblk(14,je)) +             &
        p_grf%coeff_ubcintp_e12(12,je)*                &
        p_vn_in(iidx(15,je),iblk(15,je)) 

      ! child edge 3
      vn_aux(1,je,3) = p_grf%coeff_ubcintp_e34(1,je) * &
        p_vn_in(iidx(1,je),iblk(1,je)) +               &
        p_grf%coeff_ubcintp_e34(2,je) *                &
        p_vn_in(iidx(2,je),iblk(2,je)) +               &
        p_grf%coeff_ubcintp_e34(3,je)*                 &
        p_vn_in(iidx(3,je),iblk(3,je)) +               &
        p_grf%coeff_ubcintp_e34(4,je)*                 &
        p_vn_in(iidx(4,je),iblk(4,je)) +               &
        p_grf%coeff_ubcintp_e34(5,je) *                &
        p_vn_in(iidx(5,je),iblk(5,je)) 

      ! child edge 4
      IF (p_pp%edges%refin_ctrl(iidx(1,je),iblk(1,je)) == -1) CYCLE
      vn_aux(1,je,4) = p_grf%coeff_ubcintp_e34(6,je) * &
        p_vn_in(iidx(1,je),iblk(1,je)) +               &
        p_grf%coeff_ubcintp_e34(7,je) *                &
        p_vn_in(iidx(6,je),iblk(6,je)) +               &
        p_grf%coeff_ubcintp_e34(8,je)*                 &
        p_vn_in(iidx(7,je),iblk(7,je)) +               &
        p_grf%coeff_ubcintp_e34(9,je)*                 &
        p_vn_in(iidx(8,je),iblk(8,je)) +               &
        p_grf%coeff_ubcintp_e34(10,je) *               &
        p_vn_in(iidx(9,je),iblk(9,je)) 

    ENDDO
  ENDDO
!$OMP END DO


  ! Store results in p_vn_out

  IF (my_process_is_mpi_seq()) THEN

!$OMP DO PRIVATE (jb,je,nlen,nshift,ib,ie) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = 1, nblks_ubcintp
      IF (jb == nblks_ubcintp) THEN
        nlen = npromz_ubcintp
      ELSE
        nlen = nproma_ubcintp
      ENDIF
      nshift = (jb-1)*nproma_ubcintp

      DO je = nshift+1, nshift+nlen
        ie = iidx(1,je)
        ib = iblk(1,je)

        p_vn_out(icheidx(ie,ib,1),1,icheblk(ie,ib,1)) = vn_aux(1,je,1)
        p_vn_out(icheidx(ie,ib,2),1,icheblk(ie,ib,2)) = vn_aux(1,je,2)
        p_vn_out(icheidx(ie,ib,3),1,icheblk(ie,ib,3)) = vn_aux(1,je,3)
        p_vn_out(icheidx(ie,ib,4),1,icheblk(ie,ib,4)) = vn_aux(1,je,4)

      ENDDO
    ENDDO
!$OMP END DO NOWAIT

  ENDIF
!$OMP END PARALLEL

  IF (my_process_is_mpi_parallel()) THEN

    nsendtot = SUM(p_pc%comm_pat_interpol_vec_ubc(1:4)%n_send)
    nrecvtot = SUM(p_pc%comm_pat_interpol_vec_ubc(1:4)%n_recv)
    CALL exchange_data_grf(p_pc%comm_pat_interpol_vec_ubc,1,1,nsendtot,nrecvtot, &
                           RECV1=p_vn_out,SEND1=vn_aux)

  ENDIF


END SUBROUTINE interpol_vec_ubc


!-------------------------------------------------------------------------
!
!>
!! Performs interpolation of scalar upper boundary condition fields from parent 
!! cells to child cells using the 2D gradient at the cell center.
!!
!! @par Revision History
!! Developed  by Guenther Zaengl, DWD (2011-03-16)
!!
SUBROUTINE interpol_scal_ubc(p_pp, p_pc, p_grf, nfields, f3din, f3dout, llimit_nneg)

  !
  TYPE(t_patch), TARGET, INTENT(in) :: p_pp
  TYPE(t_patch), TARGET, INTENT(in) :: p_pc

  ! Indices of source points and interpolation coefficients
  TYPE(t_gridref_single_state),   TARGET, INTENT(IN)    ::  p_grf

  ! number of fields provided on input (middle index of I/O fields for efficiency optimization)
  INTEGER, INTENT(IN) :: nfields

  ! logical switch: if present and true, limit horizontal gradient so as to avoid negative values
  LOGICAL, INTENT(IN), OPTIONAL :: llimit_nneg

  ! input scalar fields at cell points
  REAL(wp), INTENT(IN) :: f3din(:,:,:) ! dim: (nproma,nfields,nblks_c)

  ! interpolated scalar output fields
  REAL(wp), INTENT(INOUT) :: f3dout(:,:,:)

  INTEGER :: jb, jc, jn        ! loop indices
  INTEGER :: nproma_ubcintp, nblks_ubcintp, npromz_ubcintp, nlen, nshift
  INTEGER :: ib, ic
  INTEGER :: nsendtot, nrecvtot ! for MPI communication call

  LOGICAL :: l_limit_nneg     ! local variable corresponding to llimit_nneg

  ! Auxiliary fields
  REAL(wp), DIMENSION(nfields,p_grf%npoints_ubcintp_c) :: &
    grad_x, grad_y, maxval_neighb, minval_neighb, val_ctr
  REAL(wp) :: h_aux(nfields,p_grf%npoints_ubcintp_c,4)
  REAL(wp) :: limfac1, limfac2, limfac, min_expval, max_expval, epsi, ovsht_fac, r_ovsht_fac, &
              relaxed_minval, relaxed_maxval

  ! Pointers to index fields
  INTEGER, DIMENSION(:,:),   POINTER :: iidx, iblk
  INTEGER, DIMENSION(:,:,:), POINTER :: ichcidx, ichcblk

!-----------------------------------------------------------------------

  ! Check if gradient limiting is required
  IF (PRESENT(llimit_nneg)) THEN
    l_limit_nneg = llimit_nneg
  ELSE
    l_limit_nneg = .FALSE.
  ENDIF

  ! Compute values for dynamic nproma blocking
  nproma_ubcintp = MIN(nproma,256)
  nblks_ubcintp  = INT(p_grf%npoints_ubcintp_c/nproma_ubcintp)
  npromz_ubcintp = MOD(p_grf%npoints_ubcintp_c,nproma_ubcintp)
  IF (npromz_ubcintp > 0) THEN
    nblks_ubcintp = nblks_ubcintp + 1
  ELSE
    npromz_ubcintp = nproma_ubcintp
  ENDIF

  epsi = 1.e-75_wp
  ovsht_fac = 1.0_wp ! factor of allowed overshooting
  r_ovsht_fac = 1._wp/ovsht_fac
 
  ! Pointers to index lists for gradient computation
  iidx => p_grf%idxlist_ubcintp_c
  iblk => p_grf%blklist_ubcintp_c

  ! child cell indices and blocks for non-MPI parent-to-child communication
  ichcidx => p_pp%cells%child_idx
  ichcblk => p_pp%cells%child_blk


!$OMP PARALLEL

!$OMP DO PRIVATE (jb,nlen,nshift,jn,jc,limfac1,limfac2,limfac, &
!$OMP   maxval_neighb,minval_neighb,relaxed_minval,relaxed_maxval) ICON_OMP_DEFAULT_SCHEDULE
  DO jb = 1, nblks_ubcintp
    IF (jb == nblks_ubcintp) THEN
      nlen = npromz_ubcintp
    ELSE
      nlen = nproma_ubcintp
    ENDIF
    nshift = (jb-1)*nproma_ubcintp

#ifdef __LOOP_EXCHANGE
    DO jc = nshift+1, nshift+nlen
      DO jn = 1, nfields
#else
!CDIR NOLOOPCHG
    DO jn = 1, nfields
      DO jc = nshift+1, nshift+nlen
#endif

        val_ctr(jn,jc) = f3din(iidx(1,jc),jn,iblk(1,jc))
        grad_x(jn,jc) =  &
          p_grf%coeff_ubcintp_c(1,1,jc)*f3din(iidx(1,jc),jn,iblk(1,jc)) + &
          p_grf%coeff_ubcintp_c(2,1,jc)*f3din(iidx(2,jc),jn,iblk(2,jc)) + &
          p_grf%coeff_ubcintp_c(3,1,jc)*f3din(iidx(3,jc),jn,iblk(3,jc)) + &
          p_grf%coeff_ubcintp_c(4,1,jc)*f3din(iidx(4,jc),jn,iblk(4,jc)) + &
          p_grf%coeff_ubcintp_c(5,1,jc)*f3din(iidx(5,jc),jn,iblk(5,jc)) + &
          p_grf%coeff_ubcintp_c(6,1,jc)*f3din(iidx(6,jc),jn,iblk(6,jc)) + &
          p_grf%coeff_ubcintp_c(7,1,jc)*f3din(iidx(7,jc),jn,iblk(7,jc)) + &
          p_grf%coeff_ubcintp_c(8,1,jc)*f3din(iidx(8,jc),jn,iblk(8,jc)) + &
          p_grf%coeff_ubcintp_c(9,1,jc)*f3din(iidx(9,jc),jn,iblk(9,jc)) + &
          p_grf%coeff_ubcintp_c(10,1,jc)*f3din(iidx(10,jc),jn,iblk(10,jc))
        grad_y(jn,jc) =  &
          p_grf%coeff_ubcintp_c(1,2,jc)*f3din(iidx(1,jc),jn,iblk(1,jc)) + &
          p_grf%coeff_ubcintp_c(2,2,jc)*f3din(iidx(2,jc),jn,iblk(2,jc)) + &
          p_grf%coeff_ubcintp_c(3,2,jc)*f3din(iidx(3,jc),jn,iblk(3,jc)) + &
          p_grf%coeff_ubcintp_c(4,2,jc)*f3din(iidx(4,jc),jn,iblk(4,jc)) + &
          p_grf%coeff_ubcintp_c(5,2,jc)*f3din(iidx(5,jc),jn,iblk(5,jc)) + &
          p_grf%coeff_ubcintp_c(6,2,jc)*f3din(iidx(6,jc),jn,iblk(6,jc)) + &
          p_grf%coeff_ubcintp_c(7,2,jc)*f3din(iidx(7,jc),jn,iblk(7,jc)) + &
          p_grf%coeff_ubcintp_c(8,2,jc)*f3din(iidx(8,jc),jn,iblk(8,jc)) + &
          p_grf%coeff_ubcintp_c(9,2,jc)*f3din(iidx(9,jc),jn,iblk(9,jc)) + &
          p_grf%coeff_ubcintp_c(10,2,jc)*f3din(iidx(10,jc),jn,iblk(10,jc))
        maxval_neighb(jn,jc) =                 &
          MAX(f3din(iidx(1,jc),jn,iblk(1,jc)), &
          f3din(iidx(2,jc),jn,iblk(2,jc)),     &
          f3din(iidx(3,jc),jn,iblk(3,jc)),     &
          f3din(iidx(4,jc),jn,iblk(4,jc)),     &
          f3din(iidx(5,jc),jn,iblk(5,jc)),     &
          f3din(iidx(6,jc),jn,iblk(6,jc)),     &
          f3din(iidx(7,jc),jn,iblk(7,jc)),     &
          f3din(iidx(8,jc),jn,iblk(8,jc)),     &
          f3din(iidx(9,jc),jn,iblk(9,jc)),     &
          f3din(iidx(10,jc),jn,iblk(10,jc)))
        minval_neighb(jn,jc) =                 &
          MIN(f3din(iidx(1,jc),jn,iblk(1,jc)), &
          f3din(iidx(2,jc),jn,iblk(2,jc)),     &
          f3din(iidx(3,jc),jn,iblk(3,jc)),     &
          f3din(iidx(4,jc),jn,iblk(4,jc)),     &
          f3din(iidx(5,jc),jn,iblk(5,jc)),     &
          f3din(iidx(6,jc),jn,iblk(6,jc)),     &
          f3din(iidx(7,jc),jn,iblk(7,jc)),     &
          f3din(iidx(8,jc),jn,iblk(8,jc)),     &
          f3din(iidx(9,jc),jn,iblk(9,jc)),     &
          f3din(iidx(10,jc),jn,iblk(10,jc)))
      ENDDO
    ENDDO

#ifdef __LOOP_EXCHANGE
    DO jc = nshift+1, nshift+nlen
      DO jn = 1, nfields
#else
!CDIR NOLOOPCHG
    DO jn = 1, nfields
      DO jc = nshift+1, nshift+nlen
#endif
        min_expval = MIN(grad_x(jn,jc)*p_grf%dist_pc2cc_ubc(1,1,jc) + &
                         grad_y(jn,jc)*p_grf%dist_pc2cc_ubc(1,2,jc),  &
                         grad_x(jn,jc)*p_grf%dist_pc2cc_ubc(2,1,jc) + &
                         grad_y(jn,jc)*p_grf%dist_pc2cc_ubc(2,2,jc),  &
                         grad_x(jn,jc)*p_grf%dist_pc2cc_ubc(3,1,jc) + &
                         grad_y(jn,jc)*p_grf%dist_pc2cc_ubc(3,2,jc),  &
                         grad_x(jn,jc)*p_grf%dist_pc2cc_ubc(4,1,jc) + &
                         grad_y(jn,jc)*p_grf%dist_pc2cc_ubc(4,2,jc),  &
                         -1.e-80_wp )
        max_expval = MAX(grad_x(jn,jc)*p_grf%dist_pc2cc_ubc(1,1,jc) + &
                         grad_y(jn,jc)*p_grf%dist_pc2cc_ubc(1,2,jc),  &
                         grad_x(jn,jc)*p_grf%dist_pc2cc_ubc(2,1,jc) + &
                         grad_y(jn,jc)*p_grf%dist_pc2cc_ubc(2,2,jc),  &
                         grad_x(jn,jc)*p_grf%dist_pc2cc_ubc(3,1,jc) + &
                         grad_y(jn,jc)*p_grf%dist_pc2cc_ubc(3,2,jc),  &
                         grad_x(jn,jc)*p_grf%dist_pc2cc_ubc(4,1,jc) + &
                         grad_y(jn,jc)*p_grf%dist_pc2cc_ubc(4,2,jc),  &
                         1.e-80_wp )

        limfac1 = 1._wp
        limfac2 = 1._wp
        ! Allow a limited amount of over-/undershooting in the downscaled fields
        IF (minval_neighb(jn,jc) > 0._wp) THEN
          relaxed_minval = r_ovsht_fac*minval_neighb(jn,jc)
        ELSE
          relaxed_minval = ovsht_fac*minval_neighb(jn,jc)
        ENDIF
        IF (maxval_neighb(jn,jc) > 0._wp) THEN
          relaxed_maxval = ovsht_fac*maxval_neighb(jn,jc)
        ELSE
          relaxed_maxval = r_ovsht_fac*maxval_neighb(jn,jc)
        ENDIF

        IF (val_ctr(jn,jc) + min_expval < relaxed_minval-epsi) THEN
          limfac1 = ABS((relaxed_minval-val_ctr(jn,jc))/min_expval)
        ENDIF
        IF (val_ctr(jn,jc) + max_expval > relaxed_maxval+epsi) THEN
          limfac2 = ABS((relaxed_maxval-val_ctr(jn,jc))/max_expval)
        ENDIF
        limfac = MIN(limfac1,limfac2)

        grad_x(jn,jc) = grad_x(jn,jc)*limfac
        grad_y(jn,jc) = grad_y(jn,jc)*limfac

      ENDDO
    ENDDO

    IF (l_limit_nneg) THEN
#ifdef __LOOP_EXCHANGE
      DO jc = nshift+1, nshift+nlen
        DO jn = 1, nfields
#else
!CDIR NOLOOPCHG
      DO jn = 1, nfields
        DO jc = nshift+1, nshift+nlen
#endif

          h_aux(jn,jc,1) = MAX(0._wp, val_ctr(jn,jc)    + &
            grad_x(jn,jc)*p_grf%dist_pc2cc_ubc(1,1,jc)  + &
            grad_y(jn,jc)*p_grf%dist_pc2cc_ubc(1,2,jc))
          h_aux(jn,jc,2) = MAX(0._wp, val_ctr(jn,jc)    + &
            grad_x(jn,jc)*p_grf%dist_pc2cc_ubc(2,1,jc)  + &
            grad_y(jn,jc)*p_grf%dist_pc2cc_ubc(2,2,jc))
          h_aux(jn,jc,3) = MAX(0._wp, val_ctr(jn,jc)    + &
            grad_x(jn,jc)*p_grf%dist_pc2cc_ubc(3,1,jc)  + &
            grad_y(jn,jc)*p_grf%dist_pc2cc_ubc(3,2,jc))
          h_aux(jn,jc,4) = MAX(0._wp, val_ctr(jn,jc)    + &
            grad_x(jn,jc)*p_grf%dist_pc2cc_ubc(4,1,jc)  + &
            grad_y(jn,jc)*p_grf%dist_pc2cc_ubc(4,2,jc))

        ENDDO
      ENDDO
        ELSE
#ifdef __LOOP_EXCHANGE
      DO jc = nshift+1, nshift+nlen
        DO jn = 1, nfields
#else
!CDIR NOLOOPCHG
      DO jn = 1, nfields
        DO jc = nshift+1, nshift+nlen
#endif

          h_aux(jn,jc,1) = val_ctr(jn,jc)              + &
            grad_x(jn,jc)*p_grf%dist_pc2cc_ubc(1,1,jc) + &
            grad_y(jn,jc)*p_grf%dist_pc2cc_ubc(1,2,jc)
          h_aux(jn,jc,2) = val_ctr(jn,jc)              + &
            grad_x(jn,jc)*p_grf%dist_pc2cc_ubc(2,1,jc) + &
            grad_y(jn,jc)*p_grf%dist_pc2cc_ubc(2,2,jc)
          h_aux(jn,jc,3) = val_ctr(jn,jc)              + &
            grad_x(jn,jc)*p_grf%dist_pc2cc_ubc(3,1,jc) + &
            grad_y(jn,jc)*p_grf%dist_pc2cc_ubc(3,2,jc)
          h_aux(jn,jc,4) = val_ctr(jn,jc)              + &
            grad_x(jn,jc)*p_grf%dist_pc2cc_ubc(4,1,jc) + &
            grad_y(jn,jc)*p_grf%dist_pc2cc_ubc(4,2,jc)

        ENDDO
      ENDDO
    ENDIF

  ENDDO ! blocks
!$OMP END DO


  ! Store results in p_out

  IF (my_process_is_mpi_seq()) THEN

!$OMP DO PRIVATE (jb,nlen,nshift,jn,jc,ic,ib) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = 1, nblks_ubcintp
      IF (jb == nblks_ubcintp) THEN
        nlen = npromz_ubcintp
      ELSE
        nlen = nproma_ubcintp
      ENDIF
      nshift = (jb-1)*nproma_ubcintp

#ifdef __LOOP_EXCHANGE
      DO jc = nshift+1, nshift+nlen
        ic = iidx(1,jc)
        ib = iblk(1,jc)
        DO jn = 1, nfields
#else
!CDIR NOLOOPCHG
      DO jn = 1, nfields
!CDIR NODEP
        DO jc = nshift+1, nshift+nlen
          ic = iidx(1,jc)
          ib = iblk(1,jc)
#endif

          f3dout(ichcidx(ic,ib,1),jn,ichcblk(ic,ib,1)) = h_aux(jn,jc,1)
          f3dout(ichcidx(ic,ib,2),jn,ichcblk(ic,ib,2)) = h_aux(jn,jc,2)
          f3dout(ichcidx(ic,ib,3),jn,ichcblk(ic,ib,3)) = h_aux(jn,jc,3)
          f3dout(ichcidx(ic,ib,4),jn,ichcblk(ic,ib,4)) = h_aux(jn,jc,4)

        ENDDO
      ENDDO
    ENDDO
!$OMP END DO NOWAIT

  ENDIF
!$OMP END PARALLEL

IF (my_process_is_mpi_parallel()) THEN

  nsendtot = SUM(p_pc%comm_pat_interpol_scal_ubc(1:4)%n_send)
  nrecvtot = SUM(p_pc%comm_pat_interpol_scal_ubc(1:4)%n_recv)


  CALL exchange_data_grf(p_pc%comm_pat_interpol_scal_ubc,1,nfields,nsendtot, &
                         nrecvtot,RECV1=f3dout,SEND1=h_aux)

ENDIF

END SUBROUTINE interpol_scal_ubc


!-------------------------------------------------------------------------
END MODULE mo_grf_ubcintp

