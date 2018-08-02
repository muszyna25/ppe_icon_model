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
USE mo_model_domain,        ONLY: t_patch
USE mo_parallel_config,     ONLY: nproma
USE mo_communication,       ONLY: exchange_data_grf

USE mo_grf_intp_data_strc


IMPLICIT NONE

PRIVATE

PUBLIC :: interpol_vec_ubc, interpol_scal_ubc

CONTAINS

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
  TYPE(t_patch), INTENT(in) :: p_pp
  TYPE(t_patch), INTENT(inout) :: p_pc

  ! input: delta vn at interface level
  REAL(wp), INTENT(IN) :: p_vn_in(:,:) ! dim: (nproma,nblks_e)

  ! Indices of source points and interpolation coefficients
  TYPE(t_gridref_single_state),   TARGET, INTENT(IN) ::  p_grf

  ! reconstructed edge-normal wind component
  REAL(wp),INTENT(INOUT) :: p_vn_out(nproma,1,p_pc%nblks_e)


  INTEGER :: jb, je                    ! loop indices
  INTEGER :: nproma_ubcintp, nblks_ubcintp, npromz_ubcintp, nlen, nshift

  REAL(wp) :: vn_aux(1,p_grf%npoints_ubcintp_e,4)

  ! Pointers to index fields/lists
  INTEGER,  DIMENSION(:,:),   POINTER :: iidx, iblk

  !-----------------------------------------------------------------------
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
    nlen = MERGE(nproma_ubcintp, npromz_ubcintp, jb /= nblks_ubcintp)
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
!$OMP END PARALLEL

  ! Store results in p_vn_out

  CALL exchange_data_grf(p_pc%comm_pat_coll_interpol_vec_ubc,1,1, &
    &                    RECV1=p_vn_out,SEND1=vn_aux)

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
  TYPE(t_patch), INTENT(in) :: p_pp
  TYPE(t_patch), INTENT(inout) :: p_pc

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

  LOGICAL :: l_limit_nneg     ! local variable corresponding to llimit_nneg

  ! Auxiliary fields
  REAL(wp), DIMENSION(nfields,p_grf%npoints_ubcintp_c) :: &
    grad_x, grad_y, maxval_neighb, minval_neighb, val_ctr
  REAL(wp) :: h_aux(nfields,p_grf%npoints_ubcintp_c,4)
  REAL(wp) :: limfac1, limfac2, limfac, min_expval, max_expval, epsi, ovsht_fac, r_ovsht_fac, &
              relaxed_minval, relaxed_maxval

  ! Pointers to index fields
  INTEGER, DIMENSION(:,:),   POINTER :: iidx, iblk

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

!$OMP PARALLEL

!$OMP DO PRIVATE (jb,nlen,nshift,jn,jc,limfac1,limfac2,limfac,min_expval,max_expval, &
!$OMP   relaxed_minval,relaxed_maxval) ICON_OMP_DEFAULT_SCHEDULE
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
!$OMP END PARALLEL

  ! Store results in p_out

  CALL exchange_data_grf(p_pc%comm_pat_coll_interpol_scal_ubc,1,nfields, &
    &                    RECV1=f3dout,SEND1=h_aux)

END SUBROUTINE interpol_scal_ubc


!-------------------------------------------------------------------------
END MODULE mo_grf_ubcintp

