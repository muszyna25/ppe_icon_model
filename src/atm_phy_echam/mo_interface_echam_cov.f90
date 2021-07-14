!>
!! @brief Subroutine interface_echam_cov calls the cloud colver scheme.
!!
!! @author Hui Wan, MPI-M
!! @author Marco Giorgetta, MPI-M
!!
!! @par Revision History
!!  Original version from ECHAM6 (revision 2028)
!!  Modified for ICOHAM by Hui Wan and Marco Giorgetta (2010)
!!  Modified for ICONAM by Marco Giorgetta (2014)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!

MODULE mo_interface_echam_cov

  USE mo_kind                ,ONLY: wp

  USE mo_echam_phy_memory    ,ONLY: t_echam_phy_field, prm_field
  USE mo_echam_sfc_indices   ,ONLY: nsfc_type, iwtr, iice
  
  USE mo_timer               ,ONLY: ltimer, timer_start, timer_stop, timer_cov

  USE mo_run_config          ,ONLY: iqv, iqc, iqi
  USE mo_cover               ,ONLY: cover
  
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: interface_echam_cov

CONTAINS

  SUBROUTINE interface_echam_cov(jg,jb,jcs,jce ,&
       &                         nproma,nlev   ) 

    ! Arguments
    !
    INTEGER                 ,INTENT(in) :: jg,jb,jcs,jce
    INTEGER                 ,INTENT(in) :: nproma,nlev

    ! Pointers
    !
    TYPE(t_echam_phy_field) ,POINTER    :: field

    ! Local variables
    !
    INTEGER  :: nlevp1, jc
    INTEGER  :: itype(nproma) !< type of convection
    REAL(wp) :: zfrw (nproma) !< cell area fraction of open water
    REAL(wp) :: zfri (nproma) !< cell area fraction of ice covered water

    IF (ltimer) call timer_start(timer_cov)

    field => prm_field(jg)

    nlevp1 = nlev+1

    !$ACC DATA PRESENT( field% rtype, field%frac_tile ) &
    !$ACC       CREATE( itype, zfrw, zfri )

    !$ACC PARALLEL DEFAULT(PRESENT)
    !$ACC LOOP GANG VECTOR
    DO jc = jcs, jce
      itype(jc) = NINT(field%rtype(jc,jb))
    END DO
    !$ACC END PARALLEL

    !$ACC PARALLEL DEFAULT(PRESENT)
    !$ACC LOOP GANG VECTOR
    DO jc = 1, nproma
      IF (iwtr.LE.nsfc_type) THEN
         zfrw(jc) = field%frac_tile(jc,jb,iwtr)
      ELSE
         zfrw(jc) = 0.0_wp
      END IF

      IF (iice.LE.nsfc_type) THEN
         zfri(jc) = field%frac_tile(jc,jb,iice)
      ELSE
         zfri(jc) = 0.0_wp
      END IF
    END DO
    !$ACC END PARALLEL

    CALL cover(    jg,                        &! in
         &         jb,                        &! in
         &         jcs, jce, nproma,          &! in
         &         nlev, nlevp1,              &! in
         &         itype,                     &! in
         &         zfrw(:),                   &! in
         &         zfri(:),                   &! in
         &         field% zf(:,:,jb),         &! in
         &         field% phalf(:,:,jb),      &! in
         &         field% pfull(:,:,jb),      &! in
         &         field%  ta(:,:,jb),        &! in    tm1
         &         field%  qtrc(:,:,jb,iqv),  &! in    qm1
         &         field%  qtrc(:,:,jb,iqc),  &! in    xlm1
         &         field%  qtrc(:,:,jb,iqi),  &! in    xim1
         &         field%  aclc(:,:,jb),      &! out   (for "radiation" and "vdiff_down")
         &         field% rintop(:,  jb)     ) ! out   (for output)

    !$ACC END DATA

    NULLIFY(field)

    IF (ltimer) call timer_stop(timer_cov)

  END SUBROUTINE interface_echam_cov

END MODULE mo_interface_echam_cov
