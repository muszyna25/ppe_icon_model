!>
!! @brief Module to collect routines that write radiation diagnostics on
!!        output variables
!!
!! @remarks
!!   collect routines here that write variables of the radiation part onto
!!   output variables contained in mo_echam_phy_memory.f90
!!
!!
!! @author Sebastian Rast, MPI-M, Hamburg (2015-06-19): Original source
!!
!! $ID: n/a$
!!
!! @par Origin
!! Old rad_aero_diag of mo_radiation relocated here because it is needed in
!! in the psrad radiation also.
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!
MODULE mo_rad_diag

  USE mo_kind,                 ONLY: wp
  USE mo_echam_phy_memory,     ONLY: prm_field

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: rad_aero_diag

CONTAINS

!>
!! SUBROUTINE rad_aero_diag writes actual aerosol optical properties to output stream
!!
!! @author J.S.Rast, MPI-Met Hamburg
!!
!! @par Revision History
!! Origianl Source by J.S.Rast, MPI-Met Hamburg, (2013-08-30)
!!
!!
SUBROUTINE rad_aero_diag (                kg              , &
      & kb              ,kcs             ,kce             , &
      & kbdim           ,klev            ,kpband          , &
      & kpsw            ,paer_tau_lw_vr  ,paer_tau_sw_vr  , &
      & paer_piz_sw_vr  ,paer_cg_sw_vr                      )

      INTEGER, INTENT(in)    :: kg      ! domain index
      INTEGER, INTENT(in)    :: kb      ! block index
      INTEGER, INTENT(in)    :: kcs     ! actual block length (start index)
      INTEGER, INTENT(in)    :: kce     ! actual block length (end index)
      INTEGER, INTENT(in)    :: kbdim   ! declaration block length
      INTEGER, INTENT(in)    :: klev    ! levels
      INTEGER, INTENT(in)    :: kpband  ! number of lw bands
      INTEGER, INTENT(in)    :: kpsw    ! number of sw bands
      REAL(wp), INTENT(in)   :: paer_tau_lw_vr(kbdim,klev,kpband) ! aod thermal wavelengths
      REAL(wp), INTENT(in)   :: paer_tau_sw_vr(kbdim,klev,kpsw)   ! aod solar wavelengths
      REAL(wp), INTENT(in)   :: paer_piz_sw_vr(kbdim,klev,kpsw)   ! ssa solar wavelengths
      REAL(wp), INTENT(in)   :: paer_cg_sw_vr(kbdim,klev,kpsw)    ! asy solar wavelengths

      prm_field(kg)%aer_aod_9731(kcs:kce,1:klev,kb) = &
                   paer_tau_lw_vr(kcs:kce,klev:1:-1,7)
      prm_field(kg)%aer_aod_533 (kcs:kce,1:klev,kb) = &
                   paer_tau_sw_vr(kcs:kce,klev:1:-1,10)
      prm_field(kg)%aer_ssa_533 (kcs:kce,1:klev,kb) = &
                   paer_piz_sw_vr(kcs:kce,klev:1:-1,10)
      prm_field(kg)%aer_asy_533 (kcs:kce,1:klev,kb) = &
                   paer_cg_sw_vr(kcs:kce,klev:1:-1,10)
      prm_field(kg)%aer_aod_2325(kcs:kce,1:klev,kb) = &
                   paer_tau_sw_vr(kcs:kce,klev:1:-1,3)
      prm_field(kg)%aer_ssa_2325(kcs:kce,1:klev,kb) = &
                   paer_piz_sw_vr(kcs:kce,klev:1:-1,3)
      prm_field(kg)%aer_asy_2325(kcs:kce,1:klev,kb) = &
                   paer_cg_sw_vr(kcs:kce,klev:1:-1,3)
END SUBROUTINE rad_aero_diag

END MODULE mo_rad_diag
