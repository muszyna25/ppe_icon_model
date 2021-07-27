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
SUBROUTINE rad_aero_diag (                                  &
      & kcs             ,kce             ,kbdim           , &
      & klev            ,kpband          ,kpsw            , &
      & paer_tau_lw_vr  ,paer_tau_sw_vr  ,paer_piz_sw_vr  , &
      & paer_cg_sw_vr                                     , &
      & aer_aod_533     ,aer_ssa_533     ,aer_asy_533     , &
      & aer_aod_2325    ,aer_ssa_2325    ,aer_asy_2325    , &
      & aer_aod_9731                                      , &
      & opt_use_acc                                       )

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
      REAL(wp), INTENT(inout)  :: &
           & aer_aod_533 (kbdim,klev), & ! aod at 533 nm      
           & aer_ssa_533 (kbdim,klev), & ! ssa at 533 nm
           & aer_asy_533 (kbdim,klev), & ! asy at 533 nm
           & aer_aod_2325(kbdim,klev), & ! aod at 2325 nm      
           & aer_ssa_2325(kbdim,klev), & ! ssa at 2325 nm
           & aer_asy_2325(kbdim,klev), & ! asy at 2325 nm
           & aer_aod_9731(kbdim,klev)    ! aod at 9731 nm
      LOGICAL, INTENT(IN), OPTIONAL  :: opt_use_acc

      LOGICAL :: use_acc   = .FALSE.  ! Default: no acceleration
      IF (PRESENT(opt_use_acc)) use_acc = opt_use_acc

      !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF (use_acc)
      aer_aod_9731(kcs:kce,1:klev) = &
           paer_tau_lw_vr(kcs:kce,klev:1:-1,7)
      aer_aod_533 (kcs:kce,1:klev) = &
           paer_tau_sw_vr(kcs:kce,klev:1:-1,10)
      aer_ssa_533 (kcs:kce,1:klev) = &
           paer_piz_sw_vr(kcs:kce,klev:1:-1,10)
      aer_asy_533 (kcs:kce,1:klev) = &
           paer_cg_sw_vr(kcs:kce,klev:1:-1,10)
      aer_aod_2325(kcs:kce,1:klev) = &
           paer_tau_sw_vr(kcs:kce,klev:1:-1,3)
      aer_ssa_2325(kcs:kce,1:klev) = &
           paer_piz_sw_vr(kcs:kce,klev:1:-1,3)
      aer_asy_2325(kcs:kce,1:klev) = &
           paer_cg_sw_vr(kcs:kce,klev:1:-1,3)
      aer_aod_9731(kcs:kce,1:klev) = &
           paer_tau_lw_vr(kcs:kce,klev:1:-1,7)
      !$ACC END KERNELS

END SUBROUTINE rad_aero_diag

END MODULE mo_rad_diag
