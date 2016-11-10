!>
!!  This module contains parameters, subroutines and
!!  functions used in the Held-Suarez test of the 3D
!!  hydrostatic dynamical core.
!!
!! @par See also
!! Held, I. M., and M. J. Suarez, 1994: A proposal for the intercomparison
!!   of the dynamical cores of atmospheric general circulation models.
!!   Bulletin of the Amer. Meteor. Soc., vol. 75, pp. 1825-1830.
!! Wan, H., M. A. Giorgetta and L. Bonaventura, 2008: Ensemble Held-Suarez
!!   test with a spectral transform model: variability, sensitivity, and
!!   convergence. Mon. Wea. Rev., Vol. 136, 1075-1992.
!!
!! @par Revision History
!!  Original implementation in ECHAM5 by Hui Wan, MPI-M (2005-07)
!!  Adaptation for the ICOHDC by Hui Wan, MPI-M (2008-05-30)
!!  Further adaptation for the restructured code by Hui Wan, MPI-M (2009-02-03)
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
MODULE mo_hs_test

  USE mo_kind,                ONLY: wp
  USE mo_physical_constants,  ONLY: rd, cpd, cvd
  USE mo_physical_constants,  ONLY: rdaylen

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: held_suarez_forcing_vn
  PUBLIC :: held_suarez_forcing_temp

  ! Module parameters

  REAL(wp), PARAMETER :: HSsigb    = 0.7_wp
  REAL(wp), PARAMETER :: HSvcoeff1 = 1._wp/(1._wp-HSsigb)
  REAL(wp), PARAMETER :: HSvcoeff2 = -HSsigb/(1._wp-HSsigb)

  REAL(wp), PARAMETER :: HSkf    = 1._wp/ rdaylen
  REAL(wp), PARAMETER :: HSka    = 1._wp/(rdaylen*40._wp)
  REAL(wp), PARAMETER :: HSks    = 1._wp/(rdaylen* 4._wp )

  REAL(wp), PARAMETER :: HSdty   = 60._wp
  REAL(wp), PARAMETER :: HSdthz  = 10._wp

  REAL(wp), PARAMETER :: HScappa = rd/cpd
  REAL(wp), PARAMETER :: HSp0    = 1.E5_wp
  REAL(wp), PARAMETER :: HSt0    = 200._wp
  REAL(wp), PARAMETER :: HSt1    = 315._wp

CONTAINS
  !>
  !! Linear damping of velocity in the 'boundary layer'.
  !!
  !! @par Revision History
  !! Original version for ECHAM5 by Hui Wan, MPI-M (2005-07)
  !! Adaptation for the ICOHDC by Hui Wan, MPI-M (2008-05-30):
  !! Further adaptation for the restructured code by Hui Wan, MPI-M (2009-02-03)
  !!
  SUBROUTINE held_suarez_forcing_vn( pvn, psigma,           & !in
                                   & nlev, nproma, is, ie,  & !in
                                   & fvn_hs )                 !out

    INTEGER, INTENT(IN) :: nlev            !< number of vertical layers
    INTEGER, INTENT(IN) :: nproma, is, ie  !<

    REAL(wp),INTENT(IN) :: pvn      ( nproma, nlev ) !< normal velocity
    REAL(wp),INTENT(IN) :: psigma   ( nproma, nlev ) !< sigma = pres/pres_sfc

    REAL(wp),INTENT(INOUT) :: fvn_hs( nproma, nlev ) !< forcing on velocity

    REAL(wp) :: ztmp( nproma, nlev )
    !---

    ztmp  (is:ie,:) = psigma(is:ie,:)*HSvcoeff1 + HSvcoeff2
    fvn_hs(is:ie,:) = -HSkf*MAX( 0._wp,ztmp(is:ie,:) )*pvn(is:ie,:)

  END SUBROUTINE held_suarez_forcing_vn
  !-------------
  !>
  !! Newtonian cooling.
  !!
  !! @par Revision History
  !! Original version for ECHAM5 by Hui Wan, MPI-M (2005-07)
  !! Adaptation for the ICOHDC by Hui Wan, MPI-M (2008-05-30):
  !! Further adaptation for the restructured code by Hui Wan (2009-02-03)
  !!
  SUBROUTINE held_suarez_forcing_temp( ptemp_mc, ppres_mc,   &! in
                                     & psigma, plat,         &! in
                                     & nlev, nproma, is, ie, &! in
                                     & fT_hs,                &! out
                                     & opt_ekinh, opt_ldissip_heat)

    INTEGER, INTENT(IN) :: nlev                   !< number of vertical layers
    INTEGER, INTENT(IN) :: nproma, is, ie
    REAL(wp),INTENT(IN) :: ptemp_mc (nproma,nlev) !< temperature in Kelvin
    REAL(wp),INTENT(IN) :: ppres_mc (nproma,nlev) !< pressure in Pa
    REAL(wp),INTENT(IN) :: psigma   (nproma,nlev) !< sigma = pres/pres_sfc
    REAL(wp),INTENT(IN) :: plat     (nproma)      !< latitide in radians

    REAL(wp),INTENT(IN), OPTIONAL :: opt_ekinh(nproma,nlev) !< kinetic energy
    LOGICAL, OPTIONAL  :: opt_ldissip_heat        !< dissipative heating or not

    REAL(wp),INTENT(OUT):: fT_hs (nproma,nlev) !< forcing on temperature

    INTEGER :: jk  !vertical layer index

    REAL(wp) :: zsinlat2(nproma), zcoslat2(nproma), zcoslat4(nproma)
    REAL(wp) :: zsigma0(nproma), zTempEq(nproma), ztmp(nproma), kT_hs(nproma)
    LOGICAL  :: l_friheat

    !------------------------------

    ! check, whether dissipative heating should be performed
    IF ( PRESENT(opt_ldissip_heat) ) THEN
      l_friheat=opt_ldissip_heat
    ELSE
      l_friheat=.FALSE.
    ENDIF

    ! latitude related parameters

    zsinlat2(is:ie) = SIN(plat(is:ie))**2
    zcoslat2(is:ie) = 1._wp - zsinlat2(is:ie)
    zcoslat4(is:ie) = zcoslat2(is:ie)**2

    DO jk = 1,nlev

       ! equilibrium temperature

       zsigma0(is:ie) = ppres_mc(is:ie,jk)/HSp0
       zTempEq(is:ie) = HSt1 - HSdty*zsinlat2(is:ie)               &
                      & - HSdthz*LOG(zsigma0(is:ie))*zcoslat2(is:ie)
       zTempEq(is:ie) = zTempEq(is:ie) * zsigma0(is:ie)**HScappa
       zTempEq(is:ie) = MAX( 200._wp, zTempEq(is:ie) )

       ! pressure-dependent coefficient

       ztmp (is:ie)   = psigma(is:ie,jk)*HSvcoeff1 + HSvcoeff2
       kT_hs(is:ie)   = HSka + (HSks-HSka)                       &
                      & *zcoslat4(is:ie)*MAX(0._wp,ztmp(is:ie))

       ! Newtonian cooling

       fT_hs(is:ie,jk) = -kT_hs(is:ie)                        &
                       & *( ptemp_mc(is:ie,jk)-zTempEq(is:ie) )

    ENDDO !vertical layer loop

    IF (l_friheat) THEN
      DO jk = 1,nlev
        ztmp (is:ie) = psigma(is:ie,jk)*HSvcoeff1 + HSvcoeff2
        fT_hs(is:ie,jk) = fT_hs(is:ie,jk) &
        & + HSkf*MAX( 0._wp,ztmp(is:ie))*2.0_wp*opt_ekinh(is:ie,jk)/cvd
      ENDDO
    ENDIF

  END SUBROUTINE held_suarez_forcing_temp
  !-------------

END MODULE mo_hs_test
