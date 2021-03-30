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

    REAL(wp) :: ztmp
    INTEGER  :: i, jk
    !---

    !$ACC DATA PRESENT( pvn, psigma, fvn_hs )
    !$ACC PARALLEL
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO jk=1,nlev
      DO i = is, ie
        ztmp = psigma(i,jk)*HSvcoeff1 + HSvcoeff2
        fvn_hs(i,jk) = -HSkf*MAX( 0._wp,ztmp )*pvn(i,jk)
      ENDDO
    ENDDO
    !$ACC END PARALLEL
    !$ACC END DATA
   

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
    INTEGER :: i   ! nproma index
    REAL(wp) :: zsinlat2(nproma), zcoslat2(nproma), zcoslat4(nproma)
    REAL(wp) :: zsigma0, zTempEq, ztmp, kT_hs
    LOGICAL  :: l_friheat

    !------------------------------

    ! check, whether dissipative heating should be performed
    IF ( PRESENT(opt_ldissip_heat) ) THEN
      l_friheat=opt_ldissip_heat
    ELSE
      l_friheat=.FALSE.
    ENDIF

    ! latitude related parameters

    !$ACC DATA CREATE( zsinlat2, zcoslat2, zcoslat4 ) PRESENT( ptemp_mc, ppres_mc, psigma, plat, fT_hs )
    !$ACC DATA PRESENT( opt_ekinh ) IF( PRESENT( opt_ekinh ) )
    !$ACC PARALLEL
    !$ACC LOOP GANG VECTOR
    DO i=is, ie
      zsinlat2(i) = SIN(plat(i))**2
      zcoslat2(i) = 1._wp - zsinlat2(i)
      zcoslat4(i) = zcoslat2(i)**2
    ENDDO
    !$ACC END PARALLEL
   
    !$ACC PARALLEL
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO jk = 1,nlev

      DO i=is,ie
        ! equilibrium temperature

        zsigma0 = ppres_mc(i,jk)/HSp0
        zTempEq = HSt1 - HSdty*zsinlat2(i)               &
                       & - HSdthz*LOG(zsigma0)*zcoslat2(i)
        zTempEq = zTempEq * zsigma0**HScappa
        zTempEq = MAX( 200._wp, zTempEq )

        ! pressure-dependent coefficient

        ztmp   = psigma(i,jk)*HSvcoeff1 + HSvcoeff2
        kT_hs  = HSka + (HSks-HSka)*zcoslat4(i)*MAX(0._wp,ztmp)

        ! Newtonian cooling

        fT_hs(i,jk) = -kT_hs *( ptemp_mc(i,jk)-zTempEq )
      ENDDO

    ENDDO !vertical layer loop
    !$ACC END PARALLEL

    IF (l_friheat) THEN
      !$ACC PARALLEL
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO jk = 1,nlev
        DO i = is,ie
          ztmp = psigma(i,jk)*HSvcoeff1 + HSvcoeff2
          fT_hs(i,jk) = fT_hs(i,jk) + HSkf*MAX( 0._wp,ztmp)*2.0_wp*opt_ekinh(i,jk)/cvd
        ENDDO
     ENDDO
     !$ACC END PARALLEL
  ENDIF

  !$ACC END DATA
  !$ACC END DATA

  END SUBROUTINE held_suarez_forcing_temp
  !-------------

END MODULE mo_hs_test
