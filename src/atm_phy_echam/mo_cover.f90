#ifdef __xlC__
@PROCESS HOT
@PROCESS XLF90(NOSIGNEDZERO)
#endif
#include "fsel.inc"

!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
!! @brief Module diagnoses cloud cover for current timestep
!!
!! @remarks
!!     In ECHAM4 all cloud cover calculations were performed
!!     in routine CLOUD.  The routine diagnosed the cloud cover
!!     from the previous timestep using m1 variables, but then
!!     used a "first guess" calculation of T and Q to give a
!!     cloud cover estimate for the next timestep.  This meant
!!     that the radiation scheme used cloud cover values that were
!!     from a different timestep to the temperature and water vapour
!!     values, and also that the T and Q values were anyway
!!     preliminary.  Finally, the cover calculation was performed
!!     twice each timestep, when one calculation suffices.
!!
!!     This scheme calculates cover diagnostically and is called
!!     once at the beginning of each timestep.  It uses the
!!     standard relative humidity calculation from Lohmann and
!!     Roeckner (96).
!!
!! @references.
!!     Diagnostic CC scheme: Lohmann and Roeckner 96, Clim. Dyn.
!!
!! @author A. Tompkins    MPI-Hamburg        2000
!!         K. Ketelesen   NEC,         April 2002
!!         L. Kornblueh   MPI-Hamburg, April 2002
!!
!! @par Revision History
!!    v2: first working version
!!    v5: lookup table added
!!    v8: zriddr and functions replaced for vectorization
!!    v9: optimizations, longer vector loop and less indirect addressing
!!       - introduction of additional arrays
!!       - change structure if "beta function scheme" IF block
!!       - scattered loops "ictit" and "ictdg" are collected over "kproma" and "klev"
!!       - intoducing additional arrays to hold data in "ictit" and "ictdg"
!!         addressing scheme
!! - Taken from ECHAM6.2, wrapped in module and modified for ICON
!!   by Monika Esch, MPI-M (2013-11)
!! - updated to ECHAM.6.3 and unified for use in ICON; removed Tompkins scheme
!!   by Monika Esch, MPI-M (2015-05)
!!
!
MODULE mo_cover

  USE mo_kind,                 ONLY : wp
  USE mo_physical_constants,   ONLY : vtmpc1, cpd, grav
  USE mo_echam_convect_tables, ONLY : prepare_ua_index_spline,lookup_ua_eor_uaw_spline
#ifndef __ICON__
  USE mo_echam_cloud_params,   ONLY : jbmin, jbmax, csatsc, crt, crs, nex, nadd, cinv
#else
  USE mo_echam_cloud_config,   ONLY: echam_cloud_config
#endif
#ifdef _PROFILE
  USE mo_profile,              ONLY : trace_start, trace_stop
#endif

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: cover

#ifdef __ICON__
  ! to simplify access to components of echam_cloud_config
  INTEGER , POINTER :: jbmin, jbmax, nex, nadd
  REAL(wp), POINTER :: csatsc, crt, crs, cinv
#endif


CONTAINS
  !>
  !!
  SUBROUTINE cover (         kproma,   kbdim, ktdia, klev, klevp1                  & !in
    &                      , ktype,    pfrw,     pfri                              & !in
    &                      , paphm1,   papm1,    pgeo                              & !in
    &                      , ptm1,     pqm1,     pxim1                             & !in
    &                      , paclc                                                 & !inout
    &                      , printop                                               & !out
    &              )
    !---------------------------------------------------------------------------------
    !
    INTEGER, INTENT(IN)    :: kbdim, klevp1, klev, kproma, ktdia
    INTEGER, INTENT(IN)    ::  &
      & ktype(kbdim)          !< type of convection
    REAL(wp),INTENT(IN)    ::  &
      & pfrw(kbdim)         ,&!< water mask
      & pfri(kbdim)           !< ice mask
    REAL(wp),INTENT(IN)    ::  &
      & paphm1(kbdim,klevp1),&!< pressure at half levels                   (n-1)
      & papm1(kbdim,klev)   ,&!< pressure at full levels                   (n-1)
      & pgeo(kbdim,klev)    ,&!<
      & pqm1(kbdim,klev)    ,&!< specific humidity                         (n-1)
      & ptm1(kbdim,klev)    ,&!< temperature                               (n-1)
      & pxim1(kbdim,klev)     !< cloud ice                                 (n-1)
    REAL(wp),INTENT(INOUT) ::  &
      & paclc(kbdim,klev)     !< cloud cover
    REAL(wp),INTENT(OUT)   ::  &
      & printop(kbdim)

    INTEGER :: jl, jk, jb, jbn
    INTEGER :: locnt, nl, ilev
    REAL(wp):: zdtdz, zcor, zrhc, zsat, zqr
    INTEGER :: itv1(kproma*klev), itv2(kproma*klev)

    !
    !   Temporary arrays
    !
    REAL(wp)   ::  zdtmin(kbdim), za(kbdim)
    !
    !   Pointers and counters for iteration and diagnostic loop:
    !
    INTEGER :: nphase
    !
    !   variables required for zriddr iteration scheme:
    !
    REAL(wp) :: zjk, zgam
    REAL(wp) :: zqsm1(kproma*klev)
    REAL(wp) :: ua(kproma)

    LOGICAL :: lao, lao1

    REAL(wp) :: zpapm1i(kbdim,klev),     ztmp(kproma*klev)

    REAL(wp) :: zknvb(kbdim),            zphase(kbdim)

    INTEGER :: knvb(kbdim), loidx(kproma*klev)

#ifdef __ICON__
    ! to simplify access to components of echam_cloud_config
    jbmin  => echam_cloud_config% jbmin
    jbmax  => echam_cloud_config% jbmax
    csatsc => echam_cloud_config% csatsc
    crs    => echam_cloud_config% crs
    crt    => echam_cloud_config% crt
    nex    => echam_cloud_config% nex
    nadd   => echam_cloud_config% nadd
    cinv   => echam_cloud_config% cinv
#endif

#ifdef _PROFILE
    CALL trace_start ('cover', 9)
#endif
    !
    !   Initialize variables
    !
    DO jl = 1,kproma
      zdtmin(jl) = -cinv * grav/cpd   ! fraction of dry adiabatic lapse rate
      zknvb(jl)  = 1.0_wp
      printop(jl)= 0.0_wp
    END DO
    !
    DO jk = ktdia,klev
      DO jl = 1,kproma
         zpapm1i(jl,jk) = SWDIV_NOCHK(1._wp,papm1(jl,jk))
      END DO
    END DO
    !
    !       1.3   Checking occurrence of low-level inversion
    !             (below 2000 m, sea points only, no convection)
    !
    locnt = 0
    DO jl = 1,kproma
      IF (pfrw(jl).GT.0.5_wp.AND.pfri(jl).LT.1.e-12_wp.AND.ktype(jl).EQ.0) THEN
        locnt = locnt + 1
        loidx(locnt) = jl
      END IF
    END DO

    IF (locnt.GT.0) THEN
      DO jk = klev,jbmin,-1

!IBM* ASSERT(NODEPS)
!IBM* novector
        DO nl = 1,locnt
          jl = loidx(nl)
          ztmp(nl) = (ptm1(jl,jk-1)-ptm1(jl,jk))*grav/(pgeo(jl,jk-1)-pgeo(jl,jk))
        END DO

        zjk = REAL(jk,wp)
!IBM* ASSERT(NODEPS)
        DO nl = 1,locnt
          jl = loidx(nl)
          zdtdz       = MIN(0.0_wp, ztmp(nl))
          zknvb(jl)   = FSEL(zdtmin(jl)-zdtdz,zknvb(jl),zjk)
          zdtmin(jl)  = MAX(zdtdz,zdtmin(jl))
        END DO
      END DO
    END IF
    knvb(1:kproma) = INT(zknvb(1:kproma))
    !
    !       1.   Calculate the saturation mixing ratio
    !
    IF (ktdia < klev+1) THEN

      DO jk = ktdia,klev

        CALL prepare_ua_index_spline('cover (2)',kproma,ptm1(1,jk),itv1(1),          &
                                         za(1),pxim1(1,jk),nphase,zphase,itv2)
        CALL lookup_ua_eor_uaw_spline(kproma,itv1(1),za(1),nphase,itv2(1),ua(1))

!IBM* novector
        DO jl = 1,kproma
          zqsm1(jl) = MIN(ua(jl)*zpapm1i(jl,jk),0.5_wp)
          zcor      = 1._wp/(1._wp-vtmpc1*zqsm1(jl))
          zqsm1(jl) = zqsm1(jl)*zcor
        END DO
        !
        !       Threshold relative humidity, qsat and cloud cover
        !       This is from cloud, and is the original calculation for
        !       cloud cover, based on relative humidity
        !       (Lohmann and Roeckner, Clim. Dyn.  96)
        !
        DO jl = 1,kproma
        !
          zrhc=crt+(crs-crt)*EXP(1._wp-(paphm1(jl,klevp1)/papm1(jl,jk))**nex)
          zsat=1._wp
          jb=knvb(jl)
          jbn=jb+nadd                  ! mo_echam_cloud_params: nadd=0 except for T31
          lao=(jb.GE.jbmin .AND. jb.LE.jbmax)
          lao1=(jk.EQ.jb .OR. jk.EQ.jbn)
          ilev=klev
          IF (lao .AND. lao1) THEN
          !  ilev=klevp1-jb
            ilev=100
            printop(jl)=REAL(ilev,wp)
            zdtdz = (ptm1(jl,jb-1)-ptm1(jl,jb))*grav/(pgeo(jl,jb-1)-pgeo(jl,jb))
            zgam  = MAX(0.0_wp,-zdtdz*cpd/grav)
            zsat  = MIN(1.0_wp,csatsc+zgam)
          END IF
          zqr=pqm1(jl,jk)/(zqsm1(jl)*zsat)
          paclc(jl,jk)=(zqr-zrhc)/(1.0_wp-zrhc)
          paclc(jl,jk)=MAX(MIN(paclc(jl,jk),1.0_wp),0.0_wp)
          paclc(jl,jk)=1._wp-SQRT(1._wp-paclc(jl,jk))
        END DO !jl
      END DO  !jk
    END IF

#ifdef _PROFILE
    CALL trace_stop ('cover', 9)
#endif
    !
    !
  END SUBROUTINE cover

END MODULE mo_cover
