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
  USE mo_echam_cld_config,     ONLY : echam_cld_config

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: cover

CONTAINS
  !>
  !!
  SUBROUTINE cover (         jg                                                    & !in
    &                      , jb                                                    & !in
    &                      , jcs,      kproma,   kbdim, klev, klevp1               & !in
    &                      , ktype,    pfrw,     pfri                              & !in
    &                      , zf                                                    & !in
    &                      , paphm1,   papm1                                       & !in
    &                      , ptm1,     pqm1,     pxim1                             & !in
    &                      , paclc                                                 & !inout
    &                      , printop                                               & !out
    &              )
    !---------------------------------------------------------------------------------
    !
    INTEGER, INTENT(in)    :: jg
    INTEGER, INTENT(in)    :: jb                    !< number of block
    INTEGER, INTENT(in)    :: kbdim, klevp1, klev, jcs, kproma
    INTEGER, INTENT(in)    :: ktype(kbdim)          !< type of convection
    REAL(wp),INTENT(in)    :: pfrw(kbdim)         ,&!< water mask
      &                       pfri(kbdim)           !< ice mask
    REAL(wp),INTENT(in)    :: zf(kbdim,klev)      ,&!< geometric height thickness [m]
      &                       paphm1(kbdim,klevp1),&!< pressure at half levels
      &                       papm1(kbdim,klev)   ,&!< pressure at full levels
      &                       pqm1(kbdim,klev)    ,&!< specific humidity
      &                       ptm1(kbdim,klev)    ,&!< temperature
      &                       pxim1(kbdim,klev)     !< cloud ice
    REAL(wp),INTENT(out)   :: paclc(kbdim,klev)     !< cloud cover
    REAL(wp),INTENT(out)   :: printop(kbdim)

    INTEGER :: jl, jk, jbm
    INTEGER :: locnt, nl, ilev
    REAL(wp):: zdtdz, zcor, zrhc, zsat, zqr
    INTEGER :: itv1(kproma), itv2(kproma)

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

    ! Shortcuts to components of echam_cld_config
    !
    INTEGER , POINTER :: jks, jbmin, jbmax, nex
    REAL(wp), POINTER :: csatsc, crt, crs, cinv
    !
    jks    => echam_cld_config(jg)% jks
    jbmin  => echam_cld_config(jg)% jbmin
    jbmax  => echam_cld_config(jg)% jbmax
    csatsc => echam_cld_config(jg)% csatsc
    crs    => echam_cld_config(jg)% crs
    crt    => echam_cld_config(jg)% crt
    nex    => echam_cld_config(jg)% nex
    cinv   => echam_cld_config(jg)% cinv

    !
    !   Initialize variables
    !
    DO jk = 1,jks-1
      DO jl = jcs,kproma
         paclc(jl,jk) = 0.0_wp
      END DO
    END DO
    !
    DO jl = jcs,kproma
      zdtmin(jl) = -cinv * grav/cpd   ! fraction of dry adiabatic lapse rate
      zknvb(jl)  = 1.0_wp
      printop(jl)= 0.0_wp
    END DO
    !
    DO jk = jks,klev
      DO jl = jcs,kproma
         zpapm1i(jl,jk) = SWDIV_NOCHK(1._wp,papm1(jl,jk))
      END DO
    END DO
    !
    !       1.3   Checking occurrence of low-level inversion
    !             (below 2000 m, sea points only, no convection)
    !
    locnt = jcs-1
    DO jl = jcs,kproma
      IF (pfrw(jl).GT.0.5_wp.AND.pfri(jl).LT.1.e-12_wp.AND.ktype(jl).EQ.0) THEN
        locnt = locnt + 1
        loidx(locnt) = jl
      END IF
    END DO

    IF (locnt.GT.jcs-1) THEN
      DO jk = klev,jbmin,-1

!IBM* ASSERT(NODEPS)
!IBM* novector
        DO nl = jcs,locnt
          jl = loidx(nl)
          ztmp(nl) = (ptm1(jl,jk-1)-ptm1(jl,jk))/(zf(jl,jk-1)-zf(jl,jk))
        END DO

        zjk = REAL(jk,wp)
!IBM* ASSERT(NODEPS)
        DO nl = jcs,locnt
          jl = loidx(nl)
          zdtdz       = MIN(0.0_wp, ztmp(nl))
          zknvb(jl)   = FSEL(zdtmin(jl)-zdtdz,zknvb(jl),zjk)
          zdtmin(jl)  = MAX(zdtdz,zdtmin(jl))
        END DO
      END DO
    END IF
    knvb(jcs:kproma) = INT(zknvb(jcs:kproma))
    !
    !       1.   Calculate the saturation mixing ratio
    !
    IF (jks < klev+1) THEN

      DO jk = jks,klev

        CALL prepare_ua_index_spline(jg,'cover (2)',jcs,kproma,ptm1(1,jk),itv1(1),   &
                                         za(1),pxim1(1,jk),nphase,zphase,itv2,       &
                                         klev=jk,kblock=jb,kblock_size=kbdim)
        ! output: itv1=idx, za=zalpha, nphase, zphase, itv2
        CALL lookup_ua_eor_uaw_spline(jcs,kproma,itv1(1),za(1),nphase,itv2(1),ua(1))
        ! output: ua

!IBM* novector

        DO jl = jcs,kproma
          zqsm1(jl) = MIN(ua(jl)*zpapm1i(jl,jk),0.5_wp)
          zcor      = 1._wp/(1._wp-vtmpc1*zqsm1(jl))
          zqsm1(jl) = zqsm1(jl)*zcor  ! qsat
        END DO
        !
        !       Threshold relative humidity, qsat and cloud cover
        !       This is from cloud, and is the original calculation for
        !       cloud cover, based on relative humidity
        !       (Lohmann and Roeckner, Clim. Dyn.  96)
        !
        DO jl = jcs,kproma
        !
          zrhc=crt+(crs-crt)*EXP(1._wp-(paphm1(jl,klevp1)/papm1(jl,jk))**nex)
          zsat=1._wp
          jbm=knvb(jl)
          lao=(jbm.GE.jbmin .AND. jbm.LE.jbmax)
          lao1=(jk.EQ.jbm)
          ilev=klev
          IF (lao .AND. lao1) THEN
          !  ilev=klevp1-jbm
            ilev=100
            printop(jl)=REAL(ilev,wp)
            zdtdz = (ptm1(jl,jbm-1)-ptm1(jl,jbm))/(zf(jl,jk-1)-zf(jl,jk))
            zgam  = MAX(0.0_wp,-zdtdz*cpd/grav)
            zsat  = MIN(1.0_wp,csatsc+zgam)
          END IF
          zqr=pqm1(jl,jk)/(zqsm1(jl)*zsat)      ! r in Lohmann-scheme (grid-mean rel hum)
          paclc(jl,jk)=(zqr-zrhc)/(1.0_wp-zrhc) ! b_o in Lohman-scheme
          paclc(jl,jk)=MAX(MIN(paclc(jl,jk),1.0_wp),0.0_wp)
          paclc(jl,jk)=1._wp-SQRT(1._wp-paclc(jl,jk))
        END DO !jl
      END DO  !jk
    END IF

    !
    !
  END SUBROUTINE cover

END MODULE mo_cover
