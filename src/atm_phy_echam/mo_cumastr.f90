#ifdef __xlC__
@PROCESS STRICT
#endif
#include "fsel.inc"

!>
!! @brief Master routine - provides interface for cumulus parameterization
!!
!! @remarks
!!     *cumastr*
!!     This routine computes the physical tendencies of the prognostic variables
!!     t,q,u and v due to convective processes. Processes considered are: convective
!!     fluxes, formation of precipitation, evaporation of falling rain below cloud
!!     base, saturated cumulus downdrafts.
!!     It takes its input from the long-term storage t,q,u,v,phi and p and moisture
!!     tendencies. It returns its output to the same space
!!     1. modified tendencies of model variables
!!     2. rates of convective precipitation (for surface)
!!     Method: Parameterization is done using a massflux-scheme.
!!        (1) Define constants and parameters
!!        (2) Specify values (t,q,qs...) at half levels and
!!            initialize updraft- and downdraft-values in 'cuini'
!!        (3) Calculate cloud base in 'cubase'
!!            and specify cloud base massflux from pbl moisture budget
!!        (4) Do cloud ascent in 'cuasc' in absence of downdrafts
!!        (5) Do downdraft calculations:
!!              (A) Determine values at lfs in 'cudlfs'
!!              (B) Determine moist descent in 'cuddraf'
!!              (C) Recalculate cloud base massflux considering the
!!                  effect of cumulus-downdrafts
!!        (6) Do final cloud ascent in 'cuasc'
!!        (7) Do final adjustments to convective fluxes in 'cuflx',
!!            do evaporation in subcloud layer
!!        (8) Calculate increments of t and q in 'cudtdq'
!!        (9) Calculate increments of u and v in 'cududv'
!!
!! @author M. Tiedtke, ECMWF,    1986/1987/1989
!!         M. Tiedtke, ECMWF,    Dec 1989
!!
!! @references.
!!       Tiedtke, 1989: Mon. Wea. Rev., 117, 1779-1800
!!
!! @par Revision History
!! - Taken from ECHAM6.3, unified for ICON/ECHAM by Monika Esch, MPI-M (2015-06)
!!
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_cumastr

  USE mo_kind,                 ONLY: wp
  USE mo_echam_convect_tables, ONLY: prepare_ua_index_spline,lookup_ua_spline, lookup_ubc
  USE mo_physical_constants,   ONLY: grav, alv, als, tmelt, vtmpc1, rd
  USE mo_echam_cnv_config,     ONLY: echam_cnv_config, cevapcu
  USE mo_cuinitialize,         ONLY: cuini, cubase
  USE mo_cuascent,             ONLY: cuasc
  USE mo_cudescent,            ONLY: cudlfs, cuddraf
  USE mo_cufluxdts,            ONLY: cuflx, cudtdq, cududv

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: cumastr

CONTAINS
  !>
  !!
  SUBROUTINE cumastr(  jg, jb,                                            &
    &                  jcs, kproma,   kbdim,                              &
    &                  klev,     klevp1,   klevm1,                        &
    &                  pdtime,                                            &
    &                  pzf,      pzh,                                     &
    &                  pmref,                                             &
    &                  pten,     pqen,     pxen,     puen,     pven,      &
    &                  ktrac,    ldland,                                  &
    &                  pxten,                                             &
    &                  pverv,    pqhfla,                                  &
    &                  papp1,    paphp1,                                  &
    &                  pgeo,     pgeoh,                                   &
    &                  pqte,                                              &
    &                  pthvsig,                                           &
    &                  ktype,    kctop,                                   &
    &                  prsfc,    pssfc,                                   &
    &                  pcon_dtrl,pcon_dtri,pcon_iqte,                     &
    &                  pq_cnv,   pvom_cnv, pvol_cnv, pqte_cnv, pxtte_cnv, &
    &                  pxtecl,   pxteci,                                  &
    &                  ptop                                               )
    !
    INTEGER, INTENT(IN)   :: jg, jb
    INTEGER, INTENT(IN)   :: jcs, kproma, kbdim, klev, klevp1, ktrac, klevm1
    REAL(wp),INTENT(IN)   :: pdtime
    REAL(wp),INTENT(IN)   :: pzf(kbdim,klev),         pzh(kbdim,klevp1)
    REAL(wp),INTENT(IN)   :: pmref(kbdim,klev)

    REAL(wp),INTENT(IN)   :: pten(kbdim,klev),        pqen(kbdim,klev),        &
      &                      pxen(kbdim,klev),        pxten(kbdim,klev,ktrac), &
      &                      puen(kbdim,klev),        pven(kbdim,klev),        &
      &                      pverv(kbdim,klev),       pqhfla(kbdim),           &
      &                      papp1(kbdim,klev),       paphp1(kbdim,klevp1),    &
      &                      pgeo(kbdim,klev),        pgeoh(kbdim,klevp1),     &
      &                      pqte(kbdim,klev),                                 &
      &                      pthvsig(kbdim)

    INTEGER, INTENT(OUT)  :: ktype(kbdim),         kctop(kbdim)
    REAL(wp),INTENT(OUT)  :: pcon_dtrl(kbdim),     pcon_dtri(kbdim)
    REAL(wp),INTENT(OUT)  :: pcon_iqte(kbdim)
    REAL(wp),INTENT(OUT)  :: pq_cnv(kbdim,klev)
    REAL(wp),INTENT(OUT)  :: pvom_cnv(kbdim,klev), pvol_cnv(kbdim,klev)
    REAL(wp),INTENT(OUT)  :: pqte_cnv(kbdim,klev), pxtte_cnv(kbdim,klev,ktrac)
    REAL(wp),INTENT(OUT)  :: prsfc(kbdim),         pssfc(kbdim)
    REAL(wp),INTENT(OUT)  :: pxtecl(kbdim,klev),   pxteci(kbdim,klev)
    REAL(wp),INTENT(OUT)  :: ptop(kbdim)

    !
    REAL(wp):: ptu(kbdim,klev),         pqu(kbdim,klev),                   &
      &        plu(kbdim,klev),         plude(kbdim,klev),                 &
      &        pmfu(kbdim,klev),        pmfd(kbdim,klev)
    INTEGER :: kcbot(kbdim)
    REAL(wp):: pqude(kbdim,klev)
    REAL(wp):: zqsen(kbdim,klev),       zcpen(kbdim,klev)
    REAL(wp):: ztenh(kbdim,klev),       zqenh(kbdim,klev),                 &
      &        zxenh(kbdim,klev),       zalvsh(kbdim,klev),                &
      ! zalvsh: latent heat of vaporisation/sublimation defined at half levels
      &        zqsenh(kbdim,klev),                                         &
      &        ztd(kbdim,klev),         zqd(kbdim,klev),                   &
      &        zmfus(kbdim,klev),       zmfds(kbdim,klev),                 &
      &        zmfuq(kbdim,klev),       zmfdq(kbdim,klev),                 &
      &        zdmfup(kbdim,klev),      zdmfdp(kbdim,klev),                &
      &        zmful(kbdim,klev),       zrfl(kbdim),                       &
      &        zuu(kbdim,klev),         zvu(kbdim,klev),                   &
      &        zud(kbdim,klev),         zvd(kbdim,klev)
    REAL(wp):: zcpcu(kbdim,klev)
    REAL(wp):: zentr(kbdim),            zhcbase(kbdim),                    &
      &        zmfub(kbdim),            zmfub1(kbdim),                     &
      &        zktype(kbdim),           zldcum(kbdim),                     &
      &        zcpcui(kbdim,klev),      zkcbot(kbdim),                     &
      &        zictop0(kbdim),          ztmp1(kbdim),                      &
      &        ztmp2(kbdim),            ztmp3(kbdim),                      &
      &        za(kbdim)
    REAL(wp):: zsfl(kbdim),             zdpmel(kbdim,klev)
    REAL(wp):: zcape(kbdim),            zheat(kbdim)
    REAL(wp):: ua(kbdim), dua(kbdim), ub(kbdim)
    REAL(wp):: zhmin(kbdim),            zihmin(kbdim)
    REAL(wp):: zhhatt(kbdim,klev),      zdqpbl(kbdim)
    REAL(wp):: zdqcv(kbdim)
    INTEGER :: ihmin(kbdim), ilo1(kbdim), ldidx(kbdim), loidx(kbdim)
    INTEGER :: ilab(kbdim,klev),        idtop(kbdim),                      &
      &        ictop0(kbdim),           ilwmin(kbdim)
    INTEGER :: itopec2(kbdim)
    REAL(wp):: zxtu(kbdim,klev,ktrac),                                     &
      &        zxtenh(kbdim,klev,ktrac),zxtd(kbdim,klev,ktrac),            &
      &        zmfuxt(kbdim,klev,ktrac),zmfdxt(kbdim,klev,ktrac)
    LOGICAL :: loddraf(kbdim),          ldland(kbdim)
    LOGICAL :: ldcum(kbdim)
    LOGICAL :: llo1, llo1mask(kbdim)
    !
    INTEGER :: nl, jl, jk, ikb, jt, itopm2, locnt, ldcnt
    REAL(wp):: zcons, zqumqe, zdqmin, zmfmax, zalvs, zalvdcp, zqalv       &
      &      , zhsat, zes, zcor, zqsat, zdqsdt, zgam, zzz, zhhat           &
      &      , zbi, zroi, zdz, zdhdz, zdepth, zfac, zrh, zeps              &
      &      , zjk, zhelp, zlo1, zpaphp1i, ztenhi, zkctop, zpbmpt, za1, za2
    !
    !  INTRINSIC FUNCTIONS
    INTRINSIC MIN, MAX
    !
    !  Executable statements

    ! Shortcuts to components of echam_cnv_config
    !
    LOGICAL           :: lmfdd, lmfdudv
    REAL(wp)          :: entrpen, entrscv, cmfdeps, cmftau
    REAL(wp), POINTER :: zevapcu(:)
    !
    lmfdd   = echam_cnv_config(jg)% lmfdd
    lmfdudv = echam_cnv_config(jg)% lmfdudv
    entrscv = echam_cnv_config(jg)% entrscv
    entrpen = echam_cnv_config(jg)% entrpen
    cmfdeps = echam_cnv_config(jg)% cmfdeps
    cmftau  = echam_cnv_config(jg)% cmftau
    !
    zevapcu(1:klev) => cevapcu(1:klev,jg)

    !$ACC DATA PRESENT( pxten ) IF( ktrac > 0 )
    !$ACC DATA PRESENT( pxtte_cnv ) IF( ktrac> 0 )
    !$ACC DATA PRESENT( pzf, pzh, pmref, pten, pqen, pxen, puen, pven,        &
    !$ACC               ldland, pverv, pqhfla, papp1, paphp1, pgeo,           &
    !$ACC               pgeoh, pqte, pthvsig, ktype, kctop, prsfc, pssfc,     &
    !$ACC               pcon_dtrl, pcon_dtri, pcon_iqte, pq_cnv, pvom_cnv,    &
    !$ACC               pvol_cnv, pqte_cnv, pxtecl, pxteci, ptop, zevapcu )   &
    !$ACC       CREATE( ptu, pqu, plu, plude, pmfu, pmfd, kcbot, pqude, zqsen,&
    !$ACC               zcpen, ztenh, zqenh, zxenh, zalvsh, zqsenh, ztd, zqd, &
    !$ACC               zmfus, zmfds, zmfuq, zmfdq, zdmfup, zdmfdp, zmful,    &
    !$ACC               zrfl, zuu, zvu, zud, zvd, zcpcu, zentr, zhcbase,      &
    !$ACC               zmfub, zmfub1, zktype, zldcum, zcpcui, zkcbot,        &
    !$ACC               zictop0, ztmp1, ztmp2, ztmp3, za, zsfl, zdpmel, zcape,&
    !$ACC               zheat, ua, dua, ub, zhmin, zihmin, zhhatt, zdqpbl,    &
    !$ACC               zdqcv, ihmin, ilo1, ldidx, loidx, ilab, idtop, ictop0,&
    !$ACC               ilwmin, itopec2, loddraf, ldland, ldcum, llo1mask )
    !$ACC DATA  CREATE( zxtu, zxtenh, zxtd, zmfuxt, zmfdxt ) IF( ktrac > 0 )

    !-----------------------------------------------------------------------
    !
    !     1.           Specify constants and parameters
    !                  --------------------------------
    !
    zcons=1._wp/pdtime
    !
    !----------------------------------------------------------------------
    !
    !*    2.           Initialize values at vertical grid points in 'cuini'
    !                  ---------------------------------------------------
    !
    CALL cuini(jg, jb,                                                   &
      &        jcs, kproma, kbdim, klev, klevp1, klevm1,                 &
      &        pten,     pqen,     zqsen,    pxen,     puen,     pven,   &
      &        ktrac,                                                    &
      &        pxten,    zxtenh,   zxtu,     zxtd,     zmfuxt,   zmfdxt, &
      &        pverv,    papp1,    pgeo,     paphp1,   pgeoh,            &
      &        ztenh,    zqenh,    zqsenh,   zxenh,    ilwmin,           &
      &        ptu,      pqu,      ztd,      zqd,                        &
      &        zuu,      zvu,      zud,      zvd,                        &
      &        pmfu,     pmfd,     zmfus,    zmfds,                      &
      &        zmfuq,    zmfdq,    zdmfup,   zdmfdp,                     &
      &        zcpen,    zcpcu,    zalvsh,                               &
      &        zdpmel,   plu,      plude,    pqude,    ilab             )
    !
    !-----------------------------------------------------------------------
    !
    !*    3.0          Cloud base calculations
    !                  -----------------------
    !
    !*             (A) Determine cloud base values in 'cubase'
    !                  ---------------------------------------
    !
    CALL cubase(jg, jb,                                                  &
      &         jcs, kproma,   kbdim,    klev,     klevp1,    klevm1,    &
      &         ztenh,    zqenh,    pgeoh,    paphp1,    pthvsig,        &
      &         ptu,      pqu,      plu,                                 &
      &         puen,     pven,     zuu,      zvu,                       &
      &         zcpcu,                                                   &
      &         ldcum,    kcbot,    ilab)
    !
    !*             (B) Determine total moisture convergence and
    !*                 then decide on type of cumulus convection
    !                  -----------------------------------------
    !
    !$ACC PARALLEL DEFAULT(PRESENT)
    !$ACC LOOP GANG VECTOR
    DO jl = jcs, kproma
      zkcbot(jl) = REAL(kcbot(jl),wp)
    END DO
    !$ACC END PARALLEL

    jk=1
    !$ACC PARALLEL DEFAULT(PRESENT)
    !$ACC LOOP GANG VECTOR
    DO jl=jcs,kproma
      zdqpbl(jl)=0.0_wp
      zdqcv(jl)=pqte(jl,jk)*pmref(jl,jk)
      idtop(jl)=0
    END DO
    !$ACC END PARALLEL
    !$ACC PARALLEL DEFAULT(PRESENT)
    !$ACC LOOP SEQ
    DO jk=2,klev
      zjk = REAL(jk,wp)
      !$ACC LOOP GANG VECTOR PRIVATE( zhelp )
      DO jl=jcs,kproma
        zhelp      = pmref(jl,jk)
        zdqcv(jl)  = zdqcv(jl)+pqte(jl,jk)*zhelp
        zdqpbl(jl) = zdqpbl(jl) + FSEL(zjk - zkcbot(jl),pqte(jl,jk),0._wp)*zhelp
#ifdef __ibmdbg__
        PRINT '(A6,I3,I4,E18.10,L3)','zdqpbl',jk,jl,zdqpbl(jl),(FSEL(zjk - zkcbot(jl),1._wp,0._wp)==1._wp)
#endif
      END DO
    END DO
    !$ACC END PARALLEL
    !
    !*             (C) Determine moisture supply for boundary layer and determine
    !*                 cloud base massflux ignoring the effects of downdrafts
    !*                 at this stage
    !                  ---------------------------------------------------------------
    !
    !$ACC PARALLEL DEFAULT(PRESENT)
    !$ACC LOOP GANG VECTOR
    DO jl = jcs, kproma
      zldcum(jl) = MERGE(1._wp,0._wp,ldcum(jl))
    END DO
    !$ACC END PARALLEL
    !$ACC PARALLEL DEFAULT(PRESENT)
    !$ACC LOOP GANG VECTOR
    DO jl = 1, kbdim
      ktype(jl) = 0
    END DO
    !$ACC END PARALLEL

!DIR$ IVDEP
#ifdef _CRAYFTN
!DIR$ CONCURRENT
#endif
!IBM* ASSERT(NODEPS)
!IBM* NOVECTOR
    !$ACC PARALLEL DEFAULT(PRESENT)
    !$ACC LOOP GANG VECTOR PRIVATE( ikb, zqumqe, zdqmin, zlo1, zmfmax, zhelp )
    DO jl=jcs,kproma
      ikb=kcbot(jl)
      zqumqe=pqu(jl,ikb)+plu(jl,ikb)-zqenh(jl,ikb)
      zdqmin=MAX(0.01_wp*zqenh(jl,ikb),1.e-10_wp)
      zlo1 = FSEL(-zdqpbl(jl),0._wp,1._wp)
      zlo1 = FSEL(zdqmin - zqumqe,0._wp,zlo1) * zldcum(jl)
      zmfub(jl)=FSEL(-zlo1,0.01_wp,zdqpbl(jl)/MAX(zqumqe,zdqmin))
      zmfmax=pmref(jl,ikb-1)*zcons
      zmfub(jl)=MIN(zmfub(jl),zmfmax)
      zldcum(jl) = zlo1
      zhelp = MAX(0._wp,-1.1_wp*pqhfla(jl))
      zktype(jl) = FSEL(zhelp - zdqcv(jl),2._wp,1._wp)
      zentr(jl)  = FSEL(zhelp - zdqcv(jl),entrscv,entrpen)
      ktype(jl) = INT(zktype(jl))
#ifdef __ibmdbg__
      PRINT '(A6,I3,I4,2 E18.10,I3,L3)','zmfub',ikb,jl,zmfub(jl),zentr(jl),ktype(jl),(zlo1==1._wp)
#endif
    END DO
    !$ACC END PARALLEL
    !$ACC PARALLEL DEFAULT(PRESENT)
    !$ACC LOOP GANG VECTOR
    DO jl = jcs, kproma
      ldcum(jl) = (zldcum(jl).GT.0._wp)
    END DO
    !$ACC END PARALLEL
    !
    !-----------------------------------------------------------------------
    !*    4.0          Determine cloud ascent for entraining plume
    !                  -------------------------------------------
    !
    !*   (A) Estimate cloud height for entrainment/detrainment calculations in
    !*       cuasc (max. possible cloud height for non-entraining plume,
    !*       following A.-S.,1974)
    !        -------------------------------------------------------------------------
    !
!DIR$ IVDEP
    !$ACC PARALLEL DEFAULT(PRESENT)
    !$ACC LOOP GANG VECTOR PRIVATE( ikb, zalvs )
    DO jl=jcs,kproma
      ikb=kcbot(jl)
      zalvs=FSEL(tmelt-ptu(jl,ikb),als,alv)
      zhcbase(jl) = zcpcu(jl,ikb)*ptu(jl,ikb) + pgeoh(jl,ikb) + zalvs*pqu(jl,ikb)
      zictop0(jl) = zkcbot(jl)-1._wp
    END DO
    !$ACC END PARALLEL
    DO jk=klev,1,-1
      !$ACC PARALLEL DEFAULT(PRESENT)
      !$ACC LOOP GANG VECTOR
      DO jl = jcs, kproma
        zcpcui(jl,jk) = 1._wp/zcpcu(jl,jk)
      END DO
      !$ACC END PARALLEL
      IF (jk <= klevm1 .AND. jk >= 3) THEN
        ! mpuetz: too few instructions (FP dependencies)
        CALL prepare_ua_index_spline(jg,'cumastr',jcs,kproma,ztenh(:,jk),loidx(:),za(:))
        CALL lookup_ua_spline(jcs,kproma,loidx(:),za(:),ua(:),dua(:))
        CALL lookup_ubc(jcs,kproma,ztenh(:,jk),ub(:))
        zjk = REAL(jk,wp)
!IBM* NOVECTOR
        !$ACC PARALLEL DEFAULT(PRESENT)
        !$ACC LOOP GANG VECTOR PRIVATE( zalvs, zalvdcp, zqalv, zhsat, zpaphp1i, zes, zcor, &
        !$ACC                           zqsat, zdqsdt, zgam, zzz, zhhat, zlo1 )
        DO jl=jcs,kproma
          ! mpuetz: move some of these into the previous loop
          zalvs=FSEL(tmelt - ztenh(jl,jk),als,alv)
          zalvdcp=zalvs*zcpcui(jl,jk)
          zqalv=1._wp/zalvs
          zhsat=zcpcu(jl,jk)*ztenh(jl,jk)+pgeoh(jl,jk)+zalvs*zqsenh(jl,jk)
          zpaphp1i = 1._wp/paphp1(jl,jk)
          zes=ua(jl)*zpaphp1i
          zes=MIN(0.5_wp,zes)
          zcor=1._wp/(1._wp-vtmpc1*zes)
          zqsat=zes*zcor
          zdqsdt=zpaphp1i*zcor**2*dua(jl)
          zgam=FSEL(zes - 0.4_wp,zqsat*zcor*ub(jl),zalvdcp*zdqsdt)
          zzz=zcpcu(jl,jk)*ztenh(jl,jk)*vtmpc1
          zhhat=zhsat-((zzz+zgam*zzz)/(1._wp+zgam*zzz*zqalv))*             &
            &  MAX(zqsenh(jl,jk)-zqenh(jl,jk),0._wp)
          zhhatt(jl,jk)=zhhat
          zlo1 = FSEL(zjk - zictop0(jl),0._wp,1._wp)
          zlo1 = FSEL(zhhat - zhcbase(jl),0._wp,zlo1)
          zictop0(jl) = FSEL(-zlo1,zictop0(jl),zjk)
#ifdef __ibmdbg__
          PRINT '(A6,I3,I4,3 E18.10,I3)','ictop',jk,jl,zes,zqsdt,zhhatt(jl,jk),INT(zictop0(jl))
#endif
        END DO
        !$ACC END PARALLEL
      END IF
    END DO
    !$ACC PARALLEL DEFAULT(PRESENT)
    !$ACC LOOP GANG VECTOR
    DO jl = jcs, kproma
      ictop0(jl) = INT(zictop0(jl))
    END DO
    !$ACC END PARALLEL
    !!
    !!     DEEP CONVECTION IF CLOUD DEPTH > 200 HPA, ELSE SHALLOW
    !!     (CLOUD DEPTH FROM NON-ENTRAINIG PLUME)
    !!
    !  DO jl=jcs,kproma
    !     ktype(jl)=MERGE(1,2,                                              &
    !                paphp1(jl,kcbot(jl))-paphp1(jl,ictop0(jl)).gt.2.E4_wp)
    !     zentr(jl)=MERGE(entrpen,entrscv,ktype(jl).eq.1)
    !  ENDDO
    !!
    !              Find lowest possible org. detrainment level
    !              -------------------------------------------
    !
    ldcnt = jcs-1
    !$ACC PARALLEL DEFAULT(PRESENT)
    !$ACC LOOP GANG VECTOR
    DO jl=jcs,kproma
      zhmin(jl)=0._wp
      zihmin(jl)=0._wp
      llo1mask(jl)=ldcum(jl).AND.ktype(jl).EQ.1
      IF(llo1mask(jl)) THEN
        zihmin(jl)=zkcbot(jl)
      ENDIF
    ENDDO
    !$ACC END PARALLEL
    !$ACC UPDATE HOST( llo1mask )
    DO jl=jcs,kproma
      IF(llo1mask(jl)) THEN
        ldcnt = ldcnt + 1
        ldidx(ldcnt) = jl
      ENDIF
    ENDDO
    !$ACC UPDATE DEVICE( ldidx )
    !
    zbi=1._wp/25._wp
    DO jk=klev,1,-1
      ! mpuetz: compute the update criterion
      zjk = REAL(jk,wp)
!IBM* ASSERT(NODEPS)
      !$ACC PARALLEL DEFAULT(PRESENT)
      !$ACC LOOP GANG VECTOR PRIVATE( jl, zlo1 )
      DO nl = jcs,ldcnt
        jl = ldidx(nl)
        zlo1 = FSEL(zjk - zkcbot(jl),0._wp,1._wp)
        zlo1 = FSEL(zjk - zictop0(jl),zlo1,0._wp)
        zlo1 = FSEL(-ABS(zihmin(jl)-zkcbot(jl)),zlo1,0._wp)
        ilo1(nl) = INT(zlo1)
      END DO
      !$ACC END PARALLEL
      ! mpuetz: compute the indices of elements to be updated
      locnt = 0
      !$ACC UPDATE HOST( ilo1, ldidx )
      DO nl = jcs,ldcnt
        IF (ilo1(nl).GT.0) THEN
          locnt = locnt + 1
          loidx(locnt) = ldidx(nl)
        END IF
      END DO
      !$ACC UPDATE DEVICE( loidx )
!IBM* ASSERT(NODEPS)
      !$ACC PARALLEL DEFAULT(PRESENT)
      !$ACC LOOP GANG VECTOR PRIVATE( jl, zalvs, ikb, zroi, zdz, za1, za2, zdhdz, zdepth )
      DO nl=1,locnt
        jl = loidx(nl)
        zalvs = FSEL(tmelt - ztenh(jl,jk),als,alv)
      !  zalvs=MERGE(alv,als,ztenh(jl,jk)>tmelt)
        ikb   = kcbot(jl)
        zroi  = SWDIV_NOCHK(rd*ztenh(jl,jk)*(1._wp+vtmpc1*zqenh(jl,jk)),paphp1(jl,jk))
        zdz   = pmref(jl,jk-1)*zroi
        za1 = (zcpen(jl,jk-1)*pten(jl,jk-1) - zcpen(jl,jk)*pten(jl,jk)               &
             + zalvs*(pqen(jl,jk-1) - pqen(jl,jk))+(pgeo(jl,jk-1)-pgeo(jl,jk)))
        za2 = pzf(jl,jk-1)-pzf(jl,jk)
        zdhdz = SWDIV_NOCHK(za1, za2)
        zdepth    = pzh(jl,jk)-pzh(jl,ikb)
        ztmp1(nl) = zalvs
        ztmp2(nl) = zdz*zdhdz
        ztmp3(nl) = 1._wp+zdepth*zbi
      END DO
      !$ACC END PARALLEL
      !$ACC PARALLEL DEFAULT(PRESENT)
      !$ACC LOOP GANG VECTOR
      DO nl = 1, locnt
        ztmp3(nl) = SQRT(ztmp3(nl))
      END DO
      !$ACC END PARALLEL
!IBM* ASSERT(NODEPS)
      !$ACC PARALLEL DEFAULT(PRESENT)
      !$ACC LOOP GANG VECTOR PRIVATE( jl, zalvs, zfac, zrh )
      DO nl=1,locnt
        jl = loidx(nl)
        zalvs     = ztmp1(nl)
        zfac      = ztmp3(nl)
        zhmin(jl) = zhmin(jl) + zfac*ztmp2(nl)
        zrh       =-zalvs*(zqsenh(jl,jk)-zqenh(jl,jk))*zfac
        zihmin(jl) = FSEL(zrh - zhmin(jl),zihmin(jl),zjk)
      !  IF(zhmin(jl).GT.zrh) ihmin(jl)=jk
      ENDDO
      !$ACC END PARALLEL
    ENDDO
!IBM* ASSERT(NODEPS)
    !$ACC PARALLEL DEFAULT(PRESENT)
    !$ACC LOOP GANG VECTOR PRIVATE( jl )
    DO nl=jcs,ldcnt
      jl = ldidx(nl)
#ifdef __ibmdbg__
      PRINT '(A6,I4,I4,I4)','ihmin',jl,INT(zihmin(jl)),ictop0(jl)
#endif
      zihmin(jl) = FSEL(zihmin(jl)-zictop0(jl),zihmin(jl),zictop0(jl))
    !    IF(ihmin(jl).LT.ictop0(jl)) ihmin(jl)=ictop0(jl)
    ENDDO
    !$ACC END PARALLEL
    !$ACC PARALLEL DEFAULT(PRESENT)
    !$ACC LOOP GANG VECTOR
    DO jl = jcs, kproma
      ihmin(jl) = INT(zihmin(jl))
    END DO
    !$ACC END PARALLEL
    !
    !*         (B) Do ascent in 'cuasc' in absence of downdrafts
    !              ---------------------------------------------
    !
    CALL cuasc(jg, jb,                                                   &
      &        jcs,  kproma, kbdim, klev, klevp1, klevm1,                &
      &        pzf,      pzh,      pmref,                                &
      &        ztenh,    zqenh,    puen,     pven,                       &
      &        ktrac,                                                    &
      &        pdtime,                                                   &
      &        zxtenh,   pxten,    zxtu,     zmfuxt,                     &
      &        pten,     pqen,     zqsen,                                &
      &        pgeo,     pgeoh,    paphp1,   pthvsig,                    &
      &        pqte,               pverv,    ilwmin,                     &
      &        ldcum,    ldland,   ktype,    ilab,                       &
      &        ptu,      pqu,      plu,      zuu,      zvu,              &
      &        pmfu,     zmfub,    zentr,                                &
      &        zmfus,    zmfuq,                                          &
      &        zmful,    plude,    pqude,    zdmfup,                     &
      &        ihmin,    zhhatt,   zhcbase,  zqsenh,                     &
      &        zcpen,    zcpcu,                                          &
      &        kcbot,    kctop,    ictop0                                &
      &       )
    !
    !*     (C) Check cloud depth and change entrainment rate accordingly
    !          Calculate precipitation rate (for downdraft calculation)
    !          ---------------------------------------------------------
    !
!DIR$ IVDEP
    !$ACC PARALLEL DEFAULT(PRESENT)
    !$ACC LOOP GANG VECTOR PRIVATE( zpbmpt )
    DO jl=jcs,kproma
      zpbmpt=paphp1(jl,kcbot(jl))-paphp1(jl,kctop(jl))
      IF(ldcum(jl).AND.ktype(jl).EQ.1.AND.zpbmpt.LT.2.e4_wp) ktype(jl)=2  ! cloud thickness < 200hPa --> shallow conv.
      IF(ldcum(jl)) ictop0(jl)=kctop(jl)
      IF(ktype(jl).EQ.2) zentr(jl)=entrscv
      zrfl(jl)=zdmfup(jl,1)
    END DO
    !$ACC END PARALLEL
    !$ACC PARALLEL DEFAULT(PRESENT)
    !$ACC LOOP SEQ
    DO jk=2,klev
      !$ACC LOOP GANG VECTOR
      DO jl=jcs,kproma
        zrfl(jl)=zrfl(jl)+zdmfup(jl,jk)
      END DO
    END DO
    !$ACC END PARALLEL
    ! mpuetz: must recompute ldidx() since ktype could have changed
    ldcnt = jcs-1 
    !$ACC PARALLEL DEFAULT(PRESENT)
    !$ACC LOOP GANG VECTOR
    DO jl=jcs,kproma
      llo1mask(jl)=ldcum(jl).AND.ktype(jl).EQ.1
    ENDDO
    !$ACC END PARALLEL
    !$ACC UPDATE HOST( llo1mask )
    DO jl=jcs,kproma
      IF(llo1mask(jl)) THEN
        ldcnt = ldcnt + 1
        ldidx(ldcnt) = jl
      ENDIF
    ENDDO
    !$ACC UPDATE DEVICE( ldidx )
    !
    !-----------------------------------------------------------------------
    !*    5.0          Cumulus downdraft calculations
    !                  ------------------------------
    !
    IF(lmfdd) THEN
      !
      !*             (A) Determine lfs in 'cudlfs'
      !                  -------------------------
      CALL cudlfs(jg,       jb,                                         &
        &         jcs,     kproma,   kbdim,    klev,     klevp1,        &
        &         ztenh,    zqenh,    puen,     pven,                   &
        &         ktrac,                                                &
        &         zxtenh,   zxtu,     zxtd,     zmfdxt,                 &
        &         pgeoh,    paphp1,                                     &
        &         ptu,      pqu,      zuu,      zvu,                    &
        &         ldcum,    kcbot,    kctop,    zmfub,    zrfl,         &
        &         ztd,      zqd,      zud,      zvd,                    &
        &         pmfd,     zmfds,    zmfdq,    zdmfdp,                 &
        &         zcpcu,                                                &
        &         idtop,    loddraf)
      !
      !*            (B)  Determine downdraft t,q and fluxes in 'cuddraf'
      !                  -----------------------------------------------
      CALL cuddraf(jg,       jb,                                        &
        &          jcs,     kproma,   kbdim,    klev,     klevp1,       &
        &          pmref,                                               &
        &          ztenh,    zqenh,    puen,     pven,                  &
        &          ktrac,                                               &
        &          zxtenh,   zxtd,     zmfdxt,                          &
        &          pgeoh,    paphp1,   zrfl,                            &
        &          ztd,      zqd,      zud,      zvd,                   &
        &          pmfd,     zmfds,    zmfdq,    zdmfdp,                &
        &          zcpcu,                                               &
        &          loddraf)
      !
    END IF
    !
    !*    5.5      Recalculate cloud base massflux from a cape closure for deep
    !*             convection (ktype=1) and by pbl equilibrium taking downdrafts
    !*             into account for shallow convection (ktype=2)
    !              -------------------------------------------------------------
    ! mpuetz: cuasc can modify kctop !!
    !
    !$ACC PARALLEL DEFAULT(PRESENT)
    !$ACC LOOP GANG VECTOR
    DO jl = jcs, kproma
      zkcbot(jl) = REAL(kcbot(jl),wp)
    END DO
    !$ACC END PARALLEL
    !
    !$ACC PARALLEL DEFAULT(PRESENT)
    !$ACC LOOP GANG VECTOR
    DO jl=jcs,kproma
      zheat(jl)=0._wp
      zcape(jl)=0._wp
      zmfub1(jl)=zmfub(jl)
    ! zhelp(jl) =0._wp
    END DO
    !$ACC END PARALLEL
    !
    DO jk=1,klev
      zjk = REAL(jk,wp)
!IBM* ASSERT(NODEPS)
      !$ACC PARALLEL DEFAULT(PRESENT)
      !$ACC LOOP GANG VECTOR PRIVATE( jl, zkctop, zlo1 )
      DO nl = jcs,ldcnt
        jl = ldidx(nl)
        zkctop = REAL(kctop(jl),wp)
        zlo1 = FSEL(zkcbot(jl) - zjk,1._wp,0._wp)* FSEL(zkctop - zjk,0._wp,1._wp)
        ilo1(nl) = INT(zlo1)
      END DO
      !$ACC END PARALLEL
      locnt = 0
      !$ACC UPDATE HOST( ldidx, kcbot, kctop )
      DO nl = jcs,ldcnt
        jl = ldidx(nl)
        IF(jk.LE.kcbot(jl).AND.jk.GT.kctop(jl)) THEN
      !  IF (ilo1(nl).GT.0) THEN
          locnt = locnt + 1
          loidx(locnt) = jl
        END IF
      END DO
      !$ACC UPDATE DEVICE( loidx )
#ifdef __ibmdbg__
      PRINT '(A6,I3,2 I4)','cumas1',jk,ldcnt,locnt,INT(zkcbot(jl))
#endif
      ! mpuetz: there is reuse of zro, maybe we can calulate once and store
!IBM* ASSERT(NODEPS)
      !$ACC PARALLEL DEFAULT(PRESENT)
      !$ACC LOOP GANG VECTOR PRIVATE( jl, zroi, ztenhi, zdz )
      DO nl=1,locnt
        jl = loidx(nl)
        zroi = SWDIV_NOCHK(rd*ztenh(jl,jk)*(1._wp+vtmpc1*zqenh(jl,jk)),paphp1(jl,jk))
        ztenhi = SWDIV_NOCHK(1._wp,ztenh(jl,jk))
        zdz    = pmref(jl,jk-1)*zroi
        zheat(jl) = zheat(jl) + grav*((pten(jl,jk-1)-pten(jl,jk) + grav*zdz*zcpcui(jl,jk))*ztenhi  &
          &                           +vtmpc1*(pqen(jl,jk-1)-pqen(jl,jk)))                         &
          &                         *(pmfu(jl,jk)+pmfd(jl,jk))*zroi
        zcape(jl) = zcape(jl) + grav*((ptu(jl,jk)   -ztenh(jl,jk))*ztenhi                          &
          &                           +vtmpc1*(pqu(jl,jk)-zqenh(jl,jk))-plu(jl,jk))                &
          &                         *zdz
#ifdef __ibmdbg__
        PRINT '(A6,I3,I4,5 E18.10,I4)','zheat',jk,jl,zheat(jl),zcape(jl),zroi,ztenhi,zdz,kcbot(jl)
        PRINT '(A6,I3,I4,6 E18.10)','zcape',jk,jl,pten(jl,jk-1)-pten(jl,jk),pqen(jl,jk-1)-pqen(jl,jk),pmfu(jl,jk)+pmfd(jl,jk),ptu(jl,jk),pqu(jl,jk),plu(jl,jk)
#endif
      ENDDO
      !$ACC END PARALLEL
    ENDDO
    !
    !  DO jl=jcs,kproma
    !     DO jk=2,klev
    !        zro=paphp1(jl,jk)/(rd*ztenh(jl,jk)*                         &
    !             (1._wp+vtmpc1*zqenh(jl,jk)))
    !        zdz=pmref(jl,jk-1)/zro
    !        zhelp(jl)=zhelp(jl) +                               &
    !             (grav*(ptu(jl,jk)-ztenh(jl,jk))/ztenh(jl,jk)     &
    !             +grav*vtmpc1*(pqu(jl,jk)-zqenh(jl,jk))      &
    !             -grav*plu(jl,jk) ) * zdz
    !     ENDDO
    !     if (zhelp(jl).lt.0._wp) zhelp(jl)=0._wp
    !  ENDDO
    !$ACC PARALLEL DEFAULT(PRESENT)
    !$ACC LOOP GANG VECTOR PRIVATE( ikb, zmfmax )
    DO jl=jcs,kproma
      IF(ldcum(jl).AND.ktype(jl).EQ.1) THEN
        ikb=kcbot(jl)
        zmfub1(jl) = SWDIV_NOCHK((zcape(jl)*zmfub(jl)),(zheat(jl)*cmftau))
        zmfub1(jl) = MAX(zmfub1(jl),0.001_wp)
        zmfmax     = pmref(jl,ikb-1)*zcons
        zmfub1(jl) = MIN(zmfub1(jl),zmfmax)
      ENDIF
    ENDDO
    !$ACC END PARALLEL
    !
    !*      Recalculate convective fluxes due to effect of downdrafts on boundary
    !*      layer moisture budget for shallow convection (ktype=2)
    !       ---------------------------------------------------------------------
    ldcnt = jcs-1
    !$ACC UPDATE HOST( ktype )
    DO jl=jcs,kproma
      llo1=ktype(jl).EQ.2
      IF(llo1) THEN
        ldcnt = ldcnt + 1
        ldidx(ldcnt) = jl
      ENDIF
    ENDDO
    !$ACC UPDATE DEVICE( ldidx )
!DIR$ IVDEP
!IBM* ASSERT(NODEPS)
    !$ACC PARALLEL DEFAULT(PRESENT)
    !$ACC LOOP GANG VECTOR PRIVATE( jl, ikb, llo1, zeps, zqumqe, zdqmin, zmfmax )
    DO nl=jcs,ldcnt
      jl = ldidx(nl)
      ikb=kcbot(jl)
      llo1=pmfd(jl,ikb).LT.0._wp.AND.loddraf(jl)
      zeps=MERGE(cmfdeps,0._wp,llo1)
      zqumqe=pqu(jl,ikb)+plu(jl,ikb)-zeps*zqd(jl,ikb)-(1._wp-zeps)*zqenh(jl,ikb)
      zdqmin=MAX(0.01_wp*zqenh(jl,ikb),1.e-10_wp)
      zmfmax=pmref(jl,ikb-1)*zcons
      llo1=zdqpbl(jl).GT.0._wp.AND.zqumqe.GT.zdqmin.AND.ldcum(jl)                    &
        &                                          .AND.zmfub(jl).LT.zmfmax
      zmfub1(jl)=MERGE(zdqpbl(jl)/MAX(zqumqe,zdqmin),zmfub(jl),llo1)
      zmfub1(jl)=MERGE(zmfub1(jl),zmfub(jl),                                         &
        &                          ABS(zmfub1(jl)-zmfub(jl)).LT.0.2_wp*zmfub(jl))
#ifdef __ibmdbg__
      PRINT '(A6,I4,E18.10)','zmfub1',jl,zmfub1(jl)
#endif
    END DO
    !$ACC END PARALLEL
    ldcnt = jcs-1
    !$ACC UPDATE HOST( ldcum )
    DO jl=jcs,kproma
      IF (ldcum(jl)) THEN
        ldcnt = ldcnt + 1
        ldidx(ldcnt) = jl
      ENDIF
    ENDDO
    !$ACC UPDATE DEVICE( ldidx )
    DO jk=1,klev
!IBM* ASSERT(NODEPS)
      !$ACC PARALLEL DEFAULT(PRESENT)
      !$ACC LOOP GANG VECTOR PRIVATE( jl, zfac )
      DO nl=jcs,ldcnt
        jl = ldidx(nl)
        ztmp1(nl) = SWDIV_NOCHK(zmfub1(jl),MAX(zmfub(jl),1.e-10_wp))
        zfac      = ztmp1(nl)
        pmfd(jl,jk)   = pmfd(jl,jk)*zfac
        zmfds(jl,jk)  = zmfds(jl,jk)*zfac
        zmfdq(jl,jk)  = zmfdq(jl,jk)*zfac
        zdmfdp(jl,jk) = zdmfdp(jl,jk)*zfac
#ifdef __ibmdbg__
        PRINT '(A6,I3,I4,E18.10)','zfac',jk,jl,ztmp1(nl)
#endif
      END DO
      !$ACC END PARALLEL
!IBM* unroll(4)
      !$ACC PARALLEL DEFAULT(PRESENT) IF( ktrac > 0 )
      !$ACC LOOP SEQ
      DO jt=1,ktrac
!IBM* ASSERT(NODEPS)
        !$ACC LOOP GANG VECTOR PRIVATE( jl, zfac )
        DO nl=jcs,ldcnt
          jl = ldidx(nl)
          zfac             = ztmp1(nl)
          zmfdxt(jl,jk,jt) = zmfdxt(jl,jk,jt)*zfac
        END DO
      END DO
      !$ACC END PARALLEL
      !
    END DO
    !
    !*       New values of cloud base mass flux
    !        ----------------------------------
    !
    !$ACC PARALLEL DEFAULT(PRESENT)
    !$ACC LOOP GANG VECTOR PRIVATE( jl )
    DO nl=jcs,ldcnt
      jl = ldidx(nl)
      zmfub(jl) = zmfub1(jl)
    END DO
    !$ACC END PARALLEL
    !
    !-----------------------------------------------------------------------
    !*    6.0   Determine final cloud ascent for entraining plume
    !*          for penetrative convection (type=1),
    !*          for shallow to medium convection (type=2)
    !*          and for mid-level convection (type=3).
    !           --------------------------------------------------
    !
    CALL cuasc(jg,     jb,                                               &
      &        jcs, kproma, kbdim, klev, klevp1, klevm1,                 &
      &        pzf,      pzh,      pmref,                                &
      &        ztenh,    zqenh,    puen,     pven,                       &
      &        ktrac,                                                    &
      &        pdtime,                                                   &
      &        zxtenh,   pxten,    zxtu,     zmfuxt,                     &
      &        pten,     pqen,     zqsen,                                &
      &        pgeo,     pgeoh,    paphp1,   pthvsig,                    &
      &        pqte,               pverv,    ilwmin,                     &
      &        ldcum,    ldland,   ktype,    ilab,                       &
      &        ptu,      pqu,      plu,      zuu,      zvu,              &
      &        pmfu,     zmfub,    zentr,                                &
      &        zmfus,    zmfuq,                                          &
      &        zmful,    plude,    pqude,    zdmfup,                     &
      &        ihmin,    zhhatt,   zhcbase,  zqsenh,                     &
      &        zcpen,    zcpcu,                                          &
      &        kcbot,    kctop,    ictop0                                &
      &       )

    !-----------------------------------------------------------------------
    !*    7.0      Determine final convective fluxes in 'cuflx'
    !              --------------------------------------------
    !
    CALL cuflx(jcs,      kproma,   kbdim,    klev,     klevp1,           &
      &        pmref,                                                    &
      &        pqen,     zqsen,    ztenh,    zqenh,                      &
      &        ktrac,                                                    &
      &        pdtime,                                                   &
      &        zevapcu,                                                  &
      &        zxtenh,   zmfuxt,   zmfdxt,                               &
      &        paphp1,   pgeoh,                                          &
      &        kcbot,    kctop,    idtop,                                &
      &        ktype,    loddraf,  ldcum,                                &
      &        pmfu,     pmfd,     zmfus,    zmfds,                      &
      &        zmfuq,    zmfdq,    zmful,                                &
      &        zdmfup,   zdmfdp,   zrfl,                                 &
      &        zcpcu,                                                    &
      &        pten,     zsfl,     zdpmel,   itopm2                     )

    !$ACC PARALLEL DEFAULT(PRESENT)
    !$ACC LOOP GANG VECTOR
    DO jl = 1, kbdim
      prsfc(jl) = zrfl(jl)
      pssfc(jl) = zsfl(jl)
    END DO
    !$ACC END PARALLEL

    !-----------------------------------------------------------------------
    !
    !*    8.0      Update tendencies for t and q in subroutine'cudtdq'
    !              ---------------------------------------------------
    !
    CALL cudtdq(jcs, kproma, kbdim, klev, itopm2, ldcum, ktrac,          &
      &         pmref,    pten,                                          &
      &         zmfuxt,   zmfdxt,                                        &
      &         zmfus,    zmfds,    zmfuq,    zmfdq,                     &
      &         zmful,    zdmfup,   zdmfdp,   plude,                     &
      &         zdpmel,                                                  &
      &         zalvsh,                                                  &
      &         pcon_dtrl,pcon_dtri,pcon_iqte,                           &
      &         pq_cnv,   pqte_cnv, pxtte_cnv,                           &
      &         pxtecl,   pxteci                                        )
    !
    !-----------------------------------------------------------------------
    !
    !*    9.0      Update tendencies for u and v in subroutine'cududv'
    !              ---------------------------------------------------
    !
    IF(lmfdudv) THEN
      CALL cududv(jcs,      kproma,   kbdim,    klev,     klevp1,       &
        &         itopm2,   ktype,    kcbot,    paphp1,   ldcum,        &
        &         pmref,    puen,     pven,     pvom_cnv, pvol_cnv,     &
        &         zuu,      zud,      zvu,      zvd,                    &
        &         pmfu,     pmfd)
      !
    END IF
    !
    !-----------------------------------------------------------------------
    !
    !*   10.0      Pressure altitude of convective cloud tops
    !              ------------------------------------------
    !
    !$ACC PARALLEL DEFAULT(PRESENT)
    !$ACC LOOP GANG VECTOR
    DO jl=jcs,kproma
      itopec2(jl)=klevp1
    END DO
    !$ACC END PARALLEL
    !
    !$ACC PARALLEL DEFAULT(PRESENT)
    !$ACC LOOP SEQ
    DO jk=1,klev-4
      !$ACC LOOP GANG VECTOR
      DO jl=jcs,kproma
        IF(ilab(jl,jk).EQ.2 .AND. itopec2(jl).EQ.klevp1) THEN
          itopec2(jl)=jk
        END IF
      END DO
    END DO
    !$ACC END PARALLEL
    !
    !$ACC PARALLEL DEFAULT(PRESENT)
    !$ACC LOOP GANG VECTOR
    DO jl = 1, kbdim
      ptop(jl)=99999._wp
    END DO
    !$ACC END PARALLEL
    !
    !$ACC PARALLEL DEFAULT(PRESENT)
    !$ACC LOOP GANG VECTOR
    DO jl=jcs,kproma
      IF(itopec2(jl).EQ.1) THEN
        ptop(jl)=(paphp1(jl,1)+paphp1(jl,2))*0.5_wp
      ELSE IF(itopec2(jl).NE.klevp1) THEN
        ptop(jl)=paphp1(jl,itopec2(jl))
      ELSE
        ptop(jl)=99999._wp
      END IF
    END DO
    !$ACC END PARALLEL
    !
    !$ACC END DATA
    !$ACC END DATA
    !$ACC END DATA
    !$ACC END DATA
    !
    !-----------------------------------------------------------------------
    !
  END SUBROUTINE cumastr

END MODULE mo_cumastr
