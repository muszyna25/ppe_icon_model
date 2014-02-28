MODULE mo_cuinitialize
#ifdef __xlC__
@PROCESS HOT
#endif
 
USE mo_kind,               ONLY: wp
USE mo_physical_constants, ONLY: rd, cpd, vtmpc1, alv, als, tmelt
USE mo_echam_conv_constants,ONLY: lmfdudv, cbfac, cminbuoy, cmaxbuoy
USE mo_cuadjust,           ONLY: cuadjtq

IMPLICIT NONE

PRIVATE

PUBLIC :: cuini, cubase

CONTAINS 

SUBROUTINE cuini(kproma, kbdim, klev, klevp1, klevm1,                  &
           pten,     pqen,     pqsen,    pxen,     puen,     pven,     &
           ptven,    ktrac,                                            &
           pxten,    pxtenh,   pxtu,     pxtd,     pmfuxt,   pmfdxt,   &
           pverv,    pgeo,     paphp1,   pgeoh,                        &
           ptenh,    pqenh,    pqsenh,   pxenh,    klwmin,             &
           ptu,      pqu,      ptd,      pqd,                          &
           puu,      pvu,      pud,      pvd,                          &
           pmfu,     pmfd,     pmfus,    pmfds,                        &
           pmfuq,    pmfdq,    pdmfup,   pdmfdp,                       &
           pcpen,    pcpcu,    palvsh,                                 &
           pdpmel,   plu,      plude,    pqude,    klab                )
!
!          M.TIEDTKE         E.C.M.W.F.     12/89
!
!          PURPOSE
!          -------
!
!          THIS ROUTINE INTERPOLATES LARGE-SCALE FIELDS OF T,Q ETC.
!          TO HALF LEVELS (I.E. GRID FOR MASSFLUX SCHEME),
!          DETERMINES LEVEL OF MAXIMUM VERTICAL VELOCITY
!          AND INITIALIZES VALUES FOR UPDRAFTS AND DOWNDRAFTS
!
!          INTERFACE
!          ---------
!          THIS ROUTINE IS CALLED FROM *CUMASTR*.
!
!          METHOD.
!          --------
!          FOR EXTRAPOLATION TO HALF LEVELS SEE TIEDTKE(1989)
!
!          EXTERNALS
!          ---------
!          *CUADJTQ* TO SPECIFY QS AT HALF LEVELS
!

INTEGER, INTENT (IN) :: kproma, kbdim, klev, klevp1, klevm1, ktrac

REAL(wp):: pten(kbdim,klev),          pqen(kbdim,klev),                &
           puen(kbdim,klev),          pven(kbdim,klev),                &
           pqsen(kbdim,klev),         pverv(kbdim,klev),               &
           pgeo(kbdim,klev),          pgeoh(kbdim,klev),               &
           paphp1(kbdim,klevp1),      ptenh(kbdim,klev),               &
           pxenh(kbdim,klev),         pxen(kbdim,klev),                &
           ptven(kbdim,klev),         palvsh(kbdim,klev),              &
           pqenh(kbdim,klev),         pqsenh(kbdim,klev)
REAL(wp):: pcpen(kbdim,klev),         pcpcu(kbdim,klev)
!
REAL(wp):: ptu(kbdim,klev),           pqu(kbdim,klev),                 &
           ptd(kbdim,klev),           pqd(kbdim,klev),                 &
           puu(kbdim,klev),           pud(kbdim,klev),                 &
           pvu(kbdim,klev),           pvd(kbdim,klev),                 &
           pmfu(kbdim,klev),          pmfd(kbdim,klev),                &
           pmfus(kbdim,klev),         pmfds(kbdim,klev),               &
           pmfuq(kbdim,klev),         pmfdq(kbdim,klev),               &
           pdmfup(kbdim,klev),        pdmfdp(kbdim,klev),              &
           plu(kbdim,klev),           plude(kbdim,klev),               &
           pqude(kbdim,klev)
REAL(wp):: pdpmel(kbdim,klev)
INTEGER :: klab(kbdim,klev),          klwmin(kbdim)
!
REAL(wp):: zwmax(kbdim)
REAL(wp):: zph(kbdim)
INTEGER :: loidx(kbdim)
REAL(wp):: pxten(kbdim,klev,ktrac),   pxtenh(kbdim,klev,ktrac),        &
           pxtu(kbdim,klev,ktrac),    pxtd(kbdim,klev,ktrac),          &
           pmfuxt(kbdim,klev,ktrac),  pmfdxt(kbdim,klev,ktrac)

INTEGER :: jk, jl, jt, ik, icall
REAL(wp):: zarg, zcpm, zzs
LOGICAL :: llo1

!
!  INTRINSIC FUNCTIONS
INTRINSIC MAX, MIN
!
!----------------------------------------------------------------------
!
!*    1.           SPECIFY LARGE SCALE PARAMETERS AT HALF LEVELS
!*                 ADJUST TEMPERATURE FIELDS IF STATICLY UNSTABLE
!*                 FIND LEVEL OF MAXIMUM VERTICAL VELOCITY
!                  ----------------------------------------------
!
  DO 105 jl=1,kproma
     zarg=paphp1(jl,klevp1)/paphp1(jl,klev)
     pgeoh(jl,klev)=rd*ptven(jl,klev)*LOG(zarg)
105 END DO
  DO 107 jk=klevm1,2,-1
     DO 106 jl=1,kproma
        zarg=paphp1(jl,jk+1)/paphp1(jl,jk)
        pgeoh(jl,jk)=pgeoh(jl,jk+1)+rd*ptven(jl,jk)*LOG(zarg)
106  END DO
107 END DO
  DO 130 jk=2,klev
!IBM* NOVECTOR
     DO 110 jl=1,kproma
        zcpm=(pcpen(jl,jk)+pcpen(jl,jk-1))*0.5_wp
        pcpcu(jl,jk)=zcpm
        ptenh(jl,jk)=(MAX(pcpen(jl,jk-1)*pten(jl,jk-1)+pgeo(jl,jk-1),  &
                   pcpen(jl,jk)*pten(jl,jk)+pgeo(jl,jk))-pgeoh(jl,jk))
        ptenh(jl,jk) = ptenh(jl,jk)/zcpm
        pqsenh(jl,jk)=pqsen(jl,jk-1)
        zph(jl)=paphp1(jl,jk)
        loidx(jl)=jl
110  END DO
!
     DO 1104 jt=1,ktrac
        DO 1102 jl=1,kproma
           pxtenh(jl,jk,jt)=(pxten(jl,jk,jt)+pxten(jl,jk-1,jt))*0.5_wp
1102    END DO
1104 END DO
!
!
     ik=jk
     icall=0
     CALL cuadjtq(kproma, kbdim, klev, ik,                         &
          zph,      ptenh,    pqsenh,   loidx, kproma, icall)
!
     DO 120 jl=1,kproma
        pxenh(jl,jk)=(pxen(jl,jk)+pxen(jl,jk-1))*0.5_wp
        pqenh(jl,jk)=MIN(pqen(jl,jk-1),pqsen(jl,jk-1))                 &
                          +(pqsenh(jl,jk)-pqsen(jl,jk-1))
        pqenh(jl,jk)=MAX(pqenh(jl,jk),0._wp)
120  END DO
130 END DO
!
!IBM* NOVECTOR
  DO 140 jl=1,kproma
     ptenh(jl,klev)=pcpen(jl,klev)*pten(jl,klev)+pgeo(jl,klev)-pgeoh(jl,klev)
     ptenh(jl,klev)=ptenh(jl,klev)/pcpen(jl,klev)
     pxenh(jl,klev)=pxen(jl,klev)
     pqenh(jl,klev)=pqen(jl,klev)
     pcpcu(jl,1)=pcpen(jl,1)
     ptenh(jl,1)=pten(jl,1)
     pxenh(jl,1)=pxen(jl,1)
     pqenh(jl,1)=pqen(jl,1)
     pgeoh(jl,1)=pgeo(jl,1)
     klwmin(jl)=klev
     zwmax(jl)=0._wp
140 END DO
!
  DO 1404 jt=1,ktrac
     DO 1402 jl=1,kproma
        pxtenh(jl,klev,jt)=pxten(jl,klev,jt)
        pxtenh(jl,1,jt)=pxten(jl,1,jt)
1402 END DO
1404 END DO
!
!
  DO 160 jk=klevm1,2,-1
!IBM* NOVECTOR
     DO 150 jl=1,kproma
        zzs=MAX(pcpcu(jl,jk)*ptenh(jl,jk)+pgeoh(jl,jk),                &
                      pcpcu(jl,jk+1)*ptenh(jl,jk+1)+pgeoh(jl,jk+1))
        ptenh(jl,jk)=(zzs-pgeoh(jl,jk))/pcpcu(jl,jk)
150  END DO
160 END DO
!
  DO 170 jk=1,klev
! NOVECTOR ?
     DO 165 jl=1,kproma
        llo1 = (ptenh(jl,jk)-tmelt) .GT. 0.0_wp
        palvsh(jl,jk) = MERGE(alv,als,llo1)
165 END DO
170 END DO
!
  DO 190 jk=klev,3,-1
!DIR$ IVDEP
!OCL NOVREC
     DO 180 jl=1,kproma
        IF(pverv(jl,jk).LT.zwmax(jl)) THEN
           zwmax(jl)=pverv(jl,jk)
           klwmin(jl)=jk
        END IF
180  END DO
190 END DO
!
!
!-----------------------------------------------------------------------
!*    2.0          INITIALIZE VALUES FOR UPDRAFTS AND DOWNDRAFTS
!*                 ---------------------------------------------
!
  DO 230 jk=1,klev
     ik=jk-1
     IF(jk.EQ.1) ik=1
     DO 220 jl=1,kproma
        ptu(jl,jk)=ptenh(jl,jk)
        ptd(jl,jk)=ptenh(jl,jk)
        pqu(jl,jk)=pqenh(jl,jk)
        pqd(jl,jk)=pqenh(jl,jk)
        plu(jl,jk)=0._wp
        puu(jl,jk)=puen(jl,ik)
        pud(jl,jk)=puen(jl,ik)
        pvu(jl,jk)=pven(jl,ik)
        pvd(jl,jk)=pven(jl,ik)
        pmfu(jl,jk)=0._wp
        pmfd(jl,jk)=0._wp
        pmfus(jl,jk)=0._wp
        pmfds(jl,jk)=0._wp
        pmfuq(jl,jk)=0._wp
        pmfdq(jl,jk)=0._wp
        pdmfup(jl,jk)=0._wp
        pdmfdp(jl,jk)=0._wp
        pdpmel(jl,jk)=0._wp
        plude(jl,jk)=0._wp
        pqude(jl,jk)=0._wp
        klab(jl,jk)=0
220  END DO
!
     DO 2204 jt=1,ktrac
        DO 2202 jl=1,kproma
           pxtu(jl,jk,jt)=pxtenh(jl,jk,jt)
           pxtd(jl,jk,jt)=pxtenh(jl,jk,jt)
           pmfuxt(jl,jk,jt)=0._wp
           pmfdxt(jl,jk,jt)=0._wp
2202    END DO
2204 END DO
!
230 END DO
!
  RETURN
END SUBROUTINE cuini

SUBROUTINE cubase(   kproma, kbdim, klev, klevp1, klevm1,              &
           ptenh,    pqenh,    pgeoh,    paph,   pthvsig,              &
           ptu,      pqu,      plu,                                    &
           puen,     pven,     puu,      pvu,                          &
           pcpcu,                                                      &
           ldcum,    kcbot,    klab)
!
!          THIS ROUTINE CALCULATES CLOUD BASE VALUES (T AND Q)
!          FOR CUMULUS PARAMETERIZATION
!
!          M.TIEDTKE         E.C.M.W.F.     7/86 MODIF. 12/89
!
!          PURPOSE.
!          --------
!          TO PRODUCE CLOUD BASE VALUES FOR CU-PARAMETRIZATION
!
!          INTERFACE
!          ---------
!          THIS ROUTINE IS CALLED FROM *CUMASTR*.
!          INPUT ARE ENVIRONM. VALUES OF T,Q,P,PHI AT HALF LEVELS.
!          IT RETURNS CLOUD BASE VALUES AND FLAGS AS FOLLOWS;
!                 KLAB=1 FOR SUBCLOUD LEVELS
!                 KLAB=2 FOR CONDENSATION LEVEL
!
!          METHOD.
!          --------
!          LIFT SURFACE AIR DRY-ADIABATICALLY TO CLOUD BASE
!          (NON ENTRAINING PLUME,I.E.CONSTANT MASSFLUX)
!
!          EXTERNALS
!          ---------
!          *CUADJTQ* FOR ADJUSTING T AND Q DUE TO CONDENSATION IN ASCENT
!

INTEGER, INTENT (IN) :: kproma, kbdim, klev, klevp1, klevm1

REAL(wp):: ptenh(kbdim,klev),       pqenh(kbdim,klev),                 &
           pgeoh(kbdim,klev),       paph(kbdim,klevp1),                &
           pthvsig(kbdim)
!
REAL(wp):: ptu(kbdim,klev),         pqu(kbdim,klev),                   &
           plu(kbdim,klev)
REAL(wp):: puen(kbdim,klev),        pven(kbdim,klev),                  &
           puu(kbdim,klev),         pvu(kbdim,klev)
REAL(wp):: pcpcu(kbdim,klev)
INTEGER :: klab(kbdim,klev),        kcbot(kbdim)
LOGICAL :: ldcum(kbdim)
!
REAL(wp):: zqold(kbdim)
REAL(wp):: zph(kbdim)
INTEGER :: loidx(kbdim)

INTEGER :: jl, jk, nl, is, ik, ikb, icall
REAL(wp):: zbuo, zz, zlift
!
!
!----------------------------------------------------------------------
!
!     1.           INITIALIZE VALUES AT LIFTING LEVEL
!                  ----------------------------------
!
  DO 110 jl=1,kproma
     klab(jl,klev)=1
     kcbot(jl)=klevm1
     ldcum(jl)=.FALSE.
     puu(jl,klev)=puen(jl,klev)*(paph(jl,klevp1)-paph(jl,klev))
     pvu(jl,klev)=pven(jl,klev)*(paph(jl,klevp1)-paph(jl,klev))
110 END DO
!
!
!----------------------------------------------------------------------
!
!     2.0          DO ASCENT IN SUBCLOUD LAYER,
!                  CHECK FOR EXISTENCE OF CONDENSATION LEVEL,
!                  ADJUST T,Q AND L ACCORDINGLY IN *CUADJTQ*,
!                  CHECK FOR BUOYANCY AND SET FLAGS
!                  -------------------------------------
!
  DO 290 jk=klevm1,2,-1
     is=0
     DO 210 jl=1,kproma
        IF (klab(jl,jk+1).EQ.1) THEN
           is = is + 1
           loidx(is) = jl
        END IF
        zph(jl)=paph(jl,jk)
210  END DO
     IF(is.EQ.0) CYCLE !GOTO 290
!IBM* ASSERT(NODEPS)
     DO 220 nl=1,is
        jl = loidx(nl)
        zlift=MAX(cminbuoy,MIN(cmaxbuoy,pthvsig(jl)*cbfac))
        zlift=MIN(zlift,1.0_wp)
        pqu(jl,jk)=pqu(jl,jk+1)
        ptu(jl,jk)=(pcpcu(jl,jk+1)*ptu(jl,jk+1)+pgeoh(jl,jk+1)      &
             -pgeoh(jl,jk))/pcpcu(jl,jk)
        zbuo=ptu(jl,jk)*(1._wp+vtmpc1*pqu(jl,jk))-ptenh(jl,jk)      &
             *(1._wp+vtmpc1*pqenh(jl,jk))+zlift
        IF(zbuo.GT.0._wp) klab(jl,jk)=1
        zqold(jl)=pqu(jl,jk)
220  END DO
!
     ik=jk
     icall=1
     CALL cuadjtq(kproma, kbdim, klev, ik,                             &
                      zph,    ptu,   pqu,  loidx,   is,   icall)
!
!DIR$ IVDEP
!OCL NOVREC
     DO 240 nl=1,is
        jl = loidx(nl)
        IF(pqu(jl,jk).LT.zqold(jl)) THEN
           klab(jl,jk)=2
           zlift=MAX(cminbuoy,MIN(cmaxbuoy,pthvsig(jl)*cbfac))
           zlift=MIN(zlift,1.0_wp)
           plu(jl,jk)=plu(jl,jk)+zqold(jl)-pqu(jl,jk)
           zbuo=ptu(jl,jk)*(1._wp+vtmpc1*pqu(jl,jk)-plu(jl,jk))-       &
                        ptenh(jl,jk)*(1._wp+vtmpc1*pqenh(jl,jk))+zlift
           IF(zbuo.GT.0.) THEN
              kcbot(jl)=jk
              ldcum(jl)=.TRUE.
           END IF
        END IF
240  END DO
!
!             CALCULATE AVERAGES OF U AND V FOR SUBCLOUD ARA,.
!             THE VALUES WILL BE USED TO DEFINE CLOUD BASE VALUES.
!
     IF(lmfdudv) THEN
        DO 250 jl=1,kproma
           IF(jk.GE.kcbot(jl)) THEN
              puu(jl,klev)=puu(jl,klev)+                               &
                             puen(jl,jk)*(paph(jl,jk+1)-paph(jl,jk))
              pvu(jl,klev)=pvu(jl,klev)+                               &
                             pven(jl,jk)*(paph(jl,jk+1)-paph(jl,jk))
           END IF
250     END DO
     END IF
!
290 END DO
!
!
  IF(lmfdudv) THEN
     DO 310 jl=1,kproma
        IF(ldcum(jl)) THEN
           ikb=kcbot(jl)
           zz=1._wp/(paph(jl,klevp1)-paph(jl,ikb))
           puu(jl,klev)=puu(jl,klev)*zz
           pvu(jl,klev)=pvu(jl,klev)*zz
        ELSE
           puu(jl,klev)=puen(jl,klevm1)
           pvu(jl,klev)=pven(jl,klevm1)
        END IF
310  END DO
  END IF
!
  RETURN
END SUBROUTINE cubase
END MODULE mo_cuinitialize
