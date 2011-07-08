!>
!! Module contains subroutine cumastr -- the master routine for cumulus convection.
!!
!! @author M. Tiedtke, ECMWF
!!
!! @par Revision History
!! - M. Tiedtke, ECMWF, 1986/1987/1989
!! - Taken from ECHAM6, wrapped in module and modified for ICON
!!   by Hui Wan, MPI-M (2010-08)
!!
!! @par Copyright
!! 2002-2010 by DWD and MPI-M
!! This software is provided for non-commercial use only.
!! See the LICENSE and the WARRANTY conditions.
!!
!! @par License
!! The use of ICON is hereby granted free of charge for an unlimited time,
!! provided the following rules are accepted and applied:
!! <ol>
!! <li> You may use or modify this code for your own non commercial and non
!!    violent purposes.
!! <li> The code may not be re-distributed without the consent of the authors.
!! <li> The copyright notice and statement of authorship must appear in all
!!    copies.
!! <li> You accept the warranty conditions (see WARRANTY).
!! <li> In case you intend to use the code commercially, we oblige you to sign
!!    an according license agreement with DWD and MPI-M.
!! </ol>
!!
!! @par Warranty
!! This code has been tested up to a certain level. Defects and weaknesses,
!! which may be included in the code, do not establish any warranties by the
!! authors.
!! The authors do not make any warranty, express or implied, or assume any
!! liability or responsibility for the use, acquisition or application of this
!! software.
!!
#ifndef __xlC__
#define FSEL(a,b,c) MERGE(b,c,(a) >= 0._dp)
#define SWDIV_NOCHK(a,b) ((a)/(b))
#endif

MODULE mo_cumastr

  USE mo_kind,                ONLY: dp

#ifdef __ICON__
  USE mo_physical_constants,  ONLY: g=>grav, alv, als, tmelt, vtmpc1, rd
  USE mo_echam_conv_nml,      ONLY: entrpen, entrscv, cmfdeps
#else
  USE mo_control,             ONLY: nn
  USE mo_constants,           ONLY: g, alv, als, tmelt, vtmpc1, rd
  USE mo_cumulus_flux,        ONLY: entrpen, entrscv, cmfdeps
  USE mo_tracer_processes,    ONLY: xt_conv_massfix
#endif

  USE mo_cuini,               ONLY: cuini
  USE mo_cuasc,               ONLY: cuasc
  USE mo_cubase,              ONLY: cubase
  USE mo_cuflx,               ONLY: cuflx
  USE mo_cudtdq,              ONLY: cudtdq
  USE mo_cududv,              ONLY: cududv
  USE mo_cudlfs,              ONLY: cudlfs
  USE mo_cuddraf,             ONLY: cuddraf

#ifdef __ibmspline__
  USE mo_convect_tables, ONLY: prepare_ua_index_spline, lookup_ua_spline, lookup_ubc
#else
  USE mo_convect_tables, ONLY: prepare_ua_index, lookup_ua, lookup_ubc
#endif


  IMPLICIT NONE
  PRIVATE
  PUBLIC :: cumastr

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

CONTAINS
  !>
  !!
  SUBROUTINE cumastr(  ncvmicro, lmfdudv, lmfdd, lmfmid, dlev, cmftau,    &
                       cmfctop, cprcon, cminbuoy, &
                       pdtime, ptime_step_len,                            &
                       kproma, kbdim, klev, klevp1, klevm1, ilab,         &
!!$                       krow,                                              &
                       papp1,                                             &
                       pten,     pqen,     pxen,     puen,     pven,      &
                       ptven,    ktrac,    ldland,                        &
                       pxten,    pxtu,                                    &
!!$                       pxtte,                                             &
                       pverv,    pqsen,    pqhfla,                        &
                       paphp1,   pgeo,                                    &
                       ptte,     pqte,     pvom,     pvol,                &
                       prsfc,    pssfc,    paprc,    paprs,    pxtec,     &
                       pqtec,    pqude,                                   &
                       ldcum,    ktype,    kcbot,    kctop,               &
                       ptu,      pqu,      plu,      plude,               &
                       pmfu,     pmfd,     prain,    pthvsig,             &
                       pcvcbot,  pwcape,                                  &! for CDNC/IC
                       pxtecl,   pxteci,   pxtecnl,  pxtecni,             &! for CDNC/IC
!!$                       ptkem1,                                            &
                       ptte_cnv, pvom_cnv, pvol_cnv, pqte_cnv, pxtte_cnv )
!
!**** *CUMASTR*  MASTER ROUTINE FOR CUMULUS MASSFLUX-SCHEME
!
!     M.TIEDTKE      E.C.M.W.F.     1986/1987/1989
!
!     PURPOSE
!     -------
!
!          THIS ROUTINE COMPUTES THE PHYSICAL TENDENCIES OF THE
!     PROGNOSTIC VARIABLES T,Q,U AND V DUE TO CONVECTIVE PROCESSES.
!     PROCESSES CONSIDERED ARE: CONVECTIVE FLUXES, FORMATION OF
!     PRECIPITATION, EVAPORATION OF FALLING RAIN BELOW CLOUD BASE,
!     SATURATED CUMULUS DOWNDRAFTS.
!
!**   INTERFACE.
!     ----------
!
!          *CUMASTR* IS CALLED FROM *CUCALL*
!     THE ROUTINE TAKES ITS INPUT FROM THE LONG-TERM STORAGE
!     T,Q,U,V,PHI AND P AND MOISTURE TENDENCIES.
!     IT RETURNS ITS OUTPUT TO THE SAME SPACE
!      1.MODIFIED TENDENCIES OF MODEL VARIABLES
!      2.RATES OF CONVECTIVE PRECIPITATION
!        (USED IN SUBROUTINE SURF)
!
!     METHOD
!     -------
!
!     PARAMETERIZATION IS DONE USING A MASSFLUX-SCHEME.
!        (1) DEFINE CONSTANTS AND PARAMETERS
!        (2) SPECIFY VALUES (T,Q,QS...) AT HALF LEVELS AND
!            INITIALIZE UPDRAFT- AND DOWNDRAFT-VALUES IN 'CUINI'
!        (3) CALCULATE CLOUD BASE IN 'CUBASE'
!            AND SPECIFY CLOUD BASE MASSFLUX FROM PBL MOISTURE BUDGET
!        (4) DO CLOUD ASCENT IN 'CUASC' IN ABSENCE OF DOWNDRAFTS
!        (5) DO DOWNDRAFT CALCULATIONS:
!              (A) DETERMINE VALUES AT LFS IN 'CUDLFS'
!              (B) DETERMINE MOIST DESCENT IN 'CUDDRAF'
!              (C) RECALCULATE CLOUD BASE MASSFLUX CONSIDERING THE
!                  EFFECT OF CU-DOWNDRAFTS
!        (6) DO FINAL CLOUD ASCENT IN 'CUASC'
!        (7) DO FINAL ADJUSMENTS TO CONVECTIVE FLUXES IN 'CUFLX',
!            DO EVAPORATION IN SUBCLOUD LAYER
!        (8) CALCULATE INCREMENTS OF T AND Q IN 'CUDTDQ'
!        (9) CALCULATE INCREMENTS OF U AND V IN 'CUDUDV'
!
!     EXTERNALS.
!     ----------
!
!       CUINI:  INITIALIZES VALUES AT VERTICAL GRID USED IN CU-PARAMETR.
!       CUBASE: CLOUD BASE CALCULATION FOR PENETR.AND SHALLOW CONVECTION
!       CUASC:  CLOUD ASCENT FOR ENTRAINING PLUME
!       CUDLFS: DETERMINES VALUES AT LFS FOR DOWNDRAFTS
!       CUDDRAF:DOES MOIST DESCENT FOR CUMULUS DOWNDRAFTS
!       CUFLX:  FINAL ADJUSTMENTS TO CONVECTIVE FLUXES (ALSO IN PBL)
!       CUDQDT: UPDATES TENDENCIES FOR T AND Q
!       CUDUDV: UPDATES TENDENCIES FOR U AND V
!
!     SWITCHES.
!     --------
!
!          LMFPEN=.T.   PENETRATIVE CONVECTION IS SWITCHED ON
!          LMFSCV=.T.   SHALLOW CONVECTION IS SWITCHED ON
!          LMFMID=.T.   MIDLEVEL CONVECTION IS SWITCHED ON
!          LMFDD=.T.    CUMULUS DOWNDRAFTS SWITCHED ON
!          LMFDUDV=.T.  CUMULUS FRICTION SWITCHED ON
!
!     MODEL PARAMETERS (DEFINED IN SUBROUTINE CUPARAM)
!     ------------------------------------------------
!     ENTRPEN    ENTRAINMENT RATE FOR PENETRATIVE CONVECTION
!     ENTRSCV    ENTRAINMENT RATE FOR SHALLOW CONVECTION
!     ENTRMID    ENTRAINMENT RATE FOR MIDLEVEL CONVECTION
!     ENTRDD     ENTRAINMENT RATE FOR CUMULUS DOWNDRAFTS
!     CMFCTOP    RELATIVE CLOUD MASSFLUX AT LEVEL ABOVE NONBUOYANCY LEVEL
!     CMFCMAX    MAXIMUM MASSFLUX VALUE ALLOWED FOR
!     CMFCMIN    MINIMUM MASSFLUX VALUE (FOR SAFETY)
!     CMFDEPS    FRACTIONAL MASSFLUX FOR DOWNDRAFTS AT LFS
!     CPRCON     COEFFICIENT FOR CONVERSION FROM CLOUD WATER TO RAIN
!
!     REFERENCE.
!     ----------
!
!          PAPER ON MASSFLUX SCHEME (TIEDTKE,1989)
!
!
INTEGER, INTENT (IN) :: ncvmicro
LOGICAL, INTENT (IN) :: lmfdudv, lmfdd, lmfmid 
REAL(dp),INTENT (IN) :: dlev, cmftau, cmfctop, cprcon, cminbuoy
INTEGER, INTENT (IN) :: kproma, kbdim, klev, klevp1, ktrac, klevm1
!---Included for in-cloud scavenging (Philip Stier, 19/01/06):----------
!!$INTEGER, INTENT (IN) :: krow
!---End Included for scavenging-----------------------------------------
REAL(dp),INTENT(IN) :: pdtime
REAL(dp),INTENT(IN) :: ptime_step_len

REAL(dp):: pten(kbdim,klev),        pqen(kbdim,klev),                  &
           pxen(kbdim,klev),        ptven(kbdim,klev),                 &
           puen(kbdim,klev),        pven(kbdim,klev),                  &
           ptte(kbdim,klev),        pqte(kbdim,klev),                  &
           pvom(kbdim,klev),        pvol(kbdim,klev),                  &
           pqsen(kbdim,klev),       pgeo(kbdim,klev),                  &
           paphp1(kbdim,klevp1),                                       &
           pverv(kbdim,klev)
REAL(dp),INTENT(OUT) :: pvom_cnv(kbdim,klev), pvol_cnv(kbdim,klev)
REAL(dp),INTENT(OUT) :: ptte_cnv(kbdim,klev), pqte_cnv(kbdim,klev)
REAL(dp),INTENT(OUT) ::pxtte_cnv(kbdim,klev,ktrac)
REAL(dp):: ptu(kbdim,klev),         pqu(kbdim,klev),                   &
           plu(kbdim,klev),         plude(kbdim,klev),                 &
           pmfu(kbdim,klev),        pmfd(kbdim,klev),                  &
           paprc(kbdim),            paprs(kbdim),                      &
           prsfc(kbdim),            pssfc(kbdim),                      &
           prain(kbdim),            pqhfla(kbdim)
REAL(dp):: pthvsig(kbdim)
INTEGER :: kcbot(kbdim),            kctop(kbdim),                      &
           ktype(kbdim)
REAL(dp):: pxtec(kbdim,klev),       pqtec(kbdim,klev),                 &
           pqude(kbdim,klev)
REAL(dp):: ztenh(kbdim,klev),       zqenh(kbdim,klev),                 &
           zxenh(kbdim,klev),                                          &
           zgeoh(kbdim,klev),       zqsenh(kbdim,klev),                &
           ztd(kbdim,klev),         zqd(kbdim,klev),                   &
           zmfus(kbdim,klev),       zmfds(kbdim,klev),                 &
           zmfuq(kbdim,klev),       zmfdq(kbdim,klev),                 &
           zdmfup(kbdim,klev),      zdmfdp(kbdim,klev),                &
           zmful(kbdim,klev),       zrfl(kbdim),                       &
           zuu(kbdim,klev),         zvu(kbdim,klev),                   &
           zud(kbdim,klev),         zvd(kbdim,klev)
REAL(dp):: zcpen(kbdim,klev),       zcpcu(kbdim,klev)
REAL(dp):: zentr(kbdim),            zhcbase(kbdim),                    &
           zmfub(kbdim),            zmfub1(kbdim),                     &
           zktype(kbdim),           zldcum(kbdim),                     &
           zcpcui(kbdim,klev),      zkcbot(kbdim),                     &
           zictop0(kbdim),          ztmp1(kbdim),                      &
           ztmp2(kbdim),            ztmp3(kbdim)
#ifdef __ibmspline__
REAL(dp):: za(kbdim)
#endif
REAL(dp):: zsfl(kbdim),             zdpmel(kbdim,klev)
REAL(dp):: zcape(kbdim),            zheat(kbdim)
REAL(dp):: ua(kbdim), dua(kbdim), ub(kbdim)
REAL(dp):: zhmin(kbdim),            zihmin(kbdim)
REAL(dp):: zhhatt(kbdim,klev),      zdqpbl(kbdim)
REAL(dp):: zdqcv(kbdim)
INTEGER :: ihmin(kbdim), ilo1(kbdim), ldidx(kbdim), loidx(kbdim)
INTEGER :: ilab(kbdim,klev),        idtop(kbdim),                      &
           ictop0(kbdim),           ilwmin(kbdim)
REAL(dp):: pxten(kbdim,klev,ktrac),                                    &
!!$           pxtte(kbdim,klev,ktrac),                                    &
           pxtu(kbdim,klev,ktrac),                                     &
           zxtenh(kbdim,klev,ktrac),zxtd(kbdim,klev,ktrac),            &
           zmfuxt(kbdim,klev,ktrac),zmfdxt(kbdim,klev,ktrac)
LOGICAL :: loddraf(kbdim),          ldland(kbdim)
LOGICAL :: ldcum(kbdim)
LOGICAL :: llo1
!---Included for scavenging in xtwetdep (Philip Stier, 19/02/04, UL, 28.3.07):-------
REAL(dp):: papp1(kbdim,klev)
REAL(dp):: zmwc(kbdim,klev),        zmrateprecip(kbdim,klev)
REAL(dp):: zmratesnow(kbdim,klev)
!---End Included for scavenging-----------------------------------------
!--- Included for prognostic CDNC/IC scheme ----------------------------
REAL(dp):: pcvcbot(kbdim),          pwcape(kbdim)
REAL(dp):: plul(kbdim,klev),        plui(kbdim,klev),                   &
           pludel(kbdim,klev),      pludei(kbdim,klev),                 &
           pxtecl(kbdim,klev),      pxteci(kbdim,klev),                 &
           pxtecnl(kbdim,klev),     pxtecni(kbdim,klev),                &
           zmfull(kbdim,klev),      zmfuli(kbdim,klev)
!!$REAL(dp) :: ptkem1(kbdim,klev)
!--- End Included for CDNC/IC ------------------------------------------
!
INTEGER :: nl, jl, jk, ikb, jt, itopm2, locnt, ldcnt
REAL(dp):: zcons2, ztau, zqumqe, zdqmin, zmfmax, zalvs, zalvdcp, zqalv &
         , zhsat, zes, zcor, zqsat, zqst1, zdqsdt, zgam, zzz, zhhat    &
         , zb, zbi,      zroi, zdz, zdhdz, zdepth, zfac, zrh, zeps     &
         , zjk, zhelp, zlo1, zpaphp1i, ztenhi, zkctop, zpbmpt, za1, za2
INTEGER :: icuasc
!
!  INTRINSIC FUNCTIONS
INTRINSIC MIN, MAX
!
!  Executable statements

!!$#ifdef __ICON__
!!$#else
!!$  IF (lconvmassfix) THEN
!!$     CALL xt_conv_massfix(kproma,         kbdim,         klev,              &
!!$                          klevp1,         ktrac,         krow,              &
!!$                          papp1,          paphp1,        pxtte,             &
!!$                          .TRUE.  )
!!$  END IF
!!$#endif

!-----------------------------------------------------------------------
!
!     1.           SPECIFY CONSTANTS AND PARAMETERS
!                  --------------------------------
!
!100 CONTINUE
!
  zcons2=1._dp/(g*ptime_step_len)
#ifdef __ICON__
  ztau=cmftau
#else
  ztau=MIN(3._dp*3600._dp,7200._dp*63._dp/nn)
#endif
  pwcape(:)=0._dp
!
!----------------------------------------------------------------------
!
!*    2.           INITIALIZE VALUES AT VERTICAL GRID POINTS IN 'CUINI'
!                  ---------------------------------------------------
!
!200 CONTINUE
  CALL cuini(ncvmicro, kproma, kbdim, klev, klevp1, klevm1,            &
             pten,     pqen,     pqsen,    pxen,     puen,     pven,   &
             ptven,    ktrac,                                          &
             pxten,    zxtenh,   pxtu,     zxtd,     zmfuxt,   zmfdxt, &
             pverv,    pgeo,     paphp1,   zgeoh,                      &
             ztenh,    zqenh,    zqsenh,   zxenh,    ilwmin,           &
             ptu,      pqu,      ztd,      zqd,                        &
             zuu,      zvu,      zud,      zvd,                        &
             pmfu,     pmfd,     zmfus,    zmfds,                      &
             zmfuq,    zmfdq,    zdmfup,   zdmfdp,                     &
             zcpen,    zcpcu,                                          &
             zdpmel,   plu,      plude,    pqude,    ilab,             &
!--- Included for prognostic CDNC/IC scheme ----------------------------
             plul,     plui,     pludel,   pludei                       )
!--- End Included for CDNC/IC ------------------------------------------
!
!-----------------------------------------------------------------------
!
!*    3.0          CLOUD BASE CALCULATIONS
!                  -----------------------
!
!300 CONTINUE
!
!*             (A) DETERMINE CLOUD BASE VALUES IN 'CUBASE'
!                  ---------------------------------------
!
  CALL cubase(lmfdudv,  cminbuoy, kproma, kbdim, klev, klevp1, klevm1, &
              ztenh,    zqenh,    zgeoh,    paphp1,    pthvsig,        &
              ptu,      pqu,      plu,                                 &
              puen,     pven,     zuu,      zvu,                       &
              zcpcu,                                                   &
              ldcum,    kcbot,    ilab)
!
!*             (B) DETERMINE TOTAL MOISTURE CONVERGENCE AND
!*                 THEN DECIDE ON TYPE OF CUMULUS CONVECTION
!                  -----------------------------------------
!
  zkcbot(1:kproma) = REAL(kcbot(1:kproma),dp)

  jk=1
  DO 310 jl=1,kproma
     zdqpbl(jl)=0.0_dp
     zdqcv(jl)=pqte(jl,jk)*(paphp1(jl,jk+1)-paphp1(jl,jk))
     idtop(jl)=0
310 END DO
  DO 320 jk=2,klev
     zjk = REAL(jk,dp)
     DO 315 jl=1,kproma
        zhelp      = paphp1(jl,jk+1)-paphp1(jl,jk)
        zdqcv(jl)  = zdqcv(jl)+pqte(jl,jk)*zhelp
        zdqpbl(jl) = zdqpbl(jl) + FSEL(zjk - zkcbot(jl),pqte(jl,jk),0._dp)*zhelp
#ifdef __ibmdbg__
        PRINT '(A6,I3,I4,E18.10,L3)','zdqpbl',jk,jl,zdqpbl(jl),(FSEL(zjk - zkcbot(jl),1._dp,0._dp)==1._dp)
#endif
315  END DO
320 END DO
!
!*             (C) DETERMINE MOISTURE SUPPLY FOR BOUNDARY LAYER
!*                 AND DETERMINE CLOUD BASE MASSFLUX IGNORING
!*                 THE EFFECTS OF DOWNDRAFTS AT THIS STAGE
!                  ------------------------------------------
!
  zldcum(1:kproma) = MERGE(1._dp,0._dp,ldcum(1:kproma))

!DIR$ IVDEP
!DIR$ CONCURRENT
!CDIR NODEP
!IBM* ASSERT(NODEPS)
!IBM* NOVECTOR
  DO 340 jl=1,kproma
     ikb=kcbot(jl)
     zqumqe=pqu(jl,ikb)+plu(jl,ikb)-zqenh(jl,ikb)
     zdqmin=MAX(0.01_dp*zqenh(jl,ikb),1.e-10_dp)
     zlo1 = FSEL(-zdqpbl(jl),0._dp,1._dp)
     zlo1 = FSEL(zdqmin - zqumqe,0._dp,zlo1) * zldcum(jl)
     zmfub(jl)=FSEL(-zlo1,0.01_dp,(zdqpbl(jl)/(g*MAX(zqumqe,zdqmin))))
     zmfmax=(paphp1(jl,ikb)-paphp1(jl,ikb-1))*zcons2
     zmfub(jl)=MIN(zmfub(jl),zmfmax)
     zldcum(jl) = zlo1
     zhelp = MAX(0._dp,-1.1_dp*pqhfla(jl)*g)
     zktype(jl) = FSEL(zhelp - zdqcv(jl),2._dp,1._dp)
     zentr(jl)  = FSEL(zhelp - zdqcv(jl),entrscv,entrpen)
     ktype(jl) = INT(zktype(jl))
#ifdef __ibmdbg__
     PRINT '(A6,I3,I4,2 E18.10,I3,L3)','zmfub',ikb,jl,zmfub(jl),zentr(jl),ktype(jl),(zlo1==1._dp)
#endif
340 END DO
  ldcum(1:kproma) = (zldcum(1:kproma).GT.0._dp)
!
!-----------------------------------------------------------------------
!*    4.0          DETERMINE CLOUD ASCENT FOR ENTRAINING PLUME
!                  -------------------------------------------
!
!400 CONTINUE
!
!*             (A) ESTIMATE CLOUD HEIGHT FOR ENTRAINMENT/DETRAINMENT
!*                 CALCULATIONS IN CUASC (MAX.POSSIBLE CLOUD HEIGHT
!*                 FOR NON-ENTRAINING PLUME, FOLLOWING A.-S.,1974)
!                  -------------------------------------------------
!
!DIR$ IVDEP
  DO 410 jl=1,kproma
     ikb=kcbot(jl)
     zalvs=FSEL(tmelt-ptu(jl,ikb),als,alv)
     zhcbase(jl) = zcpcu(jl,ikb)*ptu(jl,ikb) + zgeoh(jl,ikb) + zalvs*pqu(jl,ikb)
     zictop0(jl) = zkcbot(jl)-1._dp
410 END DO
  DO 430 jk=klev,1,-1
     zcpcui(1:kproma,jk) = 1._dp/zcpcu(1:kproma,jk)
     IF (jk <= klevm1 .AND. jk >= 3) THEN
        ! mpuetz: too few instructions (FP dependencies)
#ifdef __ibmspline__
        CALL prepare_ua_index_spline('cumastr',kproma,ztenh(1,jk),loidx(1),za(1))
        CALL lookup_ua_spline(kproma,loidx(1),za(1),ua(1),dua(1))
#else
        CALL prepare_ua_index('cumastr',kproma,ztenh(1,jk),loidx(1))
        CALL lookup_ua(kproma,loidx(1),ua(1),dua(1))
#endif
        CALL lookup_ubc(kproma,ztenh(1,jk),ub(1))
        zjk = REAL(jk,dp)
!IBM* NOVECTOR
        DO 421 jl=1,kproma
           ! mpuetz: move some of these into the previous loop
           zalvs=FSEL(tmelt - ztenh(jl,jk),als,alv)
           zalvdcp=zalvs*zcpcui(jl,jk)
           zqalv=1._dp/zalvs
           zhsat=zcpcu(jl,jk)*ztenh(jl,jk)+zgeoh(jl,jk)+zalvs*zqsenh(jl,jk)
           zpaphp1i = 1._dp/paphp1(jl,jk)
           zes=ua(jl)*zpaphp1i
           zes=MIN(0.5_dp,zes)
           zcor=1._dp/(1._dp-vtmpc1*zes)
           zqsat=zes*zcor
#ifdef __ibmspline2__
           zdqsdt=zpaphp1i*zcor**2*dua(jl)
#else
           zqst1=(ua(jl) + 0.001_dp*dua(jl))*zpaphp1i
           zqst1=MIN(0.5_dp,zqst1)
           zqst1=zqst1/(1._dp-vtmpc1*zqst1)
           zdqsdt=(zqst1-zqsat)*1000._dp
#endif
           zgam=FSEL(zes - 0.4_dp,zqsat*zcor*ub(jl),zalvdcp*zdqsdt)
           zzz=zcpcu(jl,jk)*ztenh(jl,jk)*vtmpc1
           zhhat=zhsat-((zzz+zgam*zzz)/(1._dp+zgam*zzz*zqalv))*             &
                MAX(zqsenh(jl,jk)-zqenh(jl,jk),0._dp)
           zhhatt(jl,jk)=zhhat
           zlo1 = FSEL(zjk - zictop0(jl),0._dp,1._dp)
           zlo1 = FSEL(zhhat - zhcbase(jl),0._dp,zlo1)
           zictop0(jl) = FSEL(-zlo1,zictop0(jl),zjk)
#ifdef __ibmdbg__
           PRINT '(A6,I3,I4,3 E18.10,I3)','ictop',jk,jl,zes,zqsdt,zhhatt(jl,jk),INT(zictop0(jl))
#endif
421     END DO
     END IF
430 END DO
  ictop0(1:kproma) = INT(zictop0(1:kproma))
!!
!!     DEEP CONVECTION IF CLOUD DEPTH > 200 HPA, ELSE SHALLOW
!!     (CLOUD DEPTH FROM NON-ENTRAINIG PLUME)
!!
!  DO jl=1,kproma
!     ktype(jl)=MERGE(1,2,                                              &
!                paphp1(jl,kcbot(jl))-paphp1(jl,ictop0(jl)).gt.2.E4_dp)
!     zentr(jl)=MERGE(entrpen,entrscv,ktype(jl).eq.1)
!  ENDDO
!!
!
!                  FIND LOWEST POSSIBLE ORG. DETRAINMENT LEVEL
!                  -------------------------------------------
!
  ldcnt = 0
  DO jl=1,kproma
     zhmin(jl)=0._dp
     zihmin(jl)=0._dp
     llo1=ldcum(jl).AND.ktype(jl).EQ.1
     IF(llo1) THEN
        zihmin(jl)=zkcbot(jl)
        ldcnt = ldcnt + 1
        ldidx(ldcnt) = jl
     ENDIF
  ENDDO
!
  zb=25._dp
  zbi=1._dp/(zb*g)
  DO jk=klev,1,-1

     ! mpuetz: compute the update criterion

     zjk = REAL(jk,dp)
!CDIR NODEP
!IBM* ASSERT(NODEPS)
     DO nl = 1,ldcnt
        jl = ldidx(nl)
        zlo1 = FSEL(zjk - zkcbot(jl),0._dp,1._dp)
        zlo1 = FSEL(zjk - zictop0(jl),zlo1,0._dp)
        zlo1 = FSEL(-ABS(zihmin(jl)-zkcbot(jl)),zlo1,0._dp)
        ilo1(nl) = INT(zlo1)
     END DO

     ! mpuetz: compute the indices of elements to be updated

     locnt = 0
     DO nl = 1,ldcnt
        IF (ilo1(nl).GT.0) THEN
           locnt = locnt + 1
           loidx(locnt) = ldidx(nl)
        END IF
     END DO

!CDIR NODEP
!IBM* ASSERT(NODEPS)
     DO nl=1,locnt
        jl = loidx(nl)

        zalvs = FSEL(tmelt - ztenh(jl,jk),als,alv)
!        zalvs=MERGE(alv,als,ztenh(jl,jk)>tmelt)
        ikb   = kcbot(jl)
        zroi  = SWDIV_NOCHK(rd*ztenh(jl,jk)*(1._dp+vtmpc1*zqenh(jl,jk)),paphp1(jl,jk))
        zdz   = (paphp1(jl,jk)-paphp1(jl,jk-1))*zroi/g

        za1 = (zcpen(jl,jk-1)*pten(jl,jk-1)        &
             - zcpen(jl,jk)*pten(jl,jk)            &
             + zalvs*(pqen(jl,jk-1) - pqen(jl,jk)) &
             + (pgeo(jl,jk-1)-pgeo(jl,jk)) )*g
        za2 = pgeo(jl,jk-1)-pgeo(jl,jk)
        zdhdz = SWDIV_NOCHK(za1, za2)

        zdepth    = zgeoh(jl,jk)-zgeoh(jl,ikb)

        ztmp1(nl) = zalvs
        ztmp2(nl) = zdz*zdhdz
        ztmp3(nl) = 1._dp+zdepth*zbi
     END DO

     ztmp3(1:locnt) = SQRT(ztmp3(1:locnt))

!CDIR NODEP
!IBM* ASSERT(NODEPS)
     DO nl=1,locnt
        jl = loidx(nl)
        zalvs     = ztmp1(nl)
        zfac      = ztmp3(nl)

        zdepth    = zgeoh(jl,jk)-zgeoh(jl,ikb)
        zhmin(jl) = zhmin(jl) + zfac*ztmp2(nl)
        zrh       =-zalvs*(zqsenh(jl,jk)-zqenh(jl,jk))*zfac

        zihmin(jl) = FSEL(zrh - zhmin(jl),zihmin(jl),zjk)
!        IF(zhmin(jl).GT.zrh) ihmin(jl)=jk
     ENDDO
  ENDDO
!
!CDIR NODEP
!IBM* ASSERT(NODEPS)
  DO nl=1,ldcnt
     jl = ldidx(nl)
#ifdef __ibmdbg__
     PRINT '(A6,I4,I4,I4)','ihmin',jl,INT(zihmin(jl)),ictop0(jl)
#endif
     zihmin(jl) = FSEL(zihmin(jl)-zictop0(jl),zihmin(jl),zictop0(jl))
!        IF(ihmin(jl).LT.ictop0(jl)) ihmin(jl)=ictop0(jl)
  ENDDO
  ihmin(1:kproma) = INT(zihmin(1:kproma))
!
!*             (B) DO ASCENT IN 'CUASC'IN ABSENCE OF DOWNDRAFTS
!                  --------------------------------------------
!
  icuasc=1
!
  CALL cuasc(ncvmicro, lmfdudv,  lmfmid, dlev, cmfctop, cprcon, cminbuoy, &
             pdtime, ptime_step_len,                                   &
             kproma, kbdim, klev, klevp1, klevm1,                      &
             ztenh,    zqenh,    puen,     pven,                       &
             ktrac,                                                    &
             zxtenh,   pxten,    pxtu,     zmfuxt,                     &
             pten,     pqen,     pqsen,                                &
             pgeo,     zgeoh,    paphp1,   pthvsig,                    &
             pqte,     pverv,    ilwmin,                               &
             ldcum,    ldland,   ktype,    ilab,                       &
             ptu,      pqu,      plu,      zuu,      zvu,              &
             pmfu,     zmfub,    zentr,                                &
             zmfus,    zmfuq,                                          &
             zmful,    plude,    pqude,    zdmfup,                     &
             ihmin,    zhhatt,   zhcbase,  zqsenh,                     &
             zcpen,    zcpcu,                                          &
             kcbot,    kctop,    ictop0,                               &
!---Included for scavenging in xtwetdep (Philip Stier, 19/02/04, UL, 28.3.07):-------
             zmwc,     zmrateprecip,  zmratesnow,                      &
!---End Included for scavenging-----------------------------------------
!--- Included for prognostic CDNC/IC scheme ----------------------------
             pludel,   pludei,   pxtecnl,  pxtecni,                    &
             plul,     plui,     papp1,                                &
             zmfull,   zmfuli                                          )
!!$             pwcape,   ptkem1,   krow,     icuasc                      )
!--- End Included for CDNC/IC ------------------------------------------
!
!*         (C) CHECK CLOUD DEPTH AND CHANGE ENTRAINMENT RATE ACCORDINGLY
!              CALCULATE PRECIPITATION RATE (FOR DOWNDRAFT CALCULATION)
!              ---------------------------------------------------------
!
!DIR$ IVDEP
  DO 440 jl=1,kproma
     zpbmpt=paphp1(jl,kcbot(jl))-paphp1(jl,kctop(jl))
     IF(ldcum(jl).AND.ktype(jl).EQ.1.AND.zpbmpt.LT.2.e4_dp) ktype(jl)=2
     IF(ldcum(jl)) ictop0(jl)=kctop(jl)
     IF(ktype(jl).EQ.2) zentr(jl)=entrscv
     zrfl(jl)=zdmfup(jl,1)
440 END DO
  DO 460 jk=2,klev
     DO 450 jl=1,kproma
        zrfl(jl)=zrfl(jl)+zdmfup(jl,jk)
450  END DO
460 END DO

  ! mpuetz: must recompute ldidx() since ktype could have changed

  ldcnt = 0
  DO jl=1,kproma
     llo1=ldcum(jl).AND.ktype(jl).EQ.1
     IF(llo1) THEN
        ldcnt = ldcnt + 1
        ldidx(ldcnt) = jl
     ENDIF
  ENDDO
!
!-----------------------------------------------------------------------
!*    5.0          CUMULUS DOWNDRAFT CALCULATIONS
!                  ------------------------------
!
!500 CONTINUE
!
  IF(lmfdd) THEN
!
!*             (A) DETERMINE LFS IN 'CUDLFS'
!                  -------------------------
!
     CALL cudlfs(ncvmicro, lmfdudv, lmfdd, kproma,kbdim, klev, klevp1, &
                 ztenh,    zqenh,    puen,     pven,                   &
                 ktrac,                                                &
                 zxtenh,   pxtu,     zxtd,     zmfdxt,                 &
                 zgeoh,    paphp1,                                     &
                 ptu,      pqu,      zuu,      zvu,                    &
                 ldcum,    kcbot,    kctop,    zmfub,    zrfl,         &
                 ztd,      zqd,      zud,      zvd,                    &
                 pmfd,     zmfds,    zmfdq,    zdmfdp,                 &
                 zcpcu,                                                &
                 idtop,    loddraf,plui)
!
!*            (B)  DETERMINE DOWNDRAFT T,Q AND FLUXES IN 'CUDDRAF'
!                  -----------------------------------------------
!
     CALL cuddraf(ncvmicro, lmfdudv,  kproma,   kbdim,  klev, klevp1,  &
                  ztenh,    zqenh,    puen,     pven,                  &
                  ktrac,                                               &
                  zxtenh,   zxtd,     zmfdxt,                          &
                  zgeoh,    paphp1,   zrfl,                            &
                  ztd,      zqd,      zud,      zvd,                   &
                  pmfd,     zmfds,    zmfdq,    zdmfdp,                &
                  zcpcu,                                               &
                  loddraf,plui)
!
  END IF
!
!*    5.5          RECALCULATE CLOUD BASE MASSFLUX FROM A
!*                 CAPE CLOSURE FOR DEEP CONVECTION (KTYPE=1)
!*                 AND BY PBL EQUILIBRUM TAKING DOWNDRAFTS INTO
!*                 ACCOUNT FOR SHALLOW CONVECTION (KTYPE=2)
!                  -------------------------------------------
!
  ! mpuetz: cuasc can modify kctop !!

  zkcbot(1:kproma) = REAL(kcbot(1:kproma),dp)

  DO jl=1,kproma
     zheat(jl)=0._dp
     zcape(jl)=0._dp
     zmfub1(jl)=zmfub(jl)
!     zhelp(jl) =0._dp
  ENDDO
!
  DO jk=1,klev
     zjk = REAL(jk,dp)
!CDIR NODEP
!IBM* ASSERT(NODEPS)
     DO nl = 1,ldcnt
        jl = ldidx(nl)
        zkctop = REAL(kctop(jl),dp)
        zlo1 = FSEL(zkcbot(jl) - zjk,1._dp,0._dp) &
             * FSEL(zkctop - zjk,0._dp,1._dp)
        ilo1(nl) = INT(zlo1)
     END DO

     locnt = 0
!CDIR NODEP
     DO nl = 1,ldcnt
        jl = ldidx(nl)
        IF(jk.LE.kcbot(jl).AND.jk.GT.kctop(jl)) THEN
!        IF (ilo1(nl).GT.0) THEN
           locnt = locnt + 1
           loidx(locnt) = jl
        END IF
     END DO

#ifdef __ibmdbg__
     PRINT '(A6,I3,2 I4)','cumas1',jk,ldcnt,locnt,INT(zkcbot(jl))
#endif

     ! mpuetz: there is reuse of zro, maybe we can calulate once and store

!CDIR NODEP
!IBM* ASSERT(NODEPS)
     DO nl=1,locnt
        jl = loidx(nl)

        zroi      = SWDIV_NOCHK(rd*ztenh(jl,jk)*(1._dp+vtmpc1*zqenh(jl,jk)),paphp1(jl,jk))
        ztenhi    = SWDIV_NOCHK(1._dp,ztenh(jl,jk))
        zdz       = (paphp1(jl,jk)-paphp1(jl,jk-1))*zroi/g
        zheat(jl) = zheat(jl) +                                           &
             (  (pten(jl,jk-1)-pten(jl,jk) + g*zdz*zcpcui(jl,jk))*ztenhi  &
             +  vtmpc1*(pqen(jl,jk-1)-pqen(jl,jk)) )                      &
             * (g*(pmfu(jl,jk)+pmfd(jl,jk)))*zroi
        zcape(jl) = zcape(jl) +                                           &
             ( g*(ptu(jl,jk)-ztenh(jl,jk))*ztenhi                         &
             + g*vtmpc1*(pqu(jl,jk)-zqenh(jl,jk))                         &
             - g*plu(jl,jk) ) * zdz

#ifdef __ibmdbg__
        PRINT '(A6,I3,I4,5 E18.10,I4)','zheat',jk,jl,zheat(jl),zcape(jl),zroi,ztenhi,zdz,kcbot(jl)
        PRINT '(A6,I3,I4,6 E18.10)','zcape',jk,jl,pten(jl,jk-1)-pten(jl,jk),pqen(jl,jk-1)-pqen(jl,jk),pmfu(jl,jk)+pmfd(jl,jk),ptu(jl,jk),pqu(jl,jk),plu(jl,jk)
#endif
    ENDDO
  ENDDO
!

!  DO jl=1,kproma
!     DO jk=2,klev
!        zro=paphp1(jl,jk)/(rd*ztenh(jl,jk)*                         &
!             (1._dp+vtmpc1*zqenh(jl,jk)))
!        zdz=(paphp1(jl,jk)-paphp1(jl,jk-1))/(g*zro)
!        zhelp(jl)=zhelp(jl) +                               &
!             (g*(ptu(jl,jk)-ztenh(jl,jk))/ztenh(jl,jk)     &
!             +g*vtmpc1*(pqu(jl,jk)-zqenh(jl,jk))      &
!             -g*plu(jl,jk) ) * zdz
!     ENDDO
!     if (zhelp(jl).lt.0._dp) zhelp(jl)=0._dp
!  ENDDO

  DO jl=1,kproma
     !--- Included for prognostic CDNC/IC scheme ----------------------------
!     IF(ldcum(jl)) THEN
!         pcvcbot(jl)=0._dp
!         pwcape(jl)=0._dp
!         IF (zcape(jl).GT.0._dp) THEN
!            pcvcbot(jl)=REAL(kcbot(jl),dp)
!            pwcape(jl)=0.5_dp*SQRT(zcape(jl))
!         ENDIF
!      ENDIF

!S.K. changed to
     pcvcbot(jl)=0._dp
     pwcape(jl)=0._dp
     IF(ldcum(jl)) THEN
        IF (zcape(jl).GT.0._dp) THEN
           pcvcbot(jl)=REAL(kcbot(jl),dp)
           pwcape(jl)=2._dp*SQRT(zcape(jl))
        ENDIF
     ENDIF
!S.K. end change
     IF(ldcum(jl).AND.ktype(jl).EQ.1) THEN
       ikb=kcbot(jl)
       zmfub1(jl) = SWDIV_NOCHK((zcape(jl)*zmfub(jl)),(zheat(jl)*ztau))
       zmfub1(jl) = MAX(zmfub1(jl),0.001_dp)
       zmfmax     = (paphp1(jl,ikb)-paphp1(jl,ikb-1))*zcons2
       zmfub1(jl) = MIN(zmfub1(jl),zmfmax)
     ENDIF
  ENDDO
!
!*                 RECALCULATE CONVECTIVE FLUXES DUE TO EFFECT OF
!*                 DOWNDRAFTS ON BOUNDARY LAYER MOISTURE BUDGET
!*                 FOR SHALLOW CONVECTION (KTYPE=2)
!                  --------------------------------------------
!
  ldcnt = 0
  DO jl=1,kproma
     llo1=ktype(jl).EQ.2
     IF(llo1) THEN
        ldcnt = ldcnt + 1
        ldidx(ldcnt) = jl
     ENDIF
  ENDDO
!DIR$ IVDEP
!CDIR NODEP
!IBM* ASSERT(NODEPS)
  DO 520 nl=1,ldcnt
     jl = ldidx(nl)
     ikb=kcbot(jl)
     llo1=pmfd(jl,ikb).LT.0._dp.AND.loddraf(jl)
     zeps=MERGE(cmfdeps,0._dp,llo1)
     zqumqe=pqu(jl,ikb)+plu(jl,ikb)-                                &
          zeps*zqd(jl,ikb)-(1._dp-zeps)*zqenh(jl,ikb)
     zdqmin=MAX(0.01_dp*zqenh(jl,ikb),1.e-10_dp)
     zmfmax=(paphp1(jl,ikb)-paphp1(jl,ikb-1))*zcons2
     llo1=zdqpbl(jl).GT.0._dp.AND.zqumqe.GT.zdqmin.AND.ldcum(jl)    &
          .AND.zmfub(jl).LT.zmfmax
     zmfub1(jl)=MERGE(zdqpbl(jl)/(g*MAX(zqumqe,zdqmin)),            &
          zmfub(jl),llo1)
     zmfub1(jl)=MERGE(zmfub1(jl),zmfub(jl),                         &
          ABS(zmfub1(jl)-zmfub(jl)).LT.0.2_dp*zmfub(jl))
#ifdef __ibmdbg__
     PRINT '(A6,I4,E18.10)','zmfub1',jl,zmfub1(jl)
#endif
520 END DO
  ldcnt = 0
  DO jl=1,kproma
     IF (ldcum(jl)) THEN
        ldcnt = ldcnt + 1
        ldidx(ldcnt) = jl
     ENDIF
  ENDDO
  DO 540 jk=1,klev
!CDIR NODEP
!IBM* ASSERT(NODEPS)
     DO 530 nl=1,ldcnt
        jl = ldidx(nl)

        ztmp1(nl) = SWDIV_NOCHK(zmfub1(jl),MAX(zmfub(jl),1.e-10_dp))
        zfac      = ztmp1(nl)
        pmfd(jl,jk)   = pmfd(jl,jk)*zfac
        zmfds(jl,jk)  = zmfds(jl,jk)*zfac
        zmfdq(jl,jk)  = zmfdq(jl,jk)*zfac
        zdmfdp(jl,jk) = zdmfdp(jl,jk)*zfac
#ifdef __ibmdbg__
        PRINT '(A6,I3,I4,E18.10)','zfac',jk,jl,ztmp1(nl)
#endif
530  END DO
!
!IBM* unroll(4)
     DO 5304 jt=1,ktrac
!CDIR NODEP
!IBM* ASSERT(NODEPS)
        DO 5302 nl=1,ldcnt
           jl = ldidx(nl)
           zfac             = ztmp1(nl)
           zmfdxt(jl,jk,jt) = zmfdxt(jl,jk,jt)*zfac
5302    END DO
5304 END DO
!
540 END DO
!
!*                 NEW VALUES OF CLOUD BASE MASS FLUX
!                  ----------------------------------
!
  DO 550 nl=1,ldcnt
     jl = ldidx(nl)
     zmfub(jl) = zmfub1(jl)
550 END DO
!
!-----------------------------------------------------------------------
!*    6.0          DETERMINE FINAL CLOUD ASCENT FOR ENTRAINING PLUME
!*                 FOR PENETRATIVE CONVECTION (TYPE=1),
!*                 FOR SHALLOW TO MEDIUM CONVECTION (TYPE=2)
!*                 AND FOR MID-LEVEL CONVECTION (TYPE=3).
!                  --------------------------------------------------
  icuasc=2
!
!600 CONTINUE
  CALL cuasc(ncvmicro, lmfdudv, lmfmid, dlev, cmfctop, cprcon, cminbuoy,&
             pdtime, ptime_step_len,        &
             kproma, kbdim, klev, klevp1, klevm1,                      &
             ztenh,    zqenh,    puen,     pven,                       &
             ktrac,                                                    &
             zxtenh,   pxten,    pxtu,     zmfuxt,                     &
             pten,     pqen,     pqsen,                                &
             pgeo,     zgeoh,    paphp1,   pthvsig,                    &
             pqte,     pverv,    ilwmin,                               &
             ldcum,    ldland,   ktype,    ilab,                       &
             ptu,      pqu,      plu,      zuu,      zvu,              &
             pmfu,     zmfub,    zentr,                                &
             zmfus,    zmfuq,                                          &
             zmful,    plude,    pqude,    zdmfup,                     &
             ihmin,    zhhatt,   zhcbase,  zqsenh,                     &
             zcpen,    zcpcu,                                          &
             kcbot,    kctop,    ictop0,                               &
!---Included for scavenging in xtwetdep (Philip Stier, 19/02/04, UL, 28.3.07):-------
             zmwc,     zmrateprecip,  zmratesnow,                      &
!---End Included for scavenging-----------------------------------------
!--- Included for prognostic CDNC/IC scheme ----------------------------
             pludel,   pludei,   pxtecnl,  pxtecni,                    &
             plul,     plui,     papp1,                                &
             zmfull,   zmfuli                                          )
!!$             pwcape,   ptkem1,   krow,     icuasc                      )
!--- End Included for CDNC/IC ------------------------------------------
!
!-----------------------------------------------------------------------
!*    7.0          DETERMINE FINAL CONVECTIVE FLUXES IN 'CUFLX'
!                  ------------------------------------------
!
!700 CONTINUE

  CALL cuflx(ncvmicro,  &
             ptime_step_len, kproma,   kbdim,    klev,     klevp1,     &
             pqen,     pqsen,    ztenh,    zqenh,                      &
             ktrac,                                                    &
!---Included for scavenging in xtwetdep (Philip Stier, 28/03/01, UL, 28.3.07):-------
!!$             krow,                                                     &
!!$             pxtte,    pxtu,     ptu,                                  &
!!$             zmwc,     zmrateprecip,   zmratesnow,                     &
!---End Included for scavenging-----------------------------------------
             zxtenh,   zmfuxt,   zmfdxt,                               &
             paphp1,   zgeoh,                                          &
             kcbot,    kctop,    idtop,                                &
             ktype,    loddraf,  ldcum,                                &
             pmfu,     pmfd,     zmfus,    zmfds,                      &
             zmfuq,    zmfdq,    zmful,                                &
             zdmfup,   zdmfdp,   zrfl,     prain,                      &
             zcpcu,                                                    &
             pten,     zsfl,     zdpmel,   itopm2,                     &
!--- Included for prognostic CDNC/IC scheme ----------------------------
             zmfull,   zmfuli)
!!$             plul,     plui)
!--- End Included for CDNC/IC ------------------------------------------
!-----------------------------------------------------------------------
!
!*    8.0          UPDATE TENDENCIES FOR T AND Q IN SUBROUTINE CUDTDQ
!                  --------------------------------------------------
!
!800 CONTINUE
  CALL cudtdq(ncvmicro, &
              pdtime, kproma, kbdim, klev, klevp1, itopm2, ldcum, ktrac, &
!--- Included for dust emissions (Philip Stier 23/01/06)-----------------
!!$              krow,                                                    &
!--- End Included for dust emissions in ---------------------------------
              paphp1,   pten,     ptte,     pqte,                      &
!!$              pxtte,                                                   &
              pxtec,                                                   &
!!$              zmfuxt,   zmfdxt,                                        &
              zmfus,    zmfds,    zmfuq,    zmfdq,                     &
              zmful,    zdmfup,   zdmfdp,   plude,                     &
              zdpmel,   zrfl,     zsfl,                                &
              zcpen,    pqtec,    pqude,                               &
              prsfc,    pssfc,    paprc,    paprs,                     &
!--- Included for prognostic CDNC/IC scheme ----------------------------
              pludel,   pludei,   pxtecl,   pxteci,                    &
              zmfull,   zmfuli,                                        &
!!$              plui,                                                    &
!--- End Included for CDNC/IC ------------------------------------------
              ptte_cnv, pqte_cnv, pxtte_cnv                        )
!
!-----------------------------------------------------------------------
!
!*    9.0          UPDATE TENDENCIES FOR U AND U IN SUBROUTINE CUDUDV
!                  --------------------------------------------------
!
!900 CONTINUE
  IF(lmfdudv) THEN
     CALL cududv(kproma,   kbdim,    klev,     klevp1,                 &
                 itopm2,   ktype,    kcbot,    paphp1,   ldcum,        &
                 puen,     pven,     pvom,     pvol,                   &
                 zuu,      zud,      zvu,      zvd,                    &
                 pmfu,     pmfd,     pvom_cnv, pvol_cnv )
!
  END IF
!
!!$#ifdef __ICON__
!!$#else
!!$  IF (lconvmassfix) THEN
!!$     CALL xt_conv_massfix(kproma,         kbdim,         klev,              &
!!$                          klevp1,         ktrac,         krow,              &
!!$                          papp1,          paphp1,        pxtte,             &
!!$                          .FALSE.  )
!!$  END IF
!!$#endif
!
!1000 CONTINUE
!
    RETURN
  END SUBROUTINE cumastr
  !-------------

END MODULE mo_cumastr
