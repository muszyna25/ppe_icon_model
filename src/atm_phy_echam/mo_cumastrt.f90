!>
!! Module contains subroutine cumastrt -- master routine for
!! the Tiedtke cumulus convection scheme.
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
MODULE mo_cumastrt

  USE mo_kind,                ONLY: dp

#ifdef __ICON__
  USE mo_physical_constants,  ONLY: g=>grav, alv, als, tmelt, vtmpc1
  USE mo_echam_conv_nml,      ONLY: entrpen, entrscv,      &
                                  & lmfdd, cmfdeps, lmfdudv
#else
  USE mo_constants,           ONLY: g, alv, als, tmelt, vtmpc1
  USE mo_cumulus_flux,        ONLY: entrpen, entrscv, lmfdd, cmfdeps, lmfdudv
  USE mo_tracer_processes,    ONLY: xt_conv_massfix
  USE mo_param_switches,      ONLY: lconvmassfix
#endif

  USE mo_cuini,               ONLY: cuini
  USE mo_cuasct,              ONLY: cuasct
  USE mo_cubase,              ONLY: cubase
  USE mo_cuflx,               ONLY: cuflx
  USE mo_cudtdq,              ONLY: cudtdq
  USE mo_cududv,              ONLY: cududv
  USE mo_cudlfs,              ONLY: cudlfs
  USE mo_cuddraf,             ONLY: cuddraf
  USE mo_convect_tables,      ONLY: tlucua,                &! table a
                                  & tlucub,                &! table b
                                  & jptlucu1, jptlucu2,    &
                                  & lookuperror

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: cumastrt

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

CONTAINS
  !>
  !!
  SUBROUTINE cumastrt( ncvmicro, pdtime, ptime_step_len,                    &
                       kproma, kbdim, klev, klevp1, klevm1,                 &
                       ilab,                                                &
!0                     krow,                                                &
!0                     papp1,                                               &!in-cloud scavenging
                       pten,     pqen,     pxen,     puen,     pven,        &
                       ptven,    ktrac,    ldland,                          &
                       pxten,    pxtu,                                      &
!0                     pxtte,                                               &
                       pverv,    pqsen,    pqhfla,                          &
                       paphp1,   pgeo,                                      &
                       ptte,     pqte,     pvom,     pvol,                  &
                       prsfc,    pssfc,    paprc,    paprs,    pxtec,       &
                       pqtec,    pqude,                                     &
                       ldcum,    ktype,    kcbot,    kctop,                 &
                       ptu,      pqu,      plu,      plude,                 &
                       pmfu,     pmfd,     prain,    pthvsig,               &
                       pxtecl,   pxteci,                                    &!prognostic CDNC/IC
                       ptte_cnv, pvom_cnv, pvol_cnv, pqte_cnv, pxtte_cnv  )
!
!**** *CUMASTRT*  MASTER ROUTINE FOR CUMULUS MASSFLUX-SCHEME
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
!          *CUMASTRT* IS CALLED FROM *CUCALL* IN CASE ICONV .EQ. 2
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
!       CUASCT: CLOUD ASCENT FOR ENTRAINING PLUME
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
!     MODEL PARAMETERS (DEFINED IN MODULE mo_cumulus_flux)
!     ----------------------------------------------------
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
INTEGER, INTENT (IN) :: kproma, kbdim, klev, klevp1, ktrac, klevm1
!---Included for in-cloud scavenging (Philip Stier, 19/01/06):----------
!!$INTEGER, INTENT (IN) :: krow
!---End Included for scavenging-----------------------------------------
REAL(dp),INTENT(IN) :: pdtime
REAL(dp),INTENT(IN) :: ptime_step_len
!
REAL(dp):: pten(kbdim,klev),        pqen(kbdim,klev),                  &
           pxen(kbdim,klev),        ptven(kbdim,klev),                 &
           puen(kbdim,klev),        pven(kbdim,klev),                  &
           ptte(kbdim,klev),        pqte(kbdim,klev),                  &
           pvom(kbdim,klev),        pvol(kbdim,klev),                  &
           pqsen(kbdim,klev),       pgeo(kbdim,klev),                  &
           paphp1(kbdim,klevp1),                                       &
           pverv(kbdim,klev)
REAL(dp),INTENT(OUT) :: pvom_cnv(kbdim,klev),pvol_cnv(kbdim,klev)
REAL(dp),INTENT(OUT) :: ptte_cnv(kbdim,klev),pqte_cnv(kbdim,klev)
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
           zdqpbl(kbdim),           zdqcv(kbdim)
REAL(dp):: zsfl(kbdim),             zdpmel(kbdim,klev)
INTEGER :: ilab(kbdim,klev),        idtop(kbdim),                      &
           ictop0(kbdim),           ilwmin(kbdim)
REAL(dp):: pxten(kbdim,klev,ktrac),                                    &
!!$           pxtte(kbdim,klev,ktrac),                                    &
           pxtu(kbdim,klev,ktrac),                                     &
           zxtenh(kbdim,klev,ktrac),zxtd(kbdim,klev,ktrac),            &
           zmfuxt(kbdim,klev,ktrac),zmfdxt(kbdim,klev,ktrac)
LOGICAL :: loddraf(kbdim),          ldland(kbdim)
LOGICAL :: ldcum(kbdim)
LOGICAL :: llo1, lo
!---Included for scavenging in xtwetdep (Philip Stier, 19/02/04):-------
!!$REAL(dp):: papp1(kbdim,klev)
REAL(dp):: zmwc(kbdim,klev),        zmrateprecip(kbdim,klev)
!!$REAL(dp):: zmratesnow(kbdim,klev)
!---End Included for scavenging-----------------------------------------
!--- Included for prognostic CDNC/IC scheme ----------------------------
REAL(dp):: plul(kbdim,klev),        plui(kbdim,klev),                   &
           pludel(kbdim,klev),      pludei(kbdim,klev),                 &
           pxtecl(kbdim,klev),      pxteci(kbdim,klev),                 &
           zmfull(kbdim,klev),      zmfuli(kbdim,klev)
!--- End Included for CDNC/IC ------------------------------------------
!
INTEGER :: jl, jk, ikb, it, it1, jt, itopm2
REAL(dp):: zcons2, zqumqe, zdqmin, zmfmax, zalvs, zalvdcp, zqalv       &
         , zhsat, zes, zcor, zqsat, zqst1, zdqsdt, zgam, zzz, zhhat    &
         , zfac, zpbmpt, zeps
LOGICAL :: lookupoverflow
!
!  INTRINSIC FUNCTIONS
INTRINSIC MIN, MAX
!
!  Executable statements

  lookupoverflow = .FALSE.

!!$#ifdef __ICON__
!!$#else
!!$  IF(lconvmassfix) THEN
!!$     CALL xt_conv_massfix(kproma,         kbdim,            klev,    &
!!$                          klevp1,         ktrac,            krow,    &
!!$                          papp1,          paphp1,           pxtte,   &
!!$                          .TRUE.                                     )
!!$  END IF
!!$#endif
!---------------------------------------------------------------------
!
!     1.           SPECIFY CONSTANTS AND PARAMETERS
!                  --------------------------------
!
!100 CONTINUE
!
  zcons2=1._dp/(g*ptime_step_len)
!
!---------------------------------------------------------------------
!
!*    2.           INITIALIZE VALUES AT VERTICAL GRID POINTS IN 'CUINI'
!                  ---------------------------------------------------
!
!200 CONTINUE
  CALL cuini(ncvmicro, kproma,   kbdim,    klev,     klevp1,   klevm1, &
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
!---------------------------------------------------------------------
!
!*    3.0          CLOUD BASE CALCULATIONS
!                  -----------------------
!
!300 CONTINUE
!
!*             (A) DETERMINE CLOUD BASE VALUES IN 'CUBASE'
!                  ---------------------------------------
!
  CALL cubase(kproma, kbdim, klev, klevp1, klevm1,                     &
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
  jk=1
  DO 310 jl=1,kproma
     zdqpbl(jl)=0.0_dp
     zdqcv(jl)=pqte(jl,jk)*(paphp1(jl,jk+1)-paphp1(jl,jk))
     idtop(jl)=0
310 END DO
  DO 320 jk=2,klev
     DO 315 jl=1,kproma
        zdqcv(jl)=zdqcv(jl)+pqte(jl,jk)*(paphp1(jl,jk+1)-paphp1(jl,jk))
        IF(jk.GE.kcbot(jl)) zdqpbl(jl)=zdqpbl(jl)+pqte(jl,jk)          &
                             *(paphp1(jl,jk+1)-paphp1(jl,jk))
315  END DO
320 END DO
!
!*             (C) DETERMINE MOISTURE SUPPLY FOR BOUNDARY LAYER
!*                 AND DETERMINE CLOUD BASE MASSFLUX IGNORING
!*                 THE EFFECTS OF DOWNDRAFTS AT THIS STAGE
!                  ------------------------------------------
!
!DIR$ IVDEP
  DO 340 jl=1,kproma
     ikb=kcbot(jl)
     zqumqe=pqu(jl,ikb)+plu(jl,ikb)-zqenh(jl,ikb)
     zdqmin=MAX(0.01_dp*zqenh(jl,ikb),1.e-10_dp)
     llo1=zdqpbl(jl).GT.0._dp.AND.zqumqe.GT.zdqmin.AND.ldcum(jl)
     zmfub(jl)=MERGE(zdqpbl(jl)/(g*MAX(zqumqe,zdqmin)),0.01_dp,llo1)
     zmfmax=(paphp1(jl,ikb)-paphp1(jl,ikb-1))*zcons2
     zmfub(jl)=MIN(zmfub(jl),zmfmax)
     IF(.NOT.llo1) ldcum(jl)=.FALSE.
     ktype(jl)=MERGE(1,2,zdqcv(jl).GT.MAX(0._dp,-1.1_dp*pqhfla(jl)*g))
     zentr(jl)=MERGE(entrpen,entrscv,ktype(jl).EQ.1)
340 END DO
!
!---------------------------------------------------------------------
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
     zalvs=MERGE(alv,als,ptu(jl,ikb)>tmelt)
     zhcbase(jl)=zcpcu(jl,ikb)*ptu(jl,ikb)+zgeoh(jl,ikb)+zalvs*pqu(jl,ikb)
     ictop0(jl)=kcbot(jl)-1
410 END DO
  DO 430 jk=klevm1,3,-1
!DIR$ IVDEP
     DO 420 jl=1,kproma
        zalvs=MERGE(alv,als,ztenh(jl,jk)>tmelt)
        zalvdcp=zalvs/zcpcu(jl,jk)
        zqalv=1._dp/zalvs
        zhsat=zcpcu(jl,jk)*ztenh(jl,jk)+zgeoh(jl,jk)+zalvs*zqsenh(jl,jk)
        it = NINT(ztenh(jl,jk)*1000._dp)
        IF (it<jptlucu1 .OR. it>jptlucu2) lookupoverflow = .TRUE.
        it = MAX(MIN(it,jptlucu2),jptlucu1)
        zes=tlucua(it)/paphp1(jl,jk)
        zes=MIN(0.5_dp,zes)
        LO=zes<0.4_dp
        zcor=1._dp/(1._dp-vtmpc1*zes)
        zqsat=zes*zcor
        it1=it+1
        it1 = MAX(MIN(it1,jptlucu2),jptlucu1)
        zqst1=tlucua(it1)/paphp1(jl,jk)
        zqst1=MIN(0.5_dp,zqst1)
        zqst1=zqst1/(1._dp-vtmpc1*zqst1)
        zdqsdt=(zqst1-zqsat)*1000._dp
        zgam=MERGE(zalvdcp*zdqsdt,zqsat*zcor*tlucub(it),LO)
        zzz=zcpcu(jl,jk)*ztenh(jl,jk)*vtmpc1
        zhhat=zhsat-(zzz+zgam*zzz)/(1._dp+zgam*zzz*zqalv)*             &
                       MAX(zqsenh(jl,jk)-zqenh(jl,jk),0._dp)
        IF(jk.LT.ictop0(jl).AND.zhcbase(jl).GT.zhhat) ictop0(jl)=jk
420  END DO
430 END DO
!
     IF (lookupoverflow) CALL lookuperror ('cumastrt')
!!
!!     DEEP CONVECTION IF CLOUD DEPTH > 200 HPA, ELSE SHALLOW
!!     (CLOUD DEPTH FROM NON-ENTRAINIG PLUME)
!!
!  DO jl=1,kproma
!     ktype(jl)=MERGE(1,2,                                             &
!                paphp1(jl,kcbot(jl))-paphp1(jl,ictop0(jl)).gt.2.E4_dp)
!     zentr(jl)=MERGE(entrpen,entrscv,ktype(jl).eq.1)
!  ENDDO
!!
!*             (B) DO ASCENT IN 'CUASCT' IN ABSENCE OF DOWNDRAFTS
!                  ----------------------------------------------
!
  CALL cuasct(ptime_step_len, kproma, kbdim, klev, klevp1, klevm1,     &
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
             zcpen,    zcpcu,                                          &
             kcbot,    kctop,    ictop0,                               &
!---Included for scavenging in xtwetdep (Philip Stier, 19/02/04):-------
             zmwc,     zmrateprecip,                                   &
!---End Included for scavenging-----------------------------------------
!--- Included for prognostic CDNC/IC scheme ----------------------------
             plul,     plui,     zmfull,   zmfuli)
!--- End Included for CDNC/IC ------------------------------------------
!
!*        (C) CHECK CLOUD DEPTH AND CHANGE ENTRAINMENT RATE ACCORDINGLY
!             CALCULATE PRECIPITATION RATE (FOR DOWNDRAFT CALCULATION)
!             ---------------------------------------------------------
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
!
!---------------------------------------------------------------------
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
     CALL cudlfs(ncvmicro, kproma,   kbdim,    klev,     klevp1,       &
                 ztenh,    zqenh,    puen,     pven,                   &
                 ktrac,                                                &
                 zxtenh,   pxtu,     zxtd,     zmfdxt,                 &
                 zgeoh,    paphp1,                                     &
                 ptu,      pqu,      zuu,      zvu,                    &
                 ldcum,    kcbot,    kctop,    zmfub,    zrfl,         &
                 ztd,      zqd,      zud,      zvd,                    &
                 pmfd,     zmfds,    zmfdq,    zdmfdp,                 &
                 zcpcu,                                                &
                 idtop,    loddraf,  plui)
!
!*            (B)  DETERMINE DOWNDRAFT T,Q AND FLUXES IN 'CUDDRAF'
!                 -----------------------------------------------
!
     CALL cuddraf(ncvmicro, kproma,   kbdim,    klev,     klevp1,      &
                  ztenh,    zqenh,    puen,     pven,                  &
                  ktrac,                                               &
                  zxtenh,   zxtd,     zmfdxt,                          &
                  zgeoh,    paphp1,   zrfl,                            &
                  ztd,      zqd,      zud,      zvd,                   &
                  pmfd,     zmfds,    zmfdq,    zdmfdp,                &
                  zcpcu,                                               &
                  loddraf,  plui)
!
!*            (C)  RECALCULATE CONVECTIVE FLUXES DUE TO EFFECT OF
!*                 DOWNDRAFTS ON BOUNDARY LAYER MOISTURE BUDGET
!                  --------------------------------------------
!
!DIR$ IVDEP
  DO 520 jl=1,kproma
     IF(loddraf(jl)) THEN
        ikb=kcbot(jl)
        llo1=pmfd(jl,ikb).LT.0._dp
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
                           (ktype(jl).EQ.1.OR.ktype(jl).EQ.2) .AND.    &
                         ABS(zmfub1(jl)-zmfub(jl)).LT.0.2_dp*zmfub(jl))
     END IF
520 END DO
  DO 540 jk=1,klev
     DO 530 jl=1,kproma
        IF(loddraf(jl)) THEN
           zfac=zmfub1(jl)/MAX(zmfub(jl),1.e-10_dp)
           pmfd(jl,jk)=pmfd(jl,jk)*zfac
           zmfds(jl,jk)=zmfds(jl,jk)*zfac
           zmfdq(jl,jk)=zmfdq(jl,jk)*zfac
           zdmfdp(jl,jk)=zdmfdp(jl,jk)*zfac
        END IF
530  END DO
!
     DO 5304 jt=1,ktrac
        DO 5302 jl=1,kproma
           IF(loddraf(jl)) THEN
              zfac=zmfub1(jl)/MAX(zmfub(jl),1.e-10_dp)
              zmfdxt(jl,jk,jt)=zmfdxt(jl,jk,jt)*zfac
           ENDIF
5302    END DO
5304 END DO
!
540 END DO
!
!*                 NEW VALUES OF CLOUD BASE MASS FLUX
!                  ----------------------------------
  DO 550 jl=1,kproma
     IF(loddraf(jl)) zmfub(jl)=zmfub1(jl)
550 END DO
!
  END IF
!
!---------------------------------------------------------------------
!*    6.0          DETERMINE FINAL CLOUD ASCENT FOR ENTRAINING PLUME
!*                 FOR PENETRATIVE CONVECTION (TYPE=1),
!*                 FOR SHALLOW TO MEDIUM CONVECTION (TYPE=2)
!*                 AND FOR MID-LEVEL CONVECTION (TYPE=3).
!                  ----------------------------------------------------
!600 CONTINUE
!
  CALL cuasct(ptime_step_len, kproma,  kbdim,    klev,     klevp1,   klevm1,  &
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
             zcpen,    zcpcu,                                          &
             kcbot,    kctop,    ictop0,                               &
!---Included for scavenging in xtwetdep (Philip Stier, 19/02/04):-------
             zmwc,     zmrateprecip,                                   &
!---End Included for scavenging-----------------------------------------
!--- Included for prognostic CDNC/IC scheme ----------------------------
             plul,     plui,     zmfull,   zmfuli)
!--- End Included for CDNC/IC ------------------------------------------
!
!---------------------------------------------------------------------
!*    7.0          DETERMINE FINAL CONVECTIVE FLUXES IN 'CUFLX'
!                  ------------------------------------------
!
!700 CONTINUE
  CALL cuflx(ncvmicro, ptime_step_len, kproma, kbdim, klev, klevp1,    &
             pqen,     pqsen,    ztenh,    zqenh,                      &
             ktrac,                                                    &
!---Included for scavenging in xtwetdep (Philip Stier, 28/03/01):-------
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
!---------------------------------------------------------------------
!
!*    8.0          UPDATE TENDENCIES FOR T AND Q IN SUBROUTINE CUDTDQ
!                  --------------------------------------------------
!
!800 CONTINUE
  CALL cudtdq(ncvmicro, pdtime,                                        &
              kproma, kbdim, klev, klevp1, itopm2, ldcum, ktrac,       &
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
              ptte_cnv, pqte_cnv, pxtte_cnv                     )
!
!---------------------------------------------------------------------
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
                 pmfu,     pmfd,     pvom_cnv, pvol_cnv         )
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
!
  RETURN
  END SUBROUTINE cumastrt
  !-------------

END MODULE mo_cumastrt
