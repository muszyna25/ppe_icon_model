!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
MODULE mo_cumastr
#ifdef __xlC__
@PROCESS HOT
#endif
#ifndef __xlC__
#define FSEL(a,b,c) MERGE(b,c,(a) >= 0._wp)
#define SWDIV_NOCHK(a,b) ((a)/(b))
#endif

  USE mo_kind,               ONLY: wp
!  USE mo_control,            ONLY: nn
  USE mo_physical_constants, ONLY: grav, alv, als, tmelt, vtmpc1, rd, cpd, cpv
!  USE mo_time_control,       ONLY: time_step_len
!  USE mo_param_switches,     ONLY: iconv, lconvmassfix
  USE mo_echam_convect_tables,     ONLY: prepare_ua_index_spline,lookup_ua_spline, &
                                   lookup_ubc
!  USE mo_tracer_processes,   ONLY: xt_conv_massfix
  USE mo_echam_conv_constants, ONLY: entrpen, entrscv, lmfdd, cmfdeps, lmfdudv, cmftau
  USE mo_cuinitialize,       ONLY: cuini, cubase
  USE mo_cuascent,           ONLY: cuasc, cuasct
  USE mo_cudescent,          ONLY: cudlfs, cuddraf
  USE mo_cufluxdts,          ONLY: cuflx, cudtdq, cududv
  USE mo_echam_phy_memory,   ONLY: t_echam_phy_field, prm_field

  IMPLICIT NONE

PRIVATE

PUBLIC :: cucall

CONTAINS 

SUBROUTINE cucall(   iconv,                                          &! in
                     nmctop, cevapcu,                                &! in
                     kproma, kbdim, klev, klevp1, klevm1,            &
                     ktrac,                                          &
                     krow,                                           &
                     time_step_len,                                  &! in
                     ldland,                                         &! in
                     ptm1,     pum1,     pvm1,                       &! in
                     pverv,                                          &! in
                     pgeo,                                           &! in
                     pqm1,     pxlm1,    pxim1, pxtm1,               &! in
                     pcd, pcv,                                       &! in
                     pxlte,    pxite,                                &! in
                     papp1,    paphp1,                               &! in
                     pqhfla,                                         &! in
                     pthvsig,                                        &! in
                     ptte,     pvom,     pvol,                       &! inout
                     pqte,     pxtte,                                &! inout
                     pqtec,                                          &! inout
                     pch_concloud, pcon_dtrl, pcon_dtri,             &! inout
                     pcon_iqte,                                      &! inout
                     pxtecl,   pxteci,                               &
                     prsfc,    pssfc,                                &! out
                     ptopmax,                                        &! inout
                     ktype,                                          &! inout
                     ilab,                                           &! out
!                     pcvcbot,  pwcape,                               &! out
                     ptte_cnv, pvom_cnv, pvol_cnv, pqte_cnv,         &! out
                     pxtte_cnv                               )        ! out

!
!
!          *CUCALL* - MASTER ROUTINE - PROVIDES INTERFACE FOR:
!                     *CUMASTR* (CUMULUS PARAMETERIZATION)
!
!           M.TIEDTKE      E.C.M.W.F.     12/1989
!
!**   PURPOSE.
!     --------
!
!          *CUCALL* - INTERFACE FOR *CUMASTR*:
!                     PROVIDES INPUT FOR CUMASTR
!                     RECEIVES UPDATED TENDENCIES, PRECIPITATION.
!
!**   INTERFACE.
!     ----------
!
!          *CUCALL* IS CALLED FROM *PHYSC*
!
!     EXTERNALS.
!     ----------
!
!          CUMASTR, CUMASTRT OR CUMASTRH
!

  INTEGER, INTENT (IN) :: iconv, klev, klevm1, klevp1, kproma, kbdim, ktrac
  REAL(wp),INTENT (IN) :: cevapcu(:)
  INTEGER, INTENT (IN) :: nmctop
  INTEGER, INTENT (IN) :: krow
  REAL(wp),INTENT (IN) :: time_step_len, pcd, pcv

  REAL(wp)::  ptm1(kbdim,klev),         pqm1(kbdim,klev),              &
              pum1(kbdim,klev),         pvm1(kbdim,klev),              &
              ptte(kbdim,klev),         pqte(kbdim,klev),              &
              pvom(kbdim,klev),         pvol(kbdim,klev),              &
              pverv(kbdim,klev),        pgeo(kbdim,klev),              &
              papp1(kbdim,klev),        paphp1(kbdim,klevp1)
  REAL(wp)::  prsfc(kbdim),             pssfc(kbdim)
  REAL(wp)::  pthvsig(kbdim)
  INTEGER ::  ktype(kbdim)
  REAL(wp)::  pqhfla(kbdim)
  REAL(wp)::  ptopmax(kbdim)
  INTEGER ::  ilab(kbdim,klev)
  REAL(wp)::  pxtec(kbdim,klev),        pqtec(kbdim,klev)
  REAL(wp)::  pxtecl(kbdim,klev),       pxteci(kbdim,klev)
  REAL(wp)::  pxlm1(kbdim,klev),        pxim1(kbdim,klev),             &
              pxlte(kbdim,klev),        pxite(kbdim,klev)
  
  REAL(wp),INTENT(INOUT) :: pch_concloud(kbdim), pcon_dtrl(kbdim), pcon_dtri(kbdim)
  REAL(wp),INTENT(INOUT) :: pcon_iqte(kbdim)
  REAL(wp),INTENT(INOUT) :: ptte_cnv(kbdim,klev)                               
  REAL(wp),INTENT(INOUT) :: pvom_cnv(kbdim,klev), pvol_cnv(kbdim,klev)         
  REAL(wp),INTENT(INOUT) :: pqte_cnv(kbdim,klev), pxtte_cnv(kbdim,klev,ktrac)  

  REAL(wp)::  ztp1(kbdim,klev),         zqp1(kbdim,klev),              &
              zxp1(kbdim,klev),         ztvp1(kbdim,klev),             &
              zup1(kbdim,klev),         zvp1(kbdim,klev),              &
              ztu(kbdim,klev),          zqu(kbdim,klev),               &
              zlu(kbdim,klev),          zlude(kbdim,klev),             &
              zqude(kbdim,klev),                                       &
              zcpq(kbdim,klev),         zcq(kbdim,klev),               &
              zmfu(kbdim,klev),         zmfd(kbdim,klev),              &
              zqsat(kbdim,klev),        zrain(kbdim),                  &
              za(kbdim),                ua(kbdim)

  INTEGER ::  itopec2(kbdim),           idx(kbdim)
  INTEGER ::  icbot(kbdim),             ictop(kbdim)
  REAL(wp)::  zxtp1(kbdim,klev,ktrac),  zxtu(kbdim,klev,ktrac),        &
              pxtm1(kbdim,klev,ktrac),  pxtte(kbdim,klev,ktrac)
  REAL(wp)::  ztopmax(kbdim)
  LOGICAL ::  locum(kbdim),             ldland(kbdim)

  !  Local scalars: 
  REAL(wp):: ztmst, zxlp1, zxip1
  INTEGER :: ilevmin, jk, jl, jt

  TYPE(t_echam_phy_field),   POINTER :: field

  field  => prm_field(1)

  !  Executable statements 

!
!-----------------------------------------------------------------------
!*    1.           CALCULATE T,Q AND QS AT MAIN LEVELS
!*                 -----------------------------------
!
!
  ztmst=time_step_len
  DO 120 jk=1,klev
!IBM* NOVECTOR
     DO jl=1,kproma
        ztp1(jl,jk)=ptm1(jl,jk)+ptte(jl,jk)*ztmst
        zqp1(jl,jk)=MAX(0._wp,pqm1(jl,jk)+pqte(jl,jk)*ztmst)
        zxlp1=pxlm1(jl,jk)+pxlte(jl,jk)*ztmst
        zxip1=pxim1(jl,jk)+pxite(jl,jk)*ztmst
        zxp1(jl,jk)=MAX(0._wp,zxlp1+zxip1)
        ztvp1(jl,jk)=ztp1(jl,jk)*(1._wp+vtmpc1*zqp1(jl,jk)-zxp1(jl,jk))
        zup1(jl,jk)=pum1(jl,jk)+pvom(jl,jk)*ztmst
        zvp1(jl,jk)=pvm1(jl,jk)+pvol(jl,jk)*ztmst
     END DO

     CALL prepare_ua_index_spline('cucall',kproma,ztp1(1,jk),idx(1),za(1))
     CALL lookup_ua_spline(kproma,idx(1),za(1),ua(1))

!IBM* NOVECTOR
     DO jl=1,kproma
        zqsat(jl,jk)=ua(jl)/papp1(jl,jk)
        zqsat(jl,jk)=MIN(0.5_wp,zqsat(jl,jk))
        zqsat(jl,jk)=zqsat(jl,jk)/(1._wp-vtmpc1*zqsat(jl,jk))
        zcpq(jl,jk)=cpd+(cpv-cpd)*MAX(pqm1(jl,jk),0.0_wp) ! cp of moist air for comp. of fluxes
        zcq (jl,jk)=pcd+(pcv-pcd)*MAX(pqm1(jl,jk),0.0_wp) ! cp or cv of moist air for temp tendency
     END DO

     DO 1104 jt=1,ktrac
        DO 1102 jl=1,kproma
           zxtp1(jl,jk,jt)=pxtm1(jl,jk,jt)+pxtte(jl,jk,jt)*ztmst
1102    END DO
1104 END DO

120 END DO
  DO 130 jl=1,kproma
     zrain(jl)=0._wp
     locum(jl)=.FALSE.
130 END DO

! temporary fields for debugging.  
  DO jk=1,klev
     DO jl=1,kproma
        field%debug_3d_1(jl,jk,krow)=ptm1(jl,jk)        
        field%debug_3d_2(jl,jk,krow)=pqm1(jl,jk)        
        field%debug_3d_2b(jl,jk,krow)=pxlm1(jl,jk)        
        field%debug_3d_2c(jl,jk,krow)=pxim1(jl,jk)        
        field%debug_3d_3(jl,jk,krow)=ptte(jl,jk)        
        field%debug_3d_4(jl,jk,krow)=pqte(jl,jk)        
        field%debug_3d_5(jl,jk,krow)=ztp1(jl,jk)        
        field%debug_3d_6(jl,jk,krow)=zqp1(jl,jk)        
        field%debug_3d_7(jl,jk,krow)=zqsat(jl,jk)        
        field%debug_3d_8(jl,jk,krow)=pverv(jl,jk)        
        field%debug_3d_9(jl,jk,krow)=paphp1(jl,jk)        
        field%debug_3d_10(jl,jk,krow)=papp1(jl,jk)        
        field%debug_3d_11(jl,jk,krow)=pgeo(jl,jk)    
     END DO
  END DO
  DO jl=1,kproma
     field%debug_2d_1(jl,krow)=pqhfla(jl)   
     field%debug_2d_2(jl,krow)=pthvsig(jl)        
  END DO
!
!
!-----------------------------------------------------------------------
!
!*    2.     CALL 'CUMASTR'(MASTER-ROUTINE FOR CUMULUS PARAMETERIZATION)
!*           -----------------------------------------------------------
!
!
  SELECT CASE (iconv)
  CASE(1)
     CALL cumastr(kproma, kbdim, klev, klevp1, klevm1, ilab,           &
!---Included for in-cloud scavenging (Philip Stier, 19/01/06):----------
                  krow,     papp1,                                     &
!---End Included for scavenging-----------------------------------------
                  ztp1,     zqp1,     zxp1,     zup1,   zvp1,          &
                  ztvp1,    ktrac,    ldland,                          &
                  zxtp1,    zxtu,     pxtte,                           &
                  pverv,    zqsat,    pqhfla,                          &
                  paphp1,   pgeo,                                      &
                  ptte,     pqte,     pvom,     pvol,                  &
                  prsfc,    pssfc,    pxtec,                           &
                  pqtec,    zqude,    zcpq,     zcq,                   &
                  pch_concloud, pcon_dtrl, pcon_dtri,                  &
                  pcon_iqte,                                           &
                  locum,    ktype,    icbot,    ictop,                 &
                  ztu,      zqu,      zlu,      zlude,                 &
                  zmfu,     zmfd,     zrain,    pthvsig,               &
                  pxtecl,   pxteci,   nmctop,   cevapcu, time_step_len,&
                  ptte_cnv, pvom_cnv, pvol_cnv, pqte_cnv,pxtte_cnv     )
  CASE(2)
     CALL cumastrt(kproma, kbdim, klev, klevp1, klevm1, ilab,          &
!---Included for in-cloud scavenging (Philip Stier, 19/01/06):----------
                  krow,     papp1,                                     &
!---End Included for scavenging-----------------------------------------
                  ztp1,     zqp1,     zxp1,     zup1,   zvp1,          &
                  ztvp1,    ktrac,    ldland,                          &
                  zxtp1,    zxtu,     pxtte,                           &
                  pverv,    zqsat,    pqhfla,                          &
                  paphp1,   pgeo,                                      &
                  ptte,     pqte,     pvom,     pvol,                  &
                  prsfc,    pssfc,    pxtec,                           &
                  pqtec,    zqude,    zcpq,     zcq,                   &
                  pch_concloud, pcon_dtrl, pcon_dtri,                  &
                  pcon_iqte,                                           &
                  locum,    ktype,    icbot,    ictop,                 &
                  ztu,      zqu,      zlu,      zlude,                 &
                  zmfu,     zmfd,     zrain,    pthvsig,               &
                  pxtecl,   pxteci,   nmctop,   cevapcu, time_step_len,&
                  ptte_cnv, pvom_cnv, pvol_cnv, pqte_cnv,pxtte_cnv     )
  CASE(3)
     CALL cumastrh(kproma, kbdim, klev, klevp1, klevm1, ilab,          &
!---Included for in-cloud scavenging (Philip Stier, 19/01/06):----------
                  krow,     papp1,                                     &
!---End Included for scavenging-----------------------------------------
                  ztp1,     zqp1,     zxp1,     zup1,   zvp1,          &
                  ztvp1,    ktrac,    ldland,                          &
                  zxtp1,    zxtu,     pxtte,                           &
                  pverv,    zqsat,    pqhfla,                          &
                  paphp1,   pgeo,                                      &
                  ptte,     pqte,     pvom,     pvol,                  &
                  prsfc,    pssfc,    pxtec,                           &
                  pqtec,    zqude,    zcpq,     zcq,                   &
                  pch_concloud, pcon_dtrl, pcon_dtri,                  &
                  pcon_iqte,                                           &
                  locum,    ktype,    icbot,    ictop,                 &
                  ztu,      zqu,      zlu,      zlude,                 &
                  zmfu,     zmfd,     zrain,    pthvsig,               &
                  pxtecl,   pxteci,   nmctop,   cevapcu, time_step_len,&
                  ptte_cnv, pvom_cnv, pvol_cnv, pqte_cnv,pxtte_cnv     )
  END SELECT
!
!
! ------------------------------------------------------------------
!
!*     3.     PRESSURE ALTITUDE OF CONVECTIVE CLOUD TOPS.
!             -------- -------- -- ---------- ----- -----
!
  ilevmin=klev-4
!
  DO 301 jl=1,kproma
     itopec2(jl)=klevp1
301 END DO
!
  DO 303 jk=1,ilevmin
     DO 302 jl=1,kproma
        IF(ilab(jl,jk).EQ.2 .AND. itopec2(jl).EQ.klevp1) THEN
           itopec2(jl)=jk
        END IF
302  END DO
303 END DO
!
  ztopmax(1:kproma) = ptopmax(1:kproma)

  DO 304 jl=1,kproma
     IF(itopec2(jl).EQ.1) THEN
        ptopmax(jl)=papp1(jl,1)
     ELSE IF(itopec2(jl).NE.klevp1) THEN
        ptopmax(jl)=paphp1(jl,itopec2(jl))
     ELSE
        ptopmax(jl)=99999._wp
     END IF
     ptopmax(jl)=MIN(ptopmax(jl),ztopmax(jl))
304 END DO
!
!
!---------------------------------------------------------------------
!
  RETURN
END SUBROUTINE cucall

SUBROUTINE cumastr(  kproma, kbdim, klev, klevp1, klevm1, ilab,         &
                      krow,    papp1,                                   &
           pten,     pqen,     pxen,     puen,     pven,                &
           ptven,    ktrac,    ldland,                                  &
           pxten,    pxtu,     pxtte,                                   &
           pverv,    pqsen,    pqhfla,                                  &
           paphp1,   pgeo,                                              &
           ptte,     pqte,     pvom,     pvol,                          &
           prsfc,    pssfc,    pxtec,                                   &
           pqtec,    pqude,    zcpen,    zcen,                          &
           pch_concloud, pcon_dtrl, pcon_dtri, pcon_iqte,               &
           ldcum,    ktype,    kcbot,    kctop,                         &
           ptu,      pqu,      plu,      plude,                         &
           pmfu,     pmfd,     prain,    pthvsig,                       &
           pxtecl,   pxteci,   nmctop,   cevapcu,  time_step_len,       &        
           ptte_cnv, pvom_cnv, pvol_cnv, pqte_cnv, pxtte_cnv     )
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
INTEGER, INTENT (IN) :: kproma, kbdim, klev, klevp1, ktrac, klevm1
REAL(wp),INTENT(IN)  :: cevapcu(:), time_step_len
INTEGER, INTENT(IN)  :: nmctop
REAL(wp),INTENT(INOUT):: pch_concloud(kbdim)
REAL(wp),INTENT(INOUT):: pcon_dtrl(kbdim), pcon_dtri(kbdim)
REAL(wp),INTENT(INOUT):: pcon_iqte(kbdim)
REAL(wp),INTENT(INOUT):: ptte_cnv(kbdim,klev)                               
REAL(wp),INTENT(INOUT):: pvom_cnv(kbdim,klev), pvol_cnv(kbdim,klev)         
REAL(wp),INTENT(INOUT):: pqte_cnv(kbdim), pxtte_cnv(kbdim,klev,ktrac)  
!---Included for in-cloud scavenging (Philip Stier, 19/01/06):----------
INTEGER, INTENT (IN) :: krow
!---End Included for scavenging-----------------------------------------
!
REAL(wp):: pten(kbdim,klev),        pqen(kbdim,klev),                  &
           pxen(kbdim,klev),        ptven(kbdim,klev),                 &
           puen(kbdim,klev),        pven(kbdim,klev),                  &
           ptte(kbdim,klev),        pqte(kbdim,klev),                  &
           pvom(kbdim,klev),        pvol(kbdim,klev),                  &
           pqsen(kbdim,klev),       pgeo(kbdim,klev),                  &
           paphp1(kbdim,klevp1),                                       &
           pverv(kbdim,klev)                                            
REAL(wp):: ptu(kbdim,klev),         pqu(kbdim,klev),                   &
           plu(kbdim,klev),         plude(kbdim,klev),                 &
           pmfu(kbdim,klev),        pmfd(kbdim,klev),                  &
           prsfc(kbdim),            pssfc(kbdim),                      &
           prain(kbdim),            pqhfla(kbdim)
REAL(wp):: pthvsig(kbdim)
INTEGER :: kcbot(kbdim),            kctop(kbdim),                      &
           ktype(kbdim)
REAL(wp):: pxtec(kbdim,klev),       pqtec(kbdim,klev),                 &
           pxtecl(kbdim,klev),      pxteci(kbdim,klev),                &
           pqude(kbdim,klev)
REAL(wp):: ztenh(kbdim,klev),       zqenh(kbdim,klev),                 &
           zxenh(kbdim,klev),       zalvsh(kbdim,klev),                &
! zalvsh: latent heat of vaporisation/sublimation defined at half levels
           zgeoh(kbdim,klev),       zqsenh(kbdim,klev),                &
           ztd(kbdim,klev),         zqd(kbdim,klev),                   &
           zmfus(kbdim,klev),       zmfds(kbdim,klev),                 &
           zmfuq(kbdim,klev),       zmfdq(kbdim,klev),                 &
           zdmfup(kbdim,klev),      zdmfdp(kbdim,klev),                &
           zmful(kbdim,klev),       zrfl(kbdim),                       &
           zuu(kbdim,klev),         zvu(kbdim,klev),                   &
           zud(kbdim,klev),         zvd(kbdim,klev)
REAL(wp):: zcpen(kbdim,klev),       zcen(kbdim,klev),                  &
           zcpcu(kbdim,klev)
REAL(wp):: zentr(kbdim),            zhcbase(kbdim),                    &
           zmfub(kbdim),            zmfub1(kbdim),                     &
           zktype(kbdim),           zldcum(kbdim),                     &
           zcpcui(kbdim,klev),      zkcbot(kbdim),                     &
           zictop0(kbdim),          ztmp1(kbdim),                      &
           ztmp2(kbdim),            ztmp3(kbdim),                      &
           za(kbdim)
REAL(wp):: zsfl(kbdim),             zdpmel(kbdim,klev)
REAL(wp):: zcape(kbdim),            zheat(kbdim)
REAL(wp):: ua(kbdim), dua(kbdim), ub(kbdim)
REAL(wp):: zhmin(kbdim),            zihmin(kbdim)
REAL(wp):: zhhatt(kbdim,klev),      zdqpbl(kbdim)
REAL(wp):: zdqcv(kbdim)
INTEGER :: ihmin(kbdim), ilo1(kbdim), ldidx(kbdim), loidx(kbdim)
INTEGER :: ilab(kbdim,klev),        idtop(kbdim),                      &
           ictop0(kbdim),           ilwmin(kbdim)
REAL(wp):: pxten(kbdim,klev,ktrac), pxtte(kbdim,klev,ktrac),           &
           pxtu(kbdim,klev,ktrac),                                     &
           zxtenh(kbdim,klev,ktrac),zxtd(kbdim,klev,ktrac),            &
           zmfuxt(kbdim,klev,ktrac),zmfdxt(kbdim,klev,ktrac)
LOGICAL :: loddraf(kbdim),          ldland(kbdim)
LOGICAL :: ldcum(kbdim)
LOGICAL :: llo1
!---Included for scavenging in wetdep_interface (Philip Stier, 19/02/04, UL, 28.3.07):-------
REAL(wp):: papp1(kbdim,klev)
REAL(wp):: zmwc(kbdim,klev),        zmrateprecip(kbdim,klev)
!---End Included for scavenging-----------------------------------------
!
INTEGER :: nl, jl, jk, ikb, jt, itopm2, locnt, ldcnt
REAL(wp):: zcons2, ztau, zqumqe, zdqmin, zmfmax, zalvs, zalvdcp, zqalv &
         , zhsat, zes, zcor, zqsat, zdqsdt, zgam, zzz, zhhat           &
         , zb, zbi, zroi, zdz, zdhdz, zdepth, zfac, zrh, zeps          &
         , zjk, zhelp, zlo1, zpaphp1i, ztenhi, zkctop, zpbmpt, za1, za2
INTEGER :: icuasc
!
!  INTRINSIC FUNCTIONS
INTRINSIC MIN, MAX
!
!  Executable statements

!  IF (lconvmassfix) THEN
!     CALL xt_conv_massfix(kproma,         kbdim,         klev,              &
!                          klevp1,         ktrac,         krow,              &
!                          papp1,          paphp1,        pxtte,             &
!                          .TRUE.  )
!  END IF

!-----------------------------------------------------------------------
!
!     1.           SPECIFY CONSTANTS AND PARAMETERS
!                  --------------------------------
!
  zcons2=1._wp/(grav*time_step_len)
  ztau=cmftau
!  ztau=MIN(3._wp*3600._wp,7200._wp*63._wp/nn)
!
!----------------------------------------------------------------------
!
!*    2.           INITIALIZE VALUES AT VERTICAL GRID POINTS IN 'CUINI'
!                  ---------------------------------------------------
!
  CALL cuini(kproma, kbdim, klev, klevp1, klevm1,                      &
             pten,     pqen,     pqsen,    pxen,     puen,     pven,   &
             ptven,    ktrac,                                          &
             pxten,    zxtenh,   pxtu,     zxtd,     zmfuxt,   zmfdxt, &
             pverv,    pgeo,     paphp1,   zgeoh,                      &
             ztenh,    zqenh,    zqsenh,   zxenh,    ilwmin,           &
             ptu,      pqu,      ztd,      zqd,                        &
             zuu,      zvu,      zud,      zvd,                        &
             pmfu,     pmfd,     zmfus,    zmfds,                      &
             zmfuq,    zmfdq,    zdmfup,   zdmfdp,                     &
             zcpen,    zcpcu,    zalvsh,                               &
             zdpmel,   plu,      plude,    pqude,    ilab             )
!
!-----------------------------------------------------------------------
!
!*    3.0          CLOUD BASE CALCULATIONS
!                  -----------------------
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
  zkcbot(1:kproma) = REAL(kcbot(1:kproma),wp)

  jk=1
  DO 310 jl=1,kproma
     zdqpbl(jl)=0.0_wp
     zdqcv(jl)=pqte(jl,jk)*(paphp1(jl,jk+1)-paphp1(jl,jk))
     idtop(jl)=0
310 END DO
  DO 320 jk=2,klev
     zjk = REAL(jk,wp)
     DO 315 jl=1,kproma
        zhelp      = paphp1(jl,jk+1)-paphp1(jl,jk)
        zdqcv(jl)  = zdqcv(jl)+pqte(jl,jk)*zhelp
        zdqpbl(jl) = zdqpbl(jl) + FSEL(zjk - zkcbot(jl),pqte(jl,jk),0._wp)*zhelp
#ifdef __ibmdbg__
        PRINT '(A6,I3,I4,E18.10,L3)','zdqpbl',jk,jl,zdqpbl(jl),(FSEL(zjk - zkcbot(jl),1._wp,0._wp)==1._wp)
#endif
315  END DO
320 END DO
!
!*             (C) DETERMINE MOISTURE SUPPLY FOR BOUNDARY LAYER
!*                 AND DETERMINE CLOUD BASE MASSFLUX IGNORING
!*                 THE EFFECTS OF DOWNDRAFTS AT THIS STAGE
!                  ------------------------------------------
!
  zldcum(1:kproma) = MERGE(1._wp,0._wp,ldcum(1:kproma))

!DIR$ IVDEP
!DIR$ CONCURRENT
!IBM* ASSERT(NODEPS)
!IBM* NOVECTOR
  DO 340 jl=1,kproma
     ikb=kcbot(jl)
     zqumqe=pqu(jl,ikb)+plu(jl,ikb)-zqenh(jl,ikb)
     zdqmin=MAX(0.01_wp*zqenh(jl,ikb),1.e-10_wp)
     zlo1 = FSEL(-zdqpbl(jl),0._wp,1._wp)
     zlo1 = FSEL(zdqmin - zqumqe,0._wp,zlo1) * zldcum(jl)
     zmfub(jl)=FSEL(-zlo1,0.01_wp,(zdqpbl(jl)/(grav*MAX(zqumqe,zdqmin))))
     zmfmax=(paphp1(jl,ikb)-paphp1(jl,ikb-1))*zcons2
     zmfub(jl)=MIN(zmfub(jl),zmfmax)
     zldcum(jl) = zlo1
     zhelp = MAX(0._wp,-1.1_wp*pqhfla(jl)*grav)
     zktype(jl) = FSEL(zhelp - zdqcv(jl),2._wp,1._wp)
     zentr(jl)  = FSEL(zhelp - zdqcv(jl),entrscv,entrpen)
     ktype(jl) = INT(zktype(jl))
#ifdef __ibmdbg__
     PRINT '(A6,I3,I4,2 E18.10,I3,L3)','zmfub',ikb,jl,zmfub(jl),zentr(jl),ktype(jl),(zlo1==1._wp)
#endif
340 END DO
  ldcum(1:kproma) = (zldcum(1:kproma).GT.0._wp)
!
!-----------------------------------------------------------------------
!*    4.0          DETERMINE CLOUD ASCENT FOR ENTRAINING PLUME
!                  -------------------------------------------
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
     zictop0(jl) = zkcbot(jl)-1._wp
410 END DO
  DO 430 jk=klev,1,-1
     zcpcui(1:kproma,jk) = 1._wp/zcpcu(1:kproma,jk)
     IF (jk <= klevm1 .AND. jk >= 3) THEN
        ! mpuetz: too few instructions (FP dependencies)
        CALL prepare_ua_index_spline('cumastr',kproma,ztenh(1,jk),loidx(1),za(1))
        CALL lookup_ua_spline(kproma,loidx(1),za(1),ua(1),dua(1))
        CALL lookup_ubc('cumastr',kproma,ztenh(1,jk),ub(1))
        zjk = REAL(jk,wp)
!IBM* NOVECTOR
        DO 421 jl=1,kproma
           ! mpuetz: move some of these into the previous loop
           zalvs=FSEL(tmelt - ztenh(jl,jk),als,alv)
           zalvdcp=zalvs*zcpcui(jl,jk)
           zqalv=1._wp/zalvs
           zhsat=zcpcu(jl,jk)*ztenh(jl,jk)+zgeoh(jl,jk)+zalvs*zqsenh(jl,jk)
           zpaphp1i = 1._wp/paphp1(jl,jk)
           zes=ua(jl)*zpaphp1i
           zes=MIN(0.5_wp,zes)
           zcor=1._wp/(1._wp-vtmpc1*zes)
           zqsat=zes*zcor
           zdqsdt=zpaphp1i*zcor**2*dua(jl)
           zgam=FSEL(zes - 0.4_wp,zqsat*zcor*ub(jl),zalvdcp*zdqsdt)
           zzz=zcpcu(jl,jk)*ztenh(jl,jk)*vtmpc1
           zhhat=zhsat-((zzz+zgam*zzz)/(1._wp+zgam*zzz*zqalv))*             &
                MAX(zqsenh(jl,jk)-zqenh(jl,jk),0._wp)
           zhhatt(jl,jk)=zhhat
           zlo1 = FSEL(zjk - zictop0(jl),0._wp,1._wp)
           zlo1 = FSEL(zhhat - zhcbase(jl),0._wp,zlo1)
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
!                paphp1(jl,kcbot(jl))-paphp1(jl,ictop0(jl)).gt.2.E4_wp)
!     zentr(jl)=MERGE(entrpen,entrscv,ktype(jl).eq.1)
!  ENDDO
!!
!
!                  FIND LOWEST POSSIBLE ORG. DETRAINMENT LEVEL
!                  -------------------------------------------
!
  ldcnt = 0
  DO jl=1,kproma
     zhmin(jl)=0._wp
     zihmin(jl)=0._wp
     llo1=ldcum(jl).AND.ktype(jl).EQ.1
     IF(llo1) THEN
        zihmin(jl)=zkcbot(jl)
        ldcnt = ldcnt + 1
        ldidx(ldcnt) = jl
     ENDIF
  ENDDO
!
  zb=25._wp
  zbi=1._wp/(zb*grav)
  DO jk=klev,1,-1

     ! mpuetz: compute the update criterion

     zjk = REAL(jk,wp)
!IBM* ASSERT(NODEPS)
     DO nl = 1,ldcnt
        jl = ldidx(nl)
        zlo1 = FSEL(zjk - zkcbot(jl),0._wp,1._wp)
        zlo1 = FSEL(zjk - zictop0(jl),zlo1,0._wp)
        zlo1 = FSEL(-ABS(zihmin(jl)-zkcbot(jl)),zlo1,0._wp)
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

!IBM* ASSERT(NODEPS)
     DO nl=1,locnt
        jl = loidx(nl)

        zalvs = FSEL(tmelt - ztenh(jl,jk),als,alv)
!        zalvs=MERGE(alv,als,ztenh(jl,jk)>tmelt)
        ikb   = kcbot(jl)
        zroi  = SWDIV_NOCHK(rd*ztenh(jl,jk)*(1._wp+vtmpc1*zqenh(jl,jk)),paphp1(jl,jk))
        zdz   = (paphp1(jl,jk)-paphp1(jl,jk-1))*zroi/grav

        za1 = (zcpen(jl,jk-1)*pten(jl,jk-1)        &
             - zcpen(jl,jk)*pten(jl,jk)            &
             + zalvs*(pqen(jl,jk-1) - pqen(jl,jk)) &
             + (pgeo(jl,jk-1)-pgeo(jl,jk)) )*grav
        za2 = pgeo(jl,jk-1)-pgeo(jl,jk)
        zdhdz = SWDIV_NOCHK(za1, za2)

        zdepth    = zgeoh(jl,jk)-zgeoh(jl,ikb)

        ztmp1(nl) = zalvs
        ztmp2(nl) = zdz*zdhdz
        ztmp3(nl) = 1._wp+zdepth*zbi
     END DO

     ztmp3(1:locnt) = SQRT(ztmp3(1:locnt))

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
  CALL cuasc(kproma, kbdim, klev, klevp1, klevm1,                      &
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
             kcbot,    kctop,    ictop0,   nmctop,                     &
!---Included for scavenging in wetdep_interface (Philip Stier, 19/02/04, UL, 28.3.07):-------
             zmwc,     zmrateprecip,       time_step_len               )
!---End Included for scavenging-----------------------------------------
!
!*         (C) CHECK CLOUD DEPTH AND CHANGE ENTRAINMENT RATE ACCORDINGLY
!              CALCULATE PRECIPITATION RATE (FOR DOWNDRAFT CALCULATION)
!              ---------------------------------------------------------
!
!DIR$ IVDEP
  DO 440 jl=1,kproma
     zpbmpt=paphp1(jl,kcbot(jl))-paphp1(jl,kctop(jl))
     IF(ldcum(jl).AND.ktype(jl).EQ.1.AND.zpbmpt.LT.2.e4_wp) ktype(jl)=2
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
  IF(lmfdd) THEN
!
!*             (A) DETERMINE LFS IN 'CUDLFS'
!                  -------------------------
!
     CALL cudlfs(kproma,   kbdim,    klev,     klevp1,                 &
                 ztenh,    zqenh,    puen,     pven,                   &
                 ktrac,                                                &
                 zxtenh,   pxtu,     zxtd,     zmfdxt,                 &
                 zgeoh,    paphp1,                                     &
                 ptu,      pqu,      zuu,      zvu,                    &
                 ldcum,    kcbot,    kctop,    zmfub,    zrfl,         &
                 ztd,      zqd,      zud,      zvd,                    &
                 pmfd,     zmfds,    zmfdq,    zdmfdp,                 &
                 zcpcu,                                                &
                 idtop,    loddraf)
!
!*            (B)  DETERMINE DOWNDRAFT T,Q AND FLUXES IN 'CUDDRAF'
!                  -----------------------------------------------
!
     CALL cuddraf(kproma,   kbdim,    klev,     klevp1,                &
                  ztenh,    zqenh,    puen,     pven,                  &
                  ktrac,                                               &
                  zxtenh,   zxtd,     zmfdxt,                          &
                  zgeoh,    paphp1,   zrfl,                            &
                  ztd,      zqd,      zud,      zvd,                   &
                  pmfd,     zmfds,    zmfdq,    zdmfdp,                &
                  zcpcu,                                               &
                  loddraf)
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

  zkcbot(1:kproma) = REAL(kcbot(1:kproma),wp)

  DO jl=1,kproma
     zheat(jl)=0._wp
     zcape(jl)=0._wp
     zmfub1(jl)=zmfub(jl)
!     zhelp(jl) =0._wp
  ENDDO
!
  DO jk=1,klev
     zjk = REAL(jk,wp)
!IBM* ASSERT(NODEPS)
     DO nl = 1,ldcnt
        jl = ldidx(nl)
        zkctop = REAL(kctop(jl),wp)
        zlo1 = FSEL(zkcbot(jl) - zjk,1._wp,0._wp) &
             * FSEL(zkctop - zjk,0._wp,1._wp)
        ilo1(nl) = INT(zlo1)
     END DO

     locnt = 0
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

!IBM* ASSERT(NODEPS)
     DO nl=1,locnt
        jl = loidx(nl)

        zroi      = SWDIV_NOCHK(rd*ztenh(jl,jk)*(1._wp+vtmpc1*zqenh(jl,jk)),paphp1(jl,jk))
        ztenhi    = SWDIV_NOCHK(1._wp,ztenh(jl,jk))
        zdz       = (paphp1(jl,jk)-paphp1(jl,jk-1))*zroi/grav
        zheat(jl) = zheat(jl) +                                           &
             (  (pten(jl,jk-1)-pten(jl,jk) + grav*zdz*zcpcui(jl,jk))*ztenhi  &
             +  vtmpc1*(pqen(jl,jk-1)-pqen(jl,jk)) )                      &
             * (grav*(pmfu(jl,jk)+pmfd(jl,jk)))*zroi
        zcape(jl) = zcape(jl) +                                           &
             ( grav*(ptu(jl,jk)-ztenh(jl,jk))*ztenhi                         &
             + grav*vtmpc1*(pqu(jl,jk)-zqenh(jl,jk))                         &
             - grav*plu(jl,jk) ) * zdz

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
!             (1._wp+vtmpc1*zqenh(jl,jk)))
!        zdz=(paphp1(jl,jk)-paphp1(jl,jk-1))/(grav*zro)
!        zhelp(jl)=zhelp(jl) +                               &
!             (grav*(ptu(jl,jk)-ztenh(jl,jk))/ztenh(jl,jk)     &
!             +grav*vtmpc1*(pqu(jl,jk)-zqenh(jl,jk))      &
!             -grav*plu(jl,jk) ) * zdz
!     ENDDO
!     if (zhelp(jl).lt.0._wp) zhelp(jl)=0._wp
!  ENDDO

  DO jl=1,kproma
     IF(ldcum(jl).AND.ktype(jl).EQ.1) THEN
       ikb=kcbot(jl)
       zmfub1(jl) = SWDIV_NOCHK((zcape(jl)*zmfub(jl)),(zheat(jl)*ztau))
       zmfub1(jl) = MAX(zmfub1(jl),0.001_wp)
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
!IBM* ASSERT(NODEPS)
  DO 520 nl=1,ldcnt
     jl = ldidx(nl)
     ikb=kcbot(jl)
     llo1=pmfd(jl,ikb).LT.0._wp.AND.loddraf(jl)
     zeps=MERGE(cmfdeps,0._wp,llo1)
     zqumqe=pqu(jl,ikb)+plu(jl,ikb)-                                &
          zeps*zqd(jl,ikb)-(1._wp-zeps)*zqenh(jl,ikb)
     zdqmin=MAX(0.01_wp*zqenh(jl,ikb),1.e-10_wp)
     zmfmax=(paphp1(jl,ikb)-paphp1(jl,ikb-1))*zcons2
     llo1=zdqpbl(jl).GT.0._wp.AND.zqumqe.GT.zdqmin.AND.ldcum(jl)    &
          .AND.zmfub(jl).LT.zmfmax
     zmfub1(jl)=MERGE(zdqpbl(jl)/(grav*MAX(zqumqe,zdqmin)),            &
          zmfub(jl),llo1)
     zmfub1(jl)=MERGE(zmfub1(jl),zmfub(jl),                         &
          ABS(zmfub1(jl)-zmfub(jl)).LT.0.2_wp*zmfub(jl))
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
!IBM* ASSERT(NODEPS)
     DO 530 nl=1,ldcnt
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
530  END DO
!
!IBM* unroll(4)
     DO 5304 jt=1,ktrac
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
  CALL cuasc(kproma, kbdim, klev, klevp1, klevm1,                      &
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
             kcbot,    kctop,    ictop0,   nmctop,                     &
!---Included for scavenging in wetdep_interface (Philip Stier, 19/02/04, UL, 28.3.07):-------
             zmwc,     zmrateprecip,       time_step_len               )
!---End Included for scavenging-----------------------------------------
!-----------------------------------------------------------------------
!*    7.0          DETERMINE FINAL CONVECTIVE FLUXES IN 'CUFLX'
!                  ------------------------------------------
!
  CALL cuflx(kproma,   kbdim,    klev,     klevp1,                     &
             pqen,     pqsen,    ztenh,    zqenh,                      &
             ktrac,                                                    &
!---Included for scavenging in wetdep_interface (Philip Stier, 28/03/01, UL, 28.3.07):-------
             krow,                                                     &
             pxtte,    pxtu,     ptu,                                  &
             zmwc,     zmrateprecip,                                   &
!---End Included for scavenging-----------------------------------------
             zxtenh,   zmfuxt,   zmfdxt,                               &
             paphp1,   zgeoh,                                          &
             kcbot,    kctop,    idtop,                                &
             ktype,    loddraf,  ldcum,                                &
             pmfu,     pmfd,     zmfus,    zmfds,                      &
             zmfuq,    zmfdq,    zmful,                                &
             zdmfup,   zdmfdp,   zrfl,     prain,                      &
             zcpcu,                                                    &
             pten,     zsfl,     zdpmel,   itopm2,  cevapcu,           &
             time_step_len                                             )
!-----------------------------------------------------------------------
!
!*    8.0          UPDATE TENDENCIES FOR T AND Q IN SUBROUTINE CUDTDQ
!                  --------------------------------------------------
!
  CALL cudtdq(kproma, kbdim, klev, klevp1, itopm2, ldcum, ktrac,       &
!--- Included for dust emissions (Philip Stier 23/01/06)-----------------
              krow,                                                    &
!--- End Included for dust emissions in ---------------------------------
              paphp1,   pten,     ptte,     pqte,                      &
              pxtte,    pxtec,    zmfuxt,   zmfdxt,                    &
              zmfus,    zmfds,    zmfuq,    zmfdq,                     &
              zmful,    zdmfup,   zdmfdp,   plude,                     &
              zdpmel,   zrfl,     zsfl,                                &
              zcen,     zalvsh,   pqtec,    pqude,                     &
              prsfc,    pssfc,                                         &
              pch_concloud, pcon_dtrl, pcon_dtri,                      &
              pxtecl,   pxteci, pcon_iqte,                             &
              ptte_cnv, pqte_cnv, pxtte_cnv                        )
!
!-----------------------------------------------------------------------
!
!*    9.0          UPDATE TENDENCIES FOR U AND U IN SUBROUTINE CUDUDV
!                  --------------------------------------------------
!
  IF(lmfdudv) THEN
     CALL cududv(kproma,   kbdim,    klev,     klevp1,                 &
                 itopm2,   ktype,    kcbot,    paphp1,   ldcum,        &
                 puen,     pven,     pvom,     pvol,                   &
                 zuu,      zud,      zvu,      zvd,                    &
                 pmfu,     pmfd,     pvom_cnv, pvol_cnv                )
!
  END IF
!
!  IF (lconvmassfix) THEN
!     CALL xt_conv_massfix(kproma,         kbdim,         klev,              &
!                          klevp1,         ktrac,         krow,              &
!                          papp1,          paphp1,        pxtte,             &
!                          .FALSE.  )
!  END IF
!
  RETURN
END SUBROUTINE cumastr

SUBROUTINE cumastrt( kproma, kbdim, klev, klevp1, klevm1, ilab,        &
!---Included for in-cloud scavenging (Philip Stier, 19/01/06):----------
                     krow,    papp1,                                   &
!---End Included for scavenging-----------------------------------------
           pten,     pqen,     pxen,     puen,     pven,               &
           ptven,    ktrac,    ldland,                                 &
           pxten,    pxtu,     pxtte,                                  &
           pverv,    pqsen,    pqhfla,                                 &
           paphp1,   pgeo,                                             &
           ptte,     pqte,     pvom,     pvol,                         &
           prsfc,    pssfc,    pxtec,                                  &
           pqtec,    pqude,    zcpen,    zcen,                         &
           pch_concloud, pcon_dtrl, pcon_dtri, pcon_iqte,              &
           ldcum,    ktype,    kcbot,    kctop,                        &
           ptu,      pqu,      plu,      plude,                        &
           pmfu,     pmfd,     prain,    pthvsig,                      &
           pxtecl,   pxteci,   nmctop,   cevapcu,  time_step_len,      &        
           ptte_cnv, pvom_cnv, pvol_cnv, pqte_cnv, pxtte_cnv     )
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
INTEGER, INTENT (IN) :: kproma, kbdim, klev, klevp1, ktrac, klevm1
REAL(wp),INTENT(IN)  :: cevapcu(:), time_step_len
INTEGER, INTENT(IN)  :: nmctop
REAL(wp),INTENT(INOUT):: pch_concloud(kbdim)
REAL(wp),INTENT(INOUT):: pcon_dtrl(kbdim), pcon_dtri(kbdim)
REAL(wp),INTENT(INOUT):: pcon_iqte(kbdim)
REAL(wp),INTENT(INOUT):: ptte_cnv(kbdim,klev)                               
REAL(wp),INTENT(INOUT):: pvom_cnv(kbdim,klev), pvol_cnv(kbdim,klev)         
REAL(wp),INTENT(INOUT):: pqte_cnv(kbdim), pxtte_cnv(kbdim,klev,ktrac)  
!---Included for in-cloud scavenging (Philip Stier, 19/01/06):----------
INTEGER, INTENT (IN) :: krow
!---End Included for scavenging-----------------------------------------
!
REAL(wp):: pten(kbdim,klev),        pqen(kbdim,klev),                  &
           pxen(kbdim,klev),        ptven(kbdim,klev),                 &
           puen(kbdim,klev),        pven(kbdim,klev),                  &
           ptte(kbdim,klev),        pqte(kbdim,klev),                  &
           pvom(kbdim,klev),        pvol(kbdim,klev),                  &
           pqsen(kbdim,klev),       pgeo(kbdim,klev),                  &
           paphp1(kbdim,klevp1),                                       &
           pverv(kbdim,klev)
REAL(wp):: ptu(kbdim,klev),         pqu(kbdim,klev),                   &
           plu(kbdim,klev),         plude(kbdim,klev),                 &
           pmfu(kbdim,klev),        pmfd(kbdim,klev),                  &
           prsfc(kbdim),            pssfc(kbdim),                      &
           prain(kbdim),            pqhfla(kbdim)
REAL(wp):: pthvsig(kbdim)
INTEGER :: kcbot(kbdim),            kctop(kbdim),                      &
           ktype(kbdim)
REAL(wp):: pxtec(kbdim,klev),       pqtec(kbdim,klev),                 &
           pxtecl(kbdim,klev),      pxteci(kbdim,klev),                &  
           pqude(kbdim,klev)
REAL(wp):: ztenh(kbdim,klev),       zqenh(kbdim,klev),                 &
           zxenh(kbdim,klev),       zalvsh(kbdim,klev),                &
! zalvsh: latent heat of vaporisation/sublimation defined at half levels
           zgeoh(kbdim,klev),       zqsenh(kbdim,klev),                &
           ztd(kbdim,klev),         zqd(kbdim,klev),                   &
           zmfus(kbdim,klev),       zmfds(kbdim,klev),                 &
           zmfuq(kbdim,klev),       zmfdq(kbdim,klev),                 &
           zdmfup(kbdim,klev),      zdmfdp(kbdim,klev),                &
           zmful(kbdim,klev),       zrfl(kbdim),                       &
           zuu(kbdim,klev),         zvu(kbdim,klev),                   &
           zud(kbdim,klev),         zvd(kbdim,klev)
REAL(wp):: zcpen(kbdim,klev),       zcen(kbdim,klev),                  &
           zcpcu(kbdim,klev)
REAL(wp):: zentr(kbdim),            zhcbase(kbdim),                    &
           zmfub(kbdim),            zmfub1(kbdim),                     &
           zcpcui(kbdim,klev),      zictop0(kbdim),                    &
           zdqpbl(kbdim),           zdqcv(kbdim)
REAL(wp):: zsfl(kbdim),             zdpmel(kbdim,klev)
REAL(wp):: ua(kbdim), dua(kbdim), ub(kbdim), za(kbdim)
REAL(wp):: zhhatt(kbdim,klev)
INTEGER :: ilab(kbdim,klev),        idtop(kbdim),                      &
           ictop0(kbdim),           ilwmin(kbdim), loidx(kbdim)
REAL(wp):: pxten(kbdim,klev,ktrac), pxtte(kbdim,klev,ktrac),           &
           pxtu(kbdim,klev,ktrac),                                     &
           zxtenh(kbdim,klev,ktrac),zxtd(kbdim,klev,ktrac),            &
           zmfuxt(kbdim,klev,ktrac),zmfdxt(kbdim,klev,ktrac)
LOGICAL :: loddraf(kbdim),          ldland(kbdim)
LOGICAL :: ldcum(kbdim)
LOGICAL :: llo1
!---Included for scavenging in wetdep_interface (Philip Stier, 19/02/04):-------
REAL(wp):: papp1(kbdim,klev)
REAL(wp):: zmwc(kbdim,klev),        zmrateprecip(kbdim,klev)
!---End Included for scavenging-----------------------------------------
!
INTEGER :: jl, jk, ikb, jt, itopm2
REAL(wp):: zcons2, zqumqe, zdqmin, zmfmax, zalvs, zalvdcp, zqalv       &
         , zhsat, zes, zcor, zqsat, zdqsdt, zgam, zzz, zhhat           &
         , zfac, zpbmpt, zeps, zjk, zpaphp1i, zlo1
!
!  INTRINSIC FUNCTIONS
INTRINSIC MIN, MAX
!
!  Executable statements

!  IF(lconvmassfix) THEN
!     CALL xt_conv_massfix(kproma,         kbdim,            klev,    &
!                          klevp1,         ktrac,            krow,    &
!                          papp1,          paphp1,           pxtte,   &
!                          .TRUE.                                     )
!  END IF
!---------------------------------------------------------------------
!
!     1.           SPECIFY CONSTANTS AND PARAMETERS
!                  --------------------------------
!
  zcons2=1._wp/(grav*time_step_len)
!
!---------------------------------------------------------------------
!
!*    2.           INITIALIZE VALUES AT VERTICAL GRID POINTS IN 'CUINI'
!                  ---------------------------------------------------
!
  CALL cuini(kproma,   kbdim,    klev,     klevp1,   klevm1,           &
             pten,     pqen,     pqsen,    pxen,     puen,     pven,   &
             ptven,    ktrac,                                          &
             pxten,    zxtenh,   pxtu,     zxtd,     zmfuxt,   zmfdxt, &
             pverv,    pgeo,     paphp1,   zgeoh,                      &
             ztenh,    zqenh,    zqsenh,   zxenh,    ilwmin,           &
             ptu,      pqu,      ztd,      zqd,                        &
             zuu,      zvu,      zud,      zvd,                        &
             pmfu,     pmfd,     zmfus,    zmfds,                      &
             zmfuq,    zmfdq,    zdmfup,   zdmfdp,                     &
             zcpen,    zcpcu,    zalvsh,                               &
             zdpmel,   plu,      plude,    pqude,    ilab             )
!
!---------------------------------------------------------------------
!
!*    3.0          CLOUD BASE CALCULATIONS
!                  -----------------------
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
     zdqpbl(jl)=0.0_wp
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
     zdqmin=MAX(0.01_wp*zqenh(jl,ikb),1.e-10_wp)
     llo1=zdqpbl(jl).GT.0._wp.AND.zqumqe.GT.zdqmin.AND.ldcum(jl)
     zmfub(jl)=MERGE(zdqpbl(jl)/(grav*MAX(zqumqe,zdqmin)),0.01_wp,llo1)
     zmfmax=(paphp1(jl,ikb)-paphp1(jl,ikb-1))*zcons2
     zmfub(jl)=MIN(zmfub(jl),zmfmax)
     IF(.NOT.llo1) ldcum(jl)=.FALSE.
     ktype(jl)=MERGE(1,2,zdqcv(jl).GT.MAX(0._wp,-1.1_wp*pqhfla(jl)*grav))
     zentr(jl)=MERGE(entrpen,entrscv,ktype(jl).EQ.1)
340 END DO
!
!---------------------------------------------------------------------
!*    4.0          DETERMINE CLOUD ASCENT FOR ENTRAINING PLUME
!                  -------------------------------------------
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
  DO 430 jk=klev,1,-1
     zcpcui(1:kproma,jk) = 1._wp/zcpcu(1:kproma,jk)
     IF (jk <= klevm1 .AND. jk >= 3) THEN
        ! mpuetz: too few instructions (FP dependencies)
        CALL prepare_ua_index_spline('cumastr',kproma,ztenh(1,jk),loidx(1),za(1))
        CALL lookup_ua_spline(kproma,loidx(1),za(1),ua(1),dua(1))
        CALL lookup_ubc('cumastr',kproma,ztenh(1,jk),ub(1))
        zjk = REAL(jk,wp)
!IBM* NOVECTOR
        DO 421 jl=1,kproma
           ! mpuetz: move some of these into the previous loop
           zalvs=FSEL(tmelt - ztenh(jl,jk),als,alv)
           zalvdcp=zalvs*zcpcui(jl,jk)
           zqalv=1._wp/zalvs
           zhsat=zcpcu(jl,jk)*ztenh(jl,jk)+zgeoh(jl,jk)+zalvs*zqsenh(jl,jk)
           zpaphp1i = 1._wp/paphp1(jl,jk)
           zes=ua(jl)*zpaphp1i
           zes=MIN(0.5_wp,zes)
           zcor=1._wp/(1._wp-vtmpc1*zes)
           zqsat=zes*zcor
           zdqsdt=zpaphp1i*zcor**2*dua(jl)
           zgam=FSEL(zes - 0.4_wp,zqsat*zcor*ub(jl),zalvdcp*zdqsdt)
           zzz=zcpcu(jl,jk)*ztenh(jl,jk)*vtmpc1
           zhhat=zhsat-((zzz+zgam*zzz)/(1._wp+zgam*zzz*zqalv))*             &
                MAX(zqsenh(jl,jk)-zqenh(jl,jk),0._wp)
           zhhatt(jl,jk)=zhhat
           zlo1 = FSEL(zjk - zictop0(jl),0._wp,1._wp)
           zlo1 = FSEL(zhhat - zhcbase(jl),0._wp,zlo1)
           zictop0(jl) = FSEL(-zlo1,zictop0(jl),zjk)
421     END DO
     END IF
430 END DO
!
!!
!!     DEEP CONVECTION IF CLOUD DEPTH > 200 HPA, ELSE SHALLOW
!!     (CLOUD DEPTH FROM NON-ENTRAINIG PLUME)
!!
!  DO jl=1,kproma
!     ktype(jl)=MERGE(1,2,                                             &
!                paphp1(jl,kcbot(jl))-paphp1(jl,ictop0(jl)).gt.2.E4_wp)
!     zentr(jl)=MERGE(entrpen,entrscv,ktype(jl).eq.1)
!  ENDDO
!!
!*             (B) DO ASCENT IN 'CUASCT' IN ABSENCE OF DOWNDRAFTS
!                  ----------------------------------------------
!
  CALL cuasct(kproma, kbdim, klev, klevp1, klevm1,                     &
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
             kcbot,    kctop,    ictop0,   nmctop,                     &
!---Included for scavenging in wetdep_interface (Philip Stier, 19/02/04):-------
             zmwc,     zmrateprecip,       time_step_len               )
!---End Included for scavenging-----------------------------------------
!
!*        (C) CHECK CLOUD DEPTH AND CHANGE ENTRAINMENT RATE ACCORDINGLY
!             CALCULATE PRECIPITATION RATE (FOR DOWNDRAFT CALCULATION)
!             ---------------------------------------------------------
!
!DIR$ IVDEP
  DO 440 jl=1,kproma
     zpbmpt=paphp1(jl,kcbot(jl))-paphp1(jl,kctop(jl))
     IF(ldcum(jl).AND.ktype(jl).EQ.1.AND.zpbmpt.LT.2.e4_wp) ktype(jl)=2
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
  IF(lmfdd) THEN
!
!*             (A) DETERMINE LFS IN 'CUDLFS'
!                  -------------------------
!
     CALL cudlfs(kproma,   kbdim,    klev,     klevp1,                 &
                 ztenh,    zqenh,    puen,     pven,                   &
                 ktrac,                                                &
                 zxtenh,   pxtu,     zxtd,     zmfdxt,                 &
                 zgeoh,    paphp1,                                     &
                 ptu,      pqu,      zuu,      zvu,                    &
                 ldcum,    kcbot,    kctop,    zmfub,    zrfl,         &
                 ztd,      zqd,      zud,      zvd,                    &
                 pmfd,     zmfds,    zmfdq,    zdmfdp,                 &
                 zcpcu,                                                &
                 idtop,    loddraf)
!
!*            (B)  DETERMINE DOWNDRAFT T,Q AND FLUXES IN 'CUDDRAF'
!                 -----------------------------------------------
!
     CALL cuddraf(kproma,   kbdim,    klev,     klevp1,                &
                  ztenh,    zqenh,    puen,     pven,                  &
                  ktrac,                                               &
                  zxtenh,   zxtd,     zmfdxt,                          &
                  zgeoh,    paphp1,   zrfl,                            &
                  ztd,      zqd,      zud,      zvd,                   &
                  pmfd,     zmfds,    zmfdq,    zdmfdp,                &
                  zcpcu,                                               &
                  loddraf)
!
!*            (C)  RECALCULATE CONVECTIVE FLUXES DUE TO EFFECT OF
!*                 DOWNDRAFTS ON BOUNDARY LAYER MOISTURE BUDGET
!                  --------------------------------------------
!
!DIR$ IVDEP
  DO 520 jl=1,kproma
     IF(loddraf(jl)) THEN
        ikb=kcbot(jl)
        llo1=pmfd(jl,ikb).LT.0._wp
        zeps=MERGE(cmfdeps,0._wp,llo1)
        zqumqe=pqu(jl,ikb)+plu(jl,ikb)-                                &
                    zeps*zqd(jl,ikb)-(1._wp-zeps)*zqenh(jl,ikb)
        zdqmin=MAX(0.01_wp*zqenh(jl,ikb),1.e-10_wp)
        zmfmax=(paphp1(jl,ikb)-paphp1(jl,ikb-1))*zcons2
        llo1=zdqpbl(jl).GT.0._wp.AND.zqumqe.GT.zdqmin.AND.ldcum(jl)    &
                    .AND.zmfub(jl).LT.zmfmax
        zmfub1(jl)=MERGE(zdqpbl(jl)/(grav*MAX(zqumqe,zdqmin)),            &
                               zmfub(jl),llo1)
        zmfub1(jl)=MERGE(zmfub1(jl),zmfub(jl),                         &
                           (ktype(jl).EQ.1.OR.ktype(jl).EQ.2) .AND.    &
                         ABS(zmfub1(jl)-zmfub(jl)).LT.0.2_wp*zmfub(jl))
     END IF
520 END DO
  DO 540 jk=1,klev
     DO 530 jl=1,kproma
        IF(loddraf(jl)) THEN
           zfac=zmfub1(jl)/MAX(zmfub(jl),1.e-10_wp)
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
              zfac=zmfub1(jl)/MAX(zmfub(jl),1.e-10_wp)
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
  CALL cuasct(kproma,  kbdim,    klev,     klevp1,   klevm1,           &
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
             kcbot,    kctop,    ictop0,   nmctop,                     &
!---Included for scavenging in wetdep_interface (Philip Stier, 19/02/04):-------
             zmwc,     zmrateprecip,       time_step_len               )
!---End Included for scavenging-----------------------------------------
!
!---------------------------------------------------------------------
!*    7.0          DETERMINE FINAL CONVECTIVE FLUXES IN 'CUFLX'
!                  ------------------------------------------
!
  CALL cuflx(kproma,   kbdim,    klev,     klevp1,                     &
             pqen,     pqsen,    ztenh,    zqenh,                      &
             ktrac,                                                    &
!---Included for scavenging in wetdep_interface (Philip Stier, 28/03/01):-------
             krow,                                                     &
             pxtte,    pxtu,     ptu,                                  &
             zmwc,     zmrateprecip,                                   &
!---End Included for scavenging-----------------------------------------
             zxtenh,   zmfuxt,   zmfdxt,                               &
             paphp1,   zgeoh,                                          &
             kcbot,    kctop,    idtop,                                &
             ktype,    loddraf,  ldcum,                                &
             pmfu,     pmfd,     zmfus,    zmfds,                      &
             zmfuq,    zmfdq,    zmful,                                &
             zdmfup,   zdmfdp,   zrfl,     prain,                      &
             zcpcu,                                                    &
             pten,     zsfl,     zdpmel,   itopm2,   cevapcu,          &
             time_step_len                                             )
!---------------------------------------------------------------------
!
!*    8.0          UPDATE TENDENCIES FOR T AND Q IN SUBROUTINE CUDTDQ
!                  --------------------------------------------------
!
  CALL cudtdq(kproma, kbdim, klev, klevp1, itopm2, ldcum, ktrac,       &
!--- Included for dust emissions (Philip Stier 23/01/06)-----------------
              krow,                                                    &
!--- End Included for dust emissions in ---------------------------------
              paphp1,   pten,     ptte,     pqte,                      &
              pxtte,    pxtec,    zmfuxt,   zmfdxt,                    &
              zmfus,    zmfds,    zmfuq,    zmfdq,                     &
              zmful,    zdmfup,   zdmfdp,   plude,                     &
              zdpmel,   zrfl,     zsfl,                                &
              zcen,     zalvsh,   pqtec,    pqude,                     &
              prsfc,    pssfc,                                         &
              pch_concloud, pcon_dtrl, pcon_dtri,                      &
              pxtecl,   pxteci, pcon_iqte,                             &
              ptte_cnv, pqte_cnv, pxtte_cnv                            )
!
!---------------------------------------------------------------------
!
!*    9.0          UPDATE TENDENCIES FOR U AND U IN SUBROUTINE CUDUDV
!                  --------------------------------------------------
!
  IF(lmfdudv) THEN
     CALL cududv(kproma,   kbdim,    klev,     klevp1,                 &
                 itopm2,   ktype,    kcbot,    paphp1,   ldcum,        &
                 puen,     pven,     pvom,     pvol,                   &
                 zuu,      zud,      zvu,      zvd,                    &
                 pmfu,     pmfd,     pvom_cnv, pvol_cnv                )
!
  END IF
!
!  IF (lconvmassfix) THEN
!     CALL xt_conv_massfix(kproma,         kbdim,         klev,              &
!                          klevp1,         ktrac,         krow,              &
!                          papp1,          paphp1,        pxtte,             &
!                          .FALSE.  )
!  END IF
!
  RETURN
END SUBROUTINE cumastrt

SUBROUTINE cumastrh( kproma, kbdim, klev, klevp1, klevm1, ilab,        &
                     krow,    papp1,                                   &
           pten,     pqen,     pxen,     puen,     pven,               &
           ptven,    ktrac,    ldland,                                 &
           pxten,    pxtu,     pxtte,                                  &
           pverv,    pqsen,    pqhfla,                                 &
           paphp1,   pgeo,                                             &
           ptte,     pqte,     pvom,     pvol,                         &
           prsfc,    pssfc,    pxtec,                                  &
           pqtec,    pqude,    zcpen,    zcen,                         &
           pch_concloud, pcon_dtrl, pcon_dtri, pcon_iqte,              &
           ldcum,    ktype,    kcbot,    kctop,                        &
           ptu,      pqu,      plu,      plude,                        &
           pmfu,     pmfd,     prain,    pthvsig,                      &
           pxtecl,   pxteci,   nmctop,   cevapcu,  time_step_len,      &        
           ptte_cnv, pvom_cnv, pvol_cnv, pqte_cnv, pxtte_cnv     )
!
!**** *CUMASTRH*  MASTER ROUTINE FOR CUMULUS MASSFLUX-SCHEME
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
!          *CUMASTRH* IS CALLED FROM *CUCALL* IN CASE ICONV .EQ. 3
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
!       CUASCT:  CLOUD ASCENT FOR ENTRAINING PLUME
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
INTEGER, INTENT (IN) :: kproma, kbdim, klev, klevp1, ktrac, klevm1
REAL(wp),INTENT(IN)  :: cevapcu(:), time_step_len
INTEGER, INTENT(IN)  :: nmctop
REAL(wp),INTENT(INOUT):: pch_concloud(kbdim)
REAL(wp),INTENT(INOUT):: pcon_dtrl(kbdim), pcon_dtri(kbdim)
REAL(wp),INTENT(INOUT):: pcon_iqte(kbdim)
REAL(wp),INTENT(INOUT):: ptte_cnv(kbdim,klev)                               
REAL(wp),INTENT(INOUT):: pvom_cnv(kbdim,klev), pvol_cnv(kbdim,klev)         
REAL(wp),INTENT(INOUT):: pqte_cnv(kbdim), pxtte_cnv(kbdim,klev,ktrac)  
!---Included for in-cloud scavenging (Philip Stier, 19/01/06):----------
INTEGER, INTENT (IN) :: krow
!---End Included for scavenging-----------------------------------------
!
REAL(wp):: pten(kbdim,klev),        pqen(kbdim,klev),                  &
           pxen(kbdim,klev),        ptven(kbdim,klev),                 &
           puen(kbdim,klev),        pven(kbdim,klev),                  &
           ptte(kbdim,klev),        pqte(kbdim,klev),                  &
           pvom(kbdim,klev),        pvol(kbdim,klev),                  &
           pqsen(kbdim,klev),       pgeo(kbdim,klev),                  &
           paphp1(kbdim,klevp1),                                       &
           pverv(kbdim,klev)
REAL(wp):: ptu(kbdim,klev),         pqu(kbdim,klev),                   &
           plu(kbdim,klev),         plude(kbdim,klev),                 &
           pmfu(kbdim,klev),        pmfd(kbdim,klev),                  &
           prsfc(kbdim),            pssfc(kbdim),                      &
           prain(kbdim),            pqhfla(kbdim)
REAL(wp):: pthvsig(kbdim)
INTEGER :: kcbot(kbdim),            kctop(kbdim),                      &
           ktype(kbdim)
REAL(wp):: pxtec(kbdim,klev),       pqtec(kbdim,klev),                 &
           pxtecl(kbdim,klev),      pxteci(kbdim,klev),                &
           pqude(kbdim,klev)
REAL(wp):: ztenh(kbdim,klev),       zqenh(kbdim,klev),                 &
           zxenh(kbdim,klev),       zalvsh(kbdim,klev),                &
! zalvsh: latent heat of vaporisation/sublimation defined at half levels
           zgeoh(kbdim,klev),       zqsenh(kbdim,klev),                &
           ztd(kbdim,klev),         zqd(kbdim,klev),                   &
           zmfus(kbdim,klev),       zmfds(kbdim,klev),                 &
           zmfuq(kbdim,klev),       zmfdq(kbdim,klev),                 &
           zdmfup(kbdim,klev),      zdmfdp(kbdim,klev),                &
           zmful(kbdim,klev),       zrfl(kbdim),                       &
           zuu(kbdim,klev),         zvu(kbdim,klev),                   &
           zud(kbdim,klev),         zvd(kbdim,klev)
REAL(wp):: zcpen(kbdim,klev),       zcen(kbdim,klev),                  &
           zcpcu(kbdim,klev)
REAL(wp):: zentr(kbdim),            zhcbase(kbdim),                    &
           zmfub(kbdim),            zmfub1(kbdim),                     &
           zcpcui(kbdim,klev),      zictop0(kbdim),                    &
           zdqpbl(kbdim),           zdqcv(kbdim)
REAL(wp):: zsfl(kbdim),             zdpmel(kbdim,klev)
REAL(wp):: zcape(kbdim),            zheat(kbdim)
REAL(wp):: ua(kbdim), dua(kbdim), ub(kbdim), za(kbdim)
REAL(wp):: zhhatt(kbdim,klev)
INTEGER :: ilab(kbdim,klev),        idtop(kbdim),                      &
           ictop0(kbdim),           ilwmin(kbdim), loidx(kbdim)
REAL(wp):: pxten(kbdim,klev,ktrac), pxtte(kbdim,klev,ktrac),           &
           pxtu(kbdim,klev,ktrac),                                     &
           zxtenh(kbdim,klev,ktrac),zxtd(kbdim,klev,ktrac),            &
           zmfuxt(kbdim,klev,ktrac),zmfdxt(kbdim,klev,ktrac)
LOGICAL :: loddraf(kbdim),          ldland(kbdim)
LOGICAL :: ldcum(kbdim)
LOGICAL :: llo1
!---Included for scavenging in wetdep_interface (Philip Stier, 19/02/04):-------
REAL(wp):: papp1(kbdim,klev)
REAL(wp):: zmwc(kbdim,klev),        zmrateprecip(kbdim,klev)
!---End Included for scavenging-----------------------------------------
!
INTEGER :: jl, jk, ikb, jt, itopm2
REAL(wp):: zcons2, ztau, zqumqe, zdqmin, zmfmax, zalvs, zalvdcp, zqalv &
         , zhsat, zes, zcor, zqsat, zdqsdt, zgam, zzz, zhhat           &
         , zro, zdz, zfac, zpbmpt, zeps, zjk, zpaphp1i, zlo1
!
!  INTRINSIC FUNCTIONS
INTRINSIC MIN, MAX
!
!  Executable statements

!  IF(lconvmassfix) THEN
!     CALL xt_conv_massfix(kproma,         kbdim,            klev,    &
!                          klevp1,         ktrac,            krow,    &
!                          papp1,          paphp1,           pxtte,   &
!                          .TRUE.                                     )
!  END IF

!-----------------------------------------------------------------------
!
!     1.           SPECIFY CONSTANTS AND PARAMETERS
!                  --------------------------------
!
  zcons2=1._wp/(grav*time_step_len)
  ztau=cmftau
!  ztau=MIN(3._wp*3600._wp,7200._wp*63._wp/nn)
!
!-----------------------------------------------------------------------
!
!*    2.           INITIALIZE VALUES AT VERTICAL GRID POINTS IN 'CUINI'
!                  ---------------------------------------------------
!
  CALL cuini(kproma,   kbdim,    klev,     klevp1,   klevm1,           &
             pten,     pqen,     pqsen,    pxen,     puen,     pven,   &
             ptven,    ktrac,                                          &
             pxten,    zxtenh,   pxtu,     zxtd,     zmfuxt,   zmfdxt, &
             pverv,    pgeo,     paphp1,   zgeoh,                      &
             ztenh,    zqenh,    zqsenh,   zxenh,    ilwmin,           &
             ptu,      pqu,      ztd,      zqd,                        &
             zuu,      zvu,      zud,      zvd,                        &
             pmfu,     pmfd,     zmfus,    zmfds,                      &
             zmfuq,    zmfdq,    zdmfup,   zdmfdp,                     &
             zcpen,    zcpcu,    zalvsh,                               &
             zdpmel,   plu,      plude,    pqude,    ilab             )
!-----------------------------------------------------------------------
!
!*    3.0          CLOUD BASE CALCULATIONS
!                  -----------------------
!
!*             (A) DETERMINE CLOUD BASE VALUES IN 'CUBASE'
!                  ---------------------------------------
!
  CALL cubase(kproma,   kbdim,    klev,     klevp1,   klevm1,          &
              ztenh,    zqenh,    zgeoh,    paphp1,   pthvsig,         &
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
     zdqpbl(jl)=0.0_wp
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
     zdqmin=MAX(0.01_wp*zqenh(jl,ikb),1.e-10_wp)
     llo1=zdqpbl(jl).GT.0._wp.AND.zqumqe.GT.zdqmin.AND.ldcum(jl)
     zmfub(jl)=MERGE(zdqpbl(jl)/(grav*MAX(zqumqe,zdqmin)),0.01_wp,llo1)
     zmfmax=(paphp1(jl,ikb)-paphp1(jl,ikb-1))*zcons2
     zmfub(jl)=MIN(zmfub(jl),zmfmax)
     IF(.NOT.llo1) ldcum(jl)=.FALSE.
     ktype(jl)=MERGE(1,2,zdqcv(jl).GT.MAX(0._wp,-1.1_wp*pqhfla(jl)*grav))
     zentr(jl)=MERGE(entrpen,entrscv,ktype(jl).EQ.1)
340 END DO
!
!-----------------------------------------------------------------------
!*    4.0          DETERMINE CLOUD ASCENT FOR ENTRAINING PLUME
!                  -------------------------------------------
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
     zhcbase(jl)=zcpcu(jl,ikb)*ptu(jl,ikb)+zgeoh(jl,ikb)+              &
                                                      zalvs*pqu(jl,ikb)
     ictop0(jl)=kcbot(jl)-1
410 END DO
  DO 430 jk=klev,1,-1
     zcpcui(1:kproma,jk) = 1._wp/zcpcu(1:kproma,jk)
     IF (jk <= klevm1 .AND. jk >= 3) THEN
        ! mpuetz: too few instructions (FP dependencies)
        CALL prepare_ua_index_spline('cumastr',kproma,ztenh(1,jk),loidx(1),za(1))
        CALL lookup_ua_spline(kproma,loidx(1),za(1),ua(1),dua(1))
        CALL lookup_ubc('cumastr',kproma,ztenh(1,jk),ub(1))
        zjk = REAL(jk,wp)
!IBM* NOVECTOR
        DO 421 jl=1,kproma
           ! mpuetz: move some of these into the previous loop
           zalvs=FSEL(tmelt - ztenh(jl,jk),als,alv)
           zalvdcp=zalvs*zcpcui(jl,jk)
           zqalv=1._wp/zalvs
           zhsat=zcpcu(jl,jk)*ztenh(jl,jk)+zgeoh(jl,jk)+zalvs*zqsenh(jl,jk)
           zpaphp1i = 1._wp/paphp1(jl,jk)
           zes=ua(jl)*zpaphp1i
           zes=MIN(0.5_wp,zes)
           zcor=1._wp/(1._wp-vtmpc1*zes)
           zqsat=zes*zcor
           zdqsdt=zpaphp1i*zcor**2*dua(jl)
           zgam=FSEL(zes - 0.4_wp,zqsat*zcor*ub(jl),zalvdcp*zdqsdt)
           zzz=zcpcu(jl,jk)*ztenh(jl,jk)*vtmpc1
           zhhat=zhsat-((zzz+zgam*zzz)/(1._wp+zgam*zzz*zqalv))*             &
                MAX(zqsenh(jl,jk)-zqenh(jl,jk),0._wp)
           zhhatt(jl,jk)=zhhat
           zlo1 = FSEL(zjk - zictop0(jl),0._wp,1._wp)
           zlo1 = FSEL(zhhat - zhcbase(jl),0._wp,zlo1)
           zictop0(jl) = FSEL(-zlo1,zictop0(jl),zjk)
421     END DO
     END IF
430 END DO
!
!!
!!     DEEP CONVECTION IF CLOUD DEPTH > 200 HPA, ELSE SHALLOW
!!     (CLOUD DEPTH FROM NON-ENTRAINIG PLUME)
!!
!  DO jl=1,kproma
!     ktype(jl)=MERGE(1,2,                                               &
!                paphp1(jl,kcbot(jl))-paphp1(jl,ictop0(jl)).gt.2.E4_wp)
!     zentr(jl)=MERGE(entrpen,entrscv,ktype(jl).eq.1)
!  ENDDO
!!
!*             (B) DO ASCENT IN 'CUASCT' IN ABSENCE OF DOWNDRAFTS
!                  ----------------------------------------------
!
  CALL cuasct(kproma,  kbdim,    klev,     klevp1,   klevm1,           &
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
             kcbot,    kctop,    ictop0,   nmctop,                     &
!---Included for scavenging in wetdep_interface (Philip Stier, 19/02/04):-------
             zmwc,     zmrateprecip,       time_step_len               )
!---End Included for scavenging-----------------------------------------
!
!*         (C) CHECK CLOUD DEPTH AND CHANGE ENTRAINMENT RATE ACCORDINGLY
!              CALCULATE PRECIPITATION RATE (FOR DOWNDRAFT CALCULATION)
!              ---------------------------------------------------------
!
!DIR$ IVDEP
  DO 440 jl=1,kproma
     zpbmpt=paphp1(jl,kcbot(jl))-paphp1(jl,kctop(jl))
     IF(ldcum(jl).AND.ktype(jl).EQ.1.AND.zpbmpt.LT.2.e4_wp) ktype(jl)=2
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
!-----------------------------------------------------------------------
!*    5.0          CUMULUS DOWNDRAFT CALCULATIONS
!                  ------------------------------
!
  IF(lmfdd) THEN
!
!*             (A) DETERMINE LFS IN 'CUDLFS'
!                  -------------------------
!
     CALL cudlfs(kproma,   kbdim,    klev,     klevp1,                 &
                 ztenh,    zqenh,    puen,     pven,                   &
                 ktrac,                                                &
                 zxtenh,   pxtu,     zxtd,     zmfdxt,                 &
                 zgeoh,    paphp1,                                     &
                 ptu,      pqu,      zuu,      zvu,                    &
                 ldcum,    kcbot,    kctop,    zmfub,    zrfl,         &
                 ztd,      zqd,      zud,      zvd,                    &
                 pmfd,     zmfds,    zmfdq,    zdmfdp,                 &
                 zcpcu,                                                &
                 idtop,    loddraf)
!
!*            (B)  DETERMINE DOWNDRAFT T,Q AND FLUXES IN 'CUDDRAF'
!                  -----------------------------------------------
!
     CALL cuddraf(kproma,   kbdim,    klev,     klevp1,                &
                  ztenh,    zqenh,    puen,     pven,                  &
                  ktrac,                                               &
                  zxtenh,   zxtd,     zmfdxt,                          &
                  zgeoh,    paphp1,   zrfl,                            &
                  ztd,      zqd,      zud,      zvd,                   &
                  pmfd,     zmfds,    zmfdq,    zdmfdp,                &
                  zcpcu,                                               &
                  loddraf)
!
  END IF
!
!*    5.5          RECALCULATE CLOUD BASE MASSFLUX FROM A
!*                 CAPE CLOSURE FOR DEEP CONVECTION (KTYPE=1)
!*                 AND BY PBL EQUILIBRUM TAKING DOWNDRAFTS INTO
!*                 ACCOUNT FOR SHALLOW CONVECTION (KTYPE=2)
!                  -------------------------------------------
!
  DO jl=1,kproma
     zheat(jl)=0._wp
     zcape(jl)=0._wp
     zmfub1(jl)=zmfub(jl)
  ENDDO
!
  DO jk=1,klev
     DO jl=1,kproma
        llo1=ldcum(jl).AND.ktype(jl).EQ.1
        IF(llo1.AND.jk.LE.kcbot(jl).AND.jk.GT.kctop(jl)) THEN
           ikb=kcbot(jl)
           zro=paphp1(jl,jk)/(rd*ztenh(jl,jk)*                         &
                                           (1._wp+vtmpc1*zqenh(jl,jk)))
           zdz=(paphp1(jl,jk)-paphp1(jl,jk-1))/(grav*zro)
           zheat(jl)=zheat(jl) +                                       &
                (  (pten(jl,jk-1)-pten(jl,jk) + grav*zdz/zcpcu(jl,jk))    &
                     /ztenh(jl,jk)                                     &
                    +  vtmpc1*(pqen(jl,jk-1)-pqen(jl,jk))  ) *         &
                       (grav*(pmfu(jl,jk)+pmfd(jl,jk)))/zro
           zcape(jl)=zcape(jl) +                                       &
                         (grav*(ptu(jl,jk)-ztenh(jl,jk))/ztenh(jl,jk)     &
                              +grav*vtmpc1*(pqu(jl,jk)-zqenh(jl,jk))      &
                              -grav*plu(jl,jk) ) * zdz
        ENDIF
     ENDDO
  ENDDO
!

  DO jl=1,kproma
     IF(ldcum(jl).AND.ktype(jl).EQ.1) THEN
        ikb=kcbot(jl)
        zmfub1(jl)=(zcape(jl)*zmfub(jl))/(zheat(jl)*ztau)
        zmfub1(jl)=MAX(zmfub1(jl),0.001_wp)
        zmfmax=(paphp1(jl,ikb)-paphp1(jl,ikb-1))*zcons2
        zmfub1(jl)=MIN(zmfub1(jl),zmfmax)
     ENDIF
  ENDDO
!
!*                 RECALCULATE CONVECTIVE FLUXES DUE TO EFFECT OF
!*                 DOWNDRAFTS ON BOUNDARY LAYER MOISTURE BUDGET
!*                 FOR SHALLOW CONVECTION (KTYPE=2)
!                  --------------------------------------------
!
!DIR$ IVDEP
  DO 520 jl=1,kproma
     IF(ktype(jl).EQ.2) THEN
        ikb=kcbot(jl)
        llo1=pmfd(jl,ikb).LT.0._wp.AND.loddraf(jl)
        zeps=MERGE(cmfdeps,0._wp,llo1)
        zqumqe=pqu(jl,ikb)+plu(jl,ikb)-                                &
                    zeps*zqd(jl,ikb)-(1._wp-zeps)*zqenh(jl,ikb)
        zdqmin=MAX(0.01_wp*zqenh(jl,ikb),1.e-10_wp)
        zmfmax=(paphp1(jl,ikb)-paphp1(jl,ikb-1))*zcons2
        llo1=zdqpbl(jl).GT.0._wp.AND.zqumqe.GT.zdqmin.AND.ldcum(jl)    &
                    .AND.zmfub(jl).LT.zmfmax
        zmfub1(jl)=MERGE(zdqpbl(jl)/(grav*MAX(zqumqe,zdqmin)),            &
                               zmfub(jl),llo1)
        zmfub1(jl)=MERGE(zmfub1(jl),zmfub(jl),                         &
                         ABS(zmfub1(jl)-zmfub(jl)).LT.0.2_wp*zmfub(jl))
     END IF
520 END DO
  DO 540 jk=1,klev
     DO 530 jl=1,kproma
        IF(ldcum(jl)) THEN
           zfac=zmfub1(jl)/MAX(zmfub(jl),1.e-10_wp)
           pmfd(jl,jk)=pmfd(jl,jk)*zfac
           zmfds(jl,jk)=zmfds(jl,jk)*zfac
           zmfdq(jl,jk)=zmfdq(jl,jk)*zfac
           zdmfdp(jl,jk)=zdmfdp(jl,jk)*zfac
        END IF
530  END DO
!
     DO 5304 jt=1,ktrac
        DO 5302 jl=1,kproma
           IF(ldcum(jl)) THEN
              zfac=zmfub1(jl)/MAX(zmfub(jl),1.e-10_wp)
              zmfdxt(jl,jk,jt)=zmfdxt(jl,jk,jt)*zfac
           ENDIF
5302    END DO
5304 END DO
!
540 END DO
!
!*                 NEW VALUES OF CLOUD BASE MASS FLUX
!                  ----------------------------------
!
  DO 550 jl=1,kproma
     IF(ldcum(jl)) zmfub(jl)=zmfub1(jl)
550 END DO
!
!-----------------------------------------------------------------------
!*    6.0          DETERMINE FINAL CLOUD ASCENT FOR ENTRAINING PLUME
!*                 FOR PENETRATIVE CONVECTION (TYPE=1),
!*                 FOR SHALLOW TO MEDIUM CONVECTION (TYPE=2)
!*                 AND FOR MID-LEVEL CONVECTION (TYPE=3).
!                  --------------------------------------------------
!
  CALL cuasct(kproma,  kbdim,    klev,     klevp1,   klevm1,           &
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
             kcbot,    kctop,    ictop0,   nmctop,                     &
!---Included for scavenging in wetdep_interface (Philip Stier, 19/02/04):-------
             zmwc,     zmrateprecip,       time_step_len               )
!---End Included for scavenging-----------------------------------------
!
!-----------------------------------------------------------------------
!*    7.0          DETERMINE FINAL CONVECTIVE FLUXES IN 'CUFLX'
!                  ------------------------------------------
!
  CALL cuflx(kproma,   kbdim,    klev,     klevp1,                     &
             pqen,     pqsen,    ztenh,    zqenh,                      &
             ktrac,                                                    &
!---Included for scavenging in wetdep_interface (Philip Stier, 28/03/01):-------
             krow,                                                     &
             pxtte,    pxtu,     ptu,                                  &
             zmwc,     zmrateprecip,                                   &
!---End Included for scavenging-----------------------------------------
             zxtenh,   zmfuxt,   zmfdxt,                               &
             paphp1,   zgeoh,                                          &
             kcbot,    kctop,    idtop,                                &
             ktype,    loddraf,  ldcum,                                &
             pmfu,     pmfd,     zmfus,    zmfds,                      &
             zmfuq,    zmfdq,    zmful,                                &
             zdmfup,   zdmfdp,   zrfl,     prain,                      &
             zcpcu,                                                    &
             pten,     zsfl,     zdpmel,   itopm2,   cevapcu,          &
             time_step_len                                             )
!-----------------------------------------------------------------------
!
!*    8.0          UPDATE TENDENCIES FOR T AND Q IN SUBROUTINE CUDTDQ
!                  --------------------------------------------------
!
  CALL cudtdq(kproma, kbdim, klev, klevp1, itopm2, ldcum, ktrac,       &
!--- Included for dust emissions (Philip Stier 23/01/06)-----------------
              krow,                                                    &
!--- End Included for dust emissions in ---------------------------------
              paphp1,   pten,     ptte,     pqte,                      &
              pxtte,    pxtec,    zmfuxt,   zmfdxt,                    &
              zmfus,    zmfds,    zmfuq,    zmfdq,                     &
              zmful,    zdmfup,   zdmfdp,   plude,                     &
              zdpmel,   zrfl,     zsfl,                                &
              zcen,     zalvsh,   pqtec,    pqude,                     &
              prsfc,    pssfc,                                         &
              pch_concloud, pcon_dtrl, pcon_dtri,                      &
              pxtecl,   pxteci, pcon_iqte,                             &
              ptte_cnv, pqte_cnv, pxtte_cnv                            )
!
!-----------------------------------------------------------------------
!
!*    9.0          UPDATE TENDENCIES FOR U AND U IN SUBROUTINE CUDUDV
!                  --------------------------------------------------
!
  IF(lmfdudv) THEN
     CALL cududv(kproma,   kbdim,    klev,     klevp1,                 &
                 itopm2,   ktype,    kcbot,    paphp1,   ldcum,        &
                 puen,     pven,     pvom,     pvol,                   &
                 zuu,      zud,      zvu,      zvd,                    &
                 pmfu,     pmfd,     pvom_cnv, pvol_cnv                )
!
  END IF
!
!  IF (lconvmassfix) THEN
!     CALL xt_conv_massfix(kproma,         kbdim,         klev,              &
!                          klevp1,         ktrac,         krow,              &
!                          papp1,          paphp1,        pxtte,             &
!                          .FALSE.  )
!  END IF
!
  RETURN
END SUBROUTINE cumastrh
END MODULE mo_cumastr
