!>
!! Prepare surface boundary conditions for T and Q
!!
!! @author Martin Koehler, DWD
!!
!! @par Revision History
!! Imported from IFS by Martin Koehler  (starting 2012-4-30)
!!   (IFS cycle CY37R2)
!!
!!-----------------------------------------------------------------------------
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!-----------------------------------------------------------------------------

MODULE mo_vsurf
 
  PUBLIC :: vsurf

CONTAINS

SUBROUTINE VSURF(KIDIA,KFDIA,KLON,KLEVS,KTILE,&
 & KTVL,KTVH,&
 & PTMLEV, PQMLEV  ,PAPHMS,&
 & PTSKM1M,PWSAM1M,PTSAM1M,KSOTY,&
 & PSRFD ,PRAQ  ,PQSAM ,&
 & PQS   ,PDQS  ,&
 & PWETB ,PCPTS ,PWETL, PWETH, PWETHS)  

! USE PARKIND1  ,ONLY : JPIM     ,JPRB
! USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
! 
! USE YOS_VEG   ,ONLY : RCEPSW   ,RVROOTSA ,RVRSMIN  ,RVHSTR   ,RVLAI
! USE YOS_CST   ,ONLY : RCPD     ,RETV     ,RLVTT    ,RLSTT    ,RTT
! USE YOS_THF   ,ONLY : R2ES     ,R3LES    ,R3IES    ,R4LES    ,R4IES &
!                   &  ,R5LES    ,R5IES    ,RVTMP2
! USE YOS_SOIL  ,ONLY : RWCAP    ,RWPWP    ,RQWEVAP  &
!                   &  ,RWCAPM   ,RWPWPM   ,RQWEVAPM , RWRESTM &
!                   &  ,RTF1     ,RTF2     ,RTF3     ,RTF4 &
!                   &  ,LEVGEN
! USE YOS_EXC   ,ONLY : LEOCSA

!ICON definitions:
USE mo_kind         ,ONLY : JPRB=>wp ,JPIM=>i4
USE mo_cuparameters ,ONLY : lhook    ,dr_hook  ,&
      & RCPD     ,RETV     ,RVTMP2                            !yomcst (& yos_cst)
USE mo_edmf_param   ,ONLY : &
!dmk  & RWCAP    ,RWPWP    ,RQWEVAP  ,&                       ! (& yos_soil)
      & RWCAPM   ,RWPWPM   ,RQWEVAPM ,RWRESTM  ,&             ! (&  ---    )
      & RTF1     ,RTF2     ,RTF3     ,RTF4     ,&             ! (&  ---    )
!dmk  & LEVGEN   ,&                                           !yoephy (& yos_soil)
      & LEOCSA   ,&                                           !yoephy (& yos_exc)
      & RCEPSW   ,RVROOTSA ,RVRSMIN  ,RVHSTR   ,RVLAI  ,&     ! (& yos_veg)
      & FOEEW    ,FOEDESU

!mk: attention, LEVGEN=.TRUE. assumed to simplify code (new soil theory)
 

! #ifdef DOC
!     ------------------------------------------------------------------

!**   *VSURF* - PREPARES SURFACE BOUNDARY CONDITION FOR T AND Q

!     DERIVED FROM VDIFF (CY34) BY
!     A.C.M. BELJAARS       E.C.M.W.F.    18-1-90
!     Modified P.VITERBO AND A.C.M. BELJAARS  E.C.M.W.F.    16-3-93
!     Modified ACM Beljaars  26-03-99  Tiling of the surface
!     P. Viterbo     24-05-2004     Change surface units
!     P. Viterbo ECMWF 12/05/2005 Externalize SURF
!                     (based on VDFSURF)
!     G. Balsamo ECMWF 22/05/2006   Evaporative fraction f(soil)
!     G. Balsamo ECMWF 03/07/2006   Add soil type
!     E. Dutra/G. Balsamo 01/05/2008   Add lake tile

!     PURPOSE
!     -------

!     PREPARE SURFACE BOUNDARY CONDITION FOR Q AND T, E.G. FRACTIONAL
!     SURFACE COVER (SNOW AND VEGETATION), SATURATION SPECIFIC HUMIDITY
!     AT THE SURFACE, RELATIVE HUMIDITY OVER BARE LAND AND THE STOMATAL
!     RESISTANCE.

!     INTERFACE
!     ---------

!     *VSURF* IS CALLED BY *SURFEXCDRIVER*

!     INPUT PARAMETERS (INTEGER):

!     *KIDIA*        START POINT
!     *KFDIA*        END POINT
!     *KLON*         NUMBER OF GRID POINTS PER PACKET
!     *KLEVS*        NUMBER OF SOIL LAYERS
!     *KTILE*        TILE INDEX
!     *KTVL*         VEGETATION TYPE FOR LOW VEGETATION FRACTION
!     *KTVH*         VEGETATION TYPE FOR HIGH VEGETATION FRACTION
!     *KSOTY*        SOIL TYPE (1-7) 

!     INPUT PARAMETERS (REAL):

!     *PTMLEV*      TEMPERATURE AT T-1, lowest model level
!     *PQMLEV*      SPECIFIC HUMIDITY AT T-1, lowest model level
!     *PAPHMS*      PRESSURE AT T-1, surface
!     *PTSKM1M*      SURFACE TEMPERATURE
!     *PWSAM1M*      SOIL MOISTURE ALL LAYERS                   M**3/M**3
!     *PTSAM1M*      SOIL TEMPERATURE ALL LAYERS  
!     *PSRFD*        DOWNWARD SHORT WAVE RADIATION FLUX AT SURFACE
!     *PRAQ*         PRELIMINARY AERODYNAMIC RESISTANCE

!     OUTPUT PARAMETERS (REAL):

!     *PQSAM*        SPECIFIC HUMIDITY AT THE SURFACE
!     *PQS*          SATURATION Q AT SURFACE
!     *PDQS*         DERIVATIVE OF SATURATION Q-CURVE AT SURFACE T
!     *PWETB*        BARE SOIL RESISTANCE
!     *PCPTS*        DRY STATIC ENRGY AT SURFACE
!     *PWETL*        CANOPY RESISTANCE LOW VEGETATION
!     *PWETH*        CANOPY RESISTANCE HIGH VEGETATION, SNOW FREE
!     *PWETHS*       CANOPY RESISTANCE HIGH VEGETATION WITH SNOW

!     METHOD
!     ------

!     SEE DOCUMENTATION

!     ------------------------------------------------------------------
! #endif

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEVS 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KTILE 
INTEGER(KIND=JPIM),INTENT(IN)    :: KTVL(:) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KTVH(:) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KSOTY(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTMLEV(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQMLEV(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPHMS(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSKM1M(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PWSAM1M(:,:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSAM1M(:,:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSRFD(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRAQ(:) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PQSAM(:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PQS(:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDQS(:) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PWETB(:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCPTS(:) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PWETL(:) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PWETH(:) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PWETHS(:) 

!*    LOCAL STORAGE
!     ----- -------

REAL(KIND=JPRB) ::      ZLIQ(KLON,KLEVS)

INTEGER(KIND=JPIM) :: JK, JL, JS

REAL(KIND=JPRB) ::   &
 & ZCOR, ZEPSF3, ZF, ZF1H, ZF1L, ZF2H, ZF2L, ZF2B, &
 & ZF3H, ZF3L, ZHSTRH, ZHSTRL, ZLAIH, ZLAIL, &
 & ZQSAIR, ZROOT1H, ZROOT1L, ZROOT2H, ZROOT2L, &
 & ZROOT3H, ZROOT3L, ZROOT4H, ZROOT4L, ZRSMINH, &
 & ZRSMINL, ZRSMINB, ZSRFL, ZWROOTH, ZWROOTL, &
 & ZQWEVAP, ZWPWP,&
 & ZSALIN
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

! #include "fcsttre.h"

!     ------------------------------------------------------------------

!*       1.     INITIALIZE CONSTANTS
!               ---------- ----------

IF (LHOOK) CALL DR_HOOK('VSURF_MOD:VSURF',0,ZHOOK_HANDLE)

ZEPSF3=0.00001_JPRB ! security value for exponential sat-deficit dependence
ZRSMINB=50._JPRB  ! bare soil minimum resistance

!     ------------------------------------------------------------------

!          2.    PREPARE SURFACE BOUNDARY CONDITION
!                ------- ------- -------- ---------

!*         2.1   RELATIVE HUMIDITY OVER THE BARE LAND PART

!                BARE SOIL RESISTANCE IS COMPUTED FOR KTILE=4

!*         2.2   SATURATION PARAMETERS,

IF (LEOCSA .AND. KTILE  ==  1) THEN
  ZSALIN=0.98_JPRB
ELSE
  ZSALIN=1.0_JPRB
ENDIF

DO JL=KIDIA,KFDIA
  PQS(JL)=FOEEW(PTSKM1M(JL))/PAPHMS(JL)
  ZCOR=ZSALIN/(1.0_JPRB-RETV  *PQS(JL))
  PQS(JL)=PQS(JL)*ZCOR
  PDQS(JL)=PQS(JL)*ZCOR*FOEDESU(PTSKM1M(JL))
ENDDO

!*         2.3   DEFINITION OF THE STOMATAL RESISTANCE AND BARE SOIL RES
!*               DOES WORK FOR TYPE 4, 6 AND 8 WHEN ROUTINE IS CALLED FOR 
!*               TYPE 4

IF (KTILE  ==  4) THEN

!                Compute first liquid fraction of soil water to 
!                be used later in stress functions
!          CONTRIBUTION TO APPARENT ENERGY, TAKING INTO ACCOUNT
!          FREEZING/MELTING OF SOIL WATER.

  DO JK=1,KLEVS
    DO JL=KIDIA,KFDIA
      IF(PTSAM1M(JL,JK) < RTF1.AND.PTSAM1M(JL,JK) > RTF2) THEN
        ZF=0.5_JPRB*(1.0_JPRB-SIN(RTF4*(PTSAM1M(JL,JK)-RTF3)))
      ELSEIF (PTSAM1M(JL,JK) <= RTF2) THEN
        ZF=1.0_JPRB
      ELSE
        ZF=0.0_JPRB
      ENDIF
!dmk  IF (LEVGEN) THEN
        JS=KSOTY(JL)
        ZLIQ(JL,JK)=MAX(RWPWPM(JS),MIN(RWCAPM(JS),PWSAM1M(JL,JK)*(1.0_JPRB-ZF)))
!dmk  ELSE
!       ZLIQ(JL,JK)=MAX(RWPWP,MIN(RWCAP,PWSAM1M(JL,JK)*(1.0_JPRB-ZF)))
!xxx  ENDIF
    ENDDO
  ENDDO

  DO JL=KIDIA,KFDIA
!           minimal stomatal resistance : ZRSMIN
    ZRSMINL=RVRSMIN(KTVL(JL))
    ZRSMINH=RVRSMIN(KTVH(JL))

!           leaf area index  : ZLAI
    ZLAIL=RVLAI(KTVL(JL))
    ZLAIH=RVLAI(KTVH(JL))

!           soil moisture stress function : F2
    ZROOT1L=RVROOTSA(1,KTVL(JL))
    ZROOT2L=RVROOTSA(2,KTVL(JL))
    ZROOT3L=RVROOTSA(3,KTVL(JL))
    ZROOT4L=RVROOTSA(4,KTVL(JL))
    ZROOT1H=RVROOTSA(1,KTVH(JL))
    ZROOT2H=RVROOTSA(2,KTVH(JL))
    ZROOT3H=RVROOTSA(3,KTVH(JL))
    ZROOT4H=RVROOTSA(4,KTVH(JL))

    ZWROOTL=ZLIQ(JL,1)*ZROOT1L+&
     & ZLIQ(JL,2)*ZROOT2L+&
     & ZLIQ(JL,3)*ZROOT3L+&
     & ZLIQ(JL,4)*ZROOT4L  
    ZWROOTH=ZLIQ(JL,1)*ZROOT1H+&
     & ZLIQ(JL,2)*ZROOT2H+&
     & ZLIQ(JL,3)*ZROOT3H+&
     & ZLIQ(JL,4)*ZROOT4H  
!dmk IF (LEVGEN) THEN
       JS=KSOTY(JL)
       ZWPWP=RWPWPM(JS)
       ZQWEVAP=RQWEVAPM(JS)
!      bare ground evaporation scales with residual moisture in VG
!      ZQWEVAP=1._JPRB/(RWCAPM(JS)-RWRESTM(JS))
!      ZF2B=MAX(RCEPSW,MIN(1.0_JPRB,(ZLIQ(JL,1)-RWRESTM(JS))*ZQWEVAP))
!dmk ELSE
!      ZWPWP=RWPWP
!      ZQWEVAP=RQWEVAP
!xxx ENDIF
    ZF2B=MAX(RCEPSW,MIN(1.0_JPRB,(ZLIQ(JL,1)-ZWPWP)*ZQWEVAP))
    ZF2L=MAX(RCEPSW,MIN(1.0_JPRB,(ZWROOTL-ZWPWP)*ZQWEVAP))
    ZF2H=MAX(RCEPSW,MIN(1.0_JPRB,(ZWROOTH-ZWPWP)*ZQWEVAP))

!           radiation stress function (proposed by Alan Betts): ZF1 
    ZSRFL=PSRFD(JL)/250._JPRB
    ZF1L=1.0_JPRB/MAX(1.0_JPRB,0.81_JPRB*(1.+ZSRFL)/(ZSRFL+0.05_JPRB))
    ZF1H=ZF1L

!           atmospheric moisture deficit stress function : F3
    ZHSTRL=RVHSTR(KTVL(JL))
    ZHSTRH=RVHSTR(KTVH(JL))
    ZQSAIR=FOEEW(PTMLEV(JL))/PAPHMS(JL)
    ZCOR=1.0_JPRB/(1.0_JPRB-RETV  *ZQSAIR)
    ZQSAIR=ZQSAIR*ZCOR
    ZF3L=EXP(-ZHSTRL*(ZQSAIR-PQMLEV(JL)))
    ZF3H=EXP(-ZHSTRH*(ZQSAIR-PQMLEV(JL)))
    ZF3L=MAX(ZEPSF3,MIN(1.0_JPRB,ZF3L))
    ZF3H=MAX(ZEPSF3,MIN(1.0_JPRB,ZF3H))

    PWETL(JL)=ZRSMINL/ZLAIL/ZF1L/ZF2L/ZF3L
    PWETH(JL)=ZRSMINH/ZLAIH/ZF1H/ZF2H/ZF3H
    PWETHS(JL)=PWETH(JL)
    PWETB(JL)=ZRSMINB/ZF2B

  ENDDO
ENDIF

IF (KTILE == 4) THEN
  DO JL=KIDIA,KFDIA
    IF (PQMLEV(JL) > PQS(JL)) THEN
      PWETL(JL)=0.0_JPRB
    ENDIF
  ENDDO
ELSEIF (KTILE == 6) THEN
  DO JL=KIDIA,KFDIA
    IF (PQMLEV(JL) > PQS(JL)) THEN
      PWETH(JL)=0.0_JPRB
    ENDIF
  ENDDO
ELSEIF (KTILE == 7) THEN
  DO JL=KIDIA,KFDIA
    IF (PQMLEV(JL) > PQS(JL)) THEN
      PWETHS(JL)=0.0_JPRB
    ENDIF
  ENDDO
ELSEIF (KTILE == 8) THEN
  DO JL=KIDIA,KFDIA
    IF (PQMLEV(JL) > PQS(JL)) THEN
      PWETB(JL)=0.0_JPRB
    ENDIF
  ENDDO
ENDIF

!*         2.4   APPARENT SURFACE HUMIDITY

IF (KTILE  ==  1.OR. KTILE  ==  2.OR. KTILE  ==  3.OR. KTILE  ==  5 .OR. KTILE == 9 ) THEN   
  DO JL=KIDIA,KFDIA
    PQSAM(JL)=PQS(JL)
  ENDDO
ELSEIF (KTILE  ==  8) THEN
  DO JL=KIDIA,KFDIA
    PQSAM(JL)=PQS(JL)+(PQMLEV(JL)-PQS(JL))*PWETB(JL)/(PWETB(JL)+PRAQ(JL))
  ENDDO
ELSEIF (KTILE  ==  4) THEN
  DO JL=KIDIA,KFDIA
    PQSAM(JL)=PQS(JL)+(PQMLEV(JL)-PQS(JL))*PWETL(JL)/(PWETL(JL)+PRAQ(JL))
  ENDDO
ELSEIF (KTILE == 6) THEN ! I.E. HIGH VEGETATION, SNOW FREE
  DO JL=KIDIA,KFDIA
    PQSAM(JL)=PQS(JL)+(PQMLEV(JL)-PQS(JL))*PWETH(JL)/(PWETH(JL)+PRAQ(JL))
  ENDDO
ELSE ! I.E. HIGH VEGETATION WITH SNOW (7)
  DO JL=KIDIA,KFDIA
    PQSAM(JL)=PQS(JL)+(PQMLEV(JL)-PQS(JL))*PWETHS(JL)/(PWETHS(JL)+PRAQ(JL))
  ENDDO
ENDIF

!*         2.5   DRY STATIC ENERGY AT THE SURFACE

DO JL=KIDIA,KFDIA
  PCPTS(JL)=PTSKM1M(JL)*RCPD*(1.0_JPRB+RVTMP2*PQSAM(JL))
ENDDO

IF (LHOOK) CALL DR_HOOK('VSURF_MOD:VSURF',1,ZHOOK_HANDLE)
END SUBROUTINE VSURF


END MODULE MO_VSURF
