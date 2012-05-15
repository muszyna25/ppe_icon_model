!>
!! Computes evapotranspiration.
!!
!! @author Martin Koehler, DWD
!!
!! @par Revision History
!! Imported from IFS by Martin Koehler  (starting 2012-5-07)
!!   (IFS cycle CY36R1)
!!
!!-----------------------------------------------------------------------------
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
!!-----------------------------------------------------------------------------

MODULE mo_vevap

  PUBLIC :: vevap

CONTAINS

SUBROUTINE VEVAP(KIDIA,KFDIA,KLON,PTMST,PRVDIFTS,KTILE,&
 & PWLMX ,PTMLEV  ,PQMLEV  ,PAPHMS, PTSKM1M,PTSAM1M,&
 & PQS   ,PCFQ  ,PWETB  ,PWETL,PWETH,PWETHS,&
 & PCPTS ,PCSAT ,PCAIR ,PCSNW )  

! USE PARKIND1  ,ONLY : JPIM     ,JPRB
! USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
!
! USE YOS_CST  , ONLY : RG       ,RD       ,RCPD     ,RETV
! USE YOS_THF  , ONLY : RVTMP2
! USE YOS_VEG  , ONLY : RLHAERO, RLHAEROS

!ICON definitions:
USE mo_kind         ,ONLY : JPRB=>wp ,JPIM=>i4
USE mo_cuparameters ,ONLY : lhook    ,dr_hook  ,&           !yomcst  (& yos_exc)
      & RG       ,RD       ,RCPD     ,RETV     ,&           !yomcst  (& yos_cst)
      & RVTMP2                                              !yoethf  (& yos_thf)
USE mo_edmf_param   ,ONLY : &
      & RLHAERO, RLHAEROS                                   !yos_veg

!     ------------------------------------------------------------------

!**   *VEVAP* - COMPUTE EQUIVALENT EVAPOTRANSPIRATION EFFICIENCY

!     DERIVED FROM VDIFF (CY34) BY
!     A.C.M. BELJAARS       E.C.M.W.F.     18/01/90.

!     OBUKHOV-L UPDATE      ACMB           26/03/90.
!     (MAINLY TECHNICAL; TO MAKE CODE MORE READABLE)
!     Tiling of land surface ACMB          26/03/99.
!     Change surface units  P Viterbo      24/05/2004
!     Move to SURF library  P Viterbo      15/05/2005
!          (based on VDFEVAP)

!     PURPOSE
!     -------

!     COMPUTE EQUIVALENT EVAPOTRANSPIRATION EFFICIENCY

!     INTERFACE
!     ---------

!     *VEVAP* IS CALLED BY *SURFEXCDRIVER*

!     INPUT PARAMETERS (INTEGER):

!     *KIDIA*        START POINT
!     *KFDIA*        END POINT
!     *KLON*         NUMBER OF GRID POINTS PER PACKET
!     *KTILE*        TILE INDEX

!     INPUT PARAMETERS (REAL):

!     *PTMST*        TIME STEP
!     *PRVDIFTS*     Semi-implicit factor for vertical diffusion discretization
!     *PWLMX*        MAXIMUM INTERCEPTION LAYER CAPACITY
!     *PTMLEV*       TEMPERATURE AT T-1, lowest atmospheric level
!     *PQMLEV*       SPECIFIC HUMUDITY AT T-1, lowest atmospheric level
!     *PAPHMS*       PRESSURE AT T-1, surface
!     *PTSKM1M*      SKIN TEMPERATURE
!     *PTSAM1M*      SURFACE TEMPERATURE
!     *PQS*          SATURATION Q AT SURFACE
!     *PCFQ*         PROP. TO EXCH. COEFF. FOR MOISTURE(C-STAR IN DOC.)
!                    (SURFACE LAYER ON;Y)
!     *PWETB*        BARE SOIL RESISTANCE
!     *PWETL*        STOMATAL RESISTANCE LOW VEGETATION
!     *PWETH*        STOMATAL RESISTANCE HIGH VEGETATION, SNOW FREE
!     *PWETHS*       STOMATAL RESISTANCE HIGH VEGETATION WITH SNOW

!     OUTPUT PARAMETERS (REAL):

!     *PCPTS*        DRY STATIC ENRGY AT SURFACE
!     *PCSAT*        MULTIPLICATION FACTOR FOR QS AT SURFACE
!                    FOR SURFACE FLUX COMPUTATION
!     *PCAIR*        MULTIPLICATION FACTOR FOR Q AT LOWEST MODEL LEVEL
!                    FOR SURFACE FLUX COMPUTATION
!     *PCSNW*        MULTIPLICATION FACTOR FOR MOISTURE FLUX
!                    COMPUTATION FROM SNOW THROUGH CANOPY (TILE 7)

!     METHOD
!     ------

!     SEE DOCUMENTATION

!     ------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTMST 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRVDIFTS
INTEGER(KIND=JPIM),INTENT(IN)    :: KTILE 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PWLMX(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTMLEV(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQMLEV(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPHMS(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSKM1M(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSAM1M(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQS(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCFQ(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PWETB(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PWETL(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PWETH(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PWETHS(:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCPTS(:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCSAT(:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCAIR(:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCSNW(:) 
INTEGER(KIND=JPIM) :: JL

REAL(KIND=JPRB) :: ZCONS, ZCONS12, ZZWET, ZZWETS, ZRAS, ZCONS16, ZEP, ZEMAX
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!*    LOCAL STORAGE
!     ----- -------

!     ------------------------------------------------------------------

!*       1.     INITIALIZE CONSTANTS
!               ---------- ----------

! aerodynamic resistance for moisture transport from 
! top of high vegetation canopy to underlying snow (s/m)
! The value 1200. below is an estimate of rho*Cp at the surface

IF (LHOOK) CALL DR_HOOK('VEVAP_MOD:VEVAP',0,ZHOOK_HANDLE)
ZCONS12=PRVDIFTS*PTMST*RG/RD
ZCONS16=1./(RG*PTMST*PRVDIFTS)

!     ------------------------------------------------------------------

!          2.    COMPUTE EQUIVALENT EFFICIENCY FOR EVAPORATION
!                ------- ---------- ---------- --- -----------

PCSAT(KIDIA:KFDIA)=1.0_JPRB
PCAIR(KIDIA:KFDIA)=1.0_JPRB

!      interception reservoir

IF(KTILE == 3)THEN
  DO JL=KIDIA,KFDIA
    ZEP=ZCONS16*PCFQ(JL)*(PQMLEV(JL)-PQS(JL))
    ZEMAX=-PWLMX(JL)/PTMST
    IF (ZEP < 0.0_JPRB .AND. ZEMAX < 0.0_JPRB .AND. ZEP < ZEMAX) THEN
      PCAIR(JL)=MIN(1.0_JPRB,ZEMAX/ZEP)
      PCSAT(JL)=PCAIR(JL)
    ENDIF
  ENDDO
ENDIF

!      LOW VEGETATION

IF (KTILE  ==  4) THEN
  DO JL=KIDIA,KFDIA
    IF (PQMLEV(JL)  <=  PQS(JL)) THEN
      ZCONS=ZCONS12*PAPHMS(JL)/&
       & (PTMLEV(JL)*(1.0_JPRB+RETV*PQMLEV(JL)))  
      ZZWET=PWETL(JL)/ZCONS
      PCSAT(JL)=1.0_JPRB/(1.0_JPRB+PCFQ(JL)*ZZWET)
      PCAIR(JL)=1.0_JPRB/(1.0_JPRB+PCFQ(JL)*ZZWET)
    ENDIF
  ENDDO
ENDIF

!      HIGH VEGETATION (NO SNOW)

IF (KTILE  ==  6) THEN
  DO JL=KIDIA,KFDIA
    IF (PQMLEV(JL)  <=  PQS(JL)) THEN
      ZCONS=ZCONS12*PAPHMS(JL)/&
       & (PTMLEV(JL)*(1.0_JPRB+RETV*PQMLEV(JL)))  
      ZZWET=PWETH(JL)/ZCONS
      PCSAT(JL)=1.0_JPRB/(1.0_JPRB+PCFQ(JL)*ZZWET)
      PCAIR(JL)=PCSAT(JL)
    ENDIF
  ENDDO
ENDIF

!      HIGH VEGETATION (WITH UNDERLYING SNOW)

IF (KTILE == 7) THEN
  PCSNW(KIDIA:KFDIA)=0.0_JPRB
  DO JL=KIDIA,KFDIA
    IF (PQMLEV(JL) <= PQS(JL)) THEN
      IF(PTSKM1M(JL) > PTSAM1M(JL)) THEN
        ZRAS=1200._JPRB/RLHAEROS  
      ELSE
        ZRAS=1200._JPRB/RLHAERO  
      ENDIF
      ZCONS=ZCONS12*PAPHMS(JL)/&
       & (PTMLEV(JL)*(1.0_JPRB+RETV*PQMLEV(JL)))  
      ZZWET=PWETHS(JL)/ZCONS
      ZZWETS=ZRAS/ZCONS
      PCSAT(JL)=1.0_JPRB/(1.0_JPRB+ZZWET*PCFQ(JL)+ZZWET/ZZWETS)
      PCAIR(JL)=PCSAT(JL)
      IF (PWETHS(JL) > 1.0_JPRB) &
       & PCSNW(JL)=1.0_JPRB/(1.0_JPRB+ZZWETS*PCFQ(JL)+ZZWETS/ZZWET)  
    ENDIF
  ENDDO
ENDIF

!      BARE SOIL

IF (KTILE  ==  8) THEN
  DO JL=KIDIA,KFDIA
    IF (PQMLEV(JL)  <=  PQS(JL)) THEN
      ZCONS=ZCONS12*PAPHMS(JL)/&
       & (PTMLEV(JL)*(1.0_JPRB+RETV*PQMLEV(JL)))  
      ZZWET=PWETB(JL)/ZCONS
      PCSAT(JL)=1.0_JPRB/(1.0_JPRB+PCFQ(JL)*ZZWET)
      PCAIR(JL)=PCSAT(JL)
    ENDIF
  ENDDO
ENDIF

DO JL=KIDIA,KFDIA
  PCPTS(JL)=PTSKM1M(JL)*RCPD*(1.0_JPRB+RVTMP2*&
   & (PCSAT(JL)*PQS(JL)+(1.0_JPRB-PCAIR(JL))*PQMLEV(JL)))  
ENDDO

IF (LHOOK) CALL DR_HOOK('VEVAP_MOD:VEVAP',1,ZHOOK_HANDLE)
END SUBROUTINE VEVAP

END MODULE MO_VEVAP
