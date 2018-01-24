!>
!! COMPUTES THE SURFACE (2 M) TEMPERATURE AND HUMIDITY
!! WITH STABILITY AS FUNCTION OF OBUKHOV-L
!!
!! @author Martin Koehler, DWD
!!
!! @par Revision History
!! Imported from IFS by Martin Koehler  (starting 2012-4-30)
!!   (IFS cycle CY36R1)
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

MODULE mo_sppcfl
 
  PUBLIC :: sppcfl

CONTAINS

SUBROUTINE SPPCFL(KIDIA, KFDIA, KLON &
 & , PUMLEV, PVMLEV, PQMLEV, PGEOMLEV, PCPTS, PCPTGZLEV &
 & , PAPHMS, PZ0MM, PZ0HM, PZ0QM, PZDL, PQSA &
 & , PBLEND, PFBLEND, PUCURR, PVCURR &
 ! OUTPUTS
 & , PU10, PV10, PT2, PD2, PQ2  )  

! USE PARKIND1  ,ONLY : JPIM     ,JPRB
! USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
! 
! USE YOS_CST   , ONLY : RG       ,RD       ,RCPD     ,RETV     ,&
!  & RLVTT    ,RLSTT    ,RTT  
! USE YOS_THF   , ONLY : R2ES     ,R3LES    ,R4LES    ,&
!  & RVTMP2
! USE YOS_EXCS  , ONLY : RCHBA    ,RCHBB    ,RCHBD    ,RCHB23A  ,&
!  & RCHBBCD  ,RCHBCD   ,RCHETA   ,RCHETB   ,RCHBHDL  ,&
!  & RCDHALF  ,RCDHPI2

!ICON definitions:
USE mo_kind         ,ONLY : JPRB=>wp ,JPIM=>i4
USE mo_cuparameters ,ONLY : lhook    ,dr_hook  ,&           !yomcst  (& yos_exc)
      & RG       ,RD       ,RCPD     ,RETV     ,&           !yomcst  (& yos_cst)
      & RTT      ,&                                         ! -
      & R2ES     ,R3LES    ,R4LES    ,RVTMP2   ,&           !yoethf (& yos_thf)
      & RCHBHDL                                             ! -
USE mo_edmf_param   ,ONLY : &
      & PSIHU    ,PSIHS    ,PSIMU    ,PSIMS                 !fcsvdfs.h

!     ------------------------------------------------------------------

!**   *SPPCFL* - COMPUTES THE SURFACE (2 M) TEMPERATURE AND HUMIDITY
!                WITH STABILITY AS FUNCTION OF OBUKHOV-L

!     P. VITERBO         E.C.M.W.F.    08/10/93. (MODIFIED FROM VDFT2M)

!     Modified   P. Viterbo ECMWF 12/05/2005 Externalize SURF
!                                   (based on vdfppcfl)

!     PURPOSE
!     -------

!     COMPUTE WIND OCMPONENTS, TEMPERATURE AND DEWPOINT TEMPERATURE
!     AT SCREEN LEVEL HEIGHT

!     INTERFACE
!     ---------

!     *SPPCFL* IS CALLED BY *SURFPP*

!     INPUT PARAMETERS (INTEGER):

!     *KIDIA*        START POINT
!     *KFDIA*        END POINT
!     *KLON*         NUMBER OF GRID POINTS PER PACKET

!     INPUT PARAMETERS (REAL):

!     *PUMLEV*       U-COMPONENT WIND AT T-1, lowest atmospheric level
!     *PVMLEV*       V-COMPONENT WIND AT T-1, lowest atmospheric level
!     *PQMLEV*       SPECIFIC HUMIDITY AT T-1, lowest atmospheric level
!     *PAPHMS*       surface pressure
!     *PGEOMLEV*     GEOPOTENTIAL AT T-1, lowest atmospheric level
!     *PCPTS*        DRY STATIC ENRGY AT SURFACE
!     *PCPTGZLEV*    DRY STATIC ENERGY    AT T-1, lowest atmospheric level
!     *PZ0MM*        AERODYNAMIC ROUGHNESS LENGTH
!     *PZ0HM*        ROUGHNESS LENGTH FOR TEMPERATURE
!     *PZ0QM*        ROUGHNESS LENGTH FOR MOISTURE
!     *PZDL*         ZNLEV+Z0M DEVIDED BY OBUKHOV LENGTH
!     *PQSA*         APPARENT SURFACE HUMIDITY
!     *PBLEND*       HEIGHT FROM WHICH WIND SPEED IS INTERPOLATED TO 10 M
!     *PFBLEND*      WIDN SPEED AT HEIGHT PBLEND
!     *PUCURR*       U COMPONENT OF OCEAN SURFACE CURRENT
!     *PVCURR*       V COMPONENT OF OCEAN SURFACE CURRENT

!     OUTPUT PARAMETERS (REAL):

!     *PU10*         U-COMPONENT WIND AT 10 M
!     *PV10*         V-COMPONENT WIND AT 10 M
!     *PT2*          TEMPERATURE AT 2 M
!     *PD2*          DEW POINT TEMPERATURE AT 2 M
!     *PQ2*          SPECIFIC HUMIDITY AT 2 M

!     METHOD
!     ------

!     UNIVERSAL FUNCTIONS ARE USED TO INTERPOLATE BETWEEN SURF. AND NLEV

!     ------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUMLEV(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVMLEV(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQMLEV(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGEOMLEV(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCPTS(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCPTGZLEV(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPHMS(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PZ0MM(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PZ0HM(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PZ0QM(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PZDL(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQSA(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PBLEND(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PFBLEND(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUCURR(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVCURR(:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PU10(:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PV10(:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PT2(:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PD2(:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PQ2(:) 

!*    LOCAL STORAGE
!     ----- -------

INTEGER(KIND=JPIM) :: JL

REAL(KIND=JPRB) :: Z10DNL, Z10M, Z10MP, Z2M, Z2MP,&
 & ZAPH2M, ZCPT2M, ZCVM3, ZCVM4, ZDL, &
 & ZF1, ZFRAC, ZHTQ, ZL, ZNLEV, ZPRH0, ZPRH1, &
 & ZPRH2, ZPRM0, ZPRM10, ZPRM2, ZPRQ0, &
 & ZWIND, ZZQM1  
REAL(KIND=JPRB) :: ZHOOK_HANDLE

! #include "fcsvdfs.h" ! replaced by use statements

!     ------------------------------------------------------------------

!*       1.   INITIALIZE CONSTANTS
!             ---------- ----------

!     THIS SPECIFIES THE HEIGHT FOR U,V (10 M) AND T,Q (2 M)

IF (LHOOK) CALL DR_HOOK('SPPCFL_MOD:SPPCFL',0,ZHOOK_HANDLE)
Z10M=10._JPRB
Z2M=2.0_JPRB
ZHTQ=Z2M*RG

!     ------------------------------------------------------------------

!        2.   COMPUTE  WIND COMPONENTS AND TEMPERATURE
!             ----------------------------------------

DO JL=KIDIA,KFDIA
  ZNLEV=PGEOMLEV(JL)/RG+PZ0MM(JL)
  ZDL=PZDL(JL)
  ZL=ZNLEV/PZDL(JL)

  IF (PZDL(JL)  >  RCHBHDL) THEN
    ZNLEV=RCHBHDL*ZL+MAX(PZ0MM(JL),PZ0HM(JL))
    ZDL=ZNLEV/ZL
  ENDIF

  Z2MP =Z2M+PZ0MM(JL)
  IF (Z2MP/ZL  >  RCHBHDL) THEN
    Z2MP=RCHBHDL*ZL+MAX(PZ0MM(JL),PZ0HM(JL))
  ENDIF

  Z10MP =Z10M+PZ0MM(JL)
  IF (Z10MP/ZL  >  RCHBHDL) THEN
    Z10MP=RCHBHDL*ZL+MAX(PZ0MM(JL),PZ0HM(JL))
  ENDIF

  ZWIND =PBLEND(JL)+PZ0MM(JL)
  IF (ZWIND/ZL  >  RCHBHDL) THEN
    ZWIND=RCHBHDL*ZL+MAX(PZ0MM(JL),PZ0HM(JL))
  ENDIF

  IF (PZDL(JL)  >  0.0_JPRB) THEN
    ZPRH2=PSIHS(ZDL*Z2MP/ZNLEV)
    ZPRH1=PSIHS(ZDL)
    ZPRH0=PSIHS(ZDL*PZ0HM(JL)/ZNLEV)

    ZPRQ0=PSIHS(ZDL*PZ0QM(JL)/ZNLEV)

    ZPRM2 =PSIMS(ZDL*ZWIND/ZNLEV)
    ZPRM10=PSIMS(ZDL*Z10MP/ZNLEV)
    ZPRM0 =PSIMS(ZDL*PZ0MM(JL)/ZNLEV)
  ELSE
    ZPRH2=PSIHU(ZDL*Z2MP/ZNLEV)
    ZPRH1=PSIHU(ZDL)
    ZPRH0=PSIHU(ZDL*PZ0HM(JL)/ZNLEV)

    ZPRQ0=PSIHU(ZDL*PZ0QM(JL)/ZNLEV)

    ZPRM2 =PSIMU(ZDL*ZWIND/ZNLEV)
    ZPRM10=PSIMU(ZDL*Z10MP/ZNLEV)
    ZPRM0 =PSIMU(ZDL*PZ0MM(JL)/ZNLEV)
  ENDIF

!     10 M WIND COMPONENT

  Z10DNL= (LOG(Z10MP/PZ0MM(JL))-ZPRM10+ZPRM0)&
   & /(LOG(ZWIND/PZ0MM(JL))-ZPRM2+ZPRM0)  
  ZF1=MAX(0.1_JPRB,SQRT((PUMLEV(JL)-PUCURR(JL))**2&
   & +(PVMLEV(JL)-PVCURR(JL))**2))
  PU10(JL)=(PUMLEV(JL)-PUCURR(JL))*Z10DNL*PFBLEND(JL)/ZF1+PUCURR(JL)
  PV10(JL)=(PVMLEV(JL)-PVCURR(JL))*Z10DNL*PFBLEND(JL)/ZF1+PVCURR(JL)

!     2M TEMPERATURE AND 2M SPECIFIC HUMIDITY

  ZZQM1=MAX(1.E-12_JPRB,PQMLEV(JL))
  ZCPT2M=PCPTS(JL)+ (PCPTGZLEV(JL)-PCPTS(JL))&
   & *(LOG(Z2MP /PZ0HM(JL))-ZPRH2+ZPRH0)&
   & /(LOG(ZNLEV/PZ0HM(JL))-ZPRH1+ZPRH0)  
  PQ2(JL)=PQSA(JL) + (ZZQM1-PQSA(JL))&
   & *(LOG(Z2MP /PZ0QM(JL))-ZPRH2+ZPRQ0)&
   & /(LOG(ZNLEV/PZ0QM(JL))-ZPRH1+ZPRQ0)  

!     APPROXIMATE MOISTURE CORRECTION IN CP WITH QNLEV

  PT2(JL)=(ZCPT2M-RG*Z2M)/( RCPD*(1.0_JPRB+RVTMP2*ZZQM1) )
ENDDO

!     ------------------------------------------------------------------

!        3.   COMPUTE  DEW POINT TEMPERATURE
!             ------------------------------

DO JL=KIDIA,KFDIA
  ZZQM1=MAX(1.E-12_JPRB,PQMLEV(JL))
  ZAPH2M=PAPHMS(JL)*(1.0_JPRB-ZHTQ/(RD*PT2(JL)*(1.0_JPRB+RETV*ZZQM1)))
! Note the use ofs saturation with respect to WATER only, in conformity to
!   WMO observation reporting standards
  ZCVM3=R3LES
  ZCVM4=R4LES
  ZFRAC=LOG(ZAPH2M*PQ2(JL)/(R2ES*(1.0_JPRB+RETV*PQ2(JL))))/ZCVM3
  PD2(JL)=(RTT-ZFRAC*ZCVM4)/(1.0_JPRB-ZFRAC)

!     LIMIT DEW POINT TEMPERATURE < TEMPERATURE

  PD2(JL)=MIN(PT2(JL),PD2(JL))
ENDDO

IF (LHOOK) CALL DR_HOOK('SPPCFL_MOD:SPPCFL',1,ZHOOK_HANDLE)
END SUBROUTINE SPPCFL


END MODULE MO_SPPCFL
