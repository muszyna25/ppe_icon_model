!>
!! Calculation of exchange coefficients between surface and lowest model level
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

MODULE mo_vexcs
 
  PUBLIC :: vexcs

CONTAINS

SUBROUTINE VEXCS(KIDIA,KFDIA,KLON,KITT,K_VMASS,LDINIT,PTMST,PRVDIFTS,&
 & PUMLEV,PVMLEV,PTMLEV,PQMLEV,PAPHMS,PGEOMLEV,PCPTGZLEV,&
 & PCPTS,PQSAM,PZ0MM,PZ0HM,PZ0QM,PZDL,PBUOM,PUCURR,PVCURR,&
 & PCFM,PCFH,PCFQ,PKM,PKH,PCM,PCH)  

! USE PARKIND1  ,ONLY : JPIM     ,JPRB
! USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
! 
! USE YOS_CST   , ONLY : RG       ,RD       ,RETV
! USE YOS_THF   , ONLY : RVTMP2
! USE YOS_EXC   , ONLY : RKAP     ,REPDU2   ,RPARZI
! USE YOS_EXCS  , ONLY : JPRITBL  ,RITBL    ,RCHBA    ,RCHBB    ,&
!  & RCHBD    ,RCHB23A  ,RCHBBCD  ,RCHBCD   ,RCHETA   ,&
!  & RCHETB   ,RCHBHDL  ,RCDHALF  ,RCDHPI2  ,DRITBL  

!ICON definitions:
USE mo_kind         ,ONLY : JPRB=>wp ,JPIM=>i4
USE mo_cuparameters ,ONLY : lhook    ,dr_hook  ,&           !yomcst  (& yos_exc)
      & RG       ,RD       ,RETV     ,&                     !yomcst  (& yos_cst)   
      & RVTMP2   ,&                                         !yoethf  (& yos_thf)
      & RKAP     ,REPDU2   ,RPARZI   ,&                     !yoevdf  (& yos_exc)
      & JPRITBL  ,RITBL    ,&                               !yoevdfs (& yos_excs)
      & RCHBHDL  ,DRITBL   ,&                     ! -
      & VDIV     ,VEXP     ,VREC     ,VLOG
USE mo_edmf_param   ,ONLY : &
      & PSIHU    ,PSIHS    ,PSIMU    ,PSIMS                 !fcsvdfs.h


! #ifdef DOC
!     ------------------------------------------------------------------

!**   *VEXCS* - DETERMINES THE EXCHANGE COEFFICIENTS BETWEEN THE
!                 SURFACE AND THE LOWEST MODEL LEVEL WITH HELP OF
!                 STABILITY AS FUNCTION OF OBUKHOV-L.

!     Original A.C.M. BELJAARS       E.C.M.W.F.    26/03/90.
!     Modified A.C.M. BELJAARS   26/03/99 Tiling of the land surface.
!     Modified J. HAGUE          13/01/03 MASS Vector Functions       
!     Modified   P. Viterbo ECMWF 12/05/2005 Externalize SURF
!                                   (based on vdfexcs)

!     PURPOSE
!     -------

!     DETERMINE EXCHANGE COEFFICIENTS BETWEEN THE SURFACE AND THE
!     LOWEST MODEL LEVEL

!     INTERFACE
!     ---------

!     *VEXCS* IS CALLED BY *SURFEXCDRIVER*

!     INPUT PARAMETERS (INTEGER):

!     *KIDIA*        START POINT
!     *KFDIA*        END POINT
!     *KLON*         NUMBER OF GRID POINTS PER PACKET
!     *KITT*         NUMBER OF ITERATIONS (BACK SUBST. IS ALWAYS DONE) [#]
!     *K_VMASS*      Controls the use of vector functions in the IBM scientific
!                     library. Set K_VMASS=0 to use standard functions

!     INPUT PARAMETER (LOGICAL):

!     *LDINIT*       IF .T. : ROUTINE GENERATES ITS OWN INITIAL GUESS  [#]
!                    IF .F. : ROUTINE STARTS WITH PZDL  AS INITIAL GUESS

!     INPUT PARAMETERS (REAL):

!     *PTMST*       TIME STEP
!     *PRVDIFTS*    Semi-implicit factor for vertical diffusion discretization
!     *PUMLEV*      X-VELOCITY COMPONENT AT T-1, lowest model level
!     *PVMLEV*      Y-VELOCITY COMPONENT AT T-1, lowest model level
!     *PTMLEV*      TEMPERATURE AT T-1, lowest model level
!     *PQMLEV*      SPECIFIC HUMUDITY AT T-1, lowest model level
!     *PAPHMS*      PRESSURE AT T-1, surface
!     *PGEOMLEV*    GEOPOTENTIAL AT T-1, lowest model level
!     *PCPTGZLEV*    DRY STATIC ENERGY, LOWEST MODEL LEVEL
!     *PCPTS*        DRY STATIC ENERGY AT THE SURFACE
!     *PQSAM*        SPECIFIC HUMIDITY AT THE SURFACE
!     *PZ0MM*        AERODYNAMIC ROUGHNESS LENGTH
!     *PZ0HM*        ROUGHNESS LENGTH FOR TEMPERATURE
!     *PZ0QM*        ROUGHNESS LENGTH FOR MOISTURE
!     *PZDL*         ZNLEV DEVIDED BY OBUKHOV LENGTH                 [#]
!     *PBUOM*        BUOYANCE FLUX AT THE SURFACE
!     *PUCURR*       OCEAN CURRENT X-COMPONENT
!     *PVCURR*       OCEAN CURRENT Y-COMPONENT

!     OUTPUT PARAMETERS (REAL):

!     *PCFM*         PROP. TO EXCH. COEFF. FOR MOMENTUM (C-STAR IN DOC.)
!                    AT LOWEST MODEL LEVEL; CALLED WITH ZCFM(1,KLEV)
!     *PCFH*         PROP. TO EXCH. COEFF. FOR HEAT     (C-STAR IN DOC.)
!                    AT LOWEST MODEL LEVEL; CALLED WITH ZCFM(1,KLEV)
!     *PCFQ*         PROP. TO EXCH. COEFF. FOR MOISTURE (C-STAR IN DOC.)
!                    (THIS ARRAY APPLIES TO BOTTOM LAYER ONLY)
!     *PKM*          CM*U  IN SURFACE LAYER                        (M/S)
!     *PKH*          CH*U  IN SURFACE LAYER                        (M/S)
!     *PCM*          CM    IN SURFACE LAYER                        (1)
!     *PCH*          CH    IN SURFACE LAYER                        (1)

!     REMARK: [#] UNUSED PARAMETERS IN TANGENT LINEAR AND ADJOINT VERSIONS
!     ------

!     METHOD
!     ------

!     THE ALGEBRAIC RELATION BETWEEN Z/L AND THE RICHARDSON NUMBER
!     IS SOLVED ITERATIVELY. THE STABILITY FUNCTIONS ARE THE SO-CALLED
!     PROFILE PSI-FUNCTIONS.
!     THE INITIAL GUESS (E.G. FROM PREVIOUS TIMESTEP) IS BACK-
!     SUBSTITUTED TO OBTAIN A SECOND APPROXIMATION. FURTHER ITERATION
!     IS DONE BY LINEAR INTER(EXTRA)POLATION (NEWTON'S  METHOD WITH
!     THE DERIVATIVE FROM SUCCESSIVE APPROXIMATIONS). IF NO INITIAL
!     GUESS IS PROVIDED, THE ROUTINE CAN PRODUCE ITS OWN. IN THE LATTER
!     CASE LDINIT=.T. HAS TO BE SPECIFIED AND IT IS RECOMMENDED TO
!     CHOOSE KITT=3. WITH INITIAL GUESSES FROM THE PREVIOUS TIME STEP
!     0 OR 1 ITERATION SHOULD BE SUFFICIENT.

!     ------------------------------------------------------------------
! #endif

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KITT 
INTEGER(KIND=JPIM),INTENT(IN)    :: K_VMASS
LOGICAL           ,INTENT(IN)    :: LDINIT 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTMST 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRVDIFTS
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUMLEV(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVMLEV(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTMLEV(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQMLEV(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPHMS(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGEOMLEV(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCPTGZLEV(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCPTS(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQSAM(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PZ0MM(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PZ0HM(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PZ0QM(:) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PZDL(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PBUOM(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUCURR(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVCURR(:)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PCFM(:) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PCFH(:) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PCFQ(:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PKM(:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PKH(:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCM(:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCH(:) 
!*            LOCAL STORAGE
!             ----- -------
INTEGER(KIND=JPIM) :: JLEN

REAL(KIND=JPRB) ::  Z1DZ0M(KLON+K_VMASS)  ,Z1DZ0H(KLON+K_VMASS),&
 & Z1DZ0Q(KLON+K_VMASS)  ,ZXLNM(KLON+K_VMASS),&
 & Z1DZ0MD(KLON+K_VMASS) ,Z1DZ0HD(KLON+K_VMASS) ,Z1DZ0QD(KLON+K_VMASS),&
 & ZXLNH(KLON+K_VMASS)   ,ZXLNQ(KLON+K_VMASS),&
 & ZETA2(KLON)           ,ZETA3(KLON) ,&
 & ZF2(KLON)             ,ZF3(KLON)             ,ZDU2(KLON+K_VMASS),&
 & ZNLEV(KLON),&
 & Z1DZ1D(KLON+K_VMASS)  ,ZRICLS(KLON+K_VMASS)  
REAL(KIND=JPRB) :: ZTMP1(KFDIA-KIDIA+1+K_VMASS)
REAL(KIND=JPRB) :: ZTMP2(KFDIA-KIDIA+1+K_VMASS)
REAL(KIND=JPRB) :: ZTMP3(KFDIA-KIDIA+1+K_VMASS)
REAL(KIND=JPRB) :: ZTMP4(KFDIA-KIDIA+1+K_VMASS)

INTEGER(KIND=JPIM) :: IRIB, JIT, JL

REAL(KIND=JPRB) :: ZA, ZAUX1, ZAUX2, ZB, &
 & ZCON1, ZCON2, ZCONS12, &
 & ZDRORO, ZETA, ZETA1, ZHU, ZIPBL, &
 & ZL, ZPRH, ZPRH0, ZPRH1, ZPRM, ZPRM0, ZPRM1, &
 & ZPRQ, ZPRQ0, ZRIB1, ZTPFAC1, ZUABS, ZWST2, &
 & ZX2  
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!             INCLUDE STABILITY FUNCTIONS
!             ------- --------- ---------

! #include "fcsvdfs.h" ! replaced by use statements

!     ------------------------------------------------------------------

!*       1.   INITIALIZE CONSTANTS
!             ---------- ----------

IF (LHOOK) CALL DR_HOOK('VEXCS_MOD:VEXCS',0,ZHOOK_HANDLE)
ZTPFAC1=PRVDIFTS

ZCONS12=ZTPFAC1*PTMST*RG/RD
ZCON1  =RVTMP2-RETV
ZCON2  =2.0_JPRB/3._JPRB

!     PBL HEIGHT FOR W* - EFFECT

ZIPBL=RPARZI

!     ------------------------------------------------------------------

!*       2.   COMPUTATION OF BASIC QUANTITIES
!             ----------- -- ----- ----------

IF(K_VMASS <= 0) THEN ! Not using Vector MASS

  DO JL=KIDIA,KFDIA

!             (W*)**2, WIND SHEAR, RICHARDSON NUMBER,

    IF (PBUOM(JL)  <=  0.0_JPRB) THEN
      ZWST2=0.0_JPRB
    ELSE
      ZWST2=(PBUOM(JL)*ZIPBL)**ZCON2
    ENDIF
    ZDU2(JL)=MAX(REPDU2,(PUMLEV(JL)-PUCURR(JL))**2&
     & +(PVMLEV(JL)-PVCURR(JL))**2+ZWST2)
    ZDRORO=2.0_JPRB*(PCPTGZLEV(JL)-PCPTS(JL))&
     & /(PCPTGZLEV(JL)+PCPTS(JL)-PGEOMLEV(JL))&
     & -ZCON1*(PQMLEV(JL)-MIN(PQSAM(JL),0.1_JPRB))
    ZRICLS(JL)=PGEOMLEV(JL)*ZDRORO/ZDU2(JL)

!             COMMON FACTORS IN NEUTRAL FORMULAE AND
!             DRAG COEFFICIENTS.

    ZNLEV (JL)=PGEOMLEV(JL)/RG+PZ0MM(JL)
    Z1DZ0M(JL)=ZNLEV(JL)/PZ0MM(JL)
    Z1DZ0H(JL)=ZNLEV(JL)/PZ0HM(JL)
    Z1DZ0Q(JL)=ZNLEV(JL)/PZ0QM(JL)
    Z1DZ1D(JL)=Z1DZ0M(JL)/(Z1DZ0M(JL)-1.0_JPRB)
    ZXLNM(JL) =LOG(Z1DZ0M(JL))
    ZXLNH(JL) =LOG(Z1DZ0H(JL))
    ZXLNQ(JL) =LOG(Z1DZ0Q(JL))
  ENDDO

ELSE ! Using Vector VMASS

  JLEN=KFDIA-KIDIA+K_VMASS-MOD(KFDIA-KIDIA,K_VMASS)
  DO JL=KIDIA,KFDIA
    ZTMP1(JL-KIDIA+1)=PCPTGZLEV(JL)+PCPTS(JL)-PGEOMLEV(JL)
    ZTMP3(JL-KIDIA+1)=ABS(PBUOM(JL)*ZIPBL)
  ENDDO 

  IF(KFDIA-KIDIA+1 /= JLEN) THEN
    ZTMP1(KFDIA-KIDIA+2:JLEN)=1.0_JPRB
    ZTMP3(KFDIA-KIDIA+2:JLEN)=1.0_JPRB
  ENDIF
  CALL VREC(ZTMP2,ZTMP1,JLEN)
  CALL VLOG(ZTMP4,ZTMP3,JLEN)
  DO JL=KIDIA,KFDIA
    ZTMP3(JL-KIDIA+1)=ZTMP4(JL-KIDIA+1)*ZCON2
  ENDDO
  CALL VEXP(ZTMP4,ZTMP3,JLEN)

  DO JL=KIDIA,KFDIA
!             (W*)**2, WIND SHEAR, RICHARDSON NUMBER,
    IF (PBUOM(JL)  <=  0.0_JPRB) THEN
      ZWST2=0.0_JPRB
    ELSE
      ZWST2=ZTMP4(JL-KIDIA+1)
    ENDIF
    ZDU2(JL)=MAX(REPDU2,(PUMLEV(JL)-PUCURR(JL))**2&
     & +(PVMLEV(JL)-PVCURR(JL))**2+ZWST2)
    ZDRORO=2.0_JPRB*(PCPTGZLEV(JL)-PCPTS(JL))*ZTMP2(JL-KIDIA+1)&
     &- ZCON1*(PQMLEV(JL)-MIN(PQSAM(JL),0.1_JPRB))
    ZTMP1(JL-KIDIA+1)=PGEOMLEV(JL)*ZDRORO

!             COMMON FACTORS IN NEUTRAL FORMULAE AND
!             DRAG COEFFICIENTS.

    ZNLEV (JL)=PGEOMLEV(JL)*(1.0_JPRB/RG)+PZ0MM(JL)
    Z1DZ0M(JL)=ZNLEV(JL)/PZ0MM(JL)
    Z1DZ0H(JL)=ZNLEV(JL)/PZ0HM(JL)
    Z1DZ0Q(JL)=ZNLEV(JL)/PZ0QM(JL)
    ZTMP2(JL-KIDIA+1)=Z1DZ0M(JL)-1.0_JPRB
  ENDDO

  IF(KFDIA-KIDIA+1 /= JLEN) THEN
    Z1DZ0M(KFDIA+1:JLEN)=1.0_JPRB
    Z1DZ0H(KFDIA+1:JLEN)=1.0_JPRB
    Z1DZ0Q(KFDIA+1:JLEN)=1.0_JPRB
    ZDU2  (KFDIA+1:JLEN)=1.0_JPRB
    ZDU2  (KFDIA+1:JLEN)=1.0_JPRB
    ZTMP2  (KFDIA-KIDIA+2:JLEN)=1.0_JPRB
  ENDIF
  CALL VDIV(ZRICLS(KIDIA),ZTMP1,ZDU2(KIDIA),JLEN)
  CALL VDIV(Z1DZ1D(KIDIA),Z1DZ0M(KIDIA),ZTMP2,JLEN)
  CALL VLOG(ZXLNM(KIDIA),Z1DZ0M(KIDIA),JLEN)  
  CALL VLOG(ZXLNH(KIDIA),Z1DZ0H(KIDIA),JLEN)
  CALL VLOG(ZXLNQ(KIDIA),Z1DZ0Q(KIDIA),JLEN)
ENDIF

!*       3. EXCHANGE COEFFIENTS DEPENDING UPON MONIN-OBUKHOV LENGTH
!           -------- ---------- --------- ---- ------------- ------

!           SOLVE Z/L ITERATIVELY
!           ----- --- -----------

IF(K_VMASS > 0) THEN
  CALL VREC(Z1DZ0MD(KIDIA),Z1DZ0M(KIDIA),JLEN)
  CALL VREC(Z1DZ0HD(KIDIA),Z1DZ0H(KIDIA),JLEN)
ELSE 
  DO JL=KIDIA,KFDIA
    Z1DZ0MD(JL)=1.0_JPRB/Z1DZ0M(JL)
    Z1DZ0HD(JL)=1.0_JPRB/Z1DZ0H(JL)
  ENDDO
ENDIF

!        3.1  FIRST AND SECOND GUESS TO INITIALIZE
!             ----- --- ------ ----- -- ----------

IF (LDINIT) THEN
  DO JL=KIDIA,KFDIA

!             FIRST GUESS: INTERPLOLATE ETA-TABLE FOR POS. RI
!             -----------

    IF (ZRICLS(JL)  >  0.0_JPRB) THEN
      ZRIB1=ZRICLS(JL)*ZXLNM(JL)/ZXLNH(JL)
      IRIB=INT(ZRIB1/DRITBL)+1
      IF (IRIB  >=  JPRITBL) THEN
        ZETA1=RITBL(JPRITBL)
      ELSE
        ZX2  = IRIB*DRITBL
        ZA   = (ZX2-ZRIB1)/DRITBL
        ZB   = 1.0_JPRB-ZA
        ZETA1= ZA*RITBL(IRIB) + ZB*RITBL(IRIB+1)
      ENDIF
    ELSE
      ZETA1=ZRICLS(JL)*ZXLNM(JL)/ZXLNH(JL)
    ENDIF

!             SECOND GUESS: WITH LOGARITHMIC TERMS
!             ------ -----

    ZETA2(JL)=ZETA1*ZXLNM(JL)*Z1DZ1D(JL)

  ENDDO
ELSE
  DO JL=KIDIA,KFDIA
    ZETA2(JL)=PZDL(JL)
  ENDDO
ENDIF

!        3.2  BACK SUBSTITUTION
!             ---- ------------

DO JL=KIDIA,KFDIA
  IF (ZETA2(JL)  >  0.0_JPRB) THEN
    ZPRH1=PSIHS(ZETA2(JL))
    ZPRH0=PSIHS(ZETA2(JL)*Z1DZ0HD(JL))
    ZPRM1=PSIMS(ZETA2(JL))
    ZPRM0=PSIMS(ZETA2(JL)*Z1DZ0MD(JL))
  ELSE
    ZPRH1=PSIHU(ZETA2(JL))
    ZPRH0=PSIHU(ZETA2(JL)*Z1DZ0HD(JL))
    ZPRM1=PSIMU(ZETA2(JL))
    ZPRM0=PSIMU(ZETA2(JL)*Z1DZ0MD(JL))
  ENDIF
  ZPRH=ZXLNH(JL)-ZPRH1+ZPRH0
  ZPRM=ZXLNM(JL)-ZPRM1+ZPRM0
  ZHU =ZPRM**2*Z1DZ1D(JL)/ZPRH
  ZETA3(JL) =ZRICLS(JL)*ZHU
  ZF2(JL)  =ZRICLS(JL)-ZETA2(JL)/ZHU
ENDDO

!        3.3  ITERATE KITT TIMES
!             ------- ---- -----

DO JIT=1, KITT
  DO JL=KIDIA,KFDIA
    IF (ZETA3(JL)  >  0.0_JPRB) THEN
      ZPRH1=PSIHS(ZETA3(JL))
      ZPRH0=PSIHS(ZETA3(JL)*Z1DZ0HD(JL))
      ZPRM1=PSIMS(ZETA3(JL))
      ZPRM0=PSIMS(ZETA3(JL)*Z1DZ0MD(JL))
    ELSE
      ZPRH1=PSIHU(ZETA3(JL))
      ZPRH0=PSIHU(ZETA3(JL)*Z1DZ0HD(JL))
      ZPRM1=PSIMU(ZETA3(JL))
      ZPRM0=PSIMU(ZETA3(JL)*Z1DZ0MD(JL))
    ENDIF
    ZPRH   =ZXLNH(JL)-ZPRH1+ZPRH0
    ZPRM   =ZXLNM(JL)-ZPRM1+ZPRM0
    ZF3(JL)=ZRICLS(JL)-ZETA3(JL)*ZPRH/(Z1DZ1D(JL)*ZPRM**2)
    IF (ABS(ZF3(JL)-ZF2(JL))  >  1.E-25_JPRB) THEN
      ZETA=ZETA2(JL)-ZF2(JL)*(ZETA3(JL)-ZETA2(JL))/(ZF3(JL)-ZF2(JL))
    ELSE
      ZETA=ZETA3(JL)
    ENDIF

!             SAFETY PROVISION FOR DIVERGING ITERATIONS

    IF (ZETA*ZETA3(JL)  <  0.0_JPRB) ZETA=ZETA3(JL)*0.5_JPRB
    ZETA2(JL)=ZETA3(JL)
    ZF2(JL)  =ZF3(JL)
    ZETA3(JL)=ZETA
  ENDDO
ENDDO

!        3.4  COPY ZNLEV/L TO OUTPUT ARRAY
!             ---- ------- -- ------ -----

DO JL=KIDIA,KFDIA
  PZDL(JL)=ZETA3(JL)
ENDDO

!        4.   COMPUTE EXCHANGE COEFFICIENTS
!             ------- -------- ------------

DO JL=KIDIA,KFDIA

!             ASSUME MAXIMUM NOCTURNAL BOUNDARY LAYER HEIGHT
!             EQUAL TO RCHBHDL*L (5 TIMES OBUKHOV LENGTH)
!             AND UNIFORM PROFILES FROM H TO ZNLEV

  IF (ZETA3(JL)  >  RCHBHDL) THEN
    ZL=ZNLEV(JL)/ZETA3(JL)
    ZNLEV(JL)=RCHBHDL*ZL+MAX(PZ0MM(JL),PZ0QM(JL))
    ZETA3(JL)=ZNLEV(JL)/ZL
    Z1DZ0M(JL)=ZNLEV(JL)/PZ0MM(JL)
    Z1DZ0H(JL)=ZNLEV(JL)/PZ0HM(JL)
    Z1DZ0Q(JL)=ZNLEV(JL)/PZ0QM(JL)
    ZXLNM(JL) =LOG(Z1DZ0M(JL))
    ZXLNH(JL) =LOG(Z1DZ0H(JL))
    ZXLNQ(JL) =LOG(Z1DZ0Q(JL))
  ENDIF
ENDDO

IF(K_VMASS > 0) THEN
  CALL VREC(Z1DZ0MD(KIDIA),Z1DZ0M(KIDIA),JLEN)
  CALL VREC(Z1DZ0HD(KIDIA),Z1DZ0H(KIDIA),JLEN)
  CALL VREC(Z1DZ0QD(KIDIA),Z1DZ0Q(KIDIA),JLEN)
ELSE
  DO JL=KIDIA,KFDIA
    Z1DZ0MD(JL)=1.0_JPRB/Z1DZ0M(JL)
    Z1DZ0HD(JL)=1.0_JPRB/Z1DZ0H(JL)
    Z1DZ0QD(JL)=1.0_JPRB/Z1DZ0Q(JL)
  ENDDO
ENDIF

DO JL=KIDIA,KFDIA
  IF (ZETA3(JL)  >  0.0_JPRB) THEN
    ZPRM1=PSIMS(ZETA3(JL))
    ZPRM0=PSIMS(ZETA3(JL)*Z1DZ0MD(JL))
    ZPRH1=PSIHS(ZETA3(JL))
    ZPRH0=PSIHS(ZETA3(JL)*Z1DZ0HD(JL))
    ZPRQ0=PSIHS(ZETA3(JL)*Z1DZ0QD(JL))
  ELSE
    ZPRM1=PSIMU(ZETA3(JL))
    ZPRM0=PSIMU(ZETA3(JL)*Z1DZ0MD(JL))
    ZPRH1=PSIHU(ZETA3(JL))
    ZPRH0=PSIHU(ZETA3(JL)*Z1DZ0HD(JL))
    ZPRQ0=PSIHU(ZETA3(JL)*Z1DZ0QD(JL))
  ENDIF
  ZPRM   =ZXLNM(JL)-ZPRM1+ZPRM0
  ZPRH   =ZXLNH(JL)-ZPRH1+ZPRH0
  ZPRQ   =ZXLNQ(JL)-ZPRH1+ZPRQ0
  ZAUX1  =ZCONS12*PAPHMS(JL)/&
   & (PTMLEV(JL)*(1.0_JPRB+RETV*PQMLEV(JL)))  ! g * rho * alpha * dt
  ZUABS  =SQRT(ZDU2(JL))                      ! |U|
  ZAUX2  =ZAUX1*ZUABS*RKAP**2                 ! |U| * g * rho * alpha * dt * kappa^2

  PCFM(JL)=ZAUX2/(ZPRM**2)                    ! C* = C * |U| * g * rho * alpha * dt
  PCFH(JL)=ZAUX2/(ZPRM*ZPRH)
  PCFQ(JL)=ZAUX2/(ZPRM*ZPRQ)

!             U*CH FOR POSTPROCESSING ONLY  [m/s]

  PKM(JL)=PCFM(JL)/ZAUX1
  PKH(JL)=PCFH(JL)/ZAUX1

!             CM and CH FOR TERRA           [1]

  PCM(JL)=PCFM(JL)/ZAUX1/ZUABS
  PCH(JL)=PCFH(JL)/ZAUX1/ZUABS


ENDDO

IF (LHOOK) CALL DR_HOOK('VEXCS_MOD:VEXCS',1,ZHOOK_HANDLE)
END SUBROUTINE VEXCS


END MODULE MO_VEXCS
