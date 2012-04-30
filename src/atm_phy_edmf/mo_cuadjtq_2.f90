!>
!! saturation adjustment for EDMF DUALM
!!
!! @author Martin Koehler, DWD
!!
!! @par Revision History
!! Imported from IFS by Martin Koehler  (starting 2011-9-30)
!!   (IFS cycle CY36R1_DUALM_M8b)
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

MODULE mo_cuadjtq_2
 
  PUBLIC :: cuadjtq2

CONTAINS

!! #ifdef RS6K
!! @PROCESS HOT NOSTRICT
!! #endif
SUBROUTINE CUADJTQ2 &
 & (KIDIA,    KFDIA,    KLON,    KLEV,     KK,&
 &  PSP,      PT,       PQ,      LDFLAG,   KCALL)  

!          M.TIEDTKE         E.C.M.W.F.     12/89

!          MODIFICATIONS
!          -------------
!          D.SALMOND         CRAY(UK))      12/8/91
!          J.J. MORCRETTE    ECMWF          92-09-18   Update to Cy44
!          J.F. MAHFOUF      ECMWF          96-06-11   Smoothing option
!          D.SALMOND & M.HAMRUD ECMWF       99-06-04   Optimisation
!          J.HAGUE                          03-01-13   MASS Vector Functions
!          J.HAGUE                          03-07-07   More MASS V.F.
!        M.Hamrud              01-Oct-2003 CY28 Cleaning
!        J.Hague & D.Salmond   22-Nov-2005 Optimisations 

!          PURPOSE.
!          --------
!          TO PRODUCE T,Q AND L VALUES FOR CLOUD ASCENT

!          INTERFACE
!          ---------
!          THIS ROUTINE IS CALLED FROM SUBROUTINES:
!              *COND*     (T AND Q AT CONDENSATION LEVEL)
!              *CUBASE*   (T AND Q AT CONDENSATION LEVEL)
!              *CUASC*    (T AND Q AT CLOUD LEVELS)
!              *CUINI*    (ENVIRONMENTAL T AND QS VALUES AT HALF LEVELS)
!              *CUSTRAT*  (T AND Q AT CONDENSATION LEVEL)
!          INPUT ARE UNADJUSTED T AND Q VALUES,
!          IT RETURNS ADJUSTED VALUES OF T AND Q

!     PARAMETER     DESCRIPTION                                   UNITS
!     ---------     -----------                                   -----
!     INPUT PARAMETERS (INTEGER):

!    *KIDIA*        START POINT
!    *KFDIA*        END POINT
!    *KLON*         NUMBER OF GRID POINTS PER PACKET
!    *KLEV*         NUMBER OF LEVELS
!    *KK*           LEVEL
!    *KCALL*        DEFINES CALCULATION AS
!                      KCALL=0  ENV. T AND QS IN*CUINI*
!                      KCALL=1  CONDENSATION IN UPDRAFTS  (E.G. CUBASE, CUASC)
!                      KCALL=2  EVAPORATION IN DOWNDRAFTS (E.G. CUDLFS,CUDDRAF)

!     INPUT PARAMETERS (LOGICAL):

!    *LDLAND*       LAND-SEA MASK (.TRUE. FOR LAND POINTS)

!     INPUT PARAMETERS (REAL):

!    *PSP*          PRESSURE                                        PA

!     UPDATED PARAMETERS (REAL):

!    *PT*           TEMPERATURE                                     K
!    *PQ*           SPECIFIC HUMIDITY                             KG/KG

!          EXTERNALS   
!          ---------
!          3 LOOKUP TABLES ( TLUCUA, TLUCUB, TLUCUC )
!          FOR CONDENSATION CALCULATIONS.
!          THE TABLES ARE INITIALISED IN *SUPHEC*.

!----------------------------------------------------------------------

! USE PARKIND1  ,ONLY : JPIM     ,JPRB
! USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
! USE YOMCST   , ONLY : RETV     ,RLVTT    ,RLSTT    ,RTT
! USE YOETHF   , ONLY : R2ES     ,R3LES    ,R3IES    ,R4LES    ,&
!  & R4IES    ,R5LES    ,R5IES    ,R5ALVCP  ,R5ALSCP  ,&
!  & RALVDCP  ,RALSDCP  ,RTWAT    ,RTICE    ,RTICECU  ,&
!  & RTWAT_RTICE_R      ,RTWAT_RTICECU_R  
! USE YOEPHLI  , ONLY : LPHYLIN  ,RLPTRC   ,RLPAL1   ,RLPAL2
! USE YOMJFH   , ONLY : N_VMASS

!ICON definitions:
USE mo_kind         ,ONLY : JPRB=>wp ,JPIM=>i4
USE mo_cuparameters ,ONLY : lhook    ,dr_hook  ,&
                & RETV     ,RLVTT    ,RLSTT    ,RTT      ,&           !yomcst
                & R2ES     ,R3LES    ,R3IES    ,R4LES    ,&           !yoethf
                & R4IES    ,R5LES    ,R5IES    ,R5ALVCP  ,R5ALSCP  ,& ! -
                & RALVDCP  ,RALSDCP  ,RTWAT    ,RTICE    ,RTICECU  ,& ! -
                & RTWAT_RTICE_R      ,RTWAT_RTICECU_R    ,&           ! -
                & LPHYLIN  ,RLPTRC   ,RLPAL1   ,RLPAL2   ,&           !yoephli
                & vdiv     ,vexp     ,vrec
USE mo_edmf_param   ,ONLY : &
                & N_VMASS  ,&                                         !yomjfh
                & FOEALFCU ,FOEEWMCU ,FOEDEMCU ,FOELDCPMCU         ,& !fcttre.h
                & FOEALFA  ,FOEEWM   ,FOEDEM   ,FOELDCPM              ! -

IMPLICIT NONE


INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KK 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSP(KLON) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PT(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PQ(KLON,KLEV) 
LOGICAL           ,INTENT(IN)    :: LDFLAG(KLON) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KCALL 
INTEGER(KIND=JPIM) :: JL, JLEN, JJ, JJJ, JINDEX(KLON), IKNDEX(KLON), IKK

REAL(KIND=JPRB) :: ZTMP0(KFDIA-KIDIA+1)
REAL(KIND=JPRB) :: ZTMP1(2*(KFDIA-KIDIA+1))
REAL(KIND=JPRB) :: ZTMP3(2*(KFDIA-KIDIA+1))
! REAL(KIND=JPRB) :: ZTMP4(KFDIA-KIDIA+1)
REAL(KIND=JPRB) :: ZTMP5(2*(KFDIA-KIDIA+1))
! REAL(KIND=JPRB) :: ZTMP6(KFDIA-KIDIA+1)
REAL(KIND=JPRB) :: ZTEMP1(2*KLON)
REAL(KIND=JPRB) :: ZTEMP2(2*KLON)
REAL(KIND=JPRB) :: ZTEMP3(2*KLON)
REAL(KIND=JPRB) :: ZTEMP4(2*KLON)
REAL(KIND=JPRB) :: ZTEMP5(2*KLON)
REAL(KIND=JPRB) :: ZTEMP6(2*KLON)

REAL(KIND=JPRB) :: Z1S, Z2S, ZCOND,ZCOND1, ZCOR, ZFOEEWI, ZFOEEWL,&
 & ZOEALFA, ZQMAX, ZQSAT, ZTARG, ZQP
REAL(KIND=JPRB) :: ZL, ZI, ZF
REAL(KIND=JPRB) :: ZPT, ZPQ

!DIR$ VFUNCTION EXPHF
! #include "fcttre.h" ! replaced by use statements

!     STATEMENT FUNCTIONS
!REAL_B :: FOEALFAJ,FOEDEMJ,FOELDCPMJ,FOEEWMJ

REAL(KIND=JPRB) :: MINJ, MAXJ, X, Y
REAL(KIND=JPRB) :: ZHOOK_HANDLE

MINJ(X,Y) = Y - 0.5_JPRB*(ABS(X-Y)-(X-Y))
MAXJ(X,Y) = Y + 0.5_JPRB*(ABS(X-Y)+(X-Y))

!----------------------------------------------------------------------

!     1.           DEFINE CONSTANTS
!                  ----------------

IF (LHOOK) CALL DR_HOOK('CUADJTQ',0,ZHOOK_HANDLE)

IF(N_VMASS >  0) THEN
  JLEN=KFDIA-KIDIA+1
ENDIF

ZQMAX=0.5_JPRB

!*********************************************
IF (.NOT.LPHYLIN) THEN
!*********************************************                 

!     2.           CALCULATE CONDENSATION AND ADJUST T AND Q ACCORDINGLY
!                  -----------------------------------------------------


  IF (KCALL == 1 ) THEN

!DIR$    IVDEP
!OCL NOVREC
    JJJ=0
    DO JL=KIDIA,KFDIA
      IF(LDFLAG(JL)) THEN
        JJJ = JJJ + 1
        JINDEX(JJJ) = JL
      ENDIF
    ENDDO

    IF(N_VMASS <= 0)  THEN  ! Don't use Vector MASS

      DO JL=KIDIA,KFDIA
        IF(LDFLAG(JL)) THEN
          ZQP    =1.0_JPRB/PSP(JL)
!         ZQSAT=FOEEWMCU(PT(JL,KK))*ZQP    
          ! FOEEWMCU ( PTARE ) = R2ES *&
          !  &(FOEALFCU(PTARE)*EXP(R3LES*(PTARE-RTT)/(PTARE-R4LES))+&
          !  &(1.0_JPRB-FOEALFCU(PTARE))*EXP(R3IES*(PTARE-RTT)/(PTARE-R4IES)))
          ZL=1.0_JPRB/(PT(JL,KK)-R4LES)
          ZI=1.0_JPRB/(PT(JL,KK)-R4IES)
          ZQSAT=R2ES *(FOEALFCU(PT(JL,KK))*EXP(R3LES*(PT(JL,KK)-RTT)*ZL)+&
            &(1.0_JPRB-FOEALFCU(PT(JL,KK)))*EXP(R3IES*(PT(JL,KK)-RTT)*ZI))
          ZQSAT=ZQSAT*ZQP
          ZQSAT=MIN(0.5_JPRB,ZQSAT)
          ZCOR=1.0_JPRB-RETV*ZQSAT
!         ZCOND=(PQ(JL,KK)*ZCOR**2-ZQSAT*ZCOR)/(ZCOR**2+ZQSAT*FOEDEMCU(PT(JL,KK)))
          ! FOEDEMCU ( PTARE )=FOEALFCU(PTARE)*R5ALVCP*(1.0_JPRB/(PTARE-R4LES)**2)+&
          !   &(1.0_JPRB-FOEALFCU(PTARE))*R5ALSCP*(1.0_JPRB/(PTARE-R4IES)**2)
          ZF=FOEALFCU(PT(JL,KK))*R5ALVCP*ZL**2 + &
            &(1.0_JPRB-FOEALFCU(PT(JL,KK)))*R5ALSCP*ZI**2
          ZCOND=(PQ(JL,KK)*ZCOR**2-ZQSAT*ZCOR)/(ZCOR**2+ZQSAT*ZF)
!         ZCOND=MAX(ZCOND,0.0_JPRB)
          IF(ZCOND > 0.0_JPRB)THEN
            PT(JL,KK)=PT(JL,KK)+FOELDCPMCU(PT(JL,KK))*ZCOND
            PQ(JL,KK)=PQ(JL,KK)-ZCOND
!           ZQSAT=FOEEWMCU(PT(JL,KK))*ZQP    
            ZL=1.0_JPRB/(PT(JL,KK)-R4LES)
            ZI=1.0_JPRB/(PT(JL,KK)-R4IES)
            ZQSAT=R2ES *(FOEALFCU(PT(JL,KK))*EXP(R3LES*(PT(JL,KK)-RTT)*ZL)+&
              &(1.0_JPRB-FOEALFCU(PT(JL,KK)))*EXP(R3IES*(PT(JL,KK)-RTT)*ZI))
            ZQSAT=ZQSAT*ZQP
            ZQSAT=MINJ(0.5_JPRB,ZQSAT)
            ZCOR=1.0_JPRB-RETV*ZQSAT
!           ZCOND1=(PQ(JL,KK)*ZCOR**2-ZQSAT*ZCOR)/(ZCOR**2+ZQSAT*FOEDEMCU(PT(JL,KK)))
            ZF=FOEALFCU(PT(JL,KK))*R5ALVCP*ZL**2 + &
              &(1.0_JPRB-FOEALFCU(PT(JL,KK)))*R5ALSCP*ZI**2
            ZCOND1=(PQ(JL,KK)*ZCOR**2-ZQSAT*ZCOR)/(ZCOR**2+ZQSAT*ZF)
            IF(ZCOND ==  0.0_JPRB)ZCOND1=0.0_JPRB
            PT(JL,KK)=PT(JL,KK)+FOELDCPMCU(PT(JL,KK))*ZCOND1
            PQ(JL,KK)=PQ(JL,KK)-ZCOND1
          ENDIF
        ENDIF
      ENDDO

    ELSE ! Use Vector Mass

      DO JJ=1,JJJ
        JL = JINDEX(JJ)
        ZL=1.0_JPRB/(PT(JL,KK)-R4LES)
        ZI=1.0_JPRB/(PT(JL,KK)-R4IES)
        ZTEMP1(2*JJ-1) = (R3LES*(PT(JL,KK)-RTT)*ZL)
        ZTEMP1(2*JJ  ) = (R3IES*(PT(JL,KK)-RTT)*ZI)
        ZTEMP2(2*JJ-1) = ZI
        ZTEMP2(2*JJ  ) = ZL
      ENDDO
      DO JJ=1,2*JJJ
        ZTEMP1 (JJ)= EXP(ZTEMP1(JJ))
      ENDDO
      IKK=0
      DO JJ=1,JJJ
        JL = JINDEX(JJ)
        ZI= ZTEMP2(2*JJ-1)
        ZL= ZTEMP2(2*JJ  )
        ZQP    =1.0_JPRB/PSP(JL)
        ZQSAT=R2ES *(FOEALFCU(PT(JL,KK))*(ZTEMP1(2*JJ-1))+&
          &(1.0_JPRB-FOEALFCU(PT(JL,KK)))*(ZTEMP1(2*JJ  )))
        ZQSAT=ZQSAT*ZQP
        ZQSAT=MIN(0.5_JPRB,ZQSAT)
        ZCOR=1.0_JPRB-RETV*ZQSAT
        ZF=FOEALFCU(PT(JL,KK))*R5ALVCP*ZL**2 + &
          &(1.0_JPRB-FOEALFCU(PT(JL,KK)))*R5ALSCP*ZI**2
        ZCOND=(PQ(JL,KK)*ZCOR**2-ZQSAT*ZCOR)/(ZCOR**2+ZQSAT*ZF)
        ZPT=PT(JL,KK)+FOELDCPMCU(PT(JL,KK))*ZCOND
        ZPQ=PQ(JL,KK)-ZCOND
        ZL=1.0_JPRB/(ZPT-R4LES)
        ZI=1.0_JPRB/(ZPT-R4IES)
        IF(ZCOND > 0.0_JPRB)THEN
          IKK=IKK+1
          IKNDEX(IKK)=JL
          ZTEMP3(2*IKK-1) = R3LES*(ZPT-RTT)*ZL
          ZTEMP3(2*IKK  ) = R3IES*(ZPT-RTT)*ZI
          ZTEMP4(2*IKK-1) = ZI
          ZTEMP4(2*IKK  ) = ZL
          ZTEMP5(2*IKK-1) = ZPT
          ZTEMP5(2*IKK  ) = ZPQ
          ZTEMP6(2*IKK  ) = ZQP
        ENDIF
      ENDDO
      DO JJ=1,2*IKK
        ZTEMP3(JJ) = EXP(ZTEMP3(JJ))
      ENDDO
      DO JJ=1,IKK
        JL = IKNDEX(JJ)
        ZI= ZTEMP4(2*JJ-1)
        ZL= ZTEMP4(2*JJ  )
        ZPT = ZTEMP5(2*JJ-1)
        ZPQ = ZTEMP5(2*JJ  )
        ZQP = ZTEMP6(2*JJ  )
          ZQSAT=R2ES *(FOEALFCU(ZPT)*(ZTEMP3(2*JJ-1)  )+&
            &(1.0_JPRB-FOEALFCU(ZPT))*(ZTEMP3(2*JJ)  ))
          ZQSAT=ZQSAT*ZQP
          ZQSAT=MIN(0.5_JPRB,ZQSAT)
          ZCOR=1.0_JPRB-RETV*ZQSAT
          ZF=FOEALFCU(ZPT)*R5ALVCP*ZL**2 + &
            &(1.0_JPRB-FOEALFCU(ZPT))*R5ALSCP*ZI**2
          ZCOND1=(ZPQ*ZCOR**2-ZQSAT*ZCOR)/(ZCOR**2+ZQSAT*ZF)
        ZPT = ZPT+FOELDCPMCU(ZPT)*ZCOND1
        PT(JL,KK) = ZPT
        PQ(JL,KK)=ZPQ-ZCOND1
      ENDDO
    ENDIF

  ENDIF

  IF(KCALL == 2) THEN

!DIR$    IVDEP
!OCL NOVREC
    DO JL=KIDIA,KFDIA
      IF(LDFLAG(JL)) THEN
        ZQP    =1.0_JPRB/PSP(JL)
        ZQSAT=FOEEWMCU(PT(JL,KK))*ZQP    
        ZQSAT=MIN(0.5_JPRB,ZQSAT)
        ZCOR=1.0_JPRB/(1.0_JPRB-RETV  *ZQSAT)
        ZQSAT=ZQSAT*ZCOR
        ZCOND=(PQ(JL,KK)-ZQSAT)/(1.0_JPRB+ZQSAT*ZCOR*FOEDEMCU(PT(JL,KK)))
        ZCOND=MIN(ZCOND,0.0_JPRB)
        PT(JL,KK)=PT(JL,KK)+FOELDCPMCU(PT(JL,KK))*ZCOND
        PQ(JL,KK)=PQ(JL,KK)-ZCOND
        ZQSAT=FOEEWMCU(PT(JL,KK))*ZQP    
        ZQSAT=MIN(0.5_JPRB,ZQSAT)
        ZCOR=1.0_JPRB/(1.0_JPRB-RETV  *ZQSAT)
        ZQSAT=ZQSAT*ZCOR
        ZCOND1=(PQ(JL,KK)-ZQSAT)/(1.0_JPRB+ZQSAT*ZCOR*FOEDEMCU(PT(JL,KK)))
        IF(ZCOND == 0.0_JPRB)ZCOND1=MIN(ZCOND1,0.0_JPRB)
        PT(JL,KK)=PT(JL,KK)+FOELDCPMCU(PT(JL,KK))*ZCOND1
        PQ(JL,KK)=PQ(JL,KK)-ZCOND1
      ENDIF
    ENDDO

  ENDIF

  IF(KCALL == 0) THEN

!DIR$    IVDEP
!OCL NOVREC

    DO JL=KIDIA,KFDIA
      ZQP    =1.0_JPRB/PSP(JL)
      ZQSAT=FOEEWM(PT(JL,KK))*ZQP    
      ZQSAT=MIN(0.5_JPRB,ZQSAT)
      ZCOR=1.0_JPRB/(1.0_JPRB-RETV  *ZQSAT)
      ZQSAT=ZQSAT*ZCOR
      ZCOND1=(PQ(JL,KK)-ZQSAT)/(1.0_JPRB+ZQSAT*ZCOR*FOEDEM(PT(JL,KK)))
      PT(JL,KK)=PT(JL,KK)+FOELDCPM(PT(JL,KK))*ZCOND1
      PQ(JL,KK)=PQ(JL,KK)-ZCOND1
      ZQSAT=FOEEWM(PT(JL,KK))*ZQP    
      ZQSAT=MIN(0.5_JPRB,ZQSAT)
      ZCOR=1.0_JPRB/(1.0_JPRB-RETV  *ZQSAT)
      ZQSAT=ZQSAT*ZCOR
      ZCOND1=(PQ(JL,KK)-ZQSAT)/(1.0_JPRB+ZQSAT*ZCOR*FOEDEM(PT(JL,KK)))
      PT(JL,KK)=PT(JL,KK)+FOELDCPM(PT(JL,KK))*ZCOND1
      PQ(JL,KK)=PQ(JL,KK)-ZCOND1
    ENDDO

  ENDIF

  IF(KCALL == 4 )THEN

!DIR$    IVDEP
!OCL NOVREC
    DO JL=KIDIA,KFDIA
      IF(LDFLAG(JL)) THEN
        ZQP    =1.0_JPRB/PSP(JL)
        ZQSAT=FOEEWM(PT(JL,KK))*ZQP    
        ZQSAT=MIN(0.5_JPRB,ZQSAT)
        ZCOR=1.0_JPRB/(1.0_JPRB-RETV  *ZQSAT)
        ZQSAT=ZQSAT*ZCOR
        ZCOND=(PQ(JL,KK)-ZQSAT)/(1.0_JPRB+ZQSAT*ZCOR*FOEDEM(PT(JL,KK)))
        PT(JL,KK)=PT(JL,KK)+FOELDCPM(PT(JL,KK))*ZCOND
        PQ(JL,KK)=PQ(JL,KK)-ZCOND
        ZQSAT=FOEEWM(PT(JL,KK))*ZQP    
        ZQSAT=MIN(0.5_JPRB,ZQSAT)
        ZCOR=1.0_JPRB/(1.0_JPRB-RETV  *ZQSAT)
        ZQSAT=ZQSAT*ZCOR
        ZCOND1=(PQ(JL,KK)-ZQSAT)/(1.0_JPRB+ZQSAT*ZCOR*FOEDEM(PT(JL,KK)))
        PT(JL,KK)=PT(JL,KK)+FOELDCPM(PT(JL,KK))*ZCOND1
        PQ(JL,KK)=PQ(JL,KK)-ZCOND1
      ENDIF
    ENDDO
  ENDIF 

  IF(KCALL == 5) THEN  ! Same as 4 but with LDFLAG all true

!DIR$    IVDEP
!OCL NOVREC
    IF(N_VMASS <= 0)  THEN ! Not using Vector MASS
      DO JL=KIDIA,KFDIA
        ZQP    =1.0_JPRB/PSP(JL)
        ZQSAT=FOEEWM(PT(JL,KK))*ZQP    
        ZQSAT=MIN(0.5_JPRB,ZQSAT)
        ZCOR=1.0_JPRB/(1.0_JPRB-RETV  *ZQSAT)
        ZQSAT=ZQSAT*ZCOR
        ZCOND=(PQ(JL,KK)-ZQSAT)/(1.0_JPRB+ZQSAT*ZCOR*FOEDEM(PT(JL,KK)))
        PT(JL,KK)=PT(JL,KK)+FOELDCPM(PT(JL,KK))*ZCOND
        PQ(JL,KK)=PQ(JL,KK)-ZCOND
        ZQSAT=FOEEWM(PT(JL,KK))*ZQP    
        ZQSAT=MIN(0.5_JPRB,ZQSAT)
        ZCOR=1.0_JPRB/(1.0_JPRB-RETV  *ZQSAT)
        ZQSAT=ZQSAT*ZCOR
        ZCOND1=(PQ(JL,KK)-ZQSAT)/(1.0_JPRB+ZQSAT*ZCOR*FOEDEM(PT(JL,KK)))
        PT(JL,KK)=PT(JL,KK)+FOELDCPM(PT(JL,KK))*ZCOND1
        PQ(JL,KK)=PQ(JL,KK)-ZCOND1
      ENDDO
    ELSE ! Using Vector VMASS
      DO JL=KIDIA,KFDIA
!       ZTMP1(JL-KIDIA+1)=R3LES*(PT(JL,KK)-RTT)
!       ZTMP2(JL-KIDIA+1)=R3IES*(PT(JL,KK)-RTT)
!       ZTMP3(JL-KIDIA+1)=PT(JL,KK)-R4LES
!       ZTMP4(JL-KIDIA+1)=PT(JL,KK)-R4IES
        ZTMP1(2*(JL-KIDIA)+1)=R3LES*(PT(JL,KK)-RTT)
        ZTMP1(2*(JL-KIDIA)+2)=R3IES*(PT(JL,KK)-RTT)
        ZTMP3(2*(JL-KIDIA)+1)=PT(JL,KK)-R4LES
        ZTMP3(2*(JL-KIDIA)+2)=PT(JL,KK)-R4IES
      ENDDO
!     CALL VDIV(ZTMP5,ZTMP1,ZTMP3,JLEN)
      CALL VDIV(ZTMP5,ZTMP1,ZTMP3,2*JLEN)
!     CALL VDIV(ZTMP6,ZTMP2,ZTMP4,JLEN)
!     CALL VEXP(ZTMP1,ZTMP5,JLEN)
      CALL VEXP(ZTMP1,ZTMP5,2*JLEN)
!     CALL VEXP(ZTMP2,ZTMP6,JLEN)
!     CALL VREC(ZTMP5,ZTMP3,JLEN)
      CALL VREC(ZTMP5,ZTMP3,2*JLEN)
!     CALL VREC(ZTMP6,ZTMP4,JLEN)
      DO JL=KIDIA,KFDIA
        ZQP    =1.0_JPRB/PSP(JL)
!       ZQSAT=R2ES*(FOEALFA(PT(JL,KK))*ZTMP1(JL-KIDIA+1)+&
!        & (1.0_JPRB-FOEALFA(PT(JL,KK)))*ZTMP2(JL-KIDIA+1))*ZQP  
        ZQSAT=R2ES*(FOEALFA(PT(JL,KK))*ZTMP1(2*(JL-KIDIA)+1)+&
         & (1.0_JPRB-FOEALFA(PT(JL,KK)))*ZTMP1(2*(JL-KIDIA)+2))*ZQP  
        ZQSAT=MINJ(0.5_JPRB,ZQSAT)
        ZCOR=1.0_JPRB-RETV*ZQSAT
!       ZCOND=(PQ(JL,KK)*ZCOR**2-ZQSAT*ZCOR)/(ZCOR**2+ZQSAT*FOEDEM(PT(JL,KK)))
        ! FOEDEM(PTARE) = FOEALFA(PTARE)*R5ALVCP*(1.0_JPRB/(PTARE-R4LES)**2)+&
        !   &(1.0_JPRB-FOEALFA(PTARE))*R5ALSCP*(1.0_JPRB/(PTARE-R4IES)**2)
!       ZF = FOEALFA(PT(JL,KK))*R5ALVCP*(ZTMP5(JL-KIDIA+1)**2)+&
!         &(1.0_JPRB-FOEALFA(PT(JL,KK)))*R5ALSCP*(ZTMP6(JL-KIDIA+1))**2)
        ZF = FOEALFA(PT(JL,KK))*R5ALVCP*(ZTMP5(2*(JL-KIDIA)+1)**2)+&
          &(1.0_JPRB-FOEALFA(PT(JL,KK)))*R5ALSCP*(ZTMP5(2*(JL-KIDIA)+2)**2)
        ZCOND=(PQ(JL,KK)*ZCOR**2-ZQSAT*ZCOR)/(ZCOR**2+ZQSAT*ZF)
        PT(JL,KK)=PT(JL,KK)+FOELDCPM(PT(JL,KK))*ZCOND
        PQ(JL,KK)=PQ(JL,KK)-ZCOND
        ZTMP0(JL-KIDIA+1)=ZQP
!       ZTMP1(JL-KIDIA+1)=R3LES*(PT(JL,KK)-RTT)
!       ZTMP2(JL-KIDIA+1)=R3IES*(PT(JL,KK)-RTT)
!       ZTMP3(JL-KIDIA+1)=PT(JL,KK)-R4LES
!       ZTMP4(JL-KIDIA+1)=PT(JL,KK)-R4IES
        ZTMP1(2*(JL-KIDIA)+1)=R3LES*(PT(JL,KK)-RTT)
        ZTMP1(2*(JL-KIDIA)+2)=R3IES*(PT(JL,KK)-RTT)
        ZTMP3(2*(JL-KIDIA)+1)=PT(JL,KK)-R4LES
        ZTMP3(2*(JL-KIDIA)+2)=PT(JL,KK)-R4IES
      ENDDO
!     CALL VDIV(ZTMP5,ZTMP1,ZTMP3,JLEN)
      CALL VDIV(ZTMP5,ZTMP1,ZTMP3,2*JLEN)
!     CALL VDIV(ZTMP6,ZTMP2,ZTMP4,JLEN)
!     CALL VEXP(ZTMP1,ZTMP5,JLEN)
      CALL VEXP(ZTMP1,ZTMP5,2*JLEN)
!     CALL VEXP(ZTMP2,ZTMP6,JLEN)
!     CALL VREC(ZTMP5,ZTMP3,JLEN)
      CALL VREC(ZTMP5,ZTMP3,2*JLEN)
!     CALL VREC(ZTMP6,ZTMP4,JLEN)
      DO JL=KIDIA,KFDIA
        ZQP  = ZTMP0(JL-KIDIA+1)
!       ZQSAT=R2ES*(FOEALFA(PT(JL,KK))*ZTMP1(JL-KIDIA+1)+&
!        & (1.0_JPRB-FOEALFA(PT(JL,KK)))*ZTMP2(JL-KIDIA+1))*ZQP  
        ZQSAT=R2ES*(FOEALFA(PT(JL,KK))*ZTMP1(2*(JL-KIDIA)+1)+&
         & (1.0_JPRB-FOEALFA(PT(JL,KK)))*ZTMP1(2*(JL-KIDIA)+2))*ZQP  
        ZQSAT=MINJ(0.5_JPRB,ZQSAT)
        ZCOR=1.0_JPRB-RETV*ZQSAT
!       ZCOND1=(PQ(JL,KK)*ZCOR**2-ZQSAT*ZCOR)/(ZCOR**2+ZQSAT*FOEDEM(PT(JL,KK)))
        ! FOEDEM(PTARE) = FOEALFA(PTARE)*R5ALVCP*(1.0_JPRB/(PTARE-R4LES)**2)+&
        !   &(1.0_JPRB-FOEALFA(PTARE))*R5ALSCP*(1.0_JPRB/(PTARE-R4IES)**2)
!       ZF = FOEALFA(PT(JL,KK))*R5ALVCP*(ZTMP5(JL-KIDIA+1)**2)+&
!         &(1.0_JPRB-FOEALFA(PT(JL,KK)))*R5ALSCP*(ZTMP6(JL-KIDIA+1)**2)
        ZF = FOEALFA(PT(JL,KK))*R5ALVCP*(ZTMP5(2*(JL-KIDIA)+1)**2)+&
          &(1.0_JPRB-FOEALFA(PT(JL,KK)))*R5ALSCP*(ZTMP5(2*(JL-KIDIA)+2)**2)
        ZCOND1=(PQ(JL,KK)*ZCOR**2-ZQSAT*ZCOR)/(ZCOR**2+ZQSAT*ZF)
        PT(JL,KK)=PT(JL,KK)+FOELDCPM(PT(JL,KK))*ZCOND1
        PQ(JL,KK)=PQ(JL,KK)-ZCOND1
      ENDDO
    ENDIF
  ENDIF 

  IF(KCALL == 3) THEN 
!DIR$    IVDEP !OCL NOVREC 
    IF(N_VMASS <=  0)  THEN ! Not using Vector MASS
      DO JL=KIDIA,KFDIA
        ZQP    =1.0_JPRB/PSP(JL)
        ZQSAT=FOEEWMCU(PT(JL,KK))*ZQP
        ZQSAT=MIN(0.5_JPRB,ZQSAT)
        ZCOR=1.0_JPRB/(1.0_JPRB-RETV  *ZQSAT)
        ZQSAT=ZQSAT*ZCOR
        ZCOND1=(PQ(JL,KK)-ZQSAT)/(1.0_JPRB+ZQSAT*ZCOR*FOEDEMCU(PT(JL,KK)))
        PT(JL,KK)=PT(JL,KK)+FOELDCPMCU(PT(JL,KK))*ZCOND1
        PQ(JL,KK)=PQ(JL,KK)-ZCOND1
        ZQSAT=FOEEWMCU(PT(JL,KK))*ZQP
        ZQSAT=MIN(0.5_JPRB,ZQSAT)
        ZCOR=1.0_JPRB/(1.0_JPRB-RETV  *ZQSAT)
        ZQSAT=ZQSAT*ZCOR
        ZCOND1=(PQ(JL,KK)-ZQSAT)/(1.0_JPRB+ZQSAT*ZCOR*FOEDEMCU(PT(JL,KK)))
        PT(JL,KK)=PT(JL,KK)+FOELDCPMCU(PT(JL,KK))*ZCOND1
        PQ(JL,KK)=PQ(JL,KK)-ZCOND1
      ENDDO
    ELSE
      DO JL=KIDIA,KFDIA 
!       ZTMP1(JL-KIDIA+1)=R3LES*(PT(JL,KK)-RTT) 
!       ZTMP2(JL-KIDIA+1)=R3IES*(PT(JL,KK)-RTT) 
!       ZTMP3(JL-KIDIA+1)=PT(JL,KK)-R4LES 
!       ZTMP4(JL-KIDIA+1)=PT(JL,KK)-R4IES 
        ZTMP1(2*(JL-KIDIA)+1)=R3LES*(PT(JL,KK)-RTT)
        ZTMP1(2*(JL-KIDIA)+2)=R3IES*(PT(JL,KK)-RTT)
        ZTMP3(2*(JL-KIDIA)+1)=PT(JL,KK)-R4LES 
        ZTMP3(2*(JL-KIDIA)+2)=PT(JL,KK)-R4IES 
      ENDDO 
!     CALL VDIV(ZTMP5,ZTMP1,ZTMP3,JLEN)
!     CALL VDIV(ZTMP6,ZTMP2,ZTMP4,JLEN)
!     CALL VEXP(ZTMP1,ZTMP5,JLEN)
!     CALL VEXP(ZTMP2,ZTMP6,JLEN)
      CALL VDIV(ZTMP5,ZTMP1,ZTMP3,2*JLEN)
      CALL VEXP(ZTMP1,ZTMP5,2*JLEN)
      DO JL=KIDIA,KFDIA
        ZQP    =1.0_JPRB/PSP(JL)
        ZQSAT=R2ES*(FOEALFCU(PT(JL,KK))*ZTMP1(2*(JL-KIDIA)+1)+&
         & (1.0_JPRB-FOEALFCU(PT(JL,KK)))*ZTMP1(2*(JL-KIDIA)+2))*ZQP  
        ZQSAT=MINJ(0.5_JPRB,ZQSAT)
        ZCOR=1.0_JPRB-RETV*ZQSAT
        ZCOND1=(PQ(JL,KK)*ZCOR**2-ZQSAT*ZCOR)/(ZCOR**2+ZQSAT*FOEDEMCU(PT(JL,KK)))
        PT(JL,KK)=PT(JL,KK)+FOELDCPMCU(PT(JL,KK))*ZCOND1
        PQ(JL,KK)=PQ(JL,KK)-ZCOND1
        ZTMP0(JL-KIDIA+1)=ZQP
!       ZTMP1(JL-KIDIA+1)=R3LES*(PT(JL,KK)-RTT)
!       ZTMP2(JL-KIDIA+1)=R3IES*(PT(JL,KK)-RTT)
!       ZTMP3(JL-KIDIA+1)=PT(JL,KK)-R4LES
!       ZTMP4(JL-KIDIA+1)=PT(JL,KK)-R4IES
        ZTMP1(2*(JL-KIDIA)+1)=R3LES*(PT(JL,KK)-RTT)
        ZTMP1(2*(JL-KIDIA)+2)=R3IES*(PT(JL,KK)-RTT)
        ZTMP3(2*(JL-KIDIA)+1)=PT(JL,KK)-R4LES
        ZTMP3(2*(JL-KIDIA)+2)=PT(JL,KK)-R4IES
      ENDDO
!     CALL VDIV(ZTMP5,ZTMP1,ZTMP3,JLEN)
!     CALL VDIV(ZTMP6,ZTMP2,ZTMP4,JLEN)
!     CALL VEXP(ZTMP1,ZTMP5,JLEN)
!     CALL VEXP(ZTMP2,ZTMP6,JLEN)
      CALL VDIV(ZTMP5,ZTMP1,ZTMP3,2*JLEN)
      CALL VEXP(ZTMP1,ZTMP5,2*JLEN)
      DO JL=KIDIA,KFDIA
        ZQP  = ZTMP0(JL-KIDIA+1)
        ZQSAT=R2ES*(FOEALFCU(PT(JL,KK))*ZTMP1(2*(JL-KIDIA)+1)+&
         & (1.0_JPRB-FOEALFCU(PT(JL,KK)))*ZTMP1(2*(JL-KIDIA)+2))*ZQP  
        ZQSAT=MINJ(0.5_JPRB,ZQSAT)
        ZCOR=1.0_JPRB-RETV*ZQSAT
        ZCOND1=(PQ(JL,KK)*ZCOR**2-ZQSAT*ZCOR)/(ZCOR**2+ZQSAT*FOEDEMCU(PT(JL,KK)))
        PT(JL,KK)=PT(JL,KK)+FOELDCPMCU(PT(JL,KK))*ZCOND1
        PQ(JL,KK)=PQ(JL,KK)-ZCOND1
      ENDDO
    ENDIF

  ENDIF
!*********************************************
ELSE
!*********************************************                 

!     2.           CALCULATE CONDENSATION AND ADJUST T AND Q ACCORDINGLY
!                  -----------------------------------------------------

  IF (KCALL == 1 ) THEN

!DIR$    IVDEP
!OCL NOVREC
    DO JL=KIDIA,KFDIA
      IF(LDFLAG(JL)) THEN
        ZQP    =1.0_JPRB/PSP(JL)
        ZTARG=PT(JL,KK)
        ZOEALFA=0.5_JPRB*(TANH(RLPAL1*(ZTARG-RLPTRC))+1.0_JPRB)
        ZFOEEWL=R2ES*EXP(R3LES*(ZTARG-RTT)/(ZTARG-R4LES))
        ZFOEEWI=R2ES*EXP(R3IES*(ZTARG-RTT)/(ZTARG-R4IES))
        ZQSAT=ZQP    *(ZOEALFA*ZFOEEWL+(1.0_JPRB-ZOEALFA)*ZFOEEWI)
        Z1S=TANH(RLPAL2*(ZQSAT-ZQMAX))
        ZQSAT=0.5_JPRB*((1.0_JPRB-Z1S)*ZQSAT+(1.0_JPRB+Z1S)*ZQMAX)

        ZCOR=1.0_JPRB/(1.0_JPRB-RETV  *ZQSAT)
        ZQSAT=ZQSAT*ZCOR

        Z2S=    ZOEALFA *R5ALVCP*(1.0_JPRB/(ZTARG-R4LES)**2)+&
         & (1.0_JPRB-ZOEALFA)*R5ALSCP*(1.0_JPRB/(ZTARG-R4IES)**2)  
        ZCOND=(PQ(JL,KK)-ZQSAT)/(1.0_JPRB+ZQSAT*ZCOR*Z2S)

        ZCOND=MAX(ZCOND,0.0_JPRB)

        IF(ZCOND /= 0.0_JPRB) THEN

          PT(JL,KK)=PT(JL,KK)+&
           & (ZOEALFA*RALVDCP+(1.0_JPRB-ZOEALFA)*RALSDCP)*ZCOND  
          PQ(JL,KK)=PQ(JL,KK)-ZCOND
          ZTARG=PT(JL,KK)
          ZOEALFA=0.5_JPRB*(TANH(RLPAL1*(ZTARG-RLPTRC))+1.0_JPRB)
          ZFOEEWL=R2ES*EXP(R3LES*(ZTARG-RTT)/(ZTARG-R4LES))
          ZFOEEWI=R2ES*EXP(R3IES*(ZTARG-RTT)/(ZTARG-R4IES))
          ZQSAT=ZQP    *(ZOEALFA*ZFOEEWL+(1.0_JPRB-ZOEALFA)*ZFOEEWI)
          Z1S=TANH(RLPAL2*(ZQSAT-ZQMAX))
          ZQSAT=0.5_JPRB*((1.0_JPRB-Z1S)*ZQSAT+(1.0_JPRB+Z1S)*ZQMAX)
  
          ZCOR=1.0_JPRB/(1.0_JPRB-RETV  *ZQSAT)
          ZQSAT=ZQSAT*ZCOR
  
          Z2S=    ZOEALFA *R5ALVCP*(1.0_JPRB/(ZTARG-R4LES)**2)+&
           & (1.0_JPRB-ZOEALFA)*R5ALSCP*(1.0_JPRB/(ZTARG-R4IES)**2)  
          ZCOND1=(PQ(JL,KK)-ZQSAT)/(1.0_JPRB+ZQSAT*ZCOR*Z2S)
  
          PT(JL,KK)=PT(JL,KK)+(ZOEALFA*RALVDCP+(1.0_JPRB-ZOEALFA)*RALSDCP)*ZCOND1
  
          PQ(JL,KK)=PQ(JL,KK)-ZCOND1
        ENDIF
      ENDIF
    ENDDO

  ENDIF

  IF(KCALL == 2) THEN

!DIR$    IVDEP
!OCL NOVREC
    DO JL=KIDIA,KFDIA
      IF(LDFLAG(JL)) THEN
        ZQP    =1.0_JPRB/PSP(JL)

        ZTARG=PT(JL,KK)
        ZOEALFA=0.5_JPRB*(TANH(RLPAL1*(ZTARG-RLPTRC))+1.0_JPRB)
        ZFOEEWL=R2ES*EXP(R3LES*(ZTARG-RTT)/(ZTARG-R4LES))
        ZFOEEWI=R2ES*EXP(R3IES*(ZTARG-RTT)/(ZTARG-R4IES))
        ZQSAT=ZQP    *(ZOEALFA*ZFOEEWL+(1.0_JPRB-ZOEALFA)*ZFOEEWI)
        Z1S=TANH(RLPAL2*(ZQSAT-ZQMAX))
        ZQSAT=0.5_JPRB*((1.0_JPRB-Z1S)*ZQSAT+(1.0_JPRB+Z1S)*ZQMAX)

        ZCOR=1.0_JPRB/(1.0_JPRB-RETV  *ZQSAT)
        ZQSAT=ZQSAT*ZCOR

        Z2S=    ZOEALFA *R5ALVCP*(1.0_JPRB/(ZTARG-R4LES)**2)+&
         & (1.0_JPRB-ZOEALFA)*R5ALSCP*(1.0_JPRB/(ZTARG-R4IES)**2)  
        ZCOND=(PQ(JL,KK)-ZQSAT)/(1.0_JPRB+ZQSAT*ZCOR*Z2S)

        ZCOND=MIN(ZCOND,0.0_JPRB)

        IF(ZCOND /= 0.0_JPRB) THEN

          PT(JL,KK)=PT(JL,KK)+&
           & (ZOEALFA*RALVDCP+(1.0_JPRB-ZOEALFA)*RALSDCP)*ZCOND  
          PQ(JL,KK)=PQ(JL,KK)-ZCOND
          ZTARG=PT(JL,KK)
          ZOEALFA=0.5_JPRB*(TANH(RLPAL1*(ZTARG-RLPTRC))+1.0_JPRB)
          ZFOEEWL=R2ES*EXP(R3LES*(ZTARG-RTT)/(ZTARG-R4LES))
          ZFOEEWI=R2ES*EXP(R3IES*(ZTARG-RTT)/(ZTARG-R4IES))
          ZQSAT=ZQP    *(ZOEALFA*ZFOEEWL+(1.0_JPRB-ZOEALFA)*ZFOEEWI)
          Z1S=TANH(RLPAL2*(ZQSAT-ZQMAX))
          ZQSAT=0.5_JPRB*((1.0_JPRB-Z1S)*ZQSAT+(1.0_JPRB+Z1S)*ZQMAX)
  
          ZCOR=1.0_JPRB/(1.0_JPRB-RETV  *ZQSAT)
          ZQSAT=ZQSAT*ZCOR
  
          Z2S=    ZOEALFA *R5ALVCP*(1.0_JPRB/(ZTARG-R4LES)**2)+&
           & (1.0_JPRB-ZOEALFA)*R5ALSCP*(1.0_JPRB/(ZTARG-R4IES)**2)  
          ZCOND1=(PQ(JL,KK)-ZQSAT)/(1.0_JPRB+ZQSAT*ZCOR*Z2S)
  
          PT(JL,KK)=PT(JL,KK)+(ZOEALFA*RALVDCP+(1.0_JPRB-ZOEALFA)*RALSDCP)*ZCOND1
  
          PQ(JL,KK)=PQ(JL,KK)-ZCOND1
        ENDIF
      ENDIF
    ENDDO

  ENDIF

  IF(KCALL == 0) THEN

!DIR$    IVDEP
!OCL NOVREC
    DO JL=KIDIA,KFDIA
      ZQP    =1.0_JPRB/PSP(JL)

      ZTARG=PT(JL,KK)
      ZOEALFA=0.5_JPRB*(TANH(RLPAL1*(ZTARG-RLPTRC))+1.0_JPRB)
      ZFOEEWL=R2ES*EXP(R3LES*(ZTARG-RTT)/(ZTARG-R4LES))
      ZFOEEWI=R2ES*EXP(R3IES*(ZTARG-RTT)/(ZTARG-R4IES))
      ZQSAT=ZQP    *(ZOEALFA*ZFOEEWL+(1.0_JPRB-ZOEALFA)*ZFOEEWI)
      Z1S=TANH(RLPAL2*(ZQSAT-ZQMAX))
      ZQSAT=0.5_JPRB*((1.0_JPRB-Z1S)*ZQSAT+(1.0_JPRB+Z1S)*ZQMAX)

      ZCOR=1.0_JPRB/(1.0_JPRB-RETV  *ZQSAT)
      ZQSAT=ZQSAT*ZCOR

      Z2S=    ZOEALFA *R5ALVCP*(1.0_JPRB/(ZTARG-R4LES)**2)+&
       & (1.0_JPRB-ZOEALFA)*R5ALSCP*(1.0_JPRB/(ZTARG-R4IES)**2)  
      ZCOND1=(PQ(JL,KK)-ZQSAT)/(1.0_JPRB+ZQSAT*ZCOR*Z2S)

      PT(JL,KK)=PT(JL,KK)+(ZOEALFA*RALVDCP+(1.0_JPRB-ZOEALFA)*RALSDCP)*ZCOND1

      PQ(JL,KK)=PQ(JL,KK)-ZCOND1

      ZTARG=PT(JL,KK)
      ZOEALFA=0.5_JPRB*(TANH(RLPAL1*(ZTARG-RLPTRC))+1.0_JPRB)
      ZFOEEWL=R2ES*EXP(R3LES*(ZTARG-RTT)/(ZTARG-R4LES))
      ZFOEEWI=R2ES*EXP(R3IES*(ZTARG-RTT)/(ZTARG-R4IES))
      ZQSAT=ZQP    *(ZOEALFA*ZFOEEWL+(1.0_JPRB-ZOEALFA)*ZFOEEWI)
      Z1S=TANH(RLPAL2*(ZQSAT-ZQMAX))
      ZQSAT=0.5_JPRB*((1.0_JPRB-Z1S)*ZQSAT+(1.0_JPRB+Z1S)*ZQMAX)

      ZCOR=1.0_JPRB/(1.0_JPRB-RETV  *ZQSAT)
      ZQSAT=ZQSAT*ZCOR

      Z2S=    ZOEALFA *R5ALVCP*(1.0_JPRB/(ZTARG-R4LES)**2)+&
       & (1.0_JPRB-ZOEALFA)*R5ALSCP*(1.0_JPRB/(ZTARG-R4IES)**2)  
      ZCOND1=(PQ(JL,KK)-ZQSAT)/(1.0_JPRB+ZQSAT*ZCOR*Z2S)

      PT(JL,KK)=PT(JL,KK)+(ZOEALFA*RALVDCP+(1.0_JPRB-ZOEALFA)*RALSDCP)*ZCOND1

      PQ(JL,KK)=PQ(JL,KK)-ZCOND1
    ENDDO

  ENDIF

  IF(KCALL == 4) THEN

!DIR$    IVDEP
!OCL NOVREC
    DO JL=KIDIA,KFDIA
      IF(LDFLAG(JL)) THEN
        ZQP    =1.0_JPRB/PSP(JL)

        ZTARG=PT(JL,KK)
        ZOEALFA=0.5_JPRB*(TANH(RLPAL1*(ZTARG-RLPTRC))+1.0_JPRB)
        ZFOEEWL=R2ES*EXP(R3LES*(ZTARG-RTT)/(ZTARG-R4LES))
        ZFOEEWI=R2ES*EXP(R3IES*(ZTARG-RTT)/(ZTARG-R4IES))
        ZQSAT=ZQP    *(ZOEALFA*ZFOEEWL+(1.0_JPRB-ZOEALFA)*ZFOEEWI)
        Z1S=TANH(RLPAL2*(ZQSAT-ZQMAX))
        ZQSAT=0.5_JPRB*((1.0_JPRB-Z1S)*ZQSAT+(1.0_JPRB+Z1S)*ZQMAX)
        
        ZCOR=1.0_JPRB/(1.0_JPRB-RETV  *ZQSAT)
        ZQSAT=ZQSAT*ZCOR
        
        Z2S=    ZOEALFA *R5ALVCP*(1.0_JPRB/(ZTARG-R4LES)**2)+&
         & (1.0_JPRB-ZOEALFA)*R5ALSCP*(1.0_JPRB/(ZTARG-R4IES)**2)  
        ZCOND=(PQ(JL,KK)-ZQSAT)/(1.0_JPRB+ZQSAT*ZCOR*Z2S)
        
        PT(JL,KK)=PT(JL,KK)+(ZOEALFA*RALVDCP+(1.0_JPRB-ZOEALFA)*RALSDCP)*ZCOND
        
        PQ(JL,KK)=PQ(JL,KK)-ZCOND
        
        ZTARG=PT(JL,KK)
        ZOEALFA=0.5_JPRB*(TANH(RLPAL1*(ZTARG-RLPTRC))+1.0_JPRB)
        ZFOEEWL=R2ES*EXP(R3LES*(ZTARG-RTT)/(ZTARG-R4LES))
        ZFOEEWI=R2ES*EXP(R3IES*(ZTARG-RTT)/(ZTARG-R4IES))
        ZQSAT=ZQP    *(ZOEALFA*ZFOEEWL+(1.0_JPRB-ZOEALFA)*ZFOEEWI)
        Z1S=TANH(RLPAL2*(ZQSAT-ZQMAX))
        ZQSAT=0.5_JPRB*((1.0_JPRB-Z1S)*ZQSAT+(1.0_JPRB+Z1S)*ZQMAX)
        
        ZQSAT=MIN(ZQMAX,ZQSAT)
        ZCOR=1.0_JPRB/(1.0_JPRB-RETV  *ZQSAT)
        ZQSAT=ZQSAT*ZCOR
        
        Z2S=    ZOEALFA *R5ALVCP*(1.0_JPRB/(ZTARG-R4LES)**2)+&
         & (1.0_JPRB-ZOEALFA)*R5ALSCP*(1.0_JPRB/(ZTARG-R4IES)**2)  
        ZCOND1=(PQ(JL,KK)-ZQSAT)/(1.0_JPRB+ZQSAT*ZCOR*Z2S)
        
        PT(JL,KK)=PT(JL,KK)+(ZOEALFA*RALVDCP+(1.0_JPRB-ZOEALFA)*RALSDCP)*ZCOND1

        PQ(JL,KK)=PQ(JL,KK)-ZCOND1
      ENDIF
    ENDDO

  ENDIF

!*********************************************
ENDIF
!*********************************************                 

IF (LHOOK) CALL DR_HOOK('CUADJTQ',1,ZHOOK_HANDLE)
END SUBROUTINE CUADJTQ2


END MODULE mo_cuadjtq_2
