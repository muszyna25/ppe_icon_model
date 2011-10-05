!>
!! Diagnostic PBL height for EDMF DUALM
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

MODULE mo_vdfdpbl
 
  PUBLIC :: vdfdpbl

CONTAINS

!! !OPTIONS XOPT(HSFUN)
SUBROUTINE VDFDPBL(KIDIA,KFDIA,KLON,KLEV,&
 & PUM1,PVM1,PTM1,PQM1,PGEOM1,&
 & PKMFL,PKHFL,PKQFL,PDHPBL)  
!     ------------------------------------------------------------------

!**   *VDFDPBL* - VDFDPBL (Diagnostic PBL height) determines  
!                 PBL height for diagnostic purposes only

!     A.C.M. BELJAARS       E.C.M.W.F.    17/02/1998.

!     PURPOSE
!     -------

!     Determine PBL height

!     INTERFACE
!     ---------

!     *VDFDPBL* IS CALLED BY *VDFMAIN*

!     INPUT PARAMETERS (INTEGER):

!     *KIDIA*        START POINT
!     *KFDIA*        END POINT
!     *KLEV*         NUMBER OF LEVELS
!     *KLON*         NUMBER OF GRID POINTS PER PACKET

!     INPUT PARAMETERS (REAL):

!     *PUM1*         X-VELOCITY COMPONENT AT T-1
!     *PVM1*         Y-VELOCITY COMPONENT AT T-1
!     *PTM1*         TEMPERATURE AT T-1
!     *PQM1*         SPECIFIC HUMUDITY AT T-1
!     *PGEOM1*       GEOPOTENTIAL AT T-1
!     *PKMFL*        KINEMATIC MOMENTUM FLUX                
!     *PKHFL*        KINEMATIC HEAT FLUX                    
!     *PKQFL*        KINEMATIC MOISTURE FLUX           

!     OUTPUT PARAMETERS (REAL):

!     *PDHPBL*        Boundary layer height                  m

!     METHOD
!     ------

!     Troen and Mahrt method using bulk Richardson criterion
!     (see documentation)

!     ------------------------------------------------------------------

! USE PARKIND1  ,ONLY : JPIM     ,JPRB
! USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
! USE YOMCST   , ONLY : RG       ,RCPD     ,RETV
! USE YOEVDF   , ONLY : REPDU2

!ICON definitions:
USE mo_kind         ,ONLY : JPRB=>wp ,JPIM=>i4
USE mo_cuparameters ,ONLY : lhook    ,dr_hook  ,&
                & RG       ,RCPD     ,RETV     ,& !yomcst
                & REPDU2                          !yoevdf
USE mo_edmf_param   ,ONLY : &
                & REPUST                          !yos_exc

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGEOM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PKMFL(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PKHFL(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PKQFL(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDHPBL(KLON) 
!*            LOCAL STORAGE
!             ----- -------

LOGICAL :: LLDONE(KLON)
REAL(KIND=JPRB) ::    ZRI(KLON),ZDU2(KLON),ZUST(KLON),&
 & ZSVBOT(KLON),ZSVBOTP(KLON),ZSVINC(KLON),ZKHVFL(KLON),&
 & ZWS(KLON)  

INTEGER(KIND=JPIM) :: ILEVM1, ITOT, JIT, JK, JL

REAL(KIND=JPRB) :: ZBUST, ZCONS13, ZCONS14, ZCONS15, ZDRORO,&
 & ZEPS, ZPAR, ZPAR1, ZPARZI, ZRICRI, ZRILEV, &
 & ZSV,ZREPUST
REAL(KIND=JPRB) :: ZHOOK_HANDLE

INTERFACE
! #include "surf_inq.h"
END INTERFACE

!     ------------------------------------------------------------------

!*       1.     INITIALIZE CONSTANTS
!               ---------- ----------

IF (LHOOK) CALL DR_HOOK('VDFDPBL',0,ZHOOK_HANDLE)
ZPAR   = 8.5_JPRB
ZPAR1  = 0.6_JPRB
ZCONS14= RCPD*ZPAR
ZCONS15= ZPAR1*RG
ZEPS   = 1.E-10_JPRB
ZBUST  = 100._JPRB
ZRICRI = 0.25_JPRB
ZPARZI = 1000._JPRB
ZEPS   = 0.5_JPRB*RCPD

!MK: convert to use statement
!  CALL SURF_INQ(PREPUST=ZREPUST)
ZREPUST = REPUST

ZCONS13=1.0_JPRB/3._JPRB

ILEVM1 = KLEV-1

!        1.    PREPARE SURFACE PARAMETERS
!              ------- ------- ----------

DO JL=KIDIA,KFDIA
  ZUST(JL)=MAX(SQRT(PKMFL(JL)),ZREPUST)
  ZKHVFL(JL)=PKHFL(JL)+RETV*PTM1(JL,KLEV)*PKQFL(JL)
  ZSVBOT(JL)=RCPD*PTM1(JL,KLEV)*(1.0_JPRB+RETV*PQM1(JL,KLEV))+PGEOM1(JL,KLEV)
  PDHPBL(JL)=ZPARZI
ENDDO

!        2.    DO 1 OR 2 ITERATIONS ON PBL HEIGHT
!              -- - -- - ---------- -- --- ------

!*************************************
DO JIT=1,2
!*************************************

!        2.1   UPDATE VELOCITY SCALE AND EXCESS TEMPERATURE
!              ------ -------- ----- --- ------ -----------

  DO JL=KIDIA,KFDIA
    LLDONE(JL)=.FALSE.
    ZRI(JL)=0.0_JPRB
    IF (ZKHVFL(JL)  <  0.0_JPRB) THEN
      ZWS(JL)=(ZUST(JL)**3 &
       & -ZCONS15*ZKHVFL(JL)*PDHPBL(JL)/PTM1(JL,KLEV))**ZCONS13  
      ZSVINC(JL)=-ZCONS14*ZKHVFL(JL)/ZWS(JL)
    ELSE
      ZWS(JL)=ZUST(JL)
      ZSVINC(JL)=0.0_JPRB
    ENDIF
    ZSVBOTP(JL)=ZSVBOT(JL)+ZSVINC(JL)+ZEPS
  ENDDO

!        2.2    VERTICAL SCAN TO DETERMINE MIXED LAYER DEPTH
!               -------- ---- -- --------- ----- ----- -----

!***
  DO JK=ILEVM1,1,-1
!***

    ITOT=KFDIA-KIDIA+1
    DO JL=KIDIA,KFDIA
      IF (.NOT. LLDONE(JL)) THEN
        ZSV=RCPD*PTM1(JL,JK)*(1.0_JPRB+RETV*PQM1(JL,JK))+PGEOM1(JL,JK)
        ZDU2(JL)=MAX(REPDU2, (PUM1(JL,JK)-PUM1(JL,KLEV))**2 &
         & +(PVM1(JL,JK)-PVM1(JL,KLEV))**2 &
         & +ZBUST*ZUST(JL)**2 )  
        ZDRORO=2.0_JPRB*(ZSV-ZSVBOTP(JL))&
         & /(ZSV+ZSVBOT(JL)-PGEOM1(JL,JK)-PGEOM1(JL,KLEV))  
        ZRILEV=(PGEOM1(JL,JK)-PGEOM1(JL,KLEV))*ZDRORO/ZDU2(JL)
        IF (ZRILEV  >  ZRICRI) THEN
          PDHPBL(JL)=( (ZRILEV-ZRICRI)*PGEOM1(JL,JK+1)&
           & +(ZRICRI-ZRI(JL))*PGEOM1(JL,JK) )/&
           & ((ZRILEV-ZRI(JL))*RG)  
          LLDONE(JL)=.TRUE.
          ITOT=ITOT-1
        ELSE
          ZRI(JL)=ZRILEV
        ENDIF
      ELSE
        ITOT=ITOT-1
      ENDIF
    ENDDO
    IF (ITOT  <=  0) EXIT
!***
  ENDDO
!***

!*************************************
ENDDO
!*************************************

IF (LHOOK) CALL DR_HOOK('VDFDPBL',1,ZHOOK_HANDLE)
END SUBROUTINE VDFDPBL


END MODULE mo_vdfdpbl
