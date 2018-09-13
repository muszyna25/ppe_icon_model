!>
!! Saturation routine for EDMF DUALM
!!
!! @author Martin Koehler, DWD
!!
!! @par Revision History
!! Imported from IFS by Martin Koehler  (starting 2011-9-30)
!!   (IFS cycle CY36R1_DUALM_M8b)
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

MODULE mo_vdfthermo
 
  PUBLIC :: vdfthermo

CONTAINS

!! !OPTIONS XOPT(HSFUN)
SUBROUTINE VDFTHERMO(thl,qt,qs,ql,pf,gz)

!----------------------------------------------------------------------!
!                                                                      !
!           Sommeria-Deardorff (JAS, 1977) saturation routine          !
!           -------------------------------------------------          !
!                                                                      !
! 4-step iteration                                                     !
!                                                                      !
! qt has to be in g/kg, pf in Pa!                                      !
! qs, ql are also returned in g/kg!                                    !
!                                                                      !
! If the input qt <=0., it is assumed the dry saturation curve         !
! (qsl) is  asked for. In that case, no iteration is required.         !
!                                                                      !
! This routine is also used in the KNMI LES model.                     !
!                                                                      !
!                                                                      !
!     PARAMETER     DESCRIPTION                          UNITS         !
!     ---------     -----------                          -----         !
!     INPUT PARAMETERS (INTEGER):                                      !
!                                                                      !
!                                                                      !
!     INPUT PARAMETERS (REAL):                                         !
!                                                                      !
!      *PF*  pressure                                                  !
!      *GZ*  geopotential                                              !
!                                                                      !
!     OUTPUT PARAMETERS (INTEGER):                                     !
!                                                                      !
!      *QS*  saturation specific humidity                              !
!      *QL*  condensate                                                !
!                                                                      !
! Roel Neggers, ECMWF    22-02-2006                                    !
!                                                                      !
!----------------------------------------------------------------------!
      
! USE PARKIND1  ,ONLY : JPIM     ,JPRB
! USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
! USE YOMCST   , ONLY : RD, RV, RCPD, RETV, RTT, RLVTT
! USE YOETHF   , ONLY : R2ES,  R3LES,  R4LES  

!ICON definitions:
USE mo_kind         ,ONLY : JPRB=>wp ,JPIM=>i4
USE mo_cuparameters ,ONLY : lhook    ,dr_hook ,&
                   & RD    ,RV       ,RCPD    ,RETV ,RTT ,RLVTT ,& !yomcst
                   & R2ES  ,R3LES    ,R4LES                        !yoethf

IMPLICIT NONE


!*         0.1    GLOBAL VARIABLES

REAL(KIND=JPRB)   ,INTENT(IN)    :: thl,qt,pf,gz
REAL(KIND=JPRB)   ,INTENT(OUT)   :: qs,ql

!*         0.2    LOCAL VARIABLES
      
REAL(KIND=JPRB) :: tlin,tl,qtin,qttemp,es,qsl,b1

REAL(KIND=JPRB) :: tmelt, es0, at, bt, cp, rdry, rlv

INTEGER(KIND=JPIM) :: l,liter
  
REAL(KIND=JPRB) :: ZHOOK_HANDLE

LOGICAL :: LLKNMI



!     ------------------------------------------------------------------
!
!*         1.     INITIALIZE VARIABLES
!                 --------------------
!

IF (LHOOK) CALL DR_HOOK('VDFTHERMO',0,ZHOOK_HANDLE)

LLKNMI = .FALSE.
!LLKNMI = .TRUE.


!-- KNMI LES constants ---
tmelt  = 273.16_JPRB   !RTT   = 273.16
es0    = 610.78_JPRB   !R2ES  = 611.21*RD/RV
at     = 17.27_JPRB    !R3LES = 17.502
bt     = 35.86_JPRB    !R4LES = 32.19
cp     = 1004._JPRB    !RCPD  = 1004.7088578
rdry   = 287._JPRB     !RD    = 287.05967
rlv    = 2.5e6_JPRB    !RLVTT = 2500800.
!       = 0.378_JPRB    !      = RETV * RD / RV


!-- initialize --
tlin   = thl - gz/RCPD
qtin   = qt / 1000._JPRB

qttemp = qtin
tl     = tlin      


      
!     ------------------------------------------------------------------
!
!*         2.     ITERATION
!                 ---------
!      
      
liter = 4
IF (qt.le.0._JPRB) THEN 
  liter=1
ENDIF
        
DO l=1,liter

  IF (LLKNMI) THEN
  
    !-- KNMI LES--
    es  = es0*exp(at*(tl-tmelt)/(tl-bt))
    qsl = 0.622_JPRB*es/(pf-0.378_JPRB*es)
    b1  = 0.622_JPRB*rlv/(rdry*tl)*(rlv/(cp*tl))
    qs  = qsl*(1._JPRB+b1*qttemp)/(1._JPRB+b1*qsl)
  
  ELSE

    !-- using ECMWF constants --  
    es  = R2ES*exp(R3LES*(tl-RTT)/(tl-R4LES))
    qsl = es/(pf-RETV*es)
    b1  = (RD/RV)*RLVTT/(RD*tl)*(RLVTT/(RCPD*tl))
    qs  = qsl*(1._JPRB+b1*qttemp)/(1._JPRB+b1*qsl)
  
  ENDIF
      
  ql = dim(qtin-qs,0._JPRB)
        
  !-- set point to this guess of {theta,qsat} --
  qttemp = max(qsl,qs)  
!  tl     = tlin + rlv * ql / cp
  tl     = tlin + RLVTT * ql / RCPD
        
ENDDO
      
qs = qttemp * 1000._JPRB
ql = ql * 1000._JPRB
      

IF (LHOOK) CALL DR_HOOK('VDFTHERMO',1,ZHOOK_HANDLE)

RETURN
END SUBROUTINE VDFTHERMO

END MODULE mo_vdfthermo
