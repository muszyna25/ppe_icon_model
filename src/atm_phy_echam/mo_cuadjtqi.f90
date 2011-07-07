!>
!! Module contains subroutine cuadjtqi.
!!
!! @author M. Tiedtke, ECMWF
!! @author D. Salmond, Cray UK
!!
!! @par Revision History
!! M. Tiedtke, ECMWF (1989)
!! D. Salmond, Cray UK (1991)
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
MODULE mo_cuadjtqi

  USE mo_kind,               ONLY: dp

#ifdef __ICON__
  USE mo_physical_constants, ONLY: vtmpc1, tmelt
  USE mo_echam_cloud_params, ONLY: cthomi, csecfrl
#else
  USE mo_constants,          ONLY: vtmpc1, tmelt
  USE mo_cloud,              ONLY: cthomi, csecfrl
#endif

  USE mo_convect_tables,     ONLY: tlucua,   & ! table a
                                   tlucub,   & ! table b
                                   tlucuc,   & ! table c
                                   jptlucu1, jptlucu2, &
                                   lookuperror, lookupoverflow &
                                   ,tlucuaw,tlucubw,tlucucw

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: cuadjtqi

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

CONTAINS
  !>
  !!
SUBROUTINE cuadjtqi( ncvmicro, kproma, kbdim, klev, kk,   &
           pp,       pt,       pq,       ldflag,   kcall, &
!--- Included for prognostic CDNC/IC scheme (Ulrike Lohmann, 11/02/2007)
           plui)
!-End Included for prognostic CDNC/IC scheme (Ulrike Lohmann, 11/02/2007)
!
!          M.TIEDTKE         E.C.M.W.F.     12/89
!          D.SALMOND         CRAY(UK))      12/8/91
!
!          PURPOSE.
!          --------
!          TO PRODUCE T,Q AND L VALUES FOR CLOUD ASCENT
!
!          INTERFACE
!          ---------
!          THIS ROUTINE IS CALLED FROM SUBROUTINES:
!              *CUBASE*   (T AND Q AT CONDENSTION LEVEL)
!              *CUASC*    (T AND Q AT CLOUD LEVELS)
!              *CUINI*    (ENVIRONMENTAL T AND QS VALUES AT HALF LEVELS)
!          INPUT ARE UNADJUSTED T AND Q VALUES,
!          IT RETURNS ADJUSTED VALUES OF T AND Q
!          NOTE: INPUT PARAMETER KCALL DEFINES CALCULATION AS
!               KCALL=0    ENV. T AND QS IN*CUINI*
!               KCALL=1  CONDENSATION IN UPDRAFTS  (E.G. CUBASE, CUASC)
!               KCALL=2  EVAPORATION IN DOWNDRAFTS (E.G. CUDLFS,CUDDRAF)
!
!          EXTERNALS
!          ---------
!          3 LOOKUP TABLES ( TLUCUA, TLUCUB, TLUCUC )
!          FOR CONDENSATION CALCULATIONS.
!          THE TABLES ARE INITIALISED IN *SETPHYS*.
!
  !  Scalar arguments with intent(In):
  INTEGER, INTENT (IN) :: ncvmicro     !< 0 or 1. Scheme for convective microphysics
  INTEGER, INTENT (IN) :: kcall, kk, klev, kproma, kbdim

  !  Array arguments with intent(In):
  REAL(dp), INTENT (IN) :: pp(kbdim)
  LOGICAL, INTENT (IN) :: ldflag(kbdim)

  !  Array arguments with intent(InOut):
  REAL(dp), INTENT (INOUT) :: pq(kbdim,klev), pt(kbdim,klev)

  !  Local scalars:
  REAL(dp):: zcond1, zqst1, zdqsdt, zqsat, zes, zcor, zlcdqsdt
  INTEGER :: isum, jl, it, it1
  LOGICAL :: LO

  !  Local arrays:
  REAL(dp):: zcond(kbdim)
!--- Included for prognostic CDNC/IC scheme (Ulrike Lohmann, 11/02/2007)
LOGICAL   :: lo2
REAL(dp)  :: plui(kbdim,klev),ztlucub,ztlucuc
!-------------------------------------end-----------------------------


  !  Executable statements

  lookupoverflow = .FALSE.
  zcond(:)=0._dp
!
!----------------------------------------------------------------------
!
!     2.           CALCULATE CONDENSATION AND ADJUST T AND Q ACCORDINGLY
!                  -----------------------------------------------------
!
  IF (kcall.EQ.1 ) THEN

     isum=0
!DIR$ IVDEP
!OCL NOVREC
     DO 210 jl=1,kproma
        IF(ldflag(jl)) THEN
           it = NINT(pt(jl,kk)*1000._dp)
           IF (it<jptlucu1 .OR. it>jptlucu2) lookupoverflow = .TRUE.
           it = MAX(MIN(it,jptlucu2),jptlucu1)
!--- Included for prognostic CDNC/IC scheme (Ulrike Lohmann, 11/02/2007)
           IF (ncvmicro>0) THEN
              lo2=(pt(jl,kk).LE.cthomi).OR.                              &
                   ((pt(jl,kk).LT.tmelt).AND.(plui(jl,kk).GT.csecfrl))
              zes=MERGE(tlucua(it)/pp(jl),tlucuaw(it)/pp(jl),lo2)
           ELSE
              zes=tlucua(it)/pp(jl)
           ENDIF
!-------------------------------end----------------------------------
           zes=MIN(0.5_dp,zes)
           LO=zes<0.4_dp
           zcor=1._dp/(1._dp-vtmpc1*zes)
           zqsat=zes*zcor
           it1=it+1
           it1 = MAX(MIN(it1,jptlucu2),jptlucu1)
           IF (ncvmicro>0) THEN
              lo2=(pt(jl,kk).LE.cthomi).OR.                              &
                   ((pt(jl,kk).LT.tmelt).AND.(plui(jl,kk).GT.csecfrl))
              zqst1=MERGE(tlucua(it1)/pp(jl),tlucuaw(it1)/pp(jl),lo2)
              zqst1=MIN(0.5_dp,zqst1)
              zqst1=zqst1/(1._dp-vtmpc1*zqst1)
              zdqsdt=(zqst1-zqsat)*1000._dp
              ztlucub=MERGE(tlucub(it),tlucubw(it),lo2)
              ztlucuc=MERGE(tlucuc(it),tlucucw(it),lo2)
              zlcdqsdt=MERGE(zdqsdt*ztlucuc,zqsat*zcor*ztlucub,LO)
              zcond(jl)=(pq(jl,kk)-zqsat)/(1._dp+zlcdqsdt)
              zcond(jl)=MAX(zcond(jl),0._dp)
              pt(jl,kk)=pt(jl,kk)+ztlucuc*zcond(jl)
           ELSE
              zqst1=tlucua(it1)/pp(jl)
              zqst1=MIN(0.5_dp,zqst1)
              zqst1=zqst1/(1._dp-vtmpc1*zqst1)
              zdqsdt=(zqst1-zqsat)*1000._dp
              zlcdqsdt=MERGE(zdqsdt*tlucuc(it),zqsat*zcor*tlucub(it),LO)
              zcond(jl)=(pq(jl,kk)-zqsat)/(1._dp+zlcdqsdt)
              zcond(jl)=MAX(zcond(jl),0._dp)
              pt(jl,kk)=pt(jl,kk)+tlucuc(it)*zcond(jl)
           ENDIF
           pq(jl,kk)=pq(jl,kk)-zcond(jl)
           IF(ABS(zcond(jl)) > 0._dp) isum=isum+1
        END IF
210  END DO

     IF (lookupoverflow) CALL lookuperror ('cuadjtqi (1) ')

     IF(isum.EQ.0) go to 230

!DIR$ IVDEP
!OCL NOVREC
     DO 220 jl=1,kproma
        IF(ldflag(jl)) THEN
           IF(ABS(zcond(jl)) > 0._dp) THEN
           it = NINT(pt(jl,kk)*1000._dp)
           IF (it<jptlucu1 .OR. it>jptlucu2) lookupoverflow = .TRUE.
           it = MAX(MIN(it,jptlucu2),jptlucu1)
!--- Included for prognostic CDNC/IC scheme (Ulrike Lohmann, 11/02/2007)
           IF (ncvmicro>0) THEN
              lo2=(pt(jl,kk).LE.cthomi).OR.                              &
                   ((pt(jl,kk).LT.tmelt).AND.(plui(jl,kk).GT.csecfrl))
              zes=MERGE(tlucua(it)/pp(jl),tlucuaw(it)/pp(jl),lo2)
           ELSE
              zes=tlucua(it)/pp(jl)
           ENDIF
!-------------------------------end----------------------------------
           zes=MIN(0.5_dp,zes)
           LO=zes<0.4_dp
           zcor=1._dp/(1._dp-vtmpc1*zes)
           zqsat=zes*zcor
           it1=it+1
           it1 = MAX(MIN(it1,jptlucu2),jptlucu1)
           IF (ncvmicro>0) THEN
              lo2=(pt(jl,kk).LE.cthomi).OR.                              &
                   ((pt(jl,kk).LT.tmelt).AND.(plui(jl,kk).GT.csecfrl))
              zqst1=MERGE(tlucua(it1)/pp(jl),tlucuaw(it1)/pp(jl),lo2)
              zqst1=MIN(0.5_dp,zqst1)
              zqst1=zqst1/(1._dp-vtmpc1*zqst1)
              zdqsdt=(zqst1-zqsat)*1000._dp
              ztlucub=MERGE(tlucub(it),tlucubw(it),lo2)
              ztlucuc=MERGE(tlucuc(it),tlucucw(it),lo2)
              zlcdqsdt=MERGE(zdqsdt*ztlucuc,zqsat*zcor*ztlucub,LO)
              zcond1=(pq(jl,kk)-zqsat)/(1._dp+zlcdqsdt)
              pt(jl,kk)=pt(jl,kk)+ztlucuc*zcond1
           ELSE
              zqst1=tlucua(it1)/pp(jl)
              zqst1=MIN(0.5_dp,zqst1)
              zqst1=zqst1/(1._dp-vtmpc1*zqst1)
              zdqsdt=(zqst1-zqsat)*1000._dp
              zlcdqsdt=MERGE(zdqsdt*tlucuc(it),zqsat*zcor*tlucub(it),LO)
              zcond1=(pq(jl,kk)-zqsat)/(1._dp+zlcdqsdt)
              pt(jl,kk)=pt(jl,kk)+tlucuc(it)*zcond1
           ENDIF
           pq(jl,kk)=pq(jl,kk)-zcond1
           END IF
        END IF
220  END DO

     IF (lookupoverflow) CALL lookuperror ('cuadjtqi (2) ')

230  CONTINUE

  END IF

  IF(kcall.EQ.2) THEN

     isum=0
!DIR$ IVDEP
!OCL NOVREC
     DO 310 jl=1,kproma
        IF(ldflag(jl)) THEN
           it = NINT(pt(jl,kk)*1000._dp)
           IF (it<jptlucu1 .OR. it>jptlucu2) lookupoverflow = .TRUE.
           it = MAX(MIN(it,jptlucu2),jptlucu1)
!--- Included for prognostic CDNC/IC scheme (Ulrike Lohmann, 11/02/2007)
           IF (ncvmicro>0) THEN
              lo2=(pt(jl,kk).LE.cthomi).OR.                              &
                   ((pt(jl,kk).LT.tmelt).AND.(plui(jl,kk).GT.csecfrl))
              zes=MERGE(tlucua(it)/pp(jl),tlucuaw(it)/pp(jl),lo2)
           ELSE
              zes=tlucua(it)/pp(jl)
           ENDIF
!-------------------------------end----------------------------------
           zes=MIN(0.5_dp,zes)
           LO=zes<0.4_dp
           zcor=1._dp/(1._dp-vtmpc1*zes)
           zqsat=zes*zcor
           it1=it+1
           it1 = MAX(MIN(it1,jptlucu2),jptlucu1)
           IF (ncvmicro>0) THEN
              lo2=(pt(jl,kk).LE.cthomi).OR.                              &
                   ((pt(jl,kk).LT.tmelt).AND.(plui(jl,kk).GT.csecfrl))
              zqst1=MERGE(tlucua(it1)/pp(jl),tlucuaw(it1)/pp(jl),lo2)
              zqst1=MIN(0.5_dp,zqst1)
              zqst1=zqst1/(1._dp-vtmpc1*zqst1)
              zdqsdt=(zqst1-zqsat)*1000._dp
              ztlucub=MERGE(tlucub(it),tlucubw(it),lo2)
              ztlucuc=MERGE(tlucuc(it),tlucucw(it),lo2)
              zlcdqsdt=MERGE(zdqsdt*ztlucuc,zqsat*zcor*ztlucub,LO)
              zcond(jl)=(pq(jl,kk)-zqsat)/(1._dp+zlcdqsdt)
              zcond(jl)=MIN(zcond(jl),0._dp)
              pt(jl,kk)=pt(jl,kk)+ztlucuc*zcond(jl)
           ELSE
              zqst1=tlucua(it1)/pp(jl)
              zqst1=MIN(0.5_dp,zqst1)
              zqst1=zqst1/(1._dp-vtmpc1*zqst1)
              zdqsdt=(zqst1-zqsat)*1000._dp
              zlcdqsdt=MERGE(zdqsdt*tlucuc(it),zqsat*zcor*tlucub(it),LO)
              zcond(jl)=(pq(jl,kk)-zqsat)/(1._dp+zlcdqsdt)
              zcond(jl)=MIN(zcond(jl),0._dp)
              pt(jl,kk)=pt(jl,kk)+tlucuc(it)*zcond(jl)
           ENDIF
           pq(jl,kk)=pq(jl,kk)-zcond(jl)
           IF(ABS(zcond(jl)) > 0._dp) isum=isum+1
        END IF
310  END DO

     IF (lookupoverflow) CALL lookuperror ('cuadjtqi (3) ')

     IF(isum.EQ.0) go to 330

!DIR$ IVDEP
!OCL NOVREC
     DO 320 jl=1,kproma
        IF(ldflag(jl) .AND. ABS(zcond(jl)) > 0._dp) THEN
           it = NINT(pt(jl,kk)*1000._dp)
           IF (it<jptlucu1 .OR. it>jptlucu2) lookupoverflow = .TRUE.
           it = MAX(MIN(it,jptlucu2),jptlucu1)
!--- Included for prognostic CDNC/IC scheme (Ulrike Lohmann, 11/02/2007)
           IF (ncvmicro>0) THEN
              lo2=(pt(jl,kk).LE.cthomi).OR.                              &
                   ((pt(jl,kk).LT.tmelt).AND.(plui(jl,kk).GT.csecfrl))
              zes=MERGE(tlucua(it)/pp(jl),tlucuaw(it)/pp(jl),lo2)
           ELSE
              zes=tlucua(it)/pp(jl)
           ENDIF
!-------------------------------end----------------------------------
           zes=MIN(0.5_dp,zes)
           LO=zes<0.4_dp
           zcor=1._dp/(1._dp-vtmpc1*zes)
           zqsat=zes*zcor
           it1=it+1
           it1 = MAX(MIN(it1,jptlucu2),jptlucu1)
           IF (ncvmicro>0) THEN
              lo2=(pt(jl,kk).LE.cthomi).OR.                              &
                   ((pt(jl,kk).LT.tmelt).AND.(plui(jl,kk).GT.csecfrl))
              zqst1=MERGE(tlucua(it1)/pp(jl),tlucuaw(it1)/pp(jl),lo2)
              zqst1=MIN(0.5_dp,zqst1)
              zqst1=zqst1/(1._dp-vtmpc1*zqst1)
              zdqsdt=(zqst1-zqsat)*1000._dp
              ztlucub=MERGE(tlucub(it),tlucubw(it),lo2)
              ztlucuc=MERGE(tlucuc(it),tlucucw(it),lo2)
              zlcdqsdt=MERGE(zdqsdt*ztlucuc,zqsat*zcor*ztlucub,LO)
              zcond1=(pq(jl,kk)-zqsat)/(1._dp+zlcdqsdt)
              pt(jl,kk)=pt(jl,kk)+ztlucuc*zcond1
           ELSE
              zqst1=tlucua(it1)/pp(jl)
              zqst1=MIN(0.5_dp,zqst1)
              zqst1=zqst1/(1._dp-vtmpc1*zqst1)
              zdqsdt=(zqst1-zqsat)*1000._dp
              zlcdqsdt=MERGE(zdqsdt*tlucuc(it),zqsat*zcor*tlucub(it),LO)
              zcond1=(pq(jl,kk)-zqsat)/(1._dp+zlcdqsdt)
              pt(jl,kk)=pt(jl,kk)+tlucuc(it)*zcond1
           ENDIF
           pq(jl,kk)=pq(jl,kk)-zcond1
        END IF
320  END DO

     IF (lookupoverflow) CALL lookuperror ('cuadjtqi (4) ')

330  CONTINUE

  END IF

  IF(kcall.EQ.0) THEN

     isum=0
!DIR$ IVDEP
!OCL NOVREC
     DO 410 jl=1,kproma
        it = NINT(pt(jl,kk)*1000._dp)
        IF (it<jptlucu1 .OR. it>jptlucu2) lookupoverflow = .TRUE.
        it = MAX(MIN(it,jptlucu2),jptlucu1)
!--- Included for prognostic CDNC/IC scheme (Ulrike Lohmann, 11/02/2007)
           IF (ncvmicro>0) THEN
              lo2=(pt(jl,kk).LE.cthomi).OR.                              &
                   ((pt(jl,kk).LT.tmelt).AND.(plui(jl,kk).GT.csecfrl))
              zes=MERGE(tlucua(it)/pp(jl),tlucuaw(it)/pp(jl),lo2)
           ELSE
              zes=tlucua(it)/pp(jl)
           ENDIF
!-------------------------------end----------------------------------
        zes=MIN(0.5_dp,zes)
        LO=zes<0.4_dp
        zcor=1._dp/(1._dp-vtmpc1*zes)
        zqsat=zes*zcor
        it1=it+1
        it1 = MAX(MIN(it1,jptlucu2),jptlucu1)
           IF (ncvmicro>0) THEN
              lo2=(pt(jl,kk).LE.cthomi).OR.                              &
                   ((pt(jl,kk).LT.tmelt).AND.(plui(jl,kk).GT.csecfrl))
              zqst1=MERGE(tlucua(it1)/pp(jl),tlucuaw(it1)/pp(jl),lo2)
              zqst1=MIN(0.5_dp,zqst1)
              zqst1=zqst1/(1._dp-vtmpc1*zqst1)
              zdqsdt=(zqst1-zqsat)*1000._dp
              ztlucub=MERGE(tlucub(it),tlucubw(it),lo2)
              ztlucuc=MERGE(tlucuc(it),tlucucw(it),lo2)
              zlcdqsdt=MERGE(zdqsdt*ztlucuc,zqsat*zcor*ztlucub,LO)
              zcond(jl)=(pq(jl,kk)-zqsat)/(1._dp+zlcdqsdt)
              pt(jl,kk)=pt(jl,kk)+ztlucuc*zcond(jl)
           ELSE
              zqst1=tlucua(it1)/pp(jl)
              zqst1=MIN(0.5_dp,zqst1)
              zqst1=zqst1/(1._dp-vtmpc1*zqst1)
              zdqsdt=(zqst1-zqsat)*1000._dp
              zlcdqsdt=MERGE(zdqsdt*tlucuc(it),zqsat*zcor*tlucub(it),LO)
              zcond(jl)=(pq(jl,kk)-zqsat)/(1._dp+zlcdqsdt)
              pt(jl,kk)=pt(jl,kk)+tlucuc(it)*zcond(jl)
           ENDIF
        pq(jl,kk)=pq(jl,kk)-zcond(jl)
        IF(ABS(zcond(jl)) > 0._dp) isum=isum+1
410  END DO

     IF (lookupoverflow) CALL lookuperror ('cuadjtqi (5) ')

     IF(isum.EQ.0) go to 430

!DIR$ IVDEP
!OCL NOVREC
     DO 420 jl=1,kproma
        it = NINT(pt(jl,kk)*1000._dp)
        IF (it<jptlucu1 .OR. it>jptlucu2) lookupoverflow = .TRUE.
        it = MAX(MIN(it,jptlucu2),jptlucu1)
!--- Included for prognostic CDNC/IC scheme (Ulrike Lohmann, 11/02/2007)
           IF (ncvmicro>0) THEN
              lo2=(pt(jl,kk).LE.cthomi).OR.                              &
                   ((pt(jl,kk).LT.tmelt).AND.(plui(jl,kk).GT.csecfrl))
              zes=MERGE(tlucua(it)/pp(jl),tlucuaw(it)/pp(jl),lo2)
           ELSE
              zes=tlucua(it)/pp(jl)
           ENDIF
!-------------------------------end----------------------------------
        zes=MIN(0.5_dp,zes)
        LO=zes<0.4_dp
        zcor=1._dp/(1._dp-vtmpc1*zes)
        zqsat=zes*zcor
        it1=it+1
        it1 = MAX(MIN(it1,jptlucu2),jptlucu1)
           IF (ncvmicro>0) THEN
              lo2=(pt(jl,kk).LE.cthomi).OR.                              &
                   ((pt(jl,kk).LT.tmelt).AND.(plui(jl,kk).GT.csecfrl))
              zqst1=MERGE(tlucua(it1)/pp(jl),tlucuaw(it1)/pp(jl),lo2)
              zqst1=MIN(0.5_dp,zqst1)
              zqst1=zqst1/(1._dp-vtmpc1*zqst1)
              zdqsdt=(zqst1-zqsat)*1000._dp
              ztlucub=MERGE(tlucub(it),tlucubw(it),lo2)
              ztlucuc=MERGE(tlucuc(it),tlucucw(it),lo2)
              zlcdqsdt=MERGE(zdqsdt*ztlucuc,zqsat*zcor*ztlucub,LO)
              zcond1=(pq(jl,kk)-zqsat)/(1._dp+zlcdqsdt)
              pt(jl,kk)=pt(jl,kk)+ztlucuc*zcond1
           ELSE
              zqst1=tlucua(it1)/pp(jl)
              zqst1=MIN(0.5_dp,zqst1)
              zqst1=zqst1/(1._dp-vtmpc1*zqst1)
              zdqsdt=(zqst1-zqsat)*1000._dp
              zlcdqsdt=MERGE(zdqsdt*tlucuc(it),zqsat*zcor*tlucub(it),LO)
              zcond1=(pq(jl,kk)-zqsat)/(1._dp+zlcdqsdt)
              pt(jl,kk)=pt(jl,kk)+tlucuc(it)*zcond1
           ENDIF
        pq(jl,kk)=pq(jl,kk)-zcond1
420  END DO

     IF (lookupoverflow) CALL lookuperror ('cuadjtqi (6) ')

430  CONTINUE

  END IF

  IF(kcall.EQ.4) THEN

!DIR$ IVDEP
!OCL NOVREC
     DO 510 jl=1,kproma
        it = NINT(pt(jl,kk)*1000._dp)
        IF (it<jptlucu1 .OR. it>jptlucu2) lookupoverflow = .TRUE.
        it = MAX(MIN(it,jptlucu2),jptlucu1)
!--- Included for prognostic CDNC/IC scheme (Ulrike Lohmann, 11/02/2007)
           IF (ncvmicro>0) THEN
              lo2=(pt(jl,kk).LE.cthomi).OR.                              &
                   ((pt(jl,kk).LT.tmelt).AND.(plui(jl,kk).GT.csecfrl))
              zes=MERGE(tlucua(it)/pp(jl),tlucuaw(it)/pp(jl),lo2)
           ELSE
              zes=tlucua(it)/pp(jl)
           ENDIF
!-------------------------------end----------------------------------
        zes=MIN(0.5_dp,zes)
        LO=zes<0.4_dp
        zcor=1._dp/(1._dp-vtmpc1*zes)
        zqsat=zes*zcor
        it1=it+1
        it1 = MAX(MIN(it1,jptlucu2),jptlucu1)
           IF (ncvmicro>0) THEN
              lo2=(pt(jl,kk).LE.cthomi).OR.                              &
                   ((pt(jl,kk).LT.tmelt).AND.(plui(jl,kk).GT.csecfrl))
              zqst1=MERGE(tlucua(it1)/pp(jl),tlucuaw(it1)/pp(jl),lo2)
              zqst1=MIN(0.5_dp,zqst1)
              zqst1=zqst1/(1._dp-vtmpc1*zqst1)
              zdqsdt=(zqst1-zqsat)*1000._dp
              ztlucub=MERGE(tlucub(it),tlucubw(it),lo2)
              ztlucuc=MERGE(tlucuc(it),tlucucw(it),lo2)
              zlcdqsdt=MERGE(zdqsdt*ztlucuc,zqsat*zcor*ztlucub,LO)
              zcond(jl)=(pq(jl,kk)-zqsat)/(1._dp+zlcdqsdt)
              pt(jl,kk)=pt(jl,kk)+ztlucuc*zcond(jl)
           ELSE
              zqst1=tlucua(it1)/pp(jl)
              zqst1=MIN(0.5_dp,zqst1)
              zqst1=zqst1/(1._dp-vtmpc1*zqst1)
              zdqsdt=(zqst1-zqsat)*1000._dp
              zlcdqsdt=MERGE(zdqsdt*tlucuc(it),zqsat*zcor*tlucub(it),LO)
              zcond(jl)=(pq(jl,kk)-zqsat)/(1._dp+zlcdqsdt)
              pt(jl,kk)=pt(jl,kk)+tlucuc(it)*zcond(jl)
           ENDIF
        pq(jl,kk)=pq(jl,kk)-zcond(jl)
510  END DO

     IF (lookupoverflow) CALL lookuperror ('cuadjtqi (7) ')

!DIR$ IVDEP
!OCL NOVREC
     DO 520 jl=1,kproma
        it = NINT(pt(jl,kk)*1000._dp)
        IF (it<jptlucu1 .OR. it>jptlucu2) lookupoverflow = .TRUE.
        it = MAX(MIN(it,jptlucu2),jptlucu1)
!--- Included for prognostic CDNC/IC scheme (Ulrike Lohmann, 11/02/2007)
           IF (ncvmicro>0) THEN
              lo2=(pt(jl,kk).LE.cthomi).OR.                              &
                   ((pt(jl,kk).LT.tmelt).AND.(plui(jl,kk).GT.csecfrl))
              zes=MERGE(tlucua(it)/pp(jl),tlucuaw(it)/pp(jl),lo2)
           ELSE
              zes=tlucua(it)/pp(jl)
           ENDIF
!-------------------------------end----------------------------------
        zes=MIN(0.5_dp,zes)
        LO=zes<0.4_dp
        zcor=1._dp/(1._dp-vtmpc1*zes)
        zqsat=zes*zcor
        it1=it+1
        it1 = MAX(MIN(it1,jptlucu2),jptlucu1)
           IF (ncvmicro>0) THEN
              lo2=(pt(jl,kk).LE.cthomi).OR.                              &
                   ((pt(jl,kk).LT.tmelt).AND.(plui(jl,kk).GT.csecfrl))
              zqst1=MERGE(tlucua(it1)/pp(jl),tlucuaw(it1)/pp(jl),lo2)
              zqst1=MIN(0.5_dp,zqst1)
              zqst1=zqst1/(1._dp-vtmpc1*zqst1)
              zdqsdt=(zqst1-zqsat)*1000._dp
              ztlucub=MERGE(tlucub(it),tlucubw(it),lo2)
              ztlucuc=MERGE(tlucuc(it),tlucucw(it),lo2)
              zlcdqsdt=MERGE(zdqsdt*ztlucuc,zqsat*zcor*ztlucub,LO)
              zcond1=(pq(jl,kk)-zqsat)/(1._dp+zlcdqsdt)
              pt(jl,kk)=pt(jl,kk)+ztlucuc*zcond1
           ELSE
              zqst1=tlucua(it1)/pp(jl)
              zqst1=MIN(0.5_dp,zqst1)
              zqst1=zqst1/(1._dp-vtmpc1*zqst1)
              zdqsdt=(zqst1-zqsat)*1000._dp
              zlcdqsdt=MERGE(zdqsdt*tlucuc(it),zqsat*zcor*tlucub(it),LO)
              zcond1=(pq(jl,kk)-zqsat)/(1._dp+zlcdqsdt)
              pt(jl,kk)=pt(jl,kk)+tlucuc(it)*zcond1
           ENDIF
        pq(jl,kk)=pq(jl,kk)-zcond1
520  END DO

     IF (lookupoverflow) CALL lookuperror ('cuadjtqi (8) ')

  END IF

    RETURN
  END SUBROUTINE cuadjtqi

END MODULE mo_cuadjtqi
