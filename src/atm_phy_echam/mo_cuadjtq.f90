!>
!! Module contains subroutines cuadjtq and cuadjtq_idx.
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
#ifndef __xlC__
#define SWDIV_NOCHK(a,b) (a)/(b)
#define FSEL(a,b,c) MERGE(b,c,(a) >= 0._dp)
#endif

MODULE mo_cuadjtq

  USE mo_kind,               ONLY: dp

#ifdef __ICON__
  USE mo_physical_constants, ONLY: vtmpc1
#else
  USE mo_constants,          ONLY: vtmpc1
#endif

  USE mo_convect_tables,     ONLY: tlucua,   & ! table a
                                   tlucub,   & ! table b
                                   tlucuc,   & ! table c
                                   jptlucu1, jptlucu2, &
                                   lookuperror
#ifdef __ibmspline__
  USE mo_convect_tables, ONLY: prepare_ua_index_spline, lookup_ua_spline, &
    &                          lookup_ua_list_spline, lookup_ubc,         &
    &                          lookup_ubc_list
#else
  USE mo_convect_tables, ONLY: prepare_ua_index, lookup_ua, lookup_ua_list, &
    &                          lookup_ubc, lookup_ubc_list
#endif
#ifdef _PROFILE
  USE mo_profile,        ONLY: trace_start, trace_stop
#endif

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: cuadjtq
  PUBLIC :: cuadjtq_idx

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

CONTAINS
  !>
  !!
SUBROUTINE cuadjtq_idx(  kproma, kbdim, klev, kk,             &
           pp,       pt,       pq,       ldidx, ldcnt,  kcall)
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
  INTEGER, INTENT (IN) :: kcall, kk, klev, kproma, kbdim

  !  Array arguments with intent(In):
  REAL(dp), INTENT (IN) :: pp(kbdim)
  INTEGER, INTENT (IN) :: ldidx(kbdim)
  INTEGER, INTENT (IN) :: ldcnt

  !  Array arguments with intent(InOut):
  REAL(dp), INTENT (INOUT) :: pq(kbdim,klev), pt(kbdim,klev)

  !  Local scalars:
  REAL(dp):: zcond1, zqst1, zdqsdt, zqsat, zes, zcor, zlcdqsdt
  INTEGER :: jl, nl, nsum

  !  Local arrays:
  REAL(dp) :: zcond(kbdim),zppi(kbdim)
#ifdef __ibmspline__
  REAL(dp) :: za(kbdim)
#endif
  REAL(dp) :: ua(kbdim),dua(kbdim),ub(kbdim),uc(kbdim)
  INTEGER  :: idx(kbdim),ncond(kbdim)


  !  Executable statements

#ifdef _PROFILE
  CALL trace_start ('cuadjtq_idx', 32)
#endif
!
!----------------------------------------------------------------------
!
!     2.           CALCULATE CONDENSATION AND ADJUST T AND Q ACCORDINGLY
!                  -----------------------------------------------------
!
  IF (kcall >= 0.AND. kcall <= 2 ) THEN

     CALL lookup_ubc_list(kproma,ldcnt,ldidx(1),pt(1,kk),ub(1),uc(1))
#ifdef __ibmspline__
     CALL lookup_ua_list_spline('cuadjtq_idx (1)',ldcnt,ldidx(1),pt(1,kk),ua(1),dua(1))
#else
     CALL lookup_ua_list('cuadjtq_idx (1)',kproma,ldcnt,ldidx(1),pt(1,kk),ua(1),dua(1))
#endif

!DIR$ IVDEP
!OCL NOVREC
!CDIR NODEP
!IBM* ASSERT(NODEPS)
     DO 111 nl=1,ldcnt
        jl = ldidx(nl)
        zppi(jl)=1._dp/pp(jl)
111  END DO

     IF (kcall.EQ.0) THEN

        ! mpuetz: 40% of cuadtjq (lookup of tlucua und tlucub !!)

!DIR$ IVDEP
!OCL NOVREC
!CDIR NODEP
!IBM ASSERT(NODEPS)
        DO 112 nl=1,ldcnt
           jl = ldidx(nl)

           zes  = ua(nl)*zppi(jl)
           zes  = MIN(0.5_dp,zes)
           zcor = SWDIV_NOCHK(1._dp,(1._dp-vtmpc1*zes))
           zqsat= zes*zcor

#ifdef __ibmspline2__
           zdqsdt = zppi(jl)*zcor**2*dua(nl)
#else
           zqst1=(ua(nl)+0.001_dp*dua(nl))*zppi(jl)
           zqst1=MIN(0.5_dp,zqst1)
           zqst1=SWDIV_NOCHK(zqst1,(1._dp-vtmpc1*zqst1))
           zdqsdt=(zqst1 - zqsat)*1000._dp
#endif
           zlcdqsdt  = FSEL(zes-0.4_dp,zqsat*zcor*ub(nl),zdqsdt*uc(nl))

           zcond(jl) = SWDIV_NOCHK((pq(jl,kk)-zqsat),(1._dp+zlcdqsdt))

           pt(jl,kk) = pt(jl,kk) + uc(nl)*zcond(jl)
           pq(jl,kk) = pq(jl,kk) - zcond(jl)

           ! mpuetz: zcond is almost always > 0
           ncond(nl) = INT(FSEL(-ABS(zcond(jl)),0._dp,1._dp))
112     END DO

     ELSE IF (kcall.EQ.1) THEN

!DIR$ IVDEP
!OCL NOVREC
!CDIR NODEP
!IBM ASSERT(NODEPS)
        DO 212 nl=1,ldcnt
           jl = ldidx(nl)

           zes  = ua(nl)*zppi(jl)
           zes  = MIN(0.5_dp,zes)
           zcor = SWDIV_NOCHK(1._dp,(1._dp-vtmpc1*zes))
           zqsat= zes*zcor

#ifdef __ibmspline2__
           zdqsdt = zppi(jl)*zcor**2*dua(nl)
#else
           zqst1=(ua(nl)+0.001_dp*dua(nl))*zppi(jl)
           zqst1=MIN(0.5_dp,zqst1)
           zqst1=SWDIV_NOCHK(zqst1,(1._dp-vtmpc1*zqst1))
           zdqsdt=(zqst1 - zqsat)*1000._dp
#endif
           zlcdqsdt  = FSEL(zes-0.4_dp,zqsat*zcor*ub(nl),zdqsdt*uc(nl))

           zcond(jl) = SWDIV_NOCHK((pq(jl,kk)-zqsat),(1._dp+zlcdqsdt))
           zcond(jl) = MAX(zcond(jl),0._dp)

           pt(jl,kk) = pt(jl,kk) + uc(nl)*zcond(jl)
           pq(jl,kk) = pq(jl,kk) - zcond(jl)
           ! mpuetz: zcond is almost always > 0
           ncond(nl) = INT(FSEL(-ABS(zcond(jl)),0._dp,1._dp))
212     END DO

     ELSE

!DIR$ IVDEP
!OCL NOVREC
!CDIR NODEP
!IBM* ASSERT(NODEPS)
        DO 312 nl=1,ldcnt
           jl = ldidx(nl)

           zes  = ua(nl)*zppi(jl)
           zes  = MIN(0.5_dp,zes)
           zcor = SWDIV_NOCHK(1._dp,(1._dp-vtmpc1*zes))
           zqsat= zes*zcor

#ifdef __ibmspline2__
           zdqsdt = zppi(jl)*zcor**2*dua(nl)
#else
           zqst1=(ua(nl)+0.001_dp*dua(nl))*zppi(jl)
           zqst1=MIN(0.5_dp,zqst1)
           zqst1=SWDIV_NOCHK(zqst1,(1._dp-vtmpc1*zqst1))
           zdqsdt=(zqst1 - zqsat)*1000._dp
#endif
           zlcdqsdt  = FSEL(zes-0.4_dp,zqsat*zcor*ub(nl),zdqsdt*uc(nl))
           zcond(jl) = SWDIV_NOCHK((pq(jl,kk)-zqsat),(1._dp+zlcdqsdt))

           zcond(jl) = MIN(zcond(jl),0._dp)

           pt(jl,kk) = pt(jl,kk) + uc(nl)*zcond(jl)
           pq(jl,kk) = pq(jl,kk) - zcond(jl)
           ! mpuetz: zcond is almost always > 0
           ncond(nl) = INT(FSEL(-ABS(zcond(jl)),0._dp,1._dp))
312     END DO

     END IF

     nsum = 1
     DO nl=1,ldcnt
        idx(nsum) = ldidx(nl)
        nsum = nsum + ncond(nl)
     END DO
     nsum = nsum - 1

#ifdef __ibmdbg__
     print *,'cuad(',kcall,')',ldcnt,nsum,ldcnt
#endif

     IF (nsum > 0) THEN

        CALL lookup_ubc_list(kproma,nsum,idx(1),pt(1,kk),ub(1),uc(1))
#ifdef __ibmspline__
        CALL lookup_ua_list_spline('cuadjtq_idx (2)',nsum,idx(1),pt(1,kk),ua(1),dua(1))
#else
        CALL lookup_ua_list('cuadjtq_idx (2)',kproma,nsum,idx(1),pt(1,kk),ua(1),dua(1))
#endif

!DIR$ IVDEP
!OCL NOVREC
!CDIR NODEP
!IBM* ASSERT(NODEPS)
!IBM* UNROLL(3)
        DO 116 nl=1,nsum
           jl = idx(nl)

           zes  = ua(nl)*zppi(jl)
           zes  = MIN(0.5_dp,zes)
           zcor = SWDIV_NOCHK(1._dp,(1._dp-vtmpc1*zes))
           zqsat= zes*zcor

#ifdef __ibmspline2__
           zdqsdt = zppi(jl)*zcor**2*dua(nl)
#else
           zqst1=(ua(nl)+0.001_dp*dua(nl))*zppi(jl)
           zqst1=MIN(0.5_dp,zqst1)
           zqst1=SWDIV_NOCHK(zqst1,(1._dp-vtmpc1*zqst1))
           zdqsdt=(zqst1 - zqsat)*1000._dp
#endif
           zlcdqsdt = FSEL(zes-0.4_dp,zqsat*zcor*ub(nl),zdqsdt*uc(nl))
           zcond1   = SWDIV_NOCHK((pq(jl,kk)-zqsat),(1._dp+zlcdqsdt))

           pt(jl,kk) = pt(jl,kk) + uc(nl)*zcond1
           pq(jl,kk) = pq(jl,kk) - zcond1
116     END DO

     END IF

  ELSE IF (kcall.EQ.4) THEN

     ! no indirect addressing

     zppi(1:kproma)=1._dp/pp(jl)

#ifdef __ibmspline__
     CALL prepare_ua_index_spline('cuadjtq_idx (3)',kproma,pt(1,kk),idx(1),za(1))
     CALL lookup_ua_spline(kproma,idx(1),za(1),ua(1),dua(1))
#else
     CALL prepare_ua_index('cuadjtq_idx (3)',kproma,pt(1,kk),idx(1))
     CALL lookup_ua(kproma,idx(1),ua(1),dua(1))
#endif
     CALL lookup_ubc(kproma,pt(1,kk),ub(1),uc(1))

     DO 412 jl=1,kproma

        zes  = ua(jl)*zppi(jl)
        zes  = MIN(0.5_dp,zes)
        zcor = SWDIV_NOCHK(1._dp,(1._dp-vtmpc1*zes))
        zqsat= zes*zcor

#ifdef __ibmspline2__
           zdqsdt = zppi(jl)*zcor**2*dua(nl)
#else
           zqst1=(ua(nl)+0.001_dp*dua(nl))*zppi(jl)
           zqst1=MIN(0.5_dp,zqst1)
           zqst1=SWDIV_NOCHK(zqst1,(1._dp-vtmpc1*zqst1))
           zdqsdt=(zqst1 - zqsat)*1000._dp
#endif
        zlcdqsdt  = FSEL(zes-0.4_dp,zqsat*zcor*ub(jl),zdqsdt*uc(jl))
        zcond(jl) = SWDIV_NOCHK((pq(jl,kk)-zqsat),(1._dp+zlcdqsdt))

        pt(jl,kk) = pt(jl,kk) + uc(jl)*zcond(jl)
        pq(jl,kk) = pq(jl,kk) - zcond(jl)
412  END DO

#ifdef __ibmspline__
     CALL prepare_ua_index_spline('cuadjtq_idx (4)',kproma,pt(1,kk),idx(1),za(1))
     CALL lookup_ua_spline(kproma,idx(1),za(1),ua(1),dua(1))
#else
     CALL prepare_ua_index('cuadjtq_idx (4)',kproma,pt(1,kk),idx(1))
     CALL lookup_ua(kproma,idx(1),ua(1),dua(1))
#endif
     CALL lookup_ubc(kproma,pt(1,kk),ub(1),uc(1))

     DO 413 jl=1,kproma

        zes  = ua(jl)*zppi(jl)
        zes  = MIN(0.5_dp,zes)
        zcor = SWDIV_NOCHK(1._dp,(1._dp-vtmpc1*zes))
        zqsat= zes*zcor

#ifdef __ibmspline2__
           zdqsdt = zppi(jl)*zcor**2*dua(nl)
#else
           zqst1=(ua(nl)+0.001_dp*dua(nl))*zppi(jl)
           zqst1=MIN(0.5_dp,zqst1)
           zqst1=SWDIV_NOCHK(zqst1,(1._dp-vtmpc1*zqst1))
           zdqsdt=(zqst1 - zqsat)*1000._dp
#endif
        zlcdqsdt = FSEL(zes-0.4_dp,zqsat*zcor*ub(jl),zdqsdt*uc(jl))
        zcond1   = SWDIV_NOCHK((pq(jl,kk)-zqsat),(1._dp+zlcdqsdt))

        pt(jl,kk)=pt(jl,kk) + uc(jl)*zcond1
        pq(jl,kk)=pq(jl,kk) - zcond1
413  END DO

  END IF

#ifdef _PROFILE
  CALL trace_stop ('cuadjtq_idx', 32)
#endif

    RETURN
  END SUBROUTINE cuadjtq_idx
  !-------------
  !>
  !!
  SUBROUTINE cuadjtq(  kproma, kbdim, klev, kk,             &
             pp,       pt,       pq,       ldflag,   kcall)
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
  INTEGER, INTENT (IN) :: kcall, kk, klev, kproma, kbdim

  !  Array arguments with intent(In):
  REAL(dp), INTENT (IN) :: pp(kbdim)
  LOGICAL, INTENT (IN) :: ldflag(kbdim)

  !  Array arguments with intent(InOut):
  REAL(dp), INTENT (INOUT) :: pq(kbdim,klev), pt(kbdim,klev)

  !  Local scalars:
  REAL(dp):: zcond1, zqst1, zdqsdt, zqsat, zes, zcor, zlcdqsdt
  INTEGER :: isum, jl, it, it1
  LOGICAL :: LO, lookupoverflow

  !  Local arrays:
  REAL(dp):: zcond(kbdim)


  !  Executable statements
!
#ifdef _PROFILE
  CALL trace_start ('cuadjtq', 30)
#endif
!

  lookupoverflow = .FALSE.

  zcond = 0.0_dp
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
           zes=tlucua(it)/pp(jl)
           zes=MIN(0.5_dp,zes)
           LO=zes<0.4_dp
           zcor=1._dp/(1._dp-vtmpc1*zes)
           zqsat=zes*zcor
           it1=it+1
           it1 = MAX(MIN(it1,jptlucu2),jptlucu1)
           zqst1=tlucua(it1)/pp(jl)
           zqst1=MIN(0.5_dp,zqst1)
           zqst1=zqst1/(1._dp-vtmpc1*zqst1)
           zdqsdt=(zqst1-zqsat)*1000._dp
           zlcdqsdt=MERGE(zdqsdt*tlucuc(it),zqsat*zcor*tlucub(it),LO)
           zcond(jl)=(pq(jl,kk)-zqsat)/(1._dp+zlcdqsdt)
           zcond(jl)=MAX(zcond(jl),0._dp)
           pt(jl,kk)=pt(jl,kk)+tlucuc(it)*zcond(jl)
           pq(jl,kk)=pq(jl,kk)-zcond(jl)
           IF(ABS(zcond(jl)) > 0._dp) isum=isum+1
        END IF
210  END DO

     IF (lookupoverflow) CALL lookuperror ('cuadjtq (1) ')

     IF(isum.EQ.0) go to 230

!DIR$ IVDEP
!OCL NOVREC
     DO 220 jl=1,kproma
        IF(ldflag(jl)) THEN
           IF(ABS(zcond(jl)) > 0._dp) THEN
           it = NINT(pt(jl,kk)*1000._dp)
           IF (it<jptlucu1 .OR. it>jptlucu2) lookupoverflow = .TRUE.
           it = MAX(MIN(it,jptlucu2),jptlucu1)
           zes=tlucua(it)/pp(jl)
           zes=MIN(0.5_dp,zes)
           LO=zes<0.4_dp
           zcor=1._dp/(1._dp-vtmpc1*zes)
           zqsat=zes*zcor
           it1=it+1
           it1 = MAX(MIN(it1,jptlucu2),jptlucu1)
           zqst1=tlucua(it1)/pp(jl)
           zqst1=MIN(0.5_dp,zqst1)
           zqst1=zqst1/(1._dp-vtmpc1*zqst1)
           zdqsdt=(zqst1-zqsat)*1000._dp
           zlcdqsdt=MERGE(zdqsdt*tlucuc(it),zqsat*zcor*tlucub(it),LO)
           zcond1=(pq(jl,kk)-zqsat)/(1._dp+zlcdqsdt)
           pt(jl,kk)=pt(jl,kk)+tlucuc(it)*zcond1
           pq(jl,kk)=pq(jl,kk)-zcond1
           END IF
        END IF
220  END DO

     IF (lookupoverflow) CALL lookuperror ('cuadjtq (2) ')

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
           zes=tlucua(it)/pp(jl)
           zes=MIN(0.5_dp,zes)
           LO=zes<0.4_dp
           zcor=1._dp/(1._dp-vtmpc1*zes)
           zqsat=zes*zcor
           it1=it+1
           it1 = MAX(MIN(it1,jptlucu2),jptlucu1)
           zqst1=tlucua(it1)/pp(jl)
           zqst1=MIN(0.5_dp,zqst1)
           zqst1=zqst1/(1._dp-vtmpc1*zqst1)
           zdqsdt=(zqst1-zqsat)*1000._dp
           zlcdqsdt=MERGE(zdqsdt*tlucuc(it),zqsat*zcor*tlucub(it),LO)
           zcond(jl)=(pq(jl,kk)-zqsat)/(1._dp+zlcdqsdt)
           zcond(jl)=MIN(zcond(jl),0._dp)
           pt(jl,kk)=pt(jl,kk)+tlucuc(it)*zcond(jl)
           pq(jl,kk)=pq(jl,kk)-zcond(jl)
           IF(ABS(zcond(jl)) > 0._dp) isum=isum+1
        END IF
310  END DO

     IF (lookupoverflow) CALL lookuperror ('cuadjtq (3) ')

     IF(isum.EQ.0) go to 330

!DIR$ IVDEP
!OCL NOVREC
     DO 320 jl=1,kproma
        IF(ldflag(jl) .AND. ABS(zcond(jl)) > 0._dp) THEN
           it = NINT(pt(jl,kk)*1000._dp)
           IF (it<jptlucu1 .OR. it>jptlucu2) lookupoverflow = .TRUE.
           it = MAX(MIN(it,jptlucu2),jptlucu1)
           zes=tlucua(it)/pp(jl)
           zes=MIN(0.5_dp,zes)
           LO=zes<0.4_dp
           zcor=1._dp/(1._dp-vtmpc1*zes)
           zqsat=zes*zcor
           it1=it+1
           it1 = MAX(MIN(it1,jptlucu2),jptlucu1)
           zqst1=tlucua(it1)/pp(jl)
           zqst1=MIN(0.5_dp,zqst1)
           zqst1=zqst1/(1._dp-vtmpc1*zqst1)
           zdqsdt=(zqst1-zqsat)*1000._dp
           zlcdqsdt=MERGE(zdqsdt*tlucuc(it),zqsat*zcor*tlucub(it),LO)
           zcond1=(pq(jl,kk)-zqsat)/(1._dp+zlcdqsdt)
           pt(jl,kk)=pt(jl,kk)+tlucuc(it)*zcond1
           pq(jl,kk)=pq(jl,kk)-zcond1
        END IF
320  END DO

     IF (lookupoverflow) CALL lookuperror ('cuadjtq (4) ')

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
        zes=tlucua(it)/pp(jl)
        zes=MIN(0.5_dp,zes)
        LO=zes<0.4_dp
        zcor=1._dp/(1._dp-vtmpc1*zes)
        zqsat=zes*zcor
        it1=it+1
        it1 = MAX(MIN(it1,jptlucu2),jptlucu1)
        zqst1=tlucua(it1)/pp(jl)
        zqst1=MIN(0.5_dp,zqst1)
        zqst1=zqst1/(1._dp-vtmpc1*zqst1)
        zdqsdt=(zqst1-zqsat)*1000._dp
        zlcdqsdt=MERGE(zdqsdt*tlucuc(it),zqsat*zcor*tlucub(it),LO)
        zcond(jl)=(pq(jl,kk)-zqsat)/(1._dp+zlcdqsdt)
        pt(jl,kk)=pt(jl,kk)+tlucuc(it)*zcond(jl)
        pq(jl,kk)=pq(jl,kk)-zcond(jl)
        IF(ABS(zcond(jl)) > 0._dp) isum=isum+1
410  END DO

     IF (lookupoverflow) CALL lookuperror ('cuadjtq (5) ')

     IF(isum.EQ.0) go to 430

!DIR$ IVDEP
!OCL NOVREC
     DO 420 jl=1,kproma
        it = NINT(pt(jl,kk)*1000._dp)
        IF (it<jptlucu1 .OR. it>jptlucu2) lookupoverflow = .TRUE.
        it = MAX(MIN(it,jptlucu2),jptlucu1)
        zes=tlucua(it)/pp(jl)
        zes=MIN(0.5_dp,zes)
        LO=zes<0.4_dp
        zcor=1._dp/(1._dp-vtmpc1*zes)
        zqsat=zes*zcor
        it1=it+1
        it1 = MAX(MIN(it1,jptlucu2),jptlucu1)
        zqst1=tlucua(it1)/pp(jl)
        zqst1=MIN(0.5_dp,zqst1)
        zqst1=zqst1/(1._dp-vtmpc1*zqst1)
        zdqsdt=(zqst1-zqsat)*1000._dp
        zlcdqsdt=MERGE(zdqsdt*tlucuc(it),zqsat*zcor*tlucub(it),LO)
        zcond1=(pq(jl,kk)-zqsat)/(1._dp+zlcdqsdt)
        pt(jl,kk)=pt(jl,kk)+tlucuc(it)*zcond1
        pq(jl,kk)=pq(jl,kk)-zcond1
420  END DO

     IF (lookupoverflow) CALL lookuperror ('cuadjtq (6) ')

430  CONTINUE

  END IF

  IF(kcall.EQ.4) THEN

!DIR$ IVDEP
!OCL NOVREC
     DO 510 jl=1,kproma
        it = NINT(pt(jl,kk)*1000._dp)
        IF (it<jptlucu1 .OR. it>jptlucu2) lookupoverflow = .TRUE.
        it = MAX(MIN(it,jptlucu2),jptlucu1)
        zes=tlucua(it)/pp(jl)
        zes=MIN(0.5_dp,zes)
        LO=zes<0.4_dp
        zcor=1._dp/(1._dp-vtmpc1*zes)
        zqsat=zes*zcor
        it1=it+1
        it1 = MAX(MIN(it1,jptlucu2),jptlucu1)
        zqst1=tlucua(it1)/pp(jl)
        zqst1=MIN(0.5_dp,zqst1)
        zqst1=zqst1/(1._dp-vtmpc1*zqst1)
        zdqsdt=(zqst1-zqsat)*1000._dp
        zlcdqsdt=MERGE(zdqsdt*tlucuc(it),zqsat*zcor*tlucub(it),LO)
        zcond(jl)=(pq(jl,kk)-zqsat)/(1._dp+zlcdqsdt)
        pt(jl,kk)=pt(jl,kk)+tlucuc(it)*zcond(jl)
        pq(jl,kk)=pq(jl,kk)-zcond(jl)
510  END DO

     IF (lookupoverflow) CALL lookuperror ('cuadjtq (7) ')

!DIR$ IVDEP
!OCL NOVREC
     DO 520 jl=1,kproma
        it = NINT(pt(jl,kk)*1000._dp)
        IF (it<jptlucu1 .OR. it>jptlucu2) lookupoverflow = .TRUE.
        it = MAX(MIN(it,jptlucu2),jptlucu1)
        zes=tlucua(it)/pp(jl)
        zes=MIN(0.5_dp,zes)
        LO=zes<0.4_dp
        zcor=1._dp/(1._dp-vtmpc1*zes)
        zqsat=zes*zcor
        it1=it+1
        it1 = MAX(MIN(it1,jptlucu2),jptlucu1)
        zqst1=tlucua(it1)/pp(jl)
        zqst1=MIN(0.5_dp,zqst1)
        zqst1=zqst1/(1._dp-vtmpc1*zqst1)
        zdqsdt=(zqst1-zqsat)*1000._dp
        zlcdqsdt=MERGE(zdqsdt*tlucuc(it),zqsat*zcor*tlucub(it),LO)
        zcond1=(pq(jl,kk)-zqsat)/(1._dp+zlcdqsdt)
        pt(jl,kk)=pt(jl,kk)+tlucuc(it)*zcond1
        pq(jl,kk)=pq(jl,kk)-zcond1
520  END DO

     IF (lookupoverflow) CALL lookuperror ('cuadjtq (8) ')

  END IF
!
#ifdef _PROFILE
  CALL trace_stop ('cuadjtq', 30)
#endif

  END SUBROUTINE cuadjtq

END MODULE mo_cuadjtq
