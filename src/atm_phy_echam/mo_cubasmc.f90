!>
!! Module contains subroutine cubasmc.
!!
!! @author M. Tiedtke, ECMWF
!!
!! @par Revision History
!! M. Tiedtke, ECMWF (1989)
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
MODULE mo_cubasmc

  USE mo_kind,               ONLY : dp
#ifdef __ICON__
  USE mo_physical_constants, ONLY : g=>grav
  USE mo_echam_conv_nml,     ONLY : entrmid, cmfcmin, cmfcmax
#else
  USE mo_constants,     ONLY : g
  USE mo_cumulus_flux,  ONLY : entrmid, cmfcmin, cmfcmax
#endif

  IMPLICIT NONE
  PRIVATE
  PUBLIC  :: cubasmc

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

CONTAINS
  !>
  !!
SUBROUTINE cubasmc(  lmfdudv, &
           kproma, kbdim, klev, kk, klab,                              &
           pten,     pqen,     pqsen,    puen,     pven,               &
           ktrac,                                                      &
           pxten,    pxtu,     pmfuxt,                                 &
           pverv,    pgeo,     pgeoh,    ldcum,    ktype,              &
           pmfu,     pmfub,    pentr,    kcbot,                        &
           ptu,      pqu,      plu,      puu,      pvu,                &
           pmfus,    pmfuq,    pmful,    pdmfup,   pmfuu,              &
           pcpen,                                                      &
           pmfuv,                                                      &
!--------------------------------added by Junhua Zhang---------------------------
           plul,     plui,     pmfull,   pmfuli)
!------------------------------------end-----------------------------------------
!
!          M.TIEDTKE         E.C.M.W.F.     12/89
!
!          PURPOSE.
!          --------
!          THIS ROUTINE CALCULATES CLOUD BASE VALUES
!          FOR MIDLEVEL CONVECTION
!
!          INTERFACE
!          ---------
!
!          THIS ROUTINE IS CALLED FROM *CUASC*.
!          INPUT ARE ENVIRONMENTAL VALUES T,Q ETC
!          IT RETURNS CLOUDBASE VALUES FOR MIDLEVEL CONVECTION
!
!          METHOD.
!          --------
!          S. TIEDTKE (1989)
!
!          EXTERNALS
!          ---------
!          NONE
!
!
LOGICAL, INTENT (IN) :: lmfdudv
INTEGER, INTENT (IN) :: kbdim, klev, ktrac, kproma, kk
REAL(dp) :: pten(kbdim,klev),        pqen(kbdim,klev),                 &
            puen(kbdim,klev),        pven(kbdim,klev),                 &
            pqsen(kbdim,klev),       pverv(kbdim,klev),                &
            pgeo(kbdim,klev),        pgeoh(kbdim,klev)
!
REAL(dp) :: ptu(kbdim,klev),         pqu(kbdim,klev),                  &
            puu(kbdim,klev),         pvu(kbdim,klev),                  &
            plu(kbdim,klev),         pmfu(kbdim,klev),                 &
            pmfub(kbdim),            pentr(kbdim),                     &
            pmfus(kbdim,klev),       pmfuq(kbdim,klev),                &
            pmful(kbdim,klev),       pdmfup(kbdim,klev),               &
            pmfuu(kbdim),            pmfuv(kbdim)
REAL(dp) :: pcpen(kbdim,klev)
INTEGER  :: ktype(kbdim),            kcbot(kbdim),                     &
            klab(kbdim,klev)
LOGICAL  :: ldcum(kbdim)
!
REAL(dp) :: pxten(kbdim,klev,ktrac), pxtu(kbdim,klev,ktrac),           &
            pmfuxt(kbdim,klev,ktrac)
LOGICAL  :: llo3(kbdim)
INTEGER  :: jl, jt
REAL(dp) :: zzzmb
!
!
!-----------------------added by Junhua Zhang for Micro---------------
REAL(dp) :: plul(kbdim,klev),           plui(kbdim,klev)
REAL(dp) :: pmfull(kbdim,klev),         pmfuli(kbdim,klev)
!-------------------------------------end-----------------------------
!
!----------------------------------------------------------------------
!
!*    1.           CALCULATE ENTRAINMENT AND DETRAINMENT RATES
!                  -------------------------------------------
!
!100 CONTINUE
!DIR$ IVDEP
!OCL NOVREC
  DO 150 jl=1,kproma
     llo3(jl)=.FALSE.
     IF(.NOT.ldcum(jl) .AND. klab(jl,kk+1) .EQ. 0                       &
              .AND. pqen(jl,kk)    .GT. 0.90_dp*pqsen(jl,kk)            &
              .AND. pverv(jl,kk)   .LT. 0.0_dp                          &
              .AND. pgeoh(jl,kk)/g .GT. 1500.0_dp) THEN
        llo3(jl)=.TRUE.
        ptu(jl,kk+1)=(pcpen(jl,kk)*pten(jl,kk)                         &
                        +pgeo(jl,kk)-pgeoh(jl,kk+1))/pcpen(jl,kk)
        pqu(jl,kk+1)=pqen(jl,kk)
        plu(jl,kk+1)=0._dp
!-----------------------added by Junhua Zhang for Micro---------------
        plul(jl,kk+1)=0._dp
        plui(jl,kk+1)=0._dp
!-------------------------------------end-----------------------------
        zzzmb=MAX(cmfcmin,-pverv(jl,kk)/g)
        zzzmb=MIN(zzzmb,cmfcmax)
        pmfub(jl)=zzzmb
        pmfu(jl,kk+1)=pmfub(jl)
        pmfus(jl,kk+1)=pmfub(jl)*(pcpen(jl,kk)*ptu(jl,kk+1)            &
                                        +pgeoh(jl,kk+1))
        pmfuq(jl,kk+1)=pmfub(jl)*pqu(jl,kk+1)
        pmful(jl,kk+1)=0._dp
!-----------------------added by Junhua Zhang for Micro---------------
        pmfull(jl,kk+1)=0._dp
        pmfuli(jl,kk+1)=0._dp
!-------------------------------------end-----------------------------
        pdmfup(jl,kk+1)=0._dp
        kcbot(jl)=kk
        klab(jl,kk+1)=1
        ktype(jl)=3
        pentr(jl)=entrmid
        IF(lmfdudv) THEN
           puu(jl,kk+1)=puen(jl,kk)
           pvu(jl,kk+1)=pven(jl,kk)
           pmfuu(jl)=pmfub(jl)*puu(jl,kk+1)
           pmfuv(jl)=pmfub(jl)*pvu(jl,kk+1)
        END IF
     END IF
150 END DO
!DIR$ IVDEP
!OCL NOVREC
  DO 1504 jt=1,ktrac
     DO 1502 jl=1,kproma
        IF (llo3(jl)) THEN
           pxtu(jl,kk+1,jt)=pxten(jl,kk,jt)
           pmfuxt(jl,kk+1,jt)=pmfub(jl)*pxtu(jl,kk+1,jt)
        ENDIF
1502 END DO
1504 END DO
!
!
    RETURN
  END SUBROUTINE cubasmc

END MODULE mo_cubasmc
