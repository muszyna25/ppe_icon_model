!>
!! Module contains subroutine cubase.
!!
!! @author M. Tiedtke, ECMWF
!!
!! @par Revision History
!! M. Tiedtke, ECMWF (1986,1989)
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
MODULE mo_cubase

  USE mo_kind,               ONLY: dp
  USE mo_cuadjtq,            ONLY: cuadjtq_idx
#ifdef __ICON__
  USE mo_physical_constants, ONLY: vtmpc1
  USE mo_echam_conv_nml,     ONLY: lmfdudv, cbfac, cminbuoy, cmaxbuoy
#else
  USE mo_constants,          ONLY: vtmpc1
  USE mo_cumulus_flux,       ONLY: lmfdudv, cbfac, cminbuoy, cmaxbuoy
#endif

  IMPLICIT NONE
  PRIVATE
  PUBLIC  :: cubase

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

CONTAINS
  !>
  !!
SUBROUTINE cubase(   kproma, kbdim, klev, klevp1, klevm1,              &
           ptenh,    pqenh,    pgeoh,    paph,   pthvsig,              &
           ptu,      pqu,      plu,                                    &
           puen,     pven,     puu,      pvu,                          &
           pcpcu,                                                      &
           ldcum,    kcbot,    klab)
!
!          THIS ROUTINE CALCULATES CLOUD BASE VALUES (T AND Q)
!          FOR CUMULUS PARAMETERIZATION
!
!          M.TIEDTKE         E.C.M.W.F.     7/86 MODIF. 12/89
!
!          PURPOSE.
!          --------
!          TO PRODUCE CLOUD BASE VALUES FOR CU-PARAMETRIZATION
!
!          INTERFACE
!          ---------
!          THIS ROUTINE IS CALLED FROM *CUMASTR*.
!          INPUT ARE ENVIRONM. VALUES OF T,Q,P,PHI AT HALF LEVELS.
!          IT RETURNS CLOUD BASE VALUES AND FLAGS AS FOLLOWS;
!                 KLAB=1 FOR SUBCLOUD LEVELS
!                 KLAB=2 FOR CONDENSATION LEVEL
!
!          METHOD.
!          --------
!          LIFT SURFACE AIR DRY-ADIABATICALLY TO CLOUD BASE
!          (NON ENTRAINING PLUME,I.E.CONSTANT MASSFLUX)
!
!          EXTERNALS
!          ---------
!          *CUADJTQ* FOR ADJUSTING T AND Q DUE TO CONDENSATION IN ASCENT
!
INTEGER, INTENT (IN) :: kproma, kbdim, klev, klevp1, klevm1

REAL(dp):: ptenh(kbdim,klev),       pqenh(kbdim,klev),                 &
           pgeoh(kbdim,klev),       paph(kbdim,klevp1),                &
           pthvsig(kbdim)
!
REAL(dp):: ptu(kbdim,klev),         pqu(kbdim,klev),                   &
           plu(kbdim,klev)
REAL(dp):: puen(kbdim,klev),        pven(kbdim,klev),                  &
           puu(kbdim,klev),         pvu(kbdim,klev)
REAL(dp):: pcpcu(kbdim,klev)
INTEGER :: klab(kbdim,klev),        kcbot(kbdim)
LOGICAL :: ldcum(kbdim)
!
REAL(dp):: zqold(kbdim)
REAL(dp):: zph(kbdim)
INTEGER :: loidx(kbdim)

INTEGER :: jl, jk, nl, is, ik, ikb, icall
REAL(dp):: zbuo, zz, zlift
!
!
!----------------------------------------------------------------------
!
!     1.           INITIALIZE VALUES AT LIFTING LEVEL
!                  ----------------------------------
!
!100 CONTINUE
  DO 110 jl=1,kproma
     klab(jl,klev)=1
     kcbot(jl)=klevm1
     ldcum(jl)=.FALSE.
     puu(jl,klev)=puen(jl,klev)*(paph(jl,klevp1)-paph(jl,klev))
     pvu(jl,klev)=pven(jl,klev)*(paph(jl,klevp1)-paph(jl,klev))
110 END DO
!
!
!----------------------------------------------------------------------
!
!     2.0          DO ASCENT IN SUBCLOUD LAYER,
!                  CHECK FOR EXISTENCE OF CONDENSATION LEVEL,
!                  ADJUST T,Q AND L ACCORDINGLY IN *CUADJTQ*,
!                  CHECK FOR BUOYANCY AND SET FLAGS
!                  -------------------------------------
!
!200 CONTINUE
  DO 290 jk=klevm1,2,-1
     is=0
     DO 210 jl=1,kproma
        IF (klab(jl,jk+1).EQ.1) THEN
           is = is + 1
           loidx(is) = jl
        END IF
        zph(jl)=paph(jl,jk)
210  END DO
     IF(is.EQ.0) CYCLE !GOTO 290
!CDIR NODEP
!IBM* ASSERT(NODEPS)
     DO 220 nl=1,is
        jl = loidx(nl)
        zlift=MAX(cminbuoy,MIN(cmaxbuoy,pthvsig(jl)*cbfac))
        zlift=MIN(zlift,1.0_dp)
        pqu(jl,jk)=pqu(jl,jk+1)
        ptu(jl,jk)=(pcpcu(jl,jk+1)*ptu(jl,jk+1)+pgeoh(jl,jk+1)      &
             -pgeoh(jl,jk))/pcpcu(jl,jk)
        zbuo=ptu(jl,jk)*(1._dp+vtmpc1*pqu(jl,jk))-ptenh(jl,jk)      &
             *(1._dp+vtmpc1*pqenh(jl,jk))+zlift
        IF(zbuo.GT.0._dp) klab(jl,jk)=1
        zqold(jl)=pqu(jl,jk)
220  END DO
!
     ik=jk
     icall=1
     CALL cuadjtq_idx(kproma, kbdim, klev, ik,                             &
                      zph,    ptu,   pqu,  loidx,   is,   icall)
!
!DIR$ IVDEP
!OCL NOVREC
!CDIR NODEP
     DO 240 nl=1,is
        jl = loidx(nl)
        IF(pqu(jl,jk).LT.zqold(jl)) THEN
           klab(jl,jk)=2
           zlift=MAX(cminbuoy,MIN(cmaxbuoy,pthvsig(jl)*cbfac))
           zlift=MIN(zlift,1.0_dp)
           plu(jl,jk)=plu(jl,jk)+zqold(jl)-pqu(jl,jk)
           zbuo=ptu(jl,jk)*(1._dp+vtmpc1*pqu(jl,jk)-plu(jl,jk))-       &
                        ptenh(jl,jk)*(1._dp+vtmpc1*pqenh(jl,jk))+zlift
           IF(zbuo.GT.0._dp) THEN
              kcbot(jl)=jk
              ldcum(jl)=.TRUE.
           END IF
        END IF
240  END DO
!
!             CALCULATE AVERAGES OF U AND V FOR SUBCLOUD ARA,.
!             THE VALUES WILL BE USED TO DEFINE CLOUD BASE VALUES.
!
     IF(lmfdudv) THEN
        DO 250 jl=1,kproma
           IF(jk.GE.kcbot(jl)) THEN
              puu(jl,klev)=puu(jl,klev)+                               &
                             puen(jl,jk)*(paph(jl,jk+1)-paph(jl,jk))
              pvu(jl,klev)=pvu(jl,klev)+                               &
                             pven(jl,jk)*(paph(jl,jk+1)-paph(jl,jk))
           END IF
250     END DO
     END IF
!
290 END DO
!
!
  IF(lmfdudv) THEN
     DO 310 jl=1,kproma
        IF(ldcum(jl)) THEN
           ikb=kcbot(jl)
           zz=1._dp/(paph(jl,klevp1)-paph(jl,ikb))
           puu(jl,klev)=puu(jl,klev)*zz
           pvu(jl,klev)=pvu(jl,klev)*zz
        ELSE
           puu(jl,klev)=puen(jl,klevm1)
           pvu(jl,klev)=pven(jl,klevm1)
        END IF
310  END DO
  END IF
!
    RETURN
  END SUBROUTINE cubase

END MODULE mo_cubase
