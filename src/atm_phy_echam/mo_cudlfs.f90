!>
!! Module contains subroutine cudlfs
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
MODULE mo_cudlfs

  USE mo_kind,               ONLY : dp
  USE mo_cuadjtq,            ONLY : cuadjtq
  USE mo_cuadjtqi,           ONLY : cuadjtqi
#ifdef __ICON__
  USE mo_physical_constants, ONLY : vtmpc1
  USE mo_echam_conv_params,  ONLY : lmfdudv, lmfdd, cmfdeps, ncvmicro
#else
  USE mo_constants,          ONLY : vtmpc1
  USE mo_cumulus_flux,       ONLY : lmfdudv, lmfdd, cmfdeps
  USE mo_param_switches,     ONLY : ncvmicro
#endif

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: cudlfs

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

CONTAINS
  !>
  !!
SUBROUTINE cudlfs(   kproma, kbdim, klev, klevp1,                      &
           ptenh,    pqenh,    puen,     pven,                         &
           ktrac,                                                      &
           pxtenh,   pxtu,     pxtd,     pmfdxt,                       &
           pgeoh,    paphp1,                                           &
           ptu,      pqu,      puu,      pvu,                          &
           ldcum,    kcbot,    kctop,    pmfub,    prfl,               &
           ptd,      pqd,      pud,      pvd,                          &
           pmfd,     pmfds,    pmfdq,    pdmfdp,                       &
           pcpcu,                                                      &
           kdtop,    lddraf,                                           &
!-----------------------added by Junhua Zhang for Micro---------------
           plui)
!-------------------------------------end-----------------------------
!
!          THIS ROUTINE CALCULATES LEVEL OF FREE SINKING FOR
!          CUMULUS DOWNDRAFTS AND SPECIFIES T,Q,U AND V VALUES
!
!          M.TIEDTKE         E.C.M.W.F.    12/86 MODIF. 12/89
!
!          PURPOSE.
!          --------
!          TO PRODUCE LFS-VALUES FOR CUMULUS DOWNDRAFTS
!          FOR MASSFLUX CUMULUS PARAMETERIZATION
!
!          INTERFACE
!          ---------
!          THIS ROUTINE IS CALLED FROM *CUMASTR*.
!          INPUT ARE ENVIRONMENTAL VALUES OF T,Q,U,V,P,PHI
!          AND UPDRAFT VALUES T,Q,U AND V AND ALSO
!          CLOUD BASE MASSFLUX AND CU-PRECIPITATION RATE.
!          IT RETURNS T,Q,U AND V VALUES AND MASSFLUX AT LFS.
!
!          METHOD.
!          --------
!          CHECK FOR NEGATIVE BUOYANCY OF AIR OF EQUAL PARTS OF
!          MOIST ENVIRONMENTAL AIR AND CLOUD AIR.
!
!          EXTERNALS
!          ---------
!          *CUADJTQ* FOR CALCULATING WET BULB T AND Q AT LFS
!
INTEGER, INTENT (IN) :: kbdim, klev, ktrac, kproma, klevp1
REAL(dp) :: ptenh(kbdim,klev),       pqenh(kbdim,klev),                &
            puen(kbdim,klev),        pven(kbdim,klev),                 &
            pgeoh(kbdim,klev),       paphp1(kbdim,klevp1),             &
            ptu(kbdim,klev),         pqu(kbdim,klev),                  &
            puu(kbdim,klev),         pvu(kbdim,klev),                  &
            pmfub(kbdim),            prfl(kbdim)
!
REAL(dp) :: ptd(kbdim,klev),         pqd(kbdim,klev),                  &
            pud(kbdim,klev),         pvd(kbdim,klev),                  &
            pmfd(kbdim,klev),        pmfds(kbdim,klev),                &
            pmfdq(kbdim,klev),       pdmfdp(kbdim,klev)
REAL(dp) :: pcpcu(kbdim,klev)
INTEGER  :: kcbot(kbdim),            kctop(kbdim),                     &
            kdtop(kbdim)
LOGICAL  :: ldcum(kbdim),            lddraf(kbdim)
!
REAL(dp) :: ztenwb(kbdim,klev),      zqenwb(kbdim,klev),               &
            zcond(kbdim)
REAL(dp) :: zph(kbdim)
LOGICAL  :: llo2(kbdim)
REAL(dp) :: pxtenh(kbdim,klev,ktrac),pxtu(kbdim,klev,ktrac),           &
            pxtd(kbdim,klev,ktrac),  pmfdxt(kbdim,klev,ktrac)
LOGICAL  :: llo3(kbdim)
INTEGER  :: jl, jk, ke, is, ik, icall, jt
REAL(DP) :: zttest, zqtest, zbuo, zmftop
!
!-----------------------added by Junhua Zhang for Micro---------------
REAL(dp) :: plui(kbdim,klev)
!-------------------------------------end-----------------------------
!
!----------------------------------------------------------------------
!
!     1.           SET DEFAULT VALUES FOR DOWNDRAFTS
!                  ---------------------------------
!
!100 CONTINUE
  DO 110 jl=1,kproma
     lddraf(jl)=.FALSE.
     kdtop(jl)=klevp1
110 END DO
!
  IF(.NOT.lmfdd) RETURN !go to 300
!
!
!----------------------------------------------------------------------
!
!     2.           DETERMINE LEVEL OF FREE SINKING BY
!                  DOING A SCAN FROM TOP TO BASE OF CUMULUS CLOUDS
!
!                  FOR EVERY POINT AND PROCEED AS FOLLOWS:
!
!                    (1) DETEMINE WET BULB ENVIRONMENTAL T AND Q
!                    (2) DO MIXING WITH CUMULUS CLOUD AIR
!                    (3) CHECK FOR NEGATIVE BUOYANCY
!
!                  THE ASSUMPTION IS THAT AIR OF DOWNDRAFTS IS MIXTURE
!                  OF 50% CLOUD AIR + 50% ENVIRONMENTAL AIR AT WET BULB
!                  TEMPERATURE (I.E. WHICH BECAME SATURATED DUE TO
!                  EVAPORATION OF RAIN AND CLOUD WATER)
!                  ----------------------------------------------------
!
!200 CONTINUE
!
  ke=klev-3
  DO 290 jk=3,ke
!
!
!     2.1          CALCULATE WET-BULB TEMPERATURE AND MOISTURE
!                  FOR ENVIRONMENTAL AIR IN *CUADJTQ*
!                  -------------------------------------------
!
!210  CONTINUE
     is=0
     DO 212 jl=1,kproma
        ztenwb(jl,jk)=ptenh(jl,jk)
        zqenwb(jl,jk)=pqenh(jl,jk)
        zph(jl)=paphp1(jl,jk)
        llo2(jl)=ldcum(jl).AND.prfl(jl).GT.0._dp.AND..NOT.lddraf(jl)   &
                          .AND.(jk.LT.kcbot(jl).AND.jk.GT.kctop(jl))
        is=is+MERGE(1,0,llo2(jl))
212  END DO
     IF(is.EQ.0) CYCLE !go to 290
!
     ik=jk
     icall=2
     IF (ncvmicro .GT. 0) THEN
        CALL cuadjtqi(kproma, kbdim, klev, ik,                         &
                  zph,      ztenwb,   zqenwb,   llo2,     icall,       &
!-----------------------added by Junhua Zhang for Micro---------------
                  plui)
!-------------------------------------end-----------------------------
     ELSE
        CALL cuadjtq(kproma, kbdim, klev, ik,                          &
                  zph,      ztenwb,   zqenwb,   llo2,     icall)
     ENDIF
!
!
!     2.2          DO MIXING OF CUMULUS AND ENVIRONMENTAL AIR
!                  AND CHECK FOR NEGATIVE BUOYANCY.
!                  THEN SET VALUES FOR DOWNDRAFT AT LFS.
!                  ----------------------------------------
!
!220  CONTINUE
!DIR$ IVDEP
!OCL NOVREC
     DO 222 jl=1,kproma
        llo3(jl)=.FALSE.
        IF(llo2(jl)) THEN
           zttest=0.5_dp*(ptu(jl,jk)+ztenwb(jl,jk))
           zqtest=0.5_dp*(pqu(jl,jk)+zqenwb(jl,jk))
           zbuo=zttest*(1._dp+vtmpc1*zqtest)-                          &
                         ptenh(jl,jk)*(1._dp+vtmpc1*pqenh(jl,jk))
           zcond(jl)=pqenh(jl,jk)-zqenwb(jl,jk)
           zmftop=-cmfdeps*pmfub(jl)
           IF(zbuo.LT.0._dp.AND.prfl(jl).GT.10._dp*zmftop*zcond(jl))   &
                                                                  THEN
              llo3(jl)=.TRUE.
              kdtop(jl)=jk
              lddraf(jl)=.TRUE.
              ptd(jl,jk)=zttest
              pqd(jl,jk)=zqtest
              pmfd(jl,jk)=zmftop
              pmfds(jl,jk)=pmfd(jl,jk)*(pcpcu(jl,jk)*ptd(jl,jk)        &
                                               +pgeoh(jl,jk))
              pmfdq(jl,jk)=pmfd(jl,jk)*pqd(jl,jk)
              pdmfdp(jl,jk-1)=-0.5_dp*pmfd(jl,jk)*zcond(jl)
              prfl(jl)=prfl(jl)+pdmfdp(jl,jk-1)
           END IF
        END IF
222  END DO
!
     DO 2224 jt=1,ktrac
        DO 2222 jl=1,kproma
           IF(llo3(jl)) THEN
              pxtd(jl,jk,jt)=0.5_dp*(pxtu(jl,jk,jt)+pxtenh(jl,jk,jt))
              pmfdxt(jl,jk,jt)=pmfd(jl,jk)*pxtd(jl,jk,jt)
           ENDIF
2222    END DO
2224 END DO
!
     IF(lmfdudv) THEN
        DO 224 jl=1,kproma
           IF(pmfd(jl,jk).LT.0._dp) THEN
              pud(jl,jk)=0.5_dp*(puu(jl,jk)+puen(jl,jk-1))
              pvd(jl,jk)=0.5_dp*(pvu(jl,jk)+pven(jl,jk-1))
           END IF
224     END DO
     END IF
!
290 END DO
!
!300 CONTINUE
!  RETURN
  END SUBROUTINE cudlfs

END MODULE mo_cudlfs
