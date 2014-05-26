!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
MODULE mo_cudescent

USE mo_kind,              ONLY : wp
USE mo_physical_constants,ONLY : grav, rd, vtmpc1
USE mo_echam_conv_constants,ONLY : lmfdudv, lmfdd, cmfdeps, cmfcmin, entrdd
USE mo_cuadjust,          ONLY : cuadjtq
!
IMPLICIT NONE

PRIVATE

PUBLIC :: cudlfs, cuddraf

CONTAINS

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
           kdtop,    lddraf                                           )
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
!
INTEGER, INTENT (IN) :: kbdim, klev, ktrac, kproma, klevp1
REAL(wp) :: ptenh(kbdim,klev),       pqenh(kbdim,klev),                &
            puen(kbdim,klev),        pven(kbdim,klev),                 &
            pgeoh(kbdim,klev),       paphp1(kbdim,klevp1),             &
            ptu(kbdim,klev),         pqu(kbdim,klev),                  &
            puu(kbdim,klev),         pvu(kbdim,klev),                  &
            pmfub(kbdim),            prfl(kbdim)
!
REAL(wp) :: ptd(kbdim,klev),         pqd(kbdim,klev),                  &
            pud(kbdim,klev),         pvd(kbdim,klev),                  &
            pmfd(kbdim,klev),        pmfds(kbdim,klev),                &
            pmfdq(kbdim,klev),       pdmfdp(kbdim,klev)
REAL(wp) :: pcpcu(kbdim,klev)
INTEGER  :: kcbot(kbdim),            kctop(kbdim),                     &
            kdtop(kbdim)
LOGICAL  :: ldcum(kbdim),            lddraf(kbdim)
!
REAL(wp) :: ztenwb(kbdim,klev),      zqenwb(kbdim,klev),               &
            zcond(kbdim)
REAL(wp) :: zph(kbdim)
INTEGER  :: loidx(kbdim)
LOGICAL  :: llo2(kbdim)
REAL(wp) :: pxtenh(kbdim,klev,ktrac),pxtu(kbdim,klev,ktrac),           &
            pxtd(kbdim,klev,ktrac),  pmfdxt(kbdim,klev,ktrac)
LOGICAL  :: llo3(kbdim)
INTEGER  :: jl, jk, ke, is, ik, icall, jt
REAL(wp) :: zttest, zqtest, zbuo, zmftop
!
!----------------------------------------------------------------------
!
!     1.           SET DEFAULT VALUES FOR DOWNDRAFTS
!                  ---------------------------------
!
  DO jl = 1, kproma
     lddraf(jl) = .FALSE.
     kdtop(jl) = klevp1
 END DO
!
  IF (.NOT.lmfdd) RETURN
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
  ke=klev-3
  level: DO jk=3,ke
!
!
!     2.1          CALCULATE WET-BULB TEMPERATURE AND MOISTURE
!                  FOR ENVIRONMENTAL AIR IN *CUADJTQ*
!                  -------------------------------------------
!
    is = 0
    DO jl = 1, kproma
      llo2(jl) = .FALSE.
      IF (ldcum(jl) .AND. prfl(jl) > 0.0_wp .AND. .NOT. lddraf(jl) .AND. (jk < kcbot(jl) .AND. jk > kctop(jl))) THEN
        is = is+1
        loidx(is) = jl
        llo2(jl) = .TRUE.
      ENDIF
    ENDDO
    IF (is == 0) CYCLE level
     DO 212 jl=1,kproma
        ztenwb(jl,jk)=ptenh(jl,jk)
        zqenwb(jl,jk)=pqenh(jl,jk)
        zph(jl)=paphp1(jl,jk)
212  END DO
!
     ik=jk
     icall=2
     CALL cuadjtq( kproma, kbdim, klev, ik, zph, ztenwb, zqenwb, loidx, is, icall)
!
!
!     2.2          DO MIXING OF CUMULUS AND ENVIRONMENTAL AIR
!                  AND CHECK FOR NEGATIVE BUOYANCY.
!                  THEN SET VALUES FOR DOWNDRAFT AT LFS.
!                  ----------------------------------------
!
!DIR$ IVDEP
!OCL NOVREC
     DO 222 jl=1,kproma
        llo3(jl)=.FALSE.
        IF(llo2(jl)) THEN
           zttest=0.5_wp*(ptu(jl,jk)+ztenwb(jl,jk))
           zqtest=0.5_wp*(pqu(jl,jk)+zqenwb(jl,jk))
           zbuo=zttest*(1._wp+vtmpc1*zqtest)-                          &
                         ptenh(jl,jk)*(1._wp+vtmpc1*pqenh(jl,jk))
           zcond(jl)=pqenh(jl,jk)-zqenwb(jl,jk)
           zmftop=-cmfdeps*pmfub(jl)
           IF(zbuo.LT.0._wp.AND.prfl(jl).GT.10._wp*zmftop*zcond(jl))   &
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
              pdmfdp(jl,jk-1)=-0.5_wp*pmfd(jl,jk)*zcond(jl)
              prfl(jl)=prfl(jl)+pdmfdp(jl,jk-1)
           END IF
        END IF
222  END DO
!
     DO 2224 jt=1,ktrac
        DO 2222 jl=1,kproma
           IF(llo3(jl)) THEN
              pxtd(jl,jk,jt)=0.5_wp*(pxtu(jl,jk,jt)+pxtenh(jl,jk,jt))
              pmfdxt(jl,jk,jt)=pmfd(jl,jk)*pxtd(jl,jk,jt)
           ENDIF
2222    END DO
2224 END DO
!
     IF(lmfdudv) THEN
        DO 224 jl=1,kproma
           IF(pmfd(jl,jk).LT.0._wp) THEN
              pud(jl,jk)=0.5_wp*(puu(jl,jk)+puen(jl,jk-1))
              pvd(jl,jk)=0.5_wp*(pvu(jl,jk)+pven(jl,jk-1))
           END IF
224     END DO
     END IF
!
   END DO level
   !
END SUBROUTINE cudlfs

SUBROUTINE cuddraf(  kproma, kbdim, klev, klevp1,                      &
           ptenh,    pqenh,    puen,     pven,                         &
           ktrac,                                                      &
           pxtenh,   pxtd,     pmfdxt,                                 &
           pgeoh,    paphp1,   prfl,                                   &
           ptd,      pqd,      pud,      pvd,                          &
           pmfd,     pmfds,    pmfdq,    pdmfdp,                       &
           pcpcu,                                                      &
           lddraf                                                     )
!
!          THIS ROUTINE CALCULATES CUMULUS DOWNDRAFT DESCENT
!
!          M.TIEDTKE         E.C.M.W.F.    12/86 MODIF. 12/89
!
!          PURPOSE.
!          --------
!          TO PRODUCE THE VERTICAL PROFILES FOR CUMULUS DOWNDRAFTS
!          (I.E. T,Q,U AND V AND FLUXES)
!
!          INTERFACE
!          ---------
!
!          THIS ROUTINE IS CALLED FROM *CUMASTR*.
!          INPUT IS T,Q,P,PHI,U,V AT HALF LEVELS.
!          IT RETURNS FLUXES OF S,Q AND EVAPORATION RATE
!          AND U,V AT LEVELS WHERE DOWNDRAFT OCCURS
!
!          METHOD.
!          --------
!          CALCULATE MOIST DESCENT FOR ENTRAINING/DETRAINING PLUME BY
!          A) MOVING AIR DRY-ADIABATICALLY TO NEXT LEVEL BELOW AND
!          B) CORRECTING FOR EVAPORATION TO OBTAIN SATURATED STATE.
!
!          EXTERNALS
!          ---------
!          *CUADJTQ* FOR ADJUSTING T AND Q DUE TO EVAPORATION IN
!          SATURATED DESCENT
!
!          REFERENCE
!          ---------
!          (TIEDTKE,1989)
!
!
INTEGER, INTENT (IN) :: kbdim, klev, ktrac, kproma, klevp1
!
REAL(wp) :: ptenh(kbdim,klev),       pqenh(kbdim,klev),                &
            puen(kbdim,klev),        pven(kbdim,klev),                 &
            pgeoh(kbdim,klev),       paphp1(kbdim,klevp1)
!
REAL(wp) :: ptd(kbdim,klev),         pqd(kbdim,klev),                  &
            pud(kbdim,klev),         pvd(kbdim,klev),                  &
            pmfd(kbdim,klev),        pmfds(kbdim,klev),                &
            pmfdq(kbdim,klev),       pdmfdp(kbdim,klev),               &
            prfl(kbdim)
REAL(wp) :: pcpcu(kbdim,klev)
LOGICAL  :: lddraf(kbdim)
!
REAL(wp) :: zdmfen(kbdim),           zdmfde(kbdim),                    &
            zcond(kbdim)
REAL(wp) :: zph(kbdim)
LOGICAL  :: llo2(kbdim)
INTEGER  :: loidx(kbdim)
REAL(wp) :: pxtenh(kbdim,klev,ktrac),pxtd(kbdim,klev,ktrac),           &
            pmfdxt(kbdim,klev,ktrac)
LOGICAL  :: llo1
INTEGER  :: jk, is, jl, itopde, jt, ik, icall
REAL(wp) :: zentr, zseen, zqeen, zsdde, zqdde, zmfdsk, zmfdqk, zxteen  &
          , zxtdde, zmfdxtk, zbuo, zdmfdp, zmfduk, zmfdvk
!
!----------------------------------------------------------------------
!
!     1.           CALCULATE MOIST DESCENT FOR CUMULUS DOWNDRAFT BY
!                     (A) CALCULATING ENTRAINMENT RATES, ASSUMING
!                         LINEAR DECREASE OF MASSFLUX IN PBL
!                     (B) DOING MOIST DESCENT - EVAPORATIVE COOLING
!                         AND MOISTENING IS CALCULATED IN *CUADJTQ*
!                     (C) CHECKING FOR NEGATIVE BUOYANCY AND
!                         SPECIFYING FINAL T,Q,U,V AND DOWNWARD FLUXES
!                    -------------------------------------------------
!
level: DO jk = 3, klev
  is = 0
  DO jl = 1, kproma     
    llo2(jl) = .FALSE.
    IF (lddraf(jl) .AND. pmfd(jl,jk-1) < 0.0_wp) THEN
      is = is+1
      loidx(is) = jl
      llo2(jl) = .TRUE.
    ENDIF
  ENDDO

  DO jl = 1, kproma
    zph(jl) = paphp1(jl,jk)
  END DO

     DO 122 jl=1,kproma
        IF(llo2(jl)) THEN
           zentr=entrdd*pmfd(jl,jk-1)*rd*ptenh(jl,jk-1)/               &
                          (grav*paphp1(jl,jk-1))*                         &
                          (paphp1(jl,jk)-paphp1(jl,jk-1))
           zdmfen(jl)=zentr
           zdmfde(jl)=zentr
        END IF
122  END DO
     itopde=klev-2
     IF(jk.GT.itopde) THEN
        DO 124 jl=1,kproma
           IF(llo2(jl)) THEN
              zdmfen(jl)=0._wp
              zdmfde(jl)=pmfd(jl,itopde)*                              &
                             (paphp1(jl,jk)-paphp1(jl,jk-1))/          &
                             (paphp1(jl,klevp1)-paphp1(jl,itopde))
           END IF
124     END DO
     END IF
!
     DO 126 jl=1,kproma
        IF(llo2(jl)) THEN
           pmfd(jl,jk)=pmfd(jl,jk-1)+zdmfen(jl)-zdmfde(jl)
           zseen=(pcpcu(jl,jk-1)*ptenh(jl,jk-1)+pgeoh(jl,jk-1))        &
                                                      *zdmfen(jl)
           zqeen=pqenh(jl,jk-1)*zdmfen(jl)
           zsdde=(pcpcu(jl,jk-1)*ptd(jl,jk-1)+pgeoh(jl,jk-1))*zdmfde(jl)
           zqdde=pqd(jl,jk-1)*zdmfde(jl)
           zmfdsk=pmfds(jl,jk-1)+zseen-zsdde
           zmfdqk=pmfdq(jl,jk-1)+zqeen-zqdde
           pqd(jl,jk)=zmfdqk*(1._wp/MIN(-cmfcmin,pmfd(jl,jk)))
           ptd(jl,jk)=(zmfdsk*(1._wp/MIN(-cmfcmin,pmfd(jl,jk)))-       &
                                   pgeoh(jl,jk))/pcpcu(jl,jk)
           ptd(jl,jk)=MIN(400._wp,ptd(jl,jk))
           ptd(jl,jk)=MAX(100._wp,ptd(jl,jk))
           zcond(jl)=pqd(jl,jk)
        END IF
126  END DO
!
     DO 1264 jt=1,ktrac
        DO 1262 jl=1,kproma
           IF(llo2(jl)) THEN
              zxteen=pxtenh(jl,jk-1,jt)*zdmfen(jl)
              zxtdde=pxtd(jl,jk-1,jt)*zdmfde(jl)
              zmfdxtk=pmfdxt(jl,jk-1,jt)+zxteen-zxtdde
              pxtd(jl,jk,jt)=zmfdxtk*(1._wp/MIN(-cmfcmin,pmfd(jl,jk)))
           ENDIF
1262    END DO
1264 END DO
!
!
     ik=jk
     icall=2
     CALL cuadjtq(kproma, kbdim, klev, ik, zph, ptd, pqd, loidx, is, icall) 
!
!
     DO 150 jl=1,kproma
        IF(llo2(jl)) THEN
           zcond(jl)=zcond(jl)-pqd(jl,jk)
           zbuo=ptd(jl,jk)*(1._wp+vtmpc1*pqd(jl,jk))-                  &
                       ptenh(jl,jk)*(1._wp+vtmpc1*pqenh(jl,jk))
           llo1=zbuo.LT.0._wp.AND.                                     &
                             (prfl(jl)-pmfd(jl,jk)*zcond(jl).GT.0._wp)
           pmfd(jl,jk)=MERGE(pmfd(jl,jk),0._wp,llo1)
           pmfds(jl,jk)=(pcpcu(jl,jk)*ptd(jl,jk)+pgeoh(jl,jk))         &
                                                    *pmfd(jl,jk)
           pmfdq(jl,jk)=pqd(jl,jk)*pmfd(jl,jk)
           zdmfdp=-pmfd(jl,jk)*zcond(jl)
           pdmfdp(jl,jk-1)=zdmfdp
           prfl(jl)=prfl(jl)+zdmfdp
        END IF
150  END DO
!
     DO 1504 jt=1,ktrac
        DO 1502 jl=1,kproma
           IF(llo2(jl)) THEN
              pmfdxt(jl,jk,jt)=pxtd(jl,jk,jt)*pmfd(jl,jk)
           ENDIF
1502    END DO
1504 END DO
!
     IF(lmfdudv) THEN
        DO 160 jl=1,kproma
           IF(llo2(jl).AND.pmfd(jl,jk).LT.0._wp) THEN
              zmfduk=pmfd(jl,jk-1)*pud(jl,jk-1)+                       &
                              zdmfen(jl)*puen(jl,jk-1)-                &
                              zdmfde(jl)*pud(jl,jk-1)
              zmfdvk=pmfd(jl,jk-1)*pvd(jl,jk-1)+                       &
                              zdmfen(jl)*pven(jl,jk-1)-                &
                              zdmfde(jl)*pvd(jl,jk-1)
              pud(jl,jk)=zmfduk*(1._wp/MIN(-cmfcmin,pmfd(jl,jk)))
              pvd(jl,jk)=zmfdvk*(1._wp/MIN(-cmfcmin,pmfd(jl,jk)))
           END IF
160     END DO
     END IF
!
   END DO level

END SUBROUTINE cuddraf

END MODULE mo_cudescent
