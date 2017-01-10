!>
!! @brief
!!       *cuflx*
!!        This routine does the final calculation of convective fluxes in the
!!        cloud layer and in the subcloud layer
!!        *cudtdq*
!!        This routine updates t and q tendencies, precipitation rates, does
!!        global diagnostics
!!        *cududv*
!!        This routine updates u and v tendencies, does global diagnostic of
!!        dissipation
!!
!! @author M. Tiedtke, ECMWF,    Dec 1989
!!
!! @par Revision History
!! - Taken from ECHAM6.3, unified for ICON/ECHAM by Monika Esch, MPI-M (2015-06)
!!
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_cufluxdts

  USE mo_kind,               ONLY: wp
  USE mo_physical_constants, ONLY: grav, alv, als, alf, tmelt, cpd, vtmpc2
  !
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: cuflx, cudtdq, cududv

CONTAINS
  !>
  !!
  !!++mgs: zcucov and zdpevap now 1d vectors
  !!
  SUBROUTINE cuflx(    kproma, kbdim, klev, klevp1,                                  &
    &        pqen,     pqsen,    ptenh,    pqenh,                                    &
    &        ktrac,                                                                  &
    &        ptime_step_len,                                                         &
    &        cevapcu,                                                                &
    &        pxtenh,   pmfuxt,   pmfdxt,                                             &
    &        paphp1,   pgeoh,                                                        &
    &        kcbot,    kctop,    kdtop,                                              &
    &        ktype,    lddraf,   ldcum,                                              &
    &        pmfu,     pmfd,     pmfus,    pmfds,                                    &
    &        pmfuq,    pmfdq,    pmful,                                              &
    &        pdmfup,   pdmfdp,   prfl,                                               &
    &        pcpcu,                                                                  &
    &        pten,     psfl,     pdpmel,   ktopm2                                   )
    !
    INTEGER, INTENT (IN) :: kproma, kbdim, klev, klevp1, ktrac
    INTEGER, INTENT (OUT):: ktopm2
    REAL(wp),INTENT (IN) :: ptime_step_len
    REAL(wp),INTENT (IN) :: cevapcu(klev)
    REAL(wp):: pqen(kbdim,klev),        pqsen(kbdim,klev),                           &
      &        ptenh(kbdim,klev),       pqenh(kbdim,klev),                           &
      &        paphp1(kbdim,klevp1),    pgeoh(kbdim,klev)
    REAL(wp):: pmfu(kbdim,klev),        pmfd(kbdim,klev),                            &
      &        pmfus(kbdim,klev),       pmfds(kbdim,klev),                           &
      &        pmfuq(kbdim,klev),       pmfdq(kbdim,klev),                           &
      &        pdmfup(kbdim,klev),      pdmfdp(kbdim,klev),                          &
      &        pmful(kbdim,klev),                                                    &
      &        prfl(kbdim)
    REAL(wp):: pten(kbdim,klev),        pdpmel(kbdim,klev),                          &
      &        psfl(kbdim)
    REAL(wp):: pcpcu(kbdim,klev)
    INTEGER :: kcbot(kbdim),            kctop(kbdim),                                &
      &        kdtop(kbdim),            ktype(kbdim)
    LOGICAL :: lddraf(kbdim),           ldcum(kbdim)
    REAL(wp):: pxtenh(kbdim,klev,ktrac),                                             &
      &        pmfuxt(kbdim,klev,ktrac),pmfdxt(kbdim,klev,ktrac)
    !
    INTEGER :: jl, jk, jt, ikb
    !++mgs
    REAL(wp):: zcons1, zcons2,         ztmelp2, zzp, zfac, zsnmlt                    &
      &      , zrfl, zrnew, zrmin, zrfln, zdrfl, zrsum
    REAL(wp):: zpsubcl(kbdim), zcucov(kbdim), zdpevap(kbdim)
    !--mgs
    !
    !*             Specify constants
    !
    zcons1=cpd/(alf*grav*ptime_step_len)
    zcons2=1._wp/(grav*ptime_step_len)
    ztmelp2=tmelt+2._wp
    !
    !*    1.0          Determine final convection fluxes
    !                  ---------------------------------
    !
    !  itop=klev
    DO jl=1,kproma
      ! itop=MIN(itop,kctop(jl))
      IF(.NOT.ldcum(jl).OR.kdtop(jl).LT.kctop(jl)) lddraf(jl)=.FALSE.
      IF(.NOT.ldcum(jl)) ktype(jl)=0
    END DO
    ktopm2=1
    DO jk=ktopm2,klev
!DIR$ IVDEP
!OCL NOVREC
      DO jl=1,kproma
        IF(ldcum(jl).AND.jk.GE.kctop(jl)-1) THEN
          pmfus(jl,jk)=pmfus(jl,jk)-pmfu(jl,jk)*                                     &
            &                  (pcpcu(jl,jk)*ptenh(jl,jk)+pgeoh(jl,jk))
          pmfuq(jl,jk)=pmfuq(jl,jk)-pmfu(jl,jk)*pqenh(jl,jk)
          IF(lddraf(jl).AND.jk.GE.kdtop(jl)) THEN
            pmfds(jl,jk)=pmfds(jl,jk)-pmfd(jl,jk)*                                   &
              &                (pcpcu(jl,jk)*ptenh(jl,jk)+pgeoh(jl,jk))
            pmfdq(jl,jk)=pmfdq(jl,jk)-pmfd(jl,jk)*pqenh(jl,jk)
          ELSE
            pmfd(jl,jk)=0._wp
            pmfds(jl,jk)=0._wp
            pmfdq(jl,jk)=0._wp
            pdmfdp(jl,jk-1)=0._wp
          END IF
        END IF
      END DO
      !
      DO jt=1,ktrac
        DO jl=1,kproma
          IF(ldcum(jl).AND.jk.GE.kctop(jl)-1) THEN
            pmfuxt(jl,jk,jt)=pmfuxt(jl,jk,jt)-pmfu(jl,jk)*pxtenh(jl,jk,jt)
            IF(lddraf(jl).AND.jk.GE.kdtop(jl)) THEN
              pmfdxt(jl,jk,jt)=pmfdxt(jl,jk,jt)-pmfd(jl,jk)*pxtenh(jl,jk,jt)
            ELSE
              pmfdxt(jl,jk,jt)=0._wp
            ENDIF
          ELSE
            pmfuxt(jl,jk,jt)=0._wp
            pmfdxt(jl,jk,jt)=0._wp
          ENDIF
        END DO
      END DO
      !
    END DO
    DO jk=ktopm2,klev
!DIR$ IVDEP
!OCL NOVREC
      DO jl=1,kproma
        IF(ldcum(jl).AND.jk.GT.kcbot(jl)) THEN
          ikb=kcbot(jl)
          zzp=((paphp1(jl,klevp1)-paphp1(jl,jk))/(paphp1(jl,klevp1)-paphp1(jl,ikb)))
          zzp=MERGE(zzp**2,zzp,ktype(jl).EQ.3)
          pmfu(jl,jk)=pmfu(jl,ikb)*zzp
          pmfus(jl,jk)=pmfus(jl,ikb)*zzp
          pmfuq(jl,jk)=pmfuq(jl,ikb)*zzp
          pmful(jl,jk)=pmful(jl,ikb)*zzp
        END IF
      END DO
      !
      DO jt=1,ktrac
!DIR$ IVDEP
!OCL NOVREC
        DO jl=1,kproma
          IF(ldcum(jl).AND.jk.GT.kcbot(jl)) THEN
            ikb=kcbot(jl)
            zzp=(paphp1(jl,klevp1)-paphp1(jl,jk))/(paphp1(jl,klevp1)-paphp1(jl,ikb))
            zzp=MERGE(zzp**2,zzp,ktype(jl).EQ.3)
            pmfuxt(jl,jk,jt)=pmfuxt(jl,ikb,jt)*zzp
          ENDIF
        END DO
      END DO
    END DO
    !
    !*    2.  Calculate rain/snow fall rates
    !*        Calculate melting of snow
    !*        Calculate evaporation of precipitation
    !         --------------------------------------
    !
    prfl(:)=0._wp
    psfl(:)=0._wp
      
    DO jk=ktopm2,klev
      DO jl=1,kproma
        IF(ldcum(jl)) THEN
          IF(pten(jl,jk).GT.tmelt) THEN
            prfl(jl)=prfl(jl)+pdmfup(jl,jk)+pdmfdp(jl,jk)
            IF(psfl(jl).GT.0._wp.AND.pten(jl,jk).GT.ztmelp2) THEN
              zfac=zcons1*(1._wp+vtmpc2*pqen(jl,jk))*(paphp1(jl,jk+1)-paphp1(jl,jk))
              zsnmlt=MIN(psfl(jl),zfac*(pten(jl,jk)-ztmelp2))
              pdpmel(jl,jk)=zsnmlt
              psfl(jl)=psfl(jl)-zsnmlt
              prfl(jl)=prfl(jl)+zsnmlt
            END IF
          ELSE
            psfl(jl)=psfl(jl)+pdmfup(jl,jk)+pdmfdp(jl,jk)
          END IF
        END IF
      END DO
    END DO
    DO jl=1,kproma
      prfl(jl)=MAX(prfl(jl),0._wp)
      psfl(jl)=MAX(psfl(jl),0._wp)
      zpsubcl(jl)=prfl(jl)+psfl(jl)
    END DO
    DO jk=ktopm2,klev
      zdpevap(1:kproma) = 0._wp
      zcucov(1:kproma) = 0.05_wp

      DO jl=1,kproma
        IF(ldcum(jl).AND.jk.GE.kcbot(jl).AND.zpsubcl(jl).GT.1.e-20_wp) THEN
          zrfl=zpsubcl(jl)
          zrnew=(MAX(0._wp,SQRT(zrfl/zcucov(jl))-                                    &
            &    cevapcu(jk)*(paphp1(jl,jk+1)-paphp1(jl,jk))*                        &
            &    MAX(0._wp,pqsen(jl,jk)-pqen(jl,jk))))**2*zcucov(jl)
          zrmin=zrfl-zcucov(jl)*MAX(0._wp,0.8_wp*pqsen(jl,jk)-pqen(jl,jk))           &
            &   *zcons2*(paphp1(jl,jk+1)-paphp1(jl,jk))
          zrnew=MAX(zrnew,zrmin)
          zrfln=MAX(zrnew,0._wp)
          zdrfl=MIN(0._wp,zrfln-zrfl)
          zdpevap(jl) = -zdrfl
          pdmfup(jl,jk)=pdmfup(jl,jk)+zdrfl
          zpsubcl(jl)=zrfln
        END IF
      END DO

    END DO

    !!baustelle!! (?)
    DO jl=1,kproma
      zrsum=prfl(jl)+psfl(jl)
      zdpevap(jl)=zpsubcl(jl)-zrsum
      prfl(jl)=prfl(jl)+zdpevap(jl)*prfl(jl)*(1._wp/MAX(1.e-20_wp,zrsum))
      psfl(jl)=psfl(jl)+zdpevap(jl)*psfl(jl)*(1._wp/MAX(1.e-20_wp,zrsum))
    END DO
    !
  END SUBROUTINE cuflx
  !>
  !!
  SUBROUTINE cudtdq(kproma, kbdim, klev, klevp1, ktopm2, ldcum, ktrac,               &
    &               paphp1,   pten,                                                  &
    &               pmfuxt,   pmfdxt,                                                &
    &               pmfus,    pmfds,    pmfuq,    pmfdq,                             &
    &               pmful,    pdmfup,   pdmfdp,   plude,                             &
    &               pdpmel,   prfl,     psfl,                                        &
    &               palvsh,                                                          &
    &               pcon_dtrl,pcon_dtri,pcon_iqte,                                   &
    &               pq_cnv,   pqte_cnv, pxtte_cnv,                                   &
    &               pxtecl,   pxteci                                                )
    !
    INTEGER, INTENT (IN) :: kproma, kbdim, klev, klevp1, ktopm2, ktrac
    REAL(wp),INTENT(IN)  :: prfl (kbdim), psfl(kbdim)
    REAL(wp),INTENT(OUT) :: pcon_dtrl(kbdim), pcon_dtri(kbdim)
    REAL(wp),INTENT(OUT) :: pcon_iqte(kbdim) ! integrated qv tendency
    REAL(wp),INTENT(OUT) :: pq_cnv(kbdim,klev)
    REAL(wp),INTENT(OUT) :: pqte_cnv(kbdim,klev), pxtte_cnv(kbdim,klev,ktrac)
    REAL(wp),INTENT(OUT) :: pxtecl(kbdim,klev), pxteci(kbdim,klev)
    
    LOGICAL  llo1
    !
    REAL(wp) :: pten(kbdim,klev),        paphp1(kbdim,klevp1)
    REAL(wp) :: pmfus(kbdim,klev),       pmfds(kbdim,klev),                          &
      &         pmfuq(kbdim,klev),       pmfdq(kbdim,klev),                          &
      &         pmful(kbdim,klev),       plude(kbdim,klev),                          &
      &         pdmfup(kbdim,klev),      pdmfdp(kbdim,klev)
    REAL(wp) :: pdpmel(kbdim,klev)
    REAL(wp) :: palvsh(kbdim,klev)
    LOGICAL  :: ldcum(kbdim)
    REAL(wp) :: zqte_cnv(kbdim,klev)
    REAL(wp) :: pmfuxt(kbdim,klev,ktrac),pmfdxt(kbdim,klev,ktrac)
    REAL(wp) :: zxtecl(kbdim,klev),zxteci(kbdim,klev)
    INTEGER  :: jl, jk, jt
    REAL(wp) :: zalv
    pq_cnv   (:,:)    = 0.0_wp
    zqte_cnv (:,:)    = 0.0_wp
    pqte_cnv (:,:)    = 0.0_wp
    pxtte_cnv(:,:,:)  = 0.0_wp
    zxtecl   (:,:)    = 0.0_wp
    zxteci   (:,:)    = 0.0_wp
    pxtecl   (:,:)    = 0.0_wp
    pxteci   (:,:)    = 0.0_wp
    pcon_dtrl(:)      = 0.0_wp
    pcon_dtri(:)      = 0.0_wp
    pcon_iqte(:)      = 0.0_wp
    !----------------------------------------------------------------------
    !
    !*    2.0          Incrementation of t and q tendencies
    !                  ------------------------------------
    !
    DO jk=ktopm2,klev
      !
      IF(jk.LT.klev) THEN
        DO jl=1,kproma
          IF(ldcum(jl)) THEN
            llo1=(pten(jl,jk)-tmelt).GT.0._wp
            zalv=MERGE(alv,als,llo1)
            
            pq_cnv(jl,jk)   =   pmfus(jl,jk+1)-pmfus(jl,jk)+                        &
              &                 pmfds(jl,jk+1)-pmfds(jl,jk)-                        &
              &                 alf*pdpmel(jl,jk)-                                  &
              &                 palvsh(jl,jk+1)*pmful(jl,jk+1)+                     &
              &                 palvsh(jl,jk)*pmful(jl,jk)+                         &
              &                 zalv*(plude(jl,jk)+pdmfup(jl,jk)+pdmfdp(jl,jk))
            
            zqte_cnv(jl,jk) =   pmfuq(jl,jk+1)-pmfuq(jl,jk)+                        &
              &                 pmfdq(jl,jk+1)-pmfdq(jl,jk)+                        &
              &                 pmful(jl,jk+1)-pmful(jl,jk)-                        &
              &                 plude(jl,jk)-                                       &
              &                 (pdmfup(jl,jk)+pdmfdp(jl,jk))
            pqte_cnv(jl,jk) =  (grav/(paphp1(jl,jk+1)-paphp1(jl,jk)))*zqte_cnv(jl,jk)
            
            zxteci(jl,jk)   =  MERGE(0.0_wp,plude(jl,jk),llo1)
            pxteci(jl,jk)   =  MAX  (0.0_wp,(grav/(paphp1(jl,jk+1)-paphp1(jl,jk)))*zxteci(jl,jk))
            
            zxtecl(jl,jk)   =  MERGE(plude(jl,jk),0.0_wp,llo1)
            pxtecl(jl,jk)   =  MAX  (0.0_wp,(grav/(paphp1(jl,jk+1)-paphp1(jl,jk)))*zxtecl(jl,jk))
          ENDIF
        END DO
        !
          DO jt=1,ktrac
              DO jl=1,kproma
                IF(ldcum(jl)) THEN
                  pxtte_cnv(jl,jk,jt) = (grav/(paphp1(jl,jk+1)-paphp1(jl,jk)))      &
                    &                  *(pmfuxt(jl,jk+1,jt)-pmfuxt(jl,jk,jt)        &
                    &                   +pmfdxt(jl,jk+1,jt)-pmfdxt(jl,jk,jt))
                ENDIF
              END DO
          END DO
      !
      ELSE
        DO jl=1,kproma
          IF(ldcum(jl)) THEN
            llo1=(pten(jl,jk)-tmelt).GT.0._wp
            zalv=MERGE(alv,als,llo1)
            
            pq_cnv(jl,jk)   = -(pmfus(jl,jk)+pmfds(jl,jk)+alf*pdpmel(jl,jk)-        &
              &                 palvsh(jl,jk)*pmful(jl,jk)-                         &
              &                 zalv*(pdmfup(jl,jk)+pdmfdp(jl,jk)+plude(jl,jk)))
            
            zqte_cnv(jl,jk) = -(pmfuq(jl,jk)+pmfdq(jl,jk)+plude(jl,jk)+             &
              &                 (pmful(jl,jk)+pdmfup(jl,jk)+pdmfdp(jl,jk)))
            pqte_cnv(jl,jk) =  (grav/(paphp1(jl,jk+1)-paphp1(jl,jk)))*zqte_cnv(jl,jk)
            
            zxteci(jl,jk)   =  MERGE(0.0_wp,plude(jl,jk),llo1)
            pxteci(jl,jk)   =  MAX  (0.0_wp,(grav/(paphp1(jl,jk+1)-paphp1(jl,jk)))*zxteci(jl,jk))
            
            zxtecl(jl,jk)   =  MERGE(plude(jl,jk),0.0_wp,llo1)
            pxtecl(jl,jk)   =  MAX  (0.0_wp,(grav/(paphp1(jl,jk+1)-paphp1(jl,jk)))*zxtecl(jl,jk))
          END IF
        END DO
        !
          DO jt=1,ktrac
              DO jl=1,kproma
                IF(ldcum(jl)) THEN
                  pxtte_cnv(jl,jk,jt) =-(grav/(paphp1(jl,jk+1)-paphp1(jl,jk)))      &
                    &                  *(pmfuxt(jl,jk,jt)+pmfdxt(jl,jk,jt))
                ENDIF
              END DO
          END DO
      !
      END IF
      !
    END DO
    !
    !---------------------------------------------------------------------
    !
    !      3.          Update surface fields
    !                  ---------------------
    !
    ! do we need to account for the surface, or for the top 2 layers?
    DO jk=ktopm2,klev
      DO jl=1,kproma
        ! water vapor
        pcon_iqte(jl)=pcon_iqte(jl)+zqte_cnv(jl,jk)
        ! detrained liquid water
        pcon_dtrl(jl)=pcon_dtrl(jl)+zxtecl(jl,jk)
        ! detrained ice
        pcon_dtri(jl)=pcon_dtri(jl)+zxteci(jl,jk)
      END DO
    END DO
    !
  END SUBROUTINE cudtdq
  !>
  !!
  SUBROUTINE cududv(   kproma,   kbdim,    klev,     klevp1,                         &
    &        ktopm2,   ktype,    kcbot,    paphp1,   ldcum,                          &
    &        puen,     pven,                                                         &
    &        pvom_cnv, pvol_cnv,                                                     &
    &        puu,      pud,      pvu,      pvd,                                      &
    &        pmfu,     pmfd)
    !
    INTEGER, INTENT (IN) :: kproma, kbdim, klev, klevp1, ktopm2
    REAL(wp):: puen(kbdim,klev),        pven(kbdim,klev),                            &
      &        pvom_cnv(kbdim,klev),    pvol_cnv(kbdim,klev),                        &
      &        paphp1(kbdim,klevp1)
    REAL(wp):: puu(kbdim,klev),         pud(kbdim,klev),                             &
      &        pvu(kbdim,klev),         pvd(kbdim,klev),                             &
      &        pmfu(kbdim,klev),        pmfd(kbdim,klev)
    INTEGER :: ktype(kbdim),            kcbot(kbdim)
    LOGICAL :: ldcum(kbdim)
    !
    REAL(wp):: zmfuu(kbdim,klev),       zmfdu(kbdim,klev),                           &
      &        zmfuv(kbdim,klev),       zmfdv(kbdim,klev)
    INTEGER :: jl, jk, ik, ikb
    REAL(wp):: zzp, zdudt, zdvdt

    pvom_cnv(:,:) = 0._wp
    pvol_cnv(:,:) = 0._wp
    !
    !----------------------------------------------------------------------
    !
    !*    1.0          Calculate fluxes and update u and v tendencies
    !                  ----------------------------------------------
    !
    IF(ktopm2.EQ.1) THEN
      DO jk=2,klev
        ik=jk-1
        DO jl=1,kproma
          IF(ldcum(jl)) THEN
            zmfuu(jl,jk)=pmfu(jl,jk)*(puu(jl,jk)-puen(jl,ik))
            zmfuv(jl,jk)=pmfu(jl,jk)*(pvu(jl,jk)-pven(jl,ik))
            zmfdu(jl,jk)=pmfd(jl,jk)*(pud(jl,jk)-puen(jl,ik))
            zmfdv(jl,jk)=pmfd(jl,jk)*(pvd(jl,jk)-pven(jl,ik))
          END IF
        END DO
      END DO
      DO jl=1,kproma
        IF(ldcum(jl)) THEN
          zmfuu(jl,1)=zmfuu(jl,2)
          zmfuv(jl,1)=zmfuv(jl,2)
          zmfdu(jl,1)=zmfdu(jl,2)
          zmfdv(jl,1)=zmfdv(jl,2)
        END IF
      END DO
    ELSE
      DO jk=ktopm2,klev
        ik=jk-1
        DO jl=1,kproma
          IF(ldcum(jl)) THEN
            zmfuu(jl,jk)=pmfu(jl,jk)*(puu(jl,jk)-puen(jl,ik))
            zmfuv(jl,jk)=pmfu(jl,jk)*(pvu(jl,jk)-pven(jl,ik))
            zmfdu(jl,jk)=pmfd(jl,jk)*(pud(jl,jk)-puen(jl,ik))
            zmfdv(jl,jk)=pmfd(jl,jk)*(pvd(jl,jk)-pven(jl,ik))
          END IF
        END DO
      END DO
    END IF
    DO jk=ktopm2,klev
!DIR$ IVDEP
!OCL NOVREC
      DO jl=1,kproma
        IF(ldcum(jl).AND.jk.GT.kcbot(jl)) THEN
          ikb=kcbot(jl)
          zzp=((paphp1(jl,klevp1)-paphp1(jl,jk))/(paphp1(jl,klevp1)-paphp1(jl,ikb)))
          zzp=MERGE(zzp**2,zzp,ktype(jl).EQ.3)
          zmfuu(jl,jk)=zmfuu(jl,ikb)*zzp
          zmfuv(jl,jk)=zmfuv(jl,ikb)*zzp
          zmfdu(jl,jk)=zmfdu(jl,ikb)*zzp
          zmfdv(jl,jk)=zmfdv(jl,ikb)*zzp
        END IF
      END DO
    END DO
    !
    DO jk=ktopm2,klev
      !
      IF(jk.LT.klev) THEN
        DO jl=1,kproma
          IF(ldcum(jl)) THEN
            zdudt=(grav/(paphp1(jl,jk+1)-paphp1(jl,jk)))*                            &
              &         (zmfuu(jl,jk+1)-zmfuu(jl,jk)+zmfdu(jl,jk+1)-zmfdu(jl,jk))
            zdvdt=(grav/(paphp1(jl,jk+1)-paphp1(jl,jk)))*                            &
              &         (zmfuv(jl,jk+1)-zmfuv(jl,jk)+zmfdv(jl,jk+1)-zmfdv(jl,jk))
            pvom_cnv(jl,jk)=zdudt
            pvol_cnv(jl,jk)=zdvdt
          END IF
        END DO
        !
      ELSE
        DO jl=1,kproma
          IF(ldcum(jl)) THEN
            zdudt=-(grav/(paphp1(jl,jk+1)-paphp1(jl,jk)))*                           &
              &            (zmfuu(jl,jk)+zmfdu(jl,jk))
            zdvdt=-(grav/(paphp1(jl,jk+1)-paphp1(jl,jk)))*                           &
              &            (zmfuv(jl,jk)+zmfdv(jl,jk))
            pvom_cnv(jl,jk)=zdudt
            pvol_cnv(jl,jk)=zdvdt
          END IF
        END DO
      END IF
      !
    END DO
  END SUBROUTINE cududv

END MODULE mo_cufluxdts
