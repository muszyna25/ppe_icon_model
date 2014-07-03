!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
MODULE mo_cufluxdts

USE mo_kind,               ONLY: wp
USE mo_physical_constants, ONLY: grav, alv, als, alf, tmelt, rd, cpd, vtmpc2
!USE mo_physc2,             ONLY: cevapcu
!USE mo_time_control,       ONLY: time_step_len
!USE mo_submodel,           ONLY: lanysubmodel, lham    ! ### explicit submodel dependency (zcucov)
!USE mo_submodel_interface, ONLY: cuflx_subm
!USE mo_tracdef,            ONLY : trlist
!USE mo_time_control, ONLY : delta_time
!USE mo_vphysc,       ONLY: set_vphysc_var
!++ for cfmip diagnostics
!USE mo_memory_cfdiag, ONLY : locfdiag, mc, mcu, mcd, imc, smc, dmc
!USE mo_geoloc,        ONLY : gboxarea
!USE mo_cosp_offline, ONLY: locospoffl, cospoffl_ccrain, cospoffl_ccsnow
!++-end-for cfmip diagnostics----------------------------  
!
IMPLICIT NONE

PRIVATE

PUBLIC :: cuflx, cudtdq, cududv

CONTAINS

!++mgs: zcucov and zdpevap now 1d vectors
SUBROUTINE cuflx(    kproma, kbdim, klev, klevp1,                      &
           pqen,     pqsen,    ptenh,    pqenh,                        &
           ktrac,                                                      &
!---Included for scavenging in wetdep_interface (Philip Stier, 28/03/01, UL, 2803.07):-------
           krow,                                                       &
           pxtte,    pxtu,     ptu,                                    &
           pmwc,     pmrateprecip,                                     &
!---End Included for scavenging-----------------------------------------
           pxtenh,   pmfuxt,   pmfdxt,                                 &
           paphp1,   pgeoh,                                            &
           kcbot,    kctop,    kdtop,                                  &
           ktype,    lddraf,   ldcum,                                  &
           pmfu,     pmfd,     pmfus,    pmfds,                        &
           pmfuq,    pmfdq,    pmful,                                  &
           pdmfup,   pdmfdp,   prfl,     prain,                        &
           pcpcu,                                                      &
           pten,     psfl,     pdpmel,   ktopm2,  cevapcu, time_step_len)
!
!          M.TIEDTKE         E.C.M.W.F.     7/86 MODIF. 12/89
!
!          PURPOSE
!          -------
!
!          THIS ROUTINE DOES THE FINAL CALCULATION OF CONVECTIVE
!          FLUXES IN THE CLOUD LAYER AND IN THE SUBCLOUD LAYER
!
!          INTERFACE
!          ---------
!          THIS ROUTINE IS CALLED FROM *CUMASTR*.
!
!          EXTERNALS
!          ---------
!          NONE
!
!
INTEGER, INTENT (IN) :: kproma, kbdim, klev, klevp1, ktrac
INTEGER, INTENT (OUT):: ktopm2
REAL(wp),INTENT(IN)  :: cevapcu(klev)
REAL(wp), INTENT(in) :: time_step_len 
!
REAL(wp):: pqen(kbdim,klev),        pqsen(kbdim,klev),                 &
           ptenh(kbdim,klev),       pqenh(kbdim,klev),                 &
           paphp1(kbdim,klevp1),    pgeoh(kbdim,klev)
!
REAL(wp):: pmfu(kbdim,klev),        pmfd(kbdim,klev),                  &
           pmfus(kbdim,klev),       pmfds(kbdim,klev),                 &
           pmfuq(kbdim,klev),       pmfdq(kbdim,klev),                 &
           pdmfup(kbdim,klev),      pdmfdp(kbdim,klev),                &
           pmful(kbdim,klev),                                          &
           prfl(kbdim),             prain(kbdim)
REAL(wp):: pten(kbdim,klev),        pdpmel(kbdim,klev),                &
           psfl(kbdim)
REAL(wp):: pcpcu(kbdim,klev)
INTEGER :: kcbot(kbdim),            kctop(kbdim),                      &
           kdtop(kbdim),            ktype(kbdim)
LOGICAL :: lddraf(kbdim),           ldcum(kbdim)
REAL(wp):: pxtenh(kbdim,klev,ktrac),                                   &
           pmfuxt(kbdim,klev,ktrac),pmfdxt(kbdim,klev,ktrac)
!
!
REAL(wp) :: afac(kbdim) ! for cfmip diagnostics
!
INTEGER :: jl, jk, jt, ikb
!++mgs
!!REAL(wp):: zcons1, zcons2, zcucov, ztmelp2, zzp, zfac, zsnmlt          &
!!         , zrfl, zrnew, zrmin, zrfln, zdrfl, zrsum, zdpevap,zwu
!!REAL(wp):: zpsubcl(kbdim)
REAL(wp):: zcons1, zcons2,         ztmelp2, zzp, zfac, zsnmlt          &
         , zrfl, zrnew, zrmin, zrfln, zdrfl, zrsum,         zwu
REAL(wp):: zpsubcl(kbdim), zpsubcl_sav(kbdim), zcucov(kbdim), zdpevap(kbdim)
!--mgs
!---Included for scavenging in wetdep_interface (Philip Stier, 28/03/01):-------
!
INTEGER, INTENT(IN) :: krow

REAL(wp):: pxtte(kbdim,klev,ktrac), pxtu(kbdim,klev,ktrac),            &
           ptu(kbdim,klev)

REAL(wp):: zmlwc(kbdim,klev),       zmiwc(kbdim,klev),                 &
           zmratepr(kbdim,klev),    zmrateps(kbdim,klev),              &
           zfrain(kbdim,klev),      zfsnow(kbdim,klev),                &
           zdpg(kbdim,klev),        zfevapr(kbdim,klev),               &
           zfsubls(kbdim,klev),     zaclc(kbdim,klev),                 &
           zmsnowacl(kbdim,klev),   zrhou(kbdim,klev)

REAL(wp):: pmwc(kbdim,klev),        pmrateprecip(kbdim,klev)

REAL(wp):: ztc,         zzfac

REAL(wp):: zcaa, &   !(Constants for partitioning of cloud water 
           zcab      ! into liquid and solid part)

REAL(wp):: zalpha    ! Fraction of cloud water in liquid phase
                        ! =>  (1-zalpha)      "   in solid  phase

!---End Included for scavenging-----------------------------------------
!
!*             SPECIFY CONSTANTS
!

  zcons1=cpd/(alf*grav*time_step_len)
  zcons2=1._wp/(grav*time_step_len)
!!mgs!!  zcucov=0.05_wp            !! ++mgs: now 1-d array
  ztmelp2=tmelt+2._wp
!
!---Included for scavenging in wetdep_interface (Philip Stier, 28/03/01):-------
!
!IF(lanysubmodel) THEN
!
!  zcaa = 0.0059_wp   ! (Constants for partitioning of cloud water 
!  zcab = 0.003102_wp !  into liquid and solid part)
!
!  zmlwc(1:kproma,:)     = 0._wp
!  zmiwc(1:kproma,:)     = 0._wp
!  zmratepr(1:kproma,:)  = 0._wp
!  zmrateps(1:kproma,:)  = 0._wp
!  zfrain(1:kproma,:)    = 0._wp
!  zfsnow(1:kproma,:)    = 0._wp
!  zdpg(1:kproma,:)      = 0._wp
!  zfevapr(1:kproma,:)   = 0._wp
!  zfsubls(1:kproma,:)   = 0._wp
!  zmsnowacl(1:kproma,:) = 0._wp
!
!  zwu            = 2.0_wp
!
!  !--- Set cloud cover to 1 below cloud top and bottom:
!
!  DO jk=1, klev
!     DO jl=1, kproma
!        IF (jk>=kctop(jl) .AND. jk<=kcbot(jl)) THEN
!           zaclc(jl,jk) = 1._wp
!        ELSE
!           zaclc(jl,jk) = 0._wp
!        END IF
!     END DO
!  END DO
!END IF
!
!---End Included for scavenging-----------------------------------------
!
!
!
!*    1.0          DETERMINE FINAL CONVECTIVE FLUXES
!                  ---------------------------------
!
!  itop=klev
  DO 110 jl=1,kproma
!     itop=MIN(itop,kctop(jl))
     IF(.NOT.ldcum(jl).OR.kdtop(jl).LT.kctop(jl)) lddraf(jl)=.FALSE.
     IF(.NOT.ldcum(jl)) ktype(jl)=0
110 END DO
  ktopm2=1
  DO 120 jk=ktopm2,klev
!DIR$ IVDEP
!OCL NOVREC
     DO 115 jl=1,kproma
        IF(ldcum(jl).AND.jk.GE.kctop(jl)-1) THEN
           pmfus(jl,jk)=pmfus(jl,jk)-pmfu(jl,jk)*                      &
                               (pcpcu(jl,jk)*ptenh(jl,jk)+pgeoh(jl,jk))
           pmfuq(jl,jk)=pmfuq(jl,jk)-pmfu(jl,jk)*pqenh(jl,jk)
           IF(lddraf(jl).AND.jk.GE.kdtop(jl)) THEN
              pmfds(jl,jk)=pmfds(jl,jk)-pmfd(jl,jk)*                   &
                               (pcpcu(jl,jk)*ptenh(jl,jk)+pgeoh(jl,jk))
              pmfdq(jl,jk)=pmfdq(jl,jk)-pmfd(jl,jk)*pqenh(jl,jk)
           ELSE
              pmfd(jl,jk)=0._wp
              pmfds(jl,jk)=0._wp
              pmfdq(jl,jk)=0._wp
              pdmfdp(jl,jk-1)=0._wp
           END IF
        END IF
115  END DO
!
     DO 1154 jt=1,ktrac
        DO 1152 jl=1,kproma
           IF(ldcum(jl).AND.jk.GE.kctop(jl)-1) THEN
              pmfuxt(jl,jk,jt)=pmfuxt(jl,jk,jt)                        &
                                     -pmfu(jl,jk)*pxtenh(jl,jk,jt)
              IF(lddraf(jl).AND.jk.GE.kdtop(jl)) THEN
                 pmfdxt(jl,jk,jt)=pmfdxt(jl,jk,jt)                     &
                                     -pmfd(jl,jk)*pxtenh(jl,jk,jt)
              ELSE
                 pmfdxt(jl,jk,jt)=0._wp
              ENDIF
           ELSE
              pmfuxt(jl,jk,jt)=0._wp
              pmfdxt(jl,jk,jt)=0._wp
           ENDIF
1152    END DO
1154 END DO
!
120 END DO
  DO 130 jk=ktopm2,klev
!DIR$ IVDEP
!OCL NOVREC
     DO 125 jl=1,kproma
        IF(ldcum(jl).AND.jk.GT.kcbot(jl)) THEN
           ikb=kcbot(jl)
           zzp=((paphp1(jl,klevp1)-paphp1(jl,jk))/                     &
                         (paphp1(jl,klevp1)-paphp1(jl,ikb)))
           zzp=MERGE(zzp**2,zzp,ktype(jl).EQ.3)
           pmfu(jl,jk)=pmfu(jl,ikb)*zzp
           pmfus(jl,jk)=pmfus(jl,ikb)*zzp
           pmfuq(jl,jk)=pmfuq(jl,ikb)*zzp
           pmful(jl,jk)=pmful(jl,ikb)*zzp
        END IF
125  END DO
!
     DO 1254 jt=1,ktrac
!DIR$ IVDEP
!OCL NOVREC
        DO 1252 jl=1,kproma
           IF(ldcum(jl).AND.jk.GT.kcbot(jl)) THEN
              ikb=kcbot(jl)
              zzp=(paphp1(jl,klevp1)-paphp1(jl,jk))/                   &
                            (paphp1(jl,klevp1)-paphp1(jl,ikb))
              zzp=MERGE(zzp**2,zzp,ktype(jl).EQ.3)
              pmfuxt(jl,jk,jt)=pmfuxt(jl,ikb,jt)*zzp
           ENDIF
1252    END DO
1254 END DO
!
130 END DO
!
!
!*    2.            CALCULATE RAIN/SNOW FALL RATES
!*                  CALCULATE MELTING OF SNOW
!*                  CALCULATE EVAPORATION OF PRECIP
!                   -------------------------------
!
  DO 210 jl=1,kproma
     prfl(jl)=0._wp
     psfl(jl)=0._wp
     prain(jl)=0._wp
210 END DO
!++jsr interface for scavenging
!IF (lanysubmodel) THEN
!   DO jk=ktopm2,klev
!!---Included for scavenging in wetdep_interface (Philip Stier, 28/03/01):-------
!      zdpg(1:kproma,jk)=(paphp1(1:kproma,jk+1)-paphp1(1:kproma,jk))/grav
!!---End Included for scavenging-----------------------------------------
!   END DO
!   DO jk=ktopm2,klev
!      DO jl=1,kproma
!         IF (ldcum(jl)) THEN
!            ztc=ptu(jl,jk)-tmelt
!            zzfac=MERGE(1._wp,0._wp,ztc<0._wp)
!            zalpha=(1._wp-zzfac)+zzfac*(zcaa+(1._wp-zcaa)*exp(-zcab*ztc**2))
!      
!            zmlwc(jl,jk)=zalpha*pmwc(jl,jk)
!            zmiwc(jl,jk)=(1._wp-zalpha)*pmwc(jl,jk)
!            
!            zmratepr(jl,jk)=zalpha*pmrateprecip(jl,jk)
!            zmrateps(jl,jk)=(1._wp-zalpha)*pmrateprecip(jl,jk)
!         END IF
!      END DO
!   END DO
!END IF
!!---End Included for scavenging-----------------------------------------
!--jsr interface for scavenging
      
  DO 220 jk=ktopm2,klev
     DO 215 jl=1,kproma
        IF(ldcum(jl)) THEN
           prain(jl)=prain(jl)+pdmfup(jl,jk)
           IF(pten(jl,jk).GT.tmelt) THEN
              prfl(jl)=prfl(jl)+pdmfup(jl,jk)+pdmfdp(jl,jk)
              IF(psfl(jl).GT.0._wp.AND.pten(jl,jk).GT.ztmelp2) THEN
                 zfac=zcons1*(1._wp+vtmpc2*pqen(jl,jk))                &
                             *(paphp1(jl,jk+1)-paphp1(jl,jk))
                 zsnmlt=MIN(psfl(jl),zfac*(pten(jl,jk)-ztmelp2))
                 pdpmel(jl,jk)=zsnmlt
                 psfl(jl)=psfl(jl)-zsnmlt
                 prfl(jl)=prfl(jl)+zsnmlt
              END IF
           ELSE
              psfl(jl)=psfl(jl)+pdmfup(jl,jk)+pdmfdp(jl,jk)
           END IF
        END IF
215  END DO
!     IF (lanysubmodel) THEN
!!---Included for scavenging in wetdep_interface (Philip Stier, 28/03/01):-------
!        zfrain(1:kproma,jk)=prfl(1:kproma)
!        zfsnow(1:kproma,jk)=psfl(1:kproma)
!!---End Included for scavenging-----------------------------------------      
!     END IF
220 END DO
  DO 230 jl=1,kproma
     prfl(jl)=MAX(prfl(jl),0._wp)
     psfl(jl)=MAX(psfl(jl),0._wp)
     zpsubcl(jl)=prfl(jl)+psfl(jl)
230 END DO
!++jsr scavenging interface
!++mgs ### new code
!!IF (lanysubmodel) THEN
!!   IF (lham) THEN
!!      DO jk=ktopm2,klev
!!         DO jl=1,kproma
!!            IF(ldcum(jl).AND.jk.GE.kcbot(jl).AND.zpsubcl(jl).GT.1.e-20_wp) THEN
!!               zrhou(jl,jk)=(paphp1(jl,jk+1)+paphp1(jl,jk))* &
!!                            0.5_wp/(ptu(jl,jk)*rd)
!!               zcucov=pmfu(jl,jk)/(zwu*zrhou(jl,jk))
!!               zrfl=zpsubcl(jl)
!!               zrnew=(MAX(0._wp,SQRT(zrfl/zcucov)-                         &
!!                      cevapcu(jk)*(paphp1(jl,jk+1)-paphp1(jl,jk))*   &
!!                      MAX(0._wp,pqsen(jl,jk)-pqen(jl,jk))))**2*zcucov
!!               zrmin=zrfl-zcucov*MAX(0._wp,0.8_wp*pqsen(jl,jk)-pqen(jl,jk))&
!!                     *zcons2*(paphp1(jl,jk+1)-paphp1(jl,jk))
!!               zrnew=MAX(zrnew,zrmin)
!!               zrfln=MAX(zrnew,0._wp)
!!               zdrfl=MIN(0._wp,zrfln-zrfl)
!!               pdmfup(jl,jk)=pdmfup(jl,jk)+zdrfl
!!!---Included for scavenging in wetdep_interface (Philip Stier, 28/03/01):-------
!!               zdpevap=-zdrfl
!!               zfevapr(jl,jk)=zdpevap*prfl(jl)/zpsubcl(jl)
!!               zfsubls(jl,jk)=zdpevap*psfl(jl)/zpsubcl(jl)
!!!---End Included for scavenging-----------------------------------------  
!!               zpsubcl(jl)=zrfln
!!            END IF
!!         END DO
!!      END DO
!!   ELSE
!!      DO jk=ktopm2,klev
!!         DO jl=1,kproma
!!            IF(ldcum(jl).AND.jk.GE.kcbot(jl).AND.zpsubcl(jl).GT.1.e-20_wp) THEN
!!               zrhou(jl,jk)=(paphp1(jl,jk+1)+paphp1(jl,jk))* &
!!                            0.5_wp/(ptu(jl,jk)*rd)
!!               zrfl=zpsubcl(jl)
!!               zrnew=(MAX(0._wp,SQRT(zrfl/zcucov)-                         &
!!                      cevapcu(jk)*(paphp1(jl,jk+1)-paphp1(jl,jk))*   &
!!                      MAX(0._wp,pqsen(jl,jk)-pqen(jl,jk))))**2*zcucov
!!               zrmin=zrfl-zcucov*MAX(0._wp,0.8_wp*pqsen(jl,jk)-pqen(jl,jk))&
!!                     *zcons2*(paphp1(jl,jk+1)-paphp1(jl,jk))
!!               zrnew=MAX(zrnew,zrmin)
!!               zrfln=MAX(zrnew,0._wp)
!!               zdrfl=MIN(0._wp,zrfln-zrfl)
!!               pdmfup(jl,jk)=pdmfup(jl,jk)+zdrfl
!!!---Included for scavenging in wetdep_interface (Philip Stier, 28/03/01):-------
!!               zdpevap=-zdrfl
!!               zfevapr(jl,jk)=zdpevap*prfl(jl)/zpsubcl(jl)
!!               zfsubls(jl,jk)=zdpevap*psfl(jl)/zpsubcl(jl)
!!!---End Included for scavenging-----------------------------------------  
!!               zpsubcl(jl)=zrfln
!!            END IF
!!         END DO
!!      END DO
!!   END IF
!!ELSE
!!  DO 240 jk=ktopm2,klev
!!     DO 235 jl=1,kproma
!!       IF(ldcum(jl).AND.jk.GE.kcbot(jl).AND.zpsubcl(jl).GT.1.e-20_wp)  &
!!                                                                  THEN
!!           zrfl=zpsubcl(jl)
!!           zrnew=(MAX(0._wp,SQRT(zrfl/zcucov)-                         &
!!                        cevapcu(jk)*(paphp1(jl,jk+1)-paphp1(jl,jk))*   &
!!                        MAX(0._wp,pqsen(jl,jk)-pqen(jl,jk))))**2*zcucov
!!           zrmin=zrfl-zcucov*MAX(0._wp,0.8_wp*pqsen(jl,jk)-pqen(jl,jk))&
!!                        *zcons2*(paphp1(jl,jk+1)-paphp1(jl,jk))
!!           zrnew=MAX(zrnew,zrmin)
!!           zrfln=MAX(zrnew,0._wp)
!!           zdrfl=MIN(0._wp,zrfln-zrfl)
!!           pdmfup(jl,jk)=pdmfup(jl,jk)+zdrfl
!!           zpsubcl(jl)=zrfln
!!       END IF
!!235  END DO
!!240 END DO
!!END IF
!  IF (lanysubmodel) THEN
!    DO jk=ktopm2,klev
!      zrhou(1:kproma,jk)=(paphp1(1:kproma,jk+1)+paphp1(1:kproma,jk))* &
!                          0.5_wp/(ptu(1:kproma,jk)*rd)
!    END DO
!  END IF

  DO jk=ktopm2,klev

    zdpevap(1:kproma) = 0._wp
!    IF (lanysubmodel) zpsubcl_sav(1:kproma) = zpsubcl(1:kproma)
!    ! ### value of zcucov depends on lham => explicit submodel dependence !
!    IF (lham) THEN
!      zcucov(1:kproma) = pmfu(1:kproma,jk)/(zwu*zrhou(1:kproma,jk))
!    ELSE
      zcucov(1:kproma) = 0.05_wp
!    END IF

    DO jl=1,kproma
      IF(ldcum(jl).AND.jk.GE.kcbot(jl).AND.zpsubcl(jl).GT.1.e-20_wp) THEN
        zrfl=zpsubcl(jl)
        zrnew=(MAX(0._wp,SQRT(zrfl/zcucov(jl))-                           &
               cevapcu(jk)*(paphp1(jl,jk+1)-paphp1(jl,jk))*               &
               MAX(0._wp,pqsen(jl,jk)-pqen(jl,jk))))**2*zcucov(jl)
        zrmin=zrfl-zcucov(jl)*MAX(0._wp,0.8_wp*pqsen(jl,jk)-pqen(jl,jk))  &
              *zcons2*(paphp1(jl,jk+1)-paphp1(jl,jk))
        zrnew=MAX(zrnew,zrmin)
        zrfln=MAX(zrnew,0._wp)
        zdrfl=MIN(0._wp,zrfln-zrfl)
        zdpevap(jl) = -zdrfl
        pdmfup(jl,jk)=pdmfup(jl,jk)+zdrfl
        zpsubcl(jl)=zrfln
      END IF
    END DO

!---Included for scavenging in wetdep_interface (Philip Stier, 28/03/01):-------
!    IF (lanysubmodel) THEN
!      DO jl=1,kproma
!        IF (zpsubcl_sav(jl) > 1.e-20_wp) THEN
!          zfevapr(jl,jk)=zdpevap(jl)*prfl(jl)/zpsubcl_sav(jl)
!          zfsubls(jl,jk)=zdpevap(jl)*psfl(jl)/zpsubcl_sav(jl)
!        ELSE
!          zfevapr(jl,jk) = 0._wp
!          zfsubls(jl,jk) = 0._wp
!        END IF
!      END DO
!    END IF
!---End Included for scavenging-----------------------------------------
  END DO

!--mgs ### new code

!  IF (lanysubmodel) THEN
!    CALL cuflx_subm(kbdim,  kproma,    klev,     ktopm2,   & ! dimensions
!                    krow,                                  & ! longitude (kproma-block) index
!                    pxtenh, pxtu, zrhou,                   & ! tracers
!                    pmfu,   pmfuxt,                        & ! convective fluxes and corresp. mmr
!                    zmlwc,  zmiwc,     zmratepr, zmrateps, & ! cloud properties
!                    zfrain, zfsnow,    zfevapr,  zfsubls,  & !   "       "
!                    zaclc,  zmsnowacl,                     & !   "       "
!                    ptu,    zdpg,                          & ! thermodynamic quantities
!                    pxtte                                  ) 
!  END IF

!!baustelle!! (?)
  DO 250 jl=1,kproma
     zrsum=prfl(jl)+psfl(jl)
     zdpevap(jl)=zpsubcl(jl)-zrsum
     prfl(jl)=prfl(jl)+zdpevap(jl)*prfl(jl)*(1._wp/MAX(1.e-20_wp,zrsum))
     psfl(jl)=psfl(jl)+zdpevap(jl)*psfl(jl)*(1._wp/MAX(1.e-20_wp,zrsum))
250 END DO
!
!  IF ( locfdiag ) THEN
!    DO jk=ktopm2, klev
!      IF (lham) THEN
!        zcucov(1:kproma) = pmfu(1:kproma,jk)/(zwu*zrhou(1:kproma,jk))
!      ELSE
!        zcucov(1:kproma) = 0.05_wp
!      END IF
!
!       afac(1:kproma)=zcucov(1:kproma)/gboxarea(1:kproma)
!
!        imc(1:kproma,jk,krow)  = afac(1:kproma) * (pmfu(1:kproma,jk) + pmfd(1:kproma,jk))  
!  
!        DO jl=1,kproma
!          IF ( ktype(jl).EQ.1 ) THEN
!             dmc (jl,jk,krow)  = dmc(jl,jk,krow) + imc(jl,jk,krow) * delta_time
!          END IF
!          IF ( ktype(jl).EQ.2 ) THEN
!             smc (jl,jk,krow)  = smc(jl,jk,krow) + imc(jl,jk,krow) * delta_time
!          END IF
!        END DO        
!       mc (1:kproma,jk,krow)  = mc(1:kproma,jk,krow) + imc(1:kproma,jk,krow) * delta_time
!       mcu (1:kproma,jk,krow) = mcu(1:kproma,jk,krow) + afac(1:kproma) * &
!                                (pmfu(1:kproma,jk)) * delta_time
!       mcd (1:kproma,jk,krow) = mcd(1:kproma,jk,krow) + afac(1:kproma) * &
!                                (pmfd(1:kproma,jk)) * delta_time
!       
!
!
!     END DO
!  END IF ! locfdiag

!  IF ( locospoffl ) THEN
!     DO jk=ktopm2,klev
!        DO jl=1,kproma
!           IF(ldcum(jl)) THEN
!              IF(pten(jl,jk).GT.tmelt) THEN
!                 cospoffl_ccrain(jl,jk,krow) = pdmfup(jl,jk)+pdmfdp(jl,jk)
!                 cospoffl_ccsnow(jl,jk,krow) = 0._wp
!              ELSE
!                 cospoffl_ccsnow(jl,jk,krow) = pdmfup(jl,jk)+pdmfdp(jl,jk)
!                 cospoffl_ccrain(jl,jk,krow) = 0._wp
!              END IF
!           END IF
!        END DO
!     END DO
!  ENDIF !locospoffl


  RETURN
END SUBROUTINE cuflx

SUBROUTINE cudtdq(kproma, kbdim, klev, klevp1, ktopm2, ldcum, ktrac,   &
                  krow,                                                &
                  paphp1,   pten,     ptte,     pqte,                  &
                  pxtte,    pxtec,    pmfuxt,   pmfdxt,                &
                  pmfus,    pmfds,    pmfuq,    pmfdq,                 &
                  pmful,    pdmfup,   pdmfdp,   plude,                 &
                  pdpmel,   prfl,     psfl,                            &
                  pcpen,    palvsh,   pqtec,    pqude,                 &
                  prsfc,    pssfc,                                     &
                  pch_concloud, pcon_dtrl, pcon_dtri,                  &
                  pxtecl,   pxteci, pcon_iqte,                         &
                  ptte_cnv, pqte_cnv, pxtte_cnv                        )
!
!
!**** *CUDTDQ* - UPDATES T AND Q TENDENCIES, PRECIPITATION RATES
!                DOES GLOBAL DIAGNOSTICS
!
!          M.TIEDTKE         E.C.M.W.F.     7/86 MODIF. 12/89
!
!**   INTERFACE.
!     ----------
!
!          *CUDTDQ* IS CALLED FROM *CUMASTR*
!
!
INTEGER, INTENT (IN) :: kproma, kbdim, klev, klevp1, ktopm2, ktrac
INTEGER, INTENT (IN) :: krow
! pch_concloud, pcon_dtrl, and pcon_dtri track the convective heating and 
! the detrained liquid and ice. All these are only for diag purposes
REAL(wp),INTENT(INOUT) :: pch_concloud(kbdim), pcon_dtrl(kbdim), pcon_dtri(kbdim)
REAL(wp),INTENT(INOUT) :: pcon_iqte(kbdim) ! integrated qv tendency
REAL(wp),INTENT(INOUT) :: ptte_cnv(kbdim,klev)                               
REAL(wp),INTENT(INOUT) :: pqte_cnv(kbdim,klev), pxtte_cnv(kbdim,klev,ktrac)
REAL(wp) :: zqte_cnv(kbdim,klev) 
LOGICAL  llo1
!
REAL(wp) :: ptte(kbdim,klev),        pqte(kbdim,klev),                 &
            pten(kbdim,klev),        paphp1(kbdim,klevp1),             &
            prsfc(kbdim),            pssfc(kbdim)
REAL(wp) :: pxtecl(kbdim,klev),      pxteci(kbdim,klev)
REAL(wp) :: pmfus(kbdim,klev),       pmfds(kbdim,klev),                &
            pmfuq(kbdim,klev),       pmfdq(kbdim,klev),                &
            pmful(kbdim,klev),       plude(kbdim,klev),                &
            pdmfup(kbdim,klev),      pdmfdp(kbdim,klev),               &
            pqtec(kbdim,klev),       pqude(kbdim,klev),                &
            pxtec(kbdim,klev),       prfl(kbdim)
REAL(wp) :: pdpmel(kbdim,klev),      psfl(kbdim)
REAL(wp) :: pcpen(kbdim,klev),       palvsh(kbdim,klev)
LOGICAL  :: ldcum(kbdim)
!
REAL(wp) :: zmelt(kbdim), zcpten(kbdim,klev) 
REAL(wp) :: zteci(kbdim,klev), ztecl(kbdim,klev)
REAL(wp) :: zsheat(kbdim)
REAL(wp) :: pxtte(kbdim,klev,ktrac), pmfuxt(kbdim,klev,ktrac),         &
            pmfdxt(kbdim,klev,ktrac)
!
REAL(wp) :: zrcpm ! reciprocal value of specific heat of moist air

INTEGER  :: jl, jk, jt
REAL(wp) :: zdiagt, zalv, zdtdt, zdqdt, zdxtdt
!
   ptte_cnv(1:kproma,:)   = 0._wp
   pqte_cnv(1:kproma,:)   = 0._wp
  pxtte_cnv(1:kproma,:,:) = 0._wp
!
!----------------------------------------------------------------------
!
!*    1.0          SPECIFY PARAMETERS
!                  ------------------
!
!  zdiagt=delta_time
!
!----------------------------------------------------------------------
!
!*    2.0          INCREMENTATION OF T AND Q TENDENCIES
!                  ------------------------------------
!
  DO 210 jl=1,kproma
     zmelt(jl)=0._wp
     zsheat(jl)=0._wp
210 END DO
!
  pch_concloud(1:kproma)=0._wp
  pcon_dtrl(1:kproma)=0._wp
  pcon_dtri(1:kproma)=0._wp
  pcon_iqte(1:kproma)=0._wp
  zcpten(1:kproma,ktopm2:klev)=0._wp
  zqte_cnv(1:kproma,ktopm2:klev)=0._wp
! why aren't these zeroed for all levels?!
  zteci(1:kproma,ktopm2:klev)=0._wp
  ztecl(1:kproma,ktopm2:klev)=0._wp

  DO 250 jk=ktopm2,klev
!
     IF(jk.LT.klev) THEN
        DO 220 jl=1,kproma
           IF(ldcum(jl)) THEN
              llo1=(pten(jl,jk)-tmelt).GT.0._wp
              zalv=MERGE(alv,als,llo1)
              zrcpm=1._wp/pcpen(jl,jk)
              zdtdt=(grav/(paphp1(jl,jk+1)-paphp1(jl,jk)))*zrcpm*      &
                                  (pmfus(jl,jk+1)-pmfus(jl,jk)+        &
                                   pmfds(jl,jk+1)-pmfds(jl,jk)-        &
                                   alf*pdpmel(jl,jk)-                  &
                                   palvsh(jl,jk+1)*pmful(jl,jk+1)+     &
                                   palvsh(jl,jk)*pmful(jl,jk)+         &
                                   zalv*(plude(jl,jk)+pdmfup(jl,jk)+pdmfdp(jl,jk)))
              ptte(jl,jk)=ptte(jl,jk)+zdtdt
              zcpten(jl,jk)=zdtdt*pcpen(jl,jk)
              ptte_cnv(jl,jk)=zdtdt
              zdqdt=(grav/(paphp1(jl,jk+1)-paphp1(jl,jk)))*               &
                                  (pmfuq(jl,jk+1)-pmfuq(jl,jk)+        &
                                   pmfdq(jl,jk+1)-pmfdq(jl,jk)+        &
                                   pmful(jl,jk+1)-pmful(jl,jk)-        &
                                   plude(jl,jk)-                       &
                                  (pdmfup(jl,jk)+pdmfdp(jl,jk)))
              pqte(jl,jk)=pqte(jl,jk)+zdqdt
              zqte_cnv(jl,jk)=zdqdt
              pxtec(jl,jk)=(grav/(paphp1(jl,jk+1)-                        &
                            paphp1(jl,jk)))*plude(jl,jk)
              pxteci(jl,jk)=MERGE(0.0_wp,pxtec(jl,jk),llo1)
              zteci(jl,jk)=pxteci(jl,jk)
              pxtecl(jl,jk)=MERGE(pxtec(jl,jk),0.0_wp,llo1)
              ztecl(jl,jk)=pxtecl(jl,jk)
              pqtec(jl,jk)=(grav/(paphp1(jl,jk+1)-                        &
                            paphp1(jl,jk)))*pqude(jl,jk)
              zsheat(jl)=zsheat(jl)+zalv*(pdmfup(jl,jk)+pdmfdp(jl,jk))
              zmelt(jl)=zmelt(jl)+pdpmel(jl,jk)
           ENDIF
220     END DO
!
!        IF (trlist% anyconv /= 0) THEN
!           DO 2204 jt=1,ktrac
!              IF (trlist% ti(jt)% nconv == 1) THEN
!                DO 2202 jl=1,kproma
!                   IF(ldcum(jl)) THEN
!                     zdxtdt=(grav/(paphp1(jl,jk+1)-paphp1(jl,jk)))        &
!                                 *(pmfuxt(jl,jk+1,jt)-pmfuxt(jl,jk,jt) &
!                                  +pmfdxt(jl,jk+1,jt)-pmfdxt(jl,jk,jt))
!                     pxtte(jl,jk,jt)=pxtte(jl,jk,jt)+zdxtdt
!                   ENDIF
!2202            END DO
!              ENDIF
!2204       END DO
!        ENDIF
!
!
     ELSE
        DO 230 jl=1,kproma
           IF(ldcum(jl)) THEN
              llo1=(pten(jl,jk)-tmelt).GT.0._wp
              zalv=MERGE(alv,als,llo1)
              zrcpm=1._wp/pcpen(jl,jk)
              zdtdt=-(grav/(paphp1(jl,jk+1)-paphp1(jl,jk)))*zrcpm*     &
                     (pmfus(jl,jk)+pmfds(jl,jk)+alf*pdpmel(jl,jk)-     &
                      palvsh(jl,jk)*pmful(jl,jk)-                      &
                      zalv*(pdmfup(jl,jk)+pdmfdp(jl,jk)+plude(jl,jk)))
              ptte(jl,jk)=ptte(jl,jk)+zdtdt
              zcpten(jl,jk)=zdtdt*pcpen(jl,jk)
              ptte_cnv(jl,jk)=zdtdt
              zdqdt=-(grav/(paphp1(jl,jk+1)-paphp1(jl,jk)))*              &
                        (pmfuq(jl,jk)+pmfdq(jl,jk)+plude(jl,jk)+       &
                        (pmful(jl,jk)+pdmfup(jl,jk)+pdmfdp(jl,jk)))
              pqte(jl,jk)=pqte(jl,jk)+zdqdt
              zqte_cnv(jl,jk)=zdqdt
              pxtec(jl,jk)=(grav/(paphp1(jl,jk+1)-paphp1(jl,jk)))         &
                           *plude(jl,jk)
              pxteci(jl,jk)=MERGE(0.0_wp,pxtec(jl,jk),llo1)
              zteci(jl,jk)=pxteci(jl,jk)
              pxtecl(jl,jk)=MERGE(pxtec(jl,jk),0.0_wp,llo1)
              ztecl(jl,jk)=pxtecl(jl,jk)
              pqtec(jl,jk)=(grav/(paphp1(jl,jk+1)-paphp1(jl,jk)))         &
                           *pqude(jl,jk)
              zsheat(jl)=zsheat(jl)+zalv*(pdmfup(jl,jk)+pdmfdp(jl,jk))
              zmelt(jl)=zmelt(jl)+pdpmel(jl,jk)
           END IF
230     END DO
!
!        IF (trlist% anyconv /= 0) THEN
!           DO 2304 jt=1,ktrac
!              IF (trlist% ti(jt)% nconv == 1) THEN
!                DO 2302 jl=1,kproma
!                   IF(ldcum(jl)) THEN
!                      zdxtdt=-(grav/(paphp1(jl,jk+1)-paphp1(jl,jk)))      &
!                             *(pmfuxt(jl,jk,jt)+pmfdxt(jl,jk,jt))
!                      pxtte(jl,jk,jt)=pxtte(jl,jk,jt)+zdxtdt
!                   ENDIF
!2302            END DO
!              END IF
!2304       END DO
!        ENDIF
!
     END IF
!
250 END DO
!
!
!---------------------------------------------------------------------
!
!      3.          UPDATE SURFACE FIELDS
!                  ---------------------
!
  DO 310 jl=1,kproma
     prsfc(jl)=prfl(jl)
     pssfc(jl)=psfl(jl)
310 END DO

! diagnostic computations
  DO jk=ktopm2,klev
    DO jl=1,kproma
      pch_concloud(jl)=pch_concloud(jl)+ &
      &   zcpten(jl,jk)*(paphp1(jl,jk+1)-paphp1(jl,jk))/grav
    END DO
  END DO
  DO jl=1,kproma
     pch_concloud(jl)=pch_concloud(jl)-(alv*prsfc(jl)+als*pssfc(jl))
  END DO
! do we need to account for the surface, or for the top 2 layers? 
  DO jk=ktopm2,klev
    DO jl=1,kproma 
      ! water vapor 
      pcon_iqte(jl)=pcon_iqte(jl)+   &
      &   zqte_cnv(jl,jk)*(paphp1(jl,jk+1)-paphp1(jl,jk))/grav
      ! detrained liquid water
      pcon_dtrl(jl)=pcon_dtrl(jl)+ &
      &   ztecl(jl,jk)*(paphp1(jl,jk+1)-paphp1(jl,jk))/grav
      ! detrained ice
      pcon_dtri(jl)=pcon_dtri(jl)+ &
      &   zteci(jl,jk)*(paphp1(jl,jk+1)-paphp1(jl,jk))/grav
    END DO
  END DO
!
! calculate and store convective accumulated precipitation (mm)

!  IF (lanysubmodel) THEN
!    CALL set_vphysc_var &
!         (kproma,        klev,           krow,      &
!         prflconv=prfl, psflconv=psfl)
!  ENDIF

  RETURN
END SUBROUTINE cudtdq

SUBROUTINE cududv(   kproma,   kbdim,    klev,     klevp1,             &
           ktopm2,   ktype,    kcbot,    paphp1,   ldcum,              &
           puen,     pven,     pvom,     pvol,                         &
           puu,      pud,      pvu,      pvd,                          &
           pmfu,     pmfd,     pvom_cnv, pvol_cnv                      )
!
!
!**** *CUDUDV* - UPDATES U AND V TENDENCIES,
!                DOES GLOBAL DIAGNOSTIC OF DISSIPATION
!
!          M.TIEDTKE         E.C.M.W.F.     7/86 MODIF. 12/89
!
!**   INTERFACE.
!     ----------
!
!          *CUDUDV* IS CALLED FROM *CUMASTR*
!
!
INTEGER, INTENT (IN)   :: kproma, kbdim, klev, klevp1, ktopm2
REAL(wp),INTENT(INOUT) :: pvom_cnv(kbdim,klev), pvol_cnv(kbdim,klev)     
!
REAL(wp):: puen(kbdim,klev),        pven(kbdim,klev),                  &
           pvol(kbdim,klev),        pvom(kbdim,klev),                  &
           paphp1(kbdim,klevp1)
REAL(wp):: puu(kbdim,klev),         pud(kbdim,klev),                   &
           pvu(kbdim,klev),         pvd(kbdim,klev),                   &
           pmfu(kbdim,klev),        pmfd(kbdim,klev)
INTEGER :: ktype(kbdim),            kcbot(kbdim)
LOGICAL :: ldcum(kbdim)
!
REAL(wp):: zmfuu(kbdim,klev),       zmfdu(kbdim,klev),                 &
           zmfuv(kbdim,klev),       zmfdv(kbdim,klev)
!
INTEGER :: jl, jk, ik, ikb
REAL(wp):: zzp, zdudt, zdvdt
!
  pvom_cnv(1:kproma,:) = 0._wp
  pvol_cnv(1:kproma,:) = 0._wp
!
!----------------------------------------------------------------------
!
!*    1.0          CALCULATE FLUXES AND UPDATE U AND V TENDENCIES
!                  ----------------------------------------------
!
  IF(ktopm2.EQ.1) THEN
    DO 107 jk=2,klev
       ik=jk-1
       DO 106 jl=1,kproma
          IF(ldcum(jl)) THEN
             zmfuu(jl,jk)=pmfu(jl,jk)*(puu(jl,jk)-puen(jl,ik))
             zmfuv(jl,jk)=pmfu(jl,jk)*(pvu(jl,jk)-pven(jl,ik))
             zmfdu(jl,jk)=pmfd(jl,jk)*(pud(jl,jk)-puen(jl,ik))
             zmfdv(jl,jk)=pmfd(jl,jk)*(pvd(jl,jk)-pven(jl,ik))
          END IF
106    END DO
107  END DO
    DO 105 jl=1,kproma
      IF(ldcum(jl)) THEN
        zmfuu(jl,1)=zmfuu(jl,2)
        zmfuv(jl,1)=zmfuv(jl,2)
        zmfdu(jl,1)=zmfdu(jl,2)
        zmfdv(jl,1)=zmfdv(jl,2)
      END IF
105 END DO
  ELSE
    DO 120 jk=ktopm2,klev
       ik=jk-1
       DO 110 jl=1,kproma
          IF(ldcum(jl)) THEN
             zmfuu(jl,jk)=pmfu(jl,jk)*(puu(jl,jk)-puen(jl,ik))
             zmfuv(jl,jk)=pmfu(jl,jk)*(pvu(jl,jk)-pven(jl,ik))
             zmfdu(jl,jk)=pmfd(jl,jk)*(pud(jl,jk)-puen(jl,ik))
             zmfdv(jl,jk)=pmfd(jl,jk)*(pvd(jl,jk)-pven(jl,ik))
          END IF
110    END DO
120  END DO
  END IF
  DO 140 jk=ktopm2,klev
!DIR$ IVDEP
!OCL NOVREC
     DO 130 jl=1,kproma
        IF(ldcum(jl).AND.jk.GT.kcbot(jl)) THEN
           ikb=kcbot(jl)
           zzp=((paphp1(jl,klevp1)-paphp1(jl,jk))/                     &
                         (paphp1(jl,klevp1)-paphp1(jl,ikb)))
           zzp=MERGE(zzp**2,zzp,ktype(jl).EQ.3)
           zmfuu(jl,jk)=zmfuu(jl,ikb)*zzp
           zmfuv(jl,jk)=zmfuv(jl,ikb)*zzp
           zmfdu(jl,jk)=zmfdu(jl,ikb)*zzp
           zmfdv(jl,jk)=zmfdv(jl,ikb)*zzp
        END IF
130  END DO
140 END DO
!
  DO 190 jk=ktopm2,klev
!
     IF(jk.LT.klev) THEN
        DO 160 jl=1,kproma
           IF(ldcum(jl)) THEN
              zdudt=(grav/(paphp1(jl,jk+1)-paphp1(jl,jk)))*               &
                          (zmfuu(jl,jk+1)-zmfuu(jl,jk)+                &
                           zmfdu(jl,jk+1)-zmfdu(jl,jk))
              zdvdt=(grav/(paphp1(jl,jk+1)-paphp1(jl,jk)))*               &
                          (zmfuv(jl,jk+1)-zmfuv(jl,jk)+                &
                           zmfdv(jl,jk+1)-zmfdv(jl,jk))
              pvom(jl,jk)=pvom(jl,jk)+zdudt
              pvol(jl,jk)=pvol(jl,jk)+zdvdt
              pvom_cnv(jl,jk)=zdudt
              pvol_cnv(jl,jk)=zdvdt
           END IF
160     END DO
!
     ELSE
        DO 170 jl=1,kproma
           IF(ldcum(jl)) THEN
              zdudt=-(grav/(paphp1(jl,jk+1)-paphp1(jl,jk)))*              &
                           (zmfuu(jl,jk)+zmfdu(jl,jk))
              zdvdt=-(grav/(paphp1(jl,jk+1)-paphp1(jl,jk)))*              &
                           (zmfuv(jl,jk)+zmfdv(jl,jk))
              pvom(jl,jk)=pvom(jl,jk)+zdudt
              pvol(jl,jk)=pvol(jl,jk)+zdvdt
              pvom_cnv(jl,jk)=zdudt
              pvol_cnv(jl,jk)=zdvdt
           END IF
170     END DO
     END IF
!
190 END DO
!
!
  RETURN
END SUBROUTINE cududv

END MODULE mo_cufluxdts
