!>
!! Module contains subroutine cuflx.
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
MODULE mo_cuflx

  USE mo_kind,                ONLY: dp

#ifdef __ICON__
  USE mo_physical_constants,  ONLY: g=>grav, alf, cpd, tmelt, vtmpc2
  USE mo_echam_conv_nml,      ONLY: cevapcu, ncvmicro
#else
  USE mo_constants,           ONLY: g, alf, cpd, tmelt, vtmpc2, rd
  USE mo_physc2,              ONLY: cevapcu
  USE mo_submodel,            ONLY: lanysubmodel, lham
  USE mo_submodel_interface,  ONLY: cuflx_subm
  USE mo_param_switches,      ONLY: ncvmicro
#endif

  IMPLICIT NONE
  PRIVATE
  PUBLIC  :: cuflx

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

CONTAINS
  !>
  !!
!++mgs: zcucov and zdpevap now 1d vectors
SUBROUTINE cuflx(  ptime_step_len,  kproma, kbdim, klev, klevp1,       &
           pqen,     pqsen,    ptenh,    pqenh,                        &
           ktrac,                                                      &
!---Included for scavenging in xtwetdep (Philip Stier, 28/03/01, UL, 2803.07):-------
!!$           krow,                                                       &
!!$           pxtte,    pxtu,     ptu,                                    &
!!$           pmwc,     pmrateprecip,  pmratesnow,                        &
!---End Included for scavenging-----------------------------------------
           pxtenh,   pmfuxt,   pmfdxt,                                 &
           paphp1,   pgeoh,                                            &
           kcbot,    kctop,    kdtop,                                  &
           ktype,    lddraf,   ldcum,                                  &
           pmfu,     pmfd,     pmfus,    pmfds,                        &
           pmfuq,    pmfdq,    pmful,                                  &
           pdmfup,   pdmfdp,   prfl,     prain,                        &
           pcpcu,                                                      &
           pten,     psfl,     pdpmel,   ktopm2,                       &
!-----------------------added by Junhua Zhang for Micro---------------
           pmfull,pmfuli)
!!$           plul,     plui)
!-------------------------------------end-----------------------------
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
REAL(dp),INTENT(IN) :: ptime_step_len

INTEGER, INTENT (IN) :: kproma, kbdim, klev, klevp1, ktrac
INTEGER, INTENT (OUT):: ktopm2
!
REAL(dp):: pqen(kbdim,klev),        pqsen(kbdim,klev),                 &
           ptenh(kbdim,klev),       pqenh(kbdim,klev),                 &
           paphp1(kbdim,klevp1),    pgeoh(kbdim,klev)
!
REAL(dp):: pmfu(kbdim,klev),        pmfd(kbdim,klev),                  &
           pmfus(kbdim,klev),       pmfds(kbdim,klev),                 &
           pmfuq(kbdim,klev),       pmfdq(kbdim,klev),                 &
           pdmfup(kbdim,klev),      pdmfdp(kbdim,klev),                &
           pmful(kbdim,klev),                                          &
           prfl(kbdim),             prain(kbdim)
REAL(dp):: pten(kbdim,klev),        pdpmel(kbdim,klev),                &
           psfl(kbdim)
REAL(dp):: pcpcu(kbdim,klev)
INTEGER :: kcbot(kbdim),            kctop(kbdim),                      &
           kdtop(kbdim),            ktype(kbdim)
LOGICAL :: lddraf(kbdim),           ldcum(kbdim)
REAL(dp):: pxtenh(kbdim,klev,ktrac),                                   &
           pmfuxt(kbdim,klev,ktrac),pmfdxt(kbdim,klev,ktrac)
!
!----------added by Junhua Zhang, Ulrike Lohmann for Micro---------------

!!$REAL(dp):: plul(kbdim,klev),        plui(kbdim,klev)
REAL(dp):: pmfull(kbdim,klev),      pmfuli(kbdim,klev)
!-------------------------------------end-----------------------------
!
INTEGER :: jl, jk, jt, ikb
!++mgs
!!REAL(dp):: zcons1, zcons2, zcucov, ztmelp2, zzp, zfac, zsnmlt          &
!!         , zrfl, zrnew, zrmin, zrfln, zdrfl, zrsum, zdpevap,zwu
!!REAL(dp):: zpsubcl(kbdim)
REAL(dp):: zcons1, zcons2,         ztmelp2, zzp, zfac, zsnmlt          &
         , zrfl, zrnew, zrmin, zrfln, zdrfl, zrsum
!!$REAL(dp):: zwu
REAL(dp):: zpsubcl(kbdim)
!!$REAL(dp):: zpsubcl_sav(kbdim)
REAL(dp):: zcucov(kbdim), zdpevap(kbdim)
!--mgs
!---Included for scavenging in xtwetdep (Philip Stier, 28/03/01):-------
!
!!$INTEGER, INTENT(IN) :: krow
!!$
!!$REAL(dp):: pxtte(kbdim,klev,ktrac), pxtu(kbdim,klev,ktrac),            &
!!$           ptu(kbdim,klev)
!!$
!!$REAL(dp):: zmlwc(kbdim,klev),       zmiwc(kbdim,klev),                 &
!!$           zmratepr(kbdim,klev),    zmrateps(kbdim,klev),              &
!!$           zfrain(kbdim,klev),      zfsnow(kbdim,klev),                &
!!$           zdpg(kbdim,klev),        zfevapr(kbdim,klev),               &
!!$           zfsubls(kbdim,klev),     zaclc(kbdim,klev),                 &
!!$           zmsnowacl(kbdim,klev),   zrhou(kbdim,klev)
!!$
!!$REAL(dp):: pmratesnow(kbdim,klev),            &
!!$           pmwc(kbdim,klev),        pmrateprecip(kbdim,klev)
!!$
!!$REAL(dp):: ztc,         zzfac
!!$
!!$REAL(dp):: zcaa, &   !(Constants for partitioning of cloud water
!!$           zcab      ! into liquid and solid part)
!!$
!!$REAL(dp):: zalpha    ! Fraction of cloud water in liquid phase
!!$                     ! =>  (1-zalpha)      "   in solid  phase

!---End Included for scavenging-----------------------------------------
!
!*             SPECIFY CONSTANTS
!

  zcons1=cpd/(alf*g*ptime_step_len)
  zcons2=1._dp/(g*ptime_step_len)
!!mgs!!  zcucov=0.05_dp            !! ++mgs: now 1-d array
  ztmelp2=tmelt+2._dp
!
!---Included for scavenging in xtwetdep (Philip Stier, 28/03/01):-------
!
!!$#ifdef __ICON__
!!$#else
!!$IF(lanysubmodel) THEN
!!$
!!$  zcaa = 0.0059_dp   ! (Constants for partitioning of cloud water
!!$  zcab = 0.003102_dp !  into liquid and solid part)
!!$
!!$  zmlwc(1:kproma,:)     = 0._dp
!!$  zmiwc(1:kproma,:)     = 0._dp
!!$  zmratepr(1:kproma,:)  = 0._dp
!!$  zmrateps(1:kproma,:)  = 0._dp
!!$  zfrain(1:kproma,:)    = 0._dp
!!$  zfsnow(1:kproma,:)    = 0._dp
!!$  zdpg(1:kproma,:)      = 0._dp
!!$  zfevapr(1:kproma,:)   = 0._dp
!!$  zfsubls(1:kproma,:)   = 0._dp
!!$  zmsnowacl(1:kproma,:) = 0._dp
!!$
!!$  zwu            = 2.0_dp
!!$
!!$  !--- Set cloud cover to 1 below cloud top and bottom:
!!$
!!$  DO jk=1, klev
!!$     DO jl=1, kproma
!!$        IF (jk>=kctop(jl) .AND. jk<=kcbot(jl)) THEN
!!$           zaclc(jl,jk) = 1._dp
!!$        ELSE
!!$           zaclc(jl,jk) = 0._dp
!!$        END IF
!!$     END DO
!!$  END DO
!!$END IF
!!$#endif
!
!---End Included for scavenging-----------------------------------------
!
!
!
!*    1.0          DETERMINE FINAL CONVECTIVE FLUXES
!                  ---------------------------------
!
!100 CONTINUE
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
              pmfd(jl,jk)=0._dp
              pmfds(jl,jk)=0._dp
              pmfdq(jl,jk)=0._dp
              pdmfdp(jl,jk-1)=0._dp
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
                 pmfdxt(jl,jk,jt)=0._dp
              ENDIF
           ELSE
              pmfuxt(jl,jk,jt)=0._dp
              pmfdxt(jl,jk,jt)=0._dp
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
!           pmful(jl,jk)=pmful(jl,ikb)*zzp
!-------------Added by Junhua Zhang for CONV Micro--------------
           IF (ncvmicro>0) THEN
              pmfull(jl,jk)=pmfull(jl,ikb)*zzp
              pmfuli(jl,jk)=pmfuli(jl,ikb)*zzp
           ELSE
              pmful(jl,jk)=pmful(jl,ikb)*zzp
           ENDIF
!----------------------------End--------------------------------
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
!200 CONTINUE
  DO 210 jl=1,kproma
     prfl(jl)=0._dp
     psfl(jl)=0._dp
     prain(jl)=0._dp
210 END DO

!!$#ifdef __ICON__
!!$#else
!!$!++jsr interface for scavenging
!!$IF (lanysubmodel) THEN
!!$   DO jk=ktopm2,klev
!!$!---Included for scavenging in xtwetdep (Philip Stier, 28/03/01):-------
!!$      zdpg(1:kproma,jk)=(paphp1(1:kproma,jk+1)-paphp1(1:kproma,jk))/g
!!$!---End Included for scavenging-----------------------------------------
!!$   END DO
!!$   IF (ncvmicro>0) THEN
!!$      DO jk=ktopm2,klev
!!$         WHERE (ldcum(1:kproma))
!!$!---Included for scavenging in xtwetdep (Philip Stier, 28/03/01):-------
!!$            zmlwc(1:kproma,jk)=plul(1:kproma,jk)
!!$            zmiwc(1:kproma,jk)=plui(1:kproma,jk)
!!$
!!$            zmratepr(1:kproma,jk)=pmrateprecip(1:kproma,jk)
!!$            zmrateps(1:kproma,jk)=pmratesnow(1:kproma,jk)
!!$         END WHERE
!!$      END DO
!!$   ELSE
!!$      DO jk=ktopm2,klev
!!$         DO jl=1,kproma
!!$            IF (ldcum(jl)) THEN
!!$               ztc=ptu(jl,jk)-tmelt
!!$               zzfac=MERGE(1._dp,0._dp,ztc<0._dp)
!!$               zalpha=(1._dp-zzfac)+zzfac*(zcaa+(1._dp-zcaa)*exp(-zcab*ztc**2))
!!$
!!$               zmlwc(jl,jk)=zalpha*pmwc(jl,jk)
!!$               zmiwc(jl,jk)=(1._dp-zalpha)*pmwc(jl,jk)
!!$
!!$               zmratepr(jl,jk)=zalpha*pmrateprecip(jl,jk)
!!$               zmrateps(jl,jk)=(1._dp-zalpha)*pmrateprecip(jl,jk)
!!$            END IF
!!$         END DO
!!$      END DO
!!$   END IF
!!$END IF
!!$!---End Included for scavenging-----------------------------------------
!!$#endif
!--jsr interface for scavenging

  DO 220 jk=ktopm2,klev
     DO 215 jl=1,kproma
        IF(ldcum(jl)) THEN
           prain(jl)=prain(jl)+pdmfup(jl,jk)
           IF(pten(jl,jk).GT.tmelt) THEN
              prfl(jl)=prfl(jl)+pdmfup(jl,jk)+pdmfdp(jl,jk)
              IF(psfl(jl).GT.0._dp.AND.pten(jl,jk).GT.ztmelp2) THEN
                 zfac=zcons1*(1._dp+vtmpc2*pqen(jl,jk))                &
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

#ifdef __ICON__
#else
     IF (lanysubmodel) THEN
!---Included for scavenging in xtwetdep (Philip Stier, 28/03/01):-------
        zfrain(1:kproma,jk)=prfl(1:kproma)
        zfsnow(1:kproma,jk)=psfl(1:kproma)
!---End Included for scavenging-----------------------------------------
     END IF
#endif

220 END DO
  DO 230 jl=1,kproma
     prfl(jl)=MAX(prfl(jl),0._dp)
     psfl(jl)=MAX(psfl(jl),0._dp)
     zpsubcl(jl)=prfl(jl)+psfl(jl)
230 END DO
!++jsr scavenging interface
!++mgs ### new code
!!IF (lanysubmodel) THEN
!!   IF (lham) THEN
!!      DO jk=ktopm2,klev
!!         DO jl=1,kproma
!!            IF(ldcum(jl).AND.jk.GE.kcbot(jl).AND.zpsubcl(jl).GT.1.e-20_dp) THEN
!!               zrhou(jl,jk)=(paphp1(jl,jk+1)+paphp1(jl,jk))* &
!!                            0.5_dp/(ptu(jl,jk)*rd)
!!               zcucov=pmfu(jl,jk)/(zwu*zrhou(jl,jk))
!!               zrfl=zpsubcl(jl)
!!               zrnew=(MAX(0._dp,SQRT(zrfl/zcucov)-                         &
!!                      cevapcu(jk)*(paphp1(jl,jk+1)-paphp1(jl,jk))*   &
!!                      MAX(0._dp,pqsen(jl,jk)-pqen(jl,jk))))**2*zcucov
!!               zrmin=zrfl-zcucov*MAX(0._dp,0.8_dp*pqsen(jl,jk)-pqen(jl,jk))&
!!                     *zcons2*(paphp1(jl,jk+1)-paphp1(jl,jk))
!!               zrnew=MAX(zrnew,zrmin)
!!               zrfln=MAX(zrnew,0._dp)
!!               zdrfl=MIN(0._dp,zrfln-zrfl)
!!               pdmfup(jl,jk)=pdmfup(jl,jk)+zdrfl
!!!---Included for scavenging in xtwetdep (Philip Stier, 28/03/01):-------
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
!!            IF(ldcum(jl).AND.jk.GE.kcbot(jl).AND.zpsubcl(jl).GT.1.e-20_dp) THEN
!!               zrhou(jl,jk)=(paphp1(jl,jk+1)+paphp1(jl,jk))* &
!!                            0.5_dp/(ptu(jl,jk)*rd)
!!               zrfl=zpsubcl(jl)
!!               zrnew=(MAX(0._dp,SQRT(zrfl/zcucov)-                         &
!!                      cevapcu(jk)*(paphp1(jl,jk+1)-paphp1(jl,jk))*   &
!!                      MAX(0._dp,pqsen(jl,jk)-pqen(jl,jk))))**2*zcucov
!!               zrmin=zrfl-zcucov*MAX(0._dp,0.8_dp*pqsen(jl,jk)-pqen(jl,jk))&
!!                     *zcons2*(paphp1(jl,jk+1)-paphp1(jl,jk))
!!               zrnew=MAX(zrnew,zrmin)
!!               zrfln=MAX(zrnew,0._dp)
!!               zdrfl=MIN(0._dp,zrfln-zrfl)
!!               pdmfup(jl,jk)=pdmfup(jl,jk)+zdrfl
!!!---Included for scavenging in xtwetdep (Philip Stier, 28/03/01):-------
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
!!       IF(ldcum(jl).AND.jk.GE.kcbot(jl).AND.zpsubcl(jl).GT.1.e-20_dp)  &
!!                                                                  THEN
!!           zrfl=zpsubcl(jl)
!!           zrnew=(MAX(0._dp,SQRT(zrfl/zcucov)-                         &
!!                        cevapcu(jk)*(paphp1(jl,jk+1)-paphp1(jl,jk))*   &
!!                        MAX(0._dp,pqsen(jl,jk)-pqen(jl,jk))))**2*zcucov
!!           zrmin=zrfl-zcucov*MAX(0._dp,0.8_dp*pqsen(jl,jk)-pqen(jl,jk))&
!!                        *zcons2*(paphp1(jl,jk+1)-paphp1(jl,jk))
!!           zrnew=MAX(zrnew,zrmin)
!!           zrfln=MAX(zrnew,0._dp)
!!           zdrfl=MIN(0._dp,zrfln-zrfl)
!!           pdmfup(jl,jk)=pdmfup(jl,jk)+zdrfl
!!           zpsubcl(jl)=zrfln
!!       END IF
!!235  END DO
!!240 END DO
!!END IF
#ifdef __ICON__
#else
  IF (lanysubmodel) THEN
    DO jk=ktopm2,klev
      zrhou(1:kproma,jk)=(paphp1(1:kproma,jk+1)+paphp1(1:kproma,jk))* &
                          0.5_dp/(ptu(1:kproma,jk)*rd)
    END DO
  END IF
#endif

  DO jk=ktopm2,klev

    zdpevap(1:kproma) = 0._dp

!!$#ifdef __ICON__
    zcucov(1:kproma) = 0.05_dp

!!$#else
!!$    IF (lanysubmodel) zpsubcl_sav(1:kproma) = zpsubcl(1:kproma)
!!$    ! ### value of zcucov depends on lham => explicit submodel dependence !
!!$    IF (lham) THEN
!!$      zcucov(1:kproma) = pmfu(1:kproma,jk)/(zwu*zrhou(1:kproma,jk))
!!$    ELSE
!!$      zcucov(1:kproma) = 0.05_dp
!!$    END IF
!!$#endif

    DO jl=1,kproma
      IF(ldcum(jl).AND.jk.GE.kcbot(jl).AND.zpsubcl(jl).GT.1.e-20_dp) THEN
        zrfl=zpsubcl(jl)
        zrnew=(MAX(0._dp,SQRT(zrfl/zcucov(jl))-                           &
               cevapcu(jk)*(paphp1(jl,jk+1)-paphp1(jl,jk))*               &
               MAX(0._dp,pqsen(jl,jk)-pqen(jl,jk))))**2*zcucov(jl)
        zrmin=zrfl-zcucov(jl)*MAX(0._dp,0.8_dp*pqsen(jl,jk)-pqen(jl,jk))  &
              *zcons2*(paphp1(jl,jk+1)-paphp1(jl,jk))
        zrnew=MAX(zrnew,zrmin)
        zrfln=MAX(zrnew,0._dp)
        zdrfl=MIN(0._dp,zrfln-zrfl)
        zdpevap(jl) = -zdrfl
        pdmfup(jl,jk)=pdmfup(jl,jk)+zdrfl
        zpsubcl(jl)=zrfln
      END IF
    END DO

#ifdef __ICON__
#else
!---Included for scavenging in xtwetdep (Philip Stier, 28/03/01):-------
    IF (lanysubmodel) THEN
      DO jl=1,kproma
        IF (zpsubcl_sav(jl) > 1.e-20_dp) THEN
          zfevapr(jl,jk)=zdpevap(jl)*prfl(jl)/zpsubcl_sav(jl)
          zfsubls(jl,jk)=zdpevap(jl)*psfl(jl)/zpsubcl_sav(jl)
        ELSE
          zfevapr(jl,jk) = 0._dp
          zfsubls(jl,jk) = 0._dp
        END IF
      END DO
    END IF
!---End Included for scavenging-----------------------------------------
#endif

  END DO !jk=ktopm2,klev

!!$#ifdef __ICON__
!!$#else
!!$!--mgs ### new code
!!$
!!$  IF (lanysubmodel) THEN
!!$    CALL cuflx_subm(kbdim,  kproma,    klev,     ktopm2,   & ! dimensions
!!$                    krow,                                  & ! longitude (kproma-block) index
!!$                    pxtenh, pxtu, zrhou,                   & ! tracers
!!$                    pmfu,   pmfuxt,                        & ! convective fluxes and corresp. mmr
!!$                    zmlwc,  zmiwc,     zmratepr, zmrateps, & ! cloud properties
!!$                    zfrain, zfsnow,    zfevapr,  zfsubls,  & !   "       "
!!$                    zaclc,  zmsnowacl,                     & !   "       "
!!$                    ptu,    zdpg,                          & ! thermodynamic quantities
!!$                    pxtte                                  )
!!$  END IF
!!$#endif

!!baustelle!! (?)
  DO 250 jl=1,kproma
     zrsum=prfl(jl)+psfl(jl)
     zdpevap(jl)=zpsubcl(jl)-zrsum
     prfl(jl)=prfl(jl)+zdpevap(jl)*prfl(jl)*(1._dp/MAX(1.e-20_dp,zrsum))
     psfl(jl)=psfl(jl)+zdpevap(jl)*psfl(jl)*(1._dp/MAX(1.e-20_dp,zrsum))
250 END DO
!
    RETURN
  END SUBROUTINE cuflx

END MODULE mo_cuflx
