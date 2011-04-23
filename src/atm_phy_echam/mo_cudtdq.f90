!>
!! Module contains subroutine cudtdq.
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
MODULE mo_cudtdq

  USE mo_kind,               ONLY: dp

#ifdef __ICON__
  USE mo_physical_constants, ONLY: alv, als, alf, tmelt, g=>grav
  USE mo_echam_conv_params,  ONLY: ncvmicro
#else
  USE mo_constants,          ONLY: alv, als, alf, tmelt, g
  USE mo_tracdef,            ONLY: trlist
  USE mo_param_switches,     ONLY: ncvmicro
  USE mo_submodel,           ONLY: lanysubmodel
  USE mo_vphysc,             ONLY: set_vphysc_var
#endif

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: cudtdq

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

CONTAINS
  !>
  !!
SUBROUTINE cudtdq(pdelta_time,                                         &
                  kproma, kbdim, klev, klevp1, ktopm2, ldcum, ktrac,   &
!!$                  krow,                                                &
                  paphp1,   pten,     ptte,     pqte,                  &
!!$                  pxtte,                                               &
                  pxtec,                                               &
!!$                  pmfuxt,   pmfdxt,                                    &
                  pmfus,    pmfds,    pmfuq,    pmfdq,                 &
                  pmful,    pdmfup,   pdmfdp,   plude,                 &
                  pdpmel,   prfl,     psfl,                            &
                  pcpen,    pqtec,    pqude,                           &
                  prsfc,    pssfc,    paprc,    paprs,                 &
                  pludel,pludei,pxtecl,pxteci,pmfull,pmfuli,           &!for CDNC/IC
!!$                  plui,                                                &!for CDNC/IC
                  ptte_cnv, pqte_cnv, pxtte_cnv                    )
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
REAL(dp),INTENT(IN) :: pdelta_time
INTEGER, INTENT(IN) :: kproma, kbdim, klev, klevp1, ktopm2, ktrac
!!$INTEGER, INTENT(IN) :: krow
LOGICAL  llo1
!
REAL(dp) :: ptte(kbdim,klev),        pqte(kbdim,klev),                 &
            pten(kbdim,klev),        paphp1(kbdim,klevp1),             &
            paprc(kbdim),            paprs(kbdim),                     &
            prsfc(kbdim),            pssfc(kbdim)
REAL(dp) :: pmfus(kbdim,klev),       pmfds(kbdim,klev),                &
            pmfuq(kbdim,klev),       pmfdq(kbdim,klev),                &
            pmful(kbdim,klev),       plude(kbdim,klev),                &
            pdmfup(kbdim,klev),      pdmfdp(kbdim,klev),               &
            pqtec(kbdim,klev),       pqude(kbdim,klev),                &
            pxtec(kbdim,klev),       prfl(kbdim)
REAL(dp) :: pdpmel(kbdim,klev),      psfl(kbdim)
REAL(dp) :: pcpen(kbdim,klev)
LOGICAL  :: ldcum(kbdim)
!
REAL(dp) :: zmelt(kbdim)
REAL(dp) :: zsheat(kbdim)
!!$REAL(dp) :: pxtte(kbdim,klev,ktrac)
!!$REAL(dp) :: pmfuxt(kbdim,klev,ktrac), pmfdxt(kbdim,klev,ktrac)
!
REAL(dp) :: zrcpm ! reciprocal value of specific heat of moist air

REAL(dp),INTENT(OUT) ::  ptte_cnv(kbdim,klev), pqte_cnv(kbdim,klev)
REAL(dp),INTENT(OUT) :: pxtte_cnv(kbdim,klev,ktrac)

!----Included for CDNC/IC scheme (Ulrike Lohmann 11/02/07)-----------
REAL(dp) ::  pludel(kbdim,klev),       pludei(kbdim,klev),              &
             pmfull(kbdim,klev),       pmfuli(kbdim,klev),              &
             pxtecl(kbdim,klev),       pxteci(kbdim,klev)
!!$REAL(dp) ::  plui(kbdim,klev)
!----End Included for CDNC/IC scheme (Ulrike Lohmann 11/02/07)-----------

INTEGER  :: jl, jk
!!$INTEGER  :: jt
REAL(dp) :: zdiagt, zalv, zdtdt, zdqdt
!!$REAL(dp) :: zdxtdt

   ptte_cnv(1:kproma,:)   = 0._dp
   pqte_cnv(1:kproma,:)   = 0._dp
  pxtte_cnv(1:kproma,:,:) = 0._dp
!
!----------------------------------------------------------------------
!
!*    1.0          SPECIFY PARAMETERS
!                  ------------------
!
!100 CONTINUE
  zdiagt=pdelta_time
!
!----------------------------------------------------------------------
!
!*    2.0          INCREMENTATION OF T AND Q TENDENCIES
!                  ------------------------------------
!
!200 CONTINUE
  DO 210 jl=1,kproma
     zmelt(jl)=0._dp
     zsheat(jl)=0._dp
210 END DO
!
  DO 250 jk=ktopm2,klev
!
     IF(jk.LT.klev) THEN
        DO 220 jl=1,kproma
           IF(ldcum(jl)) THEN
!----Included for CDNC/IC scheme (Ulrike Lohmann 11/02/07)---------
            IF (ncvmicro>0) THEN
!              llo1=(pten(jl,jk)-tmelt).GT.0._dp.OR.(pten(jl,jk)>cthomi.AND.plui(jl,jk)<csecfrl)
              llo1=(pten(jl,jk)-tmelt).GT.0._dp
              zalv=MERGE(alv,als,llo1)
              zrcpm=1._dp/pcpen(jl,jk)
              zdtdt=(g/(paphp1(jl,jk+1)-paphp1(jl,jk)))*zrcpm*        &
                                  (pmfus(jl,jk+1)-pmfus(jl,jk)+       &
                                   pmfds(jl,jk+1)-pmfds(jl,jk)-       &
                                   alf*pdpmel(jl,jk)-                 &
                                   alv*(pmfull(jl,jk+1)-pmfull(jl,jk)- &
                                   pludel(jl,jk))-                      &
                                   als*(pmfuli(jl,jk+1)-pmfuli(jl,jk)- &
                                   pludei(jl,jk))                      &
                                   +zalv*(pdmfup(jl,jk)+pdmfdp(jl,jk)))
              ptte(jl,jk)=ptte(jl,jk)+zdtdt
              ptte_cnv(jl,jk)=zdtdt
              zdqdt=(g/(paphp1(jl,jk+1)-paphp1(jl,jk)))*              &
                                  (pmfuq(jl,jk+1)-pmfuq(jl,jk)+       &
                                   pmfdq(jl,jk+1)-pmfdq(jl,jk)+       &
                                   pmfull(jl,jk+1)-pmfull(jl,jk)-       &
                                   pludel(jl,jk)+                      &
                                   pmfuli(jl,jk+1)-pmfuli(jl,jk)-       &
                                   pludei(jl,jk)                      &
                                  -(pdmfup(jl,jk)+pdmfdp(jl,jk)))
              pqte(jl,jk)=pqte(jl,jk)+zdqdt
              pqte_cnv(jl,jk)=zdqdt
              pxtecl(jl,jk)=(g/(paphp1(jl,jk+1)-                       &
                            paphp1(jl,jk)))*pludel(jl,jk)
              pxteci(jl,jk)=(g/(paphp1(jl,jk+1)-                       &
                            paphp1(jl,jk)))*pludei(jl,jk)
              pqtec(jl,jk)=(g/(paphp1(jl,jk+1)-                       &
                            paphp1(jl,jk)))*pqude(jl,jk)
              zsheat(jl)=zsheat(jl)+zalv*(pdmfup(jl,jk)+pdmfdp(jl,jk))
              zmelt(jl)=zmelt(jl)+pdpmel(jl,jk)
           ELSE
              llo1=(pten(jl,jk)-tmelt).GT.0._dp
              zalv=MERGE(alv,als,llo1)
              zrcpm=1._dp/pcpen(jl,jk)
              zdtdt=(g/(paphp1(jl,jk+1)-paphp1(jl,jk)))*zrcpm*         &
                                  (pmfus(jl,jk+1)-pmfus(jl,jk)+        &
                                   pmfds(jl,jk+1)-pmfds(jl,jk)-        &
                                   alf*pdpmel(jl,jk)-                  &
                                   zalv*(pmful(jl,jk+1)-pmful(jl,jk)-  &
                                   plude(jl,jk)-                       &
                                  (pdmfup(jl,jk)+pdmfdp(jl,jk))))
              ptte(jl,jk)=ptte(jl,jk)+zdtdt
              ptte_cnv(jl,jk)=zdtdt
              zdqdt=(g/(paphp1(jl,jk+1)-paphp1(jl,jk)))*               &
                                  (pmfuq(jl,jk+1)-pmfuq(jl,jk)+        &
                                   pmfdq(jl,jk+1)-pmfdq(jl,jk)+        &
                                   pmful(jl,jk+1)-pmful(jl,jk)-        &
                                   plude(jl,jk)-                       &
                                  (pdmfup(jl,jk)+pdmfdp(jl,jk)))
              pqte(jl,jk)=pqte(jl,jk)+zdqdt
              pqte_cnv(jl,jk)=zdqdt
              pxtec(jl,jk)=(g/(paphp1(jl,jk+1)-                        &
                            paphp1(jl,jk)))*plude(jl,jk)
              pqtec(jl,jk)=(g/(paphp1(jl,jk+1)-                        &
                            paphp1(jl,jk)))*pqude(jl,jk)
              zsheat(jl)=zsheat(jl)+zalv*(pdmfup(jl,jk)+pdmfdp(jl,jk))
              zmelt(jl)=zmelt(jl)+pdpmel(jl,jk)
            ENDIF
!---End Included for CDNC/IC scheme (Ulrike Lohmann 11/02/07)---------

           END IF
220     END DO
!
!!$#ifdef __ICON__
!!$#else
!!$        IF (trlist% anyconv /= 0) THEN
!!$           DO 2204 jt=1,ktrac
!!$              IF (trlist% ti(jt)% nconv == 1) THEN
!!$                DO 2202 jl=1,kproma
!!$                   IF(ldcum(jl)) THEN
!!$                     zdxtdt=(g/(paphp1(jl,jk+1)-paphp1(jl,jk)))        &
!!$                                 *(pmfuxt(jl,jk+1,jt)-pmfuxt(jl,jk,jt) &
!!$                                  +pmfdxt(jl,jk+1,jt)-pmfdxt(jl,jk,jt))
!!$                     pxtte(jl,jk,jt)=pxtte(jl,jk,jt)+zdxtdt
!!$                     pxtte_cnv(jl,jk,jt)=zdxtdt
!!$                   ENDIF
!!$2202            END DO
!!$              ENDIF
!!$2204       END DO
!!$        ENDIF
!!$#endif
!
!
     ELSE
        DO 230 jl=1,kproma
           IF(ldcum(jl)) THEN
!----Included for CDNC/IC scheme (Ulrike Lohmann 11/02/07)---------
            IF (ncvmicro>0) THEN
!              llo1=(pten(jl,jk)-tmelt).GT.0._dp.OR.(pten(jl,jk)>cthomi.AND.plui(jl,jk)<csecfrl)
              llo1=(pten(jl,jk)-tmelt).GT.0._dp
              zalv=MERGE(alv,als,llo1)
              zrcpm=1._dp/pcpen(jl,jk)
              zdtdt=-(g/(paphp1(jl,jk+1)-paphp1(jl,jk)))*zrcpm*       &
                     (pmfus(jl,jk)+pmfds(jl,jk)+alf*pdpmel(jl,jk)-    &
                                 alv*(pmfull(jl,jk)+pludel(jl,jk))-   &
                                 als*(pmfuli(jl,jk)+pludei(jl,jk))    &
                                -zalv*(pdmfdp(jl,jk)+pdmfup(jl,jk)))
              ptte(jl,jk)=ptte(jl,jk)+zdtdt
              ptte_cnv(jl,jk)=zdtdt
              zdqdt=-(g/(paphp1(jl,jk+1)-paphp1(jl,jk)))*             &
                        (pmfuq(jl,jk)+pmfdq(jl,jk)+pludel(jl,jk)+pludei(jl,jk)+     &
                        (pmfull(jl,jk)+pmfuli(jl,jk)+pdmfup(jl,jk)+pdmfdp(jl,jk)))
              pqte(jl,jk)=pqte(jl,jk)+zdqdt
              pqte_cnv(jl,jk)=zdqdt
              pxtecl(jl,jk)=(g/(paphp1(jl,jk+1)-paphp1(jl,jk)))        &
                           *pludel(jl,jk)
              pxteci(jl,jk)=(g/(paphp1(jl,jk+1)-paphp1(jl,jk)))        &
                           *pludei(jl,jk)
              pqtec(jl,jk)=(g/(paphp1(jl,jk+1)-paphp1(jl,jk)))        &
                           *pqude(jl,jk)
              zsheat(jl)=zsheat(jl)+zalv*(pdmfup(jl,jk)+pdmfdp(jl,jk))
              zmelt(jl)=zmelt(jl)+pdpmel(jl,jk)
            ELSE
              llo1=(pten(jl,jk)-tmelt).GT.0._dp
              zalv=MERGE(alv,als,llo1)
              zrcpm=1._dp/pcpen(jl,jk)
              zdtdt=-(g/(paphp1(jl,jk+1)-paphp1(jl,jk)))*zrcpm*        &
                     (pmfus(jl,jk)+pmfds(jl,jk)+alf*pdpmel(jl,jk)-     &
                                 zalv*(pmful(jl,jk)+pdmfup(jl,jk)      &
                                +pdmfdp(jl,jk)+plude(jl,jk)))
              ptte(jl,jk)=ptte(jl,jk)+zdtdt
              ptte_cnv(jl,jk)=zdtdt
              zdqdt=-(g/(paphp1(jl,jk+1)-paphp1(jl,jk)))*              &
                        (pmfuq(jl,jk)+pmfdq(jl,jk)+plude(jl,jk)+       &
                        (pmful(jl,jk)+pdmfup(jl,jk)+pdmfdp(jl,jk)))
              pqte(jl,jk)=pqte(jl,jk)+zdqdt
              pqte_cnv(jl,jk)=zdqdt
              pxtec(jl,jk)=(g/(paphp1(jl,jk+1)-paphp1(jl,jk)))         &
                           *plude(jl,jk)
              pqtec(jl,jk)=(g/(paphp1(jl,jk+1)-paphp1(jl,jk)))         &
                           *pqude(jl,jk)
              zsheat(jl)=zsheat(jl)+zalv*(pdmfup(jl,jk)+pdmfdp(jl,jk))
              zmelt(jl)=zmelt(jl)+pdpmel(jl,jk)
            END IF
          END IF
230     END DO
!
!!$#ifdef __ICON__
!!$#else
!!$        IF (trlist% anyconv /= 0) THEN
!!$           DO 2304 jt=1,ktrac
!!$              IF (trlist% ti(jt)% nconv == 1) THEN
!!$                DO 2302 jl=1,kproma
!!$                   IF(ldcum(jl)) THEN
!!$                      zdxtdt=-(g/(paphp1(jl,jk+1)-paphp1(jl,jk)))      &
!!$                             *(pmfuxt(jl,jk,jt)+pmfdxt(jl,jk,jt))
!!$                      pxtte(jl,jk,jt)=pxtte(jl,jk,jt)+zdxtdt
!!$                      pxtte_cnv(jl,jk,jt)=zdxtdt
!!$                   ENDIF
!!$2302            END DO
!!$              END IF
!!$2304       END DO
!!$        ENDIF
!!$#endif
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
!300 CONTINUE
  DO 310 jl=1,kproma
     prsfc(jl)=prfl(jl)
     pssfc(jl)=psfl(jl)
     paprc(jl)=paprc(jl)+zdiagt*(prfl(jl)+psfl(jl))
     paprs(jl)=paprs(jl)+zdiagt*psfl(jl)
310 END DO
!
! calculate and store convective accumulated precipitation (mm)

!!$#ifdef __ICON__
!!$#else
!!$  IF (lanysubmodel) CALL set_vphysc_var &
!!$                    (kproma,        klev,           krow,      &
!!$!++mgs
!!$                     prflconv=prfl, psflconv=psfl)
!!$!--mgs
!!$#endif

    RETURN
  END SUBROUTINE cudtdq

END MODULE mo_cudtdq
