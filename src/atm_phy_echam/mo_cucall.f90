!>
!! Module contains the subroutine CUCALL, the top level routine for
!! cumulus convection.
!!
!! *CUCALL* - INTERFACE FOR *CUMASTR*:
!!            PROVIDES INPUT FOR CUMASTR
!!            RECEIVES UPDATED TENDENCIES, PRECIPITATION.
!!
!! *CUCALL* IS CALLED FROM *PHYSC*
!!
!! @author M. Tiedtke, ECMWF
!!
!! @par Revision History
!! - Original version by M. Tiedtke, ECMWF (1989-12)
!! - Taken from ECHAM6, wrapped in module and modified for ICON
!!   by Hui Wan, MPI-M (2010-08)
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
#ifdef __xlC__
@PROCESS HOT
#endif

MODULE mo_cucall

  USE mo_kind,               ONLY: dp
  USE mo_exception,          ONLY: finish

  USE mo_physical_constants, ONLY: vtmpc1    ! vtmpc1=rv/rd-1

  USE mo_convect_tables,     ONLY: prepare_ua_index,lookup_ua
  USE mo_cumastr,            ONLY: cumastr

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: cucall

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

CONTAINS
  !>
  !!
  SUBROUTINE cucall( iconv, lmfdudv, lmfdd, lmfmid, dlev,            &! in
                     cmftau, cmfctop, cprcon, cminbuoy, entrpen,     &! in
                     nmctop, cevapcu,                                &! in
                     kproma, kbdim, klev, klevp1, klevm1,            &! in
                     ktrac,                                          &! in
!0                   krow,                                           &! in
                     pdtime, ptime_step_len,                         &! in
                     ldland,                                         &! in
                     ptm1,     pum1,     pvm1,                       &! in
                     pverv,                                          &! in
                     pgeo,                                           &! in
                     pqm1,     pxlm1,    pxim1, pxtm1,               &! in
                     pxlte,    pxite,                                &! in
                     papp1,    paphp1,                               &! in
                     pqhfla,                                         &! in
!0                   ptkem1,                                         &! in
                     pthvsig,                                        &! in
                     ptte,     pvom,     pvol,                       &! inout
                     pqte,     pxtte,                                &! inout
                     pqtec,    pxtec,                                &! inout
                     prsfc,    pssfc,                                &! out
                     ptopmax,                                        &! inout
                     ktype,                                          &! inout
                     ilab,                                           &! out
                     pcvcbot,  pwcape,                               &! out
                     ptte_cnv, pvom_cnv, pvol_cnv, pqte_cnv,         &! out
                     pxtte_cnv                               )        ! out

    INTEGER, INTENT(IN) :: iconv
    LOGICAL, INTENT(IN) :: lmfdudv, lmfdd, lmfmid 
    REAL(dp),INTENT(IN) :: dlev, cmftau, cmfctop, cprcon, cminbuoy, entrpen
    INTEGER, INTENT(IN) :: nmctop 
    INTEGER, INTENT(IN) :: klev, klevm1, klevp1, kproma, kbdim, ktrac
    REAL(dp),INTENT(IN) :: cevapcu(:) 
!0  INTEGER, INTENT(IN) :: krow
    REAL(dp),INTENT(IN) :: pdtime
    REAL(dp),INTENT(IN) :: ptime_step_len
    LOGICAL, INTENT(IN) :: ldland(kbdim)
    REAL(dp),INTENT(IN) :: ptm1(kbdim,klev)
    REAL(dp),INTENT(IN) :: pum1(kbdim,klev),      pvm1(kbdim,klev)
    REAL(dp),INTENT(IN) :: pverv(kbdim,klev),     pgeo(kbdim,klev)
    REAL(dp),INTENT(IN) :: pqm1(kbdim,klev)
    REAL(dp),INTENT(IN) :: pxlm1(kbdim,klev),     pxim1(kbdim,klev)
    REAL(dp),INTENT(IN) :: pxtm1(kbdim,klev,ktrac)
    REAL(dp),INTENT(IN) :: pxlte(kbdim,klev),     pxite(kbdim,klev)
    REAL(dp),INTENT(IN) :: papp1(kbdim,klev),    paphp1(kbdim,klevp1)
    REAL(dp),INTENT(IN) :: pqhfla(kbdim)
!0  REAL(dp),INTENT(IN) :: ptkem1(kbdim,klev)
    REAL(dp),INTENT(IN) :: pthvsig(kbdim)

    REAL(dp),INTENT(INOUT) :: ptte(kbdim,klev)
    REAL(dp),INTENT(INOUT) :: pvom(kbdim,klev),      pvol(kbdim,klev)
    REAL(dp),INTENT(INOUT) :: pqte(kbdim,klev)
    REAL(dp),INTENT(INOUT) :: pxtte(kbdim,klev,ktrac)
    REAL(dp),INTENT(INOUT) :: pxtec  (kbdim,klev),  pqtec  (kbdim,klev)

    REAL(dp),INTENT(INOUT) :: prsfc(kbdim),         pssfc(kbdim)  ! OUT
    REAL(dp),INTENT(INOUT) :: ptopmax(kbdim)
    INTEGER, INTENT(INOUT) :: ktype(kbdim)

    INTEGER, INTENT(INOUT) :: ilab(kbdim,klev)                                       ! OUT
    REAL(dp),INTENT(INOUT) :: pcvcbot(kbdim),       pwcape (kbdim)                   ! OUT
    REAL(dp),INTENT(INOUT) :: ptte_cnv(kbdim,klev)                                   ! OUT
    REAL(dp),INTENT(INOUT) :: pvom_cnv(kbdim,klev), pvol_cnv(kbdim,klev)             ! OUT
    REAL(dp),INTENT(INOUT) :: pqte_cnv(kbdim,klev), pxtte_cnv(kbdim,klev,ktrac)      ! OUT

    CHARACTER(LEN=*),PARAMETER :: routine = 'mo_cucall:cucall'

    REAL(dp)::  ztp1(kbdim,klev),         zqp1(kbdim,klev),              &
                zxp1(kbdim,klev),         ztvp1(kbdim,klev),             &
                zup1(kbdim,klev),         zvp1(kbdim,klev),              &
                ztu(kbdim,klev),          zqu(kbdim,klev),               &
                zlu(kbdim,klev),          zlude(kbdim,klev),             &
                zqude(kbdim,klev),                                       &
                zmfu(kbdim,klev),         zmfd(kbdim,klev),              &
                zqsat(kbdim,klev),        zrain(kbdim),                  &
                ua(kbdim)

    INTEGER ::  itopec2(kbdim),           idx(kbdim)
    INTEGER ::  icbot(kbdim),             ictop(kbdim)
    REAL(dp)::  zxtp1(kbdim,klev,ktrac),  zxtu(kbdim,klev,ktrac)
    REAL(dp)::  ztopmax(kbdim)
    LOGICAL ::  locum(kbdim)

    !  Local scalars:
    REAL(dp):: ztmst, zxlp1, zxip1
    INTEGER :: ilevmin, jk, jl, jt

    IF (iconv/=1) THEN
      !
      CALL finish(TRIM(routine),'iconv /=1 not supported in ICON')
      !
    END IF

!-----------------------------------------------------------------------
!*    1.           CALCULATE T,Q AND QS AT MAIN LEVELS
!*                 -----------------------------------
!
!
!100 CONTINUE
  ztmst=ptime_step_len
  DO 120 jk=1,klev
!IBM* NOVECTOR
     DO jl=1,kproma
        ztp1(jl,jk)=ptm1(jl,jk)+ptte(jl,jk)*ztmst
        zqp1(jl,jk)=MAX(0._dp,pqm1(jl,jk)+pqte(jl,jk)*ztmst)
        zxlp1=pxlm1(jl,jk)+pxlte(jl,jk)*ztmst
        zxip1=pxim1(jl,jk)+pxite(jl,jk)*ztmst
        zxp1(jl,jk)=MAX(0._dp,zxlp1+zxip1)
        ztvp1(jl,jk)=ztp1(jl,jk)*(1._dp+vtmpc1*zqp1(jl,jk)-zxp1(jl,jk))
        zup1(jl,jk)=pum1(jl,jk)+pvom(jl,jk)*ztmst
        zvp1(jl,jk)=pvm1(jl,jk)+pvol(jl,jk)*ztmst
     END DO

     CALL prepare_ua_index('cucall',kproma,ztp1(1,jk),idx(1))
     CALL lookup_ua(kproma,idx(1),ua(1))

!IBM* NOVECTOR
     DO jl=1,kproma
        zqsat(jl,jk)=ua(jl)/papp1(jl,jk)
        zqsat(jl,jk)=MIN(0.5_dp,zqsat(jl,jk))
        zqsat(jl,jk)=zqsat(jl,jk)/(1._dp-vtmpc1*zqsat(jl,jk))
     END DO

     DO 1104 jt=1,ktrac
        DO 1102 jl=1,kproma
           zxtp1(jl,jk,jt)=pxtm1(jl,jk,jt)+pxtte(jl,jk,jt)*ztmst
1102    END DO
1104 END DO

120 END DO
  DO 130 jl=1,kproma
     zrain(jl)=0._dp
     locum(jl)=.FALSE.
130 END DO
!
!
!-----------------------------------------------------------------------
!
!*    2.     CALL 'CUMASTR'(MASTER-ROUTINE FOR CUMULUS PARAMETERIZATION)
!*           -----------------------------------------------------------
!
!
!200 CONTINUE
     CALL cumastr(lmfdudv, lmfdd, lmfmid, dlev, cmftau,                &
                  cmfctop, cprcon, cminbuoy, entrpen, nmctop, cevapcu, &
                  pdtime, ptime_step_len,                              &
                  kproma, kbdim, klev, klevp1, klevm1, ilab,           &
                  papp1,                                               &
                  ztp1,     zqp1,     zxp1,     zup1,   zvp1,          &
                  ztvp1,    ktrac,    ldland,                          &
                  zxtp1,    zxtu,                                      &
                  pverv,    zqsat,    pqhfla,                          &
                  paphp1,   pgeo,                                      &
                  ptte,     pqte,     pvom,     pvol,                  &
                  prsfc,    pssfc,    pxtec,                           &
                  pqtec,    zqude,                                     &
                  locum,    ktype,    icbot,    ictop,                 &
                  ztu,      zqu,      zlu,      zlude,                 &
                  zmfu,     zmfd,     zrain,    pthvsig,               &
!--- Included for prognostic CDNC/IC scheme ----------------------------
                  pcvcbot,  pwcape,                                    &
!--- End Included for CDNC/IC ------------------------------------------
                  ptte_cnv, pvom_cnv, pvol_cnv, pqte_cnv, pxtte_cnv    )
!
!
! ------------------------------------------------------------------
!
!*     3.     PRESSURE ALTITUDE OF CONVECTIVE CLOUD TOPS.
!             -------- -------- -- ---------- ----- -----
!
!300 CONTINUE
!
  ilevmin=klev-4
!
  DO 301 jl=1,kproma
     itopec2(jl)=klevp1
301 END DO
!
  DO 303 jk=1,ilevmin
     DO 302 jl=1,kproma
        IF(ilab(jl,jk).EQ.2 .AND. itopec2(jl).EQ.klevp1) THEN
           itopec2(jl)=jk
        END IF
302  END DO
303 END DO
!
  ztopmax(1:kproma) = ptopmax(1:kproma)

  DO 304 jl=1,kproma
     IF(itopec2(jl).EQ.1) THEN
        ptopmax(jl)=papp1(jl,1)
     ELSE IF(itopec2(jl).NE.klevp1) THEN
        ptopmax(jl)=paphp1(jl,itopec2(jl))
     ELSE
        ptopmax(jl)=99999._dp
     END IF
     ptopmax(jl)=MIN(ptopmax(jl),ztopmax(jl))
304 END DO
!
!
!---------------------------------------------------------------------
!
  END SUBROUTINE cucall
  !-------------

END MODULE mo_cucall
