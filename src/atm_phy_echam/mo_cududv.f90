!>
!! Module contains subroutine cududv.
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
MODULE mo_cududv

  USE mo_kind,               ONLY: dp
#ifdef __ICON__
  USE mo_physical_constants, ONLY: g=>grav
#else
  USE mo_constants,          ONLY: g
#endif

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: cududv

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

CONTAINS
  !>
  !!
SUBROUTINE cududv(   kproma,   kbdim,    klev,     klevp1,             &
           ktopm2,   ktype,    kcbot,    paphp1,   ldcum,              &
           puen,     pven,     pvom,     pvol,                         &
           puu,      pud,      pvu,      pvd,                          &
           pmfu,     pmfd,     pvom_cnv, pvol_cnv)
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
INTEGER, INTENT (IN) :: kproma, kbdim, klev, klevp1, ktopm2
!
REAL(dp):: puen(kbdim,klev),        pven(kbdim,klev),                  &
           pvol(kbdim,klev),        pvom(kbdim,klev),                  &
           paphp1(kbdim,klevp1)
REAL(dp):: puu(kbdim,klev),         pud(kbdim,klev),                   &
           pvu(kbdim,klev),         pvd(kbdim,klev),                   &
           pmfu(kbdim,klev),        pmfd(kbdim,klev)
REAL(dp),INTENT(OUT) :: pvom_cnv(kbdim,klev), pvol_cnv(kbdim,klev)
INTEGER :: ktype(kbdim),            kcbot(kbdim)
LOGICAL :: ldcum(kbdim)
!
REAL(dp):: zmfuu(kbdim,klev),       zmfdu(kbdim,klev),                 &
           zmfuv(kbdim,klev),       zmfdv(kbdim,klev)
!
INTEGER :: jl, jk, ik, ikb
REAL(dp):: zzp, zdudt, zdvdt
!

  pvom_cnv(1:kproma,:) = 0._dp
  pvol_cnv(1:kproma,:) = 0._dp

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
              zdudt=(g/(paphp1(jl,jk+1)-paphp1(jl,jk)))*               &
                          (zmfuu(jl,jk+1)-zmfuu(jl,jk)+                &
                           zmfdu(jl,jk+1)-zmfdu(jl,jk))
              zdvdt=(g/(paphp1(jl,jk+1)-paphp1(jl,jk)))*               &
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
              zdudt=-(g/(paphp1(jl,jk+1)-paphp1(jl,jk)))*              &
                           (zmfuu(jl,jk)+zmfdu(jl,jk))
              zdvdt=-(g/(paphp1(jl,jk+1)-paphp1(jl,jk)))*              &
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

END MODULE mo_cududv
