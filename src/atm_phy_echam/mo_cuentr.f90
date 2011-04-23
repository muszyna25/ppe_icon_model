!>
!! Module contains subroutine cuentr.
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
#ifdef __xlC__
@PROCESS HOT
#else
#define FSEL(a,b,c) MERGE(b,c,(a) >= 0._dp)
#endif

MODULE mo_cuentr

  USE mo_kind,               ONLY: dp
#ifdef __ICON__
  USE mo_physical_constants, ONLY: g=>grav, rd, vtmpc1
  USE mo_echam_conv_params,  ONLY: centrmax, cmfcmin
#else
  USE mo_constants,          ONLY: g, rd, vtmpc1
  USE mo_cumulus_flux,       ONLY: centrmax, cmfcmin
#endif
#ifdef _PROFILE
  USE mo_profile,        ONLY: trace_start, trace_stop
#endif

  IMPLICIT NONE
  PRIVATE
  PUBLIC  :: cuentr

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

CONTAINS
  !>
  !!
SUBROUTINE cuentr(   kproma, kbdim, klev, klevp1, kk,                  &
           ptenh,    pqenh,    pqte,     paphp1,                       &
           klwmin,   ldcum,    ktype,    kcbot,    kctop0,             &
           ppbase,   pmfu,     pentr,    podetr,                       &
           khmin,    pgeoh,                                            &
           pdmfen,   pdmfde)
!
!          M.TIEDTKE         E.C.M.W.F.     12/89
!
!          PURPOSE.
!          --------
!          THIS ROUTINE CALCULATES ENTRAINMENT/DETRAINMENT RATES
!          FOR UPDRAFTS IN CUMULUS PARAMETERIZATION
!
!          INTERFACE
!          ---------
!
!          THIS ROUTINE IS CALLED FROM *CUASC*.
!          INPUT ARE ENVIRONMENTAL VALUES T,Q ETC
!          AND UPDRAFT VALUES T,Q ETC
!          IT RETURNS ENTRAINMENT/DETRAINMENT RATES
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
INTEGER, INTENT (IN) :: kbdim, klev, klevp1, kproma, kk
!
REAL(dp) :: ptenh(kbdim,klev),       pqenh(kbdim,klev),                &
            paphp1(kbdim,klevp1),                                      &
            pmfu(kbdim,klev),        pqte(kbdim,klev),                 &
            pentr(kbdim),            ppbase(kbdim)
REAL(dp) :: podetr(kbdim,klev)
REAL(dp) :: pgeoh (kbdim,klev)
INTEGER  :: khmin (kbdim)
INTEGER  :: klwmin(kbdim),           ktype(kbdim),                     &
            kcbot(kbdim),            kctop0(kbdim)
LOGICAL  :: ldcum(kbdim)
!
REAL(dp) :: pdmfen(kbdim),           pdmfde(kbdim)
!
LOGICAL  :: llo1,llo2
!
INTEGER  :: jl, ikt, ikh, iklwmin, n, ncnt
REAL(dp) :: zrg, zpmid, zentr, zentest, zzmzk, ztmzk, zorgde, zarg
REAL(dp) :: zrrho(kbdim),zdprho(kbdim)
INTEGER  :: icond1(kbdim),icond2(kbdim),icond3(kbdim),idx(kbdim)
!
!----------------------------------------------------------------------
!
!*    1.           CALCULATE ENTRAINMENT AND DETRAINMENT RATES
!                  -------------------------------------------
!
!100 CONTINUE
!
!
!*    1.1          SPECIFY ENTRAINMENT RATES FOR SHALLOW CLOUDS
!                  --------------------------------------------
!
!110 CONTINUE
!
!
!*    1.2          SPECIFY ENTRAINMENT RATES FOR DEEP CLOUDS
!                  -----------------------------------------
!
!120 CONTINUE

#ifdef _PROFILE
  CALL trace_start ('cuentr', 41)
#endif
!
  zrg=1._dp/g
!IBM* NOVECTOR
  DO 125 jl=1,kproma
     ppbase(jl) = paphp1(jl,kcbot(jl))
     zrrho(jl)  = (rd*ptenh(jl,kk+1)*(1._dp+vtmpc1*pqenh(jl,kk+1)))    &
                    / paphp1(jl,kk+1)
     zdprho(jl) = (paphp1(jl,kk+1)-paphp1(jl,kk))*zrg
     zpmid      = 0.5_dp*(ppbase(jl)+paphp1(jl,kctop0(jl)))
     icond1(jl) = INT(FSEL(zpmid-paphp1(jl,kk),0._dp,1._dp))
     icond2(jl) = INT(FSEL(0.2e5_dp - (ppbase(jl)-paphp1(jl,kk)),0._dp,1._dp))
     icond3(jl) = INT(FSEL(1.e-5_dp-pqenh(jl,kk+1),0._dp,1._dp))
     pdmfde(jl)=0._dp
     pdmfen(jl)=0._dp
     podetr(jl,kk)=0._dp
125 END DO

  ncnt = 0
  DO jl=1,kproma
     llo1=kk.LT.kcbot(jl).AND.ldcum(jl)
     IF (llo1) THEN
        ncnt = ncnt+1
        idx(ncnt) = jl
     END IF
  END DO

!CDIR NODEP
!IBM* ASSERT(NODEPS)
  DO 126 n=1,ncnt
     jl = idx(n)
     zentr=pentr(jl)*pmfu(jl,kk+1)*zdprho(jl)*zrrho(jl)
     pdmfde(jl)=zentr

     llo2=ktype(jl).EQ.2.AND.(icond2(jl).GT.0.OR.icond1(jl).GT.0)
     IF (llo2) pdmfen(jl)=zentr

     iklwmin=MAX(klwmin(jl),kctop0(jl)+2)
     llo2=ktype(jl).EQ.3 .AND. kk .GE. iklwmin
     IF(llo2) pdmfen(jl)=zentr

     IF(llo2 .AND. icond3(jl).GT.0) THEN
        pmfu(jl,kk+1) = MAX(pmfu(jl,kk+1),cmfcmin)
        zentest = MAX(pqte(jl,kk),0._dp)/pqenh(jl,kk+1)
        zentest = MIN(centrmax,zentest/(pmfu(jl,kk+1)*zrrho(jl)))
        pdmfen(jl) = zentr+zentest*pmfu(jl,kk+1)*zrrho(jl)*zdprho(jl)
     ENDIF

     llo2=ktype(jl).EQ.1 .AND.(kk.GE.iklwmin.OR.icond1(jl).GT.0)
     IF(llo2) pdmfen(jl)=zentr
126 END DO
!
!    organized detrainment, detrainment starts at khmin
!
  IF (kproma > 0) THEN
! kproma>0: mpuetz's workaround to avoid combining of two loops
!CDIR NODEP
!IBM* ASSERT(NODEPS)
     DO 127 n=1,ncnt
        jl = idx(n)
        llo2=ktype(jl).EQ.1
        IF(llo2.AND.kk.LE.khmin(jl).AND.kk.GE.kctop0(jl)) THEN
           ikt=kctop0(jl)
           ikh=khmin(jl)
           IF(ikh.GT.ikt) THEN
              zzmzk  =-(pgeoh(jl,ikh)-pgeoh(jl,kk))*zrg
              ztmzk  =-(pgeoh(jl,ikh)-pgeoh(jl,ikt))*zrg
              zarg  =3.1415_dp*(zzmzk/ztmzk)*0.5_dp
              zorgde=TAN(zarg)*3.1415_dp*0.5_dp/ztmzk
              zdprho(jl)=(paphp1(jl,kk+1)-paphp1(jl,kk))*(zrg*zrrho(jl))
              podetr(jl,kk)=MIN(zorgde,centrmax)*pmfu(jl,kk+1)*zdprho(jl)
           ENDIF
        ENDIF
127  END DO
  END IF
!
#ifdef _PROFILE
  CALL trace_stop ('cuentr', 41)
#endif
!
  END SUBROUTINE cuentr

END MODULE mo_cuentr
