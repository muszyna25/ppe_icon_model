!>
!! Module contains subroutine cuentrt.
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
MODULE mo_cuentrt

  USE mo_kind,               ONLY: dp
#ifdef __ICON__
  USE mo_physical_constants, ONLY: g=>grav, rd, vtmpc1
  USE mo_echam_conv_nml,     ONLY: centrmax, cmfcmin
#else
  USE mo_constants,          ONLY: g, rd, vtmpc1
  USE mo_cumulus_flux,       ONLY: centrmax, cmfcmin
#endif

  IMPLICIT NONE
  PRIVATE
  PUBLIC  :: cuentrt

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

CONTAINS
  !>
  !!
SUBROUTINE cuentrt(  kproma, kbdim, klev, klevp1, kk,                  &
           ptenh,    pqenh,    pqte,     paphp1,                       &
           klwmin,   ldcum,    ktype,    kcbot,    kctop0,             &
           ppbase,   pmfu,     pentr,                                  &
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
INTEGER  :: klwmin(kbdim),           ktype(kbdim),                     &
            kcbot(kbdim),            kctop0(kbdim)
LOGICAL  :: ldcum(kbdim)
!
REAL(dp) :: pdmfen(kbdim),           pdmfde(kbdim)
!
LOGICAL  :: llo1,llo2
!
INTEGER  :: jl, iklwmin
REAL(dp) :: zrg, zrrho, zdprho, zpmid, zentr, zentest

!
!----------------------------------------------------------------------
!
!*    1.           CALCULATE ENTRAINMENT AND DETRAINMENT RATES
!                  -------------------------------------------
!
!
!
!*    1.1          SPECIFY ENTRAINMENT RATES FOR SHALLOW CLOUDS
!                  --------------------------------------------
!

!
!*    1.2          SPECIFY ENTRAINMENT RATES FOR DEEP CLOUDS
!                  -----------------------------------------
!
  zrg=1._dp/g
  DO 125 jl=1,kproma
     ppbase(jl)=paphp1(jl,kcbot(jl))
     zrrho=(rd*ptenh(jl,kk+1)*(1._dp+vtmpc1*pqenh(jl,kk+1)))           &
                    /paphp1(jl,kk+1)
     zdprho=(paphp1(jl,kk+1)-paphp1(jl,kk))*zrg
     zpmid=0.5_dp*(ppbase(jl)+paphp1(jl,kctop0(jl)))
     zentr=pentr(jl)*pmfu(jl,kk+1)*zdprho*zrrho
     llo1=kk.LT.kcbot(jl).AND.ldcum(jl)
     pdmfde(jl)=MERGE(zentr,0._dp,llo1)
     llo2=llo1.AND.ktype(jl).EQ.2.AND.                                 &
                   (ppbase(jl)-paphp1(jl,kk).LT.0.2e5_dp.OR.           &
                    paphp1(jl,kk).GT.zpmid)
     pdmfen(jl)=MERGE(zentr,0._dp,llo2)
     iklwmin=MAX(klwmin(jl),kctop0(jl)+2)
     llo2=llo1.AND.(ktype(jl).EQ.1 .OR. ktype(jl).EQ.3) .AND.          &
                   (kk.GE.iklwmin.OR.paphp1(jl,kk).GT.zpmid)
     IF(llo2) pdmfen(jl)=zentr
     IF(llo2 .AND. pqenh(jl,kk+1).GT.1.E-5_dp) THEN
        pmfu(jl,kk+1) = MAX(pmfu(jl,kk+1),cmfcmin)
        zentest = MAX(pqte(jl,kk),0._dp)/pqenh(jl,kk+1)
        zentest = MIN(centrmax,zentest/(pmfu(jl,kk+1)*zrrho))
        pdmfen(jl) = zentr+zentest*pmfu(jl,kk+1)*zrrho*zdprho
     ENDIF
125 END DO
!
    RETURN
  END SUBROUTINE cuentrt

END MODULE mo_cuentrt
