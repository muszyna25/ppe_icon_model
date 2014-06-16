! this is a workaround to initialize the echam physics packages for use
! in iconam
! configure_echam_convection sets the value for two parameters:
! echam_conv_config%nmctop
! echam_conv_config%cevapcu 
!
! nmctop was computed useing vct_a and vct_b which we need to avoid when using
! the sleve coordinates.  this value is below computed using a US standard
! atmospheric profile. 

!>
!!  Initialize the physical schemes at start time
!!
!! @par Revision History
!! First implementations by Kristina Froehlich, DWD, 2010-070-20
!! Include initialitation of SST for APE experiments by P. Ripodas, DWD,2010-11
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

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_nhecham_conv_config

  USE mo_kind,                ONLY: wp
  USE mo_impl_constants,      ONLY: SUCCESS
  USE mo_physical_constants,  ONLY: grav, rd, p0sl_bg, dtdz_standardatm
  USE mo_exception,           ONLY: message, finish ,message_text, print_value
  USE mo_vertical_coord_table,ONLY: vct_a, vct
  USE mo_echam_conv_config,   ONLY: echam_conv_config

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

  PUBLIC  :: config_nhecham_conv

CONTAINS


SUBROUTINE config_nhecham_conv( nlev ) 

  INTEGER, INTENT(IN) :: nlev ! number of full levels
  REAL(wp)            :: pref(nlev)
  REAL(wp)            :: cceta(nlev)
  REAL(wp)            :: ztmp

  INTEGER :: jk, istat
!  INTEGER, INTENT(IN) :: nshift ! shift with respect to global grid

! Reference atmosphere parameters
  REAL(wp), PARAMETER :: htropo = 11000._wp       ! [m]    tropopause height
  REAL(wp), PARAMETER :: t00    = 288.15_wp       ! [m]    temperature at sea level
  REAL(wp) :: ttropo, ptropo, temp, zfull

  CHARACTER(LEN=*),PARAMETER :: &
           routine = 'mo_echam_conv_config:config_echam_convection'

!--------------------------------------------------------------
!>reference pressure according to U.S. standard atmosphere
! (with the caveat that the stratosphere is assumed isothermal, which does not hurt
!  because pref is used for determining model level indices referring to pressures
!  >= 60 hPa)
!--------------------------------------------------------------
ttropo = t00 + dtdz_standardatm*htropo
! p0sl_bg is sea level pressure
ptropo = p0sl_bg*(ttropo/t00)**(-grav/(rd*dtdz_standardatm))
DO jk = nlev, 1, -1
  zfull = 0.5_wp*(vct_a(jk) + vct_a(jk+1))
  IF (zfull < htropo) THEN
    temp = t00 + dtdz_standardatm*zfull
    pref(jk) = p0sl_bg*(temp/t00)**(-grav/(rd*dtdz_standardatm))
  ELSE
    pref(jk) = ptropo*EXP(-grav*(zfull-htropo)/(rd*ttropo))
  ENDIF
END DO 
  
! search for the 300 hPa level

CALL message('','')
CALL print_value('lowest model level for cloud base '//&
                 'of mid level convection: nmctop = ', &
                 echam_conv_config%nmctop)
CALL message('','')
DO jk = 1, nlev
  echam_conv_config%nmctop = jk
  IF(pref(jk).GE.30000.0_wp) THEN
    write(*,*) "LOOP EXIT; PREF(JK) AT JK=",jk," IS = ",pref(jk)
    write(*,*) "TROPO H[m] = ",htropo,"; temp[K] at sea level is = ",t00
    EXIT
  ENDIF
END DO

!--------------------------------------
! Set evaporation coefficient for kuo0
!--------------------------------------

! ...from echam_conv_config...
ALLOCATE( echam_conv_config%cevapcu(nlev),STAT=istat )
IF (istat/=SUCCESS) CALL finish(TRIM(routine),'allocation of cevapcu failed')

DO jk = 1,nlev
   cceta(jk) = pref(jk)/p0sl_bg
   ztmp = 1.E3_wp/(38.3_wp*0.293_wp)*SQRT(cceta(jk))
   echam_conv_config%cevapcu(jk) = 1.93E-6_wp*261._wp*SQRT(ztmp)*0.5_wp/grav
   write(*,*) "VALUE OF ECHAM_CONV_CONFIG IS ", echam_conv_config%cevapcu(jk)
END DO

!------------------------------------------------------------------------
! Print the configuration on stdio
!------------------------------------------------------------------------
CALL message('','')
CALL message('','------- configuration of the ECHAM convection scheme --------')

SELECT CASE (echam_conv_config% iconv)
CASE(1); CALL message('','--- iconv = 1 -> Convection: Nordeng (default)')
CASE(2); CALL message('','--- iconv = 2 -> Convection: Tiedtke')
CASE(3); CALL message('','--- iconv = 3 -> Convection: Hybrid')
CASE default
  WRITE(message_text,'(a,i0,a)') 'iconv = ',echam_conv_config% iconv, &
                                 ' is not supported'
  CALL finish(TRIM(routine),message_text)
END SELECT

SELECT CASE(echam_conv_config% ncvmicro)
CASE (0); CALL message('','--- ncvmicro = 0')
CASE DEFAULT
  CALL finish(TRIM(routine),'ncvmicro > 0 not yet supported in ICON')
END SELECT

CALL print_value(' lmfpen   ', echam_conv_config% lmfpen)
!CALL print_value(' lmfmid   ', echam_conv_config% lmfmid)
!CALL print_value(' lmfdd    ', echam_conv_config% lmfdd)
!CALL print_value(' lmfdudv  ', echam_conv_config% lmfdudv)

!CALL print_value(' cmftau   ', echam_conv_config% cmftau)
!CALL print_value(' cmfctop  ', echam_conv_config% cmfctop)
!CALL print_value(' cprcon   ', echam_conv_config% cprcon)
!CALL print_value(' cminbuoy ', echam_conv_config% cminbuoy)
!CALL print_value(' entrpen  ', echam_conv_config% entrpen)
!CALL print_value(' dlev     ', echam_conv_config% dlev)

CALL print_value(' ECHAM_CONV_CONFIG%NMCTOP ', echam_conv_config%nmctop)
CALL message('','---------------------------')
CALL message('','')

END SUBROUTINE config_nhecham_conv

END MODULE mo_nhecham_conv_config
