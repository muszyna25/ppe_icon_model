! $RCSfile$
! $Revision$ $Date$

!>
!! <Short description of module for listings and indices>
!! Description:
!!**** *CUCALCLPI*  ROUITINE FOR LPI COMPUTATION
!!
!!  @author  Guido Schroeder               02/02/2021    
!!
!!     PURPOSE
!!     -------

!!     THIS ROUTINE COMPUTES THE LIGHTNIN POTENTIAL INDEX AS IN 
!!     LYNN AND YAIR (2010). IT ALSO MAKES USE OF SOME IDEAS OF LOPEZ (2016)
!!     FOR THE LIGHTNING FLASH DENSITY COMPUTATION.


!!     METHOD
!!     -------

!!     IN CONTRAST TO THE ORIGINAL PAPER OF LYNN AND YAIR THE
!!     UPDRAFT PROFILE OF THE CONVECION SCHEME IS USED TO 
!!     COMPUTE THE LPI.
!!     THE UPDRAFT TEMPERATURE, LIQUID AND FROZEN WATER AS WELL AS
!!     THE KINETIC ENERGY OF THE UPDRAFT ARE USED.

!! @par Revision History
!! Initial implementation into ICON  by Guido Schroeder, DWD (2021-02-02)
!!
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_cucalclpi
  
  USE mo_kind   ,ONLY: JPRB=>wp     , &
    &                  jpim=>i4
  
  USE mo_cuparameters , ONLY :                                   &
    & rg, rd, rcpd
  
  USE mo_cufunctions, ONLY: foealfa, foeldcpm

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: cucalclpi, cucalcmlpi

CONTAINS

SUBROUTINE cucalclpi(klon, klev, ktype, ztu, zlu, zkineu, zmflxs        &
      &            , zten, pap, zdgeoh, ldland, lpi)

! Code Description:
!     PARAMETER     DESCRIPTION                                   UNITS
!     ---------     -----------                                   -----
!     INPUT PARAMETERS (INTEGER):

!    *KLON*         NUMBER OF GRID POINTS PER PACKET
!    *KLEV*         NUMBER OF LEVELS
!    *KTYPE*        TYPE OF CONVECTION
!                       1 = PENETRATIVE CONVECTION
!                       2 = SHALLOW CONVECTION
!                       3 = MIDLEVEL CONVECTION

!     INPUT PARAMETERS (LOGICAL)

!    *LDLAND*       LAND SEA MASK (.TRUE. FOR LAND)

!     INPUT PARAMETERS (REAL)

!    *ZTU*          TEMPERATURE IN UPDRAFTS                         K
!    *ZLU*          LIQUID WATER CONTENT IN UPDRAFTS              KG/KG
!    *ZKINEU*       KINETIC ENERGY IN UPDRATFS                    M2/S2
!    *ZMFLXS*       CONVECTIVE SNOW FLUX                          KG/(M2*S)
!    *ZTEN*         ENVIRONMENT TEMPERATURE (T+1)                 K
!    *PAP*          PRESSURE ON FULL LEVELS                       PA
!    *zdgeoh*       geopot thickness on full levels               M2/S2

!    OUTPUT PARAMETERS (REAL):

!    *LPI*          LIGHTNING POTENTIAL INDEX AS IN LYNN AND YAIR (2010) J/KG

IMPLICIT NONE

INTEGER(KIND=jpim),INTENT(in)  :: klon
INTEGER(KIND=jpim),INTENT(in)  :: klev
INTEGER(KIND=jpim),INTENT(in)  :: ktype(klon)
REAL(KIND=jprb)   ,INTENT(in)  :: ztu(klon,klev) 
REAL(KIND=jprb)   ,INTENT(in)  :: zlu(klon,klev) 
REAL(KIND=jprb)   ,INTENT(in)  :: zmflxs(klon,klev+1) 
REAL(KIND=jprb)   ,INTENT(in)  :: zten(klon,klev) 
REAL(KIND=jprb)   ,INTENT(in)  :: pap(klon,klev) 
REAL(KIND=jprb)   ,INTENT(in)  :: zdgeoh(klon,klev)
REAL(KIND=jprb)   ,INTENT(in)  :: zkineu(klon,klev)
LOGICAL           ,INTENT(in)  :: ldland(klon) 
REAL(KIND=jprb)   ,INTENT(out) :: lpi(klon)

REAL(KIND=jprb)   :: zrho(klon,klev)
INTEGER (KIND=jpim)   :: kland(klon)

! Parameters as in Lopez 2016
! coefficient to split solid water mass flux in snow and graupel
! beta(1) for land, beta(2) for sea - Takahashi 2006 (mentioned
! in Lopez 2016)
REAL(KIND=jprb), PARAMETER :: beta(2)= [0.7_jprb , 0.45_jprb ]  
REAL(KIND=jprb), PARAMETER :: Vgraup=3.0_jprb ! 3   m/s fall speed for graupel
REAL(KIND=jprb), PARAMETER :: Vsnow=0.5_jprb  ! 0.5 m/s fall speed for snow

! Auxiliary fields as in Lopez 2016 and Lynn and Yair 2010
REAL(KIND=jprb) :: zqIce, zqLiquid, zqGraup, zqSnow, zEps, zQi, zdz

REAL(KIND=jprb) :: unitVolume(klon) ! the volume used for the integral - region
                                  ! beween 0 and -20 celsius

INTEGER(KIND=jpim) :: jk, jl

! Make land sea mask for beta
  kland=2_jpim
  WHERE (ldland) kland=1_jpim

  zrho=pap/(zten*Rd+1E-10_jprb)

  ! Now compute the LPI (Lynn and Yair 2010)
  LPI=0.0_jprb
  unitVolume=0.0_jprb
  DO JK = 1, klev
    DO JL = 1, klon
      ! only for deep convection LPI is computed
      IF (ktype(JL) == 1) THEN
        IF (ztu(JL,JK) <=273.15_jprb .AND. ztu(JL,JK) >=273.15_jprb-20._jprb) THEN
        ! compute graupel and snow mixing ratio from the massflux
        ! of frozen precip - splitting it into snow and graupel with beta.
        ! See Eq. 1 and 2 in Lopez (2016) 
        ! "A lightning parameterization for the ECMWF integrated forecasting system"
          zqGraup=beta(kland(JL))*zmflxs(JL, JK)/zrho(JL,JK)/Vgraup
          zqSnow=(1-beta(kland(JL)))*zmflxs(JL,JK)/zrho(JL,JK)/Vsnow
          zqLiquid=zlu(JL,JK)*foealfa(ztu(JL,JK))
          zqIce=zlu(JL,JK)*(1._jprb-foealfa(ztu(JL,JK)))
        ! eq. 3 in Lynn and Yair (2010)
          zQI=zqGraup*(SQRT(zqSnow*zqGraup)/(zqSnow+zqGraup+1E-20_jprb)+&
                    SQRT(zqIce*zqGraup)/(zqIce+zqGraup+1E-20_jprb))
        ! eq. 2 in Lynn and Yair (2010)
          zeps=2._jprb*sqrt(zQI*zqLiquid)/(zQi+zqLiquid+1E-20_jprb)
          zdz=zdgeoh(jl,jk)/rg
          LPI(JL) = LPI(JL)+zeps*2._jprb*zkineu(JL, JK)*zdz
          unitVolume(JL) = unitVolume(JL)+zdz
        ENDIF
      ENDIF
    ENDDO
    !ENDDO
  ENDDO

  lpi=lpi/(unitVolume+1E-20_jprb)

!!  ! We use 15 km as vertial unit volume - after all it is just a scaling factor
!!  ! Note - In the paper the LPI is an average oer the unit volume. If we
!!  ! applied that here properly we would need to multiply with sigma
!!  ! (the updraft area proportion in the cloud). That would mean the numbers
!!  ! get smaller the larger the grid spacing. That is not necessarily
!!  ! desirable.
!!  ! In conclusion: The LPI here is representative for the updraft in the cloud
!!  ! only.
!!  lpi=lpi/15000._jprb 

END SUBROUTINE cucalclpi
 
SUBROUTINE CUCALCMLPI(klon, klev, lpi, zten, zqen, pap, paph, koi, mlpi)
! Code Description:
! Computes a modified LPI using LPI as in Lynn and Yair and KOI.
!
!     PARAMETER     DESCRIPTION                                   UNITS
!     ---------     -----------                                   -----
!     INPUT PARAMETERS (INTEGER):

!    *KLON*         NUMBER OF GRID POINTS PER PACKET
!    *KLEV*         NUMBER OF LEVELS
!    *KTYPE*        TYPE OF CONVECTION
!                       1 = PENETRATIVE CONVECTION
!                       2 = SHALLOW CONVECTION
!                       3 = MIDLEVEL CONVECTION


!     INPUT PARAMETERS (REAL)

!    *LPI*          LIGHTNING POTENTIAL INDEX AS IN LYNN AND YAIR (2010) J/KG
!    *ZTEN*         ENVIRONMENT TEMPERATURE (T+1)                 K
!    *ZQEN*         ENVIRONMENT SPEC. HUMIDITY (T+1)              KG/KG
!    *PAP*          PRESSURE ON FULL LEVELS                       PA
!    *PAPH*         PRESSURE ON HALF LEVELS                       PA

!    OUTPUT PARAMETERS (REAL):

!    *KOI*          CONVECTION INDEX (VERT. GRADIENT OF EQUIALENT POT.TEMP) K
!    *MLPI*         MODIFIED LIGHTNING POTENTIAL INDEX (COMBINATIN WITH KOI) J/KG

IMPLICIT NONE
  
INTEGER(KIND=jpim),INTENT(in) :: klon
INTEGER(KIND=jpim),INTENT(in) :: klev
REAL(KIND=jprb),INTENT(in)    :: zten(klon,klev) 
REAL(KIND=jprb),INTENT(in)    :: zqen(klon,klev) 
REAL(KIND=jprb),INTENT(in)    :: pap(klon,klev)
REAL(KIND=jprb),INTENT(in)    :: paph(klon,klev+1)
REAL(KIND=jprb),INTENT(in)    :: lpi(klon)
REAL(KIND=jprb),INTENT(out)    :: mlpi(klon)
REAL(KIND=jprb),INTENT(out)    :: koi(klon)
! Equivalent temperatures in 600 and 900 hPa
! Note - the average of 500 to 700 hPa is computed for 600hPa
!        the average below 800 hPa is computed for 900 hPa
REAL(KIND=jprb) :: te
! Equivalent potential temperatures in 600 and 900 hPa
REAL(KIND=jprb) :: thetae600(klon), thetae900(klon), thetae(klon, klev)
! The total pressure difference for the integral
REAL(KIND=jprb) :: deltap600(klon), deltap900(klon)

!* Coefficients for the MLPI formula
!* They were computed in an optimization process
REAL(KIND=jprb) :: fa(klon),fb(klon)
REAL(KIND=jprb), PARAMETER :: fe=0.391040135975501_jprb
REAL(KIND=jprb), PARAMETER :: fd=1.26663689775482_jprb
REAL(KIND=jprb), PARAMETER :: fg=4.728883_jprb !!2.58697023913754_jprb
REAL(KIND=jprb), PARAMETER :: fh=0.625209762557336_jprb
REAL(KIND=jprb), PARAMETER :: fi=1.825632_jprb !! 3.0022767963889_jprb
REAL(KIND=jprb), PARAMETER :: fj=-2.185377_jprb !!-1.3288895756172_jprb
REAL(KIND=jprb), PARAMETER :: fk=0.309079838495018_jprb
REAL(KIND=jprb), PARAMETER :: fbmax=14.04617_jprb !!2.80923302891783_jprb

INTEGER(KIND=jpim) :: jk, jl
  ! compute equivalent potential temperature
  DO JK = 1_jpim, klev
    DO JL = 1_jpim, klon
       te = zten(JL,JK)+foeldcpm(zten(JL,JK)+1E-20)*zqen(JL,JK)
       thetae(JL,JK)=te*(100000.D0/(pap(JL,JK)+1E-10))**(rd/rcpd)
    ENDDO
  ENDDO

  ! Now compute KOI
  thetae900=0.0_jprb
  deltap900=0.0_jprb
  thetae600=0.0_jprb
  deltap600=0.0_jprb
  DO JK = 1_jpim, klev
    DO JL = 1_jpim, klon
      IF (pap(JL,JK) > 80000.D0) THEN
        thetae900(JL) = thetae900(JL)                                   &
 &                       +(paph(JL,JK+1)-paph(JL,JK))*thetae(JL,JK)
        deltap900(JL) =  deltap900(JL)+(paph(JL,JK+1)-paph(JL,JK))
      ELSE IF (pap(JL,JK) < 70000.D0 .AND. pap(JL,JK) > 50000.D0)  THEN
        thetae600(JL) = thetae600(JL)                                   &
 &                       +(paph(JL,JK+1)-paph(JL,JK))*thetae(JL,JK)
        deltap600(JL) =  deltap600(JL)+(paph(JL,JK+1)-paph(JL,JK))
      ENDIF
    ENDDO
  ENDDO
  thetae900 = thetae900/(deltap900+1E-20)
  thetae600 = thetae600/(deltap600+1E-20)
  KOI=thetae600-thetae900
  ! Over mountains KOI cannot be properly computed - here we set KOI to zero to
  ! leave LPI unchanged.
  WHERE (deltap900 < 1E-20_jprb) KOI=0_jprb

  ! Compute the modified LPI
  fa = fg*LPI**fh
  ! we require a >= b
  fb = min(fa, fbmax*(1_jprb+tanh(fi*(LPI**fk+fj)))/2_jprb)
  MLPI= fb+(1_jprb+tanh(-fe*(KOI+fd)))/2_jprb*(fa-fb)
                     

END SUBROUTINE CUCALCMLPI

END MODULE mo_cucalclpi
