MODULE mo_bgc_constants

!- Description:
!
!  This module contains basic constants and derived constants
!

  USE mo_kind, ONLY: wp

 
  IMPLICIT NONE

  PUBLIC
   


   REAL(wp), PARAMETER :: tmelt = 273.15_wp        ! melting temperature of ice/snow
   
   ! Conversion factors
   REAL(wp), PARAMETER:: kilo = 1.e3_wp
   REAL(wp), PARAMETER:: s2year = 3.1536e7_wp ! seconds to years
   REAL(wp), PARAMETER:: c2gtc = 12._wp*1.e-12_wp ! kmolC to GtC 
   REAL(wp), PARAMETER:: n2tgn = 14._wp*1e-9_wp ! kmolN to TgN
  !>  GRAVITATIONAL ACCELERATION (CONSTANT 9.81 M/S**2)
   REAL(wp), PARAMETER :: g = 9.81_wp
  !>  OCEAN REFERENCE DENSITY [kg m-3]
  REAL(wp), PARAMETER :: rhoref_water = 1025.0_wp
  REAL(wp), PARAMETER:: rmnit=14._wp
   REAL(wp), PARAMETER:: cmh2ms = 6.9722e-07_wp ! cm/hr to m/sec * a constant for piston velocity

  !     -----------------------------------------------------------------
  !*         1. SET HALF PRECISION CONSTANTS
  !             --- ---- --------- ---------
  !
  REAL(wp), PARAMETER ::   SMICR = 1.E-6_wp, PERC  = 1.E-2_wp, THIRD = 1._wp/3._wp
  !
  !     -----------------------------------------------------------------
  !*         3. SET CONVERSION FACTOR SALINITY -> CHLORINITY
  !             ------ ------- -- ---- ------- -- ----------
  !             (AFTER WOOSTER ET AL., 1969)
  !
  REAL(wp), PARAMETER ::  SALCHL = 1._wp/1.80655_wp
  !
  !     -----------------------------------------------------------------
  !*         6. SET COEFFICIENTS FOR APPARENT SOLUBILITY EQUILIBRIUM
  !             OF CALCITE (Mucci 1983, in Z & W-G,2001)
  !             -- ------- ------- ----- --- ----------- -----------
  !         log10ksp = akcc1 + akcc2*t + akcc3/t + akcc4*LOG10(t)         &
  !               + (akcc5 +akcc6*t + akcc7/t) * sqrt(s)                  &
  !               + akcc8*s + akcc9*(s**(1.5))
  REAL(wp), PARAMETER ::  akcc1 = -171.9065_wp,   akcc2 = -0.077993_wp, &
       akcc3 = 2839.319_wp, akcc4 = 71.595_wp, akcc5 = -0.77712_wp,     &
       akcc6 = 0.0028426_wp, akcc7 = 178.34_wp, akcc8 = -0.07711_wp,    &
       akcc9 = 0.0041249_wp
  !
  !     -----------------------------------------------------------------
  !*        6A. SET FRACTION OF ARAGONITE IN BIOGENIC CACO3 PARTICLES
  !             --- -------- --------- -- -------- ----- ---------
  !
  REAL(wp), PARAMETER :: arafra = 0._wp  ! simplification to use solubility product for calcite and
  !               ! adjust it in case of aragonit production
  !
  !     -----------------------------------------------------------------
  !*        6B. FRACTION OF CALCITE IN BIOGENIC CACO3 PARTICLES
  !             -------- ------- -- -------- ----- ---------
  !
  REAL(wp), PARAMETER ::   calfra = 1._wp - arafra
  !
  !
  !     -----------------------------------------------------------------
  !*         7. FACTOR TO GET APPARENT SOLUBILITY PRODUCT FOR
  !             ARAGONIT BY MULTIPLICATION WITH APPARENT SOLUBILITY
  !             PRODUCT FOR CALCIT (CF. BERNER, 1976,
  !             OR BROECKER ET AL., 1982)
  !             -- -------- -- ---- ----- -- ------ --- ------- -----
  !
  REAL(wp), PARAMETER :: aracal = arafra * 1.45_wp + calfra   ! simplification instead of calculating
                                       ! solubility product for aragonit
  !
  !     -----------------------------------------------------------------
  !*         8. SET COEFFICIENTS FOR SEAWATER PRESSURE CORRECTION OF
  !             (AFTER ZEEBE and WOLF-GLADROW (2001)
  !             ------- -------- --- ---------- ----- --- -------- -
  !             in general the pressure dependency for dissoziation contants is given by
  !             ln (K_p / K_0) = - delta V/RT *P + 0.5 delta kappa/RT *P**2
  !             with delta V corresponds to molal volume  := delta V = dev + devt* TC + devt2*TC**2
  !              and delta kappa to compressibility; here neglected
  !             Thus  K_p = K_0 * exp (-(dekv + devkt*TC + devkt2*TC**2)*CP)
  !             with CP = P/RT
  !             K_p = K_0 * exp( -devk - devkt*TC - devkt2*TC**2)*CP
  !             Note that in table A.11.1 all orginal devk (1. column) are negative; herefore, the sign
  !             already is included in the definition of the constant
  !             devt2 for carbon, calcite and aragonit is equal 0.0
  !
  REAL(wp), PARAMETER :: devk1  = 25.50_wp, devk1t = 0.1271_wp, devk2  = 15.82_wp, &
       devk2t = -0.0219_wp, devkb  = 29.48_wp, devkbt = 0.1622_wp,                 &
       devkbt2= 2.6080_wp*1.E-3_wp,  devkw  = 25.60_wp,  devkwt = 0.2324_wp,       &
       devkwt2= -3.6246_wp*1.E-3_wp
  !
  !     -----------------------------------------------------------------
  !*         9. SET COEFFICIENTS FOR PRESSURE CORRECTION OF SOLUBILITY
  !             PRODUCT OF CACO3 CORRESPONDING TO ARAGONITE/CALCITE
  !             RATIO
  !             -----------------------------------------------------------------
  !
  !             Millero (1995)
  !             this formulation is  only for molal volume changes,
  !             it neglects compressibility change
  !
  REAL(wp), PARAMETER ::  devkst = 0.5304_wp
  REAL(wp), PARAMETER ::  devks  = 46_wp * arafra + 48.76_wp * calfra

  !     -----------------------------------------------------------------
  !*        11. SET GAS CONSTANT  for pressure correction after Millero(1995)
  !             --- --------- --- --------
  !             in cm**3 bar / mol K (from Zeebe and Wolf-Gladrow)
  !
  REAL(wp), PARAMETER ::   rgas = 83.131_wp
  !
  !     -----------------------------------------------------------------
  !*        12. SET BORON CONCENTRATION IN SEA WATER
  !*            IN G/KG PER O/OO CL ACCORDING
  !             TO RILEY AND SKIRROW, 1965 (P.250)
  !             -- ----- --- -------- ---- ------- -
  !
  REAL(wp), PARAMETER ::   bor1 = 0.00023_wp
  !
  !     -----------------------------------------------------------------
  !*        13. SET INVERSE OF ATOMIC WEIGHT OF BORON [G**-1]
  !             (USED TO CONVERT SPECIFIC TOTAL BORAT INTO CONCENTRATIONS)
  !             ----- -- ------- -------- ----- ----- ---- ---------------
  !
  REAL(wp), PARAMETER ::   bor2 = 1._wp/10.82_wp

  !
  !     -----------------------------------------------------------------
  !*        15. SET COEFF. FOR 1. DISSOC. OF CARBONIC ACID
  !             DOE (1994)
  !             ------------------------------------------
  !
  REAL(wp), PARAMETER ::  c10 = 2.83655_wp, c11 = -2307.1266_wp, &
       c12 = -1.5529413_wp, c13 = 0.207608410_wp, c14 = 4.0484_wp, &
       c15 = 0.0846834_wp, c16 = -0.00654208_wp, c17 = 0.001005_wp

  !     -----------------------------------------------------------------
  !*        16. SET COEFF. FOR 2. DISSOC. OF CARBONIC ACID
  !             DOE (1994)
  !             ------------------------------------------
  !
  REAL(wp), PARAMETER :: c20 = -9.226508_wp,  c21 = -3351.6106_wp,    &
       c22 = -0.2005743_wp, c23 = 0.106901773_wp,   c24 = 23.9722_wp, &
       c25 = 0.1130822_wp, c26 = -0.00846934_wp,  c27 = 0.001005_wp

  !
  !     -----------------------------------------------------------------
  !*        17. SET COEFF. FOR 1. DISSOC. OF BORIC ACID
  !             DOE (1994) (based on by Dickson 1990)
  !             ---------------------------------------
  !
  REAL(wp), PARAMETER :: cb0 = -8966.90_wp,  cb1 = -2890.53_wp,  cb2 = -77.942_wp, &
       cb3 = 1.728_wp, cb4 = -0.0996_wp, cb5 = 148.0248_wp, cb6 = 137.1942_wp,      &
       cb7 = 1.62142_wp, cb8 = 24.4344_wp, cb9 = 25.085_wp, cb10 =0.2474_wp,        &
       cb11 =0.053105_wp
  !
  !     -----------------------------------------------------------------
  !*        18. SET COEFF. FOR ION PRODUCT OF WATER
  !             DOE (1994)
  !             ---------------------------------------
  !
  REAL(wp), PARAMETER ::  cw0 = 148.96502_wp, cw1 = -13847.26_wp, &
       cw2 = -23.6521_wp, cw3 = 118.67_wp, cw4 = -5.977_wp,       &
       cw5 = 1.0495_wp, cw6 = -0.01615_wp

  !
  !     -----------------------------------------------------------------
  !*        19. SET VOLUMETRIC SOLUBILITY CONSTANTS FOR CO2 IN ML/L
  !             WEISS, R. F. (1974)
  !             CARBON DIOXIDE IN WATER AND SEAWATER: THE SOLUBILITY OF A
  !             NON IDEAL GAS. MARINE CHEMISTRY, VOL. 2, 203-215.
  !     -----------------------------------------------------------------
  REAL(wp), PARAMETER :: c00 = 9345.17_wp,  c01 = -60.2409_wp, c02 = 23.3585_wp, &
       c03 = 0.023517_wp, c04 = -0.00023656_wp, c05 = 0.0047036_wp
  !
  !     -----------------------------------------------------------------
  !*        20. SET VOLUMETRIC SOLUBILITY CONSTANTS FOR O2 IN ML/L
  !             (WEISS, 1970)
  !             ------- ------ --------- --------- --- -- -- ----
  !
  REAL(wp), PARAMETER ::  ox0 = -173.4292_wp, ox1 = 249.6339_wp, ox2 = 143.3483_wp, &
       ox3 = -21.8492_wp, ox4 = -0.033096_wp, ox5 = 0.014259_wp, ox6 = -0.0017_wp
  !
  !     -----------------------------------------------------------------
  !*            SET VOLUMETRIC SOLUBILITY CONSTANTS FOR N2 IN ML/L
  !             WEISS, R. F. (1970) THE SOLUBILITY OF NITROGEN
  !             OXYGEN AND ARGON IN WATER AND SEAWATER.
  !             DEEP-SEA RESEARCH, VOL. 17, 721-735.
  !     -----------------------------------------------------------------
  !
  REAL(wp), PARAMETER ::  an0 = -172.4965_wp, an1 = 248.4262_wp, an2 = 143.0738_wp, &
       an3 = -21.7120_wp, an4 = -0.049781_wp, an5 = 0.025018_wp, an6 = -0.0034861_wp
  !
  !      Constants for laughing gas solubility
  !      (WEISS, 1974, MARINE CHEMISTRY)
  !      --------------------------------------
  !
  REAL(wp), PARAMETER ::  a1 = -62.7062_wp, a2 = 97.3066_wp, a3 = 24.1406_wp, &
       b1 = -0.058420_wp, b2 = 0.033193_wp, b3 = -0.0051313_wp

  REAL(wp), PARAMETER :: atn2o = 3.e-7_wp

  REAL(wp), PARAMETER :: ten   = 10._wp
  ! INVERS OF NORMAL MOLAL VOLUME OF AN IDEAL GAS [CM**3]
  REAL(wp), PARAMETER :: oxyco = 1._wp/22414.4_wp
  !
  REAL(wp), PARAMETER :: invlog10 = 1.0_wp/log(10.0_wp)

END MODULE 
