!>
!! Parameters for Ritter-Geleyn radiation, taken from DWD's GME/COSMO
!!
!! @author <Thorsten Reinhardt, AGeoBw, Offenbach>
!! @author <Bodo Ritter, DWD, Offenbach>
!!
!!
!! @par Revision History
!! Initial Release by Thorsten Reinhardt (2010-12-02)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!

MODULE mo_radiation_rg_par

  USE mo_kind,                 ONLY: wp
  
  IMPLICIT NONE

  PUBLIC


!      INCLUDE "gme_par_rad.h"
!     *PARAMETERS* FOR SPECTRAL RESOLUTION OF RADIATIVE TRANSFER CODE
!
  REAL(wp), PARAMETER :: zepmu0 = 1.E-06_wp ! security constant for cos (zenith angle)
  
  INTEGER, PARAMETER :: &! Parameters for spectral resolution of radiative transfer code
    & JPSOL =3, &        ! Number of solar spectral intervals
    & JPTHER=5, &        ! Number of thermal spectral intervals
    & JPSPEC=8           ! JPSOL+JPTHER! Total number of spectral intervals

      INTEGER ncgas(jpspec,3)  ! number of coefficients for each
                               ! interval and gas (maximum=7)
      INTEGER nfast(jpspec)  ! control variable for choice between
                             ! ESFT/FESFT in each interval
!     *gme_com_gas.h : absorption properties of atmospheric gases
!                  (..,1 = h2o ; ..,2 = co2 ; ..,3 = o3)
      REAL(wp) :: coai  (7,jpspec,3)    ! weigthing coefficients
      REAL(wp) :: cobi  (7,jpspec,3)   ! absorption coefficients
      REAL(wp) :: coali (7,jpspec,3)    ! pressure    correction coefficients
      REAL(wp) :: cobti (7,jpspec,3)    ! temperature correction coefficients

!     *COMGRZ*  : Limits of spectral intervals (for information only)
!      REAL(wp) :: grenze(2,2,jpspec)
  !             WMIN1      , WMAX1       , WMIN2       , WMAX2
!  DATA grenze /  1.5300_wp ,   4.6420_wp , 999._wp     ,  0._wp    , &
!                 0.7000_wp ,   1.5300_wp , 999._wp     ,  0._wp    , &
!                 0.2451_wp ,   0.7000_wp , 999._wp     ,  0._wp    , &
!                20.0000_wp , 104.5150_wp , 999._wp     ,  0._wp    , &
!                12.5000_wp ,  20.0000_wp , 999._wp     ,  0._wp    , &
!                 8.3333_wp ,   9.0090_wp ,  10.3093_wp , 12.5000_wp, &
!                 9.0090_wp ,  10.3093_wp , 999._wp     ,  0._wp    , &
!                 4.6420_wp ,   8.3333_wp , 999._wp     ,  0._wp      /

!     DATEN ZU DEN OPTISCHEN EIGENSCHAFTEN FUER FLUESSIGWASSER
!     FUER ALLE 8 SPEKTRALINTERVALLE (SOLAR+THERMISCH).
!
!
!
!    SPEKTRALE INTEGRATION MIT NLA LSF UND
!                     EXTINKTION = ABSORPTION + STREUUNG
!    BEI AUFBEREITUNG DER KOEFFIZIENTEN, DANN EINFACHSTREUALBEDO
!
  REAL  (wp)  ::  &
    zlwe(4,jpspec), &  ! 
    zlww(2,jpspec), &  !
    zlwg(2,jpspec), &  !
    zlwemn(jpspec), &  ! minimum values of the extinction coefficients in Pa(H2O)**-1
    zlwemx(jpspec)     ! maximum values of the extinction coefficients in Pa(H2O)**-1
  
      DATA zlwe /                                                     &
     &      -23.014052_wp, .173026_wp, .811865_wp, .000453_wp, &
     &      -28.122596_wp, .172211_wp, .705673_wp, .000457_wp, &
     &      -28.162592_wp, .198665_wp, .810637_wp, .000550_wp, &
     &     -170.687770_wp, .498371_wp, .356225_wp, .001330_wp, &
     &      -68.573703_wp, .263182_wp, .568143_wp, .000776_wp, &
     &     -122.833213_wp, .297599_wp, .306486_wp, .000976_wp, &
     &     -192.594948_wp, .440659_wp, .317142_wp, .001027_wp, &
     &      -62.018469_wp, .281406_wp, .732715_wp, .000611_wp/
      DATA zlww /                           &
     &               .989679_wp,  -22.291412_wp,  &
     &               .999529_wp,     .020875_wp,  &
     &               .999999_wp,     .000000_wp,  &
     &               .302657_wp,  102.711916_wp,  &
     &               .337398_wp,   80.596716_wp,  &
     &               .449308_wp,   52.823880_wp,  &
     &               .686930_wp,  -29.876242_wp,  &
     &               .804203_wp, -103.022685_wp/
      DATA zlwg /                           &
     &               .804992_wp,   17.901033_wp,  &
     &               .814785_wp,   14.204375_wp,  &
     &               .843955_wp,    8.306586_wp,  &
     &               .279400_wp,  124.179115_wp,  &
     &               .499491_wp,  131.635343_wp,  &
     &               .696708_wp,   75.061613_wp,  &
     &               .704732_wp,   77.778408_wp,  &
     &               .784672_wp,   38.002913_wp/
!
!     MINDEST- UND MAXIMALWERTE DER EXTINKTIONSKOEFFIZIENTEN IN DER
!                    EINHEIT PA(H2O)**-1
!
      DATA zlwemn/ 5.000_wp, 5.000_wp, 4.930_wp, 5.800_wp, &
     &             5.400_wp, 5.200_wp, 5.500_wp, 5.500_wp/
      DATA zlwemx/32.500_wp,32.500_wp,31.360_wp,18.600_wp, &
     &            24.500_wp,18.200_wp,20.200_wp,32.400_wp/


!     Aerosol optical properties for 8 spectral intervals
!
!     *ZAEA(ISPEC,IAER)* : VERHAELTNIS DER OPTISCHEN DICKE FUER ABSORP-
!                          TION IM INTERVALL *ISPEC* ZUR GESAMTEN OP-
!                          TISCHEN DICKE BEI 0.55M*1.E-06 FUER DEN
!                          AEROSOLTYP *IAER* berechnet nach OPAC
!
!     *ZAES(ISPEC,IAER)* : ANALOG FUER DIE OPTISCHE DICKE FUER STREUUNG
!
!     *ZAEG(ISPEC,IAER)* : ASYMETRIEFAKTOR FUER DEN JEWEILIGEN AEROSOL-
!                          TYP IM BETRACHTETEN INTERVALL
!     IAER     AEROSOLTYP 
!        1     CONTINENTAL (GME,EZMW, Tanre--> GADS)
!        2     MARITIME clean (           "        )
!        4     DESERT         (           "        )
!        3     URBAN          (           "        )
!        5     STRATOSPHAERISCHES 'BACKGROUND'-AEROSOL (GME)
! Initialization of zaea, zaes, zaef moved to mo_aerosol_util.

!
!
!     Optische Eigenschaften fuer Eiswolken fuer alle Spektralintervalle        
!
!      Datei:  /f/for2rit/pcrad/datiw   
!
!      Erstellungsmethode:  Spektrale Mittelung mit gewichtetem, nicht
!                           linearem LSF; Rekombination von Extinktion,
!                           Streuung und Absorption nach Mittelung durch
!                           Extinktion = Streuung + Absorption
!

      REAL(wp) :: ziwe(4,jpspec),ziweMN(jpspec),ziwemx(jpspec),  &
     &                     ziww(2,jpspec),ziwg(2,jpspec)
      
      DATA ziwe /                                                     &
     &    16.726535_wp, 0.007465_wp, 1.354626_wp, 0.000112_wp, &
     &    17.531261_wp, 0.003949_wp, 0.669605_wp, 0.000058_wp, &
     &    17.698999_wp, 0.003657_wp, 0.625067_wp, 0.000055_wp, &
     &    19.592746_wp, 0.008644_wp, 1.153213_wp, 0.000101_wp, &
     &    18.990998_wp, 0.006743_wp, 0.997361_wp, 0.000080_wp, &
     &    18.482156_wp, 0.004408_wp, 0.693883_wp, 0.000060_wp, &
     &    18.603168_wp, 0.005260_wp, 0.813026_wp, 0.000064_wp, &
     &    18.437818_wp, 0.004378_wp, 0.692678_wp, 0.000057_wp/
!
!     Koeffizienten fuer logarithmischen Fit
!
      DATA ziww /                           &
     &              0.694631_wp,   -0.022160_wp,  &
     &              0.998669_wp,   -0.000107_wp,  &
     &              0.999993_wp,    0.000000_wp,  &
     &              0.289966_wp,   -0.033855_wp,  &
     &              0.555820_wp,    0.004491_wp,  &
     &              0.554495_wp,    0.004904_wp,  &
     &              0.375319_wp,   -0.017168_wp,  &
     &              0.485290_wp,   -0.004358_wp/
      DATA ziwg /                           &
     &              0.976960_wp,    0.007938_wp,  &
     &              0.914842_wp,    0.003334_wp,  &
     &              0.900536_wp,    0.001797_wp,  &
     &              1.134025_wp,    0.032141_wp,  &
     &              1.053136_wp,    0.015721_wp,  &
     &              1.010632_wp,    0.006844_wp,  &
     &              1.063545_wp,    0.011772_wp,  &
     &              1.035725_wp,    0.008755_wp/
!     
!     Koeffizienten fuer linearen Fit
!
!     DATA ZIWW /                           &
!    &              0.910929,  -82.979127,  &
!    &              0.999711,   -0.388518,  &
!    &              0.999997,   -0.005074,  &
!    &              0.612437, -105.929725,  &
!    &              0.518593,    5.803072,  &
!    &              0.505939,   19.703743,  &
!    &              0.549304,  -68.644080,  &
!    &              0.534472,  -25.724857/
!     DATA ZIWG /                           &
!    &              0.897605,   32.198971,  &
!    &              0.881326,   13.716064,  &
!    &              0.882406,    7.478366,  &
!    &              0.816268,  126.067881,  &
!    &              0.894497,   63.721895,  &
!    &              0.939750,   31.333999,  &
!    &              0.947264,   44.064562,  &
!    &              0.944844,   40.010780/
!
!     MINDEST- UND MAXIMALWERTE DER EXTINKTIONSKOEFFIZIENTEN IN DER
!                    EINHEIT PA(H2O)**-1
!
      DATA ziwemn/ 2.000_wp, 2.000_wp, 2.000_wp, 2.000_wp, &
     &             2.000_wp, 2.000_wp, 2.000_wp, 2.000_wp/
      DATA ziwemx/30.000_wp,30.000_wp,30.000_wp,30.000_wp, &
     &            30.000_wp,30.000_wp,30.000_wp,30.000_wp/
  
CONTAINS


  SUBROUTINE rad_aibi
!
!      Purpose
!      -------
!
!      *rad_aibi* provides the data statements for the calculation
!                 of the absorption properties of atmospheric gases

!      Interface
!      ---------
!
!      *rad_aibi* has to be called once at the beginning of a model run
!                 in order to convert the raw coffeicients from the data
!                 variables to those entities used in the radiation code
!                 through a corresponding common block
!      Method
!      ------
!
!      - the limits of the spectral intervals are provided (for information
!        only) through data statements
!      - the raw coefficients for each gas are also provided via data
!        statements
!      - skaling of absorption coefficients from the individual reference
!        temperature and pressure to unified conditions, i.e.
!        reference temperature = 273.15 K
!        reference pressure    = 1013.25 hPa
!
!=======================================================================
!
! Current Code Owner: DWD, B. Ritter
!    phone: +49-69-8062-2703, fax: +49-69-8062-3721
!    email: bodo.ritter@dwd.de
!
!=======================================================================
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! 1.1        2000/01/28 B. Ritter
!  Initial Release
! 1.14       2002/01/16 Helmut P. Frank
!  Introduce KIND-notation for all reals with real2kind.pl
! V2_24        2010/04/28 Michael Gertz
!  Adaptions to SVN   
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!==============================================================================

!      INCLUDE "gme_par_rad.h"
!      INCLUDE "gme_commpi.h"
!C
      INTEGER   i,jc,jg,js,jpabsk,jpgas
      PARAMETER ( jpabsk = 7, jpgas = 3)
!C
!C--------------------------------------------------------------------
!C     ACHTUNG.
!C     -------
!C     DIE MAXIMAL ZULAESSIGE ZAHL VON ABSORPTIONSKOEFFIZIENTEN
!C     *JPABSK* MUSS IM COMMON-BLOCK *COMGAS*  BERUECKSICHTIGT WERDEN.
!C--------------------------------------------------------------------
!C
       INTEGER icgas(jpspec,jpgas)
       REAL(wp) :: zai (jpabsk,jpspec,jpgas), zbi (jpabsk,jpspec,jpgas)
       REAL(wp) :: zali(jpabsk,jpspec,jpgas), zbti(jpabsk,jpspec,jpgas)
       REAL(wp) :: pgas       (jpspec,jpgas), tgas       (jpspec,jpgas)
!C
!      INCLUDE "gme_com_gas.h"
!      INCLUDE "gme_com_grz.h"
!      INCLUDE "gme_dat_grz.h"
!      INCLUDE "gme_dat_gas.h"

!
!C*   STEUERPARAMETER ZUR AUSWAHL VON ESFT (0) ODER FESFT (1)
!C*   IN DEN EINZELNEN SPEKTRALINTERVALLEN.
!C
!=======================================================================
!
! Current Code Owner: DWD, B. Ritter
!    phone: +49-69-8062-2703, fax: +49-69-8062-3721
!    email: bodo.ritter@dwd.de
!
!=======================================================================
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! 1.1        2000/01/28 B. Ritter
!  Initial Release
! 1.14       2002/01/16 Helmut P. Frank
!  Introduce KIND-notation for all reals with real2kind.pl
! 1.15       2002/02/15 A. Mueller  
!  No more DATA initialisation  but simple setting of NFAST
! V2_24        2010/04/28 Michael Gertz
!  Adaptions to SVN     
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!
!=======================================================================
!
!C    DATA NFAST / 0, 0, 0, 0, 0, 0, 0, 0/
!C
!C                                                                       
!C*** KOEFFIZIENTEN FUER H2O          INTERVALL: 1                      
!C    REFERENZDRUCK     : 101325.00 PA                                  
!C    REFERENZTEMPERATUR:    281.70 K                                   
      DATA ICGAS(1,1)/ 7/                                               
      DATA PGAS(1,1)/101325.000_wp/                                        
      DATA TGAS(1,1)/   281.700_wp/                                        
      DATA (ZAI(I,1,1),I=1,7)/                                         &
     &    .114190E-01_wp,  .600200E-01_wp,  .111201E+00_wp,  .123340E+00_wp,   &
     &    .902500E-01_wp,  .199632E+00_wp,  .404139E+00_wp/                    
      DATA (ZBI(I,1,1),I=1,7)/                                         &
     &    .209894E+02_wp,  .208930E+01_wp,  .184502E+00_wp,  .217771E-01_wp,   &
     &    .279254E-02_wp,  .463447E-03_wp, 0.000000E+00_wp/                    
      DATA (ZALI(I,1,1),I=1,7)/                                        &
     &    .285370E-01_wp,  .688620E+00_wp,  .766031E+00_wp,  .833136E+00_wp,   &
     &    .773491E+00_wp,  .768818E+00_wp,  .100000E+01_wp/                    
      DATA (ZBTI(I,1,1),I=1,7)/                                        &
     &    .473006E+00_wp, -.468639E+00_wp, -.599601E+00_wp, -.162223E+01_wp,   &
     &   -.176002E+01_wp, -.153131E+01_wp,  .100000E+01_wp/                    
!C                                                                       
!C*** KOEFFIZIENTEN FUER H2O          INTERVALL: 2                      
!C    REFERENZDRUCK     : 101325.00 PA                                  
!C    REFERENZTEMPERATUR:    281.70 K                                   
      DATA ICGAS(2,1)/ 7/                                               
      DATA PGAS(2,1)/101325.000_wp/                                        
      DATA TGAS(2,1)/   281.700_wp/                                        
      DATA (ZAI(I,2,1),I=1,7)/                                         &
     &    .201500E-02_wp,  .268530E-01_wp,  .598920E-01_wp,  .907740E-01_wp,   &
     &    .102284E+00_wp,  .217298E+00_wp,  .500884E+00_wp/                    
      DATA (ZBI(I,2,1),I=1,7)/                                         &
     &    .508159E+01_wp,  .519996E+00_wp,  .465586E-01_wp,  .891251E-02_wp,   &
     &    .159221E-02_wp,  .374973E-03_wp, 0.000000E+00_wp/                    
      DATA (ZALI(I,2,1),I=1,7)/                                        &
     &   -.482300E-02_wp,  .529161E+00_wp,  .587751E+00_wp,  .756567E+00_wp,   &
     &    .774607E+00_wp,  .733883E+00_wp,  .100000E+01_wp/                    
      DATA (ZBTI(I,2,1),I=1,7)/                                        &
     &    .499755E+00_wp, -.529716E+00_wp, -.177970E-01_wp, -.746447E+00_wp,   &
     &   -.106191E+00_wp, -.727589E+00_wp,  .100000E+01_wp/                    
!C                                                                      
!C*** KOEFFIZIENTEN FUER H2O          INTERVALL: 3                      
!C    REFERENZDRUCK     : 101325.00 PA                                  
!C    REFERENZTEMPERATUR:    281.70 K                                   
      DATA ICGAS(3,1)/ 3/                                               
      DATA PGAS(3,1)/101325.000_wp/                                        
      DATA TGAS(3,1)/   281.700_wp/                                        
      DATA (ZAI(I,3,1),I=1,7)/                                         &
     &    .566900E-02_wp,  .346720E-01_wp,  .959659E+00_wp, 0.000000E+00_wp,   &
     &   0.000000E+00_wp, 0.000000E+00_wp, 0.000000E+00_wp/                    
      DATA (ZBI(I,3,1),I=1,7)/                                         &
     &    .716144E-03_wp,  .256449E-03_wp, 0.000000E+00_wp, 0.000000E+00_wp,   &
     &   0.000000E+00_wp, 0.000000E+00_wp, 0.000000E+00_wp/                    
      DATA (ZALI(I,3,1),I=1,7)/                                        &
     &   -.281669E+00_wp,  .611673E+00_wp,  .100000E+01_wp, 0.000000E+00_wp,   &
     &   0.000000E+00_wp, 0.000000E+00_wp, 0.000000E+00_wp/                    
      DATA (ZBTI(I,3,1),I=1,7)/                                        &
     &    .418657E+00_wp,  .405230E-01_wp,  .100000E+01_wp, 0.000000E+00_wp,   &
     &   0.000000E+00_wp, 0.000000E+00_wp, 0.000000E+00_wp/                    
!C                                                                      
!C*** KOEFFIZIENTEN FUER H2O          INTERVALL: 4                      
!C    REFERENZDRUCK     :  50662.50 PA                                  
!C    REFERENZTEMPERATUR:    255.80 K                                   
      DATA ICGAS(4,1)/ 7/                                               
      DATA PGAS(4,1)/ 50662.500_wp/                                        
      DATA TGAS(4,1)/   255.800_wp/                                        
      DATA (ZAI(I,4,1),I=1,7)/                                         &
     &    .641200E-02_wp,  .362630E-01_wp,  .147064E+00_wp,  .285387E+00_wp,   &
     &    .246376E+00_wp,  .226899E+00_wp,  .515980E-01_wp/                    
      DATA (ZBI(I,4,1),I=1,7)/                                         &
     &    .298538E+04_wp,  .139959E+03_wp,  .152405E+02_wp,  .144212E+01_wp,   &
     &    .183654E+00_wp,  .283139E-01_wp,  .409261E-02_wp/                    
      DATA (ZALI(I,4,1),I=1,7)/                                        &
     &    .183780E-01_wp,  .410557E+00_wp,  .808897E+00_wp,  .897332E+00_wp,   &
     &    .932149E+00_wp,  .978389E+00_wp,  .100000E+01_wp/                    
      DATA (ZBTI(I,4,1),I=1,7)/                                        &
     &    .413777E+00_wp, -.663704E+00_wp, -.953789E+00_wp, -.111883E+01_wp,   &
     &   -.156269E+01_wp, -.330557E+01_wp,  .100000E+01_wp/                    
!C                                                                      
!C*** KOEFFIZIENTEN FUER H2O          INTERVALL: 5                      
!C    REFERENZDRUCK     :  86126.25 PA                                  
!C    REFERENZTEMPERATUR:    281.70 K                                   
      DATA ICGAS(5,1)/ 7/                                               
      DATA PGAS(5,1)/ 86126.250_wp/                                        
      DATA TGAS(5,1)/   281.700_wp/                                        
      DATA (ZAI(I,5,1),I=1,7)/                                         &
     &    .147700E-02_wp,  .345020E-01_wp,  .865590E-01_wp,  .144237E+00_wp,   &
     &    .218089E+00_wp,  .339440E+00_wp,  .175697E+00_wp/                    
      DATA (ZBI(I,5,1),I=1,7)/                                         &
     &    .126765E+02_wp,  .149624E+01_wp,  .147571E+00_wp,  .368129E-01_wp,   &
     &    .792501E-02_wp,  .208930E-02_wp, 0.000000E+00_wp/                    
      DATA (ZALI(I,5,1),I=1,7)/                                        &
     &   -.414300E-02_wp,  .504464E+00_wp,  .670985E+00_wp,  .920940E+00_wp,   &
     &    .889089E+00_wp,  .966028E+00_wp,  .100000E+01_wp/                    
      DATA (ZBTI(I,5,1),I=1,7)/                                        &
     &    .454691E+00_wp, -.423980E+01_wp, -.340869E+01_wp, -.410896E+01_wp,   &
     &   -.268068E+01_wp, -.250967E+01_wp,  .100000E+01_wp/                    
!C                                                                      
!C*** KOEFFIZIENTEN FUER H2O          INTERVALL: 6                      
!C    REFERENZDRUCK     :  86126.25 PA                                  
!C    REFERENZTEMPERATUR:    281.70 K                                   
      DATA ICGAS(6,1)/ 4/                                               
      DATA PGAS(6,1)/ 86126.250_wp/                                        
      DATA TGAS(6,1)/   281.700_wp/                                        
      DATA (ZAI(I,6,1),I=1,7)/                                         &
     &    .653200E-02_wp,  .700040E-01_wp,  .243768E+00_wp,  .679696E+00_wp,   &
     &   0.000000E+00_wp, 0.000000E+00_wp, 0.000000E+00_wp/                    
      DATA (ZBI(I,6,1),I=1,7)/                                         &
     &    .632412E+00_wp,  .473151E-02_wp,  .163305E-02_wp, 0.000000E+00_wp,   &
     &   0.000000E+00_wp, 0.000000E+00_wp, 0.000000E+00_wp/                    
      DATA (ZALI(I,6,1),I=1,7)/                                        &
     &    .794801E+00_wp,  .306898E+00_wp,  .100000E+01_wp,  .100000E+01_wp,   &
     &   0.000000E+00_wp, 0.000000E+00_wp, 0.000000E+00_wp/                    
      DATA (ZBTI(I,6,1),I=1,7)/                                        &
     &   -.100000E+02_wp, -.219711E+01_wp, -.369325E+01_wp,  .100000E+01_wp,   &
     &   0.000000E+00_wp, 0.000000E+00_wp, 0.000000E+00_wp/                    
!C                                                                      
!C*** KOEFFIZIENTEN FUER H2O          INTERVALL: 7                      
!C    REFERENZDRUCK     :  86126.25 PA                                  
!C    REFERENZTEMPERATUR:    281.70 K                                   
      DATA ICGAS(7,1)/ 3/                                               
      DATA PGAS(7,1)/ 86126.250_wp/                                        
      DATA TGAS(7,1)/   281.700_wp/                                        
      DATA (ZAI(I,7,1),I=1,7)/                                         &
     &    .138610E-01_wp,  .226595E+00_wp,  .759544E+00_wp, 0.000000E+00_wp,   &
     &   0.000000E+00_wp, 0.000000E+00_wp, 0.000000E+00_wp/                    
      DATA (ZBI(I,7,1),I=1,7)/                                         &
     &    .425598E-02_wp,  .155239E-02_wp, 0.000000E+00_wp, 0.000000E+00_wp,   &
     &   0.000000E+00_wp, 0.000000E+00_wp, 0.000000E+00_wp/                    
      DATA (ZALI(I,7,1),I=1,7)/                                        &
     &   -.736171E+00_wp,  .805828E+00_wp,  .100000E+01_wp, 0.000000E+00_wp,   &
     &   0.000000E+00_wp, 0.000000E+00_wp, 0.000000E+00_wp/                    
      DATA (ZBTI(I,7,1),I=1,7)/                                        &
     &    .308301E+00_wp, -.267573E+01_wp,  .100000E+01_wp, 0.000000E+00_wp,   &
     &   0.000000E+00_wp, 0.000000E+00_wp, 0.000000E+00_wp/                    
!C                                                                      
!C*** KOEFFIZIENTEN FUER H2O          INTERVALL: 8                      
!C    REFERENZDRUCK     :  75993.75 PA                                  
!C    REFERENZTEMPERATUR:    281.70 K                                   
      DATA ICGAS(8,1)/ 7/                                               
      DATA PGAS(8,1)/ 75993.750_wp/                                        
      DATA TGAS(8,1)/   281.700_wp/                                        
      DATA (ZAI(I,8,1),I=1,7)/                                         &
     &    .181840E-01_wp,  .106586E+00_wp,  .237611E+00_wp,  .241085E+00_wp,   &
     &    .157304E+00_wp,  .178767E+00_wp,  .604640E-01_wp/                    
      DATA (ZBI(I,8,1),I=1,7)/                                         &
     &    .822243E+02_wp,  .979490E+01_wp,  .905733E+00_wp,  .140281E+00_wp,   &
     &    .193197E-01_wp,  .320627E-02_wp, 0.000000E+00_wp/                    
      DATA (ZALI(I,8,1),I=1,7)/                                        &
     &   -.126888E+00_wp,  .701873E+00_wp,  .834941E+00_wp,  .920550E+00_wp,   &
     &    .849506E+00_wp,  .931957E+00_wp,  .100000E+01_wp/                    
      DATA (ZBTI(I,8,1),I=1,7)/                                        &
     &    .384580E+00_wp, -.187972E+01_wp, -.226834E+01_wp, -.475940E+01_wp,   &
     &   -.589531E+01_wp, -.395962E+01_wp,  .100000E+01_wp/                    
!C                                                                      
!C*** KOEFFIZIENTEN FUER CO2 + ...    INTERVALL: 1                      
!C    REFERENZDRUCK     :  86126.25 PA                                  
!C    REFERENZTEMPERATUR:    255.80 K                                   
      DATA ICGAS(1,2)/ 6/                                               
      DATA PGAS(1,2)/ 86126.250_wp/                                        
      DATA TGAS(1,2)/   255.800_wp/                                        
      DATA (ZAI(I,1,2),I=1,7)/                                         &
     &    .592000E-02_wp,  .667700E-02_wp,  .423020E-01_wp,  .732310E-01_wp,   &
     &    .140143E+00_wp,  .731727E+00_wp, 0.000000E+00_wp/                    
      DATA (ZBI(I,1,2),I=1,7)/                                         &
     &    .760326E+02_wp,  .480839E+01_wp,  .391742E+00_wp,  .133968E-01_wp,   &
     &    .355631E-02_wp, 0.000000E+00_wp, 0.000000E+00_wp/                    
      DATA (ZALI(I,1,2),I=1,7)/                                        &
     &    .659071E+00_wp,  .240858E+00_wp,  .694157E+00_wp,  .424843E+00_wp,   &
     &    .694262E+00_wp,  .100000E+01_wp, 0.000000E+00_wp/                    
      DATA (ZBTI(I,1,2),I=1,7)/                                        &
     &    .467048E+00_wp,  .395422E+00_wp, -.902210E+00_wp, -.557526E+00_wp,   &
     &   -.473196E+00_wp,  .100000E+01_wp, 0.000000E+00_wp/                    
!C                                                                      
!C*** KOEFFIZIENTEN FUER CO2 + ...    INTERVALL: 2                      
!C    REFERENZDRUCK     :  86126.25 PA                                  
!C    REFERENZTEMPERATUR:    255.80 K                                   
      DATA ICGAS(2,2)/ 3/                                               
      DATA PGAS(2,2)/ 86126.250_wp/                                        
      DATA TGAS(2,2)/   255.800_wp/                                        
      DATA (ZAI(I,2,2),I=1,7)/                                         &
     &    .278000E-02_wp,  .197330E-01_wp,  .977487E+00_wp, 0.000000E+00_wp,   &
     &   0.000000E+00_wp, 0.000000E+00_wp, 0.000000E+00_wp/                    
      DATA (ZBI(I,2,2),I=1,7)/                                         &
     &    .169434E+00_wp,  .103753E-01_wp, 0.000000E+00_wp, 0.000000E+00_wp,   &
     &   0.000000E+00_wp, 0.000000E+00_wp, 0.000000E+00_wp/                    
      DATA (ZALI(I,2,2),I=1,7)/                                        &
     &    .138563E+00_wp,  .831359E+00_wp,  .100000E+01_wp, 0.000000E+00_wp,   &
     &   0.000000E+00_wp, 0.000000E+00_wp, 0.000000E+00_wp/                    
      DATA (ZBTI(I,2,2),I=1,7)/                                        &
     &    .475293E+00_wp, -.496213E+00_wp,  .100000E+01_wp, 0.000000E+00_wp,   &
     &   0.000000E+00_wp, 0.000000E+00_wp, 0.000000E+00_wp/                    
!C                                                                      
!C*** KOEFFIZIENTEN FUER CO2 + ...    INTERVALL: 3                      
!C    REFERENZDRUCK     :  86126.25 PA                                  
!C    REFERENZTEMPERATUR:    255.80 K                                   
      DATA ICGAS(3,2)/ 2/                                               
      DATA PGAS(3,2)/ 86126.250_wp/                                        
      DATA TGAS(3,2)/   255.800_wp/                                        
      DATA (ZAI(I,3,2),I=1,7)/                                         &
     &    .306100E-02_wp,  .996939E+00_wp, 0.000000E+00_wp, 0.000000E+00_wp,   &
     &   0.000000E+00_wp, 0.000000E+00_wp, 0.000000E+00_wp/                    
      DATA (ZBI(I,3,2),I=1,7)/                                         &
     &    .101625E+00_wp, 0.000000E+00_wp, 0.000000E+00_wp, 0.000000E+00_wp,   &
     &   0.000000E+00_wp, 0.000000E+00_wp, 0.000000E+00_wp/                    
      DATA (ZALI(I,3,2),I=1,7)/                                        &
     &    .100000E+01_wp,  .100000E+01_wp, 0.000000E+00_wp, 0.000000E+00_wp,   &
     &   0.000000E+00_wp, 0.000000E+00_wp, 0.000000E+00_wp/                    
      DATA (ZBTI(I,3,2),I=1,7)/                                        &
     &   -.100670E+01_wp,  .100000E+01_wp, 0.000000E+00_wp, 0.000000E+00_wp,   &
     &   0.000000E+00_wp, 0.000000E+00_wp, 0.000000E+00_wp/                    
!C                                                                      
!C*** KOEFFIZIENTEN FUER CO2 + ...    INTERVALL: 4                      
!C    REFERENZDRUCK     :  60795.00 PA                                  
!C    REFERENZTEMPERATUR:    255.80 K                                   
      DATA ICGAS(4,2)/ 0/                                               
      DATA PGAS(4,2)/ 60795.000_wp/                                        
      DATA TGAS(4,2)/   255.800_wp/                                        
      DATA (ZAI(I,4,2),I=1,7)/                                         &
     &    .335800E-02_wp,  .996642E+00_wp, 0.000000E+00_wp, 0.000000E+00_wp,   &
     &   0.000000E+00_wp, 0.000000E+00_wp, 0.000000E+00_wp/                    
      DATA (ZBI(I,4,2),I=1,7)/                                         &
     &    .247172E+00_wp, 0.000000E+00_wp, 0.000000E+00_wp, 0.000000E+00_wp,   &
     &   0.000000E+00_wp, 0.000000E+00_wp, 0.000000E+00_wp/                    
      DATA (ZALI(I,4,2),I=1,7)/                                        &
     &    .100000E+01_wp,  .100000E+01_wp, 0.000000E+00_wp, 0.000000E+00_wp,   &
     &   0.000000E+00_wp, 0.000000E+00_wp, 0.000000E+00_wp/                    
      DATA (ZBTI(I,4,2),I=1,7)/                                        &
     &   -.807310E+00_wp,  .100000E+01_wp, 0.000000E+00_wp, 0.000000E+00_wp,   &
     &   0.000000E+00_wp, 0.000000E+00_wp, 0.000000E+00_wp/                    
!C                                                                      
!C*** KOEFFIZIENTEN FUER CO2 + ...    INTERVALL: 5                      
!C    REFERENZDRUCK     :  10132.50 PA                                  
!C    REFERENZTEMPERATUR:    229.90 K                                   
      DATA ICGAS(5,2)/ 7/                                               
      DATA PGAS(5,2)/ 10132.500_wp/                                        
      DATA TGAS(5,2)/   229.900_wp/                                        
      DATA (ZAI(I,5,2),I=1,7)/                                         &
     &    .452500E-02_wp,  .321420E-01_wp,  .659180E-01_wp,  .101074E+00_wp,   &
     &    .107224E+00_wp,  .186663E+00_wp,  .502454E+00_wp/                    
      DATA (ZBI(I,5,2),I=1,7)/                                         &
     &    .299226E+03_wp,  .364754E+02_wp,  .271644E+01_wp,  .570164E+00_wp,   &
     &    .100231E+00_wp,  .224388E-01_wp, 0.000000E+00_wp/                    
      DATA (ZALI(I,5,2),I=1,7)/                                        &
     &    .466819E+00_wp,  .319510E+00_wp,  .596734E+00_wp,  .751216E+00_wp,   &
     &    .708519E+00_wp,  .744381E+00_wp,  .100000E+01_wp/                    
      DATA (ZBTI(I,5,2),I=1,7)/                                        &
     &    .358348E+00_wp, -.739332E+00_wp, -.183599E+01_wp, -.289470E+01_wp,   &
     &   -.214575E+01_wp, -.585028E+01_wp,  .100000E+01_wp/                    
!C                                                                      
!C*** KOEFFIZIENTEN FUER CO2 + ...    INTERVALL: 6                      
!C    REFERENZDRUCK     :  50662.50 PA                                  
!C    REFERENZTEMPERATUR:    255.80 K                                   
      DATA ICGAS(6,2)/ 3/                                               
      DATA PGAS(6,2)/ 50662.500_wp/                                        
      DATA TGAS(6,2)/   255.800_wp/                                        
      DATA (ZAI(I,6,2),I=1,7)/                                         &
     &    .119551E+00_wp,  .899140E-01_wp,  .790535E+00_wp, 0.000000E+00_wp,   &
     &   0.000000E+00_wp, 0.000000E+00_wp, 0.000000E+00_wp/                    
      DATA (ZBI(I,6,2),I=1,7)/                                         &
     &    .305492E-02_wp,  .148936E-02_wp, 0.000000E+00_wp, 0.000000E+00_wp,   &
     &   0.000000E+00_wp, 0.000000E+00_wp, 0.000000E+00_wp/                    
      DATA (ZALI(I,6,2),I=1,7)/                                        &
     &    .783365E+00_wp, -.113116E+00_wp,  .100000E+01_wp, 0.000000E+00_wp,   &
     &   0.000000E+00_wp, 0.000000E+00_wp, 0.000000E+00_wp/                    
      DATA (ZBTI(I,6,2),I=1,7)/                                        &
     &   -.447333E+01_wp,  .296352E+01_wp,  .100000E+01_wp, 0.000000E+00_wp,   &
     &   0.000000E+00_wp, 0.000000E+00_wp, 0.000000E+00_wp/                    
!C                                                                      
!C*** KOEFFIZIENTEN FUER CO2 + ...    INTERVALL: 7                      
!C    REFERENZDRUCK     :  50662.50 PA                                  
!C    REFERENZTEMPERATUR:    255.80 K                                   
      DATA ICGAS(7,2)/ 3/                                               
      DATA PGAS(7,2)/ 50662.500_wp/                                        
      DATA TGAS(7,2)/   255.800_wp/                                        
      DATA (ZAI(I,7,2),I=1,7)/                                         &
     &    .577890E-01_wp,  .321750E-01_wp,  .910036E+00_wp, 0.000000E+00_wp,   &
     &   0.000000E+00_wp, 0.000000E+00_wp, 0.000000E+00_wp/                    
      DATA (ZBI(I,7,2),I=1,7)/                                         &
     &    .650130E-02_wp,  .309030E-02_wp, 0.000000E+00_wp, 0.000000E+00_wp,   &
     &   0.000000E+00_wp, 0.000000E+00_wp, 0.000000E+00_wp/                    
      DATA (ZALI(I,7,2),I=1,7)/                                        &
     &    .295465E+00_wp,  .930860E-01_wp,  .100000E+01_wp, 0.000000E+00_wp,   &
     &   0.000000E+00_wp, 0.000000E+00_wp, 0.000000E+00_wp/                    
      DATA (ZBTI(I,7,2),I=1,7)/                                        &
     &   -.562957E+01_wp, -.984577E+01_wp,  .100000E+01_wp, 0.000000E+00_wp,   &
     &   0.000000E+00_wp, 0.000000E+00_wp, 0.000000E+00_wp/                    
!C                                                                       
!C*** KOEFFIZIENTEN FUER CO2 + ...    INTERVALL: 8                      
!C    REFERENZDRUCK     :  50662.50 PA                                  
!C    REFERENZTEMPERATUR:    255.80 K                                   
      DATA ICGAS(8,2)/ 4/                                               
      DATA PGAS(8,2)/ 50662.500_wp/                                        
      DATA TGAS(8,2)/   255.800_wp/                                        
      DATA (ZAI(I,8,2),I=1,7)/                                         &
     &    .317000E-02_wp,  .127109E+00_wp,  .114118E+00_wp,  .755604E+00_wp,   &
     &   0.000000E+00_wp, 0.000000E+00_wp, 0.000000E+00_wp/                    
      DATA (ZBI(I,8,2),I=1,7)/                                         &
     &    .174181E+02_wp,  .495450E-01_wp,  .165196E-01_wp, 0.000000E+00_wp,   &
     &   0.000000E+00_wp, 0.000000E+00_wp, 0.000000E+00_wp/                    
      DATA (ZALI(I,8,2),I=1,7)/                                        &
     &    .511300E-02_wp,  .252848E+00_wp,  .851104E+00_wp,  .100000E+01_wp,   &
     &   0.000000E+00_wp, 0.000000E+00_wp, 0.000000E+00_wp/                    
      DATA (ZBTI(I,8,2),I=1,7)/                                        &
     &    .495222E+00_wp,  .445084E+00_wp,  .117957E+01_wp,  .100000E+01_wp,   &
     &   0.000000E+00_wp, 0.000000E+00_wp, 0.000000E+00_wp/                    
!C                                                                      
!C*** KOEFFIZIENTEN FUER O3           INTERVALL: 1                      
!C    REFERENZDRUCK     :3039.75PA                                      
!C    REFERENZTEMPERATUR:    229.90 K                                   
      DATA ICGAS(1,3)/ 0/                                               
      DATA PGAS(1,3)/3039.75_wp/                                           
      DATA TGAS(1,3)/   229.900_wp/                                        
      DATA (ZAI(I,1,3),I=1,7)/                                         &
     &    .306000E-03_wp,  .999694E+00_wp, 0.000000E+00_wp, 0.000000E+00_wp,   &
     &   0.000000E+00_wp, 0.000000E+00_wp, 0.000000E+00_wp/                    
      DATA (ZBI(I,1,3),I=1,7)/                                         &
     &    .409261E+02_wp, 0.000000E+00_wp, 0.000000E+00_wp, 0.000000E+00_wp,   &
     &   0.000000E+00_wp, 0.000000E+00_wp, 0.000000E+00_wp/                    
      DATA (ZALI(I,1,3),I=1,7)/                                        &
     &    .618332E+00_wp,  .100000E+01_wp, 0.000000E+00_wp, 0.000000E+00_wp,   &
     &   0.000000E+00_wp, 0.000000E+00_wp, 0.000000E+00_wp/                    
      DATA (ZBTI(I,1,3),I=1,7)/                                        &
     &   -.974847E+00_wp,  .100000E+01_wp, 0.000000E+00_wp, 0.000000E+00_wp,   &
     &   0.000000E+00_wp, 0.000000E+00_wp, 0.000000E+00_wp/                    
!C                                                                       
!C*** KOEFFIZIENTEN FUER O3           INTERVALL: 2                      
!C    REFERENZDRUCK     :3039.75PA                                      
!C    REFERENZTEMPERATUR:    229.90 K                                   
      DATA ICGAS(2,3)/ 0/                                               
      DATA PGAS(2,3)/3039.75_wp/                                           
      DATA TGAS(2,3)/   229.900_wp/                                        
      DATA (ZAI(I,2,3),I=1,7)/                                         &
     &    .154800E-02_wp,  .998452E+00_wp, 0.000000E+00_wp, 0.000000E+00_wp,   &
     &   0.000000E+00_wp, 0.000000E+00_wp, 0.000000E+00_wp/                    
      DATA (ZBI(I,2,3),I=1,7)/                                         &
     &    .395367E+02_wp, 0.000000E+00_wp, 0.000000E+00_wp, 0.000000E+00_wp,   &
     &   0.000000E+00_wp, 0.000000E+00_wp, 0.000000E+00_wp/                    
      DATA (ZALI(I,2,3),I=1,7)/                                        &
     &    .592629E+00_wp,  .100000E+01_wp, 0.000000E+00_wp, 0.000000E+00_wp,   &
     &   0.000000E+00_wp, 0.000000E+00_wp, 0.000000E+00_wp/                    
      DATA (ZBTI(I,2,3),I=1,7)/                                        &
     &   -.106087E+01_wp,  .100000E+01_wp, 0.000000E+00_wp, 0.000000E+00_wp,   &
     &   0.000000E+00_wp, 0.000000E+00_wp, 0.000000E+00_wp/                    
!C                                                                      
!C*** KOEFFIZIENTEN FUER O3           INTERVALL: 3                      
!C    REFERENZDRUCK     :3039.75PA                                      
!C    REFERENZTEMPERATUR:    229.90 K                                   
      DATA ICGAS(3,3)/ 5/                                               
      DATA PGAS(3,3)/3039.75_wp/                                           
      DATA TGAS(3,3)/   229.900_wp/                                        
      DATA (ZAI(I,3,3),I=1,7)/                                         &
     &    .564000E-03_wp,  .108690E-01_wp,  .124320E-01_wp,  .184417E+00_wp,   &
     &    .791718E+00_wp, 0.000000E+00_wp, 0.000000E+00_wp/                    
      DATA (ZBI(I,3,3),I=1,7)/                                         &
     &    .191426E+05_wp,  .579429E+03_wp,  .717794E+02_wp,  .187068E+01_wp,   &
     &   0.000000E+00_wp, 0.000000E+00_wp, 0.000000E+00_wp/                    
      DATA (ZALI(I,3,3),I=1,7)/                                        &
     &   -.204400E-02_wp,  .776840E-01_wp, -.229667E+00_wp,  .994500E-01_wp,   &
     &    .100000E+01_wp, 0.000000E+00_wp, 0.000000E+00_wp/                    
      DATA (ZBTI(I,3,3),I=1,7)/                                        &
     &    .499912E+00_wp,  .485463E+00_wp,  .464581E+00_wp, -.254634E+00_wp,   &
     &    .100000E+01_wp, 0.000000E+00_wp, 0.000000E+00_wp/                    
!C                                                                       
!C*** KOEFFIZIENTEN FUER O3           INTERVALL: 4                      
!C    REFERENZDRUCK     :3039.75PA                                      
!C    REFERENZTEMPERATUR:    229.90 K                                   
      DATA ICGAS(4,3)/ 0/                                               
      DATA PGAS(4,3)/3039.75_wp/                                           
      DATA TGAS(4,3)/   229.900_wp/                                        
      DATA (ZAI(I,4,3),I=1,7)/                                         &
     &    .540000E-04_wp,  .999946E+00_wp, 0.000000E+00_wp, 0.000000E+00_wp,   &
     &   0.000000E+00_wp, 0.000000E+00_wp, 0.000000E+00_wp/                    
      DATA (ZBI(I,4,3),I=1,7)/                                         &
     &    .210378E+03_wp, 0.000000E+00_wp, 0.000000E+00_wp, 0.000000E+00_wp,   &
     &   0.000000E+00_wp, 0.000000E+00_wp, 0.000000E+00_wp/                    
      DATA (ZALI(I,4,3),I=1,7)/                                        &
     &    .490324E+00_wp,  .100000E+01_wp, 0.000000E+00_wp, 0.000000E+00_wp,   &
     &   0.000000E+00_wp, 0.000000E+00_wp, 0.000000E+00_wp/                    
      DATA (ZBTI(I,4,3),I=1,7)/                                        &
     &    .500000E+00_wp,  .100000E+01_wp, 0.000000E+00_wp, 0.000000E+00_wp,   &
     &   0.000000E+00_wp, 0.000000E+00_wp, 0.000000E+00_wp/                    
!C                                                                      
!C*** KOEFFIZIENTEN FUER O3           INTERVALL: 5                      
!C    REFERENZDRUCK     :3039.75PA                                      
!C    REFERENZTEMPERATUR:    229.90 K                                   
      DATA ICGAS(5,3)/ 2/                                               
      DATA PGAS(5,3)/3039.75_wp/                                           
      DATA TGAS(5,3)/   229.900_wp/                                        
      DATA (ZAI(I,5,3),I=1,7)/                                         &
     &    .587700E-02_wp,  .994123E+00_wp, 0.000000E+00_wp, 0.000000E+00_wp,   &
     &   0.000000E+00_wp, 0.000000E+00_wp, 0.000000E+00_wp/                    
      DATA (ZBI(I,5,3),I=1,7)/                                         &
     &    .223357E+03_wp, 0.000000E+00_wp, 0.000000E+00_wp, 0.000000E+00_wp,   &
     &   0.000000E+00_wp, 0.000000E+00_wp, 0.000000E+00_wp/                    
      DATA (ZALI(I,5,3),I=1,7)/                                        &
     &    .551312E+00_wp,  .100000E+01_wp, 0.000000E+00_wp, 0.000000E+00_wp,   &
     &   0.000000E+00_wp, 0.000000E+00_wp, 0.000000E+00_wp/                    
      DATA (ZBTI(I,5,3),I=1,7)/                                        &
     &   -.140025E+01_wp,  .100000E+01_wp, 0.000000E+00_wp, 0.000000E+00_wp,   &
     &   0.000000E+00_wp, 0.000000E+00_wp, 0.000000E+00_wp/                    
!C                                                                       
!C*** KOEFFIZIENTEN FUER O3           INTERVALL: 6                      
!C    REFERENZDRUCK     :3039.75PA                                      
!C    REFERENZTEMPERATUR:    229.90 K                                   
      DATA ICGAS(6,3)/ 0/                                               
      DATA PGAS(6,3)/3039.75_wp/                                           
      DATA TGAS(6,3)/   229.900_wp/                                        
      DATA (ZAI(I,6,3),I=1,7)/                                         &
     &    .154100E-02_wp,  .998459E+00_wp, 0.000000E+00_wp, 0.000000E+00_wp,   &
     &   0.000000E+00_wp, 0.000000E+00_wp, 0.000000E+00_wp/                    
      DATA (ZBI(I,6,3),I=1,7)/                                         &
     &    .221820E+03_wp, 0.000000E+00_wp, 0.000000E+00_wp, 0.000000E+00_wp,   &
     &   0.000000E+00_wp, 0.000000E+00_wp, 0.000000E+00_wp/                    
      DATA (ZALI(I,6,3),I=1,7)/                                        &
     &    .546048E+00_wp,  .100000E+01_wp, 0.000000E+00_wp, 0.000000E+00_wp,   &
     &   0.000000E+00_wp, 0.000000E+00_wp, 0.000000E+00_wp/                    
      DATA (ZBTI(I,6,3),I=1,7)/                                        &
     &   -.273183E+01_wp,  .100000E+01_wp, 0.000000E+00_wp, 0.000000E+00_wp,   &
     &   0.000000E+00_wp, 0.000000E+00_wp, 0.000000E+00_wp/                    
!C                                                                      
!C*** KOEFFIZIENTEN FUER O3           INTERVALL: 7                      
!C    REFERENZDRUCK     :  10132.50 PA                                  
!C    REFERENZTEMPERATUR:    204.00 K                                   
      DATA ICGAS(7,3)/ 7/                                               
      DATA PGAS(7,3)/ 10132.500_wp/                                        
      DATA TGAS(7,3)/   204.000_wp/                                        
      DATA (ZAI(I,7,3),I=1,7)/                                         &
     &    .220500E-02_wp,  .523500E-02_wp,  .951500E-02_wp,  .578800E-01_wp,   &
     &    .277389E+00_wp,  .643850E-01_wp,  .583391E+00_wp/                    
      DATA (ZBI(I,7,3),I=1,7)/                                         &
     &    .434510E+03_wp,  .299916E+03_wp,  .121339E+03_wp,  .827942E+02_wp,   &
     &    .157398E+02_wp,  .615177E+01_wp, 0.000000E+00_wp/                    
      DATA (ZALI(I,7,3),I=1,7)/                                        &
     &    .224000E-03_wp,  .100500E-02_wp,  .571600E-02_wp,  .508760E-01_wp,   &
     &    .524641E+00_wp,  .896800E-01_wp,  .100000E+01_wp/                    
      DATA (ZBTI(I,7,3),I=1,7)/                                        &
     &    .320370E+01_wp,  .130031E+01_wp, -.332851E+01_wp,  .105177E+01_wp,   &
     &   -.561714E+00_wp, -.357670E+01_wp,  .100000E+01_wp/                    
!C                                                                      
!C*** KOEFFIZIENTEN FUER O3           INTERVALL: 8                      
!C    REFERENZDRUCK     :3039.75 PA                                     
!C    REFERENZTEMPERATUR:    229.90 K                                   
      DATA ICGAS(8,3)/ 0/                                               
      DATA PGAS(8,3)/3039.75_wp/                                           
      DATA TGAS(8,3)/   229.900_wp/                                        
      DATA (ZAI(I,8,3),I=1,7)/                                         &
     &    .397000E-03_wp,  .999603E+00_wp, 0.000000E+00_wp, 0.000000E+00_wp,   &
     &   0.000000E+00_wp, 0.000000E+00_wp, 0.000000E+00_wp/                    
      DATA (ZBI(I,8,3),I=1,7)/                                         &
     &    .230675E+03_wp, 0.000000E+00_wp, 0.000000E+00_wp, 0.000000E+00_wp,   &
     &   0.000000E+00_wp, 0.000000E+00_wp, 0.000000E+00_wp/                    
      DATA (ZALI(I,8,3),I=1,7)/                                        &
     &    .564371E+00_wp,  .100000E+01_wp, 0.000000E+00_wp, 0.000000E+00_wp,   &
     &   0.000000E+00_wp, 0.000000E+00_wp, 0.000000E+00_wp/                    
      DATA (ZBTI(I,8,3),I=1,7)/                                        &
     &   -.479075E+01_wp,  .100000E+01_wp, 0.000000E+00_wp, 0.000000E+00_wp,   &
     &   0.000000E+00_wp, 0.000000E+00_wp, 0.000000E+00_wp/                    
!C

      NFAST(:) = 1

!      IF (myproc == 0) THEN
!      PRINT *,'======================================================'
!      PRINT *,'=     Radiative transfer calculations employ data    ='
!      PRINT *,'= for gaseous absortion provided in routine rad_aibi ='
!      PRINT *,'======================================================'
!      ENDIF
 
      ! Copy of data to common block variables 
 
      DO  jg=1,3
        DO  js=1,jpspec
        ncgas(js,jg)    =icgas(js,jg)
          DO  jc=1,7
          COAI (jc,js,jg)  = ZAI(jc,js,jg)
!C        COBI (jc,js,jg)  = ZBI(jc,js,jg)
          COALI(jc,js,jg) = ZALI(jc,js,jg)
          COBTI(jc,js,jg) = ZBTI(jc,js,jg)
          END DO
        END DO
      END DO
 
!     Include reference pressure and temperature in *cobi*
 
      DO jg=1,3
        DO js=1,jpspec
          DO jc=1,ncgas(js,jg)
      COBI(jc,js,jg) = ZBI(jc,js,jg) * (1._wp/PGAS(js,jg))**COALI(jc,js,jg)  &
     &                               * (   TGAS(js,jg))**COBTI(jc,js,jg)
          END DO
        END DO
      END DO
 
!     SECURITY SETTINGS, IF A GAS SHALL NOT BE CONSIDERED IN AN
!                        INTERVAL WHERE THE ESFT WILL BE USED.
      DO JG=1,3
        DO JS=1,jpspec
         IF ((NFAST(JS).EQ.0).AND.(ICGAS(JS,JG).EQ.0)) THEN
         NCGAS     (JS,JG) = 1
         COAI    (1,JS,JG) = 1.00_wp
         COBI    (1,JS,JG) = 0.00_wp
         COALI   (1,JS,JG) = 1.00_wp
         COBTI   (1,JS,JG) = 1.00_wp
         END IF
        END DO
      END DO
      RETURN
  END SUBROUTINE rad_aibi

  !>
  !!  Subroutine aerdis is simplified version from COSMO model (version 4.16).
  !!
  !! @par Revision History
  !! Initial Release by Thorsten Reinhardt, AGeoBw, Offenbach (2011-02-28)
  !!      
  SUBROUTINE aerdis ( klevp1, kbdim, jcs, jce, petah,  pvdaes, pvdael, pvdaeu, pvdaed )
    
    !------------------------------------------------------------------------------
    !
    ! Description:
    !
    ! The module procedure aerdis provides parameters for the vertical distribution
    ! of aerosols (based on the original code of J.F. Geleyn (ECMWF, 4.11.82).
    !
    ! The routine computes the values PVDAE* (* = s, l, u or d for sea, land
    ! urban or desert) of a surfach-normalised vertical distribution of aerosols'
    ! optical depth from the argument petah (vertical coordinate) at klevp1 levels.
    ! It also sets values for non-geograpically weighted total optical depths (at
    ! 55 micrometer wavelength) paeopn for the same four types and similar optical
    ! depths diveded by pressure for bachground well-mixed aerosols of three types
    ! p**bga (** = tr, vo or st for tropospheric, volcanic (stratosperic ashes) or
    ! stratosperic (sulfuric type)). It finally sets values for the power to be
    ! applied to a temperature ratio smaller than two in order to obtain an index
    ! one in the stratosphere and zero in the troposphere with a relatively smooth
    ! transistion (ptrpt), as well as for adsorption coefficients fo water to the
    ! three type of troposperic aerosols (paeadk) with a minimum value ( in the 
    ! whole atmosphere) for the sum of the products paeadk by the optical depths
    ! divided by pressure thickness: paeadm. 
    !
    ! Method:
    !
    ! Straightforward, equivalent heights are given in meters (8434 for the
    ! atmosphere) and tropospheric and stratospheric pressure boundary values
    ! are set at 101325 and 19330 Pascal. 
    !
    !------------------------------------------------------------------------------
    
    ! Subroutine arguments:
    ! --------------------
    
    ! Input data
    ! ----------
    INTEGER, INTENT (IN) ::  &
      & klevp1,         &           ! number of model layer interfaces
      & kbdim,          &
      & jcs,            &
      & jce

    REAL    (wp), INTENT (IN) ::  &
      petah(kbdim,klevp1)    ! normalized vertical coordinate at half levels

    ! Output data
    ! -----------
    REAL    (wp), INTENT (OUT) ::  &
      pvdaes(kbdim,klevp1), & ! normalized vertical distribution (sea)
      pvdael(kbdim,klevp1), & ! normalized vertical distribution (land)
      pvdaeu(kbdim,klevp1), & ! normalized vertical distribution (urban)
      pvdaed(kbdim,klevp1)    ! normalized vertical distrubution (desert)

    ! Local parameters:
    ! -------------
    REAL (wp), PARAMETER  ::  &
      zhss = 8434.0_wp/1000.0_wp ,  & !
      zhsl = 8434.0_wp/1000.0_wp ,  & !
      zhsu = 8434.0_wp/1000.0_wp ,  & !
      zhsd = 8434.0_wp/3000.0_wp      !

    INTEGER :: jc,jk
    REAL(wp) :: log_eta

    !- End of header
    !==============================================================================

    !------------------------------------------------------------------------------
    ! Begin Subroutine aerdis              
    !------------------------------------------------------------------------------

    DO jc=jcs,jce
      pvdaes(jc,1) = 0.0_wp
      pvdael(jc,1) = 0.0_wp
      pvdaeu(jc,1) = 0.0_wp
      pvdaed(jc,1) = 0.0_wp
    ENDDO

!!$  IF(petah(1).NE.0._wp) THEN
!!$     pvdaes(1) = petah(1)**zhss
!!$     pvdael(1) = petah(1)**zhsl
!!$     pvdaeu(1) = petah(1)**zhsu
!!$     pvdaed(1) = petah(1)**zhsd
!!$  END IF

    DO jk=2,klevp1
      DO jc=jcs,jce
        log_eta       = LOG(petah(jc,jk))
        pvdaes(jc,jk) = EXP(zhss*log_eta) ! petah(jc,jk)**zhss
        pvdael(jc,jk) = pvdaes(jc,jk)     ! petah(jc,jk)**zhsl; zhsl is the same as zhss
        pvdaeu(jc,jk) = pvdaes(jc,jk)     ! petah(jc,jk)**zhsu; zhsu is the same as zhss
        pvdaed(jc,jk) = EXP(zhsd*log_eta) ! petah(jc,jk)**zhsd
      ENDDO
    ENDDO

    !------------------------------------------------------------------------------
    ! End of the subroutine 
    !------------------------------------------------------------------------------

  END SUBROUTINE aerdis
  
END MODULE mo_radiation_rg_par

