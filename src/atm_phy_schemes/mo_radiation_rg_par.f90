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

MODULE mo_radiation_rg_par

  USE mo_impl_constants,       ONLY: min_rlcell
  USE mo_impl_constants_grf,   ONLY: grf_bdywidth_c
  USE mo_kind,                 ONLY: wp
  USE mo_loopindices,          ONLY: get_indices_c
  USE mo_model_domain,         ONLY: t_patch
  
  IMPLICIT NONE

  PUBLIC

!  CHARACTER(len=*), PARAMETER :: version = '$Id$'

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
      REAL(wp) :: grenze(2,2,jpspec)
  !             WMIN1      , WMAX1       , WMIN2       , WMAX2
  DATA grenze /  1.5300_wp ,   4.6420_wp , 999._wp     ,  0._wp    , &
                 0.7000_wp ,   1.5300_wp , 999._wp     ,  0._wp    , &
                 0.2451_wp ,   0.7000_wp , 999._wp     ,  0._wp    , &
                20.0000_wp , 104.5150_wp , 999._wp     ,  0._wp    , &
                12.5000_wp ,  20.0000_wp , 999._wp     ,  0._wp    , &
                 8.3333_wp ,   9.0090_wp ,  10.3093_wp , 12.5000_wp, &
                 9.0090_wp ,  10.3093_wp , 999._wp     ,  0._wp    , &
                 4.6420_wp ,   8.3333_wp , 999._wp     ,  0._wp      /

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


  REAL  (wp)              ::           &
  zaea(jpspec,5), &  ! ratio of optical thickness for the absorption in spectral
                     ! interval jpspec  and total optical thickness at 0.55m*1.E-06 
                     ! for an aerosoltyp specified by second array index
  zaes(jpspec,5), &  ! analog for the optical thickness of scattering 
  zaeg(jpspec,5), zaef(jpspec,5)
      
      DATA ZAEA /  &
       .0345_wp, .0511_wp, .0847_wp, .0336_wp, &
       .0499_wp, .0364_wp, .0382_wp, .0260_wp, &
       .0457_wp, .0018_wp, .0015_wp, .1361_wp, &
       .2346_wp, .1177_wp, .0684_wp, .0808_wp, &
       .0707_wp, .0689_wp, .1557_wp, .1258_wp, &
       .1588_wp, .1973_wp, .2766_wp, .1134_wp, &
       .0597_wp, .1077_wp, .2095_wp, .0299_wp, &
       .0456_wp, .0358_wp, .0377_wp, .0304_wp, &
       .0103_wp, .000016_wp,.0000_wp,.0087_wp, &
       .0238_wp, .0511_wp, .0734_wp, .0809_wp /

      DATA ZAES /  &
       .1030_wp, .3977_wp,1.0680_wp, .0084_wp, &
       .0142_wp, .0191_wp, .0234_wp, .0140_wp, &
       .7894_wp, .9734_wp,1.0110_wp, .0307_wp, &
       .0531_wp, .0546_wp, .0839_wp, .2142_wp, &
       .7157_wp, .8698_wp, .8604_wp, .0645_wp, &
       .0781_wp, .1256_wp, .2317_wp, .1409_wp, &
       .0859_wp, .3442_wp, .9496_wp, .0067_wp, &
       .0113_wp, .0153_wp, .0187_wp, .0113_wp, &
       .0467_wp, .3854_wp,1.1008_wp, .0000_wp, &
       .00005_wp,.0004_wp, .0006_wp, .0006_wp /

      DATA ZAEG /  &
       .6562_wp, .6614_wp, .7109_wp, .5043_wp, &
       .6486_wp, .6814_wp, .6489_wp, .7799_wp, &
       .8105_wp, .7906_wp, .7947_wp, .4374_wp, &
       .5203_wp, .7076_wp, .7246_wp, .7535_wp, &
       .6932_wp, .6962_wp, .7402_wp, .4029_wp, &
       .5587_wp, .5618_wp, .4520_wp, .7120_wp, &
       .6462_wp, .6510_wp, .6955_wp, .5041_wp, &
       .6482_wp, .6805_wp, .6477_wp, .7753_wp, &
       .3751_wp, .6353_wp, .7259_wp, .0037_wp, &
       .0083_wp, .0177_wp, .0201_wp, .0332_wp /


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
      NFAST = 1
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
  !! Initializes aerosol (climatologically).
  !!
  !! Taken from COSMO model (version 4.16), adaptated for better vectorization.
  !!
  !! @par Revision History
  !! Initial Release by Thorsten Reinhardt, AGeoBw, Offenbach (2011-02-28)
  !!      
  SUBROUTINE init_aerosol (kbdim,pt_patch,aersea,aerlan,aerurb,aerdes)
    
    INTEGER,            INTENT(in)    :: kbdim

    TYPE(t_patch),      INTENT(in)    :: pt_patch    ! Patch

    REAL(wp),           INTENT(out)   :: &
      & aersea(kbdim,pt_patch%nblks_c),  &
      & aerlan(kbdim,pt_patch%nblks_c),  &
      & aerurb(kbdim,pt_patch%nblks_c),  &
      & aerdes(kbdim,pt_patch%nblks_c)

    ! Local arrays and scalars:
    ! -------------------------
    INTEGER                   ::  &
      jc, jb,                     & ! loop indices
      jzj, jzm1, jzm2, jzm, jzn,  & ! indices for Legendre coefficients
      imn, imnc, imns, jmm, jnn

    INTEGER :: rl_start, rl_end
    INTEGER :: i_startblk, i_endblk    !> blocks
    INTEGER :: i_startidx, i_endidx    !< slices
    INTEGER :: i_nchdom                !< domain index

    REAL (wp)              ::  &
                                ! arrays for the T10 distrubution of
      zaesc(66) , zaess (55) , & ! sea    type aerosols                     
      zaelc(66) , zaels (55) , & ! land   type aerosols
      zaeuc(66) , zaeus (55) , & ! urban  type aerosols
      zaedc(66) , zaeds (55) , & ! desert type aerosols
      zfaes(kbdim,21) , zfael (kbdim,21) , & ! coefficients for spectral
      zfaeu(kbdim,21) , zfaed (kbdim,21) , & ! expansion
      zalp (kbdim,66) ,              & !
      zsinphi(kbdim)   , zcosphi(kbdim)    , & !
      zm, z2m, zre1, ze1, ze2, & !
      zf1m(kbdim), zf2m(kbdim), zn, zn2,     & !
      zsin1, zsin2, zsin3, zsin4, zsin5,  & ! 
      zsin6, zsin7, zsin8, zsin9, zsin10, & !
      zcos1, zcos2, zcos3, zcos4, zcos5,  & ! 
      zcos6, zcos7, zcos8, zcos9, zcos10    !

    !------------------------------------------------------------------------------
    ! Section 0: Data for the Fourier coefficients of the four aerosol types          
    !------------------------------------------------------------------------------

    DATA zaesc/ &
      +.6688E+00_wp,-.1172E+00_wp,-.1013E+00_wp,+.1636E-01_wp,-.3699E-01_wp,+.1775E-01_wp,   &
      -.9635E-02_wp,+.1290E-02_wp,+.4681E-04_wp,-.9106E-04_wp,+.9355E-04_wp,                 &
      -.7076E-01_wp,-.1782E-01_wp,+.1856E-01_wp,+.1372E-01_wp,+.8210E-04_wp,+.2149E-02_wp,   &
      +.4856E-03_wp,+.2231E-03_wp,+.1824E-03_wp,+.1960E-05_wp,                               &
      +.2057E-01_wp,+.2703E-01_wp,+.2424E-01_wp,+.9716E-02_wp,+.1312E-02_wp,-.8846E-03_wp,   &
      -.3347E-03_wp,+.6231E-04_wp,+.6397E-04_wp,                                             &
      -.3341E-02_wp,-.1295E-01_wp,-.4598E-02_wp,+.3242E-03_wp,+.8122E-03_wp,-.2975E-03_wp,   &
      -.7757E-04_wp,+.7793E-04_wp,                                                           &
      +.4455E-02_wp,-.1584E-01_wp,-.2551E-02_wp,+.1174E-02_wp,+.1335E-04_wp,+.5112E-04_wp,   &
      +.5605E-04_wp,                                                                         &
      +.7412E-04_wp,+.1857E-02_wp,-.1917E-03_wp,+.4460E-03_wp,+.1767E-04_wp,-.5281E-04_wp,   &
      -.5043E-03_wp,+.2467E-03_wp,-.2497E-03_wp,-.2377E-04_wp,-.3954E-04_wp,                 &
      +.2666E-03_wp,-.8186E-03_wp,-.1441E-03_wp,-.1904E-04_wp,                               &
      +.3337E-03_wp,-.1696E-03_wp,-.2503E-04_wp,                                             &
      +.1239E-03_wp,-.9983E-04_wp,                                                           &
      -.5283E-04_wp  /

    DATA zaess/                                                                              &
      -.3374E-01_wp,-.3247E-01_wp,-.1012E-01_wp,+.6002E-02_wp,+.5190E-02_wp,+.7784E-03_wp,   &
      -.1090E-02_wp,+.3294E-03_wp,+.1719E-03_wp,-.5866E-05_wp,                               &
      -.4124E-03_wp,-.3742E-01_wp,-.5054E-02_wp,+.3430E-02_wp,+.5513E-03_wp,-.6235E-03_wp,   &
      +.2892E-03_wp,-.9730E-04_wp,+.7078E-04_wp,                                             &
      -.3300E-01_wp,+.5104E-03_wp,-.2156E-02_wp,-.3194E-02_wp,-.5079E-03_wp,-.5517E-03_wp,   &
      +.4632E-04_wp,+.5369E-04_wp,                                                           &
      -.2731E-01_wp,+.5126E-02_wp,+.2241E-02_wp,-.5789E-03_wp,-.3048E-03_wp,-.1774E-03_wp,   &
      +.1946E-05_wp,                                                                         &
      -.8247E-02_wp,+.2338E-02_wp,+.1021E-02_wp,+.1575E-04_wp,+.2612E-05_wp,+.1995E-04_wp,   &
      -.1319E-02_wp,+.1384E-02_wp,-.4159E-03_wp,-.2337E-03_wp,+.5764E-04_wp,                 &
      +.1495E-02_wp,-.3727E-03_wp,+.6075E-04_wp,-.4642E-04_wp,                               &
      +.5368E-03_wp,-.7619E-04_wp,+.3774E-04_wp,                                             &
      +.1206E-03_wp,-.4104E-06_wp,                                                           &
      +.2158E-04_wp  /

    DATA zaelc/                                                                              &
      +.1542E+00_wp,+.8245E-01_wp,-.1879E-03_wp,+.4864E-02_wp,-.5527E-02_wp,-.7966E-02_wp,   &
      -.2683E-02_wp,-.2011E-02_wp,-.8889E-03_wp,-.1058E-03_wp,-.1614E-04_wp,                 &
      +.4206E-01_wp,+.1912E-01_wp,-.9476E-02_wp,-.6780E-02_wp,+.1767E-03_wp,-.5422E-03_wp,   &
      -.7753E-03_wp,-.2106E-03_wp,-.9870E-04_wp,-.1721E-04_wp,                               &
      -.9536E-02_wp,-.9580E-02_wp,-.1050E-01_wp,-.5747E-02_wp,-.1282E-02_wp,+.2248E-03_wp,   &
      +.1694E-03_wp,-.4782E-04_wp,-.2441E-04_wp,                                             &
      +.5781E-03_wp,+.6212E-02_wp,+.1921E-02_wp,-.1102E-02_wp,-.8145E-03_wp,+.2497E-03_wp,   &
      +.1539E-03_wp,-.2538E-04_wp,                                                           &
      -.3993E-02_wp,+.9777E-02_wp,+.4837E-03_wp,-.1304E-02_wp,+.2417E-04_wp,-.1370E-04_wp,   &
      -.3731E-05_wp,                                                                         &
      +.1922E-02_wp,-.5167E-03_wp,+.4295E-03_wp,-.1888E-03_wp,+.2427E-04_wp,+.4012E-04_wp,   &
      +.1529E-02_wp,-.2120E-03_wp,+.8166E-04_wp,+.2579E-04_wp,+.3488E-04_wp,                 &
      +.2140E-03_wp,+.2274E-03_wp,-.3447E-05_wp,-.1075E-04_wp,                               &
      -.1018E-03_wp,+.2864E-04_wp,+.3442E-04_wp,                                             &
      -.1002E-03_wp,+.7117E-04_wp,                                                           &
      +.2045E-04_wp  /

    DATA zaels/                                                                              &
      +.1637E-01_wp,+.1935E-01_wp,+.1080E-01_wp,+.2784E-02_wp,+.1606E-03_wp,+.1860E-02_wp,   &
      +.1263E-02_wp,-.2707E-03_wp,-.2290E-03_wp,-.9761E-05_wp,                               &
      -.7317E-02_wp,+.2465E-01_wp,+.6799E-02_wp,-.1913E-02_wp,+.1382E-02_wp,+.6691E-03_wp,   &
      +.1414E-03_wp,+.3527E-04_wp,-.5210E-04_wp,                                             &
      +.1873E-01_wp,+.2977E-02_wp,+.4650E-02_wp,+.2509E-02_wp,+.3680E-03_wp,+.1481E-03_wp,   &
      -.6594E-04_wp,-.5634E-04_wp,                                                           &
      +.1592E-01_wp,-.1875E-02_wp,-.1093E-02_wp,+.3022E-03_wp,+.2625E-03_wp,+.3252E-04_wp,   &
      -.3803E-04_wp,                                                                         &
      +.4218E-02_wp,-.1843E-02_wp,-.1351E-02_wp,-.2952E-03_wp,-.8171E-05_wp,-.1473E-04_wp,   &
      +.9076E-03_wp,-.1057E-02_wp,+.2676E-03_wp,+.1307E-03_wp,-.3628E-04_wp,                 &
      -.9158E-03_wp,+.4335E-03_wp,+.2927E-04_wp,+.6602E-04_wp,                               &
      -.3570E-03_wp,+.5760E-04_wp,-.3465E-04_wp,                                             &
      -.8535E-04_wp,-.2011E-04_wp,                                                           &
      +.6612E-06_wp  /  

    DATA zaeuc/                                                                              &
      +.8005E-01_wp,+.7095E-01_wp,+.2014E-01_wp,-.1412E-01_wp,-.2425E-01_wp,-.1332E-01_wp,   &
      -.2904E-02_wp,+.5068E-03_wp,+.9369E-03_wp,+.4114E-03_wp,+.7549E-04_wp,                 &
      +.1922E-01_wp,+.2534E-01_wp,+.2088E-01_wp,+.1064E-01_wp,+.1063E-02_wp,-.2526E-02_wp,   &
      -.2091E-02_wp,-.9660E-03_wp,-.2030E-03_wp,+.3865E-04_wp,                               &
      -.9900E-02_wp,-.5964E-02_wp,+.2223E-02_wp,+.4941E-02_wp,+.3277E-02_wp,+.1038E-02_wp,   &
      -.1480E-03_wp,-.2844E-03_wp,-.1208E-03_wp,                                             &
      +.3999E-02_wp,+.6282E-02_wp,+.2813E-02_wp,+.1475E-02_wp,+.4571E-03_wp,-.1349E-03_wp,   &
      -.9011E-04_wp,-.1936E-04_wp,                                                           &
      +.1994E-02_wp,+.3540E-02_wp,+.8837E-03_wp,+.1992E-03_wp,+.3092E-04_wp,-.7979E-04_wp,   &
      -.2664E-04_wp,                                                                         &
      -.5006E-04_wp,+.6447E-03_wp,+.5550E-03_wp,+.1197E-03_wp,+.6657E-04_wp,+.1488E-04_wp,   &
      -.9141E-04_wp,-.2896E-03_wp,-.1561E-03_wp,-.6524E-04_wp,-.1559E-04_wp,                 &
      -.1082E-03_wp,-.4126E-03_wp,-.1732E-03_wp,-.8286E-04_wp,                               &
      -.1993E-04_wp,+.3850E-04_wp,+.2870E-04_wp,                                             &
      +.4493E-04_wp,+.4721E-04_wp,                                                           &
      +.1338E-04_wp  /

    DATA zaeus/                                                                              &
      +.6646E-02_wp,+.8373E-02_wp,+.5463E-02_wp,+.4554E-02_wp,+.3301E-02_wp,+.5725E-03_wp,   &
      -.7482E-03_wp,-.6222E-03_wp,-.2603E-03_wp,-.5127E-04_wp,                               &
      -.3849E-04_wp,+.9741E-02_wp,+.8190E-02_wp,+.5712E-02_wp,+.3039E-02_wp,+.5290E-03_wp,   &
      -.2044E-03_wp,-.2309E-03_wp,-.1160E-03_wp,                                             &
      +.9160E-02_wp,+.1286E-01_wp,+.1170E-01_wp,+.5491E-02_wp,+.1393E-02_wp,-.6288E-04_wp,   &
      -.2715E-03_wp,-.1047E-03_wp,                                                           &
      +.4873E-02_wp,+.3545E-02_wp,+.3069E-02_wp,+.1819E-02_wp,+.6947E-03_wp,+.1416E-03_wp,   &
      -.1538E-04_wp,                                                                         &
      -.4351E-03_wp,-.1907E-02_wp,-.5774E-03_wp,-.2247E-03_wp,+.5345E-04_wp,+.9052E-04_wp,   &
      -.3972E-04_wp,-.9665E-04_wp,+.7912E-04_wp,-.1094E-04_wp,-.6776E-05_wp,                 &
      +.2724E-03_wp,+.1973E-03_wp,+.6837E-04_wp,+.4313E-04_wp,                               &
      -.7174E-05_wp,+.8527E-05_wp,-.2160E-05_wp,                                             &
      -.7852E-04_wp,+.3453E-06_wp,                                                           &
      -.2402E-05_wp  /

    DATA zaedc/                                                                              &
      +.2840E-01_wp,+.1775E-01_wp,-.1069E-01_wp,-.1553E-01_wp,-.3299E-02_wp,+.3583E-02_wp,   &
      +.2274E-02_wp,+.5767E-04_wp,-.3678E-03_wp,-.1050E-03_wp,+.2133E-04_wp,                 &
      +.2326E-01_wp,+.1566E-01_wp,-.3130E-02_wp,-.8253E-02_wp,-.2615E-02_wp,+.1247E-02_wp,   &
      +.1059E-02_wp,+.1196E-03_wp,-.1303E-03_wp,-.5094E-04_wp,                               &
      +.1185E-01_wp,+.7238E-02_wp,-.1562E-02_wp,-.3665E-02_wp,-.1182E-02_wp,+.4678E-03_wp,   &
      +.4448E-03_wp,+.8307E-04_wp,-.3468E-04_wp,                                             &
      +.5273E-02_wp,+.3037E-02_wp,-.4014E-03_wp,-.1202E-02_wp,-.4647E-03_wp,+.5148E-04_wp,   &
      +.1014E-03_wp,+.2996E-04_wp,                                                           &
      +.2505E-02_wp,+.1495E-02_wp,+.2438E-03_wp,-.1223E-03_wp,-.7669E-04_wp,-.1638E-04_wp,   &
      +.1869E-05_wp,                                                                         &
      +.1094E-02_wp,+.6131E-03_wp,+.1508E-03_wp,+.1765E-04_wp,+.1360E-05_wp,-.7998E-06_wp,   &
      +.4475E-03_wp,+.2737E-03_wp,+.6430E-04_wp,-.6759E-05_wp,-.6761E-05_wp,                 &
      +.1992E-03_wp,+.1531E-03_wp,+.4828E-04_wp,+.5103E-06_wp,                               &
      +.7454E-04_wp,+.5917E-04_wp,+.2152E-04_wp,                                             &
      +.9300E-05_wp,+.9790E-05_wp,                                                           &
      -.8853E-05_wp /

    DATA zaeds/                                                                              &
      +.9815E-02_wp,+.8436E-02_wp,+.1087E-02_wp,-.2717E-02_wp,-.1755E-02_wp,-.1559E-03_wp,   &
      +.2367E-03_wp,+.8808E-04_wp,+.2001E-05_wp,-.1244E-05_wp,                               &
      +.1041E-01_wp,+.8039E-02_wp,+.1005E-02_wp,-.1981E-02_wp,-.1090E-02_wp,+.1595E-05_wp,   &
      +.1787E-03_wp,+.4644E-04_wp,-.1052E-04_wp,                                             &
      +.6593E-02_wp,+.3983E-02_wp,-.1527E-03_wp,-.1235E-02_wp,-.5078E-03_wp,+.3649E-04_wp,   &
      +.1005E-03_wp,+.3182E-04_wp,                                                           &
      +.3225E-02_wp,+.1672E-02_wp,-.7752E-04_wp,-.4312E-03_wp,-.1872E-03_wp,-.1666E-04_wp,   &
      +.1872E-04_wp,                                                                         &
      +.1133E-02_wp,+.5643E-03_wp,+.7747E-04_wp,-.2980E-04_wp,-.2092E-04_wp,-.8590E-05_wp,   &
      +.2988E-03_wp,+.6714E-04_wp,-.6249E-05_wp,+.1052E-04_wp,+.8790E-05_wp,                 &
      +.1569E-03_wp,-.1175E-04_wp,-.3033E-04_wp,-.9777E-06_wp,                               &
      +.1101E-03_wp,+.6827E-05_wp,-.1023E-04_wp,                                             &
      +.4231E-04_wp,+.4905E-05_wp,                                                           &
      +.6229E-05_wp  /


    !------------------------------------------------------------------------------
    ! Section 2: Calculation of the inverse Legendre and Fourier transformation  
    !------------------------------------------------------------------------------

    zaea=RESHAPE((/0.0477_wp, 0.0875_wp,  0.1198_wp, 0.0458_wp, &
      0.0387_wp, 0.0439_wp,  0.0599_wp, 0.0396_wp, &
      0.0381_wp, 0.0129_wp,  0.0130_wp, 0.1304_wp, &
      0.1757_wp, 0.0949_wp,  0.0653_wp, 0.0795_wp, &
      0.0962_wp, 0.2046_wp,  0.4116_wp, 0.0169_wp, &
      0.0204_wp, 0.0263_wp,  0.0348_wp, 0.0361_wp, &
      0.0030_wp, 0.0271_wp,  0.0613_wp, 0.0118_wp, &
      0.0160_wp, 0.0231_wp,  0.0287_wp, 0.0127_wp, &
      0.0103_wp, 0.000016_wp,0.0000_wp, 0.0087_wp, &
      0.0238_wp, 0.0511_wp,  0.0734_wp, 0.0809_wp/),(/8,5/))

    zaes=RESHAPE((/0.1407_wp, 0.4256_wp,  1.0066_wp, 0.0279_wp, &
      0.0391_wp, 0.0445_wp,  0.0485_wp, 0.0362_wp, &
      0.6746_wp, 0.8761_wp,  1.0139_wp, 0.0443_wp, &
      0.0624_wp, 0.0921_wp,  0.1491_wp, 0.2327_wp, &
      0.0605_wp, 0.2761_wp,  0.7449_wp, 0.0023_wp, &
      0.0034_wp, 0.0051_wp,  0.0065_wp, 0.0045_wp, &
      0.0284_wp, 0.5524_wp,  0.9683_wp, 0.0001_wp, &
      0.0004_wp, 0.0024_wp,  0.0049_wp, 0.0030_wp, &
      0.0467_wp, 0.3854_wp,  1.1008_wp, 0.0000_wp, &
      0.00005_wp,0.0004_wp,  0.0006_wp, 0.0006_wp/),(/8,5/))

    zaeg=RESHAPE((/0.6989_wp, 0.6329_wp,  0.6418_wp, 0.6243_wp, &
      0.7299_wp, 0.7430_wp,  0.7086_wp, 0.8569_wp, &
      0.7833_wp, 0.7575_wp,  0.7456_wp, 0.4997_wp, &
      0.6130_wp, 0.7440_wp,  0.7426_wp, 0.7590_wp, &
      0.5753_wp, 0.5867_wp,  0.5957_wp, 0.6027_wp, &
      0.6766_wp, 0.6117_wp,  0.5439_wp, 0.6905_wp, &
      0.5170_wp, 0.6674_wp,  0.7004_wp, 0.0340_wp, &
      0.0570_wp, 0.1289_wp,  0.1597_wp, 0.1906_wp, &
      0.3751_wp, 0.6353_wp,  0.7259_wp, 0.0037_wp, &
      0.0083_wp, 0.0177_wp,  0.0201_wp, 0.0332_wp/),(/8,5/))

    i_nchdom  = MAX(1,pt_patch%n_childdom)

    !in order to account for mesh refinement
    rl_start = grf_bdywidth_c+1
    rl_end   = min_rlcell

    i_startblk = pt_patch%cells%start_blk(rl_start,1)
    i_endblk   = pt_patch%cells%end_blk(rl_end,i_nchdom)

    zalp (:,1) = 1.0_wp

    DO jb = i_startblk, i_endblk

      CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
        &                         i_startidx, i_endidx, rl_start, rl_end)

      DO jc=i_startidx,i_endidx

        ! Calculation of the values zalp for the sine of latitude (zsinphi) of the
        ! normalized Legendre associated functions. The limit wave number is 10.

        zsinphi(jc)  = SIN(pt_patch%cells%center(jc,jb)%lat )
        zcosphi(jc)  = SQRT(1._wp-zsinphi(jc)**2)

        zf1m(jc)     = SQRT(3.0_wp)
        !        zalp (1) = 1.0_wp
        zalp (jc,2) = zf1m(jc)*zsinphi(jc)

      ENDDO

      jzj      = 2

      wave_number_loop : DO jzm1 = 1, 11

        jzm  = jzm1-1
        zm   = REAL(jzm,wp)
        z2m  = zm + zm
        zre1 = SQRT(z2m+3.0_wp)
        ze1  = 1.0_wp/zre1
        IF (jzm.NE.0) THEN
          jzj       = jzj + 1
          DO jc=i_startidx,i_endidx
            zf2m(jc)      = zf1m(jc)*zcosphi(jc)/SQRT(z2m)
            zf1m(jc)      = zf2m(jc)*zre1
            zalp(jc,jzj) = zf2m(jc)
          ENDDO
          IF(jzm ==10) CYCLE wave_number_loop
          jzj       = jzj + 1
          DO jc=i_startidx,i_endidx
            zalp(jc,jzj) = zf1m(jc)*zsinphi(jc)
          ENDDO
          IF(jzm1==10) CYCLE wave_number_loop
        ENDIF
        jzm2 = jzm+2
        DO jzn = jzm2, 10
          zn        = REAL(jzn,wp)
          zn2       = zn**2
          ze2       = SQRT( (4.0_wp*zn2-1.0_wp)/(zn2-zm**2) )
          jzj       = jzj+1
          DO jc=i_startidx,i_endidx
            zalp(jc,jzj) = ze2*(zsinphi(jc)*zalp(jc,jzj-1)-ze1*zalp(jc,jzj-2))
          ENDDO
          ze1       = 1.0_wp/ze2
        ENDDO
      ENDDO wave_number_loop


      ! Legendre transform of aerosols

      DO jc=i_startidx,i_endidx
        zfaes(jc,:) = 0.0_wp
        zfael(jc,:) = 0.0_wp
        zfaeu(jc,:) = 0.0_wp
        zfaed(jc,:) = 0.0_wp
      ENDDO !jc

      imn  = 0
      imnc = 0
      imns = 0

      DO jmm = 1, 11
        imn  = imn  + 1
        DO jnn = jmm, 11
          imnc       = imnc + 1
          DO jc=i_startidx,i_endidx
            zfaes(jc,imn) = zfaes(jc,imn)+zalp(jc,imnc)*zaesc(imnc)
            zfael(jc,imn) = zfael(jc,imn)+zalp(jc,imnc)*zaelc(imnc)
            zfaeu(jc,imn) = zfaeu(jc,imn)+zalp(jc,imnc)*zaeuc(imnc)
            zfaed(jc,imn) = zfaed(jc,imn)+zalp(jc,imnc)*zaedc(imnc)
          ENDDO
        ENDDO
        IF(jmm.NE.1) THEN
          imn  = imn+1
          DO jnn = jmm, 11
            imns       = imns + 1
            DO jc=i_startidx,i_endidx
              zfaes(jc,imn) = zfaes(jc,imn)+zalp(jc,imns+11)*zaess(imns)
              zfael(jc,imn) = zfael(jc,imn)+zalp(jc,imns+11)*zaels(imns)
              zfaeu(jc,imn) = zfaeu(jc,imn)+zalp(jc,imns+11)*zaeus(imns)
              zfaed(jc,imn) = zfaed(jc,imn)+zalp(jc,imns+11)*zaeds(imns)
            ENDDO
          ENDDO
        ENDIF
      ENDDO

      DO jc=i_startidx,i_endidx

        ! Inverse Fourier transformation

        zcos1   = COS(pt_patch%cells%center(jc,jb)%lon )
        zsin1   = SIN(pt_patch%cells%center(jc,jb)%lon )
        zcos2   = zcos1*zcos1 - zsin1*zsin1
        zsin2   = zsin1*zcos1 + zcos1*zsin1
        zcos3   = zcos2*zcos1 - zsin2*zsin1
        zsin3   = zsin2*zcos1 + zcos2*zsin1
        zcos4   = zcos3*zcos1 - zsin3*zsin1
        zsin4   = zsin3*zcos1 + zcos3*zsin1
        zcos5   = zcos4*zcos1 - zsin4*zsin1
        zsin5   = zsin4*zcos1 + zcos4*zsin1
        zcos6   = zcos5*zcos1 - zsin5*zsin1
        zsin6   = zsin5*zcos1 + zcos5*zsin1
        zcos7   = zcos6*zcos1 - zsin6*zsin1
        zsin7   = zsin6*zcos1 + zcos6*zsin1
        zcos8   = zcos7*zcos1 - zsin7*zsin1
        zsin8   = zsin7*zcos1 + zcos7*zsin1
        zcos9   = zcos8*zcos1 - zsin8*zsin1
        zsin9   = zsin8*zcos1 + zcos8*zsin1
        zcos10  = zcos9*zcos1 - zsin9*zsin1
        zsin10  = zsin9*zcos1 + zcos9*zsin1

        aersea(jc,jb) = zfaes(jc,1) + 2._wp* ( zfaes(jc,2 ) * zcos1 + zfaes(jc,3 ) * zsin1  &
          + zfaes(jc,4 ) * zcos2 + zfaes(jc,5 ) * zsin2    &
          + zfaes(jc,6 ) * zcos3 + zfaes(jc,7 ) * zsin3    &
          + zfaes(jc,8 ) * zcos4 + zfaes(jc,9 ) * zsin4    &
          + zfaes(jc,10) * zcos5 + zfaes(jc,11) * zsin5    &
          + zfaes(jc,12) * zcos6 + zfaes(jc,13) * zsin6    &
          + zfaes(jc,14) * zcos7 + zfaes(jc,15) * zsin7    &
          + zfaes(jc,16) * zcos8 + zfaes(jc,17) * zsin8    &
          + zfaes(jc,18) * zcos9 + zfaes(jc,19) * zsin9    &
          + zfaes(jc,20) * zcos10+ zfaes(jc,21) * zsin10 )

        aerlan(jc,jb) = zfael(jc,1) + 2._wp* ( zfael(jc,2 ) * zcos1 + zfael(jc,3 ) * zsin1  &
          + zfael(jc,4 ) * zcos2 + zfael(jc,5 ) * zsin2    &
          + zfael(jc,6 ) * zcos3 + zfael(jc,7 ) * zsin3    &
          + zfael(jc,8 ) * zcos4 + zfael(jc,9 ) * zsin4    &
          + zfael(jc,10) * zcos5 + zfael(jc,11) * zsin5    &
          + zfael(jc,12) * zcos6 + zfael(jc,13) * zsin6    &
          + zfael(jc,14) * zcos7 + zfael(jc,15) * zsin7    &
          + zfael(jc,16) * zcos8 + zfael(jc,17) * zsin8    &
          + zfael(jc,18) * zcos9 + zfael(jc,19) * zsin9    &
          + zfael(jc,20) * zcos10+ zfael(jc,21) * zsin10 )

        aerurb(jc,jb) = zfaeu(jc,1) + 2._wp* ( zfaeu(jc,2 ) * zcos1 + zfaeu(jc,3 ) * zsin1  &
          + zfaeu(jc,4 ) * zcos2 + zfaeu(jc,5 ) * zsin2    &
          + zfaeu(jc,6 ) * zcos3 + zfaeu(jc,7 ) * zsin3    &
          + zfaeu(jc,8 ) * zcos4 + zfaeu(jc,9 ) * zsin4    &
          + zfaeu(jc,10) * zcos5 + zfaeu(jc,11) * zsin5    &
          + zfaeu(jc,12) * zcos6 + zfaeu(jc,13) * zsin6    &
          + zfaeu(jc,14) * zcos7 + zfaeu(jc,15) * zsin7    &
          + zfaeu(jc,16) * zcos8 + zfaeu(jc,17) * zsin8    &
          + zfaeu(jc,18) * zcos9 + zfaeu(jc,19) * zsin9    &
          + zfaeu(jc,20) * zcos10+ zfaeu(jc,21) * zsin10 )

        aerdes(jc,jb) = zfaed(jc,1) + 2._wp* ( zfaed(jc,2 ) * zcos1 + zfaed(jc,3 ) * zsin1  &
          + zfaed(jc,4 ) * zcos2 + zfaed(jc,5 ) * zsin2    &
          + zfaed(jc,6 ) * zcos3 + zfaed(jc,7 ) * zsin3    &
          + zfaed(jc,8 ) * zcos4 + zfaed(jc,9 ) * zsin4    &
          + zfaed(jc,10) * zcos5 + zfaed(jc,11) * zsin5    &
          + zfaed(jc,12) * zcos6 + zfaed(jc,13) * zsin6    &
          + zfaed(jc,14) * zcos7 + zfaed(jc,15) * zsin7    &
          + zfaed(jc,16) * zcos8 + zfaed(jc,17) * zsin8    &
          + zfaed(jc,18) * zcos9 + zfaed(jc,19) * zsin9    &
          + zfaed(jc,20) * zcos10+ zfaed(jc,21) * zsin10 )

        aersea(jc,jb) = MAX( 0.0_wp, MIN( 1.0_wp, aersea(jc,jb) ) )
        aerlan(jc,jb) = MAX( 0.0_wp, MIN( 1.0_wp, aerlan(jc,jb) ) )
        aerurb(jc,jb) = MAX( 0.0_wp, MIN( 1.0_wp, aerurb(jc,jb) ) )
        aerdes(jc,jb) = MAX( 0.0_wp, MIN( 1.0_wp, aerdes(jc,jb) ) )

      ENDDO !jc

    ENDDO !jb 

  END SUBROUTINE init_aerosol

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
        pvdaes(jc,jk) = petah(jc,jk)**zhss
        pvdael(jc,jk) = petah(jc,jk)**zhsl
        pvdaeu(jc,jk) = petah(jc,jk)**zhsu
        pvdaed(jc,jk) = petah(jc,jk)**zhsd
      ENDDO
    ENDDO

    !------------------------------------------------------------------------------
    ! End of the subroutine 
    !------------------------------------------------------------------------------

  END SUBROUTINE aerdis
  
END MODULE mo_radiation_rg_par

