! $RCSfile$
! $Revision$ $Date$
!
!>
!! <Short description of module for listings and indices>
!!
!! *fcttrm* - Functions used by sucst for inital settings
!! *SUCST * - Routine to initialize the constants of the model.
!! *SUCUMF* - THIS ROUTINE DEFINES DISPOSABLE PARAMETERS FOR MASSFLUX SCHEME
!! *SUPHLI* - SET-UP ROUTINE FOR COMMON BLOCK *YOEPHLI*
!! *SUVDF*  - SET-UP ROUTINE FOR COMMON BLOCK *YOEVDF*
!! *SUVDFS* - INITIALIZES COMMON BLOCK *YOEVDFS*

!! @author    J.F. MAHFOUF         E.C.M.W.F.      1996-06-23
!! @author A.C.M. BELJAARS         E.C.M.W.F.      1989-11-02
!! @author A.C.M. BELJAARS         E.C.M.W.F.      1990-03-26
!!
!!
!! @par Revision History
!! Implementation into GME/ICON by Kristina Froehlich, DWD (2010-06-11)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_cuparameters

#ifdef __ICON__
  USE mo_kind,       ONLY: jpim=>i4, jprb=>wp
  USE mo_exception,  ONLY: message_text, message
  USE mo_physical_constants, ONLY: rday => rdaylen
  USE mo_nwp_parameters,  ONLY: t_phy_params
  USE mo_nwp_tuning_config, ONLY: tune_entrorg, tune_rhebc_land, tune_rhebc_ocean, tune_rcucov, &
    tune_texc, tune_qexc, tune_rhebc_land_trop, tune_rhebc_ocean_trop, tune_rcucov_trop
#endif

#ifdef __GME__
!  USE parkind1,      ONLY: jpim,     jprb
#endif
  
!  USE yomhook,       ONLY: lhook,    dr_hook
  
  IMPLICIT NONE

  PRIVATE


  
  !*    Common of physical constants of IFS
  ! A1.0 Fundamental constants
  REAL(KIND=jprb) :: rpi
  REAL(KIND=jprb) :: rclum
  REAL(KIND=jprb) :: rhpla
  REAL(KIND=jprb) :: rkbol
  REAL(KIND=jprb) :: rnavo
  ! A1.1 Astronomical constants
#ifdef __GME__
  REAL(KIND=jprb) :: rday
#endif
  REAL(KIND=jprb) :: rea
  REAL(KIND=jprb) :: repsm
  REAL(KIND=jprb) :: rsiyea
  REAL(KIND=jprb) :: rsiday
  REAL(KIND=jprb) :: romega
  ! A1.2 Geoide
  REAL(KIND=jprb) :: ra
  REAL(KIND=jprb) :: rg
  REAL(KIND=jprb) :: r1sa
  ! A1.3 Radiation
  REAL(KIND=jprb) :: rsigma
  REAL(KIND=jprb) :: ri0
  ! A1.4 Thermodynamic gas phase
  REAL(KIND=jprb) :: r
  REAL(KIND=jprb) :: rmd
  REAL(KIND=jprb) :: rmv
  REAL(KIND=jprb) :: rmo3
  REAL(KIND=jprb) :: rd
  REAL(KIND=jprb) :: rv
  REAL(KIND=jprb) :: rcpd
  REAL(KIND=jprb) :: rcpv
  REAL(KIND=jprb) :: rcvd
  REAL(KIND=jprb) :: rcvv
  REAL(KIND=jprb) :: rkappa
  REAL(KIND=jprb) :: retv
  ! A1.5,6 Thermodynamic liquid,solid phases
  REAL(KIND=jprb) :: rcw
  REAL(KIND=jprb) :: rcs
  ! A1.7 Thermodynamic transition of phase
  REAL(KIND=jprb) :: rlvtt
  REAL(KIND=jprb) :: rlstt
  REAL(KIND=jprb) :: rlvzer
  REAL(KIND=jprb) :: rlszer
  REAL(KIND=jprb) :: rlmlt
  REAL(KIND=jprb) :: rtt
  REAL(KIND=jprb) :: ratm
  REAL(KIND=jprb) :: rdt
  ! A1.8 Curve of saturation
  REAL(KIND=jprb) :: restt
  REAL(KIND=jprb) :: ralpw
  REAL(KIND=jprb) :: rbetw
  REAL(KIND=jprb) :: rgamw
  REAL(KIND=jprb) :: ralps
  REAL(KIND=jprb) :: rbets
  REAL(KIND=jprb) :: rgams
  REAL(KIND=jprb) :: ralpd
  REAL(KIND=jprb) :: rbetd
  REAL(KIND=jprb) :: rgamd


  ! Description:
  !     ------------------------------------------------------------------
  !*    ** *YOEVDF* CONTAINS CONSTANTS NEEDED BY *VDF....*
  !     ------------------------------------------------------------------
  !
  REAL(KIND=jprb) :: rlam
  REAL(KIND=jprb) :: rkap
  REAL(KIND=jprb) :: rchar
  REAL(KIND=jprb) :: rvdifts
  REAL(KIND=jprb) :: rz0ice
  REAL(KIND=jprb) :: repdu2
  REAL(KIND=jprb) :: repust
  REAL(KIND=jprb) :: rsez0h
  REAL(KIND=jprb) :: rsez0q
  REAL(KIND=jprb) :: rnum
  REAL(KIND=jprb) :: rnuh
  REAL(KIND=jprb) :: rnuq
  REAL(KIND=jprb) :: rentr
  REAL(KIND=jprb) :: rpar
  REAL(KIND=jprb) :: rpar1
  REAL(KIND=jprb) :: rparsrf
  REAL(KIND=jprb) :: rparzi
  REAL(KIND=jprb) :: rlamsk
  LOGICAL lelwdd
  
  !**   ** *YOEDFS* CONTAINS STABILITY FUNCTION TABLES FOR *VDF...*
  
  !     A.C.M. BELJAARS   E.C.M.W.F.       26/03/90.
  
  !      NAME      TYPE        PURPOSE
  !      ----      ----        -------
  
  !     *RCHBA*     REAL       *CONSTANT A IN *HOLTSLAG AND *DEBRUIN
  !                            FUNCTIONS FOR STABLE SITUATIONS
  !     *RCHBB*     REAL       *CONSTANT B IN *HB* FUNCTIONS
  !     *RCHBC*     REAL       *CONSTANT C IN *HB* FUNCTIONS
  !     *RCHBD*     REAL       *CONSTANT D IN *HB* FUNCTIONS
  !     *RCHB23A    REAL       2./3.*A IN *HB* FUNCTIONS
  !     *RCHBBCD    REAL       B*C/D IN *HB* FUNCTIONS
  !     *RCHBCD     REAL       C/D IN *HB* FUNCTIONS
  !     *RCHETA*    REAL       CONSTANT IN THE *HOGSTROM *ELLISON *TURNER
  !                            FUNCTIONS FOR STABLY STRATIFIED TURBULENCE
  !     *RCHETB*    REAL       CONSTANT IN THE *HET* FUNCTIONS
  !     *RCHETC*    REAL       CONSTANT IN *HET* FUNCTIONS
  !     *RCHBHDL*   REAL       MAXIM ZNLEV/L FOR STABLE BOUNDARY LAYER
  !     *RCDHALF    REAL       CONSTANT IN *DYER AND *HICKS FORMULAE
  !                            FOR UNSTABLE SITUATIONS
  !     *RCDHPI2    REAL       PI/2.
  !     *RIMAX*     REAL       *MAXIMIM RICHARDSON NUMBER TABULATED
  !     *DRITBL*    REAL       *INCREMENT OF THE RICHARDSON NUMBER
  !                            BETWEEN TABULATED VALUES.
  !     *DRI26*     REAL       DRITBL**2/6.
  !     *RITBL*     REAL ARRAY *TABULATED ETA-VALUES (Z/L) AS A FUNCTION
  !                            OF THE RICHARDSON NUMBER FOR STABLE CASES.
  !     *ARITBL*    REAL ARRAY *SECOND DERIVATIVES OF TABULATED FUNCTION
  !                            FOR SPLINE INTERPOLATION.
  !     ------------------------------------------------------------------
  
  INTEGER(KIND=jpim), PARAMETER :: jpritbl=101
  REAL(KIND=jprb) :: ritbl(jpritbl)
  REAL(KIND=jprb) :: aritbl(jpritbl)
  REAL(KIND=jprb) :: rchba
  REAL(KIND=jprb) :: rchbb
  REAL(KIND=jprb) :: rchbc
  REAL(KIND=jprb) :: rchbd
  REAL(KIND=jprb) :: rchb23a
  REAL(KIND=jprb) :: rchbbcd
  REAL(KIND=jprb) :: rchbcd
  REAL(KIND=jprb) :: rcheta
  REAL(KIND=jprb) :: rchetb
  REAL(KIND=jprb) :: rchetc
  REAL(KIND=jprb) :: rchbhdl
  REAL(KIND=jprb) :: rcdhalf
  REAL(KIND=jprb) :: rcdhpi2
  REAL(KIND=jprb) :: rimax
  REAL(KIND=jprb) :: dritbl
  REAL(KIND=jprb) :: dri26
  
  !     J.-J. MORCRETTE                   91/07/14  ADAPTED TO I.F.S.
  
  !      NAME     TYPE      PURPOSE
  !      ----     ----      -------
  
  !     *R__ES*   REAL      *CONSTANTS USED FOR COMPUTATION OF SATURATION
  !                         MIXING RATIO OVER LIQUID WATER(*R_LES*) OR
  !                         ICE(*R_IES*).
  !     *RVTMP2*  REAL      *RVTMP2=RCPV/RCPD-1.
  !     *RHOH2O*  REAL      *DENSITY OF LIQUID WATER.   (RATM/100.)
  !     *R5ALVCP* REAL      *R5LES*RLVTT/RCPD
  !     *R5ALSCP* REAL      *R5IES*RLSTT/RCPD
  !     *RALVDCP* REAL      *RLVTT/RCPD
  !     *RALSDCP* REAL      *RLSTT/RCPD
  !     *RALFDCP* REAL      *RLMLT/RCPD
  !     *RTWAT*   REAL      *RTWAT=RTT
  !     *RTBER*   REAL      *RTBER=RTT-0.05
  !     *RTBERCU  REAL      *RTBERCU=RTT-5.0
  !     *RTICE*   REAL      *RTICE=RTT-0.1
  !     *RTICECU* REAL      *RTICECU=RTT-23.0
  !     *RTWAT_RTICE_R*   REAL      *RTWAT_RTICE_R=1./(RTWAT-RTICE)
  !     *RTWAT_RTICECU_R* REAL      *RTWAT_RTICECU_R=1./(RTWAT-RTICECU)
  
  !       ----------------------------------------------------------------
  REAL(KIND=jprb) :: r2es
  REAL(KIND=jprb) :: r3les
  REAL(KIND=jprb) :: r3ies
  REAL(KIND=jprb) :: r4les
  REAL(KIND=jprb) :: r4ies
  REAL(KIND=jprb) :: r5les
  REAL(KIND=jprb) :: r5ies
  REAL(KIND=jprb) :: rvtmp2
  REAL(KIND=jprb) :: rhoh2o
  REAL(KIND=jprb) :: r5alvcp
  REAL(KIND=jprb) :: r5alscp
  REAL(KIND=jprb) :: ralvdcp
  REAL(KIND=jprb) :: ralsdcp
  REAL(KIND=jprb) :: ralfdcp
  REAL(KIND=jprb) :: rtwat
  REAL(KIND=jprb) :: rtber
  REAL(KIND=jprb) :: rtbercu
  REAL(KIND=jprb) :: rtice
  REAL(KIND=jprb) :: rticecu
  REAL(KIND=jprb) :: rtwat_rtice_r
  REAL(KIND=jprb) :: rtwat_rticecu_r
  
  ! LEPCLD : LOGICAL : TURN THE PROGNOSTIC CLOUD SCHEME ON
  LOGICAL lepcld
  
  !*     *YOEPHLI* CONTAINS CONSTANTS NEEDED BY
  !     THE LINEARIZED PHYSICS
  
  
  !     J.F. MAHFOUF        E.C.M.W.F.    23/06/96
  
  
  !     NAME        TYPE     DESCRIPTION
  !     ----        ----     -----------
  
  !     *RLPTRC*    REAL     CRITICAL TEMPERATURE FOR MIXED PHASE PROPERTIES
  !                          OF WATER
  !     *RLPAL1*    REAL     SMOOTHING COEFFICIENT
  !     *RLPAL2*    REAL     SMOOTHING COEFFICIENT
  !     *RLPBB*     REAL     CONSTANT FROM THE LOUIS ET AL. FORMULATION
  !     *RLPCC*     REAL     CONSTANT FROM THE LOUIS ET AL. FORMULATION
  !     *RLPDD*     REAL     CONSTANT FROM THE LOUIS ET AL. FORMULATION
  !     *RLPMIXL*   REAL     PSEUDO DEPTH OF THE PLANETARY BOUNDARY LAYER
  !     *RLPBETA*   REAL     REDUCTION FACTOR OF THE ASYMPTOTIC MIXING LENGTH
  !     *RLPDRAG*   REAL     COEFFICIENT FOR THE ESTIMATION OF SURFACE DRAG
  !     *RLPEVAP*   REAL     FRACTION OF POSSIBLE RAINFALL EVAPORATION
  !     *RLPP00*    REAL     PRESSURE ABOVE WHICH RADIATION IS NOT APPLIED
  !     *LPHYLIN*   LOGICAL  TRUE WHEN LINEARIZED PHYSICS IS ACTIVATED
  !     *LENOPERT   LOGICAL  TRUE WHEN NO PERTURBATION IS REQUIRED
  !                          FOR SURFACE ARRAYS
  !     *LRAISANEN   LOGICAL  TRUE WHEN RAISANEN OVERLAP SCHEME IS
  !                           ACTIVATED
  !     ------------------------------------------------------------------
  LOGICAL lphylin
  LOGICAL lenopert
  LOGICAL lraisanen
  
  REAL(KIND=jprb) :: rlptrc
  REAL(KIND=jprb) :: rlpal1
  REAL(KIND=jprb) :: rlpal2
  REAL(KIND=jprb) :: rlpbb
  REAL(KIND=jprb) :: rlpcc
  REAL(KIND=jprb) :: rlpdd
  REAL(KIND=jprb) :: rlpmixl
  REAL(KIND=jprb) :: rlpbeta
  REAL(KIND=jprb) :: rlpdrag
  REAL(KIND=jprb) :: rlpevap
  REAL(KIND=jprb) :: rlpp00
  
  !*    ** *YOECUMF* - PARAMETERS FOR CUMULUS MASSFLUX SCHEME
  
  !     M.TIEDTKE       E. C. M. W. F.      18/1/89

  !     NAME      TYPE      PURPOSE
  !     ----      ----      -------
  
  !     LMFPEN    LOGICAL  TRUE IF PENETRATIVE CONVECTION IS SWITCHED ON
  !     LMFSCV    LOGICAL  TRUE IF SHALLOW     CONVECTION IS SWITCHED ON
  !     LMFMID    LOGICAL  TRUE IF MIDLEVEL    CONVECTION IS SWITCHED ON
  !     LMFDD     LOGICAL  TRUE IF CUMULUS DOWNDRAFT      IS SWITCHED ON
  !     LMFIT     LOGICAL  TRUE IF UPDRAUGHT ITERATION
  !     LMFDUDV   LOGICAL  TRUE IF CUMULUS FRICTION       IS SWITCHED ON
  !     LMFSMOOTH LOGICAL  TRUE IF MASS FLUXES TOP/BOTTOM SMOOTHED FOR TRACER TRANSPORT
  !     LMFTRAC   LOGICAL  TRUE IF CONVECTIVE TRACER TRANSPORT IS SWITCHED ON
  !     LMFWSTAR  LOGICAL  TRUE IF GRANT W* CLOSURE USED FOR SHALLOW
  !     ENTRORG   REAL     COEFFICIENT FOR ORGANIZED ENTRAINMENT DEEP CONVECTION
  !     ENTRDD    REAL     ENTRAINMENT RATE FOR CUMULUS DOWNDRAFTS
  !     DETRPEN   REAL     DETRAINMENT RATE FOR PENETRATIVE CONVECTION
  !     RMFCMAX   REAL     MAXIMUM MASSFLUX VALUE ALLOWED FOR
  !     RMFCMIN   REAL     MINIMUM MASSFLUX VALUE (FOR SAFETY)
  !     RMFDEPS   REAL     FRACTIONAL MASSFLUX FOR DOWNDRAFTS AT LFS
  !     RDEPTHS   REAL     MAXIMUM ALLOWED CLOUD THICKNESS FOR SHALLOW
  !     RPRCON    REAL     COEFFICIENTS FOR DETERMINING CONVERSION
  !                        FROM CLOUD WATER TO RAIN
  !     RCPECONS  REAL     COEFFICIENT FOR RAIN EVAPORATION BELOW CLOUD
  !     RCUCOV    REAL     CONVECTIVE CLOUD COVER FOR RAIN EVPORATION
  !     RTAUMEL   REAL     TIME CONSTANT FOR MELTING
  !     RHEBC     REAL     REL. HUMIDITY BELOW CLOUD FOR WHICH EVAPORATION STARTS
  !     RTAU      REAL     ADJUSTMENT TIME SCALE IN CAPE CLOSURE
  !     ICAPDCYCL INTEGER  DIURNAL CYCLE MODULATION FOR CAPE
  !                        0:NO  1:SURFACE SENS HEAT FLUX 2: SUBCLOUD CAPE
  !     RTAU0     REAL     PRORTIONALITY FACTOR USED FOR RCAPDCYCL TIME SCALE
  !     RMFCFL    REAL     MULTIPLE OF CFL STABILITY CRITERIUM
  !     RMFLIC    REAL     USE CFL (1) OR ABSOLUT MASS FLUX LIMIT (0)
  !     RMFLIA    REAL     VALUE OF ABSOLUT MASS FLUX LIMIT
  !     RMFSOLCT  REAL     SOLVER FOR MASSFLUX ADVECTION EQUATION FOR TRACERS
  !                        0 : EXPLICIT  0-1 : SEMI-IMPLICIT >=1 : FULLY IMPLICIT
  !     RMFSOLUV  REAL     SOLVER FOR MASSFLUX ADVECTION EQUATION FOR MOMENTUM
  !     RMFSOLTQ  REAL     SOLVER FOR MASSFLUX ADVECTION EQUATION FOR T AND Q
  !     RUVPER    REAL     UPDRAIGHT VELOCITY PERTURBATION AT KLEV FOR IMPLICIT
  !     NJKT1, NJKT2, NJKT3-5 INTEGER  LEVEL LIMITS FOR CUBASEN/CUDDR
  !     ----------------------------------------------------------------
  
  ! REAL(KIND=jprb) :: entrorg
  REAL(KIND=jprb) :: entshalp
  REAL(KIND=jprb) :: entstpc1
  REAL(KIND=jprb) :: entstpc2
  REAL(KIND=jprb) :: entrdd
  REAL(KIND=jprb) :: detrpen
  ! REAL(KIND=jprb) :: rmfcfl -> moved into phy_params because it is resolution-dependent
  REAL(KIND=jprb) :: rmflic
  REAL(KIND=jprb) :: rmflia
  REAL(KIND=jprb) :: rmflmax
  REAL(KIND=jprb) :: rmfsoluv
  REAL(KIND=jprb) :: rmfsoltq
  REAL(KIND=jprb) :: rmfsolct
  REAL(KIND=jprb) :: rmfcmax
  REAL(KIND=jprb) :: rmfcmin
  REAL(KIND=jprb) :: rmfdeps
  REAL(KIND=jprb) :: rdepths
  REAL(KIND=jprb) :: rprcon
  ! REAL(KIND=jprb) :: rtau -> moved into phy_params because it is resolution-dependent
  ! REAL(KIND=jprb) :: rtau0 -> moved into phy_params because it is resolution-dependent
  INTEGER         :: icapdcycl
  REAL(KIND=jprb) :: rcpecons
  ! REAL(KIND=jprb) :: rcucov
  REAL(KIND=jprb) :: rtaumel
  ! REAL(KIND=jprb) :: rhebc
  REAL(KIND=jprb) :: ruvper

  LOGICAL :: lmfdd
  LOGICAL :: lmfit
  LOGICAL :: lmfdudv
  LOGICAL :: LMFUVDIS
  LOGICAL :: lmfsmooth
  LOGICAL :: lmftrac
  LOGICAL :: lmfwstar
  ! INTEGER(KIND=jpim) :: njkt1, njkt2, njkt3, njkt4, njkt5

  !     -----------------------------------------------------------------
  !     ** YOECLDP - CONTROL PARAMETERS FOR PROGNOSTIC CLOUD SCHEME
  !     -----------------------------------------------------------------

  !     * E.C.M.W.F. PHYSICS PACKAGE *

  !     C. JAKOB        E.C.M.W.F.          94/02/07
  ! Code Description:
  !      NAME     TYPE      PURPOSE
  !      ----     ----      -------

  !     *RAMID*   REAL      BASE VALUE FOR CALCULATION OF RELATIVE
  !                         HUMIDITY THRESHOLD FOR ONSET OF STRATIFORM
  !                         CONDENSATION (TIEDTKE, 1993, EQUATION 24)
  !     *RCLDIFF* REAL      DIFFUSION-COEFFICIENT FOR EVAPORATION BY
  !                         TURBULENT MIXING (IBID., EQU. 30)
  !     *RCLCRIT* REAL      BASE VALUE OF CRITICAL CLOUD WATER CONTENT
  !                         FOR CONVERSION TO RAIN (SUNDQUIST, 1988)
  !     *RKCONV*  REAL      BASE VALUE FOR CONVERSION COEFFICIENT (IBID.)
  !     *RPRC1*   REAL      COALESCENCE CONSTANT (IBID.)
  !     *RPRC2*   REAL      BERGERON-FINDEISEN CONSTANT (IBID.)
  !     *RCLDMAX* REAL      MAXIMUM CLOUD WATER CONTENT
  !     *RPECONS* REAL      EVAPORATION CONSTANT AFTER KESSLER
  !                         (TIEDTKE, 1993, EQU.35)
  !     *RENTRTU* REAL      ENTRAINMENT PARAMETER (DEARDOFF, 1976)
  !     *RENTRRA* REAL      ENTRAINMENT PARAMETER FOR RADIATION EFFECT
  !     *RAMIN*   REAL      LIMIT FOR A
  !     *RLMIN*   REAL      LIMIT FOR L
  !     *RSATQ*   REAL      LIMIT FOR SATURATION CHECK
  !     *RASMICE*   REAL    COEFFICIENT FOR "SMALL ICE" CALCULATION
  !                         (MACFARQUHAR AND HEYMSFIELD, JAS, 1997)
  !     *RBSMICE*   REAL    COEFFICIENT FOR "SMALL ICE" CALCULATION
  !                         (MACFARQUHAR AND HEYMSFIELD, JAS, 1997)

  REAL(KIND=jprb) :: ramid
  REAL(KIND=jprb) :: rcldiff
  REAL(KIND=jprb) :: rclcrit
  REAL(KIND=jprb) :: rkconv
  REAL(KIND=jprb) :: rprc1
  REAL(KIND=jprb) :: rprc2
  REAL(KIND=jprb) :: rcldmax
  REAL(KIND=jprb) :: rpecons
  REAL(KIND=jprb) :: rentrtu
  REAL(KIND=jprb) :: rentrra
  REAL(KIND=jprb) :: ramin
  REAL(KIND=jprb) :: rlmin
  REAL(KIND=jprb) :: rsatq
  REAL(KIND=jprb) :: rasmice
  REAL(KIND=jprb) :: rbsmice

  ! yomhook
  LOGICAL:: LHOOK=.FALSE.  
  !
  
  !yoecld
  PUBLIC :: rlmin
  !yomcst
  PUBLIC :: retv     ,rlvtt    ,rlstt    ,rtt       ,&
          & rg       ,rcpd     ,rd       ,rv        ,&
          & rkappa   ,ratm     ,rpi      ,rlmlt     ,&
          & rcvd     ,rsigma
  !yoecumf
  PUBLIC :: entshalp ,entstpc1 ,entstpc2            ,&
          & rprcon   ,rmfcmax  ,rmfcmin   ,detrpen  ,&
          & lmfdd    ,lmfdudv  , rdepths            ,&
          & lmfit    ,rmflic                       ,&
          & rmflia   ,rmfsoluv ,rmflmax            ,&
          & ruvper   ,rmfsoltq ,rmfsolct ,&
          & lmfsmooth,lmfwstar ,LMFUVDIS ,lmftrac  ,&
          & entrdd   ,& ! njkt1                    ,&
        ! & njkt2    ,njkt3    ,njkt4    ,njkt5    ,&
          & rcpecons ,rtaumel  ,& ! rcucov, rhebc  ,&
          & rmfdeps, icapdcycl
  !yoephli
  PUBLIC :: lphylin  ,rlptrc   ,rlpal1   ,rlpal2
  !yoephy
  PUBLIC :: lepcld
  !yoevdf
  PUBLIC :: rkap     ,rvdifts  ,repdu2   ,repust   ,rz0ice  ,&
          & rnum     ,rnuh     ,rnuq     ,rparzi   ,lelwdd
  !yoethf
  PUBLIC :: r2es     ,r3les    ,r3ies    ,r4les    ,r4ies   ,&
          & r5les    ,r5ies    ,rvtmp2   ,rhoh2o   ,r5alvcp ,&
          & r5alscp  ,ralvdcp  ,ralsdcp  ,ralfdcp  ,&
          & rtwat    ,rtber    ,rtbercu  ,rtice    ,rticecu ,&
          & rtwat_rtice_r      ,rtwat_rticecu_r
  !yoevdfs
  PUBLIC :: jpritbl  ,ritbl    ,aritbl   ,rcdhalf  ,&
          & rcdhpi2  ,rcheta   ,rchetb   ,rchbb    ,&
          & rchbcd   ,rchbd    ,rchb23a  ,rchbbcd  ,&
          & rchba    ,rchbhdl  ,rimax    ,dritbl   ,dri26
      
  PUBLIC :: phihu    ,phimu    ,phims    ,phihs

  PUBLIC :: sucst
  PUBLIC :: sucumf
  PUBLIC :: su_yoethf
  PUBLIC :: sucldp
  PUBLIC :: suphli
  PUBLIC :: suvdf
  PUBLIC :: suvdfs
  PUBLIC :: dr_hook
  PUBLIC :: lhook
  PUBLIC :: vdiv, vexp, vrec, vlog
  
CONTAINS
  
  ! fcttrm.h
  !!     ABSOLUTE THERMODYNAMICAL FUNCTIONS .
  !
  !!     RLV : LATENT HEAT OF VAPOURISATION
  !!     RLS : LATENT HEAT OF SUBLIMATION
  !!     RLF : LATENT HEAT OF FUSION
  !!     ESW : SATURATION IN PRESENCE OF WATER
  !!     ESS : SATURATION IN PRESENCE OF ICE
  !!     ES  : SATURATION (IF T>RTT THEN WATER ; IF T<RTT THEN ICE)
  !!        INPUT (FOR ALL SIX FUNCTIONS) : PTARG = TEMPERATURE .
  
  ELEMENTAL FUNCTION rlv(ptarg)
    REAL(KIND=jprb)             :: rlv
    REAL(KIND=jprb), INTENT(in) :: ptarg
    rlv = rlvtt+(rcpv-rcw)*(ptarg-rtt)
  END FUNCTION rlv
  
  ELEMENTAL FUNCTION rls(ptarg)
    REAL(KIND=jprb)             :: rls
    REAL(KIND=jprb), INTENT(in) :: ptarg
    rls = rlstt+(rcpv-rcs)*(ptarg-rtt)
  END FUNCTION rls
  
  ELEMENTAL FUNCTION rlf(ptarg)
    REAL(KIND=jprb)             :: rlf
    REAL(KIND=jprb), INTENT(in) :: ptarg
    rlf = rls(ptarg)-rlv(ptarg)
  END FUNCTION rlf
  
  ELEMENTAL FUNCTION esw(ptarg)
    REAL(KIND=jprb)             :: esw
    REAL(KIND=jprb), INTENT(in) :: ptarg
    esw = EXP(ralpw-rbetw/ptarg-rgamw*LOG(ptarg))
  END FUNCTION esw
  
  ELEMENTAL FUNCTION ess(ptarg)
    REAL(KIND=jprb)             :: ess
    REAL(KIND=jprb), INTENT(in) :: ptarg
    ess = EXP(ralps-rbets/ptarg-rgams*LOG(ptarg))
  END FUNCTION ess
  
  ELEMENTAL FUNCTION es (ptarg)
    REAL(KIND=jprb)             :: es
    REAL(KIND=jprb), INTENT(in) :: ptarg
    es   = EXP(                                                           &
      & (ralpw+ralpd*MAX(0.0_JPRB,SIGN(1.0_JPRB,rtt-ptarg)))           &
      & -(rbetw+rbetd*MAX(0.0_JPRB,SIGN(1.0_JPRB,rtt-ptarg)))/ptarg    &
      & -(rgamw+rgamd*MAX(0.0_JPRB,SIGN(1.0_JPRB,rtt-ptarg)))*LOG(ptarg))
  END FUNCTION es
  
  !! Orbit of the earth
  
  ELEMENTAL FUNCTION rteta(ptime)
    REAL(KIND=jprb)             :: rteta
    REAL(KIND=jprb), INTENT(in) :: ptime
    rteta = ptime/(rday*365.25_JPRB)
  END FUNCTION rteta
  
  ELEMENTAL FUNCTION rel(pteta)
    REAL(KIND=jprb)             :: rel
    REAL(KIND=jprb), INTENT(in) :: pteta
    rel = 1.7535_JPRB+6.283076_JPRB*pteta
  END FUNCTION rel
  
  PURE FUNCTION rem(pteta)
    REAL(KIND=jprb)             :: rem
    REAL(KIND=jprb), INTENT(in) :: pteta
    rem = 6.240075_JPRB+6.283020_JPRB*pteta
  END FUNCTION rem
  
  ELEMENTAL FUNCTION rrs(pteta)
    REAL(KIND=jprb)             :: rrs
    REAL(KIND=jprb), INTENT(in) :: pteta
    rrs = rea*(1.0001_JPRB-0.0163_JPRB*SIN(rel(pteta))&
      & +0.0037_JPRB*COS(rel(pteta)))
  END FUNCTION rrs
  
  !! Relative movement Sun/Earth
  PURE FUNCTION rlls (pteta)
    REAL(KIND=jprb)             :: rlls
    REAL(KIND=jprb), INTENT(in) :: pteta
    rlls = 4.8951_JPRB+6.283076_JPRB*pteta
  END FUNCTION rlls
  
  ELEMENTAL FUNCTION rllls (pteta)
    REAL(KIND=jprb)             :: rllls
    REAL(KIND=jprb), INTENT(in) :: pteta
    rllls = 4.8952_JPRB+6.283320_JPRB*pteta-0.0075_JPRB*SIN(rel(pteta))&
      & -0.0326_JPRB*COS(rel(pteta))-0.0003_JPRB*SIN(2.0_JPRB*rel(pteta))  &
      & +0.0002_JPRB*COS(2.0_JPRB*rel(pteta))
  END FUNCTION rllls
  
  ELEMENTAL FUNCTION rds (pteta)
    REAL(KIND=jprb)             :: rds
    REAL(KIND=jprb), INTENT(in) :: pteta
    rds = ASIN(SIN(repsm)*SIN(rllls(pteta)))
  END FUNCTION rds
  
  ELEMENTAL FUNCTION ret (pteta)
    REAL(KIND=jprb)             :: ret
    REAL(KIND=jprb), INTENT(in) :: pteta
    ret = 591.8_JPRB*SIN(2.0_JPRB*rlls(pteta))-459.4_JPRB*SIN(rem(pteta))&
      & +39.5_JPRB*SIN(rem(pteta))*COS(2.0_JPRB*rlls(pteta))            &
      & -12.7_JPRB*SIN(4._jprb*rlls(pteta))-4.8_JPRB*SIN(2.0_JPRB*rem(pteta))
  END FUNCTION ret
  !    -------------------------------------------------------------
  
  PURE FUNCTION ndd(kgrdat)
    INTEGER(KIND=jpim)             :: ndd
    INTEGER(KIND=jpim), INTENT(in) :: kgrdat
    ndd = MOD(kgrdat,100)
  END FUNCTION ndd
  
  ELEMENTAL FUNCTION nmm(kgrdat)
    INTEGER(KIND=jpim)             :: nmm
    INTEGER(KIND=jpim), INTENT(in) :: kgrdat
    nmm  =MOD((kgrdat-ndd(kgrdat))/100,100)
  END FUNCTION nmm
  
  PURE FUNCTION nccaa(kgrdat)
    INTEGER(KIND=jpim)             :: nccaa
    INTEGER(KIND=jpim), INTENT(in) :: kgrdat
    nccaa = kgrdat/10000
  END FUNCTION nccaa
  
  ELEMENTAL FUNCTION naa(kgrdat)
    INTEGER(KIND=jpim)             :: naa
    INTEGER(KIND=jpim), INTENT(in) :: kgrdat
    naa = MOD(nccaa(kgrdat),100)
  END FUNCTION naa
  
  ELEMENTAL FUNCTION namd(kgrdat)
    INTEGER(KIND=jpim)              :: namd
    INTEGER(KIND=jpim), INTENT(in)  :: kgrdat
    namd = MOD(kgrdat,1000000)
  END FUNCTION namd
  
  ELEMENTAL FUNCTION ncth(ksec)
    INTEGER(KIND=jpim)              :: ncth
    INTEGER(KIND=jpim), INTENT(in)  :: ksec
    ncth = ksec/3600
  END FUNCTION ncth
  
  ELEMENTAL FUNCTION ncent(kgrdat)
    INTEGER(KIND=jpim)              :: ncent
    INTEGER(KIND=jpim), INTENT(in)  :: kgrdat
    ncent = nccaa(kgrdat)/100+MIN(naa(kgrdat),1)
  END FUNCTION ncent
  
  ELEMENTAL FUNCTION nyearc(kgrdat)
    INTEGER(KIND=jpim)             :: nyearc
    INTEGER(KIND=jpim), INTENT(in) :: kgrdat
    nyearc = naa(kgrdat)+100*(1-MIN(naa(kgrdat),1))
  END FUNCTION nyearc
  
  ELEMENTAL FUNCTION nconstruct_date(kcent,kyearc,kmonth,kday)
    INTEGER(KIND=jpim)             :: nconstruct_date
    INTEGER(KIND=jpim) , INTENT(in):: kcent,kyearc,kmonth,kday
    nconstruct_date = (kcent-1)*10**6+kyearc*10**4+kmonth*10**2+kday
  END FUNCTION nconstruct_date
  
  ELEMENTAL FUNCTION nzzaa(kaaaa,kmm)
    INTEGER(KIND=jpim)             :: nzzaa
    INTEGER(KIND=jpim), INTENT(in) :: kaaaa,kmm
    nzzaa = kaaaa-( (1-SIGN(1,kmm-3))/2 )
  END FUNCTION nzzaa
  
  ELEMENTAL FUNCTION nzzmm(kmm)
    INTEGER(KIND=jpim)             :: nzzmm
    INTEGER(KIND=jpim), INTENT(in) :: kmm
    nzzmm = kmm+6*(1-SIGN(1,kmm-3))
  END FUNCTION nzzmm
  
  ELEMENTAL FUNCTION rjudat(kaaaa,kmm,kdd)
    REAL(KIND=jprb)                ::  rjudat
    INTEGER(KIND=jpim), INTENT(in) ::  kaaaa,kmm,kdd
    rjudat = 1720994.5_JPRB  &
      & + REAL(2-nzzaa(kaaaa,kmm)/100 + (nzzaa(kaaaa,kmm)/100)/4 &
      & + INT(365.25_JPRB*REAL(nzzaa(kaaaa,kmm),jprb))           &
      & + INT(30.601_JPRB*REAL(nzzmm(kmm)+1,jprb))               &
      & +  kdd,jprb)
  END FUNCTION rjudat
  
  ELEMENTAL FUNCTION rtime(kaaaa,kmm,kdd,kss)
    REAL(KIND=jprb)                ::  rtime
    INTEGER(KIND=jpim), INTENT(in) ::  kaaaa,kmm,kdd,kss
    rtime = (rjudat(kaaaa,kmm,kdd)-2451545._jprb)&
      & *  rday+REAL(kss,jprb)
  END FUNCTION rtime
  
  
  ! fcvdfs.h
  !     ------------------------------------------------------------------
  !     *FCVDFS** CONTAINS STATEMENT FUNCTIONS DESCRIBING STAB. FUNCT.
  !
  !          *THE STABILITY FUNCTIONS ARE THE SO-CALLED *PHI* AND
  !     *PSI*-FUNCTIONS. THE *PSI*-FUNCTIONS GIVE THE STABILITY
  !     CORRECTIONS IN THE LOGARITHMIC PROFILES FOR
  !     WIND, DRY STATIC ENERGY AND SPECIFIC HUMIDITY. THE FUNCTIONS
  !     DEPEND ON THE RATIO OF HEIGHT AND *OBUKHOV LENGTH (*ETA*).
  !          FOR THE UNSTABLE BOUNDARY LAYER, THE *DYER AND *HICKS
  !     FORMULATIONS ARE USED (CF. *DYER, 1974; *HOGSTROM, 1988). IN
  !     STABLE SITUATIONS, THE EMPIRICAL FORMS, PROPOSED BY *HOLTSLAG
  !     AND *DEBRUIN ARE USED WITH A MODIFICATION TO SATISFY A CRITICAL
  !     FLUX-*RICHARDSON NUMBER FOR LARGE *ETA*.
  !          THE *PHI* AND *PSI* FUNCTIONS ARE INTERRELATED. THE *PSI*
  !     FUNCTIONS CAN BE DERIVED FROM THE *PHI* FUNCTIONS BY INTEGRATION
  !     OF (1.-PHI)/ETA OR *PHI* FROM *PSI* BY COMPUTING
  !     (1.-ETA*DPSI/DETA) (SEE ALSO *HAUGEN, 1973; WORKSHOP ON
  !     MICROMETEOROLOGY, P. 77).
  !     ------------------------------------------------------------------


  !        *PHI AND *PSI FUNCTIONS FOR UNSTABLE SITUATIONS ACCORDING
  !        TO *DYER AND *HICKS

  ELEMENTAL FUNCTION phihu(peta)
    REAL(KIND=jprb)             :: phihu
    REAL(KIND=jprb), INTENT(in) :: peta
    phihu = 1.0_JPRB/     SQRT(1.0_JPRB-RCDHALF*PETA)
  END FUNCTION phihu

  ELEMENTAL FUNCTION phimu(peta)
    REAL(KIND=jprb)             :: phimu
    REAL(KIND=jprb), INTENT(in) :: peta
    phimu = 1.0_JPRB/SQRT(SQRT(1.0_JPRB-RCDHALF*PETA))
  END FUNCTION phimu

  !        *PHI AND *PSI FUNCTIONS FOR UNSTABLE SITUATIONS ACCORDING
  !        TO HOGSTROM FOR MOMENTUM AND DERIVED FROM THE ELLISON AND
  !        TURNER RELATION FOR THE RATIO OF PHIM AMD PHIH.
  
  ELEMENTAL FUNCTION phims(peta)
    REAL(KIND=jprb)             :: phims
    REAL(KIND=jprb), INTENT(in) :: peta
    phims = 1.0_JPRB+rcheta*peta
  END FUNCTION phims
  
  ELEMENTAL FUNCTION phihs(peta)
    REAL(KIND=jprb)             :: phihs
    REAL(KIND=jprb), INTENT(in) :: peta
    phihs = (1.0_JPRB+rchetb*peta)**2
  END FUNCTION phihs


!------------------------------------------------------------------------------
 
 
  SUBROUTINE sucst(kulout,kdat,ksss,kprintlev)
    !>
    !! Description:
    !!**** *SUCST * - Routine to initialize the constants of the model.
    
    !!     Purpose.
    !!     --------
    !!           Initialize and print the common YOMCST + initialize
    !!         date and time of YOMRIP.
    
    !!        Explicit arguments :
    !!        --------------------
    
    !!        KULOUT  - logical unit for the output
    !!        KDAT    - date in the form AAAAMMDD
    !!        KSSS    - number of seconds in the day
    !!        KPRINTLEV - printing level
    
    !!     Reference.
    !!     ----------
    !!        ECMWF Research Department documentation of the IFS
    
    !!     Author.
    !!     -------
    !!        Mats Hamrud and Philippe Courtier  *ECMWF*
    
    !!     Modifications.
    !!     --------------
    !!        Original : 87-10-15
    !!        Additions : 90-07-30 (J.-F. Geleyn)
    !!                    91-11-15 (M. Deque)
    !!                    96-08-12 M.Hamrud - Reduce printing
    !!    ------------------------------------------------------------------
    !!
    
    !USE PARKIND1  ,ONLY : JPIM     ,JPRB
    
    !USE YOMCST   , ONLY : RPI      ,RCLUM    ,RHPLA    ,RKBOL    ,&
    !            &RNAVO    ,RDAY     ,REA      ,REPSM    ,RSIYEA   ,&
    !            &RSIDAY   ,ROMEGA   ,RA       ,RG       ,R1SA     ,&
    !            &RSIGMA   ,RI0      ,R        ,RMD      ,RMV      ,&
    !            &RMO3     ,RD       ,RV       ,RCPD     ,RCPV     ,&
    !            &RCVD     ,RCVV     ,RKAPPA   ,RETV     ,RCW      ,&
    !            &RCS      ,RLVTT    ,RLSTT    ,RLVZER   ,RLSZER   ,&
    !            &RLMLT    ,RTT      ,RATM     ,RDT      ,RESTT    ,&
    !            &RALPW    ,RBETW    ,RGAMW    ,RALPS    ,RBETS    ,&
    !            &RGAMS    ,RALPD    ,RBETD    ,RGAMD
    !USE YOMRIP   , ONLY : RTIMST   ,RTIMTR
    
    !KF 'USE' instead of 'include'
    !#include "fctast.h"
    !#include "fcttrm.h"
    !#include "fcttim.h"
    
    
    !     DUMMY INTEGER SCALARS
    INTEGER(KIND=jpim) :: kdat
    INTEGER(KIND=jpim) :: kprintlev
    INTEGER(KIND=jpim) :: ksss
    INTEGER(KIND=jpim) :: kulout
    
    
    !     LOCAL INTEGER SCALARS
    INTEGER(KIND=jpim) :: ia, id, idat, im, isss, j
    
    !     LOCAL REAL SCALARS
    REAL(KIND=jprb) :: zde, zet, zju, zrs, zrsrel, zteta, zti
    
    !      -----------------------------------------------------------------
    
    !*       1.    DEFINE FUNDAMENTAL CONSTANTS.
    !              -----------------------------
    
    rpi=2._jprb*ASIN(1._jprb)
    rclum=299792458._jprb
    rhpla=6.6260755E-34_JPRB
    rkbol=1.380658E-23_JPRB
    rnavo=6.0221367E+23_JPRB
    
    !     ------------------------------------------------------------------
    
    !*       2.    DEFINE ASTRONOMICAL CONSTANTS.
    !              ------------------------------
#ifdef __GME__    
    rday=86400._jprb
#endif
    rea=149597870000._jprb
    repsm=0.409093_JPRB
    
    rsiyea=365.25_JPRB*rday*2._jprb*rpi/6.283076_JPRB
    rsiday=rday/(1._jprb+rday/rsiyea)
    romega=2._jprb*rpi/rsiday
    
    idat=kdat
    isss=ksss
    id=ndd(idat)
    im=nmm(idat)
    ia=nccaa(idat)
    zju=rjudat(ia,im,id)
    zti=rtime(ia,im,id,isss)
    !RTIMST=ZTI
    !RTIMTR=ZTI
    zteta=rteta(zti)
    zrs=rrs(zteta)
    zde=rds(zteta)
    zet=ret(zteta)
    zrsrel=zrs/rea
    
    !     ------------------------------------------------------------------
    
    !*       3.    DEFINE GEOIDE.
    !              --------------
    
    rg=9.80665_JPRB
    ra=6371229._jprb
    r1sa=REAL(1._jprb/REAL(ra,KIND(1._jprb)),KIND(r1sa))
    
    !     ------------------------------------------------------------------
    
    !*       4.    DEFINE RADIATION CONSTANTS.
    !              ---------------------------
    
    rsigma=2._jprb * rpi**5 * rkbol**4 /(15._jprb* rclum**2 * rhpla**3)
    ri0=1370._jprb
    
    !     ------------------------------------------------------------------
    
    !*       5.    DEFINE THERMODYNAMIC CONSTANTS, GAS PHASE.
    !              ------------------------------------------
    
    r=rnavo*rkbol
    rmd=28.9644_JPRB
    rmv=18.0153_JPRB
    rmo3=47.9942_JPRB
    rd=1000._jprb*r/rmd
    rv=1000._jprb*r/rmv
    rcpd=3.5_JPRB*rd
    rcvd=rcpd-rd
    rcpv=4._jprb *rv
    rcvv=rcpv-rv
    rkappa=rd/rcpd
    retv=rv/rd-1._jprb
    
    !     ------------------------------------------------------------------
    
    !*       6.    DEFINE THERMODYNAMIC CONSTANTS, LIQUID PHASE.
    !              ---------------------------------------------
    
    rcw=4218._jprb
    
    !     ------------------------------------------------------------------
    
    !*       7.    DEFINE THERMODYNAMIC CONSTANTS, SOLID PHASE.
    !              --------------------------------------------
    
    rcs=2106._jprb
    
    !     ------------------------------------------------------------------
    
    !*       8.    DEFINE THERMODYNAMIC CONSTANTS, TRANSITION OF PHASE.
    !              ----------------------------------------------------
    
    rtt=273.16_JPRB
    rdt=11.82_JPRB
    rlvtt=2.5008E+6_JPRB
    rlstt=2.8345E+6_JPRB
    rlvzer=rlvtt+rtt*(rcw-rcpv)
    rlszer=rlstt+rtt*(rcs-rcpv)
    rlmlt=rlstt-rlvtt
    ratm=100000._jprb
    
    !     ------------------------------------------------------------------
    
    !*       9.    SATURATED VAPOUR PRESSURE.
    !              --------------------------
    
    restt=611.14_JPRB
    rgamw=(rcw-rcpv)/rv
    rbetw=rlvtt/rv+rgamw*rtt
    ralpw=LOG(restt)+rbetw/rtt+rgamw*LOG(rtt)
    rgams=(rcs-rcpv)/rv
    rbets=rlstt/rv+rgams*rtt
    ralps=LOG(restt)+rbets/rtt+rgams*LOG(rtt)
    rgamd=rgams-rgamw
    rbetd=rbets-rbetw
    ralpd=ralps-ralpw
    
    !     ------------------------------------------------------------------
    
    !*      10.    PRINTS
    
    IF (kprintlev >= 1) THEN
      WRITE(kulout,'(''0*** Constants of the ICM   ***'')')
      WRITE(kulout,'('' *** Fundamental constants ***'')')
      WRITE(kulout,'(''           PI = '',E13.7,'' -'')')rpi
      WRITE(kulout,'(''            c = '',E13.7,''m s-1'')')rclum
      WRITE(kulout,'(''            h = '',E13.7,''J s'')')rhpla
      WRITE(kulout,'(''            K = '',E13.7,''J K-1'')')rkbol
      WRITE(kulout,'(''            N = '',E13.7,''mol-1'')')rnavo
      WRITE(kulout,'('' *** Astronomical constants ***'')')
      WRITE(kulout,'(''          day = '',E13.7,'' s'')')rday
      WRITE(kulout,'('' half g. axis = '',E13.7,'' m'')')rea
      WRITE(kulout,'('' mean anomaly = '',E13.7,'' -'')')repsm
      WRITE(kulout,'('' sideral year = '',E13.7,'' s'')')rsiyea
      WRITE(kulout,'(''  sideral day = '',E13.7,'' s'')')rsiday
      WRITE(kulout,'(''        omega = '',E13.7,'' s-1'')')romega
      
      WRITE(kulout,'('' The initial date of the run is :'')')
      WRITE(kulout,'(1X,I8,1X,I5,5X,I4,1X,I2,1X,I2)')idat,isss,ia,im,id
      WRITE(kulout,'('' The Julian date is : '',F11.2)') zju
      WRITE(kulout,'('' Time of the model  : '',F15.2,'' s'')')zti
      WRITE(kulout,'('' Distance Earth-Sun : '',E13.7,'' m'')')zrs
      WRITE(kulout,'('' Relative Dist. E-S : '',E13.7,'' m'')')zrsrel
      WRITE(kulout,'('' Declination        : '',F12.5)') zde
      WRITE(kulout,'('' Eq. of time        : '',F12.5,'' s'')')zet
      WRITE(kulout,'('' ***         Geoide         ***'')')
      WRITE(kulout,'(''      Gravity = '',E13.7,'' m s-2'')')rg
      WRITE(kulout,'('' Earth radius = '',E13.7,'' m'')')ra
      WRITE(kulout,'('' Inverse E.R. = '',E13.7,'' m'')')r1sa
      WRITE(kulout,'('' ***        Radiation       ***'')')
      WRITE(kulout,'('' Stefan-Bol.  = '',E13.7,'' W m-2 K-4'')')  rsigma
      WRITE(kulout,'('' Solar const. = '',E13.7,'' W m-2'')')ri0
      WRITE(kulout,'('' *** Thermodynamic, gas     ***'')')
      WRITE(kulout,'('' Perfect gas  = '',e13.7)') r
      WRITE(kulout,'('' Dry air mass = '',e13.7)') rmd
      WRITE(kulout,'('' Vapour  mass = '',e13.7)') rmv
      WRITE(kulout,'('' Ozone   mass = '',e13.7)') rmo3
      WRITE(kulout,'('' Dry air cst. = '',e13.7)') rd
      WRITE(kulout,'('' Vapour  cst. = '',e13.7)') rv
      WRITE(kulout,'(''         Cpd  = '',e13.7)') rcpd
      WRITE(kulout,'(''         Cvd  = '',e13.7)') rcvd
      WRITE(kulout,'(''         Cpv  = '',e13.7)') rcpv
      WRITE(kulout,'(''         Cvv  = '',e13.7)') rcvv
      WRITE(kulout,'(''      Rd/Cpd  = '',e13.7)') rkappa
      WRITE(kulout,'(''     Rv/Rd-1  = '',e13.7)') retv
      WRITE(kulout,'('' *** Thermodynamic, liquid  ***'')')
      WRITE(kulout,'(''         Cw   = '',E13.7)') rcw
      WRITE(kulout,'('' *** thermodynamic, solid   ***'')')
      WRITE(kulout,'(''         Cs   = '',E13.7)') rcs
      WRITE(kulout,'('' *** Thermodynamic, trans.  ***'')')
      WRITE(kulout,'('' Fusion point  = '',E13.7)') rtt
      WRITE(kulout,'('' RTT-Tx(ew-ei) = '',E13.7)') rdt
      WRITE(kulout,'(''        RLvTt  = '',E13.7)') rlvtt
      WRITE(kulout,'(''        RLsTt  = '',E13.7)') rlstt
      WRITE(kulout,'(''        RLv0   = '',E13.7)') rlvzer
      WRITE(kulout,'(''        RLs0   = '',E13.7)') rlszer
      WRITE(kulout,'(''        RLMlt  = '',E13.7)') rlmlt
      WRITE(kulout,'('' Normal press. = '',E13.7)') ratm
      WRITE(kulout,'('' Latent heat :  '')')
      WRITE(kulout,'(10(1X,E10.4))') (10._jprb*REAL(j,jprb),j=-4,4)
      WRITE(kulout,'(10(1X,E10.4))') (rlv(rtt+10._jprb*REAL(j,jprb)),j=-4,4)
      WRITE(kulout,'(10(1X,E10.4))') (rls(rtt+10._jprb*REAL(j,jprb)),j=-4,4)
      WRITE(kulout,'('' *** Thermodynamic, satur.  ***'')')
      WRITE(kulout,'('' Fusion point = '',E13.7)') rtt
      WRITE(kulout,'(''      es(Tt)  = '',e13.7)') restt
      WRITE(kulout,'('' es(T) :  '')')
      WRITE(kulout,'(10(1X,E10.4))') (10._jprb*REAL(j,jprb),j=-4,4)
      WRITE(kulout,'(10(1X,E10.4))') (esw(rtt+10._jprb*REAL(j,jprb)),j=-4,4)
      WRITE(kulout,'(10(1X,E10.4))') (ess(rtt+10._jprb*REAL(j,jprb)),j=-4,4)
      WRITE(kulout,'(10(1X,E10.4))') (es (rtt+10._jprb*REAL(j,jprb)),j=-4,4)
    ENDIF

    !RETURN
  END SUBROUTINE sucst
  

!------------------------------------------------------------------------------

  
  SUBROUTINE sucumf(rsltn,klev,pmean,phy_params,lshallow_only)

!     THIS ROUTINE DEFINES DISPOSABLE PARAMETERS FOR MASSFLUX SCHEME

!          M.TIEDTKE         E.C.M.W.F.    2/89

!          INTERFACE
!          ---------

!          THIS ROUTINE IS CALLED FROM *INIPHY*

!          MODIFICATIONS
!          -------------
!          P. Bechtold 2003-2013       Cleaning and revision of entrainment rates
!                                      options implicit, tracers, perturb, stand atmos
!                                      adding scaling factors for different planet
!                                      (modified gravity)
!                                      add options for diurnal cycle over land
!          P. Lopez, ECMWF (Oct 2007)  Put reading of NAMCUMF back in.
!          R. Forbes, May 2008         Changed factor in RTAUMEL from 
!                                      1.5 to 0.66
!          N. Semane+P.Bechtold     04-10-2012 Add RCORIOI/RPLRG/RPLDARE/RHOUR/RCVRFACTOR for small planet
!          T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!-----------------------------------------------------------------------

!USE PARKIND1 , ONLY : JPIM, JPRB
!USE YOMHOOK  , ONLY : LHOOK, DR_HOOK
!USE YOMCST   , ONLY : RG
!USE YOMLUN   , ONLY : NULOUT, NULNAM
!USE YOMDIMV  , ONLY : YRDIMV
!USE YOECUMF  , ONLY : ENTRORG   ,ENTRDD   ,&
! & ENTSHALP ,DETRPEN  ,RMFCMAX  ,RMFCMIN  ,RMFDEPS  ,RDEPTHS  ,&
! & ENTSTPC1 ,ENTSTPC2 ,&
! & RPRCON   ,RTAU     ,RTAUA    ,RTAU0     ,RCAPDCYCL,&
! & RCUCOV   ,RCPECONS ,RTAUMEL  ,RHEBC    ,RCVRFACTOR,&
! & LMFPEN   ,LMFSCV   ,LMFMID   ,LMFSMOOTH,LMFWSTAR ,LMFUVDIS,&
! & LMFDD    ,LMFDUDV  ,LMFCUCA  ,LMFPROFP ,&
! & RMFSOLUV ,RMFSOLTQ ,RMFSOLCT,&
! & RMFCFL   ,RMFLIC   ,RMFLIA   ,RUVPER   ,RBASE0   ,RMINCIN ,&
! & NJKT1    ,NJKT2    ,NJKT3    ,NJKT4    ,NJKT5    ,NJKT6
!USE YOMSTA  , ONLY : YRSTA
!USE YOMDYNCORE, ONLY : RPLRADI, RPLRG, RPLDARE

IMPLICIT NONE
!INCLUDE "mpif.h"
!INCLUDE "gme_commpi.h"

! NFLEVG : number of levels in grid point space
INTEGER(KIND=jpim) :: nflevg
INTEGER(KIND=jpim), INTENT(in) :: klev
REAL(KIND=jprb)   , INTENT(in) :: rsltn
REAL(KIND=jprb)   , INTENT(in) :: pmean(klev)
TYPE(t_phy_params), INTENT(inout) :: phy_params
LOGICAL           , INTENT(in) :: lshallow_only
!* change to operations

#ifdef __GME__
INTEGER(KIND=jpim) :: nulout=6
#endif

INTEGER(KIND=jpim) :: jlev
!INTEGER(KIND=JPIM) :: myrank,ierr,size
REAL(KIND=jprb) :: zhook_handle, zres_thresh, zres_thresh_trop, zfac, ztrans_end
!-----------------------------------------------------------------------

IF (lhook) CALL dr_hook('SUCUMF',0,zhook_handle)

nflevg=klev

!     1.           SPECIFY PARAMETERS FOR MASSFLUX-SCHEME
!                  --------------------------------------

!     DETRPEN: AVERAGE DETRAINMENT RATE FOR PENETRATIVE CONVECTION (1/M)
!     -------

detrpen=0.75E-4_JPRB  

!         NOTA:SHALLOW/DEEP ENTRAINMENT RATES ARE 
!              VERTICALLY SCALED BY FUNCTION  (qs/qsb)**3

!     ENTRORG: ENTRAINMENT FOR POSITIVELY BUOYANT DEEP/SHALLOW CONVECTION 1/(M)
!     -------

!ENTRORG=1.75E-3_JPRB     !40r3 default
!ENTRORG=1.9E-3_JPRB      !value tuned for ICON before changing to resolution-dependent setting
! ** entrorg is now set via the tuning namelist and depends on model resolution (see below) **

!     ENTSHALP: SHALLOW ENTRAINMENT DEFINED AS ENTSHALP*ENTRORG
!     --------

ENTSHALP=2.0_JPRB

!     ENTSTPC1,2: SHALLOW ENTRAINMENT CONSTANTS FOR TRIGGER TEST PARCEL ONLY
!     ----------

ENTSTPC1=0.55_JPRB
ENTSTPC2=1.E-4_JPRB
!ENTSTPC1=0.8_JPRB        !40r3 default
!ENTSTPC2=2.E-4_JPRB      !40r3 default

!     ENTRDD: AVERAGE ENTRAINMENT RATE FOR DOWNDRAFTS
!     ------

entrdd =2.0E-4_JPRB
!entrdd =3.0E-4_JPRB      !40r3 default

!     RMFCMAX:   MAXIMUM MASSFLUX VALUE ALLOWED FOR UPDRAFTS ETC
!     -------

rmfcmax=1._jprb

!     RMFCMIN:   MINIMUM MASSFLUX VALUE (FOR SAFETY)
!     -------

rmfcmin=1.e-10_JPRB

!     RMFDEPS:   FRACTIONAL MASSFLUX FOR DOWNDRAFTS AT LFS
!     -------

RMFDEPS=0.30_JPRB

!     RDEPTHS:   MAXIMUM ALLOWED SHALLOW CLOUD DEPTH (Pa)
!     -------

rdepths=2.e4_jprb

!     RPRCON:    COEFFICIENTS FOR DETERMINING CONVERSION FROM CLOUD WATER
!     ------

rprcon =1.4E-3_JPRB

!                COEFFICIENTS FOR RAIN EVAPORATION BELOW CLOUD
!                AND MELTING
!                ---------------------------------------------
!     RCPECONS:  KESSLER COEFFICIENT
!     RCUCOV:    ASSUMED CONVECTIVE CLOUD COVER
!     RTAUMEL:   MELTING TIME SCALE
!     RHEBC:     CRITICAL RELATIVE HUMIDITY BELOW CLOUD  FOR EVAPORATION

rcpecons=5.44E-4_JPRB/rg
rtaumel=5._jprb*3.6E3_JPRB*1.5_JPRB

!
! resolution-dependent setting of rhebc for mesh sizes below the threshold given by zres_thresh
zres_thresh      = 20.0E3_JPRB   ! 20 km
zres_thresh_trop = 12.5E3_JPRB   ! 12.5 km for tropics
ztrans_end       = 1.0E3_JPRB    ! 1 km - end of transition range

phy_params%rhebc_land       = tune_rhebc_land
phy_params%rhebc_ocean      = tune_rhebc_ocean
phy_params%rcucov           = tune_rcucov
phy_params%rhebc_land_trop  = tune_rhebc_land_trop
phy_params%rhebc_ocean_trop = tune_rhebc_ocean_trop
phy_params%rcucov_trop      = tune_rcucov_trop

!
IF (rsltn < zres_thresh) THEN
  phy_params%rhebc_land  = tune_rhebc_land  + (1._JPRB-tune_rhebc_land )*LOG(zres_thresh/rsltn)/LOG(zres_thresh/ztrans_end)
  phy_params%rhebc_ocean = tune_rhebc_ocean + (1._JPRB-tune_rhebc_ocean)*LOG(zres_thresh/rsltn)/LOG(zres_thresh/ztrans_end)
  !
  phy_params%rcucov      = tune_rcucov      + (1._JPRB-tune_rcucov)*(LOG(zres_thresh/rsltn)/LOG(zres_thresh/ztrans_end))**2
  !
  ! no one should use the convection scheme at resolutions finer than ztrans_end, but to be safe...
  phy_params%rhebc_land  = MIN(1._JPRB, phy_params%rhebc_land)
  phy_params%rhebc_ocean = MIN(1._JPRB, phy_params%rhebc_ocean)
  phy_params%rcucov      = MIN(1._JPRB, phy_params%rcucov)
ENDIF

IF (rsltn < zres_thresh_trop) THEN
  phy_params%rhebc_land_trop  = tune_rhebc_land_trop  + (1._JPRB-tune_rhebc_land_trop )* &
                                LOG(zres_thresh_trop/rsltn)/LOG(zres_thresh_trop/ztrans_end)
  phy_params%rhebc_ocean_trop = tune_rhebc_ocean_trop + (1._JPRB-tune_rhebc_ocean_trop)* &
                                LOG(zres_thresh_trop/rsltn)/LOG(zres_thresh_trop/ztrans_end)
  !
  phy_params%rcucov_trop      = tune_rcucov_trop      + (1._JPRB-tune_rcucov_trop)* &
                                (LOG(zres_thresh_trop/rsltn)/LOG(zres_thresh_trop/ztrans_end))**2
  !
  ! no one should use the convection scheme at resolutions finer than ztrans_end, but to be safe...
  phy_params%rhebc_land_trop  = MIN(1._JPRB, phy_params%rhebc_land_trop)
  phy_params%rhebc_ocean_trop = MIN(1._JPRB, phy_params%rhebc_ocean_trop)
  phy_params%rcucov_trop      = MIN(1._JPRB, phy_params%rcucov_trop)
ENDIF


! tuning parameter for organized entrainment of deep convection
phy_params%entrorg = tune_entrorg + 1.2E-4_JPRB*LOG(zres_thresh/rsltn)


! resolution-dependent settings for 'excess values' of temperature and QV used for convection triggering (test parcel ascent)

! This factor is 1 for dx = 20 km or coarser and 0 for dx = ztrans_end or finer
zfac = MIN(1._JPRB,LOG(MAX(1._JPRB,rsltn/ztrans_end))/LOG(zres_thresh/ztrans_end))

phy_params%texc = zfac*tune_texc   ! K
phy_params%qexc = zfac*tune_qexc   ! relative perturbation of grid-scale QV


!     SET ADJUSTMENT TIME SCALE FOR CAPE CLOSURE AS A FUNCTION
!     OF MODEL RESOLUTION

!     Cy32r1 and earlier:
!     convective adjustment time TAU=RTAU
!     RTAU IS 20 MINUTES FOR RESOLUTIONS HIGHER THAN TL319
!     RTAU IS 10 MINUTES FOR RESOLUTIONS HIGHER THAN TL511
!     RTAU IS 1 HOUR FOR ANY OTHER RESOLUTION
!     --------------------------------------------------------

!     from Cy32r3 onward:
!     CONVECTIVE ADJUSTMENT TIME TAU=Z_cld/W_cld*rtau
!     WHERE RTAU (unitless) NOW ONLY REPRESENTS THE RESOLUTION DEPENDENT PART

!phy_params%tau=1.0_JPRB+264.0_JPRB/REAL(ksmax,jprb)

! Basic resolution-dependent setting (tuned for ICON, denominator originally was 75 km)
phy_params%tau = 1.0_JPRB + rsltn/120.e3_jprb

! Set upper limit
phy_params%tau=MIN(3.0_JPRB,phy_params%tau)

! Increase adjustment time scale at resolutions below 10 km
IF (rsltn < 10.e3_jprb) phy_params%tau = phy_params%tau + (LOG(10.e3_jprb/rsltn))**2

! ** CAPE correction to improve diurnal cycle of convection ** (set now in mo_nwp_phy_nml)
! icapdcycl = 0! 0= no CAPE diurnal cycle correction (IFS default prior to cy40r1, i.e. 2013-11-19)
               ! 1=    CAPE - surface buoyancy flux (intermediate testing option)
               ! 2=    CAPE - subcloud CAPE (IFS default starting with cy40r1)
               ! 3=    Apply CAPE modification of (2) over land only, with additional restriction to the tropics

phy_params%tau0 = 1.0_jprb
IF (icapdcycl >= 2) phy_params%tau0 = 1.0_jprb/phy_params%tau

!     LOGICAL SWITCHES
!     ----------------

phy_params%lmfscv  =.TRUE.   ! shallow convection
IF (lshallow_only) THEN
  phy_params%lmfmid  =.FALSE.   ! mid-level convection
  phy_params%lmfpen  =.FALSE.   ! deep convection
ELSE
  phy_params%lmfmid  =.TRUE.   ! mid-level convection
  phy_params%lmfpen  =.TRUE.   ! deep convection
ENDIF

lmfdd   =.TRUE.   ! use downdrafts
lmfit   =.FALSE.  ! updraught iteration or not
LMFUVDIS=.TRUE.   ! use kinetic energy dissipation (addit T-tendency)
lmfdudv =.TRUE.   ! use convective momentum transport
!*UPG add to operations
lmftrac =.TRUE.   ! convective chemical tracer transport
lepcld  =.TRUE.   ! produce detrained cloud water/ice
! to reuse in prognostic cloud scheme

!     RMFCFL:     MASSFLUX MULTIPLE OF CFL STABILITY CRITERIUM
!     -------
phy_params%mfcfl = 2._JPRB*MIN(2._JPRB,1._JPRB + 2.5e-5_JPRB*rsltn)

rmflic=1.0_JPRB   ! use CFL mass flux limit (1) or absolut limit (0)
rmflia=0.0_JPRB   ! value of absolut mass flux limit
rmflmax=1.75_jprb ! mass flux limit following a suggestion by P. Bechtold [kg/(m**2s)]


!     MASSFLUX SOLVERs FOR MOMEMTUM AND TRACERS
!     0: EXPLICIT 0-1 SEMI-IMPLICIT >=1: IMPLICIT
!     -------------------------------------------

rmfsoluv=1.0_JPRB  ! mass flux solver for momentum
rmfsoltq=1.0_JPRB  ! mass flux solver for T and q
rmfsolct=1.0_JPRB  ! mass flux solver for chemical tracers
lmfsmooth=.FALSE.  ! Smoothing of mass fluxes top/bottom for Tracers
lmfwstar=.FALSE.   ! Grant w* closure for shallow convection

!     UPDRAUGHT VELOCITY PERTURBATION FOR IMPLICIT (M/S)
!     --------------------------------------------------

ruvper=0.3_JPRB

phy_params%kcon1=2
phy_params%kcon2=2
DO jlev=nflevg,2,-1
  ! IF(STPRE(JLEV) > 350.E2_JPRB)NJKT1=JLEV
  ! IF(STPRE(JLEV) >  60.E2_JPRB)NJKT2=JLEV
  ! IF(STPRE(JLEV) > 950.E2_JPRB)NJKT3=JLEV
  ! IF(STPRE(JLEV) > 850.E2_JPRB)NJKT4=JLEV
  ! IF(STPRE(JLEV) > 500.E2_JPRB)NJKT5=JLEV
  !  IF(PMEAN(JLEV)/PMEAN(KLEV)*1013.E2 > 350.E2_JPRB)NJKT1=JLEV
  !  IF(PMEAN(JLEV)/PMEAN(KLEV)*1013.E2 >  60.E2_JPRB)NJKT2=JLEV
  !  IF(PMEAN(JLEV)/PMEAN(KLEV)*1013.E2 > 950.E2_JPRB)NJKT3=JLEV
  !  IF(PMEAN(JLEV)/PMEAN(KLEV)*1013.E2 > 850.E2_JPRB)NJKT4=JLEV
  !  IF(PMEAN(JLEV)/PMEAN(KLEV)*1013.E2 > 500.E2_JPRB)NJKT5=JLEV
  IF(pmean(jlev) > 350.e2_jprb) phy_params%kcon1=jlev
  IF(pmean(jlev) >  60.e2_jprb) phy_params%kcon2=jlev
ENDDO

#ifdef __GME__
WRITE(6,*)'SUCUMF: NJKT1=',njkt1,' NJKT2=',njkt2,' NJKT3=',njkt3,' RESOLUTION=',rsltn
!WRITE(6,*)'SUCUMF: KSMAX=',KSMAX
WRITE(UNIT=nulout,FMT='('' COMMON YOECUMF '')')
WRITE(UNIT=nulout,FMT='('' LMFMID = '',L5 &
  & ,'' LMFDD = '',L5,'' LMFDUDV = '',L5 &
  & ,'' RTAU = '',E12.5,'' s-1'')') &
  & phy_params%lmfmid,lmfdd,lmfdudv,rtau
#endif

#ifdef __ICON__
CALL message('mo_cuparameters, sucumf', 'NJKT1, NJKT2, KSMAX')
WRITE(message_text,'(2i7,E12.5)') phy_params%kcon1, phy_params%kcon2, rsltn 
CALL message('mo_cuparameters, sucumf ', TRIM(message_text))
CALL message('mo_cuparameters, sucumf', 'LMFMID, LMFDD, LMFDUDV, RTAU, ENTRORG, TEXC, QEXC')
WRITE(message_text,'(4x,l6,l6,l6,F8.4,E11.4,2F8.5)')phy_params%lmfmid,lmfdd,lmfdudv,phy_params%tau,&
  phy_params%entrorg,phy_params%texc,phy_params%qexc
CALL message('mo_cuparameters, sucumf ', TRIM(message_text))
CALL message('mo_cuparameters, sucumf', 'RHEBC_LND, RHEBC_LND_TROP, RHEBC_OCE, RHEBC_OCE_TROP, RCUCOV, RCUCOV_TROP')
WRITE(message_text,'(4x,6F8.4)') phy_params%rhebc_land,phy_params%rhebc_land_trop,phy_params%rhebc_ocean, &
  phy_params%rhebc_ocean_trop,phy_params%rcucov,phy_params%rcucov_trop
CALL message('mo_cuparameters, sucumf ', TRIM(message_text))
#endif

IF (lhook) CALL dr_hook('SUCUMF',1,zhook_handle)

  END SUBROUTINE sucumf


!------------------------------------------------------------------------------


  SUBROUTINE su_yoethf
    !
    !USE PARKIND1  ,ONLY : JPIM     ,JPRB

    !USE YOETHF, ONLY : R2ES,    R3LES,  R3IES,  R4LES, R4IES,  &
    !&                   R5LES,   R5IES,  RVTMP2 ,RHOH2O,R5ALVCP,&
    !&                   R5ALSCP, RALVDCP,RALSDCP,RALFDCP,       &
    !&                   RTWAT,   RTBER,  RTBERCU,RTICE, RTICECU,&
    !&                   RTWAT_RTICE_R, RTWAT_RTICECU_R
    !USE YOMCST, ONLY : RD, RV, RCPD, RCPV, RATM, RTT, RLVTT, &
    !                   RLSTT, RLMLT

    IMPLICIT NONE

!   rvtmp2=rcpv/rcpd-1._jprb   !use cp,moist
    rvtmp2=0._jprb             !neglect cp,moist (as in IFS)
    rhoh2o=ratm/100._jprb
    !IFS original (Bolton(1980) for water and Buck(1981) for ice)
    r2es=611.21_JPRB*rd/rv
    r3les=17.502_JPRB
    r3ies=22.587_JPRB
    r4les=32.19_JPRB
    r4ies=-0.7_JPRB
    !GME values KF
    !R2ES=610.78_JPRB*RD/RV ! =B1
    !R3LES=17.269388_JPRB   ! =B2_w
    !R3IES=21.8745584_JPRB  ! =B2_i
    !R4LES= 35.86_JPRB      ! =B4_w
    !R4IES= 7.66_JPRB       ! =B4_i
    r5les=r3les*(rtt-r4les)
    r5ies=r3ies*(rtt-r4ies)
    r5alvcp=r5les*rlvtt/rcpd
    r5alscp=r5ies*rlstt/rcpd
    ralvdcp=rlvtt/rcpd
    ralsdcp=rlstt/rcpd
    ralfdcp=rlmlt/rcpd
    rtwat=rtt
    rtber=rtt-5._jprb
    rtbercu=rtt-5.0_JPRB
    rtice=rtt-23._jprb
    rticecu=rtt-23._jprb
    rtwat_rtice_r=1._jprb/(rtwat-rtice)
    rtwat_rticecu_r=1._jprb/(rtwat-rticecu)

  END SUBROUTINE su_yoethf


!------------------------------------------------------------------------------


  SUBROUTINE sucldp
    !
    ! Description:
    !**** *SUCLDP*   - INITIALIZE COMMON YOECLD CONTROLLING *CLOUDSC*

    !     PURPOSE.
    !     --------
    !           INITIALIZE YOECLDP

    !**   INTERFACE.
    !     ----------
    !        CALL *SUCLDP* FROM *SUPHEC*
    !              ------        ------

    !        EXPLICIT ARGUMENTS :
    !        --------------------
    !        NONE

    !        IMPLICIT ARGUMENTS :
    !        --------------------
    !        COMMON YOECLDP

    !     METHOD.
    !     -------
    !        SEE DOCUMENTATION

    !     EXTERNALS.
    !     ----------
    !        NONE

    !     REFERENCE.
    !     ----------
    !        ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE
    !     "INTEGRATED FORECASTING SYSTEM"

    !     AUTHOR.
    !     -------
    !        C.JAKOB   *ECMWF*

    !     MODIFICATIONS.
    !     --------------
    !        ORIGINAL : 94-02-07
    !     ------------------------------------------------------------------

    !
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Current Code Owner: DWD, Kristina Froehlich
    !    kristina.froehlich@dwd.de
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !
    ! History:
    ! Version    Date       Name
    ! ---------- ---------- ----
    ! @VERSION@  @DATE@     K. Froehlich
    !  Initial release
    !
    ! Code Description:
    ! Language: Fortran 90.
    ! Software Standards: "European Standards for Writing and
    ! Documenting Exchangeable Fortran 90 Code".
    !=======================================================================
    !

    !USE PARKIND1  ,ONLY : JPIM     ,JPRB
    !USE YOMCST   , ONLY : RG
    !USE YOECLDP  , ONLY : RAMID    ,RCLDIFF  ,RCLCRIT  ,RKCONV   ,&
    !            &RPRC1    ,RPRC2    ,RCLDMAX  ,RPECONS  ,RTAUMEL  ,&
    !            &RENTRTU  ,RENTRRA  ,RAMIN    ,RLMIN    ,RASMICE  ,&
    !            &RBSMICE, RSATQ


    !*       1.    SET VALUES
    !              ----------

    IMPLICIT NONE
    ramid=0.8_JPRB
    rcldiff=2.e-6_JPRB  !1e-6 pre 25r1 default
    rcldiff=3.e-6_jprb

    rclcrit=0.3E-3_JPRB
    rkconv=1.e-4_JPRB
    rprc1=100._jprb
    rprc2=0.5_JPRB
    rcldmax=5.e-3_JPRB

    rpecons=5.44E-4_JPRB/rg

    rentrtu=0.5_JPRB
    rentrra=0.5_JPRB

    ramin=1.e-8_JPRB
    rlmin=1.e-8_JPRB

    rasmice=0.252_JPRB
    rbsmice=0.837_JPRB

    rsatq=1.0E-7_JPRB

    !RETURN
  END SUBROUTINE sucldp


!------------------------------------------------------------------------------


  SUBROUTINE suphli
    !>
    !! Description:

    !!     ------------------------------------------------------------------
    !!**   *SUPHLI* IS THE SET-UP ROUTINE FOR COMMON BLOCK *YOEPHLI*

    !!     J.F. MAHFOUF         E.C.M.W.F.      96/06/23

    !!     PURPOSE
    !!     -------
    !!          THIS ROUTINE INITIALIZES THE CONSTANTS IN COMMON BLOCK
    !!     *YOEPHLI*

    !!     METHOD.
    !!     -------
    !!         INITIALIZATION OF THE CONSTANTS USED IN THE LINEARIZED
    !!         PHYSICS

    !!     EXTERNALS.
    !!     ----------

    !!        NONE


    !USE PARKIND1  ,ONLY : JPIM     ,JPRB

    !USE YOETHF   , ONLY : RTWAT    ,RTICE
    !USE YOEPHLI  , ONLY : LPHYLIN   ,LENOPERT ,LRAISANEN,&
    !            &RLPTRC   ,RLPAL1   ,RLPAL2   ,RLPBB    ,&
    !            &RLPCC    ,RLPDD    ,RLPMIXL  ,RLPBETA  ,&
    !            &RLPDRAG  ,RLPEVAP  ,RLPP00
    !USE YOERAD   , ONLY : LRRTM     , NOVLP
    !USE YOMLUN   , ONLY : NULNAM
    !USE YOMTLEVOL, ONLY : LTLEVOL

    IMPLICIT NONE


    !#include "namtlevol.h"

    !     ------------------------------------------------------------------

    !*         1.     SET LOGICAL TO SWICH ON LINEARIZED PHYSICS
    !                 ------------------------------------------

    lphylin = .FALSE.
    !LTLEVOL = .FALSE.

    !CALL POSNAM(NULNAM,'NAMTLEVOL')
    !READ(NULNAM,NAMTLEVOL)

    !*         1.1 No perturbation of surface arrays
    !          -------------------------------------

    lenopert = .TRUE.

    !*         2.     SET CONSTANTS RELATED TO WATER MIXED PHASE
    !                 ------------------------------------------


    rlptrc=rtice+(rtwat-rtice)/SQRT(2._jprb)
    rlpal1=0.15_JPRB
    rlpal2=20._jprb


    !*         3.     SET CONSTANTS RELATED TO VERTICAL DIFFUSION
    !                 -------------------------------------------


    !     CONSTANTS OF THE LOUIS FORMULATION

    rlpbb=5._jprb
    rlpcc=5._jprb
    rlpdd=5._jprb

    !     PSEUDO DEPTH OF THE BOUNDARY LAYER

    rlpmixl=4000._jprb

    !     REDUCTION FACTOR OF THE ASYMPTOTIC MIXING LENGTH

    rlpbeta=0.2_JPRB


    !*         4.     SET CONSTANTS RELATED TO GRAVITY WAVE DRAG
    !                 ------------------------------------------

    rlpdrag=0._jprb


    !*         5.     SET CONSTANTS RELATED TO RAINFALL EVAPORATION
    !                 ---------------------------------------------

    rlpevap=0._jprb


    !*         6.     SET CONSTANTS RELATED TO RADIATION
    !                 ----------------------------------

    !     Pressure level above which long-wave cooling is not applied

    rlpp00=30000._jprb

    !     Using Raisanen overlap scheme

    !*UPG Oper
    !IF (LRRTM .AND. (NOVLP.EQ.1)) THEN
    !  LRAISANEN = .TRUE.
    !ELSE
    lraisanen = .FALSE.
    !ENDIF


!    PRINT*, 'SUPHLI', rlptrc
    !RETURN
  END SUBROUTINE suphli


!------------------------------------------------------------------------------


  SUBROUTINE suvdf
    !>
    !! Description:
    !!     ------------------------------------------------------------------

    !!**   *SUVDF* IS THE SET-UP ROUTINE FOR COMMON BLOCK *YOEVDF*

    !!     A.C.M. BELJAARS         E.C.M.W.F.       2/11/89

    !!     PURPOSE
    !!     -------
    !!          THIS ROUTINE INITIALIZES THE CONSTANTS IN COMMON BLOCK
    !!     *YOEVDF*

    !!     EXTERNALS.
    !     ----------

    !!     NONE.

    !!     MODIFICATIONS
    !!     -------------
    !!     J.-J. MORCRETTE         E.C.M.W.F.      91/07/14
    !     ------------------------------------------------------------------

    !USE PARKIND1  ,ONLY : JPIM     ,JPRB

    !USE YOEVDF   , ONLY : RLAM     ,RKAP     ,RCHAR    ,RVDIFTS  ,&
    !            &RZ0ICE   ,REPDU2   ,REPUST   ,RSEZ0H   ,&
    !            &RSEZ0Q   ,RNUM     ,RNUH     ,RNUQ     ,RENTR    ,&
    !            &RPAR     ,RPAR1    ,RPARSRF  ,RPARZI   ,RLAMSK   ,&
    !            &LELWDD


    IMPLICIT NONE


    !     LOCAL REAL SCALARS
    REAL(KIND=jprb) :: cepz0o, znu


    !     ------------------------------------------------------------------

    !*         1.     SET FIRST SET OF CONSTANTS
    !                 --------------------------


    rlam   =150._jprb
    rkap   =0.4_JPRB
    rchar  =0.018_JPRB
    rvdifts=1.5_JPRB
    rz0ice =0.001_JPRB

    !     ------------------------------------------------------------------

    !*         2.      SET OTHER CONSTANTS
    !                  -------------------


    cepz0o=2._jprb
    repdu2 =(0.1_JPRB)**2
    repust=0.0001_JPRB
    rsez0h=1.4E-5_JPRB
    rsez0q=1.3E-4_JPRB

    !     KINEMATIC VISCOSITY OF AIR

    znu   =1.5E-5_JPRB
    rnum  =0.11_JPRB*znu
    rnuh  =0.40_JPRB*znu
    rnuq  =0.62_JPRB*znu

    !     ENTRAINMENT PARAMETRIZATION

    rentr=0.20_JPRB
    rpar=2._jprb
    rpar1=0.6_JPRB
    rparsrf=0.1_JPRB
    rparzi=1000._jprb

    !     COMPUTATION OF SKIN TEMPERATURE

    rlamsk=15._jprb
    lelwdd=.TRUE.

    !     ------------------------------------------------------------------

    !RETURN
  END SUBROUTINE suvdf


!------------------------------------------------------------------------------


  SUBROUTINE suvdfs
    !>
    !! Description:
    !!**   *SUBROUTINE* *SUVDFS* INITIALIZES COMMON BLOCK *YOEVDFS*

    !!     A. BELJAARS   E.C.M.W.F.   26/03/90


    !!      PURPOSE
    !!      -------

    !!           *SUBROUTINE *SUVDFS* INITIALIZES THE CONSTANTS NEEDED BY
    !!      THE EMPIRICAL STABILITY FUNCTIONS AND SETS UP THE TABLE THAT
    !!      GIVES THE STABILITY PARAMETER ETA (HEIGHT DEVIDED BY
    !!      *OBUKHOV LENGTH) AS A FUNCTION OF THE GRADIENT *RICHARDSON NUMBER
    !!      FOR STABLE SITUATIONS. ALSO THE SECOND DERIVATIVES ARE
    !!      TABULATED FOR LATER USE BY SPLINE INTERPOLATION.
    !!           *INIVDFS* CHECKS WETHER THE *PHI* AND *PSI* EXPRESSIONS ARE
    !!     CONSISTENT (BY MEANS OF NUMERICAL DIFFERENTIATION). IF THE
    !!      FUNCTIONS ARE NOT CONSISTENT, THE ROUTINE ABORTS.


    !!      METHOD
    !!      ------

    !!           *THE ALGEBRAIC EQUATION TO BE SOLVED IS
    !!      RI=(PHIH/PHIM**2)*ETA, WHERE *PHIH* AND *PHIM* ARE THE GRADIENT
    !!      STABILITY FUNCTIONS FOR HEAT AND MOMENTUM RESPECTIVELY.
    !!           *TO SOLVE THE IMPLICIT ALGEBRAIC EQUATION FOR *ETA* AS A
    !!      FUNCTION OF *RI*, *NEWTON'S METHOD IS EMPLOYED WITH A FIXED
    !!      NUMBER OF ITERATIONS. THE NECESSARY DERIVATIVES ARE EVALUATED
    !!     NUMERICALLY.
    !!           *AFTER COMPLETION OF THE FUNCTION TABLE, THE SECOND
    !!      DERIVATIVES ARE DERIVED FOR LATER USE BY THE SPLINE
    !!      INTERPOLATION.


    !!      EXTERNALS
    !!      ---------

    !!           *STATEMENT FUNCTIONS *PHIMS*, *PHIHS*, *PSIMS* AND *PSIHS*


    !!      REFERENCE
    !!      ---------

    !!           *SEE *PRESS ET AL. (1986; NUMERICAL RECIPES - THE ART OF
    !!      SCIENTIFIC COMPUTING) FOR DETAILS ON THE SPLINE INTERPOLATION.

    !      -----------------------------------------------------------------
    !
    !USE PARKIND1  ,ONLY : JPIM     ,JPRB

    !KF already defined in module
    !USE YOEVDFS  , ONLY : JPRITBL  ,RITBL    ,ARITBL   ,RCHBA    ,&
    !            &RCHBB    ,RCHBC    ,RCHBD    ,RCHB23A  ,RCHBBCD  ,&
    !            &RCHBCD   ,RCHETA   ,RCHETB   ,RCHETC   ,RCHBHDL  ,&
    !            &RCDHALF  ,RCDHPI2  ,RIMAX    ,DRITBL   ,DRI26

    !KF use instead of includes
    !USE  mo_cufunctions , ONLY : PHIHS,PHIMS

    IMPLICIT NONE

    REAL(KIND=jprb) :: zu(jpritbl)

    !     LOCAL INTEGER SCALARS
    INTEGER(KIND=jpim) :: itmax, jit, jjp

    !     LOCAL REAL SCALARS
    REAL(KIND=jprb) :: zdrv1, zdrvn, zeps, zeta, zetam, zetap, zfxm,&
      & zfxp, zp, zqn, zrib, zun


    !      INSERT STATEMENT FUNCTIONS (PSI-FUNCTIONS)

    !#include "fcvdfs.h"

    !*     1. INITIALIZE CONSTANTS IN COMMON BLOCK
    !         ---------- --------- -- ------ -----


    !     1.1 CONSTANTS RELATED TO ETA(RI)-TABLE

    rimax   = 10._jprb
    dritbl  = rimax/REAL(jpritbl-1,jprb)
    dri26   = dritbl**2/6._jprb
    ritbl(1)= 0._jprb

    !     1.2 CONSTANTS FOR THE ELLISON TURNER FUNCTIONS (STABLE)

    rcheta=5._jprb
    rchetb=4._jprb
    rchetc=8._jprb
    rchbhdl = 5._jprb

    !     1.21 CONSTANTS FOR HOLTSLAG AND DEBRUIN FUNCTIONS (STABLE PSI)

    rchba   = 1._jprb
    rchbb   = 2._jprb/3._jprb
    rchbc   = 5._jprb
    rchbd   = 0.35_JPRB
    rchb23a = (2._jprb/3._jprb)*rchba
    rchbbcd = rchbb*rchbc/rchbd
    rchbcd  = rchbc/rchbd

    !     1.3 CONSTANTS FOR DYER AND HICKS EXPRESSIONS (UNSTABLE)

    rcdhalf = 16._jprb
    rcdhpi2 = 2._jprb*ATAN(1._jprb)

    !*    3. LOOP OVER TABLE INDEX
    !        ---- ---- ----- -----

    DO jjp=2, jpritbl
      zrib=REAL(jjp-1,jprb)*dritbl

      !     3.1 INITIAL GUESS

      zeta=ritbl(jjp-1)
      IF (zrib  <  0.5_JPRB) THEN
        itmax=5
      ELSE
        itmax=3
      ENDIF

      !     3.2 NEWTON'S ITERATION LOOP WITH DERIVATIVE FROM FINITE DIFFERENCE

      DO jit=1, itmax
        zeps=(zeta+1._jprb)*0.001_JPRB
        zetap=zeta+zeps
        zetam=zeta-zeps
        zfxp=zrib-zetap*phihs(zetap)/phims(zetap)**2
        zfxm=zrib-zetam*phihs(zetam)/phims(zetam)**2
        zeta=zeta-zeps*(zfxp+zfxm)/(zfxp-zfxm)
      ENDDO

      !     3.3 STORE RESULT

      ritbl(jjp)=zeta

    ENDDO

    !     4. COMPUTE SECOND DERIVATIVES FROM TABULATED RESULTS
    !        ------- ------ ----------- ---- --------- -------


    !     4.1 DERIVATIVE AT RI=0. AND RI=RIMAX

    zdrv1=1._jprb
    zdrvn=(ritbl(jpritbl)-ritbl(jpritbl-1))/dritbl

    aritbl(1)=-0.5_JPRB
    zu(1)=(3._jprb/dritbl)*((ritbl(2)-ritbl(1))/dritbl-zdrv1)

    !     4.2. DECOMPOSITION LOOP OF TRIDIAGONAL MATRIX

    DO jjp=2,jpritbl-1
      zp=0.5_JPRB*aritbl(jjp-1)+2._jprb
      aritbl(jjp)=-0.5_JPRB/zp
      zu(jjp)=(6._jprb*((ritbl(jjp+1)-ritbl(jjp))/dritbl &
        & -(ritbl(jjp)-ritbl(jjp-1))&
        & /dritbl)/(2._jprb*dritbl)-0.5_JPRB*zu(jjp-1))/zp
    ENDDO

    zqn=0.5_JPRB
    zun=(3._jprb/dritbl)*(zdrvn-(ritbl(jpritbl)-ritbl(jpritbl-1))/dritbl)
    aritbl(jpritbl)=(zun-zqn*zu(jpritbl-1))/(zqn*aritbl(jpritbl-1)+1._jprb)

    !     4.3 BACK SUBSTITUTION

    DO jjp=jpritbl-1,1,-1
      aritbl(jjp)=aritbl(jjp)*aritbl(jjp+1)+zu(jjp)
    ENDDO

    !RETURN
  END SUBROUTINE suvdfs


!------------------------------------------------------------------------------


  SUBROUTINE DR_HOOK(CH,K1,P1)
    !USE PARKIND1  ,ONLY : JPIM     ,JPRB
    IMPLICIT NONE
    CHARACTER::CH
    INTEGER(KIND=JPIM)::K1
    REAL(KIND=JPRB)::P1
  END SUBROUTINE DR_HOOK 

  
!------------------------------------------------------------------------------
! Dummy routines vdiv, vexp, vrec (only use when V_MASS>0)


  SUBROUTINE vdiv(p1,p2,p3,k1)
    !USE PARKIND1  ,ONLY : JPIM     ,JPRB
    IMPLICIT NONE
    INTEGER(KIND=jpim)::k1
    REAL(KIND=jprb)::p1(k1),p2(k1),p3(k1)
  END SUBROUTINE vdiv

  
!------------------------------------------------------------------------------

  
  SUBROUTINE vexp(p1,p2,k1)
    
    !USE PARKIND1  ,ONLY : JPIM     ,JPRB
    IMPLICIT NONE
    INTEGER(KIND=jpim)::k1
    REAL(KIND=jprb)::p1(k1),p2(k1)
  END SUBROUTINE vexp

  
!------------------------------------------------------------------------------

  
  SUBROUTINE vrec(p1,p2,k1)
    
    !USE PARKIND1  ,ONLY : JPIM     ,JPRB
    IMPLICIT NONE
    INTEGER(KIND=jpim)::k1
    REAL(KIND=jprb)::p1(k1),p2(k1)
  END SUBROUTINE vrec


!------------------------------------------------------------------------------

  
  SUBROUTINE vlog(p1,p2,k1)
    
    !USE PARKIND1  ,ONLY : JPIM     ,JPRB
    IMPLICIT NONE
    INTEGER(KIND=jpim)::k1
    REAL(KIND=jprb)::p1(k1),p2(k1)
  END SUBROUTINE vlog


END MODULE mo_cuparameters

