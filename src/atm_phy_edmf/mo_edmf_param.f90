!>
!! Modules and associated setup routines for EDMF DUALM:
!! yoephy  & su0phy
!! yomsekf & susekf
!! yomjfh  & sujfh (not needed)
!! yomct0  & suct0
!! yomlun  & sulun
!!
!! @author Martin Koehler, DWD
!!
!! @par Revision History
!! Imported from IFS by Martin Koehler  (starting 2011-9-30)
!!   (IFS cycle CY36R1_DUALM_M8b)
!!
!!-----------------------------------------------------------------------------
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!-----------------------------------------------------------------------------

MODULE mo_edmf_param
 
  USE mo_kind,       ONLY: jpim=>i4, jprb=>wp

  IMPLICIT NONE

  PUBLIC

! NTILES for EDMF

  INTEGER(KIND=JPIM) :: NTILES_EDMF=8
  
! Type of EDMF setup: surface layer and TERRA calling
  INTEGER(KIND=JPIM) :: EDMF_CONF=2  ! 1: EDMF surface layer and TERRA called within EDMF
                                     ! 2: TURBTRAN surface layer and TERRA called from NWP interface


!------------------------------------------------------------------------------
! yoephy.F90 (part)
!     -----------------------------------------------------------------
!*    ** *YOEPHY* - SWITCHES RELATED TO DIABATIC PROCESSES
!     -----------------------------------------------------------------
!  NAME     TYPE     PURPOSE
!  ----  :  ----   : ---------------------------------------------------
! LVDFTRAC: LOGICAL: TURN TRACER TRANSPORT BY VERTICAL DIFFUSION ON 
! LEOCWA : LOGICAL : WARM OCEAN LAYER PARAMETRIZATION
! LEOCCO : LOGICAL : COOL OCEAN SKIN PARAMETRIZATION
! LEFLAKE : LOGICAL  : LAKE MODEL FLAKE


  LOGICAL :: LVDFTRAC
  LOGICAL :: LEOCWA
  LOGICAL :: LEOCCO
  LOGICAL :: LEOCSA    ! TRUE if SALINTY EFFECT ON SATURATION AT OCEAN SURFACE active
  LOGICAL :: LEVGEN    ! TRUE IF VAN GENUCHTEN HYDRO IS ACTIVATED
  LOGICAL :: LEFLAKE
  REAL (KIND = JPRB):: RH_ICE_MIN_FLK   ! Minimum ice thickness [m]


!------------------------------------------------------------------------------
! yomsekf.F90 (part)
!     ------------------------------------------------------------------
!*    Variables for the simplfied EKF soil moisture analysis.
!*    ------------------------------------------------------------------
!     Control variables for the EKF
!     -----------------------------
!     N_SEKF_PT  :    Number of perturbation run
!
!     Logicals
!     --------
!     LUSEKF_REF     : True if the SEKF soil moisture analysis is used
!                      (nconf = 302)
!     LUSE_JATM      : True is J computed from full 3D perturbed run

  INTEGER(KIND=JPIM) :: N_SEKF_PT
  LOGICAL :: LUSEKF_REF
  LOGICAL :: LUSE_JATM  ! switch for offline jacobians


!------------------------------------------------------------------------------
! yomct0.F90 (part)
!     ------------------------------------------------------------------
!*    Control variables for the job - constant within job
!========== ECMWF Single Column Model =========================================
! LSCMEC  : .T. = ECMWF Single Column Model
! LSFCFLX : .T. = forcing with surface fluxes (latent and sensible).
! REXTSHF : externally supplied sensible heat flux [W/m^2]
! REXTLHF : externally supplied latent   heat flux [W/m^2]
! LROUGH  : .T. = surface roughness length is externally specified 
! REXTZ0M : externally supplied roughness length for momentum [m]
! REXTZ0H : externally supplied roughness length for heat [m]

! * ECMWF Single Column Model:
LOGICAL :: LSCMEC
LOGICAL :: LSFCFLX
REAL(KIND=JPRB) :: REXTSHF
REAL(KIND=JPRB) :: REXTLHF
LOGICAL :: LROUGH
REAL(KIND=JPRB) :: REXTZ0M
REAL(KIND=JPRB) :: REXTZ0H


!------------------------------------------------------------------------------
! yomjfk.F90 (part)
!------------------------------------------------------------------------------
! Use of MASS library
!  N_VMASS: < or = 0 if not using MASS library vector routines
!           > 0      if using MASS library vector routines

INTEGER(KIND=JPIM) :: N_VMASS=0


!------------------------------------------------------------------------------
! yomlun.F90 (part)
!------------------------------------------------------------------------------
INTEGER(KIND=JPIM) :: NULERR=0


!------------------------------------------------------------------------------
! yos_exc.F90, suvexc_mod.F90, surf_inq.F90
!------------------------------------------------------------------------------
REAL(KIND=JPRB) :: REPUST=0.0001_JPRB ! MINIMUM FRICTION VELOCITY (SECURITY PARAMETER)


!------------------------------------------------------------------------------
! yos_veg.F90 (part)
!------------------------------------------------------------------------------
REAL(KIND=JPRB),ALLOCATABLE :: RVLAI(:)      ! LEAF AREA INDEX
REAL(KIND=JPRB),ALLOCATABLE :: RVROOTSA(:,:) ! PERCENTAGE OF ROOTS IN EACH SOIL LAYER
REAL(KIND=JPRB),ALLOCATABLE :: RVLAMSK(:)    ! Unstable SKIN LAYER CONDUCT. FOR EACH TILE
REAL(KIND=JPRB),ALLOCATABLE :: RVLAMSKS(:)   ! Stable SKIN LAYER CONDUCT. FOR EACH TILE
REAL(KIND=JPRB),ALLOCATABLE :: RVTRSR(:)     ! TRANSMISSION OF NET SOLAR RAD. 
                                             ! THROUGH VEG.
REAL(KIND=JPRB),ALLOCATABLE :: RVZ0M(:)      ! ROUGHNESS LENGTH FOR MOMENTUM
REAL(KIND=JPRB),ALLOCATABLE :: RVZ0H(:)      ! ROUGHNESS LENGTH FOR HEAT 
REAL(KIND=JPRB),ALLOCATABLE :: RVRSMIN(:)    ! MIN STOMATAL RESISTANCE FOR EACH VEG. TYPE (S/M)
REAL(KIND=JPRB),ALLOCATABLE :: RVHSTR(:)     ! HUMIDITY STRESS FUNCTION PARAMETER (M/S kgkg-1)
REAL(KIND=JPRB) :: RCEPSW                    ! MINIMUM RELATIVE HUMIDITY
REAL(KIND=JPRB) :: RLHAERO, RLHAEROS



!------------------------------------------------------------------------------
! yos_soil.F90 (part)
!------------------------------------------------------------------------------
REAL(KIND=JPRB),ALLOCATABLE :: RWCAPM(:)    !     RWCAP IN VAN GENUCHTEN
REAL(KIND=JPRB),ALLOCATABLE :: RWPWPM(:)    !     RWPWP IN VAN GENUCHTEN
REAL(KIND=JPRB),ALLOCATABLE :: RQWEVAPM(:)  !     RQWEVAP IN VAN GENUCHTEN
REAL(KIND=JPRB),ALLOCATABLE :: RWRESTM(:)   !     RWRST IN VAN GENUCHTEN
REAL(KIND=JPRB) :: RTF1            ! UPPER TEMPERATURE FOR SOIL WATER FREEZING
REAL(KIND=JPRB) :: RTF2            ! LOWER TEMPERATURE FOR SOIL WATER FREEZING
REAL(KIND=JPRB) :: RTF3            ! COEFFICIENT FOR SOIL WATER FREEZING FUNCTION
REAL(KIND=JPRB) :: RTF4            ! COEFFICIENT FOR SOIL WATER FREEZING FUNCTION


!------------------------------------------------------------------------------

  PUBLIC :: FOEALFA  ,FOEEWM   ,FOEDEM   ,FOELDCPM  , &
          & FOEALFCU ,FOEEWMCU ,FOEDEMCU ,FOELDCPMCU, &
          & FOEEW    ,FOEDESU  ,                      &
          & PSIHU    ,PSIMU    ,PSIHS    ,PSIMS     , &
          & RVZ0M    ,RVZ0H    ,RVLAMSK  ,RVLAMSKS  , &
          & RVTRSR   ,                                &
          & RWCAPM   ,RWPWPM   ,RQWEVAPM ,RWRESTM   , &
          & RTF1     ,RTF2     ,RTF3     ,RTF4      , &
          & RLHAERO  ,RLHAEROS
  PUBLIC :: suct0, su0phy, susekf, susveg, abort_surf

CONTAINS

!------------------------------------------------------------------------------
! fcttrm.h
!       This COMDECK includes the Thermodynamical functions for the cy39
!       ECMWF Physics package.
!       Consistent with YOMCST Basic physics constants, assuming the
!       partial pressure of water vapour is given by a first order
!       Taylor expansion of Qs(T) w.r.t. to Temperature, using constants
!       in YOETHF
!       Two sets of functions are available. In the first set only the
!       cases water or ice are distinguished by temperature.  This set 
!       consists of the functions FOEDELTA,FOEEW,FOEDE and FOELH.
!       The second set considers, besides the two cases water and ice 
!       also a mix of both for the temperature range RTICE < T < RTWAT.
!       This set contains FOEALFA,FOEEWM,FOEDEM,FOELDCPM and FOELHM.
!       FKOOP modifies the ice saturation mixing ratio for homogeneous 
!       nucleation

!       Depending on the consideration of mixed phases either the first 
!       set (e.g. surface, post-processing) or the second set 
!       (e.g. clouds, condensation, convection) should be used.
!------------------------------------------------------------------------------

ELEMENTAL FUNCTION foealfa(ptare)
  USE mo_cuparameters ,ONLY   :  rtwat, rtice, rtwat_rtice_r
  REAL(KIND=jprb)             :: foealfa
  REAL(KIND=jprb), INTENT(in) :: ptare
  FOEALFA = MIN(1.0_JPRB,((MAX(RTICE,MIN(RTWAT,PTARE))-RTICE)&
          & *RTWAT_RTICE_R)**2)
END FUNCTION foealfa

ELEMENTAL FUNCTION foeewm(ptare)
  USE mo_cuparameters ,ONLY   :  r2es, r3les, rtt, r4les, r3ies, r4ies
  REAL(KIND=jprb)             :: foeewm
  REAL(KIND=jprb), INTENT(in) :: ptare
  FOEEWM = R2ES *&
         &(FOEALFA(PTARE)*EXP(R3LES*(PTARE-RTT)/(PTARE-R4LES))+&
         &(1.0_JPRB-FOEALFA(PTARE))*EXP(R3IES*(PTARE-RTT)/(PTARE-R4IES)))
END FUNCTION foeewm

ELEMENTAL FUNCTION FOEDEM(ptare)
  USE mo_cuparameters ,ONLY   :  r5alvcp, r4les, r5alscp, r4ies
  REAL(KIND=jprb)             :: FOEDEM
  REAL(KIND=jprb), INTENT(in) :: ptare
  FOEDEM = FOEALFA(PTARE)*R5ALVCP*(1.0_JPRB/(PTARE-R4LES)**2)+&
             &(1.0_JPRB-FOEALFA(PTARE))*R5ALSCP*(1.0_JPRB/(PTARE-R4IES)**2)
END FUNCTION FOEDEM

ELEMENTAL FUNCTION FOELDCPM(ptare)
  USE mo_cuparameters ,ONLY   :  ralvdcp, ralsdcp
  REAL(KIND=jprb)             :: FOELDCPM
  REAL(KIND=jprb), INTENT(in) :: ptare
  FOELDCPM = FOEALFA(PTARE)*RALVDCP+&
            &(1.0_JPRB-FOEALFA(PTARE))*RALSDCP
END FUNCTION FOELDCPM

ELEMENTAL FUNCTION foealfcu(ptare)
  USE mo_cuparameters ,ONLY   :  rtwat, rticecu, rtwat_rticecu_r
  REAL(KIND=jprb)             :: foealfcu
  REAL(KIND=jprb), INTENT(in) :: ptare
  FOEALFCU = MIN(1.0_JPRB,((MAX(RTICECU,MIN(RTWAT,PTARE))&
           &-RTICECU)*RTWAT_RTICECU_R)**2) 
END FUNCTION foealfcu

ELEMENTAL FUNCTION FOEEWMCU(ptare)
  USE mo_cuparameters ,ONLY   :  r2es, r3les, rtt, r4les, r3ies, r4ies
  REAL(KIND=jprb)             :: FOEEWMCU
  REAL(KIND=jprb), INTENT(in) :: ptare
  FOEEWMCU = R2ES *&
     &(FOEALFCU(PTARE)*EXP(R3LES*(PTARE-RTT)/(PTARE-R4LES))+&
  &(1.0_JPRB-FOEALFCU(PTARE))*EXP(R3IES*(PTARE-RTT)/(PTARE-R4IES)))
END FUNCTION FOEEWMCU

ELEMENTAL FUNCTION FOEDEMCU(ptare)
  USE mo_cuparameters ,ONLY   :  r5alvcp, r4les, r5alscp, r4ies 
  REAL(KIND=jprb)             :: FOEDEMCU
  REAL(KIND=jprb), INTENT(in) :: ptare
  FOEDEMCU = FOEALFCU(PTARE)*R5ALVCP*(1.0_JPRB/(PTARE-R4LES)**2)+&
             &(1.0_JPRB-FOEALFCU(PTARE))*R5ALSCP*(1.0_JPRB/(PTARE-R4IES)**2)
END FUNCTION FOEDEMCU

ELEMENTAL FUNCTION FOELDCPMCU(ptare)
  USE mo_cuparameters ,ONLY   :  ralvdcp, ralsdcp
  REAL(KIND=jprb)             :: FOELDCPMCU
  REAL(KIND=jprb), INTENT(in) :: ptare
  FOELDCPMCU = FOEALFCU(PTARE)*RALVDCP+&
            &(1.0_JPRB-FOEALFCU(PTARE))*RALSDCP
END FUNCTION FOELDCPMCU


!------------------------------------------------------------------------------


!     ------------------------------------------------------------------
! fcsttre.h
!     This COMDECK includes the Thermodynamical functions for the cy39
!       ECMWF Physics package.
!       Consistent with YOMCST Basic physics constants, assuming the
!       partial pressure of water vapour is given by a first order
!       Taylor expansion of Qs(T) w.r.t. to Temperature, using constants
!       in YOETHF
!       Two sets of functions are available. In the first set only the
!       cases water or ice are distinguished by temperature.  This set 
!       consists of the functions FOEDELTA,FOEEW,FOEDE and FOELH.
!       The second set considers, besides the two cases water and ice 
!       also a mix of both for the temperature range RTICE < T < RTWAT.
!       This set contains FOEALFA,FOEEWM,FOEDEM,FOELDCPM and FOELHM.

!       Depending on the consideration of mixed phases either the first 
!       set (e.g. surface, post-processing) or the second set 
!       (e.g. clouds, condensation, convection) should be used.
!     ------------------------------------------------------------------
!     *****************************************************************
!                NO CONSIDERATION OF MIXED PHASES
!     *****************************************************************

!                  FOEDELTA = 1    water
!                  FOEDELTA = 0    ice

!     THERMODYNAMICAL FUNCTIONS .

!     Pressure of water vapour at saturation
!        INPUT : PTARE = TEMPERATURE

ELEMENTAL FUNCTION foedelta(ptare)
  USE mo_cuparameters ,ONLY   :  rtt
  REAL(KIND=JPRB)             :: FOEDELTA
  REAL(KIND=jprb), INTENT(in) :: ptare
  FOEDELTA = MAX (0.0_JPRB,SIGN(1.0_JPRB,PTARE-RTT))
END FUNCTION foedelta

ELEMENTAL FUNCTION foeew(ptare)
  USE mo_cuparameters ,ONLY   :  r2es, r3les, rtt, r4les, r3ies, r4ies
  REAL(KIND=jprb)             :: foeew
  REAL(KIND=jprb), INTENT(in) :: ptare
  FOEEW = R2ES*EXP (&
    & (R3LES*FOEDELTA(PTARE)+R3IES*(1.0_JPRB-FOEDELTA(PTARE)))*(PTARE-RTT)&
    & / (PTARE-(R4LES*FOEDELTA(PTARE)+R4IES*(1.0_JPRB-FOEDELTA(PTARE)))))
END FUNCTION foeew

ELEMENTAL FUNCTION foedesu(ptare)
  USE mo_cuparameters ,ONLY   :  r4les, r4ies, r5les, r5ies
  REAL(KIND=jprb)             :: FOEDESU
  REAL(KIND=jprb), INTENT(in) :: ptare
  FOEDESU = &
    &(FOEDELTA(PTARE)*R5LES+(1.0_JPRB-FOEDELTA(PTARE))*R5IES)&
    &/ (PTARE-(R4LES*FOEDELTA(PTARE)+R4IES*(1.0_JPRB-FOEDELTA(PTARE))))**2
END FUNCTION foedesu


!------------------------------------------------------------------------------


!     ------------------------------------------------------------------
! fcsvdfs.h
!
!     *FCVDFS** CONTAINS STATEMENT FUNCTIONS DESCRIBING STAB. FUNCT.
!
!     A.C.M. BELJAARS    E.C.M.W.F.      26/03/90.
!
!     ------------------------------------------------------------------
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
!
!
!     ------------------------------------------------------------------
!
!        *PHI AND *PSI FUNCTIONS FOR UNSTABLE SITUATIONS ACCORDING
!        TO *DYER AND *HICKS
!        3/2/03 Modification, J.Hague  X**1.5 -> X*SQRT(X)
!
!     ------------------------------------------------------------------

ELEMENTAL FUNCTION psihu(peta)
  USE mo_cuparameters ,ONLY   :  rcdhalf
  REAL(KIND=JPRB)             :: PSIHU
  REAL(KIND=JPRB), INTENT(in) :: PETA
  PSIHU= 2.0_JPRB*LOG((1.0_JPRB+     SQRT(1.0_JPRB-RCDHALF*PETA))*0.5_JPRB )
END FUNCTION psihu

ELEMENTAL FUNCTION psimu(peta)
  USE mo_cuparameters ,ONLY   :  rcdhalf, rcdhpi2
  REAL(KIND=JPRB)             :: PSIMU
  REAL(KIND=JPRB), INTENT(in) :: PETA
  PSIMU= LOG((1.0_JPRB+SQRT( SQRT(1.0_JPRB-RCDHALF*PETA)))**2 &
        &*(1.0_JPRB+         SQRT(1.0_JPRB-RCDHALF*PETA) ) *0.125_JPRB )&
        &-2.0_JPRB*ATAN(SQRT(SQRT(1.0_JPRB-RCDHALF*PETA)))&
        &+ RCDHPI2
END FUNCTION psimu

!        *PHI AND *PSI FUNCTIONS FOR UNSTABLE SITUATIONS ACCORDING
!        TO HOGSTROM FOR MOMENTUM AND DERIVED FROM THE ELLISON AND
!        TURNER RELATION FOR THE RATIO OF PHIM AMD PHIH.

!     PSI FUNCTIONS NOT COMPATIBLE WITH PHI FUNCTIONS IN THE STABLE CASE

ELEMENTAL FUNCTION psihs(peta)
  USE mo_cuparameters ,ONLY   :  rchbb, rchbcd, rchbd, rchb23a, rchbbcd
  REAL(KIND=JPRB)             :: PSIHS
  REAL(KIND=JPRB), INTENT(in) :: PETA
  PSIHS= -RCHBB*(PETA-RCHBCD)*EXP  (-RCHBD*PETA)&
        &-(1.0_JPRB+RCHB23A*PETA)*SQRT(1.0_JPRB+RCHB23A*PETA)-RCHBBCD+1.0_JPRB
END FUNCTION psihs

ELEMENTAL FUNCTION psims(peta)
  USE mo_cuparameters ,ONLY   :  rchbb, rchbcd, rchbd, rchba, rchbbcd
  REAL(KIND=JPRB)             :: PSIMS
  REAL(KIND=JPRB), INTENT(in) :: PETA
  PSIMS= -RCHBB*(PETA-RCHBCD)*EXP  (-RCHBD*PETA)&
        &-RCHBA*PETA - RCHBBCD
END FUNCTION psims


!------------------------------------------------------------------------------


SUBROUTINE suct0
!**** *SUCT0*   - Routine to initialize level 0 control common

!     Purpose.
!     --------
!           Initialize level 0 control commons

! ECMWF Single Column Model off by default
LSCMEC=.FALSE.

END SUBROUTINE suct0


!------------------------------------------------------------------------------


SUBROUTINE su0phy
!**** *SU0PHY*   - Initialize common YOxPHY controlling physics

!     Purpose.
!     --------
!           Initialize YOxPHY, the common that includes the
!           basic switches for the physics of the model.

LVDFTRAC=.TRUE.
LEOCWA=.TRUE.
LEOCCO=.TRUE.
LEOCSA=.TRUE.
LEVGEN=.TRUE.
LEFLAKE=.FALSE.
RH_ICE_MIN_FLK  = 1.0E-9_JPRB          ! Minimum ice thickness [m]

END SUBROUTINE su0phy


!------------------------------------------------------------------------------


SUBROUTINE susekf
!**** *SUSEKF*  - Routine to setup SEKF

!     Purpose.
!     --------
!           Setup the discriptors for initializing the SEKF

N_SEKF_PT = 0_JPIM
LUSEKF_REF=.FALSE.
LUSE_JATM = .False.  ! Switch to offline jacobians

END SUBROUTINE susekf


!------------------------------------------------------------------------------


SUBROUTINE SUSVEG
!**   *SUSVEG* IS THE SET-UP ROUTINE FOR COMMON BLOCK *YOS_VEG*

!     PURPOSE
!     -------
!          THIS ROUTINE INITIALIZES THE CONSTANTS IN COMMON BLOCK
!     *YOS_VEG*

!USE mo_lnd_nwp_config,       ONLY: NTILES=>ntiles_total
USE mo_srfrootfr     ,       ONLY: srfrootfr
USE mo_cuparameters  ,       ONLY: RETV

INTEGER(KIND=JPIM) :: NVTYPES, IVTYPES, NCSS, JS
REAL(KIND=JPRB)    :: ZLARGE , ZSNOW
REAL(KIND=JPRB), ALLOCATABLE :: RDAW(:)

! Number of vegetation types
NVTYPES=20
IVTYPES = MAX(NVTYPES,20)

NCSS = 4

!     CONSTANTS DEFINING SOIL DISCRETIZATION

IF(.NOT.ALLOCATED(RDAW)) ALLOCATE(RDAW(NCSS))
IF (NCSS >= 1) RDAW(1)=0.07_JPRB
IF (NCSS >= 2) RDAW(2)=0.21_JPRB
IF (NCSS >= 3) RDAW(3)=0.72_JPRB
IF (NCSS >= 4) RDAW(4)=1.89_JPRB

! Leaf area index

IF(.NOT.ALLOCATED(RVLAI)) ALLOCATE(RVLAI(0:IVTYPES))
RVLAI(1)=3._JPRB     ! Crops, Mixed Farming
RVLAI(2)=2._JPRB     ! Short Grass
RVLAI(3)=5._JPRB     ! Evergreen Needleleaf Trees
RVLAI(4)=5._JPRB     ! Deciduous Needleleaf Trees
RVLAI(5)=5._JPRB     ! Deciduous Broadleaf Trees
RVLAI(6)=6._JPRB     ! Evergreen Broadleaf Trees
RVLAI(7)=2._JPRB     ! Tall Grass
RVLAI(8)=0.5_JPRB    ! Desert
RVLAI(9)=1.0_JPRB    ! Tundra
RVLAI(10)=3._JPRB    ! Irrigated Crops
RVLAI(11)=0.5_JPRB   ! Semidesert
RVLAI(12)=0.0_JPRB   ! Ice Caps and Glaciers
RVLAI(13)=4._JPRB    ! Bogs and Marshes
RVLAI(14)=0.0_JPRB   ! Inland Water
RVLAI(15)=0.0_JPRB   ! Ocean
RVLAI(16)=3._JPRB    ! Evergreen Shrubs
RVLAI(17)=1.5_JPRB   ! Deciduous Shrubs
RVLAI(18)=5._JPRB    ! Mixed Forest/woodland
RVLAI(19)=2.5_JPRB   ! Interrupted Forest
RVLAI(20)=4._JPRB    ! Water and Land Mixtures
RVLAI(0)=RVLAI(8)

! Skin layer conductivity

ZSNOW=7._JPRB
ZLARGE=1.E10_JPRB
RLHAERO=15._JPRB
RLHAEROS=10._JPRB

! Unstable skin layer conductivity

IF(.NOT.ALLOCATED(RVLAMSK)) ALLOCATE(RVLAMSK(NTILES_EDMF))
IF (NTILES_EDMF >= 1) RVLAMSK(1)=ZLARGE           !  Water
IF (NTILES_EDMF >= 2) RVLAMSK(2)=58._JPRB         !  Sea ice
IF (NTILES_EDMF >= 3) RVLAMSK(3)=10._JPRB         !  Wet skin
IF (NTILES_EDMF >= 4) RVLAMSK(4)=10._JPRB         !  Low veg. (Aerodyn.+Rad.+Soil)
IF (NTILES_EDMF >= 5) RVLAMSK(5)=ZSNOW            !  Snow
IF (NTILES_EDMF >= 5) RVLAMSK(6)=RLHAERO+5._JPRB  !  High veg. (Aerodyn.+Rad.)
IF (NTILES_EDMF >= 7) RVLAMSK(7)=RLHAERO+5._JPRB  !  Snow under veg (Aerodyn.+Rad.)
IF (NTILES_EDMF >= 8) RVLAMSK(8)=15._JPRB         !  Bare soil
IF (NTILES_EDMF >= 9) RVLAMSK(9)=ZLARGE           !  LAKES = WATER  

! Stable skin layer conductivity

IF(.NOT.ALLOCATED(RVLAMSKS)) ALLOCATE(RVLAMSKS(NTILES_EDMF))
IF (NTILES_EDMF >= 1) RVLAMSKS(1)=ZLARGE            !  Water
IF (NTILES_EDMF >= 2) RVLAMSKS(2)=58._JPRB          !  Sea ice
IF (NTILES_EDMF >= 3) RVLAMSKS(3)=10._JPRB          !  Wet skin
IF (NTILES_EDMF >= 4) RVLAMSKS(4)=10._JPRB          !  Low veg. (Aerodyn.+Rad.+Soil)
IF (NTILES_EDMF >= 5) RVLAMSKS(5)=ZSNOW             !  Snow
IF (NTILES_EDMF >= 5) RVLAMSKS(6)=RLHAEROS+5._JPRB  !  High veg. (Aerodyn.+Rad.)
IF (NTILES_EDMF >= 7) RVLAMSKS(7)=RLHAEROS+5._JPRB  !  Snow under veg (Aerodyn.+Rad.+Snow
IF (NTILES_EDMF >= 8) RVLAMSKS(8)=15._JPRB          !  Bare soil
IF (NTILES_EDMF >= 9) RVLAMSKS(9)=ZLARGE            !  LAKES = WATER  

! Transmission of net solar rad. through vegetation

IF(.NOT.ALLOCATED(RVTRSR)) ALLOCATE(RVTRSR(NTILES_EDMF))
IF (NTILES_EDMF >= 1) RVTRSR(1)=0.00_JPRB     ! Ocean (SSR transmission does not apply)
IF (NTILES_EDMF >= 2) RVTRSR(2)=0.00_JPRB     ! Ice   (SSR transmission does not apply)
IF (NTILES_EDMF >= 3) RVTRSR(3)=0.05_JPRB     ! Wet skin (this is a compromise)
IF (NTILES_EDMF >= 4) RVTRSR(4)=0.05_JPRB     ! Low vegetation
IF (NTILES_EDMF >= 5) RVTRSR(5)=0.00_JPRB     ! Snow on vegetation (SSR transmission does not apply)
IF (NTILES_EDMF >= 5) RVTRSR(6)=0.03_JPRB     ! High veg.
IF (NTILES_EDMF >= 7) RVTRSR(7)=0.03_JPRB     ! Snow under veg.
IF (NTILES_EDMF >= 8) RVTRSR(8)=0.00_JPRB     ! Bare soil (SSR transmission does not apply)
IF (NTILES_EDMF >= 9) RVTRSR(9)=0.00_JPRB     ! LAKES (SSR transmission does not apply) 

! Root fraction

IF(.NOT.ALLOCATED(RVROOTSA)) ALLOCATE(RVROOTSA(NCSS,0:NVTYPES))
IF (NCSS >= 1) CALL SRFROOTFR(NCSS,NVTYPES,RDAW,RVROOTSA(:,1:NVTYPES))
DO JS=1,NCSS
  RVROOTSA(JS,0)=RVROOTSA(JS,8)
ENDDO

! Set other constants

RCEPSW  =1.E-3_JPRB

! Minimum stomatal resitance for each vegetation type (s/m)

!                 The additinal correction coefficient (1.30 etc.) accounts
!                 for the decreased radiation stress function at full 
!                 light saturation.

IF(.NOT.ALLOCATED(RVRSMIN)) ALLOCATE (RVRSMIN(0:IVTYPES))
RVRSMIN(1)=180._JPRB    ! Crops, Mixed Farming
RVRSMIN(2)=110._JPRB    ! Short Grass
RVRSMIN(3)=500._JPRB    ! Evergreen Needleleaf Trees
RVRSMIN(4)=500._JPRB    ! Deciduous Needleleaf Trees
RVRSMIN(5)=175._JPRB    ! Deciduous Broadleaf Trees
RVRSMIN(6)=240._JPRB    ! Evergreen Broadleaf Trees
RVRSMIN(7)=100._JPRB    ! Tall Grass
RVRSMIN(8)=250._JPRB    ! Desert
RVRSMIN(9)=80._JPRB     ! Tundra
RVRSMIN(10)=180._JPRB   ! Irrigated Crops
RVRSMIN(11)=150._JPRB   ! Semidesert
RVRSMIN(12)=0.0_JPRB    ! Ice Caps and Glaciers
RVRSMIN(13)=240._JPRB   ! Bogs and Marshes
RVRSMIN(14)=0.0_JPRB    ! Inland Water
RVRSMIN(15)=0.0_JPRB    ! Ocean
RVRSMIN(16)=225._JPRB   ! Evergreen Shrubs
RVRSMIN(17)=225._JPRB   ! Deciduous Shrubs
RVRSMIN(18)=250._JPRB   ! Mixed Forest/woodland
RVRSMIN(19)=175._JPRB   ! Interrupted Forest
RVRSMIN(20)=150._JPRB   ! Water and Land Mixtures
RVRSMIN(0)=RVRSMIN(8)

! Parameter in humidity stress function (m/s mbar, converted to M/S kgkg-1)
IF(.NOT.ALLOCATED(RVHSTR)) ALLOCATE (RVHSTR(0:IVTYPES))
RVHSTR(1)=0.0_JPRB     ! Crops, Mixed Farming
RVHSTR(2)=0.0_JPRB     ! Short Grass
RVHSTR(3)=0.03_JPRB    ! Evergreen Needleleaf Trees
RVHSTR(4)=0.03_JPRB    ! Deciduous Needleleaf Trees
RVHSTR(5)=0.03_JPRB    ! Deciduous Broadleaf Trees
RVHSTR(6)=0.03_JPRB    ! Evergreen Broadleaf Trees
RVHSTR(7)=0.0_JPRB     ! Tall Grass
RVHSTR(8)=0.0_JPRB     ! Desert
RVHSTR(9)=0.0_JPRB     ! Tundra
RVHSTR(10)=0.0_JPRB    ! Irrigated Crops
RVHSTR(11)=0.0_JPRB    ! Semidesert
RVHSTR(12)=0.0_JPRB    ! Ice Caps and Glaciers
RVHSTR(13)=0.0_JPRB    ! Bogs and Marshes
RVHSTR(14)=0.0_JPRB    ! Inland Water
RVHSTR(15)=0.0_JPRB    ! Ocean
RVHSTR(16)=0.0_JPRB    ! Evergreen Shrubs
RVHSTR(17)=0.0_JPRB    ! Deciduous Shrubs
RVHSTR(18)=0.03_JPRB   ! Mixed Forest/woodland
RVHSTR(19)=0.03_JPRB   ! Interrupted Forest
RVHSTR(20)=0.0_JPRB    ! Water and Land Mixtures
RVHSTR(0)=RVHSTR(8)
RVHSTR(:)=1013.25_JPRB*(RETV+1)*RVHSTR(:)

! Roughness length for momentum (Mahfouf et al. 1995)

IF(.NOT.ALLOCATED(RVZ0M)) ALLOCATE (RVZ0M(0:IVTYPES))
RVZ0M(1)=0.15_JPRB     ! Crops, Mixed Farming
RVZ0M(2)=0.02_JPRB     ! Short Grass
RVZ0M(3)=2.00_JPRB     ! Evergreen Needleleaf Trees
RVZ0M(4)=2.00_JPRB     ! Deciduous Needleleaf Trees
RVZ0M(5)=2.00_JPRB     ! Deciduous Broadleaf Trees
RVZ0M(6)=2.00_JPRB     ! Evergreen Broadleaf Trees
RVZ0M(7)=0.10_JPRB     ! Tall Grass
RVZ0M(8)=0.013_JPRB    ! Desert                    # Masson et al.
RVZ0M(9)=0.05_JPRB     ! Tundra
RVZ0M(10)=0.15_JPRB    ! Irrigated Crops           # Crops type 1
RVZ0M(11)=0.05_JPRB    ! Semidesert 
RVZ0M(12)=0.0013_JPRB  ! Ice Caps and Glaciers     # Mason et al. 
RVZ0M(13)=0.05_JPRB    ! Bogs and Marshes
RVZ0M(14)=0.0001_JPRB  ! Inland Water              # Not used but needs value here
RVZ0M(15)=0.0001_JPRB  ! Ocean                     # Not used but needs value here
RVZ0M(16)=0.10_JPRB    ! Evergreen Shrubs
RVZ0M(17)=0.10_JPRB    ! Deciduous Shrubs
RVZ0M(18)=2.00_JPRB    ! Mixed Forest/woodland
RVZ0M(19)=0.50_JPRB    ! Interrupted Forest        # New value invented here
RVZ0M(20)=0.02_JPRB    ! Water and Land Mixtures   # Not used but needs value here
RVZ0M(0)=RVZ0M(8)      !                           # Bare soil value

! Roughness length for heat

!IF(.NOT.ALLOCATED(RVZ0H)) ALLOCATE (RVZ0H(0:IVTYPES))
!RVZ0H(1)=RVZ0M( 1)/10._JPRB     ! Crops, Mixed Farming
!RVZ0H(2)=RVZ0M( 2)/10._JPRB     ! Short Grass
!RVZ0H(3)=RVZ0M( 3)              ! Evergreen Needleleaf Trees
!RVZ0H(4)=RVZ0M( 4)              ! Deciduous Needleleaf Trees
!RVZ0H(5)=RVZ0M( 5)              ! Deciduous Broadleaf Trees
!RVZ0H(6)=RVZ0M( 6)              ! Evergreen Broadleaf Trees
!RVZ0H(7)=RVZ0M( 7)/10._JPRB     ! Tall Grass
!RVZ0H(8)=RVZ0M( 8)/10._JPRB     ! Desert 
!RVZ0H(9)=RVZ0M( 9)/10._JPRB     ! Tundra
!RVZ0H(10)=RVZ0M(10)/10._JPRB    ! Irrigated Crops 
!RVZ0H(11)=RVZ0M(11)/10._JPRB    ! Semidesert 
!RVZ0H(12)=RVZ0M(12)/10._JPRB    ! Ice Caps and Glaciers 
!RVZ0H(13)=RVZ0M(13)/10._JPRB    ! Bogs and Marshes
!RVZ0H(14)=RVZ0M(14)/10._JPRB    ! Inland Water  
!RVZ0H(15)=RVZ0M(15)/10._JPRB    ! Ocean 
!RVZ0H(16)=RVZ0M(16)/10._JPRB    ! Evergreen Shrubs
!RVZ0H(17)=RVZ0M(17)/10._JPRB    ! Deciduous Shrubs
!RVZ0H(18)=RVZ0M(18)             ! Mixed Forest/woodland
!RVZ0H(19)=RVZ0M(19)/10._JPRB    ! Interrupted Forest 
!RVZ0H(20)=RVZ0M(20)/10._JPRB    ! Water and Land Mixtures 
!RVZ0H(0)=RVZ0H(8)        

!xmk: change z0's to update to cy28r2 (small z0/100)
!     AND lower tree values (~2m) by factor 10
!     ???

IF(.NOT.ALLOCATED(RVZ0H)) ALLOCATE (RVZ0H(0:IVTYPES))
RVZ0H(1)=RVZ0M( 1)/100._JPRB     ! Crops, Mixed Farming
RVZ0H(2)=RVZ0M( 2)/100._JPRB     ! Short Grass
RVZ0H(3)=RVZ0M( 3)/10._JPRB      ! Evergreen Needleleaf Trees
RVZ0H(4)=RVZ0M( 4)/10._JPRB      ! Deciduous Needleleaf Trees
RVZ0H(5)=RVZ0M( 5)/10._JPRB      ! Deciduous Broadleaf Trees
RVZ0H(6)=RVZ0M( 6)/10._JPRB      ! Evergreen Broadleaf Trees
RVZ0H(7)=RVZ0M( 7)/100._JPRB     ! Tall Grass
RVZ0H(8)=RVZ0M( 8)/100._JPRB     ! Desert 
RVZ0H(9)=RVZ0M( 9)/100._JPRB     ! Tundra
RVZ0H(10)=RVZ0M(10)/100._JPRB    ! Irrigated Crops 
RVZ0H(11)=RVZ0M(11)/100._JPRB    ! Semidesert 
RVZ0H(12)=RVZ0M(12)/10._JPRB     ! Ice Caps and Glaciers 
RVZ0H(13)=RVZ0M(13)/100._JPRB    ! Bogs and Marshes
RVZ0H(14)=RVZ0M(14)/10._JPRB     ! Inland Water  
RVZ0H(15)=RVZ0M(15)/10._JPRB     ! Ocean 
RVZ0H(16)=RVZ0M(16)/100._JPRB    ! Evergreen Shrubs
RVZ0H(17)=RVZ0M(17)/100._JPRB    ! Deciduous Shrubs
RVZ0H(18)=RVZ0M(18)/10._JPRB     ! Mixed Forest/woodland
RVZ0H(19)=RVZ0M(19)/100._JPRB    ! Interrupted Forest 
RVZ0H(0)=RVZ0H(8)        
RVZ0H(20)=RVZ0M(20)/10._JPRB     ! Water and Land Mixtures 

END SUBROUTINE SUSVEG


!------------------------------------------------------------------------------


SUBROUTINE SUSSOIL
!**   *SUSSSOIL* IS THE SET-UP ROUTINE FOR COMMON BLOCK *YOESOIL*

!     PURPOSE
!     -------
!          THIS ROUTINE INITIALIZES THE CONSTANTS IN COMMON BLOCK
!     *YOESOIL*

USE mo_cuparameters ,ONLY : RTT, RPI

IMPLICIT NONE

INTEGER(KIND=JPIM) :: NSOTY, JS

REAL(KIND=JPRB), ALLOCATABLE :: ZMVGALPHA(:),ZNFAC(:),ZMFAC(:),&
 & ZWSATM(:),ZWCAPM(:),ZWRES(:),ZWPWPM(:)

REAL(KIND=JPRB) :: ZPSIPWP, ZPSICAP, ZFAC

  NSOTY=7
  IF(.NOT.ALLOCATED(ZMVGALPHA)) ALLOCATE (ZMVGALPHA(0:NSOTY))
  IF(.NOT.ALLOCATED(ZNFAC)) ALLOCATE (ZNFAC(0:NSOTY))
  IF(.NOT.ALLOCATED(ZMFAC)) ALLOCATE (ZMFAC(0:NSOTY))
  IF(.NOT.ALLOCATED(ZWSATM)) ALLOCATE (ZWSATM(0:NSOTY))
  IF(.NOT.ALLOCATED(ZWCAPM)) ALLOCATE (ZWCAPM(0:NSOTY))
  IF(.NOT.ALLOCATED(ZWRES)) ALLOCATE (ZWRES(0:NSOTY))
  IF(.NOT.ALLOCATED(ZWPWPM)) ALLOCATE (ZWPWPM(0:NSOTY))
  IF(.NOT.ALLOCATED(RWCAPM)) ALLOCATE (RWCAPM(0:NSOTY))
  IF(.NOT.ALLOCATED(RWPWPM)) ALLOCATE (RWPWPM(0:NSOTY))
  IF(.NOT.ALLOCATED(RQWEVAPM)) ALLOCATE (RQWEVAPM(0:NSOTY))
  IF(.NOT.ALLOCATED(RWRESTM)) ALLOCATE (RWRESTM(0:NSOTY))
  
  ZMVGALPHA(1:NSOTY)=(/3.83_JPRB,3.14_JPRB,0.83_JPRB,3.67_JPRB,2.65_JPRB,1.300_JPRB,3.14_JPRB/)
  ZNFAC(1:NSOTY)=(/1.3774_JPRB,1.1804_JPRB,1.2539_JPRB,1.1012_JPRB,1.1033_JPRB,1.2039_JPRB,&
                &1.1804_JPRB/)
  ZWSATM(1:NSOTY)=(/0.403_JPRB,0.439_JPRB,0.430_JPRB,0.520_JPRB,0.614_JPRB,0.766_JPRB,0.439_JPRB/)
  ZWRES(1:NSOTY)=(/0.025_JPRB,0.010_JPRB,0.010_JPRB,0.010_JPRB,0.010_JPRB,0.010_JPRB,0.010_JPRB/)

  ZMVGALPHA(0)=0.0_JPRB
  ZNFAC(0)=0.0_JPRB
  ZWSATM(0)=0.0_JPRB
  ZWRES(0)=0.0_JPRB
  RWCAPM(0)=0.0_JPRB
  RWPWPM(0)=0.0_JPRB
  RQWEVAPM(0)=0.0_JPRB

  ZFAC=1000._JPRB*100._JPRB/(1000._JPRB*9.8_JPRB)
  ZPSIPWP=-15._JPRB*ZFAC   ! Classic permanent wilting point
  IF (LEVGEN) THEN
    ZPSICAP=-0.10_JPRB*ZFAC  ! Value valid for sandy soil in the literature---> larger AWC 
                             ! favourably compared to FC obs for medium soils (Ukraine/Russia)
                             ! produce invariant equilibrium after rescaling from TESSEL
                             ! enlarge the range of action for SM analysis (active if PWP<W<CAP)
  ELSE
    ZPSICAP=-0.33_JPRB*ZFAC  ! Value mostly used in literature (Hillel, 1998)
  ENDIF

  DO JS=1,NSOTY
    ZMFAC(JS)=1._JPRB-(1._JPRB/ZNFAC(JS))
    ZWCAPM(JS)=ZWRES(JS)+(ZWSATM(JS)-ZWRES(JS)) &
 &           *(1._JPRB/(1._JPRB+((ABS(ZMVGALPHA(JS)*ZPSICAP)) &
 &           **ZNFAC(JS))))**ZMFAC(JS)
    RWCAPM(JS)=ZWCAPM(JS)
    ZWPWPM(JS)=ZWRES(JS)+(ZWSATM(JS)-ZWRES(JS)) &
 &           *(1._JPRB/(1._JPRB+((ABS(ZMVGALPHA(JS)*ZPSIPWP)) &
 &           **ZNFAC(JS))))**ZMFAC(JS)
    RWPWPM(JS)=ZWPWPM(JS)
    RQWEVAPM(JS)=1._JPRB/(RWCAPM(JS)-RWPWPM(JS))
  ENDDO

!     CONSTANTS FOR SOIL WATER FREEZING

RTF1=RTT+1.0_JPRB
RTF2=RTT-3._JPRB
RTF3=0.5_JPRB*(RTF1+RTF2)
RTF4=RPI/(RTF1-RTF2)

END SUBROUTINE SUSSOIL


!------------------------------------------------------------------------------


SUBROUTINE ABORT_SURF(text)
  USE mo_io_units, ONLY: nerr
  USE mo_mpi, ONLY: abort_mpi
  CHARACTER (len=*), INTENT(in) :: text
  WRITE (nerr,'(a)')  TRIM(text)
  CALL abort_mpi
END SUBROUTINE ABORT_SURF



END MODULE mo_edmf_param
