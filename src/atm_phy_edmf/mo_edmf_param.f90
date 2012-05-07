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
!!-----------------------------------------------------------------------------

MODULE mo_edmf_param
 
  USE mo_kind,       ONLY: jpim=>i4, jprb=>wp

  IMPLICIT NONE


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

  PUBLIC

  LOGICAL :: LVDFTRAC
  LOGICAL :: LEOCWA
  LOGICAL :: LEOCCO
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
REAL(KIND=JPRB),ALLOCATABLE :: RVZ0M(:)      ! ROUGHNESS LENGTH FOR MOMENTUM
REAL(KIND=JPRB),ALLOCATABLE :: RVZ0H(:)      ! ROUGHNESS LENGTH FOR HEAT 
REAL(KIND=JPRB),ALLOCATABLE :: RVLAMSK(:)    ! Unstable SKIN LAYER CONDUCT. FOR EACH TILE
REAL(KIND=JPRB),ALLOCATABLE :: RVLAMSKS(:)   ! Stable SKIN LAYER CONDUCT. FOR EACH TILE
REAL(KIND=JPRB),ALLOCATABLE :: RVTRSR(:)     ! TRANSMISSION OF NET SOLAR RAD. 
                                             ! THROUGH VEG.


!------------------------------------------------------------------------------

  PUBLIC :: FOEALFA  ,FOEEWM   ,FOEDEM   ,FOELDCPM  , &
          & FOEALFCU ,FOEEWMCU ,FOEDEMCU ,FOELDCPMCU, &
          & FOEEW    ,FOEDESU  ,                      &
          & PSIHU    ,PSIMU    ,PSIHS    ,PSIMS     , &
          & RVZ0M    ,RVZ0H    ,RVLAMSK  ,RVLAMSKS  , &
          & RVTRSR
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
  USE mo_cuparameters ,ONLY   :  r3les, rtt, r4les, r4ies, r5les, r5ies
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
LEOCWA=.FALSE.
LEOCCO=.FALSE.
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

USE mo_lnd_nwp_config,       ONLY: NTILES=>nsfc_subs

INTEGER(KIND=JPIM) :: NVTYPES, IVTYPES
REAL(KIND=JPRB)    :: ZLARGE , ZSNOW
REAL(KIND=JPRB)    :: RLHAERO, RLHAEROS

! Number of vegetation types
NVTYPES=20
IVTYPES = MAX(NVTYPES,20)


! Skin layer conductivity

ZSNOW=7._JPRB
ZLARGE=1.E10_JPRB
RLHAERO=15._JPRB
RLHAEROS=10._JPRB

! Unstable skin layer conductivity

IF(.NOT.ALLOCATED(RVLAMSK)) ALLOCATE(RVLAMSK(NTILES))
IF (NTILES >= 1) RVLAMSK(1)=ZLARGE           !  Water
IF (NTILES >= 2) RVLAMSK(2)=58._JPRB         !  Sea ice
IF (NTILES >= 3) RVLAMSK(3)=10._JPRB         !  Wet skin
IF (NTILES >= 4) RVLAMSK(4)=10._JPRB         !  Low veg. (Aerodyn.+Rad.+Soil)
IF (NTILES >= 5) RVLAMSK(5)=ZSNOW            !  Snow
IF (NTILES >= 5) RVLAMSK(6)=RLHAERO+5._JPRB  !  High veg. (Aerodyn.+Rad.)
IF (NTILES >= 7) RVLAMSK(7)=RLHAERO+5._JPRB  !  Snow under veg (Aerodyn.+Rad.)
IF (NTILES >= 8) RVLAMSK(8)=15._JPRB         !  Bare soil
IF (NTILES >= 9) RVLAMSK(9)=ZLARGE           !  LAKES = WATER  

! Stable skin layer conductivity

IF(.NOT.ALLOCATED(RVLAMSKS)) ALLOCATE(RVLAMSKS(NTILES))
IF (NTILES >= 1) RVLAMSKS(1)=ZLARGE            !  Water
IF (NTILES >= 2) RVLAMSKS(2)=58._JPRB          !  Sea ice
IF (NTILES >= 3) RVLAMSKS(3)=10._JPRB          !  Wet skin
IF (NTILES >= 4) RVLAMSKS(4)=10._JPRB          !  Low veg. (Aerodyn.+Rad.+Soil)
IF (NTILES >= 5) RVLAMSKS(5)=ZSNOW             !  Snow
IF (NTILES >= 5) RVLAMSKS(6)=RLHAEROS+5._JPRB  !  High veg. (Aerodyn.+Rad.)
IF (NTILES >= 7) RVLAMSKS(7)=RLHAEROS+5._JPRB  !  Snow under veg (Aerodyn.+Rad.+Snow
IF (NTILES >= 8) RVLAMSKS(8)=15._JPRB          !  Bare soil
IF (NTILES >= 9) RVLAMSKS(9)=ZLARGE            !  LAKES = WATER  

! Transmission of net solar rad. through vegetation

IF(.NOT.ALLOCATED(RVTRSR)) ALLOCATE(RVTRSR(NTILES))
IF (NTILES >= 1) RVTRSR(1)=0.00_JPRB     ! Ocean (SSR transmission does not apply)
IF (NTILES >= 2) RVTRSR(2)=0.00_JPRB     ! Ice   (SSR transmission does not apply)
IF (NTILES >= 3) RVTRSR(3)=0.05_JPRB     ! Wet skin (this is a compromise)
IF (NTILES >= 4) RVTRSR(4)=0.05_JPRB     ! Low vegetation
IF (NTILES >= 5) RVTRSR(5)=0.00_JPRB     ! Snow on vegetation (SSR transmission does not apply)
IF (NTILES >= 5) RVTRSR(6)=0.03_JPRB     ! High veg.
IF (NTILES >= 7) RVTRSR(7)=0.03_JPRB     ! Snow under veg.
IF (NTILES >= 8) RVTRSR(8)=0.00_JPRB     ! Bare soil (SSR transmission does not apply)
IF (NTILES >= 9) RVTRSR(9)=0.00_JPRB     ! LAKES (SSR transmission does not apply) 

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

IF(.NOT.ALLOCATED(RVZ0H)) ALLOCATE (RVZ0H(0:IVTYPES))
RVZ0H(1)=RVZ0M( 1)/10._JPRB     ! Crops, Mixed Farming
RVZ0H(2)=RVZ0M( 2)/10._JPRB     ! Short Grass
RVZ0H(3)=RVZ0M( 3)              ! Evergreen Needleleaf Trees
RVZ0H(4)=RVZ0M( 4)              ! Deciduous Needleleaf Trees
RVZ0H(5)=RVZ0M( 5)              ! Deciduous Broadleaf Trees
RVZ0H(6)=RVZ0M( 6)              ! Evergreen Broadleaf Trees
RVZ0H(7)=RVZ0M( 7)/10._JPRB     ! Tall Grass
RVZ0H(8)=RVZ0M( 8)/10._JPRB     ! Desert 
RVZ0H(9)=RVZ0M( 9)/10._JPRB     ! Tundra
RVZ0H(10)=RVZ0M(10)/10._JPRB    ! Irrigated Crops 
RVZ0H(11)=RVZ0M(11)/10._JPRB    ! Semidesert 
RVZ0H(12)=RVZ0M(12)/10._JPRB    ! Ice Caps and Glaciers 
RVZ0H(13)=RVZ0M(13)/10._JPRB    ! Bogs and Marshes
RVZ0H(14)=RVZ0M(14)/10._JPRB    ! Inland Water  
RVZ0H(15)=RVZ0M(15)/10._JPRB    ! Ocean 
RVZ0H(16)=RVZ0M(16)/10._JPRB    ! Evergreen Shrubs
RVZ0H(17)=RVZ0M(17)/10._JPRB    ! Deciduous Shrubs
RVZ0H(18)=RVZ0M(18)             ! Mixed Forest/woodland
RVZ0H(19)=RVZ0M(19)/10._JPRB    ! Interrupted Forest 
RVZ0H(20)=RVZ0M(20)/10._JPRB    ! Water and Land Mixtures 
RVZ0H(0)=RVZ0H(8)        

END SUBROUTINE SUSVEG


!------------------------------------------------------------------------------


SUBROUTINE ABORT_SURF(text)
  USE mo_io_units, ONLY: nerr
  USE mo_mpi, ONLY: p_abort
  CHARACTER (len=*), INTENT(in) :: text
  WRITE (nerr,'(a)')  TRIM(text)
  CALL p_abort
END SUBROUTINE ABORT_SURF



END MODULE mo_edmf_param
