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

! * ECMWF Single Column Model:
LOGICAL :: LSCMEC
LOGICAL :: LSFCFLX
REAL(KIND=JPRB) :: REXTSHF
REAL(KIND=JPRB) :: REXTLHF


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

  PUBLIC :: FOEALFA  ,FOEEWM   ,FOEDEM   ,FOELDCPM, &
          & FOEALFCU ,FOEEWMCU ,FOEDEMCU ,FOELDCPMCU
  PUBLIC :: suct0, su0phy, susekf

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


END MODULE mo_edmf_param
