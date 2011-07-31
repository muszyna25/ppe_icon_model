!>
!!
!! Ritter-Geleyn radiation parameterization from GME (provided by Bodo Ritter, DWD, Offenbach)
!!
!! @author Thorsten Reinhardt, AGeoBw, Offenbach (Implementation into ICON)
!! @author Bodo Ritter, DWD, Offenbach (original code)
!!
!! @par Revision History
!! Initial release by Thorsten Reinhardt, AGeoBw, 2010-12-02
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
MODULE mo_radiation_rg

  USE mo_kind, ONLY : wp

  USE mo_exception,          ONLY: message, message_text  

  USE mo_physical_constants, ONLY: g     => grav  , & !! acceleration due to gravity
    &                              r_d   => rd, &
    &                              r_v   => rv, &
    &                              Rvd_m_o =>  vtmpc1, &
    &                              psig  => stbo

  USE mo_radiation_rg_par, ONLY : &
  &                               jpsol,jpther,jpspec, &
  &                               ncgas, nfast, &
  &                               coai, cobi, coali, cobti, &
  &                               zlwe, zlww, zlwg, zlwemn, zlwemx, &
  &                               zaea, zaes, zaef, zaeg, &
  &                               ziwe, ziwemn, ziwemx, ziww, ziwg!, &
    
  IMPLICIT NONE

  PRIVATE

  PUBLIC  :: fesft

!!$  CHARACTER(len=*), PARAMETER :: version = '$Id$'

CONTAINS

      SUBROUTINE fesft(                                            &
       pti   ,pdp  ,pclc_in ,pqv  ,pqvs  ,pqcwc,pqiwc,pduco2,pduo3,&
       paeq1 ,paeq2,paeq3,paeq4,paeq5,                             &
       papre_in ,psmu0,palso,palth,                                &
       psct ,                                                      &
       kig1s ,kig1e,ki3s,ki3e,                                     &
       ki1sc ,ki1ec,                                               &
       lsolar,lthermal,lcrf ,                                      &
       &       pflt  ,pfls  )
!!$       pfltf,pflsf,pflpar,                            &
!!$       pflsp,pflsd,pflsu       )


! 
!**** *fesft* - organizes the radiative transfer calculations
 
!     Author:    B.Ritter April 1996 based on earlier versions of
!                                    the same routine
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
! 1.14       2002/01/16 D. Majewski; H. Frank
!  Rename pqlwc to pqcwc and pdulwc to pducwc and prholwc to prhocwc
!  Introduce KIND-notation for all reals with real2kind.pl
!  Include 'gme_comconst.h' to set R_d, R_v. Replace ZRVDM1 by Rvd_m_o,
!  RD by R_d.
! V2_23        2010/04/27 J. Helmert
!  Add output of radiation flux components         
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
!     Purpose and Method
!     ------------------
!
!     This routine organizes the radiative transfer calculations by
!     calling a set of dedicated routines for the calculation of
!     basic optical properties (opt_th/opt_so), the derivation of
!     layer coefficients (coe_th/coe_so) for an implicit delta-two-
!     stream scheme (cf.Ritter and Geleyn, 1992) and the inversion
!     (inv_th/inv_so) of the corresponding system matrix. These 
!     operations are performed seperately for thermal and solar parts
!     of the spectrum and are embedded in loops over the various
!     spectral intervals. Within each interval, a data-controlled
!     decision is taken, whether the so-called ESFT or FESFT approach
!     is used for the handling of gaseous absorption (cf. Ritter and
!     Geleyn, 1992).
!     Controlled by the logical input variable LCRF, the calculation
!     of radiative fluxes in cloud-free conditions can be done in
!     addition to the results for the given atmospheric cloud structure.
!     (not implemented yet)
!     Before the actual flux calculation starts, some preliminary steps
!     provide utility arrays which are applicable to all spectral inter-
!     vals (e.g. cloud geometry factors, integrated layer water content,
!     etc.)
!==============================================================================

      INTEGER, INTENT(in) :: kig1s,kig1e
!     horizontal dimensions for arrays defined on the computational domain
      INTEGER , INTENT(in) ::ki1sc,ki1ec
!     vertical array limits
      INTEGER, INTENT(in) :: ki3s ,ki3e

      ! Temperature at layer boundaries, pressure thickness of layers,
      ! cloud cover in each layer, water vapour and saturation water
      ! vapour mixing ratio, liquid water and ice mixing ratio
      ! layer CO2 content and O3 content
      REAL(wp), INTENT(in) :: pti   (kig1s:kig1e,ki3s:ki3e+1)            ! (K)
      REAL(wp), INTENT(in) :: pdp   (kig1s:kig1e,ki3s:ki3e)              ! (Pa)
      REAL(wp), INTENT(in) :: pclc_in(kig1s:kig1e,ki3s:ki3e )! (1) clc input 
      REAL(wp), INTENT(in) :: pqv   (kig1s:kig1e,ki3s:ki3e)  ! (kg/kg)
      REAL(wp), INTENT(in) :: pduco2(kig1s:kig1e,ki3s:ki3e)              ! (Pa CO2)  
      REAL(wp), INTENT(in) :: pduo3 (kig1s:kig1e,ki3s:ki3e)              ! (Pa O3)  
!!$      REAL(wp) :: pti   (:,:)       ! (K)
!!$      REAL(wp) :: pdp   (:,:)         ! (Pa)
!!$      REAL(wp) :: pclc  (:,:)        ! (1)
!!$      REAL(wp) :: pqv   (:,:)         ! (kg/kg)
!!$      REAL(wp) :: pduco2(:,:)         ! (Pa CO2)
!!$      REAL(wp) :: pduo3 (:,:)         ! (Pa O3)  
!


      
!     The following three arrays are local in *physics*
!
      REAL(wp), INTENT(in) :: pqvs  (kig1s:kig1e,ki3s:ki3e)         ! (kg/kg) 
      REAL(wp), INTENT(in) :: pqcwc (kig1s:kig1e,ki3s:ki3e)         ! (kg/kg) 
      REAL(wp), INTENT(in) :: pqiwc (kig1s:kig1e,ki3s:ki3e)         ! (kg/kg)
!!$      REAL(wp) :: pqvs  (:,:)         ! (kg/kg) 
!!$      REAL(wp) :: pqcwc (:,:)         ! (kg/kg) 
!!$      REAL(wp) :: pqiwc (:,:)         ! (kg/kg) 
      

      ! Aerorsole optical depth at 0.55 micrometer for five types
      REAL(wp), INTENT(in) :: paeq1  (kig1s:kig1e,ki3s:ki3e)         ! (1)
      REAL(wp), INTENT(in) :: paeq2  (kig1s:kig1e,ki3s:ki3e)         ! (1)
      REAL(wp), INTENT(in) :: paeq3  (kig1s:kig1e,ki3s:ki3e)         ! (1)
      REAL(wp), INTENT(in) :: paeq4  (kig1s:kig1e,ki3s:ki3e)         ! (1)
      REAL(wp), INTENT(in) :: paeq5  (kig1s:kig1e,ki3s:ki3e)         ! (1)
!!$      REAL(wp) :: paeq1  (:,:)         ! (1)
!!$      REAL(wp) :: paeq2  (:,:)         ! (1)
!!$      REAL(wp) :: paeq3  (:,:)         ! (1)
!!$      REAL(wp) :: paeq4  (:,:)         ! (1)
!!$      REAL(wp) :: paeq5  (:,:)         ! (1)
      

      ! Surface pressure, cosine of zenith angle, thermal and solar
      ! surface albedo
      REAL(wp), INTENT(in) :: papre_in(kig1s:kig1e)              ! (Pa) input surface pressure
      REAL(wp), INTENT(in) :: psmu0  (kig1s:kig1e)               ! (1)  
      REAL(wp), INTENT(in) :: palth  (kig1s:kig1e)               ! (1)  
      REAL(wp), INTENT(in) :: palso  (kig1s:kig1e)               ! (1)
!!$      REAL(wp) :: papre  (:)               ! (Pa)
!!$      REAL(wp) :: psmu0  (:)               ! (1)  
!!$      REAL(wp) :: palth  (:)               ! (1)  
!!$      REAL(wp) :: palso  (:)               ! (1)  
      
!
      REAL(wp), INTENT(in) :: psct          ! solar constant (at time of year)
!

      LOGICAL, INTENT(in) :: lsolar     ! control switch for solar calculations
      LOGICAL, INTENT(in) :: lthermal   ! control switch for thermal calculations
      LOGICAL, INTENT(in) :: lcrf       ! control switch for cloud-free calcul.


!     Output data

      ! Thermal and solar radiative fluxes at each layer boundary
      ! dito for cloud-free conditions (TOA and surface only)
      ! surface flux of photosynthetic active radiation 
      REAL(wp), INTENT(out) :: pflt   (kig1s:kig1e,ki3s:ki3e+1) ! (W/m**2)
      REAL(wp), INTENT(out) :: pfls   (kig1s:kig1e,ki3s:ki3e+1) ! (W/m**2)
!!$      REAL(wp) :: pflt   (:,:) ! (W/m**2)
!!$      REAL(wp) :: pfls   (:,:) ! (W/m**2)
      
      REAL(wp) :: pfltf  (kig1s:kig1e,2)           ! (W/m**2)
      REAL(wp) :: pflsf  (kig1s:kig1e,2)           ! (W/m**2)

      ! Solar and thermal radiation flux components
      REAL(wp) :: pflpar (kig1s:kig1e)             ! (W/m**2)
      REAL(wp) :: pflsp  (kig1s:kig1e)             ! (W/m**2)
      REAL(wp) :: pflsd  (kig1s:kig1e)             ! (W/m**2)
      REAL(wp) :: pflsu  (kig1s:kig1e)             ! (W/m**2)

!     Arrays local to *fesft* or required for communication with
!     subroutines called from *fesft*

      REAL(wp) :: pclc  (kig1s:kig1e,ki3s:ki3e ) ! (1)  internal clc ( 0. < clc < 1. )
      REAL(wp) :: papre(kig1s:kig1e)             ! (Pa) internal surface pressure

!     'Grey' and gaseous fluxes for individual spectral intervals
      REAL(wp) :: zflux  (kig1s:kig1e,ki3s:ki3e+1) ! (W/m**2)
      REAL(wp) :: zfluxi (kig1s:kig1e,ki3s:ki3e+1) ! 1./(W/m**2)
      REAL(wp) :: zfluxu (kig1s:kig1e,ki3s:ki3e+1) ! (W/m**2)
      REAL(wp) :: zfluxui(kig1s:kig1e,ki3s:ki3e+1) ! 1./(W/m**2)
      REAL(wp) :: zfluxd (kig1s:kig1e,ki3s:ki3e+1) ! (W/m**2)
      REAL(wp) :: zfluxdi(kig1s:kig1e,ki3s:ki3e+1) ! 1./(W/m**2)
      REAL(wp) :: zfgas  (kig1s:kig1e,ki3s:ki3e+1) ! (W/m**2)
      REAL(wp) :: zfgasu (kig1s:kig1e,ki3s:ki3e+1) ! (W/m**2)
      REAL(wp) :: zfgasd (kig1s:kig1e,ki3s:ki3e+1) ! (W/m**2)

!     Black body radiation at layer boundaries
      REAL(wp) :: pbbr   (kig1s:kig1e,ki3s:ki3e+1) ! (W/m**2)
!     Solar flux at TOA                              
      REAL(wp) :: pflpt  (kig1s:kig1e)               ! 
!     Solar surface albedo for parallel radiation
      REAL(wp) :: palp   (kig1s:kig1e)               ! 
!     Inverse of cosine of zenith angle               
      REAL(wp) :: pqsmu0 (kig1s:kig1e)               ! 
!     Logarithm of layer mean temperature and pressure
      REAL(wp) :: palogt (kig1s:kig1e,ki3s:ki3e) ! ln T
      REAL(wp) :: palogp (kig1s:kig1e,ki3s:ki3e) ! ln p
!     pressure at one level                             
      REAL(wp) :: papra  (kig1s:kig1e)               ! (Pa)
!     layer water vapour content (cloudy and cloud-free)
      REAL(wp) :: pduh2oc(kig1s:kig1e,ki3s:ki3e) ! (Pa H2O-vapour) 
      REAL(wp) :: pduh2of(kig1s:kig1e,ki3s:ki3e) ! (Pa H2O-vapour) 
!     cloudy part of layer liquid water and ice content density
      REAL(wp) :: pducwc (kig1s:kig1e,ki3s:ki3e) ! (Pa H2O-liquid)
      REAL(wp) :: pduiwc (kig1s:kig1e,ki3s:ki3e) ! (Pa H2O-ice)
      REAL(wp) :: prhocwc(kig1s:kig1e,ki3s:ki3e) ! (kg/m**3)
      REAL(wp) :: prhoiwc(kig1s:kig1e,ki3s:ki3e) ! (kg/m**3)
!     water vapour e-type contribution (cloudy and cloud-free)
      REAL(wp) :: zduetpc(kig1s:kig1e,ki3s:ki3e)   !
      REAL(wp) :: zduetpf(kig1s:kig1e,ki3s:ki3e)   !
      
!     layer mean temperature, water vapour mixing ratio, utility arrays

      REAL(wp) :: ztm    (kig1s:kig1e)               !
      REAL(wp) :: zwv    
      REAL(wp) :: zcpo   (kig1s:kig1e)               !
      REAL(wp) :: zcpn   (kig1s:kig1e)               !
      REAL(wp) :: zcmo   (kig1s:kig1e)               !
      REAL(wp) :: zcmn   (kig1s:kig1e)               !

!     Output data from opt_th/opt_so

      REAL(wp) :: podac (kig1s:kig1e,ki3s:ki3e) ! absorption optical depth
      REAL(wp) :: podaf (kig1s:kig1e,ki3s:ki3e) ! in cloudy and free part
      REAL(wp) :: podsc (kig1s:kig1e,ki3s:ki3e) ! scattering optical depth
      REAL(wp) :: podsf (kig1s:kig1e,ki3s:ki3e) ! in cloudy and free part
      REAL(wp) :: pbsfc (kig1s:kig1e,ki3s:ki3e) ! backscattering fraction 
      REAL(wp) :: pbsff (kig1s:kig1e,ki3s:ki3e) ! in cloudy and free part
      REAL(wp) :: pusfc (kig1s:kig1e,ki3s:ki3e) ! upscattering   fraction 
      REAL(wp) :: pusff (kig1s:kig1e,ki3s:ki3e) ! in cloudy and free part


!     cloud geometry factors

      REAL(wp) :: pca1  (kig1s:kig1e,ki3s:ki3e) !
      REAL(wp) :: pcb1  (kig1s:kig1e,ki3s:ki3e) !
      REAL(wp) :: pcc1  (kig1s:kig1e,ki3s:ki3e) !
      REAL(wp) :: pcd1  (kig1s:kig1e,ki3s:ki3e) !
      REAL(wp) :: pca2  (kig1s:kig1e,ki3s:ki3e) !
      REAL(wp) :: pcb2  (kig1s:kig1e,ki3s:ki3e) !
      REAL(wp) :: pcc2  (kig1s:kig1e,ki3s:ki3e) !
      REAL(wp) :: pcd2  (kig1s:kig1e,ki3s:ki3e) !

!     fluxes calculated in inv_th/inv_so
      REAL(wp) :: pflfd  (kig1s:kig1e,ki3s:ki3e+1) ! (W/m**2)
      REAL(wp) :: pflfu  (kig1s:kig1e,ki3s:ki3e+1) ! (W/m**2)
      REAL(wp) :: pflfp  (kig1s:kig1e,ki3s:ki3e+1) ! (W/m**2)
      REAL(wp) :: pflcd  (kig1s:kig1e,ki3s:ki3e+1) ! (W/m**2)
      REAL(wp) :: pflcu  (kig1s:kig1e,ki3s:ki3e+1) ! (W/m**2)
      REAL(wp) :: pflcp  (kig1s:kig1e,ki3s:ki3e+1) ! (W/m**2)
 
!     Local variables for data arrays and execution control
 
      REAL(wp) :: solant (JPSOL)
      REAL(wp) :: planck (3,jpther)
      REAL(wp) :: zketyp (jpther) ,ztetyp (jpther)
      REAL(wp) :: zketypR(jpther) ,ztetypR(jpther)
      REAL(wp) :: zketypA(jpther) ,ztetypA(jpther)

!     REAL(wp) :: RD,RV,ZRVD,ZRVDM1,ztvm,zet,zteref,zzetyp,zzet,zaiprod
      REAL(wp) :: ZRVD,ztvm,zet,zteref,zzetyp,zzet,zaiprod
      REAL(wp) :: zcoai,zcobi
      REAL(wp) :: zepai    ! Could be used to save computing time for
                    ! 'unimportant' gaseous absorption coefficients
      INTEGER icgas  (3)
      INTEGER ifirst, icrf, igase, igasm1, igasz
      INTEGER j1,j3          ! loop indices over spatial dimensions
      INTEGER jc,jh2o,jco2,jo3  ! loop indices over gaseous coefficients
      INTEGER jspec,jspect      ! loop indices over spectrum
      INTEGER jg,jjg            ! loop indices over gases      
      INTEGER icc,ih2o,ico2,io3 ! loop limit for gaseous absorption
!
!     Debug switches for higher level subroutines
!
      LOGICAL, PARAMETER :: ldebug         = .FALSE. ! debug control switch
      LOGICAL, PARAMETER :: ldebug_th      = .FALSE. ! debug control switch for thermal     
      LOGICAL, PARAMETER :: ldebug_so      = .FALSE. ! debug control switch for solar       
      LOGICAL, PARAMETER :: ldebug_opt_th  = .FALSE. ! debug control switch for opt_th      
      LOGICAL, PARAMETER :: ldebug_opt_so  = .FALSE. ! debug control switch for opt_so      
      LOGICAL, PARAMETER :: ldebug_inv_th  = .FALSE. ! debug control switch for inv_th      
      LOGICAL, PARAMETER :: ldebug_inv_so  = .FALSE. ! debug control switch for inv_so      

!     Security variables   
 
      REAL(wp) :: zeplw  ! Minimum liquid water content to avoid numerical 
                                  ! problems in *opt_th* and *opt_so*

      REAL(wp) :: zepflx ! Minimum 'grey' flux to avoid 1./0.

 !     REAL(wp) :: pflfpmn,pflfdmn,pflfumn,pflcpmn,pflcdmn,pflcumn
 !KF unusedvar     REAL(wp) :: pflfpmx,pflfdmx,pflfumx,pflcpmx,pflcdmx,pflcumx

      INTEGER j1b        ! index of debug-point
 
!     data statements 
 
      DATA IFIRST/0/
      DATA j1b /1/
 
!      INCLUDE "gme_dat_sol.h"
!     Fraction of solar energy at TOA contained in individual
!     solar spectral intervals based on data from LABS and
!     NECKEL (1970/1984):
 
      DATA SOLANT/ 0.12888167_wp, 0.41683156_wp, 0.45428677_wp/
       
!      INCLUDE "gme_dat_ther.h"
!     *PLANCK* : coefficients for the description of the fraction of
!                the total black body radiation contained in an thermal
!                spectral interval as a function (2.order polynomial)
!                of temperature:
!
!                F(T) =  PLANCK(1,ISPEC)
!                      + PLANCK(2,ISPEC) * T
!                      + PLANCK(3,ISPEC) * T**2
!
      DATA PLANCK/                                             &
     &        0.157656E+01_wp, -0.711486E-02_wp,  0.908220E-05_wp, &
     &       -0.761337E-01_wp,  0.339014E-02_wp, -0.703246E-05_wp, &
     &       -0.353624E+00_wp,  0.321131E-02_wp, -0.472513E-05_wp, &
     &       -0.180726E+00_wp,  0.148131E-02_wp, -0.195189E-05_wp, &
     &        0.343422E-01_wp, -0.971936E-03_wp,  0.463714E-05_wp /
 
!     *ZKETYP*   e-type continuum-coefficient for all spectral inter-
!                vals (PA (H2O)**-2) at 296 K
!                *R*  following ROBERTS ET AL. 1976
!                *A*  implicitly derived from the AFGL spectral data
!     *ZTETYP*   constant for the temperature dependancy of e-type   
!                absorption for all intervals
!     *ZTEREF  : reference temperaure 
!
      DATA ZKETYPR/ 0.0_wp , 0.418E-06_wp , 0.859E-06_wp ,       &
     &        0.594E-06_wp , 0.767E-07_wp /   
!    & 14.900E-06 , 3.530E-06 , 0.859E-06 , 0.594E-06 , 0.433E-06 /
 
      DATA ZTETYPR/ 0._wp ,    0._wp, 1800._wp,                  &
     &           1800._wp , 1800._wp /
 
      DATA ZKETYPA/ 0.8426E-05_wp, 0.8982E-06_wp, 0.5489E-06_wp, &
     &              0.4743E-06_wp, 0.7040E-06_wp/
 
      DATA ZTETYPA/ 1365.55_wp, 1544.38_wp, 1699.06_wp,          &
     &              1724.39_wp, 1668.94_wp   /
 
      DATA ZTEREF /  296.0_wp/
 


      ZRVD   = R_v/R_d
 
!      zeplw  = 1.0E-12_wp
      zepflx = 1.0E-12_wp
      zepflx = 1.0E-08_wp
      zepai  = 0.0_wp
 
      ICRF  = 0

      IF (ldebug) THEN
      PRINT *,' **** FESFT  ******************************'
      PRINT *,' **** debug point : ',j1b
      PRINT *,' '
      PRINT *,' **** kig1s,kig1e : ',kig1s,kig1e
      PRINT *,' '
      PRINT *,' **** ki3s,ki3e   : ',ki3s,ki3e
      PRINT *,' '
      PRINT *,' solar   calculations: ',lsolar 
      PRINT *,' thermal calculations: ',lthermal
        WRITE(message_text,'(a,2E15.5)') ' max/min pti(:,:)   = ',&
             & MAXVAL(pti(:,:) ), MINVAL( pti(:,:) )
        CALL message('', TRIM(message_text))

        WRITE(message_text,'(a,2E15.5)') ' max/min pdp(:,:)   = ',&
             & MAXVAL(pdp(:,:) ), MINVAL( pdp(:,:) )
        CALL message('', TRIM(message_text))

        WRITE(message_text,'(a,2E15.5)') ' max/min pqv(:,:)   = ',&
             & MAXVAL(pqv(:,:) ), MINVAL( pqv(:,:) )
        CALL message('', TRIM(message_text))

        WRITE(message_text,'(a,2E15.5)') ' max/min palso(:)   = ',&
             & MAXVAL(palso(:) ), MINVAL( palso(:) )
        CALL message('', TRIM(message_text))

        WRITE(message_text,'(a,2E15.5)') ' max/min palth(:)   = ',&
             & MAXVAL(palth(:) ), MINVAL( palth(:) )
        CALL message('', TRIM(message_text))

        WRITE(message_text,'(a,2E15.5)') ' max/min papre(:)   = ',&
             & MAXVAL(papre(:) ), MINVAL( papre(:) )
        CALL message('', TRIM(message_text))   

        WRITE(message_text,'(a,2E15.5)') ' psct = ',&
             & psct
        CALL message('', TRIM(message_text))
      END IF
      
      DO j3 = ki3s,ki3e
        DO j1 = ki1sc,ki1ec
          pclc(j1,j3) = MAX(pclc_in(j1,j3),1.E-9_wp)       !to avoid cloud fraction = 0.0
          pclc(j1,j3) = MIN(pclc   (j1,j3),1._wp-1.E-9_wp) !to avoid cloud fraction = 1.0
        END DO
      END DO
      
!     Preset output arrays
 
      DO j3 = ki3s,ki3e+1
          DO j1 = ki1sc,ki1ec
            pflt  (j1,j3)=0._wp
            pfls  (j1,j3)=0._wp
          END DO
      END DO     ! End of vertical loop
  
      DO j1 = ki1sc,ki1ec
        pfltf (j1,1)=0._wp
        pfltf (j1,2)=0._wp
        pflsf (j1,1)=0._wp
        pflsf (j1,2)=0._wp
        pflpar(j1 ) =0._wp
        pflsp (j1 ) =0._wp
        pflsd (j1 ) =0._wp
        pflsu (j1 ) =0._wp
        papre (j1 ) = papre_in(j1)
      END DO

!      Choice of e-type-absorption and temperature correction coefficients
       DO jspec=1,jpther     
       zketyp(jspec) = zketypA(jspec)
       ztetyp(jspec) = ztetypA(jspec)
       END DO

!     cloud geometry factors
 
      IF (ldebug) THEN
      PRINT *,' FESFT   some input   at debug point'
      PRINT *,' ==================================='
      PRINT *,' clc        qv      qvs     pti        '
      DO j3=ki3s,ki3e
      PRINT 1001,pclc (j1b,j3),pqv(j1b,j3),pqvs(j1b,j3),pti(j1b,j3)
      END DO
!     STOP 'test'
      END IF

!     first part for top layer
        DO j1 = ki1sc,ki1ec
        pca2(j1,ki3s)= 1._wp-pclc(j1,ki3s)
        pcd2(j1,ki3s)= 1._wp
        zcpn(j1)  = MAX(pclc(j1,ki3s),pclc(j1,ki3s+1))
        zcmn(j1)  = MIN(pclc(j1,ki3s),pclc(j1,ki3s+1))
        pca2(j1,ki3s+1)= (1._wp-zcpn(j1))/pca2(j1,ki3s)
        pcd2(j1,ki3s+1)= zcmn(j1)/pclc(j1,ki3s)
        END DO
 
!     first part for inner layers
 
      DO j3=ki3s+1,ki3e-1
          DO j1 = ki1sc,ki1ec
          zcpo(j1)=zcpn(j1)
          zcmo(j1)=zcmn(j1)
          END DO
          DO j1 = ki1sc,ki1ec
          zcpn(j1)      = MAX (pclc(j1,j3),pclc(j1,j3+1))
          zcmn(j1)      = MIN (pclc(j1,j3),pclc(j1,j3+1))
          pca2(j1,j3+1) = (1._wp-zcpn(j1))/(1._wp-pclc(j1,j3))
          pca1(j1,j3-1) = (1._wp-zcpo(j1))/(1._wp-pclc(j1,j3))
          pcd2(j1,j3+1) = zcmn(j1)/pclc(j1,j3)
          pcd1(j1,j3-1) = zcmo(j1)/pclc(j1,j3)
          END DO
      END DO
 
!     first part for lowest layer
 
      DO j1 = ki1sc,ki1ec
      pca1(j1,ki3e-1) = (1._wp-zcpn(j1))/(1._wp-pclc(j1,ki3e))
      pcd1(j1,ki3e-1) = zcmn(j1)/pclc(j1,ki3e)
      pca1(j1,ki3e)  = 1._wp
      pcd1(j1,ki3e)  = 1._wp
      END DO
 
!     second part of geometry factors

      DO j3 = ki3s,ki3e
          DO j1 = ki1sc,ki1ec
          pcb1(j1,j3)=1._wp-pca1(j1,j3)
          pcc1(j1,j3)=1._wp-pcd1(j1,j3)
          pcb2(j1,j3)=1._wp-pca2(j1,j3)
          pcc2(j1,j3)=1._wp-pcd2(j1,j3)
          END DO
      END DO
 
!     Optically relevant layer constituents
!      (Note: CO2 and O3 amounts are provided by calling routine)
        DO j1 = ki1sc,ki1ec
        papra(j1) = papre(j1)   ! surface pressure
        END DO
 
!     water vapour, liquid water and ice content, logarithm of layer
!     mean temperature and pressure, absorber amount for e-type 
!     absorption
 
      DO j3=ki3e,ki3s,-1    ! Bottom to top
          DO j1 = ki1sc,ki1ec
          ztm   (j1)    = 0.5_wp*(pti(j1,j3)+pti(j1,j3+1))
          papra (j1)    =   papra(j1)-0.5_wp*pdp(j1,j3)
          palogt(j1,j3) = LOG  (ztm  (j1))
          palogp(j1,j3) = LOG  (papra(j1))
          END DO
 
!     cloud-free:  water vapour and e-type absorber amount
          DO j1 = ki1sc,ki1ec
          zwv            = MAX(                                    &
          (pqv(j1,j3)-pclc(j1,j3)*pqvs(j1,j3))            &
             /(1._wp-pclc(j1,j3))       , 0._wp)
          pduh2of(j1,j3) = pdp(j1,j3)*zwv
          zduetpf(j1,j3) = pduh2of(j1,j3)*pduh2of(j1,j3)  &
                             *papra(j1)*zrvd/pdp(j1,j3)
          END DO
 
!     cloudy:  water vapour, e-type absorber amount, liquid water and ice
          DO j1 = ki1sc,ki1ec
      pducwc(j1,j3)  = pdp(j1,j3)*(pqcwc(j1,j3)/pclc(j1,j3))
      pducwc(j1,j3)  = MAX(pducwc(j1,j3),0._wp)
      pduiwc(j1,j3)  = pdp(j1,j3)*(pqiwc(j1,j3)/pclc(j1,j3))
      pduiwc(j1,j3)  = MAX(pduiwc(j1,j3),0._wp)
      pduh2oc(j1,j3) = pdp(j1,j3)*pqvs(j1,j3)
      zduetpc(j1,j3) = pduh2oc(j1,j3)*pduh2oc(j1,j3)           &
                          *papra(j1)*zrvd/pdp(j1,j3)
!     ztvm              = ztm(j1)*(1._wp+zrvdm1*pqvs(j1,j3))
      ztvm              = ztm(j1)*(1._wp+Rvd_m_o*pqvs(j1,j3))
      prhocwc(j1,j3) =(pqcwc(j1,j3)/pclc(j1,j3))*papra(j1)  &
                               /(R_d*ztvm)
      prhoiwc(j1,j3) =(pqiwc(j1,j3)/pclc(j1,j3))*papra(j1)  &
                               /(R_d*ztvm)
!     Secure minium for ice density for use in empirical function with LOG !
      prhoiwc(j1,j3) = MAX (prhoiwc(j1,j3),1.E-06_wp) 
      
      papra(j1) = papra(j1)-0.5_wp*pdp(j1,j3)
          END DO
      END DO     ! End of vertical loop

!     Identify *papre* with top of model pressure (for Rayleigh scattering)
 
        DO j1 = ki1sc,ki1ec
        papre(j1) = papra(j1)
        END DO
 
    1 CONTINUE  ! Address for backward jump to perform cloud-free calculations

!     Thermal radiative flux calculations 
!     ===================================

      IF (lthermal) THEN

!     Loop over thermal spectral intervals

! ===============================================================
      DO 2 jspec=jpsol+1,jpspec
! ===============================================================
      
      jspect = jspec - jpsol
 
!     Black body radiation at layer boundaries in spectral interval
 
      DO j3 = ki3s,ki3e+1
          DO j1 = ki1sc,ki1ec
          pbbr(j1,j3)= (           planck(1,jspect) +    &
                       pti(j1,j3)*(planck(2,jspect) +    &
                       pti(j1,j3)* planck(3,jspect)  ) ) &
                         *psig*(pti(j1,j3)**2)**2
          END DO
      END DO     ! End of vertical loop
 
!     Optical properties of non-gaseous constituents        
 
      IF (ldebug) THEN
      PRINT *,' FESFT    Call to opt_th for jspec: ',jspec
      PRINT *,' ========================================='
      PRINT *,' prhocwc   pducwc    prhoiwc   pduiwc '
      DO j3=ki3s,ki3e
      PRINT 1001,prhocwc(j1b,j3),pducwc(j1b,j3),  &
                 prhoiwc(j1b,j3),pduiwc(j1b,j3)
      END DO
 1001 format (1x,4f10.5)
      END IF
      CALL opt_th(prhocwc,pducwc,prhoiwc,pduiwc,                 &
                  paeq1  ,paeq2 ,paeq3  ,paeq4 , paeq5 ,         &
                  jspec  ,kig1s, kig1e ,ki1sc ,ki1ec   ,         &
                    ki3s ,ki3e  ,ldebug_opt_th,                  &
                  podac  ,podaf ,podsc  ,podsf , pbsfc ,pbsff  )
 
!     Addition of e-type contribution
 
      IF (zketyp(jspect).NE.0._wp) THEN
      zet   = 1._wp/EXP(ztetyp(jspect)/zteref)
      DO j3 = ki3s,ki3e
          DO j1 = ki1sc,ki1ec
         ztm(j1) = 0.5_wp*(pti(j1,j3)+pti(j1,j3+1))
         zzetyp     = zet*EXP(ztetyp(jspect)/ztm(j1))*zketyp(jspect)
         podaf(j1,j3)=podaf(j1,j3)+zduetpf(j1,j3)*zzetyp
         podac(j1,j3)=podac(j1,j3)+zduetpc(j1,j3)*zzetyp
          END DO
      END DO     ! End of vertical loop
         END IF
 
!     Selection of ESFT or FESFT method for interval considered
!--------------------------------------------------------------
      IF (nfast(jspec).EQ.0) THEN           ! ESFT method
!--------------------------------------------------------------
      DO j3 = ki3s,ki3e+1
          DO j1 = ki1sc,ki1ec
          zflux(j1,j3)=0.0_wp  ! Preset flux in spectral interval
          END DO
      END DO     ! End of vertical loop
 
!     Loop over various absorption coefficients of each of the three gases
 
      IH2O = NCGAS(JSPEC,1)
      ICO2 = NCGAS(JSPEC,2)
      IO3  = NCGAS(JSPEC,3)
!
      DO jh2o=1,IH2O    ! Loop over H2O coefficients
      DO jco2=1,ICO2    ! Loop over CO2 coefficients
      DO jo3=1,IO3      ! Loop over O3  coefficients
!
      zaiprod = coai(jh2o,jspec,1)*coai(jco2,jspec,2)*coai(jo3,jspec,3)
!
      IF (icrf.EQ.0) THEN      ! partially cloudy atmosphere
      CALL inv_th (                                                    &
        pclc   ,pca1   ,pca2  ,pcb1  ,pcb2  ,pcc1  ,pcc2 ,pcd1 ,pcd2 , &
        pduh2oc,pduh2of,pduco2,pduo3 ,palogp,palogt,                   &
        podsc  ,podsf  ,podac ,podaf ,pbsfc ,pbsff ,                   &
        pbbr   ,palth,                                                 &
        jspec  ,jh2o   ,jco2  ,jo3   ,                                 &
        kig1s  ,kig1e  ,              ki3s ,ki3e ,                     &
        ki1sc  ,ki1ec  ,                                               &
        ldebug_inv_th ,                                                &
        pflcu  ,pflfu  ,pflcd ,pflfd )
      ELSE                     ! 'cloud-free' atmosphere    
      PRINT *,' CRF not yet implemented'
!     CALL INVLONF(JSPEC,jh2o,jco2,jo3,IRLON,IRLEV,IRLVP1,            &
!    & PBBR,PALTE,                                                    &
!    & pduh2of,PDUCO2,PDUO3,palogp,palogt,PODSF,                      &
!    & PODAF,PBSFF,                                                   &
!    & pflfu,pflfd )
      END IF
 
!     Incrementation of flux in spectral interval

      DO j3 = ki3s,ki3e+1
          DO j1 = ki1sc,ki1ec
          zflux(j1,j3) = zflux(j1,j3)+zaiprod*         &
                           (pflfu(j1,j3)+pflcu(j1,j3)  &
                           -pflfd(j1,j3)-pflcd(j1,j3))
          END DO
      END DO     ! End of vertical loop

      END DO           ! Loop over O3 absorption coefficients
      END DO           ! Loop over CO2 absorption coefficients
      END DO           ! Loop over H2O absorption coefficients
!--------------------------------------------------------------
      ELSE                                  ! FESFT method
!--------------------------------------------------------------
      igase = 0
      DO jg=1,3
      icgas (jg) = 0
      IF (NCGAS(jspec,jg).GT.1) THEN
      igase     = igase + 1
      END IF
      END DO
      igasm1    = igase -1
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
               IF (igase.le.1) THEN  !(no 'grey' fluxes required)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
               DO j3 = ki3s,ki3e+1
                   DO j1 = ki1sc,ki1ec
                   zfluxi(j1,j3) = 1._wp
                   END DO
               END DO
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
               ELSE   ! more than 1 gas --> 'grey' fluxes required
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF (icrf.EQ.0) THEN    ! partially cloudy atmosphere
      CALL inv_th (                                                     &
        pclc   ,pca1   ,pca2  ,pcb1  ,pcb2  ,pcc1  ,pcc2 ,pcd1 ,pcd2 ,  &
        pduh2oc,pduh2of,pduco2,pduo3 ,palogp,palogt,                    &
        podsc  ,podsf  ,podac ,podaf ,pbsfc ,pbsff ,                    &
        pbbr   ,palth,                                                  &
        jspec  ,0      ,0     ,0     ,                                  &
        kig1s  ,kig1e  ,              ki3s ,ki3e ,                      &
        ki1sc  ,ki1ec  ,                                                &
        ldebug_inv_th ,                                                 &
        pflcu  ,pflfu  ,pflcd ,pflfd )
      ELSE                     ! 'cloud-free' atmosphere    
      PRINT *,' CRF not yet implemented'
!     CALL INVLONF(JSPEC,0,0,0,IRLON,IRLEV,IRLVP1,                     &
!    & PBBR,PALTE,                                                     &
!    & pduh2of,PDUCO2,PDUO3,palogp,palogt,PODSF,                       &
!    & PODAF,PBSFF,                                                    &
!    & pflfu,pflfd )
      END IF

!     Storage of 'grey' fluxes and their inverse (**igasm1)

      DO j3 = ki3s,ki3e+1
          DO j1 = ki1sc,ki1ec
          zflux (j1,j3) = pflfu(j1,j3)+pflcu(j1,j3)   &
                            -pflfd(j1,j3)-pflcd(j1,j3)
          zfluxi(j1,j3) = 1._wp/MIN( -ZEPFLX, zflux(j1,j3) )**IGASM1
          END DO
      END DO
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
               END IF    ! No.of relevant gases
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      igasz = 0
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      DO jg=3,1,-1   !     Loop over gases
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 
      IF (ncgas(jspec,jg).GT.1) THEN   ! include gas only, if necessary
      igasz = igasz + 1
!
      DO jjg = 1,3
      ICGAS(jjg) = 0   !   Set absorption coefficient index for all
                       !   gases to zero
      END DO   
 
      DO j3 = ki3s,ki3e+1
          DO j1 = ki1sc,ki1ec
          zfgas(j1,j3) = 0.0_wp   ! Preset 'gaseous' flux
          END DO
      END DO
!
      icc = NCGAS(jspec,jg)       ! No.of relevant coefficients    
!
! -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
      DO jc = icc,1,-1        ! Loop over absorption coefficients
! -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
      zcoai = coai(jc,jspec,jg)
      zcobi = cobi(jc,jspec,jg)
 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF ( ((zcoai.GE.zepai).AND.(zcobi.GT.0.0_wp))  &
         .OR. (igase.EQ.1) ) THEN    ! Solve linear system, if necessary
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      icgas(jg) = jc
      IF (icrf.EQ.0) THEN  ! partially cloudy atmosphere
      CALL inv_th (                                                     &
        pclc   ,pca1   ,pca2  ,pcb1  ,pcb2  ,pcc1  ,pcc2 ,pcd1 ,pcd2 ,  &
        pduh2oc,pduh2of,pduco2,pduo3 ,palogp,palogt,                    &
        podsc  ,podsf  ,podac ,podaf ,pbsfc ,pbsff ,                    &
        pbbr   ,palth,                                                  &
        jspec  ,icgas(1),icgas(2),icgas(3),                             &
        kig1s  ,kig1e  ,              ki3s ,ki3e ,                      &
        ki1sc  ,ki1ec  ,                                                &
        ldebug_inv_th ,                                                 &
        pflcu  ,pflfu  ,pflcd ,pflfd )
      ELSE                     ! 'cloud-free' atmosphere    
      PRINT *,' CRF not yet implemented'
!     CALL INVLONF(JSPEC,ICGAS(1),ICGAS(2),ICGAS(3),                   &
!    &             IRLON,IRLEV,IRLVP1,                                 &
!    & PBBR,PALTE,                                                     &
!    & pduh2of,PDUCO2,PDUO3,palogp,palogt,PODSF,                       &
!    & PODAF,PBSFF,                                                    &
!    & pflfu,pflfd )
      END IF
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ELSE               ! use 'grey' fluxes directly
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      DO j3 = ki3s,ki3e+1
          DO j1 = ki1sc,ki1ec
          pflfu(j1,j3) = zflux(j1,j3)
          pflfd(j1,j3) = 0.0_wp
          pflcu(j1,j3) = 0.0_wp
          pflcd(j1,j3) = 0.0_wp
          END DO
      END DO
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      END IF            ! Necessity to calculate fluxes
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 
      DO j3 = ki3s,ki3e+1
          DO j1 = ki1sc,ki1ec
          zfgas(j1,j3) = zfgas(j1,j3) + ZCOAI*          &
                           (pflfu(j1,j3)+pflcu(j1,j3)   &
                           -pflfd(j1,j3)-pflcd(j1,j3))
          END DO
      END DO
!
!?????????????????????????????????????????????????????????????????
      IF (ldebug_th) THEN
      PRINT *,' FESFT in debug mode for thermal fluxes'
      PRINT *,' only one interval/coefficient considered '
      PRINT *,'zfgas(j1b,ki3s): ',zfgas(j1b,ki3s)
      GO TO 22    
      END IF
!?????????????????????????????????????????????????????????????????
! -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
      END DO      ! Loop over absorption coefficients
! -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
 
!     Combination of inverse of 'grey' fluxes and actual gaseous flux
 
      DO j3 = ki3s,ki3e+1
          DO j1 = ki1sc,ki1ec
          zfluxi(j1,j3) = zfluxi(j1,j3)*zfgas(j1,j3)
          END DO
      END DO
 
      IF (igasz.EQ.igasm1) THEN    !     Avoid unphysical pseudo-transmission
      DO j3 = ki3s,ki3e+1
          DO j1 = ki1sc,ki1ec
          zfluxi(j1,j3) = MIN( 1._wp, MAX( 0._wp, zfluxi(j1,j3) ) )
          END DO
      END DO
      END IF
    
          END IF                   ! Test, whether gas needs to be included
    
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      END DO          ! End of loop over gases
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 

      DO j3 = ki3s,ki3e+1
          DO j1 = ki1sc,ki1ec
          zflux(j1,j3) = zfluxi(j1,j3)  ! Store FESFT result in zflux
          END DO
      END DO
!--------------------------------------------------------------
      END IF                             ! ESFT/FESFT-Selection 
!--------------------------------------------------------------

!     Addition of flux for spectral interval to total thermal flux
 
      IF (ICRF.EQ.0) THEN     ! Add flux at all levels
      DO j3 = ki3s,ki3e+1
          DO j1 = ki1sc,ki1ec
          PFLT(j1,j3) = PFLT(j1,j3) + zflux(j1,j3)
          END DO
      END DO
      ELSE                    ! Add flux only at TOA and surface
        DO j1 = ki1sc,ki1ec
        PFLTF(j1,1) = PFLTF(j1,1) + zflux(j1,ki3s  )
        PFLTF(j1,2) = PFLTF(j1,2) + zflux(j1,ki3e+1)
        END DO
      END IF
 
! ===============================================================
    2 CONTINUE                          ! End of spectral loop
   22 CONTINUE                          ! Exit branch for debug mode

      END IF                            ! switch for inclusion/exclusion
                                        ! of thermal calculations

! ===============================================================




! ###############################################################
!     Solar flux calculations
      IF (LSOLAR) THEN     ! if solar calculations are required
! ###############################################################

!     Inverse of cosine of zenith angle and surface albedo for
!     parallel radiation

        DO j1 = ki1sc,ki1ec
        pqsmu0(j1) = 1._wp/psmu0(j1)
        palp  (j1) = (1._wp+0.5_wp*(psmu0(j1)*(1._wp/palso(j1)-1._wp)))  &
                       /(1._wp +          (psmu0(j1)*(1._wp/palso(j1)-1._wp)))**2
        END DO

! ===============================================================
!     Loop over solar spectral intervals
      DO 3 jspec=1,jpsol
! ===============================================================

      DO j3 = ki3s,ki3e+1
          DO j1 = ki1sc,ki1ec
          zflux (j1,j3) = 0.0_wp  ! Preset flux in spectral interval
          zfluxd(j1,j3) = 0.0_wp
          zfluxu(j1,j3) = 0.0_wp
          END DO
      END DO     ! End of vertical loop
 
!     Upper boundary condition and reference pressure for Rayleigh sc.

        DO j1 = ki1sc,ki1ec
        pflpt(j1) = psct*solant(jspec)*psmu0(j1)
        papra(j1) = papre(j1)*pqsmu0(j1)
        END DO
 
!     Optical properties of non-gaseous constituents        
 
      IF (ldebug) THEN
      PRINT *,' FESFT    Call to opt_so for jspec: ',jspec
      END IF
      CALL       opt_so(prhocwc,pducwc,prhoiwc,pduiwc,                 &
                        paeq1  ,paeq2 ,paeq3  ,paeq4 , paeq5 ,         &
                        pdp    ,papra ,psmu0  ,pqsmu0,                 &
                        jspec  ,kig1s, kig1e ,ki1sc ,ki1ec  ,          &
                                ki3s ,ki3e  ,ldebug_opt_so,            &
                        podac  ,podaf ,podsc  ,podsf , pbsfc ,pbsff ,  &
                        pusfc  ,pusff  )

!     Selection of ESFT or FESFT method for interval considered
!--------------------------------------------------------------
      IF (nfast(jspec).EQ.0) THEN           ! ESFT method
!--------------------------------------------------------------

      IH2O = ncgas(jspec,1)
      ICO2 = ncgas(jspec,2)
      IO3  = ncgas(jspec,3)
 
      DO jh2o=1,IH2O    ! Loop over H2O coefficients
      DO jco2=1,ICO2    ! Loop over CO2 coefficients
      DO jo3 =1,IO3     ! Loop over O3  coefficients
 
      zaiprod=coai(jh2o,jspec,1)*coai(jco2,jspec,2)*coai(jo3,jspec,3)
!
      IF (icrf.EQ.0) THEN      ! partially cloudy atmosphere
      CALL       inv_so (                                              &
        pclc   ,pca1   ,pca2  ,pcb1  ,pcb2  ,pcc1  ,pcc2 ,pcd1 ,pcd2 , &
        pflpt  ,psmu0  ,pqsmu0,palp  ,palso ,                          &
        pduh2oc,pduh2of,pduco2,pduo3 ,palogp,palogt,                   &
        podsc  ,podsf  ,podac ,podaf ,pbsfc ,pbsff ,pusfc,pusff,       &
        jspec  ,jh2o   ,jco2  ,jo3   ,                                 &
        kig1s  ,kig1e  ,              ki3s ,ki3e ,                     &
        ki1sc  ,ki1ec  ,                                               &
        ldebug_inv_so ,                                                &
        pflcu  ,pflfu  ,pflcd ,pflfd ,pflcp ,pflfp)
 
      ELSE                     ! cloud-free calculation
      PRINT *,' CRF-Code not implemented yet'
!     CALL INVSHOF(JSPEC,jh2o,jco2,jo3,IRLON,JLA,JLE,IRLEV,IRLVP1,  &
!    &     PFLPT,pflfu,pflfd,PFLFP,                                 &
!    &     pduh2of,PDUCO2,PDUO3,palogp,palogt,                      &
!    &     PODSF  ,PODAF ,PBSFF ,PUSFF,                             &
!    &     PSMU0  ,PQSMU0 ,PALP  ,PALSO     )
 
      END IF
 
!     Incrementation of flux in spectral interval

      DO j3 = ki3s,ki3e+1
          DO j1 = ki1sc,ki1ec
          zflux (j1,j3)= zflux (j1,j3)                           &
                           +zaiprod*(pflfp(j1,j3)+pflcp(j1,j3) )
          zfluxd(j1,j3)= zfluxd(j1,j3)                           &
                           +zaiprod*(pflfd(j1,j3)+pflcd(j1,j3) )
          zfluxu(j1,j3)= zfluxu(j1,j3)                           &
                           +zaiprod*(pflfu(j1,j3)+pflcu(j1,j3) )
          END DO
      END DO     ! End of vertical loop
 
      END DO           ! Loop over O3 absorption coefficients
      END DO           ! Loop over CO2 absorption coefficients
      END DO           ! Loop over H2O absorption coefficients

!--------------------------------------------------------------
      ELSE                                  ! FESFT method
!--------------------------------------------------------------

      igase = 0
      DO jg=1,3
      ICGAS(jg) = 0
      IF (NCGAS(JSPEC,jg).GT.1) THEN
      igase = igase + 1
      END IF
      END DO
      igasm1    = igase -1
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
               IF (igase.le.1) THEN  !(no 'grey' fluxes required)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
               DO j3 = ki3s,ki3e+1
                   DO j1 = ki1sc,ki1ec
                   zflux  (j1,j3) = 1.0_wp
                   zfluxd (j1,j3) = 1.0_wp
                   zfluxu (j1,j3) = 1.0_wp
                   zfluxi (j1,j3) = 1.00_wp
                   zfluxdI(j1,j3) = 1.00_wp
                   zfluxuI(j1,j3) = 1.00_wp
                   END DO
               END DO     ! End of vertical loop
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
               ELSE   ! more than 1 gas --> 'grey' fluxes required
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF (icrf.EQ.0) THEN      ! partially cloudy atmosphere
      CALL       inv_so (                                              &
        pclc   ,pca1   ,pca2  ,pcb1  ,pcb2  ,pcc1  ,pcc2 ,pcd1 ,pcd2 , &
        pflpt  ,psmu0  ,pqsmu0,palp  ,palso ,                          &
        pduh2oc,pduh2of,pduco2,pduo3 ,palogp,palogt,                   &
        podsc  ,podsf  ,podac ,podaf ,pbsfc ,pbsff ,pusfc,pusff,       &
        jspec  ,0      ,0     ,0     ,                                 &
        kig1s  ,kig1e  ,              ki3s ,ki3e ,                     &
        ki1sc  ,ki1ec  ,                                               &
        ldebug_inv_so ,                                                &
        pflcu  ,pflfu  ,pflcd ,pflfd ,pflcp ,pflfp)
      ELSE                     ! cloud-free calculation
      PRINT *,' CRF-Code not implemented yet !!!'
      STOP 'no-crf'
!     CALL INVSHOF(JSPEC,0,0,0,IRLON,JLA,JLE,IRLEV,IRLVP1,  &
!    &     PFLPT,pflfu,pflfd,PFLFP,                         &
!    &     pduh2of,PDUCO2,PDUO3,palogp,palogt,              &
!    &     PODSF  ,PODAF ,PBSFF ,PUSFF,                     &
!    &     PSMU0  ,PQSMU0 ,PALP  ,PALSO     )
      END IF

!     Storage of 'grey' fluxes and their inverse (**igasm1)

      DO j3 = ki3s,ki3e+1
          DO j1 = ki1sc,ki1ec
          zfluxi (j1,j3) = 1._wp                           &
             /MAX(pflfp(j1,j3)+pflcp(j1,j3),ZEPFLX)**IGASM1
          zfluxdi(j1,j3) = 1._wp                           &
             /MAX(pflfd(j1,j3)+pflcd(j1,j3),ZEPFLX)**IGASM1
          zfluxui(j1,j3) = 1._wp                           &
            /MAX(pflfu(j1,j3)+pflcu(j1,j3),ZEPFLX)**IGASM1
          zflux  (j1,j3) = pflfp(j1,j3)+pflcp(j1,j3)
          zfluxd (j1,j3) = pflfd(j1,j3)+pflcd(j1,j3)
          zfluxu (j1,j3) = pflfu(j1,j3)+pflcu(j1,j3)
          END DO
      END DO     ! End of vertical loop
      
          IF (ldebug_so) THEN
          PRINT *,' FESFT in debug mode for solar   fluxes'
          PRINT *,' Grey fluxes   '
          DO j3=ki3s,ki3e+1
          PRINT *,'par/down/up    : ',zflux (j1b,j3)    &
                                     ,zfluxd(j1b,j3)    &
                                     ,zfluxu(j1b,j3),j3
          END DO
          END IF
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
               END IF    ! No.of relevant gases
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      igasz = 0
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      DO jg=3,1,-1   !     Loop over gases
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 
      IF (ncgas(jspec,jg).GT.1) THEN   ! include gas only, IF necessary
   
      igasz = igasz + 1
!
      DO jjg = 1,3
      ICGAS(jjg) = 0
      END DO 
 
      DO j3 = ki3s,ki3e+1
          DO j1 = ki1sc,ki1ec
          zfgas (j1,j3) = 0.0_wp    ! Initialize 'gaseous' fluxes
          zfgasD(j1,j3) = 0.0_wp
          zfgasU(j1,j3) = 0.0_wp
          END DO
      END DO     ! End of vertical loop
 
! -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
      icc = NCGAS(jspec,jg)
      DO jc = icc,1,-1        ! Loop over absorption coefficients
! -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
      zcoai = coai(jc,jspec,jg)
      zcobi = cobi(jc,jspec,jg)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF ( ((zcoai.ge.zepai).and.(zcobi.gt.0.0_wp))           &
         .OR. (igase.EQ.1) ) THEN    ! Solve linear system, if necessary
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ICGAS(jg) = jc
      IF (icrf.EQ.0) THEN      ! partially cloudy atmosphere
      CALL       inv_so (                                               &
        pclc   ,pca1   ,pca2  ,pcb1  ,pcb2  ,pcc1  ,pcc2 ,pcd1 ,pcd2 ,  &
        pflpt  ,psmu0  ,pqsmu0,palp  ,palso ,                           &
        pduh2oc,pduh2of,pduco2,pduo3 ,palogp,palogt,                    &
        podsc  ,podsf  ,podac ,podaf ,pbsfc ,pbsff ,pusfc,pusff,        &
        jspec  ,icgas(1),icgas(2),icgas(3),                             &
        kig1s  ,kig1e  ,              ki3s ,ki3e ,                      &
        ki1sc  ,ki1ec  ,                                                &
        ldebug_inv_so ,                                                 &
        pflcu  ,pflfu  ,pflcd ,pflfd ,pflcp ,pflfp)
      ELSE                     ! cloud-free clculations
      PRINT *,'crf-code not yet implemented'
      STOP 'no-crf'
!     CALL INVSHOF(JSPEC,ICGAS(1),ICGAS(2),ICGAS(3),IRLON,JLA,JLE,  &
!    &     IRLON,JLA,JLE,IRLEV,IRLVP1,                              &
!    &     PFLPT,pflfu,pflfd,PFLFP,                                 &
!    &     pduh2of,PDUCO2,PDUO3,palogp,palogt,                      &
!    &     PODSF  ,PODAF ,PBSFF ,PUSFF,                             &
!    &     PSMU0  ,PQSMU0 ,PALP  ,PALSO     )
      END IF
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ELSE               ! use 'grey' fluxes directly
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      DO j3 = ki3s,ki3e+1
          DO j1 = ki1sc,ki1ec
          pflfp(j1,j3) = zflux (j1,j3)
          pflcp(j1,j3) = 0.0_wp
          pflfd(j1,j3) = zfluxd(j1,j3)
          pflcd(j1,j3) = 0.0_wp
          pflfu(j1,j3) = zfluxu(j1,j3)
          pflcu(j1,j3) = 0.0_wp
          END DO
      END DO     ! End of vertical loop
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      END IF            ! Necessity to calculate fluxes
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 
      DO j3 = ki3s,ki3e+1
          DO j1 = ki1sc,ki1ec
          zfgas (j1,j3)= zfgas (j1,j3)                          &
                           +zcoai*( pflfp(j1,j3)+pflcp(j1,j3))
          zfgasd(j1,j3)= zfgasd(j1,j3)                          &
                           +zcoai*( pflfd(j1,j3)+pflcd(j1,j3))
          zfgasu(j1,j3)= zfgasu(j1,j3)                          &
                           +zcoai*( pflfu(j1,j3)+pflcu(j1,j3))
          END DO
      END DO     ! End of vertical loop
 
!?????????????????????????????????????????????????????????????????
      IF (ldebug_so) THEN
      PRINT *,' FESFT in debug mode for solar   fluxes'
      PRINT *,' only one interval/coefficient considered '
      PRINT *,' zcoai = ',zcoai
      DO j3=ki3s,ki3e+1
      PRINT *,'zfgas(j1b,j3): ',zfgas(j1b,j3)
      END DO
      GO TO 33    
      END IF
!?????????????????????????????????????????????????????????????????
! -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
! -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
      END DO      ! Loop over absorption coefficients
! -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
 
!     Combination of inverse of 'grey' fluxes and actual gaseous flux
 
      DO j3 = ki3s,ki3e+1
          DO j1 = ki1sc,ki1ec
          zfluxi (j1,j3) = zfluxi (j1,j3)*zfgas (j1,j3)
          zfluxdi(j1,j3) = zfluxdi(j1,j3)*zfgasd(j1,j3)
          zfluxui(j1,j3) = zfluxui(j1,j3)*zfgasu(j1,j3)
          END DO
      END DO     ! End of vertical loop
!
      IF (igasz.EQ.igasm1) THEN    !     Avoid unphysical pseudo-transmission
      DO j3 = ki3s,ki3e+1
          DO j1 = ki1sc,ki1ec
          zfluxi (j1,j3) = MIN( 1._wp, MAX( 0._wp, zfluxi (j1,j3) ) )
          zfluxdi(j1,j3) = MIN( 1._wp, MAX( 0._wp, zfluxdi(j1,j3)) )
          zfluxui(j1,j3) = MIN( 1._wp, MAX( 0._wp, zfluxui(j1,j3)) )
          END DO
      END DO     ! End of vertical loop
      END IF

          END IF                   ! Test, whether gas needs to be included

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      END DO          ! End of loop over gases
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
      DO j3 = ki3s,ki3e+1
          DO j1 = ki1sc,ki1ec
          zflux (j1,j3) = zfluxi (j1,j3)   ! Store flux, so that it
          zfluxd(j1,j3) = zfluxdi(j1,j3)   ! can be used in the same
          zfluxu(j1,j3) = zfluxui(j1,j3)   ! way as with ESFT-method
          END DO
      END DO     ! End of vertical loop
 
!--------------------------------------------------------------
      END IF                             ! ESFT/FESFT-Selection 
!--------------------------------------------------------------

!     Addition of flux for spectral interval to total solar flux
 
      IF (ICRF.EQ.0) THEN     ! Add flux at all levels

      DO j3 = ki3s,ki3e+1
          DO j1 = ki1sc,ki1ec
          pfls(j1,j3) = pfls  (j1,j3)                    &
                          +zflux (j1,j3)+zfluxd(j1,j3)   &
                          -zfluxu(j1,j3)
          END DO
      END DO     ! End of vertical loop

          DO j1 = ki1sc,ki1ec
            pflsp(j1) = pflsp (j1)  + zflux  (j1,ki3e+1)
            pflsd(j1) = pflsd (j1)  + zfluxd (j1,ki3e+1)
            pflsu(j1) = pflsu (j1)  + zfluxu (j1,ki3e+1)
          END DO
        IF (jspec.EQ.3) THEN   ! Photosynthetic active radiation
            DO j1 = ki1sc,ki1ec
            pflpar(j1) = pflpar(j1)                               &
                           +zflux (j1,ki3e+1)+ zfluxd(j1,ki3e+1)  &
                           -zfluxu(j1,ki3e+1)
            END DO
        END IF
      ELSE                    ! Add flux only at TOA and surface
        DO j1 = ki1sc,ki1ec
        pflsf(j1,1)  = pflsf (j1,1)                        &
                         +zflux (j1,ki3s)+zfluxd(j1,ki3s)  &
                         -zfluxu(j1,ki3s)
        pflsf(j1,2)  = pflsf (j1,2)                            &
                         +zflux (j1,ki3e+1)+zfluxd(j1,ki3e+1)  &
                         -zfluxu(j1,ki3e+1)
        END DO
      END IF
 
! ===============================================================
    3 CONTINUE                       ! End of solar spectral loop
   33 CONTINUE                       ! Exit branch for debug mode
! ===============================================================
 
! ###############################################################
      END IF     ! Test, whether solar calculation or not
! ###############################################################
 
!     Repeat calculations for cloud-free fluxes IF switch for CRF
!     is set to .true. and cloud-free fluxes have not yet been
!     computed

      IF (.not.lcrf) THEN
        DO j1=ki1sc,ki1ec
        pflsf(j1,1)  = pfls(j1,ki3s  )
        pflsf(j1,2)  = pfls(j1,ki3e+1)
        pfltf(j1,1)  = pflt(j1,ki3s  )
        pfltf(j1,2)  = pflt(j1,ki3e+1)
        END DO
      ELSE IF (icrf.EQ.0) THEN  ! Branch to cloud-free calculations
      ICRF = 1                  ! only once 
      GO TO 1
      END IF

!!$        WRITE(message_text,'(a,2E15.5)') ' max/min pflt(:,:)   = ',&
!!$             & MAXVAL(pflt(:,:) ), MINVAL( pflt(:,:) )
!!$        CALL message('', TRIM(message_text))
!!$        WRITE(message_text,'(a,2E15.5)') ' max/min pfls(:,:)   = ',&
!!$             & MAXVAL(pfls(:,:) ), MINVAL( pfls(:,:) )
!!$        CALL message('', TRIM(message_text))
      
      RETURN
      END SUBROUTINE fesft

      SUBROUTINE opt_th(prholwc,pdulwc,prhoiwc,pduiwc,                 &
     &                  paeq1  ,paeq2 ,paeq3  ,paeq4 , paeq5 ,         &
     &                  kspec  ,kig1s, kig1e, ki1sc ,ki1ec  ,          &
     &                          ki3sc ,ki3ec  ,ldebug,                 &
     &                  podac  ,podaf ,podsc  ,podsf , pbsfc ,pbsff  )
!
 
!**** *opt_th* - calculates the optical properties of the non-gaseous
!                constituents for one spectral
!                interval in the thermal part of the spectrum
 
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
!  Include gme_comconst.h for the gravity constant
! V2_24        2010/04/28 Michael Gertz
!  Adaptions to SVN       
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!=======================================================================
!     Purpose:   Calculation of (absorption and scattering) optical
!                depth and backward scattered fraction of diffuse
!                radiation (excluding the contribution by gaseous
!                constituents)

!     METHOD
!     ------
 
!     - determination of optical properies (i.e. extinction coefficient,
!       single scattering albedo and asymetry factor of the phase function)
!       using approximate relations for cloud water and cloud ice and
!       combination of five type of aerosols
!
!     - calculation of optical depth (scattering and absorption) and back-
!       scattered fraction suitable for use in an implicit delta-two-stream
!       scheme
!       
!==============================================================================

      !     Input data
      INTEGER kig1s,kig1e
      INTEGER ki1sc          ! start index for first  array dimension
      INTEGER ki1ec          ! end   index for first  array dimension
      INTEGER ki3sc          ! start index for third  array dimension
      INTEGER ki3ec          ! end   index for third  array dimension

!     Liquid and ice water density and content within for the cloudy
!     part of each layer

      REAL(wp) :: prholwc(kig1s:kig1e,ki3sc:ki3ec) 
      REAL(wp) :: prhoiwc(kig1s:kig1e,ki3sc:ki3ec) 
      REAL(wp) :: pdulwc (kig1s:kig1e,ki3sc:ki3ec) 
      REAL(wp) :: pduiwc (kig1s:kig1e,ki3sc:ki3ec) 

!     Aerosole contents (optical depths at 0.55 micrometer) for 5 types
!
      REAL(wp) :: paeq1  (kig1s:kig1e,ki3sc:ki3ec) 
      REAL(wp) :: paeq2  (kig1s:kig1e,ki3sc:ki3ec) 
      REAL(wp) :: paeq3  (kig1s:kig1e,ki3sc:ki3ec) 
      REAL(wp) :: paeq4  (kig1s:kig1e,ki3sc:ki3ec) 
      REAL(wp) :: paeq5  (kig1s:kig1e,ki3sc:ki3ec) 
 
!
      INTEGER kspec          ! selected spectral interval
      LOGICAL ldebug         ! debug control switch       

!     Output data

      REAL(wp) :: podac (kig1s:kig1e,ki3sc:ki3ec) ! absorption optical depth
      REAL(wp) :: podaf (kig1s:kig1e,ki3sc:ki3ec) ! in cloudy and free part
      REAL(wp) :: podsc (kig1s:kig1e,ki3sc:ki3ec) ! scattering optical depth
      REAL(wp) :: podsf (kig1s:kig1e,ki3sc:ki3ec) ! in cloudy and free part
      REAL(wp) :: pbsfc (kig1s:kig1e,ki3sc:ki3ec) ! backscattering fraction 
      REAL(wp) :: pbsff (kig1s:kig1e,ki3sc:ki3ec) ! in cloudy and free part


!     Local arrays and variables

!     optical properties (absorption, scattering, backscatter fraction)
!     for liquid water, ice and total aerosole

      REAL(wp) :: zlwoda (kig1s:kig1e) !
      REAL(wp) :: zlwods (kig1s:kig1e) ! 
      REAL(wp) :: zlwb0  (kig1s:kig1e) ! 
      REAL(wp) :: ziwoda (kig1s:kig1e) ! 
      REAL(wp) :: ziwods (kig1s:kig1e) ! 
      REAL(wp) :: ziwb0  (kig1s:kig1e) ! 
      REAL(wp) :: zaeoda (kig1s:kig1e) ! 
      REAL(wp) :: zaeods (kig1s:kig1e) ! 
      REAL(wp) :: zaeb0  (kig1s:kig1e) ! 
      
!     individual optical properties of liquid water and ice 

      REAL(wp) :: z_lwe, z_iwe                 ! extinction coefficient
      REAL(wp) :: z_lww, z_iww                 ! single scattering coefficient
      REAL(wp) :: z_lwg, z_iwg                 ! asymetry factor
      REAL(wp) :: z_lwf, z_iwf                 ! forward scattered fraction 

      INTEGER j1,j3                 ! loop indices over spatial dimensions
      INTEGER ja                    ! local loop index
      INTEGER j1b                   ! indices of debug point


!     Spectral properties 

!      REAL(wp) :: zaea(jpspec,5),zaes(jpspec,5),                 &
!     &                     zaeg(jpspec,5),zaef(jpspec,5)
!      REAL(wp) :: zlwe(4,jpspec),zlwemn(jpspec),zlwemx(jpspec),  &
!     &                     zlww(2,jpspec),zlwg(2,jpspec)
!      REAL(wp) :: ziwe(4,jpspec),ziweMN(jpspec),ziwemx(jpspec),  &
!     &                     ziww(2,jpspec),ziwg(2,jpspec)

!    
!     Local utility variables
!
      REAL(wp) :: z1dg       ! 1./g
      REAL(wp) :: z1d8       ! 1./8
      REAL(wp) :: zzg

      REAL(wp) :: ZEPOPD     ! Security constant for optical depth
      REAL(wp) :: ZEPSSA     ! Security constant for single scattering albedo
!      INCLUDE "gme_dat_lw.h"
!      INCLUDE "gme_dat_iw.h"
!      INCLUDE "gme_dat_aer.h"
!      INCLUDE "gme_comconst.h"

      DATA ZEPOPD /1.E-06_wp/
      DATA ZEPSSA /1.E-06_wp/
      DATA j1b/1/
 
!     Z1DG     = 1.0_wp/9.80665_wp
      Z1DG     = 1.0_wp/G
      Z1D8     = 0.125_wp    
      IF (ldebug) THEN
      PRINT *,' **** opt-th   start ********************'
      PRINT *,' **** debug point : ',j1b
      END IF
 
      DO ja=1,5        
      zaef(kspec,ja)  = zaeg(kspec,ja)**2 ! forward sc.fraction f.aerosols
      END DO   
      IF (ldebug) THEN
      DO ja=1,5
      PRINT *,'ja, zaef(kspec,ja): ',ja,zaef(kspec,ja)
      END DO
      END IF
      IF (ldebug) PRINT *,' In opt-th   vor vertical loop'
!     Vertical loop
 
      DO j3=ki3sc,ki3ec
 
      IF (ldebug) PRINT *,' In opt-th   j3 = ',j3  

!     Optical properties of liquid water and ice as function of the spec
!     liquid water and ice content
 
          DO j1=ki1sc,ki1ec

!     Liquid water

      z_lwg    = zlwg(1,kspec) + zlwg(2,kspec)*prholwc(j1,j3)
      z_lwg    = MAX (0.0_wp,MIN(1.0_wp,z_lwg))
      z_lwf    = z_lwg*z_lwg
      z_lww    = zlww(1,kspec) + zlww(2,kspec)*prholwc(j1,j3)
      z_lww    = MAX(zepssa,MIN(1.0_wp,z_lww)) 
      z_lwe    = z1dg * (zlwe(1,kspec) + zlwe(2,kspec)/           &
     &           (zlwe(3,kspec)*prholwc(j1,j3)+zlwe(4,kspec)))
      z_lwe    = MAX(zlwemn(kspec),MIN(zlwemx(kspec),z_lwe))

      zlwoda(j1)= z_lwe *pdulwc(j1,j3) * (1._wp-z_lww) 
      zlwods(j1)= z_lwe *pdulwc(j1,j3) *     z_lww  * (1._wp-z_lwf)
      zlwb0 (j1)= z1d8*(4._wp+z_lwg)/(1._wp+z_lwg)
 
!     Ice
 
!     z_iwg    = ziwg(1,kspec) + ziwg(2,kspec)*prhoiwc(j1,j3)
      z_iwg    = ziwg(1,kspec) + ziwg(2,kspec)*LOG(prhoiwc(j1,j3))
      z_iwg    = MAX(0.0_wp,MIN(1.0_wp,z_iwg))
      z_iwf    = z_iwg*z_iwg
!     z_iww    = ziww(1,kspec) + ziww(2,kspec)*prhoiwc(j1,j3)
      z_iww    = ziww(1,kspec) + ziww(2,kspec)*LOG(prhoiwc(j1,j3))
      z_iww    = MAX(1.E-12_wp,MIN(1.0_wp,z_iww))
      z_iwe    = z1dg * (ziwe(1,kspec) + ziwe(2,kspec)/            &
     &           (ziwe(3,kspec)*prhoiwc(j1,j3)+ziwe(4,kspec)))
      z_iwe    = MAX(ziwemn(kspec),MIN(ziwemx(kspec),z_iwe  ))
 
      ziwoda(j1) = z_iwe * pduiwc(j1,j3)*(1.0_wp-z_iww)
      ziwods(j1) = z_iwe * pduiwc(j1,j3)*     z_iww  *(1._wp-z_iwf)
      ziwb0 (j1) = z1d8*(4._wp+z_iwg)/(1._wp+z_iwg)
          END DO

      IF (ldebug) THEN
      PRINT *,' prholwc (j1b) :',prholwc (j1b,j3)
      PRINT *,' pdulwc  (j1b) :',pdulwc  (j1b,j3)
      PRINT *,' zlwoda  (j1b) :',zlwoda  (j1b)
      PRINT *,' zlwods  (j1b) :',zlwods  (j1b)
      PRINT *,' zlwb0   (j1b) :',zlwb0   (j1b)
      PRINT *,' z_lwg                 :',z_lwg                
      PRINT *,' z_lwf                 :',z_lwf                
      PRINT *,' z_lww                 :',z_lww                
      PRINT *,' z_lwe                 :',z_lwe                
      PRINT *,' prhoiwc (j1b) :',prhoiwc (j1b,j3)
      PRINT *,' pduiwc  (j1b) :',pduiwc  (j1b,j3)
      PRINT *,' ziwoda  (j1b) :',ziwoda  (j1b)
      PRINT *,' ziwods  (j1b) :',ziwods  (j1b)
      PRINT *,' ziwb0   (j1b) :',ziwb0   (j1b)
      END IF
 
!     Aerosoles
 
          DO j1=ki1sc,ki1ec
      zaeoda(j1) = paeq1(j1,j3)*zaea(kspec,1)  &
     &               +paeq2(j1,j3)*zaea(kspec,2)  &
     &               +paeq3(j1,j3)*zaea(kspec,3)  &
     &               +paeq4(j1,j3)*zaea(kspec,4)  &
     &               +paeq5(j1,j3)*zaea(kspec,5)
      zaeods(j1) =                                         &
     &  ( paeq1(j1,j3)*zaes(kspec,1)*(1._wp-zaef(kspec,1)) )  &
     & +( paeq2(j1,j3)*zaes(kspec,2)*(1._wp-zaef(kspec,2)) )  &
     & +( paeq3(j1,j3)*zaes(kspec,3)*(1._wp-zaef(kspec,3)) )  &
     & +( paeq4(j1,j3)*zaes(kspec,4)*(1._wp-zaef(kspec,4)) )  &
     & +( paeq5(j1,j3)*zaes(kspec,5)*(1._wp-zaef(kspec,5)) )
!
      zzg=                                                                &
     &((paeq1(j1,j3)*zaes(kspec,1)*(1._wp-zaef(kspec,1)))*zaeg(kspec,1)   &
     &+(paeq2(j1,j3)*zaes(kspec,2)*(1._wp-zaef(kspec,2)))*zaeg(kspec,2)   &
     &+(paeq3(j1,j3)*zaes(kspec,3)*(1._wp-zaef(kspec,3)))*zaeg(kspec,3)   &
     &+(paeq4(j1,j3)*zaes(kspec,4)*(1._wp-zaef(kspec,4)))*zaeg(kspec,4)   &
     &+(paeq5(j1,j3)*zaes(kspec,5)*(1._wp-zaef(kspec,5)))*zaeg(kspec,5))  &
     &   / MAX( zaeods(j1),zepopd)
      zaeb0 (j1) = z1d8*(4._wp+zzg)/(1._wp+zzg)
          END DO

      IF (ldebug) THEN
      PRINT *,' zaeoda  (j1b) :',zaeoda  (j1b)
      PRINT *,' zaeods  (j1b) :',zaeods  (j1b)
      PRINT *,' zaeb0   (j1b) :',zaeb0   (j1b)
      END IF
 
!     Combined effects
 
!     a) cloud free part of each layer

          DO j1=ki1sc,ki1ec
          podaf(j1,j3)=MAX(zaeoda(j1),zepopd)   
          podsf(j1,j3)=MAX(zaeods(j1),zepopd)
          pbsff(j1,j3)=    zaeb0 (j1)
          END DO
      IF (ldebug) THEN
      PRINT *,' podaf   (j1b,j3) :',podaf   (j1b,j3)
      PRINT *,' podsf   (j1b,j3) :',podsf   (j1b,j3)
      PRINT *,' pbsff   (j1b,j3) :',pbsff   (j1b,j3)
      END IF
 
!     b) cloudy part of each layer

          DO j1=ki1sc,ki1ec
      podac(j1,j3)=MAX(zlwoda (j1)                          &
     &                  + ziwoda (j1) + zaeoda(j1) ,zepopd)
      podsc(j1,j3)=MAX(zlwods (j1)                          &
     &                  + ziwods (j1) + zaeods(j1) ,zepopd)
      podsc(j1,j3)=MIN(podsc(j1,j3),                        &
     &             (1._wp-zepssa)*(podac(j1,j3)+podsc(j1,j3)))
 
      pbsfc(j1,j3)= (zlwb0 (j1)*zlwods (j1)   &
     &                + ziwb0 (j1)*ziwods (j1)   &
     &                + zaeb0 (j1)*zaeods (j1)) / podsc(j1,j3)
          END DO
      IF (ldebug) THEN
      PRINT *,' podac   (j1b,j3) :',podac   (j1b,j3)
      PRINT *,' podsc   (j1b,j3) :',podsc   (j1b,j3)
      PRINT *,' pbsfc   (j1b,j3) :',pbsfc   (j1b,j3)
      END IF
 
      END DO       ! End of vertical loop
 
      RETURN
    END SUBROUTINE opt_th

    
    SUBROUTINE inv_th (                                               &
     &  pclc   ,pca1   ,pca2  ,pcb1  ,pcb2  ,pcc1  ,pcc2 ,pcd1 ,pcd2 ,  &
     &  pduh2oc,pduh2of,pduco2,pduo3 ,palogp,palogt,                    &
     &  podsc  ,podsf  ,podac ,podaf ,pbsfc ,pbsff ,                    &
     &  pbbr   ,palth,                                                  &
     &  kspec  ,kh2o   ,kco2  ,ko3   ,                                  &
     &  kig1s  ,kig1e  ,              ki3s ,ki3e ,                      &
     &  ki1sc  ,ki1ec  ,                                                &
     &  ldebug ,                                                        &
     &  pflcu  ,pflfu  ,pflcd ,pflfd )
!
!     *inv_th* solves the linear system of equations for thermal fluxes
!
!     Purpose
!     -------
!
!     The routine solves a linear equation system for thermal fluxes using
!     a Gaussian elimination-backsubstitution algorithm dedicated to the
!     specific structure of the system matrix
!
!     Method
!     ------
!
!     - setting of the RHS of the system using the layer boundary black
!       body radiation and allowing for partial cloud cover in each layer
!     - solution of the equation system including the lower boundary
!       condition
!     - matrix coefficients are calculated in the course of the elimination
!       step for one layer at a time through a call to routine *coe_th*
!     - the final result, i.e. the so-called black body flux differences
!       (cf.Ritter and Geleyn, 1992) are stored separately for cloudy and
!       cloud-free part of each layer boundary
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
!
!=======================================================================

!     Input data 
!
      INTEGER kig1s,kig1e,ki3s,ki3e   ! global array bounds
      INTEGER ki1sc,ki1ec             ! computational domain


      REAL(wp) :: pclc  (kig1s:kig1e,ki3s:ki3e) ! cloud cover
      REAL(wp) :: pca1  (kig1s:kig1e,ki3s:ki3e) ! cloud geometry factor  
      REAL(wp) :: pca2  (kig1s:kig1e,ki3s:ki3e) ! cloud geometry factor  
      REAL(wp) :: pcb1  (kig1s:kig1e,ki3s:ki3e) ! cloud geometry factor  
      REAL(wp) :: pcb2  (kig1s:kig1e,ki3s:ki3e) ! cloud geometry factor  
      REAL(wp) :: pcc1  (kig1s:kig1e,ki3s:ki3e) ! cloud geometry factor  
      REAL(wp) :: pcc2  (kig1s:kig1e,ki3s:ki3e) ! cloud geometry factor  
      REAL(wp) :: pcd1  (kig1s:kig1e,ki3s:ki3e) ! cloud geometry factor  
      REAL(wp) :: pcd2  (kig1s:kig1e,ki3s:ki3e) ! cloud geometry factor  

      REAL(wp) :: pbbr  (kig1s:kig1e,ki3s:ki3e+1) ! black body radiation   
      REAL(wp) :: palth (kig1s:kig1e)             ! surface albedo
 
!     Input arrays to be passed to *coe_th*
 
!     layer gas contents (cloudy and cloud-free, if distinction necessary)
      REAL(wp) :: pduh2oc(kig1s:kig1e,ki3s:ki3e) ! h2o-vapour cloudy      
      REAL(wp) :: pduh2of(kig1s:kig1e,ki3s:ki3e) ! h2o-vapour cloud-free  
      REAL(wp) :: pduco2 (kig1s:kig1e,ki3s:ki3e) ! co2
      REAL(wp) :: pduo3  (kig1s:kig1e,ki3s:ki3e) ! o3 
!     optical properties of 'grey' constituents (cloudy and cloud-free)
      REAL(wp) :: podsc  (kig1s:kig1e,ki3s:ki3e) ! scattering optical depth
      REAL(wp) :: podsf  (kig1s:kig1e,ki3s:ki3e) ! scattering optical depth
      REAL(wp) :: podac  (kig1s:kig1e,ki3s:ki3e) ! absorption optical depth
      REAL(wp) :: podaf  (kig1s:kig1e,ki3s:ki3e) ! absorption optical depth
      REAL(wp) :: pbsfc  (kig1s:kig1e,ki3s:ki3e) ! backscatter fraction
      REAL(wp) :: pbsff  (kig1s:kig1e,ki3s:ki3e) ! backscatter fraction
!     
      REAL(wp) :: palogp (kig1s:kig1e,ki3s:ki3e) ! ln(p)
      REAL(wp) :: palogt (kig1s:kig1e,ki3s:ki3e) ! ln(T)
 
!
      INTEGER kspec          ! selected spectral interval
      INTEGER kh2o           ! table index for h2o-vapour absorption
      INTEGER kco2           ! table index for co2        absorption
      INTEGER ko3            ! table index for o3         absorption
      LOGICAL ldebug         ! debug control switch       

!     Output data
 
      REAL(wp) :: pflcu (kig1s:kig1e,ki3s:ki3e+1) ! flux up   cloudy
      REAL(wp) :: pflfu (kig1s:kig1e,ki3s:ki3e+1) ! flux up   cloud-free 
      REAL(wp) :: pflcd (kig1s:kig1e,ki3s:ki3e+1) ! flux down cloudy     
      REAL(wp) :: pflfd (kig1s:kig1e,ki3s:ki3e+1) ! flux down cloud-free 


!     Local arrays and variables

!     layer properties calculated in *coe_th*

      REAL(wp) :: pa1c   (kig1s:kig1e) !
      REAL(wp) :: pa1f   (kig1s:kig1e) ! 
      REAL(wp) :: pa2c   (kig1s:kig1e) ! 
      REAL(wp) :: pa2f   (kig1s:kig1e) ! 
      REAL(wp) :: pa3c   (kig1s:kig1e) ! 
      REAL(wp) :: pa3f   (kig1s:kig1e) ! 

!     Utitlity arrays
      REAL(wp) :: ztu1 (kig1s:kig1e,ki3s:ki3e) !
      REAL(wp) :: ztu2 (kig1s:kig1e,ki3s:ki3e) !
      REAL(wp) :: ztu3 (kig1s:kig1e,ki3s:ki3e) !
      REAL(wp) :: ztu4 (kig1s:kig1e,ki3s:ki3e) !
      REAL(wp) :: ztu5 (kig1s:kig1e,ki3s:ki3e) !
      REAL(wp) :: ztu6 (kig1s:kig1e,ki3s:ki3e) !
      REAL(wp) :: ztu7 (kig1s:kig1e,ki3s:ki3e) !
      REAL(wp) :: ztu8 (kig1s:kig1e,ki3s:ki3e) !
      REAL(wp) :: ztu9 (kig1s:kig1e,ki3s:ki3e) !

!     Utility variables
      REAL(wp) :: ztd1 ,ztd2 ,ztd3 ,ztd4 ,ztd5 ,ztd6 , ztd7
      REAL(wp) :: ztds1,ztds2,ztds3,ztus1
      REAL(wp) :: zt1,zt2,zt3,zs1,zs2,zs3 ! Temporaries for optimization
      
      INTEGER j1,j3       ! loop indices over spatial dimensions
      INTEGER j1b        ! indices for debug point
      LOGICAL ldebug_coe_th  ! Debug switch for *coe_th*

      DATA j1b /1/

!     ------------------------------------------------------------------
!
      ldebug_coe_th = .FALSE. 

!     Upper boundary condition
 
        DO j1=ki1sc,ki1ec
        pflfd(j1,ki3s)=pbbr(j1,ki3s)
        pflcd(j1,ki3s)=0._wp
        END DO
#ifdef DEBUG
      IF (ldebug) THEN
        PRINT *,' *** INV_TH **************************'
        PRINT *,' *** debug point : ',j1b
        PRINT *,'pflfd(j1b,ki3s) : ',pflfd(j1b,ki3s)
        PRINT *,'pflcd(j1b,ki3s) : ',pflcd(j1b,ki3s)
      END IF
#endif
 
!     Determine effects of first layer in *coe_th*
      CALL       coe_th (                                &
     & pduh2oc,pduh2of,pduco2 ,pduo3  ,palogp ,palogt ,  &
     & podsc  ,podsf  ,podac  ,podaf  ,pbsfc  ,pbsff  ,  &
     & ki3s  ,kspec  ,kh2o   ,kco2   ,ko3    ,           &
     & kig1s ,kig1e, ki1sc  ,ki1ec  ,  ki3s  ,ki3e,      &
     & ldebug_coe_th ,                                   &
     & pa1c   ,pa1f   ,pa2c   ,pa2f   ,pa3c   ,pa3f)
 
!     Set RHS
 
        DO j1=ki1sc,ki1ec
        pflfu(j1,ki3s)=(1._wp-pclc(j1,ki3s))*pa3f(j1)*       &
     &                              (pbbr(j1,ki3s)-pbbr(j1,ki3s+1))
        pflcu(j1,ki3s)=    pclc(j1,ki3s) *pa3c(j1)*    &
     &                     (pbbr(j1,ki3s)-pbbr(j1,ki3s+1))
        pflfd(j1,ki3s+1)=-pflfu(j1,ki3s)
        pflcd(j1,ki3s+1)=-pflcu(j1,ki3s)
        END DO
 
#ifdef DEBUG
      IF (ldebug) THEN
        PRINT *,'pflfd(j1b,ki3s+1) : ',pflfd(j1b,ki3s+1)
      END IF
#endif
!
!     Elimination for first layer
        DO j1=ki1sc,ki1ec
        pflfu(j1,ki3s)  =pflfu(j1,ki3s  ) + pa2f(j1)*    &
     &                      (pca2 (j1,ki3s)*pflfd(j1,ki3s))
        pflfd(j1,ki3s+1)=pflfd(j1,ki3s+1) + pa1f(j1)*    &
     &                      (pca2 (j1,ki3s)*pflfd(j1,ki3s))
        pflcu(j1,ki3s)  =pflcu(j1,ki3s  ) + pa2c(j1)*    &
     &                      (pcb2 (j1,ki3s)*pflfd(j1,ki3s))
        pflcd(j1,ki3s+1)=pflcd(j1,ki3s+1) + pa1c(j1)*    &
     &                      (pcb2 (j1,ki3s)*pflfd(j1,ki3s))
        END DO
#ifdef DEBUG
      IF (ldebug) THEN
        PRINT *,' after elimination'
        PRINT *,'pflfd(j1b,ki3s+1) : ',pflfd(j1b,ki3s+1)
      END IF
#endif
 
!     Store some utitlity variables for first layer
        DO j1=ki1sc,ki1ec
        ztu1(j1,ki3s)=0._wp
        ztu2(j1,ki3s)=pca1(j1,ki3s)*pa1f(j1)
        ztu3(j1,ki3s)=pcc1(j1,ki3s)*pa1f(j1)
        ztu4(j1,ki3s)=pcb1(j1,ki3s)*pa1c(j1)
        ztu5(j1,ki3s)=pcd1(j1,ki3s)*pa1c(j1)
        ztu6(j1,ki3s)=pca1(j1,ki3s)*pa2f(j1)
        ztu7(j1,ki3s)=pcc1(j1,ki3s)*pa2f(j1)
        ztu8(j1,ki3s)=pcb1(j1,ki3s)*pa2c(j1)
        ztu9(j1,ki3s)=pcd1(j1,ki3s)*pa2c(j1)
        END DO
 
!     Vertical loop
 
      DO j3=ki3s+1,ki3e
!     Determine effect of the layer in *coe_th*
      ldebug_coe_th = .false.
      CALL       coe_th (                                &
     & pduh2oc,pduh2of,pduco2 ,pduo3  ,palogp ,palogt ,  &
     & podsc  ,podsf  ,podac  ,podaf  ,pbsfc  ,pbsff  ,  &
     & j3     ,kspec  ,kh2o   ,kco2   ,ko3    ,          &
     & kig1s ,kig1e, ki1sc  ,ki1ec  ,  ki3s  ,ki3e,      &
     & ldebug_coe_th ,                                   &
     & pa1c   ,pa1f   ,pa2c   ,pa2f   ,pa3c   ,pa3f)
 
!     Set RHS
 
          DO j1=ki1sc,ki1ec
          pflfu(j1,j3  )=(1._wp-pclc(j1,j3))*pa3f(j1)*     &
     &                                (pbbr(j1,j3)-pbbr(j1,j3+1))
          pflcu(j1,j3  )=    pclc(j1,j3) *pa3c(j1)*  &
     &                      (pbbr(j1,j3)-pbbr(j1,j3+1))
          pflfd(j1,j3+1)=-pflfu(j1,j3)
          pflcd(j1,j3+1)=-pflcu(j1,j3)
          END DO
 
#ifdef DEBUG
      IF (ldebug) THEN
        PRINT *,' in vertical loop j3=',j3
        PRINT *,'pflfd(j1b,j3+1) : ',pflfd(j1b,j3+1)
      END IF
#endif
!
!     Elimination and storage of utility variables
 
          DO j1=ki1sc,ki1ec

          zt1=( pca2(j1,j3)*ztu6 (j1,j3-1)  &
     &         +pcc2(j1,j3)*ztu8 (j1,j3-1))
          zt2=( pca2(j1,j3)*ztu7 (j1,j3-1)  &
     &         +pcc2(j1,j3)*ztu9 (j1,j3-1))
          zt3=( pca2(j1,j3)*pflfd(j1,j3)    &
     &         +pcc2(j1,j3)*pflcd(j1,j3))
          zs1=( pcb2(j1,j3)*ztu6 (j1,j3-1)  &
     &         +pcd2(j1,j3)*ztu8 (j1,j3-1))
          zs2=( pcb2(j1,j3)*ztu7 (j1,j3-1)  &
     &         +pcd2(j1,j3)*ztu9 (j1,j3-1))
          zs3=( pcb2(j1,j3)*pflfd(j1,j3)    &
     &         +pcd2(j1,j3)*pflcd(j1,j3))

          ztd1 = 1._wp/(1._wp-pa2f(j1)*zt1)
          pflfu(j1,j3) = ztd1*( pflfu(j1,j3) +  &
     &                  pa2f(j1)*zt3)
          ztu1 (j1,j3) = ztd1*                     &
     &                  pa2f(j1)*zt2
          ztu2 (j1,j3) =(ztd1*pa1f(j1))*pca1(j1,j3)
          ztu3 (j1,j3) =(ztd1*pa1f(j1))*pcc1(j1,j3)
          ztd2 = pa2c(j1)*zs1
          ztd3 = 1._wp/(1._wp-pa2c(j1)*zs2                     &
     &                 -ztd2*ztu1(j1,j3))
          pflcu(j1,j3) = ztd3*( pflcu(j1,j3) +        &
     &                  pa2c(j1)*zs3                     &
     &                 +ztd2*pflfu(j1,j3))
          ztu4 (j1,j3) = ztd3*( pa1c(j1)*pcb1(j1,j3)         &
     &                            +ztd2*ztu2(j1,j3))
          ztu5 (j1,j3) = ztd3*( pa1c(j1)*pcd1(j1,j3)         &
     &                            +ztd2*ztu3(j1,j3))
          ztd4 = pa1f(j1)*zt1
          ztd5 = pa1f(j1)*zt2
          pflfd(j1,j3+1) = pflfd(j1,j3+1)+pa1f(j1)*zt3       &
     &           + ztd4*pflfu(j1,j3) + ztd5*pflcu(j1,j3)
          ztu6 (j1,j3) = pa2f(j1)*pca1(j1,j3)                &
     &                     +ztd4*ztu2(j1,j3)+ztd5*ztu4(j1,j3)
          ztu7 (j1,j3) = pa2f(j1)*pcc1(j1,j3)                &
     &                     +ztd4*ztu3(j1,j3)+ztd5*ztu5(j1,j3)
          ztd6 = pa1c(j1)*zs1
          ztd7 = pa1c(j1)*zs2
          pflcd(j1,j3+1) = pflcd(j1,j3+1)+pa1c(j1)*zs3       &
     &           + ztd6*pflfu(j1,j3) + ztd7*pflcu(j1,j3)
          ztu8(j1,j3) = pa2c(j1)*pcb1(j1,j3)                 &
     &                    +ztd6*ztu2(j1,j3)+ztd7*ztu4(j1,j3)
          ztu9(j1,j3) = pa2c(j1)*pcd1(j1,j3)                 &
     &                    +ztd6*ztu3(j1,j3)+ztd7*ztu5(j1,j3)
          END DO
!
#ifdef DEBUG
      IF (ldebug) THEN
        PRINT *,' after elimination in vertical loop j3=',j3
        PRINT *,'pflfd(j1b,j3+1) : ',pflfd(j1b,j3+1)
      END IF
#endif
 
      END DO     ! End of vertical loop over layers
 
!     Elimination and backsubstitution at the surface
        DO j1=ki1sc,ki1ec
        ztds1    =1._wp/(1._wp-palth(j1)*ztu6(j1,ki3e))
        pflfu(j1,ki3e+1)= ztds1 *palth(j1)*pflfd(j1,ki3e+1)
        ztus1    =ztds1 *palth(j1)*ztu7(j1,ki3e)
        ztds2    =palth(j1)*ztu8(j1,ki3e)
        ztds3    =1._wp/(1._wp-palth(j1)*ztu9(j1,ki3e)-ztds2*ztus1)
        pflcu(j1,ki3e+1)=ztds3*( palth(j1)*pflcd(j1,ki3e+1)   &
     &                              +ztds2       *pflfu(j1,ki3e+1))
        pflfu(j1,ki3e+1)=pflfu(j1,ki3e+1)                        &
     &                            +ztus1*pflcu(j1,ki3e+1)
        END DO

!     Layer-by-layer backsubstitution
 
      DO j3 =ki3e,ki3s,-1
          DO j1=ki1sc,ki1ec
          pflcd(j1,j3+1) = pflcd(j1,j3+1)                  &
     &                       +ztu8 (j1,j3)*pflfu(j1,j3+1)  &
     &                       +ztu9 (j1,j3)*pflcu(j1,j3+1)
          pflfd(j1,j3+1) = pflfd(j1,j3+1)                  &
     &                       +ztu6 (j1,j3)*pflfu(j1,j3+1)  &
     &                       +ztu7 (j1,j3)*pflcu(j1,j3+1)
          pflcu(j1,j3  ) = pflcu(j1,j3)                    &
     &                       +ztu4 (j1,j3)*pflfu(j1,j3+1)  &
     &                       +ztu5 (j1,j3)*pflcu(j1,j3+1)
          pflfu(j1,j3  ) = pflfu(j1,j3)                    &
     &                       +ztu2 (j1,j3)*pflfu(j1,j3+1)  &
     &                       +ztu3 (j1,j3)*pflcu(j1,j3+1)  &
     &                       +ztu1 (j1,j3)*pflcu(j1,j3)
          END DO
!
#ifdef DEBUG
      IF (ldebug) THEN
        PRINT *,' after backsubst.  in vertical loop j3=',j3
        PRINT *,'pflfd(j1b,j3+1) : ',pflfd(j1b,j3+1)
      END IF
#endif
!
      END DO    ! End of the loop over layers
 
      RETURN
    END SUBROUTINE inv_th

    SUBROUTINE opt_so(prholwc,pdulwc,prhoiwc,pduiwc,                &
     &                  paeq1  ,paeq2 ,paeq3  ,paeq4 , paeq5 ,        &
     &                  pdp    ,papra ,psmu0  ,pqsmu0,                &
     &                  kspec  ,kig1s, kig1e,ki1sc ,ki1ec  ,          &
     &                          ki3sc ,ki3ec  ,ldebug,                &
     &                  podac  ,podaf ,podsc  ,podsf , pbsfc ,pbsff , &
     &                  pusfc  ,pusff  )
!

!**** *opt_so* - calculates the optical properties of the non-gaseous
!                constituents for one spectral
!                interval in the solar part of the spectrum

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
!  Include gme_comconst.h for the gravity constant
! V2_24        2010/04/28 Michael Gertz
!  Adaptions to SVN       
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!=======================================================================
!     Purpose:   Calculation of (absorption and scattering) optical
!                depth and backward scattered fraction for diffuse
!                and upward scattered fraction for direct solar
!                radiation (excluding the contribution by gaseous
!                constituents)

!     METHOD
!     ------

!     - determination of optical properies (i.e. extinction coefficient,
!       single scattering albedo and asymetry factor of the phase function)
!       using approximate relations for Rayleigh effects, cloud water,
!       cloud ice and a combination of five type of aerosols
!
!     - calculation of optical depth (scattering and absorption), back-
!       scattered and upscattered fraction of radiation suitable for use
!       in an implicit delta-two-stream scheme
!     
!==============================================================================

      !     Input data
      INTEGER kig1s, kig1e
      INTEGER ki1sc          ! start index for first  array dimension
      INTEGER ki1ec          ! end   index for first  array dimension
      INTEGER ki3sc          ! start index for third  array dimension
      INTEGER ki3ec          ! end   index for third  array dimension

!     Liquid and ice water density and content within for the cloudy
!     part of each layer

      REAL(wp) :: prholwc(kig1s:kig1e,ki3sc:ki3ec) 
      REAL(wp) :: prhoiwc(kig1s:kig1e,ki3sc:ki3ec) 
      REAL(wp) :: pdulwc (kig1s:kig1e,ki3sc:ki3ec) 
      REAL(wp) :: pduiwc (kig1s:kig1e,ki3sc:ki3ec) 

!     Aerosole contents (optical depths at 0.55 micrometer) for 5 types
!
      REAL(wp) :: paeq1  (kig1s:kig1e,ki3sc:ki3ec) 
      REAL(wp) :: paeq2  (kig1s:kig1e,ki3sc:ki3ec) 
      REAL(wp) :: paeq3  (kig1s:kig1e,ki3sc:ki3ec) 
      REAL(wp) :: paeq4  (kig1s:kig1e,ki3sc:ki3ec) 
      REAL(wp) :: paeq5  (kig1s:kig1e,ki3sc:ki3ec) 

!     pressure thickness of layers and mean layer pressure (TOA on input)
      REAL(wp) :: pdp    (kig1s:kig1e,ki3sc:ki3ec) 
      REAL(wp) :: papra  (kig1s:kig1e) 
 
!     zenith angle and it's inverse
      REAL(wp) :: psmu0  (kig1s:kig1e) 
      REAL(wp) :: pqsmu0 (kig1s:kig1e) 
!
      INTEGER kspec          ! selected spectral interval
      LOGICAL ldebug         ! debug control switch       

!     Output data

      REAL(wp) :: podac (kig1s:kig1e,ki3sc:ki3ec) ! absorption optical depth
      REAL(wp) :: podaf (kig1s:kig1e,ki3sc:ki3ec) ! in cloudy and free part
      REAL(wp) :: podsc (kig1s:kig1e,ki3sc:ki3ec) ! scattering optical depth
      REAL(wp) :: podsf (kig1s:kig1e,ki3sc:ki3ec) ! in cloudy and free part
      REAL(wp) :: pbsfc (kig1s:kig1e,ki3sc:ki3ec) ! backscattering fraction 
      REAL(wp) :: pbsff (kig1s:kig1e,ki3sc:ki3ec) ! in cloudy and free part
      REAL(wp) :: pusfc (kig1s:kig1e,ki3sc:ki3ec) ! upward scattered fraction
      REAL(wp) :: pusff (kig1s:kig1e,ki3sc:ki3ec) ! in cloudy and free part


!     Local arrays and variables

!     optical properties (absorption, scattering, back- and upscatter fraction)
!     for liquid water, ice and total aerosole
!     for Rayleigh scattering, only the optical depth is stored as array, since
!     back- and upscatter fractions are constant (=0.5)

      REAL(wp) :: zlwoda (kig1s:kig1e) !
      REAL(wp) :: zlwods (kig1s:kig1e) ! 
      REAL(wp) :: zlwb0  (kig1s:kig1e) ! 
      REAL(wp) :: zlwb   (kig1s:kig1e) ! 
      REAL(wp) :: ziwoda (kig1s:kig1e) ! 
      REAL(wp) :: ziwods (kig1s:kig1e) ! 
      REAL(wp) :: ziwb0  (kig1s:kig1e) ! 
      REAL(wp) :: ziwb   (kig1s:kig1e) ! 
      REAL(wp) :: zaeoda (kig1s:kig1e) ! 
      REAL(wp) :: zaeods (kig1s:kig1e) ! 
      REAL(wp) :: zaeb0  (kig1s:kig1e) ! 
      REAL(wp) :: zaeb   (kig1s:kig1e) ! 
      REAL(wp) :: zraods (kig1s:kig1e) ! 
      
!     individual optical properties of liquid water and ice 

      REAL(wp) :: z_lwe, z_iwe         ! extinction coefficient
      REAL(wp) :: z_lww, z_iww         ! single scattering coefficient
      REAL(wp) :: z_lwg, z_iwg         ! asymetry factor
      REAL(wp) :: z_lwf, z_iwf         ! forward scattered fraction 

      INTEGER j1,j3                  ! loop indices over spatial dimensions
      INTEGER ja                    ! local loop index
      INTEGER j1b                   ! debug point index


!     Spectral properties (aerosols, liquid water, ice and Rayleigh scattering)

!      REAL(wp) :: zaea(jpspec,5),zaes(jpspec,5),                &
!     &                     zaeg(jpspec,5),zaef(jpspec,5)
!      REAL(wp) :: zlwe(4,jpspec),zlwemn(jpspec),zlwemx(jpspec), &
!     &                     zlww(2,jpspec),zlwg(2,jpspec)
!      REAL(wp) :: ziwe(4,jpspec),ziweMN(jpspec),ziwemx(jpspec), &
!     &                     ziww(2,jpspec),ziwg(2,jpspec)
      REAL(wp) :: zrsc(  jpsol)

!    
!     Local utility variables
!
      REAL(wp) :: z1dg       ! 1./g
      REAL(wp) :: z1d8       ! 1./8
      REAL(wp) :: zzg
      INTEGER isp     ! (=kspec, but shorter for notation purposes)

      REAL(wp) :: ZEPOPD     ! Security constant for optical depth
      REAL(wp) :: ZEPSSA     ! Security constant for single scattering albedo

!      INCLUDE "gme_dat_lw.h"
!      INCLUDE "gme_dat_iw.h"
!      INCLUDE "gme_dat_aer.h"
!      INCLUDE "gme_dat_ray.h"
!      INCLUDE "gme_comconst.h"


!     *zrsc(ispec)* Rayleigh scattering coefficients in solar
!                   spectral intervals
 
      DATA zrsc /0.59776370E-08_wp, 0.13266702E-06_wp,  &
     &           0.20634412E-05_wp/
 

      
      DATA zepopd /1.E-06_wp/
      DATA zepssa /1.E-06_wp/
      DATA j1b / 1/
 
      isp      = kspec
!     z1dg     = 1.0_wp/9.80665_wp
      z1dg     = 1.0_wp/G
      z1d8     = 0.125_wp    
#ifdef DEBUG
      IF (ldebug) THEN
      PRINT *,' **** opt-so   start **********************'
      PRINT *,' **** debug point : ',j1b
      PRINT *,' **** interval    : ',isp    
      END IF
#endif
 
      DO ja=1,5        
      zaef(isp,ja)  = zaeg(isp,ja)**2 ! forward sc.fraction f.aerosols
      END DO   
#ifdef DEBUG
      IF (ldebug) THEN
      DO ja=1,5
      PRINT *,'ja, zaef(isp,ja): ',ja,zaef(isp,ja)
      END DO
      END IF
#endif
!     Vertical loop
 
      DO j3=ki3sc,ki3ec
 
!     Optical properties of liquid water and ice as function of the specific
!     liquid water and ice content
 
          DO j1=ki1sc,ki1ec
 
!         liquid water effects
            
          z_lwg = zlwg(1,isp) + zlwg(2,isp)*prholwc(j1,j3) !T.R.
          z_lwg = MAX(0.0_wp,MIN(1.0_wp,z_lwg))
          z_lwf = z_lwg*z_lwg 
          z_lww = zlww(1,isp) + zlww(2,isp)*prholwc(j1,j3)
          z_lww = MAX(zepssa,MIN(1.0_wp-zepssa,z_lww))
          z_lwe = z1dg * (zlwe(1,isp) + zlwe(2,isp)/          &
     &            (zlwe(3,isp)*prholwc(j1,j3)+zlwe(4,isp)))
          z_lwe = MAX(zlwemn(isp),MIN(zlwemx(isp),z_lwe))
 
          zlwoda(j1) = z_lwe*pdulwc(j1,j3)*(1.0_wp-z_lww)
          zlwods(j1) = z_lwe*pdulwc(j1,j3)*     z_lww *(1._wp-z_lwf)
          zlwb0 (j1) = z1d8*(4._wp+z_lwg)/(1._wp+z_lwg) 
          zlwb  (j1) = 0.5_wp-0.75_wp*psmu0(j1)*z_lwg/(1._wp+z_lwg)
 
!         ice water effects
 
!         z_iwg = ziwg(1,isp) + ziwg(2,isp)*prhoiwc(j1,j3)
          z_iwg = ziwg(1,isp) + ziwg(2,isp)*LOG(prhoiwc(j1,j3))
          z_iwg = MAX(0.0_wp,MIN(1.0_wp,z_iwg))
          z_iwf = z_iwg*z_iwg
!         z_iww = ziww(1,isp) + ziww(2,isp)*prhoiwc(j1,j3)
          z_iww = ziww(1,isp) + ziww(2,isp)*LOG(prhoiwc(j1,j3))
          z_iww = MAX(zepssa,MIN(1.0_wp-zepssa,z_iww))
          z_iwe = z1dg * (ziwe(1,isp) + ziwe(2,isp)/              &
     &               (ziwe(3,isp)*prhoiwc(j1,j3)+ziwe(4,isp)))
          z_iwe = MAX(ziwemn(isp),MIN(ziwemx(isp),z_iwe  ))
!
          ziwoda(j1) = z_iwe*pduiwc(j1,j3)*(1.0_wp-z_iww)
          ziwods(j1) = z_iwe*pduiwc(j1,j3)*     z_iww *(1._wp-z_iwf)
          ziwb0 (j1) = z1d8*(4._wp+z_iwg)/(1._wp+z_iwg)
          ziwb  (j1) = 0.5_wp-0.75_wp*psmu0(j1)*z_iwg/(1._wp+z_iwg)
!          END DO
        END DO
#ifdef DEBUG
      IF (ldebug) THEN
      PRINT *,' prholwc (j1b) :',prholwc (j1b,j3)   
      PRINT *,' pdulwc  (j1b) :',pdulwc  (j1b,j3)
      PRINT *,' zlwoda  (j1b) :',zlwoda  (j1b)
      PRINT *,' zlwods  (j1b) :',zlwods  (j1b)
      PRINT *,' zlwb0   (j1b) :',zlwb0   (j1b)
      PRINT *,' zlwb    (j1b) :',zlwb    (j1b)
      PRINT *,' z_lwg             :',z_lwg                
      PRINT *,' z_lwf             :',z_lwf                
      PRINT *,' z_lww             :',z_lww                
      PRINT *,' z_lwe             :',z_lwe                
      PRINT *,' prhoiwc (j1b) :',prhoiwc (j1b,j3)
      PRINT *,' pduiwc  (j1b) :',pduiwc  (j1b,j3)
      PRINT *,' ziwoda  (j1b) :',ziwoda  (j1b)
      PRINT *,' ziwods  (j1b) :',ziwods  (j1b)
      PRINT *,' ziwb0   (j1b) :',ziwb0   (j1b)
      PRINT *,' ziwb    (j1b) :',ziwb    (j1b)
      END IF
#endif

!     Optical properties of five aerosol types combined
 
          DO j1=ki1sc,ki1ec
          zaeoda(j1)=                                            &
     &    paeq1(j1,j3)*zaea(isp,1)+paeq2(j1,j3)*zaea(isp,2)+  &
     &    paeq3(j1,j3)*zaea(isp,3)+paeq4(j1,j3)*zaea(isp,4)+  &
     &    paeq5(j1,j3)*zaea(isp,5)

          zaeods(j1)=                                            &
     &      ( paeq1(j1,j3)*zaes(isp,1)*(1._wp-zaef(isp,1)) )        &
     &     +( paeq2(j1,j3)*zaes(isp,2)*(1._wp-zaef(isp,2)) )        &
     &     +( paeq3(j1,j3)*zaes(isp,3)*(1._wp-zaef(isp,3)) )        &
     &     +( paeq4(j1,j3)*zaes(isp,4)*(1._wp-zaef(isp,4)) )        &
     &     +( paeq5(j1,j3)*zaes(isp,5)*(1._wp-zaef(isp,5)) )
 
          zzg=(                                                         &
     &     (paeq1(j1,j3)*zaes(isp,1)*(1._wp-zaef(isp,1)))*zaeg(isp,1)   &
     &    +(paeq2(j1,j3)*zaes(isp,2)*(1._wp-zaef(isp,2)))*zaeg(isp,2)   &
     &    +(paeq3(j1,j3)*zaes(isp,3)*(1._wp-zaef(isp,3)))*zaeg(isp,3)   &
     &    +(paeq4(j1,j3)*zaes(isp,4)*(1._wp-zaef(isp,4)))*zaeg(isp,4)   &
     &    +(paeq5(j1,j3)*zaes(isp,5)*(1._wp-zaef(isp,5)))*zaeg(isp,5))  &
     &    / MAX( zaeods(j1),zepopd)

          zaeb0(j1) = z1d8*(4._wp+zzg)/(1._wp+zzg)
          zaeb (j1) = 0.5_wp-0.75_wp*psmu0(j1)*zzg/(1._wp+zzg)

          END DO
#ifdef DEBUG
      IF (ldebug) THEN
      PRINT *,' zaeoda  (j1b) :',zaeoda  (j1b)
      PRINT *,' zaeods  (j1b) :',zaeods  (j1b)
      PRINT *,' zaeb0   (j1b) :',zaeb0   (j1b)
      PRINT *,' zaeb    (j1b) :',zaeb    (j1b)
      END IF
#endif

!     Optical thickness for Rayleigh scattering

          DO j1=ki1sc,ki1ec
          zraods(j1)  = zrsc(isp)*pdp(j1,j3)/(1._wp+zrsc(isp)*   &
     &                   (papra(j1)+0.5_wp*pdp(j1,j3)*pqsmu0(j1)))
          papra(j1)  = papra(j1)+pdp(j1,j3)*pqsmu0(j1)
          END DO
#ifdef DEBUG
      IF (ldebug) THEN
!     PRINT *,' Rayleigh coefficient:',zrsc(isp)
!     PRINT *,' Papra  (j1b)    :',papra   (j1b)
!     PRINT *,' Pqsmu0 (j1b)    :',pqsmu0  (j1b)
!     PRINT *,' Pdp    (j1b)    :',pdp     (j1b,j3)
      PRINT *,' zraods (j1b)    :',zraods  (j1b)
      END IF
#endif
!     ------------------------------------------------------------------
 
!     Linear combination of individual contributions 
 
!     a) cloud free part of layer
 
          DO j1=ki1sc,ki1ec
          podaf(j1,j3)=MAX( zaeoda(j1)                ,zepopd)
          podsf(j1,j3)=MAX( zaeods(j1) + zraods(j1),zepopd)
          podsf(j1,j3)=MIN( podsf (j1,j3),                         &
     &            (1._wp-zepssa)*(podaf(j1,j3)+podsf(j1,j3)))
 
          pbsff(j1,j3)=(zaeb0(j1)*zaeods(j1)+0.5_wp*zraods(j1)) &
     &                             / podsf(j1,j3)
          pusff(j1,j3)=(zaeb (j1)*zaeods(j1)+0.5_wp*zraods(j1)) &
     &                             / podsf(j1,j3)
!         pbsff(j1,j3)= MAX (pbsff(j1,j3),0.00001)
 
!     b) cloudy part of layer
 
          podac(j1,j3)= zlwoda(j1)+ziwoda(j1)+zaeoda(j1)
          podsc(j1,j3)= zlwods(j1)+ziwods(j1)+zaeods(j1)  &
     &                   + zraods(j1)
          podac(j1,j3)= MAX(podac(j1,j3), zepopd)
          podsc(j1,j3)= MAX(podsc(j1,j3), zepopd)
          podsc(j1,j3)=MIN( podsc (j1,j3),                      &
     &            (1._wp-zepssa)*(podac(j1,j3)+podsc(j1,j3)))
 
          pbsfc(j1,j3)= (zlwb0(j1)*zlwods(j1)   &
     &                     +ziwb0(j1)*ziwods(j1)   &
     &                     +zaeb0(j1)*zaeods(j1)   &
     &                     +  0.5_wp*zraods(j1))  &
     &                       / podsc(j1,j3)
          pusfc(j1,j3)= (zlwb (j1)*zlwods(j1)   &
     &                     +ziwb (j1)*ziwods(j1)   &
     &                     +zaeb (j1)*zaeods(j1)   &
     &                     +  0.5_wp*zraods(j1))  &
     &                       / podsc(j1,j3)
!         pbsfc(j1,j3)= MAX (pbsfc(j1,j3),0.00001)
          END DO
 
      END DO    ! End of vertical loop
#ifdef DEBUG
      IF (ldebug) THEN
      PRINT *,' podaf   (j1b,1) :',podaf   (j1b,1)
      PRINT *,' podsf   (j1b,1) :',podsf   (j1b,1)
      PRINT *,' pbsff   (j1b,1) :',pbsff   (j1b,1)
      PRINT *,' pusff   (j1b,1) :',pusff   (j1b,1)
      PRINT *,' podac   (j1b,1) :',podac   (j1b,1)
      PRINT *,' podsc   (j1b,1) :',podsc   (j1b,1)
      PRINT *,' pbsfc   (j1b,1) :',pbsfc   (j1b,1)
      PRINT *,' pusfc   (j1b,1) :',pusfc   (j1b,1)
      PRINT *,' -------------------------------------------------'
      END IF
#endif
 
      RETURN
    END SUBROUTINE opt_so

   SUBROUTINE inv_so (                                               &
     &  pclc   ,pca1   ,pca2  ,pcb1  ,pcb2  ,pcc1  ,pcc2 ,pcd1 ,pcd2 ,  &
     &  pflpt  ,psmu0  ,pqsmu0,palp  ,palso ,                           &
     &  pduh2oc,pduh2of,pduco2,pduo3 ,palogp,palogt,                    &
     &  podsc  ,podsf  ,podac ,podaf ,pbsfc ,pbsff ,pusfc,pusff,        &
     &  kspec  ,kh2o   ,kco2  ,ko3   ,                                  &
     &  kig1s  ,kig1e  ,              ki3s ,ki3e ,                      &
     &  ki1sc  ,ki1ec  ,                                                &
     &  ldebug ,                                                        &
     &  pflcu  ,pflfu  ,pflcd ,pflfd ,pflcp ,pflfp)
!
 
!     *inv_th* solves the linear system of equations for solar fluxes
 
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
!=======================================================================
 
!     Purpose
!     -------
!
!     The routine solves a linear equation system for solar fluxes using
!     a Gaussian elimination-backsubstitution algorithm dedicated to the
!     specific structure of the system matrix
 
!     Method
!     ------
!
!     - setting of the RHS of the system using the parallel solar radiation
!       at the top of the atmosphere and allowing for partial cloud cover
!     - solution of the equation system including the lower boundary
!       condition
!     - matrix coefficients are calculated in the course of the elimination
!       step for one layer at a time through a call to routine *coe_so*
!     - the final result, i.e. upward and downward diffuse and parallel 
!       solar fluxes are stored seperately for cloudy and cloud-free parts
!       of each layer boundary
!==============================================================================

!     Input data 
 
      INTEGER kig1s,kig1e,ki3s,ki3e   ! global array bounds
      INTEGER ki1sc,ki1ec             ! computational domain

      REAL(wp) :: pclc  (kig1s:kig1e,ki3s:ki3e) ! cloud cover
      REAL(wp) :: pca1  (kig1s:kig1e,ki3s:ki3e) ! cloud geometry factor  
      REAL(wp) :: pca2  (kig1s:kig1e,ki3s:ki3e) ! cloud geometry factor  
      REAL(wp) :: pcb1  (kig1s:kig1e,ki3s:ki3e) ! cloud geometry factor  
      REAL(wp) :: pcb2  (kig1s:kig1e,ki3s:ki3e) ! cloud geometry factor  
      REAL(wp) :: pcc1  (kig1s:kig1e,ki3s:ki3e) ! cloud geometry factor  
      REAL(wp) :: pcc2  (kig1s:kig1e,ki3s:ki3e) ! cloud geometry factor  
      REAL(wp) :: pcd1  (kig1s:kig1e,ki3s:ki3e) ! cloud geometry factor  
      REAL(wp) :: pcd2  (kig1s:kig1e,ki3s:ki3e) ! cloud geometry factor  

      REAL(wp) :: pflpt (kig1s:kig1e) ! parallel solar flux at TOA
      REAL(wp) :: palp  (kig1s:kig1e) ! surface albedo for parallel
      REAL(wp) :: palso (kig1s:kig1e) ! and for diffuse radiation  
 
!     Input arrays to be passed to *coe_so*
 
!     layer gas contents (cloudy and cloud-free, if distinction necessary)
      REAL(wp) :: pduh2oc(kig1s:kig1e,ki3s:ki3e) ! h2o-vapour cloudy      
      REAL(wp) :: pduh2of(kig1s:kig1e,ki3s:ki3e) ! h2o-vapour cloud-free  
      REAL(wp) :: pduco2 (kig1s:kig1e,ki3s:ki3e) ! co2
      REAL(wp) :: pduo3  (kig1s:kig1e,ki3s:ki3e) ! o3 
!     optical properties of 'grey' constituents (cloudy and cloud-free)
      REAL(wp) :: podsc  (kig1s:kig1e,ki3s:ki3e) ! scattering optical depth
      REAL(wp) :: podsf  (kig1s:kig1e,ki3s:ki3e) ! scattering optical depth
      REAL(wp) :: podac  (kig1s:kig1e,ki3s:ki3e) ! absorption optical depth
      REAL(wp) :: podaf  (kig1s:kig1e,ki3s:ki3e) ! absorption optical depth
      REAL(wp) :: pbsfc  (kig1s:kig1e,ki3s:ki3e) ! backscatter fraction
      REAL(wp) :: pbsff  (kig1s:kig1e,ki3s:ki3e) ! backscatter fraction
      REAL(wp) :: pusfc  (kig1s:kig1e,ki3s:ki3e) ! upscatter   fraction
      REAL(wp) :: pusff  (kig1s:kig1e,ki3s:ki3e) ! upscatter   fraction
!     
      REAL(wp) :: palogp (kig1s:kig1e,ki3s:ki3e) ! ln(p)
      REAL(wp) :: palogt (kig1s:kig1e,ki3s:ki3e) ! ln(T)
!
      REAL(wp) :: psmu0 (kig1s:kig1e) ! cosine of zenith angle
      REAL(wp) :: pqsmu0(kig1s:kig1e) ! 1._wp/cosine of zenith angle

!
      INTEGER kspec          ! selected spectral interval
      INTEGER kh2o           ! table index for h2o-vapour absorption
      INTEGER kco2           ! table index for co2        absorption
      INTEGER ko3            ! table index for o3         absorption
      LOGICAL ldebug         ! debug control switch       

!     Output data
 
      REAL(wp) :: pflcu (kig1s:kig1e,ki3s:ki3e+1) ! flux up   cloudy
      REAL(wp) :: pflfu (kig1s:kig1e,ki3s:ki3e+1) ! flux up   cloud-free 
      REAL(wp) :: pflcd (kig1s:kig1e,ki3s:ki3e+1) ! flux down cloudy     
      REAL(wp) :: pflfd (kig1s:kig1e,ki3s:ki3e+1) ! flux down cloud-free 
      REAL(wp) :: pflcp (kig1s:kig1e,ki3s:ki3e+1) ! flux par. cloudy     
      REAL(wp) :: pflfp (kig1s:kig1e,ki3s:ki3e+1) ! flux par. cloud-free 

!     Local arrays and variables

!     layer properties calculated in *coe_so*

      REAL(wp) :: pa1c   (kig1s:kig1e) !
      REAL(wp) :: pa1f   (kig1s:kig1e) ! 
      REAL(wp) :: pa2c   (kig1s:kig1e) ! 
      REAL(wp) :: pa2f   (kig1s:kig1e) ! 
      REAL(wp) :: pa3c   (kig1s:kig1e) ! 
      REAL(wp) :: pa4f   (kig1s:kig1e) ! 
      REAL(wp) :: pa4c   (kig1s:kig1e) ! 
      REAL(wp) :: pa5f   (kig1s:kig1e) ! 
      REAL(wp) :: pa5c   (kig1s:kig1e) ! 
      REAL(wp) :: pa3f   (kig1s:kig1e) ! 

!     Utility arrays
      REAL(wp) :: ztu1 (kig1s:kig1e,ki3s:ki3e) !
      REAL(wp) :: ztu2 (kig1s:kig1e,ki3s:ki3e) !
      REAL(wp) :: ztu3 (kig1s:kig1e,ki3s:ki3e) !
      REAL(wp) :: ztu4 (kig1s:kig1e,ki3s:ki3e) !
      REAL(wp) :: ztu5 (kig1s:kig1e,ki3s:ki3e) !
      REAL(wp) :: ztu6 (kig1s:kig1e,ki3s:ki3e) !
      REAL(wp) :: ztu7 (kig1s:kig1e,ki3s:ki3e) !
      REAL(wp) :: ztu8 (kig1s:kig1e,ki3s:ki3e) !
      REAL(wp) :: ztu9 (kig1s:kig1e,ki3s:ki3e) !

!     Utility variables
      REAL(wp) :: ztd1 ,ztd2 ,ztd3 ,ztd4 ,ztd5 ,ztd6 , ztd7
      REAL(wp) :: ztds1,ztds2,ztds3,ztus1
      
      INTEGER j1,j3       ! loop indices over spatial dimensions
      INTEGER j1b        ! index of debug-point
      LOGICAL ldebug_coe_so  ! debug switch for *coe_so*
      DATA j1b /1/
!     ------------------------------------------------------------------
      ldebug_coe_so = .false. 
!
!
!     Upper boundary condition
 
        DO j1=ki1sc,ki1ec
        pflfp(j1,ki3s) = PFLPT(j1)
        pflcp(j1,ki3s) = 0._wp
        pflfd(j1,ki3s) = 0._wp
        pflcd(j1,ki3s) = 0._wp
        END DO
#ifdef DEBUG
      IF (ldebug) THEN
      PRINT *,' *** INV_SO **************************'
      PRINT *,' *** Debug point: ',j1b
      PRINT *,'pflfp(j1b,ki3s) : ',pflfp(j1b,ki3s)
      PRINT *,'pflcp(j1b,ki3s) : ',pflcp(j1b,ki3s)
      PRINT *,'pflfp(j1b,ki3s) : ',pflfp(j1b,ki3s)
      PRINT *,'pflcd(j1b,ki3s) : ',pflcd(j1b,ki3s)
      END IF
#endif

!     Determine effects of first layer in *coe_so*
      CALL       coe_so (                                               &
     & pduh2oc,pduh2of,pduco2 ,pduo3  ,palogp ,palogt ,                 &
     & podsc  ,podsf  ,podac  ,podaf  ,pbsfc  ,pbsff  ,pusfc ,pusff ,   &
     & psmu0  ,pqsmu0 ,                                                 &
     & ki3s  ,kspec  ,kh2o   ,kco2   ,ko3    ,                          &
     & kig1s ,kig1e, ki1sc  ,ki1ec ,                ki3s  ,ki3e,        &
     & ldebug_coe_so ,                                                  &
     & pa1c   ,pa1f   ,pa2c   ,pa2f   ,pa3c   ,pa3f ,                   &
     & pa4c   ,pa4f   ,pa5c   ,pa5f )
 
!     Top layer elimination
 
        DO j1=ki1sc,ki1ec
        pflfu(j1,ki3s  )= pa3f(j1)                         &
     &                       *pca2(j1,ki3s)*pflfp(j1,ki3s)
        pflfp(j1,ki3s+1)= pa1f(j1)                         &
     &                       *pca2(j1,ki3s)*pflfp(j1,ki3s)
        pflfd(j1,ki3s+1)= pa2f(j1)                         &
     &                       *pca2(j1,ki3s)*pflfp(j1,ki3s)
        pflcu(j1,ki3s  )= pa3c(j1)                         &
     &                       *pcb2(j1,ki3s)*pflfp(j1,ki3s)
        pflcp(j1,ki3s+1)= pa1c(j1)                         &
     &                       *pcb2(j1,ki3s)*pflfp(j1,ki3s)
        pflcd(j1,ki3s+1)= pa2c(j1)                         &
     &                       *pcb2(j1,ki3s)*pflfp(j1,ki3s)
        END DO
#ifdef DEBUG
      IF (ldebug) THEN
      PRINT *,' *** INV_SO **************************'
      PRINT *,'pflfu(j1b,ki3s)  : ',pflfu(j1b,ki3s)
      PRINT *,'pflcu(j1b,ki3s)  : ',pflcu(j1b,ki3s)
      PRINT *,'pflfd(j1b,ki3s+1): ',pflfd(j1b,ki3s+1)
      PRINT *,'pflcd(j1b,ki3s+1): ',pflcd(j1b,ki3s+1)
      PRINT *,'pa1f (j1b)        : ',pa1f (j1b)        
      PRINT *,'pa1c (j1b)        : ',pa1c (j1b)        
      PRINT *,'pa2f (j1b)        : ',pa2f (j1b)        
      PRINT *,'pa2c (j1b)        : ',pa2c (j1b)        
      PRINT *,'pa3f (j1b)        : ',pa3f (j1b)        
      PRINT *,'pa3c (j1b)        : ',pa3c (j1b)        
!     STOP 'test'
      END IF
#endif
 
!     Storage of utility arrays for the top layer
 
        DO j1=ki1sc,ki1ec
        ztu1(j1,1) =0._wp
        ztu2(j1,1) = pca1(j1,1)*pa4f(j1)
        ztu3(j1,1) = pcc1(j1,1)*pa4f(j1)
        ztu4(j1,1) = pcb1(j1,1)*pa4c(j1)
        ztu5(j1,1) = pcd1(j1,1)*pa4c(j1)
        ztu6(j1,1) = pca1(j1,1)*pa5f(j1)
        ztu7(j1,1) = pcc1(j1,1)*pa5f(j1)
        ztu8(j1,1) = pcb1(j1,1)*pa5c(j1)
        ztu9(j1,1) = pcd1(j1,1)*pa5c(j1)
        END DO
 
!     Suczessive layer-by-layer elimination
 
      DO j3=ki3s+1,ki3e         ! Loop over vertical

!     Determine effects of layer in *coe_so*
      CALL       coe_so (                                              &
     & pduh2oc,pduh2of,pduco2 ,pduo3  ,palogp ,palogt ,                &
     & podsc  ,podsf  ,podac  ,podaf  ,pbsfc  ,pbsff  ,pusfc ,pusff ,  &
     & psmu0  ,pqsmu0 ,                                                &
     & j3     ,kspec  ,kh2o   ,kco2   ,ko3    ,                        &
     & kig1s ,kig1e,  ki1sc  ,ki1ec  , ki3s  ,ki3e,                    &
     & ldebug_coe_so ,                                                 &
     & pa1c   ,pa1f   ,pa2c   ,pa2f   ,pa3c   ,pa3f ,                  &
     & pa4c   ,pa4f   ,pa5c   ,pa5f )
 
!     Elimination
 
          DO j1=ki1sc,ki1ec
          ztd1 = 1._wp/(1._wp-pa5f(j1)*(pca2(j1,j3)*ztu6(j1,j3-1)  &
     &                                            +pcc2(j1,j3)*ztu8(j1,j3-1)))
          pflfu(j1,j3)= ztd1*(                             &
     &     pa3f(j1)*( pca2(j1,j3)*pflfp(j1,j3)       &
     &                  +pcc2(j1,j3)*pflcp(j1,j3) )     &
     &    +pa5f(j1)*( pca2(j1,j3)*pflfd(j1,j3)       &
     &                  +pcc2(j1,j3)*pflcd(j1,j3) )   )
          ztu1 (j1,j3)= ztd1*pa5f(j1)*                  &
     &                 ( pca2(j1,j3)*ztu7(j1,j3-1)      &
     &                  +pcc2(j1,j3)*ztu9(j1,j3-1))
          ztu2(j1,j3) =  ztd1*pa4f(j1)*pca1(j1,j3)
          ztu3(j1,j3) =  ztd1*pa4f(j1)*pcc1(j1,j3)
          ztd2 = pa5c(j1)*( pcb2(j1,j3)*ztu6(j1,j3-1)          &
     &                        +pcd2(j1,j3)*ztu8(j1,j3-1) )
          ztd3 = 1._wp/( 1._wp-pa5c(j1)*(pcb2(j1,j3)*ztu7(j1,j3-1)   &
     &                                             +pcd2(j1,j3)*ztu9(j1,j3-1))  &
     &                                -ztd2*ztu1(j1,j3) )
          pflcu(j1,j3)= ztd3 *(                             &
     &     pa3c(j1)*( pcb2(j1,j3)*pflfp(j1,j3)        &
     &                  +pcd2(j1,j3)*pflcp(j1,j3) )      &
     &    +pa5c(j1)*( pcb2(j1,j3)*pflfd(j1,j3)        &
     &                  +pcd2(j1,j3)*pflcd(j1,j3) )      &
     &    +ztd2*pflfu(j1,j3)   )
          ztu4(j1,j3) = ztd3 *( pa4c(j1)*pcb1(j1,j3)  &
     &                            +ztd2       *ztu2(j1,j3) )
          ztu5(j1,j3) = ztd3 *( pa4c(j1)*pcd1(j1,j3)  &
     &                            +ztd2       *ztu3(j1,j3) )
          pflfp(j1,j3+1)=pa1f(j1)*(pca2(j1,j3)*pflfp(j1,j3)  &
     &                                  +pcc2(j1,j3)*pflcp(j1,j3))
          pflcp(j1,j3+1)=pa1c(j1)*(pcb2(j1,j3)*pflfp(j1,j3)  &
     &                                  +pcd2(j1,j3)*pflcp(j1,j3))
          ztd4 = pa4f(j1)*( pca2(j1,j3)*ztu6(j1,j3-1)           &
     &                        +pcc2(j1,j3)*ztu8(j1,j3-1) )
          ztd5 = pa4f(j1)*( pca2(j1,j3)*ztu7(j1,j3-1)           &
     &                        +pcc2(j1,j3)*ztu9(j1,j3-1) )
          pflfd(j1,j3+1)=pa2f(j1)*(pca2(j1,j3)*pflfp(j1,j3)  &
     &                                  +pcc2(j1,j3)*pflcp(j1,j3)) &
     &                     +pa4f(j1)*(pca2(j1,j3)*pflfd(j1,j3)  &
     &                                  +pcc2(j1,j3)*pflcd(j1,j3)) &
     &                     +ztd4*pflfu(j1,j3)                         &
     &                     +ztd5*pflcu(j1,j3)
          ztu6(j1,j3) = pa5f(j1)*pca1(j1,j3)  &
     &                    +ztd4*ztu2(j1,j3)         &
     &                    +ztd5*ztu4(j1,j3)
          ztu7(j1,j3) = pa5f(j1)*pcc1(j1,j3)  &
     &                    +ztd4*ztu3(j1,j3)         &
     &                    +ztd5*ztu5(j1,j3)
          ztd6 = pa4c(j1)*( pcb2(j1,j3)*ztu6(j1,j3-1)  &
     &                        +pcd2(j1,j3)*ztu8(j1,j3-1) )
          ztd7 = pa4c(j1)*( pcb2(j1,j3)*ztu7(j1,j3-1)  &
     &                        +pcd2(j1,j3)*ztu9(j1,j3-1) )
          pflcd(j1,j3+1)=pa2c(j1)*(pcb2(j1,j3)*pflfp(j1,j3)  &
     &                                  +pcd2(j1,j3)*pflcp(j1,j3)) &
     &                     +pa4c(j1)*(pcb2(j1,j3)*pflfd(j1,j3)  &
     &                                  +pcd2(j1,j3)*pflcd(j1,j3)) &
     &                     +ztd6*pflfu(j1,j3)                         &
     &                     +ztd7*pflcu(j1,j3)
          ztu8(j1,j3) = pa5c(j1)*pcb1(j1,j3)  &
     &                    +ztd6*ztu2(j1,j3)         &
     &                    +ztd7*ztu4(j1,j3)
          ztu9(j1,j3) = pa5c(j1)*pcd1(j1,j3)  &
     &                    +ztd6*ztu3(j1,j3)         &
     &                    +ztd7*ztu5(j1,j3)
          END DO
#ifdef DEBUG
      IF (ldebug) THEN
      PRINT *,' inv_so j3=',j3
      PRINT *,'pflfu(j1b,j3)  : ',pflfu(j1b,j3)
      PRINT *,'pflcu(j1b,j3)  : ',pflcu(j1b,j3)
      PRINT *,'pflfd(j1b,j3+1): ',pflfd(j1b,j3+1)
      PRINT *,'pflcd(j1b,j3+1): ',pflcd(j1b,j3+1)
      PRINT *,'pflfp(j1b,j3+1): ',pflfp(j1b,j3+1)
      PRINT *,'pflcp(j1b,j3+1): ',pflcp(j1b,j3+1)
      PRINT *,'ztu1 (j1b,j3)  : ',ztu1 (j1b,j3)
      PRINT *,'ztu2 (j1b,j3)  : ',ztu2 (j1b,j3)
      PRINT *,'ztu3 (j1b,j3)  : ',ztu3 (j1b,j3)
      PRINT *,'ztu4 (j1b,j3)  : ',ztu4 (j1b,j3)
      PRINT *,'ztu5 (j1b,j3)  : ',ztu5 (j1b,j3)
      PRINT *,'ztu6 (j1b,j3)  : ',ztu6 (j1b,j3)
      PRINT *,'ztu7 (j1b,j3)  : ',ztu7 (j1b,j3)
      PRINT *,'ztu8 (j1b,j3)  : ',ztu8 (j1b,j3)
      PRINT *,'ztu9 (j1b,j3)  : ',ztu9 (j1b,j3)
      PRINT *,' .....'
      END IF
#endif
      END DO       ! Vertical loop
 
!     Elimination and back-substitution at surface
 
        DO j1=ki1sc,ki1ec
        ztds1  =1._wp/(1._wp-palso(j1)*ztu6(j1,ki3e))
        pflfu(j1,ki3e+1)= ztds1 *(palp (j1)*pflfp(j1,ki3e+1)   &
     &                              +palso(j1)*pflfd(j1,ki3e+1))
        ztus1  =               ztds1*palso(j1)*ztu7(j1,ki3e)
        ztds2  =                     palso(j1)*ztu8(j1,ki3e)
        ztds3  =1._wp/(1._wp-palso(j1)*ztu9(j1,ki3e)-ztds2*ztus1)
        pflcu(j1,ki3e+1)= ztds3 *(palp (j1)*pflcp(j1,ki3e+1)   &
     &                              +palso(j1)*pflcd(j1,ki3e+1)   &
     &                              +ztds2       *pflfu(j1,ki3e+1))
        pflfu(j1,ki3e+1)=       pflfu(j1,ki3e+1)                  & 
     &                      +ztus1*pflcu(j1,ki3e+1)
        END DO
#ifdef DEBUG
      IF (ldebug) THEN
      PRINT *,' inv_so surface ------------------------------'
      PRINT *,'pflfu(j1b,ki3e+1) : ',pflfu(j1b,ki3e+1)
      PRINT *,'pflcu(j1b,ki3e+1) : ',pflcu(j1b,ki3e+1)
      PRINT *,'palp (j1b): ',palp (j1b)     
      PRINT *,'palso(j1b): ',palso(j1b)     
      PRINT *,'ztds1               ',ztds1                  
      PRINT *,'ztds2               ',ztds2                  
      PRINT *,'ztds3               ',ztds3                  
      PRINT *,'ztus1               ',ztus1                  
      PRINT *,' .....'
      END IF
#endif

!     Layer-by-layer backsubstitution

#ifdef DEBUG
      IF (ldebug) THEN
      PRINT *,' inv_so BACKSUBSTITUTION'
      END IF
#endif
      DO j3=ki3e,ki3s,-1
          DO j1=ki1sc,ki1ec
          pflcd(j1,j3+1) =                pflcd(j1,j3+1)  &
     &                       +ztu8(j1,j3)*pflfu(j1,j3+1)  &
     &                       +ztu9(j1,j3)*pflcu(j1,j3+1)

          pflfd(j1,j3+1) =                pflfd(j1,j3+1)  &
     &                       +ztu6(j1,j3)*pflfu(j1,j3+1)  &
     &                       +ztu7(j1,j3)*pflcu(j1,j3+1)

          pflcu(j1,j3)   =                pflcu(j1,j3)    &
     &                       +ztu4(j1,j3)*pflfu(j1,j3+1)  &
     &                       +ztu5(j1,j3)*pflcu(j1,j3+1)

          pflfu(j1,j3)   =                pflfu(j1,j3)    &
     &                       +ztu2(j1,j3)*pflfu(j1,j3+1)  &
     &                       +ztu3(j1,j3)*pflcu(j1,j3+1)  &
     &                       +ztu1(j1,j3)*pflcu(j1,j3)
          END DO
#ifdef DEBUG
      IF (ldebug) THEN
      PRINT *,' inv_so j3=',j3
      PRINT *,'pflfu(j1b,j3)  : ',pflfu(j1b,j3)
      PRINT *,'pflcu(j1b,j3)  : ',pflcu(j1b,j3)
      PRINT *,'pflfd(j1b,j3+1): ',pflfd(j1b,j3+1)
      PRINT *,'pflcd(j1b,j3+1): ',pflcd(j1b,j3+1)
      END IF
#endif
      END DO  
 
      RETURN
    END SUBROUTINE inv_so

    SUBROUTINE coe_th (                                 &
       pduh2oc,pduh2of,pduco2 ,pduo3  ,palogp ,palogt ,   &
       podsc  ,podsf  ,podac  ,podaf  ,pbsfc  ,pbsff  ,   &
       ki3    ,kspec  ,kh2o   ,kco2   ,ko3    ,           &
       kig1s ,kig1e,ki1sc  ,ki1ec  ,ki3sc  ,ki3ec,     &
       ldebug ,                                           &
       pa1c   ,pa1f   ,pa2c   ,pa2f   ,pa3c   ,pa3f)
!
!**** *coe_th* - calculates the optical effects of atmospheric
!                layers on thermal radiation based on basic
!                optical properties of non-gaseous constituents
!                and gaseous absorption coefficients selected
!                through the corresponding control variables
!                in the argument list

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
!  zargli, ztsec given as constants (same values as in coe_so.f90)
! 1.17       2002/05/03 A. Mueller 
!  Minor correction
! 1.22       2003/02/12 A. Mueller 
!  Avoid one EXP calculation
! V2_20        2009/04/06 H. Frank
!  Restructuring for NEC SX. Replace ztsec by exp(-zargli)
! V2_24        2010/04/28 Michael Gertz
!  Adaptions to SVN
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!=======================================================================
!     Purpose
!     -----
!
!     This routine computes layer effects (transmissivity, reflectivity
!     and emmisivity) in the thermal part of the radiative spectrum
!     both for the cloud-free and the cloudy part of a model layer.
!     The calculation is based on the implicit delt-two-stream equations
!     (cf. Ritter and Geleyn, 1992) and uses basic optical properties
!     (i.e. absorption and scattering optical depth and backscattered
!     fraction for non-gaseous atmospheric constituents as well as 
!     gaseous absorption properties) as input. 

!     Method
!     ------
!
!     - addition of individual gaseous absorption effects to the optical
!       properties of the non-gaseous constituents
!     - determination of layer effects (cf. Zdunkowski et al., 1982, 1986
!       and Ritter and Geleyn, 1992)
!==============================================================================

!     Input data

!     opticall relevant gas quantities (Pa)

      INTEGER ki3         ! vertical layer considered
      INTEGER kspec       ! spectral interval considered
      INTEGER kh2o        ! table index for h2o absorption properties
      INTEGER kco2        ! table index for co2 absorption properties
      INTEGER ko3         ! table index for o3  absorption properties

      !
      INTEGER kig1s ,kig1e
      INTEGER ki1sc          ! start index for first  array dimension
      INTEGER ki1ec          ! end   index for first  array dimension
      INTEGER ki3sc          ! start index for third  array dimension
      INTEGER ki3ec          ! end   index for third  array dimension

      REAL(wp) :: pduh2oc(kig1s:kig1e,ki3sc:ki3ec) ! h2o inside cloud
      REAL(wp) :: pduh2of(kig1s:kig1e,ki3sc:ki3ec) ! h2o out of cloud
      REAL(wp) :: pduco2 (kig1s:kig1e,ki3sc:ki3ec) ! co2 content 
      REAL(wp) :: pduo3  (kig1s:kig1e,ki3sc:ki3ec) ! o3  content 

!     Logarithm of layer mean temperature and pressure

      REAL(wp) :: palogt (kig1s:kig1e,ki3sc:ki3ec) ! ln T
      REAL(wp) :: palogp (kig1s:kig1e,ki3sc:ki3ec) ! ln p
    
!     Optical properties of non-gaseous constituents (..c=cloudy; ..f=free)  

      REAL(wp) :: podsc  (kig1s:kig1e,ki3sc:ki3ec) ! 
      REAL(wp) :: podsf  (kig1s:kig1e,ki3sc:ki3ec) ! 
      REAL(wp) :: podac  (kig1s:kig1e,ki3sc:ki3ec) ! 
      REAL(wp) :: podaf  (kig1s:kig1e,ki3sc:ki3ec) ! 
      REAL(wp) :: pbsfc  (kig1s:kig1e,ki3sc:ki3ec) ! 
      REAL(wp) :: pbsff  (kig1s:kig1e,ki3sc:ki3ec) ! 
 
!
      LOGICAL ldebug         ! debug control switch       

!     Output data

      REAL(wp) :: pa1c  (kig1s:kig1e) ! transmissivity in cloud
      REAL(wp) :: pa1f  (kig1s:kig1e) ! transmissivity cloud-free  
      REAL(wp) :: pa2c  (kig1s:kig1e) ! reflectivity   in cloud 
      REAL(wp) :: pa2f  (kig1s:kig1e) ! reflectivity   cloud-free 
      REAL(wp) :: pa3c  (kig1s:kig1e) ! emissivity     in cloud
      REAL(wp) :: pa3f  (kig1s:kig1e) ! emissivity     cloud-free


!     Local arrays and variables

      REAL(wp) :: zodgf  (kig1s:kig1e) !
      REAL(wp) :: zodgc  (kig1s:kig1e) !
      REAL(wp) :: zod1f  (kig1s:kig1e) ! 
      REAL(wp) :: zod1c  (kig1s:kig1e) ! 
      REAL(wp) :: zod2f  (kig1s:kig1e) ! 
      REAL(wp) :: zod2c  (kig1s:kig1e) ! 
      
!     Security parameters (machine dependent) and utility variables

      REAL(wp), PARAMETER :: zargli = 70.0_wp         ! argument limit for EXP 
      REAL(wp), PARAMETER :: zodmax = 1.E+06_wp ! maximum allowed optical depth
      REAL(wp) :: zeps                         !
      REAL(wp) :: ztau                         !
      REAL(wp) :: zrho                         !
      REAL(wp) :: ztmp                         !

!     Diffusivity factors for gases and other constituents

      REAL(wp), PARAMETER :: zudiff = 2.0_wp     ! gases
      REAL(wp) ::            zangfa

      INTEGER :: j1,j3     ! loop indices over spatial dimensions
!
!      INCLUDE "gme_com_gas.h"
! 
      INTEGER :: j1b        ! index of debug-point
      j1b = ki1sc

      zangfa = EXP( 0.5_wp)
!
!=======================================================================
!
      j3     = ki3

      IF (ldebug) THEN
        PRINT *,'**** coe_th ******************************'
        PRINT *,'**** debug point : ',j1b
        PRINT *,'**** coe_th kspec=',kspec
        PRINT *,'**** coe_th j3   =',j3   
        PRINT *,'**** coe_th kh2o =',kh2o 
        PRINT *,'**** coe_th kco2 =',kco2 
        PRINT *,'**** coe_th ko3  =',ko3  
        PRINT *,'**** pduh2of(j1b,j3)=',pduh2of(ki1sc,j3)
        PRINT *,'**** pduh2oc(j1b,j3)=',pduh2oc(ki1sc,j3)
        PRINT *,'**** pduco2 (j1b,j3)=',pduco2 (ki1sc,j3)
        PRINT *,'**** pduo3  (j1b,j3)=',pduo3  (ki1sc,j3)
      END IF

!     Optical depth of gases

        DO j1=ki1sc,ki1ec
          zodgf(j1) = 0.0_wp     ! Initialisation
        END DO
 
      IF ( kco2 >0) THEN     ! Include CO2 contribution
          DO j1=ki1sc,ki1ec
            zodgf(j1)=zodgf(j1) + pduco2(j1,j3)*    &
                   ( cobi(kco2,kspec,2)*                     &
               EXP (coali(kco2,kspec,2)*palogp(j1,j3)     &
                   -cobti(kco2,kspec,2)*palogt(j1,j3)))
          END DO
      END IF                  ! CO2
      IF (ldebug) PRINT *,'**** zodgf(CO2)        =',zodgf(j1b)
 
      IF ( ko3 >0) THEN  ! Include O3 contribution
          DO j1=ki1sc,ki1ec
            zodgf(j1)=zodgf(j1) + pduo3 (j1,j3)*    &
                   ( cobi(ko3 ,kspec,3)*                     &
               EXP (coali(ko3 ,kspec,3)*palogp(j1,j3)     &
                   -cobti(ko3 ,kspec,3)*palogt(j1,j3)))
          END DO
      END IF
      IF (ldebug) PRINT *,'**** zodgf(CO2+O3)     =',zodgf(j1b)
 
!     Cloudy = cloud free for CO2 and O3 :
 
        DO j1=ki1sc,ki1ec
          zodgc(j1) = zodgf(j1)
        END DO
 
      IF ( kh2o >0) THEN  ! Include H2O contribution
          DO j1=ki1sc,ki1ec
            ztmp = (cobi(kh2o,kspec,1)*                    &
              EXP (coali(kh2o,kspec,1)*palogp(j1,j3)    &
                  -cobti(kh2o,kspec,1)*palogt(j1,j3)))
            zodgf(j1)=zodgf(j1) + pduh2of(j1,j3)*ztmp
            zodgc(j1)=zodgc(j1) + pduh2oc(j1,j3)*ztmp
          END DO
      END IF
      IF (ldebug) THEN
        PRINT *,'**** zodgf(CO2+O3+H2O) =',zodgf(j1b)
        PRINT *,'**** zodgc(CO2+O3+H2O) =',zodgc(j1b)
      END IF
 
        DO j1=ki1sc,ki1ec
          zodgf (j1)=MIN(zodgf(j1),zodmax)
          zodgc (j1)=MIN(zodgc(j1),zodmax)
        END DO
! 
!     Pseudo-optical depth:
! 
        DO j1=ki1sc,ki1ec
!         in cloud-free part of layer
          zod2f(j1) = zudiff*pbsff(j1,j3)*podsf(j1,j3)
          zod1f(j1) = zod2f(j1) + zudiff*podaf(j1,j3)
          zod1f(j1) = zod1f(j1) + zangfa*zodgf(j1)

!         in cloudy part of layer
          zod2c(j1) = zudiff*pbsfc(j1,j3)*podsc(j1,j3)
          zod1c(j1) = zod2c(j1) + zudiff*podac(j1,j3)
          zod1c(j1) = zod1c(j1) + zangfa*zodgc(j1)
        END DO
 
!     Layer coefficients
 
        DO j1=ki1sc,ki1ec
!
!         in cloud-free part of layer
!
          zeps = SQRT(zod1f(j1)*zod1f(j1)-zod2f(j1)*zod2f(j1))
          ztau = EXP( -MIN( zeps, zargli) )
          zrho = zod2f(j1)/(zod1f(j1)+zeps)
          pa1f(j1)=ztau*(1._wp-(zrho**2))*(1._wp/(1._wp-(zrho**2)*(ztau**2)))
          pa2f(j1)=zrho*(1._wp-(ztau**2))*(1._wp/(1._wp-(zrho**2)*(ztau**2)))
          pa3f(j1)=(1._wp-pa1f(j1)+pa2f(j1))  &
                              /(zod1f(j1)+zod2f(j1))
! 
!         in cloudy part of layer
! 
          zeps = SQRT(zod1c(j1)*zod1c(j1)-zod2c(j1)*zod2c(j1))
          ztau = EXP( -MIN( zeps, zargli) )
          zrho = zod2c(j1)/(zod1c(j1)+zeps)
          pa1c(j1)=ztau*(1._wp-(zrho**2))*(1._wp/(1._wp-(zrho**2)*(ztau**2)))
          pa2c(j1)=zrho*(1._wp-(ztau**2))*(1._wp/(1._wp-(zrho**2)*(ztau**2)))
          pa3c(j1)=(1._wp-pa1c(j1)+pa2c(j1))  &
                              /(zod1c(j1)+zod2c(j1))
        END DO

      IF (ldebug) THEN
        PRINT *,'**** nach securit auf optical depth '
        PRINT *,'**** zodgf(CO2+O3+H2O) =',zodgf(j1b)
        PRINT *,'**** zodgc(CO2+O3+H2O) =',zodgc(j1b)
        PRINT *,'**** cloud-free zod1f (j1b)=',zod1f (j1b)
        PRINT *,'**** cloud-free zod2f (j1b)=',zod2f (j1b)
        PRINT *,'**** cloud-free pa1f (j1b)=',pa1f (j1b)
        PRINT *,'**** cloud-free pa2f (j1b)=',pa2f (j1b)
        PRINT *,'**** cloud-free pa3f (j1b)=',pa3f (j1b)
      END IF
 
      RETURN
    END SUBROUTINE coe_th

    SUBROUTINE coe_so (                                              &
       pduh2oc,pduh2of,pduco2 ,pduo3  ,palogp ,palogt ,                &
       podsc  ,podsf  ,podac  ,podaf  ,pbsfc  ,pbsff  ,pusfc ,pusff ,  &
       psmu0  ,pqsmu0 ,                                                &
       ki3    ,kspec  ,kh2o   ,kco2   ,ko3    ,                        &
       kig1s ,kig1e, ki1sc  ,ki1ec  ,ki3sc  ,ki3ec,                    &
       ldebug ,                                                        &
       pa1c   ,pa1f   ,pa2c   ,pa2f   ,pa3c   ,pa3f ,                  &
       pa4c   ,pa4f   ,pa5c   ,pa5f )
!
!**** *coe_so* - calculates the optical effects of atmospheric
!                layers on solar radiation based on basic
!                optical properties of non-gaseous constituents
!                and gaseous absorption coefficients selected
!                through the corresponding control variables
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
!  Introduce KIND-notation for all reals with real2kind.pl,
!  zargli, ztsec given as constants
! 1.17       2002/05/03 A. Mueller
!  Code cleanup
! 1.22       2003/02/12 A. Mueller
!  Avoid one EXP calculation
! V2_20        2009/04/06 H. Frank
!  Restructuring for NEC SX. Replace ztsec by exp(-zargli)
! V2_24        2010/04/28 Michael Gertz
!  Adaptions to SVN      
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!=======================================================================
!     Purpose
!     -----
!
!     This routine computes layer effects (transmissivity, reflectivity)
!     for diffuse and direct solar radiation both for the cloudy and the
!     cloud-free part of a model layer.
!     The calculation is based on the implicit delt-two-stream equations
!     (cf. Ritter and Geleyn, 1992) and uses basic optical properties
!     (i.e. absorption and scattering optical depth, backscattered and
!     upscattered fraction for non-gaseous atmospheric constituents and
!     gaseous absorption properties) as input.

!     Method
!     ------
!
!     - addition of individual gaseous absorption effects to the optical
!       properties of the non-gaseous constituents
!       (optical depth multiplied by alpha1 to alpha4)
!
!     - determination of layer effects (cf. Zdunkowski et al., 1982, 1986
!       and Ritter and Geleyn, 1992)
!     - the resonance case for those effects related to the direct solar
!       radiation is avoided by a small displacement of the local inverse
!       of the cosine of the zenith angle (if necessary)
!==============================================================================

!     Input data

!     opticall relevant gas quantities (Pa)

      INTEGER ki3         ! vertical layer considered
      INTEGER kspec       ! spectral interval considered
      INTEGER kh2o        ! table index for h2o absorption properties
      INTEGER kco2        ! table index for co2 absorption properties
      INTEGER ko3         ! table index for o3  absorption properties

      !
      INTEGER kig1s ,kig1e
      INTEGER ki1sc          ! start index for first  array dimension
      INTEGER ki1ec          ! end   index for first  array dimension
      INTEGER ki3sc          ! start index for third  array dimension
      INTEGER ki3ec          ! end   index for third  array dimension

      REAL(wp) :: pduh2oc(kig1s:kig1e,ki3sc:ki3ec) ! h2o inside cloud
      REAL(wp) :: pduh2of(kig1s:kig1e,ki3sc:ki3ec) ! h2o out of cloud
      REAL(wp) :: pduco2 (kig1s:kig1e,ki3sc:ki3ec) ! co2 content
      REAL(wp) :: pduo3  (kig1s:kig1e,ki3sc:ki3ec) ! o3  content

!     Logarithm of layer mean temperature and pressure

      REAL(wp) :: palogt (kig1s:kig1e,ki3sc:ki3ec) ! ln T
      REAL(wp) :: palogp (kig1s:kig1e,ki3sc:ki3ec) ! ln p

!     Optical properties of non-gaseous constituents (..c=cloudy; ..f=free)

      REAL(wp) :: podsc  (kig1s:kig1e,ki3sc:ki3ec) !
      REAL(wp) :: podsf  (kig1s:kig1e,ki3sc:ki3ec) !
      REAL(wp) :: podac  (kig1s:kig1e,ki3sc:ki3ec) !
      REAL(wp) :: podaf  (kig1s:kig1e,ki3sc:ki3ec) !
      REAL(wp) :: pbsfc  (kig1s:kig1e,ki3sc:ki3ec) !
      REAL(wp) :: pbsff  (kig1s:kig1e,ki3sc:ki3ec) !
      REAL(wp) :: pusfc  (kig1s:kig1e,ki3sc:ki3ec) !
      REAL(wp) :: pusff  (kig1s:kig1e,ki3sc:ki3ec) !

      REAL(wp) :: psmu0  (kig1s:kig1e) ! cosine of zenith angle
      REAL(wp) :: pqsmu0 (kig1s:kig1e) ! inverse of cosine ...
!
      LOGICAL ldebug         ! debug control switch

!     Output data

      REAL(wp) :: pa1c  (kig1s:kig1e) ! direct radiation transmis-
      REAL(wp) :: pa1f  (kig1s:kig1e) ! sivity cloudy/cloud-free
      REAL(wp) :: pa2c  (kig1s:kig1e) ! direct radition downward
      REAL(wp) :: pa2f  (kig1s:kig1e) ! scattering cloudy/cloud-free
      REAL(wp) :: pa3c  (kig1s:kig1e) ! direct radiation upward
      REAL(wp) :: pa3f  (kig1s:kig1e) ! scattering cloudy/cloud-free
      REAL(wp) :: pa4c  (kig1s:kig1e) ! diffuse flux transmissivity
      REAL(wp) :: pa4f  (kig1s:kig1e) ! cloudy/cloud-free
      REAL(wp) :: pa5c  (kig1s:kig1e) ! diffuse flux reflectivity
      REAL(wp) :: pa5f  (kig1s:kig1e) ! cloudy/cloud-free

!     Local arrays and variables

      REAL(wp) :: zodgf  (kig1s:kig1e)
      REAL(wp) :: zodgc  (kig1s:kig1e)
      REAL(wp) :: zod1f  (kig1s:kig1e)
      REAL(wp) :: zod1c  (kig1s:kig1e)
      REAL(wp) :: zod2f  (kig1s:kig1e)
      REAL(wp) :: zod2c  (kig1s:kig1e)
      REAL(wp) :: zod3f  (kig1s:kig1e)
      REAL(wp) :: zod3c  (kig1s:kig1e)
      REAL(wp) :: zod4f  (kig1s:kig1e)
      REAL(wp) :: zod4c  (kig1s:kig1e)
      REAL(wp) :: zod5f  (kig1s:kig1e)
      REAL(wp) :: zod5c  (kig1s:kig1e)
!
      REAL(wp) :: ztmp   !temporary varable for EXP

!     Security parameters (machine dependent) and utility variables

      REAL(wp), PARAMETER :: zargli = 70.0_wp   ! argument limit for EXP
      REAL(wp), PARAMETER :: zodmax = 1.E+06_wp ! maximum allowed optical depth
      REAL(wp), PARAMETER :: zepres = 100._wp*EPSILON(1._wp) ! for resonance case avoidance
      REAL(wp) :: zeps
      REAL(wp) :: ze,zm,zg1,zg2,ze1mwf,zmu0if

!     Diffusivity factors for gases and other constituents

      REAL(wp), PARAMETER :: zudiff = 2.0_wp         ! gases
      REAL(wp) ::            zangfa

      INTEGER :: j1,j3       ! loop indices over spatial dimensions
!
!      INCLUDE "gme_com_gas.h"
!
      INTEGER :: j1b        ! index of debug-point
      j1b = ki1sc

      zangfa = EXP( 0.5_wp)
!
!=======================================================================
!
      j3     = ki3

      IF (ldebug) THEN
        PRINT *,'**** coe_so ******************************'
        PRINT *,'**** debug point index : ',j1b
        PRINT *,'**** coe_so kspec=',kspec
        PRINT *,'**** coe_so j3   =',j3
        PRINT *,'**** coe_so kh2o =',kh2o
        PRINT *,'**** coe_so kco2 =',kco2
        PRINT *,'**** coe_so ko3  =',ko3
        PRINT *,'**** pduh2of(j1b,j3)=',pduh2of(j1b,j3)
        PRINT *,'**** pduh2oc(j1b,j3)=',pduh2oc(j1b,j3)
        PRINT *,'**** pduco2 (j1b,j3)=',pduco2 (j1b,j3)
        PRINT *,'**** pduo3  (j1b,j3)=',pduo3  (j1b,j3)
        PRINT *,'**** psmu0  (j1b)   =',psmu0  (j1b)
      END IF

!     Optical depth of gases

        DO j1=ki1sc,ki1ec
          zodgf(j1) = 0.0_wp     ! Initialisation
        END DO

      IF ( kco2 >0) THEN     ! Include CO2 contribution

          DO j1=ki1sc,ki1ec
            zodgf(j1)=zodgf(j1) + pduco2(j1,j3)*    &
                   ( cobi(kco2,kspec,2)*                     &
               EXP (coali(kco2,kspec,2)*palogp(j1,j3)     &
                   -cobti(kco2,kspec,2)*palogt(j1,j3)))
          END DO
      END IF                  ! CO2
      IF (ldebug) PRINT *,'**** zodgf(CO2)        =',zodgf(j1b)

      IF ( ko3 >0) THEN  ! Include O3 contribution
          DO j1=ki1sc,ki1ec
            zodgf(j1)=zodgf(j1) + pduo3 (j1,j3)*    &
                   ( cobi(ko3 ,kspec,3)*                     &
               EXP (coali(ko3 ,kspec,3)*palogp(j1,j3)     &
                   -cobti(ko3 ,kspec,3)*palogt(j1,j3)))
          END DO
      END IF
      IF (ldebug) PRINT *,'**** zodgf(CO2+O3)     =',zodgf(j1b)
!
!     Cloudy = cloud free for CO2 and O3 :
!
        DO j1=ki1sc,ki1ec
          zodgc(j1) = zodgf(j1)
        END DO

      IF ( kh2o >0) THEN  ! Include H2O contribution
          DO j1=ki1sc,ki1ec
            ztmp =   ( cobi(kh2o,kspec,1)*                      &
                 EXP (coali(kh2o,kspec,1)*palogp(j1,j3)      &
                     -cobti(kh2o,kspec,1)*palogt(j1,j3)))
             zodgf(j1)=zodgf(j1) + pduh2of(j1,j3)*ztmp
             zodgc(j1)=zodgc(j1) + pduh2oc(j1,j3)*ztmp
          END DO
      END IF
      IF (ldebug) THEN
        PRINT *,'**** zodgf(CO2+O3+H2O) =',zodgf(j1b)
        PRINT *,'**** zodgc(CO2+O3+H2O) =',zodgc(j1b)
      END IF

        DO j1=ki1sc,ki1ec
          zodgf (j1)=MIN(zodgf(j1),zodmax)
          zodgc (j1)=MIN(zodgc(j1),zodmax)
        END DO
      IF (ldebug) THEN
        PRINT *,'**** after security on optical depth '
        PRINT *,'**** zodgf(CO2+O3+H2O) =',zodgf(j1b)
        PRINT *,'**** zodgc(CO2+O3+H2O) =',zodgc(j1b)
      END IF
!
!     Pseudo-optical depth
!
        DO j1=ki1sc,ki1ec
!
!         in cloud-free part of layer
!
          zod2f(j1) = zudiff*pbsff(j1,j3)*podsf(j1,j3)
          zod1f(j1) = zod2f(j1) + zudiff*podaf(j1,j3) + zangfa*zodgf(j1)
          zod3f(j1) = pusff(j1,j3)*podsf(j1,j3)
          zod4f(j1) = podsf(j1,j3) - zod3f(j1)
          zod5f(j1) = podsf(j1,j3) + podaf(j1,j3) + zodgf(j1)
!
!         in cloudy part of layer
!
          zod2c(j1) = zudiff*pbsfc(j1,j3)*podsc(j1,j3)
          zod1c(j1) = zod2c(j1) + zudiff*podac(j1,j3) + zangfa*zodgc(j1)
          zod3c(j1) = pusfc(j1,j3)*podsc(j1,j3)
          zod4c(j1) = podsc(j1,j3) - zod3c(j1)
          zod5c(j1) = podsc(j1,j3) + podac(j1,j3) + zodgc(j1)
        END DO
      IF (ldebug) THEN
        PRINT *,'**** cloud-free zod1f (j1b)    =',zod1f (j1b)
        PRINT *,'**** cloud-free zod2f (j1b)    =',zod2f (j1b)
        PRINT *,'**** cloud-free zod3f (j1b)    =',zod3f (j1b)
        PRINT *,'**** cloud-free zod4f (j1b)    =',zod4f (j1b)
        PRINT *,'**** cloud-free zod5f (j1b)    =',zod5f (j1b)
        PRINT *,'**** cloudy     zod1c (j1b)    =',zod1c (j1b)
        PRINT *,'**** cloudy     zod2c (j1b)    =',zod2c (j1b)
        PRINT *,'**** cloudy     zod3c (j1b)    =',zod3c (j1b)
        PRINT *,'**** cloudy     zod4c (j1b)    =',zod4c (j1b)
        PRINT *,'**** cloudy     zod5c (j1b)    =',zod5c (j1b)
      END IF
!
!     Layer coefficients:
!
        DO j1=ki1sc,ki1ec
!
!     Layer coefficients in cloud-free part of layer
!
          zeps=SQRT(zod1f(j1)*zod1f(j1)-zod2f(j1)*zod2f(j1))
          ze = EXP( -MIN( zeps, zargli) )
          zm = zod2f(j1)/(zod1f(j1)+zeps)
          pa4f(j1)=ze*(1._wp-(zm**2))*(1._wp/(1._wp-(zm**2)*(ze**2)))
          pa5f(j1)=zm*(1._wp-(ze**2))*(1._wp/(1._wp-(zm**2)*(ze**2)))

          ze1mwf = zeps / zod5f(j1)
          zmu0if = ze1mwf + SIGN ( MAX(ABS(pqsmu0(j1)-ze1mwf),zepres),  &
                                  (pqsmu0(j1)-ze1mwf) )
          zod3f(j1)= zod3f(j1)*zmu0if
          zod4f(j1)= zod4f(j1)*zmu0if
          zod5f(j1)= zod5f(j1)*zmu0if
          pa1f(j1) = EXP( -MIN( zod5f(j1), zargli) )
          zg1 = ( zod3f(j1)*(zod5f(j1)-zod1f(j1))   &
                 -zod2f(j1)* zod4f(j1))                &
                /(zod5f(j1)* zod5f(j1) - zeps*zeps)
          zg2 =-( zod4f(j1)*(zod5f(j1)+zod1f(j1))   &
                 +zod2f(j1)* zod3f(j1))                &
                /(zod5f(j1)* zod5f(j1) - zeps*zeps)
          pa2f(j1) = zg2*(pa1f(j1)-pa4f(j1))        &
                        -zg1*pa5f(j1)*pa1f(j1)
          pa3f(j1) = zg1*(1._wp-pa4f(j1)*pa1f(j1))  &
                       -zg2*pa5f(j1)

!     Layer coefficients in cloudy part of layer

          zeps=SQRT(zod1c(j1)*zod1c(j1)-zod2c(j1)*zod2c(j1))
          ze = EXP( -MIN( zeps, zargli) )
          zm = zod2c(j1)/(zod1c(j1)+zeps)
          pa4c(j1)=ze*(1._wp-(zm**2))*(1._wp/(1._wp-(zm**2)*(ze**2)))
          pa5c(j1)=zm*(1._wp-(ze**2))*(1._wp/(1._wp-(zm**2)*(ze**2)))

          ze1mwf = zeps / zod5c(j1)
          zmu0if = ze1mwf + SIGN ( MAX(ABS(pqsmu0(j1)-ze1mwf),zepres)  &
                                  ,(pqsmu0(j1)-ze1mwf) )
          zod3c(j1)= zod3c(j1)*zmu0if
          zod4c(j1)= zod4c(j1)*zmu0if
          zod5c(j1)= zod5c(j1)*zmu0if
          pa1c(j1) = EXP( -MIN( zod5c(j1), zargli) )
          zg1 = ( zod3c(j1)*(zod5c(j1)-zod1c(j1))   &
                 -zod2c(j1)*zod4c(j1))                &
                /(zod5c(j1)*zod5c(j1) - zeps*zeps)
          zg2 =-( zod4c(j1)*(zod5c(j1)+zod1c(j1))   &
                 +zod2c(j1)*zod3c(j1))                &
                /(zod5c(j1)*zod5c(j1) - zeps*zeps)
          pa2c(j1) = zg2*(pa1c(j1)-pa4c(j1))     &
                        -zg1*pa5c(j1)*pa1c(j1)
          pa3c(j1) = zg1*(1._wp-pa4c(j1)*pa1c(j1))  &
                       -zg2*pa5c(j1)
        END DO

      IF (ldebug) THEN
        PRINT *,'**** cloud-free pa1f (j1b)    =',pa1f (j1b)
        PRINT *,'**** cloud-free pa2f (j1b)    =',pa2f (j1b)
        PRINT *,'**** cloud-free pa3f (j1b)    =',pa3f (j1b)
        PRINT *,'**** cloud-free pa4f (j1b)    =',pa4f (j1b)
        PRINT *,'**** cloud-free pa5f (j1b)    =',pa5f (j1b)
        PRINT *,'**** cloudy     pa1c (j1b)    =',pa1c (j1b)
        PRINT *,'**** cloudy     pa2c (j1b)    =',pa2c (j1b)
        PRINT *,'**** cloudy     pa3c (j1b)    =',pa3c (j1b)
        PRINT *,'**** cloudy     pa4c (j1b)    =',pa4c (j1b)
        PRINT *,'**** cloudy     pa5c (j1b)    =',pa5c (j1b)
      END IF

      RETURN
    END SUBROUTINE coe_so

    
END MODULE mo_radiation_rg

