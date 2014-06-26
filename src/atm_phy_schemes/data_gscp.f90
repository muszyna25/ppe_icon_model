!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!+ Data module for variables of the grid scale parameterization
!------------------------------------------------------------------------------

MODULE data_gscp

!------------------------------------------------------------------------------
!
! Description:
!  This module contains variables that are used in the grid scale 
!  parameterizations (Microphysics). 
!
! Current Code Owner: DWD, Axel Seifert
!  phone:  +49  69  8062 2729
!  fax:    +49  69  8062 3721
!  email:  ulrich.schaettler@dwd.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! 3.22       2007/01/24 Axel Seifert
!  Initial Release
! V4_5         2008/09/10 Ulrich Schaettler
!  Added variables mu_rain and cloud_num, which are now Namelist variables
! V4_13        2010/05/11 Michael Gertz
!  Adaptions to SVN
! V4_14        2010/06/14 Axel Seifert
!  Introduced v0snow as global variable
! V4_20        2011/08/31 Axel Seifert
!  Moved some global variables from src_gscp to data_gscp
! V4_21        2011/12/06 Axel Seifert
!  Additional variable rain_n0_factor
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!==============================================================================
!
! Declarations:
!
! Modules used:
#ifdef __COSMO__
USE data_parameters, ONLY :   &
    ireals,    & ! KIND-type parameter for real variables
    iintegers    ! KIND-type parameter for standard integer variables
#endif

#ifdef __ICON__
USE mo_kind,               ONLY: ireals=>wp     , &
                                 iintegers=>i4
#endif


!==============================================================================

IMPLICIT NONE


PUBLIC

!==============================================================================


! Variables for hydci_pp
! ----------------------

  REAL (KIND=ireals) ::           &
    ccsrim,    & !
    ccsagg,    & !
    ccsdep,    & !
    ccsvel,    & !
    ccsvxp,    & !
    ccslam,    & !
    ccslxp,    & !
    ccsaxp,    & !
    ccsdxp,    & !
    ccshi1,    & !
    ccdvtp,    & !
    ccidep,    & !
    ccswxp,    & !
    zconst,    & !
    zcev,      & !
    zbev,      & !
    zcevxp,    & !
    zbevxp,    & !
    zvzxp,     & !
    zvz0r,     & !
  !CK>
    vtxexp,    & !!
    kc_c1,     & !!
    kc_c2        !!  

! Variables for hydci_pp, hydci_pp_ice and hydci_pp_gr
! --------------------------------------------------------

  REAL (KIND=ireals) ::              &
    v0snow         = 20.0_ireals,    & ! factor in the terminal velocity for snow
    rain_n0_factor = 1.0_ireals,     & ! COSMO_EU default
    mu_rain,                         & ! 
    mu_snow
  
!CK>
! constants needed for cloud ice sedimentation
! Small hex plates from Mitchell (1996)
REAL(KIND=ireals), PARAMETER :: &  
  kc_alpha =  0.5870086, & !..alf, CGS is 0.00739
  kc_beta  =  2.45,      & !..exponent  in mass-size relation
  kc_gamma =  0.120285,  & !..gam, CGS is 0.24
  kc_sigma =  1.85    ,  & !..exponent  in area-size relation
  do_i     =  5.83    ,  & ! coefficients for drag correction
  co_i     =  0.6          ! coefficients for turbulence correction
!CK<

#ifdef __COSMO__
  REAL (KIND=ireals) ::              &
    cloud_num = 5.00e+08_ireals        ! cloud droplet number concentration
#else
    REAL (KIND=ireals) ::            &
    cloud_num = 200.00e+06_ireals      ! cloud droplet number concentration
#endif

    
#ifdef __COSMO__
    PRINT*,"define mu_rain and mu_snow in data_gscp.f90. Stop."
    STOP
#endif

! Parameters for Segal-Khain parameterization (aerosol-microphysics coupling)
REAL(KIND=ireals), PARAMETER :: r2_fix    = 0.03_ireals , & ! Parameters for simplified lookup table computation; 
                                lsigs_fix = 0.3_ireals      ! relevant for r2_lsigs_are_fixed = .TRUE.

LOGICAL,           PARAMETER :: r2_lsigs_are_fixed = .TRUE. , & !
                                lincloud           = .FALSE.    ! ignore in-cloud nucleation


!==============================================================================

END MODULE data_gscp
