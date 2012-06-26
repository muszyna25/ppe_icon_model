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
    zvz0r

! Variables for hydci_pp, hydci_pp_ice and hydci_pp_gr
! --------------------------------------------------------

  REAL (KIND=ireals) ::              &
    v0snow         = 20.0_ireals,    & ! factor in the terminal velocity for snow
    mu_rain        = 0.0_ireals,     & ! COSMO_EU default
    rain_n0_factor = 1.0_ireals        ! COSMO_EU default


#ifdef __COSMO__
  REAL (KIND=ireals) ::              &
    cloud_num = 5.00e+08_ireals        ! cloud droplet number concentration
#else
    REAL (KIND=ireals) ::            &
    cloud_num = 200.00e+06_ireals      ! cloud droplet number concentration
#endif




!==============================================================================

END MODULE data_gscp
