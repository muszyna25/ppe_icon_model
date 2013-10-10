! =====================
!  Standing alone sea ice
!  Version 2, based on version 1, with several new features added
!  Most important are true VP solver and FCT advection
!  Questions to S. Danilov (dynamics) and Q. Wang and R. Timmermann
! (thermodynamics). Many updates and corrections
!  to version 1 are due to R. Timmermann
! ======================

! This file contails all modules, subroutine feim_step (ice model
! step) and main program. The modules here are derived from
! FESOM, the coupled version relies on full FESOM modules.


MODULE mo_ice_atm_forcing
USE mo_ice_data_types
USE mo_kind,    ONLY: wp
implicit none
PUBLIC
save

  REAL(wp), ALLOCATABLE, DIMENSION(:)         :: u_wind, v_wind
  REAL(wp), ALLOCATABLE, DIMENSION(:)         :: shortwave, longwave
  REAL(wp), ALLOCATABLE, DIMENSION(:)         :: Tair, Tdew
  REAL(wp), ALLOCATABLE, DIMENSION(:)         :: Pair, precipitation
  REAL(wp), ALLOCATABLE, DIMENSION(:)         :: evaporation
  REAL(wp), ALLOCATABLE, DIMENSION(:)         :: shum, cloudiness  !RTnew


END MODULE mo_ice_atm_forcing

!=====================================================================
