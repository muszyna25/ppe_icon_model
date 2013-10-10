MODULE mo_ice_dyn_parms
!introduced by Ralph Timmermann, 17.9.2004
USE mo_kind,    ONLY: wp
implicit none
PUBLIC
REAL(wp), parameter :: cdwin = 2.25e-3_wp   ! drag coefficient atmosphere - ice FR284
REAL(wp), parameter :: cdwat = 5.00e-3_wp   ! drag coefficient ocean - ice      FR284
REAL(wp), parameter :: cdao  = 1.20e-3_wp   ! drag coefficient atmosphere - ocean !FR284

END MODULE mo_ice_dyn_parms
!==============================================================================
