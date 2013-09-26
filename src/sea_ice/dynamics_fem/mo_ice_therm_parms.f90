MODULE mo_ice_therm_parms
USE mo_kind,    ONLY: wp
implicit none
PUBLIC
REAL(wp), parameter :: rhoair=  1.3_wp    ! Air density
REAL(wp), parameter :: rhowat= 1025._wp   ! Water density
REAL(wp), parameter :: rhoice=  910._wp   ! Ice density
REAL(wp), parameter :: rhosno=  290._wp   ! Snow density

REAL(wp), parameter :: sice = 5.0_wp      ! Ice salinity 3.2--5.0 ppt.

integer, parameter :: iclasses=7     ! Number of ice thickness gradations for ice growth calcs.
REAL(wp), parameter :: h0=1.0_wp     ! Lead, closing parameter 0.5 [m] standard
REAL(wp), parameter :: Saterm=0.5_wp   ! Sa - term parameter 0.5 [m] standard
REAL(wp), parameter :: hmin= 0.05_wp   ! Cut-off ice thickness
REAL(wp), parameter :: Armin=0.15_wp   ! Minimum ice concentration

REAL(wp), parameter :: tmelt=273.16_wp ! 0 deg C expressed in K
REAL(wp), parameter :: cc=4.20E6_wp    ! Volumetr. heat cap. of water [J/m**3/K](cc = rhoice*cp)
REAL(wp), parameter :: cl=3.02E8_wp    ! Volumetric latent heat of sea ice [J/m**3]
REAL(wp), parameter :: q0=1._wp/cl     ! 1/volumetric heat of fusion
REAL(wp), parameter :: cpair=1004._wp  ! Specific heat of air [J/(kg * K)]
REAL(wp), parameter :: clhw=2.500E6_wp ! Specific latent heat [J/kg]: water -> water vapor;
REAL(wp), parameter :: clhi=2.834E6_wp !                              sea ice -> water vapor
REAL(wp), parameter :: cdsens=1.75E-3_wp  ! Bulk sensible heat transfer coefficient, SIOM standard
REAL(wp), parameter :: cdlat =1.75E-3_wp  ! Bulk latent heat transfer coefficients,  SIOM standard

!RT:
REAL(wp), parameter :: gamma_t=1.e-4_wp   ! heat transfer coefficient ocean -> ice
!RT-

REAL(wp), parameter :: qsw=17.2694_wp  ! Constants for latent heat fluxes
REAL(wp), parameter :: tqw=237.3_wp    ! over water
REAL(wp), parameter :: qsi=21.8746_wp  ! and ice
REAL(wp), parameter :: tqi=265.5_wp

REAL(wp), parameter :: d1=rhoair*cpair*cdsens  ! coefficients for bulk formulas for heat fluxes
REAL(wp), parameter :: d2w=rhoair*clhw*cdlat
REAL(wp), parameter :: d2i=rhoair*clhi*cdlat

REAL(wp), parameter :: Cd_thrm_i_o=1.0e-3_wp  ! Ocean-Ice thermoconductivity coefficient
!rt REAL(wp), parameter :: H_ML=30.            ! Ocean Mixed layer depth

REAL(wp), parameter :: emiss=0.97_wp   !rt was 0.99     ! Emissivity of Snow/Ice

REAL(wp), parameter :: boltzmann=5.67E-8_wp   ! S. Boltzmann const.*longw. emissivity
REAL(wp), parameter :: d3=boltzmann*emiss  ! SIOM standard (MH)
REAL(wp), parameter :: con   = 2.1656_wp ! Thermal conductivities: ice;
REAL(wp), parameter :: consn = 0.31_wp   !                         snow

!REAL(wp), parameter :: albsn=  0.81         ! Albedo: frozen snow
!REAL(wp), parameter :: albsnm=    0.77         !         melting snow
!REAL(wp), parameter :: albi= 0.70         !         frozen ice
!REAL(wp), parameter :: albm= 0.68         !         melting ice
!REAL(wp), parameter :: albw= 0.10         !         open water

!++RT BRIOS:
REAL(wp), parameter :: albsn= 0.85_wp      ! Albedo: frozen snow
REAL(wp), parameter :: albsnm=0.75_wp      !         melting snow
REAL(wp), parameter :: albi=  0.75_wp      !         frozen ice
REAL(wp), parameter :: albm=  0.66_wp      !         melting ice
REAL(wp), parameter :: albw=  0.10_wp      !         open water
!-RT--

END MODULE mo_ice_therm_parms

!=======================================================================
