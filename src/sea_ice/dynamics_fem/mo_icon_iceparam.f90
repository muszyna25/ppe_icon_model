
!
MODULE mo_ice_iceparam
  !
  ! Ice specific parameters
  !
  USE mo_kind,    ONLY: wp
  USE mo_physical_constants,    ONLY: Cd_io
  IMPLICIT NONE
  PUBLIC
  SAVE

  REAL(wp), ALLOCATABLE, DIMENSION(:)    :: coriolis_nod2D
 ! ice model parameters:
  ! RHEOLOGY
  REAL(wp), PARAMETER  :: Pstar = 27500._wp         !15000 !was 20000.   ![N/m^2]
  REAL(wp), PARAMETER  :: ellipse =2._wp            !
  REAL(wp), PARAMETER  :: c_pressure =20.0_wp       !
  REAL(wp)             :: delta_min=2.0e-9_wp       ! [s^(-1)]
  REAL(wp)             :: Clim_evp=615._wp          ! kg/m^2
  REAL(wp)             :: zeta_min=4.0e+8_wp        ! kg/s
  INTEGER                  :: evp_rheol_steps=120       ! EVP rheology
                                                        ! cybcycling steps
  REAL(wp)             :: ice_gamma_fct=0.25_wp     ! smoothing parameter
                                                        ! in ice fct advection
  REAL(wp)             :: ice_diff=100.0_wp         ! diffusion to stabilize
                                                        ! ice advection
  LOGICAL                  :: ice_VP_rheology=.true.    ! VP=.true., EVP=.false.
  REAL(wp)             :: ice_VP_soltol=1.0e-6_wp
  REAL(wp)             :: Tevp_inv
  REAL(wp)             :: theta_io                  ! rotation angle
                                                        ! (ice-ocean), available
              ! in EVP
  INTEGER                  :: ice_advection         ! type of ice advection
                                                        ! scheme
  INTEGER, PARAMETER       :: ice_FCT=1
  INTEGER, PARAMETER       :: ice_BE=0
  REAL(wp), PARAMETER  :: C_d_io= Cd_io !0.0055_wp!sd  was 0.003    ! ice-water drag coefficient
  REAL(wp), PARAMETER  :: rho_ice=900._wp
  REAL(wp), PARAMETER  :: rho_snow=300._wp
END MODULE mo_ice_iceparam

!=============================================================================
