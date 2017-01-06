!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
MODULE mo_ssodrag

  ! Description:
  !
  ! Set up parameters for gravity wave drag calculations
  !
  ! Authors:
  !           Martin Miller, ECMWF, Jan 1990
  !           Francois Lott, LMD,   Jul 1999  
  !           Elisa Manzini, MPI,   Aug 2000
  !
  ! References: 
  !     Lott, 1999: Alleviation of stationary biases in a GCM through...
  !                 Monthly Weather Review, 127, pp 788-801.

  USE mo_kind                 ,ONLY: wp
  USE mo_exception            ,ONLY: message, print_value

  USE mo_run_config           ,ONLY: nlevp1
  USE mo_vertical_coord_table ,ONLY: vct

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: sugwd  ,                          &
    &       nktopg , ntop  ,                  &
    &       gpicmea, gstd  , gkdrag, gkwake , &
    &       gfrcrit, grcrit, gklift, grahilo, & 
    &       gsigcr , gssec , gtsec , gvsec

  CHARACTER(len=*), PARAMETER :: thismodule = 'mo_ssodrag'

  INTEGER  :: nktopg    ! Security value for blocked flow level
  INTEGER  :: ntop = 1  ! An estimate to qualify the upper levels of
                        ! the model where one wants to impose strees
                        ! profiles
  !
  ! Parameters depending on model resolution
  !
  REAL(wp) :: gpicmea   ! (PEAK-mean) threshold for activation of scheme
  REAL(wp) :: gstd      ! Standard deviation threshold for activation of scheme
  REAL(wp) :: gkdrag    ! Gravity wave drag coefficient                  (G  in (3), LOTT 1999)
  REAL(wp) :: gkwake    ! Bluff-body drag coefficient for low level wake (Cd in (2), LOTT 1999)

  !      SET_UP THE "TUNABLE PARAMETERS" OF THE VARIOUS SSO SCHEMES

  REAL(wp), PARAMETER :: gfrcrit = 0.5_wp      ! Critical Non-dimensional mountain Height
  !                                              (HNC in (1), LOTT 1999)
  REAL(wp), PARAMETER :: grcrit  = 0.25_wp     ! Critical Richardson Number 
  !                                              (Ric, end of first column p791, LOTT 1999)
  REAL(wp), PARAMETER :: gklift  = 0.00_wp     ! Mountain Lift coefficient
  !                                              (Cl in (4), LOTT 1999)
  REAL(wp), PARAMETER :: grahilo = 1.00_wp     ! Set-up the trapped waves fraction
  !                                              (Beta , end of first column, LOTT 1999)
  REAL(wp), PARAMETER :: ghmax   = 10000.0_wp  ! Not used
  REAL(wp), PARAMETER :: gvcrit  = 0.1_wp      ! no documentation

  !       SET_UP  VALUES OF SECURITY PARAMETERS

  REAL(wp), PARAMETER :: gsigcr = 0.80_wp      ! Security value for blocked flow depth
  REAL(wp), PARAMETER :: gssec  = 0.0001_wp    ! Security min value for low-level B-V frequency
  REAL(wp), PARAMETER :: gtsec  = 0.00001_wp   ! Security min value for anisotropy and GW stress.
  REAL(wp), PARAMETER :: gvsec  = 0.10_wp      ! Security min value for ulow

CONTAINS
  !======================================================================
  SUBROUTINE sugwd(klev)

  CHARACTER(len=*), PARAMETER :: routine = 'ssodrag'

  !  Scalar arguments with intent(In):
  INTEGER, INTENT (IN) :: klev

  ! local scalar
  INTEGER  :: jk
  REAL(wp) :: zsigt, zpm1r, zpr

  !          SET THE VALUES OF THE PARAMETERS
  !

  ! Define the mask for the SSO parameterization:
    gpicmea = 40.0_wp ! only where  (peak - mean height) is typically > 1st layer depth
    gstd    = 10.0_wp ! only where SSO slope, asymmetry and orientation are defined by EXTPAR

  ! Define the tuning parameters for SSO drag. These values depend on:
  ! (1) the resolution of the topography DATA used to compute the SSO parameters, and
  ! (2) the model resolution.
    gkdrag  = 0.05_wp ! gravity wave drag
    gkwake  = 0.00_wp ! low level blocking

  ! Compute globally valid layer index nktopg from the vertical coordinate tables that
  ! describe the hybrid sigma coordinate as follows:
  ! - Pressure sigma coord.: ph(jk) = a(jk) + p_sfc*b(jk) [Pa]
  ! - Height   sigma coord.: zh(jk) = a(jk) + z_sfc*b(jk) [m]
  !
  ! a(:) and b(:) are obtained from vct(:) as follows:
  ! - a(1:nlev+1) = vct(       1:       nlev+1) and 
  ! - b(1:nlev+1) = vct(nlevp1+1:nlevp1+nlev+1) with nlevp1=nlev+1
  !
  ! The vertical arrays are ordered from the top of the model (tom) to the surface (sfc):
  ! - ph(1) = p_tom, ph(nlev+1) = p_sfc
  ! - zh(1) = z_tom, zh(nlev+1) = z_sfc
  !
  ! The a(:) coefficients at the top of the model, where the sigma portion b(:) is zero,
  ! can be used to distinguish the pressure and height sigma grid:
  ! - pressure sigma grid : a(1) < a(2)
  ! - height   sigma grid : a(1) > a(2)

  IF (vct(1) < vct(2)) THEN
    !
    ! pressure sigma grid
    !
    zpr   = 80000.0_wp ! Pa   , surface pressure
    zsigt = 0.94_wp    ! Pa/Pa, sigma for blocked flow depth (0.94 * 800 hPa = 750 hPa)
    !
    DO 110 jk=klev,1,-1
      !
      ! full level pressure pf(jk) = (ph(jk)+ph(jk+1)/2 for p_sfc=zpr
      zpm1r = 0.5_wp*(vct(jk)+vct(jk+1)+zpr*(vct(nlevp1+jk)+vct(nlevp1+jk+1)))
      !
      ! full level sigma(jk) = pf(jk)/p_sfc for p_sfc=zpr
      zpm1r = zpm1r/zpr
      !
      ! Find highest full level with sigma(jk) >= zsigt
      IF (zpm1r >= zsigt) THEN
        nktopg=jk
      END IF
      !
110 END DO
    !
  ELSE
    !
    ! height sigma grid
    !
    zpr   =  1950._wp ! m (800 hPa in International Standard Atmosphere)
    zsigt =  2460._wp ! m (750 hPa in International Standard Atmosphere)
    !
    DO 120 jk=klev,1,-1
      !
      ! full level height zf(jk) = (zh(jk)+zh(jk+1)/2 for z_sfc=zpr
      zpm1r = 0.5_wp*(vct(jk)+vct(jk+1)+zpr*(vct(nlevp1+jk)+vct(nlevp1+jk+1)))
      !
      ! Find highest full level with zf(jk) <= zsigt 
      IF (zpm1r <= zsigt) THEN
        nktopg=jk
      END IF
      !
120 END DO
    !
  END IF

  CALL message(thismodule//':'//routine,'*** SSO essential constants ***')
  CALL print_value('Gravity wave drag coeff.    =',gkdrag )
  CALL print_value('Trapped/total wave drag     =',grahilo)
  CALL print_value('Critical Richardson number  =',grcrit )
  CALL print_value('Critical Froude number      =',gfrcrit)
  CALL print_value('Low level wake bluff coeff. =',gkwake )
  CALL print_value('Low level lift  coeff.      =',gklift )
  CALL message(thismodule//':'//routine,'-------------------------------')

  END SUBROUTINE sugwd
  !======================================================================
END MODULE mo_ssodrag
