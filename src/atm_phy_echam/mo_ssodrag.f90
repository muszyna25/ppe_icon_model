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

  USE mo_kind,      ONLY: wp
  USE mo_exception, ONLY: finish

  IMPLICIT NONE

  ! nombre de vrais traceurs
  INTEGER, PARAMETER :: nqmx=2
  INTEGER, PARAMETER :: nbtr=nqmx-2+1/(nqmx-1)

  INTEGER :: nktopg ! Security value for blocked flow level
  INTEGER :: ntop = 1   ! An estimate to qualify the upper levels of
  !                       the model where one wants to impose strees
  !                       profiles
  INTEGER ::  nstra ! no documentation
  !
  ! Parameters depending on model resolution
  !
  REAL(wp) :: gpicmea   ! (PEAK-mean) threshold for activation of scheme
  REAL(wp) :: gstd      ! Standard deviation threshold for activation of scheme
  REAL(wp) :: gkdrag    ! Gravity wave drag coefficient (G in (3),      LOTT 1999)
  REAL(wp) :: gkwake    ! Bluff-body drag coefficient for low level wake 
  !                    (Cd in (2), LOTT 1999)

  !      SET_UP THE "TUNABLE PARAMETERS" OF THE VARIOUS SSO SCHEMES

  REAL(wp), PARAMETER :: gfrcrit = 0.5_wp  ! Critical Non-dimensional mountain 
  !                                    Height (HNC in (1), LOTT 1999)
  REAL(wp), PARAMETER :: grcrit  = 0.25_wp ! Critical Richardson Number (Ric, End of
  !                                   first column p791 of LOTT 1999)
  REAL(wp), PARAMETER :: gklift  = 0.00_wp ! Mountain Lift coefficient
  !                                   (Cl in (4),     LOTT 1999)
  REAL(wp), PARAMETER :: grahilo = 1.00_wp ! Set-up the trapped waves fraction
  !                                  (Beta , End of first column,  LOTT 1999)
  REAL(wp), PARAMETER :: ghmax   = 10000.0_wp  ! Not used
  REAL(wp), PARAMETER :: gvcrit  = 0.1_wp     ! no documentation

  !       SET_UP  VALUES OF SECURITY PARAMETERS

  REAL(wp), PARAMETER :: gsigcr = 0.80_wp    ! Security value for blocked flow depth
  REAL(wp), PARAMETER :: gssec  = 0.0001_wp  ! Security min value for low-level B-V 
  !                                     frequency
  REAL(wp), PARAMETER :: gtsec  = 0.00001_wp ! Security min value for anisotropy 
  !                                     and GW stress.
  REAL(wp), PARAMETER :: gvsec  = 0.10_wp    ! Security min value for ulow

  PUBLIC :: sugwd
  PUBLIC :: nktopg, ntop, gstd, gkdrag,gkwake,gfrcrit,grcrit 
  PUBLIC :: gklift, gsigcr, gssec, gtsec, gvsec, grahilo,gpicmea

CONTAINS
  !======================================================================
  SUBROUTINE sugwd(klon,klev)

!  USE mo_control, ONLY: nvclev, vct, nn

  !  Scalar arguments with intent(In):
  INTEGER, INTENT (IN) :: klon, klev

!! ATTENTION
    INTEGER:: nvclev,  nn
    REAL(wp):: vct(klev)

  ! local scalar
  INTEGER :: jk
  REAL(wp):: zstra, zsigt, zpm1r, zpr

  !          SET THE VALUES OF THE PARAMETERS
  !
  IF (nn == 31) THEN
    gpicmea = 500.0_wp
    gstd    = 200.0_wp
    gkdrag  = 0.8_wp
    gkwake  = 0.8_wp
  ELSE IF(nn == 42) THEN
    gpicmea = 400.0_wp
    gstd    = 100.0_wp
    gkdrag  = 0.8_wp
    gkwake  = 0.8_wp
  ELSE IF(nn == 63) THEN
    gpicmea = 400.0_wp
    gstd    = 100.0_wp
    IF(klev == 95) THEN
      gkdrag  = 0.25_wp
      gkwake  = 0.25_wp
    ELSE
      gkdrag  = 0.5_wp
      gkwake  = 0.5_wp
    END IF
  ELSE IF(nn == 127) THEN
    gpicmea = 200.0_wp
    gstd    =  50.0_wp
    gkdrag  = 0.5_wp
    gkwake  = 0.5_wp
  ELSE IF(nn == 255) THEN
    gpicmea = 100.0_wp
    gstd    =  10.0_wp
    gkdrag  = 0.8_wp
    gkwake  = 0.8_wp
  ELSE
    CALL finish ('mo_ssodrag', 'Truncation not supported.')
  ENDIF

  ! PRINT *,' DANS SUGWD NLEV=',klev

  !
  zpr=80000.0_wp
  zstra=0.001_wp
  zsigt=0.94_wp
  !old  ZSIGT=0.85_wp
  !
  DO 110 jk=klev,1,-1
    zpm1r = 0.5_wp*(vct(jk)+vct(jk+1)+zpr*(vct(nvclev+jk)+vct(nvclev+jk+1)))
    zpm1r = zpm1r/zpr

    IF(zpm1r >= zsigt)THEN
      nktopg=jk
    ENDIF
    IF(zpm1r >= zstra)THEN
      nstra=jk-1
    ENDIF
110 END DO

!    PRINT *,' DANS SUGWD nktopg=', nktopg
!    PRINT *,' DANS SUGWD nstra=', nstra
!    PRINT *,' DANS SUGWD ntop=', ntop
    !


!    WRITE(unit=6,fmt='('' *** SSO essential constants ***'')')
!    WRITE(unit=6,fmt='('' *** SPECIFIED IN SUGWD ***'')')
!    WRITE(unit=6,fmt='('' Gravity wave ct '',E13.7,'' '')')gkdrag
!    WRITE(unit=6,fmt='('' Trapped/total wave dag '',E13.7,'' '')')    &
!         grahilo
!    WRITE(unit=6,fmt='('' Critical Richardson   = '',E13.7,'' '')')   &
!         grcrit
!    WRITE(unit=6,fmt='('' Critical Froude'',e13.7)') gfrcrit
!    WRITE(unit=6,fmt='('' Low level Wake bluff cte'',e13.7)') gkwake
!    WRITE(unit=6,fmt='('' Low level lift  cte'',e13.7)') gklift

  END SUBROUTINE sugwd
  !======================================================================
END MODULE mo_ssodrag
