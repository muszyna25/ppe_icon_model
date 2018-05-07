MODULE mo_psrad_setup

  USE mo_psrad_general, ONLY: wp

  PUBLIC :: psrad_basic_setup

  INTEGER, PARAMETER :: n_ref = 59, tropopause_ref = 13, jp_mult = 5
  REAL(wp), DIMENSION(n_ref) :: log_pref
  REAL(wp) :: jp_coeff

  REAL(wp), PARAMETER :: stpfac_orig  = 296._wp/1013._wp
  REAL(wp) :: stpfac, pressure_scale, droplet_scale



  ! Reference pressures are chosen such that the ln of the first pressure
  ! has only a few non-zero digits (i.e. ln(PREF(1)) = 6.96000) and
  ! each subsequent ln(pressure) differs from the previous one by 0.2.
  REAL(WP), PARAMETER :: jp_fudge = 0.04_wp, &
    dlog_pref = 0.2_wp, &
    log_pref_tropopause_orig = 4.56_wp, &
    ! These are the temperatures associated with the respective pressures
    tref(n_ref) = (/ &
      2.9420e+02_wp, 2.8799e+02_wp, 2.7894e+02_wp, 2.6925e+02_wp, &
      2.5983e+02_wp, 2.5017e+02_wp, 2.4077e+02_wp, 2.3179e+02_wp, &
      2.2306e+02_wp, 2.1578e+02_wp, 2.1570e+02_wp, 2.1570e+02_wp, &
      2.1570e+02_wp, 2.1706e+02_wp, 2.1858e+02_wp, 2.2018e+02_wp, &
      2.2174e+02_wp, 2.2328e+02_wp, 2.2479e+02_wp, 2.2655e+02_wp, &
      2.2834e+02_wp, 2.3113e+02_wp, 2.3401e+02_wp, 2.3703e+02_wp, &
      2.4022e+02_wp, 2.4371e+02_wp, 2.4726e+02_wp, 2.5085e+02_wp, &
      2.5457e+02_wp, 2.5832e+02_wp, 2.6216e+02_wp, 2.6606e+02_wp, &
      2.6999e+02_wp, 2.7340e+02_wp, 2.7536e+02_wp, 2.7568e+02_wp, &
      2.7372e+02_wp, 2.7163e+02_wp, 2.6955e+02_wp, 2.6593e+02_wp, &
      2.6211e+02_wp, 2.5828e+02_wp, 2.5360e+02_wp, 2.4854e+02_wp, &
      2.4348e+02_wp, 2.3809e+02_wp, 2.3206e+02_wp, 2.2603e+02_wp, &
      2.2000e+02_wp, 2.1435e+02_wp, 2.0887e+02_wp, 2.0340e+02_wp, &
      1.9792e+02_wp, 1.9290e+02_wp, 1.8809e+02_wp, 1.8329e+02_wp, &
      1.7849e+02_wp, 1.7394e+02_wp, 1.7212e+02_wp /) 

  PUBLIC :: log_pref, dlog_pref, tref, n_ref, &
    tropopause_ref, jp_coeff, jp_mult, jp_fudge, &
    stpfac, stpfac_orig, pressure_scale, droplet_scale

CONTAINS

  SUBROUTINE psrad_basic_setup(upwards_, klev, &
    pressure_scale_, droplet_scale_, &
    zinhoml1_, zinhoml2_, zinhoml3_, zinhomi_)
    USE mo_psrad_cloud_optics, ONLY: setup_cloud_optics
    USE mo_psrad_flat_data, ONLY: setup_flat_data
    USE mo_psrad_general, ONLY: finish_cb, default_finish, &
      dummy4, upwards, jTOA, jSFC, jINC, jABOVE, jBELOW

    IMPLICIT NONE

    LOGICAL, INTENT(IN) :: upwards_
    INTEGER, INTENT(IN) :: klev
    REAL(wp), INTENT(IN) :: pressure_scale_, droplet_scale_, &
      zinhoml1_, zinhoml2_, zinhoml3_, zinhomi_
    INTEGER :: i
    REAL(wp) :: log_factor

    ALLOCATE(dummy4(1,1,1,1))
    finish_cb => default_finish

    upwards = upwards_
    pressure_scale = pressure_scale_
    droplet_scale = droplet_scale_
    IF (upwards) THEN
      jTOA = klev
      jSFC = 1
      jINC = 1
      jABOVE = 1
      jBELOW = 0
    ELSE
      jTOA = 1
      jSFC = klev
      jINC = -1
      jABOVE = 0
      jBELOW = 1
    ENDIF

    CALL setup_cloud_optics(droplet_scale, &
      zinhoml1_, zinhoml2_, zinhoml3_, zinhomi_)
    CALL setup_flat_data(pressure_scale, droplet_scale)

    stpfac = stpfac_orig / pressure_scale
    log_factor = LOG(pressure_scale)
    log_pref = (/( &
      log_pref_tropopause_orig + log_factor - i * dlog_pref, &
      i = 1 - tropopause_ref, n_ref - tropopause_ref )/)
    jp_coeff = 1 + (log_pref(1) + jp_fudge) / dlog_pref

  END SUBROUTINE psrad_basic_setup


END MODULE mo_psrad_setup
