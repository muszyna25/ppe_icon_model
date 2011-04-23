!     path:      $Source: /storm/rc1/cvsroot/rc/rrtmg_lw/src/rrtmg_lw_rtrnmr.f90,v $
!     author:    $Author: mike $
!     revision:  $Revision: 1.7 $
!     created:   $Date: 2009/11/12 20:52:26 $
!
#ifdef __xlC__
@process hot
@process nostrict
#endif

MODULE mo_lrtm_rtrnmr

  !  --------------------------------------------------------------------------
  ! |                                                                          |
  ! |  Copyright 2002-2009, Atmospheric & Environmental Research, Inc. (AER).  |
  ! |  This software may be used, copied, or redistributed as long as it is    |
  ! |  not sold and this copyright notice is reproduced on each copy made.     |
  ! |  This model is provided as is without any express or implied warranties. |
  ! |                       (http://www.rtweb.aer.com/)                        |
  ! |                                                                          |
  !  --------------------------------------------------------------------------

  ! ------- Modules -------

  USE mo_kind,             ONLY : wp

#ifdef __ICON__
  USE mo_math_constants,   ONLY : api => pi
#else
  USE mo_constants,        ONLY : api
#endif

  USE mo_lrtm_par,         ONLY : nbndlw, delwave, ngs
#ifdef __TEST_lrtm_rtrnmr
  USE mo_timer,            ONLY: timer_lrtm_1, timer_lrtm_2, timer_start, timer_stop
  USE mo_run_nml,          ONLY: ltimer
  USE mo_exception,        ONLY: finish
#endif

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: ntbl, bpade, tau_tbl, exp_tbl, tfn_tbl, &
    & lrtm_rtrnmr

  REAL(wp), PARAMETER :: fluxfac = 2.0e+04_wp * api

  INTEGER,  PARAMETER :: ntbl = 10000
  REAL(wp), PARAMETER :: tblint = 10000.0_wp
  REAL(wp), DIMENSION(0:ntbl) :: tau_tbl
  REAL(wp), DIMENSION(0:ntbl) :: exp_tbl
  REAL(wp), DIMENSION(0:ntbl) :: tfn_tbl
  REAL(wp)                    :: bpade

CONTAINS

  !-----------------------------------------------------------------------------
#ifdef __TEST_lrtm_rtrnmr
  SUBROUTINE lrtm_rtrnmr_1(&
#else
  SUBROUTINE lrtm_rtrnmr(&
#endif
    & kproma, nlayers, istart, iend, iout, semiss, ncbands, &
    & cldfrac, taucloud, planklay, planklev, plankbnd, &
    & pwvcm, fracs, taut, &
    & totuflux, totdflux, fnet, &
    & totuclfl, totdclfl, fnetc, &
    & idrv, dplankbnd_dt, dtotuflux_dt, dtotuclfl_dt )
    !-----------------------------------------------------------------------------
    !
    !  Original version:   E. J. Mlawer, et al. RRTM_V3.0
    !  Revision for GCMs:  Michael J. Iacono; October, 2002
    !  Revision for F90:  Michael J. Iacono; June, 2006
    !  Revision for dFdT option: M. J. Iacono and E. J. Mlawer, November 2009
    !
    !  This program calculates the upward fluxes, downward fluxes, and
    !  heating rates for an arbitrary clear or cloudy atmosphere.  The input
    !  to this program is the atmospheric profile, all Planck function
    !  information, and the cloud fraction by layer.  A variable diffusivity
    !  angle (SECDIFF) is used for the angle integration.  Bands 2-3 and 5-9
    !  use a value for SECDIFF that varies from 1.50 to 1.80 as a function of
    !  the column water vapor, and other bands use a value of 1.66.  The Gaussian
    !  weight appropriate to this angle (WTDIFF=0.5) is applied here.  Note that
    !  use of the emissivity angle for the flux integration can cause errors of
    !  1 to 4 W/m2 within cloudy layers.
    !  Clouds are treated with a maximum-random cloud overlap method.
    !  This subroutine also provides the optional capability to calculate
    !  the derivative of upward flux respect to surface temperature using
    !  the pre-tabulated derivative of the Planck function with respect to
    !  temperature integrated over each spectral band.
    !***************************************************************************

    ! ------- Declarations -------

    ! ----- Input -----
    INTEGER, INTENT(in) :: kproma          ! number of columns
    INTEGER, INTENT(in) :: nlayers         ! total number of layers
    INTEGER, INTENT(in) :: istart          ! beginning band of calculation
    INTEGER, INTENT(in) :: iend            ! ending band of calculation
    INTEGER, INTENT(in) :: iout            ! output option flag

    ! Atmosphere
    !    Dimensions: (0:nlayers)
    REAL(wp), INTENT(in) :: pwvcm(:)           ! precipitable water vapor (cm)
    REAL(wp), INTENT(in) :: semiss(:,:)        ! lw surface emissivity
    !    Dimensions: (nbndlw)
    REAL(wp), INTENT(in) :: planklay(:,:,:)      !
    !    Dimensions: (nlayers,nbndlw)
    REAL(wp), INTENT(in) :: planklev(:,0:,:)     !
    !    Dimensions: (0:nlayers,nbndlw)
    REAL(wp), INTENT(in) :: plankbnd(:,:)        !
    !    Dimensions: (nbndlw)
    REAL(wp), INTENT(in) :: fracs(:,:,:)         !
    !    Dimensions: (kproma,nlayers,ngptw)
    REAL(wp), INTENT(in) :: taut(:,:,:)          ! gaseous + aerosol optical depths
    !    Dimensions: (kproma,nlayers,ngptlw)

    ! Clouds
    INTEGER, INTENT(in) :: ncbands(:)          ! number of cloud spectral bands
    REAL(wp), INTENT(in) :: cldfrac(:,:)       ! layer cloud fraction
    !    Dimensions: (nlayers)
    REAL(wp), INTENT(in) :: taucloud(:,:,:)      ! layer cloud optical depth
    !    Dimensions: (nlayers,nbndlw)
    INTEGER, INTENT(in) :: idrv            ! flag for calculation of dF/dt from
    ! Planck derivative [0=off, 1=on]
    REAL(wp), INTENT(in) :: dplankbnd_dt(:,:)    ! derivative of Planck function wrt temp
    !    Dimensions: (nbndlw)

    ! ----- Output -----
    REAL(wp), INTENT(out) :: totuflux(:,0:)      ! upward longwave flux (w/m2)
    !    Dimensions: (0:nlayers)
    REAL(wp), INTENT(out) :: totdflux(:,0:)      ! downward longwave flux (w/m2)
    !    Dimensions: (0:nlayers)
    REAL(wp), INTENT(out) :: fnet(:,0:)          ! net longwave flux (w/m2)
    !    Dimensions: (0:nlayers)
    REAL(wp), INTENT(out) :: totuclfl(:,0:)      ! clear sky upward longwave flux (w/m2)
    !    Dimensions: (0:nlayers)
    REAL(wp), INTENT(out) :: totdclfl(:,0:)      ! clear sky downward longwave flux (w/m2)
    !    Dimensions: (0:nlayers)
    REAL(wp), INTENT(out) :: fnetc(:,0:)         ! clear sky net longwave flux (w/m2)
    !    Dimensions: (0:nlayers)
    REAL(wp), INTENT(out) :: dtotuflux_dt(:,0:)  ! change in upward longwave flux (w/m2/k)
    ! with respect to surface temperature
    !    Dimensions: (0:nlayers)
    REAL(wp), INTENT(out) :: dtotuclfl_dt(:,0:)  ! change in upward longwave flux (w/m2/k)
    ! with respect to surface temperature
    !    Dimensions: (0:nlayers)

    ! Vectorized version implemented by Guenther Zaengl, DWD
    ! See above for variable descriptions
    ! ----- Local -----
    ! Declarations for radiative transfer
    REAL(wp) :: atot(kproma,nlayers)
    REAL(wp) :: atrans(kproma,nlayers)
    REAL(wp) :: bbugas(kproma,nlayers)
    REAL(wp) :: bbutot(kproma,nlayers)
    REAL(wp) :: clrurad(kproma,0:nlayers)
    REAL(wp) :: clrdrad(kproma,0:nlayers)
    REAL(wp) :: uflux
    REAL(wp) :: dflux
    REAL(wp) :: urad(kproma,0:nlayers)
    REAL(wp) :: drad(kproma,0:nlayers)
    REAL(wp) :: uclfl
    REAL(wp) :: dclfl

    REAL(wp) :: secdiff(kproma,nbndlw)          ! secant of diffusivity angle
    REAL(wp) :: a0(nbndlw),a1(nbndlw),a2(nbndlw)! diffusivity angle adjustment coefficients
    REAL(wp) :: wtdiff, rec_6, dplankup(kproma,nlayers), dplankdn(kproma,nlayers)
    REAL(wp) :: radld(kproma), radclrd(kproma), plfrac
    REAL(wp) :: odepth(kproma), odtot(kproma), odepth_rec, odtot_rec, gassrc(kproma), ttot
    REAL(wp) :: tblind, tfactot, bbd(kproma), bbdtot(kproma), tfacgas, transc, tausfac
    REAL(wp) :: rad0, reflect, radlu(kproma), radclru(kproma)

    REAL(wp) :: duflux_dt
    REAL(wp) :: duclfl_dt
    REAL(wp) :: d_urad_dt(kproma,0:nlayers)
    REAL(wp) :: d_clrurad_dt(kproma,0:nlayers)
    REAL(wp) :: d_rad0_dt, d_radlu_dt(kproma), d_radclru_dt(kproma)

    INTEGER :: icldlyr(kproma,nlayers)             ! flag for cloud in layer
    INTEGER :: ibnd, ib, iband, lay, lev    ! loop indices
    INTEGER :: igc                          ! g-point interval counter
    INTEGER :: iclddn(kproma)               ! flag for cloud in down path
    INTEGER :: ittot, itgas, itr            ! lookup table indices
    INTEGER :: ipat(16,0:2)
    INTEGER :: jl

    INTEGER :: ibv(kproma)

    ! Declarations for cloud overlap adjustment
    REAL(wp) :: faccld1(kproma,nlayers+1),faccld2(kproma,nlayers+1)
    REAL(wp) :: facclr1(kproma,nlayers+1),facclr2(kproma,nlayers+1)
    REAL(wp) :: faccmb1(kproma,nlayers+1),faccmb2(kproma,nlayers+1)
    REAL(wp) :: faccld1d(kproma,0:nlayers),faccld2d(kproma,0:nlayers)
    REAL(wp) :: facclr1d(kproma,0:nlayers),facclr2d(kproma,0:nlayers)
    REAL(wp) :: faccmb1d(kproma,0:nlayers),faccmb2d(kproma,0:nlayers)

    REAL(wp) :: fmax, fmin, rat1(kproma), rat2(kproma)
    REAL(wp), DIMENSION(kproma) :: clrradd, cldradd, clrradu, cldradu, oldclr, oldcld, &
      & rad, cldsrc, radmod

    INTEGER :: istcld(kproma,nlayers+1),istcldd(kproma,0:nlayers)

    ! Local variables for cloud / no cloud index lists
    INTEGER, DIMENSION(kproma,nlayers) :: icld_ind,iclear_ind
    INTEGER :: icld, iclear, n_cloudpoints(nlayers), n_clearpoints(nlayers)
    INTEGER :: icld1, npoints1, npoints2, npoints3, &
      &        ilist1(kproma), ilist2(kproma), ilist3(kproma)

    ! These arrays indicate the spectral 'region' (used in the
    ! calculation of ice cloud optical depths) corresponding
    ! to each spectral band.  See cldprop.f for more details.
    DATA ipat /1,1,1,1,1,1,1,1,1, 1, 1, 1, 1, 1, 1, 1, &
      & 1,2,3,3,3,4,4,4,5, 5, 5, 5, 5, 5, 5, 5, &
      & 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16/

    ! This secant and weight corresponds to the standard diffusivity
    ! angle.  This initial value is redefined below for some bands.
    DATA wtdiff /0.5_wp/
    DATA rec_6 /0.166667_wp/

    ! Reset diffusivity angle for Bands 2-3 and 5-9 to vary (between 1.50
    ! and 1.80) as a function of total column water vapor.  The function
    ! has been defined to minimize flux and cooling rate errors in these bands
    ! over a wide range of precipitable water values.
    DATA a0 / 1.66_wp,  1.55_wp,  1.58_wp,  1.66_wp, &
      & 1.54_wp, 1.454_wp,  1.89_wp,  1.33_wp, &
      & 1.668_wp,  1.66_wp,  1.66_wp,  1.66_wp, &
      & 1.66_wp,  1.66_wp,  1.66_wp,  1.66_wp /
    DATA a1 / 0.00_wp,  0.25_wp,  0.22_wp,  0.00_wp, &
      & 0.13_wp, 0.446_wp, -0.10_wp,  0.40_wp, &
      & -0.006_wp,  0.00_wp,  0.00_wp,  0.00_wp, &
      & 0.00_wp,  0.00_wp,  0.00_wp,  0.00_wp /
    DATA a2 / 0.00_wp, -12.0_wp, -11.7_wp,  0.00_wp, &
      & -0.72_wp,-0.243_wp,  0.19_wp,-0.062_wp, &
      & 0.414_wp,  0.00_wp,  0.00_wp,  0.00_wp, &
      & 0.00_wp,  0.00_wp,  0.00_wp,  0.00_wp /

    DO ibnd = 1,nbndlw
      IF (ibnd.EQ.1 .OR. ibnd.EQ.4 .OR. ibnd.GE.10) THEN
        DO jl = 1, kproma  ! loop over columns
          secdiff(jl,ibnd) = 1.66_wp
        ENDDO
      ELSE
        DO jl = 1, kproma  ! loop over columns
          secdiff(jl,ibnd) = a0(ibnd) + a1(ibnd)*EXP(a2(ibnd)*pwvcm(jl))
          IF (secdiff(jl,ibnd) .GT. 1.80_wp) secdiff(jl,ibnd) = 1.80_wp
          IF (secdiff(jl,ibnd) .LT. 1.50_wp) secdiff(jl,ibnd) = 1.50_wp
        ENDDO
      ENDIF
    ENDDO

!CDIR BEGIN COLLAPSE
    faccld1(:,:) = 0.0_wp
    faccld2(:,:) = 0.0_wp
    facclr1(:,:) = 0.0_wp
    facclr2(:,:) = 0.0_wp
    faccmb1(:,:) = 0.0_wp
    faccmb2(:,:) = 0.0_wp
    faccld1d(:,:) = 0.0_wp
    faccld2d(:,:) = 0.0_wp
    facclr1d(:,:) = 0.0_wp
    facclr2d(:,:) = 0.0_wp
    faccmb1d(:,:) = 0.0_wp
    faccmb2d(:,:) = 0.0_wp
    urad(:,:)     = 0.0_wp
    drad(:,:)     = 0.0_wp
    totuflux(:,:) = 0.0_wp
    totdflux(:,:) = 0.0_wp
    clrurad(:,:)  = 0.0_wp
    clrdrad(:,:)  = 0.0_wp
    totuclfl(:,:) = 0.0_wp
    totdclfl(:,:) = 0.0_wp
!CDIR END

    IF (idrv .EQ. 1) THEN
!CDIR BEGIN COLLAPSE
      d_urad_dt(:,:)    = 0.0_wp
      d_clrurad_dt(:,:) = 0.0_wp
      dtotuflux_dt(:,:) = 0.0_wp
      dtotuclfl_dt(:,:) = 0.0_wp
!CDIR END
    ENDIF

    DO lay = 1, nlayers

      icld   = 0
      iclear = 0

      DO jl = 1, kproma  ! loop over columns
        IF (cldfrac(jl,lay) .GE. 1.e-6_wp) THEN
          icldlyr(jl,lay) = 1
          icld = icld + 1
          icld_ind(icld,lay) = jl
        ELSE
          icldlyr(jl,lay) = 0
          iclear = iclear + 1
          iclear_ind(iclear,lay) = jl
        ENDIF
      ENDDO

      n_cloudpoints(lay) = icld
      n_clearpoints(lay) = iclear

    ENDDO

    ! Maximum/Random cloud overlap parameter

    istcld(:,1) = 1
    istcldd(:,nlayers) = 1

    DO lev = 1, nlayers
      IF (n_cloudpoints(lev) == kproma) THEN ! all points are cloudy
        DO jl = 1, kproma ! Thus, direct addressing can be used
          ! Maximum/random cloud overlap
          istcld(jl,lev+1) = 0
          IF (lev .EQ. nlayers) THEN
            faccld1(jl,lev+1) = 0._wp
            faccld2(jl,lev+1) = 0._wp
            facclr1(jl,lev+1) = 0._wp
            facclr2(jl,lev+1) = 0._wp
            faccmb1(jl,lev+1) = 0._wp
            faccmb2(jl,lev+1) = 0._wp
          ELSEIF (cldfrac(jl,lev+1) .GE. cldfrac(jl,lev)) THEN
            faccld1(jl,lev+1) = 0._wp
            faccld2(jl,lev+1) = 0._wp
            IF (istcld(jl,lev) .EQ. 1) THEN
              facclr1(jl,lev+1) = 0._wp
              facclr2(jl,lev+1) = 0._wp
              IF (cldfrac(jl,lev) .LT. 1._wp) facclr2(jl,lev+1) = &
                & (cldfrac(jl,lev+1)-cldfrac(jl,lev))/(1._wp-cldfrac(jl,lev))
              facclr2(jl,lev) = 0._wp
              faccld2(jl,lev) = 0._wp
            ELSE
              fmax = MAX(cldfrac(jl,lev),cldfrac(jl,lev-1))
              IF (cldfrac(jl,lev+1) .GT. fmax) THEN
                facclr1(jl,lev+1) = rat2(jl)
                facclr2(jl,lev+1) = (cldfrac(jl,lev+1)-fmax)/(1._wp-fmax)
              ELSEIF (cldfrac(jl,lev+1) .LT. fmax) THEN
                facclr1(jl,lev+1) = (cldfrac(jl,lev+1)-cldfrac(jl,lev))/ &
                  & (cldfrac(jl,lev-1)-cldfrac(jl,lev))
                facclr2(jl,lev+1) = 0._wp
              ELSE
                facclr1(jl,lev+1) = rat2(jl)
                facclr2(jl,lev+1) = 0._wp
              ENDIF
            ENDIF
            IF (facclr1(jl,lev+1).GT.0._wp .OR. facclr2(jl,lev+1).GT.0._wp) THEN
              rat1(jl) = 1._wp
              rat2(jl) = 0._wp
            ELSE
              rat1(jl) = 0._wp
              rat2(jl) = 0._wp
            ENDIF
          ELSE
            facclr1(jl,lev+1) = 0._wp
            facclr2(jl,lev+1) = 0._wp
            IF (istcld(jl,lev) .EQ. 1) THEN
              faccld1(jl,lev+1) = 0._wp
              faccld2(jl,lev+1) = (cldfrac(jl,lev)-cldfrac(jl,lev+1))/cldfrac(jl,lev)

              facclr2(jl,lev) = 0._wp
              faccld2(jl,lev) = 0._wp
            ELSE
              fmin = MIN(cldfrac(jl,lev),cldfrac(jl,lev-1))
              IF (cldfrac(jl,lev+1) .LE. fmin) THEN
                faccld1(jl,lev+1) = rat1(jl)
                faccld2(jl,lev+1) = (fmin-cldfrac(jl,lev+1))/fmin
              ELSE
                faccld1(jl,lev+1) = (cldfrac(jl,lev)-cldfrac(jl,lev+1))/(cldfrac(jl,lev)-fmin)
                faccld2(jl,lev+1) = 0._wp
              ENDIF
            ENDIF
            IF (faccld1(jl,lev+1).GT.0._wp .OR. faccld2(jl,lev+1).GT.0._wp) THEN
              rat1(jl) = 0._wp
              rat2(jl) = 1._wp
            ELSE
              rat1(jl) = 0._wp
              rat2(jl) = 0._wp
            ENDIF
          ENDIF
          IF (lev == 1) THEN
            faccmb1(jl,lev+1) = 0._wp
            faccmb2(jl,lev+1) = faccld1(jl,lev+1) * facclr2(jl,lev)
          ELSE
            faccmb1(jl,lev+1) = facclr1(jl,lev+1) * faccld2(jl,lev) * cldfrac(jl,lev-1)
            faccmb2(jl,lev+1) = faccld1(jl,lev+1) * facclr2(jl,lev) * (1._wp - cldfrac(jl,lev-1))
          ENDIF
        ENDDO
      ELSE IF (n_clearpoints(lev) == kproma) THEN ! all points are clear
        istcld(1:kproma,lev+1) = 1
      ELSE ! use index list for the case that not all points are cloudy
!CDIR NODEP,VOVERTAKE,VOB
        DO icld = 1, n_cloudpoints(lev)
          jl = icld_ind(icld,lev)
          ! Maximum/random cloud overlap
          istcld(jl,lev+1) = 0
          IF (lev .EQ. nlayers) THEN
            faccld1(jl,lev+1) = 0._wp
            faccld2(jl,lev+1) = 0._wp
            facclr1(jl,lev+1) = 0._wp
            facclr2(jl,lev+1) = 0._wp
            faccmb1(jl,lev+1) = 0._wp
            faccmb2(jl,lev+1) = 0._wp
          ELSEIF (cldfrac(jl,lev+1) .GE. cldfrac(jl,lev)) THEN
            faccld1(jl,lev+1) = 0._wp
            faccld2(jl,lev+1) = 0._wp
            IF (istcld(jl,lev) .EQ. 1) THEN
              facclr1(jl,lev+1) = 0._wp
              facclr2(jl,lev+1) = 0._wp
              IF (cldfrac(jl,lev) .LT. 1._wp) facclr2(jl,lev+1) = &
                & (cldfrac(jl,lev+1)-cldfrac(jl,lev))/(1._wp-cldfrac(jl,lev))
              facclr2(jl,lev) = 0._wp
              faccld2(jl,lev) = 0._wp
            ELSE
              fmax = MAX(cldfrac(jl,lev),cldfrac(jl,lev-1))
              IF (cldfrac(jl,lev+1) .GT. fmax) THEN
                facclr1(jl,lev+1) = rat2(jl)
                facclr2(jl,lev+1) = (cldfrac(jl,lev+1)-fmax)/(1._wp-fmax)
              ELSEIF (cldfrac(jl,lev+1) .LT. fmax) THEN
                facclr1(jl,lev+1) = (cldfrac(jl,lev+1)-cldfrac(jl,lev))/ &
                  & (cldfrac(jl,lev-1)-cldfrac(jl,lev))
                facclr2(jl,lev+1) = 0._wp
              ELSE
                facclr1(jl,lev+1) = rat2(jl)
                facclr2(jl,lev+1) = 0._wp
              ENDIF
            ENDIF
            IF (facclr1(jl,lev+1).GT.0._wp .OR. facclr2(jl,lev+1).GT.0._wp) THEN
              rat1(jl) = 1._wp
              rat2(jl) = 0._wp
            ELSE
              rat1(jl) = 0._wp
              rat2(jl) = 0._wp
            ENDIF
          ELSE
            facclr1(jl,lev+1) = 0._wp
            facclr2(jl,lev+1) = 0._wp
            IF (istcld(jl,lev) .EQ. 1) THEN
              faccld1(jl,lev+1) = 0._wp
              faccld2(jl,lev+1) = (cldfrac(jl,lev)-cldfrac(jl,lev+1))/cldfrac(jl,lev)

              facclr2(jl,lev) = 0._wp
              faccld2(jl,lev) = 0._wp
            ELSE
              fmin = MIN(cldfrac(jl,lev),cldfrac(jl,lev-1))
              IF (cldfrac(jl,lev+1) .LE. fmin) THEN
                faccld1(jl,lev+1) = rat1(jl)
                faccld2(jl,lev+1) = (fmin-cldfrac(jl,lev+1))/fmin
              ELSE
                faccld1(jl,lev+1) = (cldfrac(jl,lev)-cldfrac(jl,lev+1))/(cldfrac(jl,lev)-fmin)
                faccld2(jl,lev+1) = 0._wp
              ENDIF
            ENDIF
            IF (faccld1(jl,lev+1).GT.0._wp .OR. faccld2(jl,lev+1).GT.0._wp) THEN
              rat1(jl) = 0._wp
              rat2(jl) = 1._wp
            ELSE
              rat1(jl) = 0._wp
              rat2(jl) = 0._wp
            ENDIF
          ENDIF
          IF (lev == 1) THEN
            faccmb1(jl,lev+1) = 0._wp
            faccmb2(jl,lev+1) = faccld1(jl,lev+1) * facclr2(jl,lev)
          ELSE
            faccmb1(jl,lev+1) = facclr1(jl,lev+1) * faccld2(jl,lev) * cldfrac(jl,lev-1)
            faccmb2(jl,lev+1) = faccld1(jl,lev+1) * facclr2(jl,lev) * (1._wp - cldfrac(jl,lev-1))
          ENDIF
        ENDDO
!CDIR NODEP,VOVERTAKE,VOB
        DO iclear = 1, n_clearpoints(lev)
          jl = iclear_ind(iclear,lev)
          istcld(jl,lev+1) = 1
        ENDDO
      ENDIF

    ENDDO

    DO lev = nlayers, 1, -1

      IF (n_cloudpoints(lev) == kproma) THEN ! all points are cloudy
        DO jl = 1, kproma ! Thus, direct addressing can be used
          istcldd(jl,lev-1) = 0
          IF (lev .EQ. 1) THEN
            faccld1d(jl,lev-1) = 0._wp
            faccld2d(jl,lev-1) = 0._wp
            facclr1d(jl,lev-1) = 0._wp
            facclr2d(jl,lev-1) = 0._wp
            faccmb1d(jl,lev-1) = 0._wp
            faccmb2d(jl,lev-1) = 0._wp
          ELSEIF (cldfrac(jl,lev-1) .GE. cldfrac(jl,lev)) THEN
            faccld1d(jl,lev-1) = 0._wp
            faccld2d(jl,lev-1) = 0._wp
            IF (istcldd(jl,lev) .EQ. 1) THEN
              facclr1d(jl,lev-1) = 0._wp
              facclr2d(jl,lev-1) = 0._wp
              IF (cldfrac(jl,lev) .LT. 1._wp) facclr2d(jl,lev-1) = &
                & (cldfrac(jl,lev-1)-cldfrac(jl,lev))/(1._wp-cldfrac(jl,lev))
              facclr2d(jl,lev) = 0._wp
              faccld2d(jl,lev) = 0._wp
            ELSE
              fmax = MAX(cldfrac(jl,lev),cldfrac(jl,lev+1))
              IF (cldfrac(jl,lev-1) .GT. fmax) THEN
                facclr1d(jl,lev-1) = rat2(jl)
                facclr2d(jl,lev-1) = (cldfrac(jl,lev-1)-fmax)/(1._wp-fmax)
              ELSEIF (cldfrac(jl,lev-1) .LT. fmax) THEN
                facclr1d(jl,lev-1) = (cldfrac(jl,lev-1)-cldfrac(jl,lev))/ &
                  & (cldfrac(jl,lev+1)-cldfrac(jl,lev))
                facclr2d(jl,lev-1) = 0._wp
              ELSE
                facclr1d(jl,lev-1) = rat2(jl)
                facclr2d(jl,lev-1) = 0._wp
              ENDIF
            ENDIF
            IF (facclr1d(jl,lev-1).GT.0._wp .OR. facclr2d(jl,lev-1).GT.0._wp)THEN
              rat1(jl) = 1._wp
              rat2(jl) = 0._wp
            ELSE
              rat1(jl) = 0._wp
              rat2(jl) = 0._wp
            ENDIF
          ELSE
            facclr1d(jl,lev-1) = 0._wp
            facclr2d(jl,lev-1) = 0._wp
            IF (istcldd(jl,lev) .EQ. 1) THEN
              faccld1d(jl,lev-1) = 0._wp
              faccld2d(jl,lev-1) = (cldfrac(jl,lev)-cldfrac(jl,lev-1))/cldfrac(jl,lev)
              facclr2d(jl,lev) = 0._wp
              faccld2d(jl,lev) = 0._wp
            ELSE
              fmin = MIN(cldfrac(jl,lev),cldfrac(jl,lev+1))
              IF (cldfrac(jl,lev-1) .LE. fmin) THEN
                faccld1d(jl,lev-1) = rat1(jl)
                faccld2d(jl,lev-1) = (fmin-cldfrac(jl,lev-1))/fmin
              ELSE
                faccld1d(jl,lev-1) = (cldfrac(jl,lev)-cldfrac(jl,lev-1))/(cldfrac(jl,lev)-fmin)
                faccld2d(jl,lev-1) = 0._wp
              ENDIF
            ENDIF
            IF (faccld1d(jl,lev-1).GT.0._wp .OR. faccld2d(jl,lev-1).GT.0._wp)THEN
              rat1(jl) = 0._wp
              rat2(jl) = 1._wp
            ELSE
              rat1(jl) = 0._wp
              rat2(jl) = 0._wp
            ENDIF
          ENDIF
          IF (lev == nlayers) THEN
            faccmb1d(jl,lev-1) = 0._wp
            faccmb2d(jl,lev-1) = faccld1d(jl,lev-1) * facclr2d(jl,lev)
          ELSE
            faccmb1d(jl,lev-1) = facclr1d(jl,lev-1) * faccld2d(jl,lev) * cldfrac(jl,lev+1)
            faccmb2d(jl,lev-1) = faccld1d(jl,lev-1) * facclr2d(jl,lev) * (1._wp-cldfrac(jl,lev+1))
          ENDIF
        ENDDO
      ELSE IF (n_clearpoints(lev) == kproma) THEN ! all points are clear
        istcldd(1:kproma,lev-1) = 1
      ELSE  ! use index list for the case that not all points are cloudy
!CDIR NODEP,VOVERTAKE,VOB
        DO icld = 1, n_cloudpoints(lev)
          jl = icld_ind(icld,lev)

          istcldd(jl,lev-1) = 0
          IF (lev .EQ. 1) THEN
            faccld1d(jl,lev-1) = 0._wp
            faccld2d(jl,lev-1) = 0._wp
            facclr1d(jl,lev-1) = 0._wp
            facclr2d(jl,lev-1) = 0._wp
            faccmb1d(jl,lev-1) = 0._wp
            faccmb2d(jl,lev-1) = 0._wp
          ELSEIF (cldfrac(jl,lev-1) .GE. cldfrac(jl,lev)) THEN
            faccld1d(jl,lev-1) = 0._wp
            faccld2d(jl,lev-1) = 0._wp
            IF (istcldd(jl,lev) .EQ. 1) THEN
              facclr1d(jl,lev-1) = 0._wp
              facclr2d(jl,lev-1) = 0._wp
              IF (cldfrac(jl,lev) .LT. 1._wp) facclr2d(jl,lev-1) = &
                & (cldfrac(jl,lev-1)-cldfrac(jl,lev))/(1._wp-cldfrac(jl,lev))
              facclr2d(jl,lev) = 0._wp
              faccld2d(jl,lev) = 0._wp
            ELSE
              fmax = MAX(cldfrac(jl,lev),cldfrac(jl,lev+1))
              IF (cldfrac(jl,lev-1) .GT. fmax) THEN
                facclr1d(jl,lev-1) = rat2(jl)
                facclr2d(jl,lev-1) = (cldfrac(jl,lev-1)-fmax)/(1._wp-fmax)
              ELSEIF (cldfrac(jl,lev-1) .LT. fmax) THEN
                facclr1d(jl,lev-1) = (cldfrac(jl,lev-1)-cldfrac(jl,lev))/ &
                  & (cldfrac(jl,lev+1)-cldfrac(jl,lev))
                facclr2d(jl,lev-1) = 0._wp
              ELSE
                facclr1d(jl,lev-1) = rat2(jl)
                facclr2d(jl,lev-1) = 0._wp
              ENDIF
            ENDIF
            IF (facclr1d(jl,lev-1).GT.0._wp .OR. facclr2d(jl,lev-1).GT.0._wp)THEN
              rat1(jl) = 1._wp
              rat2(jl) = 0._wp
            ELSE
              rat1(jl) = 0._wp
              rat2(jl) = 0._wp
            ENDIF
          ELSE
            facclr1d(jl,lev-1) = 0._wp
            facclr2d(jl,lev-1) = 0._wp
            IF (istcldd(jl,lev) .EQ. 1) THEN
              faccld1d(jl,lev-1) = 0._wp
              faccld2d(jl,lev-1) = (cldfrac(jl,lev)-cldfrac(jl,lev-1))/cldfrac(jl,lev)
              facclr2d(jl,lev) = 0._wp
              faccld2d(jl,lev) = 0._wp
            ELSE
              fmin = MIN(cldfrac(jl,lev),cldfrac(jl,lev+1))
              IF (cldfrac(jl,lev-1) .LE. fmin) THEN
                faccld1d(jl,lev-1) = rat1(jl)
                faccld2d(jl,lev-1) = (fmin-cldfrac(jl,lev-1))/fmin
              ELSE
                faccld1d(jl,lev-1) = (cldfrac(jl,lev)-cldfrac(jl,lev-1))/(cldfrac(jl,lev)-fmin)
                faccld2d(jl,lev-1) = 0._wp
              ENDIF
            ENDIF
            IF (faccld1d(jl,lev-1).GT.0._wp .OR. faccld2d(jl,lev-1).GT.0._wp)THEN
              rat1(jl) = 0._wp
              rat2(jl) = 1._wp
            ELSE
              rat1(jl) = 0._wp
              rat2(jl) = 0._wp
            ENDIF
          ENDIF
          IF (lev == nlayers) THEN
            faccmb1d(jl,lev-1) = 0._wp
            faccmb2d(jl,lev-1) = faccld1d(jl,lev-1) * facclr2d(jl,lev)
          ELSE
            faccmb1d(jl,lev-1) = facclr1d(jl,lev-1) * faccld2d(jl,lev) * cldfrac(jl,lev+1)
            faccmb2d(jl,lev-1) = faccld1d(jl,lev-1) * facclr2d(jl,lev) * (1._wp-cldfrac(jl,lev+1))
          ENDIF
        ENDDO
!CDIR NODEP,VOVERTAKE,VOB
        DO iclear = 1, n_clearpoints(lev)
          jl = iclear_ind(iclear,lev)
          istcldd(jl,lev-1) = 1
        ENDDO
      ENDIF

    ENDDO

    igc = 1
    ! Loop over frequency bands.
    DO iband = istart, iend

      DO jl = 1, kproma
        IF (ncbands(jl) .EQ. 1) THEN
          ibv(jl) = ipat(iband,0)
        ELSEIF (ncbands(jl) .EQ.  5) THEN
          ibv(jl) = ipat(iband,1)
        ELSEIF (ncbands(jl) .EQ. 16) THEN
          ibv(jl) = ipat(iband,2)
        ENDIF
      ENDDO

      DO lev = 1, nlayers
        DO jl = 1, kproma
          dplankup(jl,lev) = planklev(jl,lev,iband) - planklay(jl,lev,iband)
          dplankdn(jl,lev) = planklev(jl,lev-1,iband) - planklay(jl,lev,iband)
        ENDDO
      ENDDO

      ! Reinitialize g-point counter for each band if output for each band is requested.
      IF (iout.GT.0.AND.iband.GE.2) igc = ngs(iband-1)+1

      ! Loop over g-channels.
      1000    CONTINUE

      ! Radiative transfer starts here.
      radld(:) = 0._wp
      radclrd(:) = 0._wp
      iclddn(:) = 0

      ! Downward radiative transfer loop.
      DO lev = nlayers, 1, -1
        IF (n_cloudpoints(lev) == kproma) THEN ! all points are cloudy

          DO jl = 1, kproma ! Thus, direct addressing can be used

            ib = ibv(jl)
            plfrac = fracs(jl,lev,igc)
            odepth(jl) = MAX(0.0_wp, secdiff(jl,iband) * taut(jl,lev,igc))

            iclddn(jl) = 1
            odtot(jl) = odepth(jl) + secdiff(jl,ib) * taucloud(jl,lev,ib)
            IF (odtot(jl) .LT. 0.06_wp) THEN
              atrans(jl,lev) = odepth(jl) - 0.5_wp*odepth(jl)*odepth(jl)
              odepth_rec = rec_6*odepth(jl)
              gassrc(jl) = plfrac*(planklay(jl,lev,iband) &
                         + dplankdn(jl,lev)*odepth_rec)*atrans(jl,lev)

              atot(jl,lev) = odtot(jl) - 0.5_wp*odtot(jl)*odtot(jl)
              odtot_rec = rec_6*odtot(jl)
              bbdtot(jl) =  plfrac * (planklay(jl,lev,iband)+dplankdn(jl,lev)*odtot_rec)
              bbd(jl) = plfrac*(planklay(jl,lev,iband)+dplankdn(jl,lev)*odepth_rec)

              bbugas(jl,lev) =  plfrac * (planklay(jl,lev,iband)+dplankup(jl,lev)*odepth_rec)
              bbutot(jl,lev) =  plfrac * (planklay(jl,lev,iband)+dplankup(jl,lev)*odtot_rec)
            ELSEIF (odepth(jl) .LE. 0.06_wp) THEN
              atrans(jl,lev) = odepth(jl) - 0.5_wp*odepth(jl)*odepth(jl)
              odepth_rec = rec_6*odepth(jl)
              gassrc(jl) = plfrac*(planklay(jl,lev,iband) &
                         + dplankdn(jl,lev)*odepth_rec)*atrans(jl,lev)

              tblind = odtot(jl)/(bpade+odtot(jl))
              ittot = INT(tblint*tblind + 0.5_wp)
              tfactot = tfn_tbl(ittot)
              bbdtot(jl) = plfrac * (planklay(jl,lev,iband) + tfactot*dplankdn(jl,lev))
              bbd(jl) = plfrac*(planklay(jl,lev,iband)+dplankdn(jl,lev)*odepth_rec)
              atot(jl,lev) = 1._wp - exp_tbl(ittot)

              bbugas(jl,lev) = plfrac * (planklay(jl,lev,iband) + dplankup(jl,lev)*odepth_rec)
              bbutot(jl,lev) = plfrac * (planklay(jl,lev,iband) + tfactot * dplankup(jl,lev))
            ELSE
              tblind = odepth(jl)/(bpade+odepth(jl))
              itgas = INT(tblint*tblind+0.5_wp)
              odepth(jl) = tau_tbl(itgas)
              atrans(jl,lev) = 1._wp - exp_tbl(itgas)
              tfacgas = tfn_tbl(itgas)
              gassrc(jl) = atrans(jl,lev) * plfrac * (planklay(jl,lev,iband) &
                                                   + tfacgas*dplankdn(jl,lev))

              tblind = odtot(jl)/(bpade+odtot(jl))
              ittot = INT(tblint*tblind + 0.5_wp)
              tfactot = tfn_tbl(ittot)
              bbdtot(jl) = plfrac * (planklay(jl,lev,iband) + tfactot*dplankdn(jl,lev))
              bbd(jl) = plfrac*(planklay(jl,lev,iband)+tfacgas*dplankdn(jl,lev))
              atot(jl,lev) = 1._wp - exp_tbl(ittot)

              bbugas(jl,lev) = plfrac * (planklay(jl,lev,iband) + tfacgas * dplankup(jl,lev))
              bbutot(jl,lev) = plfrac * (planklay(jl,lev,iband) + tfactot * dplankup(jl,lev))
            ENDIF

            IF (istcldd(jl,lev) .EQ. 1) THEN
              cldradd(jl) = cldfrac(jl,lev) * radld(jl)
              clrradd(jl) = radld(jl) - cldradd(jl)
              oldcld(jl) = cldradd(jl)
              oldclr(jl) = clrradd(jl)
              rad(jl) = 0._wp
            ENDIF
            ttot = 1._wp - atot(jl,lev)
            cldsrc(jl) = bbdtot(jl) * atot(jl,lev)
            cldradd(jl) = cldradd(jl) * ttot + cldfrac(jl,lev) * cldsrc(jl)
            clrradd(jl) = clrradd(jl) * (1._wp-atrans(jl,lev)) + &
              & (1._wp-cldfrac(jl,lev))*gassrc(jl)
            radld(jl) = cldradd(jl) + clrradd(jl)
            drad(jl,lev-1) = drad(jl,lev-1) + radld(jl)

            radmod(jl) = rad(jl) * &
              & (facclr1d(jl,lev-1) * (1._wp-atrans(jl,lev)) + &
              & faccld1d(jl,lev-1) *  ttot) - &
              & faccmb1d(jl,lev-1) * gassrc(jl) + &
              & faccmb2d(jl,lev-1) * cldsrc(jl)

            oldcld(jl) = cldradd(jl) - radmod(jl)
            oldclr(jl) = clrradd(jl) + radmod(jl)
            rad(jl) = -radmod(jl) + facclr2d(jl,lev-1)*oldclr(jl) -&
              &  faccld2d(jl,lev-1)*oldcld(jl)
            cldradd(jl) = cldradd(jl) + rad(jl)
            clrradd(jl) = clrradd(jl) - rad(jl)
          ENDDO

        ELSE IF (n_clearpoints(lev) == kproma) THEN ! all points are clear

          DO jl = 1, kproma ! Thus, direct addressing can be used

            plfrac = fracs(jl,lev,igc)
            odepth(jl) = MAX(0.0_wp, secdiff(jl,iband) * taut(jl,lev,igc))

            IF (odepth(jl) .LE. 0.06_wp) THEN
              atrans(jl,lev) = odepth(jl)-0.5_wp*odepth(jl)*odepth(jl)
              odepth_rec = rec_6*odepth(jl)
              bbd(jl) = plfrac*(planklay(jl,lev,iband)+dplankdn(jl,lev)*odepth_rec)
              bbugas(jl,lev) = plfrac*(planklay(jl,lev,iband)+dplankup(jl,lev)*odepth_rec)
            ELSE
              tblind = odepth(jl)/(bpade+odepth(jl))
              itr = INT(tblint*tblind+0.5_wp)
              transc = exp_tbl(itr)
              atrans(jl,lev) = 1._wp-transc
              tausfac = tfn_tbl(itr)
              bbd(jl) = plfrac*(planklay(jl,lev,iband)+tausfac*dplankdn(jl,lev))
              bbugas(jl,lev) = plfrac * (planklay(jl,lev,iband) + tausfac * dplankup(jl,lev))
            ENDIF
            radld(jl) = radld(jl) + (bbd(jl)-radld(jl))*atrans(jl,lev)
            drad(jl,lev-1) = drad(jl,lev-1) + radld(jl)
          ENDDO

        ELSE ! both cloudy and clear points are in the vector

          npoints1 = 0
          npoints2 = 0
          npoints3 = 0
 
          ! Cloudy layer
!CDIR NODEP,VOVERTAKE,VOB
          DO icld = 1, n_cloudpoints(lev)
            jl = icld_ind(icld,lev)

            ib = ibv(jl)
            odepth(jl) = MAX(0.0_wp, secdiff(jl,iband) * taut(jl,lev,igc))

            iclddn(jl) = 1
            odtot(jl) = odepth(jl) + secdiff(jl,ib) * taucloud(jl,lev,ib)

            IF (odtot(jl) .LT. 0.06_wp) THEN
              npoints1 = npoints1 + 1
              ilist1(npoints1) = jl
            ELSEIF (odepth(jl) .LE. 0.06_wp) THEN
              npoints2 = npoints2 + 1
              ilist2(npoints2) = jl
            ELSE
              npoints3 = npoints3 + 1
              ilist3(npoints3) = jl
            ENDIF
          ENDDO

!CDIR NODEP,VOVERTAKE,VOB
          DO icld1 = 1, npoints1
            jl = ilist1(icld1)

            plfrac = fracs(jl,lev,igc)
            atrans(jl,lev) = odepth(jl) - 0.5_wp*odepth(jl)*odepth(jl)
            odepth_rec = rec_6*odepth(jl)
            gassrc(jl) = plfrac*(planklay(jl,lev,iband) &
                       + dplankdn(jl,lev)*odepth_rec)*atrans(jl,lev)

            atot(jl,lev) = odtot(jl) - 0.5_wp*odtot(jl)*odtot(jl)
            odtot_rec = rec_6*odtot(jl)
            bbdtot(jl) =  plfrac * (planklay(jl,lev,iband)+dplankdn(jl,lev)*odtot_rec)
            bbd(jl) = plfrac*(planklay(jl,lev,iband)+dplankdn(jl,lev)*odepth_rec)

            bbugas(jl,lev) =  plfrac * (planklay(jl,lev,iband)+dplankup(jl,lev)*odepth_rec)
            bbutot(jl,lev) =  plfrac * (planklay(jl,lev,iband)+dplankup(jl,lev)*odtot_rec)
          ENDDO

!CDIR NODEP,VOVERTAKE,VOB
          DO icld1 = 1, npoints2
            jl = ilist2(icld1)

            plfrac = fracs(jl,lev,igc)
            atrans(jl,lev) = odepth(jl) - 0.5_wp*odepth(jl)*odepth(jl)
            odepth_rec = rec_6*odepth(jl)
            gassrc(jl) = plfrac*(planklay(jl,lev,iband) &
                       + dplankdn(jl,lev)*odepth_rec)*atrans(jl,lev)

            tblind = odtot(jl)/(bpade+odtot(jl))
            ittot = INT(tblint*tblind + 0.5_wp)
            tfactot = tfn_tbl(ittot)
            bbdtot(jl) = plfrac * (planklay(jl,lev,iband) + tfactot*dplankdn(jl,lev))
            bbd(jl) = plfrac*(planklay(jl,lev,iband)+dplankdn(jl,lev)*odepth_rec)
            atot(jl,lev) = 1._wp - exp_tbl(ittot)

            bbugas(jl,lev) = plfrac * (planklay(jl,lev,iband) + dplankup(jl,lev)*odepth_rec)
            bbutot(jl,lev) = plfrac * (planklay(jl,lev,iband) + tfactot * dplankup(jl,lev))
          ENDDO

!CDIR NODEP,VOVERTAKE,VOB
          DO icld1 = 1, npoints3
            jl = ilist3(icld1)

            plfrac = fracs(jl,lev,igc)
            tblind = odepth(jl)/(bpade+odepth(jl))
            itgas = INT(tblint*tblind+0.5_wp)
            odepth(jl) = tau_tbl(itgas)
            atrans(jl,lev) = 1._wp - exp_tbl(itgas)
            tfacgas = tfn_tbl(itgas)
            gassrc(jl) = atrans(jl,lev) * plfrac * (planklay(jl,lev,iband) &
                                                 + tfacgas*dplankdn(jl,lev))

            tblind = odtot(jl)/(bpade+odtot(jl))
            ittot = INT(tblint*tblind + 0.5_wp)
            tfactot = tfn_tbl(ittot)
            bbdtot(jl) = plfrac * (planklay(jl,lev,iband) + tfactot*dplankdn(jl,lev))
            bbd(jl) = plfrac*(planklay(jl,lev,iband)+tfacgas*dplankdn(jl,lev))
            atot(jl,lev) = 1._wp - exp_tbl(ittot)

            bbugas(jl,lev) = plfrac * (planklay(jl,lev,iband) + tfacgas * dplankup(jl,lev))
            bbutot(jl,lev) = plfrac * (planklay(jl,lev,iband) + tfactot * dplankup(jl,lev))
          ENDDO

!CDIR NODEP,VOVERTAKE,VOB
          DO icld = 1, n_cloudpoints(lev)
            jl = icld_ind(icld,lev)

            IF (istcldd(jl,lev) .EQ. 1) THEN
              cldradd(jl) = cldfrac(jl,lev) * radld(jl)
              clrradd(jl) = radld(jl) - cldradd(jl)
              oldcld(jl) = cldradd(jl)
              oldclr(jl) = clrradd(jl)
              rad(jl) = 0._wp
            ENDIF
            ttot = 1._wp - atot(jl,lev)
            cldsrc(jl) = bbdtot(jl) * atot(jl,lev)
            cldradd(jl) = cldradd(jl) * ttot + cldfrac(jl,lev) * cldsrc(jl)
            clrradd(jl) = clrradd(jl) * (1._wp-atrans(jl,lev)) + &
              & (1._wp-cldfrac(jl,lev))*gassrc(jl)
            radld(jl) = cldradd(jl) + clrradd(jl)
            drad(jl,lev-1) = drad(jl,lev-1) + radld(jl)

            radmod(jl) = rad(jl) * &
              & (facclr1d(jl,lev-1) * (1._wp-atrans(jl,lev)) + &
              & faccld1d(jl,lev-1) *  ttot) - &
              & faccmb1d(jl,lev-1) * gassrc(jl) + &
              & faccmb2d(jl,lev-1) * cldsrc(jl)

            oldcld(jl) = cldradd(jl) - radmod(jl)
            oldclr(jl) = clrradd(jl) + radmod(jl)
            rad(jl) = -radmod(jl) + facclr2d(jl,lev-1)*oldclr(jl) -&
              &  faccld2d(jl,lev-1)*oldcld(jl)
            cldradd(jl) = cldradd(jl) + rad(jl)
            clrradd(jl) = clrradd(jl) - rad(jl)
          ENDDO

          ! Clear layer
!CDIR NODEP,VOVERTAKE,VOB
          DO iclear = 1, n_clearpoints(lev)
            jl = iclear_ind(iclear,lev)

            plfrac = fracs(jl,lev,igc)
            odepth(jl) = MAX(0.0_wp, secdiff(jl,iband) * taut(jl,lev,igc))

            IF (odepth(jl) .LE. 0.06_wp) THEN
              atrans(jl,lev) = odepth(jl)-0.5_wp*odepth(jl)*odepth(jl)
              odepth_rec = rec_6*odepth(jl)
              bbd(jl) = plfrac*(planklay(jl,lev,iband)+dplankdn(jl,lev)*odepth_rec)
              bbugas(jl,lev) = plfrac*(planklay(jl,lev,iband)+dplankup(jl,lev)*odepth_rec)
            ELSE
              tblind = odepth(jl)/(bpade+odepth(jl))
              itr = INT(tblint*tblind+0.5_wp)
              transc = exp_tbl(itr)
              atrans(jl,lev) = 1._wp-transc
              tausfac = tfn_tbl(itr)
              bbd(jl) = plfrac*(planklay(jl,lev,iband)+tausfac*dplankdn(jl,lev))
              bbugas(jl,lev) = plfrac * (planklay(jl,lev,iband) + tausfac * dplankup(jl,lev))
            ENDIF
            radld(jl) = radld(jl) + (bbd(jl)-radld(jl))*atrans(jl,lev)
            drad(jl,lev-1) = drad(jl,lev-1) + radld(jl)
          ENDDO

        ENDIF

        !  Set clear sky stream to total sky stream as long as layers
        !  remain clear.  Streams diverge when a cloud is reached (iclddn=1),
        !  and clear sky stream must be computed separately from that point.
        DO jl = 1, kproma
          IF (iclddn(jl).EQ.1) THEN
            radclrd(jl) = radclrd(jl) + (bbd(jl)-radclrd(jl)) * atrans(jl,lev)
            clrdrad(jl,lev-1) = clrdrad(jl,lev-1) + radclrd(jl)
          ELSE
            radclrd(jl) = radld(jl)
            clrdrad(jl,lev-1) = drad(jl,lev-1)
          ENDIF
        ENDDO
      ENDDO

      DO jl = 1, kproma
        ! Spectral emissivity & reflectance
        !  Include the contribution of spectrally varying longwave emissivity
        !  and reflection from the surface to the upward radiative transfer.
        !  Note: Spectral and Lambertian reflection are identical for the
        !  diffusivity angle flux integration used here.
        !  Note: The emissivity is applied to plankbnd and dplankbnd_dt when
        !  they are defined in subroutine setcoef.

        rad0 = fracs(jl,1,igc) * plankbnd(jl,iband)

        !  Add in reflection of surface downward radiance.
        reflect = 1._wp - semiss(jl,iband)
        radlu(jl) = rad0 + reflect * radld(jl)
        radclru(jl) = rad0 + reflect * radclrd(jl)

        ! Upward radiative transfer loop.

        urad(jl,0) = urad(jl,0) + radlu(jl)
        clrurad(jl,0) = clrurad(jl,0) + radclru(jl)
      ENDDO
      IF (idrv .EQ. 1) THEN
        DO jl = 1, kproma
          d_rad0_dt = fracs(jl,1,igc) * dplankbnd_dt(jl,iband)
          d_radlu_dt(jl) = d_rad0_dt
          d_urad_dt(jl,0) = d_urad_dt(jl,0) + d_radlu_dt(jl)
          d_radclru_dt(jl) = d_rad0_dt
          d_clrurad_dt(jl,0) = d_clrurad_dt(jl,0) + d_radclru_dt(jl)
        ENDDO
      ENDIF

      DO lev = 1, nlayers
        IF (n_cloudpoints(lev) == kproma) THEN ! all points are cloudy

          DO jl = 1, kproma ! Thus, direct addressing can be used

            gassrc(jl) = bbugas(jl,lev) * atrans(jl,lev)
            IF (istcld(jl,lev) .EQ. 1) THEN
              cldradu(jl) = cldfrac(jl,lev) * radlu(jl)
              clrradu(jl) = radlu(jl) - cldradu(jl)
              oldcld(jl) = cldradu(jl)
              oldclr(jl) = clrradu(jl)
              rad(jl) = 0._wp
            ENDIF
            ttot = 1._wp - atot(jl,lev)
            cldsrc(jl) = bbutot(jl,lev) * atot(jl,lev)
            cldradu(jl) = cldradu(jl) * ttot + cldfrac(jl,lev) * cldsrc(jl)
            clrradu(jl) = clrradu(jl) * (1.0_wp-atrans(jl,lev))+(1._wp-cldfrac(jl,lev))*gassrc(jl)
            ! Total sky radiance
            radlu(jl) = cldradu(jl) + clrradu(jl)
            urad(jl,lev) = urad(jl,lev) + radlu(jl)
            radmod(jl) = rad(jl) * &
              & (facclr1(jl,lev+1)*(1.0_wp-atrans(jl,lev))+ &
              & faccld1(jl,lev+1) *  ttot) - &
              & faccmb1(jl,lev+1) * gassrc(jl) + &
              & faccmb2(jl,lev+1) * cldsrc(jl)
            oldcld(jl) = cldradu(jl) - radmod(jl)
            oldclr(jl) = clrradu(jl) + radmod(jl)
            rad(jl) = -radmod(jl) + facclr2(jl,lev+1)*oldclr(jl) - faccld2(jl,lev+1)*oldcld(jl)
            cldradu(jl) = cldradu(jl) + rad(jl)
            clrradu(jl) = clrradu(jl) - rad(jl)
          ENDDO
          IF (idrv .EQ. 1) THEN
            DO jl = 1, kproma
              d_radlu_dt(jl) = d_radlu_dt(jl) * cldfrac(jl,lev) * (1.0_wp - atot(jl,lev)) + &
                & d_radlu_dt(jl) * (1.0_wp - cldfrac(jl,lev)) * (1.0_wp - atrans(jl,lev))
              d_urad_dt(jl,lev) = d_urad_dt(jl,lev) + d_radlu_dt(jl)
            ENDDO
          ENDIF

        ELSE IF (n_clearpoints(lev) == kproma) THEN ! all points are clear

          DO jl = 1, kproma ! thus, direct addressing can be used
            radlu(jl) = radlu(jl) + (bbugas(jl,lev)-radlu(jl))*atrans(jl,lev)
            urad(jl,lev) = urad(jl,lev) + radlu(jl)
          ENDDO
          IF (idrv .EQ. 1) THEN
            DO jl = 1, kproma
              d_radlu_dt(jl) = d_radlu_dt(jl) * (1.0_wp - atrans(jl,lev))
              d_urad_dt(jl,lev) = d_urad_dt(jl,lev) + d_radlu_dt(jl)
            ENDDO
          ENDIF

        ELSE ! both cloudy and clear points are in the vector

          ! Cloudy layer
!CDIR NODEP,VOVERTAKE,VOB
          DO icld = 1, n_cloudpoints(lev)
            jl = icld_ind(icld,lev)

            gassrc(jl) = bbugas(jl,lev) * atrans(jl,lev)
            IF (istcld(jl,lev) .EQ. 1) THEN
              cldradu(jl) = cldfrac(jl,lev) * radlu(jl)
              clrradu(jl) = radlu(jl) - cldradu(jl)
              oldcld(jl) = cldradu(jl)
              oldclr(jl) = clrradu(jl)
              rad(jl) = 0._wp
            ENDIF
            ttot = 1._wp - atot(jl,lev)
            cldsrc(jl) = bbutot(jl,lev) * atot(jl,lev)
            cldradu(jl) = cldradu(jl) * ttot + cldfrac(jl,lev) * cldsrc(jl)
            clrradu(jl) = clrradu(jl) * (1.0_wp-atrans(jl,lev))+(1._wp-cldfrac(jl,lev))*gassrc(jl)
            ! Total sky radiance
            radlu(jl) = cldradu(jl) + clrradu(jl)
            urad(jl,lev) = urad(jl,lev) + radlu(jl)
            radmod(jl) = rad(jl) * &
              & (facclr1(jl,lev+1)*(1.0_wp-atrans(jl,lev))+ &
              & faccld1(jl,lev+1) *  ttot) - &
              & faccmb1(jl,lev+1) * gassrc(jl) + &
              & faccmb2(jl,lev+1) * cldsrc(jl)
            oldcld(jl) = cldradu(jl) - radmod(jl)
            oldclr(jl) = clrradu(jl) + radmod(jl)
            rad(jl) = -radmod(jl) + facclr2(jl,lev+1)*oldclr(jl) - faccld2(jl,lev+1)*oldcld(jl)
            cldradu(jl) = cldradu(jl) + rad(jl)
            clrradu(jl) = clrradu(jl) - rad(jl)
          ENDDO

          ! Clear layer
!CDIR NODEP,VOVERTAKE,VOB
          DO iclear = 1, n_clearpoints(lev)
            jl = iclear_ind(iclear,lev)

            radlu(jl) = radlu(jl) + (bbugas(jl,lev)-radlu(jl))*atrans(jl,lev)
            urad(jl,lev) = urad(jl,lev) + radlu(jl)
          ENDDO

          IF (idrv .EQ. 1) THEN
            DO jl = 1, kproma
              IF (icldlyr(jl,lev) .EQ. 1) THEN
                d_radlu_dt(jl) = d_radlu_dt(jl) * cldfrac(jl,lev) * (1.0_wp - atot(jl,lev)) + &
                  & d_radlu_dt(jl) * (1.0_wp - cldfrac(jl,lev)) * (1.0_wp - atrans(jl,lev))
                d_urad_dt(jl,lev) = d_urad_dt(jl,lev) + d_radlu_dt(jl)
              ELSE
                d_radlu_dt(jl) = d_radlu_dt(jl) * (1.0_wp - atrans(jl,lev))
                d_urad_dt(jl,lev) = d_urad_dt(jl,lev) + d_radlu_dt(jl)
              ENDIF
            ENDDO
          ENDIF

        ENDIF
        !  Set clear sky stream to total sky stream as long as all layers
        !  are clear (iclddn=0).  Streams must be calculated separately at
        !  all layers when a cloud is present (iclddn=1), because surface
        !  reflectance is different for each stream.
        DO jl = 1, kproma
          IF (iclddn(jl).EQ.1) THEN
            radclru(jl) = radclru(jl) + (bbugas(jl,lev)-radclru(jl))*atrans(jl,lev)
            clrurad(jl,lev) = clrurad(jl,lev) + radclru(jl)
          ELSE
            radclru(jl) = radlu(jl)
            clrurad(jl,lev) = urad(jl,lev)
          ENDIF
        ENDDO
        IF (idrv .EQ. 1) THEN
          DO jl = 1, kproma
            IF (iclddn(jl).EQ.1) THEN
              d_radclru_dt(jl) = d_radclru_dt(jl) * (1.0_wp - atrans(jl,lev))
              d_clrurad_dt(jl,lev) = d_clrurad_dt(jl,lev) + d_radclru_dt(jl)
            ELSE
              d_radclru_dt(jl) = d_radlu_dt(jl)
              d_clrurad_dt(jl,lev) = d_urad_dt(jl,lev)
            ENDIF
          ENDDO
        ENDIF
      ENDDO

      ! Increment g-point counter
      igc = igc + 1
      ! Return to continue radiative transfer for all g-channels in present band
      IF (igc .LE. ngs(iband)) GO TO 1000

      ! Process longwave output from band.
      ! Calculate upward, downward, and net flux.
      DO lev = nlayers, 0, -1
        DO jl = 1, kproma
          uflux = urad(jl,lev)*wtdiff
          dflux = drad(jl,lev)*wtdiff
          urad(jl,lev) = 0.0_wp
          drad(jl,lev) = 0.0_wp
          totuflux(jl,lev) = totuflux(jl,lev) + uflux * delwave(iband)
          totdflux(jl,lev) = totdflux(jl,lev) + dflux * delwave(iband)
          uclfl = clrurad(jl,lev)*wtdiff
          dclfl = clrdrad(jl,lev)*wtdiff
          clrurad(jl,lev) = 0.0_wp
          clrdrad(jl,lev) = 0.0_wp
          totuclfl(jl,lev) = totuclfl(jl,lev) + uclfl * delwave(iband)
          totdclfl(jl,lev) = totdclfl(jl,lev) + dclfl * delwave(iband)
        ENDDO
      ENDDO

      ! Calculate total change in upward flux wrt surface temperature
      IF (idrv .EQ. 1) THEN
        DO lev = nlayers, 0, -1
          DO jl = 1, kproma
            duflux_dt = d_urad_dt(jl,lev) * wtdiff
            d_urad_dt(jl,lev) = 0.0_wp
            dtotuflux_dt(jl,lev) = dtotuflux_dt(jl,lev) + duflux_dt * delwave(iband) * fluxfac

            duclfl_dt = d_clrurad_dt(jl,lev) * wtdiff
            d_clrurad_dt(jl,lev) = 0.0_wp
            dtotuclfl_dt(jl,lev) = dtotuclfl_dt(jl,lev) + duclfl_dt * delwave(iband) * fluxfac
          ENDDO
        ENDDO
      ENDIF

      ! End spectral band loop
    ENDDO

    ! Calculate fluxes at surface
    DO jl = 1, kproma  ! loop over columns
      totuflux(jl,0) = totuflux(jl,0) * fluxfac
      totdflux(jl,0) = totdflux(jl,0) * fluxfac
      fnet(jl,0) = totuflux(jl,0) - totdflux(jl,0)

      totuclfl(jl,0) = totuclfl(jl,0) * fluxfac
      totdclfl(jl,0) = totdclfl(jl,0) * fluxfac
      fnetc(jl,0) = totuclfl(jl,0) - totdclfl(jl,0)
    ENDDO

    ! Calculate fluxes at model levels
    DO lev = 1, nlayers
      DO jl = 1, kproma  ! loop over columns
        totuflux(jl,lev) = totuflux(jl,lev) * fluxfac
        totdflux(jl,lev) = totdflux(jl,lev) * fluxfac
        fnet(jl,lev) = totuflux(jl,lev) - totdflux(jl,lev)
        totuclfl(jl,lev) = totuclfl(jl,lev) * fluxfac
        totdclfl(jl,lev) = totdclfl(jl,lev) * fluxfac
        fnetc(jl,lev) = totuclfl(jl,lev) - totdclfl(jl,lev)
      ENDDO
    ENDDO

    ! end vectorized version
#ifdef __TEST_lrtm_rtrnmr
  END SUBROUTINE lrtm_rtrnmr_1
#else
  END SUBROUTINE lrtm_rtrnmr
#endif


#ifdef __TEST_lrtm_rtrnmr
    !-----------------------------------------------------------------------------
  SUBROUTINE lrtm_rtrnmr(kproma, nlayers, istart, iend, iout, semiss, ncbands, &
    & cldfrac, taucloud, planklay, planklev, plankbnd, &
    & pwvcm, fracs, taut, &
    & totuflux, totdflux, fnet, &
    & totuclfl, totdclfl, fnetc, &
    & idrv, dplankbnd_dt, dtotuflux_dt, dtotuclfl_dt )

    ! ------- Declarations -------

    ! ----- Input -----
    INTEGER, INTENT(in) :: kproma          ! number of columns
    INTEGER, INTENT(in) :: nlayers         ! total number of layers
    INTEGER, INTENT(in) :: istart          ! beginning band of calculation
    INTEGER, INTENT(in) :: iend            ! ending band of calculation
    INTEGER, INTENT(in) :: iout            ! output option flag

    ! Atmosphere
    !    Dimensions: (0:nlayers)
    REAL(wp), INTENT(in) :: pwvcm(:)           ! precipitable water vapor (cm)
    REAL(wp), INTENT(in) :: semiss(:,:)        ! lw surface emissivity
    !    Dimensions: (nbndlw)
    REAL(wp), INTENT(in) :: planklay(:,:,:)      !
    !    Dimensions: (nlayers,nbndlw)
    REAL(wp), INTENT(in) :: planklev(:,0:,:)     !
    !    Dimensions: (0:nlayers,nbndlw)
    REAL(wp), INTENT(in) :: plankbnd(:,:)        !
    !    Dimensions: (nbndlw)
    REAL(wp), INTENT(in) :: fracs(:,:,:)         !
    !    Dimensions: (kproma,nlayers,ngptw)
    REAL(wp), INTENT(in) :: taut(:,:,:)          ! gaseous + aerosol optical depths
    !    Dimensions: (kproma,nlayers,ngptlw)

    ! Clouds
    INTEGER, INTENT(in) :: ncbands(:)          ! number of cloud spectral bands
    REAL(wp), INTENT(in) :: cldfrac(:,:)       ! layer cloud fraction
    !    Dimensions: (nlayers)
    REAL(wp), INTENT(in) :: taucloud(:,:,:)      ! layer cloud optical depth
    !    Dimensions: (nlayers,nbndlw)
    INTEGER, INTENT(in) :: idrv            ! flag for calculation of dF/dt from
    ! Planck derivative [0=off, 1=on]
    REAL(wp), INTENT(in) :: dplankbnd_dt(:,:)    ! derivative of Planck function wrt temp
    !    Dimensions: (nbndlw)

    ! ----- Output -----
    REAL(wp), INTENT(out) :: totuflux(:,0:)      ! upward longwave flux (w/m2)
    !    Dimensions: (0:nlayers)
    REAL(wp), INTENT(out) :: totdflux(:,0:)      ! downward longwave flux (w/m2)
    !    Dimensions: (0:nlayers)
    REAL(wp), INTENT(out) :: fnet(:,0:)          ! net longwave flux (w/m2)
    !    Dimensions: (0:nlayers)
    REAL(wp), INTENT(out) :: totuclfl(:,0:)      ! clear sky upward longwave flux (w/m2)
    !    Dimensions: (0:nlayers)
    REAL(wp), INTENT(out) :: totdclfl(:,0:)      ! clear sky downward longwave flux (w/m2)
    !    Dimensions: (0:nlayers)
    REAL(wp), INTENT(out) :: fnetc(:,0:)         ! clear sky net longwave flux (w/m2)
    !    Dimensions: (0:nlayers)
    REAL(wp), INTENT(out) :: dtotuflux_dt(:,0:)  ! change in upward longwave flux (w/m2/k)
    ! with respect to surface temperature
    !    Dimensions: (0:nlayers)
    REAL(wp), INTENT(out) :: dtotuclfl_dt(:,0:)  ! change in upward longwave flux (w/m2/k)
    ! with respect to surface temperature
    !    Dimensions: (0:nlayers)


    ! tmp variables for testing
    REAL(wp), POINTER :: tmp_totuflux(:,:)      ! upward longwave flux (w/m2)
    !    Dimensions: (0:nlayers)
    REAL(wp), POINTER:: tmp_totdflux(:,:)      ! downward longwave flux (w/m2)
    !    Dimensions: (0:nlayers)
    REAL(wp), POINTER:: tmp_fnet(:,:)          ! net longwave flux (w/m2)
    !    Dimensions: (0:nlayers)
    REAL(wp), POINTER:: tmp_totuclfl(:,:)      ! clear sky upward longwave flux (w/m2)
    !    Dimensions: (0:nlayers)
    REAL(wp), POINTER:: tmp_totdclfl(:,:)      ! clear sky downward longwave flux (w/m2)
    !    Dimensions: (0:nlayers)
    REAL(wp), POINTER:: tmp_fnetc(:,:)         ! clear sky net longwave flux (w/m2)
    !    Dimensions: (0:nlayers)
    REAL(wp), POINTER:: tmp_dtotuflux_dt(:,:)  ! change in upward longwave flux (w/m2/k)
    ! with respect to surface temperature
    !    Dimensions: (0:nlayers)
    REAL(wp), POINTER:: tmp_dtotuclfl_dt(:,:)  ! change in upward longwave flux (w/m2/k)
    ! with respect to surface temperature

    ALLOCATE(tmp_totuflux(lbound(totuflux,1):ubound(totuflux,1), &
                          lbound(totuflux,2):ubound(totuflux,2)))
    ALLOCATE(tmp_totdflux(lbound(totdflux,1):ubound(totdflux,1), &
                          lbound(totdflux,2):ubound(totdflux,2)))
    ALLOCATE(tmp_fnet(lbound(fnet,1):ubound(fnet,1), &
                      lbound(fnet,2):ubound(fnet,2)))
    ALLOCATE(tmp_totuclfl(lbound(totuclfl,1):ubound(totuclfl,1), &
                          lbound(totuclfl,2):ubound(totuclfl,2)))
    ALLOCATE(tmp_totdclfl(lbound(totdclfl,1):ubound(totdclfl,1), &
                          lbound(totdclfl,2):ubound(totdclfl,2)))
    ALLOCATE(tmp_fnetc(lbound(fnetc,1):ubound(fnetc,1), &
                       lbound(fnetc,2):ubound(fnetc,2)))
    ALLOCATE(tmp_dtotuflux_dt(lbound(dtotuflux_dt,1):ubound(dtotuflux_dt,1), &
                              lbound(dtotuflux_dt,2):ubound(dtotuflux_dt,2)))
    ALLOCATE(tmp_dtotuclfl_dt(lbound(dtotuclfl_dt,1):ubound(dtotuclfl_dt,1), &
                              lbound(dtotuclfl_dt,2):ubound(dtotuclfl_dt,2)))

    tmp_totuflux(:,:) = 0._wp
    tmp_totdflux(:,:) = 0._wp
    tmp_fnet(:,:) = 0._wp
    tmp_totuclfl(:,:) = 0._wp
    tmp_totdclfl(:,:) = 0._wp
    tmp_fnetc(:,:) = 0._wp
    tmp_dtotuflux_dt(:,:) = 0._wp
    tmp_dtotuclfl_dt(:,:) = 0._wp

    totuflux(:,:) = 0._wp
    totdflux(:,:) = 0._wp
    fnet(:,:) = 0._wp
    totuclfl(:,:) = 0._wp
    totdclfl(:,:) = 0._wp
    fnetc(:,:) = 0._wp
    dtotuflux_dt(:,:) = 0._wp
    dtotuclfl_dt(:,:) = 0._wp


    CALL timer_start(timer_lrtm_2)
    CALL lrtm_rtrnmr_1(kproma, nlayers, istart, iend, iout, semiss, ncbands, &
    & cldfrac, taucloud, planklay, planklev, plankbnd, &
    & pwvcm, fracs, taut, &
    & totuflux, totdflux, fnet, &
    & totuclfl, totdclfl, fnetc, &
    & idrv, dplankbnd_dt, dtotuflux_dt, dtotuclfl_dt )
    CALL timer_stop(timer_lrtm_2)

    CALL timer_start(timer_lrtm_1)
    CALL lrtm_rtrnmr_2(kproma, nlayers, istart, iend, iout, semiss, ncbands, &
    & cldfrac, taucloud, planklay, planklev, plankbnd, &
    & pwvcm, fracs, taut, &
    & tmp_totuflux, tmp_totdflux, tmp_fnet, &
    & tmp_totuclfl, tmp_totdclfl, tmp_fnetc, &
    & idrv, dplankbnd_dt, tmp_dtotuflux_dt, tmp_dtotuclfl_dt )
    CALL timer_stop(timer_lrtm_1)


    IF (MAXVAL(ABS(tmp_totuflux-totuflux)) > 0.0_wp) &
      CALL finish('lrtm_rtrnmr','Results not identical for totuflux')
    IF (MAXVAL(ABS(tmp_totdflux-totdflux)) > 0.0_wp) &
      CALL finish('lrtm_rtrnmr','Results not identical for totdflux')
    IF (MAXVAL(ABS(tmp_fnet-fnet)) > 0.0_wp) &
      CALL finish('lrtm_rtrnmr','Results not identical for fnet')
    IF (MAXVAL(ABS(tmp_totuclfl-totuclfl)) > 0.0_wp) &
      CALL finish('lrtm_rtrnmr','Results not identical for totuclfl')
    IF (MAXVAL(ABS(tmp_totdclfl-totdclfl)) > 0.0_wp) &
      CALL finish('lrtm_rtrnmr','Results not identical for totdclfl')
    IF (MAXVAL(ABS(tmp_fnetc-fnetc)) > 0.0_wp) &
      CALL finish('lrtm_rtrnmr','Results not identical for fnetc')
    IF (MAXVAL(ABS(tmp_dtotuflux_dt-dtotuflux_dt)) > 0.0_wp) &
      CALL finish('lrtm_rtrnmr','Results not identical for dtotuflux_dt')
    IF (MAXVAL(ABS(tmp_dtotuclfl_dt-dtotuclfl_dt)) > 0.0_wp) &
      CALL finish('lrtm_rtrnmr','Results not identical for dtotuclfl_dt')

    DEALLOCATE(tmp_totuflux)
    DEALLOCATE(tmp_totdflux)
    DEALLOCATE(tmp_fnet)
    DEALLOCATE(tmp_totuclfl)
    DEALLOCATE(tmp_totdclfl)
    DEALLOCATE(tmp_fnetc)
    DEALLOCATE(tmp_dtotuflux_dt)
    DEALLOCATE(tmp_dtotuclfl_dt)
  END SUBROUTINE lrtm_rtrnmr
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  ! Optimized lrtm_rtrnmr with extensive gather/scatter
  ! Leonidas Linardakis, MPI-M, Oct 2010
  ! Under development !!!
  !-----------------------------------------------------------------------------
  SUBROUTINE lrtm_rtrnmr_2(kproma, nlayers, istart, iend, iout, semiss, ncbands, &
    & cldfrac, taucloud, planklay, planklev, plankbnd, &
    & pwvcm, fracs, taut, &
    & totuflux, totdflux, fnet, &
    & totuclfl, totdclfl, fnetc, &
    & idrv, dplankbnd_dt, dtotuflux_dt, dtotuclfl_dt )
    !-----------------------------------------------------------------------------
    !
    !  Original version:   E. J. Mlawer, et al. RRTM_V3.0
    !  Revision for GCMs:  Michael J. Iacono; October, 2002
    !  Revision for F90:  Michael J. Iacono; June, 2006
    !  Revision for dFdT option: M. J. Iacono and E. J. Mlawer, November 2009
    !
    !  This program calculates the upward fluxes, downward fluxes, and
    !  heating rates for an arbitrary clear or cloudy atmosphere.  The input
    !  to this program is the atmospheric profile, all Planck function
    !  information, and the cloud fraction by layer.  A variable diffusivity
    !  angle (SECDIFF) is used for the angle integration.  Bands 2-3 and 5-9
    !  use a value for SECDIFF that varies from 1.50 to 1.80 as a function of
    !  the column water vapor, and other bands use a value of 1.66.  The Gaussian
    !  weight appropriate to this angle (WTDIFF=0.5) is applied here.  Note that
    !  use of the emissivity angle for the flux integration can cause errors of
    !  1 to 4 W/m2 within cloudy layers.
    !  Clouds are treated with a maximum-random cloud overlap method.
    !  This subroutine also provides the optional capability to calculate
    !  the derivative of upward flux respect to surface temperature using
    !  the pre-tabulated derivative of the Planck function with respect to
    !  temperature integrated over each spectral band.
    !***************************************************************************

    ! ------- Declarations -------

    ! ----- Input -----
    INTEGER, INTENT(in) :: kproma          ! number of columns
    INTEGER, INTENT(in) :: nlayers         ! total number of layers
    INTEGER, INTENT(in) :: istart          ! beginning band of calculation
    INTEGER, INTENT(in) :: iend            ! ending band of calculation
    INTEGER, INTENT(in) :: iout            ! output option flag

    ! Atmosphere
    !    Dimensions: (0:nlayers)
    REAL(wp), INTENT(in) :: pwvcm(:)           ! precipitable water vapor (cm)
    REAL(wp), INTENT(in) :: semiss(:,:)        ! lw surface emissivity
    !    Dimensions: (nbndlw)
    REAL(wp), INTENT(in) :: planklay(:,:,:)      !
    !    Dimensions: (nlayers,nbndlw)
    REAL(wp), INTENT(in) :: planklev(:,0:,:)     !
    !    Dimensions: (0:nlayers,nbndlw)
    REAL(wp), INTENT(in) :: plankbnd(:,:)        !
    !    Dimensions: (nbndlw)
    REAL(wp), INTENT(in) :: fracs(:,:,:)         !
    !    Dimensions: (kproma,nlayers,ngptw)
    REAL(wp), INTENT(in) :: taut(:,:,:)          ! gaseous + aerosol optical depths
    !    Dimensions: (kproma,nlayers,ngptlw)

    ! Clouds
    INTEGER, INTENT(in) :: ncbands(:)          ! number of cloud spectral bands
    REAL(wp), INTENT(in) :: cldfrac(:,:)       ! layer cloud fraction
    !    Dimensions: (nlayers)
    REAL(wp), INTENT(in) :: taucloud(:,:,:)      ! layer cloud optical depth
    !    Dimensions: (nlayers,nbndlw)
    INTEGER, INTENT(in) :: idrv            ! flag for calculation of dF/dt from
    ! Planck derivative [0=off, 1=on]
    REAL(wp), INTENT(in) :: dplankbnd_dt(:,:)    ! derivative of Planck function wrt temp
    !    Dimensions: (nbndlw)

    ! ----- Output -----
    REAL(wp), INTENT(out) :: totuflux(:,0:)      ! upward longwave flux (w/m2)
    !    Dimensions: (0:nlayers)
    REAL(wp), INTENT(out) :: totdflux(:,0:)      ! downward longwave flux (w/m2)
    !    Dimensions: (0:nlayers)
    REAL(wp), INTENT(out) :: fnet(:,0:)          ! net longwave flux (w/m2)
    !    Dimensions: (0:nlayers)
    REAL(wp), INTENT(out) :: totuclfl(:,0:)      ! clear sky upward longwave flux (w/m2)
    !    Dimensions: (0:nlayers)
    REAL(wp), INTENT(out) :: totdclfl(:,0:)      ! clear sky downward longwave flux (w/m2)
    !    Dimensions: (0:nlayers)
    REAL(wp), INTENT(out) :: fnetc(:,0:)         ! clear sky net longwave flux (w/m2)
    !    Dimensions: (0:nlayers)
    REAL(wp), INTENT(out) :: dtotuflux_dt(:,0:)  ! change in upward longwave flux (w/m2/k)
    ! with respect to surface temperature
    !    Dimensions: (0:nlayers)
    REAL(wp), INTENT(out) :: dtotuclfl_dt(:,0:)  ! change in upward longwave flux (w/m2/k)
    ! with respect to surface temperature
    !    Dimensions: (0:nlayers)

    ! Vectorized version implemented by Guenther Zaengl, DWD
    ! See above for variable descriptions
    ! ----- Local -----
    ! Declarations for radiative transfer
    REAL(wp) :: atot(kproma,nlayers)
    REAL(wp) :: atrans(kproma,nlayers)
    REAL(wp) :: bbugas(kproma,nlayers)
    REAL(wp) :: bbutot(kproma,nlayers)
    REAL(wp) :: clrurad(kproma,0:nlayers)
    REAL(wp) :: clrdrad(kproma,0:nlayers)
    REAL(wp) :: uflux
    REAL(wp) :: dflux
    REAL(wp) :: urad(kproma,0:nlayers)
    REAL(wp) :: drad(kproma,0:nlayers)
    REAL(wp) :: uclfl
    REAL(wp) :: dclfl

    REAL(wp) :: secdiff(kproma,nbndlw)          ! secant of diffusivity angle
    REAL(wp) :: a0(nbndlw),a1(nbndlw),a2(nbndlw)! diffusivity angle adjustment coefficients
    REAL(wp) :: wtdiff, rec_6
    REAL(wp) :: radld(kproma), radclrd(kproma), plfrac, blay, dplankup, dplankdn
    REAL(wp) :: odepth(kproma), odtot(kproma), odepth_rec, odtot_rec, gassrc(kproma), ttot
    REAL(wp) :: tblind, tfactot, bbd(kproma), bbdtot(kproma), tfacgas, transc, tausfac
    REAL(wp) :: rad0, reflect, radlu(kproma), radclru(kproma)

    REAL(wp) :: duflux_dt
    REAL(wp) :: duclfl_dt
    REAL(wp) :: d_urad_dt(kproma,0:nlayers)
    REAL(wp) :: d_clrurad_dt(kproma,0:nlayers)
    REAL(wp) :: d_rad0_dt, d_radlu_dt(kproma), d_radclru_dt(kproma)

    INTEGER :: icldlyr(kproma,nlayers)             ! flag for cloud in layer
    INTEGER :: ibnd, ib, iband, lay, lev    ! loop indices
    INTEGER :: igc                          ! g-point interval counter
    INTEGER :: iclddn(kproma)               ! flag for cloud in down path
    INTEGER :: ittot, itgas, itr            ! lookup table indices
    INTEGER :: ipat(16,0:2)
    INTEGER :: jl

    INTEGER :: ibv(kproma)

    ! Declarations for cloud overlap adjustment
    REAL(wp) :: faccld1(kproma,nlayers+1),faccld2(kproma,nlayers+1)
    REAL(wp) :: facclr1(kproma,nlayers+1),facclr2(kproma,nlayers+1)
    REAL(wp) :: faccmb1(kproma,nlayers+1),faccmb2(kproma,nlayers+1)
    REAL(wp) :: faccld1d(kproma,0:nlayers),faccld2d(kproma,0:nlayers)
    REAL(wp) :: facclr1d(kproma,0:nlayers),facclr2d(kproma,0:nlayers)
    REAL(wp) :: faccmb1d(kproma,0:nlayers),faccmb2d(kproma,0:nlayers)

    REAL(wp) :: fmax, fmin, rat1(kproma), rat2(kproma)
    REAL(wp), DIMENSION(kproma) :: clrradd, cldradd, clrradu, cldradu, oldclr, oldcld, &
      & rad, cldsrc, radmod



    INTEGER :: istcld(kproma,nlayers+1),istcldd(kproma,0:nlayers)

    ! Local variables for cloud / no cloud index lists
    INTEGER, DIMENSION(kproma,nlayers) :: icld_ind,iclear_ind
    INTEGER :: icld, iclear, n_cloudpoints(nlayers), n_clearpoints(nlayers)

    INTEGER, DIMENSION(kproma,nlayers) :: clear_to_cloud_idx, cloud_to_cloud_idx
    INTEGER :: could_cloud, clear_cloud, clear_to_cloud_cnt(nlayers), cloud_to_cloud_cnt(nlayers)
    REAL(wp) :: cldfrac_diff, cldfrac_diff_plus, cldfrac_diff_min, r_diff_ge_0
    LOGICAL :: diff_ge_0, diff_gt_0, diff_lt_0

    ! These arrays indicate the spectral 'region' (used in the
    ! calculation of ice cloud optical depths) corresponding
    ! to each spectral band.  See cldprop.f for more details.
    DATA ipat /1,1,1,1,1,1,1,1,1, 1, 1, 1, 1, 1, 1, 1, &
      & 1,2,3,3,3,4,4,4,5, 5, 5, 5, 5, 5, 5, 5, &
      & 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16/

    ! This secant and weight corresponds to the standard diffusivity
    ! angle.  This initial value is redefined below for some bands.
    DATA wtdiff /0.5_wp/
    DATA rec_6 /0.166667_wp/

    ! Reset diffusivity angle for Bands 2-3 and 5-9 to vary (between 1.50
    ! and 1.80) as a function of total column water vapor.  The function
    ! has been defined to minimize flux and cooling rate errors in these bands
    ! over a wide range of precipitable water values.
    DATA a0 / 1.66_wp,  1.55_wp,  1.58_wp,  1.66_wp, &
      & 1.54_wp, 1.454_wp,  1.89_wp,  1.33_wp, &
      & 1.668_wp,  1.66_wp,  1.66_wp,  1.66_wp, &
      & 1.66_wp,  1.66_wp,  1.66_wp,  1.66_wp /
    DATA a1 / 0.00_wp,  0.25_wp,  0.22_wp,  0.00_wp, &
      & 0.13_wp, 0.446_wp, -0.10_wp,  0.40_wp, &
      & -0.006_wp,  0.00_wp,  0.00_wp,  0.00_wp, &
      & 0.00_wp,  0.00_wp,  0.00_wp,  0.00_wp /
    DATA a2 / 0.00_wp, -12.0_wp, -11.7_wp,  0.00_wp, &
      & -0.72_wp,-0.243_wp,  0.19_wp,-0.062_wp, &
      & 0.414_wp,  0.00_wp,  0.00_wp,  0.00_wp, &
      & 0.00_wp,  0.00_wp,  0.00_wp,  0.00_wp /

    DO ibnd = 1,nbndlw
      IF (ibnd.EQ.1 .OR. ibnd.EQ.4 .OR. ibnd.GE.10) THEN
        DO jl = 1, kproma  ! loop over columns
          secdiff(jl,ibnd) = 1.66_wp
        ENDDO
      ELSE
        DO jl = 1, kproma  ! loop over columns
          secdiff(jl,ibnd) = a0(ibnd) + a1(ibnd)*EXP(a2(ibnd)*pwvcm(jl))
          IF (secdiff(jl,ibnd) .GT. 1.80_wp) secdiff(jl,ibnd) = 1.80_wp
          IF (secdiff(jl,ibnd) .LT. 1.50_wp) secdiff(jl,ibnd) = 1.50_wp
        ENDDO
      ENDIF
    ENDDO

    faccld1(:,:)  = 0.0_wp
    faccld2(:,:)  = 0.0_wp
    facclr1(:,:)  = 0.0_wp
    facclr2(:,:)  = 0.0_wp
    faccmb1(:,:)  = 0.0_wp
    faccmb2(:,:)  = 0.0_wp
    faccld1d(:,:) = 0.0_wp
    faccld2d(:,:) = 0.0_wp
    facclr1d(:,:) = 0.0_wp
    facclr2d(:,:) = 0.0_wp
    faccmb1d(:,:) = 0.0_wp
    faccmb2d(:,:) = 0.0_wp
    urad(:,:)     = 0.0_wp
    drad(:,:)     = 0.0_wp
    totuflux(:,:) = 0.0_wp
    totdflux(:,:) = 0.0_wp
    clrurad(:,:)  = 0.0_wp
    clrdrad(:,:)  = 0.0_wp
    totuclfl(:,:) = 0.0_wp
    totdclfl(:,:) = 0.0_wp

    IF (idrv .EQ. 1) THEN
      d_urad_dt(:,:)    = 0.0_wp
      d_clrurad_dt(:,:) = 0.0_wp
      dtotuflux_dt(:,:) = 0.0_wp
      dtotuclfl_dt(:,:) = 0.0_wp
    ENDIF

    !---------------------------------------------
    ! fill cloud, non-cloud points for each layer
    clear_to_cloud_cnt(:) = 0
    cloud_to_cloud_cnt(:) = 0
    DO lay = 1, nlayers

      icld   = 0
      iclear = 0

      DO jl = 1, kproma  ! loop over columns
        IF (cldfrac(jl,lay) .GE. 1.e-6_wp) THEN
          icldlyr(jl,lay) = 1
          icld = icld + 1
          icld_ind(icld,lay) = jl
        ELSE
          icldlyr(jl,lay) = 0
          iclear = iclear + 1
          iclear_ind(iclear,lay) = jl
        ENDIF
      ENDDO

      n_cloudpoints(lay) = icld
      n_clearpoints(lay) = iclear

    ENDDO
    !---------------------------------------------

    ! Maximum/Random cloud overlap parameter

    ! fill istcld
    istcld(:,1) = 1

    lev = 1
    DO iclear = 1, n_clearpoints(lev)
      jl = iclear_ind(iclear,lev)
      istcld(jl,lev+1) = 1
    ENDDO
    DO icld = 1, n_cloudpoints(lev)
      jl = icld_ind(icld,lev)
      istcld(jl,lev+1) = 0
      clear_to_cloud_idx(icld, lev) = jl
    ENDDO
    clear_to_cloud_cnt(lev) = n_cloudpoints(lev)
    cloud_to_cloud_cnt(lev) = 0

    DO lev = 2, nlayers
      DO iclear = 1, n_clearpoints(lev)
        jl = iclear_ind(iclear,lev)
        istcld(jl,lev+1) = 1
      ENDDO

      could_cloud = 0
      clear_cloud = 0
      DO icld = 1, n_cloudpoints(lev)
        jl = icld_ind(icld,lev)
        istcld(jl,lev+1) = 0
        IF (istcld(jl,lev) == 0) THEN
          could_cloud = could_cloud + 1
          cloud_to_cloud_idx(could_cloud, lev) = jl
        ELSE
          clear_cloud = clear_cloud + 1
          clear_to_cloud_idx(clear_cloud, lev) = jl
        ENDIF
      ENDDO
      clear_to_cloud_cnt(lev) = clear_cloud
      cloud_to_cloud_cnt(lev) = could_cloud
    ENDDO

    !---------------------------------------------
    DO lev = 1, nlayers-1
      IF (clear_to_cloud_cnt(lev) == kproma) THEN
        DO jl = 1, kproma
          cldfrac_diff = cldfrac(jl,lev+1) - cldfrac(jl,lev)

          diff_gt_0 =  cldfrac_diff .GT. 0._wp
          diff_lt_0 =  cldfrac_diff .lT. 0._wp

        !  cldfrac_diff_plus = MAX(cldfrac_diff, 0._wp)
        !  cldfrac_diff_min  = MAX(-cldfrac_diff, 0._wp)

          faccld2(jl,lev+1) = MERGE(&
            &  -cldfrac_diff/cldfrac(jl,lev), &
            &  0._wp,  &
            &  diff_lt_0)
          facclr2(jl,lev+1) = MERGE(&
                  & cldfrac_diff/(1._wp-cldfrac(jl,lev)), &
                  & 0._wp, &
  !                 & cldfrac(jl,lev) .LT. 1._wp .AND. diff_gt_0)
                  diff_gt_0)

          rat2(jl) = MERGE(1._wp, 0._wp, faccld2(jl,lev+1) > 0._wp)
          rat1(jl) = MERGE(1._wp, 0._wp, facclr2(jl,lev+1) > 0._wp)
        ENDDO
      ELSE
        DO clear_cloud = 1, clear_to_cloud_cnt(lev)
          jl = clear_to_cloud_idx(clear_cloud, lev)

          cldfrac_diff = cldfrac(jl,lev+1) - cldfrac(jl,lev)

          diff_gt_0 =  cldfrac_diff .GT. 0._wp
          diff_lt_0 =  cldfrac_diff .lT. 0._wp

         ! cldfrac_diff_plus = MAX(cldfrac_diff, 0._wp)
         ! cldfrac_diff_min  = MAX(-cldfrac_diff, 0._wp)

          faccld2(jl,lev+1) = MERGE(&
            &  -cldfrac_diff/cldfrac(jl,lev), &
            &  0._wp,  &
            &  diff_lt_0)

          facclr2(jl,lev+1) = MERGE(&
                  & cldfrac_diff/(1._wp-cldfrac(jl,lev)), &
                  & 0._wp, &
  !                 & cldfrac(jl,lev) .LT. 1._wp .AND. diff_gt_0)
                  diff_gt_0)

          rat2(jl) = MERGE(1._wp, 0._wp, faccld2(jl,lev+1) > 0._wp)
          rat1(jl) = MERGE(1._wp, 0._wp, facclr2(jl,lev+1) > 0._wp)

        ENDDO
      ENDIF
      !----------------------------------------

      !----------------------------------------
      IF (cloud_to_cloud_cnt(lev) == kproma) THEN
        DO jl=1,kproma

          cldfrac_diff = cldfrac(jl,lev+1) - cldfrac(jl,lev)
          diff_ge_0 =  cldfrac_diff .GE. 0._wp

          fmax = MAX(cldfrac(jl,lev),cldfrac(jl,lev-1))
          fmin = MIN(cldfrac(jl,lev),cldfrac(jl,lev-1))
          r_diff_ge_0 = MERGE(1._wp, 0._wp, diff_ge_0)

  !        write(0,*) 'facclr2...'
          facclr2(jl,lev+1) = MERGE(&
              & (cldfrac(jl,lev+1)-fmax)/(1._wp-fmax),  &
              &  0._wp, &
              &  cldfrac(jl,lev+1) .GT. fmax .AND. diff_ge_0)

  !        write(0,*) 'facclr1...'
          facclr1(jl,lev+1) = MERGE(&
              & cldfrac_diff/(cldfrac(jl,lev-1)-cldfrac(jl,lev)),&
              & rat2(jl)*r_diff_ge_0, &
              & cldfrac(jl,lev+1) .LT. fmax .AND. diff_ge_0)

  !        write(0,*) 'faccld1...'
          faccld1(jl,lev+1) = MERGE(&
            & -cldfrac_diff/(cldfrac(jl,lev)-fmin), &
            & rat1(jl) * (1._wp-r_diff_ge_0),  &
            & cldfrac(jl,lev+1) .GT. fmin .AND. .NOT. diff_ge_0)

  !        write(0,*) 'faccld2...'
          faccld2(jl,lev+1) = MERGE(&
            & (fmin-cldfrac(jl,lev+1))/fmin, &
            & 0._wp, &
            & cldfrac(jl,lev+1) .LE. fmin .AND. .NOT. diff_ge_0)

          rat1(jl) = MERGE(1._wp, 0._wp, &
            & facclr1(jl,lev+1) > 0._wp .OR. facclr2(jl,lev+1) > 0._wp)
          rat2(jl) = MERGE(1._wp,  0._wp, &
            & faccld1(jl,lev+1) > 0._wp .OR. faccld2(jl,lev+1) > 0._wp)

          faccmb1(jl,lev+1) = facclr1(jl,lev+1) * faccld2(jl,lev) * cldfrac(jl,lev-1)
          faccmb2(jl,lev+1) = faccld1(jl,lev+1) * facclr2(jl,lev) * (1._wp - cldfrac(jl,lev-1))

        ENDDO
      ELSE
        DO could_cloud = 1, cloud_to_cloud_cnt(lev)
          jl = cloud_to_cloud_idx(could_cloud, lev)

          cldfrac_diff = cldfrac(jl,lev+1) - cldfrac(jl,lev)
          diff_ge_0 =  cldfrac_diff .GE. 0._wp

          fmax = MAX(cldfrac(jl,lev),cldfrac(jl,lev-1))
          fmin = MIN(cldfrac(jl,lev),cldfrac(jl,lev-1))
          r_diff_ge_0 = MERGE(1._wp, 0._wp, diff_ge_0)

  !        write(0,*) 'facclr2...'
          facclr2(jl,lev+1) = MERGE(&
              & (cldfrac(jl,lev+1)-fmax)/(1._wp-fmax),  &
              &  0._wp, &
              &  cldfrac(jl,lev+1) .GT. fmax .AND. diff_ge_0)

  !        write(0,*) 'facclr1...'
          facclr1(jl,lev+1) = MERGE(&
              & cldfrac_diff/(cldfrac(jl,lev-1)-cldfrac(jl,lev)),&
              & rat2(jl)*r_diff_ge_0, &
              & cldfrac(jl,lev+1) .LT. fmax .AND. diff_ge_0)

  !        write(0,*) 'faccld1...'
          faccld1(jl,lev+1) = MERGE(&
            & -cldfrac_diff/(cldfrac(jl,lev)-fmin), &
            & rat1(jl) * (1._wp-r_diff_ge_0),  &
            & cldfrac(jl,lev+1) .GT. fmin .AND. .NOT. diff_ge_0)

  !        write(0,*) 'faccld2...'
          faccld2(jl,lev+1) = MERGE(&
            & (fmin-cldfrac(jl,lev+1))/fmin, &
            & 0._wp, &
            & cldfrac(jl,lev+1) .LE. fmin .AND. .NOT. diff_ge_0)

          rat1(jl) = MERGE(1._wp, 0._wp, &
            & facclr1(jl,lev+1) > 0._wp .OR. facclr2(jl,lev+1) > 0._wp)
          rat2(jl) = MERGE(1._wp,  0._wp, &
            & faccld1(jl,lev+1) > 0._wp .OR. faccld2(jl,lev+1) > 0._wp)

          faccmb1(jl,lev+1) = facclr1(jl,lev+1) * faccld2(jl,lev) * cldfrac(jl,lev-1)
          faccmb2(jl,lev+1) = faccld1(jl,lev+1) * facclr2(jl,lev) * (1._wp - cldfrac(jl,lev-1))

        ENDDO
      ENDIF

    ENDDO

    !---------------------------
    ! process last layer
    ! this is already 0
!     lev = nlayers
!     DO icld = 1, n_cloudpoints(lev)
!       jl = icld_ind(icld,lev)
!       faccld1(jl,lev+1) = 0._wp
!       faccld2(jl,lev+1) = 0._wp
!       facclr1(jl,lev+1) = 0._wp
!       facclr2(jl,lev+1) = 0._wp
!       faccmb1(jl,lev+1) = 0._wp
!       faccmb2(jl,lev+1) = 0._wp
!       faccmb1(jl,lev+1) = 0._wp
!       faccmb2(jl,lev+1) = 0._wp
!     ENDDO
    !---------------------------

    istcldd(:,nlayers) = 1

    DO lev = nlayers, 1, -1

      IF (n_cloudpoints(lev) == kproma) THEN ! all points are cloudy
        DO jl = 1, kproma ! Thus, direct addressing can be used
          istcldd(jl,lev-1) = 0
          IF (lev .EQ. 1) THEN
            faccld1d(jl,lev-1) = 0._wp
            faccld2d(jl,lev-1) = 0._wp
            facclr1d(jl,lev-1) = 0._wp
            facclr2d(jl,lev-1) = 0._wp
            faccmb1d(jl,lev-1) = 0._wp
            faccmb2d(jl,lev-1) = 0._wp
          ELSEIF (cldfrac(jl,lev-1) .GE. cldfrac(jl,lev)) THEN
            faccld1d(jl,lev-1) = 0._wp
            faccld2d(jl,lev-1) = 0._wp
            IF (istcldd(jl,lev) .EQ. 1) THEN
              facclr1d(jl,lev-1) = 0._wp
              facclr2d(jl,lev-1) = 0._wp
              IF (cldfrac(jl,lev) .LT. 1._wp) facclr2d(jl,lev-1) = &
                & (cldfrac(jl,lev-1)-cldfrac(jl,lev))/(1._wp-cldfrac(jl,lev))
              facclr2d(jl,lev) = 0._wp
              faccld2d(jl,lev) = 0._wp
            ELSE
              fmax = MAX(cldfrac(jl,lev),cldfrac(jl,lev+1))
              IF (cldfrac(jl,lev-1) .GT. fmax) THEN
                facclr1d(jl,lev-1) = rat2(jl)
                facclr2d(jl,lev-1) = (cldfrac(jl,lev-1)-fmax)/(1._wp-fmax)
              ELSEIF (cldfrac(jl,lev-1) .LT. fmax) THEN
                facclr1d(jl,lev-1) = (cldfrac(jl,lev-1)-cldfrac(jl,lev))/ &
                  & (cldfrac(jl,lev+1)-cldfrac(jl,lev))
                facclr2d(jl,lev-1) = 0._wp
              ELSE
                facclr1d(jl,lev-1) = rat2(jl)
                facclr2d(jl,lev-1) = 0._wp
              ENDIF
            ENDIF
            IF (facclr1d(jl,lev-1).GT.0._wp .OR. facclr2d(jl,lev-1).GT.0._wp)THEN
              rat1(jl) = 1._wp
              rat2(jl) = 0._wp
            ELSE
              rat1(jl) = 0._wp
              rat2(jl) = 0._wp
            ENDIF
          ELSE
            facclr1d(jl,lev-1) = 0._wp
            facclr2d(jl,lev-1) = 0._wp
            IF (istcldd(jl,lev) .EQ. 1) THEN
              faccld1d(jl,lev-1) = 0._wp
              faccld2d(jl,lev-1) = (cldfrac(jl,lev)-cldfrac(jl,lev-1))/cldfrac(jl,lev)
              facclr2d(jl,lev) = 0._wp
              faccld2d(jl,lev) = 0._wp
            ELSE
              fmin = MIN(cldfrac(jl,lev),cldfrac(jl,lev+1))
              IF (cldfrac(jl,lev-1) .LE. fmin) THEN
                faccld1d(jl,lev-1) = rat1(jl)
                faccld2d(jl,lev-1) = (fmin-cldfrac(jl,lev-1))/fmin
              ELSE
                faccld1d(jl,lev-1) = (cldfrac(jl,lev)-cldfrac(jl,lev-1))/(cldfrac(jl,lev)-fmin)
                faccld2d(jl,lev-1) = 0._wp
              ENDIF
            ENDIF
            IF (faccld1d(jl,lev-1).GT.0._wp .OR. faccld2d(jl,lev-1).GT.0._wp)THEN
              rat1(jl) = 0._wp
              rat2(jl) = 1._wp
            ELSE
              rat1(jl) = 0._wp
              rat2(jl) = 0._wp
            ENDIF
          ENDIF
          IF (lev == nlayers) THEN
            faccmb1d(jl,lev-1) = 0._wp
            faccmb2d(jl,lev-1) = faccld1d(jl,lev-1) * facclr2d(jl,lev)
          ELSE
            faccmb1d(jl,lev-1) = facclr1d(jl,lev-1) * faccld2d(jl,lev) * cldfrac(jl,lev+1)
            faccmb2d(jl,lev-1) = faccld1d(jl,lev-1) * facclr2d(jl,lev) * (1._wp-cldfrac(jl,lev+1))
          ENDIF
        ENDDO
      ELSE IF (n_clearpoints(lev) == kproma) THEN ! all points are clear
        istcldd(1:kproma,lev-1) = 1
      ELSE  ! use index list for the case that not all points are cloudy
        DO icld = 1, n_cloudpoints(lev)
          jl = icld_ind(icld,lev)

          istcldd(jl,lev-1) = 0
          IF (lev .EQ. 1) THEN
            faccld1d(jl,lev-1) = 0._wp
            faccld2d(jl,lev-1) = 0._wp
            facclr1d(jl,lev-1) = 0._wp
            facclr2d(jl,lev-1) = 0._wp
            faccmb1d(jl,lev-1) = 0._wp
            faccmb2d(jl,lev-1) = 0._wp
          ELSEIF (cldfrac(jl,lev-1) .GE. cldfrac(jl,lev)) THEN
            faccld1d(jl,lev-1) = 0._wp
            faccld2d(jl,lev-1) = 0._wp
            IF (istcldd(jl,lev) .EQ. 1) THEN
              facclr1d(jl,lev-1) = 0._wp
              facclr2d(jl,lev-1) = 0._wp
              IF (cldfrac(jl,lev) .LT. 1._wp) facclr2d(jl,lev-1) = &
                & (cldfrac(jl,lev-1)-cldfrac(jl,lev))/(1._wp-cldfrac(jl,lev))
              facclr2d(jl,lev) = 0._wp
              faccld2d(jl,lev) = 0._wp
            ELSE
              fmax = MAX(cldfrac(jl,lev),cldfrac(jl,lev+1))
              IF (cldfrac(jl,lev-1) .GT. fmax) THEN
                facclr1d(jl,lev-1) = rat2(jl)
                facclr2d(jl,lev-1) = (cldfrac(jl,lev-1)-fmax)/(1._wp-fmax)
              ELSEIF (cldfrac(jl,lev-1) .LT. fmax) THEN
                facclr1d(jl,lev-1) = (cldfrac(jl,lev-1)-cldfrac(jl,lev))/ &
                  & (cldfrac(jl,lev+1)-cldfrac(jl,lev))
                facclr2d(jl,lev-1) = 0._wp
              ELSE
                facclr1d(jl,lev-1) = rat2(jl)
                facclr2d(jl,lev-1) = 0._wp
              ENDIF
            ENDIF
            IF (facclr1d(jl,lev-1).GT.0._wp .OR. facclr2d(jl,lev-1).GT.0._wp)THEN
              rat1(jl) = 1._wp
              rat2(jl) = 0._wp
            ELSE
              rat1(jl) = 0._wp
              rat2(jl) = 0._wp
            ENDIF
          ELSE
            facclr1d(jl,lev-1) = 0._wp
            facclr2d(jl,lev-1) = 0._wp
            IF (istcldd(jl,lev) .EQ. 1) THEN
              faccld1d(jl,lev-1) = 0._wp
              faccld2d(jl,lev-1) = (cldfrac(jl,lev)-cldfrac(jl,lev-1))/cldfrac(jl,lev)
              facclr2d(jl,lev) = 0._wp
              faccld2d(jl,lev) = 0._wp
            ELSE
              fmin = MIN(cldfrac(jl,lev),cldfrac(jl,lev+1))
              IF (cldfrac(jl,lev-1) .LE. fmin) THEN
                faccld1d(jl,lev-1) = rat1(jl)
                faccld2d(jl,lev-1) = (fmin-cldfrac(jl,lev-1))/fmin
              ELSE
                faccld1d(jl,lev-1) = (cldfrac(jl,lev)-cldfrac(jl,lev-1))/(cldfrac(jl,lev)-fmin)
                faccld2d(jl,lev-1) = 0._wp
              ENDIF
            ENDIF
            IF (faccld1d(jl,lev-1).GT.0._wp .OR. faccld2d(jl,lev-1).GT.0._wp)THEN
              rat1(jl) = 0._wp
              rat2(jl) = 1._wp
            ELSE
              rat1(jl) = 0._wp
              rat2(jl) = 0._wp
            ENDIF
          ENDIF
          IF (lev == nlayers) THEN
            faccmb1d(jl,lev-1) = 0._wp
            faccmb2d(jl,lev-1) = faccld1d(jl,lev-1) * facclr2d(jl,lev)
          ELSE
            faccmb1d(jl,lev-1) = facclr1d(jl,lev-1) * faccld2d(jl,lev) * cldfrac(jl,lev+1)
            faccmb2d(jl,lev-1) = faccld1d(jl,lev-1) * facclr2d(jl,lev) * (1._wp-cldfrac(jl,lev+1))
          ENDIF
        ENDDO
        DO iclear = 1, n_clearpoints(lev)
          jl = iclear_ind(iclear,lev)
          istcldd(jl,lev-1) = 1
        ENDDO
      ENDIF

    ENDDO

    igc = 1
    ! Loop over frequency bands.
    DO iband = istart, iend

      ! Reinitialize g-point counter for each band if output for each band is requested.
      IF (iout.GT.0.AND.iband.GE.2) igc = ngs(iband-1)+1

      ! Loop over g-channels.
1000    CONTINUE

      ! Radiative transfer starts here.
      radld(:) = 0._wp
      radclrd(:) = 0._wp
      iclddn(:) = 0

      DO jl = 1, kproma
        IF (ncbands(jl) .EQ. 1) THEN
          ibv(jl) = ipat(iband,0)
        ELSEIF (ncbands(jl) .EQ.  5) THEN
          ibv(jl) = ipat(iband,1)
        ELSEIF (ncbands(jl) .EQ. 16) THEN
          ibv(jl) = ipat(iband,2)
        ENDIF
      ENDDO

      ! Downward radiative transfer loop.
      DO lev = nlayers, 1, -1
        IF (n_cloudpoints(lev) == kproma) THEN ! all points are cloudy

          DO jl = 1, kproma ! Thus, direct addressing can be used

            ib = ibv(jl)

            plfrac = fracs(jl,lev,igc)
            blay = planklay(jl,lev,iband)
            dplankup = planklev(jl,lev,iband) - blay
            dplankdn = planklev(jl,lev-1,iband) - blay
            odepth(jl) = secdiff(jl,iband) * taut(jl,lev,igc)

            IF (odepth(jl) .LT. 0.0_wp) odepth(jl) = 0.0_wp

            iclddn(jl) = 1
            odtot(jl) = odepth(jl) + secdiff(jl,ib) * taucloud(jl,lev,ib)
            IF (odtot(jl) .LT. 0.06_wp) THEN
              atrans(jl,lev) = odepth(jl) - 0.5_wp*odepth(jl)*odepth(jl)
              odepth_rec = rec_6*odepth(jl)
              gassrc(jl) = plfrac*(blay+dplankdn*odepth_rec)*atrans(jl,lev)

              atot(jl,lev) = odtot(jl) - 0.5_wp*odtot(jl)*odtot(jl)
              odtot_rec = rec_6*odtot(jl)
              bbdtot(jl) =  plfrac * (blay+dplankdn*odtot_rec)
              bbd(jl) = plfrac*(blay+dplankdn*odepth_rec)

              bbugas(jl,lev) =  plfrac * (blay+dplankup*odepth_rec)
              bbutot(jl,lev) =  plfrac * (blay+dplankup*odtot_rec)
            ELSEIF (odepth(jl) .LE. 0.06_wp) THEN
              atrans(jl,lev) = odepth(jl) - 0.5_wp*odepth(jl)*odepth(jl)
              odepth_rec = rec_6*odepth(jl)
              gassrc(jl) = plfrac*(blay+dplankdn*odepth_rec)*atrans(jl,lev)

              odtot(jl) = odepth(jl) + secdiff(jl,ib) * taucloud(jl,lev,ib)
              tblind = odtot(jl)/(bpade+odtot(jl))
              ittot = INT(tblint*tblind + 0.5_wp)
              tfactot = tfn_tbl(ittot)
              bbdtot(jl) = plfrac * (blay + tfactot*dplankdn)
              bbd(jl) = plfrac*(blay+dplankdn*odepth_rec)
              atot(jl,lev) = 1._wp - exp_tbl(ittot)

              bbugas(jl,lev) = plfrac * (blay + dplankup*odepth_rec)
              bbutot(jl,lev) = plfrac * (blay + tfactot * dplankup)
            ELSE
              tblind = odepth(jl)/(bpade+odepth(jl))
              itgas = INT(tblint*tblind+0.5_wp)
              odepth(jl) = tau_tbl(itgas)
              atrans(jl,lev) = 1._wp - exp_tbl(itgas)
              tfacgas = tfn_tbl(itgas)
              gassrc(jl) = atrans(jl,lev) * plfrac * (blay + tfacgas*dplankdn)

              odtot(jl) = odepth(jl) + secdiff(jl,ib) * taucloud(jl,lev,ib)
              tblind = odtot(jl)/(bpade+odtot(jl))
              ittot = INT(tblint*tblind + 0.5_wp)
              tfactot = tfn_tbl(ittot)
              bbdtot(jl) = plfrac * (blay + tfactot*dplankdn)
              bbd(jl) = plfrac*(blay+tfacgas*dplankdn)
              atot(jl,lev) = 1._wp - exp_tbl(ittot)

              bbugas(jl,lev) = plfrac * (blay + tfacgas * dplankup)
              bbutot(jl,lev) = plfrac * (blay + tfactot * dplankup)
            ENDIF

            IF (istcldd(jl,lev) .EQ. 1) THEN
              cldradd(jl) = cldfrac(jl,lev) * radld(jl)
              clrradd(jl) = radld(jl) - cldradd(jl)
              oldcld(jl) = cldradd(jl)
              oldclr(jl) = clrradd(jl)
              rad(jl) = 0._wp
            ENDIF
            ttot = 1._wp - atot(jl,lev)
            cldsrc(jl) = bbdtot(jl) * atot(jl,lev)
            cldradd(jl) = cldradd(jl) * ttot + cldfrac(jl,lev) * cldsrc(jl)
            clrradd(jl) = clrradd(jl) * (1._wp-atrans(jl,lev)) + &
              & (1._wp-cldfrac(jl,lev))*gassrc(jl)
            radld(jl) = cldradd(jl) + clrradd(jl)
            drad(jl,lev-1) = drad(jl,lev-1) + radld(jl)

            radmod(jl) = rad(jl) * &
              & (facclr1d(jl,lev-1) * (1._wp-atrans(jl,lev)) + &
              & faccld1d(jl,lev-1) *  ttot) - &
              & faccmb1d(jl,lev-1) * gassrc(jl) + &
              & faccmb2d(jl,lev-1) * cldsrc(jl)

            oldcld(jl) = cldradd(jl) - radmod(jl)
            oldclr(jl) = clrradd(jl) + radmod(jl)
            rad(jl) = -radmod(jl) + facclr2d(jl,lev-1)*oldclr(jl) - faccld2d(jl,lev-1)*oldcld(jl)
            cldradd(jl) = cldradd(jl) + rad(jl)
            clrradd(jl) = clrradd(jl) - rad(jl)
          ENDDO

        ELSE IF (n_clearpoints(lev) == kproma) THEN ! all points are clear

          DO jl = 1, kproma ! Thus, direct addressing can be used

            ib = ibv(jl)

            plfrac = fracs(jl,lev,igc)
            blay = planklay(jl,lev,iband)
            dplankup = planklev(jl,lev,iband) - blay
            dplankdn = planklev(jl,lev-1,iband) - blay
            odepth(jl) = secdiff(jl,iband) * taut(jl,lev,igc)

            IF (odepth(jl) .LT. 0.0_wp) odepth(jl) = 0.0_wp

            IF (odepth(jl) .LE. 0.06_wp) THEN
              atrans(jl,lev) = odepth(jl)-0.5_wp*odepth(jl)*odepth(jl)
              odepth(jl) = rec_6*odepth(jl)
              bbd(jl) = plfrac*(blay+dplankdn*odepth(jl))
              bbugas(jl,lev) = plfrac*(blay+dplankup*odepth(jl))
            ELSE
              tblind = odepth(jl)/(bpade+odepth(jl))
              itr = INT(tblint*tblind+0.5_wp)
              transc = exp_tbl(itr)
              atrans(jl,lev) = 1._wp-transc
              tausfac = tfn_tbl(itr)
              bbd(jl) = plfrac*(blay+tausfac*dplankdn)
              bbugas(jl,lev) = plfrac * (blay + tausfac * dplankup)
            ENDIF
            radld(jl) = radld(jl) + (bbd(jl)-radld(jl))*atrans(jl,lev)
            drad(jl,lev-1) = drad(jl,lev-1) + radld(jl)
          ENDDO

        ELSE ! both cloudy and clear points are in the vector; use index list
          ! or if clause
!! CDIR NODEP
          ! Cloudy layer
          !          do icld = 1, n_cloudpoints(lev)
          !            jl = icld_ind(icld,lev)
          DO jl = 1, kproma
            IF (icldlyr(jl,lev) .EQ. 1) THEN
              ib = ibv(jl)

              plfrac = fracs(jl,lev,igc)
              blay = planklay(jl,lev,iband)
              dplankup = planklev(jl,lev,iband) - blay
              dplankdn = planklev(jl,lev-1,iband) - blay
              odepth(jl) = secdiff(jl,iband) * taut(jl,lev,igc)

              IF (odepth(jl) .LT. 0.0_wp) odepth(jl) = 0.0_wp

              iclddn(jl) = 1
              odtot(jl) = odepth(jl) + secdiff(jl,ib) * taucloud(jl,lev,ib)
              IF (odtot(jl) .LT. 0.06_wp) THEN
                atrans(jl,lev) = odepth(jl) - 0.5_wp*odepth(jl)*odepth(jl)
                odepth_rec = rec_6*odepth(jl)
                gassrc(jl) = plfrac*(blay+dplankdn*odepth_rec)*atrans(jl,lev)

                atot(jl,lev) = odtot(jl) - 0.5_wp*odtot(jl)*odtot(jl)
                odtot_rec = rec_6*odtot(jl)
                bbdtot(jl) =  plfrac * (blay+dplankdn*odtot_rec)
                bbd(jl) = plfrac*(blay+dplankdn*odepth_rec)

                bbugas(jl,lev) =  plfrac * (blay+dplankup*odepth_rec)
                bbutot(jl,lev) =  plfrac * (blay+dplankup*odtot_rec)
              ELSEIF (odepth(jl) .LE. 0.06_wp) THEN
                atrans(jl,lev) = odepth(jl) - 0.5_wp*odepth(jl)*odepth(jl)
                odepth_rec = rec_6*odepth(jl)
                gassrc(jl) = plfrac*(blay+dplankdn*odepth_rec)*atrans(jl,lev)

                odtot(jl) = odepth(jl) + secdiff(jl,ib) * taucloud(jl,lev,ib)
                tblind = odtot(jl)/(bpade+odtot(jl))
                ittot = INT(tblint*tblind + 0.5_wp)
                tfactot = tfn_tbl(ittot)
                bbdtot(jl) = plfrac * (blay + tfactot*dplankdn)
                bbd(jl) = plfrac*(blay+dplankdn*odepth_rec)
                atot(jl,lev) = 1._wp - exp_tbl(ittot)

                bbugas(jl,lev) = plfrac * (blay + dplankup*odepth_rec)
                bbutot(jl,lev) = plfrac * (blay + tfactot * dplankup)
              ELSE
                tblind = odepth(jl)/(bpade+odepth(jl))
                itgas = INT(tblint*tblind+0.5_wp)
                odepth(jl) = tau_tbl(itgas)
                atrans(jl,lev) = 1._wp - exp_tbl(itgas)
                tfacgas = tfn_tbl(itgas)
                gassrc(jl) = atrans(jl,lev) * plfrac * (blay + tfacgas*dplankdn)

                odtot(jl) = odepth(jl) + secdiff(jl,ib) * taucloud(jl,lev,ib)
                tblind = odtot(jl)/(bpade+odtot(jl))
                ittot = INT(tblint*tblind + 0.5_wp)
                tfactot = tfn_tbl(ittot)
                bbdtot(jl) = plfrac * (blay + tfactot*dplankdn)
                bbd(jl) = plfrac*(blay+tfacgas*dplankdn)
                atot(jl,lev) = 1._wp - exp_tbl(ittot)

                bbugas(jl,lev) = plfrac * (blay + tfacgas * dplankup)
                bbutot(jl,lev) = plfrac * (blay + tfactot * dplankup)
              ENDIF

              IF (istcldd(jl,lev) .EQ. 1) THEN
                cldradd(jl) = cldfrac(jl,lev) * radld(jl)
                clrradd(jl) = radld(jl) - cldradd(jl)
                oldcld(jl) = cldradd(jl)
                oldclr(jl) = clrradd(jl)
                rad(jl) = 0._wp
              ENDIF
              ttot = 1._wp - atot(jl,lev)
              cldsrc(jl) = bbdtot(jl) * atot(jl,lev)
              cldradd(jl) = cldradd(jl) * ttot + cldfrac(jl,lev) * cldsrc(jl)
              clrradd(jl) = clrradd(jl) * (1._wp-atrans(jl,lev)) + &
                & (1._wp-cldfrac(jl,lev))*gassrc(jl)
              radld(jl) = cldradd(jl) + clrradd(jl)
              drad(jl,lev-1) = drad(jl,lev-1) + radld(jl)

              radmod(jl) = rad(jl) * &
                & (facclr1d(jl,lev-1) * (1._wp-atrans(jl,lev)) + &
                & faccld1d(jl,lev-1) *  ttot) - &
                & faccmb1d(jl,lev-1) * gassrc(jl) + &
                & faccmb2d(jl,lev-1) * cldsrc(jl)

              oldcld(jl) = cldradd(jl) - radmod(jl)
              oldclr(jl) = clrradd(jl) + radmod(jl)
              rad(jl) = -radmod(jl) + facclr2d(jl,lev-1)*oldclr(jl) &
                & - faccld2d(jl,lev-1)*oldcld(jl)
              cldradd(jl) = cldradd(jl) + rad(jl)
              clrradd(jl) = clrradd(jl) - rad(jl)
              !          enddo
            ELSE
              ! !CDIR NODEP
              ! Clear layer
              !          do iclear = 1, n_clearpoints(lev)
              !            jl = iclear_ind(iclear,lev)

              ib = ibv(jl)

              plfrac = fracs(jl,lev,igc)
              blay = planklay(jl,lev,iband)
              dplankup = planklev(jl,lev,iband) - blay
              dplankdn = planklev(jl,lev-1,iband) - blay
              odepth(jl) = secdiff(jl,iband) * taut(jl,lev,igc)

              IF (odepth(jl) .LT. 0.0_wp) odepth(jl) = 0.0_wp

              IF (odepth(jl) .LE. 0.06_wp) THEN
                atrans(jl,lev) = odepth(jl)-0.5_wp*odepth(jl)*odepth(jl)
                odepth(jl) = rec_6*odepth(jl)
                bbd(jl) = plfrac*(blay+dplankdn*odepth(jl))
                bbugas(jl,lev) = plfrac*(blay+dplankup*odepth(jl))
              ELSE
                tblind = odepth(jl)/(bpade+odepth(jl))
                itr = INT(tblint*tblind+0.5_wp)
                transc = exp_tbl(itr)
                atrans(jl,lev) = 1._wp-transc
                tausfac = tfn_tbl(itr)
                bbd(jl) = plfrac*(blay+tausfac*dplankdn)
                bbugas(jl,lev) = plfrac * (blay + tausfac * dplankup)
              ENDIF
              radld(jl) = radld(jl) + (bbd(jl)-radld(jl))*atrans(jl,lev)
              drad(jl,lev-1) = drad(jl,lev-1) + radld(jl)
            ENDIF
          ENDDO

        ENDIF
        !          endif
        !  Set clear sky stream to total sky stream as long as layers
        !  remain clear.  Streams diverge when a cloud is reached (iclddn=1),
        !  and clear sky stream must be computed separately from that point.
        DO jl = 1, kproma
          IF (iclddn(jl).EQ.1) THEN
            radclrd(jl) = radclrd(jl) + (bbd(jl)-radclrd(jl)) * atrans(jl,lev)
            clrdrad(jl,lev-1) = clrdrad(jl,lev-1) + radclrd(jl)
          ELSE
            radclrd(jl) = radld(jl)
            clrdrad(jl,lev-1) = drad(jl,lev-1)
          ENDIF
        ENDDO
      ENDDO

      DO jl = 1, kproma
        ! Spectral emissivity & reflectance
        !  Include the contribution of spectrally varying longwave emissivity
        !  and reflection from the surface to the upward radiative transfer.
        !  Note: Spectral and Lambertian reflection are identical for the
        !  diffusivity angle flux integration used here.
        !  Note: The emissivity is applied to plankbnd and dplankbnd_dt when
        !  they are defined in subroutine setcoef.

        rad0 = fracs(jl,1,igc) * plankbnd(jl,iband)

        !  Add in reflection of surface downward radiance.
        reflect = 1._wp - semiss(jl,iband)
        radlu(jl) = rad0 + reflect * radld(jl)
        radclru(jl) = rad0 + reflect * radclrd(jl)

        ! Upward radiative transfer loop.

        urad(jl,0) = urad(jl,0) + radlu(jl)
        clrurad(jl,0) = clrurad(jl,0) + radclru(jl)
      ENDDO
      IF (idrv .EQ. 1) THEN
        DO jl = 1, kproma
          d_rad0_dt = fracs(jl,1,igc) * dplankbnd_dt(jl,iband)
          d_radlu_dt(jl) = d_rad0_dt
          d_urad_dt(jl,0) = d_urad_dt(jl,0) + d_radlu_dt(jl)
          d_radclru_dt(jl) = d_rad0_dt
          d_clrurad_dt(jl,0) = d_clrurad_dt(jl,0) + d_radclru_dt(jl)
        ENDDO
      ENDIF

      DO lev = 1, nlayers
        IF (n_cloudpoints(lev) == kproma) THEN ! all points are cloudy

          DO jl = 1, kproma ! Thus, direct addressing can be used

            gassrc(jl) = bbugas(jl,lev) * atrans(jl,lev)
            IF (istcld(jl,lev) .EQ. 1) THEN
              cldradu(jl) = cldfrac(jl,lev) * radlu(jl)
              clrradu(jl) = radlu(jl) - cldradu(jl)
              oldcld(jl) = cldradu(jl)
              oldclr(jl) = clrradu(jl)
              rad(jl) = 0._wp
            ENDIF
            ttot = 1._wp - atot(jl,lev)
            cldsrc(jl) = bbutot(jl,lev) * atot(jl,lev)
            cldradu(jl) = cldradu(jl) * ttot + cldfrac(jl,lev) * cldsrc(jl)
            clrradu(jl) = clrradu(jl) * (1.0_wp-atrans(jl,lev))+(1._wp-cldfrac(jl,lev))*gassrc(jl)
            ! Total sky radiance
            radlu(jl) = cldradu(jl) + clrradu(jl)
            urad(jl,lev) = urad(jl,lev) + radlu(jl)
            radmod(jl) = rad(jl) * &
              & (facclr1(jl,lev+1)*(1.0_wp-atrans(jl,lev))+ &
              & faccld1(jl,lev+1) *  ttot) - &
              & faccmb1(jl,lev+1) * gassrc(jl) + &
              & faccmb2(jl,lev+1) * cldsrc(jl)
            oldcld(jl) = cldradu(jl) - radmod(jl)
            oldclr(jl) = clrradu(jl) + radmod(jl)
            rad(jl) = -radmod(jl) + facclr2(jl,lev+1)*oldclr(jl) - faccld2(jl,lev+1)*oldcld(jl)
            cldradu(jl) = cldradu(jl) + rad(jl)
            clrradu(jl) = clrradu(jl) - rad(jl)
            IF (idrv .EQ. 1) THEN
              d_radlu_dt(jl) = d_radlu_dt(jl) * cldfrac(jl,lev) * (1.0_wp - atot(jl,lev)) + &
                & d_radlu_dt(jl) * (1.0_wp - cldfrac(jl,lev)) * (1.0_wp - atrans(jl,lev))
              d_urad_dt(jl,lev) = d_urad_dt(jl,lev) + d_radlu_dt(jl)
            ENDIF
          ENDDO

        ELSE IF (n_clearpoints(lev) == kproma) THEN ! all points are clear

          DO jl = 1, kproma ! thus, direct addressing can be used
            radlu(jl) = radlu(jl) + (bbugas(jl,lev)-radlu(jl))*atrans(jl,lev)
            urad(jl,lev) = urad(jl,lev) + radlu(jl)
            IF (idrv .EQ. 1) THEN
              d_radlu_dt(jl) = d_radlu_dt(jl) * (1.0_wp - atrans(jl,lev))
              d_urad_dt(jl,lev) = d_urad_dt(jl,lev) + d_radlu_dt(jl)
            ENDIF
          ENDDO

        ELSE ! both cloudy and clear points are in the vector; use index list
          ! or if clause
          ! !CDIR NODEP
          ! Cloudy layer
          !          do icld = 1, n_cloudpoints(lev)
          !            jl = icld_ind(icld,lev)
          DO jl = 1, kproma
            IF (icldlyr(jl,lev) .EQ. 1) THEN

              gassrc(jl) = bbugas(jl,lev) * atrans(jl,lev)
              IF (istcld(jl,lev) .EQ. 1) THEN
                cldradu(jl) = cldfrac(jl,lev) * radlu(jl)
                clrradu(jl) = radlu(jl) - cldradu(jl)
                oldcld(jl) = cldradu(jl)
                oldclr(jl) = clrradu(jl)
                rad(jl) = 0._wp
              ENDIF
              ttot = 1._wp - atot(jl,lev)
              cldsrc(jl) = bbutot(jl,lev) * atot(jl,lev)
              cldradu(jl) = cldradu(jl) * ttot + cldfrac(jl,lev) * cldsrc(jl)
              clrradu(jl) = clrradu(jl) * (1.0_wp-atrans(jl,lev)) &
                & +(1._wp-cldfrac(jl,lev))*gassrc(jl)
              ! Total sky radiance
              radlu(jl) = cldradu(jl) + clrradu(jl)
              urad(jl,lev) = urad(jl,lev) + radlu(jl)
              radmod(jl) = rad(jl) * &
                & (facclr1(jl,lev+1)*(1.0_wp-atrans(jl,lev))+ &
                & faccld1(jl,lev+1) *  ttot) - &
                & faccmb1(jl,lev+1) * gassrc(jl) + &
                & faccmb2(jl,lev+1) * cldsrc(jl)
              oldcld(jl) = cldradu(jl) - radmod(jl)
              oldclr(jl) = clrradu(jl) + radmod(jl)
              rad(jl) = -radmod(jl) + facclr2(jl,lev+1)*oldclr(jl) - &
                & faccld2(jl,lev+1)*oldcld(jl)
              cldradu(jl) = cldradu(jl) + rad(jl)
              clrradu(jl) = clrradu(jl) - rad(jl)
              IF (idrv .EQ. 1) THEN
                d_radlu_dt(jl) = d_radlu_dt(jl) * cldfrac(jl,lev) * (1.0_wp - atot(jl,lev)) + &
                  & d_radlu_dt(jl) * (1.0_wp - cldfrac(jl,lev)) * (1.0_wp - atrans(jl,lev))
                d_urad_dt(jl,lev) = d_urad_dt(jl,lev) + d_radlu_dt(jl)
              ENDIF
              !          enddo
            ELSE
              ! !CDIR NODEP
              ! Clear layer
              !          do iclear = 1, n_clearpoints(lev)
              !            jl = iclear_ind(iclear,lev)

              radlu(jl) = radlu(jl) + (bbugas(jl,lev)-radlu(jl))*atrans(jl,lev)
              urad(jl,lev) = urad(jl,lev) + radlu(jl)
              IF (idrv .EQ. 1) THEN
                d_radlu_dt(jl) = d_radlu_dt(jl) * (1.0_wp - atrans(jl,lev))
                d_urad_dt(jl,lev) = d_urad_dt(jl,lev) + d_radlu_dt(jl)
              ENDIF
            ENDIF
          ENDDO

        ENDIF
        !  Set clear sky stream to total sky stream as long as all layers
        !  are clear (iclddn=0).  Streams must be calculated separately at
        !  all layers when a cloud is present (iclddn=1), because surface
        !  reflectance is different for each stream.
        DO jl = 1, kproma
          IF (iclddn(jl).EQ.1) THEN
            radclru(jl) = radclru(jl) + (bbugas(jl,lev)-radclru(jl))*atrans(jl,lev)
            clrurad(jl,lev) = clrurad(jl,lev) + radclru(jl)
          ELSE
            radclru(jl) = radlu(jl)
            clrurad(jl,lev) = urad(jl,lev)
          ENDIF
        ENDDO
        IF (idrv .EQ. 1) THEN
          DO jl = 1, kproma
            IF (iclddn(jl).EQ.1) THEN
              d_radclru_dt(jl) = d_radclru_dt(jl) * (1.0_wp - atrans(jl,lev))
              d_clrurad_dt(jl,lev) = d_clrurad_dt(jl,lev) + d_radclru_dt(jl)
            ELSE
              d_radclru_dt(jl) = d_radlu_dt(jl)
              d_clrurad_dt(jl,lev) = d_urad_dt(jl,lev)
            ENDIF
          ENDDO
        ENDIF
      ENDDO

      ! Increment g-point counter
      igc = igc + 1
      ! Return to continue radiative transfer for all g-channels in present band
      IF (igc .LE. ngs(iband)) GO TO 1000

      ! Process longwave output from band.
      ! Calculate upward, downward, and net flux.
      DO lev = nlayers, 0, -1
        DO jl = 1, kproma
          uflux = urad(jl,lev)*wtdiff
          dflux = drad(jl,lev)*wtdiff
          urad(jl,lev) = 0.0_wp
          drad(jl,lev) = 0.0_wp
          totuflux(jl,lev) = totuflux(jl,lev) + uflux * delwave(iband)
          totdflux(jl,lev) = totdflux(jl,lev) + dflux * delwave(iband)
          uclfl = clrurad(jl,lev)*wtdiff
          dclfl = clrdrad(jl,lev)*wtdiff
          clrurad(jl,lev) = 0.0_wp
          clrdrad(jl,lev) = 0.0_wp
          totuclfl(jl,lev) = totuclfl(jl,lev) + uclfl * delwave(iband)
          totdclfl(jl,lev) = totdclfl(jl,lev) + dclfl * delwave(iband)
        ENDDO
      ENDDO

      ! Calculate total change in upward flux wrt surface temperature
      IF (idrv .EQ. 1) THEN
        DO lev = nlayers, 0, -1
          DO jl = 1, kproma
            duflux_dt = d_urad_dt(jl,lev) * wtdiff
            d_urad_dt(jl,lev) = 0.0_wp
            dtotuflux_dt(jl,lev) = dtotuflux_dt(jl,lev) + duflux_dt * delwave(iband) * fluxfac

            duclfl_dt = d_clrurad_dt(jl,lev) * wtdiff
            d_clrurad_dt(jl,lev) = 0.0_wp
            dtotuclfl_dt(jl,lev) = dtotuclfl_dt(jl,lev) + duclfl_dt * delwave(iband) * fluxfac
          ENDDO
        ENDDO
      ENDIF

      ! End spectral band loop
    ENDDO

    ! Calculate fluxes at surface
    DO jl = 1, kproma  ! loop over columns
      totuflux(jl,0) = totuflux(jl,0) * fluxfac
      totdflux(jl,0) = totdflux(jl,0) * fluxfac
      fnet(jl,0) = totuflux(jl,0) - totdflux(jl,0)

      totuclfl(jl,0) = totuclfl(jl,0) * fluxfac
      totdclfl(jl,0) = totdclfl(jl,0) * fluxfac
      fnetc(jl,0) = totuclfl(jl,0) - totdclfl(jl,0)
    ENDDO

    ! Calculate fluxes at model levels
    DO lev = 1, nlayers
      DO jl = 1, kproma  ! loop over columns
        totuflux(jl,lev) = totuflux(jl,lev) * fluxfac
        totdflux(jl,lev) = totdflux(jl,lev) * fluxfac
        fnet(jl,lev) = totuflux(jl,lev) - totdflux(jl,lev)
        totuclfl(jl,lev) = totuclfl(jl,lev) * fluxfac
        totdclfl(jl,lev) = totdclfl(jl,lev) * fluxfac
        fnetc(jl,lev) = totuclfl(jl,lev) - totdclfl(jl,lev)
      ENDDO
    ENDDO

  END SUBROUTINE lrtm_rtrnmr_2
  !-----------------------------------------------------------------------------
#endif

END MODULE mo_lrtm_rtrnmr

