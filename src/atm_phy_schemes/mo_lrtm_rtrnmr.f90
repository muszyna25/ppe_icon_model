!     path:      $Source: /storm/rc1/cvsroot/rc/rrtmg_lw/src/rrtmg_lw_rtrnmr.f90,v $
!     author:    $Author: mike $
!     revision:  $Revision: 1.7 $
!     created:   $Date: 2009/11/12 20:52:26 $
!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
#if defined __xlC__ && !defined NOXLFPROCESS
@PROCESS HOT
@PROCESS NOSTRICT
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

  USE mo_math_constants,   ONLY : pi

  USE mo_lrtm_par,         ONLY : nbndlw, delwave, ngs
  USE mo_lrtm_setup,       ONLY : ntbl, bpade, tau_tbl, exp_tbl, tfn_tbl

  IMPLICIT NONE

  PRIVATE


  PUBLIC :: lrtm_rtrnmr

  REAL(wp), PARAMETER :: fluxfac = 2.0e+04_wp * pi

  REAL(wp), PARAMETER :: tblint = REAL(ntbl, wp)

CONTAINS

  !-----------------------------------------------------------------------------
  SUBROUTINE lrtm_rtrnmr(&
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
    REAL(wp) :: dplankup(kproma,nlayers), dplankdn(kproma,nlayers)
    REAL(wp), PARAMETER :: rec_6 = 0.166667_wp, &
      ! diffusivity angle adjustment coefficients
      a0(nbndlw) = (/ 1.66_wp,  1.55_wp,  1.58_wp,  1.66_wp, &
      & 1.54_wp, 1.454_wp,  1.89_wp,  1.33_wp, &
      & 1.668_wp,  1.66_wp,  1.66_wp,  1.66_wp, &
      & 1.66_wp,  1.66_wp,  1.66_wp,  1.66_wp /),&
      a1(nbndlw) = (/ 0.00_wp,  0.25_wp,  0.22_wp,  0.00_wp, &
      & 0.13_wp, 0.446_wp, -0.10_wp,  0.40_wp, &
      & -0.006_wp,  0.00_wp,  0.00_wp,  0.00_wp, &
      & 0.00_wp,  0.00_wp,  0.00_wp,  0.00_wp /),&
      a2(nbndlw) = (/ 0.00_wp, -12.0_wp, -11.7_wp,  0.00_wp, &
      & -0.72_wp,-0.243_wp,  0.19_wp,-0.062_wp, &
      & 0.414_wp,  0.00_wp,  0.00_wp,  0.00_wp, &
      & 0.00_wp,  0.00_wp,  0.00_wp,  0.00_wp /), &
    ! This secant and weight corresponds to the standard diffusivity
    ! angle.  This initial value is redefined below for some bands.
      wtdiff = 0.5_wp
    REAL(wp) :: radld(kproma), radclrd(kproma), plfrac
    REAL(wp) :: odepth(kproma), odtot(kproma), odepth_rec, odtot_rec, &
         gassrc(kproma), ttot, odepth_temp
    REAL(wp) :: tfactot, bbd(kproma), bbdtot(kproma), tfacgas
    REAL(wp) :: rad0, reflect, radlu(kproma), radclru(kproma)

    REAL(wp) :: duflux_dt
    REAL(wp) :: duclfl_dt
    REAL(wp) :: d_urad_dt(kproma,0:nlayers)
    REAL(wp) :: d_clrurad_dt(kproma,0:nlayers)
    REAL(wp) :: d_rad0_dt, d_radlu_dt(kproma), d_radclru_dt(kproma)

    LOGICAL :: lcldlyr(kproma,nlayers)             ! flag for cloud in layer
    INTEGER :: ibnd, ib, iband, lay, lev    ! loop indices
    INTEGER :: igc                          ! g-point interval counter
    LOGICAL :: iclddn(kproma)               ! flag for cloud in down path
    INTEGER :: ittot, itgas            ! lookup table indices
    INTEGER, PARAMETER :: ipat(16,0:2) = RESHAPE((/&
      ! These arrays indicate the spectral 'region' (used in the
      ! calculation of ice cloud optical depths) corresponding
      ! to each spectral band.  See cldprop.f for more details.
      & 1,1,1,1,1,1,1,1,1, 1, 1, 1, 1, 1, 1, 1, &
      & 1,2,3,3,3,4,4,4,5, 5, 5, 5, 5, 5, 5, 5, &
      & 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16/), (/16, 3/))

    INTEGER :: jl

    INTEGER :: ibv(kproma)

    ! Declarations for cloud overlap adjustment
    REAL(wp) :: faccld1(kproma,nlayers+1),faccld2(kproma,nlayers+1)
    REAL(wp) :: facclr1(kproma,nlayers+1),facclr2(kproma,nlayers+1)
    REAL(wp) :: faccmb1(kproma,nlayers+1),faccmb2(kproma,nlayers+1)
    REAL(wp) :: faccld1d(kproma,0:nlayers),faccld2d(kproma,0:nlayers)
    REAL(wp) :: facclr1d(kproma,0:nlayers),facclr2d(kproma,0:nlayers)
    REAL(wp) :: faccmb1d(kproma,0:nlayers),faccmb2d(kproma,0:nlayers)
    !--------------------------------------------------------------------------
    ! Maximum/Random cloud overlap variables
    ! for upward radiative transfer
    !  facclr2  fraction of clear radiance from previous layer that needs to
    !           be switched to cloudy stream
    !  facclr1  fraction of the radiance that had been switched in the previous
    !           layer from cloudy to clear that needs to be switched back to
    !           cloudy in the current layer
    !  faccld2  fraction of cloudy radiance from previous layer that needs to
    !           be switched to clear stream
    !  faccld1  fraction of the radiance that had been switched in the previous
    !           layer from clear to cloudy that needs to be switched back to
    !           clear in the current layer
    ! for downward radiative transfer
    !  facclr2d fraction of clear radiance from previous layer that needs to
    !           be switched to cloudy stream
    !  facclr1d fraction of the radiance that had been switched in the previous
    !           layer from cloudy to clear that needs to be switched back to
    !           cloudy in the current layer
    !  faccld2d fraction of cloudy radiance from previous layer that needs to
    !           be switched to clear stream
    !  faccld1d fraction of the radiance that had been switched in the previous
    !           layer from clear to cloudy that needs to be switched back to
    !           clear in the current layer
    !--------------------------------------------------------------------------

    REAL(wp) :: odepth_rec_or_tfacgas, odtot_rec_or_tfactot, cldsrc, &
         oldclr, oldcld, clrradd_temp, cldradd_temp, radmod, &
         odepth_rec_or_tausfac
    REAL(wp), DIMENSION(kproma) :: clrradd, cldradd, clrradu, cldradu, &
      & rad

    LOGICAL :: istcld(kproma,nlayers+1),istcldd(kproma,0:nlayers), &
         branch_od1, branch_od2

    ! ------- Definitions -------
    ! input
    !    nlayers                      ! number of model layers
    !    ngptlw                       ! total number of g-point subintervals
    !    nbndlw                       ! number of longwave spectral bands
    !    ncbands                      ! number of spectral bands for clouds
    !    secdiff                      ! diffusivity angle
    !    wtdiff                       ! weight for radiance to flux conversion
    !    pavel                        ! layer pressures (mb)
    !    tavel                        ! layer temperatures (k)
    !    tz                           ! level (interface) temperatures(mb)
    !    tbound                       ! surface temperature (k)
    !    cldfrac                      ! layer cloud fraction
    !    taucloud                     ! layer cloud optical depth
    !    lcldlyr                      ! flag for cloudy layers
    !    iclddn                       ! flag for cloud in column at any layer
    !    semiss                       ! surface emissivities for each band
    !    reflect                      ! surface reflectance
    !    bpade                        ! 1/(pade constant)
    !    tau_tbl                      ! clear sky optical depth look-up table
    !    exp_tbl                      ! exponential look-up table for transmittance
    !    tfn_tbl                      ! tau transition function look-up table

    ! local
    !    atrans                       ! gaseous absorptivity
    !    atot                         ! combined gaseous and cloud absorptivity
    !    odclr                        ! clear sky (gaseous) optical depth
    !    odcld                        ! cloud optical depth
    !    odtot                        ! optical depth of gas and cloud
    !    tfacgas                      ! gas-only pade factor, used for planck fn
    !    tfactot                      ! gas and cloud pade factor, used for planck fn
    !    bbdgas                       ! gas-only planck function for downward rt
    !    bbugas                       ! gas-only planck function for upward rt
    !    bbdtot                       ! gas and cloud planck function for downward rt
    !    bbutot                       ! gas and cloud planck function for upward calc.
    !    gassrc                       ! source radiance due to gas only
    !    radlu                        ! spectrally summed upward radiance
    !    radclru                      ! spectrally summed clear sky upward radiance
    !    urad                         ! upward radiance by layer
    !    clrurad                      ! clear sky upward radiance by layer
    !    radld                        ! spectrally summed downward radiance
    !    radclrd                      ! spectrally summed clear sky downward radiance
    !    drad                         ! downward radiance by layer
    !    clrdrad                      ! clear sky downward radiance by layer
    !    d_radlu_dt                   ! spectrally summed upward radiance
    !    d_radclru_dt                 ! spectrally summed clear sky upward radiance
    !    d_urad_dt                    ! upward radiance by layer
    !    d_clrurad_dt                 ! clear sky upward radiance by layer

    ! output
    !    totuflux                     ! upward longwave flux (w/m2)
    !    totdflux                     ! downward longwave flux (w/m2)
    !    fnet                         ! net longwave flux (w/m2)
    !    totuclfl                     ! clear sky upward longwave flux (w/m2)
    !    totdclfl                     ! clear sky downward longwave flux (w/m2)
    !    fnetc                        ! clear sky net longwave flux (w/m2)
    !    dtotuflux_dt                 ! change in upward longwave flux (w/m2/k)
    !                                 ! with respect to surface temperature
    !    dtotuclfl_dt                 ! change in clear sky upward longwave flux (w/m2/k)
    !

    ! Local variables for cloud / no cloud index lists
    INTEGER, DIMENSION(kproma,nlayers) :: icld_ind,iclear_ind
    INTEGER :: icld, iclear, n_cloudpoints(nlayers), n_clearpoints(nlayers)


    ! Reset diffusivity angle for Bands 2-3 and 5-9 to vary (between 1.50
    ! and 1.80) as a function of total column water vapor.  The function
    ! has been defined to minimize flux and cooling rate errors in these bands
    ! over a wide range of precipitable water values.

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
          lcldlyr(jl,lay) = .TRUE.
          icld = icld + 1
          icld_ind(icld,lay) = jl
        ELSE
          lcldlyr(jl,lay) = .FALSE.
          iclear = iclear + 1
          iclear_ind(iclear,lay) = jl
        ENDIF
      ENDDO

      n_cloudpoints(lay) = icld
      n_clearpoints(lay) = iclear

    ENDDO

    ! Maximum/Random cloud overlap parameter

    istcld(:,1) = .TRUE.
    istcldd(:,nlayers) = .TRUE.

    CALL compute_radiative_transfer(kproma, nlayers, 1, cldfrac, &
      n_cloudpoints, n_clearpoints, icld_ind,iclear_ind, &
      1, nlayers, 1, &
      faccld1, faccld2, facclr1, facclr2, faccmb1, faccmb2, istcld)

    CALL compute_radiative_transfer(kproma, nlayers, 0, cldfrac, &
      n_cloudpoints, n_clearpoints, icld_ind,iclear_ind, &
      nlayers, 1, -1, &
      faccld1d, faccld2d, facclr1d, facclr2d, faccmb1d, faccmb2d, istcldd)


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
      iclddn(:) = .FALSE.
      clrradd = 0._wp
      cldradd = 0._wp

      ! Downward radiative transfer loop.
      DO lev = nlayers, 1, -1
        IF (n_clearpoints(lev) /= kproma) THEN
          DO jl = 1, kproma ! Thus, direct addressing can be used
            ib = ibv(jl)
            plfrac = fracs(jl,lev,igc)
            odepth_temp = MAX(0.0_wp, secdiff(jl,iband) * taut(jl,lev,igc))
            branch_od2 = odepth_temp .LE. 0.06_wp

            iclddn(jl) = iclddn(jl) .OR. lcldlyr(jl,lev)
            IF (lcldlyr(jl,lev)) THEN
              odtot(jl) = odepth_temp + secdiff(jl,ib) * taucloud(jl,lev,ib)
              branch_od1 = odtot(jl) .LT. 0.06_wp
              itgas = INT(tblint * odepth_temp/(bpade+odepth_temp) + 0.5_wp)
              tfacgas = tfn_tbl(itgas)
              ittot = MERGE(0, INT(tblint * odtot(jl)/(bpade+odtot(jl)) + 0.5_wp), branch_od1)
              tfactot = MERGE(0.0_wp, tfn_tbl(ittot), branch_od1)
              odepth(jl) = MERGE(odepth_temp, tau_tbl(itgas), branch_od1 .OR. branch_od2)

              odepth_rec = rec_6*odepth(jl)
              odtot_rec = MERGE(rec_6*odtot(jl), 0.0_wp, branch_od1)
              odepth_rec_or_tfacgas = MERGE(odepth_rec, tfacgas, branch_od1 .OR. branch_od2)
              odtot_rec_or_tfactot = MERGE(odtot_rec, tfactot, branch_od1)

              atot(jl,lev) = MERGE(odtot(jl) - 0.5_wp*odtot(jl)*odtot(jl), &
                &               1._wp - exp_tbl(ittot), branch_od1)

              atrans(jl,lev) = MERGE(odepth(jl) - 0.5_wp*odepth(jl)*odepth(jl), &
                &           1._wp - exp_tbl(itgas), branch_od1 .OR. branch_od2)
              bbdtot(jl) = plfrac * (planklay(jl,lev,iband) + odtot_rec_or_tfactot * dplankdn(jl,lev))
              bbd(jl) = plfrac * (planklay(jl,lev,iband) + odepth_rec_or_tfacgas * dplankdn(jl,lev))
              gassrc(jl) = plfrac * (planklay(jl,lev,iband) &
                + odepth_rec_or_tfacgas * dplankdn(jl,lev)) * atrans(jl,lev)
              bbugas(jl,lev) = plfrac * (planklay(jl,lev,iband) &
                + odepth_rec_or_tfacgas * dplankup(jl,lev))
              bbutot(jl,lev) = plfrac * (planklay(jl,lev,iband) &
                + odtot_rec_or_tfactot * dplankup(jl,lev))

              ttot = 1._wp - atot(jl,lev)
              cldsrc = bbdtot(jl) * atot(jl,lev)
              clrradd_temp = MERGE(radld(jl) - cldfrac(jl,lev) * radld(jl), &
                & clrradd(jl), istcldd(jl,lev)) * (1._wp-atrans(jl,lev)) + &
                & (1._wp-cldfrac(jl,lev))*gassrc(jl)
              cldradd_temp = MERGE(cldfrac(jl,lev) * radld(jl), cldradd(jl), &
                istcldd(jl,lev)) * ttot + cldfrac(jl,lev) * cldsrc
              radld(jl) = cldradd_temp + clrradd_temp

              radmod = MERGE(0._wp, rad(jl), istcldd(jl,lev)) * &
                & (facclr1d(jl,lev-1) * (1._wp-atrans(jl,lev)) + &
                & faccld1d(jl,lev-1) *  ttot) - &
                & faccmb1d(jl,lev-1) * gassrc(jl) + &
                & faccmb2d(jl,lev-1) * cldsrc

              oldcld = cldradd_temp - radmod
              oldclr = clrradd_temp + radmod
              rad(jl) = -radmod + facclr2d(jl,lev-1)*oldclr -&
                &  faccld2d(jl,lev-1)*oldcld
              cldradd(jl) = cldradd_temp + rad(jl)
              clrradd(jl) = clrradd_temp - rad(jl)
            ELSE
              ! only needed for NAG
              atot(jl, lev) = 0.0_wp
              bbutot(jl, lev) = 0.0_wp

              odepth_rec_or_tausfac = MERGE(rec_6*odepth_temp, &
                tfn_tbl(INT(tblint*odepth_temp/(bpade+odepth_temp)+0.5_wp)), &
                branch_od2)
              atrans(jl,lev) = MERGE(odepth_temp-0.5_wp*odepth_temp*odepth_temp, &
                1._wp - exp_tbl(INT(tblint*odepth_temp/(bpade+odepth_temp) + 0.5_wp)), &
                branch_od2)
              odepth(jl) = odepth_temp
              bbd(jl) = plfrac * (planklay(jl,lev,iband) &
                + odepth_rec_or_tausfac * dplankdn(jl,lev))
              bbugas(jl,lev) = plfrac * (planklay(jl,lev,iband) &
                + odepth_rec_or_tausfac * dplankup(jl,lev))
              radld(jl) = radld(jl) + (bbd(jl)-radld(jl))*atrans(jl,lev)

            END IF
            drad(jl,lev-1) = drad(jl,lev-1) + radld(jl)
          ENDDO

        ELSE

          DO jl = 1, kproma ! Thus, direct addressing can be used

            plfrac = fracs(jl,lev,igc)
            odepth_temp = MAX(0.0_wp, secdiff(jl,iband) * taut(jl,lev,igc))

            branch_od2 = odepth_temp .LE. 0.06_wp

            ! only needed for NAG
            atot(jl, lev) = 0.0_wp
            bbutot(jl, lev) = 0.0_wp

            odepth_rec_or_tausfac = MERGE(rec_6*odepth_temp, &
              tfn_tbl(INT(tblint*odepth_temp/(bpade+odepth_temp)+0.5_wp)), &
              branch_od2)
            atrans(jl,lev) = MERGE(odepth_temp-0.5_wp*odepth_temp*odepth_temp, &
              1._wp - exp_tbl(INT(tblint*odepth_temp/(bpade+odepth_temp) + 0.5_wp)), &
              branch_od2)
            odepth(jl) = odepth_temp
            bbd(jl) = plfrac * (planklay(jl,lev,iband) &
              + odepth_rec_or_tausfac * dplankdn(jl,lev))
            bbugas(jl,lev) = plfrac * (planklay(jl,lev,iband) &
              + odepth_rec_or_tausfac * dplankup(jl,lev))
            radld(jl) = radld(jl) + (bbd(jl)-radld(jl))*atrans(jl,lev)
            drad(jl,lev-1) = drad(jl,lev-1) + radld(jl)
          ENDDO

        ENDIF

        !  Set clear sky stream to total sky stream as long as layers
        !  remain clear.  Streams diverge when a cloud is reached (iclddn=1),
        !  and clear sky stream must be computed separately from that point.
        DO jl = 1, kproma
          IF (iclddn(jl)) THEN
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

          DO jl = 1, kproma
            gassrc(jl) = bbugas(jl,lev) * atrans(jl,lev)
            IF (istcld(jl,lev)) THEN
              cldradu(jl) = cldfrac(jl,lev) * radlu(jl)
              clrradu(jl) = radlu(jl) - cldradu(jl)
              rad(jl) = 0._wp
            ENDIF
            ttot = 1._wp - atot(jl,lev)
            cldsrc = bbutot(jl,lev) * atot(jl,lev)
            cldradu(jl) = cldradu(jl) * ttot + cldfrac(jl,lev) * cldsrc
            clrradu(jl) = clrradu(jl) * (1.0_wp-atrans(jl,lev))+(1._wp-cldfrac(jl,lev))*gassrc(jl)
            ! Total sky radiance
            ! Clear layer if lcldlyr(jl,lev) == .true.
            radlu(jl) = MERGE(cldradu(jl) + clrradu(jl), &
                 radlu(jl) + (bbugas(jl,lev)-radlu(jl))*atrans(jl,lev), &
                 lcldlyr(jl,lev))

            urad(jl,lev) = urad(jl,lev) + radlu(jl)
            radmod = rad(jl) * &
              & (facclr1(jl,lev+1)*(1.0_wp-atrans(jl,lev))+ &
              & faccld1(jl,lev+1) *  ttot) - &
              & faccmb1(jl,lev+1) * gassrc(jl) + &
              & faccmb2(jl,lev+1) * cldsrc
            oldcld = cldradu(jl) - radmod
            oldclr = clrradu(jl) + radmod
            rad(jl) = -radmod + facclr2(jl,lev+1)*oldclr - faccld2(jl,lev+1)*oldcld
            cldradu(jl) = cldradu(jl) + rad(jl)
            clrradu(jl) = clrradu(jl) - rad(jl)
          ENDDO

        !  Set clear sky stream to total sky stream as long as all layers
        !  are clear (iclddn=true).  Streams must be calculated separately at
        !  all layers when a cloud is present (iclddn=false), because surface
        !  reflectance is different for each stream.
        DO jl = 1, kproma
          radclru(jl) = MERGE(radclru(jl) + (bbugas(jl,lev)-radclru(jl))*atrans(jl,lev), &
               radlu(jl), iclddn(jl))
          clrurad(jl,lev) = MERGE(clrurad(jl,lev) + radclru(jl), urad(jl,lev), &
               iclddn(jl))
        ENDDO

        IF (idrv .EQ. 1) THEN
          DO jl = 1, kproma
            d_radlu_dt(jl) = d_radlu_dt(jl) * MERGE(cldfrac(jl,lev) * (1.0_wp - atot(jl,lev)) + &
                 & (1.0_wp - cldfrac(jl,lev)) * (1.0_wp - atrans(jl,lev)), &
                 & (1.0_wp - atrans(jl,lev)), lcldlyr(jl,lev))
            d_urad_dt(jl,lev) = d_urad_dt(jl,lev) + d_radlu_dt(jl)
            d_radclru_dt(jl) = MERGE(d_radclru_dt(jl) * (1.0_wp - atrans(jl,lev)), &
                 d_radlu_dt(jl), iclddn(jl))
            d_clrurad_dt(jl,lev) = MERGE(d_clrurad_dt(jl,lev) + d_radclru_dt(jl), &
                 d_urad_dt(jl,lev), iclddn(jl))
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
  END SUBROUTINE lrtm_rtrnmr

  SUBROUTINE compute_radiative_transfer(kproma, nlayers, ofs, cldfrac, &
       n_cloudpoints, n_clearpoints, icld_ind,iclear_ind, &
       start_lev, end_lev, lev_incr, &
       faccld1, faccld2, facclr1, facclr2, faccmb1, faccmb2, istcld)
    INTEGER, INTENT(in) :: kproma          ! number of columns
    INTEGER, INTENT(in) :: nlayers         ! total number of layers
    INTEGER, intent(in) :: ofs             ! layer offset to use for
                                           ! references to facc*
    INTEGER, INTENT(in) :: n_cloudpoints(nlayers), n_clearpoints(nlayers)
    INTEGER, DIMENSION(kproma,nlayers) :: icld_ind,iclear_ind
    INTEGER, INTENT(in) :: start_lev, end_lev, lev_incr


    REAL(wp), INTENT(inout) :: faccld1(kproma,nlayers+1), &
         faccld2(kproma,nlayers+1), &
         facclr1(kproma,nlayers+1), &
         facclr2(kproma,nlayers+1), &
         faccmb1(kproma,nlayers+1), &
         faccmb2(kproma,nlayers+1)
    LOGICAL, INTENT(out) :: istcld(kproma,nlayers+1)
    REAL(wp), INTENT(in) :: cldfrac(:,:)       ! layer cloud fraction

    REAL(wp) :: fmax, fmin, rat1(kproma), rat2(kproma)
    INTEGER :: lev, jl, olev, icld, iclear

    DO lev = start_lev, end_lev, lev_incr
      olev = lev + ofs
      IF (n_cloudpoints(lev) == kproma) THEN ! all points are cloudy
        DO jl = 1, kproma ! Thus, direct addressing can be used
          ! Maximum/random cloud overlap
          istcld(jl,olev) = .FALSE.
          IF (lev .EQ. end_lev) THEN
            faccld1(jl,olev) = 0._wp
            faccld2(jl,olev) = 0._wp
            facclr1(jl,olev) = 0._wp
            facclr2(jl,olev) = 0._wp
            faccmb1(jl,olev) = 0._wp
            faccmb2(jl,olev) = 0._wp
          ELSEIF (cldfrac(jl,olev) .GE. cldfrac(jl,lev)) THEN
            faccld1(jl,olev) = 0._wp
            faccld2(jl,olev) = 0._wp
            IF (istcld(jl,lev+1-ofs)) THEN
              facclr1(jl,olev) = 0._wp
              facclr2(jl,olev) = 0._wp
              IF (cldfrac(jl,lev) .LT. 1._wp) facclr2(jl,olev) = &
                & (cldfrac(jl,olev)-cldfrac(jl,lev))/(1._wp-cldfrac(jl,lev))
              facclr2(jl,lev+1-ofs) = 0._wp
              faccld2(jl,lev+1-ofs) = 0._wp
            ELSE
              fmax = MAX(cldfrac(jl,lev),cldfrac(jl,lev-lev_incr))
              IF (cldfrac(jl,olev) .GT. fmax) THEN
                facclr1(jl,olev) = rat2(jl)
                facclr2(jl,olev) = (cldfrac(jl,olev)-fmax)/(1._wp-fmax)
              ELSEIF (cldfrac(jl,olev) .LT. fmax) THEN
                facclr1(jl,olev) = (cldfrac(jl,olev)-cldfrac(jl,lev))/ &
                  & (cldfrac(jl,lev-lev_incr)-cldfrac(jl,lev))
                facclr2(jl,olev) = 0._wp
              ELSE
                facclr1(jl,olev) = rat2(jl)
                facclr2(jl,olev) = 0._wp
              ENDIF
            ENDIF
            rat1(jl) = MERGE(1._wp, 0._wp, &
                 facclr1(jl,olev) > 0._wp .OR. facclr2(jl,olev) > 0._wp)
            rat2(jl) = 0._wp
          ELSE
            facclr1(jl,olev) = 0._wp
            facclr2(jl,olev) = 0._wp
            IF (istcld(jl,lev+1-ofs)) THEN
              faccld1(jl,olev) = 0._wp
              faccld2(jl,olev) = (cldfrac(jl,lev)-cldfrac(jl,olev))/cldfrac(jl,lev)

              facclr2(jl,lev+1-ofs) = 0._wp
              faccld2(jl,lev+1-ofs) = 0._wp
            ELSE
              fmin = MIN(cldfrac(jl,lev),cldfrac(jl,lev-lev_incr))
              IF (cldfrac(jl,olev) .LE. fmin) THEN
                faccld1(jl,olev) = rat1(jl)
                faccld2(jl,olev) = (fmin-cldfrac(jl,olev))/fmin
              ELSE
                faccld1(jl,olev) = (cldfrac(jl,lev)-cldfrac(jl,olev))/(cldfrac(jl,lev)-fmin)
                faccld2(jl,olev) = 0._wp
              ENDIF
            ENDIF
            rat2(jl) = MERGE(1._wp, 0._wp, &
                 faccld1(jl,olev) > 0._wp .OR. faccld2(jl,olev) > 0._wp)
            rat1(jl) = 0._wp
          ENDIF
          IF (lev == start_lev) THEN
            faccmb1(jl,olev) = 0._wp
            faccmb2(jl,olev) = faccld1(jl,olev) * facclr2(jl,lev)
          ELSE
            faccmb1(jl,olev) = facclr1(jl,olev) * faccld2(jl,lev) * cldfrac(jl,lev-lev_incr)
            faccmb2(jl,olev) = faccld1(jl,olev) * facclr2(jl,lev) * (1._wp - cldfrac(jl,lev-lev_incr))
          ENDIF
        ENDDO
      ELSE IF (n_clearpoints(lev) == kproma) THEN ! all points are clear
        istcld(1:kproma,olev) = .TRUE.
      ELSE ! use index list for the case that not all points are cloudy
!CDIR NODEP,VOVERTAKE,VOB
        DO icld = 1, n_cloudpoints(lev)
          jl = icld_ind(icld,lev)
          ! Maximum/random cloud overlap
          istcld(jl,olev) = .FALSE.
          IF (lev .EQ. end_lev) THEN
            faccld1(jl,olev) = 0._wp
            faccld2(jl,olev) = 0._wp
            facclr1(jl,olev) = 0._wp
            facclr2(jl,olev) = 0._wp
            faccmb1(jl,olev) = 0._wp
            faccmb2(jl,olev) = 0._wp
          ELSEIF (cldfrac(jl,olev) .GE. cldfrac(jl,lev)) THEN
            faccld1(jl,olev) = 0._wp
            faccld2(jl,olev) = 0._wp
            IF (istcld(jl,lev+1-ofs)) THEN
              facclr1(jl,olev) = 0._wp
              facclr2(jl,olev) = 0._wp
              IF (cldfrac(jl,lev) .LT. 1._wp) facclr2(jl,olev) = &
                & (cldfrac(jl,olev)-cldfrac(jl,lev))/(1._wp-cldfrac(jl,lev))
              facclr2(jl,lev+1-ofs) = 0._wp
              faccld2(jl,lev+1-ofs) = 0._wp
            ELSE
              fmax = MAX(cldfrac(jl,lev),cldfrac(jl,lev-lev_incr))
              IF (cldfrac(jl,olev) .GT. fmax) THEN
                facclr1(jl,olev) = rat2(jl)
                facclr2(jl,olev) = (cldfrac(jl,olev)-fmax)/(1._wp-fmax)
              ELSEIF (cldfrac(jl,olev) .LT. fmax) THEN
                facclr1(jl,olev) = (cldfrac(jl,olev)-cldfrac(jl,lev))/ &
                  & (cldfrac(jl,lev-lev_incr)-cldfrac(jl,lev))
                facclr2(jl,olev) = 0._wp
              ELSE
                facclr1(jl,olev) = rat2(jl)
                facclr2(jl,olev) = 0._wp
              ENDIF
            ENDIF
            rat1(jl) = MERGE(1._wp, 0._wp, &
              facclr1(jl,olev) > 0._wp .OR. facclr2(jl,olev) > 0._wp)
            rat2(jl) = 0._wp
          ELSE
            facclr1(jl,olev) = 0._wp
            facclr2(jl,olev) = 0._wp
            IF (istcld(jl,lev+1-ofs)) THEN
              faccld1(jl,olev) = 0._wp
              faccld2(jl,olev) = (cldfrac(jl,lev)-cldfrac(jl,olev))/cldfrac(jl,lev)

              facclr2(jl,lev+1-ofs) = 0._wp
              faccld2(jl,lev+1-ofs) = 0._wp
            ELSE
              fmin = MIN(cldfrac(jl,lev),cldfrac(jl,lev-lev_incr))
              IF (cldfrac(jl,olev) .LE. fmin) THEN
                faccld1(jl,olev) = rat1(jl)
                faccld2(jl,olev) = (fmin-cldfrac(jl,olev))/fmin
              ELSE
                faccld1(jl,olev) = (cldfrac(jl,lev)-cldfrac(jl,olev))/(cldfrac(jl,lev)-fmin)
                faccld2(jl,olev) = 0._wp
              ENDIF
            ENDIF
            rat2(jl) = MERGE(1._wp, 0._wp, &
              faccld1(jl,olev) > 0._wp .OR. faccld2(jl,olev) > 0._wp)
            rat1(jl) = 0._wp
          ENDIF
          IF (lev == start_lev) THEN
            faccmb1(jl,olev) = 0._wp
            faccmb2(jl,olev) = faccld1(jl,olev) * facclr2(jl,lev)
          ELSE
            faccmb1(jl,olev) = facclr1(jl,olev) * faccld2(jl,lev) * cldfrac(jl,lev-lev_incr)
            faccmb2(jl,olev) = faccld1(jl,olev) * facclr2(jl,lev) * (1._wp - cldfrac(jl,lev-lev_incr))
          ENDIF
        ENDDO
!CDIR NODEP,VOVERTAKE,VOB
        DO iclear = 1, n_clearpoints(lev)
          jl = iclear_ind(iclear,lev)
          istcld(jl,olev) = .TRUE.
        ENDDO
      ENDIF

    ENDDO

  END SUBROUTINE compute_radiative_transfer

END MODULE mo_lrtm_rtrnmr
!
! Local Variables:
! f90-continuation-indent: 2
! End:
!
