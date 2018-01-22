!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
!! @brief Module to provide interface to rrtmg longwave radiation
!!
!! @remarks
!!   This module contains routines that provide the interface between ECHAM
!!   and the AER RRTMG radiation code.  Mostly it organizes and calculates the 
!!   information necessary to call the radiative transfer solvers.
!!
!! @author Bjorn Stevens, MPI-M, Hamburg (2010-08)
!!         Robert Pincus, U. Colorado, visiting MPI-M, Hamburg (2011-07)
!!
!! $ID: n/a$
!!
!! @par Origin
!!   Major segments of this code combines and rewrites (for the ICON standard) 
!!   code previously provided by AER and copyrighed by them.  The authors of the
!!   original AER code are: Eli J. Mlawer, Jennifer S. Delamere, Michael J. 
!!   Iacono and Shepard A. Clough with acknowledgments to Steven J. Taubman, 
!!   Karen Cady-Pereira, Patrick D. Brown, Ronald E. Farren, Luke Chen, Robert 
!!   Bergstrom. The rewrites were designed to better interface with the structure
!!   of the ICON family of models and elements of the ICON programming standard.
!!
!! @par Copyright
!!   The AER copyright
!!   on the original code is as follows: Copyright 2002-2009, Atmospheric and
!!   Environmental Research, Inc. (AER). This software may be used, copied, or
!!   redistributed as long as it is not sold and this copyright notice is
!!   reproduced on each copy made.  This model is provided as is without any
!!   express or implied warranties. (http://www.rtweb.aer.com/)               
!! 
!
MODULE mo_psrad_lrtm_driver

  USE mo_psrad_general
  USE mo_psrad_radiation_parameters, ONLY: i_overlap
  USE mo_psrad_radiation_parameters, ONLY: l_do_sep_clear_sky, rad_undef
  USE mo_psrad_lrtm_kgs, ONLY: ngpt, use_fixed_secdiff, &
    precipitable_vapor_factor
  USE mo_psrad_lrtm_gas_optics, ONLY: gas_optics_lw
  USE mo_psrad_lrtm_solver,  ONLY: lrtm_solver, fill_secdiff
  USE mo_psrad_cld_sampling, ONLY: sample_cld_state
#ifdef PSRAD_TIMING
  USE mo_timer, ONLY: ltimer, timer_start, timer_stop, &
    timer_sample_cloud_lw, timer_gas_optics_lw
#endif
#ifdef PSRAD_DEBUG
  USE mo_psrad_matlab
#endif
#ifdef PSRAD_DEVEL
#ifdef _OPENMP
  USE omp_lib, ONLY: omp_get_thread_num
#endif
#endif
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: lrtm

CONTAINS
  !-----------------------------------------------------------------------------
  !>
  !! @brief Prepares information for radiation call
  !! 
  !! @remarks: This program is the driver subroutine for the longwave radiative
  !! transfer routine.  This routine is adapted from the AER LW RRTMG_LW model
  !! that itself has been adapted from RRTM_LW for improved efficiency.  Our 
  !! routine does the spectral integration externally (the solver is explicitly
  !! called for each g-point, so as to facilitate sampling of g-points
  !! This routine:
  !!    1) calls INATM to read in the atmospheric profile from GCM;
  !!       all layering in RRTMG is ordered from surface to toa. 
  !!    2) calls COEFFS to calculate various quantities needed for 
  !!       the radiative transfer algorithm.  This routine is called only once for
  !!       any given thermodynamic state, i.e., it does not change if clouds chanege
  !!    3) calls TAUMOL to calculate gaseous optical depths for each 
  !!       of the 16 spectral bands, this is updated band by band.
  !!    4) calls SOLVER (for both clear and cloudy profiles) to perform the
  !!       radiative transfer calculation with a maximum-random cloud
  !!       overlap method, or calls RTRN to use random cloud overlap.
  !!    5) passes the necessary fluxes and cooling rates back to GCM
  !!
  !
  SUBROUTINE lrtm(kproma, kbdim, klev, play, psfc, tlay, tlev, tsfc, &
    gases, wx, coldry, emis, cldfr, taucld, tauaer, rnseeds, &
    ratio, scaleminorn2, fac, laytrop, iabs, &
    jp, indminor, h2o_factor, h2o_fraction, h2o_index, &
    colbrd, minorfrac, scaleminor, &
    flux_up, flux_dn, flux_up_clr, flux_dn_clr)

    USE mo_psrad_general, ONLY: ngptlw, ngas, ih2o, jTOA, jSFC, jINC 
#ifdef PSRAD_DEVEL
    USE mo_psrad_dump, ONLY: store
#endif
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: &
      kbdim, & ! Maximum block length
      kproma, & ! Number of horizontal columns
      klev ! Number of model layers

    REAL(wp), INTENT(IN) :: &
      play(KBDIM,klev), & ! Layer pressures [hPa, mb] (KBDIM,klev)
      psfc(KBDIM), & ! Surface pressure [hPa, mb] (KBDIM)
      tlay(KBDIM,klev), & ! Layer temperatures [K] (KBDIM,klev)
      tlev(KBDIM,klev+1), & ! Interface temperatures [K] (KBDIM,klev+1)
      tsfc(KBDIM), & ! Surface temperature [K] (KBDIM)
      wx(KBDIM,klev,ncfc), & ! CFC type gas volume mixing ratios
      coldry(KBDIM,klev), & ! Column dry amount
      emis(KBDIM,nbndlw), & ! Surface emissivity  (KBDIM,nbndlw)
      cldfr(KBDIM,klev), & ! Cloud fraction  (KBDIM,klev)
      taucld(KBDIM,klev,nbndlw), & ! Coud optical depth (KBDIM,klev,nbndlw)
      tauaer(KBDIM,klev,nbndlw)    ! Aerosol optical depth (KBDIM,klev,nbndlw)

    REAL(wp), TARGET, INTENT(IN) :: &
      gases(:,:,:), ratio(:,:,:,:), scaleminorn2(:,:)
    REAL(wp), INTENT(IN) :: fac(KBDIM,2,2,klev)
    INTEGER, INTENT(IN) :: laytrop (KBDIM), & ! tropopause layer index
      iabs(KBDIM,2,2,klev)
    INTEGER, DIMENSION(KBDIM,klev), INTENT(IN) :: jp, &
      indminor
    REAL(wp), DIMENSION(KBDIM,klev,2), INTENT(IN) :: h2o_factor,h2o_fraction
    INTEGER, DIMENSION(KBDIM,klev,2), INTENT(IN) :: h2o_index
    REAL(wp), DIMENSION(KBDIM,klev), INTENT(IN) :: colbrd, &
      minorfrac, scaleminor

    !INTEGER, INTENT(in) :: icld(KBDIM,klev)
    ! Seeds for random number generator (KBDIM,:) 
    INTEGER, INTENT(INOUT) :: rnseeds(:,:) 

    ! Longwave fluxes in [W/m2]
    REAL(wp), DIMENSION(KBDIM,klev+1), INTENT(INOUT) :: & 
      flux_up, & ! Total sky upward
      flux_dn, & ! Total sky downward
      flux_up_clr, & ! Clear sky upward
      flux_dn_clr ! Clearr sky downward

    REAL(wp) :: taug(KBDIM,klev,ngptlw), & ! gas optical depth 
      vapor(KBDIM), & ! precipitable water vapor [cm]
      secdiff(KBDIM), & ! diffusivity angle for RT calculation 
      P_lay(KBDIM,klev), & ! Planck function at mid-layer
      P_lev(KBDIM,klev+1), & ! Planck function at level interfaces
      P_sfc(KBDIM) ! Planck function at surface
    REAL(wp), DIMENSION(KBDIM,klev,ngptlw) :: &
      fracs, & ! Planck fraction per g-point
      taut ! gaseous + aerosol optical depths for all columns
    REAL(wp), DIMENSION(KBDIM,klev+1) :: & 
      partial_clear_dn, & ! < gpoint clearsky downward flux
      partial_clear_up, & ! < gpoint clearsky downward flux
      partial_full_dn, & ! < gpoint fullsky downward flux
      partial_full_up ! < gpoint fullsky downward flux

    ! -----------------
    ! Variables for gas optics calculations

         
    ! Normalized CFC amounts (molecules/cm^2) 
    !REAL(wp) :: wx_loc(KBDIM,SIZE(wx, 2), SIZE(wx, 3)) 


    INTEGER :: i, j, ig, gpt_in, band
    !REAL(wp), PARAMETER :: cldmin = 1.e-20_wp ! minimum val for clouds

    ! Variables for sampling strategy 
    REAL(WP) :: clear_sky_scaling(KBDIM) !, gpt_scaling 
    LOGICAL :: &
      col_cld_mask(KBDIM,ngptlw) ! cloud mask for each column
    LOGICAL :: &
      cld_mask(KBDIM,klev,ngptlw) ! cloud mask in each cell


    ! Cloud optical depth is only saved for the band associated 
    ! with this g-point
    ! We sample clouds first because we may want to adjust water vapor based 
    ! on presence/absence of clouds
    ! BUG: should pass icld flag instead of cldfr
#ifdef PSRAD_DEVEL
    INTEGER :: thread_id
#ifdef _OPENMP
    thread_id = omp_get_thread_num() + 1
#else
    thread_id = 1
#endif
#endif

#ifdef PSRAD_TIMING
    IF (ltimer) CALL timer_start(timer_sample_cloud_lw)
#endif
    CALL sample_cld_state(kproma, KBDIM, klev, &
      ngptlw, rnseeds, i_overlap, cldfr, cld_mask)
#ifdef PSRAD_TIMING
    IF (ltimer) CALL timer_stop(timer_sample_cloud_lw)
#endif
#ifdef PSRAD_DEVEL
    store(thread_id)%cld_mask_lw = 0
    WHERE(cld_mask(1:kproma,:,:))
      store(thread_id)%cld_mask_lw(1:kproma,:,:) = 1
    END WHERE
#endif

    ! Cloud masks for sorting out clear skies - by cell and by column

    ! FIXME: code for the l_do_sep_clear_sky = .false. case is poorly written
    ! and doesn't vectorize well - but is not currently used...
    IF(.not. l_do_sep_clear_sky) THEN
      ! Are any layers cloudy? 
      col_cld_mask(1:kproma,1:ngptlw) = &
        ANY(cld_mask(1:kproma,:,1:ngptlw), DIM=2)
      ! Clear-sky scaling is gpt_scaling/frac_clr or 0 if all samples 
      ! are cloudy 
      !clear_sky_scaling = gpt_scaling * MERGE( &
      clear_sky_scaling(1:kproma) = MERGE( &
          REAL(ngptlw, KIND=wp) / &
            REAL(ngptlw - COUNT(col_cld_mask, DIM=2), KIND=wp), &
          0._wp, ANY(.not. col_cld_mask(1:kproma,:), DIM=2))
    END IF
         
      !  Loop over g-points calculating gas optical properties. 
#ifdef PSRAD_TIMING
    IF (ltimer) CALL timer_start(timer_gas_optics_lw)
#endif
    CALL gas_optics_lw(kproma, KBDIM, klev, jp, fac, iabs, &
      laytrop, gases, h2o_factor, h2o_fraction, h2o_index, &
      play, wx, coldry, colbrd, ratio, &
      minorfrac, scaleminor, scaleminorn2, indminor, fracs, taug)
#ifdef PSRAD_TIMING
    IF (ltimer) CALL timer_stop(timer_gas_optics_lw)
#endif
#ifdef PSRAD_DEVEL
    store(thread_id)%laytrop_lw(1:kproma) = laytrop(1:kproma)
    store(thread_id)%laytrop_lw(kproma+1:kbdim) = 0
    store(thread_id)%tau_lw(1:kproma,:,:) = taug(1:kproma,:,:)
    store(thread_id)%tau_lw(kproma+1:kbdim,:,:) = 0.0
    store(thread_id)%fracs_lw(1:kproma,:,:) = fracs(1:kproma,:,:)
    store(thread_id)%fracs_lw(kproma+1:kbdim,:,:) = 0.0
#endif

    ig = 0
    DO band = 1,nbndlw
    DO gpt_in = 1, ngpt(band)
      ig = ig + 1
      DO i = 1, kproma  
        ! Gas concentrations in colxx variables are normalized by 1.e-20_wp 
        ! in lrtm_coeffs; CFC gas concentrations (wx) need the same 
        ! normalization. Per Eli Mlawer the k values used in gas optics 
        ! tables have been multiplied by 1e20
        DO j = 1, klev
          taut(i,j,ig) = taug(i,j,ig) + tauaer(i,j,band)
        ENDDO
      END DO
    END DO
    END DO

    ! Compute radiative transfer.
    !
    flux_up = 0.0_wp
    flux_dn = 0.0_wp
    flux_up_clr = 0.0_wp
    flux_dn_clr = 0.0_wp
    
        
    ! FIXME: this calculation seems wrong: relating the sum over 
    ! the whole column instead of a per-cell relation then summing?
    ! Precipitable water vapor in each column - this can affect the 
    ! integration angle secdiff
    vapor(1:kproma) = SUM(gases(1:kproma,jSFC:jTOA:jINC,ih2o), DIM=2)

    vapor(1:kproma) = precipitable_vapor_factor * (vapor(1:kproma) / &
      (SUM(coldry(1:kproma,jSFC:jTOA:jINC), DIM=2) + &
      vapor(1:kproma))) * psfc(1:kproma)

    ! Compute radiative transfer for each set of samples
    ig = 0
    DO band = 1,nbndlw
      IF (use_fixed_secdiff(band) == 0) THEN
        CALL fill_secdiff(kproma, KBDIM, band, vapor, secdiff)
      ENDIF
      ! Planck function in each band at layers and boundaries
      CALL fast_planck(kproma, KBDIM, klev, band, tlay, P_lay)
      CALL fast_planck(kproma, KBDIM, klev+1, band, tlev, P_lev)
      CALL fast_planck(kproma, KBDIM, 1, band, tsfc, P_sfc)
      DO gpt_in = 1, ngpt(band)
        ig = ig + 1

        IF(l_do_sep_clear_sky) THEN
          ! Remove clouds and do second RT calculation
          CALL lrtm_solver(kproma, KBDIM, klev, taut(:,:,ig), &
            P_lay, P_lev, fracs(:,:,ig), &
          use_fixed_secdiff(band), secdiff, &
            P_sfc, emis(:,band), partial_clear_up, partial_clear_dn)
          DO j = 1, klev+1
          DO i = 1, kproma
            flux_up_clr(i,j) = flux_up_clr(i,j) + partial_clear_up(i,j) 
            flux_dn_clr(i,j) = flux_dn_clr(i,j) + partial_clear_dn(i,j) 
          ENDDO
          ENDDO
        ENDIF
        
        DO j = 1, klev
        DO i = 1, kproma
          IF (cld_mask(i,j,ig)) THEN
            taut(i,j,ig) = taut(i,j,ig) + taucld(i,j,band)
          ENDIF
        ENDDO
        ENDDO
        ! All sky fluxes
        CALL lrtm_solver(kproma, KBDIM, klev, taut(:,:,ig), &
          P_lay, P_lev, fracs(:,:,ig), &
          use_fixed_secdiff(band), secdiff, &
          P_sfc, emis(:,band), partial_full_up, partial_full_dn)

        DO j = 1, klev+1
        DO i = 1, kproma
          flux_up(i,j) = flux_up(i,j) + partial_full_up(i,j)
          flux_dn(i,j) = flux_dn(i,j) + partial_full_dn(i,j)
        ENDDO
        ENDDO
                               
        ! Clear-sky fluxes
        IF (.not. l_do_sep_clear_sky) THEN
          ! Accumulate fluxes by excluding cloudy subcolumns, 
          ! weighting to account for smaller sample size
          DO j = 1, klev+1
          DO i = 1, kproma
            IF (.not. col_cld_mask(i,ig)) THEN
              flux_up_clr(i,j) = flux_up_clr(i,j) + &
                partial_full_up(i,j) * clear_sky_scaling(i)
              flux_dn_clr(i,j) = flux_dn_clr(i,j) + &
                partial_full_dn(i,j) * clear_sky_scaling(i)
            ENDIF
          END DO 
          END DO 
        END IF 
      END DO
    END DO

    ! If computing clear-sky fluxes from samples, flag any columns 
    ! where all samples were cloudy
    IF(.not. l_do_sep_clear_sky) THEN 
       DO i = 1, kproma
          IF(ALL(col_cld_mask(i,:))) THEN
             flux_up_clr(i,:) = rad_undef
             flux_dn_clr(i,:) = rad_undef
          END IF
       END DO
    END IF
  END SUBROUTINE lrtm

  SUBROUTINE fast_planck(kproma, kbdim, klev, band, T, F)
    USE mo_psrad_lrtm_kgs, ONLY: P => totplanck
    ! Compute the blackbody emission in a given band as a 
    ! function of temperature
    INTEGER,  INTENT(IN) :: kproma, kbdim, klev, band
    REAL(wp), INTENT(IN) :: T(KBDIM,klev)
    REAL(wp), INTENT(INOUT) :: F(KBDIM,klev)
    
    INTEGER  :: i, j, idx
    REAL(WP) :: frc
    
    DO j = 1, klev
    DO i = 1, kproma
      frc = T(i,j) - 159._wp
      idx = MIN(MAX(1, INT(frc)),180)
      frc = frc - idx
      F(i,j) = P(idx,band) + frc * (P(idx+1,band) - P(idx,band))
    ENDDO
    ENDDO
  END SUBROUTINE fast_planck

END MODULE mo_psrad_lrtm_driver

