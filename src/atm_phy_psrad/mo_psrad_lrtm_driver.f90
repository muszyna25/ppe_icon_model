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
  USE mo_psrad_radiation_parameters, &
                             ONLY : i_overlap, l_do_sep_clear_sky, rad_undef
  USE mo_psrad_lrtm_setup,   ONLY: ngb, delwave
  USE mo_psrad_lrtm_kgs, ONLY: totplanck
  USE mo_psrad_rrtm_coeffs,  ONLY: lrtm_coeffs
  USE mo_psrad_lrtm_gas_optics, ONLY: gas_optics_lw
  USE mo_psrad_lrtm_solver,  ONLY: lrtm_solver, find_secdiff
  USE mo_psrad_cld_sampling, ONLY: sample_cld_state
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
    wkl, wx, coldry, emis, cldfr, taucld, tauaer, rnseeds, uflx, dflx, &
    uflxc, dflxc)

    use mo_psrad_general, only: ngptlw, ngas, ih2o, ico2, nreact
    INTEGER, INTENT(in) :: &
      kbdim, & !< Maximum block length
      kproma, & !< Number of horizontal columns
      klev !< Number of model layers

    REAL(wp), INTENT(in) :: &
      play(KBDIM,klev), & !< Layer pressures [hPa, mb] (KBDIM,klev)
      psfc(KBDIM), & !< Surface pressure [hPa, mb] (KBDIM)
      tlay(KBDIM,klev), & !< Layer temperatures [K] (KBDIM,klev)
      tlev(KBDIM,klev+1), & !< Interface temperatures [K] (KBDIM,klev+1)
      tsfc(KBDIM), & !< Surface temperature [K] (KBDIM)
      wkl(KBDIM,klev,ngas), & !< Gas volume mixing ratios
      wx(KBDIM,ncfc,klev), & !< CFC type gas volume mixing ratios
      coldry(KBDIM,klev), & !< Column dry amount
      emis(KBDIM,nbndlw), & !< Surface emissivity  (KBDIM,nbndlw)
      cldfr(KBDIM,klev), & !< Cloud fraction  (KBDIM,klev)
      taucld(KBDIM,klev,nbndlw), & !< Coud optical depth (KBDIM,klev,nbndlw)
      tauaer(KBDIM,klev,nbndlw)    !< Aerosol optical depth (KBDIM,klev,nbndlw)

    !INTEGER, INTENT(in) :: icld(KBDIM,klev)
    ! Seeds for random number generator (KBDIM,:) 
    INTEGER, INTENT(INOUT) :: rnseeds(:,:) 

    ! Longwave fluxes in [W/m2]
    REAL(wp), DIMENSION(KBDIM,klev+1), INTENT(out) :: & 
      uflx, & ! Total sky upward
      dflx, & ! Total sky downward
      uflxc, & ! Clear sky upward
      dflxc ! Clearr sky downward

    REAL(wp) :: taug(KBDIM,klev,ngptlw), & ! gas optical depth 
      pwvcm(KBDIM), & ! precipitable water vapor [cm]
      secdiff(KBDIM), & ! diffusivity angle for RT calculation 
      planklay(KBDIM,klev,nbndlw), & ! Planck function at mid-layer
      planklev(KBDIM,klev+1,nbndlw), & ! Planck function at level interfaces
      plankbnd(KBDIM,nbndlw) ! Planck function at surface
    REAL(wp), DIMENSION(KBDIM,klev,ngptlw) :: &
      fracs, & ! Planck fraction per g-point
      taut, & ! gaseous + aerosol optical depths for all columns
      tautot !< cloud + gaseous + aerosol optical depths for all columns
    REAL(wp), DIMENSION(KBDIM,klev+1) :: & 
      zgpcd, & ! < gpoint clearsky downward flux
      zgpcu, & ! < gpoint clearsky downward flux
      zgpfd, & ! < gpoint fullsky downward flux
      zgpfu ! < gpoint fullsky downward flux

    ! -----------------
    ! Variables for gas optics calculations
    INTEGER :: laytrop (KBDIM), & !< tropopause layer index
      iabs(KBDIM,2,2,klev)
    INTEGER, DIMENSION(KBDIM,klev) :: jp, jp1, jt, jt1, &
      indminor

    REAL(wp), DIMENSION(KBDIM,klev) :: colbrd, &
      minorfrac, scaleminor, scaleminorn2, wbrodl
    REAL(wp), DIMENSION(KBDIM,klev,2) :: h2o_factor,h2o_fraction
    INTEGER, DIMENSION(KBDIM,klev,2) :: h2o_index
         
    ! Normalized CFC amounts (molecules/cm^2) 
    REAL(wp) :: wx_loc(KBDIM,SIZE(wx, 2), SIZE(wx, 3)) 

    REAL(wp) :: fac(KBDIM,2,2,klev), gases(KBDIM,klev,ngas), &
      ratio(KBDIM,2,klev,nreact) 

    INTEGER :: jl, jk, ig, ib

    !REAL(wp), PARAMETER :: cldmin = 1.e-20_wp ! minimum val for clouds

    ! Variables for sampling strategy 
    REAL(WP) :: clrSky_scaling(1:KBDIM) !, gpt_scaling 
    REAL(WP) :: smp_tau(KBDIM, klev, ngptlw) 
    LOGICAL :: &
      cldMask(KBDIM,klev,ngptlw), & ! cloud mask in each cell
      colcldMask(KBDIM,ngptlw) ! cloud mask for each column

    ! Cloud optical depth is only saved for the band associated 
    ! with this g-point
    ! We sample clouds first because we may want to adjust water vapor based 
    ! on presence/absence of clouds
    ! BUG: should pass icld flag instead of cldfr
    CALL sample_cld_state(kproma, KBDIM, klev, ngptlw, rnseeds(:,:), &
      i_overlap, cldfr(:,:), cldMask(:,:,:))
!IBM* ASSERT(NODEPS)
    DO ig = 1, ngptlw
      ib = ngb(ig)
      DO jl = 1, kproma
        smp_tau(jl,:,ig) = &
          MERGE(taucld(jl,:,ib), 0._wp, cldMask(jl,:,ig)) 
      END DO
    END DO
    IF (kproma /= kbdim) THEN
      smp_tau(kproma+1:kbdim,:,:) = 0
    ENDIF
    
    ! Cloud masks for sorting out clear skies - by cell and by column
    IF(.not. l_do_sep_clear_sky) THEN
      ! Are any layers cloudy? 
      colcldMask(:,1:ngptlw) = &
        ANY(cldMask(:,:,1:ngptlw), DIM=2)
      ! Clear-sky scaling is gpt_scaling/frac_clr or 0 if all samples 
      ! are cloudy 
      !clrSky_scaling(:) = gpt_scaling * MERGE( &
      clrSky_scaling(:) = MERGE( &
          REAL(ngptlw,KIND=wp) / &
            REAL(ngptlw - COUNT(colCldMask(:,:),DIM=2),KIND=wp), &
          0._wp, ANY(.not. colCldMask(:,:),DIM=2))
    END IF

    ! Calculate information needed by the radiative transfer routine
    ! that is specific to this atmosphere, especially some of the 
    ! coefficients and indices needed to compute the optical depths
    ! by interpolating data from stored reference atmospheres. 
    ! The coefficients are functions of temperature and pressure and 
    ! remain the same for all g-point samples.
    ! If gas concentrations, temperatures, or pressures vary with sample (ig) 
    ! the coefficients need to be calculated inside the loop over samples

    ! Broadening gases -- the number of molecules per cm^2 of all gases 
    ! not specified explicitly (water is excluded) 
    wbrodl(:,:) = coldry(:,:) - &
      SUM(wkl(:,:,(/ico2,io3,in2o,ich4,io2/)), DIM=3)
    CALL lrtm_coeffs(KBDIM, klev, play, tlay, coldry, wkl, wbrodl, &
      laytrop, jp, jp1, jt, jt1, iabs, gases, colbrd, fac, ratio, &
      h2o_factor, h2o_fraction, h2o_index, minorfrac, &
      scaleminor, scaleminorn2, indminor)
         
      !  Loop over g-points calculating gas optical properties. 
!IBM* ASSERT(NODEPS)
    wx_loc(:,:,:) = 1.e-20_wp * wx(:,:,:)
    CALL gas_optics_lw(KBDIM, klev, play, wx_loc, &
      coldry, laytrop, jp1, iabs, gases, colbrd, fac, ratio, &
      h2o_factor, h2o_fraction, h2o_index, &
      minorfrac, scaleminor, scaleminorn2, indminor,fracs,taug)

    DO ig = 1, ngptlw
      ib   = ngb(ig) 
      DO jl = 1, kbdim
        ! Gas concentrations in colxx variables are normalized by 1.e-20_wp 
        ! in lrtm_coeffs; CFC gas concentrations (wx) need the same 
        ! normalization. Per Eli Mlawer the k values used in gas optics 
        ! tables have been multiplied by 1e20
        DO jk = 1, klev
          taut (jl,jk,ig) = taug(jl,jk,ig) + tauaer(jl,jk,ib)
        ENDDO
      END DO
    END DO
    ! All-sky optical depth. Mask for 0 cloud optical depth? 
    tautot(:,:,:) = taut(:,:,:) + smp_tau(:,:,:) 
    
    ! Compute radiative transfer.
    !
    uflx (:,:) = 0.0_wp
    dflx (:,:) = 0.0_wp
    uflxc(:,:) = 0.0_wp
    dflxc(:,:) = 0.0_wp
    
    ! Planck function in each band at layers and boundaries
!IBM* ASSERT(NODEPS)
    DO ig = 1, nbndlw
      planklay(:,:,ig) = planckFunction(tlay(:,:),ig) 
      planklev(:,:,ig) = planckFunction(tlev(:,:),ig) 
      plankbnd(:,ig) = planckFunction(tsfc(:),ig) 
    END DO
        
    ! Precipitable water vapor in each column - this can affect the 
    ! integration angle secdiff
    pwvcm(:) = ((amw * SUM(wkl(:,:,ih2o), DIM=2)) / & 
      (amd * SUM(coldry(:,:) + wkl(:,:,ih2o), DIM=2))) * & 
      (1.e3_wp * psfc(:)) / (1.e2_wp * grav)

    ! Compute radiative transfer for each set of samples
    DO ig = 1, ngptlw
      ib = ngb(ig)
      secdiff(:) = find_secdiff(ib, pwvcm(:))
      
      ! All sky fluxes
      CALL lrtm_solver(KBDIM, klev, tautot(:,:,ig), &
        planklay(:,:,ib), planklev(:,:,ib), fracs(:,:,ig), secdiff, &
        plankbnd(:,ib), emis(:,ib), zgpfu, zgpfd)
       
      uflx(:,:) = uflx (:,:) + zgpfu(:,:)
      dflx(:,:) = dflx (:,:) + zgpfd(:,:) 
	                           
      ! Clear-sky fluxes
      IF(l_do_sep_clear_sky) THEN
        ! Remove clouds and do second RT calculation
        CALL lrtm_solver(KBDIM, klev, taut(:,:,ig), &
          planklay(:,:,ib), planklev(:,:,ib), fracs(:,:,ig), secdiff, &
          plankbnd(:,ib), emis(:,ib), zgpcu, zgpcd)
        uflxc(:,:) = uflxc(:,:) + zgpcu(:,:) 
        dflxc(:,:) = dflxc(:,:) + zgpcd(:,:) 
      ELSE
        ! Accumulate fluxes by excluding cloudy subcolumns, 
        ! weighting to account for smaller sample size
!IBM* ASSERT(NODEPS)
        DO jk = 0, klev
          uflxc(:,jk) = uflxc(:,jk) + MERGE(0._wp, &
            zgpfu(:,jk) * clrSky_scaling(:), & 
            colCldMask(:,ig))
          dflxc(:,jk) = dflxc(:,jk) + MERGE(0._wp, & 
            zgpfd(:,jk) * clrSky_scaling(:), & 
            colCldMask(:,ig)) 
        END DO 
      END IF 
    END DO

    ! If computing clear-sky fluxes from samples, flag any columns 
    ! where all samples were cloudy
    IF(.not. l_do_sep_clear_sky) THEN 
!IBM* ASSERT(NODEPS)
       DO jl = 1, kproma
          IF(ALL(colCldMask(jl,:))) THEN
             uflxc(jl,:) = rad_undef
             dflxc(jl,:) = rad_undef
          END IF
       END DO
    END IF
  END SUBROUTINE lrtm

  !----------------------------------------------------------------------------
  ELEMENTAL FUNCTION planckFunction(temp, band)
    ! Compute the blackbody emission in a given band as a 
    ! function of temperature
    REAL(WP), INTENT(IN) :: temp
    INTEGER,  INTENT(IN) :: band 
    REAL(WP) :: planckFunction
    
    INTEGER  :: index
    REAL(WP) :: fraction 
    
    index = MIN(MAX(1, INT(temp - 159._wp)),180)
    fraction = temp - 159._wp - float(index)
    
    planckFunction = totplanck(index, band) + &
      fraction * (totplanck(index+1, band) - totplanck(index, band))
    planckFunction = planckFunction * delwave(band)
  END FUNCTION planckFunction

END MODULE mo_psrad_lrtm_driver

