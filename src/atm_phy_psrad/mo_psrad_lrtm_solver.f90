!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
!! @brief Module to compute longwave radiative transfer (no scattering)
!!
!! @remarks
!!   This module contains routines that compute longwave radiative transfer. 
!!   The algorithms follow those used by RRTMG but have been simplified and vectorized
!!   accross columns. 
!!
!! @author Robert Pincus, U. Colorado, visiting MPI-M, Hamburg (2011-07)
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
MODULE mo_psrad_lrtm_solver
  USE mo_psrad_general,             ONLY : wp, pi, nbndlw
  USE mo_psrad_fastmath,   ONLY : transmit, tautrans 

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: lrtm_solver, find_secdiff

CONTAINS
  SUBROUTINE lrtm_solver(jcs, kproma, kbdim, klev, tau, layer_f, level_f, &
    weights, secant, surface_f, surface_eps, flux_up, flux_down)

    ! Compute IR (no scattering) radiative transfer for a set of columns 
    ! Based on AER code RRTMG_LW_RTNMC, including approximations used there
    ! Layers are ordered from botton to top (i.e. tau(1) is tau in lowest layer)
    ! Computes all-sky RT given a total optical thickness in each layer
    !
    INTEGER,  INTENT(in) :: jcs, kproma, kbdim, klev
      
    REAL(wp), DIMENSION(KBDIM, klev), INTENT(in) :: &
      tau, & ! Longwave optical thickness
      layer_f, & ! Planck function at layer centers
      weights ! Fraction of total Planck function for this g-point 
    REAL(wp), INTENT(in) :: &
      ! Planck function at layer edges, level i is the top of layer i
      level_f(KBDIM,klev+1)
    REAL(wp), DIMENSION(KBDIM), INTENT(in) :: &
      surface_f, & ! Planck function at surface
      surface_eps, & ! Surface emissivity
      secant ! secant of integration angle - 
      !...depends on band, column water vapor
    REAL(wp), DIMENSION(KBDIM,klev+1), INTENT(out) :: &
      flux_up, & ! Fluxes at the interfaces 
      flux_down
                                         
    INTEGER  :: jk
      
    REAL(wp) :: &
      layer_transmissivity(KBDIM,klev), & ! Layer transmissivity
      tfn(KBDIM,klev), & ! TFN_TBL 
                ! Tau transition function; i.e. the transition of the Planck
                ! function from that for the mean layer temperature to that for
                ! the layer boundary temperature as a function of optical depth.
                ! The "linear in tau" method is used to make the table.
      tmp(KBDIM) 

    ! This secant and weight corresponds to the standard diffusivity 
    ! angle.  The angle is redefined for some bands.
    REAL(wp), PARAMETER :: wtdiff = 0.5_wp!, rec_6 = 0.166667_wp
    REAL(wp) :: fudge
    fudge = 2.0e+04_wp * pi

    ! Plank function derivatives and total emission for linear-in-tau 
    !approximation
    DO jk = 1, klev
      ! Weight optical depth by 1/cos(diffusivity angle), which depends on band 
      ! This will be used to compute layer transmittance
      !BUG? Bjorn: you set this to 0!
      !layer_effective_ir_optical_depth(:,jk) = max(1.e-9_wp,...
      !tmp == layer_effective_ir_optical_depth      
      tmp(jcs:kproma) = max(0.0_wp, secant(jcs:kproma) * tau(jcs:kproma,jk))
      tfn(jcs:kproma,jk) = tautrans(tmp(jcs:kproma), kproma-jcs+1)
      layer_transmissivity(jcs:kproma,jk) = transmit(tmp(jcs:kproma), kproma-jcs+1)
    END DO         

    ! Downward radiative transfer
    ! Radiance down at propagation angle
    flux_down(:, klev+1) = 0. ! Upper boundary condition - no downwelling IR
    DO jk = klev, 1, -1
      ! Interpolated downward emission 
      tmp(jcs:kproma) = weights(jcs:kproma,jk) * &
        (layer_f(jcs:kproma,jk) + (level_f(jcs:kproma,jk) - layer_f(jcs:kproma,jk)) * tfn(jcs:kproma,jk))
      flux_down(jcs:kproma,jk) = flux_down(jcs:kproma,jk+1) + layer_transmissivity(jcs:kproma,jk) * &
        (tmp(jcs:kproma) - flux_down(jcs:kproma,jk+1)) 
    END DO 

    ! Surface contribution, including reflection
    flux_up(jcs:kproma, 1) = weights(jcs:kproma, 1) * surface_eps(jcs:kproma) * surface_f(jcs:kproma) + &
      (1._wp - surface_eps(jcs:kproma)) * flux_down(jcs:kproma, 1)
    
    ! Upward radiative transfer
    ! Radiance up at propagation angle
    DO jk = 1, klev
      ! Interpolated upward emission 
      tmp(jcs:kproma) = weights(jcs:kproma,jk) * &
        (layer_f(jcs:kproma,jk) + (level_f(jcs:kproma,jk+1) - layer_f(jcs:kproma,jk)) * tfn(jcs:kproma,jk))
      flux_up(jcs:kproma,jk+1) = flux_up(jcs:kproma,jk) * &
        (1._wp - layer_transmissivity(jcs:kproma,jk)) + &
        layer_transmissivity(jcs:kproma,jk) * tmp(jcs:kproma)
    END DO 
    
    ! Covert intensities at diffusivity angles (radiance) to fluxes
    flux_up(jcs:kproma,:) = flux_up(jcs:kproma,:) * wtdiff * fudge
    flux_down(jcs:kproma,:) = flux_down(jcs:kproma,:) * wtdiff * fudge
    
  END SUBROUTINE lrtm_solver

  ELEMENTAL FUNCTION find_secdiff(iband, pwvcm) 
    INTEGER,  INTENT(in) :: iband ! RRTMG LW band number
    REAL(wp),  INTENT(in) :: pwvcm ! Precipitable water vapor (cm) 
    REAL(wp) :: find_secdiff

    ! Compute diffusivity angle for Bands 2-3 and 5-9 to vary (between 1.50
    ! and 1.80) as a function of total column water vapor.  The function
    ! has been defined to minimize flux and cooling rate errors in these bands
    ! over a wide range of precipitable water values.
    REAL(wp), DIMENSION(nbndlw), PARAMETER :: &
      a0 = (/ 1.66_wp,  1.55_wp,  1.58_wp,  1.66_wp, &
        1.54_wp, 1.454_wp,  1.89_wp,  1.33_wp, &
        1.668_wp, 1.66_wp,  1.66_wp,  1.66_wp, &
        1.66_wp,  1.66_wp,  1.66_wp,  1.66_wp /), &
      a1 = (/ 0.00_wp,  0.25_wp,  0.22_wp,  0.00_wp, &
        0.13_wp, 0.446_wp, -0.10_wp,  0.40_wp, &
        -0.006_wp, 0.00_wp,  0.00_wp,  0.00_wp, &
        0.00_wp,  0.00_wp,  0.00_wp,  0.00_wp /), &
      a2 = (/ 0.00_wp, -12.0_wp, -11.7_wp,  0.00_wp, &
        -0.72_wp,-0.243_wp,  0.19_wp,-0.062_wp, &
        0.414_wp, 0.00_wp,  0.00_wp,  0.00_wp, &
        0.00_wp,  0.00_wp,  0.00_wp, 0.00_wp /)

    if (iband == 1 .or. iband == 4 .or. iband >= 10) then
      find_secdiff = 1.66_wp
    else
      find_secdiff = MAX(1.5_wp, &
        MIN(a0(iband) + a1(iband) * exp(a2(iband)*pwvcm), 1.8_wp))
    endif

  END FUNCTION find_secdiff
  ! -------------------------------------------------------------------------------
END MODULE mo_psrad_lrtm_solver
