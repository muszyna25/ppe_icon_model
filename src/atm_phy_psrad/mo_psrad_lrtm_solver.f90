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
#include "psrad_fastmath.inc"
MODULE mo_psrad_lrtm_solver
  USE mo_psrad_general, ONLY : wp, nbndlw

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: lrtm_solver, fill_secdiff
  REAL(wp), PARAMETER :: rec_6  = 1._wp/6._wp, &
    fixed_secdiff = 1.66_wp
  
CONTAINS

! Compute IR (no scattering) radiative transfer for a set of columns 
! Based on AER code RRTMG_LW_RTNMC, including approximations used there
! Layers are ordered from botton to top (i.e. tau(1) is tau in lowest layer)
! Computes all-sky RT given a total optical thickness in each layer
  SUBROUTINE lrtm_solver(kproma, kbdim, klev, tau, &
    P_lay, P_lev, fracs, use_fixed_secdiff, secdiff, P_sfc, emis, & 
    flux_up, flux_dn)
    USE mo_psrad_general, ONLY : jTOA, jSFC, jINC, jABOVE, jBELOW
    INTEGER, INTENT(IN) :: kproma, kbdim, klev
    INTEGER, INTENT(IN) :: use_fixed_secdiff
      
    REAL(wp), DIMENSION(KBDIM,klev), INTENT(in) :: &
      tau, & ! Longwave optical thickness
      P_lay, & ! Planck function at layer centers
      fracs ! Fraction of total Planck function for this g-point 
      
    REAL(wp), INTENT(IN) :: &
      P_lev(KBDIM,klev+1) ! Planck function at layer edges, level i is the top of layer i
      
    REAL(wp), DIMENSION(KBDIM), INTENT(IN) :: &
      P_sfc, & ! Planck function at surface
      emis, & ! Surface emissivity
      secdiff ! secant of integration angle - depends on band, column water vapor
    REAL(wp), INTENT(INOUT) :: &
      flux_up(KBDIM,klev+1), & ! Fluxes at the interfaces 
      flux_dn(KBDIM,klev+1) 
                                         
    ! -----------
    INTEGER  :: i, j
      
    REAL(wp), DIMENSION(KBDIM,klev) :: &
      trans ! Layer transmissivity

    REAL(wp) :: &
      tmp, tfn, &
      effective_tau   ! Effective IR optical depth of layer

    ! Level 1 is closest to the ground
    !rad_dn -> flux_dn(KBDIM,klev+1), & ! Radiance down at propagation angle
    flux_dn(:, jTOA+jABOVE) = 0. ! Upper boundary condition - no downwelling IR

    !DO j = klev, 1, -1
    DO j = jTOA, jSFC, -jINC
    DO i = 1, kproma
      IF (use_fixed_secdiff == 1) THEN
        effective_tau = max(1.e-9_wp, fixed_secdiff * tau(i,j))
      ELSE
        effective_tau = max(1.e-9_wp, secdiff(i) * tau(i,j))
      ENDIF
      INV_EXPON(effective_tau, tmp)
      trans(i,j) = 1.0_wp - tmp
      tfn = effective_tau * rec_6
      IF (effective_tau > 1.0e-3_wp) THEN
        tfn = 1.0_wp - 2.0_wp * (1._wp/effective_tau - tmp/trans(i,j))
      ENDIF

      flux_up(i,j+jABOVE) = trans(i,j) * fracs(i,j) * (P_lay(i,j) + &
        (P_lev(i,j+jABOVE) - P_lay(i,j)) * tfn)
      flux_dn(i,j+jBELOW) = flux_dn(i,j+jABOVE) + (fracs(i,j) * &
        (P_lay(i,j) + (P_lev(i,j+jBELOW) - P_lay(i,j)) * tfn) - &
        flux_dn(i,j+jABOVE)) * trans(i,j)
    END DO 
    END DO 

    DO i = 1, kproma
      flux_up(i, jSFC+jBELOW) = fracs(i,jSFC) * emis(i) * P_sfc(i) + &
        (1._wp - emis(i)) * flux_dn(i,jSFC+jBELOW)
    END DO 

    !DO j = 1, klev
    DO j = jSFC, jTOA, jINC
    DO i = 1, kproma
      flux_up(i,j+jABOVE) = flux_up(i,j+jABOVE) + &
        flux_up(i,j+jBELOW) * (1._wp - trans(i,j))
    END DO 
    END DO 
    
  END SUBROUTINE lrtm_solver
  ! -------------------------------------------------------------------------------
  SUBROUTINE fill_secdiff(kproma, band, pwvcm, secdiff) 
    INTEGER, INTENT(IN) :: kproma, band
    REAL(wp), INTENT(IN) :: pwvcm(:) ! Precipitable water vapor (cm) 
    REAL(wp), INTENT(INOUT) :: secdiff(:)

    ! Compute diffusivity angle for Bands 2-3 and 5-9 to vary (between 1.50
    ! and 1.80) as a function of total column water vapor.  The function
    ! has been defined to minimize flux and cooling rate errors in these bands
    ! over a wide range of precipitable water values.
    REAL(wp), DIMENSION(3,nbndlw), PARAMETER :: a = RESHAPE((/ &
      1.66_wp, 0.00_wp, 0.00_wp, &
      1.55_wp, 0.25_wp, -12.0_wp, &
      1.58_wp, 0.22_wp, -11.7_wp, &
      1.66_wp, 0.00_wp, 0.00_wp, &
      1.54_wp, 0.13_wp, -0.72_wp, &
      1.454_wp, 0.446_wp, -0.243_wp, &
      1.89_wp, -0.10_wp, 0.19_wp, &
      1.33_wp, 0.40_wp, -0.062_wp, &
      1.668_wp, -0.006_wp, 0.414_wp, &
      1.66_wp, 0.00_wp, 0.00_wp, &
      1.66_wp, 0.00_wp, 0.00_wp, &
      1.66_wp, 0.00_wp, 0.00_wp, &
      1.66_wp, 0.00_wp, 0.00_wp, &
      1.66_wp, 0.00_wp, 0.00_wp, &
      1.66_wp, 0.00_wp, 0.00_wp, &
      1.66_wp, 0.00_wp, 0.00_wp/), SHAPE=(/3,nbndlw/))

    INTEGER :: i

    DO i = 1, kproma
      secdiff(i) = MAX(1.5_wp, &
        MIN(a(1,band) + a(2,band) * exp(a(3,band)*pwvcm(i)), 1.8_wp))
    ENDDO

  END SUBROUTINE fill_secdiff

END MODULE mo_psrad_lrtm_solver
