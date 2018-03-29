!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
!! @brief Module to do shortwave radiative transfer calculation
!!
!! @remarks
!!   This module contains routines to compute shortwave radiative transfer 
!!     given optical properties. Two-stream methods provide
!!     layer transmittance and reflectance; adding computing flux profiles
!! 
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
!!   Many of the comments are in the original. 
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
MODULE mo_psrad_srtm_solver

  USE mo_psrad_general, ONLY: wp

  IMPLICIT NONE
  PRIVATE 
  
#include "psrad_fastmath.inc"

  REAL(WP), PARAMETER :: &
    ! Min single scattering albedo for conservative approx.
    w0Min = 0.9999995_wp

  PUBLIC :: srtm_solver_tr, srtm_reftra_ec

CONTAINS

! Solver using optical properties of each layer as input
! Following are comments from the AER code from which this is derived
! Purpose: Contains spectral loop to compute the shortwave radiative fluxes, 
!          using the two-stream method of H. Barker. 
! Method:
!    Adapted from two-stream model of H. Barker;
!    Two-stream model options (selected with kmodts in rrtmg_sw_reftra.F90):
!        1: Eddington, 2: PIFM, Zdunkowski et al., 3: discret ordinates
! Modifications:
! Original: H. Barker
! Revision: Merge with RRTMG_SW: J.-J.Morcrette, ECMWF, Feb 2003
! Revision: Add adjustment for Earth/Sun distance : MJIacono, AER, Oct 2003
! Revision: Bug fix for use of PALBP and PALBD: MJIacono, AER, Nov 2003
! Revision: Bug fix to apply delta scaling to clear sky: AER, Dec 2004
!           if routine cldprop is used; delta scaling can be applied by 
!           switching code below if cldprop is not used to get cloud 
!           properties. 
!           AER, Jan 2005
! Revision: Uniform formatting for RRTMG: MJIacono, AER, Jul 2006 
! Revision: Use exponential lookup table for transmittance: MJIacono, AER, 
!           Aug 2007 

! Solver using optical properties and pre-computed direct/diffuse
!   reflectanc/transmittance to compute fluxes layer-by-layer 
  SUBROUTINE srtm_solver_tr(kproma, kbdim, klev, &
    albdif, albdir, Rc, Rd, Tc, Td, Tb, &
    flux_dn, flux_up, direct_flux)

    USE mo_psrad_general, ONLY: jTOA, jSFC, jINC, jABOVE, jBELOW

    IMPLICIT NONE

    INTEGER,  INTENT(IN) :: kproma, kbdim, klev
    REAL(wp), DIMENSION(KBDIM), INTENT(IN) :: &
      albdif, & ! surface albedo (diffuse)
      albdir ! surface albedo (direct)
    ! Reflectance/Transmittance, collimated/diffuse 
    REAL(wp), DIMENSION(KBDIM,klev+1), INTENT(INOUT) :: &
      Rc, Rd, & 
      Tc, Td, &
      Tb ! Direct beam transmission
    REAL(wp), INTENT(INOUT) :: &
      flux_dn(KBDIM,klev+1), flux_up(KBDIM,klev+1), & ! Flux down, up
      direct_flux(KBDIM)

    INTEGER  :: i, jk
    REAL(wp), DIMENSION(KBDIM,klev+1) :: &
      accTb ! Direct beam transmittance total to each layer

    ! Surface values
    DO i = 1, kproma
      Rc(i,jSFC+jBELOW) = albdir(i)
      Rd(i,jSFC+jBELOW) = albdif(i)
      Tc(i,jSFC+jBELOW) = 0.0_wp
      Td(i,jSFC+jBELOW) = 0.0_wp
    ENDDO

    Tb(1:kproma,jSFC+jBELOW) = 0.0_wp
    ! Accumulated direct beam transmission 
    accTb(1:kproma,jTOA+jABOVE) = 1._wp
    DO jk=jTOA+jABOVE,jSFC+jABOVE, -jINC
    DO i = 1, kproma
      accTb(i,jk-jINC) = Tb(i,jk)*accTb(i,jk)
    ENDDO
    ENDDO

    DO i = 1, kproma
      direct_flux(i) = accTb(i,jSFC+jBELOW)
    ENDDO

    ! Vertical quadrature for cloudy fluxes
    CALL adding(kproma, KBDIM, klev, &
      Rc, Rd, Tc, Td, &      
      Tb, accTb, &
      albdir, albdif, flux_dn, flux_up) 
  END SUBROUTINE srtm_solver_tr

! "vertical quadrature" i.e. flux profiles based on layer-by-layer direct
! and diffuse reflectance and transmittance
! Equations are developed in doi:10.1002/qj.49712555316. 
! Modifications.
! Original: H. Barker
! Revision: Integrated with rrtmg_sw, J.-J. Morcrette, ECMWF, Oct 2002
! Revision: Reformatted for consistency with rrtmg_lw: MJIacono, AER, Jul 2006

  SUBROUTINE adding(kproma, kbdim, klev, Rc, Rd, Tc, Td, &
       Tb, accTb, albdir, albdif, flux_dn, flux_up)
    USE mo_psrad_general, ONLY: jTOA, jSFC, jINC, jABOVE, jBELOW
    INTEGER, INTENT (in) :: kproma, kbdim, klev
    ! klev+1 indicates surface
    REAL(wp), DIMENSION(KBDIM,klev+1), INTENT(IN) :: &
      Rc, Rd, &  ! direct amd diffuse beam reflectivity
      Tc, Td, &  ! direct and diffuse beam transmissivity
      Tb, accTb  ! direct beam transmission, by layer and accumulated
    ! Surface albedo for direct and diffuse incidence
    REAL(wp), INTENT(IN) :: albdir(KBDIM), albdif(KBDIM) 
    REAL(wp), DIMENSION(KBDIM,klev+1), INTENT(INOUT) :: &
      flux_dn, & ! downwelling flux (W/m2)
      flux_up ! upwelling   flux (W/m2)

    INTEGER  :: i, jk
    REAL(wp) :: zreflect(KBDIM)
    REAL(wp), DIMENSION(KBDIM,klev+1) :: ztdn, zrup, zrupd, zrdnd

    ! Link lowest layer with surface
    DO i = 1, kproma
      zrup (i,jSFC+jBELOW) = albdir(i)
      zrupd(i,jSFC+jBELOW) = albdif(i)
      zreflect(i) = 1._wp / (1._wp - Rd(i,jSFC+jBELOW) * Rd(i,jSFC+jABOVE))
      zrup (i,jSFC+jABOVE) = Rc(i,jSFC+jABOVE) +&
        (Td(i,jSFC+jABOVE) * ((Tc(i,jSFC+jABOVE) - Tb(i,jSFC+jABOVE)) * Rd(i,jSFC+jBELOW) + & 
         Tb (i,jSFC+jABOVE) * Rc(i,jSFC+jBELOW))) * zreflect(i)
      zrupd(i,jSFC+jABOVE) = Rd(i,jSFC+jABOVE) + &
        Td(i,jSFC+jABOVE) * Td(i,jSFC+jABOVE) * Rd(i,jSFC+jBELOW) * zreflect(i)
    ENDDO

    ! Pass from bottom to top 
    DO jk = jSFC+jINC, jTOA, jINC
      DO i = 1, kproma
        zreflect(i) = 1._wp / (1._wp - zrupd(i,jk+jBELOW) * Rd(i,jk+jABOVE))
        zrup(i,jk+jABOVE) = Rc(i,jk+jABOVE) + (Td(i,jk+jABOVE) * &
          ((Tc(i,jk+jABOVE) - Tb(i,jk+jABOVE)) * zrupd(i,jk+jBELOW) + & 
          Tb(i,jk+jABOVE) * zrup(i,jk+jBELOW))) * zreflect(i)
        zrupd(i,jk+jABOVE) = Rd(i,jk+jABOVE) + &
          Td(i,jk+jABOVE) * Td(i,jk+jABOVE) * zrupd(i,jk+jBELOW) * &
          zreflect(i)
      ENDDO
    ENDDO

    ! Upper boundary conditions
    DO i = 1, kproma
      ztdn (i,jTOA+jABOVE) = 1._wp
      zrdnd(i,jTOA+jABOVE) = 0._wp
      ztdn (i,JTOA+jBELOW) = Tc(i,jTOA+jABOVE)
      zrdnd(i,JTOA+jBELOW) = Rd(i,jTOA+jABOVE)
    ENDDO

    ! Pass from top to bottom
    !DO jk = klev,2,-1
    DO jk = jTOA-jINC, jSFC, -jINC
      DO i = 1, kproma
        zreflect(i) = 1._wp / (1._wp - Rd(i,jk+jABOVE) * zrdnd(i,jk+jABOVE))
        ztdn(i,jk+jBELOW) = accTb(i,jk+jABOVE) * Tc(i,jk+jABOVE) + &
          (Td(i,jk+jABOVE) * ((ztdn(i,jk+jABOVE) - accTb(i,jk+jABOVE)) + &
           accTb(i,jk+jABOVE) * Rc(i,jk+jABOVE) * zrdnd(i,jk+jABOVE))) * zreflect(i)
        zrdnd(i,jk+jBELOW) = Rd(i,jk+jABOVE) + &
           Td(i,jk+jABOVE) * Td(i,jk+jABOVE) * zrdnd(i,jk+jABOVE) * zreflect(i)
      ENDDO
    ENDDO

    ! Up and down-welling fluxes at levels
    DO jk = 1,klev+1
      DO i = 1, kproma
        zreflect(i) = 1._wp / (1._wp - zrdnd(i,jk) * zrupd(i,jk))
        flux_up(i,jk) = (accTb(i,jk) * zrup(i,jk) + &
          (ztdn(i,jk) - accTb(i,jk)) * zrupd(i,jk)) * zreflect(i)
        flux_dn(i,jk) = accTb(i,jk) + &
          (ztdn(i,jk) - accTb(i,jk) + &
            accTb(i,jk) * zrup(i,jk) * zrdnd(i,jk)) * zreflect(i)
      END DO
    END DO
  END SUBROUTINE adding

!! @brief Two-stream solutions to direct and diffuse reflectance and transmittance for a layer
!!    with optical depth tau, single scattering albedo w0, and asymmetery parameter g.  
!!
!! @ remarks Equations are developed in Meador and Weaver, 1980, 
!!    doi:10.1175/1520-0469(1980)037<0630:TSATRT>2.0.CO;2
!**** *SRTM_REFTRA* - REFLECTIVITY AND TRANSMISSIVITY
!           COMPUTES THE REFLECTIVITY AND TRANSMISSIVITY OF A CLEAR OR 
!     CLOUDY LAYER USING A CHOICE OF VARIOUS APPROXIMATIONS.
!        EXPLICIT ARGUMENTS :
!        --------------------
! INPUTS
!      KMODTS  = 1 EDDINGTON (JOSEPH ET AL., 1976)
!              = 2 PIFM (ZDUNKOWSKI ET AL., 1980)
!              = 3 DISCRETE ORDINATES (LIOU, 1973)
!      LDRTCHK = .T. IF CLOUDY
!              = .F. IF CLEAR-SKY
!      PGG     = ASSYMETRY FACTOR
!      PRMUZ   = COSINE SOLAR ZENITH ANGLE
!      PRMUZI  = INVERSE COSINE SOLAR ZENITH ANGLE
!      PTAU    = OPTICAL THICKNESS
!      PW      = SINGLE SCATTERING ALBEDO

! OUTPUTS
!      PREF    : COLLIMATED BEAM REFLECTIVITY
!      PREFD   : DIFFUSE BEAM REFLECTIVITY 
!      PTRA    : COLLIMATED BEAM TRANSMISSIVITY
!      PTRAD   : DIFFUSE BEAM TRANSMISSIVITY
!     METHOD.
!     -------
!          STANDARD DELTA-EDDINGTON, P.I.F.M., OR D.O.M. LAYER CALCULATIONS.
!     AUTHOR.
!        JEAN-JACQUES MORCRETTE  *ECMWF*
!     MODIFICATIONS.
!        ORIGINAL : 03-02-27
!        M.Hamrud   01-Oct-2003      CY28 Cleaning
!        Mike Iacono, AER, Mar 2004: bug fix 
!        D.Salmond  31-Oct-2007 Vector version in the style of RRTM from Meteo France & NEC
!        M. Puetz   20-Apr-2010 Gather/Scatter for better pipelining on scalar cpus
  SUBROUTINE srtm_reftra_ec(kproma, kbdim, klev, &
    ignore_mask, mask, &
    cos_mu0, inv_cos_mu0, tau, asymm, omega, Rc, Rd, Tc, Td, Tb)

    USE mo_psrad_general, ONLY: jABOVE

    IMPLICIT NONE

    INTEGER,INTENT(IN) :: kproma, kbdim, klev
    LOGICAL, INTENT(IN) :: ignore_mask, mask(KBDIM,klev) 
    REAL(wp), DIMENSION(KBDIM), INTENT(IN) :: cos_mu0, inv_cos_mu0
    REAL(wp), DIMENSION(KBDIM,klev), INTENT(INOUT) :: tau, asymm, omega
    REAL(wp), DIMENSION(KBDIM,klev+1), INTENT(INOUT) :: &
      Rc, Rd, Tc, Td, Tb

    INTEGER :: jk, ic, jc, ict, icc, icn, ns, pass
    INTEGER, DIMENSION(KBDIM) :: idxt, idxc, idxn, shuffle
    REAL(wp), DIMENSION(KBDIM) :: zgamma1, zgamma2, zgamma3, &
      zrk, zem1, zem2, zep1, zep2, za, za1, za2, zemm, zdend, &
      zgt, zr1, zr2, zr3, zr4, zr5, zrk2, zrkg, zrp, &
      zrpp, zt1, zt2, zt3, ztemp, zw, zg, zg3, &
      stau, scos, s_invcos, wf
    REAL(wp) :: zgamma4
    REAL(wp), PARAMETER :: replog=1e-12     ! epsilon for lapace transformation
    INTEGER,  PARAMETER :: kmodts=2

    REAL(wp) :: zsr3
    zsr3=SQRT(3._wp)

    IF (ignore_mask) THEN
      idxt(1:kproma) = (/(ic, ic = 1,kproma)/)
      ict = kproma
    ENDIF

    DO jk = 1, klev
      IF (.not. ignore_mask) THEN
        ict = 0
        DO ic = 1, kproma
          IF (mask(ic,jk)) THEN
            ict = ict + 1
            idxt(ict) = ic
          ENDIF
        ENDDO
      END IF
      icc = 0
      icn = 0
      DO ic = 1,ict
        zw(ic) = omega(idxt(ic),jk)
        zg(ic) = asymm(idxt(ic),jk)
        stau(ic) = tau(idxt(ic), jk)
        IF (zw(ic) >= w0Min) THEN
          icc = icc + 1
          idxc(icc) = idxt(ic)
        ELSE
          icn = icn + 1
          idxn(icn) = idxt(ic)
        ENDIF

        wf(ic) = zw(ic)* zg(ic)*zg(ic)
        stau(ic) = (1.0_wp - wf(ic)) * stau(ic)
        zw(ic) = (zw(ic) - wf(ic)) / (1.0_wp - wf(ic))
        zg(ic) = zg(ic) / (1.0_wp + zg(ic))

        omega(idxt(ic),jk) = zw(ic)
        asymm(idxt(ic),jk) = zg(ic)
        tau(idxt(ic),jk) = stau(ic)
      ENDDO

      shuffle(1:icc) = idxc(1:icc)
      ns=icc
      DO pass = 1,2
        DO ic = 1,ns
          zw(ic) = omega(shuffle(ic),jk)
          zg(ic) = asymm(shuffle(ic),jk)
          zg3(ic) = 3.0_wp * zg(ic)
          stau(ic) = tau(shuffle(ic),jk)
          s_invcos(ic) = inv_cos_mu0(shuffle(ic))
          scos(ic) = cos_mu0(shuffle(ic))
        ENDDO
        !-- GENERAL TWO-STREAM EXPRESSIONS
        IF (kmodts == 1) THEN
          DO ic = 1,ns
            zgamma1(ic)= (7._wp - zw(ic) * (4._wp + zg3(ic))) * 0.25_wp
            zgamma2(ic)=-(1._wp - zw(ic) * (4._wp - zg3(ic))) * 0.25_wp
            zgamma3(ic)= (2._wp - zg3(ic) * scos(ic) ) * 0.25_wp
          ENDDO
        ELSEIF (kmodts == 2) THEN  
          DO ic = 1,ns
            zgamma1(ic)= (8._wp - zw(ic) * (5._wp + zg3(ic))) * 0.25_wp
            zgamma2(ic)=  3._wp *(zw(ic) * (1._wp - zg(ic) )) * 0.25_wp
            zgamma3(ic)= (2._wp - zg3(ic) * scos(ic) ) * 0.25_wp
          ENDDO
        ELSEIF (kmodts == 3) THEN  
          DO ic = 1,ns
            zgamma1(ic)= zsr3 * (2._wp - zw(ic) * (1._wp + zg(ic))) * 0.5_wp
            zgamma2(ic)= zsr3 * zw(ic) * (1._wp - zg(ic) ) * 0.5_wp
            zgamma3(ic)= (1._wp - zsr3 * zg(ic) * scos(ic) ) * 0.5_wp
          ENDDO
        ENDIF

        IF (pass == 1) THEN
          !-- conservative scattering
          !-- Homogeneous reflectance and transmittance
!DIR$ NOVECTOR
          DO ic = 1,ns
            zem2(ic) = -MIN(stau(ic) * s_invcos(ic),500._wp)
            zem2(ic) = EXP(zem2(ic))
            za(ic)  = zgamma1(ic) * scos(ic) 
            za1(ic) = za(ic) - zgamma3(ic)
            zgt(ic) = zgamma1(ic) * stau(ic)
          ENDDO

          DO ic = 1,ns
            ! collimated beam
            ztemp(ic)=1.0_wp/(1._wp + zgt(ic))
            zrp(ic) = &
              (zgt(ic) - za1(ic) * (1._wp - zem2(ic))) * ztemp(ic)
            Rc(shuffle(ic),jk+jABOVE) = zrp(ic)
            Tc(shuffle(ic),jk+jABOVE) = 1._wp - zrp(ic)
          ENDDO

          DO ic = 1,ns
            ! isotropic incidence
            zrp(ic) = zgt(ic) * ztemp(ic)
            Rd(shuffle(ic),jk+jABOVE) = zrp(ic)
            Td(shuffle(ic),jk+jABOVE) = 1._wp - zrp(ic)
          ENDDO
        ELSE
          !-- non-conservative scattering
          !-- Homogeneous reflectance and transmittance
          DO ic = 1,ns
            zrk(ic) = SQRT(MAX(zgamma1(ic)**2 - zgamma2(ic)**2,replog))
          ENDDO

          DO ic = 1,ns
            zep1(ic) = MIN(zrk(ic) * stau(ic), 500._wp)
            zep1(ic) = EXP(zep1(ic))
            zem1(ic) = 1.0_wp/zep1(ic)
          ENDDO

          DO ic = 1,ns
            zep2(ic) = MIN(stau(ic) * s_invcos(ic),500._wp)
            zep2(ic) = EXP(zep2(ic))
            zem2(ic) = 1.0_wp/zep2(ic)
          ENDDO


          DO ic = 1,ns
            zgamma4 = 1._wp - zgamma3(ic)
            za1(ic) = zgamma1(ic) * zgamma4 + &
              zgamma2(ic) * zgamma3(ic)
            za2(ic) = zgamma1(ic) * zgamma3(ic) + &
              zgamma2(ic) * zgamma4
          ENDDO

          DO ic = 1,ns
            zrp(ic) = zrk(ic) * scos(ic)
            zrk2(ic) = 2._wp * zrk(ic)
            zrpp(ic) = 1._wp - zrp(ic)*zrp(ic)
            zrkg(ic) = zrk(ic) + zgamma1(ic)
          ENDDO
          DO ic = 1,ns
            zr1(ic) = (1._wp - zrp(ic)) * (za2(ic) + zrk(ic) * zgamma3(ic))
            zr2(ic) = (1._wp + zrp(ic)) * (za2(ic) - zrk(ic) * zgamma3(ic))
            zr3(ic) = zrk2(ic) * (zgamma3(ic) - za2(ic) * scos(ic))
            zr4(ic) = zrpp(ic) * zrkg(ic)
            zr5(ic) = zrpp(ic) * (zrk(ic) - zgamma1(ic))
          ENDDO
          DO ic = 1,ns
            zgamma4 = 1._wp - zgamma3(ic)
            zt1(ic) = (1._wp + zrp(ic)) * (za1(ic) + zrk(ic) * zgamma4)
            zt2(ic) = (1._wp - zrp(ic)) * (za1(ic) - zrk(ic) * zgamma4)
            zt3(ic) = zrk2(ic) * (zgamma4 + za1(ic) * scos(ic))
            !zt4(ic) = zr4(ic)
            !zt5(ic) = zr5(ic)
          ENDDO

          ! collimated beam
          DO ic = 1,ns
            jc = shuffle(ic)
            Rc(jc,jk+jABOVE) = zw(ic)  * (zr1(ic)*zep1(ic) - &
              zr2(ic)*zem1(ic) - zr3(ic)*zem2(ic)) / &
              (zr4(ic)*zep1(ic) + zr5(ic)*zem1(ic))
              
            Tc(jc,jk+jABOVE) = &
              zem2(ic) * (1._wp - zw(ic) * (zt1(ic)*zep1(ic) - &
              zt2(ic)*zem1(ic) - zt3(ic)*zep2(ic)) / &
              (zr4(ic)*zep1(ic) + zr5(ic)*zem1(ic)))
          ENDDO

          ! diffuse beam
          DO ic = 1,ns
            zemm(ic) = zem1(ic)*zem1(ic)
            zdend(ic) = 1._wp / (zrkg(ic) + (zrk(ic) - zgamma1(ic)) * zemm(ic))
            jc = shuffle(ic)
            Rd(jc,jk+jABOVE) =  zgamma2(ic) * (1._wp - zemm(ic)) * zdend(ic)
            Td(jc,jk+jABOVE) =  zrk2(ic)*zem1(ic)*zdend(ic)
          ENDDO
        ENDIF
        DO ic = 1,ns
          Tb(shuffle(ic),jk+jABOVE) = zem2(ic)
        ENDDO

        shuffle(1:icn) = idxn(1:icn)
        ns=icn
      ENDDO
    ENDDO
  END SUBROUTINE srtm_reftra_ec

END MODULE mo_psrad_srtm_solver
