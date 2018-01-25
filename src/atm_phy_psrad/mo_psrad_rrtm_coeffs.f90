!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
!! @brief Module to provide coefficients for use in longwave solver
!!
!! @remarks
!!   For a given atmosphere, calculate the indices and fractions related to the 
!!   pressure and temperature interpolations. Also calculate the values of the 
!!   integrated Planck functions for each band at the level and layer temperatures.
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
MODULE mo_psrad_rrtm_coeffs

  USE mo_psrad_general, ONLY : wp, preflog, tref, &
    ngas, ih2o, ico2, ich4, io3, in2o, ico, &
    nreact, ih2oco2, ih2oo3, ih2on2o, ih2och4, in2oco2, io3co2
  USE mo_psrad_lrtm_kgs, ONLY : chi_mls
  IMPLICIT NONE

  PRIVATE

  INTEGER, PARAMETER :: reaction(3,nreact) = RESHAPE((/ &
    ih2oco2, ih2o, ico2, &
    ih2oo3, ih2o, io3, & ! Needed only in lower atmos (plog > 4.56_wp) 
    ih2on2o, ih2o, in2o, &
    ih2och4, ih2o, ich4, &
    in2oco2, in2o, ico2, &
    io3co2, io3, ico2/), & ! Needed only in upper atmos (plog <= 4.56_wp) 
    SHAPE=(/3,nreact/))

  PUBLIC :: lrtm_coeffs, srtm_coeffs
  
  REAL(wp), PARAMETER :: stpfac  = 296._wp/1013._wp

CONTAINS
! --------------------------------------------------------------------------------------------
  SUBROUTINE lrtm_coeffs(jcs, kproma, kbdim, klev, play, tlay, coldry, wkl, &
    wbroad, laytrop, jp, jp1, jt, jt1, iabs, gases, colbrd, fac, ratio, &
    h2o_factor, h2o_fraction, h2o_index, &
    minorfrac, scaleminor, scaleminorn2, indminor)

    use mo_psrad_general, only: ngas

    INTEGER, INTENT(in) :: jcs, kproma, kbdim, klev
    REAL(wp), INTENT(in) :: &
      play(KBDIM,klev), & ! layer pressures (mb) 
      tlay(KBDIM,klev), & ! layer temperatures (K)
      coldry(KBDIM,klev), & ! dry air column density (mol/cm2)
      wbroad(KBDIM,klev), & ! broadening gas column density (mol/cm2)
      wkl(KBDIM,klev,ngas) !< molecular amounts (mol/cm-2) (klev,ngas)

    INTEGER, INTENT(out) :: laytrop(KBDIM), & !< tropopause layer index
      iabs(KBDIM,2,2,klev)
    INTEGER, DIMENSION(KBDIM,klev), INTENT(out) :: jp, jp1, jt, jt1, &
      indminor
    REAL(wp), DIMENSION(KBDIM,klev), INTENT(out) :: colbrd, &
      minorfrac, scaleminor, &
      scaleminorn2
    REAL(wp), INTENT(out) :: fac(KBDIM,2,2,klev), gases(KBDIM,klev,ngas), &
      ratio(KBDIM,2,klev,nreact)
    REAL(wp), DIMENSION(KBDIM,klev,2), INTENT(out) :: h2o_factor,h2o_fraction
    INTEGER, DIMENSION(KBDIM,klev,2), INTENT(out) :: h2o_index
    
    INTEGER  :: r, jk
    REAL(wp) :: colmol(KBDIM,klev), factor(KBDIM,klev) 

    CALL srtm_coeffs(jcs, kproma, KBDIM, klev, play, tlay, coldry, wkl, laytrop, &
      jp, jp1, jt, jt1, iabs, gases, colmol, fac, &
      h2o_factor(:,:,1), h2o_fraction(:,:,1), h2o_index(:,:,1), &
      h2o_factor(:,:,2), h2o_fraction(:,:,2), h2o_index(:,:,2))
    
    colbrd(jcs:kproma,:) = 1.e-20_wp * wbroad(jcs:kproma,:)
    gases(jcs:kproma,:,ico) = MERGE(1.e-20_wp * wkl(jcs:kproma,:,ico),  &
                                    1.e-32_wp * coldry(jcs:kproma,:), &
                                    wkl(jcs:kproma,:,ico) > 0._wp)
    
    ! Water vapor continuum broadening factors are used differently 
    ! in LW and SW? 
    h2o_factor(jcs:kproma,:,1) = h2o_factor(jcs:kproma,:,1) * gases(jcs:kproma,:,ih2o)
    h2o_factor(jcs:kproma,:,2) = h2o_factor(jcs:kproma,:,2) * gases(jcs:kproma,:,ih2o)
    
    !  Setup reference ratio to be used in calculation of binary species 
    ! parameter.
    ratio = 1.e20
    DO r = 1,nreact
      DO jk = 1,klev
        WHERE (chi_mls(reaction(3,r), jp(jcs:kproma,jk)) /= 0)
          ratio(jcs:kproma,1,jk,reaction(1,r)) = chi_mls(reaction(2,r), jp(jcs:kproma,jk)) /&
            chi_mls(reaction(3,r), jp(jcs:kproma,jk))
        ENDWHERE
        WHERE (chi_mls(reaction(3,r), jp1(jcs:kproma,jk)) /= 0)
          ratio(jcs:kproma,2,jk,reaction(1,r)) = chi_mls(reaction(2,r), jp1(jcs:kproma,jk)) /&
            chi_mls(reaction(3,r), jp1(jcs:kproma,jk))
        ENDWHERE
      END DO
    END DO

    !  Set up factors needed to separately include the minor gases
    !  in the calculation of absorption coefficient
    scaleminor(jcs:kproma,:) = play(jcs:kproma,:) / &
      tlay(jcs:kproma,:)
    scaleminorn2(jcs:kproma,:) = scaleminor(jcs:kproma,:) * &
      (wbroad(jcs:kproma,:) / &
      (coldry(jcs:kproma,:)+wkl(jcs:kproma,:,ih2o)))
    factor(jcs:kproma,:) = (tlay(jcs:kproma,:)-180.8_wp) / 7.2_wp
    indminor(jcs:kproma,:) = MIN(18, MAX(1, INT(factor(jcs:kproma,:))))
    minorfrac(jcs:kproma,:) = (tlay(jcs:kproma,:)-180.8_wp)/7.2_wp - &
      FLOAT(indminor(jcs:kproma,:))

  END SUBROUTINE lrtm_coeffs    

  SUBROUTINE srtm_coeffs(jcs, kproma, kbdim, klev, play, tlay, coldry, wkl, &
    laytrop, jp, jp1, jt, jt1, iabs, gases, colmol, fac, selffac, selffrac, &
    indself, forfac, forfrac, indfor)
    INTEGER, INTENT(in) :: jcs, kproma, kbdim, klev
    REAL(wp), INTENT(in) :: &
      play(KBDIM,klev), & ! layer pressures (mb) 
      tlay(KBDIM,klev), & ! layer temperatures (K)
      coldry(KBDIM,klev), & ! dry air column density (mol/cm2)
      wkl(KBDIM,klev,ngas) ! molecular amounts (mol/cm-2) 

    ! Output Dimensions kproma, klev unless otherwise specified
    INTEGER, INTENT(out) :: laytrop(KBDIM), & !< tropopause layer index
      iabs(KBDIM,2,2,klev)
    INTEGER, DIMENSION(KBDIM,klev), INTENT(out) :: jp, jp1, jt, jt1, &
      indself, indfor
    REAL(wp), DIMENSION(KBDIM,klev), INTENT(out) :: colmol, selffac, &
      selffrac, forfac, forfrac
    REAL(wp), INTENT(out) :: fac(KBDIM,2,2,klev), gases(KBDIM,klev,ngas)
    INTEGER :: jk
    REAL(wp), DIMENSION(KBDIM,klev) :: plog, ft, ft1, water, scalefac, &
      factor 
    REAL(wp), DIMENSION(KBDIM) :: fp
    !  Find the two reference pressures on either side of the
    !  layer pressure.  Store them in JP and JP1.  Store in FP the
    !  fraction of the difference (in ln(pressure)) between these
    !  two values that the layer pressure lies.
    plog(jcs:kproma,:) = LOG(play(jcs:kproma,:))
    jp (jcs:kproma,:) = MIN(58, &
      MAX(1,INT(36._wp - 5*(plog(jcs:kproma,:)+0.04_wp))))
    jp1(jcs:kproma,:) = jp(jcs:kproma,:) + 1
    !  Determine, for each reference pressure (JP and JP1), which
    !  reference temperature (these are different for each  
    !  reference pressure) is nearest the layer temperature but does
    !  not exceed it.  Store these indices in JT and JT1, resp.
    !  Store in FT (resp. FT1) the fraction of the way between JT
    !  (JT1) and the next highest reference temperature that the 
    !  layer temperature falls.
    DO jk = 1, klev
      jt(jcs:kproma,jk) = MIN(4,&
        MAX(1,INT(3._wp + (tlay(jcs:kproma,jk) - tref(jp(jcs:kproma,jk)))/15._wp)))
      jt1(jcs:kproma,jk) = MIN(4,&
        MAX(1,INT(3._wp + (tlay(jcs:kproma,jk) - tref(jp1(jcs:kproma,jk)))/15._wp)))
    END DO 
    DO jk = 1, klev
      ft(jcs:kproma,jk) = ((tlay(jcs:kproma,jk)-tref(jp(jcs:kproma,jk)))/15._wp) - &
        float(jt(jcs:kproma,jk)-3)
      ft1(jcs:kproma,jk) = ((tlay(jcs:kproma,jk)-tref(jp1(jcs:kproma,jk)))/15._wp) - &
        float(jt1(jcs:kproma,jk)-3)
    END DO 

    iabs(jcs:kproma,1,1,:) = MIN(11,jp(jcs:kproma,:)-1)*5+(jt(jcs:kproma,:)-1)
    iabs(jcs:kproma,2,1,:) = MIN(12,jp(jcs:kproma,:))*5+(jt1(jcs:kproma,:)-1)
    iabs(jcs:kproma,1,2,:) = MAX(0,jp(jcs:kproma,:)-13)*5+(jt(jcs:kproma,:)-1)
    iabs(jcs:kproma,2,2,:) = MAX(1,jp(jcs:kproma,:)-12)*5+(jt1(jcs:kproma,:)-1)

    water(jcs:kproma,:) = wkl(jcs:kproma,:,ih2o)/coldry(jcs:kproma,:)
    scalefac(jcs:kproma,:) = play(jcs:kproma,:) * stpfac / &
      tlay(jcs:kproma,:)

    !  We have now isolated the layer ln pressure and temperature,
    !  between two reference pressures and two reference temperatures 
    !  (for each reference pressure).  We multiply the pressure 
    !  fraction FP with the appropriate temperature fractions to get 
    !  the factors that will be needed for the interpolation that yields
    !  the optical depths (performed in routines TAUGBn for band n).`
    !
    do jk = 1, klev
      fp(jcs:kproma) = 5._wp *(preflog(jp(jcs:kproma,jk)) - plog(jcs:kproma,jk))
      fac(jcs:kproma,2,2,jk) = fp(jcs:kproma) * ft1(jcs:kproma,jk)
      fac(jcs:kproma,1,2,jk) = fp(jcs:kproma) * (1._wp - ft1(jcs:kproma,jk))
      fp(jcs:kproma) = 1. - fp(jcs:kproma)
      fac(jcs:kproma,2,1,jk) = fp(jcs:kproma) * ft(jcs:kproma,jk)
      fac(jcs:kproma,1,1,jk) = fp(jcs:kproma) * (1._wp - ft(jcs:kproma,jk))
    end do 
    
    ! Tropopause defined in terms of pressure (~100 hPa)
    ! We're looking for the first layer (counted from the bottom) at which 
    ! the pressure reaches or falls below this value
    laytrop(jcs:kproma) = COUNT(plog(jcs:kproma,:) > 4.56_wp, DIM = 2) 

    ! Calculate needed column amounts. Only a few ratios are used in the 
    ! upper atmosphere but masking may be less efficient
    ! BUG: if size mismatch cannot be determined at compile time,
    ! no exception is thrown at runtime either!
    gases(jcs:kproma,:,ih2o:in2o) = MERGE(1.e-20_wp * wkl(jcs:kproma,:,ih2o:in2o),  &
      SPREAD(1.e-32_wp * coldry(jcs:kproma,:),3,in2o-ih2o+1), &
      wkl(jcs:kproma,:,ih2o:in2o) > 0._wp)
    colmol(jcs:kproma,:) = 1.e-20_wp * coldry(jcs:kproma,:) + gases(jcs:kproma,:,ih2o)

    ! Interpolation coefficients 
    forfac(jcs:kproma,:) = scalefac(jcs:kproma,:) / &
      (1._wp+water(jcs:kproma,:))
    ! Set up factors needed to separately include the water vapor
    ! self-continuum in the calculation of absorption coefficient.
    selffac(jcs:kproma,:)  = water(jcs:kproma,:) * &
      forfac(jcs:kproma,:)

    ! If the pressure is less than ~100mb, perform a different set of species
    ! interpolations.
    factor(jcs:kproma,:) = (332.0_wp-tlay(jcs:kproma,:)) / 36.0_wp
    indfor(jcs:kproma,:) = MERGE(3, MIN(2, &
      MAX(1, INT(factor(jcs:kproma,:)))), plog(jcs:kproma,:) <= 4.56_wp)
            
    forfrac(jcs:kproma,:) = MERGE( &
      (tlay(jcs:kproma,:)-188.0_wp)/36.0_wp - 1.0_wp, &
      factor(jcs:kproma,:) - FLOAT(indfor(jcs:kproma,:)), &
      plog(jcs:kproma,:) <= 4.56_wp)

    ! In RRTMG code, this calculation is done only in the lower atmosphere 
    ! (plog > 4.56) 
    factor(jcs:kproma,:) = (tlay(jcs:kproma,:)-188.0_wp)/7.2_wp
    indself(jcs:kproma,:) = MIN(9, MAX(1, INT(factor(jcs:kproma,:))-7))
    selffrac(jcs:kproma,:) = factor(jcs:kproma,:) - &
      float(indself(jcs:kproma,:) + 7)
  END SUBROUTINE srtm_coeffs
END MODULE mo_psrad_rrtm_coeffs

