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
  SUBROUTINE lrtm_coeffs(kbdim, klev, play, tlay, coldry, wkl, &
    wbroad, laytrop, jp, jp1, jt, jt1, iabs, gases, colbrd, fac, ratio, &
    h2o_factor, h2o_fraction, h2o_index, &
    minorfrac, scaleminor, scaleminorn2, indminor)

    use mo_psrad_general, only: ngas

    INTEGER, INTENT(in) :: kbdim, klev
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

    CALL srtm_coeffs(KBDIM, klev, play, tlay, coldry, wkl, laytrop, &
      jp, jp1, jt, jt1, iabs, gases, colmol, fac, &
      h2o_factor(:,:,1), h2o_fraction(:,:,1), h2o_index(:,:,1), &
      h2o_factor(:,:,2), h2o_fraction(:,:,2), h2o_index(:,:,2))
    
    colbrd(:,:) = 1.e-20_wp * wbroad(:,:)
    gases(:,:,ico) = MERGE(1.e-20_wp * wkl(:,:,ico),  &
                                    1.e-32_wp * coldry(:,:), &
                                    wkl(:,:,ico) > 0._wp)
    
    ! Water vapor continuum broadening factors are used differently 
    ! in LW and SW? 
    h2o_factor(:,:,1) = h2o_factor(:,:,1) * gases(:,:,ih2o)
    h2o_factor(:,:,2) = h2o_factor(:,:,2) * gases(:,:,ih2o)
    
    !  Setup reference ratio to be used in calculation of binary species 
    ! parameter.
    ratio = 1e20
    DO r = 1,nreact
      DO jk = 1,klev
        WHERE (chi_mls(reaction(3,r), jp(:,jk)) /= 0)
          ratio(:,1,jk,reaction(1,r)) = chi_mls(reaction(2,r), jp(:,jk)) /&
            chi_mls(reaction(3,r), jp(:,jk))
        ENDWHERE
        WHERE (chi_mls(reaction(3,r), jp1(:,jk)) /= 0)
          ratio(:,2,jk,reaction(1,r)) = chi_mls(reaction(2,r), jp1(:,jk)) /&
            chi_mls(reaction(3,r), jp1(:,jk))
        ENDWHERE
      END DO
    END DO

    !  Set up factors needed to separately include the minor gases
    !  in the calculation of absorption coefficient
    scaleminor(:,:) = play(:,:) / &
      tlay(:,:)
    scaleminorn2(:,:) = scaleminor(:,:) * &
      (wbroad(:,:) / &
      (coldry(:,:)+wkl(:,:,ih2o)))
    factor(:,:) = (tlay(:,:)-180.8_wp) / 7.2_wp
    indminor(:,:) = MIN(18, MAX(1, INT(factor(:,:))))
    minorfrac(:,:) = (tlay(:,:)-180.8_wp)/7.2_wp - &
      FLOAT(indminor(:,:))

  END SUBROUTINE lrtm_coeffs    

  SUBROUTINE srtm_coeffs(kbdim, klev, play, tlay, coldry, wkl, &
    laytrop, jp, jp1, jt, jt1, iabs, gases, colmol, fac, selffac, selffrac, &
    indself, forfac, forfrac, indfor)
    INTEGER, INTENT(in) :: kbdim, klev
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
    plog(:,:) = LOG(play(:,:))
    jp (:,:) = MIN(58, &
      MAX(1,INT(36._wp - 5*(plog(:,:)+0.04_wp))))
    jp1(:,:) = jp(:,:) + 1
    !  Determine, for each reference pressure (JP and JP1), which
    !  reference temperature (these are different for each  
    !  reference pressure) is nearest the layer temperature but does
    !  not exceed it.  Store these indices in JT and JT1, resp.
    !  Store in FT (resp. FT1) the fraction of the way between JT
    !  (JT1) and the next highest reference temperature that the 
    !  layer temperature falls.
    DO jk = 1, klev
      jt(:,jk) = MIN(4,&
        MAX(1,INT(3._wp + (tlay(:,jk) - tref(jp(:,jk)))/15._wp)))
      jt1(:,jk) = MIN(4,&
        MAX(1,INT(3._wp + (tlay(:,jk) - tref(jp1(:,jk)))/15._wp)))
    END DO 
    DO jk = 1, klev
      ft(:,jk) = ((tlay(:,jk)-tref(jp(:,jk)))/15._wp) - &
        float(jt(:,jk)-3)
      ft1(:,jk) = ((tlay(:,jk)-tref(jp1(:,jk)))/15._wp) - &
        float(jt1(:,jk)-3)
    END DO 

    iabs(:,1,1,:) = MIN(11,jp(:,:)-1)*5+(jt(:,:)-1)
    iabs(:,2,1,:) = MIN(12,jp(:,:))*5+(jt1(:,:)-1)
    iabs(:,1,2,:) = MAX(0,jp(:,:)-13)*5+(jt(:,:)-1)
    iabs(:,2,2,:) = MAX(1,jp(:,:)-12)*5+(jt1(:,:)-1)

    water(:,:) = wkl(:,:,ih2o)/coldry(:,:)
    scalefac(:,:) = play(:,:) * stpfac / &
      tlay(:,:)

    !  We have now isolated the layer ln pressure and temperature,
    !  between two reference pressures and two reference temperatures 
    !  (for each reference pressure).  We multiply the pressure 
    !  fraction FP with the appropriate temperature fractions to get 
    !  the factors that will be needed for the interpolation that yields
    !  the optical depths (performed in routines TAUGBn for band n).`
    !
    do jk = 1, klev
      fp(:) = 5._wp *(preflog(jp(:,jk)) - plog(:,jk))
      fac(:,2,2,jk) = fp(:) * ft1(:,jk)
      fac(:,1,2,jk) = fp(:) * (1._wp - ft1(:,jk))
      fp(:) = 1. - fp(:)
      fac(:,2,1,jk) = fp(:) * ft(:,jk)
      fac(:,1,1,jk) = fp(:) * (1._wp - ft(:,jk))
    end do 
    
    ! Tropopause defined in terms of pressure (~100 hPa)
    ! We're looking for the first layer (counted from the bottom) at which 
    ! the pressure reaches or falls below this value
    laytrop(:) = COUNT(plog(:,:) > 4.56_wp, DIM = 2) 

    ! Calculate needed column amounts. Only a few ratios are used in the 
    ! upper atmosphere but masking may be less efficient
    ! BUG: if size mismatch cannot be determined at compile time,
    ! no exception is thrown at runtime either!
    gases(:,:,ih2o:in2o) = MERGE(1.e-20_wp * wkl(:,:,ih2o:in2o),  &
      SPREAD(1.e-32_wp * coldry(:,:),3,in2o-ih2o+1), &
      wkl(:,:,ih2o:in2o) > 0._wp)
    colmol(:,:) = 1.e-20_wp * coldry(:,:) + gases(:,:,ih2o)

    ! Interpolation coefficients 
    forfac(:,:) = scalefac(:,:) / &
      (1._wp+water(:,:))
    ! Set up factors needed to separately include the water vapor
    ! self-continuum in the calculation of absorption coefficient.
    selffac(:,:)  = water(:,:) * &
      forfac(:,:)

    ! If the pressure is less than ~100mb, perform a different set of species
    ! interpolations.
    factor(:,:) = (332.0_wp-tlay(:,:)) / 36.0_wp
    indfor(:,:) = MERGE(3, MIN(2, &
      MAX(1, INT(factor(:,:)))), plog(:,:) <= 4.56_wp)
            
    forfrac(:,:) = MERGE( &
      (tlay(:,:)-188.0_wp)/36.0_wp - 1.0_wp, &
      factor(:,:) - FLOAT(indfor(:,:)), &
      plog(:,:) <= 4.56_wp)

    ! In RRTMG code, this calculation is done only in the lower atmosphere 
    ! (plog > 4.56) 
    factor(:,:) = (tlay(:,:)-188.0_wp)/7.2_wp
    indself(:,:) = MIN(9, MAX(1, INT(factor(:,:))-7))
    selffrac(:,:) = factor(:,:) - &
      float(indself(:,:) + 7)
  END SUBROUTINE srtm_coeffs
END MODULE mo_psrad_rrtm_coeffs

