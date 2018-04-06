!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!

! Bands 16->29 renumbered 1->14
! Band (gpt_range): wave number range (low key species; high key species)
! 16/1 (1-6): 2600-3250 cm-1 (low - h2o,ch4; high - ch4)
! 17/2 (7-18): 3250-4000 cm-1 (low - h2o,co2; high - h2o,co2)
! 18/3 (19-26): 4000-4650 cm-1 (low - h2o,ch4; high - ch4)
! 19/4 (27-34): 4650-5150 cm-1 (low - h2o,co2; high - co2)
! 20/5 (35-44): 5150-6150 cm-1 (low - h2o; high - h2o)
! 21/6 (45-54): 6150-7700 cm-1 (low - h2o,co2; high - h2o,co2)
! 22/7 (55-56): 6150-7700 cm-1 (low - h2o,co2; high - h2o,co2)
! In this band the ratio of total O2 band intensity (lines 
! and Mate continuum) to O2 band intensity (line only) is 1.6 and is used
! to adjust the optical depths since the k's include only lines.  This is
! done by multiplying swght1 by 1.6 and abs_ab(:,2,:) by 1.6 to account for 
! this difference in below and above 100 hPa respectively.  Also note that
! the minor gas does not have a g-point dependent absorption because the
! o2 abosrption is a continuum effect

! 23/8 (57-66): 8050-12850 cm-1 (low - h2o; high - nothing)
!  Average Giver et al. correction factor for this band is 1.029
!  and is implemented by multiplying the abosrption coefficients

! 24/9 (67-74): 12850-16000 cm-1 (low - h2o,o2; high - o2)
! 25/10 (75-80): 16000-22650 cm-1 (low - h2o; high - nothing)
! 26/11 (81-86): 22650-29000 cm-1 (low - nothing; high - nothing)
! 27/12 (87-94):  29000-38000 cm-1 (low - o3; high - o3)
! Kurucz solar source function
! The values in sfluxref were obtained using the "low resolution"
! version of the Kurucz solar source function.  For unknown reasons,
! the total irradiance in this band differs from the corresponding
! total in the "high-resolution" version of the Kurucz function.
! Therefore, these values are scaled below by the factor SCALEKUR.

! 28/13 (95-100): 38000-50000 cm-1 (low - o3,o2; high - o3,o2)
! 29/14 (101-112): 820-2600 cm-1 (low - h2o; high - co2)

MODULE mo_psrad_srtm_kgs

  USE mo_psrad_general

  IMPLICIT NONE

  PUBLIC

  INTEGER, PARAMETER :: &
    ngpt(nbndsw) = (/6,12,8,8,10, 10,2,10,8,6, 6,8,6,12/), &
    nsp(2,nbndsw) = RESHAPE((/& 
      9,9,9,9,1, 9,9,1,9,1, 0,1,9,1, &
      1,5,1,1,1, 5,1,0,1,0, 0,1,5,1/), &
      SHAPE=(/2,nbndsw/), ORDER=(/2,1/)), &
    fracs_mult(2,nbndsw) = MAX(0,nsp-1), &
    minor_species(2,nbndsw) = RESHAPE((/ &
      0,0,0,0,ich4, 0,io2,0,io3,io3, 0,0,0,ico2, &
      0,0,0,0,ich4, 0,io2,0,io3,io3, 0,0,0,ih2o/), &
      SHAPE=(/2,nbndsw/), ORDER=(/2,1/)), &
    major_species(2,2,nbndsw) = RESHAPE((/& 
      ih2o,ich4, ich4,0, & ! Band 1
      ih2o,ico2, ih2o,ico2, & ! Band 2
      ih2o,ich4, ich4,0, & ! Band 3
      ih2o,ico2, ico2,0, & ! Band 4
      ih2o,0, ih2o,0, & ! Band 5
      ih2o,ico2, ih2o,ico2, & ! Band 6
      ih2o,io2, io2,0, & ! Band 7
      ih2o,0, 0,0, & ! Band 8
      ih2o,io2, io2,0, & ! Band 9
      ih2o,0, 0,0, & ! Band 10
      0,0, 0,0, & ! Band 11
      io3,0, io3,0, & ! Band 12
      io3,io2, io3,io2, & ! Band 13
      ih2o,0, ico2,0/), & ! Band 14
      SHAPE=(/2,2,nbndsw/)), &
    h2o_absorption_flag(2,nbndsw) = RESHAPE((/& 
      1,1,1,1,1, 1,1,1,1,0, 0,0,0,1, &
      0,1,0,0,1, 1,0,0,0,0, 0,0,0,0/), &
      SHAPE=(/2,nbndsw/), ORDER=(/2,1/)), &
    nh2oref(2,nbndsw) = RESHAPE((/& 
      10,10,10,10,10, 10,10,10,10,0, 0,0,0,10, &
      3,4,3,3,4,      4,3,3,3,0,     0,0,0,4/), &
      SHAPE=(/2,nbndsw/), ORDER=(/2,1/)), &
    nsfluxref(nbndsw) = (/1,5,9,9,1, 9,9,1,9,1, 1,1,5,1/), &
    rayl_type(2,nbndsw) = RESHAPE((/ &
      0,0,0,0,0, 0,0,1,9,1, 1,1,0,0, &
      0,0,0,0,0, 0,0,1,1,1, 1,1,0,0/), &
      SHAPE=(/2,nbndsw/), ORDER=(/2,1/)), &
    minor_species_missing_data(2,nbndsw) = RESHAPE((/ &
      0,0,0,0,0, 0,0,0,0,0, 0,0,0,0, &
      0,0,0,0,1, 0,0,0,0,0, 0,0,0,0/), &
      SHAPE=(/2,nbndsw/), ORDER=(/2,1/))


  ! Shortwave spectral band limits (wavenumbers)
  REAL(wp), PARAMETER :: wavenum1(nbndsw) = (/ &
       2600._wp, 3250._wp, 4000._wp, 4650._wp, 5150._wp, 6150._wp, 7700._wp, &
       8050._wp,12850._wp,16000._wp,22650._wp,29000._wp,38000._wp,  820._wp/)
  REAL(wp), PARAMETER :: wavenum2(nbndsw) = (/ &
       3250., 4000., 4650., 5150., 6150., 7700., 8050., &
       12850.,16000.,22650.,29000.,38000.,50000., 2600./), &
    delwave(nbndsw)  = (/ &
       650.,  750.,  650.,  500., 1000., 1550.,  350., &
       4800., 3150., 6650., 6350., 9000.,12000., 1780./)

END MODULE mo_psrad_srtm_kgs

