!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!

! Band (gpt_range): wave number range (low key species; high key species)
! 1(1 - 10): 10-250 cm-1 (low - h2o; high - h2o)
! 2(11 - 22): 250-500 cm-1 (low - h2o; high - h2o) - 
! 3(23 - 38): 500-630 cm-1 (low - h2o,co2; high - h2o,co2)
! 4(39 - 52): 630-700 cm-1 (low - h2o,co2; high - o3,co2)
! 5(53 - 68): 700-820 cm-1 (low - h2o,co2; high - o3,co2)
! 6(69 - 76): 820-980 cm-1 (low - h2o; high - nothing)
! 7(77 - 88): 980-1080 cm-1 (low - h2o,o3; high - o3)
! 8(89 - 96): 1080-1180 cm-1 (low (i.e.>~300mb) - h2o; high - o3)
! 9(97 - 108): 1180-1390 cm-1 (low - h2o,ch4; high - ch4)
! 10(109 - 114): 1390-1480 cm-1 (low - h2o; high - h2o)
! 11(115 - 122): 1480-1800 cm-1 (low - h2o; high - h2o)
! 12(123 - 130): 1800-2080 cm-1 (low - h2o,co2; high - nothing)
! 13(131 - 134): 2080-2250 cm-1 (low - h2o,n2o; high - nothing)
! 14(135 - 136): 2250-2380 cm-1 (low - co2; high - co2)
! 15(137 - 138): 2380-2600 cm-1 (low - n2o,co2; high - nothing)
! 16(139 - 140): 2600-3000 cm-1 (low - h2o,ch4; high - nothing)

! band 1:  10-350 cm-1 (low key - h2o; low minor - n2)
!                      (high key - h2o; high minor - n2)
!   old:   10-250 cm-1 (low - h2o; high - h2o)
! Minor gas mapping levels:
!     lower - n2, p = 142.5490 mbar, t = 215.70 k
!     upper - n2, p = 142.5490 mbar, t = 215.70 k
! band 2:  350-500 cm-1 (low key - h2o; high key - h2o)
!   old:   250-500 cm-1 (low - h2o; high - h2o)
! band 3:  500-630 cm-1 (low key - h2o,co2; low minor - n2o)
!                       (high key - h2o,co2; high minor - n2o)
!    old:  500-630 cm-1 (low - h2o,co2; high - h2o,co2)
! Minor gas mapping levels:
!     lower - n2o, p = 706.272 mbar, t = 278.94 k
!     upper - n2o, p = 95.58 mbar, t = 215.7 k
! band 4:  630-700 cm-1 (low key - h2o,co2; high key - o3,co2)
!    old:  630-700 cm-1 (low - h2o,co2; high - o3,co2)
! band 5:  700-820 cm-1 (low key - h2o,co2; low minor - o3, ccl4)
!                       (high key - o3,co2)
!    old:  700-820 cm-1 (low - h2o,co2; high - o3,co2)
! Minor gas mapping level :
!     lower - o3, p = 317.34 mbar, t = 240.77 k
!     lower - ccl4
!band 6:  820-980 cm-1 (low key - h2o; low minor - co2)
!                      (high key - nothing; high minor - cfc11, cfc12)
!   old:  820-980 cm-1 (low - h2o; high - nothing)
! NOTE: code has cfc11/12 in lower atmosphere as well.
! Minor gas mapping level:
!     lower - co2, p = 706.2720 mb, t = 294.2 k
!     upper - cfc11, cfc12
! band 7:  98051080 cm-1 (low key - h2o,o3; low minor - co2)
!                        (high key - o3; high minor - co2)
!    old:  980-1080 cm-1 (low - h2o,o3; high - o3)
! Minor gas mapping level :
!     lower - co2, p = 706.2620 mbar, t= 278.94 k
!     upper - co2, p = 12.9350 mbar, t = 234.01 k
! band 8:  1080-1180 cm-1 (low key - h2o; low minor - co2,o3,n2o)
!                         (high key - o3; high minor - co2, n2o)
! NOTE: Code contained cfc12 and cfc22 as well!!!
!    old:  1080-1180 cm-1 (low (i.e.>~300mb) - h2o; high - o3)
! Minor gas mapping level:
!     lower - co2, p = 1053.63 mb, t = 294.2 k
!     lower - o3,  p = 317.348 mb, t = 240.77 k
!     lower - n2o, p = 706.2720 mb, t= 278.94 k
!     lower - cfc12,cfc11
!     upper - co2, p = 35.1632 mb, t = 223.28 k
!     upper - n2o, p = 8.716e-2 mb, t = 226.03 k
! band 9:  1180-1390 cm-1 (low key - h2o,ch4; low minor - n2o)
!                         (high key - ch4; high minor - n2o)
!    old:  1180-1390 cm-1 (low - h2o,ch4; high - ch4)
! Minor gas mapping level :
!     lower - n2o, p = 706.272 mbar, t = 278.94 k
!     upper - n2o, p = 95.58 mbar, t = 215.7 k
! band 10:  1390-1480 cm-1 (low key - h2o; high key - h2o)
!     old:  1390-1480 cm-1 (low - h2o; high - h2o)
! band 11:  1480-1800 cm-1 (low - h2o; low minor - o2)
!                          (high key - h2o; high minor - o2)
!     old:  1480-1800 cm-1 (low - h2o; low minor - o2)
!                          (high key - h2o; high minor - o2)
! Minor gas mapping level :
!     lower - o2, p = 706.2720 mbar, t = 278.94 k
!     upper - o2, p = 4.758820 mbarm t = 250.85 k
! band 12:  1800-2080 cm-1 (low - h2o,co2; high - nothing)
!     old:  1800-2080 cm-1 (low - h2o,co2; high - nothing)
! band 13:  2080-2250 cm-1 (low key - h2o,n2o; high minor - o3 minor)
!     old:  2080-2250 cm-1 (low - h2o,n2o; high - nothing)
! NOTE: comment did not reflect code, see minor species for lower atm
! Minor gas mapping levels :
!     lower - co2, p = 1053.63 mb, t = 294.2 k
!     lower - co, p = 706 mb, t = 278.94 k
!     upper - o3, p = 95.5835 mb, t = 215.7 k
! band 14:  2250-2380 cm-1 (low - co2; high - co2)
!     old:  2250-2380 cm-1 (low - co2; high - co2)
! band 15:  2380-2600 cm-1 (low - n2o,co2; low minor - n2)
!                          (high - nothing)
!     old:  2380-2600 cm-1 (low - n2o,co2; high - nothing)
! Minor gas mapping level : 
!     Lower - Nitrogen Continuum, P = 1053., T = 294.
! band 16:  2600-3250 cm-1 (low key- h2o,ch4; high key - ch4)
!     old:  2600-3000 cm-1 (low - h2o,ch4; high - nothing)


MODULE mo_psrad_lrtm_kgs

  USE mo_psrad_general

  IMPLICIT NONE

  PUBLIC

  REAL(wp), PARAMETER :: &
    delwave(nbndlw) = (/& ! Spectral band width in wavenumbers
      340._wp, 150._wp, 130._wp,  70._wp, 120._wp, 160._wp, &
      100._wp, 100._wp, 210._wp,  90._wp, 320._wp, 280._wp, &
      170._wp, 130._wp, 220._wp, 650._wp/)

  REAL(wp), PARAMETER :: wavenum1(nbndlw) = (/ & !< Spectral band lower boundary in wavenumbers
    10._wp, 350._wp, 500._wp, 630._wp, 700._wp, 820._wp, &
    980._wp,1080._wp,1180._wp,1390._wp,1480._wp,1800._wp, &
    2080._wp,2250._wp,2380._wp,2600._wp/)
  REAL(wp), PARAMETER :: wavenum2(nbndlw) = (/ & !< Spectral band upper boundary in wavenumbers
    350._wp, 500._wp, 630._wp, 700._wp, 820._wp, 980._wp, &
    1080._wp,1180._wp,1390._wp,1480._wp,1800._wp,2080._wp, &
    2250._wp,2380._wp,2600._wp,3250._wp/)

  REAL(wp) :: chi_mls(ngas,59), planck_ratio(ngas,ngas,59), &
    ! planck function for each band, scaled by delwave
    totplanck(181,nbndlw), & 
    totplanck16(181) ! for band 16 NOTE: not used

  INTEGER, PARAMETER :: &
    ngpt(nbndlw) = (/10,12,16,14,16, 8,12,8,12,6, 8,8,4,2,2, 2/),&
    nsp(2,nbndlw) = RESHAPE((/& ! Number of reference lower/upper atmospheres
      1,1,9,9,9, 1,9,1,9,1, 1,9,9,1,9, 9, &
      1,1,5,5,5, 0,1,1,1,1, 1,0,0,1,0, 1/), &
      SHAPE=(/2,nbndlw/), ORDER=(/2,1/)), &
    nsp_fraction(2,nbndlw) = RESHAPE((/& ! Number of reference lower/upper atmospheres
      1,1,9,9,9, 1,9,1,9,1, 1,9,9,1,9, 9, &
      1,1,5,5,5, 1,1,1,1,1, 1,0,1,1,0, 1/), &
      SHAPE=(/2,nbndlw/), ORDER=(/2,1/)), &
    fracs_mult(2,nbndlw) = MAX(0,nsp_fraction-1), &
    !NOTE: Arrays differ only by nsp_*(2,13) and (2,6)
    nh2oref(2,nbndlw) = RESHAPE((/& 
      10,10,10,10,10, 10,10,10,10,10, 10,10,10,10,10, 10, &
      4,4,4,4,4,      4,4,4,4,4,     4,4,4,4,4, 4/), &
      SHAPE=(/2,nbndlw/), ORDER=(/2,1/)), &
    skip_atmosphere(2,nbndlw) = RESHAPE((/ &
      0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0, &
      0,0,0,0,0, 0,0,0,0,0, 0,1,0,0,1, 0/), &
      SHAPE=(/2,nbndlw/), ORDER=(/2,1/)), &
    copy_planck_fraction_from_other_band(nbndlw) = &
      (/0,0,0,0,0, 2,0,0,0,0, 0,0,0,0,0, 0/), &
    h2o_absorption_flag(2,2,nbndlw) = RESHAPE((/ &
      1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, & ! self,lower
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, & ! self,upper
      1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, & ! foreign, lower
      1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0/), & !foreign, upper
      SHAPE = (/2,2,nbndlw/), ORDER = (/3,2,1/)), &
    pressure_dependent_tau_correction(2,nbndlw) = RESHAPE((/ &
      1, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/), &
      SHAPE=(/2,nbndlw/), ORDER=(/2,1/)), &
    major_species(3,2,nbndlw) = RESHAPE((/ &
      ih2o,0,0, ih2o,0,0, & ! Band 1
      ih2o,0,0, ih2o,0,0, & ! Band 2
      ih2o,ico2,ih2oco2, ih2o,ico2,ih2oco2, & ! Band 3
      ih2o,ico2,ih2oco2, io3,ico2,io3co2, & ! Band 4
      ih2o,ico2,ih2oco2, io3,ico2,io3co2, & ! Band 5
      ih2o,0,0, 0,0,0, & ! Band 6
      ih2o,io3,ih2oo3, io3,0,0, & ! Band 7
      ih2o,0,0, io3,0,0, & ! Band 8
      ih2o,ich4,ih2och4,  ich4,0,0, & ! Band 9
      ih2o,0,0, ih2o,0,0, & ! Band 10
      ih2o,0,0, ih2o,0,0, & ! Band 11
      ih2o,ico2,ih2oco2, 0,0,0, & ! Band 12
      ih2o,in2o,ih2on2o, 0,0,0, & ! Band 13
      ico2,0,0, ico2,0,0, & ! Band 14
      in2o,ico2,in2oco2, 0,0,0, & ! Band 15
      ih2o,ich4,ih2och4, ich4,0,0/), & ! Band 16
      SHAPE=(/3,2,nbndlw/)), &
    planck_fraction_interpolation_layer(2,nbndlw) = RESHAPE((/ &
      0,0,9,11,5,   0,3,0,9,0, 0,10,5,0,1, 6, &
      0,0,13,13,43, 0,0,0,0,0, 0,0,0,0,0,  0/), &
      SHAPE = (/2,nbndlw/), ORDER=(/2,1/)), &
    stratosphere_fudge_idx(nbndlw) = &
      (/ 0,0,0,1,0, 0,2,0,0,0, 0,0,0,0,0, 0/), &
    use_fixed_secdiff(nbndlw) = &
      (/1,0,0,1,0, 0,0,0,0,1, 1,1,1,1,1, 1/)

  REAL(wp), PARAMETER :: &
    minor_species_fudge(4,6) = RESHAPE((/ &
      1.5_wp, 0.5_wp, 0.65_wp, 0.0_wp, &
      3.0_wp, 2.0_wp, 0.77_wp, 0.0_wp, &
      3.0_wp, 3.0_wp, 0.79_wp, 0.0_wp, &
      3.0_wp, 2.0_wp, 0.79_wp, 0.0_wp, &
      3.0_wp, 2.0_wp, 0.65_wp, 0.0_wp, &
      3.0_wp, 2.0_wp, 0.68_wp, 3.55e-4_wp/), SHAPE=(/4,6/)), &
    stratosphere_fudge(ngpt_orig,2) = RESHAPE((/ &
      1.,1.,1.,1.,1.,  1.,1.,0.92,0.88,1.07,     1.1,0.99,0.88,0.943,0., 0., &
      1.,1.,1.,1.,1.,  0.92,0.88,1.07,1.1,0.99,  0.855,1.,0.,0.,0., 0./), &
      SHAPE=(/ngpt_orig,2/))

  TYPE :: T !TMinorInfo => T for brevity
    INTEGER :: gas, M, N, scaling_type, fudge_index, interpolation_layer
  END TYPE

  TYPE(T), PARAMETER :: O = T(0,0,0,0,0,0), &
    minor_species(max_minor_species,2,nbndlw) = RESHAPE((/ &
     T(in2,19,0,1,0,0),O,O,O,O, & ! Band1, lower
     T(in2,19,0,1,0,0),O,O,O,O, & ! Band1, upper

     O,O,O,O,O, & !Band 2, lower
     O,O,O,O,O, & !Band 2, upper, and so on

     T(in2o,9,19,4,1,3),O,O,O,O, &
     T(in2o,5,19,4,1,13),O,O,O,O, &

     O,O,O,O,O, &
     O,O,O,O,O, &

     T(io3,9,19,0,0,7), T(iccl4, 1, 0,0,0,0),O,O,O, & !Band 5, lower
     O,O,O,O,O, & !Band 5, upper

     T(ico2,19,0,4,2,0), T(icfc11,1,0,0,0,0), T(icfc12,1,0,0,0,0), &
     O,O, & ! Band 6, lower
     T(icfc11,1,0,0,0,0), T(icfc12,1,0,0,0,0), O,O,O, & !Band 6, upper

     T(ico2,9,19,4,3,3),O,O,O,O, & !CHECK
     T(ico2,19,0,4,4,0),O,O,O,O, &

     T(ico2,19,0,4,5,0), T(io3,19,0,0,0,0), T(in2o,19,0,0,0,0), &
     T(icfc12,1,0,0,0,0), T(icfc22,1,0,0,0,0), & !Band 8, lower
     T(ico2,19,0,4,5,0), T(in2o,19,0,0,0,0), T(icfc12,1,0,0,0,0), &
     T(icfc22,1,0,0,0,0), O, & !Band8, upper

     T(in2o,9,19,4,1,3),O,O,O,O, &
     T(in2o,19,0,4,1,0),O,O,O,O, &

     O,O,O,O,O, &
     O,O,O,O,O, &

     T(io2,19,0,2,0,0),O,O,O,O, & !Band 11
     T(io2,19,0,2,0,0),O,O,O,O, &

     O,O,O,O,O, &
     O,O,O,O,O, &

     T(ico2,9,19,4,6,1),T(ico,9,19,0,0,3), O,O,O, &
     T(io3,19,0,0,0,0),O,O,O,O, &

     O,O,O,O,O, & !Band 14
     O,O,O,O,O, &

     T(in2,9,19,3,0,1),O,O,O,O, & 
     O,O,O,O,O, &

     O,O,O,O,O, &
     O,O,O,O,O/), &
     SHAPE=(/max_minor_species,2,nbndlw/))

  REAL(wp) :: precipitable_vapor_factor
  REAL(wp) :: pa1, pa2, pa3, pb1, pb2, pc1, pc2, pc3

  PUBLIC :: lrtm_count_flat
#ifdef PSRAD_WITH_LEGACY
  PUBLIC :: lrtm_pack_data
#else
  PUBLIC :: lrtm_unpack_data
#endif

CONTAINS

  SUBROUTINE lrtm_count_flat(n_lrtm)
    INTEGER, INTENT(OUT) :: n_lrtm
    n_lrtm = ngas * 59 + & ! chi_mls
      181 * nbndlw + & ! totplanck
      ngas * ngas * 59 + & ! planck_ratio
      9 !  precipitable_vapor_factor, pa1, pa2, pa3, pb1, pb2, pc1, pc2, pc3
  END SUBROUTINE lrtm_count_flat

#ifdef PSRAD_WITH_LEGACY

  SUBROUTINE lrtm_pack_data(lrtm_data)
    REAL(wp), INTENT(OUT) :: lrtm_data(:)
    INTEGER :: offset, s, e, l
    
    offset = 0

    l = ngas * 59
    s = offset + 1
    e = offset + l
    lrtm_data(s:e) = RESHAPE(chi_mls, SHAPE = (/l/))
    offset = offset + l

    l = 181 * nbndlw
    s = offset + 1
    e = offset + l
    lrtm_data(s:e) = RESHAPE(totplanck, SHAPE = (/l/))
    offset = offset + l

    l = ngas * ngas * 59
    s = offset + 1
    e = offset + l
    lrtm_data(s:e) = RESHAPE(planck_ratio, SHAPE = (/l/))
    offset = offset + l

    l = 9
    s = offset + 1
    e = offset + l
    lrtm_data(s:e) = (/precipitable_vapor_factor, &
      pa1, pa2, pa3, pb1, pb2, pc1, pc2, pc3/)
  
  END SUBROUTINE lrtm_pack_data
    
#else

  SUBROUTINE lrtm_unpack_data(lrtm_data)
    REAL(wp), INTENT(IN) :: lrtm_data(:)
    INTEGER :: offset, s, e, l
    
    offset = 0

    l = ngas * 59
    s = offset + 1
    e = offset + l
    chi_mls = RESHAPE(lrtm_data(s:e), SHAPE = (/ngas,59/))
    offset = offset + l

    l = 181 * nbndlw
    s = offset + 1
    e = offset + l
    totplanck = RESHAPE(lrtm_data(s:e), SHAPE = (/181,nbndlw/))
    offset = offset + l

    l = ngas * ngas * 59
    s = offset + 1
    e = offset + l
    planck_ratio = RESHAPE(lrtm_data(s:e), SHAPE = (/ngas,ngas,59/))
    offset = offset + l

    precipitable_vapor_factor = lrtm_data(offset + 1)
    pa1 = lrtm_data(offset + 2)
    pa2 = lrtm_data(offset + 3)
    pa3 = lrtm_data(offset + 4)
    pb1 = lrtm_data(offset + 5)
    pb2 = lrtm_data(offset + 6)
    pc1 = lrtm_data(offset + 7)
    pc2 = lrtm_data(offset + 8)
    pc3 = lrtm_data(offset + 9)
  
  END SUBROUTINE lrtm_unpack_data

#endif

END MODULE mo_psrad_lrtm_kgs
