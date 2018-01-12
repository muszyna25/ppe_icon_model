!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!

MODULE mo_psrad_lrtm_kgs

  USE mo_psrad_general ,ONLY : wp, nbndlw, &
    maxperband, max_minor_species, &
    ngas, ncfc
  IMPLICIT NONE

  PUBLIC

  SAVE

  REAL(wp) :: chi_mls(ngas,59), planck_ratio(ngas,ngas,59), &
    totplanck(181,nbndlw), & ! planck function for each band
    totplanck16(181) ! for band 16

  INTEGER, PARAMETER :: no  = 16, & ! original abs coefficients, all bands
    ng(nbndlw) = (/10, 12, 16, 14, 16, 8, 12, 8, 12, 6, 8, 8, 4, 2, 2, 2/), &
    nsp_fraction(2,nbndlw) = RESHAPE((/& ! Number of reference lower/upper atmospheres
      1,1,9,9,9, 1,9,1,9,1, 1,9,9,1,9,9, &
      1,1,5,5,5, 1,1,1,1,1, 1,0,1,1,0,1/), &
      SHAPE=(/2,nbndlw/), ORDER=(/2,1/)), &
    nsp_species(2,nbndlw) = RESHAPE((/& ! Number of reference lower/upper atmospheres
      1,1,9,9,9, 1,9,1,9,1, 1,9,9,1,9,9, &
      1,1,5,5,5, 0,1,1,1,1, 1,0,0,1,0,1/), &
      SHAPE=(/2,nbndlw/), ORDER=(/2,1/)), &
    !NOTE: Arrays differ only by nsp_*(2,13) and (2,6)
    nsp_species_broken_16(2,nbndlw) = RESHAPE((/& 
    ! Number of reference lower/upper atmospheres for gas_optics, 
    ! differs by nsp_*(2,16)
      1,1,9,9,9, 1,9,1,9,1, 1,9,9,1,9,9, &
      1,1,5,5,5, 0,1,1,1,1, 1,0,0,1,0,0/), &
      SHAPE=(/2,nbndlw/), ORDER=(/2,1/)), &
    npressure(2) = (/13,47/), &
    nh2oref(2) = (/10, 4/)! #self, #foreign

  TYPE ptr1
    REAL(wp), POINTER :: v(:)
  END TYPE ptr1
  TYPE ptr2
    REAL(wp), POINTER :: v(:,:)
  END TYPE ptr2
  TYPE ptr3
    REAL(wp), POINTER :: v(:,:,:)
  END TYPE ptr3
  TYPE ptr4
    REAL(wp), POINTER :: v(:,:,:,:)
  END TYPE ptr4

  TYPE(ptr2), DIMENSION(2,nbndlw) :: kmajor
  TYPE(ptr3), DIMENSION(2,nbndlw) :: kmajor3_in
  TYPE(ptr4), DIMENSION(2,nbndlw) :: kmajor4_in
  TYPE(ptr2), DIMENSION(ngas,2,nbndlw) :: kgas2, kgas2_delta, kgas2_in
  TYPE(ptr3), DIMENSION(ngas,2,nbndlw) :: kgas3, &
    kgas3_delta, kgas3_in
  INTEGER :: kgas2_list(ngas,2,nbndlw), kgas3_list(2,ngas,2,nbndlw)

  TYPE(ptr1), DIMENSION(ncfc,nbndlw) :: cfc, cfc_in
  INTEGER :: cfc_list(ncfc,nbndlw) ! 1/2 - reduce with/without weighing

  REAL(wp), DIMENSION(MAXVAL(nh2oref),maxperband,2,nbndlw) :: &
    h2oref_in, h2oref, h2oref_delta

  REAL(wp), DIMENSION(maxperband,2,nbndlw) :: planck_fraction1_in, &
    planck_fraction1
  TYPE(ptr2), DIMENSION(2,nbndlw) :: planck_fraction2_in, &
    planck_fraction2, planck_fraction2_delta

  INTEGER :: planck_fraction_interpolation_layer(2,nbndlw)

  INTEGER :: skip_atmosphere(2,nbndlw)
  INTEGER :: missing_planck_fraction(nbndlw)

  INTEGER :: key_species(2,2,nbndlw), &
    minor_species(max_minor_species,2,nbndlw), &
    minor_species_scale(max_minor_species,2,nbndlw), &
    minor_species_interpolation_layer(max_minor_species,2,nbndlw)
  REAL(wp) :: minor_species_fudge(4,max_minor_species,2,nbndlw)
  INTEGER :: n_key_species(2,nbndlw), n_minor_species(2,nbndlw)
  INTEGER :: reaction_table(ngas,ngas)
  REAL(wp) :: stratosphere_fudge(maxperband,nbndlw)
  INTEGER :: stratosphere_fudge_flag(nbndlw)

  INTEGER :: h2o_absorption_flag(2,2,nbndlw)
  INTEGER :: pressure_dependent_tau_correction(2,nbndlw)


!! @brief LRTM Band 1 Parameters: 10-250 cm-1 (low - h2o; high - h2o)
!! @brief LRTM Band 2 Parameters 250-500 cm-1 (low - h2o; high - h2o)
!! @brief LRTM Band 3 Parameters 500-630 cm-1 (low - h2o,co2; high - h2o,co2)
!! @brief LRTM Band 4 Parameters 630-700 cm-1 (low - h2o,co2; high - o3,co2)
!! @brief LRTM Band 5 Parameters 700-820 cm-1 (low - h2o,co2; high - o3,co2)
!! @brief LRTM Band 6 Parameters 820-980 cm-1 (low - h2o; high - nothing)
!! @brief LRTM Band 7 Parameters 980-1080 cm-1 (low - h2o,o3; high - o3)
!! @brief LRTM Band 8 Parameters 1080-1180 cm-1 (low (i.e.>~300mb) - h2o; high - o3)
!! @brief LRTM Band 9 parameters 1180-1390 cm-1 (low - h2o,ch4; high - ch4)
!! @brief LRTM Band 10 Parameters 1390-1480 cm-1 (low - h2o; high - h2o)
!! @brief LRTM Band 11 Parameters 1480-1800 cm-1 (low - h2o; high - h2o)
!! @brief LRTM Band 12 Parameters 1800-2080 cm-1 (low - h2o,co2; high - nothing)
!! @brief LRTM Band 13 Parameters 2080-2250 cm-1 (low - h2o,n2o; high - nothing)
!! @brief LRTM Band 14 Parameters 2250-2380 cm-1 (low - co2; high - co2)
!! @brief LRTM Band 15 Parameters 2380-2600 cm-1 (low - n2o,co2; high - nothing)
!! @brief LRTM Band 16 Parameters 2600-3000 cm-1 (low - h2o,ch4; high - nothing)

END MODULE mo_psrad_lrtm_kgs
