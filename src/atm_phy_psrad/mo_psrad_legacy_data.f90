MODULE mo_psrad_legacy_data

#ifdef PSRAD_WITH_LEGACY
  USE mo_psrad_general, ONLY: wp, nbndlw, nbndsw, max_minor_species

  PUBLIC

  TYPE ptr1
    REAL(wp), ALLOCATABLE :: v(:)
  END TYPE ptr1
  TYPE ptr2
    REAL(wp), ALLOCATABLE :: v(:,:)
  END TYPE ptr2
  TYPE ptr3
    REAL(wp), ALLOCATABLE :: v(:,:,:)
  END TYPE ptr3
  TYPE ptr4
    REAL(wp), ALLOCATABLE :: v(:,:,:,:)
  END TYPE ptr4


  TYPE(ptr2), DIMENSION(2,nbndlw) :: lw_kmajor_red
  TYPE(ptr2), DIMENSION(2,nbndlw) :: lw_h2oref_red
  TYPE(ptr2), DIMENSION(max_minor_species,2,nbndlw) :: lw_kgas_red
  TYPE(ptr2), DIMENSION(2,nbndlw) :: lw_planck_red

  TYPE(ptr2), DIMENSION(2,nbndsw) :: sw_kmajor_red
  TYPE(ptr2), DIMENSION(2,nbndsw) :: sw_h2oref_red
  TYPE(ptr1), DIMENSION(2,nbndsw) :: sw_kgas_red
  REAL(wp) :: sw_rayl0(nbndsw)
  TYPE(ptr1) :: sw_rayl1_red(nbndsw)
  TYPE(ptr2) :: sw_rayl2_red(2,nbndsw)
  TYPE(ptr2), DIMENSION(nbndsw) :: sw_sfluxref_red
#endif

END MODULE mo_psrad_legacy_data

