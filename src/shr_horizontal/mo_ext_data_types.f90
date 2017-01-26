!>
!! External data type definition.
!!
!! External data type definition.
!!
!! @author Daniel Reinert, DWD
!! @author Hermann Asensio, DWD
!!
!!
!! @par Revision History
!! Initial revision by Daniel Reinert, DWD (2012-03-22)
!! Separated from former module mo_ext_data
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_ext_data_types

  USE mo_kind,               ONLY: wp
  USE mo_fortran_tools,      ONLY: t_ptr_2d3d, t_ptr_i2d3d 
  USE mo_linked_list,        ONLY: t_var_list

  IMPLICIT NONE


  PRIVATE

  PUBLIC :: t_external_data
  PUBLIC :: t_external_atmos
  PUBLIC :: t_external_atmos_td
  PUBLIC :: t_external_ocean
  PUBLIC :: t_external_bgc



  !>
  !! atmosphere external data class
  !!
  !! atmosphere external data class
  !!
  TYPE :: t_external_atmos

    !
    ! *** Topography ***
    REAL(wp), POINTER ::   &   !< topographic height at cell centers      [m]
      &  topography_c(:,:)     ! index1=1,nproma, index2=1,nblks_c

    REAL(wp), POINTER ::   &   !< geopotential (S)                        [m**2/s**2]
      &  fis(:,:)              ! index1=1,nproma, index2=1,nblks_c

    REAL(wp), POINTER ::   &   !< gradient of topographic height at cell centers [ ]
      &  grad_topo(:,:,:)      ! index1=1,2, index2=1,nproma, index3=1,nblks_c


    !
    ! *** Land-Sea-Mask ***

    LOGICAL, POINTER  ::   &   !< land-sea-mask for cell centers          [ ]
      &  llsm_atm_c(:,:)       ! .TRUE. if landpoint
                               ! index1=1,nproma, index2=1,nblks_c

    LOGICAL, POINTER  ::   &   !< mask function for lake points           [ ]
      &  llake_c(:,:)          ! .TRUE. if lake point
                               ! index1=1,nproma, index2=1,nblks_c

    INTEGER, POINTER  ::   &   !< land-sea-mask for cell centers          [ ]
      &  lsm_ctr_c(:,:)        !  index1=1,nproma, index2=1,nblks_c

    REAL(wp), POINTER  ::   &  !< elevation at cell centers               [m]
      &  elevation_c(:,:)

    REAL(wp), POINTER ::   &   !< fraction land in a grid element         [ ]
      &  fr_land(:,:)          ! 0. for water, 1.0 indicates 100% land
                               ! index1=1,nproma, index2=1,nblks_c

    REAL(wp), POINTER ::    &  !< fraction land glacier in a grid element [ ]
      &  fr_glac(:,:)          ! 1.0 indicates 100% glacier
                               ! index1=1,nproma, index2=1,nblks_c    

    REAL(wp), POINTER ::   &   !< fraction land in a grid element         [ ]
      &  fr_land_smt(:,:)      !  = smoothed fr_land

    REAL(wp), POINTER ::   &   !< fraction land glacier in a grid element [ ]
      &  fr_glac_smt(:,:)      ! = smoothed fr_glac


    !  
    ! *** roughness length ***
    REAL(wp), POINTER ::   &   !< surface roughness                       [m]
      &  z0(:,:)               ! index1=1,nproma, index2=1,nblks_c


    !
    ! *** FLake ***
    REAL(wp), POINTER ::   &   !< fraction of fresh water                 [ ]
      &  fr_lake(:,:)          ! as partition of total area of the
                               ! grid element
                               ! index1=1,nproma, index2=1,nblks_c

    REAL(wp), POINTER ::   &   !< lake depth                              [m]
      &  depth_lk(:,:)         ! index1=1,nproma, index2=1,nblks_c

    REAL(wp), POINTER ::   &   !< wind fetch over lake                    [m]
      &  fetch_lk(:,:)         ! index1=1,nproma, index2=1,nblks_c

    REAL(wp), POINTER ::   &   !< depth of the thermally active layer of  [m]
      &  dp_bs_lk(:,:)         !< bottom sediments
                               ! index1=1,nproma, index2=1,nblks_c

    REAL(wp), POINTER ::   &   !< climatological temperature at the       [K]
      &  t_bs_lk(:,:)          !< bottom of the thermally active layer of sediments
                               ! index1=1,nproma, index2=1,nblks_c

    REAL(wp), POINTER ::   &   !< attenuation coefficient of the lake     [m^-1]
      &  gamso_lk(:,:)         !< water with respect to solar radiation
                               ! index1=1,nproma, index2=1,nblks_c



    !
    ! *** subgrid scale orography ***
    REAL(wp), POINTER ::   &   !< standard deviation of sub-grid scale orography [m]
      &  sso_stdh(:,:)         ! index1=1,nproma, index2=1,nblks_c

    REAL(wp), POINTER ::   &   !< "raw" sso_stdh without correction for orography smoothing [m]
      &  sso_stdh_raw(:,:)     ! index1=1,nproma, index2=1,nblks_c

    REAL(wp), POINTER ::   &   !< anisotropy of sub-grid scale orography  [ ]
      &  sso_gamma(:,:)        ! index1=1,nproma, index2=1,nblks_c

    REAL(wp), POINTER ::   &   !< angle betw. principal axis of orography and E [rad]
      &  sso_theta(:,:)        ! index1=1,nproma, index2=1,nblks_c

    REAL(wp), POINTER ::   &   !< mean slope of sub-grid scale orography  [ ]
      &  sso_sigma(:,:)        ! index1=1,nproma, index2=1,nblks_c


    !
    ! *** vegetation parameters ***
    REAL(wp), POINTER ::   &   !< ground fraction covered by plants (vegetation period)  [ ]
      & plcov_mx(:,:)          ! index1=1,nproma, index2=1,nblks_c
    REAL(wp), POINTER ::   &   !< ground fraction covered by plants (vegetation period)  [ ]
      & plcov(:,:)             ! index1=1,nproma, index2=1,nblks_c
    REAL(wp), POINTER ::   &   !< ground fraction covered by plants (vegetation period)  [ ]
      & plcov_t(:,:,:)         ! index1=1,nproma, index2=1,nblks_c, ntiles_total

    REAL(wp), POINTER ::   &   !< max leaf area index (vegetation period)        [ ]
      &  lai_mx(:,:)           ! index1=1,nproma, index2=1,nblks_c
    REAL(wp), POINTER ::   &   !< leaf area index (vegetation period)            [ ]
      &  lai(:,:)              ! index1=1,nproma, index2=1,nblks_c
    REAL(wp), POINTER ::   &   !< surface area index (vegetation period)         [ ]
      &  sai(:,:)              ! index1=1,nproma, index2=1,nblks_c
    REAL(wp), POINTER ::   &   !< transpiration area index (vegetation period)   [ ]
      &  tai(:,:)              ! index1=1,nproma, index2=1,nblks_c
    REAL(wp), POINTER ::   &   !< earth area (evaporative surface area)          [ ] 
      &  eai(:,:)              ! index1=1,nproma, index2=1,nblks_c
    REAL(wp), POINTER ::   &   !< surface area index (vegetation period)         [ ]
      &  sai_t(:,:,:)          ! index1=1,nproma, index2=1,nblks_c, ntiles_total
    REAL(wp), POINTER ::   &   !< transpiration area index (vegetation period)   [ ]
      &  tai_t(:,:,:)          ! index1=1,nproma, index2=1,nblks_c, ntiles_total
    REAL(wp), POINTER ::   &   !< earth area (evaporative surface area)          [ ]
      &  eai_t(:,:,:)          !< index (vegetation period)
                               ! index1=1,nproma, index2=1,nblks_c, ntiles_total

    REAL(wp), POINTER ::   &   !< root depth                              [m]
      &  rootdp(:,:)           ! index1=1,nproma, index2=1,nblks_c
    REAL(wp), POINTER ::   &   !< root depth                              [m]
      &  rootdp_t(:,:,:)       ! index1=1,nproma, index2=1,nblks_c, ntiles_total

    REAL(wp), POINTER ::   &   !< ground fraction covered by evergreen forest [ ]
      &  for_e(:,:)            ! index1=1,nproma, index2=1,nblks_c

    REAL(wp), POINTER ::   &   !< ground fraction covered by deciduous forest [ ]
      &  for_d(:,:)            ! index1=1,nproma, index2=1,nblks_c

    REAL(wp), POINTER ::   &   !< minimum value of stomata resistance     [ s/m ]
      &  rsmin(:,:)            ! index1=1,nproma, index2=1,nblks_c
    REAL(wp), POINTER ::   &   !< minimum value of stomata resistance     [ s/m ]
      &  rsmin2d_t(:,:,:)      ! index1=1,nproma, index2=1,nblks_c, ntiles_total

    REAL(wp), POINTER ::   &   !< annual maximum NDVI                     [ ]
      &  ndvi_max(:,:)         ! index1=1,nproma, index2=1,nblks_c

    REAL(wp), POINTER ::   &   !< proportion of actual value/maximum 
      &  ndviratio(:,:)        !< normalized differential vegetation index [ ]
                               !< for starting time of model integration
                               !< (derived from atm_td%ndvi_mrat)
                               ! index1=1,nproma, index2=1,nblks_c


    !
    ! *** soil parameters ***
    INTEGER, POINTER  ::   &   !< soil texture, keys 0-9                  [ ]
      &  soiltyp(:,:)          ! index1=1,nproma, index2=1,nblks_c
    ! soiltyp_t refers to the land point index list
    ! this field is dimensioned with ntiles_total even though this appears to be
    ! unnecessary in order not to disturb the runtime optimization
    INTEGER, POINTER  ::   &   !< soil texture, keys 0-9                  [ ]
      &  soiltyp_t(:,:,:)      ! index1=1,nproma, index2=1,nblks_c, index3=ntiles_total


    REAL(wp), POINTER ::   &   !< Near surface temperature (climatological mean)  [ K ]
      &  t_cl(:,:)             !  used as climatological layer (deepest layer) of T_SO
                               ! index1=1,nproma, index2=1,nblks_c


    ! *** radiation parameters ***
    REAL(wp), POINTER ::   &   !< ozone mixing ratio                        [ kg kg^-1 ]
      &  o3(:,:,:)             ! index1=1,nproma, index2=nlev,index3=1,nblks_c

    REAL(wp), POINTER ::   &   !< longwave surface emissivity             [ ]
      &  emis_rad(:,:)         ! index1=1,nproma, index2=1,nblks_c

    REAL(wp), POINTER ::   &   !< shortwave broadband albedo for diffuse radiation  [1]
      &  alb_dif(:,:)          !< (0.3 - 5.0 um)
                               ! index1=1,nproma, index2=1,nblks_c

    REAL(wp), POINTER ::   &   !< UV visible broadband albedo for diffuse radiation [1]
      &  albuv_dif(:,:)        !< (0.3 - 0.7 um)
                               ! index1=1,nproma, index2=1,nblks_c

    REAL(wp), POINTER ::   &   !< Near IR broadband albedo for diffuse radiation    [1]
      &  albni_dif(:,:)        !< (0.7 - 5.0 um)
                               ! index1=1,nproma, index2=1,nblks_c



    ! *** flow control parameters for tile approach ***
    REAL(wp), POINTER ::  &    !< Landuse class fraction                  [ ]
      & lu_class_fraction(:,:,:) ! index1=1,nproma, index2=1,nblks_c, index3=1,nclass_lu

    INTEGER, POINTER ::  &    !< Static land point index list for each block  [ ]
      & idx_lst_lp(:,:)       ! index1=1,nproma, index2=1,nblks_c
    INTEGER, POINTER ::  &    !< Land point count per block       [ ]
      & lp_count(:)           ! index1=1,nblks_c
    INTEGER, POINTER ::  &    !< Static sea point index list for each block   [ ]
      & idx_lst_sp(:,:)       ! index1=1,nproma, index2=1,nblks_c
    INTEGER, POINTER ::  &    !< Sea point count per block        [ ]
      & sp_count(:)           ! index1=1,nblks_c
    INTEGER, POINTER ::  &    !< static lake point index list for each block  [ ]
      & idx_lst_fp(:,:)       ! index1=1,nproma, index2=1,nblks_c
    INTEGER, POINTER ::  &    !< Lake point count per block        [ ]
      & fp_count(:)           ! index1=1,nblks_c

    INTEGER, POINTER ::  &    !< Static grid point index list for each block and tile [ ]
      & idx_lst_lp_t(:,:,:)   ! index1=1,nproma, index2=1,nblks_c, index3=ntiles_total
    INTEGER, POINTER ::  &    !< Corresponding grid point count per block and tile index      [ ]
      & lp_count_t(:,:)       ! index1=1,nblks_c, index2=ntiles_total

    INTEGER, POINTER ::  &    !< Land cover class for each tile index  [ ]
      & lc_class_t(:,:,:)     ! index1=1,nproma, index2=1,nblks_c, index3=ntiles_total
    REAL(wp), POINTER ::  &   !< Normalized land cover fraction for each tile index  [ ]
      & lc_frac_t(:,:,:)      ! index1=1,nproma, index2=1,nblks_c, index3=ntiles_total

    INTEGER, POINTER ::  &    !< Dynamic grid point index list (if lsnowtile=true) for each block and tile [ ]
      & idx_lst_t(:,:,:)      ! index1=1,nproma, index2=1,nblks_c, index3=ntiles_total
    INTEGER, POINTER ::  &    !< Corresponding grid point count per block and tile index      [ ]
      & gp_count_t(:,:)       ! index1=1,nblks_c, index2=ntiles_total

    INTEGER, POINTER ::  &    !< Snowtile flag field [ ]
      & snowtile_flag_t(:,:,:)! index1=1,nproma, index2=1,nblks_c, index3=ntiles_total
                              ! -1: no separation between snow tile and snow-free tile
                              !  0: inactive
                              !  1: active
                              !  2: newly activated; initialization from corresponding tile required
    REAL(wp), POINTER ::  &   !< Actual area fraction for each tile index  [ ]
      & frac_t(:,:,:)         ! index1=1,nproma, index2=1,nblks_c, index3=ntiles_total

    REAL(wp), POINTER ::  &   !< Inverse of fr_land derived from actual land tile fractions [ ]
      & inv_frland_from_tiles(:,:) ! needed for aggregation of land-only fields. 
                              ! Approximately equal to inverse of  
                              ! fr_land (extpar) + fr_lake(extpar, where fr_lake<frlake_thrhld)
                              ! index1=1,nproma, index2=1,nblks_c

    ! Sub-lists for sea points (idx_lst_sp), in order to distinguish between ice-covered and open
    ! sea points.
    ! 
    INTEGER, POINTER ::  &    !< Dynamic sea water point index list for each block and tile [ ]
      & idx_lst_spw(:,:)      ! index1=1,nproma, index2=1,nblks_c
    INTEGER, POINTER ::  &    !< Corresponding grid point count per block                   [ ]
      & spw_count(:)          ! index1=1,nblks_c
    INTEGER, POINTER ::  &    !< Dynamic sea ice point index list for each block and tile   [ ]
      & idx_lst_spi(:,:)      ! index1=1,nproma, index2=1,nblks_c
    INTEGER, POINTER ::  &    !< Corresponding grid point count per block                   [ ]
      & spi_count(:)          ! index1=1,nblks_c

 

    ! *** storage for lookup table data for each landuse class ***
    ! (needed to simplify switching between GLC2000 and Globcover2009, which
    !  have a different order of the landcover classes)
    REAL(wp), POINTER ::  &    !< Land-cover related roughness length   [m]
      & z0_lcc(:)              ! index1=1,23
    REAL(wp), POINTER ::  &    !< Minimum land-cover related roughness length   [m]
      & z0_lcc_min(:)          ! index1=1,23
    REAL(wp), POINTER ::  &    !< Maximum plant cover fraction for each land-cover class  [ ]
      & plcovmax_lcc(:)        ! index1=1,23
    REAL(wp), POINTER ::  &    !< Maximum leaf area index for each land-cover class  [ ]
      & laimax_lcc(:)          ! index1=1,23
    REAL(wp), POINTER ::  &    !< Maximum root depth for each land-cover class  [ ]
      & rootdmax_lcc(:)        ! index1=1,23
    REAL(wp), POINTER ::  &    !< Minimum stomata resistance for each land-cover class  [ ]
      & stomresmin_lcc(:)      ! index1=1,23
    REAL(wp), POINTER ::  &    !< Albedo in case of snow cover for each land-cover class  [ ]
      & snowalb_lcc(:)         ! index1=1,23
    LOGICAL, POINTER ::  &     !< Existence of separate snow tiles for land-cover class [ ]
      & snowtile_lcc(:)        ! index1=1,23
    INTEGER :: i_lc_snow_ice   !< Land-cover classification index for snow and ice
    INTEGER :: i_lc_water      !< Land-cover classification index for water
    INTEGER :: i_lc_urban      !< Land-cover classification index for urban / artificial surface
    INTEGER :: i_lc_shrub_eg   !< Shrub cover Grassland/Forest // Evergreen. Currently used by ART only
    INTEGER :: i_lc_shrub      !< Closed to open Shrubland (deciduous). Currently used by ART only
    INTEGER :: i_lc_grass      !< Grassland//herbaceous. Currently used by ART only
    INTEGER :: i_lc_bare_soil  !< Land-cover classification index for bare soil. Currently used by ART only
    INTEGER :: i_lc_sparse     !< Land-cover classification index for sparse vergetation. Currently used by ART only

    ! for output purposes.
    TYPE(t_ptr_i2d3d), ALLOCATABLE :: lc_class_t_ptr(:)
    TYPE(t_ptr_2d3d),  ALLOCATABLE :: frac_t_ptr(:)
    TYPE(t_ptr_2d3d),  ALLOCATABLE :: plcov_t_ptr(:)

  END TYPE t_external_atmos



  !>
  !! atmosphere external data class (time dependent)
  !!
  !! Contains auxiliary time dependent versions of some external atmospheric data
  !! fields which are already defined in external_atmos. These fields will
  !! be used to store e.g. montly means from which updated external data can be
  !! derived. The updated interpolated fields can be copied into the time independent
  !! counterparts, which are defined in external_atmos.
  !!
  TYPE :: t_external_atmos_td

    ! *** radiation parameters ***
    REAL(wp), POINTER ::   &   !< ozone mixing ratio    [ ]
      &  o3(:,:,:,:)           ! index1=1,nproma, index2=nlev_o3,
                               ! index3=1,nblks_c, index4=1,ntimes

    REAL(wp),POINTER::  &
      &   zf  (:),      &      !full levels of ozone gemometric height
      &   pfoz(:),      &      !full levels of of ozone pressure
      &   phoz(:)              !half levels of ozone pressure field 
    !
    ! *** radiation parameters ***
    REAL(wp), POINTER ::   &   !< aerosol optical thickness of black carbon    [ ]
      &  aer_bc(:,:,:)         ! index1=1,nproma, index2=1,nblks_c, index3=1,ntimes

    REAL(wp), POINTER ::   &   !< aerosol optical thickness of ambient aerosol [ ]
      &  aer_dust(:,:,:)       ! index1=1,nproma, index2=1,nblks_c, index3=1,ntimes

    REAL(wp), POINTER ::   &   !< aerosol optical thickness of particulate     [ ]
      &  aer_org(:,:,:)        !< organic_matter_ambient_aerosol             
                               ! index1=1,nproma, index2=1,nblks_c, index3=1,ntimes

    REAL(wp), POINTER ::   &   !< aerosol optical thickness of sulfate aerosol [ ]
      &  aer_so4(:,:,:)        ! index1=1,nproma, index2=1,nblks_c, index3=1,ntimes

    REAL(wp), POINTER ::   &   !< aerosol optical thickness of seasalt aerosol [ ]
      &  aer_ss(:,:,:)         ! index1=1,nproma, index2=1,nblks_c, index3=1,ntimes

    REAL(wp), POINTER ::   &   !< shortwave broadband albedo for diffuse radiation  [1]
      &  alb_dif(:,:,:)        !< (0.3 - 5.0 um)
                               ! index1=1,nproma, index2=1,nblks_c, index3=1,ntimes

    REAL(wp), POINTER ::   &   !< UV visible broadband albedo for diffuse radiation [1]
      &  albuv_dif(:,:,:)      !< (0.3 - 0.7 um)
                               ! index1=1,nproma, index2=1,nblks_c, index3=1,ntimes

    REAL(wp), POINTER ::   &   !< Near IR broadband albedo for diffuse radiation    [1]
      &  albni_dif(:,:,:)      !< (0.7 - 5.0 um)
                               ! index1=1,nproma, index2=1,nblks_c, index3=1,ntimes


    !
    ! *** vegetation parameters ***
    REAL(wp), POINTER ::   &   !< (monthly) proportion of actual value/maximum 
      &  ndvi_mrat(:,:,:)      !< normalized differential vegetation index   [ ]
                               ! index1=1,nproma, index2=1,nblks_c, index3=1,ntimes
    !
    ! ***SST and sea ice fraction
    REAL(wp), POINTER ::   &   !< (monthly) SST
      &  sst_m(:,:,:)          ! index1=1,nproma, index2=1,nblks_c, index3=1,ntimes
    REAL(wp), POINTER ::   &   !< (monthly) sea ice fraction
      &  fr_ice_m(:,:,:)       ! index1=1,nproma, index2=1,nblks_c, index3=1,ntimes

  END TYPE t_external_atmos_td


  TYPE :: t_external_bgc
   REAL(wp), POINTER ::   &   !< (monthly) dust deposition
      &  dust(:,:,:)   
   REAL(wp), POINTER ::   &   !< (monthly) nitrogen deposition
      &  nitro(:,:,:)   
  END TYPE t_external_bgc

  !>
  !! ocean external data class
  !!
  !! ocean external data class
  !!
  TYPE :: t_external_ocean

    ! ocean topography <=> bathymetric height used in the ocean 
    ! cell centers and edges only
    !
    REAL(wp), POINTER ::   &  !<  bathymetric height at cell centers [m]
      &  bathymetry_c(:,:)    !  index1=1,nproma, index2=1,nblks_c

    REAL(wp), POINTER ::   &  !< topographic height at cell edges    [m]
      &  bathymetry_e(:,:)    !  index1=1,nproma, index2=1,nblks_e

  ! REAL(wp), POINTER ::   &  !< topographic height at cell vertices [m]
  !   &  bathymetry_v(:,:)    !  index1=1,nproma, index2=1,nblks_v

    ! *** Land-Sea-Mask ***
    INTEGER, POINTER  ::   &  !< land-sea-mask for cell centers          [ ]
      &  lsm_ctr_c(:,:)       !  index1=1,nproma, index2=1,nblks_c
    INTEGER, POINTER ::    &  !< land-sea-mask for cell edges
      &  lsm_ctr_e(:,:)       !  index1=1,nproma, index2=1,nblks_e

    ! OMIP/NCEP type flux-forcing fluxes on cell centers. no_of_fluxes=3
    !
    REAL(wp), POINTER ::   &       !< omip monthly/daily mean forcing fluxes
      &  flux_forc_mon_c(:,:,:,:)  !  index1=nproma, index2=time, index3=nblks_c, index4=no_of_fluxes

  END TYPE t_external_ocean



!  !>
!  !! ocean external data class (time dependent)
!  !!
!  !! This data type contains auxiliary time dependent versions of
!  !! some external oceanic data fields already defined in external_ocean. These
!  !! fields will be used to store e.g. montly means from which interpolated external
!  !! data can be derived. The updated fields are copied into the time independent
!  !! counterparts which are defined in external_ocean.
!  !!
!  TYPE :: external_ocean_td
!  END TYPE external_ocean_td



  !>
  !! External data class including lists
  !!
  !! External data class including lists
  !!
  TYPE :: t_external_data

    TYPE(t_external_atmos)    :: atm
    TYPE(t_var_list)          :: atm_list

    TYPE(t_external_atmos_td) :: atm_td
    TYPE(t_var_list)          :: atm_td_list

    TYPE(t_external_ocean)    :: oce
    TYPE(t_var_list)          :: oce_list

    TYPE(t_external_bgc)      :: bgc
    TYPE(t_var_list)          :: bgc_list
!    TYPE(t_external_ocean_td) :: oce_td
!    TYPE(t_var_list), POINTER :: oce_td_list

  END TYPE t_external_data



END MODULE mo_ext_data_types

