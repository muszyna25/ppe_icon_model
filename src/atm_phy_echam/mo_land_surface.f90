MODULE mo_land_surface

!  USE mo_jsbach_grid,   ONLY: grid_type, domain_type, kstart, kend, nidx
!  USE mo_jsbach_lctlib, ONLY: lctlib_type
  USE mo_kind,          ONLY: dp=>wp

  IMPLICIT NONE

!  PUBLIC ::  init_land_surface
!  PUBLIC ::  init_albedo
  PUBLIC ::  update_albedo_echam5
!  PUBLIC ::  update_albedo
!  PUBLIC ::  update_albedo_temp
!  PUBLIC ::  update_albedo_snow_age
!  PUBLIC ::  land_surface_diagnostics
!  PUBLIC ::  scale_cover_fract
!  PUBLIC ::  define_masks

!!$  TYPE land_surface_type
!!$     INTEGER              :: ntiles
!!$     INTEGER, POINTER, DIMENSION(:,:)     :: & !! (nland, ntiles)
!!$          cover_type                           !! Index into LctLibrary
!!$     REAL(dp),    POINTER, DIMENSION(:,:) :: & !! (nland, ntiles)
!!$          cover_type_real,                   & !! cover_type converted to REAL for netCDF files
!!$          cover_fract,                       & !! Fraction of coverage for each land cover type
!!$          cover_fract_pot,                   & !! Fraction of coverage for potential natural vegetation 
!!$          albedo,                            & !! Surface albedo
!!$          albedo_vis,                        & !! Surface albedo in the visible range
!!$          albedo_nir,                        & !! Surface albedo in the NIR range
!!$          veg_ratio,                         & !! This is the actual fraction of a tile covered by a vegetation canopy, i.e. this
!!$                                               !! value includes the actual leaf area index (i.e. foliar projective cover)
!!$          box_veg_ratio                        !! in contrast to veg_ratio box_veg_ratio takes into account actual cover fractions
!!$     REAL(dp),    POINTER, DIMENSION(:,:) :: & !! (nland, month)
!!$          veg_fract
!!$     REAL(dp), POINTER, DIMENSION(:)      :: & !! (nland)
!!$          veg_ratio_max,                     & !! maximal fraction of grid box covered by leaves, assuming infinite LAI
!!$          rock_fract,                        & !! fraction of grid box not suitable for plants
!!$          elevation,                         & !! Mean land grid cell elevation (m)
!!$          oro_std_dev,                       & !! Standard deviation of orography
!!$          forest_fract,                      & !! Forest fraction (TEMPORARY!)
!!$          area,                              & !! Area of land fraction of grid cell (m^2)
!!$          swdown_acc,                        & !! accumulated solar downward radiation
!!$          swdown_reflect_acc                   !! accumulated reflected solar radiation
!!$     LOGICAL, POINTER, DIMENSION(:,:) ::     & !! (nland,  ntiles)
!!$          is_bare_soil,     &     !! Tile is bare soil (but not glacier)?
!!$          is_vegetation,    &     !! Tile is vegetation (includes all natural and anthropogenically used types)
!!$          is_C4vegetation,  &     !! Tile is vegetation with C4 photosynthesis
!!$          is_naturalVeg,    &     !! Tile is natural vegetation
!!$          is_forest,        &     !! Tile is forest
!!$          is_grass,         &     !! Tile is covered with grassland          
!!$          is_pasture,       &     !! Tile is covered with pasture
!!$          is_crop,          &     !! Tile is covered with crops
!!$          is_glacier,       &     !! Tile is glacier?
!!$          is_lake,          &     !! Tile is lake?
!!$          is_present              !! Is tile present, i.e. should it be handled by the land model?
!!$  END TYPE land_surface_type
!!$  PUBLIC ::  land_surface_type
!!$
!!$  TYPE land_surface_diag_type
!!$     REAL(dp), POINTER, DIMENSION(:)     :: & !! (nland)
!!$          albedo
!!$  END TYPE land_surface_diag_type
!!$
!!$  TYPE albedo_options_type
!!$     LOGICAL :: useAlbedocanopy
!!$     LOGICAL :: UseSnowAge             !! Use snow age within albedo model?
!!$  END TYPE albedo_options_type
!!$  PUBLIC :: albedo_options_type
!!$  TYPE(albedo_options_type), SAVE :: albedo_options
!!$
!!$  TYPE albedo_params_type          !vg: Some arrays of the land_surface_type should go in here
!!$     REAL(dp) :: dummy             !vg: just to have an entry to enable compilation 
!!$  END TYPE albedo_params_type
!!$  PUBLIC :: albedo_params_type
!!$  TYPE(albedo_params_type), SAVE :: albedo_params

  !! The following parameter fract_small is the smallest allowed value for cover_fract and veg_ratio_max
!!$  REAL(dp), PARAMETER :: fract_small = 1.e-10_dp       !! very small fraction (e.g. minimum value of cover_fract)

  !! Instead, the parameter box_fract_small below is used as an estimate of the smallest box-cover-fraction of PFTs -- this  
  !! should not be confused with fract_small, the smallest value for cover_fract, which refers to canopy-area instead of   
  !! box-area. The relation between these two is obtained from
  !!         box-cover-fraction=veg_ratio_max * (1-exp(-a*LAI_max)) * cover_fract. 
  !! Since the clumping term is even for a small LAI_max=0.5 larger than 0.1, and because veg_ratio_max and cover_fract are
  !! prepared in the initial files such that they are larger than fract_small, box-cover-fraction > 0.1*fract_small*fract_small. 
  !! This explains the choice of box_fract_small as
!!$  real(dp),parameter :: box_fract_small = 0.1*fract_small*fract_small !! Estimate of smallest value for box-cover-fractions of PFTs

!!$  PUBLIC :: fract_small,box_fract_small

!!$  REAL(dp), POINTER, SAVE :: init_cover_fract(:,:)
!!$  REAL(dp), POINTER, SAVE :: init_cover_fract_pot(:,:)
!!$  REAL(dp), POINTER, SAVE :: init_veg_ratio_max(:)
!!$  PUBLIC :: init_cover_fract, init_cover_fract_pot, init_veg_ratio_max

  PRIVATE

!!$  INTEGER, SAVE :: nlct = -1        !! Number of land cover types
!!$  INTEGER, SAVE :: ntiles = -1      !! Maximum number of land cover types actually used (<= nclt)

!!$  TYPE(land_surface_diag_type), SAVE    :: land_surface_diag

!!$  LOGICAL, SAVE :: module_initialized = .FALSE.

CONTAINS

  !----------------------------------------------------------------------------
  SUBROUTINE update_albedo_echam5(kidx, ntiles, is_glacier, forest_fract, &
          surface_temperature, snow_fract, albedo_background, canopy_snow_fract, lai, &
          albedo_echam5)

! This routine is coded for ECHAM5 compatibility!

    USE mo_physical_constants,     ONLY: tmelt

    INTEGER, INTENT(in) :: kidx
    INTEGER, INTENT(in) :: ntiles
!!$    TYPE(lctlib_type),       INTENT(in)    :: lctlib
!!$    TYPE(land_surface_type), INTENT(inout) :: land_surface

    LOGICAL, INTENT(in), DIMENSION(kidx,ntiles) :: &
         is_glacier
    REAL(dp), INTENT(in), DIMENSION(kidx) ::   &
         forest_fract,        &
         surface_temperature
    REAL(dp), INTENT(in), DIMENSION(kidx, ntiles) ::   &
         snow_fract,          &                    ! Fraction of snow covered ground
         albedo_background,   &                    ! background albedo
         canopy_snow_fract,   &                    ! Fraction of snow covered canopy
         lai

    REAL(dp), INTENT(out), DIMENSION(kidx, ntiles) ::   &
         albedo_echam5

    ! Local variables
    REAL(dp)    :: min_temp_snow_albedo             ! Temperature threshold below which maximum snow albedo is used
    REAL(dp)    :: snow_albedo_ground(kidx,ntiles)  ! Temperature dependend snow albedo over ground
    REAL(dp)    :: sky_view_fract(kidx,ntiles)      ! Fraction of bare ground below canopy for albedo calculation
    REAL(dp)    :: forest_view_fract(kidx,ntiles)   ! Factor for albedo from canopy (uses forest fraction and sky view factor)
    INTEGER :: i, j, itile
    REAL(dp), PARAMETER :: SkyViewFactor = 1.0_dp        !! Constant in albedo calculation
    REAL(dp), PARAMETER :: AlbedoCanopySnow = 0.20_dp    !! Albedo of snow covered canopy
    REAL(dp), PARAMETER :: AlbedoSnowMin = 0.40_dp       !! Albedo of at tmelt
    REAL(dp), PARAMETER :: AlbedoSnowMax = 0.80_dp       !! Albedo of snow at low temperature
!!$    INTEGER :: ctype(kidx)

    min_temp_snow_albedo = tmelt - 5.0_dp
    snow_albedo_ground = 0._dp

    do itile=1,ntiles
      do i=1,kidx
        ! Temperature dependend snow albedo over ground
        IF (surface_temperature(i) >= tmelt) THEN
          snow_albedo_ground(i,itile) = AlbedoSnowMin
        ELSE IF (surface_temperature(i) < min_temp_snow_albedo) THEN
          snow_albedo_ground(i,itile) = AlbedoSnowMax
        ELSE
         snow_albedo_ground(i,itile) = AlbedoSnowMin + &
                     (tmelt - surface_temperature(i)) * &
                     (AlbedoSnowMax - AlbedoSnowMin) / &
                     (tmelt - min_temp_snow_albedo)
        END IF
      END DO
    ENDDO

    sky_view_fract(:,:) = EXP(-SkyViewFactor * MAX(lai(:,:),2.0_dp)) ! the value 2.0 should be replaced by StemArea

    ! calculate fraction for which albedo is computed from canopy
    forest_view_fract(:,:) = SPREAD(forest_fract(:), DIM=2, NCOPIES=ntiles) * &
                             (1._dp - sky_view_fract(:,:))

    WHERE (is_glacier(:,:))
       albedo_echam5(:,:) = snow_albedo_ground
    ELSEWHERE
       ! albedo = weighted mean of albedo of ground below canopy and albedo of canopy
       albedo_echam5(:,:) = &
            MAX(((1._dp - forest_view_fract(:,:)) * (snow_fract(:,:) * &
            snow_albedo_ground(:,:) + (1._dp - snow_fract(:,:)) * albedo_background(:,:)) + &
            forest_view_fract(:,:) * (canopy_snow_fract(:,:) * AlbedoCanopySnow + &
            (1._dp - canopy_snow_fract(:,:)) * albedo_background(:,:))), albedo_background(:,:))
    END WHERE

  END SUBROUTINE update_albedo_echam5

  !=================================================================================================

!!$  SUBROUTINE update_albedo(lctlib,                                   &
!!$          l_inquire, l_trigrad, l_trigradm1, l_start,                &
!!$          l_standalone,                                              &
!!$          nidx, ntiles, cover_type, is_glacier, is_forest,           &
!!$          veg_ratio_max, radiation_net_vis, radiation_net_nir,       &
!!$          snow_age, czenith,                                         &
!!$          surface_temperature, snow_fract,                           &
!!$          background_albedo_soil_vis, background_albedo_soil_nir,    &
!!$          background_albedo_veg_vis, background_albedo_veg_nir,      &
!!$          lai, canopy_snow_fract,                                    &
!!$          albedo_vis, albedo_nir, albedo)
!!$
!!$    TYPE(lctlib_type),       INTENT(in)    :: lctlib
!!$    LOGICAL, INTENT(in)  ::  l_inquire              ! trigger model initialisation
!!$    LOGICAL, INTENT(in)  ::  l_trigrad              ! trigger full radiation time step
!!$    LOGICAL, INTENT(in)  ::  l_trigradm1            ! trigger one time step before full radiation time step
!!$    LOGICAL, INTENT(in)  ::  l_start                ! trigger first time step after model initialisation
!!$    LOGICAL, INTENT(in)  ::  l_standalone           ! model runs without atmosphere
!!$    INTEGER, INTENT(in)  ::  nidx                   ! domain
!!$    INTEGER, INTENT(in)  ::  ntiles                 ! number of tiles
!!$    INTEGER, INTENT(in), DIMENSION(nidx,ntiles) ::          &
!!$         cover_type                                 ! number of plant functional type (PFT)
!!$    LOGICAL, INTENT(in), DIMENSION(nidx,ntiles) ::          &
!!$         is_glacier,                   &            ! glacier flag
!!$         is_forest                                  ! forest flag
!!$    REAL(dp), INTENT(in), DIMENSION(nidx)  ::                   &
!!$         veg_ratio_max,                &            ! maximal fraction of the grid boy covered by vegetation
!!$         radiation_net_vis,            &            ! net solar radiation in the visible range [W/m2]
!!$         radiation_net_nir,            &            ! net solar radiation in the NIR range [W/m2]
!!$         snow_age,                     &            ! non-dimensional age of snow
!!$         czenith                                    ! cosine of zenith angle
!!$    REAL(dp), INTENT(in), DIMENSION(nidx,ntiles)  ::            &
!!$         surface_temperature,          &            !
!!$         snow_fract,                   &            ! fraction of snow covered ground
!!$         background_albedo_soil_vis,   &            ! background albedo vegetation NIR
!!$         background_albedo_soil_nir,   &            ! background albedo soil NIR
!!$         background_albedo_veg_vis,    &            ! background albedo vegetation visible
!!$         background_albedo_veg_nir,    &            ! background albedo soil visible
!!$         lai,                          &            ! leaf area index
!!$         canopy_snow_fract                          ! fraction of snow covered canopy (forest)
!!$    REAL(dp), INTENT(inout), DIMENSION(nidx,ntiles)  ::         &
!!$         albedo_vis,                   &            ! albedo of the visible range
!!$         albedo_nir,                   &            ! albedo of the NIR range
!!$         albedo                                     ! albedo of the whole solar range
!!$
!!$    IF (albedo_options%UseSnowAge) THEN
!!$       CALL update_albedo_snow_age(lctlib,                           &
!!$          l_inquire, l_trigrad, l_trigradm1, l_start,                &
!!$          l_standalone,                                              &
!!$          nidx, ntiles, cover_type, is_glacier, is_forest,           &
!!$          veg_ratio_max, radiation_net_vis, radiation_net_nir,       &
!!$          snow_age, czenith,                                         &
!!$          surface_temperature, snow_fract,                           &
!!$          background_albedo_soil_vis, background_albedo_soil_nir,    &
!!$          background_albedo_veg_vis, background_albedo_veg_nir,      &
!!$          lai, canopy_snow_fract,                                    &
!!$          albedo_vis, albedo_nir, albedo)
!!$    ELSE
!!$       CALL update_albedo_temp(lctlib,                               &
!!$          l_inquire, l_trigrad, l_trigradm1, l_start,                &
!!$          l_standalone,                                              &
!!$          nidx, ntiles, cover_type, is_glacier, is_forest,           &
!!$          veg_ratio_max, radiation_net_vis, radiation_net_nir,       &
!!$          surface_temperature, snow_fract,                           &
!!$          background_albedo_soil_vis, background_albedo_soil_nir,    &
!!$          background_albedo_veg_vis, background_albedo_veg_nir,      &
!!$          lai, canopy_snow_fract,                                    &
!!$          albedo_vis, albedo_nir, albedo)
!!$    END IF
!!$
!!$  END SUBROUTINE update_albedo
!!$
!!$  !=================================================================================================
!!$
!!$  SUBROUTINE update_albedo_temp(lctlib,                              &
!!$          l_inquire, l_trigrad, l_trigradm1, l_start,                &
!!$          l_standalone,                                              &
!!$          nidx, ntiles, cover_type, is_glacier, is_forest,           &
!!$          veg_ratio_max, radiation_net_vis, radiation_net_nir,       &
!!$          surface_temperature, snow_fract,                           &
!!$          background_albedo_soil_vis, background_albedo_soil_nir,    &
!!$          background_albedo_veg_vis, background_albedo_veg_nir,      &
!!$          lai, canopy_snow_fract,                                    &
!!$          albedo_vis, albedo_nir, albedo)
!!$
!!$    USE mo_physical_constants,     ONLY : tmelt
!!$
!!$    TYPE(lctlib_type),       INTENT(in)    :: lctlib
!!$    LOGICAL, INTENT(in)  ::  l_inquire              ! trigger model initialisation
!!$    LOGICAL, INTENT(in)  ::  l_trigrad              ! trigger full radiation time step
!!$    LOGICAL, INTENT(in)  ::  l_trigradm1            ! trigger one time step before full radiation time step
!!$    LOGICAL, INTENT(in)  ::  l_start                ! trigger first time step after model initialisation
!!$    LOGICAL, INTENT(in)  ::  l_standalone           ! model runs without atmosphere
!!$    INTEGER, INTENT(in)  ::  nidx                   ! domain
!!$    INTEGER, INTENT(in)  ::  ntiles                 ! number of tiles
!!$    INTEGER, INTENT(in), DIMENSION(nidx,ntiles) ::          &
!!$         cover_type                                 ! number of plant functional type (PFT)
!!$    LOGICAL, INTENT(in), DIMENSION(nidx,ntiles) ::          &
!!$         is_glacier,                   &            ! glacier flag
!!$         is_forest                                  ! forest flag
!!$    REAL(dp), INTENT(in), DIMENSION(nidx)  ::                   &
!!$         veg_ratio_max,                &            ! maximal fraction of the grid boy covered by vegetation
!!$         radiation_net_vis,            &            ! net solar radiation in the visible range [W/m2]
!!$         radiation_net_nir                          ! net solar radiation in the NIR range [W/m2]
!!$    REAL(dp), INTENT(in), DIMENSION(nidx,ntiles)  ::            &
!!$         surface_temperature,          &            !
!!$         snow_fract,                   &            ! fraction of snow covered ground
!!$         background_albedo_soil_vis,   &            ! background albedo vegetation NIR
!!$         background_albedo_soil_nir,   &            ! background albedo soil NIR
!!$         background_albedo_veg_vis,    &            ! background albedo vegetation visible
!!$         background_albedo_veg_nir,    &            ! background albedo soil visible
!!$         lai,                          &            ! leaf area index
!!$         canopy_snow_fract                          ! fraction of snow covered canopy (forest)
!!$    REAL(dp), INTENT(inout), DIMENSION(nidx,ntiles)  ::         &
!!$         albedo_vis,                   &            ! albedo of the visible range
!!$         albedo_nir,                   &            ! albedo of the NIR range
!!$         albedo                                     ! albedo of the whole solar range
!!$
!!$    ! parameters
!!$    REAL(dp), PARAMETER :: AlbedoCanopySnow = 0.25_dp    !! Albedo of snow covered canopy
!!$    REAL(dp), PARAMETER :: AlbedoGlacierVisMin = 0.78_dp !! Albedo of glacier in the visible range at the melting point
!!$    REAL(dp), PARAMETER :: AlbedoGlacierVisMax = 0.9_dp  !! Albedo of glacier in the visible range at hard frost
!!$    REAL(dp), PARAMETER :: AlbedoGlacierNirMin = 0.44_dp !! Albedo of glacier in the NIR range at at the melting point
!!$    REAL(dp), PARAMETER :: AlbedoGlacierNirMax = 0.8_dp  !! Albedo of glacier in the NIR range at hard frost
!!$
!!$    ! Local variables
!!$    REAL(dp)  ::  background_albedo_vis                ! albedo without snow in the visible range
!!$    REAL(dp)  ::  background_albedo_nir                ! albedo without snow in the NIR range
!!$    REAL(dp)  ::  background_albedo_canopy_vis         ! albedo without snow of the canopy in the visible range
!!$    REAL(dp)  ::  background_albedo_canopy_nir         ! albedo without snow of the canopy in the NIR range
!!$    REAL(dp)  ::  fraction_down_vis(nidx,ntiles)       ! fraction of solar downward radiation in the visible range
!!$    REAL(dp)  ::  min_temp_snow_albedo                 ! temperature threshold below which maximum snow albedo is used
!!$    REAL(dp)  ::  snow_albedo_soil_vis(nidx,ntiles)    ! albedo(temp.) of snow (covering the soil) in the visible range
!!$    REAL(dp)  ::  snow_albedo_soil_nir(nidx,ntiles)    ! albedo(temp.) of snow (covering the soil) in the NIR range
!!$    REAL(dp)  ::  sky_view_fract                       ! fraction of bare ground below canopy for albedo calculation
!!$    REAL(dp)  ::  sky_view_fract_stem                  ! fraction added to snow covered canopy due to stem area
!!$    LOGICAL  ::  l_rad(nidx)                           ! flag to indicate if it is day or night
!!$    INTEGER  ::  i,itile
!!$
!!$    min_temp_snow_albedo = tmelt - 5.0_dp
!!$    snow_albedo_soil_vis(:,:) = 0.0_dp
!!$    snow_albedo_soil_nir(:,:) = 0.0_dp
!!$
!!$    IF (l_standalone .AND. l_start) THEN
!!$       albedo(:,:) = 0.2_dp
!!$       albedo_vis(:,:) = 0.1_dp
!!$       albedo_nir(:,:) = 0.3_dp
!!$    END IF
!!$
!!$    l_rad(:) = .FALSE.
!!$    IF (l_inquire) THEN
!!$      l_rad(:) = .TRUE.
!!$    ELSE
!!$      DO i = 1,nidx
!!$        IF (radiation_net_vis(i) + radiation_net_nir(i) > 1.0e-09_dp) l_rad(i) = .TRUE.
!!$      END DO
!!$    END IF
!!$
!!$    ! calculate albedo of snow
!!$    IF (l_standalone .OR. l_trigradm1 .OR. l_inquire) THEN
!!$    DO itile = 1,ntiles
!!$      DO i = 1,nidx
!!$        IF (surface_temperature(i,itile) >= tmelt) THEN
!!$          snow_albedo_soil_vis(i,itile) = lctlib%AlbedoSnowVisMin(cover_type(i,itile))
!!$          snow_albedo_soil_nir(i,itile) = lctlib%AlbedoSnowNirMin(cover_type(i,itile))
!!$        ELSE IF (surface_temperature(i,itile) < min_temp_snow_albedo) THEN
!!$          snow_albedo_soil_vis(i,itile) = lctlib%AlbedoSnowVisMax(cover_type(i,itile))
!!$          snow_albedo_soil_nir(i,itile) = lctlib%AlbedoSnowNirMax(cover_type(i,itile))
!!$        ELSE
!!$          snow_albedo_soil_vis(i,itile) = lctlib%AlbedoSnowVisMin(cover_type(i,itile)) +             &
!!$            (tmelt - surface_temperature(i,itile)) * (lctlib%AlbedoSnowVisMax(cover_type(i,itile)) - &
!!$            lctlib%AlbedoSnowVisMin(cover_type(i,itile))) / (tmelt - min_temp_snow_albedo)
!!$          snow_albedo_soil_nir(i,itile) = lctlib%AlbedoSnowNirMin(cover_type(i,itile)) +             &
!!$            (tmelt - surface_temperature(i,itile)) * (lctlib%AlbedoSnowNirMax(cover_type(i,itile)) - &
!!$            lctlib%AlbedoSnowNirMin(cover_type(i,itile))) / (tmelt - min_temp_snow_albedo)
!!$        END IF
!!$      END DO
!!$    END DO
!!$
!!$    ! calculate new glacier albedo (visible and NIR) only at grid points with solar radiation or at model initialisation
!!$    DO itile = 1,ntiles
!!$      DO i = 1,nidx
!!$        IF (l_rad(i)) THEN
!!$          IF (surface_temperature(i,itile) >= tmelt) THEN
!!$            albedo_vis(i,itile) = AlbedoGlacierVisMin
!!$            albedo_nir(i,itile) = AlbedoGlacierNirMin          
!!$          ELSE IF (surface_temperature(i,itile) < min_temp_snow_albedo) THEN
!!$            albedo_vis(i,itile) = AlbedoGlacierVisMax
!!$            albedo_nir(i,itile) = AlbedoGlacierNirMax
!!$          ELSE
!!$            albedo_vis(i,itile) = AlbedoGlacierVisMin + &
!!$              (tmelt - surface_temperature(i,itile)) * (AlbedoGlacierVisMax - AlbedoGlacierVisMin) / &
!!$              (tmelt - min_temp_snow_albedo)
!!$            albedo_nir(i,itile) = AlbedoGlacierNirMin + &
!!$              (tmelt - surface_temperature(i,itile)) * (AlbedoGlacierNirMax - AlbedoGlacierNirMin) / &
!!$              (tmelt - min_temp_snow_albedo)
!!$          END IF
!!$        END IF
!!$      END DO
!!$    END DO
!!$    END IF
!!$    ! diagnose albedo of glaciers (whole solar spectral range) only at grid points with solar radiation or at model initialisation
!!$    IF (l_standalone .OR. l_trigrad .AND. .NOT. l_inquire) THEN
!!$    DO itile = 1,ntiles
!!$      DO i = 1,nidx
!!$        IF (l_rad(i) .AND. is_glacier(i,itile)) THEN
!!$          IF (l_standalone) THEN
!!$            fraction_down_vis(i,itile) = radiation_net_vis(i) / (radiation_net_vis(i) + radiation_net_nir(i))
!!$          ELSE
!!$            fraction_down_vis(i,itile) = (radiation_net_vis(i) / (1.0_dp - albedo_vis(i,itile))) /   &
!!$                                         ((radiation_net_vis(i) / (1.0_dp - albedo_vis(i,itile))) +  &
!!$                                         (radiation_net_nir(i) / (1.0_dp - albedo_nir(i,itile))))
!!$          ENDIF
!!$          albedo(i,itile) = fraction_down_vis(i,itile) * albedo_vis(i,itile) +  &
!!$                            (1.0_dp - fraction_down_vis(i,itile)) * albedo_nir(i,itile)
!!$        END IF
!!$      END DO
!!$    END DO
!!$    END IF
!!$    IF (l_inquire) THEN
!!$    DO itile = 1,ntiles
!!$      DO i = 1,nidx
!!$        IF (is_glacier(i,itile)) THEN
!!$          albedo(i,itile) = 0.5_dp * albedo_vis(i,itile) + 0.5_dp * albedo_nir(i,itile)
!!$        END IF
!!$      END DO
!!$    END DO
!!$    END IF
!!$    ! calculate new albedo (visible and NIR) only at grid points with solar radiation or at model initialisation
!!$    IF (l_standalone .OR. l_trigradm1 .OR. l_inquire) THEN
!!$    DO itile = 1,ntiles
!!$      DO i = 1,nidx
!!$        ! albedo of the canopy (vegetation)
!!$        IF (albedo_options%useAlbedocanopy) THEN
!!$          background_albedo_canopy_vis = background_albedo_veg_vis(i,itile)
!!$          background_albedo_canopy_nir = background_albedo_veg_nir(i,itile)
!!$        ELSE
!!$          background_albedo_canopy_vis = lctlib%AlbedoCanopyVIS(cover_type(i,itile))
!!$          background_albedo_canopy_nir = lctlib%AlbedoCanopyNIR(cover_type(i,itile))
!!$        END IF
!!$        ! fraction for which albedo is computed from canopy
!!$        sky_view_fract = 1.0_dp - (veg_ratio_max(i) * &
!!$                         (1.0_dp - EXP(-0.5_dp * lai(i,itile))))
!!$        ! fraction which is added to the (snow covered) canopy fraction due to stem area
!!$        sky_view_fract_stem = sky_view_fract - (1.0_dp - (veg_ratio_max(i) * &
!!$                              (1.0_dp - EXP(-0.5_dp * (lai(i,itile) +        &
!!$                              lctlib%StemArea(cover_type(i,itile)))))))
!!$
!!$        IF (l_rad(i) .AND. .NOT. is_glacier(i,itile)) THEN
!!$          background_albedo_vis = (1.0_dp - sky_view_fract) * background_albedo_canopy_vis + &
!!$                                  sky_view_fract * background_albedo_soil_vis(i,itile)
!!$          background_albedo_nir = (1.0_dp - sky_view_fract) * background_albedo_canopy_nir + &
!!$                                  sky_view_fract * background_albedo_soil_nir(i,itile)
!!$          ! albedo of forests = weighted mean of albedo of ground below canopy and albedo of canopy
!!$          IF (is_forest(i,itile)) THEN
!!$            albedo_vis(i,itile) =                                                                                 &
!!$              MAX((sky_view_fract - sky_view_fract_stem) * (snow_fract(i,itile) * snow_albedo_soil_vis(i,itile) + &
!!$              (1.0_dp - snow_fract(i,itile)) * background_albedo_soil_vis(i,itile)) +                             &
!!$              sky_view_fract_stem  * (canopy_snow_fract(i,itile) * AlbedoCanopySnow +                             &
!!$              (1.0_dp - canopy_snow_fract(i,itile)) * background_albedo_soil_vis(i,itile)) +                      &
!!$              (1.0_dp - sky_view_fract) * (canopy_snow_fract(i,itile) * AlbedoCanopySnow +                        &
!!$              (1.0_dp - canopy_snow_fract(i,itile)) * background_albedo_canopy_vis),                              &
!!$              background_albedo_vis)
!!$            albedo_nir(i,itile) =                                                                                 &
!!$              MAX((sky_view_fract - sky_view_fract_stem) * (snow_fract(i,itile) * snow_albedo_soil_nir(i,itile) + &
!!$              (1.0_dp - snow_fract(i,itile)) * background_albedo_soil_nir(i,itile)) +                             &
!!$              sky_view_fract_stem * (canopy_snow_fract(i,itile) * AlbedoCanopySnow +                              &
!!$              (1.0_dp - canopy_snow_fract(i,itile)) * background_albedo_soil_nir(i,itile)) +                      &
!!$              (1.0_dp - sky_view_fract) * (canopy_snow_fract(i,itile) * AlbedoCanopySnow +                        &
!!$              (1.0_dp - canopy_snow_fract(i,itile)) * background_albedo_canopy_nir),                              &
!!$              background_albedo_nir)
!!$          ELSE
!!$            albedo_vis(i,itile) =                                              &
!!$              MAX(snow_fract(i,itile) * snow_albedo_soil_vis(i,itile) +        &
!!$              (1.0_dp - snow_fract(i,itile)) * background_albedo_vis,          &
!!$              background_albedo_vis)
!!$            albedo_nir(i,itile) =                                              &
!!$              MAX(snow_fract(i,itile) * snow_albedo_soil_nir(i,itile) +        &
!!$              (1.0_dp - snow_fract(i,itile)) * background_albedo_nir,          &
!!$              background_albedo_nir)
!!$          END IF
!!$        END IF
!!$      END DO
!!$    END DO
!!$    END IF
!!$    ! diagnose albedo (whole solar spectral range) only at grid points with solar radiation or at model initialisation
!!$    IF (l_standalone .OR. l_trigrad .AND. .NOT. l_inquire) THEN
!!$    DO itile = 1,ntiles
!!$      DO i = 1,nidx
!!$        IF (l_rad(i) .AND. .NOT. is_glacier(i,itile)) THEN
!!$          IF (l_standalone) THEN
!!$            fraction_down_vis(i,itile) = radiation_net_vis(i) / (radiation_net_vis(i) + radiation_net_nir(i))
!!$          ELSE
!!$            fraction_down_vis(i,itile) = (radiation_net_vis(i) / (1.0_dp - albedo_vis(i,itile))) /  &
!!$                                         ((radiation_net_vis(i) / (1.0_dp - albedo_vis(i,itile))) + &
!!$                                         (radiation_net_nir(i) / (1.0_dp - albedo_nir(i,itile))))
!!$          END IF
!!$          albedo(i,itile) = fraction_down_vis(i,itile) * albedo_vis(i,itile) +                    &
!!$                            (1.0_dp - fraction_down_vis(i,itile)) * albedo_nir(i,itile)
!!$        END IF
!!$      END DO
!!$    END DO
!!$    END IF
!!$    IF (l_inquire) THEN
!!$    DO itile = 1,ntiles
!!$      DO i = 1,nidx
!!$        IF (.NOT. is_glacier(i,itile)) THEN
!!$          albedo(i,itile) = 0.5_dp * albedo_vis(i,itile) + 0.5_dp * albedo_nir(i,itile)
!!$        END IF
!!$      END DO
!!$    END DO
!!$    END IF
!!$
!!$  END SUBROUTINE update_albedo_temp
!!$
!!$  !=================================================================================================
!!$
!!$  SUBROUTINE update_albedo_snow_age(lctlib,                          &
!!$          l_inquire, l_trigrad, l_trigradm1, l_start,                &
!!$          l_standalone,                                              &
!!$          nidx, ntiles, cover_type, is_glacier, is_forest,           &
!!$          veg_ratio_max, radiation_net_vis, radiation_net_nir,       &
!!$          snow_age, czenith,                                         &
!!$          surface_temperature, snow_fract,                           &
!!$          background_albedo_soil_vis, background_albedo_soil_nir,    &
!!$          background_albedo_veg_vis, background_albedo_veg_nir,      &
!!$          lai, canopy_snow_fract,                                    &
!!$          albedo_vis, albedo_nir, albedo)
!!$
!!$    USE mo_physical_constants,     ONLY : tmelt
!!$
!!$    TYPE(lctlib_type),       INTENT(in)    :: lctlib
!!$    LOGICAL, INTENT(in)  ::  l_inquire              ! trigger model initialisation
!!$    LOGICAL, INTENT(in)  ::  l_trigrad              ! trigger full radiation time step
!!$    LOGICAL, INTENT(in)  ::  l_trigradm1            ! trigger one time step before full radiation time step
!!$    LOGICAL, INTENT(in)  ::  l_start                ! trigger first time step after model initialisation
!!$    LOGICAL, INTENT(in)  ::  l_standalone           ! model runs without atmosphere
!!$    INTEGER, INTENT(in)  ::  nidx                   ! domain
!!$    INTEGER, INTENT(in)  ::  ntiles                 ! number of tiles
!!$    INTEGER, INTENT(in), DIMENSION(nidx,ntiles) ::          &
!!$         cover_type                                 ! number of plant functional type (PFT)
!!$    LOGICAL, INTENT(in), DIMENSION(nidx,ntiles) ::          &
!!$         is_glacier,                   &            ! glacier flag
!!$         is_forest                                  ! forest flag
!!$    REAL(dp), INTENT(in), DIMENSION(nidx)  ::                   &
!!$         veg_ratio_max,                &            ! maximal fraction of the grid boy covered by vegetation
!!$         radiation_net_vis,            &            ! net solar radiation in the visible range [W/m2]
!!$         radiation_net_nir,            &            ! net solar radiation in the NIR range [W/m2]
!!$         snow_age,                     &            ! non-dimensional age of snow
!!$         czenith                                    ! cosine of zenith angle
!!$    REAL(dp), INTENT(in), DIMENSION(nidx,ntiles)  ::            &
!!$         surface_temperature,          &            !
!!$         snow_fract,                   &            ! fraction of snow covered ground
!!$         background_albedo_soil_vis,   &            ! background albedo vegetation NIR
!!$         background_albedo_soil_nir,   &            ! background albedo soil NIR
!!$         background_albedo_veg_vis,    &            ! background albedo vegetation visible
!!$         background_albedo_veg_nir,    &            ! background albedo soil visible
!!$         lai,                          &            ! leaf area index
!!$         canopy_snow_fract                          ! fraction of snow covered canopy (forest)
!!$    REAL(dp), INTENT(inout), DIMENSION(nidx,ntiles)  ::         &
!!$         albedo_vis,                   &            ! albedo of the visible range
!!$         albedo_nir,                   &            ! albedo of the NIR range
!!$         albedo                                     ! albedo of the whole solar range
!!$
!!$    ! parameters
!!$    REAL(dp), PARAMETER :: AlbedoCanopySnow = 0.20_dp    !! Albedo of snow covered canopy
!!$    REAL(dp), PARAMETER :: AlbedoGlacierVisMin = 0.8_dp  !! Albedo of glacier in the visible range at the melting point
!!$    REAL(dp), PARAMETER :: AlbedoGlacierVisMax = 0.9_dp  !! Albedo of glacier in the visible range at hard frost
!!$    REAL(dp), PARAMETER :: AlbedoGlacierNirMin = 0.6_dp  !! Albedo of glacier in the NIR range at at the melting point
!!$    REAL(dp), PARAMETER :: AlbedoGlacierNirMax = 0.8_dp  !! Albedo of glacier in the NIR range at hard frost
!!$    REAL(dp), PARAMETER :: TempAlbedoGlacierMax = 1.0_dp !! Maximum glacier albedo at this temperature below melting point of H2O
!!$    REAL(dp), PARAMETER :: AlbedoSnowVisMax = 0.95_dp    !! Maximum albedo of fresh snow in the visible range
!!$    REAL(dp), PARAMETER :: AlbedoSnowNirMax = 0.65_dp    !! Maximum albedo of fresh snow in the NIR range
!!$    REAL(dp), PARAMETER :: AlbedoSnowVisAge = 0.2_dp     !! Maximal rel. reduction of snow albedo by aging in the visible range
!!$    REAL(dp), PARAMETER :: AlbedoSnowNirAge = 0.5_dp     !! Maximal rel. reduction of snow albedo by aging in the NIR range
!!$    REAL(dp), PARAMETER :: AlbedoSnowAngle = 0.4_dp      !! Maximal rel. reduction of snow absorption by large solar zenith angle
!!$    REAL(dp), PARAMETER :: ZenithAngleFactor = 2._dp     !! Factor in solar zenith angle dependence of snow albedo
!!$                                                         !! (the increase of snow albedo is the higher this factor
!!$    REAL(dp), PARAMETER :: SkyViewFactor = 1.0_dp        !! Constant in calculating the sky view fraction depending on
!!$                                                         !! lai and stem area
!!$    ! Local variables
!!$    REAL(dp)  ::  background_albedo_vis                ! albedo without snow in the visible range
!!$    REAL(dp)  ::  background_albedo_nir                ! albedo without snow in the NIR range
!!$    REAL(dp)  ::  background_albedo_canopy_vis         ! albedo without snow of the canopy in the visible range
!!$    REAL(dp)  ::  background_albedo_canopy_nir         ! albedo without snow of the canopy in the NIR range
!!$    REAL(dp)  ::  fraction_down_vis(nidx,ntiles)       ! fraction of solar downward radiation in the visible range
!!$    REAL(dp)  ::  min_temp_glacier_albedo              ! temperature threshold below which maximum glacier albedo is used
!!$    REAL(dp)  ::  snow_albedo_soil_vis(nidx,ntiles)    ! albedo(temp.) of snow (covering the soil) in the visible range
!!$    REAL(dp)  ::  snow_albedo_soil_nir(nidx,ntiles)    ! albedo(temp.) of snow (covering the soil) in the NIR range
!!$    REAL(dp)  ::  snow_age_factor                      ! snow aging factor
!!$    REAL(dp)  ::  snow_albedo_angle_factor             ! function of solar zenith angle (for albedo of snow on land)
!!$    REAL(dp)  ::  sky_view_fract                       ! fraction of bare ground below canopy for albedo calculation
!!$    REAL(dp)  ::  sky_view_fract_stem                  ! fraction added to snow covered canopy due to stem area
!!$    LOGICAL  ::  l_rad(nidx)                           ! flag to indicate if it is day or night
!!$    INTEGER  ::  i,itile
!!$
!!$    min_temp_glacier_albedo = tmelt - TempAlbedoGlacierMax ! ztalb
!!$    snow_albedo_soil_vis(:,:) = 0.0_dp
!!$    snow_albedo_soil_nir(:,:) = 0.0_dp
!!$
!!$    IF (l_standalone .AND. l_start) THEN
!!$       albedo(:,:) = 0.2_dp
!!$       albedo_vis(:,:) = 0.1_dp
!!$       albedo_nir(:,:) = 0.3_dp
!!$    END IF
!!$
!!$    l_rad(:) = .FALSE.
!!$    IF (l_inquire) THEN
!!$      l_rad(:) = .TRUE.
!!$    ELSE
!!$      DO i = 1,nidx
!!$        IF (radiation_net_vis(i) + radiation_net_nir(i) > 1.0e-09_dp) l_rad(i) = .TRUE.
!!$      END DO
!!$    END IF
!!$
!!$    ! calculate albedo of snow
!!$    IF (l_standalone .OR. l_trigradm1 .OR. l_inquire) THEN
!!$    DO itile = 1,ntiles
!!$      DO i = 1,nidx
!!$        IF (snow_fract(i,itile) > EPSILON(1._dp)) THEN
!!$          snow_age_factor = snow_age(i) / (1._dp + snow_age(i))
!!$          IF (.NOT. l_inquire) THEN
!!$            IF (czenith(i) < 0.5_dp) THEN
!!$              snow_albedo_angle_factor = ((1.0_dp + ZenithAngleFactor)/ &
!!$                                         (1.0_dp + 2.0_dp * ZenithAngleFactor * &
!!$                                         czenith(i)) - 1.0_dp) / ZenithAngleFactor
!!$            ELSE
!!$              snow_albedo_angle_factor = 0._dp
!!$            END IF
!!$          ELSE
!!$            snow_albedo_angle_factor = 0._dp
!!$          END IF
!!$          snow_albedo_soil_vis(i,itile) = AlbedoSnowVisMax * &
!!$            (1._dp - AlbedoSnowVisAge * snow_age_factor)
!!$          snow_albedo_soil_vis(i,itile) = snow_albedo_soil_vis(i,itile) + &
!!$            AlbedoSnowAngle * snow_albedo_angle_factor * &
!!$            (1._dp - snow_albedo_soil_vis(i,itile))
!!$          snow_albedo_soil_nir(i,itile) = AlbedoSnowNirMax * &
!!$            (1._dp - AlbedoSnowNirAge * snow_age_factor)
!!$          snow_albedo_soil_nir(i,itile) = snow_albedo_soil_nir(i,itile) + &
!!$            AlbedoSnowAngle * snow_albedo_angle_factor * &
!!$            (1._dp - snow_albedo_soil_nir(i,itile))
!!$        ELSE
!!$          snow_albedo_soil_vis(i,itile) = 0._dp
!!$          snow_albedo_soil_nir(i,itile) = 0._dp
!!$        END IF
!!$      END DO
!!$    END DO
!!$
!!$    ! calculate new glacier albedo (visible and NIR) only at grid points with solar radiation or at model initialisation
!!$    DO itile = 1,ntiles
!!$      DO i = 1,nidx
!!$        IF (l_rad(i)) THEN
!!$          IF (surface_temperature(i,itile) >= tmelt) THEN
!!$            albedo_vis(i,itile) = AlbedoGlacierVisMin
!!$            albedo_nir(i,itile) = AlbedoGlacierNirMin          
!!$          ELSE IF (surface_temperature(i,itile) < min_temp_glacier_albedo) THEN
!!$            albedo_vis(i,itile) = AlbedoGlacierVisMax
!!$            albedo_nir(i,itile) = AlbedoGlacierNirMax
!!$          ELSE
!!$            albedo_vis(i,itile) = AlbedoGlacierVisMin + &
!!$              (tmelt - surface_temperature(i,itile)) * (AlbedoGlacierVisMax - AlbedoGlacierVisMin) / &
!!$              TempAlbedoGlacierMax
!!$            albedo_nir(i,itile) = AlbedoGlacierNirMin + &
!!$              (tmelt - surface_temperature(i,itile)) * (AlbedoGlacierNirMax - AlbedoGlacierNirMin) / &
!!$              TempAlbedoGlacierMax
!!$          END IF
!!$        END IF
!!$      END DO
!!$    END DO
!!$    END IF
!!$
!!$    ! diagnose albedo of glaciers (whole solar spectral range) only at grid points with solar radiation or at model initialisation
!!$    IF (l_standalone .OR. l_trigrad .AND. .NOT. l_inquire) THEN
!!$    DO itile = 1,ntiles
!!$      DO i = 1,nidx
!!$        IF (l_rad(i) .AND. is_glacier(i,itile)) THEN
!!$          IF (l_standalone) THEN
!!$            fraction_down_vis(i,itile) = radiation_net_vis(i) / (radiation_net_vis(i) + radiation_net_nir(i))
!!$          ELSE
!!$            fraction_down_vis(i,itile) = (radiation_net_vis(i) / (1.0_dp - albedo_vis(i,itile))) /   &
!!$                                         ((radiation_net_vis(i) / (1.0_dp - albedo_vis(i,itile))) +  &
!!$                                         (radiation_net_nir(i) / (1.0_dp - albedo_nir(i,itile))))
!!$          END IF 
!!$          albedo(i,itile) = fraction_down_vis(i,itile) * albedo_vis(i,itile) +  &
!!$                            (1.0_dp - fraction_down_vis(i,itile)) * albedo_nir(i,itile)
!!$        END IF
!!$      END DO
!!$    END DO
!!$    END IF
!!$    IF (l_inquire) THEN
!!$    DO itile = 1,ntiles
!!$      DO i = 1,nidx
!!$        IF (is_glacier(i,itile)) THEN
!!$          albedo(i,itile) = 0.5_dp * albedo_vis(i,itile) + 0.5_dp * albedo_nir(i,itile)
!!$        END IF
!!$      END DO
!!$    END DO
!!$    END IF
!!$
!!$    ! calculate new albedo (visible and NIR) only at grid points with solar radiation or at model initialisation
!!$    IF (l_standalone .OR. l_trigradm1 .OR. l_inquire) THEN
!!$    DO itile = 1,ntiles
!!$      DO i = 1,nidx
!!$        ! albedo of the canopy (vegetation)
!!$        IF (albedo_options%useAlbedocanopy) THEN
!!$          background_albedo_canopy_vis = background_albedo_veg_vis(i,itile)
!!$          background_albedo_canopy_nir = background_albedo_veg_nir(i,itile)
!!$        ELSE
!!$          background_albedo_canopy_vis = lctlib%AlbedoCanopyVIS(cover_type(i,itile))
!!$          background_albedo_canopy_nir = lctlib%AlbedoCanopyNIR(cover_type(i,itile))
!!$        END IF
!!$        ! fraction for which albedo is computed from canopy
!!$        sky_view_fract = 1.0_dp - (veg_ratio_max(i) * &
!!$                         (1.0_dp - EXP(-SkyViewFactor * lai(i,itile))))
!!$        ! fraction which is added to the (snow covered) canopy fraction due to stem area
!!$        sky_view_fract_stem = sky_view_fract - (1.0_dp - (veg_ratio_max(i) * &
!!$                              (1.0_dp - EXP(-SkyViewFactor *                 &
!!$                              (lai(i,itile) + lctlib%StemArea(cover_type(i,itile)))))))
!!$
!!$        IF (l_rad(i) .AND. .NOT. is_glacier(i,itile)) THEN
!!$          background_albedo_vis = (1.0_dp - sky_view_fract) * background_albedo_canopy_vis + &
!!$                                  sky_view_fract * background_albedo_soil_vis(i,itile)
!!$          background_albedo_nir = (1.0_dp - sky_view_fract) * background_albedo_canopy_nir + &
!!$                                  sky_view_fract * background_albedo_soil_nir(i,itile)
!!$          ! albedo of forests = weighted mean of albedo of ground below canopy and albedo of canopy
!!$          IF (is_forest(i,itile)) THEN
!!$            albedo_vis(i,itile) =                                                                                 &
!!$              MAX((sky_view_fract - sky_view_fract_stem) * (snow_fract(i,itile) * snow_albedo_soil_vis(i,itile) + &
!!$              (1.0_dp - snow_fract(i,itile)) * background_albedo_soil_vis(i,itile)) +                             &
!!$              sky_view_fract_stem  * (canopy_snow_fract(i,itile) * AlbedoCanopySnow +                             &
!!$              (1.0_dp - canopy_snow_fract(i,itile)) * background_albedo_soil_vis(i,itile)) +                      &
!!$              (1.0_dp - sky_view_fract) * (canopy_snow_fract(i,itile) * AlbedoCanopySnow +                        &
!!$              (1.0_dp - canopy_snow_fract(i,itile)) * background_albedo_canopy_vis),                              &
!!$              background_albedo_vis)
!!$            albedo_nir(i,itile) =                                                                                 &
!!$              MAX((sky_view_fract - sky_view_fract_stem) * (snow_fract(i,itile) * snow_albedo_soil_nir(i,itile) + &
!!$              (1.0_dp - snow_fract(i,itile)) * background_albedo_soil_nir(i,itile)) +                             &
!!$              sky_view_fract_stem * (canopy_snow_fract(i,itile) * AlbedoCanopySnow +                              &
!!$              (1.0_dp - canopy_snow_fract(i,itile)) * background_albedo_soil_nir(i,itile)) +                      &
!!$              (1.0_dp - sky_view_fract) * (canopy_snow_fract(i,itile) * AlbedoCanopySnow +                        &
!!$              (1.0_dp - canopy_snow_fract(i,itile)) * background_albedo_canopy_nir),                              &
!!$              background_albedo_nir)
!!$          ELSE
!!$            albedo_vis(i,itile) =                                              &
!!$              MAX(snow_fract(i,itile) * snow_albedo_soil_vis(i,itile) +        &
!!$              (1.0_dp - snow_fract(i,itile)) * background_albedo_vis,          &
!!$              background_albedo_vis)
!!$            albedo_nir(i,itile) =                                              &
!!$              MAX(snow_fract(i,itile) * snow_albedo_soil_nir(i,itile) +        &
!!$              (1.0_dp - snow_fract(i,itile)) * background_albedo_nir,          &
!!$              background_albedo_nir)
!!$          END IF
!!$        END IF
!!$      END DO
!!$    END DO
!!$    END IF
!!$    ! diagnose albedo (whole solar spectral range) only at grid points with solar radiation or at model initialisation
!!$    IF (l_standalone .OR. l_trigrad .AND. .NOT. l_inquire) THEN
!!$    DO itile = 1,ntiles
!!$      DO i = 1,nidx
!!$        IF (l_rad(i) .AND. .NOT. is_glacier(i,itile)) THEN
!!$          IF (l_standalone) THEN
!!$            fraction_down_vis(i,itile) = radiation_net_vis(i) / ( radiation_net_vis(i) + radiation_net_nir(i) )
!!$          ELSE
!!$            fraction_down_vis(i,itile) = (radiation_net_vis(i) / (1.0_dp - albedo_vis(i,itile))) /  &
!!$                                         ((radiation_net_vis(i) / (1.0_dp - albedo_vis(i,itile))) + &
!!$                                         (radiation_net_nir(i) / (1.0_dp - albedo_nir(i,itile))))
!!$          END IF
!!$          albedo(i,itile) = fraction_down_vis(i,itile) * albedo_vis(i,itile) +                    &
!!$                            (1.0_dp - fraction_down_vis(i,itile)) * albedo_nir(i,itile)
!!$        END IF
!!$      END DO
!!$    END DO
!!$    END IF
!!$    IF (l_inquire) THEN
!!$    DO itile = 1,ntiles
!!$      DO i = 1,nidx
!!$        IF (.NOT. is_glacier(i,itile)) THEN
!!$          albedo(i,itile) = 0.5_dp * albedo_vis(i,itile) + 0.5_dp * albedo_nir(i,itile)
!!$        END IF
!!$      END DO
!!$    END DO
!!$    END IF
!!$
!!$  END SUBROUTINE update_albedo_snow_age

END MODULE mo_land_surface
