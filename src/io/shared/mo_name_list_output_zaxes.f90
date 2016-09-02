!>
!! Module handling the specification of vertical axes for the output
!! module.
!!
!! @par Revision History
!! Initial implementation by R. Johanni, taken from io_vlist module.
!! Moved to a separate module: F. Prill, DWD (2014-08-12)
!! Added vertical level selection: F. Prill, DWD (2014-08-15)
!!
!! --------------------------------
!!    Details of the implementation
!! --------------------------------
!!
!! Axis creation:
!!
!! The CDI routines for the creation of the vertical axes are called
!! for each output file separately in "mo_name_list_output_init",
!! subroutine "setup_output_vlist". 
!!
!! For model level output, all axes (atmospheric, surface, soil etc.)
!! are created.  For pressure, height, and isentropic level output
!! only the respective axis is created together with a surface axis
!! for the grid information.
!!
!! Note: The vertical axes for the output of the ocean do not support
!!       level selection yet.
!!
!! Derived data type "t_level_selection":
!!
!! Vertical levels are selected via objects of the derived data type
!! "t_level_selection" (module "mo_name_list_output_types").  If such
!! a data object has been initialized for the output file, then the
!! vertical axis definitions below create CDI axis objects for the
!! selected levels only.
!! The write routines in the module "mo_name_list_output" skip levels
!! which are not part of the selection.
!!
!! Furthermore, if a "t_level_selection" object is present, then the
!! memory windows for the asynchronous one-sided MPI communication are
!! adjusted to the reduced level size
!! "output_file%level_selection%n_selected".
!!
!! Creation of level selection data:
!!
!! During the setup phase, level selection objects are created in two
!! different situations:
!!
!! a) Definition of namelist parameter "m_levels"
!!
!!    The user may specify levels and/or level ranges in the form of a
!!    (string) namelist parameter. This string is parsed with the help
!!    of the module "mo_util_string_parse" and converted into a
!!    "t_level_selection" object.
!!
!! b) Definition of namelist parameters "p_levels", "h_levels", "i_levels"
!!
!!    Each output namelist "output_nml" may specify its own range of
!!    pressure levels for output (or for height, isentropic levels as
!!    well). However, the internal post-processing routines perform
!!    vertical interpolation for the union set of all requested
!!    pressure levels of a specific domain at once (merging the
!!    information from several "output_nml" namelists. Afterwards, the
!!    write routine for the output files merely copies the respective
!!    levels. For this purpose, the described "t_level_selection"
!!    mechanism is applied as well.
!!
MODULE mo_name_list_output_zaxes

  USE ISO_C_BINDING,                        ONLY: C_SIGNED_CHAR
  USE mo_cdi,                               ONLY: CDI_UNDEFID, ZAXIS_DEPTH_BELOW_SEA, ZAXIS_GENERIC, ZAXIS_SURFACE, &
                                                & ZAXIS_ISENTROPIC, ZAXIS_ALTITUDE, ZAXIS_PRESSURE, ZAXIS_CLOUD_BASE, &
                                                & ZAXIS_CLOUD_TOP, ZAXIS_DEPTH_BELOW_LAND, ZAXIS_HEIGHT, ZAXIS_HYBRID, &
                                                & ZAXIS_HYBRID_HALF, ZAXIS_ISOTHERM_ZERO, ZAXIS_LAKE_BOTTOM, ZAXIS_MEANSEA, &
                                                & ZAXIS_MIX_LAYER, ZAXIS_REFERENCE, ZAXIS_SEDIMENT_BOTTOM_TW, ZAXIS_SNOW, &
                                                & ZAXIS_ATMOSPHERE,  &
                                                & ZAXIS_TOA, zaxisCreate, zaxisDefNumber, zaxisDefUUID, zaxisDefLevels, &
                                                & zaxisDefLbounds, zaxisDefUbounds, zaxisDefVct, zaxisDefUnits, zaxisDefNlevRef
  USE mo_cdi_constants,                     ONLY: ZA_depth_below_sea, ZA_depth_below_sea_half, ZA_GENERIC_ICE, ZA_surface, &
                                                & ZA_isentropic, ZA_altitude, ZA_pressure, ZA_cloud_base, ZA_cloud_top, &
                                                & ZA_depth_below_land, ZA_depth_below_land_p1, ZA_depth_runoff_g, &
                                                & ZA_depth_runoff_s, ZA_height_10m, ZA_height_2m, ZA_hybrid, ZA_hybrid_half, &
                                                & ZA_hybrid_half_hhl, ZA_isotherm_zero, ZA_lake_bottom, ZA_lake_bottom_half, &
                                                & ZA_meansea, ZA_mix_layer, ZA_pressure_0, ZA_pressure_400, ZA_pressure_800, &
                                                & ZA_ATMOSPHERE, ZA_PRES_FL_SFC_200, ZA_PRES_FL_200_350, ZA_PRES_FL_350_550,      &
                                                & ZA_PRES_FL_SFC_100, ZA_PRES_FL_100_245, ZA_PRES_FL_245_390, ZA_PRES_FL_390_530, &
                                                & ZA_reference, ZA_reference_half, ZA_reference_half_hhl, &
                                                & ZA_sediment_bottom_tw_half, ZA_snow, ZA_snow_half, ZA_toa, ZA_OCEAN_SEDIMENT
  USE mo_kind,                              ONLY: wp, dp
  USE mo_impl_constants,                    ONLY: zml_soil, SUCCESS
  USE mo_var_list_element,                  ONLY: level_type_ml, level_type_pl, level_type_hl,    &
    &                                             level_type_il
  USE mo_exception,                         ONLY: finish
  USE mo_name_list_output_types,            ONLY: t_output_file, t_level_selection
  USE mo_util_vgrid_types,                  ONLY: vgrid_buffer
  USE mo_vertical_coord_table,              ONLY: vct
  USE mo_math_utilities,                    ONLY: set_zlev, t_value_set, find_values_in_set
  USE mo_run_config,                        ONLY: num_lev
  USE mo_util_string_parse,                 ONLY: util_do_parse_intlist
#ifndef __NO_ICON_ATMO__
  USE mo_nh_pzlev_config,                   ONLY: nh_pzlev_config
  USE mo_nonhydrostatic_config,             ONLY: ivctype
  USE mo_lnd_nwp_config,                    ONLY: nlev_snow
#endif
#ifndef __NO_ICON_OCEAN__
  USE mo_ocean_nml,                         ONLY: n_zlev, dzlev_m,lhamocc
  USE mo_sedmnt,                            ONLY: ks, ksp, dzsed


#endif

  IMPLICIT NONE

  PRIVATE

  ! subroutines
  PUBLIC :: setup_ml_axes_atmo
  PUBLIC :: setup_pl_axis_atmo
  PUBLIC :: setup_hl_axis_atmo
  PUBLIC :: setup_il_axis_atmo
  PUBLIC :: setup_zaxes_oce
  PUBLIC :: create_mipz_level_selections
  PUBLIC :: deallocate_level_selection

  INTERFACE create_level_selection
    MODULE PROCEDURE create_level_selection_str
    MODULE PROCEDURE create_level_selection_set
  END INTERFACE

  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_name_list_output_zaxes'

CONTAINS

  ! --------------------------------------------------------------------------------------
  !> Setup of vertical axes for output module: Atmosphere component, model levels
  !
  SUBROUTINE setup_ml_axes_atmo(of)
    TYPE(t_output_file),     INTENT(INOUT)        :: of                  !< output file meta-data
    ! local variables
    REAL(dp), ALLOCATABLE           :: levels(:), lbounds(:), ubounds(:)
    INTEGER                         :: k, nlev, nlevp1, znlev_soil

#ifndef __NO_ICON_ATMO__

    of%cdiZaxisID(:) = CDI_UNDEFID ! not all are set

    nlev   = num_lev(of%log_patch_id)
    nlevp1 = num_lev(of%log_patch_id) + 1

    ! introduce temporary variable znlev_soil, since global variable
    ! nlev_soil is unknown to the I/O-Processor.
    znlev_soil = SIZE(zml_soil)


    ! --------------------------------------------------------------------------------------
    ! Definitions for single levels --------------------------------------------------------
    ! --------------------------------------------------------------------------------------

    ! surface level
    CALL define_single_level_axis(of, ZA_surface, ZAXIS_SURFACE)

    ! CLOUD BASE LEVEL
    CALL define_single_level_axis(of, ZA_cloud_base, ZAXIS_CLOUD_BASE)

    ! CLOUD TOP LEVEL
    CALL define_single_level_axis(of, ZA_cloud_top, ZAXIS_CLOUD_TOP)

    ! LEVEL of 0\deg C isotherm
    CALL define_single_level_axis(of, ZA_isotherm_zero, ZAXIS_ISOTHERM_ZERO)

    ! Define axis for output on mean sea level
    CALL define_single_level_axis(of, ZA_meansea, ZAXIS_MEANSEA)

    ! Specific soil axis for Runoff_s
    CALL define_single_layer_axis(of, ZA_depth_runoff_s, ZAXIS_DEPTH_BELOW_LAND, 0._dp, 0.1_dp, &
      &                           unit="m")

    ! Specific soil axis for Runoff_g
    CALL define_single_layer_axis(of, ZA_depth_runoff_g, ZAXIS_DEPTH_BELOW_LAND, 0.1_dp, 1._dp, &
      &                           unit="m")

    ! Specified height level above ground: 2m
    CALL define_single_level_axis(of, ZA_height_2m, ZAXIS_HEIGHT, opt_level_value=2._dp)

    ! Specified height level above ground: 10m
    CALL define_single_level_axis(of, ZA_height_10m, ZAXIS_HEIGHT, opt_level_value=10._dp)

    ! Top of atmosphere
    CALL define_single_level_axis(of, ZA_toa, ZAXIS_TOA)

    ! Bottom of sediment layer penetrated by thermal wave (interface, i.e. only typeOfFirstFixedSurface)
    CALL define_single_level_axis(of, ZA_sediment_bottom_tw_half, ZAXIS_SEDIMENT_BOTTOM_TW, opt_unit="m")

    ! Lake bottom half (interface, i.e. only typeOfFirstFixedSurface)
    CALL define_single_level_axis(of, ZA_lake_bottom_half, ZAXIS_LAKE_BOTTOM, opt_unit="m")

    ! for having ice variable in the atmosphere (like AMIP)
    of%cdiZaxisID(ZA_GENERIC_ICE) = zaxisCreate(ZAXIS_GENERIC, 1)


    ! --------------------------------------------------------------------------------------
    ! Definitions for single layers --------------------------------------------------------
    ! --------------------------------------------------------------------------------------

    ! Isobaric surface 800 hPa (layer)
    !DR Note that for this particular axis in GME and COSMO, 
    !DR   scaleFactorOfSecondFixedSurface = scaledValueOfSecondFixedSurface = missing
    !DR whereas in ICON
    !DR   scaleFactorOfSecondFixedSurface = scaledValueOfSecondFixedSurface = 0
    !DR This is because, CDI does not allow layer-type vertical axis without setting 
    !DR scaleFactorOfSecondFixedSurface and scaledValueOfSecondFixedSurface.
    !
    CALL define_single_layer_axis(of, ZA_pressure_800, ZAXIS_PRESSURE, 800._dp,   0._dp, "hPa")
    ! Isobaric surface 400 hPa (layer)
    CALL define_single_layer_axis(of, ZA_pressure_400, ZAXIS_PRESSURE, 400._dp, 800._dp, "hPa")
    ! Isobaric surface 0 hPa (layer)
    CALL define_single_layer_axis(of, ZA_pressure_0,   ZAXIS_PRESSURE,   0._dp, 400._dp, "hPa")

    ! Specific vertical axis for Lake-model ------------------------------------------------
    ! Lake bottom (we define it as a layer in order to be able to re-set
    ! either the first- or secondFixedSurfaces if necessary)
    CALL define_single_layer_axis(of, ZA_lake_bottom, ZAXIS_LAKE_BOTTOM, 0._dp, 0._dp, "m")

    ! Mixing layer (we define it as a layer in order to be able to re-set
    ! either the first- or secondFixedSurfaces if necessary)
    CALL define_single_layer_axis(of, ZA_mix_layer, ZAXIS_MIX_LAYER, 1._dp, 0._dp, "m")

    ! Volcanic ash products - Maximum total mass concentration in flight level range
    !                         defined by pressure layers
    CALL define_single_layer_axis(of, ZA_PRES_FL_SFC_200, ZAXIS_PRESSURE, 465.00_dp, 1013.25_dp, "hPa")
    CALL define_single_layer_axis(of, ZA_PRES_FL_200_350, ZAXIS_PRESSURE, 240.00_dp,  465.00_dp, "hPa")
    CALL define_single_layer_axis(of, ZA_PRES_FL_350_550, ZAXIS_PRESSURE,  91.00_dp,  240.00_dp, "hPa")
    CALL define_single_layer_axis(of, ZA_PRES_FL_SFC_100, ZAXIS_PRESSURE, 700.00_dp, 1013.25_dp, "hPa")
    CALL define_single_layer_axis(of, ZA_PRES_FL_100_245, ZAXIS_PRESSURE, 385.00_dp,  700.00_dp, "hPa")
    CALL define_single_layer_axis(of, ZA_PRES_FL_245_390, ZAXIS_PRESSURE, 200.00_dp,  385.00_dp, "hPa")
    CALL define_single_layer_axis(of, ZA_PRES_FL_390_530, ZAXIS_PRESSURE, 100.00_dp,  200.00_dp, "hPa")
    ! Volcanic ash products - Colummn integrated total mass concentration (entire atmosphere)
    CALL define_single_level_axis(of, ZA_ATMOSPHERE, ZAXIS_ATMOSPHERE)

    ! --------------------------------------------------------------------------------------
    ! Definitions for reference grids (ZAXIS_REFERENCE) ------------------------------------
    ! --------------------------------------------------------------------------------------

    ! REFERENCE
    CALL define_vertical_axis(of, ZA_reference, ZAXIS_REFERENCE, nlev,                   &
      &                           levels           = (/ ( REAL(k,dp),   k=1,nlevp1 ) /), &
      &                           opt_set_bounds   = .TRUE.,                             &
      &                           opt_number       = get_numberOfVgridUsed(ivctype),     &
      &                           opt_uuid         = vgrid_buffer(of%log_patch_id)%uuid%DATA )
    ! Define number of half levels for z-axis 
    CALL zaxisDefNlevRef(of%cdiZaxisID(ZA_reference),nlevp1)

    ! REFERENCE HALF
    CALL define_vertical_axis(of, ZA_reference_half, ZAXIS_REFERENCE, nlevp1,        &
      &                           levels   = (/ ( REAL(k,dp),   k=1,nlevp1 ) /),     &
      &                           opt_number = get_numberOfVgridUsed(ivctype),       &
      &                           opt_uuid   = vgrid_buffer(of%log_patch_id)%uuid%DATA )
    ! Define number of half levels for z-axis 
    CALL zaxisDefNlevRef(of%cdiZaxisID(ZA_reference_half),nlevp1)

    ! REFERENCE (special version for HHL)
    CALL define_vertical_axis(of, ZA_reference_half_hhl, ZAXIS_REFERENCE, nlevp1,        &
      &                           levels           = (/ ( REAL(k,dp),   k=1,nlevp1 ) /), &
      &                           opt_set_bounds   = .TRUE.,                             &
      &                           opt_set_ubounds_value = 0._dp,                         &
      &                           opt_number       = get_numberOfVgridUsed(ivctype),     &
      &                           opt_uuid         = vgrid_buffer(of%log_patch_id)%uuid%DATA )
    ! Define number of half levels for z-axis 
    CALL zaxisDefNlevRef(of%cdiZaxisID(ZA_reference_half_hhl),nlevp1)


    ! --------------------------------------------------------------------------------------
    ! Definitions for hybrid z-axes (ZAXIS_HYBRID) -----------------------------------------
    ! --------------------------------------------------------------------------------------

    ! HYBRID_LAYER
    !
    CALL define_vertical_axis(of, ZA_hybrid, ZAXIS_HYBRID, nlev,               &
      &                           levels = (/ ( REAL(k,dp),   k=1,nlevp1 ) /), &
      &                           opt_set_bounds  = .TRUE.)
    CALL zaxisDefVct(of%cdiZaxisID(ZA_hybrid), 2*nlevp1, vct(1:2*nlevp1))

    ! HYBRID
    !
    ! Note: "ZAXIS_HYBRID_HALF" is deprecated and will soon be
    ! removed from the CDI (in principle its use should be simply
    ! replaced by ZAXIS_HALF, as long as lbounds and ubounds are set
    ! correctly).
    CALL define_vertical_axis(of, ZA_hybrid_half, ZAXIS_HYBRID_HALF, nlevp1,   &
      &                           levels = (/ ( REAL(k,dp),   k=1,nlevp1 ) /))
    CALL zaxisDefVct(of%cdiZaxisID(ZA_hybrid_half), 2*nlevp1, vct(1:2*nlevp1))

    ! HYBRID (special version for HHL)
    !
    CALL define_vertical_axis(of, ZA_hybrid_half_hhl, ZAXIS_HYBRID_HALF, nlevp1, &
      &                           levels = (/ ( REAL(k,dp),   k=1,nlevp1 ) /),   &
      &                           opt_set_bounds = .TRUE.,                       &
      &                           opt_set_ubounds_value  = 0._dp)
    CALL zaxisDefVct(of%cdiZaxisID(ZA_hybrid_half_hhl), 2*nlevp1, vct(1:2*nlevp1))


    ! --------------------------------------------------------------------------------------
    ! Axes for soil model (ZAXIS_DEPTH_BELOW_LAND) -----------------------------------------
    ! --------------------------------------------------------------------------------------

    of%cdiZaxisID(ZA_depth_below_land_p1) = &
      & zaxisCreate(ZAXIS_DEPTH_BELOW_LAND, znlev_soil+1)
    ALLOCATE(levels(znlev_soil+1))
    levels(1) = 0._dp
    DO k = 1, znlev_soil
      levels(k+1) = REAL(zml_soil(k)*1000._wp,dp)  ! in mm
    END DO
    CALL zaxisDefLevels(of%cdiZaxisID(ZA_depth_below_land_p1), levels)
    CALL zaxisDefUnits(of%cdiZaxisID(ZA_depth_below_land_p1), "mm")
    DEALLOCATE(levels)

    of%cdiZaxisID(ZA_depth_below_land) = zaxisCreate(ZAXIS_DEPTH_BELOW_LAND, znlev_soil)
    ALLOCATE(lbounds(znlev_soil), ubounds(znlev_soil), levels(znlev_soil))
    lbounds(1) = 0._dp   ! surface
    DO k = 2, znlev_soil
      lbounds(k)   = REAL((zml_soil(k-1) + (zml_soil(k-1) - lbounds(k-1))),dp)
    ENDDO
    DO k = 1, znlev_soil
      ubounds(k) = REAL((zml_soil(k) + (zml_soil(k) - lbounds(k))),dp)
      levels(k)  = REAL(zml_soil(k)*1000._wp,dp)
    ENDDO
    ubounds(:) = ubounds(:) * 1000._dp        ! in mm
    lbounds(:) = lbounds(:) * 1000._dp        ! in mm

    CALL zaxisDefLbounds(of%cdiZaxisID(ZA_depth_below_land), lbounds) !necessary for GRIB2
    CALL zaxisDefUbounds(of%cdiZaxisID(ZA_depth_below_land), ubounds) !necessary for GRIB2
    CALL zaxisDefLevels (of%cdiZaxisID(ZA_depth_below_land), levels)  !necessary for NetCDF
    CALL zaxisDefUnits  (of%cdiZaxisID(ZA_depth_below_land), "mm")
    DEALLOCATE(lbounds, ubounds, levels)


    ! --------------------------------------------------------------------------------------
    ! Axes for multi-layer snow model (ZAXIS_SNOW) -----------------------------------------
    ! --------------------------------------------------------------------------------------

    ! SNOW-layer axis (for multi-layer snow model)
    of%cdiZaxisID(ZA_snow) = zaxisCreate(ZAXIS_SNOW, nlev_snow)
    ALLOCATE(levels(nlev_snow), lbounds(nlev_snow), ubounds(nlev_snow))
    DO k = 1, nlev_snow
      lbounds(k) = REAL(k,dp)
      levels(k)  = REAL(k,dp)
    ENDDO
    DO k = 1, nlev_snow
      ubounds(k) = REAL(k+1,dp)
    ENDDO

    CALL zaxisDefLbounds(of%cdiZaxisID(ZA_snow), lbounds) !necessary for GRIB2
    CALL zaxisDefUbounds(of%cdiZaxisID(ZA_snow), ubounds) !necessary for GRIB2
    CALL zaxisDefLevels(of%cdiZaxisID(ZA_snow), levels)   !necessary for NetCDF
    DEALLOCATE(levels, lbounds, ubounds)

    of%cdiZaxisID(ZA_snow_half) = zaxisCreate(ZAXIS_SNOW, nlev_snow+1)
    ALLOCATE(levels(nlev_snow+1))
    DO k = 1, nlev_snow+1
      levels(k) = REAL(k,dp)
    END DO
    CALL zaxisDefLevels(of%cdiZaxisID(ZA_snow_half), levels)
    DEALLOCATE(levels)

#endif
    ! #ifndef __NO_ICON_ATMO__

  END SUBROUTINE setup_ml_axes_atmo


  ! --------------------------------------------------------------------------------------
  !> Setup of vertical axis for output module: Pressure levels
  !
  SUBROUTINE setup_pl_axis_atmo(of)
    TYPE(t_output_file), INTENT(INOUT) :: of      !< output file meta-data
#ifndef __NO_ICON_ATMO__
    ! local variables:
    TYPE (t_value_set), POINTER :: levels
    of%cdiZaxisID(:) = CDI_UNDEFID ! not all are set

    ! surface level (required, e.g., for RLON, RLAT)
    CALL define_single_level_axis(of, ZA_surface, ZAXIS_SURFACE)
    ! p-axis
    !
    levels => nh_pzlev_config(of%log_patch_id)%plevels
    CALL define_vertical_axis(of, ZA_pressure, ZAXIS_PRESSURE,      &
      &                       SIZE(levels%values),                  &
      &                       levels = levels%values,               &
      &                       opt_set_vct_as_levels = .TRUE.)
#endif
    ! #ifndef __NO_ICON_ATMO__
  END SUBROUTINE setup_pl_axis_atmo


  ! --------------------------------------------------------------------------------------
  !> Setup of vertical axis for output module: Height levels
  !
  SUBROUTINE setup_hl_axis_atmo(of)
    TYPE(t_output_file), INTENT(INOUT) :: of      !< output file meta-data
#ifndef __NO_ICON_ATMO__
    ! local variables:
    TYPE (t_value_set), POINTER :: levels
    of%cdiZaxisID(:) = CDI_UNDEFID ! not all are set

    ! surface level (required, e.g., for RLON, RLAT)
    CALL define_single_level_axis(of, ZA_surface, ZAXIS_SURFACE)
    ! Altitude above mean sea level
    !
    levels => nh_pzlev_config(of%log_patch_id)%zlevels
    CALL define_vertical_axis(of, ZA_altitude, ZAXIS_ALTITUDE,      &
      &                       SIZE(levels%values),                  &
      &                       levels = levels%values,               &
      &                       opt_set_vct_as_levels = .TRUE.)
#endif
    ! #ifndef __NO_ICON_ATMO__
  END SUBROUTINE setup_hl_axis_atmo


  ! --------------------------------------------------------------------------------------
  !> Setup of vertical axis for output module: Isentropes
  !
  SUBROUTINE setup_il_axis_atmo(of)
    TYPE(t_output_file), INTENT(INOUT) :: of      !< output file meta-data
#ifndef __NO_ICON_ATMO__
    ! local variables:
    TYPE (t_value_set), POINTER :: levels
    of%cdiZaxisID(:) = CDI_UNDEFID ! not all are set

    ! surface level (required, e.g., for RLON, RLAT)
    CALL define_single_level_axis(of, ZA_surface, ZAXIS_SURFACE)
    ! i-axis (isentropes)
    !
    levels => nh_pzlev_config(of%log_patch_id)%ilevels
    CALL define_vertical_axis(of, ZA_isentropic, ZAXIS_ISENTROPIC,  &
      &                       SIZE(levels%values),                  &
      &                       levels = levels%values,               &
      &                       opt_set_vct_as_levels = .TRUE.)
#endif
    ! #ifndef __NO_ICON_ATMO__
  END SUBROUTINE setup_il_axis_atmo


  ! --------------------------------------------------------------------------------------
  !> Setup of vertical axes for output module: Ocean component
  !
  SUBROUTINE setup_zaxes_oce(of)
    TYPE(t_output_file), INTENT(INOUT) :: of  !< output file meta-data
    ! local variables
    REAL(dp), ALLOCATABLE           :: levels(:)
    REAL(wp), ALLOCATABLE           :: levels_i(:), levels_m(:)
    REAL(wp), ALLOCATABLE           :: levels_s(:), levels_sp(:)
    INTEGER                         :: nzlevp1

#ifndef __NO_ICON_OCEAN__
    of%cdiZaxisID(:) = CDI_UNDEFID ! not all are set

    ! surface level
 !  of%cdiZaxisID(ZA_surface) = zaxisCreate(ZAXIS_SURFACE, 1)
 !  ALLOCATE(levels(1))
 !  levels(1) = 0.0_dp
 !  CALL zaxisDefLevels(of%cdiZaxisID(ZA_surface), levels)
 !  DEALLOCATE(levels)

 !  of%cdiZaxisID(ZA_depth_below_sea)      = zaxisCreate(ZAXIS_DEPTH_BELOW_SEA, n_zlev)
 !  of%cdiZaxisID(ZA_depth_below_sea_half) = zaxisCreate(ZAXIS_DEPTH_BELOW_SEA, nzlevp1)
    
    ALLOCATE(levels_i(n_zlev+1))
    ALLOCATE(levels_m(n_zlev))
    CALL set_zlev(levels_i, levels_m, n_zlev, dzlev_m)
!   CALL zaxisDefLevels(of%cdiZaxisID(ZA_DEPTH_BELOW_SEA), REAL(levels_m,dp))
!   CALL zaxisDefLevels(of%cdiZaxisID(ZA_DEPTH_BELOW_SEA_HALF), REAL(levels_i,dp))
!   of%cdiZaxisID(ZA_GENERIC_ICE) = zaxisCreate(ZAXIS_GENERIC, 1)

    CALL define_single_level_axis(of, ZA_surface, ZAXIS_SURFACE)
    CALL define_vertical_axis(of, ZA_depth_below_sea,      ZAXIS_DEPTH_BELOW_SEA, n_zlev  , levels = REAL(levels_m,wp))
    CALL define_vertical_axis(of, ZA_depth_below_sea_half, ZAXIS_DEPTH_BELOW_SEA, n_zlev+1, levels = REAL(levels_i,wp))
    CALL define_single_level_axis(of, ZA_GENERIC_ICE, ZAXIS_GENERIC)

    DEALLOCATE(levels_i)
    DEALLOCATE(levels_m)
    of%cdiZaxisID(ZA_GENERIC_ICE) = zaxisCreate(ZAXIS_GENERIC, 1)
    if(lhamocc)then
    ! ocean sediment
    of%cdiZaxisID(ZA_OCEAN_SEDIMENT) = zaxisCreate(ZAXIS_GENERIC, ks)
    ALLOCATE(levels_s(ks))
    ALLOCATE(levels_sp(ksp))

    CALL set_zlev(levels_sp, levels_s, ks, dzsed*1000._dp)
    CALL zaxisDefLevels(of%cdiZaxisID(ZA_OCEAN_SEDIMENT), REAL(levels_s,dp))
    DEALLOCATE(levels_s)
    DEALLOCATE(levels_sp)
    endif
#endif

  END SUBROUTINE setup_zaxes_oce


  ! --------------------------------------------------------------------------------------
  !> Utility function: defines z-axis with a single level
  !
  SUBROUTINE define_single_level_axis(of, za_type, cdi_type, opt_level_value, opt_unit)
    TYPE(t_output_file), INTENT(INOUT) :: of       !< output file meta-data
    INTEGER,             INTENT(IN)    :: za_type  !< ICON-internal axis ID (see mo_cdi_constants)
    INTEGER,             INTENT(IN)    :: cdi_type !< CDI-internal axis name

    REAL(dp),         INTENT(IN), OPTIONAL :: opt_level_value   !< level value
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: opt_unit          !< axis unit
    ! local variables
    REAL(dp), ALLOCATABLE           :: levels(:)

    of%cdiZaxisID(za_type) = zaxisCreate(cdi_type, 1)
    ALLOCATE(levels(1))
    IF (PRESENT(opt_level_value)) THEN
      levels(1) = opt_level_value
    ELSE
      levels(1) = 0.0_dp
    END IF
    CALL zaxisDefLevels(of%cdiZaxisID(za_type), levels)
    IF (PRESENT(opt_unit)) THEN
      CALL zaxisDefUnits(of%cdiZaxisID(za_type), TRIM(opt_unit))
    END IF
    DEALLOCATE(levels)
  END SUBROUTINE define_single_level_axis


  ! --------------------------------------------------------------------------------------
  !> Utility function: defines z-axis with a single *layer*
  !
  SUBROUTINE define_single_layer_axis(of, za_type, cdi_type, lbound, ubound, unit)
    TYPE(t_output_file), INTENT(INOUT) :: of             !< output file meta-data
    INTEGER,             INTENT(IN)    :: za_type        !< ICON-internal axis ID (see mo_cdi_constants)
    INTEGER,             INTENT(IN)    :: cdi_type       !< CDI-internal axis name
    REAL(dp),            INTENT(IN)    :: lbound, ubound !< lower and upper bound
    CHARACTER(LEN=*),    INTENT(IN)    :: unit           !< axis unit
    ! local variables
    REAL(dp), ALLOCATABLE           :: levels(:), lbounds(:), ubounds(:)

    of%cdiZaxisID(za_type) = zaxisCreate(cdi_type, 1)
    ALLOCATE(lbounds(1), ubounds(1), levels(1))
    lbounds(1) = lbound
    ubounds(1) = ubound
    levels(1)  = lbound
    CALL zaxisDefLbounds(of%cdiZaxisID(za_type), lbounds) !necessary for GRIB2
    CALL zaxisDefUbounds(of%cdiZaxisID(za_type), ubounds) !necessary for GRIB2
    CALL zaxisDefLevels(of%cdiZaxisID(za_type), levels)
    CALL zaxisDefUnits(of%cdiZaxisID(za_type), TRIM(unit))
    DEALLOCATE(lbounds, ubounds, levels)
  END SUBROUTINE define_single_layer_axis


  ! --------------------------------------------------------------------------------------
  !> Utility function: defines z-axis with given levels, lower and upper bounds.
  !
  SUBROUTINE define_vertical_axis(of, za_type, cdi_type, in_nlevs, levels,     &
    &                             opt_set_bounds, opt_set_ubounds_value,       &
    &                             opt_name, opt_number, opt_uuid, opt_set_vct_as_levels)
    TYPE(t_output_file),     INTENT(INOUT) :: of             !< output file meta-data
    INTEGER,                 INTENT(IN)    :: za_type        !< ICON-internal axis ID (see mo_cdi_constants)
    INTEGER,                 INTENT(IN)    :: cdi_type       !< CDI-internal axis name 
    INTEGER,                 INTENT(IN)    :: in_nlevs       !< no. of levels (if no selection)
    REAL(dp),                INTENT(IN)    :: levels(:)      !< axis levels
    LOGICAL,                 INTENT(IN), OPTIONAL :: opt_set_bounds        !< Flag. Set lower/upper bounds if .TRUE.
    REAL(dp),                INTENT(IN), OPTIONAL :: opt_set_ubounds_value !< Explicit value for ubounds
    CHARACTER(*),            INTENT(IN), OPTIONAL :: opt_name              !< Name of the zaxis    
    INTEGER,                 INTENT(IN), OPTIONAL :: opt_number            !< numberOfVGridUsed
    INTEGER(KIND = C_SIGNED_CHAR), INTENT(IN), OPTIONAL :: opt_uuid(16)          !< UUID of vertical grid
    LOGICAL,                 INTENT(IN), OPTIONAL :: opt_set_vct_as_levels !< set VCT to level values
    ! local variables
    CHARACTER(*), PARAMETER :: routine = TRIM(modname)//":define_vertical_axis"
    REAL(dp), ALLOCATABLE :: axis_levels(:), ubounds(:)     !< level list, upper bounds

    INTEGER :: nlev                                         !< no. of axis levels
    LOGICAL :: set_bounds, select_levs
    INTEGER :: ierrstat, i, isrc_lev

    ! First, build a level list. If no "of%level_selection" was
    ! provided, this level lists is simply a copy of the input
    ! parameter "levels(:)":
    select_levs = ASSOCIATED(of%level_selection)
    IF (.NOT. select_levs) THEN
      nlev = in_nlevs
    ELSE
      ! count the number of feasible levels that have been selected:
      nlev = 0
      DO i=1,of%level_selection%n_selected
        IF (of%level_selection%global_idx(i) <= SIZE(levels)) THEN
          nlev = nlev + 1
        END IF
      END DO
    END IF
    ALLOCATE(axis_levels(nlev),STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')
    IF (.NOT. select_levs) THEN
      axis_levels(1:nlev) = levels(1:nlev)
    ELSE
      axis_levels(1:nlev) = levels(of%level_selection%global_idx(1:nlev))
    END IF

    ! create vertical axis object
    of%cdiZaxisID(za_type) = zaxisCreate(cdi_type, nlev)
    CALL zaxisDefLevels(of%cdiZaxisID(za_type), axis_levels ) !necessary for NetCDF

    IF (PRESENT(opt_name)) CALL zaxisDefName(of%cdiZaxisID(za_type), TRIM(opt_name))

    set_bounds = .FALSE.
    IF (PRESENT(opt_set_bounds)) set_bounds = opt_set_bounds

    IF (set_bounds) THEN
      ! note: the axes in ICON use the "axis_levels" values as lower bounds
      CALL zaxisDefLbounds(of%cdiZaxisID(za_type), axis_levels) !necessary for GRIB2

      ALLOCATE(ubounds(nlev))
      IF (PRESENT(opt_set_ubounds_value)) THEN
        ubounds(:) = opt_set_ubounds_value
      ELSE
        ubounds(:) = 0._dp
        IF (.NOT. select_levs) THEN
          ubounds(1:nlev) = levels(2:(nlev+1))
        ELSE
          DO i=1,nlev
            isrc_lev = of%level_selection%global_idx(i)+1
            IF (isrc_lev <= SIZE(levels)) THEN
              ubounds(i) = levels(isrc_lev)  
            END IF
          END DO
        END IF
      END IF
      CALL zaxisDefUbounds(of%cdiZaxisID(za_type), ubounds) !necessary for GRIB2
      DEALLOCATE(ubounds)
    END IF

    ! Set numberOfVGridUsed
    IF (PRESENT(opt_number)) THEN
      ! Dependent on the algorithm chosen to generate the vertical grid (ivctype)
      CALL zaxisDefNumber(of%cdiZaxisID(za_type), opt_number )
    END IF
    ! Write vertical grid UUID
    IF (PRESENT(opt_uuid)) THEN
      CALL zaxisDefUUID (of%cdiZaxisID(za_type), opt_uuid ) !uuidOfVGrid
    END IF
    ! Set Vertical Coordinate Table (VCT) to level values
    IF (PRESENT(opt_set_vct_as_levels)) THEN
      IF (opt_set_vct_as_levels)  CALL zaxisDefVct(of%cdiZaxisID(za_type), nlev, axis_levels)
    END IF

    DEALLOCATE(axis_levels, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed.')
  END SUBROUTINE define_vertical_axis


  ! --------------------------------------------------------------------------------------
  !> FUNCTION get_numberOfVGridUsed
  !  Depending on the vertical axis chosen for ICON (ivctype), it gives back the value for
  !  the GRIB2-key 'numberOVGridUsed'. Here, we adhere to the COSMO implementation:
  !
  !  |       Description               |  ivctype  |  numberOfVGridUsed  |
  !  =====================================================================
  !  | height based hybrid Gal-Chen    |    1      |       2             |
  !  | height based SLEVE              |    2      |       4             |
  !
  !
  FUNCTION get_numberOfVgridUsed(ivctype)
    INTEGER                 :: get_numberOfVgridUsed
    INTEGER, INTENT(IN)     :: ivctype
    CHARACTER(*), PARAMETER :: routine = TRIM(modname)//":get_numberOfVgridUsed"

    get_numberOfVgridUsed = -1
    SELECT CASE(ivctype)
      CASE(1)
        get_numberOfVgridUsed = 2
      CASE(2)
        get_numberOfVgridUsed = 4
      CASE DEFAULT
        CALL finish(routine, "invalid ivctype! Must be 1 or 2")
    END SELECT
  END FUNCTION get_numberOfVgridUsed


  ! --------------------------------------------------------------------------------------
  !> Creates a level selection data object from a selection that is
  !  described as a character string, e.g. "1,5...10,15". See the
  !  module "mo_util_string_parse" for a detailed description of valid
  !  selection strings.
  ! 
  SUBROUTINE create_level_selection_str(selection_str, nlevs, of, opt_nlev_value)
    CHARACTER(LEN=*),        INTENT(IN)    :: selection_str    !< selection described as string
    INTEGER,                 INTENT(IN)    :: nlevs            !< total no. of levels
    TYPE(t_output_file),     INTENT(INOUT) :: of               !< output file meta-data containing level selection
    INTEGER, INTENT(IN), OPTIONAL :: opt_nlev_value            !< number to substitute for "N"/"nlev"
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = '::create_level_selection_str'
    INTEGER :: int_list(0:nlevs) ! note: 0-lower bound required
    INTEGER :: ierrstat, nlev_value

    IF (TRIM(selection_str) == "") RETURN ! do nothing

    ALLOCATE(of%level_selection)
    ALLOCATE(of%level_selection%s(nlevs), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')
    ! parse character string into a LOGICAL array:
    nlev_value = nlevs
    IF (PRESENT(opt_nlev_value)) nlev_value = opt_nlev_value
    CALL util_do_parse_intlist(selection_str, nlev_value, int_list, ierrstat)
    of%level_selection%s(1:nlevs) = (int_list(1:nlevs) == 1)
    IF (ierrstat /= 0) CALL finish(routine, 'Parsing of level selection failed.')
    ! count no. of selected levels
    of%level_selection%n_selected = COUNT(of%level_selection%s)
    ! get mapping: global level -> local level index
    ALLOCATE(of%level_selection%global_idx(of%level_selection%n_selected), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')
    CALL get_level_indices_in_selection(of%level_selection%s, of%level_selection%n_selected, &
      &                                 of%level_selection%global_idx)
    ! get mapping: local level -> global level index
    ALLOCATE(of%level_selection%local_idx(nlevs), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')
    CALL get_selection_index(of%level_selection%s, of%level_selection%local_idx)
  END SUBROUTINE create_level_selection_str


  ! --------------------------------------------------------------------------------------
  !  > Creates a level selection data object from a selection that is
  !    described as list of REAL(wp) that partly match the levels
  !    contained in a "t_value_set" object.
  !
  !  @note A list of level values starting with a negative value means
  !        that no level selection is to be created (this is the
  !        default).
  ! 
  SUBROUTINE create_level_selection_set(selection, value_set, of)
    REAL(wp),                INTENT(IN)    :: selection(:)     !< selected levels
    TYPE (t_value_set),      INTENT(IN)    :: value_set        !< total set of available levels
    TYPE(t_output_file),     INTENT(INOUT) :: of               !< output file meta-data containing level selection
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = '::create_level_selection_set'
    INTEGER :: ierrstat, nlevs

    IF (selection(1) < 0._wp)  RETURN ! do nothing
    ! count the no. of levels
    DO nlevs=1,SIZE(selection)
      IF (selection(nlevs) < 0._wp) EXIT
    END DO
    nlevs = nlevs - 1
    ! allocate the level selection object
    ALLOCATE(of%level_selection)
    ALLOCATE(of%level_selection%s(value_set%nvalues), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')
    ! determine the selected levels
    CALL find_values_in_set(nlevs, selection, value_set, of%level_selection%s)

    ! count no. of selected levels
    of%level_selection%n_selected = COUNT(of%level_selection%s)
    ! get mapping: global level -> local level index
    ALLOCATE(of%level_selection%global_idx(of%level_selection%n_selected), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')
    CALL get_level_indices_in_selection(of%level_selection%s, of%level_selection%n_selected, &
      &                                 of%level_selection%global_idx)
    ! get mapping: local level -> global level index
    ALLOCATE(of%level_selection%local_idx(value_set%nvalues), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')
    CALL get_selection_index(of%level_selection%s, of%level_selection%local_idx)
  END SUBROUTINE create_level_selection_set


  ! --------------------------------------------------------------------------------------
  !> Frees a level selection data object.
  !
  SUBROUTINE deallocate_level_selection(selection)
    TYPE(t_level_selection), INTENT(INOUT) :: selection !< level selection object
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = '::deallocate_level_selection'
    INTEGER :: ierrstat
    DEALLOCATE(selection%s, selection%global_idx, selection%local_idx, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed.')
  END SUBROUTINE deallocate_level_selection


  ! --------------------------------------------------------------------------------------
  !> Utility function: Assume that we have N vertical levels, out of
  !  which only a few levels are selected. This subroutine takes this
  !  selection as input in the form of a LOGICAL array s(1...N), where
  !  "s(i)=.TRUE." means that level "i" is selected. As an output, we
  !  get an integer list idx(1...n_selected) containing the selected
  !  level indices.
  !
  SUBROUTINE get_level_indices_in_selection(s, n_selected, idx)
    LOGICAL, INTENT(IN)  :: s(:)       !< level selection (LOGICAL array)
    INTEGER, INTENT(OUT) :: n_selected !< no. of selected levels
    INTEGER, INTENT(OUT) :: idx(:)     !< selected level indices
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = '::get_level_indices_in_selection'
    INTEGER :: nlevs, i

    nlevs = SIZE(s)
    n_selected = 0
    DO i=1,nlevs
      IF (s(i)) THEN
        n_selected = n_selected + 1
        IF (SIZE(idx) < n_selected)  CALL finish(routine, "Dimension mismatch!")
        idx(n_selected) = i
      END IF
    END DO
  END SUBROUTINE get_level_indices_in_selection


  ! --------------------------------------------------------------------------------------
  !> Utility function: Assume that we have N vertical levels, out of
  !  which only a few levels are selected. This subroutine takes this
  !  selection as input in the form of a LOGICAL array s(1...N), where
  !  "s(i)=.TRUE." means that level "i" is selected. As an output, we
  !  get an integer list idx(1...N) containing the local index in the
  !  list of selected level indices (i.e. an integer number in the
  !  range 1...n_selected).
  !
  SUBROUTINE get_selection_index(s, idx)
    LOGICAL, INTENT(IN)  :: s(:)       !< level selection (LOGICAL array)
    INTEGER, INTENT(OUT) :: idx(:)     !< selected level indices
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = '::get_selection_index'
    INTEGER :: nlevs, i, n_selected

    nlevs = SIZE(s)
    n_selected = 0
    DO i=1,nlevs
      IF (s(i)) THEN
        n_selected = n_selected + 1
        idx(i)     = n_selected
      ELSE
        idx(i)     = 0
      END IF
    END DO
  END SUBROUTINE get_selection_index


  !------------------------------------------------------------------------------------------------
  ! Loop over output file (p_of) and create the "selection" from the
  ! union set of vertical model/pressure/isentropic/height (MIPZ)
  ! levels:
  SUBROUTINE create_mipz_level_selections(output_file)
    TYPE(t_output_file), TARGET, INTENT(INOUT) :: output_file(:)
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = '::create_mipz_level_selections'
    TYPE (t_output_file), POINTER   :: p_of
    INTEGER :: i

#ifndef __NO_ICON_ATMO__
    DO i=1,SIZE(output_file)
      p_of => output_file(i)
      SELECT CASE(p_of%ilev_type)
      CASE (level_type_ml) 
        ! note: "ml" comprises both full and half levels
        CALL create_level_selection(p_of%name_list%m_levels,              &
          &                         num_lev(p_of%log_patch_id)+1, p_of,   &
          &                         opt_nlev_value = num_lev(p_of%log_patch_id))
      CASE (level_type_pl) 
        CALL create_level_selection(p_of%name_list%p_levels, &
          &                         nh_pzlev_config(p_of%log_patch_id)%plevels, p_of)
      CASE (level_type_hl) 
        CALL create_level_selection(p_of%name_list%z_levels, &
          &                         nh_pzlev_config(p_of%log_patch_id)%zlevels, p_of)
      CASE (level_type_il) 
        CALL create_level_selection(p_of%name_list%i_levels, &
          &                         nh_pzlev_config(p_of%log_patch_id)%ilevels, p_of)
      CASE DEFAULT
        CALL finish(routine, "Internal error!")
      END SELECT
    END DO
#endif

  END SUBROUTINE create_mipz_level_selections

END MODULE mo_name_list_output_zaxes

