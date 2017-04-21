!>
!! @brief configuration setup for NWP land scheme TERRA
!!
!! configuration setup for NWP land scheme TERRA
!! <Details of procedures are documented below with their definitions.>
!! <Include any applicable external references inline as module::procedure,>
!! <external_procedure(), or by using @see.>
!! <Don't forget references to literature.>
!!
!! @author Kristina Froehlich, MPI-M (2011-07-14)
!! @author <name, affiliation>
!!
!!
!! @par Revision History
!! <Description of activity> by <name, affiliation> (<YYYY-MM-DD>)
!!
!! Modifications by Dmitrii Mironov, DWD (2016-08-04)
!! - Logical switch is introduced to allow the use of 
!!   a rate equation for the sea-ice albedo.
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_lnd_nwp_config

  USE mo_kind,               ONLY: wp
  USE mo_impl_constants,     ONLY: zml_soil, dzsoil, GLOBCOVER2009, GLC2000
  USE mo_var_metadata_types, ONLY: CLASS_TILE, CLASS_TILE_LAND
  USE mo_io_units,           ONLY: filename_max
  USE mo_exception,          ONLY: finish
  USE mo_util_cdi,           ONLY: trivial_tileinfo, trivial_tileId
  USE mo_util_string,        ONLY: int2string

  IMPLICIT NONE

  PRIVATE

  ! FUNCTIONS/SUBROUTINES
  PUBLIC :: getNumberOfTiles
  PUBLIC :: select_tile
  PUBLIC :: configure_lnd_nwp
  PUBLIC :: convert_luc_ICON2GRIB
  PUBLIC :: get_tile_suffix
  PUBLIC :: find_tile_id

  ! TYPES
  PUBLIC :: t_GRIB2_tile
  PUBLIC :: t_GRIB2_att
  PUBLIC :: t_tile

  ! VARIABLES
  PUBLIC :: nlev_soil, nlev_snow, ibot_w_so, ntiles_total, ntiles_lnd, ntiles_water
  PUBLIC :: frlnd_thrhld, frlndtile_thrhld, frlake_thrhld, frsea_thrhld
  PUBLIC :: lseaice, lprog_albsi, llake, lmelt, lmelt_var, lmulti_snow, lsnowtile, max_toplaydepth
  PUBLIC :: itype_trvg, itype_evsl, itype_lndtbl, l2lay_rho_snow
  PUBLIC :: itype_root, itype_heatcond, itype_interception, &
             itype_hydbound, idiag_snowfrac, cwimax_ml, c_soil, c_soil_urb
  PUBLIC :: lstomata, l2tls, lana_rho_snow 
  PUBLIC :: isub_water, isub_lake, isub_seaice
  PUBLIC :: sstice_mode, sst_td_filename, ci_td_filename
  PUBLIC :: tiles


  !--------------------------------------------------------------------------
  ! Basic configuration setup for NWP land
  !--------------------------------------------------------------------------

!  TYPE t_nwp_lnd_config

  ! namelist variables
  INTEGER ::  nlev_snow          !< number of snow layers
  INTEGER ::  ntiles_lnd         !< number of static land surface types
  REAL(wp)::  frlnd_thrhld       !< fraction threshold for creating a land grid point
  REAL(wp)::  frlndtile_thrhld   !< fraction threshold for retaining the respective 
                                 !< tile for a grid point
  REAL(wp)::  frlake_thrhld      !< fraction threshold for creating a lake grid point
  REAL(wp)::  frsea_thrhld       !< fraction threshold for creating a sea grid point
  INTEGER ::  itype_trvg         !< type of vegetation transpiration parameterization
  INTEGER ::  itype_evsl         !< type of parameterization of bare soil evaporation
  INTEGER ::  itype_lndtbl       !< choice of table for associating surface parameters to land-cover classes
  INTEGER ::  itype_root         !< type of root density distribution
  INTEGER ::  itype_heatcond     !< type of soil heat conductivity
  INTEGER ::  itype_interception !< type of plant interception
  REAL(wp)::  cwimax_ml          !< scaling parameter for maximum interception storage
  REAL(wp)::  c_soil             !< surface area density of the (evaporative) soil surface
  REAL(wp)::  c_soil_urb         !< surface area density of the (evaporative) soil surface, urban areas
  INTEGER ::  itype_hydbound     !< type of hydraulic lower boundary condition
  INTEGER ::  idiag_snowfrac     !< method for diagnosis of snow-cover fraction

  LOGICAL ::  lseaice     !> forecast with sea-ice model
  LOGICAL ::  lprog_albsi !> sea-ice albedo is computed prognostically from a rate equation
  LOGICAL ::  llake       !! forecast with lake model FLake
  LOGICAL ::  lmelt       !! soil model with melting process
  LOGICAL ::  lmelt_var   !! freezing temperature dependent on water content
  LOGICAL ::  lmulti_snow !! run the multi-layer snow model
  LOGICAL ::  l2lay_rho_snow !! use two-layer snow density for single-layer snow model
  REAL(wp)::  max_toplaydepth !< maximum depth of uppermost snow layer for multi-layer snow scheme
  LOGICAL ::  lstomata    !! map of minimum stomata resistance
  LOGICAL ::  l2tls       !! forecast with 2-TL integration scheme
  LOGICAL ::  lana_rho_snow !! if .TRUE., take rho_snow-values from analysis file 
  LOGICAL ::  lsnowtile   !! if .TRUE., snow is considered as a separate tile
 
  INTEGER ::  sstice_mode      !< set if SST and sea ice cover are read from the analysis
                               !< and kept constant or read from external data files 
                               !< and updated regularly in run time
  CHARACTER(LEN=filename_max) :: sst_td_filename, ci_td_filename

  ! derived variables
  INTEGER ::  nlev_soil    !< number of soil layers (based on zml_soil in impl_constants)
  INTEGER ::  ibot_w_so    !< number of hydrological active soil layers 
  INTEGER ::  isub_water   !< (open) water points tile number
  INTEGER ::  isub_lake    !< lake points tile number
  INTEGER ::  isub_seaice  !< seaice tile number
  INTEGER ::  ntiles_total !< total number of TILES
  INTEGER ::  ntiles_water !< number of extra tiles for ocean, seaice and lakes

!  END TYPE t_nwp_lnd_config

  !>
  !!
!  TYPE(t_nwp_lnd_config) :: nwp_lnd_config(max_dom)


   TYPE t_GRIB2_tile
     INTEGER :: tileIndex                    ! tile Index (1,...,numberOfTiles)
     INTEGER :: numberOfTileAttributes       ! number of used tile attributes
   END TYPE t_GRIB2_tile

   TYPE t_GRIB2_att
     INTEGER     :: tileAttribute
   END TYPE t_GRIB2_att

   TYPE t_tile
     TYPE(t_GRIB2_tile) :: GRIB2_tile
     TYPE(t_GRIB2_att)  :: GRIB2_att
   END TYPE t_tile

   TYPE(t_tile), ALLOCATABLE :: tiles(:)

   CHARACTER(LEN = *), PARAMETER :: modname = "mo_lnd_nwp_config"

CONTAINS

  !>
  !! setup components of the NWP land scheme depending on its namelist
  !!
  !! Setup of additional nwp-land control variables depending on the 
  !! land-NAMELIST and potentially other namelists. This routine is 
  !! called, after all namelists have been read and a synoptic consistency 
  !! check has been done.
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2011-08-01)
  !!
  SUBROUTINE configure_lnd_nwp()
  !
    ! local variables
    INTEGER  :: kso              ! soil loop index
    REAL(wp) :: depth_hl         ! half level depth

    CHARACTER(len=*), PARAMETER::  &
      &  routine = 'mo_lnd_nwp_config: configure_lnd_nwp'
    !-----------------------------------------------------------------------

    ! number of soil layers
    ! Note that this number must be consistent with the number of entries 
    ! in zml_soil. zml_soil provides soil layer full level heights (confirmed).
    nlev_soil = SIZE(zml_soil)


    ! number of hydraulical active soil layers ibot_w_so
    !
    ! currently, we take all those layers which are located completely 
    ! above 2.5m soil depth
    DO kso=1,nlev_soil
      depth_hl = zml_soil(kso) + 0.5_wp*dzsoil(kso)
      IF (depth_hl<=2.5_wp) ibot_w_so=kso
    ENDDO
    ! make sure that ibot_w_so>=2
    ibot_w_so = MAX(2, ibot_w_so)


    IF (ntiles_lnd == 1) THEN ! Reset options that can be used in combination with tile approach
      lsnowtile     = .FALSE.
      frlnd_thrhld  = 0.5_wp
      frlake_thrhld = 0.5_wp
      frsea_thrhld  = 0.5_wp
    ENDIF

    ! number of tiles; ntiles_lnd is set by the namelist variable ntiles
    IF(lsnowtile) THEN
      ntiles_total = 2*ntiles_lnd
    ELSE
      ntiles_total = ntiles_lnd
    END IF

    IF (ntiles_total == 1) THEN 
      ! no tile approach, thus no extra tile index for water points
      ntiles_water = 0
    ELSE
      ! extra tiles for water points
      ntiles_water = 3 ! one for open sea points
                       ! another one for lake points
                       ! another one for seaice
    ENDIF


    !
    ! ASCII art for internal tile structure, showing the configuration for:
    ! - ntiles_lnd = 3
    ! - lsnowtile  = .TRUE.
    !
    !------------------------------------------------------------------!
    !   1  |   2  |   3  ||   4  |   5  |   6  ||   7  |   8  |    9   !
    ! land | land | land || snow | snow | snow || open | lake | seaice !
    !      |      |      ||      |      |      || sea  |      |        !
    !------------------------------------------------------------------!
    !
    !<---  ntiles_lnd --->                      <---- ntiles_water ---->
    !
    !<-------------- ntiles_total ------------>
    !
    !<---------------- ntiles_total + ntiles_water ------------------->
    !

    ! ASCII art for corresponding GRIB2 tile structure:
    !
    !-----------------------------------------------------------------!
    !    land    |    land    |    land    |     oce    |    lake     !
    !      1     |      2     |      3     |      4     |     5       !  <- GRIB2 Tile ID
    !     / \    |     / \    |     / \    |     / \    |     |       !
    !****/***\***|****/***\***|****/***\***|****/***\***|*****|*******!  *******************
    ! umod | snw | umod | snw | umod | snw | umod | ice |   undef     !  <- Attribute
    !  1   |  4  |   2  |  5  |   3  |  6  |   7  |  9  |     8       !  <- Internal Tile ID
    !-----------------------------------------------------------------!     


    ! (open) water points tile number
    isub_water  = MAX(1,ntiles_total + ntiles_water - 2)

    ! lake points tile number
    isub_lake   = MAX(1,ntiles_total + ntiles_water - 1)

    ! sea-ice tile number
    isub_seaice = ntiles_total + ntiles_water


    ! Setup tile meta-information required for tile I/O.
    CALL setup_tile_metainfo (tiles, ntiles_lnd, ntiles_total, ntiles_water)

  END SUBROUTINE configure_lnd_nwp


  !>
  !! Provides number of tiles
  !!
  !! Provides number of tiles, as it is required for GRIB2 encoding.
  !! It differs from the ICON internal counting rules for tiles, in that 
  !! snowtiles and the sea-ice tile are not considered as separate tiles.
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2015-01-22)
  !!
  FUNCTION getNumberOfTiles (class)  RESULT (numberOfTiles)

    INTEGER, INTENT(IN) :: class           ! CLASS_TILE
                                           !   variable contains land and water tiles 
                                           ! CLASS_TILE_LAND
                                           !   variable contains only land tiles
    INTEGER             :: numberOfTiles   ! number of tiles for GRIB2 encoding

    SELECT CASE(class)
    CASE(CLASS_TILE)
      ! snow tiles and sea-ice tile are not taken into account.
      numberOfTiles = MAX(1,ntiles_lnd + ntiles_water - 1)
    CASE(CLASS_TILE_LAND)
      ! water tiles, snow tiles and sea-ice tile are not taken into account.
      numberOfTiles = ntiles_lnd
    CASE DEFAULT
      CALL finish( 'mo_lnd_nwp_config:getNumberOfTiles', 'class not supported' )

    END SELECT

  END FUNCTION getNumberOfTiles


  !>
  !! Setup tile meta-information
  !!
  !! Setup tile meta-information required for tile I/O. 
  !! For each ICON-tile the following keys are defined:
  !! - tileIndex
  !! - numberOfTileAttributes
  !! - tileAttribute
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2015-01-23)
  !!
  SUBROUTINE setup_tile_metainfo (tiles, ntiles_lnd, ntiles_total, ntiles_water)

    TYPE(t_tile), ALLOCATABLE, INTENT(INOUT) :: tiles(:)
    INTEGER                  , INTENT(IN)    :: ntiles_lnd 
    INTEGER                  , INTENT(IN)    :: ntiles_total
    INTEGER                  , INTENT(IN)    :: ntiles_water

    ! Local
    INTEGER :: i                 ! tile loop counter
    INTEGER :: nat_lnd, nat_oce  ! number of attributes for land tiles and ocean tile
    INTEGER :: istat             ! error flag

    ! tile attributes
    INTEGER, PARAMETER :: UNDEF = 0  ! undefined
    INTEGER, PARAMETER :: UNMOD = 1  ! unmodified
    INTEGER, PARAMETER :: SNOW  = 2  ! snow-covered
    INTEGER, PARAMETER :: SEAICE= 4  ! sea ice covered

    !-----------------------------------------------------------------------

    ! allocate tile array
    ALLOCATE(tiles(ntiles_total + ntiles_water), STAT=istat)
    IF (istat /= 0) THEN
      CALL finish('mo_lnd_nwp_config:setup_tile_metainfo', &
        &      'allocation of array tiles failed')
    ENDIF


    IF (ntiles_total > 1) THEN  ! tile approach is used

      ! define number of attributes for land tiles
      ! can be unmodified or snow covered
      IF (lsnowtile) THEN
        nat_lnd = 2
      ELSE
        nat_lnd = 1
      ENDIF

      ! define number of attributes for ocean tile
      ! can be unmodified or sea ice covered
      IF (ntiles_water > 0) THEN
        nat_oce = 2
      ELSE
        nat_oce = 0
      ENDIF


      ! fill in tile meta information
      DOTILES: DO i=1, SIZE(tiles)
 
        IF (i<=ntiles_lnd) THEN  ! dominating land tiles
  
          tiles(i)%GRIB2_tile%tileIndex              = i
          tiles(i)%GRIB2_tile%numberOfTileAttributes = nat_lnd
          tiles(i)%GRIB2_att%tileAttribute           = MERGE(UNMOD, UNDEF, nat_lnd>1)
  
        ELSE IF ( (i>ntiles_lnd) .AND. (i<= ntiles_total) ) THEN ! corresponding snow tiles (if present)
  
          tiles(i)%GRIB2_tile%tileIndex              = i - ntiles_lnd
          tiles(i)%GRIB2_tile%numberOfTileAttributes = nat_lnd
          tiles(i)%GRIB2_att%tileAttribute           = SNOW
  
        ELSE IF ( i == isub_water ) THEN ! ocean tiles (unfrozen)
  
          tiles(i)%GRIB2_tile%tileIndex              = ntiles_lnd + 1
          tiles(i)%GRIB2_tile%numberOfTileAttributes = nat_oce
          tiles(i)%GRIB2_att%tileAttribute           = UNMOD
  
        ELSE IF ( i == isub_lake ) THEN  ! lake tile
  
          tiles(i)%GRIB2_tile%tileIndex              = ntiles_lnd + 2
          tiles(i)%GRIB2_tile%numberOfTileAttributes = 1
          tiles(i)%GRIB2_att%tileAttribute           = UNDEF
  
        ELSE IF ( i == isub_seaice ) THEN ! sea-ice tile
  
          tiles(i)%GRIB2_tile%tileIndex              = ntiles_lnd + 1
          tiles(i)%GRIB2_tile%numberOfTileAttributes = nat_oce
          tiles(i)%GRIB2_att%tileAttribute           = SEAICE

        ELSE
          CALL finish('mo_lnd_nwp_config:setup_tile_metainfo', &
            &      ' failed')
        END IF
  
      ENDDO DOTILES

    ELSE IF (ntiles_total == 1) THEN  ! no tile apporach
      ! set trivial tile info
      tiles(1)%GRIB2_tile%tileIndex                = trivial_tileinfo%idx
      tiles(1)%GRIB2_tile%numberOfTileAttributes   = 0
      tiles(1)%GRIB2_att%tileAttribute             = trivial_tileinfo%att
    ELSE
      CALL finish('mo_lnd_nwp_config:setup_tile_metainfo', &
        &    ' invalid number of tiles ntiles_total')
    ENDIF  ! ntiles_total > 1

  END SUBROUTINE setup_tile_metainfo


  !>
  !! Select tile
  !!
  !! Select specific tile and corresponding tile information.
  !! Since the ICON internal tile nomenclature and structure differes 
  !! from the GRIB2 tile nomenclature, this function provides the GRIB2 
  !! tileIndex, given the internal tile ID (tile number) as input.
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2015-01-23)
  !!
  FUNCTION select_tile (tileID_int)  RESULT (selected_tile)

    INTEGER, INTENT(IN) :: tileID_int

    TYPE(t_tile) :: selected_tile
    CHARACTER(LEN = *), PARAMETER :: routine = modname//":select_tile"

    !-----------------------------------------------------------------------

    IF(tileID_int == trivial_tileId) then
        selected_tile%GRIB2_tile%tileIndex = trivial_tileinfo%idx
        selected_tile%GRIB2_tile%numberOfTileAttributes = 0
        selected_tile%GRIB2_att%tileAttribute = trivial_tileinfo%att
    ELSE IF(tileID_int < 1 .OR. tileID_int > SIZE(tiles, 1)) THEN
        CALL finish(routine, "illegal tileID "//int2string(tileID_int))
    ELSE
        selected_tile =  tiles(tileID_int)
    END IF

  END FUNCTION select_tile


  !>
  !! Convert tile/attribute pair into ICON-internal varname suffix
  !!
  !! Convert tile/attribute pair into ICON-internal varname suffix '_t_X'
  !! Returns the empty string for the trivial_tileinfo.
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2015-03-05)
  !!
  FUNCTION get_tile_suffix(tileIdx, tileAtt) RESULT (tileSuffix)
    INTEGER, VALUE :: tileIdx        ! GRIB2 Tile ID
    INTEGER, VALUE :: tileAtt        ! GRIB2 Tile Attribute

    INTEGER :: tileId
    CHARACTER(LEN=2) :: tileIdx_str
    CHARACTER(LEN=5) :: tileSuffix

    !special case: trivial tileinfo => empty suffix
    IF(tileIdx == trivial_tileinfo%idx .OR. tileAtt == trivial_tileinfo%att) THEN
        tileSuffix =''
        RETURN
    END IF

    !get the tile ID
    tileId = find_tile_id(tileIdx, tileAtt)
    IF(tileId == 0) THEN
        CALL finish('mo_lnd_nwp_config:get_tile_suffix', &
        &           'Invalid tile/attribute pair ('//int2string(tileIdx)//', '//int2string(tileAtt)//')' )
    END IF

    !convert the tile ID to the corresponding suffix string
    WRITE(tileIdx_str,'(i2)') tileId
    tileSuffix = '_t_'//TRIM(ADJUSTL(tileIdx_str))
  END FUNCTION get_tile_suffix

  !>
  !! Find a given tile index-attribute pair within the tiles array AND RETURN its index.
  !! This results IN a unique, ICON-internal ID for the tile.
  !!
  !! Returns zero on failure to find the given tile.
  !!
  FUNCTION find_tile_id(tileIndex, tileAttribute) RESULT(resultVar)
    INTEGER, VALUE :: tileIndex, tileAttribute  !GRIB2 tile index AND attribute values
    INTEGER :: resultVar

    DO resultVar = 1, SIZE(tiles, 1)
        IF(tiles(resultVar)%GRIB2_tile%tileIndex == tileIndex) THEN
            IF(tiles(resultVar)%GRIB2_att%tileAttribute == tileAttribute) THEN
                RETURN
            END IF
        END IF
    END DO
    resultVar = -1
  END FUNCTION find_tile_id


  !>
  !! Given the internal land use class index, provide the official GRIB2 index.
  !!
  !! Given the ICON-internal land use class index, provide the official GRIB2 
  !! index according to GRIB2 table 4.2.43 (or vice versa).
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2015-01-26)
  !!
  FUNCTION convert_luc_ICON2GRIB(lc_datbase,iluc_in,opt_linverse)  RESULT (iluc_out)

    INTEGER          , INTENT(IN) :: lc_datbase       !< landuse class database
    INTEGER          , INTENT(IN) :: iluc_in          !< landuse class for grid point i
    LOGICAL, OPTIONAL, INTENT(IN) :: opt_linverse     !< given the GRIB2 landsuse class, 
                                                      !< provide the internal one

    ! Local
    LOGICAL :: linverse
    INTEGER :: iluc_out          !< landuse class index output
                                 !< Default: according to GRIB2 table 4.243

    INTEGER :: iluc_GLOBCOVER2009(23)  ! num_lcc
    INTEGER :: iluc_GLC2000(23)        ! num_lcc
    !
    ! maps GLOBCOVER2009 landuse classes onto tile class table 4.243
    DATA iluc_GLOBCOVER2009 / 24, &  ! 1 :irrigated croplands
                            & 25, &  ! 2 :rainfed croplands 
                            & 26, &  ! 3 :mosaic cropland (50-70%) - vegetation (20-50%)
                            & 27, &  ! 4 :mosaic vegetation (50-70%) - cropland (20-50%)
                            & 28, &  ! 5 :closed broadleaved evergreen forest           
                            & 2 , &  ! 6 :closed broadleaved deciduous forest           
                            & 3 , &  ! 7 :open broadleaved deciduous forest             
                            & 29, &  ! 8 :closed needleleaved evergreen forest          
                            & 30, &  ! 9 :open needleleaved deciduous forest            
                            & 31, &  ! 10:mixed broadleaved and needleleaved forest     
                            & 32, &  ! 11:mosaic shrubland (50-70%) - grassland (20-50%)
                            & 33, &  ! 12:mosaic grassland (50-70%) - shrubland (20-50%)
                            & 34, &  ! 13:closed to open shrubland                      
                            & 25, &  ! 14:closed to open herbaceous vegetation          
                            & 35, &  ! 15:sparse vegetation                             
                            & 36, &  ! 16:closed to open forest regularly flooded        
                            & 37, &  ! 17:closed forest or shrubland permanently flooded
                            & 38, &  ! 18:closed to open grassland regularly flooded    
                            & 22, &  ! 19:artificial surfaces                           
                            & 19, &  ! 20:bare areas                                     
                            & 20, &  ! 21:water bodies                                   
                            & 21, &  ! 22:permanent snow and ice                         
                            & 39  /  ! 23:undefined
    !
    ! maps GLC2000 landuse classes onto tile class table 4.243
    DATA iluc_GLC2000       / 1 , &  ! 1 :evergreen broadleaf forest
                            & 2 , &  ! 2 :deciduous broadleaf closed forest 
                            & 3 , &  ! 3 :deciduous broadleaf open forest
                            & 4 , &  ! 4 :evergreen needleleaf forest
                            & 5 , &  ! 5 :deciduous needleleaf forest           
                            & 6 , &  ! 6 :mixed leaf trees           
                            & 7 , &  ! 7 :fresh water flooded trees             
                            & 8 , &  ! 8 :saline water flooded trees          
                            & 9 , &  ! 9 :mosaic tree / natural vegetation            
                            & 10, &  ! 10:burnt tree cover     
                            & 11, &  ! 11:evergreen shrubs closed-open
                            & 12, &  ! 12:decidous shrubs closed-open
                            & 13, &  ! 13:herbaceous vegetation closed-open                      
                            & 14, &  ! 14:sparse herbaceous or grass          
                            & 15, &  ! 15:flooded shrubs or herbaceous                             
                            & 16, &  ! 16:cultivated & managed areas        
                            & 17, &  ! 17:mosaic crop / tree / natural vegetation
                            & 18, &  ! 18:mosaic crop / shrub / grass    
                            & 19, &  ! 19:bare areas                    
                            & 20, &  ! 20:water                                     
                            & 21, &  ! 21:snow & ice                                   
                            & 22, &  ! 22:artificial surface                         
                            & 39  /  ! 23:undefined
    !-----------------------------------------------------------------------

    IF (PRESENT(opt_linverse)) THEN
      linverse = opt_linverse
    ELSE
      linverse = .FALSE.
    ENDIF


    ! Special treatment for undefined points
    IF (iluc_in == -1) THEN
      iluc_out = 39                !undefined
      RETURN
    ENDIF
    IF (iluc_in == 0) THEN
      iluc_out = -9999       ! warning
      RETURN
    ENDIF


    SELECT CASE(lc_datbase)
    CASE(GLOBCOVER2009)
      IF (.NOT. linverse) THEN
        ! internal to GRIB2
        iluc_out = iluc_GLOBCOVER2009(iluc_in)
      ELSE
        ! GRIB2 to internal
!!$        iluc_out = find_loc(iluc_in,iluc_GLOBCOVER2009)
      ENDIF

    CASE(GLC2000)

      IF (.NOT. linverse) THEN
        ! internal to GRIB2
        iluc_out = iluc_GLC2000(iluc_in)
      ELSE
        ! GRIB2 to internal
!!$        iluc_out = find_loc(iluc_in,iluc_GLC2000)
      ENDIF

    END SELECT 


  END FUNCTION convert_luc_ICON2GRIB

END MODULE mo_lnd_nwp_config
