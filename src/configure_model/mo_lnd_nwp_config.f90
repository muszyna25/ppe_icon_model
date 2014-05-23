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
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_lnd_nwp_config

  USE mo_kind,            ONLY: wp
  USE mo_exception,       ONLY: finish
  USE mo_impl_constants,  ONLY: MAX_NTRACER, SUCCESS, MAX_CHAR_LENGTH, &
    &                           max_dom, zml_soil
  USE mo_model_domain,    ONLY: t_patch
  USE mo_io_units,        ONLY: filename_max

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: nlev_soil, nlev_snow, ntiles_total, ntiles_lnd, ntiles_water
  PUBLIC :: frlnd_thrhld, frlndtile_thrhld, frlake_thrhld, frsea_thrhld
  PUBLIC :: lseaice,  llake, lmelt, lmelt_var, lmulti_snow, lsnowtile, max_toplaydepth
  PUBLIC :: itype_trvg, itype_evsl, itype_lndtbl
  PUBLIC :: itype_root, itype_heatcond, itype_interception, &
             itype_hydbound, idiag_snowfrac
  PUBLIC :: lstomata,   l2tls, lana_rho_snow 
  PUBLIC :: isub_water, isub_lake, isub_seaice
  PUBLIC :: sstice_mode, sst_td_filename, ci_td_filename

  PUBLIC :: configure_lnd_nwp

  CHARACTER(len=*),PARAMETER,PRIVATE :: version = '$Id$'

  !--------------------------------------------------------------------------
  ! Basic configuration setup for NWP land
  !--------------------------------------------------------------------------

!  TYPE t_nwp_lnd_config

  ! namelist variables
  INTEGER ::  nlev_snow          !< number of snow layers
  INTEGER ::  ntiles_total       !< total number of TILES
  INTEGER ::  ntiles_lnd         !< number of static land surface types
  INTEGER ::  ntiles_water       !< number of extra tiles for ocean, seaice and lakes
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
  INTEGER ::  itype_hydbound     !< type of hydraulic lower boundary condition
  INTEGER ::  idiag_snowfrac     !< method for diagnosis of snow-cover fraction

  LOGICAL ::  lseaice     !> forecast with sea-ice model
  LOGICAL ::  llake       !! forecast with lake model FLake
  LOGICAL ::  lmelt       !! soil model with melting process
  LOGICAL ::  lmelt_var   !! freezing temperature dependent on water content
  LOGICAL ::  lmulti_snow !! run the multi-layer snow model
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
  INTEGER ::  nlev_soil   !< number of soil layers (based on zml_soil in impl_constants)
  INTEGER ::  isub_water  !< (open) water points tile number
  INTEGER ::  isub_lake   !< lake points tile number
  INTEGER ::  isub_seaice !< seaice tile number

!  END TYPE t_nwp_lnd_config

  !>
  !!
!  TYPE(t_nwp_lnd_config) :: nwp_lnd_config(max_dom)

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
    CHARACTER(len=*), PARAMETER::  &
      &  routine = 'mo_lnd_nwp_config: configure_lnd_nwp'
    !-----------------------------------------------------------------------

    ! number of soil layers
    ! Note that this number must be consistent with the number of entries 
    ! in zml_soil. zml_soil provides soil layer full level heights (confirmed).
    nlev_soil = SIZE(zml_soil)

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
    ! ASCII art for tile approach, showing the configuration for:
    ! - ntiles_lnd = 3
    ! - lsnowtile  = .TRUE.
    !
    !------------------------------------------------------------------!
    !      |      |      ||      |      |      ||      |      |        !
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

    ! (open) water points tile number
    isub_water  = MAX(1,ntiles_total + ntiles_water - 2)

    ! lake points tile number
    isub_lake   = MAX(1,ntiles_total + ntiles_water - 1)

    ! sea-ice tile number
    isub_seaice = ntiles_total + ntiles_water

  END SUBROUTINE configure_lnd_nwp

END MODULE mo_lnd_nwp_config
