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
!! @par Copyright
!! 2002-2011 by DWD and MPI-M
!! This software is provided for non-commercial use only.
!! See the LICENSE and the WARRANTY conditions.
!!
!! @par License
!! The use of ICON is hereby granted free of charge for an unlimited time,
!! provided the following rules are accepted and applied:
!! <ol>
!! <li> You may use or modify this code for your own non commercial and non
!!    violent purposes.
!! <li> The code may not be re-distributed without the consent of the authors.
!! <li> The copyright notice and statement of authorship must appear in all
!!    copies.
!! <li> You accept the warranty conditions (see WARRANTY).
!! <li> In case you intend to use the code commercially, we oblige you to sign
!!    an according license agreement with DWD and MPI-M.
!! </ol>
!!
!! @par Warranty
!! This code has been tested up to a certain level. Defects and weaknesses,
!! which may be included in the code, do not establish any warranties by the
!! authors.
!! The authors do not make any warranty, express or implied, or assume any
!! liability or responsibility for the use, acquisition or application of this
!! software.
!!
MODULE mo_lnd_nwp_config

  USE mo_kind,            ONLY: wp
  USE mo_exception,       ONLY: finish
  USE mo_impl_constants,  ONLY: MAX_NTRACER, SUCCESS, MAX_CHAR_LENGTH, &
    &                           max_dom, zml_soil
  USE mo_model_domain,    ONLY: t_patch

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: nlev_soil, nlev_snow, ntiles_total, ntiles_lnd, ntiles_water, nlists_water
  PUBLIC :: frlnd_thrhld, frlndtile_thrhld, frlake_thrhld, frsea_thrhld
  PUBLIC :: lseaice,  llake, lmelt, lmelt_var, lmulti_snow, lsnowtile 
  PUBLIC :: itype_gscp, itype_trvg ,    itype_evsl, itype_tran 
  PUBLIC :: itype_root, itype_heatcond, itype_hydbound, idiag_snowfrac
  PUBLIC :: lstomata,   l2tls, lana_rho_snow, itype_subs 
  PUBLIC :: isub_water, isub_seaice

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
  INTEGER ::  ntiles_water       !< number of extra tiles for ocean and lakes
  INTEGER ::  nlists_water       !< number of extra index lists for ocean and lakes
  REAL(wp)::  frlnd_thrhld       !< fraction threshold for creating a land grid point
  REAL(wp)::  frlndtile_thrhld   !< fraction threshold for retaining the respective 
                                 !< tile for a grid point
  REAL(wp)::  frlake_thrhld      !< fraction threshold for creating a lake grid point
  REAL(wp)::  frsea_thrhld       !< fraction threshold for creating a sea grid point
  INTEGER ::  itype_gscp         !< type of grid-scale precipitation physics
  INTEGER ::  itype_trvg         !< type of vegetation transpiration parameterization
  INTEGER ::  itype_evsl         !< type of parameterization of bare soil evaporation
  INTEGER ::  itype_tran         !< type of surface to atmospher transfer
  INTEGER ::  itype_root         !< type of root density distribution
  INTEGER ::  itype_heatcond     !< type of soil heat conductivity
  INTEGER ::  itype_hydbound     !< type of hydraulic lower boundary condition
  INTEGER ::  itype_subs         !< type of subscale surface treatment =1 MOSAIC, =2 TILE 
  INTEGER ::  idiag_snowfrac     !< method for diagnosis of snow-cover fraction

  LOGICAL ::  lseaice     !> forecast with sea-ice model
  LOGICAL ::  llake       !! forecast with lake model FLake
  LOGICAL ::  lmelt       !! soil model with melting process
  LOGICAL ::  lmelt_var   !! freezing temperature dependent on water content
  LOGICAL ::  lmulti_snow !! run the multi-layer snow model
  LOGICAL ::  lstomata    !! map of minimum stomata resistance
  LOGICAL ::  l2tls       !! forecast with 2-TL integration scheme
  LOGICAL ::  lana_rho_snow !! if .TRUE., take rho_snow-values from analysis file 
  LOGICAL ::  lsnowtile   !! if .TRUE., snow is considered as a separate tile


  ! derived variables
  INTEGER ::  nlev_soil   !< number of soil layers (based on zml_soil in impl_constants)
  INTEGER ::  isub_water  !< (open) water points tile number
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
  SUBROUTINE configure_lnd_nwp(p_patch, n_dom, nproma)
  !
    TYPE(t_patch), INTENT(IN)  :: p_patch(:) 
    INTEGER      , INTENT(IN)  :: n_dom      !< number of model domains
    INTEGER      , INTENT(IN)  :: nproma
    INTEGER                    :: jg, isubs  !< loop index 
    INTEGER                    :: ist        !< status

    CHARACTER(len=*), PARAMETER::  &
      &  routine = 'mo_lnd_nwp_config: configure_lnd_nwp'
    !-----------------------------------------------------------------------

    ! number of soil layers
    ! Note that this number must be consistent with the number of entries 
    ! in zml_soil. zml_soil provides soil layer full level heights.
    nlev_soil = SIZE(zml_soil)-1  !< currently 7

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
      nlists_water = 0
    ELSE
      ! extra tiles for water points
      ntiles_water = 2 ! one for lake and ocean points
                       ! another one for seaice

      nlists_water = 3 ! one for ocean points
                       ! one for lake points
                       ! one for sea-ice points 
    ENDIF

    ! (open) water points tile number
    isub_water  = MAX(1,ntiles_total + ntiles_water - 1)

    ! sea-ice tile number
    isub_seaice = ntiles_total + ntiles_water

  END SUBROUTINE configure_lnd_nwp

END MODULE mo_lnd_nwp_config
