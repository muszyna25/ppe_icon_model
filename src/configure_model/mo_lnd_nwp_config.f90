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

  PUBLIC :: nlev_soil, nztlev ,nlev_snow ,nsfc_subs, nsfc_snow
  PUBLIC :: lseaice,  llake, lmelt , lmelt_var ,   lmulti_snow 
  PUBLIC :: itype_gscp, itype_trvg ,    itype_evsl, itype_tran 
  PUBLIC :: itype_root, itype_heatcond, itype_hydbound  
  PUBLIC :: lstomata,   l2tls, lana_rho_snow, itype_subs 
  PUBLIC :: t_tiles, p_tiles

  PUBLIC :: configure_lnd_nwp
  PUBLIC :: destruct_tiles_arrays

  CHARACTER(len=*),PARAMETER,PRIVATE :: version = '$Id$'

  !--------------------------------------------------------------------------
  ! Basic configuration setup for NWP land
  !--------------------------------------------------------------------------

!  TYPE t_nwp_lnd_config

  ! namelist variables
  INTEGER ::  nztlev             !< time integration scheme
  INTEGER ::  nlev_snow          !< number of snow layers
  INTEGER ::  nsfc_subs          !< number of TILES
  INTEGER ::  nsfc_snow          !< number of static surface types 
                                 !< which can have snow as a tile
  INTEGER ::  itype_gscp         !< type of grid-scale precipitation physics
  INTEGER ::  itype_trvg         !< type of vegetation transpiration parameterization
  INTEGER ::  itype_evsl         !< type of parameterization of bare soil evaporation
  INTEGER ::  itype_tran         !< type of surface to atmospher transfer
  INTEGER ::  itype_root         !< type of root density distribution
  INTEGER ::  itype_heatcond     !< type of soil heat conductivity
  INTEGER ::  itype_hydbound     !< type of hydraulic lower boundary condition
  INTEGER ::  itype_subs         !< type of subscale surface treatment =1 MOSAIC, =2 TILE 


  LOGICAL ::  lseaice     !> forecast with sea ice model
  LOGICAL ::  llake       !! forecast with lake model FLake
  LOGICAL ::  lmelt       !! soil model with melting process
  LOGICAL ::  lmelt_var   !! freezing temperature dependent on water content
  LOGICAL ::  lmulti_snow !! run the multi-layer snow model
  LOGICAL ::  lstomata    !! map of minimum stomata resistance
  LOGICAL ::  l2tls       !! forecast with 2-TL integration scheme
  LOGICAL ::  lana_rho_snow !! if .TRUE., take rho_snow-values from analysis file 


  ! derived variables
  INTEGER ::  nlev_soil   !< number of soil layers (based on zml_soil in impl_constants)

!  END TYPE t_nwp_lnd_config

  !>
  !!
!  TYPE(t_nwp_lnd_config) :: nwp_lnd_config(max_dom)

  TYPE t_tiles
    INTEGER, ALLOCATABLE :: length(:)    !< dim: nblks_c
    INTEGER, ALLOCATABLE :: corrsp(:,:)  !< dim: nproma,nblks_c

    LOGICAL :: snow_tile                 !< whether it is a snow tile
    LOGICAL :: snowfree_tile             !< whether it is a snow-free tile -  
                                         !< a counterpart to a snow tile
    INTEGER :: conjunct                  !< index of a counterpart for a given 
                                         !< tile in the array of tiles
                                         !< (snow-free for snow tile, snow for 
                                         !< snow-free tile, itself for tiles-surface 
                                         !< types for which no explicit snow tile is 
                                         !< considered
    LOGICAL :: lake_tile                 !< whether it is a lake tile

  END TYPE t_tiles

  TYPE(t_tiles), TARGET, ALLOCATABLE :: p_tiles(:,:) 

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



    ! setup tile arrays
    !
    ALLOCATE (p_tiles(n_dom, nsfc_subs), stat=ist)
    IF (ist /= success) THEN
      CALL finish(TRIM(routine),'allocation for p_tiles failed')
    ENDIF
    CALL construct_tiles_arrays (p_patch, p_tiles, n_dom, nproma)

    DO jg = 1, n_dom
      DO isubs = 1, nsfc_subs - nsfc_snow
        p_tiles(jg,isubs)%snow_tile     = .FALSE.
        p_tiles(jg,isubs)%snowfree_tile = .FALSE.
        p_tiles(jg,isubs)%conjunct      = isubs
      END DO

      DO isubs = nsfc_subs - nsfc_snow + 1, nsfc_subs-1, 2
        p_tiles(jg,isubs  )%snow_tile     = .FALSE.
        p_tiles(jg,isubs+1)%snow_tile     = .TRUE.
        p_tiles(jg,isubs  )%snowfree_tile = .TRUE.
        p_tiles(jg,isubs+1)%snowfree_tile = .FALSE.
        p_tiles(jg,isubs  )%conjunct      = isubs+1
        p_tiles(jg,isubs+1)%conjunct      = isubs
      END DO
!!$    p_tiles(jg,:)%lake_tile = .FALSE.
!!$    IF(nsfc_subs .NE. nsfc_snow) THEN       !temporary
!!$      p_tiles(jg,2)%lake_tile = .TRUE.
!!$    END IF

!!$    DO ns = 1, nsfc_subs
!!$      pt_tiles%length(ns,jb) = 0
!!$
!!$      DO jc = i_startidx, i_endidx
!!$
!!$        IF(subsfrac(jc,1,ns) > frac_thres) THEN
!!$          pt_tiles%length(ns,jb) = pt_tiles%length(ns,jb) + 1
!!$          pt_tiles%corrsp(pt_tiles%length(ns,jb),ns,jb) = jc
!!$        END IF
!!$      END DO
!!$
!!$    END DO
    ENDDO  ! jg


  END SUBROUTINE configure_lnd_nwp


  !-------------------------------------------------------------------------
  !>
  !! Construction of tiles arrays.
  !!
  !! Initialization of components with zero.
  !!
  !! @par Revision History
  !! Initial release by Ekaterina Machulskaya, DWD (2011-07-26)
  !!
  !!
  SUBROUTINE construct_tiles_arrays (p_patch, p_tiles, n_dom, nproma)
! 
    TYPE(t_patch), TARGET, INTENT(IN)   :: p_patch(:)    

    ! arrays of correspondences between tiles and grid points
    TYPE(t_tiles), TARGET, INTENT(INOUT):: p_tiles(:,:) 
    INTEGER, INTENT(IN) :: n_dom
    INTEGER, INTENT(IN) :: nproma  

    INTEGER :: nblks_c, & ! number of cell blocks to allocate
      &        jg     , & ! index of domain
      &        ns     , & ! index of tile
      &        ist        ! status 
                 
!-----------------------------------------------------------------------

    DO jg = 1, n_dom
      DO ns = 1, nsfc_subs
  
      !determine size of arrays
      nblks_c = p_patch(jg)%nblks_c
  
      ! length
      ALLOCATE(p_tiles(jg,ns)%length(nblks_c), STAT = ist)
      IF (ist/=SUCCESS)THEN
        CALL finish('mo_lnd_state:construct_tiles_arrays', &
                    'allocation for length of arrays failed')
      ENDIF
      p_tiles(jg,ns)%length(:) = 0

      ! corrsp
      ALLOCATE(p_tiles(jg,ns)%corrsp(nproma,nblks_c), STAT = ist)
      IF (ist/=SUCCESS)THEN
        CALL finish('mo_lnd_state:construct_tiles_arrays', &
                    'allocation for correspondences arrays failed')
      ENDIF
      p_tiles(jg,ns)%corrsp(:,:) = 0

      ENDDO !nsfc_subs
    ENDDO !ndom

  END SUBROUTINE construct_tiles_arrays

  !-------------------------------------------------------------------------
  !>
  !! Destruction of tiles arrays.
  !!
  !! @par Revision History
  !! Initial release by Ekaterina Machulskaya, DWD (2011-07-26)
  !!
  !!
  SUBROUTINE destruct_tiles_arrays (p_tiles)

    ! arrays of correspondences between tiles and grid points
    TYPE(t_tiles), TARGET, INTENT(INOUT):: p_tiles(:,:)

    INTEGER :: jg     , & ! index of domain
      &        ns     , & ! index of tile
      &        ist        ! status

    INTEGER :: n_dom
                 
    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = 'mo_nwp_lnd_state:destruct_tiles_arrays'
!-----------------------------------------------------------------------

    ! get number of model domains
    n_dom = UBOUND(p_tiles,1)

    DO jg = 1, n_dom
      DO ns = 1, nsfc_subs

        DEALLOCATE(p_tiles(jg,ns)%length, STAT=ist)
          IF(ist/=SUCCESS)THEN
            CALL finish (TRIM(routine),  &
              &  'deallocation of grid points arrays length failed')
          ENDIF
        DEALLOCATE(p_tiles(jg,ns)%corrsp, STAT=ist)
          IF(ist/=SUCCESS)THEN
            CALL finish (TRIM(routine),  &
              &  'deallocation of arrays  th the correspondence between &
              &   grid points and one-dimensional arrays failed')
          ENDIF
      ENDDO !nsfc_subs
    ENDDO !ndom

  END SUBROUTINE destruct_tiles_arrays

END MODULE mo_lnd_nwp_config
