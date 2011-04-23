!>
!! Definition, allocation/deallocation and reading of external datasets
!!
!! This module contains the type-declaration for the external datasets.
!! including subroutines for memory allocation/deallocation and reading.
!!
!! @author Daniel Reinert, DWD
!! @author Hermann Asensio, DWD
!!
!!
!! @par Revision History
!! Initial revision by Daniel Reinert, DWD (2010-07-12)
!! Modification by Hermann Asensio, DWD (2010-07-16)
!!  - add miscellaneous variables for external parameters
!!
!! @par Copyright
!! 2002-2010 by DWD and MPI-M
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
MODULE mo_ext_data

  USE mo_kind
  USE mo_io_units,           ONLY: filename_max
  USE mo_run_nml,            ONLY: nproma, i_cell_type, itopo, locean, &
    &                              iforcing, inwp, fac_smooth_topo, n_iter_smooth_topo
  USE mo_model_domain,       ONLY: t_patch
  USE mo_impl_constants,     ONLY: SUCCESS, max_char_length
  USE mo_exception,          ONLY: message, finish
  USE mo_model_domimp_setup, ONLY: reshape_real, reshape_int
  USE mo_grid_nml,           ONLY: n_dom
  USE mo_interpolation,      ONLY: t_int_state, cells2verts_scalar
  USE mo_math_operators,     ONLY: nabla4_scalar!, nabla2_scalar
  USE mo_loopindices,        ONLY: get_indices_c
  USE mo_sync,               ONLY: SYNC_C, SYNC_V, sync_patch_array
  USE mo_mpi,                ONLY: p_pe, p_io, p_bcast
  USE mo_parallel_nml,       ONLY: p_test_run, p_comm_work_test, p_comm_work
  USE mo_communication,      ONLY: idx_no, blk_no

!  USE mo_impl_constants_grf, ONLY: grf_bdywidth_c

  IMPLICIT NONE

  ! required for testing/reading topography
  INCLUDE 'netcdf.inc'

  PRIVATE

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

  PUBLIC :: t_external_data
  PUBLIC :: t_external_atmos
  PUBLIC :: t_external_atmos_td
  PUBLIC :: t_external_ocean
!  PUBLIC :: external_ocean_td
  PUBLIC :: ext_data

  PUBLIC :: init_ext_data
  PUBLIC :: construct_ext_data_atm  ! PUBLIC attribute only necessary for postpro.f90
  PUBLIC :: destruct_ext_data_atm
  PUBLIC :: destruct_ext_data_oce


  !>
  !! atmosphere external data class
  !!
  !! atmosphere external data class
  !!
  TYPE :: t_external_atmos

    ! *** Topography ***
    REAL(wp), ALLOCATABLE ::   &  !< topographic height at cell centers      [m]
      &  topography_c(:,:)        ! index1=1,nproma, index2=1,nblks_c

    REAL(wp), ALLOCATABLE ::   &  !< smoothed topographic height at cell centers [m]
      &  topography_smt_c(:,:)    ! index1=1,nproma, index2=1,nblks_c

    REAL(wp), ALLOCATABLE ::   &  !< topographic height at cell edges        [m]
      &  topography_e(:,:)        ! index1=1,nproma, index2=1,nblks_e

    REAL(wp), ALLOCATABLE ::   &  !< topographic height at cell vertices     [m]
      &  topography_v(:,:)        ! index1=1,nproma, index2=1,nblks_v

    REAL(wp), ALLOCATABLE ::   &  !< smoothed topographic height at vertices [m]
      &  topography_smt_v(:,:)    ! index1=1,nproma, index2=1,nblks_v

    REAL(wp), ALLOCATABLE ::   &  !< geometric height times grav        [m**2/s**2]
      &  fis(:,:)                 ! index1=1,nproma, index2=1,nblks_c


    ! *** Land-Sea-Mask ***
    INTEGER, ALLOCATABLE ::    &  !< land-sea-mask for cell centers          [ ]
      &  lsm_atm_c(:,:)           ! index1=1,nproma, index2=1,nblks_c

    INTEGER, ALLOCATABLE ::    &  !< land-sea-mask for cell edges            [ ]
      &  lsm_atm_e(:,:)           ! index1=1,nproma, index2=1,nblks_e

    INTEGER, ALLOCATABLE ::    &  !< land-sea-mask for cell vertices         [ ]
      &  lsm_atm_v(:,:)           ! index1=1,nproma, index2=1,nblks_v


    REAL(wp), ALLOCATABLE ::    &  !< fraction land in a grid element        [ ]
      &  fr_land(:,:)              ! 0. for water, 1.0 indicates 100% land
                                   ! index1=1,nproma, index2=1,nblks_c

    REAL(wp), ALLOCATABLE ::    &  !< fraction land glacier in a grid element [ ]
      &  fr_glac(:,:)              ! 1.0 indicates 100% glacier
                                   ! index1=1,nproma, index2=1,nblks_c    
    REAL(wp), ALLOCATABLE ::    &  !< fraction land in a grid element        [ ]
      &  fr_land_smt(:,:)          !  = smoothed fr_land

    REAL(wp), ALLOCATABLE ::    &  !< fraction land glacier in a grid element [ ]
      &  fr_glac_smt(:,:)          ! = smoothed fr_glac

    
   ! *** roghness length ***
    REAL(wp), ALLOCATABLE ::    &  !< surface roughness                      [m]
      &  z0(:,:)                   ! index1=1,nproma, index2=1,nblks_c


    ! *** FLake ***
    REAL(wp), ALLOCATABLE ::    &  !< fraction of fresh water                [ ]
      &  fr_lake(:,:)              ! as partition of total area of the
                                   ! grid element
                                   ! index1=1,nproma, index2=1,nblks_c

    REAL(wp), ALLOCATABLE ::    &  !< lake depth                             [m]
      &  depth_lk(:,:)             ! index1=1,nproma, index2=1,nblks_c

    ! *** subgrid scale orography ***
    REAL(wp), ALLOCATABLE ::    &  !< standard deviation of sub-grid scale orography [m]
      &  sso_stdh(:,:)             ! index1=1,nproma, index2=1,nblks_c

    REAL(wp), ALLOCATABLE ::    &  !< anisotropy of sub-grid scale orography  [ ]
      &  sso_gamma(:,:)            ! index1=1,nproma, index2=1,nblks_c

    REAL(wp), ALLOCATABLE ::    &  !< angle betw. principal axis of orography and E [rad]
      &  sso_theta(:,:)            ! index1=1,nproma, index2=1,nblks_c

    REAL(wp), ALLOCATABLE ::    &  !< mean slope of sub-grid scale orography [ ]
      &  sso_sigma(:,:)            ! index1=1,nproma, index2=1,nblks_c


    ! *** vegetation parameters ***

    REAL(wp), ALLOCATABLE ::    &  !< ground fraction covered by plants (vegetation period)  [ ]
      & plcov_mx(:,:)              ! index1=1,nproma, index2=1,nblks_c


    REAL(wp), ALLOCATABLE ::    &  !< leaf area index (vegetation period)   [ ]
      &  lai_mx(:,:)               ! index1=1,nproma, index2=1,nblks_c


    REAL(wp), ALLOCATABLE ::    &  !< root depth  [ m ]
      &  root_dp(:,:)              ! index1=1,nproma, index2=1,nblks_c


    REAL(wp), ALLOCATABLE ::    &  !< ground fraction covered by evergreen forest   [ ]
      &  forest_e(:,:)                ! index1=1,nproma, index2=1,nblks_c


    REAL(wp), ALLOCATABLE ::    &  !< ground fraction covered by deciduous forest   [ ]
      &  forest_d(:,:)                ! index1=1,nproma, index2=1,nblks_c


    REAL(wp), ALLOCATABLE ::    &  !< minimum value of stomata resistance  [ s/m ]
      &  plant_res_min(:,:)        ! index1=1,nproma, index2=1,nblks_c


  !  REAL(wp), ALLOCATABLE ::   &  !<  ratio of NDVI to annual maximum NDVI [ ]
  !    &  ndvi_ratio(:,:)          ! index1=1,nproma, index2=1,nblks_c


    ! *** soil parameters ***
    INTEGER, ALLOCATABLE ::   &  !<  soil texture, keys 0-9       []
      &  soiltyp(:,:)            ! index1=1,nproma, index2=1,nblks_c

    REAL(wp), ALLOCATABLE ::   &  !< Near surface temperature (climatological mean)  [ K ]
      &  t_clim(:,:)              !  used as climatological layer (deepest layer) of T_SO
                                  ! index1=1,nproma, index2=1,nblks_c


    ! *** parameters for radiation ***
    REAL(wp), ALLOCATABLE ::   &  !<         []
      &  opt_thick(:,:)           ! index1=1,nproma, index2=1,nblks_c

  END TYPE t_external_atmos



  !>
  !! atmosphere external data class (time dependent)
  !!
  !! Contains auxiliary time dependent versions of some external atmospheric data
  !! fields already which are already defined in external_atmos. These fields will
  !! be used to store e.g. montly means from which updated external data can be
  !! derived. The updated interpolated fields can be copied into the time independent
  !! counterparts, which are defined in external_atmos.
  !!
  TYPE :: t_external_atmos_td

    REAL(wp), ALLOCATABLE ::   &  !<  (monthly) ratio of NDVI to annual maximum NDVI [ ]
      &  ndvi_ratio_td(:,:,:)     ! index1=1,nproma, index2=1,nblks_c, index3=1,ntimes

    REAL(wp), ALLOCATABLE ::   &  !<         []
      &  opt_thick_td(:,:,:)      ! index1=1,nproma, index2=1,nblks_c, index3=1,ntimes

  END TYPE t_external_atmos_td



  !>
  !! ocean external data class
  !!
  !! ocean external data class
  !!
  TYPE :: t_external_ocean

    REAL(wp), ALLOCATABLE ::   &  !<  topographic height at cell centers  [m]
      &  bathymetry_c(:,:)        !  index1=1,nproma, index2=1,nblks_c

    REAL(wp), ALLOCATABLE ::   &  !<  topographic height at cell edges    [m]
      &  bathymetry_e(:,:)        ! index1=1,nproma, index2=1,nblks_e

    REAL(wp), ALLOCATABLE ::   &  !<  topographic height at cell vertices [m]
      &  bathymetry_v(:,:)        ! index1=1,nproma, index2=1,nblks_v

    INTEGER, ALLOCATABLE ::    &  !< land-sea-mask for cell centers
      &  lsm_oce_c(:,:,:)         ! index1=1,nproma, index2=1,n_zlev, index3=1,nblks_c

    INTEGER, ALLOCATABLE ::    &  !< land-sea-mask for cell edges
      &  lsm_oce_e(:,:,:)         ! index1=1,nproma, index2=1,n_zlev, index3=1,nblks_e

    INTEGER, ALLOCATABLE ::    &  !< land-sea-mask for cell vertices
      &  lsm_oce_v(:,:,:)         ! index1=1,nproma, index2=1,n_zlev, index3=1,nblks_v

    !  Pointer to array that contains indices of boundary cells
    !  and number of boundary cells. A cell is a boundary cell, if
    !  not all of its surrounding cells in the primal grid are of
    !  the same type (land/sea).
    !  The indices are local, i.e. with respect to the patch.
    !  Ocean cells have maximal 2 boundary edges
    !
    ! index1=1,n_land_sea_boundary_c
    INTEGER, ALLOCATABLE :: land_sea_boundary_idx_c(:)
    INTEGER, ALLOCATABLE :: land_sea_boundary_blk_c(:)
    INTEGER :: n_land_sea_boundary_c

    !
    !  Pointer to array that contains indices of edges at lateral boundaries
    !  and number of boundary edges. Edge is a boundary edge, if its two
    !  adjacent cells are not of the same type (land/sea).
    !  The indices are local, i.e. with respect to the patch.
    !
    ! index1=1,n_land_sea_boundary_e
    INTEGER, ALLOCATABLE :: land_sea_boundary_idx_e(:)
    INTEGER, ALLOCATABLE :: land_sea_boundary_blk_e(:)
    INTEGER :: n_land_sea_boundary_e

    !
    !  Pointer to array that contains indices of ocean vertices and
    !  number of boundary vertices. Vertex is boundary vertex if
    !  not all of its surrounding cells in the dual grid are of
    !  the same type.
    !  The indices are local, i.e. with respect to the patch
    !  Ocean vertices have cells of different type among its neighbor cells.
    !
    ! index1=1,n_land_sea_boundary_v
    INTEGER, ALLOCATABLE :: land_sea_boundary_idx_v(:)
    INTEGER, ALLOCATABLE :: land_sea_boundary_blk_v(:)
    INTEGER :: n_land_sea_boundary_v

  END TYPE t_external_ocean



!  !>
!  !! ocean external data class (time dependent)
!  !!
!  !! This data type contains additional auxiliary time dependent versions of
!  !! some external oceanic data fields already defined in external_ocean. These
!  !! fields will be used to store e.g. montly means from which interpolated external
!  !! data can be derived. The updated fields are copied into the time independent
!  !! counterparts which are defined in external_ocean.
!  !!
!  TYPE :: external_ocean_td
!  END TYPE external_ocean_td



  !>
  !! External data class
  !!
  !! External data class
  !!
  TYPE :: t_external_data

    TYPE(t_external_atmos)    :: atm
    TYPE(t_external_atmos_td) :: atm_td

    TYPE(t_external_ocean) :: oce
    TYPE(t_external_ocean) :: oce_td

  END TYPE t_external_data


  TYPE(t_external_data), ALLOCATABLE :: &
    &  ext_data(:)  ! n_dom

!-------------------------------------------------------------------------

CONTAINS


  !-------------------------------------------------------------------------
  !>
  !! Init external data for atmosphere and ocean
  !!
  !! Init external data for atmosphere and ocean.
  !! 1. Memory is allocated.
  !! 2. External data are read in from netCDF file (optional)
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2010-07-16)
  !!
  SUBROUTINE init_ext_data (p_patch, p_int_state, ext_data)

    TYPE(t_patch), INTENT(IN)            :: p_patch(:)
    TYPE(t_int_state), INTENT(IN)        :: p_int_state(n_dom)
    TYPE(t_external_data), INTENT(INOUT) :: ext_data(:)

    INTEGER :: jg

    CHARACTER(len=max_char_length), PARAMETER :: &
      routine = 'mo_ext_data:init_ext_data'

!-------------------------------------------------------------------------

    ! Allocate memory for atmospheric external data
    CALL construct_ext_data_atm (p_patch, ext_data)

! #Daniel#: This still has to be coded
    ! Allocate memory for oceanic external data
!DR    IF ( locean ) THEN
!DR      CALL construct_ext_data_oce (p_patch, ext_data)
!DR    ENDIF


    ! Check, whether external data should be read from file
    ! currently this is done via 'itopo' (provisional).
    SELECT CASE(itopo)

    CASE(0) ! do not read external data
            ! topography from analytical functions

      CALL message( TRIM(routine),'Running with analytical topography' )

    CASE(1) ! read external data from netcdf dataset

      IF ( .NOT.locean ) THEN
        CALL message( TRIM(routine),'Running with atmosphere topography' )

        CALL read_ext_data_atm (p_patch, ext_data)
        IF (n_iter_smooth_topo > 0) THEN
          DO jg = 1, n_dom
            CALL smooth_topography (p_patch(jg), p_int_state(jg))
          ENDDO
        ENDIF
      ELSE
        CALL finish(TRIM(routine),&
          &    'OCEAN tried to read in atmospheric topography data')

!DR        CALL read_ext_data_atm (p_patch, ext_data)
!DR        CALL read_ext_data_oce (p_patch, ext_data)
      ENDIF

    CASE DEFAULT

      CALL finish( TRIM(routine), 'topography selection not supported' )

    END SELECT


  END SUBROUTINE init_ext_data



  !-------------------------------------------------------------------------
  !>
  !! Allocate memory for atmospheric external data
  !!
  !! Allocate memory for atmospheric external data
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2010-07-12)
  !!
  SUBROUTINE construct_ext_data_atm (p_patch, ext_data)

    TYPE(t_patch),          INTENT(IN)    :: p_patch(:)
    TYPE(t_external_data),  INTENT(INOUT) :: ext_data(:)

    CHARACTER(len=max_char_length), PARAMETER :: &
      routine = 'mo_ext_data:construct_ext_data_atm'

    INTEGER :: nblks_c, nblks_e, nblks_v
    INTEGER :: ntimes
    INTEGER :: jg
    INTEGER :: ist  ! status variable
!-------------------------------------------------------------------------

    DO jg = 1,n_dom

      ! values for the blocking
      nblks_c = p_patch(jg)%nblks_c
      nblks_e = p_patch(jg)%nblks_e
      nblks_v = p_patch(jg)%nblks_v


      !topographic height at cells, edges and vertices
      !
      !cells
      ALLOCATE(ext_data(jg)%atm%topography_c(nproma,nblks_c),STAT=ist)
      IF (ist /= success) THEN
        CALL finish (routine,'allocating topography_c failed')
      ENDIF
      !
      !smoothed at cells
      ALLOCATE(ext_data(jg)%atm%topography_smt_c(nproma,nblks_c),STAT=ist)
      IF (ist /= success) THEN
        CALL finish (routine,'allocating topography_smt_c failed')
      ENDIF
      !
      !edges
      ALLOCATE(ext_data(jg)%atm%topography_e(nproma,nblks_e),STAT=ist)
      IF (ist /= success) THEN
        CALL finish (routine,'allocating topography_e failed')
      ENDIF
      !
      !vertices
      ALLOCATE(ext_data(jg)%atm%topography_v(nproma,nblks_v),STAT=ist)
      IF (ist /= success) THEN
        CALL finish (routine,'allocating topography_v failed')
      ENDIF
      !
      !smoothed at vertices
      ALLOCATE(ext_data(jg)%atm%topography_smt_v(nproma,nblks_v),STAT=ist)
      IF (ist /= success) THEN
        CALL finish (routine,'allocating topography_smt_v failed')
      ENDIF


      !
      ! !land sea mask for cells, edges and vertices
      !
      !cells
      ALLOCATE(ext_data(jg)%atm%lsm_atm_c(nproma,nblks_c),STAT=ist)
      IF (ist /= success) THEN
        CALL finish (routine,'allocating lsm_atm_c failed')
      ENDIF
      !
      !edges
      ALLOCATE(ext_data(jg)%atm%lsm_atm_e(nproma,nblks_e),STAT=ist)
      IF (ist /= success) THEN
        CALL finish (routine,'allocating lsm_atm_e failed')
      ENDIF
      !
      !vertices
      ALLOCATE(ext_data(jg)%atm%lsm_atm_v(nproma,nblks_v),STAT=ist)
      IF (ist /= success) THEN
        CALL finish (routine,'allocating lsm_atm_v failed')
      ENDIF

      ! fr_land
      ALLOCATE(ext_data(jg)%atm%fr_land(nproma,nblks_c),STAT=ist)
      IF (ist /= success) THEN
        CALL finish (routine,'allocating fr_land failed')
      ENDIF

      ! maybe the next three (fr_glac, fr_land_smt, fr_glac_smt)
      ! should be moved into corresponding if block
      ! fr_glac
      ALLOCATE(ext_data(jg)%atm%fr_glac(nproma,nblks_c),STAT=ist)
      IF (ist /= success) THEN
        CALL finish (routine,'allocating fr_glac failed')
      ENDIF

      ! fr_land_smt
      ALLOCATE(ext_data(jg)%atm%fr_land_smt(nproma,nblks_c),STAT=ist)
      IF (ist /= success) THEN
        CALL finish (routine,'allocating fr_land_smt failed')
      ENDIF

      ! fr_glac_smt
      ALLOCATE(ext_data(jg)%atm%fr_glac_smt(nproma,nblks_c),STAT=ist)
      IF (ist /= success) THEN
        CALL finish (routine,'allocating fr_glac_smt failed')
      ENDIF      

      !
      ! !geometric height times grav
      !
      ALLOCATE(ext_data(jg)%atm%fis(nproma,nblks_c),STAT=ist)
      IF (ist /= success) THEN
        CALL finish (routine,'allocating fis failed')
      ENDIF


      ! Several IF-statements are necessary in order to allocate only
      ! those fields which are necessary for the chosen
      ! parameterizations.

      ! external parameter for NWP forcing
      IF (iforcing == inwp) THEN

        ! *** roughness length ***
        ALLOCATE(ext_data(jg)%atm%z0(nproma,nblks_c),STAT=ist)
        IF (ist /= success) THEN
          CALL finish (routine,'allocating z0 failed')
        ENDIF

        ! *** FLake ***
        ALLOCATE(ext_data(jg)%atm%fr_lake(nproma,nblks_c),STAT=ist)
        IF (ist /= success) THEN
          CALL finish (routine,'allocating fr_lake failed')
        ENDIF

        ALLOCATE(ext_data(jg)%atm%depth_lk(nproma,nblks_c),STAT=ist)
        IF (ist /= success) THEN
          CALL finish (routine,'allocating depth_lk failed')
        ENDIF

        ! *** subgrid scale orography ***
        ALLOCATE(ext_data(jg)%atm%sso_stdh(nproma,nblks_c),STAT=ist)
        IF (ist /= success) THEN
          CALL finish (routine,'allocating sso_stdh failed')
        ENDIF

        ALLOCATE(ext_data(jg)%atm%sso_gamma(nproma,nblks_c),STAT=ist)
        IF (ist /= success) THEN
          CALL finish (routine,'allocating sso_gamma failed')
        ENDIF

        ALLOCATE(ext_data(jg)%atm%sso_theta(nproma,nblks_c),STAT=ist)
        IF (ist /= success) THEN
          CALL finish (routine,'allocating sso_theta failed')
        ENDIF

        ALLOCATE(ext_data(jg)%atm%sso_sigma(nproma,nblks_c),STAT=ist)
        IF (ist /= success) THEN
          CALL finish (routine,'allocating sso_sigma failed')
        ENDIF

        ! *** vegetation parameters ***
        ALLOCATE(ext_data(jg)%atm%plcov_mx(nproma,nblks_c),STAT=ist)
        IF (ist /= success) THEN
          CALL finish (routine,'allocating plcov_mx failed')
        ENDIF

        ALLOCATE(ext_data(jg)%atm%lai_mx(nproma,nblks_c),STAT=ist)
        IF (ist /= success) THEN
          CALL finish (routine,'allocating lai_mx failed')
        ENDIF

        ALLOCATE(ext_data(jg)%atm%root_dp(nproma,nblks_c),STAT=ist)
        IF (ist /= success) THEN
          CALL finish (routine,'allocating root_dp failed')
        ENDIF

        ALLOCATE(ext_data(jg)%atm%forest_e(nproma,nblks_c),STAT=ist)
        IF (ist /= success) THEN
          CALL finish (routine,'allocating forest_e failed')
        ENDIF

        ALLOCATE(ext_data(jg)%atm%forest_d(nproma,nblks_c),STAT=ist)
        IF (ist /= success) THEN
          CALL finish (routine,'allocating forest_d failed')
        ENDIF

        ALLOCATE(ext_data(jg)%atm%plant_res_min(nproma,nblks_c),STAT=ist)
        IF (ist /= success) THEN
          CALL finish (routine,'allocating plant_res_min failed')
        ENDIF

        ! *** soil parameters ***
        ALLOCATE(ext_data(jg)%atm%soiltyp(nproma,nblks_c),STAT=ist)
        IF (ist /= success) THEN
          CALL finish (routine,'allocating soiltyp failed')
        ENDIF

        ALLOCATE(ext_data(jg)%atm%t_clim(nproma,nblks_c),STAT=ist)
        IF (ist /= success) THEN
          CALL finish (routine,'allocating t_clim failed')
        ENDIF

        ! *** time dependent external parameters ***
        ! monthly mean data
        ntimes = 12

        ALLOCATE(ext_data(jg)%atm_td%ndvi_ratio_td(nproma,nblks_c,ntimes),STAT=ist)
        IF (ist /= success) THEN
          CALL finish (routine,'allocating ndvi_ratio_td failed')
        ENDIF


      ENDIF ! iforcing



! #Hermann#
! Echam physics external parameter have to be alllocated and initialized as well

! #Daniel#
! time dependent fields need to be allocated and initialized as well.
!

      !
      ! Initialize all fields with zero
      !
      ext_data(jg)%atm%topography_c(:,:)    = 0._wp
      ext_data(jg)%atm%topography_smt_c(:,:)= 0._wp
      ext_data(jg)%atm%topography_e(:,:)    = 0._wp
      ext_data(jg)%atm%topography_v(:,:)    = 0._wp
      ext_data(jg)%atm%topography_smt_v(:,:)= 0._wp
      ext_data(jg)%atm%lsm_atm_c(:,:)       = 0
      ext_data(jg)%atm%lsm_atm_e(:,:)       = 0
      ext_data(jg)%atm%lsm_atm_v(:,:)       = 0
      ext_data(jg)%atm%fr_land(:,:)         = 0._wp
      ext_data(jg)%atm%fr_glac(:,:)         = 0._wp
      ext_data(jg)%atm%fr_land_smt(:,:)     = 0._wp
      ext_data(jg)%atm%fr_glac_smt(:,:)     = 0._wp
      ext_data(jg)%atm%fis(:,:)             = 0._wp

      ! external parameter for NWP forcing
      IF (iforcing == inwp) THEN
        ext_data(jg)%atm%z0(:,:)            = 0._wp
        ext_data(jg)%atm%fr_lake(:,:)       = 0._wp
        ext_data(jg)%atm%depth_lk(:,:)      = 0._wp
        ext_data(jg)%atm%sso_stdh(:,:)      = 0._wp
        ext_data(jg)%atm%sso_gamma(:,:)     = 0._wp
        ext_data(jg)%atm%sso_theta(:,:)     = 0._wp
        ext_data(jg)%atm%sso_sigma(:,:)     = 0._wp
        ext_data(jg)%atm%plcov_mx(:,:)      = 0._wp
        ext_data(jg)%atm%lai_mx(:,:)        = 0._wp
        ext_data(jg)%atm%root_dp(:,:)       = 0._wp
        ext_data(jg)%atm%forest_e(:,:)      = 0._wp
        ext_data(jg)%atm%forest_d(:,:)      = 0._wp
        ext_data(jg)%atm%plant_res_min(:,:) = 0._wp
        ext_data(jg)%atm%t_clim(:,:)        = 0._wp
        ext_data(jg)%atm%soiltyp(:,:)       = 0

        ext_data(jg)%atm_td%ndvi_ratio_td(:,:,:)= 0._wp
      ENDIF ! iforcing

    END DO

  END SUBROUTINE construct_ext_data_atm



  !-------------------------------------------------------------------------
  !>
  !! Allocate memory for oceanic external data
  !!
  !! Allocate memory for oceanic external data
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2010-07-12)
  !!
  SUBROUTINE construct_ext_data_oce (p_patch, ext_data)

    TYPE(t_patch),           TARGET,INTENT(IN)    :: p_patch(:)  ! please do not remove
    TYPE(t_external_data),   TARGET,INTENT(INOUT) :: ext_data(:) ! please do not remove

    CHARACTER(len=max_char_length), PARAMETER :: &
      routine = 'mo_ext_data:construct_ext_data_oce'

!DR    INTEGER :: nblks_c, nblks_e, nblks_v
!DR    INTEGER :: jg
!DR    INTEGER :: ist  ! status variable
!-------------------------------------------------------------------------

    CALL finish(TRIM(routine),'sorry, not implemented yet')

!    DO jg = 1,n_dom

!      ! values for the blocking
!      nblks_c = p_patch(jg)%nblks_c
!      nblks_e = p_patch(jg)%nblks_e
!      nblks_v = p_patch(jg)%nblks_v

!     ALLOCATE(ext_data(jg)%oce%land_sea_boundary_idx_c(,),  &
!       &      ext_data(jg)%land_sea_boundary_blk_c(,),  &
!       &        STAT=ist)
!     IF(ist/=SUCCESS)THEN
!       CALL finish  (routine,'allocating land_sea_boundary_c failed')
!     ENDIF

!     ALLOCATE(ext_data(jg)%oce%land_sea_boundary_idx_e(,),  &
!       &      ext_data(jg)%oce%land_sea_boundary_blk_e(,),  &
!       &        STAT=ist)
!     IF(ist/=SUCCESS)THEN
!       CALL finish  (routine,'allocating land_sea_boundary_e failed')
!     ENDIF

!     ALLOCATE(ext_data(jg)%oce%land_sea_boundary_idx_v(,),  &
!       &      ext_data(jg)%oce%land_sea_boundary_idx_v(,),  &
!       &        STAT=ist)
!     IF(ist/=SUCCESS)THEN
!       CALL finish  (routine,'allocating land_sea_boundary_v failed')
!     ENDIF

!   ENDDO

  END SUBROUTINE construct_ext_data_oce



  !-------------------------------------------------------------------------
  !>
  !! Deallocate memory for atmospheric external data
  !!
  !! Deallocate memory for atmospheric external data
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2010-07-12)
  !!
  SUBROUTINE destruct_ext_data_atm (ext_data)

    TYPE(t_external_data), INTENT(INOUT) :: ext_data(:)

    CHARACTER(len=max_char_length), PARAMETER :: &
      routine = 'mo_ext_data:destruct_ext_data_atm'

    INTEGER :: jg
    INTEGER :: ist     ! status variable

!-------------------------------------------------------------------------

    DO jg = 1,n_dom

      !topographic height at cells, edges and vertices
      !
      !cells
      DEALLOCATE(ext_data(jg)%atm%topography_c,STAT=ist)
      IF (ist /= success) THEN
        CALL finish (routine,'deallocating topography_c failed')
      ENDIF
      !
      !smoothed at cells
      DEALLOCATE(ext_data(jg)%atm%topography_smt_c,STAT=ist)
      IF (ist /= success) THEN
        CALL finish (routine,'deallocating topography_smt_c failed')
      ENDIF
      !
      !edges
      DEALLOCATE(ext_data(jg)%atm%topography_e,STAT=ist)
      IF (ist /= success) THEN
        CALL finish (routine,'deallocating topography_e failed')
      ENDIF
      !
      !vertices
      DEALLOCATE(ext_data(jg)%atm%topography_v,STAT=ist)
      IF (ist /= success) THEN
        CALL finish (routine,'deallocating topography_v failed')
      ENDIF
      !
      !smoothed at vertices
      DEALLOCATE(ext_data(jg)%atm%topography_smt_v,STAT=ist)
      IF (ist /= success) THEN
        CALL finish (routine,'deallocating topography_smt_v failed')
      ENDIF


      !
      ! !land sea mask for cells, edges and vertices
      !
      !cells
      DEALLOCATE(ext_data(jg)%atm%lsm_atm_c,STAT=ist)
      IF (ist /= success) THEN
        CALL finish (routine,'deallocating lsm_atm_c failed')
      ENDIF
      !
      !edges
      DEALLOCATE(ext_data(jg)%atm%lsm_atm_e,STAT=ist)
      IF (ist /= success) THEN
        CALL finish (routine,'deallocating lsm_atm_e failed')
      ENDIF
      !
      !vertices
      DEALLOCATE(ext_data(jg)%atm%lsm_atm_v,STAT=ist)
      IF (ist /= success) THEN
        CALL finish (routine,'deallocating lsm_atm_v failed')
      ENDIF


      !
      !fr_land
      DEALLOCATE(ext_data(jg)%atm%fr_land,STAT=ist)
      IF (ist /= success) THEN
        CALL finish (routine,'deallocating fr_land failed')
      ENDIF
      !fr_glac
      DEALLOCATE(ext_data(jg)%atm%fr_glac,STAT=ist)
      IF (ist /= success) THEN
        CALL finish (routine,'deallocating fr_glac failed')
      ENDIF
      !fr_land_smt
      DEALLOCATE(ext_data(jg)%atm%fr_land_smt,STAT=ist)
      IF (ist /= success) THEN
        CALL finish (routine,'deallocating fr_land_smt failed')
      ENDIF
      !fr_glac_smt
      DEALLOCATE(ext_data(jg)%atm%fr_glac_smt,STAT=ist)
      IF (ist /= success) THEN
        CALL finish (routine,'deallocating fr_glac_smt failed')
      ENDIF

      !
      ! geometric height times grav
      DEALLOCATE(ext_data(jg)%atm%fis,STAT=ist)
      IF (ist /= success) THEN
        CALL finish (routine,'deallocating fis failed')
      ENDIF


      !
      ! external parameter for NWP forcing
      !
      IF (iforcing == inwp) THEN

        ! *** roghness length ***
        DEALLOCATE(ext_data(jg)%atm%z0,STAT=ist)
        IF (ist /= success) THEN
          CALL finish (routine,'deallocating z0 failed')
        ENDIF

        ! *** FLake ***
        DEALLOCATE(ext_data(jg)%atm%fr_lake,STAT=ist)
        IF (ist /= success) THEN
          CALL finish (routine,'deallocating fr_lake failed')
        ENDIF

        DEALLOCATE(ext_data(jg)%atm%depth_lk,STAT=ist)
        IF (ist /= success) THEN
          CALL finish (routine,'deallocating depth_lk failed')
        ENDIF

        ! *** subgrid scale orography ***
        DEALLOCATE(ext_data(jg)%atm%sso_stdh,STAT=ist)
        IF (ist /= success) THEN
          CALL finish (routine,'deallocating sso_stdh failed')
        ENDIF

        DEALLOCATE(ext_data(jg)%atm%sso_gamma,STAT=ist)
        IF (ist /= success) THEN
          CALL finish (routine,'deallocating sso_gamma failed')
        ENDIF

        DEALLOCATE(ext_data(jg)%atm%sso_sigma,STAT=ist)
        IF (ist /= success) THEN
          CALL finish (routine,'deallocating sso_sigma failed')
        ENDIF

        ! *** vegetation parameters ***
        DEALLOCATE(ext_data(jg)%atm%plcov_mx,STAT=ist)
        IF (ist /= success) THEN
          CALL finish (routine,'deallocating plcov_mx failed')
        ENDIF

        DEALLOCATE(ext_data(jg)%atm%lai_mx,STAT=ist)
        IF (ist /= success) THEN
          CALL finish (routine,'deallocating lai_mx failed')
        ENDIF

        DEALLOCATE(ext_data(jg)%atm%root_dp,STAT=ist)
        IF (ist /= success) THEN
          CALL finish (routine,'deallocating root_dp failed')
        ENDIF

        DEALLOCATE(ext_data(jg)%atm%forest_e,STAT=ist)
        IF (ist /= success) THEN
          CALL finish (routine,'deallocating forest_e failed')
        ENDIF

        DEALLOCATE(ext_data(jg)%atm%forest_d,STAT=ist)
        IF (ist /= success) THEN
          CALL finish (routine,'deallocating forest_d failed')
        ENDIF

        DEALLOCATE(ext_data(jg)%atm%plant_res_min,STAT=ist)
        IF (ist /= success) THEN
          CALL finish (routine,'deallocating plant_res_min failed')
        ENDIF

        ! *** soil parameters ***
        DEALLOCATE(ext_data(jg)%atm%soiltyp,STAT=ist)
        IF (ist /= success) THEN
          CALL finish (routine,'deallocating soiltyp failed')
        ENDIF

        DEALLOCATE(ext_data(jg)%atm%t_clim,STAT=ist)
        IF (ist /= success) THEN
          CALL finish (routine,'deallocating t_clim failed')
        ENDIF

        ! *** time dependent external parameters ***

        DEALLOCATE(ext_data(jg)%atm_td%ndvi_ratio_td,STAT=ist)
        IF (ist /= success) THEN
          CALL finish (routine,'deallocating ndvi_ratio_td failed')
        ENDIF


      ENDIF ! iforcing

! #Hermann#
! Echam physics external parameter have to be dealllocated as well

! #Daniel#
! time dependent fields need to be deallocated as well.
!

    END DO

  END SUBROUTINE destruct_ext_data_atm


  !-------------------------------------------------------------------------
  !>
  !! Deallocate memory for oceanic external data
  !!
  !! Deallocate memory for oceanic external data
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2010-07-12)
  !!
  SUBROUTINE destruct_ext_data_oce (ext_data)

    TYPE(t_external_data), TARGET, INTENT(INOUT) :: ext_data(:) ! please do not remove

    CHARACTER(len=max_char_length), PARAMETER :: &
      routine = 'mo_ext_data:destruct_ext_data_oce'

!DR    INTEGER :: jg
!DR    INTEGER :: ist     ! status variable

!-------------------------------------------------------------------------

    CALL finish(TRIM(routine),'sorry, not implemented yet')

!    DO jg = 1,n_dom

!     DEALLOCATE(ext_data(jg)%oce%land_sea_boundary_idx_c,  &
!       &        ext_data(jg)%oce%land_sea_boundary_blk_c,  &
!       &        STAT=ist)
!     IF(ist/=SUCCESS)THEN
!       CALL finish  (routine,'deallocating land_sea_boundary_c failed')
!     ENDIF

!     DEALLOCATE(ext_data(jg)%oce%land_sea_boundary_idx_e,  &
!       &        ext_data(jg)%oce%land_sea_boundary_blk_e,  &
!       &        STAT=ist)
!     IF(ist/=SUCCESS)THEN
!       CALL finish  (routine,'deallocating land_sea_boundary_e failed')
!     ENDIF

!     DEALLOCATE(ext_data(jg)%oce%land_sea_boundary_idx_v,  &
!       &        ext_data(jg)%oce%land_sea_boundary_idx_v,  &
!       &        STAT=ist)
!     IF(ist/=SUCCESS)THEN
!       CALL finish  (routine,'deallocating for land_sea_boundary_v failed')
!     ENDIF

!   ENDDO

  END SUBROUTINE destruct_ext_data_oce



  !-------------------------------------------------------------------------
  !>
  !! Read atmospheric external data
  !!
  !! Read atmospheric external data from netcdf
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2010-07-14)
  !!
  SUBROUTINE read_ext_data_atm (p_patch, ext_data)

    TYPE(t_patch), INTENT(IN)            :: p_patch(:)
    TYPE(t_external_data), INTENT(INOUT) :: ext_data(:)

    CHARACTER(len=max_char_length), PARAMETER :: &
      routine = 'mo_ext_data:read_ext_data_atm'

    CHARACTER(filename_max) :: topo_file  !< file name for reading in

    LOGICAL :: l_exist
    INTEGER :: jg
    INTEGER :: i_lev, no_cells, no_verts
    INTEGER :: ncid, dimid


!-------------------------------------------------------------------------

    DO jg = 1,n_dom

      i_lev = p_patch(jg)%level

      IF(p_pe == p_io) THEN
        !
        ! generate file name
        !
        WRITE(topo_file,'(a)') 'external_parameter_icon.nc'

        INQUIRE (FILE=topo_file, EXIST=l_exist)
        IF (.NOT.l_exist) THEN
          CALL finish(TRIM(routine),'Topography file is not found.')
        ENDIF

        !
        ! open file
        !
        CALL nf(nf_open(TRIM(topo_file), NF_NOWRITE, ncid))

        !
        ! get number of cells and vertices
        !
        CALL nf(nf_inq_dimid(ncid, 'cell', dimid))
        IF (i_cell_type == 3) THEN ! triangular grid
          CALL nf(nf_inq_dimlen(ncid, dimid, no_cells))
        ELSEIF (i_cell_type == 6) THEN ! hexagonal grid
          CALL nf(nf_inq_dimlen(ncid, dimid, no_verts))
        ENDIF

        CALL nf(nf_inq_dimid(ncid, 'vertex', dimid))
        IF (i_cell_type == 3) THEN ! triangular grid
          CALL nf(nf_inq_dimlen(ncid, dimid, no_verts))
        ELSEIF (i_cell_type == 6) THEN ! hexagonal grid
          CALL nf(nf_inq_dimlen(ncid, dimid, no_cells))
        ENDIF

        !
        ! check the number of cells and verts
        !
        IF(p_patch(jg)%n_patch_cells_g /= no_cells) THEN
          CALL finish(TRIM(ROUTINE),&
          & 'Number of patch cells and cells in topography file do not match.')
        ENDIF
        IF(p_patch(jg)%n_patch_verts_g /= no_verts) THEN
          CALL finish(TRIM(ROUTINE),&
          & 'Number of patch verts and verts in topography file do not match.')
        ENDIF
      ENDIF


      !-------------------------------------------------------
      !
      ! Read topography for triangle centers and vertices
      !
      !-------------------------------------------------------

      ! triangle center

      IF (i_cell_type == 3) THEN     ! triangular grid
        CALL read_netcdf_data (ncid, 'topography_c', p_patch(jg)%n_patch_cells_g,       &
          &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
          &                     ext_data(jg)%atm%topography_c)
      ELSEIF (i_cell_type == 6) THEN ! hexagonal grid
        CALL read_netcdf_data (ncid, 'topography_c', p_patch(jg)%n_patch_verts_g,       &
          &                     p_patch(jg)%n_patch_verts, p_patch(jg)%verts%glb_index, &
          &                     ext_data(jg)%atm%topography_v)
      ENDIF

      ! triangle vertex
      IF (i_cell_type == 3) THEN     ! triangular grid
        CALL read_netcdf_data (ncid, 'topography_v', p_patch(jg)%n_patch_verts_g,       &
          &                     p_patch(jg)%n_patch_verts, p_patch(jg)%verts%glb_index, &
          &                     ext_data(jg)%atm%topography_v)
      ELSEIF (i_cell_type == 6) THEN ! hexagonal grid
        CALL read_netcdf_data (ncid, 'topography_v', p_patch(jg)%n_patch_cells_g,       &
          &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
          &                     ext_data(jg)%atm%topography_c)
      ENDIF

      !other external parameters on triangular grid
      IF (i_cell_type == 3) THEN     ! triangular grid

        CALL read_netcdf_data (ncid, 'fr_land', p_patch(jg)%n_patch_cells_g,              &
          &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index,   &
          &                     ext_data(jg)%atm%fr_land)
        
        IF (iforcing == inwp) THEN
          CALL read_netcdf_data (ncid, 'plcov_mx', p_patch(jg)%n_patch_cells_g,           &
            &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
            &                     ext_data(jg)%atm%plcov_mx)

          CALL read_netcdf_data (ncid, 'lai_mx', p_patch(jg)%n_patch_cells_g,             &
            &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
            &                     ext_data(jg)%atm%lai_mx)

          CALL read_netcdf_data (ncid, 'rs_min', p_patch(jg)%n_patch_cells_g,             &
            &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
            &                     ext_data(jg)%atm%plant_res_min)

          CALL read_netcdf_data (ncid, 'for_d', p_patch(jg)%n_patch_cells_g,              &
            &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
            &                     ext_data(jg)%atm%forest_d)

          CALL read_netcdf_data (ncid, 'for_e', p_patch(jg)%n_patch_cells_g,              &
            &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
            &                     ext_data(jg)%atm%forest_e)

          CALL read_netcdf_data (ncid, 'root', p_patch(jg)%n_patch_cells_g,               &
            &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
            &                     ext_data(jg)%atm%root_dp)

          CALL read_netcdf_data (ncid, 'z0', p_patch(jg)%n_patch_cells_g,                 &
            &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
            &                     ext_data(jg)%atm%z0)

          CALL read_netcdf_data (ncid, 'sso_stdh', p_patch(jg)%n_patch_cells_g,           &
            &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
            &                     ext_data(jg)%atm%sso_stdh)

          CALL read_netcdf_data (ncid, 'sso_theta', p_patch(jg)%n_patch_cells_g,          &
            &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
            &                     ext_data(jg)%atm%sso_theta)

          CALL read_netcdf_data (ncid, 'sso_gamma', p_patch(jg)%n_patch_cells_g,          &
            &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
            &                     ext_data(jg)%atm%sso_gamma)


          CALL read_netcdf_data (ncid, 'sso_sigma', p_patch(jg)%n_patch_cells_g,          &
            &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
            &                     ext_data(jg)%atm%sso_sigma)


          CALL read_netcdf_data (ncid, 'tem_clim', p_patch(jg)%n_patch_cells_g,           &
            &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
            &                     ext_data(jg)%atm%t_clim)


          CALL read_netcdf_data_int (ncid, 'soiltype_fao', p_patch(jg)%n_patch_cells_g,   &
            &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
            &                     ext_data(jg)%atm%soiltyp)

        ENDIF ! (iforcing == inwp)
        
!      ELSEIF (i_cell_type == 6) THEN ! hexagonal grid

      ENDIF

      !-------------------------------------------------------
      !
      ! Read ...
      !
      !-------------------------------------------------------


      !
      ! close file
      !
      IF(p_pe == p_io) CALL nf(nf_close(ncid))

    ENDDO


!    CALL finish(TRIM(routine),'End of topo-test implementation.')

  END SUBROUTINE read_ext_data_atm
  !-------------------------------------------------------------------------



  !-------------------------------------------------------------------------
  !>
  !! Read dataset from netcdf file
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2010-07-14)
  !! Adapted for parallel runs by Rainer Johanni (2010-12-07)
  !!
  SUBROUTINE read_netcdf_data (ncid, varname, glb_arr_len, loc_arr_len, glb_index, var_out)

    CHARACTER(len=*), INTENT(IN)  ::  &  !< Var name of field to be read
      &  varname

    INTEGER, INTENT(IN) :: ncid          !< id of netcdf file
    INTEGER, INTENT(IN) :: glb_arr_len   !< length of 1D field (global)
    INTEGER, INTENT(IN) :: loc_arr_len   !< length of 1D field (local)
    INTEGER, INTENT(IN) :: glb_index(:)  !< Index mapping local to global

    REAL(wp), INTENT(INOUT) :: &         !< output field
      &  var_out(:,:)

    INTEGER :: varid, mpi_comm, j, jl, jb
    REAL(wp):: z_dummy_array(glb_arr_len)!< local dummy array
  !-------------------------------------------------------------------------

    ! Get var ID
    IF(p_pe==p_io) CALL nf(nf_inq_varid(ncid, TRIM(varname), varid))

    IF(p_test_run) THEN
      mpi_comm = p_comm_work_test
    ELSE
      mpi_comm = p_comm_work
    ENDIF

    ! I/O PE reads and broadcasts data

    IF(p_pe==p_io) CALL nf(nf_get_var_double(ncid, varid, z_dummy_array(:)))
    CALL p_bcast(z_dummy_array, p_io, mpi_comm)

    var_out(:,:) = 0._wp

    ! Set var_out from global data
    DO j = 1, loc_arr_len

      jb = blk_no(j) ! Block index in distributed patch
      jl = idx_no(j) ! Line  index in distributed patch

      var_out(jl,jb) = z_dummy_array(glb_index(j))
    ENDDO

  END SUBROUTINE read_netcdf_data



  !-------------------------------------------------------------------------
  !>
  !! Read dataset from netcdf file
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2010-07-14)
  !! Adapted for parallel runs by Rainer Johanni (2010-12-07)
  !!
  SUBROUTINE read_netcdf_data_int (ncid, varname, glb_arr_len, loc_arr_len, glb_index, var_out)

    CHARACTER(len=*), INTENT(IN)  ::  &  !< Var name of field to be read
      &  varname

    INTEGER, INTENT(IN) :: ncid          !< id of netcdf file
    INTEGER, INTENT(IN) :: glb_arr_len   !< length of 1D field (global)
    INTEGER, INTENT(IN) :: loc_arr_len   !< length of 1D field (local)
    INTEGER, INTENT(IN) :: glb_index(:)  !< Index mapping local to global

    INTEGER, INTENT(INOUT) :: &          !< output field
      &  var_out(:,:)

    INTEGER :: varid, mpi_comm, j, jl, jb
    INTEGER :: z_dummy_array(glb_arr_len)!< local dummy array
  !-------------------------------------------------------------------------

    ! Get var ID
    IF(p_pe==p_io) CALL nf(nf_inq_varid(ncid, TRIM(varname), varid))

    IF(p_test_run) THEN
      mpi_comm = p_comm_work_test
    ELSE
      mpi_comm = p_comm_work
    ENDIF

    ! I/O PE reads and broadcasts data

    IF(p_pe==p_io) CALL nf(nf_get_var_int(ncid, varid, z_dummy_array(:)))
    CALL p_bcast(z_dummy_array, p_io, mpi_comm)

    var_out(:,:) = 0

    ! Set var_out from global data
    DO j = 1, loc_arr_len

      jb = blk_no(j) ! Block index in distributed patch
      jl = idx_no(j) ! Line  index in distributed patch

      var_out(jl,jb) = z_dummy_array(glb_index(j))
    ENDDO

  END SUBROUTINE read_netcdf_data_int


  !-------------------------------------------------------------------------


  SUBROUTINE nf(status)

    INTEGER, INTENT(in) :: status

    IF (status /= nf_noerr) THEN
      CALL finish('mo_ext_data netCDF error', nf_strerror(status))
    ENDIF

  END SUBROUTINE nf


   SUBROUTINE smooth_topography (p_patch, p_int)

    TYPE(t_patch), INTENT(IN)            :: p_patch
    TYPE(t_int_state), INTENT(IN)        :: p_int

    ! local variables
    INTEGER  :: jg, jb, jc, iter, il
    INTEGER  :: i_startblk, nblks_c, i_startidx, i_endidx
    REAL(wp) :: z_topo(nproma,1,p_patch%nblks_c),z_nabla4_topo(nproma,1,p_patch%nblks_c), &
      &         z_topo_old(nproma,1,p_patch%nblks_c)
    REAL(wp) :: z_topo_v(nproma,1,p_patch%nblks_v)
    REAL(wp) :: zmaxtop,zmintop,z_topo_new


    jg = p_patch%id

    z_topo_v(:,1,:) = ext_data(jg)%atm%topography_v(:,:)
    z_nabla4_topo(:,1,:) = 0._wp

    i_startblk = p_patch%cells%start_blk(2,1)
    nblks_c    = p_patch%nblks_c

!    write(0,*) 'n_iter_smooth_topo=',n_iter_smooth_topo
!    write(0,*) 'fac_smooth_topo=',fac_smooth_topo

    DO iter = 1, n_iter_smooth_topo
      z_topo(:,1,:)   = ext_data(jg)%atm%topography_c(:,:)
      z_topo_old(:,1,:) = z_topo(:,1,:)

      CALL nabla4_scalar(z_topo, p_patch, p_int, z_nabla4_topo, &
        & opt_slev=1, opt_elev=1  )

      DO jb = i_startblk,nblks_c

        CALL get_indices_c(p_patch, jb, i_startblk, nblks_c, &
                           i_startidx, i_endidx, 2)

        DO jc = i_startidx, i_endidx

          !Limiter to avoid amplification of local extrema

          ! compute maximum (zmaxtop) and minimum (zmintop) of neighbor cells topography
          !set zmaxtop and zmintop to first neighbor's value
          zmaxtop = z_topo(p_patch%cells%neighbor_idx(jc,jb,1),1, &
              &          p_patch%cells%neighbor_blk(jc,jb,1))
          zmintop = zmaxtop
          DO il=2,i_cell_type

            IF ( z_topo(p_patch%cells%neighbor_idx(jc,jb,il),1, &
              &          p_patch%cells%neighbor_blk(jc,jb,il)) > &
              & zmaxtop ) THEN
              zmaxtop = z_topo(p_patch%cells%neighbor_idx(jc,jb,il),1, &
                &          p_patch%cells%neighbor_blk(jc,jb,il) )
            ENDIF
            IF ( z_topo(p_patch%cells%neighbor_idx(jc,jb,il),1, &
              &          p_patch%cells%neighbor_blk(jc,jb,il)) < &
              & zmintop ) THEN
              zmintop = z_topo(p_patch%cells%neighbor_idx(jc,jb,il),1, &
                &          p_patch%cells%neighbor_blk(jc,jb,il) )
            ENDIF
           !zmaxtop and zmintop are now max resp min of all neighbors
          ENDDO

          z_topo_new = z_topo (jc,1,jb) -  &
              fac_smooth_topo * z_nabla4_topo(jc,1,jb) * &
              p_patch%cells%area(jc,jb)*p_patch%cells%area(jc,jb)

          !If it was a local maximum in the old field, dont make it higher
          IF ( zmaxtop < z_topo (jc,1,jb) ) THEN
             IF ( z_nabla4_topo(jc,1,jb) < 0.0_wp ) CYCLE
          ENDIF
          !If it was a local minimum in the old field, dont make it lower:
          IF ( zmintop > z_topo (jc,1,jb) ) THEN
             IF ( z_nabla4_topo(jc,1,jb) > 0.0_wp ) CYCLE
          ENDIF

          !If it became a local maximum in the new field with regard to old neighbors, avoid it:
          IF (( zmaxtop < z_topo_new  ) .AND. &
            & ( z_nabla4_topo(jc,1,jb) < 0.0_wp )) THEN
            ext_data(jg)%atm%topography_c(jc,jb) = &
               & MAX(z_topo_old(jc,1,jb),zmaxtop)
            CYCLE
          ENDIF
          !If it became a local minimum in the new field with regard to old neighbors, avoid it:
          IF (( zmintop > z_topo_new  ) .AND. &
            & ( z_nabla4_topo(jc,1,jb) > 0.0_wp )) THEN
            ext_data(jg)%atm%topography_c(jc,jb) = &
               & MIN(z_topo_old(jc,1,jb),zmintop)
            CYCLE
          ENDIF

          ext_data(jg)%atm%topography_c(jc,jb)=z_topo_new

        ENDDO

      ENDDO

      z_topo(:,1,:)   = ext_data(jg)%atm%topography_c(:,:)
      CALL sync_patch_array(SYNC_C, p_patch, z_topo)
      ext_data(jg)%atm%topography_c(:,:)=z_topo(:,1,:)

      DO jb = i_startblk,nblks_c

        CALL get_indices_c(p_patch, jb, i_startblk, nblks_c, &
                           i_startidx, i_endidx, 2)
        DO jc = i_startidx, i_endidx

          !Limiter to avoid amplification of local extrema
          zmaxtop = z_topo(p_patch%cells%neighbor_idx(jc,jb,1),1, &
              &          p_patch%cells%neighbor_blk(jc,jb,1))
          zmintop = zmaxtop
          DO il=2,i_cell_type

            IF ( z_topo(p_patch%cells%neighbor_idx(jc,jb,il),1, &
              &          p_patch%cells%neighbor_blk(jc,jb,il)) > &
              & zmaxtop ) THEN
              zmaxtop = z_topo(p_patch%cells%neighbor_idx(jc,jb,il),1, &
                &          p_patch%cells%neighbor_blk(jc,jb,il) )
            ENDIF
            IF ( z_topo(p_patch%cells%neighbor_idx(jc,jb,il),1, &
              &          p_patch%cells%neighbor_blk(jc,jb,il)) < &
              & zmintop ) THEN
              zmintop = z_topo(p_patch%cells%neighbor_idx(jc,jb,il),1, &
                &          p_patch%cells%neighbor_blk(jc,jb,il) )
            ENDIF
           !zmaxtop and zmintop are now mx resp min of all neighbors
          ENDDO

           !If it became e local maximum in the new field, avoid it
           IF ( ( zmaxtop < z_topo(jc,1,jb) ) .AND. &
             & ( z_nabla4_topo(jc,1,jb) < 0.0_wp ) ) THEN
             ext_data(jg)%atm%topography_c(jc,jb) = &
               & MAX(z_topo_old(jc,1,jb),zmaxtop)
           !If it became e local minimum in the new field, avoid it
           ELSEIF ( ( zmintop > z_topo(jc,1,jb) ) .AND. &
               & ( z_nabla4_topo(jc,1,jb) > 0.0_wp ) ) THEN
             ext_data(jg)%atm%topography_c(jc,jb) = &
               & MIN(z_topo_old(jc,1,jb),zmintop)
           ENDIF

        ENDDO

      ENDDO

      z_topo(:,1,:)   = ext_data(jg)%atm%topography_c(:,:)
      CALL sync_patch_array(SYNC_C, p_patch, z_topo)
      ext_data(jg)%atm%topography_c(:,:)=z_topo(:,1,:)

    ENDDO !iter


    ! Interpolate smooth topography from cells to vertices
    z_topo(:,1,:)   = ext_data(jg)%atm%topography_c(:,:)
    CALL cells2verts_scalar(z_topo,p_patch,p_int%cells_aw_verts,z_topo_v,1,1)

    CALL sync_patch_array(SYNC_V,p_patch,z_topo_v)

     ext_data(jg)%atm%topography_v(:,:) = z_topo_v(:,1,:)

!    PRINT *,'Max topo_c=',MAXVAL(ext_data(jg)%atm%topography_c(:,:))
!    PRINT *,'Min topo_c=',MINVAL(ext_data(jg)%atm%topography_c(:,:))
!    PRINT *,'Max topo_v=',MAXVAL(ext_data(jg)%atm%topography_v(:,:))
!    PRINT *,'Min topo_v=',MINVAL(ext_data(jg)%atm%topography_v(:,:))

  END SUBROUTINE smooth_topography


END MODULE mo_ext_data

