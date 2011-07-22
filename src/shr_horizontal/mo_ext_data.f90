!>
!! Definition, allocation/deallocation and reading of external datasets
!!
!! This module contains the type-declaration for the external datasets,
!! including memory allocation/deallocation and reading.
!!
!! @author Daniel Reinert, DWD
!! @author Hermann Asensio, DWD
!!
!!
!! @par Revision History
!! Initial revision by Daniel Reinert, DWD (2010-07-12)
!! Modification by Hermann Asensio, DWD (2010-07-16)
!!  - add miscellaneous variables for external parameters
!! Modification by Daniel Reinert, DWD (2011-05-03)
!! - Memory allocation method changed from explicit allocation to Luis'
!!   infrastructure
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
  USE mo_parallel_config,  ONLY: nproma
  USE mo_impl_constants,     ONLY: inwp, ihs_ocean
  USE mo_run_config,         ONLY: iforcing
  USE mo_extpar_config,      ONLY: itopo, fac_smooth_topo, n_iter_smooth_topo
  USE mo_dynamics_config,    ONLY: iequations
  USE mo_radiation_config,   ONLY: irad_o3
  USE mo_model_domain,       ONLY: t_patch
  USE mo_impl_constants,     ONLY: MAX_CHAR_LENGTH
  USE mo_exception,          ONLY: message, message_text, finish
  USE mo_grid_config,        ONLY: n_dom
  USE mo_interpolation,      ONLY: t_int_state, cells2verts_scalar
  USE mo_math_operators,     ONLY: nabla4_scalar
  USE mo_loopindices,        ONLY: get_indices_c
  USE mo_sync,               ONLY: SYNC_C, SYNC_V, sync_patch_array
  USE mo_mpi,                ONLY: p_pe, p_io, p_bcast, p_comm_work_test, p_comm_work
  USE mo_parallel_config,  ONLY: p_test_run
  USE mo_communication,      ONLY: idx_no, blk_no
  USE mo_linked_list,        ONLY: t_var_list
  USE mo_var_list,           ONLY: default_var_list_settings, &
    &                              add_var,                   &
    &                              new_var_list,              &
    &                              delete_var_list
  USE mo_cf_convention
  USE mo_grib2
  USE mo_cdi_constants


  IMPLICIT NONE

  INTEGER::  nlev_pres, nlev_height, nmonths ! 

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
  PUBLIC :: construct_ext_data  ! PUBLIC attribute only necessary for postpro.f90
  PUBLIC :: destruct_ext_data
  PUBLIC :: smooth_topography, read_netcdf_data

  INTERFACE read_netcdf_data
    MODULE PROCEDURE read_netcdf_2d
    MODULE PROCEDURE read_netcdf_2d_int
    MODULE PROCEDURE read_netcdf_3d
    MODULE PROCEDURE read_netcdf_4d
  END INTERFACE read_netcdf_data


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

    REAL(wp), POINTER ::   &   !< smoothed topographic height at cell centers [m]
      &  topography_smt_c(:,:) ! index1=1,nproma, index2=1,nblks_c

    REAL(wp), POINTER ::   &   !< topographic height at cell edges        [m]
      &  topography_e(:,:)     ! index1=1,nproma, index2=1,nblks_e

    REAL(wp), POINTER ::   &   !< topographic height at cell vertices     [m]
      &  topography_v(:,:)     ! index1=1,nproma, index2=1,nblks_v

    REAL(wp), POINTER ::   &   !< smoothed topographic height at vertices [m]
      &  topography_smt_v(:,:) ! index1=1,nproma, index2=1,nblks_v

    REAL(wp), POINTER ::   &   !< geopotential (S)                        [m**2/s**2]
      &  fis(:,:)              ! index1=1,nproma, index2=1,nblks_c


    !
    ! *** Land-Sea-Mask ***
    INTEGER, POINTER  ::   &   !< land-sea-mask for cell centers          [ ]
      &  lsm_atm_c(:,:)        ! index1=1,nproma, index2=1,nblks_c

    INTEGER, POINTER  ::   &   !< land-sea-mask for cell edges            [ ]
      &  lsm_atm_e(:,:)        ! index1=1,nproma, index2=1,nblks_e

    INTEGER, POINTER  ::   &   !< land-sea-mask for cell vertices         [ ]
      &  lsm_atm_v(:,:)        ! index1=1,nproma, index2=1,nblks_v


    REAL(wp), POINTER ::   &   !< fraction land in a grid element         [ ]
      &  fr_land(:,:)          ! 0. for water, 1.0 indicates 100% land
                               ! index1=1,nproma, index2=1,nblks_c

    REAL(wp), POINTER ::    &  !< fraction land glacier in a grid element [ ]
      &  fr_glac(:,:)          ! 1.0 indicates 100% glacier
                               ! index1=1,nproma, index2=1,nblks_c    

    REAL(wp), POINTER ::   &   !< fraction sea ice cover in a grid element [ ]
      &  fr_ice(:,:)           ! 1.0 indicates 100% ice
                               ! index1=1,nproma, index2=1,nblks_c 
   
    REAL(wp), POINTER ::   &   !< fraction land in a grid element         [ ]
      &  fr_land_smt(:,:)      !  = smoothed fr_land

    REAL(wp), POINTER ::   &   !< fraction sea ice cover in a grid element [ ]
      &  fr_ice_smt(:,:)       ! = smoothed fr_ice

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


    !
    ! *** subgrid scale orography ***
    REAL(wp), POINTER ::   &   !< standard deviation of sub-grid scale orography [m]
      &  sso_stdh(:,:)         ! index1=1,nproma, index2=1,nblks_c

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

    REAL(wp), POINTER ::   &   !< leaf area index (vegetation period)     [ ]
      &  lai_mx(:,:)           ! index1=1,nproma, index2=1,nblks_c

    REAL(wp), POINTER ::   &   !< root depth                              [m]
      &  rootdp(:,:)           ! index1=1,nproma, index2=1,nblks_c

    REAL(wp), POINTER ::   &   !< ground fraction covered by evergreen forest [ ]
      &  for_e(:,:)            ! index1=1,nproma, index2=1,nblks_c

    REAL(wp), POINTER ::   &   !< ground fraction covered by deciduous forest [ ]
      &  for_d(:,:)            ! index1=1,nproma, index2=1,nblks_c

    REAL(wp), POINTER ::   &   !< urban area fraction                     [ ]
      &  urban(:,:)            ! index1=1,nproma, index2=1,nblks_c

    REAL(wp), POINTER ::   &   !< minimum value of stomata resistance     [ s/m ]
      &  rsmin(:,:)            ! index1=1,nproma, index2=1,nblks_c

    REAL(wp), POINTER ::   &   !< annual maximum NDVI                     [ ]
      &  ndvi_max(:,:)         ! index1=1,nproma, index2=1,nblks_c



    !
    ! *** soil parameters ***
    INTEGER, POINTER  ::   &   !< soil texture, keys 0-9                  [ ]
      &  soiltyp(:,:)          ! index1=1,nproma, index2=1,nblks_c

    INTEGER, POINTER  ::   &   !< soil texture, keys 0-9                  [ ]
      &  soiltyp_frac(:,:,:)   ! index1=1,nproma, index2=1,nblks_c, index3=1,nsfc_subs

    REAL(wp), POINTER ::   &   !< Near surface temperature (climatological mean)  [ K ]
      &  t_cl(:,:)             !  used as climatological layer (deepest layer) of T_SO
                               ! index1=1,nproma, index2=1,nblks_c

    REAL(wp), POINTER ::   &   !< longwave surface emissivity             [ ]
      &  emis_rad(:,:)         ! index1=1,nproma, index2=1,nblks_c

    REAL(wp), POINTER ::  &    !< Landuse class fraction                  [ ]
      & lu_class_fraction(:,:,:) ! index1=1,nproma, index2=1,nblks_c, index3=1,nclass_lu

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
    REAL(wp), POINTER ::   &   !< aerosol optical thickness of black carbon    [ ]
      &  o3(:,:,:,:)           ! index1=1,nproma, index2=nlev_pres,
                               ! index3=1,nblks_c, index4=1,ntimes
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



    !
    ! *** vegetation parameters ***
    REAL(wp), POINTER ::   &   !< normalized difference vegetation index [ ]
      &  ndvi(:,:,:)           !< (monthly mean)
                               ! index1=1,nproma, index2=1,nblks_c, index3=1,ntimes

    REAL(wp), POINTER ::   &   !< (monthly) proportion of actual value/maximum 
      &  ndvi_mrat(:,:)        !< normalized differential vegetation index   [ ]
                               ! index1=1,nproma, index2=1,nblks_c

  END TYPE t_external_atmos_td



  !>
  !! ocean external data class
  !!
  !! ocean external data class
  !!
  TYPE :: t_external_ocean

    REAL(wp), POINTER ::   &   !< topographic height at cell centers  [m]
      &  bathymetry_c(:,:)     !  index1=1,nproma, index2=1,nblks_c

    REAL(wp), POINTER ::   &   !< topographic height at cell edges    [m]
      &  bathymetry_e(:,:)     ! index1=1,nproma, index2=1,nblks_e

    REAL(wp), POINTER ::   &   !< topographic height at cell vertices [m]
      &  bathymetry_v(:,:)     ! index1=1,nproma, index2=1,nblks_v

    INTEGER, POINTER ::    &   !< land-sea-mask for cell centers
      &  lsm_oce_c(:,:,:)      ! index1=1,nproma, index2=1,n_zlev, index3=1,nblks_c

    INTEGER, POINTER ::    &   !< land-sea-mask for cell edges
      &  lsm_oce_e(:,:,:)      ! index1=1,nproma, index2=1,n_zlev, index3=1,nblks_e

    INTEGER, POINTER ::    &   !< land-sea-mask for cell vertices
      &  lsm_oce_v(:,:,:)      ! index1=1,nproma, index2=1,n_zlev, index3=1,nblks_v

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

!    TYPE(t_external_ocean_td) :: oce_td
!    TYPE(t_var_list), POINTER :: oce_td_list

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
  !! 1. Build data structure, including field lists and 
  !!    memory allocation.
  !! 2. External data are read in from netCDF file or set analytically
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

    !-------------------------------------------------------------------------
    !  1.  inquire external files for their data structure
    !-------------------------------------------------------------------------

    IF(irad_o3 == 3) THEN 
      CALL inquire_external_files(p_patch)
    ENDIF

    !------------------------------------------------------------------
    !  2.  construct external fields for the model
    !------------------------------------------------------------------

    ! top-level procedure for building data structures for 
    ! external data.
    CALL construct_ext_data(p_patch, ext_data)

    !-------------------------------------------------------------------------
    !  3.  read the data into the fields
    !-------------------------------------------------------------------------

    ! Check, whether external data should be read from file
    ! 
    SELECT CASE(itopo)

    CASE(0) ! do not read external data
            ! topography from analytical functions

      CALL message( TRIM(routine),'Running with analytical topography' )
      IF(irad_o3 == 3) THEN
        CALL message( TRIM(routine),'external ozone data required' )
        CALL read_ext_data_atm (p_patch, ext_data)
      ENDIF


      IF(iforcing == inwp) THEN
      !
      ! initalize external data with meaningful data, in the case that they 
      ! are not read in from file.
      DO jg = 1, n_dom
        ext_data(jg)%atm%fr_land(:,:)     = 1._wp   ! land fraction
        ext_data(jg)%atm%fr_land_smt(:,:) = 1._wp   ! land fraction (smoothed)
        ext_data(jg)%atm%plcov_mx(:,:)    = 0.5_wp  ! plant cover
        ext_data(jg)%atm%lai_mx(:,:)      = 3._wp   ! max Leaf area index
        ext_data(jg)%atm%rootdp(:,:)      = 1._wp   ! root depth
        ext_data(jg)%atm%rsmin(:,:)       = 150._wp ! minimal stomata resistence
        ext_data(jg)%atm%soiltyp(:,:)     = 8       ! soil type
      END DO

    ENDIF

    CASE(1) ! read external data from netcdf dataset

      IF ( iforcing/=ihs_ocean) THEN
        CALL message( TRIM(routine),'Start reading external data from file' )

        CALL read_ext_data_atm (p_patch, ext_data)
        IF (n_iter_smooth_topo > 0) THEN
          DO jg = 1, n_dom
            CALL smooth_topography (p_patch(jg), p_int_state(jg),  &
                                    ext_data(jg)%atm%topography_c, &
                                    ext_data(jg)%atm%topography_v)
          ENDDO
        ENDIF

       CALL message( TRIM(routine),'Finished reading external data' )

      ELSE
        CALL finish(TRIM(routine),&
          &    'OCEAN tried to read in atmospheric topography data')

!DR       CALL read_ext_data_oce (p_patch, ext_data)
      ENDIF

    CASE DEFAULT

      CALL finish( TRIM(routine), 'topography selection not supported' )

    END SELECT


  END SUBROUTINE init_ext_data



  !-------------------------------------------------------------------------
  !>
  !! Top-level procedure for building external data structure
  !!
  !! Top-level procedure for building external data structure
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2010-07-12)
  !!
  SUBROUTINE construct_ext_data (p_patch, ext_data)

    TYPE(t_patch),          INTENT(IN)    :: p_patch(:)
    TYPE(t_external_data),  INTENT(INOUT) :: ext_data(:)

    INTEGER :: jg

    CHARACTER(len=MAX_CHAR_LENGTH) :: listname

    CHARACTER(len=max_char_length), PARAMETER :: &
      routine = 'mo_ext_data:construct_ext_data'

!-------------------------------------------------------------------------


    CALL message (TRIM(routine), 'Construction of data structure for ' // &
      &                          'external data started')

    DO jg = 1, n_dom

      !
      ! Build external data list for constant in time atmospheric fields
      WRITE(listname,'(a,i2.2)') 'ext_data_atm_D',jg
      CALL new_ext_data_atm_list(p_patch(jg), ext_data(jg)%atm,       &
        &                        ext_data(jg)%atm_list, TRIM(listname))


      IF (iforcing==inwp) THEN
        !       
        ! Build external data list for time-dependent atmospheric fields
        WRITE(listname,'(a,i2.2)') 'ext_data_atm_td_D',jg
        CALL new_ext_data_atm_td_list(p_patch(jg), ext_data(jg)%atm_td,       &
          &                           ext_data(jg)%atm_td_list, TRIM(listname))
      END IF

      !
      IF (iforcing==ihs_ocean) THEN
        ! Build external data list for constant in time oceanic fields
        WRITE(listname,'(a,i2.2)') 'ext_data_oce_D',jg
        CALL new_ext_data_oce_list(p_patch(jg), ext_data(jg)%oce,       &
          &                        ext_data(jg)%oce_list, TRIM(listname))
      ENDIF
      !
      ! Build external data list for time-dependent oceanic fields
      ! ### to be done ###

    END DO

    CALL message (TRIM(routine), 'Construction of data structure for ' // &
      &                          'external data finished')

  END SUBROUTINE construct_ext_data




  !-------------------------------------------------------------------------
  !
  !
  !>
  !! Allocation of atmospheric external data structure
  !!
  !! Allocation of atmospheric external data structure (constant in time 
  !! elements).
  !!
  !! Initialization of elements with zero.
  !!
  !! @par Revision History
  !! Initial release by Daniel Reinert (2011-05-03)
  !! Statements that assign initial value added by Hui Wan (MPI-M, 2011-05-30)
  !!
  SUBROUTINE new_ext_data_atm_list ( p_patch, p_ext_atm, p_ext_atm_list, &
    &                                listname)
!
    TYPE(t_patch), TARGET, INTENT(IN)     :: & !< current patch
      &  p_patch

    TYPE(t_external_atmos), INTENT(INOUT) :: & !< current external data structure
      &  p_ext_atm 

    TYPE(t_var_list) :: p_ext_atm_list !< current external data list

    CHARACTER(len=*), INTENT(IN)      :: & !< list name
      &  listname

    TYPE(t_cf_var)    :: cf_desc
    TYPE(t_grib2_var) :: grib2_desc

    INTEGER :: nblks_c, &    !< number of cell blocks to allocate
      &        nblks_e, &    !< number of edge blocks to allocate
      &        nblks_v       !< number of vertex blocks to allocate

    INTEGER :: shape2d_c(2), shape2d_e(2), shape2d_v(2) !, shape3d_c(3)

    INTEGER :: ientr         !< "entropy" of horizontal slice
    !--------------------------------------------------------------

    !determine size of arrays
    nblks_c = p_patch%nblks_c
    nblks_e = p_patch%nblks_e
    nblks_v = p_patch%nblks_v


    ientr = 16   ! "entropy" of horizontal slice

    ! predefined array shapes
    shape2d_c = (/ nproma, nblks_c /)
    shape2d_e = (/ nproma, nblks_e /)
    shape2d_v = (/ nproma, nblks_v /)
!DR shape3d_c = (/ nproma, nblks_c, nsfc_subs /)


    !
    ! Register a field list and apply default settings
    !
    CALL new_var_list( p_ext_atm_list, TRIM(listname) )
    CALL default_var_list_settings( p_ext_atm_list,            &
                                  & lrestart=.TRUE.,           &
                                  & restart_type=FILETYPE_NC2  )


    ! topography height at cell center
    !
    ! topography_c  p_ext_atm%topography_c(nproma,nblks_c)
    cf_desc    = t_cf_var('surface_height', 'm', &
      &                   'geometric height of the earths surface above sea level')
    grib2_desc = t_grib2_var( 2, 0, 7, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_ext_atm_list, 'topography_c', p_ext_atm%topography_c,      &
      &           GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d_c )


  IF (iequations==3) THEN

    ! smoothed topography height at cell center
    !
    ! topography_smt_c  p_ext_atm%topography_smt_c(nproma,nblks_c)
    cf_desc    = t_cf_var('smoothed_surface_height', 'm', &
      &                   'smoothed geometric height of the earths surface above sea level')
    grib2_desc = t_grib2_var( 2, 0, 7, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_ext_atm_list, 'topography_smt_c', p_ext_atm%topography_smt_c, &
      &           GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d_c )


    ! topography height at edge midpoint
    !
    ! topography_e  p_ext_atm%topography_e(nproma,nblks_e)
    cf_desc    = t_cf_var('surface_height', 'm', &
      &                   'geometric height of the earths surface above sea level')
    grib2_desc = t_grib2_var( 2, 0, 7, ientr, GRID_REFERENCE, GRID_EDGE)
    CALL add_var( p_ext_atm_list, 'topography_e', p_ext_atm%topography_e, &
      &           GRID_UNSTRUCTURED_EDGE, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d_e)


    ! topography height at vertex
    !
    ! topography_v  p_ext_atm%topography_v(nproma,nblks_v)
    cf_desc    = t_cf_var('surface_height', 'm', &
      &                   'geometric height of the earths surface above sea level')
    grib2_desc = t_grib2_var( 2, 0, 7, ientr, GRID_REFERENCE, GRID_VERTEX)
    CALL add_var( p_ext_atm_list, 'topography_v', p_ext_atm%topography_v, &
      &           GRID_UNSTRUCTURED_VERT, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d_v )


    ! smoothed topography height at vertex
    !
    ! topography_smt_v  p_ext_atm%topography_smt_v(nproma,nblks_v)
    cf_desc    = t_cf_var('smoothed_surface_height', 'm', &
      &                   'smoothed geometric height of the earths surface above sea level')
    grib2_desc = t_grib2_var( 2, 0, 7, ientr, GRID_REFERENCE, GRID_VERTEX)
    CALL add_var( p_ext_atm_list, 'topography_smt_v', p_ext_atm%topography_smt_v, &
      &           GRID_UNSTRUCTURED_VERT, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d_v )


    ! land sea mask for cells
    !
    ! lsm_atm_c    p_ext_atm%lsm_atm_c(nproma,nblks_c)
    cf_desc    = t_cf_var('land_sea_mask_(cell)', '-', &
      &                   'land sea mask (cell)')
    grib2_desc = t_grib2_var( 2, 0, 0, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_ext_atm_list, 'lsm_atm_c', p_ext_atm%lsm_atm_c, &
      &           GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, &
      &          grib2_desc, ldims=shape2d_c, lrestart=.FALSE.  )


    ! land sea mask for edges
    !
    ! lsm_atm_e    p_ext_atm%lsm_atm_e(nproma,nblks_e)
    cf_desc    = t_cf_var('land_sea_mask_(edge)', '-', &
      &                   'land sea mask (edge)')
    grib2_desc = t_grib2_var( 2, 0, 0, ientr, GRID_REFERENCE, GRID_EDGE)
    CALL add_var( p_ext_atm_list, 'lsm_atm_e', p_ext_atm%lsm_atm_e, &
      &           GRID_UNSTRUCTURED_EDGE, ZAXIS_SURFACE, cf_desc, &
      &           grib2_desc, ldims=shape2d_e, lrestart=.FALSE.  )


    ! land fraction
    !
    ! fr_land      p_ext_atm%fr_land(nproma,nblks_c)
    cf_desc    = t_cf_var('land_area_fraction', '-', 'Fraction land')
    grib2_desc = t_grib2_var( 2, 0, 0, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_ext_atm_list, 'fr_land', p_ext_atm%fr_land, &
      &           GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc,&
      &            grib2_desc, ldims=shape2d_c )


    ! glacier fraction
    !
    ! fr_glac      p_ext_atm%fr_glac(nproma,nblks_c)
    cf_desc    = t_cf_var('glacier_area_fraction', '-', 'Fraction glacier')
    grib2_desc = t_grib2_var( 2, 0, 192, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_ext_atm_list, 'fr_glac', p_ext_atm%fr_glac, &
      &           GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d_c )


    ! maybe the next three (ice, fr_land_smt, fr_ice_smt)
    ! should be moved into corresponding if block

    ! sea Ice fraction
    !
    ! fr_ice       p_ext_atm%fr_ice(nproma,nblks_c)
    cf_desc    = t_cf_var('Sea_ice_fraction', '-', 'Sea ice fraction')
    grib2_desc = t_grib2_var( 10, 2, 0, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_ext_atm_list, 'fr_ice', p_ext_atm%fr_ice, &
      &           GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d_c )


    ! land fraction (smoothed)
    !
    ! fr_land_smt  p_ext_atm%fr_land_smt(nproma,nblks_c)
    cf_desc    = t_cf_var('land_area_fraction_(smoothed)', '-', &
      &                   'land area fraction (smoothed)')
    grib2_desc = t_grib2_var( 2, 0, 0, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_ext_atm_list, 'fr_land_smt', p_ext_atm%fr_land_smt, &
      &           GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d_c )



    ! glacier area fraction (smoothed)
    !
    ! fr_glac_smt  p_ext_atm%fr_glac_smt(nproma,nblks_c)
    cf_desc    = t_cf_var('glacier_area_fraction_(smoothed)', '-', &
      &                   'glacier area fraction (smoothed)')
    grib2_desc = t_grib2_var( 2, 0, 192, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_ext_atm_list, 'fr_glac_smt', p_ext_atm%fr_glac_smt, &
      &           GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d_c )



    ! sea Ice fraction (smoothed)
    !
    ! fr_ice_smt  p_ext_atm%fr_glac_smt(nproma,nblks_c)
    cf_desc    = t_cf_var('Sea_ice_fraction (smoothed)', '-', &
      &                   'Sea ice fraction (smoothed)')
    grib2_desc = t_grib2_var( 10, 2, 0, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_ext_atm_list, 'fr_ice_smt', p_ext_atm%fr_ice_smt, &
      &           GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d_c )


    ! geopotential (s)
    !
    ! fis          p_ext_atm%fis(nproma,nblks_c)
    cf_desc    = t_cf_var('Geopotential_(s)', 'm2 s-2', &
      &                   'Geopotential (s)')
    grib2_desc = t_grib2_var( 0, 3, 4, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_ext_atm_list, 'fis', p_ext_atm%fis, &
      &           GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d_c )


    ! external parameter for NWP forcing
    IF (iforcing == inwp) THEN

      ! Several IF-statements are necessary in order to allocate only
      ! those fields which are necessary for the chosen
      ! parameterizations.


      ! roughness length
      !
      ! z0           p_ext_atm%z0(nproma,nblks_c)
      cf_desc    = t_cf_var('roughtness_length', 'm', 'roughtness length')
      grib2_desc = t_grib2_var( 2, 0, 1, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'z0', p_ext_atm%z0, &
        &           GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d_c )


      ! fraction lake
      !
      ! fr_lake      p_ext_atm%fr_lake(nproma,nblks_c)
      cf_desc    = t_cf_var('fraction_lake', '-', 'fraction lake')
      grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'fr_lake', p_ext_atm%fr_lake, &
        &           GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d_c )


      ! lake depth
      !
      ! depth_lk     p_ext_atm%depth_lk(nproma,nblks_c)
      cf_desc    = t_cf_var('lake_depth', '-', 'lake depth')
      grib2_desc = t_grib2_var( 192, 228, 7, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'depth_lk', p_ext_atm%depth_lk, &
        &           GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d_c )



      !--------------------------------
      ! sub-gridscale orography
      !--------------------------------

      ! Standard deviation of sub-grid scale orography
      !
      ! sso_stdh     p_ext_atm%sso_stdh(nproma,nblks_c)
      cf_desc    = t_cf_var('standard_deviation_of_height', 'm',&
        &                   'Standard deviation of sub-grid scale orography')
      grib2_desc = t_grib2_var( 0, 3, 20, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'sso_stdh', p_ext_atm%sso_stdh, &
        &           GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d_c )



      ! Anisotropy of sub-gridscale orography
      !
      ! sso_gamma    p_ext_atm%sso_gamma(nproma,nblks_c)
      cf_desc    = t_cf_var('anisotropy_factor', '-',&
        &                   'Anisotropy of sub-gridscale orography')
      grib2_desc = t_grib2_var( 0, 3, 20, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'sso_gamma', p_ext_atm%sso_gamma, &
        &           GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d_c )



      ! Angle of sub-gridscale orography
      !
      ! sso_theta    p_ext_atm%sso_theta(nproma,nblks_c)
      cf_desc    = t_cf_var('angle_of_principal_axis', 'radians',&
        &                   'Angle of sub-gridscale orography')
      grib2_desc = t_grib2_var( 0, 3, 21, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'sso_theta', p_ext_atm%sso_theta, &
        &           GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d_c )



      ! Slope of sub-gridscale orography
      !
      ! sso_sigma    p_ext_atm%sso_sigma(nproma,nblks_c)
      cf_desc    = t_cf_var('slope_of_terrain', '-',&
        &                   'Slope of sub-gridscale orography')
      grib2_desc = t_grib2_var( 0, 3, 22, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'sso_sigma', p_ext_atm%sso_sigma, &
        &           GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d_c )





      !--------------------------------
      ! vegetation parameters
      !--------------------------------

      ! Plant covering degree in the vegetation phase
      !
      ! plcov_mx     p_ext_atm%plcov_mx(nproma,nblks_c)
      cf_desc    = t_cf_var('vegetation_area_fraction_vegetation_period', '-',&
        &                   'Plant covering degree in the vegetation phase')
      grib2_desc = t_grib2_var( 2, 0, 4, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'plcov_mx', p_ext_atm%plcov_mx, &
        &           GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d_c )



      ! Max Leaf area index
      !
      ! lai_mx       p_ext_atm%lai_mx(nproma,nblks_c)
      cf_desc    = t_cf_var('leaf_area_index_vegetation_period', '-',&
        &                   'Leaf Area Index Maximum')
      grib2_desc = t_grib2_var( 2, 0, 28, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'lai_mx', p_ext_atm%lai_mx, &
        &           GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d_c )



      ! root depth of vegetation
      !
      ! rootdp      p_ext_atm%rootdp(nproma,nblks_c)
      cf_desc    = t_cf_var('root_depth_of_vegetation', 'm',&
        &                   'root depth of vegetation')
      grib2_desc = t_grib2_var( 2, 0, 32, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'rootdp', p_ext_atm%rootdp, &
        &           GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d_c )



      ! evergreen forest
      !
      ! for_e        p_ext_atm%for_e(nproma,nblks_c)
      cf_desc    = t_cf_var('fraction_of_evergreen_forest_cover', '-',&
        &                   'Fraction of evergreen forest')
      grib2_desc = t_grib2_var( 2, 0, 29, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'for_e', p_ext_atm%for_e, &
        &           GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d_c )



      ! deciduous forest
      !
      ! for_d     p_ext_atm%for_d(nproma,nblks_c)
      cf_desc    = t_cf_var('fraction_of_deciduous_forest_cover', '-',&
        &                   'Fraction of deciduous forest')
      grib2_desc = t_grib2_var( 2, 0, 30, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'for_d', p_ext_atm%for_d, &
        &           GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d_c )



      ! urban area fraction
      !
      ! urban        p_ext_atm%urban(nproma,nblks_c)
      cf_desc    = t_cf_var('fraction_of_urban_areas', '-',&
        &                   'urban area fraction')
      grib2_desc = t_grib2_var( 2, 0, 30, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'urban', p_ext_atm%urban, &
        &           GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d_c )


      ! Minimal stomata resistence
      !
      ! rsmin        p_ext_atm%rsmin(nproma,nblks_c)
      cf_desc    = t_cf_var('RSMIN', 's m-1', 'Minimal stomata resistence')
      grib2_desc = t_grib2_var( 2, 0, 16, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'rsmin', p_ext_atm%rsmin, &
        &           GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d_c )



      ! NDVI yearly maximum
      !
      ! ndvi_max        p_ext_atm%ndvi_max(nproma,nblks_c)
      cf_desc    = t_cf_var('normalized_difference_vegetation_index', '-', &
        &                   'NDVI yearly maximum')
      grib2_desc = t_grib2_var( 2, 0, 31, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'ndvi_max', p_ext_atm%ndvi_max, &
        &           GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d_c )



      !--------------------------------
      ! soil parameters
      !--------------------------------

      ! soil type
      !
      ! soiltyp      p_ext_atm%soiltyp(nproma,nblks_c)
      cf_desc    = t_cf_var('soil_type', '-','soil type')
      grib2_desc = t_grib2_var( 2, 3, 0, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'soiltyp', p_ext_atm%soiltyp, &
        &           GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, &
        &           grib2_desc, ldims=shape2d_c, lrestart=.FALSE.   )


!DR      ! soil texture, keys 0-9
!DR      !
!DR      ! soiltyp_frac    p_ext_atm%soiltyp_frac(nproma,nblks_c,nsfc_subs)
!DR      cf_desc    = t_cf_var('soil_texture', '-','soil texture')
!DR      grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
!DR      CALL add_var( p_ext_atm_list, 'soiltyp_frac', p_ext_atm%soiltyp_frac, &
!DR        &           GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape3d_c )


      ! Climat. temperature
      !
      ! t_cl         p_ext_atm%t_cl(nproma,nblks_c)
      cf_desc    = t_cf_var('soil_temperature', 'K',                  &
        &                   'CRU near surface temperature climatology')
      grib2_desc = t_grib2_var( 2, 3, 18, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_ext_atm_list, 't_cl', p_ext_atm%t_cl, &
        &           GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d_c )


      ! longwave surface emissivity
      !
      ! emis_rad     p_ext_atm%emis_rad(nproma,nblks_c)
      cf_desc    = t_cf_var('emis_rad', '-', 'longwave surface emissivity')
      grib2_desc = t_grib2_var( 2, 3, 196, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'emis_rad', p_ext_atm%emis_rad, &
        &           GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, &
        &           grib2_desc, ldims=shape2d_c)


    ENDIF ! iforcing = nwp
  ENDIF ! iequations = 3


! #Hermann#
! Echam physics external parameter have to be allocated and initialized as well

  END SUBROUTINE new_ext_data_atm_list



  !-------------------------------------------------------------------------
  !
  !
  !>
  !! Allocation of atmospheric external data structure (time dependent)
  !!
  !! Allocation of atmospheric external data structure (time dependent  
  !! elements).
  !!
  !! Initialization of elements with zero.
  !!
  !! @par Revision History
  !! Initial release by Daniel Reinert (2011-05-03)
  !!
  SUBROUTINE new_ext_data_atm_td_list ( p_patch, p_ext_atm_td, &
    &                               p_ext_atm_td_list, listname)
!
    TYPE(t_patch), TARGET, INTENT(IN)     :: & !< current patch
      &  p_patch

    TYPE(t_external_atmos_td), INTENT(INOUT) :: & !< current external data structure
      &  p_ext_atm_td 

    TYPE(t_var_list) :: p_ext_atm_td_list  !< current external data list

    CHARACTER(len=*), INTENT(IN)      :: & !< list name
      &  listname

    TYPE(t_cf_var)    :: cf_desc
    TYPE(t_grib2_var) :: grib2_desc

    INTEGER :: nblks_c      !< number of cell blocks to allocate

    INTEGER :: shape3d_c(3)
    INTEGER :: shape4d_c(4)

    INTEGER :: ientr         !< "entropy" of horizontal slice
    INTEGER :: ntimes        !< number of time slices
    !--------------------------------------------------------------

    !determine size of arrays
    nblks_c = p_patch%nblks_c

    ientr  = 16   ! "entropy" of horizontal slice
    ntimes = 12   ! number of time slices

    ! predefined array shapes
    shape3d_c = (/ nproma, nblks_c, ntimes /)
    shape4d_c = (/ nproma, nlev_pres, nblks_c, nmonths /) 

    !
    ! Register a field list and apply default settings
    !
    CALL new_var_list( p_ext_atm_td_list, TRIM(listname) )
    CALL default_var_list_settings( p_ext_atm_td_list,         &
                                  & lrestart=.FALSE.,           &
                                  & restart_type=FILETYPE_NC2  )


    !--------------------------------
    ! radiation parameters
    !--------------------------------


    ! ozone on pressure levels
    ! ATTENTION: a GRIB2 number will go to 
    ! the ozone mass mixing ratio...
    !
    IF(irad_o3 == 3) THEN 
      ! o3       p_ext_atm_td%o3(nproma,nlev_pres,nblks_c,nmonths)
      cf_desc    = t_cf_var('O3', 'mole mole^-1',   &
        &                   'mole_fraction_of_ozone_in_air')
      grib2_desc = t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_ext_atm_td_list, 'O3', p_ext_atm_td%O3, &
        &           GRID_UNSTRUCTURED_CELL, ZAXIS_PRESSURE, cf_desc, &
        &           grib2_desc, ldims=shape4d_c )
    END IF


    ! Black carbon aerosol
    !
    ! aer_bc       p_ext_atm%aer_bc(nproma,nblks_c,ntimes)
    cf_desc    = t_cf_var('aerosol optical thickness of black carbon', '-',   &
      &                   'atmosphere_absorption_optical_thickness_due_to_' //&
      &                   'black_carbon_ambient_aerosol')
    grib2_desc = t_grib2_var( 0, 13, 195, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_ext_atm_td_list, 'aer_bc', p_ext_atm_td%aer_bc, &
      &           GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc,  &
      &            grib2_desc, ldims=shape3d_c )


    ! Dust aerosol
    !
    ! aer_dust     p_ext_atm%aer_dust(nproma,nblks_c,ntimes)
    cf_desc    = t_cf_var('aot_dust', '-', &
      &                   'atmosphere absorption optical thickness due '//  &
      &                   'to dust ambient aerosol')
    grib2_desc = t_grib2_var( 0, 13, 193, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_ext_atm_td_list, 'aer_dust', p_ext_atm_td%aer_dust, &
      &           GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, &
      &           ldims=shape3d_c )


    ! Organic aerosol
    !
    ! aer_org      p_ext_atm%aer_org(nproma,nblks_c,ntimes)
    cf_desc    = t_cf_var('aot_org', '-', &
      &                   'atmosphere absorption optical thickness due '//  &
      &                   'to particulate organic matter ambient aerosol')
    grib2_desc = t_grib2_var( 0, 13, 194, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_ext_atm_td_list, 'aer_org', p_ext_atm_td%aer_org, &
      &           GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc,&
      &           ldims=shape3d_c )


    ! Sulfate aerosol
    !
    ! aer_so4      p_ext_atm%aer_so4(nproma,nblks_c,ntimes)
    cf_desc    = t_cf_var('aot_so4', '-', &
      &                   'atmosphere absorption optical thickness due '//  &
      &                   'to sulfate_ambient_aerosol')
    grib2_desc = t_grib2_var( 0, 13, 192, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_ext_atm_td_list, 'aer_so4', p_ext_atm_td%aer_so4, &
      &           GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc,&
      &           ldims=shape3d_c )


    ! Seasalt aerosol
    !
    ! aer_ss       p_ext_atm%aer_ss(nproma,nblks_c,ntimes)
    cf_desc    = t_cf_var('aot_ss', '-', &
      &                   'atmosphere absorption optical thickness due '//  &
      &                   'to seasalt_ambient_aerosol')
    grib2_desc = t_grib2_var( 0, 13, 196, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_ext_atm_td_list, 'aer_ss', p_ext_atm_td%aer_ss, &
      &           GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc,&
      &           ldims=shape3d_c )


    !--------------------------------
    ! vegetation parameters
    !--------------------------------

    ! monthly mean normalized difference vegetation index
    !
    ! ndvi         p_ext_atm%ndvi(nproma,nblks_c,ntimes)
    cf_desc    = t_cf_var('normalized_difference_vegetation_index', '-', &
      &                   'monthly mean NDVI')
    grib2_desc = t_grib2_var( 2, 0, 217, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_ext_atm_td_list, 'ndvi', p_ext_atm_td%ndvi, &
      &           GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc,&
      &           ldims=shape3d_c )


    ! (monthly) proportion of actual value/maximum NDVI
    !
    ! ndvi_mrat     p_ext_atm%ndvi_mrat(nproma,nblks_c,ntimes)
    cf_desc    = t_cf_var('normalized_difference_vegetation_index', '-', &
      &                   '(monthly) proportion of actual value/maximum ' // &
      &                   'normalized differential vegetation index')
    grib2_desc = t_grib2_var( 2, 0, 192, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_ext_atm_td_list, 'ndvi_mrat', p_ext_atm_td%ndvi_mrat, &
      &           GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc,&
      &           ldims=shape3d_c)


  END SUBROUTINE new_ext_data_atm_td_list



  !-------------------------------------------------------------------------
  !
  !
  !>
  !! Allocation of oceanic external data structure
  !!
  !! Allocation of oceanic external data structure (constant in time 
  !! elements).
  !!
  !! Initialization of elements with zero.
  !!
  !! @par Revision History
  !! Initial release by Daniel Reinert (2011-06-24)
  !!
  SUBROUTINE new_ext_data_oce_list ( p_patch, p_ext_oce, p_ext_oce_list, &
    &                                listname)
!
    TYPE(t_patch), TARGET, INTENT(IN)   :: & !< current patch
      &  p_patch

    TYPE(t_external_ocean), INTENT(INOUT) :: & !< current external data structure
      &  p_ext_oce 

    TYPE(t_var_list) :: p_ext_oce_list !< current external data list

    CHARACTER(len=*), INTENT(IN)  :: & !< list name
      &  listname

    TYPE(t_cf_var)    :: cf_desc
    TYPE(t_grib2_var) :: grib2_desc

    INTEGER :: nblks_c, &    !< number of cell blocks to allocate
      &        nblks_e       !< number of edge blocks to allocate

    INTEGER :: shape2d_c(2), shape2d_e(2)

    INTEGER :: ientr         !< "entropy" of horizontal slice
    !--------------------------------------------------------------

    !determine size of arrays
    nblks_c = p_patch%nblks_c
    nblks_e = p_patch%nblks_e


    ientr = 16   ! "entropy" of horizontal slice

    ! predefined array shapes
    shape2d_c = (/ nproma, nblks_c /)
    shape2d_e = (/ nproma, nblks_e /)

    !
    ! Register a field list and apply default settings
    !
    CALL new_var_list( p_ext_oce_list, TRIM(listname) )
    CALL default_var_list_settings( p_ext_oce_list,            &
                                  & lrestart=.TRUE.,           &
                                  & restart_type=FILETYPE_NC2  )



    ! bathymetric height at cell center
    !
    ! bathymetry_c  p_ext_oce%bathymetry_c(nproma,nblks_c)
    cf_desc    = t_cf_var('Model bathymetry at cell center', 'm', &
      &                   'Model bathymetry')
    grib2_desc = t_grib2_var( 192, 140, 219, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_ext_oce_list, 'bathymetry_c', p_ext_oce%bathymetry_c,      &
      &           GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d_c )


    ! bathymetric height at cell edge
    !
    ! bathymetry_e  p_ext_oce%bathymetry_e(nproma,nblks_e)
    cf_desc    = t_cf_var('Model bathymetry at cell edge', 'm', &
      &                   'Model bathymetry')
    grib2_desc = t_grib2_var( 192, 140, 219, ientr, GRID_REFERENCE, GRID_EDGE)
    CALL add_var( p_ext_oce_list, 'bathymetry_e', p_ext_oce%bathymetry_e,      &
      &           GRID_UNSTRUCTURED_EDGE, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d_e )

  END SUBROUTINE new_ext_data_oce_list



  !-------------------------------------------------------------------------
  !>
  !! Destruct external data data structure and lists
  !!
  !! Destruct external data data structure and lists
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2011-05-04)
  !!
  SUBROUTINE destruct_ext_data

    INTEGER :: jg
    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
      routine = 'mo_ext_data:destruct_ext_data'
!-------------------------------------------------------------------------

    CALL message (TRIM(routine), 'Destruction of data structure for' // &
      &                          'external data started')

    DO jg = 1,n_dom

      ! Delete list of constant in time atmospheric elements
      CALL delete_var_list( ext_data(jg)%atm_list )

      IF (iforcing==inwp) THEN
      ! Delete list of time-dependent atmospheric elements
      CALL delete_var_list( ext_data(jg)%atm_td_list )
      END IF

      IF (iforcing==ihs_ocean) THEN
        ! Delete list of constant in time oceanic elements
        CALL delete_var_list( ext_data(jg)%oce_list )
      ENDIF

      ! Delete list of time-dependent oceanic elements
      ! ### to be added if necessary ###

    ENDDO

    CALL message (TRIM(routine), 'Destruction of data structure for' // &
      &                          'external data finished')

  END SUBROUTINE destruct_ext_data




 !-------------------------------------------------------------------------
 !-------------------------------------------------------------------------

  SUBROUTINE inquire_external_files(p_patch)

    !-------------------------------------------------------
    !
    ! open netcdf files and investigate the data structure  
    ! of the external parameters
    !
    !-------------------------------------------------------

    TYPE(t_patch), INTENT(IN)            :: p_patch(:)
    INTEGER :: no_cells
    INTEGER :: ncid, dimid
    INTEGER :: jg

    LOGICAL :: l_exist

    CHARACTER(len=max_char_length), PARAMETER :: &
      routine = 'mo_ext_data: inquire_external_files'

    CHARACTER(filename_max) :: ozone_file  !< file name for reading in

    ! default values for nlev_pres and nmonths
    nlev_pres = 1
    nmonths   = 1


    DO jg= 1,n_dom

       IF(irad_o3 == 3) THEN

       IF(p_pe == p_io) THEN

        WRITE(ozone_file,'(a,i2.2,a)') 'o3_icon_DOM',jg,'.nc'

        INQUIRE (FILE=ozone_file, EXIST=l_exist)
        IF (.NOT.l_exist) THEN
          CALL finish(TRIM(routine),'ozone file of domain is not found.')
        ENDIF

        !
        ! open file: do I have to enhance the number of ncid?
        !
        CALL nf(nf_open(TRIM(ozone_file), NF_NOWRITE, ncid))

        ! get number of cells in triangles and hexagons
        !

        !triangles
        IF (p_patch(jg)%cell_type == 3) THEN ! triangular grid
           CALL nf(nf_inq_dimid (ncid, 'cell', dimid))
           CALL nf(nf_inq_dimlen(ncid, dimid, no_cells))
       ENDIF
       
       !hexagons
        IF (p_patch(jg)%cell_type == 6) THEN ! hexagonal grid
           CALL nf(nf_inq_dimid (ncid, 'vertex', dimid))
           CALL nf(nf_inq_dimlen(ncid, dimid, no_cells))
        ENDIF

        !
        ! check the number of cells and verts
        !
        IF(p_patch(jg)%n_patch_cells_g /= no_cells) THEN
          CALL finish(TRIM(ROUTINE),&
          & 'Number of patch cells and cells in ozone file do not match.')
        ENDIF
          
      ! check the vertical structure
        CALL nf(nf_inq_dimid (ncid, 'plev', dimid))
        CALL nf(nf_inq_dimlen(ncid, dimid, nlev_pres))
        
        CALL message(TRIM(ROUTINE),message_text)
        WRITE(message_text,'(A,I4)')  &
           & 'Number of pressure levels in ozone file = ', &
           & nlev_pres

        ! check the time structure
        CALL nf(nf_inq_dimid (ncid, 'time', dimid))
        CALL nf(nf_inq_dimlen(ncid, dimid, nmonths))
        
        CALL message(TRIM(ROUTINE),message_text)
        WRITE(message_text,'(A,I4)')  &
           & 'Number of months in ozone file = ', &
           & nmonths

     ENDIF ! pe

  ENDIF !o3

  IF(p_pe == p_io) CALL nf(nf_close(ncid))

  ENDDO ! ndom

END SUBROUTINE inquire_external_files


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

    CHARACTER(filename_max) :: extpar_file !< file name for reading in
    CHARACTER(filename_max) :: ozone_file  !< file name for reading in

    LOGICAL :: l_exist
    INTEGER :: jg
    INTEGER :: i_lev, no_cells, no_verts
    INTEGER :: ncid, dimid


!-------------------------------------------------------------------------

    IF(itopo == 1 ) THEN
    DO jg = 1,n_dom

      i_lev = p_patch(jg)%level

      IF(p_pe == p_io) THEN
        !
        ! generate file name
        !
!         WRITE(extpar_file,'(a,i1,a,i1,a,i1,a)') &
!           & 'extpar_R',nroot,'B0',start_lev,'_DOM0',jg,'.nc'
        WRITE(extpar_file,'(a,a)') &
          & 'extpar_',TRIM(p_patch(jg)%grid_filename)
        CALL message("read_ext_data_atm, extpar_file=",extpar_file)
        
        INQUIRE (FILE=extpar_file, EXIST=l_exist)
        IF (.NOT.l_exist) THEN
          CALL finish(TRIM(routine),'external data file is not found.')
        ENDIF

        !
        ! open file
        !
        CALL nf(nf_open(TRIM(extpar_file), NF_NOWRITE, ncid))

        !
        ! get number of cells and vertices
        !
        CALL nf(nf_inq_dimid(ncid, 'cell', dimid))
        IF (p_patch(jg)%cell_type == 3) THEN ! triangular grid
          CALL nf(nf_inq_dimlen(ncid, dimid, no_cells))
        ELSEIF (p_patch(jg)%cell_type == 6) THEN ! hexagonal grid
          CALL nf(nf_inq_dimlen(ncid, dimid, no_verts))
        ENDIF

        CALL nf(nf_inq_dimid(ncid, 'vertex', dimid))
        IF (p_patch(jg)%cell_type == 3) THEN ! triangular grid
          CALL nf(nf_inq_dimlen(ncid, dimid, no_verts))
        ELSEIF (p_patch(jg)%cell_type == 6) THEN ! hexagonal grid
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


      !
      ! topography
      !
      IF (p_patch(jg)%cell_type == 3) THEN     ! triangular grid

        ! triangle center
        CALL read_netcdf_data (ncid, 'topography_c', p_patch(jg)%n_patch_cells_g,       &
          &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
          &                     ext_data(jg)%atm%topography_c)

        ! triangle vertex
        CALL read_netcdf_data (ncid, 'topography_v', p_patch(jg)%n_patch_verts_g,       &
          &                     p_patch(jg)%n_patch_verts, p_patch(jg)%verts%glb_index, &
          &                     ext_data(jg)%atm%topography_v)

      ELSEIF (p_patch(jg)%cell_type == 6) THEN ! hexagonal grid

        ! As extpar "knows" only the triangular grid, cells and vertices need to be switched here
        ! triangle center
        CALL read_netcdf_data (ncid, 'topography_c', p_patch(jg)%n_patch_verts_g,       &
          &                     p_patch(jg)%n_patch_verts, p_patch(jg)%verts%glb_index, &
          &                     ext_data(jg)%atm%topography_v)

        ! triangle vertex
        CALL read_netcdf_data (ncid, 'topography_v', p_patch(jg)%n_patch_cells_g,       &
          &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
          &                     ext_data(jg)%atm%topography_c)

      ENDIF


      !
      ! other external parameters on triangular grid
      !
      IF (p_patch(jg)%cell_type == 3) THEN     ! triangular grid

        CALL read_netcdf_data (ncid, 'FR_LAND', p_patch(jg)%n_patch_cells_g,              &
          &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index,   &
          &                     ext_data(jg)%atm%fr_land)

        CALL read_netcdf_data (ncid, 'ICE', p_patch(jg)%n_patch_cells_g,                &
          &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
          &                     ext_data(jg)%atm%fr_ice)


        IF (iforcing == inwp) THEN

          CALL read_netcdf_data (ncid, 'PLCOV_MX', p_patch(jg)%n_patch_cells_g,           &
            &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
            &                     ext_data(jg)%atm%plcov_mx)

          CALL read_netcdf_data (ncid, 'LAI_MX', p_patch(jg)%n_patch_cells_g,             &
            &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
            &                     ext_data(jg)%atm%lai_mx)

          CALL read_netcdf_data (ncid, 'ROOTDP', p_patch(jg)%n_patch_cells_g,             &
            &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
            &                     ext_data(jg)%atm%rootdp)

          CALL read_netcdf_data (ncid, 'RSMIN', p_patch(jg)%n_patch_cells_g,             &
            &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
            &                     ext_data(jg)%atm%rsmin)

          CALL read_netcdf_data (ncid, 'FOR_D', p_patch(jg)%n_patch_cells_g,              &
            &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
            &                     ext_data(jg)%atm%for_d)

          CALL read_netcdf_data (ncid, 'FOR_E', p_patch(jg)%n_patch_cells_g,              &
            &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
            &                     ext_data(jg)%atm%for_e)

          CALL read_netcdf_data (ncid, 'URBAN', p_patch(jg)%n_patch_cells_g,              &
            &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
            &                     ext_data(jg)%atm%urban)

          CALL read_netcdf_data (ncid, 'Z0', p_patch(jg)%n_patch_cells_g,                 &
            &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
            &                     ext_data(jg)%atm%z0)

          CALL read_netcdf_data (ncid, 'NDVI_MAX', p_patch(jg)%n_patch_cells_g,           &
            &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
            &                     ext_data(jg)%atm%ndvi_max)

          CALL read_netcdf_data (ncid, 'SOILTYP', p_patch(jg)%n_patch_cells_g,        &
            &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
            &                     ext_data(jg)%atm%soiltyp)

          CALL read_netcdf_data (ncid, 'EMIS_RAD', p_patch(jg)%n_patch_cells_g,           &
            &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
            &                     ext_data(jg)%atm%emis_rad)

          CALL read_netcdf_data (ncid, 'T_CL', p_patch(jg)%n_patch_cells_g,               &
            &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
            &                     ext_data(jg)%atm%t_cl)

          CALL read_netcdf_data (ncid, 'SSO_STDH', p_patch(jg)%n_patch_cells_g,           &
            &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
            &                     ext_data(jg)%atm%sso_stdh)

          CALL read_netcdf_data (ncid, 'SSO_THETA', p_patch(jg)%n_patch_cells_g,          &
            &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
            &                     ext_data(jg)%atm%sso_theta)

          CALL read_netcdf_data (ncid, 'SSO_GAMMA', p_patch(jg)%n_patch_cells_g,          &
            &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
            &                     ext_data(jg)%atm%sso_gamma)

          CALL read_netcdf_data (ncid, 'SSO_SIGMA', p_patch(jg)%n_patch_cells_g,          &
            &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
            &                     ext_data(jg)%atm%sso_sigma)

          CALL read_netcdf_data (ncid, 'FR_LAKE', p_patch(jg)%n_patch_cells_g,            &
            &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
            &                     ext_data(jg)%atm%fr_lake)

          CALL read_netcdf_data (ncid, 'DEPTH_LK', p_patch(jg)%n_patch_cells_g,           &
            &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
            &                     ext_data(jg)%atm%depth_lk)
        ENDIF ! (iforcing == inwp)
        
      ELSEIF (p_patch(jg)%cell_type == 6) THEN ! hexagonal grid

        CALL finish(TRIM(ROUTINE),&
        & 'Hexagonal grid is not supported, yet.')

      ENDIF

      !
      ! close file
      !
      IF(p_pe == p_io) CALL nf(nf_close(ncid))

    ENDDO

    ENDIF ! itopo


      !-------------------------------------------------------
      ! Read ozone
      !-------------------------------------------------------

    IF(irad_o3 == 3) THEN
      DO jg = 1,n_dom
        IF(p_pe == p_io) THEN


          !
          ! open file
          !
          WRITE(ozone_file,'(a,I2.2,a)') 'o3_icon_DOM',jg,'.nc'
          CALL nf(nf_open(TRIM(ozone_file), NF_NOWRITE, ncid))

          IF (p_patch(jg)%cell_type == 3) THEN     ! triangular grid
            CALL read_netcdf_data (ncid, 'O3', & ! &
              &                    p_patch(jg)%n_patch_cells_g,  &
              &                    p_patch(jg)%n_patch_cells,    &
              &                    p_patch(jg)%cells%glb_index,  & 
              &                    nlev_pres,  nmonths,          &
              &                    ext_data(jg)%atm_td%O3)
          ELSEIF (p_patch(jg)%cell_type == 6) THEN ! hexagonal grid
            CALL read_netcdf_data (ncid, 'O3', & 
              &                    p_patch(jg)%n_patch_verts_g,  &
              &                    p_patch(jg)%n_patch_verts,    & 
              &                    p_patch(jg)%verts%glb_index,  &
              &                    nlev_pres, nmonths,           &
              &                    ext_data(jg)%atm_td%O3)
          ENDIF
        ENDIF ! pe

        !
        ! close file
        !
        IF(p_pe == p_io) CALL nf(nf_close(ncid))

      ENDDO ! ndom
    ENDIF ! irad_o3

!    write(0,*)'try to give a number of ozone'
!    write(0,*)'any value jan', ext_data(1)%atm_td%O3(1,1,1,1)
!    write(0,*)'maxval dec',MAXVAL( ext_data(1)%atm_td%O3(:,nlev_pres,:,nmonths))
 

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
  SUBROUTINE read_netcdf_2d (ncid, varname, glb_arr_len, loc_arr_len, glb_index, var_out)

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

  END SUBROUTINE read_netcdf_2d



  !-------------------------------------------------------------------------
  !>
  !! Read dataset from netcdf file
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2010-07-14)
  !! Adapted for parallel runs by Rainer Johanni (2010-12-07)
  !!
  SUBROUTINE read_netcdf_2d_int (ncid, varname, glb_arr_len, loc_arr_len, glb_index, var_out)

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

  END SUBROUTINE read_netcdf_2d_int

 !-------------------------------------------------------------------------
  !>
  !! Read 3D (inlcuding height) dataset from netcdf file
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2010-07-14)
  !! Adapted for parallel runs by Rainer Johanni (2010-12-07)
  !! Adapted for 3 D by Guenther Zaengl, DWD (2011-07-11)
  !!
  SUBROUTINE read_netcdf_3d (ncid, varname, glb_arr_len, &
       &                     loc_arr_len, glb_index, &
       &                     nlevs,      var_out)

    CHARACTER(len=*), INTENT(IN)  ::  &  !< Var name of field to be read
      &  varname

    INTEGER, INTENT(IN) :: ncid          !< id of netcdf file
    INTEGER, INTENT(IN) :: nlevs         !< vertical levels of netcdf file
    INTEGER, INTENT(IN) :: glb_arr_len   !< length of 1D field (global)
    INTEGER, INTENT(IN) :: loc_arr_len   !< length of 1D field (local)
    INTEGER, INTENT(IN) :: glb_index(:)  !< Index mapping local to global

    REAL(wp), INTENT(INOUT) :: &         !< output field
      &  var_out(:,:,:)

    INTEGER :: varid, mpi_comm, j, jl, jb, jk
    REAL(wp):: z_dummy_array(glb_arr_len,nlevs)!< local dummy array
  !-------------------------------------------------------------------------

    ! Get var ID
    IF(p_pe==p_io) CALL nf(nf_inq_varid(ncid, TRIM(varname), varid))

    IF(p_test_run) THEN
      mpi_comm = p_comm_work_test
    ELSE
      mpi_comm = p_comm_work
    ENDIF

    ! I/O PE reads and broadcasts data

    IF(p_pe==p_io) CALL nf(nf_get_var_double(ncid, varid, z_dummy_array(:,:)))
    CALL p_bcast(z_dummy_array, p_io, mpi_comm)

    var_out(:,:,:) = 0._wp

    ! Set var_out from global data
     DO jk = 1, nlevs
       DO j = 1, loc_arr_len

         jb = blk_no(j) ! Block index in distributed patch
         jl = idx_no(j) ! Line  index in distributed patch
               
         var_out(jl,jk,jb) = z_dummy_array(glb_index(j),jk)

       ENDDO
     ENDDO

  END SUBROUTINE read_netcdf_3d


 !-------------------------------------------------------------------------
  !>
  !! Read 4D (inlcuding height and time) dataset from netcdf file
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2010-07-14)
  !! Adapted for parallel runs by Rainer Johanni (2010-12-07)
  !! Adapted for 4 D by Kristina Froehlich, MPI-M (2011-06-16)
  !!
  SUBROUTINE read_netcdf_4d (ncid, varname, glb_arr_len, &
       &                     loc_arr_len, glb_index, &
       &                     nlevs, ntime,      var_out)

    CHARACTER(len=*), INTENT(IN)  ::  &  !< Var name of field to be read
      &  varname

    INTEGER, INTENT(IN) :: ncid          !< id of netcdf file
    INTEGER, INTENT(IN) :: nlevs         !< vertical levels of netcdf file
    INTEGER, INTENT(IN) :: ntime         !< time levels of netcdf file
    INTEGER, INTENT(IN) :: glb_arr_len   !< length of 1D field (global)
    INTEGER, INTENT(IN) :: loc_arr_len   !< length of 1D field (local)
    INTEGER, INTENT(IN) :: glb_index(:)  !< Index mapping local to global

    REAL(wp), INTENT(INOUT) :: &         !< output field
      &  var_out(:,:,:,:)

    INTEGER :: varid, mpi_comm, j, jl, jb, jk, jt
    REAL(wp):: z_dummy_array(glb_arr_len,nlevs,ntime)!< local dummy array
  !-------------------------------------------------------------------------

    ! Get var ID
    IF(p_pe==p_io) CALL nf(nf_inq_varid(ncid, TRIM(varname), varid))

    IF(p_test_run) THEN
      mpi_comm = p_comm_work_test
    ELSE
      mpi_comm = p_comm_work
    ENDIF

    ! I/O PE reads and broadcasts data

    IF(p_pe==p_io) CALL nf(nf_get_var_double(ncid, varid, z_dummy_array(:,:,:)))
    CALL p_bcast(z_dummy_array, p_io, mpi_comm)

    var_out(:,:,:,:) = 0._wp

    ! Set var_out from global data
    DO jt = 1, ntime
        DO jk = 1, nlevs
           DO j = 1, loc_arr_len

             jb = blk_no(j) ! Block index in distributed patch
             jl = idx_no(j) ! Line  index in distributed patch
               
             var_out(jl,jk,jb,jt) = z_dummy_array(glb_index(j),jk,jt)

          ENDDO
       ENDDO
    ENDDO

  END SUBROUTINE read_netcdf_4d


  !-------------------------------------------------------------------------

  SUBROUTINE nf(status)

    INTEGER, INTENT(in) :: status

    IF (status /= nf_noerr) THEN
      CALL finish('mo_ext_data netCDF error', nf_strerror(status))
    ENDIF

  END SUBROUTINE nf

  !-------------------------------------------------------------------------

   SUBROUTINE smooth_topography (p_patch, p_int, topography_c, topography_v)

    TYPE(t_patch), INTENT(IN)       :: p_patch
    TYPE(t_int_state), INTENT(IN)   :: p_int
    REAL(wp), INTENT(INOUT)         :: topography_c(:,:), topography_v(:,:)

    ! local variables
    INTEGER  :: jg, jb, jc, iter, il
    INTEGER  :: i_startblk, nblks_c, i_startidx, i_endidx
    REAL(wp) :: z_topo(nproma,1,p_patch%nblks_c),z_nabla4_topo(nproma,1,p_patch%nblks_c), &
      &         z_topo_old(nproma,1,p_patch%nblks_c)
    REAL(wp) :: z_topo_v(nproma,1,p_patch%nblks_v)
    REAL(wp) :: zmaxtop,zmintop,z_topo_new


    jg = p_patch%id

    z_topo_v(:,1,:) = topography_v(:,:)
    z_nabla4_topo(:,1,:) = 0._wp

    i_startblk = p_patch%cells%start_blk(2,1)
    nblks_c    = p_patch%nblks_c

!    write(0,*) 'n_iter_smooth_topo=',n_iter_smooth_topo
!    write(0,*) 'fac_smooth_topo=',fac_smooth_topo

    DO iter = 1, n_iter_smooth_topo
      z_topo(:,1,:)   = topography_c(:,:)
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
          DO il=2,p_patch%cell_type

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
          IF (( zmaxtop < z_topo_new  ) .AND. ( z_nabla4_topo(jc,1,jb) < 0.0_wp )) THEN
            topography_c(jc,jb) = MAX(z_topo_old(jc,1,jb),zmaxtop)
            CYCLE
          ENDIF
          !If it became a local minimum in the new field with regard to old neighbors, avoid it:
          IF (( zmintop > z_topo_new  ) .AND. ( z_nabla4_topo(jc,1,jb) > 0.0_wp )) THEN
            topography_c(jc,jb) = MIN(z_topo_old(jc,1,jb),zmintop)
            CYCLE
          ENDIF

          topography_c(jc,jb)=z_topo_new

        ENDDO

      ENDDO

      z_topo(:,1,:)   = topography_c(:,:)
      CALL sync_patch_array(SYNC_C, p_patch, z_topo)
      topography_c(:,:)=z_topo(:,1,:)

      DO jb = i_startblk,nblks_c

        CALL get_indices_c(p_patch, jb, i_startblk, nblks_c, &
                           i_startidx, i_endidx, 2)
        DO jc = i_startidx, i_endidx

          !Limiter to avoid amplification of local extrema
          zmaxtop = z_topo(p_patch%cells%neighbor_idx(jc,jb,1),1, &
              &          p_patch%cells%neighbor_blk(jc,jb,1))
          zmintop = zmaxtop
          DO il=2,p_patch%cell_type

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

           !If it became a local maximum in the new field, avoid it
           IF ( ( zmaxtop < z_topo(jc,1,jb) ) .AND. ( z_nabla4_topo(jc,1,jb) < 0.0_wp ) ) THEN
             topography_c(jc,jb) = MAX(z_topo_old(jc,1,jb),zmaxtop)
           !If it became a local minimum in the new field, avoid it
           ELSEIF (( zmintop > z_topo(jc,1,jb) ) .AND. ( z_nabla4_topo(jc,1,jb) > 0.0_wp )) THEN
             topography_c(jc,jb) = MIN(z_topo_old(jc,1,jb),zmintop)
           ENDIF

        ENDDO

      ENDDO

      z_topo(:,1,:)   = topography_c(:,:)
      CALL sync_patch_array(SYNC_C, p_patch, z_topo)
      topography_c(:,:)=z_topo(:,1,:)

    ENDDO !iter


    ! Interpolate smooth topography from cells to vertices
    z_topo(:,1,:)   = topography_c(:,:)
    CALL cells2verts_scalar(z_topo,p_patch,p_int%cells_aw_verts,z_topo_v,1,1)

    CALL sync_patch_array(SYNC_V,p_patch,z_topo_v)

    topography_v(:,:) = z_topo_v(:,1,:)


  END SUBROUTINE smooth_topography


END MODULE mo_ext_data

