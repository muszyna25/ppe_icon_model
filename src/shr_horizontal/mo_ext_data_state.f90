!>
!! Allocation/deallocation and reading of external datasets
!!
!! This module contains routines for setting up the external data state for 
!! atmosphere and ocean.
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
!! Modification by Daniel Reinert, DWD (2012-02-23)
!! - Routine smooth_topography moved to a new module named mo_smooth_topo
!! Modification by Daniel Reinert, DWD (2012-03-22)
!! - Type declaration moved to new module mo_ext_data_types
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

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_ext_data_state

  USE mo_kind,               ONLY: wp
  USE mo_io_units,           ONLY: filename_max
  USE mo_parallel_config,    ONLY: nproma
  USE mo_impl_constants,     ONLY: inwp, iecham, ildf_echam, io3_clim, io3_ape, &
    &                              ihs_ocean, ihs_atm_temp, ihs_atm_theta, inh_atmosphere, &
    &                              max_char_length, min_rlcell_int,                        &
    &                              VINTP_METHOD_LIN
  USE mo_physical_constants, ONLY: ppmv2gg, zemiss_def
  USE mo_run_config,         ONLY: iforcing
  USE mo_ocean_nml,          ONLY: iforc_oce, iforc_type, iforc_len
  USE mo_impl_constants_grf, ONLY: grf_bdywidth_c
  USE mo_lnd_nwp_config,     ONLY: ntiles_total, ntiles_lnd, ntiles_water, lsnowtile, frlnd_thrhld, &
                                   frlndtile_thrhld, frlake_thrhld, frsea_thrhld, isub_water,       &
                                   sstice_mode
  USE mo_extpar_config,      ONLY: itopo, l_emiss, extpar_filename, generate_filename
  USE mo_time_config,        ONLY: time_config
  USE mo_dynamics_config,    ONLY: iequations
  USE mo_radiation_config,   ONLY: irad_o3, irad_aero
  USE mo_echam_phy_config,   ONLY: echam_phy_config
  USE mo_smooth_topo,        ONLY: smooth_topography
  USE mo_model_domain,       ONLY: t_patch
  USE mo_exception,          ONLY: message, message_text, finish
  USE mo_grid_config,        ONLY: n_dom, nroot, dynamics_grid_filename
  USE mo_intp_data_strc,     ONLY: t_int_state
  USE mo_loopindices,        ONLY: get_indices_c
  USE mo_mpi,                ONLY: my_process_is_stdio, p_io, p_bcast, &
    &                              p_comm_work_test, p_comm_work
  USE mo_sync,               ONLY: global_sum_array
  USE mo_parallel_config,    ONLY: p_test_run
  USE mo_linked_list,        ONLY: t_var_list
  USE mo_ext_data_types,     ONLY: t_external_data, t_external_atmos,    &
    &                              t_external_atmos_td, t_external_ocean
  USE mo_var_list,           ONLY: default_var_list_settings, &
    &                              add_var, add_ref,          &
    &                              new_var_list,              &
    &                              delete_var_list, create_vert_interp_metadata
  USE mo_master_nml,         ONLY: model_base_dir
  USE mo_cf_convention,      ONLY: t_cf_var
  USE mo_grib2,              ONLY: t_grib2_var
  USE mo_util_netcdf,        ONLY: read_netcdf_data, read_netcdf_lu, nf
  USE mo_util_string,        ONLY: t_keyword_list,  &
    &                              associate_keyword, with_keywords
  USE mo_phyparam_soil,      ONLY: c_lnd, c_soil
  USE mo_datetime,           ONLY: t_datetime, month2hour
  USE mo_cdi_constants,      ONLY: GRID_UNSTRUCTURED_CELL, GRID_UNSTRUCTURED_EDGE, &
    &                              GRID_UNSTRUCTURED_VERT, GRID_REFERENCE,         &
    &                              GRID_CELL, GRID_EDGE, GRID_VERTEX, ZA_SURFACE,  &
    &                              ZA_HYBRID, ZA_PRESSURE, DATATYPE_FLT32,         &
    &                              DATATYPE_PACK16, FILETYPE_NC2, TSTEP_CONSTANT

  IMPLICIT NONE

  ! required for reading external data
  INCLUDE 'netcdf.inc'

  PRIVATE

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

  INTEGER::  nlev_o3, nmonths

  CHARACTER(len=6)  :: levelname
  CHARACTER(len=4)  :: zlevelname
  CHARACTER(len=6)  :: cellname
  CHARACTER(len=5)  :: o3name
  CHARACTER(len=20) :: o3unit

  ! Number of landcover classes provided by external parameter data
  ! Needs to be changed into a variable if landcover classifications 
  ! with a different number of classes become available
  INTEGER, PARAMETER :: num_lcc = 23, n_param_lcc = 7

  INTEGER, ALLOCATABLE :: nclass_lu(:)  !< number of landuse classes
                                        !< dim: n_dom
  INTEGER, ALLOCATABLE :: nmonths_ext(:)!< number of months in external data file
                                        !< dim: n_dom

  PUBLIC :: ext_data
  PUBLIC :: nmonths
  PUBLIC :: nlev_o3

  PUBLIC :: init_ext_data  
  PUBLIC :: init_ext_data_oce
  PUBLIC :: init_index_lists
  PUBLIC :: destruct_ext_data


  TYPE(t_external_data),TARGET, ALLOCATABLE :: &
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
    TYPE(t_int_state), INTENT(IN)        :: p_int_state(:)
    TYPE(t_external_data), INTENT(INOUT) :: ext_data(:)


    INTEGER :: jg

    CHARACTER(len=max_char_length), PARAMETER :: &
      routine = 'mo_ext_data:init_ext_data'

    !-------------------------------------------------------------------------
    CALL message (TRIM(routine), 'Start')

    !-------------------------------------------------------------------------
    !  1.  inquire external files for their data structure
    !-------------------------------------------------------------------------

    ALLOCATE(nclass_lu(n_dom))
    ! Set default value for nclass_lu. Will be overwritten, if external data 
    ! are read from file
    nclass_lu(1:n_dom) = 1

    ALLOCATE(nmonths_ext(n_dom))
    ! Set default value for nmonths_ext. Will be overwritten, if external data 
    ! are read from file
    nmonths_ext(1:n_dom) = 1

    CALL inquire_external_files(p_patch)

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

    !-------------------------------------------------------------------------
    ! Now for atmosphere only:
    !-------------------------------------------------------------------------

    IF (iequations/=ihs_ocean ) THEN

!!$    !-------------------------------------------------------------------------
!!$    !Ozone and aerosols
!!$    !-------------------------------------------------------------------------
!!$
!!$      IF((irad_o3 == io3_clim) .OR. (irad_o3 == io3_ape) ) THEN
!!$        CALL message( TRIM(routine),'external ozone data required' )
!!$        CALL read_ext_data_atm (p_patch, ext_data, nlev_o3)   ! read ozone
!!$      ENDIF


    SELECT CASE(itopo)

    CASE(0) ! do not read external data  (except land-sea mask for JSBACH and ozone)
            ! topography from analytical functions 

    !-------------------------------------------------------------------------
    ! surface/vegetation  parameter
    !-------------------------------------------------------------------------

      !
      ! initalize external data with meaningful data, in the case that they 
      ! are not read in from file.
      SELECT CASE ( iforcing )
      CASE ( inwp )
        DO jg = 1, n_dom
          ext_data(jg)%atm%emis_rad(:,:)    = zemiss_def ! longwave surface emissivity
          ext_data(jg)%atm%fr_land(:,:)     = 0._wp      ! land fraction
          ext_data(jg)%atm%fr_land_smt(:,:) = 0._wp      ! land fraction (smoothed)
          ext_data(jg)%atm%fr_glac_smt(:,:) = 0._wp      ! glacier fraction (smoothed)
          ext_data(jg)%atm%llsm_atm_c(:,:)  = .FALSE.    ! land-sea mask
          ext_data(jg)%atm%plcov_mx(:,:)    = 0.5_wp     ! plant cover
          ext_data(jg)%atm%lai_mx(:,:)      = 3._wp      ! max Leaf area index
          ext_data(jg)%atm%rootdp(:,:)      = 1._wp      ! root depth
          ext_data(jg)%atm%rsmin(:,:)       = 150._wp    ! minimal stomata resistence
          ext_data(jg)%atm%soiltyp(:,:)     = 8          ! soil type
          ext_data(jg)%atm%z0(:,:)          = 0.001_wp   ! roughness length
        END DO
      CASE ( iecham, ildf_echam)
        DO jg = 1, n_dom
          ext_data(jg)%atm%emis_rad(:,:)    = zemiss_def ! longwave surface emissivity
        END DO
      END SELECT

      CALL message( TRIM(routine),'Running with analytical topography' )

      ! call read_ext_data_atm to read land-sea mask for JSBACH and to read O3
      IF (echam_phy_config%ljsbach .OR. irad_o3 == io3_clim .OR. irad_o3 == io3_ape) THEN
        CALL read_ext_data_atm (p_patch, ext_data, nlev_o3)
      END IF 
    CASE(1) ! read external data from netcdf dataset

      CALL message( TRIM(routine),'Start reading external data from file' )

      CALL read_ext_data_atm (p_patch, ext_data, nlev_o3)
      DO jg = 1, n_dom
        CALL smooth_topography (p_patch(jg), p_int_state(jg),  &
                                ext_data(jg)%atm%topography_c, &
                                ext_data(jg)%atm%topography_v, &
                                ext_data(jg)%atm%sso_stdh      )
      ENDDO

      CALL message( TRIM(routine),'Finished reading external data' )

      ! Get interpolated ndvi_mrat. Interpolation is done in time, based 
      ! on ini_datetime. For NWP applications it is assumed, that 
      ! ndvi_mrat is constant in time.
      !
      SELECT CASE ( iforcing )
      CASE ( inwp )
        CALL interpol_ndvi_time (p_patch, ext_data, time_config%ini_datetime)
      END SELECT

    CASE DEFAULT

      CALL finish( TRIM(routine), 'topography selection not supported' )

    END SELECT

    ! read external ocean data from ocean gridfile (no itopo needed)
    ELSE
      CALL read_ext_data_oce (p_patch, ext_data)
    END IF


  END SUBROUTINE init_ext_data


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
  SUBROUTINE init_ext_data_oce (p_patch, ext_data)

    TYPE(t_patch), INTENT(IN)            :: p_patch(:)
    TYPE(t_external_data), INTENT(INOUT) :: ext_data(:)


    INTEGER :: jg

    CHARACTER(len=max_char_length), PARAMETER :: &
      routine = 'mo_ext_data:init_ext_data'

    !-------------------------------------------------------------------------
    CALL message (TRIM(routine), 'Start')

    !-------------------------------------------------------------------------
    !  1.  inquire external files for their data structure
    !-------------------------------------------------------------------------

    ALLOCATE(nclass_lu(n_dom))
    ! Set default value for nclass_lu. Will be overwritten, if external data 
    ! are read from file
    nclass_lu(1:n_dom) = 1

    ALLOCATE(nmonths_ext(n_dom))
    ! Set default value for nmonths_ext. Will be overwritten, if external data 
    ! are read from file
    nmonths_ext(1:n_dom) = 1

    CALL inquire_external_files(p_patch)

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

      CALL read_ext_data_oce (p_patch, ext_data)
  END SUBROUTINE init_ext_data_oce

  !-------------------------------------------------------------------------
  !>
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

    IF (iequations/=ihs_ocean) THEN  ! atmosphere model ---------------------

      ! Build external data list for constant-in-time fields for the atm model
      write(*,*) 'create new external data list for atmosphere'
      DO jg = 1, n_dom
        WRITE(listname,'(a,i2.2)') 'ext_data_atm_D',jg
        CALL new_ext_data_atm_list(p_patch(jg), ext_data(jg)%atm,       &
          &                        ext_data(jg)%atm_list, TRIM(listname))
      END DO

      ! Build external data list for time-dependent fields
      IF (iforcing > 1 ) THEN ! further distinction is made inside
      DO jg = 1, n_dom
          WRITE(listname,'(a,i2.2)') 'ext_data_atm_td_D',jg
          CALL new_ext_data_atm_td_list(p_patch(jg), ext_data(jg)%atm_td,       &
            &                           ext_data(jg)%atm_td_list, TRIM(listname))
      END DO
      END IF

    ELSE ! iequations==ihs_ocean ------------------------------------------

      write(0,*) 'create new external data list for ocean'
      ! Build external data list for constant-in-time fields for the ocean model
      DO jg = 1, n_dom
        WRITE(listname,'(a,i2.2)') 'ext_data_oce_D',jg
        CALL new_ext_data_oce_list(p_patch(jg), ext_data(jg)%oce,       &
          &                        ext_data(jg)%oce_list, TRIM(listname))
      END DO ! jg = 1,n_dom

      ! Build external data list for time-dependent fields
      ! ### to be done ###

    ENDIF ! atmosphere or ocean ------

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

    INTEGER :: jg

    INTEGER :: nlev          !< number of vertical levels

    INTEGER :: nblks_c, &    !< number of cell blocks to allocate
      &        nblks_e, &    !< number of edge blocks to allocate
      &        nblks_v       !< number of vertex blocks to allocate

    INTEGER :: shape2d_c(2), shape2d_e(2), shape2d_v(2)
    INTEGER :: shape3d_c(3), shape3d_mon_c(3)
    INTEGER :: shape3d_sfc(3), shape3d_nt(3), shape3d_ntw(3)

    INTEGER :: ibits         !< "entropy" of horizontal slice

    INTEGER          :: jsfc
    CHARACTER(LEN=1) :: csfc

    !--------------------------------------------------------------

    !determine size of arrays
    nblks_c = p_patch%nblks_c
    nblks_e = p_patch%nblks_e
    nblks_v = p_patch%nblks_v

    ! get patch ID
    jg = p_patch%id
    ibits = DATATYPE_PACK16   ! "entropy" of horizontal slice

    ! number of vertical levels
    nlev = p_patch%nlev


    ! predefined array shapes
    shape2d_c  = (/ nproma, nblks_c /)
    shape2d_e  = (/ nproma, nblks_e /)
    shape2d_v  = (/ nproma, nblks_v /)
    shape3d_c  = (/ nproma, nlev, nblks_c       /)
    shape3d_mon_c = (/ nproma, 12, nblks_c       /)
    shape3d_sfc= (/ nproma, nblks_c, nclass_lu(jg) /) 
    shape3d_nt = (/ nproma, nblks_c, ntiles_total     /) 
    shape3d_ntw = (/ nproma, nblks_c, ntiles_total + ntiles_water /) 


    !
    ! Register a field list and apply default settings
    !
    CALL new_var_list( p_ext_atm_list, TRIM(listname), patch_id=p_patch%id )
    CALL default_var_list_settings( p_ext_atm_list,            &
                                  & lrestart=.FALSE.,          &
                                  & restart_type=FILETYPE_NC2  )


    ! topography height at cell center
    !
    ! topography_c  p_ext_atm%topography_c(nproma,nblks_c)
    cf_desc    = t_cf_var('surface_height', 'm', &
      &                   'geometric height of the earths surface above sea level', DATATYPE_FLT32)
    grib2_desc = t_grib2_var( 0, 3, 6, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_ext_atm_list, 'topography_c', p_ext_atm%topography_c,  &
      &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,             &
      &           grib2_desc, ldims=shape2d_c, loutput=.TRUE.,             &
      &           isteptype=TSTEP_CONSTANT )

    IF (echam_phy_config%ljsbach) THEN
    ! atmosphere land-sea-mask at surface on cell centers
    !
    ! lsm_ctr_c  p_ext_atm%lsm_ctr_c(nproma,nblks_c)
    cf_desc    = t_cf_var('Atmosphere model land-sea-mask at cell center', '-2/-1/1/2', &
      &                   'Atmosphere model land-sea-mask', DATATYPE_FLT32)
    grib2_desc = t_grib2_var( 192, 140, 219, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_ext_atm_list, 'lsm_ctr_c', p_ext_atm%lsm_ctr_c,        &
      &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,             &
                  grib2_desc, ldims=shape2d_c )
    ! elevation p_ext_atm%elevation_c(nproma,nblks_c)
    cf_desc    = t_cf_var('elevation at cell center', 'm', &
      &                   'elevation', DATATYPE_FLT32)
    grib2_desc = t_grib2_var( 192, 140, 219, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_ext_atm_list, 'elevation_c', p_ext_atm%elevation_c,        &
      &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,          &
                  grib2_desc, ldims=shape2d_c )
    ! SST  p_ext_atm%sst(nproma,nblks_c)
    cf_desc    = t_cf_var('SST', 'K', &
      &                   'SST as prescribed for atmosphere simulations', DATATYPE_FLT32)
    grib2_desc = t_grib2_var( 192, 140, 219, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_ext_atm_list, 'sst', p_ext_atm%sst,                    &
      &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,             &
                  grib2_desc, ldims=shape2d_c )
    ! monthly SST  p_ext_atm%sst_mon(nproma,nblks_c)
    cf_desc    = t_cf_var('SST monthly', 'K', &
      &                   'monthly SST as prescribed for atmosphere simulations', DATATYPE_FLT32)
    grib2_desc = t_grib2_var( 192, 140, 219, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_ext_atm_list, 'sst_mon', p_ext_atm%sst_mon,        &
      &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,          &
                  grib2_desc, ldims=shape3d_mon_c )
    ! albedo_vis_soil  p_ext_atm%albedo_vis_soil(nproma,nblks_c)
    cf_desc    = t_cf_var('soil surface albedo visible', '', &
      &                   'soil surface albedo in the visible range', DATATYPE_FLT32)
    grib2_desc = t_grib2_var( 192, 140, 219, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_ext_atm_list, 'albedo_vis_soil', p_ext_atm%albedo_vis_soil,        &
      &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,             &
                  grib2_desc, ldims=shape2d_c )

    ! albedo_nir_soil  p_ext_atm%albedo_nir_soil(nproma,nblks_c)
    cf_desc    = t_cf_var('soil surface albedo NIR', '', &
      &                   'soil surface albedo in the NIR range', DATATYPE_FLT32)
    grib2_desc = t_grib2_var( 192, 140, 219, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_ext_atm_list, 'albedo_nir_soil', p_ext_atm%albedo_nir_soil,        &
      &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,             &
                  grib2_desc, ldims=shape2d_c )

    ! albedo_vis_canopy  p_ext_atm%albedo_vis_canopy(nproma,nblks_c)
    cf_desc    = t_cf_var('canopy reflectance visible', '', &
      &                   'canopy reflectance in the visible range', DATATYPE_FLT32)
    grib2_desc = t_grib2_var( 192, 140, 219, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_ext_atm_list, 'albedo_vis_canopy', p_ext_atm%albedo_vis_canopy,    &
      &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,             &
                  grib2_desc, ldims=shape2d_c )

    ! albedo_nir_canopy  p_ext_atm%albedo_nir_canopy(nproma,nblks_c)
    cf_desc    = t_cf_var('canopy reflectance NIR', '', &
      &                   'canopy reflectance in the NIR range', DATATYPE_FLT32)
    grib2_desc = t_grib2_var( 192, 140, 219, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_ext_atm_list, 'albedo_nir_canopy', p_ext_atm%albedo_nir_canopy,    &
      &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,             &
                  grib2_desc, ldims=shape2d_c )

    ! albedo_background  p_ext_atm%albedo_background(nproma,nblks_c)
    cf_desc    = t_cf_var('background albedo', '', &
      &                   'background albedo for echam5 albedo scheme', DATATYPE_FLT32)
    grib2_desc = t_grib2_var( 192, 140, 219, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_ext_atm_list, 'albedo_background', p_ext_atm%albedo_background,    &
      &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,             &
                  grib2_desc, ldims=shape2d_c )

    ! forest_fract  p_ext_atm%forest_fract(nproma,nblks_c)
    cf_desc    = t_cf_var('forest fraction', '', &
      &                   'forest fraction for echam5 albedo scheme', DATATYPE_FLT32)
    grib2_desc = t_grib2_var( 192, 140, 219, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_ext_atm_list, 'forest_fract', p_ext_atm%forest_fract,   &
      &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,              &
                  grib2_desc, ldims=shape2d_c )
    END IF

    ! ozone mixing ratio
    !
    ! o3            p_ext_atm%o3(nproma,nlev,nblks_c)
    cf_desc    = t_cf_var('ozone mixing ratio', 'kg kg-1', &
      &                   'ozone mixing ratio', DATATYPE_FLT32)
    grib2_desc = t_grib2_var( 0, 14, 1, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_ext_atm_list, 'o3', p_ext_atm%o3,                      &
      &           GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc,              &
      &           grib2_desc, ldims=shape3d_c, loutput=.TRUE. )


  SELECT CASE ( iequations )
  CASE ( inh_atmosphere )
    ! smoothed topography height at cell center
    !
    ! topography_smt_c  p_ext_atm%topography_smt_c(nproma,nblks_c)
    cf_desc    = t_cf_var('smoothed_surface_height', 'm', &
      &                   'smoothed geometric height of the earths surface above sea level', &
      &                   DATATYPE_FLT32)
    grib2_desc = t_grib2_var( 2, 0, 7, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_ext_atm_list, 'topography_smt_c', p_ext_atm%topography_smt_c, &
      &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,                    &
      &           grib2_desc, ldims=shape2d_c, loutput=.FALSE.,                   &
      &           isteptype=TSTEP_CONSTANT )


    ! topography height at edge midpoint
    !
    ! topography_e  p_ext_atm%topography_e(nproma,nblks_e)
    cf_desc    = t_cf_var('surface_height', 'm', &
      &                   'geometric height of the earths surface above sea level', DATATYPE_FLT32)
    grib2_desc = t_grib2_var( 0, 3, 6, ibits, GRID_REFERENCE, GRID_EDGE)
    CALL add_var( p_ext_atm_list, 'topography_e', p_ext_atm%topography_e, &
      &           GRID_UNSTRUCTURED_EDGE, ZA_SURFACE, cf_desc,            &
      &           grib2_desc, ldims=shape2d_e, loutput=.FALSE.,           &
      &           isteptype=TSTEP_CONSTANT )


    ! topography height at vertex
    !
    ! topography_v  p_ext_atm%topography_v(nproma,nblks_v)
    cf_desc    = t_cf_var('surface_height', 'm', &
      &                   'geometric height of the earths surface above sea level', DATATYPE_FLT32)
    grib2_desc = t_grib2_var( 0, 3, 6, ibits, GRID_REFERENCE, GRID_VERTEX)
    CALL add_var( p_ext_atm_list, 'topography_v', p_ext_atm%topography_v, &
      &           GRID_UNSTRUCTURED_VERT, ZA_SURFACE, cf_desc,            &
      &           grib2_desc, ldims=shape2d_v, loutput=.FALSE.,           &
      &           isteptype=TSTEP_CONSTANT )


    ! smoothed topography height at vertex
    !
    ! topography_smt_v  p_ext_atm%topography_smt_v(nproma,nblks_v)
    cf_desc    = t_cf_var('smoothed_surface_height', 'm', &
      &                   'smoothed geometric height of the earths surface above sea level', &
      &                   DATATYPE_FLT32)
    grib2_desc = t_grib2_var( 2, 0, 7, ibits, GRID_REFERENCE, GRID_VERTEX)
    CALL add_var( p_ext_atm_list, 'topography_smt_v', p_ext_atm%topography_smt_v, &
      &           GRID_UNSTRUCTURED_VERT, ZA_SURFACE, cf_desc,                    &
      &           grib2_desc, ldims=shape2d_v, loutput=.FALSE.,                   &
      &           isteptype=TSTEP_CONSTANT )



    ! land sea mask for cells (LOGICAL)
    ! Note: Here "loutput" is set to .FALSE. since the output
    !       scheme operates on REAL model variables only and
    !       throws an error on this.
    !
    ! llsm_atm_c    p_ext_atm%llsm_atm_c(nproma,nblks_c)
    cf_desc    = t_cf_var('land_sea_mask_(cell)', '-', &
      &                   'land sea mask (cell)', DATATYPE_FLT32)
    grib2_desc = t_grib2_var( 2, 0, 0, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_ext_atm_list, 'llsm_atm_c', p_ext_atm%llsm_atm_c, &
      &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,        &
      &           grib2_desc, ldims=shape2d_c, loutput=.FALSE.,       &
      &           isteptype=TSTEP_CONSTANT )

    ! land fraction
    !
    ! fr_land      p_ext_atm%fr_land(nproma,nblks_c)
    cf_desc    = t_cf_var('land_area_fraction', '-', 'Fraction land', DATATYPE_FLT32)
    grib2_desc = t_grib2_var( 2, 0, 0, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_ext_atm_list, 'fr_land', p_ext_atm%fr_land,   &
      &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,    &
      &           grib2_desc, ldims=shape2d_c, loutput=.TRUE.,    &
      &           isteptype=TSTEP_CONSTANT )


    ! glacier fraction
    !
    ! fr_glac      p_ext_atm%fr_glac(nproma,nblks_c)
    cf_desc    = t_cf_var('glacier_area_fraction', '-', 'Fraction glacier', DATATYPE_FLT32)
    grib2_desc = t_grib2_var( 2, 0, 192, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_ext_atm_list, 'fr_glac', p_ext_atm%fr_glac,   &
      &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,    &
      &           grib2_desc, ldims=shape2d_c, loutput=.FALSE. )


    ! maybe the next three (ice, fr_land_smt, fr_ice_smt)
    ! should be moved into corresponding if block

    ! sea Ice fraction
    !
    ! fr_ice       p_ext_atm%fr_ice(nproma,nblks_c)
    cf_desc    = t_cf_var('Sea_ice_fraction', '-', 'Sea ice fraction', DATATYPE_FLT32)
    grib2_desc = t_grib2_var( 10, 2, 0, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_ext_atm_list, 'fr_ice', p_ext_atm%fr_ice,     &
      &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,    &
      &           grib2_desc, ldims=shape2d_c, loutput=.FALSE. )


    ! land fraction (smoothed)
    !
    ! fr_land_smt  p_ext_atm%fr_land_smt(nproma,nblks_c)
    cf_desc    = t_cf_var('land_area_fraction_(smoothed)', '-', &
      &                   'land area fraction (smoothed)', DATATYPE_FLT32)
    grib2_desc = t_grib2_var( 2, 0, 0, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_ext_atm_list, 'fr_land_smt', p_ext_atm%fr_land_smt, &
      &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,          &
      &           grib2_desc, ldims=shape2d_c, loutput=.FALSE.,         &
      &           isteptype=TSTEP_CONSTANT )


    ! glacier area fraction (smoothed)
    !
    ! fr_glac_smt  p_ext_atm%fr_glac_smt(nproma,nblks_c)
    cf_desc    = t_cf_var('glacier_area_fraction_(smoothed)', '-', &
      &                   'glacier area fraction (smoothed)', DATATYPE_FLT32)
    grib2_desc = t_grib2_var( 2, 0, 192, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_ext_atm_list, 'fr_glac_smt', p_ext_atm%fr_glac_smt, &
      &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,          &
      &           grib2_desc, ldims=shape2d_c, loutput=.FALSE. )


    ! sea Ice fraction (smoothed)
    !
    ! fr_ice_smt  p_ext_atm%fr_ice_smt(nproma,nblks_c)
    cf_desc    = t_cf_var('Sea_ice_fraction (smoothed)', '-', &
      &                   'Sea ice fraction (smoothed)', DATATYPE_FLT32)
    grib2_desc = t_grib2_var( 10, 2, 0, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_ext_atm_list, 'fr_ice_smt', p_ext_atm%fr_ice_smt, &
      &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,        &
      &           grib2_desc, ldims=shape2d_c, loutput=.FALSE. )


    ! geopotential (s)
    !
    ! fis          p_ext_atm%fis(nproma,nblks_c)
    cf_desc    = t_cf_var('Geopotential_(s)', 'm2 s-2', &
      &                   'Geopotential (s)', DATATYPE_FLT32)
    grib2_desc = t_grib2_var( 0, 3, 4, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_ext_atm_list, 'fis', p_ext_atm%fis,           &
      &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,    &
      &           grib2_desc, ldims=shape2d_c, loutput=.TRUE. )



    SELECT CASE ( iforcing )
    CASE ( inwp )
      ! external parameter for NWP forcing

      ! roughness length
      !
      ! z0           p_ext_atm%z0(nproma,nblks_c)
      cf_desc    = t_cf_var('roughtness_length', 'm', 'roughtness length', DATATYPE_FLT32)
      grib2_desc = t_grib2_var( 2, 0, 1, ibits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'z0', p_ext_atm%z0,             &
        &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,    &
        &           grib2_desc, ldims=shape2d_c, loutput=.FALSE. )


      ! fraction lake
      !
      ! fr_lake      p_ext_atm%fr_lake(nproma,nblks_c)
      cf_desc    = t_cf_var('fraction_lake', '-', 'fraction lake', DATATYPE_FLT32)
      grib2_desc = t_grib2_var( 1, 2, 2, ibits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'fr_lake', p_ext_atm%fr_lake,   &
        &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,    &
        &           grib2_desc, ldims=shape2d_c, loutput=.FALSE.,   &
        &           isteptype=TSTEP_CONSTANT )


      ! lake depth
      !
      ! depth_lk     p_ext_atm%depth_lk(nproma,nblks_c)
      cf_desc    = t_cf_var('lake_depth', '-', 'lake depth', DATATYPE_FLT32)
      grib2_desc = t_grib2_var( 1, 2, 0, ibits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'depth_lk', p_ext_atm%depth_lk, &
        &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,    &
        &           grib2_desc, ldims=shape2d_c, loutput=.FALSE.,   &
        &           isteptype=TSTEP_CONSTANT )




      !--------------------------------
      ! sub-gridscale orography
      !--------------------------------

      ! Standard deviation of sub-grid scale orography
      !
      ! sso_stdh     p_ext_atm%sso_stdh(nproma,nblks_c)
      cf_desc    = t_cf_var('standard_deviation_of_height', 'm',    &
        &                   'Standard deviation of sub-grid scale orography', DATATYPE_FLT32)
      grib2_desc = t_grib2_var( 0, 3, 20, ibits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'sso_stdh', p_ext_atm%sso_stdh, &
        &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,    &
        &           grib2_desc, ldims=shape2d_c, loutput=.FALSE.,   &
        &           isteptype=TSTEP_CONSTANT )



      ! Anisotropy of sub-gridscale orography
      !
      ! sso_gamma    p_ext_atm%sso_gamma(nproma,nblks_c)
      cf_desc    = t_cf_var('anisotropy_factor', '-',&
        &                   'Anisotropy of sub-gridscale orography', DATATYPE_FLT32)
      grib2_desc = t_grib2_var( 0, 3, 20, ibits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'sso_gamma', p_ext_atm%sso_gamma, &
        &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,      &
        &           grib2_desc, ldims=shape2d_c, loutput=.FALSE.,     &
        &           isteptype=TSTEP_CONSTANT )



      ! Angle of sub-gridscale orography
      !
      ! sso_theta    p_ext_atm%sso_theta(nproma,nblks_c)
      cf_desc    = t_cf_var('angle_of_principal_axis', 'radians',&
        &                   'Angle of sub-gridscale orography', DATATYPE_FLT32)
      grib2_desc = t_grib2_var( 0, 3, 21, ibits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'sso_theta', p_ext_atm%sso_theta, &
        &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,      &
        &           grib2_desc, ldims=shape2d_c, loutput=.FALSE.,     &
        &           isteptype=TSTEP_CONSTANT )



      ! Slope of sub-gridscale orography
      !
      ! sso_sigma    p_ext_atm%sso_sigma(nproma,nblks_c)
      cf_desc    = t_cf_var('slope_of_terrain', '-',&
        &                   'Slope of sub-gridscale orography', DATATYPE_FLT32)
      grib2_desc = t_grib2_var( 0, 3, 22, ibits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'sso_sigma', p_ext_atm%sso_sigma, &
        &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,      &
        &           grib2_desc, ldims=shape2d_c, loutput=.FALSE.,     &
        &           isteptype=TSTEP_CONSTANT )





      !--------------------------------
      ! vegetation parameters
      !--------------------------------

      ! Plant covering degree in the vegetation phase
      !
      ! plcov_mx     p_ext_atm%plcov_mx(nproma,nblks_c)
      cf_desc    = t_cf_var('vegetation_area_fraction_vegetation_period', '-',&
        &                   'Plant covering degree in the vegetation phase', DATATYPE_FLT32)
      grib2_desc = t_grib2_var( 2, 0, 4, ibits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'plcov_mx', p_ext_atm%plcov_mx, &
        &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,    &
        &           grib2_desc, ldims=shape2d_c, loutput=.FALSE. )

      ! plcov_t     p_ext_atm%plcov_t(nproma,nblks_c,ntiles_total)
      cf_desc    = t_cf_var('vegetation_area_fraction_vegetation_period', '-',&
        &                   'Plant covering degree in the vegetation phase', DATATYPE_FLT32)
      grib2_desc = t_grib2_var( 2, 0, 4, ibits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'plcov_t', p_ext_atm%plcov_t, &
        &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,  &
        &           grib2_desc, ldims=shape3d_nt, loutput=.FALSE. )


      ! Max Leaf area index
      !
      ! lai_mx       p_ext_atm%lai_mx(nproma,nblks_c)
      cf_desc    = t_cf_var('leaf_area_index_vegetation_period', '-',&
        &                   'Leaf Area Index Maximum', DATATYPE_FLT32)
      grib2_desc = t_grib2_var( 2, 0, 28, ibits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'lai_mx', p_ext_atm%lai_mx,     &
        &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,    &
        &           grib2_desc, ldims=shape2d_c, loutput=.FALSE. )

      ! sai_t       p_ext_atm%sai_t(nproma,nblks_c,ntiles_total+ntiles_water)
      cf_desc    = t_cf_var('surface_area_index_vegetation_period', '-',&
        &                   'Surface Area Index', DATATYPE_FLT32)
      grib2_desc = t_grib2_var( 2, 0, 28, ibits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'sai_t', p_ext_atm%sai_t,     &
        &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,  &
        &           grib2_desc, ldims=shape3d_ntw, loutput=.FALSE. )

      ! tai_t       p_ext_atm%tai_t(nproma,nblks_c,ntiles_total)
      cf_desc    = t_cf_var('transpiration_area_index_vegetation_period', '-',&
        &                   'Transpiration Area Index', DATATYPE_FLT32)
      grib2_desc = t_grib2_var( 2, 0, 28, ibits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'tai_t', p_ext_atm%tai_t,     &
        &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,  &
        &           grib2_desc, ldims=shape3d_nt, loutput=.FALSE. )

      ! eai_t       p_ext_atm%eai_t(nproma,nblks_c,ntiles_total)
      cf_desc    = t_cf_var('evaporative_surface_area_index_vegetation_period', '-',&
        &                   'Earth Area (evaporative surface area) Index', DATATYPE_FLT32)
      grib2_desc = t_grib2_var( 2, 0, 28, ibits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'eai_t', p_ext_atm%eai_t,     &
        &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,  &
        &           grib2_desc, ldims=shape3d_nt, loutput=.FALSE. )


      ! root depth of vegetation
      !
      ! rootdp      p_ext_atm%rootdp(nproma,nblks_c)
      cf_desc    = t_cf_var('root_depth_of_vegetation', 'm',&
        &                   'root depth of vegetation', DATATYPE_FLT32)
      grib2_desc = t_grib2_var( 2, 0, 32, ibits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'rootdp', p_ext_atm%rootdp,     &
        &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,    &
        &           grib2_desc, ldims=shape2d_c, loutput=.FALSE. )

      ! rootdp_t      p_ext_atm%rootdp_t(nproma,nblks_c,ntiles_total)
      cf_desc    = t_cf_var('root_depth_of_vegetation', 'm',&
        &                   'root depth of vegetation', DATATYPE_FLT32)
      grib2_desc = t_grib2_var( 2, 0, 32, ibits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'rootdp_t', p_ext_atm%rootdp_t, &
        &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,    &
        &           grib2_desc, ldims=shape3d_nt, loutput=.FALSE. )

      ! evergreen forest
      !
      ! for_e        p_ext_atm%for_e(nproma,nblks_c)
      cf_desc    = t_cf_var('fraction_of_evergreen_forest_cover', '-',&
        &                   'Fraction of evergreen forest', DATATYPE_FLT32)
      grib2_desc = t_grib2_var( 2, 0, 29, ibits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'for_e', p_ext_atm%for_e,       &
        &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,    &
        &           grib2_desc, ldims=shape2d_c, loutput=.FALSE. )
 


      ! deciduous forest
      !
      ! for_d     p_ext_atm%for_d(nproma,nblks_c)
      cf_desc    = t_cf_var('fraction_of_deciduous_forest_cover', '-',&
        &                   'Fraction of deciduous forest', DATATYPE_FLT32)
      grib2_desc = t_grib2_var( 2, 0, 30, ibits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'for_d', p_ext_atm%for_d,       &
        &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,    &
        &           grib2_desc, ldims=shape2d_c )



      ! urban area fraction
      !
      ! urban        p_ext_atm%urban(nproma,nblks_c)
      cf_desc    = t_cf_var('fraction_of_urban_areas', '-',&
        &                   'urban area fraction', DATATYPE_FLT32)
      grib2_desc = t_grib2_var( 2, 0, 30, ibits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'urban', p_ext_atm%urban,       &
        &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,    &
        &           grib2_desc, ldims=shape2d_c )


      ! Minimal stomata resistence
      !
      ! rsmin        p_ext_atm%rsmin(nproma,nblks_c)
      cf_desc    = t_cf_var('RSMIN', 's m-1', 'Minimal stomata resistence', DATATYPE_FLT32)
      grib2_desc = t_grib2_var( 2, 0, 16, ibits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'rsmin', p_ext_atm%rsmin,       &
        &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,    &
        &           grib2_desc, ldims=shape2d_c )

      ! rsmin2d_t        p_ext_atm%rsmin2d_t(nproma,nblks_c,ntiles_total)
      cf_desc    = t_cf_var('RSMIN', 's m-1', 'Minimal stomata resistence', DATATYPE_FLT32)
      grib2_desc = t_grib2_var( 2, 0, 16, ibits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'rsmin2d_t', p_ext_atm%rsmin2d_t,       &
        &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,            &
        &           grib2_desc, ldims=shape3d_nt, loutput=.FALSE. )

      ! NDVI yearly maximum
      !
      ! ndvi_max        p_ext_atm%ndvi_max(nproma,nblks_c)
      cf_desc    = t_cf_var('normalized_difference_vegetation_index', '-', &
        &                   'NDVI yearly maximum', DATATYPE_FLT32)
      grib2_desc = t_grib2_var( 2, 0, 31, ibits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'ndvi_max', p_ext_atm%ndvi_max, &
        &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,    &
        &           grib2_desc, ldims=shape2d_c, loutput=.FALSE.  )

      ! proportion of actual value/maximum NDVI (at ini_datetime)
      !
      ! ndvi_mrat        p_ext_atm%ndvi_mrat(nproma,nblks_c)
      cf_desc    = t_cf_var('normalized_difference_vegetation_index', '-',     &
        &                   '(monthly) proportion of actual value/maximum ' // &
        &                   'NDVI (at init time)', DATATYPE_FLT32)
      grib2_desc = t_grib2_var( 2, 0, 192, ibits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'ndvi_mrat', p_ext_atm%ndvi_mrat, &
        &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,      &
        &           grib2_desc, ldims=shape2d_c, loutput=.FALSE.  )

      ! Control fields for tile approach
      ! idx_lst_lp          p_ext_atm%idx_lst_lp(nproma,nblks_c)
      cf_desc    = t_cf_var('land point index list', '-', &
        &                   'land point index list', DATATYPE_FLT32)
      grib2_desc = t_grib2_var( 255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'idx_lst_lp', p_ext_atm%idx_lst_lp, &
        &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,        &
        &           grib2_desc, ldims=shape2d_c, loutput=.FALSE. )

      ! idx_lst_sp          p_ext_atm%idx_lst_sp(nproma,nblks_c)
      cf_desc    = t_cf_var('sea point index list', '-', &
        &                   'sea point index list', DATATYPE_FLT32)
      grib2_desc = t_grib2_var( 255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'idx_lst_sp', p_ext_atm%idx_lst_sp, &
        &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,        &
        &           grib2_desc, ldims=shape2d_c, loutput=.FALSE. )

      ! idx_lst_fp          p_ext_atm%idx_lst_sp(nproma,nblks_c)
      cf_desc    = t_cf_var('lake point index list', '-', &
        &                   'lake point index list', DATATYPE_FLT32)
      grib2_desc = t_grib2_var( 255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'idx_lst_fp', p_ext_atm%idx_lst_fp, &
        &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,        &
        &           grib2_desc, ldims=shape2d_c, loutput=.FALSE. )

      ! idx_lst_lp_t        p_ext_atm%idx_lst_lp_t(nproma,nblks_c,ntiles_total)
      cf_desc    = t_cf_var('static land tile point index list', '-', &
        &                   'static land tile point index list', DATATYPE_FLT32)
      grib2_desc = t_grib2_var( 255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'idx_lst_lp_t', p_ext_atm%idx_lst_lp_t, &
        &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,            &
        &           grib2_desc, ldims=shape3d_nt, loutput=.FALSE. )

      ! idx_lst_t        p_ext_atm%idx_lst_t(nproma,nblks_c,ntiles_total)
      cf_desc    = t_cf_var('dynamic land tile point index list', '-', &
        &                   'dynamic land tile point index list', DATATYPE_FLT32)
      grib2_desc = t_grib2_var( 255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'idx_lst_t', p_ext_atm%idx_lst_t, &
        &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,      &
        &           grib2_desc, ldims=shape3d_nt, loutput=.FALSE. )


      ! idx_lst_spw      p_ext_atm%idx_lst_spw(nproma,nblks_c)
      cf_desc    = t_cf_var('sea water point index list', '-', &
        &                   'sea water point index list', DATATYPE_FLT32)
      grib2_desc = t_grib2_var( 255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'idx_lst_spw', p_ext_atm%idx_lst_spw, &
        &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,          &
        &           grib2_desc, ldims=shape2d_c, loutput=.FALSE. )

      ! idx_lst_spi      p_ext_atm%idx_lst_spi(nproma,nblks_c)
      cf_desc    = t_cf_var('sea ice point index list', '-', &
        &                   'sea ice point index list', DATATYPE_FLT32)
      grib2_desc = t_grib2_var( 255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'idx_lst_spi', p_ext_atm%idx_lst_spi, &
        &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,          &
        &           grib2_desc, ldims=shape2d_c, loutput=.FALSE. )


      ! snowtile_flag_t   p_ext_atm%snowtile_flag_t(nproma,nblks_c,ntiles_total)
      ! -1: no separation between snow tile and snow-free tile
      !  0: inactive
      !  1: active
      !  2: newly activated; initialization from corresponding tile required
      cf_desc    = t_cf_var('flag of activity', '-', &
        &                   'flag of activity', DATATYPE_FLT32)
      grib2_desc = t_grib2_var( 255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'snowtile_flag_t', p_ext_atm%snowtile_flag_t, &
        &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,                  &
        &           grib2_desc, ldims=shape3d_nt, loutput=.FALSE. )

      ! not sure if these dimensions are supported by add_var...
      ALLOCATE(p_ext_atm%lp_count(nblks_c), p_ext_atm%gp_count_t(nblks_c,ntiles_total), &
               p_ext_atm%lp_count_t(nblks_c,ntiles_total) )
      ALLOCATE(p_ext_atm%sp_count (nblks_c),p_ext_atm%fp_count (nblks_c))

      ! allocate grid point counts per block for dynamic ocean ice/water point 
      ! index lists
      ALLOCATE(p_ext_atm%spw_count(nblks_c),p_ext_atm%spi_count(nblks_c))


      ! lc_class_t        p_ext_atm%lc_class_t(nproma,nblks_c,ntiles_total+ntiles_water)
      cf_desc    = t_cf_var('tile point land cover class list', '-', &
        &                   'tile point land cover class list', DATATYPE_FLT32)
      grib2_desc = t_grib2_var( 255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'lc_class_t', p_ext_atm%lc_class_t, &
        &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,        &
        &           grib2_desc, ldims=shape3d_ntw, loutput=.FALSE. )

      ! lc_frac_t        p_ext_atm%lc_frac_t(nproma,nblks_c,ntiles_total+ntiles_water)
      cf_desc    = t_cf_var('tile point land cover fraction list', '-', &
        &                   'tile point land cover fraction list', DATATYPE_FLT32)
      grib2_desc = t_grib2_var( 255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'lc_frac_t', p_ext_atm%lc_frac_t, &
        &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,      &
        &           grib2_desc, ldims=shape3d_ntw, loutput=.FALSE. )

      ! frac_t        p_ext_atm%frac_t(nproma,nblks_c,ntiles_total+ntiles_water)
      cf_desc    = t_cf_var('tile point area fraction list', '-', &
        &                   'tile point area fraction list', DATATYPE_FLT32)
      grib2_desc = t_grib2_var( 255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'frac_t', p_ext_atm%frac_t,   &
        &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,  &
        &           grib2_desc, ldims=shape3d_ntw, loutput=.FALSE., lcontainer=.TRUE.)

      ! fill the separate variables belonging to the container frac_t
      ALLOCATE(p_ext_atm%frac_t_ptr(ntiles_total+ntiles_water))
      DO jsfc = 1,ntiles_total + ntiles_water
      WRITE(csfc,'(i1)') jsfc
      CALL add_ref( p_ext_atm_list, 'frac_t', 'frac_t_'//TRIM(ADJUSTL(csfc)),  &
        &           p_ext_atm%frac_t_ptr(jsfc)%p_2d,                           &
        &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,   &
        &           ldims=shape2d_c, loutput=.TRUE. )
      ENDDO



      ! Storage for table values - not sure if these dimensions are supported by add_var
      ! The dimension (num_lcc) is currently hard-wired to 23
       ALLOCATE(p_ext_atm%z0_lcc(num_lcc),         & ! Land-cover related roughness length
                p_ext_atm%plcovmax_lcc(num_lcc),   & ! Maximum plant cover fraction for each land-cover class
                p_ext_atm%laimax_lcc(num_lcc),     & ! Maximum leaf area index for each land-cover class
                p_ext_atm%rootdmax_lcc(num_lcc),   & ! Maximum root depth each land-cover class
                p_ext_atm%stomresmin_lcc(num_lcc), & ! Minimum stomata resistance for each land-cover class
                p_ext_atm%snowalb_lcc(num_lcc),    & ! Albedo in case of snow cover for each land-cover class
                p_ext_atm%snowtile_lcc(num_lcc)    ) ! Specification of snow tiles for land-cover class


      !--------------------------------
      ! soil parameters
      !--------------------------------

      ! soil type
      !
      ! soiltyp      p_ext_atm%soiltyp(nproma,nblks_c)
      cf_desc    = t_cf_var('soil_type', '-','soil type', DATATYPE_FLT32)
      grib2_desc = t_grib2_var( 2, 3, 0, ibits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'soiltyp', p_ext_atm%soiltyp,   &
        &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,    &
        &           grib2_desc, ldims=shape2d_c, loutput=.FALSE. )

      ! soiltyp_t      p_ext_atm%soiltyp_t(nproma,nblks_c,ntiles_total)
      cf_desc    = t_cf_var('soil_type', '-','soil type', DATATYPE_FLT32)
      grib2_desc = t_grib2_var( 2, 3, 0, ibits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'soiltyp_t', p_ext_atm%soiltyp_t,   &
        &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,        &
        &           grib2_desc, ldims=shape3d_nt, loutput=.FALSE. )


      ! Climat. temperature
      !
      ! t_cl         p_ext_atm%t_cl(nproma,nblks_c)
      cf_desc    = t_cf_var('soil_temperature', 'K',                  &
        &                   'CRU near surface temperature climatology', DATATYPE_FLT32)
      grib2_desc = t_grib2_var( 2, 3, 18, ibits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_ext_atm_list, 't_cl', p_ext_atm%t_cl,         &
        &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,    &
        &           grib2_desc, ldims=shape2d_c, loutput=.FALSE.  )


      ! longwave surface emissivity
      !
      ! emis_rad     p_ext_atm%emis_rad(nproma,nblks_c)
      cf_desc    = t_cf_var('emis_rad', '-', 'longwave surface emissivity', DATATYPE_FLT32)
      grib2_desc = t_grib2_var( 2, 3, 196, ibits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'emis_rad', p_ext_atm%emis_rad, &
        &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,    &
        &           grib2_desc, ldims=shape2d_c, loutput=.FALSE. )


      ! landuse class fraction
      !
      ! lu_class_fraction    p_ext_atm%lu_class_fraction(nproma,nblks_c,nclass_lu)
      cf_desc    = t_cf_var('lu_class_fraction', '-', 'landuse class fraction', DATATYPE_FLT32)
      grib2_desc = t_grib2_var( 255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'lu_class_fraction', p_ext_atm%lu_class_fraction, &
        &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,    &
        &           grib2_desc, ldims=shape3d_sfc, loutput=.FALSE. )

    CASE ( iecham, ildf_echam )

      ! longwave surface emissivity
      !
      ! emis_rad     p_ext_atm%emis_rad(nproma,nblks_c)
      cf_desc    = t_cf_var('emis_rad', '-', 'longwave surface emissivity', DATATYPE_FLT32)
      grib2_desc = t_grib2_var( 2, 3, 196, ibits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'emis_rad', p_ext_atm%emis_rad, &
        &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,    &
        &           grib2_desc, ldims=shape2d_c, loutput=.FALSE. )

    END SELECT ! iforcing

  CASE ( ihs_atm_temp, ihs_atm_theta )

    SELECT CASE ( iforcing )
    CASE ( iecham, ildf_echam)
      ! longwave surface emissivity
      !
      ! emis_rad     p_ext_atm%emis_rad(nproma,nblks_c)
      cf_desc    = t_cf_var('emis_rad', '-', 'longwave surface emissivity', DATATYPE_FLT32)
      grib2_desc = t_grib2_var( 2, 3, 196, ibits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'emis_rad', p_ext_atm%emis_rad, &
        &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,    &
        &           grib2_desc, ldims=shape2d_c, loutput=.FALSE. )     
    END SELECT
        
  END SELECT ! iequations


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
    INTEGER :: jg           !< patch ID

    INTEGER :: shape3d_c(3)
    INTEGER :: shape3d_ape(3)
    INTEGER :: shape4d_c(4)
    INTEGER :: shape3d_sstice(3)

    INTEGER :: ibits         !< "entropy" of horizontal slice

    CHARACTER(len=max_char_length), PARAMETER :: &
      routine = 'mo_ext_data_state:new_ext_data_atm_td_list'
    !--------------------------------------------------------------

    !determine size of arrays
    nblks_c = p_patch%nblks_c

    ! get patch ID
    jg = p_patch%id

    ibits  = 16   ! "entropy" of horizontal slice

    ! predefined array shapes
    shape3d_c   = (/ nproma, nblks_c, nmonths_ext(jg)  /)
    shape4d_c   = (/ nproma, nlev_o3, nblks_c, nmonths /) 
    shape3d_ape = (/ nproma, nlev_o3, nblks_c          /) 

    IF ( sstice_mode > 1 ) THEN
     SELECT CASE (sstice_mode)
     CASE(2)
      shape3d_sstice = (/ nproma, nblks_c, 12 /)  
     CASE(3)
      shape3d_sstice = (/ nproma, nblks_c,  2 /)  
     CASE(4)
      CALL finish (TRIM(routine), 'sstice_mode=4  not implemented!')
     CASE DEFAULT
      CALL finish (TRIM(routine), 'sstice_mode not valid!')
     END SELECT
    END IF

    !
    ! Register a field list and apply default settings
    !
    CALL new_var_list( p_ext_atm_td_list, TRIM(listname), patch_id=p_patch%id )
    CALL default_var_list_settings( p_ext_atm_td_list,         &
                                  & lrestart=.FALSE.,          &
                                  & loutput=.FALSE.,           &
                                  & restart_type=FILETYPE_NC2  )


    !-------------------------------- 
    ! radiation parameters
    !--------------------------------


    ! ozone on pressure levels
    ! ATTENTION: a GRIB2 number will go to 
    ! the ozone mass mixing ratio...
    !
    IF(irad_o3 == io3_clim .OR. irad_o3 == io3_ape) THEN 

      WRITE(0,*) 'generate ext ozone field'

      ! o3  main height level from read-in file
      cf_desc    = t_cf_var('O3_zf', 'm',   &
        &                   'ozone geometric height level', DATATYPE_FLT32)
      grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_ext_atm_td_list, 'O3_zf', p_ext_atm_td%zf,  &
        &           GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc,   &
        &           grib2_desc, ldims=(/nlev_o3/), loutput=.FALSE.  )

      ! o3  main pressure level from read-in file
      cf_desc    = t_cf_var('O3_pf', 'Pa',   &
        &                   'ozone main pressure level', DATATYPE_FLT32)
      grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_ext_atm_td_list, 'O3_pf', p_ext_atm_td%pfoz, &
        &           GRID_UNSTRUCTURED_CELL, ZA_PRESSURE, cf_desc,  &
        &           grib2_desc, ldims=(/nlev_o3/), loutput=.FALSE.  )

      ! o3  intermediate pressure level
      cf_desc    = t_cf_var('O3_ph', 'Pa',   &
        &                   'ozone intermediate pressure level', DATATYPE_FLT32)
      grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_ext_atm_td_list, 'O3_ph', p_ext_atm_td%phoz, &
        &           GRID_UNSTRUCTURED_CELL, ZA_PRESSURE, cf_desc,  &
        &           grib2_desc, ldims=(/nlev_o3+1/), loutput=.FALSE.  )

        ! o3       p_ext_atm_td%o3(nproma,nlev_o3,nblks_c,nmonths)
        cf_desc    = t_cf_var('O3', TRIM(o3unit),   &
          &                   'mole_fraction_of_ozone_in_air', DATATYPE_FLT32)
        grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
        CALL add_var( p_ext_atm_td_list, 'O3', p_ext_atm_td%O3, &
          &           GRID_UNSTRUCTURED_CELL, ZA_PRESSURE, cf_desc, &
          &           grib2_desc, ldims=shape4d_c, loutput=.FALSE.  )

    END IF ! irad_o3

    IF(iforcing == inwp) THEN

    ! Black carbon aerosol
    !
    ! aer_bc       p_ext_atm%aer_bc(nproma,nblks_c,ntimes)
    cf_desc    = t_cf_var('aerosol optical thickness of black carbon', '-',   &
      &                   'atmosphere_absorption_optical_thickness_due_to_' //&
      &                   'black_carbon_ambient_aerosol', DATATYPE_FLT32)
    grib2_desc = t_grib2_var( 0, 13, 195, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_ext_atm_td_list, 'aer_bc', p_ext_atm_td%aer_bc, &
      &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,      &
      &            grib2_desc, ldims=shape3d_c, loutput=.FALSE. )


    ! Dust aerosol
    !
    ! aer_dust     p_ext_atm%aer_dust(nproma,nblks_c,ntimes)
    cf_desc    = t_cf_var('aot_dust', '-', &
      &                   'atmosphere absorption optical thickness due '//  &
      &                   'to dust ambient aerosol', DATATYPE_FLT32)
    grib2_desc = t_grib2_var( 0, 13, 193, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_ext_atm_td_list, 'aer_dust', p_ext_atm_td%aer_dust, &
      &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
      &           ldims=shape3d_c, loutput=.FALSE. )


    ! Organic aerosol
    !
    ! aer_org      p_ext_atm%aer_org(nproma,nblks_c,ntimes)
    cf_desc    = t_cf_var('aot_org', '-', &
      &                   'atmosphere absorption optical thickness due '//  &
      &                   'to particulate organic matter ambient aerosol', DATATYPE_FLT32)
    grib2_desc = t_grib2_var( 0, 13, 194, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_ext_atm_td_list, 'aer_org', p_ext_atm_td%aer_org,     &
      &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,&
      &           ldims=shape3d_c, loutput=.FALSE. )


    ! Sulfate aerosol
    !
    ! aer_so4      p_ext_atm%aer_so4(nproma,nblks_c,ntimes)
    cf_desc    = t_cf_var('aot_so4', '-', &
      &                   'atmosphere absorption optical thickness due '//  &
      &                   'to sulfate_ambient_aerosol', DATATYPE_FLT32)
    grib2_desc = t_grib2_var( 0, 13, 192, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_ext_atm_td_list, 'aer_so4', p_ext_atm_td%aer_so4, &
      &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,&
      &           ldims=shape3d_c, loutput=.FALSE. )


    ! Seasalt aerosol
    !
    ! aer_ss       p_ext_atm%aer_ss(nproma,nblks_c,ntimes)
    cf_desc    = t_cf_var('aot_ss', '-', &
      &                   'atmosphere absorption optical thickness due '//  &
      &                   'to seasalt_ambient_aerosol', DATATYPE_FLT32)
    grib2_desc = t_grib2_var( 0, 13, 196, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_ext_atm_td_list, 'aer_ss', p_ext_atm_td%aer_ss, &
      &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,&
      &           ldims=shape3d_c, loutput=.FALSE. )


    !--------------------------------
    ! vegetation parameters
    !--------------------------------

    ! (monthly) proportion of actual value/maximum NDVI
    !
    ! ndvi_mrat     p_ext_atm%ndvi_mrat(nproma,nblks_c,ntimes)
    cf_desc    = t_cf_var('normalized_difference_vegetation_index', '-', &
      &                   '(monthly) proportion of actual value/maximum ' // &
      &                   'normalized differential vegetation index', DATATYPE_FLT32)
    grib2_desc = t_grib2_var( 2, 0, 192, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_ext_atm_td_list, 'ndvi_mrat', p_ext_atm_td%ndvi_mrat, &
      &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,&
      &           ldims=shape3d_c, loutput=.FALSE. )
    !--------------------------------
    !SST and sea ice fraction
    !--------------------------------
    IF ( sstice_mode > 1 ) THEN
     ! sst_m     p_ext_atm_td%sst_m(nproma,nblks_c,ntimes)
     cf_desc    = t_cf_var('sst_m', 'K', &
       &                   '(monthly) sea surface temperature '  &
       &                   , DATATYPE_FLT32)
     grib2_desc = t_grib2_var(192 ,128 , 34, ibits, GRID_REFERENCE, GRID_CELL)
     CALL add_var( p_ext_atm_td_list, 'sst_m', p_ext_atm_td%sst_m, &
       &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,&
       &           ldims=shape3d_sstice, loutput=.FALSE. )

     ! fr_ice_m     p_ext_atm_td%fr_ice_m(nproma,nblks_c,ntimes)
     cf_desc    = t_cf_var('fr_ice_m', '(0-1)', &
       &                   '(monthly) sea ice fraction '  &
       &                   , DATATYPE_FLT32)
     grib2_desc = t_grib2_var( 192,128 ,31 , ibits, GRID_REFERENCE, GRID_CELL)
     CALL add_var( p_ext_atm_td_list, 'fr_ice_m', p_ext_atm_td%fr_ice_m, &
       &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,&
       &           ldims=shape3d_sstice, loutput=.FALSE. )


    ENDIF ! sstice_mode

    ENDIF ! inwp

  END SUBROUTINE new_ext_data_atm_td_list



  !-------------------------------------------------------------------------
  !
  !
  !>
  !! Allocation of oceanic external data structure
  !!
  !! Allocation of oceanic external data structure used for elements that are
  !! stationary in time.
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

    INTEGER :: shape2d_c(2), shape2d_e(2), shape4d_c(4)
    INTEGER :: idim_omip

    INTEGER :: ibits         !< "entropy" of horizontal slice

    !--------------------------------------------------------------

    !determine size of arrays
    nblks_c = p_patch%nblks_c
    nblks_e = p_patch%nblks_e


    ibits = 16   ! "entropy" of horizontal slice

    ! predefined array shapes
    shape2d_c = (/ nproma, nblks_c /)
    shape2d_e = (/ nproma, nblks_e /)

    ! OMIP/NCEP or other flux forcing data on cell centers: 3, 5 or 12 variables, iforc_len data sets
    ! for type of forcing see mo_oce_bulk
    IF (iforc_type == 1 ) idim_omip =  3    !  stress (x, y) and SST
    IF (iforc_type == 2 ) idim_omip = 13    !  OMIP type forcing
    IF (iforc_type == 3 ) idim_omip =  5    !  stress (x, y), SST, net heat and freshwater
    IF (iforc_type == 4 ) idim_omip =  9    !  stress (x, y), SST, and 6 parts of net fluxes
    IF (iforc_type == 5 ) idim_omip = 13    !  NCEP type forcing - time dependent read in mo_oce_bulk
    shape4d_c = (/ nproma, iforc_len, nblks_c, idim_omip /)

    !
    ! Register a field list and apply default settings
    !
    CALL new_var_list( p_ext_oce_list, TRIM(listname), patch_id=p_patch%id )
    CALL default_var_list_settings( p_ext_oce_list,            &
                                  & lrestart=.FALSE.,          &
                                  & restart_type=FILETYPE_NC2, &
                                  & model_type='oce'  )

    ! bathymetric height at cell center
    !
    ! bathymetry_c  p_ext_oce%bathymetry_c(nproma,nblks_c)
    cf_desc    = t_cf_var('Model bathymetry at cell center', 'm', &
      &                   'Model bathymetry', DATATYPE_FLT32)
    grib2_desc = t_grib2_var( 192, 140, 219, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_ext_oce_list, 'bathymetry_c', p_ext_oce%bathymetry_c,      &
      &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d_c )


    ! bathymetric height at cell edge
    !
    ! bathymetry_e  p_ext_oce%bathymetry_e(nproma,nblks_e)
    cf_desc    = t_cf_var('Model bathymetry at edge', 'm', &
      &                   'Model bathymetry', DATATYPE_FLT32)
    grib2_desc = t_grib2_var( 192, 140, 219, ibits, GRID_REFERENCE, GRID_EDGE)
    CALL add_var( p_ext_oce_list, 'bathymetry_e', p_ext_oce%bathymetry_e,      &
      &           GRID_UNSTRUCTURED_EDGE, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d_e )

    ! ocean land-sea-mask at surface on cell centers
    !
    ! lsm_ctr_c  p_ext_oce%lsm_ctr_c(nproma,nblks_c)
    cf_desc    = t_cf_var('Ocean model land-sea-mask at cell center', '-2/-1/1/2', &
      &                   'Ocean model land-sea-mask', DATATYPE_FLT32)
    grib2_desc = t_grib2_var( 192, 140, 219, ibits, GRID_REFERENCE, GRID_CELL)
    !#slo-2011-08-08# does not compile yet?
    CALL add_var( p_ext_oce_list, 'lsm_ctr_c', p_ext_oce%lsm_ctr_c, &
      &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d_c )

    ! ocean land-sea-mask at surface on cell edge
    !
    cf_desc    = t_cf_var('Ocean model land-sea-mask at cell edge', '-2/0/2', &
      &                   'Ocean model land-sea-mask', DATATYPE_FLT32)
    grib2_desc = t_grib2_var( 192, 140, 219, ibits, GRID_REFERENCE, GRID_EDGE)
    CALL add_var( p_ext_oce_list, 'lsm_ctr_e', p_ext_oce%lsm_ctr_e,      &
      &           GRID_UNSTRUCTURED_EDGE, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d_e )

    ! omip forcing data on cell edge
    !
    IF (iforc_oce == 12) THEN
      cf_desc    = t_cf_var('Ocean model OMIP forcing data at cell edge', 'Pa, K', &
        &                   'OMIP forcing data', DATATYPE_FLT32)
      grib2_desc = t_grib2_var( 192, 140, 219, ibits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_ext_oce_list, 'flux_forc_mon_c', p_ext_oce%flux_forc_mon_c,  &
        &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape4d_c )
    END IF

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

    INTEGER :: jg, errstat
    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
      routine = 'mo_ext_data:destruct_ext_data'
    !-------------------------------------------------------------------------

    CALL message (TRIM(routine), 'Destruction of data structure for ' // &
      &                          'external data started')

    IF (iequations/=ihs_ocean) THEN  ! atmosphere model ------------------

      DO jg = 1,n_dom
        ! Delete list of constant in time atmospheric elements
        CALL delete_var_list( ext_data(jg)%atm_list )
      ENDDO

      IF (iforcing > 1 ) THEN
      DO jg = 1,n_dom
        ! Delete list of time-dependent atmospheric elements
        CALL delete_var_list( ext_data(jg)%atm_td_list )
      ENDDO
      END IF

    ELSE ! iequations==ihs_ocean ------------------------------------------

      DO jg = 1,n_dom
        ! Delete list of constant in time oceanic elements
        CALL delete_var_list( ext_data(jg)%oce_list )
      ENDDO

        ! Delete list of time-dependent oceanic elements
        ! ### to be added if necessary ###

    END IF

    DEALLOCATE(nclass_lu, STAT=errstat)
    IF (errstat /= 0)  &
      CALL finish (TRIM(routine), 'Error in DEALLOCATE operation!')

    CALL message (TRIM(routine), 'Destruction of data structure for ' // &
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

    TYPE(t_patch), INTENT(IN)      :: p_patch(:)

    INTEGER :: jg, mpi_comm
    INTEGER :: no_cells, no_verts
    INTEGER :: ncid, dimid

    LOGICAL :: l_exist

    CHARACTER(len=max_char_length), PARAMETER :: &
      routine = 'mo_ext_data:inquire_external_files'

    CHARACTER(filename_max) :: extpar_file !< file name for reading in
    CHARACTER(filename_max) :: ozone_file  !< file name for reading in

!--------------------------------------------------------------------------

    IF(p_test_run) THEN
      mpi_comm = p_comm_work_test
    ELSE
      mpi_comm = p_comm_work
    ENDIF

    DO jg= 1,n_dom

      !------------------------------------------------!
      ! 1. Check validity of external parameter file   !
      !------------------------------------------------!
      IF (itopo == 1 .AND. iequations /=ihs_ocean ) THEN

        IF( my_process_is_stdio()) THEN
          !
          ! generate file name
          extpar_file = generate_filename(extpar_filename,                   &
            &                             model_base_dir,                    &
            &                             TRIM(p_patch(jg)%grid_filename))
          CALL message("read_ext_data_atm, extpar_file=",extpar_file)
        
          INQUIRE (FILE=extpar_file, EXIST=l_exist)
          IF (.NOT.l_exist) THEN
            CALL finish(TRIM(routine),'external data file is not found.')
          ENDIF

          !
          ! open file
          !
          CALL nf(nf_open(TRIM(extpar_file), NF_NOWRITE, ncid), routine)

          !
          ! get number of cells and vertices
          !
          CALL nf(nf_inq_dimid(ncid, 'cell', dimid), routine)
          IF (p_patch(jg)%cell_type == 3) THEN ! triangular grid
            CALL nf(nf_inq_dimlen(ncid, dimid, no_cells), routine)
          ELSEIF (p_patch(jg)%cell_type == 6) THEN ! hexagonal grid
            CALL nf(nf_inq_dimlen(ncid, dimid, no_verts), routine)
          ENDIF

          CALL nf(nf_inq_dimid(ncid, 'vertex', dimid), routine)
          IF (p_patch(jg)%cell_type == 3) THEN ! triangular grid
            CALL nf(nf_inq_dimlen(ncid, dimid, no_verts), routine)
          ELSEIF (p_patch(jg)%cell_type == 6) THEN ! hexagonal grid
            CALL nf(nf_inq_dimlen(ncid, dimid, no_cells), routine)
          ENDIF

          !
          ! check the number of cells and verts
          !
          IF(p_patch(jg)%n_patch_cells_g /= no_cells) THEN
            CALL finish(TRIM(ROUTINE),&
            & 'Number of patch cells and cells in extpar file do not match.')
          ENDIF
          IF(p_patch(jg)%n_patch_verts_g /= no_verts) THEN
            CALL finish(TRIM(ROUTINE),&
            & 'Number of patch verts and verts in extpar file do not match.')
          ENDIF

          !
          ! get the number of landuse classes
          !
          CALL nf(nf_inq_dimid (ncid, 'nclass_lu', dimid), routine)
          CALL nf(nf_inq_dimlen(ncid, dimid, nclass_lu(jg)), routine)


          ! get time dimension from external data file
          CALL nf(nf_inq_dimid (ncid, 'time', dimid), routine)
          CALL nf(nf_inq_dimlen(ncid, dimid, nmonths_ext(jg)), routine)

          WRITE(message_text,'(A,I4)')  &
            & 'Number of months in external data file = ', &
            & nmonths_ext(jg)
          CALL message(TRIM(ROUTINE),message_text)

          !
          ! close file
          !
          CALL nf(nf_close(ncid), routine)

        ENDIF ! my_process_is_stdio()

        ! broadcast nclass_lu from I-Pe to WORK Pes
        CALL p_bcast(nclass_lu(jg), p_io, mpi_comm)

        ! broadcast nmonths from I-Pe to WORK Pes
        CALL p_bcast(nmonths_ext(jg), p_io, mpi_comm)

      ENDIF


      !------------------------------------------------!
      ! 2. Check validity of ozone file                !
      !------------------------------------------------!

      ! default values for nlev_o3 and nmonths
      nlev_o3 = 1
      nmonths   = 1

      O3 : IF((irad_o3 == io3_clim) .OR. (irad_o3 == io3_ape )) THEN

        IF(p_test_run) THEN
          mpi_comm = p_comm_work_test
        ELSE
          mpi_comm = p_comm_work
        ENDIF

        IF(irad_o3 == io3_ape ) THEN
          levelname = 'level'
          zlevelname = 'zlev'
          cellname  = 'ncells'
          o3name    = 'O3'
          o3unit    = 'g/g'
        ELSE ! o3_clim
          levelname= 'plev'
          cellname = 'cell'
          o3name   = 'O3'
          o3unit   = 'g/g' !this unit ozon will have after being read out
                           ! and converted from ppmv
        ENDIF

        IF_IO : IF(my_process_is_stdio()) THEN
          
          WRITE(ozone_file,'(a,i2.2,a)') 'o3_icon_DOM',jg,'.nc'

          ! Note resolution assignment is done per script by symbolic links

          INQUIRE (FILE=ozone_file, EXIST=l_exist)
          IF (.NOT.l_exist) THEN
            WRITE(0,*) 'DOMAIN=',jg
            CALL finish(TRIM(routine),'ozone file of domain is not found.')
          ENDIF

          !
          ! open file
          !
          CALL nf(nf_open(TRIM(ozone_file), NF_NOWRITE, ncid), routine)

          WRITE(0,*)'open ozone file'
          ! get number of cells in triangles and hexagons
          !

          !triangles
          IF (p_patch(jg)%cell_type == 3) THEN ! triangular grid
            CALL nf(nf_inq_dimid (ncid, TRIM(cellname), dimid), routine)
            CALL nf(nf_inq_dimlen(ncid, dimid, no_cells), routine)
  
            WRITE(0,*)'number of cells are', no_cells
          !
          ! check the number of cells and verts
          !
            IF(p_patch(jg)%n_patch_cells_g /= no_cells) THEN
            CALL finish(TRIM(ROUTINE),&
              & 'Number of patch cells and cells in ozone file do not match.')
          ENDIF
        ENDIF
       
          !hexagons
          IF (p_patch(jg)%cell_type == 6) THEN ! hexagonal grid
            CALL nf(nf_inq_dimid (ncid, TRIM(cellname), dimid), routine)
            CALL nf(nf_inq_dimlen(ncid, dimid, no_cells), routine)

            WRITE(0,*)'number of hexcells_o3 are', no_cells
            WRITE(0,*)'number of hexverts are', p_patch(jg)%n_patch_verts_g
            WRITE(0,*)'number of hexcellsmo are', p_patch(jg)%n_patch_cells_g 
          !
          ! check the number of cells and verts
          !
            IF(p_patch(jg)%n_patch_cells_g /= no_cells) THEN
!          IF(p_patch(jg)%n_patch_verts_g /= no_cells) THEN
            CALL finish(TRIM(ROUTINE),&
            & 'Number of patch cells and cells in ozone file do not match.')
          ENDIF
          ENDIF

          ! check the time structure
          CALL nf(nf_inq_dimid (ncid, 'time', dimid), routine)
          CALL nf(nf_inq_dimlen(ncid, dimid, nmonths), routine)

          WRITE(message_text,'(A,I4)')  &
            & 'Number of months in ozone file = ', &
            & nmonths
          CALL message(TRIM(ROUTINE),message_text)

          ! check the vertical structure

!          SELECT CASE (iequations)
!          CASE(ihs_atm_temp,ihs_atm_theta)

            CALL nf(nf_inq_dimid (ncid,TRIM(levelname), dimid), routine)
            CALL nf(nf_inq_dimlen(ncid, dimid, nlev_o3), routine)

            WRITE(message_text,'(A,I4)')  &
              & 'Number of pressure levels in ozone file = ', nlev_o3
            CALL message(TRIM(ROUTINE),message_text)

!          CASE(inh_atmosphere)
!            CALL nf(nf_inq_dimid (ncid,TRIM(zlevelname), dimid), routine)
!            CALL nf(nf_inq_dimlen(ncid, dimid, nlev_o3), routine)
!            WRITE(message_text,'(A,I4)')  &
!              & 'Number of height levels in ozone file = ', nlev_o3
!            CALL message(TRIM(ROUTINE),message_text)
!          END SELECT
  
          !
          ! close file
          !
          CALL nf(nf_close(ncid), routine)

        END IF IF_IO ! pe

        CALL p_bcast(nlev_o3, p_io, mpi_comm)      
        CALL p_bcast(nmonths,   p_io, mpi_comm)      

      END IF O3 !o3

    ENDDO ! ndom

  END SUBROUTINE inquire_external_files

  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !>
  !! Read atmospheric external data
  !!
  !! Read atmospheric external data from netcdf
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2010-07-14)
  !!
  SUBROUTINE read_ext_data_atm (p_patch, ext_data, nlev_o3)

    TYPE(t_patch), INTENT(IN)            :: p_patch(:)
    TYPE(t_external_data), INTENT(INOUT) :: ext_data(:)
    INTEGER, INTENT(IN)                  :: nlev_o3

    CHARACTER(len=max_char_length), PARAMETER :: &
      routine = 'mo_ext_data:read_ext_data_atm'

    CHARACTER(filename_max) :: extpar_file !< file name for reading in
    CHARACTER(filename_max) :: ozone_file  !< file name for reading in
    CHARACTER(len=max_char_length) :: rawdata_attr

    INTEGER :: jg, jc, jb, i, mpi_comm, ilu
    INTEGER :: i_lev,jk
    INTEGER :: ncid, varid

    INTEGER :: rl_start, rl_end
    INTEGER :: i_startblk, i_endblk   !> blocks
    INTEGER :: i_startidx, i_endidx   !< slices
    INTEGER :: i_nchdom               !< domain index

    INTEGER :: i_lctype(n_dom) ! stores the landcover classification used for the external parameter data
                               ! 1: GLC2000, 2: Globcover2009

    REAL(wp):: zdummy_o3lev(nlev_o3) ! will be used for pressure and height levels

    REAL(wp), DIMENSION(num_lcc*n_param_lcc):: lu_glc2000   ! < lookup table landuse class GLC2000
    REAL(wp), DIMENSION(num_lcc*n_param_lcc):: lu_gcv2009   ! < lookup table landuse class GlobCover2009

!                    z0         pcmx      laimx rd      rsmin      snowalb snowtile
!
 DATA lu_glc2000 /   1.00_wp,  0.8_wp,  5.0_wp, 1.0_wp, 250.0_wp,  0.38_wp,-1._wp, & ! evergreen broadleaf forest   
                 &   1.00_wp,  0.9_wp,  6.0_wp, 1.0_wp, 150.0_wp,  0.31_wp,-1._wp, & ! deciduous broadleaf closed forest
                 &   0.15_wp,  0.8_wp,  4.0_wp, 2.0_wp, 150.0_wp,  0.31_wp,-1._wp, & ! deciduous broadleaf open   forest
                 &   1.00_wp,  0.8_wp,  5.0_wp, 0.6_wp, 150.0_wp,  0.27_wp,-1._wp, & ! evergreen needleleaf forest   
                 &   1.00_wp,  0.9_wp,  5.0_wp, 0.6_wp, 150.0_wp,  0.33_wp,-1._wp, & ! deciduous needleleaf forest
                 &   1.00_wp,  0.9_wp,  5.0_wp, 0.8_wp, 150.0_wp,  0.29_wp,-1._wp, & ! mixed leaf trees            
                 &   1.00_wp,  0.8_wp,  5.0_wp, 1.0_wp, 150.0_wp,  -1.0_wp,-1._wp, & ! fresh water flooded trees
                 &   1.00_wp,  0.8_wp,  5.0_wp, 1.0_wp, 150.0_wp,  -1.0_wp,-1._wp, & ! saline water flooded trees
                 &   0.20_wp,  0.8_wp,  2.5_wp, 1.0_wp, 150.0_wp,  -1.0_wp, 1._wp, & ! mosaic tree / natural vegetation
                 &   0.05_wp,  0.5_wp,  0.6_wp, 0.3_wp, 150.0_wp,  -1.0_wp, 1._wp, & ! burnt tree cover
                 &   0.20_wp,  0.8_wp,  3.0_wp, 1.0_wp, 120.0_wp,  -1.0_wp, 1._wp, & ! evergreen shrubs closed-open
                 &   0.15_wp,  0.8_wp,  1.5_wp, 2.0_wp, 120.0_wp,  -1.0_wp, 1._wp, & ! decidous shrubs closed-open
                 &   0.03_wp,  0.9_wp,  3.1_wp, 0.6_wp,  40.0_wp,  -1.0_wp, 1._wp, & ! herbaceous vegetation closed-open
                 &   0.05_wp,  0.5_wp,  0.6_wp, 0.3_wp,  40.0_wp,  -1.0_wp, 1._wp, & ! sparse herbaceous or grass 
                 &   0.05_wp,  0.8_wp,  2.0_wp, 0.4_wp,  40.0_wp,  -1.0_wp,-1._wp, & ! flooded shrubs or herbaceous
                 &   0.07_wp,  0.9_wp,  3.3_wp, 1.0_wp, 120.0_wp,  -1.0_wp, 1._wp, & ! cultivated & managed areas
                 &   0.25_wp,  0.8_wp,  3.0_wp, 1.0_wp, 120.0_wp,  -1.0_wp, 1._wp, & ! mosaic crop / tree / natural vegetation
                 &   0.07_wp,  0.9_wp,  3.5_wp, 1.0_wp, 100.0_wp,  -1.0_wp, 1._wp, & ! mosaic crop / shrub / grass
                 &   0.05_wp,  0.05_wp, 0.6_wp, 0.3_wp, 120.0_wp,  -1.0_wp, 1._wp, & ! bare areas                       
                 &   0.0002_wp,0.0_wp,  0.0_wp, 0.0_wp, 120.0_wp,  -1.0_wp,-1._wp, & ! water
                 &   0.01_wp,  0.0_wp,  0.0_wp, 0.0_wp, 120.0_wp,  -1.0_wp,-1._wp, & ! snow & ice 
                 &   1.00_wp,  0.2_wp,  1.0_wp, 0.6_wp, 120.0_wp,  -1.0_wp,-1._wp, & ! artificial surface  
                 &   0.00_wp,  0.0_wp,  0.0_wp, 0.0_wp,  40.0_wp,  -1.0_wp,-1._wp / ! undefined           

 DATA lu_gcv2009 /   0.07_wp,  0.9_wp,  3.3_wp, 1.0_wp, 120.0_wp,  -1.0_wp, 1._wp, & ! irrigated croplands                           
                 &   0.07_wp,  0.9_wp,  3.3_wp, 1.0_wp, 120.0_wp,  -1.0_wp, 1._wp, & ! rainfed croplands                             
                 &   0.25_wp,  0.8_wp,  3.0_wp, 1.0_wp, 120.0_wp,  -1.0_wp, 1._wp, & ! mosaic cropland (50-70%) - vegetation (20-50%)
                 &   0.07_wp,  0.9_wp,  3.5_wp, 1.0_wp, 100.0_wp,  -1.0_wp, 1._wp, & ! mosaic vegetation (50-70%) - cropland (20-50%)
                 &   1.00_wp,  0.8_wp,  5.0_wp, 1.0_wp, 250.0_wp,  0.38_wp,-1._wp, & ! closed broadleaved evergreen forest           
                 &   1.00_wp,  0.9_wp,  6.0_wp, 1.0_wp, 150.0_wp,  0.31_wp,-1._wp, & ! closed broadleaved deciduous forest           
                 &   0.15_wp,  0.8_wp,  4.0_wp, 2.0_wp, 150.0_wp,  0.31_wp,-1._wp, & ! open broadleaved deciduous forest             
                 &   1.00_wp,  0.8_wp,  5.0_wp, 0.6_wp, 150.0_wp,  0.27_wp,-1._wp, & ! closed needleleaved evergreen forest          
                 &   1.00_wp,  0.9_wp,  5.0_wp, 0.6_wp, 150.0_wp,  0.33_wp,-1._wp, & ! open needleleaved deciduous forest            
                 &   1.00_wp,  0.9_wp,  5.0_wp, 0.8_wp, 150.0_wp,  0.29_wp,-1._wp, & ! mixed broadleaved and needleleaved forest     
                 &   0.20_wp,  0.8_wp,  2.5_wp, 1.0_wp, 150.0_wp,  -1.0_wp, 1._wp, & ! mosaic shrubland (50-70%) - grassland (20-50%)
                 &   0.20_wp,  0.8_wp,  2.5_wp, 1.0_wp, 150.0_wp,  -1.0_wp, 1._wp, & ! mosaic grassland (50-70%) - shrubland (20-50%)
                 &   0.15_wp,  0.8_wp,  2.5_wp, 1.5_wp, 120.0_wp,  -1.0_wp, 1._wp, & ! closed to open shrubland                      
                 &   0.03_wp,  0.9_wp,  3.1_wp, 0.6_wp,  40.0_wp,  -1.0_wp, 1._wp, & ! closed to open herbaceous vegetation          
                 &   0.05_wp,  0.5_wp,  0.6_wp, 0.3_wp,  40.0_wp,  -1.0_wp, 1._wp, & ! sparse vegetation                             
                 &   1.00_wp,  0.8_wp,  5.0_wp, 1.0_wp, 150.0_wp,  -1.0_wp,-1._wp, & ! closed to open forest regulary flooded        
                 &   1.00_wp,  0.8_wp,  5.0_wp, 1.0_wp, 150.0_wp,  -1.0_wp,-1._wp, & ! closed forest or shrubland permanently flooded
                 &   0.05_wp,  0.8_wp,  2.0_wp, 1.0_wp,  40.0_wp,  -1.0_wp,-1._wp, & ! closed to open grassland regularly flooded    
                 &   1.00_wp,  0.2_wp,  1.6_wp, 0.6_wp, 120.0_wp,  -1.0_wp,-1._wp, & ! artificial surfaces                           
                 &   0.05_wp,  0.05_wp, 0.6_wp, 0.3_wp, 120.0_wp,  -1.0_wp, 1._wp, & ! bare areas                                    
                 &   0.0002_wp,0.0_wp,  0.0_wp, 0.0_wp, 120.0_wp,  -1.0_wp,-1._wp, & ! water bodies                                  
                 &   0.01_wp,  0.0_wp,  0.0_wp, 0.0_wp, 120.0_wp,  -1.0_wp,-1._wp, & ! permanent snow and ice                        
                 &   0.00_wp,  0.0_wp,  0.0_wp, 0.0_wp, 250.0_wp,  -1.0_wp,-1._wp  / !undefined                                  





    !----------------------------------------------------------------------

    IF(p_test_run) THEN
      mpi_comm = p_comm_work_test
    ELSE
      mpi_comm = p_comm_work
    ENDIF

    ! Read land-sea mask, sst, land surface albedo, and forest fraction if JSBACH is used

    IF (echam_phy_config%ljsbach) THEN
    DO jg = 1,n_dom

      i_lev = p_patch(jg)%level

      IF(my_process_is_stdio()) &
        CALL nf(nf_open(TRIM(p_patch(jg)%grid_filename), NF_NOWRITE, ncid), routine)

      IF (p_patch(jg)%cell_type == 3) THEN     ! triangular grid

        ! get land-sea-mask on cells, integer marks are:
        ! inner sea (-2), boundary sea (-1, cells and vertices), boundary (0, edges),
        ! boundary land (1, cells and vertices), inner land (2)
        CALL read_netcdf_data (ncid, 'cell_sea_land_mask', p_patch(jg)%n_patch_cells_g, &
          &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
          &                     ext_data(jg)%atm%lsm_ctr_c)

        ! get topography [m]
        CALL read_netcdf_data (ncid, 'cell_elevation', p_patch(jg)%n_patch_cells_g, &
          &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
          &                     ext_data(jg)%atm%elevation_c)
        ext_data(jg)%atm%elevation_c(:,:) = MAX(0._wp, ext_data(jg)%atm%elevation_c(:,:))
        ext_data(jg)%atm%topography_c(:,:) = 0._wp

      ENDIF

      IF( my_process_is_stdio()) CALL nf(nf_close(ncid), routine)

      IF(my_process_is_stdio()) CALL nf(nf_open('sst.nc', NF_NOWRITE, ncid), routine)

      IF (p_patch(jg)%cell_type == 3) THEN     ! triangular grid      

        CALL read_netcdf_data (ncid, 'SST', p_patch(jg)%n_patch_cells_g, &
          &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
          &                     ext_data(jg)%atm%sst)      

      END IF

      IF( my_process_is_stdio()) CALL nf(nf_close(ncid), routine)

      IF(my_process_is_stdio()) CALL nf(nf_open('sst_mon.nc', NF_NOWRITE, ncid), routine)

      IF (p_patch(jg)%cell_type == 3) THEN     ! triangular grid      

        CALL read_netcdf_data (ncid, 'SST', p_patch(jg)%n_patch_cells_g, &
          &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
          &                     12,ext_data(jg)%atm%sst_mon)      

      END IF

      IF( my_process_is_stdio()) CALL nf(nf_close(ncid), routine)

      IF(my_process_is_stdio()) THEN
        !
        ! generate file name and open file
        extpar_file = generate_filename(extpar_filename,                  &
          &                             model_base_dir,                   &
          &                             TRIM(p_patch(jg)%grid_filename))

        CALL nf(nf_open(TRIM(extpar_file), NF_NOWRITE, ncid), routine)
!!$ TR      CALL nf(nf_open('/pf/m/m212070/jsbach_R2B04_v1.nc', NF_NOWRITE, ncid), routine)
      END IF

      IF (p_patch(jg)%cell_type == 3) THEN     ! triangular grid

        ! get land surface albedo on cells
        CALL read_netcdf_data (ncid, 'albedo_soil_vis', p_patch(jg)%n_patch_cells_g, &
          &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
          &                     ext_data(jg)%atm%albedo_vis_soil)
        CALL read_netcdf_data (ncid, 'albedo_soil_nir', p_patch(jg)%n_patch_cells_g, &
          &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
          &                     ext_data(jg)%atm%albedo_nir_soil)
        CALL read_netcdf_data (ncid, 'albedo_veg_vis', p_patch(jg)%n_patch_cells_g, &
          &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
          &                     ext_data(jg)%atm%albedo_vis_canopy)
        CALL read_netcdf_data (ncid, 'albedo_veg_nir', p_patch(jg)%n_patch_cells_g, &
          &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
          &                     ext_data(jg)%atm%albedo_nir_canopy)
        CALL read_netcdf_data (ncid, 'albedo', p_patch(jg)%n_patch_cells_g, &
          &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
          &                     ext_data(jg)%atm%albedo_background)
        CALL read_netcdf_data (ncid, 'forest_fract', p_patch(jg)%n_patch_cells_g, &
          &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
          &                     ext_data(jg)%atm%forest_fract)
      ENDIF

      IF( my_process_is_stdio()) CALL nf(nf_close(ncid), routine)

    END DO
    END IF

    IF(itopo == 1 ) THEN
      DO jg = 1,n_dom

        i_lev = p_patch(jg)%level

        IF(my_process_is_stdio()) THEN

          !
          ! generate file name and open file
          extpar_file = generate_filename(extpar_filename,                  &
            &                             model_base_dir,                   &
            &                             TRIM(p_patch(jg)%grid_filename))
          CALL nf(nf_open(TRIM(extpar_file), NF_NOWRITE, ncid), routine)

          ! Check which data source has been used to generate the external perameters
          CALL nf(nf_get_att_text(ncid, nf_global, 'rawdata', rawdata_attr), routine)

          IF (INDEX(rawdata_attr,'GLC2000') /= 0) THEN
            i_lctype(jg) = 1
          ELSE IF (INDEX(rawdata_attr,'GLOBCOVER2009') /= 0) THEN
            i_lctype(jg) = 2
          ELSE
            CALL finish(TRIM(ROUTINE),'Unknown landcover data source')
          ENDIF

        ENDIF

        ! Broadcast i_lctype from IO-PE to others
        CALL p_bcast(i_lctype(jg), p_io, mpi_comm)

        ! Preset parameter fields with the correct table values
        ilu = 0
        IF (i_lctype(jg) == 1) THEN
          ext_data(jg)%atm%i_lc_snow_ice = 21
          ext_data(jg)%atm%i_lc_water    = 20
          DO i = 1, num_lcc*n_param_lcc, n_param_lcc
            ilu=ilu+1
            ext_data(jg)%atm%z0_lcc(ilu)          = lu_glc2000(i  )  ! Land-cover related roughness length
            ext_data(jg)%atm%plcovmax_lcc(ilu)    = lu_glc2000(i+1)  ! Maximum plant cover fraction for each land-cover class
            ext_data(jg)%atm%laimax_lcc(ilu)      = lu_glc2000(i+2)  ! Maximum leaf area index for each land-cover class
            ext_data(jg)%atm%rootdmax_lcc(ilu)    = lu_glc2000(i+3)  ! Maximum root depth for each land-cover class
            ext_data(jg)%atm%stomresmin_lcc(ilu)  = lu_glc2000(i+4)  ! Minimum stomata resistance for each land-cover class
            ext_data(jg)%atm%snowalb_lcc(ilu)     = lu_glc2000(i+5)  ! Albedo in case of snow cover for each land-cover class
            ext_data(jg)%atm%snowtile_lcc(ilu)    = MERGE(.TRUE.,.FALSE.,lu_glc2000(i+6)>0._wp) ! Existence of snow tiles for land-cover class
          ENDDO
        ELSE IF (i_lctype(jg) == 2) THEN
          i_lctype(jg) = 2
          ext_data(jg)%atm%i_lc_snow_ice = 22
          ext_data(jg)%atm%i_lc_water    = 21
          DO i = 1, num_lcc*n_param_lcc, n_param_lcc
            ilu=ilu+1
            ext_data(jg)%atm%z0_lcc(ilu)          = lu_gcv2009(i  )  ! Land-cover related roughness length
            ext_data(jg)%atm%plcovmax_lcc(ilu)    = lu_gcv2009(i+1)  ! Maximum plant cover fraction for each land-cover class
            ext_data(jg)%atm%laimax_lcc(ilu)      = lu_gcv2009(i+2)  ! Maximum leaf area index for each land-cover class
            ext_data(jg)%atm%rootdmax_lcc(ilu)    = lu_gcv2009(i+3)  ! Maximum root depth for each land-cover class
            ext_data(jg)%atm%stomresmin_lcc(ilu)  = lu_gcv2009(i+4)  ! Minimum stomata resistance for each land-cover class
            ext_data(jg)%atm%snowalb_lcc(ilu)     = lu_gcv2009(i+5)  ! Albedo in case of snow cover for each land-cover class
            ext_data(jg)%atm%snowtile_lcc(ilu)    = MERGE(.TRUE.,.FALSE.,lu_gcv2009(i+6)>0._wp) ! Existence of snow tiles for land-cover class
          ENDDO
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


          SELECT CASE ( iforcing )
          CASE ( inwp )

            CALL read_netcdf_data (ncid, 'PLCOV_MX', p_patch(jg)%n_patch_cells_g,           &
              &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
              &                     ext_data(jg)%atm%plcov_mx)

            CALL read_netcdf_data (ncid, 'LAI_MX', p_patch(jg)%n_patch_cells_g,             &
              &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
              &                     ext_data(jg)%atm%lai_mx)

            CALL read_netcdf_data (ncid, 'ROOTDP', p_patch(jg)%n_patch_cells_g,             &
              &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
              &                     ext_data(jg)%atm%rootdp)

            CALL read_netcdf_data (ncid, 'RSMIN', p_patch(jg)%n_patch_cells_g,              &
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

            CALL read_netcdf_data (ncid, 'SOILTYP', p_patch(jg)%n_patch_cells_g,            &
              &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
              &                     ext_data(jg)%atm%soiltyp)

            CALL read_netcdf_lu (ncid, 'LU_CLASS_FRACTION', p_patch(jg)%n_patch_cells_g,    &
              &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
              &                     nclass_lu(jg), ext_data(jg)%atm%lu_class_fraction )

            CALL read_netcdf_data (ncid, 'ICE', p_patch(jg)%n_patch_cells_g,                &
              &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
              &                     ext_data(jg)%atm%fr_ice)

           
            IF ( l_emiss ) THEN
              CALL read_netcdf_data (ncid, 'EMIS_RAD', p_patch(jg)%n_patch_cells_g,           &
                &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
                &                     ext_data(jg)%atm%emis_rad)
            ELSE
              ext_data(jg)%atm%emis_rad(:,:)= zemiss_def
            ENDIF

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


            !!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! Read time dependent data
            !!!!!!!!!!!!!!!!!!!!!!!!!!!
            IF ( irad_aero == 6 ) THEN

              CALL read_netcdf_data (ncid, nmonths_ext(jg), 'AER_SS', &
                &                    p_patch(jg)%n_patch_cells_g,     &
                &                    p_patch(jg)%n_patch_cells,       &
                &                    p_patch(jg)%cells%glb_index,     & 
                &                    ext_data(jg)%atm_td%aer_ss)

              CALL read_netcdf_data (ncid, nmonths_ext(jg), 'AER_DUST',&
                &                    p_patch(jg)%n_patch_cells_g,      &
                &                    p_patch(jg)%n_patch_cells,        &
                &                    p_patch(jg)%cells%glb_index,      & 
                &                    ext_data(jg)%atm_td%aer_dust)

              CALL read_netcdf_data (ncid, nmonths_ext(jg), 'AER_ORG', &
                &                    p_patch(jg)%n_patch_cells_g,      &
                &                    p_patch(jg)%n_patch_cells,        &
                &                    p_patch(jg)%cells%glb_index,      & 
                &                    ext_data(jg)%atm_td%aer_org)

              CALL read_netcdf_data (ncid, nmonths_ext(jg), 'AER_SO4', &
                &                    p_patch(jg)%n_patch_cells_g,      &
                &                    p_patch(jg)%n_patch_cells,        &
                &                    p_patch(jg)%cells%glb_index,      & 
                &                    ext_data(jg)%atm_td%aer_so4)

              CALL read_netcdf_data (ncid, nmonths_ext(jg), 'AER_BC',  &
                &                    p_patch(jg)%n_patch_cells_g,      &
                &                    p_patch(jg)%n_patch_cells,        &
                &                    p_patch(jg)%cells%glb_index,      & 
                &                    ext_data(jg)%atm_td%aer_bc)

            ENDIF


            CALL read_netcdf_data (ncid, nmonths_ext(jg), 'NDVI_MRAT', &
              &                    p_patch(jg)%n_patch_cells_g,        &
              &                    p_patch(jg)%n_patch_cells,          &
              &                    p_patch(jg)%cells%glb_index,        &
              &                    ext_data(jg)%atm_td%ndvi_mrat)
           

          CASE ( iecham, ildf_echam )

            IF ( l_emiss ) THEN
              CALL read_netcdf_data (ncid, 'EMIS_RAD', p_patch(jg)%n_patch_cells_g,           &
                &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
                &                     ext_data(jg)%atm%emis_rad)
            ELSE
              ext_data(jg)%atm%emis_rad(:,:)= zemiss_def
            ENDIF

          END SELECT ! iforcing



          !
          ! derived external parameter fields
          !

          ! land sea mask at cell centers (LOGICAL)
          !
          i_nchdom  = MAX(1,p_patch(jg)%n_childdom)

          rl_start = 1
          rl_end   = min_rlcell_int

          i_startblk = p_patch(jg)%cells%start_blk(rl_start,1)
          i_endblk   = p_patch(jg)%cells%end_blk(rl_end,i_nchdom)

          DO jb = i_startblk, i_endblk
            CALL get_indices_c(p_patch(jg), jb, i_startblk, i_endblk, &
              &                i_startidx, i_endidx, rl_start, rl_end)

            ! Loop starts with 1 instead of i_startidx 
            ! because the start index is missing in RRTM
            DO jc = 1,i_endidx
              IF (ext_data(jg)%atm%fr_land(jc,jb) > 0.5_wp) THEN
                ext_data(jg)%atm%llsm_atm_c(jc,jb) = .TRUE.  ! land point
              ELSE
                ext_data(jg)%atm%llsm_atm_c(jc,jb) = .FALSE.  ! water point
              ENDIF
            ENDDO
          ENDDO


          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !DR !!!! Quick fix, as long as a routine for computing smoothed external 
          ! parameter fields is missing. Just copy.
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !
          DO jb = i_startblk, i_endblk
            CALL get_indices_c(p_patch(jg), jb, i_startblk, i_endblk, &
              &                i_startidx, i_endidx, rl_start, rl_end)
            ! Loop starts with 1 instead of i_startidx 
            ! because the start index is missing in RRTM
            DO jc = 1,i_endidx
              ext_data(jg)%atm%fr_land_smt(jc,jb) = ext_data(jg)%atm%fr_land(jc,jb)
              ext_data(jg)%atm%fr_ice_smt(jc,jb)  = ext_data(jg)%atm%fr_ice(jc,jb)
            ENDDO

            ! glacier fraction (still missing): QUICK FIX
            DO jc = 1,i_endidx
              IF ( ext_data(jg)%atm%soiltyp(jc,jb) == 1 ) THEN   ! 1: ice
                ext_data(jg)%atm%fr_glac(jc,jb) = 1._wp
              ENDIF
              ext_data(jg)%atm%fr_glac_smt(jc,jb) = ext_data(jg)%atm%fr_glac(jc,jb)
            ENDDO
          ENDDO


        ELSEIF (p_patch(jg)%cell_type == 6) THEN ! hexagonal grid

          CALL finish(TRIM(ROUTINE),&
            & 'Hexagonal grid is not supported, yet.')

        ENDIF

        !
        ! close file
        !
        IF( my_process_is_stdio()) CALL nf(nf_close(ncid), routine)

      ENDDO  ! jg

    ENDIF ! itopo


    !-------------------------------------------------------
    ! Read ozone
    !-------------------------------------------------------

    IF((irad_o3 == io3_clim) .OR. (irad_o3 == io3_ape)) THEN

      IF(p_test_run) THEN
        mpi_comm = p_comm_work_test
      ELSE
        mpi_comm = p_comm_work
      ENDIF

      DO jg = 1,n_dom

        IF(my_process_is_stdio()) THEN
          ! open file
          !
          WRITE(ozone_file,'(a,I2.2,a)') 'o3_icon_DOM',jg,'.nc'
          CALL nf(nf_open(TRIM(ozone_file), NF_NOWRITE, ncid), routine)
          WRITE(0,*)'read ozone levels'

          !          SELECT CASE (iequations)
          !          CASE(ihs_atm_temp,ihs_atm_theta)

          CALL nf(nf_inq_varid(ncid, TRIM(levelname), varid), routine)
          CALL nf(nf_get_var_double(ncid, varid, zdummy_o3lev(:)), routine)

          !          CASE(inh_atmosphere)
          !            CALL nf(nf_inq_varid(ncid, TRIM(zlevelname), varid), routine)
          !            CALL nf(nf_get_var_double(ncid, varid, zdummy_o3lev(:)), routine)
          !          END SELECT

          !
        ENDIF ! pe

        CALL p_bcast(zdummy_o3lev(:), p_io, mpi_comm)      

        !         SELECT CASE (iequations)
        !         CASE(ihs_atm_temp,ihs_atm_theta)

        DO jk=1,nlev_o3
          ext_data(jg)%atm_td%pfoz(jk)=zdummy_o3lev(jk)
        ENDDO

        ! define half levels of ozone pressure grid
        ! upper boundary: ph =      0.Pa -> extrapolation of uppermost value
        ! lower boundary: ph = 125000.Pa -> extrapolation of lowermost value
        ext_data(jg)%atm_td%phoz(1)           = 0._wp
        ext_data(jg)%atm_td%phoz(2:nlev_o3) = (ext_data(jg)%atm_td%pfoz(1:nlev_o3-1) &
          &                                   +  ext_data(jg)%atm_td%pfoz(2:nlev_o3))*.5_wp
        ext_data(jg)%atm_td%phoz(nlev_o3+1) = 125000._wp

        DO i=1,nlev_o3
          WRITE(0,*) 'full/half level press ozone ', i, ext_data(jg)%atm_td%pfoz(i),&
            &                                           ext_data(jg)%atm_td%phoz(i+1)
        ENDDO

        !         CASE(inh_atmosphere)
        !           DO jk=1,nlev_o3
        !             ext_data(jg)%atm_td%zf(jk)=zdummy_o3lev(jk)
        !           ENDDO
        !       END SELECT

        ! we have 2 different ozone files for hexagons and triangels at the moment
        ! therefore ozone data are both stored at the cell centers

        !         IF (p_patch(jg)%cell_type == 3) THEN     ! triangular grid

        CALL read_netcdf_data (ncid, TRIM(o3name), & ! &
          &                    p_patch(jg)%n_patch_cells_g,  &
          &                    p_patch(jg)%n_patch_cells,    &
          &                    p_patch(jg)%cells%glb_index,  & 
          &                    nlev_o3,  nmonths,          &
          &                    ext_data(jg)%atm_td%O3)

        !        ELSEIF (p_patch(jg)%cell_type == 6) THEN ! hexagonal grid
        !
        !          CALL read_netcdf_data (ncid, TRIM(o3name), & 
        !            &                    p_patch(jg)%n_patch_verts_g,  &
        !            &                    p_patch(jg)%n_patch_verts,    & 
        !            &                    p_patch(jg)%verts%glb_index,  &
        !            &                    nlev_o3, nmonths,           &
        !            &                    ext_data(jg)%atm_td%O3)
        !
        !         ENDIF ! patches


        WRITE(0,*)'MAX/min o3 ppmv',MAXVAL(ext_data(jg)%atm_td%O3(:,:,:,:)),&
          &                        MINVAL(ext_data(jg)%atm_td%O3(:,:,:,:))

        ! convert from ppmv to g/g only in case of APE ozone
        IF(irad_o3 == io3_ape) &
          &         ext_data(jg)%atm_td%O3(:,:,:,:)= ext_data(jg)%atm_td%O3(:,:,:,:)*ppmv2gg

        WRITE(0,*)'MAX/min o3 g/g',MAXVAL(ext_data(jg)%atm_td%O3(:,:,:,:)),&
          &                        MINVAL(ext_data(jg)%atm_td%O3(:,:,:,:))

        ! close file
        IF(my_process_is_stdio()) CALL nf(nf_close(ncid), routine)

      ENDDO ! ndom
    END IF ! irad_o3

!
! Read time dependent SST and ICE Fraction  
!
   IF (sstice_mode == 2) THEN

      IF(p_test_run) THEN
        mpi_comm = p_comm_work_test
      ELSE
        mpi_comm = p_comm_work
      ENDIF

      DO jg = 1,n_dom
       !Read the climatological values for SST and ice cover

        CALL finish(TRIM(routine),'sstice_mode == 2 not yet implemented .')
      END DO ! ndom

   END IF ! sstice_mode


  END SUBROUTINE read_ext_data_atm
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !>
  !! Read ocean external data
  !!
  !! Read ocean external data from netcdf
  !!
  !! @par Revision History
  !! Initial revision by Stephan Lorenz, MPI (2011-06-17)
  !!
  SUBROUTINE read_ext_data_oce (p_patch, ext_data)

    TYPE(t_patch), INTENT(IN)            :: p_patch(:)
    TYPE(t_external_data), INTENT(INOUT) :: ext_data(:)

    CHARACTER(len=max_char_length), PARAMETER :: &
      routine = 'mo_ext_data:read_ext_data_oce'

    CHARACTER(filename_max) :: grid_file   !< file name for reading in
    CHARACTER(filename_max) :: omip_file   !< file name for reading in

    LOGICAL :: l_exist
    INTEGER :: jg, i_lev, i_cell_type, no_cells, no_verts, no_tst
    INTEGER :: ncid, dimid
    INTEGER :: mpi_comm

    !REAL(wp):: z_flux(nproma, 12,p_patch(1)%nblks_c)
    REAL(wp):: z_flux(nproma,iforc_len,p_patch(1)%nblks_c)
    TYPE (t_keyword_list), POINTER :: keywords => NULL()

!-------------------------------------------------------------------------

    CALL message (TRIM(routine), 'start')

!-------------------------------------------------------------------------

    !  READ OCEAN BATHYMETRY

    !-------------------------------------------------------------------------

    !DO jg = 1,n_dom
    jg = 1

    i_lev       = p_patch(jg)%level
    i_cell_type = p_patch(jg)%cell_type

    IF(my_process_is_stdio()) THEN
      !
      ! bathymetry and lsm are read from the general ICON grid file
      ! bathymetry and lsm integrated into grid file by grid generator
      !WRITE (bathy_file,'(a,i0,a,i2.2,a)') 'iconR',nroot,'B',i_lev, '-grid.nc'
      !write(*,*) 'bathy_file = ',TRIM(bathy_file)
      !write(*,*) 'dynamics_grid_filename = ', TRIM(dynamics_grid_filename(1))
      CALL associate_keyword("<path>", TRIM(model_base_dir), keywords)
      grid_file = TRIM(with_keywords(keywords, TRIM(dynamics_grid_filename(jg))))

      INQUIRE (FILE=grid_file, EXIST=l_exist)
      IF (.NOT.l_exist) THEN
        CALL finish(TRIM(routine),'Grid file for reading bathymetry not found.')
      ENDIF

      !
      ! open file
      !
      CALL nf(nf_open(TRIM(grid_file), NF_NOWRITE, ncid), routine)

      !
      ! get number of cells and vertices
      !
      CALL nf(nf_inq_dimid(ncid, 'cell', dimid), routine)
      IF (i_cell_type == 3) THEN ! triangular grid
        CALL nf(nf_inq_dimlen(ncid, dimid, no_cells), routine)
      ELSEIF (i_cell_type == 6) THEN ! hexagonal grid
        CALL nf(nf_inq_dimlen(ncid, dimid, no_verts), routine)
      ENDIF

      CALL nf(nf_inq_dimid(ncid, 'vertex', dimid), routine)
      IF (i_cell_type == 3) THEN ! triangular grid
        CALL nf(nf_inq_dimlen(ncid, dimid, no_verts), routine)
      ELSEIF (i_cell_type == 6) THEN ! hexagonal grid
        CALL nf(nf_inq_dimlen(ncid, dimid, no_cells), routine)
      ENDIF

      !
      ! check the number of cells and verts
      !
      IF(p_patch(jg)%n_patch_cells_g /= no_cells) THEN
        CALL finish(TRIM(ROUTINE),&
        & 'Number of patch cells and cells in bathymetry file do not match.')
      ENDIF
      IF(p_patch(jg)%n_patch_verts_g /= no_verts) THEN
        CALL finish(TRIM(ROUTINE),&
        & 'Number of patch verts and verts in bathymetry file do not match.')
      ENDIF

      WRITE(message_text,'(3(a,i6))') 'No of cells =', no_cells, &
        &                           '  no of edges =', p_patch(jg)%n_patch_edges_g, &
        &                           '  no of verts =', no_verts
      CALL message( TRIM(routine),TRIM(message_text))
    ENDIF


    !-------------------------------------------------------
    !
    ! Read bathymetry for triangle centers and edges
    !
    !-------------------------------------------------------

    ! triangle center and edges

    IF (i_cell_type == 3) THEN     ! triangular grid

      ! These arrays are not included in standard icon-grid, but they are
      ! created by "create_ocean_grid"
      CALL read_netcdf_data (ncid, 'cell_elevation', p_patch(jg)%n_patch_cells_g,     &
        &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
        &                     ext_data(jg)%oce%bathymetry_c)

      CALL read_netcdf_data (ncid, 'edge_elevation', p_patch(jg)%n_patch_edges_g,     &
        &                     p_patch(jg)%n_patch_edges, p_patch(jg)%edges%glb_index, &
        &                     ext_data(jg)%oce%bathymetry_e)

      ! get land-sea-mask on cells, integer marks are:
      ! inner sea (-2), boundary sea (-1, cells and vertices), boundary (0, edges),
      ! boundary land (1, cells and vertices), inner land (2)
      CALL read_netcdf_data (ncid, 'cell_sea_land_mask', p_patch(jg)%n_patch_cells_g, &
        &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
        &                     ext_data(jg)%oce%lsm_ctr_c)

      CALL read_netcdf_data (ncid, 'edge_sea_land_mask', p_patch(jg)%n_patch_edges_g, &
        &                    p_patch(jg)%n_patch_edges, p_patch(jg)%edges%glb_index,  &
        &                    ext_data(jg)%oce%lsm_ctr_e)

    ENDIF

    !
    ! close file
    !
    IF(my_process_is_stdio()) CALL nf(nf_close(ncid), routine)

    !ENDDO ! jg

    CALL message( TRIM(routine),'Ocean bathymetry for external data read' )

    !-------------------------------------------------------------------------

    !  READ OMIP FORCING

    !-------------------------------------------------------------------------

    IF (iforc_type .NE. 5 .AND. iforc_oce == 12) THEN

    !DO jg = 1,n_dom
      jg = 1

      i_lev       = p_patch(jg)%level
      i_cell_type = p_patch(jg)%cell_type

      IF(my_process_is_stdio()) THEN
        !
        WRITE (omip_file,'(a,i0,a,i2.2,a)') 'iconR',nroot,'B',i_lev, '-flux.nc'

        !omip_file=TRIM('/pool/data/ICON/external/iconR2B04-flux.nc')
        !omip_file='/scratch/local1/m212053/ICON/trunk/icon-dev/grids/omip4icon-R2B02-monmean.nc'
        CALL message( TRIM(routine),'Ocean OMIP forcing flux file is: '//TRIM(omip_file) )
        INQUIRE (FILE=omip_file, EXIST=l_exist)
        IF (.NOT.l_exist) THEN
          write(*,*)'FORCING FILE: ',TRIM(omip_file)
          CALL finish(TRIM(routine),'OMIP forcing flux file is not found - ABORT')
        ENDIF

        !
        ! open file
        !
        CALL nf(nf_open(TRIM(omip_file), NF_NOWRITE, ncid), routine)
        CALL message( TRIM(routine),'Ocean OMIP flux file opened for read' )

        !
        ! get and check number of cells in OMIP data
        !
     !  CALL nf(nf_inq_dimid(ncid, 'cell', dimid), routine)  !  workaround for r2b2 omip.daily
        CALL nf(nf_inq_dimid (ncid, 'ncells', dimid), routine)
        CALL nf(nf_inq_dimlen(ncid, dimid, no_cells), routine)

        IF(p_patch(jg)%n_patch_cells_g /= no_cells) THEN
          CALL finish(TRIM(ROUTINE),&
          & 'Number of patch cells and cells in OMIP flux file do not match - ABORT')
        ENDIF

        !
        ! get number of timesteps
        !
        CALL nf(nf_inq_dimid (ncid, 'time', dimid), routine)
        CALL nf(nf_inq_dimlen(ncid, dimid, no_tst), routine)
        !
        ! check
        !
        WRITE(message_text,'(A,I6,A)')  'Ocean OMIP flux file contains',no_tst,' data sets'
        CALL message( TRIM(routine), TRIM(message_text) )
        IF(no_tst /= iforc_len ) THEN
          CALL finish(TRIM(ROUTINE),&
          & 'Number of forcing timesteps is not equal iforc_len specified in namelist - ABORT')
        ENDIF
      ENDIF
      IF(p_test_run) THEN
        mpi_comm = p_comm_work_test
      ELSE
        mpi_comm = p_comm_work
      ENDIF
      CALL p_bcast(no_tst, p_io, mpi_comm)

      !-------------------------------------------------------
      !
      ! Read complete OMIP data for triangle centers
      !
      !-------------------------------------------------------

      ! provide OMIP fluxes for wind stress forcing
      ! 1:  'stress_x': zonal wind stress       [m/s]
      ! 2:  'stress_y': meridional wind stress  [m/s]
      ! 3:  'SST"     : sea surface temperature [K]

      ! zonal wind stress
      CALL read_netcdf_data (ncid, 'stress_x', p_patch(jg)%n_patch_cells_g,          &
        &                    p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
        &                    no_tst, z_flux(:,:,:))
      ext_data(jg)%oce%flux_forc_mon_c(:,:,:,1) = z_flux(:,:,:)

      ! meridional wind stress
      CALL read_netcdf_data (ncid, 'stress_y', p_patch(jg)%n_patch_cells_g,          &
        &                    p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
        &                    no_tst, z_flux)
      ext_data(jg)%oce%flux_forc_mon_c(:,:,:,2) = z_flux(:,:,:)

      ! SST
      CALL read_netcdf_data (ncid, 'SST', p_patch(jg)%n_patch_cells_g,           &
        &                    p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
        &                    no_tst, z_flux)
      ext_data(jg)%oce%flux_forc_mon_c(:,:,:,3) = z_flux(:,:,:)

      IF (iforc_type == 2) THEN

      ! Read complete OMIP data sets for focing ocean model
      ! 4:  tafo(:,:),   &  ! 2 m air temperature                              [C]
      ! 5:  ftdew(:,:),  &  ! 2 m dew-point temperature                        [K]
      ! 6:  fu10(:,:) ,  &  ! 10 m wind speed                                  [m/s]
      ! 7:  fclou(:,:),  &  ! Fractional cloud cover
      ! 8:  pao(:,:),    &  ! Surface atmospheric pressure                     [hPa]
      ! 9:  fswr(:,:),   &  ! Incoming surface solar radiation                 [W/m]

        ! 2m-temperature
        CALL read_netcdf_data (ncid, 'temp_2m', p_patch(jg)%n_patch_cells_g,           &
          &                    p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
          &                    no_tst, z_flux)
        ext_data(jg)%oce%flux_forc_mon_c(:,:,:,4) = z_flux(:,:,:)
     
        ! 2m dewpoint temperature
        CALL read_netcdf_data (ncid, 'dpt_temp_2m', p_patch(jg)%n_patch_cells_g,       &
          &                    p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
          &                    no_tst, z_flux)
        ext_data(jg)%oce%flux_forc_mon_c(:,:,:,5) = z_flux(:,:,:)
     
        ! Scalar wind
        CALL read_netcdf_data (ncid, 'scalar_wind', p_patch(jg)%n_patch_cells_g,       &
          &                    p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
          &                    no_tst, z_flux)
        ext_data(jg)%oce%flux_forc_mon_c(:,:,:,6) = z_flux(:,:,:)
     
        ! cloud cover
        CALL read_netcdf_data (ncid, 'cloud', p_patch(jg)%n_patch_cells_g,             &
          &                    p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
          &                    no_tst, z_flux)
        ext_data(jg)%oce%flux_forc_mon_c(:,:,:,7) = z_flux(:,:,:)
     
        ! sea level pressure
        CALL read_netcdf_data (ncid, 'pressure', p_patch(jg)%n_patch_cells_g,          &
          &                    p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
          &                    no_tst, z_flux)
        ext_data(jg)%oce%flux_forc_mon_c(:,:,:,8) = z_flux(:,:,:)
     
        ! total solar radiation
        CALL read_netcdf_data (ncid, 'tot_solar', p_patch(jg)%n_patch_cells_g,         &
          &                    p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
          &                    no_tst, z_flux)
        ext_data(jg)%oce%flux_forc_mon_c(:,:,:,9) = z_flux(:,:,:)
     
        ! precipitation
        CALL read_netcdf_data (ncid, 'precip', p_patch(jg)%n_patch_cells_g,            &
          &                    p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
          &                    no_tst, z_flux)
        ext_data(jg)%oce%flux_forc_mon_c(:,:,:,10) = z_flux(:,:,:)
     
        ! evaporation
        CALL read_netcdf_data (ncid, 'evap', p_patch(jg)%n_patch_cells_g,              &
          &                    p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
          &                    no_tst, z_flux)
        ext_data(jg)%oce%flux_forc_mon_c(:,:,:,11) = z_flux(:,:,:)
     
        ! runoff
        CALL read_netcdf_data (ncid, 'runoff', p_patch(jg)%n_patch_cells_g,            &
          &                    p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
          &                    no_tst, z_flux)
        ext_data(jg)%oce%flux_forc_mon_c(:,:,:,12) = z_flux(:,:,:)

      END IF

      ! provide heat and freshwater flux for focing ocean model
      IF (iforc_type == 3) THEN

        ! net surface heat flux
        CALL read_netcdf_data (ncid, 'net_hflx', p_patch(jg)%n_patch_cells_g,          &
          &                    p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
          &                    no_tst, z_flux)
        ext_data(jg)%oce%flux_forc_mon_c(:,:,:,4) = z_flux(:,:,:)
     
        ! surface freshwater flux
        CALL read_netcdf_data (ncid, 'net_fflx', p_patch(jg)%n_patch_cells_g,          &
          &                    p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
          &                    no_tst, z_flux)
        ext_data(jg)%oce%flux_forc_mon_c(:,:,:,5) = z_flux(:,:,:)

      END IF

      ! provide 4 parts of heat and 2 parts of freshwater flux for focing ocean model
      IF (iforc_type == 4) THEN

        ! surface short wave heat flux
        CALL read_netcdf_data (ncid, 'swflxsfc_avg', p_patch(jg)%n_patch_cells_g,      &
          &                    p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
          &                    no_tst, z_flux)
        ext_data(jg)%oce%flux_forc_mon_c(:,:,:,4) = z_flux(:,:,:)
     
        ! surface long wave heat flux
        CALL read_netcdf_data (ncid, 'lwflxsfc_avg', p_patch(jg)%n_patch_cells_g,      &
          &                    p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
          &                    no_tst, z_flux)
        ext_data(jg)%oce%flux_forc_mon_c(:,:,:,5) = z_flux(:,:,:)
     
        ! surface sensible heat flux
        CALL read_netcdf_data (ncid, 'shflx_avg',    p_patch(jg)%n_patch_cells_g,      &
          &                    p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
          &                    no_tst, z_flux)
        ext_data(jg)%oce%flux_forc_mon_c(:,:,:,6) = z_flux(:,:,:)
     
        ! surface latent heat flux
        CALL read_netcdf_data (ncid, 'lhflx_avg',    p_patch(jg)%n_patch_cells_g,      &
          &                    p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
          &                    no_tst, z_flux)
        ext_data(jg)%oce%flux_forc_mon_c(:,:,:,7) = z_flux(:,:,:)
     
        ! total precipiation
        CALL read_netcdf_data (ncid, 'precip', p_patch(jg)%n_patch_cells_g,            &
          &                    p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
          &                    no_tst, z_flux)
        ext_data(jg)%oce%flux_forc_mon_c(:,:,:,8) = z_flux(:,:,:)
     
        ! evaporation
        CALL read_netcdf_data (ncid, 'evap'  , p_patch(jg)%n_patch_cells_g,            &
          &                    p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
          &                    no_tst, z_flux)
        ext_data(jg)%oce%flux_forc_mon_c(:,:,:,9) = z_flux(:,:,:)

      END IF

      !
      ! close file
      !
      IF(my_process_is_stdio()) CALL nf(nf_close(ncid), routine)

    !ENDDO

      CALL message( TRIM(routine),'Ocean OMIP fluxes for external data read' )

    END IF ! iforc_oce=12 and iforc_type.ne.5

  END SUBROUTINE read_ext_data_oce
  !-------------------------------------------------------------------------



  SUBROUTINE init_index_lists (p_patch, ext_data)

    TYPE(t_patch), INTENT(IN)            :: p_patch(:)
    TYPE(t_external_data), INTENT(INOUT) :: ext_data(:)

    INTEGER :: i_lu, jb,jc, jg, n_lu, i_count, i_count_flk, ic, jt, jt_in
    INTEGER :: i_count_sea             ! number of sea points
    INTEGER :: rl_start, rl_end
    INTEGER :: i_startblk, i_endblk    !> blocks
    INTEGER :: i_startidx, i_endidx    !< slices
    INTEGER :: i_nchdom                !< domain index
    LOGICAL  :: tile_mask(num_lcc) = .true. 
    REAL(wp) :: tile_frac(num_lcc), sum_frac
    INTEGER  :: lu_subs, it_count(ntiles_total)
    INTEGER  :: npoints, npoints_sea, npoints_lake
    INTEGER  :: i_lc_water

    REAL(wp), POINTER  ::  &  !< pointer to proportion of actual value/maximum
      &  ptr_ndvi_mrat(:,:)   !< NDVI (for starting time of model integration)

    CHARACTER(len=max_char_length), PARAMETER :: &
      routine = 'mo_ext_data:init_index_lists'
    !-------------------------------------------------------------------------

    WRITE(message_text,'(a,i4)') &
      &  'Index list generation - number of land tiles: ', ntiles_lnd
    CALL message('', TRIM(message_text))
    WRITE(message_text,'(a,i4)')  'Total number of tiles: ', ntiles_total
    CALL message('', TRIM(message_text))


    DO jg = 1, n_dom 

       n_lu = nclass_lu(jg)

       ptr_ndvi_mrat => ext_data(jg)%atm%ndvi_mrat(:,:)

       i_nchdom  = MAX(1,p_patch(jg)%n_childdom)

       ! exclude the boundary interpolation zone of nested domains
       rl_start = grf_bdywidth_c+1
       rl_end   = min_rlcell_int

       i_startblk = p_patch(jg)%cells%start_blk(rl_start,1)
       i_endblk   = p_patch(jg)%cells%end_blk(rl_end,i_nchdom)
     
       i_lc_water = ext_data(jg)%atm%i_lc_water

       ! Initialization of index list counts - moved here in order to avoid uninitialized elements
       ! along nest boundaries
       ext_data(jg)%atm%lp_count(:) = 0
       ext_data(jg)%atm%sp_count(:) = 0
       ext_data(jg)%atm%fp_count(:) = 0

       ext_data(jg)%atm%spi_count(:) = 0
       ext_data(jg)%atm%spw_count(:) = 0  

       ext_data(jg)%atm%gp_count_t(:,:) = 0
       ext_data(jg)%atm%lp_count_t(:,:) = 0
      

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,i_lu,i_startidx,i_endidx,i_count,i_count_sea,i_count_flk,tile_frac,&
!$OMP            tile_mask,lu_subs,sum_frac,it_count,ic,jt,jt_in ) ICON_OMP_DEFAULT_SCHEDULE
       DO jb=i_startblk, i_endblk

         CALL get_indices_c(p_patch(jg), jb, i_startblk, i_endblk, &
            & i_startidx, i_endidx, rl_start, rl_end)

         i_count                       = 0   ! counter for land points
         i_count_sea                   = 0   ! counter for sea points
         i_count_flk                   = 0   ! counter for lake points

         it_count(:)                   = 0 ! counter for tiles

         DO jc = i_startidx, i_endidx
           IF (ext_data(jg)%atm%fr_land(jc,jb)> frlnd_thrhld) THEN ! searching for land-points 
             i_count=i_count+1
             ext_data(jg)%atm%idx_lst_lp(i_count,jb) = jc  ! write index of land-points

             tile_frac(:)= ext_data(jg)%atm%lu_class_fraction(jc,jb,:)
             tile_mask(:)=.true.
             tile_mask(i_lc_water)=.false. ! exclude water points
             ext_data(jg)%atm%lc_class_t(jc,jb,:) = -1    ! dummy value for undefined points

             ext_data(jg)%atm%lp_count(jb) = i_count

             IF (ntiles_lnd == 1) THEN 

               ! i_lu=1 contains grid-box mean values from EXTPAR!
               ext_data(jg)%atm%rootdp_t (jc,jb,1)  = ext_data(jg)%atm%rootdp(jc,jb)
               ext_data(jg)%atm%plcov_t  (jc,jb,1)  = ptr_ndvi_mrat(jc,jb)  &
                 &                                       * ext_data(jg)%atm%plcov_mx(jc,jb)
               ext_data(jg)%atm%tai_t    (jc,jb,1)  = ptr_ndvi_mrat(jc,jb)**2  &
                 & * ext_data(jg)%atm%plcov_mx(jc,jb)*ext_data(jg)%atm%lai_mx(jc,jb)
               ext_data(jg)%atm%sai_t    (jc,jb,1)  = c_lnd+ext_data(jg)%atm%tai_t(jc,jb,1)
               ext_data(jg)%atm%eai_t    (jc,jb,1)  = c_soil      
               ext_data(jg)%atm%rsmin2d_t(jc,jb,1)  = ext_data(jg)%atm%rsmin(jc,jb)
               ext_data(jg)%atm%soiltyp_t(jc,jb,1)  = ext_data(jg)%atm%soiltyp(jc,jb)
               ext_data(jg)%atm%lc_frac_t(jc,jb,1)  = 1._wp
               ext_data(jg)%atm%lc_class_t(jc,jb,1) = MAXLOC(tile_frac,1,tile_mask)
               !  Workaround for GLC2000 hole below 60 deg S
               IF (ext_data(jg)%atm%lc_class_t(jc,jb,1) <= 0) &
                 ext_data(jg)%atm%lc_class_t(jc,jb,1) = ext_data(jg)%atm%i_lc_snow_ice

               ! static index list and corresponding counter
               ext_data(jg)%atm%idx_lst_lp_t(i_count,jb,1)  = jc
               ext_data(jg)%atm%lp_count_t(jb,1)            = i_count

               ! initialize dynamic index list (in case of lsnowtile=true) with the same values
               ext_data(jg)%atm%idx_lst_t(i_count,jb,1) = jc
               ext_data(jg)%atm%gp_count_t(jb,1)        = i_count

               ! initialize snowtile flag with 1 if the tile is eligible for separate treatment of
               ! a snow-covered and a snow-free part, otherwise with -1
               IF (ext_data(jg)%atm%snowtile_lcc(ext_data(jg)%atm%lc_class_t(jc,jb,1))) THEN
                 ext_data(jg)%atm%snowtile_flag_t(jc,jb,1)  = 1
               ELSE
                 ext_data(jg)%atm%snowtile_flag_t(jc,jb,1)  = -1
               ENDIF

             ELSE    
               ext_data(jg)%atm%lc_frac_t(jc,jb,:)  = 0._wp ! to be really safe

               DO i_lu = 1, ntiles_lnd
                 lu_subs = MAXLOC(tile_frac,1,tile_mask)
                 IF (tile_frac(lu_subs) >= frlndtile_thrhld) THEN
                   it_count(i_lu)    = it_count(i_lu) + 1
                   tile_mask(lu_subs)= .FALSE.

                   ! static index list and corresponding counter
                   ext_data(jg)%atm%idx_lst_lp_t(it_count(i_lu),jb,i_lu) = jc
                   ext_data(jg)%atm%lp_count_t(jb,i_lu)                  = it_count(i_lu)

                   ! initialize dynamic index list (in case of lsnowtile=true) with the same values
                   ext_data(jg)%atm%idx_lst_t(it_count(i_lu),jb,i_lu) = jc
                   ext_data(jg)%atm%gp_count_t(jb,i_lu)               = it_count(i_lu)

                   ! initialize snowtile flag with 1 if the tile is eligible for separate treatment of
                   ! a snow-covered and a snow-free part, otherwise with -1
                   IF (ext_data(jg)%atm%snowtile_lcc(lu_subs)) THEN
                     ext_data(jg)%atm%snowtile_flag_t(jc,jb,i_lu)  = 1
                   ELSE
                     ext_data(jg)%atm%snowtile_flag_t(jc,jb,i_lu)  = -1
                   ENDIF

                   ext_data(jg)%atm%lc_frac_t(jc,jb,i_lu)  = tile_frac(lu_subs)
                   ext_data(jg)%atm%lc_class_t(jc,jb,i_lu) = lu_subs
                 ELSE
                   EXIT ! no more land cover classes exceeding the threshold
                 ENDIF

               END DO

               IF (ext_data(jg)%atm%fr_land(jc,jb) < 0.5_wp) THEN
                 ! fix for non-dominant land points: reset soil type to sandy loam ...
                 ext_data(jg)%atm%soiltyp(jc,jb) = 4
                 ! ... and reset ndvi_mrat to 0.5
                 ptr_ndvi_mrat(jc,jb) = 0.5_wp
               ENDIF

               sum_frac = SUM(ext_data(jg)%atm%lc_frac_t(jc,jb,1:ntiles_lnd))

               DO i_lu = 1, ntiles_lnd 

                 IF ( sum_frac < 1.e-10_wp) THEN !  Workaround for GLC2000 hole below 60 deg S
                   IF (i_lu == 1) THEN
                     it_count(i_lu)    = it_count(i_lu) + 1
                     ! static index list and corresponding counter
                     ext_data(jg)%atm%idx_lst_lp_t(it_count(i_lu),jb,i_lu) = jc
                     ext_data(jg)%atm%lp_count_t(jb,i_lu)               = it_count(i_lu)

                     ! initialize dynamic index list (in case of lsnowtile=true) with the same values
                     ext_data(jg)%atm%idx_lst_t(it_count(i_lu),jb,i_lu) = jc
                     ext_data(jg)%atm%gp_count_t(jb,i_lu)               = it_count(i_lu)

                     ! the snowtile flag is initialized with -1 here because the snow/ice class is not
                     ! eligible for separate consideration of a snow-free and a snow-covered part
                     ext_data(jg)%atm%snowtile_flag_t(jc,jb,i_lu)         = -1

                     ext_data(jg)%atm%lc_class_t(jc,jb,i_lu) = ext_data(jg)%atm%i_lc_snow_ice
                     ext_data(jg)%atm%lc_frac_t(jc,jb,i_lu)  = ext_data(jg)%atm%fr_land(jc,jb)
                   ELSE
                     ext_data(jg)%atm%lc_class_t(jc,jb,i_lu) = -1
                     ext_data(jg)%atm%lc_frac_t(jc,jb,i_lu)  = 0._wp
                   ENDIF 
                 END IF

                 lu_subs = ext_data(jg)%atm%lc_class_t(jc,jb,i_lu)
                 IF (lu_subs < 0) CYCLE

                 ! root depth
                 ext_data(jg)%atm%rootdp_t (jc,jb,i_lu)  = ext_data(jg)%atm%rootdmax_lcc(lu_subs)
                 ! plant cover
                 ext_data(jg)%atm%plcov_t  (jc,jb,i_lu)  =    &
                   &         ptr_ndvi_mrat(jc,jb) * ext_data(jg)%atm%plcovmax_lcc(lu_subs)
                 ! max leaf area index
                 ext_data(jg)%atm%tai_t    (jc,jb,i_lu)  =    &
                   &  ptr_ndvi_mrat(jc,jb)**2 *ext_data(jg)%atm%plcovmax_lcc(lu_subs)  &
                   &     * ext_data(jg)%atm%laimax_lcc(lu_subs)

                 ! max leaf area index
                 ext_data(jg)%atm%sai_t    (jc,jb,i_lu)  = c_lnd+ ext_data(jg)%atm%tai_t (jc,jb,i_lu)

                 ! max leaf area index
                 ext_data(jg)%atm%eai_t    (jc,jb,i_lu)  = c_soil

                 ! minimal stomata resistence
                 ext_data(jg)%atm%rsmin2d_t(jc,jb,i_lu)  = ext_data(jg)%atm%stomresmin_lcc(lu_subs)
    
                 ! soil type
                 ext_data(jg)%atm%soiltyp_t(jc,jb,i_lu)  = ext_data(jg)%atm%soiltyp(jc,jb)
               END DO
             END IF ! ntiles
           ENDIF

           jt = ntiles_total + MIN(1,ntiles_water)
           !
           ! searching for lake-points
           !
           IF (ext_data(jg)%atm%fr_lake(jc,jb) >= frlake_thrhld) THEN 
             i_count_flk=i_count_flk+1
             ext_data(jg)%atm%idx_lst_fp(i_count_flk,jb) = jc  ! write index of lake-points
             ext_data(jg)%atm%fp_count(jb) = i_count_flk
             ! set land-cover class
             ext_data(jg)%atm%lc_class_t(jc,jb,jt)  = ext_data(jg)%atm%i_lc_water
             ! set also area fractions
             ext_data(jg)%atm%lc_frac_t(jc,jb,jt)  = ext_data(jg)%atm%fr_lake(jc,jb)
           ELSE IF (1._wp-ext_data(jg)%atm%fr_land(jc,jb) >= frsea_thrhld) THEN
             ! 
             ! searching for sea points 
             !
             i_count_sea=i_count_sea + 1
             ext_data(jg)%atm%idx_lst_sp(i_count_sea,jb) = jc  ! write index of sea-points
             ext_data(jg)%atm%sp_count(jb) = i_count_sea
             ! set land-cover class
             ext_data(jg)%atm%lc_class_t(jc,jb,jt) = ext_data(jg)%atm%i_lc_water
             ! set also area fractions
             ext_data(jg)%atm%lc_frac_t(jc,jb,jt)  = 1._wp-ext_data(jg)%atm%fr_land(jc,jb)
           ENDIF

         END DO ! jc



         ! Normalize land-cover fractions over the land tile indices plus the first water index to a sum of 1
         IF (ntiles_lnd > 1) THEN
           DO jc = i_startidx, i_endidx
             sum_frac = SUM(ext_data(jg)%atm%lc_frac_t(jc,jb,1:ntiles_lnd)) + &
                            ext_data(jg)%atm%lc_frac_t(jc,jb,ntiles_total+MIN(1,ntiles_water))

             DO jt = 1, ntiles_total + MIN(1,ntiles_water)
               ext_data(jg)%atm%lc_frac_t(jc,jb,jt) = ext_data(jg)%atm%lc_frac_t(jc,jb,jt) / sum_frac
             ENDDO
           ENDDO
         ELSE ! overwrite fractional settings over water points if tile approach is turned off
           DO jc = i_startidx, i_endidx
             ext_data(jg)%atm%lc_frac_t(jc,jb,1) = 1
           ENDDO
         ENDIF

         IF (lsnowtile) THEN ! copy static external data fields to snow tile grid points
           DO jt = ntiles_lnd+1, ntiles_total

             jt_in = jt - ntiles_lnd
             ext_data(jg)%atm%lp_count_t(jb,jt)     = ext_data(jg)%atm%lp_count_t(jb,jt_in)
             ext_data(jg)%atm%idx_lst_lp_t(:,jb,jt) = ext_data(jg)%atm%idx_lst_lp_t(:,jb,jt_in)

!CDIR NODEP
             DO ic = 1, ext_data(jg)%atm%lp_count_t(jb,jt)
               jc = ext_data(jg)%atm%idx_lst_lp_t(ic,jb,jt)
               ext_data(jg)%atm%rootdp_t(jc,jb,jt)   = ext_data(jg)%atm%rootdp_t(jc,jb,jt_in)
               ext_data(jg)%atm%plcov_t(jc,jb,jt)    = ext_data(jg)%atm%plcov_t(jc,jb,jt_in)
               ext_data(jg)%atm%tai_t(jc,jb,jt)      = ext_data(jg)%atm%tai_t(jc,jb,jt_in)
               ext_data(jg)%atm%sai_t(jc,jb,jt)      = ext_data(jg)%atm%sai_t(jc,jb,jt_in)
               ext_data(jg)%atm%eai_t(jc,jb,jt)      = ext_data(jg)%atm%eai_t(jc,jb,jt_in)
               ext_data(jg)%atm%rsmin2d_t(jc,jb,jt)  = ext_data(jg)%atm%rsmin2d_t(jc,jb,jt_in)
               ext_data(jg)%atm%soiltyp_t(jc,jb,jt)  = ext_data(jg)%atm%soiltyp_t(jc,jb,jt_in)
               ext_data(jg)%atm%lc_class_t(jc,jb,jt) = ext_data(jg)%atm%lc_class_t(jc,jb,jt_in)
               ext_data(jg)%atm%lc_frac_t(jc,jb,jt)  = ext_data(jg)%atm%lc_frac_t(jc,jb,jt_in)
             ENDDO

           ENDDO
         ENDIF



         ! Initialize frac_t with lc_frac_t on all static grid points
         ! Recall: frac_t differs from lc_frac_t in the presence of time-dependent sub-lists 
         !         (snow tiles or sea ice)
         ! In this case, frac_t refers to the time-dependent sub-tiles.
         ! ** Aggregation operations always must use frac_t **
         DO jt = 1, ntiles_lnd
           DO jc = i_startidx, i_endidx
             ext_data(jg)%atm%frac_t(jc,jb,jt)  = ext_data(jg)%atm%lc_frac_t(jc,jb,jt)
           ENDDO
         ENDDO
         DO jc = i_startidx, i_endidx
           ext_data(jg)%atm%frac_t(jc,jb,isub_water)  = ext_data(jg)%atm%lc_frac_t(jc,jb,isub_water)
         ENDDO
         ! frac_t(jc,jb,isub_seaice) is set in init_sea_lists

       END DO !jb
!$OMP END DO NOWAIT
!$OMP END PARALLEL

         ! Some useful diagnostics
       npoints = SUM(ext_data(jg)%atm%lp_count(i_startblk:i_endblk))
       npoints = global_sum_array(npoints)
       WRITE(message_text,'(a,i3,a,i10)') 'Number of land points in domain',jg,':', npoints
       CALL message('', TRIM(message_text))
       npoints_sea = SUM(ext_data(jg)%atm%sp_count(i_startblk:i_endblk))
       npoints_sea = global_sum_array(npoints_sea)
       WRITE(message_text,'(a,i3,a,i10)') 'Number of sea points in domain',jg,':', npoints_sea
       CALL message('', TRIM(message_text))
       npoints_lake = SUM(ext_data(jg)%atm%fp_count(i_startblk:i_endblk))
       npoints_lake = global_sum_array(npoints_lake)
       WRITE(message_text,'(a,i3,a,i10)') 'Number of lake points in domain',jg,':', npoints_lake
       CALL message('', TRIM(message_text))


       DO i_lu = 1, ntiles_lnd
         npoints = SUM(ext_data(jg)%atm%gp_count_t(i_startblk:i_endblk,i_lu))
         npoints = global_sum_array(npoints)
         WRITE(message_text,'(a,i2,a,i10)') 'Number of points in tile',i_lu,':',npoints
         CALL message('', TRIM(message_text))
       ENDDO

    END DO  !jg

  

  END SUBROUTINE init_index_lists



  !-------------------------------------------------------------------------
  !>
  !! Get ndvi_mrat for starting time of model integration
  !!
  !! Get ndvi_mrat for starting time of model integration.
  !! Linear interpolation in time. 
  !!
  !! @par Revision History
  !! Initial revision by Juergen Helmert, DWD (2012-04-17)
  !!
  SUBROUTINE interpol_ndvi_time(p_patch, ext_data, datetime)

    TYPE(t_patch),         INTENT(IN)    :: p_patch(:)
    TYPE(t_external_data), INTENT(INOUT) :: ext_data(:)
    TYPE(t_datetime),      INTENT(IN)    :: datetime


    INTEGER :: jg, mo1, mo2
    INTEGER :: jc, jb               !< loop index
    INTEGER :: i_startblk, i_endblk, i_nchdom
    INTEGER :: rl_start, rl_end
    INTEGER :: i_startidx, i_endidx
    REAL(wp):: zw1, zw2

    CHARACTER(len=max_char_length), PARAMETER :: &
      routine = 'mo_ext_data: interpol_ndvi_time'

    !---------------------------------------------------------------

    ! Find the 2 nearest months m1, m2 and the weights pw1, pw2 
    ! to the actual date and time
    CALL month2hour( datetime, mo1, mo2, zw2 )

    zw1 = 1._wp - zw2

    ! Get interpolated ndvi_mrat
    DO jg = 1, n_dom

      i_nchdom  = MAX(1,p_patch(jg)%n_childdom)

      ! exclude the boundary interpolation zone of nested domains
      rl_start = grf_bdywidth_c+1
      rl_end   = min_rlcell_int

      i_startblk = p_patch(jg)%cells%start_blk(rl_start,1)
      i_endblk   = p_patch(jg)%cells%end_blk(rl_end,i_nchdom)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx)
      DO jb=i_startblk, i_endblk

        CALL get_indices_c(p_patch(jg), jb, i_startblk, i_endblk, &
           & i_startidx, i_endidx, rl_start, rl_end)

        DO jc = i_startidx, i_endidx
          ext_data(jg)%atm%ndvi_mrat(jc,jb) = zw1*ext_data(jg)%atm_td%ndvi_mrat(jc,jb,mo1) & 
            &                               + zw2*ext_data(jg)%atm_td%ndvi_mrat(jc,jb,mo2)
        ENDDO
      ENDDO
!$OMP END DO
!$OMP END PARALLEL

    END DO ! jg

  END SUBROUTINE interpol_ndvi_time


END MODULE mo_ext_data_state

