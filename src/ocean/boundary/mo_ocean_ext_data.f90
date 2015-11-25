!>
!! Allocation/deallocation and reading of ocean external datasets
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
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_ocean_ext_data

  USE mo_kind,               ONLY: wp
  USE mo_io_units,           ONLY: filename_max
  USE mo_parallel_config,    ONLY: nproma
  USE mo_impl_constants,     ONLY: max_char_length, LAND
  USE mo_math_constants,     ONLY: dbl_eps
  USE mo_ocean_nml,          ONLY: iforc_oce, &
    & forcing_timescale, &
    & forcing_windstress_u_type, &
    & forcing_windstress_v_type, &
    & forcing_fluxes_type,       &
    & OMIP_FluxFromFile,         &
    & use_omip_windstress, use_omip_fluxes, use_omip_forcing
  USE mo_sea_ice_nml,        ONLY: i_ice_dyn
   
  USE mo_model_domain,       ONLY: t_patch
  USE mo_exception,          ONLY: message, message_text, finish
  USE mo_grid_config,        ONLY: n_dom, nroot, dynamics_grid_filename
  USE mo_mpi,                ONLY: my_process_is_stdio, p_io, p_bcast, &
    &                              p_comm_work_test, p_comm_work
  USE mo_sync,               ONLY: global_sum_array
  USE mo_parallel_config,    ONLY: p_test_run
  USE mo_linked_list,        ONLY: t_var_list
  USE mo_ext_data_types,     ONLY: t_external_data, t_external_atmos,    &
    &                              t_external_atmos_td, t_external_ocean
  USE mo_var_list,           ONLY: default_var_list_settings,   &
    &                              add_var, add_ref,            &
    &                              new_var_list,                &
    &                              delete_var_list
  USE mo_master_config,      ONLY: getModelBaseDir
  USE mo_cf_convention,      ONLY: t_cf_var
  USE mo_grib2,              ONLY: t_grib2_var, grib2_var
  USE mo_read_interface,     ONLY: openInputFile, closeFile, on_cells, &
    &                              t_stream_id, nf, read_2D, read_2D_int, &
    &                              read_3D, on_edges
  USE mo_util_string,        ONLY: t_keyword_list,  &
    &                              associate_keyword, with_keywords
  USE mo_datetime,           ONLY: t_datetime, month2hour
  USE mo_cdi,                ONLY: DATATYPE_FLT32, DATATYPE_PACK16, TSTEP_CONSTANT, &
    &                              TSTEP_MAX, TSTEP_AVG, GRID_UNSTRUCTURED
  USE mo_cdi_constants,      ONLY: GRID_UNSTRUCTURED_CELL, GRID_UNSTRUCTURED_EDGE, &
    &                              GRID_UNSTRUCTURED_VERT,                         &
    &                              GRID_CELL, GRID_EDGE, GRID_VERTEX, ZA_SURFACE,  &
    &                              ZA_HYBRID, ZA_PRESSURE, ZA_HEIGHT_2M

  USE mo_master_config,      ONLY: isRestart


  IMPLICIT NONE

  ! required for reading external data
  INCLUDE 'netcdf.inc'

  PRIVATE

  PUBLIC :: construct_ocean_ext_data
  PUBLIC :: destruct_ocean_ext_data

  PUBLIC :: ext_data

  TYPE(t_external_data),TARGET, ALLOCATABLE :: ext_data(:)  ! n_dom

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
!<Optimize:inUse>
  SUBROUTINE construct_ocean_ext_data (p_patch, ext_data)

    TYPE(t_patch), INTENT(IN)            :: p_patch(:)
    TYPE(t_external_data), INTENT(INOUT) :: ext_data(:)

    INTEGER :: jg
    CHARACTER(len=MAX_CHAR_LENGTH) :: listname

    CHARACTER(len=max_char_length), PARAMETER :: &
      routine = 'mo_ext_data:init_ext_data'

    !-------------------------------------------------------------------------
    CALL message (TRIM(routine), 'Start')

    !-------------------------------------------------------------------------
    !  1.  inquire external files for their data structure
    !-------------------------------------------------------------------------

    !------------------------------------------------------------------
    !  2.  construct external fields for the model
    !------------------------------------------------------------------

    ! top-level procedure for building data structures for 
    ! external data.
    CALL message (TRIM(routine), 'Construction of data structure for ' // &
      &                          'external data started')


    ! write(0,*) 'create new external data list for ocean'
    ! Build external data list for constant-in-time fields for the ocean model
    DO jg = 1, n_dom
      WRITE(listname,'(a,i2.2)') 'ext_data_oce_D',jg
      CALL new_ext_data_oce_list(p_patch(jg), ext_data(jg)%oce, ext_data(jg)%oce_list, TRIM(listname))
    END DO ! jg = 1,n_dom

    ! Build external data list for time-dependent fields
    ! ### to be done ###

    CALL message (TRIM(routine), 'Construction of data structure for ' // &
      &                          'external data finished')


    !-------------------------------------------------------------------------
    !  3.  read the data into the fields
    !-------------------------------------------------------------------------

    ! Check, whether external data should be read from file
    !
    CALL read_ext_data_oce (p_patch, ext_data)

  END SUBROUTINE construct_ocean_ext_data
  !-------------------------------------------------------------------------


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
!<Optimize:inUse>
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

    LOGICAL :: use_windstress_only, use_full_file

    !--------------------------------------------------------------

    !determine size of arrays
    nblks_c = p_patch%alloc_cell_blocks
    nblks_e = p_patch%nblks_e


    ibits = 16   ! "entropy" of horizontal slice

    ! predefined array shapes
    shape2d_c = (/ nproma, nblks_c /)
    shape2d_e = (/ nproma, nblks_e /)

    ! OMIP/NCEP or other flux forcing data on cell centers: 3, 5 or 12 variables, forcing_timescale data sets
    ! for type of forcing see mo_oce_bulk
    idim_omip = 0
    use_windstress_only = (&
      & (forcing_windstress_u_type == 1 .OR. forcing_windstress_u_type == 5) .AND. &
      & (forcing_windstress_v_type == 1 .OR. forcing_windstress_v_type == 5) .AND. &
      & (forcing_fluxes_type       >  100 )                                        &
      & )
    use_full_file = ( &
      & (forcing_windstress_u_type == 1 .OR. forcing_windstress_u_type == 5) .AND. &
      & (forcing_windstress_v_type == 1 .OR. forcing_windstress_v_type == 5) .AND. &
      & (forcing_fluxes_type       == 1 .OR. forcing_fluxes_type       == 5)       &
      & )

    IF ( use_windstress_only ) idim_omip =  3    !  stress (x, y) and SST
    IF ( use_full_file )       idim_omip = 14    !  OMIP type forcing

    shape4d_c = (/ nproma, forcing_timescale, nblks_c, idim_omip /)

    !
    ! Register a field list and apply default settings
    !
    CALL new_var_list( p_ext_oce_list, TRIM(listname), patch_id=p_patch%id )
    CALL default_var_list_settings( p_ext_oce_list,            &
                                  & lrestart=.FALSE.,          &
                                 & model_type='oce'  )

    ! bathymetric height at cell center
    !
    ! bathymetry_c  p_ext_oce%bathymetry_c(nproma,nblks_c)
    cf_desc    = t_cf_var('Model bathymetry at cell center', 'm', &
      &                   'Model bathymetry', DATATYPE_FLT32)
    grib2_desc = grib2_var( 192, 140, 219, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_ext_oce_list, 'bathymetry_c', p_ext_oce%bathymetry_c,      &
      &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d_c )


    ! bathymetric height at cell edge
    !
    ! bathymetry_e  p_ext_oce%bathymetry_e(nproma,nblks_e)
    cf_desc    = t_cf_var('Model bathymetry at edge', 'm', &
      &                   'Model bathymetry', DATATYPE_FLT32)
    grib2_desc = grib2_var( 192, 140, 219, ibits, GRID_UNSTRUCTURED, GRID_EDGE)
    CALL add_var( p_ext_oce_list, 'bathymetry_e', p_ext_oce%bathymetry_e,      &
      &           GRID_UNSTRUCTURED_EDGE, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d_e )

    ! ocean land-sea-mask at surface on cell centers
    !
    ! lsm_ctr_c  p_ext_oce%lsm_ctr_c(nproma,nblks_c)
    cf_desc    = t_cf_var('Ocean model land-sea-mask at cell center', '-2/-1/1/2', &
      &                   'Ocean model land-sea-mask', DATATYPE_FLT32)
    grib2_desc = grib2_var( 192, 140, 219, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    !#slo-2011-08-08# does not compile yet?
    CALL add_var( p_ext_oce_list, 'lsm_ctr_c', p_ext_oce%lsm_ctr_c, &
      &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d_c )

    ! ocean land-sea-mask at surface on cell edge
    !
    cf_desc    = t_cf_var('Ocean model land-sea-mask at cell edge', '-2/0/2', &
      &                   'Ocean model land-sea-mask', DATATYPE_FLT32)
    grib2_desc = grib2_var( 192, 140, 219, ibits, GRID_UNSTRUCTURED, GRID_EDGE)
    CALL add_var( p_ext_oce_list, 'lsm_ctr_e', p_ext_oce%lsm_ctr_e,      &
      &           GRID_UNSTRUCTURED_EDGE, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d_e )

    ! omip forcing data on cell edge
    !
    IF (iforc_oce == OMIP_FluxFromFile) THEN
      cf_desc    = t_cf_var('Ocean model OMIP forcing data at cell edge', 'Pa, K', &
        &                   'OMIP forcing data', DATATYPE_FLT32)
      grib2_desc = grib2_var( 192, 140, 219, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_ext_oce_list, 'flux_forc_mon_c', p_ext_oce%flux_forc_mon_c,  &
        &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape4d_c )
    END IF

  END SUBROUTINE new_ext_data_oce_list
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !>
  !! Destruct external data data structure and lists
  !!
  !! Destruct external data data structure and lists
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2011-05-04)
  !!
!<Optimize:inUse>
  SUBROUTINE destruct_ocean_ext_data

    INTEGER :: jg, errstat
    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
      routine = 'mo_ext_data:destruct_ocean_ext_data'
    !-------------------------------------------------------------------------

    CALL message (TRIM(routine), 'Destruction of data structure for ' // &
      &                          'external data started')

    DO jg = 1,n_dom
      ! Delete list of constant in time oceanic elements
      CALL delete_var_list( ext_data(jg)%oce_list )
    ENDDO

    CALL message (TRIM(routine), 'Destruction of data structure for ' // &
      &                          'external data finished')

  END SUBROUTINE destruct_ocean_ext_data
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
!<Optimize:inUse>
  SUBROUTINE read_ext_data_oce (p_patch, ext_data)

    TYPE(t_patch), TARGET, INTENT(IN)            :: p_patch(:)
    TYPE(t_external_data), INTENT(INOUT) :: ext_data(:)

    CHARACTER(len=max_char_length), PARAMETER :: &
      routine = 'mo_ext_data:read_ext_data_oce'

    CHARACTER(filename_max) :: grid_file   !< file name for reading in
    CHARACTER(filename_max) :: omip_file   !< file name for reading in

    LOGICAL :: l_exist
    INTEGER :: jg, i_lev, no_cells, no_verts, no_tst
    INTEGER :: ncid, dimid
    TYPE(t_stream_id) :: stream_id
    INTEGER :: mpi_comm

    !REAL(wp):: z_flux(nproma, 12,p_patch(1)%nblks_c)
    REAL(wp):: z_flux(nproma,forcing_timescale,p_patch(1)%alloc_cell_blocks)
    TYPE (t_keyword_list), POINTER :: keywords => NULL()

    CALL message (TRIM(routine), 'start')

    !-------------------------------------------------------------------------
    !  READ OCEAN BATHYMETRY
    !-------------------------------------------------------------------------

    jg = 1

    i_lev       = p_patch(jg)%level
    z_flux(:,:,:) = 0.0_wp

    CALL associate_keyword("<path>", TRIM(getModelBaseDir()), keywords)
    grid_file = TRIM(with_keywords(keywords, TRIM(dynamics_grid_filename(jg))))

    IF(my_process_is_stdio()) THEN
      !
      ! bathymetry and lsm are read from the general ICON grid file
      ! bathymetry and lsm integrated into grid file by grid generator
      !WRITE (bathy_file,'(a,i0,a,i2.2,a)') 'iconR',nroot,'B',i_lev, '-grid.nc'
      !write(*,*) 'bathy_file = ',TRIM(bathy_file)
      !write(*,*) 'dynamics_grid_filename = ', TRIM(dynamics_grid_filename(1))

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
      CALL nf(nf_inq_dimlen(ncid, dimid, no_cells), routine)

      CALL nf(nf_inq_dimid(ncid, 'vertex', dimid), routine)
      CALL nf(nf_inq_dimlen(ncid, dimid, no_verts), routine)

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
      
      CALL nf(nf_close(ncid), routine)
    ENDIF

    stream_id = openInputFile(grid_file, p_patch(jg))

    !-------------------------------------------------------
    !
    ! Read bathymetry for triangle centers and edges
    !
    !-------------------------------------------------------
    ! These arrays are not included in standard icon-grid, but they are
    ! created by "create_ocean_grid"
    ! first initialise everything as land
    ! we need to do this since dummy entities may exist in the arryas but not in the grid data
    ext_data(jg)%oce%bathymetry_c(:,:) = 99999999.0_wp
    ext_data(jg)%oce%bathymetry_e(:,:) = 99999999.0_wp
    ext_data(jg)%oce%lsm_ctr_c(:,:)    = LAND
    ext_data(jg)%oce%lsm_ctr_e(:,:)    = LAND
     
    CALL read_2D(stream_id, on_cells, 'cell_elevation', &
      &          ext_data(jg)%oce%bathymetry_c)

    CALL read_2D(stream_id, on_edges, 'edge_elevation', &
      &          ext_data(jg)%oce%bathymetry_e)

    ! get land-sea-mask on cells, integer marks are:
    ! inner sea (-2), boundary sea (-1, cells and vertices), boundary (0, edges),
    ! boundary land (1, cells and vertices), inner land (2)
    CALL read_2D_int(stream_id, on_cells, 'cell_sea_land_mask', &
      &              ext_data(jg)%oce%lsm_ctr_c)
    CALL read_2D_int(stream_id, on_edges, 'edge_sea_land_mask', &
      &              ext_data(jg)%oce%lsm_ctr_e)

    !
    ! close file
    !
    CALL closeFile(stream_id)

    !ENDDO ! jg

    CALL message( TRIM(routine),'Ocean bathymetry for external data read' )

    !-------------------------------------------------------------------------

    !  READ OMIP FORCING

    !-------------------------------------------------------------------------

    IF ( use_omip_forcing .AND. iforc_oce == OMIP_FluxFromFile) THEN

    !DO jg = 1,n_dom
      jg = 1

      i_lev       = p_patch(jg)%level

!       WRITE (omip_file,'(a,i0,a,i2.2,a)') 'iconR',nroot,'B',i_lev, '-flux.nc'
      omip_file='ocean-flux.nc'

      !omip_file=TRIM('/pool/data/ICON/external/iconR2B04-flux.nc')
      !omip_file='/scratch/local1/m212053/ICON/trunk/icon-dev/grids/omip4icon-R2B02-monmean.nc'
      CALL message( TRIM(routine),'Ocean OMIP forcing flux file is: '//TRIM(omip_file) )

      IF(my_process_is_stdio()) THEN
        !
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
        IF(no_tst /= forcing_timescale ) THEN
          CALL finish(TRIM(ROUTINE),&
          & 'Number of forcing timesteps is not equal forcing_timescale specified in namelist - ABORT')
        ENDIF

        CALL nf(nf_close(ncid), routine)
      ENDIF

      stream_id = openInputFile(omip_file, p_patch(jg))
      
      IF(p_test_run) THEN
        mpi_comm = p_comm_work_test
      ELSE
        mpi_comm = p_comm_work
      ENDIF
      no_tst = forcing_timescale
      ! write(0,*) "no_tst=", no_tst
      ! CALL p_bcast(no_tst, p_io, mpi_comm)

      !-------------------------------------------------------
      !
      ! Read complete OMIP data for triangle centers
      !
      !-------------------------------------------------------

      ! provide OMIP fluxes for wind stress forcing
      ! 1:  'stress_x': zonal wind stress       [m/s]
      ! 2:  'stress_y': meridional wind stress  [m/s]
      ! 3:  'SST"     : sea surface temperature [K]
!     IF ( use_omip_windstress ) THEN
      ! zonal wind stress
      CALL read_3D(stream_id, on_cells, 'stress_x', z_flux(:,:,:))
      ext_data(jg)%oce%flux_forc_mon_c(:,:,:,1) = z_flux(:,:,:)

      ! meridional wind stress
      CALL read_3D(stream_id, on_cells, 'stress_y', z_flux)
      ext_data(jg)%oce%flux_forc_mon_c(:,:,:,2) = z_flux(:,:,:)

      ! SST
      CALL read_3D(stream_id, on_cells, 'SST', z_flux)
      ext_data(jg)%oce%flux_forc_mon_c(:,:,:,3) = z_flux(:,:,:)
!     ENDIF
      IF ( use_omip_fluxes ) THEN

      ! Read complete OMIP data sets for focing ocean model
      !  - names are used in type t_atmos_for_ocean in mo_se_ice_types
      !  4:  tafo(:,:),   &  ! 2 m air temperature                              [C]
      !  5:  ftdew(:,:),  &  ! 2 m dew-point temperature                        [K]
      !  6:  fu10(:,:) ,  &  ! 10 m wind speed                                  [m/s]
      !  7:  fclou(:,:),  &  ! Fractional cloud cover
      !  8:  pao(:,:),    &  ! Surface atmospheric pressure                     [hPa]
      !  9:  fswr(:,:),   &  ! Incoming surface solar radiation                 [W/m]
      ! 10:  precip(:,:), &  ! precipitation rate                               [m/s]
      ! 11:  evap  (:,:), &  ! evaporation   rate                               [m/s]
      ! 12:  runoff(:,:)     ! river runoff  rate                               [m/s]
      ! 13:  u_wind_10m   &  ! zonal wind speed                                 [m/s]
      ! 14:  v_wind_10m   &  ! meridional wind speed                            [m/s]

        ! 2m-temperature
        CALL read_3D(stream_id, on_cells, 'temp_2m', z_flux)
        ext_data(jg)%oce%flux_forc_mon_c(:,:,:,4) = z_flux(:,:,:)
     
        ! 2m dewpoint temperature
        CALL read_3D(stream_id, on_cells, 'dpt_temp_2m', z_flux)
        ext_data(jg)%oce%flux_forc_mon_c(:,:,:,5) = z_flux(:,:,:)
     
        ! Scalar wind
        CALL read_3D(stream_id, on_cells, 'scalar_wind', z_flux)
        ext_data(jg)%oce%flux_forc_mon_c(:,:,:,6) = z_flux(:,:,:)
     
        ! cloud cover
        CALL read_3D(stream_id, on_cells, 'cloud', z_flux)
        ext_data(jg)%oce%flux_forc_mon_c(:,:,:,7) = z_flux(:,:,:)
     
        ! sea level pressure
        CALL read_3D(stream_id, on_cells, 'pressure', z_flux)
        ext_data(jg)%oce%flux_forc_mon_c(:,:,:,8) = z_flux(:,:,:)
     
        ! total solar radiation
        CALL read_3D(stream_id, on_cells, 'tot_solar', z_flux)
        ext_data(jg)%oce%flux_forc_mon_c(:,:,:,9) = z_flux(:,:,:)
     
        ! precipitation
        CALL read_3D(stream_id, on_cells, 'precip', z_flux)
        ext_data(jg)%oce%flux_forc_mon_c(:,:,:,10) = z_flux(:,:,:)
     
        ! evaporation
        CALL read_3D(stream_id, on_cells, 'evap', z_flux)
        ext_data(jg)%oce%flux_forc_mon_c(:,:,:,11) = z_flux(:,:,:)
     
        ! runoff
        CALL read_3D(stream_id, on_cells, 'runoff', z_flux)
        ext_data(jg)%oce%flux_forc_mon_c(:,:,:,12) = z_flux(:,:,:)

!        IF (i_ice_dyn==1) THEN
!          ! zonal wind speed
!          CALL read_3D(stream_id, on_cells, 'u_wind_10m', z_flux)
!          ext_data(jg)%oce%flux_forc_mon_c(:,:,:,13) = z_flux(:,:,:)
!
!          ! meridional wind speed
!          CALL read_3D(stream_id, on_cells, 'v_wind_10m', z_flux)
!          ext_data(jg)%oce%flux_forc_mon_c(:,:,:,14) = z_flux(:,:,:)
!        ENDIF


      END IF

    ! TODO not needed at the moment, disabled through the restructuring of the forcing
    ! ! provide heat and freshwater flux for focing ocean model
    ! IF (iforc_type == 3) THEN

    !   ! net surface heat flux
    !   CALL read_3D(stream_id, on_cells, 'net_hflx', z_flux)
    !   ext_data(jg)%oce%flux_forc_mon_c(:,:,:,4) = z_flux(:,:,:)
    !
    !   ! surface freshwater flux
    !   CALL read_3D(stream_id, on_cells, 'net_fflx', z_flux)
    !   ext_data(jg)%oce%flux_forc_mon_c(:,:,:,5) = z_flux(:,:,:)

    ! END IF

    ! ! provide 4 parts of heat and 2 parts of freshwater flux for focing ocean model
    ! IF (iforc_type == 4) THEN

    !   ! surface short wave heat flux
    !   CALL read_3D(stream_id, on_cells, 'swflxsfc_avg', z_flux)
    !   ext_data(jg)%oce%flux_forc_mon_c(:,:,:,4) = z_flux(:,:,:)
    !
    !   ! surface long wave heat flux
    !   CALL read_3D(stream_id, on_cells, 'lwflxsfc_avg', z_flux)
    !   ext_data(jg)%oce%flux_forc_mon_c(:,:,:,5) = z_flux(:,:,:)
    !
    !   ! surface sensible heat flux
    !   CALL read_3D(stream_id, on_cells, 'shflx_avg', z_flux)
    !   ext_data(jg)%oce%flux_forc_mon_c(:,:,:,6) = z_flux(:,:,:)
    !
    !   ! surface latent heat flux
    !   CALL read_3D(stream_id, on_cells, 'lhflx_avg', z_flux)
    !   ext_data(jg)%oce%flux_forc_mon_c(:,:,:,7) = z_flux(:,:,:)
    !
    !   ! total precipiation
    !   CALL read_3D(stream_id, on_cells, 'precip', z_flux)
    !   ext_data(jg)%oce%flux_forc_mon_c(:,:,:,8) = z_flux(:,:,:)
    !
    !   ! evaporation
    !   CALL read_3D(stream_id, on_cells, 'evap', z_flux)
    !   ext_data(jg)%oce%flux_forc_mon_c(:,:,:,9) = z_flux(:,:,:)

    ! END IF

      !
      ! close file
      !
      CALL closeFile(stream_id)

    !ENDDO

      CALL message( TRIM(routine),'Ocean OMIP fluxes for external data read' )

    END IF ! iforc_oce=OMIP_FluxFromFile and iforc_type.ne.5

  END SUBROUTINE read_ext_data_oce
  !-------------------------------------------------------------------------


END MODULE mo_ocean_ext_data
