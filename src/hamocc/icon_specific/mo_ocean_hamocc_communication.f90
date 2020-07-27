!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_ocean_hamocc_communication
#ifndef NOMPI
#ifdef HAVE_YAXT
#define USE_OCEAN_HAMOCC_COMMUNICATION
#endif
#endif

  USE iso_c_binding,               ONLY: c_ptr

#ifndef USE_OCEAN_HAMOCC_COMMUNICATION
  USE mo_exception,                ONLY: finish
#endif

!  USE mo_kind,                     ONLY: wp
  USE mo_master_control,           ONLY: ocean_process, hamocc_process
  USE mo_parallel_config,          ONLY: nproma
  USE mo_model_domain,             ONLY: t_patch
  USE mo_grid_config,              ONLY: n_dom, n_dom_start
  USE mo_run_config,               ONLY: num_lev


#ifdef USE_OCEAN_HAMOCC_COMMUNICATION
  USE mpi
  USE mo_mpi,                      ONLY: get_mpi_work_intercomm, p_pe_work, &
                                         get_my_global_mpi_communicator
  USE yaxt,                        ONLY: xt_initialize, xt_initialized, &
                                         xt_redist, xt_redist_collection_new, &
                                         xt_redist_p2p_ext_new, xt_redist_p2p_new, &
                                         xt_offset_ext, xt_redist_p2p_off_new, &
                                         xt_redist_delete, xt_redist_s_exchange, &
                                         xt_int_kind, xt_xmap, xt_xmap_delete, &
                                         xt_xmap_dist_dir_intercomm_new, &
                                         xt_idxvec_new, xt_idxempty_new, &
                                         xt_idxlist_delete, xt_idxlist
  USE mo_communication,            ONLY: blk_no, idx_no
#endif

  IMPLICIT NONE

  PRIVATE

#ifdef USE_OCEAN_HAMOCC_COMMUNICATION
  INTEGER :: ocean_hamocc_intercomm
  TYPE(xt_redist) :: exchange_redist_ocean_2_hamocc
  TYPE(xt_redist) :: exchange_redist_hamocc_2_ocean
#endif

  PUBLIC :: setup_ocean_2_hamocc_communication, &
            setup_hamocc_2_ocean_communication, &
            exchange_data_ocean_2_hamocc, &
            exchange_data_hamocc_2_ocean, &
            free_ocean_hamocc_communication

CONTAINS

#ifdef USE_OCEAN_HAMOCC_COMMUNICATION
  SUBROUTINE generate_redists(is_ocean, p_patch, no_of_levels, &
                              exchange_redist_ocean_2_hamocc, &
                              exchange_redist_hamocc_2_ocean)

    LOGICAL, INTENT(IN)             :: is_ocean
    TYPE(t_patch), INTENT(IN)       :: p_patch
    INTEGER, INTENT(IN)             :: no_of_levels

    TYPE(xt_redist), INTENT(OUT) :: exchange_redist_ocean_2_hamocc
    TYPE(xt_redist), INTENT(OUT) :: exchange_redist_hamocc_2_ocean

    INTEGER :: i, j, dt_wp_nbndsw

    INTEGER(xt_int_kind), ALLOCATABLE :: global_index_cells_3d_halfLevels(:,:), global_index_edges_3d(:,:)
    INTEGER(xt_int_kind) :: n_patch_cells_g, n_patch_edges_g

    INTEGER :: num_unmasked_cells, num_unmasked_edges
    INTEGER :: idx, blk
    INTEGER, ALLOCATABLE :: offsets_cells(:), offsets_edges(:)
    TYPE(xt_offset_ext), ALLOCATABLE :: extents_cells_3d(:), extents_cells_3d_halfLevels(:), &
      extents_edges_3d(:) 

    TYPE(xt_idxlist) :: idxlist_cells_2d, idxlist_cells_3d, idxlist_cells_3d_halfLevels, &
                        idxlist_empty, idxlist_edges_3d

    TYPE(xt_xmap) :: ocean_2_hamocc_xmap_cells_2d, hamocc_2_ocean_xmap_cells_2d
    TYPE(xt_xmap) :: ocean_2_hamocc_xmap_cells_3d, hamocc_2_ocean_xmap_cells_3d
    TYPE(xt_xmap) :: ocean_2_hamocc_xmap_edges_3d
    TYPE(xt_xmap) :: ocean_2_hamocc_xmap_cells_3d_halfLevels

    TYPE(xt_redist) :: ocean_2_hamocc_redist_cells_2d
    TYPE(xt_redist) :: ocean_2_hamocc_redist_cells_3d
    TYPE(xt_redist) :: ocean_2_hamocc_redist_cells_3d_halfLevels
    TYPE(xt_redist) :: ocean_2_hamocc_redist_edges_3d
    TYPE(xt_redist) :: hamocc_2_ocean_redist_cells_2d
    TYPE(xt_redist) :: hamocc_2_ocean_redist_cells_3d
 
    INTEGER :: ierr

    num_unmasked_cells = &
      COUNT(p_patch%cells%decomp_info%owner_local == p_pe_work)
    n_patch_cells_g = INT(p_patch%n_patch_cells_g, xt_int_kind)
    num_unmasked_edges = &
      COUNT(p_patch%edges%decomp_info%owner_local == p_pe_work)
    n_patch_edges_g = INT(p_patch%n_patch_edges_g, xt_int_kind)

    ALLOCATE(global_index_cells_3d_halfLevels(no_of_levels+1, num_unmasked_cells))
    ALLOCATE(global_index_edges_3d(no_of_levels, num_unmasked_edges))

    
    j = 1
    DO i = 1, p_patch%n_patch_cells
      IF (p_patch%cells%decomp_info%owner_local(i) == p_pe_work) THEN
        global_index_cells_3d_halfLevels(1,j) = &
          INT(p_patch%cells%decomp_info%glb_index(i), xt_int_kind)
        j = j + 1
      END IF
    END DO
    DO i = 2, no_of_levels+1
      global_index_cells_3d_halfLevels(i,:) = &
        global_index_cells_3d_halfLevels(1,:) + n_patch_cells_g * INT(i - 1, xt_int_kind)
    END DO
    
    j = 1
    DO i = 1, p_patch%n_patch_edges
      IF (p_patch%edges%decomp_info%owner_local(i) == p_pe_work) THEN
        global_index_edges_3d(1,j) = &
          INT(p_patch%edges%decomp_info%glb_index(i), xt_int_kind)
        j = j + 1
      END IF
    END DO
    DO i = 2, no_of_levels
      global_index_edges_3d(i,:) = &
        global_index_edges_3d(1,:) + n_patch_cells_g * INT(i - 1, xt_int_kind)
    END DO

    idxlist_empty               = xt_idxempty_new()
    idxlist_cells_2d            = xt_idxvec_new(global_index_cells_3d_halfLevels(1,:))
    idxlist_cells_3d            = xt_idxvec_new(global_index_cells_3d_halfLevels(1:no_of_levels,:))
    idxlist_cells_3d_halfLevels = xt_idxvec_new(global_index_cells_3d_halfLevels(:,:))
    idxlist_edges_3d            = xt_idxvec_new(global_index_edges_3d(:,:))   

    DEALLOCATE(global_index_cells_3d_halfLevels, global_index_edges_3d)
    
    ALLOCATE(offsets_cells(num_unmasked_cells), &
             offsets_edges(num_unmasked_edges), &
             extents_cells_3d(num_unmasked_cells), &
             extents_cells_3d_halfLevels(num_unmasked_cells), &
             extents_edges_3d(num_unmasked_edges))
    j = 1
    DO i = 1, p_patch%n_patch_cells
      IF (p_patch%cells%decomp_info%owner_local(i) == p_pe_work) THEN
        idx = idx_no(i) - 1
        blk = blk_no(i) - 1
        offsets_cells(j) = i - 1
        extents_cells_3d(j)%start = idx + blk * nproma * no_of_levels
        extents_cells_3d(j)%size = no_of_levels
        extents_cells_3d(j)%stride = nproma
        extents_cells_3d_halfLevels(j)%start = idx + blk * nproma * (no_of_levels + 1)
        extents_cells_3d_halfLevels(j)%size = no_of_levels + 1
        extents_cells_3d_halfLevels(j)%stride = nproma
        j = j + 1;
      END IF
    END DO
    j = 1
    DO i = 1, p_patch%n_patch_edges
      IF (p_patch%edges%decomp_info%owner_local(i) == p_pe_work) THEN
        idx = idx_no(i) - 1
        blk = blk_no(i) - 1
        offsets_edges(j) = i - 1
        extents_edges_3d(j)%start = idx + blk * nproma * no_of_levels
        extents_edges_3d(j)%size = no_of_levels
        extents_edges_3d(j)%stride = nproma
        j = j + 1;
      END IF
    END DO
    
    
    IF (is_ocean) THEN
       ocean_2_hamocc_xmap_cells_2d = &
        xt_xmap_dist_dir_intercomm_new( &
          idxlist_cells_2d, idxlist_empty, ocean_hamocc_intercomm, p_patch%comm)
      ocean_2_hamocc_xmap_cells_3d = &
        xt_xmap_dist_dir_intercomm_new( &
          idxlist_cells_3d, idxlist_empty, ocean_hamocc_intercomm, p_patch%comm)
      ocean_2_hamocc_xmap_cells_3d_halfLevels = &
        xt_xmap_dist_dir_intercomm_new( &
          idxlist_cells_3d_halfLevels, idxlist_empty, ocean_hamocc_intercomm, p_patch%comm)
      ocean_2_hamocc_xmap_edges_3d = &
        xt_xmap_dist_dir_intercomm_new( &
          idxlist_edges_3d, idxlist_empty, ocean_hamocc_intercomm, p_patch%comm)
      hamocc_2_ocean_xmap_cells_2d = &
        xt_xmap_dist_dir_intercomm_new( &
          idxlist_empty, idxlist_cells_2d, ocean_hamocc_intercomm, p_patch%comm)
      hamocc_2_ocean_xmap_cells_3d = &
        xt_xmap_dist_dir_intercomm_new( &
          idxlist_empty, idxlist_cells_3d, ocean_hamocc_intercomm, p_patch%comm)
   
    ELSE
    
      ocean_2_hamocc_xmap_cells_2d = &
        xt_xmap_dist_dir_intercomm_new( &
          idxlist_empty, idxlist_cells_2d, ocean_hamocc_intercomm, p_patch%comm)
      ocean_2_hamocc_xmap_cells_3d = &
        xt_xmap_dist_dir_intercomm_new( &
          idxlist_empty, idxlist_cells_3d, ocean_hamocc_intercomm, p_patch%comm)
      ocean_2_hamocc_xmap_cells_3d_halfLevels = &
        xt_xmap_dist_dir_intercomm_new( &
          idxlist_empty, idxlist_cells_3d_halfLevels, ocean_hamocc_intercomm, p_patch%comm)
      ocean_2_hamocc_xmap_edges_3d = &
        xt_xmap_dist_dir_intercomm_new( &
          idxlist_empty, idxlist_edges_3d, ocean_hamocc_intercomm, p_patch%comm)
      hamocc_2_ocean_xmap_cells_2d = &
        xt_xmap_dist_dir_intercomm_new( &
          idxlist_cells_2d, idxlist_empty, ocean_hamocc_intercomm, p_patch%comm)
      hamocc_2_ocean_xmap_cells_3d = &
        xt_xmap_dist_dir_intercomm_new( &
          idxlist_cells_3d, idxlist_empty, ocean_hamocc_intercomm, p_patch%comm)
    END IF

    ocean_2_hamocc_redist_cells_2d = &
      xt_redist_p2p_off_new( &
        ocean_2_hamocc_xmap_cells_2d, offsets_cells, offsets_cells, MPI_DOUBLE_PRECISION)
    hamocc_2_ocean_redist_cells_2d = &
      xt_redist_p2p_off_new( &
        hamocc_2_ocean_xmap_cells_2d, offsets_cells, offsets_cells, MPI_DOUBLE_PRECISION)
    ocean_2_hamocc_redist_cells_3d = &
      xt_redist_p2p_ext_new( &
        ocean_2_hamocc_xmap_cells_3d, extents_cells_3d, extents_cells_3d, MPI_DOUBLE_PRECISION)
    hamocc_2_ocean_redist_cells_3d = &
      xt_redist_p2p_ext_new( &
        hamocc_2_ocean_xmap_cells_3d, extents_cells_3d, extents_cells_3d, MPI_DOUBLE_PRECISION)
    ocean_2_hamocc_redist_cells_3d_halfLevels = &
      xt_redist_p2p_ext_new( &
        ocean_2_hamocc_xmap_cells_3d_halfLevels, extents_cells_3d_halfLevels, extents_cells_3d_halfLevels, &
        MPI_DOUBLE_PRECISION)

    exchange_redist_ocean_2_hamocc = &
      xt_redist_collection_new((/ocean_2_hamocc_redist_cells_2d, & ! top_dilution_coeff
                                 ocean_2_hamocc_redist_cells_2d, & ! h_old
                                 ocean_2_hamocc_redist_cells_2d, & ! h_new
                                 ocean_2_hamocc_redist_cells_2d, & ! ice_concentration_sum
                                 ocean_2_hamocc_redist_cells_3d, & ! temperature
                                 ocean_2_hamocc_redist_cells_3d, & ! salinity
                                 ocean_2_hamocc_redist_cells_3d_halfLevels, & !ver_diffusion_coeff
!                                  ocean_2_hamocc_redist_edges_3d, & ! hor_diffusion_coeff
                                 ocean_2_hamocc_redist_cells_2d, & ! short_wave_flux
                                 ocean_2_hamocc_redist_cells_2d, & ! wind10m
                                 ocean_2_hamocc_redist_cells_2d, & ! co2 mixing ratio
                                 ocean_2_hamocc_redist_edges_3d, & ! transport: mass_flux_e
                                 ocean_2_hamocc_redist_edges_3d, & ! transport: vn
                                 ocean_2_hamocc_redist_cells_3d_halfLevels & ! transport:w
                                 /), & !
                               ocean_hamocc_intercomm)

    exchange_redist_hamocc_2_ocean = &
      xt_redist_collection_new((/hamocc_2_ocean_redist_cells_2d/), & !  co2_flux
                                ocean_hamocc_intercomm)
      ! hamocc_2_ocean_redist_cells_3d,   & !  swr_fraction, in the future
      
    ! clean up
    CALL xt_redist_delete(ocean_2_hamocc_redist_cells_2d)
    CALL xt_redist_delete(ocean_2_hamocc_redist_cells_3d)
    CALL xt_redist_delete(ocean_2_hamocc_redist_cells_3d_halfLevels)
    CALL xt_redist_delete(ocean_2_hamocc_redist_edges_3d)
    CALL xt_redist_delete(hamocc_2_ocean_redist_cells_2d)
    CALL xt_redist_delete(hamocc_2_ocean_redist_cells_3d)


    CALL xt_xmap_delete(ocean_2_hamocc_xmap_cells_2d)
    CALL xt_xmap_delete(ocean_2_hamocc_xmap_cells_3d)
    CALL xt_xmap_delete(ocean_2_hamocc_xmap_cells_3d_halfLevels)
    CALL xt_xmap_delete(ocean_2_hamocc_xmap_edges_3d)
    CALL xt_xmap_delete(hamocc_2_ocean_xmap_cells_2d)
    CALL xt_xmap_delete(hamocc_2_ocean_xmap_cells_3d)

    CALL xt_idxlist_delete(idxlist_cells_2d)
    CALL xt_idxlist_delete(idxlist_cells_3d)
    CALL xt_idxlist_delete(idxlist_cells_3d_halfLevels)
    CALL xt_idxlist_delete(idxlist_edges_3d)
    CALL xt_idxlist_delete(idxlist_empty)

  END SUBROUTINE generate_redists
#endif

#ifdef USE_OCEAN_HAMOCC_COMMUNICATION
  SUBROUTINE setup_communication(is_ocean, p_patch, no_of_levels)

    LOGICAL, INTENT(IN) :: is_ocean
    TYPE(t_patch), INTENT(IN)       :: p_patch
    INTEGER, INTENT(in) :: no_of_levels

    IF (.NOT. xt_initialized()) &
      CALL xt_initialize(get_my_global_mpi_communicator())

    CALL generate_redists(is_ocean, p_patch, no_of_levels, &
                            exchange_redist_ocean_2_hamocc, &
                            exchange_redist_hamocc_2_ocean)

  END SUBROUTINE setup_communication
#endif

  SUBROUTINE setup_ocean_2_hamocc_communication(p_patch, no_of_levels)
    TYPE(t_patch), INTENT(IN)       :: p_patch
    INTEGER, INTENT(in) :: no_of_levels

#ifdef USE_OCEAN_HAMOCC_COMMUNICATION
    ocean_hamocc_intercomm = get_mpi_work_intercomm(hamocc_process)

    CALL setup_communication(.TRUE., p_patch, no_of_levels)
#endif

  END SUBROUTINE setup_ocean_2_hamocc_communication

  SUBROUTINE setup_hamocc_2_ocean_communication(p_patch, no_of_levels)
    TYPE(t_patch), INTENT(IN)       :: p_patch
    INTEGER, INTENT(in) :: no_of_levels

#ifdef USE_OCEAN_HAMOCC_COMMUNICATION
    ocean_hamocc_intercomm = get_mpi_work_intercomm(ocean_process)

    CALL setup_communication(.FALSE., p_patch, no_of_levels)
#endif

  END SUBROUTINE setup_hamocc_2_ocean_communication

  SUBROUTINE free_ocean_hamocc_communication
 
#ifdef USE_OCEAN_HAMOCC_COMMUNICATION
   INTEGER :: ierr

   CALL xt_redist_delete(exchange_redist_ocean_2_hamocc)
   CALL xt_redist_delete(exchange_redist_hamocc_2_ocean)

   CALL MPI_COMM_FREE(ocean_hamocc_intercomm, ierr)
#endif

  END SUBROUTINE free_ocean_hamocc_communication

  SUBROUTINE exchange_data_ocean_2_hamocc( &
    &  top_dilution_coeff, &
    &  h_old, &
    &  h_new, &
    &  ice_concentration_sum, &
    &  temperature, &
    &  salinity, &
    &  ver_diffusion_coeff, &
    &  short_wave_flux, &
    &  wind10m, &
    &  co2_mixing_ratio, &
    &  mass_flux_e, &
    &  vn, &
    &  w)


    !INTEGER,INTENT(INOUT)  :: &
    TYPE(c_ptr), INTENT(IN) :: &
      &  top_dilution_coeff, &
      &  h_old, &
      &  h_new, &
      &  ice_concentration_sum, &
      &  temperature, &
      &  salinity, &
      &  ver_diffusion_coeff, &
      &  short_wave_flux, &
      &  wind10m, &
      &  co2_mixing_ratio, &
      &  mass_flux_e, &
      &  vn, &
      &  w


    TYPE(c_ptr) :: src_data_cptr(13), dst_data_cptr(13)

    src_data_cptr( 1) =  top_dilution_coeff
    src_data_cptr( 2) =  h_old
    src_data_cptr( 3) =  h_new
    src_data_cptr( 4) =  ice_concentration_sum
    src_data_cptr( 5) =  temperature
    src_data_cptr( 6) =  salinity
    src_data_cptr( 7) =  ver_diffusion_coeff
    src_data_cptr( 8) =  short_wave_flux
    src_data_cptr( 9) =  wind10m
    src_data_cptr(10) =  co2_mixing_ratio
    src_data_cptr(11) =  mass_flux_e
    src_data_cptr(12) =  vn
    src_data_cptr(13) =  w

    dst_data_cptr = src_data_cptr

#ifdef USE_OCEAN_HAMOCC_COMMUNICATION
    CALL xt_redist_s_exchange(exchange_redist_ocean_2_hamocc, &
                              src_data_cptr, dst_data_cptr)
#else
    CALL finish("exchange_data_ocean_2_hamocc", " Requires the YAXT library")
#endif


  END SUBROUTINE exchange_data_ocean_2_hamocc


  SUBROUTINE exchange_data_hamocc_2_ocean(co2_flux)

    TYPE(c_ptr), INTENT(IN)  :: co2_flux

    TYPE(c_ptr) :: src_data_cptr(1), dst_data_cptr(1)
   
    src_data_cptr( 1) = co2_flux
    dst_data_cptr = src_data_cptr

#ifdef USE_OCEAN_HAMOCC_COMMUNICATION
    CALL xt_redist_s_exchange(exchange_redist_hamocc_2_ocean, &
                              src_data_cptr, dst_data_cptr)
#else
    CALL finish("exchange_data_ocean_2_hamocc", " Requires the YAXT library")
#endif

  END SUBROUTINE exchange_data_hamocc_2_ocean

END MODULE mo_ocean_hamocc_communication
