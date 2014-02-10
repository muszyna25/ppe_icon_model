!>
!! Module handling the collection of grid information for synchronous
!! and asynchronous output.
!!
!! @author R. Johanni, F. Prill
!!
MODULE mo_name_list_output_gridinfo

  USE mo_cdi_constants          ! We need all
  USE mo_kind,                              ONLY: wp
  USE mo_parallel_config,                   ONLY: nproma
  USE mo_exception,                         ONLY: finish
  USE mo_model_domain,                      ONLY: t_patch
  USE mo_math_utilities,                    ONLY: t_geographical_coordinates,               &
    &                                             check_orientation
  USE mo_communication,                     ONLY: t_comm_gather_pattern, exchange_data
  USE mo_grib2,                             ONLY: t_grib2_var
  USE mo_lonlat_grid,                       ONLY: t_lon_lat_grid, compute_lonlat_specs,     &
    &                                             rotate_latlon_grid
  USE mo_intp_data_strc,                    ONLY: lonlat_grid_list
  USE mo_math_constants,                    ONLY: pi_180, pi
  USE mo_impl_constants,                    ONLY: SUCCESS, min_rlcell_int, min_rledge_int,  &
    &                                             min_rlvert
  USE mo_mpi,                               ONLY: p_comm_work_2_io,                         &
    &                                             my_process_is_mpi_test, my_process_is_io, &
    &                                             my_process_is_mpi_workroot, p_bcast
  USE mo_master_control,                    ONLY: my_process_is_ocean
  USE mo_gribout_config,                    ONLY: gribout_config, t_gribout_config
  USE mo_loopindices,                       ONLY: get_indices_c, get_indices_e, get_indices_v
  USE mo_util_cdi,                          ONLY: set_additional_GRIB2_keys

  USE mo_name_list_output_types,            ONLY: t_patch_info, t_grid_info, t_output_file, &
    &                                             REMAP_NONE, REMAP_REGULAR_LATLON,         &
    &                                             ILATLON, ICELL, IEDGE, IVERT, IRLAT,      &
    &                                             IRLON

  IMPLICIT NONE

  PRIVATE

  !-----------------------------------------------------------------
  ! include NetCDF headers (direct NetCDF library calls are required
  ! for output of grid information).
  !-----------------------------------------------------------------
  INCLUDE 'netcdf.inc'

  ! constants
  PUBLIC :: GRID_INFO_NONE, GRID_INFO_FILE, GRID_INFO_BCAST
  ! subroutines
  PUBLIC :: deallocate_all_grid_info
  PUBLIC :: collect_all_grid_info
  PUBLIC :: set_grid_info_netcdf
  PUBLIC :: set_grid_info_grb2
  PUBLIC :: copy_grid_info
  PUBLIC :: bcast_grid_info
  PUBLIC :: write_grid_info_grb2

  ! constants: modes how to collect grid information (for output)
  INTEGER, PARAMETER :: GRID_INFO_NONE  = 1
  INTEGER, PARAMETER :: GRID_INFO_FILE  = 2
  INTEGER, PARAMETER :: GRID_INFO_BCAST = 3


  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_name_list_output_gridinfo'

CONTAINS

  SUBROUTINE deallocate_grid_info(grid_info)
    TYPE(t_grid_info),  INTENT(INOUT) :: grid_info
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::deallocate_grid_info"
    ! local variables
    INTEGER :: ierrstat

    IF (ALLOCATED(grid_info%lon)) THEN
      DEALLOCATE(grid_info%lon,  grid_info%lat,  &
        &        grid_info%lonv, grid_info%latv, &
        &        STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed.')
    END IF
    IF (ALLOCATED(grid_info%log_dom_index)) THEN
      DEALLOCATE(grid_info%log_dom_index,        &
        &        STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed.')
    END IF
  END SUBROUTINE deallocate_grid_info


  SUBROUTINE deallocate_all_grid_info(patch_info)
    TYPE(t_patch_info),   INTENT(INOUT) :: patch_info
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::deallocate_grid_info"
    ! local variables
    INTEGER :: ierrstat

    CALL deallocate_grid_info(patch_info%cells%grid_info)
    CALL deallocate_grid_info(patch_info%edges%grid_info)
    CALL deallocate_grid_info(patch_info%verts%grid_info)
  END SUBROUTINE deallocate_all_grid_info


  SUBROUTINE collect_all_grid_info(p_patch, patch_info)
    TYPE(t_patch),        INTENT(IN)    :: p_patch 
    TYPE(t_patch_info),   INTENT(INOUT) :: patch_info 
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::collect_all_grid_info"
    INTEGER                             :: idom_log, ierrstat, max_cell_connectivity
    REAL(wp), ALLOCATABLE               :: lonv(:,:,:), latv(:,:,:)

    ! logical domain ID
    idom_log = patch_info%log_patch_id
    max_cell_connectivity = p_patch%cells%max_connectivity
    !-- collect domain data on working PE 0
    ! --cells
    ALLOCATE(lonv(nproma, p_patch%nblks_c, max_cell_connectivity), &
      &      latv(nproma, p_patch%nblks_c, max_cell_connectivity), &
      &      STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')
    CALL cf_1_1_grid_cells(p_patch, lonv, latv)
    CALL collect_grid_info(patch_info%nblks_glb_c,              &
      &                    p_patch%nblks_c,                     &
      &                    p_patch%cells%center,                &
      &                    lonv, latv,                          &
      &                    patch_info%cells%grid_info,          &
      &                    max_cell_connectivity,                    &
      &                    patch_info%p_pat_c)
    DEALLOCATE(lonv, latv, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed.')

    !-- edges
    ALLOCATE(lonv(nproma, p_patch%nblks_e, 4), &
      &      latv(nproma, p_patch%nblks_e, 4), &
      &      STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')
    CALL cf_1_1_grid_edges(p_patch, lonv, latv)
    CALL collect_grid_info(patch_info%nblks_glb_e,              &
      &                    p_patch%nblks_e,                     &
      &                    p_patch%edges%center,                &
      &                    lonv, latv,                          &
      &                    patch_info%edges%grid_info,          &
      &                    4,                                   &
      &                    patch_info%p_pat_e)
    DEALLOCATE(lonv, latv, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed.')

    !-- verts
    ALLOCATE(lonv(nproma, p_patch%nblks_v, 9-max_cell_connectivity), &
      &      latv(nproma, p_patch%nblks_v, 9-max_cell_connectivity), &
      &      STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')

    IF (my_process_is_ocean()) THEN
      CALL cf_1_1_grid_verts_ocean(p_patch, lonv, latv)
    ELSE
      CALL cf_1_1_grid_verts(p_patch, lonv, latv)
    ENDIF

    CALL collect_grid_info(patch_info%nblks_glb_v,              &
      &                    p_patch%nblks_v,                     &
      &                    p_patch%verts%vertex,                &
      &                    lonv, latv,                          &
      &                    patch_info%verts%grid_info,          &
      &                    9-max_cell_connectivity,                  &
      &                    patch_info%p_pat_v)
    DEALLOCATE(lonv, latv, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed.')
  END SUBROUTINE collect_all_grid_info


  !------------------------------------------------------------------------------------------------
  !> SUBROUTINE collect_grid_info
  !
  !  Prepare the output of grid information: For each physical/logical
  !  patch we must collect the geographical locations of cells, edges,
  !  and vertices first on working PE 0 - from where it will be
  !  broadcasted to the pure I/O PEs.
  !
  SUBROUTINE collect_grid_info(nblks_glb, nblks_loc, in_lonlat, lonv, latv, out_lonlat, &
    &                          dim3, p_pat)
    INTEGER,                          INTENT(IN)    :: nblks_glb,     &  ! global number of blocks
      &                                                nblks_loc         ! local  number of blocks
    TYPE(t_geographical_coordinates), INTENT(IN)    :: in_lonlat(:,:)
    REAL(wp),                         INTENT(IN)    :: lonv(:,:,:), latv(:,:,:)
    TYPE(t_grid_info),                INTENT(INOUT) :: out_lonlat
    INTEGER,                          INTENT(IN)    :: dim3
    TYPE(t_comm_gather_pattern),  POINTER           :: p_pat
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine  =  modname//"::collect_grid_info"
    INTEGER               :: ierrstat, jb, jc, idim
    REAL(wp), ALLOCATABLE :: r_tmp_lon(:,:), r_tmp_lat(:,:)

    ! skip this on test PE...
    IF (my_process_is_mpi_test()) RETURN

    ! allocate destination (on work root)
    IF ( my_process_is_mpi_workroot() ) THEN
      ALLOCATE(out_lonlat%lon (nproma*nblks_glb),      &
        &      out_lonlat%lat (nproma*nblks_glb),      &
        &      out_lonlat%lonv(dim3,nproma*nblks_glb), &
        &      out_lonlat%latv(dim3,nproma*nblks_glb), &
        &      STAT=ierrstat)
    ELSE
      ALLOCATE(out_lonlat%lon(1), out_lonlat%lat (1),            &
        &      out_lonlat%lonv(dim3,1), out_lonlat%latv(dim3,1), &
        &      STAT=ierrstat)
    END IF
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')

    ! allocate temporary fields:
    IF ( my_process_is_mpi_workroot()) THEN
      ALLOCATE(r_tmp_lon (nproma, nblks_glb), r_tmp_lat (nproma, nblks_glb), &
        &      STAT=ierrstat)
    ELSE
      ALLOCATE(r_tmp_lon (nproma, nblks_loc), r_tmp_lat (nproma, nblks_loc), &
        &      STAT=ierrstat)
    ENDIF
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')

    !-- part 1: exchange lon/lat coordinates:
    IF (.NOT. my_process_is_mpi_test()) THEN

      ! copy coordinates into sender array:
      DO jb=1,nblks_loc
        DO jc=1,nproma
          r_tmp_lon(jc,jb) = in_lonlat(jc,jb)%lon
          r_tmp_lat(jc,jb) = in_lonlat(jc,jb)%lat
        END DO
      END DO

      CALL exchange_data(r_tmp_lon(:,:), out_lonlat%lon(:), p_pat)
      CALL exchange_data(r_tmp_lat(:,:), out_lonlat%lat(:), p_pat)

      !-- part 2: exchange vertex lon/lat coordinates:
      DO idim=1,dim3

        CALL exchange_data(lonv(1:nproma,1:nblks_loc,idim), &
          &                out_lonlat%lonv(idim,:), p_pat)
        CALL exchange_data(latv(1:nproma,1:nblks_loc,idim), &
          &                out_lonlat%latv(idim,:), p_pat)

      END DO ! idim

    ENDIF !  (.NOT. my_process_is_mpi_test())

    ! clean up
    DEALLOCATE(r_tmp_lon, r_tmp_lat, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed.')

  END SUBROUTINE collect_grid_info


  !------------------------------------------------------------------------------------------------
  !> Reformat vertex lon/lat coordinates "clonv", "clatv" into a CF1.1 compatible form.
  !
  !  based on SUBROUTINE mo_gridrefinement::write_patch(p)
  !
  SUBROUTINE cf_1_1_grid_cells(p_patch, lonv, latv)
    TYPE(t_patch),      INTENT(IN)    :: p_patch
    REAL(wp),           INTENT(INOUT) :: lonv(:,:,:), latv(:,:,:)
    ! local variables
    INTEGER :: jc, jb, j, iidx, iblk,                  &
      &        rl_start, rl_end, i_startblk, i_endblk, &
      &        i_startidx, i_endidx, i_nchdom
    INTEGER :: max_cell_connectivity

    rl_start   = 1
    rl_end     = min_rlcell_int
    i_startblk = p_patch%cells%start_blk(rl_start,1)
    i_nchdom   = MAX(1,p_patch%n_childdom)
    i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)
    max_cell_connectivity = p_patch%cells%max_connectivity

    lonv(:,:,:) = 0.0_wp
    latv(:,:,:) = 0.0_wp
    DO j = 1, max_cell_connectivity
      DO jb = i_startblk, i_endblk
        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
          i_startidx, i_endidx, rl_start, rl_end)
        DO jc = i_startidx, i_endidx
          iidx = p_patch%cells%vertex_idx(jc,jb,j)
          iblk = p_patch%cells%vertex_blk(jc,jb,j)
          lonv(jc,jb,j) = p_patch%verts%vertex(iidx,iblk)%lon
          latv(jc,jb,j) = p_patch%verts%vertex(iidx,iblk)%lat
        END DO
      END DO
    END DO
    WHERE (ABS(lonv(:,:,:)) < EPSILON(0.0_wp))
      lonv(:,:,:) = 0.0_wp
    END WHERE
    WHERE (ABS(latv(:,:,:)) < EPSILON(0.0_wp))
      latv(:,:,:) = 0.0_wp
    END WHERE
    DO j = 1, max_cell_connectivity
      DO jb = i_startblk, i_endblk
        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
          i_startidx, i_endidx, rl_start, rl_end)
        DO jc = i_startidx, i_endidx
          IF (ABS(latv(jc,jb,j)) > 0.5_wp*pi-EPSILON(0.0_wp)) THEN
            lonv(jc,jb,j) = p_patch%cells%center(jc,jb)%lon
          ENDIF
        END DO
      END DO
    END DO

  END SUBROUTINE cf_1_1_grid_cells


  !------------------------------------------------------------------------------------------------
  !> Reformat vertex lon/lat coordinates "elonv", "elatv" into a CF1.1 compatible form.
  !
  !  based on SUBROUTINE mo_gridrefinement::write_patch(p)
  !
  SUBROUTINE cf_1_1_grid_edges(p_patch, lonv, latv)
    TYPE(t_patch),      INTENT(IN)    :: p_patch
    REAL(wp),           INTENT(INOUT) :: lonv(:,:,:), latv(:,:,:)
    ! local variables
    INTEGER  :: jc, jb, j, iidx, iblk,                 &
      &        rl_start, rl_end, i_startblk, i_endblk, &
      &        i_startidx, i_endidx, i_nchdom
    REAL(wp) :: swap(4)

    rl_start   = 1
    rl_end     = min_rledge_int
    i_startblk = p_patch%edges%start_blk(rl_start,1)
    i_nchdom   = MAX(1,p_patch%n_childdom)
    i_endblk   = p_patch%edges%end_blk(rl_end,i_nchdom)

    lonv(:,:,:) = 0.0_wp
    latv(:,:,:) = 0.0_wp
    DO jb = i_startblk, i_endblk
      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
        &                i_startidx, i_endidx, rl_start, rl_end)
      DO jc = i_startidx, i_endidx
        iidx = p_patch%edges%vertex_idx(jc,jb,1)
        iblk = p_patch%edges%vertex_blk(jc,jb,1)
        lonv(jc,jb,1) = p_patch%verts%vertex(iidx,iblk)%lon
        latv(jc,jb,1) = p_patch%verts%vertex(iidx,iblk)%lat

        iidx = p_patch%edges%vertex_idx(jc,jb,2)
        iblk = p_patch%edges%vertex_blk(jc,jb,2)
        lonv(jc,jb,3) = p_patch%verts%vertex(iidx,iblk)%lon
        latv(jc,jb,3) = p_patch%verts%vertex(iidx,iblk)%lat

        IF (p_patch%edges%cell_idx(jc,jb,1) > 0) THEN
          iidx = p_patch%edges%cell_idx(jc,jb,1)
          iblk = p_patch%edges%cell_blk(jc,jb,1)
          lonv(jc,jb,4) = p_patch%cells%center(iidx,iblk)%lon
          latv(jc,jb,4) = p_patch%cells%center(iidx,iblk)%lat
        ELSE
          lonv(jc,jb,4) = 0._wp
          latv(jc,jb,4) = 0._wp
        ENDIF

        IF (p_patch%edges%cell_idx(jc,jb,2) > 0) THEN
          iidx = p_patch%edges%cell_idx(jc,jb,2)
          iblk = p_patch%edges%cell_blk(jc,jb,2)
          lonv(jc,jb,2) = p_patch%cells%center(iidx,iblk)%lon
          latv(jc,jb,2) = p_patch%cells%center(iidx,iblk)%lat
        ELSE
          lonv(jc,jb,2) = 0._wp
          latv(jc,jb,2) = 0._wp
        END IF
      END DO
    END DO

    WHERE (ABS(lonv(:,:,:)) < EPSILON(0.0_wp))
      lonv(:,:,:) = 0.0_wp
    END WHERE
    WHERE (ABS(latv(:,:,:)) < EPSILON(0.0_wp))
      latv(:,:,:) = 0.0_wp
    END WHERE

    DO j = 1, 4
      DO jb = i_startblk, i_endblk
        CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
          &                i_startidx, i_endidx, rl_start, rl_end)
        DO jc = i_startidx, i_endidx
          IF ( ABS(latv(jc,jb,j)) > 0.5_wp*pi-EPSILON(0.0_wp)) THEN
            lonv(jc,jb,j) = p_patch%edges%center(jc,jb)%lon
          ENDIF
        ENDDO
      ENDDO
    END DO
    DO jb = i_startblk, i_endblk
      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
        &                i_startidx, i_endidx, rl_start, rl_end)
      DO jc = i_startidx, i_endidx
        IF ( check_orientation(p_patch%edges%center(jc,jb)%lon, &
          &                    lonv(jc,jb,:), latv(jc,jb,:), 4) < 0 ) THEN
          swap(1:4) = lonv(jc,jb,4:1:-1)
          lonv(jc,jb,:) = swap(:)
          swap(1:4) = latv(jc,jb,4:1:-1)
          latv(jc,jb,:) = swap(:)
        END IF
      END DO
    ENDDO

  END SUBROUTINE cf_1_1_grid_edges


  !------------------------------------------------------------------------------------------------
  !> Reformat vertex lon/lat coordinates "vlonv", "vlatv" into a CF1.1 compatible form for ocean
  ! This is needed for general connectivity grids, the original routine  cf_1_1_grid_verts will not work
  !
  ! This is also sloppy for boundary vertices, as they seem not to be defined
  !
  !  based on SUBROUTINE cf_1_1_grid_verts
  !
  SUBROUTINE cf_1_1_grid_verts_ocean(patch_2D, lonv, latv)
    TYPE(t_patch),      INTENT(IN)    :: patch_2D
    REAL(wp),           INTENT(INOUT) :: lonv(:,:,:), latv(:,:,:)
    ! local variables
    INTEGER :: jc, jb, j, iidx, iblk,                  &
      &        rl_start, rl_end, i_startblk, i_endblk, &
      &        i_startidx, i_endidx, i_nchdom
    INTEGER :: last_valid_cell

    rl_start   = 2
    rl_end     = min_rlvert
    i_nchdom   = MAX(1,patch_2D%n_childdom)
    i_startblk = patch_2D%verts%start_blk(rl_start,1)
    i_endblk   = patch_2D%verts%end_blk(rl_end,i_nchdom)

    lonv(:,:,:) = 0.0_wp
    latv(:,:,:) = 0.0_wp
    DO jb = i_startblk, i_endblk
      CALL get_indices_v(patch_2D, jb, i_startblk, i_endblk, &
        &                i_startidx, i_endidx, rl_start, rl_end)
      DO jc = i_startidx, i_endidx
        last_valid_cell = 0
        DO j = 1,patch_2D%verts%max_connectivity
          IF ((patch_2D%verts%cell_idx(jc,jb,j) > 0)) THEN
            last_valid_cell = j
            iidx = patch_2D%verts%cell_idx(jc,jb,j)
            iblk = patch_2D%verts%cell_blk(jc,jb,j)
            lonv(jc,jb, j) = patch_2D%cells%center(iidx,iblk)%lon
            latv(jc,jb, j) = patch_2D%cells%center(iidx,iblk)%lat
          ELSE
            iidx = patch_2D%verts%cell_idx(jc,jb,last_valid_cell)
            iblk = patch_2D%verts%cell_blk(jc,jb,last_valid_cell)
            lonv(jc,jb, j) = patch_2D%cells%center(iidx,iblk)%lon
            latv(jc,jb, j) = patch_2D%cells%center(iidx,iblk)%lat
          ENDIF
        ENDDO
      ENDDO
    END DO
  END SUBROUTINE cf_1_1_grid_verts_ocean
  !------------------------------------------------------------------------------------------------



  !------------------------------------------------------------------------------------------------
  !> Reformat vertex lon/lat coordinates "vlonv", "vlatv" into a CF1.1 compatible form.
  !
  !  based on SUBROUTINE mo_gridrefinement::write_patch(p)
  !
  SUBROUTINE cf_1_1_grid_verts(p_patch, lonv, latv)
    TYPE(t_patch),      INTENT(IN)    :: p_patch
    REAL(wp),           INTENT(INOUT) :: lonv(:,:,:), latv(:,:,:)
    ! local variables
    INTEGER :: jc, jb, j, iidx, iblk,                  &
      &        rl_start, rl_end, i_startblk, i_endblk, &
      &        i_startidx, i_endidx, i_nchdom
    INTEGER :: max_cell_connectivity

    rl_start   = 2
    rl_end     = min_rlvert
    i_nchdom   = MAX(1,p_patch%n_childdom)
    i_startblk = p_patch%verts%start_blk(rl_start,1)
    i_endblk   = p_patch%verts%end_blk(rl_end,i_nchdom)
    max_cell_connectivity = p_patch%cells%max_connectivity

    lonv(:,:,:) = 0.0_wp
    latv(:,:,:) = 0.0_wp
    DO j = 1,(9-max_cell_connectivity)
      DO jb = i_startblk, i_endblk
        CALL get_indices_v(p_patch, jb, i_startblk, i_endblk, &
          &                i_startidx, i_endidx, rl_start, rl_end)
        DO jc = i_startidx, i_endidx
          IF ((p_patch%verts%cell_idx(jc,jb,j) == 0) .AND. &
            & (p_patch%verts%refin_ctrl(jc,jb) /= 1)) THEN
            iidx = p_patch%verts%cell_idx(jc,jb,5)
            iblk = p_patch%verts%cell_blk(jc,jb,5)
            lonv(jc,jb,7-j) = p_patch%cells%center(iidx,iblk)%lon
            latv(jc,jb,7-j) = p_patch%cells%center(iidx,iblk)%lat
          ELSE IF ((p_patch%verts%cell_idx(jc,jb,j) < 0) .OR. &
            &      (p_patch%verts%refin_ctrl(jc,jb) == 1)) THEN
            lonv(jc,jb,7-j) = 0._wp
            latv(jc,jb,7-j) = 0._wp
          ELSE
            iidx = p_patch%verts%cell_idx(jc,jb,j)
            iblk = p_patch%verts%cell_blk(jc,jb,j)
            lonv(jc,jb,7-j) = p_patch%cells%center(iidx,iblk)%lon
            latv(jc,jb,7-j) = p_patch%cells%center(iidx,iblk)%lat
          ENDIF
        ENDDO
      ENDDO
    END DO
  END SUBROUTINE cf_1_1_grid_verts


  !------------------------------------------------------------------------------------------------
  !> Sets the grid information in output file
  !
  SUBROUTINE set_grid_info_netcdf(gridID, grid_info)
    INTEGER,             INTENT(IN) :: gridID
    TYPE(t_grid_info),   INTENT(IN) :: grid_info
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::set_grid_info_netcdf"

    CALL gridDefXvals  (gridID, grid_info%lon)
    CALL gridDefYvals  (gridID, grid_info%lat)
    CALL gridDefXbounds(gridID, grid_info%lonv)
    CALL gridDefYbounds(gridID, grid_info%latv)
  END SUBROUTINE set_grid_info_netcdf


  !------------------------------------------------------------------------------------------------
  !> Declaration of the grid information (RLAT/RLON) in output file, GRIB2 format.
  !
  SUBROUTINE set_grid_info_grb2(of)
    TYPE (t_output_file), INTENT(INOUT) :: of
    ! local variables:
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::set_grid_info_grb"
    CHARACTER(LEN=4), PARAMETER :: grid_coord_name(2) = (/ "RLON", "RLAT" /)
    TYPE (t_grib2_var), PARAMETER :: grid_coord_grib2(2) = (/  &
      ! geographical longitude RLON
      & t_grib2_var(               0,   &  ! discipline
      &                          191,   &  ! category
      &                            2,   &  ! number
      &              DATATYPE_PACK16,   &  ! bits
      &               GRID_REFERENCE,   &  ! gridtype
      &                    GRID_CELL ), &  ! subgridtype
      ! geographical latitude RLAT
      & t_grib2_var(               0,   &  ! discipline
      &                          191,   &  ! category
      &                            1,   &  ! number
      &              DATATYPE_PACK16,   &  ! bits
      &               GRID_REFERENCE,   &  ! gridtype
      &                    GRID_CELL )  &  ! subgridtype
      /)

    INTEGER :: igrid,i,vlistID,idx(3),gridID(3),zaxisID

    vlistID = of%cdiVlistID
    zaxisID = of%cdiZaxisID(ZA_surface)

    SELECT CASE(of%name_list%remap)
    CASE (REMAP_NONE)
      idx(:)    = (/ ICELL, IEDGE, IVERT /)
      gridID(:) = (/ of%cdiCellGridID, of%cdiEdgeGridID, of%cdiVertGridID /)
      DO igrid=1,3
        DO i=1,2 ! for longitude, latitude:
          of%cdi_grb2(idx(igrid),i) = vlistDefVar(vlistID, gridID(igrid), zaxisID, TSTEP_CONSTANT)
          CALL vlistDefVarDatatype(vlistID,  of%cdi_grb2(idx(igrid),i), grid_coord_grib2(i)%bits)
          CALL vlistDefVarTsteptype(vlistID, of%cdi_grb2(idx(igrid),i), TSTEP_CONSTANT)
          CALL vlistDefVarName(vlistID,      of%cdi_grb2(idx(igrid),i), TRIM(grid_coord_name(i)))

          ! Set GRIB2 Triplet
          CALL vlistDefVarParam( vlistID, of%cdi_grb2(idx(igrid),i),  &
            &  cdiEncodeParam(grid_coord_grib2(i)%number,             &
            &                 grid_coord_grib2(i)%category,           &
            &                 grid_coord_grib2(i)%discipline) )

          ! GRIB2 Quick hack: Set additional GRIB2 keys
          CALL set_additional_GRIB2_keys(vlistID, of%cdi_grb2(idx(igrid),i),   &
            &                            gribout_config(of%phys_patch_id), 0 )
        END DO
      END DO

    CASE (REMAP_REGULAR_LATLON)
      ! for longitude, latitude:
      DO i=1,2
        of%cdi_grb2(ILATLON,i) = vlistDefVar(vlistID, of%cdiLonLatGridID, zaxisID, TSTEP_CONSTANT)
        CALL vlistDefVarDatatype(vlistID,  of%cdi_grb2(ILATLON,i), grid_coord_grib2(i)%bits)
        CALL vlistDefVarTsteptype(vlistID, of%cdi_grb2(ILATLON,i), TSTEP_CONSTANT)
        CALL vlistDefVarName(vlistID,      of%cdi_grb2(ILATLON,i), TRIM(grid_coord_name(i)))
        ! Set GRIB2 Triplet
        CALL vlistDefVarParam( vlistID, of%cdi_grb2(ILATLON,i),   &
          &  cdiEncodeParam(grid_coord_grib2(i)%number,            &
          &                 grid_coord_grib2(i)%category,          &
          &                 grid_coord_grib2(i)%discipline) )

        ! GRIB2 Quick hack: Set additional GRIB2 keys
        CALL set_additional_GRIB2_keys(vlistID, of%cdi_grb2(ILATLON,i),      &
          &                            gribout_config(of%phys_patch_id), 0 )
      END DO

    CASE DEFAULT
      CALL finish(routine, "Unsupported grid type.")
    END SELECT

  END SUBROUTINE set_grid_info_grb2


  !------------------------------------------------------------------------------------------------
  !
  ! Copies the grid information from grid file to output file
  !
  SUBROUTINE copy_grid_info(of, patch_info)
    TYPE(t_patch_info),   INTENT(INOUT) :: patch_info (:)
    TYPE (t_output_file), INTENT(IN) :: of

    INTEGER :: ncid, dimid, varid
    INTEGER :: i_nc, i_ne, i_nv, max_cell_connectivity
    INTEGER :: i_dom

    REAL(wp), ALLOCATABLE :: clon(:), clat(:), clonv(:,:), clatv(:,:)
    REAL(wp), ALLOCATABLE :: elon(:), elat(:), elonv(:,:), elatv(:,:)
    REAL(wp), ALLOCATABLE :: vlon(:), vlat(:), vlonv(:,:), vlatv(:,:)

    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::copy_grid_info"

    i_dom = of%phys_patch_id
    ! Please note: The following is more or less a copy from mo_io_vlist with adaptions
    ! to the data structures used here.
    ! Unfortunately it seems necessary to open the gridfile for reading the information
    ! since it is not read and stored during patch input.

    !---------------------------------------------------------------------------
    ! Open grid file, read dimensions and make a cross check if they match.
    ! This is just for safety and could be skipped, of course.

    CALL nf(nf_open(TRIM(patch_info(i_dom)%grid_filename), NF_NOWRITE, ncid))

    CALL nf(nf_inq_dimid(ncid, 'nv', dimid))
    CALL nf(nf_inq_dimlen(ncid, dimid, max_cell_connectivity))
    !
    ! SELECT CASE (max_cell_connectivity)
    ! CASE (3)
    CALL nf(nf_inq_dimid(ncid, 'cell', dimid))
    ! CASE (6)
    !  CALL nf(nf_inq_dimid(ncid, 'vertex', dimid))
    ! END SELECT
    CALL nf(nf_inq_dimlen(ncid, dimid, i_nc))
    !
    CALL nf(nf_inq_dimid(ncid, 'edge', dimid))
    CALL nf(nf_inq_dimlen(ncid, dimid, i_ne))
    !
    ! SELECT CASE (max_cell_connectivity)
    ! CASE (3)
      CALL nf(nf_inq_dimid(ncid, 'vertex', dimid))
    ! CASE (6)
      CALL nf(nf_inq_dimid(ncid, 'cell', dimid))
    ! END SELECT

    CALL nf(nf_inq_dimlen(ncid, dimid, i_nv))

    IF(i_nc /= patch_info(i_dom)%cells%n_log) &
      CALL finish(routine,'Number of cells differs in '//TRIM(patch_info(i_dom)%grid_filename))
    IF(i_ne /= patch_info(i_dom)%edges%n_log) &
      CALL finish(routine,'Number of edges differs in '//TRIM(patch_info(i_dom)%grid_filename))
    IF(i_nv /= patch_info(i_dom)%verts%n_log) &
      CALL finish(routine,'Number of verts differs in '//TRIM(patch_info(i_dom)%grid_filename))
    !
    !---------------------------------------------------------------------------
    ! cell grid

    ! SELECT CASE (max_cell_connectivity)
    ! CASE (3)
    CALL nf(nf_inq_varid(ncid, 'clon', varid))
    ! CASE (6)
    !   CALL nf(nf_inq_varid(ncid, 'vlon', varid))
    ! END SELECT

    ALLOCATE(clon(i_nc))
    CALL nf(nf_get_var_double(ncid, varid, clon))
    CALL reorder1(patch_info(i_dom)%cells%n_glb, &
      &           patch_info(i_dom)%cells%grid_info%log_dom_index, clon)
    CALL gridDefXvals(of%cdiCellGridID, clon)
    DEALLOCATE(clon)


    ! SELECT CASE (max_cell_connectivity)
    ! CASE (3)
    CALL nf(nf_inq_varid(ncid, 'clat', varid))
    ! CASE (6)
    !  CALL nf(nf_inq_varid(ncid, 'vlat', varid))
    ! END SELECT

    ALLOCATE(clat(i_nc))
    CALL nf(nf_get_var_double(ncid, varid, clat))
    CALL reorder1(patch_info(i_dom)%cells%n_glb, &
      &           patch_info(i_dom)%cells%grid_info%log_dom_index, clat)

    CALL gridDefYvals(of%cdiCellGridID, clat)
    DEALLOCATE(clat)


    ! SELECT CASE (max_cell_connectivity)
    ! CASE (3)
    CALL nf(nf_inq_varid(ncid, 'clon_vertices', varid))
    ! CASE (6)
    ! CALL nf(nf_inq_varid(ncid, 'vlon_vertices', varid))
    ! END SELECT

    ALLOCATE(clonv(max_cell_connectivity, i_nc))
    CALL nf(nf_get_var_double(ncid, varid, clonv))
    CALL reorder2(patch_info(i_dom)%cells%n_glb, &
      &           patch_info(i_dom)%cells%grid_info%log_dom_index, clonv)

    CALL gridDefXbounds(of%cdiCellGridID, clonv)
    DEALLOCATE(clonv)


    ! SELECT CASE (max_cell_connectivity)
    ! CASE (3)
    CALL nf(nf_inq_varid(ncid, 'clat_vertices', varid))
    !CASE (6)
    !  CALL nf(nf_inq_varid(ncid, 'vlat_vertices', varid))
    ! END SELECT

    ALLOCATE(clatv(max_cell_connectivity, i_nc))
    CALL nf(nf_get_var_double(ncid, varid, clatv))
    CALL reorder2(patch_info(i_dom)%cells%n_glb, &
      &           patch_info(i_dom)%cells%grid_info%log_dom_index, clatv)

    CALL gridDefYbounds(of%cdiCellGridID, clatv)
    DEALLOCATE(clatv)

    !-------------------------------------------------------------------------
    ! edge grid

    ALLOCATE(elon(i_ne))
    CALL nf(nf_inq_varid(ncid, 'elon', varid))
    CALL nf(nf_get_var_double(ncid, varid, elon))
    CALL reorder1(patch_info(i_dom)%edges%n_glb, &
      &           patch_info(i_dom)%edges%grid_info%log_dom_index, elon)

    CALL gridDefXvals(of%cdiEdgeGridID, elon)
    DEALLOCATE(elon)

    ALLOCATE(elat(i_ne))
    CALL nf(nf_inq_varid(ncid, 'elat', varid))
    CALL nf(nf_get_var_double(ncid, varid, elat))
    CALL reorder1(patch_info(i_dom)%edges%n_glb, &
      &           patch_info(i_dom)%edges%grid_info%log_dom_index, elat)

    CALL gridDefYvals(of%cdiEdgeGridID, elat)
    DEALLOCATE(elat)

    ALLOCATE(elonv(4, i_ne))
    CALL nf(nf_inq_varid(ncid, 'elon_vertices', varid))
    CALL nf(nf_get_var_double(ncid, varid, elonv))
    CALL reorder2(patch_info(i_dom)%edges%n_glb, &
      &           patch_info(i_dom)%edges%grid_info%log_dom_index, elonv)

    CALL gridDefXbounds(of%cdiEdgeGridID, elonv)
    DEALLOCATE(elonv)

    ALLOCATE(elatv(4, i_ne))
    CALL nf(nf_inq_varid(ncid, 'elat_vertices', varid))
    CALL nf(nf_get_var_double(ncid, varid, elatv))
    CALL reorder2(patch_info(i_dom)%edges%n_glb, &
      &           patch_info(i_dom)%edges%grid_info%log_dom_index, elatv)

    CALL gridDefYbounds(of%cdiEdgeGridID, elatv)
    DEALLOCATE(elatv)

    !-------------------------------------------------------------------------
    ! vertex grid

    ! SELECT CASE (max_cell_connectivity)
    ! CASE (3)
    CALL nf(nf_inq_varid(ncid, 'vlon', varid))
    ! CASE (6)
    !  CALL nf(nf_inq_varid(ncid, 'clon', varid))
    !END SELECT

    ALLOCATE(vlon(i_nv))
    CALL nf(nf_get_var_double(ncid, varid, vlon))
    CALL reorder1(patch_info(i_dom)%verts%n_glb, &
      &           patch_info(i_dom)%verts%grid_info%log_dom_index, vlon)

    CALL gridDefXvals(of%cdiVertGridID, vlon)
    DEALLOCATE(vlon)

    ! SELECT CASE (max_cell_connectivity)
    ! CASE (3)
      CALL nf(nf_inq_varid(ncid, 'vlat', varid))
    ! CASE (6)
    !  CALL nf(nf_inq_varid(ncid, 'clat', varid))
    !END SELECT

    ALLOCATE(vlat(i_nv))
    CALL nf(nf_get_var_double(ncid, varid, vlat))
    CALL reorder1(patch_info(i_dom)%verts%n_glb, &
      &           patch_info(i_dom)%verts%grid_info%log_dom_index, vlat)

    CALL gridDefYvals(of%cdiVertGridID, vlat)
    DEALLOCATE(vlat)

    !IF(max_cell_connectivity==3) THEN
      CALL nf(nf_inq_varid(ncid, 'vlon_vertices', varid))
    !ELSEIF(max_cell_connectivity==6) THEN
    !  CALL nf(nf_inq_varid(ncid, 'clon_vertices', varid))
    !ENDIF

    ALLOCATE(vlonv(9-max_cell_connectivity, i_nv))
    CALL nf(nf_get_var_double(ncid, varid, vlonv))
    CALL reorder2(patch_info(i_dom)%verts%n_glb, &
      &           patch_info(i_dom)%verts%grid_info%log_dom_index, vlonv)

    CALL gridDefXbounds(of%cdiVertGridID, vlonv)
    DEALLOCATE(vlonv)

    ! IF(max_cell_connectivity==3) THEN
    CALL nf(nf_inq_varid(ncid, 'vlat_vertices', varid))
    ! ELSEIF(max_cell_connectivity==6) THEN
    !  CALL nf(nf_inq_varid(ncid, 'clat_vertices', varid))
    !ENDIF

    ALLOCATE(vlatv(9-max_cell_connectivity, i_nv))
    CALL nf(nf_get_var_double(ncid, varid, vlatv))
    CALL reorder2(patch_info(i_dom)%verts%n_glb, &
      &           patch_info(i_dom)%verts%grid_info%log_dom_index, vlatv)

    CALL gridDefYbounds(of%cdiVertGridID, vlatv)
    DEALLOCATE(vlatv)

    !-------------------------------------------------------------------------

    ! Close NetCDF file, it is not needed any more
    CALL nf(nf_close(ncid))

  CONTAINS

    SUBROUTINE nf(status)
      INTEGER, INTENT(in) :: status

      IF (status /= nf_noerr) THEN
        CALL finish(routine, 'NetCDF Error: '//nf_strerror(status))
      ENDIF
    END SUBROUTINE nf

    ! reorder1: get the physical patch points from the logical patch
    ! Note that this works within the array as long as idx is monotonically increasing
    SUBROUTINE reorder1(n, idx, array)
      INTEGER, INTENT(IN)     :: n, idx(:)
      REAL(wp), INTENT(INOUT) :: array(:)
      INTEGER :: i

      DO i = 1, n
        array(i) = array(idx(i))
      ENDDO
    END SUBROUTINE reorder1

    ! reorder2: same as reorder1 for 2D array
    SUBROUTINE reorder2(n, idx, array)
      INTEGER, INTENT(IN)     :: n, idx(:)
      REAL(wp), INTENT(INOUT) :: array(:,:)
      INTEGER :: i

      DO i = 1, n
        array(:,i) = array(:,idx(i))
      ENDDO
    END SUBROUTINE reorder2

  END SUBROUTINE copy_grid_info

  
  SUBROUTINE bcast_grid_info(patch_info, bcast_root)
    TYPE(t_patch_info), INTENT(INOUT) :: patch_info
    INTEGER,            intent(IN)    :: bcast_root
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::bcast_grid_info"
    INTEGER :: ierrstat, dim_c, dim_e, dim_v

    dim_c = patch_info%max_cell_connectivity
    dim_e =                  4
    dim_v = 9 - patch_info%max_cell_connectivity

    IF(my_process_is_io()) THEN
      ALLOCATE(patch_info%cells%grid_info%lon(nproma*patch_info%nblks_glb_c),         &
        &      patch_info%cells%grid_info%lat(nproma*patch_info%nblks_glb_c),         &
        &      patch_info%cells%grid_info%lonv(nproma*patch_info%nblks_glb_c, dim_c), &
        &      patch_info%cells%grid_info%latv(nproma*patch_info%nblks_glb_c, dim_c), &
        !
        &      patch_info%edges%grid_info%lon(nproma*patch_info%nblks_glb_e),         &
        &      patch_info%edges%grid_info%lat(nproma*patch_info%nblks_glb_e),         &
        &      patch_info%edges%grid_info%lonv(nproma*patch_info%nblks_glb_e, dim_e), &
        &      patch_info%edges%grid_info%latv(nproma*patch_info%nblks_glb_e, dim_e), &
        !
        &      patch_info%verts%grid_info%lon(nproma*patch_info%nblks_glb_v),         &
        &      patch_info%verts%grid_info%lat(nproma*patch_info%nblks_glb_v),         &
        &      patch_info%verts%grid_info%lonv(nproma*patch_info%nblks_glb_v, dim_v), &
        &      patch_info%verts%grid_info%latv(nproma*patch_info%nblks_glb_v, dim_v), &
        &      STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')
    END IF

    ! cells
    CALL p_bcast(patch_info%cells%grid_info%lon,  bcast_root, p_comm_work_2_io)
    CALL p_bcast(patch_info%cells%grid_info%lat,  bcast_root, p_comm_work_2_io)
    CALL p_bcast(patch_info%cells%grid_info%lonv, bcast_root, p_comm_work_2_io)
    CALL p_bcast(patch_info%cells%grid_info%latv, bcast_root, p_comm_work_2_io)
    ! edges
    CALL p_bcast(patch_info%edges%grid_info%lon,  bcast_root, p_comm_work_2_io)
    CALL p_bcast(patch_info%edges%grid_info%lat,  bcast_root, p_comm_work_2_io)
    CALL p_bcast(patch_info%edges%grid_info%lonv, bcast_root, p_comm_work_2_io)
    CALL p_bcast(patch_info%edges%grid_info%latv, bcast_root, p_comm_work_2_io)
    ! vertices
    CALL p_bcast(patch_info%verts%grid_info%lon,  bcast_root, p_comm_work_2_io)
    CALL p_bcast(patch_info%verts%grid_info%lat,  bcast_root, p_comm_work_2_io)
    CALL p_bcast(patch_info%verts%grid_info%lonv, bcast_root, p_comm_work_2_io)
    CALL p_bcast(patch_info%verts%grid_info%latv, bcast_root, p_comm_work_2_io)
  END SUBROUTINE bcast_grid_info


  !------------------------------------------------------------------------------------------------
  !> Writes the grid information in output file, GRIB2 format.
  !
! #ifndef __NO_ICON_ATMO__
  SUBROUTINE write_grid_info_grb2(of, patch_info)
    TYPE (t_output_file), INTENT(INOUT)           :: of
    TYPE(t_patch_info),   INTENT(IN),   TARGET    :: patch_info (:)
    ! local variables:
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::write_grid_info_grb2"

    TYPE t_grid_info_ptr
      TYPE (t_grid_info), POINTER :: ptr
    END TYPE t_grid_info_ptr

    INTEGER                        :: errstat, idom, igrid, idx(3), isize(3), idom_log
    TYPE (t_lon_lat_grid), POINTER :: grid
    REAL(wp), ALLOCATABLE          :: rotated_pts(:,:,:), r_out_dp(:,:), r_out_dp_1D(:)
    TYPE(t_grid_info_ptr)          :: gptr(3)

    ! skip this on test PE...
    IF (my_process_is_mpi_test()) RETURN

    SELECT CASE(of%name_list%remap)
    CASE (REMAP_NONE)
      idom     = of%phys_patch_id
      idom_log = patch_info(idom)%log_patch_id
      idx(:)    = (/ ICELL, IEDGE, IVERT /)
      isize(:)  = (/ patch_info(idom_log)%cells%n_glb, &
        &            patch_info(idom_log)%edges%n_glb, &
        &            patch_info(idom_log)%verts%n_glb /)
      gptr(1)%ptr => patch_info(idom_log)%cells%grid_info
      gptr(2)%ptr => patch_info(idom_log)%edges%grid_info
      gptr(3)%ptr => patch_info(idom_log)%verts%grid_info
      DO igrid=1,3
        ! allocate data buffer:
        ALLOCATE(r_out_dp_1D(isize(igrid)), stat=errstat)
        IF (errstat /= SUCCESS) CALL finish(routine, 'ALLOCATE failed!')
        ! write RLON, RLAT
        r_out_dp_1D(:) = gptr(igrid)%ptr%lon(1:isize(igrid)) / pi_180
        CALL streamWriteVar(of%cdiFileID, of%cdi_grb2(idx(igrid),IRLON), r_out_dp_1D, 0)
        r_out_dp_1D(:) = gptr(igrid)%ptr%lat(1:isize(igrid)) / pi_180
        CALL streamWriteVar(of%cdiFileID, of%cdi_grb2(idx(igrid),IRLAT), r_out_dp_1D, 0)
        ! clean up
        DEALLOCATE(r_out_dp_1D, stat=errstat)
        IF (errstat /= SUCCESS) CALL finish(routine, 'DEALLOCATE failed!')
      END DO

    CASE (REMAP_REGULAR_LATLON)
      ! allocate data buffer:
      grid => lonlat_grid_list(of%name_list%lonlat_id)%grid
      ! compute some entries of lon-lat grid specification:
      CALL compute_lonlat_specs(grid)
      ALLOCATE(rotated_pts(grid%lon_dim, grid%lat_dim, 2), &
        &      r_out_dp(grid%lon_dim,grid%lat_dim), stat=errstat)
      IF (errstat /= SUCCESS) CALL finish(routine, 'ALLOCATE failed!')
      ! compute grid points of rotated lon/lat grid
      CALL rotate_latlon_grid(grid, rotated_pts)
      ! write RLON, RLAT
      r_out_dp(:,:) = rotated_pts(:,:,1) / pi_180
      CALL streamWriteVar(of%cdiFileID, of%cdi_grb2(ILATLON,IRLON), r_out_dp, 0)
      r_out_dp(:,:) = rotated_pts(:,:,2) / pi_180
      CALL streamWriteVar(of%cdiFileID, of%cdi_grb2(ILATLON,IRLAT), r_out_dp, 0)
      ! clean up
      DEALLOCATE(rotated_pts, r_out_dp, stat=errstat)
      IF (errstat /= SUCCESS) CALL finish(routine, 'DEALLOCATE failed!')

    CASE DEFAULT
      CALL finish(routine, "Unsupported grid type.")
    END SELECT

  END SUBROUTINE write_grid_info_grb2
! #endif

END MODULE mo_name_list_output_gridinfo
