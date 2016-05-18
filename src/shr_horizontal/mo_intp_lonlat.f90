#ifdef __xlC__
  @PROCESS smp=noopt
  @PROCESS noopt
#endif
#ifdef __PGI
  !pgi$g opt=1
#endif

!>
!! Contains the implementation of interpolation onto regular grids.
!!
!! @par Revision History
!! Moved from mo_intp_rbf_coeffs : 2012-03-20, F. Prill (DWD)
!! Modified by Anurag Dipankar, MPIM, 2012-12-28
!!-Replaced usage of ptr_int%cart_edge%coord with ptr_patch%edges%cartesian_center which is now calculated
!! within the grid_generator. The ptr_int variable is not calculated anymore.
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
!! TODO[FP] Move, if possible, lonlat interpolation data structure from
!!          "mo_intp_data_strc" into this module.
!!
  MODULE mo_intp_lonlat
    !-------------------------------------------------------------------------
    !
    !    ProTeX FORTRAN source: Style 2
    !    modified for ICON project, DWD/MPI-M 2006
    !
    !-------------------------------------------------------------------------
    !
!$  USE OMP_LIB
    USE mo_kind,                ONLY: wp
    USE mo_exception,           ONLY: message, message_text, finish
    USE mo_impl_constants,      ONLY: SUCCESS, min_rlcell_int,                                &
      &                               HINTP_TYPE_NONE, HINTP_TYPE_LONLAT_RBF,                 &
      &                               HINTP_TYPE_LONLAT_NNB, HINTP_TYPE_LONLAT_BCTR,          &
      &                               SCALE_MODE_TABLE, SCALE_MODE_AUTO, SCALE_MODE_PRESET
    USE mo_impl_constants_grf,  ONLY: grf_bdywidth_c
    USE mo_model_domain,        ONLY: t_patch
    USE mo_run_config,          ONLY: ltimer
    USE mo_grid_config,         ONLY: n_dom, grid_sphere_radius, is_plane_torus
    USE mo_timer,               ONLY: timer_start, timer_stop, timer_lonlat_setup
    USE mo_math_utilities,      ONLY: gc2cc, gvec2cvec, arc_length_v,                         &
      &                               t_cartesian_coordinates,                                &
      &                               t_geographical_coordinates
    USE mo_physical_constants,  ONLY: earth_radius
    USE mo_math_utility_solvers, ONLY: solve_chol_v, choldec_v
    USE mo_lonlat_grid,         ONLY: t_lon_lat_grid, latlon_compute_area_weights
    USE mo_parallel_config,     ONLY: nproma, p_test_run
    USE mo_io_config,           ONLY: lnetcdf_flt64_output
    USE mo_loopindices,         ONLY: get_indices_c, get_indices_e
    USE mo_intp_data_strc,      ONLY: t_int_state, t_lon_lat_intp, n_lonlat_grids,            &
      &                               lonlat_grid_list, n_lonlat_grids, MAX_LONLAT_GRIDS
    USE mo_interpol_config,     ONLY: rbf_vec_dim_c, rbf_c2grad_dim, rbf_vec_kern_ll,         &
      &                               rbf_vec_scale_ll, rbf_dim_c2l,                          &
      &                               l_mono_c2l, rbf_scale_mode_ll, support_baryctr_intp
    USE mo_gnat_gridsearch,     ONLY: gnat_init_grid, gnat_destroy, t_gnat_tree,              &
      &                               gnat_query_containing_triangles,                        &
      &                               gnat_merge_distributed_queries, gk, SKIP_NODE,          &
      &                               INVALID_NODE, gnat_recursive_proximity_query
    USE mo_mpi,                 ONLY: get_my_mpi_work_id,                                     &
      &                               p_min, p_comm_work,                                     &
      &                               my_process_is_mpi_test, p_max, p_send,                  &
      &                               p_recv, process_mpi_all_test_id,                        &
      &                               process_mpi_all_workroot_id, p_pe,                      &
      &                               my_process_is_stdio
    USE mo_communication,       ONLY: blk_no, idx_no, setup_comm_gather_pattern, idx_1d
    USE mo_lonlat_grid,         ONLY: rotate_latlon_grid
    USE mo_cf_convention,       ONLY: t_cf_var
    USE mo_grib2,               ONLY: t_grib2_var, grib2_var
    USE mo_cdi,                 ONLY: TSTEP_CONSTANT, GRID_UNSTRUCTURED,                      &
      &                               DATATYPE_PACK16, DATATYPE_FLT32, DATATYPE_FLT64
    USE mo_cdi_constants,       ONLY: GRID_REGULAR_LONLAT, GRID_CELL, ZA_SURFACE
    USE mo_nonhydro_state,      ONLY: p_nh_state_lists
    USE mo_var_list,            ONLY: add_var
    USE mo_var_metadata,        ONLY: create_hor_interp_metadata
    USE mo_linked_list,         ONLY: t_list_element
    USE mo_sync,                ONLY: SYNC_C, sync_idx, sync_patch_array
    USE mo_rbf_errana,          ONLY: estimate_rbf_parameter
    USE mo_delaunay_types,      ONLY: t_point_list, OPERATOR(/), t_triangulation, OPERATOR(+)
    USE mo_intp_lonlat_baryctr, ONLY: baryctr_interpol_lonlat, setup_barycentric_intp_lonlat, &
      &                               setup_barycentric_intp_lonlat_repartition,              &
      &                               compute_auxiliary_triangulation
    USE mo_util_string,         ONLY: int2string

    IMPLICIT NONE

    !> level of output verbosity (for debugging purposes)
    INTEGER, PARAMETER  :: dbg_level = 0

    CHARACTER(LEN=*), PARAMETER :: modname = 'mo_intp_lonlat'

    PRIVATE

    PUBLIC :: init_lonlat_grid_list
    PUBLIC :: destroy_lonlat_grid_list
    PUBLIC :: compute_lonlat_intp_coeffs
    PUBLIC :: compute_lonlat_area_weights
    PUBLIC :: rbf_vec_interpol_lonlat
    PUBLIC :: interpol_lonlat

    INTERFACE interpol_lonlat
      MODULE PROCEDURE interpol_lonlat_int
      MODULE PROCEDURE interpol_lonlat_real
    END INTERFACE

    INTERFACE nnb_interpol_lonlat
      MODULE PROCEDURE nnb_interpol_lonlat_real
      MODULE PROCEDURE nnb_interpol_lonlat_int
    END INTERFACE

  CONTAINS

#include "intp_functions.inc"

    !===============================================================
    ! SETUP OF LON-LAT GRID LIST

    !---------------------------------------------------------------
    !> Setup of lon-lat registry
    !
    SUBROUTINE init_lonlat_grid_list()
      INTEGER :: i
      CHARACTER(*), PARAMETER :: routine = modname//"::init_lonlat_grid_list"

      IF (dbg_level > 5) CALL message(routine, "Enter")

      ! not much to do yet...
      n_lonlat_grids = 0
      DO i=1,MAX_LONLAT_GRIDS
        lonlat_grid_list(i)%l_dom(:)         = .FALSE.
        lonlat_grid_list(i)%l_initialized(:) = .FALSE.
      END DO
      IF (dbg_level > 5) CALL message(routine, "Done")
    END SUBROUTINE init_lonlat_grid_list


    !---------------------------------------------------------------
    !> Frees all data allocated by lon-lat grid list
    !
    SUBROUTINE destroy_lonlat_grid_list()
      INTEGER :: i, jg
      CHARACTER(*), PARAMETER :: routine = modname//"::destroy_lonlat_grid_list"

      IF (dbg_level > 5) CALL message(routine, "Enter")
      DO i=1, n_lonlat_grids
        DO jg=1,n_dom
          IF (lonlat_grid_list(i)%l_initialized(jg)) THEN
            CALL deallocate_int_state_lonlat(lonlat_grid_list(i)%intp(jg))
            lonlat_grid_list(i)%l_initialized(jg) = .FALSE.
          END IF
        END DO
      END DO
      IF (dbg_level > 5) CALL message(routine, "Done")

    END SUBROUTINE destroy_lonlat_grid_list


    !---------------------------------------------------------------
    !> Resize arrays with lon-lat interpolation data from global to
    !  local size. This can be performed as soon as the setup phase is
    !  finished.
    SUBROUTINE resize_lonlat_state(ptr_int_lonlat, nblks_local, nlocal_pts)
      TYPE (t_lon_lat_intp), TARGET, INTENT(INOUT) :: ptr_int_lonlat
      INTEGER,                       INTENT(IN)    :: nblks_local, nlocal_pts

      ! local variables
      CHARACTER(*), PARAMETER :: routine = modname//"::resize_lonlat_state"
      INTEGER, ALLOCATABLE  :: tmp_tri_idx(:,:,:), tmp_global_idx(:)
      INTEGER               :: errstat

      ! perform a triangle copy for the two fields:
      IF (dbg_level > 5) CALL message(routine, "Enter")

      ! first allocate temporary storage and copy fields:
      ALLOCATE(tmp_tri_idx(2, nproma, nblks_local),   &
        &      tmp_global_idx(nlocal_pts), STAT=errstat )
      IF (errstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed')
      tmp_tri_idx(:,:,1:nblks_local) = ptr_int_lonlat%tri_idx(:,:,1:nblks_local)
      tmp_global_idx(1:nlocal_pts)   = ptr_int_lonlat%global_idx(1:nlocal_pts)

      ! resize tri_idx and global_idx:
      DEALLOCATE(ptr_int_lonlat%tri_idx, ptr_int_lonlat%global_idx, stat=errstat)
      IF (errstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed')
      ALLOCATE(ptr_int_lonlat%tri_idx(2, nproma, nblks_local),   &
        &      ptr_int_lonlat%global_idx(nlocal_pts), STAT=errstat )
      IF (errstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed')

      ! copy contents back:
      ptr_int_lonlat%tri_idx(:,:,1:nblks_local) = tmp_tri_idx(:,:,1:nblks_local)
      ptr_int_lonlat%global_idx(1:nlocal_pts)   = tmp_global_idx(1:nlocal_pts)

      ! clean up
      DEALLOCATE(tmp_tri_idx, tmp_global_idx, stat=errstat)
      IF (errstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed')

      IF (dbg_level > 5) CALL message(routine, "Done")

    END SUBROUTINE resize_lonlat_state


    !---------------------------------------------------------------
    !> Setup of lon-lat interpolation
    !
    SUBROUTINE compute_lonlat_intp_coeffs(p_patch, p_int_state)
      TYPE(t_patch),        INTENT(IN)    :: p_patch(:)
      TYPE(t_int_state),    INTENT(INOUT) :: p_int_state(:)
      ! local variables
      CHARACTER(*), PARAMETER :: routine = modname//"::compute_lonlat_intp_coeffs"
      LOGICAL,      PARAMETER :: lrepartition = .TRUE.

      INTEGER                :: i, jg, n_points, my_id, nthis_local_pts
      TYPE(t_gnat_tree)      :: gnat
      TYPE (t_triangulation) :: tri
      TYPE (t_point_list)    :: points
!$    DOUBLE PRECISION       :: time_s_total, toc

      ! -----------------------------------------------------------

      IF (dbg_level > 5) CALL message(routine, "Enter")
      DO i=1, n_lonlat_grids
        lonlat_grid_list(i)%l_dom(n_dom+1:) = .FALSE.
      END DO

      ! allocate and compute coefficients needed for lon-lat
      ! interpolation:
      DO jg=1,n_dom
        IF (.NOT. ANY(lonlat_grid_list(1:n_lonlat_grids)%l_dom(jg))) CYCLE

        ! build GNAT data structure
        IF (dbg_level > 1) CALL message(routine, "build GNAT data structure")
        CALL gnat_init_grid(gnat, p_patch(jg))

        DO i=1, n_lonlat_grids
          IF (lonlat_grid_list(i)%l_dom(jg)) THEN
            IF (ltimer) CALL timer_start(timer_lonlat_setup)
            ! compute weights for Radial Basis Function (RBF)
            ! interpolation:
            CALL rbf_setup_interpol_lonlat_grid(      &
              &         lonlat_grid_list(i)%grid,     &
              &         gnat,                         &
              &         p_patch(jg),                  &
              &         lonlat_grid_list(i)%intp(jg), &
              &         p_int_state(jg) )

            ! compute weights for barycentric interpolation.
            ! 
            ! note: this requires a previous call to
            ! "rbf_setup_interpol_lonlat_grid" where the interpolation
            ! data structure is allocated.
            !
            ! --------------------------------------------------------------------------
            !
            ! IMPORTANT: Currently, barycentric interpolation supported only
            !            - if the namelist parameter "interpol_nml/support_baryctr_intp"
            !              has been set to .TRUE.
            !
            ! --------------------------------------------------------------------------

            IF (support_baryctr_intp) THEN
!$            time_s_total = omp_get_wtime()
              IF ((jg > 1) .OR. (.NOT. lrepartition)) THEN
                CALL compute_auxiliary_triangulation(p_patch(jg), tri, points)
                CALL setup_barycentric_intp_lonlat(tri, points, lonlat_grid_list(i)%intp(jg))
                CALL tri%destructor()
                CALL points%destructor()
              ELSE
                ! for global domains there exists an alternative
                ! algorithm which repartitions the grid mass points
                ! before triangulating them:
                CALL setup_barycentric_intp_lonlat_repartition(p_patch(jg), lonlat_grid_list(i)%intp(jg))
              END IF
!$            toc = omp_get_wtime() - time_s_total
!$            IF (my_process_is_stdio() .and. (dbg_level > 5)) THEN
!$              WRITE (0,*) trim(routine), " :: total elapsed time: ", toc
!$            END IF
            END IF

            ! mask out nest boundary points
            IF (jg > 1) THEN 
              CALL mask_out_boundary( p_patch(jg), lonlat_grid_list(i)%intp(jg) )
            END IF

            IF (ltimer) CALL timer_stop(timer_lonlat_setup)
            lonlat_grid_list(i)%l_initialized(jg) = .TRUE.
          END IF
        END DO ! i

        ! clean up
        CALL gnat_destroy(gnat)
      END DO ! jg

      DO i=1, n_lonlat_grids
        ! set communication pattern for gathering lon-lat grids:
        DO jg=1,n_dom
          IF (lonlat_grid_list(i)%l_dom(jg)) THEN
            ! n_points     : total number of points in the RECEIVER array
            n_points = lonlat_grid_list(i)%grid%lon_dim * &
              &        lonlat_grid_list(i)%grid%lat_dim

            my_id = get_my_mpi_work_id()
            nthis_local_pts = lonlat_grid_list(i)%intp(jg)%nthis_local_pts

            ! resize global data structures, after the setup only
            ! local information is needed:
            CALL resize_lonlat_state(lonlat_grid_list(i)%intp(jg),     &
              &                      (nthis_local_pts - 1)/nproma + 1, &
              &                      nthis_local_pts)

            CALL setup_comm_gather_pattern(n_points,                      &
              (/(my_id, i = 1, nthis_local_pts)/),                        &
              lonlat_grid_list(i)%intp(jg)%global_idx(1:nthis_local_pts), &
              lonlat_grid_list(i)%p_pat(jg),.true.)

          END IF
        END DO ! jg

      END DO ! i
      IF (dbg_level > 5) CALL message(routine, "Done")

    END SUBROUTINE compute_lonlat_intp_coeffs


    !---------------------------------------------------------------
    !> Adds a special metrics variable containing the area weights of
    !  the regular lon-lat grid.
    !
    !  @note This new variable is time-constant!
    !
    SUBROUTINE compute_lonlat_area_weights()
      ! local variables
      CHARACTER(*), PARAMETER :: routine = modname//"::compute_lonlat_area_weights"
      TYPE(t_cf_var)       :: cf_desc
      TYPE(t_grib2_var)    :: grib2_desc
      INTEGER              :: var_shape(3), nblks_lonlat, i, jg, ierrstat, &
        &                     i_lat, idx_glb, jc, jb, j, datatype_flt
      TYPE (t_lon_lat_grid), POINTER :: grid
      TYPE (t_lon_lat_intp), POINTER :: ptr_int_lonlat
!DR      REAL(wp),              POINTER :: area_weights(:), p_dummy(:,:,:)
!DR !!! Using the POINTER attribute for area_weights(:) mysteriously leads
!DR !!! to a bus error on NEC SX9 (tested with compiler revision 450). However,
!DR !!! using the ALLOCATABLE attribute, instead, works.
      REAL(wp), ALLOCATABLE          :: area_weights(:)
      REAL(wp),              POINTER :: p_dummy(:,:,:)
      TYPE(t_list_element),  POINTER :: new_element

      ! define NetCDF output precision
      IF ( lnetcdf_flt64_output ) THEN
        datatype_flt = DATATYPE_FLT64
      ELSE
        datatype_flt = DATATYPE_FLT32
      ENDIF

      ! Add area weights
      DO i=1, n_lonlat_grids
        DO jg=1,n_dom
          IF (lonlat_grid_list(i)%l_dom(jg)) THEN

            grid           => lonlat_grid_list(i)%grid
            ptr_int_lonlat => lonlat_grid_list(i)%intp(jg)

            nblks_lonlat   =  (ptr_int_lonlat%nthis_local_pts - 1)/nproma + 1
            var_shape = (/ nproma, 1, nblks_lonlat /)
            cf_desc    = t_cf_var('aw', '1', 'area weights for regular lat-lon grid', datatype_flt)
            grib2_desc = grib2_var(0, 191, 193, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL)

            ALLOCATE(area_weights(grid%lat_dim), STAT=ierrstat)
            IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')

            CALL add_var( p_nh_state_lists(jg)%diag_list,                             &
              &           "aw", p_dummy,                                        &
              &           GRID_REGULAR_LONLAT, ZA_SURFACE, cf_desc, grib2_desc, &
              &           ldims=var_shape, lrestart=.FALSE.,                    &
              &           loutput=.TRUE., new_element=new_element,              &
              &           hor_interp=create_hor_interp_metadata(                &
              &             hor_intp_type=HINTP_TYPE_NONE ),                    &
              &           isteptype=TSTEP_CONSTANT )
            ! link this new variable to the lon-lat grid:
            new_element%field%info%hor_interp%lonlat_id = i
            ! compute area weights:
!CDIR NOIEXPAND
            CALL latlon_compute_area_weights(grid, earth_radius, area_weights)
            ! for each local lon-lat point on this PE:
            DO j=1, ptr_int_lonlat%nthis_local_pts
              ! determine block, index
              jb = (j-1)/nproma + 1
              jc = j - (jb-1)*nproma
              ! determine latitude index:
              idx_glb = ptr_int_lonlat%global_idx(j)
              i_lat   = (idx_glb-1)/grid%lon_dim + 1
              ! set area weight:
              p_dummy(jc,1,jb) = area_weights(i_lat)
            END DO

            DEALLOCATE(area_weights, STAT=ierrstat)
            IF (ierrstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed.')
          END IF
        END DO
      END DO

    END SUBROUTINE compute_lonlat_area_weights


    !-------------------------------------------------------------------------
    !
    !
    !> Allocation of components of interpolation state using an arbitrary grid.
    !!
    !! @par Revision History
    !! Initial implementation: F. Prill, DWD, 2011
    !! Introduced arbitrary grid, Rainer Johanni (2010-11-25)
    !!
    SUBROUTINE allocate_int_state_lonlat_grid( nblks_c, nblks_lonlat, ptr_int_lonlat )
      !
      INTEGER,                       INTENT(IN)    :: nblks_lonlat, nblks_c
      TYPE (t_lon_lat_intp), TARGET, INTENT(INOUT) :: ptr_int_lonlat

      CHARACTER(*), PARAMETER :: routine = modname//"::allocate_int_state_lonlat"
      INTEGER :: ist

      ! ----------------------------------------------------------------------

      ! allocate memory only when needed.
      ALLOCATE ( &
        &  ptr_int_lonlat%rbf_vec_coeff(rbf_vec_dim_c, 2, nproma, nblks_lonlat),     &
        &  ptr_int_lonlat%rbf_vec_idx(rbf_vec_dim_c, nproma, nblks_lonlat),          &
        &  ptr_int_lonlat%rbf_vec_blk(rbf_vec_dim_c, nproma, nblks_lonlat),          &
        &  ptr_int_lonlat%rbf_vec_stencil(nproma, nblks_lonlat),                     &
        &  STAT=ist )
      IF (ist /= SUCCESS) THEN
        CALL finish (routine, 'allocation for rbf lon-lat coeffs failed')
      ENDIF

      ptr_int_lonlat%rbf_vec_idx       = 0
      ptr_int_lonlat%rbf_vec_blk       = 0
      ptr_int_lonlat%rbf_vec_stencil   = 0
      ptr_int_lonlat%rbf_vec_coeff     = 0._wp

      ALLOCATE ( ptr_int_lonlat%rbf_c2l_coeff(rbf_dim_c2l, nproma, nblks_lonlat),         &
        &  ptr_int_lonlat%rbf_c2l_idx(rbf_dim_c2l, nproma, nblks_c),                      &
        &  ptr_int_lonlat%rbf_c2l_blk(rbf_dim_c2l, nproma, nblks_c),                      &
        &  ptr_int_lonlat%rbf_c2lr_idx(rbf_dim_c2l, nproma, nblks_lonlat),                &
        &  ptr_int_lonlat%rbf_c2lr_blk(rbf_dim_c2l, nproma, nblks_lonlat),                &
        &  ptr_int_lonlat%rbf_c2l_stencil(nproma, nblks_c),                               &
        &  ptr_int_lonlat%rbf_c2lr_stencil(nproma, nblks_lonlat),                         &
        STAT=ist )
      IF (ist /= SUCCESS) &
        CALL finish (routine, 'allocation for rbf lon-lat coeffs failed')
      
      ptr_int_lonlat%rbf_c2l_coeff    = 0._wp
      ptr_int_lonlat%rbf_c2l_idx      = 0
      ptr_int_lonlat%rbf_c2l_blk      = 0
      ptr_int_lonlat%rbf_c2lr_idx     = 0
      ptr_int_lonlat%rbf_c2lr_blk     = 0
      ptr_int_lonlat%rbf_c2l_stencil  = 0
      ptr_int_lonlat%rbf_c2lr_stencil = 0

      ! --- nearest-neighbor interpolation
      !
      !     (we need interpolation weights also here to blank out points)

      ALLOCATE (ptr_int_lonlat%nnb_coeff(1, nproma, nblks_lonlat), STAT=ist )
      IF (ist /= SUCCESS) THEN
        CALL finish (routine, 'allocation for nnb lon-lat coeffs failed')
      ENDIF
      ptr_int_lonlat%nnb_coeff = 1._wp

      ! --- barycentric interpolation

      ALLOCATE ( &
        &  ptr_int_lonlat%baryctr_coeff(3, nproma, nblks_lonlat), &
        &  ptr_int_lonlat%baryctr_idx(3, nproma, nblks_lonlat),   &
        &  ptr_int_lonlat%baryctr_blk(3, nproma, nblks_lonlat),   &
        &  STAT=ist )
      IF (ist /= SUCCESS) THEN
        CALL finish (routine, 'allocation for barycentric lon-lat coeffs failed')
      ENDIF
      ptr_int_lonlat%baryctr_coeff = 0._wp
      ptr_int_lonlat%baryctr_idx   = -1
      ptr_int_lonlat%baryctr_blk   = -1

      ! -- array with lon-lat coordinates
      ALLOCATE(ptr_int_lonlat%ll_coord(nproma, nblks_lonlat), stat=ist)
      IF (ist /= SUCCESS)  CALL finish (routine, 'Allocation of array with lon-lat coordinates!')

    END SUBROUTINE allocate_int_state_lonlat_grid


    !-------------------------------------------------------------------------
    !
    !>
    !! Deallocation of components of a single 2d lon-lat interpolation state.
    !!
    SUBROUTINE deallocate_int_state_lonlat( ptr_int_lonlat )

      TYPE (t_lon_lat_intp), TARGET, INTENT(INOUT) :: ptr_int_lonlat
      CHARACTER(*), PARAMETER :: routine = modname//"::deallocate_int_state_lonlat"
      INTEGER :: ist

      !-----------------------------------------------------------------------

      DEALLOCATE (ptr_int_lonlat%rbf_vec_coeff,           &
        &         ptr_int_lonlat%rbf_vec_idx,             &
        &         ptr_int_lonlat%rbf_vec_blk,             &
        &         ptr_int_lonlat%rbf_vec_stencil,         &
        &         ptr_int_lonlat%tri_idx,                 &
        &         ptr_int_lonlat%global_idx,              &
        &         STAT=ist )
      IF (ist /= SUCCESS) THEN
        CALL finish (routine, 'deallocation for lon-lat coefficients failed')
      ENDIF

      DEALLOCATE ( ptr_int_lonlat%rbf_c2l_coeff,    &
        &          ptr_int_lonlat%rbf_c2l_idx,      &
        &          ptr_int_lonlat%rbf_c2l_blk,      &
        &          ptr_int_lonlat%rbf_c2lr_idx,     &
        &          ptr_int_lonlat%rbf_c2lr_blk,     &
        &          ptr_int_lonlat%rbf_c2l_stencil,  &
        &          ptr_int_lonlat%rbf_c2lr_stencil, &
        &          STAT=ist )
      IF (ist /= SUCCESS) &
        CALL finish (routine, 'deallocation for rbf lon-lat coeffs failed')

      IF (ALLOCATED(ptr_int_lonlat%nnb_coeff)) THEN
        DEALLOCATE (ptr_int_lonlat%nnb_coeff, STAT=ist)
        IF (ist /= SUCCESS) THEN
          CALL finish (routine, 'deallocation for nnb lon-lat coefficients failed')
        ENDIF
      END IF

      IF (ALLOCATED(ptr_int_lonlat%baryctr_idx)) THEN
        DEALLOCATE (ptr_int_lonlat%baryctr_idx, ptr_int_lonlat%baryctr_blk, STAT=ist)
        IF (ist /= SUCCESS) THEN
          CALL finish (routine, 'deallocation for barycentric stencil indices failed')
        ENDIF
      END IF
      IF (ALLOCATED(ptr_int_lonlat%baryctr_coeff)) THEN
        DEALLOCATE (ptr_int_lonlat%baryctr_coeff, STAT=ist)
        IF (ist /= SUCCESS) THEN
          CALL finish (routine, 'deallocation for barycentric lon-lat coefficients failed')
        ENDIF
      END IF

      ! -- array with lon-lat coordinates
      DEALLOCATE(ptr_int_lonlat%ll_coord, stat=ist)
      IF (ist /= SUCCESS)  CALL finish (routine, 'Deallocation of array with lon-lat coordinates!')

    END SUBROUTINE deallocate_int_state_lonlat


    !>
    !! This routine creates the index list needed for RBF reconstruction
    !! of scalar values at cell centers to lon-lat points.
    !!
    SUBROUTINE rbf_c2l_index( ptr_patch, ptr_int, ptr_int_lonlat )

      ! patch on which computation is performed
      TYPE(t_patch), TARGET, INTENT(in)    :: ptr_patch
      TYPE (t_int_state),   INTENT(IN)     :: ptr_int
      TYPE(t_lon_lat_intp), INTENT(inout)  :: ptr_int_lonlat

      CHARACTER(*), PARAMETER :: routine = modname//"::rbf_c2l_index"
      INTEGER :: ilc_n(3), ibc_n(3)       ! line and block index for neighbors of
      ! direct neighbors
      INTEGER :: ilv(3), ibv(3)           ! vertex line and block indices
      INTEGER :: ilc_v(3,6), ibc_v(3,6)   ! cell line and block indices
      ! around each of the three vertices
      INTEGER :: jb                       ! loop index blocks
      INTEGER :: jc                       ! loop index cells
      INTEGER :: jj                       ! loop index
      INTEGER :: jec                      ! loop index
      INTEGER :: jtri                     ! loop index
      INTEGER :: cnt                      ! counter
      INTEGER :: i_startblk               ! start block
      INTEGER :: i_startidx               ! start index
      INTEGER :: i_endidx                 ! end index
      INTEGER :: rl_start, rl_end         ! refinement control start level
      INTEGER :: i_nchdom, i_endblk

      REAL(wp) :: z_stencil(UBOUND(ptr_int_lonlat%rbf_c2l_stencil,1), &
        &                   UBOUND(ptr_int_lonlat%rbf_c2l_stencil,2))

      !--------------------------------------------------------------------

      ! stencil size: 4  = nearest neighbor,
      !               13 = vertex stencil,
      !               rbf_c2grad_dim = edge stencil

      IF ((rbf_dim_c2l == 4) .OR. (rbf_dim_c2l == 13)) THEN

        rl_start = 2
        rl_end   = min_rlcell_int
        ptr_int_lonlat%rbf_c2l_idx(:,:,:) = 0

        ! The start block depends on the width of the stencil
        i_nchdom   = MAX(1,ptr_patch%n_childdom)
        i_startblk = ptr_patch%cells%start_blk(rl_start,1)
        i_endblk   = ptr_patch%cells%end_blk(rl_end,i_nchdom)

!$OMP PARALLEL

        ! The stencil consists of 12 cells surrounding the control
        ! volume that share a vertex with the control volume, plus the
        ! control volume itself.
        !
        ! Note: At pentagon points the size of the stencil reduces to 12.

!$OMP DO PRIVATE(jb,jc,jec,jj,jtri,i_startidx,i_endidx,cnt,ilv,ibv, &
!$OMP            ilc_v,ibc_v,ilc_n,ibc_n)
        DO jb = i_startblk, i_endblk

          CALL get_indices_c(ptr_patch, jb,        &
            &                i_startblk, i_endblk, &
            &                i_startidx, i_endidx, &
            &                rl_start, rl_end)

          DO jc = i_startidx, i_endidx

            IF(.NOT. ptr_patch%cells%decomp_info%owner_mask(jc,jb)) CYCLE

            cnt = 1

            ! get get line and block indices of cell vertices
            ilv(1:3) = ptr_patch%cells%vertex_idx(jc,jb,1:3)
            ibv(1:3) = ptr_patch%cells%vertex_blk(jc,jb,1:3)

            ! for each vertex: get all the cells which share this vertex
            DO jj = 1,3
              ilc_v(jj,:)=ptr_patch%verts%cell_idx(ilv(jj),ibv(jj),:)
              ibc_v(jj,:)=ptr_patch%verts%cell_blk(ilv(jj),ibv(jj),:)
            ENDDO

            ! 1. add the 3 direct neighbors to the stencil
            DO jec = 1, 3
              ! get line and block indices of direct neighbors
              ilc_n(jec) = ptr_patch%cells%neighbor_idx(jc,jb,jec)
              ibc_n(jec) = ptr_patch%cells%neighbor_blk(jc,jb,jec)

              ptr_int_lonlat%rbf_c2l_idx(cnt,jc,jb) = ilc_n(jec)
              ptr_int_lonlat%rbf_c2l_blk(cnt,jc,jb) = ibc_n(jec)

              cnt = cnt + 1
            ENDDO

            ! 2. loop over the vertices and add all the cells
            !    that are no direct neighbors and not our CV.

            IF (rbf_dim_c2l == 13) THEN
              DO jj = 1,3   ! loop over vertices
                DO jtri=1,6 ! loop over cells around each vertex

                  IF (.NOT.( (ilc_v(jj,jtri) == ilc_n(1) .AND. ibc_v(jj,jtri) == ibc_n(1))  &
                    &  .OR.  (ilc_v(jj,jtri) == ilc_n(2) .AND. ibc_v(jj,jtri) == ibc_n(2))  &
                    &  .OR.  (ilc_v(jj,jtri) == ilc_n(3) .AND. ibc_v(jj,jtri) == ibc_n(3))  &
                    &  .OR.  (ilc_v(jj,jtri) == jc       .AND. ibc_v(jj,jtri) == jb)        &
                    &  .OR.  (ilc_v(jj,jtri) == 0        .AND. ibc_v(jj,jtri) == 0 ) ) ) THEN

                    ptr_int_lonlat%rbf_c2l_idx(cnt,jc,jb) = ilc_v(jj,jtri)
                    ptr_int_lonlat%rbf_c2l_blk(cnt,jc,jb) = ibc_v(jj,jtri)

                    cnt = cnt + 1
                  ENDIF
                ENDDO
              ENDDO
            END IF

            ! 3. Finally, add the control volume itself
            ptr_int_lonlat%rbf_c2l_idx(cnt,jc,jb) = jc
            ptr_int_lonlat%rbf_c2l_blk(cnt,jc,jb) = jb

            ptr_int_lonlat%rbf_c2l_stencil(jc,jb) = cnt

          ENDDO ! loop over cells
        ENDDO ! loop over blocks
!$OMP END DO NOWAIT
!$OMP END PARALLEL

        DO cnt = 1, rbf_dim_c2l
          CALL sync_idx(SYNC_C, SYNC_C, ptr_patch, ptr_int_lonlat%rbf_c2l_idx(cnt,:,:), &
            &                                      ptr_int_lonlat%rbf_c2l_blk(cnt,:,:))
        ENDDO
        z_stencil(:,:) = REAL(ptr_int_lonlat%rbf_c2l_stencil(:,:),wp)
        CALL sync_patch_array(SYNC_C,ptr_patch,z_stencil)
        ptr_int_lonlat%rbf_c2l_stencil(:,:) = NINT(z_stencil(:,:))

      ELSE IF (rbf_dim_c2l == rbf_c2grad_dim) THEN

        ! copy stencil info from "ptr_int%rbf_c2grad_idx":
        ptr_int_lonlat%rbf_c2l_idx(:,:,:)   = ptr_int%rbf_c2grad_idx
        ptr_int_lonlat%rbf_c2l_blk(:,:,:)   = ptr_int%rbf_c2grad_blk
        ptr_int_lonlat%rbf_c2l_stencil(:,:) = rbf_dim_c2l

      ELSE

        CALL finish(routine, "Unknown stencil!")

      END IF

    END SUBROUTINE rbf_c2l_index


    !===============================================================
    ! COMPUTATION OF LON-LAT INTERPOLATION COEFFICIENTS

    !-------------------------------------------------------------------------
    !> Computes the RBF coefficients for lon-lat grid points.
    !-------------------------------------------------------------------------
    !
    ! This routine is based on mo_intp_rbf_coeffs::rbf_vec_compute_coeff_cell()
    !
    ! @par Revision History
    !      Initial implementation  by  F. Prill, DWD (2011-08)
    !
    SUBROUTINE rbf_vec_compute_coeff_lonlat( ptr_patch, ptr_int_lonlat,  &
      &                                      nblks_lonlat, npromz_lonlat, rbf_shape_param )

      ! Input parameters
      TYPE(t_patch),                 INTENT(IN)    :: ptr_patch
      TYPE (t_lon_lat_intp),         INTENT(INOUT) :: ptr_int_lonlat
      INTEGER,                       INTENT(IN)    :: nblks_lonlat, npromz_lonlat ! blocking info
      REAL(wp),                      INTENT(IN)    :: rbf_shape_param
      ! Local parameters
      CHARACTER(*), PARAMETER :: routine = modname//"::rbf_vec_compute_coeff_lonlat"
      REAL(wp)                         :: cc_e1(3), cc_e2(3), cc_c(nproma,3)  ! coordinates of edge midpoints
      TYPE(t_cartesian_coordinates)    :: cc_center                  ! cartes. coordinates of lon-lat points
      REAL(wp), DIMENSION (nproma,3)   :: z_nx1, z_nx2, z_nx3        ! 3d  normal velocity vectors at edge midpoints
      REAL(wp)                         :: z_lon, z_lat,            & ! longitude and latitude
        &                                 z_norm,                  & ! norm of velocity vectors
        &                                 z_dist,                  & ! distance between data points
        &                                 z_nxprod                   ! scalar prod of normal veloc.
      REAL(wp),ALLOCATABLE             :: z_rbfmat(:,:,:),         & ! RBF interpolation matrix
        &                                 z_diag(:,:),             & ! diagonal of Cholesky decomp.
        &                                 z_rbfval(:,:),           & ! RBF function value
        &                                 z_rhs1(:,:),             & ! right hand side of linear
        &                                 z_rhs2(:,:)                ! interpolation system in 2d
      INTEGER                          :: &
        &                                 jc, jb,                  & ! integer over lon-lat grid points
        &                                 je1, je2,                & ! integer over edges
        &                                 ile1, ibe1, ile2, ibe2,  & ! edge indices
        &                                 ist,                     & ! return value of array allocation
        &                                 i_startidx, i_endidx,    & ! start/end index
        &                                 istencil(nproma)           ! actual number of stencil points
      REAL(wp)                         :: checksum_u,checksum_v      ! to check if sum of interpolation coefficients is correct
      TYPE(t_geographical_coordinates) :: grid_point

      !--------------------------------------------------------------------

      CALL message(routine, '')
      IF (ptr_patch%n_patch_cells == 0) RETURN;

!$OMP PARALLEL PRIVATE (z_rbfmat,z_diag,z_rbfval,z_rhs1,z_rhs2, ist,   &
!$OMP                   grid_point),                                   &
!$OMP          SHARED  (nproma, rbf_vec_dim_c, nblks_lonlat,           &
!$OMP                   npromz_lonlat, ptr_int_lonlat, ptr_patch,      &
!$OMP                   rbf_vec_kern_ll)

      ALLOCATE( z_rbfmat(nproma,rbf_vec_dim_c,rbf_vec_dim_c),  &
        z_diag(nproma,rbf_vec_dim_c),                  &
        z_rbfval(nproma,rbf_vec_dim_c),                &
        z_rhs1(nproma,rbf_vec_dim_c),                  &
        z_rhs2(nproma,rbf_vec_dim_c), STAT=ist )
      IF (ist /= SUCCESS) THEN
        CALL finish (routine, 'allocation for working arrays failed')
      ENDIF

!$OMP DO PRIVATE (jb,jc,i_startidx,i_endidx,je1,je2,istencil,       &
!$OMP             ist,ile1,ibe1,cc_e1,z_lon,z_lat,z_norm,           &
!$OMP             z_nx1,ile2,ibe2,cc_e2,cc_c,z_nx2,z_nxprod,z_dist, &
!$OMP             cc_center,z_nx3,checksum_u,checksum_v)
      BLOCKS: DO jb = 1, nblks_lonlat

        i_startidx = 1
        i_endidx   = nproma
        if (jb == nblks_lonlat) i_endidx = npromz_lonlat

        ! for each cell, build the vector RBF interpolation matrix
        DO je1 = 1, rbf_vec_dim_c
          DO je2 = 1, je1
            DO jc = i_startidx, i_endidx

              ! Get actual number of stencil points
              istencil(jc) = ptr_int_lonlat%rbf_vec_stencil(jc,jb)

              ! paranoia:
              IF ( (je1 > istencil(jc)) .OR. (je2 > istencil(jc)) ) CYCLE

              ! line and block indices for each edge je1 and je2 of RBF stencil
              ile1 = ptr_int_lonlat%rbf_vec_idx(je1,jc,jb)
              ibe1 = ptr_int_lonlat%rbf_vec_blk(je1,jc,jb)
              ile2 = ptr_int_lonlat%rbf_vec_idx(je2,jc,jb)
              ibe2 = ptr_int_lonlat%rbf_vec_blk(je2,jc,jb)

              ! get Cartesian coordinates and orientation vectors
              cc_e1(:) = ptr_patch%edges%cartesian_center(ile1,ibe1)%x(:)
              cc_e2(:) = ptr_patch%edges%cartesian_center(ile2,ibe2)%x(:)
              !
              z_nx1(jc,:) = ptr_patch%edges%primal_cart_normal(ile1,ibe1)%x(:)
              z_nx2(jc,:) = ptr_patch%edges%primal_cart_normal(ile2,ibe2)%x(:)

              ! compute dot product of normal vectors and distance between edge midpoints
              z_nxprod = DOT_PRODUCT(z_nx1(jc,:),z_nx2(jc,:))
              z_dist   = arc_length_v(cc_e1,cc_e2)

              ! set up interpolation matrix
              IF      (rbf_vec_kern_ll == 1) THEN
                z_rbfmat(jc,je1,je2) = z_nxprod * gaussi(z_dist,rbf_shape_param)
              ELSE IF (rbf_vec_kern_ll == 3) THEN
                z_rbfmat(jc,je1,je2) = z_nxprod * inv_multiq(z_dist,rbf_shape_param)
              ENDIF

              IF (je1 > je2) z_rbfmat(jc,je2,je1) = z_rbfmat(jc,je1,je2)

            END DO
          END DO
        END DO

        ! apply Cholesky decomposition to matrix
        !
!CDIR NOIEXPAND
#ifdef __SX__
        CALL choldec_v(i_startidx,i_endidx,istencil,rbf_vec_dim_c,z_rbfmat,z_diag)
#else
        CALL choldec_v(i_startidx,i_endidx,istencil,              z_rbfmat,z_diag)
#endif

        DO jc = i_startidx, i_endidx

          !
          ! Solve immediately for coefficients
          !
          ! convert coordinates to cartesian vector
          !
          grid_point = ptr_int_lonlat%ll_coord(jc,jb)
          cc_center = gc2cc(grid_point)
          cc_c(jc,1:3) = cc_center%x(1:3)

          z_lon = grid_point%lon
          z_lat = grid_point%lat

          ! Zonal wind component
          CALL gvec2cvec(1._wp,0._wp, z_lon,z_lat, &
            &            z_nx1(jc,1),z_nx1(jc,2),z_nx1(jc,3))

          z_norm = SQRT( DOT_PRODUCT(z_nx1(jc,:),z_nx1(jc,:)) )
          z_nx1(jc,:)  = 1._wp/z_norm * z_nx1(jc,:)

          ! Meridional wind component
          CALL gvec2cvec(0._wp,1._wp, z_lon,z_lat, &
            &            z_nx2(jc,1),z_nx2(jc,2),z_nx2(jc,3))

          z_norm = SQRT( DOT_PRODUCT(z_nx2(jc,:),z_nx2(jc,:)) )
          z_nx2(jc,:)  = 1._wp/z_norm * z_nx2(jc,:)

        END DO

        !
        ! set up right hand side for interpolation system
        !
        DO je2 = 1, rbf_vec_dim_c
          DO jc = i_startidx, i_endidx

            IF (je2 > istencil(jc)) CYCLE

            ! get indices and coordinates of edge midpoints and compute distance
            ! to cell center
            ile2   = ptr_int_lonlat%rbf_vec_idx(je2,jc,jb)
            ibe2   = ptr_int_lonlat%rbf_vec_blk(je2,jc,jb)

            cc_e2(:)  = ptr_patch%edges%cartesian_center(ile2,ibe2)%x(:)

            z_dist = arc_length_v(cc_c(jc,:), cc_e2)

            ! get Cartesian orientation vector
            z_nx3(jc,:) = ptr_patch%edges%primal_cart_normal(ile2,ibe2)%x(:)

            IF (rbf_vec_kern_ll == 1) THEN
              z_rbfval(jc,je2) = gaussi(z_dist,rbf_shape_param)
            ELSE IF (rbf_vec_kern_ll == 3) THEN
              z_rbfval(jc,je2) = inv_multiq(z_dist,rbf_shape_param)
            ENDIF

            ! compute projection on target vector orientation
            z_rhs1(jc,je2) = z_rbfval(jc,je2) * &
              &                     DOT_PRODUCT(z_nx1(jc,:),z_nx3(jc,:))
            z_rhs2(jc,je2) = z_rbfval(jc,je2) * &
              &                     DOT_PRODUCT(z_nx2(jc,:),z_nx3(jc,:))

          END DO
        END DO

        ! compute vector coefficients
!CDIR NOIEXPAND
#ifdef __SX__
        CALL solve_chol_v(i_startidx, i_endidx, istencil, rbf_vec_dim_c, z_rbfmat,  &
          &               z_diag, z_rhs1, ptr_int_lonlat%rbf_vec_coeff(:,1,:,jb))
#else
        CALL solve_chol_v(i_startidx, i_endidx, istencil,                z_rbfmat,  &
          &               z_diag, z_rhs1, ptr_int_lonlat%rbf_vec_coeff(:,1,:,jb))
#endif
!CDIR NOIEXPAND
#ifdef __SX__
        CALL solve_chol_v(i_startidx, i_endidx, istencil, rbf_vec_dim_c, z_rbfmat,  &
          &               z_diag, z_rhs2, ptr_int_lonlat%rbf_vec_coeff(:,2,:,jb))
#else
        CALL solve_chol_v(i_startidx, i_endidx, istencil,                z_rbfmat,  &
          &               z_diag, z_rhs2, ptr_int_lonlat%rbf_vec_coeff(:,2,:,jb))
#endif

        DO jc = i_startidx, i_endidx

          ! Ensure that sum of interpolation coefficients is correct
          checksum_u = 0._wp
          checksum_v = 0._wp

          DO je1 = 1, istencil(jc)
            ile1   = ptr_int_lonlat%rbf_vec_idx(je1,jc,jb)
            ibe1   = ptr_int_lonlat%rbf_vec_blk(je1,jc,jb)

            ! get Cartesian orientation vector
            z_nx3(jc,:) = ptr_patch%edges%primal_cart_normal(ile1,ibe1)%x(:)

            checksum_u = checksum_u + ptr_int_lonlat%rbf_vec_coeff(je1,1,jc,jb)* &
              DOT_PRODUCT(z_nx1(jc,:),z_nx3(jc,:))
            checksum_v = checksum_v + ptr_int_lonlat%rbf_vec_coeff(je1,2,jc,jb)* &
              DOT_PRODUCT(z_nx2(jc,:),z_nx3(jc,:))
          ENDDO

          DO je1 = 1, istencil(jc)
            ptr_int_lonlat%rbf_vec_coeff(je1,1,jc,jb) = &
              ptr_int_lonlat%rbf_vec_coeff(je1,1,jc,jb) / checksum_u
          END DO
          DO je1 = 1, istencil(jc)
            ptr_int_lonlat%rbf_vec_coeff(je1,2,jc,jb) = &
              ptr_int_lonlat%rbf_vec_coeff(je1,2,jc,jb) / checksum_v
          END DO

        END DO ! jc

      END DO BLOCKS
!$OMP END DO

      DEALLOCATE( z_rbfmat, z_diag, z_rbfval, z_rhs1, z_rhs2,  STAT=ist )
      IF (ist /= SUCCESS) THEN
        CALL finish (routine, 'deallocation for working arrays failed')
      ENDIF
!$OMP END PARALLEL

    END SUBROUTINE rbf_vec_compute_coeff_lonlat


    !-------------------------------------------------------------------------
    !> Computes the RBF coefficients for lon-lat grid points.
    !
    !  We directly interpolate from cell centers to lon-lat points,
    !  otherwise we do gradient interpolation and reconstruction.
    !  -------------------------------------------------------------------------
    !
    ! @par Revision History
    !      Initial implementation  by  F. Prill, DWD (2012-06-13)
    !
    SUBROUTINE rbf_compute_coeff_c2l( ptr_patch, ptr_int_lonlat, &
      &                               nblks_lonlat, npromz_lonlat, rbf_shape_param )

      ! Input parameters
      TYPE(t_patch),                 INTENT(IN)    :: ptr_patch
      TYPE (t_lon_lat_intp),         INTENT(INOUT) :: ptr_int_lonlat
      INTEGER,                       INTENT(IN)    :: nblks_lonlat, npromz_lonlat ! blocking info
      REAL(wp),                      INTENT(IN)    :: rbf_shape_param
      ! Local parameters
      CHARACTER(*), PARAMETER :: routine = modname//"::rbf_compute_coeff_c2l"
      REAL(wp)                         :: cc_c(nproma,3)             ! coordinates of cell centers
      TYPE(t_cartesian_coordinates)    :: cc_center, cc_1, cc_2      ! temporary variables
      REAL(wp)                         :: z_dist                     ! distance between data points
      REAL(wp),ALLOCATABLE             :: z_rbfmat(:,:,:),         & ! RBF interpolation matrix
        &                                 z_diag(:,:),             & ! diagonal of Cholesky decomp.
        &                                 z_rbfval(:,:)              ! RBF function value
      INTEGER                          :: jc, jb,                  & ! integer over lon-lat grid points
        &                                 je1, je2,                & ! integer over edges
        &                                 ilc1, ibc1, ilc2, ibc2,  & ! cell indices
        &                                 ist,                     & ! return value of array allocation
        &                                 i_startidx, i_endidx,    & ! start/end index
        &                                 istencil(nproma),        & ! actual number of stencil points
        &                                 jb_cell, jc_cell
      REAL(wp)                         :: checksum                   ! to check if sum of interpolation coefficients is correct
      TYPE(t_geographical_coordinates) :: grid_point

      !--------------------------------------------------------------------

      CALL message(routine, '')
      IF (ptr_patch%n_patch_cells == 0) RETURN;

!$OMP PARALLEL PRIVATE (z_rbfmat,z_diag,z_rbfval, ist)

      ALLOCATE( z_rbfmat(nproma,rbf_dim_c2l,rbf_dim_c2l),    &
        &       z_diag(nproma,rbf_dim_c2l),                  &
        &       z_rbfval(nproma,rbf_dim_c2l),                &
        &       STAT=ist )
      IF (ist /= SUCCESS) &
        & CALL finish (routine, 'allocation for working arrays failed')

!$OMP DO PRIVATE (jb,jc,i_startidx,i_endidx,je1,je2,istencil, grid_point, &
!$OMP             ist,ilc1,ibc1, ilc2,ibc2,cc_1, cc_2,cc_c,z_dist,        &
!$OMP             cc_center, checksum, jb_cell, jc_cell )
      BLOCKS: DO jb = 1, nblks_lonlat

        i_startidx = 1
        i_endidx   = nproma
        if (jb == nblks_lonlat) i_endidx = npromz_lonlat

        ! for each cell, build the vector RBF interpolation matrix
        z_rbfmat(:,:,:) = 0._wp
        DO je1 = 1, rbf_dim_c2l
          DO je2 = 1, je1
            DO jc = i_startidx, i_endidx

              jc_cell = ptr_int_lonlat%tri_idx(1,jc, jb)
              jb_cell = ptr_int_lonlat%tri_idx(2,jc, jb)
              ! Get actual number of stencil points
              istencil(jc) = ptr_int_lonlat%rbf_c2l_stencil(jc_cell,jb_cell)
              ! paranoia:
              IF ( (je1 > istencil(jc)) .OR. (je2 > istencil(jc)) ) CYCLE

              ! line and block indices for each cell of RBF stencil
              ilc1 = ptr_int_lonlat%rbf_c2l_idx(je1,jc_cell,jb_cell)
              ibc1 = ptr_int_lonlat%rbf_c2l_blk(je1,jc_cell,jb_cell)
              ilc2 = ptr_int_lonlat%rbf_c2l_idx(je2,jc_cell,jb_cell)
              ibc2 = ptr_int_lonlat%rbf_c2l_blk(je2,jc_cell,jb_cell)

              ! get Cartesian coordinates and orientation vectors
              cc_1 = gc2cc(ptr_patch%cells%center(ilc1,ibc1))
              cc_2 = gc2cc(ptr_patch%cells%center(ilc2,ibc2))
              ! compute distance between cell centers
              z_dist = arc_length_v(cc_1%x(:),cc_2%x(:))

              ! set up interpolation matrix
              IF      (rbf_vec_kern_ll == 1) THEN
                z_rbfmat(jc,je1,je2) = gaussi(z_dist,rbf_shape_param)
              ELSE IF (rbf_vec_kern_ll == 3) THEN
                z_rbfmat(jc,je1,je2) = inv_multiq(z_dist,rbf_shape_param)
              ENDIF

              IF (je1 > je2) z_rbfmat(jc,je2,je1) = z_rbfmat(jc,je1,je2)

            END DO
          END DO
        END DO

        ! apply Cholesky decomposition to matrix
        !
!CDIR NOIEXPAND
#ifdef __SX__
        CALL choldec_v(i_startidx,i_endidx,istencil,rbf_dim_c2l,z_rbfmat,z_diag)
#else
        CALL choldec_v(i_startidx,i_endidx,istencil,            z_rbfmat,z_diag)
#endif

        ! compute RHS for coefficient computation
        DO jc = i_startidx, i_endidx

          grid_point = ptr_int_lonlat%ll_coord(jc,jb)
          ! convert coordinates to cartesian vector
          cc_center = gc2cc(grid_point)
          cc_c(jc,1:3) = cc_center%x(1:3)

        END DO

        !
        ! set up right hand side for interpolation system
        !
        DO je2 = 1, rbf_dim_c2l
          DO jc = i_startidx, i_endidx

            IF (je2 > istencil(jc)) CYCLE

            jc_cell = ptr_int_lonlat%tri_idx(1,jc, jb)
            jb_cell = ptr_int_lonlat%tri_idx(2,jc, jb)

            ! get indices and coordinates of cell centers and compute
            ! distance to lon-lat point:
            ilc2   = ptr_int_lonlat%rbf_c2l_idx(je2,jc_cell,jb_cell)
            ibc2   = ptr_int_lonlat%rbf_c2l_blk(je2,jc_cell,jb_cell)
            cc_1 = gc2cc(ptr_patch%cells%center(ilc2,ibc2))

            z_dist = arc_length_v(cc_c(jc,:), cc_1%x(:))

            IF (rbf_vec_kern_ll == 1) THEN
              z_rbfval(jc,je2) = gaussi(z_dist,rbf_shape_param)
            ELSE IF (rbf_vec_kern_ll == 3) THEN
              z_rbfval(jc,je2) = inv_multiq(z_dist,rbf_shape_param)
            ENDIF

          END DO
        END DO

        ! compute vector coefficients
!CDIR NOIEXPAND
#ifdef __SX__
        CALL solve_chol_v(i_startidx, i_endidx, istencil, rbf_dim_c2l, z_rbfmat,  &
          &               z_diag, z_rbfval, ptr_int_lonlat%rbf_c2l_coeff(:,:,jb))
#else
        CALL solve_chol_v(i_startidx, i_endidx, istencil,              z_rbfmat,  &
          &               z_diag, z_rbfval, ptr_int_lonlat%rbf_c2l_coeff(:,:,jb))
#endif

        DO jc = i_startidx, i_endidx

          ! Ensure that sum of interpolation coefficients is correct
          checksum = 0._wp
          DO je1 = 1, istencil(jc)
            checksum = checksum + ptr_int_lonlat%rbf_c2l_coeff(je1,jc,jb)
          ENDDO
          DO je1 = 1, istencil(jc)
            ptr_int_lonlat%rbf_c2l_coeff(je1,jc,jb) = &
              ptr_int_lonlat%rbf_c2l_coeff(je1,jc,jb) / checksum
          END DO
        END DO ! jc

      END DO BLOCKS
!$OMP END DO

      DEALLOCATE( z_rbfmat, z_diag, z_rbfval, STAT=ist )
      IF (ist /= SUCCESS) &
        CALL finish (routine, 'deallocation for working arrays failed')
!$OMP END PARALLEL

    END SUBROUTINE rbf_compute_coeff_c2l


    !-------------------------------------------------------------------------
    !> Mask out nest boundary zone.
    !  -------------------------------------------------------------------------
    !
    SUBROUTINE mask_out_boundary( ptr_patch, ptr_int_lonlat )

      ! Input parameters
      TYPE(t_patch),                 INTENT(IN)    :: ptr_patch
      TYPE (t_lon_lat_intp), TARGET, INTENT(INOUT) :: ptr_int_lonlat
      ! Local parameters
      CHARACTER(*), PARAMETER :: routine      = modname//"::mask_out_boundary"
      INTEGER,      PARAMETER :: intp_types(3) = (/ HINTP_TYPE_LONLAT_RBF, &
        &                                           HINTP_TYPE_LONLAT_NNB, &
        &                                           HINTP_TYPE_LONLAT_BCTR /)

      INTEGER  :: last_bdry_index, rl_start, rl_end, i_nchdom, i_startblk, i_endblk, &
        &         i_startidx, i_endidx, local_idx, glb_idx, jc, jb, je1,             &
        &         istencil, istencil2, iidx, iblk, itype, hintp_type, nblks_lonlat,  &
        &         npromz_lonlat
      REAL(wp) :: area
      REAL(wp), DIMENSION(:,:,:), POINTER :: ptr_coeff

      !--------------------------------------------------------------------

      IF (dbg_level > 1)  CALL message(routine,'')
      IF (ptr_patch%n_patch_cells == 0) RETURN;

      ! set local values for "nblks" and "npromz"
      nblks_lonlat  = (ptr_int_lonlat%nthis_local_pts - 1)/nproma + 1
      npromz_lonlat = ptr_int_lonlat%nthis_local_pts - (nblks_lonlat-1)*nproma

      ! communicate the largest global index of the lateral boundary cells.
      ! 
      ! note that this implicitly assumes that the indices of the nest
      ! boundary region are stored at the beginning of the solution
      ! vector in a contiguous fashion.
      last_bdry_index = 0

      rl_start   = 1
      rl_end     = grf_bdywidth_c
      i_nchdom   = MAX(1,ptr_patch%n_childdom)
      i_startblk = ptr_patch%cells%start_blk(rl_start,1)
      i_endblk   = ptr_patch%cells%end_blk(rl_end,i_nchdom)
      glb_idx    = 0
      DO jb=i_startblk,i_endblk
        CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
          i_startidx, i_endidx, rl_start, rl_end)
        DO jc=i_startidx, i_endidx
          local_idx = idx_1d(jc,jb)
          glb_idx   = MAX(glb_idx, ptr_patch%cells%decomp_info%glb_index(local_idx))
        END DO
      END DO
      last_bdry_index = p_max(glb_idx, p_comm_work)

      DO itype=1,SIZE(intp_types)
        hintp_type = intp_types(itype)
        IF (dbg_level > 1)  CALL message(routine, "interpolation type: "//int2string(hintp_type))

        ! depending on the interpolation type: set a pointer to the
        ! coefficient vector
        NULLIFY(ptr_coeff)
        SELECT CASE(hintp_type)
        CASE (HINTP_TYPE_LONLAT_RBF)
          ptr_coeff => ptr_int_lonlat%rbf_c2l_coeff
        CASE (HINTP_TYPE_LONLAT_NNB)
          ptr_coeff => ptr_int_lonlat%nnb_coeff
        CASE (HINTP_TYPE_LONLAT_BCTR)
          IF (.NOT. support_baryctr_intp)  CYCLE
          ptr_coeff => ptr_int_lonlat%baryctr_coeff
        CASE DEFAULT
          CALL finish(routine, "Internal error!")
        END SELECT

!$OMP PARALLEL DO PRIVATE(jb, i_startidx, i_endidx, jc, istencil, &
!$OMP                     istencil2, area, je1, iidx, iblk, local_idx, glb_idx)
        BLOCKS: DO jb = 1, nblks_lonlat

          i_startidx = 1
          i_endidx   = nproma
          IF (jb == nblks_lonlat) i_endidx = npromz_lonlat

          ! loop over all lonlat points and their input stencils
          DO jc = i_startidx, i_endidx

            ! Get actual number of stencil points
            istencil = 0
            SELECT CASE(hintp_type)
            CASE (HINTP_TYPE_LONLAT_RBF)
              istencil  = ptr_int_lonlat%rbf_c2lr_stencil(jc,jb)
            CASE (HINTP_TYPE_LONLAT_NNB)
              istencil  = 1
            CASE (HINTP_TYPE_LONLAT_BCTR)
              istencil  = 3
            END SELECT
            istencil2 = istencil
            ! first, multiply interpolation weights by total "area" (the
            ! sum of weights in the stencil). we perform this step
            ! though we may assume normalized coefficients (area=1).
            area = 0._wp
            DO je1 = 1, istencil
              area = area + ptr_coeff(je1,jc,jb)
            ENDDO
            DO je1 = 1, istencil
              ptr_coeff(je1,jc,jb) = ptr_coeff(je1,jc,jb) * area
            END DO
            ! remove stencil points from the nest boundary region
            ! together with their coefficients
            DO je1 = 1, istencil
              SELECT CASE(hintp_type)
              CASE (HINTP_TYPE_LONLAT_RBF)
                iidx = ptr_int_lonlat%rbf_c2lr_idx(je1,jc,jb)
                iblk = ptr_int_lonlat%rbf_c2lr_blk(je1,jc,jb)
              CASE (HINTP_TYPE_LONLAT_NNB)
                iidx = ptr_int_lonlat%tri_idx(1,jc,jb)
                iblk = ptr_int_lonlat%tri_idx(2,jc,jb)
              CASE (HINTP_TYPE_LONLAT_BCTR)
                iidx = ptr_int_lonlat%baryctr_idx(je1,jc,jb)
                iblk = ptr_int_lonlat%baryctr_blk(je1,jc,jb)
              END SELECT

              local_idx = idx_1d(iidx,iblk)
              glb_idx   = ptr_patch%cells%decomp_info%glb_index(local_idx)
              IF (glb_idx <= last_bdry_index) THEN
                area = area - ptr_coeff(je1,jc,jb)
                ptr_coeff(je1,jc,jb) = 0._wp
                istencil2 = istencil2 - 1
                IF (istencil2 == 0)  area = 1._wp ! avoid zero divide
              END IF
            END DO
            ! re-normalize weights
            DO je1 = 1, istencil
              ptr_coeff(je1,jc,jb) = ptr_coeff(je1,jc,jb) / area
            END DO
          END DO ! jc

        END DO BLOCKS
!$OMP END PARALLEL DO
      END DO
      IF (dbg_level > 1)  CALL message(routine, "done.")

    END SUBROUTINE mask_out_boundary


    ! -----------------------------------------------------------------------------------
    !> Utility function: distance computation
    !
    PURE FUNCTION dist_p(p1, p2)
      REAL(gk)              :: dist_p
      REAL(gk), INTENT(IN)  :: p1(2), p2(2)
      ! local variables
      REAL(gk) :: val
      
      ! spherical distance:
      val = SIN(p1(2))*SIN(p2(2)) + COS(p1(2))*COS(p2(2))*COS(p1(1)-p2(1))
      dist_p = ACOS( MIN(1._gk, MAX(-1._gk, val)) )
    END FUNCTION dist_p


    ! -----------------------------------------------------------------------------------
    !> Utility function: Divide-and-Conquer algorithm to determine
    !  points in lon-lat grid which are "far off" the current patch of
    !  the triangular grid.
    !
    RECURSIVE SUBROUTINE flag_ll_points(rotated_pts, s_lon, e_lon, s_lat, e_lat, &
      &                                 pts_flags, recursion_depth, max_recursion, gnat, min_radius)
      REAL(wp), INTENT(IN)    :: rotated_pts(:,:,:)  ! dim (lon,lat,1:2)
      INTEGER,  INTENT(IN)    :: s_lon, e_lon, s_lat, e_lat
      INTEGER,  INTENT(INOUT) :: pts_flags(:,:)
      INTEGER,  INTENT(IN)    :: recursion_depth, max_recursion
      TYPE(t_gnat_tree), INTENT(IN) :: gnat
      REAL(wp), INTENT(IN)    :: min_radius
      ! local variables
      INTEGER  :: m_lon, m_lat, c_lon(4), c_lat(4), v_lon(4,4), v_lat(4,4), &
        &         s_lon2(4), s_lat2(4), e_lon2(4), e_lat2(4), icirc, ivert, &
        &         i1, i2
      REAL(gk) :: p1(4,2), p2(2), radius(4)
      LOGICAL  :: is_pt_inside(4)

      IF (recursion_depth == max_recursion) RETURN
      if ((s_lon == e_lon) .or. (s_lat == e_lat)) return

      ! determine mid longitude and latitude in current range, thus
      ! dividing the regular grid into 4 subgrids:
      m_lon = (e_lon + s_lon)/2
      m_lat = (e_lat + s_lat)/2

      ! nomenclature
      !
      ! *-------*-------*   - e_lat
      ! |       |       |
      ! |   3   |   4   |
      ! |       |       |
      ! *-------*-------*   - m_lat
      ! |       |       |
      ! |   1   |   2   |
      ! |       |       |
      ! *-------*-------*   - s_lat
      !
      ! |       |       |
      !  s_lon   m_lon   e_lon

      s_lon2(1:4) = (/ s_lon, m_lon, s_lon, m_lon /)
      e_lon2(1:4) = (/ m_lon, e_lon, m_lon, e_lon /)
      s_lat2(1:4) = (/ s_lat, s_lat, m_lat, m_lat /)
      e_lat2(1:4) = (/ m_lat, m_lat, e_lat, e_lat /)
      ! for each of the four subblocks, determine centre point and
      ! vertex indices:
      DO icirc=1,4
        c_lon(icirc) = (s_lon2(icirc) + e_lon2(icirc)) / 2
        c_lat(icirc) = (s_lat2(icirc) + e_lat2(icirc)) / 2
        ! vertex numbering counter-clockwise
        v_lon(icirc,1:4) = (/ s_lon2(icirc), e_lon2(icirc), e_lon2(icirc), s_lon2(icirc) /)
        v_lat(icirc,1:4) = (/ s_lat2(icirc), s_lat2(icirc), e_lat2(icirc), e_lat2(icirc) /)
      END DO
      ! determine radii of circles overlapping each subblock:
      radius(:) = 0._gk
      DO icirc=1,4
        p1(icirc,1:2) = rotated_pts(c_lon(icirc),c_lat(icirc),1:2)
        DO ivert=1,4
          p2(1:2) = rotated_pts(v_lon(icirc,ivert),v_lat(icirc,ivert),1:2)
          radius(icirc) = MAX(radius(icirc), dist_p(p1(icirc,1:2), p2))
        END DO
      END DO
      radius(1:4) = MAX(min_radius, 1.25_gk * radius(1:4)) ! safety margin
      ! test, if there exists a triangle center inside the circles:
      DO icirc=1,4
        is_pt_inside(icirc) = .FALSE.
        CALL gnat_recursive_proximity_query(gnat, gnat%gnat_tree, p1(icirc,1:2), radius(icirc), is_pt_inside(icirc))
      END DO
      ! for each of the four subblocks:
      DO icirc=1,4
        IF (.NOT. is_pt_inside(icirc)) THEN
          ! if patch of triangular grid does not intersect circle:
          ! mark corresponding points in regular grid with "SKIP_NODE"
          DO i2=s_lat2(icirc),e_lat2(icirc)
            DO i1=s_lon2(icirc),e_lon2(icirc)
              pts_flags(i1,i2) = SKIP_NODE
            END DO
          END DO
        ELSE
          ! if patch of triangular grid does intersect circle:
          ! recursion and further subdivision:
          CALL flag_ll_points(rotated_pts, s_lon2(icirc), e_lon2(icirc), s_lat2(icirc), e_lat2(icirc), &
            &                 pts_flags, (recursion_depth+1), max_recursion, gnat, min_radius)          
        END IF
      END DO
    END SUBROUTINE flag_ll_points


    !-------------------------------------------------------------------------
    !> Setup routine for RBF reconstruction at lon-lat grid points for
    !  an arbitrary grid.
    !
    ! @par Revision History
    !      Initial implementation  by  F. Prill, DWD (2011-08)
    !      Changed for abritrary grids by Rainer Johanni (2011-11)
    !
    ! @note   This subroutine assumes that the GNAT data structure for
    !         this patch has already been initialized.
    !
    ! Note:   Two variables are created globally: "tri_idx" and "global_idx"
    !         This is necessary for the splitting of the lon-lat grid
    !         points over the PEs and for creating the corresponding
    !         communication patterns. However, after this initial phase,
    !         these data structures are resized to the local sizes.
    !
    SUBROUTINE rbf_setup_interpol_lonlat_grid(grid, gnat, ptr_patch, ptr_int_lonlat, ptr_int)

      TYPE (t_lon_lat_grid), INTENT(INOUT)         :: grid
      TYPE (t_gnat_tree),    INTENT(INOUT)         :: gnat
      ! data structure containing grid info:
      TYPE(t_patch), TARGET, INTENT(IN)            :: ptr_patch
      ! Indices of source points and interpolation coefficients
      TYPE (t_lon_lat_intp), TARGET, INTENT(INOUT) :: ptr_int_lonlat
      TYPE (t_int_state),    TARGET, INTENT(INOUT) :: ptr_int
      ! Local Parameters:
      CHARACTER(*), PARAMETER :: routine = modname//"::rbf_setup_interpol_lonlat"
      ! Flag: .TRUE., if we want to erase values outside local domains
      LOGICAL,      PARAMETER :: l_cutoff_local_domains = .TRUE.
      ! Flag: .TRUE., if we want raw performance measurements
      LOGICAL,      PARAMETER :: l_measure_time         = .FALSE.
      ! Zero threshold constant for equality test
      REAL(wp),     PARAMETER :: ZERO_TOL               = 1.e-15_wp

      REAL(wp), ALLOCATABLE            :: rotated_pts(:,:,:)
      REAL(gk), ALLOCATABLE            :: in_points(:,:,:)
      INTEGER,  ALLOCATABLE            :: pts_flags(:,:)
      REAL(gk), ALLOCATABLE            :: min_dist(:,:)       ! minimal distance
      INTEGER                          :: jb, jc, i_startidx, i_endidx,           &
        &                                 nblks_lonlat, npromz_lonlat, errstat,   &
        &                                 jb_lonlat, jc_lonlat,                   &
        &                                 rl_start, rl_end, i_startblk,           &
        &                                 i, j, row_idx(2), ibeg_idx, ibeg_blk,   &
        &                                 iend_idx, iend_blk, tri_idx_pole,       &
        &                                 ibeg_glb
      REAL(gk)                         :: min_dist_pole
      REAL(wp)                         :: point(2),                               &
        &                                 max_dist, start_radius, rbf_shape_param
!$    REAL                             :: time1
      LOGICAL                          :: l_grid_is_unrotated, l_grid_contains_poles


      !-----------------------------------------------------------------------

      IF (dbg_level > 1) THEN
        WRITE(message_text,*) "SETUP : interpol_lonlat"
        CALL message(routine, message_text)
      END IF
      IF (ptr_patch%geometry_info%cell_type == 6) &
        &   CALL finish(routine, "Lon-lat interpolation not yet implemented for cell_type == 6!")

      nblks_lonlat  = grid%nblks
      npromz_lonlat = grid%npromz

      ! ------------------------------------------
      ! check, if we have a global, unrotated grid
      ! ------------------------------------------

      l_grid_is_unrotated   = ((ABS(90._wp - grid%north_pole(2)) < ZERO_TOL)       .AND.  &
        &                       ABS( 0._wp - grid%north_pole(1)) < ZERO_TOL)
      l_grid_contains_poles = ((ABS(90._wp - ABS(grid%reg_lat_def(1))) < ZERO_TOL) .AND.  &
        &                      (ABS(90._wp - ABS(grid%reg_lat_def(3))) < ZERO_TOL))

      ! ---------------------------------------------------------
      ! preliminaries: make a sensible guess for maximum distance
      ! ---------------------------------------------------------

      ! make a sensible guess for maximum distance (just taking a
      ! "randomly chosen" edge length)
      rl_start   = 2
      rl_end     = min_rlcell_int
      i_startblk = ptr_patch%cells%start_blk(rl_start,1)
      CALL get_indices_e(ptr_patch, i_startblk, i_startblk, i_startblk, &
        &                i_startidx, i_endidx, rl_start, rl_end)
      max_dist = 3._gk * &
        & REAL(ptr_patch%edges%primal_edge_length(i_startidx,i_startblk)/grid_sphere_radius, gk)
      ! for MPI-independent behaviour: determine global max.
      max_dist = p_max(max_dist, comm=p_comm_work)
      IF(p_test_run) THEN
        IF(.NOT. my_process_is_mpi_test()) THEN
          ! Send to test PE
          CALL p_send(max_dist, process_mpi_all_test_id, 1)
        ELSE
          ! Receive result from parallel worker PEs
          CALL p_recv(max_dist, process_mpi_all_workroot_id, 1)
        END IF
      END IF

      ! -------------------------------------------
      ! compute grid points of rotated lon/lat grid
      ! -------------------------------------------

      ALLOCATE(in_points(nproma, nblks_lonlat, 2),            &
        &      rotated_pts(grid%lon_dim, grid%lat_dim, 2),    &
        &      pts_flags(grid%lon_dim, grid%lat_dim),         &
        &      min_dist(nproma, nblks_lonlat),                &
        &      stat=errstat)
      IF (errstat /= SUCCESS) CALL finish (routine, 'allocation for working arrays failed')

      IF (dbg_level > 1) CALL message(routine, "rotate lon-lat grid points")
      ! compute grid points (in radians): "rotated_pts(lon,lat,1:2)"
      ! contains the longitude and latitude of the point with grid
      ! index lon/lat.
      CALL rotate_latlon_grid( grid, rotated_pts )
      in_points(:,:,1) = RESHAPE(REAL(rotated_pts(:,:,1), gk), (/ nproma, nblks_lonlat /), (/ 0._gk /) )
      in_points(:,:,2) = RESHAPE(REAL(rotated_pts(:,:,2), gk), (/ nproma, nblks_lonlat /), (/ 0._gk /) )

      ! --------------------------------------------------------------
      ! discard lon/lat grid points which are "far off"
      ! --------------------------------------------------------------

      pts_flags(:,:) = INVALID_NODE
      CALL flag_ll_points(rotated_pts, 1,grid%lon_dim, 1,grid%lat_dim, pts_flags, 1,5, gnat, max_dist)

      !-- if we have a global, unrotated grid: skip all points of the
      !   first and last latitude row (except one):
      IF (l_grid_is_unrotated .AND. l_grid_contains_poles) THEN
        pts_flags(2:,           1) = SKIP_NODE
        pts_flags(2:,grid%lat_dim) = SKIP_NODE
      END IF

      ! --------------------------------------------------------------
      ! proximity query, build a list of cell indices that contain the
      ! lon/lat grid points
      ! --------------------------------------------------------------

      ! allocate global arrays for distributed computation:
      nblks_lonlat  = (grid%total_dim - 1)/nproma + 1
      npromz_lonlat = grid%total_dim - (nblks_lonlat-1)*nproma
      ALLOCATE(ptr_int_lonlat%tri_idx(2, nproma, nblks_lonlat),   &
        &      ptr_int_lonlat%global_idx(grid%total_dim), STAT=errstat )
      IF (errstat /= SUCCESS) CALL finish (routine, 'allocation for lon-lat point distribution failed')
      ptr_int_lonlat%tri_idx(1,:,:) = RESHAPE(pts_flags(:,:), (/ nproma, nblks_lonlat /), (/ INVALID_NODE /) )
      ptr_int_lonlat%tri_idx(2,:,:) = 0

!$    time1 = REAL(omp_get_wtime())

      ! Perform query. Note that for distributed patches we receive a
      ! local list of "in_points" actually located on this portion of the
      ! domain.
      IF (dbg_level > 1) CALL message(routine, "proximity query")
      start_radius = grid_sphere_radius
      ! work-around: increased search radius on torus
      if (is_plane_torus)  start_radius = 10._wp
      CALL gnat_query_containing_triangles(gnat, ptr_patch, in_points(:,:,:),      &
        &                 nproma, nblks_lonlat, npromz_lonlat, start_radius,       &
        &                 p_test_run, ptr_int_lonlat%tri_idx(:,:,:), min_dist(:,:))

!$    IF (l_measure_time) WRITE (0,*) "elapsed time (query): ",  REAL(omp_get_wtime()) - time1

      !-- if we have a global, unrotated grid: copy the distance
      !   result for all points of the first and last latitude row:
      IF (l_grid_is_unrotated .AND. l_grid_contains_poles) THEN
        row_idx = (/ 1, grid%lat_dim /)
        DO i=1,2
          ! compute indices for pole row
          ibeg_glb = (row_idx(i)-1)*grid%lon_dim
          ibeg_idx = idx_no(ibeg_glb + 1)
          ibeg_blk = blk_no(ibeg_glb + 1)
          iend_idx = idx_no(ibeg_glb + grid%lon_dim)
          iend_blk = blk_no(ibeg_glb + grid%lon_dim)
          ! copy indices and distances:
          DO j=1,2
            tri_idx_pole = ptr_int_lonlat%tri_idx(j,ibeg_idx,ibeg_blk)
            IF (iend_blk > ibeg_blk) THEN
              ptr_int_lonlat%tri_idx(j, (ibeg_idx+1):nproma,                  ibeg_blk) = tri_idx_pole
              ptr_int_lonlat%tri_idx(j,                   :, (ibeg_blk+1):(iend_blk-1)) = tri_idx_pole
              ptr_int_lonlat%tri_idx(j,          1:iend_idx,                  iend_blk) = tri_idx_pole
            ELSE
              ptr_int_lonlat%tri_idx(j, (ibeg_idx+1):iend_idx, ibeg_blk) = tri_idx_pole
            END IF
          END DO
          min_dist_pole = min_dist(ibeg_idx,ibeg_blk)
          IF (iend_blk > ibeg_blk) THEN
            min_dist((ibeg_idx+1):nproma,                  ibeg_blk) = min_dist_pole
            min_dist(                  :, (ibeg_blk+1):(iend_blk-1)) = min_dist_pole
            min_dist(         1:iend_idx,                  iend_blk) = min_dist_pole
          ELSE
            min_dist((ibeg_idx+1):iend_idx, ibeg_blk) = min_dist_pole
          END IF
        END DO
      END IF

      CALL gnat_merge_distributed_queries(ptr_patch, grid%total_dim, nproma, grid%nblks, &
        &                 min_dist, ptr_int_lonlat%tri_idx(:,:,:), in_points(:,:,:),     &
        &                 ptr_int_lonlat%global_idx(:), ptr_int_lonlat%nthis_local_pts)

!$    if (l_measure_time) WRITE (0,*) "elapsed time: ",  REAL(omp_get_wtime()) - time1

      ! set local values for "nblks" and "npromz"
      nblks_lonlat  = (ptr_int_lonlat%nthis_local_pts - 1)/nproma + 1
      npromz_lonlat = ptr_int_lonlat%nthis_local_pts - (nblks_lonlat-1)*nproma

      ! After the lon-lat points have been distributed, we know
      ! exactly for each PE the number of points
      CALL allocate_int_state_lonlat_grid( ptr_patch%nblks_c, nblks_lonlat, ptr_int_lonlat)

      ! store the geographical coordinates for the lon-lat points
      !
!$OMP PARALLEL DO PRIVATE(jb,jc,i_startidx,i_endidx,point)
      DO jb = 1,nblks_lonlat
        i_startidx = 1
        i_endidx   = nproma
        IF (jb == nblks_lonlat) i_endidx = npromz_lonlat

        DO jc = i_startidx, i_endidx
          point(1:2) = REAL(in_points(jc,jb,1:2), wp)
          ptr_int_lonlat%ll_coord(jc,jb)%lon = point(1)
          ptr_int_lonlat%ll_coord(jc,jb)%lat = point(2)
        END DO
      END DO
!$OMP END PARALLEL DO

      ! clean up
      DEALLOCATE(in_points, pts_flags, rotated_pts, min_dist, stat=errstat)
      IF (errstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE of working arrays failed')


      !-- edge indices of stencil
      ! are available through previous calls to SR rbf_vec_index_cell()
      ! and rbf_c2grad_index()

      ! -----------------------------------------------------
      ! copy neighbor indices from standard RBF interpolation
      ! -----------------------------------------------------

      IF (dbg_level > 1) CALL message(routine, "copy neighbor indices from standard RBF interpolation")

!$OMP PARALLEL SHARED(nblks_lonlat, nproma, npromz_lonlat,   &
!$OMP                 ptr_int_lonlat, ptr_int),              &
!$OMP          DEFAULT (NONE)
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jc_lonlat,jb_lonlat), SCHEDULE(runtime)
      DO jb_lonlat = 1,nblks_lonlat

        i_startidx = 1
        i_endidx   = nproma
        IF (jb_lonlat == nblks_lonlat) i_endidx = npromz_lonlat

        DO jc_lonlat = i_startidx, i_endidx

          jc = ptr_int_lonlat%tri_idx(1,jc_lonlat, jb_lonlat)
          jb = ptr_int_lonlat%tri_idx(2,jc_lonlat, jb_lonlat)

          ptr_int_lonlat%rbf_vec_stencil(jc_lonlat,jb_lonlat) = ptr_int%rbf_vec_stencil_c(jc,jb)
          ptr_int_lonlat%rbf_vec_idx(:,jc_lonlat,jb_lonlat)   = ptr_int%rbf_vec_idx_c(:,jc,jb)
          ptr_int_lonlat%rbf_vec_blk(:,jc_lonlat,jb_lonlat)   = ptr_int%rbf_vec_blk_c(:,jc,jb)
        ENDDO

      ENDDO
!$OMP END DO
!$OMP END PARALLEL

      IF (dbg_level > 1) CALL message(routine, "compute lon-lat interpolation coefficients")

      ! -------------------------------------------------------------------------
      ! compute interpolation coefficients for RBF interpolation of vector fields
      ! -------------------------------------------------------------------------

      SELECT CASE (rbf_scale_mode_ll)
      CASE (SCALE_MODE_TABLE) 
        rbf_shape_param = rbf_vec_scale_ll(MAX(ptr_patch%id,1))
      CASE (SCALE_MODE_AUTO)
        ! if no shape parameter has been set: compute an estimate 
        CALL estimate_rbf_parameter(nblks_lonlat, npromz_lonlat, ptr_patch%edges%center,              &
          &                         ptr_int_lonlat%rbf_vec_idx, ptr_int_lonlat%rbf_vec_blk,           &
          &                         ptr_int_lonlat%rbf_vec_stencil, rbf_vec_dim_c,                    &
          &                         ptr_int_lonlat%global_idx, rbf_shape_param)
        rbf_shape_param          = p_min(rbf_shape_param, comm=p_comm_work)
        ptr_int_lonlat%rbf_scale = rbf_shape_param
        IF (my_process_is_stdio()) THEN
          WRITE(0,*) routine, ": auto-estimated shape_param = ", rbf_shape_param
        END IF
      CASE (SCALE_MODE_PRESET)
        IF (ptr_int_lonlat%rbf_scale <= 0._wp) THEN
          CALL finish(routine, "Explicitly presetting RBF shape parameter... invalid value!")
        END IF
        rbf_shape_param = ptr_int_lonlat%rbf_scale
      CASE DEFAULT
        CALL finish(routine, "Unknown value for rbf_scale_mode_ll!")
      END SELECT
      ! fallback solution:
      IF (rbf_shape_param <= 0._wp) &
        &  rbf_shape_param = rbf_vec_scale_ll(MAX(ptr_patch%id,1))
      WRITE(message_text,*) ": chosen rbf_shape_param = ", rbf_shape_param
      CALL message(routine, TRIM(message_text))

      CALL rbf_vec_compute_coeff_lonlat( ptr_patch, ptr_int_lonlat, nblks_lonlat, &
        &                                npromz_lonlat, rbf_shape_param )

      ! -------------------------------------------------------------------------
      ! compute interpolation coefficients for RBF interpolation of scalar values
      ! -------------------------------------------------------------------------

      CALL rbf_c2l_index( ptr_patch, ptr_int, ptr_int_lonlat )
      ! Compute reordered index lists to avoid nested indirect addressing at runtime
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jb_lonlat,jc_lonlat), SCHEDULE(runtime)
      DO jb_lonlat = 1,nblks_lonlat

        i_startidx = 1
        i_endidx   = nproma
        IF (jb_lonlat == nblks_lonlat) i_endidx = npromz_lonlat
        
        DO jc_lonlat = i_startidx, i_endidx
          jc = ptr_int_lonlat%tri_idx(1,jc_lonlat, jb_lonlat)
          jb = ptr_int_lonlat%tri_idx(2,jc_lonlat, jb_lonlat)
          
          ptr_int_lonlat%rbf_c2lr_idx(:,jc_lonlat,jb_lonlat)     = ptr_int_lonlat%rbf_c2l_idx(:,jc,jb)
          ptr_int_lonlat%rbf_c2lr_blk(:,jc_lonlat,jb_lonlat)     = ptr_int_lonlat%rbf_c2l_blk(:,jc,jb)
          ptr_int_lonlat%rbf_c2lr_stencil(jc_lonlat,jb_lonlat)   = ptr_int_lonlat%rbf_c2l_stencil(jc,jb)
        ENDDO
        
      ENDDO
!$OMP END DO
!$OMP END PARALLEL

      SELECT CASE (rbf_scale_mode_ll)
      CASE (SCALE_MODE_TABLE) 
        rbf_shape_param = rbf_vec_scale_ll(MAX(ptr_patch%id,1))
      CASE (SCALE_MODE_AUTO)
        ! if no shape parameter has been set: compute an estimate 
        CALL estimate_rbf_parameter(nblks_lonlat, npromz_lonlat, ptr_patch%cells%center,              &
          &                         ptr_int_lonlat%rbf_c2lr_idx, ptr_int_lonlat%rbf_c2lr_blk,         &
          &                         ptr_int_lonlat%rbf_c2lr_stencil, rbf_dim_c2l,                     &
          &                         ptr_int_lonlat%global_idx, rbf_shape_param)
        rbf_shape_param = p_min(rbf_shape_param, comm=p_comm_work)
        ptr_int_lonlat%rbf_scale = rbf_shape_param
        IF (my_process_is_stdio()) THEN
          WRITE(0,*) routine, ": auto-estimated shape_param = ", rbf_shape_param
        END IF
      CASE (SCALE_MODE_PRESET)
        rbf_shape_param = ptr_int_lonlat%rbf_scale
      CASE DEFAULT
        CALL finish(routine, "Unknown value for rbf_scale_mode_ll!")
      END SELECT
      ! fallback solution:
      IF (rbf_shape_param <= 0._wp) &
        &  rbf_shape_param = rbf_vec_scale_ll(MAX(ptr_patch%id,1))
      WRITE(message_text,*) ": chosen rbf_shape_param = ", rbf_shape_param
      CALL message(routine, TRIM(message_text))
      ! compute coefficients:
      CALL rbf_compute_coeff_c2l( ptr_patch, ptr_int_lonlat, &
          &                         nblks_lonlat, npromz_lonlat, rbf_shape_param )

    END SUBROUTINE rbf_setup_interpol_lonlat_grid


    !===============================================================
    ! APPLY LON-LAT INTERPOLATION

    !-------------------------------------------------------------------------
    !> Performs vector RBF reconstruction at lon-lat grid points.
    !
    ! This routine is based on mo_intp_rbf::rbf_vec_interpol_cell()
    !
    ! @par Revision History
    !      Initial implementation  by  F. Prill, DWD (2011-08)
    !
    SUBROUTINE rbf_vec_interpol_lonlat( p_vn_in, ptr_int, &
      &                                 grad_x, grad_y,              &
      &                                 nblks_lonlat, npromz_lonlat, &
      &                                 opt_slev, opt_elev)

      ! INPUT PARAMETERS

      ! input normal components of (velocity) vectors at edge midpoints
      REAL(wp),                   INTENT(IN) :: p_vn_in(:,:,:)                ! dim: (nproma,nlev,nblks_e)
      ! Indices of source points and interpolation coefficients
      TYPE (t_lon_lat_intp), TARGET, INTENT(IN) :: ptr_int

      ! reconstructed x/y-components of velocity vector
      REAL(wp),         INTENT(INOUT)        :: grad_x(:,:,:), grad_y(:,:,:)  ! dim: (nproma,nlev,nblks_lonlat)

      INTEGER,          INTENT(IN)           :: nblks_lonlat, npromz_lonlat   ! blocking info

      INTEGER,          INTENT(in), OPTIONAL :: opt_slev, opt_elev            ! optional vertical start/end level

      ! LOCAL VARIABLES
      INTEGER :: slev, elev,                 & ! vertical start and end level
        &        i_startidx, i_endidx,       & ! start/end index
        &        jc, jb, jk                    ! integer over lon-lat points, levels

      INTEGER,  DIMENSION(:,:,:),   POINTER :: iidx, iblk
      REAL(wp), DIMENSION(:,:,:,:), POINTER :: ptr_coeff

      !-----------------------------------------------------------------------

      slev = 1
      elev = UBOUND(p_vn_in,2)
      ! check optional arguments
      IF ( PRESENT(opt_slev) ) slev = opt_slev
      IF ( PRESENT(opt_elev) ) elev = opt_elev

      iidx => ptr_int%rbf_vec_idx
      iblk => ptr_int%rbf_vec_blk
      ptr_coeff => ptr_int%rbf_vec_coeff

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk,jc), SCHEDULE(runtime)
      DO jb = 1,nblks_lonlat

        i_startidx = 1
        i_endidx   = nproma
        IF (jb == nblks_lonlat) i_endidx = npromz_lonlat

#ifdef __LOOP_EXCHANGE
        DO jc = i_startidx, i_endidx
          DO jk = slev, elev
#else
!CDIR UNROLL=2
        DO jk = slev, elev
          DO jc = i_startidx, i_endidx
#endif

            grad_x(jc,jk,jb) =                                               &
              ptr_coeff(1,1,jc,jb)*p_vn_in(iidx(1,jc,jb),jk,iblk(1,jc,jb)) + &
              ptr_coeff(2,1,jc,jb)*p_vn_in(iidx(2,jc,jb),jk,iblk(2,jc,jb)) + &
              ptr_coeff(3,1,jc,jb)*p_vn_in(iidx(3,jc,jb),jk,iblk(3,jc,jb)) + &
              ptr_coeff(4,1,jc,jb)*p_vn_in(iidx(4,jc,jb),jk,iblk(4,jc,jb)) + &
              ptr_coeff(5,1,jc,jb)*p_vn_in(iidx(5,jc,jb),jk,iblk(5,jc,jb)) + &
              ptr_coeff(6,1,jc,jb)*p_vn_in(iidx(6,jc,jb),jk,iblk(6,jc,jb)) + &
              ptr_coeff(7,1,jc,jb)*p_vn_in(iidx(7,jc,jb),jk,iblk(7,jc,jb)) + &
              ptr_coeff(8,1,jc,jb)*p_vn_in(iidx(8,jc,jb),jk,iblk(8,jc,jb)) + &
              ptr_coeff(9,1,jc,jb)*p_vn_in(iidx(9,jc,jb),jk,iblk(9,jc,jb))

            grad_y(jc,jk,jb) =                                               &
              ptr_coeff(1,2,jc,jb)*p_vn_in(iidx(1,jc,jb),jk,iblk(1,jc,jb)) + &
              ptr_coeff(2,2,jc,jb)*p_vn_in(iidx(2,jc,jb),jk,iblk(2,jc,jb)) + &
              ptr_coeff(3,2,jc,jb)*p_vn_in(iidx(3,jc,jb),jk,iblk(3,jc,jb)) + &
              ptr_coeff(4,2,jc,jb)*p_vn_in(iidx(4,jc,jb),jk,iblk(4,jc,jb)) + &
              ptr_coeff(5,2,jc,jb)*p_vn_in(iidx(5,jc,jb),jk,iblk(5,jc,jb)) + &
              ptr_coeff(6,2,jc,jb)*p_vn_in(iidx(6,jc,jb),jk,iblk(6,jc,jb)) + &
              ptr_coeff(7,2,jc,jb)*p_vn_in(iidx(7,jc,jb),jk,iblk(7,jc,jb)) + &
              ptr_coeff(8,2,jc,jb)*p_vn_in(iidx(8,jc,jb),jk,iblk(8,jc,jb)) + &
              ptr_coeff(9,2,jc,jb)*p_vn_in(iidx(9,jc,jb),jk,iblk(9,jc,jb))

          ENDDO
        ENDDO

      ENDDO
!$OMP END DO
!$OMP END PARALLEL

    END SUBROUTINE rbf_vec_interpol_lonlat


    !-------------------------------------------------------------------------
    !> Performs nearest neighbor interpolation, REAL implementation
    !
    ! @par Revision History
    !      Initial implementation  by  F. Prill, DWD (2013-02)
    !
    SUBROUTINE nnb_interpol_lonlat_real( p_cell_in, ptr_int,                   &
      &                                  p_out, nblks_lonlat, npromz_lonlat,   &
      &                                  opt_slev, opt_elev)
      ! INPUT PARAMETERS
      !
      ! input cell-based variable for which gradient at cell center is computed
      REAL(wp),                      INTENT(IN)    :: p_cell_in(:,:,:) ! dim: (nproma,nlev,nblks_c)
      ! Indices of source points and interpolation coefficients
      TYPE (t_lon_lat_intp), TARGET, INTENT(IN)    :: ptr_int
      ! reconstructed scalar value at lon-lat point
      REAL(wp),                      INTENT(INOUT) :: p_out(:,:,:)     ! dim: (nproma,nlev,nblks_lonlat)
      ! lon-lat grid blocking info
      INTEGER,                       INTENT(IN)    :: nblks_lonlat, npromz_lonlat
      ! optional vertical start/end level
      INTEGER,                       INTENT(IN), OPTIONAL :: opt_slev, opt_elev

      ! LOCAL VARIABLES
      INTEGER :: slev, elev,                 & ! vertical start and end level
        &        i_startidx, i_endidx,       & ! start/end index
        &        jc, jb, jk                    ! integer over lon-lat points, levels

      slev = 1
      elev = UBOUND(p_cell_in,2)
      ! check optional arguments
      IF ( PRESENT(opt_slev) ) slev = opt_slev
      IF ( PRESENT(opt_elev) ) elev = opt_elev

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk,jc), SCHEDULE(runtime)
      DO jb = 1,nblks_lonlat
        i_startidx = 1
        i_endidx   = nproma
        IF (jb == nblks_lonlat) i_endidx = npromz_lonlat
#ifdef __LOOP_EXCHANGE
          DO jc = i_startidx, i_endidx
            DO jk = slev, elev
#else
!CDIR UNROLL=3
          DO jk = slev, elev
            DO jc = i_startidx, i_endidx
#endif
              p_out(jc,jk,jb) = ptr_int%nnb_coeff(1,jc,jb) * &
                &  p_cell_in(ptr_int%tri_idx(1,jc,jb), jk, ptr_int%tri_idx(2,jc,jb))
            ENDDO
          ENDDO
      END DO
!$OMP END DO
!$OMP END PARALLEL
    END SUBROUTINE nnb_interpol_lonlat_real


    !-------------------------------------------------------------------------
    !> Performs nearest neighbor interpolation, INTEGER implementation
    !
    ! @par Revision History
    !      Initial implementation  by  F. Prill, DWD (2013-02)
    !
    SUBROUTINE nnb_interpol_lonlat_int( p_cell_in, ptr_int,                   &
      &                                 p_out, nblks_lonlat, npromz_lonlat,   &
      &                                 opt_slev, opt_elev)
      ! INPUT PARAMETERS
      !
      ! input cell-based variable for which gradient at cell center is computed
      INTEGER,                       INTENT(IN)    :: p_cell_in(:,:,:) ! dim: (nproma,nlev,nblks_c)
      ! Indices of source points and interpolation coefficients
      TYPE (t_lon_lat_intp), TARGET, INTENT(IN)    :: ptr_int
      ! reconstructed scalar value at lon-lat point
      INTEGER,                       INTENT(INOUT) :: p_out(:,:,:)     ! dim: (nproma,nlev,nblks_lonlat)
      ! lon-lat grid blocking info
      INTEGER,                       INTENT(IN)    :: nblks_lonlat, npromz_lonlat
      ! optional vertical start/end level
      INTEGER,                       INTENT(IN), OPTIONAL :: opt_slev, opt_elev

      ! LOCAL VARIABLES
      INTEGER :: slev, elev,                 & ! vertical start and end level
        &        i_startidx, i_endidx,       & ! start/end index
        &        jc, jb, jk                    ! integer over lon-lat points, levels

      slev = 1
      elev = UBOUND(p_cell_in,2)
      ! check optional arguments
      IF ( PRESENT(opt_slev) ) slev = opt_slev
      IF ( PRESENT(opt_elev) ) elev = opt_elev

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk,jc), SCHEDULE(runtime)
      DO jb = 1,nblks_lonlat
        i_startidx = 1
        i_endidx   = nproma
        IF (jb == nblks_lonlat) i_endidx = npromz_lonlat
#ifdef __LOOP_EXCHANGE
          DO jc = i_startidx, i_endidx
            DO jk = slev, elev
#else
!CDIR UNROLL=3
          DO jk = slev, elev
            DO jc = i_startidx, i_endidx
#endif
              p_out(jc,jk,jb) = INT(ptr_int%nnb_coeff(1,jc,jb)) *  &
                & p_cell_in(ptr_int%tri_idx(1,jc,jb), jk, ptr_int%tri_idx(2,jc,jb))
            ENDDO
          ENDDO
      END DO
!$OMP END DO
!$OMP END PARALLEL
    END SUBROUTINE nnb_interpol_lonlat_int
   

    !-------------------------------------------------------------------------
    !> Performs vector RBF reconstruction at lon-lat grid points.
    !  This routine takes a cell-based variable as input and reconstructs
    !  the gradient (of the FD approximation) at the lon-lat grid points.
    !
    !  More precisely, the routine "rbf_interpol_c2grad_lonlat" combines
    !  the computation of gradients from scalar variables,
    !  CALL grad_fd_norm( p_cell_in(:,:,:), &
    !    &                ptr_patch, grad_norm_psi_e )
    !
    !  and the application of interpolation coefficients,
    !  CALL rbf_vec_interpol_lonlat( grad_norm_psi_e, ptr_patch, ptr_int, grad_x,  &
    !    &                           grad_y, ptr_int%tri_idx(:,:,:),        &
    !    &                           nblks_lonlat, npromz_lonlat )
    !
    ! @par Revision History
    !      Initial implementation  by  F. Prill, DWD (2011-08)
    !      based on "rbf_interpol_c2grad"
    !
    SUBROUTINE rbf_interpol_c2l( p_cell_in, ptr_int,                   &
      &                          p_out, nblks_lonlat, npromz_lonlat,   &
      &                          opt_slev, opt_elev)
      ! !INPUT PARAMETERS
      !
      ! input cell-based variable for which gradient at cell center is computed
      REAL(wp),                  INTENT(IN) :: p_cell_in(:,:,:) ! dim: (nproma,nlev,nblks_c)
      ! Indices of source points and interpolation coefficients
      TYPE (t_lon_lat_intp), TARGET, INTENT(IN) :: ptr_int
      ! reconstructed scalar value at lon-lat point
      REAL(wp),INTENT(INOUT) :: p_out(:,:,:) ! dim: (nproma,nlev,nblks_lonlat)

      ! lon-lat grid blocking info
      INTEGER,                   INTENT(IN) :: nblks_lonlat, npromz_lonlat
      ! optional vertical start/end level
      INTEGER,                   INTENT(IN), OPTIONAL :: opt_slev, opt_elev

      ! Local variables
      INTEGER :: slev, elev,               &  ! vertical start and end level
        &        jc, jb, jk,               &  ! integer over lon-lat points, levels
        &        i_startidx, i_endidx         ! start/end index

      REAL(wp) :: vmin, vmax
      INTEGER,  DIMENSION(:,:,:), POINTER :: iidx, iblk
      REAL(wp), DIMENSION(:,:,:), POINTER :: ptr_coeff

      !-----------------------------------------------------------------------

      slev = 1
      elev = UBOUND(p_cell_in,2)
      ! check optional arguments
      IF ( PRESENT(opt_slev) ) slev = opt_slev
      IF ( PRESENT(opt_elev) ) elev = opt_elev

      iidx => ptr_int%rbf_c2lr_idx
      iblk => ptr_int%rbf_c2lr_blk

      ptr_coeff => ptr_int%rbf_c2l_coeff

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk,jc,vmin,vmax), SCHEDULE(runtime)

      DO jb = 1,nblks_lonlat

        i_startidx = 1
        i_endidx   = nproma
        IF (jb == nblks_lonlat) i_endidx = npromz_lonlat

        ! we have to duplicate code here for different stencil sizes,
        ! otherwise we would break vectorization...
        SELECT CASE(rbf_dim_c2l)

        CASE(4)

#ifdef __LOOP_EXCHANGE
          DO jc = i_startidx, i_endidx
            DO jk = slev, elev
#else
!CDIR UNROLL=3
          DO jk = slev, elev
            DO jc = i_startidx, i_endidx
#endif

              p_out(jc,jk,jb) =                                                   &
                ptr_coeff(1 ,jc,jb)*p_cell_in(iidx(1 ,jc,jb),jk,iblk(1 ,jc,jb)) + &
                ptr_coeff(2 ,jc,jb)*p_cell_in(iidx(2 ,jc,jb),jk,iblk(2 ,jc,jb)) + &
                ptr_coeff(3 ,jc,jb)*p_cell_in(iidx(3 ,jc,jb),jk,iblk(3 ,jc,jb)) + &
                ptr_coeff(4 ,jc,jb)*p_cell_in(iidx(4 ,jc,jb),jk,iblk(4 ,jc,jb))

              ! monotonicity can be enforced by demanding that the interpolated
              ! value is not higher or lower than the stencil point values.

              ! Cf. the "lmono" implementation in the GME:
              ! D. Majewski, "Documentation of the new global model (GME)
              !               of the DWD" (1996)
              IF (l_mono_c2l) THEN

                vmin = MIN(                                     &
                  p_cell_in(iidx(1 ,jc,jb),jk,iblk(1 ,jc,jb)) , &
                  p_cell_in(iidx(2 ,jc,jb),jk,iblk(2 ,jc,jb)) , &
                  p_cell_in(iidx(3 ,jc,jb),jk,iblk(3 ,jc,jb)) , &
                  p_cell_in(iidx(4 ,jc,jb),jk,iblk(4 ,jc,jb))   )

                vmax = MAX(                                     &
                  p_cell_in(iidx(1 ,jc,jb),jk,iblk(1 ,jc,jb)) , &
                  p_cell_in(iidx(2 ,jc,jb),jk,iblk(2 ,jc,jb)) , &
                  p_cell_in(iidx(3 ,jc,jb),jk,iblk(3 ,jc,jb)) , &
                  p_cell_in(iidx(4 ,jc,jb),jk,iblk(4 ,jc,jb))   )

                p_out(jc,jk,jb) = MAX( MIN(p_out(jc,jk,jb), vmax), vmin )
              END IF

            ENDDO
          ENDDO

        CASE(10)

#ifdef __LOOP_EXCHANGE
          DO jc = i_startidx, i_endidx
            DO jk = slev, elev
#else
          DO jk = slev, elev
            DO jc = i_startidx, i_endidx
#endif

              p_out(jc,jk,jb) =                                                   &
                ptr_coeff(1 ,jc,jb)*p_cell_in(iidx(1 ,jc,jb),jk,iblk(1 ,jc,jb)) + &
                ptr_coeff(2 ,jc,jb)*p_cell_in(iidx(2 ,jc,jb),jk,iblk(2 ,jc,jb)) + &
                ptr_coeff(3 ,jc,jb)*p_cell_in(iidx(3 ,jc,jb),jk,iblk(3 ,jc,jb)) + &
                ptr_coeff(4 ,jc,jb)*p_cell_in(iidx(4 ,jc,jb),jk,iblk(4 ,jc,jb)) + &
                ptr_coeff(5 ,jc,jb)*p_cell_in(iidx(5 ,jc,jb),jk,iblk(5 ,jc,jb)) + &
                ptr_coeff(6 ,jc,jb)*p_cell_in(iidx(6 ,jc,jb),jk,iblk(6 ,jc,jb)) + &
                ptr_coeff(7 ,jc,jb)*p_cell_in(iidx(7 ,jc,jb),jk,iblk(7 ,jc,jb)) + &
                ptr_coeff(8 ,jc,jb)*p_cell_in(iidx(8 ,jc,jb),jk,iblk(8 ,jc,jb)) + &
                ptr_coeff(9 ,jc,jb)*p_cell_in(iidx(9 ,jc,jb),jk,iblk(9 ,jc,jb)) + &
                ptr_coeff(10,jc,jb)*p_cell_in(iidx(10,jc,jb),jk,iblk(10,jc,jb))

              ! monotonicity can be enforced by demanding that the interpolated
              ! value is not higher or lower than the stencil point values.

              ! Cf. the "lmono" implementation in the GME:
              ! D. Majewski, "Documentation of the new global model (GME)
              !               of the DWD" (1996)
              IF (l_mono_c2l) THEN

                vmin = MIN(                                     &
                  p_cell_in(iidx(1 ,jc,jb),jk,iblk(1 ,jc,jb)) , &
                  p_cell_in(iidx(2 ,jc,jb),jk,iblk(2 ,jc,jb)) , &
                  p_cell_in(iidx(3 ,jc,jb),jk,iblk(3 ,jc,jb)) , &
                  p_cell_in(iidx(4 ,jc,jb),jk,iblk(4 ,jc,jb)) , &
                  p_cell_in(iidx(5 ,jc,jb),jk,iblk(5 ,jc,jb)) , &
                  p_cell_in(iidx(6 ,jc,jb),jk,iblk(6 ,jc,jb)) , &
                  p_cell_in(iidx(7 ,jc,jb),jk,iblk(7 ,jc,jb)) , &
                  p_cell_in(iidx(8 ,jc,jb),jk,iblk(8 ,jc,jb)) , &
                  p_cell_in(iidx(9 ,jc,jb),jk,iblk(9 ,jc,jb)) , &
                  p_cell_in(iidx(10,jc,jb),jk,iblk(10,jc,jb))   )

                vmax = MAX(                                     &
                  p_cell_in(iidx(1 ,jc,jb),jk,iblk(1 ,jc,jb)) , &
                  p_cell_in(iidx(2 ,jc,jb),jk,iblk(2 ,jc,jb)) , &
                  p_cell_in(iidx(3 ,jc,jb),jk,iblk(3 ,jc,jb)) , &
                  p_cell_in(iidx(4 ,jc,jb),jk,iblk(4 ,jc,jb)) , &
                  p_cell_in(iidx(5 ,jc,jb),jk,iblk(5 ,jc,jb)) , &
                  p_cell_in(iidx(6 ,jc,jb),jk,iblk(6 ,jc,jb)) , &
                  p_cell_in(iidx(7 ,jc,jb),jk,iblk(7 ,jc,jb)) , &
                  p_cell_in(iidx(8 ,jc,jb),jk,iblk(8 ,jc,jb)) , &
                  p_cell_in(iidx(9 ,jc,jb),jk,iblk(9 ,jc,jb)) , &
                  p_cell_in(iidx(10,jc,jb),jk,iblk(10,jc,jb))   )

                p_out(jc,jk,jb) = MAX( MIN(p_out(jc,jk,jb), vmax), vmin )
              END IF

            ENDDO
          ENDDO

        CASE(13)

#ifdef __LOOP_EXCHANGE
          DO jc = i_startidx, i_endidx
            DO jk = slev, elev
#else
          DO jk = slev, elev
            DO jc = i_startidx, i_endidx
#endif

              p_out(jc,jk,jb) =                                                   &
                ptr_coeff(1 ,jc,jb)*p_cell_in(iidx(1 ,jc,jb),jk,iblk(1 ,jc,jb)) + &
                ptr_coeff(2 ,jc,jb)*p_cell_in(iidx(2 ,jc,jb),jk,iblk(2 ,jc,jb)) + &
                ptr_coeff(3 ,jc,jb)*p_cell_in(iidx(3 ,jc,jb),jk,iblk(3 ,jc,jb)) + &
                ptr_coeff(4 ,jc,jb)*p_cell_in(iidx(4 ,jc,jb),jk,iblk(4 ,jc,jb)) + &
                ptr_coeff(5 ,jc,jb)*p_cell_in(iidx(5 ,jc,jb),jk,iblk(5 ,jc,jb)) + &
                ptr_coeff(6 ,jc,jb)*p_cell_in(iidx(6 ,jc,jb),jk,iblk(6 ,jc,jb)) + &
                ptr_coeff(7 ,jc,jb)*p_cell_in(iidx(7 ,jc,jb),jk,iblk(7 ,jc,jb)) + &
                ptr_coeff(8 ,jc,jb)*p_cell_in(iidx(8 ,jc,jb),jk,iblk(8 ,jc,jb)) + &
                ptr_coeff(9 ,jc,jb)*p_cell_in(iidx(9 ,jc,jb),jk,iblk(9 ,jc,jb)) + &
                ptr_coeff(10,jc,jb)*p_cell_in(iidx(10,jc,jb),jk,iblk(10,jc,jb)) + &
                ptr_coeff(11,jc,jb)*p_cell_in(iidx(11,jc,jb),jk,iblk(11,jc,jb)) + &
                ptr_coeff(12,jc,jb)*p_cell_in(iidx(12,jc,jb),jk,iblk(12,jc,jb)) + &
                ptr_coeff(13,jc,jb)*p_cell_in(iidx(13,jc,jb),jk,iblk(13,jc,jb))

              ! monotonicity can be enforced by demanding that the interpolated
              ! value is not higher or lower than the stencil point values.

              ! Cf. the "lmono" implementation in the GME:
              ! D. Majewski, "Documentation of the new global model (GME)
              !               of the DWD" (1996)
              IF (l_mono_c2l) THEN

                vmin = MIN(                                     &
                  p_cell_in(iidx(1 ,jc,jb),jk,iblk(1 ,jc,jb)) , &
                  p_cell_in(iidx(2 ,jc,jb),jk,iblk(2 ,jc,jb)) , &
                  p_cell_in(iidx(3 ,jc,jb),jk,iblk(3 ,jc,jb)) , &
                  p_cell_in(iidx(4 ,jc,jb),jk,iblk(4 ,jc,jb)) , &
                  p_cell_in(iidx(5 ,jc,jb),jk,iblk(5 ,jc,jb)) , &
                  p_cell_in(iidx(6 ,jc,jb),jk,iblk(6 ,jc,jb)) , &
                  p_cell_in(iidx(7 ,jc,jb),jk,iblk(7 ,jc,jb)) , &
                  p_cell_in(iidx(8 ,jc,jb),jk,iblk(8 ,jc,jb)) , &
                  p_cell_in(iidx(9 ,jc,jb),jk,iblk(9 ,jc,jb)) , &
                  p_cell_in(iidx(10,jc,jb),jk,iblk(10,jc,jb)) , &
                  p_cell_in(iidx(11,jc,jb),jk,iblk(11,jc,jb)) , &
                  p_cell_in(iidx(12,jc,jb),jk,iblk(12,jc,jb)) , &
                  p_cell_in(iidx(13,jc,jb),jk,iblk(13,jc,jb))   )

                vmax = MAX(                                     &
                  p_cell_in(iidx(1 ,jc,jb),jk,iblk(1 ,jc,jb)) , &
                  p_cell_in(iidx(2 ,jc,jb),jk,iblk(2 ,jc,jb)) , &
                  p_cell_in(iidx(3 ,jc,jb),jk,iblk(3 ,jc,jb)) , &
                  p_cell_in(iidx(4 ,jc,jb),jk,iblk(4 ,jc,jb)) , &
                  p_cell_in(iidx(5 ,jc,jb),jk,iblk(5 ,jc,jb)) , &
                  p_cell_in(iidx(6 ,jc,jb),jk,iblk(6 ,jc,jb)) , &
                  p_cell_in(iidx(7 ,jc,jb),jk,iblk(7 ,jc,jb)) , &
                  p_cell_in(iidx(8 ,jc,jb),jk,iblk(8 ,jc,jb)) , &
                  p_cell_in(iidx(9 ,jc,jb),jk,iblk(9 ,jc,jb)) , &
                  p_cell_in(iidx(10,jc,jb),jk,iblk(10,jc,jb)) , &
                  p_cell_in(iidx(11,jc,jb),jk,iblk(11,jc,jb)) , &
                  p_cell_in(iidx(12,jc,jb),jk,iblk(12,jc,jb)) , &
                  p_cell_in(iidx(13,jc,jb),jk,iblk(13,jc,jb))   )

                p_out(jc,jk,jb) = MAX( MIN(p_out(jc,jk,jb), vmax), vmin )
              END IF

            ENDDO
          ENDDO

        END SELECT

      ENDDO
!$OMP END DO
!$OMP END PARALLEL

    END SUBROUTINE rbf_interpol_c2l


    !-------------------------------------------------------------------------
    !> REAL fields: Driver routine for RBF reconstruction of
    !  cell-based variables at lon-lat grid points.
    !
    ! @par Revision History
    !      Initial implementation  by  F. Prill, DWD (2011-08)
    !
    SUBROUTINE interpol_lonlat_real( name, p_cell_in, ptr_int, p_lonlat_out, nblks_lonlat, &
      &                              npromz_lonlat, hintp_type)
      CHARACTER(LEN=*),      INTENT(IN)           :: name                !< variable name
      ! input cell-based variable for which gradient at cell center is computed
      REAL(wp),              INTENT(IN)           :: p_cell_in(:,:,:)    ! dim: (nproma,nlev,nblks_c)
      ! Indices of source points and interpolation coefficients
      TYPE (t_lon_lat_intp), TARGET, INTENT(IN)   :: ptr_int
      ! output lon-lat-based variable
      REAL(wp),              INTENT(INOUT)        :: p_lonlat_out(:,:,:) ! dim: (nproma,nlev,nblks_lonlat)
      ! lon-lat grid blocking info
      INTEGER,               INTENT(IN)           :: nblks_lonlat, npromz_lonlat
      ! optional: horizontal interpolation type
      INTEGER,               INTENT(IN)           :: hintp_type

      ! Local Parameters:
      CHARACTER(*), PARAMETER :: routine = modname//"::interpol_lonlat_real"
      INTEGER  :: slev, elev

      !-----------------------------------------------------------------------

      slev = 1
      elev = UBOUND(p_cell_in,2)

      !-- apply interpolation coefficients
      IF (dbg_level > 1) THEN
        WRITE(message_text,*) "PE #", p_pe, ": apply interpolation coefficients"
        CALL message(routine, TRIM(message_text))
      END IF

      SELECT CASE(hintp_type)
      CASE (HINTP_TYPE_LONLAT_RBF)

        ! ---------------------------------------------------------------
        ! RBF interpolation
        ! ---------------------------------------------------------------

        CALL rbf_interpol_c2l( p_cell_in(:,:,:), ptr_int, p_lonlat_out(:,:,:), &
          &                    nblks_lonlat, npromz_lonlat, slev, elev)

      CASE (HINTP_TYPE_LONLAT_NNB)

        ! ---------------------------------------------------------------
        ! Nearest-neighbor interpolation
        ! ---------------------------------------------------------------

        CALL nnb_interpol_lonlat( p_cell_in(:,:,:), ptr_int, p_lonlat_out(:,:,:), &
          &                       nblks_lonlat, npromz_lonlat, slev, elev)

      CASE (HINTP_TYPE_LONLAT_BCTR)

        ! ---------------------------------------------------------------
        ! Barycentric interpolation
        ! ---------------------------------------------------------------

        CALL baryctr_interpol_lonlat( p_cell_in(:,:,:), ptr_int, p_lonlat_out(:,:,:), &
          &                           nblks_lonlat, npromz_lonlat, slev, elev)

      CASE DEFAULT
        CALL finish(routine, "Internal error with variable "//TRIM(name))

      END SELECT

    END SUBROUTINE interpol_lonlat_real


    !-------------------------------------------------------------------------
    !> INTEGER fields: Driver routine for RBF reconstruction of
    !  cell-based variables at lon-lat grid points.
    !
    ! @par Revision History
    !      Initial implementation  by  F. Prill, DWD (2013-02)
    !
    SUBROUTINE interpol_lonlat_int( name, p_cell_in, ptr_int, p_lonlat_out, nblks_lonlat, &
      &                             npromz_lonlat, hintp_type)
      CHARACTER(LEN=*),      INTENT(IN)           :: name                !< variable name
      ! input cell-based variable for which gradient at cell center is computed
      INTEGER,               INTENT(IN)           :: p_cell_in(:,:,:)    ! dim: (nproma,nlev,nblks_c)
      ! Indices of source points and interpolation coefficients
      TYPE (t_lon_lat_intp), TARGET, INTENT(IN)   :: ptr_int
      ! output lon-lat-based variable
      INTEGER,               INTENT(INOUT)        :: p_lonlat_out(:,:,:) ! dim: (nproma,nlev,nblks_lonlat)
      ! lon-lat grid blocking info
      INTEGER,               INTENT(IN)           :: nblks_lonlat, npromz_lonlat
      ! optional: horizontal interpolation type
      INTEGER,               INTENT(IN)           :: hintp_type

      ! Local Parameters:
      CHARACTER(*), PARAMETER :: routine = modname//"::interpol_lonlat_int"
      INTEGER :: slev, elev

      slev = 1
      elev = UBOUND(p_cell_in,2)

      !-- apply interpolation coefficients
      IF (dbg_level > 1) THEN
        WRITE(message_text,*) "PE #", p_pe, ": apply interpolation coefficients"
        CALL message(routine, message_text)
      END IF

      SELECT CASE(hintp_type)
      
      CASE (HINTP_TYPE_LONLAT_NNB)

        ! ---------------------------------------------------------------
        ! Nearest-neighbor interpolation
        ! ---------------------------------------------------------------

        CALL nnb_interpol_lonlat( p_cell_in(:,:,:), ptr_int, p_lonlat_out(:,:,:), &
          &                       nblks_lonlat, npromz_lonlat, slev, elev)

      CASE DEFAULT
        CALL finish(routine, "Internal error with variable "//TRIM(name))

      END SELECT

    END SUBROUTINE interpol_lonlat_int

  END MODULE mo_intp_lonlat
