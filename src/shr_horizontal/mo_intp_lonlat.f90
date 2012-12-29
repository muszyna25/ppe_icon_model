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
  !! @par Copyright
  !! 2002-2007 by DWD and MPI-M
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
    USE mo_kind,                ONLY: wp
    USE mo_exception,           ONLY: message, message_text, finish
    USE mo_impl_constants,      ONLY: SUCCESS, min_rlcell_int, max_dom, &
      &                               HINTP_TYPE_NONE
    USE mo_model_domain,        ONLY: t_patch
    USE mo_run_config,          ONLY: ltimer
    USE mo_grid_config,         ONLY: n_dom, grid_sphere_radius
    USE mo_timer,               ONLY: timer_start, timer_stop, &
      &                               timers_level,            &
      &                               timer_lonlat_setup
    USE mo_math_utilities,      ONLY: gc2cc, gvec2cvec, arc_length_v, &
      &                               t_cartesian_coordinates,        &
      &                               t_geographical_coordinates,     &
      &                               latlon_compute_area_weights
    USE mo_math_utility_solvers, ONLY: solve_chol_v, choldec_v
    USE mo_lonlat_grid,         ONLY: t_lon_lat_grid,                         &
      &                               compute_lonlat_blocking,                &
      &                               compute_lonlat_specs
    USE mo_parallel_config,     ONLY: nproma, p_test_run
    USE mo_loopindices,         ONLY: get_indices_c, get_indices_e, get_indices_v
    USE mo_intp_data_strc,      ONLY: t_int_state, t_lon_lat_intp, n_lonlat_grids, &
      &                               lonlat_grid_list, n_lonlat_grids, MAX_LONLAT_GRIDS
    USE mo_interpol_config,     ONLY: rbf_vec_dim_c, rbf_c2grad_dim, rbf_vec_kern_ll,   &
      &                               rbf_vec_scale_ll, rbf_dim_c2l, l_intp_c2l,        &
      &                               l_mono_c2l
    USE mo_gnat_gridsearch,     ONLY: gnat_init_grid, gnat_destroy, gnat_tree,&
      &                               gnat_query_containing_triangles,        &
      &                               gnat_merge_distributed_queries, gk
    USE mo_math_utilities,      ONLY: rotate_latlon_grid
    USE mo_mpi,                 ONLY: my_process_is_mpi_workroot,                 &
      &                               get_my_mpi_work_id, p_n_work,               &
      &                               p_max, get_my_mpi_work_communicator,        &
      &                               my_process_is_mpi_seq, p_comm_work,         &
      &                               my_process_is_mpi_test, p_max, p_send,      &
      &                               p_recv, process_mpi_all_test_id,            &
      &                               process_mpi_all_workroot_id, p_pe
    USE mo_communication,       ONLY: idx_1d, blk_no, idx_no, &
      &                               setup_comm_pattern
    USE mo_lonlat_grid,         ONLY: t_lon_lat_grid
    USE mo_cf_convention,       ONLY: t_cf_var
    USE mo_grib2,               ONLY: t_grib2_var
    USE mo_cdi_constants,       ONLY: GRID_REGULAR_LONLAT, GRID_REFERENCE,  &
      &                               GRID_CELL, ZA_SURFACE, TIME_CONSTANT, &
      &                               TSTEP_CONSTANT, DATATYPE_PACK16, DATATYPE_FLT32
    USE mo_nonhydro_state,      ONLY: p_nh_state
    USE mo_var_list,            ONLY: add_var, create_hor_interp_metadata
    USE mo_linked_list,         ONLY: t_list_element
    USE mo_sync,                ONLY: SYNC_C, sync_idx, sync_patch_array

    IMPLICIT NONE

    !> level of output verbosity
    INTEGER, PARAMETER  :: dbg_level = 0

    PRIVATE

    PUBLIC :: init_lonlat_grid_list
    PUBLIC :: destroy_lonlat_grid_list
    PUBLIC :: compute_lonlat_intp_coeffs
    PUBLIC :: compute_lonlat_area_weights
    PUBLIC :: rbf_vec_compute_coeff_lonlat
    PUBLIC :: rbf_compute_coeff_c2grad_lonlat
    PUBLIC :: rbf_setup_interpol_lonlat_grid
    PUBLIC :: rbf_vec_interpol_lonlat_nl
    PUBLIC :: rbf_interpol_c2grad_lonlat_nl
    PUBLIC :: rbf_interpol_lonlat_nl
    PUBLIC :: allocate_int_state_lonlat_grid
    PUBLIC :: deallocate_int_state_lonlat

  CONTAINS

#include "intp_functions.inc"

    !===============================================================
    ! SETUP OF LON-LAT GRID LIST

    !---------------------------------------------------------------
    !> Setup of lon-lat registry
    !
    SUBROUTINE init_lonlat_grid_list()
      INTEGER :: i
      CHARACTER(*), PARAMETER :: routine = &
        &  TRIM("mo_intp_lonlat:init_lonlat_grid_list")

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
      CHARACTER(*), PARAMETER :: routine = &
        &  TRIM("mo_intp_lonlat:destroy_lonlat_grid_list")
      
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
      CHARACTER(*), PARAMETER :: routine = &
        &  TRIM("mo_intp_lonlat:resize_lonlat_state")
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
      CHARACTER(*), PARAMETER :: routine = &
        &  TRIM("mo_intp_lonlat:compute_lonlat_intp_coeffs")
      INTEGER              :: i, j, jg, n_points,  &
        &                     ist, glb_idx, my_id, nthis_local_pts
      INTEGER, ALLOCATABLE :: glb_owner(:), local_index(:)
      TYPE(t_lon_lat_grid), POINTER  :: grid

      IF (dbg_level > 5) CALL message(routine, "Enter")
      DO i=1, n_lonlat_grids

        ! compute some entries of lon-lat grid specification:
        CALL compute_lonlat_specs(lonlat_grid_list(i)%grid)
 
        CALL compute_lonlat_blocking(lonlat_grid_list(i)%grid, nproma)
        
        ! allocate and compute coefficients needed for lon-lat
        ! interpolation:
        lonlat_grid_list(i)%l_dom(n_dom+1:) = .FALSE.
        DO jg=1,n_dom
          IF (lonlat_grid_list(i)%l_dom(jg)) THEN

            IF (ltimer) CALL timer_start(timer_lonlat_setup)
            CALL rbf_setup_interpol_lonlat_grid(      &
              &         lonlat_grid_list(i)%grid,     &
              &         p_patch(jg),                  &
              &         lonlat_grid_list(i)%intp(jg), &
              &         p_int_state(jg) )
            IF (ltimer) CALL timer_stop(timer_lonlat_setup)
            lonlat_grid_list(i)%l_initialized(jg) = .TRUE.
          END IF
        END DO ! jg
        
        ! set communication pattern for gathering lon-lat grids:
        DO jg=1,n_dom
          IF (lonlat_grid_list(i)%l_dom(jg)) THEN
            ! n_points     : total number of points in the RECEIVER array
            grid => lonlat_grid_list(i)%grid
            n_points = grid%lon_dim * grid%lat_dim
            ALLOCATE( glb_owner(n_points), local_index(n_points), STAT=ist )
            IF (ist /= SUCCESS) &
              CALL finish (routine, 'allocation for working arrays failed')

            my_id = get_my_mpi_work_id()
            nthis_local_pts = lonlat_grid_list(i)%intp(jg)%nthis_local_pts

            ! glb_owner    : owner PE number of every point in the RECEIVER array
            glb_owner(:) = -1
            ! local_index  : local index in the SENDER array.
            ! positive on all positions owner by sender
            local_index(:) = -1
            DO j=1,nthis_local_pts
              glb_idx = lonlat_grid_list(i)%intp(jg)%global_idx(j)
              glb_owner  (glb_idx) = my_id
              local_index(glb_idx) = j
            END DO

            ! on work root PE we have to gather owner ranks for all points
            glb_owner = p_max(glb_owner, comm=get_my_mpi_work_communicator(), root=0)

            IF (.NOT. my_process_is_mpi_seq()) THEN
              IF (my_process_is_mpi_workroot()) THEN
                CALL setup_comm_pattern(n_points, glb_owner, local_index=local_index, &
                  &                     p_pat=lonlat_grid_list(i)%p_pat(jg))
              ELSE
                ! We don't want to receive any data, i.e. the number of
                ! lon-lat points is 0 and owner/global index are
                ! dummies!
                glb_owner(:) = -1
                CALL setup_comm_pattern(0, glb_owner, local_index=local_index, &
                  &                     p_pat=lonlat_grid_list(i)%p_pat(jg))
              END IF
            END IF

            DEALLOCATE( glb_owner, local_index, STAT=ist )
            IF (ist /= SUCCESS) &
              CALL finish (routine, 'deallocation for working arrays failed')

            ! resize global data structures, after the setup only
            ! local information is needed:
            CALL resize_lonlat_state(lonlat_grid_list(i)%intp(jg),     &
              &                      (nthis_local_pts - 1)/nproma + 1, &
              &                      nthis_local_pts)
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
      CHARACTER(*), PARAMETER :: routine = &
        &  TRIM("mo_intp_lonlat:compute_lonlat_area_weights")
      TYPE(t_cf_var)       :: cf_desc
      TYPE(t_grib2_var)    :: grib2_desc
      INTEGER              :: var_shape(3), nblks_lonlat, i, jg, ierrstat, &
        &                     i_lat, idx_glb, jc, jb, j
      TYPE (t_lon_lat_grid), POINTER :: grid
      TYPE (t_lon_lat_intp), POINTER :: ptr_int_lonlat
!DR      REAL(wp),              POINTER :: area_weights(:), p_dummy(:,:,:)
!DR !!! Using the POINTER attribute for area_weights(:) mysteriously leads
!DR !!! to a bus error on NEC SX9 (tested with compiler revision 450). However, 
!DR !!! using the ALLOCATABLE attribute, instead, works.
      REAL(wp), ALLOCATABLE          :: area_weights(:)
      REAL(wp),              POINTER :: p_dummy(:,:,:)
      TYPE(t_list_element),  POINTER :: new_element

      ! Add area weights
      DO i=1, n_lonlat_grids
        DO jg=1,n_dom
          IF (lonlat_grid_list(i)%l_dom(jg)) THEN
            
            grid           => lonlat_grid_list(i)%grid
            ptr_int_lonlat => lonlat_grid_list(i)%intp(jg)

            nblks_lonlat   =  (ptr_int_lonlat%nthis_local_pts - 1)/nproma + 1
            var_shape = (/ nproma, 1, nblks_lonlat /)
            cf_desc    = t_cf_var('gw', '1', 'area weights', DATATYPE_FLT32)
            grib2_desc = t_grib2_var(0, 0, 0, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL)

            ALLOCATE(area_weights(grid%lat_dim), STAT=ierrstat)
            IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')
            
            CALL add_var( p_nh_state(jg)%diag_list,                             &
              &           "gw", p_dummy,                                        &
              &           GRID_REGULAR_LONLAT, ZA_SURFACE, cf_desc, grib2_desc, &
              &           ldims=var_shape, lrestart=.FALSE.,                    &
              &           loutput=.TRUE., new_element=new_element,              &
              &           hor_interp=create_hor_interp_metadata(                &
              &             hor_intp_type=HINTP_TYPE_NONE ),                    &
              &           isteptype=TSTEP_CONSTANT,                             &
              &           cdiTimeID=TIME_CONSTANT )                ! deprecated
            ! link this new variable to the lon-lat grid:
            new_element%field%info%hor_interp%lonlat_id = i
            ! compute area weights:
            CALL latlon_compute_area_weights(grid, area_weights)
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

      CHARACTER(*), PARAMETER :: routine = TRIM("mo_intp_state:allocate_int_state_lonlat")
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

      ALLOCATE(ptr_int_lonlat%rdist(2, nproma, nblks_lonlat), stat=ist)
      IF (ist /= SUCCESS)  CALL finish (routine, 'allocation for working arrays failed')

      ptr_int_lonlat%rdist             = 0._wp
      ptr_int_lonlat%rbf_vec_idx       = 0
      ptr_int_lonlat%rbf_vec_blk       = 0
      ptr_int_lonlat%rbf_vec_stencil   = 0
      ptr_int_lonlat%rbf_vec_coeff     = 0._wp

      IF (l_intp_c2l) THEN
        ALLOCATE ( ptr_int_lonlat%rbf_c2l_coeff(rbf_dim_c2l, nproma, nblks_lonlat),         &
          &  ptr_int_lonlat%rbf_c2l_idx(rbf_dim_c2l, nproma, nblks_c),                      &
          &  ptr_int_lonlat%rbf_c2l_blk(rbf_dim_c2l, nproma, nblks_c),                      &
          &  ptr_int_lonlat%rbf_c2lr_idx(rbf_dim_c2l, nproma, nblks_lonlat),                      &
          &  ptr_int_lonlat%rbf_c2lr_blk(rbf_dim_c2l, nproma, nblks_lonlat),                      &
          &  ptr_int_lonlat%rbf_c2l_stencil(nproma, nblks_c),                               &
          STAT=ist )
        IF (ist /= SUCCESS) &
          CALL finish (routine, 'allocation for rbf lon-lat coeffs failed')

        ptr_int_lonlat%rbf_c2l_coeff   = 0._wp
        ptr_int_lonlat%rbf_c2l_idx     = 0
        ptr_int_lonlat%rbf_c2l_blk     = 0
        ptr_int_lonlat%rbf_c2lr_idx    = 0
        ptr_int_lonlat%rbf_c2lr_blk    = 0
        ptr_int_lonlat%rbf_c2l_stencil = 0
      ELSE
        ALLOCATE ( &
          &  ptr_int_lonlat%rbf_c2grad_coeff(rbf_c2grad_dim, 2, nproma, nblks_lonlat), &
          &  ptr_int_lonlat%rbf_c2grad_idx(rbf_c2grad_dim, nproma, nblks_lonlat),      &
          &  ptr_int_lonlat%rbf_c2grad_blk(rbf_c2grad_dim, nproma, nblks_lonlat),      &
          &  ptr_int_lonlat%cell_vert_dist(nproma, 3, 2, nblks_lonlat),                &
          &  STAT=ist )
        IF (ist /= SUCCESS) THEN
          CALL finish (routine, 'allocation for rbf lon-lat coeffs failed')
        ENDIF
        ptr_int_lonlat%rbf_c2grad_coeff  = 0._wp
        ptr_int_lonlat%rbf_c2grad_idx    = 0
        ptr_int_lonlat%rbf_c2grad_blk    = 0
        ptr_int_lonlat%cell_vert_dist    = 0._wp
      END IF

    END SUBROUTINE allocate_int_state_lonlat_grid

    
    !-------------------------------------------------------------------------
    !
    !>
    !! Deallocation of components of a single 2d lon-lat interpolation state.
    !!
    SUBROUTINE deallocate_int_state_lonlat( ptr_int_lonlat )

      TYPE (t_lon_lat_intp), TARGET, INTENT(INOUT) :: ptr_int_lonlat
      CHARACTER(*), PARAMETER :: routine = TRIM("mo_intp_state:deallocate_int_state_lonlat")
      INTEGER :: ist

      !-----------------------------------------------------------------------

      DEALLOCATE (ptr_int_lonlat%rbf_vec_coeff,           &
        &         ptr_int_lonlat%rbf_vec_idx,             &
        &         ptr_int_lonlat%rbf_vec_blk,             &
        &         ptr_int_lonlat%rbf_vec_stencil,         &
        &         ptr_int_lonlat%rdist,                   &
        &         ptr_int_lonlat%tri_idx,                 &
        &         ptr_int_lonlat%global_idx,              &
        &         STAT=ist )
      IF (ist /= SUCCESS) THEN
        CALL finish (routine, 'deallocation for lon-lat coefficients failed')
      ENDIF

      IF (l_intp_c2l) THEN
        DEALLOCATE ( ptr_int_lonlat%rbf_c2l_coeff,   &
          &          ptr_int_lonlat%rbf_c2l_idx,     &
          &          ptr_int_lonlat%rbf_c2l_blk,     &
          &          ptr_int_lonlat%rbf_c2lr_idx,    &
          &          ptr_int_lonlat%rbf_c2lr_blk,    &
          &          ptr_int_lonlat%rbf_c2l_stencil, &
          &          STAT=ist )
        IF (ist /= SUCCESS) &
          CALL finish (routine, 'deallocation for rbf lon-lat coeffs failed')
      ELSE
        DEALLOCATE (ptr_int_lonlat%rbf_c2grad_coeff,        &
          &         ptr_int_lonlat%rbf_c2grad_idx,          &
          &         ptr_int_lonlat%rbf_c2grad_blk,          &
          &         ptr_int_lonlat%cell_vert_dist,          &
          &         STAT=ist )
        IF (ist /= SUCCESS) THEN
          CALL finish (routine, 'deallocation for lon-lat coefficients failed')
        ENDIF
      END IF

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

      CHARACTER(*), PARAMETER :: routine = TRIM("mo_intp_lonlat:rbf_c2l_index")
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
      INTEGER :: nblks_c
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

        ! values for the blocking
        nblks_c  = ptr_patch%nblks_int_c

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

            IF(.NOT. ptr_patch%cells%owner_mask(jc,jb)) CYCLE

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
    !      Initial implementation  by  F.Prill, DWD (2011-08)
    !
    SUBROUTINE rbf_vec_compute_coeff_lonlat( ptr_patch, ptr_int, ptr_int_lonlat, lon_lat_points, &
      &                                      nblks_lonlat, npromz_lonlat )

      ! Input parameters
      TYPE(t_patch),                 INTENT(IN)    :: ptr_patch
      TYPE (t_int_state),    TARGET, INTENT(INOUT) :: ptr_int
      TYPE (t_lon_lat_intp),         INTENT(INOUT) :: ptr_int_lonlat
      REAL(gk),                      INTENT(IN)    :: lon_lat_points(:,:,:)
      INTEGER,                       INTENT(IN)    :: nblks_lonlat, npromz_lonlat ! blocking info
      ! Local parameters
      CHARACTER(*), PARAMETER :: routine = TRIM("mo_intp_lonlat:rbf_vec_compute_coeff_lonlat")
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
        &                                 istencil(nproma),        & ! actual number of stencil points
        &                                 jg
      REAL(wp)                         :: checksum_u,checksum_v      ! to check if sum of interpolation coefficients is correct
      TYPE(t_geographical_coordinates) :: grid_point

      !--------------------------------------------------------------------

      CALL message(routine, '')
      IF (ptr_patch%n_patch_cells == 0) RETURN;

      jg = ptr_patch%id

!$OMP PARALLEL PRIVATE (z_rbfmat,z_diag,z_rbfval,z_rhs1,z_rhs2, ist,   &
!$OMP                   grid_point),                                   &
!$OMP          SHARED  (nproma, rbf_vec_dim_c, nblks_lonlat,           &
!$OMP                   npromz_lonlat, ptr_int_lonlat, ptr_patch,      &
!$OMP                   rbf_vec_kern_ll, jg, rbf_vec_scale_ll )

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
                z_rbfmat(jc,je1,je2) = z_nxprod * gaussi(z_dist,rbf_vec_scale_ll(MAX(jg,1)))
              ELSE IF (rbf_vec_kern_ll == 3) THEN
                z_rbfmat(jc,je1,je2) = z_nxprod * inv_multiq(z_dist,rbf_vec_scale_ll(MAX(jg,1)))
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
          grid_point%lon = REAL( lon_lat_points(jc, jb,1), wp)
          grid_point%lat = REAL( lon_lat_points(jc, jb,2), wp)
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
              z_rbfval(jc,je2) = gaussi(z_dist,rbf_vec_scale_ll(MAX(jg,1)))
            ELSE IF (rbf_vec_kern_ll == 3) THEN
              z_rbfval(jc,je2) = inv_multiq(z_dist,rbf_vec_scale_ll(MAX(jg,1)))
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


    !> Computed combined coefficient/finite difference matrix for lon-lat interpolation.
    !! This routine computes the coefficients needed for reconstructing the gradient
    !! of a cell-based variable at the cell center. The operations performed here
    !! combine taking the centered-difference gradient (like in grad_fd_norm) at
    !! the 9 edges entering into the usual 9-point RBF vector reconstruction at
    !! cell centers, followed by applying the RBF reconstruction.
    !!
    !! @note Use of a Gradient Limiter:
    !!       --------------------------
    !!       In order to avoid over- and undershoots we suggest to
    !!       employ a gradient limiter (not yet implemented). Because
    !!       the gradient must be limited wrt. all lon-lat points
    !!       falling into a given cell, and because the number of
    !!       these "target points" is not known a priori, it would be
    !!       efficient to limit the gradient wrt. the cell vertices.
    !!
    !! @par Revision History
    !! based on
    !!   mo_intp_rbf_coeffs::rbf_compute_coeff_c2grad
    !! developed and tested by Guenther Zaengl, DWD (2009-12-15)
    !! Restructuring F. Prill, DWD (2011-08-18)
    !!
    SUBROUTINE rbf_compute_coeff_c2grad_lonlat (ptr_patch, ptr_int,          &
      &                                         nblks_lonlat, npromz_lonlat, &
      &                                         ptr_int_lonlat)

      TYPE(t_patch),     INTENT(in)    :: ptr_patch                   ! patch on which computation is performed
      TYPE (t_int_state),    TARGET, INTENT(INOUT) :: ptr_int
      INTEGER,           INTENT(IN)    :: nblks_lonlat, npromz_lonlat ! lon-lat grid blocking info
      TYPE (t_lon_lat_intp), TARGET, INTENT(inout) :: ptr_int_lonlat  ! interpolation state

      ! local variables
      INTEGER  :: je, jcc,                                   &
        &         jc, jb,                                    &  ! integer over lon-lat points
        &         i_startidx, i_endidx,                      &
        &         ile, ibe, ilc, ibc, ilcc, ibcc, jc_c, jb_c
      REAL(wp) :: inv_length, coeff(2),                      &
        &         aux_coeff(nproma,rbf_c2grad_dim,2)
      REAL(wp), DIMENSION(:,:,:,:), POINTER :: ptr_coeff

      !--------------------------------------------------------------------

      ptr_coeff => ptr_int%rbf_vec_coeff_c

      ! loop through all patch cells (and blocks)

!$OMP PARALLEL 
!$OMP DO PRIVATE(jb,je,jcc,jc,i_startidx,i_endidx,ile,ibe,ilc,ibc, &
!$OMP            aux_coeff, inv_length, coeff, ilcc, ibcc, jc_c, jb_c)
      DO jb = 1,nblks_lonlat

        i_startidx = 1
        i_endidx   = nproma
        IF (jb == nblks_lonlat) i_endidx = npromz_lonlat

        aux_coeff(:,:,:) = 0._wp

        DO je = 1, rbf_vec_dim_c
          DO jcc = 1, rbf_c2grad_dim
!CDIR NODEP
            DO jc = i_startidx, i_endidx

              ile = ptr_int_lonlat%rbf_vec_idx(je,jc,jb)
              ibe = ptr_int_lonlat%rbf_vec_blk(je,jc,jb)

              ilcc = ptr_int_lonlat%rbf_c2grad_idx(jcc,jc,jb)
              ibcc = ptr_int_lonlat%rbf_c2grad_blk(jcc,jc,jb)

              inv_length = ptr_patch%edges%inv_dual_edge_length(ile,ibe)*grid_sphere_radius

              jc_c = ptr_int_lonlat%tri_idx(1,jc, jb)
              jb_c = ptr_int_lonlat%tri_idx(2,jc, jb)
              coeff(1:2) = ptr_coeff(je,1:2,jc_c,jb_c) * inv_length

              ! every edge (ile,ibe) of the stencil for cell (jc,jb)
              ! contributes to two entries of the combined coeff/finite
              ! difference matrix block "aux_coeff" (with the sign "-/+1")

              ilc = ptr_patch%edges%cell_idx(ile,ibe,1)
              ibc = ptr_patch%edges%cell_blk(ile,ibe,1)

              IF (ilcc == ilc .AND. ibcc == ibc) THEN
                aux_coeff(jc,jcc,1:2) = aux_coeff(jc,jcc,1:2) - coeff(1:2)
              ELSE
                ilc = ptr_patch%edges%cell_idx(ile,ibe,2)
                ibc = ptr_patch%edges%cell_blk(ile,ibe,2)
                IF (ilcc == ilc .AND. ibcc == ibc) THEN
                  aux_coeff(jc,jcc,1:2) = aux_coeff(jc,jcc,1:2) + coeff(1:2)
                END IF
              END IF
            END DO ! jc
          END DO ! jcc
        ENDDO !je

        FORALL (jcc = 1:rbf_c2grad_dim)
          ptr_int_lonlat%rbf_c2grad_coeff(jcc,1,i_startidx:i_endidx,jb) = &
            & aux_coeff(i_startidx:i_endidx,jcc,1)
          ptr_int_lonlat%rbf_c2grad_coeff(jcc,2,i_startidx:i_endidx,jb) = &
            & aux_coeff(i_startidx:i_endidx,jcc,2)
        END FORALL

      END DO !block loop
!$OMP END DO
!$OMP END PARALLEL

    END SUBROUTINE rbf_compute_coeff_c2grad_lonlat


    !-------------------------------------------------------------------------
    !> Computes the RBF coefficients for lon-lat grid points.
    !
    !  We directly interpolate from cell centers to lon-lat points,
    !  otherwise we do gradient interpolation and reconstruction.
    !  -------------------------------------------------------------------------
    !
    ! @par Revision History
    !      Initial implementation  by  F.Prill, DWD (2012-06-13)
    !
    SUBROUTINE rbf_compute_coeff_c2l( ptr_patch, ptr_int_lonlat, lon_lat_points, &
      &                               nblks_lonlat, npromz_lonlat )

      ! Input parameters
      TYPE(t_patch),                 INTENT(IN)    :: ptr_patch
      TYPE (t_lon_lat_intp),         INTENT(INOUT) :: ptr_int_lonlat
      REAL(gk),                      INTENT(IN)    :: lon_lat_points(:,:,:)
      INTEGER,                       INTENT(IN)    :: nblks_lonlat, npromz_lonlat ! blocking info
      ! Local parameters
      CHARACTER(*), PARAMETER :: routine = TRIM("mo_intp_lonlat:rbf_compute_coeff_c2l")
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
        &                                 jg, jb_cell, jc_cell
      REAL(wp)                         :: checksum                   ! to check if sum of interpolation coefficients is correct
      TYPE(t_geographical_coordinates) :: grid_point

      !--------------------------------------------------------------------

      CALL message(routine, '')
      IF (ptr_patch%n_patch_cells == 0) RETURN;

      jg = ptr_patch%id

!OMP PARALLEL PRIVATE (z_rbfmat,z_diag,z_rbfval, ist, grid_point),    &
!OMP          SHARED  (nproma, rbf_dim_c2l, nblks_lonlat,             &
!OMP                   npromz_lonlat, ptr_int_lonlat, ptr_patch,      &
!OMP                   rbf_vec_kern_ll, jg, rbf_vec_scale_ll )

      ALLOCATE( z_rbfmat(nproma,rbf_dim_c2l,rbf_dim_c2l),    &
        &       z_diag(nproma,rbf_dim_c2l),                  &
        &       z_rbfval(nproma,rbf_dim_c2l),                &
        &       STAT=ist )
      IF (ist /= SUCCESS) &
        & CALL finish (routine, 'allocation for working arrays failed')

!OMP DO PRIVATE (jb,jc,i_startidx,i_endidx,je1,je2,istencil,       &
!OMP             ist,ilc1,ibc1, ilc2,ibc2,cc_1, cc_2,cc_c,z_dist,  &
!OMP             cc_center, checksum, jb_cell, jc_cell, stencil )
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
                z_rbfmat(jc,je1,je2) = gaussi(z_dist,rbf_vec_scale_ll(MAX(jg,1)))
              ELSE IF (rbf_vec_kern_ll == 3) THEN
                z_rbfmat(jc,je1,je2) = inv_multiq(z_dist,rbf_vec_scale_ll(MAX(jg,1)))
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

          grid_point%lon = REAL( lon_lat_points(jc, jb,1), wp)
          grid_point%lat = REAL( lon_lat_points(jc, jb,2), wp)
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
              z_rbfval(jc,je2) = gaussi(z_dist,rbf_vec_scale_ll(MAX(jg,1)))
            ELSE IF (rbf_vec_kern_ll == 3) THEN
              z_rbfval(jc,je2) = inv_multiq(z_dist,rbf_vec_scale_ll(MAX(jg,1)))
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
!OMP END DO

      DEALLOCATE( z_rbfmat, z_diag, z_rbfval, STAT=ist )
      IF (ist /= SUCCESS) &
        CALL finish (routine, 'deallocation for working arrays failed')
!OMP END PARALLEL
      
    END SUBROUTINE rbf_compute_coeff_c2l



    !-------------------------------------------------------------------------
    !> Setup routine for RBF reconstruction at lon-lat grid points for an arbitrary grid.
    ! 
    ! @par Revision History
    !      Initial implementation  by  F.Prill, DWD (2011-08)
    !      Changed for abritrary grids by Rainer Johanni (2011-11)
    !
    ! Note:   Two variables are created globally: "tri_idx" and "global_idx"
    !         This is necessary for the splitting of the lon-lat grid
    !         points over the PEs and for creating the corresponding
    !         communication patterns. However, after this initial phase,
    !         these data structures are resized to the local sizes.
    !
    ! @todo   For global grids, one should add only one point of the
    !         first and last latitude rows to the interpolation point
    !         list. Then, during interpolation, the interpolation
    !         result should be copied on to the other points with
    !         +-90deg.
    !
    SUBROUTINE rbf_setup_interpol_lonlat_grid(grid, ptr_patch, ptr_int_lonlat, ptr_int)

      TYPE (t_lon_lat_grid), INTENT(INOUT)         :: grid
      ! data structure containing grid info:
      TYPE(t_patch), TARGET, INTENT(IN)            :: ptr_patch
      ! Indices of source points and interpolation coefficients
      TYPE (t_lon_lat_intp), TARGET, INTENT(INOUT) :: ptr_int_lonlat
      TYPE (t_int_state),    TARGET, INTENT(INOUT) :: ptr_int
      ! Local Parameters:
      CHARACTER(*), PARAMETER :: routine = TRIM("mo_intp_rbf_coeffs:rbf_setup_interpol_lonlat")
      REAL(wp), ALLOCATABLE            :: rotated_pts(:,:,:)
      REAL(gk), ALLOCATABLE            :: in_points(:,:,:)
      REAL(gk), ALLOCATABLE            :: min_dist(:,:)       ! minimal distance
      INTEGER                          :: jb, jc, i_startidx, i_endidx,         &
        &                                 nblks_lonlat, npromz_lonlat, errstat, &
        &                                 jb_lonlat, jc_lonlat,                 &
        &                                 rl_start, rl_end, i_nchdom, i_startblk
      LOGICAL                          :: l_cutoff_local_domains
      LOGICAL, ALLOCATABLE             :: l_cutoff(:,:)
      TYPE(t_geographical_coordinates) :: cell_center, lonlat_pt
      TYPE(t_cartesian_coordinates)    :: p1, p2
      REAL(wp)                         :: point(2), z_norm, z_nx1(3), z_nx2(3), &
        &                                 max_dist

      !-----------------------------------------------------------------------

      ! Flag: .TRUE., if we want to erase values outside local domains
      l_cutoff_local_domains = .TRUE.

      IF (dbg_level > 1) THEN
        WRITE(message_text,*) "SETUP : rbf_interpol_lonlat"
        CALL message(routine, message_text)
      END IF
      IF (ptr_patch%cell_type == 6) THEN
        CALL finish(routine, "Lon-lat interpolation not yet implemented for cell_type == 6!")
      END IF

      nblks_lonlat  = grid%nblks
      npromz_lonlat = grid%npromz

      !-- compute grid points of rotated lon/lat grid

      ALLOCATE(in_points(nproma, nblks_lonlat, 2),            &
        &      rotated_pts(grid%lon_dim, grid%lat_dim, 2),    &
        &      min_dist(nproma, nblks_lonlat),                &
        &      stat=errstat)
      IF (errstat /= SUCCESS) THEN
        CALL finish (routine, 'allocation for working arrays failed')
      ENDIF

      IF (dbg_level > 1) &
        CALL message(routine, "rotate lon-lat grid points")
      CALL rotate_latlon_grid( grid, rotated_pts )

      in_points(:,:,1) = RESHAPE(REAL(rotated_pts(:,:,1), gk), &
        &                        (/ nproma, nblks_lonlat /), (/ 0._gk /) )
      in_points(:,:,2) = RESHAPE(REAL(rotated_pts(:,:,2), gk), &
        &                        (/ nproma, nblks_lonlat /), (/ 0._gk /) )

      !-- proximity query, build a list of cell indices that contain the
      !   lon/lat grid points

      ! build GNAT data structure
      IF (dbg_level > 1) &
        CALL message(routine, "build GNAT data structure")
      CALL gnat_init_grid(ptr_patch)

      ! allocate global arrays for distributed computation:
      nblks_lonlat  = (grid%total_dim - 1)/nproma + 1
      npromz_lonlat = grid%total_dim - (nblks_lonlat-1)*nproma
      ALLOCATE(ptr_int_lonlat%tri_idx(2, nproma, nblks_lonlat),   &
        &      ptr_int_lonlat%global_idx(grid%total_dim), STAT=errstat )
      IF (errstat /= SUCCESS) &
        CALL finish (routine, 'allocation for lon-lat point distribution failed')
      ptr_int_lonlat%tri_idx(:,:,:) = 0

      ! Perform query. Note that for distributed patches we receive a
      ! local list of "in_points" actually located on this portion of the
      ! domain.
      IF (dbg_level > 1) &
        CALL message(routine, "proximity query")
      CALL gnat_query_containing_triangles(ptr_patch, gnat_tree, in_points(:,:,:),       &
        &                   nproma, nblks_lonlat, npromz_lonlat, grid_sphere_radius,     &
        &                   p_test_run, ptr_int_lonlat%tri_idx(:,:,:), min_dist(:,:))
      CALL gnat_merge_distributed_queries(ptr_patch, grid%total_dim, nproma, grid%nblks, &
        &                 min_dist, ptr_int_lonlat%tri_idx(:,:,:), in_points(:,:,:),     &
        &                 ptr_int_lonlat%global_idx(:), ptr_int_lonlat%nthis_local_pts)
        
      ! set local values for "nblks" and "npromz"
      nblks_lonlat  = (ptr_int_lonlat%nthis_local_pts - 1)/nproma + 1
      npromz_lonlat = ptr_int_lonlat%nthis_local_pts - (nblks_lonlat-1)*nproma
      ! clean up
      CALL gnat_destroy()

      ! After the lon-lat points have been distributed, we know
      ! exactly for each PE the number of points
      CALL allocate_int_state_lonlat_grid( ptr_patch%nblks_c, nblks_lonlat, ptr_int_lonlat)

      !-- edge indices of stencil
      ! are available through previous calls to SR rbf_vec_index_cell()
      ! and rbf_c2grad_index()

      !-- copy neighbor indices from standard RBF interpolation
      IF (dbg_level > 1) &
        CALL message(routine, "copy neighbor indices from standard RBF interpolation")

!$OMP PARALLEL SHARED(nblks_lonlat, nproma, npromz_lonlat,   &
!$OMP                 ptr_int_lonlat, ptr_int, l_intp_c2l),  &
!$OMP          DEFAULT (NONE)
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc), SCHEDULE(runtime)
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

        IF (.NOT.l_intp_c2l) THEN
          DO jc_lonlat = i_startidx, i_endidx

            jc = ptr_int_lonlat%tri_idx(1,jc_lonlat, jb_lonlat)
            jb = ptr_int_lonlat%tri_idx(2,jc_lonlat, jb_lonlat)

            ptr_int_lonlat%rbf_c2grad_idx(:,jc_lonlat,jb_lonlat) = ptr_int%rbf_c2grad_idx(:,jc,jb)
            ptr_int_lonlat%rbf_c2grad_blk(:,jc_lonlat,jb_lonlat) = ptr_int%rbf_c2grad_blk(:,jc,jb)
            ptr_int_lonlat%cell_vert_dist(jc_lonlat,1:3,1:2,jb_lonlat) = ptr_int%cell_vert_dist(jc,1:3,1:2,jb)
          ENDDO
        ENDIF

      ENDDO
!$OMP END DO
!$OMP END PARALLEL

      !-- compute interpolation coefficients
      IF (dbg_level > 1) &
        CALL message(routine, "compute lon-lat interpolation coefficients")
      CALL rbf_vec_compute_coeff_lonlat( ptr_patch, ptr_int, ptr_int_lonlat,  &
        &                                in_points, nblks_lonlat, npromz_lonlat )

      
      IF (l_intp_c2l) THEN
        CALL rbf_c2l_index( ptr_patch, ptr_int, ptr_int_lonlat )
        CALL rbf_compute_coeff_c2l( ptr_patch, ptr_int_lonlat,  &
          &                         in_points, nblks_lonlat, npromz_lonlat )

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

            ptr_int_lonlat%rbf_c2lr_idx(:,jc_lonlat,jb_lonlat) = ptr_int_lonlat%rbf_c2l_idx(:,jc,jb)   
            ptr_int_lonlat%rbf_c2lr_blk(:,jc_lonlat,jb_lonlat) = ptr_int_lonlat%rbf_c2l_blk(:,jc,jb)
          ENDDO

        ENDDO
!$OMP END DO
!$OMP END PARALLEL

      END IF

      ! compute distances (x0i - xc)
      DO jb=1,nblks_lonlat
        i_startidx = 1
        i_endidx   = nproma
        IF (jb == nblks_lonlat) i_endidx = npromz_lonlat

        DO jc=i_startidx,i_endidx
          point(:)    = REAL(in_points(jc,jb,:), wp)
          lonlat_pt%lon = point(1)
          lonlat_pt%lat = point(2)
          cell_center = ptr_patch%cells%center(ptr_int_lonlat%tri_idx(1,jc,jb),      &
            &                                  ptr_int_lonlat%tri_idx(2,jc,jb))

          ! convert to Cartesian coordinate system:
          p1 = gc2cc(lonlat_pt)
          p2 = gc2cc(cell_center)
          p1%x(:) = p1%x(:) - p2%x(:)

          ! Zonal component
          CALL gvec2cvec(1._wp,0._wp, point(1), point(2), &
            &            z_nx1(1),z_nx1(2),z_nx1(3))
          z_norm = SQRT( DOT_PRODUCT(z_nx1(:),z_nx1(:)) )
          z_nx1(:)  = 1._wp/z_norm * z_nx1(:)
          ! Meridional component
          CALL gvec2cvec(0._wp, 1._wp, point(1), point(2), &
            &            z_nx2(1),z_nx2(2),z_nx2(3))
          z_norm = SQRT( DOT_PRODUCT(z_nx2(:),z_nx2(:)) )
          z_nx2(:)  = 1._wp/z_norm * z_nx2(:)
          ! projection, scale with earth radius
          ptr_int_lonlat%rdist(1, jc, jb) = DOT_PRODUCT(p1%x(:), z_nx1)
          ptr_int_lonlat%rdist(2, jc, jb) = DOT_PRODUCT(p1%x(:), z_nx2)
        END DO
      END DO

      IF (.NOT.l_intp_c2l) THEN
        CALL rbf_compute_coeff_c2grad_lonlat (ptr_patch, ptr_int, nblks_lonlat, &
          &                                   npromz_lonlat, ptr_int_lonlat)
      ENDIF

      IF ( l_cutoff_local_domains ) THEN
        ! make a sensible guess for maximum distance (just taking a
        ! "randomly chosen" edge length)
        rl_start = 2
        rl_end   = min_rlcell_int
        i_nchdom   = MAX(1,ptr_patch%n_childdom)
        i_startblk = ptr_patch%cells%start_blk(rl_start,1)
        CALL get_indices_e(ptr_patch, i_startblk,  &
          &                i_startblk, i_startblk, &
          &                i_startidx, i_endidx, &
          &                rl_start, rl_end)
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

        ! build a logical array, .true. where lon-lat point is outside of
        ! the current patch:
        ALLOCATE(l_cutoff(nproma, nblks_lonlat), stat=errstat)
        IF (errstat /= SUCCESS) CALL finish (routine, 'allocation for working arrays failed')
        l_cutoff(:,:) = (MAX(ABS(ptr_int_lonlat%rdist(1,:,:)), &
          &                  ABS(ptr_int_lonlat%rdist(2,:,:))) >= max_dist)
        
        WHERE (l_cutoff(:,:))
          ptr_int_lonlat%rdist(1,:,:) = 0._wp
          ptr_int_lonlat%rdist(2,:,:) = 0._wp
        END WHERE
        
        ! clean up
        DEALLOCATE(l_cutoff, stat=errstat)
        IF (errstat /= SUCCESS)  CALL finish (routine, 'DEALLOCATE of working arrays failed')
      END IF
      
      DEALLOCATE(in_points, rotated_pts, min_dist, stat=errstat)
      IF (errstat /= SUCCESS) &
        CALL finish (routine, 'DEALLOCATE of working arrays failed')

    END SUBROUTINE rbf_setup_interpol_lonlat_grid


    !===============================================================
    ! APPLY LON-LAT INTERPOLATION

    !-------------------------------------------------------------------------
    !> Performs vector RBF reconstruction at lon-lat grid points.
    !
    ! This routine is based on mo_intp_rbf::rbf_vec_interpol_cell()
    ! 
    ! @par Revision History
    !      Initial implementation  by  F.Prill, DWD (2011-08)
    !
    SUBROUTINE rbf_vec_interpol_lonlat_nl( p_vn_in, ptr_int, &
      &                                    grad_x, grad_y,              &
      &                                    nblks_lonlat, npromz_lonlat, &
      &                                    opt_slev, opt_elev)

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

    END SUBROUTINE rbf_vec_interpol_lonlat_nl


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
    !      Initial implementation  by  F.Prill, DWD (2011-08)
    !      based on "rbf_interpol_c2grad"
    !
    SUBROUTINE rbf_interpol_c2grad_lonlat_nl( p_cell_in, ptr_int,            &
      &                                       grad_x, grad_y,                &
      &                                       nblks_lonlat, npromz_lonlat,   &
      &                                       opt_slev, opt_elev)

      ! !INPUT PARAMETERS
      !
      ! input cell-based variable for which gradient at cell center is computed
      REAL(wp),                  INTENT(IN) :: p_cell_in(:,:,:) ! dim: (nproma,nlev,nblks_c)
      ! Indices of source points and interpolation coefficients
      TYPE (t_lon_lat_intp), TARGET, INTENT(IN) :: ptr_int

      ! reconstructed zonal (x) component of gradient vector
      REAL(wp),INTENT(INOUT) :: grad_x(:,:,:) ! dim: (nproma,nlev,nblks_lonlat)
      ! reconstructed zonal (x) component of gradient vector
      REAL(wp),INTENT(INOUT) :: grad_y(:,:,:) ! dim: (nproma,nlev,nblks_lonlat)

      ! lon-lat grid blocking info
      INTEGER,                   INTENT(IN) :: nblks_lonlat, npromz_lonlat
      ! optional vertical start/end level
      INTEGER,                   INTENT(IN), OPTIONAL :: opt_slev, opt_elev

      ! Local variables
      INTEGER :: slev, elev,               &  ! vertical start and end level
        &          jc, jb, jk,             &  ! integer over lon-lat points, levels
        &          i_startidx, i_endidx       ! start/end index

      ! maximum and minimum difference between external stencil points and local point
      REAL(wp), DIMENSION(nproma,UBOUND(p_cell_in,2)) :: maxdif, mindif

      ! maximum and minimum extrapolation increments and reduction factor for gradient
      REAL(wp) :: maxextr, minextr, redfac

      INTEGER,  DIMENSION(:,:,:),   POINTER :: iidx, iblk
      REAL(wp), DIMENSION(:,:,:,:), POINTER :: ptr_coeff, ptr_c2v_dist

      !-----------------------------------------------------------------------

      slev = 1
      elev = UBOUND(p_cell_in,2)
      ! check optional arguments
      IF ( PRESENT(opt_slev) ) slev = opt_slev
      IF ( PRESENT(opt_elev) ) elev = opt_elev

      iidx => ptr_int%rbf_c2grad_idx
      iblk => ptr_int%rbf_c2grad_blk

      ptr_coeff => ptr_int%rbf_c2grad_coeff
      ptr_c2v_dist => ptr_int%cell_vert_dist

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk,jc,maxdif,mindif,maxextr,minextr,redfac), SCHEDULE(runtime)
      DO jb = 1,nblks_lonlat

        i_startidx = 1
        i_endidx   = nproma
        IF (jb == nblks_lonlat) i_endidx = npromz_lonlat

#ifdef __LOOP_EXCHANGE
        DO jc = i_startidx, i_endidx
          DO jk = slev, elev
#else
        DO jk = slev, elev
          DO jc = i_startidx, i_endidx
#endif

            grad_x(jc,jk,jb) = ptr_coeff(1, 1,jc,jb) * p_cell_in(iidx(1 ,jc,jb),jk,iblk(1 ,jc,jb)) &
              &              + ptr_coeff(2, 1,jc,jb) * p_cell_in(iidx(2 ,jc,jb),jk,iblk(2 ,jc,jb)) &
              &              + ptr_coeff(3, 1,jc,jb) * p_cell_in(iidx(3 ,jc,jb),jk,iblk(3 ,jc,jb)) &
              &              + ptr_coeff(4, 1,jc,jb) * p_cell_in(iidx(4 ,jc,jb),jk,iblk(4 ,jc,jb)) &
              &              + ptr_coeff(5, 1,jc,jb) * p_cell_in(iidx(5 ,jc,jb),jk,iblk(5 ,jc,jb)) &
              &              + ptr_coeff(6, 1,jc,jb) * p_cell_in(iidx(6 ,jc,jb),jk,iblk(6 ,jc,jb)) &
              &              + ptr_coeff(7, 1,jc,jb) * p_cell_in(iidx(7 ,jc,jb),jk,iblk(7 ,jc,jb)) &
              &              + ptr_coeff(8, 1,jc,jb) * p_cell_in(iidx(8 ,jc,jb),jk,iblk(8 ,jc,jb)) &
              &              + ptr_coeff(9, 1,jc,jb) * p_cell_in(iidx(9 ,jc,jb),jk,iblk(9 ,jc,jb)) &
              &              + ptr_coeff(10,1,jc,jb) * p_cell_in(iidx(10,jc,jb),jk,iblk(10,jc,jb))

            grad_y(jc,jk,jb) = ptr_coeff(1, 2,jc,jb) * p_cell_in(iidx(1 ,jc,jb),jk,iblk(1 ,jc,jb)) &
              &              + ptr_coeff(2, 2,jc,jb) * p_cell_in(iidx(2 ,jc,jb),jk,iblk(2 ,jc,jb)) &
              &              + ptr_coeff(3, 2,jc,jb) * p_cell_in(iidx(3 ,jc,jb),jk,iblk(3 ,jc,jb)) &
              &              + ptr_coeff(4, 2,jc,jb) * p_cell_in(iidx(4 ,jc,jb),jk,iblk(4 ,jc,jb)) &
              &              + ptr_coeff(5, 2,jc,jb) * p_cell_in(iidx(5 ,jc,jb),jk,iblk(5 ,jc,jb)) &
              &              + ptr_coeff(6, 2,jc,jb) * p_cell_in(iidx(6 ,jc,jb),jk,iblk(6 ,jc,jb)) &
              &              + ptr_coeff(7, 2,jc,jb) * p_cell_in(iidx(7 ,jc,jb),jk,iblk(7 ,jc,jb)) &
              &              + ptr_coeff(8, 2,jc,jb) * p_cell_in(iidx(8 ,jc,jb),jk,iblk(8 ,jc,jb)) &
              &              + ptr_coeff(9, 2,jc,jb) * p_cell_in(iidx(9 ,jc,jb),jk,iblk(9 ,jc,jb)) &
              &              + ptr_coeff(10,2,jc,jb) * p_cell_in(iidx(10,jc,jb),jk,iblk(10,jc,jb))


            IF (l_mono_c2l) THEN ! prepare input for gradient limiter
              maxdif(jc,jk) = MAX(p_cell_in(iidx(2 ,jc,jb),jk,iblk(2 ,jc,jb)), &
                                  p_cell_in(iidx(3 ,jc,jb),jk,iblk(3 ,jc,jb)), &
                                  p_cell_in(iidx(4 ,jc,jb),jk,iblk(4 ,jc,jb)), &
                                  p_cell_in(iidx(5 ,jc,jb),jk,iblk(5 ,jc,jb)), &
                                  p_cell_in(iidx(6 ,jc,jb),jk,iblk(6 ,jc,jb)), &
                                  p_cell_in(iidx(7 ,jc,jb),jk,iblk(7 ,jc,jb)), &
                                  p_cell_in(iidx(8 ,jc,jb),jk,iblk(8 ,jc,jb)), &
                                  p_cell_in(iidx(9 ,jc,jb),jk,iblk(9 ,jc,jb)), &
                                  p_cell_in(iidx(10,jc,jb),jk,iblk(10,jc,jb))) &
                                - p_cell_in(iidx(1 ,jc,jb),jk,iblk(1 ,jc,jb))

              mindif(jc,jk) = MIN(p_cell_in(iidx(2 ,jc,jb),jk,iblk(2 ,jc,jb)), &
                                  p_cell_in(iidx(3 ,jc,jb),jk,iblk(3 ,jc,jb)), &
                                  p_cell_in(iidx(4 ,jc,jb),jk,iblk(4 ,jc,jb)), &
                                  p_cell_in(iidx(5 ,jc,jb),jk,iblk(5 ,jc,jb)), &
                                  p_cell_in(iidx(6 ,jc,jb),jk,iblk(6 ,jc,jb)), &
                                  p_cell_in(iidx(7 ,jc,jb),jk,iblk(7 ,jc,jb)), &
                                  p_cell_in(iidx(8 ,jc,jb),jk,iblk(8 ,jc,jb)), &
                                  p_cell_in(iidx(9 ,jc,jb),jk,iblk(9 ,jc,jb)), &
                                  p_cell_in(iidx(10,jc,jb),jk,iblk(10,jc,jb))) &
                                - p_cell_in(iidx(1 ,jc,jb),jk,iblk(1 ,jc,jb))
            ENDIF

          ENDDO
        ENDDO

        IF (l_mono_c2l) THEN ! apply gradient limiter in order to avoid overshoots

          DO jk = slev, elev
            DO jc = i_startidx, i_endidx

              maxextr = MAX(grad_x(jc,jk,jb)*ptr_c2v_dist(jc,1,1,jb)+grad_y(jc,jk,jb)*ptr_c2v_dist(jc,1,2,jb), &
                            grad_x(jc,jk,jb)*ptr_c2v_dist(jc,2,1,jb)+grad_y(jc,jk,jb)*ptr_c2v_dist(jc,2,2,jb), &
                            grad_x(jc,jk,jb)*ptr_c2v_dist(jc,3,1,jb)+grad_y(jc,jk,jb)*ptr_c2v_dist(jc,3,2,jb) )

              minextr = MIN(grad_x(jc,jk,jb)*ptr_c2v_dist(jc,1,1,jb)+grad_y(jc,jk,jb)*ptr_c2v_dist(jc,1,2,jb), &
                            grad_x(jc,jk,jb)*ptr_c2v_dist(jc,2,1,jb)+grad_y(jc,jk,jb)*ptr_c2v_dist(jc,2,2,jb), &
                            grad_x(jc,jk,jb)*ptr_c2v_dist(jc,3,1,jb)+grad_y(jc,jk,jb)*ptr_c2v_dist(jc,3,2,jb) )

              ! In the case of a local extremum, the gradient must be set to zero in order to avoid over-/undershoots
              IF (maxdif(jc,jk) <= 0._wp .OR. mindif(jc,jk) >= 0._wp) THEN
                redfac = 0._wp
              ! If the gradients are zero anyway, no limitation is needed
              ELSE IF (maxextr <= 0._wp .OR. minextr >= 0._wp) THEN
                redfac = 1._wp
              ELSE ! compute reduction factor for gradient
                redfac = MIN(1._wp, maxdif(jc,jk)/maxextr, mindif(jc,jk)/minextr)
              ENDIF

              grad_x(jc,jk,jb) = grad_x(jc,jk,jb)*redfac
              grad_y(jc,jk,jb) = grad_y(jc,jk,jb)*redfac

            ENDDO
          ENDDO

        ENDIF  

      ENDDO
!$OMP END DO
!$OMP END PARALLEL

    END SUBROUTINE rbf_interpol_c2grad_lonlat_nl


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
    !      Initial implementation  by  F.Prill, DWD (2011-08)
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
    
    !> Driver routine for RBF reconstruction of cell-based variables at
    !  lon-lat grid points.
    !
    ! @par Revision History
    !      Initial implementation  by  F.Prill, DWD (2011-08)
    !
    SUBROUTINE rbf_interpol_lonlat_nl( p_cell_in, ptr_int,            &
      &                                p_lonlat_out,                  &
      &                                nblks_lonlat, npromz_lonlat,   &
      &                                opt_slev, opt_elev)
      ! input cell-based variable for which gradient at cell center is computed
      REAL(wp),              INTENT(IN)           :: p_cell_in(:,:,:)    ! dim: (nproma,nlev,nblks_c)
      ! Indices of source points and interpolation coefficients
      TYPE (t_lon_lat_intp), TARGET, INTENT(IN)   :: ptr_int
      ! output lon-lat-based variable
      REAL(wp),              INTENT(INOUT)        :: p_lonlat_out(:,:,:) ! dim: (nproma,nlev,nblks_lonlat)
      ! lon-lat grid blocking info
      INTEGER,                   INTENT(IN)       :: nblks_lonlat, npromz_lonlat
      ! optional vertical start/end level
      INTEGER,                   INTENT(IN), OPTIONAL :: opt_slev, opt_elev

      ! Local Parameters:
      CHARACTER(*), PARAMETER :: routine = TRIM("mo_intp_rbf:rbf_interpol_lonlat")
      INTEGER  :: jb, jk, jc, i_startidx, i_endidx, slev, elev
      REAL(wp) :: grad_x(nproma, SIZE(p_cell_in,2), SIZE(p_lonlat_out,3)), &
        &         grad_y(nproma, SIZE(p_cell_in,2), SIZE(p_lonlat_out,3))

      !-----------------------------------------------------------------------

      slev = 1
      elev = UBOUND(p_cell_in,2)
      ! check optional arguments
      IF ( PRESENT(opt_slev) ) slev = opt_slev
      IF ( PRESENT(opt_elev) ) elev = opt_elev

      !-- apply interpolation coefficients
      IF (dbg_level > 1) THEN
        WRITE(message_text,*) "PE #", p_pe, ": apply interpolation coefficients"
        CALL message(routine, message_text)
      END IF

      IF (.NOT. l_intp_c2l) THEN

        CALL rbf_interpol_c2grad_lonlat_nl( p_cell_in(:,:,:), ptr_int,   &
          &                                 grad_x, grad_y,              &
          &                                 nblks_lonlat, npromz_lonlat, &
          &                                 slev, elev)

        !-- reconstruct scalar from gradient information
        IF (dbg_level > 1) THEN
          WRITE(message_text,*) "PE #", p_pe, ": reconstruct scalar from gradient information"
          CALL message(routine, message_text)
        END IF

        ! simple linear reconstruction
        ! given: zonal, meridional gradients d_1/2 in lon-lat grid points (x_0i, y_0i)
        !        and scalar values in cell centers (x_c, y_c)

        ! extrapolate: f(x_0i) = f(x_c) + (x_0i-x_c)*d_1 + (y_0i - y_c)*d_2

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk,jc), SCHEDULE(runtime)
        DO jb=1,nblks_lonlat
          i_startidx = 1
          i_endidx   = nproma
          IF (jb == nblks_lonlat) i_endidx = npromz_lonlat
!CDIR UNROLL=5
          DO jk=slev,elev
            DO jc=i_startidx,i_endidx

              p_lonlat_out(jc,jk,jb) = p_cell_in(ptr_int%tri_idx(1,jc,jb), jk, &
                &                                ptr_int%tri_idx(2,jc,jb))     &
                &           +  ptr_int%rdist(1,jc,jb) * grad_x(jc,jk,jb)       &
                &           +  ptr_int%rdist(2,jc,jb) * grad_y(jc,jk,jb)

            END DO
          END DO
        END DO
!$OMP END DO
!$OMP END PARALLEL

      ELSE

        CALL rbf_interpol_c2l( p_cell_in(:,:,:), ptr_int,   &
          &                    p_lonlat_out(:,:,:),         &
          &                    nblks_lonlat, npromz_lonlat, &
          &                    slev, elev)

      END IF

    END SUBROUTINE rbf_interpol_lonlat_nl

  END MODULE mo_intp_lonlat
