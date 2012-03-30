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
    USE mo_impl_constants,      ONLY: SUCCESS, min_rlcell_int, max_dom
    USE mo_model_domain,        ONLY: t_patch
    USE mo_run_config,          ONLY: ltimer
    USE mo_grid_config,         ONLY: n_dom
    USE mo_timer,               ONLY: timer_start, timer_stop, &
      &                               timers_level,            &
      &                               timer_lonlat_setup
    USE mo_math_utilities,      ONLY: gc2cc, gvec2cvec, solve_chol_v, choldec_v, &
      &                               arc_length_v, t_cartesian_coordinates,     &
      &                               t_geographical_coordinates
    USE mo_lonlat_grid,         ONLY: t_lon_lat_grid,                         &
      &                               compute_lonlat_blocking,                &
      &                               compute_lonlat_specs
    USE mo_parallel_config,     ONLY: nproma
    USE mo_loopindices,         ONLY: get_indices_c, get_indices_e, get_indices_v
    USE mo_intp_data_strc,      ONLY: t_int_state, t_lon_lat_intp, n_lonlat_grids, &
      &                               lonlat_grid_list, n_lonlat_grids, MAX_LONLAT_GRIDS
    USE mo_interpol_config,     ONLY: rbf_vec_dim_c, rbf_vec_dim_e, rbf_vec_dim_v,     &
      &                               rbf_c2grad_dim, rbf_vec_kern_c, rbf_vec_kern_e,  &
      &                               rbf_vec_kern_v, rbf_vec_scale_c, rbf_vec_scale_e,&
      &                               rbf_vec_scale_v
    USE mo_gnat_gridsearch,     ONLY: gnat_init_grid, gnat_destroy, gnat_tree,&
      &                               gnat_query_containing_triangles,        &
      &                               gnat_merge_distributed_queries, gk
    USE mo_math_utilities,      ONLY: rotate_latlon_grid
    USE mo_physical_constants,  ONLY: re
    USE mo_mpi,                 ONLY: p_gather_field, my_process_is_mpi_workroot, &
      &                               get_my_mpi_work_id, p_n_work
    USE mo_communication,       ONLY: idx_1d, blk_no, idx_no, &
      &                               setup_comm_pattern
    USE mo_mpi,                 ONLY: p_pe
    USE mo_lonlat_grid,         ONLY: t_lon_lat_grid
    USE mo_intp_state,          ONLY: allocate_int_state_lonlat_grid,    &
      &                               deallocate_int_state_lonlat

    IMPLICIT NONE

    !> level of output verbosity
    INTEGER, PARAMETER  :: dbg_level = 0

    PRIVATE

    PUBLIC :: init_lonlat_grid_list
    PUBLIC :: destroy_lonlat_grid_list
    PUBLIC :: compute_lonlat_intp_coeffs
    PUBLIC :: rbf_vec_compute_coeff_lonlat
    PUBLIC :: rbf_compute_coeff_c2grad_lonlat
    PUBLIC :: rbf_setup_interpol_lonlat_grid
    PUBLIC :: rbf_vec_interpol_lonlat
    PUBLIC :: rbf_interpol_c2grad_lonlat
    PUBLIC :: rbf_vec_interpol_lonlat_nl
    PUBLIC :: rbf_interpol_c2grad_lonlat_nl
    PUBLIC :: rbf_interpol_lonlat_nl

  CONTAINS

#include "intp_functions.f90.inc"

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
    !> Setup of lon-lat interpolation
    !
    SUBROUTINE compute_lonlat_intp_coeffs(p_patch, p_int_state)
      TYPE(t_patch),        INTENT(IN)    :: p_patch(:)
      TYPE(t_int_state),    INTENT(INOUT) :: p_int_state(:)
      ! local variables
      CHARACTER(*), PARAMETER :: routine = &
        &  TRIM("mo_intp_lonlat:compute_lonlat_intp_coeffs")
      INTEGER              :: i, j, jg, n_points, iloc, ist
      INTEGER, ALLOCATABLE :: owner(:), local_index(:)
      TYPE(t_lon_lat_grid), POINTER  :: grid
      
      IF (dbg_level > 5) CALL message(routine, "Enter")
      DO i=1, n_lonlat_grids

        ! compute some entries of lon-lat grid specification:
        CALL compute_lonlat_specs(lonlat_grid_list(i)%grid)
        CALL compute_lonlat_blocking(lonlat_grid_list(i)%grid, nproma)
        
        ! allocate and compute coefficients needed for lon-lat
        ! interpolation:
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
            ALLOCATE( owner(n_points), local_index(n_points), STAT=ist )
            IF (ist /= SUCCESS) &
              CALL finish (routine, 'allocation for working arrays failed')
            
            ! owner        : owner PE number of every point in the RECEIVER array
            ! local_index  : local index in the SENDER array.
            owner(:)       = lonlat_grid_list(i)%intp(jg)%owner(:)            
            local_index(:) = -1
            iloc = 0
            DO j=1,n_points
              IF (owner(j) == get_my_mpi_work_id()) THEN
                iloc = iloc + 1
                local_index(j) = iloc
              END IF
            END DO

            IF (my_process_is_mpi_workroot()) THEN
              CALL setup_comm_pattern(n_points, owner, local_index=local_index, &
                &                     p_pat=lonlat_grid_list(i)%p_pat(jg))
            ELSE
              ! We don't want to receive any data, i.e. the number of
              ! lon-lat points is 0 and owner/global index are
              ! dummies!
              owner(:) = -1
              CALL setup_comm_pattern(0, owner, local_index=local_index, &
                &                     p_pat=lonlat_grid_list(i)%p_pat(jg))
            END IF

            DEALLOCATE( owner, local_index, STAT=ist )
            IF (ist /= SUCCESS) &
              CALL finish (routine, 'deallocation for working arrays failed')
          END IF
        END DO ! jg

      END DO ! i
      IF (dbg_level > 5) CALL message(routine, "Done")

    END SUBROUTINE compute_lonlat_intp_coeffs


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
      CHARACTER(*), PARAMETER :: routine = TRIM("mo_interpolation:rbf_vec_compute_coeff_lonlat")
      REAL(wp)                         :: cc_e1(3), cc_e2(3), cc_c(nproma,3)  ! coordinates of edge midpoints
      TYPE(t_cartesian_coordinates)    :: cc_center                  ! coordinates of cell centers
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
        &                                 i_rcstartlev,            & ! refinement control start level
        &                                 istencil(nproma),        & ! actual number of stencil points
        &                                 jg
      REAL(wp)                         :: checksum_u,checksum_v      ! to check if sum of interpolation coefficients is correct
      TYPE(t_geographical_coordinates) :: grid_point

      !--------------------------------------------------------------------

      CALL message(routine, '')
      IF (ptr_patch%n_patch_cells == 0) RETURN;

      jg = ptr_patch%id

      i_rcstartlev = 2

!$OMP PARALLEL PRIVATE (z_rbfmat,z_diag,z_rbfval,z_rhs1,z_rhs2, ist,   &
!$OMP                   grid_point),                                   &
!$OMP          SHARED  (nproma, rbf_vec_dim_c, nblks_lonlat,           &
!$OMP                   npromz_lonlat, ptr_int_lonlat, ptr_patch,      &
!$OMP                   rbf_vec_kern_c, jg, rbf_vec_scale_c )

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
              cc_e1(:) = ptr_int%cart_edge_coord(ile1,ibe1,:)
              cc_e2(:) = ptr_int%cart_edge_coord(ile2,ibe2,:)
              !
              z_nx1(jc,:) = ptr_patch%edges%primal_cart_normal(ile1,ibe1)%x(:)
              z_nx2(jc,:) = ptr_patch%edges%primal_cart_normal(ile2,ibe2)%x(:)

              ! compute dot product of normal vectors and distance between edge midpoints
              z_nxprod = DOT_PRODUCT(z_nx1(jc,:),z_nx2(jc,:))
              z_dist   = arc_length_v(cc_e1,cc_e2)

              ! set up interpolation matrix
              IF      (rbf_vec_kern_c == 1) THEN
                z_rbfmat(jc,je1,je2) = z_nxprod * gaussi(z_dist,rbf_vec_scale_c(MAX(jg,1)))
              ELSE IF (rbf_vec_kern_c == 3) THEN
                z_rbfmat(jc,je1,je2) = z_nxprod * inv_multiq(z_dist,rbf_vec_scale_c(MAX(jg,1)))
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
          ! convert coordinates of cell center to cartesian vector
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

            cc_e2(:)  = ptr_int%cart_edge_coord(ile2,ibe2,:)

            z_dist = arc_length_v(cc_c(jc,:), cc_e2)

            ! get Cartesian orientation vector
            z_nx3(jc,:) = ptr_patch%edges%primal_cart_normal(ile2,ibe2)%x(:)

            IF (rbf_vec_kern_c == 1) THEN
              z_rbfval(jc,je2) = gaussi(z_dist,rbf_vec_scale_c(MAX(jg,1)))
            ELSE IF (rbf_vec_kern_c == 3) THEN
              z_rbfval(jc,je2) = inv_multiq(z_dist,rbf_vec_scale_c(MAX(jg,1)))
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
    !! @par Revision History
    !! based on
    !!   mo_intp_rbf_coeffs::rbf_compute_coeff_c2grad
    !! developed and tested by Guenther Zaengl, DWD (2009-12-15)
    !! Restructuring F. Prill, DWD (2011-08-18)
    !!
    SUBROUTINE rbf_compute_coeff_c2grad_lonlat (ptr_patch,                   &
      &                                         nblks_lonlat, npromz_lonlat, &
      &                                         ptr_int_lonlat)

      TYPE(t_patch),     INTENT(in)    :: ptr_patch                   ! patch on which computation is performed
      INTEGER,           INTENT(IN)    :: nblks_lonlat, npromz_lonlat ! lon-lat grid blocking info
      TYPE (t_lon_lat_intp), TARGET, INTENT(inout) :: ptr_int_lonlat  ! interpolation state

      ! local variables
      INTEGER  :: je, jcc,                                   &
        &         jc, jb,                                    &  ! integer over lon-lat points
        &         i_startidx, i_endidx,                      &
        &         ile, ibe, ilc, ibc, ilcc, ibcc
      REAL(wp) :: inv_length, coeff(2),                      &
        &         aux_coeff(nproma,rbf_c2grad_dim,2)
      REAL(wp), DIMENSION(:,:,:,:), POINTER :: ptr_coeff

      !--------------------------------------------------------------------

      ptr_coeff => ptr_int_lonlat%rbf_vec_coeff

      ! loop through all patch cells (and blocks)

!$OMP PARALLEL 
!$OMP DO PRIVATE(jb,je,jcc,jc,i_startidx,i_endidx,ile,ibe,ilc,ibc,&
!$OMP            aux_coeff, inv_length, coeff, ilcc, ibcc)
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

              inv_length = ptr_patch%edges%inv_dual_edge_length(ile,ibe)
              coeff(1:2) = ptr_coeff(je,1:2,jc,jb) * inv_length

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
    !> Setup routine for RBF reconstruction at lon-lat grid points for an arbitrary grid.
    ! 
    ! @par Revision History
    !      Initial implementation  by  F.Prill, DWD (2011-08)
    !      Changed for abritrary grids by Rainer Johanni (2011-11)
    !
    ! TODO[FP] Two variables are created globally: "tri_idx" and "owner"
    !          This is necessary for the splitting of the lon-lat grid
    !          points over the PEs and for creating the corresponding
    !          communication patterns. However, after this initial phase,
    !          these data structures should be resized to the local sizes.
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
        &                                 block_size, glb_index, i,             &
        &                                 rl_start, rl_end, i_nchdom, i_startblk
      LOGICAL                          :: l_cutoff_local_domains
      LOGICAL, ALLOCATABLE             :: l_cutoff(:,:)
      INTEGER                          :: array_shape_2d(2), array_shape_3d(3), &
        &                                 array_shape_4d(4)
      TYPE(t_geographical_coordinates) :: cell_center, lonlat_pt
      TYPE(t_cartesian_coordinates)    :: p1, p2
      REAL(wp)                         :: point(2), z_norm, z_nx1(3), z_nx2(3), &
        &                                 max_dist
      INTEGER                          :: ithis_local_pts

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
        &      ptr_int_lonlat%owner(grid%total_dim), STAT=errstat )
      IF (errstat /= SUCCESS) &
        CALL finish (routine, 'allocation for lon-lat point distribution failed')
      ptr_int_lonlat%tri_idx(:,:,:) = 0

      ! Perform query. Note that for distributed patches we receive a
      ! local list of "in_points" actually located on this portion of the
      ! domain.
      IF (dbg_level > 1) &
        CALL message(routine, "proximity query")
      CALL gnat_query_containing_triangles(ptr_patch, gnat_tree, in_points(:,:,:), &
        &                   nproma, nblks_lonlat, npromz_lonlat,    &
        &                   ptr_int_lonlat%tri_idx(:,:,:), min_dist(:,:))
      CALL gnat_merge_distributed_queries(ptr_patch, grid%total_dim, nproma, grid%nblks, &
        &                 min_dist, ptr_int_lonlat%tri_idx(:,:,:), in_points(:,:,:),     &
        &                 ptr_int_lonlat%owner(:), ptr_int_lonlat%nthis_local_pts)
        
      ! set local values for "nblks" and "npromz"
      nblks_lonlat  = (ptr_int_lonlat%nthis_local_pts - 1)/nproma + 1
      npromz_lonlat = ptr_int_lonlat%nthis_local_pts - (nblks_lonlat-1)*nproma
      ! clean up
      CALL gnat_destroy()

      ! After the lon-lat points have been distributed, we know
      ! exactly for each PE the number of points
      CALL allocate_int_state_lonlat_grid(nblks_lonlat, ptr_int_lonlat)

      !-- edge indices of stencil
      ! are available through previous calls to SR rbf_vec_index_cell()
      ! and rbf_c2grad_index()

      !-- copy neighbor indices from standard RBF interpolation
      IF (dbg_level > 1) &
        CALL message(routine, "copy neighbor indices from standard RBF interpolation")

!$OMP PARALLEL SHARED(nblks_lonlat, nproma, npromz_lonlat, &
!$OMP                 ptr_int_lonlat, ptr_int),            &
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

          ptr_int_lonlat%rbf_c2grad_idx(:,jc_lonlat,jb_lonlat) = ptr_int%rbf_c2grad_idx(:,jc,jb)
          ptr_int_lonlat%rbf_c2grad_blk(:,jc_lonlat,jb_lonlat) = ptr_int%rbf_c2grad_blk(:,jc,jb)

        ENDDO
      ENDDO
!$OMP END DO
!$OMP END PARALLEL

      !-- compute interpolation coefficients
      IF (dbg_level > 1) &
        CALL message(routine, "compute lon-lat interpolation coefficients")
      CALL rbf_vec_compute_coeff_lonlat( ptr_patch, ptr_int, ptr_int_lonlat,  &
        &                                in_points, nblks_lonlat, npromz_lonlat )

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
          ptr_int_lonlat%rdist(1, jc, jb) = re*DOT_PRODUCT(p1%x(:), z_nx1)
          ptr_int_lonlat%rdist(2, jc, jb) = re*DOT_PRODUCT(p1%x(:), z_nx2)
        END DO
      END DO

      CALL rbf_compute_coeff_c2grad_lonlat (ptr_patch, nblks_lonlat, &
        &                                   npromz_lonlat, ptr_int_lonlat)

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
        max_dist = 3._gk * REAL(ptr_patch%edges%primal_edge_length(i_startidx,i_startblk)/re, gk)
        
        ! build a logical array, .true. where lon-lat point is outside of
        ! the current patch:
        ALLOCATE(l_cutoff(nproma, nblks_lonlat), stat=errstat)
        IF (errstat /= SUCCESS) CALL finish (routine, 'allocation for working arrays failed')
        l_cutoff(:,:) = (MAX(ABS(ptr_int_lonlat%rdist(1,:,:)), &
          &                  ABS(ptr_int_lonlat%rdist(2,:,:))) >= max_dist)
        
        WHERE (l_cutoff(:,:))
          ptr_int_lonlat%rdist(1,:,:) = 0.
          ptr_int_lonlat%rdist(2,:,:) = 0.
        END WHERE
        
        ! clean up
        DEALLOCATE(l_cutoff, stat=errstat)
        IF (errstat /= SUCCESS)  CALL finish (routine, 'DEALLOCATE of working arrays failed')
      END IF
      
      DEALLOCATE(in_points, rotated_pts, min_dist, stat=errstat)
      IF (errstat /= SUCCESS) THEN
        CALL finish (routine, 'DEALLOCATE of working arrays failed')
      ENDIF

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
    SUBROUTINE rbf_vec_interpol_lonlat( p_vn_in, ptr_patch, ptr_int, &
      &                                 grad_x, grad_y,              &
      &                                 nblks_lonlat, npromz_lonlat, &
      &                                 opt_slev, opt_elev)

      ! INPUT PARAMETERS

      TYPE(t_patch),      TARGET, INTENT(IN) :: ptr_patch                     ! patch on which computation is performed
      ! input normal components of (velocity) vectors at edge midpoints
      REAL(wp),                   INTENT(IN) :: p_vn_in(:,:,:)                ! dim: (nproma,nlev,nblks_e)
      ! Indices of source points and interpolation coefficients
      TYPE (t_lon_lat_intp), TARGET, INTENT(IN) :: ptr_int

      INTEGER,          INTENT(IN)           :: nblks_lonlat, npromz_lonlat   ! blocking info

      INTEGER,          INTENT(in), OPTIONAL :: opt_slev, opt_elev            ! optional vertical start/end level

      ! reconstructed x/y-components of velocity vector
      REAL(wp),         INTENT(INOUT)        :: grad_x(:,:,:), grad_y(:,:,:)  ! dim: (nproma,nlev,nblks_lonlat)

      ! LOCAL VARIABLES
      INTEGER :: slev, elev,                 & ! vertical start and end level
        &        i_startidx, i_endidx,       & ! start/end index
        &        i,                          &
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

!CDIR EXPAND=9
            grad_x(jc,jk,jb) =  &
              SUM( (/ ( ptr_coeff(i,1,jc,jb) *      &
              &         p_vn_in(iidx(i,jc,jb),jk,iblk(i,jc,jb)) , &
              &    i=1,9 ) /) )
!CDIR EXPAND=9
            grad_y(jc,jk,jb) =  &
              SUM( (/ ( ptr_coeff(i,2,jc,jb) *      &
              &         p_vn_in(iidx(i,jc,jb),jk,iblk(i,jc,jb)) , &
              &    i=1,9 ) /) )
            
          ENDDO
        ENDDO
        
      ENDDO
!$OMP END DO
!$OMP END PARALLEL

    END SUBROUTINE rbf_vec_interpol_lonlat


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
    SUBROUTINE rbf_interpol_c2grad_lonlat( p_cell_in, ptr_patch, ptr_int, &
      &                                    nblks_lonlat, npromz_lonlat,   &
      &                                    grad_x, grad_y,                &
      &                                    opt_slev, opt_elev)

      ! !INPUT PARAMETERS
      !
      !  patch on which computation is performed
      !
      TYPE(t_patch), TARGET,     INTENT(in) :: ptr_patch
      ! input cell-based variable for which gradient at cell center is computed
      REAL(wp),                  INTENT(IN) :: p_cell_in(:,:,:) ! dim: (nproma,nlev,nblks_c)
      ! Indices of source points and interpolation coefficients
      TYPE (t_lon_lat_intp), TARGET, INTENT(IN) :: ptr_int
      ! lon-lat grid blocking info
      INTEGER,                   INTENT(IN) :: nblks_lonlat, npromz_lonlat

      ! reconstructed zonal (x) component of gradient vector
      REAL(wp),INTENT(INOUT) :: grad_x(:,:,:) ! dim: (nproma,nlev,nblks_lonlat)
      ! reconstructed zonal (x) component of gradient vector
      REAL(wp),INTENT(INOUT) :: grad_y(:,:,:) ! dim: (nproma,nlev,nblks_lonlat)

      ! optional vertical start/end level
      INTEGER,                   INTENT(IN), OPTIONAL :: opt_slev, opt_elev

      ! Local variables
      INTEGER :: slev, elev,               &  ! vertical start and end level
        &          jc, jb, jk,               &  ! integer over lon-lat points, levels
        &          i_startidx, i_endidx,     &  ! start/end index
        &          i

      INTEGER,  DIMENSION(:,:,:),   POINTER :: iidx, iblk
      REAL(wp), DIMENSION(:,:,:,:), POINTER :: ptr_coeff

      !-----------------------------------------------------------------------

      slev = 1
      elev = UBOUND(p_cell_in,2)
      ! check optional arguments
      IF ( PRESENT(opt_slev) ) slev = opt_slev
      IF ( PRESENT(opt_elev) ) elev = opt_elev

      iidx => ptr_int%rbf_c2grad_idx
      iblk => ptr_int%rbf_c2grad_blk

      ptr_coeff => ptr_int%rbf_c2grad_coeff

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

!CDIR EXPAND=10
            grad_x(jc,jk,jb) =  SUM( &
              & (/ ( ptr_coeff(i, 1, jc, jb) * &
              &      p_cell_in(iidx(i,jc,jb), jk,              &
              &                iblk(i,jc,jb)) , i=1,10 ) /) )
!CDIR EXPAND=10
            grad_y(jc,jk,jb) =  SUM( &
              & (/ ( ptr_coeff(i, 2, jc, jb) * &
              &      p_cell_in(iidx(i,jc,jb), jk,              &
              &                iblk(i,jc,jb)) , i=1,10 ) /) )
            
          ENDDO
        ENDDO
        
      ENDDO
!$OMP END DO
!$OMP END PARALLEL
      
    END SUBROUTINE rbf_interpol_c2grad_lonlat


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
        &        i,                          &
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

!CDIR EXPAND=9
            grad_x(jc,jk,jb) =  &
              SUM( (/ ( ptr_coeff(i,1,jc,jb) *      &
              &         p_vn_in(iidx(i,jc,jb),jk,iblk(i,jc,jb)) , &
              &    i=1,9 ) /) )
!CDIR EXPAND=9
            grad_y(jc,jk,jb) =  &
              SUM( (/ ( ptr_coeff(i,2,jc,jb) *      &
              &         p_vn_in(iidx(i,jc,jb),jk,iblk(i,jc,jb)) , &
              &    i=1,9 ) /) )
            
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
        &          jc, jb, jk,               &  ! integer over lon-lat points, levels
        &          i_startidx, i_endidx,     &  ! start/end index
        &          i

      INTEGER,  DIMENSION(:,:,:),   POINTER :: iidx, iblk
      REAL(wp), DIMENSION(:,:,:,:), POINTER :: ptr_coeff

      !-----------------------------------------------------------------------

      slev = 1
      elev = UBOUND(p_cell_in,2)
      ! check optional arguments
      IF ( PRESENT(opt_slev) ) slev = opt_slev
      IF ( PRESENT(opt_elev) ) elev = opt_elev

      iidx => ptr_int%rbf_c2grad_idx
      iblk => ptr_int%rbf_c2grad_blk

      ptr_coeff => ptr_int%rbf_c2grad_coeff

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

!CDIR EXPAND=10
            grad_x(jc,jk,jb) =  SUM( &
              & (/ ( ptr_coeff(i, 1, jc, jb) * &
              &      p_cell_in(iidx(i,jc,jb), jk,              &
              &                iblk(i,jc,jb)) , i=1,10 ) /) )
!CDIR EXPAND=10
            grad_y(jc,jk,jb) =  SUM( &
              & (/ ( ptr_coeff(i, 2, jc, jb) * &
              &      p_cell_in(iidx(i,jc,jb), jk,              &
              &                iblk(i,jc,jb)) , i=1,10 ) /) )
            
          ENDDO
        ENDDO
        
      ENDDO
!$OMP END DO
!$OMP END PARALLEL

    END SUBROUTINE rbf_interpol_c2grad_lonlat_nl


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
      DO jb=1,nblks_lonlat
        DO jk=slev,elev

          i_startidx = 1
          i_endidx   = nproma
          IF (jb == nblks_lonlat) i_endidx = npromz_lonlat
          FORALL (jc=i_startidx:i_endidx)

            p_lonlat_out(jc,jk,jb) = p_cell_in(ptr_int%tri_idx(1,jc,jb), jk,  &
              &                                ptr_int%tri_idx(2,jc,jb))      &
              &           +  ptr_int%rdist(1,jc,jb) * grad_x(jc,jk,jb)        &
              &           +  ptr_int%rdist(2,jc,jb) * grad_y(jc,jk,jb)

          END FORALL
        END DO
      END DO

    END SUBROUTINE rbf_interpol_lonlat_nl

  END MODULE mo_intp_lonlat
