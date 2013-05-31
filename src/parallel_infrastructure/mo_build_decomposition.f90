MODULE mo_build_decomposition
  USE mo_complete_subdivision
  USE mo_setup_subdivision
  USE mo_ext_decompose_patches
  USE mo_sync,                ONLY: sync_patch_array, sync_idx,disable_sync_checks, &
  &                                 enable_sync_checks, SYNC_E, SYNC_V
  USE mo_grid_config,         ONLY: grid_sphere_radius, n_dom, n_dom_start
  USE mo_mpi
  USE mo_kind
  USE mo_math_utilities
  USE mo_loopindices,         ONLY: get_indices_e
  USE mo_impl_constants
  USE mo_model_domain,        ONLY: p_patch,t_patch_3D,t_patch
  USE mo_model_domimp_patches
  USE mo_parallel_config,     ONLY: p_test_run, l_test_openmp, num_io_procs, division_method
  USE mo_impl_constants,      ONLY: success, MAX_DOM
  USE mo_exception,           ONLY: finish, message, message_text, get_filename_noext

#ifndef __ICON_OCEAN_ONLY__
  USE mo_dump_restore,        ONLY: restore_patches_netcdf
#endif

  IMPLICIT NONE

  PUBLIC :: build_decomposition
  PUBLIC :: complete_patchinfo_oce

  CONTAINS

  SUBROUTINE build_decomposition(nlev,nlevp1,num_lev,num_levp1,nshift,&
      &                          l_is_ocean,l_restore_states, p_patch_3D)
    INTEGER, INTENT(IN) :: nlev,nlevp1,num_lev(MAX_DOM),num_levp1(MAX_DOM),nshift(MAX_DOM)
    LOGICAL, INTENT(IN) :: l_is_ocean,l_restore_states
    TYPE(t_patch_3D), POINTER, OPTIONAL :: p_patch_3D

    TYPE(t_patch), ALLOCATABLE :: p_patch_global(:)
    INTEGER :: error_status,jg
    CHARACTER(LEN=*), PARAMETER :: routine = 'build_decomposition'

    ALLOCATE(p_patch(n_dom_start:n_dom), stat=error_status)
    IF (error_status/=success) THEN
      CALL finish(TRIM(routine), 'allocation of patch failed')
    ENDIF


#ifndef __ICON_OCEAN_ONLY__
    IF(l_restore_states) THEN
      ! Before the restore set p_comm_input_bcast to null
      CALL set_comm_input_bcast(null_comm_type)
      IF( .NOT. my_process_is_mpi_test()) THEN
        !CALL restore_patches_netcdf( p_patch, .TRUE. )
        !CALL set_patch_communicators(p_patch)
        !The 3D-ocean version of previous calls
        CALL restore_patches_netcdf( p_patch, .TRUE. )
        CALL set_patch_communicators(p_patch)

      ELSE
        CALL import_basic_patches(p_patch,num_lev,num_levp1,nshift)
        CALL disable_sync_checks
        CALL complete_patches( p_patch )
        CALL enable_sync_checks
      ENDIF
      ! After the restore is done set p_comm_input_bcast in the
      ! same way as it would be set in parallel_nml_setup when
      ! no restore is wanted:
      CALL set_comm_input_bcast()

    ELSE
#endif
      ! Please note: ldump_dd/lread_dd not (yet?) implemented
      IF(my_process_is_mpi_parallel()) THEN

        IF (division_method(1) > 100) THEN
          ! use ext decomposition library driver
          ALLOCATE(p_patch_global(n_dom_start:n_dom))
          CALL import_basic_patches(p_patch_global,num_lev,num_levp1,nshift)
          CALL ext_decompose_patches(p_patch, p_patch_global)
          DEALLOCATE(p_patch_global)
          IF (l_is_ocean) THEN
            CALL complete_parallel_setup_oce(p_patch)
          ELSE
            CALL complete_parallel_setup
          END IF

        ELSE

          ! use internal decomposition
          !The 3D-ocean version of previous calls
          ALLOCATE(p_patch_global(n_dom_start:n_dom))
          CALL import_basic_patches(p_patch_global,num_lev,num_levp1,nshift)
          IF (l_is_ocean) THEN
            CALL decompose_domain_oce(p_patch,p_patch_global)
          ELSE
            CALL decompose_domain(p_patch_global)
          END IF
          DEALLOCATE(p_patch_global)
          IF (l_is_ocean) THEN
            CALL complete_parallel_setup_oce(p_patch)
          ELSE
            CALL complete_parallel_setup
          END IF
!          CALL finish(routine, "Old decomposition not available for the ocean" )
        ENDIF

      ELSE

        !The 3D-ocean version of previous calls 
        CALL import_basic_patches(p_patch,num_lev,num_levp1,nshift) 

      ENDIF

      ! Complete information which is not yet read or calculated
      CALL complete_patches(p_patch)

#ifndef __ICON_OCEAN_ONLY__
    ENDIF
#endif

    DO jg = n_dom_start, n_dom
      CALL complete_patchinfo_oce(p_patch(jg))
      !The 3D-ocean version of previous calls
      ! Note: this apperas to be problematic for removing the land points
      !CALL complete_patchinfo(p_patch(jg))
    END DO
    !--------------------------------------------        
    ! Setup the information for the physical patches
    CALL setup_phys_patches

    ! In case of a test run: Copy processor splitting to test PE
    !IF(p_test_run) CALL copy_processor_splitting(p_patch)
    !The 3D-ocean version of previous calls 
    IF(p_test_run) CALL copy_processor_splitting(p_patch)
    !--------------------------------------------------------------------------------
    ! 5. Construct interpolation state, compute interpolation coefficients.
    !--------------------------------------------------------------------------------

    !-------------------------------------------------------------------
    ! 7. Finalize domain decomposition
    !-------------------------------------------------------------------   
    IF (my_process_is_mpi_parallel() .AND. .NOT.l_restore_states) THEN

      !CALL finalize_decomposition()
      !The 3D-ocean version of previous calls 
      IF (l_is_ocean) THEN
        CALL finalize_decomposition_oce(p_patch)
      ELSE
        CALL finalize_decomposition
      END IF
    ENDIF

    ! set the horizontal attribute pointer of the 3D patch
    IF (PRESENT(p_patch_3D)) THEN
      ALLOCATE(p_patch_3D, stat=error_status)
      IF (error_status/=SUCCESS) THEN
        CALL finish(TRIM(routine), 'allocation of patch_3D failed')
      ENDIF
      p_patch_3D%p_patch_2D => p_patch
    END IF
  END SUBROUTINE build_decomposition
  !----------------------------------------------------------------------------
  !
  !
  !>
  !! Computes the local orientation of the edge primal normal and dual normal.
  !!
  !! Computes the local orientation of the edge primal normal and dual normal
  !! at the location of the cell centers and vertices.
  !! Moreover, the Cartesian orientation vectors of the edge primal normals
  !! are stored for use in the RBF initialization routines, and inverse
  !! primal and dual edge lengths are computed
  !!
  !! Note: Not clear if all the included calclulations are needed
  !!
  !! @par Revision History
  !!  developed by Guenther Zaengl, 2009-03-31
  !!
  SUBROUTINE complete_patchinfo_oce( ptr_patch)
  !

  !
  !  patch on which computation is performed
  !
  TYPE(t_patch), TARGET, INTENT(inout) :: ptr_patch

  !

  INTEGER :: jb, je!, jc
  INTEGER :: rl_start, rl_end
  INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx, i_nchdom

  INTEGER :: ilc1, ibc1, ilv1, ibv1, ilc2, ibc2, ilv2, ibv2, &
             ilv3, ibv3, ilv4, ibv4!, ile1, ibe1

  REAL(wp) :: z_nu, z_nv, z_lon, z_lat, z_nx1(3), z_nx2(3), z_norm

  TYPE(t_cartesian_coordinates) :: cc_ev3, cc_ev4

  !-----------------------------------------------------------------------

  i_nchdom   = MAX(1,ptr_patch%n_childdom)

  ! !$OMP PARALLEL  PRIVATE(rl_start,rl_end,i_startblk,i_endblk)
  rl_start = 1
  rl_end = min_rlcell

  ! values for the blocking
  i_startblk = ptr_patch%cells%start_blk(rl_start,1)
  i_endblk   = ptr_patch%cells%end_blk(rl_end,i_nchdom)


  rl_start = 1
  rl_end = min_rledge

  ! values for the blocking
  i_startblk = ptr_patch%edges%start_blk(rl_start,1)
  i_endblk   = ptr_patch%edges%end_blk(rl_end,i_nchdom)
  !
  ! First step: compute Cartesian coordinates and Cartesian vectors on full domain
  ! this is needed to vectorize RBF initialization; the existing field carrying
  ! the Cartesian orientation vectors (primal_cart_normal) did not work for that
  ! because it is a derived data type
  ! In addition, the fields for the inverse primal and dual edge lengths are
  ! initialized here.
  !
  ! !$OMP DO PRIVATE(jb,i_startidx,i_endidx,je) ICON_OMP_DEFAULT_SCHEDULE
  DO jb = i_startblk, i_endblk

    CALL get_indices_e(ptr_patch, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, rl_start, rl_end)

    DO je =  i_startidx, i_endidx

      IF(.NOT.ptr_patch%edges%owner_mask(je,jb)) CYCLE

      ! compute Cartesian coordinates (needed for RBF initialization)
      ptr_patch%edges%inv_primal_edge_length(je,jb) = &
        1._wp/ptr_patch%edges%primal_edge_length(je,jb)

    ENDDO

  END DO !block loop
  ! !$OMP END DO

  rl_start = 2
  rl_end = min_rledge

  ! Second step: computed projected orientation vectors and related information
  i_startblk = ptr_patch%edges%start_blk(rl_start,1)
  i_endblk   = ptr_patch%edges%end_blk(rl_end,i_nchdom)

  ! Initialization of lateral boundary points
  IF (ptr_patch%id > 1) THEN
  ! !$OMP WORKSHARE
    ptr_patch%edges%inv_dual_edge_length(:,1:i_startblk)    = 0._wp
    ptr_patch%edges%vertex_idx(:,1:i_startblk,3)            = 0
    ptr_patch%edges%vertex_idx(:,1:i_startblk,4)            = 0
    ptr_patch%edges%vertex_blk(:,1:i_startblk,3)            = 0
    ptr_patch%edges%vertex_blk(:,1:i_startblk,4)            = 0
    ptr_patch%edges%inv_vert_vert_length(:,1:i_startblk)    = 0._wp
    ptr_patch%edges%primal_normal_cell(:,1:i_startblk,:)%v1 = 0._wp
    ptr_patch%edges%dual_normal_cell  (:,1:i_startblk,:)%v1 = 0._wp
    ptr_patch%edges%primal_normal_vert(:,1:i_startblk,:)%v1 = 0._wp
    ptr_patch%edges%dual_normal_vert  (:,1:i_startblk,:)%v1 = 0._wp
    ptr_patch%edges%primal_normal_cell(:,1:i_startblk,:)%v2 = 0._wp
    ptr_patch%edges%dual_normal_cell  (:,1:i_startblk,:)%v2 = 0._wp
    ptr_patch%edges%primal_normal_vert(:,1:i_startblk,:)%v2 = 0._wp
    ptr_patch%edges%dual_normal_vert  (:,1:i_startblk,:)%v2 = 0._wp
  ! !$OMP END WORKSHARE
  ENDIF
  !
  ! loop through all patch edges
  !
  ! !$OMP DO PRIVATE(jb,i_startidx,i_endidx,je,ilc1,ibc1,ilv1,ibv1,ilc2,ibc2,ilv2, &
  ! !$OMP            ibv2,ilv3,ibv3,ilv4,ibv4,z_nu,z_nv,z_lon,z_lat,z_nx1,z_nx2,   &
  ! !$OMP            cc_ev3,cc_ev4,z_norm) ICON_OMP_DEFAULT_SCHEDULE
  DO jb = i_startblk, i_endblk

    CALL get_indices_e(ptr_patch, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, rl_start, rl_end)

    DO je =  i_startidx, i_endidx

      IF(.NOT.ptr_patch%edges%owner_mask(je,jb)) CYCLE

      ! compute inverse dual edge length (undefined for refin_ctrl=1)

      ptr_patch%edges%inv_dual_edge_length(je,jb) = &
        1._wp/ptr_patch%edges%dual_edge_length(je,jb)

      ! compute edge-vertex indices (and blocks) 3 and 4, which
      ! are the outer vertices of cells 1 and 2, respectively,
      ! and the inverse length bewtween vertices 3 and 4

      ilc1 = ptr_patch%edges%cell_idx(je,jb,1)
      ibc1 = ptr_patch%edges%cell_blk(je,jb,1)
      ilc2 = ptr_patch%edges%cell_idx(je,jb,2)
      ibc2 = ptr_patch%edges%cell_blk(je,jb,2)

      ilv1 = ptr_patch%edges%vertex_idx(je,jb,1)
      ibv1 = ptr_patch%edges%vertex_blk(je,jb,1)
      ilv2 = ptr_patch%edges%vertex_idx(je,jb,2)
      ibv2 = ptr_patch%edges%vertex_blk(je,jb,2)

      IF ((ptr_patch%cells%vertex_idx(ilc1,ibc1,1) /= &
           ptr_patch%edges%vertex_idx(je,jb,1) .OR.  &
           ptr_patch%cells%vertex_blk(ilc1,ibc1,1) /= &
           ptr_patch%edges%vertex_blk(je,jb,1)) .AND.  &
          (ptr_patch%cells%vertex_idx(ilc1,ibc1,1) /= &
           ptr_patch%edges%vertex_idx(je,jb,2) .OR.  &
           ptr_patch%cells%vertex_blk(ilc1,ibc1,1) /= &
           ptr_patch%edges%vertex_blk(je,jb,2)) )        THEN

        ptr_patch%edges%vertex_idx(je,jb,3) = ptr_patch%cells%vertex_idx(ilc1,ibc1,1)
        ptr_patch%edges%vertex_blk(je,jb,3) = ptr_patch%cells%vertex_blk(ilc1,ibc1,1)

      ELSE IF ((ptr_patch%cells%vertex_idx(ilc1,ibc1,2) /= &
                ptr_patch%edges%vertex_idx(je,jb,1) .OR.  &
                ptr_patch%cells%vertex_blk(ilc1,ibc1,2) /= &
                ptr_patch%edges%vertex_blk(je,jb,1)) .AND.  &
               (ptr_patch%cells%vertex_idx(ilc1,ibc1,2) /= &
                ptr_patch%edges%vertex_idx(je,jb,2) .OR.  &
                ptr_patch%cells%vertex_blk(ilc1,ibc1,2) /= &
                ptr_patch%edges%vertex_blk(je,jb,2)) )        THEN

        ptr_patch%edges%vertex_idx(je,jb,3) = ptr_patch%cells%vertex_idx(ilc1,ibc1,2)
        ptr_patch%edges%vertex_blk(je,jb,3) = ptr_patch%cells%vertex_blk(ilc1,ibc1,2)

      ELSE IF ((ptr_patch%cells%vertex_idx(ilc1,ibc1,3) /= &
                ptr_patch%edges%vertex_idx(je,jb,1) .OR.  &
                ptr_patch%cells%vertex_blk(ilc1,ibc1,3) /= &
                ptr_patch%edges%vertex_blk(je,jb,1)) .AND.  &
               (ptr_patch%cells%vertex_idx(ilc1,ibc1,3) /= &
                ptr_patch%edges%vertex_idx(je,jb,2) .OR.  &
                ptr_patch%cells%vertex_blk(ilc1,ibc1,3) /= &
                ptr_patch%edges%vertex_blk(je,jb,2)) )        THEN

        ptr_patch%edges%vertex_idx(je,jb,3) = ptr_patch%cells%vertex_idx(ilc1,ibc1,3)
        ptr_patch%edges%vertex_blk(je,jb,3) = ptr_patch%cells%vertex_blk(ilc1,ibc1,3)

      ENDIF

      IF ((ptr_patch%cells%vertex_idx(ilc2,ibc2,1) /= &
           ptr_patch%edges%vertex_idx(je,jb,1) .OR.  &
           ptr_patch%cells%vertex_blk(ilc2,ibc2,1) /= &
           ptr_patch%edges%vertex_blk(je,jb,1)) .AND.  &
          (ptr_patch%cells%vertex_idx(ilc2,ibc2,1) /= &
           ptr_patch%edges%vertex_idx(je,jb,2) .OR.  &
           ptr_patch%cells%vertex_blk(ilc2,ibc2,1) /= &
           ptr_patch%edges%vertex_blk(je,jb,2)) )        THEN

        ptr_patch%edges%vertex_idx(je,jb,4) = ptr_patch%cells%vertex_idx(ilc2,ibc2,1)
        ptr_patch%edges%vertex_blk(je,jb,4) = ptr_patch%cells%vertex_blk(ilc2,ibc2,1)

      ELSE IF ((ptr_patch%cells%vertex_idx(ilc2,ibc2,2) /= &
                ptr_patch%edges%vertex_idx(je,jb,1) .OR.  &
                ptr_patch%cells%vertex_blk(ilc2,ibc2,2) /= &
                ptr_patch%edges%vertex_blk(je,jb,1)) .AND.  &
               (ptr_patch%cells%vertex_idx(ilc2,ibc2,2) /= &
                ptr_patch%edges%vertex_idx(je,jb,2) .OR.  &
                ptr_patch%cells%vertex_blk(ilc2,ibc2,2) /= &
                ptr_patch%edges%vertex_blk(je,jb,2)) )        THEN

        ptr_patch%edges%vertex_idx(je,jb,4) = ptr_patch%cells%vertex_idx(ilc2,ibc2,2)
        ptr_patch%edges%vertex_blk(je,jb,4) = ptr_patch%cells%vertex_blk(ilc2,ibc2,2)

      ELSE IF ((ptr_patch%cells%vertex_idx(ilc2,ibc2,3) /= &
                ptr_patch%edges%vertex_idx(je,jb,1) .OR.  &
                ptr_patch%cells%vertex_blk(ilc2,ibc2,3) /= &
                ptr_patch%edges%vertex_blk(je,jb,1)) .AND.  &
               (ptr_patch%cells%vertex_idx(ilc2,ibc2,3) /= &
                ptr_patch%edges%vertex_idx(je,jb,2) .OR.  &
                ptr_patch%cells%vertex_blk(ilc2,ibc2,3) /= &
                ptr_patch%edges%vertex_blk(je,jb,2)) )        THEN

        ptr_patch%edges%vertex_idx(je,jb,4) = ptr_patch%cells%vertex_idx(ilc2,ibc2,3)
        ptr_patch%edges%vertex_blk(je,jb,4) = ptr_patch%cells%vertex_blk(ilc2,ibc2,3)

      ENDIF

      ilv3 = ptr_patch%edges%vertex_idx(je,jb,3)
      ibv3 = ptr_patch%edges%vertex_blk(je,jb,3)
      ilv4 = ptr_patch%edges%vertex_idx(je,jb,4)
      ibv4 = ptr_patch%edges%vertex_blk(je,jb,4)

  !     IF ( ilv3 > 0 .AND. ilv4 > 0) THEN
  !      write(*,*) "cells:", ilc1, ibc1, ilc2, ibc2
  !      write(*,*) "cell1 vertex idx :", ptr_patch%cells%vertex_idx(ilc1,ibc1,:)
  !      write(*,*) "cell1 vertex blk :", ptr_patch%cells%vertex_blk(ilc1,ibc1,:)
  !      write(*,*) "cell2 vertex idx :", ptr_patch%cells%vertex_idx(ilc2,ibc2,:)
  !      write(*,*) "cell2 vertex blk :", ptr_patch%cells%vertex_blk(ilc2,ibc2,:)
  !      write(*,*) "chosen vertexes:", ilv3,ibv3, ilv4,ibv4
        cc_ev3 = gc2cc(ptr_patch%verts%vertex(ilv3,ibv3))
        cc_ev4 = gc2cc(ptr_patch%verts%vertex(ilv4,ibv4))

        ! inverse length bewtween vertices 3 and 4
        ptr_patch%edges%inv_vert_vert_length(je,jb) = 1._wp / &
          & (grid_sphere_radius * arc_length(cc_ev3,cc_ev4))
  !      ENDIF

      ! next step: compute projected orientation vectors for cells and vertices
      ! bordering to each edge (incl. vertices 3 and 4 intorduced above)

      ! transform orientation vectors at local edge center to Cartesian space
      z_lon = ptr_patch%edges%center(je,jb)%lon
      z_lat = ptr_patch%edges%center(je,jb)%lat

      ! transform primal normal to cartesian vector z_nx1
      z_nx1(:)=ptr_patch%edges%primal_cart_normal(je,jb)%x(:)

      ! transform dual normal to cartesian vector z_nx2
      z_nx2(:)=ptr_patch%edges%dual_cart_normal(je,jb)%x(:)

      ! get location of cell 1

      z_lon = ptr_patch%cells%center(ilc1,ibc1)%lon
      z_lat = ptr_patch%cells%center(ilc1,ibc1)%lat

      ! compute local primal and dual normals at cell 1

      CALL cvec2gvec(z_nx1(1),z_nx1(2),z_nx1(3),z_lon,z_lat,z_nu,z_nv)
      z_norm = SQRT(z_nu*z_nu+z_nv*z_nv)

      ptr_patch%edges%primal_normal_cell(je,jb,1)%v1 = z_nu/z_norm
      ptr_patch%edges%primal_normal_cell(je,jb,1)%v2 = z_nv/z_norm

      CALL cvec2gvec(z_nx2(1),z_nx2(2),z_nx2(3),z_lon,z_lat,z_nu,z_nv)
      z_norm = SQRT(z_nu*z_nu+z_nv*z_nv)

      ptr_patch%edges%dual_normal_cell(je,jb,1)%v1 = z_nu/z_norm
      ptr_patch%edges%dual_normal_cell(je,jb,1)%v2 = z_nv/z_norm

      ! get location of cell 2

      z_lon = ptr_patch%cells%center(ilc2,ibc2)%lon
      z_lat = ptr_patch%cells%center(ilc2,ibc2)%lat

      ! compute local primal and dual normals at cell 2

      CALL cvec2gvec(z_nx1(1),z_nx1(2),z_nx1(3),z_lon,z_lat,z_nu,z_nv)
      z_norm = SQRT(z_nu*z_nu+z_nv*z_nv)

      ptr_patch%edges%primal_normal_cell(je,jb,2)%v1 = z_nu/z_norm
      ptr_patch%edges%primal_normal_cell(je,jb,2)%v2 = z_nv/z_norm

      CALL cvec2gvec(z_nx2(1),z_nx2(2),z_nx2(3),z_lon,z_lat,z_nu,z_nv)
      z_norm = SQRT(z_nu*z_nu+z_nv*z_nv)

      ptr_patch%edges%dual_normal_cell(je,jb,2)%v1 = z_nu/z_norm
      ptr_patch%edges%dual_normal_cell(je,jb,2)%v2 = z_nv/z_norm

      ! get location of vertex 1

      z_lon = ptr_patch%verts%vertex(ilv1,ibv1)%lon
      z_lat = ptr_patch%verts%vertex(ilv1,ibv1)%lat

      ! compute local primal and dual normals at vertex 1

      CALL cvec2gvec(z_nx1(1),z_nx1(2),z_nx1(3),z_lon,z_lat,z_nu,z_nv)
      z_norm = SQRT(z_nu*z_nu+z_nv*z_nv)

      ptr_patch%edges%primal_normal_vert(je,jb,1)%v1 = z_nu/z_norm
      ptr_patch%edges%primal_normal_vert(je,jb,1)%v2 = z_nv/z_norm

      CALL cvec2gvec(z_nx2(1),z_nx2(2),z_nx2(3),z_lon,z_lat,z_nu,z_nv)
      z_norm = SQRT(z_nu*z_nu+z_nv*z_nv)

      ptr_patch%edges%dual_normal_vert(je,jb,1)%v1 = z_nu/z_norm
      ptr_patch%edges%dual_normal_vert(je,jb,1)%v2 = z_nv/z_norm

      ! get location of vertex 2

      z_lon = ptr_patch%verts%vertex(ilv2,ibv2)%lon
      z_lat = ptr_patch%verts%vertex(ilv2,ibv2)%lat

      ! compute local primal and dual normals at vertex 2

      CALL cvec2gvec(z_nx1(1),z_nx1(2),z_nx1(3),z_lon,z_lat,z_nu,z_nv)
      z_norm = SQRT(z_nu*z_nu+z_nv*z_nv)

      ptr_patch%edges%primal_normal_vert(je,jb,2)%v1 = z_nu/z_norm
      ptr_patch%edges%primal_normal_vert(je,jb,2)%v2 = z_nv/z_norm

      CALL cvec2gvec(z_nx2(1),z_nx2(2),z_nx2(3),z_lon,z_lat,z_nu,z_nv)
      z_norm = SQRT(z_nu*z_nu+z_nv*z_nv)

      ptr_patch%edges%dual_normal_vert(je,jb,2)%v1 = z_nu/z_norm
      ptr_patch%edges%dual_normal_vert(je,jb,2)%v2 = z_nv/z_norm

      ! get location of vertex 3

      z_lon = ptr_patch%verts%vertex(ilv3,ibv3)%lon
      z_lat = ptr_patch%verts%vertex(ilv3,ibv3)%lat

      ! compute local primal and dual normals at vertex 3

      CALL cvec2gvec(z_nx1(1),z_nx1(2),z_nx1(3),z_lon,z_lat,z_nu,z_nv)
      z_norm = SQRT(z_nu*z_nu+z_nv*z_nv)

      ptr_patch%edges%primal_normal_vert(je,jb,3)%v1 = z_nu/z_norm
      ptr_patch%edges%primal_normal_vert(je,jb,3)%v2 = z_nv/z_norm

      CALL cvec2gvec(z_nx2(1),z_nx2(2),z_nx2(3),z_lon,z_lat,z_nu,z_nv)
      z_norm = SQRT(z_nu*z_nu+z_nv*z_nv)

      ptr_patch%edges%dual_normal_vert(je,jb,3)%v1 = z_nu/z_norm
      ptr_patch%edges%dual_normal_vert(je,jb,3)%v2 = z_nv/z_norm

      ! get location of vertex 4

      z_lon = ptr_patch%verts%vertex(ilv4,ibv4)%lon
      z_lat = ptr_patch%verts%vertex(ilv4,ibv4)%lat

      ! compute local primal and dual normals at vertex 2

      CALL cvec2gvec(z_nx1(1),z_nx1(2),z_nx1(3),z_lon,z_lat,z_nu,z_nv)
      z_norm = SQRT(z_nu*z_nu+z_nv*z_nv)

      ptr_patch%edges%primal_normal_vert(je,jb,4)%v1 = z_nu/z_norm
      ptr_patch%edges%primal_normal_vert(je,jb,4)%v2 = z_nv/z_norm

      CALL cvec2gvec(z_nx2(1),z_nx2(2),z_nx2(3),z_lon,z_lat,z_nu,z_nv)
      z_norm = SQRT(z_nu*z_nu+z_nv*z_nv)

      ptr_patch%edges%dual_normal_vert(je,jb,4)%v1 = z_nu/z_norm
      ptr_patch%edges%dual_normal_vert(je,jb,4)%v2 = z_nv/z_norm

    ENDDO

  END DO !block loop
  ! !$OMP END DO NOWAIT

  ! !$OMP END PARALLEL

    ! primal_normal_cell must be sync'd before next loop,
    ! so do a sync for all above calculated quantities

    CALL sync_patch_array(SYNC_E,ptr_patch,ptr_patch%edges%inv_primal_edge_length)

    CALL sync_idx(SYNC_E,SYNC_V,ptr_patch,ptr_patch%edges%vertex_idx(:,:,3), &
                                        & ptr_patch%edges%vertex_blk(:,:,3))
    CALL sync_idx(SYNC_E,SYNC_V,ptr_patch,ptr_patch%edges%vertex_idx(:,:,4), &
                                        & ptr_patch%edges%vertex_blk(:,:,4))

    CALL sync_patch_array(SYNC_E,ptr_patch,ptr_patch%edges%inv_dual_edge_length)
    CALL sync_patch_array(SYNC_E,ptr_patch,ptr_patch%edges%inv_vert_vert_length)

    CALL sync_patch_array(SYNC_E,ptr_patch,ptr_patch%edges%primal_normal_cell(:,:,1)%v1)
    CALL sync_patch_array(SYNC_E,ptr_patch,ptr_patch%edges%primal_normal_cell(:,:,2)%v1)
    CALL sync_patch_array(SYNC_E,ptr_patch,ptr_patch%edges%primal_normal_vert(:,:,1)%v1)
    CALL sync_patch_array(SYNC_E,ptr_patch,ptr_patch%edges%primal_normal_vert(:,:,2)%v1)
    CALL sync_patch_array(SYNC_E,ptr_patch,ptr_patch%edges%primal_normal_vert(:,:,3)%v1)
    CALL sync_patch_array(SYNC_E,ptr_patch,ptr_patch%edges%primal_normal_vert(:,:,4)%v1)

    CALL sync_patch_array(SYNC_E,ptr_patch,ptr_patch%edges%dual_normal_cell(:,:,1)%v1)
    CALL sync_patch_array(SYNC_E,ptr_patch,ptr_patch%edges%dual_normal_cell(:,:,2)%v1)
    CALL sync_patch_array(SYNC_E,ptr_patch,ptr_patch%edges%dual_normal_vert(:,:,1)%v1)
    CALL sync_patch_array(SYNC_E,ptr_patch,ptr_patch%edges%dual_normal_vert(:,:,2)%v1)
    CALL sync_patch_array(SYNC_E,ptr_patch,ptr_patch%edges%dual_normal_vert(:,:,3)%v1)
    CALL sync_patch_array(SYNC_E,ptr_patch,ptr_patch%edges%dual_normal_vert(:,:,4)%v1)

    CALL sync_patch_array(SYNC_E,ptr_patch,ptr_patch%edges%primal_normal_cell(:,:,1)%v2)
    CALL sync_patch_array(SYNC_E,ptr_patch,ptr_patch%edges%primal_normal_cell(:,:,2)%v2)
    CALL sync_patch_array(SYNC_E,ptr_patch,ptr_patch%edges%primal_normal_vert(:,:,1)%v2)
    CALL sync_patch_array(SYNC_E,ptr_patch,ptr_patch%edges%primal_normal_vert(:,:,2)%v2)
    CALL sync_patch_array(SYNC_E,ptr_patch,ptr_patch%edges%primal_normal_vert(:,:,3)%v2)
    CALL sync_patch_array(SYNC_E,ptr_patch,ptr_patch%edges%primal_normal_vert(:,:,4)%v2)

    CALL sync_patch_array(SYNC_E,ptr_patch,ptr_patch%edges%dual_normal_cell(:,:,1)%v2)
    CALL sync_patch_array(SYNC_E,ptr_patch,ptr_patch%edges%dual_normal_cell(:,:,2)%v2)
    CALL sync_patch_array(SYNC_E,ptr_patch,ptr_patch%edges%dual_normal_vert(:,:,1)%v2)
    CALL sync_patch_array(SYNC_E,ptr_patch,ptr_patch%edges%dual_normal_vert(:,:,2)%v2)
    CALL sync_patch_array(SYNC_E,ptr_patch,ptr_patch%edges%dual_normal_vert(:,:,3)%v2)
    CALL sync_patch_array(SYNC_E,ptr_patch,ptr_patch%edges%dual_normal_vert(:,:,4)%v2)


  !!$OMP PARALLEL  PRIVATE(rl_start,rl_end,i_startblk,i_endblk)


  !!$OMP END PARALLEL

  END SUBROUTINE complete_patchinfo_oce

END MODULE mo_build_decomposition
