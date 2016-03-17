!>
!! Contains the implementation of interpolation needed for the FEM ice model.
!!
!! @par Revision History
!! Developed  by Einar Olason (2013).
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
!!
!----------------------------
#include "omp_definitions.inc"
!----------------------------
!-----------------------------------------------------------------------
!
!  ! 1) Various initialization routines for the FEM sea-ice model
!  ! 2) Sync routines used in the EVP solver
!
!-----------------------------------------------------------------------
MODULE mo_ice_fem_init
  !-------------------------------------------------------------------------
  !
  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: finish
  USE mo_exception,           ONLY: message
  USE mo_impl_constants,      ONLY: max_char_length
  USE mo_impl_constants,      ONLY: success
!  USE mo_mpi,                 ONLY: get_my_global_mpi_id
!  USE mo_util_dbg_prnt,       ONLY: dbg_print

  USE mo_parallel_config,     ONLY: nproma
  USE mo_model_domain,        ONLY: t_patch, t_patch_3D
  USE mo_grid_subset,         ONLY: t_subset_range, get_index_range
  USE mo_math_constants,      ONLY: rad2deg, deg2rad


  IMPLICIT NONE

  PUBLIC  :: ice_fem_init_vel
  PUBLIC  :: init_fem_wgts
  PUBLIC  :: destruct_fem_wgts
  PUBLIC  :: ice_fem_grid_init
  PUBLIC  :: ice_fem_grid_post
  PUBLIC  :: ice_fem_update_vel_for_restart

  PUBLIC  :: exchange_nod2D
  PUBLIC  :: exchange_nod2Di
  PUBLIC  :: exchange_elem2D

  PRIVATE :: basisfunctions_nod
  PRIVATE :: calc_f_rot

  TYPE(t_patch), POINTER, PRIVATE :: fem_patch

  REAL(wp), ALLOCATABLE, PUBLIC   :: c2v_wgt(:,:,:)
  REAL(wp), ALLOCATABLE, PUBLIC  :: rot_mat(:,:,:,:), rot_mat_3D(:,:)
  ! Longitude and latitude of the north pole of the rotated grid
  ! These values put the north pole on the Indonesian/Malasian island Kalimantan and the south pole
  ! in north-west Brazil, near the border to Venezuela and Colombia.
  REAL(wp), PARAMETER, PUBLIC :: pollon = 114._wp*deg2rad, pollat = 0._wp
  !  the North pole is centered on Greenland (40◦W, 75◦N) and the South pole on Antarctica
  !  was used to test consistency of the results with respect to the change in the rotated grid
  !  REAL(wp), PARAMETER, PRIVATE :: pollon = -40._wp*deg2rad, pollat = 75._wp*deg2rad

CONTAINS

  !-------------------------------------------------------------------------
  !
  !> Constructor of weights for FEM sea-ice model
  !! We calculate only c2v_wgt (PRIVATE), which is used by cells2verts_scalar
  !!
  !! @par Revision History
  !! Developed by Einar Olason, MPI-M (2013-06-05)
  !! Modified by Vladimir Lapin, MPI-M (2015-11-04)
  !<Optimize:inUse>
  SUBROUTINE init_fem_wgts(p_patch_3D)
    TYPE(t_patch_3D), TARGET, INTENT(IN)    :: p_patch_3D

    !Local variables
    CHARACTER(LEN=max_char_length), PARAMETER :: routine = 'mo_ice_fem_utils:init_fem_wgts'

    ! Patch
    TYPE(t_patch), POINTER :: p_patch
    INTEGER,  DIMENSION(:,:,:),   POINTER :: iidx, iblk

    ! Indexing
    INTEGER :: ist, nblks_v
    INTEGER :: jb, jv, ji
    INTEGER :: i_startidx_v, i_endidx_v
    INTEGER :: cell_index, cell_block

    ! The denominator
    REAL(wp) :: rdeno, rsum

    !-------------------------------------------------------------------------
    CALL message(TRIM(routine), 'start' )

    p_patch => p_patch_3D%p_patch_2D(1)
    nblks_v = p_patch%nblks_v

    ALLOCATE(c2v_wgt(nproma,6,nblks_v),STAT=ist)
    IF (ist /= SUCCESS) THEN
      CALL finish (routine,'allocating c2v_wgt failed')
    ENDIF
    c2v_wgt=0.0_wp

    ! Weights for the cells2verts_scalar call
!     i_startblk = p_patch%verts%all%start_block
!     i_endblk   = p_patch%verts%all%end_block
    ! Indexing of neighbours
    iidx => p_patch%verts%cell_idx
    iblk => p_patch%verts%cell_blk
    !
    ! loop through all patch edges
    !
!ICON_OMP_PARALLEL_DO PRIVATE(jb,i_startidx_v,i_endidx_v,jv,ji,rsum,     &
!ICON_OMP        rdeno,cell_index,cell_block) ICON_OMP_DEFAULT_SCHEDULE
  DO jb =  p_patch%verts%in_domain%start_block, p_patch%verts%in_domain%end_block
    CALL get_index_range(p_patch%verts%in_domain, jb, i_startidx_v, i_endidx_v)

    DO jv =  i_startidx_v, i_endidx_v

        rsum = 0.0_wp

        DO ji = 1, 6
    !    in case of pentagons and boundary nodes:
    !    cell_index = 0 when use_duplicated_connectivity = .false. (always the case for the ocean)
    !    cell_index = dummy values (duplicates) when use_duplicated_connectivity = .true.
          cell_index = iidx(jv,jb,ji)
          cell_block = iblk(jv,jb,ji)
          IF (cell_index > 0)                                      &
            & rsum = rsum +                                        &
            &      p_patch_3D%wet_c  (cell_index,1,cell_block) *   &
            &      p_patch%cells%area(cell_index,  cell_block)

        ENDDO

        rdeno =  1._wp/MAX( TINY(rdeno), rsum )

        DO ji = 1, 6
          cell_index = iidx(jv,jb,ji)
          cell_block = iblk(jv,jb,ji)
          IF (cell_index > 0)                                 &
            & c2v_wgt(jv,ji,jb) =                             &
            &   p_patch_3D%wet_c  (cell_index,1,cell_block) * &
            &   p_patch%cells%area(cell_index,  cell_block) * rdeno

        ENDDO

      ENDDO
 
    ENDDO
!ICON_OMP_END_PARALLEL_DO

    CALL message (TRIM(routine), 'end')        

  END SUBROUTINE init_fem_wgts
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !
  !> Destructor of patch for FEM sea-ice model, deallocates c2v_wgt
  !!
  !! @par Revision History
  !! Developed by Einar Olason, MPI-M (2013-06-05)
  !
  SUBROUTINE destruct_fem_wgts

    !Local variables
    INTEGER :: ist
    CHARACTER(LEN=max_char_length), PARAMETER :: routine = 'mo_ice_fem_utils:destruct_fem_wgts'
    !-------------------------------------------------------------------------
    CALL message(TRIM(routine), 'start' )

    DEALLOCATE(c2v_wgt,STAT=ist)
    IF (ist /= SUCCESS) THEN
      CALL finish (routine,'deallocating c2v_wgt failed')
    ENDIF

    DEALLOCATE(rot_mat,STAT=ist)
    IF (ist /= SUCCESS) THEN
      CALL finish (routine,'deallocating rot_mat failed')
    ENDIF

    DEALLOCATE(rot_mat_3D,STAT=ist)
    IF (ist /= SUCCESS) THEN
      CALL finish (routine,'deallocating rot_mat_3D failed')
    ENDIF

    CALL message (TRIM(routine), 'end')        

  END SUBROUTINE destruct_fem_wgts
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !
  !> Initialise the FEM grid
  !! This replaces the routine Sergey used to read the grid from file. We give values to
  !! coord_nod2D, nod2D, index_nod2D, elem2D_nodes and elem2D from the ICON grid.
  !!
  !! @par Revision History
  !! Developed by Einar Olason, MPI-M (2013-08-05)
  !
  SUBROUTINE ice_fem_grid_init(p_patch_3D)

  USE mo_math_utilities,      ONLY: rotate_latlon, rotate_latlon_vec!, disp_new_vect

    ! For the AWI FEM
    USE mo_ice_mesh,           ONLY: coord_nod2D, nod2D, index_nod2D
    USE mo_ice_elements,       ONLY: elem2D_nodes, elem2D, nod2D_elems
    USE mo_physical_constants, ONLY: earth_angular_velocity
    USE mo_ice_mesh,           ONLY: coriolis_nod2D
   ! USE mo_mpi
!    USE mo_sync,                ONLY: SYNC_C, sync_patch_array, sync_patch_array_mult

    TYPE(t_patch_3D), TARGET, INTENT(in) :: p_patch_3D

    ! local variables
    CHARACTER(LEN=max_char_length), PARAMETER :: routine = 'mo_ice_fem_utils:ice_fem_grid_init'

    ! Patch
    TYPE(t_patch), POINTER :: p_patch

    ! Indexing
    INTEGER :: k, jb, jc, jv
    INTEGER :: i_startidx_c, i_endidx_c
    INTEGER :: i_startidx_v, i_endidx_v
    INTEGER :: cell_index, cell_block

    ! Temporary variables/buffers and flags
    INTEGER :: verts(nproma,p_patch_3D%p_patch_2D(1)%nblks_v)
    INTEGER :: elems(nproma,p_patch_3D%p_patch_2D(1)%nblks_e)
    INTEGER :: buffy(nproma*p_patch_3D%p_patch_2D(1)%nblks_v)
    INTEGER :: ist
    REAL(wp) :: lat, lon, cos_d, sin_d
    REAL(wp) :: cos_lonp, sin_lonp, cos_latp, sin_latp

    ! Masking for the halo-points
    INTEGER :: halo_mask(nproma, p_patch_3D%p_patch_2D(1)%nblks_v)

   ! REAL(wp) :: test(nproma, p_patch_3D%p_patch_2D(1)%alloc_cell_blocks, 3)
   ! INTEGER :: my_mpi_id

!--------------------------------------------------------------------------------------------------

    p_patch => p_patch_3D%p_patch_2D(1)
    fem_patch => p_patch
    ! my_mpi_id = get_my_global_mpi_id()

    ! nod2D is the number of nodes = number of vertices
    nod2D = p_patch%n_patch_verts
    ALLOCATE(coord_nod2D(2,nod2D), index_nod2D(nod2D),STAT=ist)
    IF (ist /= SUCCESS) THEN
      CALL finish (routine,'allocating coord_nod2D failed')
    ENDIF
    ALLOCATE(coriolis_nod2D(nod2D),STAT=ist)
    IF (ist /= SUCCESS) THEN
      CALL finish (routine,'allocating coriolis_nod2D failed')
    ENDIF
    ALLOCATE(rot_mat(nproma,p_patch%nblks_v,2,2),STAT=ist)
    IF (ist /= SUCCESS) THEN
      CALL finish (routine,'allocating rot_mat failed')
    ENDIF
    ALLOCATE(rot_mat_3D(3,3),STAT=ist)
    IF (ist /= SUCCESS) THEN
      CALL finish (routine,'allocating rot_mat_3D failed')
    ENDIF


    ! Go through all the vertices and assign coordinates and land mask to coord_nod2D and
    ! index_nod2D
    k=0
    DO jb = 1,p_patch%nblks_v
      CALL get_index_range(p_patch%verts%all, jb, i_startidx_v, i_endidx_v) 
      DO jv = i_startidx_v,i_endidx_v
        k=k+1
        verts(jv,jb) = k
        lat = p_patch%verts%vertex(jv,jb)%lat
        lon = p_patch%verts%vertex(jv,jb)%lon

        ! Calculate the rotation matrix once
        ! see p.27 at http://www.cosmo-model.org/content/model/documentation/core/cosmoDyncsNumcs.pdf
        CALL rotate_latlon_vec( lon, lat, pollon, pollat, sin_d, cos_d )
        rot_mat(jv,jb,1,:) = (/ cos_d, sin_d /)
        rot_mat(jv,jb,2,:) = (/-sin_d, cos_d /)

        ! Use a rotated grid with the north pole at (pollon,pollat)
        CALL rotate_latlon(lat, lon, pollat, pollon)

!==========================================================================
        coriolis_nod2D(k) = calc_f_rot(lat,lon,pollat,earth_angular_velocity)

        ! x-coords in degrees
        coord_nod2D(1,k) = lon*rad2deg
        ! y-coords in degrees
        coord_nod2D(2,k) = lat*rad2deg
        ! border (and land) points (0 or 1)
        index_nod2D(k)   =  MIN(1, MAX(0, 1+p_patch_3D%surface_vertex_sea_land_mask(jv,jb)))
      ENDDO
    ENDDO

    ! the rotation matrix is a modified version of r2*r1 from disp_new_vect
    cos_lonp=COS(pollon)
    sin_lonp=SIN(pollon)
    cos_latp=COS(pollat)
    sin_latp=SIN(pollat)

    rot_mat_3D(1,:) =  (/  cos_lonp*sin_latp,  sin_lonp*sin_latp, -cos_latp /)
    rot_mat_3D(2,:) =  (/ -sin_lonp ,      cos_lonp,       0.0_wp /)
    rot_mat_3D(3,:) =  (/  cos_lonp*cos_latp,  sin_lonp*cos_latp,  sin_latp /)

    ! elem2D is the number of elements = number of cells
   ! test(:,:,:) = 0.0_wp
    elem2D = p_patch%n_patch_cells
    ALLOCATE(elem2D_nodes(3,elem2D),STAT=ist)
    IF (ist /= SUCCESS) THEN
      CALL finish (routine,'allocating coord_nod2D failed')
    ENDIF

    ! Establish the connectivity between elements and nodes. The array elem2D_nodes lists the
    ! nodes each element is made up of.
    k=0
    DO jb = p_patch%cells%all%start_block, p_patch%cells%all%end_block
      CALL get_index_range(p_patch%cells%all, jb, i_startidx_c, i_endidx_c) 
      DO jc = i_startidx_c,i_endidx_c
        k=k+1
        elem2D_nodes(1,k) =   &
          &   verts(p_patch%cells%vertex_idx(jc,jb,1),p_patch%cells%vertex_blk(jc,jb,1))
        elem2D_nodes(2,k) =   &
          &   verts(p_patch%cells%vertex_idx(jc,jb,2),p_patch%cells%vertex_blk(jc,jb,2))
        elem2D_nodes(3,k) =   &
          &   verts(p_patch%cells%vertex_idx(jc,jb,3),p_patch%cells%vertex_blk(jc,jb,3))

!         test(jc, jb, 1) = &
!           & REAL(p_patch%verts%glb_index(verts(p_patch%cells%vertex_idx(jc,jb,1),p_patch%cells%vertex_blk(jc,jb,1))), wp)
!         test(jc, jb, 2) = &
!           & REAL(p_patch%verts%glb_index(verts(p_patch%cells%vertex_idx(jc,jb,2),p_patch%cells%vertex_blk(jc,jb,2))), wp)
!         test(jc, jb, 3) = &
!           & REAL(p_patch%verts%glb_index(verts(p_patch%cells%vertex_idx(jc,jb,3),p_patch%cells%vertex_blk(jc,jb,3))), wp)
!
!         write(0,*) my_mpi_id, "::", p_patch%cells%glb_index(k), ":",  test(jc, jb, 1), test(jc, jb, 2), test(jc, jb, 3)

      ENDDO
    ENDDO

    ! allocate verts->cells connectivity index for the restructured solver loop
    ALLOCATE(nod2D_elems(6,nod2D),STAT=ist) ! 6 neighbouring cells by default, can be less for pentagons and boundary nodes
    IF (ist /= SUCCESS) THEN
      CALL finish (routine,'allocating nod2D_elems failed')
    ENDIF
    ! initialize with zeros
    nod2D_elems = 0

    k=0
    DO jb = p_patch%cells%all%start_block, p_patch%cells%all%end_block
      CALL get_index_range(p_patch%cells%all, jb, i_startidx_c, i_endidx_c)
      DO jc = i_startidx_c,i_endidx_c
        k=k+1
        elems(jc,jb) = k
      ENDDO
    ENDDO
    ! Establish the connectivity between elements and nodes.
    ! Array nod2D_elems lists the elements that share a given node.
    k=0
    DO jb = p_patch%verts%all%start_block, p_patch%verts%all%end_block
      CALL get_index_range(p_patch%verts%all, jb, i_startidx_v, i_endidx_v)
      DO jv = i_startidx_v, i_endidx_v
        k=k+1
        DO jc = 1, 6
          cell_index = p_patch%verts%cell_idx(jv,jb,jc)
          cell_block = p_patch%verts%cell_blk(jv,jb,jc)

          IF (cell_index > 0)                                      &
             nod2D_elems(jc,k) = elems(cell_index,cell_block)

        END DO
      END DO
    END DO

!    write(0,*) "sync test(jc, jb, 1)..."
!    CALL sync_patch_array(SYNC_C, p_patch, test(:, :, 1))
!    write(0,*) "sync test(jc, jb, 2)..."
!    CALL sync_patch_array(SYNC_C, p_patch, test(:, :, 2))
!    write(0,*) "sync test(jc, jb, 3)..."
!    CALL sync_patch_array(SYNC_C, p_patch, test(:, :, 3))
!    write(0,*) "sync test( done."


    ! Mask out halopoints
    ! This is neccesary so that Sergey's EVP routine doesn't try to calculate velocities on the
    ! halo-points. In his code ponts with index_nod2D=0 are used for calculation and the velocity
    ! at points with index_nod2D=1 (land borders) is set to zero (no-slip). The velocity at points
    ! with index_nod2D>1 will therefore not be calculated and not set to zero.
    halo_mask = 0
    DO jb = p_patch%verts%not_owned%start_block, p_patch%verts%not_owned%end_block
      CALL get_index_range(p_patch%verts%not_owned, jb, i_startidx_v, i_endidx_v)
        DO jv = i_startidx_v, i_endidx_v
          halo_mask(jv,jb) = 2
        ENDDO
    ENDDO

    buffy = RESHAPE(halo_mask(:,:), SHAPE(buffy))
    index_nod2D = index_nod2D+buffy(1:SIZE(index_nod2D))

  END SUBROUTINE ice_fem_grid_init

  !-------------------------------------------------------------------------
  !
  !> Copy basis functions arrays (bafux, bafuy) to (bafux_nod, bafuy_nod)
  !! which are indexed for each node, rather than for each element
  !! Those arrays are necessary for the OMP-parallelized version of stress2rhs_omp
  !! Where we loop over the nodes to update the rhs value in the momentum equatoins.
  !!
  !! @par Revision History
  !! Developed by Vladimir Lapin, MPI-M (2015-11-30)
  !
  SUBROUTINE basisfunctions_nod

    USE mo_ice_mesh,           ONLY: nod2D
    USE mo_ice_elements,       ONLY: bafux, bafuy, bafux_nod, bafuy_nod,    &
                                   & elem2D_nodes, nod2D_elems

    ! Indexing
    INTEGER :: je, jn, k
    INTEGER :: ist
    ! Temporary var
    INTEGER ::  elem, elnodes(3)
    CHARACTER(LEN=max_char_length), PARAMETER :: routine = 'mo_ice_fem_utils:basisfunctions_nod'
    LOGICAL :: found

    ! allocate bafux_nod, bafuy_nod
    ALLOCATE(bafux_nod(6,nod2D),bafuy_nod(6,nod2D),STAT=ist) ! 6 neighbouring cells by default, can be less for pentagons and boundary nodes
    IF (ist /= SUCCESS) THEN
      CALL finish (routine,'allocating bafux_nod/bafuy_nod failed')
    ENDIF
    ! initialize with zeros
    bafux_nod = 0.0_wp
    bafuy_nod = 0.0_wp

!ICON_OMP_PARALLEL_DO PRIVATE(jn,je,elem,elnodes,k) ICON_OMP_DEFAULT_SCHEDULE
    DO jn=1, nod2D
         DO je=1,6
            elem=nod2D_elems(je,jn)

            IF (elem > 0) THEN
              elnodes=elem2D_nodes(:,elem)

              ! find the node for the giveN element elem, that matches jn
              ! just brute force search
              found = .false.
              DO k=1,3
                IF (jn==elnodes(k)) THEN
                    found = .true.
                    exit
                ENDIF
              END DO

              IF (found) THEN
                  bafux_nod(je,jn) = bafux(k,elem)
                  bafuy_nod(je,jn) = bafuy(k,elem)
              ELSE
                CALL finish(routine, "FEM element-vertex connectivity inconsistency")
              ENDIF

            END IF
         END DO
    END DO
!ICON_OMP_END_PARALLEL_DO

  END SUBROUTINE basisfunctions_nod

  !-------------------------------------------------------------------------
  !> Calculates the Coriolis parameters at a point on the rotated grid with
  !> the new pole at (pollon, pollat).
  !-------------------------------------------------------------------------
  !!
  !! @par Revision History
  !! Developed by Vladimir Lapin, MPI-M (2015-08-11)
  !
  function calc_f_rot(lat, lon, pollat, omega)

  implicit none

  real(wp), intent(in)  :: lat, lon, pollat, omega
  real(wp)              :: calc_f_rot

  calc_f_rot = 2._wp*omega*(-cos(lon)*cos(lat)*cos(pollat)+sin(lat)*sin(pollat))

  end function calc_f_rot

  !-------------------------------------------------------------------------
  !
  !> Synchronisation routine for the FEM
  !! This replaces the routine Sergey used to synchronize arrays across CPUs. It is just a wrapper
  !! around sync_patch_array using the PRIVATE variable fem_patch as patch and reshaping the
  !! variable before and after syncing.
  !!
  !! @par Revision History
  !! Developed by Einar Olason, MPI-M (2013-08-05)
  !
  SUBROUTINE exchange_nod2D(u_ice)

    USE mo_sync,                ONLY: SYNC_V, sync_patch_array!, sync_patch_array_mult

    REAL(wp), INTENT(INOUT) :: u_ice(fem_patch%n_patch_verts)

    ! local variables
    ! Temporary variables/buffers and pad
    REAL(wp) :: u(nproma, fem_patch%nblks_v)
    REAL(wp) :: buffy(nproma*fem_patch%nblks_v)
    REAL(wp), ALLOCATABLE :: pad(:)

  !-------------------------------------------------------------------------

    ! Reshape and copy ice velocities to ICON variables
    ALLOCATE(pad(nproma*fem_patch%nblks_v - fem_patch%n_patch_verts))
    pad = -9999._wp
    u=RESHAPE(u_ice, SHAPE(u), pad)
    DEALLOCATE(pad)

    CALL sync_patch_array(SYNC_V, fem_patch, u(:,:))

    ! Reshape and copy ice velocity to FEM model variables
    buffy = RESHAPE(u, SHAPE(buffy))
    u_ice = buffy(1:SIZE(u_ice))

  END SUBROUTINE exchange_nod2D

  !-------------------------------------------------------------------------
  !
  !> Synchronisation routine for the FEM having an integer array as input
  !! For testing purposes only
  !!
  !! @par Revision History
  !! Developed by Einar Olason, MPI-M (2013-08-05)
  !
  SUBROUTINE exchange_nod2Di(index)
    USE mo_sync,                ONLY: SYNC_V, sync_patch_array!, sync_patch_array_mult

    INTEGER, INTENT(INOUT) :: index(fem_patch%n_patch_verts)

    ! local variables
    ! Temporary variables/buffers and pad
    INTEGER :: i(nproma, fem_patch%nblks_v)
    INTEGER :: buffy(nproma*fem_patch%nblks_v)
    INTEGER, ALLOCATABLE :: pad(:)

  !-------------------------------------------------------------------------

    ! Reshape and copy ice velocities to ICON variables
    ALLOCATE(pad(nproma*fem_patch%nblks_v - fem_patch%n_patch_verts))
    pad = -9999
    i=RESHAPE(index, SHAPE(i), pad)
    DEALLOCATE(pad)

    CALL sync_patch_array(SYNC_V, fem_patch, i(:,:))

    ! Reshape and copy ice velocity to FEM model variables
    buffy = RESHAPE(i, SHAPE(buffy))
    index = buffy(1:SIZE(index))

  END SUBROUTINE exchange_nod2Di

  !-------------------------------------------------------------------------
  !
  !> Synchronisation routine for the FEM elements
  !! For testing purposes only
  !!
  !! @par Revision History
  !! Developed by Einar Olason, MPI-M (2013-08-05)
  !
  SUBROUTINE exchange_elem2D(sigma)
    USE mo_sync,                ONLY: SYNC_C, sync_patch_array!, sync_patch_array_mult

    REAL(wp), INTENT(INOUT) :: sigma(fem_patch%n_patch_cells)

    ! local variables
    ! Temporary variables/buffers and pad
    REAL(wp) :: s(nproma, fem_patch%alloc_cell_blocks)
    REAL(wp) :: buffy(nproma*fem_patch%alloc_cell_blocks)
    REAL(wp), ALLOCATABLE :: pad(:)

  !-------------------------------------------------------------------------

    ! Reshape and copy fem variable to ICON variable
    ALLOCATE(pad(nproma*fem_patch%alloc_cell_blocks - fem_patch%n_patch_cells))
    pad = -9999._wp
    s=RESHAPE(sigma, SHAPE(s), pad)
    DEALLOCATE(pad)

    CALL sync_patch_array(SYNC_C, fem_patch, s(:,:))

    ! Reshape and copy ice velocity to FEM model variables
    buffy = RESHAPE(s, SHAPE(buffy))
    sigma = buffy(1:SIZE(sigma))

  END SUBROUTINE exchange_elem2D

  !-------------------------------------------------------------------------
  !
  !> Post-initialisation work for the FEM grid
  !! In order for the FEM model to pass the p-test we need to sort the lists myList_elem2D and
  !! myList_nod2D according to the global index. This is neccesary so that summation in stress2rhs
  !! (in ice_dyn_fem/ice_evp.f90) is always done in the same order on all CPUs.
  !! In addition to this we also need to set lmass_matrix to the dual area size, since lmass_matrix
  !! does not synchronise across CPUs the way it is calculated here. We also set voltriangle to the
  !! cell area, since the cell area is more accurately calculated by ICON than the FEM code.
  !!
  !! @par Revision History
  !! Developed by Einar Olason, MPI-M (2013-08-05)
  !
  SUBROUTINE ice_fem_grid_post(p_patch)
    USE mo_ice_parsup,          ONLY: myList_nod2D, myList_elem2D
    USE mo_ice,                 ONLY: lmass_matrix
    USE mo_ice_elements,        ONLY: voltriangle
    USE mo_ice_mesh,            ONLY: cos_elem2D, sin_elem2D, coriolis_nod2D!, nod2D, elem2D

    USE mo_util_sort,           ONLY: quicksort

    TYPE(t_patch), TARGET, INTENT(IN)    :: p_patch

    ! Local variables
    ! Patch and ranges
    TYPE(t_subset_range), POINTER :: all_verts
    ! Indexing
    INTEGER :: i_startidx_v, i_endidx_v
    INTEGER :: i_startidx_c, i_endidx_c
    INTEGER :: k, jb, jv, jc

    ! Global lists
    INTEGER :: globList_nod2D(p_patch%n_patch_verts)
    INTEGER :: globList_elem2D(p_patch%n_patch_cells)

    ! Temporary variables/buffers
    REAL(wp) :: buffy_v(nproma*p_patch%nblks_v)
    REAL(wp) :: buffy_c(nproma*p_patch%alloc_cell_blocks), buffy_c2(nproma*p_patch%alloc_cell_blocks)

  !-------------------------------------------------------------------------

    ! Copy basis functions arrays (bafux, bafuy) to (bafux_nod, bafuy_nod)
    ! Used by the OMP-parallelized loop in stress2rhs_omp
    CALL basisfunctions_nod

!!!! check that input myList_nod2D==(1:nod2D) and myList_elem2D==(1:elem2D)
!
!  write(0,*) "max(myList_nod2D - 1:nod2D)", maxval(myList_nod2D - (/(k, k=1, nod2D)/) ), minval(myList_nod2D - (/(k, k=1, nod2D)/) )
!  write(0,*) "max(myList_elem2D - 1:elem2D)", maxval(myList_elem2D - (/(k, k=1, elem2D)/) ), minval(myList_elem2D - (/(k, k=1, elem2D)/) )

    all_verts => p_patch%verts%all
    
    ! Sort the list myList_nod2D using the global index
    k=0
    DO jb = 1,p_patch%nblks_v
      CALL get_index_range(all_verts, jb, i_startidx_v, i_endidx_v)
      DO jv = i_startidx_v,i_endidx_v
        k=k+1
        globList_nod2D(k) = p_patch%verts%decomp_info%glb_index(k)
      ENDDO
    ENDDO

    CALL quicksort(globList_nod2D,myList_nod2D)

    ! Sort the list myList_elem2D using the global index
    k=0
    DO jb = p_patch%cells%all%start_block, p_patch%cells%all%end_block
      CALL get_index_range(p_patch%cells%all, jb, i_startidx_c, i_endidx_c) 
      DO jc = i_startidx_c,i_endidx_c
        k=k+1
        globList_elem2D(k) = p_patch%cells%decomp_info%glb_index(k)
      ENDDO
    ENDDO

    CALL quicksort(globList_elem2D,myList_elem2D)

!    k=0
!    DO jb = p_patch%cells%all%start_block, p_patch%cells%all%end_block
!      CALL get_index_range(p_patch%cells%all, jb, i_startidx_c, i_endidx_c)
!      DO jc = i_startidx_c,i_endidx_c
!        k=k+1
!        write(0,*) get_my_global_mpi_id(), ":: myList_elem2D(k)=", myList_elem2D(k), " globList_elem2D(k)=", globList_elem2D(k),  " vs k=",  k
!      ENDDO
!    ENDDO
!   CALL finish ('check finished')
!

    ! Re-create lmass_matrix as one third of the sum of areas of triangles that share vertex n
    ! Or in ICON-terms as the dual area

    buffy_v = RESHAPE(p_patch%verts%dual_area(:,:), SHAPE(buffy_v))
    lmass_matrix = buffy_v(1:SIZE(lmass_matrix))

    ! voltriangle should also be the cell area
    buffy_c = RESHAPE(p_patch%cells%area(:,:), SHAPE(buffy_c))
    voltriangle = buffy_c(1:SIZE(voltriangle))

    ! coriolis_nod2D, cos_elem2D, sin_elem2D need to be reset because of grid rotation

    ! IS NOT NECESSARY, already handled by the previously called intitialization routines
    ! But for later restructuring -- keep this way, instead of recalculating
    buffy_v = RESHAPE(p_patch%verts%f_v(:,:), SHAPE(buffy_v))
    coriolis_nod2D = buffy_v(1:SIZE(coriolis_nod2D))

    ! Calc cos&sin for metrics terms were previously incorrectly calculated on the original grid
    !    buffy_c = RESHAPE(p_patch%cells%center(:,:)%lat, SHAPE(buffy_c))
    !    cos_elem2D = COS(buffy_c(1:SIZE(cos_elem2D)))
    !    buffy_c = RESHAPE(p_patch%cells%center(:,:)%lat, SHAPE(buffy_c))
    !    sin_elem2D = SIN(buffy_c(1:SIZE(sin_elem2D)))

    ! Insted, one can either use the formulas below
    ! buffy_c = RESHAPE(p_patch%cells%center(:,:)%lat, SHAPE(buffy_c))
    ! buffy_c2 = RESHAPE(p_patch%cells%center(:,:)%lon, SHAPE(buffy_c))

    ! sin_elem2D = SIN(buffy_c(1:SIZE(sin_elem2D)))*SIN(pollat) + COS(buffy_c(1:SIZE(sin_elem2D)))*COS(pollat)*COS(buffy_c2(1:SIZE(sin_elem2D))-pollon)
    ! cos_elem2D = COS(ASIN( SIN(buffy_c(1:SIZE(sin_elem2D)))*SIN(pollat) + COS(buffy_c(1:SIZE(sin_elem2D)))*COS(pollat)*COS(buffy_c2(1:SIZE(sin_elem2D))-pollon) ))

    ! or simply let the fem init routines to handle this (how it is currently done)

  END SUBROUTINE ice_fem_grid_post

  !-------------------------------------------------------------------------
  !
  !> Initialize u_ice, v_ice in case of a restart
  !!
  !! @par Revision History
  !! Developed by Einar Olason, MPI-M (2013-08-05)
  !
  SUBROUTINE ice_fem_init_vel(p_patch, p_ice)
    USE mo_ice,                 ONLY: u_ice, v_ice
    USE mo_sea_ice_types,       ONLY: t_sea_ice

    TYPE(t_patch), TARGET, INTENT(IN)    :: p_patch
    TYPE(t_sea_ice),       INTENT(IN)    :: p_ice

    ! Local variables
    ! Patch and ranges
    TYPE(t_subset_range), POINTER :: all_verts
    ! Indexing
    INTEGER :: i_startidx_v, i_endidx_v
    INTEGER :: jk, jb, jv

  !-------------------------------------------------------------------------

    all_verts => p_patch%verts%all

    jk=0
    DO jb = all_verts%start_block, all_verts%end_block
      CALL get_index_range(all_verts, jb, i_startidx_v, i_endidx_v)
      DO jv = i_startidx_v, i_endidx_v
        jk=jk+1
    !------------------------------------------------------------------------------------------
    ! Also, for some reason, ice free drift used to be set here. Not sure what was the point...
    !------------------------------------------------------------------------------------------
!          ! Set the ice speed to free drift speed where concentration is less than 0.01
!          IF ( a_ice(jk) <= 0.01_wp ) THEN
!            u_ice(jk) = 0._wp; v_ice(jk) = 0._wp; u_change = 0._wp
!            ! TODO: Change u_change to speed_change
!            DO WHILE ( u_change > 1e-6_wp )
!              u_change = SQRT(u_ice(jk)**2+v_ice(jk)**2)
!              delu = SQRT( (u_w(jk)-u_ice(jk))**2 + (v_w(jk)-v_ice(jk))**2 )
!              u_ice(jk) = stress_atmice_x(jk)/( Cd_io*rho_ref*delu )
!              v_ice(jk) = stress_atmice_y(jk)/( Cd_io*rho_ref*delu )
!              u_change = ABS(u_change-SQRT(u_ice(jk)**2+v_ice(jk)**2))
!            ENDDO
!          ELSE
!            ! Strictly speaking p_ice%u_prog and p_ice%v_prog are only used by the restart files
!            ! now, so this does not need to be done every timestep, only after restart file is
!            ! read.
            u_ice(jk) = p_ice%u_prog(jv,jb)
            v_ice(jk) = p_ice%v_prog(jv,jb)
!          ENDIF
      END DO
    END DO

  END SUBROUTINE ice_fem_init_vel

  !-------------------------------------------------------------------------
  !
  !> Update p_ice%u_prog with last u_ice, v_ice values before writing a restart file
  !!
  !! @par Revision History
  !! Developed by Einar Olason, MPI-M (2013-08-05)
  !
  SUBROUTINE ice_fem_update_vel_for_restart(p_patch, p_ice)
    USE mo_ice,                 ONLY: u_ice, v_ice
    USE mo_sea_ice_types,       ONLY: t_sea_ice

    TYPE(t_patch), TARGET, INTENT(in)    :: p_patch
    TYPE(t_sea_ice),       INTENT(inout)    :: p_ice

    ! Local variables
    ! Patch and ranges
    TYPE(t_subset_range), POINTER :: all_verts
    ! Indexing
    INTEGER :: i_startidx_v, i_endidx_v
    INTEGER :: jk, jb, jv

  !-------------------------------------------------------------------------

    all_verts => p_patch%verts%all

    jk=0
    DO jb = all_verts%start_block, all_verts%end_block
      CALL get_index_range(all_verts, jb, i_startidx_v, i_endidx_v)
      DO jv = i_startidx_v, i_endidx_v
        jk=jk+1
!        ! Strictly speaking p_ice%u_prog and p_ice%v_prog are only used by the restart files now,
!        ! so this does not need to be done every timestep, only before restart file is written.
        p_ice%u_prog(jv,jb) = u_ice(jk)
        p_ice%v_prog(jv,jb) = v_ice(jk)

      END DO
    END DO

  END SUBROUTINE ice_fem_update_vel_for_restart

END MODULE mo_ice_fem_init