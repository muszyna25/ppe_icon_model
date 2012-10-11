!>
!!               The module <i>mo_model_import_domain</i>.
!!
!!               The module <i>mo_model_import_domain</i>
!! provides functionality to import information about the models computational
!! domain. This information is read from several files that were generated by
!! the patch generator programm. The data types describing the model domain are
!! contained in <i>mo_domain_model</i>.
!!
!! @par Revision History
!! Initial version  by: Peter Korn,  MPI-M, Hamburg, June 2005
!! Modification by Thomas Heinze (2006-02-21):
!! - renamed m_modules to mo_modules
!! Modification by Thomas Heinze (2006-09-20):
!! - added routine grid_and_patch_diagnosis
!! Modification by Pilar Ripodas, DWD, (2007-01-31)
!! - addapted to the new TYPE grid_edges (system_orientation added)
!! Modification by Peter Korn,  MPI-M, (2006-12)
!! - implementation of topography and boundary treatment, i.e.
!!   initialization of the grid & patch components that carry
!!   information about topography and the lateral boundaries of
!!   the domain; this is not related to patch boundaries.
!!   topography can either be computed by analytical l,eans or
!!   by reading from database files.
!! Modification by Hui Wan, MPI-M, (2007-02-23)
!! - Subroutine <i>init_import</i> was changed to <i>setup_grid</i>.
!!   Namelist hierarchy_ini was renamed to <i>grid_ctl</i>, moved
!!   from <i>mo_io_utilities</i> to this module and now read from
!!   an external file in subroutine <i>setup_grid</i>.
!! - Some changes in <i>init_ocean_patch_component</i> after
!!   discussion with Peter.
!! - Calculation of the min. primal edge length was added to
!!   <i>import_patches</i>. However, shouldn't it be an array with
!!   one element for each patch, rather than a scalar?
!! Modification by P. Ripodas, DWD, (2007-03-14):
!! - Now the output of "import_patches" is the min_dual_edge_lenght
!!   instead of the min_primal_edge_lenght. It will be used to set
!!   the horizontal diffusion parameter. Now it is done as it was
!!   in the prototype.
!! Modification by Almut Gassmann, MPI-M (2007-04)
!! - removed loptimize to make compatible with new grid generator
!! - removed itoa for good programming style
!! - reorganized patch input to be compatible with the new patch generator
!! - cleaning up "destruct_patches"
!! Modification by Almut Gassmann, MPI-M (2007-04-13)
!! - remove grid type and perform related adaptations
!!   (grid information comes now inside a patch)
!! - changed subroutine name form setup_grid to setup_files
!! Modified by Hui Wan, MPI-M, (2008-04-04)
!!  - topography_file_dir renamed topo_file_dir
!!  - for the hydro_atmos, control variable testtype renamed ctest_name.
!! Modified by Almut Gassmann, MPI-M, (2008-04-23)
!!  - itopo distinguishes now shallow water (itopo=1) orography function
!!    from hydro_atmos orography function (itopo=2)
!! Modification by Jochen Foerstner, DWD, (2008-07-16)
!!  - new fields in the derived type for the edges:
!!    grid_edges%primal_cart_normal (Cartesian normal to edge),
!!    grid_edges%quad_idx, grid_edges%quad_area and grid_edges%quad_orientation
!!    (indices of edges and area of the quadrilateral formed by two adjacent cells)
!!    up to now these new fields are initialized in the new routines
!!    calculate_primal_cart_normal and init_quad_twoadjcells
!!    rather than read from a grid/patch file.
!! Modification by Almut Gassmann, MPI-M, (2008-09-21)
!!  - remove reference to mask and height files, they are never used
!!  - use cell_type to distinguish cells as triangles or hexagons
!! Modification by Almut Gassmann, MPI-M (2008-10-30)
!!  - add subroutine init_coriolis to initialize Coriolis parameter
!! Modification by Constantin Junk, MPI-M (2011-04-05)
!! - ...
!! Modification by Daniel Reinert, DWD (2012-04-11)
!!  - new routine which initializes the butterfly data structure
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
!!

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_model_domimp_setup
  !-------------------------------------------------------------------------
  USE mo_kind,               ONLY: wp
  USE mo_exception,          ONLY: finish, warning
  USE mo_model_domain,       ONLY: t_patch
!   USE mo_physical_constants, ONLY: earth_angular_velocity
  USE mo_parallel_config,    ONLY: nproma
  USE mo_math_utilities,     ONLY: gvec2cvec, t_cartesian_coordinates
  USE mo_math_constants,     ONLY: pi_2
  USE mo_loopindices,        ONLY: get_indices_e
  USE mo_grid_config,        ONLY: corio_lat, grid_angular_velocity
  USE mo_sync,               ONLY: sync_c, sync_e, sync_patch_array, sync_idx
  USE mo_grid_subset,        ONLY: fill_subset,t_subset_range, get_index_range
  USE mo_mpi,                ONLY: work_mpi_barrier, get_my_mpi_work_id, my_process_is_mpi_seq
  USE mo_impl_constants,     ONLY: halo_levels_ceiling
  
  IMPLICIT NONE
  
  PRIVATE
  
  PUBLIC :: reshape_int
  PUBLIC :: reshape_real
  PUBLIC :: calculate_cart_normal
  PUBLIC :: init_quad_twoadjcells
  PUBLIC :: init_coriolis
  PUBLIC :: set_verts_phys_id
  PUBLIC :: init_butterfly_idx
  
  PUBLIC :: fill_grid_subsets
  
  CHARACTER(LEN=*), PARAMETER :: version = '$Id$'
  
  
  !-------------------------------------------------------------------------
  
CONTAINS

  !-------------------------------------------------------------------------
  !>
  !! This routine calculates the edge area of the patch
  !!
  SUBROUTINE calculate_edge_area(ptr_patch)
    TYPE(t_patch), TARGET, INTENT(inout) :: ptr_patch

    INTEGER                       :: je, jb, istart_e, iend_e
    TYPE(t_subset_range), POINTER :: all_edges
    !-------------------------------------------------------------------------
    all_edges => ptr_patch%edges%all

    ! a) the control volume associated to each edge is defined as the
    ! quadrilateral whose edges are the primal edge and the associated dual edge
    !----------------------------------------------------------------------------
    ! loop over all blocks and edges

!$OMP PARALLEL DO PRIVATE(jb,je) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = all_edges%start_block, all_edges%end_block
      CALL get_index_range(all_edges, jb, istart_e, iend_e)
      DO je = istart_e, iend_e
        ptr_patch%edges%area_edge(je,jb) =  &
            &    ptr_patch%edges%primal_edge_length(je,jb)  &
            &  * ptr_patch%edges%dual_edge_length(je,jb)
      END DO
    END DO
!$OMP END PARALLEL DO
  END SUBROUTINE calculate_edge_area
  
  !-------------------------------------------------------------------------
  !>
  !! This routine calculates the Cartesian components of the edge primal normals.
  !!
  !! Therefore the given zonal and meridional
  !! components of the edge primal normals are used.
  !! This information should be provided by the grid generator and
  !! read in the routine read_patch.
  !!
  !! @par Revision History
  !! Initial release by Jochen Foerstner (2008-05-19)
  !! Modifiaction by A. Gassmann(2010-09-05)
  !! - added also tangential normal, and generalize to lplane
  !!
  SUBROUTINE calculate_cart_normal( lplane, p_patch )
    !
    LOGICAL,INTENT(in) :: lplane
    TYPE(t_patch), TARGET, INTENT(inout) :: p_patch  ! patch on a specific level
    
    TYPE(t_cartesian_coordinates) :: z_vec
    REAL(wp) :: z_lon, z_lat, z_u, z_v  ! location and components of normal
    REAL(wp) :: z_norm                  ! norm of Cartesian normal
    
    INTEGER :: nlen, nblks_e, npromz_e
    INTEGER :: jb, je                   ! loop indices
    
    !-----------------------------------------------------------------------
    
    ! values for the blocking
    nblks_e  = p_patch%nblks_e
    npromz_e = p_patch%npromz_e
    
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,je,nlen,z_lon,z_lat,z_u,z_v,z_norm,z_vec) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = 1, nblks_e
      
      IF (jb /= nblks_e) THEN
        nlen = nproma
      ELSE
        nlen = npromz_e
      ENDIF
      
      ! loop over edges
      DO je = 1, nlen
        
        ! location of edge midpoint
        IF (.NOT.lplane) THEN
          z_lon = p_patch%edges%center(je,jb)%lon
          z_lat = p_patch%edges%center(je,jb)%lat
        ELSE
          z_lon = -pi_2 !-90
          z_lat =  pi_2 !+90
        ENDIF
        
        ! zonal and meridional component of primal normal
        z_u = p_patch%edges%primal_normal(je,jb)%v1
        z_v = p_patch%edges%primal_normal(je,jb)%v2
        
        ! calculate Cartesian components of primal normal
        CALL gvec2cvec( z_u, z_v, z_lon, z_lat, z_vec%x(1), z_vec%x(2), z_vec%x(3) )
        
        ! compute unit normal to edge je
        z_norm = SQRT( DOT_PRODUCT(z_vec%x(1:3),z_vec%x(1:3)) )
        z_vec%x(1:3) = 1._wp / z_norm * z_vec%x(1:3)
        
        ! save the values in the according type structure of the patch
        p_patch%edges%primal_cart_normal(je,jb)%x(1) = z_vec%x(1)
        p_patch%edges%primal_cart_normal(je,jb)%x(2) = z_vec%x(2)
        p_patch%edges%primal_cart_normal(je,jb)%x(3) = z_vec%x(3)
        
        ! zonal and meridional component of dual normal
        z_u = p_patch%edges%dual_normal(je,jb)%v1
        z_v = p_patch%edges%dual_normal(je,jb)%v2
        
        ! calculate Cartesian components of primal normal
        CALL gvec2cvec( z_u, z_v, z_lon, z_lat, z_vec%x(1), z_vec%x(2), z_vec%x(3) )
        
        ! compute unit normal to edge je
        z_norm = SQRT( DOT_PRODUCT(z_vec%x(1:3),z_vec%x(1:3)) )
        z_vec%x(1:3) = 1._wp / z_norm * z_vec%x(1:3)
        
        ! save the values in the according type structure of the patch
        p_patch%edges%dual_cart_normal(je,jb)%x(1) = z_vec%x(1)
        p_patch%edges%dual_cart_normal(je,jb)%x(2) = z_vec%x(2)
        p_patch%edges%dual_cart_normal(je,jb)%x(3) = z_vec%x(3)
        
      END DO
      
    END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
    
  END SUBROUTINE calculate_cart_normal
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  !>
  !! This routine initializes the data for the quadrilateral cells.
  !!
  !! This routine initializes the data for the quadrilateral cells
  !! formed by the two adjacent cells of an edge.
  !! The four edge indices of this qudrilateral and its area are stored
  !! in the derived type for the edges.
  !! This information should be provided by the grid generator and
  !! read in the routine read_patch.
  !!
  !! @par Revision History
  !! Initial release by Jochen Foerstner (2008-07-16)
  !! Modification by Guenther Zaengl (2009-05-08): vectorization
  !! Also suitable for hexagons: Almut Gassmann (2009-10-01)
  !!
  SUBROUTINE init_quad_twoadjcells( p_patch )
    !
    TYPE(t_patch), TARGET, INTENT(inout) :: p_patch  ! patch on a specific level
    
    INTEGER :: nblks_e
    INTEGER :: je, jb                   ! loop index
    INTEGER :: iie, iie1(nproma), iie2(nproma), ierror(nproma)
    INTEGER :: ilc1, ibc1, ilc2, ibc2                   ! cell indices
    INTEGER :: ile1, ibe1, ile2, ibe2, ile3, ibe3       ! edge indices
    INTEGER :: ilv1, ibv1, ilv2, ibv2, ilev, ibev, jjev ! vertex indices
    INTEGER :: i_startblk, i_startidx, i_endidx
    
    
    !-----------------------------------------------------------------------
    
    ! values for the blocking
    nblks_e  = p_patch%nblks_e
    
    ! Quad cells cannot be computed along the lateral edges of a nested domain
    i_startblk = p_patch%edges%start_blk(2,1)
    
!$OMP PARALLEL
    
    ! Initialize array elements along nest boundaries with zero
    IF (p_patch%id > 1) THEN
!$OMP WORKSHARE
      p_patch%edges%quad_area(:,1:i_startblk)            = 0._wp
      p_patch%edges%quad_idx(:,1:i_startblk,1:4)         = 0
      p_patch%edges%quad_blk(:,1:i_startblk,1:4)         = 0
      p_patch%edges%quad_orientation(:,1:i_startblk,1:4) = 0._wp
!$OMP END WORKSHARE
    ENDIF
    
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,je,iie,ilc1,ibc1,ilc2,ibc2,ile1,ibe1,ile2,ibe2,&
!$OMP            ile3,ibe3,iie1,iie2,ierror,ilv1,ibv1,ilv2,ibv2,ilev,ibev,&
!$OMP jjev) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, nblks_e
      
      CALL get_indices_e(p_patch, jb, i_startblk, nblks_e, &
        & i_startidx, i_endidx, 2)
      
      IF (p_patch%cell_type == 3) THEN ! vectorized code for triangular grid
        
        ierror = 0
        
        DO je = i_startidx, i_endidx
          
          IF(.NOT.p_patch%edges%owner_mask(je,jb)) CYCLE
          !
          ! get global indices of the edges of the two neighboring cells
          !
          ilc1 = p_patch%edges%cell_idx(je,jb,1)
          ibc1 = p_patch%edges%cell_blk(je,jb,1)
          ilc2 = p_patch%edges%cell_idx(je,jb,2)
          ibc2 = p_patch%edges%cell_blk(je,jb,2)
          
          ! sum area of the two adjacent cells
          p_patch%edges%quad_area(je,jb) = &
            & p_patch%cells%area(ilc1,ibc1) + p_patch%cells%area(ilc2,ibc2)
          
          ! Indices of edges adjacent to cell 1
          ile1 = p_patch%cells%edge_idx(ilc1,ibc1,1)
          ibe1 = p_patch%cells%edge_blk(ilc1,ibc1,1)
          ile2 = p_patch%cells%edge_idx(ilc1,ibc1,2)
          ibe2 = p_patch%cells%edge_blk(ilc1,ibc1,2)
          ile3 = p_patch%cells%edge_idx(ilc1,ibc1,3)
          ibe3 = p_patch%cells%edge_blk(ilc1,ibc1,3)
          
          IF (je == ile1 .AND. jb == ibe1) THEN
            iie1(je) = 2
            iie2(je) = 3
          ELSE IF (je == ile2 .AND. jb == ibe2) THEN
            iie1(je) = 3
            iie2(je) = 1
          ELSE IF (je == ile3 .AND. jb == ibe3) THEN
            iie1(je) = 1
            iie2(je) = 2
          ELSE
            ierror(je) = ierror(je) + 1
          ENDIF
          
          p_patch%edges%quad_idx(je,jb,1) = p_patch%cells%edge_idx(ilc1,ibc1,iie1(je))
          p_patch%edges%quad_blk(je,jb,1) = p_patch%cells%edge_blk(ilc1,ibc1,iie1(je))
          p_patch%edges%quad_orientation(je,jb,1) =  &
            & p_patch%cells%edge_orientation(ilc1,ibc1,iie1(je))
          
          p_patch%edges%quad_idx(je,jb,2) = p_patch%cells%edge_idx(ilc1,ibc1,iie2(je))
          p_patch%edges%quad_blk(je,jb,2) = p_patch%cells%edge_blk(ilc1,ibc1,iie2(je))
          p_patch%edges%quad_orientation(je,jb,2) =  &
            & p_patch%cells%edge_orientation(ilc1,ibc1,iie2(je))
          
          ! Indices of edges adjacent to cell 2
          ile1 = p_patch%cells%edge_idx(ilc2,ibc2,1)
          ibe1 = p_patch%cells%edge_blk(ilc2,ibc2,1)
          ile2 = p_patch%cells%edge_idx(ilc2,ibc2,2)
          ibe2 = p_patch%cells%edge_blk(ilc2,ibc2,2)
          ile3 = p_patch%cells%edge_idx(ilc2,ibc2,3)
          ibe3 = p_patch%cells%edge_blk(ilc2,ibc2,3)
          
          IF (je == ile1 .AND. jb == ibe1) THEN
            iie1(je) = 2
            iie2(je) = 3
          ELSE IF (je == ile2 .AND. jb == ibe2) THEN
            iie1(je) = 3
            iie2(je) = 1
          ELSE IF (je == ile3 .AND. jb == ibe3) THEN
            iie1(je) = 1
            iie2(je) = 2
          ELSE
            ierror(je) = ierror(je) + 1
          ENDIF
          
          p_patch%edges%quad_idx(je,jb,3) = p_patch%cells%edge_idx(ilc2,ibc2,iie1(je))
          p_patch%edges%quad_blk(je,jb,3) = p_patch%cells%edge_blk(ilc2,ibc2,iie1(je))
          p_patch%edges%quad_orientation(je,jb,3) =  &
            & p_patch%cells%edge_orientation(ilc2,ibc2,iie1(je))
          
          p_patch%edges%quad_idx(je,jb,4) = p_patch%cells%edge_idx(ilc2,ibc2,iie2(je))
          p_patch%edges%quad_blk(je,jb,4) = p_patch%cells%edge_blk(ilc2,ibc2,iie2(je))
          p_patch%edges%quad_orientation(je,jb,4) =  &
            & p_patch%cells%edge_orientation(ilc2,ibc2,iie2(je))
          
        ENDDO
        IF (MAXVAL(ierror) > 0) THEN
          WRITE(0,*) "Number of errors", SUM(ierror(1:nproma))
          CALL finish ('mo_model_domain_import:init_quad_twoadjcells',  &
            & 'edge-cell index relationships are apparently incorrect')
        ENDIF
        
      ELSE ! cell_type == 6
        
        DO je = i_startidx, i_endidx
          
          IF(.NOT.p_patch%edges%owner_mask(je,jb)) CYCLE
          
          iie = 0
          !
          ! get global indices of the edges of the two neighboring verts
          !
          ilv1 = p_patch%edges%vertex_idx(je,jb,1)
          ibv1 = p_patch%edges%vertex_blk(je,jb,1)
          ilv2 = p_patch%edges%vertex_idx(je,jb,2)
          ibv2 = p_patch%edges%vertex_blk(je,jb,2)
          
          ! sum area of the two adjacent verts
          p_patch%edges%quad_area(je,jb) = &
            & p_patch%verts%dual_area(ilv1,ibv1) + p_patch%verts%dual_area(ilv2,ibv2)
          
          DO jjev = 1, 3
            
            ilev = p_patch%verts%edge_idx(ilv1,ibv1,jjev)
            ibev = p_patch%verts%edge_blk(ilv1,ibv1,jjev)
            
            ! test if edge is not the current one, then store line and
            ! block indices and orientation of quad edge
            IF ( ilev /= je .OR. ibev /= jb ) THEN
              iie = iie+1
              p_patch%edges%quad_idx(je,jb,iie) = ilev
              p_patch%edges%quad_blk(je,jb,iie) = ibev
              p_patch%edges%quad_orientation(je,jb,iie) =  &
                & p_patch%verts%edge_orientation(ilv1,ibv1,jjev)
            ENDIF
            
          ENDDO
          
          DO jjev = 1, 3
            
            ilev = p_patch%verts%edge_idx(ilv2,ibv2,jjev)
            ibev = p_patch%verts%edge_blk(ilv2,ibv2,jjev)
            
            ! test if edge is not the current one, then store line and
            ! block indices and orientation of quad edge
            IF ( ilev /= je .OR. ibev /= jb ) THEN
              iie = iie+1
              p_patch%edges%quad_idx(je,jb,iie) = ilev
              p_patch%edges%quad_blk(je,jb,iie) = ibev
              p_patch%edges%quad_orientation(je,jb,iie) =  &
                & p_patch%verts%edge_orientation(ilv2,ibv2,jjev)
            ENDIF
            
          ENDDO
          
          IF (iie /= 4) THEN
            PRINT *, "ERROR ==>  iie = ", iie,  " /= 4 "
            CALL finish ('mo_model_domain_import:init_quad_twoadjcells',  &
              & 'wrong number of edge indices for quad')
          ENDIF
          
        END DO ! end edge loop
        
      ENDIF ! cell_type
      
    END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
    
    CALL sync_patch_array(sync_e,p_patch,p_patch%edges%quad_area)
    DO iie = 1, 4
      CALL sync_patch_array(sync_e,p_patch,p_patch%edges%quad_orientation(:,:,iie))
      
      CALL sync_idx(sync_e,sync_e,p_patch,p_patch%edges%quad_idx(:,:,iie), &
        & p_patch%edges%quad_blk(:,:,iie), &
        & opt_remap=.FALSE.)
    ENDDO
    
    
  END SUBROUTINE init_quad_twoadjcells
  !-------------------------------------------------------------------------
  

  !-------------------------------------------------------------------------
  !>
  !! This routine initializes the butterfly index field.
  !!
  !! This routine initializes the butterfly index field, formed by the 
  !! 4 cells sharing the 2 vertices which bound a given edge.
  !! The 4 cell indices are stored in the derived type for the edges.
  !!
  !!     -------------------------------
  !!       \          /\          /
  !!        \ (2,2)  /  \  (2,1) /             (edge_neighbor,cell)
  !!         \      /    \      /
  !!          \    /      \    /
  !!           \  /   2    \  /             ^  : edge normal
  !!            \/          \/             /|\
  !!       ----2-============-1---------    |
  !!            /\          /\
  !!           /  \   1    /  \
  !!          /    \      /    \
  !!         /      \    /      \
  !!        /  (1,2) \  /  (1,1) \
  !!       /          \/          \
  !!     -------------------------------
  !!
  !!
  !! Storage order:  butterfly_idx(nproma,nblks_c,edge_neighbor,cell)
  !!                 butterfly_blk(nproma,nblks_c,edge_neighbor,cell)
  !!
  !! The cells are subdivided into two classes: Cells which are 
  !! neighbors of edge-neighbor 1 and cells which are neighbors of edge-
  !! neighbor 2. These cells are then numbered according to the number 
  !! of the edge-vertex they share. (see ASCII-ART)
  !!
  !! @par Revision History
  !! Initial release by Daniel Reinert (2012-04-05)
  !!
  SUBROUTINE init_butterfly_idx( p_patch )
    !
    TYPE(t_patch), TARGET, INTENT(inout) :: p_patch  ! patch on a specific level

    INTEGER :: cnt               ! counter
    INTEGER :: nblks_e
    INTEGER :: rl_start
    INTEGER :: je, jen, jn, jv, jb                   ! loop index
    INTEGER :: iln, ibn                              ! cell indices (neighbors)
    INTEGER :: iie, ije
    INTEGER :: il_bf1(2), ib_bf1(2)
    INTEGER :: il_bf2(2), ib_bf2(2)
    INTEGER :: ilc1, ibc1, ilc2, ibc2                ! cell indices
    INTEGER :: ilv1, ibv1, ilv2, ibv2, ilv, ibv      ! vertex indices
    INTEGER :: i_startblk, i_startidx, i_endidx
    
    
    !-----------------------------------------------------------------------
    
    ! values for the blocking
    nblks_e  = p_patch%nblks_e
    
    rl_start = 3
    i_startblk = p_patch%edges%start_blk(rl_start,1)
    
!$OMP PARALLEL
    
    ! Initialize array elements along nest boundaries with zero
    IF (p_patch%id > 1) THEN
!$OMP WORKSHARE
      p_patch%edges%butterfly_idx(:,1:i_startblk,:,:)  = 0
      p_patch%edges%butterfly_blk(:,1:i_startblk,:,:)  = 0
!$OMP END WORKSHARE
    ENDIF
    
!$OMP DO PRIVATE(jb,je,i_startidx,i_endidx,ilc1,ilc2,ibc1,ibc2, &
!$OMP            cnt,jen,iln,ibn,il_bf1,ib_bf1,il_bf2,ib_bf2,   &
!$OMP            ilv1,ilv2,ibv1,ibv2,jn,jv,ilv,ibv)
    DO jb = i_startblk, nblks_e
      
      CALL get_indices_e(p_patch, jb, i_startblk, nblks_e, &
        & i_startidx, i_endidx, rl_start)
      
        
      DO je = i_startidx, i_endidx
          
        IF(.NOT.p_patch%edges%owner_mask(je,jb)) CYCLE

        !
        ! get global indices of the two neighboring cells
        !
        ilc1 = p_patch%edges%cell_idx(je,jb,1)
        ibc1 = p_patch%edges%cell_blk(je,jb,1)
        ilc2 = p_patch%edges%cell_idx(je,jb,2)
        ibc2 = p_patch%edges%cell_blk(je,jb,2)
          

        ! get global indices of the two cells neighboring 
        ! each edge-neighbor, while excluding the other 
        ! edge neighbor.
        !
        ! neighbor 1 (ilc1,ibc1)
        cnt = 0
        DO jen = 1,3
          iln= p_patch%cells%neighbor_idx(ilc1,ibc1,jen)
          ibn= p_patch%cells%neighbor_blk(ilc1,ibc1,jen)
          IF ((iln == ilc2) .AND. (ibn == ibc2)) CYCLE
          cnt = cnt + 1
          il_bf1(cnt) = iln
          ib_bf1(cnt) = ibn
        ENDDO


        ! neighbor 2 (ilc2,ibc2)
        cnt = 0
        DO jen = 1,3
          iln= p_patch%cells%neighbor_idx(ilc2,ibc2,jen)
          ibn= p_patch%cells%neighbor_blk(ilc2,ibc2,jen)
          IF ((iln == ilc1) .AND. (ibn == ibc1)) CYCLE
          cnt = cnt + 1
          il_bf2(cnt) = iln
          ib_bf2(cnt) = ibn
        ENDDO


        ! Now, the global indices of the 4 cells are known, and 
        ! they are already distinguished with respect to the edge 
        ! neighbor. However, we do not know yet, whether they share 
        ! vertex 1 or 2.

        ! Indices of vertices bounding edge je
        !
        ilv1= p_patch%edges%vertex_idx(je,jb,1)
        ibv1= p_patch%edges%vertex_blk(je,jb,1)
        ilv2= p_patch%edges%vertex_idx(je,jb,2)
        ibv2= p_patch%edges%vertex_blk(je,jb,2)


        ! check neighbors of neighbor 1, whether they share 
        ! vertex 1 or 2
        DO jn= 1, 2   ! loop over 2 neighbors of neighbor 1
          !
          ! get vertex indices of neighbor jn
          DO jv=1,3
            ilv = p_patch%cells%vertex_idx(il_bf1(jn),ib_bf1(jn),jv)
            ibv = p_patch%cells%vertex_blk(il_bf1(jn),ib_bf1(jn),jv)

            IF ( ilv == ilv1 .AND. ibv == ibv1 ) THEN
              p_patch%edges%butterfly_idx(je,jb,1,1) = il_bf1(jn)
              p_patch%edges%butterfly_blk(je,jb,1,1) = ib_bf1(jn)
            ELSE IF ( ilv == ilv2 .AND. ibv == ibv2 ) THEN
              p_patch%edges%butterfly_idx(je,jb,1,2) = il_bf1(jn)
              p_patch%edges%butterfly_blk(je,jb,1,2) = ib_bf1(jn)
            ENDIF
          ENDDO
        ENDDO

        ! check neighbors of neighbor 2, whether they share 
        ! vertex 1 or 2
        DO jn= 1, 2   ! loop over 2 neighbors of neighbor 2
          !
          ! get vertex indices of neighbor jn
          DO jv=1,3
            ilv = p_patch%cells%vertex_idx(il_bf2(jn),ib_bf2(jn),jv)
            ibv = p_patch%cells%vertex_blk(il_bf2(jn),ib_bf2(jn),jv)

            IF ( ilv == ilv1 .AND. ibv == ibv1 ) THEN
              p_patch%edges%butterfly_idx(je,jb,2,1) = il_bf2(jn)
              p_patch%edges%butterfly_blk(je,jb,2,1) = ib_bf2(jn)
            ELSE IF ( ilv == ilv2 .AND. ibv == ibv2 ) THEN
              p_patch%edges%butterfly_idx(je,jb,2,2) = il_bf2(jn)
              p_patch%edges%butterfly_blk(je,jb,2,2) = ib_bf2(jn)
            ENDIF
          ENDDO
        ENDDO       

      END DO  ! je      
    END DO  ! jb
!$OMP END DO
!$OMP END PARALLEL

    DO iie = 1, 2
      DO ije = 1, 2
        CALL sync_idx(sync_e,sync_c,p_patch,            &
          &  p_patch%edges%butterfly_idx(:,:,iie,ije),  &
          &  p_patch%edges%butterfly_blk(:,:,iie,ije),  &
          &  opt_remap=.FALSE.)
      ENDDO
    ENDDO
    
    
  END SUBROUTINE init_butterfly_idx



  !-------------------------------------------------------------------------
  !>
  !! Initializes the Coriolis components of the grid with analytical values.
  !!
  !! @par Revision History
  !! Developed by Almut Gassmann, MPI-M (2008-10-30).
  !! Modification by Daniel Reinert, DWD (2010-07-13),
  !! - Allocation of f_c, f_e, f_v has been moved to read_patch
  !!
  SUBROUTINE init_coriolis( lcorio, lplane, p_patch )
    
    LOGICAL,INTENT(in) :: lcorio, lplane
    TYPE(t_patch), TARGET, INTENT(inout) :: p_patch ! patch on specific level
    
    INTEGER :: nlen,                         &
      & nblks_c,  nblks_e,  nblks_v,  &
      & npromz_c, npromz_e, npromz_v
    INTEGER :: jc, je, jv, jb
    
    REAL(wp) :: zlat
    
    !-----------------------------------------------------------------------
    
    ! values for the blocking
    nblks_c  = p_patch%nblks_c
    npromz_c = p_patch%npromz_c
    nblks_e  = p_patch%nblks_e
    npromz_e = p_patch%npromz_e
    nblks_v  = p_patch%nblks_v
    npromz_v = p_patch%npromz_v
    
    IF (lcorio .AND. .NOT. lplane) THEN
      
      DO jb = 1, nblks_c
        IF (jb /= nblks_c) THEN
          nlen = nproma
        ELSE
          nlen = npromz_c
        ENDIF
        DO jc = 1, nlen
          zlat = p_patch%cells%center(jc,jb)%lat
          p_patch%cells%f_c(jc,jb) = 2._wp * grid_angular_velocity * SIN(zlat)
        END DO
      END DO
      
      DO jb = 1, nblks_e
        IF (jb /= nblks_e) THEN
          nlen = nproma
        ELSE
          nlen = npromz_e
        ENDIF
        DO je = 1, nlen
          zlat = p_patch%edges%center(je,jb)%lat
          p_patch%edges%f_e(je,jb) = 2._wp * grid_angular_velocity * SIN(zlat)
        END DO
      END DO
      
      DO jb = 1, nblks_v
        IF (jb /= nblks_v) THEN
          nlen = nproma
        ELSE
          nlen = npromz_v
        ENDIF
        DO jv = 1, nlen
          zlat = p_patch%verts%vertex(jv,jb)%lat
          p_patch%verts%f_v(jv,jb) = 2._wp * grid_angular_velocity * SIN(zlat)
        END DO
      END DO
      
    ELSEIF (lcorio .AND. lplane) THEN
      
      p_patch%cells%f_c(:,:) = 2._wp*grid_angular_velocity*SIN(corio_lat)
      p_patch%edges%f_e(:,:) = 2._wp*grid_angular_velocity*SIN(corio_lat)
      p_patch%verts%f_v(:,:) = 2._wp*grid_angular_velocity*SIN(corio_lat)
      
    ELSE
      grid_angular_velocity = 0.0_wp
      p_patch%cells%f_c(:,:) = 0.0_wp
      p_patch%edges%f_e(:,:) = 0.0_wp
      p_patch%verts%f_v(:,:) = 0.0_wp
      
    ENDIF
    
  END SUBROUTINE init_coriolis
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  !>
  !! Sets phys_id for verts since this is not read from input
  !!
  !! @par Revision History
  !! Developed by Rainer Johanni (2010-12-05).
  !!
  SUBROUTINE set_verts_phys_id( p_patch )
    
    TYPE(t_patch), INTENT(inout) :: p_patch ! patch on specific level
    
    INTEGER :: nlen, jb, jv, ji, ilc, ibc
    
    !-----------------------------------------------------------------------
    
    IF (p_patch%cell_type == 3) THEN
      
      DO jb = 1, p_patch%nblks_v
        IF (jb /= p_patch%nblks_v) THEN
          nlen = nproma
        ELSE
          nlen = p_patch%npromz_v
        ENDIF
        DO jv = 1, nlen
          ! Just go over all neighboring cells of a vertex and take the phys_id
          ! of the first valid neighbor cell
          DO ji = 1, 6
            ilc = p_patch%verts%cell_idx(jv,jb,ji)
            ibc = p_patch%verts%cell_blk(jv,jb,ji)
            IF(ilc > 0) THEN
              p_patch%verts%phys_id(jv,jb) = p_patch%cells%phys_id(ilc, ibc)
              EXIT
            ENDIF
          END DO
        END DO
      END DO
      
    ELSE
      
      ! There is no distinction between physical and logical patches for hex grids
      p_patch%verts%phys_id(:,:) = p_patch%id
      
    ENDIF
    
  END SUBROUTINE set_verts_phys_id
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  !>
  !!               reshape an integer array for the blocking.
  !!
  !! @par Revision History
  !! Developed  by  Jochen Foerstner, DWD (2008-10-21)
  !!
  SUBROUTINE reshape_int( p_int_array_in, nblks, npromz,  &
    & p_reshaped_int_array_out )
    
    
    ! input array
    INTEGER, INTENT(in) :: p_int_array_in(:)
    ! number of blocks
    INTEGER, INTENT(in) :: nblks
    ! chunk length
    INTEGER, INTENT(in) :: npromz
    
    ! output array
    INTEGER, INTENT(inout) :: p_reshaped_int_array_out(:,:)
    
    INTEGER :: nlen
    INTEGER :: jl, jb
    INTEGER :: il
    
    !-----------------------------------------------------------------------
    
    DO jb = 1, nblks
      
      IF (jb /= nblks) THEN
        nlen = nproma
      ELSE
        nlen = npromz
        DO jl = npromz+1, nproma
          p_reshaped_int_array_out(jl,nblks) = 0
        END DO
      END IF
      
      DO jl = 1, nlen
        il = jl + (jb-1)*nproma
        p_reshaped_int_array_out(jl,jb) = p_int_array_in(il)
      END DO
      
    END DO
    
  END SUBROUTINE reshape_int
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  !>
  !!               reshape a real array for the blocking.
  !!
  !!
  !! @par Revision History
  !! Developed  by  Jochen Foerstner, DWD (2008-10-21)
  !!
  SUBROUTINE reshape_real( p_real_array_in, nblks, npromz,  &
    & p_reshaped_real_array_out )
    
    ! input array
    REAL(wp), INTENT(in):: p_real_array_in(:)
    ! number of blocks
    INTEGER, INTENT(in) :: nblks
    ! chunk length
    INTEGER, INTENT(in) :: npromz
    
    ! output array
    REAL(wp), INTENT(inout) :: p_reshaped_real_array_out(:,:)
    
    INTEGER :: nlen
    INTEGER :: jl, jb
    INTEGER :: il
    
    !-----------------------------------------------------------------------
    
    DO jb = 1, nblks
      
      IF (jb /= nblks) THEN
        nlen = nproma
      ELSE
        nlen = npromz
        DO jl = npromz+1, nproma
          p_reshaped_real_array_out(jl,nblks) = 0._wp
        END DO
      END IF
      
      DO jl = 1, nlen
        il = jl + (jb-1)*nproma
        p_reshaped_real_array_out(jl,jb) = p_real_array_in(il)
      END DO
      
    END DO
    
  END SUBROUTINE reshape_real
  !-------------------------------------------------------------------------
  
  !-----------------------------------------------------------------------------
  !>
  !! Fills the grid's subsets
  SUBROUTINE fill_grid_subsets(p_patch)
    TYPE(t_patch), INTENT(inout) :: p_patch
    
    INTEGER :: j, jb, jl, jg
    
    !--------------------------------------------------------------------------------
    ! aliasing the halo_level to decomp_domain
    p_patch%cells%halo_level => p_patch%cells%decomp_domain
    p_patch%edges%halo_level => p_patch%edges%decomp_domain
    p_patch%verts%halo_level => p_patch%verts%decomp_domain
    !--------------------------------------------------------------------------------
    ! make sure the levels are correct when running sequentially
    IF (my_process_is_mpi_seq()) THEN
      p_patch%cells%halo_level(:,:) = 0
      p_patch%cells%halo_level(p_patch%npromz_c + 1 :nproma, p_patch%nblks_c) = -1
      p_patch%edges%halo_level(:,:) = 0
      p_patch%edges%halo_level(p_patch%npromz_e + 1 :nproma, p_patch%nblks_e) = -1
      p_patch%verts%halo_level(:,:) = 0
      p_patch%verts%halo_level(p_patch%npromz_v + 1 :nproma, p_patch%nblks_v) = -1
    ENDIF
    !--------------------------------------------------------------------------------
    ! fill cell subsets
    !    fill_subset(range subset,  halo levels, start level, end level)
    CALL fill_subset(p_patch%cells%all, p_patch, p_patch%cells%halo_level, &
      & 0, halo_levels_ceiling)
    p_patch%cells%all%is_in_domain   = .false.
    
    CALL fill_subset(p_patch%cells%owned, p_patch, p_patch%cells%halo_level, 0, 0)
    p_patch%cells%owned%is_in_domain = .true.
    
    CALL fill_subset(p_patch%cells%in_domain, p_patch, p_patch%cells%halo_level, 0, 0)
    p_patch%cells%in_domain%is_in_domain = .true.
    
    CALL fill_subset(p_patch%cells%not_owned, p_patch, p_patch%cells%halo_level,&
      & 1, halo_levels_ceiling)
    p_patch%cells%not_owned%is_in_domain = .false.
    
    CALL fill_subset(p_patch%cells%not_in_domain,  p_patch, &
      & p_patch%cells%halo_level, 2, halo_levels_ceiling)
    p_patch%cells%not_in_domain%is_in_domain = .false.
    
    CALL fill_subset(p_patch%cells%one_edge_in_domain, p_patch, p_patch%cells%halo_level, 1, 1)
    p_patch%cells%one_edge_in_domain%is_in_domain = .false.
    
    IF (p_patch%cells%in_domain%no_of_holes > 0) &
      CALL warning("p_patch%cells%in_domain", "no_of_holes > 0")
    IF (p_patch%cells%one_edge_in_domain%no_of_holes > 0) &
      CALL warning("p_patch%cells%one_edge_in_domain", "no_of_holes > 0")
    
    ! fill edge subsets
    !    fill_subset(range subset,  halo levels, start level, end level)
    CALL fill_subset(p_patch%edges%all, p_patch, p_patch%edges%halo_level, 0, halo_levels_ceiling)
    p_patch%edges%all%is_in_domain   = .false.
    
    CALL fill_subset(p_patch%edges%owned, p_patch,p_patch%edges%halo_level, 0, 0)
    p_patch%edges%owned%is_in_domain = .true.
    
    CALL fill_subset(p_patch%edges%in_domain, p_patch, p_patch%edges%halo_level, 0, 1)
    p_patch%edges%in_domain%is_in_domain   = .true.
    
    CALL fill_subset(p_patch%edges%not_owned, p_patch, p_patch%edges%halo_level, &
      & 1, halo_levels_ceiling)
    p_patch%edges%not_owned%is_in_domain   = .false.
    
    CALL fill_subset(p_patch%edges%not_in_domain, p_patch, &
      & p_patch%edges%halo_level, 2, halo_levels_ceiling)
    p_patch%edges%not_in_domain%is_in_domain   = .false.
    
    ! fill vertex subsets
    !    fill_subset(range subset,  halo levels, start level, end level)
    CALL fill_subset(p_patch%verts%all,  p_patch, p_patch%verts%halo_level, &
      & 0, halo_levels_ceiling)
    p_patch%verts%all%is_in_domain   = .false.
    
    CALL fill_subset(p_patch%verts%owned, p_patch, p_patch%verts%halo_level, 0, 0)
    p_patch%verts%owned%is_in_domain   = .true.
    
    CALL fill_subset(p_patch%verts%in_domain, p_patch, p_patch%verts%halo_level, 0, 1)
    p_patch%verts%in_domain%is_in_domain   = .true.
    
    CALL fill_subset(p_patch%verts%not_owned, p_patch, p_patch%verts%halo_level, 1, halo_levels_ceiling)
    p_patch%verts%not_owned%is_in_domain   = .false.
    
    CALL fill_subset(p_patch%verts%not_in_domain, p_patch, &
      & p_patch%verts%halo_level, 2, halo_levels_ceiling)
    p_patch%verts%not_in_domain%is_in_domain   = .false.
    
    
      
    ! write some info:
    !     IF (get_my_mpi_work_id() == 1) CALL work_mpi_barrier()
    !     write(0,*) get_my_mpi_work_id(), "egdes global_index:", p_patch%edges%glb_index
    !     write(0,*) get_my_mpi_work_id(), "egdes halo_level:", p_patch%edges%halo_level
    !     write(0,*) get_my_mpi_work_id(), "egdes owner:", p_patch%edges%owner_local
    !     IF (get_my_mpi_work_id() == 0) CALL work_mpi_barrier()
    
  END SUBROUTINE fill_grid_subsets
  !-----------------------------------------------------------------------------
  
END MODULE mo_model_domimp_setup
