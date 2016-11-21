!>
!! @page pagegraphgeneratorf901 Main routines for new graph/grid/patch generator
!!
!! @author
!!     Luca Bonaventura, Luis Kornblueh, Peter Korn, Uwe Schulzweida,
!!     (MPI-M)
!!
!! @author
!!     Thomas Heinze, Pilar Ripodas, Guenther Zaengl
!!     (DWD)
!!
!! @author
!!     Peter Sanders
!!     (Uni Karlsruhe)
!!
!!
!! @date 2009-12-14 07:35:04 UTC
!!

!>
!!               This is the main program for the ICON grid generator.
!!
!!               Starting from the original icosahedron, a finer  <i>root</i>
!!               grid is produced, which is the basis for subsequent dyadic
!!               refinements that allow to reach any desired resolution.
!!               The grid generator uses an internal representation of the
!!               data that is defined in the module <i>mo_base_datatypes.</i>
!!               The grid construction proceeds in four steps:
!!
!!               <ul>
!!               <li>  a tree structure is defined by  <i>generate_tree</i>,
!!                     which allows to identify each triangular cell in the
!!                     grid hierarchy and its relationship to coarser cells
!!                     to which it belongs and to finer cells that belong to
!!                     it.
!!               <li>  neighbourhood relationships among cells, edges and
!!                     vertices are defined in <i>generate_graph</i>
!!               <li>  the geometric quantities are computed by
!!                     <i>generate_geometry</i>
!!               <li>  the computed quantities and the neighbourhood
!!                     relationships are stored onto a separate <i>grid</i>
!!                     structure, defined in <i>mo_grid</i> that is used by
!!                     the numerical subroutines.
!!               </ul>
!!
!!               The <i>grid</i> structure represents the data as one
!!               dimensional vectors and uses an addressing that inherits
!!               the properties of the internal tree representation, which
!!               ensures data locality.
!!
!!               The data are stored in graphmap netCDF files.
!!
!! @par Revision History
!!  Inital brainstorming by Luis Kornblueh, Peter Sanders, Uwe Schulzweida
!!  and Luca Bonaventura, (2004)
!! @par
!!  Development of the indexing procedure for triangular cells
!!  by Luis Kornblueh (2004).
!!  Development of the inherited full graph structure
!!  by Luca Bonaventura (2005).
!!  Modification by Thomas Heinze, DWD, (2007-05-04)
!!  - split grid generator in two parts, one for creating graph (this file)
!!    another for calculating the grid information (grid_generator)
!!  Almut Gassmann, MPI-M (2009-01-15)
!!  - introduce planar option
!!
!!  prepare_gridref developed by Guenther Zaengl, DWD (2009-07-22)
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
MODULE mo_global_grid_generator
  !-------------------------------------------------------------------------
  !
  !            Almut Gassmann
  !
  !-------------------------------------------------------------------------
  USE ISO_C_BINDING,         ONLY: C_DOUBLE

  USE mo_exception, ONLY: message, message_text, open_log, finish
  USE mo_topology,  ONLY: init_topology_graph,  init_topology_grid, &
    & generate_tree,          &
    & generate_graph,         &
    & write_graph,            &
    & destruct_topology_graph, destruct_topology_grid, read_graph
  USE mo_io_graph,  ONLY: init_graphgen
  USE mo_geometry,  ONLY: generate_geometry
  USE mo_io_units,  ONLY: filename_max
  USE mo_grid
  USE mo_grid_levels
  USE mo_io_grid
  USE mo_gridref
  USE mo_util_uuid,  ONLY: t_uuid, uuid_generate, uuid_unparse, uuid_string_length
  USE mo_impl_constants, ONLY: max_dom

  IMPLICIT NONE
  
  PRIVATE
  
  !--------------------------------------------------------------
  !  Public subroutines
  PUBLIC global_graph_generator, global_grid_generator, prepare_gridref
CONTAINS
  
  SUBROUTINE global_graph_generator()
    
    INTEGER :: grid_levels  ! number of grid levels
    INTEGER :: nroot        ! number of initial sub division of basic icosahedron
    LOGICAL :: lplane       ! planar and not spherical grid construction
    
    !--------------------------------------------------------------------
    !BOC
    
    ! Switch logging on
    
    CALL open_log ('graph_generator.log')      ! setup log file
    
    CALL message ('', '')
    CALL message ('', 'Initialize graph generation ...')
    
    ! Read namelist to initialize graph generator    
    CALL init_graphgen(grid_levels, nroot, lplane)
    
    CALL message ('', 'Run graph generation ...')
    
    !  Build graph on basic icosahedron    
    CALL init_topology_graph (grid_levels, nroot, lplane)
    
    !  Generate tree structure    
    CALL message ('', 'Generate tree structure ...')
    
    CALL generate_tree (lplane)
    
    !  Fill graph tree    
    CALL message ('', 'Fill tree structure ...')
    
    CALL generate_graph(grid_levels, lplane)
    
    !  Store graph in GRAPHMAP file    
    CALL message ('', 'Write graphs ...')
    
    CALL write_graph(nroot, grid_levels, lplane)
    
    CALL message ('',  'Finished graph generation.')
    
    CALL destruct_topology_graph (grid_levels)
    
    CALL message ('', '')
    
  END SUBROUTINE global_graph_generator
  
  
  SUBROUTINE global_grid_generator()
    
    INTEGER :: grid_levels  ! number of grid levels
    INTEGER :: nroot        ! number of initial sub division of basic icosahedron
    LOGICAL :: lplane       ! planar and not spherical grid construction
    LOGICAL :: lread_graph  ! switch to decide if graph is read from file or computed from scratch
                            ! (the latter replaces running the graph generator before the grid generator)

    !--------------------------------------------------------------------
    ! Switch logging on    
    CALL open_log ('grid_generator.log')      ! setup log file
    
    CALL message ('', '')
    CALL message ('', 'Initialize grid generation ...')
    
    ! Read values to initialize grid generator
    
    CALL init_gridgen(grid_levels, nroot, lplane, lread_graph)
    
    CALL message ('', 'Run grid generation ...')
    
    !  Allocate everything needed
    
    CALL init_topology_grid (grid_levels, nroot, lplane)
    
    !  Generate tree structure
    
    CALL message ('', 'Generate tree structure ...')
    
    CALL generate_tree (lplane)

    IF (lread_graph) THEN
    
      ! Read graph from precomputed file
    
      CALL message ('', 'Read graphs ...')
    
      CALL read_graph(grid_levels, lplane)

    ELSE

      !  Fill graph tree    
      CALL message ('', 'Fill tree structure ...')
    
      CALL generate_graph(grid_levels, lplane)

    ENDIF
    
    !  NOTE: the grid optimization calls are all done in
    !   generate_geometry
    
    CALL generate_geometry(lplane)
    
    CALL init_grid(lplane)
    
    CALL message ('', 'Write grids ...')
    
    CALL output_grid(lplane)
    
    CALL statistics(lplane)
    
    CALL message ('', 'Finished grid generation.')
    
    ! cleanup section
    
    CALL destruct_grids
    
    CALL destruct_topology_grid (grid_levels)
    
    CALL message ('', '')
    
  END SUBROUTINE global_grid_generator
  
  
  SUBROUTINE prepare_gridref()
    
    CHARACTER(LEN=filename_max) :: filename
    
    INTEGER :: i, ilev, idom, ip, iplev, ic, in, iclev, istartlev, n_dom_start
    
    TYPE(t_patch), POINTER :: p_patch=>NULL()
    TYPE(t_patch), POINTER :: pc_patch=>NULL() ! for child level
    TYPE(t_patch), POINTER :: pp_patch=>NULL() ! for parent level
    TYPE(t_cell_list),POINTER :: local_cell_list=>NULL()
    TYPE(t_cell_list),POINTER :: parent_cell_list=>NULL()
    
    TYPE(t_patch), TARGET :: global_parent

    TYPE(t_uuid)                      :: uuid
    CHARACTER(len=uuid_string_length) :: uuid_grid(0:max_dom), uuid_parent(0:max_dom), &
                                         uuid_child(0:max_dom,5)

    LOGICAL :: child_exist, parent_exist
    
    CALL message ('', '')
    WRITE(message_text,'(a)')  'Preparation program for grid refinement'
    CALL message ('', TRIM(message_text))
    
    CALL init_gridref
    CALL check_namelist_gridref
    
    ! write_hierarchy = 2 enforces for writing output with mapped,
    ! parent-child index relationships back to grid level 0
    ! This is useful for search algorithms (external data etc.)
    IF (write_hierarchy == 2) THEN
      istartlev = 0
      n_dom_start = 0
    ELSE IF (write_hierarchy == 1) THEN
      istartlev = start_lev - 1
      n_dom_start = 0
    ELSE
      istartlev = start_lev
      n_dom_start = 1
    ENDIF
    !
    !  read grid quantities from files produced by grid generator
    !
    ALLOCATE (grid_on_level(istartlev:end_lev))
    
    DO ilev = istartlev, end_lev
      
      WRITE(filename,'(a,i0,a,i2.2,a)')'iconR', grid_root, 'B', ilev, '-grid.nc'

      CALL input_grid(grid_on_level(ilev), filename)

      ! Initialize refin_ctrl marker field
      grid_on_level(ilev)%cells%refin_ctrl(:) = 0
      grid_on_level(ilev)%edges%refin_ctrl(:) = 0
      grid_on_level(ilev)%verts%refin_ctrl(:) = 0
      
    ENDDO
    !
    ! allocate patch hierarchy
    !
    ALLOCATE (hierarchy(n_dom))
    
    ! allocate cell index list fields
    ALLOCATE( ptr_cell_list(n_dom), ptr_inv_cell_list(n_dom) )
    ALLOCATE( ptr_parent_cell_list(n_dom), ptr_inv_parent_cell_list(n_dom) )
    ALLOCATE( ptr_phys_cell_list(n_phys_dom), ptr_phys_inv_cell_list(n_phys_dom) )
    ALLOCATE( ptr_phys_parent_cell_list(n_phys_dom), ptr_phys_inv_parent_cell_list(n_phys_dom) )
    
    WRITE(message_text,'(a)') 'Start prepare_gridref'
    CALL message ('', TRIM(message_text))
    
    IF (n_dom > 1) THEN
      CALL setup_index_lists(grid_on_level(start_lev:end_lev), 1)
    ENDIF
    
    ! Merge nested domains located at the same grid level to a
    ! smaller number of logical grids if requested in the namelist
    CALL merge_nested_domains(grid_on_level(start_lev:end_lev))
    
    p_patch=>hierarchy(1)
    
    CALL create_global_domain(p_patch,grid_on_level(start_lev))
    
    DO idom = 2, n_dom
      
      ilev = grid_level(idom)
      ip = parent_grid_id(idom)
      iplev = grid_level(ip)
      
      p_patch          => hierarchy(idom)
      local_cell_list  => ptr_cell_list(idom)
      parent_cell_list => ptr_parent_cell_list(idom)
      
      CALL create_local_domain(p_patch,grid_on_level(ilev), &
        & local_cell_list%ip, idom)
      
    ENDDO
    
    DO idom = n_dom, 1, -1
      
      p_patch  => hierarchy(idom)
      ilev = grid_level(idom)
      ip = parent_grid_id(idom)
      
      IF (ip > 0) THEN
        iplev = grid_level(ip)
        pp_patch => hierarchy(ip)
      ELSE
        iplev = ilev
        pp_patch => hierarchy(idom)
      ENDIF
            
      IF (idom > 1) THEN
        parent_exist = .TRUE.
      ELSE
        parent_exist = .FALSE.
      ENDIF

      IF (n_childdom(idom) > 5) THEN
        WRITE(message_text,'(a)')  'At most 5 child domains per parent grid are admitted'
        CALL message ('', TRIM(message_text))
      ENDIF

      IF (n_childdom(idom) > 0) THEN
        DO in = 1, n_childdom(idom)
          ic = child_id(idom,in)
          iclev = grid_level(ic)
          pc_patch => hierarchy(ic)
          child_exist = .TRUE.
          CALL complete_index_lists(p_patch,pp_patch,pc_patch,ic,     &
            & grid_on_level(ilev),grid_on_level(iplev),grid_on_level(iclev),&
            & child_exist,parent_exist)
        ENDDO
      ELSE
        pc_patch => hierarchy(idom)
        child_exist = .FALSE.
        iclev = ilev
        CALL complete_index_lists(p_patch,pp_patch,pc_patch,1, &
          & grid_on_level(ilev),grid_on_level(iplev),grid_on_level(iclev),&
          & child_exist,parent_exist)
      ENDIF
      
    ENDDO
    
    ! If requested, create global parent patch.
    ! This is required for computing physics parameterizations on a reduced grid
    IF(n_dom_start == 0) THEN
      p_patch  => hierarchy(1)
      CALL create_global_parent_domain(p_patch, global_parent, grid_on_level(start_lev-1))
    ENDIF
    
    ! Generate uuid's where needed
    DO idom = n_dom_start, n_dom
      i = idom+1-n_dom_start
      IF (TRIM(uuid_sourcefile(i)) /= 'EMPTY') THEN
        CALL read_uuid(uuid_sourcefile(i), uuid_grid(idom))
      ELSE
        ! UUID is generated as fingerprint of "clon" field:
        CALL uuid_generate(REAL(p_patch%cells%center(:)%lon, C_DOUBLE), SIZE(p_patch%cells%center), uuid)
        CALL uuid_unparse(uuid, uuid_grid(idom))
      ENDIF
      IF (idom==1 .AND. n_dom_start==0) THEN
        uuid_child(0,1)  = uuid_grid(1)
        uuid_child(0,2:) = ''
        uuid_parent(1)   = uuid_grid(0)
        uuid_parent(0)   = uuid_grid(0)
      ELSE IF (idom==1) THEN
        uuid_parent(1)   = uuid_grid(1)
      ENDIF
    ENDDO

    ! Map uuid's between parent and child grids
    DO idom = n_dom, 1, -1
      IF (n_childdom(idom) > 0) THEN
        DO in = 1, n_childdom(idom)
          ic = child_id(idom,in)
          uuid_child(idom,in) = uuid_grid(ic)
          uuid_parent(ic)  = uuid_grid(idom)
        ENDDO
          uuid_child(idom,n_childdom(idom)+1:) = ''
      ELSE
        uuid_child(idom,:) = ''
      ENDIF
    ENDDO

    DO idom = n_dom_start, n_dom
      
      IF (idom==0) THEN
        p_patch => global_parent
      ELSE
        p_patch => hierarchy(idom)
      ENDIF
      
      IF (l_plot .AND. idom>0) THEN
        ilev = grid_level(idom)
        CALL plot_local_domains(p_patch,grid_on_level(ilev))
      ENDIF
      
      CALL write_patch(p_patch, uuid_grid(idom), uuid_parent(idom), uuid_child(idom,:))
      
    ENDDO
    
    ! Write output back to grid level zero if write_hierarchy=2,
    ! by default, write output only for the parent grid of the global grid
    DO ilev = istartlev , start_lev - 2
      
      WRITE(filename,'(a,i0,a,i2.2,a)') 'iconR', grid_root, 'B', ilev, '.nc'
      
      CALL write_grid(grid_on_level(ilev),filename)
      
    ENDDO
            
    WRITE(message_text,'(a)')  'Done prepare_gridref'
    CALL message ('', TRIM(message_text))
    
    ! cleanup section
    DO i= istartlev,end_lev
      CALL destruct_grid (grid_on_level(i))
    ENDDO
    
    DEALLOCATE (grid_on_level)
    DEALLOCATE (hierarchy)
    DEALLOCATE (ptr_cell_list,ptr_inv_cell_list,ptr_parent_cell_list,&
      & ptr_inv_parent_cell_list)
    
    CALL message ('', '')
    
  END SUBROUTINE prepare_gridref
  
END MODULE mo_global_grid_generator

