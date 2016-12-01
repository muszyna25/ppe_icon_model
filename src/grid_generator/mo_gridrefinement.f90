!>
!!  This module replaces the former mo_hierarchy and allows for generating.
!!
!!  This module replaces the former mo_hierarchy and allows for generating
!!  refined model domains with (relatively) arbitrary domain configurations.
!!  In contrast to the previous version, multiple nests on the same grid level
!!  are possible.
!!
!! @par Revision History
!! Initial version of mo_hierarchy by:
!! @par
!! Luca Bonaventura,  MPI-M, Hamburg, January 2005
!! @par
!! Implementing boundary exchange for parallel patches
!! Almut Gassmann,  MPI-M, Hamburg, February 2007
!! @par
!! Almut Gassmann,  MPI-M, Hamburg, April 2007
!! - Providing grid information on patches so that the global grid becomes
!!   obsolete
!! - removed edges%neighbor_index (not used)
!! @par
!! Guenther Zaengl, DWD, 2008-10-10
!! Implement preparations for mesh refinement, mainly
!! - compute parent/child index relations for edges
!! - map parent/child index relations for cells from global grid to patch
!! - compute refin_ctrl variable for flow control
!! - reorder index lists according to refin_ctrl value for efficiency
!!   (developed between March and September 2008)
!! Guenther Zaengl, DWD, 2008-10-15
!! Change subroutine write_output from binary to netcdf. Patch and halo
!! cells are not included because the parallelization concept needs to be
!! reconsidered anyway. The old output routine is kept under the name
!! 'write_output_old' for the time being.
!! @par
!! Initial version of mo_gridref by
!!   Guenther Zaengl, DWD, 2009-07-22
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_gridref
  !-------------------------------------------------------------------------
  !
  !    ProTeX FORTRAN source: Style 2
  !    modified for ICON project, DWD/MPI-M 2005
  !
  !-------------------------------------------------------------------------
  !
  !
  !
  USE mo_kind,               ONLY: wp
  USE mo_math_constants,     ONLY: rad2deg,pi,pi2
  USE mo_physical_constants, ONLY: earth_radius
  USE mo_exception,          ONLY: message_text, message, finish
  USE mo_io_units,           ONLY: ngmt, nnml, filename_max
  USE mo_namelist,           ONLY: position_nml, open_nml, positioned
  USE mo_grid,               ONLY: t_grid,                                  &
       &                           t_grid_cells, t_grid_vertices, t_grid_edges
  USE mo_math_utilities,     ONLY: gvec2cvec, t_cartesian_coordinates,      &
       &                           gc2cc, arc_length
  USE mo_grid_levels,        ONLY: nf
  USE mo_math_utilities,     ONLY: rotate_latlon, check_orientation
  USE mo_impl_constants,     ONLY: min_rlcell, max_rlcell, min_rlvert,    &
       &                           max_rlvert, min_rledge, max_rledge,    &
       &                           min_rlcell_int, min_rlvert_int,        &
       &                           min_rledge_int, max_dom, max_phys_dom
  USE mo_util_uuid_types,    ONLY: uuid_string_length

  IMPLICIT NONE

  PRIVATE

  INCLUDE 'netcdf.inc'

  ! !modules interface objects
  PUBLIC t_patch
  PUBLIC hierarchy
  PUBLIC grid_root
  PUBLIC start_lev
  PUBLIC n_dom, n_phys_dom
  PUBLIC init_gridref
  PUBLIC check_namelist_gridref
  PUBLIC deallocate_patch
  PUBLIC create_global_domain
  PUBLIC create_local_domain
  PUBLIC create_global_parent_domain
  PUBLIC complete_index_lists
  PUBLIC map_indices
  PUBLIC setup_index_lists
  PUBLIC merge_nested_domains
  PUBLIC plot_local_domains
  PUBLIC write_patch, write_reduced_output
  PUBLIC t_cell_list
  PUBLIC ptr_cell_list, ptr_parent_cell_list, ptr_inv_cell_list, &
    & ptr_inv_parent_cell_list
  PUBLIC ptr_phys_cell_list, ptr_phys_parent_cell_list, ptr_phys_inv_cell_list, &
    & ptr_phys_inv_parent_cell_list
  PUBLIC grid_level, n_childdom, parent_grid_id, child_id, max_childdom, &
    & end_lev, l_plot, write_hierarchy, uuid_sourcefile

  ! Namelist variables
  INTEGER :: grid_root, start_lev, n_dom, parent_id(max_phys_dom-1), &
    & bdy_indexing_depth, logical_id(max_phys_dom-1), n_phys_dom, write_hierarchy
    
  LOGICAL :: l_circ, l_rotate, l_plot, lsep_gridref_info
  CHARACTER(LEN=filename_max) :: uuid_sourcefile(max_dom+1)
  REAL(wp), DIMENSION(max_phys_dom-1) :: radius, center_lon, center_lat, &
    & hwidth_lon, hwidth_lat

  ! Additional global variables needed to simplify flow control
  INTEGER, DIMENSION(max_phys_dom) :: phys_grid_level,n_phys_childdom,logical_grid_id
  INTEGER, DIMENSION(max_dom) :: grid_level,n_childdom,parent_grid_id, &
    & n_physdom

  INTEGER :: phys_child_id(max_phys_dom,max_phys_dom)
  INTEGER :: child_id(max_dom,max_dom),physdom_list(max_dom,max_phys_dom)
  INTEGER :: max_childdom, end_lev

  INTEGER, DIMENSION (max_dom)      :: n_parent_count, n_patch_count
  INTEGER, DIMENSION (max_phys_dom) :: n_phys_parent_count, n_phys_patch_count

  !---------------------------------------
  !
  ! type to handle cell lists
  !
  TYPE t_cell_list
    INTEGER,POINTER :: ip(:)
  END TYPE t_cell_list

  ! Logical and physical cell lists
  TYPE(t_cell_list), POINTER :: ptr_cell_list(:)                 => NULL()
  TYPE(t_cell_list), POINTER :: ptr_parent_cell_list(:)          => NULL()
  TYPE(t_cell_list), POINTER :: ptr_inv_cell_list(:)             => NULL()
  TYPE(t_cell_list), POINTER :: ptr_inv_parent_cell_list(:)      => NULL()
  TYPE(t_cell_list), POINTER :: ptr_phys_cell_list(:)            => NULL()
  TYPE(t_cell_list), POINTER :: ptr_phys_parent_cell_list(:)     => NULL()
  TYPE(t_cell_list), POINTER :: ptr_phys_inv_cell_list(:)        => NULL()
  TYPE(t_cell_list), POINTER :: ptr_phys_inv_parent_cell_list(:) => NULL()

  !----------------------------------------

  TYPE t_patch
    !
    ! ! pointer to global grid on which patch lives
    !
    TYPE(t_grid),   POINTER :: pgrid
    !
    ! ! level in grid hierarchy on which patch lives
    !
    INTEGER :: level
    !
    ! domain ID of current domain
    INTEGER :: id
    !
    ! domain ID of parent domain
    INTEGER :: parent_id
    !
    ! list of number of physical domains contained in each logical domain
    INTEGER :: n_phys_dom(max_dom)
    !
    ! list of child domain ID's
    INTEGER :: child_id(max_dom)
    !
    ! actual number of child domains
    INTEGER :: n_childdom
    !
    ! maximum number of child domains
    INTEGER :: max_childdom
    !
    ! number of parent cells overlapping with the present nest
    INTEGER :: n_parent_cells
    !
    !
    ! ! number of patch items (total no. of cells, edges, vertices)
    !
    INTEGER :: ncells
    INTEGER :: nedges
    INTEGER :: nverts
    !
    ! ! index arrays for cells /edges/vertices inside a patch
    !
    INTEGER, POINTER :: ic(:) => NULL()
    INTEGER, POINTER :: ie(:) => NULL()
    INTEGER, POINTER :: iv(:) => NULL()
    !
    ! ! index arrays for global cells and edges
    !
    INTEGER, POINTER :: iic(:) => NULL()
    INTEGER, POINTER :: iie(:) => NULL()
    INTEGER, POINTER :: iiv(:) => NULL()
    !
    !
    ! ! pointers to patches descending from present patch
    !
    TYPE(t_patch),  POINTER :: children(:)
    !
    ! ! grid information on the patch
    !
    TYPE(t_grid_cells)    :: cells
    TYPE(t_grid_edges)    :: edges
    TYPE(t_grid_vertices) :: verts

  END TYPE t_patch
  !
  !
  ! ! patch object, global to the hierarchy module
  !
  TYPE(t_patch), ALLOCATABLE, TARGET :: hierarchy(:)


  INTEGER, PARAMETER, PRIVATE ::max_length = 1024 ! max length of a char array
  !

  !-------------------------------------------------------------------------
  ! meta data descriptors
  CHARACTER(len=*), PARAMETER :: uri_pathname = 'http://icon-downloads.mpimet.mpg.de/grids/public/'
  INTEGER :: number_of_grid_used(max_dom+1)! index for reference to xml table for grib2, 0 for private use
  INTEGER :: centre                        ! centre running the grid generator: 78 - edzw (DWD), 252 - MPIM
  INTEGER :: subcentre                     ! subcentre to be assigned by centre, usually 0
  INTEGER :: outname_style                 ! Naming convention
  !
!  INTEGER :: annotate_level          ! which grid level to annotate
  !
  !-------------------------------------------------------------------------

CONTAINS

  !-------------------------------------------------------------------------
  !
  !
  !>
  !!                Reads the namelist for grid refinement setup.
  !!
  !!
  SUBROUTINE init_gridref
    INTEGER :: i_status, grid_root_d, start_lev_d, n_dom_d, parent_id_d(max_phys_dom-1), &
      & bdy_indexing_depth_d, n_phys_dom_d, write_hierarchy_d, logical_id_d(max_phys_dom-1)
    REAL(wp), DIMENSION(max_phys_dom-1) :: radius_d, center_lon_d, center_lat_d, &
      & hwidth_lon_d, hwidth_lat_d
    LOGICAL :: l_circ_d, l_rotate_d, l_plot_d, lsep_gridref_info_d
    INTEGER :: i, j
    !--------------------------------------------------------------------
    !BOC
    !
    ! initialize with grid parameters from files
    !
    NAMELIST /gridref_ini/ grid_root, start_lev, n_dom, parent_id, l_plot,       &
      & l_circ, l_rotate, radius, center_lon, center_lat, n_phys_dom,            &
      & hwidth_lon, hwidth_lat, write_hierarchy, bdy_indexing_depth, logical_id, &
      & lsep_gridref_info, uuid_sourcefile
    NAMELIST /gridref_metadata/ number_of_grid_used, centre, subcentre, outname_style ! , &
!      & annotate_level

    ! set default values for gridref_ini
    grid_root  = 2
    start_lev  = 4
    n_dom      = 2
    n_phys_dom = n_dom

    ! Note: the namelist arrays parent_id, radius etc. are not defined for the
    ! global domain, so the parent_id(1), i.e. the parent ID of the first nest,
    ! is always 1
    DO i = 1, max_phys_dom-1
      parent_id(i) = i
      logical_id(i) = i + 1
    ENDDO

    l_circ     = .false.
    l_rotate   = .false.
    l_plot     = .false.
    uuid_sourcefile(:) = 'EMPTY'

    lsep_gridref_info = .false.

    bdy_indexing_depth = 12
    write_hierarchy    = 1

    radius     = 30.0_wp
    center_lat = 90.0_wp
    center_lon = 30.0_wp
    hwidth_lat = 20.0_wp
    hwidth_lon = 20.0_wp

    ! set default values for grid_options
    number_of_grid_used(:) = 0
    centre = 65535                ! missing value from WMO common code table C-11
    subcentre = 0
    outname_style = 1             ! Default naming convention
!    annotate_level = start_lev    ! level to annotate, to allow definition of higher level nests

    ! copy default values to "default" variables
    grid_root_d  = grid_root
    start_lev_d  = start_lev
    n_dom_d      = n_dom
    n_phys_dom_d = n_phys_dom
    parent_id_d  = parent_id
    logical_id_d = logical_id
    l_plot_d     = l_plot
    l_circ_d     = l_circ
    l_rotate_d   = l_rotate
    radius_d     = radius
    center_lat_d = center_lat
    center_lon_d = center_lon
    hwidth_lat_d = hwidth_lat
    hwidth_lon_d = hwidth_lon
    bdy_indexing_depth_d = bdy_indexing_depth
    write_hierarchy_d    = write_hierarchy
    lsep_gridref_info_d  = lsep_gridref_info

    ! Set logical_id and n_phys_dom temporarily back to -1 in order to check
    ! if they are specified in the namelist
    logical_id = -1
    n_phys_dom = -1

    ! read namelist
    CALL open_nml('NAMELIST_GRIDREF')

    CALL position_nml('gridref_ini',STATUS=i_status)
    IF (i_status == positioned) THEN
      READ (nnml,gridref_ini)
    ENDIF

    CALL position_nml('gridref_metadata', STATUS=i_status)
    IF (i_status == positioned) THEN
      READ (nnml,gridref_metadata)
    ENDIF

    CLOSE(nnml)

    ! Restore logical_id back to the physical domain ID if not specified
    IF (MAXVAL(logical_id) == -1) logical_id = logical_id_d

    ! Set n_phys_dom = n_dom if not specified
    IF (n_phys_dom == -1) n_phys_dom = n_dom

    ! Avoid unreasonably small values of bdy_indexing_depth
    bdy_indexing_depth = MAX(10,max_rlcell+1,bdy_indexing_depth)

    ! Define additional variables needed for flow control

    ! NOTE: the namelist variables are defined on physical model domains;
    ! the flow control variables refer to logical domains unless indicated by "phys"

    ! a) physical domains

    n_phys_childdom(:) = 0
    phys_child_id(:,:) = 0
    phys_grid_level(1) = start_lev
    logical_grid_id(1) = 1

    ! Determine the grid level of each domain; the "i-1" on the rhs is because
    ! parent_id is not defined for the global domain
    DO i = 2, n_phys_dom
      phys_grid_level(i) = phys_grid_level(parent_id(i-1)) + 1
      logical_grid_id(i) = logical_id(i-1)
    ENDDO

    DO i = 1, n_phys_dom-1
      DO j = i+1,n_phys_dom
        IF (i == parent_id(j-1)) THEN
          n_phys_childdom(i) = n_phys_childdom(i)+1
          phys_child_id(i,n_phys_childdom(i)) = j
        ENDIF
      ENDDO
    ENDDO

    end_lev = MAXVAL(phys_grid_level(1:n_phys_dom))

    ! b) logical domains

    n_childdom(:)  = 0
    child_id(:,:) = 0
    grid_level(1)  = start_lev
    max_childdom   = 0
    n_physdom(:) = 0
    physdom_list(:,:) = 0

    ! list holding the phyiscal domain IDs for a given logical domain
    physdom_list(1,1) = 1 ! for the global domain, logical and physical are always the same
    n_physdom(1)      = 1 ! which implies that there is exactly 1 physical domain
    DO i = 2, n_dom
      DO j = 2, n_phys_dom
        IF (logical_grid_id(j) == i) THEN
          n_physdom(i) = n_physdom(i) + 1
          physdom_list(i,n_physdom(i)) = j
        ENDIF
      ENDDO
    ENDDO

    ! parent_grid_id converts the information of parent_id to the logical space,
    ! and its first element refers to the global domain
    parent_grid_id(1) = 0
    DO i = 2, n_dom
      parent_grid_id(i) = logical_grid_id(parent_id(physdom_list(i,1)-1))
      ! for safety, we check if all physical domains merged into a logical domain
      ! have the same parent
      DO j = 2, n_physdom(i)
        IF (parent_grid_id(i) /= logical_grid_id(parent_id(physdom_list(i,j)-1))) THEN
          CALL finish  ('init_gridref', &
            &           'merged domains must have the same parent')
        ENDIF
      ENDDO
    ENDDO

    DO i = 2, n_dom
      grid_level(i) = grid_level(parent_grid_id(i)) + 1
    ENDDO

    DO i = 1, n_dom-1
      DO j = i+1,n_dom
        IF (i == parent_grid_id(j)) THEN
          n_childdom(i) = n_childdom(i)+1
          child_id(i,n_childdom(i)) = j
          max_childdom = MAX(1,max_childdom,n_childdom(i))
        ENDIF
      ENDDO
    ENDDO

    ! This is needed for the special case that patch0 is to be created without a nest
    max_childdom = MAX(1,max_childdom)


    ! write control output of namelist variables

    CALL message ('', '')
    WRITE(message_text,'(a)')'Namelist Group: gridref_ini'
    CALL message ('', message_text)
    WRITE(message_text,'(a)')'-----------------------------'
    CALL message ('', message_text)
    WRITE(message_text,'(t7,a,t28,a,t43,a,t59,a )') 'Variable', 'Actual Value','Default Value'
    CALL message ('', message_text)
    WRITE(message_text,'(t8,a,t28,i12,t44,i12,t60,a3)') 'grid_root', grid_root, grid_root_d
    CALL message ('', message_text)
    WRITE(message_text,'(t8,a,t28,i12,t44,i12,t60,a3)') 'start_lev', start_lev, start_lev_d
    CALL message ('', message_text)
    WRITE(message_text,'(t8,a,t28,i12,t44,i12,t60,a3)') 'n_dom', n_dom, n_dom_d
    CALL message ('', message_text)
    WRITE(message_text,'(t8,a,t28,i12,t44,i12,t60,a3)') 'n_phys_dom', n_phys_dom, n_phys_dom_d
    CALL message ('', message_text)
    WRITE(message_text,'(t8,a,t28,l12,t44,l12,t60,a3)') 'l_circ',l_circ, l_circ_d
    CALL message ('', message_text)
    WRITE(message_text,'(t8,a,t28,l12,t44,l12,t60,a3)') 'l_rotate',l_rotate, l_rotate_d
    CALL message ('', message_text)

    WRITE(message_text,'(t8,a,t28,i12,t44,i12,t60,a3)') 'bdy_indexing_depth', &
                                                         bdy_indexing_depth, bdy_indexing_depth_d
    CALL message ('', message_text)
    WRITE(message_text,'(t8,a,t28,i12,t44,i12,t60,a3)') 'write_hierarchy', &
                                                         write_hierarchy, write_hierarchy_d
    CALL message ('', message_text)
    WRITE(message_text,'(t8,a,t28,l12,t44,l12,t60,a3)') 'lsep_gridref_info',lsep_gridref_info, lsep_gridref_info_d
    CALL message ('', message_text)

    DO i = 2, n_phys_dom
      j = i - 1
      WRITE(message_text,'(a)')''
      CALL message ('', message_text)
      WRITE(message_text,'(a,i3)')'Domain-specific settings for nest', j
      CALL message ('', message_text)
      WRITE(message_text,'(a)')'--------------------------------'
      CALL message ('', message_text)
      WRITE(message_text,'(t7,a,t28,a,t43,a,t59,a )') &
        & 'Variable', 'Actual Value', 'Default Value'
      CALL message ('', message_text)
      WRITE(message_text,'(t8,a,t28,i12,t44,i12,t60,a3)') &
        & 'parent_id',parent_id(j), parent_id_d(j)
      CALL message ('', message_text)
      WRITE(message_text,'(t8,a,t28,i12,t44,i12,t60,a3)') &
        & 'logical_id',logical_id(j), logical_id_d(j)
      CALL message ('', message_text)
      WRITE(message_text,'(t8,a,t28,f12.4,t44,f12.4,t60,a3)') &
        & 'radius', radius(j), radius_d(j)
      CALL message ('', message_text)
      WRITE(message_text,'(t8,a,t28,f12.4,t44,f12.4,t60,a3)') &
        & 'center_lat', center_lat(j), center_lat_d(j)
      CALL message ('', message_text)
      WRITE(message_text,'(t8,a,t28,f12.4,t44,f12.4,t60,a3)') &
        & 'center_lon', center_lon(j), center_lon_d(j)
      CALL message ('', message_text)
      WRITE(message_text,'(t8,a,t28,f12.4,t44,f12.4,t60,a3)') &
        & 'hwidth_lat', hwidth_lat(j), hwidth_lat_d(j)
      CALL message ('', message_text)
      WRITE(message_text,'(t8,a,t28,f12.4,t44,f12.4,t60,a3)') &
        & 'hwidth_lon', hwidth_lon(j), hwidth_lon_d(j)
      CALL message ('', message_text)
      CALL message ('', '')
    ENDDO

  END SUBROUTINE init_gridref
  !EOC
  !-------------------------------------------------------------------------
  !
  !
  !>
  !!                Checks the namelist settings for grid refinement.
  !!
  !!
  !! @par Revision History
  !!
  SUBROUTINE check_namelist_gridref
    !--------------------------------------------------------------------
    !BOC

    IF((grid_root /= 2) .and. (grid_root /= 3) .and. (grid_root /= 5) .and. (grid_root /= 7)) THEN
      CALL finish  ('check_namelist_gridref', 'grid_root must be either 2, 3, 5 or 7')
    ENDIF
    IF(start_lev<1) THEN
      CALL finish  ('check_namelist_gridref', 'startlev must be >=1')
    ENDIF
    IF(parent_id(1) /= 1) THEN
      CALL finish  ('check_namelist_gridref', 'parent_id(1) must be 1')
    ENDIF

    IF((n_dom > max_dom)) THEN
      CALL finish  ('check_namelist_gridref', &
                    'the maximum allowed number of model domains is max_dom')
    ENDIF

  END SUBROUTINE check_namelist_gridref


  !
  !EOC
  !-------------------------------------------------------------------------
  !
  !
  !>
  !!  Creates the grid for the global model domain.
  !!
  !! Based on the former
  !!    make_global_patch
  !!
  SUBROUTINE create_global_domain(p,g)
    !
    ! !REVISION HISTORY (make_global_patch):
    !  Luca Bonaventura, MPI-M, Hamburg, March 2005
    !  Almut Gassmann, MPI-M, Hamburg, April 2007
    !  Guenther Zaengl, DWD, July 2008:
    !    reorder index list if refined meshes are present
    !
    !  Guenther Zaengl, DWD, 2009-07-22:
    !    Initial version of create_global_domain
    !
    TYPE(t_grid) , TARGET, INTENT(in)  :: g
    TYPE(t_patch) :: p
    ! local variables:
    INTEGER, POINTER :: list(:), itmpe(:), itmpv(:)
    INTEGER :: irv, i, j, k, g_cell, g_edge, g_vert, l2g, i1, icd, nest_id
    !-------------------------------------------------------------------------
    !BOC

    p%pgrid=>g

    p%level=g%level

    p%ncells=g%ncells
    p%nedges=g%nedges
    p%nverts=g%nverts

    ! Allocate all fields in p

    CALL allocate_patch(p)

    p%id = 1
    p%parent_id = 0
    p%child_id(:) = child_id(1,:)
    p%n_childdom  = n_childdom(1)
    p%n_phys_dom  = n_physdom(1)
    p%max_childdom = max_childdom

    ALLOCATE(list(p%ncells))  ! Reordered index list for cells

    p%cells%indlist(1:max_rlcell,1) = 1
    p%cells%indlist(min_rlcell:0,1) = p%ncells+1
    p%cells%indlist(1:max_rlcell,2) = 0
    p%cells%indlist(min_rlcell:0,2) = p%ncells

    p%cells%start_idx(1:max_rlcell,1) = 1
    p%cells%start_idx(min_rlcell:0,1) = p%ncells+1
    p%cells%end_idx(1:max_rlcell,1)   = 0
    p%cells%end_idx(min_rlcell:0,1)   = p%ncells

    ! The following loop used to accomplish the reordering of nest overlap points
    ! (deactivated since Aug. 2013 for restructured nesting flow control)
    k = 0
    irv = 0
    DO i = 1, p%ncells
      IF (g%cells%refin_ctrl(i) <= irv) THEN
        k = k + 1
        list(k) = i
        p%cells%indlist(irv,1) = MIN(p%cells%indlist(irv,1),k)
        p%cells%indlist(irv,2) = k
        p%cells%start_idx(irv,1) = MIN(p%cells%start_idx(irv,1),k)
        p%cells%end_idx(irv,1) = k
      ENDIF
    ENDDO

    IF (k /= p%ncells) CALL finish ('create_global_domain', &
      & 'Error in reordering global cell index list')

    !  allocate temporary  global index array for edges
    !
    ALLOCATE( itmpe(p%nedges) )

    p%edges%indlist(1:max_rledge,1) = 1
    p%edges%indlist(min_rledge:0,1) = p%nedges+1
    p%edges%indlist(1:max_rledge,2) = 0
    p%edges%indlist(min_rledge:0,2) = p%nedges

    p%edges%start_idx(1:max_rledge,1) = 1
    p%edges%start_idx(min_rledge:0,1) = p%nedges+1
    p%edges%end_idx(1:max_rledge,1) = 0
    p%edges%end_idx(min_rledge:0,1) = p%nedges

    ! The following loop used to accomplish the reordering of nest overlap points
    ! (deactivated since Aug. 2013 for restructured nesting flow control)
    k = 0
    irv = 0
    DO i = 1, p%nedges
      IF (g%edges%refin_ctrl(i) <= irv) THEN
        k = k + 1
        itmpe(k) = i
        p%edges%indlist(irv,1) = MIN(p%edges%indlist(irv,1),k)
        p%edges%indlist(irv,2) = k
        p%edges%start_idx(irv,1) = MIN(p%edges%start_idx(irv,1),k)
        p%edges%end_idx(irv,1) = k
      ENDIF
    ENDDO

    IF (k /= p%nedges) CALL finish ('create_global_domain', &
      & 'Error in reordering global edge index list')

    !  allocate global index array for vertices
    !
    ALLOCATE( itmpv(p%nverts) )

    p%verts%indlist(1:max_rlvert,1) = 1
    p%verts%indlist(min_rlvert:0,1) = p%nverts+1
    p%verts%indlist(1:max_rlvert,2) = 0
    p%verts%indlist(min_rlvert:0,2) = p%nverts

    p%verts%start_idx(1:max_rlvert,1) = 1
    p%verts%start_idx(min_rlvert:0,1) = p%nverts+1
    p%verts%end_idx(1:max_rlvert,1) = 0
    p%verts%end_idx(min_rlvert:0,1) = p%nverts

    ! The following loop used to accomplish the reordering of nest overlap points
    ! (deactivated since Aug. 2013 for restructured nesting flow control)
    k = 0
    irv = 0
    DO i = 1, p%nverts
      IF (g%verts%refin_ctrl(i) <= irv) THEN
        k = k + 1
        itmpv(k) = i
        p%verts%indlist(irv,1) = MIN(p%verts%indlist(irv,1),k)
        p%verts%indlist(irv,2) = k
        p%verts%start_idx(irv,1) = MIN(p%verts%start_idx(irv,1),k)
        p%verts%end_idx(irv,1) = k
      ENDIF
    ENDDO

    IF (k /= p%nverts) CALL finish ('create_global_domain', &
      & 'Error in reordering global vertex index list')

    !
    !  allocate and fill arrays holding original global indices for edges,
    !  cells, vertices
    !
    ALLOCATE(p%ic(1:p%ncells),p%iic(1:p%ncells))
    ALLOCATE(p%ie(1:p%nedges),p%iie(1:p%nedges))
    ALLOCATE(p%iv(1:p%nverts),p%iiv(1:p%nverts))

    p%ic(1:p%ncells)=list (1:p%ncells)
    p%ie(1:p%nedges)=itmpe(1:p%nedges)
    p%iv(1:p%nverts)=itmpv(1:p%nverts)

    DO k=1,p%ncells
      p%iic(p%ic(k)) = k
    ENDDO

    DO k=1,p%nedges
      p%iie(p%ie(k)) = k
    ENDDO

    DO k=1,p%nverts
      p%iiv(p%iv(k)) = k
    ENDDO

    !   deallocate auxiliary arrays
    !
    DEALLOCATE(list, itmpv, itmpe)

    ! Map grid relationships to reordered index lists (copied from make_local_patch)
    DO j = 1, p%ncells

      p%cells%idx(j) = j
      g_cell = p%ic(j)
      p%cells%parent_index(j) = 0
      ! Child indices are set later in subroutine complete_local_patch
      p%cells%child_index(j,1:4) = 0

      l2g = g%cells%neighbor_index(g_cell,1)
      p%cells%neighbor_index(j     ,1) = p%iic(l2g)
      l2g = g%cells%neighbor_index(g_cell,2)
      p%cells%neighbor_index(j     ,2) = p%iic(l2g)
      l2g = g%cells%neighbor_index(g_cell,3)
      p%cells%neighbor_index(j     ,3) = p%iic(l2g)

      l2g = g%cells%edge_index(g_cell,1)
      p%cells%edge_index(j     ,1) = p%iie(l2g)
      l2g = g%cells%edge_index(g_cell,2)
      p%cells%edge_index(j     ,2) = p%iie(l2g)
      l2g = g%cells%edge_index(g_cell,3)
      p%cells%edge_index(j     ,3) = p%iie(l2g)

      l2g = g%cells%vertex_index(g_cell,1)
      p%cells%vertex_index(j     ,1) = p%iiv(l2g)
      l2g = g%cells%vertex_index(g_cell,2)
      p%cells%vertex_index(j     ,2) = p%iiv(l2g)
      l2g = g%cells%vertex_index(g_cell,3)
      p%cells%vertex_index(j     ,3) = p%iiv(l2g)

      p%cells%edge_orientation(j,1) = g%cells%edge_orientation(g_cell,1)
      p%cells%edge_orientation(j,2) = g%cells%edge_orientation(g_cell,2)
      p%cells%edge_orientation(j,3) = g%cells%edge_orientation(g_cell,3)
      p%cells%center          (j  ) = g%cells%center          (g_cell  )
      p%cells%area            (j  ) = g%cells%area            (g_cell  )
      p%cells%refin_ctrl      (j  ) = g%cells%refin_ctrl      (g_cell  )
      p%cells%child_id        (j  ) = g%cells%child_id        (g_cell  )
      p%cells%phys_id         (j  ) = g%cells%phys_id         (g_cell  )

    ENDDO

    DO j = 1, p%nedges

      p%edges%idx(j) = j

      g_edge = p%ie(j)

      l2g = g%edges%cell_index(g_edge,1)
      p%edges%cell_index(j     ,1) = p%iic(l2g)
      l2g = g%edges%cell_index(g_edge,2)
      p%edges%cell_index(j     ,2) = p%iic(l2g)

      l2g = g%edges%vertex_index(g_edge,1)
      p%edges%vertex_index(j     ,1) = p%iiv(l2g)
      l2g = g%edges%vertex_index(g_edge,2)
      p%edges%vertex_index(j     ,2) = p%iiv(l2g)

      p%edges%system_orientation(j) = g%edges%system_orientation(g_edge)
      p%edges%center            (j) = g%edges%center            (g_edge)
      p%edges%primal_normal     (j) = g%edges%primal_normal     (g_edge)
      p%edges%dual_normal       (j) = g%edges%dual_normal       (g_edge)
      p%edges%primal_edge_length(j) = g%edges%primal_edge_length(g_edge)
      p%edges%dual_edge_length  (j) = g%edges%dual_edge_length  (g_edge)
      p%edges%edge_vert_length(j,:) = g%edges%edge_vert_length(g_edge,:)
      p%edges%edge_cell_length(j,:) = g%edges%edge_cell_length(g_edge,:)
      p%edges%refin_ctrl        (j) = g%edges%refin_ctrl        (g_edge)
      p%edges%child_id          (j) = g%edges%child_id          (g_edge)
      p%edges%phys_id           (j) = g%edges%phys_id           (g_edge)

      p%edges%parent_index(j)       = 0
      ! Child indices are set later in subroutine complete_local_patch
      p%edges%child_index(j,1:4)    = 0

    ENDDO

    DO j = 1, p%nverts

      p%verts%idx(j) = j

      g_vert = p%iv(j)
      !
      ! note: the dummy edges, neighbor_vertices and cells get the index 0
      !
      DO i1 = 1, 6

        l2g = g%verts%neighbor_index(g_vert,i1)
        IF (l2g > 0) THEN
          p%verts%neighbor_index(j     ,i1) = p%iiv(l2g)
        ELSE
          p%verts%neighbor_index(j     ,i1) = 0
        ENDIF

        l2g = g%verts%cell_index(g_vert,i1)
        IF (l2g > 0) THEN
          p%verts%cell_index(j     ,i1) = p%iic(l2g)
        ELSE
          p%verts%cell_index(j     ,i1) = 0
        ENDIF

        l2g = g%verts%edge_index(g_vert,i1)
        IF (l2g > 0) THEN
          p%verts%edge_index(j     ,i1) = p%iie(l2g)
        ELSE
          p%verts%edge_index(j     ,i1) = 0
        ENDIF

      ENDDO

      p%verts%edge_orientation(j,1:6) = g%verts%edge_orientation(g_vert,1:6)
      p%verts%dual_area       (j    ) = g%verts%dual_area       (g_vert    )
      p%verts%vertex          (j    ) = g%verts%vertex          (g_vert    )
      p%verts%refin_ctrl      (j    ) = g%verts%refin_ctrl      (g_vert    )
      p%verts%child_id        (j    ) = g%verts%child_id        (g_vert    )
      p%verts%phys_id         (j    ) = g%verts%phys_id         (g_vert    )

    ENDDO

  END SUBROUTINE create_global_domain

  !EOC
  !-------------------------------------------------------------------------
  !
  !
  !>
  !!  Allocates the fields of the patch datatype.
  !!
  !!
  !! @par Revision History
  !!  Rainer Johanni, 2010-10-28
  !!
  SUBROUTINE allocate_patch(p)
    !

    TYPE(t_patch) :: p

    ! p%ncells, p%nedges, p%nverts must already be set
    ! p%ic, p%iic, p%ie, p%iie, p%iv, p%iiv are not allocated here!

    !
    ! allocate fields for grid attributes
    !
    ! cells
    !------
    ALLOCATE(p%cells%idx(p%ncells))
    ALLOCATE(p%cells%parent_index(p%ncells))
    ALLOCATE(p%cells%child_index(p%ncells,4))
    ALLOCATE(p%cells%child_id(p%ncells))
    ALLOCATE(p%cells%phys_id(p%ncells))
    ALLOCATE(p%cells%neighbor_index(p%ncells,3))
    ALLOCATE(p%cells%edge_index(p%ncells,3))
    ALLOCATE(p%cells%vertex_index(p%ncells,3))
    ALLOCATE(p%cells%edge_orientation(p%ncells,3))
    ALLOCATE(p%cells%center(p%ncells))
    ALLOCATE(p%cells%area(p%ncells))
    ALLOCATE(p%cells%refin_ctrl(p%ncells))

    ALLOCATE(p%cells%indlist(min_rlcell:max_rlcell,2))
    ALLOCATE(p%cells%start_idx(min_rlcell:max_rlcell,1))
    ALLOCATE(p%cells%end_idx(min_rlcell:max_rlcell,1))

    ! edges
    !------
    ALLOCATE(p%edges%idx(p%nedges))
    ALLOCATE(p%edges%cell_index(p%nedges,2))
    ALLOCATE(p%edges%vertex_index(p%nedges,2))
    ALLOCATE(p%edges%system_orientation(p%nedges))
    ALLOCATE(p%edges%center(p%nedges))
    ALLOCATE(p%edges%primal_normal(p%nedges))
    ALLOCATE(p%edges%dual_normal(p%nedges))
    ALLOCATE(p%edges%primal_edge_length(p%nedges))
    ALLOCATE(p%edges%dual_edge_length(p%nedges))
    ALLOCATE(p%edges%edge_vert_length(p%nedges,2))
    ALLOCATE(p%edges%edge_cell_length(p%nedges,2))
    ALLOCATE(p%edges%refin_ctrl(p%nedges))
    ALLOCATE(p%edges%parent_index(p%nedges))
    ALLOCATE(p%edges%child_index(p%nedges,4))
    ALLOCATE(p%edges%child_id(p%nedges))
    ALLOCATE(p%edges%phys_id(p%nedges))

    ALLOCATE(p%edges%indlist(min_rledge:max_rledge,2))
    ALLOCATE(p%edges%start_idx(min_rledge:max_rledge,1))
    ALLOCATE(p%edges%end_idx(min_rledge:max_rledge,1))

    ! vertices
    !---------
    ALLOCATE(p%verts%idx(p%nverts))
    ALLOCATE(p%verts%neighbor_index(p%nverts,6))
    ALLOCATE(p%verts%cell_index(p%nverts,6))
    ALLOCATE(p%verts%edge_index(p%nverts,6))
    ALLOCATE(p%verts%edge_orientation(p%nverts,6))
    ALLOCATE(p%verts%dual_area(p%nverts))
    ALLOCATE(p%verts%vertex(p%nverts))
    ALLOCATE(p%verts%refin_ctrl(p%nverts))
    ALLOCATE(p%verts%child_id(p%nverts))
    ALLOCATE(p%verts%phys_id(p%nverts))

    ALLOCATE(p%verts%indlist(min_rlvert:max_rlvert,2))
    ALLOCATE(p%verts%start_idx(min_rlvert:max_rlvert,1))
    ALLOCATE(p%verts%end_idx(min_rlvert:max_rlvert,1))

  END SUBROUTINE allocate_patch

  !-------------------------------------------------------------------------
  !
  !
  !>
  !!  Deallocates the fields of the patch datatype.
  !!
  !!
  !! @par Revision History
  !!  Th.Heinze, DWD, Offenbach, 2006-10-25
  !!  Almut Gassmann, MPI-M, Hamburg, 2007-02-08
  !!
  SUBROUTINE deallocate_patch(p)
    !

    TYPE(t_patch) :: p
    !-------------------------------------------------------------------------
    !BOC

    ! deallocate fields for grid attributes
    !
    ! cells
    !------
    DEALLOCATE(p%cells%idx)
    DEALLOCATE(p%cells%parent_index)
    DEALLOCATE(p%cells%child_index)
    DEALLOCATE(p%cells%child_id)
    DEALLOCATE(p%cells%phys_id)
    DEALLOCATE(p%cells%neighbor_index)
    DEALLOCATE(p%cells%edge_index)
    DEALLOCATE(p%cells%vertex_index)
    DEALLOCATE(p%cells%edge_orientation)
    DEALLOCATE(p%cells%center)
    DEALLOCATE(p%cells%area)
    DEALLOCATE(p%cells%refin_ctrl)
    DEALLOCATE(p%cells%indlist)
    DEALLOCATE(p%cells%start_idx)
    DEALLOCATE(p%cells%end_idx)

    ! edges
    !------
    DEALLOCATE(p%edges%idx)
    DEALLOCATE(p%edges%cell_index)
    DEALLOCATE(p%edges%vertex_index)
    DEALLOCATE(p%edges%system_orientation)
    DEALLOCATE(p%edges%center)
    DEALLOCATE(p%edges%primal_normal)
    DEALLOCATE(p%edges%dual_normal)
    DEALLOCATE(p%edges%primal_edge_length)
    DEALLOCATE(p%edges%dual_edge_length)
    DEALLOCATE(p%edges%edge_vert_length)
    DEALLOCATE(p%edges%edge_cell_length)
    DEALLOCATE(p%edges%refin_ctrl)
    DEALLOCATE(p%edges%indlist)
    DEALLOCATE(p%edges%parent_index)
    DEALLOCATE(p%edges%child_index)
    DEALLOCATE(p%edges%child_id)
    DEALLOCATE(p%edges%phys_id)
    DEALLOCATE(p%edges%start_idx)
    DEALLOCATE(p%edges%end_idx)

    ! vertices
    !---------
    DEALLOCATE(p%verts%idx)
    DEALLOCATE(p%verts%neighbor_index)
    DEALLOCATE(p%verts%cell_index)
    DEALLOCATE(p%verts%edge_index)
    DEALLOCATE(p%verts%edge_orientation)
    DEALLOCATE(p%verts%dual_area)
    DEALLOCATE(p%verts%vertex)
    DEALLOCATE(p%verts%refin_ctrl)
    DEALLOCATE(p%verts%indlist)
    DEALLOCATE(p%verts%start_idx)
    DEALLOCATE(p%verts%end_idx)
    DEALLOCATE(p%verts%child_id)
    DEALLOCATE(p%verts%phys_id)

    IF( ASSOCIATED(p%ic) ) DEALLOCATE(p%ic)
    IF( ASSOCIATED(p%ie) ) DEALLOCATE(p%ie)
    IF( ASSOCIATED(p%iv) ) DEALLOCATE(p%iv)

    IF( ASSOCIATED(p%iic) ) DEALLOCATE(p%iic)
    IF( ASSOCIATED(p%iie) ) DEALLOCATE(p%iie)
    IF( ASSOCIATED(p%iiv) ) DEALLOCATE(p%iiv)

    NULLIFY(p%children)
    NULLIFY(p%pgrid)


  END SUBROUTINE deallocate_patch

  !EOC
  !-------------------------------------------------------------------------
  !
  !
  !>
  !!  Replaces the former routine make_local_patch.
  !!
  !!  Maps (nearly) all elements of the patch datatype from the global
  !!  grid to the refined grid(s) defined in setup_index_lists.
  !!  The exception is the mapping of the parent/childd index relationship,
  !!  which is done in complete_index_lists (see below).
  !!
  SUBROUTINE create_local_domain(p, g, listc, idom)
    !
    ! !REVISION HISTORY (make_local_patch):
    !  Luca Bonaventura, MPI-M, Hamburg, March-August 2005
    !  Peter Korn, MPI-M, Hamburg, February 2006
    !  Almut Gassmann, MPI-M, Hamburg, February 2007
    !  - corrected and modified large parts
    !  Almut Gassmann, MPI-M, Hamburg, April 2007
    !  - added grid information on patch
    !  Guenther Zaengl, DWD, March-June 2008
    !  - complete fields needed for mesh refinement control
    !
    !  Guenther Zaengl, DWD, 2009-07-22:
    !    Initial version of create_local_domain
    !
    ! patch to be constructed
    !
    TYPE(t_patch) :: p
    !
    ! grid on which patch lives
    !
    TYPE(t_grid) , TARGET, INTENT(in)  :: g
    !
    !   list of cells produced by the make_lists_local_patch
    !   routines or by any other preprocessing software
    !   or even by hand.......
    !
    INTEGER, POINTER :: listc(:)
    ! Current domain ID
    INTEGER, INTENT(in) :: idom

    INTEGER :: ji, j, jj
    INTEGER :: jj1, jj2, jj3
    INTEGER :: i1
    INTEGER :: ie_count, iv_count, ic_count
    INTEGER, POINTER :: itmpe(:), itmpe2(:)
    INTEGER, POINTER :: itmpv(:), itmpv2(:)
    INTEGER, POINTER :: itmpc(:)

    INTEGER :: nc_int, ne_int, nv
    INTEGER :: ist, istat
    INTEGER :: g_cell, g_edge, g_vert, l2g, irv, icd, nest_id
    LOGICAL, ALLOCATABLE, DIMENSION(:) :: l_counted

    !-------------------------------------------------------------------------
    !BOC
    !
    !   set pointer to grid on which patch lives
    !
    p%pgrid => g

    p%level=g%level

    ! More patch attributes
    p%id = idom
    p%parent_id = parent_grid_id(idom)
    p%child_id(:) = child_id(idom,:)
    p%n_childdom  = n_childdom(idom)
    p%n_phys_dom  = n_physdom(idom)
    p%max_childdom = max_childdom

    !
    !   number of cells belonging to patch is known by size of
    !   array listc, which is provided by the make_list routines;
    !   it is a temporary solution and to be improved in the future
    !
    p%ncells = n_patch_count(idom)

    IF (idom > 1) THEN
      p%n_parent_cells = n_parent_count(idom)
    ELSE
      p%n_parent_cells = 0
    ENDIF
    !
    !  allocate local index arrays for edges/cells/vertices
    !  and initialize with default value;
    !
    ALLOCATE(p%iie(1:g%nedges))
    p%iie=-1

    ALLOCATE(p%iic(1:g%ncells))
    p%iic=-1

    ALLOCATE(p%iiv(1:g%nverts))
    p%iiv=-1

    ALLOCATE(l_counted(1:g%ncells))
    l_counted = .false.

    !
    ! initialize patch variables to preliminary values. All halos will be written in
    ! one vector.
    !

    p%nedges = 0
    ie_count = 0
    ic_count = 0

    !loop over triangles in patch list
    !(listc contains no of triangles that belong to the patch)
    DO ji=1,p%ncells

      !get indices of 3 neighboring triangles of triangle ji in patch list
      jj1 = g%cells%neighbor_index(listc(ji),1)
      jj2 = g%cells%neighbor_index(listc(ji),2)
      jj3 = g%cells%neighbor_index(listc(ji),3)

      !check if neighbor triangles are already counted
      IF (.not.l_counted(jj1)) p%nedges = p%nedges+1
      IF (.not.l_counted(jj2)) p%nedges = p%nedges+1
      IF (.not.l_counted(jj3)) p%nedges = p%nedges+1

      l_counted(listc(ji)) = .true.

    ENDDO
    !
    !  allocate temporary global index arrays for edges
    !
    ALLOCATE( itmpe(p%nedges),itmpe2(p%nedges),      &
      & p%edges%indlist(min_rledge:max_rledge,2),    &
      & p%edges%start_idx(min_rledge:max_rledge,1),  &
      & p%edges%end_idx(min_rledge:max_rledge,1)     )
    !
    !
    !  fill global index array for edges
    !
    ie_count = 0
    ic_count = 0
    l_counted = .false.

    !proceed like in the edge counting loop above:
    !get neighbor triangles of triangle in patch list
    !check if neighbors belong to patch list
    !if this is true then store the index of edge in temporary array.
    DO ji=1,p%ncells

      jj1 = g%cells%neighbor_index(listc(ji),1)
      jj2 = g%cells%neighbor_index(listc(ji),2)
      jj3 = g%cells%neighbor_index(listc(ji),3)


      IF (.not.l_counted(jj1)) THEN
        ie_count = ie_count+1
        itmpe(ie_count) = g%cells%edge_index(listc(ji),1)
      ENDIF

      IF (.not.l_counted(jj2)) THEN
        ie_count = ie_count+1
        itmpe(ie_count) = g%cells%edge_index(listc(ji),2)
      ENDIF

      IF (.not.l_counted(jj3)) THEN
        ie_count = ie_count+1
        itmpe(ie_count) = g%cells%edge_index(listc(ji),3)
      ENDIF

      l_counted(listc(ji)) = .true.

    ENDDO

    p%edges%indlist(1:max_rledge,1) = p%nedges
    p%edges%indlist(min_rledge:0,1) = p%nedges+1
    p%edges%indlist(1:max_rledge,2) = 0
    p%edges%indlist(min_rledge:0,2) = p%nedges

    p%edges%start_idx(1:max_rledge,1) = p%nedges
    p%edges%start_idx(min_rledge:0,1) = p%nedges+1
    p%edges%end_idx(1:max_rledge,1) = 0
    p%edges%end_idx(min_rledge:0,1) = p%nedges

    ! Reorder edge index list according to their value of refin_ctrl
    ! The indlist field contains the corresponding start and end indices
    jj1 = 0
    DO irv = 1, max_rledge
      DO ji = 1, p%nedges
        jj = itmpe(ji)
        IF (g%edges%refin_ctrl(jj) == irv) THEN
          jj1 = jj1 + 1
          itmpe2(jj1) = jj
          p%edges%indlist(irv,1) = MIN(p%edges%indlist(irv,1),jj1)
          p%edges%indlist(irv,2) = jj1
          p%edges%start_idx(irv,1) = MIN(p%edges%start_idx(irv,1),jj1)
          p%edges%end_idx(irv,1) = jj1
        ENDIF
      ENDDO
    ENDDO

    ! The following loop used to accomplish the reordering of nest overlap points
    ! (deactivated since Aug. 2013 for restructured nesting flow control)
    irv = 0
    DO ji = 1, p%nedges
      jj = itmpe(ji)
      IF (g%edges%refin_ctrl(jj) <= irv .OR. g%edges%refin_ctrl(jj) > max_rledge) THEN
        jj1 = jj1 + 1
        itmpe2(jj1) = jj
        p%edges%indlist(irv,1) = MIN(p%edges%indlist(irv,1),jj1)
        p%edges%indlist(irv,2) = jj1
        p%edges%start_idx(irv,1) = MIN(p%edges%start_idx(irv,1),jj1)
        p%edges%end_idx(irv,1) = jj1
      ENDIF
    ENDDO


    IF (jj1 /= p%nedges) CALL finish ('create_local_domain', &
      & 'Error in restructuring edge index list')

    DO ji = 1, p%nedges
      itmpe(ji) = itmpe2(ji)
    ENDDO

    DEALLOCATE (itmpe2,l_counted)

    ALLOCATE(l_counted(1:g%nverts))
    l_counted = .false.

    !
    !  count vertices belonging to patch
    !
    p%nverts=0

    DO ji=1,p%nedges

      !get adjacent vertices
      jj1 = g%edges%vertex_index(itmpe(ji),1)
      jj2 = g%edges%vertex_index(itmpe(ji),2)

      IF (.not.l_counted(jj1)) p%nverts = p%nverts+1
      IF (.not.l_counted(jj2)) p%nverts = p%nverts+1

      l_counted(jj1) = .true.
      l_counted(jj2) = .true.

    ENDDO
    !
    !  allocate global index arrays for vertices
    !
    ALLOCATE(itmpv(p%nverts), itmpv2(p%nverts),      &
      & p%verts%indlist(min_rlvert:max_rlvert,2),    &
      & p%verts%start_idx(min_rlvert:max_rlvert,1),  &
      & p%verts%end_idx(min_rlvert:max_rlvert,1)     )
    !
    !  fill global index array for vertices
    !
    iv_count=0
    l_counted = .false.

    DO ji=1,p%nedges

      jj1 = g%edges%vertex_index(itmpe(ji),1)
      jj2 = g%edges%vertex_index(itmpe(ji),2)

      IF (.not.l_counted(jj1)) THEN
        iv_count=iv_count+1
        itmpv(iv_count) = jj1
      ENDIF

      IF (.not.l_counted(jj2)) THEN
        iv_count = iv_count+1
        itmpv(iv_count) = jj2
      ENDIF

      l_counted(jj1) = .true.
      l_counted(jj2) = .true.

    ENDDO

    p%verts%indlist(1:max_rlvert,1) = p%nverts
    p%verts%indlist(min_rlvert:0,1) = p%nverts+1
    p%verts%indlist(1:max_rlvert,2) = 0
    p%verts%indlist(min_rlvert:0,2) = p%nverts

    p%verts%start_idx(1:max_rlvert,1) = p%nverts
    p%verts%start_idx(min_rlvert:0,1) = p%nverts+1
    p%verts%end_idx(1:max_rlvert,1) = 0
    p%verts%end_idx(min_rlvert:0,1) = p%nverts

    ! Reorder vertex index list according to their value of refin_ctrl
    ! The indlist field contains the corresponding start and end indices
    jj1 = 0
    DO irv = 1, max_rlvert
      DO ji = 1, p%nverts
        jj = itmpv(ji)
        IF (g%verts%refin_ctrl(jj) == irv) THEN
          jj1 = jj1 + 1
          itmpv2(jj1) = jj
          p%verts%indlist(irv,1) = MIN(p%verts%indlist(irv,1),jj1)
          p%verts%indlist(irv,2) = jj1
          p%verts%start_idx(irv,1) = MIN(p%verts%start_idx(irv,1),jj1)
          p%verts%end_idx(irv,1) = jj1
        ENDIF
      ENDDO
    ENDDO

    ! The following loop used to accomplish the reordering of nest overlap points
    ! (deactivated since Aug. 2013 for restructured nesting flow control)
    irv = 0
    DO ji = 1, p%nverts
      jj = itmpv(ji)
      IF (g%verts%refin_ctrl(jj) <= irv .OR. g%verts%refin_ctrl(jj) > max_rlvert) THEN
        jj1 = jj1 + 1
        itmpv2(jj1) = jj
        p%verts%indlist(irv,1) = MIN(p%verts%indlist(irv,1),jj1)
        p%verts%indlist(irv,2) = jj1
        p%verts%start_idx(irv,1) = MIN(p%verts%start_idx(irv,1),jj1)
        p%verts%end_idx(irv,1) = jj1
      ENDIF
    ENDDO


    IF (jj1 /= p%nverts)  CALL finish ('create_local_domain', &
      & 'Error in restructuring vertex index list')

    DO ji = 1, p%nverts
      itmpv(ji) = itmpv2(ji)
    ENDDO

    DEALLOCATE (itmpv2)

    ! Add list of start and end indices for cells
    ALLOCATE(p%cells%indlist(min_rlcell:max_rlcell,2),itmpc(p%ncells), &
      & p%cells%start_idx(min_rlcell:max_rlcell,1),    &
      & p%cells%end_idx(min_rlcell:max_rlcell,1)       )

    p%cells%indlist(1:max_rlcell,1) = p%ncells
    p%cells%indlist(min_rlcell:0,1) = p%ncells+1
    p%cells%indlist(1:max_rlcell,2) = 0
    p%cells%indlist(min_rlcell:0,2) = p%ncells

    p%cells%start_idx(1:max_rlcell,1) = p%ncells
    p%cells%start_idx(min_rlcell:0,1) = p%ncells+1
    p%cells%end_idx(1:max_rlcell,1) = 0
    p%cells%end_idx(min_rlcell:0,1) = p%ncells

    jj1 = 0
    DO irv = 1, max_rlcell
      DO ji = 1, p%ncells
        jj = listc(ji)
        IF (g%cells%refin_ctrl(jj) == irv) THEN
          jj1 = jj1 + 1
          itmpc(jj1) = jj
          p%cells%indlist(irv,1) = MIN(p%cells%indlist(irv,1),jj1)
          p%cells%indlist(irv,2) = jj1
          p%cells%start_idx(irv,1) = MIN(p%cells%start_idx(irv,1),jj1)
          p%cells%end_idx(irv,1) = jj1
        ENDIF
      ENDDO
    ENDDO

    ! The following loop used to accomplish the reordering of nest overlap points
    ! (deactivated since Aug. 2013 for restructured nesting flow control)
    irv = 0
    DO ji = 1, p%ncells
      jj = listc(ji)
      IF (g%cells%refin_ctrl(jj) <= irv .OR. g%cells%refin_ctrl(jj) > max_rlcell) THEN
        jj1 = jj1 + 1
        itmpc(jj1) = jj
        p%cells%indlist(irv,1) = MIN(p%cells%indlist(irv,1),jj1)
        p%cells%indlist(irv,2) = jj1
        p%cells%start_idx(irv,1) = MIN(p%cells%start_idx(irv,1),jj1)
        p%cells%end_idx(irv,1) = jj1
      ENDIF
    ENDDO

    !
    !  allocate final  global index arrays for edges, cells, vertices
    !
    ALLOCATE(p%ic(1:p%ncells))
    ALLOCATE(p%ie(1:p%nedges))
    ALLOCATE(p%iv(1:p%nverts))

    !
    !   fill final  global index arrays for edges, cells, vertices:
    !   items belonging to patch are copied from listc and temporary
    !   arrays, items belonging to external halo are copied from halo
    !   arrays
    !
    p%ic(1:p%ncells)=itmpc(1:p%ncells)
    p%ie(1:p%nedges)=itmpe(1:p%nedges)
    p%iv(1:p%nverts)=itmpv(1:p%nverts)



    DO ji=1,p%ncells
      p%iic(p%ic(ji)) = ji
    ENDDO

    DO ji=1,p%nedges
      p%iie(p%ie(ji)) = ji
    ENDDO

    DO ji=1,p%nverts
      p%iiv(p%iv(ji)) = ji
    ENDDO


    !
    !   deallocate auxiliary arrays
    !

    DEALLOCATE(itmpc,itmpe,itmpv)
    DEALLOCATE(l_counted)


    !
    ! Add grid information to the patch
    !
    !
    ! Fill in the cells
    !
    nc_int = p%ncells

    ist = 0
    ALLOCATE(p%cells%idx             (nc_int  ),stat=istat)
    ist = ist + ABS(istat)
    ALLOCATE(p%cells%parent_index    (nc_int  ),stat=istat)
    ist = ist + ABS(istat)
    ALLOCATE(p%cells%child_index     (nc_int,4),stat=istat)
    ist = ist + ABS(istat)
    ALLOCATE(p%cells%child_id        (nc_int  ),stat=istat)
    ist = ist + ABS(istat)
    ALLOCATE(p%cells%phys_id         (nc_int  ),stat=istat)
    ist = ist + ABS(istat)
    ALLOCATE(p%cells%neighbor_index  (nc_int,3),stat=istat)
    ist = ist + ABS(istat)
    ALLOCATE(p%cells%edge_index      (nc_int,3),stat=istat)
    ist = ist + ABS(istat)
    ALLOCATE(p%cells%vertex_index    (nc_int,3),stat=istat)
    ist = ist + ABS(istat)
    ALLOCATE(p%cells%edge_orientation(nc_int,3),stat=istat)
    ist = ist + ABS(istat)
    ALLOCATE(p%cells%center          (nc_int  ),stat=istat)
    ist = ist + ABS(istat)
    ALLOCATE(p%cells%area            (nc_int  ),stat=istat)
    ist = ist + ABS(istat)
    ALLOCATE(p%cells%refin_ctrl      (nc_int  ),stat=istat)
    ist = ist + ABS(istat)

    IF (ist>0) THEN
      WRITE (message_text,'(a,i2)') 'Problem in allocating local patch cells'
      CALL finish ('create_local_domain', message_text)
    ENDIF

    DO j = 1, nc_int

      p%cells%idx(j) = j

      g_cell = p%ic(j)

      ! These indices are set later in subroutine complete_local_patch
      p%cells%parent_index(j) = 0

      p%cells%child_index(j,1:4) = 0

      l2g = g%cells%neighbor_index(g_cell,1)
      p%cells%neighbor_index(j     ,1) = p%iic(l2g)
      l2g = g%cells%neighbor_index(g_cell,2)
      p%cells%neighbor_index(j     ,2) = p%iic(l2g)
      l2g = g%cells%neighbor_index(g_cell,3)
      p%cells%neighbor_index(j     ,3) = p%iic(l2g)

      l2g = g%cells%edge_index(g_cell,1)
      p%cells%edge_index(j     ,1) = p%iie(l2g)
      l2g = g%cells%edge_index(g_cell,2)
      p%cells%edge_index(j     ,2) = p%iie(l2g)
      l2g = g%cells%edge_index(g_cell,3)
      p%cells%edge_index(j     ,3) = p%iie(l2g)

      l2g = g%cells%vertex_index(g_cell,1)
      p%cells%vertex_index(j     ,1) = p%iiv(l2g)
      l2g = g%cells%vertex_index(g_cell,2)
      p%cells%vertex_index(j     ,2) = p%iiv(l2g)
      l2g = g%cells%vertex_index(g_cell,3)
      p%cells%vertex_index(j     ,3) = p%iiv(l2g)

      p%cells%refin_ctrl      (j  ) = g%cells%refin_ctrl      (g_cell  )
      p%cells%child_id        (j  ) = g%cells%child_id        (g_cell  )
      p%cells%phys_id         (j  ) = g%cells%phys_id         (g_cell  )

      p%cells%edge_orientation(j,1) = g%cells%edge_orientation(g_cell,1)
      p%cells%edge_orientation(j,2) = g%cells%edge_orientation(g_cell,2)
      p%cells%edge_orientation(j,3) = g%cells%edge_orientation(g_cell,3)
      p%cells%center          (j  ) = g%cells%center          (g_cell  )
      p%cells%area            (j  ) = g%cells%area            (g_cell  )

    ENDDO
    !
    ! Fill in the edges
    !
    ne_int = p%nedges

    ist = 0
    ALLOCATE(p%edges%idx               (ne_int  ),stat=istat)
    ist = ist + ABS(istat)
    ALLOCATE(p%edges%cell_index        (ne_int,2),stat=istat)
    ist = ist + ABS(istat)
    ALLOCATE(p%edges%vertex_index      (ne_int,2),stat=istat)
    ist = ist + ABS(istat)
    ALLOCATE(p%edges%system_orientation(ne_int  ),stat=istat)
    ist = ist + ABS(istat)
    ALLOCATE(p%edges%center            (ne_int  ),stat=istat)
    ist = ist + ABS(istat)
    ALLOCATE(p%edges%primal_normal     (ne_int  ),stat=istat)
    ist = ist + ABS(istat)
    ALLOCATE(p%edges%dual_normal       (ne_int  ),stat=istat)
    ist = ist + ABS(istat)
    ALLOCATE(p%edges%primal_edge_length(ne_int  ),stat=istat)
    ist = ist + ABS(istat)
    ALLOCATE(p%edges%dual_edge_length  (ne_int  ),stat=istat)
    ist = ist + ABS(istat)
    ALLOCATE(p%edges%edge_vert_length  (ne_int,2),stat=istat)
    ist = ist + ABS(istat)
    ALLOCATE(p%edges%edge_cell_length  (ne_int,2),stat=istat)
    ist = ist + ABS(istat)
    ALLOCATE(p%edges%refin_ctrl        (ne_int  ),stat=istat)
    ist = ist + ABS(istat)
    ALLOCATE(p%edges%parent_index      (ne_int  ),stat=istat)
    ist = ist + ABS(istat)
    ALLOCATE(p%edges%child_index       (ne_int,4),stat=istat)
    ist = ist + ABS(istat)
    ALLOCATE(p%edges%child_id          (ne_int  ),stat=istat)
    ist = ist + ABS(istat)
    ALLOCATE(p%edges%phys_id           (ne_int  ),stat=istat)
    ist = ist + ABS(istat)

    IF (ist>0) THEN
      WRITE (message_text,'(a)') 'Problem in allocating local patch edges'
      CALL finish ('create_local_domain', message_text)
    ENDIF

    DO j = 1, ne_int

      p%edges%idx(j) = j

      g_edge = p%ie(j)

      l2g = g%edges%cell_index(g_edge,1)
      p%edges%cell_index(j     ,1) = p%iic(l2g)
      l2g = g%edges%cell_index(g_edge,2)
      p%edges%cell_index(j     ,2) = p%iic(l2g)


      l2g = g%edges%vertex_index(g_edge,1)
      p%edges%vertex_index(j     ,1) = p%iiv(l2g)
      l2g = g%edges%vertex_index(g_edge,2)
      p%edges%vertex_index(j     ,2) = p%iiv(l2g)

      p%edges%refin_ctrl        (j) = g%edges%refin_ctrl        (g_edge)
      p%edges%child_id          (j) = g%edges%child_id          (g_edge)
      p%edges%phys_id           (j) = g%edges%phys_id           (g_edge)
      ! These indices are set later in subroutine complete_local_patch
      p%edges%parent_index(j)       = 0
      p%edges%child_index(j,1:4)    = 0

      p%edges%system_orientation(j) = g%edges%system_orientation(g_edge)
      p%edges%center            (j) = g%edges%center            (g_edge)
      p%edges%primal_normal     (j) = g%edges%primal_normal     (g_edge)
      p%edges%dual_normal       (j) = g%edges%dual_normal       (g_edge)
      p%edges%primal_edge_length(j) = g%edges%primal_edge_length(g_edge)
      p%edges%dual_edge_length  (j) = g%edges%dual_edge_length  (g_edge)
      p%edges%edge_vert_length(j,:) = g%edges%edge_vert_length(g_edge,:)
      p%edges%edge_cell_length(j,:) = g%edges%edge_cell_length(g_edge,:)


    ENDDO
    !
    ! Fill in the verts
    !
    nv = p%nverts

    ist = 0
    ALLOCATE(p%verts%idx             (nv  ),stat=istat)
    ist = ist + ABS(istat)
    ALLOCATE(p%verts%neighbor_index  (nv,6),stat=istat)
    ist = ist + ABS(istat)
    ALLOCATE(p%verts%cell_index      (nv,6),stat=istat)
    ist = ist + ABS(istat)
    ALLOCATE(p%verts%edge_index      (nv,6),stat=istat)
    ist = ist + ABS(istat)
    ALLOCATE(p%verts%edge_orientation(nv,6),stat=istat)
    ist = ist + ABS(istat)
    ALLOCATE(p%verts%vertex          (nv  ),stat=istat)
    ist = ist + ABS(istat)
    ALLOCATE(p%verts%dual_area       (nv  ),stat=istat)
    ist = ist + ABS(istat)
    ALLOCATE(p%verts%refin_ctrl      (nv  ),stat=istat)
    ist = ist + ABS(istat)
    ALLOCATE(p%verts%child_id        (nv  ),stat=istat)
    ist = ist + ABS(istat)
    ALLOCATE(p%verts%phys_id         (nv  ),stat=istat)
    ist = ist + ABS(istat)

    IF (ist>0) THEN
      WRITE (message_text,'(a)') 'Problem in allocating local patch verts'
      CALL finish ('create_local_domain', message_text)
    ENDIF

    DO j = 1, nv

      p%verts%idx(j) = j

      g_vert = p%iv(j)
      !
      ! note: the dummy edges, neighbor_vertices and cells get the index 0
      !
      DO i1 = 1, 6

        l2g = g%verts%neighbor_index(g_vert,i1)
        IF (l2g > 0) THEN
          p%verts%neighbor_index(j     ,i1) = p%iiv(l2g)
        ELSE
          p%verts%neighbor_index(j     ,i1) = 0
        ENDIF

        l2g = g%verts%cell_index(g_vert,i1)
        IF (l2g > 0) THEN
          p%verts%cell_index(j     ,i1) = p%iic(l2g)
        ELSE
          p%verts%cell_index(j     ,i1) = 0
        ENDIF

        l2g = g%verts%edge_index(g_vert,i1)
        IF (l2g > 0) THEN
          p%verts%edge_index(j     ,i1) = p%iie(l2g)
        ELSE
          p%verts%edge_index(j     ,i1) = 0
        ENDIF

      ENDDO

      p%verts%edge_orientation(j,1:6) = g%verts%edge_orientation(g_vert,1:6)
      p%verts%dual_area       (j    ) = g%verts%dual_area       (g_vert    )
      p%verts%vertex          (j    ) = g%verts%vertex          (g_vert    )
      p%verts%refin_ctrl      (j    ) = g%verts%refin_ctrl      (g_vert    )
      p%verts%child_id        (j    ) = g%verts%child_id        (g_vert    )
      p%verts%phys_id         (j    ) = g%verts%phys_id         (g_vert    )

    ENDDO

  END SUBROUTINE create_local_domain

  !-------------------------------------------------------------------------
  !
  !
  !>
  !!  This subroutine defines the parent index and child index fields for the.
  !!
  !!  This subroutine defines the parent index and child index fields for the
  !!  patches generated in make_local_patch. A separate routine is required for
  !!  this because the patches need be created for all refinement levels before the
  !!  parent/child indices can be transferred from the former global grid to the
  !!  patches.
  !!
  !! @par Revision History
  !! Initial version: Guenther Zaengl, DWD, March 2008
  !! @par
  !! Guenther Zaengl, DWD, 2009-07-22:
  !!   Updated version for generalized grid refinement
  !!
  SUBROUTINE complete_index_lists(p, pp, pc, ichid, g, gp, gc, child_exist, parent_exist)
    !

    !
    ! patches of current level, parent level and child level
    !
    TYPE(t_patch) :: p, pp, pc
    !
    ! grids on which patches live
    !
    TYPE(t_grid) , TARGET, INTENT(in)  :: g, gp, gc
    !
    LOGICAL, INTENT(in) :: child_exist,parent_exist

    INTEGER, INTENT(in) :: ichid    ! child domain ID

    INTEGER :: j, k, ig, j1, igpar, igchi(4)

    INTEGER, ALLOCATABLE, DIMENSION(:) :: index_map

    ! Set indices for parent cells and edges if a parent level exists
    IF (parent_exist) THEN

      ALLOCATE (index_map(gp%ncells))
      index_map = 0
      DO j1 = 1, pp%ncells
        index_map(pp%ic(j1)) = j1
      ENDDO

      DO j = 1, p%ncells
        ig = p%ic(j)
        igpar = g%cells%parent_index(ig)
        p%cells%parent_index(j) = index_map(igpar)
      ENDDO

      DEALLOCATE (index_map)
      ALLOCATE (index_map(gp%nedges))
      index_map = 0

      p%edges%parent_index(:) = -1 ! Initialization

      DO j1 = 1, pp%nedges
        index_map(pp%ie(j1)) = j1
      ENDDO

      DO j = 1, p%nedges
        ig = p%ie(j)
        igpar = g%edges%parent_index(ig)
        p%edges%parent_index(j) = index_map(igpar)
      ENDDO
      DEALLOCATE (index_map)

    ENDIF

    ! Set indices for child cells and edges if a child level exists
    IF (child_exist) THEN

      ALLOCATE (index_map(gc%ncells))
      index_map = 0
      DO j1 = 1, pc%ncells
        index_map(pc%ic(j1)) = j1
      ENDDO

      DO j = 1, p%ncells
        ig = p%ic(j)
        igchi(1:4) = g%cells%child_index(ig,1:4)
        DO k = 1, 4
          IF (index_map(igchi(k)) /= 0)  &
            & p%cells%child_index(j,k) = index_map(igchi(k))
        ENDDO
      ENDDO
      DEALLOCATE (index_map)

      ! The edge child index requires a special treatment because the third and
      ! fourth can be negative, denoting that the orientation of the edge of
      ! the inner small triangle is opposite to that of the main triangle edge
      ALLOCATE (index_map(0:gc%nedges))
      index_map = 0
      DO j1 = 1, pc%nedges
        index_map(pc%ie(j1)) = j1
      ENDDO

      DO j = 1, p%nedges
        ig = p%ie(j)
        igchi(1:4) = g%edges%child_index(ig,1:4)
        DO k = 1, 4
          IF (index_map(ABS(igchi(k))) /= 0) &
            & p%edges%child_index(j,k) = index_map(ABS(igchi(k)))*sign(1,igchi(k))
        ENDDO
      ENDDO
      DEALLOCATE (index_map)

    ENDIF

    IF (parent_exist) THEN

      ! Consistency check 1a: For child cells flagged with 1 or 2,
      ! the parent cells should have -1

      DO j = 1, p%ncells
        IF ((p%cells%refin_ctrl(j) == 1) .or. (p%cells%refin_ctrl(j) == 2)) THEN
          j1 = p%cells%parent_index(j)
          IF (pp%cells%refin_ctrl(j1) /= -1)  CALL finish &
            & ('complete_index_lists', 'refin_ctrl consistency error 1a')

          ! Consistency check 1b: For child cells flagged with 3 or 4,
          ! the parent cells should have -2

        ELSE IF ((p%cells%refin_ctrl(j) == 3) .or. (p%cells%refin_ctrl(j) == 4)) THEN
          j1 = p%cells%parent_index(j)
          IF (pp%cells%refin_ctrl(j1) /= -2) CALL finish &
            & ('complete_index_lists', 'refin_ctrl consistency error 1b')

          ! Consistency check 1c: For child cells flagged with 5 or 6,
          ! the parent cells should have -3

        ELSE IF ((p%cells%refin_ctrl(j) == 5) .or. (p%cells%refin_ctrl(j) == 6)) THEN
          j1 = p%cells%parent_index(j)
          IF (pp%cells%refin_ctrl(j1) /= -3) CALL finish &
            & ('complete_index_lists', 'refin_ctrl consistency error 1c')

          ! Consistency check 1d: For child cells flagged with 7 or 8,
          ! the parent cells should have -4

        ELSE IF ((p%cells%refin_ctrl(j) == 7) .or. (p%cells%refin_ctrl(j) == 8)) THEN
          j1 = p%cells%parent_index(j)
          IF (pp%cells%refin_ctrl(j1) /= -4) CALL finish &
            & ('complete_index_lists', 'refin_ctrl consistency error 1d')
        ENDIF
      ENDDO

    ENDIF


    IF (child_exist) THEN

      ! Consistency check 2a: For parent cells flagged with -1,
      ! the child cells should have 1 or 2

      DO j = 1, p%ncells
        IF ((p%cells%refin_ctrl(j) == -1) .and. (p%cells%child_id(j) == ichid)) THEN
          DO k=1,4
            j1 = p%cells%child_index(j,k)
            IF ((pc%cells%refin_ctrl(j1) < 1).or.(pc%cells%refin_ctrl(j1) > 2)) &
              & CALL finish ('complete_index_lists', &
              & 'refin_ctrl consistency error 2a')
          ENDDO

          ! Consistency check 2b: For parent cells flagged with -2,
          ! the child cells should have 3 or 4

        ELSE IF ((p%cells%refin_ctrl(j) == -2) .and. (p%cells%child_id(j) == ichid)) THEN
          DO k=1,4
            j1 = p%cells%child_index(j,k)
            IF ((pc%cells%refin_ctrl(j1) < 3).or.(pc%cells%refin_ctrl(j1) > 4)) &
              & CALL finish ('complete_index_lists', &
              & 'refin_ctrl consistency error 2b')
          ENDDO

        ENDIF
      ENDDO

      ! Consistency check 3: Are child edges set correctly?

      DO j = 1, p%nedges
        IF ((p%edges%refin_ctrl(j) == -1) .and. (p%edges%child_id(j) == ichid)) THEN
          DO k=1,3
            j1 = p%edges%child_index(j,k)
            IF (j1 == 0) CALL finish ('complete_index_lists', &
              & 'child edges at refin_ctrl=-1 incorrect')
          ENDDO
        ELSE IF ((p%edges%refin_ctrl(j) < -1) .and. (p%edges%child_id(j) == ichid)) THEN
          DO k=1,4
            j1 = p%edges%child_index(j,k)
            IF (j1 == 0) CALL finish ('complete_index_lists', &
              & 'child edges at refin_ctrl<-1 incorrect')
          ENDDO
        ENDIF
        DO k=1,2
          j1 = p%edges%child_index(j,k)
          IF (j1 < 0) CALL finish ('complete_index_lists', &
            & 'negative child edge index at parent edge')
        ENDDO
      ENDDO

    ENDIF

  END SUBROUTINE complete_index_lists


  !-------------------------------------------------------------------------
  !
  !
  !>
  !!  This subroutine maps the parent/child index relationships between the global.
  !!
  !!  This subroutine maps the parent/child index relationships between the global
  !!  computational grid and the next lower grid level. This is needed for search
  !!  algorithms such as used in external parameter programs.
  !!
  !! @par Revision History
  !! Initial version: Guenther Zaengl, DWD, 2009-07-27
  !! @par
  !! patch at global level
  !!
  SUBROUTINE map_indices(p, gp)
    !
    TYPE(t_patch) :: p
    !
    ! grids at global level and its paret grid
    !
    TYPE(t_grid) :: gp
    !
    INTEGER :: j, k

    DO j = 1, gp%ncells
      DO k = 1, 4
        gp%cells%child_index(j,k) = p%iic(gp%cells%child_index(j,k))
        p%cells%parent_index(gp%cells%child_index(j,k)) = j
      ENDDO
    ENDDO

  END SUBROUTINE map_indices

  !-------------------------------------------------------------------------
  !
  !>
  !! Creates the parent patch for the global domain (for calculations on reduced grid)
  !! and sets the parent/child relations for the cells in pp AND p.
  !! The parent/child relations for the edges are currently NOT set, since they
  !! are not present.
  !!
  !! @par Revision History
  !! Initial version: Rainer Johanni, 2010-10-29

  SUBROUTINE create_global_parent_domain(p, pp, gp)

    TYPE(t_patch) :: p, pp
    TYPE(t_grid), INTENT(IN), TARGET :: gp

    INTEGER :: j, k


    pp%pgrid=>gp

    pp%level=gp%level

    pp%ncells=gp%ncells
    pp%nedges=gp%nedges
    pp%nverts=gp%nverts

    CALL allocate_patch(pp)

    pp%id = 0
    pp%parent_id = -1
    pp%child_id(:) = 0
    pp%child_id(1) = p%id ! i.e. 1
    pp%n_childdom  = 1
    pp%n_phys_dom  = 1
    pp%max_childdom = max_childdom

    ! cells
    !------

    pp%cells%indlist(:,:) = 0 ! not needed

    pp%cells%start_idx(:,1) = 1
    pp%cells%end_idx(:,1)   = 0
    pp%cells%start_idx(min_rlcell:min_rlcell_int-1,1) = pp%ncells+1
    pp%cells%end_idx(min_rlcell:min_rlcell_int,1) = pp%ncells

    DO j = 1, pp%ncells
      pp%cells%idx(j) = j
      DO k = 1, 4
        pp%cells%child_index(j,k) = p%iic(gp%cells%child_index(j,k))
        p%cells%parent_index(pp%cells%child_index(j,k)) = j
      ENDDO
    ENDDO

    pp%cells%parent_index    (:  ) = -1
    pp%cells%child_id        (:  ) = p%id ! i.e. 1
    pp%cells%neighbor_index  (:,:) = gp%cells%neighbor_index  (:,:)
    pp%cells%edge_index      (:,:) = gp%cells%edge_index      (:,:)
    pp%cells%vertex_index    (:,:) = gp%cells%vertex_index    (:,:)
    pp%cells%edge_orientation(:,:) = gp%cells%edge_orientation(:,:)
    pp%cells%center          (:  ) = gp%cells%center          (:  )
    pp%cells%area            (:  ) = gp%cells%area            (:  )
    pp%cells%refin_ctrl      (:  ) = min_rlcell_int
    pp%cells%phys_id         (:  ) = 0 ! not needed

    ! edges
    !------

    pp%edges%indlist(:,:) = 0 ! not needed

    pp%edges%start_idx(:,1) = 1
    pp%edges%end_idx(:,1)   = 0
    pp%edges%start_idx(min_rledge:min_rledge_int-1,1) = pp%nedges+1
    pp%edges%end_idx(min_rledge:min_rledge_int,1) = pp%nedges

    DO j = 1, pp%nedges
      pp%edges%idx(j) = j
    ENDDO

    p%edges%parent_index       (:)   = 0 ! set below
    pp%edges%child_index       (:,:) = 0 ! set below

    pp%edges%parent_index      (:)   = -1
    pp%edges%child_id          (:)   = p%id ! i.e. 1
    pp%edges%cell_index        (:,:) = gp%edges%cell_index        (:,:)
    pp%edges%vertex_index      (:,:) = gp%edges%vertex_index      (:,:)
    pp%edges%system_orientation(:)   = gp%edges%system_orientation(:)
    pp%edges%center            (:)   = gp%edges%center            (:)
    pp%edges%primal_normal     (:)   = gp%edges%primal_normal     (:)
    pp%edges%dual_normal       (:)   = gp%edges%dual_normal       (:)
    pp%edges%primal_edge_length(:)   = gp%edges%primal_edge_length(:)
    pp%edges%dual_edge_length  (:)   = gp%edges%dual_edge_length  (:)
    pp%edges%edge_vert_length  (:,:) = gp%edges%edge_vert_length  (:,:)
    pp%edges%edge_cell_length  (:,:) = gp%edges%edge_cell_length  (:,:)
    pp%edges%refin_ctrl        (:)   = min_rledge_int
    pp%edges%phys_id           (:)   = 0 ! not needed

    ! Now set parent/child relationship
    CALL get_edge_parent_child_relations(p, pp)


    ! vertices
    !---------

    pp%verts%indlist(:,:) = 0 ! not needed

    pp%verts%start_idx(:,1) = 1
    pp%verts%end_idx(:,1)   = 0
    pp%verts%start_idx(min_rlvert:min_rlvert_int-1,1) = pp%nverts+1
    pp%verts%end_idx(min_rlvert:min_rlvert_int,1) = pp%nverts

    DO j = 1, pp%nverts
      pp%verts%idx(j) = j
    ENDDO

    pp%verts%child_id        (:)   = -1

    pp%verts%neighbor_index  (:,:) = gp%verts%neighbor_index  (:,:)
    pp%verts%cell_index      (:,:) = gp%verts%cell_index      (:,:)
    pp%verts%edge_index      (:,:) = gp%verts%edge_index      (:,:)
    pp%verts%edge_orientation(:,:) = gp%verts%edge_orientation(:,:)
    pp%verts%dual_area       (:)   = gp%verts%dual_area       (:)
    pp%verts%vertex          (:)   = gp%verts%vertex          (:)
    pp%verts%refin_ctrl      (:)   = min_rlvert_int
    pp%verts%phys_id         (:)   = 0 ! not needed

  END SUBROUTINE create_global_parent_domain

  !-------------------------------------------------------------------------
  !
  !> Calculate edge parent child relations for two patches (child and parent)
  !! Sets child_index in pp and parent_index in p
  !! For consistency, the same method as in setup_index_lists() is used
  !! but it is applied to the whole patch
  !!
  !! @par Revision History
  !! Initial version: Rainer Johanni, 2010-10-29

  SUBROUTINE get_edge_parent_child_relations(p, pp)

    TYPE(t_patch), TARGET, INTENT(INOUT) :: p, pp

    INTEGER :: i, j, k, i1, i2, i_edp(3), i_edc(4,3), nc, nct, nct2

    REAL(wp) :: auxp(4), auxc(4), mindist, dist(3,3,4), orient(3,3,4), auxpc(3,3), auxcc(4,3,3)

    TYPE(t_cartesian_coordinates) :: aux_epc, aux_ecc

    TYPE(t_grid_cells),POINTER :: cc =>null()
    TYPE(t_grid_cells),POINTER :: ccp =>null()
    TYPE(t_grid_edges),POINTER :: ce =>null()
    TYPE(t_grid_edges),POINTER :: cep =>null()

    ! Calculate parent and child indices for edges

    cc=>p%cells
    ce=>p%edges
    ccp=>pp%cells
    cep=>pp%edges

    ! Logic of the child edge indices: For a given edge, children 1 and 2 are
    ! the sub-triangle edges coinciding with the parent edge; children 3 and 4
    ! are the edges of the inner sub-triangles having (roughly) the same
    ! orientation as the parent edge

    DO i = 1, pp%ncells
      DO j = 1, 3
        i_edp(j) = ccp%edge_index(i,j)
        auxp(1) = cep%center(i_edp(j))%lat
        auxp(2) = cep%center(i_edp(j))%lon
        auxp(3) = cep%primal_normal(i_edp(j))%v1
        auxp(4) = cep%primal_normal(i_edp(j))%v2
        CALL gvec2cvec(auxp(3),auxp(4),auxp(2),auxp(1),&
          & auxpc(j,1),auxpc(j,2),auxpc(j,3))
      ENDDO
      DO k = 1, 4
        i1 = ccp%child_index(i,k)
        DO j = 1, 3
          i_edc(k,j) = cc%edge_index(i1,j)
          auxc(1) = ce%center(i_edc(k,j))%lat
          auxc(2) = ce%center(i_edc(k,j))%lon
          auxc(3) = ce%primal_normal(i_edc(k,j))%v1
          auxc(4) = ce%primal_normal(i_edc(k,j))%v2
          CALL gvec2cvec(auxc(3),auxc(4),auxc(2),auxc(1),&
            & auxcc(k,j,1),auxcc(k,j,2),auxcc(k,j,3))
        ENDDO
      ENDDO
      DO i1 = 1, 3
        DO k = 1, 4
          DO i2 = 1, 3
            aux_epc         = gc2cc(cep%center(i_edp(i1)))
            aux_ecc         = gc2cc(ce%center(i_edc(k,i2)))
            dist(i1,i2,k)   = arc_length(aux_epc,aux_ecc)
            orient(i1,i2,k) = DOT_PRODUCT(auxpc(i1,:),auxcc(k,i2,:))* &
              & ABS(DOT_PRODUCT(auxpc(i1,:),auxcc(k,i2,:)))/ &
              & DOT_PRODUCT(auxpc(i1,:),auxpc(i1,:))/ &
              & DOT_PRODUCT(auxcc(k,i2,:),auxcc(k,i2,:))
          ENDDO
        ENDDO
      ENDDO
      mindist = MINVAL(dist)
      nct = 0
      nct2 = 0
      DO i1 = 1, 3
        nc = 0
        DO k = 1, 4
          DO i2 = 1, 3
            IF ( (dist(i1,i2,k) <= 1.5_wp*mindist) .and. &
              & (orient(i1,i2,k) > 0.99_wp) ) THEN
              nc = nc+1
              nct = nct+1
              cep%child_index(i_edp(i1),nc) = i_edc(k,i2)
              ce%parent_index(i_edc(k,i2)) = i_edp(i1)
              IF (nc > 2) CALL finish ('get_edge_parent_child_relations', &
                & 'Error 1 in calculation of parent/child edge indices')
            ELSE IF ( (dist(i1,i2,k) > 1.5_wp*mindist) .and. &
              & (orient(i1,i2,k) > 0.99_wp) ) THEN
              nct2 = nct2+1
              IF (cep%child_index(i_edp(i1),3) == 0) THEN
                cep%child_index(i_edp(i1),3) = i_edc(k,i2)
                ce%parent_index(i_edc(k,i2)) = i_edp(i1)
              ELSE IF (cep%child_index(i_edp(i1),3) /= i_edc(k,i2)) THEN
                cep%child_index(i_edp(i1),4) = i_edc(k,i2)
                ce%parent_index(i_edc(k,i2)) = i_edp(i1)
              ENDIF
            ELSE IF ( (dist(i1,i2,k) > 1.5_wp*mindist) .and.  &
              & (orient(i1,i2,k) < -0.99_wp) ) THEN
              nct2 = nct2+1
              IF (cep%child_index(i_edp(i1),3) == 0) THEN
                cep%child_index(i_edp(i1),3) = i_edc(k,i2)
                ce%parent_index(i_edc(k,i2)) = i_edp(i1)
              ELSE IF (cep%child_index(i_edp(i1),3) /= i_edc(k,i2)) THEN
                cep%child_index(i_edp(i1),4) = i_edc(k,i2)
                ce%parent_index(i_edc(k,i2)) = i_edp(i1)
              ENDIF
            ENDIF
          ENDDO
        ENDDO
      ENDDO
      ! Check sum of child edges
      IF (nct /= 6) CALL finish ('get_edge_parent_child_relations', &
        & 'Error 2 in calculation of parent/child edge indices')
      ! There are actually 3 inner edges, but they are counted twice
      IF (nct2 /= 6) CALL finish ('get_edge_parent_child_relations', &
        & 'Error 3 in calculation of parent/child edge indices')
    ENDDO

  END SUBROUTINE get_edge_parent_child_relations

  !-------------------------------------------------------------------------
  !BOP
  !
  ! !IROUTINE: write_patch
  !
  ! !SUBROUTINE INTERFACE:
  SUBROUTINE write_patch(p, uuid_grid, uuid_parent, uuid_child)
    !
    ! !DESCRIPTION:
    !
    !  Writes patch data to file. New version for netcdf output
    !
    ! !REVISION HISTORY:
    !  Guenther Zaengl, DWD, October 2008
    !  USES:
    TYPE(t_patch), POINTER, INTENT(INOUT) :: p

    CHARACTER(len=uuid_string_length), INTENT(IN) :: uuid_grid, uuid_parent, uuid_child(:)

    !---local variables
    CHARACTER(LEN=filename_max) :: filename, filename_grfinfo
    CHARACTER(LEN=1) :: child_id

    TYPE(t_grid), POINTER :: g =>null()
    INTEGER :: i, j, i_lev

!!$    INTEGER, PARAMETER :: patch_unit = 89

    ! Copied from output_grid...

    INTEGER :: i_nc, i_ne, i_nv
    INTEGER :: dom_id, parent_dom_id
    INTEGER :: old_mode
    INTEGER :: ncid
    INTEGER :: dimids(2)

    INTEGER :: dim_ncell, dim_nvertex, dim_nedge
    INTEGER :: dim_two, dim_cell_refine, dim_edge_refine, dim_vert_refine
    INTEGER :: dim_nvertex_per_cell, dim_ncells_per_edge, dim_nedges_per_vertex, &
      & dim_nchilds_per_cell, dim_list, dim_nchdom

    INTEGER :: varid_clon, varid_clat, varid_clonv, varid_clatv
    INTEGER :: varid_vlon, varid_vlat, varid_vlonv, varid_vlatv
    INTEGER :: varid_elon, varid_elat, varid_elonv, varid_elatv

    INTEGER :: varid_carea, varid_varea, varid_earea

    INTEGER :: varid1, varid2, varid3, varid4, varid5, varid6, varid7, varid8
    INTEGER :: varid9, varid10, varid11, varid12, varid13, varid14, varid15, varid16
    INTEGER :: varid17, varid18, varid19, varid20, varid21, varid22, varid23, varid24
    INTEGER :: varid25, varid26, varid27, varid28, varid29, varid30, varid31, varid32
    INTEGER :: varid33, varid34, varid35, varid36, varid37, varid38, varid39, varid40
    INTEGER :: varid41, varid42, varid43, varid44, varid45, varid46, varid47, varid48
    INTEGER :: varid49, varid50

    INTEGER :: ifs2icon_cell, ifs2icon_edge, ifs2icon_vertex

    INTEGER :: istat
    CHARACTER(LEN=256) :: command_line

    INTEGER :: command_line_len

    REAL(wp), ALLOCATABLE :: zv2d(:,:), zv2dx(:,:), zv2dy(:,:)

    REAL(wp) :: swap(4)
    INTEGER :: str_idx, end_idx
    INTEGER :: current_number_of_grid_used

    !EOP
    !-------------------------------------------------------------------------
    !BOC

    g=>p%pgrid

    i_nc = p%ncells
    i_ne = p%nedges
    i_nv = p%nverts
    i_lev = p%level

    dom_id        = p%id
    parent_dom_id = p%parent_id

    IF (write_hierarchy == 0) THEN ! no radiation grid
      current_number_of_grid_used = number_of_grid_used(dom_id)
    ELSE IF (write_hierarchy == 1) THEN ! radiation grid is generated - number of entries in number_of_grid_used must be n_dom+1
      current_number_of_grid_used = number_of_grid_used(dom_id+1)
    ELSE ! the special mode for writing the grid hierarchy back to level 0 does not allow for setting the number_of_grid_used
      current_number_of_grid_used = 0
    ENDIF


    IF ( outname_style==1 ) THEN   ! Default naming convention
      WRITE(filename,'(a,i0,2(a,i2.2),a)') &
        & 'iconR', grid_root, 'B', p%level, '_DOM', dom_id, '.nc'
      IF (lsep_gridref_info) THEN ! generate separate files for information on parent-child connectivities
      WRITE(filename_grfinfo,'(a,i0,2(a,i2.2),a)') &
        & 'iconR', grid_root, 'B', p%level, '_DOM', dom_id, '-grfinfo.nc'
      ENDIF
    ELSE                           ! DWD naming convention
      IF (dom_id==0) THEN
        WRITE(filename,'(a,i4.4,2(a,i2.2),a)') &
          & 'icon_grid_', current_number_of_grid_used, '_R', grid_root, 'B', p%level, '_R.nc'
      ELSE IF (dom_id==1) THEN
        WRITE(filename,'(a,i4.4,2(a,i2.2),a)') &
          & 'icon_grid_', current_number_of_grid_used, '_R', grid_root, 'B', p%level, '_G.nc'
      ELSE
        WRITE(filename,'(a,i4.4,3(a,i2.2),a)') &
          & 'icon_grid_', current_number_of_grid_used, '_R', grid_root, 'B', p%level, '_N', dom_id, '.nc'
      ENDIF
      IF (lsep_gridref_info) THEN ! generate separate files for information on parent-child connectivities
        IF (dom_id==0) THEN
          WRITE(filename_grfinfo,'(a,i4.4,2(a,i2.2),a)') &
            & 'icon_grid_', current_number_of_grid_used, '_R', grid_root, 'B', p%level, '_R-grfinfo.nc'
        ELSE IF (dom_id==1) THEN
          WRITE(filename_grfinfo,'(a,i4.4,2(a,i2.2),a)') &
            & 'icon_grid_', current_number_of_grid_used, '_R', grid_root, 'B', p%level, '_G-grfinfo.nc'
        ELSE
          WRITE(filename_grfinfo,'(a,i4.4,3(a,i2.2),a)') &
            & 'icon_grid_', current_number_of_grid_used, '_R', grid_root, 'B', p%level, '_N', dom_id, '-grfinfo.nc'
        ENDIF
      ENDIF
    ENDIF


    WRITE(message_text,'(a,a)') 'Write ICON grid file: ', TRIM(filename)
    CALL message ('', message_text)


    CALL GET_COMMAND(command_line, command_line_len, istat)

    CALL nf(nf_set_default_format(nf_format_64bit, old_mode))
    CALL nf(nf_create(TRIM(filename), nf_clobber, ncid))
    CALL nf(nf_set_fill(ncid, nf_nofill, old_mode))

    ! Set global attributes
    !
    CALL set_global_attributes(.FALSE.)

    ! Set dimensions
    !
    CALL set_dimensions()

    ! Grid variables
    !
    ! Define variables for basic grid and geometry information
    !
    CALL define_standard_variables()

    ! Define variables for parent-child connectivity information
    !
    IF (.NOT. lsep_gridref_info) CALL define_connectivity_variables()

    ! end of definition part
    !
    CALL nf(nf_enddef(ncid))

    ! Write variables for basic grid and geometry information
    !
    CALL write_standard_variables()

    ! Write variables for parent-child connectivity information
    !
    IF (.NOT. lsep_gridref_info) CALL write_connectivity_variables()

    ! close file
    CALL nf(nf_close(ncid))

    ! Generate separate file with connectivity info if requested
    IF (lsep_gridref_info) THEN

      WRITE(message_text,'(a,a)') 'Write ICON grid connectivity file: ', TRIM(filename_grfinfo)
      CALL message ('', message_text)

      CALL nf(nf_create(TRIM(filename_grfinfo), nf_clobber, ncid))
      CALL nf(nf_set_fill(ncid, nf_nofill, old_mode))

      ! Set global attributes
      !
      CALL set_global_attributes(.TRUE.)

      ! Set dimensions
      !
      CALL set_dimensions()

      ! Define variables for parent-child connectivity information
      !
      CALL define_connectivity_variables()

      CALL nf(nf_enddef(ncid))

      ! Write variables for parent-child connectivity information
      !
      CALL write_connectivity_variables()

      ! close file
      CALL nf(nf_close(ncid))

    ENDIF

    IF (current_number_of_grid_used /= 0) THEN
      CALL dump_grid_table_xml_metadata(TRIM(filename),current_number_of_grid_used)
    ENDIF

   !------------------------------------------------------------------------
    write(*,*) '---',TRIM(filename), '---'
    str_idx=LBOUND(p%verts%start_idx, 1)
    end_idx=str_idx+SIZE(p%verts%start_idx, 1)-1
    DO i=str_idx,end_idx
      write(*,*) 'verts%start_idx, end:', i, p%verts%start_idx(i,1), p%verts%end_idx(i,1)
    ENDDO

    str_idx=LBOUND(p%edges%start_idx, 1)
    end_idx=str_idx+SIZE(p%edges%start_idx, 1)-1
    DO i=str_idx,end_idx
      write(*,*) 'edges%start_idx, end:', i, p%edges%start_idx(i,1), p%edges%end_idx(i,1)
    ENDDO

    str_idx=LBOUND(p%cells%start_idx, 1)
    end_idx=str_idx+SIZE(p%cells%start_idx, 1)-1
    DO i=str_idx,end_idx
      write(*,*) 'cells%start_idx, end:', i, p%cells%start_idx(i,1), p%cells%end_idx(i,1)
    ENDDO
    write(*,*) '-------------------'
    
  CONTAINS

  SUBROUTINE set_global_attributes(lwrite_pcc)

    LOGICAL, INTENT(IN) :: lwrite_pcc
    INTEGER :: i

    CALL nf(nf_put_att_text    (ncid, nf_global, 'title', 21, 'ICON grid description'))
    CALL nf(nf_put_att_text    (ncid, nf_global, 'history', command_line_len, command_line))
    CALL nf(nf_put_att_text    (ncid, nf_global, 'institution', 59, &
      & 'Max Planck Institute for Meteorology/Deutscher Wetterdienst'))
    CALL nf(nf_put_att_text    (ncid, nf_global, 'source', 10, 'icon-dev'))
    CALL nf(nf_put_att_text    (ncid, nf_global, 'uuidOfHGrid' , uuid_string_length, TRIM(uuid_grid)))

    IF (lwrite_pcc) THEN
      CALL nf(nf_put_att_text    (ncid, nf_global, 'uuidOfParHGrid' , uuid_string_length, TRIM(uuid_parent)))
      DO i=1,5
        WRITE(child_id,'(i1)') i
        CALL nf(nf_put_att_text    (ncid, nf_global, 'uuidOfChiHGrid_'//child_id , uuid_string_length, TRIM(uuid_child(i))))
      ENDDO
    ENDIF

    IF (current_number_of_grid_used == 0) THEN
      CALL message('','number_of_grid_used is 0 and cannot be added to the ICON master grid table')
      CALL nf(nf_put_att_int   (ncid, nf_global, 'number_of_grid_used', nf_int, 1, 0))
      CALL nf(nf_put_att_text  (ncid, nf_global, 'ICON_grid_file_uri' , 7, 'private'))
    ELSE
      CALL nf(nf_put_att_int   (ncid, nf_global, 'number_of_grid_used', nf_int, 1, current_number_of_grid_used))
      CALL nf(nf_put_att_text  (ncid, nf_global, 'ICON_grid_file_uri' , &
           &                                      LEN_TRIM(uri_pathname//TRIM(filename)), TRIM(uri_pathname//TRIM(filename))))
    ENDIF
    CALL nf(nf_put_att_int     (ncid, nf_global, 'centre', nf_int, 1, centre))
    CALL nf(nf_put_att_int     (ncid, nf_global, 'subcentre', nf_int, 1, subcentre))
    CALL nf(nf_put_att_int     (ncid, nf_global, 'outname_style', nf_int, 1, outname_style))
    CALL nf(nf_put_att_text    (ncid, nf_global, 'grid_mapping_name' , 18, 'lat_long_on_sphere'))
    CALL nf(nf_put_att_text    (ncid, nf_global, 'crs_id' , 28, 'urn:ogc:def:cs:EPSG:6.0:6422'))
    CALL nf(nf_put_att_text    (ncid, nf_global, 'crs_name', 30,'Spherical 2D Coordinate System'))
    CALL nf(nf_put_att_text    (ncid, nf_global, 'ellipsoid_name' , 6, 'Sphere'))
    CALL nf(nf_put_att_double  (ncid, nf_global, 'semi_major_axis' , nf_double, 1, earth_radius))
    CALL nf(nf_put_att_double  (ncid, nf_global, 'inverse_flattening' , nf_double, 1, 0.0_wp))
    CALL nf(nf_put_att_int     (ncid, nf_global, 'grid_level', nf_int, 1, i_lev))
    CALL nf(nf_put_att_int     (ncid, nf_global, 'grid_root', nf_int, 1, grid_root))
    CALL nf(nf_put_att_int     (ncid, nf_global, 'grid_ID', nf_int, 1, dom_id))
    CALL nf(nf_put_att_int     (ncid, nf_global, 'parent_grid_ID', nf_int, 1, parent_dom_id))
    CALL nf(nf_put_att_int     (ncid, nf_global, 'max_childdom', nf_int, 1, max_childdom))

  END SUBROUTINE set_global_attributes


  SUBROUTINE set_dimensions

    CALL nf(nf_def_dim(ncid, 'cell',   i_nc, dim_ncell))
    CALL nf(nf_def_dim(ncid, 'vertex', i_nv, dim_nvertex))
    CALL nf(nf_def_dim(ncid, 'edge',   i_ne, dim_nedge))
    !
    CALL nf(nf_def_dim(ncid, 'nc',    2, dim_ncells_per_edge))
    CALL nf(nf_def_dim(ncid, 'nv',    3, dim_nvertex_per_cell))
    CALL nf(nf_def_dim(ncid, 'ne',    6, dim_nedges_per_vertex))
    !
    CALL nf(nf_def_dim(ncid, 'no',    4, dim_nchilds_per_cell))

    ! Dimensions for refinement
    CALL nf(nf_def_dim(ncid, 'two_grf',        2, dim_two))
    CALL nf(nf_def_dim(ncid, 'max_chdom', 1, dim_nchdom)) ! was originally max_childdom; kept for backward 
    dim_list = max_rlcell-min_rlcell+1                    ! compatibility of old grid files
    CALL nf(nf_def_dim(ncid, 'cell_grf',dim_list, dim_cell_refine))
    dim_list = max_rledge-min_rledge+1
    CALL nf(nf_def_dim(ncid, 'edge_grf',dim_list, dim_edge_refine))
    dim_list = max_rlvert-min_rlvert+1
    CALL nf(nf_def_dim(ncid, 'vert_grf',dim_list, dim_vert_refine))

  END SUBROUTINE set_dimensions


  SUBROUTINE define_standard_variables
    ! cell part:
    !
    CALL nf(nf_def_var(ncid, 'clon', nf_double, 1, dim_ncell, varid_clon))
    CALL nf(nf_put_att_text(ncid, varid_clon, 'long_name', 16, 'center longitude'))
    CALL nf(nf_put_att_text(ncid, varid_clon, 'units', 6, 'radian'))
    CALL nf(nf_put_att_text(ncid, varid_clon, 'standard_name', 14, 'grid_longitude'))
    CALL nf(nf_put_att_text(ncid, varid_clon, 'bounds', 13, 'clon_vertices'))
    !
    CALL nf(nf_def_var(ncid, 'clat', nf_double, 1, dim_ncell, varid_clat))
    CALL nf(nf_put_att_text(ncid, varid_clat, 'long_name', 15, 'center latitude'))
    CALL nf(nf_put_att_text(ncid, varid_clat, 'units', 6, 'radian'))
    CALL nf(nf_put_att_text(ncid, varid_clat, 'standard_name', 13, 'grid_latitude'))
    CALL nf(nf_put_att_text(ncid, varid_clat, 'bounds', 13, 'clat_vertices'))
    !
    dimids = (/ dim_nvertex_per_cell, dim_ncell /)
    CALL nf(nf_def_var(ncid, 'clon_vertices', nf_double, 2, dimids, varid_clonv))
    CALL nf(nf_put_att_text(ncid, varid_clonv, 'units', 6, 'radian'))
    !
    dimids = (/ dim_nvertex_per_cell, dim_ncell /)
    CALL nf(nf_def_var(ncid, 'clat_vertices', nf_double, 2, dimids, varid_clatv))
    CALL nf(nf_put_att_text(ncid, varid_clatv, 'units', 6, 'radian'))
    !
    ! vertex part:
    !
    CALL nf(nf_def_var(ncid, 'vlon', nf_double, 1, dim_nvertex, varid_vlon))
    CALL nf(nf_put_att_text(ncid, varid_vlon, 'long_name', 16, 'vertex longitude'))
    CALL nf(nf_put_att_text(ncid, varid_vlon, 'units', 6, 'radian'))
    CALL nf(nf_put_att_text(ncid, varid_vlon, 'standard_name', 14, 'grid_longitude'))
    CALL nf(nf_put_att_text(ncid, varid_vlon, 'bounds', 13, 'vlon_vertices'))
    !
    CALL nf(nf_def_var(ncid, 'vlat', nf_double, 1, dim_nvertex, varid_vlat))
    CALL nf(nf_put_att_text(ncid, varid_vlat, 'long_name', 15, 'vertex latitude'))
    CALL nf(nf_put_att_text(ncid, varid_vlat, 'units', 6, 'radian'))
    CALL nf(nf_put_att_text(ncid, varid_vlat, 'standard_name', 13, 'grid_latitude'))
    CALL nf(nf_put_att_text(ncid, varid_vlat, 'bounds', 13, 'vlat_vertices'))
    !
    dimids = (/ dim_nedges_per_vertex, dim_nvertex /)
    CALL nf(nf_def_var(ncid, 'vlon_vertices', nf_double, 2, dimids, varid_vlonv))
    CALL nf(nf_put_att_text(ncid, varid_vlonv, 'units', 6, 'radian'))
    !
    dimids = (/ dim_nedges_per_vertex, dim_nvertex /)
    CALL nf(nf_def_var(ncid, 'vlat_vertices', nf_double, 2, dimids, varid_vlatv))
    CALL nf(nf_put_att_text(ncid, varid_vlatv, 'units', 6, 'radian'))
    !
    ! edge part:
    !
    CALL nf(nf_def_var(ncid, 'elon', nf_double, 1, dim_nedge, varid_elon))
    CALL nf(nf_put_att_text(ncid, varid_elon, 'long_name', 23, 'edge midpoint longitude'))
    CALL nf(nf_put_att_text(ncid, varid_elon, 'units', 6, 'radian'))
    CALL nf(nf_put_att_text(ncid, varid_elon, 'standard_name', 14, 'grid_longitude'))
    CALL nf(nf_put_att_text(ncid, varid_elon, 'bounds', 13, 'elon_vertices'))
    !
    CALL nf(nf_def_var(ncid, 'elat', nf_double, 1, dim_nedge, varid_elat))
    CALL nf(nf_put_att_text(ncid, varid_elat, 'long_name', 22, 'edge midpoint latitude') )
    CALL nf(nf_put_att_text(ncid, varid_elat, 'units', 6, 'radian'))
    CALL nf(nf_put_att_text(ncid, varid_elat, 'standard_name', 13, 'grid_latitude'))
    CALL nf(nf_put_att_text(ncid, varid_elat, 'bounds', 13, 'elat_vertices'))
    !
    dimids = (/ dim_nchilds_per_cell, dim_nedge /)
    CALL nf(nf_def_var(ncid, 'elon_vertices', nf_double, 2, dimids, varid_elonv))
    CALL nf(nf_put_att_text(ncid, varid_elonv, 'units', 6, 'radian'))
    !
    dimids = (/ dim_nchilds_per_cell, dim_nedge /)
    CALL nf(nf_def_var(ncid, 'elat_vertices', nf_double, 2, dimids, varid_elatv))
    CALL nf(nf_put_att_text(ncid, varid_elatv, 'units', 6, 'radian'))
    !
    !---------------------------------------------------------------------
    !
    ! ifs2icon variables 
    !
    CALL nf(nf_def_var(ncid, 'ifs2icon_cell_grid', nf_double, 1, dim_ncell, ifs2icon_cell))
    CALL nf(nf_put_att_text(ncid, ifs2icon_cell, 'long_name', 17, 'ifs to icon cells'))
    CALL nf(nf_put_att_text(ncid, ifs2icon_cell, 'coordinates', 9, 'clon clat'))
    !
    CALL nf(nf_def_var(ncid, 'ifs2icon_edge_grid', nf_double, 1, dim_nedge, ifs2icon_edge))
    CALL nf(nf_put_att_text(ncid, ifs2icon_edge, 'long_name', 16, 'ifs to icon edge'))
    CALL nf(nf_put_att_text(ncid, ifs2icon_edge, 'coordinates', 9, 'elon elat'))
    !
    CALL nf(nf_def_var(ncid, 'ifs2icon_vertex_grid', nf_double, 1, dim_nvertex, ifs2icon_vertex))
    CALL nf(nf_put_att_text(ncid, ifs2icon_vertex, 'long_name', 18, 'ifs to icon vertex'))
    CALL nf(nf_put_att_text(ncid, ifs2icon_vertex, 'coordinates', 9, 'vlon vlat'))
    !
    ! test variables (areas)
    !
    CALL nf(nf_def_var(ncid, 'cell_area', nf_double, 1, dim_ncell, varid_carea))
    CALL nf(nf_put_att_text(ncid, varid_carea, 'long_name', 17, 'area of grid cell'))
    CALL nf(nf_put_att_text(ncid, varid_carea, 'units', 9, 'steradian'))
    CALL nf(nf_put_att_text(ncid, varid_carea, 'standard_name', 4, 'area'))
    CALL nf(nf_put_att_text(ncid, varid_carea, 'coordinates', 9, 'clon clat'))
    !
    CALL nf(nf_def_var(ncid, 'dual_area', nf_double, 1, dim_nvertex, varid_varea))
    CALL nf(nf_put_att_text(ncid, varid_varea, 'long_name', 40, &
      & 'areas of dual hexagonal/pentagonal cells'))
    CALL nf(nf_put_att_text(ncid, varid_varea, 'units', 9, 'steradian'))
    CALL nf(nf_put_att_text(ncid, varid_varea, 'standard_name', 4, 'area'))
    CALL nf(nf_put_att_text(ncid, varid_varea, 'coordinates', 9, 'vlon vlat'))
    !
    CALL nf(nf_def_var(ncid, 'quadrilateral_area', nf_double, 1, dim_nedge, varid_earea))
    CALL nf(nf_put_att_text(ncid, varid_earea, 'long_name', 69, &
      & 'edge area given by the quadrilateral of adjacent centers and vertices'))
    CALL nf(nf_put_att_text(ncid, varid_earea, 'units', 9, 'steradian'))
    CALL nf(nf_put_att_text(ncid, varid_earea, 'standard_name', 4, 'area'))
    CALL nf(nf_put_att_text(ncid, varid_earea, 'coordinates', 9, 'elon elat'))
    !
    !---------------------------------------------------------------------
    !
    ! private grid information
    !
    CALL nf(nf_def_var(ncid, 'lon_cell_centre', nf_double, 1, dim_ncell, varid1))
    CALL nf(nf_put_att_text(ncid, varid1, 'long_name', 24, 'longitude of cell centre'))
    CALL nf(nf_put_att_text(ncid, varid1, 'units', 6, 'radian'))
    CALL nf(nf_put_att_text(ncid, varid1, 'cdi', 6, 'ignore'))
    !
    CALL nf(nf_def_var(ncid, 'lat_cell_centre', nf_double, 1, dim_ncell, varid2))
    CALL nf(nf_put_att_text(ncid, varid2, 'long_name', 23, 'latitude of cell centre'))
    CALL nf(nf_put_att_text(ncid, varid2, 'units', 6, 'radian'))
    CALL nf(nf_put_att_text(ncid, varid2, 'cdi', 6, 'ignore'))
    !
    CALL nf(nf_def_var(ncid, 'longitude_vertices', nf_double, 1, dim_nvertex, varid3))
    CALL nf(nf_put_att_text(ncid, varid3, 'long_name', 21, 'longitude of vertices'))
    CALL nf(nf_put_att_text(ncid, varid3, 'units', 6, 'radian'))
    CALL nf(nf_put_att_text(ncid, varid3, 'cdi', 6, 'ignore'))
    !
    CALL nf(nf_def_var(ncid, 'latitude_vertices', nf_double, 1, dim_nvertex, varid4))
    CALL nf(nf_put_att_text(ncid, varid4, 'long_name', 20, 'latitude of vertices'))
    CALL nf(nf_put_att_text(ncid, varid4, 'units', 6, 'radian'))
    CALL nf(nf_put_att_text(ncid, varid4, 'cdi', 6, 'ignore'))
    !
    CALL nf(nf_def_var(ncid, 'lon_edge_centre', nf_double, 1, dim_nedge, varid5))
    CALL nf(nf_put_att_text(ncid, varid5, 'long_name', 28, 'longitudes of edge midpoints'))
    CALL nf(nf_put_att_text(ncid, varid5, 'units', 6, 'radian'))
    CALL nf(nf_put_att_text(ncid, varid5, 'cdi', 6, 'ignore'))
    !
    CALL nf(nf_def_var(ncid, 'lat_edge_centre', nf_double, 1, dim_nedge, varid6))
    CALL nf(nf_put_att_text(ncid, varid6, 'long_name', 27, 'latitudes of edge midpoints'))
    CALL nf(nf_put_att_text(ncid, varid6, 'units', 6, 'radian'))
    CALL nf(nf_put_att_text(ncid, varid6, 'cdi', 6, 'ignore'))
    !
    dimids = (/ dim_ncell, dim_nvertex_per_cell /)
    CALL nf(nf_def_var(ncid, 'edge_of_cell', nf_int, 2, dimids, varid7))
    CALL nf(nf_put_att_text(ncid, varid7, 'long_name', 29, 'edges of each triangular cell'))
    CALL nf(nf_put_att_text(ncid, varid7, 'cdi', 6, 'ignore'))
    !
    dimids = (/ dim_ncell, dim_nvertex_per_cell /)
    CALL nf(nf_def_var(ncid, 'vertex_of_cell', nf_int, 2, dimids, varid8))
    CALL nf(nf_put_att_text(ncid, varid8, 'long_name', 32, 'vertices of each triangular cell'))
    CALL nf(nf_put_att_text(ncid, varid8, 'cdi', 6, 'ignore'))
    !
    dimids = (/ dim_nedge, dim_ncells_per_edge /)
    CALL nf(nf_def_var(ncid, 'adjacent_cell_of_edge', nf_int, 2, dimids, varid9))
    CALL nf(nf_put_att_text(ncid, varid9, 'long_name', 27, 'cells adjacent to each edge'))
    CALL nf(nf_put_att_text(ncid, varid9, 'cdi', 6, 'ignore'))
    !
    dimids = (/ dim_nedge, dim_ncells_per_edge /)
    CALL nf(nf_def_var(ncid, 'edge_vertices', nf_int, 2, dimids, varid10))
    CALL nf(nf_put_att_text(ncid, varid10, 'long_name', 35,'vertices at the end of of each edge'))
    CALL nf(nf_put_att_text(ncid, varid10, 'cdi', 6, 'ignore'))
    !
    dimids = (/ dim_nvertex, dim_nedges_per_vertex /)
    CALL nf(nf_def_var(ncid, 'cells_of_vertex', nf_int, 2, dimids, varid11))
    CALL nf(nf_put_att_text(ncid, varid11, 'long_name', 24, 'cells around each vertex'))
    CALL nf(nf_put_att_text(ncid, varid11, 'cdi', 6, 'ignore'))
    !
    dimids = (/ dim_nvertex, dim_nedges_per_vertex /)
    CALL nf(nf_def_var(ncid, 'edges_of_vertex', nf_int, 2, dimids, varid12))
    CALL nf(nf_put_att_text(ncid, varid12, 'long_name', 24, 'edges around each vertex'))
    CALL nf(nf_put_att_text(ncid, varid12, 'cdi', 6, 'ignore'))
    !
    dimids = (/ dim_nvertex, dim_nedges_per_vertex /)
    CALL nf(nf_def_var(ncid, 'vertices_of_vertex', nf_int, 2, dimids, varid13))
    CALL nf(nf_put_att_text(ncid, varid13, 'long_name', 27, 'vertices around each vertex'))
    CALL nf(nf_put_att_text(ncid, varid13, 'cdi', 6, 'ignore'))
    !
    CALL nf(nf_def_var(ncid, 'cell_area_p', nf_double, 1, dim_ncell, varid14))
    CALL nf(nf_put_att_text(ncid, varid14, 'long_name', 17, 'area of grid cell'))
    CALL nf(nf_put_att_text(ncid, varid14, 'units', 9, 'steradian'))
    CALL nf(nf_put_att_text(ncid, varid14, 'cdi', 6, 'ignore'))
    !
    CALL nf(nf_def_var(ncid, 'dual_area_p', nf_double, 1, dim_nvertex, varid15))
    CALL nf(nf_put_att_text(ncid, varid15, 'long_name', 40, &
      & 'areas of dual hexagonal/pentagonal cells'))
    CALL nf(nf_put_att_text(ncid, varid15, 'units', 9, 'steradian'))
    CALL nf(nf_put_att_text(ncid, varid15, 'cdi', 6, 'ignore'))
    !
    CALL nf(nf_def_var(ncid, 'edge_length', nf_double, 1, dim_nedge, varid16))
    CALL nf(nf_put_att_text(ncid, varid16, 'long_name',36,'lengths of edges of triangular cells'))
    CALL nf(nf_put_att_text(ncid, varid16, 'units', 6, 'radian'))
    CALL nf(nf_put_att_text(ncid, varid16, 'cdi', 6, 'ignore'))
    !
    dimids = (/ dim_nedge, dim_ncells_per_edge /)
    CALL nf(nf_def_var(ncid, 'edge_cell_distance', nf_double, 2, dimids, varid40))
    CALL nf(nf_put_att_text(ncid, varid40, 'long_name', 63, &
      & 'distances between edge midpoint and adjacent triangle midpoints'))
    CALL nf(nf_put_att_text(ncid, varid40, 'units', 1, 'm'))
    CALL nf(nf_put_att_text(ncid, varid40, 'cdi', 6, 'ignore'))
    !
    CALL nf(nf_def_var(ncid, 'dual_edge_length', nf_double, 1, dim_nedge, varid17))
    CALL nf(nf_put_att_text(ncid, varid17, 'long_name', 70, &
      & 'lengths of dual edges (distances between triangular cell circumcenters)'))
    CALL nf(nf_put_att_text(ncid, varid17, 'units', 6, 'radian'))
    CALL nf(nf_put_att_text(ncid, varid17, 'cdi', 6, 'ignore'))
    !
    dimids = (/ dim_nedge, dim_ncells_per_edge /)
    CALL nf(nf_def_var(ncid, 'edge_vert_distance', nf_double, 2, dimids, varid18))
    CALL nf(nf_put_att_text(ncid, varid18, 'long_name', 57, &
      & 'distances between edge midpoint and vertices of that edge'))
    CALL nf(nf_put_att_text(ncid, varid18, 'units', 1, 'm'))
    CALL nf(nf_put_att_text(ncid, varid18, 'cdi', 6, 'ignore'))
    !
    CALL nf(nf_def_var(ncid, 'zonal_normal_primal_edge', nf_double, 1, dim_nedge, varid19))
    CALL nf(nf_put_att_text(ncid, varid19, 'long_name', 40, &
      & 'zonal component of normal to primal edge'))
    CALL nf(nf_put_att_text(ncid, varid19, 'units', 6, 'radian'))
    CALL nf(nf_put_att_text(ncid, varid19, 'cdi', 6, 'ignore'))
    !
    CALL nf(nf_def_var(ncid, 'meridional_normal_primal_edge', nf_double, 1, dim_nedge, varid20))
    CALL nf(nf_put_att_text(ncid, varid20, 'long_name', 45, &
      & 'meridional component of normal to primal edge'))
    CALL nf(nf_put_att_text(ncid, varid20, 'units', 6, 'radian'))
    CALL nf(nf_put_att_text(ncid, varid20, 'cdi', 6, 'ignore'))
    !
    CALL nf(nf_def_var(ncid, 'zonal_normal_dual_edge', nf_double, 1, dim_nedge, varid21))
    CALL nf(nf_put_att_text(ncid, varid21, 'long_name', 38, &
      & 'zonal component of normal to dual edge'))
    CALL nf(nf_put_att_text(ncid, varid21, 'units', 6, 'radian'))
    CALL nf(nf_put_att_text(ncid, varid21, 'cdi', 6, 'ignore'))
    !
    CALL nf(nf_def_var(ncid, 'meridional_normal_dual_edge', nf_double, 1, dim_nedge, varid22))
    CALL nf(nf_put_att_text(ncid, varid22, 'long_name', 43, &
      & 'meridional component of normal to dual edge'))
    CALL nf(nf_put_att_text(ncid, varid22, 'units', 6, 'radian'))
    CALL nf(nf_put_att_text(ncid, varid22, 'cdi', 6, 'ignore'))
    !
    dimids = (/ dim_ncell, dim_nvertex_per_cell /)
    CALL nf(nf_def_var(ncid, 'orientation_of_normal', nf_int, 2, dimids, varid23))
    CALL nf(nf_put_att_text(ncid, varid23, 'long_name', 48, &
      & 'orientations of normals to triangular cell edges'))
    CALL nf(nf_put_att_text(ncid, varid23, 'cdi', 6, 'ignore'))
    !
    CALL nf(nf_def_var(ncid, 'cell_index', nf_int, 1, dim_ncell, varid24))
    CALL nf(nf_put_att_text(ncid, varid24, 'long_name', 10, 'cell index'))
    CALL nf(nf_put_att_text(ncid, varid24, 'cdi', 6, 'ignore'))
    !
    dimids = (/ dim_ncell, dim_nvertex_per_cell /)
    CALL nf(nf_def_var(ncid, 'neighbor_cell_index', nf_int, 2, dimids, varid26))
    CALL nf(nf_put_att_text(ncid, varid26, 'long_name', 19, 'cell neighbor index'))
    CALL nf(nf_put_att_text(ncid, varid26, 'cdi', 6, 'ignore'))
    !
    CALL nf(nf_def_var(ncid, 'edge_index', nf_int, 1, dim_nedge, varid28))
    CALL nf(nf_put_att_text(ncid, varid28, 'long_name', 10, 'edge index'))
    CALL nf(nf_put_att_text(ncid, varid28, 'cdi', 6, 'ignore'))
    !
    CALL nf(nf_def_var(ncid, 'vertex_index', nf_int, 1, dim_nvertex, varid29))
    CALL nf(nf_put_att_text(ncid, varid29, 'long_name', 14, 'vertices index'))
    CALL nf(nf_put_att_text(ncid, varid29, 'cdi', 6, 'ignore'))
    !
    dimids = (/ dim_nvertex, dim_nedges_per_vertex /)
    CALL nf(nf_def_var(ncid, 'edge_orientation', nf_int, 2, dimids, varid30))
    CALL nf(nf_put_att_text(ncid, varid30, 'long_name', 16, 'edge orientation'))
    CALL nf(nf_put_att_text(ncid, varid30, 'cdi', 6, 'ignore'))
    !
    CALL nf(nf_def_var(ncid, 'edge_system_orientation', nf_int, 1, dim_nedge, varid31))
    CALL nf(nf_put_att_text(ncid, varid31, 'long_name', 23, 'edge system orientation'))
    CALL nf(nf_put_att_text(ncid, varid31, 'cdi', 6, 'ignore'))

  END SUBROUTINE define_standard_variables


  SUBROUTINE define_connectivity_variables
    !
    !
    CALL nf(nf_def_var(ncid, 'parent_cell_index', nf_int, 1, dim_ncell, varid25))
    CALL nf(nf_put_att_text(ncid, varid25, 'long_name', 17, 'parent cell index'))
    CALL nf(nf_put_att_text(ncid, varid25, 'cdi', 6, 'ignore'))
    !
    dimids = (/ dim_ncell, dim_nchilds_per_cell /)
    CALL nf(nf_def_var(ncid, 'child_cell_index', nf_int, 2, dimids, varid27))
    CALL nf(nf_put_att_text(ncid, varid27, 'long_name', 16, 'child cell index'))
    CALL nf(nf_put_att_text(ncid, varid27, 'cdi', 6, 'ignore'))
    !
    CALL nf(nf_def_var(ncid, 'child_cell_id', nf_int, 1, dim_ncell, varid41))
    CALL nf(nf_put_att_text(ncid, varid41, 'long_name', 23, 'domain ID of child cell'))
    CALL nf(nf_put_att_text(ncid, varid41, 'cdi', 6, 'ignore'))
    !
    CALL nf(nf_def_var(ncid, 'refin_c_ctrl', nf_int, 1, dim_ncell, varid32))
    CALL nf(nf_put_att_text(ncid, varid32, 'long_name', 33, 'refinement control flag for cells'))
    CALL nf(nf_put_att_text(ncid, varid32, 'cdi', 6, 'ignore'))
    !
    dimids = (/ dim_cell_refine, dim_two /)
    CALL nf(nf_def_var(ncid, 'index_c_list', nf_int, 2, dimids, varid33))
    CALL nf(nf_put_att_text(ncid, varid33, 'long_name', 73, &
      & 'list of start and end indices for each refinement control level for cells'))
    CALL nf(nf_put_att_text(ncid, varid33, 'cdi', 6, 'ignore'))
    !
    dimids = (/ dim_cell_refine, dim_nchdom /)
    CALL nf(nf_def_var(ncid, 'start_idx_c', nf_int, 2, dimids, varid43))
    CALL nf(nf_put_att_text(ncid, varid43, 'long_name', 65, &
      & 'list of start indices for each refinement control level for cells'))
    CALL nf(nf_put_att_text(ncid, varid43, 'cdi', 6, 'ignore'))
    !
    CALL nf(nf_def_var(ncid, 'end_idx_c', nf_int, 2, dimids, varid44))
    CALL nf(nf_put_att_text(ncid, varid44, 'long_name', 63, &
      & 'list of end indices for each refinement control level for cells'))
    CALL nf(nf_put_att_text(ncid, varid44, 'cdi', 6, 'ignore'))
    !
    CALL nf(nf_def_var(ncid, 'refin_e_ctrl', nf_int, 1, dim_nedge, varid34))
    CALL nf(nf_put_att_text(ncid, varid34, 'long_name', 33, 'refinement control flag for edges'))
    CALL nf(nf_put_att_text(ncid, varid34, 'cdi', 6, 'ignore'))
    !
    dimids = (/ dim_edge_refine, dim_two /)
    CALL nf(nf_def_var(ncid, 'index_e_list', nf_int, 2, dimids, varid35))
    CALL nf(nf_put_att_text(ncid, varid35, 'long_name', 73, &
      & 'list of start and end indices for each refinement control level for edges'))
    CALL nf(nf_put_att_text(ncid, varid35, 'cdi', 6, 'ignore'))
    !
    dimids = (/ dim_edge_refine, dim_nchdom /)
    CALL nf(nf_def_var(ncid, 'start_idx_e', nf_int, 2, dimids, varid45))
    CALL nf(nf_put_att_text(ncid, varid45, 'long_name', 65, &
      & 'list of start indices for each refinement control level for edges'))
    CALL nf(nf_put_att_text(ncid, varid45, 'cdi', 6, 'ignore'))
    !
    CALL nf(nf_def_var(ncid, 'end_idx_e', nf_int, 2, dimids, varid46))
    CALL nf(nf_put_att_text(ncid, varid46, 'long_name', 63, &
      & 'list of end indices for each refinement control level for edges'))
    CALL nf(nf_put_att_text(ncid, varid46, 'cdi', 6, 'ignore'))
    !
    CALL nf(nf_def_var(ncid, 'refin_v_ctrl', nf_int, 1, dim_nvertex, varid36))
    CALL nf(nf_put_att_text(ncid, varid36, 'long_name',36,'refinement control flag for vertices'))
    CALL nf(nf_put_att_text(ncid, varid36, 'cdi', 6, 'ignore'))
    !
    dimids = (/ dim_vert_refine, dim_two /)
    CALL nf(nf_def_var(ncid, 'index_v_list', nf_int, 2, dimids, varid37))
    CALL nf(nf_put_att_text(ncid, varid37, 'long_name', 76, &
      & 'list of start and end indices for each refinement control level for vertices'))
    CALL nf(nf_put_att_text(ncid, varid37, 'cdi', 6, 'ignore'))
    !
    dimids = (/ dim_vert_refine, dim_nchdom /)
    CALL nf(nf_def_var(ncid, 'start_idx_v', nf_int, 2, dimids, varid47))
    CALL nf(nf_put_att_text(ncid, varid47, 'long_name', 68, &
      & 'list of start indices for each refinement control level for vertices'))
    CALL nf(nf_put_att_text(ncid, varid47, 'cdi', 6, 'ignore'))
    !
    CALL nf(nf_def_var(ncid, 'end_idx_v', nf_int, 2, dimids, varid48))
    CALL nf(nf_put_att_text(ncid, varid48, 'long_name', 66, &
      & 'list of end indices for each refinement control level for vertices'))
    CALL nf(nf_put_att_text(ncid, varid48, 'cdi', 6, 'ignore'))
    !
    CALL nf(nf_def_var(ncid, 'parent_edge_index', nf_int, 1, dim_nedge, varid38))
    CALL nf(nf_put_att_text(ncid, varid38, 'long_name', 17, 'parent edge index'))
    CALL nf(nf_put_att_text(ncid, varid38, 'cdi', 6, 'ignore'))
    !
    dimids = (/ dim_nedge, dim_nchilds_per_cell /)
    CALL nf(nf_def_var(ncid, 'child_edge_index', nf_int, 2, dimids, varid39))
    CALL nf(nf_put_att_text(ncid, varid39, 'long_name', 16, 'child edge index'))
    CALL nf(nf_put_att_text(ncid, varid39, 'cdi', 6, 'ignore'))
    !
    CALL nf(nf_def_var(ncid, 'child_edge_id', nf_int, 1, dim_nedge, varid42))
    CALL nf(nf_put_att_text(ncid, varid42, 'long_name', 23, 'domain ID of child edge'))
    CALL nf(nf_put_att_text(ncid, varid42, 'cdi', 6, 'ignore'))
    !
    CALL nf(nf_def_var(ncid, 'phys_cell_id', nf_int, 1, dim_ncell, varid49))
    CALL nf(nf_put_att_text(ncid, varid49, 'long_name', 26, 'physical domain ID of cell'))
    CALL nf(nf_put_att_text(ncid, varid49, 'cdi', 6, 'ignore'))
    !
    CALL nf(nf_def_var(ncid, 'phys_edge_id', nf_int, 1, dim_nedge, varid50))
    CALL nf(nf_put_att_text(ncid, varid50, 'long_name', 26, 'physical domain ID of edge'))
    CALL nf(nf_put_att_text(ncid, varid50, 'cdi', 6, 'ignore'))
    !
  END SUBROUTINE define_connectivity_variables



  SUBROUTINE write_standard_variables
    !
    !--------------------------------------------------------------------------------------
    !
    ! cell part:
    !
    ! Transpose of index array necessary for CF-1.1 Convention
    !
    ALLOCATE(zv2dx(3, i_nc), zv2dy(3, i_nc))
    DO j = 1, 3
      DO i = 1, i_nc
        zv2dx(j,i) = p%verts%vertex(p%cells%vertex_index(i,j))%lon
      ENDDO
    ENDDO
    DO j = 1, 3
      DO i = 1, i_nc
        zv2dy(j,i) = p%verts%vertex(p%cells%vertex_index(i,j))%lat
      ENDDO
    ENDDO
    WHERE (ABS(zv2dx(:,:)) < EPSILON(0.0_wp))
    zv2dx(:,:) = 0.0_wp
    endwhere
    WHERE (ABS(zv2dy(:,:)) < EPSILON(0.0_wp))
    zv2dy(:,:) = 0.0_wp
    endwhere
    DO j = 1, i_nc
      DO i = 1, 3
        IF (ABS(zv2dy(i,j)) > 0.5_wp*pi-EPSILON(0.0_wp)) THEN
          zv2dx(i,j) = p%cells%center(j)%lon
        ENDIF
      ENDDO
    ENDDO
    CALL nf(nf_put_var_double(ncid, varid_clatv, zv2dy))
    CALL nf(nf_put_var_double(ncid, varid_clonv, zv2dx))
    DEALLOCATE(zv2dx, zv2dy)

    CALL nf(nf_put_var_double(ncid, varid_clon, p%cells%center(:)%lon))
    CALL nf(nf_put_var_double(ncid, varid_clat, p%cells%center(:)%lat))

    !
    ! vertex part:
    !
    CALL nf(nf_put_var_double(ncid, varid_vlon,  p%verts%vertex(:)%lon))
    CALL nf(nf_put_var_double(ncid, varid_vlat,  p%verts%vertex(:)%lat))
    !
    ! Transpose of index array necessary for CF-1.1 Convention
    !
    ALLOCATE(zv2d(6, i_nv))
    DO j = 1, 6
      DO i = 1, i_nv
        IF ((p%verts%cell_index(i,j) == 0) .and. (p%verts%refin_ctrl(i) /= 1)) THEN
          zv2d(7-j,i) = p%cells%center(p%verts%cell_index(i,5))%lon
        ELSE IF ((p%verts%cell_index(i,j) < 0) .or. (p%verts%refin_ctrl(i) == 1)) THEN
          zv2d(7-j,i) = 0._wp
        ELSE
          zv2d(7-j,i) = p%cells%center(p%verts%cell_index(i,j))%lon
        ENDIF
      ENDDO
    ENDDO
    CALL nf(nf_put_var_double(ncid, varid_vlonv, zv2d))
    DO j = 1, 6
      DO i = 1, i_nv
        IF ((p%verts%cell_index(i,j) == 0) .and. (p%verts%refin_ctrl(i) /= 1)) THEN
          zv2d(7-j,i) = p%cells%center(p%verts%cell_index(i,5))%lat
        ELSE IF ((p%verts%cell_index(i,j) < 0) .or. (p%verts%refin_ctrl(i) == 1)) THEN
          zv2d(7-j,i) = 0._wp
        ELSE
          zv2d(7-j,i) = p%cells%center(p%verts%cell_index(i,j))%lat
        ENDIF
      ENDDO
    ENDDO
    CALL nf(nf_put_var_double(ncid, varid_vlatv, zv2d))
    DEALLOCATE (zv2d)
    !
    ! edge part:
    !
    CALL nf(nf_put_var_double(ncid, varid_elon,  p%edges%center(:)%lon))
    CALL nf(nf_put_var_double(ncid, varid_elat,  p%edges%center(:)%lat))
    !
    ! Transpose of index array necessary for CF-1.1 Convention
    !
    ALLOCATE(zv2dx(4, i_ne), zv2dy(4, i_ne))
    DO i = 1, i_ne
      zv2dx(1,i) = p%verts%vertex(p%edges%vertex_index(i,1))%lon
    ENDDO
    DO i = 1, i_ne
      zv2dx(3,i) = p%verts%vertex(p%edges%vertex_index(i,2))%lon
    ENDDO
    DO i = 1, i_ne
      IF (p%edges%cell_index(i,1) > 0) THEN
        zv2dx(4,i) = p%cells%center(p%edges%cell_index(i,1))%lon
      ELSE
        zv2dx(4,i) = 0._wp
      ENDIF
    ENDDO
    DO i = 1, i_ne
      IF (p%edges%cell_index(i,2) > 0) THEN
        zv2dx(2,i) = p%cells%center(p%edges%cell_index(i,2))%lon
      ELSE
        zv2dx(2,i) = 0._wp
      ENDIF
    ENDDO
    DO i = 1, i_ne
      zv2dy(1,i) = p%verts%vertex(p%edges%vertex_index(i,1))%lat
    ENDDO
    DO i = 1, i_ne
      zv2dy(3,i) = p%verts%vertex(p%edges%vertex_index(i,2))%lat
    ENDDO
    DO i = 1, i_ne
      IF (p%edges%cell_index(i,1) > 0) THEN
        zv2dy(4,i) = p%cells%center(p%edges%cell_index(i,1))%lat
      ELSE
        zv2dy(4,i) = 0._wp
      ENDIF
    ENDDO
    DO i = 1, i_ne
      IF (p%edges%cell_index(i,2) > 0) THEN
        zv2dy(2,i) = p%cells%center(p%edges%cell_index(i,2))%lat
      ELSE
        zv2dy(2,i) = 0._wp
      ENDIF
    ENDDO
    WHERE (ABS(zv2dx(:,:)) < EPSILON(0.0_wp))
    zv2dx(:,:) = 0.0_wp
    endwhere
    WHERE (ABS(zv2dy(:,:)) < EPSILON(0.0_wp))
    zv2dy(:,:) = 0.0_wp
    endwhere
    DO j = 1, i_ne
      DO i = 1, 4
        IF ( ABS(zv2dy(i,j)) > 0.5_wp*pi-EPSILON(0.0_wp)) THEN
          zv2dx(i,j) = p%edges%center(j)%lon
        ENDIF
      ENDDO
    ENDDO
    DO i = 1, i_ne
      IF ( check_orientation(p%edges%center(i)%lon, &
        & zv2dx(:,i),zv2dy(:,i),4) < 0 ) THEN
        swap(1:4) = zv2dx(4:1:-1,i)
        zv2dx(:,i) = swap(:)
        swap(1:4) = zv2dy(4:1:-1,i)
        zv2dy(:,i) = swap(:)
      ENDIF
      IF (check_orientation(p%edges%center(i)%lon, &
        & zv2dx(:,i),zv2dy(:,i),4) < 0) THEN
      ENDIF
    ENDDO
    !
    CALL nf(nf_put_var_double(ncid, varid_elonv, zv2dx))
    CALL nf(nf_put_var_double(ncid, varid_elatv, zv2dy))
    !
    DEALLOCATE (zv2dx, zv2dy)
    !
    !--------------------------------------------------------------------------------------
    !
    CALL nf(nf_put_var_double(ncid, varid_carea, p%cells%area))
    CALL nf(nf_put_var_double(ncid, varid_varea, p%verts%dual_area))
    !
    !--------------------------------------------------------------------------------------
    !
    CALL nf(nf_put_var_double(ncid, varid1,  p%cells%center(:)%lon))
    CALL nf(nf_put_var_double(ncid, varid2,  p%cells%center(:)%lat))
    CALL nf(nf_put_var_double(ncid, varid3,  p%verts%vertex(:)%lon))
    CALL nf(nf_put_var_double(ncid, varid4,  p%verts%vertex(:)%lat))
    CALL nf(nf_put_var_double(ncid, varid5,  p%edges%center(:)%lon))
    CALL nf(nf_put_var_double(ncid, varid6,  p%edges%center(:)%lat))
    CALL nf(nf_put_var_int   (ncid, varid7,  p%cells%edge_index))
    CALL nf(nf_put_var_int   (ncid, varid8,  p%cells%vertex_index))
    CALL nf(nf_put_var_int   (ncid, varid9,  p%edges%cell_index))
    CALL nf(nf_put_var_int   (ncid, varid10, p%edges%vertex_index))
    CALL nf(nf_put_var_int   (ncid, varid11, p%verts%cell_index))
    CALL nf(nf_put_var_int   (ncid, varid12, p%verts%edge_index))
    CALL nf(nf_put_var_int   (ncid, varid13, p%verts%neighbor_index))
    CALL nf(nf_put_var_double(ncid, varid14, p%cells%area))
    CALL nf(nf_put_var_double(ncid, varid15, p%verts%dual_area))
    CALL nf(nf_put_var_double(ncid, varid16, p%edges%primal_edge_length))
    CALL nf(nf_put_var_double(ncid, varid17, p%edges%dual_edge_length))
    CALL nf(nf_put_var_double(ncid, varid18, p%edges%edge_vert_length))
    CALL nf(nf_put_var_double(ncid, varid40, p%edges%edge_cell_length))
    CALL nf(nf_put_var_double(ncid, varid19, p%edges%primal_normal(:)%v1))
    CALL nf(nf_put_var_double(ncid, varid20, p%edges%primal_normal(:)%v2))
    CALL nf(nf_put_var_double(ncid, varid21, p%edges%dual_normal(:)%v1))
    CALL nf(nf_put_var_double(ncid, varid22, p%edges%dual_normal(:)%v2) )
    CALL nf(nf_put_var_int   (ncid, varid23, p%cells%edge_orientation))
    CALL nf(nf_put_var_int   (ncid, varid24, p%cells%idx))
    CALL nf(nf_put_var_int   (ncid, varid26, p%cells%neighbor_index))
    CALL nf(nf_put_var_int   (ncid, varid28, p%edges%idx))
    CALL nf(nf_put_var_int   (ncid, varid29, p%verts%idx))
    CALL nf(nf_put_var_int   (ncid, varid30, p%verts%edge_orientation))
    CALL nf(nf_put_var_int   (ncid, varid31, p%edges%system_orientation))

  END SUBROUTINE write_standard_variables

  SUBROUTINE write_connectivity_variables

    CALL nf(nf_put_var_int   (ncid, varid25, p%cells%parent_index))
    CALL nf(nf_put_var_int   (ncid, varid27, p%cells%child_index))
    CALL nf(nf_put_var_int   (ncid, varid41, p%cells%child_id))
    CALL nf(nf_put_var_int   (ncid, varid38, p%edges%parent_index))
    CALL nf(nf_put_var_int   (ncid, varid39, p%edges%child_index))
    CALL nf(nf_put_var_int   (ncid, varid42, p%edges%child_id))
    CALL nf(nf_put_var_int   (ncid, varid32, p%cells%refin_ctrl))
    CALL nf(nf_put_var_int   (ncid, varid33, p%cells%indlist))
    CALL nf(nf_put_var_int   (ncid, varid43, p%cells%start_idx))
    CALL nf(nf_put_var_int   (ncid, varid44, p%cells%end_idx))
    CALL nf(nf_put_var_int   (ncid, varid34, p%edges%refin_ctrl))
    CALL nf(nf_put_var_int   (ncid, varid35, p%edges%indlist))
    CALL nf(nf_put_var_int   (ncid, varid45, p%edges%start_idx))
    CALL nf(nf_put_var_int   (ncid, varid46, p%edges%end_idx))
    CALL nf(nf_put_var_int   (ncid, varid36, p%verts%refin_ctrl))
    CALL nf(nf_put_var_int   (ncid, varid37, p%verts%indlist))
    CALL nf(nf_put_var_int   (ncid, varid47, p%verts%start_idx))
    CALL nf(nf_put_var_int   (ncid, varid48, p%verts%end_idx))
    CALL nf(nf_put_var_int   (ncid, varid49, p%cells%phys_id))
    CALL nf(nf_put_var_int   (ncid, varid50, p%edges%phys_id))

  END SUBROUTINE write_connectivity_variables


  END SUBROUTINE write_patch


  !EOC
  !-------------------------------------------------------------------------
  !
  !
  !>
  !!  Writes reduced data for grid levels below the global computational domain.
  !!
  !!
  !! @par Revision History
  !!  Guenther Zaengl, DWD, 2009-07-27
  !!  USES:
  !!
  SUBROUTINE write_reduced_output(g)
    !
    TYPE(t_grid) :: g

    !---local variables
    CHARACTER(LEN=filename_max) :: filename

    INTEGER :: i_lev

!!$    INTEGER, PARAMETER :: patch_unit = 89

    ! Copied from output_grid...

    INTEGER :: i_nc, i_ne, i_nv
    INTEGER :: old_mode
    INTEGER :: ncid
    INTEGER :: dimids(2)

    INTEGER :: dim_ncell, dim_nvertex, dim_nedge
    INTEGER :: dim_nvertex_per_cell, dim_ncells_per_edge, dim_nedges_per_vertex, &
      & dim_nchilds_per_cell

    INTEGER :: varid1, varid2, varid3, varid4, varid5, varid6, varid7, varid8, varid9

    INTEGER :: istat
    CHARACTER(LEN=256) :: command_line

    INTEGER :: command_line_len

    !-------------------------------------------------------------------------
    !BOC


    i_nc = g%ncells
    i_ne = g%nedges
    i_nv = g%nverts
    i_lev = g%level


    WRITE(filename,'(a,i0,a,i2.2,a)') &
      & 'iconR', grid_root, 'B', g%level, '.nc'

    WRITE(message_text,'(a,a)') 'Write gridmap file: ', TRIM(filename)
    CALL message ('', message_text)

    CALL GET_COMMAND(command_line, command_line_len, istat)
    CALL nf(nf_set_default_format(nf_format_64bit, old_mode))
    CALL nf(nf_create(TRIM(filename), nf_clobber, ncid))
    CALL nf(nf_set_fill(ncid, nf_nofill, old_mode))

    ! Global attributes
    !
    CALL nf(nf_put_att_text    (ncid, nf_global, 'title', 21, 'ICON grid description'))
    CALL nf(nf_put_att_text    (ncid, nf_global, 'history', command_line_len, command_line))
    CALL nf(nf_put_att_text    (ncid, nf_global, 'institution', 59, &
      & 'Max Planck Institute for Meteorology/Deutscher Wetterdienst'))
    CALL nf(nf_put_att_text    (ncid, nf_global, 'source', 10, 'icon-dev'))
    CALL nf(nf_put_att_text    (ncid, nf_global, 'grid_mapping_name' , 18, 'lat_long_on_sphere'))
    CALL nf(nf_put_att_text    (ncid, nf_global, 'crs_id' , 28, 'urn:ogc:def:cs:EPSG:6.0:6422'))
    CALL nf(nf_put_att_text    (ncid, nf_global, 'crs_name',30, 'Spherical 2D Coordinate System'))
    CALL nf(nf_put_att_text    (ncid, nf_global, 'ellipsoid_name' , 6, 'Sphere'))
    CALL nf(nf_put_att_double  (ncid, nf_global, 'semi_major_axis' , nf_double, 1, earth_radius))
    CALL nf(nf_put_att_double  (ncid, nf_global, 'inverse_flattening' , nf_double, 1, 0.0_wp))
    CALL nf(nf_put_att_int     (ncid, nf_global, 'grid_level', nf_int, 1, i_lev))
    CALL nf(nf_put_att_int     (ncid, nf_global, 'grid_root', nf_int, 1, grid_root))

    ! Dimensions
    !
    CALL nf(nf_def_dim(ncid, 'cell',   i_nc, dim_ncell))
    CALL nf(nf_def_dim(ncid, 'vertex', i_nv, dim_nvertex))
    CALL nf(nf_def_dim(ncid, 'edge',   i_ne, dim_nedge))
    !
    CALL nf(nf_def_dim(ncid, 'nc',    2, dim_ncells_per_edge))
    CALL nf(nf_def_dim(ncid, 'nv',    3, dim_nvertex_per_cell))
    CALL nf(nf_def_dim(ncid, 'ne',    6, dim_nedges_per_vertex))
    !
    CALL nf(nf_def_dim(ncid, 'no',    4, dim_nchilds_per_cell))


    ! Grid variables
    !
    !---------------------------------------------------------------------
    !
    CALL nf(nf_def_var(ncid, 'lon_cell_centre', nf_double, 1, dim_ncell, varid1))
    CALL nf(nf_put_att_text(ncid, varid1, 'long_name', 24, 'longitude of cell centre'))
    CALL nf(nf_put_att_text(ncid, varid1, 'units', 6, 'radian'))
    CALL nf(nf_put_att_text(ncid, varid1, 'cdi', 6, 'ignore'))
    !
    CALL nf(nf_def_var(ncid, 'lat_cell_centre', nf_double, 1, dim_ncell, varid2))
    CALL nf(nf_put_att_text(ncid, varid2, 'long_name', 23, 'latitude of cell centre'))
    CALL nf(nf_put_att_text(ncid, varid2, 'units', 6, 'radian'))
    CALL nf(nf_put_att_text(ncid, varid2, 'cdi', 6, 'ignore'))
    !
    CALL nf(nf_def_var(ncid, 'longitude_vertices', nf_double, 1, dim_nvertex, varid3))
    CALL nf(nf_put_att_text(ncid, varid3, 'long_name', 21, 'longitude of vertices'))
    CALL nf(nf_put_att_text(ncid, varid3, 'units', 6, 'radian'))
    CALL nf(nf_put_att_text(ncid, varid3, 'cdi', 6, 'ignore'))
    !
    CALL nf(nf_def_var(ncid, 'latitude_vertices', nf_double, 1, dim_nvertex, varid4))
    CALL nf(nf_put_att_text(ncid, varid4, 'long_name', 20, 'latitude of vertices'))
    CALL nf(nf_put_att_text(ncid, varid4, 'units', 6, 'radian'))
    CALL nf(nf_put_att_text(ncid, varid4, 'cdi', 6, 'ignore'))
    !
    dimids = (/ dim_ncell, dim_nvertex_per_cell /)
    CALL nf(nf_def_var(ncid, 'vertex_of_cell', nf_int, 2, dimids, varid5))
    CALL nf(nf_put_att_text(ncid, varid5, 'long_name', 32, 'vertices of each triangular cell'))
    CALL nf(nf_put_att_text(ncid, varid5, 'cdi', 6, 'ignore'))
    !
    dimids = (/ dim_nvertex, dim_nedges_per_vertex /)
    CALL nf(nf_def_var(ncid, 'cells_of_vertex', nf_int, 2, dimids, varid6))
    CALL nf(nf_put_att_text(ncid, varid6, 'long_name', 24, 'cells around each vertex'))
    CALL nf(nf_put_att_text(ncid, varid6, 'cdi', 6, 'ignore'))
    !
    CALL nf(nf_def_var(ncid, 'parent_cell_index', nf_int, 1, dim_ncell, varid7))
    CALL nf(nf_put_att_text(ncid, varid7, 'long_name', 17, 'parent cell index'))
    CALL nf(nf_put_att_text(ncid, varid7, 'cdi', 6, 'ignore'))
    !
    dimids = (/ dim_ncell, dim_nvertex_per_cell /)
    CALL nf(nf_def_var(ncid, 'neighbor_cell_index', nf_int, 2, dimids, varid8))
    CALL nf(nf_put_att_text(ncid, varid8, 'long_name', 19, 'cell neighbor index'))
    CALL nf(nf_put_att_text(ncid, varid8, 'cdi', 6, 'ignore'))
    !
    dimids = (/ dim_ncell, dim_nchilds_per_cell /)
    CALL nf(nf_def_var(ncid, 'child_cell_index', nf_int, 2, dimids, varid9))
    CALL nf(nf_put_att_text(ncid, varid9, 'long_name', 16, 'child cell index'))
    CALL nf(nf_put_att_text(ncid, varid9, 'cdi', 6, 'ignore'))
    !
    CALL nf(nf_enddef(ncid))
    !
    !--------------------------------------------------------------------------------------
    !
    CALL nf(nf_put_var_double(ncid, varid1,  g%cells%center(:)%lon))
    CALL nf(nf_put_var_double(ncid, varid2,  g%cells%center(:)%lat))
    CALL nf(nf_put_var_double(ncid, varid3,  g%verts%vertex(:)%lon))
    CALL nf(nf_put_var_double(ncid, varid4,  g%verts%vertex(:)%lat))
    CALL nf(nf_put_var_int   (ncid, varid5,  g%cells%vertex_index))
    CALL nf(nf_put_var_int   (ncid, varid6,  g%verts%cell_index))
    CALL nf(nf_put_var_int   (ncid, varid7,  g%cells%parent_index))
    CALL nf(nf_put_var_int   (ncid, varid8,  g%cells%neighbor_index))
    CALL nf(nf_put_var_int   (ncid, varid9,  g%cells%child_index))
    !
    !--------------------------------------------------------------------------------------

    CALL nf(nf_close(ncid))

  END SUBROUTINE write_reduced_output


  !EOC
  !-------------------------------------------------------------------------
  !
  !
  !>
  !!  Plots the grid portion associated with a local patch.
  !!
  !!
  !! @par Revision History
  !!  Luca Bonaventura, MPI-M, Hamburg, August 2005
  !! @par
  !!  Cleanup by Luca Bonaventura, Polimi, February 2006:
  !!  upper grid argument removed and printing of special
  !!  grid items on fort.* files removed
  !!
  SUBROUTINE plot_local_domains(p, g)
    !


    TYPE(t_patch),POINTER :: p
    TYPE(t_grid), INTENT(in) :: g
    !    TYPE(grid), INTENT(IN) :: gup


    INTEGER ::j    =0
    INTEGER ::ist  =0
    REAL(wp) ::zlon =0.0_wp
    REAL(wp) ::zlat =0.0_wp
    CHARACTER(LEN=max_length) :: filename
    CHARACTER(LEN=max_length) :: string1,string2,string3,string4
    !-------------------------------------------------------------------------
    !BOC


    WRITE(string1,'(i2.2)') p%id

    WRITE(string2,'(30a)') TRIM(ADJUSTL('cells'))

    WRITE(string3,'(i2.2)') p%level

    WRITE(string4,'(50a)') TRIM(ADJUSTL(string2)),TRIM(ADJUSTL('_dom')),&
      & TRIM(ADJUSTL(string1)),TRIM(ADJUSTL('lev')),&
      & TRIM(ADJUSTL(string3)),TRIM(ADJUSTL('.gmt'))
    WRITE(filename,'(50a)') TRIM(ADJUSTL(string4))

    OPEN(ngmt,FILE=filename, IOSTAT=ist)
    IF(ist/=0)THEN
      WRITE(*,*)' plot_local_patch: opening file:', filename ,' failed'
      STOP
    ENDIF
    print*, 'Write plotting output in ',TRIM(filename)
    DO j=1,p%ncells
      zlon=g%verts%vertex(g%cells%vertex_index(p%ic(j),1))%lon
      zlat=g%verts%vertex(g%cells%vertex_index(p%ic(j),1))%lat
      WRITE (ngmt,'(5x,2f8.1)') zlon*rad2deg, zlat*rad2deg
      zlon=g%verts%vertex(g%cells%vertex_index(p%ic(j),2))%lon
      zlat=g%verts%vertex(g%cells%vertex_index(p%ic(j),2))%lat
      WRITE (ngmt,'(5x,2f8.1)') zlon*rad2deg, zlat*rad2deg
      zlon=g%verts%vertex(g%cells%vertex_index(p%ic(j),3))%lon
      zlat=g%verts%vertex(g%cells%vertex_index(p%ic(j),3))%lat
      WRITE (ngmt,'(5x,2f8.1)') zlon*rad2deg, zlat*rad2deg
      WRITE (ngmt,'(a)') '>'
    ENDDO

    CLOSE(ngmt)

  END SUBROUTINE plot_local_domains

  !-------------------------------------------------------------------------
  !
  !>
  !!  Writes the index lists for refined model domains.
  !!
  !! Subdomains
  !!  can be either circular or rectangular (in lat lon coordinates).
  !!  In principle, arbitrary configurations of nested domains are
  !!  possible, but the maximum number of model domains is currently
  !!  set to 10. Overlapping nests are prevented automatically.
  !!  Replaces the former routine make_lists_local_patch
  !!
  RECURSIVE SUBROUTINE setup_index_lists(gg,curr_id)
    !
    ! !REVISION HISTORY (make_lists_local_patch):
    !  Luca Bonaventura, MPI-M, Hamburg, August 2005
    !  Modified by Peter Korn, MPI-M, Hamburg
    !
    !  Guenther Zaengl, DWD, March-June 2008:
    !  Correct bug in level indexing
    !  Enhance width of boundary zone in case of multiple nesting
    !  Add calculations needed for mesh refinement control
    !
    !  Guenther Zaengl, DWD, 2009-07-22:
    !  Rename to setup_index_lists for new version supporting arbitrary
    !    nesting configurations

    INTEGER, INTENT(in) :: curr_id
    TYPE(t_grid), TARGET, INTENT(inout):: gg(start_lev:end_lev)

    !---local variables
    INTEGER :: jn            =0
    INTEGER :: i_nc          =0
    INTEGER :: i_patch_count =0
    INTEGER :: i_parent_count=0
    REAL(wp) :: long0         =0.0_wp
    REAL(wp) :: lat0          =0.0_wp
    REAL(wp) :: long1         =0.0_wp
    REAL(wp) :: lat1          =0.0_wp
    REAL(wp) :: long2         =0.0_wp
    REAL(wp) :: lat2          =0.0_wp
    REAL(wp) :: r1            =0.0_wp
    REAL(wp) :: r2            =0.0_wp
    INTEGER :: ist           =0
    INTEGER :: i, j, i1, i2, i3, jn2, irv, iref, n_neigh_count, &
      & i_edp(3), i_edc(4,3), nc, nct, nct2, k, ipc_save
    INTEGER :: ilev, nest_id, nest_lev, inest

    REAL(wp) :: auxp(4), auxc(4), mindist, dist(3,3,4), orient(3,3,4), &
      & auxpc(3,3), auxcc(4,3,3)

    REAL(wp) :: pollat, pollon

    TYPE(t_grid_cells),POINTER :: cc =>null()
    TYPE(t_grid_cells),POINTER :: ccp =>null()
    TYPE(t_grid_vertices),POINTER :: cv =>null()
    TYPE(t_grid_edges),POINTER :: ce =>null()
    TYPE(t_grid_edges),POINTER :: cep =>null()
    TYPE(t_cell_list),POINTER :: ptr_list=>null()
    TYPE(t_cell_list),POINTER :: ptr_parent_list=>null()
    TYPE(t_cell_list),POINTER :: ptr_inv_list=>null()
    TYPE(t_cell_list),POINTER :: ptr_inv_parent_list=>null()
    TYPE(t_cartesian_coordinates) :: aux_epc, aux_ecc
    !-------------------------------------------------------------------------
    !BOC

    ! Set Domain ID of global domain to 1
    IF (curr_id == 1) THEN
      ilev = phys_grid_level(curr_id)
      gg(ilev)%cells%curr_id(:) = 1
    ENDIF

    ! Loop over subdomains of current domain
    nest_loop: DO inest = 1, n_phys_childdom(curr_id)

      ilev = phys_grid_level(curr_id)
      i_nc= gg(ilev)%ncells
      nest_id = phys_child_id(curr_id,inest)
      nest_lev = ilev + 1

      ptr_list            => ptr_phys_cell_list(nest_id)
      ptr_parent_list     => ptr_phys_parent_cell_list(nest_id)
      ptr_inv_list        => ptr_phys_inv_cell_list(nest_id)
      ptr_inv_parent_list => ptr_phys_inv_parent_cell_list(nest_id)

      WRITE(message_text,'(a,3i5)')'domain ID, parent grid level, nest level, ',&
        & nest_id,ilev,nest_lev
      CALL message ('', message_text)

      ALLOCATE(ptr_list%ip(1:4*i_nc), ptr_inv_list%ip(1:4*i_nc), stat=ist)
      IF(ist/=0)THEN
        WRITE(message_text,'(a,i5)')'allocate list failed', ist
        CALL message ('', message_text)
      ENDIF
      ALLOCATE(ptr_parent_list%ip(1:i_nc), ptr_inv_parent_list%ip(1:i_nc), stat=ist)
      IF(ist/=0)THEN
        WRITE(message_text,'(a,i5)')'allocate parent list failed', ist
        CALL message ('', message_text)
      ENDIF

      cc=>gg(ilev)%cells
      cv=>gg(ilev)%verts
      i_patch_count=0
      i_parent_count=0

      ptr_list%ip(:) = 0
      ptr_parent_list%ip(:) = 0
      ptr_inv_list%ip(:) = 0
      ptr_inv_parent_list%ip(:) = 0

      IF(l_circ)THEN

        long2 = center_lon(nest_id-1)/rad2deg
        lat2  = center_lat(nest_id-1)/rad2deg
        r2    = radius(nest_id-1)/rad2deg

        DO jn = 1, i_nc

          ! This ensures sufficient distance from the boundary of the parent domain
          ! and excludes the possibility of overlapping nests
          IF (cc%refin_ctrl(jn) /= 0) CYCLE
          IF (cc%curr_id(jn) /= curr_id) CYCLE
          IF (cc%refin_ctrl(cc%neighbor_index(jn,1)) /= 0) CYCLE
          IF (cc%refin_ctrl(cc%neighbor_index(jn,2)) /= 0) CYCLE
          IF (cc%refin_ctrl(cc%neighbor_index(jn,3)) /= 0) CYCLE
          IF (cv%refin_ctrl(cc%vertex_index(jn,1)) /= 0) CYCLE
          IF (cv%refin_ctrl(cc%vertex_index(jn,2)) /= 0) CYCLE
          IF (cv%refin_ctrl(cc%vertex_index(jn,3)) /= 0) CYCLE

          long0 = cc%center(jn)%lon
          lat0  = cc%center(jn)%lat
          r1=ACOS(SIN(lat2)*sin(lat0)+COS(lat2)*cos(lat0)*cos(long0-long2))

          IF (r1<r2) THEN

            i_parent_count= i_parent_count+1
            ptr_parent_list%ip(i_parent_count) = jn
            ptr_inv_parent_list%ip(jn) = i_parent_count

            i_patch_count=i_patch_count+1
            ptr_list%ip(i_patch_count) = cc%child_index(jn,1)
            ptr_inv_list%ip(cc%child_index(jn,1)) = i_patch_count
            gg(nest_lev)%cells%curr_id(cc%child_index(jn,1)) = nest_id

            i_patch_count=i_patch_count+1
            ptr_list%ip(i_patch_count) = cc%child_index(jn,2)
            ptr_inv_list%ip(cc%child_index(jn,2)) = i_patch_count
            gg(nest_lev)%cells%curr_id(cc%child_index(jn,2)) = nest_id

            i_patch_count=i_patch_count+1
            ptr_list%ip(i_patch_count) = cc%child_index(jn,3)
            ptr_inv_list%ip(cc%child_index(jn,3)) = i_patch_count
            gg(nest_lev)%cells%curr_id(cc%child_index(jn,3)) = nest_id

            i_patch_count=i_patch_count+1
            ptr_list%ip(i_patch_count) = cc%child_index(jn,4)
            ptr_inv_list%ip(cc%child_index(jn,4)) = i_patch_count
            gg(nest_lev)%cells%curr_id(cc%child_index(jn,4)) = nest_id

          ENDIF
        ENDDO

      ELSE IF (.not.l_rotate) THEN

        long2 = center_lon(nest_id-1)/rad2deg
        lat2  = center_lat(nest_id-1)/rad2deg
        long1 = hwidth_lon(nest_id-1)/rad2deg
        lat1  = hwidth_lat(nest_id-1)/rad2deg

        DO jn = 1, i_nc

          ! This ensures sufficient distance from the boundary of the parent domain
          ! and excludes the possibility of overlapping nests
          IF (cc%refin_ctrl(jn) /= 0) CYCLE
          IF (cc%curr_id(jn) /= curr_id) CYCLE
          IF (cc%refin_ctrl(cc%neighbor_index(jn,1)) /= 0) CYCLE
          IF (cc%refin_ctrl(cc%neighbor_index(jn,2)) /= 0) CYCLE
          IF (cc%refin_ctrl(cc%neighbor_index(jn,3)) /= 0) CYCLE
          IF (cv%refin_ctrl(cc%vertex_index(jn,1)) /= 0) CYCLE
          IF (cv%refin_ctrl(cc%vertex_index(jn,2)) /= 0) CYCLE
          IF (cv%refin_ctrl(cc%vertex_index(jn,3)) /= 0) CYCLE

          long0 = cc%center(jn)%lon
          lat0  = cc%center(jn)%lat
          IF (long0-long2<-pi) THEN
            r1 = (long0+2._wp*pi-long2)
          ELSE IF (long0-long2 > pi) THEN
            r1 = (long0-2._wp*pi-long2 )
          ELSE
            r1=(long0-long2)
          ENDIF
          r2 = lat0 -lat2

          IF ((ABS(r2)<=lat1).and.(ABS(r1)<=long1)) THEN

            i_parent_count= i_parent_count+1
            ptr_parent_list%ip(i_parent_count) = jn
            ptr_inv_parent_list%ip(jn) = i_parent_count

            i_patch_count=i_patch_count+1
            ptr_list%ip(i_patch_count) = cc%child_index(jn,1)
            ptr_inv_list%ip(cc%child_index(jn,1)) = i_patch_count
            gg(nest_lev)%cells%curr_id(cc%child_index(jn,1)) = nest_id

            i_patch_count=i_patch_count+1
            ptr_list%ip(i_patch_count) = cc%child_index(jn,2)
            ptr_inv_list%ip(cc%child_index(jn,2)) = i_patch_count
            gg(nest_lev)%cells%curr_id(cc%child_index(jn,2)) = nest_id

            i_patch_count=i_patch_count+1
            ptr_list%ip(i_patch_count) = cc%child_index(jn,3)
            ptr_inv_list%ip(cc%child_index(jn,3)) = i_patch_count
            gg(nest_lev)%cells%curr_id(cc%child_index(jn,3)) = nest_id

            i_patch_count=i_patch_count+1
            ptr_list%ip(i_patch_count) = cc%child_index(jn,4)
            ptr_inv_list%ip(cc%child_index(jn,4)) = i_patch_count
            gg(nest_lev)%cells%curr_id(cc%child_index(jn,4)) = nest_id

          ENDIF
        ENDDO

      ELSE

        long2 = center_lon(nest_id-1)/rad2deg
        lat2  = center_lat(nest_id-1)/rad2deg
        long1 = hwidth_lon(nest_id-1)/rad2deg
        lat1  = hwidth_lat(nest_id-1)/rad2deg

        ! Rotate center point into the equator
        IF (lat2 >= 0._wp) THEN
          pollat = lat2 - pi2/4._wp
        ELSE
          pollat = lat2 + pi2/4._wp
        ENDIF
        pollon = long2

        CALL rotate_latlon( lat2, long2, pollat, pollon )

        DO jn = 1, i_nc

          ! This ensures sufficient distance from the boundary of the parent domain
          ! and excludes the possibility of overlapping nests
          IF (cc%refin_ctrl(jn) /= 0) CYCLE
          IF (cc%curr_id(jn) /= curr_id) CYCLE
          IF (cc%refin_ctrl(cc%neighbor_index(jn,1)) /= 0) CYCLE
          IF (cc%refin_ctrl(cc%neighbor_index(jn,2)) /= 0) CYCLE
          IF (cc%refin_ctrl(cc%neighbor_index(jn,3)) /= 0) CYCLE
          IF (cv%refin_ctrl(cc%vertex_index(jn,1)) /= 0) CYCLE
          IF (cv%refin_ctrl(cc%vertex_index(jn,2)) /= 0) CYCLE
          IF (cv%refin_ctrl(cc%vertex_index(jn,3)) /= 0) CYCLE

          long0 = cc%center(jn)%lon
          lat0  = cc%center(jn)%lat

          CALL rotate_latlon( lat0, long0, pollat, pollon )

          IF (long0-long2 < -pi) THEN
            r1 = (long0+2._wp*pi-long2)*cos(lat0)
          ELSE IF (long0-long2 > pi) THEN
            r1 = (long0-2._wp*pi-long2 )*cos(lat0)
          ELSE
            r1=(long0-long2)*cos(lat0)
          ENDIF
          r2 = lat0 -lat2

          IF ((ABS(r2)<=lat1).and.(ABS(r1)<=long1)) THEN

            i_parent_count= i_parent_count+1
            ptr_parent_list%ip(i_parent_count) = jn
            ptr_inv_parent_list%ip(jn) = i_parent_count

            i_patch_count=i_patch_count+1
            ptr_list%ip(i_patch_count) = cc%child_index(jn,1)
            ptr_inv_list%ip(cc%child_index(jn,1)) = i_patch_count
            gg(nest_lev)%cells%curr_id(cc%child_index(jn,1)) = nest_id

            i_patch_count=i_patch_count+1
            ptr_list%ip(i_patch_count) = cc%child_index(jn,2)
            ptr_inv_list%ip(cc%child_index(jn,2)) = i_patch_count
            gg(nest_lev)%cells%curr_id(cc%child_index(jn,2)) = nest_id

            i_patch_count=i_patch_count+1
            ptr_list%ip(i_patch_count) = cc%child_index(jn,3)
            ptr_inv_list%ip(cc%child_index(jn,3)) = i_patch_count
            gg(nest_lev)%cells%curr_id(cc%child_index(jn,3)) = nest_id

            i_patch_count=i_patch_count+1
            ptr_list%ip(i_patch_count) = cc%child_index(jn,4)
            ptr_inv_list%ip(cc%child_index(jn,4)) = i_patch_count
            gg(nest_lev)%cells%curr_id(cc%child_index(jn,4)) = nest_id

          ENDIF
        ENDDO
      ENDIF

      ! Fill "holes" in refinement zone, defined as triangle having two
      ! neighbors inside the refinement zone

      ipc_save = i_parent_count

      DO jn = 1, i_nc

        iref = 0
        IF (ptr_inv_parent_list%ip(jn) == 0) THEN
          i1 = ptr_inv_parent_list%ip(cc%neighbor_index(jn,1))
          i2 = ptr_inv_parent_list%ip(cc%neighbor_index(jn,2))
          i3 = ptr_inv_parent_list%ip(cc%neighbor_index(jn,3))
          IF ((i1 /= 0) .and. (i1 <= ipc_save)) iref = iref + 1
          IF ((i2 /= 0) .and. (i2 <= ipc_save)) iref = iref + 1
          IF ((i3 /= 0) .and. (i3 <= ipc_save)) iref = iref + 1
        ENDIF

        IF (iref >= 2) THEN

          IF (cc%refin_ctrl(jn) /= 0) CYCLE
          IF (cc%refin_ctrl(cc%neighbor_index(jn,1)) /= 0) CYCLE
          IF (cc%refin_ctrl(cc%neighbor_index(jn,2)) /= 0) CYCLE
          IF (cc%refin_ctrl(cc%neighbor_index(jn,3)) /= 0) CYCLE
          IF (cv%refin_ctrl(cc%vertex_index(jn,1)) /= 0) CYCLE
          IF (cv%refin_ctrl(cc%vertex_index(jn,2)) /= 0) CYCLE
          IF (cv%refin_ctrl(cc%vertex_index(jn,3)) /= 0) CYCLE

          i_parent_count= i_parent_count+1
          ptr_parent_list%ip(i_parent_count) = jn
          ptr_inv_parent_list%ip(jn) = i_parent_count

          i_patch_count=i_patch_count+1
          ptr_list%ip(i_patch_count) = cc%child_index(jn,1)
          ptr_inv_list%ip(cc%child_index(jn,1)) = i_patch_count
          gg(nest_lev)%cells%curr_id(cc%child_index(jn,1)) = nest_id

          i_patch_count=i_patch_count+1
          ptr_list%ip(i_patch_count) = cc%child_index(jn,2)
          ptr_inv_list%ip(cc%child_index(jn,2)) = i_patch_count
          gg(nest_lev)%cells%curr_id(cc%child_index(jn,2)) = nest_id

          i_patch_count=i_patch_count+1
          ptr_list%ip(i_patch_count) = cc%child_index(jn,3)
          ptr_inv_list%ip(cc%child_index(jn,3)) = i_patch_count
          gg(nest_lev)%cells%curr_id(cc%child_index(jn,3)) = nest_id

          i_patch_count=i_patch_count+1
          ptr_list%ip(i_patch_count) = cc%child_index(jn,4)
          ptr_inv_list%ip(cc%child_index(jn,4)) = i_patch_count
          gg(nest_lev)%cells%curr_id(cc%child_index(jn,4)) = nest_id
        ENDIF
      ENDDO

      n_phys_parent_count(nest_id) = i_parent_count
      n_phys_patch_count(nest_id) = i_patch_count


      ! Compute "refin_ctrl" field for mesh refinement control:
      ! Negative values mean that a nested grid is present at the location of a
      ! given cell/edge/vertex. The counting starts at -1 along the lateral boundary
      ! of the nest and then proceeds to successively more negative values with
      ! increasing distance from the lateral edge. It ends at the value of min_rlcell_int
      ! (min_rledge_int/min_rlvert_int); all points in the interior of a nest overlap area
      ! are flagged with these values.
      ! Positive values are set along the lateral boundaries of a given domain,
      ! starting with 1 for boundary cell/edges/vertices. The maximum values attained
      ! are max_rlcell/max_rledge/max_rlvert; then the refin_ctrl flag jumps to 0.
      ! For a global model domain, there are no positive refin_ctrl flags.
      !
      ! The general reasoning behind this system is that in the presence of negative
      ! refin_ctrl flag, downward interpolation and/or feedback takes place
      ! Positive flags (up to a certain value) indicate that tendencies are
      ! interpolated from the parent domain because not all operations can be
      ! executed on the local domain (because the lateral boundary is too close)
      ! They can also be used to control relaxation or filtering operations
      ! along the lateral model boundary and to specify the use of terrain blending


      ! Condition for the classification as a lateral boundary point (flag -1):
      ! A grid cell has child indices (implying that it is included in the parent
      ! list), but at least one of its vertices borders on a cell without children.

      cc=>gg(ilev)%cells
      cv=>gg(ilev)%verts

      DO i = 1, n_phys_parent_count(nest_id)
        jn = ptr_parent_list%ip(i) ! global index of cell point
        n_neigh_count = 0
        DO i1 = 1, 3
          j = cc%vertex_index(jn,i1)
          DO i2 = 1, 6
            k = cv%cell_index(j,i2)
            IF (k <= 0) THEN
              ! The 6th index is zero for pentagon points
              ! In this case, we also count the cell so as to arrive at a sum of 18
              n_neigh_count = n_neigh_count + 1
            ELSE IF (ptr_inv_parent_list%ip(k) /= 0) THEN
              ! adjacent cell lies within refinement area
              n_neigh_count = n_neigh_count + 1
            ENDIF
          ENDDO
        ENDDO
        IF (n_neigh_count < 18) THEN
          IF (cc%child_id(jn) == 0) THEN
            cc%refin_ctrl(jn) = -1
            cc%child_id(jn) = nest_id
          ELSE
            CALL finish &
              & ('setup_index_lists', 'Boundary zone too small')
          ENDIF
        ENDIF
      ENDDO


      ! For the lower (more negative) values of the refin_ctrl variable,
      ! we look for neighbors being flagged with the next higher one:

      DO irv = 2, ABS(min_rlcell_int)-1
        DO i = 1, n_phys_parent_count(nest_id)
          jn = ptr_parent_list%ip(i) ! global index of cell point
          IF (cc%refin_ctrl(jn) <= -1) CYCLE
          DO i1 = 1, 3
            j = cc%vertex_index(jn,i1)
            DO i2 = 1, 6
              k = cv%cell_index(j,i2)
              IF (k .eq. 0) EXIT
              IF (cc%refin_ctrl(k) == 1-irv) THEN
                cc%refin_ctrl(jn) = -irv
                cc%child_id(jn) = nest_id
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      ENDDO

      ! Finally consider all remaining points appearing in the parent list that
      ! have not yet been flagged. 

      DO i = 1, n_phys_parent_count(nest_id)
        jn = ptr_parent_list%ip(i) ! global index of cell point
        IF (cc%refin_ctrl(jn) == 0) THEN
          cc%refin_ctrl(jn) = min_rlcell_int
          cc%child_id(jn) = nest_id
        ENDIF
      ENDDO

      ! Now we conduct a similar search on the refined mesh, starting with a
      ! look for border grid cells

      cc=>gg(nest_lev)%cells
      cv=>gg(nest_lev)%verts

      DO i = 1, n_phys_patch_count(nest_id)
        jn = ptr_list%ip(i) ! global index of cell point
        n_neigh_count = 0
        DO i1 = 1, 3
          j = cc%vertex_index(jn,i1)
          DO i2 = 1, 6
            k = cv%cell_index(j,i2)
            IF (k <= 0) THEN
              ! The 6th index is zero for pentagon points
              ! In this case, we also count the cell so as to arrive at a sum of 18
              n_neigh_count = n_neigh_count + 1
            ELSE IF (ptr_inv_list%ip(k) /= 0) THEN
              ! adjacent cell lies within refinement area
              n_neigh_count = n_neigh_count + 1
            ENDIF
          ENDDO
        ENDDO
        IF (n_neigh_count < 18) cc%refin_ctrl(jn) = 1
      ENDDO

      ! For the higher values of the refin_ctrl variable, we look for
      ! neighbors being flagged with the next lower one:

      DO irv = 2,MAX(max_rlcell,bdy_indexing_depth)
        DO i = 1, n_phys_patch_count(nest_id)
          jn = ptr_list%ip(i) ! global index of cell point
          IF (cc%refin_ctrl(jn) >= 1) CYCLE
          DO i1 = 1, 3
            j = cc%vertex_index(jn,i1)
            DO i2 = 1, 6
              k = cv%cell_index(j,i2)
              IF (k .eq. 0) EXIT
              IF (cc%refin_ctrl(k) == irv-1) THEN
                IF (cc%refin_ctrl(jn) < 0) CALL finish &
                  & ('setup_index_lists', 'Boundary zone too small')
                cc%refin_ctrl(jn) = irv
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      ENDDO


      ! Consistency check Ia: For parent cells flagged with -1,
      ! the child cells should have 1 or 2

      DO i = 1, n_phys_parent_count(nest_id)
        jn = ptr_parent_list%ip(i)
        IF (gg(ilev)%cells%refin_ctrl(jn) == -1) THEN
          DO j = 4*(i-1)+1, 4*i
            jn2 = ptr_list%ip(j)
            IF ((gg(nest_lev)%cells%refin_ctrl(jn2) < 1) .or. &
              & (gg(nest_lev)%cells%refin_ctrl(jn2) > 2)) &
              & CALL finish ('setup_index_lists', &
              & 'refin_ctrl consistency error Ia')
          ENDDO

          ! Consistency check ib: For parent cells flagged with -2,
          ! the child cells should have 3 or 4

        ELSE IF (gg(ilev)%cells%refin_ctrl(jn) == -2) THEN
          DO j = 4*(i-1)+1, 4*i
            jn2 = ptr_list%ip(j)
            IF ((gg(nest_lev)%cells%refin_ctrl(jn2) < 3) .or. &
              & (gg(nest_lev)%cells%refin_ctrl(jn2) > 4)) &
              & CALL finish ('setup_index_lists', &
              & 'refin_ctrl consistency error Ib')
          ENDDO

        ENDIF
      ENDDO

      ! Now we compute refinement control flags for edges and vertices:
      ! For a vertex belonging to a certain cell, this flag is supposed to
      ! have the same value as cells%refin_ctrl, with priority given to the
      ! higher absolute value (which means that in the outermost row of
      ! triangles, the vertices bordering to row 2 have refin_ctrl = +/-2)
      !
      ! Edges located right at the lateral boundary of a refined mesh have a
      ! value of 1, whereas interior edges between cells with cells%refin_ctrl
      ! = 1 have a value of 2; next, edges between cells with cells%refin_ctrl
      ! = 1 and 2 have a value of 3, and so on...

      cc=>gg(nest_lev)%cells
      ce=>gg(nest_lev)%edges
      cv=>gg(nest_lev)%verts

      DO irv = 1,MAX(max_rlcell,bdy_indexing_depth)
        DO i = 1, n_phys_patch_count(nest_id)
          jn = ptr_list%ip(i) ! global index of cell point
          IF (cc%refin_ctrl(jn) == irv) THEN
            DO i1 = 1, 3
              k = cc%vertex_index(jn,i1)
              IF (cv%refin_ctrl(k) < 0) &
                & CALL finish ('setup_index_lists', &
                & 'Boundary zone too small')
              cv%refin_ctrl(k) = irv
            ENDDO
          ENDIF
        ENDDO
      ENDDO


      DO irv = 1,MAX(max_rlcell,bdy_indexing_depth-1)
        DO i = 1, n_phys_patch_count(nest_id)
          jn = ptr_list%ip(i) ! global index of cell point
          IF (cc%refin_ctrl(jn) == irv) THEN
            DO i1 = 1, 3
              j = cc%edge_index(jn,i1)
              i2 = ce%cell_index(j,1)
              i3 = ce%cell_index(j,2)
              IF (i2 /= jn) k = i2
              IF (i3 /= jn) k = i3
              IF (cc%refin_ctrl(k) == irv) THEN
                ce%refin_ctrl(j) = 2*irv
              ELSE IF (cc%refin_ctrl(k) == irv-1) THEN
                ce%refin_ctrl(j) = 2*irv-1
              ELSE IF (cc%refin_ctrl(k) == irv+1) THEN
                ce%refin_ctrl(j) = 2*irv+1
                !             ELSE IF ((irv == 1) .AND. (cc%refin_ctrl(k) == 999)) THEN
                !               ce%refin_ctrl(j) = 1  ! lateral boundary edge
              ELSE
                CALL finish ('setup_index_lists', &
                  & 'refin_ctrl consistency error IIa')
              ENDIF
            ENDDO
          ENDIF
        ENDDO
      ENDDO

      cc=>gg(ilev)%cells
      ce=>gg(ilev)%edges
      cv=>gg(ilev)%verts

      DO irv = -1,min_rlcell_int,-1
        DO i = 1, n_phys_parent_count(nest_id)
          jn = ptr_parent_list%ip(i) ! global index of cell point
          IF (cc%refin_ctrl(jn) == irv) THEN

            DO i1 = 1, 3
              k = cc%vertex_index(jn,i1)
              cv%refin_ctrl(k) = irv
              cv%child_id(k)   = nest_id
            ENDDO

            DO i1 = 1, 3
              j = cc%edge_index(jn,i1)
              i2 = ce%cell_index(j,1)
              i3 = ce%cell_index(j,2)
              IF (i2 /= jn) k = i2
              IF (i3 /= jn) k = i3
              IF (cc%refin_ctrl(k) == irv) THEN
                ce%refin_ctrl(j) = 2*irv
                ce%child_id(j)   = nest_id
              ELSE IF (cc%refin_ctrl(k) == irv-1) THEN
                ce%refin_ctrl(j) = 2*irv-1
                ce%child_id(j)   = nest_id
              ELSE IF (cc%refin_ctrl(k) == irv+1) THEN
                ce%refin_ctrl(j) = 2*irv+1
                ce%child_id(j)   = nest_id
              ELSE IF ((irv == -1) .and. cc%refin_ctrl(k) == max_rlcell) THEN
                ce%refin_ctrl(j) = -1
                ce%child_id(j)   = nest_id
              ELSE
                CALL finish ('setup_index_lists', &
                  & 'refin_ctrl consistency error IIb')
              ENDIF
            ENDDO

          ENDIF
        ENDDO
      ENDDO

      ! Calculate parent and child indices for edges

      cc=>gg(nest_lev)%cells
      ce=>gg(nest_lev)%edges
      ccp=>gg(ilev)%cells
      cep=>gg(ilev)%edges

      ! Logic of the child edge indices: For a given edge, children 1 and 2 are
      ! the sub-triangle edges coinciding with the parent edge; children 3 and 4
      ! are the edges of the inner sub-triangles having (roughly) the same
      ! orientation as the parent edge

      DO i = 1, n_phys_parent_count(nest_id)
        jn = ptr_parent_list%ip(i) ! global index of cell point
        IF (ccp%refin_ctrl(jn) < 0) THEN
          DO j = 1, 3
            i_edp(j) = ccp%edge_index(jn,j)
            auxp(1) = cep%center(i_edp(j))%lat
            auxp(2) = cep%center(i_edp(j))%lon
            auxp(3) = cep%primal_normal(i_edp(j))%v1
            auxp(4) = cep%primal_normal(i_edp(j))%v2
            CALL gvec2cvec(auxp(3),auxp(4),auxp(2),auxp(1),&
              & auxpc(j,1),auxpc(j,2),auxpc(j,3))
          ENDDO
          DO k = 1, 4
            i1 = ccp%child_index(jn,k)
            DO j = 1, 3
              i_edc(k,j) = cc%edge_index(i1,j)
              auxc(1) = ce%center(i_edc(k,j))%lat
              auxc(2) = ce%center(i_edc(k,j))%lon
              auxc(3) = ce%primal_normal(i_edc(k,j))%v1
              auxc(4) = ce%primal_normal(i_edc(k,j))%v2
              CALL gvec2cvec(auxc(3),auxc(4),auxc(2),auxc(1),&
                & auxcc(k,j,1),auxcc(k,j,2),auxcc(k,j,3))
            ENDDO
          ENDDO
          DO i1 = 1, 3
            DO k = 1, 4
              DO i2 = 1, 3
                aux_epc         = gc2cc(cep%center(i_edp(i1)))
                aux_ecc         = gc2cc(ce%center(i_edc(k,i2)))
                dist(i1,i2,k)   = arc_length(aux_epc,aux_ecc)
                orient(i1,i2,k) = DOT_PRODUCT(auxpc(i1,:),auxcc(k,i2,:))* &
                  & ABS(DOT_PRODUCT(auxpc(i1,:),auxcc(k,i2,:)))/ &
                  & DOT_PRODUCT(auxpc(i1,:),auxpc(i1,:))/ &
                  & DOT_PRODUCT(auxcc(k,i2,:),auxcc(k,i2,:))
              ENDDO
            ENDDO
          ENDDO
          mindist = MINVAL(dist)
          nct = 0
          nct2 = 0
          DO i1 = 1, 3
            nc = 0
            DO k = 1, 4
              DO i2 = 1, 3
                IF ( (dist(i1,i2,k) <= 1.5_wp*mindist) .and. &
                  & (orient(i1,i2,k) > 0.99_wp) ) THEN
                  nc = nc+1
                  nct = nct+1
                  cep%child_index(i_edp(i1),nc) = i_edc(k,i2)
                  ce%parent_index(i_edc(k,i2)) = i_edp(i1)
                  IF (nc > 2) CALL finish ('setup_index_lists', &
                    & 'Error 1 in calculation of parent/child edge indices')
                ELSE IF ( (dist(i1,i2,k) > 1.5_wp*mindist) .and. &
                  & (orient(i1,i2,k) > 0.99_wp) ) THEN
                  nct2 = nct2+1
                  IF (cep%child_index(i_edp(i1),3) == 0) THEN
                    cep%child_index(i_edp(i1),3) = i_edc(k,i2)
                    ce%parent_index(i_edc(k,i2)) = i_edp(i1)
                  ELSE IF (cep%child_index(i_edp(i1),3) /= i_edc(k,i2)) THEN
                    cep%child_index(i_edp(i1),4) = i_edc(k,i2)
                    ce%parent_index(i_edc(k,i2)) = i_edp(i1)
                  ENDIF
                ELSE IF ( (dist(i1,i2,k) > 1.5_wp*mindist) .and.  &
                  & (orient(i1,i2,k) < -0.99_wp) ) THEN
                  nct2 = nct2+1
                  IF (cep%child_index(i_edp(i1),3) == 0) THEN
                    cep%child_index(i_edp(i1),3) = i_edc(k,i2)
                    ce%parent_index(i_edc(k,i2)) = i_edp(i1)
                  ELSE IF (cep%child_index(i_edp(i1),3) /= i_edc(k,i2)) THEN
                    cep%child_index(i_edp(i1),4) = i_edc(k,i2)
                    ce%parent_index(i_edc(k,i2)) = i_edp(i1)
                  ENDIF
                ENDIF
              ENDDO
            ENDDO
          ENDDO
          ! Check sum of child edges
          IF (nct /= 6) CALL finish ('setup_index_lists', &
            & 'Error 2 in calculation of parent/child edge indices')
          ! There are actually 3 inner edges, but they are counted twice
          IF (nct2 /= 6) CALL finish ('setup_index_lists', &
            & 'Error 3 in calculation of parent/child edge indices')
        ENDIF
      ENDDO

      ! Consistency check: Are edge child indices set at all points with
      ! edges%refin_ctrl < 0?
      DO i = 1, gg(ilev)%nedges
        IF (cep%refin_ctrl(i) == -1) THEN
          DO j = 1,3
            jn = cep%child_index(i,j)
            IF (jn == 0) CALL finish ('setup_index_lists', &
              & 'Global child edge index inconsistency 1')
          ENDDO
        ELSE IF (cep%refin_ctrl(i) < -1) THEN
          DO j = 1,4
            jn = cep%child_index(i,j)
            IF (jn == 0) CALL finish ('setup_index_lists', &
              & 'Global child edge index inconsistency 2')
          ENDDO
        ENDIF
        DO j = 1,2
          jn = cep%child_index(i,j)
          IF (jn < 0) CALL finish ('setup_index_lists', &
            & 'Global child edge index inconsistency 3')
        ENDDO
      ENDDO

      WRITE(message_text,'(a,i5)')'diagnostic output for physical domain ', nest_id
      CALL message ('', message_text)
      WRITE(message_text,'(a,2i9)')'number of grid points and parent grid points:', &
        & n_phys_patch_count(nest_id),n_phys_parent_count(nest_id)
      CALL message ('', message_text)

      IF (n_phys_patch_count(nest_id) < 250 ) THEN
        CALL finish ('setup_index_lists', &
          & 'Unreasonably small size of nested domain')
      ENDIF

      IF (n_phys_childdom(nest_id) > 0) CALL setup_index_lists(gg, nest_id)

    ENDDO nest_loop


  END SUBROUTINE setup_index_lists

  !>
  !! This module merges physical domains into logical domains.
  !!
  !! If the number of logical model domains requested in the namelist
  !! is smaller than the number of physical domains, the accordingly required
  !! merging of the domains is accomplished here. For the most part, this
  !! requires concatenating the index lists defined in setup_index_lists
  !! and resetting the parent and child ID's
  !!
  !! @par Revision History
  !! Developed by Guenther Zaengl, DWD, 2010-07-05
  !!
  SUBROUTINE merge_nested_domains(gg)

    TYPE(t_grid), TARGET, INTENT(inout):: gg(start_lev:end_lev)

    !---local variables

    TYPE(t_grid_cells),POINTER :: cc =>null()
    TYPE(t_grid_cells),POINTER :: ccp =>null()
    TYPE(t_grid_edges),POINTER :: ce =>null()
    TYPE(t_grid_edges),POINTER :: cep =>null()
    TYPE(t_grid_vertices),POINTER :: cv =>null()
    TYPE(t_grid_vertices),POINTER :: cvp =>null()

    TYPE(t_cell_list),POINTER :: ptr_phys_list=>null()
    TYPE(t_cell_list),POINTER :: ptr_phys_parent_list=>null()
    TYPE(t_cell_list),POINTER :: ptr_phys_inv_list=>null()
    TYPE(t_cell_list),POINTER :: ptr_phys_inv_parent_list=>null()
    TYPE(t_cell_list),POINTER :: ptr_list=>null()
    TYPE(t_cell_list),POINTER :: ptr_parent_list=>null()
    TYPE(t_cell_list),POINTER :: ptr_inv_list=>null()
    TYPE(t_cell_list),POINTER :: ptr_inv_parent_list=>null()

    INTEGER :: idom, ilev, ilevp, i_nc, ipdom, np_patch, np_parent, jp, js, ist
    INTEGER :: ishift_patch, ishift_parent

    !-------------------------------------------------------------------------

    ! Loop over logical domains specified in namelist
    ! Nothing has to be done for the global domain because there exists only one
    dom_loop: DO idom = 2, n_dom

      ilev = grid_level(idom)
      ilevp = grid_level(parent_grid_id(idom))
      i_nc= gg(ilev)%ncells

      ! corresponding physical domain ID (the first one if more than one exist)
      ipdom = physdom_list(idom,1)

      ! Lists of logical domains that need to be filled now
      ptr_list            => ptr_cell_list(idom)
      ptr_parent_list     => ptr_parent_cell_list(idom)
      ptr_inv_list        => ptr_inv_cell_list(idom)
      ptr_inv_parent_list => ptr_inv_parent_cell_list(idom)

      ALLOCATE(ptr_list%ip(1:4*i_nc), ptr_inv_list%ip(1:4*i_nc), stat=ist)
      IF(ist/=0)THEN
        WRITE(message_text,'(a,i5)')'allocate list failed', ist
        CALL message ('', message_text)
      ENDIF
      ALLOCATE(ptr_parent_list%ip(1:i_nc), ptr_inv_parent_list%ip(1:i_nc), stat=ist)
      IF(ist/=0)THEN
        WRITE(message_text,'(a,i5)')'allocate parent list failed', ist
        CALL message ('', message_text)
      ENDIF

      ptr_list%ip(:) = 0
      ptr_parent_list%ip(:) = 0
      ptr_inv_list%ip(:) = 0
      ptr_inv_parent_list%ip(:) = 0

      cc=>gg(ilev)%cells
      ce=>gg(ilev)%edges
      cv=>gg(ilev)%verts

      ccp=>gg(ilevp)%cells
      cep=>gg(ilevp)%edges
      cvp=>gg(ilevp)%verts

      ! Lists of physical domains that already exist
      ptr_phys_list            => ptr_phys_cell_list(ipdom)
      ptr_phys_parent_list     => ptr_phys_parent_cell_list(ipdom)
      ptr_phys_inv_list        => ptr_phys_inv_cell_list(ipdom)
      ptr_phys_inv_parent_list => ptr_phys_inv_parent_cell_list(ipdom)

      ! For the first physical domain located in the current logical domain,
      ! a simple copying is sufficient
      np_patch  = n_phys_patch_count(ipdom)
      np_parent = n_phys_parent_count(ipdom)

      ptr_list%ip(1:np_patch)             = ptr_phys_list%ip(1:np_patch)
      ptr_parent_list%ip(1:np_parent)     = ptr_phys_parent_list%ip(1:np_parent)
      ptr_inv_list%ip(1:np_patch)         = ptr_phys_inv_list%ip(1:np_patch)
      ptr_inv_parent_list%ip(1:np_parent) = ptr_phys_inv_parent_list%ip(1:np_parent)

      n_patch_count(idom)  = n_phys_patch_count(ipdom)
      n_parent_count(idom) = n_phys_parent_count(ipdom)

      ! Reset child ID's of parent patch to logical value
      DO jp = 1, np_parent
        ccp%child_id(ptr_parent_list%ip(jp)) = idom
        cep%child_id(ccp%edge_index(ptr_parent_list%ip(jp),1:3)) = idom
        cvp%child_id(ccp%vertex_index(ptr_parent_list%ip(jp),1:3)) = idom
      ENDDO

      ! Set field for physical domain ID in current patch
      DO jp = 1, np_patch
        cc%phys_id(ptr_list%ip(jp)) = ipdom
        ce%phys_id(cc%edge_index(ptr_list%ip(jp),1:3)) = ipdom
        cv%phys_id(cc%vertex_index(ptr_list%ip(jp),1:3)) = ipdom
      ENDDO

      ! Loop over remaining physical domains included in current logical domain
      DO js = 2, n_physdom(idom)

        ipdom         = physdom_list(idom,js)
        ishift_patch  = n_patch_count(idom)
        ishift_parent = n_parent_count(idom)

        ! Reset pointers to physical domains
        ptr_phys_list            => ptr_phys_cell_list(ipdom)
        ptr_phys_parent_list     => ptr_phys_parent_cell_list(ipdom)
        ptr_phys_inv_list        => ptr_phys_inv_cell_list(ipdom)
        ptr_phys_inv_parent_list => ptr_phys_inv_parent_cell_list(ipdom)

        np_patch  = n_phys_patch_count(ipdom)
        np_parent = n_phys_parent_count(ipdom)

        ptr_list%ip(ishift_patch+1:ishift_patch+np_patch)             = &
          & ptr_phys_list%ip(1:np_patch)
        ptr_parent_list%ip(ishift_parent+1:ishift_parent+np_parent)   = &
          & ptr_phys_parent_list%ip(1:np_parent)

        DO jp = ishift_patch+1, ishift_patch+np_patch
          ptr_inv_list%ip(ptr_list%ip(jp)) = jp
        ENDDO

        DO jp = ishift_parent+1, ishift_parent+np_parent
          ptr_inv_parent_list%ip(ptr_parent_list%ip(jp)) = jp
        ENDDO

        n_patch_count(idom)  = n_patch_count(idom)  + n_phys_patch_count(ipdom)
        n_parent_count(idom) = n_parent_count(idom) + n_phys_parent_count(ipdom)

        ! Reset child ID's of parent patch to logical value
        DO jp = ishift_parent+1, ishift_parent+np_parent
          ccp%child_id(ptr_parent_list%ip(jp)) = idom
          cep%child_id(ccp%edge_index(ptr_parent_list%ip(jp),1:3)) = idom
          cvp%child_id(ccp%vertex_index(ptr_parent_list%ip(jp),1:3)) = idom
        ENDDO

        ! Set field for physical domain ID in current patch
        DO jp = ishift_patch+1, ishift_patch+np_patch
          cc%phys_id(ptr_list%ip(jp)) = ipdom
          ce%phys_id(cc%edge_index(ptr_list%ip(jp),1:3)) = ipdom
          cv%phys_id(cc%vertex_index(ptr_list%ip(jp),1:3)) = ipdom
        ENDDO

      ENDDO

      WRITE(message_text,'(a,i5)')'diagnostic output for logical domain ', idom
      CALL message ('', message_text)
      WRITE(message_text,'(a,2i9,i5)') &
        & 'number of grid points, parent grid points and physical domains:', &
        & n_patch_count(idom),n_parent_count(idom),n_physdom(idom)
      CALL message ('', message_text)

    ENDDO dom_loop

  END SUBROUTINE merge_nested_domains

  !-------------------------------------------------------------------------
  !>
  !!
  !! @par Revision History
  !! Luis Kornblueh, MPI-M, Hamburg, April 2013
  !! dump the xml table entry for grib2 handling of horizontal grids
  !!
  SUBROUTINE dump_grid_table_xml_metadata(filename, grid_number)

    CHARACTER(len=*), INTENT(in) :: filename
    INTEGER,          INTENT(in) :: grid_number ! = number of grid used valid for the current domain
    CHARACTER(len=255) :: uriname, uri_subcentre

    CALL message('','')
    CALL message('','---------- BEGIN XML grid table descriptor ----------')
    CALL message('','')

    ! grouping element: grid
    WRITE(message_text,'(a,i0,a,i0,a,i0,a)') &
          &          '<grid number_of_grid_used="', grid_number, &
          &              '" centre="', centre, &
          &              '" subcentre="', subcentre, &
          &              '">'
    CALL message('',message_text)

    IF (subcentre > 0) THEN
      WRITE(uri_subcentre,'(i0,a)') subcentre, '/'
    ELSE
      uri_subcentre = ''
    ENDIF
    SELECT CASE (centre)
    CASE(78)
      uriname = TRIM(uri_pathname)//'edzw/'//TRIM(uri_subcentre)
    CASE(255)
      uriname = TRIM(uri_pathname)//'mpim/'//TRIM(uri_subcentre)
    CASE DEFAULT
      uriname = TRIM(uri_pathname)//'unknown/'//TRIM(uri_subcentre)
    END SELECT

    ! contained element: uri
    WRITE(message_text,'(a)') '    <uri>'
    CALL message('',message_text)
    WRITE(message_text,'(a,a,a)') '        ', uri_pathname, TRIM(filename)
    CALL message('',message_text)
    WRITE(message_text,'(a)') '    </uri>'
    CALL message('',message_text)

    ! contained element: description
    WRITE(message_text,'(a)') '    <description>'
    CALL message('',message_text)
    WRITE(message_text,'(a)') '    </description>'
    CALL message('',message_text)

    ! finishing grouping element: grid
    WRITE(message_text,'(a)') '</grid>'
    CALL message('',message_text)

    CALL message('','')
    CALL message('','----------- END XML grid table descriptor -----------')
    CALL message('','')

  END SUBROUTINE dump_grid_table_xml_metadata

END MODULE mo_gridref

