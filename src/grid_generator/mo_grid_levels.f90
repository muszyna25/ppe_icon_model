!>
!!               The subroutines in this module compute.
!!
!!               The subroutines in this module compute
!!               all the geometric information for the grid hierarchy
!!               in the internal data representation used by the grid generator.
!!
!! @par Revision History
!! Initial version  by:
!! @par
!! Luis Kornblueh,  MPI-M, Hamburg, January 2004
!! @par
!! Luis Kornblueh, MPI-M, Hamburg, February 2004
!!  - change to compiling and running version
!!  - create level 0 subdivision, connect triangle pointers
!! Luis Kornblueh, MPI-M, Hamburg, March 2004
!!  - include parents
!!  - full working version of tree and graph representation
!! Luca Bonaventura,  MPI-M, Hamburg, October 2004
!!  - include documentation and Protex headers
!! Luca Bonaventura,  MPI-M, Hamburg, April 2005
!!  - revision to account for patches
!! Luca Bonaventura,  MPI-M, Hamburg, August 2005
!!  - revision to include grid optimization and
!!    correct computation of normal vectors to edges
!! Thomas Heinze,  DWD, Offenbach, 2005-09
!!  - cleaning up code according to latest programming guide
!!  - debugging for grid optimization: changed index of optimize_heikes
!!  Modifications by Th.Heinze, DWD (2006-10-25):
!!  - changed index to idx in TYPE declarations of edge_info, vertex_info,
!!    triangle_info, grid_cells, grid_edges and grid_vertices
!!  Modification by P. Ripodas, DWD, (2007-02):
!!  - The edge%system_orientation is calculated (system orientation
!!  of (v2-v1), (c2-c1)
!! Almut Gassmannn, MPI-M, (2007-03)
!!  - Equal area subdivision and small circle subdivision parts added
!! Luis Kornblueh, MPI-M, March 2008
!!  - removed unnecessary vertex0/1/2 references
!! Guenther Zaengl, DWD, 2008-12-08
!!  - add option for stopping grid optimization at a certain level
!!    (needed for combining spring dynamics with nesting)
!! Almut Gassmann, MPI-M (2009-01-15)
!!  - adding a planar option
!!  - cleaning and restructuring
!!  - Earth dimensions used already here for avoiding the patch generator
!!    in every case
!! Almut Gassmann, MPI-M (2009-03-02)
!!  - tuning parameter beta_spring for spring dynamics comes via namelist input
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
MODULE mo_grid_levels
  !-------------------------------------------------------------------------
  !
  !    ProTeX FORTRAN source: Style 2
  !    modified for ICON project, DWD/MPI-M 2005
  !
  !-------------------------------------------------------------------------
  USE mo_kind,               ONLY: wp
  USE mo_exception,          ONLY: message_text, message, finish

  USE mo_io_units,           ONLY: filename_max, nnml, nstat
  USE mo_namelist,           ONLY: position_nml, open_nml, positioned
  USE mo_math_constants,     ONLY: pi, rad2deg
  USE mo_physical_constants, ONLY: earth_radius
  USE mo_math_utilities,     ONLY: t_geographical_coordinates,  &
       &                           cc2gc
  USE mo_base_geometry,      ONLY:  x_rot_angle, y_rot_angle, z_rot_angle

  USE mo_grid

  USE mo_base_datatypes,     ONLY: t_spheres

  USE mo_topology,           ONLY: spheres_on_levels


  USE mo_impl_constants,     ONLY: min_rlcell, max_rlcell, &
       &                           min_rlvert, max_rlvert, &
       &                           min_rledge, max_rledge
  USE mo_util_uuid,          ONLY: t_uuid, uuid_generate, &
       &                           uuid_unparse, uuid_string_length

  USE mo_math_utilities,     ONLY: check_orientation

  IMPLICIT NONE

  PRIVATE

  INCLUDE 'netcdf.inc'


  PUBLIC :: init_grid, destruct_grids, nf
  PUBLIC :: grid_on_level,number_of_grid_levels,     &
    &       l_c_grid, itype_optimize, maxlev_optim,  &
    &       tria_arc_km, beta_spring, init_gridgen,  &
    &       check_orientation, statistics, output_grid

  LOGICAL :: lgrid_initialized = .false.
  INTEGER :: number_of_grid_levels

  LOGICAL :: l_c_grid        ! Use C-grid constraint for last refinement level
  INTEGER :: itype_optimize  ! Type of grid-optimization:
  ! = 0 : natural grid
  ! = 1 : Heikes optimization
  ! = 2 : equal area subdivision
  ! = 3 : c-grid small circle constraint
  ! = 4 : spring dynamics
  INTEGER :: maxlev_optim    ! For itype_optimize=1,4: highest level for which
  ! optimization is executed
  REAL(wp):: tria_arc_km     ! grid distance (km) for the planar grid option
  ! on the finest chosen level
  REAL(wp):: beta_spring     ! for spring dynamics: weighting of a larger or
  ! a smaller target grid length
  ! (choice: 1.11 larger target length: good for
  !          C-grid constraint,
  !          0.9  smaller target length: good for
  !          similar triangle shapes
  INTEGER :: grid_root
  ! meta data descriptors
  CHARACTER(len=*), PARAMETER :: uri_pathname = 'http://icon-downloads.mpimet.mpg.de/grids/public/'
  INTEGER :: number_of_grid_used     ! index for reference to xml table for grib2, 0 for private use
  INTEGER :: centre                  ! centre running the grid generator: 78 - edzw (DWD), 252 - MPIM
  INTEGER :: subcentre               ! subcentre to be assigned by centre, usually 0

  INTEGER :: annotate_level          ! which grid level to annotate
  !--------------------------------------------------------------
  TYPE(t_grid) , ALLOCATABLE, TARGET :: grid_on_level(:)


CONTAINS

  !-------------------------------------------------------------------------
  !>
  !!  Initializes the grid generator  run,.
  !!
  !!  Initializes the grid generator  run,
  !!  reading namelists from the <b>GRID_INI</b> file
  !!  in the directory <b>input</b>.
  !!
  SUBROUTINE init_gridgen(grid_levels,nroot,lplane,lread_graph)
    INTEGER, INTENT(out)::grid_levels, nroot
    LOGICAL, INTENT(out)::lplane, lread_graph

    INTEGER :: i_status, &
      & itype_optimize_d, nroot_d, grid_levels_d, maxlev_optim_d

    LOGICAL :: l_c_grid_d, lplane_d, lread_graph_d

    REAL(wp):: x_rot_angle_d, y_rot_angle_d, z_rot_angle_d, tria_arc_km_d, &
      & beta_spring_d

    !--------------------------------------------------------------------
    !BOC
    !
    ! ! initialize with grid parameters from files
    !
    NAMELIST /grid_ini/     grid_levels, nroot, lplane, lread_graph
    NAMELIST /grid_options/ x_rot_angle, y_rot_angle, z_rot_angle, &
      & itype_optimize, l_c_grid, maxlev_optim, &
      & beta_spring
    NAMELIST /grid_metadata/ number_of_grid_used, centre, subcentre, &
      & annotate_level
    NAMELIST /plane_options/tria_arc_km

    ! set default values for grid_options
    l_c_grid       = .false.
    itype_optimize = 4
    x_rot_angle    = 0.0_wp
    y_rot_angle    = 0.0_wp
    z_rot_angle    = 0.0_wp
    maxlev_optim   = 100
    beta_spring    = 0.90_wp

    ! set default values for grid_ini
    nroot          = 2
    grid_levels    = 4
    lplane         = .false.
    lread_graph    = .false.

    ! set default values for grid_options
    number_of_grid_used = 0
    centre = 65535                  ! missing value from WMO common code table C-11
    subcentre = 0
    annotate_level = grid_levels    ! level to annotate, to allow definition of higher level nests

    ! set default values for plane_options
    tria_arc_km    = 10.0_wp ! resolution in kilometers

    ! copy default values to "default" variables
    nroot_d          = nroot
    grid_levels_d    = grid_levels
    l_c_grid_d       = l_c_grid
    itype_optimize_d = itype_optimize
    maxlev_optim_d   = maxlev_optim
    x_rot_angle_d    = x_rot_angle
    y_rot_angle_d    = y_rot_angle
    z_rot_angle_d    = z_rot_angle
    lplane_d         = lplane
    lread_graph_d    = lread_graph
    tria_arc_km_d    = tria_arc_km
    beta_spring_d    = beta_spring

    ! read namelist
    CALL open_nml('NAMELIST_GRID')

    CALL position_nml('grid_ini', STATUS=i_status)
    IF (i_status == positioned) THEN
      READ (nnml,grid_ini)
    ENDIF

    CALL position_nml('grid_options', STATUS=i_status)
    IF (i_status == positioned) THEN
      READ (nnml,grid_options)
    ENDIF

   CALL position_nml('grid_metadata', STATUS=i_status)
    IF (i_status == positioned) THEN
      READ (nnml,grid_metadata)
    ENDIF

    IF (lplane) THEN
      CALL position_nml('plane_options', STATUS=i_status)
      IF (i_status == positioned) THEN
        READ (nnml,plane_options)
      ENDIF
    ENDIF

    CLOSE(nnml)

    IF(lplane)THEN
      itype_optimize=0
      l_c_grid=.false.
    ENDIF

    grid_root = nroot

    ! write control output of namelist variables

    CALL message('', message_text)
    WRITE(message_text,'(a)')'Namelist Group: grid_ini'
    CALL message('', message_text)
    WRITE(message_text,'(a)')'------------------------'
    CALL message('', message_text)
    WRITE(message_text,'(t7,a,t28,a,t43,a)') 'Variable', 'Actual Value', 'Default Value'
    CALL message('', message_text)
    WRITE(message_text,'(t8,a,t28,i12,t44,i12)') 'grid_levels', grid_levels, grid_levels_d
    CALL message('', message_text)
    WRITE(message_text,'(t8,a,t28,i12,t44,i12)') 'nroot' , nroot, nroot_d
    CALL message('', message_text)
    WRITE(message_text,'(t8,a,t28,l12,t44,l12)') 'lplane', lplane, lplane_d
    CALL message('', message_text)
    WRITE(message_text,'(t8,a,t28,l12,t44,l12)') 'lread_graph', lread_graph, lread_graph_d
    CALL message('', message_text)
    CALL message ('','')
    WRITE(message_text,'(a)')'Namelist Group: grid_options'
    CALL message('', message_text)
    WRITE(message_text,'(a)')'----------------------------'
    CALL message('', message_text)
    WRITE(message_text,'(t7,a,t28,a,t43,a)')'Variable', 'Actual Value', 'Default Value'
    CALL message('', message_text)
    WRITE(message_text,'(t8,a,t28,i12,t44,i12)')'itype_optimize', itype_optimize, itype_optimize_d
    CALL message('', message_text)
    WRITE(message_text,'(t8,a,t28,i12,t44,i12)')'maxlev_optim', maxlev_optim, maxlev_optim_d
    CALL message('', message_text)
    WRITE(message_text,'(t8,a,t28,l12,t44,l12)')'l_c_grid', l_c_grid, l_c_grid_d
    CALL message('', message_text)
    WRITE(message_text,'(t8,a,t28,f12.4,t44,f12.4)')'beta_spring', beta_spring, beta_spring_d
    CALL message('', message_text)
    WRITE(message_text,'(t8,a,t28,f12.4,t44,f12.4)')'x_rot_angle', x_rot_angle, x_rot_angle_d
    CALL message('', message_text)
    WRITE(message_text,'(t8,a,t28,f12.4,t44,f12.4)')'y_rot_angle', y_rot_angle, y_rot_angle_d
    CALL message('', message_text)
    WRITE(message_text,'(t8,a,t28,f12.4,t44,f12.4)')'z_rot_angle', z_rot_angle, z_rot_angle_d
    CALL message('', message_text)
    CALL message ('','')
    WRITE(message_text,'(a)')'Namelist Group: plane_options'
    CALL message('', message_text)
    WRITE(message_text,'(a)')'-----------------------------'
    CALL message('', message_text)
    WRITE(message_text,'(t8,a,t28,f12.4,t44,f12.4)')'tria_arc_km', tria_arc_km, tria_arc_km_d
    CALL message('', message_text)
    CALL message ('','')

    ! convert rotation angles do radians

    x_rot_angle= x_rot_angle/rad2deg
    y_rot_angle= y_rot_angle/rad2deg
    z_rot_angle= z_rot_angle/rad2deg

  END SUBROUTINE init_gridgen
  !-------------------------------------------------------------------------
  !>
  !!
  !! @par Revision History
  !! Luis Kornblueh, MPI-M, Hamburg, April 2013
  !! dump the xml table entry for grib2 handling of horizontal grids
  !!
  SUBROUTINE dump_grid_table_xml_metadata(filename)

    CHARACTER(len=*), INTENT(in) :: filename
    CHARACTER(len=255) :: uriname, uri_subcentre

    CALL message('','')
    CALL message('','---------- BEGIN XML grid table descriptor ----------')
    CALL message('','')

    ! grouping element: grid
    WRITE(message_text,'(a,i0,a,i0,a,i0,a)') &
          &          '<grid number_of_grid_used="', number_of_grid_used, &
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

  !-------------------------------------------------------------------------
  !>
  !!
  !! @par Revision History
  !!  Luis Kornblueh, MPI-M, Hamburg, March 2004
  !!  Luca Bonaventura, MPI-M, Hamburg, March 2005
  !!  Dual cell orientation and edge area computation
  !!  added by Luca Bonaventura, MPI-M, Hamburg, December 2005
  !!
  SUBROUTINE init_grid(lplane)
    !

    LOGICAL, INTENT(in)    :: lplane

    INTEGER :: j, jlev, ie, iee, iv1, iv2
    INTEGER, ALLOCATABLE :: nots(:),noeds(:),novs(:)
    TYPE(t_grid), POINTER :: cg
    TYPE(t_spheres), POINTER :: ct
    TYPE(t_geographical_coordinates) :: geog

    !-------------------------------------------------------------------------

    ! initialization

    number_of_grid_levels = UBOUND(spheres_on_levels,1)

    ALLOCATE(nots(0:number_of_grid_levels))
    ALLOCATE(noeds(0:number_of_grid_levels))
    ALLOCATE(novs(0:number_of_grid_levels))

    DO jlev = 0, number_of_grid_levels
      nots(jlev)   = spheres_on_levels(jlev)%no_triangles
      noeds(jlev)  = spheres_on_levels(jlev)%no_edges
      novs(jlev)   = spheres_on_levels(jlev)%no_vertices
    ENDDO

    CALL construct_grids (nots, noeds, novs)

    DO jlev = 0, number_of_grid_levels

      WRITE(message_text,'(a,i0)') 'Mapping data structure on level ', jlev
      CALL message('', message_text)

      cg => grid_on_level(jlev)
      ct => spheres_on_levels(jlev)

      cg%ncells=nots(jlev)
      cg%nverts=novs(jlev)
      cg%nedges=noeds(jlev)
      DO j = 0, nots(jlev)-1
        cg%cells%idx(j+1)               = ct%ts(j)%info%idx+1
        IF (jlev>0) THEN
          cg%cells%parent_index(j+1)    = ct%ts(j)%parent%info%idx+1
        ELSE
          cg%cells%parent_index(j+1)    = ct%ts(j)%info%parent_number+1
        ENDIF
        IF (jlev == number_of_grid_levels) THEN
          cg%cells%child_index(j+1,1)  = -1
          cg%cells%child_index(j+1,2)  = -1
          cg%cells%child_index(j+1,3)  = -1
          cg%cells%child_index(j+1,4)  = -1
        ELSE
          cg%cells%child_index(j+1,1)  = ct%ts(j)%sub_triangle0%info%idx+1
          cg%cells%child_index(j+1,2)  = ct%ts(j)%sub_triangle1%info%idx+1
          cg%cells%child_index(j+1,3)  = ct%ts(j)%sub_triangle2%info%idx+1
          cg%cells%child_index(j+1,4)  = ct%ts(j)%sub_triangle3%info%idx+1
        ENDIF
        cg%cells%neighbor_index(j+1,1)  = ct%ts(j)%neighbor0%info%idx+1
        cg%cells%neighbor_index(j+1,2)  = ct%ts(j)%neighbor1%info%idx+1
        cg%cells%neighbor_index(j+1,3)  = ct%ts(j)%neighbor2%info%idx+1
        cg%cells%edge_index(j+1,1)      = ct%ts(j)%edge0%info%idx+1
        cg%cells%edge_index(j+1,2)      = ct%ts(j)%edge1%info%idx+1
        cg%cells%edge_index(j+1,3)      = ct%ts(j)%edge2%info%idx+1
        cg%cells%edge_orientation(j+1,1)= ct%ts(j)%edge0_orientation
        cg%cells%edge_orientation(j+1,2)= ct%ts(j)%edge1_orientation
        cg%cells%edge_orientation(j+1,3)= ct%ts(j)%edge2_orientation
        cg%cells%vertex_index(j+1,1)    = ct%ts(j)%vertex0%info%idx+1
        cg%cells%vertex_index(j+1,2)    = ct%ts(j)%vertex1%info%idx+1
        cg%cells%vertex_index(j+1,3)    = ct%ts(j)%vertex2%info%idx+1
        IF (.not.lplane) THEN
          geog                         = cc2gc(ct%ts(j)%triangle_center)
          cg%cells%center(j+1)%lon     = geog%lon
          cg%cells%center(j+1)%lat     = geog%lat
        ELSE
          cg%cells%center(j+1)%lon     = ct%ts(j)%triangle_center%x(1)
          cg%cells%center(j+1)%lat     = ct%ts(j)%triangle_center%x(2)
        ENDIF
        cg%cells%area(j+1)              = ct%ts(j)%triangle_area
        ! Information on mesh refinement becomes only available in patchgen
        cg%cells%refin_ctrl(j+1)        = 0
        cg%cells%child_id(j+1)          = 0
      ENDDO

      DO j = 0, noeds(jlev)-1
        cg%edges%idx(j+1)               = ct%es(j)%info%idx+1
        cg%edges%cell_index(j+1,1)      = ct%es(j)%triangle0%info%idx+1
        cg%edges%cell_index(j+1,2)      = ct%es(j)%triangle1%info%idx+1
        cg%edges%vertex_index(j+1,1)    = ct%es(j)%vertex0%info%idx+1
        cg%edges%vertex_index(j+1,2)    = ct%es(j)%vertex1%info%idx+1
        cg%edges%system_orientation(j+1)= ct%es(j)%system_orientation
        IF (.not.lplane) THEN
          geog                         = cc2gc(ct%es(j)%edge_center)
          cg%edges%center(j+1)%lon     = geog%lon
          cg%edges%center(j+1)%lat     = geog%lat
        ELSE
          cg%edges%center(j+1)%lon     = ct%es(j)%edge_center%x(1)
          cg%edges%center(j+1)%lat     = ct%es(j)%edge_center%x(2)
        ENDIF
        cg%edges%primal_normal(j+1)     = ct%es(j)%edge_primal_normal
        cg%edges%dual_normal(j+1)       = ct%es(j)%edge_dual_normal
        cg%edges%primal_edge_length(j+1)= ct%es(j)%edge_primal_arc
        cg%edges%dual_edge_length(j+1)  = ct%es(j)%edge_dual_arc
        cg%edges%edge_vert_length(j+1,1)= ct%es(j)%edge_vert0_arc
        cg%edges%edge_vert_length(j+1,2)= ct%es(j)%edge_vert1_arc
        cg%edges%edge_cell_length(j+1,1)= ct%es(j)%edge_cell0_arc
        cg%edges%edge_cell_length(j+1,2)= ct%es(j)%edge_cell1_arc
        ! Information on mesh refinement becomes only available in patchgen
        cg%edges%refin_ctrl(j+1)        = 0
        ! Also, parent and child indices for edges are only set in patchgen
        cg%edges%parent_index(j+1)      = -1
        cg%edges%child_index(j+1,1:4)   = 0
        cg%edges%child_id(j+1)          = 0
      ENDDO

      DO j = 0, novs(jlev)-1
        cg%verts%idx(j+1)               = ct%vs(j)%info%idx+1
        cg%verts%neighbor_index(j+1,1)  = ct%vs(j)%neighbor0%info%idx+1
        cg%verts%neighbor_index(j+1,2)  = ct%vs(j)%neighbor1%info%idx+1
        cg%verts%neighbor_index(j+1,3)  = ct%vs(j)%neighbor2%info%idx+1
        cg%verts%neighbor_index(j+1,4)  = ct%vs(j)%neighbor3%info%idx+1
        cg%verts%neighbor_index(j+1,5)  = ct%vs(j)%neighbor4%info%idx+1
        cg%verts%neighbor_index(j+1,6)  = ct%vs(j)%neighbor5%info%idx+1
        cg%verts%cell_index(j+1,1)      = ct%vs(j)%triangle0%info%idx+1
        cg%verts%cell_index(j+1,2)      = ct%vs(j)%triangle1%info%idx+1
        cg%verts%cell_index(j+1,3)      = ct%vs(j)%triangle2%info%idx+1
        cg%verts%cell_index(j+1,4)      = ct%vs(j)%triangle3%info%idx+1
        cg%verts%cell_index(j+1,5)      = ct%vs(j)%triangle4%info%idx+1
        cg%verts%cell_index(j+1,6)      = ct%vs(j)%triangle5%info%idx+1
        cg%verts%edge_index(j+1,1)      = ct%vs(j)%edge0%info%idx+1
        cg%verts%edge_index(j+1,2)      = ct%vs(j)%edge1%info%idx+1
        cg%verts%edge_index(j+1,3)      = ct%vs(j)%edge2%info%idx+1
        cg%verts%edge_index(j+1,4)      = ct%vs(j)%edge3%info%idx+1
        cg%verts%edge_index(j+1,5)      = ct%vs(j)%edge4%info%idx+1
        cg%verts%edge_index(j+1,6)      = ct%vs(j)%edge5%info%idx+1
        cg%verts%dual_area(j+1)         = ct%vs(j)%area_dual_cell
        IF (.not.lplane) THEN
          geog                         = cc2gc(ct%vs(j)%vertex)
          cg%verts%vertex(j+1)%lon     = geog%lon
          cg%verts%vertex(j+1)%lat     = geog%lat
        ELSE
          cg%verts%vertex(j+1)%lon     = ct%vs(j)%vertex%x(1)
          cg%verts%vertex(j+1)%lat     = ct%vs(j)%vertex%x(2)
        ENDIF
        ! Information on mesh refinement becomes only available in patchgen
        cg%verts%refin_ctrl(j+1)        = 0
        cg%verts%child_id(j+1)          = 0
      ENDDO

      !
      !       compute dual edge orientations
      !
      DO j = 1, novs(jlev)

        DO ie=1,6

          iee = cg%verts%edge_index(j,ie)

          IF(iee==0) CYCLE

          iv1   = cg%edges%vertex_index(iee,1)
          iv2   = cg%edges%vertex_index(iee,2)

          IF(j==iv1) THEN

            IF(iv2 > iv1)THEN
              cg%verts%edge_orientation(j,ie)= cg%edges%system_orientation(iee)
            ELSEIF( iv2 < iv1)THEN
              cg%verts%edge_orientation(j,ie)=-cg%edges%system_orientation(iee)
            ENDIF

          ELSEIF(j==iv2)THEN

            IF(iv2 < iv1)THEN
              cg%verts%edge_orientation(j,ie)= cg%edges%system_orientation(iee)
            ELSEIF(iv2 > iv1)THEN
              cg%verts%edge_orientation(j,ie)=-cg%edges%system_orientation(iee)
            ENDIF

          ENDIF

        ENDDO

      ENDDO

      !
      ! Set dummy index list arrays
      !
      cg%verts%indlist(0:max_rlvert ,1) = 1
      cg%verts%indlist(min_rlvert:-1,1) = novs(jlev)+1
      cg%verts%indlist(1:max_rlvert ,2) = 0
      cg%verts%indlist(min_rlvert:0 ,2) = novs(jlev)

      cg%edges%indlist(0:max_rledge ,1) = 1
      cg%edges%indlist(min_rledge:-1,1) = noeds(jlev)+1
      cg%edges%indlist(1:max_rledge ,2) = 0
      cg%edges%indlist(min_rledge:0 ,2) = noeds(jlev)

      cg%cells%indlist(0:max_rlcell ,1) = 1
      cg%cells%indlist(min_rlcell:-1,1) = nots(jlev)+1
      cg%cells%indlist(1:max_rlcell ,2) = 0
      cg%cells%indlist(min_rlcell:0 ,2) = nots(jlev)

      cg%verts%start_idx(0:max_rlvert ,1) = 1
      cg%verts%start_idx(min_rlvert:-1,1) = novs(jlev)+1
      cg%verts%end_idx(1:max_rlvert ,1) = 0
      cg%verts%end_idx(min_rlvert:0 ,1) = novs(jlev)

      cg%edges%start_idx(0:max_rledge ,1) = 1
      cg%edges%start_idx(min_rledge:-1,1) = noeds(jlev)+1
      cg%edges%end_idx(1:max_rledge ,1) = 0
      cg%edges%end_idx(min_rledge:0 ,1) = noeds(jlev)

      cg%cells%start_idx(0:max_rlcell ,1) = 1
      cg%cells%start_idx(min_rlcell:-1,1) = nots(jlev)+1
      cg%cells%end_idx(1:max_rlcell ,1) = 0
      cg%cells%end_idx(min_rlcell:0 ,1) = nots(jlev)

    ENDDO

    DEALLOCATE (nots)
    DEALLOCATE (noeds)
    DEALLOCATE (novs)

    lgrid_initialized = .true.

  END SUBROUTINE init_grid
  !------------------------------------------------------------------------------

  !>
  !! Allocates arrays for complete grid hierarchy.
  !!
  !!
  !! @par Revision History
  !! Developed  by  Luis Kornblueh, MPI-M  (2005).
  !! Some adjustment by Luca Bonaventuta, MPI-M (2005).
  !! Modification by Thomas Heinze, DWD (2006-10-31):
  !! - replaced start_level by k_start_level in definition part
  !!
  SUBROUTINE construct_grids (k_ce,k_ed,k_ve)

    INTEGER, INTENT(in) :: k_ce(0:number_of_grid_levels), &
      & k_ed(0:number_of_grid_levels), &
      & k_ve(0:number_of_grid_levels)

    INTEGER :: j,istat

    !-----------------------------------------------------------------------

    ALLOCATE (grid_on_level(0:number_of_grid_levels),stat=istat)
    IF (istat >0) THEN
      CALL finish ('construct_grids', 'Problem in allocating grid_on_level')
    ENDIF

    DO j = 0, number_of_grid_levels

      CALL construct_grid(grid_on_level(j),k_ce(j),k_ed(j),k_ve(j))

    ENDDO

  END SUBROUTINE construct_grids
  !------------------------------------------------------------------------------

  !------------------------------------------------------------------------------
  !>
  !! Deallocates arrays for complete grid hierarchy.
  !!
  !!
  !! @par Revision History
  !! Developed  by  Luis Kornblueh, MPI-M (2005).
  !! calls to deallocation subroutines by Almut Gassmann, MPI-M (2007-04)
  !!
  SUBROUTINE destruct_grids
    !-----------------------------------------------------------------------
    INTEGER :: j,istat

    IF (ALLOCATED (grid_on_level)) THEN

      DO j = 0, number_of_grid_levels

        CALL destruct_grid(grid_on_level(j))

      ENDDO

      DEALLOCATE (grid_on_level,stat=istat)
      IF (istat /= 0 ) THEN
        WRITE (message_text,'(a,i2)') &
          & 'Problem in destructing grid_on_levels ', j
        CALL finish ('destruct_grids', message_text)
      ENDIF

    END IF

  END SUBROUTINE destruct_grids
  !-------------------------------------------------------------------------
  ! !IROUTINE: output_grid
  !
  ! !SUBROUTINE INTERFACE:
  SUBROUTINE output_grid(lplane)
    !
    ! !DESCRIPTION:
    ! Creates files  with  cross reference tables and coordinates
    ! of grid items to be read by the model. In most cases for all grid levels,
    ! the file {\bf GRIDMAP.*} is produced, where * denotes the grid level.
    ! The exception is Heikes-Randall + C grid optimization in the last step. Here
    ! only the {\bf GRIDMAP} file of the last level is created. The background
    ! is in this case the {\bf GRIDMAP} files of all the coarser levels would
    ! produce just Heikes-Randall optimized grid without the C grid optimization.
    ! Only in this case gridgen has to be run for all grid levels.\\
    ! The files are produced as sequential access, unformatted output.
    !
    ! !REVISION HISTORY:
    ! Developed and tested  by L.Bonaventura  and Th. Heinze (2004-5).
    ! Modified by Th. Heinze (DWD) (2007-10-31):
    ! - introduced seperate treatment of Heikes-Randall + C grid optimization
    !
    LOGICAL, INTENT(in) :: lplane
    !
    ! !LOCAL VARIABLES:
    TYPE(t_grid), POINTER :: gg     ! pointer to the grid
    TYPE(t_spheres), POINTER :: sphere_level

    CHARACTER(LEN=filename_max) :: output

    INTEGER :: i_loop_start       ! start index of loop
    INTEGER :: i_nc, i_ne, i_nv   ! number of cells, edges and vertices
    INTEGER :: istat              ! status of output
    INTEGER :: jgrid              ! loop index

    INTEGER :: old_mode

    INTEGER :: ncid
    INTEGER :: dimids(2)

    INTEGER :: dim_ncell, dim_nvertex, dim_nedge, dim_two, dim_cell_refine, &
      & dim_edge_refine, dim_vert_refine, dim_nvertex_per_cell,      &
      & dim_ncells_per_edge, dim_nedges_per_vertex, dim_nchilds_per_cell, &
      & dim_list, dim_nchdom

    INTEGER :: varid_clon, varid_clat, varid_clonv, varid_clatv
    INTEGER :: varid_vlon, varid_vlat, varid_vlonv, varid_vlatv
    INTEGER :: varid_vx, varid_vy, varid_vz
    INTEGER :: varid_elon, varid_elat, varid_elonv, varid_elatv

    INTEGER :: varid_carea, varid_varea

    INTEGER :: varid1, varid2, varid3, varid4, varid5, varid6, varid7, varid8, &
      & varid9, varid10, varid11, varid12, varid13, varid14, varid15,   &
      & varid16, varid17, varid18, varid19, varid20, varid21, varid22,  &
      & varid23, varid24, varid25, varid26, varid27, varid28, varid29,  &
      & varid30, varid31, varid32, varid33, varid34, varid35, varid36,  &
      & varid37, varid38, varid39, varid40, varid41, varid42, varid43,  &
      & varid44, varid45, varid46, varid47, varid48, varid251,          &
      & varid282, varid403, varid49, varid50

    INTEGER :: i, j
    INTEGER :: str_idx, end_idx

#ifdef __SX__
    INTEGER :: iargc
    CHARACTER(LEN= 32) :: arg_str
#endif
    CHARACTER(LEN=256) :: command_line
    INTEGER :: command_line_len

    REAL(wp), ALLOCATABLE :: zv2d(:,:), zv2dx(:,:), zv2dy(:,:)
    
    TYPE(t_uuid) :: uuid
    CHARACTER(len=uuid_string_length) :: uuid_string

    REAL(wp) :: rotation_vector(3)

    REAL(wp) :: swap(4)

    CHARACTER(LEN=7) :: optfix

    !EOP
    !-------------------------------------------------------------------------

    !-------------------------------------------------------------------------
    ! get unique grid file identifier for GRIB2 and updated CF-Convention
    CALL uuid_generate(uuid)
    CALL uuid_unparse(uuid, uuid_string)

    ! distinguish between gridgeneration for optimization strategies

    IF ( l_c_grid) THEN

      CALL message ('', '')
      WRITE(message_text,'(a,a)')  '!!! --- REMINDER --- !!!'
      CALL message ('', message_text)
      CALL message ('', '')
      WRITE(message_text,'(a)') 'Rerun gridgen to create GRIDMAP files on other levels!'
      CALL message ('', message_text)
      CALL message ('', '')

      i_loop_start = number_of_grid_levels

    ELSE

      i_loop_start = 0
      i_loop_start = 0

    ENDIF

    DO jgrid = i_loop_start, number_of_grid_levels

      IF (jgrid /= number_of_grid_levels .and. l_c_grid ) CYCLE

      SELECT CASE (itype_optimize)
      CASE (0)
        ! natural grid
        optfix = 'noo'
      CASE (1)
        IF (l_c_grid) THEN
          ! Heikes-Randall with additional c-grid optimization
          optfix = 'hrc'
        ELSE
          ! Heikes-Randall optimization
          optfix = 'hro'
        ENDIF
      CASE (2)
        ! equal area subdivision
        IF (l_c_grid) THEN
          optfix = 'eac'
        ELSE
          optfix = 'eas'
        ENDIF
      CASE (3)
        ! c-grid small circle constraint
        optfix = 'scc'
      CASE (4)
        ! spring dynamics
        IF (l_c_grid) THEN
          optfix = 'spc'
        ELSE
          WRITE(optfix,'(a3,f4.2)') 'spr',beta_spring
        ENDIF
      END SELECT

      IF(lplane) THEN
        WRITE (output,'(a,i0,a,i2.2,a,a,a)') &
          & 'planR', grid_root, 'B', jgrid , '-grid.nc'
      ELSE
        IF (jgrid >= maxlev_optim) THEN
          WRITE (output,'(a,i0,a,i2.2,3a,i1,a)') &
            & 'iconR', grid_root, 'B', jgrid , '-grid_', TRIM(optfix), '_M', maxlev_optim, '.nc'
        ELSE
          WRITE (output,'(a,i0,a,i2.2,a,a,a)') &
            & 'iconR', grid_root, 'B', jgrid , '-grid_', TRIM(optfix), '.nc'
        ENDIF
      ENDIF

      WRITE(message_text,'(a,a)') 'Write gridmap file: ', TRIM(output)
      CALL message ('', message_text)

      gg => grid_on_level(jgrid)
      sphere_level => spheres_on_levels(jgrid)

      i_nc = gg%ncells
      i_ne = gg%nedges
      i_nv = gg%nverts

#ifdef __SX__
      command_line = ''
      DO i = 0, iargc()
        CALL getarg(i, arg_str)
        command_line = TRIM(command_line)//TRIM(arg_str)
      ENDDO
      command_line_len = LEN_TRIM(command_line)
#else
      CALL GET_COMMAND(command_line, command_line_len, istat)
#endif
      rotation_vector = (/ x_rot_angle, y_rot_angle, z_rot_angle /)

      !----------------------------------------------------------------------
      !
      CALL nf(nf_set_default_format(nf_format_64bit, old_mode))
      CALL nf(nf_create(TRIM(output), nf_clobber, ncid))
      CALL nf(nf_set_fill(ncid, nf_nofill, old_mode))
      !
      !----------------------------------------------------------------------
      !
      ! Global attributes
      !
      CALL nf(nf_put_att_text  (ncid, nf_global, 'title', 21, 'ICON grid description'))
      CALL nf(nf_put_att_text  (ncid, nf_global, 'history', command_line_len, command_line))
      CALL nf(nf_put_att_text  (ncid, nf_global, 'institution', 59, &
        & 'Max Planck Institute for Meteorology/Deutscher Wetterdienst'))
      CALL nf(nf_put_att_text  (ncid, nf_global, 'source', 10, 'icon-dev'))
      !
      CALL nf(nf_put_att_text  (ncid, nf_global, 'uuidOfHGrid' , uuid_string_length, TRIM(uuid_string)))
      IF (jgrid == annotate_level) THEN
        IF (number_of_grid_used == 0) THEN
          CALL message('','number_of_grid_used is 0 and cannot be added to the ICON master grid table')
        ENDIF
        CALL nf(nf_put_att_int (ncid, nf_global, 'number_of_grid_used', nf_int, 1, number_of_grid_used))
        CALL nf(nf_put_att_text  (ncid, nf_global, 'ICON_grid_file_uri' , &
            &                                       LEN_TRIM(uri_pathname//output), TRIM(uri_pathname//output)))
      ELSE
        CALL nf(nf_put_att_int (ncid, nf_global, 'number_of_grid_used', nf_int, 1, 0))
        CALL nf(nf_put_att_text  (ncid, nf_global, 'ICON_grid_file_uri' , 7, 'private'))
      ENDIF
      CALL nf(nf_put_att_int   (ncid, nf_global, 'centre', nf_int, 1, centre))
      CALL nf(nf_put_att_int   (ncid, nf_global, 'subcentre', nf_int, 1, subcentre))
      CALL nf(nf_put_att_text  (ncid, nf_global, 'grid_mapping_name' , 18, 'lat_long_on_sphere'))
      CALL nf(nf_put_att_text  (ncid, nf_global, 'crs_id' , 28, 'urn:ogc:def:cs:EPSG:6.0:6422'))
      CALL nf(nf_put_att_text  (ncid, nf_global, 'crs_name',30,'Spherical 2D Coordinate System'))
      CALL nf(nf_put_att_text  (ncid, nf_global, 'ellipsoid_name' , 6, 'Sphere'))
      CALL nf(nf_put_att_double(ncid, nf_global, 'semi_major_axis' , nf_double, 1, earth_radius))
      CALL nf(nf_put_att_double(ncid, nf_global, 'inverse_flattening' , nf_double, 1, 0.0_wp))
      CALL nf(nf_put_att_int   (ncid, nf_global, 'grid_level', nf_int, 1, jgrid))
      CALL nf(nf_put_att_int   (ncid, nf_global, 'grid_root', nf_int, 1, grid_root))
      ! The following three attributes are nontrivial in the presence of grid refinement
      ! and are set here because they are checked in the ICON code
      CALL nf(nf_put_att_int     (ncid, nf_global, 'grid_ID', nf_int, 1, 1))
      CALL nf(nf_put_att_int     (ncid, nf_global, 'parent_grid_ID', nf_int, 1, 0))
      CALL nf(nf_put_att_int     (ncid, nf_global, 'max_childdom', nf_int, 1, 1))
      SELECT CASE (itype_optimize)
      CASE (0)
        CALL nf(nf_put_att_text  (ncid, nf_global, 'grid_optimization', 12, 'natural grid'))
      CASE (1)
        IF (l_c_grid) THEN
          CALL nf(nf_put_att_text(ncid, nf_global, 'grid_optimization', 50, &
            & 'Heikes-Randall with additional c-grid optimization'))
        ELSE
          CALL nf(nf_put_att_text(ncid, nf_global, 'grid_optimization', 27, &
            & 'Heikes-Randall optimization'))
        ENDIF
      CASE (2)
        CALL nf(nf_put_att_text  (ncid, nf_global, 'grid_optimization', 22, &
          & 'equal area subdivision'))
      CASE (3)
        CALL nf(nf_put_att_text  (ncid, nf_global, 'grid_optimization', 30, &
          & 'c-grid small circle constraint'))
      END SELECT
      CALL nf(nf_put_att_double  (ncid, nf_global, 'rotation_vector', nf_double,3,&
        & rotation_vector))
      !
      ! Dimensions
      !
      CALL nf(nf_def_dim(ncid, 'cell',   i_nc, dim_ncell))
      CALL nf(nf_def_dim(ncid, 'vertex', i_nv, dim_nvertex))
      CALL nf(nf_def_dim(ncid, 'edge',   i_ne, dim_nedge))
      !
      CALL nf(nf_def_dim(ncid, 'nc',        2, dim_ncells_per_edge))
      CALL nf(nf_def_dim(ncid, 'nv',        3, dim_nvertex_per_cell))
      CALL nf(nf_def_dim(ncid, 'ne',        6, dim_nedges_per_vertex))
      !
      CALL nf(nf_def_dim(ncid, 'no',        4, dim_nchilds_per_cell))

      ! Dimensions for refinement
      CALL nf(nf_def_dim(ncid, 'two_grf',   2,    dim_two))
      CALL nf(nf_def_dim(ncid, 'max_chdom', 1, dim_nchdom))
      dim_list = max_rlcell-min_rlcell+1
      CALL nf(nf_def_dim(ncid, 'cell_grf',dim_list, dim_cell_refine))
      dim_list = max_rledge-min_rledge+1
      CALL nf(nf_def_dim(ncid, 'edge_grf',dim_list, dim_edge_refine))
      dim_list = max_rlvert-min_rlvert+1
      CALL nf(nf_def_dim(ncid, 'vert_grf',dim_list, dim_vert_refine))
      !
      !---------------------------------------------------------------------
      !
      ! Grid variables
      !
      !---------------------------------------------------------------------
      !
      ! public grid information
      !
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
      CALL nf(nf_def_var(ncid, 'cartesian_x_vertices', nf_double, 1, dim_nvertex, varid_vx))
      CALL nf(nf_put_att_text(ncid, varid_vx, 'long_name', 40, &
        & 'vertex cartesian coordinate x on unit sphere'))
      CALL nf(nf_put_att_text(ncid, varid_vx, 'units', 6, 'meters'))
      !
      CALL nf(nf_def_var(ncid, 'cartesian_y_vertices', nf_double, 1, dim_nvertex, varid_vy))
      CALL nf(nf_put_att_text(ncid, varid_vy, 'long_name', 40, &
        & 'vertex cartesian coordinate y on unit sphere'))
      CALL nf(nf_put_att_text(ncid, varid_vy, 'units', 6, 'meters'))
      !
      CALL nf(nf_def_var(ncid, 'cartesian_z_vertices', nf_double, 1, dim_nvertex, varid_vz))
      CALL nf(nf_put_att_text(ncid, varid_vz, 'long_name', 40, &
        & 'vertex cartesian coordinate x on unit sphere'))
      CALL nf(nf_put_att_text(ncid, varid_vz, 'units', 6, 'meters'))
      !
      !
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
      ! test variables (areas)
      !
      CALL nf(nf_def_var(ncid, 'cell_area', nf_double, 1, dim_ncell, varid_carea))
      CALL nf(nf_put_att_text(ncid, varid_carea, 'long_name', 17, 'area of grid cell'))
      CALL nf(nf_put_att_text(ncid, varid_carea, 'units', 2, 'm2'))
      CALL nf(nf_put_att_text(ncid, varid_carea, 'standard_name', 4, 'area'))
      CALL nf(nf_put_att_text(ncid, varid_carea, 'coordinates', 9, 'clon clat'))
      !
      CALL nf(nf_def_var(ncid, 'dual_area', nf_double, 1, dim_nvertex, varid_varea))
      CALL nf(nf_put_att_text(ncid, varid_varea, 'long_name', 40, &
        & 'areas of dual hexagonal/pentagonal cells'))
      CALL nf(nf_put_att_text(ncid, varid_varea, 'units', 2, 'm2'))
      CALL nf(nf_put_att_text(ncid, varid_varea, 'standard_name', 4, 'area'))
      CALL nf(nf_put_att_text(ncid, varid_varea, 'coordinates', 9, 'vlon vlat'))
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
      CALL nf(nf_put_att_text(ncid, varid10, 'long_name',35,&
        & 'vertices at the end of of each edge'))
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
      CALL nf(nf_put_att_text(ncid, varid14, 'units', 2, 'm2'))
      CALL nf(nf_put_att_text(ncid, varid14, 'cdi', 6, 'ignore'))
      !
      CALL nf(nf_def_var(ncid, 'dual_area_p', nf_double, 1, dim_nvertex, varid15))
      CALL nf(nf_put_att_text(ncid, varid15, 'long_name', 40, &
        & 'areas of dual hexagonal/pentagonal cells'))
      CALL nf(nf_put_att_text(ncid, varid15, 'units', 2, 'm2'))
      CALL nf(nf_put_att_text(ncid, varid15, 'cdi', 6, 'ignore'))
      !
      CALL nf(nf_def_var(ncid, 'edge_length', nf_double, 1, dim_nedge, varid16))
      CALL nf(nf_put_att_text(ncid, varid16, 'long_name', 36, &
        & 'lengths of edges of triangular cells'))
      CALL nf(nf_put_att_text(ncid, varid16, 'units', 1, 'm'))
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
      CALL nf(nf_put_att_text(ncid, varid17, 'long_name', 71, &
        & 'lengths of dual edges (distances between triangular cell circumcenters)'))
      CALL nf(nf_put_att_text(ncid, varid17, 'units', 1, 'm'))
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
      CALL nf(nf_def_var(ncid, 'meridional_normal_primal_edge', nf_double, 1, &
        & dim_nedge, varid20))
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
      CALL nf(nf_def_var(ncid, 'parent_cell_index', nf_int, 1, dim_ncell, varid25))
      CALL nf(nf_put_att_text(ncid, varid25, 'long_name', 17, 'parent cell index'))
      CALL nf(nf_put_att_text(ncid, varid25, 'cdi', 6, 'ignore'))
      !
      CALL nf(nf_def_var(ncid, 'parent_cell_type', nf_int, 1, dim_ncell, varid251))
      CALL nf(nf_put_att_text(ncid, varid251, 'long_name', 17, 'parent cell type'))
      CALL nf(nf_put_att_text(ncid, varid251, 'cdi', 6, 'ignore'))
      !
      dimids = (/ dim_ncell, dim_nvertex_per_cell /)
      CALL nf(nf_def_var(ncid, 'neighbor_cell_index', nf_int, 2, dimids, varid26))
      CALL nf(nf_put_att_text(ncid, varid26, 'long_name', 19, 'cell neighbor index'))
      CALL nf(nf_put_att_text(ncid, varid26, 'cdi', 6, 'ignore'))
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
      CALL nf(nf_def_var(ncid, 'edge_index', nf_int, 1, dim_nedge, varid28))
      CALL nf(nf_put_att_text(ncid, varid28, 'long_name', 10, 'edge index'))
      CALL nf(nf_put_att_text(ncid, varid28, 'cdi', 6, 'ignore'))
      !
      !       CALL nf(nf_def_var(ncid, 'edge_parent', NF_INT, 1, dim_nedge, varid281))
      !       CALL nf(nf_put_att_text(ncid, varid281, 'long_name', 10, 'edge parent'))
      !       CALL nf(nf_put_att_text(ncid, varid281, 'cdi', 6, 'ignore'))
      !
      CALL nf(nf_def_var(ncid, 'edge_parent_type', nf_int, 1, dim_nedge, varid282))
      CALL nf(nf_put_att_text(ncid, varid282, 'long_name', 10, 'edge parent type'))
      CALL nf(nf_put_att_text(ncid, varid282, 'cdi', 6, 'ignore'))
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
      !
      ! Variables added for mesh refinement
      !
      CALL nf(nf_def_var(ncid, 'refin_c_ctrl', nf_int, 1, dim_ncell, varid32))
      CALL nf(nf_put_att_text(ncid, varid32, 'long_name', 33, &
        & 'refinement control flag for cells'))
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
      CALL nf(nf_put_att_text(ncid, varid34, 'long_name', 33, &
        & 'refinement control flag for edges'))
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
      CALL nf(nf_put_att_text(ncid, varid36, 'long_name', 36, &
        & 'refinement control flag for vertices'))
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
      !
      CALL nf(nf_def_var(ncid, 'parent_vertex_index', nf_int, 1, dim_nvertex, varid403))
      CALL nf(nf_put_att_text(ncid, varid403, 'long_name', 19, 'parent vertex index'))
      CALL nf(nf_put_att_text(ncid, varid403, 'cdi', 6, 'ignore'))
      !
    !
    CALL nf(nf_def_var(ncid, 'phys_cell_id', nf_int, 1, dim_ncell, varid49))
    CALL nf(nf_put_att_text(ncid, varid49, 'long_name', 26, 'physical domain ID of cell'))
    CALL nf(nf_put_att_text(ncid, varid49, 'cdi', 6, 'ignore'))
    !
    CALL nf(nf_def_var(ncid, 'phys_edge_id', nf_int, 1, dim_nedge, varid50))
    CALL nf(nf_put_att_text(ncid, varid50, 'long_name', 26, 'physical domain ID of edge'))
    CALL nf(nf_put_att_text(ncid, varid50, 'cdi', 6, 'ignore'))

      CALL nf(nf_enddef(ncid))
      !
      !-------------------------------------------------------------------------
      !
      ! cell part:
      !
      ! Transpose of index array necessary for CF-1.1 Convention
      !
      ALLOCATE(zv2dx(3, i_nc), zv2dy(3, i_nc))
      DO j = 1, 3
        DO i = 1, i_nc
          zv2dx(j,i) = gg%verts%vertex(gg%cells%vertex_index(i,j))%lon
        ENDDO
      ENDDO
      DO j = 1, 3
        DO i = 1, i_nc
          zv2dy(j,i) = gg%verts%vertex(gg%cells%vertex_index(i,j))%lat
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
            zv2dx(i,j) = gg%cells%center(j)%lon
          ENDIF
        ENDDO
      ENDDO
      CALL nf(nf_put_var_double(ncid, varid_clatv, zv2dy))
      CALL nf(nf_put_var_double(ncid, varid_clonv, zv2dx))
      DEALLOCATE(zv2dx, zv2dy)

      CALL nf(nf_put_var_double(ncid, varid_clon, gg%cells%center(:)%lon))
      CALL nf(nf_put_var_double(ncid, varid_clat, gg%cells%center(:)%lat))

      !
      ! vertex part:
      !
      CALL nf(nf_put_var_double(ncid, varid_vlon,  gg%verts%vertex(:)%lon))
      CALL nf(nf_put_var_double(ncid, varid_vlat,  gg%verts%vertex(:)%lat))
      !

      ! Transpose of index array necessary for CF-1.1 Convention
      !
      ALLOCATE(zv2d(6, i_nv))
      DO j = 1, 6
        DO i = 1, i_nv
          IF (gg%verts%cell_index(i,j) == 0) THEN
            zv2d(7-j,i) = gg%cells%center(gg%verts%cell_index(i,5))%lon
          ELSE
            zv2d(7-j,i) = gg%cells%center(gg%verts%cell_index(i,j))%lon
          ENDIF
        ENDDO
      ENDDO
      CALL nf(nf_put_var_double(ncid, varid_vlonv, zv2d))
      DO j = 1, 6
        DO i = 1, i_nv
          IF (gg%verts%cell_index(i,j) == 0) THEN
            zv2d(7-j,i) = gg%cells%center(gg%verts%cell_index(i,5))%lat
          ELSE
            zv2d(7-j,i) = gg%cells%center(gg%verts%cell_index(i,j))%lat
          ENDIF
        ENDDO
      ENDDO
      CALL nf(nf_put_var_double(ncid, varid_vlatv, zv2d))
      DEALLOCATE (zv2d)
      !
      ! edge part:
      !
      CALL nf(nf_put_var_double(ncid, varid_elon,  gg%edges%center(:)%lon))
      CALL nf(nf_put_var_double(ncid, varid_elat,  gg%edges%center(:)%lat))
      !
      ! Transpose of index array necessary for CF-1.1 Convention
      !
      ALLOCATE(zv2dx(4, i_ne), zv2dy(4, i_ne))
      DO i = 1, i_ne
        zv2dx(1,i) = gg%verts%vertex(gg%edges%vertex_index(i,1))%lon
      ENDDO
      DO i = 1, i_ne
        zv2dx(3,i) = gg%verts%vertex(gg%edges%vertex_index(i,2))%lon
      ENDDO
      DO i = 1, i_ne
        zv2dx(4,i) = gg%cells%center(gg%edges%cell_index(i,1))%lon
      ENDDO
      DO i = 1, i_ne
        zv2dx(2,i) = gg%cells%center(gg%edges%cell_index(i,2))%lon
      ENDDO
      DO i = 1, i_ne
        zv2dy(1,i) = gg%verts%vertex(gg%edges%vertex_index(i,1))%lat
      ENDDO
      DO i = 1, i_ne
        zv2dy(3,i) = gg%verts%vertex(gg%edges%vertex_index(i,2))%lat
      ENDDO
      DO i = 1, i_ne
        zv2dy(4,i) = gg%cells%center(gg%edges%cell_index(i,1))%lat
      ENDDO
      DO i = 1, i_ne
        zv2dy(2,i) = gg%cells%center(gg%edges%cell_index(i,2))%lat
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
            zv2dx(i,j) = gg%edges%center(j)%lon
          ENDIF
        ENDDO
      ENDDO
      DO i = 1, i_ne
        IF (check_orientation(gg%edges%center(i)%lon, &
          & zv2dx(:,i),zv2dy(:,i),4) < 0) THEN
          swap(1:4) = zv2dx(4:1:-1,i)
          zv2dx(:,i) = swap(:)
          swap(1:4) = zv2dy(4:1:-1,i)
          zv2dy(:,i) = swap(:)
        ENDIF
        IF (check_orientation(gg%edges%center(i)%lon, &
          & zv2dx(:,i),zv2dy(:,i),4) < 0) THEN
        ENDIF
      ENDDO
      !
      CALL nf(nf_put_var_double(ncid, varid_elonv, zv2dx))
      CALL nf(nf_put_var_double(ncid, varid_elatv, zv2dy))
      !
      DEALLOCATE (zv2dx, zv2dy)
      !
      !-----------------------------------------------------------------------
      !
      CALL nf(nf_put_var_double(ncid, varid_carea, gg%cells%area))
      CALL nf(nf_put_var_double(ncid, varid_varea, gg%verts%dual_area))
      !
      !-----------------------------------------------------------------------
      !
      CALL nf(nf_put_var_double(ncid, varid1,  gg%cells%center(:)%lon))
      CALL nf(nf_put_var_double(ncid, varid2,  gg%cells%center(:)%lat))
      CALL nf(nf_put_var_double(ncid, varid3,  gg%verts%vertex(:)%lon))
      CALL nf(nf_put_var_double(ncid, varid4,  gg%verts%vertex(:)%lat))
      CALL nf(nf_put_var_double(ncid, varid5,  gg%edges%center(:)%lon))

      CALL nf(nf_put_var_double(ncid, varid_vx,    sphere_level%vs(:)%vertex%x(1)))
      CALL nf(nf_put_var_double(ncid, varid_vy,    sphere_level%vs(:)%vertex%x(2)))
      CALL nf(nf_put_var_double(ncid, varid_vz,    sphere_level%vs(:)%vertex%x(3)))

      CALL nf(nf_put_var_double(ncid, varid6,  gg%edges%center(:)%lat))
      CALL nf(nf_put_var_int   (ncid, varid7,  gg%cells%edge_index))
      CALL nf(nf_put_var_int   (ncid, varid8,  gg%cells%vertex_index))
      CALL nf(nf_put_var_int   (ncid, varid9,  gg%edges%cell_index))
      CALL nf(nf_put_var_int   (ncid, varid10, gg%edges%vertex_index))
      CALL nf(nf_put_var_int   (ncid, varid11, gg%verts%cell_index))
      CALL nf(nf_put_var_int   (ncid, varid12, gg%verts%edge_index))
      CALL nf(nf_put_var_int   (ncid, varid13, gg%verts%neighbor_index))
      CALL nf(nf_put_var_double(ncid, varid14, gg%cells%area))
      CALL nf(nf_put_var_double(ncid, varid15, gg%verts%dual_area))
      CALL nf(nf_put_var_double(ncid, varid16, gg%edges%primal_edge_length))
      CALL nf(nf_put_var_double(ncid, varid17, gg%edges%dual_edge_length))
      CALL nf(nf_put_var_double(ncid, varid18, gg%edges%edge_vert_length))
      CALL nf(nf_put_var_double(ncid, varid40, gg%edges%edge_cell_length))
      CALL nf(nf_put_var_double(ncid, varid19, gg%edges%primal_normal(:)%v1))
      CALL nf(nf_put_var_double(ncid, varid20, gg%edges%primal_normal(:)%v2))
      CALL nf(nf_put_var_double(ncid, varid21, gg%edges%dual_normal(:)%v1))
      CALL nf(nf_put_var_double(ncid, varid22, gg%edges%dual_normal(:)%v2) )
      CALL nf(nf_put_var_int   (ncid, varid23, gg%cells%edge_orientation))
      CALL nf(nf_put_var_int   (ncid, varid24, gg%cells%idx))
      CALL nf(nf_put_var_int   (ncid, varid25, gg%cells%parent_index))
      CALL nf(nf_put_var_int   (ncid, varid251, gg%cells%parent_child_type))
      CALL nf(nf_put_var_int   (ncid, varid26, gg%cells%neighbor_index))
      CALL nf(nf_put_var_int   (ncid, varid27, gg%cells%child_index))
      CALL nf(nf_put_var_int   (ncid, varid41, gg%cells%child_id))
      CALL nf(nf_put_var_int   (ncid, varid28, gg%edges%idx))
      !      CALL nf(nf_put_var_int   (ncid, varid281, gg%edges%parent_index))
      CALL nf(nf_put_var_int   (ncid, varid282, gg%edges%parent_child_type))
      CALL nf(nf_put_var_int   (ncid, varid29, gg%verts%idx))
      CALL nf(nf_put_var_int   (ncid, varid30, gg%verts%edge_orientation))
      CALL nf(nf_put_var_int   (ncid, varid31, gg%edges%system_orientation))

      ! Variables added for mesh refinement
      CALL nf(nf_put_var_int   (ncid, varid38, gg%edges%parent_index))
      CALL nf(nf_put_var_int   (ncid, varid39, gg%edges%child_index))
      CALL nf(nf_put_var_int   (ncid, varid42, gg%edges%child_id))
      CALL nf(nf_put_var_int   (ncid, varid32, gg%cells%refin_ctrl))
      CALL nf(nf_put_var_int   (ncid, varid33, gg%cells%indlist))
      CALL nf(nf_put_var_int   (ncid, varid43, gg%cells%start_idx))
      CALL nf(nf_put_var_int   (ncid, varid44, gg%cells%end_idx))
      CALL nf(nf_put_var_int   (ncid, varid34, gg%edges%refin_ctrl))
      CALL nf(nf_put_var_int   (ncid, varid35, gg%edges%indlist))
      CALL nf(nf_put_var_int   (ncid, varid45, gg%edges%start_idx))
      CALL nf(nf_put_var_int   (ncid, varid46, gg%edges%end_idx))
      CALL nf(nf_put_var_int   (ncid, varid36, gg%verts%refin_ctrl))
      CALL nf(nf_put_var_int   (ncid, varid37, gg%verts%indlist))
      CALL nf(nf_put_var_int   (ncid, varid47, gg%verts%start_idx))
      CALL nf(nf_put_var_int   (ncid, varid48, gg%verts%end_idx))

      CALL nf(nf_put_var_int   (ncid, varid403,gg%verts%parent_index))
      CALL nf(nf_put_var_int   (ncid, varid49, gg%cells%phys_id))
      CALL nf(nf_put_var_int   (ncid, varid50, gg%edges%phys_id))
      !------------------------------------------------------------------------

      CALL nf(nf_close(ncid))

      IF (jgrid == annotate_level) THEN
        CALL dump_grid_table_xml_metadata(TRIM(output))
      ENDIF
    !------------------------------------------------------------------------
    !write(*,*) TRIM(output)
    !str_idx=LBOUND(gg%verts%start_idx, 1)
    !end_idx=str_idx+SIZE(gg%verts%start_idx, 1)-1
    !DO i=str_idx,end_idx
    !  write(*,*) 'verts%start_idx, end:', i, gg%verts%start_idx(i,1), gg%verts%end_idx(i,1)
    !ENDDO
    !
    !str_idx=LBOUND(gg%edges%start_idx, 1)
    !end_idx=str_idx+SIZE(gg%edges%start_idx, 1)-1
    !DO i=str_idx,end_idx
    !  write(*,*) 'edges%start_idx, end:', i, gg%edges%start_idx(i,1), gg%edges%end_idx(i,1)
    !ENDDO
    !
    !str_idx=LBOUND(gg%cells%start_idx, 1)
    !end_idx=str_idx+SIZE(gg%cells%start_idx, 1)-1
    !DO i=str_idx,end_idx
    !  write(*,*) 'cells%start_idx, end:', i, gg%cells%start_idx(i,1), gg%cells%end_idx(i,1)
    !ENDDO
    !write(*,*) '-------------------'

    ENDDO

  END SUBROUTINE output_grid

  !-------------------------------------------------------------------------
  !>
  !!  Computes statistics on the whole hierarchy of icosahedral grids.
  !!
  !!
  !! @par Revision History
  !!  Developed and tested  by L.Bonaventura  (2004)
  !!  Modified for new data structure by Th.Heinze (2005-10)
  !!
  SUBROUTINE statistics(lplane)
    !
    LOGICAL, INTENT(in) :: lplane
    TYPE(t_grid), POINTER :: ptr_gg

    REAL(wp) :: z_emx, z_emn, z_dmn, z_dmx, z_demxloc, z_demnloc, z_dmnloc, &
      & z_pemxloc, z_pemnloc, z_peratio, z_pfratio,                        &
      & z_dmxloc, z_dav, z_amn, z_amx, z_atot, z_dfratio, z_deratio,       &
      & z_coemax, z_coedge, z_colen, z_edist, z_fdist, z_farea, z_tarea,    &
      & z_pav

    INTEGER :: i_ngrids, igrid, istat, i_nc, i_ne, i_nv, ie0, iv0, jj, if1

    CHARACTER(LEN=filename_max) :: output
    CHARACTER(LEN=7) :: optfix

    !-------------------------------------------------------------------------
    !BOC

    ! mark Heikes-Randall + additional C-grid optimization in last step

    i_ngrids=number_of_grid_levels

    SELECT CASE (itype_optimize)
    CASE (0)
      ! natural grid
      optfix = 'noo'
    CASE (1)
      IF (l_c_grid) THEN
        ! Heikes-Randall with additional c-grid optimization
        optfix = 'hrc'
      ELSE
        ! Heikes-Randall optimization
        optfix = 'hro'
      ENDIF
    CASE (2)
      ! equal area subdivision
      IF (l_c_grid) THEN
        optfix = 'eac'
      ELSE
        optfix = 'eas'
      ENDIF
    CASE (3)
      ! c-grid small circle constraint
      optfix = 'scc'
    CASE (4)
      ! spring dynamics
      IF (l_c_grid) THEN
        optfix = 'spc'
      ELSE
        WRITE(optfix,'(a3,f4.2)') 'spr',beta_spring
      ENDIF
    END SELECT

    !
    ! loop over all grids
    !
    DO igrid = 0, i_ngrids

      IF(lplane) THEN
        WRITE (output,'(a,i0,a,i2.2,a,a)') &
          & 'planR', grid_root, 'B', igrid , '-statistics'
      ELSE
        IF (igrid >= maxlev_optim) THEN
          WRITE (output,'(a,i0,a,i2.2,3a,i1)') &
            & 'iconR', grid_root, 'B', igrid , '-statistics_', TRIM(optfix), '_M', maxlev_optim
        ELSE
          WRITE (output,'(a,i0,a,i2.2,a,a)') &
            & 'iconR', grid_root, 'B', igrid , '-statistics_', TRIM(optfix)
        ENDIF
      ENDIF

      OPEN(nstat,FILE=TRIM(output), STATUS='unknown',IOSTAT=istat)
      IF (istat>0) THEN
        CALL finish  ('statistics', &
          & 'error in opening grid statistics file '//TRIM(output))
      ENDIF

      ptr_gg => grid_on_level(igrid)            ! point on current grid
      i_nc= ptr_gg%ncells                       ! get current number of cells
      i_ne= ptr_gg%nedges                       ! get current number of edges
      i_nv= ptr_gg%nverts                       ! get current number of vertices
      !
      ! set defaults
      !
      z_coemax = 0._wp              ! maximum offcentering on tri C grid
      z_demnloc= 10000000.0_wp      ! minimum local distance tri cell centers
      z_demxloc= 0.0_wp             ! maximum local distance tri cell centers
      z_pemnloc= 10000000.0_wp      ! minimum edge length of a primal cell
      z_pemxloc= 0.0_wp             ! maximum edge length of a primal cell
      z_dmnloc = 10000000.0_wp      ! minimum local distance hex cell centers
      z_dmxloc = 0.0_wp             ! maximum local distance hex cell centers
      z_emn    = 10000000.0_wp      ! minimum distance tri cell centers
      z_emx    = 0.0_wp             ! maximum distance tri cell centers
      z_dmn    = 10000000.0_wp      ! minimum distance hex cell centers
      z_dmx    = 0.0_wp             ! maximum distance hex cell centers
      z_dav    = 0.0_wp             ! average resolution [rad] (triangles)
      z_pav    = 0.0_wp             ! average resolution [rad] (hexagons)
      z_amn    = 1.0e15_wp          ! minimum triangular/hexagon area
      z_amx    = 0.0_wp             ! maximum triangular/hexagon area
      z_atot   = 0._wp              ! total area of triangles/hexagons
      z_deratio= 1.0_wp             ! minimum of local ratios between
      ! min/max distances of tri cell centers
      z_peratio= 1.0_wp             ! minimum of local ratios between
      ! min/max edge lengths of tri cells
      z_dfratio= 1.0_wp             ! minimum of local ratios between
      ! min/max distances of hex cell centers
      z_pfratio= 1.0_wp             ! minimum of local ratios between
      ! min/max edge lenghts of hex cells
      !
      ! loop over all edges
      !
      DO ie0 = 1, i_ne
        z_edist=ptr_gg%edges%dual_edge_length(ie0)   ! distance tri cell centers
        z_fdist=ptr_gg%edges%primal_edge_length(ie0) ! distance hex cell centers
        z_emn = MIN(z_emn,z_edist)
        z_emx = MAX(z_emx,z_edist)
        z_dmn = MIN(z_dmn,z_fdist)
        z_dmx = MAX(z_dmx,z_fdist)
        z_dav = z_dav + z_edist                      ! sum up the dual edge lengths
        z_pav = z_pav + z_fdist                      ! sum up the primal edge lengths
        !
        ! calculate coemax
        !
        z_colen=ptr_gg%edges%edge_cell_length(ie0,1)
        z_coedge    = z_colen/z_edist
        z_coedge    = ABS(z_coedge - 0.5_wp)  ! 0.5 relates to no offcentering
        z_coemax    = MAX(z_coemax,z_coedge)
        z_colen=ptr_gg%edges%edge_cell_length(ie0,2)
        z_coedge    = z_colen/z_edist
        z_coedge    = ABS(z_coedge - 0.5_wp)  ! 0.5 relates to no offcentering
        z_coemax    = MAX(z_coemax,z_coedge)
      ENDDO
      !
      z_dav = z_dav/REAL(i_ne,wp)            ! calculate the average resolution (tri)
      z_pav = z_pav/REAL(i_ne,wp)            ! calculate the average resolution (hex)
      !
      ! loop over all cells
      !
      DO iv0 =1,i_nc
        z_demnloc=10000000.0_wp                       ! set default
        z_demxloc=0.0_wp                              ! set default
        z_pemnloc=10000000.0_wp                       ! set default
        z_pemxloc=0.0_wp                              ! set default
        DO jj=1,3                                     ! loop over all 3 edges
          ie0=ptr_gg%cells%edge_index(iv0,jj)         ! edge index

          z_edist=ptr_gg%edges%dual_edge_length(ie0)  ! dist. tri centers
          z_demnloc=MIN(z_edist,z_demnloc)
          z_demxloc=MAX(z_edist,z_demxloc)

          z_edist=ptr_gg%edges%primal_edge_length(ie0) ! triangle edge length
          z_pemnloc=MIN(z_edist,z_pemnloc)
          z_pemxloc=MAX(z_edist,z_pemxloc)
        ENDDO
        z_deratio=MIN(z_demnloc/z_demxloc,z_deratio)
        z_peratio=MIN(z_pemnloc/z_pemxloc,z_peratio)
      ENDDO
      !
      ! loop over all vertices
      !
      DO if1 =1,i_nv
        z_dmnloc=10000000.0_wp                         ! set default
        z_dmxloc=0.0_wp                                ! set default
        z_pemnloc=10000000.0_wp                        ! set default
        z_pemxloc=0.0_wp                               ! set default
        DO jj=1,5                                      ! loop over first 5 edges
          ie0=ptr_gg%verts%edge_index(if1,jj)          ! edge index

          z_fdist=ptr_gg%edges%primal_edge_length(ie0) ! dist. pent/hex centers
          z_dmnloc=MIN(z_fdist,z_dmnloc)
          z_dmxloc=MAX(z_fdist,z_dmxloc)

          z_fdist=ptr_gg%edges%dual_edge_length(ie0)   ! pent/hex edge length
          z_pemnloc=MIN(z_fdist,z_pemnloc)
          z_pemxloc=MAX(z_fdist,z_pemxloc)
        ENDDO
        !
        ! get 6th edge, if if1 relates to pentagon, the index is 0
        !
        ie0=ptr_gg%verts%edge_index(if1,6)
        !
        ! check if pentagon or hexagon
        !
        IF (ie0 > 0 ) THEN                              ! hexagon
          z_fdist=ptr_gg%edges%primal_edge_length(ie0)
          z_dmnloc=MIN(z_fdist,z_dmnloc)
          z_dmxloc=MAX(z_fdist,z_dmxloc)

          z_fdist=ptr_gg%edges%dual_edge_length(ie0)
          z_pemnloc=MIN(z_fdist,z_pemnloc)
          z_pemxloc=MAX(z_fdist,z_pemxloc)
        ENDIF
        z_dfratio=MIN(z_dmnloc/z_dmxloc,z_dfratio)
        z_pfratio=MIN(z_pemnloc/z_pemxloc,z_pfratio)
      ENDDO
      !
      ! statistics output part 1
      !
      WRITE (nstat,'(a,f12.3)') 'average distance tri cell centers [km]: ', &
        & 0.001_wp*z_dav
      WRITE (nstat,'(a,f12.3)') 'minimum distance tri cell centers [km]: ', &
        & 0.001_wp*z_emn
      WRITE (nstat,'(a,f12.3)') 'maximum distance tri cell centers [km]: ', &
        & 0.001_wp*z_emx
      WRITE (nstat,'(a,f12.3)') 'ratio of distance tri cell centers [%]: ', &
        & 100._wp*z_emn/z_emx
      WRITE (nstat,'(2a,f7.3)') 'minimum of local ratios between min/max ', &
        & 'distances of tri cell centers [%]: ',      &
        & z_deratio*100._wp
      WRITE (nstat,'(2a,f7.3)') 'minimum of local ratios between min/max ', &
        & 'triangle edge lengths [%]: ',              &
        & z_peratio*100._wp
      WRITE (nstat,*)
      WRITE (nstat,'(a,f12.3)') 'average distance hex cell centers [km]: ', &
        & 0.001_wp*z_pav
      WRITE (nstat,'(a,f12.3)') 'minimum distance hex cell centers [km]: ', &
        & 0.001_wp*z_dmn
      WRITE (nstat,'(a,f12.3)') 'maximum distance hex cell centers [km]: ', &
        & 0.001_wp*z_dmx
      WRITE (nstat,'(a,f12.3)') 'ratio of distance hex cell centers [%]: ', &
        & 100._wp*z_dmn/z_dmx
      WRITE (nstat,'(2a,f7.3)') 'minimum of local ratios between min/max ', &
        & 'distances of hex cell centers [%]: ',      &
        & z_dfratio*100._wp
      WRITE (nstat,'(2a,f7.3)') 'minimum of local ratios between min/max ', &
        & 'hexagon edge lengths [%]: ',               &
        & z_pfratio*100._wp
      WRITE (nstat,*)
      WRITE (nstat,'(a,f12.3)') 'maximum offcentering on tri C grid [%]: ', &
        & 100._wp*z_coemax
      !
      ! loop over all vertices
      !
      DO if1 = 1, i_nv
        z_farea = ptr_gg%verts%dual_area(if1)         ! area of pent/hex
        z_atot  = z_atot + z_farea
        z_amn   = MIN(z_amn,z_farea)
        z_amx   = MAX(z_amx,z_farea)
      ENDDO
      !
      ! statistics output part 2
      !
      WRITE (nstat,*)
      WRITE (nstat,'(a,f19.16,a,f19.16,a)') 'total area of hexagons = ', &
        & z_atot/earth_radius/earth_radius,' (', 4.0_wp*pi, ')'
      WRITE (nstat,'(a,f11.1)') 'minimum hexagon area         [km**2]: ', &
        & z_amn*1e-6_wp
      WRITE (nstat,'(a,f11.1)') 'maximum hexagon area         [km**2]: ', &
        & z_amx*1e-6_wp
      WRITE (nstat,'(a,f11.1)') 'average area                 [km**2]: ', &
        & z_atot*1e-6_wp/REAL(i_nv,wp)
      WRITE (nstat,'(a,f11.1)') 'ratio of minimum to maximum area [%]: ', &
        & 100._wp*z_amn/z_amx
      !
      ! set defaults
      !
      z_amn  = 1.0e15_wp
      z_amx  = 0.0_wp
      z_atot = 0.0_wp
      !
      ! loop over all cells
      !
      DO iv0 = 1 , i_nc
        z_tarea = ptr_gg%cells%area(iv0)                    ! area of triangle
        z_atot  = z_atot + z_tarea
        z_amn   = MIN(z_amn,z_tarea)
        z_amx   = MAX(z_amx,z_tarea)
      ENDDO
      !
      ! statistics output part 3
      !
      WRITE (nstat,*)
      WRITE (nstat,'(a,f19.16,a,f19.16,a)') 'total area of triangles = ', &
        & z_atot/earth_radius/earth_radius, ' (', 4._wp*pi, ')'
      WRITE (nstat,'(a,f11.1)') 'minimum triangle area        [km**2]: ', &
        & z_amn*1.e-6_wp
      WRITE (nstat,'(a,f11.1)') 'maximum triangle  area       [km**2]: ', &
        & z_amx*1.e-6_wp
      WRITE (nstat,'(a,f11.1)') 'average area                 [km**2]: ', &
        & z_atot*1e-6_wp/REAL(i_nc,wp)
      WRITE (nstat,'(a,f11.1)') 'ratio of minimum to maximum area [%]: ', &
        & 100._wp*z_amn/z_amx
      WRITE (nstat,*)
      !
    ENDDO
    !
    CLOSE(nstat)

  END SUBROUTINE statistics
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  SUBROUTINE nf(STATUS)
    INTEGER, INTENT(in) :: STATUS

    IF (STATUS /= nf_noerr) THEN
      CALL finish('mo_io_grid netCDF error', nf_strerror(STATUS))
    ENDIF

  END SUBROUTINE nf


END MODULE mo_grid_levels

