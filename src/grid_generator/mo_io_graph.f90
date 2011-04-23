!>
!!  Contains the IO subroutines for the graph.
!!
!! The subroutines
!!  <i>output_graph</i> and <i>input_graph</i> write and read
!!  graph mapping netCDF files.
!!
!! @par Revision History
!!  Developed  by Thomas Heinze, DWD (2007-05-10) to avoid cyclic dependencies
!!  in Makefile. This would have be the case if these files were contained in
!!  mo_io_grid
!! Modification by Luis Kornblueh (2008-04-21):
!! - rearrange for use of netCDF and naming scheme for map files to explain
!!   content as are root division and level. The resulting file is
!!   architectur independent and needs to be generated once only.
!! Modification by Almut Gassmann (2009-01-08)
!! - adapt for generation of planar region
!! Modification by Marco Giorgetta (2009-05-08)
!! - set default resolution to R2B04
!!
!! @par Copyright
!! 2002-2008 by DWD and MPI-M
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
MODULE mo_io_graph

  !-------------------------------------------------------------------------
  !
  !    ProTeX FORTRAN source: Style 2
  !    modified for ICON project, DWD/MPI-M 2004-5
  !
  !-------------------------------------------------------------------------
  !
  !
  !

  USE mo_io_units,           ONLY: filename_max, nnml
  USE mo_exception,          ONLY: message_text, message, finish
  USE mo_namelist,           ONLY: position_nml, open_nml, close_nml, positioned
  USE mo_base_datatypes,     ONLY: t_triangle, edge, vertex, t_spheres,   &
                                   dummy_e, dummy_c, dummy_v

  IMPLICIT NONE

  PRIVATE

  INCLUDE 'netcdf.inc'

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

  PUBLIC output_graph, input_graph, init_graphgen

  !--------------------------------------------------------------------
  !BOC

CONTAINS

  !EOC
  !-------------------------------------------------------------------------
  !
  !

  !>
  !!  Initializes the graph generator run,.
  !!
  !!  Initializes the graph generator run,
  !!  reading namelists from the <b>GRID_INI</b> file
  !!  in the directory <b>run</b>.
  !!
  SUBROUTINE init_graphgen( grid_levels, nroot, lplane)

    INTEGER, INTENT(OUT)::grid_levels, nroot
    LOGICAL, INTENT(OUT)::lplane

    INTEGER :: i_status, nroot_d, grid_levels_d
    LOGICAL :: lplane_d


    !--------------------------------------------------------------------
    !BOC
    !
    ! ! initialize with grid parameters from files

    NAMELIST /graph_ini/ grid_levels, nroot, lplane

    ! set default values for graph_ini

    nroot       = 2
    grid_levels = 4
    lplane      = .FALSE. ! graph on the sphere is the default

    ! copy default values to "default" variables

    nroot_d          = nroot
    grid_levels_d    = grid_levels
    lplane_d         = lplane

    ! read namelist

    CALL open_nml('NAMELIST_GRAPH')

    CALL position_nml ('graph_ini', lrewind=.FALSE., status=i_status)
    IF (i_status == POSITIONED) THEN
      READ (nnml,graph_ini)
    ENDIF

    CALL close_nml

  END SUBROUTINE init_graphgen

  !EOC
  !-------------------------------------------------------------------------
  !
  !
  !>
  !! Read <b>GRAPHMAP.
  !!
  !! *</b> file with cross reference tables of graph items,
  !! where * denotes the grid level.
  !! The file has sequential access, unformatted input.
  !!
  SUBROUTINE input_graph(k_root, sol, lplane)
    !

    INTEGER,       INTENT(in) :: k_root
    TYPE(t_spheres), INTENT(inout) :: sol
    LOGICAL,       INTENT(in) :: lplane


    TYPE(t_triangle), POINTER :: ct(:) ! local array of triangles
    TYPE(edge)    , POINTER :: ce(:) ! local array of edges
    TYPE(vertex)  , POINTER :: cv(:) ! local array of vertices

    ! local arrays associated with triangles
    INTEGER, ALLOCATABLE    :: tvert0(:), tvert1(:), tvert2(:)
    INTEGER, ALLOCATABLE    :: tedge0(:), tedge1(:), tedge2(:)

    ! local arrays associated with edges
    INTEGER, ALLOCATABLE    :: evert0(:), evert1(:)
    INTEGER, ALLOCATABLE    :: etria0(:), etria1(:)

    ! local arrays associated with vertices
    INTEGER, ALLOCATABLE    :: vedge0(:), vedge1(:), vedge2(:)
    INTEGER, ALLOCATABLE    :: vedge3(:), vedge4(:), vedge5(:)
    INTEGER, ALLOCATABLE    :: vtria0(:), vtria1(:), vtria2(:)
    INTEGER, ALLOCATABLE    :: vtria3(:), vtria4(:), vtria5(:)
    INTEGER, ALLOCATABLE    :: vneig0(:), vneig1(:), vneig2(:)
    INTEGER, ALLOCATABLE    :: vneig3(:), vneig4(:), vneig5(:)

    CHARACTER(len=filename_max) :: input

    INTEGER :: ist, istat
    INTEGER :: j, jgrid, i_check, i_spec_pt
    INTEGER :: i_not, i_nov, i_noe

    INTEGER :: ncid, dimid, varid
    INTEGER :: start(2), count(2)

    !-------------------------------------------------------------------------
    !BOC

    ! set dummy edge, vertex and cells(triangles)

    dummy_e%info%idx=-1
    dummy_c%info%idx=-1
    dummy_v%info%idx=-1

    dummy_e%vertex0=>dummy_v
    dummy_e%vertex1=>dummy_v
    dummy_v%triangle5=>dummy_c
    dummy_v%neighbor5=>dummy_v
    dummy_v%edge5=>dummy_e

    ! point local arrays to global one

    ct => sol%ts
    ce => sol%es
    cv => sol%vs

    ! number of triangles, vetices and edges

    i_not = sol%no_triangles
    i_nov = sol%no_vertices
    i_noe = sol%no_edges

    ! allocate more local arrays

    ALLOCATE(tvert0(0:i_not-1), STAT=istat)
    ist = istat
    ALLOCATE(tvert1(0:i_not-1), STAT=istat)
    ist = ist + istat
    ALLOCATE(tvert2(0:i_not-1), STAT=istat)
    ist = ist + istat
    ALLOCATE(tedge0(0:i_not-1), STAT=istat)
    ist = ist + istat
    ALLOCATE(tedge1(0:i_not-1), STAT=istat)
    ist = ist + istat
    ALLOCATE(tedge2(0:i_not-1), STAT=istat)
    ist = ist + istat
    ALLOCATE(evert0(0:i_noe-1), STAT=istat)
    ist = ist + istat
    ALLOCATE(evert1(0:i_noe-1), STAT=istat)
    ist = ist + istat
    ALLOCATE(etria0(0:i_noe-1), STAT=istat)
    ist = ist + istat
    ALLOCATE(etria1(0:i_noe-1), STAT=istat)
    ist = ist + istat
    ALLOCATE(vedge0(0:i_nov-1), STAT=istat)
    ist = ist + istat
    ALLOCATE(vedge1(0:i_nov-1), STAT=istat)
    ist = ist + istat
    ALLOCATE(vedge2(0:i_nov-1), STAT=istat)
    ist = ist + istat
    ALLOCATE(vedge3(0:i_nov-1), STAT=istat)
    ist = ist + istat
    ALLOCATE(vedge4(0:i_nov-1), STAT=istat)
    ist = ist + istat
    ALLOCATE(vedge5(0:i_nov-1), STAT=istat)
    ist = ist + istat
    ALLOCATE(vtria0(0:i_nov-1), STAT=istat)
    ist = ist + istat
    ALLOCATE(vtria1(0:i_nov-1), STAT=istat)
    ist = ist + istat
    ALLOCATE(vtria2(0:i_nov-1), STAT=istat)
    ist = ist + istat
    ALLOCATE(vtria3(0:i_nov-1), STAT=istat)
    ist = ist + istat
    ALLOCATE(vtria4(0:i_nov-1), STAT=istat)
    ist = ist + istat
    ALLOCATE(vtria5(0:i_nov-1), STAT=istat)
    ist = ist + istat
    ALLOCATE(vneig0(0:i_nov-1), STAT=istat)
    ist = ist + istat
    ALLOCATE(vneig1(0:i_nov-1), STAT=istat)
    ist = ist + istat
    ALLOCATE(vneig2(0:i_nov-1), STAT=istat)
    ist = ist + istat
    ALLOCATE(vneig3(0:i_nov-1), STAT=istat)
    ist = ist + istat
    ALLOCATE(vneig4(0:i_nov-1), STAT=istat)
    ist = ist + istat
    ALLOCATE(vneig5(0:i_nov-1), STAT=istat)
    ist = ist + istat

    IF (ist>0 ) THEN
      WRITE (message_text,'(a)')    'Problem in allocating local arrays'
      CALL message ('',TRIM(message_text))
      CALL finish ('mo_io_graph/input_graph', message_text)
    ENDIF

    ! open GRAPHMAP file

    jgrid = sol%level

    IF (.NOT.lplane) THEN
      WRITE (input,'(a,i0,a,i2.2,a)') 'iconR', k_root, 'B', jgrid, '-graph.nc'
    ELSE
      WRITE (input,'(a,i0,a,i2.2,a)') 'planR', k_root, 'B', jgrid, '-graph.nc'
    ENDIF

    WRITE(message_text,'(a,a)') 'Read graphmap file ', TRIM(input)
    CALL message ('', TRIM(message_text))

    !---------------------------------------------------------------------
    !
    CALL nf(nf_open(TRIM(input), NF_NOWRITE, ncid))
    !
    CALL nf(nf_get_att_int(ncid, NF_GLOBAL, 'grid_root', i_check))
    IF (i_check /= k_root) THEN
      WRITE(message_text,'(a,a)')  'The "grid_root" global attribute in the ',&
           &  'graph file differs from the "R" parameter in the filename'
      CALL finish  ('mo_io_graph/input_graph', TRIM(message_text))
    END IF
    !
    CALL nf(nf_get_att_int(ncid, NF_GLOBAL, 'grid_level', i_check))
    IF (i_check /= jgrid) THEN
      WRITE(message_text,'(a,a)')  'The "grid_level" global attribute in the ',&
           & 'graph file differs from the "B" parameter in the filename'
      CALL finish  ('mo_io_graph/input_graph', TRIM(message_text))
    ENDIF
    !
    CALL nf(nf_inq_dimid(ncid, 'cell', dimid))
    CALL nf(nf_inq_dimlen(ncid, dimid, i_check))
    IF (i_check /= i_not) THEN
      WRITE(message_text,'(a)')  'inconsistent number of triangles in GRAPHMAP file'
      CALL finish  ('mo_io_graph/input_graph', TRIM(message_text))
    ENDIF
    !
    CALL nf(nf_inq_dimid(ncid, 'edge', dimid))
    CALL nf(nf_inq_dimlen(ncid, dimid, i_check))
    IF (i_check /= i_noe) THEN
      WRITE(message_text,'(a)')  'inconsistent number of edges in GRAPHMAP file'
      CALL finish  ('mo_io_graph/input_graph', TRIM(message_text))
    ENDIF
    !
    CALL nf(nf_inq_dimid(ncid, 'vertex', dimid))
    CALL nf(nf_inq_dimlen(ncid, dimid, i_check))
    IF (i_check /= i_nov) THEN
      WRITE(message_text,'(a)')  'inconsistent number of vertices in GRAPHMAP file'
      CALL finish  ('mo_io_graph/input_graph', TRIM(message_text))
    ENDIF
    !
    !---------------------------------------------------------------------
    !
    ! Triangle vertex indices
    count = (/ i_not, 1 /)
    CALL nf(nf_inq_varid(ncid, 'triangle_vertex_index', varid))
    start = (/ 1, 1 /)
    CALL nf(nf_get_vara_int(ncid, varid, start, count, tvert0))
    start = (/ 1, 2 /)
    CALL nf(nf_get_vara_int(ncid, varid, start, count, tvert1))
    start = (/ 1, 3 /)
    CALL nf(nf_get_vara_int(ncid, varid, start, count, tvert2))
    !
    ! Triangle edge indices
    count = (/ i_not, 1 /)
    CALL nf(nf_inq_varid(ncid, 'triangle_edge_index', varid))
    start = (/ 1, 1 /)
    CALL nf(nf_get_vara_int(ncid, varid, start, count, tedge0))
    start = (/ 1, 2 /)
    CALL nf(nf_get_vara_int(ncid, varid, start, count, tedge1))
    start = (/ 1, 3 /)
    CALL nf(nf_get_vara_int(ncid, varid, start, count, tedge2))
    !
    ! Edge vertex indices
    count = (/ i_noe, 1 /)
    CALL nf(nf_inq_varid(ncid, 'edge_vertex_index', varid))
    start = (/ 1, 1 /)
    CALL nf(nf_get_vara_int(ncid, varid, start, count, evert0))
    start = (/ 1, 2 /)
    CALL nf(nf_get_vara_int(ncid, varid, start, count, evert1))
    !
    ! Edge triangle indices
    count = (/ i_noe, 1 /)
    CALL nf(nf_inq_varid(ncid, 'edge_triangle_index', varid))
    start = (/ 1, 1 /)
    CALL nf(nf_get_vara_int(ncid, varid, start, count, etria0))
    start = (/ 1, 2 /)
    CALL nf(nf_get_vara_int(ncid, varid, start, count, etria1))
    !
    ! Vertex edge indices
    count = (/ i_nov, 1 /)
    CALL nf(nf_inq_varid(ncid, 'vertex_edge_index', varid))
    start = (/ 1, 1 /)
    CALL nf(nf_get_vara_int(ncid, varid, start, count, vedge0))
    start = (/ 1, 2 /)
    CALL nf(nf_get_vara_int(ncid, varid, start, count, vedge1))
    start = (/ 1, 3 /)
    CALL nf(nf_get_vara_int(ncid, varid, start, count, vedge2))
    start = (/ 1, 4 /)
    CALL nf(nf_get_vara_int(ncid, varid, start, count, vedge3))
    start = (/ 1, 5 /)
    CALL nf(nf_get_vara_int(ncid, varid, start, count, vedge4))
    start = (/ 1, 6 /)
    CALL nf(nf_get_vara_int(ncid, varid, start, count, vedge5))
    !
    ! Vertex triangle indices
    count = (/ i_nov, 1 /)
    CALL nf(nf_inq_varid(ncid, 'vertex_triangle_index', varid))
    start = (/ 1, 1 /)
    CALL nf(nf_get_vara_int(ncid, varid, start, count, vtria0))
    start = (/ 1, 2 /)
    CALL nf(nf_get_vara_int(ncid, varid, start, count, vtria1))
    start = (/ 1, 3 /)
    CALL nf(nf_get_vara_int(ncid, varid, start, count, vtria2))
    start = (/ 1, 4 /)
    CALL nf(nf_get_vara_int(ncid, varid, start, count, vtria3))
    start = (/ 1, 5 /)
    CALL nf(nf_get_vara_int(ncid, varid, start, count, vtria4))
    start = (/ 1, 6 /)
    CALL nf(nf_get_vara_int(ncid, varid, start, count, vtria5))
    !
    ! Vertex neighbor indices
    count = (/ i_nov, 1 /)
    CALL nf(nf_inq_varid(ncid, 'vertex_neighbor_index', varid))
    start = (/ 1, 1 /)
    CALL nf(nf_get_vara_int(ncid, varid, start, count, vneig0))
    start = (/ 1, 2 /)
    CALL nf(nf_get_vara_int(ncid, varid, start, count, vneig1))
    start = (/ 1, 3 /)
    CALL nf(nf_get_vara_int(ncid, varid, start, count, vneig2))
    start = (/ 1, 4 /)
    CALL nf(nf_get_vara_int(ncid, varid, start, count, vneig3))
    start = (/ 1, 5 /)
    CALL nf(nf_get_vara_int(ncid, varid, start, count, vneig4))
    start = (/ 1, 6 /)
    CALL nf(nf_get_vara_int(ncid, varid, start, count, vneig5))
    !
    !---------------------------------------------------------------------
    !
    CALL nf(nf_close(ncid))
    !
    !---------------------------------------------------------------------

    ! set triangle index

    DO j = 0, i_not-1
      ct(j)%info%idx = j
    ENDDO

    ! point triangle vertices and edges to correct counterparts

    DO j = 0, i_not-1
      ct(j)%vertex0 => cv(tvert0(j))
      ct(j)%vertex1 => cv(tvert1(j))
      ct(j)%vertex2 => cv(tvert2(j))
      ct(j)%edge0   => ce(tedge0(j))
      ct(j)%edge1   => ce(tedge1(j))
      ct(j)%edge2   => ce(tedge2(j))
    ENDDO

    ! set edge index

    DO j = 0, i_noe-1
      ce(j)%info%idx = j
    ENDDO

    ! point edge vertices and triangles to correct counterparts

    DO j = 0, i_noe-1
      ce(j)%vertex0   => cv(evert0(j))
      ce(j)%vertex1   => cv(evert1(j))
      ce(j)%triangle0 => ct(etria0(j))
      ce(j)%triangle1 => ct(etria1(j))
    ENDDO

    ! set vertex index

    DO j = 0, i_nov-1
      cv(j)%info%idx = j
    ENDDO

    ! point vertex edges, triangles and neighbors to correct counterparts

    i_spec_pt = 0

    DO j = 0, i_nov-1
      cv(j)%edge0     => ce(vedge0(j))
      cv(j)%edge1     => ce(vedge1(j))
      cv(j)%edge2     => ce(vedge2(j))
      cv(j)%edge3     => ce(vedge3(j))
      cv(j)%edge4     => ce(vedge4(j))
      cv(j)%triangle0 => ct(vtria0(j))
      cv(j)%triangle1 => ct(vtria1(j))
      cv(j)%triangle2 => ct(vtria2(j))
      cv(j)%triangle3 => ct(vtria3(j))
      cv(j)%triangle4 => ct(vtria4(j))
      cv(j)%neighbor0 => cv(vneig0(j))
      cv(j)%neighbor1 => cv(vneig1(j))
      cv(j)%neighbor2 => cv(vneig2(j))
      cv(j)%neighbor3 => cv(vneig3(j))
      cv(j)%neighbor4 => cv(vneig4(j))

      ! special point treatment

      IF (vedge5(j) <0) THEN              ! special point found
        i_spec_pt = i_spec_pt + 1
        cv(j)%edge5     => dummy_e
        cv(j)%triangle5 => dummy_c
        cv(j)%neighbor5 => dummy_v
      ELSE                                ! no special point
        cv(j)%edge5     => ce(vedge5(j))
        cv(j)%triangle5 => ct(vtria5(j))
        cv(j)%neighbor5 => cv(vneig5(j))
      ENDIF
    ENDDO

    IF (.NOT.lplane) THEN
      IF (i_spec_pt /= 12) THEN
        WRITE (message_text,'(i0,a)') i_spec_pt, &
             ' special points found. Should be 12.'
        CALL finish ('mo_io_graph/input_graph', message_text)
      ENDIF
    ENDIF

    ! deallocate local arrays

    DEALLOCATE(tvert0, STAT=istat)
    ist = istat
    DEALLOCATE(tvert1, STAT=istat)
    ist = ist + istat
    DEALLOCATE(tvert2, STAT=istat)
    ist = ist + istat
    DEALLOCATE(tedge0, STAT=istat)
    ist = ist + istat
    DEALLOCATE(tedge1, STAT=istat)
    ist = ist + istat
    DEALLOCATE(tedge2, STAT=istat)
    ist = ist + istat
    DEALLOCATE(evert0, STAT=istat)
    ist = ist + istat
    DEALLOCATE(evert1, STAT=istat)
    ist = ist + istat
    DEALLOCATE(etria0, STAT=istat)
    ist = ist + istat
    DEALLOCATE(etria1, STAT=istat)
    ist = ist + istat
    DEALLOCATE(vedge0, STAT=istat)
    ist = ist + istat
    DEALLOCATE(vedge1, STAT=istat)
    ist = ist + istat
    DEALLOCATE(vedge2, STAT=istat)
    ist = ist + istat
    DEALLOCATE(vedge3, STAT=istat)
    ist = ist + istat
    DEALLOCATE(vedge4, STAT=istat)
    ist = ist + istat
    DEALLOCATE(vedge5, STAT=istat)
    ist = ist + istat
    DEALLOCATE(vtria0, STAT=istat)
    ist = ist + istat
    DEALLOCATE(vtria1, STAT=istat)
    ist = ist + istat
    DEALLOCATE(vtria2, STAT=istat)
    ist = ist + istat
    DEALLOCATE(vtria3, STAT=istat)
    ist = ist + istat
    DEALLOCATE(vtria4, STAT=istat)
    ist = ist + istat
    DEALLOCATE(vtria5, STAT=istat)
    ist = ist + istat
    DEALLOCATE(vneig0, STAT=istat)
    ist = ist + istat
    DEALLOCATE(vneig1, STAT=istat)
    ist = ist + istat
    DEALLOCATE(vneig2, STAT=istat)
    ist = ist + istat
    DEALLOCATE(vneig3, STAT=istat)
    ist = ist + istat
    DEALLOCATE(vneig4, STAT=istat)
    ist = ist + istat
    DEALLOCATE(vneig5, STAT=istat)
    ist = ist + istat

    IF (ist>0 ) THEN
      WRITE (message_text,'(a)')    'Problem in deallocating local arrays'
      CALL finish ('mo_io_graph/input_graph', message_text)
    ENDIF

  END SUBROUTINE input_graph

  !EOC
  !-------------------------------------------------------------------------
  !
  !
  !>
  !! Creates files  with  cross reference tables of graph items to be read by.
  !!
  !! Creates files  with  cross reference tables of graph items to be read by
  !! the grid generator. For all graphs, the file <b>GRAPHMAP.*</b>
  !! is produced, where * denotes the grid level.
  !! The files are produced as sequential access, unformatted output.
  !!
  SUBROUTINE output_graph(k_root, sol, lplane)
    !

    INTEGER,       INTENT(in) :: k_root
    TYPE(t_spheres), INTENT(in) :: sol
    LOGICAL,       INTENT(in) :: lplane


    INTEGER :: i_grid_level                    ! grid level
    INTEGER :: istat
    INTEGER :: i, i_nt, i_ne, i_nv

    CHARACTER(len=filename_max) :: output

    INTEGER :: old_mode

    INTEGER :: ncid
    INTEGER :: dimids(2), start(2), count(2)

    INTEGER :: dim_ncell, dim_nvertex, dim_nedge
    INTEGER :: dim_ntriangles_per_edge, dim_nvertex_per_triangle, dim_nedges_per_vertex

    INTEGER :: varid1, varid2, varid3, varid4, varid5, varid6, varid7

#ifdef __SX__
    INTEGER :: iargc
    CHARACTER(len= 32) :: arg_str
#endif
    CHARACTER(len=256) :: command_line
    INTEGER :: command_line_len

    INTEGER, ALLOCATABLE :: kv1d(:)

    !-------------------------------------------------------------------------
    !BOC

    i_grid_level = sol%level

    IF (.NOT.lplane) THEN
       WRITE (output,'(a,i0,a,i2.2,a)') 'iconR', k_root, 'B', &
                                         i_grid_level, '-graph.nc'
    ELSE
       WRITE (output,'(a,i0,a,i2.2,a)') 'planR', k_root, 'B', &
                                         i_grid_level, '-graph.nc'
    ENDIF

    WRITE(message_text,'(a,a)') 'Write graphmap file: ', TRIM(output)
    CALL message ('', TRIM(message_text))

#ifdef __SX__
    command_line = ''
    DO i = 0, iargc()
      CALL getarg(i, arg_str)
      command_line = TRIM(command_line)//TRIM(arg_str)
    ENDDO
    command_line_len = LEN_TRIM(command_line)
#else
    CALL get_command(command_line, command_line_len, istat)
#endif

    !---------------------------------------------------------------------
    !
    CALL nf(nf_set_default_format(NF_FORMAT_64BIT, old_mode))
    CALL nf(nf_create(output, NF_CLOBBER, ncid))
    CALL nf(nf_set_fill(ncid, NF_NOFILL, old_mode))
    !
    !---------------------------------------------------------------------
    !
    ! Global attributes
    !
    CALL nf(nf_put_att_text    (ncid, NF_GLOBAL, 'title', 22, 'ICON graph description'))
    CALL nf(nf_put_att_text    (ncid, NF_GLOBAL, 'history', command_line_len, command_line))
    CALL nf(nf_put_att_text    (ncid, NF_GLOBAL, 'institution', 59, &
         'Max Planck Institute for Meteorology/Deutscher Wetterdienst'))
    CALL nf(nf_put_att_text    (ncid, NF_GLOBAL, 'source', 10, 'icon-dev'))
    CALL nf(nf_put_att_int     (ncid, NF_GLOBAL, 'grid_level', NF_INT, 1, i_grid_level))
    CALL nf(nf_put_att_int     (ncid, NF_GLOBAL, 'grid_root', NF_INT, 1, k_root))
    !
    ! Dimensions
    !
    i_nt = sol%no_triangles
    i_ne = sol%no_edges
    i_nv = sol%no_vertices
    !
    CALL nf(nf_def_dim(ncid, 'cell',   i_nt, dim_ncell))
    CALL nf(nf_def_dim(ncid, 'vertex', i_nv, dim_nvertex))
    CALL nf(nf_def_dim(ncid, 'edge',   i_ne, dim_nedge))
    !
    CALL nf(nf_def_dim(ncid, 'nt',        2, dim_ntriangles_per_edge))
    CALL nf(nf_def_dim(ncid, 'nv',        3, dim_nvertex_per_triangle))
    CALL nf(nf_def_dim(ncid, 'ne',        6, dim_nedges_per_vertex))
    !
    !---------------------------------------------------------------------
    !
    ! Grid variables
    !
    !---------------------------------------------------------------------
    !
    ! public graph information
    !
    ! Triangle vertex indices
    dimids = (/  dim_ncell, dim_nvertex_per_triangle /)
    CALL nf(nf_def_var(ncid, 'triangle_vertex_index', NF_INT, 2, dimids, varid1))
    CALL nf(nf_put_att_text(ncid, varid1, 'long_name', 21, 'triangle vertex index'))
    CALL nf(nf_put_att_text(ncid, varid1, 'cdi', 6, 'ignore'))
    !
    ! Triangle edge indices
    dimids = (/  dim_ncell, dim_nvertex_per_triangle /)
    CALL nf(nf_def_var(ncid, 'triangle_edge_index', NF_INT, 2, dimids, varid2))
    CALL nf(nf_put_att_text(ncid, varid2, 'long_name', 19, 'triangle edge index'))
    CALL nf(nf_put_att_text(ncid, varid2, 'cdi', 6, 'ignore'))
    !
    ! Edge vertex indices
    dimids = (/  dim_nedge, dim_ntriangles_per_edge /)
    CALL nf(nf_def_var(ncid, 'edge_vertex_index', NF_INT, 2, dimids, varid3))
    CALL nf(nf_put_att_text(ncid, varid3, 'long_name', 17, 'edge vertex index'))
    CALL nf(nf_put_att_text(ncid, varid3, 'cdi', 6, 'ignore'))
    !
    ! Edge triangle indices
    dimids = (/  dim_nedge, dim_ntriangles_per_edge /)
    CALL nf(nf_def_var(ncid, 'edge_triangle_index', NF_INT, 2, dimids, varid4))
    CALL nf(nf_put_att_text(ncid, varid4, 'long_name', 19, 'edge triangle index'))
    CALL nf(nf_put_att_text(ncid, varid4, 'cdi', 6, 'ignore'))
    !
    ! Vertex edge indices
    dimids = (/  dim_nvertex, dim_nedges_per_vertex /)
    CALL nf(nf_def_var(ncid, 'vertex_edge_index', NF_INT, 2, dimids, varid5))
    CALL nf(nf_put_att_text(ncid, varid5, 'long_name', 17, 'vertex edge index'))
    CALL nf(nf_put_att_text(ncid, varid5, 'cdi', 6, 'ignore'))
    !
    ! Vertex triangle indices
    dimids = (/  dim_nvertex, dim_nedges_per_vertex /)
    CALL nf(nf_def_var(ncid, 'vertex_triangle_index', NF_INT, 2, dimids, varid6))
    CALL nf(nf_put_att_text(ncid, varid6, 'long_name', 21, 'vertex triangle index'))
    CALL nf(nf_put_att_text(ncid, varid6, 'cdi', 6, 'ignore'))
    !
    ! Vertex neighbor indices
    dimids = (/ dim_nvertex, dim_nedges_per_vertex /)
    CALL nf(nf_def_var(ncid, 'vertex_neighbor_index', NF_INT, 2, dimids, varid7))
    CALL nf(nf_put_att_text(ncid, varid7, 'long_name', 21, 'vertex neighbor index'))
    CALL nf(nf_put_att_text(ncid, varid7, 'cdi', 6, 'ignore'))
    !
    CALL nf(nf_enddef(ncid))
    !
    !---------------------------------------------------------------------
    !
    ALLOCATE(kv1d(0:i_nt-1))
    !
    count = (/ i_nt, 1 /)
    DO i = 0, i_nt-1
      kv1d(i) = sol%ts(i)%vertex0%info%idx
    ENDDO
    start = (/ 1, 1 /)
    CALL nf(nf_put_vara_int(ncid, varid1, start, count, kv1d))
    DO i = 0, i_nt-1
      kv1d(i) = sol%ts(i)%vertex1%info%idx
    ENDDO
    start = (/ 1, 2 /)
    CALL nf(nf_put_vara_int(ncid, varid1, start, count, kv1d))
    DO i = 0, i_nt-1
      kv1d(i) = sol%ts(i)%vertex2%info%idx
    ENDDO
    start = (/ 1, 3 /)
    CALL nf(nf_put_vara_int(ncid, varid1, start, count, kv1d))
    !
    count = (/ i_nt, 1 /)
    DO i = 0, i_nt-1
      kv1d(i) = sol%ts(i)%edge0%info%idx
    ENDDO
    start = (/ 1, 1 /)
    CALL nf(nf_put_vara_int(ncid, varid2, start, count, kv1d))
    DO i = 0, i_nt-1
      kv1d(i) = sol%ts(i)%edge1%info%idx
    ENDDO
    start = (/ 1, 2 /)
    CALL nf(nf_put_vara_int(ncid, varid2, start, count, kv1d))
    DO i = 0, i_nt-1
      kv1d(i) = sol%ts(i)%edge2%info%idx
    ENDDO
    start = (/ 1, 3 /)
    CALL nf(nf_put_vara_int(ncid, varid2, start, count, kv1d))
    !
    DEALLOCATE(kv1d)
    !
    !---------------------------------------------------------------------
    !
    ALLOCATE(kv1d(0:i_ne-1))
    !
    count = (/ i_ne, 1 /)
    DO i = 0, i_ne-1
      kv1d(i) = sol%es(i)%vertex0%info%idx
    ENDDO
    start = (/ 1, 1 /)
    CALL nf(nf_put_vara_int(ncid, varid3, start, count, kv1d))
    DO i = 0, i_ne-1
      kv1d(i) = sol%es(i)%vertex1%info%idx
    ENDDO
    start = (/ 1, 2 /)
    CALL nf(nf_put_vara_int(ncid, varid3, start, count, kv1d))
    !
    count = (/ i_ne, 1 /)
    DO i = 0, i_ne-1
      kv1d(i) = sol%es(i)%triangle0%info%idx
    ENDDO
    start = (/ 1, 1 /)
    CALL nf(nf_put_vara_int(ncid, varid4, start, count, kv1d))
    DO i = 0, i_ne-1
      kv1d(i) = sol%es(i)%triangle1%info%idx
    ENDDO
    start = (/ 1, 2 /)
    CALL nf(nf_put_vara_int(ncid, varid4, start, count, kv1d))
    !
    DEALLOCATE(kv1d)
    !
    !---------------------------------------------------------------------
    !
    ALLOCATE(kv1d(0:i_nv-1))
    !
    count = (/ i_nv, 1 /)
    DO i = 0, i_nv-1
      kv1d(i) = sol%vs(i)%edge0%info%idx
    ENDDO
    start = (/ 1, 1 /)
    CALL nf(nf_put_vara_int(ncid, varid5, start, count, kv1d))
    DO i = 0, i_nv-1
      kv1d(i) = sol%vs(i)%edge1%info%idx
    ENDDO
    start = (/ 1, 2 /)
    CALL nf(nf_put_vara_int(ncid, varid5, start, count, kv1d))
    DO i = 0, i_nv-1
      kv1d(i) = sol%vs(i)%edge2%info%idx
    ENDDO
    start = (/ 1, 3 /)
    CALL nf(nf_put_vara_int(ncid, varid5, start, count, kv1d))
    DO i = 0, i_nv-1
      kv1d(i) = sol%vs(i)%edge3%info%idx
    ENDDO
    start = (/ 1, 4 /)
    CALL nf(nf_put_vara_int(ncid, varid5, start, count, kv1d))
    DO i = 0, i_nv-1
      kv1d(i) = sol%vs(i)%edge4%info%idx
    ENDDO
    start = (/ 1, 5 /)
    CALL nf(nf_put_vara_int(ncid, varid5, start, count, kv1d))
    DO i = 0, i_nv-1
      kv1d(i) = sol%vs(i)%edge5%info%idx
    ENDDO
    start = (/ 1, 6 /)
    CALL nf(nf_put_vara_int(ncid, varid5, start, count, kv1d))
    !
    count = (/ i_nv, 1 /)
    DO i = 0, i_nv-1
      kv1d(i) = sol%vs(i)%triangle0%info%idx
    ENDDO
    start = (/ 1, 1 /)
    CALL nf(nf_put_vara_int(ncid, varid6, start, count, kv1d))
    DO i = 0, i_nv-1
      kv1d(i) = sol%vs(i)%triangle1%info%idx
    ENDDO
    start = (/ 1, 2 /)
    CALL nf(nf_put_vara_int(ncid, varid6, start, count, kv1d))
    DO i = 0, i_nv-1
      kv1d(i) = sol%vs(i)%triangle2%info%idx
    ENDDO
    start = (/ 1, 3 /)
    CALL nf(nf_put_vara_int(ncid, varid6, start, count, kv1d))
    start = (/ 1, 4 /)
    DO i = 0, i_nv-1
      kv1d(i) = sol%vs(i)%triangle3%info%idx
    ENDDO
    CALL nf(nf_put_vara_int(ncid, varid6, start, count, kv1d))
    DO i = 0, i_nv-1
      kv1d(i) = sol%vs(i)%triangle4%info%idx
    ENDDO
    start = (/ 1, 5 /)
    CALL nf(nf_put_vara_int(ncid, varid6, start, count, kv1d))
    DO i = 0, i_nv-1
      kv1d(i) = sol%vs(i)%triangle5%info%idx
    ENDDO
    start = (/ 1, 6 /)
    CALL nf(nf_put_vara_int(ncid, varid6, start, count, kv1d))
    !
    count = (/ i_nv, 1 /)
    DO i = 0, i_nv-1
      kv1d(i) = sol%vs(i)%neighbor0%info%idx
    ENDDO
    start = (/ 1, 1 /)
    CALL nf(nf_put_vara_int(ncid, varid7, start, count, kv1d))
    DO i = 0, i_nv-1
      kv1d(i) = sol%vs(i)%neighbor1%info%idx
    ENDDO
    start = (/ 1, 2 /)
    CALL nf(nf_put_vara_int(ncid, varid7, start, count, kv1d))
    DO i = 0, i_nv-1
      kv1d(i) = sol%vs(i)%neighbor2%info%idx
    ENDDO
    start = (/ 1, 3 /)
    CALL nf(nf_put_vara_int(ncid, varid7, start, count, kv1d))
    DO i = 0, i_nv-1
      kv1d(i) = sol%vs(i)%neighbor3%info%idx
    ENDDO
    start = (/ 1, 4 /)
    CALL nf(nf_put_vara_int(ncid, varid7, start, count, kv1d))
    DO i = 0, i_nv-1
      kv1d(i) = sol%vs(i)%neighbor4%info%idx
    ENDDO
    start = (/ 1, 5 /)
    CALL nf(nf_put_vara_int(ncid, varid7, start, count, kv1d))
    DO i = 0, i_nv-1
      kv1d(i) = sol%vs(i)%neighbor5%info%idx
    ENDDO
    start = (/ 1, 6 /)
    CALL nf(nf_put_vara_int(ncid, varid7, start, count, kv1d))
    !
    DEALLOCATE(kv1d)
    !
    !---------------------------------------------------------------------
    !
    CALL nf(nf_close(ncid))
    !
    !---------------------------------------------------------------------

  END SUBROUTINE output_graph

  !EOC

  SUBROUTINE nf(status)
    INTEGER, INTENT(in) :: status

    IF (status /= nf_noerr) THEN
      CALL finish('mo_io_graph/nf: netCDF error', nf_strerror(status))
    ENDIF

  END SUBROUTINE nf

END MODULE mo_io_graph


