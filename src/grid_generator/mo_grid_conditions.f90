!-------------------------------------------------------------------------------------
! mo_grid_conditions
!>
!! Set od geometric conditions for defining a sub-grid
!!
!! @author Leonidas Linardakis, MPI-M
!!
!! @par Revision History
!!   First implementation by Leonidas Linardakis, MPI-M, 2010-03-01
!!
!! @par Copyright
!! 2002-2010 by DWD and MPI-M
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
!-------------------------------------------------------------------------------------
MODULE mo_grid_conditions
#include "grid_definitions.inc"
  !-------------------------------------------------------------------------
  USE mo_kind,               ONLY: wp
  USE mo_math_constants,     ONLY: rad2deg,pi
  USE mo_exception,          ONLY: message_text, message, finish
  USE mo_io_units,           ONLY: nnml, filename_max
  USE mo_namelist,           ONLY: position_nml, open_nml, positioned
  USE mo_local_grid,         ONLY: t_grid_cells, new_grid, delete_grid, &
    & cut_off_grid, set_grid_creation, get_cells, t_integer_list!, undefined
  USE mo_base_geometry,      ONLY: t_cartesian_coordinates, t_geographical_coordinates, gc2cc, &
    & arc_length
!  USE mo_local_grid_geometry,ONLY: geographical_to_cartesian
  USE mo_grid_toolbox,       ONLY: smooth_boundaryfrom_cell_list, &
    & get_grid_from_cell_list
  USE mo_io_local_grid,      ONLY: read_netcdf_grid, write_netcdf_grid
  USE mo_timer

  IMPLICIT NONE

  PRIVATE

  ! !VERSION CONTROL:
  CHARACTER(LEN=*), PARAMETER :: version = '$Id$'

  PUBLIC :: cut_local_grid, cut_conditional_grid, read_grid_conditions, get_conditional_cells
  !----------------------------------------

  ! !DEFINE PARAMETERS:
  INTEGER, PARAMETER ::  max_no_of_conditions = 40
  INTEGER, PARAMETER ::  rectangle_shape = 1
  INTEGER, PARAMETER ::  circle_shape = 2

  !-------------------------------------------------------------------------
  INTEGER :: no_of_conditions = 0
  REAL(wp) :: patch_center_x(max_no_of_conditions), patch_center_y(max_no_of_conditions)
  TYPE(t_geographical_coordinates) :: patch_center_geocoord(max_no_of_conditions)
  TYPE(t_cartesian_coordinates)    :: patch_center_cartesian(max_no_of_conditions)
  REAL(wp) :: rectangle_xradious(max_no_of_conditions), rectangle_yradious(max_no_of_conditions)
  REAL(wp) :: circle_radious(max_no_of_conditions)
  INTEGER :: patch_shape(max_no_of_conditions)

  CHARACTER(LEN=filename_max) :: input_file, output_file
  !-------------------------------------------------------------------------


CONTAINS

  !-------------------------------------------------------------------------
  !   FUNCTION read_grid_conditions(param_file_name)
  !>
  !! reads the parameters for the cutting or refining a grid. Private
  FUNCTION read_grid_conditions(param_file_name) result(no_of_read_conditions)

    CHARACTER(LEN=*), INTENT(in) :: param_file_name
    INTEGER :: no_of_read_conditions

    INTEGER :: i_status

    NAMELIST /grid_geometry_conditions/ input_file, output_file, &
      & no_of_conditions, patch_shape,   &
      & patch_center_x, patch_center_y,         &
      & rectangle_xradious, rectangle_yradious, &
      & circle_radious

    ! set default values
    input_file=''
    output_file = ''
    no_of_conditions = 0
    patch_shape=0
    rectangle_xradious = 0.0_WP
    rectangle_yradious = 0.0_WP
    circle_radious = 0.0_WP
    patch_center_x = 0.0_WP
    patch_center_y = 0.0_WP
    no_of_conditions = 0
    no_of_read_conditions = 0
    
    ! read namelist
    CALL open_nml(param_file_name)
    CALL position_nml('grid_geometry_conditions',STATUS=i_status)
    IF (i_status == positioned) THEN
      READ (nnml,grid_geometry_conditions)
    ELSE
      RETURN
!       WRITE(message_text,'(a,a)') " File", param_file_name, " not POSITIONED"
!       CALL finish ('read_grid_conditions', message_text)
    ENDIF
    CLOSE(nnml)

    WRITE(message_text,'(a)') "===================================="
    CALL message ('', TRIM(message_text))
    WRITE(message_text,'(a)') "----- grid_geometry_conditions -----"
    CALL message ('', TRIM(message_text))
    IF (input_file /= '') THEN
      WRITE(message_text,'(a,a)') "inputFile=",  TRIM(input_file)
      CALL message ('', TRIM(message_text))
      WRITE(message_text,'(a,a)') "outputFile=", TRIM(output_file)
      CALL message ('', TRIM(message_text))
    ENDIF
    WRITE(message_text,'(a,i2)') "no_of_conditions=", no_of_conditions
    CALL message ('', TRIM(message_text))
    WRITE(message_text,'(a,40i2)') "patch_shapes=", patch_shape(1:no_of_conditions)
    CALL message ('', TRIM(message_text))
    WRITE(message_text,'(a,40F6.2)') "patch_center_x=", patch_center_x(1:no_of_conditions)
    CALL message ('', TRIM(message_text))
    WRITE(message_text,'(a,40F6.2)') "patch_center_y=", patch_center_y(1:no_of_conditions)
    CALL message ('', TRIM(message_text))
    WRITE(message_text,'(a,40F6.2)') "rectangle_xradious=", &
      & rectangle_xradious(1:no_of_conditions)
    CALL message ('', TRIM(message_text))
    WRITE(message_text,'(a,40F6.2)') "rectangle_yradious=", &
      & rectangle_yradious(1:no_of_conditions)
    CALL message ('', TRIM(message_text))
    WRITE(message_text,'(a,40F16.3)') "circle_radious=", &
      & circle_radious(1:no_of_conditions)
    CALL message ('', TRIM(message_text))
    WRITE(message_text,'(a)') "===================================="
    CALL message ('', TRIM(message_text))

    CALL ref_set_param()
    no_of_read_conditions = no_of_conditions
    
  END FUNCTION read_grid_conditions
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !   SUBROUTINE cut_local_grid(param_file_name)
  !>
  !! Main routine for cutting-off a sub-grid. Public
  SUBROUTINE cut_local_grid(param_file_name)

    CHARACTER(LEN=*), INTENT(in) :: param_file_name

    INTEGER :: in_grid_id, out_grid_id, tmp

    tmp=read_grid_conditions(param_file_name)

    in_grid_id = new_grid()
    CALL read_netcdf_grid(in_grid_id, input_file)

    out_grid_id = cut_conditional_grid(in_grid_id)
    IF (out_grid_id == in_grid_id .OR. out_grid_id < 1)&
      & CALL finish ('grid_cutLocalGrid', 'Failed to cut grid')

    CALL write_netcdf_grid(out_grid_id, output_file)
    CALL delete_grid(in_grid_id)
    CALL delete_grid(out_grid_id)

  END SUBROUTINE cut_local_grid
  !-------------------------------------------------------------------------



  !-------------------------------------------------------------------------
  !   SUBROUTINE cut_local_grid(param_file_name)
  !>
  !! Cuts-off the cut_grid_id from in_grid_id. Private
  INTEGER FUNCTION cut_conditional_grid(in_grid_id, param_file) result(cut_grid_id)
    INTEGER, INTENT(in) :: in_grid_id
    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: param_file
!     INTEGER, INTENT(out):: cut_grid_id
!     INTEGER, INTENT(out):: cut_status

    TYPE(t_integer_list)  :: cut_cell_list, smooth_cell_list
    INTEGER :: tmp

!    cut_status = -
    IF (PRESENT(param_file)) tmp=read_grid_conditions(param_file)
    cut_grid_id = in_grid_id

    ! IF (patch_shape == undefined) RETURN
    IF (no_of_conditions  < 1) RETURN

    CALL get_conditional_cells(in_grid_id, cut_cell_list)
    CALL smooth_boundaryfrom_cell_list(in_grid_id, cut_cell_list, smooth_cell_list)
    cut_grid_id = get_grid_from_cell_list(in_grid_id, smooth_cell_list)
    CALL set_grid_creation(cut_grid_id, cut_off_grid)

    DEALLOCATE(cut_cell_list%value)
    DEALLOCATE(smooth_cell_list%value)

    ! cut_status = 0

  END FUNCTION cut_conditional_grid
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !   SUBROUTINE ref_set_param
  !>
  !! Set the actual parameters for cutting and refining a grid. Private
  !! Changes degrees to rads, km to unit-sphere length
  SUBROUTINE ref_set_param()

    INTEGER :: i

    patch_center_x(1:no_of_conditions) = patch_center_x(1:no_of_conditions) / rad2deg
    patch_center_y(1:no_of_conditions) = patch_center_y(1:no_of_conditions) / rad2deg
    patch_center_geocoord(1:no_of_conditions)%lon = patch_center_x(1:no_of_conditions)
    patch_center_geocoord(1:no_of_conditions)%lat = patch_center_y(1:no_of_conditions)
    DO i=1,no_of_conditions
      patch_center_cartesian(i)    = gc2cc(patch_center_geocoord(i))
    ENDDO

    rectangle_xradious(1:no_of_conditions) = rectangle_xradious(1:no_of_conditions) / rad2deg
    rectangle_yradious(1:no_of_conditions) = rectangle_yradious(1:no_of_conditions) / rad2deg

    circle_radious(1:no_of_conditions) = circle_radious(1:no_of_conditions) / rad2deg

    WRITE(message_text,'(a,40F9.6)') "circle_radious=", circle_radious(1:no_of_conditions)
    CALL message ('', TRIM(message_text))

  END SUBROUTINE ref_set_param
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !   SUBROUTINE get_conditional_cells(in_grid_id, cell_list)
  !>
  !! Returns the cells in cell_list that are contained into an area. Private
  SUBROUTINE get_conditional_cells(in_grid_id, t_cell_list)

    INTEGER, INTENT(in):: in_grid_id
    TYPE(t_integer_list) :: t_cell_list

    TYPE(t_grid_cells), POINTER :: cells
    INTEGER :: no_of_input_cells,no_of_list_cells, cell_index
    ! TYPE(cartesian_coordinates), POINTER :: cartesian_cellcenter(:)
    TYPE(t_cartesian_coordinates) :: c_cellcenter
    REAL(wp) :: cell_lon,cell_lat,x_dist,y_dist

    INTEGER ::  istat, condition_index
    !-------------------------------------------------------------------------

    cells => get_cells(in_grid_id)
    no_of_input_cells = cells%no_of_allocatedcells

    ! cartesian_cellcenter => cells%cartesian_center
    !CALL geographical_to_cartesian(cells%center,no_of_input_cells, cartesian_cellcenter)

    ALLOCATE (t_cell_list%value(no_of_input_cells),stat=istat)
    IF (istat >0) THEN
      CALL finish ('get_conditional_cells', 'Problem in allocating cell_list')
    ENDIF

    no_of_list_cells=0
    DO cell_index = 1, no_of_input_cells
      cell_lon = cells%center(cell_index)%lon
      cell_lat = cells%center(cell_index)%lat
      c_cellcenter = cells%cartesian_center(cell_index)

      DO condition_index = 1, no_of_conditions

        SELECT CASE(patch_shape(condition_index))

        CASE(circle_shape)
          IF (arc_length(c_cellcenter, patch_center_cartesian(condition_index)) &
            & <= circle_radious(condition_index)) THEN
            ! cell center in circle
            ! add cell and exit loop
            no_of_list_cells=no_of_list_cells+1
            t_cell_list%value(no_of_list_cells) = cell_index
            EXIT
          ENDIF

        CASE(rectangle_shape)
          y_dist =  ABS(cell_lat - patch_center_y(condition_index))
          IF ( y_dist <= rectangle_yradious(condition_index)) THEN
            !  check the x distance
            x_dist = cell_lon - patch_center_x(condition_index)
            IF ( x_dist < -pi) THEN
              x_dist = x_dist + 2._wp*pi
            ELSE IF ( x_dist > pi) THEN
              x_dist = x_dist - 2._wp
            ENDIF
            IF (ABS(x_dist) <= rectangle_xradious(condition_index)) THEN
              ! cell center in rectangle
              ! add cell and exit loop
              no_of_list_cells=no_of_list_cells+1
              t_cell_list%value(no_of_list_cells) = cell_index
              EXIT
            ENDIF
          ENDIF

        CASE default
          CALL finish ('get_conditional_cells', 'Unkown patch_shape')

        END SELECT

      ENDDO !condition_index = 1, no_of_conditions
    ENDDO ! cell_index = 1, no_of_input_cells

    !  allocate cellList adjusted to the no_of_list_cells
    ! IF (ALLOCATED(cell_list)) DEALLOCATE (cell_list,stat=istat)

    ! fill cellList
    t_cell_list%list_size = no_of_list_cells
    ! cell_list%value(:) = tmp_cell_list%value(1:no_of_list_cells)

    ! clean-up
    ! DEALLOCATE (tmp_cell_list%value,stat=istat)
    RETURN

  END SUBROUTINE get_conditional_cells
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !   SUBROUTINE cells_center_in_rectangle(cells, tagged_cell_list)
  !>
  !! Returns in tagged_cell_list the cells whose center is in a rectangle. Private
  !   SUBROUTINE cells_center_in_rectangle(cells_id, tagged_cell_list)
  !     INTEGER, INTENT(in):: cells_id
  !     TYPE(integer_list) :: tagged_cell_list
  !
  !     TYPE(grid_cells), POINTER :: cells
  !     INTEGER :: no_of_input_cells,no_of_list_cells
  !     INTEGER :: cell_index,istat
  !     REAL(wp) :: cell_lon,cell_lat,x_dist,y_dist
  !
  !     cells => get_cells(cells_id)
  !     no_of_input_cells = cells%no_of_cells
  !
  !     WRITE(message_text,'(a)') "===================================="
  !     CALL message ('', TRIM(message_text))
  !     WRITE(message_text,'(a)')  "get cells in rectangle"
  !     CALL message ('', TRIM(message_text))
  !
  !     ALLOCATE (tagged_cell_list%value(no_of_input_cells),stat=istat)
  !     IF (istat >0) THEN
  !       CALL finish ('grid_getConditionalCellIndex', 'Problem in allocating tmpCellList')
  !     ENDIF
  !
  !     !  fill tagged_cell_list
  !     no_of_list_cells=0
  !     DO cell_index = 1, no_of_input_cells
  !       cell_lon = cells%center(cell_index)%lon
  !       cell_lat = cells%center(cell_index)%lat
  !
  !       y_dist =  ABS(cell_lat - patch_center_y)
  !       IF ( y_dist <= rectangle_yradious) THEN
  !         !  check the x distance
  !         x_dist = cell_lon - patch_center_x
  !         IF ( x_dist < -pi) THEN
  !           x_dist = x_dist + 2._wp*pi
  !         ELSE IF ( x_dist > pi) THEN
  !           x_dist = x_dist - 2._wp
  !         ENDIF
  !         IF (ABS(x_dist) <= rectangle_xradious) THEN
  !           ! cell center in rectangle
  !           no_of_list_cells=no_of_list_cells+1
  !           tagged_cell_list%value(no_of_list_cells) = cell_index
  !         ENDIF
  !       ENDIF
  !     ENDDO
  !
  !     tagged_cell_list%list_size = no_of_list_cells
  !     WRITE(message_text,'(a, i9, i9)') "input/output Cells:", no_of_input_cells, no_of_list_cells
  !     CALL message ('', TRIM(message_text))
  !
  !     RETURN
  !
  !   END SUBROUTINE cells_center_in_rectangle
  !   !-------------------------------------------------------------------------
  !
  !   !-------------------------------------------------------------------------
  !   !   SUBROUTINE cells_center_in_circle(cells, tagged_cell_list)
  !   !>
  !   !! Returns in tagged_cell_list the cells whose center is in a rectangle. Private
  !   SUBROUTINE cells_center_in_circle(cells_id, tagged_cell_list)
  !     INTEGER, INTENT(in):: cells_id
  !     TYPE(integer_list) :: tagged_cell_list
  !
  !     TYPE(grid_cells), POINTER :: cells
  !     INTEGER :: no_of_input_cells,no_of_list_cells
  !     TYPE(cartesian_coordinates), POINTER :: cartesian_cellcenter(:)
  !
  !     INTEGER :: cell_index,istat
  !
  !     cells => get_cells(cells_id)
  !     no_of_input_cells = cells%no_of_cells
  !
  !     WRITE(message_text,'(a)') "===================================="
  !     CALL message ('', TRIM(message_text))
  !     WRITE(message_text,'(a)')  "get cells in circle"
  !     CALL message ('', TRIM(message_text))
  !
  !     CALL geographical_to_cartesian(cells%center,no_of_input_cells, cartesian_cellcenter)
  !
  !     ALLOCATE (tagged_cell_list%value(no_of_input_cells),stat=istat)
  !     IF (istat >0) THEN
  !       CALL finish ('grid_getConditionalCellIndex', 'Problem in allocating tmpCellList')
  !     ENDIF
  !
  !     !  fill tagged_cell_list
  !     no_of_list_cells=0
  !     DO cell_index = 1, no_of_input_cells
  !       IF (arc_length(cartesian_cellcenter(cell_index), patch_center_cartesian) &
  !         & <= circle_radious) THEN
  !         ! cell center in circle
  !         no_of_list_cells=no_of_list_cells+1
  !         tagged_cell_list%value(no_of_list_cells) = cell_index
  !       ENDIF
  !     ENDDO
  !
  !     tagged_cell_list%list_size = no_of_list_cells
  !     WRITE(message_text,'(a, i9, i9)') "input/output Cells:", no_of_input_cells, no_of_list_cells
  !     CALL message ('', TRIM(message_text))
  !
  !     DEALLOCATE(cartesian_cellcenter)
  !
  !     RETURN
  !
  !   END SUBROUTINE cells_center_in_circle
  !-------------------------------------------------------------------------



END MODULE mo_grid_conditions

