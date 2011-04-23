!-------------------------------------------------------------------------------------
! mo_create_ocean_grid
!>
!! A collection of grid tools assocoated with the ocean grid
!!
!! NOTE : negative elevations are below sea level (sea depth)
!!
!! @author Leonidas Linardakis, MPI-M
!!
!! @par Revision History
!!   First implementation by Leonidas Linardakis, MPI-M, 2009-12-4
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
MODULE mo_zlayers
#include "grid_definitions.inc"
  !-------------------------------------------------------------------------
  USE mo_kind,                 ONLY: wp
  USE mo_exception,            ONLY: finish
! USE mo_io_units,             ONLY: filename_max, nnml
! USE mo_namelist,             ONLY: position_nml, open_nml, positioned
  USE mo_local_grid

  IMPLICIT NONE

  PRIVATE

  ! !VERSION CONTROL:
  CHARACTER(LEN=*), PARAMETER :: version = '$Id$'

  !-------------------------------------------------------------------------
  PUBLIC :: read_zlayers
  PUBLIC :: get_grid_zlayers
  PUBLIC :: count_dry_wet_cells
  PUBLIC :: count_active_layer_cells
  !-------------------------------------------------------------------------

CONTAINS



  !-------------------------------------------------------------------------
  ! get_grid_zlayers
  !>
  !!
  SUBROUTINE get_grid_zlayers(in_grid_id, vertical)
    INTEGER, INTENT(in) :: in_grid_id
    TYPE(t_vertical_ocean_structure), TARGET, INTENT(inout)  :: vertical

    TYPE(t_grid), POINTER :: in_grid
    TYPE(t_grid_cells), POINTER :: cells
    TYPE(t_grid_edges), POINTER :: edges
    REAL(wp) :: cell_depth
    INTEGER :: no_of_cells,no_of_edges
    INTEGER :: no_of_zlayers
    INTEGER :: i,j,cell1,cell2
    INTEGER :: totalcells

    in_grid => get_grid(in_grid_id)
    cells => in_grid%cells
    no_of_cells = cells%no_of_existcells
    edges => in_grid%edges
    no_of_edges = edges%no_of_existedges

    no_of_zlayers = vertical%no_of_zlayers
    in_grid%vertical_structure => vertical

    !----------------------------------------
    ! get the number of layers for each cell
    totalcells = 0
    DO i=1,no_of_cells
      cell_depth = -cells%elevation(i)
      IF (cell_depth >= vertical%layer_middle(no_of_zlayers)) THEN
         cells%no_of_zlayers(i) = no_of_zlayers
!           WRITE(*,*) i, cell_depth, no_of_zlayers, vertical%layer_middle(no_of_zlayers)
      ELSE
        DO j=2,no_of_zlayers
          IF (cell_depth < vertical%layer_middle(j)) THEN
            cells%no_of_zlayers(i) = j-1
!            WRITE(*,*) i, cell_depth, j-1, vertical%layer_middle(j-1), vertical%layer_middle(j)
            EXIT
          ENDIF
        ENDDO
      ENDIF
      ! write(*,*) i,' cell depth=', cell_depth, cells%no_of_zlayers(i)
      totalcells = totalcells + cells%no_of_zlayers(i)
    ENDDO ! i=1,no_of_cells
    write(*,*) 'get_grid_zlayers: totalcells=', totalcells
    !----------------------------------------
    ! get the number of layers for each edge
    DO i=1,no_of_edges
      cell1 = edges%get_cell_index(i,1)
      cell2 = edges%get_cell_index(i,2)
      IF (cell1 == 0 .OR. cell2 == 0) THEN
         edges%no_of_zlayers(i) = 0
         CYCLE
      ENDIF
      edges%no_of_zlayers(i) = MIN&
        & (cells%no_of_zlayers(cell1), cells%no_of_zlayers(cell2))
    ENDDO ! i=1,no_of_edges

  END SUBROUTINE get_grid_zlayers
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE read_zlayers(vertical, filename)

    TYPE(t_vertical_ocean_structure) :: vertical
    CHARACTER(LEN=*) :: filename

    INTEGER :: no_of_zlayers
    INTEGER :: i,j

!    CHARACTER(LEN=filename_max) :: verticalfilename

    WRITE (*,*) "--------------------------------"
    WRITE (*,*) "Read vertical layers"

    OPEN (500, FILE = TRIM(filename),STATUS = 'OLD')
    READ (500, *) no_of_zlayers
    WRITE (*,*) "no_of_zlayers:",no_of_zlayers
    CALL allocate_vertical_zlayers(vertical,no_of_zlayers)

    DO i=1,no_of_zlayers
      READ (500, *) j,vertical%layer_thicknes(i)
      WRITE (*,*) i,j, ":", vertical%layer_thicknes(i)
      IF (j /= i) &
        CALL finish ('read_zlayers','wrong layer index')
    ENDDO
    CLOSE(500)

    ! compute vertical parameters
    vertical%layer_bed(1) = vertical%layer_thicknes(1)
    vertical%layer_middle(1) = vertical%layer_thicknes(1) * 0.5_wp
    DO i=2,no_of_zlayers
      vertical%layer_bed(i) = vertical%layer_bed(i-1) + vertical%layer_thicknes(i)
      vertical%layer_middle(i) = vertical%layer_bed(i-1) + vertical%layer_thicknes(i) * 0.5_wp
    ENDDO

  END SUBROUTINE read_zlayers
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE count_dry_wet_cells(in_grid_id)
    INTEGER, INTENT(in) :: in_grid_id

    TYPE(t_grid), POINTER :: in_grid
    TYPE(t_grid_cells), POINTER :: cells

    TYPE(t_vertical_ocean_structure), POINTER :: vertical
    INTEGER :: no_of_zlayers, no_of_cells
    INTEGER :: i,j,total_active
    INTEGER, POINTER :: active_cells_in_layer(:)

    in_grid => get_grid(in_grid_id)
    cells => in_grid%cells
    no_of_cells = cells%no_of_existcells
    vertical => in_grid%vertical_structure
    no_of_zlayers = vertical%no_of_zlayers

    ! compute  number of dry/wet cells per layer
    ALLOCATE(active_cells_in_layer(no_of_zlayers), STAT=i)
    IF (i /= 0) THEN
      CALL finish ('count_dry_wet_cells','allocating active_cells_in_layer')
    ENDIF
    active_cells_in_layer(:) = 0
    DO i=1,no_of_cells
      active_cells_in_layer(cells%no_of_zlayers(i)) = &
        & active_cells_in_layer(cells%no_of_zlayers(i)) + 1
    ENDDO
    ! add active on all the previous ones
    DO j=no_of_zlayers-1,1,-1
      active_cells_in_layer(j) = active_cells_in_layer(j) + active_cells_in_layer(j+1)
    ENDDO

    total_active = 0
    DO j=1,no_of_zlayers
      WRITE(0,'(a10,i2, a6,f8.2, a5,i9, a5,i9, f6.2 )') &
        "Level:",j, " Depth:",vertical%layer_bed(j),&
        "Wet=",active_cells_in_layer(j),&
        "Dry=",no_of_cells-active_cells_in_layer(j),&
        real(no_of_cells-active_cells_in_layer(j))* 100.0 / real(no_of_cells)

        total_active = total_active + active_cells_in_layer(j)
    ENDDO
    WRITE(0,'(a10,i2, a5,i9, a5,i9, f6.2 )') &
      "Total:",j,"Wet=",total_active,&
      "Dry=",no_of_cells*no_of_zlayers-total_active,&
      real(no_of_cells*no_of_zlayers-total_active)* 100.0 / real(no_of_cells*no_of_zlayers)

    DEALLOCATE(active_cells_in_layer)

  END SUBROUTINE count_dry_wet_cells
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE count_active_layer_cells(in_grid_id)
    INTEGER, INTENT(in) :: in_grid_id

    TYPE(t_grid_cells), POINTER :: cells
    INTEGER :: no_of_cells
    INTEGER :: i
    INTEGER :: totalcells

    cells => get_cells(in_grid_id)
    no_of_cells = cells%no_of_existcells

    !----------------------------------------
    totalcells = 0
    DO i=1,no_of_cells
      ! write(*,*) i,' cell depth=', cell_depth, cells%no_of_zlayers(i)
      totalcells = totalcells + cells%no_of_zlayers(i)
    ENDDO ! i=1,no_of_cells
    write(*,*) 'count_active_layer_cells: totalcells=', totalcells
    !----------------------------------------
  END SUBROUTINE count_active_layer_cells

END MODULE mo_zlayers

