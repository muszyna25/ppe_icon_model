!>
!! Definition of the ICON grid for coupling
!!
!! <Description>
!! 
!! @author Rene Redler, Max-Planck Institute for Meteorology, Germany
!!
!! $Id:$
!!
!! @par Revision History
!! first implementation by Rene Redler (2010-02-13)
!!
!! @par License
!! The use of ICON is hereby granted free of charge for an unlimited time,
!! provided the following rules are accepted and applied:
!! &ltol>
!! &ltli> You may use or modify this code for your own non commercial and non
!!    violent purposes.
!! &ltli> The code may not be re-distributed without the consent of the authors.
!! &ltli> The copyright notice and statement of authorship must appear in all
!!    copies.
!! &ltli> You accept the warranty conditions (see WARRANTY).
!! &ltli> In case you intend to use the code commercially, we oblige you to sign
!!    an according license agreement with DWD and MPI-M.
!! </ol>
!!
!! @par Warranty
!!
!! This code has been tested up to a certain level. Defects and weaknesses,
!! which may be included in the code, do not establish any warranties by the
!! authors.
!! The authors do not make any warranty, express or implied, or assume any
!! liability or responsibility for the use, acquisition or application of this
!! software.
!!
MODULE mo_icon_cpl_def_grid

  USE mo_icon_cpl, ONLY : ICON_comm, t_grid, grids, &
   &                      nbr_active_grids, nbr_ICON_grids


  IMPLICIT NONE

  PRIVATE

  TYPE(t_grid), POINTER  :: gptr

  ! Return code for error handling

  INTEGER                :: ierr

  PUBLIC :: ICON_cpl_def_grid

CONTAINS

  SUBROUTINE ICON_cpl_def_grid ( comp_id, grid_shape, grid_glob_index, grid_id, ierror )

    INTEGER, INTENT(in)  :: comp_id            !<  component id
    INTEGER, INTENT(in)  :: grid_shape(2)      !<  shape of index array
    INTEGER, INTENT(in)  :: grid_glob_index & 
                        (grid_shape(1):grid_shape(2)) !<  index array

    INTEGER, INTENT(out) :: grid_id            !<  grid ID

    INTEGER, INTENT(out) :: ierror             !<  returned error code

    ! -------------------------------------------------------------------
    ! Initialise variables
    ! -------------------------------------------------------------------

    ierror = 0

    ! -------------------------------------------------------------------
    ! Initialise grid
    ! -------------------------------------------------------------------

    DO grid_id = 1, nbr_ICON_grids
       IF ( .NOT. grids(grid_id)%l_grid_status ) EXIT
    ENDDO

    IF ( grid_id > nbr_ICON_grids ) THEN
       WRITE ( * , * ) 'number of requested grids exceeds maximum of ', nbr_ICON_grids
#ifndef NOMPI
       CALL MPI_Abort ( ICON_comm, 1, ierr )
#endif
    ENDIF

    gptr => grids(grid_id)

    ! -------------------------------------------------------------------
    ! Store grid parameters and global index list
    ! -------------------------------------------------------------------

    gptr%comp_id       = comp_id
    gptr%l_grid_status = .TRUE.

    ALLOCATE ( gptr%grid_glob_index(grid_shape(1):grid_shape(2)), STAT = ierr )
    IF ( ierr > 0 ) THEN
       WRITE ( * , * ) ' Error allocating Grids '
#ifndef NOMPI
       CALL MPI_Abort ( ICON_comm, 1, ierr )
#endif
    ENDIF

    gptr%grid_shape         = grid_shape
    gptr%grid_glob_index(grid_shape(1):grid_shape(2)) = &
      &  grid_glob_index(grid_shape(1):grid_shape(2))

    nbr_active_grids = nbr_active_grids + 1

  END SUBROUTINE ICON_cpl_def_grid

END MODULE mo_icon_cpl_def_grid
