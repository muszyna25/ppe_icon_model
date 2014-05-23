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
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_icon_cpl_def_grid

  USE mo_icon_cpl, ONLY : ICON_comm, t_grid, grids, &
   &                      nbr_active_grids, nbr_ICON_grids, &
   &                      ICON_local_rank


  IMPLICIT NONE

  PRIVATE

  TYPE(t_grid), POINTER  :: gptr

  ! Return code for error handling

  INTEGER                :: ierr

  PUBLIC :: ICON_cpl_def_grid, ICON_cpl_def_location

CONTAINS

  SUBROUTINE ICON_cpl_def_grid ( grid_shape, grid_glob_index, grid_id, ierror )

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
       ierror = 1
#ifndef NOMPI
       CALL MPI_Abort ( ICON_comm, 1, ierr )
#endif
    ENDIF

    gptr => grids(grid_id)

    ! -------------------------------------------------------------------
    ! Store grid parameters and global index list
    ! -------------------------------------------------------------------

    gptr%l_grid_status = .TRUE.

    ALLOCATE ( gptr%grid_glob_index(grid_shape(1):grid_shape(2)), STAT = ierr )
    IF ( ierr > 0 ) THEN
       WRITE ( * , * ) ' Error allocating Grids '
       ierror = 1
#ifndef NOMPI
       CALL MPI_Abort ( ICON_comm, 1, ierr )
#endif
    ENDIF

    gptr%grid_shape         = grid_shape
    gptr%grid_glob_index(grid_shape(1):grid_shape(2)) = &
      &  grid_glob_index(grid_shape(1):grid_shape(2))

    nbr_active_grids = nbr_active_grids + 1

  END SUBROUTINE ICON_cpl_def_grid

  ! --------------------------------------------------------------------

  SUBROUTINE ICON_cpl_def_location ( grid_id, grid_shape, glob_index_rank, this_owner, ierror )

    INTEGER, INTENT(in)  :: grid_id            !<  grid ID
    INTEGER, INTENT(in)  :: grid_shape(2)      !<  shape of index array
    INTEGER, INTENT(in)  :: glob_index_rank & 
                        (grid_shape(1):grid_shape(2)) !<  list of ranks for each vertex
    INTEGER, INTENT(in)   :: this_owner

    INTEGER, INTENT(out) :: ierror             !<  returned error code

    ! -------------------------------------------------------------------
    ! Initialise variables
    ! -------------------------------------------------------------------

    ierror = 0

    ! -------------------------------------------------------------------
    ! Initialise grid
    ! -------------------------------------------------------------------

    IF ( .NOT. ASSOCIATED(grids) ) THEN
       WRITE ( * , * ) 'First call ICON_cpl_def_grid.'
       ierror = 1
#ifndef NOMPI
       CALL MPI_Abort ( ICON_comm, 1, ierr )
#endif
    ENDIF

    IF ( grid_id > nbr_ICON_grids .OR. .NOT. grids(grid_id)%l_grid_status ) THEN
       WRITE ( * , * ) 'Invalid grid ID ', grid_id
       WRITE ( * , * ) 'First call ICON_cpl_def_grid.'
       ierror = 1
#ifndef NOMPI
       CALL MPI_Abort ( ICON_comm, 1, ierr )
#endif
    ENDIF

    gptr => grids(grid_id)

    IF ( gptr%grid_shape(1) /= grid_shape(1) .OR. &
         gptr%grid_shape(2) /= grid_shape(2) ) THEN
       WRITE ( * , * ) 'Inconsistent grid shape ', grid_shape(1), grid_shape(2)
       WRITE ( * , * ) 'Expected: ', gptr%grid_shape(1), gptr%grid_shape(2)
       ierror = 1
#ifndef NOMPI
       CALL MPI_Abort ( ICON_comm, 1, ierr )
#endif
    ENDIF


    ! -------------------------------------------------------------------
    ! Check if this owner matches the ICON_local_rank
    ! -------------------------------------------------------------------
    IF (ICON_local_rank /= this_owner) THEN
      WRITE(0,*) "ICON_local_rank /= this_owner", ICON_local_rank, this_owner
#ifndef NOMPI
       CALL MPI_Abort ( ICON_comm, 1, ierr )
#endif
    ENDIF
    ! -------------------------------------------------------------------
    ! Store grid parameters and global index list
    ! -------------------------------------------------------------------


    ALLOCATE ( gptr%glob_index_rank(grid_shape(1):grid_shape(2)), STAT = ierr )
    IF ( ierr > 0 ) THEN
       WRITE ( * , * ) ' Error allocating rank array '
       ierror = 1
#ifndef NOMPI
       CALL MPI_Abort ( ICON_comm, 1, ierr )
#endif
    ENDIF

    gptr%glob_index_rank(grid_shape(1):grid_shape(2)) = &
         & glob_index_rank(grid_shape(1):grid_shape(2))

  END SUBROUTINE ICON_cpl_def_location

END MODULE mo_icon_cpl_def_grid
