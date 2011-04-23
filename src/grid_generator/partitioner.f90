PROGRAM partitioner

  USE mo_kind
  USE mo_io_units,  ONLY: filename_max
  USE mo_exception, ONLY: message, message_text
  USE mo_grid
  USE mo_io_grid
  USE mo_math_constants

  IMPLICIT NONE

  TYPE(t_grid) , ALLOCATABLE, TARGET :: gg(:)

  CHARACTER(LEN=filename_max) :: filename

  INTEGER, ALLOCATABLE :: xadj(:), adjcny(:)
  INTEGER, ALLOCATABLE :: vwgt(:), adjwgt(:)
#ifdef _METIS
  INTEGER :: edgecut
#endif
  INTEGER :: options(5)

  INTEGER, ALLOCATABLE :: part(:)

  INTEGER :: i, j, i_nc, i_ne

  WRITE(message_text,'(a)')  'Running partitioner'
  CALL message ('', TRIM(message_text))

  ALLOCATE (gg(0:4))

  DO i = 0, 4

    WRITE(filename,'(a,i2.2,a)')'GRIDMAP.', i, '.nc'

    CALL input_grid(gg(i), filename)

  ENDDO


  !----------------------------------------------------------------------------
  ! METIS language translation:
  !
  ! vertex: cells%center
  ! edge:   cells%neighbor_index
  ! vwgt:   weights associated with a vertex
  ! adjwgt: weights associated with an edge

  i_nc = SIZE(gg(4)%cells%center)
  i_ne = SIZE(gg(4)%edges%center)

  ALLOCATE(xadj(i_nc+1), adjcny(2*i_ne))
  ALLOCATE(vwgt(i_nc), adjwgt(2*i_ne))

  xadj(1) = 1
  DO i = 1, i_nc
    xadj(i+1) = xadj(i)+3
    adjcny(xadj(i):xadj(i+1)-1) = gg(4)%cells%neighbor_index(i,1:3)
  ENDDO

  ! weights are initially 1

  vwgt(:) = 1
  adjwgt(:) = 1

  options(1) = 0

  DO i = 1, i_nc
    ! Europe: lon = [-20,40], lat = [30,80]
    IF (rad2deg*gg(4)%cells%center(i)%lon > -20.0_dp &
      & .and. rad2deg*gg(4)%cells%center(i)%lon < 40.0_dp) THEN
      IF (rad2deg*gg(4)%cells%center(i)%lat > 30.0_dp &
        & .and. rad2deg*gg(4)%cells%center(i)%lat < 80.0_dp) THEN
        vwgt(i) = 16
      ENDIF
    ENDIF
  ENDDO

  ALLOCATE(part(i_nc))
  part(:) = 0

#ifdef _METIS
  CALL metis_partgraphkway(i_nc, xadj, adjcny, vwgt, adjwgt, &
    & 3, 1, 32, options, edgecut, part)
#endif

  DO i = 1, i_nc
    WRITE (12,'(a,i0)') '> -Z', part(i)
    DO j = 1, 3
      WRITE (12,'(2f15.6)') &
        & gg(4)%verts%vertex(gg(4)%cells%vertex_index(i,j))%lon*rad2deg, &
        & gg(4)%verts%vertex(gg(4)%cells%vertex_index(i,j))%lat*rad2deg
    ENDDO
  ENDDO

END PROGRAM partitioner
