!>
!!
!! @par Revision History
!!!
!! @par Copyright
!! 2013 by DKRZ
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
!! $Id: n/a$
!!
MODULE mo_dist_dir

  USE mo_exception, ONLY: finish
  USE mo_mpi, ONLY: p_alltoall, p_alltoallv

  IMPLICIT NONE

  TYPE t_dist_dir
    INTEGER, ALLOCATABLE :: owner(:)
    INTEGER :: global_size, local_start_index
    INTEGER :: comm, comm_size, comm_rank
  END TYPE t_dist_dir

  PUBLIC :: dist_dir_setup
  PUBLIC :: dist_dir_get_owners
  PUBLIC :: t_dist_dir

CONTAINS

  !> Function for setting up a distributed directory
  !> @param[in] owned_indices array containing all indices owned by the
  !>                          local process (indices have to be in the range
  !>                          1 : global_size)
  !> @param[in] global_size   highest global index
  !> @param[in] comm          communicator that contains all processes that
  !>                          are part of the distributed directory
  !> @param[in] comm_rank     rank of process in the provided communicator
  !> @param[in] comm_size     number of processes in the provided
  !>                          communicator
  !> @returns distributed directory
  !> @remarks this function is collective for all processes in comm
  SUBROUTINE dist_dir_setup(dist_dir, owned_indices, global_size, comm, &
    &                       comm_rank, comm_size)

    TYPE(t_dist_dir), INTENT(OUT) :: dist_dir
    INTEGER, INTENT(IN) :: owned_indices(:)
    INTEGER, INTENT(IN) :: global_size
    INTEGER, INTENT(IN) :: comm
    INTEGER, INTENT(IN) :: comm_rank
    INTEGER, INTENT(IN) :: comm_size

    INTEGER :: num_indices_per_process
    INTEGER, POINTER :: recv_buffer(:)
    INTEGER :: dummy(comm_size)
    INTEGER :: num_recv_indices_per_process(comm_size)
    INTEGER :: i, j, n

    num_indices_per_process = (global_size + comm_size - 1) / comm_size
    dist_dir%global_size = global_size
    dist_dir%local_start_index = comm_rank * num_indices_per_process + 1
    dist_dir%comm = comm
    dist_dir%comm_rank = comm_rank
    dist_dir%comm_size = comm_size

    ALLOCATE(dist_dir%owner(num_indices_per_process))

    CALL distribute_indices(owned_indices, recv_buffer, dummy, &
      &                     num_recv_indices_per_process, &
      &                     num_indices_per_process, comm, comm_size)

    dist_dir%owner(:) = -1
    n = 0
    DO i = 1, comm_size
      DO j = 1, num_recv_indices_per_process(i)
        n = n + 1
        dist_dir%owner(recv_buffer(n)) = i - 1
      END DO
    END DO

    DEALLOCATE(recv_buffer)

  END SUBROUTINE dist_dir_setup

  !> gets for each provided global index the rank of the process it
  !> belongs to
  !> @param[in] dist_dir distributed directory
  !> @param[in] indices  indices for which the owner are to be retrieved
  !> @return returns the owners of the provided indices
  !> @remark this routine is collective for all processes that belong to
  !>         dist_dir
  !> @remark all provided global indices need to be valid
  FUNCTION dist_dir_get_owners(dist_dir, indices)

    TYPE(t_dist_dir), INTENT(IN) :: dist_dir
    INTEGER, INTENT(IN) :: indices(:)
    INTEGER :: dist_dir_get_owners(SIZE(indices))

    INTEGER :: num_indices_per_process
    INTEGER :: num_send_indices_per_process(dist_dir%comm_size)
    INTEGER :: num_recv_indices_per_process(dist_dir%comm_size)
    INTEGER, POINTER :: recv_buffer(:)
    INTEGER, ALLOCATABLE :: send_buffer(:)
    INTEGER :: send_displ(dist_dir%comm_size), recv_displ(dist_dir%comm_size+1)

    INTEGER :: i, n

    IF (ANY(indices(:) < 1 .OR. indices(:) > dist_dir%global_size)) &
      CALL finish("get_owner_ranks", "invalid global index")

    recv_buffer => null()

    num_indices_per_process = (dist_dir%global_size + dist_dir%comm_size - 1) &
      &                       / dist_dir%comm_size

    CALL distribute_indices(indices, recv_buffer, &
      &                     num_recv_indices_per_process, &
      &                     num_send_indices_per_process, &
      &                     num_indices_per_process, dist_dir%comm, &
      &                     dist_dir%comm_size)

    ALLOCATE(send_buffer(SIZE(recv_buffer(:))))

    send_buffer(:) = dist_dir%owner(recv_buffer(:))

    DEALLOCATE(recv_buffer)
    ALLOCATE(recv_buffer(SIZE(indices(:))))

    send_displ(1) = 0
    recv_displ(1) = 0
    DO i = 2, dist_dir%comm_size
      send_displ(i) = send_displ(i-1) + num_send_indices_per_process(i-1)
      recv_displ(i) = recv_displ(i-1) + num_recv_indices_per_process(i-1)
    END DO
    recv_displ(dist_dir%comm_size+1) = recv_displ(dist_dir%comm_size) + &
      num_recv_indices_per_process(dist_dir%comm_size)

    CALL p_alltoallv(send_buffer, num_send_indices_per_process, send_displ, &
      &              recv_buffer, num_recv_indices_per_process, recv_displ, &
      &              dist_dir%comm)

    DO i = SIZE(indices(:)), 1, -1
      n = 2 + (indices(i)-1) / num_indices_per_process
      dist_dir_get_owners(i) = recv_buffer(recv_displ(n))
      recv_displ(n) = recv_displ(n) - 1
    END DO

    DEALLOCATE(send_buffer, recv_buffer)

  END FUNCTION dist_dir_get_owners

  SUBROUTINE distribute_indices(indices, recv_buffer, &
    &                           num_send_indices_per_process, &
    &                           num_recv_indices_per_process, &
    &                           num_indices_per_process, comm, comm_size)

    INTEGER, INTENT(IN) :: indices(:)
    INTEGER, INTENT(OUT), POINTER :: recv_buffer(:)
    INTEGER, INTENT(IN) :: comm, comm_size
    INTEGER, INTENT(OUT) :: num_send_indices_per_process(comm_size)
    INTEGER, INTENT(OUT) :: num_recv_indices_per_process(comm_size)
    INTEGER, INTENT(IN) :: num_indices_per_process

    INTEGER :: send_displ(comm_size+1), recv_displ(comm_size)
    INTEGER, ALLOCATABLE :: send_buffer(:)
    INTEGER :: i, j

    num_send_indices_per_process(:) = 0
    DO i = 1, SIZE(indices(:))
      j = 1 + (indices(i)-1) / num_indices_per_process
      num_send_indices_per_process(j) = &
        num_send_indices_per_process(j) + 1
    END DO

    call p_alltoall(num_send_indices_per_process(:), &
      &             num_recv_indices_per_process(:), comm)

    ALLOCATE(send_buffer(SUM(num_send_indices_per_process(:))), &
      &      recv_buffer(SUM(num_recv_indices_per_process(:))))

    send_displ(1:2) = 0
    recv_displ(1) = 0
    DO i = 2, comm_size
      send_displ(i+1) = send_displ(i) + num_send_indices_per_process(i-1)
      recv_displ(i)   = recv_displ(i-1) + num_recv_indices_per_process(i-1)
    END DO

    DO i = 1, SIZE(indices(:))
      j = 2 + (indices(i)-1) / num_indices_per_process
      send_displ(j) = send_displ(j) + 1
      send_buffer(send_displ(j)) = MOD(indices(i)-1, num_indices_per_process) + 1
    END DO

    CALL p_alltoallv(send_buffer, num_send_indices_per_process, send_displ, &
      &              recv_buffer, num_recv_indices_per_process, recv_displ, &
      &              comm)

    DEALLOCATE(send_buffer)

  END SUBROUTINE distribute_indices

END MODULE mo_dist_dir
