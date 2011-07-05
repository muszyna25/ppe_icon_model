!>
!!        
!! @par Revision History
!!   Revision History in mo_global_variables.f90 (r3513)
!!   Created by Leonidas Linardakis, MPI-M, 2011-05-07
!!
!! @par Copyright
!! 2002-2006 by DWD and MPI-M
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
MODULE mo_parallel_configuration

!    USE mo_kind,               ONLY: wp
!   USE mo_exception,          ONLY: message, message_text, finish

  IMPLICIT NONE

  PRIVATE
  ! Exported variables:
  PUBLIC :: nproma
#ifdef __OMP_RADIATION__
  PUBLIC :: radiation_threads, nh_stepping_threads
#endif
  PUBLIC :: n_ghost_rows,                                     &
       &    div_from_file, div_geometric, div_metis, division_method, &
       &    l_log_checks, l_fast_sum,                                 &
       &    p_test_run, l_test_openmp,                                &
       &    num_test_procs, num_work_procs, num_io_procs,             &
       &    p_test_pe, p_work_pe0, p_io_pe0,                          &
       &    p_n_work, p_pe_work, p_comm_work, p_comm_work_test,       &
       &    p_comm_work_io, p_comm_work_2_io, p_comm_input_bcast,     &
       &    pio_type, itype_comm, iorder_sendrecv
  PUBLIC :: set_nproma, get_nproma
  
  ! computing setup
  ! ---------------
  INTEGER            :: nproma              ! inner loop length/vector length

  ! Number of rows of ghost cells
  INTEGER :: n_ghost_rows


  ! Division method for area subdivision
  INTEGER, PARAMETER :: div_from_file = 0  ! Read from file
  INTEGER, PARAMETER :: div_geometric = 1  ! Geometric subdivision
  INTEGER, PARAMETER :: div_metis     = 2  ! Use Metis

  INTEGER :: division_method 

  ! Flag if checks in a verification run should be logged

  LOGICAL :: l_log_checks

  ! Flag if fast but nonreproducible sum should be used

  LOGICAL :: l_fast_sum

  ! Please note for the following variables: The default settings are for NO_MPI runs!

  ! p_test_run indicates a verification run, i.e. a run where 1 PE runs the complete
  ! model whereas the other PEs do a real parallelized run

  LOGICAL :: p_test_run

  ! if l_test_openmp is set together with p_test_run, then the verification PE uses
  ! only 1 thread. This allows for verifying the OpenMP implementation

  LOGICAL :: l_test_openmp

  ! Processor distribution:
  ! num_test_procs: 0 or 1
  ! num_work_procs: number of procs running in parallel on the model
  ! num_io_procs:   number of procs for I/O
  ! num_test_procs + num_work_procs + num_io_procs = p_nprocs

  INTEGER :: num_test_procs
  INTEGER :: num_work_procs
  INTEGER :: num_io_procs

  ! Note: p_test_pe, p_work_pe0, p_io_pe0 are identical on all PEs

  ! In a verification run, p_test_pe is the number of the PE running the complete model,
  ! otherwise it contains -1

  INTEGER :: p_test_pe     ! Number of test PE
  INTEGER :: p_work_pe0    ! Number of workgroup PE 0 within all PEs
  INTEGER :: p_io_pe0      ! Number of I/O PE 0 within all PEs (p_nprocs if no I/O PEs)

  ! Note: p_n_work, p_pe_work are NOT identical on all PEs

  ! p_n_work: Number of PEs working together:
  ! - num_work_procs for non-verification runs
  ! - num_work_procs for verification runs on pes != p_test_pe
  ! - 1              for verification runs on p_test_pe
  ! - num_io_procs   always on I/O pes

  INTEGER :: p_n_work
  INTEGER :: p_pe_work        ! PE number within work group

  ! MPI communicators
  INTEGER :: p_comm_work_io     ! Communicator spanning work group and I/O PEs
  INTEGER :: p_comm_work_2_io   ! Inter(!)communicator work PEs - I/O PEs
  INTEGER :: p_comm_input_bcast ! Communicator for broadcasts in NetCDF input

  ! Type of parallel I/O

  INTEGER :: pio_type

  ! Type of (halo) communication: 
  ! 1 = synchronous communication with local memory for exchange buffers
  ! 2 = synchronous communication with global memory for exchange buffers
  ! 3 = asynchronous communication within dynamical core with global memory 
  !     for exchange buffers (not yet implemented)
  INTEGER :: itype_comm

  ! Order of send/receive sequence in exchange routines
  ! 1 = irecv, send
  ! 2 = isend, recv
  INTEGER :: iorder_sendrecv

  INTEGER :: radiation_threads
  INTEGER :: nh_stepping_threads 

  ! MPI communicators
  INTEGER :: p_comm_work        ! Communicator for work group
  INTEGER :: p_comm_work_test   ! Communicator spanning work group and test PE



CONTAINS
  
  !-------------------------------------------------------------------------
  !>
  SUBROUTINE set_nproma(new_nproma)
    INTEGER, INTENT(in) :: new_nproma

    nproma = new_nproma

  END SUBROUTINE set_nproma
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  !>
  INTEGER FUNCTION get_nproma()

    get_nproma = nproma

  END FUNCTION get_nproma
  !-------------------------------------------------------------------------


END MODULE mo_parallel_configuration
