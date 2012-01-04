!>
!!        
!! @par Revision History
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
MODULE mo_parallel_config

  USE mo_exception,          ONLY: message, finish
  USE mo_io_units,           ONLY: filename_max

  IMPLICIT NONE

  PRIVATE
  ! Exported variables:
  PUBLIC :: nproma, openmp_threads
  PUBLIC :: radiation_ompthreads, nh_stepping_ompthreads, parallel_radiation_omp
  PUBLIC :: parallel_radiation_mpi, test_parallel_radiation

  PUBLIC :: n_ghost_rows,                                     &
       &    div_from_file, div_geometric, div_metis, division_method, &
       &    division_file_name, l_log_checks, l_fast_sum,             &
       &    p_test_run, l_test_openmp,                                &
       &    pio_type, itype_comm, iorder_sendrecv, num_io_procs,      &
       &    use_icon_comm, icon_comm_debug, max_send_recv_buffer_size, &
       &    use_dycore_barrier
       
  PUBLIC :: set_nproma, get_nproma, check_parallel_configuration
  
  ! computing setup
  ! ---------------
  INTEGER  :: nproma = 1              ! inner loop length/vector length

  ! Number of rows of ghost cells
  INTEGER :: n_ghost_rows = 1
  
  INTEGER :: openmp_threads = -1 ! < 0 means this value is not used,
                                 ! instead the system's value will be used
  
  ! Division method for area subdivision
  INTEGER, PARAMETER :: div_from_file = 0  ! Read from file
  INTEGER, PARAMETER :: div_geometric = 1  ! Geometric subdivision
  INTEGER, PARAMETER :: div_metis     = 2  ! Use Metis

  INTEGER :: division_method = 1
  CHARACTER(LEN=filename_max) :: division_file_name ! if div_from_file

  ! Flag if checks in a verification run should be logged
  LOGICAL :: l_log_checks = .false.

  ! Flag if fast but nonreproducible sum should be used
  LOGICAL :: l_fast_sum = .false.

  ! Please note for the following variables: The default settings are for NO_MPI runs!

  ! p_test_run indicates a verification run, i.e. a run where 1 PE runs the complete
  ! model whereas the other PEs do a real parallelized run
  LOGICAL :: p_test_run = .false.

  LOGICAL :: use_dycore_barrier = .false. ! put an mpi barrier before the dycore to eliminate imbalances
  
  ! if l_test_openmp is set together with p_test_run, then the verification PE uses
  ! only 1 thread. This allows for verifying the OpenMP implementation
  LOGICAL :: l_test_openmp = .false.

  LOGICAL :: parallel_radiation_omp = .false.
  LOGICAL :: parallel_radiation_mpi = .false.
  LOGICAL :: test_parallel_radiation = .false.

  LOGICAL :: use_icon_comm = .false.
  LOGICAL :: icon_comm_debug= .false.
  INTEGER :: max_send_recv_buffer_size = 131072
  ! Type of parallel I/O
  INTEGER :: pio_type = 1
  
  INTEGER :: num_io_procs = 0

  ! Type of (halo) communication: 
  ! 1 = synchronous communication with local memory for exchange buffers
  ! 2 = synchronous communication with global memory for exchange buffers
  ! 3 = asynchronous communication within dynamical core with global memory 
  !     for exchange buffers (not yet implemented)
  INTEGER :: itype_comm = 1

  ! Order of send/receive sequence in exchange routines
  ! 1 = irecv, send
  ! 2 = isend, recv
  INTEGER :: iorder_sendrecv = 1

  INTEGER :: radiation_ompthreads   = 1
  INTEGER :: nh_stepping_ompthreads = 1



CONTAINS
  
  !-------------------------------------------------------------------------
  SUBROUTINE check_parallel_configuration()
 
!   !local variables
!     INTEGER :: i_status, my_color, peer_comm, p_error
    CHARACTER(*), PARAMETER :: method_name = "check_parallel_configuration"

    !------------------------------------------------------------
    !  check the consistency of the parameters
    !------------------------------------------------------------
    IF (nproma<=0) CALL finish(TRIM(method_name),'"nproma" must be positive')

! check l_test_openmp
#ifndef _OPENMP
  IF (l_test_openmp) THEN
    CALL message(method_name, &
       & 'l_test_openmp has no effect if the model is compiled without OpenMP support')
    CALL message(method_name, &
       & '--> l_test_openmp set to .FALSE.')
    l_test_openmp = .FALSE.
  END IF
#endif

    ! check p_test_run and num_io_procs
#ifdef NOMPI
    ! Unconditionally set p_test_run to .FALSE. and num_io_procs to 0,
    ! all other variables are already set correctly
    IF (p_test_run) THEN
      CALL message(method_name, &
       & 'p_test_run has no effect if the model is compiled with the NOMPI compiler directive')
      CALL message(method_name, &
       & '--> p_test_run set to .FALSE.')
      p_test_run = .FALSE.
    END IF
    IF (num_io_procs /= 0) THEN
      CALL message(method_name, &
       & 'num_io_procs has no effect if the model is compiled with the NOMPI compiler directive')
      CALL message(method_name, &
       & '--> num_io_procs set to 0')
      num_io_procs = 0
    END IF

#else

    ! check n_ghost_rows
    IF (n_ghost_rows<1) THEN
      CALL finish(method_name, &
          & 'n_ghost_rows<1 in parallel_nml namelist is not allowed')
    END IF

    ! check division_method
    SELECT CASE (division_method)
    CASE(div_from_file, div_geometric)
      ! ok
    CASE(div_metis)
#ifdef HAVE_METIS
    ! ok
#else
      CALL finish(method_name, &
        & 'division_method=div_metis=2 in parallel_nml namelist is not allowed')
#endif
    CASE DEFAULT
      CALL finish(method_name, &
        & 'value of division_method in parallel_nml namelist is not allowed')
    END SELECT
    ! for safety only
    IF(num_io_procs < 0) num_io_procs = 0
  
#endif

  END SUBROUTINE check_parallel_configuration
  !-------------------------------------------------------------------------
  
  
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


END MODULE mo_parallel_config
