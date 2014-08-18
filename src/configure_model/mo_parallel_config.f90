!>
!!
!! @par Revision History
!!   Created by Leonidas Linardakis, MPI-M, 2011-05-07
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
MODULE mo_parallel_config

  USE mo_exception,          ONLY: message, finish
  USE mo_io_units,           ONLY: filename_max
  USE mo_impl_constants,     ONLY: max_dom, MAX_NUM_IO_PROCS
  USE mo_util_string,        ONLY: int2string

  IMPLICIT NONE

  PRIVATE
  ! Exported variables:
  PUBLIC :: nproma, openmp_threads
!   PUBLIC :: radiation_ompthreads, nh_stepping_ompthreads, parallel_radiation_omp
  PUBLIC :: parallel_radiation_mode, test_parallel_radiation

  PUBLIC :: n_ghost_rows,                                     &
       &  div_from_file, div_geometric, div_metis, division_method, &
       &  division_file_name, radiation_division_file_name, &
       &  l_log_checks, l_fast_sum, ldiv_phys_dom,                  &
       &  p_test_run, l_test_openmp,                                &
       &  pio_type, itype_comm, iorder_sendrecv, num_io_procs,      &
       &  num_restart_procs, num_pref_procs,                                       &
       &  use_icon_comm, icon_comm_debug, max_send_recv_buffer_size,&
       &  use_dycore_barrier, itype_exch_barrier, use_dp_mpi2io,    &
       &  icon_comm_method, icon_comm_openmp, max_no_of_comm_variables, &
       &  max_no_of_comm_processes, max_no_of_comm_patterns,        &
       &  sync_barrier_mode, max_mpi_message_size, use_physics_barrier, &
       &  redrad_split_factor, restart_chunk_size
  PUBLIC :: ext_div_medial, ext_div_medial_cluster, ext_div_medial_redrad, &
       & ext_div_medial_redrad_cluster, ext_div_from_file

  PUBLIC :: set_nproma, get_nproma, check_parallel_configuration, use_async_restart_output, &
       &    use_async_prefetch

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
  INTEGER, PARAMETER :: ext_div_medial = 101
  INTEGER, PARAMETER :: ext_div_medial_cluster = 102
  INTEGER, PARAMETER :: ext_div_medial_redrad = 103
  INTEGER, PARAMETER :: ext_div_medial_redrad_cluster = 104
  INTEGER, PARAMETER :: ext_div_from_file = 201

  INTEGER :: division_method(0:max_dom) = 1
  CHARACTER(LEN=filename_max) :: division_file_name(0:max_dom)! if div_from_file
  CHARACTER(LEN=filename_max) :: radiation_division_file_name(max_dom)! if parallel_radiation_mode = 1
  INTEGER :: redrad_split_factor = 6

  ! Flag if (in case of merged domains) physical domains shall be considered for
  ! computing the domain decomposition
  LOGICAL :: ldiv_phys_dom = .FALSE.

  ! Flag if checks in a verification run should be logged
  LOGICAL :: l_log_checks = .false.

  ! Flag if fast but nonreproducible sum should be used
  LOGICAL :: l_fast_sum = .false.

  ! Please note for the following variables: The default settings are for NO_MPI runs!

  ! p_test_run indicates a verification run, i.e. a run where 1 PE runs the complete
  ! model whereas the other PEs do a real parallelized run
  LOGICAL :: p_test_run = .false.

  LOGICAL :: use_dycore_barrier = .false. ! acivate an mpi barrier before the dycore
                                          ! to synchronize MPI tasks
  LOGICAL :: use_physics_barrier = .false. ! activate mpi barrier after the physics

  INTEGER :: itype_exch_barrier = 0  ! 1: put an mpi barrier at the beginning of exchange calls to synchronize MPI tasks
                                     ! 2: put an mpi barrier after MPI_WAIT to synchronize MPI tasks
                                     ! 3: 1+2

  ! if l_test_openmp is set together with p_test_run, then the verification PE uses
  ! only 1 thread. This allows for verifying the OpenMP implementation
  LOGICAL :: l_test_openmp = .false.

!   LOGICAL :: parallel_radiation_omp = .false.
  INTEGER :: parallel_radiation_mode(max_dom) = 0
  LOGICAL :: test_parallel_radiation = .false.

  LOGICAL :: use_icon_comm = .false.
  LOGICAL :: icon_comm_debug= .false.
  INTEGER :: max_send_recv_buffer_size = 262144  ! size in doubles (x8)
  INTEGER :: max_mpi_message_size      = 65536   ! size in doubles (x8)
  INTEGER :: max_no_of_comm_variables  = 64
  INTEGER :: max_no_of_comm_processes  = 64
  INTEGER :: max_no_of_comm_patterns   = 32
  INTEGER :: icon_comm_method = 1
  INTEGER :: sync_barrier_mode = 0

  LOGICAL :: icon_comm_openmp = .false.

  ! Flag whether async restart output is used, it is set in the main program:
  LOGICAL :: use_async_restart_output = .FALSE.

  ! Flag whether async prefetching is used, it is set in main program:
  INTEGER :: use_async_prefetch = 0

  ! Type of parallel I/O
  INTEGER :: pio_type = 1

  INTEGER :: num_io_procs = 0

  ! The number of PEs used for writing restart files (0 means, the worker PE0 writes)
  INTEGER :: num_restart_procs = 0

  ! The number of PEs used for async prefetching of input (0 means, the PE0 prefetches input)
  INTEGER :: num_pref_procs = 0

  ! Type of (halo) communication:
  ! 1 = synchronous communication with local memory for exchange buffers
  ! 2 = synchronous communication with global memory for exchange buffers
  ! 3 = asynchronous communication within dynamical core with global memory
  !     for exchange buffers (not yet implemented)
  INTEGER :: itype_comm = 1

  ! Order of send/receive sequence in exchange routines
  ! 1 = irecv, send
  ! 2 = isend, recv
  ! 3 = irecv, isend
  INTEGER :: iorder_sendrecv = 1

!   INTEGER :: radiation_ompthreads   = 1
!   INTEGER :: nh_stepping_ompthreads = 1

  ! Flag. Enable this flag if output fields shall be gathered in
  ! DOUBLE PRECISION. The resulting files are identical to the "normal"
  ! operation with simple REAL, since NetCDFs are written in single precision
  ! anyway by default and GRIB2 has 16 or 24 bit per data word.
  !
  LOGICAL :: use_dp_mpi2io

  ! The (asynchronous) restart is capable of writing and communicating
  ! more than one 2D slice at once
  INTEGER :: restart_chunk_size

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
#ifndef __SX__
    ! migration helper: catch nproma's that were obviously intended
    !                   for a vector machine.
    IF (nproma>256) CALL finish(TRIM(method_name),'The value of "nproma" seems to be set for a vector machine!')
#endif

    icon_comm_openmp = .false.
! check l_test_openmp
#ifndef _OPENMP
    IF (l_test_openmp) THEN
      CALL message(method_name, &
         & 'l_test_openmp has no effect if the model is compiled without OpenMP support')
      CALL message(method_name, &
         & '--> l_test_openmp set to .FALSE.')
      l_test_openmp = .FALSE.
    END IF
#else
    IF (icon_comm_method > 100) &
      & icon_comm_openmp = .true.
#endif

    ! check p_test_run, num_io_procs, num_restart_procs and num_pref_procs
#ifdef NOMPI
    ! Unconditionally set p_test_run to .FALSE. and num_io_procs to 0, num_restart_procs to 0 
    ! and num_pref_procs to 0 
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
    IF (num_restart_procs /= 0) THEN
      CALL message(method_name, &
       & 'num_restart_procs has no effect if the model is compiled with the NOMPI compiler directive')
      CALL message(method_name, &
       & '--> num_restart_procs set to 0')
      num_restart_procs = 0
    END IF
    IF (num_pref_procs /= 0) THEN
      CALL message(method_name, &
       & 'num_pref_procs has no effect if the model is compiled with the NOMPI compiler directive')
      CALL message(method_name, &
       & '--> num_pref_procs set to 0')
      num_pref_procs = 0
    END IF

#else

    ! check n_ghost_rows
    IF (n_ghost_rows<1) THEN
      CALL finish(method_name, &
          & 'n_ghost_rows<1 in parallel_nml namelist is not allowed')
    END IF

    ! check division_method
    ! this will be checked during the decomposition
!     SELECT CASE (division_method)
!     CASE(div_from_file, div_geometric, ext_div_medial, ext_div_medial_cluster, &
!       & ext_div_medial_redrad, ext_div_medial_redrad_cluster)
!       ! ok
!     CASE(div_metis)
! #ifdef HAVE_METIS
!     ! ok
! #else
!       CALL finish(method_name, &
!         & 'division_method=div_metis=2 in parallel_nml namelist is not allowed')
! #endif
!     CASE DEFAULT
!       CALL finish(method_name, &
!         & 'value of division_method in parallel_nml namelist is not allowed')
!     END SELECT
    ! for safety only
    IF (num_io_procs < 0)      num_io_procs = 0
    IF (num_restart_procs < 0) num_restart_procs = 0
    IF (num_io_procs > MAX_NUM_IO_PROCS) THEN
      CALL finish(method_name, "Namelist parameter num_io_procs chosen too large ( > "//TRIM(int2string(MAX_NUM_IO_PROCS))//")!")
    END IF
    IF(num_pref_procs < 0) num_pref_procs = 0
    IF(num_pref_procs > 1) &
         CALL finish(TRIM(method_name),'The no of prefetch processor can be zero or one, but should not be set more than one!')
    IF((num_pref_procs == 1) .AND. (use_async_prefetch /= 1)) THEN
       CALL finish(TRIM(method_name), &
                'The no of prefetch processor is set to one but the use_async_prefetch flag is not set to one!')
    ELSE IF((num_pref_procs /= 1) .AND. (use_async_prefetch == 1)) THEN
       CALL finish(TRIM(method_name), &
                'The flag use_async_prefetch is set to one but the no of prefetch processor is not set to one!')
    ENDIF
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
