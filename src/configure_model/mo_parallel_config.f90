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
  USE mo_impl_constants,     ONLY: max_dom

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
       &  num_restart_procs,                                        &
       &  use_icon_comm, icon_comm_debug, max_send_recv_buffer_size,&
       &  use_dycore_barrier, itype_exch_barrier, use_dp_mpi2io,    &
       &  icon_comm_method, icon_comm_openmp, max_no_of_comm_variables, &
       &  max_no_of_comm_processes, max_no_of_comm_patterns,        &
       &  sync_barrier_mode, max_mpi_message_size, use_physics_barrier, &
       &  redrad_split_factor
  PUBLIC :: ext_div_medial, ext_div_medial_cluster, ext_div_medial_redrad, &
       & ext_div_medial_redrad_cluster, ext_div_from_file

  PUBLIC :: set_nproma, get_nproma, check_parallel_configuration, use_async_restart_output

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

  ! Type of parallel I/O
  INTEGER :: pio_type = 1

  INTEGER :: num_io_procs = 0

  ! The number of PEs used for writing restart files (0 means, the worker PE0 writes)
  INTEGER :: num_restart_procs = 0

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
    IF (num_restart_procs /= 0) THEN
      CALL message(method_name, &
       & 'num_restart_procs has no effect if the model is compiled with the NOMPI compiler directive')
      CALL message(method_name, &
       & '--> num_restart_procs set to 0')
      num_restart_procs = 0
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
    IF(num_io_procs < 0) num_io_procs = 0
    IF(num_restart_procs < 0) num_restart_procs = 0

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
