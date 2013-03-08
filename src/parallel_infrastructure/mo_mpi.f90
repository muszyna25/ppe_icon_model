MODULE mo_mpi

  ! Comment: Please use basic WRITE to nerr for messaging in the whole
  !          MPI package to achieve proper output.

  ! actual method (MPI-2)
#ifndef NOMPI
#if !defined (__SUNPRO_F95)
  USE mpi
#endif
#endif

#ifdef _OPENMP
  USE omp_lib, ONLY: omp_get_max_threads, omp_set_num_threads
#endif

  USE mo_kind
  USE mo_io_units, ONLY: nerr
  
  IMPLICIT NONE

  PRIVATE                          ! all declarations are private

#ifndef NOMPI
#if defined (__SUNPRO_F95)
  INCLUDE "mpif.h"
#endif
#endif

  ! start/stop methods
  PUBLIC :: start_mpi
  PUBLIC :: p_stop, p_abort

  ! split the global communicator to _process_mpi_communicator
  PUBLIC :: split_global_mpi_communicator
  !The given communicator will be the all communicator for this component
!   PUBLIC :: set_process_mpi_communicator
  ! Sets the test, work, i/o communicators
  PUBLIC :: set_mpi_work_communicators
  ! Sets the p_comm_input_bcast
  PUBLIC :: set_comm_input_bcast
  ! set other parameters
  PUBLIC :: set_process_mpi_name
  
  PUBLIC :: push_glob_comm, pop_glob_comm, get_glob_proc0
  
  ! Logical functions
  PUBLIC :: run_is_global_mpi_parallel
  PUBLIC :: my_process_is_stdio, my_process_is_mpi_parallel, my_process_is_mpi_all_parallel
  PUBLIC :: my_process_is_mpi_seq, my_process_is_mpi_test, my_process_is_mpi_workroot
  PUBLIC :: my_process_is_mpi_ioroot
  PUBLIC :: my_process_is_mpi_all_seq, my_process_is_io

  ! get parameters
  PUBLIC :: get_my_mpi_all_communicator   ! the communicator for the specific component, ie the process_mpi_all_comm
  PUBLIC :: get_my_mpi_all_comm_size   ! this is the the size of the communicator for the specific component
  PUBLIC :: get_my_mpi_work_communicator   ! the communicator for the workers of this component
  PUBLIC :: get_my_mpi_work_comm_size   ! this is the the size of the workers
  
  PUBLIC :: get_mpi_all_workroot_id, get_my_global_mpi_id, get_my_mpi_all_id
  PUBLIC :: get_mpi_all_ioroot_id
  PUBLIC :: get_my_mpi_work_id
  PUBLIC :: default_comm_type, null_comm_type


  ! some public communicators
  PUBLIC :: process_mpi_all_comm
  PUBLIC :: process_mpi_all_test_id, process_mpi_all_workroot_id, &
    &       process_mpi_all_ioroot_id

  
  
  PUBLIC :: p_comm_work, p_comm_work_test
  PUBLIC :: p_comm_work_2_io, p_comm_input_bcast, p_comm_work_io
  PUBLIC :: p_communicator_a, p_communicator_b, p_communicator_d

  PUBLIC :: process_mpi_io_size
  PUBLIC :: process_mpi_stdio_id
  
  ! Main communication methods
  PUBLIC :: global_mpi_barrier
  PUBLIC :: work_mpi_barrier

  PUBLIC :: p_send, p_recv, p_sendrecv, p_bcast, p_barrier
  PUBLIC :: p_isend, p_irecv, p_wait, p_wait_any, &
    &       p_irecv_packed, p_send_packed,        &
    &       p_pack_int, p_pack_real,              &
    &       p_pack_int_1d, p_pack_real_1d,        &
    &       p_pack_string, p_pack_real_2d,        &
    &       p_unpack_int, p_unpack_real,          &
    &       p_unpack_int_1d, p_unpack_real_1d,    &
    &       p_unpack_string, p_unpack_real_2d
  PUBLIC :: p_gather, p_max, p_min, p_sum, p_global_sum, p_field_sum
  PUBLIC :: p_probe
  PUBLIC :: p_allreduce_minloc
  PUBLIC :: p_gather_field
  PUBLIC :: p_gatherv
  PUBLIC :: p_scatterv
  PUBLIC :: p_allreduce_max
  PUBLIC :: p_commit_type_struct
  PUBLIC :: p_alltoall
  PUBLIC :: p_clear_request
  PUBLIC :: p_mpi_wtime

  !----------- to be removed -----------------------------------------
  PUBLIC :: p_pe, p_io
  PUBLIC :: num_test_procs, num_work_procs,     &
       &    p_work_pe0, p_io_pe0,               &
       &    p_n_work, p_pe_work
  !     p_test_pe,
  !--------------------------------------------------------------------

#ifndef NOMPI
  PUBLIC :: MPI_INTEGER, MPI_STATUS_SIZE, MPI_SUCCESS,                     &
            MPI_INFO_NULL, MPI_ADDRESS_KIND, MPI_COMM_NULL, MPI_COMM_SELF, &
            MPI_UNDEFINED
#endif

  PUBLIC :: MPI_ANY_SOURCE

  ! real data type matching real type of MPI implementation
  PUBLIC :: p_real_dp, p_real_sp
  PUBLIC :: p_int
  PUBLIC :: p_int_i8
  PUBLIC :: p_bool
  PUBLIC :: p_address_kind
  PUBLIC :: p_real_dp_byte, p_real_sp_byte

  PUBLIC ::get_mpi_time,ellapsed_mpi_time

  ! old fashioned method (MPI-1)

!!$#ifndef NOMPI
!!$  INCLUDE 'mpif.h'
!!$#endif

  ! general run time information

#ifndef NOMPI
  INTEGER :: version, subversion   ! MPI version
 ! MPI call inherent variables
  INTEGER :: p_error                     ! MPI error number
  INTEGER :: p_status(MPI_STATUS_SIZE)   ! standard information of MPI_RECV
  
  INTEGER, ALLOCATABLE, SAVE :: p_request(:) ! request values for non blocking calls
  INTEGER :: p_irequest ! the first p_irequest values of p_request are in use
  INTEGER :: p_mrequest ! actual size of p_request
  INTEGER, PARAMETER :: p_request_alloc_size = 4096
  INTEGER, PARAMETER :: p_address_kind = MPI_ADDRESS_KIND
#else
  INTEGER, PARAMETER :: p_address_kind = i8    ! should not get touched at all
  INTEGER, PARAMETER :: MPI_COMM_NULL  = 0    
  ! dummy arguments for function calls:
  INTEGER, PARAMETER :: MPI_ANY_SOURCE = 0
#endif


  ! public parallel run information
  CHARACTER(len=64) :: global_mpi_name
  CHARACTER(len=64) :: process_mpi_name
  ! Do not change the following parameters
  INTEGER, PARAMETER :: process_mpi_stdio_id = 0
  INTEGER, PARAMETER :: process_mpi_root_id = 0
  INTEGER, PARAMETER :: default_comm_type = 1
  INTEGER, PARAMETER :: null_comm_type = 0
  
  ! communicator sets
  ! this is the global communicator
  INTEGER :: global_mpi_communicator  ! replaces MPI_COMM_WORLD 
  INTEGER :: global_mpi_size          ! total number of processes in global world
  INTEGER :: my_global_mpi_id         ! process id in global world
  LOGICAL :: is_global_mpi_parallel
  
  ! this is the communicator for the whole component (atmo/ocean/etc)
  INTEGER :: process_mpi_all_comm     ! communicator in the whole model-component
  INTEGER :: process_mpi_all_size     ! total number of processes in the whole model-component
  INTEGER :: my_process_mpi_all_id
  INTEGER :: process_mpi_all_workroot_id  ! the root process in component
  INTEGER :: process_mpi_all_ioroot_id    ! the first I/O process
  INTEGER :: process_mpi_all_test_id  ! the test process in component

  
  LOGICAL :: process_is_mpi_parallel
  LOGICAL :: process_is_stdio
  LOGICAL :: is_mpi_test_run = .false.
  LOGICAL :: is_openmp_test_run = .false.
  
  
  ! this is the local work communicator (computation, i/o, etc)
!   INTEGER :: process_mpi_local_comm     ! communicator in the work group
!   INTEGER :: process_mpi_local_size     ! total number of processes in the whole model-component
!   INTEGER :: my_process_mpi_local_id
  
  INTEGER :: my_mpi_function  ! test, or work, or i/o
  INTEGER, PARAMETER :: test_mpi_process = 1
  INTEGER, PARAMETER :: work_mpi_process = 2
  INTEGER, PARAMETER :: io_mpi_process = 3
  
  !------------------------------------------------------------
  ! Processor distribution:
  ! num_test_procs: 0 or 1
  ! num_work_procs: number of procs running in parallel on the model
  ! process_mpi_io_size:   number of procs for I/O
  ! num_test_procs + num_work_procs + process_mpi_io_size = process_mpi_all_size
  INTEGER :: num_test_procs
  INTEGER :: num_work_procs
  INTEGER :: process_mpi_io_size

  ! Note: p_test_pe, p_work_pe0, p_io_pe0 are identical on all PEs

  ! In a verification run, p_test_pe is the number of the PE running the complete model,
  ! otherwise it contains -1

  INTEGER :: p_test_pe     ! Number of test PE
  INTEGER :: p_work_pe0    ! Number of workgroup PE 0 within all PEs
  INTEGER :: p_io_pe0      ! Number of I/O PE 0 within all PEs (process_mpi_all_size if no I/O PEs)

  ! Note: p_n_work, p_pe_work are NOT identical on all PEs

  ! p_n_work: Number of PEs working together:
  ! - num_work_procs for non-verification runs
  ! - num_work_procs for verification runs on pes != p_test_pe
  ! - 1              for verification runs on p_test_pe
  ! - num_io_procs   always on I/O pes

  INTEGER :: p_n_work
  INTEGER :: p_pe_work        ! PE number within work group

  ! MPI communicators
  INTEGER :: p_comm_work        ! Communicator for work group
  INTEGER :: p_comm_work_test   ! Communicator spanning work group and test PE
  INTEGER :: p_comm_work_io     ! Communicator spanning work group and I/O PEs
  INTEGER :: p_comm_work_2_io   ! Inter(!)communicator work PEs - I/O PEs
  INTEGER :: p_comm_input_bcast ! Communicator for broadcasts in NetCDF input

  INTEGER :: p_communicator_a ! for Set A
  INTEGER :: p_communicator_b ! for Set B
  INTEGER :: p_communicator_d ! for debug node

  INTEGER :: p_pe     = 0     ! this is the PE number of this task
  INTEGER :: p_io     = 0     ! PE number of PE handling IO
  
! non blocking calls

  ! module intrinsic names

!!#ifndef NOMPI
!!LK  INTEGER :: iope                  ! PE able to do IO
!!#endif

#ifdef DEBUG
  INTEGER :: nbcast                ! counter for broadcasts for debugging
#endif

  ! MPI native transfer types

  INTEGER :: p_int     = 0
  INTEGER :: p_real    = 0
  INTEGER :: p_bool    = 0
  INTEGER :: p_char    = 0

  ! MPI transfer types

  INTEGER :: p_real_dp = 0
  INTEGER :: p_real_sp = 0
  INTEGER :: p_int_i4  = 0
  INTEGER :: p_int_i8  = 0

  ! MPI native transfer types byte size

  INTEGER :: p_int_byte     = 0
  INTEGER :: p_real_byte    = 0
!  INTEGER :: p_bool_byte    = 0
!  INTEGER :: p_char_byte    = 0

  ! MPI transfer types byte size

  INTEGER :: p_real_dp_byte = 0
  INTEGER :: p_real_sp_byte = 0
  INTEGER :: p_int_i4_byte  = 0
  INTEGER :: p_int_i8_byte  = 0

  CHARACTER(len=256) :: message_text = ''

  ! Flag if processor splitting is active
  LOGICAL, PUBLIC :: proc_split = .FALSE.

  ! communicator stack for global sums
  INTEGER, PARAMETER :: max_lev = 10 ! 2 is sufficient
  INTEGER, PUBLIC :: comm_lev = 0, glob_comm(0:max_lev), comm_proc0(0:max_lev)

  ! define generic interfaces to allow proper compiling with picky compilers
  ! like NAG f95 for clean argument checking and shortening the call sequence.

  INTERFACE p_send
     MODULE PROCEDURE p_send_char
     MODULE PROCEDURE p_send_real
     MODULE PROCEDURE p_send_int
     MODULE PROCEDURE p_send_bool
     MODULE PROCEDURE p_send_real_1d
     MODULE PROCEDURE p_send_int_1d
     MODULE PROCEDURE p_send_bool_1d
     MODULE PROCEDURE p_send_real_2d
     MODULE PROCEDURE p_send_int_2d
     MODULE PROCEDURE p_send_bool_2d
     MODULE PROCEDURE p_send_real_3d
     MODULE PROCEDURE p_send_int_3d
     MODULE PROCEDURE p_send_bool_3d
     MODULE PROCEDURE p_send_real_4d
     MODULE PROCEDURE p_send_int_4d
     MODULE PROCEDURE p_send_bool_4d
     MODULE PROCEDURE p_send_real_5d
  END INTERFACE

  INTERFACE p_isend
     MODULE PROCEDURE p_isend_char
     MODULE PROCEDURE p_isend_real
     MODULE PROCEDURE p_isend_int
     MODULE PROCEDURE p_isend_bool
     MODULE PROCEDURE p_isend_real_1d
     MODULE PROCEDURE p_isend_int_1d
     MODULE PROCEDURE p_isend_bool_1d
     MODULE PROCEDURE p_isend_real_2d
     MODULE PROCEDURE p_isend_int_2d
     MODULE PROCEDURE p_isend_bool_2d
     MODULE PROCEDURE p_isend_real_3d
     MODULE PROCEDURE p_isend_int_3d
     MODULE PROCEDURE p_isend_bool_3d
     MODULE PROCEDURE p_isend_real_4d
     MODULE PROCEDURE p_isend_int_4d
     MODULE PROCEDURE p_isend_bool_4d
     MODULE PROCEDURE p_isend_real_5d
  END INTERFACE

  INTERFACE p_recv
     MODULE PROCEDURE p_recv_char
     MODULE PROCEDURE p_recv_real
     MODULE PROCEDURE p_recv_int
     MODULE PROCEDURE p_recv_bool
     MODULE PROCEDURE p_recv_real_1d
     MODULE PROCEDURE p_recv_int_1d
     MODULE PROCEDURE p_recv_bool_1d
     MODULE PROCEDURE p_recv_real_2d
     MODULE PROCEDURE p_recv_int_2d
     MODULE PROCEDURE p_recv_bool_2d
     MODULE PROCEDURE p_recv_real_3d
     MODULE PROCEDURE p_recv_int_3d
     MODULE PROCEDURE p_recv_bool_3d
     MODULE PROCEDURE p_recv_real_4d
     MODULE PROCEDURE p_recv_int_4d
     MODULE PROCEDURE p_recv_bool_4d
     MODULE PROCEDURE p_recv_real_5d
  END INTERFACE

  INTERFACE p_irecv
     MODULE PROCEDURE p_irecv_char
     MODULE PROCEDURE p_irecv_real
     MODULE PROCEDURE p_irecv_int
     MODULE PROCEDURE p_irecv_bool
     MODULE PROCEDURE p_irecv_real_1d
     MODULE PROCEDURE p_irecv_int_1d
     MODULE PROCEDURE p_irecv_bool_1d
     MODULE PROCEDURE p_irecv_real_2d
     MODULE PROCEDURE p_irecv_int_2d
     MODULE PROCEDURE p_irecv_bool_2d
     MODULE PROCEDURE p_irecv_real_3d
     MODULE PROCEDURE p_irecv_int_3d
     MODULE PROCEDURE p_irecv_bool_3d
     MODULE PROCEDURE p_irecv_real_4d
     MODULE PROCEDURE p_irecv_int_4d
     MODULE PROCEDURE p_irecv_bool_4d
  END INTERFACE

  INTERFACE p_sendrecv
     MODULE PROCEDURE p_sendrecv_real_1d
     MODULE PROCEDURE p_sendrecv_real_2d
     MODULE PROCEDURE p_sendrecv_real_3d
     MODULE PROCEDURE p_sendrecv_real_4d
  END INTERFACE

  INTERFACE p_bcast
     MODULE PROCEDURE p_bcast_real
     MODULE PROCEDURE p_bcast_int_i4
     MODULE PROCEDURE p_bcast_int_i8
     MODULE PROCEDURE p_bcast_bool
     MODULE PROCEDURE p_bcast_real_1d
     MODULE PROCEDURE p_bcast_int_1d
     MODULE PROCEDURE p_bcast_int_i8_1d
     MODULE PROCEDURE p_bcast_bool_1d
     MODULE PROCEDURE p_bcast_real_2d
     MODULE PROCEDURE p_bcast_real_2d_single
     MODULE PROCEDURE p_bcast_int_2d
     MODULE PROCEDURE p_bcast_bool_2d
     MODULE PROCEDURE p_bcast_real_3d
     MODULE PROCEDURE p_bcast_int_3d
     MODULE PROCEDURE p_bcast_bool_3d
     MODULE PROCEDURE p_bcast_real_4d
     MODULE PROCEDURE p_bcast_int_4d
     MODULE PROCEDURE p_bcast_bool_4d
     MODULE PROCEDURE p_bcast_real_5d
     MODULE PROCEDURE p_bcast_char
     MODULE PROCEDURE p_bcast_cchar
     MODULE PROCEDURE p_bcast_char_1d
  END INTERFACE

  INTERFACE p_gather
     MODULE PROCEDURE p_gather_real_0d1d
     MODULE PROCEDURE p_gather_real_1d2d
     MODULE PROCEDURE p_gather_real_5d6d
     MODULE PROCEDURE p_gather_int
  END INTERFACE

  INTERFACE p_scatterv
    MODULE PROCEDURE p_scatterv_real1D2D
  END INTERFACE

  INTERFACE p_gatherv
    MODULE PROCEDURE p_gatherv_int
    MODULE PROCEDURE p_gatherv_real2D1D
    MODULE PROCEDURE p_gatherv_real3D1D
    MODULE PROCEDURE p_gatherv_int2D1D
  END INTERFACE

  INTERFACE p_max
     MODULE PROCEDURE p_max_0d
     MODULE PROCEDURE p_max_int_0d
     MODULE PROCEDURE p_max_1d
     MODULE PROCEDURE p_max_int_1d
     MODULE PROCEDURE p_max_2d
     MODULE PROCEDURE p_max_3d
  END INTERFACE

  INTERFACE p_min
     MODULE PROCEDURE p_min_0d
     MODULE PROCEDURE p_min_int_0d
     MODULE PROCEDURE p_min_1d
     MODULE PROCEDURE p_min_int_1d
     MODULE PROCEDURE p_min_2d
     MODULE PROCEDURE p_min_3d
  END INTERFACE

  INTERFACE p_sum
     MODULE PROCEDURE p_sum_dp_0d
     MODULE PROCEDURE p_sum_dp_1d
     MODULE PROCEDURE p_sum_i8_1d
  END INTERFACE

  INTERFACE p_global_sum
     MODULE PROCEDURE p_global_sum_1d
  END INTERFACE

  INTERFACE p_field_sum
     MODULE PROCEDURE p_field_sum_1d
  END INTERFACE

  !> generic interface for MPI communication calls
  INTERFACE p_gather_field
    MODULE PROCEDURE p_gather_field_3d
    MODULE PROCEDURE p_gather_field_4d
    MODULE PROCEDURE p_gather_field_2d_int
    MODULE PROCEDURE p_gather_field_3d_int
  END INTERFACE

  INTERFACE p_allreduce_max
    MODULE PROCEDURE p_allreduce_max_int_1d
  END INTERFACE

  INTERFACE p_alltoall
    MODULE PROCEDURE p_alltoall_int
  END INTERFACE

CONTAINS

  !------------------------------------------------------------------------------
  REAL(dp) FUNCTION get_mpi_time()
#ifdef NOMPI
    get_mpi_time = 0.0_dp
#else
    get_mpi_time = MPI_Wtime()
#endif    
  END FUNCTION get_mpi_time
  !------------------------------------------------------------------------------
  
  !------------------------------------------------------------------------------
  REAL(dp) FUNCTION ellapsed_mpi_time(start_time)
    REAL(dp), INTENT(in) :: start_time
#ifdef NOMPI
    ellapsed_mpi_time = 0.0_dp
#else
    ellapsed_mpi_time = MPI_Wtime() - start_time
#endif    
  END FUNCTION ellapsed_mpi_time
  !------------------------------------------------------------------------------

  !------------------------------------------------------------------------------
  SUBROUTINE set_process_mpi_name(name)
    CHARACTER(len=*), INTENT(in) ::name
    
    process_mpi_name = TRIM(name)

  END SUBROUTINE set_process_mpi_name
  !------------------------------------------------------------------------------

  !------------------------------------------------------------------------------
  ! Rene: this has become obsolete and should be erased.
  ! It is used to idendify globally the proccess that caused an error
  INTEGER FUNCTION get_my_global_mpi_id()
    get_my_global_mpi_id = my_global_mpi_id
  END FUNCTION get_my_global_mpi_id
  !------------------------------------------------------------------------------

  !------------------------------------------------------------------------------
  INTEGER FUNCTION get_my_mpi_all_communicator()
    get_my_mpi_all_communicator = process_mpi_all_comm
  END FUNCTION get_my_mpi_all_communicator
  !------------------------------------------------------------------------------

  !------------------------------------------------------------------------------
  INTEGER FUNCTION get_my_mpi_all_comm_size()
    get_my_mpi_all_comm_size = process_mpi_all_size
  END FUNCTION get_my_mpi_all_comm_size
  !------------------------------------------------------------------------------
  
  !------------------------------------------------------------------------------
  INTEGER FUNCTION get_my_mpi_all_id()
    get_my_mpi_all_id = my_process_mpi_all_id
  END FUNCTION get_my_mpi_all_id
  !------------------------------------------------------------------------------
  
  !------------------------------------------------------------------------------
  INTEGER FUNCTION get_mpi_all_workroot_id()
    get_mpi_all_workroot_id = process_mpi_all_workroot_id
  END FUNCTION get_mpi_all_workroot_id
  !------------------------------------------------------------------------------
  
  !------------------------------------------------------------------------------
  INTEGER FUNCTION get_mpi_all_ioroot_id()
    get_mpi_all_ioroot_id = process_mpi_all_ioroot_id
  END FUNCTION get_mpi_all_ioroot_id
  !------------------------------------------------------------------------------

  !------------------------------------------------------------------------------
  INTEGER FUNCTION get_my_mpi_work_communicator()
    get_my_mpi_work_communicator = p_comm_work
  END FUNCTION get_my_mpi_work_communicator
  !------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  INTEGER FUNCTION get_my_mpi_work_comm_size()
    get_my_mpi_work_comm_size = p_n_work
  END FUNCTION get_my_mpi_work_comm_size
  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  INTEGER FUNCTION get_my_mpi_work_id()
    get_my_mpi_work_id = p_pe_work
  END FUNCTION get_my_mpi_work_id
  !------------------------------------------------------------------------------

 

  !------------------------------------------------------------------------------
  LOGICAL FUNCTION my_process_is_stdio()
#ifdef NOMPI
    my_process_is_stdio = .TRUE.
#else
    my_process_is_stdio = process_is_stdio
#endif
  END FUNCTION my_process_is_stdio
  !------------------------------------------------------------------------------

  !------------------------------------------------------------------------------
  LOGICAL FUNCTION my_process_is_io()
     my_process_is_io = (my_mpi_function == io_mpi_process)
  END FUNCTION my_process_is_io
  !------------------------------------------------------------------------------
  
  !------------------------------------------------------------------------------
  !>
  LOGICAL FUNCTION my_process_is_mpi_ioroot()
    my_process_is_mpi_ioroot = (my_process_mpi_all_id == process_mpi_all_ioroot_id)
  END FUNCTION my_process_is_mpi_ioroot
  !------------------------------------------------------------------------------

  !------------------------------------------------------------------------------
  LOGICAL FUNCTION my_process_is_mpi_test()
    my_process_is_mpi_test = (my_mpi_function == test_mpi_process)
  END FUNCTION my_process_is_mpi_test
  !------------------------------------------------------------------------------

  !------------------------------------------------------------------------------
  !>
  ! If is mpi parallel and not a test process
  !! Note: mpi i/o processes do not count is mpi paralell work
  !!        Only computational processes are checked for running in parallel
  LOGICAL FUNCTION my_process_is_mpi_parallel()
    my_process_is_mpi_parallel = process_is_mpi_parallel
  END FUNCTION my_process_is_mpi_parallel
  !------------------------------------------------------------------------------


  !------------------------------------------------------------------------------
  !>
  ! If is mpi all parallel
  LOGICAL FUNCTION my_process_is_mpi_all_parallel()
    my_process_is_mpi_all_parallel = (process_mpi_all_size > 1)
  END FUNCTION my_process_is_mpi_all_parallel
  !------------------------------------------------------------------------------
  
  !------------------------------------------------------------------------------
  !>
  ! If is mpi all parallel
  LOGICAL FUNCTION my_process_is_mpi_all_seq()
    my_process_is_mpi_all_seq = (process_mpi_all_size <= 1)
  END FUNCTION my_process_is_mpi_all_seq
  !------------------------------------------------------------------------------
  
  !------------------------------------------------------------------------------
  !>
  !! If is not mpi work paralellel or this is a test process
  !! returns true
  !! Note: mpi i/o processes do not count is mpi paralell work
  !!        Only computational processes are checked for running in parallel
  LOGICAL FUNCTION my_process_is_mpi_seq()
    my_process_is_mpi_seq = .NOT. process_is_mpi_parallel
  END FUNCTION my_process_is_mpi_seq
  !------------------------------------------------------------------------------
 
  !------------------------------------------------------------------------------
  !>
  LOGICAL FUNCTION my_process_is_mpi_workroot()
    my_process_is_mpi_workroot = (my_process_mpi_all_id == process_mpi_all_workroot_id)
  END FUNCTION my_process_is_mpi_workroot
  !------------------------------------------------------------------------------
  
 
  !------------------------------------------------------------------------------
  !>
  ! If is not mpi paralellel or is a test process
  LOGICAL FUNCTION run_is_global_mpi_parallel()
    run_is_global_mpi_parallel = is_global_mpi_parallel
  END FUNCTION run_is_global_mpi_parallel
  !------------------------------------------------------------------------------
  
  !------------------------------------------------------------------------------
  SUBROUTINE print_info_stderr (name, text)
    CHARACTER (len=*), INTENT(in) :: name
    CHARACTER (len=*), INTENT(in) :: text

    IF (my_process_is_stdio() ) THEN
      WRITE (nerr,'(a,a,a)') TRIM(name), ": ", TRIM(text)
    ENDIF
  
  END SUBROUTINE print_info_stderr
  !------------------------------------------------------------------------------
  
  
  !------------------------------------------------------------------------------
  SUBROUTINE finish (name, text)
    CHARACTER (len=*), INTENT(in) :: name
    CHARACTER (len=*), INTENT(in) :: text

    WRITE (nerr,'(i4.4,a,a,a,a)') my_global_mpi_id, ": ", TRIM(name), ": ", TRIM(text)
    CALL p_abort
      
  END SUBROUTINE finish
  !------------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !
  !> Pushes the communicator and proc0 onto the communicator stack.
  !! The communicator stack is needed for global sums if the processor
  !! set is split among different 1st level patches.
  
  SUBROUTINE push_glob_comm(comm, proc0)
    INTEGER, INTENT(IN) :: comm, proc0
    
    ! Safety check
    IF(comm_lev>=max_lev) &
      CALL finish('push_glob_comm','max_lev exceeded')
    
    comm_lev = comm_lev+1
    glob_comm(comm_lev) = comm
    comm_proc0(comm_lev) = proc0
    
  END SUBROUTINE push_glob_comm
  
  !-------------------------------------------------------------------------
  !
  !> Pops one level of the communicator stack
  
  SUBROUTINE pop_glob_comm()
    
    ! Safety check
    IF(comm_lev<=0) &
      CALL finish('pop_glob_comm','stack empty')
    
    comm_lev = comm_lev-1
    
  END SUBROUTINE pop_glob_comm
  
  
  !-------------------------------------------------------------------------
  !
  !> @return current top of communicator stack
  
  FUNCTION get_glob_proc0()
    INTEGER :: get_glob_proc0
    get_glob_proc0 = comm_proc0(comm_lev)
  END FUNCTION get_glob_proc0
  

  !------------------------------------------------------------------------------
  !>
  !! Sets the p_comm_input_bcast
  !! If comm_flag == null_comm_type then
  !!    only the test process reads and
  !!    no broadcast takes place
  !! Otherwise
  !!    the test or the root process reads
  !!    and broadcasts to the rest
  SUBROUTINE set_comm_input_bcast (comm_flag)
    INTEGER, INTENT(in), OPTIONAL:: comm_flag

    INTEGER :: comm_type

    comm_type = default_comm_type
    
    IF (PRESENT(comm_flag)) THEN
       comm_type = comm_flag
    ENDIF

#ifndef NOMPI
    SELECT CASE(comm_type)

    CASE(null_comm_type)
      IF(my_process_is_mpi_test()) THEN
        p_comm_input_bcast = MPI_COMM_SELF ! i.e. effectively no broadcast
      ELSE
        p_comm_input_bcast = MPI_COMM_NULL ! Must not be used!
      ENDIF
    
    CASE default
    
      IF (my_process_is_io()) THEN
        ! I/O PEs never participate in reading
        p_comm_input_bcast = MPI_COMM_NULL
      ELSE      
        IF(is_mpi_test_run) THEN
          ! Test PE reads and broadcasts to workers
          p_comm_input_bcast = p_comm_work_test
        ELSE
          ! PE 0 reads and broadcasts
          p_comm_input_bcast = p_comm_work
        ENDIF
      ENDIF

    END SELECT
#endif
    
  END SUBROUTINE set_comm_input_bcast
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  SUBROUTINE set_mpi_work_communicators(p_test_run, l_test_openmp, num_io_procs)
    LOGICAL,INTENT(INOUT) :: p_test_run, l_test_openmp
    INTEGER,INTENT(INOUT) :: num_io_procs
    
!   !local variables
    INTEGER :: my_color, peer_comm, p_error
    CHARACTER(*), PARAMETER :: method_name = "set_mpi_work_communicators"


     
! check l_test_openmp
#ifndef _OPENMP
    IF (l_test_openmp) THEN
      CALL print_info_stderr(method_name, &
        & 'l_test_openmp has no effect if the model is compiled without OpenMP support')
      CALL print_info_stderr(method_name, &
        & '--> l_test_openmp set to .FALSE.')
      l_test_openmp = .FALSE.
    END IF
#endif

    ! check p_test_run and num_io_procs
#ifdef NOMPI
    ! Unconditionally set p_test_run to .FALSE. and num_io_procs to 0,
    ! all other variables are already set correctly
    IF (p_test_run) THEN
      CALL print_info_stderr(method_name, &
       & 'p_test_run has no effect if the model is compiled with the NOMPI compiler directive')
      CALL print_info_stderr(method_name, &
       & '--> p_test_run set to .FALSE.')
      p_test_run = .FALSE.
    END IF
    IF (num_io_procs /= 0) THEN
      CALL print_info_stderr(method_name, &
       & 'num_io_procs has no effect if the model is compiled with the NOMPI compiler directive')
      CALL print_info_stderr(method_name, &
       & '--> num_io_procs set to 0')
      num_io_procs = 0
    END IF

    ! set the sequential values
    p_work_pe0 = 0
    num_io_procs = 0
    num_work_procs = 1

#else

    ! A run on 1 PE is never a verification run,
    ! correct this if the user should set it differently    
    IF (process_mpi_all_size < 2) THEN
      IF (p_test_run) THEN
        CALL print_info_stderr(method_name, &
            & 'p_test_run has no effect in seq run')
        CALL print_info_stderr(method_name, &
            & '--> p_test_run set to .FALSE.')
        p_test_run = .FALSE.
      ENDIF
      IF (num_io_procs > 0) THEN
        CALL print_info_stderr(method_name, &
            & 'num_io_procs cannot be > 0 in seq run')
        CALL print_info_stderr(method_name, &
            & '--> num_io_procs set to 0')
        num_io_procs = 0
      ENDIF
    ENDIF
    IF(num_io_procs < 0) num_io_procs = 0 ! for safety only  
    
    ! -----------------------------------------
    ! Set if test
    IF(p_test_run) THEN
      num_test_procs = 1
      p_test_pe = 0
    ELSE
      num_test_procs = 0
      p_test_pe = -1
    ENDIF
        
    ! -----------------------------------------
    ! how many work processors?
    num_work_procs = process_mpi_all_size - num_test_procs - num_io_procs
    
    ! Check if there are sufficient PEs at all
    IF(num_work_procs < 1) THEN
      CALL finish(method_name, &
      & 'not enough processors for given values of p_test_run/num_io_procs')
    ELSE IF (p_test_run .AND. num_work_procs == 1) THEN
      CALL finish(method_name, &
      & 'running p_test_run with only 1 work processor does not make sense')
    ENDIF

    WRITE(message_text,'(3(a,i0))') 'Number of procs for test: ',num_test_procs, &
      & ', work: ',num_work_procs,', I/O: ',num_io_procs
    CALL print_info_stderr(method_name, message_text)

    ! Everything seems ok. Proceed to setup the communicators and ids
    ! Set up p_test_pe, p_work_pe0, p_io_pe0 which are identical on all PEs
    p_work_pe0 = num_test_procs
    p_io_pe0   = num_test_procs + num_work_procs   

    ! Set up p_n_work and p_pe_work which are NOT identical on all PEs
    IF(p_pe < p_work_pe0) THEN
      ! Test PE (if present)
      p_n_work  = 1          ! 1 PE in verification work group
      p_pe_work = 0          ! PE number within work group
    ELSE IF(p_pe < p_io_pe0) THEN
      ! Work PE
      p_n_work  = num_work_procs
      p_pe_work = p_pe - num_test_procs
    ELSE
      ! I/O PE (if present)
      p_n_work  = num_io_procs
      p_pe_work = p_pe - num_test_procs - num_work_procs
    ENDIF


    ! Set communicators
    ! =================

    ! Split communicator process_mpi_all_comm between test/work/io
    ! to get p_comm_work which is the communicator for
    ! usage WITHIN every group of the 3 different type
    IF(p_pe < p_work_pe0) THEN
      my_mpi_function = test_mpi_process
    ELSE IF(p_pe < p_io_pe0) THEN
      my_mpi_function = work_mpi_process
    ELSE
      my_mpi_function = io_mpi_process
    ENDIF

    CALL MPI_Comm_split(process_mpi_all_comm, my_mpi_function, p_pe, p_comm_work, p_error)

    ! Set p_comm_work_test, the communicator spanning work group and test PE
    IF(p_test_run) THEN
      IF(p_pe < p_io_pe0) THEN
        my_color = 1
      ELSE
        my_color = MPI_UNDEFINED ! p_comm_work_test must never be used on I/O PEs
      ENDIF

      CALL MPI_Comm_split(process_mpi_all_comm, my_color, p_pe, p_comm_work_test, p_error)
    ELSE
      ! If not a test run, p_comm_work_test must not be used at all
      p_comm_work_test = MPI_COMM_NULL
    ENDIF

    ! Set p_comm_work_io, the communicator spanning work group and I/O PEs
    IF(num_io_procs > 0) THEN
      IF(p_pe < p_work_pe0) THEN
        my_color = MPI_UNDEFINED ! p_comm_work_io must never be used on test PE
      ELSE
        my_color = 1
      ENDIF

      CALL MPI_Comm_split(process_mpi_all_comm, my_color, p_pe, p_comm_work_io, p_error)
    ELSE
      ! If no I/O PEs are present, p_comm_work_io must not be used at all
      p_comm_work_io = MPI_COMM_NULL
    ENDIF

    
!     The following is moved to set_comm_input_bcast
!     ! Set p_comm_input_bcast, the communicator for broadcasting the NetCDF input
!     IF(lrestore_states) THEN
!       ! NetCDF input is only read by the test pe and MUST NOT be broadcast
!       IF(p_pe == p_test_pe) THEN
!         p_comm_input_bcast = MPI_COMM_SELF ! i.e. effectively no broadcast
!       ELSE
!         p_comm_input_bcast = MPI_COMM_NULL ! Must not be used!
!       ENDIF
!     ELSE
!       IF(p_pe < p_io_pe0) THEN
!         IF(p_test_run) THEN
!           ! Test PE reads and broadcasts to workers
!           p_comm_input_bcast = p_comm_work_test
!         ELSE
!           ! PE 0 reads and broadcasts
!           p_comm_input_bcast = p_comm_work
!         ENDIF
!       ELSE
!         ! I/O PEs never participate in reading
!         p_comm_input_bcast = MPI_COMM_NULL
!       ENDIF
!     ENDIF

    ! Create Intercommunicator work PEs - I/O PEs

    ! From MPI-Report:
    ! Advice to users: We recommend using a dedicated peer communicator, such as a
    ! duplicate of MPI_COMM_WORLD, to avoid trouble with peer communicators.

    ! No idea what these troubles may be in reality, but let us follow the advice

    CALL MPI_Comm_dup(process_mpi_all_comm, peer_comm, p_error)

    IF(p_pe /= p_test_pe .AND. num_io_procs>0) THEN

      IF(p_pe < p_io_pe0) THEN
        CALL MPI_Intercomm_create(p_comm_work, 0, peer_comm, p_io_pe0, &
          &  1, p_comm_work_2_io, p_error)
      ELSE
        CALL MPI_Intercomm_create(p_comm_work, 0, peer_comm, p_work_pe0,&
          & 1, p_comm_work_2_io, p_error)
      ENDIF
    ELSE
      ! No Intercommunicator
      p_comm_work_2_io = MPI_COMM_NULL
    ENDIF


    ! if OpenMP is used, the test PE uses only 1 thread in order to check
    ! the correctness of the OpenMP implementation
    ! Currently the I/O PEs are also single threaded!
#ifdef _OPENMP
    IF (l_test_openmp .AND. p_pe == p_test_pe) CALL OMP_SET_NUM_THREADS(1)
    IF (p_pe >= p_io_pe0) CALL OMP_SET_NUM_THREADS(1)
#endif
  
#endif

    ! fill some derived variables
    process_mpi_all_workroot_id = p_work_pe0
    process_mpi_all_ioroot_id   = p_io_pe0
    process_mpi_all_test_id     = p_work_pe0-1
    process_mpi_io_size         = num_io_procs

    ! In case of test run, only the test process is stdio
    process_is_stdio = (my_process_mpi_all_id == process_mpi_stdio_id)
    process_is_mpi_parallel = (num_work_procs > 1)
    IF (my_process_is_mpi_test()) process_is_mpi_parallel = .false.

    ! still to be filled
!     process_mpi_local_comm  = process_mpi_all_comm
!     process_mpi_local_size  = process_mpi_all_size
!     my_process_mpi_local_id = my_process_mpi_all_id

    ! fill my  parameters
    is_mpi_test_run = p_test_run
    is_openmp_test_run = l_test_openmp

    ! fill other default
    CALL set_comm_input_bcast()    
    
  END SUBROUTINE set_mpi_work_communicators
  !-------------------------------------------------------------------------
  
  !------------------------------------------------------------------------------
  !>
  SUBROUTINE set_default_mpi_work_variables()
    
    ! fill some derived variables
    process_is_stdio = (my_process_mpi_all_id == process_mpi_stdio_id)
    process_is_mpi_parallel = (process_mpi_all_size > 1)
    my_mpi_function = work_mpi_process
    
    process_mpi_all_test_id     = -1
    process_mpi_all_workroot_id = 0
    process_mpi_io_size         = 0
    is_mpi_test_run = .false.
    is_openmp_test_run = .false.

!     process_mpi_local_comm  = process_mpi_all_comm
!     process_mpi_local_size  = process_mpi_all_size
!     my_process_mpi_local_id = my_process_mpi_all_id         

    ! set some of the old variables
    ! should be removed once the old variables are cleaned
    p_pe           = my_process_mpi_all_id
    p_io           = 0 
    num_test_procs = 0
    num_work_procs = process_mpi_all_size
    p_test_pe      = -1
    p_work_pe0     = 0
    p_io_pe0       = process_mpi_all_size    ! Number of I/O PE 0 within all PEs (process_mpi_all_size if no I/O PEs)
    p_n_work       = process_mpi_all_size
    p_pe_work      = my_process_mpi_all_id
    
    p_comm_work             = process_mpi_all_comm
    p_comm_input_bcast      = process_mpi_all_comm
    p_comm_work_io          = MPI_COMM_NULL
    p_comm_work_test        = MPI_COMM_NULL
    p_comm_work_2_io        = MPI_COMM_NULL

    ! print some info
    IF ( .NOT. process_is_mpi_parallel) THEN
      WRITE (nerr,'(a)')  'Single processor run.'
    ELSEIF (process_is_stdio) THEN
      WRITE (nerr,'(a,a,i6,a)') TRIM(process_mpi_name), &
        '  runs on ', process_mpi_all_size, ' mpi processes.'
    END IF
    
  END SUBROUTINE set_default_mpi_work_variables
  !------------------------------------------------------------------------------

  !------------------------------------------------------------------------------
  !>
  !! Splits the global communicator into this component's communicator
  !! Should be called before the component configuration
  SUBROUTINE split_global_mpi_communicator(component_no)
    INTEGER, INTENT(in) :: component_no

    INTEGER :: new_communicator
    LOGICAL             :: l_mpi_is_initialised
    CHARACTER(len=*), PARAMETER :: method_name = 'split_process_mpi_communicator'
#ifdef NOMPI
    RETURN
#else    
    !--------------------------------------------
    ! check if mpi is initialized
    CALL MPI_INITIALIZED(l_mpi_is_initialised, p_error)
    IF (p_error /= MPI_SUCCESS) THEN
      WRITE (nerr,'(a,a)') method_name, ' MPI_INITITIALIZED failed.'
      WRITE (nerr,'(a,i4)') ' Error =  ', p_error
      STOP
    END IF
    !--------------------------------------------
    ! split global_mpi_communicator 
    CALL MPI_Comm_split(global_mpi_communicator, component_no, my_global_mpi_id, &
      & new_communicator, p_error)
    IF (p_error /= MPI_SUCCESS) THEN
      WRITE (nerr,'(a,a)') method_name, ' MPI_Comm_split failed.'
      WRITE (nerr,'(a,i4)') ' Error =  ', p_error
      STOP
    END IF
    CALL set_process_mpi_communicator(new_communicator)
#endif

  END SUBROUTINE split_global_mpi_communicator
  !------------------------------------------------------------------------------

  !------------------------------------------------------------------------------
  !>
  !! Set this component's communicator
  !! The new communicator id is duplicated from the input communicator
  !! Should be called before the component configuration
  SUBROUTINE set_process_mpi_communicator(new_communicator)
    INTEGER, INTENT(in) :: new_communicator

    LOGICAL             :: l_mpi_is_initialised
    CHARACTER(len=*), PARAMETER :: method_name = 'set_process_mpi_communicator'

#ifdef NOMPI
    process_mpi_all_comm    = new_communicator
    process_mpi_all_size    = 1
    my_process_mpi_all_id   = 0
    
#else
    ! since here we define the process all communicator
    ! the work communicator is identical to the all communicator
    ! and the no test or i/o processes are present    
    CALL MPI_INITIALIZED(l_mpi_is_initialised, p_error)

    IF (p_error /= MPI_SUCCESS) THEN
      WRITE (nerr,'(a,a)') method_name, ' MPI_INITITIALIZED failed.'
      WRITE (nerr,'(a,i4)') ' Error =  ', p_error
      STOP
    END IF

    IF ( .NOT. l_mpi_is_initialised ) THEN
       WRITE (nerr,'(a,a)') method_name, &
         & ' MPI_Init or p_start needs to be called first.'
       STOP
    ENDIF

    IF ( process_mpi_all_comm /= MPI_COMM_NULL) THEN   
      ! free original communicator
      CALL MPI_COMM_FREE(process_mpi_all_comm, p_error)
      IF (p_error /= MPI_SUCCESS) THEN
        WRITE (nerr,'(a,a)') method_name, &
          & ' MPI_COMM_FREE failed. p_start needs to be called before.'
        WRITE (nerr,'(a,i4)') ' Error =  ', p_error
        CALL p_abort
      END IF
    ENDIF
    ! assign MPI communicator generated elsewhere to process_mpi_all_comm

    CALL MPI_COMM_DUP(new_communicator, process_mpi_all_comm, p_error)
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,a)') method_name,&
         & ' MPI_COMM_DUP failed for process_mpi_all_comm.'
       WRITE (nerr,'(a,i4)') ' Error =  ', p_error
       CALL p_abort
    END IF

    ! get local PE identification
    CALL MPI_COMM_RANK (process_mpi_all_comm, my_process_mpi_all_id, p_error)

    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,a)') method_name, ' MPI_COMM_RANK failed.'
       WRITE (nerr,'(a,i4)') ' Error =  ', p_error
       CALL p_abort
    ELSE
#ifdef __DEBUG__
       WRITE (nerr,'(a,a,i4,a)') method_name, ' my_process_mpi_all_id ', &
         & my_process_mpi_all_id, ' started.'
#endif
    END IF


    CALL MPI_COMM_SIZE (process_mpi_all_comm, process_mpi_all_size, p_error)

    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,a,i4,a)') method_name, ' PE: ',&
         & my_process_mpi_all_id, ' MPI_COMM_SIZE failed.'
       WRITE (nerr,'(a,i4)') ' Error =  ', p_error
       CALL p_abort
    END IF


!     CALL MPI_COMM_DUP(process_mpi_all_comm,p_comm_work,p_error)
!     IF (p_error /= MPI_SUCCESS) THEN
!        WRITE (nerr,'(a,a)') method_name, ' MPI_COMM_DUP failed for p_comm_work.'
!        WRITE (nerr,'(a,i4)') ' Error =  ', p_error
!        CALL p_abort
!     END IF
! 
!     CALL MPI_COMM_DUP(process_mpi_all_comm,p_comm_work_test,p_error)
!     IF (p_error /= MPI_SUCCESS) THEN
!        WRITE (nerr,'(a,a)') method_name, ' MPI_COMM_DUP failed for p_comm_work_test.'
!        WRITE (nerr,'(a,i4)') ' Error =  ', p_error
!        CALL p_abort
!     END IF
   
#endif

    CALL set_default_mpi_work_variables()
        
  END SUBROUTINE set_process_mpi_communicator
  !------------------------------------------------------------------------------

  !------------------------------------------------------------------------------
  SUBROUTINE start_mpi(global_name)

#ifdef _OPENMP
    USE mo_util_string, ONLY: toupper
#endif

#ifndef NOMPI
#if defined (__prism) && defined (use_comm_MPI1)
    USE mod_prism_proto, ONLY: prism_ok
#endif
#endif

    CHARACTER(len=*), INTENT(in), OPTIONAL :: global_name

    ! variables are required for determing I/O size in bytes of the defined
    ! KIND types for assigning the right MPI data types with the used kinds
    INTEGER     :: iig = 0
    INTEGER(i4) :: ii4 = 0_i4
    INTEGER(i8) :: ii8 = 0_i8
    REAL        :: rrg = 0.0
    REAL(sp)    :: rsp = 0.0_sp
    REAL(dp)    :: rdp = 0.0_dp

    CHARACTER(len=132) :: yname

    ! variables used for determing the OpenMP threads
    ! suitable as well for coupled models
#if (defined _OPENMP)
    CHARACTER(len=32) :: env_name
    CHARACTER(len=32) :: thread_num
    INTEGER :: env_threads, threads
    INTEGER :: global_no_of_threads
#ifndef NOMPI
    INTEGER :: provided
#endif
#ifndef __SX__
    ! status
    INTEGER :: istat
#else
    EXTERNAL :: getenv
#endif
#endif

#ifndef NOMPI
#if defined (__prism) && defined (use_comm_MPI1)
    INTEGER :: prism_model_number
    CHARACTER(len=132) :: prism_model_name

    EXTERNAL :: prism_abort_proto
    EXTERNAL :: prism_init_comp_proto
    EXTERNAL :: prism_get_localcomm_proto
#endif
#endif
    CHARACTER(len=*), PARAMETER :: method_name = 'start_mpi'

    
! #ifdef _OPENMP
!     global_no_of_threads = 1
! #endif

    
    ! start MPI
#ifndef NOMPI
#ifdef _OPENMP
#ifdef __MULTIPLE_MPI_THREADS
    CALL MPI_INIT_THREAD(MPI_THREAD_MULTIPLE,provided,p_error)
#else
    CALL MPI_INIT_THREAD(MPI_THREAD_FUNNELED,provided,p_error)
#endif
#else
    CALL MPI_INIT (p_error)
#endif
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,a)') method_name, ' MPI_INIT failed.'
       WRITE (nerr,'(a,i4)') ' Error =  ', p_error
       STOP
    END IF
    
#ifdef _OPENMP
    ! Check if MPI_INIT_THREAD returned at least MPI_THREAD_FUNNELED in "provided"
#ifdef __MULTIPLE_MPI_THREADS
    IF (provided < MPI_THREAD_MULTIPLE) THEN
       WRITE (nerr,'(a,a)') method_name, &
         & ' MPI_INIT_THREAD did not return desired level of thread support'
       WRITE (nerr,'(a,i0)') " provided: ", provided
       WRITE (nerr,'(a,i0)') " required: ", MPI_THREAD_MULTIPLE
       CALL MPI_Finalize(p_error)
       STOP
    END IF
#else
    IF (provided < MPI_THREAD_FUNNELED) THEN
       WRITE (nerr,'(a,a)') method_name, &
         & ' MPI_INIT_THREAD did not return desired level of thread support'
       WRITE (nerr,'(a,i0)') " provided: ", provided
       WRITE (nerr,'(a,i0)') " required: ", MPI_THREAD_FUNNELED
       CALL MPI_Finalize(p_error)
       STOP
    END IF
#endif
#endif

    process_mpi_all_comm = MPI_COMM_NULL
    IF (PRESENT(global_name)) THEN
      yname = TRIM(global_name)
    ELSE
      yname = '(unnamed)'
    END IF
    global_mpi_name = TRIM(yname)
    
    ! create communicator for this process alone before
    ! potentially joining MPI2
#if defined (__prism) && defined (use_comm_MPI1)

    prism_model_name = TRIM(yname)

    CALL prism_init_comp_proto (prism_model_number, TRIM(prism_model_name), &
         p_error)

    IF (p_error /= prism_ok) THEN
      WRITE (nerr,*) method_name, ' prism_init_comp_proto failed'
      CALL prism_abort_proto(prism_model_number, TRIM(yname),'abort1')
    ENDIF

    CALL prism_get_localcomm_proto(global_mpi_communicator, p_error)

    IF (p_error /= prism_ok) THEN
      WRITE (nerr,*) method_name, ' prism_get_localcomm_proto failed'
      CALL prism_abort_proto(prism_model_number, TRIM(yname),'abort2')
    ENDIF

#else

    CALL MPI_COMM_DUP (MPI_COMM_WORLD, global_mpi_communicator, p_error)

    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,a)') method_name, ' MPI_COMM_DUP failed.'
       WRITE (nerr,'(a,i4)') ' Error =  ', p_error
       CALL p_abort
    END IF

#endif

    ! get local PE identification
    CALL MPI_COMM_RANK (global_mpi_communicator, my_global_mpi_id, p_error)

    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,a)') method_name, ' MPI_COMM_RANK failed.'
       WRITE (nerr,'(a,i4)') ' Error =  ', p_error
       CALL p_abort
    END IF

    ! Information ...
    IF (my_global_mpi_id == 0) THEN
       WRITE (nerr,'(/,a,a,a)') ' ', &
            TRIM(yname), ' MPI interface runtime information:'
    END IF
    
#ifdef DEBUG
    WRITE (nerr,'(a,a,i4,a)') method_name, ' PE ', my_global_mpi_id, ' started.'
#endif

    ! get number of available PEs
    CALL MPI_COMM_SIZE (global_mpi_communicator, global_mpi_size, p_error)

    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,a,i4,a)') method_name ,' PE: ', my_global_mpi_id, &
       &  ' MPI_COMM_SIZE failed.'                                         
       WRITE (nerr,'(a,i4)') ' Error =  ', p_error
       CALL p_abort
    END IF

    ! for non blocking calls
    p_mrequest = p_request_alloc_size
    ALLOCATE(p_request(p_mrequest))
    p_irequest = 0
    
    ! lets check the available MPI version
    CALL MPI_GET_VERSION (version, subversion, p_error)

    IF (p_error /= MPI_SUCCESS) THEN
      WRITE (nerr,'(a,a)') method_name , ' MPI_GET_VERSION failed.'
      WRITE (nerr,'(a,i4)') ' Error =  ', p_error
      CALL p_abort
    END IF

    IF (my_global_mpi_id == 0) THEN
      WRITE (nerr,'(a,a,i0,a1,i0)') method_name, &
           '  Used MPI version: ', version, '.', subversion
    END IF

    IF (my_global_mpi_id == 0) THEN
      WRITE (nerr,'(a,a,a,i0,a)') method_name, &
        & TRIM(yname), ': Globally run on ',&
        & global_mpi_size, ' mpi processes.'
    END IF

   ! due to a possible circular dependency with mo_machine and other
    ! modules, we determine here locally the I/O size of the different
    ! kind types (assume 8 bit/byte. This is than used for determing
    ! the right MPI send/receive type parameters.
    CALL MPI_SIZEOF(iig, p_int_byte, p_error)
    CALL MPI_SIZEOF(ii4, p_int_i4_byte, p_error)
    CALL MPI_SIZEOF(ii8, p_int_i8_byte, p_error)
    CALL MPI_SIZEOF(rrg, p_real_byte, p_error)
    CALL MPI_SIZEOF(rsp, p_real_sp_byte, p_error)
    CALL MPI_SIZEOF(rdp, p_real_dp_byte, p_error)

    p_int     = MPI_INTEGER
    p_real    = MPI_REAL
    p_bool    = MPI_LOGICAL
    p_char    = MPI_CHARACTER

    CALL MPI_TYPE_MATCH_SIZE(MPI_TYPECLASS_REAL, p_real_sp_byte, p_real_sp, p_error)
    CALL MPI_TYPE_MATCH_SIZE(MPI_TYPECLASS_REAL, p_real_dp_byte, p_real_dp, p_error)
    CALL MPI_TYPE_MATCH_SIZE(MPI_TYPECLASS_INTEGER, p_int_i4_byte, p_int_i4, p_error)
    CALL MPI_TYPE_MATCH_SIZE(MPI_TYPECLASS_INTEGER, p_int_i8_byte, p_int_i8, p_error)


#ifdef __DEBUG__
    WRITE (nerr,'(/,a)')    ' MPI transfer sizes [bytes]:'
    WRITE (nerr,'(a,i4)') '  INTEGER generic:', p_int_byte
    WRITE (nerr,'(a,i4)') '  INTEGER 4 byte :', p_int_i4_byte
    WRITE (nerr,'(a,i4)') '  INTEGER 8 byte :', p_int_i8_byte
    WRITE (nerr,'(a,i4)') '  REAL generic   :', p_real_byte
    WRITE (nerr,'(a,i4)') '  REAL single    :', p_real_sp_byte
    WRITE (nerr,'(a,i4)') '  REAL double    :', p_real_dp_byte
#endif

! MPI ends here
#else

    WRITE (nerr,'(a,a)')  method_name, ' No MPI: Single processor run.'
    ! set defaults for sequential run
    global_mpi_communicator = MPI_COMM_NULL
    global_mpi_size  = 1        ! total number of processes in global world
    my_global_mpi_id = 0        ! process id in global world
    is_global_mpi_parallel = .false.
    process_mpi_name = 'uknown'    
    process_mpi_all_comm = MPI_COMM_NULL

#endif


    IF (global_mpi_size > 1) THEN
      is_global_mpi_parallel = .TRUE.
    ELSE
      global_mpi_size = 1 ! just in case we got wrong size
    ENDIF

#ifdef _OPENMP
    ! The number of threads, if varying, will be defined via
    ! namelists
    global_no_of_threads = omp_get_max_threads()
#ifndef NOMPI
    ! Make number of threads from environment available to all model PEs
!     CALL MPI_BCAST (global_no_of_threads, 1, MPI_INTEGER, 0, global_mpi_communicator, p_error)
#endif

     IF (my_global_mpi_id == 0) THEN

      IF (is_global_mpi_parallel) THEN
        WRITE (nerr,'(/,a,a)') method_name, &
          & ': Running globally hybrid OpenMP-MPI mode.'
      ELSE
        WRITE (nerr,'(/,a,a)') method_name,': Running globally OpenMP mode.'
      ENDIF
    ENDIF
    WRITE (nerr,'(a, a, i3, a, i3)') method_name,': PE:', my_global_mpi_id, &
      & ' global_no_of_threads is ', global_no_of_threads
    
#endif


    ! The number of threads, if varying, will be defined via
    ! namelists
! #ifdef _OPENMP
!     ! Expect that PE 0 did got the information of OMP_NUM_THREADS.
!     ! That might be wrong in the coupled case when the model is
!     ! started via MPI dynamic process creation. So we have to check
!     ! the environment variable too.
!     IF (my_global_mpi_id == 0) THEN
! 
!       IF (is_global_mpi_parallel) THEN
!         WRITE (nerr,'(/,a)') ' Running globally hybrid OpenMP-MPI mode.'
!       ELSE
!         WRITE (nerr,'(/,a)') ' Running globally OpenMP mode.'
!       ENDIF
! 
!       env_name = toupper(TRIM(yname)) // '_THREADS'
! #ifdef __SX__
!       CALL getenv(TRIM(env_name), thread_num)
! 
!       IF (thread_num /= ' ') THEN
! #else
!       CALL get_environment_variable(name=TRIM(env_name), value=thread_num, &
!           status=istat)
!       IF (istat == 0) THEN
! #endif
!          READ(thread_num,*) global_no_of_threads
!        ELSE
!          WRITE (nerr,'(1x,a,/)') ' Global number of OpenMP threads not given!'
!          WRITE (nerr,'(1x,a,a,a,/,1x,a,a,a)') &
!              ' Environment variable ', TRIM(env_name), ' either not set,', &
!              ' or not available to ', TRIM(yname), ' root PE.'
!          global_no_of_threads = omp_get_max_threads()
!        ENDIF
!     ENDIF ! (my_global_mpi_id == 0)
! 
! #ifndef NOMPI
!     ! Make number of threads from environment available to all model PEs
!     CALL MPI_BCAST (global_no_of_threads, 1, MPI_INTEGER, 0, global_mpi_communicator, p_error)
! #endif
!     ! Inform on OpenMP thread usage
! !     CALL OMP_SET_NUM_THREADS(threads)
! !     threads = OMP_GET_MAX_THREADS()
! 
!     IF (my_global_mpi_id == 0) THEN
!        WRITE (nerr,*)
!        WRITE (nerr,'(1x,a,i3)') ' global_no_of_threads is ', global_no_of_threads
!     ENDIF
! #endif


    ! by default, the global communicator is the process communicator
    CALL set_process_mpi_name(global_mpi_name)
    CALL set_process_mpi_communicator(global_mpi_communicator)
    
  END SUBROUTINE start_mpi
  !------------------------------------------------------------------------------

  !------------------------------------------------------------------------------
  SUBROUTINE p_stop

    ! finish MPI and clean up all PEs

#ifndef NOMPI
    ! to prevent abort due to unfinished communication
    CALL p_barrier(process_mpi_all_comm)

    CALL MPI_FINALIZE (p_error)

    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a)') ' MPI_FINALIZE failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
    process_is_mpi_parallel = .FALSE.
    DEALLOCATE(p_request)
#endif

  END SUBROUTINE p_stop
  !------------------------------------------------------------------------------

  !------------------------------------------------------------------------------
  SUBROUTINE p_abort

    ! this routine should be used instead of abort, util_abort() or STOP
    ! in all routines for proper clean up of all PEs

#ifndef __STANDALONE
    EXTERNAL :: util_exit
#endif

#ifndef NOMPI
    CALL MPI_ABORT (MPI_COMM_WORLD, 0, p_error)

    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a)') ' MPI_ABORT failed.'
       WRITE (nerr,'(a,i4)') ' Error =  ', p_error
       STOP
    END IF
#else
#ifndef __STANDALONE
    CALL util_exit(1)
#else
    STOP 'mo_mpi: p_abort ..'
#endif
#endif

  END SUBROUTINE p_abort
  !------------------------------------------------------------------------------

  ! communicator set up
  SUBROUTINE p_set_communicator (nproca, nprocb, mapmesh, debug_parallel)

    INTEGER, INTENT(in) :: nproca, nprocb
    INTEGER, INTENT(in) :: mapmesh(0:,0:)
    INTEGER, INTENT(in) :: debug_parallel

#ifndef NOMPI
    INTEGER :: all_debug_pes(SIZE(mapmesh))

    INTEGER :: group_world, group_a, group_b, group_d
    INTEGER :: p_communicator_tmp

    INTEGER :: n, members

    INTEGER :: ranks(1) = 0

    ! first set global group

    CALL MPI_COMM_GROUP (process_mpi_all_comm, group_world, p_error)

    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' PE: ', my_process_mpi_all_id, &
         & ' MPI_COMM_GROUP failed.'
       WRITE (nerr,'(a,i4)') ' Error =  ', p_error
       CALL p_abort
    END IF

    ! communicator is process_mpi_all_comm

    IF (debug_parallel >= 0 ) THEN

       CALL MPI_GROUP_INCL (group_world, 1, ranks, group_d, p_error)

       IF (p_error /= MPI_SUCCESS) THEN
          WRITE (nerr,'(a,i4,a)') ' PE: ', my_process_mpi_all_id, &
            &' MPI_GROUP_INCL failed.'
          WRITE (nerr,'(a,i4)') ' Error =  ', p_error
          CALL p_abort
       END IF

       CALL MPI_COMM_CREATE (process_mpi_all_comm, group_d, p_communicator_tmp, &
            p_error)

       IF (p_error /= MPI_SUCCESS) THEN
          WRITE (nerr,'(a,i4,a)') ' PE: ', my_process_mpi_all_id, &
            &' MPI_COMM_CREATE failed.'
          WRITE (nerr,'(a,i4)') ' Error =  ', p_error
          CALL p_abort
       END IF

       IF (my_process_mpi_all_id == 0) p_communicator_d = p_communicator_tmp

       DO n = 1, SIZE(mapmesh)
          all_debug_pes(n) = n
       END DO

       CALL MPI_GROUP_INCL (group_world, SIZE(mapmesh), all_debug_pes, &
            group_d, p_error)

       IF (p_error /= MPI_SUCCESS) THEN
          WRITE (nerr,'(a,i4,a)') ' PE: ', my_process_mpi_all_id, &
            & ' MPI_GROUP_INCL failed.'
          WRITE (nerr,'(a,i4)') ' Error =  ', p_error
          CALL p_abort
       END IF

       CALL MPI_COMM_CREATE (process_mpi_all_comm, group_d, p_communicator_tmp, &
            p_error)

       IF (p_error /= MPI_SUCCESS) THEN
          WRITE (nerr,'(a,i4,a)') ' PE: ', my_process_mpi_all_id, &
            & ' MPI_COMM_CREATE failed.'
          WRITE (nerr,'(a,i4)') ' Error =  ', p_error
          CALL p_abort
       END IF

       IF (my_process_mpi_all_id /= 0) p_communicator_d = p_communicator_tmp

    ELSE
       p_communicator_d = process_mpi_all_comm
    END IF

    DO n = 0, nproca-1
       members = nprocb
       CALL MPI_GROUP_INCL (group_world, members, mapmesh(:,n), group_a, &
            p_error)

       IF (p_error /= MPI_SUCCESS) THEN
          WRITE (nerr,'(a,i4,a)') ' PE: ', my_process_mpi_all_id, &
            & ' MPI_GROUP_INCL failed.'
          WRITE (nerr,'(a,i4)') ' Error =  ', p_error
          CALL p_abort
       END IF

       CALL MPI_COMM_CREATE (process_mpi_all_comm, group_a, p_communicator_tmp, &
            p_error)

       IF (p_error /= MPI_SUCCESS) THEN
          WRITE (nerr,'(a,i4,a)') ' PE: ', my_process_mpi_all_id, &
            & ' MPI_COMM_CREATE failed.'
          WRITE (nerr,'(a,i4)') ' Error =  ', p_error
          CALL p_abort
       END IF
       IF(p_communicator_tmp/=MPI_COMM_NULL) &
         p_communicator_a = p_communicator_tmp

    END DO

    ! create groups for set Bs

    DO n = 0, nprocb-1
       members = nproca
       CALL MPI_GROUP_INCL (group_world, members, mapmesh(n,:), group_b, &
            p_error)

       IF (p_error /= MPI_SUCCESS) THEN
          WRITE (nerr,'(a,i4,a)') ' PE: ', my_process_mpi_all_id, &
            & ' MPI_GROUP_INCL failed.'
          WRITE (nerr,'(a,i4)') ' Error =  ', p_error
          CALL p_abort
       END IF

       CALL MPI_COMM_CREATE (process_mpi_all_comm, group_b, p_communicator_tmp, &
            p_error)

       IF (p_error /= MPI_SUCCESS) THEN
          WRITE (nerr,'(a,i4,a)') ' PE: ', my_process_mpi_all_id, &
            & ' MPI_COMM_CREATE failed.'
          WRITE (nerr,'(a,i4)') ' Error =  ', p_error
          CALL p_abort
       END IF
       IF(p_communicator_tmp/=MPI_COMM_NULL) &
         p_communicator_b = p_communicator_tmp

    END DO

    CALL MPI_BARRIER (process_mpi_all_comm, p_error)

    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' PE: ', my_process_mpi_all_id, '&
         & MPI_BARRIER failed.'
       WRITE (nerr,'(a,i4)') ' Error =  ', p_error
       CALL p_abort
    END IF

    IF (debug_parallel >= 0 .AND. my_process_mpi_all_id == 0) THEN
      p_communicator_a = p_communicator_d
      p_communicator_b = p_communicator_d
    ENDIF

#ifdef DEBUG
    WRITE (nerr,'(a,i4,a,3i8)') &
         'p_set_communicator on PE ', my_process_mpi_all_id, ': ', &
         p_communicator_d, &
         p_communicator_a, &
         p_communicator_b
#endif
#endif
  END SUBROUTINE p_set_communicator

!=========================================================================

  ! send implementation

  SUBROUTINE p_send_real (t_buffer, p_destination, p_tag, p_count, comm)

    REAL (dp), INTENT(in) :: t_buffer
    INTEGER,   INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_SEND (t_buffer, p_count, p_real_dp, p_destination, p_tag, &
            p_comm, p_error)
    ELSE
       CALL MPI_SEND (t_buffer, 1, p_real_dp, p_destination, p_tag, &
            p_comm, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', my_process_mpi_all_id, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_send_real

  SUBROUTINE p_send_real_1d (t_buffer, p_destination, p_tag, p_count, comm)

    REAL (dp), INTENT(in) :: t_buffer(:)
    INTEGER,   INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_SEND (t_buffer, p_count, p_real_dp, p_destination, p_tag, &
            p_comm, p_error)
    ELSE
       CALL MPI_SEND (t_buffer, SIZE(t_buffer), p_real_dp, p_destination, p_tag, &
            p_comm, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', my_process_mpi_all_id, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_send_real_1d

  SUBROUTINE p_send_real_2d (t_buffer, p_destination, p_tag, p_count, comm)

    REAL (dp), INTENT(in) :: t_buffer(:,:)
    INTEGER,   INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_SEND (t_buffer, p_count, p_real_dp, p_destination, p_tag, &
            p_comm, p_error)
    ELSE
       CALL MPI_SEND (t_buffer, SIZE(t_buffer), p_real_dp, p_destination, p_tag, &
            p_comm, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', my_process_mpi_all_id, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_send_real_2d

  SUBROUTINE p_send_real_3d (t_buffer, p_destination, p_tag, p_count, comm)

    REAL (dp), INTENT(in) :: t_buffer(:,:,:)
    INTEGER,   INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_SEND (t_buffer, p_count, p_real_dp, p_destination, p_tag, &
            p_comm, p_error)
    ELSE
       CALL MPI_SEND (t_buffer, SIZE(t_buffer), p_real_dp, p_destination, p_tag, &
            p_comm, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', my_process_mpi_all_id, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_send_real_3d

  SUBROUTINE p_send_real_4d (t_buffer, p_destination, p_tag, p_count, comm)

    REAL (dp), INTENT(in) :: t_buffer(:,:,:,:)
    INTEGER,   INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_SEND (t_buffer, p_count, p_real_dp, p_destination, p_tag, &
            p_comm, p_error)
    ELSE
       CALL MPI_SEND (t_buffer, SIZE(t_buffer), p_real_dp, p_destination, p_tag, &
            p_comm, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', my_process_mpi_all_id, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_send_real_4d

  SUBROUTINE p_send_real_5d (t_buffer, p_destination, p_tag, p_count, comm)

    REAL (dp), INTENT(in) :: t_buffer(:,:,:,:,:)
    INTEGER,   INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_SEND (t_buffer, p_count, p_real_dp, p_destination, p_tag, &
            p_comm, p_error)
    ELSE
       CALL MPI_SEND (t_buffer, SIZE(t_buffer), p_real_dp, p_destination, p_tag, &
            p_comm, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', my_process_mpi_all_id, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_send_real_5d

  SUBROUTINE p_send_int (t_buffer, p_destination, p_tag, p_count, comm)

    INTEGER, INTENT(in) :: t_buffer
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_SEND (t_buffer, p_count, p_int, p_destination, p_tag, &
            p_comm, p_error)
    ELSE
       CALL MPI_SEND (t_buffer, 1, p_int, p_destination, p_tag, &
            p_comm, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', my_process_mpi_all_id, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_send_int

  SUBROUTINE p_send_int_1d (t_buffer, p_destination, p_tag, p_count, comm)

    INTEGER, INTENT(in) :: t_buffer(:)
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_SEND (t_buffer, p_count, p_int, p_destination, p_tag, &
            p_comm, p_error)
    ELSE
       CALL MPI_SEND (t_buffer, SIZE(t_buffer), p_int, p_destination, p_tag, &
            p_comm, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', my_process_mpi_all_id, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_send_int_1d

  SUBROUTINE p_send_int_2d (t_buffer, p_destination, p_tag, p_count, comm)

    INTEGER, INTENT(in) :: t_buffer(:,:)
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_SEND (t_buffer, p_count, p_int, p_destination, p_tag, &
            p_comm, p_error)
    ELSE
       CALL MPI_SEND (t_buffer, SIZE(t_buffer), p_int, p_destination, p_tag, &
            p_comm, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', my_process_mpi_all_id, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_send_int_2d

  SUBROUTINE p_send_int_3d (t_buffer, p_destination, p_tag, p_count, comm)

    INTEGER, INTENT(in) :: t_buffer(:,:,:)
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_SEND (t_buffer, p_count, p_int, p_destination, p_tag, &
            p_comm, p_error)
    ELSE
       CALL MPI_SEND (t_buffer, SIZE(t_buffer), p_int, p_destination, p_tag, &
            p_comm, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', my_process_mpi_all_id, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_send_int_3d

  SUBROUTINE p_send_int_4d (t_buffer, p_destination, p_tag, p_count, comm)

    INTEGER, INTENT(in) :: t_buffer(:,:,:,:)
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_SEND (t_buffer, p_count, p_int, p_destination, p_tag, &
            p_comm, p_error)
    ELSE
       CALL MPI_SEND (t_buffer, SIZE(t_buffer), p_int, p_destination, p_tag, &
            p_comm, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', my_process_mpi_all_id, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_send_int_4d


  SUBROUTINE p_send_bool (t_buffer, p_destination, p_tag, p_count, comm)

    LOGICAL, INTENT(in) :: t_buffer
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_SEND (t_buffer, p_count, p_bool, p_destination, p_tag, &
            p_comm, p_error)
    ELSE
       CALL MPI_SEND (t_buffer, 1, p_bool, p_destination, p_tag, &
            p_comm, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', my_process_mpi_all_id, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_send_bool

  SUBROUTINE p_send_bool_1d (t_buffer, p_destination, p_tag, p_count, comm)

    LOGICAL, INTENT(in) :: t_buffer(:)
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_SEND (t_buffer, p_count, p_bool, p_destination, p_tag, &
            p_comm, p_error)
    ELSE
       CALL MPI_SEND (t_buffer, SIZE(t_buffer), p_bool, p_destination, p_tag, &
            p_comm, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', my_process_mpi_all_id, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_send_bool_1d

  SUBROUTINE p_send_bool_2d (t_buffer, p_destination, p_tag, p_count, comm)

    LOGICAL, INTENT(in) :: t_buffer(:,:)
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_SEND (t_buffer, p_count, p_bool, p_destination, p_tag, &
            p_comm, p_error)
    ELSE
       CALL MPI_SEND (t_buffer, SIZE(t_buffer), p_bool, p_destination, p_tag, &
            p_comm, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', my_process_mpi_all_id, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_send_bool_2d

  SUBROUTINE p_send_bool_3d (t_buffer, p_destination, p_tag, p_count, comm)

    LOGICAL, INTENT(in) :: t_buffer(:,:,:)
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_SEND (t_buffer, p_count, p_bool, p_destination, p_tag, &
            p_comm, p_error)
    ELSE
       CALL MPI_SEND (t_buffer, SIZE(t_buffer), p_bool, p_destination, p_tag, &
            p_comm, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', my_process_mpi_all_id, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_send_bool_3d

  SUBROUTINE p_send_bool_4d (t_buffer, p_destination, p_tag, p_count, comm)

    LOGICAL, INTENT(in) :: t_buffer(:,:,:,:)
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_SEND (t_buffer, p_count, p_bool, p_destination, p_tag, &
            p_comm, p_error)
    ELSE
       CALL MPI_SEND (t_buffer, SIZE(t_buffer), p_bool, p_destination, p_tag, &
            p_comm, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', my_process_mpi_all_id, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_send_bool_4d

  SUBROUTINE p_send_char (t_buffer, p_destination, p_tag, p_count, comm)

    CHARACTER(len=*),  INTENT(in) :: t_buffer
    INTEGER,           INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_SEND (t_buffer, p_count, p_char, p_destination, p_tag, &
            p_comm, p_error)
    ELSE
       CALL MPI_SEND (t_buffer, LEN(t_buffer), p_char, p_destination, p_tag, &
            p_comm, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', my_process_mpi_all_id, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_send_char

! non-blocking sends

  SUBROUTINE p_inc_request
    INTEGER, ALLOCATABLE :: tmp(:)

#ifndef NOMPI
    p_irequest = p_irequest + 1
    IF (p_irequest > p_mrequest) THEN
      ALLOCATE(tmp(p_mrequest))
      tmp(:) = p_request(:)
      DEALLOCATE(p_request)
      ALLOCATE(p_request(p_mrequest+p_request_alloc_size))
      p_request(1:p_mrequest) = tmp(:)
      p_mrequest = p_mrequest+p_request_alloc_size
      DEALLOCATE(tmp)
    ENDIF
#endif

  END SUBROUTINE p_inc_request

  SUBROUTINE p_isend_real (t_buffer, p_destination, p_tag, p_count, comm)

    REAL (dp), INTENT(inout) :: t_buffer
    INTEGER,   INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_ISEND (t_buffer, p_count, p_real_dp, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_ISEND (t_buffer, 1, p_real_dp, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_ISEND from ', my_process_mpi_all_id, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_isend_real

  SUBROUTINE p_isend_real_1d (t_buffer, p_destination, p_tag, p_count, comm, request)

    REAL (dp), INTENT(inout) :: t_buffer(:)
    INTEGER,   INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
    INTEGER,   INTENT(INOUT), OPTIONAL :: request
#ifndef NOMPI
    INTEGER :: p_comm, icount, out_request

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    icount = SIZE(t_buffer)
    IF (PRESENT(p_count))  icount = p_count

    CALL p_inc_request
    CALL MPI_ISEND (t_buffer, icount, p_real_dp, p_destination, p_tag, &
      &             p_comm, out_request, p_error)

    IF (PRESENT(request)) THEN
      request               = out_request
    ELSE
      p_request(p_irequest) = out_request
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_ISEND from ', my_process_mpi_all_id, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_isend_real_1d

  SUBROUTINE p_isend_real_2d (t_buffer, p_destination, p_tag, p_count, comm)

    REAL (dp), INTENT(inout) :: t_buffer(:,:)
    INTEGER,   INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_ISEND (t_buffer, p_count, p_real_dp, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_ISEND (t_buffer, SIZE(t_buffer), p_real_dp, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_ISEND from ', my_process_mpi_all_id, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_isend_real_2d

  SUBROUTINE p_isend_real_3d (t_buffer, p_destination, p_tag, p_count, comm)

    REAL (dp), INTENT(inout) :: t_buffer(:,:,:)
    INTEGER,   INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_ISEND (t_buffer, p_count, p_real_dp, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_ISEND (t_buffer, SIZE(t_buffer), p_real_dp, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_ISEND from ', my_process_mpi_all_id, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_isend_real_3d

  SUBROUTINE p_isend_real_4d (t_buffer, p_destination, p_tag, p_count, comm)

    REAL (dp), INTENT(inout) :: t_buffer(:,:,:,:)
    INTEGER,   INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_ISEND (t_buffer, p_count, p_real_dp, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_ISEND (t_buffer, SIZE(t_buffer), p_real_dp, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_ISEND from ', my_process_mpi_all_id, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_isend_real_4d

  SUBROUTINE p_isend_real_5d (t_buffer, p_destination, p_tag, p_count, comm)

    REAL (dp), INTENT(inout) :: t_buffer(:,:,:,:,:)
    INTEGER,   INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_ISEND (t_buffer, p_count, p_real_dp, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_ISEND (t_buffer, SIZE(t_buffer), p_real_dp, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_ISEND from ', my_process_mpi_all_id, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_isend_real_5d

  SUBROUTINE p_isend_int (t_buffer, p_destination, p_tag, p_count, comm)

    INTEGER, INTENT(inout) :: t_buffer
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_ISEND (t_buffer, p_count, p_int, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_ISEND (t_buffer, 1, p_int, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_ISEND from ', my_process_mpi_all_id, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_isend_int

  SUBROUTINE p_isend_int_1d (t_buffer, p_destination, p_tag, p_count, comm)

    INTEGER, INTENT(inout) :: t_buffer(:)
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_ISEND (t_buffer, p_count, p_int, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_ISEND (t_buffer, SIZE(t_buffer), p_int, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_ISEND from ', my_process_mpi_all_id, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_isend_int_1d

  SUBROUTINE p_isend_int_2d (t_buffer, p_destination, p_tag, p_count, comm)

    INTEGER, INTENT(inout) :: t_buffer(:,:)
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_ISEND (t_buffer, p_count, p_int, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_ISEND (t_buffer, SIZE(t_buffer), p_int, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_ISEND from ', my_process_mpi_all_id, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_isend_int_2d

  SUBROUTINE p_isend_int_3d (t_buffer, p_destination, p_tag, p_count, comm)

    INTEGER, INTENT(inout) :: t_buffer(:,:,:)
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_ISEND (t_buffer, p_count, p_int, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_ISEND (t_buffer, SIZE(t_buffer), p_int, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_ISEND from ', my_process_mpi_all_id, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_isend_int_3d

  SUBROUTINE p_isend_int_4d (t_buffer, p_destination, p_tag, p_count, comm)

    INTEGER, INTENT(inout) :: t_buffer(:,:,:,:)
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_ISEND (t_buffer, p_count, p_int, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_ISEND (t_buffer, SIZE(t_buffer), p_int, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_ISEND from ', my_process_mpi_all_id, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_isend_int_4d


  SUBROUTINE p_isend_bool (t_buffer, p_destination, p_tag, p_count, comm)

    LOGICAL, INTENT(inout) :: t_buffer
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_ISEND (t_buffer, p_count, p_bool, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_ISEND (t_buffer, 1, p_bool, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_ISEND from ', my_process_mpi_all_id, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_isend_bool

  SUBROUTINE p_isend_bool_1d (t_buffer, p_destination, p_tag, p_count, comm)

    LOGICAL, INTENT(inout) :: t_buffer(:)
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_ISEND (t_buffer, p_count, p_bool, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_ISEND (t_buffer, SIZE(t_buffer), p_bool, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_ISEND from ', my_process_mpi_all_id, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_isend_bool_1d

  SUBROUTINE p_isend_bool_2d (t_buffer, p_destination, p_tag, p_count, comm)

    LOGICAL, INTENT(inout) :: t_buffer(:,:)
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_ISEND (t_buffer, p_count, p_bool, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_ISEND (t_buffer, SIZE(t_buffer), p_bool, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_ISEND from ', my_process_mpi_all_id, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_isend_bool_2d

  SUBROUTINE p_isend_bool_3d (t_buffer, p_destination, p_tag, p_count, comm)

    LOGICAL, INTENT(inout) :: t_buffer(:,:,:)
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_ISEND (t_buffer, p_count, p_bool, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_ISEND (t_buffer, SIZE(t_buffer), p_bool, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_ISEND from ', my_process_mpi_all_id, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_isend_bool_3d

  SUBROUTINE p_isend_bool_4d (t_buffer, p_destination, p_tag, p_count, comm)

    LOGICAL, INTENT(inout) :: t_buffer(:,:,:,:)
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_ISEND (t_buffer, p_count, p_bool, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_ISEND (t_buffer, SIZE(t_buffer), p_bool, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_ISEND from ', my_process_mpi_all_id, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_isend_bool_4d

  SUBROUTINE p_isend_char (t_buffer, p_destination, p_tag, p_count, comm)

    CHARACTER(len=*),  INTENT(inout) :: t_buffer
    INTEGER,           INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_ISEND (t_buffer, p_count, p_char, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_ISEND (t_buffer, LEN(t_buffer), p_char, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_ISEND from ', my_process_mpi_all_id, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_isend_char

  ! recv implementation

  SUBROUTINE p_recv_real (t_buffer, p_source, p_tag, p_count, comm)

    REAL (dp), INTENT(out) :: t_buffer
    INTEGER,   INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_RECV (t_buffer, p_count, p_real_dp, p_source, p_tag, &
            p_comm, p_status, p_error)
    ELSE
       CALL MPI_RECV (t_buffer, 1, p_real_dp, p_source, p_tag, &
            p_comm, p_status, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_RECV on ', my_process_mpi_all_id, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_recv_real

  SUBROUTINE p_recv_real_1d (t_buffer, p_source, p_tag, p_count, comm, displs)

    REAL (dp), INTENT(out) :: t_buffer(:)
    INTEGER,   INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm, displs
#ifndef NOMPI
    INTEGER :: p_comm, icount, idispls

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    icount = SIZE(t_buffer)
    IF (PRESENT(p_count))  icount = p_count

    idispls = 1
    IF (PRESENT(displs)) idispls = displs

    CALL MPI_RECV (t_buffer(idispls), icount, p_real_dp, p_source, p_tag, &
      &            p_comm, p_status, p_error)

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_RECV on ', my_process_mpi_all_id, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_recv_real_1d

  SUBROUTINE p_recv_real_2d (t_buffer, p_source, p_tag, p_count, comm)

    REAL (dp), INTENT(out) :: t_buffer(:,:)
    INTEGER,   INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_RECV (t_buffer, p_count, p_real_dp, p_source, p_tag, &
            p_comm, p_status, p_error)
    ELSE
       CALL MPI_RECV (t_buffer, SIZE(t_buffer), p_real_dp, p_source, p_tag, &
            p_comm, p_status, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_RECV on ', my_process_mpi_all_id, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_recv_real_2d

  SUBROUTINE p_recv_real_3d (t_buffer, p_source, p_tag, p_count, comm)

    REAL (dp), INTENT(out) :: t_buffer(:,:,:)
    INTEGER,   INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_RECV (t_buffer, p_count, p_real_dp, p_source, p_tag, &
            p_comm, p_status, p_error)
    ELSE
       CALL MPI_RECV (t_buffer, SIZE(t_buffer), p_real_dp, p_source, p_tag, &
            p_comm, p_status, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_RECV on ', my_process_mpi_all_id, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_recv_real_3d

  SUBROUTINE p_recv_real_4d (t_buffer, p_source, p_tag, p_count, comm)

    REAL (dp), INTENT(out) :: t_buffer(:,:,:,:)
    INTEGER,   INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_RECV (t_buffer, p_count, p_real_dp, p_source, p_tag, &
            p_comm, p_status, p_error)
    ELSE
       CALL MPI_RECV (t_buffer, SIZE(t_buffer), p_real_dp, p_source, p_tag, &
            p_comm, p_status, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_RECV on ', my_process_mpi_all_id, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_recv_real_4d

  SUBROUTINE p_recv_real_5d (t_buffer, p_source, p_tag, p_count, comm)

    REAL (dp), INTENT(out) :: t_buffer(:,:,:,:,:)
    INTEGER,   INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_RECV (t_buffer, p_count, p_real_dp, p_source, p_tag, &
            p_comm, p_status, p_error)
    ELSE
       CALL MPI_RECV (t_buffer, SIZE(t_buffer), p_real_dp, p_source, p_tag, &
            p_comm, p_status, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_RECV on ', my_process_mpi_all_id, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_recv_real_5d

  SUBROUTINE p_recv_int (t_buffer, p_source, p_tag, p_count, comm)

    INTEGER, INTENT(out) :: t_buffer
    INTEGER, INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_RECV (t_buffer, p_count, p_int, p_source, p_tag, &
            p_comm, p_status, p_error)
    ELSE
       CALL MPI_RECV (t_buffer, 1, p_int, p_source, p_tag, &
            p_comm, p_status, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_RECV on ', my_process_mpi_all_id, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_recv_int

  SUBROUTINE p_recv_int_1d (t_buffer, p_source, p_tag, p_count, comm)

    INTEGER, INTENT(out) :: t_buffer(:)
    INTEGER, INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_RECV (t_buffer, p_count, p_int, p_source, p_tag, &
            p_comm, p_status, p_error)
    ELSE
       CALL MPI_RECV (t_buffer, SIZE(t_buffer), p_int, p_source, p_tag, &
            p_comm, p_status, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_RECV on ', my_process_mpi_all_id, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_recv_int_1d

  SUBROUTINE p_recv_int_2d (t_buffer, p_source, p_tag, p_count, comm)

    INTEGER, INTENT(out) :: t_buffer(:,:)
    INTEGER, INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_RECV (t_buffer, p_count, p_int, p_source, p_tag, &
            p_comm, p_status, p_error)
    ELSE
       CALL MPI_RECV (t_buffer, SIZE(t_buffer), p_int, p_source, p_tag, &
            p_comm, p_status, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_RECV on ', my_process_mpi_all_id, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_recv_int_2d

  SUBROUTINE p_recv_int_3d (t_buffer, p_source, p_tag, p_count, comm)

    INTEGER, INTENT(out) :: t_buffer(:,:,:)
    INTEGER, INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_RECV (t_buffer, p_count, p_int, p_source, p_tag, &
            p_comm, p_status, p_error)
    ELSE
       CALL MPI_RECV (t_buffer, SIZE(t_buffer), p_int, p_source, p_tag, &
            p_comm, p_status, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_RECV on ', my_process_mpi_all_id, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_recv_int_3d

  SUBROUTINE p_recv_int_4d (t_buffer, p_source, p_tag, p_count, comm)

    INTEGER, INTENT(out) :: t_buffer(:,:,:,:)
    INTEGER, INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_RECV (t_buffer, p_count, p_int, p_source, p_tag, &
            p_comm, p_status, p_error)
    ELSE
       CALL MPI_RECV (t_buffer, SIZE(t_buffer), p_int, p_source, p_tag, &
            p_comm, p_status, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_RECV on ', my_process_mpi_all_id, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_recv_int_4d


  SUBROUTINE p_recv_bool (t_buffer, p_source, p_tag, p_count, comm)

    LOGICAL, INTENT(out) :: t_buffer
    INTEGER, INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_RECV (t_buffer, p_count, p_bool, p_source, p_tag, &
            p_comm, p_status, p_error)
    ELSE
       CALL MPI_RECV (t_buffer, 1, p_bool, p_source, p_tag, &
            p_comm, p_status, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_RECV on ', my_process_mpi_all_id, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_recv_bool

  SUBROUTINE p_recv_bool_1d (t_buffer, p_source, p_tag, p_count, comm)

    LOGICAL, INTENT(out) :: t_buffer(:)
    INTEGER, INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_RECV (t_buffer, p_count, p_bool, p_source, p_tag, &
            p_comm, p_status, p_error)
    ELSE
       CALL MPI_RECV (t_buffer, SIZE(t_buffer), p_bool, p_source, p_tag, &
            p_comm, p_status, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_RECV on ', my_process_mpi_all_id, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_recv_bool_1d

  SUBROUTINE p_recv_bool_2d (t_buffer, p_source, p_tag, p_count, comm)

    LOGICAL, INTENT(out) :: t_buffer(:,:)
    INTEGER, INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_RECV (t_buffer, p_count, p_bool, p_source, p_tag, &
            p_comm, p_status, p_error)
    ELSE
       CALL MPI_RECV (t_buffer, SIZE(t_buffer), p_bool, p_source, p_tag, &
            p_comm, p_status, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_RECV on ', my_process_mpi_all_id, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_recv_bool_2d

  SUBROUTINE p_recv_bool_3d (t_buffer, p_source, p_tag, p_count, comm)

    LOGICAL, INTENT(out) :: t_buffer(:,:,:)
    INTEGER, INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_RECV (t_buffer, p_count, p_bool, p_source, p_tag, &
            p_comm, p_status, p_error)
    ELSE
       CALL MPI_RECV (t_buffer, SIZE(t_buffer), p_bool, p_source, p_tag, &
            p_comm, p_status, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_RECV on ', my_process_mpi_all_id, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_recv_bool_3d

  SUBROUTINE p_recv_bool_4d (t_buffer, p_source, p_tag, p_count, comm)

    LOGICAL, INTENT(out) :: t_buffer(:,:,:,:)
    INTEGER, INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_RECV (t_buffer, p_count, p_bool, p_source, p_tag, &
            p_comm, p_status, p_error)
    ELSE
       CALL MPI_RECV (t_buffer, SIZE(t_buffer), p_bool, p_source, p_tag, &
            p_comm, p_status, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_RECV on ', my_process_mpi_all_id, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_recv_bool_4d

  SUBROUTINE p_recv_char (t_buffer, p_source, p_tag, p_count, comm)

    CHARACTER(len=*),  INTENT(out) :: t_buffer
    INTEGER,           INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_RECV (t_buffer, p_count, p_char, p_source, p_tag, &
            p_comm, p_status, p_error)
    ELSE
       CALL MPI_RECV (t_buffer, LEN(t_buffer), p_char, p_source, p_tag, &
            p_comm, p_status, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) CALL finish (routine, 'MPI function call failed')
#endif
#endif

  END SUBROUTINE p_recv_char

  ! non-blocking receives

  !================================================================================================
  ! CHARACTER SECTION -----------------------------------------------------------------------------
  ! 
  SUBROUTINE p_irecv_char (t_buffer, p_source, p_tag, p_count, comm)

    CHARACTER(len=*), INTENT(inout) :: t_buffer
    INTEGER,   INTENT(in) :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_IRECV (t_buffer, p_count, p_char, p_source, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_IRECV (t_buffer, 1, p_char, p_source, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_IRECV on ', my_process_mpi_all_id, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_irecv_char
  !================================================================================================
  ! REAL SECTION ----------------------------------------------------------------------------------
  ! 
  SUBROUTINE p_irecv_real (t_buffer, p_source, p_tag, p_count, comm)

    REAL(dp),  INTENT(inout) :: t_buffer
    INTEGER,   INTENT(in) :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_IRECV (t_buffer, p_count, p_real_dp, p_source, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_IRECV (t_buffer, 1, p_real_dp, p_source, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_IRECV on ', my_process_mpi_all_id, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_irecv_real

  SUBROUTINE p_irecv_real_1d (t_buffer, p_source, p_tag, p_count, comm)

    REAL(dp),  INTENT(inout) :: t_buffer(:)
    INTEGER,   INTENT(in) :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_IRECV (t_buffer, p_count, p_real_dp, p_source, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_IRECV (t_buffer, SIZE(t_buffer), p_real_dp, p_source, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_IRECV on ', my_process_mpi_all_id, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_irecv_real_1d

  SUBROUTINE p_irecv_real_2d (t_buffer, p_source, p_tag, p_count, comm)

    REAL(dp),  INTENT(inout) :: t_buffer(:,:)
    INTEGER,   INTENT(in) :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_IRECV (t_buffer, p_count, p_real_dp, p_source, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_IRECV (t_buffer, SIZE(t_buffer), p_real_dp, p_source, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_IRECV on ', my_process_mpi_all_id, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_irecv_real_2d

  SUBROUTINE p_irecv_real_3d (t_buffer, p_source, p_tag, p_count, comm)

    REAL(dp),  INTENT(inout) :: t_buffer(:,:,:)
    INTEGER,   INTENT(in) :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_IRECV (t_buffer, p_count, p_real_dp, p_source, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_IRECV (t_buffer, SIZE(t_buffer), p_real_dp, p_source, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_IRECV on ', my_process_mpi_all_id, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_irecv_real_3d

  SUBROUTINE p_irecv_real_4d (t_buffer, p_source, p_tag, p_count, comm)

    REAL(dp),  INTENT(inout) :: t_buffer(:,:,:,:)
    INTEGER,   INTENT(in) :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_IRECV (t_buffer, p_count, p_real_dp, p_source, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_IRECV (t_buffer, SIZE(t_buffer), p_real_dp, p_source, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_IRECV on ', my_process_mpi_all_id, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_irecv_real_4d
  !================================================================================================
  ! INTEGER SECTION -------------------------------------------------------------------------------
  ! 
  SUBROUTINE p_irecv_int (t_buffer, p_source, p_tag, p_count, comm)

    INTEGER, INTENT(inout) :: t_buffer
    INTEGER,   INTENT(in) :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_IRECV (t_buffer, p_count, p_int, p_source, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_IRECV (t_buffer, 1, p_int, p_source, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_IRECV on ', my_process_mpi_all_id, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_irecv_int

  SUBROUTINE p_irecv_int_1d (t_buffer, p_source, p_tag, p_count, comm)

    INTEGER, INTENT(inout) :: t_buffer(:)
    INTEGER,   INTENT(in) :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_IRECV (t_buffer, p_count, p_int, p_source, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_IRECV (t_buffer, SIZE(t_buffer), p_int, p_source, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_IRECV on ', my_process_mpi_all_id, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_irecv_int_1d

  SUBROUTINE p_irecv_int_2d (t_buffer, p_source, p_tag, p_count, comm)

    INTEGER, INTENT(inout) :: t_buffer(:,:)
    INTEGER,   INTENT(in) :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_IRECV (t_buffer, p_count, p_int, p_source, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_IRECV (t_buffer, SIZE(t_buffer), p_int, p_source, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_IRECV on ', my_process_mpi_all_id, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_irecv_int_2d

  SUBROUTINE p_irecv_int_3d (t_buffer, p_source, p_tag, p_count, comm)

    INTEGER, INTENT(inout) :: t_buffer(:,:,:)
    INTEGER,   INTENT(in) :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_IRECV (t_buffer, p_count, p_int, p_source, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_IRECV (t_buffer, SIZE(t_buffer), p_int, p_source, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_IRECV on ', my_process_mpi_all_id, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_irecv_int_3d

  SUBROUTINE p_irecv_int_4d (t_buffer, p_source, p_tag, p_count, comm)

    INTEGER, INTENT(inout) :: t_buffer(:,:,:,:)
    INTEGER,   INTENT(in) :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_IRECV (t_buffer, p_count, p_int, p_source, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_IRECV (t_buffer, SIZE(t_buffer), p_int, p_source, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_IRECV on ', my_process_mpi_all_id, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_irecv_int_4d

  !================================================================================================
  ! LOGICAL SECTION -------------------------------------------------------------------------------
  ! 
  SUBROUTINE p_irecv_bool (t_buffer, p_source, p_tag, p_count, comm)

    LOGICAL, INTENT(inout) :: t_buffer
    INTEGER,   INTENT(in) :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_IRECV (t_buffer, p_count, p_bool, p_source, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_IRECV (t_buffer, 1, p_bool, p_source, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_IRECV on ', my_process_mpi_all_id, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_irecv_bool

  SUBROUTINE p_irecv_bool_1d (t_buffer, p_source, p_tag, p_count, comm)

    LOGICAL, INTENT(inout) :: t_buffer(:)
    INTEGER,   INTENT(in) :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_IRECV (t_buffer, p_count, p_bool, p_source, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_IRECV (t_buffer, SIZE(t_buffer), p_bool, p_source, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_IRECV on ', my_process_mpi_all_id, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_irecv_bool_1d

  SUBROUTINE p_irecv_bool_2d (t_buffer, p_source, p_tag, p_count, comm)

    LOGICAL, INTENT(inout) :: t_buffer(:,:)
    INTEGER,   INTENT(in) :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_IRECV (t_buffer, p_count, p_bool, p_source, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_IRECV (t_buffer, SIZE(t_buffer), p_bool, p_source, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_IRECV on ', my_process_mpi_all_id, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_irecv_bool_2d

  SUBROUTINE p_irecv_bool_3d (t_buffer, p_source, p_tag, p_count, comm)

    LOGICAL, INTENT(inout) :: t_buffer(:,:,:)
    INTEGER,   INTENT(in) :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_IRECV (t_buffer, p_count, p_bool, p_source, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_IRECV (t_buffer, SIZE(t_buffer), p_bool, p_source, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_IRECV on ', my_process_mpi_all_id, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_irecv_bool_3d

  SUBROUTINE p_irecv_bool_4d (t_buffer, p_source, p_tag, p_count, comm)

    LOGICAL, INTENT(inout) :: t_buffer(:,:,:,:)
    INTEGER,   INTENT(in) :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_IRECV (t_buffer, p_count, p_bool, p_source, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_IRECV (t_buffer, SIZE(t_buffer), p_bool, p_source, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_IRECV on ', my_process_mpi_all_id, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_irecv_bool_4d

  SUBROUTINE p_pack_int (t_var, t_buffer, p_buf_size, p_pos, comm)

    INTEGER,   INTENT(IN)    :: t_var
    CHARACTER, INTENT(INOUT) :: t_buffer(:)
    INTEGER,   INTENT(IN)    :: p_buf_size
    INTEGER,   INTENT(INOUT) :: p_pos
    INTEGER, OPTIONAL, INTENT(IN) :: comm
#ifndef NOMPI
    INTEGER :: p_comm, outsize

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    outsize = p_buf_size
    CALL MPI_PACK(t_var, 1, p_int, t_buffer, outsize, p_pos, p_comm, p_error)
#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) CALL finish ("p_pack_int", 'MPI call failed')
#endif
#endif
  END SUBROUTINE p_pack_int

  SUBROUTINE p_pack_real (t_var, t_buffer, p_buf_size, p_pos, comm)

    REAL(wp),  INTENT(IN)    :: t_var
    CHARACTER, INTENT(INOUT) :: t_buffer(:)
    INTEGER,   INTENT(IN)    :: p_buf_size
    INTEGER,   INTENT(INOUT) :: p_pos
    INTEGER, OPTIONAL, INTENT(IN) :: comm
#ifndef NOMPI
    INTEGER :: p_comm, outsize

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    outsize = p_buf_size
    CALL MPI_PACK(t_var, 1, p_real_dp, t_buffer, outsize, p_pos, p_comm, p_error)
#ifdef DEBUG
   IF (p_error /= MPI_SUCCESS) CALL finish ("p_pack_real", 'MPI call failed')
#endif
#endif
  END SUBROUTINE p_pack_real

  SUBROUTINE p_pack_int_1d (t_var, p_count, t_buffer, p_buf_size, p_pos, comm)

    INTEGER,   INTENT(IN)    :: t_var(:)
    INTEGER,   INTENT(IN)    :: p_count
    CHARACTER, INTENT(INOUT) :: t_buffer(:)
    INTEGER,   INTENT(IN)    :: p_buf_size
    INTEGER,   INTENT(INOUT) :: p_pos
    INTEGER, OPTIONAL, INTENT(IN) :: comm
#ifndef NOMPI
    INTEGER :: p_comm, outsize

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    outsize = p_buf_size
    CALL MPI_PACK(t_var, p_count, p_int, t_buffer, outsize, p_pos, p_comm, p_error)
#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) CALL finish ("p_pack_int_1d", 'MPI call failed')
#endif
#endif
  END SUBROUTINE p_pack_int_1d

  SUBROUTINE p_pack_real_1d (t_var, p_count, t_buffer, p_buf_size, p_pos, comm)

    REAL(wp),  INTENT(IN)    :: t_var(:)
    INTEGER,   INTENT(IN)    :: p_count
    CHARACTER, INTENT(INOUT) :: t_buffer(:)
    INTEGER,   INTENT(IN)    :: p_buf_size
    INTEGER,   INTENT(INOUT) :: p_pos
    INTEGER, OPTIONAL, INTENT(IN) :: comm
#ifndef NOMPI
    INTEGER :: p_comm, outsize

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    outsize = p_buf_size
    CALL MPI_PACK(t_var, p_count, p_real_dp, t_buffer, outsize, p_pos, p_comm, p_error)
#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) CALL finish ("p_pack_real_1d", 'MPI call failed')
#endif
#endif
  END SUBROUTINE p_pack_real_1d
  

  SUBROUTINE p_pack_string (t_var, t_buffer, p_buf_size, p_pos, comm)

    CHARACTER(LEN=*), INTENT(IN) :: t_var
    CHARACTER, INTENT(INOUT) :: t_buffer(:)
    INTEGER,   INTENT(IN)    :: p_buf_size
    INTEGER,   INTENT(INOUT) :: p_pos
    INTEGER, OPTIONAL, INTENT(IN) :: comm
#ifndef NOMPI
    INTEGER :: p_comm, outsize, ilength
    CHARACTER(LEN=128) :: tmp_var

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    outsize = p_buf_size
    tmp_var =     TRIM(t_var)
    ilength = LEN_TRIM(t_var)
    ! first pack string length, then the character sequence
    CALL MPI_PACK(ilength,       1,  p_int, t_buffer, outsize, p_pos, p_comm, p_error)
#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) CALL finish ("p_pack_char_1d", 'MPI call failed')
#endif
    CALL MPI_PACK(tmp_var, ilength, p_char, t_buffer, outsize, p_pos, p_comm, p_error)
#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) CALL finish ("p_pack_char_1d", 'MPI call failed')
#endif
#endif
  END SUBROUTINE p_pack_string

  SUBROUTINE p_pack_real_2d (t_var, p_count, t_buffer, p_buf_size, p_pos, comm)

    REAL(wp),  INTENT(IN)    :: t_var(:,:)
    INTEGER,   INTENT(IN)    :: p_count
    CHARACTER, INTENT(INOUT) :: t_buffer(:)
    INTEGER,   INTENT(IN)    :: p_buf_size
    INTEGER,   INTENT(INOUT) :: p_pos
    INTEGER, OPTIONAL, INTENT(IN) :: comm
#ifndef NOMPI
    INTEGER :: p_comm, outsize

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    outsize = p_buf_size
    CALL MPI_PACK(t_var, p_count, p_real_dp, t_buffer, outsize, p_pos, p_comm, p_error)
#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) CALL finish ("p_pack_real_2d", 'MPI call failed')
#endif
#endif
  END SUBROUTINE p_pack_real_2d

  SUBROUTINE p_unpack_int (t_buffer, p_buf_size, p_pos, t_var, comm)

    CHARACTER, INTENT(IN)    :: t_buffer(:)
    INTEGER,   INTENT(IN)    :: p_buf_size
    INTEGER,   INTENT(INOUT) :: p_pos
    INTEGER,   INTENT(OUT)   :: t_var
    INTEGER, OPTIONAL, INTENT(IN) :: comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    CALL MPI_UNPACK(t_buffer, p_buf_size, p_pos, t_var, 1, p_int, p_comm, p_error)
#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) CALL finish ("p_unpack_int", 'MPI call failed')
#endif
#endif
  END SUBROUTINE p_unpack_int

  SUBROUTINE p_unpack_real (t_buffer, p_buf_size, p_pos, t_var, comm)

    CHARACTER, INTENT(IN)    :: t_buffer(:)
    INTEGER,   INTENT(IN)    :: p_buf_size
    INTEGER,   INTENT(INOUT) :: p_pos
    REAL(wp),  INTENT(OUT)   :: t_var
    INTEGER, OPTIONAL, INTENT(IN) :: comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    CALL MPI_UNPACK(t_buffer, p_buf_size, p_pos, t_var, 1, p_real_dp, p_comm, p_error)
#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) CALL finish ("p_unpack_real", 'MPI call failed')
#endif
#endif
  END SUBROUTINE p_unpack_real

  SUBROUTINE p_unpack_int_1d (t_buffer, p_buf_size, p_pos, t_var, p_count, comm)

    CHARACTER, INTENT(IN)    :: t_buffer(:)
    INTEGER,   INTENT(IN)    :: p_buf_size
    INTEGER,   INTENT(INOUT) :: p_pos
    INTEGER,   INTENT(INOUT) :: t_var(:)
    INTEGER,   INTENT(IN)    :: p_count
    INTEGER, OPTIONAL, INTENT(IN) :: comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    CALL MPI_UNPACK(t_buffer, p_buf_size, p_pos, t_var, p_count, p_int, p_comm, p_error)
#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) CALL finish ("p_unpack_int_1d", 'MPI call failed')
#endif
#endif
  END SUBROUTINE p_unpack_int_1d

  SUBROUTINE p_unpack_real_1d (t_buffer, p_buf_size, p_pos, t_var, p_count, comm)

    CHARACTER, INTENT(IN)    :: t_buffer(:)
    INTEGER,   INTENT(IN)    :: p_buf_size
    INTEGER,   INTENT(INOUT) :: p_pos
    REAL(wp),  INTENT(INOUT) :: t_var(:)
    INTEGER,   INTENT(IN)    :: p_count
    INTEGER, OPTIONAL, INTENT(IN) :: comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    CALL MPI_UNPACK(t_buffer, p_buf_size, p_pos, t_var, p_count, p_real_dp, p_comm, p_error)
#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) CALL finish ("p_unpack_real_1d", 'MPI call failed')
#endif
#endif
  END SUBROUTINE p_unpack_real_1d
  

  SUBROUTINE p_unpack_string (t_buffer, p_buf_size, p_pos, t_var, comm)

    CHARACTER, INTENT(IN) :: t_buffer(:)
    INTEGER,   INTENT(IN)    :: p_buf_size
    INTEGER,   INTENT(INOUT) :: p_pos
    CHARACTER(LEN=*), INTENT(INOUT) :: t_var
    INTEGER, OPTIONAL, INTENT(IN) :: comm
#ifndef NOMPI
    INTEGER :: p_comm, ilength

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF
    ! first unpack string length, then the character sequence
    CALL MPI_UNPACK(t_buffer, p_buf_size, p_pos, ilength,       1,  p_int, p_comm, p_error)
#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) CALL finish ("p_unpack_char_1d", 'MPI call failed')
#endif
    CALL MPI_UNPACK(t_buffer, p_buf_size, p_pos,   t_var, ilength, p_char, p_comm, p_error)
#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) CALL finish ("p_unpack_char_1d", 'MPI call failed')
#endif
#endif
  END SUBROUTINE p_unpack_string


  SUBROUTINE p_unpack_real_2d (t_buffer, p_buf_size, p_pos, t_var, p_count, comm)

    CHARACTER, INTENT(IN)    :: t_buffer(:)
    INTEGER,   INTENT(IN)    :: p_buf_size
    INTEGER,   INTENT(INOUT) :: p_pos
    REAL(wp),  INTENT(INOUT) :: t_var(:,:)
    INTEGER,   INTENT(IN)    :: p_count
    INTEGER, OPTIONAL, INTENT(IN) :: comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    CALL MPI_UNPACK(t_buffer, p_buf_size, p_pos, t_var, p_count, p_real_dp, p_comm, p_error)
#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) CALL finish ("p_unpack_real_2d", 'MPI call failed')
#endif
#endif
  END SUBROUTINE p_unpack_real_2d

  SUBROUTINE p_irecv_packed (t_buffer, p_source, p_tag, p_count, comm)

    CHARACTER, INTENT(INOUT) :: t_buffer(:)
    INTEGER,   INTENT(IN)    :: p_source, p_tag
    INTEGER,   INTENT(IN)    :: p_count
    INTEGER, OPTIONAL, INTENT(IN) :: comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    CALL p_inc_request
    CALL MPI_IRECV (t_buffer, p_count, MPI_PACKED, p_source, p_tag, &
      p_comm, p_request(p_irequest), p_error)

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_IRECV on ', my_process_mpi_all_id, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_irecv_packed

  SUBROUTINE p_send_packed (t_buffer, p_destination, p_tag, p_count, comm)

    CHARACTER, INTENT(in) :: t_buffer(:)
    INTEGER,   INTENT(in) :: p_destination, p_tag
    INTEGER,   INTENT(in) :: p_count
    INTEGER, OPTIONAL, INTENT(in) :: comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    CALL MPI_SEND (t_buffer, p_count, MPI_PACKED, p_destination, p_tag, &
      p_comm, p_error)

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', my_process_mpi_all_id, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_send_packed

  !
  !================================================================================================

  ! sendrecv implementation

  SUBROUTINE p_sendrecv_real_1d (sendbuf, p_dest, recvbuf, p_source, &
                                  p_tag, comm)

    REAL(dp), INTENT(in)           :: sendbuf (:)
    INTEGER,  INTENT(in)           :: p_dest
    REAL(dp), INTENT(out)          :: recvbuf (:)
    INTEGER,  INTENT(in)           :: p_source
    INTEGER,  INTENT(in)           :: p_tag
    INTEGER,  INTENT(in) ,OPTIONAL :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    CALL MPI_SENDRECV (sendbuf, SIZE(sendbuf), p_real_dp, p_dest,   p_tag, &
                        recvbuf, SIZE(recvbuf), p_real_dp, p_source, p_tag, &
                        p_comm, p_status, p_error)

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i4,a,i6,a)') ' MPI_SENDRECV by ', my_process_mpi_all_id, &
            ' to ', p_dest, ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_sendrecv_real_1d

  SUBROUTINE p_sendrecv_real_2d (sendbuf, p_dest, recvbuf, p_source, &
                                  p_tag, comm)

    REAL(dp), INTENT(in)           :: sendbuf (:,:)
    INTEGER,  INTENT(in)           :: p_dest
    REAL(dp), INTENT(out)          :: recvbuf (:,:)
    INTEGER,  INTENT(in)           :: p_source
    INTEGER,  INTENT(in)           :: p_tag
    INTEGER,  INTENT(in) ,OPTIONAL :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    CALL MPI_SENDRECV (sendbuf, SIZE(sendbuf), p_real_dp, p_dest,   p_tag, &
                        recvbuf, SIZE(recvbuf), p_real_dp, p_source, p_tag, &
                        p_comm, p_status, p_error)

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i4,a,i6,a)') ' MPI_SENDRECV by ', my_process_mpi_all_id, &
            ' to ', p_dest, ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_sendrecv_real_2d

  SUBROUTINE p_sendrecv_real_3d (sendbuf, p_dest, recvbuf, p_source, &
                                  p_tag, comm)

    REAL(dp), INTENT(in)           :: sendbuf (:,:,:)
    INTEGER,  INTENT(in)           :: p_dest
    REAL(dp), INTENT(out)          :: recvbuf (:,:,:)
    INTEGER,  INTENT(in)           :: p_source
    INTEGER,  INTENT(in)           :: p_tag
    INTEGER,  INTENT(in) ,OPTIONAL :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    CALL MPI_SENDRECV (sendbuf, SIZE(sendbuf), p_real_dp, p_dest,   p_tag, &
                        recvbuf, SIZE(recvbuf), p_real_dp, p_source, p_tag, &
                        p_comm, p_status, p_error)

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i4,a,i6,a)') ' MPI_SENDRECV by ', my_process_mpi_all_id, &
            ' to ', p_dest, ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_sendrecv_real_3d

  SUBROUTINE p_sendrecv_real_4d (sendbuf, p_dest, recvbuf, p_source, &
                                  p_tag, comm)

    REAL(dp), INTENT(in)           :: sendbuf (:,:,:,:)
    INTEGER,  INTENT(in)           :: p_dest
    REAL(dp), INTENT(out)          :: recvbuf (:,:,:,:)
    INTEGER,  INTENT(in)           :: p_source
    INTEGER,  INTENT(in)           :: p_tag
    INTEGER,  INTENT(in) ,OPTIONAL :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    CALL MPI_SENDRECV (sendbuf, SIZE(sendbuf), p_real_dp, p_dest,   p_tag, &
                        recvbuf, SIZE(recvbuf), p_real_dp, p_source, p_tag, &
                        p_comm, p_status, p_error)

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i4,a,i6,a)') ' MPI_SENDRECV by ', my_process_mpi_all_id, &
            ' to ', p_dest, ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_sendrecv_real_4d

  ! bcast implementation

  SUBROUTINE p_bcast_real (t_buffer, p_source, comm)

    REAL (dp), INTENT(inout) :: t_buffer
    INTEGER,   INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

#ifdef DEBUG
    nbcast = nbcast+1
#endif

    IF (process_mpi_all_size == 1) THEN
       RETURN
    ELSE
       CALL MPI_BCAST (t_buffer, 1, p_real_dp, p_source, &
            p_comm, p_error)
    ENDIF

#ifdef DEBUG
    WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' with broadcast number ', nbcast, ' succesfull.'

    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_bcast_real

  SUBROUTINE p_bcast_real_1d (t_buffer, p_source, comm)

    REAL (dp), INTENT(inout) :: t_buffer(:)
    INTEGER,   INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

#ifdef DEBUG
    nbcast = nbcast+1
#endif

    IF (process_mpi_all_size == 1) THEN
       RETURN
    ELSE
       CALL MPI_BCAST (t_buffer, SIZE(t_buffer), p_real_dp, p_source, &
            p_comm, p_error)
    ENDIF

#ifdef DEBUG
    WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' with broadcast number ', nbcast, ' succesfull.'

     IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_bcast_real_1d

  SUBROUTINE p_bcast_real_2d (t_buffer, p_source, comm)

    REAL (dp), INTENT(inout) :: t_buffer(:,:)
    INTEGER,   INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

#ifdef DEBUG
    nbcast = nbcast+1
#endif

    IF (process_mpi_all_size == 1) THEN
       RETURN
    ELSE
       CALL MPI_BCAST (t_buffer, SIZE(t_buffer), p_real_dp, p_source, &
            p_comm, p_error)
    ENDIF

#ifdef DEBUG
    WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' with broadcast number ', nbcast, ' succesfull.'

     IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_bcast_real_2d


  SUBROUTINE p_bcast_real_2d_single (t_buffer, p_source, comm)

    REAL (sp), INTENT(inout) :: t_buffer(:,:) ! SINGLE PRECISION
    INTEGER,   INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

#ifdef DEBUG
    nbcast = nbcast+1
#endif

    IF (process_mpi_all_size == 1) THEN
       RETURN
    ELSE
       CALL MPI_BCAST (t_buffer, SIZE(t_buffer), p_real_sp, p_source, &
            p_comm, p_error)
    ENDIF

#ifdef DEBUG
    WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' with broadcast number ', nbcast, ' succesfull.'

     IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_bcast_real_2d_single


  SUBROUTINE p_bcast_real_3d (t_buffer, p_source, comm)

    REAL (dp), INTENT(inout) :: t_buffer(:,:,:)
    INTEGER,   INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

#ifdef DEBUG
    nbcast = nbcast+1
#endif

    IF (process_mpi_all_size == 1) THEN
       RETURN
    ELSE
       CALL MPI_BCAST (t_buffer, SIZE(t_buffer), p_real_dp, p_source, &
            p_comm, p_error)
    ENDIF

#ifdef DEBUG
    WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' with broadcast number ', nbcast, ' succesfull.'

    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_bcast_real_3d

  SUBROUTINE p_bcast_real_4d (t_buffer, p_source, comm)

    REAL (dp), INTENT(inout) :: t_buffer(:,:,:,:)
    INTEGER,   INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

#ifdef DEBUG
    nbcast = nbcast+1
#endif

    IF (process_mpi_all_size == 1) THEN
       RETURN
    ELSE
       CALL MPI_BCAST (t_buffer, SIZE(t_buffer), p_real_dp, p_source, &
            p_comm, p_error)
    ENDIF

#ifdef DEBUG
    WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' with broadcast number ', nbcast, ' succesfull.'

     IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_bcast_real_4d

  SUBROUTINE p_bcast_real_5d (t_buffer, p_source, comm)

    REAL (dp), INTENT(inout) :: t_buffer(:,:,:,:,:)
    INTEGER,   INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

#ifdef DEBUG
    nbcast = nbcast+1
#endif

    IF (process_mpi_all_size == 1) THEN
       RETURN
    ELSE
       CALL MPI_BCAST (t_buffer, SIZE(t_buffer), p_real_dp, p_source, &
            p_comm, p_error)
    ENDIF

#ifdef DEBUG
    WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' with broadcast number ', nbcast, ' succesfull.'

     IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_bcast_real_5d

  SUBROUTINE p_bcast_int_i4 (t_buffer, p_source, comm)

    INTEGER (i4), INTENT(inout) :: t_buffer
    INTEGER, INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

#ifdef DEBUG
    nbcast = nbcast+1
#endif

    IF (process_mpi_all_size == 1) THEN
       RETURN
    ELSE
       CALL MPI_BCAST (t_buffer, 1, p_int_i4, p_source, &
            p_comm, p_error)
    ENDIF

#ifdef DEBUG
    WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' with broadcast number ', nbcast, ' succesfull.'

    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_bcast_int_i4

  SUBROUTINE p_bcast_int_i8 (t_buffer, p_source, comm)

    INTEGER (i8), INTENT(inout) :: t_buffer
    INTEGER, INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

#ifdef DEBUG
    nbcast = nbcast+1
#endif

    IF (process_mpi_all_size == 1) THEN
       RETURN
    ELSE
       CALL MPI_BCAST (t_buffer, 1, p_int_i8, p_source, &
            p_comm, p_error)
    ENDIF

#ifdef DEBUG
    WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' with broadcast number ', nbcast, ' succesfull.'

     IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_bcast_int_i8

  SUBROUTINE p_bcast_int_1d (t_buffer, p_source, comm)

    INTEGER, INTENT(inout) :: t_buffer(:)
    INTEGER, INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

#ifdef DEBUG
    nbcast = nbcast+1
#endif

    IF (process_mpi_all_size == 1) THEN
       RETURN
    ELSE
       CALL MPI_BCAST (t_buffer, SIZE(t_buffer), p_int, p_source, &
            p_comm, p_error)
    ENDIF

#ifdef DEBUG
    WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' with broadcast number ', nbcast, ' succesfull.'

     IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_bcast_int_1d

  SUBROUTINE p_bcast_int_i8_1d (t_buffer, p_source, comm)

    INTEGER(i8), INTENT(inout) :: t_buffer(:)
    INTEGER, INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

#ifdef DEBUG
    nbcast = nbcast+1
#endif

    IF (process_mpi_all_size == 1) THEN
       RETURN
    ELSE
       CALL MPI_BCAST (t_buffer, SIZE(t_buffer), p_int_i8, p_source, &
            p_comm, p_error)
    ENDIF

#ifdef DEBUG
    WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' with broadcast number ', nbcast, ' succesfull.'

     IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_bcast_int_i8_1d

  SUBROUTINE p_bcast_int_2d (t_buffer, p_source, comm)

    INTEGER, INTENT(inout) :: t_buffer(:,:)
    INTEGER, INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

#ifdef DEBUG
    nbcast = nbcast+1
#endif

    IF (process_mpi_all_size == 1) THEN
       RETURN
    ELSE
       CALL MPI_BCAST (t_buffer, SIZE(t_buffer), p_int, p_source, &
            p_comm, p_error)
    ENDIF

#ifdef DEBUG
    WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' with broadcast number ', nbcast, ' succesfull.'

     IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_bcast_int_2d

  SUBROUTINE p_bcast_int_3d (t_buffer, p_source, comm)

    INTEGER, INTENT(inout) :: t_buffer(:,:,:)
    INTEGER, INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

#ifdef DEBUG
    nbcast = nbcast+1
#endif

    IF (process_mpi_all_size == 1) THEN
       RETURN
    ELSE
       CALL MPI_BCAST (t_buffer, SIZE(t_buffer), p_int, p_source, &
            p_comm, p_error)
    ENDIF

#ifdef DEBUG
    WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' with broadcast number ', nbcast, ' succesfull.'

     IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_bcast_int_3d

  SUBROUTINE p_bcast_int_4d (t_buffer, p_source, comm)

    INTEGER, INTENT(inout) :: t_buffer(:,:,:,:)
    INTEGER, INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

#ifdef DEBUG
    nbcast = nbcast+1
#endif

    IF (process_mpi_all_size == 1) THEN
       RETURN
    ELSE
       CALL MPI_BCAST (t_buffer, SIZE(t_buffer), p_int, p_source, &
            p_comm, p_error)
    ENDIF

#ifdef DEBUG
    WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' with broadcast number ', nbcast, ' succesfull.'

     IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_bcast_int_4d


  SUBROUTINE p_bcast_bool (t_buffer, p_source, comm)

    LOGICAL, INTENT(inout) :: t_buffer
    INTEGER, INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

#ifdef DEBUG
    nbcast = nbcast+1
#endif

    IF (process_mpi_all_size == 1) THEN
       RETURN
    ELSE
       CALL MPI_BCAST (t_buffer, 1, p_bool, p_source, &
            p_comm, p_error)
    ENDIF

#ifdef DEBUG
    WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' with broadcast number ', nbcast, ' succesfull.'

     IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_bcast_bool

  SUBROUTINE p_bcast_bool_1d (t_buffer, p_source, comm)

    LOGICAL, INTENT(inout) :: t_buffer(:)
    INTEGER, INTENT(in) :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

#ifdef DEBUG
    nbcast = nbcast+1
#endif

    IF (process_mpi_all_size == 1) THEN
       RETURN
    ELSE
       CALL MPI_BCAST (t_buffer, SIZE(t_buffer), p_bool, p_source, &
            p_comm, p_error)
    ENDIF

#ifdef DEBUG
    WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' with broadcast number ', nbcast, ' succesfull.'

    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_bcast_bool_1d

  SUBROUTINE p_bcast_bool_2d (t_buffer, p_source, comm)

    LOGICAL, INTENT(inout) :: t_buffer(:,:)
    INTEGER, INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

#ifdef DEBUG
    nbcast = nbcast+1
#endif

    IF (process_mpi_all_size == 1) THEN
       RETURN
    ELSE
       CALL MPI_BCAST (t_buffer, SIZE(t_buffer), p_bool, p_source, &
            p_comm, p_error)
    ENDIF

#ifdef DEBUG
    WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' with broadcast number ', nbcast, ' succesfull.'

     IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_bcast_bool_2d

  SUBROUTINE p_bcast_bool_3d (t_buffer, p_source, comm)

    LOGICAL, INTENT(inout) :: t_buffer(:,:,:)
    INTEGER, INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

#ifdef DEBUG
    nbcast = nbcast+1
#endif

    IF (process_mpi_all_size == 1) THEN
       RETURN
    ELSE
       CALL MPI_BCAST (t_buffer, SIZE(t_buffer), p_bool, p_source, &
            p_comm, p_error)
    ENDIF

#ifdef DEBUG
    WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' with broadcast number ', nbcast, ' succesfull.'

    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_bcast_bool_3d

  SUBROUTINE p_bcast_bool_4d (t_buffer, p_source, comm)

    LOGICAL, INTENT(inout) :: t_buffer(:,:,:,:)
    INTEGER, INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

#ifdef DEBUG
    nbcast = nbcast+1
#endif

    IF (process_mpi_all_size == 1) THEN
       RETURN
    ELSE
       CALL MPI_BCAST (t_buffer, SIZE(t_buffer), p_bool, p_source, &
            p_comm, p_error)
    ENDIF

#ifdef DEBUG
    WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' with broadcast number ', nbcast, ' succesfull.'

     IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_bcast_bool_4d

  SUBROUTINE p_bcast_char (t_buffer, p_source, comm)

    CHARACTER(len=*),  INTENT(inout) :: t_buffer
    INTEGER,           INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

#ifdef DEBUG
    nbcast = nbcast+1
#endif

    IF (process_mpi_all_size == 1) THEN
       RETURN
    ELSE
       CALL MPI_BCAST (t_buffer, LEN(t_buffer), p_char, p_source, &
            p_comm, p_error)
    ENDIF

#ifdef DEBUG
    WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' with broadcast number ', nbcast, ' succesfull.'

     IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_bcast_char

  SUBROUTINE p_bcast_char_1d (t_buffer, p_source, comm)

    CHARACTER (*), INTENT(inout) :: t_buffer(:)
    INTEGER,       INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm
    INTEGER                      :: lexlength, flength

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

#ifdef DEBUG
    nbcast = nbcast+1
#endif

    IF (process_mpi_all_size == 1) THEN
       RETURN
    ELSE
       lexlength=LEN(t_buffer(1))
       flength=SIZE(t_buffer)
       lexlength=lexlength*flength
       CALL MPI_BCAST (t_buffer, lexlength, p_char, p_source, p_comm, p_error)
#ifdef DEBUG
       WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' with broadcast number ', nbcast, ' succesfull.'

       IF (p_error /= MPI_SUCCESS) THEN
          WRITE (nerr,'(a,i4,a)') ' MPI_BCAST from ', p_source, &
               ' failed.'
          WRITE (nerr,'(a,i4)') ' Error = ', p_error
          CALL p_abort
       END IF
#endif
    ENDIF
#endif

  END SUBROUTINE p_bcast_char_1d



  SUBROUTINE p_bcast_cchar (t_buffer, buflen, p_source, p_comm)

#ifdef __SX__
    USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_CHAR
#else
    USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_SIGNED_CHAR
#endif
    INTEGER,           INTENT(IN)    :: buflen
#ifdef __SX__
    CHARACTER(C_CHAR), INTENT(INOUT) :: t_buffer(buflen)
#else
    INTEGER(C_SIGNED_CHAR), INTENT(INOUT) :: t_buffer(buflen)
#endif

    INTEGER,           INTENT(IN)    :: p_source
    INTEGER,           INTENT(IN)    :: p_comm

#ifndef NOMPI

#ifdef DEBUG
    nbcast = nbcast+1
#endif

    IF (process_mpi_all_size == 1) THEN
       RETURN
    ELSE
       CALL MPI_BCAST (t_buffer, buflen, p_char, p_source, p_comm, p_error)
#ifdef DEBUG
       WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' with broadcast number ', nbcast, ' succesfull.'

       IF (p_error /= MPI_SUCCESS) THEN
          WRITE (nerr,'(a,i4,a)') ' MPI_BCAST from ', p_source, &
               ' failed.'
          WRITE (nerr,'(a,i4)') ' Error = ', p_error
          CALL p_abort
       END IF
#endif
    ENDIF
#endif

  END SUBROUTINE p_bcast_cchar


  ! probe implementation

  SUBROUTINE p_probe (p_tagcount, p_tagtable, p_source, &
       p_tag, p_count, comm)

    INTEGER,   INTENT(in)  :: p_tagcount, p_tagtable(:)
    INTEGER,   INTENT(out) :: p_source, p_tag, p_count
    INTEGER, OPTIONAL, INTENT(in) :: comm

#ifndef NOMPI
    INTEGER :: p_comm, i
    LOGICAL :: flag

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    p_tag = -1
    DO WHILE (p_tag == -1)
       DO i = 1, p_tagcount
          CALL MPI_IPROBE (MPI_ANY_SOURCE, p_tagtable(i), p_comm, &
               flag, p_status, p_error)
#ifdef DEBUG
          IF (p_error /= MPI_SUCCESS) THEN
             WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_IPROBE on ', my_process_mpi_all_id, &
                  ' for tag ', p_tagtable(i), ' failed.'
             WRITE (nerr,'(a,i4)') ' Error = ', p_error
             CALL p_abort
          END IF
#endif
          IF (flag) THEN
             p_source = p_status(MPI_SOURCE)
             p_tag = p_status(MPI_TAG)
             CALL MPI_GET_COUNT(p_status, p_real_dp, p_count, p_error)
#ifdef DEBUG
             IF (p_error /= MPI_SUCCESS) THEN
                WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_GET_COUNT on ', &
                  & my_process_mpi_all_id, &
                     ' for tag ', p_tag, ' from ' , p_source, ' failed.'
                WRITE (nerr,'(a,i4)') ' Error = ', p_error
                CALL p_abort
             END IF
#endif
             EXIT
          ELSE
             p_tag = -1
          END IF
       END DO
    END DO

#endif

  END SUBROUTINE p_probe
  !------------------------------------------------------

  !------------------------------------------------------
  SUBROUTINE p_wait(request)
    INTEGER, INTENT(INOUT), OPTIONAL :: request
#ifndef NOMPI
    INTEGER :: p_status_wait(MPI_STATUS_SIZE,p_irequest)

    IF (PRESENT(request)) THEN
      CALL MPI_WAIT(request,p_status_wait(:,1), p_error)
      CALL p_clear_request(request)
    ELSE
      CALL MPI_WAITALL(p_irequest, p_request, p_status_wait, p_error)
      p_irequest = 0
    END IF
#endif
  END SUBROUTINE p_wait

  SUBROUTINE p_wait_any(return_pe)

    INTEGER, INTENT(out) :: return_pe
#ifndef NOMPI
    INTEGER :: i

    CALL MPI_WAITANY(p_irequest, p_request, i, p_status, p_error)
    IF (i == MPI_UNDEFINED) THEN
      p_irequest = 0
      return_pe = -1
    ELSE
      return_pe = p_status(MPI_SOURCE)
    ENDIF
#else
    return_pe = 0
#endif

  END SUBROUTINE p_wait_any
  !------------------------------------------------------

  !------------------------------------------------------
  SUBROUTINE p_barrier (comm)
  INTEGER ,INTENT(IN) ,OPTIONAL :: comm
#ifndef NOMPI
    INTEGER :: com
    com = MPI_COMM_WORLD; IF(PRESENT(comm)) com = comm
    CALL MPI_BARRIER (com, p_error)

!#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' MPI_BARRIER on ', my_process_mpi_all_id, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
!#endif
#endif

  END SUBROUTINE p_barrier
  !------------------------------------------------------

  !------------------------------------------------------
  SUBROUTINE global_mpi_barrier 
#ifndef NOMPI
    CALL MPI_BARRIER (global_mpi_communicator, p_error)

!#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' global_mpi_barrier on ', get_my_global_mpi_id(), ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
!#endif
#endif

  END SUBROUTINE global_mpi_barrier
  !------------------------------------------------------

  !------------------------------------------------------
  SUBROUTINE work_mpi_barrier 
#ifndef NOMPI
    CALL MPI_BARRIER (p_comm_work, p_error)

!#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' global_mpi_barrier on ', get_my_global_mpi_id(), ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
!#endif
#endif

  END SUBROUTINE work_mpi_barrier
  !------------------------------------------------------

  FUNCTION p_sum_dp_0d (zfield, comm) RESULT (p_sum)

    REAL(dp)                      :: p_sum
    REAL(dp),          INTENT(in) :: zfield
    INTEGER, OPTIONAL, INTENT(in) :: comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    IF (my_process_is_mpi_parallel()) THEN
       CALL MPI_ALLREDUCE (zfield, p_sum, 1, p_real_dp, &
            MPI_SUM, p_comm, p_error)
    ELSE
       p_sum = zfield
    END IF
#else
    p_sum = zfield
#endif

  END FUNCTION p_sum_dp_0d

  FUNCTION p_sum_dp_1d (zfield, comm) RESULT (p_sum)

    REAL(dp),          INTENT(in) :: zfield(:)
    INTEGER, OPTIONAL, INTENT(in) :: comm
    REAL(dp)                      :: p_sum (SIZE(zfield))

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    IF (my_process_is_mpi_all_parallel()) THEN
       CALL MPI_ALLREDUCE (zfield, p_sum, SIZE(zfield), p_real_dp, &
            MPI_SUM, p_comm, p_error)
    ELSE
       p_sum = zfield
    END IF
#else
    p_sum = zfield
#endif

  END FUNCTION p_sum_dp_1d

  FUNCTION p_sum_i8_1d (kfield, comm) RESULT (p_sum)

    INTEGER(i8),       INTENT(in) :: kfield(:)
    INTEGER, OPTIONAL, INTENT(in) :: comm
    INTEGER(i8)                   :: p_sum (SIZE(kfield))

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    IF (my_process_is_mpi_all_parallel()) THEN
       CALL MPI_ALLREDUCE (kfield, p_sum, SIZE(kfield), p_int_i8, &
            MPI_SUM, p_comm, p_error)
    ELSE
       p_sum = kfield
    END IF
#else
    p_sum = kfield
#endif

  END FUNCTION p_sum_i8_1d

  FUNCTION p_global_sum_1d (zfield, comm) RESULT (p_sum)

    REAL(dp),          INTENT(in) :: zfield(:)
    INTEGER, OPTIONAL, INTENT(in) :: comm
    REAL(dp)                      :: p_sum
    REAL(dp)                      :: pe_sums(SIZE(zfield))

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    IF (my_process_is_mpi_all_parallel()) THEN
       CALL MPI_REDUCE (zfield, pe_sums, SIZE(zfield), p_real_dp, &
            MPI_SUM, process_mpi_root_id, p_comm, p_error)
       p_sum = SUM(pe_sums)
    ELSE
       p_sum = SUM(zfield)
    END IF
#else
    p_sum = SUM(zfield)
#endif

  END FUNCTION p_global_sum_1d

  FUNCTION p_field_sum_1d (zfield, comm) RESULT (p_sum)

    REAL(dp),          INTENT(in) :: zfield(:)
    INTEGER, OPTIONAL, INTENT(in) :: comm
    REAL(dp)                      :: p_sum (SIZE(zfield))

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    IF (my_process_is_mpi_all_parallel()) THEN
       CALL MPI_REDUCE (zfield, p_sum, SIZE(zfield), p_real_dp, &
            MPI_SUM, process_mpi_root_id, p_comm, p_error)
       IF (.NOT. my_process_is_stdio()) p_sum = 0.0_dp
    ELSE
       p_sum = zfield
    END IF
#else
    p_sum = zfield
#endif

  END FUNCTION p_field_sum_1d

  !> computes a global maximum of real numbers
  !
  ! @param[out]   proc_id  (Optional:) PE number of maximum value
  ! @param[inout] keyval   (Optional:) additional meta information
  ! @param[in]    root     (Optional:) root PE, otherwise we perform an
  !                                    ALL-TO-ALL operation
  !
  ! The parameter @p keyval can be used to communicate
  ! additional data on the maximum value, e.g., the level
  ! index where the maximum occurred.
  !
  FUNCTION p_max_0d (zfield, proc_id, keyval, comm, root) RESULT (p_max)

    REAL(dp)                         :: p_max
    REAL(dp),          INTENT(in)    :: zfield
    INTEGER, OPTIONAL, INTENT(inout) :: proc_id
    INTEGER, OPTIONAL, INTENT(inout) :: keyval
    INTEGER, OPTIONAL, INTENT(in)    :: root
    INTEGER, OPTIONAL, INTENT(in)    :: comm
#ifndef NOMPI
    INTEGER  :: p_comm
    INTEGER  :: meta_info, ikey
    REAL(dp) :: in_val(2), rcv_val(2)

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    IF (my_process_is_mpi_all_parallel()) THEN

      IF (PRESENT(proc_id) .OR. PRESENT(keyval)) THEN
        ! encode meta information
        meta_info = 0
        IF (PRESENT(keyval))  meta_info = meta_info + keyval*process_mpi_all_size
        IF (PRESENT(proc_id)) meta_info = meta_info + proc_id
        ! use MPI_MAXLOC to transfer additional data
        in_val(1) = zfield
        in_val(2) = REAL( meta_info, wp )
        IF (PRESENT(root)) THEN
          CALL MPI_REDUCE (in_val, rcv_val, 1, MPI_2DOUBLE_PRECISION, &
            MPI_MAXLOC, root, p_comm, p_error)
        ELSE
          CALL MPI_ALLREDUCE (in_val, rcv_val, 1, MPI_2DOUBLE_PRECISION, &
            MPI_MAXLOC, p_comm, p_error)
        END IF
        ! decode meta info:
        p_max = rcv_val(1)
        ikey = INT(rcv_val(2)+0.5)/process_mpi_all_size
        IF (PRESENT(keyval))  keyval  = ikey
        IF (PRESENT(proc_id)) proc_id = INT(rcv_val(2)+0.5) - ikey
      ELSE
        ! compute simple (standard) maximum
        IF (PRESENT(root)) THEN
          CALL MPI_REDUCE (zfield, p_max, 1, p_real_dp, &
            MPI_MAX, root, p_comm, p_error)
        ELSE
          CALL MPI_ALLREDUCE (zfield, p_max, 1, p_real_dp, &
            MPI_MAX, p_comm, p_error)
        END IF
     END IF
    ELSE
       p_max = zfield
    END IF
#else
    p_max = zfield
#endif

  END FUNCTION p_max_0d

  FUNCTION p_max_int_0d (zfield, comm) RESULT (p_max)

    INTEGER                          :: p_max
    INTEGER,           INTENT(in)    :: zfield
    INTEGER, OPTIONAL, INTENT(in)    :: comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    IF (my_process_is_mpi_all_parallel()) THEN
       CALL MPI_ALLREDUCE (zfield, p_max, 1, p_int, &
            MPI_MAX, p_comm, p_error)
    ELSE
       p_max = zfield
    END IF
#else
    p_max = zfield
#endif

  END FUNCTION p_max_int_0d

  !> computes a global maximum of a real number array
  !
  ! @param[out]   proc_id  (Optional:) PE number of maximum value
  ! @param[inout] keyval   (Optional:) additional meta information
  ! @param[in]    root     (Optional:) root PE, otherwise we perform an
  !                                    ALL-TO-ALL operation
  !
  ! The parameter @p keyval can be used to communicate
  ! additional data on the maximum value, e.g., the level
  ! index where the maximum occurred.
  !
  FUNCTION p_max_1d (zfield, proc_id, keyval, comm, root) RESULT (p_max)

    REAL(dp),          INTENT(in)    :: zfield(:)
    INTEGER, OPTIONAL, INTENT(inout) :: proc_id(SIZE(zfield))
    INTEGER, OPTIONAL, INTENT(inout) :: keyval(SIZE(zfield))
    INTEGER, OPTIONAL, INTENT(in)    :: root
    INTEGER, OPTIONAL, INTENT(in)    :: comm
    REAL(dp)                         :: p_max (SIZE(zfield))

#ifndef NOMPI
    INTEGER :: p_comm
    INTEGER  :: meta_info, ikey, j
    REAL(dp) :: in_val(2, SIZE(zfield)), rcv_val(2, SIZE(zfield))

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    IF (my_process_is_mpi_all_parallel()) THEN

      IF (PRESENT(proc_id) .OR. PRESENT(keyval)) THEN
        ! encode meta information
        in_val(1,:) = zfield(:)
        DO j=1, SIZE(zfield)
          meta_info = 0
          IF (PRESENT(keyval))  meta_info = meta_info + keyval(j)*process_mpi_all_size
          IF (PRESENT(proc_id)) meta_info = meta_info + proc_id(j)
          in_val(2,j) = REAL( meta_info, wp )
        END DO
        ! use MPI_MAXLOC to transfer additional data
        IF (PRESENT(root)) THEN
          CALL MPI_REDUCE (in_val, rcv_val, SIZE(zfield), MPI_2DOUBLE_PRECISION, &
            &              MPI_MAXLOC, root, p_comm, p_error)
        ELSE
          CALL MPI_ALLREDUCE (in_val, rcv_val, SIZE(zfield), MPI_2DOUBLE_PRECISION, &
            &                 MPI_MAXLOC, p_comm, p_error)
        END IF
        ! decode meta info:
        p_max = rcv_val(1,:)
        DO j=1, SIZE(zfield)
          ikey = INT(rcv_val(2,j)+0.5)/process_mpi_all_size
          IF (PRESENT(keyval))  keyval(j)  = ikey
          IF (PRESENT(proc_id)) proc_id(j) = &
            & INT(rcv_val(2,j)+0.5) - ikey*process_mpi_all_size
        END DO
      ELSE
        ! compute simple (standard) maximum
        IF (PRESENT(root)) THEN
          CALL MPI_REDUCE (zfield, p_max, SIZE(zfield), p_real_dp, &
            &              MPI_MAX, root, p_comm, p_error)
        ELSE
          CALL MPI_ALLREDUCE (zfield, p_max, SIZE(zfield), p_real_dp, &
            &                 MPI_MAX, p_comm, p_error)
        END IF
      END IF

    ELSE
       p_max = zfield
    END IF
#else
    p_max = zfield
#endif

  END FUNCTION p_max_1d

  ! Computes maximum of a 1D field of integers.
  !
  ! @param[in] root Optional root PE, otherwise we perform an
  !                 ALL-TO-ALL operation
  FUNCTION p_max_int_1d (zfield, comm, root) RESULT (p_max)

    INTEGER,           INTENT(in) :: zfield(:)
    INTEGER, OPTIONAL, INTENT(in) :: comm
    INTEGER, OPTIONAL, INTENT(in) :: root
    INTEGER                       :: p_max (SIZE(zfield))

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    IF (my_process_is_mpi_all_parallel()) THEN
      IF (PRESENT(root)) THEN
        CALL MPI_REDUCE (zfield, p_max, SIZE(zfield), p_int, &
          MPI_MAX, root, p_comm, p_error)
      ELSE
        CALL MPI_ALLREDUCE (zfield, p_max, SIZE(zfield), p_int, &
          MPI_MAX, p_comm, p_error)
      END IF ! if present(root)
    ELSE
       p_max = zfield
    END IF
#else
    p_max = zfield
#endif

  END FUNCTION p_max_int_1d

  FUNCTION p_max_2d (zfield, comm) RESULT (p_max)

    REAL(dp),          INTENT(in) :: zfield(:,:)
    INTEGER, OPTIONAL, INTENT(in) :: comm
    REAL(dp)                      :: p_max (SIZE(zfield,1),SIZE(zfield,2))

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    IF (my_process_is_mpi_all_parallel()) THEN
       CALL MPI_ALLREDUCE (zfield, p_max, SIZE(zfield), p_real_dp, &
            MPI_MAX, p_comm, p_error)
    ELSE
       p_max = zfield
    END IF
#else
    p_max = zfield
#endif

  END FUNCTION p_max_2d

  FUNCTION p_max_3d (zfield, comm) RESULT (p_max)

    REAL(dp),          INTENT(in) :: zfield(:,:,:)
    INTEGER, OPTIONAL, INTENT(in) :: comm
    REAL(dp)                      :: p_max (SIZE(zfield,1),SIZE(zfield,2)&
                                           ,SIZE(zfield,3))

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    IF (my_process_is_mpi_all_parallel()) THEN
       CALL MPI_ALLREDUCE (zfield, p_max, SIZE(zfield), p_real_dp, &
            MPI_MAX, p_comm, p_error)
    ELSE
       p_max = zfield
    END IF
#else
    p_max = zfield
#endif

  END FUNCTION p_max_3d

  FUNCTION p_min_0d (zfield, comm) RESULT (p_min)

    REAL(dp),          INTENT(in) :: zfield
    INTEGER, OPTIONAL, INTENT(in) :: comm
    REAL(dp)                      :: p_min
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    IF (my_process_is_mpi_all_parallel()) THEN
       CALL MPI_ALLREDUCE (zfield, p_min, 1, p_real_dp, &
            MPI_MIN, p_comm, p_error)
    ELSE
       p_min = zfield
    END IF
#else
    p_min = zfield
#endif

  END FUNCTION p_min_0d

  FUNCTION p_min_int_0d (zfield, comm) RESULT (p_min)

    INTEGER,           INTENT(in) :: zfield
    INTEGER, OPTIONAL, INTENT(in) :: comm
    INTEGER                       :: p_min
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    IF (my_process_is_mpi_all_parallel()) THEN
       CALL MPI_ALLREDUCE (zfield, p_min, 1, p_int, &
            MPI_MIN, p_comm, p_error)
    ELSE
       p_min = zfield
    END IF
#else
    p_min = zfield
#endif

  END FUNCTION p_min_int_0d

  FUNCTION p_min_1d (zfield, comm) RESULT (p_min)

    REAL(dp),          INTENT(in) :: zfield(:)
    INTEGER, OPTIONAL, INTENT(in) :: comm
    REAL(dp)                      :: p_min (SIZE(zfield))

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    IF (my_process_is_mpi_all_parallel()) THEN
       CALL MPI_ALLREDUCE (zfield, p_min, SIZE(zfield), p_real_dp, &
            MPI_MIN, p_comm, p_error)
    ELSE
       p_min = zfield
    END IF
#else
    p_min = zfield
#endif

  END FUNCTION p_min_1d

  FUNCTION p_min_int_1d (zfield, comm) RESULT (p_min)

    INTEGER,           INTENT(in) :: zfield(:)
    INTEGER, OPTIONAL, INTENT(in) :: comm
    INTEGER                       :: p_min (SIZE(zfield))

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    IF (my_process_is_mpi_all_parallel()) THEN
       CALL MPI_ALLREDUCE (zfield, p_min, SIZE(zfield), p_int, &
            MPI_MIN, p_comm, p_error)
    ELSE
       p_min = zfield
    END IF
#else
    p_min = zfield
#endif

  END FUNCTION p_min_int_1d

  FUNCTION p_min_2d (zfield, comm) RESULT (p_min)

    REAL(dp),          INTENT(in) :: zfield(:,:)
    INTEGER, OPTIONAL, INTENT(in) :: comm
    REAL(dp)                      :: p_min (SIZE(zfield,1),SIZE(zfield,2))

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    IF (my_process_is_mpi_all_parallel()) THEN
       CALL MPI_ALLREDUCE (zfield, p_min, SIZE(zfield), p_real_dp, &
            MPI_MIN, p_comm, p_error)
    ELSE
       p_min = zfield
    END IF
#else
    p_min = zfield
#endif

  END FUNCTION p_min_2d

  FUNCTION p_min_3d (zfield, comm) RESULT (p_min)

    REAL(dp),          INTENT(in) :: zfield(:,:,:)
    INTEGER, OPTIONAL, INTENT(in) :: comm
    REAL(dp)                      :: p_min (SIZE(zfield,1),SIZE(zfield,2)&
                                           ,SIZE(zfield,3))
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    IF (my_process_is_mpi_all_parallel()) THEN
       CALL MPI_ALLREDUCE (zfield, p_min, SIZE(zfield), p_real_dp, &
            MPI_MIN, p_comm, p_error)
    ELSE
       p_min = zfield
    END IF
#else
    p_min = zfield
#endif

  END FUNCTION p_min_3d


  ! Computes maximum of a 1D field of integers.
  !
  ! @param[in] root Optional root PE, otherwise we perform an
  !                 ALL-TO-ALL operation
  SUBROUTINE p_allreduce_max_int_1d (zfield, comm, root) 
    INTEGER, INTENT(inout) :: zfield(:)
    INTEGER, INTENT(in) :: comm
    INTEGER, INTENT(in), OPTIONAL :: root

#if !defined(NOMPI)
    INTEGER :: p_comm, p_error

    p_comm = comm
    IF (PRESENT(root)) THEN
      CALL MPI_REDUCE (MPI_IN_PLACE, zfield, SIZE(zfield), p_int, &
        MPI_MAX, root, p_comm, p_error)
    ELSE
      CALL MPI_ALLREDUCE (MPI_IN_PLACE, zfield, SIZE(zfield), p_int, &
        MPI_MAX, p_comm, p_error)
    END IF ! if present(root)
#endif
  END SUBROUTINE p_allreduce_max_int_1d


  SUBROUTINE p_gather_real_0d1d (sendbuf, recvbuf, p_dest, comm)
    REAL(dp),          INTENT(inout) :: sendbuf, recvbuf(:)
    INTEGER,           INTENT(in) :: p_dest
    INTEGER, OPTIONAL, INTENT(in) :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

     CALL MPI_GATHER(sendbuf, 1, p_real_dp, &
                     recvbuf, 1, p_real_dp, &
                     p_dest, p_comm, p_error)
#else
     recvbuf(:) = sendbuf
#endif
   END SUBROUTINE p_gather_real_0d1d

  SUBROUTINE p_gather_real_1d2d (sendbuf, recvbuf, p_dest, comm)
    REAL(dp),          INTENT(inout) :: sendbuf(:), recvbuf(:,:)
    INTEGER,           INTENT(in) :: p_dest
    INTEGER, OPTIONAL, INTENT(in) :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

     CALL MPI_GATHER(sendbuf, SIZE(sendbuf), p_real_dp, &
                     recvbuf, SIZE(sendbuf), p_real_dp, &
                     p_dest, p_comm, p_error)
#else
     recvbuf(:,LBOUND(recvbuf,2)) = sendbuf(:)
#endif
   END SUBROUTINE p_gather_real_1d2d

   SUBROUTINE p_gather_real_5d6d (sendbuf, recvbuf, p_dest, comm)

     REAL(dp)                      :: sendbuf(:,:,:,:,:), recvbuf(:,:,:,:,:,:)
     INTEGER,           INTENT(in) :: p_dest
     INTEGER, OPTIONAL, INTENT(in) :: comm

#ifndef NOMPI
     INTEGER :: p_comm

     IF (PRESENT(comm)) THEN
       p_comm = comm
     ELSE
       p_comm = process_mpi_all_comm
     ENDIF

     CALL MPI_GATHER(sendbuf, SIZE(sendbuf), p_real_dp, &
                     recvbuf, SIZE(sendbuf), p_real_dp, &
                     p_dest, p_comm, p_error)

     IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a)') ' p_gather_real_5d6d failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       STOP
     END IF

#else
     recvbuf(:,:,:,:,:,LBOUND(recvbuf,6)) = sendbuf(:,:,:,:,:)
#endif
   END SUBROUTINE p_gather_real_5d6d


   SUBROUTINE p_gather_int (sendbuf, recvbuf, p_dest, comm)
     INTEGER,           INTENT(inout) :: sendbuf, recvbuf(:)
     INTEGER,           INTENT(in) :: p_dest
     INTEGER, OPTIONAL, INTENT(in) :: comm
     
#ifndef NOMPI
     CHARACTER(*), PARAMETER :: routine = TRIM("mo_mpi:p_gather_int")
     INTEGER :: p_comm
     
     IF (PRESENT(comm)) THEN
       p_comm = comm
     ELSE
       p_comm = process_mpi_all_comm
     ENDIF
     
     CALL MPI_GATHER(sendbuf, 1, MPI_INTEGER, &
       &             recvbuf, 1, MPI_INTEGER, &
       &             p_dest, p_comm, p_error)
     IF (p_error /=  MPI_SUCCESS) CALL finish (routine, 'Error in MPI_GATHER operation!')
#else
     recvbuf(:) = sendbuf
#endif
   END SUBROUTINE p_gather_int


   SUBROUTINE p_gatherv_int (sendbuf, sendcount, recvbuf, recvcounts, &
     &                       displs, p_dest, comm)
     INTEGER,           INTENT(in)    :: sendbuf(:), sendcount
     INTEGER,           INTENT(inout) :: recvbuf(:), recvcounts(:), displs(:)
     INTEGER,           INTENT(in)    :: p_dest
     INTEGER, OPTIONAL, INTENT(in)    :: comm
     
#ifndef NOMPI
     CHARACTER(*), PARAMETER :: routine = TRIM("mo_mpi:p_gatherv_int")
     INTEGER :: p_comm
     
     IF (PRESENT(comm)) THEN
       p_comm = comm
     ELSE
       p_comm = process_mpi_all_comm
     ENDIF
     
     CALL MPI_GATHERV(sendbuf, sendcount,  MPI_INTEGER, &
       &              recvbuf, recvcounts, displs, MPI_INTEGER, &
       &              p_dest, p_comm, p_error)
     IF (p_error /=  MPI_SUCCESS) CALL finish (routine, 'Error in MPI_GATHERV operation!')
#else
     recvbuf((displs(1)+1):(displs(1)+sendcount)) = sendbuf(1:sendcount)
#endif
   END SUBROUTINE p_gatherv_int 

   
   SUBROUTINE p_gatherv_real2D1D (sendbuf, sendcount, recvbuf, recvcounts, displs, p_dest, comm)
     REAL(dp),          INTENT(IN)    :: sendbuf(:,:)
     INTEGER,           INTENT(IN)    :: sendcount
     REAL(dp),          INTENT(INOUT) :: recvbuf(:)
     INTEGER,           intent(IN)    :: recvcounts(:), displs(:)
     INTEGER,           INTENT(in)    :: p_dest
     INTEGER,           INTENT(in)    :: comm

#if !defined(NOMPI)
     CHARACTER(*), PARAMETER :: routine = TRIM("mo_mpi:p_gatherv_real2D1D")
     INTEGER :: p_comm, p_error

     p_comm = comm
     CALL MPI_GATHERV(sendbuf, sendcount, p_real_dp,   &    ! sendbuf, sendcount, sendtype
       &              recvbuf, recvcounts, displs,     &    ! recvbuf, recvcounts, displs
       &              p_real_dp, p_dest, p_comm_work, p_error)  ! recvtype, root, comm, error
     IF (p_error /=  MPI_SUCCESS) CALL finish (routine, 'Error in MPI_GATHERV operation!')
#else
     recvbuf(:) = RESHAPE(sendbuf, (/ SIZE(recvbuf) /) )
#endif
   END SUBROUTINE p_gatherv_real2D1D


  SUBROUTINE p_gatherv_int2D1D (sendbuf, sendcount, recvbuf, recvcounts, displs, p_dest, comm)
    INTEGER,           INTENT(IN)    :: sendbuf(:,:)
    INTEGER,           INTENT(IN)    :: sendcount
    INTEGER,           INTENT(INOUT) :: recvbuf(:)
    INTEGER,           INTENT(IN)    :: recvcounts(:), displs(:)
    INTEGER,           INTENT(in)    :: p_dest
    INTEGER,           INTENT(in)    :: comm

#if !defined(NOMPI)
    CHARACTER(*), PARAMETER :: routine = TRIM("mo_mpi:p_gatherv_int2D1D")
    INTEGER :: p_comm, p_error

    p_comm = comm
    CALL MPI_GATHERV(sendbuf, sendcount, p_int,       &    ! sendbuf, sendcount, sendtype
      &              recvbuf, recvcounts, displs,     &    ! recvbuf, recvcounts, displs
      &              p_int, p_dest, p_comm_work, p_error)  ! recvtype, root, comm, error
    IF (p_error /=  MPI_SUCCESS) CALL finish (routine, 'Error in MPI_GATHERV operation!')
#else
    recvbuf(:) = RESHAPE(sendbuf, (/ SIZE(recvbuf) /) )
#endif
  END SUBROUTINE p_gatherv_int2D1D


   SUBROUTINE p_gatherv_real3D1D (sendbuf, sendcount, recvbuf, recvcounts, displs, p_dest, comm)
     REAL(dp),          INTENT(IN)    :: sendbuf(:,:,:)
     INTEGER,           INTENT(IN)    :: sendcount
     REAL(dp),          INTENT(INOUT) :: recvbuf(:)
     INTEGER,           intent(IN)    :: recvcounts(:), displs(:)
     INTEGER,           INTENT(in)    :: p_dest
     INTEGER,           INTENT(in)    :: comm

#if !defined(NOMPI)
     CHARACTER(*), PARAMETER :: routine = TRIM("mo_mpi:p_gatherv_real2D1D")
     INTEGER :: p_comm, p_error

     p_comm = comm
     CALL MPI_GATHERV(sendbuf, sendcount, p_real_dp,   &    ! sendbuf, sendcount, sendtype
       &              recvbuf, recvcounts, displs,     &    ! recvbuf, recvcounts, displs
       &              p_real_dp, p_dest, p_comm_work, p_error)  ! recvtype, root, comm, error
     IF (p_error /=  MPI_SUCCESS) CALL finish (routine, 'Error in MPI_GATHERV operation!')
#else
     recvbuf(:) = RESHAPE(sendbuf, (/ SIZE(recvbuf) /) )
#endif
   END SUBROUTINE p_gatherv_real3D1D


   SUBROUTINE p_scatterv_real1D2D (sendbuf, sendcounts, displs, recvbuf, recvcount, p_dest, comm)
     REAL(dp),          INTENT(IN)    :: sendbuf(:)
     INTEGER,           INTENT(IN)    :: sendcounts(:), displs(:)
     REAL(dp),          INTENT(INOUT) :: recvbuf(:,:)
     INTEGER,           INTENT(IN)    :: recvcount
     INTEGER,           INTENT(in)    :: p_dest
     INTEGER,           INTENT(in)    :: comm

#if !defined(NOMPI)
     CHARACTER(*), PARAMETER :: routine = TRIM("mo_mpi:p_scatterv_real1D2D")
     INTEGER :: p_comm, p_error

     p_comm = comm
     CALL MPI_SCATTERV(sendbuf, sendcounts, displs,   &           ! sendbuf, sendcount, displs
       &               p_real_dp, recvbuf, recvcount, &           ! sendtype, recvbuf, recvcounts, 
       &               p_real_dp, p_dest, p_comm_work, p_error)   ! recvtype, root, comm, error
     IF (p_error /=  MPI_SUCCESS) CALL finish (routine, 'Error in MPI_SCATTERV operation!')
#else
     recvbuf(:,:) = RESHAPE(sendbuf, (/ SIZE(recvbuf,1), SIZE(recvbuf,2) /))
#endif
   END SUBROUTINE p_scatterv_real1D2D

   
   SUBROUTINE p_allreduce_minloc(array, ntotal, comm)

     INTEGER,           INTENT(IN)    :: ntotal     
     REAL,              INTENT(INOUT) :: array(2,ntotal)
     INTEGER, OPTIONAL, INTENT(IN)    :: comm
     
#ifndef NOMPI
     INTEGER :: p_comm, ierr
     CHARACTER(*), PARAMETER :: routine = TRIM("mo_mpi:p_allreduce_minloc")
     
     IF (PRESENT(comm)) THEN
       p_comm = comm
     ELSE
       p_comm = process_mpi_all_comm
     ENDIF
     CALL MPI_ALLREDUCE( MPI_IN_PLACE, array, ntotal, MPI_2REAL, MPI_MINLOC, &
       &                 p_comm, ierr )
     IF (ierr /= 0)  &
       CALL finish (routine, 'Error in MPI_ALLREDUCE operation!')
#else
     ! do nothing
#endif
   END SUBROUTINE p_allreduce_minloc


   !-------------------------------------------------------------------------
   !> Collects a 2D integer array, organized as (nproma, nblks), and returns
   !  the global array. From each PE we may receive a different number of 
   !  entries. The ordering of the array is defined by the argument "iowner".
   !
   SUBROUTINE p_gather_field_2d_int(total_dim, nlocal_pts, owner, array)

     INTEGER, INTENT(IN)     :: total_dim              ! dimension of global vector             
     INTEGER, INTENT(IN)     :: nlocal_pts(:)          ! number of points located on each patch 
     INTEGER, INTENT(IN)     :: owner(:)               ! for each point: owning process         
     INTEGER, INTENT(INOUT)  :: array(:,:)             ! gathered output data                   

     !-----------------------------------------------------------------------

#ifndef NOMPI
     ! Local Parameters:
     REAL(wp) :: r_array(1,1,UBOUND(array,1),UBOUND(array,2))

     r_array(1,1,:,:) = REAL(array(:,:),wp)
     CALL p_gather_field_4d(total_dim, nlocal_pts, owner, r_array)
     array(:,:) = INT(r_array(1,1,:,:))
#else
     ! do nothing in the sequential case
#endif

   END SUBROUTINE p_gather_field_2d_int


   !-------------------------------------------------------------------------
   !> Collects a 3D integer array, organized as (:,nproma, nblks), and returns
   !  the global array. From each PE we may receive a different number of 
   !  entries. The ordering of the array is defined by the argument "iowner". 
   !
   SUBROUTINE p_gather_field_3d_int(total_dim, nlocal_pts, owner, array)

     INTEGER,  INTENT(IN)    :: total_dim              ! dimension of global vector             
     INTEGER,  INTENT(IN)    :: nlocal_pts(:)          ! number of points located on each patch 
     INTEGER,  INTENT(IN)    :: owner(:)               ! for each point: owning process         
     INTEGER,  INTENT(INOUT) :: array(:,:,:)           ! gathered output data                   

     !-----------------------------------------------------------------------

#ifndef NOMPI
     ! Local Parameters:
     REAL(wp) :: r_array(1,UBOUND(array,1),UBOUND(array,2),UBOUND(array,3))

     r_array(1,:,:,:) = REAL(array(:,:,:),wp)
     CALL p_gather_field_4d(total_dim, nlocal_pts, owner, r_array)
     array(:,:,:) = INT(r_array(1,:,:,:))
#else
     ! do nothing in the sequential case
#endif

   END SUBROUTINE p_gather_field_3d_int


   !-------------------------------------------------------------------------
   !> Collects a 3D REAL array, organized as (:,nproma, nblks), and returns
   !  the global array. From each PE we may receive a different number of 
   !  entries. The ordering of the array is defined by the argument "iowner". 
   !
   SUBROUTINE p_gather_field_3d(total_dim, nlocal_pts, owner, array)

     INTEGER,  INTENT(IN)    :: total_dim               ! dimension of global vector             
     INTEGER,  INTENT(IN)    :: nlocal_pts(:)           ! number of points located on each patch 
     INTEGER,  INTENT(IN)    :: owner(:)                ! for each point: owning process         
     REAL(dp), INTENT(INOUT) :: array(:,:,:)            ! gathered output data                   

     !-----------------------------------------------------------------------

#ifndef NOMPI
     ! Local Parameters:
     REAL(wp) :: r_array(1,UBOUND(array,1),UBOUND(array,2),UBOUND(array,3))

     r_array(1,:,:,:) = array(:,:,:)
     CALL p_gather_field_4d(total_dim, nlocal_pts, owner, r_array)
     array(:,:,:) = r_array(1,:,:,:)

#else
     ! do nothing in the sequential case
#endif

   END SUBROUTINE p_gather_field_3d


   !-------------------------------------------------------------------------
   !> Collects a 4D REAL array, organized as (:,:,nproma, nblks), and returns
   !  the global array. From each PE we may receive a different number of 
   !  entries. The ordering of the array is defined by the argument "iowner". 
   !
   SUBROUTINE p_gather_field_4d(total_dim, nlocal_pts, owner, array)

     INTEGER,  INTENT(IN)    :: total_dim                ! dimension of global vector
     INTEGER,  INTENT(IN)    :: nlocal_pts(:)            ! number of points located on each patch
     INTEGER,  INTENT(IN)    :: owner(:)                 ! for each point: owning process
     REAL(dp), INTENT(INOUT) :: array(:,:,:,:)           ! gathered output data                   

     !-----------------------------------------------------------------------

#ifndef NOMPI
     ! Local Parameters:
     REAL(wp) :: tot_array(UBOUND(array,1),UBOUND(array,2),UBOUND(array,3)*UBOUND(array,4))
     INTEGER :: n, i, mpierr

     ! Check if result fits into output
     IF(UBOUND(tot_array,3) < total_dim) &
       CALL finish('gather_field','Output array too small')

     ! Insert my contribution into tot_array at the final location
     tot_array(:,:,:) = 0._wp
     n = 0
     DO i = 1, total_dim
       IF(owner(i) == p_pe_work) THEN
         tot_array(:,:,i) = array(:,:,MOD(n,UBOUND(array,3))+1,n/UBOUND(array,3)+1)
         n = n+1
       ENDIF
    ENDDO

    ! Plausibility check
    IF(n /= nlocal_pts(p_pe_work+1)) &
       CALL finish('gather_field','nlocal_pts inconsistent')

    ! Take the global sum over tot_array, the result goes back into array

    CALL MPI_Allreduce(tot_array, array, SIZE(tot_array), p_real_dp, MPI_SUM, p_comm_work, mpierr)
#else
     ! do nothing in the sequential case
#endif

   END SUBROUTINE p_gather_field_4d
   
   
   ! Commits a user-defined MPI type
   !
   FUNCTION p_commit_type_struct(oldtypes, blockcounts) RESULT(newtype)
     INTEGER :: newtype
     INTEGER, INTENT(IN) :: oldtypes(2), blockcounts(2)
#if !defined(NOMPI)
     INTEGER :: ierr, offsets(2), extent
     ! define structured type and commit it  
     CALL MPI_TYPE_EXTENT(oldtypes(1), extent, ierr) 
     offsets(:) = (/ 0, extent /)
     CALL MPI_TYPE_STRUCT(2, blockcounts, offsets, oldtypes, newtype, ierr) 
     CALL MPI_TYPE_COMMIT(newtype, ierr) 
#else
     newtype = 0
#endif
   END FUNCTION p_commit_type_struct


   SUBROUTINE p_alltoall_int (sendbuf, recvbuf, comm)
     INTEGER,           INTENT(inout) :: sendbuf(:), recvbuf(:)
     INTEGER,           INTENT(in) :: comm
#if !defined(NOMPI)
     CHARACTER(*), PARAMETER :: routine = TRIM("mo_mpi:p_alltoall_int")
     INTEGER :: p_comm, p_error

     p_comm = comm
     CALL MPI_ALLTOALL(sendbuf, 1, p_int, recvbuf, 1, p_int, p_comm, p_error)
     IF (p_error /=  MPI_SUCCESS) CALL finish (routine, 'Error in MPI_ALLTOALL operation!')
#else
     recvbuf(:) = sendbuf
#endif
   END SUBROUTINE p_alltoall_int


  SUBROUTINE p_clear_request(request)
    INTEGER, INTENT(INOUT) :: request
#if !defined(NOMPI)
    request = MPI_REQUEST_NULL
#endif
  END SUBROUTINE p_clear_request


  FUNCTION p_mpi_wtime()
    REAL :: p_mpi_wtime
#ifndef NOMPI
    p_mpi_wtime = MPI_Wtime()
#else
    p_mpi_wtime = 0.
#endif
  END FUNCTION p_mpi_wtime

END MODULE mo_mpi
