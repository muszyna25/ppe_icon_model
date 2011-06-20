!>
!!               Contains parameters for parallel run control.
!!
!!
!! @par Revision History
!! Initial version by Rainer Johanni, Nov 2009
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
MODULE mo_parallel_nml

  USE mo_exception,          ONLY: message, finish
  USE mo_io_units,           ONLY: nnml, nnml_output
  USE mo_namelist,           ONLY: position_nml, POSITIONED
  USE mo_master_nml,         ONLY: lrestart
  USE mo_run_nml,            ONLY: lrestore_states
  USE mo_mpi,                ONLY: p_pe, p_io, p_nprocs, p_all_comm
#ifndef NOMPI
  USE mo_mpi,                ONLY: MPI_COMM_NULL, MPI_COMM_SELF, MPI_UNDEFINED, &
    &   p_comm_work, p_comm_work_test   ! Communicator spanning work group and test PE
  USE mo_exception,          ONLY: message_text
#else
  USE mo_mpi,                ONLY:  p_comm_work, p_comm_work_test
#endif
  USE mo_io_restart_namelist,ONLY: open_tmpfile, store_and_close_namelist,   &
                                 & open_and_restore_namelist, close_tmpfile

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

  ! ------------------------------------------------------------------------
  ! 1.0 Namelist variables and auxiliary parameters
  ! ------------------------------------------------------------------------

  ! Number of rows of ghost cells
  INTEGER :: n_ghost_rows = 1


  ! Division method for area subdivision
  INTEGER, PARAMETER :: div_from_file = 0  ! Read from file
  INTEGER, PARAMETER :: div_geometric = 1  ! Geometric subdivision
  INTEGER, PARAMETER :: div_metis     = 2  ! Use Metis

  INTEGER :: division_method = div_geometric

  ! Flag if checks in a verification run should be logged

  LOGICAL :: l_log_checks = .FALSE.

  ! Flag if fast but nonreproducible sum should be used

  LOGICAL :: l_fast_sum = .FALSE.

  ! Please note for the following variables: The default settings are for NO_MPI runs!

  ! p_test_run indicates a verification run, i.e. a run where 1 PE runs the complete
  ! model whereas the other PEs do a real parallelized run

  LOGICAL :: p_test_run = .FALSE.

  ! if l_test_openmp is set together with p_test_run, then the verification PE uses
  ! only 1 thread. This allows for verifying the OpenMP implementation

  LOGICAL :: l_test_openmp = .FALSE.

  ! Processor distribution:
  ! num_test_procs: 0 or 1
  ! num_work_procs: number of procs running in parallel on the model
  ! num_io_procs:   number of procs for I/O
  ! num_test_procs + num_work_procs + num_io_procs = p_nprocs

  INTEGER :: num_test_procs = 0
  INTEGER :: num_work_procs = 1
  INTEGER :: num_io_procs = 0

  ! Note: p_test_pe, p_work_pe0, p_io_pe0 are identical on all PEs

  ! In a verification run, p_test_pe is the number of the PE running the complete model,
  ! otherwise it contains -1

  INTEGER :: p_test_pe  = -1      ! Number of test PE
  INTEGER :: p_work_pe0 =  0      ! Number of workgroup PE 0 within all PEs
  INTEGER :: p_io_pe0   =  1      ! Number of I/O PE 0 within all PEs (p_nprocs if no I/O PEs)

  ! Note: p_n_work, p_pe_work are NOT identical on all PEs

  ! p_n_work: Number of PEs working together:
  ! - num_work_procs for non-verification runs
  ! - num_work_procs for verification runs on pes != p_test_pe
  ! - 1              for verification runs on p_test_pe
  ! - num_io_procs   always on I/O pes

  INTEGER :: p_n_work=1
  INTEGER :: p_pe_work=0        ! PE number within work group

  ! MPI communicators
  INTEGER :: p_comm_work_io     ! Communicator spanning work group and I/O PEs
  INTEGER :: p_comm_work_2_io   ! Inter(!)communicator work PEs - I/O PEs
  INTEGER :: p_comm_input_bcast ! Communicator for broadcasts in NetCDF input

  ! Type of parallel I/O

  INTEGER :: pio_type = 1

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

  INTEGER :: radiation_threads = 1
  INTEGER :: nh_stepping_threads = 1
  
  NAMELIST/parallel_ctl/ n_ghost_rows, division_method, &
       &                 p_test_run,   l_test_openmp,   &
                         l_log_checks, l_fast_sum,      &
                         num_io_procs, pio_type,        &
                         itype_comm,   iorder_sendrecv, &
                         radiation_threads, nh_stepping_threads


  ! Exported routines:

  PUBLIC :: parallel_nml_setup
  
  ! Exported variables:
#ifdef __OMP_RADIATION__
  PUBLIC :: radiation_threads, nh_stepping_threads
#endif

  PUBLIC :: n_ghost_rows,                                             &
       &    div_from_file, div_geometric, div_metis, division_method, &
       &    l_log_checks, l_fast_sum,                                 &
       &    p_test_run, l_test_openmp,                                &
       &    num_test_procs, num_work_procs, num_io_procs,             &
       &    p_test_pe, p_work_pe0, p_io_pe0,                          &
       &    p_n_work, p_pe_work, p_comm_work, p_comm_work_test,       &
       &    p_comm_work_io, p_comm_work_2_io, p_comm_input_bcast,     &
       &    pio_type, itype_comm, iorder_sendrecv

  CONTAINS
  !>
  !!
  !! @par Revision History
  !! Initial version by Rainer Johanni, Nov 2009
  !! Adapted for I/O PEs, Rainer Johanni, Nov 2010
  !!
  SUBROUTINE parallel_nml_setup

    INTEGER :: istat, my_color, peer_comm, p_error
    INTEGER :: funit

    !----------------------------------------------------------------
    ! If this is a resumed integration, overwrite the defaults above
    ! by values in the previous integration.
    !----------------------------------------------------------------
    IF (lrestart) THEN
      funit = open_and_restore_namelist('parallel_ctl')
      READ(funit,NML=parallel_ctl)
      CALL close_tmpfile(funit)
    END IF

    !--------------------------------------------------------------------
    ! Read user's (new) specifications (Done so far by all MPI processes) 
    !--------------------------------------------------------------------
    CALL position_nml ('parallel_ctl', STATUS=istat)
    SELECT CASE (istat)
    CASE (POSITIONED)
      READ (nnml, parallel_ctl)
    END SELECT

  !------------------------------------------------------------
  ! 3.0 check the consistency of the parameters
  !------------------------------------------------------------

  ! check n_ghost_rows
  IF (n_ghost_rows<1) THEN
    CALL finish('mo_parallel_ctl:parallel_nml_setup', &
         & 'n_ghost_rows<1 in parallel_ctl namelist is not allowed')
  END IF

  ! check division_method
  SELECT CASE (division_method)
  CASE(div_from_file, div_geometric)
    ! ok
  CASE(div_metis)
#ifdef HAVE_METIS
    ! ok
#else
    CALL finish('mo_parallel_ctl:parallel_nml_setup', &
       & 'division_method=div_metis=2 in parallel_ctl namelist is not allowed')
#endif
  CASE DEFAULT
    CALL finish('mo_parallel_ctl:parallel_nml_setup', &
       & 'value of division_method in parallel_ctl namelist is not allowed')
  END SELECT

  ! check p_test_run and num_io_procs
#ifdef NOMPI
  ! Unconditionally set p_test_run to .FALSE. and num_io_procs to 0,
  ! all other variables are already set correctly
  IF (p_test_run) THEN
    CALL message('mo_parallel_ctl:parallel_nml_setup', &
         & 'p_test_run has no effect if the model is compiled with the NOMPI compiler directive')
    CALL message('mo_parallel_ctl:parallel_nml_setup', &
         & '--> p_test_run set to .FALSE.')
    p_test_run = .FALSE.
  END IF
  IF (num_io_procs /= 0) THEN
    CALL message('mo_parallel_ctl:parallel_nml_setup', &
         & 'num_io_procs has no effect if the model is compiled with the NOMPI compiler directive')
    CALL message('mo_parallel_ctl:parallel_nml_setup', &
         & '--> num_io_procs set to 0')
    num_io_procs = 0
  END IF
#else
  ! A run on 1 PE is never a verification run,
  ! correct this if the user should set it differently
  IF (p_test_run .AND. p_nprocs == 1) THEN
    CALL message('mo_parallel_ctl:parallel_nml_setup', &
         & 'p_test_run has no effect if p_nprocs=1')
    CALL message('mo_parallel_ctl:parallel_nml_setup', &
         & '--> p_test_run set to .FALSE.')
     p_test_run = .FALSE.
  ENDIF
  ! for safety only
  IF(num_io_procs < 0) num_io_procs = 0
#endif

  ! check l_test_openmp
#ifndef _OPENMP
  IF (l_test_openmp) THEN
    CALL message('mo_parallel_ctl:parallel_nml_setup', &
         & 'l_test_openmp has no effect if the model is compiled without OpenMP support')
    CALL message('mo_parallel_ctl:parallel_nml_setup', &
         & '--> l_test_openmp set to .FALSE.')
    l_test_openmp = .FALSE.
  END IF
#endif

  ! no checks for l_log_checks

  ! no checks for l_fast_sum

    !-----------------------------------------------------
    ! Store the namelist for restart
    !-----------------------------------------------------
    funit = open_tmpfile()
    WRITE(funit,NML=parallel_ctl)
    CALL store_and_close_namelist(funit, 'parallel_ctl')

    ! Write the final namelist to an ASCII file
    IF (p_pe == p_io) WRITE(nnml_output,nml=parallel_ctl)

!-----------------------------------------------------------------------

  ! Set dependent control variables according
  ! to the (modified) NAMELIST varaibles
  ! -----------------------------------------

#ifndef NOMPI
  ! Set up processor numbers

  IF(p_test_run) THEN
    num_test_procs = 1
  ELSE
    num_test_procs = 0
  ENDIF

  num_work_procs = p_nprocs - num_test_procs - num_io_procs

  ! Check if there are sufficient PEs at all

  IF(num_work_procs < 1) THEN
    CALL finish('mo_parallel_ctl:parallel_nml_setup', &
       & 'not enough processors for given values of p_test_run/num_io_procs')
  ELSE IF (p_test_run .AND. num_work_procs == 1) THEN
    CALL finish('mo_parallel_ctl:parallel_nml_setup', &
       & 'running p_test_run with only 1 work processor does not make sense')
  ENDIF

  WRITE(message_text,'(3(a,i0))') 'Number of procs for test: ',num_test_procs, &
    & ', work: ',num_work_procs,', I/O: ',num_io_procs

  CALL message('mo_parallel_ctl:parallel_nml_setup', message_text)

  ! Set up p_test_pe, p_work_pe0, p_io_pe0 which are identical on all PEs

  IF(p_test_run) THEN
    p_test_pe = 0
  ELSE
    p_test_pe = -1
  ENDIF

  p_work_pe0 = num_test_procs
  p_io_pe0   = num_test_procs + num_work_procs

  ! if OpenMP is used, the test PE uses only 1 thread in order to check
  ! the correctness of the OpenMP implementation
  ! Currently the I/O PEs are also single threaded!
#ifdef _OPENMP
  IF (l_test_openmp .AND. p_pe == p_test_pe) CALL OMP_SET_NUM_THREADS(1)
  IF (p_pe >= p_io_pe0) CALL OMP_SET_NUM_THREADS(1)
#endif

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

  ! Split communicator p_all_comm between test/work/io
  ! to get p_comm_work which is the communicator for
  ! usage WITHIN every group of the 3 different type

  IF(p_pe < p_work_pe0) THEN
    my_color = 1 ! Test PE
  ELSE IF(p_pe < p_io_pe0) THEN
    my_color = 2 ! Work PE
  ELSE
    my_color = 3 ! I/O PE
  ENDIF

  CALL MPI_Comm_split(p_all_comm, my_color, p_pe, p_comm_work, p_error)

  ! Set p_comm_work_test, the communicator spanning work group and test PE

  IF(p_test_run) THEN
    IF(p_pe < p_io_pe0) THEN
      my_color = 1
    ELSE
      my_color = MPI_UNDEFINED ! p_comm_work_test must never be used on I/O PEs
    ENDIF

    CALL MPI_Comm_split(p_all_comm, my_color, p_pe, p_comm_work_test, p_error)
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

    CALL MPI_Comm_split(p_all_comm, my_color, p_pe, p_comm_work_io, p_error)
  ELSE
    ! If no I/O PEs are present, p_comm_work_io must not be used at all
    p_comm_work_io = MPI_COMM_NULL
  ENDIF

  ! Set p_comm_input_bcast, the communicator for broadcasting the NetCDF input

  IF(lrestore_states) THEN
    ! NetCDF input is only read by the test pe and MUST NOT be broadcast
    IF(p_pe == p_test_pe) THEN
      p_comm_input_bcast = MPI_COMM_SELF ! i.e. effectively no broadcast
    ELSE
      p_comm_input_bcast = MPI_COMM_NULL ! Must not be used!
    ENDIF
  ELSE
    IF(p_pe < p_io_pe0) THEN
      IF(p_test_run) THEN
        ! Test PE reads and broadcasts to workers
        p_comm_input_bcast = p_comm_work_test
      ELSE
        ! PE 0 reads and broadcasts
        p_comm_input_bcast = p_comm_work
      ENDIF
    ELSE
      ! I/O PEs never participate in reading
      p_comm_input_bcast = MPI_COMM_NULL
    ENDIF
  ENDIF

  ! Create Intercommunicator work PEs - I/O PEs

  ! From MPI-Report:
  ! Advice to users: We recommend using a dedicated peer communicator, such as a
  ! duplicate of MPI_COMM_WORLD, to avoid trouble with peer communicators.

  ! No idea what these troubles may be in reality, but let us follow the advice

  CALL MPI_Comm_dup(p_all_comm, peer_comm, p_error)

  IF(p_pe /= p_test_pe .AND. num_io_procs>0) THEN

    IF(p_pe < p_io_pe0) THEN
      CALL MPI_Intercomm_create(p_comm_work, 0, peer_comm, p_io_pe0,  1, p_comm_work_2_io, p_error)
    ELSE
      CALL MPI_Intercomm_create(p_comm_work, 0, peer_comm, p_work_pe0,1, p_comm_work_2_io, p_error)
    ENDIF
  ELSE
    ! No Intercommunicator
    p_comm_work_2_io = MPI_COMM_NULL
  ENDIF
#endif

END SUBROUTINE parallel_nml_setup

END MODULE mo_parallel_nml
