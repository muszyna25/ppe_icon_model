!>
!!     Contains namelists for parallel run control.
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
  USE mo_parallel_configuration

  IMPLICIT NONE

  PRIVATE
  ! Exported routines:
  PUBLIC :: read_parallel_namelist, fill_parallel_nml_configure

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

  ! ------------------------------------------------------------------------
  ! 1.0 Namelist variables and auxiliary parameters
  ! ------------------------------------------------------------------------

  ! Number of rows of ghost cells
  INTEGER :: nml_n_ghost_rows

  INTEGER :: nml_division_method

  ! Flag if checks in a verification run should be logged

  LOGICAL :: nml_l_log_checks

  ! Flag if fast but nonreproducible sum should be used

  LOGICAL :: nml_l_fast_sum

  ! Please note for the following variables: The default settings are for NO_MPI runs!

  ! p_test_run indicates a verification run, i.e. a run where 1 PE runs the complete
  ! model whereas the other PEs do a real parallelized run

  LOGICAL :: nml_p_test_run

  ! if l_test_openmp is set together with p_test_run, then the verification PE uses
  ! only 1 thread. This allows for verifying the OpenMP implementation

  LOGICAL :: nml_l_test_openmp

  ! Processor distribution:
  ! num_test_procs: 0 or 1
  ! num_work_procs: number of procs running in parallel on the model
  ! num_io_procs:   number of procs for I/O
  ! num_test_procs + num_work_procs + num_io_procs = p_nprocs

  INTEGER :: nml_num_test_procs
  INTEGER :: nml_num_work_procs
  INTEGER :: nml_num_io_procs

  ! Note: p_test_pe, p_work_pe0, p_io_pe0 are identical on all PEs

  ! In a verification run, p_test_pe is the number of the PE running the complete model,
  ! otherwise it contains -1

  INTEGER :: nml_p_test_pe      ! Number of test PE
  INTEGER :: nml_p_work_pe0     ! Number of workgroup PE 0 within all PEs
  INTEGER :: nml_p_io_pe0       ! Number of I/O PE 0 within all PEs (p_nprocs if no I/O PEs)

  ! Note: p_n_work, p_pe_work are NOT identical on all PEs

  ! p_n_work: Number of PEs working together:
  ! - num_work_procs for non-verification runs
  ! - num_work_procs for verification runs on pes != p_test_pe
  ! - 1              for verification runs on p_test_pe
  ! - num_io_procs   always on I/O pes

  INTEGER :: nml_p_n_work
  INTEGER :: nml_p_pe_work        ! PE number within work group

  ! Type of parallel I/O
  INTEGER :: nml_pio_type

  ! Type of (halo) communication: 
  ! 1 = synchronous communication with local memory for exchange buffers
  ! 2 = synchronous communication with global memory for exchange buffers
  ! 3 = asynchronous communication within dynamical core with global memory 
  !     for exchange buffers (not yet implemented)
  INTEGER :: nml_itype_comm

  ! Order of send/receive sequence in exchange routines
  ! 1 = irecv, send
  ! 2 = isend, recv
  INTEGER :: nml_iorder_sendrecv

  INTEGER :: nml_radiation_threads
  INTEGER :: nml_nh_stepping_threads

  INTEGER :: nml_nproma    ! inner loop length/vector length

  
  NAMELIST/parallel_ctl/ nml_n_ghost_rows, nml_division_method, &
       &                 nml_p_test_run,   nml_l_test_openmp,   &
                         nml_l_log_checks, nml_l_fast_sum,      &
                         nml_num_io_procs, nml_pio_type,        &
                         nml_itype_comm,   nml_iorder_sendrecv, &
                         nml_radiation_threads, nml_nh_stepping_threads, &
                         nml_nproma


  
  ! Exported variables:
!   PUBLIC ::              nml_n_ghost_rows, nml_division_method, &
!                          nml_p_test_run,   nml_l_test_openmp,   &
!                          nml_l_log_checks, nml_l_fast_sum,      &
!                          nml_num_io_procs, nml_pio_type,        &
!                          nml_itype_comm,   nml_iorder_sendrecv, &
!                          nml_nproma
       
  CONTAINS
  !>
  !!
  !! @par Revision History
  !! Initial version by Rainer Johanni, Nov 2009
  !! Adapted for I/O PEs, Rainer Johanni, Nov 2010
  !!
  SUBROUTINE read_parallel_namelist

    INTEGER :: istat, my_color, peer_comm, p_error
    INTEGER :: funit
    CHARACTER(len=*), PARAMETER ::   &
            &  method_name = 'read_parallel_namelist'

    !--------------------------------------------
    ! set default values
    !--------------------------------------------
    ! Number of rows of ghost cells
    nml_n_ghost_rows = 1
    nml_division_method = div_geometric

    ! Flag if checks in a verification run should be logged
    nml_l_log_checks = .FALSE.

    ! Flag if fast but nonreproducible sum should be used
    nml_l_fast_sum = .FALSE.

    ! Please note for the following variables: The default settings are for NO_MPI runs!
    ! p_test_run indicates a verification run, i.e. a run where 1 PE runs the complete
    ! model whereas the other PEs do a real parallelized run
    nml_p_test_run = .FALSE.

    ! if l_test_openmp is set together with p_test_run, then the verification PE uses
    ! only 1 thread. This allows for verifying the OpenMP implementation
    nml_l_test_openmp = .FALSE.

    ! Processor distribution:
    ! num_test_procs: 0 or 1
    ! num_work_procs: number of procs running in parallel on the model
    ! num_io_procs:   number of procs for I/O
    ! num_test_procs + num_work_procs + num_io_procs = p_nprocs
    nml_num_test_procs = 0
    nml_num_work_procs = 1
    nml_num_io_procs = 0

    ! Note: p_test_pe, p_work_pe0, p_io_pe0 are identical on all PEs
    ! In a verification run, p_test_pe is the number of the PE running the complete model,
    ! otherwise it contains -1
    nml_p_test_pe  = -1      ! Number of test PE
    nml_p_work_pe0 =  0      ! Number of workgroup PE 0 within all PEs
    nml_p_io_pe0   =  1      ! Number of I/O PE 0 within all PEs (p_nprocs if no I/O PEs)

    ! Note: p_n_work, p_pe_work are NOT identical on all PEs
    ! p_n_work: Number of PEs working together:
    ! - num_work_procs for non-verification runs
    ! - num_work_procs for verification runs on pes != p_test_pe
    ! - 1              for verification runs on p_test_pe
    ! - num_io_procs   always on I/O pes
    nml_p_n_work=1
    nml_p_pe_work=0        ! PE number within work group

    ! Type of parallel I/O
    nml_pio_type = 1

    ! Type of (halo) communication: 
    ! 1 = synchronous communication with local memory for exchange buffers
    ! 2 = synchronous communication with global memory for exchange buffers
    ! 3 = asynchronous communication within dynamical core with global memory 
    !     for exchange buffers (not yet implemented)
    nml_itype_comm = 1

    ! Order of send/receive sequence in exchange routines
    ! 1 = irecv, send
    ! 2 = isend, recv
    nml_iorder_sendrecv = 1

    nml_radiation_threads = 1
    nml_nh_stepping_threads = 1

    nml_nproma = 1    ! inner loop length/vector length
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
    
    !-----------------------------------------------------
    ! Store the namelist for restart
    !-----------------------------------------------------
    funit = open_tmpfile()
    WRITE(funit,NML=parallel_ctl)
    CALL store_and_close_namelist(funit, 'parallel_ctl')
    
    !-----------------------------------------------------
    ! Write the final namelist to an ASCII file
    IF (p_pe == p_io) WRITE(nnml_output,nml=parallel_ctl)

    CALL fill_parallel_nml_configure()
    
 END SUBROUTINE read_parallel_namelist
   
   
 SUBROUTINE fill_parallel_nml_configure
    !-----------------------------------------------------
    ! fill the parallel_configuration
    n_ghost_rows        = nml_n_ghost_rows
    division_method     = nml_division_method
    p_test_run          = nml_p_test_run
    l_test_openmp       = nml_l_test_openmp
    l_log_checks        = nml_l_log_checks
    l_fast_sum          = nml_l_fast_sum
    num_io_procs        = nml_num_io_procs
    pio_type            = nml_pio_type
    itype_comm          = nml_itype_comm
    iorder_sendrecv     = nml_iorder_sendrecv
    radiation_threads   = nml_radiation_threads
    nh_stepping_threads = nml_nh_stepping_threads
    nproma              = nml_nproma

    ! check the configuration
    CALL check_parallel_configuration()

END SUBROUTINE fill_parallel_nml_configure

  !-------------------------------------------------------------------------

END MODULE mo_parallel_nml
