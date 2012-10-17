!>
!!     Contains namelists for parallel run control.
!!
!!
!! @par Revision History
!! Initial version by Rainer Johanni, Nov 2009
!!
!! @par Copyright
!! 2002-2011 by DWD and MPI-M
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

  USE mo_io_units,            ONLY: nnml, nnml_output
  USE mo_namelist,            ONLY: position_nml, POSITIONED, open_nml, close_nml
  USE mo_master_control,      ONLY: is_restart_run
  USE mo_mpi,                 ONLY: my_process_is_stdio
  USE mo_io_restart_namelist, ONLY: open_tmpfile, store_and_close_namelist,  &
    & open_and_restore_namelist, close_tmpfile
  USE mo_io_units,           ONLY: filename_max
!   USE mo_impl_constants,     ONLY: max_dom

  USE mo_parallel_config,     ONLY: &
    & config_n_ghost_rows        => n_ghost_rows,        &
    & config_division_method     => division_method,     &
    & config_division_file_name  => division_file_name,  &
    & config_ldiv_phys_dom       => ldiv_phys_dom,       &
    & config_rad_division_file_name  => radiation_division_file_name,  &
    & config_l_log_checks        => l_log_checks,        &
    & config_l_fast_sum          => l_fast_sum,          &
    & config_p_test_run          => p_test_run,          &
    & config_l_test_openmp       => l_test_openmp,       &
    & config_num_io_procs        => num_io_procs,        &
    & config_pio_type            => pio_type,            &
    & config_itype_comm          => itype_comm,          &
    & config_iorder_sendrecv     => iorder_sendrecv,     &
    & config_exch_msgsize        => exch_msgsize,        &
    & config_radiation_threads   => radiation_ompthreads,   &
    & config_nh_stepping_threads => nh_stepping_ompthreads, &
    & config_nproma              => nproma,                 &
    & config_openmp_threads      => openmp_threads,         &
    & config_parallel_radiation_omp => parallel_radiation_omp,  &
    & config_parallel_radiation_mode => parallel_radiation_mode,  &
    & config_test_parallel_radiation=> test_parallel_radiation, &
    & config_use_icon_comm       => use_icon_comm,        &
    & config_icon_comm_debug     => icon_comm_debug,        &
    & div_geometric, check_parallel_configuration,          &
    & config_max_sr_buffer_size => max_send_recv_buffer_size, &
    & config_use_dycore_barrier => use_dycore_barrier,        &
    & config_itype_exch_barrier => itype_exch_barrier,        &
    & config_use_sp_output      => use_sp_output,             &
    & config_icon_comm_method   => icon_comm_method,          &
    & config_max_no_of_comm_var => max_no_of_comm_variables,  &
    & config_max_no_of_comm_proc => max_no_of_comm_processes, &
    & config_max_no_of_comm_patt => max_no_of_comm_patterns,  &
    & config_sync_barrier_mode   => sync_barrier_mode,        &
    & config_max_mpi_message_size => max_mpi_message_size,    &
    & config_use_physics_barrier  => use_physics_barrier

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: read_parallel_namelist

  CHARACTER(len=*), PARAMETER :: version = &
    &  '$Id$'


       
  CONTAINS
  
  !-------------------------------------------------------------------------
  !>
  !! @par Revision History
  !! Initial version by Rainer Johanni, Nov 2009
  !! Adapted for I/O PEs, Rainer Johanni, Nov 2010
  !! Leonidas Linardakis, namelist restructuring, Jul 2011
  SUBROUTINE read_parallel_namelist( filename )

    ! ------------------------------------------------------------------------
    ! Number of rows of ghost cells
    INTEGER :: n_ghost_rows

    INTEGER :: division_method  ! div_from_file = 0  ! Read from file
                                ! div_geometric = 1  ! Geometric subdivision
                                ! div_metis     = 2  ! Use Metis

    CHARACTER(LEN=filename_max) :: division_file_name ! if div_from_file
    CHARACTER(LEN=filename_max) :: radiation_division_file_name ! if div_from_file

    ! Flag if (in case of merged domains) physical domains shall be considered for 
    ! computing the domain decomposition
    LOGICAL :: ldiv_phys_dom

    ! Flag if checks in a verification run should be logged
    LOGICAL :: l_log_checks

    ! Flag if fast but nonreproducible sum should be used
    LOGICAL :: l_fast_sum

    ! Please note for the following variables: The default settings are for NO_MPI runs!

    ! p_test_run indicates a verification run, i.e. a run where 1 PE runs the complete
    ! model whereas the other PEs do a real parallelized run
    LOGICAL :: p_test_run

    LOGICAL :: use_dycore_barrier ! put an mpi barrier before the dycore to synchronize MPI tasks
    LOGICAL :: use_physics_barrier
    INTEGER :: itype_exch_barrier ! 1: put an mpi barrier at the beginning of exchange calls to synchronize MPI tasks
                                  ! 2: put an mpi barrier after MPI_WAIT to synchronize MPI tasks
                                  ! 3: 1+2

    ! if l_test_openmp is set together with p_test_run, then the verification PE uses
    ! only 1 thread. This allows for verifying the OpenMP implementation
    LOGICAL :: l_test_openmp

    LOGICAL :: use_icon_comm
    LOGICAL :: icon_comm_debug
    INTEGER :: max_send_recv_buffer_size
    INTEGER :: max_mpi_message_size
    INTEGER :: icon_comm_method
    INTEGER :: max_no_of_comm_variables
    INTEGER :: max_no_of_comm_processes
    INTEGER :: max_no_of_comm_patterns
    INTEGER :: sync_barrier_mode
    
    ! Type of parallel I/O
    INTEGER :: pio_type
    INTEGER :: num_io_procs
    
    ! Type of (halo) communication: 
    ! 1 = synchronous communication with local memory for exchange buffers
    ! 2 = synchronous communication with global memory for exchange buffers
    ! 3 = asynchronous communication within dynamical core with global memory 
    !     for exchange buffers (not yet implemented)
    INTEGER :: itype_comm

    ! Order of send/receive sequence in exchange routines
    ! 1 = irecv, send
    ! 2 = isend, recv
    ! 3 = irecv, isend
    ! 4 = irecv, send with message blocking (maximum message size is given by exch_msgsize)
    INTEGER :: iorder_sendrecv

    ! maximum message size if iorder_sendrecv = 4
    INTEGER :: exch_msgsize

    !--------------------------------------------
    ! namelist for parallel radiation
    LOGICAL :: parallel_radiation_omp
    INTEGER :: parallel_radiation_mode
    LOGICAL :: test_parallel_radiation
    INTEGER :: radiation_threads
    INTEGER :: nh_stepping_threads
    !--------------------------------------------

    INTEGER :: nproma    ! inner loop length/vector length

    INTEGER :: openmp_threads  

    LOGICAL :: use_sp_output
    
    NAMELIST /parallel_nml/ n_ghost_rows,  division_method, ldiv_phys_dom, &
      & l_log_checks,      l_fast_sum,          &
      & p_test_run,        l_test_openmp,       &
      & num_io_procs,      pio_type,            &
      & itype_comm,        iorder_sendrecv,     &
      & radiation_threads, nh_stepping_threads, &
      & nproma, parallel_radiation_omp,         &
      & parallel_radiation_mode,  use_icon_comm, &
      & test_parallel_radiation, openmp_threads, &
      & icon_comm_debug, max_send_recv_buffer_size, &
      & division_file_name, radiation_division_file_name, use_dycore_barrier, &
      & use_sp_output, itype_exch_barrier, exch_msgsize, &
      & icon_comm_method, max_no_of_comm_variables,       &
      & max_no_of_comm_processes, max_no_of_comm_patterns, &
      & sync_barrier_mode, max_mpi_message_size, use_physics_barrier

    CHARACTER(LEN=*), INTENT(IN) :: filename
    INTEGER :: istat
    INTEGER :: funit
    !0!CHARACTER(len=*), PARAMETER ::   &
    !0!        &  method_name = 'mo_parallel_nml:read_parallel_namelist'

    !--------------------------------------------
    ! set default values
    !--------------------------------------------
    ! Number of rows of ghost cells
    n_ghost_rows = 1
    division_method = div_geometric
    division_file_name = ""
    radiation_division_file_name = ""
     
    ! Flag if (in case of merged domains) physical domains shall be considered for 
    ! computing the domain decomposition
    ldiv_phys_dom = .TRUE.

    ! Flag if checks in a verification run should be logged
    l_log_checks = .FALSE.

    ! Flag if fast but nonreproducible sum should be used
    l_fast_sum = .FALSE.

    ! Please note for the following variables: The default settings are for NO_MPI runs!
    ! p_test_run indicates a verification run, i.e. a run where 1 PE runs the complete
    ! model whereas the other PEs do a real parallelized run
    p_test_run = .FALSE.

    ! The barriers should be used for dedicated tests only, not for production runs
    use_dycore_barrier = config_use_dycore_barrier
    use_physics_barrier= config_use_physics_barrier
    itype_exch_barrier = 0
    
    ! if l_test_openmp is set together with p_test_run, then the verification PE uses
    ! only 1 thread. This allows for verifying the OpenMP implementation
    l_test_openmp = .FALSE.

    use_icon_comm = .FALSE.
    icon_comm_debug = .FALSE.
    max_send_recv_buffer_size = config_max_sr_buffer_size
    max_mpi_message_size      = config_max_mpi_message_size
    icon_comm_method          = config_icon_comm_method
    max_no_of_comm_variables  = config_max_no_of_comm_var
    max_no_of_comm_processes  = config_max_no_of_comm_proc
    max_no_of_comm_patterns   = config_max_no_of_comm_patt
    sync_barrier_mode = 0
          
    ! Type of parallel I/O
    pio_type = 1
    num_io_procs = 0
    
    ! Type of (halo) communication: 
    ! 1 = synchronous communication with local memory for exchange buffers
    ! 2 = synchronous communication with global memory for exchange buffers
    ! 3 = asynchronous communication within dynamical core with global memory 
    !     for exchange buffers (not yet implemented)
    itype_comm = 1

    ! Order of send/receive sequence in exchange routines
    ! 1 = irecv, send
    ! 2 = isend, recv
    ! 3 = irecv, isend
    ! 4 = irecv, send with message blocking (maximum message size is given by exch_msgsize)
    iorder_sendrecv = 1

    ! maximum message size if iorder_sendrecv = 4
    exch_msgsize = 8192  ! this means 64 KB for real*8 

    ! inner loop length/vector length
    nproma = 1
    openmp_threads = -1 ! < 0 means do not use this value
    
    ! parallel_radiation
    parallel_radiation_mode = 0
    parallel_radiation_omp = .false.
    test_parallel_radiation = .false.
    radiation_threads = 1
    nh_stepping_threads = 1

    ! output in single precision
    use_sp_output = .FALSE.

    !----------------------------------------------------------------
    ! If this is a resumed integration, overwrite the defaults above
    ! by values in the previous integration.
    !----------------------------------------------------------------
    IF (is_restart_run()) THEN
      funit = open_and_restore_namelist('parallel_nml')
      READ(funit,NML=parallel_nml)
      CALL close_tmpfile(funit)
    END IF

    !--------------------------------------------------------------------
    ! Read user's (new) specifications (Done so far by all MPI processes) 
    !--------------------------------------------------------------------
    CALL open_nml(TRIM(filename))
    CALL position_nml ('parallel_nml', STATUS=istat)
    SELECT CASE (istat)
    CASE (POSITIONED)
      READ (nnml, parallel_nml)
    END SELECT
    CALL close_nml
    
    !-----------------------------------------------------
    ! Store the namelist for restart
    !-----------------------------------------------------
    IF(my_process_is_stdio())  THEN
      funit = open_tmpfile()
      WRITE(funit,NML=parallel_nml)
      CALL store_and_close_namelist(funit, 'parallel_nml')
    ENDIF
    !-----------------------------------------------------
    ! Write the final namelist to an ASCII file
    IF (my_process_is_stdio()) WRITE(nnml_output,nml=parallel_nml)

    !-----------------------------------------------------
    ! fill_parallel_nml_configure       
    config_n_ghost_rows        = n_ghost_rows
    config_division_method     = division_method
    config_division_file_name  = division_file_name
    config_ldiv_phys_dom       = ldiv_phys_dom
    config_rad_division_file_name  = radiation_division_file_name
    config_l_log_checks        = l_log_checks
    config_l_fast_sum          = l_fast_sum
    config_p_test_run          = p_test_run
    config_l_test_openmp       = l_test_openmp
    config_num_io_procs        = num_io_procs
    config_pio_type            = pio_type
    config_itype_comm          = itype_comm
    config_iorder_sendrecv     = iorder_sendrecv
    config_exch_msgsize        = exch_msgsize
    config_radiation_threads   = radiation_threads
    config_nh_stepping_threads = nh_stepping_threads
    config_nproma              = nproma
    config_openmp_threads         = openmp_threads
    
    config_use_icon_comm       = use_icon_comm
    config_icon_comm_debug     = icon_comm_debug
    config_icon_comm_method    = icon_comm_method
    config_max_no_of_comm_var  = max_no_of_comm_variables
    config_max_no_of_comm_proc = max_no_of_comm_processes
    config_max_no_of_comm_patt = max_no_of_comm_patterns
    config_sync_barrier_mode   = sync_barrier_mode
    config_parallel_radiation_omp  = parallel_radiation_omp
    config_parallel_radiation_mode = parallel_radiation_mode
    config_test_parallel_radiation = test_parallel_radiation
    config_max_sr_buffer_size   = max_send_recv_buffer_size
    config_max_mpi_message_size = max_mpi_message_size
    config_use_dycore_barrier   = use_dycore_barrier
    config_use_physics_barrier  = use_physics_barrier
    config_itype_exch_barrier   = itype_exch_barrier
    config_use_sp_output        = use_sp_output
    !-----------------------------------------------------
    CALL check_parallel_configuration()
    
  END SUBROUTINE read_parallel_namelist
  !-------------------------------------------------------------------------

END MODULE mo_parallel_nml
