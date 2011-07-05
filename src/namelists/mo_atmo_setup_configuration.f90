!>
!! @brief Main program for the ICON atmospheric model
!!
!! @author
!!  Leonidas Linardakis (MPI-M)
!!  Hui Wan             (MPI-M)
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
MODULE mo_atmo_setup_configuration

  USE mo_exception,           ONLY: message, finish
  USE mo_mpi,                 ONLY: p_stop, p_pe, p_io, p_nprocs
  USE mo_timer,               ONLY: init_timer, print_timer
  USE mo_master_nml,          ONLY: lrestart
  USE mo_namelist,            ONLY: open_nml,  close_nml, open_nml_output, close_nml_output
  USE mo_output,              ONLY: init_output_files, close_output_files, write_output

  USE mo_parallel_nml,        ONLY: parallel_nml_setup,  get_nml_nproma, & ! process parallel run ctl. params.
    & p_comm_work_test, p_comm_input_bcast, & ! communicators
    & p_test_pe,            & !    internal parameter
    & p_comm_work,          &
    & p_test_run,           &
    & p_io_pe0                ! Number of first I/O PE

  USE mo_io_async,            ONLY: io_main_proc            ! main procedure for I/O PEs


  ! Control parameters: run control, dynamics, i/o
  !
  USE mo_global_variables,    ONLY: setup_physics           ! process forcing control parameters
  USE mo_nonhydrostatic_nml,  ONLY: ivctype,              & ! type of vertical coordinate
    & nonhydrostatic_nml_setup
  USE mo_dynamics_nml,        ONLY: dynamics_nml_setup,   &
    &                               cleanup_dyn_params 
  USE mo_diffusion_nml,       ONLY: diffusion_nml_setup
  USE mo_io_nml,              ONLY: io_nml_setup,         & ! process I/O
    & dt_data,              & !    :
    & dt_file,              & !    :
    & dt_diag,              & !    :
    & dt_checkpoint,        & !    :
    & lprepare_output         ! internal parameter
  USE mo_run_nml,             ONLY: run_nml_setup,            & ! process run control parameters
    & current_datetime,     & !    module variable
    & dtime,                & !    namelist parameter
    & nsteps,               & !    :
    & i_cell_type,          & !    :
    & ltransport,           & !    :
    & lforcing,             & !    :
    & ltestcase,            & !    :
    & ltimer,               & !    :
    & iequations,           & !    internal parameters
    & ihs_atm_temp,         & !    :
    & ihs_atm_theta,        & !    :
    & inh_atmosphere,       & !    :
    & ishallow_water,       & !    :
    & iforcing,             & !    namelist parameter
    & ildf_dry,             & !    :
    & ildf_echam,           & !    :
    & inoforcing,           & !    internal parameter
    & iheldsuarez,          & !    :
    & iecham,               & !    :
    & inwp,                 & !    :
    & ldump_states,         & ! flag if states should be dumped
    & lrestore_states         ! flag if states should be restored



  USE mo_advection_nml,       ONLY: transport_nml_setup,  & ! process transport
    & setup_transport         ! control parameters

  ! Test cases
  !
  USE mo_hydro_testcases,     ONLY: setup_testcase          ! process hyd. atm. tests ctl. params.
  USE mo_nh_testcases,        ONLY: setup_nh_testcase       ! process non-hyd. atm. test ctl. par.


  ! Horizontal grid
  !
  USE mo_model_domain_import, ONLY: grid_nml_setup,          & ! process grid control parameters
    & n_dom,                & !    :
    & n_dom_start,          & !    :
    & parent_id,            & !    :
    & import_patches,       & !
    & destruct_patches        !

  ! Horizontal interpolation
  !
  USE mo_interpol_nml,        ONLY: interpol_nml_setup   ! process interpol. ctl. params.
  USE mo_gridref_nml,         ONLY: gridref_nml_setup
  
  ! Vertical grid
  !
  USE mo_vertical_coord_table,ONLY: init_vertical_coord_table
  USE mo_vertical_grid,       ONLY: init_hybrid_coord, init_sleve_coord
    
  USE mo_echam_phy_setup,     ONLY: setup_echam_phy
  USE mo_atm_phy_nwp_nml,     ONLY: setup_nwp_phy, inwp_surface
  USE mo_lnd_nwp_nml,         ONLY: setup_nwp_lnd
 
  USE mo_impl_constants,      ONLY: SUCCESS, MAX_CHAR_LENGTH
    
  USE mo_io_restart_namelist,  ONLY: read_restart_namelists
  USE mo_io_restart_attributes,ONLY: read_restart_attributes, get_restart_attribute

  !-------------------------------------------------------------------------
  ! use configure modules
  USE mo_parallel_configuration, ONLY: set_nproma

  
  IMPLICIT NONE
  
  PRIVATE
  
  PUBLIC :: read_atmo_namelists, setup_atmo_configuration
  
CONTAINS
  
  !-------------------------------------------------------------------------
  !>
  !!
  SUBROUTINE read_atmo_namelists(namelist_filename)
    
    CHARACTER(LEN=*), INTENT(in) :: namelist_filename

    CHARACTER(*), PARAMETER :: method_name = "read_atmo_namelists"

    CHARACTER(LEN=MAX_CHAR_LENGTH) :: grid_file_name 
    INTEGER :: n_io, jg, jfile, n_file, ist, n_diag, n_chkpt
    LOGICAL :: lsuccess, l_have_output
   

    !-------------------------------------------------------------------
    ! 1. Open the atmosphere-specific namelist file and create a 
    !    new file in which all the namelist variables and their
    !    actual values used in the model run will be stored.
    !-------------------------------------------------------------------
    CALL open_nml(TRIM(namelist_filename))
    IF(p_pe == p_io) CALL open_nml_output('NAMELIST_ICON_output_atm')

    ! read namelsist
    CALL read_parallel_namelist()
        
    ! close namelist file
    CALL close_nml
    IF (p_pe == p_io) THEN
      CALL close_nml_output
    END IF
        
  END SUBROUTINE read_atmo_namelists
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  !>
  !!
  SUBROUTINE old_read_atmo_namelists(namelist_filename)
    
    CHARACTER(LEN=*), INTENT(in) :: namelist_filename

    CHARACTER(*), PARAMETER :: method_name = "read_atmo_namelists"

    CHARACTER(LEN=MAX_CHAR_LENGTH) :: grid_file_name 
    INTEGER :: n_io, jg, jfile, n_file, ist, n_diag, n_chkpt
    LOGICAL :: lsuccess, l_have_output
   

    !-------------------------------------------------------------------
    ! 1. Open the atmosphere-specific namelist file and create a 
    !    new file in which all the namelist variables and their
    !    actual values used in the model run will be stored.
    !-------------------------------------------------------------------
    CALL open_nml(TRIM(namelist_filename))
    IF(p_pe == p_io) CALL open_nml_output('NAMELIST_ICON_output_atm')

    ! The namelists ('run_nml' and 'testcase_ctl') are read in seperate
    ! subroutines. The validity of the user-specified configurations is
    ! checked right after each namelist is read.
    ! The two 'setup' subroutines above must be called before grid and patch
    ! import because during the import procedure the topography and land-sea
    ! mask will be initialized. They are related to the atmos/ocean switch
    ! as well as the selected test case.
    
    CALL run_nml_setup
    
    !-------------------------------------------------------------------
    ! parallel_nml_setup must be called after setup_run since it needs
    ! some variables read in setup_run
    
    CALL parallel_nml_setup
    
    !-------------------------------------------------------------------
    ! Initialize test case setup 
    !-------------------------------------------------------------------
    IF (ltestcase) THEN
      
      SELECT CASE (iequations)
      CASE (ishallow_water, ihs_atm_temp, ihs_atm_theta)
        CALL setup_testcase
        
      CASE (inh_atmosphere)
        CALL setup_nh_testcase
        
      CASE DEFAULT
        CALL finish( TRIM(method_name),' invalid value for iequations!' )
      END SELECT
      
    ENDIF
    
    !-------------------------------------------------------------------
    ! step 2: import computational domain of the model
    !-------------------------------------------------------------------
    ! 2a) read namelist 'grid_ctl' from the namelist file (already
    !     opened above) and set up the grid/patch configuration.
    !-------------------------------------------------------------------
    
    CALL grid_nml_setup
        
    !------------------------------------------------------------------
    ! step 3a: Read namelist 'dynamics_ctl'. This has to be done
    !          before the interpolation state is computed.
    !------------------------------------------------------------------
    CALL dynamics_nml_setup(n_dom)
    
    !------------------------------------------------------------------
    ! Read specific dynamics namelists for the nonhydrost. dynamical core
    !------------------------------------------------------------------
    CALL nonhydrostatic_nml_setup
    
    !------------------------------------------------------------------
    ! step 3a-a: ! read nwp physics namelist, ...
    !------------------------------------------------------------------
!     IF ( iforcing == inwp) THEN
!       CALL setup_nwp_phy( p_patch_global(1:) )  ! read Namelist, ...
!       IF (inwp_surface > 0) CALL setup_nwp_lnd
!     ENDIF
    
    !------------------------------------------------------------------
    ! step 3b: Read namelist 'transport_ctl',
    !------------------------------------------------------------------
    CALL transport_nml_setup
    
    !------------------------------------------------------------------
    ! step 4: construct interpolation and the grid refinement state
    !------------------------------------------------------------------
    ! - Allocate memory for the interpolation state;
    ! - Calculate interpolation coefficients.
    !------------------------------------------------------------------
    
    CALL gridref_nml_setup
    
    ! interpolation state not used for ocean model
    ! #slo# - temporarily switched on for comparison with rbf-reconstruction
    
!     CALL interpol_nml_setup(p_patch_global)
    
    !-------------------------------------------------------------------------
    ! set up horizontal diffusion after the vertical coordinate is configured
    !-------------------------------------------------------------------------
!     CALL diffusion_nml_setup(n_dom,parent_id)
    
   
    !------------------------------------------------------------------
    ! step 6: initialize output
    !------------------------------------------------------------------    
    CALL io_nml_setup
     
    !------------------------------------------------------------------
    ! Create and optionally read external data fields
    !------------------------------------------------------------------
        
    CALL close_nml
    IF (p_pe == p_io) THEN
      CALL close_nml_output
    END IF
        
  END SUBROUTINE old_read_atmo_namelists
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE setup_atmo_configuration()
  
    CALL setup_parallel_configuration()

  END SUBROUTINE setup_atmo_configuration
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  SUBROUTINE setup_parallel_configuration()
  

  !------------------------------------------------------------
  ! 3.0 check the consistency of the parameters
  !------------------------------------------------------------
  IF (nproma<=0) CALL finish(TRIM(method_name),'"nproma" must be positive')

  ! check n_ghost_rows
  IF (n_ghost_rows<1) THEN
    CALL finish(method_name, &
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
    CALL finish(method_name, &
       & 'division_method=div_metis=2 in parallel_ctl namelist is not allowed')
#endif
  CASE DEFAULT
    CALL finish(method_name, &
       & 'value of division_method in parallel_ctl namelist is not allowed')
  END SELECT

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
  ! A run on 1 PE is never a verification run,
  ! correct this if the user should set it differently
  IF (p_test_run .AND. p_nprocs == 1) THEN
    CALL message(method_name, &
         & 'p_test_run has no effect if p_nprocs=1')
    CALL message(method_name, &
         & '--> p_test_run set to .FALSE.')
     p_test_run = .FALSE.
  ENDIF
  ! for safety only
  IF(num_io_procs < 0) num_io_procs = 0
#endif

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
    CALL finish(method_name, &
       & 'not enough processors for given values of p_test_run/num_io_procs')
  ELSE IF (p_test_run .AND. num_work_procs == 1) THEN
    CALL finish(method_name, &
       & 'running p_test_run with only 1 work processor does not make sense')
  ENDIF

  WRITE(message_text,'(3(a,i0))') 'Number of procs for test: ',num_test_procs, &
    & ', work: ',num_work_procs,', I/O: ',num_io_procs

  CALL message(method_name, message_text)

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
    CALL set_nproma(get_nml_nproma())

  END SUBROUTINE setup_parallel_configuration
  !-------------------------------------------------------------------------
  
  
END MODULE mo_atmo_setup_configuration

