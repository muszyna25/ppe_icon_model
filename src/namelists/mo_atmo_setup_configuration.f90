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

  USE mo_parallel_nml,        ONLY: parallel_nml_setup,   & ! process parallel run ctl. params.
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

  ! Memory
  !
  USE mo_subdivision,         ONLY: decompose_atmo_domain,         &
    & copy_processor_splitting,      &
    & set_patch_communicators
  USE mo_dump_restore,        ONLY: dump_patch_state_netcdf,       &
    & restore_patches_netcdf,        &
    & restore_interpol_state_netcdf, &
    & restore_gridref_state_netcdf

  USE mo_icoham_dyn_memory,   ONLY: p_hydro_state
  USE mo_atmo_control,        ONLY: p_patch_global, p_patch_subdiv, p_patch,             &
    & p_nh_state, p_lnd_state,                             &
    & p_int_state_global, p_int_state_subdiv, p_int_state, &
    & p_grf_state_global, p_grf_state_subdiv, p_grf_state

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
  USE mo_intp_state,          ONLY: construct_2d_interpol_state, &
    & destruct_2d_interpol_state
  USE mo_interpolation,       ONLY: rbf_vec_interpol_cell,       &
    & edges2cells_scalar
  USE mo_gridref_nml,         ONLY: gridref_nml_setup
  USE mo_grf_interpolation,   ONLY: construct_2d_gridref_state,  &
    & destruct_2d_gridref_state
  
  ! Vertical grid
  !
  USE mo_vertical_coord_table,ONLY: init_vertical_coord_table
  USE mo_vertical_grid,       ONLY: init_hybrid_coord, init_sleve_coord
  
  ! State variables
  !
  USE mo_icoham_dyn_memory,   ONLY: destruct_icoham_dyn_state
  USE mo_nonhydro_state,      ONLY: destruct_nh_state
  
  
  ! Parameterized forcing
  !
  USE mo_echam_phy_memory,    ONLY: destruct_echam_phy_state
  USE mo_echam_phy_setup,     ONLY: setup_echam_phy
  USE mo_echam_phy_init,      ONLY: prepare_echam_phy, initcond_echam_phy, &
                                  & additional_restart_init
  USE mo_echam_phy_cleanup,   ONLY: cleanup_echam_phy
  USE mo_gmt_output,          ONLY: setup_gmt_output
  USE mo_nwp_phy_state,       ONLY: construct_nwp_phy_state,   &
    & destruct_nwp_phy_state
  USE mo_atm_phy_nwp_nml,     ONLY: setup_nwp_phy, inwp_surface
  USE mo_lnd_nwp_nml,         ONLY: setup_nwp_lnd
  USE mo_nwp_lnd_state,       ONLY: construct_nwp_lnd_state,   &
    & destruct_nwp_lnd_state
 
  USE mo_impl_constants,      ONLY: SUCCESS, MAX_CHAR_LENGTH
  
  ! Time integration
  !
  USE mo_ha_stepping,         ONLY: prepare_ha_dyn, initcond_ha_dyn, perform_ha_stepping
  USE mo_nh_stepping,         ONLY: prepare_nh_integration, perform_nh_stepping
  ! External data
  USE mo_ext_data,            ONLY: ext_data, init_ext_data, destruct_ext_data
  
  !  USE mo_nwp_phy_init,          ONLY: init_nwp_phy
  !!$  USE mo_gscp_cosmo,          ONLY: hydci_pp_init
  
  USE mo_io_restart,           ONLY: read_restart_info_file, read_restart_files
  USE mo_io_restart_namelist,  ONLY: read_restart_namelists
  USE mo_io_restart_attributes,ONLY: read_restart_attributes, get_restart_attribute

  !-------------------------------------------------------------------------
  IMPLICIT NONE
  
  PRIVATE
  
  PUBLIC :: read_atmo_namelists
  
CONTAINS
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
    IF ( iforcing == inwp) THEN
      CALL setup_nwp_phy( p_patch_global(1:) )  ! read Namelist, ...
      IF (inwp_surface > 0) CALL setup_nwp_lnd
    ENDIF
    
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
    
    CALL interpol_nml_setup(p_patch_global)
    
    !
    ! allocate type for interpolation state
    
    !-------------------------------------------------------------------------
    ! set up horizontal diffusion after the vertical coordinate is configured
    !-------------------------------------------------------------------------
    CALL diffusion_nml_setup(n_dom,parent_id)
    
    
    !------------------------------------------------------------------
    ! step 6: initialize output
    !------------------------------------------------------------------
    
    CALL io_nml_setup
 
    
    !------------------------------------------------------------------
    ! Create and optionally read external data fields
    !------------------------------------------------------------------
    ALLOCATE (ext_data(n_dom), STAT=ist)
    IF (ist /= SUCCESS) THEN
      CALL finish(TRIM(method_name),'allocation for ext_data failed')
    ENDIF
        
    CALL close_nml
    IF (p_pe == p_io) THEN
      CALL close_nml_output
    END IF
        
  END SUBROUTINE read_atmo_namelists
  
END MODULE mo_atmo_setup_configuration

