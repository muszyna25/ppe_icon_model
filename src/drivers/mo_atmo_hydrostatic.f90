!>
!! @brief workflow for the ICON atmospheric hydrostatic model
!!
!! @author
!!  Hui Wan             (MPI-M)
!!  Kristina Froehlich, MPI-M
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
MODULE mo_atmo_hydrostatic

USE mo_exception,           ONLY: message, finish
USE mo_mpi,                 ONLY: p_stop, p_pe, p_io,  &
& my_process_is_io,  my_process_is_mpi_seq, my_process_is_mpi_test, &
& my_process_is_stdio
USE mo_timer,               ONLY: init_timer, print_timer
USE mo_master_nml,          ONLY: lrestart
USE mo_namelist,            ONLY: open_nml,  close_nml, open_nml_output, close_nml_output
USE mo_output,              ONLY: init_output_files, close_output_files, write_output

USE mo_parallel_config, ONLY: p_test_run

USE mo_io_async,            ONLY: io_main_proc            ! main procedure for I/O PEs

! Control parameters: run control, dynamics, i/o
!
USE mo_global_variables,    ONLY: setup_physics           ! process forcing control parameters
USE mo_nonhydrostatic_config,ONLY: ivctype, kstart_moist, kstart_qv, l_open_ubc
!USE mo_io_config,           ONLY: dt_data, dt_file, dt_diag, dt_checkpoint 
USE mo_io_config,         ONLY:  dt_data,dt_file,dt_diag,dt_checkpoint
USE mo_dynamics_config,   ONLY: iequations
USE mo_run_config,        ONLY: &
& dtime,                & !    namelist parameter
& nsteps,               & !    :
& ltransport,           & !    :
& lforcing,             & !    :
& ltestcase,            & !    :
& ltimer,               & !    :
& iforcing,             & !    namelist parameter
& ldump_states,         & ! flag if states should be dumped
& lrestore_states,      & ! flag if states should be restored
& nlev,nlevp1,          &
& num_lev,num_levp1,    &
& iqv, nshift,          &
& lvert_nest, ntracer

USE mo_impl_constants, ONLY:&
    & ihs_atm_temp,         & !    :
    & ihs_atm_theta,        & !    :
    & inh_atmosphere,       & !    :
    & ishallow_water,       & !    :
    & ildf_dry,             & !    :
    & ildf_echam,           & !    :
    & inoforcing,           & !    :
    & iheldsuarez,          & !    :
    & iecham,               & !    :
    & inwp

USE mo_advection_nml,       ONLY: transport_nml_setup,  & ! process transport
& setup_transport         ! control parameters

! For the coupling
USE mo_impl_constants, ONLY: CELLS
USE mo_master_control, ONLY : atmo_process, is_coupled_run
USE mo_icon_cpl_init_comp, ONLY : get_my_local_comp_id
USE mo_icon_cpl_def_grid, ONLY : ICON_cpl_def_grid
USE mo_icon_cpl_def_field, ONLY : ICON_cpl_def_field
USE mo_icon_cpl_search, ONLY : ICON_cpl_search
USE mo_model_domain_import, ONLY : get_patch_global_indexes

! Test cases
!
USE mo_ha_testcases,     ONLY: read_ha_testcase_namelist, ctest_name
USE mo_nh_testcases,     ONLY: read_nh_testcase_namelist ! process non-hyd. atm. test ctl. par.

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
USE mo_nonhydro_state,      ONLY: p_nh_state
USE mo_atmo_control,        ONLY: p_patch_global, p_patch_subdiv, p_patch

! Horizontal grid
!
USE mo_grid_nml,            ONLY: read_grid_namelist
USE mo_model_domain_import, ONLY: &  !grid_nml_setup,          & ! process grid control parameters
& n_dom,                & !    :
& n_dom_start,          & !    :
!    & parent_id,            & !    :
& import_patches,       & !
& destruct_patches        !

! Horizontal interpolation
!
USE mo_intp_state,          ONLY: construct_2d_interpol_state, &
& destruct_2d_interpol_state
USE mo_interpolation,       ONLY: rbf_vec_interpol_cell,       &
& edges2cells_scalar
USE mo_grf_interpolation,   ONLY: construct_2d_gridref_state,  &
& destruct_2d_gridref_state

! Vertical grid
!
USE mo_vertical_coord_table,ONLY: init_vertical_coord_table, &
  &                               vct_a, vct_b, ceta, apzero
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
USE mo_atm_phy_nwp_config, ONLY: atm_phy_nwp_config,configure_atm_phy_nwp
USE mo_lnd_nwp_nml,         ONLY: setup_nwp_lnd
USE mo_nwp_lnd_state,       ONLY: construct_nwp_lnd_state,   &
& destruct_nwp_lnd_state, p_lnd_state

USE mo_impl_constants,      ONLY: SUCCESS, MAX_CHAR_LENGTH

! Time integration
!
USE mo_ha_stepping,         ONLY: prepare_ha_dyn, initcond_ha_dyn, perform_ha_stepping
USE mo_nh_stepping,         ONLY: prepare_nh_integration, perform_nh_stepping
! External data
USE mo_ext_data,            ONLY: ext_data, init_ext_data, destruct_ext_data

!  USE mo_nwp_phy_init,          ONLY: init_nwp_phy
!!$  USE mo_gscp_cosmo,          ONLY: hydci_pp_init


!-------------------------------------------------------------------------
! to break circular dependency

USE mo_intp_data_strc,      ONLY: p_int_state_global, p_int_state_subdiv, p_int_state
USE mo_grf_intp_data_strc,  ONLY: p_grf_state_global, p_grf_state_subdiv, p_grf_state

!-------------------------------------------------------------------------
USE mo_io_restart,           ONLY: read_restart_info_file, read_restart_files
USE mo_io_restart_namelist,  ONLY: read_restart_namelists
USE mo_io_restart_attributes,ONLY: read_restart_attributes, get_restart_attribute

USE mo_atmo_setup_configuration, ONLY: read_atmo_namelists
USE mo_atm_nml_crosscheck,       ONLY: atm_crosscheck

USE mo_time_config,        ONLY: time_config      ! variable
USE mo_dynamics_config,    ONLY: configure_dynamics  ! subroutine
USE mo_interpol_config,    ONLY: configure_interpolation 
USE mo_advection_config,   ONLY: configure_advection
USE mo_diffusion_config,   ONLY: configure_diffusion
USE mo_echam_phy_config,   ONLY: configure_echam_phy 
USE mo_echam_conv_config,  ONLY: configure_echam_convection

!-------------------------------------------------------------------------
IMPLICIT NONE
PRIVATE

PUBLIC :: atmo_hydrostatic

CONTAINS
!>
!!
  SUBROUTINE atmo_hydrostatic

    CHARACTER(*), PARAMETER :: routine = "mo_atmo_hydrostatic"
    LOGICAL :: lsuccess
    LOGICAL :: l_have_output
    INTEGER, POINTER :: grid_glob_index(:)

   
   !CALL configure_echam_phy (ltestcase, ctest_name)
   !CALL configure_echam_convection(nlev, vct_a, vct_b, ceta)


!    !------------------------------------------------------------------
!    ! Create and optionally read external data fields
!    !------------------------------------------------------------------
!    ALLOCATE (ext_data(n_dom), STAT=ist)
!    IF (ist /= SUCCESS) THEN
!      CALL finish(TRIM(routine),'allocation for ext_data failed')
!    ENDIF
!    
!    ! allocate memory for atmospheric/oceanic external data and
!    ! optionally read those data from netCDF file.
!    CALL init_ext_data (p_patch(1:), p_int_state, ext_data)
!

   
!      ALLOCATE (p_hydro_state(n_dom), stat=ist)
!      IF (ist /= success) THEN
!        CALL finish(TRIM(routine),'allocation for p_hydro_state failed')
!      ENDIF

!        IF(iforcing=iecham .OR. iforcing = ildf_echam)THEN
!         CALL setup_echam_phy
!        ENDIF


    !---------------------------------------------------------------------
    !  Perform time stepping
    !---------------------------------------------------------------------
    ! Initial conditions
    !------------------------------------------------------------------
    ! Prepare for time integration
    !------------------------------------------------------------------

!   !------------------------------------------------------------------
!   ! Initialize parameters and solvers;
!   ! Allocate memory for model state vectors.
!   !------------------------------------------------------------------
!   CALL prepare_ha_dyn( p_patch(1:) )
!   IF (iforcing==IECHAM.OR.iforcing==ILDF_ECHAM) THEN
!     CALL prepare_echam_phy( p_patch(1:) )
!   END IF
!
!   !------------------------------------------------------------------
!   ! Set initial conditions for time integration.
!   !------------------------------------------------------------------
!   IF (lrestart) THEN
!   ! This is an resumed integration. Read model state from restart file(s).
!
!     CALL read_restart_files
!     CALL message(TRIM(routine),'normal exit from read_restart_files')
!
!     ! Initialize logical variables in echam physics state.
!     ! This is necessary for now because logical arrays can not yet be
!     ! written into restart files.
!
!     IF (iforcing==IECHAM.OR.iforcing==ILDF_ECHAM) THEN                                       
!       CALL additional_restart_init( p_patch(1:) )                                            
!     END IF                                                                                   
!
!   ELSE
!   ! This is an initial run (cold start). Compute initial condition for 
!   ! test cases, or read externally given initial conditions.
!
!     CALL initcond_ha_dyn( p_patch(1:), p_int_state(1:),  &
!                         & p_grf_state(1:), p_hydro_state )
!
!     IF (iforcing==IECHAM.OR.iforcing==ILDF_ECHAM)      &
!     CALL initcond_echam_phy( p_patch(1:),p_hydro_state )
!
!   END IF ! lrestart
!                                                                                              

!    !---------------------------------------------------------
!    ! The most primitive event handling algorithm: 
!    ! compute time step interval for taking a certain action
!    !--------------------------------------------------------- 
! 
!    n_io    = NINT(dt_data/dtime)        ! write output
!    n_file  = NINT(dt_file/dtime)        ! trigger new output file
!    n_chkpt = NINT(dt_checkpoint/dtime)  ! write restart files
!    n_diag  = MAX(1,NINT(dt_diag/dtime)) ! diagnose of total integrals
!
!    !------------------------------------------------------------------
!    ! Prepare output file
!    !------------------------------------------------------------------
!    IF (.NOT.lrestart) THEN
!    ! Initialize the first output file which will contain also the 
!    ! initial conditions.
!
!      jfile = 1
!      CALL init_output_files(jfile, lclose=.FALSE.)
!
!    ELSE
!    ! No need to write out the initial condition, thus no output
!    ! during the first integration step. This run will produce
!    ! output if n_io <= integration_length. 
!
!      CALL get_restart_attribute('next_output_file',jfile)
!
!      IF (n_io.le.(nsteps-1)) THEN
!         CALL init_output_files(jfile, lclose=.FALSE.)
!         l_have_output = .TRUE.  
!      ELSE
!         l_have_output = .FALSE.
!      END IF
!
!    END IF
! 

!    !------------------------------------------------------------------
!    !  get and write out some of the inital values
!    !------------------------------------------------------------------
!    IF (.NOT.lrestart) THEN
!
!    ! diagnose u and v to have meaningful initial output
!    
!    DO jg = 1, n_dom

!        SELECT CASE (p_patch(jg)%cell_type)
!        CASE (3)
!          CALL rbf_vec_interpol_cell(p_hydro_state(jg)%prog(1)%vn,p_patch(jg), &
!            & p_int_state(jg),p_hydro_state(jg)%diag%u,p_hydro_state(jg)%diag%v)
!        CASE (6)
!          CALL edges2cells_scalar(p_hydro_state(jg)%prog(1)%vn,p_patch(jg), &
!            & p_int_state(jg)%hex_east,p_hydro_state(jg)%diag%u)
!          CALL edges2cells_scalar(p_hydro_state(jg)%prog(1)%vn,p_patch(jg), &
!            & p_int_state(jg)%hex_north,p_hydro_state(jg)%diag%v)
!        END SELECT
!
!    ENDDO
!    
!    ! Note: here the derived output variables are not yet available
!    ! (omega, divergence, vorticity)
!    CALL write_output( time_config%cur_datetime )
!    l_have_output = .TRUE.
!
!    END IF ! not lrestart
!

!    !------------------------------------------------------------------
!    ! Now start the time stepping:
!    ! The special initial time step for the three time level schemes
!    ! is executed within process_grid_level
!    !------------------------------------------------------------------

!      CALL perform_ha_stepping( p_patch(1:), p_int_state(1:), p_grf_state(1:), &
!                              & p_hydro_state, time_config%cur_datetime,       &
!                              & n_io, n_file, n_chkpt, n_diag, jfile,          &
!                              & l_have_output                                  )
!
!    IF (ltimer) CALL print_timer
!
    ! 
    !---------------------------------------------------------------------
    ! 6. Integration finished. Clean up.
    !---------------------------------------------------------------------

!   CALL message(TRIM(routine),'start to clean up')
!   
!   ! Delete state variables
!

!     CALL destruct_icoham_dyn_state
!     DEALLOCATE (p_hydro_state, STAT=ist)
!     IF (ist /= SUCCESS) THEN
!       CALL finish(TRIM(routine),'deallocation for p_hydro_state failed')
!     ENDIF
!
!   ! Delete output variable lists
!   IF (l_have_output) CALL close_output_files
!   
!
!   IF(iforcing == iecham .OR. iforcing== ildf_echam) THEN
!     CALL destruct_echam_phy_state  ! deallocate state vector
!     CALL cleanup_echam_phy         ! deallocate parameter arrays
!   ENDIF

!   ! Deallocate external data
!   CALL destruct_ext_data
!
!   ! deallocate ext_data array
!   DEALLOCATE(ext_data, stat=ist)
!   IF (ist/=success) THEN
!     CALL finish(TRIM(routine), 'deallocation of ext_data')
!   ENDIF
!   
!   CALL message(TRIM(routine),'clean-up finished')
!   

  END SUBROUTINE atmo_hydrostatic
  
END MODULE mo_atmo_hydrostatic

