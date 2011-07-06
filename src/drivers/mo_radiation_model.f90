!>
!! Main program for the ICON radiation model
!!
!!  Used for running the radiation aand the atmo model in parallel
!! 
!! @author
!!     Leonidas Linardakis, MPI-M
!!
!! @date 2011-06-16 
!!
!!
!! @par Copyright
!! 2002-2009 by DWD and MPI-M
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
!-------------------------------------------------------------------------------------
MODULE mo_radiation_model
  
  USE mo_exception,           ONLY: message, finish  ! use always
  !$  USE mo_exception,         ONLY: message_text     ! use only if compiled with OpenMP
  
  USE mo_mpi,                 ONLY: p_start, p_stop, p_pe, p_io, p_nprocs
  USE mo_timer,               ONLY: init_timer, print_timer
  USE mo_namelist,            ONLY: open_nml,  close_nml, open_nml_output, close_nml_output
  USE mo_datetime,            ONLY: t_datetime
  USE mo_output,              ONLY: init_output_files, close_output_files, write_output
  USE mo_io_vlist,            ONLY: write_vlist_oce, destruct_vlist_oce
  
  USE mo_parallel_configuration,        ONLY: & ! process parallel run ctl. params.
    & p_comm_work_test, p_comm_input_bcast, & ! communicators
    & p_test_pe,            & !    internal parameter
    & p_comm_work,          &
    & p_test_run,           &
    & p_io_pe0                ! Number of first I/O PE
  
  USE mo_io_async,            ONLY: io_main_proc            ! main procedure for I/O PEs
  
  
  ! Control parameters: run control, dynamics, i/o
  !
  USE mo_global_variables,    ONLY: setup_physics           ! process forcing control parameters
  !!$     &                              impiom,            & !    :
  USE mo_nonhydrostatic_nml,  ONLY: ivctype,              & ! type of vertical coordinate
    & nonhydrostatic_nml_setup
  USE mo_dynamics_nml,        ONLY: dynamics_nml_setup,   &
    & cleanup_dyn_params
  USE mo_diffusion_nml,       ONLY: diffusion_nml_setup
  USE mo_io_nml,              ONLY: io_nml_setup,         & ! process I/O
    & dt_data,              & !    :
    & dt_file,              & !    :
    & dt_diag,              & !    :
    & lprepare_output         ! internal parameter
!   USE mo_run_nml,             ONLY: run_nml_setup,            & ! process run control parameters
! !     & ini_datetime,         & !    namelist parameter
!     & dtime,                & !    :
!     & i_cell_type,          & !    :
!     & ltransport,           & !    :
!     & lforcing,             & !    :
!     & ltestcase,            & !    :
!     & ltimer,               & !    :
!     & iequations,           & !    internal parameters
!     & ihs_atm_temp,         & !    :
!     & ihs_atm_theta,        & !    :
!     & inh_atmosphere,       & !    :
!     & ishallow_water,       & !    :
!     & ihs_ocean,            & !
!     & iforcing,             & !    namelist parameter
!     & ildf_dry,             & !    :
!     & inoforcing,           & !    internal parameter
!     & iheldsuarez,          & !    :
!     & inwp,                 & !    :
!   !!$     &                              impiom,               & !    :
!     & ldump_states,         & ! flag if states should be dumped
!     & lrestore_states         ! flag if states should be restored
  
  
  
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
  USE mo_gmt_output,          ONLY: setup_gmt_output
  USE mo_nwp_phy_state,       ONLY: construct_nwp_phy_state,   &
    & destruct_nwp_phy_state
  USE mo_atm_phy_nwp_nml,     ONLY: setup_nwp_phy, inwp_surface
  USE mo_lnd_nwp_nml,         ONLY: setup_nwp_lnd
  USE mo_nwp_lnd_state,       ONLY: construct_nwp_lnd_state,   &
    & destruct_nwp_lnd_state
  
  !!$  USE mo_mpiom_phy_state,     ONLY: construct_mpiom_phy_state, &
  !!$    &                               destruct_mpiom_phy_state
  
  USE mo_impl_constants,      ONLY: success
  
  ! External data
  USE mo_ext_data,            ONLY: ext_data, init_ext_data, &
    & destruct_ext_data
  
  !  USE mo_nwp_phy_init,          ONLY: init_nwp_phy
  !!$  USE mo_gscp_cosmo,          ONLY: hydci_pp_init
  
  !-------------------------------------------------------------------------
  IMPLICIT NONE
  
  PRIVATE
  
  PUBLIC :: radiation_model
  !-------------------------------------------------------------------------
  
CONTAINS
  
  !-------------------------------------------------------------------------
  !>
  SUBROUTINE radiation_model(namelist_filename)
    
    CHARACTER(LEN=*), INTENT(in) :: namelist_filename

    ! local variables
    CHARACTER(*), PARAMETER :: method_name = "radiation_model"

    TYPE(t_datetime)        :: datetime    
    
    CALL message(method_name,'start model initialization.')
    
    RETURN
        
  END SUBROUTINE radiation_model
  
END MODULE mo_radiation_model

