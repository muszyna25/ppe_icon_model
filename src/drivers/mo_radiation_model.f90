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
MODULE mo_radiation_model

USE mo_exception,           ONLY: message, finish
USE mo_mpi,                 ONLY: p_stop, &
  & my_process_is_io,  my_process_is_mpi_seq, my_process_is_mpi_test, &
  & set_mpi_work_communicators, set_comm_input_bcast, null_comm_type
USE mo_timer,               ONLY: init_timer
USE mo_parallel_config,     ONLY: p_test_run, l_test_openmp, num_io_procs
USE mo_master_control,      ONLY: is_restart_run, get_my_couple_id


USE mo_io_async,            ONLY: io_main_proc            ! main procedure for I/O PEs

! Control parameters: run control, dynamics, i/o
!
USE mo_nonhydrostatic_config,ONLY: ivctype, kstart_moist, kstart_qv,    &
  &                                l_open_ubc, configure_nonhydrostatic
USE mo_dynamics_config,   ONLY: iequations
USE mo_run_config,        ONLY: configure_run, &
  & ltimer,               & !    :
  & iforcing,             & !    namelist parameter
  & nlev,nlevp1,          &
  & num_lev,num_levp1,    &
  & lvert_nest, ntracer,  &
  & nshift
  
USE mo_impl_constants, ONLY:&
  & ihs_atm_temp,         & !    :
  & ihs_atm_theta,        & !    :
  & inh_atmosphere,       & !    :
  & ishallow_water



! For the coupling
USE mo_impl_constants, ONLY: CELLS
USE mo_master_control, ONLY : is_coupled_run
USE mo_icon_cpl_init_comp, ONLY : get_my_local_comp_id
USE mo_icon_cpl_def_grid, ONLY : ICON_cpl_def_grid
USE mo_icon_cpl_def_field, ONLY : ICON_cpl_def_field
USE mo_icon_cpl_search, ONLY : ICON_cpl_search
USE mo_model_domain_import, ONLY : get_patch_global_indexes

! Memory
!
USE mo_subdivision,         ONLY: decompose_atmo_domain,         &
& copy_processor_splitting,      &
& set_patch_communicators

USE mo_atmo_control,        ONLY: p_patch_global, p_patch_subdiv, p_patch

! Horizontal grid
USE mo_grid_config,         ONLY: n_dom, n_dom_start
USE mo_model_domain_import, ONLY: import_patches, destruct_patches

! Horizontal interpolation
!
! USE mo_intp_state,          ONLY: construct_2d_interpol_state,  destruct_2d_interpol_state
! USE mo_grf_interpolation,   ONLY: construct_2d_gridref_state,   destruct_2d_gridref_state

! Vertical grid
!
USE mo_vertical_coord_table,ONLY: init_vertical_coord_table, &
  &                               vct_a, vct_b, apzero
USE mo_nh_init_utils,       ONLY: init_hybrid_coord, init_sleve_coord

! Parameterized forcing
!
USE mo_impl_constants,      ONLY: SUCCESS, MAX_CHAR_LENGTH

! Time integration
!
! External data
! USE mo_ext_data,            ONLY: ext_data, init_ext_data, destruct_ext_data

!  USE mo_nwp_phy_init,          ONLY: init_nwp_phy
!!$  USE mo_gscp_cosmo,          ONLY: hydci_pp_init


!-------------------------------------------------------------------------
! to break circular dependency
! 
! USE mo_intp_data_strc,      ONLY: p_int_state_global, p_int_state_subdiv, p_int_state
! USE mo_grf_intp_data_strc,  ONLY: p_grf_state_global, p_grf_state_subdiv, p_grf_state

!-------------------------------------------------------------------------

USE mo_read_namelists,     ONLY: read_atmo_namelists
USE mo_atm_nml_crosscheck, ONLY: atm_crosscheck

! USE mo_dynamics_config,    ONLY: configure_dynamics  ! subroutine
!USE mo_interpol_config,    ONLY: configure_interpolation 
! USE mo_interpol_config
! USE mo_advection_config,   ONLY: configure_advection
! USE mo_diffusion_config,   ONLY: configure_diffusion
! 
! USE mo_atmo_hydrostatic,    ONLY: atmo_hydrostatic 
! USE mo_atmo_nonhydrostatic, ONLY: atmo_nonhydrostatic 

!-------------------------------------------------------------------------
IMPLICIT NONE
PRIVATE

PUBLIC :: radiation_model

CONTAINS

  !-------------------------------------------------------------------------
  !>
  !!
  SUBROUTINE radiation_model(rad_namelist_filename,shr_namelist_filename)

    CHARACTER(LEN=*), INTENT(in) :: rad_namelist_filename
    CHARACTER(LEN=*), INTENT(in) :: shr_namelist_filename

    CHARACTER(LEN=MAX_CHAR_LENGTH) :: grid_file_name 
    CHARACTER(*), PARAMETER :: routine = "mo_radiation_model:atmo_model"
    LOGICAL :: lsuccess
    INTEGER :: jg

    ! For the coupling

    INTEGER, PARAMETER :: no_of_fields = 12

    CHARACTER(LEN=MAX_CHAR_LENGTH) ::  field_name(no_of_fields)
    INTEGER :: field_id(no_of_fields)
    INTEGER :: comp_id
    INTEGER :: grid_id
    INTEGER :: grid_shape(2) 
    INTEGER :: field_shape(3) 
    INTEGER :: i, error_status
    INTEGER :: no_of_entities
    INTEGER :: patch_no
    INTEGER, POINTER :: grid_glob_index(:)


    !---------------------------------------------------------------------
    ! 1.1 Read namelists (newly) specified by the user; fill the 
    !     corresponding sections of the configuration states.
    !---------------------------------------------------------------------

    CALL read_atmo_namelists(rad_namelist_filename,shr_namelist_filename)

    !---------------------------------------------------------------------
    ! 1.2 Cross-check namelist setups
    !---------------------------------------------------------------------

    CALL atm_crosscheck

    !---------------------------------------------------------------------
    ! 2. Call configure_run to finish filling the run_config state.
    !    This needs to be done very early (but anyway after atm_crosscheck)
    !    because some component of the state, e.g., num_lev, may be 
    !    modified in this subroutine which affect the following CALLs.
    !---------------------------------------------------------------------
    CALL configure_run

    !-------------------------------------------------------------------
    ! 3.1 Initialize the mpi work groups
    !-------------------------------------------------------------------
    CALL set_mpi_work_communicators(p_test_run, l_test_openmp, num_io_procs)

    !-------------------------------------------------------------------
    ! 3.2 Initialize various timers
    !-------------------------------------------------------------------
    IF (ltimer) CALL init_timer

    !------------------
    ! Next, define the horizontal and vertical grids since they are aready
    ! needed for some derived control parameters. This includes
    ! - patch import
    ! - domain decompistion
    ! - vertical coordinates
    !-------------------------------------------------------------------
    ! 4. Import patches
    !-------------------------------------------------------------------
    
    ! Check patch allocation status

    IF ( ALLOCATED(p_patch_global)) THEN
      CALL finish(TRIM(routine), 'patch already allocated')
    END IF
     
    ! Allocate patch array to start patch construction

    ALLOCATE(p_patch_global(n_dom_start:n_dom), stat=error_status)
    IF (error_status/=success) THEN
      CALL finish(TRIM(routine), 'allocation of patch failed')
    ENDIF
    
    CALL import_patches( p_patch_global,                       &
      & nlev,nlevp1,num_lev,num_levp1,nshift)

    !--------------------------------------------------------------------------------
    ! 5. Construct interpolation state, compute interpolation coefficients.
    !--------------------------------------------------------------------------------
    

    !-----------------------------------------------------------------------------
    ! 6. Construct grid refinment state, compute coefficients
    !-----------------------------------------------------------------------------
    ! For the NH model, the initialization routines called from
    ! construct_2d_gridref_state require the metric terms to be present


    !-------------------------------------------------------------------
    ! 7. Domain decomposition: 
    !    Divide patches and interpolation states for parallel runs.
    !    This is only done if the model runs really in parallel.
    !-------------------------------------------------------------------   
    IF (my_process_is_mpi_seq()) THEN
      
      ! This is a run on a single processor
      ! just set pointers
      
      p_patch => p_patch_global
      
      CALL set_patch_communicators(p_patch)
      
    ELSE
      
      CALL decompose_atmo_domain()
      
    ENDIF

    ! Note: from this point the p_patch is used
    ! In case of a test run: Copy processor splitting to test PE
    IF(p_test_run) CALL copy_processor_splitting(p_patch)
    
   !---------------------------------------------------------------------
   ! 8. Import vertical grid/ define vertical coordinate
   !---------------------------------------------------------------------
    SELECT CASE (iequations)
    
    CASE (ishallow_water)
      CALL init_vertical_coord_table(iequations, p_patch(1)%nlev)
      
    CASE (ihs_atm_temp, ihs_atm_theta)
      CALL init_vertical_coord_table(iequations, p_patch(1)%nlev)
      
    CASE (inh_atmosphere)
      IF (ivctype == 1) THEN
        CALL init_hybrid_coord(p_patch(1)%nlev)
      ELSE IF (ivctype == 2) THEN
        CALL init_sleve_coord(p_patch(1)%nlev)
      ENDIF

    CASE DEFAULT
    END SELECT


    !---------------------------------------------------------------------
    ! 9. Horizontal and vertical grid(s) are now defined.
    !    Assign values to derived variables in the configuration states
    !---------------------------------------------------------------------


    !------------------------------------------------------------------
    ! 10. Create and optionally read external data fields
    !------------------------------------------------------------------
!     ALLOCATE (ext_data(n_dom), STAT=error_status)
!     IF (error_status /= SUCCESS) THEN
!       CALL finish(TRIM(routine),'allocation for ext_data failed')
!     ENDIF
    
    ! allocate memory for atmospheric/oceanic external data and
    ! optionally read those data from netCDF file.
!     CALL init_ext_data (p_patch(1:), p_int_state(1:), ext_data)

    !---------------------------------------------------------------------
    ! 11. Do the setup for the coupled run
    !
    ! For the time being this could all go into a subroutine which is
    ! common to atmo and ocean. Does this make sense if the setup deviates
    ! too much in future.
    !---------------------------------------------------------------------
    IF ( is_coupled_run() ) THEN
 
      comp_id = get_my_couple_id ()
      patch_no = 1

      grid_shape(1)=1
      grid_shape(2)=p_patch(patch_no)%n_patch_cells

      ! CALL get_patch_global_indexes ( patch_no, CELLS, no_of_entities, grid_glob_index )
      ! should grid_glob_index become a pointer in ICON_cpl_def_grid as well?
      CALL ICON_cpl_def_grid ( &
        & comp_id, grid_shape, p_patch(patch_no)%cells%glb_index, & ! input
        & grid_id, error_status )                                   ! output
  
      field_name(1) = "SST"
      field_name(2) = "TAUX"
      field_name(3) = "TAUY"
      field_name(4) = ""
      field_name(5) = ""
      field_name(6) = ""
      field_name(7) = ""
      field_name(8) = ""
 
!       DO i = 1, no_of_fields
!          CALL ICON_cpl_def_field ( field_name(i), comp_id, grid_id, field_id(i), &
!                                  & field_shape, error_status )
!       ENDDO

      CALL ICON_cpl_search

    ENDIF

    !---------------------------------------------------------------------
    ! 12. The hydrostatic and nonhydrostatic models branch from this point
    !---------------------------------------------------------------------
!     SELECT CASE(iequations)
!     CASE(ishallow_water,ihs_atm_temp,ihs_atm_theta)
!       CALL atmo_hydrostatic
! 
!     CASE(inh_atmosphere)
!       CALL atmo_nonhydrostatic
! 
!     CASE DEFAULT
!       CALL finish( TRIM(routine),'unknown choice for iequaions.')
!     END SELECT


    !---------------------------------------------------------------------
    ! 13. Integration finished. Carry out the shared clean-up processes
    !---------------------------------------------------------------------
    ! Destruct external data state

!     CALL destruct_ext_data

    ! deallocate ext_data array
!     DEALLOCATE(ext_data, stat=error_status)
!     IF (error_status/=success) THEN
!       CALL finish(TRIM(routine), 'deallocation of ext_data')
!     ENDIF

    ! Deconstruct grid refinement state



    ! Deallocate grid patches

    CALL destruct_patches( p_patch )

    IF (my_process_is_mpi_seq()) THEN
      DEALLOCATE( p_patch_global, STAT=error_status )
    ELSE
      DEALLOCATE( p_patch_subdiv, STAT=error_status )
    ENDIF
    IF (error_status/=SUCCESS) THEN
      CALL finish(TRIM(routine),'deallocate for patch array failed')
    ENDIF

    CALL message(TRIM(routine),'clean-up finished')

  END SUBROUTINE radiation_model
  !-------------------------------------------------------------------------

END MODULE mo_radiation_model

