!>
!! @brief Main program for the prep_icon preprocessor
!!
!! @author
!!  Guenther Zaengl (DWD)
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
PROGRAM prep_icon

USE mo_exception,           ONLY: message, finish
!$  USE mo_exception,         ONLY: message_text     ! use only if compiled with OpenMP
USE mo_parallel_config,     ONLY: p_test_run, l_test_openmp, num_io_procs
USE mo_mpi,                 ONLY: start_mpi, p_stop, &
  & my_process_is_io, my_process_is_mpi_seq, set_mpi_work_communicators
USE mo_timer,               ONLY: init_timer


! Control parameters: run control, dynamics, i/o
!
USE mo_nonhydrostatic_config,ONLY: ivctype
USE mo_dynamics_config,   ONLY: iequations
USE mo_run_config,        ONLY: configure_run, &
  & ltimer,               & !    :
  & nlev,nlevp1,          &
  & num_lev,num_levp1, nshift

USE mo_impl_constants, ONLY:  inh_atmosphere

! Memory
!
USE mo_subdivision,         ONLY: decompose_atmo_domain,         &
& copy_processor_splitting,  set_patch_communicators


USE mo_atmo_control,        ONLY: p_patch_global, p_patch_subdiv, p_patch

! Horizontal grid
USE mo_grid_config,         ONLY: n_dom, n_dom_start, global_cell_type
                                  
USE mo_model_domain_import, ONLY: import_patches, destruct_patches

! Horizontal interpolation
!
USE mo_intp_state,          ONLY: construct_2d_interpol_state, &
& destruct_2d_interpol_state

USE mo_grf_interpolation,   ONLY: construct_2d_gridref_state,  &
& destruct_2d_gridref_state

! Vertical grid
!

USE mo_nh_init_utils,       ONLY: init_hybrid_coord, init_sleve_coord

! State variables
!


USE mo_impl_constants,      ONLY: SUCCESS, MAX_CHAR_LENGTH


USE mo_intp_data_strc,      ONLY: p_int_state_global, p_int_state_subdiv, p_int_state
USE mo_grf_intp_data_strc,  ONLY: p_grf_state_global, p_grf_state_subdiv, p_grf_state


USE mo_read_namelists,     ONLY: read_atmo_namelists
USE mo_atm_nml_crosscheck,       ONLY: atm_crosscheck

USE mo_time_config,        ONLY: time_config      ! variable
USE mo_dynamics_config,    ONLY: configure_dynamics  ! subroutine
USE mo_interpol_config

USE mo_atmo_nonhydrostatic, ONLY: atmo_nonhydrostatic 

! USE statements referring directly to prep_icon
USE mo_time_config, ONLY: time_config 
USE mo_prepicon_utils,      ONLY: init_prepicon, prepicon, write_prepicon_output, &
  & compute_coord_fields, init_topo_output_files, close_prepicon_output_files,    &
  & convert_variables, init_atmo_output_files, deallocate_prepicon
  
USE mo_prepicon_nml,        ONLY: i_oper_mode, prepicon_nml_setup, l_zp_out
USE mo_nh_vert_interp,      ONLY: vertical_interpolation, interpolate_to_p_and_z_levels

USE mo_extpar_config,      ONLY: itopo

!-------------------------------------------------------------------------
IMPLICIT NONE

  ! local variables
  CHARACTER(*), PARAMETER :: routine = "prep_icon"

    
  INTEGER :: ist
  INTEGER :: error_status

  !declaration of OpenMP Runtime Library Routines:
!$  INTEGER omp_get_max_threads
!$  INTEGER omp_get_num_threads
!$  INTEGER omp_get_num_procs
!$  INTEGER omp_get_thread_num
!$  LOGICAL omp_get_dynamic

!$  INTEGER :: max_threads_omp, num_procs_omp
!$  LOGICAL :: ldynamic_omp

  !--------------------------------------------------------------------
  !BOC

  !print out some information about OpenMP parallelization
!$  max_threads_omp  = omp_get_max_threads()
!$  num_procs_omp    = omp_get_num_procs()
!$  ldynamic_omp     = omp_get_dynamic()
!$  WRITE(message_text,'(A,I3,A,I3)')                &
!$    & "OpenMP:  MAX_THREADS = ", max_threads_omp,  &
!$    & ",  NUM_PROCS = ", num_procs_omp
!$  CALL message('control_model',message_text)
!$  WRITE(message_text,'(A,L3)')  &
!$    & "OpenMP:  DYNAMIC = ", ldynamic_omp
!$  CALL message('control_model',message_text)


#ifdef __INTEL_COMPILER
  ! Important on Intel: disable underflow exceptions:
  CALL disable_ufl_exception
#endif

    !-------------------------------------------------------------------
    ! Initialize MPI, this should aleays be the first call
    CALL start_mpi('PREP_ICON')

    !---------------------------------------------------------------------
    ! 1.1 Read namelists (newly) specified by the user; fill the 
    !     corresponding sections of the configuration states.
    !---------------------------------------------------------------------

    CALL read_atmo_namelists("NAMELIST_PREPICON","icon_master.namelist")

    CALL prepicon_nml_setup("NAMELIST_PREPICON")

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

    ! This is hardcoded for prep_icon because running this program is nonsense otherwise
    itopo = 1


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
    
    CALL import_patches( p_patch_global,nlev,nlevp1,num_lev,num_levp1,nshift)


    !--------------------------------------------------------------------------------
    ! 5. Construct interpolation state, compute interpolation coefficients.
    !--------------------------------------------------------------------------------
    
    CALL configure_interpolation( global_cell_type, n_dom, p_patch_global(1:)%level )

    ! Allocate array for interpolation state
    
    ALLOCATE( p_int_state_global(n_dom_start:n_dom), &
            & p_grf_state_global(n_dom_start:n_dom),STAT=error_status)
    IF (error_status /= SUCCESS) THEN
      CALL finish(TRIM(routine),'allocation for ptr_int_state failed')
    ENDIF
    

    ! Interpolation state is constructed for
    ! the full domain on every PE and divided later
    CALL construct_2d_interpol_state(p_patch_global, p_int_state_global)

    !-----------------------------------------------------------------------------
    ! 6. Construct grid refinment state, compute coefficients
    !-----------------------------------------------------------------------------
    ! For the NH model, the initialization routines called from
    ! construct_2d_gridref_state require the metric terms to be present

    IF (n_dom_start==0 .OR. n_dom > 1) THEN
      CALL construct_2d_gridref_state (p_patch_global, p_grf_state_global)
    ENDIF

    !-------------------------------------------------------------------
    ! 7. Domain decomposition: 
    !    Divide patches and interpolation states for parallel runs.
    !    This is only done if the model runs really in parallel.
    !-------------------------------------------------------------------   
    IF (my_process_is_mpi_seq()) THEN
      
      ! This is a verification run or a run on a single processor
      ! or the divided states have been read, just set pointers
      
      p_patch => p_patch_global
      p_int_state => p_int_state_global
      p_grf_state => p_grf_state_global
      
!       IF (my_process_is_mpi_seq()) THEN
!         p_patch(:)%comm = p_comm_work
!       ELSE
        CALL set_patch_communicators(p_patch)
!       ENDIF
      
    ELSE
      
      CALL decompose_atmo_domain()
      
    ENDIF

    ! Note: fro this point the p_patch is used
    ! In case of a test run: Copy processor splitting to test PE
    IF(p_test_run) CALL copy_processor_splitting(p_patch)
    


   !---------------------------------------------------------------------
   ! 8. Import vertical grid/ define vertical coordinate
   !---------------------------------------------------------------------
    SELECT CASE (iequations)
      
    CASE (inh_atmosphere)
      IF (ivctype == 1) THEN
        CALL init_hybrid_coord(p_patch(1)%nlev)
      ELSE IF (ivctype == 2) THEN
        CALL init_sleve_coord(p_patch(1)%nlev)
      ENDIF
    CASE DEFAULT
      CALL finish(TRIM(routine),'prep_icon only works for NH model')
    END SELECT


    !---------------------------------------------------------------------
    ! 9. Horizontal and vertical grid(s) are now defined.
    !    Assign values to derived variables in the configuration states
    !---------------------------------------------------------------------

    CALL configure_dynamics ( n_dom )

    ! Allocate prepicon data type
    ALLOCATE (prepicon(n_dom), stat=ist)
    IF (ist /= success) THEN
      CALL finish(routine,'allocation for prepicon failed')
    ENDIF
    
    ! allocate memory for topography and coordinate fields,
    ! read topo data from netCDF file, 
    ! optionally smooth topography data,
    ! and, in case of nesting, perform topography blending and feedback
    CALL init_prepicon (p_int_state(1:), p_grf_state(1:), prepicon)

    ! Compute the 3D coordinate fields
    CALL compute_coord_fields(p_int_state(1:), prepicon)

    IF (i_oper_mode == 2) THEN

      ! Convert (hydrostatic) IFS2ICON output variables to the NH set of prognostic variables
      CALL convert_variables(p_int_state(1:), prepicon)

    ELSE IF (i_oper_mode == 3) THEN
      ! Perform vertical interpolation from intermediate IFS2ICON grid to ICON grid
      ! and convert variables to the NH set of prognostic variables
      CALL vertical_interpolation(p_patch(1:), p_int_state(1:), p_grf_state(1:), prepicon)
    ENDIF

    IF (i_oper_mode >= 2 .AND. l_zp_out) THEN
      ! Interpolate prognostic variables to pressure and height levels for diagnostic output
      CALL interpolate_to_p_and_z_levels(p_patch(1:), p_int_state(1:), prepicon)
    ENDIF

    !------------------------------------------------------------------
    ! Write output
    !------------------------------------------------------------------
 
    IF (i_oper_mode == 1) THEN

      CALL init_topo_output_files
      CALL write_prepicon_output( time_config%ini_datetime )
      CALL close_prepicon_output_files

    ELSE

      CALL init_atmo_output_files
      CALL write_prepicon_output( time_config%ini_datetime )
      CALL close_prepicon_output_files

    ENDIF

    ! Deallocate prepicon data type
    CALL deallocate_prepicon(prepicon)
    DEALLOCATE (prepicon, stat=ist)
    IF (ist /= success) THEN
      CALL finish(routine,'deallocation for prepicon failed')
    ENDIF



    ! Deconstruct grid refinement state

    IF (n_dom > 1) THEN
      CALL destruct_2d_gridref_state( p_patch, p_grf_state )
    ENDIF

    IF (my_process_is_mpi_seq()) THEN
      DEALLOCATE (p_grf_state_global, STAT=error_status)
    ELSE
      DEALLOCATE (p_grf_state_subdiv, STAT=error_status)
    ENDIF
    IF (error_status /= SUCCESS) THEN
      CALL finish(TRIM(routine),'deallocation for ptr_grf_state failed')
    ENDIF

    ! Deallocate interpolation fields

    CALL destruct_2d_interpol_state( p_int_state )
    IF  (my_process_is_mpi_seq()) THEN
      DEALLOCATE (p_int_state_global, STAT=error_status)
    ELSE
      DEALLOCATE (p_int_state_subdiv, STAT=error_status)
    ENDIF
    IF (error_status /= SUCCESS) THEN
      CALL finish(TRIM(routine),'deallocation for ptr_int_state failed')
    ENDIF

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


  CALL p_stop

END PROGRAM prep_icon
