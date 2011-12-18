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
MODULE mo_cpl_dummy_model

USE mo_kind,                ONLY: wp
USE mo_exception,           ONLY: message, message_text, finish
USE mo_parallel_config,     ONLY: p_test_run, l_test_openmp, num_io_procs
USE mo_mpi,                 ONLY: p_stop, &
  & my_process_is_io,  my_process_is_mpi_seq, my_process_is_mpi_test, &
  & my_process_is_mpi_parallel,                                       &
  & set_mpi_work_communicators, set_comm_input_bcast, null_comm_type, &
  & global_mpi_barrier, p_pe_work
USE mo_timer,               ONLY: init_timer
USE mo_master_control,      ONLY: is_restart_run, get_my_process_name, &
                                  get_my_model_no


! Control parameters: run control, dynamics, i/o
!
USE mo_nonhydrostatic_config,ONLY: ivctype, kstart_moist, kstart_qv,    &
  &                                kend_qvsubstep, l_open_ubc,          &
  &                                configure_nonhydrostatic
USE mo_dynamics_config,   ONLY: iequations
USE mo_run_config,        ONLY: configure_run, &
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
  & ishallow_water


! For the coupling
USE mo_impl_constants, ONLY: CELLS
USE mo_coupling_config, ONLY : is_coupled_run
USE mo_icon_cpl_def_grid, ONLY : ICON_cpl_def_grid, ICON_cpl_def_location
USE mo_icon_cpl_def_field, ONLY : ICON_cpl_def_field
USE mo_icon_cpl_search, ONLY : ICON_cpl_search
USE mo_icon_cpl_exchg, ONLY : ICON_cpl_put, ICON_cpl_get
USE mo_icon_cpl_finalize,   ONLY: icon_cpl_finalize
USE mo_model_domain_import, ONLY : get_patch_global_indexes

! Memory
!
USE mo_subdivision,         ONLY: decompose_domain,         &
& complete_parallel_setup, &
& finalize_decomposition,  &
& copy_processor_splitting

USE mo_model_domain,        ONLY: t_patch, p_patch, p_patch_local_parent

! Horizontal grid
USE mo_grid_config,         ONLY: n_dom, n_dom_start, global_cell_type, &
                                  dynamics_parent_grid_id
USE mo_model_domimp_patches,ONLY: import_basic_patches, complete_patches, destruct_patches

! Horizontal interpolation
!
USE mo_intp_state,          ONLY: construct_2d_interpol_state,  destruct_2d_interpol_state, &
  &                               transfer_interpol_state
USE mo_grf_interpolation,   ONLY: construct_2d_gridref_state,   destruct_2d_gridref_state
USE mo_grf_intp_state,      ONLY: transfer_grf_state

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
USE mo_ext_data,            ONLY: ext_data, init_ext_data, destruct_ext_data

!  USE mo_nwp_phy_init,          ONLY: init_nwp_phy
!!$  USE mo_gscp_cosmo,          ONLY: hydci_pp_init


!-------------------------------------------------------------------------
! to break circular dependency

USE mo_intp_data_strc,      ONLY: p_int_state, p_int_state_local_parent
USE mo_grf_intp_data_strc,  ONLY: p_grf_state, p_grf_state_local_parent

!-------------------------------------------------------------------------
USE mo_io_restart,           ONLY: read_restart_info_file
USE mo_io_restart_namelist,  ONLY: read_restart_namelists
USE mo_io_restart_attributes,ONLY: read_restart_attributes

USE mo_read_namelists,     ONLY: read_cpl_dummy_namelists
USE mo_nml_crosscheck,       ONLY: atm_crosscheck

USE mo_dynamics_config,    ONLY: configure_dynamics  ! subroutine
!USE mo_interpol_config,    ONLY: configure_interpolation 
USE mo_interpol_config
USE mo_advection_config,   ONLY: configure_advection
USE mo_diffusion_config,   ONLY: configure_diffusion

USE mo_atmo_hydrostatic,    ONLY: atmo_hydrostatic 
USE mo_atmo_nonhydrostatic, ONLY: atmo_nonhydrostatic 

!-------------------------------------------------------------------------
IMPLICIT NONE
PRIVATE

PUBLIC :: cpl_dummy_model

CONTAINS
!>
!!
  SUBROUTINE cpl_dummy_model(cpl_dummy_namelist_filename,shr_namelist_filename)

    CHARACTER(LEN=*), INTENT(in) :: cpl_dummy_namelist_filename
    CHARACTER(LEN=*), INTENT(in) :: shr_namelist_filename

    CHARACTER(*), PARAMETER :: method_name = "mo_cpl_dummy_model:cpl_dummy_model"
    INTEGER :: jg, jgp

    TYPE(t_patch), ALLOCATABLE :: p_patch_global(:)

    ! For the coupling

    INTEGER, PARAMETER :: no_of_fields = 8

    CHARACTER(LEN=MAX_CHAR_LENGTH) ::  field_name(no_of_fields)
    INTEGER :: field_id(no_of_fields)
    INTEGER :: grid_id
    INTEGER :: grid_shape(2) 
    INTEGER :: field_shape(3) 
    INTEGER :: i, error_status
!rr INTEGER :: no_of_entities
    INTEGER :: patch_no

    !---------------------------------------------------------------------
    ! 1.1 Read namelists (newly) specified by the user; fill the 
    !     corresponding sections of the configuration states.
    !---------------------------------------------------------------------

    CALL global_mpi_barrier()
    write(0,*) TRIM(get_my_process_name()), ': Start of ', method_name
!     CALL p_stop
!     STOP
    
    CALL read_cpl_dummy_namelists(cpl_dummy_namelist_filename,shr_namelist_filename)
    write(0,*) TRIM(get_my_process_name()), ': namelist is read '
    CALL global_mpi_barrier()

    !---------------------------------------------------------------------
    ! 1.2 Cross-check namelist setups
    !---------------------------------------------------------------------

    CALL atm_crosscheck
    write(0,*) TRIM(get_my_process_name()), ': atm_crosscheck is done '
    CALL global_mpi_barrier()

    !---------------------------------------------------------------------
    ! 2. Call configure_run to finish filling the run_config state.
    !    This needs to be done very early (but anyway after atm_crosscheck)
    !    because some component of the state, e.g., num_lev, may be 
    !    modified in this subroutine which affect the following CALLs.
    !---------------------------------------------------------------------
    CALL configure_run
    write(0,*) TRIM(get_my_process_name()), ': configure_run is done '
    CALL global_mpi_barrier()

    !-------------------------------------------------------------------
    ! 3.1 Initialize the mpi work groups
    !-------------------------------------------------------------------
    CALL set_mpi_work_communicators(p_test_run, l_test_openmp, num_io_procs)
    write(0,*) TRIM(get_my_process_name()), ': set_mpi_work_communicators is done '
    CALL global_mpi_barrier()

    !-------------------------------------------------------------------
    ! 3.2 Initialize various timers
    !-------------------------------------------------------------------
    IF (ltimer) CALL init_timer
    write(0,*) TRIM(get_my_process_name()), ': init_timer is done '
    CALL global_mpi_barrier()

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
    IF ( ALLOCATED(p_patch)) THEN
      CALL finish(TRIM(method_name), 'patch already allocated')
    END IF
     
    ! Allocate patch array to start patch construction
    ALLOCATE(p_patch(n_dom_start:n_dom), stat=error_status)
    IF (error_status/=success) THEN
      CALL finish(TRIM(method_name), 'allocation of patch failed')
    ENDIF
    
    IF(my_process_is_mpi_parallel()) THEN
      ALLOCATE(p_patch_global(n_dom_start:n_dom))
      CALL import_basic_patches(p_patch_global,nlev,nlevp1,num_lev,num_levp1,nshift)
      CALL decompose_domain(p_patch_global)
      DEALLOCATE(p_patch_global)
      CALL complete_parallel_setup()
    ELSE
      CALL import_basic_patches(p_patch,nlev,nlevp1,num_lev,num_levp1,nshift)
    ENDIF

    ! Complete information which is not yet read or calculated
    CALL complete_patches( p_patch )

    ! In case of a test run: Copy processor splitting to test PE
    IF(p_test_run) CALL copy_processor_splitting(p_patch)

    !--------------------------------------------------------------------------------
    ! 5. Construct interpolation state, compute interpolation coefficients.
    !--------------------------------------------------------------------------------
    
    CALL configure_interpolation( global_cell_type, n_dom, p_patch(1:)%level )
    write(0,*) TRIM(get_my_process_name()), ': configure_interpolation is done '
    CALL global_mpi_barrier()

    ! Allocate array for interpolation state
    
    ALLOCATE( p_int_state(n_dom_start:n_dom), &
            & p_grf_state(n_dom_start:n_dom),STAT=error_status)
    IF (error_status /= SUCCESS) THEN
      CALL finish(TRIM(method_name),'allocation for ptr_int_state failed')
    ENDIF

    IF(my_process_is_mpi_parallel()) THEN
      ALLOCATE( p_int_state_local_parent(n_dom_start+1:n_dom), &
              & p_grf_state_local_parent(n_dom_start+1:n_dom), &
              & STAT=error_status)
      IF (error_status /= SUCCESS) &
        CALL finish(TRIM(method_name),'allocation for local parents failed')
    ENDIF
    
    ! Construct interpolation state
    ! Please note that for pararllel runs the divided state is constructed here
    CALL construct_2d_interpol_state(p_patch, p_int_state)
    write(0,*) TRIM(get_my_process_name()), ': construct_2d_interpol_state is done '
    CALL global_mpi_barrier()

    IF(my_process_is_mpi_parallel()) THEN
      ! Transfer interpolation state to local parent
      DO jg = n_dom_start+1, n_dom
        jgp = p_patch(jg)%parent_id
        CALL transfer_interpol_state(p_patch(jgp),p_patch_local_parent(jg), &
                                  &  p_int_state(jgp), p_int_state_local_parent(jg))
      ENDDO
    ENDIF

    !-----------------------------------------------------------------------------
    ! 6. Construct grid refinment state, compute coefficients
    !-----------------------------------------------------------------------------
    ! For the NH model, the initialization routines called from
    ! construct_2d_gridref_state require the metric terms to be present

    IF (n_dom_start==0 .OR. n_dom > 1) THEN
      ! Construct gridref state
      ! For non-parallel runs (either 1 PE only or on the test PE) this is done as usual,
      ! for parallel runs, the main part of the gridref state is constructed on the
      ! local parent with the following call
      CALL construct_2d_gridref_state (p_patch, p_grf_state)
      IF(my_process_is_mpi_parallel()) THEN
        ! Transfer gridref state from local parent to p_grf_state
        DO jg = n_dom_start+1, n_dom
          jgp = p_patch(jg)%parent_id
          CALL transfer_grf_state(p_patch(jgp), p_patch_local_parent(jg), &
                               &  p_grf_state(jgp), p_grf_state_local_parent(jg), &
                               &  p_patch(jg)%parent_child_index)
        ENDDO
      ENDIF
    ENDIF

    !-------------------------------------------------------------------
    ! 7. Finalize domain decomposition
    !-------------------------------------------------------------------   
    IF (my_process_is_mpi_parallel()) THEN

      CALL finalize_decomposition()
      write(0,*) TRIM(get_my_process_name()), ': decompose_domain is done '
      CALL global_mpi_barrier()
      
    ENDIF

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

    CALL configure_dynamics ( n_dom )
    CALL configure_diffusion( n_dom, dynamics_parent_grid_id, &
                            & nlev, vct_a, vct_b, apzero      )

    IF (iequations == inh_atmosphere) THEN
      DO jg =1,n_dom
        CALL configure_nonhydrostatic(jg, p_patch(jg)%nlev, p_patch(jg)%nshift_total)
      ENDDO
    ENDIF

    DO jg =1,n_dom
     CALL configure_advection( jg, p_patch(jg)%nlev, p_patch(1)%nlev,        &
       &                      iequations, iforcing, iqv, kstart_moist(jg),   &
       &                      kstart_qv(jg), kend_qvsubstep(jg), lvert_nest, &
       &                      l_open_ubc, ntracer ) 
    ENDDO

    !------------------------------------------------------------------
    ! 10. Create and optionally read external data fields
    !------------------------------------------------------------------
    ALLOCATE (ext_data(n_dom), STAT=error_status)
    IF (error_status /= SUCCESS) THEN
      CALL finish(TRIM(method_name),'allocation for ext_data failed')
    ENDIF
    
    ! allocate memory for atmospheric/oceanic external data and
    ! optionally read those data from netCDF file.
    CALL init_ext_data (p_patch(1:), p_int_state(1:), ext_data)

    !---------------------------------------------------------------------
    ! 11. Do the setup for the coupled run
    !
    ! For the time being this could all go into a subroutine which is
    ! common to atmo and ocean. Does this make sense if the setup deviates
    ! too much in future.
    !---------------------------------------------------------------------
    IF ( is_coupled_run() ) THEN
 
      patch_no = 1

      grid_shape(1)=1
      grid_shape(2)=p_patch(patch_no)%n_patch_cells

      ! CALL get_patch_global_indexes ( patch_no, CELLS, no_of_entities, grid_glob_index )
      ! should grid_glob_index become a pointer in ICON_cpl_def_grid as well?

       CALL ICON_cpl_def_grid ( &
         & grid_shape, p_patch(patch_no)%cells%glb_index, & ! input
         & grid_id, error_status )                          ! output

      ! Marker for internal and halo points, a list which contains the
      ! rank where the native vertices are located.
      CALL ICON_cpl_def_location ( &
        & grid_id, grid_shape, p_patch(patch_no)%cells%owner_local, & ! input
        & p_pe_work,  & ! this owner id
        & error_status )                                            ! output
      write(0,*) "p_pe_work:",  p_pe_work
      DO i=grid_shape(1),grid_shape(2)
        write(0,*) i, "global idx:", p_patch(patch_no)%cells%glb_index(i),&
         & " owner:",  p_patch(patch_no)%cells%owner_local(i)
      ENDDO
      ! Define exchange fields
      !
      ! We could also get the number of defined fields from the namelist module
      ! and assign the names according the name. In practice this is a bit dangerous
      ! as fields may require special treatment which is only known to the programmer
      ! but not to the one who possibly modifies the namelist.

      field_name(1) = "TEST1"
      field_name(2) = "TEST2"
      field_name(3) = "TEST3"
      field_name(4) = "TEST4"
      field_name(5) = "TEST5"
      field_name(6) = "TEST6"
      field_name(7) = "TEST7"
      field_name(8) = "TEST8"

      ! horizontal dimension of the exchange field

      field_shape(1) = grid_shape(1)
      field_shape(2) = grid_shape(2)

      ! number of bundles or vertical levels

      field_shape(3) = 1

      DO i = 1, no_of_fields
          CALL ICON_cpl_def_field ( field_name(i), grid_id, field_id(i), &
                                  & field_shape, error_status )
      ENDDO

      CALL ICON_cpl_search

    ENDIF

    !---------------------------------------------------------------------
    ! 12. Call dummy procedure
    ! This procedure simulates the timestepping for the sent - receive
    ! coupling process
    !---------------------------------------------------------------------

    CALL cpl_dummy_timestepping(no_of_fields, field_id, field_shape)

    !---------------------------------------------------------------------
    ! 13. Dummy Carry out the shared clean-up processes
    !---------------------------------------------------------------------
    ! Destruct external data state

    CALL destruct_ext_data

    ! deallocate ext_data array
    DEALLOCATE(ext_data, stat=error_status)
    IF (error_status/=success) THEN
      CALL finish(TRIM(method_name), 'deallocation of ext_data')
    ENDIF

    ! Deconstruct grid refinement state

    IF (n_dom > 1) THEN
      CALL destruct_2d_gridref_state( p_patch, p_grf_state )
    ENDIF

    DEALLOCATE (p_grf_state, STAT=error_status)
    IF (error_status /= SUCCESS) THEN
      CALL finish(TRIM(method_name),'deallocation for ptr_grf_state failed')
    ENDIF

    ! Deallocate interpolation fields

    CALL destruct_2d_interpol_state( p_int_state )
    DEALLOCATE (p_int_state, STAT=error_status)
    IF (error_status /= SUCCESS) THEN
      CALL finish(TRIM(method_name),'deallocation for ptr_int_state failed')
    ENDIF

    ! Deallocate grid patches

    CALL destruct_patches( p_patch )

    DEALLOCATE( p_patch, STAT=error_status )
    IF (error_status/=SUCCESS) THEN
      CALL finish(TRIM(method_name),'deallocate for patch array failed')
    ENDIF
  
    IF ( is_coupled_run() ) CALL ICON_cpl_finalize

    CALL message(TRIM(method_name),'clean-up finished')

  END SUBROUTINE cpl_dummy_model
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !> This method simulates the timestepping for the sent - receive
  ! coupling process
  SUBROUTINE cpl_dummy_timestepping(no_of_fields, field_id, field_shape)

    INTEGER, INTENT(IN) :: no_of_fields
    INTEGER, INTENT(IN) :: field_id(no_of_fields)
    INTEGER, INTENT(IN) :: field_shape(3)

    CHARACTER(*), PARAMETER :: method_name = "mo_cpl_dummy_model:cpl_dummy_model"

    INTEGER, PARAMETER :: nfld_fix = 2

    INTEGER       :: nfld, i, nb
    INTEGER       :: patch_no
    INTEGER       :: id
    INTEGER       :: ierror
    INTEGER       :: info
    INTEGER       :: model_id

    REAL(kind=wp) :: recv_field (field_shape(1):field_shape(2),field_shape(3))
    REAL(kind=wp) :: send_field (field_shape(1):field_shape(2),field_shape(3))

    CALL global_mpi_barrier()
    write(0,*) TRIM(get_my_process_name()), ': start cpl_dummy_timestepping... '
    CALL global_mpi_barrier()

    model_id = get_my_model_no()
    patch_no = 1

    ! what is the proper conversion from int to wp?

    SELECT CASE (model_id)

    CASE ( 1 )  

       DO nfld = 1, 4

          IF ( nfld == nfld_fix ) THEN
             DO nb = 1, field_shape(3)
                DO i = field_shape(1), field_shape(2)
                   send_field(i,nb) = REAL(p_patch(patch_no)%cells%glb_index(i),wp)
                ENDDO
             ENDDO
          ELSE
             DO nb = 1, field_shape(3)
                DO i = field_shape(1), field_shape(2)
                   send_field(i,nb) = real(nfld,wp)
                ENDDO
             ENDDO
          ENDIF

          id = field_id(nfld)
          CALL ICON_cpl_put ( id, field_shape, send_field, ierror )
       ENDDO

       DO nfld = 5, no_of_fields
          id = field_id(nfld)
          CALL ICON_cpl_get ( id, field_shape, recv_field, info, ierror )
       ENDDO

    CASE ( 2 )

       DO nfld = 5, no_of_fields

          IF ( nfld == 4 + nfld_fix ) THEN
             DO nb = 1, field_shape(3)
                DO i = field_shape(1), field_shape(2)
                   send_field(i,nb) = REAL(p_patch(patch_no)%cells%glb_index(i),wp)
                ENDDO
             ENDDO
          ELSE
             DO nb = 1, field_shape(3)
                DO i = field_shape(1), field_shape(2)
                   send_field(i,nb) = real(nfld,wp)
                ENDDO
             ENDDO
          ENDIF

          id = field_id(nfld)
          CALL ICON_cpl_put ( id, field_shape, send_field, ierror )

       ENDDO

       DO nfld = 1, 4

          id = field_id(nfld)
          CALL ICON_cpl_get ( id, field_shape, recv_field, info, ierror )

          IF ( info > 0 ) THEN

             ! Control received results
             ! Exchange routine only work on the inner points, not on the halo!

             IF ( nfld == nfld_fix ) THEN
                DO nb = 1, field_shape(3)
                   DO i = field_shape(1), field_shape(2)

                      IF ( p_patch(patch_no)%cells%owner_local(i) == p_pe_work ) THEN
                        IF (recv_field(i,nb) /= REAL(p_patch(patch_no)%cells%glb_index(i),wp)) THEN
                            WRITE(message_text,'(i6,a11,i6,a9,f13.4)') i, ': Expected ',    &
                                 &                    p_patch(patch_no)%cells%glb_index(i), &
                                 &                    ' but got ',                          &
                                 &                      recv_field(i,nb)
                            CALL message('mo_cpl_dummy_model', TRIM(message_text))
                        ENDIF
                      ENDIF

                   ENDDO
                ENDDO
             ELSE
                DO nb = 1, field_shape(3)
                   DO i = field_shape(1), field_shape(2)

                      IF ( p_patch(patch_no)%cells%owner_local(i) == p_pe_work ) THEN
                         IF ( recv_field(i,nb) /= real(nfld,wp) ) THEN
                            WRITE(message_text,'(i6,a11,i6,a9,f13.4)') i, ': Expected ',   &
                                 &                            nfld,                        &
                                 &                           ' but got ',                  &
                                 &                            recv_field(i,nb)
                            CALL message('mo_cpl_dummy_model', TRIM(message_text))
                         ENDIF
                      ENDIF

                   ENDDO
                ENDDO
             ENDIF

          ENDIF ! info > 0

       ENDDO

    CASE DEFAULT

       CALL message(TRIM(method_name),'more than 2 dummy comps are not supported')    
       CALL finish (TRIM(method_name),'deallocate for patch array failed')

    END SELECT

    CALL global_mpi_barrier()
    write(0,*) TRIM(get_my_process_name()), ': cpl_dummy_timestepping ends '

  END SUBROUTINE cpl_dummy_timestepping
  !-------------------------------------------------------------------------


END MODULE mo_cpl_dummy_model

