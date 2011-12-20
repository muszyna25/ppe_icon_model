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
MODULE mo_test_coupler

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: message, message_text, finish
  USE mo_mpi,                 ONLY: global_mpi_barrier, p_pe_work
  USE mo_timer,               ONLY: init_timer
  USE mo_master_control,      ONLY: is_restart_run, get_my_process_name, &
                                    get_my_model_no
  ! For the coupling
  USE mo_icon_cpl_init,       ONLY: icon_cpl_init
  USE mo_icon_cpl_init_comp,  ONLY: icon_cpl_init_comp
  USE mo_impl_constants,      ONLY : MAX_CHAR_LENGTH
! USE mo_impl_constants,      ONLY: SUCCESS
  USE mo_coupling_config,     ONLY : is_coupled_run, config_debug_coupler_level
  USE mo_icon_cpl_def_grid,   ONLY : ICON_cpl_def_grid, ICON_cpl_def_location
  USE mo_icon_cpl_def_field,  ONLY : ICON_cpl_def_field
  USE mo_icon_cpl_search,     ONLY : ICON_cpl_search
  USE mo_icon_cpl_exchg,      ONLY : ICON_cpl_put, ICON_cpl_get
  USE mo_icon_cpl_finalize,   ONLY : icon_cpl_finalize
  USE mo_model_domain_import, ONLY : get_patch_global_indexes
!   USE mo_read_namelists,      ONLY: read_cpl_dummy_namelists

  USE mo_model_domain,        ONLY:  p_patch
  
  USE mo_atmo_model,          ONLY: construct_atmo_model, destruct_atmo_model

!-------------------------------------------------------------------------
IMPLICIT NONE
PRIVATE

PUBLIC :: test_coupler

CONTAINS
!>
!!
  SUBROUTINE test_coupler(cpl_dummy_namelist_filename,shr_namelist_filename)

    CHARACTER(LEN=*), INTENT(in) :: cpl_dummy_namelist_filename
    CHARACTER(LEN=*), INTENT(in) :: shr_namelist_filename

    CHARACTER(*), PARAMETER :: method_name = "mo_test_coupler:test_coupler"

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
    
    !---------------------------------------------------------------------

    CALL construct_atmo_model(cpl_dummy_namelist_filename,shr_namelist_filename)
    
!     CALL read_cpl_dummy_namelists(cpl_dummy_namelist_filename,shr_namelist_filename)
!     write(0,*) TRIM(get_my_process_name()), ': namelist is read '
    CALL global_mpi_barrier()
    
    !---------------------------------------------------------------------
    ! 11. Do the setup for the coupled run
    !
    ! For the time being this could all go into a subroutine which is
    ! common to atmo and ocean. Does this make sense if the setup deviates
    ! too much in future.
    !---------------------------------------------------------------------
    IF ( is_coupled_run() ) THEN
 
      !------------------------------------------------------------
      CALL icon_cpl_init(debug_level=config_debug_coupler_level)
      ! Inform the coupler about what we are
      CALL icon_cpl_init_comp ( get_my_process_name(), get_my_model_no(), error_status )
      ! split the global_mpi_communicator into the components
      !------------------------------------------------------------

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
    CALL destruct_atmo_model()
  
    IF ( is_coupled_run() ) CALL ICON_cpl_finalize

    CALL message(TRIM(method_name),'clean-up finished')

  END SUBROUTINE test_coupler
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


END MODULE mo_test_coupler

