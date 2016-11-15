  !>
  !! This module contains the I/O routines for lateral boundary nudging
  !!
  !! @author M. Pondkule (DWD)
  !!
  !!
  !! @par Revision History
  !! Initial release by M. Pondkule, DWD (2013-10-31)
  !!
  !!
  !! @par Copyright
  !! 2002-2013 by DWD and MPI-M
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

MODULE mo_async_latbc

  USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_ptr, c_intptr_t, c_f_pointer

#ifndef NOMPI
    USE mpi
#endif

    ! basic modules
    USE mo_kind,                      ONLY: i8, sp
    USE mo_exception,                 ONLY: finish, message
    USE mo_mpi,                       ONLY: stop_mpi, my_process_is_io,  my_process_is_pref, &
         &                                  my_process_is_mpi_test, p_int, p_real_sp
    USE mo_parallel_config,           ONLY: nproma
    USE mo_model_domain,              ONLY: p_patch
    USE mo_model_domain,              ONLY: t_patch
    USE mo_nonhydro_types,            ONLY: t_nh_state
    USE mo_intp_data_strc,            ONLY: t_int_state
    ! MPI Communicators
    USE mo_mpi,                       ONLY: p_comm_work, p_comm_work_pref, p_comm_work_2_pref
    ! MPI Communication routines
    USE mo_mpi,                       ONLY: p_bcast
    ! MPI Process type intrinsics
    USE mo_mpi,                       ONLY: my_process_is_work
    ! MPI Process group sizes
    USE mo_mpi,                       ONLY: num_work_procs, p_n_work
    ! Processor numbers
    USE mo_mpi,                       ONLY: p_pe_work, p_work_pe0, p_comm_work_pref_compute_pe0
    USE mo_time_config,               ONLY: time_config
    USE mo_datetime,                  ONLY: t_datetime
    USE mo_async_latbc_types,         ONLY: t_patch_data, t_reorder_data, latbc_buffer
    USE mo_grid_config,               ONLY: nroot
    USE mo_async_latbc_utils,         ONLY: pref_latbc_data, prepare_pref_latbc_data, &
         &                                  compute_wait_for_async_pref, compute_shutdown_async_pref, &
         &                                  async_pref_send_handshake,  async_pref_wait_for_start
    USE mo_impl_constants,            ONLY: SUCCESS, MAX_CHAR_LENGTH, MODE_DWDANA, MODE_ICONVREMAP, &
                                            MODE_IAU_OLD, MODE_IAU
    USE mo_communication,             ONLY: idx_no, blk_no
    USE mo_nonhydro_state,            ONLY: p_nh_state
    USE mo_intp_data_strc,            ONLY: p_int_state
    USE mo_ext_data_state,            ONLY: ext_data
    USE mo_linked_list,               ONLY: t_var_list, t_list_element
    USE mo_var_metadata_types,        ONLY: t_var_metadata, VARNAME_LEN
    USE mo_var_list,                  ONLY: nvar_lists, var_lists, new_var_list, &
         &                                  collect_group
    USE mo_limarea_config,            ONLY: latbc_config, generate_filename, t_glb_indices
    USE mo_dictionary,                ONLY: t_dictionary, dict_get, dict_init, dict_loadfile, &
         &                                  dict_finalize
    USE mo_util_string,               ONLY: add_to_list, tolower
    USE mo_initicon_config,           ONLY: init_mode
    USE mo_time_config,               ONLY: time_config
    USE mo_cdi,                       ONLY: vlistInqVarZaxis , streamOpenRead, streamInqVlist, &
         &                                  vlistNvars, zaxisInqSize, vlistInqVarName,         &
         &                                  vlistInqVarGrid, streamClose, streamInqFiletype,   &
         &                                  FILETYPE_NC2, FILETYPE_NC4, FILETYPE_GRB2
    USE mo_cdi_constants,             ONLY: GRID_UNSTRUCTURED_CELL, GRID_UNSTRUCTURED_EDGE
    USE mo_io_units,                  ONLY: filename_max, nerr
    USE mo_io_util,                   ONLY: read_netcdf_int_1d
    USE mo_util_file,                 ONLY: util_filesize
    USE mo_util_cdi,                  ONLY: test_cdi_varID, cdiGetStringError

    IMPLICIT NONE

    PRIVATE

    ! subroutines
    PUBLIC :: latbc_buffer
    PUBLIC :: prefetch_input
    PUBLIC :: prefetch_main_proc
    PUBLIC :: init_prefetch
    PUBLIC :: close_prefetch

    !------------------------------------------------------------------------------------------------
    ! CONSTANTS
    !------------------------------------------------------------------------------------------------
    ! grp name in lower case letter
    CHARACTER(LEN=VARNAME_LEN), ALLOCATABLE :: StrLowCasegrp(:)

    ! maximum text lengths in this module
    INTEGER, PARAMETER :: MAX_ERROR_LENGTH  = 256

    ! NetCDF file IDs / CDI stream IDs for first guess and analysis file
    INTEGER, ALLOCATABLE :: fileID_fg(:)

    !> constant for better readability
    INTEGER, PARAMETER :: WAIT_UNTIL_FINISHED = -1
    ! common constant strings
    CHARACTER(LEN=*), PARAMETER :: ALLOCATE_FAILED   = 'ALLOCATE failed!'
    CHARACTER(LEN=*), PARAMETER :: UNKNOWN_GRID_TYPE = 'Unknown grid type!'
    CHARACTER(LEN=*), PARAMETER :: modname = 'mo_async_latbc'

    TYPE(t_patch_data), PUBLIC, SAVE, TARGET :: patch_data

    !------------------------------------------------------------------------------------------------
    ! Broadcast root for intercommunicator broadcasts form compute PEs to prefetching
    ! PE using p_comm_work_2_pref
    INTEGER :: bcast_root

  CONTAINS

    !------------------------------------------------------------------------------------------------
    !> Close all name_list files and deallocate variables
    !
    SUBROUTINE close_prefetch()

#ifndef NOMPI

      IF((.not. my_process_is_io() .AND.&
           & .not. my_process_is_pref()) .AND.&
           & .not. my_process_is_mpi_test()) THEN

         CALL compute_wait_for_async_pref()
         CALL compute_shutdown_async_pref()

      END IF

      ! deallocating patch data
      DEALLOCATE(patch_data%var_data, patch_data%cells%reorder_index, patch_data%cells%pe_own,     &
        &        patch_data%cells%pe_off, patch_data%edges%reorder_index, patch_data%edges%pe_own, &
        &        patch_data%edges%pe_off)

      ! deallocating intermediate storage latbc_buffer
      DEALLOCATE(latbc_buffer%grp_vars, latbc_buffer%hgrid, latbc_buffer%vars,                     &
        &        latbc_buffer%mapped_name, latbc_buffer%internal_name, latbc_buffer%varID,         &
        &        latbc_buffer%nlev)

      ! clean up global indices data structure.
      CALL latbc_config%global_index%finalize()

#endif
      ! NOMPI

    END SUBROUTINE close_prefetch

    !------------------------------------------------------------------------------------------------
    !> Loop over all output_name_list's, write the ones for which output is due
    !  This routine also cares about opening the output files the first time
    !  and reopening the files after a certain number of steps.
    !
    SUBROUTINE prefetch_input( datetime, p_patch, p_int_state, p_nh_state)
      TYPE(t_datetime), OPTIONAL, INTENT(INOUT) :: datetime
      TYPE(t_patch),          OPTIONAL, INTENT(IN)   :: p_patch
      TYPE(t_int_state),      OPTIONAL, INTENT(IN)   :: p_int_state
      TYPE(t_nh_state),       OPTIONAL, INTENT(INOUT):: p_nh_state  !< nonhydrostatic state on the global domain
      CHARACTER(*), PARAMETER :: method_name = "prefetch_input"

#ifndef NOMPI
      ! Set input prefetch attributes
      IF( my_process_is_work()) THEN
         CALL pref_latbc_data(patch_data, p_patch, p_nh_state, p_int_state, datetime=datetime)
      ELSE IF( my_process_is_pref()) THEN
         CALL pref_latbc_data(patch_data)
      END IF
#endif
      ! NOMPI

    END SUBROUTINE prefetch_input

    !------------------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------------------
    ! The following routines are only needed for asynchronous Input prefetching
    !-------------------------------------------------------------------------------------------------
    !> Main routine for Input Prefetcing PEs.
    !  Please note that this routine never returns.
    !
    SUBROUTINE prefetch_main_proc()
      ! local variables
      LOGICAL                         :: done
      CHARACTER(*), PARAMETER :: method_name = "prefetch_main_proc"

      ! call to initalize the prefetch processor with grid data
      CALL init_prefetch()
      ! Enter prefetch loop
      DO
         ! Wait for a message from the compute PEs to start
         CALL async_pref_wait_for_start(done)
         IF(done) EXIT ! leave loop, we are done
         ! perform input prefetching
         CALL prefetch_input()
         ! Inform compute PEs that we are done
         CALL async_pref_send_handshake()
      END DO
      ! Finalization sequence:
      CALL close_prefetch
      ! Shut down MPI
      CALL stop_mpi

      STOP

    END SUBROUTINE prefetch_main_proc

    !------------------------------------------------------------------------------------------------
    !
    ! Common helper routines for error handling and releasing of resources.
    !
    !------------------------------------------------------------------------------------------------
    !
    !  A simple helper routine to check the result of the last MPI call.
    !
    SUBROUTINE check_mpi_error(routine, mpi_call, mpi_error, l_finish)

      CHARACTER(LEN=*), INTENT(IN)    :: routine, mpi_call
      INTEGER, INTENT(IN)             :: mpi_error
      LOGICAL, INTENT(IN)             :: l_finish

      CHARACTER (LEN=MAX_ERROR_LENGTH):: error_message

#ifndef NOMPI
      IF (mpi_error /= MPI_SUCCESS) THEN
         IF (l_finish) THEN
            WRITE (error_message, '(2a,i5)')TRIM(mpi_call), &
                 &                    ' returned with error=',mpi_error
            CALL finish(routine, TRIM(error_message))
         ELSE
            WRITE (error_message, '(4a,i5)')TRIM(routine), ".", TRIM(mpi_call), &
                 &                    ' returned with error=',mpi_error
            WRITE (nerr, TRIM(error_message))
         ENDIF
      ENDIF
#endif

    END SUBROUTINE check_mpi_error

    !------------------------------------------------------------------------------------------------
    !
    !> Initialize data structures for prefetching boundary data
    !
    !  This routine is called after reading the namelists AND setting up
    !  the domains and variables.
    !
    SUBROUTINE set_patch_data()

#ifndef NOMPI
      ! local variables:
      CHARACTER(LEN=*), PARAMETER   :: routine = modname//"::set_patch_data"

      ! allocate patch data structure
      ! set number of global cells/edges/verts and patch ID

      IF (my_process_is_work() .AND. .NOT. my_process_is_pref()) THEN
         patch_data%nlev          = p_patch(1)%nlev
         patch_data%nlevp1        = p_patch(1)%nlevp1
         patch_data%level         = p_patch(1)%level
         patch_data%nblks_c       = p_patch(1)%nblks_c
         patch_data%nblks_e       = p_patch(1)%nblks_e
         patch_data%n_patch_cells = p_patch(1)%n_patch_cells
         patch_data%n_patch_edges = p_patch(1)%n_patch_edges
         patch_data%n_patch_cells_g = p_patch(1)%n_patch_cells_g
         patch_data%n_patch_edges_g = p_patch(1)%n_patch_edges_g
      END IF

      ! transfer data to prefetching PE
      CALL p_bcast(patch_data%nlev, bcast_root, p_comm_work_2_pref)
      CALL p_bcast(patch_data%nlevp1, bcast_root, p_comm_work_2_pref)
      CALL p_bcast(patch_data%level, bcast_root, p_comm_work_2_pref)
      CALL p_bcast(patch_data%nblks_c, bcast_root, p_comm_work_2_pref)
      CALL p_bcast(patch_data%nblks_e, bcast_root, p_comm_work_2_pref)
      CALL p_bcast(patch_data%n_patch_cells, bcast_root, p_comm_work_2_pref)
      CALL p_bcast(patch_data%n_patch_edges, bcast_root, p_comm_work_2_pref)
      CALL p_bcast(patch_data%n_patch_cells_g, bcast_root, p_comm_work_2_pref)
      CALL p_bcast(patch_data%n_patch_edges_g, bcast_root, p_comm_work_2_pref)

      ! ---------------------------------------------------------------------------
      ! replicate data on prefetch proc
      ! ---------------------------------------------------------------------------
      CALL replicate_data_on_pref_proc()

      IF(.NOT. my_process_is_pref()) THEN

         ! set reorder data on work PE
         CALL set_reorder_data( p_patch(1)%n_patch_cells_g, p_patch(1)%n_patch_cells, &
              p_patch(1)%cells%decomp_info%owner_mask, p_patch(1)%cells%decomp_info%glb_index, &
              patch_data%cells )

         CALL set_reorder_data( p_patch(1)%n_patch_edges_g, p_patch(1)%n_patch_edges, &
              p_patch(1)%edges%decomp_info%owner_mask, p_patch(1)%edges%decomp_info%glb_index, &
              patch_data%edges )

      ENDIF

      IF(.NOT. my_process_is_mpi_test()) THEN
         ! transfer reorder data to prefetch PE
         CALL transfer_reorder_data(patch_data%cells)
         CALL transfer_reorder_data(patch_data%edges)
      ENDIF

      ! subroutine to check whether some variable is specified
      ! in input file and setting flag for its further usage
      IF (latbc_config%itype_latbc == 1) &
           &  CALL check_variable()
#endif

    END SUBROUTINE set_patch_data


    !------------------------------------------------------------------------------------------------
    !> Replicates data (mainly the variable lists) needed for async prefetching
    !  on the prefetching procs.
    SUBROUTINE init_prefetch()

      ! local variables:
      CHARACTER(LEN=*), PARAMETER :: routine = modname//"::init_prefetch"

      ! bcast_root is not used in this case
      bcast_root = 0

#ifndef NOMPI

      ! Set broadcast root for intercommunicator broadcasts
      IF(my_process_is_pref()) THEN
         ! Root is proc 0 on the prefetch PE
         bcast_root = 0
      ELSE
         ! Special root setting for inter-communicators:
         ! The PE really sending must use MPI_ROOT, the others MPI_PROC_NULL
         IF(p_pe_work == 0) THEN
            bcast_root = MPI_ROOT
         ELSE
            bcast_root = MPI_PROC_NULL
         ENDIF
      ENDIF

      ! create and transfer patch data
      CAll set_patch_data()

      ! open and read file containing information of prefetch variables
      CALL read_init_file()

      ! initialize the memory window for communication
      CALL init_remote_memory_window

      ! --- "sparse latbc mode": read only data for boundary rows
      !
      !     this requires index information obtained from an additional
      !     grid file:
      IF (my_process_is_pref() .AND. latbc_config%lsparse_latbc) THEN
        CALL read_netcdf_int_1d(latbc_config%latbc_boundary_grid,                          &
          &                     varname1     = "global_cell_index",                        &
          &                     var1         = latbc_config%global_index%cells,            &
          &                     opt_varname2 = "global_edge_index",                        &
          &                     opt_var2     = latbc_config%global_index%edges,            &
          &                     opt_attname  = "nglobal",                                  &
          &                     opt_attvar1  = latbc_config%global_index%n_patch_cells_g,  &
          &                     opt_attvar2  = latbc_config%global_index%n_patch_edges_g)
        ! consistency checks:
        IF (latbc_config%global_index%n_patch_cells_g /= patch_data%n_patch_cells_g) THEN
          CALL finish(routine, "LatBC boundary cell list does not match in size!")
        END IF
        IF (latbc_config%global_index%n_patch_edges_g /= patch_data%n_patch_edges_g) THEN
          CALL finish(routine, "LatBC boundary edge list does not match in size!")
        END IF
      END IF

      IF( my_process_is_work()) THEN
         ! allocate input data for lateral boundary nudging
         CALL prepare_pref_latbc_data(patch_data, p_patch(1), p_int_state(1), p_nh_state(1), ext_data(1))
      ELSE IF( my_process_is_pref()) THEN
         ! allocate input data for lateral boundary nudging
         CALL prepare_pref_latbc_data(patch_data)
      ENDIF

      CALL message(routine,'Done')
#endif

    END SUBROUTINE init_prefetch


    !-------------------------------------------------------------------------------------------------
    !> open files containing first variable list and analysis
    !
    SUBROUTINE read_init_file()

#ifndef NOMPI
      ! local variables
      CHARACTER(LEN=VARNAME_LEN)              :: grp_name
      CHARACTER(LEN=VARNAME_LEN), ALLOCATABLE :: grp_vars(:)
      ! dictionary which maps prefetch variable names onto
      ! GRIB2 shortnames or NetCDF var names.
      TYPE (t_dictionary) :: latbc_varnames_dict
      CHARACTER(LEN=filename_max) :: latbc_file
      CHARACTER(*), PARAMETER :: routine = "mo_async_latbc::read_init_files"
      CHARACTER (len=MAX_CHAR_LENGTH) :: name
      CHARACTER(len=132) :: message_text
      INTEGER :: jlev, ierrstat, vlistID, nvars, varID, zaxisID, gridID, &
           &       jp, fileID_latbc, counter, filetype, ngrp_prefetch_vars
      INTEGER(KIND=i8) :: flen_latbc
      LOGICAL :: l_exist
      CHARACTER(LEN=filename_max)    :: latbc_filename
      CHARACTER(LEN=MAX_CHAR_LENGTH) :: cdiErrorText

      ! allocating buffers containing name of variables
      ALLOCATE(latbc_buffer%grp_vars(200))
      ALLOCATE(latbc_buffer%mapped_name(200))
      ALLOCATE(latbc_buffer%internal_name(200))
      ALLOCATE(grp_vars(200))
      ALLOCATE(StrLowCasegrp(200))

      ! initialising counter
      counter = 0

      !>Looks for variable groups ("group:xyz") and collects
      ! them to map prefetch variable names onto
      ! GRIB2 shortnames or NetCDF var names.
      !
      ! variables in this group are ICON, COSMO or IFS
      ! data which are read by the prefetch PE:
      grp_name = 'LATBC_PREFETCH_VARS'

      ! loop over all variables and collects the variables names
      ! corresponding to the group "grp_name"

      CALL collect_group( TRIM(grp_name), grp_vars, ngrp_prefetch_vars, loutputvars_only=.FALSE., &
           &                lremap_lonlat=.FALSE. )

      ! adding the variable 'GEOSP' to the list by add_to_list
      ! as the variable cannot be found in metadata variable list
      IF (latbc_config%itype_latbc == 1) &
           CALL add_to_list( grp_vars, ngrp_prefetch_vars, (/latbc_buffer%geop_ml_var/) , 1)

      ! read the map file into dictionary data structure
      CALL dict_init(latbc_varnames_dict, lcase_sensitive=.FALSE.)

      IF(LEN_TRIM(latbc_config%latbc_varnames_map_file) > 0) THEN
         CALL dict_loadfile(latbc_varnames_dict, TRIM(latbc_config%latbc_varnames_map_file))
      END IF

      ! allocate the number of vertical levels with the
      ! same size as number of variables
      ALLOCATE(latbc_buffer%nlev(ngrp_prefetch_vars), STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish(routine, ALLOCATE_FAILED)

      ! allocate the array for variable ID
      ! with same size as number of variables
      ALLOCATE(latbc_buffer%varID(ngrp_prefetch_vars), STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish(routine, ALLOCATE_FAILED)

      IF(my_process_is_work() .AND.  p_pe_work == p_work_pe0) THEN !!!!!!!use prefetch processor here
         jlev = patch_data%level
         ! generate file name
         latbc_filename = generate_filename(nroot, jlev, time_config%ini_datetime)
         latbc_file = TRIM(latbc_config%latbc_path)//TRIM(latbc_filename)
         INQUIRE (FILE=latbc_file, EXIST=l_exist)
         IF (.NOT.l_exist) THEN
            CALL finish(TRIM(routine),'DWD FG file not found: '//TRIM(latbc_file))
         ENDIF

         ! open file
         !
         fileID_latbc = streamOpenRead(TRIM(latbc_file))
         IF (fileID_latbc < 0) THEN
           CALL cdiGetStringError(fileID_latbc, cdiErrorText)
           WRITE(message_text,'(4a)') 'File ', TRIM(latbc_file), &
                ' cannot be opened: ', TRIM(cdiErrorText)
           CALL finish(routine, TRIM(message_text))
         ENDIF

         filetype = streamInqFiletype(fileID_latbc)

         ! Search name mapping for name in GRIB2 file
         SELECT CASE(filetype)
         CASE (FILETYPE_NC2, FILETYPE_NC4)
            IF (latbc_config%itype_latbc == 1) THEN
               DO jp= 1, ngrp_prefetch_vars
                  latbc_buffer%grp_vars(jp) = TRIM(dict_get(latbc_varnames_dict, grp_vars(jp), default=grp_vars(jp)))
               ENDDO
            ELSE
               DO jp= 1, ngrp_prefetch_vars
                  latbc_buffer%grp_vars(jp) = TRIM(grp_vars(jp))
               ENDDO
            ENDIF
         CASE (FILETYPE_GRB2)
            IF (latbc_config%itype_latbc == 1) THEN
               DO jp= 1, ngrp_prefetch_vars
                  latbc_buffer%grp_vars(jp) = TRIM(dict_get(latbc_varnames_dict, grp_vars(jp), default=grp_vars(jp)))
               ENDDO
            ELSE
               DO jp= 1, ngrp_prefetch_vars
                  latbc_buffer%grp_vars(jp) = TRIM(grp_vars(jp))
               ENDDO
            ENDIF
         CASE DEFAULT
            CALL finish(routine, "Unknown file type")
         END SELECT

         !     WRITE(0,*) 'latbc_buffer%grp_name ',  latbc_buffer%grp_vars(1:ngrp_prefetch_vars), 'ngrp_vars ', ngrp_prefetch_vars

         ! check whether the file is empty (does not work unfortunately; internal CDI error)
         flen_latbc = util_filesize(TRIM(latbc_file))
         IF (flen_latbc <= 0 ) THEN
            WRITE(message_text,'(a)') 'File '//TRIM(latbc_file)//' is empty'
            CALL message(TRIM(routine), TRIM(message_text))
            CALL finish(routine, "STOP: Empty input file")
         ENDIF

         vlistID = streamInqVlist(fileID_latbc)

         ! get the number of variables
         nvars = vlistNvars(vlistID)

         ! get the number of vertical levels for the
         ! required prefetch variables
         LOOP : DO varID=0,(nvars-1)
            CALL vlistInqVarName(vlistID, varID, name)
            !    WRITE(0,*) 'name ', name
            DO jp = 1, ngrp_prefetch_vars !latbc_buffer%ngrp_vars
               IF(tolower(name) == tolower(latbc_buffer%grp_vars(jp))) THEN
                  ! get the vertical axis ID
                  zaxisID = vlistInqVarZaxis(vlistID, varID)
                  ! get the grid ID using vlistID and varID
                  gridID = vlistInqVarGrid(vlistID, varID)
                  counter = counter + 1
                  ! get the respective vertical levels for
                  ! the respective variable
                  latbc_buffer%nlev(counter) = zaxisInqSize(zaxisID)
                  !variable ID for variable to be read from file
                  latbc_buffer%varID(counter) = varID

                  ! getting the variable name in the file
                  latbc_buffer%mapped_name(counter) = TRIM(name)
                  latbc_buffer%internal_name(counter) = &
                    TRIM(dict_get(latbc_varnames_dict, TRIM(name), linverse=.TRUE., default=TRIM(name)))
                  ! getting the variable name in lower case letter
                  StrLowCasegrp(counter) = TRIM(grp_vars(jp))

                  !       WRITE(0,*) 'mapped_name ',  latbc_buffer%mapped_name(counter)
               ENDIF
            ENDDO
         END DO LOOP

         ! closes the open dataset
         CALL streamClose(fileID_latbc)
      END IF

      CALL p_bcast(latbc_buffer%nlev(:), p_comm_work_pref_compute_pe0, p_comm_work_pref)
      CALL p_bcast(latbc_buffer%varID(:), p_comm_work_pref_compute_pe0, p_comm_work_pref)
      CALL p_bcast(latbc_buffer%mapped_name(:), p_comm_work_pref_compute_pe0, p_comm_work_pref)
      CALL p_bcast(latbc_buffer%internal_name(:), p_comm_work_pref_compute_pe0, p_comm_work_pref)
      CALL p_bcast(StrLowCasegrp(:), p_comm_work_pref_compute_pe0, p_comm_work_pref)
      CALL p_bcast(counter, p_comm_work_pref_compute_pe0, p_comm_work_pref)

      !WRITE(0,*) 'mapped_name ',  latbc_buffer%mapped_name(1:counter), 'ngrp_vars ', ngrp_prefetch_vars

      ! getting the count of number of variables in
      ! the file to be read
      latbc_buffer%ngrp_vars = counter

      ! destroy variable name dictionaries:
      CALL dict_finalize(latbc_varnames_dict)

      ! allocate the number of buffer sizes for variables
      ! with same size as number of variables
      ALLOCATE(latbc_buffer%vars(latbc_buffer%ngrp_vars), STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish(routine, ALLOCATE_FAILED)
#endif

    END SUBROUTINE read_init_file


    !-------------------------------------------------------------------------------------------------
    !> open file to determine a variable or its alternative variable is provided as input
    ! and setting flag for its further usage

    SUBROUTINE check_variable()

#ifndef NOMPI
      ! local variables
      CHARACTER(LEN=filename_max) :: latbc_file
      CHARACTER(*), PARAMETER :: routine = "mo_async_latbc::read_init_files"
      INTEGER :: jlev, fileID_latbc
      LOGICAL :: l_exist
      CHARACTER(len=132) :: message_text
      CHARACTER(LEN=filename_max)    :: latbc_filename
      CHARACTER(LEN=MAX_CHAR_LENGTH) :: cdiErrorText

      ! prefetch processor opens the file and checks if variables are present
      IF( my_process_is_work() .AND.  p_pe_work == p_work_pe0) THEN !!!!!!!use prefetch processor here
         jlev = patch_data%level
         ! generate file name
         latbc_filename = generate_filename(nroot, jlev, time_config%ini_datetime)
         latbc_file = TRIM(latbc_config%latbc_path)//TRIM(latbc_filename)
         INQUIRE (FILE=latbc_file, EXIST=l_exist)
         IF (.NOT.l_exist) THEN
            CALL finish(TRIM(routine),'DWD FG file not found: '//TRIM(latbc_file))
         ENDIF

         ! open file
         !
         fileID_latbc = streamOpenRead(TRIM(latbc_file))
         IF (fileID_latbc < 0) THEN
           CALL cdiGetStringError(fileID_latbc, cdiErrorText)
           WRITE(message_text,'(4a)') 'File ', TRIM(latbc_file), &
                ' cannot be opened: ', TRIM(cdiErrorText)
           CALL finish(routine, TRIM(message_text))
         ENDIF

         !
         ! Check if the prognostic thermodynamic variables (rho and theta_v) are provided as input
         !
         IF ((test_cdi_varID(fileID_latbc, 'RHO') /= -1 .OR. test_cdi_varID(fileID_latbc, 'DEN') /= -1) &
             .AND. test_cdi_varID(fileID_latbc, 'THETA_V') /= -1) THEN
           latbc_buffer%lthd_progvars = .true.
           CALL message(TRIM(routine),'Prognostic thermodynamic variables (rho and theta_v) are used')
         ELSE
           latbc_buffer%lthd_progvars = .false.
         ENDIF


         !
         ! Check if surface pressure (PS) or its logarithm (LNPS) is provided as input
         !
         IF (test_cdi_varID(fileID_latbc, 'PS') /= -1) THEN
            latbc_buffer%psvar = 'PS'
         ELSE IF (test_cdi_varID(fileID_latbc, 'LNPS') /= -1) THEN
            latbc_buffer%psvar = 'LNPS'
         ENDIF

         !
         ! Check if model-level surface Geopotential is provided as GEOSP or GEOP_ML
         !
         IF (test_cdi_varID(fileID_latbc, 'GEOSP') /= -1) THEN
            latbc_buffer%geop_ml_var = 'GEOSP'
         ELSE IF (test_cdi_varID(fileID_latbc, 'GEOP_ML') /= -1) THEN
            latbc_buffer%geop_ml_var = 'GEOP_ML'
         ELSE IF (.NOT. (latbc_buffer%lthd_progvars .OR. test_cdi_varID(fileID_latbc, 'HHL') /= -1 &
                  .OR. test_cdi_varID(fileID_latbc, 'Z_IFC') /= -1) ) THEN
            CALL finish(TRIM(routine),'Could not find model-level sfc geopotential')
         ENDIF

         !
         ! Check if rain water (QR) is provided as input
         !
         IF (test_cdi_varID(fileID_latbc, 'QR') /= -1) THEN
            latbc_buffer%lread_qr = .true.
         ELSE
            latbc_buffer%lread_qr = .false.
            CALL message(TRIM(routine),'Rain water (QR) not available in input data')
         ENDIF

         !
         ! Check if snow water (QS) is provided as input
         !
         IF (test_cdi_varID(fileID_latbc, 'QS') /= -1) THEN
            latbc_buffer%lread_qs = .true.
         ELSE
            latbc_buffer%lread_qs = .false.
            CALL message(TRIM(routine),'Snow water (QS) not available in input data')
         ENDIF

         !
         ! Check if surface pressure (VN) is provided as input
         !
         IF (test_cdi_varID(fileID_latbc, 'VN') /= -1) THEN
            latbc_buffer%lread_vn = .TRUE.
         ELSE
            latbc_buffer%lread_vn = .FALSE.
         ENDIF

         ! closes the open dataset
         CALL streamClose(fileID_latbc)

      ENDIF

      ! broadcast data to prefetching and compute PE's
      ! public constant: p_comm_work_pref_compute_pe0
      CALL p_bcast(latbc_buffer%psvar, p_comm_work_pref_compute_pe0, p_comm_work_pref)
      CALL p_bcast(latbc_buffer%geop_ml_var, p_comm_work_pref_compute_pe0, p_comm_work_pref)
      CALL p_bcast(latbc_buffer%lread_qs, p_comm_work_pref_compute_pe0, p_comm_work_pref)
      CALL p_bcast(latbc_buffer%lread_qr, p_comm_work_pref_compute_pe0, p_comm_work_pref)
      CALL p_bcast(latbc_buffer%lread_vn, p_comm_work_pref_compute_pe0, p_comm_work_pref)
      CALL p_bcast(latbc_buffer%lthd_progvars, p_comm_work_pref_compute_pe0, p_comm_work_pref)
#endif

    END SUBROUTINE check_variable


    !-------------------------------------------------------------------------------------------------
    !> Replicates data needed for async prefetch on the prefetch proc.
    !  ATTENTION: The data is not completely replicated, only as far as needed for prefetch.
    !
    !  This routine has to be called by all PEs (work and prefetch)
    !
    SUBROUTINE replicate_data_on_pref_proc()

#ifndef NOMPI
      ! local variables
      CHARACTER(len=*), PARAMETER   :: routine = modname//"::replicate_data_on_pref_proc"
      INTEGER                       :: info_size, i, iv, nelems, nv, n, list_info, all_var, &
           &                           all_nvars, nvars, i2
      INTEGER, ALLOCATABLE          :: info_storage(:,:)
      TYPE(t_list_element), POINTER :: element
      TYPE(t_var_metadata)          :: info
      TYPE(t_var_list)              :: p_var_list
      ! var_list_name should have at least the length of var_list names
      CHARACTER(LEN=256)            :: var_list_name

      ! There is nothing to do for the test PE:
      IF(my_process_is_mpi_test()) RETURN

      ! get the size - in default INTEGER words - which is needed to
      ! hold the contents of TYPE(t_var_metadata)
      info_size = SIZE(TRANSFER(info, (/ 0 /)))

      ! get the number of var lists
      IF(.NOT. my_process_is_pref()) nv = nvar_lists
      CALL p_bcast(nv, bcast_root, p_comm_work_2_pref)

      IF(.NOT.my_process_is_pref()) THEN
         all_nvars = 0
         DO i = 1, nvar_lists

            ! Count the number of variable entries
            element => var_lists(i)%p%first_list_element
            !   IF(element%field%info%used_dimensions(2) == 0) CYCLE
            nvars = 0
            DO
               IF(.NOT.ASSOCIATED(element)) EXIT
               !      IF(element%field%info%used_dimensions(2) == 0) CYCLE
               nvars = nvars+1
               element => element%next_list_element
            ENDDO
            all_nvars = all_nvars + nvars
         ENDDO
      ENDIF

      ! get the number of var lists
      IF(.NOT. my_process_is_pref()) all_var = all_nvars
      CALL p_bcast(all_var, bcast_root, p_comm_work_2_pref)

      IF (all_var <= 0) RETURN

      ! allocate the array of variables
      ALLOCATE(patch_data%var_data(all_var))

      i2 = 0
      ! For each var list, get its components
      DO iv = 1, nv

         ! Send name
         IF(.NOT.my_process_is_pref()) var_list_name = var_lists(iv)%p%name
         CALL p_bcast(var_list_name, bcast_root, p_comm_work_2_pref)

         IF(.NOT.my_process_is_pref()) THEN
            ! Count the number of variable entries
            element => var_lists(iv)%p%first_list_element
            nelems = 0
            DO
               IF(.NOT.ASSOCIATED(element)) EXIT
               nelems = nelems+1
               element => element%next_list_element
            ENDDO
            list_info = nelems
         ENDIF

         ! Send basic info:
         CALL p_bcast(list_info, bcast_root, p_comm_work_2_pref)

         IF(my_process_is_pref()) THEN
            nelems = list_info
            ! Create var list
            CALL new_var_list( p_var_list, var_list_name)
         ENDIF

         ! Get the binary representation of all info members of the variables
         ! of the list and send it to the receiver.
         ! Using the Fortran TRANSFER intrinsic may seem like a hack,
         ! but it has the advantage that it is completely independet of the
         ! actual declaration if TYPE(t_var_metadata).
         ! Thus members may added to or removed from TYPE(t_var_metadata)
         ! without affecting the code below and we don't have an additional
         ! cross dependency between TYPE(t_var_metadata) and this module.

         ALLOCATE(info_storage(info_size, nelems))

         IF(.NOT.my_process_is_pref()) THEN
            element => var_lists(iv)%p%first_list_element
            nelems = 0
            DO
               IF(.NOT.ASSOCIATED(element)) EXIT
               i2 = i2 + 1
               patch_data%var_data(i2)%info = element%field%info
               nelems = nelems+1
               info_storage(:,nelems) = TRANSFER(element%field%info, (/ 0 /))
               element => element%next_list_element
            ENDDO
         ENDIF

         ! Send binary representation of all info members

         CALL p_bcast(info_storage, bcast_root, p_comm_work_2_pref)

         IF(my_process_is_pref()) THEN

            ! Insert elements into var list
            p_var_list%p%first_list_element => NULL()
            element => NULL() ! Safety only
            DO n = 1, nelems
               IF(.NOT.ASSOCIATED(p_var_list%p%first_list_element)) THEN
                  ALLOCATE(p_var_list%p%first_list_element)
                  element => p_var_list%p%first_list_element
               ELSE
                  ALLOCATE(element%next_list_element)
                  element => element%next_list_element
               ENDIF

               i2 = i2 + 1
               element%next_list_element => NULL()

               ! Set info structure from binary representation in info_storage
               element%field%info = TRANSFER(info_storage(:, n), info)
               patch_data%var_data(i2)%info = element%field%info
            ENDDO
         ENDIF
         DEALLOCATE(info_storage)
      ENDDO
#endif

    END SUBROUTINE replicate_data_on_pref_proc


    !------------------------------------------------------------------------------------------------
    !> Sets the reorder_data for cells/edges/verts
    !  ATTENTION: This routine must only be called on work and test PE
    !             (i.e. not on prefetching PEs)
    !             The arguments don't make sense on the prefetching PE anyways
    !
    SUBROUTINE set_reorder_data(n_points_g, n_points, owner_mask, glb_index, p_reo)

      INTEGER, INTENT(IN) :: n_points_g      ! Global number of cells/edges/verts in logical patch
      INTEGER, INTENT(IN) :: n_points        ! Local number of cells/edges/verts in logical patch
      LOGICAL, INTENT(IN) :: owner_mask(:,:) ! owner_mask for logical patch
      INTEGER, INTENT(IN) :: glb_index(:)    ! glb_index for logical patch
      TYPE(t_reorder_data), INTENT(INOUT) :: p_reo ! Result: reorder info

#ifndef NOMPI
      ! local variables
      INTEGER :: i, n, il, ib, mpi_error, ierrstat
      LOGICAL, ALLOCATABLE :: phys_owner_mask(:) ! owner mask for physical patch
      INTEGER, ALLOCATABLE :: glbidx_own(:), glbidx_glb(:), reorder_index_log_dom(:)

      CHARACTER (LEN=MAX_ERROR_LENGTH):: error_message
      CHARACTER(LEN=*), PARAMETER :: routine = modname//"::set_reorder_data"

      ! Just for safety
      IF(my_process_is_pref()) CALL finish(routine, 'Must not be called on Prefetching PE')

      ! Set the physical patch owner mask
      ALLOCATE(phys_owner_mask(n_points))
      DO i = 1, n_points
         il = idx_no(i)
         ib = blk_no(i)
         phys_owner_mask(i) = owner_mask(il,ib)
      ENDDO

      ! Get number of owned cells/edges/verts (without halos, physical patch only)
      p_reo%n_own = COUNT(phys_owner_mask(:))

      !   WRITE(*,*) 'set_reorder_data p_pe_work ', p_pe_work , ' p_reo%n_own ', p_reo%n_own, ' n_points ', n_points

      ! Set index arrays to own cells/edges/verts
      ALLOCATE(p_reo%own_idx(p_reo%n_own), STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish(routine, ALLOCATE_FAILED)
      ALLOCATE(p_reo%own_blk(p_reo%n_own), STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish(routine, ALLOCATE_FAILED)

      ! Global index of my own points
      ALLOCATE(glbidx_own(p_reo%n_own), STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish(routine, ALLOCATE_FAILED)

      n = 0
      DO i = 1, n_points
         IF(phys_owner_mask(i)) THEN
            n = n+1
            p_reo%own_idx(n) = idx_no(i)
            p_reo%own_blk(n) = blk_no(i)
            glbidx_own(n)  = glb_index(i)
         ENDIF
      ENDDO

      ! Gather the number of own points for every PE into p_reo%pe_own
      ALLOCATE(p_reo%pe_own(0:p_n_work-1), STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish(routine, ALLOCATE_FAILED)
      ALLOCATE(p_reo%pe_off(0:p_n_work-1), STAT=ierrstat)
      IF (ierrstat  /= SUCCESS) CALL finish(routine, ALLOCATE_FAILED)

      ! make sure all are done
      ! CALL p_barrier(comm=p_comm_work)

      ! Gather the number of own points for every PE into p_reo%pe_own
      CALL MPI_Allgather(p_reo%n_own,  1, p_int, &
           p_reo%pe_own, 1, p_int, &
           p_comm_work, mpi_error)

      CALL check_mpi_error(routine, 'MPI_Allgather', mpi_error, .TRUE.)

      ! Get offset within result array
      p_reo%pe_off(0) = 0
      DO i = 1, p_n_work-1
         p_reo%pe_off(i) = p_reo%pe_off(i-1) + p_reo%pe_own(i-1)
      ENDDO

      ! Get global number of points for current (physical!) patch
      p_reo%n_glb = SUM(p_reo%pe_own(:))

      ! Get the global index numbers of the data when it is gathered on PE 0
      ! exactly in the same order as it is retrieved later during prefetching.
      ALLOCATE(glbidx_glb(p_reo%n_glb), STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish(routine, ALLOCATE_FAILED)

      ! make sure all are done
      ! CALL p_barrier(comm=p_comm_work)

      CALL MPI_Allgatherv(glbidx_own, p_reo%n_own, p_int, &
           glbidx_glb, p_reo%pe_own, p_reo%pe_off, p_int, &
           p_comm_work, mpi_error)

      CALL check_mpi_error(routine, 'MPI_Allgatherv', mpi_error, .TRUE.)

      ! Get reorder_index
      ALLOCATE(p_reo%reorder_index(p_reo%n_glb), STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish(routine, ALLOCATE_FAILED)

      ! Spans the complete logical domain
      ALLOCATE(reorder_index_log_dom(n_points_g), STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish(routine, ALLOCATE_FAILED)
      reorder_index_log_dom(:) = 0

      DO i = 1, p_reo%n_glb
         ! Reorder_index_log_dom stores where a global point in logical domain comes from.
         ! It is nonzero only at the physical patch locations
         reorder_index_log_dom(glbidx_glb(i)) = i
      ENDDO

      ! Gather the reorder index for the physical domain
      n = 0
      DO i = 1, n_points_g
         IF(reorder_index_log_dom(i)>0) THEN
            n = n+1
            p_reo%reorder_index(n) = reorder_index_log_dom(i)
         ENDIF
      ENDDO

      ! Safety check
      IF(n/=p_reo%n_glb) THEN
         WRITE (error_message, '(a,i8,a,i8)') 'Reordering failed: n=',n, &
              & ' /= p_reo%n_glb=',p_reo%n_glb
         CALL finish(routine,TRIM(error_message))
      ENDIF

      DEALLOCATE(phys_owner_mask)
      DEALLOCATE(glbidx_own)
      DEALLOCATE(glbidx_glb)
      DEALLOCATE(reorder_index_log_dom)
#endif

    END SUBROUTINE set_reorder_data


    !------------------------------------------------------------------------------------------------
    !
    ! Transfers reorder data to restart PEs.
    !
    SUBROUTINE transfer_reorder_data(p_reo)
      TYPE(t_reorder_data), INTENT(INOUT) :: p_reo

#ifndef NOMPI
      ! local variables
      INTEGER                             :: ierrstat
      CHARACTER(LEN=*), PARAMETER :: routine = modname//"::transfer_reorder_data"

      ! There is nothing to do for the test PE:
      IF(my_process_is_mpi_test()) RETURN

      !transfer the global number of points, this is not known on prefetching PE
      CALL p_bcast(p_reo%n_glb, bcast_root, p_comm_work_2_pref)

      IF(my_process_is_pref())THEN

         ! on prefetch PE: n_own=0, own_ide and own_blk are not allocated
         p_reo%n_own = 0

         ! pe_own must be allocated for num_work_procs, not for p_n_work
         ALLOCATE(p_reo%pe_own(0:num_work_procs-1), STAT=ierrstat)
         IF(ierrstat /= SUCCESS) CALL finish(routine, ALLOCATE_FAILED)

         ALLOCATE(p_reo%pe_off(0:num_work_procs-1), STAT=ierrstat)
         IF(ierrstat /= SUCCESS) CALL finish(routine, ALLOCATE_FAILED)

         ALLOCATE(p_reo%reorder_index(p_reo%n_glb), STAT=ierrstat)
         IF(ierrstat /= SUCCESS) CALL finish(routine, ALLOCATE_FAILED)
      ENDIF

      CALL p_bcast(p_reo%pe_own, bcast_root, p_comm_work_2_pref)
      CALL p_bcast(p_reo%pe_off, bcast_root, p_comm_work_2_pref)
      CALL p_bcast(p_reo%reorder_index, bcast_root, p_comm_work_2_pref)
#endif
    END SUBROUTINE  transfer_reorder_data


    !------------------------------------------------------------------------------------------------
    !> Initializes the remote memory window for asynchronous prefetch.
    !
    SUBROUTINE init_remote_memory_window

#ifndef NOMPI
      INTEGER :: ierrstat, iv, jp, nlevs, hgrid
      INTEGER (KIND=MPI_ADDRESS_KIND) :: mem_size
      LOGICAL ,ALLOCATABLE :: grp_vars_bool(:)
      CHARACTER(LEN=*), PARAMETER :: routine = modname//"::init_memory_window"

      ! There is nothing to do for the test PE:
      IF(my_process_is_mpi_test()) RETURN

      patch_data%mem_win%mpi_win = MPI_WIN_NULL

      ! Get size and offset of the data for the input
      mem_size = 0_i8

      ALLOCATE(grp_vars_bool(latbc_buffer%ngrp_vars))
      ALLOCATE(latbc_buffer%hgrid(latbc_buffer%ngrp_vars))

      DO jp = 1, latbc_buffer%ngrp_vars
         grp_vars_bool(jp) = .FALSE.
      ENDDO

      ! Go over all input variables
      DO iv = 1, SIZE(patch_data%var_data)
         DO jp = 1, latbc_buffer%ngrp_vars
            ! Use only the variables of time level 1 (".TL1") to determine memory sizes.
            IF((TRIM(StrLowCasegrp(jp)) == TRIM(patch_data%var_data(iv)%info%name)) .OR. &
                 & (TRIM(StrLowCasegrp(jp))//'.TL1' == TRIM(patch_data%var_data(iv)%info%name))) THEN

               nlevs = 0
               IF(.NOT. grp_vars_bool(jp))  THEN
                  IF (patch_data%var_data(iv)%info%ndims == 2) THEN
                     nlevs = 1
                  ELSE
                     nlevs = latbc_buffer%nlev(jp) !latbc_config%nlev_in
                  ENDIF

                  IF (nlevs == 0) CYCLE

                  SELECT CASE (patch_data%var_data(iv)%info%hgrid)

                  CASE (GRID_UNSTRUCTURED_CELL)
                     mem_size = mem_size + INT(nlevs*patch_data%cells%n_own,i8)

                     IF(my_process_is_work())THEN
                        ! allocate the buffer sizes for variables on compute processors
                        ALLOCATE(latbc_buffer%vars(jp)%buffer(nproma, nlevs, patch_data%nblks_c), STAT=ierrstat)
                        IF (ierrstat /= SUCCESS) CALL finish(routine, ALLOCATE_FAILED)
                     ENDIF

                     ! variable stored in cell center location
                     latbc_buffer%hgrid(jp) = patch_data%var_data(iv)%info%hgrid
                     hgrid = patch_data%var_data(iv)%info%hgrid

                  CASE (GRID_UNSTRUCTURED_EDGE)
                     mem_size = mem_size + INT(nlevs*patch_data%edges%n_own,i8)

                     IF(my_process_is_work())THEN
                        ! allocate the buffer sizes for variables on compute processors
                        ALLOCATE(latbc_buffer%vars(jp)%buffer(nproma, nlevs, patch_data%nblks_e), STAT=ierrstat)
                        IF (ierrstat /= SUCCESS) CALL finish(routine, ALLOCATE_FAILED)
                     ENDIF

                     ! variable stored in edge center of a cell
                     latbc_buffer%hgrid(jp) = patch_data%var_data(iv)%info%hgrid

                  CASE DEFAULT
                     CALL finish(routine,UNKNOWN_GRID_TYPE)
                  END SELECT

                  grp_vars_bool(jp) = .TRUE.

               ENDIF
            ENDIF
         ENDDO
      ENDDO ! vars

      IF (latbc_config%itype_latbc == 1) THEN
         DO jp = 1, latbc_buffer%ngrp_vars
            IF (TRIM(latbc_buffer%mapped_name(jp)) == TRIM(latbc_buffer%geop_ml_var)) THEN
               ! Memory for GEOSP variable taken as memory equivalent to 1 level of z_ifc
               ! as the variable GEOSP doesnt exist in metadata
               mem_size = mem_size + INT(1*patch_data%cells%n_own,i8)

               IF(my_process_is_work())THEN
                  ! allocate the buffer sizes for variable 'GEOSP' on compute processors
                  ALLOCATE(latbc_buffer%vars(jp)%buffer(nproma, 1, patch_data%nblks_c), STAT=ierrstat)
                  IF (ierrstat /= SUCCESS) CALL finish(routine, ALLOCATE_FAILED)
               ENDIF

               ! variable GEOSP is stored in cell center location
               latbc_buffer%hgrid(jp) = GRID_UNSTRUCTURED_CELL
            ENDIF
         ENDDO
      ENDIF

      DEALLOCATE(StrLowCasegrp)

      ! allocate amount of memory needed with MPI_Alloc_mem
      CALL allocate_mem_noncray(mem_size)
#endif

    END SUBROUTINE init_remote_memory_window


    !------------------------------------------------------------------------------------------------
    !> allocate amount of memory needed with MPI_Alloc_mem
    !
    !  @note Implementation for non-Cray pointers
    !
    SUBROUTINE allocate_mem_noncray(mem_size)

#ifdef NOMPI
      INTEGER, INTENT(IN)    :: mem_size
#else
      INTEGER (KIND=MPI_ADDRESS_KIND), INTENT(IN)    :: mem_size

      ! local variables
      CHARACTER(LEN=*), PARAMETER :: routine = modname//"::allocate_mem_noncray"
      TYPE(c_ptr)                     :: c_mem_ptr
      INTEGER                         :: mpierr
      INTEGER                         :: nbytes_real
      INTEGER (KIND=MPI_ADDRESS_KIND) :: mem_bytes

      ! Get the amount of bytes per REAL*4 variable (as used in MPI
      ! communication)
      CALL MPI_Type_extent(p_real_sp, nbytes_real, mpierr)

      ! For the IO PEs the amount of memory needed is 0 - allocate at least 1 word there:
      mem_bytes = MAX(mem_size,1_i8)*INT(nbytes_real,i8)

      ! TYPE(c_ptr) and INTEGER(KIND=MPI_ADDRESS_KIND) do NOT necessarily have the same size!!!
      ! So check if at least c_intptr_t and MPI_ADDRESS_KIND are the same, else we may get
      ! into deep, deep troubles!
      ! There is still a slight probability that TYPE(c_ptr) does not have the size indicated
      ! by c_intptr_t since the standard only requires c_intptr_t is big enough to hold pointers
      ! (so it may be bigger than a pointer), but I hope no vendor screws up its ISO_C_BINDING
      ! in such a way!!!
      ! If c_intptr_t<=0, this type is not defined and we can't do this check, of course.

      IF(c_intptr_t > 0 .AND. c_intptr_t /= MPI_ADDRESS_KIND) &
           & CALL finish(routine,'c_intptr_t /= MPI_ADDRESS_KIND, too dangerous to proceed!')

      CALL MPI_Alloc_mem(mem_bytes, MPI_INFO_NULL, c_mem_ptr, mpierr)

      ! The NEC requires a standard INTEGER array as 3rd argument for c_f_pointer,
      ! although it would make more sense to have it of size MPI_ADDRESS_KIND.

      NULLIFY(patch_data%mem_win%mem_ptr_sp)

#ifdef __SX__
      CALL C_F_POINTER(c_mem_ptr, patch_data%mem_win%mem_ptr_sp, (/ INT(mem_size) /) )
#else
      CALL C_F_POINTER(c_mem_ptr, patch_data%mem_win%mem_ptr_sp, (/ mem_size /) )
#endif
      ! Create memory window for communication
      patch_data%mem_win%mem_ptr_sp(:) = 0._sp
      CALL MPI_Win_create( patch_data%mem_win%mem_ptr_sp, mem_bytes, nbytes_real, MPI_INFO_NULL,&
           &                  p_comm_work_pref, patch_data%mem_win%mpi_win, mpierr )

      IF (mpierr /= 0) CALL finish(TRIM(routine), "MPI error!")
#endif

    END SUBROUTINE allocate_mem_noncray

END MODULE mo_async_latbc
