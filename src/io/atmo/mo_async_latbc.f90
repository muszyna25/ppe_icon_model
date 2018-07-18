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
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
!!
!! ------------------------------------------------------------------------
!! Which fields are read from the lateral boundary conditions file?
!! ------------------------------------------------------------------------
!!
!! This question is answered independently from the "init_icon"
!! namelist parameter of the initial state setup!
!!
!! - If "VN" is available, then it is read from file, otherwise "U","V".
!! - "W" is optional and read if available in the input data (note that "W" may in fact contain OMEGA).
!! - "QV", "QC", "QI" are always read
!! - "QR", "QS" are read if available
!!
!! The other fields for the lateral boundary conditions are read
!! from file, based on the following decision tree.  Note that the
!! basic distinction between input from non-hydrostatic and
!! hydrostatic models is made on the availability of the HHL field.
!!
!!                            +------------------+
!!                            |  HHL available?  |
!!                            +------------------+
!!                                     |
!!                     ______yes_______|________no________
!!                     |                                  |
!!         +--------------------------+            +------------------------+
!!         | RHO & THETA_V available? |            | PS,GEOP & T available? |
!!         +--------------------------+            +------------------------+
!!                     |                                        |
!!           ____yes___|____no______                    ___yes__|________no___________
!!           |                      |                  |                              |
!!           |               +----------------+        |                          +--------+
!! * read HHL,RHO,THETA_V    | P,T available? |     * read in PS,GEOP,T           | ERROR! | 
!! * read W if available     |                |     * read OMEGA if available     |        |
!! * ignore PS,GEOP          +----------------+     * compute P,HHL               +--------+
!! * diagnose P,T                   |               * CALL OMEGA -> W
!!                          ___yes__|___no____
!!                         |                  |
!!                         |               +--------+
!!                    * read HHL,P,T       | ERROR! |
!!                    * ignore PS,GEOP     +--------+
!!                    * read W if available
!!
!!
!! Afterwards, we 
!! - re-compute the virtual temperature (inside the vertical
!!   interpolation subroutine)
!! - (re-)compute RHO (inside the vertical interpolation subroutine)
!! - perform vertical interpolation
!!
!!
!!
!! ------------------------------------------------------------------------
!! Read-in of lateral boundary data: General Overview of the Implementation
!! ------------------------------------------------------------------------
!! 
!! Note: This short documentation focuses on the "asynchronous
!! prefetching" mode only, the old (possibly deprecated) synchronous
!! read-in of the boundary data, "mo_sync_latbc.f90", is not covered.
!! 
!! Read-in of boundary data is invoked via
!!   CALL recv_latbc_data
!! in the time loop (module "mo_nh_stepping").
!! 
!! 
!! Modules related to (asynchronous) boundary data read-in:
!! --------------------------------------------------------
!! 
!! src/io/atmo/mo_async_latbc.f90              : Setup of the boundary data read-in functionality
!! src/io/atmo/mo_async_utils.f90              : * Initialisation, allocation
!!                                               * Top level routines: "read-in" (e.g. "recv_latbc_data")
!!                                               * Top level routines: "fetch" (see explanation below)
!! src/io/atmo/mo_async_latbc_types.f90        : Declaration of data types.
!! src/io/atmo/mo_latbc_read_recv.f90          : Low level routines: 
!!                                               read-in and sending of field data via MPI
!! 
!! src/atm_dyn_iconam/mo_initicon_types.f90    : Type declarations for the final destination buffers
!! src/atm_dyn_iconam/mo_nh_nest_utilities.f90 : Actual usage of the boundary data: buffers -> tendencies
!! 
!! src/namelists/mo_limarea_nml.f90            : Namelist "limarea_nml"
!! src/configure_model/mo_limarea_config.f90   : Configuration state (filled by namelist)
!! 
!! 
!! 
!! Important notes for understanding the read-in process:
!! ------------------------------------------------------
!! 
!! * General switch: "latbc_config%itype_latbc" (INTEGER)
!! 
!!   This is set to the value LATBC_TYPE_EXT when the field data is
!!   read from an external file (default situation).
!! 
!! * General switch: "latbc_config%lsparse_latbc" (LOGICAL)
!! 
!!   Lateral boundary data can be provided either as a boundary strip
!!   or as the complete field, including the (unused) interior points.
!! 
!! * Variables which are considered for boundary data read-in are
!!   contained in the variable group LATBC_PREFETCH_VARS.
!!   Buffers are set up and allocated for group members.
!! 
!! * Most important data types: "t_latbc_data", and contained therein: "t_buffer"
!! 
!!   Defined in src/io/atmo/mo_async_latbc_types.f90 These derived
!!   data types contain the raw data, that has been read from file,
!!   together with a number of configurations settings
!!   (LOGICAL-Flags), e.g. if specific optional fields are requested.
!! 
!! * "Decision tree": 
!! 
!!   During the setup phase, the boundary data module opens the first
!!   boundary data file and analyzes its contents (subroutine
!!   "check_variables" in "mo_async_latbc", contains numerous calls of
!!   "test_cdi_varID").
!! 
!!   According to the available data, several LOGICAL flags
!!   "lread_XXX" are set in the intermediate buffer "latbc%buffer".
!!
!! * Distinction between "read-in" and "fetch"
!! 
!!   (In the module "mo_async_latbc_utils":)
!!   "read-in": A list of variables is read from file into a buffer.
!!              * executed on the read-in process.
!!              * reads all fields of the group LATBC_TYPE_EXT
!!              * the intermediate buffer may have only the size of
!!              * the boundary strip.
!!              * most important subroutine in this context:
!!                "read_latbc_data" and, therein, "CALL
!!                prefetch_cdi_3d".
!!
!!   "Fetch": Copies the data from the intermediate buffer into
!!            ICON-internal variables.
!!              * executed on the compute processes.
!!              * copy/processing of the fields depends on the
!!                "lread_XXX" flags (see above), followed by vertical
!!                interpolation.
!!              * most important subroutine in this context:
!!                compute_latbc_intp_data (mo_async_utils) and,
!!                therein, "CALL fetch_from_buffer".
!!              * Note that the "ICON-internal variables" are not the
!!                prognostic fields, but again intermediate buffers of
!!                the data type "t_init_state" (Modul
!!                "mo_initicon_types").
!!                Example: "latbc%latbc_data(tlev)%atm_in%u", which is
!!                then later processed into tendencies (in
!!                "mo_nh_nest_utilities").
!! 
!! * Recipe: How to implement additional boundary data
!! 
!!   1) Modify the "add_var" calls of the new fields in the ICON code,
!!      such that these variables become members of the group
!!      "LATBC_TYPE_EXT".
!!   2) Extend the data type "t_buffer" by a LOGICAL flag for the new
!!      field "lread_XXX".
!!   3) Subroutine "check_variables" in "mo_async_latbc": Additional
!!      test, if the new field is available in the boundary data file;
!!      set the flag "lread_XXX" accordingly.
!!   4) Extend the "t_latbc_data" data structure by a buffer for the
!!      new field, similar to the contents of
!!      "latbc%latbc_data(tlev)%atm_in".
!!   5) Subroutine "compute_latbc_icon_data": Additional call to
!!      "fetch_from_buffer", contained in an IF-condition "IF
!!      (lread_XXX) ...".
!!   6) Implement the vertical interpolation for the new field.
!!   7) Use the new data in the application code.
!! 
!! 

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_async_latbc

  USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_ptr, c_intptr_t, c_f_pointer

#ifndef NOMPI
    USE mpi
#endif

    ! basic modules
    USE mo_kind,                      ONLY: i8, sp
    USE mo_exception,                 ONLY: finish, message, message_text
    USE mo_mpi,                       ONLY: stop_mpi, my_process_is_io,  my_process_is_pref, &
         &                                  my_process_is_mpi_test, p_int, p_real_sp,        &
         &                                  process_work_pref0, p_gather
    USE mo_parallel_config,           ONLY: nproma, num_prefetch_proc
    USE mo_model_domain,              ONLY: p_patch, t_patch
    ! MPI Communicators
    USE mo_mpi,                       ONLY: p_comm_work, p_comm_work_pref, p_comm_work_2_pref
    ! MPI Communication routines
    USE mo_mpi,                       ONLY: p_bcast, p_barrier
    ! MPI Process type intrinsics
    USE mo_mpi,                       ONLY: my_process_is_work
    ! MPI Process group sizes
    USE mo_mpi,                       ONLY: num_work_procs, p_n_work
    ! Processor numbers
    USE mo_mpi,                       ONLY: p_pe_work, p_work_pe0, p_comm_work_pref_compute_pe0
    USE mo_time_config,               ONLY: time_config
    USE mo_async_latbc_types,         ONLY: t_patch_data, t_reorder_data, t_latbc_data
    USE mo_grid_config,               ONLY: nroot
    USE mo_async_latbc_utils,         ONLY: read_latbc_data, compute_init_latbc_data, async_init_latbc_data,&
         &                                  compute_wait_for_async_pref, compute_shutdown_async_pref, &
         &                                  async_pref_send_handshake,  async_pref_wait_for_start, &
         &                                  allocate_pref_latbc_data
    USE mo_impl_constants,            ONLY: SUCCESS, MAX_CHAR_LENGTH, TIMELEVEL_SUFFIX
    USE mo_cdi_constants,             ONLY: GRID_UNSTRUCTURED_CELL, GRID_UNSTRUCTURED_EDGE
    USE mo_communication,             ONLY: idx_no, blk_no
    USE mo_nonhydro_state,            ONLY: p_nh_state
    USE mo_intp_data_strc,            ONLY: p_int_state
    USE mo_ext_data_state,            ONLY: ext_data
    USE mo_linked_list,               ONLY: t_var_list, t_list_element
    USE mo_var_metadata_types,        ONLY: t_var_metadata, VARNAME_LEN
    USE mo_var_list,                  ONLY: nvar_lists, var_lists, new_var_list, &
         &                                  collect_group
    USE mo_limarea_config,            ONLY: latbc_config, generate_filename, LATBC_TYPE_EXT
    USE mo_dictionary,                ONLY: t_dictionary, dict_get, dict_init, dict_loadfile, &
         &                                  dict_finalize
    USE mo_util_string,               ONLY: add_to_list, tolower
    USE mo_time_config,               ONLY: time_config
    USE mtime,                        ONLY: datetime, OPERATOR(+)
    USE mo_cdi,                       ONLY: vlistInqVarZaxis , streamOpenRead, streamInqVlist, &
         &                                  vlistNvars, zaxisInqSize, vlistInqVarName,         &
         &                                  streamClose, streamInqFiletype,                    &
         &                                  FILETYPE_NC2, FILETYPE_NC4, FILETYPE_GRB2
    USE mo_read_interface,            ONLY: nf
    USE mo_io_units,                  ONLY: filename_max
    USE mo_io_util,                   ONLY: read_netcdf_int_1d, t_netcdf_att_int
    USE mo_util_file,                 ONLY: util_filesize
    USE mo_util_cdi,                  ONLY: test_cdi_varID, cdiGetStringError
    USE mo_latbc_read_recv,           ONLY: prefetch_proc_send, compute_data_receive
    USE mo_sync,                      ONLY: sync_patch_array, SYNC_E, SYNC_C
#ifdef YAC_coupling
    USE mo_coupling_config,           ONLY: is_coupled_run
    USE mo_io_coupling,               ONLY: construct_io_coupler, destruct_io_coupler
#endif

    IMPLICIT NONE

    INCLUDE 'netcdf.inc'

    PRIVATE

    ! subroutines
    PUBLIC :: prefetch_main_proc
    PUBLIC :: init_prefetch
    PUBLIC :: close_prefetch


    !------------------------------------------------------------------------------------------------
    ! CONSTANTS
    !------------------------------------------------------------------------------------------------

    ! common constant strings
    CHARACTER(LEN=*), PARAMETER :: modname = 'mo_async_latbc'

    ! variables in this group are ICON, COSMO or IFS data which are
    ! read by the prefetch PE:
    CHARACTER(LEN=*), PARAMETER :: LATBC_PREFETCH_VARS = 'LATBC_PREFETCH_VARS'
    
    ! max. number of LATBC group variables
    INTEGER,          PARAMETER :: MAX_NUM_GRPVARS = 200

  CONTAINS

    !------------------------------------------------------------------------------------------------
    !> Close all name_list files and deallocate variables
    !
    SUBROUTINE close_prefetch()
#ifndef NOMPI

      IF (my_process_is_work()) THEN

         CALL compute_wait_for_async_pref()
         CALL compute_shutdown_async_pref()
      END IF
#endif
    END SUBROUTINE close_prefetch

    !------------------------------------------------------------------------------------------------
    ! The following routines are only needed for asynchronous Input prefetching
    !-------------------------------------------------------------------------------------------------
    !> Main routine for Input Prefetcing PEs.
    !  Please note that this routine never returns.
    !
    SUBROUTINE prefetch_main_proc()
      ! local variables
      CHARACTER(*), PARAMETER :: routine = modname//"::prefetch_main_proc"
      LOGICAL                 :: done
      TYPE(t_latbc_data)      :: latbc
      TYPE(datetime)          :: latbc_read_datetime

#ifdef YAC_coupling
      ! The initialisation of YAC needs to be called by all (!) MPI processes
      ! in MPI_COMM_WORLD.
      ! construct_io_coupler needs to be called before init_name_list_output
      ! due to calling sequence in subroutine atmo_model for other atmosphere
      ! processes
      IF ( is_coupled_run() ) CALL construct_io_coupler ( "prefetch_input_io" )
#endif

      ! call to initalize the prefetch processor with grid data
      CALL init_prefetch(latbc)
      ! Enter prefetch loop
      DO
         ! Wait for a message from the compute PEs to start
         CALL async_pref_wait_for_start(done)
         IF(done) EXIT ! leave loop, we are done
         ! perform input prefetching
         latbc_read_datetime = latbc%mtime_last_read + latbc%delta_dtime
         CALL read_latbc_data(latbc, latbc_read_datetime)
         ! Inform compute PEs that we are done
         CALL async_pref_send_handshake()
      END DO
      ! Finalization sequence:
      CALL close_prefetch()
      ! clean up
      CALL latbc%finalize()
#ifdef YAC_coupling
      IF ( is_coupled_run() ) CALL destruct_io_coupler ( "prefetch_input_io" )
#endif
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

#ifndef NOMPI
      IF (mpi_error /= MPI_SUCCESS) THEN
        WRITE (message_text, *) TRIM(mpi_call), ' returned with error=', mpi_error
        IF (l_finish) THEN
          CALL finish(routine, message_text)
        ELSE
          CALL message(routine, message_text)
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
    SUBROUTINE set_patch_data(latbc, bcast_root, latbc_varnames_dict)
      TYPE (t_latbc_data), INTENT(INOUT) :: latbc
      INTEGER,             INTENT(IN)    :: bcast_root
      TYPE (t_dictionary), INTENT(IN)    :: latbc_varnames_dict

#ifndef NOMPI
      ! local variables:
      CHARACTER(LEN=*), PARAMETER   :: routine = modname//"::set_patch_data"

      ! allocate patch data structure
      ! set number of global cells/edges/verts and patch ID

      IF (my_process_is_work()) THEN
         latbc%patch_data%nlev          = p_patch(1)%nlev
         latbc%patch_data%nlevp1        = p_patch(1)%nlevp1
         latbc%patch_data%level         = p_patch(1)%level
         latbc%patch_data%nblks_c       = p_patch(1)%nblks_c
         latbc%patch_data%nblks_e       = p_patch(1)%nblks_e
         latbc%patch_data%n_patch_cells = p_patch(1)%n_patch_cells
         latbc%patch_data%n_patch_edges = p_patch(1)%n_patch_edges
         latbc%patch_data%n_patch_cells_g = p_patch(1)%n_patch_cells_g
         latbc%patch_data%n_patch_edges_g = p_patch(1)%n_patch_edges_g
      END IF

      ! transfer data to prefetching PE
      CALL p_bcast(latbc%patch_data%nlev, bcast_root, p_comm_work_2_pref)
      CALL p_bcast(latbc%patch_data%nlevp1, bcast_root, p_comm_work_2_pref)
      CALL p_bcast(latbc%patch_data%level, bcast_root, p_comm_work_2_pref)
      CALL p_bcast(latbc%patch_data%nblks_c, bcast_root, p_comm_work_2_pref)
      CALL p_bcast(latbc%patch_data%nblks_e, bcast_root, p_comm_work_2_pref)
      CALL p_bcast(latbc%patch_data%n_patch_cells, bcast_root, p_comm_work_2_pref)
      CALL p_bcast(latbc%patch_data%n_patch_edges, bcast_root, p_comm_work_2_pref)
      CALL p_bcast(latbc%patch_data%n_patch_cells_g, bcast_root, p_comm_work_2_pref)
      CALL p_bcast(latbc%patch_data%n_patch_edges_g, bcast_root, p_comm_work_2_pref)

      ! ---------------------------------------------------------------------------
      ! replicate data on prefetch proc
      ! ---------------------------------------------------------------------------
      IF (.NOT. my_process_is_mpi_test()) CALL replicate_data_on_pref_proc(latbc%patch_data, bcast_root)

      IF(.NOT. my_process_is_pref()) THEN

         ! set reorder data on work PE
         CALL set_reorder_data( p_patch(1)%n_patch_cells_g, p_patch(1)%n_patch_cells, &
              p_patch(1)%cells%decomp_info%owner_mask, p_patch(1)%cells%decomp_info%glb_index, &
              latbc%patch_data%cells )

         CALL set_reorder_data( p_patch(1)%n_patch_edges_g, p_patch(1)%n_patch_edges, &
              p_patch(1)%edges%decomp_info%owner_mask, p_patch(1)%edges%decomp_info%glb_index, &
              latbc%patch_data%edges )

      ENDIF

      IF(.NOT. my_process_is_mpi_test()) THEN
         ! transfer reorder data to prefetch PE
         CALL transfer_reorder_data(bcast_root, latbc%patch_data%cells)
         CALL transfer_reorder_data(bcast_root, latbc%patch_data%edges)
      ENDIF

      ! subroutine to read const (height level) data and to check
      ! whether some variable is specified in input file and setting
      ! flag for its further usage
      IF (latbc_config%itype_latbc == LATBC_TYPE_EXT)  CALL check_variables(latbc, latbc_varnames_dict)

#endif

    END SUBROUTINE set_patch_data


    !> Based on a 1D list of global indices, available on the
    !> prefetching PE: Create a local LOGICAL mask on each worker PE
    !> which is TRUE, when this entry is present in the index list.
    !
    SUBROUTINE create_latbc_mask(latbc, p_ri, glb_indices, nindices_g, hgrid_type)

      TYPE(t_latbc_data),   INTENT(INOUT), TARGET :: latbc
      TYPE(t_reorder_data), INTENT(INOUT) :: p_ri           ! reorder info data structure
      INTEGER,              INTENT(IN)    :: glb_indices(:)
      INTEGER,              INTENT(IN)    :: nindices_g     ! global no. of indices
      INTEGER,              INTENT(IN)    :: hgrid_type     ! grid type (CELL/EDGE)
#ifndef NOMPI
      ! local variables:
      CHARACTER(LEN=*), PARAMETER :: routine = modname//"::create_latbc_mask"
      INTEGER                   :: nblks, sync_type, ierrstat
      INTEGER(i8)               :: eoff
      INTEGER(i8), ALLOCATABLE  :: ioff(:)
      REAL(sp),   ALLOCATABLE   :: var_out(:,:,:), tmp_buf(:)

      IF (my_process_is_pref()) THEN
        ALLOCATE(p_ri%pe_skip(0:num_work_procs-1+num_prefetch_proc), STAT=ierrstat)
        IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")
        p_ri%pe_skip(:) = .FALSE.
      ELSE
        SELECT CASE(hgrid_type)
        CASE (GRID_UNSTRUCTURED_CELL)
          nblks     = p_patch(1)%nblks_c
          sync_type = SYNC_C
        CASE (GRID_UNSTRUCTURED_EDGE)
          nblks     = p_patch(1)%nblks_e
          sync_type = SYNC_E
        CASE DEFAULT
          CALL finish(routine, "Internal error!")
        END SELECT
      END IF

      ! we use the send/receive routines here to determine which
      ! entries are read from file.!
      ! 
      ! Set a "1.0" on all LATBC points and send them to the compute
      ! PEs ...
      IF (my_process_is_pref()) THEN
        
        ALLOCATE(ioff(0:num_work_procs+num_prefetch_proc-1), tmp_buf(nindices_g), STAT=ierrstat)
        IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")
        tmp_buf(:)              = 0.0_sp
        tmp_buf(glb_indices(:)) = 1.0_sp
        ioff                    = 0
        CALL prefetch_proc_send(latbc%patch_data, tmp_buf(:), 1, hgrid_type, ioff)
        DEALLOCATE(ioff, tmp_buf, STAT=ierrstat)
        IF (ierrstat /= SUCCESS) CALL finish(routine, "DEALLOCATE failed!")
        
      END IF ! my_process_is_pref()

      CALL p_barrier(p_comm_work_pref)
      
      IF (.NOT. my_process_is_pref()) THEN
        ALLOCATE(var_out(nproma, 1, nblks), STAT=ierrstat)
        IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")
        var_out        = 0.0_sp
        eoff           = 0_i8
        p_ri%this_skip = .FALSE.
        CALL compute_data_receive(hgrid_type, 1, var_out, eoff, latbc%patch_data)

        ! build a logical mask, which local points are read from file
        !
        ! TODO: does compute data receive fill the halo points? then
        ! the following sync is obsolete:
        CALL sync_patch_array(sync_type, p_patch(1), var_out)
        p_ri%read_mask(:,:) = (var_out(:,1,:)  > 0._sp)

        DEALLOCATE(var_out, STAT=ierrstat)
        IF (ierrstat /= SUCCESS) CALL finish(routine, "DEALLOCATE failed!")
      END IF

      IF (.NOT. my_process_is_pref()) THEN
        p_ri%this_skip = .NOT. ANY(p_ri%read_mask(:,:))
      ELSE
        p_ri%this_skip = .TRUE.
      END IF

      CALL p_gather(p_ri%this_skip, p_ri%pe_skip, process_work_pref0, p_comm_work_pref)
#endif
    END SUBROUTINE create_latbc_mask


    !------------------------------------------------------------------------------------------------
    !> Replicates data (mainly the variable lists) needed for async prefetching
    !  on the prefetching procs.
    SUBROUTINE init_prefetch(latbc)
      TYPE(t_latbc_data), INTENT(INOUT), TARGET :: latbc

      ! local variables:
      CHARACTER(LEN=*), PARAMETER :: routine = modname//"::init_prefetch"

      ! grp name in lower case letter
      CHARACTER(LEN=VARNAME_LEN), ALLOCATABLE :: StrLowCasegrp(:)

      ! Broadcast root for intercommunicator broadcasts form compute
      ! PEs to prefetching PE using p_comm_work_2_pref
      INTEGER :: bcast_root

      TYPE (t_dictionary)           :: latbc_varnames_dict
      TYPE(t_netcdf_att_int)        :: opt_att(2)            ! optional attribute values
      INTEGER                       :: ierrstat, ic, idx_c, blk_c
      LOGICAL                       :: is_pref

      ! bcast_root is not used in this case
      bcast_root = 0

#ifndef NOMPI

      is_pref = my_process_is_pref()
      ! Set broadcast root for intercommunicator broadcasts
      IF (is_pref) THEN
         ! Root is proc 0 on the prefetch PE
         bcast_root = 0
      ELSE
        ! Special root setting for inter-communicators:
        ! The PE really sending must use MPI_ROOT, the others MPI_PROC_NULL
        bcast_root = MERGE(MPI_ROOT, MPI_PROC_NULL, p_pe_work == 0)
      ENDIF

      ! read the map file into dictionary data structure
      CALL dict_init(latbc_varnames_dict, lcase_sensitive=.FALSE.)

      IF(LEN_TRIM(latbc_config%latbc_varnames_map_file) > 0) THEN
         CALL dict_loadfile(latbc_varnames_dict, TRIM(latbc_config%latbc_varnames_map_file))
      END IF

      ! create and transfer patch data
      CALL set_patch_data(latbc, bcast_root, latbc_varnames_dict)

      ! open and read file containing information of prefetch variables
      ALLOCATE(StrLowCasegrp(MAX_NUM_GRPVARS))
      IF( my_process_is_work() ) THEN
        CALL read_init_file(latbc, StrLowCasegrp, latbc_varnames_dict, p_patch(1))
      ELSE IF ( my_process_is_pref() ) THEN
        CALL read_init_file(latbc, StrLowCasegrp, latbc_varnames_dict)
      ENDIF

      ! destroy variable name dictionaries:
      CALL dict_finalize(latbc_varnames_dict)

      ! initialize the memory window for communication
      IF (.NOT. my_process_is_mpi_test()) &
           CALL init_remote_memory_window(latbc, StrLowCasegrp)

      DEALLOCATE(StrLowCasegrp)

      ! --- "sparse latbc mode": read only data for boundary rows
      !
      !     this requires index information obtained from an additional
      !     grid file:
      ALLOCATE(latbc%patch_data%cells%read_mask(nproma, latbc%patch_data%nblks_c), &
        &      latbc%patch_data%edges%read_mask(nproma, latbc%patch_data%nblks_e), STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")
      latbc%patch_data%cells%read_mask(:,:) = .TRUE.
      latbc%patch_data%edges%read_mask(:,:) = .TRUE.

      IF (latbc_config%lsparse_latbc) THEN

        CALL message(routine, "sparse LATBC read-in mode.")

        IF (is_pref) THEN
       
          opt_att(1)%var_name = "global_cell_index"
          opt_att(1)%att_name = "nglobal"
          opt_att(1)%value => latbc%global_index%n_patch_cells_g
          
          opt_att(2)%var_name = "global_edge_index"
          opt_att(2)%att_name = "nglobal"
          opt_att(2)%value => latbc%global_index%n_patch_edges_g
                   
          CALL read_netcdf_int_1d(latbc_config%latbc_boundary_grid,                   &
            &                     varname1     = "global_cell_index",                 &
            &                     var1         = latbc%global_index%cells,            &
            &                     opt_varname2 = "global_edge_index",                 &
            &                     opt_var2     = latbc%global_index%edges,            &
            &                     opt_att      = opt_att)

          ! consistency checks:
          IF (latbc%global_index%n_patch_cells_g /= latbc%patch_data%n_patch_cells_g) THEN
            CALL finish(routine, "LatBC boundary cell list does not match in size!")
          END IF
          IF (latbc%global_index%n_patch_edges_g /= latbc%patch_data%n_patch_edges_g) THEN
            CALL finish(routine, "LatBC boundary edge list does not match in size!")
          END IF
        END IF

        CALL create_latbc_mask(latbc, latbc%patch_data%cells, latbc%global_index%cells, &
          &                    latbc%global_index%n_patch_cells_g, GRID_UNSTRUCTURED_CELL)
        CALL create_latbc_mask(latbc, latbc%patch_data%edges, latbc%global_index%edges, &
          &                    latbc%global_index%n_patch_edges_g, GRID_UNSTRUCTURED_EDGE)

        ! status output
        IF (my_process_is_pref()) THEN
          WRITE (0,'(3a,i0,a,i0,a)') &
            &   " ", routine, ": prefetching PE reads ", latbc%global_index%n_patch_cells_g, &
            &   " cells and ", latbc%global_index%n_patch_edges_g, " edges."
          WRITE (0,'(3a,i0,a,i0,a)')      &
            &   " ", routine, ": ", COUNT(.NOT. latbc%patch_data%edges%pe_skip), " of ", num_work_procs, &
            &   " worker PEs are involved in the LATBC read-in."
        END IF

        IF (.NOT. my_process_is_pref()) THEN
          ! consistency check: test if all nudging points are filled by
          ! the LATBC read-in
          DO ic=1,p_nh_state(1)%metrics%nudge_c_dim
            idx_c = p_nh_state(1)%metrics%nudge_c_idx(ic)
            blk_c = p_nh_state(1)%metrics%nudge_c_blk(ic)
            IF (.NOT. latbc%patch_data%cells%read_mask(idx_c,blk_c)) THEN
              CALL finish(routine, "Nudging zone width mismatch: "//&
                &" Not all nudging points are filled by the LATBC READ-in.")
            END IF
          END DO
        END IF

      ELSE

        CALL message(routine, "non-sparse LATBC read-in mode.")

      END IF ! lsparse_latbc

      ! allocate input data for lateral boundary nudging
      IF( my_process_is_work()) THEN
        CALL compute_init_latbc_data(latbc, p_patch(1), p_int_state(1), p_nh_state(1), &
          &                          latbc%new_latbc_tlev)
      ELSE IF (is_pref) THEN
        CALL async_init_latbc_data(latbc)
      ENDIF

      CALL message(routine,'Done')
#endif

    END SUBROUTINE init_prefetch


    !-------------------------------------------------------------------------------------------------
    !> open files containing first variable list and analysis
    !
    SUBROUTINE read_init_file(latbc, StrLowCasegrp, latbc_varnames_dict, p_patch)
      TYPE (t_latbc_data),        INTENT(INOUT) :: latbc
      CHARACTER(LEN=VARNAME_LEN), INTENT(INOUT) :: StrLowCasegrp(:) !< grp name in lower case letter
      TYPE (t_dictionary),        INTENT(IN)    :: latbc_varnames_dict
      TYPE(t_patch), OPTIONAL,    INTENT(IN)    :: p_patch

      CHARACTER(*), PARAMETER                   :: routine = modname//"::read_init_files"
      LOGICAL,      PARAMETER                   :: ldebug  = .FALSE.
#ifndef NOMPI
      ! local variables
      CHARACTER(LEN=VARNAME_LEN), ALLOCATABLE :: grp_vars(:)
      ! dictionary which maps prefetch variable names onto
      ! GRIB2 shortnames or NetCDF var names.
      INTEGER                                 :: ierrstat, vlistID, nvars, varID, zaxisID,  &
        &                                        jp, fileID_latbc, counter, filetype, ngrp_prefetch_vars, &
        &                                        nlev_in, ncid
      INTEGER(KIND=i8)                        :: flen_latbc
      LOGICAL                                 :: l_exist
      CHARACTER(LEN=filename_max)             :: latbc_filename, latbc_file
      CHARACTER(LEN=MAX_CHAR_LENGTH)          :: name, cdiErrorText

      ! allocating buffers containing name of variables
      ALLOCATE(grp_vars(MAX_NUM_GRPVARS))

      !>Looks for variable groups ("group:xyz") and collects
      ! them to map prefetch variable names onto
      ! GRIB2 shortnames or NetCDF var names.

      ! loop over all variables and collects the variables names
      ! corresponding to the group "LATBC_PREFETCH_VARS"

      CALL collect_group(LATBC_PREFETCH_VARS, grp_vars, ngrp_prefetch_vars, loutputvars_only=.FALSE., &
        &                lremap_lonlat=.FALSE. )

      ! adding the variable 'GEOSP' to the list by add_to_list
      ! as the variable cannot be found in metadata variable list
      IF (latbc_config%itype_latbc == LATBC_TYPE_EXT) &
           CALL add_to_list( grp_vars, ngrp_prefetch_vars, (/latbc%buffer%geop_ml_var/) , 1)

      ! allocate the number of vertical levels and other fields with
      ! the same size as number of variables
      ALLOCATE(latbc%buffer%nlev(ngrp_prefetch_vars),        &
        &      latbc%buffer%grp_vars(ngrp_prefetch_vars),    &
        &      latbc%buffer%mapped_name(ngrp_prefetch_vars), &
        &      latbc%buffer%internal_name(ngrp_prefetch_vars), STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")

      ! allocate the array for variable ID
      ! with same size as number of variables
      ALLOCATE(latbc%buffer%varID(ngrp_prefetch_vars), STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")

      ! generate file name
      latbc_filename = generate_filename(nroot, latbc%patch_data%level, &
        &                                time_config%tc_exp_startdate, time_config%tc_exp_startdate)
      latbc_file = TRIM(latbc_config%latbc_path)//TRIM(latbc_filename)

      IF(my_process_is_work() .AND.  p_pe_work == p_work_pe0) THEN

         INQUIRE (FILE=latbc_file, EXIST=l_exist)
         IF (.NOT.l_exist) THEN
            CALL finish(routine,'LATBC file not found: '//TRIM(latbc_file))
         ENDIF

         ! open file
         !
         fileID_latbc = streamOpenRead(TRIM(latbc_file))
         IF (fileID_latbc < 0) THEN
           CALL cdiGetStringError(fileID_latbc, cdiErrorText)
           CALL finish(routine, "File "//TRIM(latbc_file)//" cannot be opened: "//TRIM(cdiErrorText))
         ENDIF

         filetype = streamInqFiletype(fileID_latbc)
         IF (.NOT. ANY(filetype == [FILETYPE_NC2, FILETYPE_NC4,FILETYPE_GRB2])) THEN
           CALL finish(routine, "Unknown file type")
         END IF

         IF (latbc_config%itype_latbc == LATBC_TYPE_EXT) THEN
           ! Search name mapping for name in file
           DO jp= 1, ngrp_prefetch_vars
             latbc%buffer%grp_vars(jp) = TRIM(dict_get(latbc_varnames_dict, grp_vars(jp), default=grp_vars(jp)))
           ENDDO
         ELSE
           DO jp= 1, ngrp_prefetch_vars
             latbc%buffer%grp_vars(jp) = TRIM(grp_vars(jp))
           ENDDO
         ENDIF

         ! check whether the file is empty (does not work unfortunately; internal CDI error)
         flen_latbc = util_filesize(TRIM(latbc_file))
         IF (flen_latbc <= 0 ) THEN
            CALL message(routine, "File "//TRIM(latbc_file)//" is empty")
            CALL finish(routine, "STOP: Empty input file")
         ENDIF

         vlistID = streamInqVlist(fileID_latbc)

         ! get the number of variables
         nvars = vlistNvars(vlistID)

         ! initialising counter
         counter = 0

         ! get the number of vertical levels for the
         ! required prefetch variables
         LOOP : DO varID=0,(nvars-1)
            CALL vlistInqVarName(vlistID, varID, name)

            IF (ldebug)  WRITE(0,*) 'name ', TRIM(name)

            DO jp = 1, ngrp_prefetch_vars !latbc%buffer%ngrp_vars
               IF(tolower(name) == tolower(latbc%buffer%grp_vars(jp))) THEN
                  ! get the vertical axis ID
                  zaxisID = vlistInqVarZaxis(vlistID, varID)

                  counter = counter + 1
                  ! get the respective vertical levels for
                  ! the respective variable
                  latbc%buffer%nlev(counter) = zaxisInqSize(zaxisID)
                  !variable ID for variable to be read from file
                  latbc%buffer%varID(counter) = varID

                  ! getting the variable name in the file
                  latbc%buffer%mapped_name(counter) = TRIM(name)
                  latbc%buffer%internal_name(counter) = &
                    TRIM(dict_get(latbc_varnames_dict, TRIM(name), linverse=.TRUE., default=TRIM(name)))
                  ! getting the variable name in lower case letter
                  StrLowCasegrp(counter) = grp_vars(jp)

                  IF (ldebug) THEN
                    WRITE(0,*) '=> internal_name ',  (latbc%buffer%internal_name(counter))
                    WRITE(0,*) '=> mapped_name   ',  (latbc%buffer%mapped_name(counter))
                  END IF
               ENDIF
            ENDDO
         END DO LOOP

         ! closes the open dataset
         CALL streamClose(fileID_latbc)
      
       END IF

      CALL p_bcast(latbc%buffer%nlev(:),          p_comm_work_pref_compute_pe0, p_comm_work_pref)
      CALL p_bcast(latbc%buffer%varID(:),         p_comm_work_pref_compute_pe0, p_comm_work_pref)
      CALL p_bcast(latbc%buffer%mapped_name(:),   p_comm_work_pref_compute_pe0, p_comm_work_pref)
      CALL p_bcast(latbc%buffer%internal_name(:), p_comm_work_pref_compute_pe0, p_comm_work_pref)
      CALL p_bcast(StrLowCasegrp(:),              p_comm_work_pref_compute_pe0, p_comm_work_pref)
      CALL p_bcast(counter,                       p_comm_work_pref_compute_pe0, p_comm_work_pref)

      !WRITE(0,*) 'mapped_name ',  latbc%buffer%mapped_name(1:counter), 'ngrp_vars ', ngrp_prefetch_vars

      ! getting the count of number of variables in
      ! the file to be read
      latbc%buffer%ngrp_vars = counter

      ! allocate the number of buffer sizes for variables
      ! with same size as number of variables
      ALLOCATE(latbc%buffer%vars(latbc%buffer%ngrp_vars), STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")

      ! allocate latbc buffer, use the maximum no. of vertical levels
      ! for any variable:
      IF (my_process_is_work()) THEN
        nlev_in = 0
        IF (p_pe_work == p_work_pe0) THEN
          ! set the maximum no. of levels to the size of the half
          ! level height field (HHL/z_ifc) minus 1.
          IF (latbc%buffer%lcompute_hhl_pres) THEN
            nlev_in = MAXVAL(latbc%buffer%nlev(1:latbc%buffer%ngrp_vars)) ! IFS input data have no vertical staggering
          ELSE
            nlev_in = MAXVAL(latbc%buffer%nlev(1:latbc%buffer%ngrp_vars))-1 ! All other supported input sources have staggering
          ENDIF
        END IF
        CALL p_bcast(nlev_in, 0, p_comm_work)
        CALL allocate_pref_latbc_data(latbc, nlev_in, p_nh_state(1), ext_data(1), p_patch)
      END IF

      ! clean up
      DEALLOCATE(grp_vars, STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish(routine, "DEALLOCATE failed!")


      ! Re-open the file to read constant fields.
      !
      IF (my_process_is_work() .AND. .NOT. my_process_is_pref()) THEN

        IF (latbc%buffer%lcompute_hhl_pres) THEN
          CALL nf(nf_open(TRIM(latbc_file), NF_NOWRITE, ncid), routine)
          CALL latbc%latbc_data_const%vct%construct(ncid, p_work_pe0, p_comm_work)
          CALL nf(nf_close(ncid), routine)
        END IF
      END IF

#endif
    END SUBROUTINE read_init_file


    !-------------------------------------------------------------------------------------------------
    !> Open the first of the LATBC files (start date) to determine if
    !  a variable or its alternative variable is provided as input and
    !  setting flag for its further usage.
    !
    SUBROUTINE check_variables(latbc, latbc_dict)
      TYPE(t_latbc_data),  INTENT(INOUT) :: latbc
      TYPE (t_dictionary), INTENT(IN)    :: latbc_dict

#ifndef NOMPI
      ! local variables
      CHARACTER(*), PARAMETER        :: routine = modname//"::check_variables"
      INTEGER                        :: fileID_latbc
      LOGICAL                        :: l_exist, lhave_ps_geop, lhave_ps, lhave_geop,  &
        &                               lhave_hhl, lhave_theta_rho, lhave_vn,          &
        &                               lhave_u, lhave_v, lhave_pres, lhave_temp
      CHARACTER(LEN=filename_max)    :: latbc_filename, latbc_file
      CHARACTER(LEN=MAX_CHAR_LENGTH) :: cdiErrorText

      ! prefetch processor opens the file and checks if variables are present
      IF( my_process_is_work() .AND.  p_pe_work == p_work_pe0) THEN !!!!!!!use prefetch processor here
         ! generate file name
         latbc_filename = generate_filename(nroot, latbc%patch_data%level, &
           &                                time_config%tc_exp_startdate, time_config%tc_exp_startdate)
         latbc_file = TRIM(latbc_config%latbc_path)//TRIM(latbc_filename)
         INQUIRE (FILE=latbc_file, EXIST=l_exist)
         IF (.NOT.l_exist) THEN
            CALL finish(routine,'LATBC file not found: '//TRIM(latbc_file))
         ENDIF

         ! open file
         !
         fileID_latbc = streamOpenRead(TRIM(latbc_file))
         IF (fileID_latbc < 0) THEN
           CALL cdiGetStringError(fileID_latbc, cdiErrorText)
           CALL finish(routine, "File "//TRIM(latbc_file)//" cannot be opened: "//TRIM(cdiErrorText))
         ENDIF

         ! Check if rain water (QR) is provided as input
         latbc%buffer%lread_qr = (test_cdi_varID(fileID_latbc, 'QR', latbc_dict) /= -1)

         ! Check if snow water (QS) is provided as input
         latbc%buffer%lread_qs = (test_cdi_varID(fileID_latbc, 'QS', latbc_dict) /= -1)


         ! --- CHECK WHICH VARIABLES ARE AVAILABLE IN THE DATA SET ---

         ! Check if vertical velocity (or OMEGA) is provided as input
         latbc%buffer%lread_w = (test_cdi_varID(fileID_latbc, 'W', latbc_dict) /= -1)

         ! Check if surface pressure (VN) is provided as input
         lhave_vn = (test_cdi_varID(fileID_latbc, 'VN', latbc_dict) /= -1)
         lhave_u  = (test_cdi_varID(fileID_latbc, 'U', latbc_dict)  /= -1)
         lhave_v  = (test_cdi_varID(fileID_latbc, 'V', latbc_dict)  /= -1)

         ! Check if the prognostic thermodynamic variables (rho and
         ! theta_v) are provided as input
         lhave_theta_rho = (test_cdi_varID(fileID_latbc, 'RHO', latbc_dict) /= -1)  .OR.  &
           &               (test_cdi_varID(fileID_latbc, 'DEN', latbc_dict) /= -1)  .AND. &
           &               (test_cdi_varID(fileID_latbc, 'THETA_V', latbc_dict) /= -1)


         ! Check if level heights are provided as input
         lhave_hhl =  .FALSE.
         IF (test_cdi_varID(fileID_latbc, 'Z_IFC', latbc_dict) /= -1) THEN
           lhave_hhl = .TRUE.
           latbc%buffer%hhl_var   = 'Z_IFC'
         END IF

         !
         ! Check if surface pressure (PS) or its logarithm (LNPS) is provided as input
         !
         lhave_ps = .FALSE.
         IF (test_cdi_varID(fileID_latbc, 'PS', latbc_dict) /= -1) THEN
            lhave_ps = .TRUE.
            latbc%buffer%psvar    = 'PS'
         ELSE IF (test_cdi_varID(fileID_latbc, 'LNPS', latbc_dict) /= -1) THEN
            lhave_ps = .TRUE.
            latbc%buffer%psvar    = 'LNPS'
         ENDIF

         !
         ! Check if model-level surface Geopotential is provided as GEOSP or GEOP_ML
         !
         lhave_geop = .FALSE.
         IF (test_cdi_varID(fileID_latbc, 'GEOSP', latbc_dict) /= -1) THEN
            lhave_geop  = .TRUE.
            latbc%buffer%geop_ml_var = 'GEOSP'
         ELSE IF (test_cdi_varID(fileID_latbc, 'GEOP_ML', latbc_dict) /= -1) THEN
            lhave_geop  = .TRUE.
            latbc%buffer%geop_ml_var = 'GEOP_ML'
         ELSE IF (.NOT. lhave_theta_rho .AND. .NOT. lhave_hhl) THEN
            CALL finish(routine,'Could not find model-level sfc geopotential')
         ENDIF
         lhave_ps_geop = (lhave_ps .and. lhave_geop)

         ! Check if pressure and temperature are available:
         lhave_pres = (test_cdi_varID(fileID_latbc, 'PRES', latbc_dict) /= -1)
         lhave_temp = (test_cdi_varID(fileID_latbc, 'TEMP', latbc_dict) /= -1)


         ! closes the open dataset
         CALL streamClose(fileID_latbc)


         ! --- DEFINE WHICH VARIABLES SHALL BE READ FROM FILE ---

         latbc%buffer%lread_hhl         = .FALSE.
         latbc%buffer%lread_theta_rho   = .FALSE.
         latbc%buffer%lread_pres        = .FALSE.
         latbc%buffer%lread_temp        = .FALSE.
         latbc%buffer%lread_ps_geop     = .FALSE.
         latbc%buffer%lconvert_omega2w  = .FALSE.
         latbc%buffer%lcompute_hhl_pres = .FALSE.

         IF (lhave_hhl) THEN

           IF (lhave_theta_rho) THEN
             latbc%buffer%lread_hhl       = .TRUE.
             latbc%buffer%lread_theta_rho = .TRUE.
             !
           ELSE IF (lhave_pres .AND. lhave_temp) THEN
             latbc%buffer%lread_hhl       = .TRUE.
             latbc%buffer%lread_pres      = .TRUE.
             latbc%buffer%lread_temp      = .TRUE.
             !
           ELSE
             CALL finish(routine, "Non-hydrostatic LATBC data set, but neither RHO+THETA_V nor P,T provided!")
           END IF

         ELSE
           
           IF (lhave_ps_geop .AND. lhave_temp) THEN
             latbc%buffer%lread_temp        = .TRUE.
             latbc%buffer%lread_ps_geop     = .TRUE.
             latbc%buffer%lconvert_omega2w  = .TRUE.
             latbc%buffer%lcompute_hhl_pres = .TRUE.
           ELSE
             CALL finish(routine, "Hydrostatic LATBC data set, but PS,GEOP,T not provided!")
           END IF

         END IF

         latbc%buffer%lread_vn  = .FALSE.
         latbc%buffer%lread_u_v = .FALSE.
         IF (lhave_vn) THEN
           latbc%buffer%lread_vn = .TRUE.
         ELSE
           IF (lhave_u .AND. lhave_v) THEN
             latbc%buffer%lread_u_v = .TRUE.
           ELSE
             CALL finish(routine, "No VN or U&V available in LATBC data set!")
           END IF
         END IF


         !
         ! Consistency checks
         ! 
         IF (latbc_config%init_latbc_from_fg .AND. .NOT. latbc%buffer%lread_hhl) THEN
           CALL finish(routine, "Init LATBC from first guess requires BCs from non-hydrostatic model!")
         END IF


         !
         ! Write some status output:
         !
         IF (latbc%buffer%lread_theta_rho) THEN
           CALL message(routine,'Prognostic thermodynamic variables (RHO and THETA_V) are read from file.')
         ENDIF

         IF (.NOT. latbc%buffer%lread_qr) THEN
           CALL message(routine,'Rain water (QR) not available in input data.')
         ENDIF

         IF (.NOT. latbc%buffer%lread_qs) THEN
            CALL message(routine,'Snow water (QS) not available in input data.')
         ENDIF

         IF (latbc%buffer%lread_hhl) THEN
            CALL message(routine,'Input levels (HHL) are read from file.')
         ELSE
            CALL message(routine,'Input levels (HHL) are computed from sfc geopotential.')
         ENDIF

         IF (.NOT. latbc%buffer%lread_w) THEN
           CALL message(routine, "Neither W nor OMEGA provided! W is set to zero at LBCs")
         ELSE IF (latbc%buffer%lconvert_omega2w) THEN
            CALL message(routine,'Compute W from OMEGA.')
         ENDIF

         IF (latbc%buffer%lread_ps_geop) THEN
           CALL message(routine,'PS and GEOP are read from file.')
         ELSE IF (lhave_ps_geop) THEN
           CALL message(routine,'PS and GEOP are ignored.')
         END IF

         IF (latbc%buffer%lcompute_hhl_pres) THEN
           CALL message(routine,'HHL and PRES are computed based on PS and GEOP.')
         END IF

         IF (latbc%buffer%lread_pres) THEN
           CALL message(routine,'PRES is read from file.')
         ELSE
           CALL message(routine,'PRES is diagnosed.')
         END IF

         IF (latbc%buffer%lread_temp) THEN
           CALL message(routine,'TEMP is read from file.')
         ELSE
           CALL message(routine,'TEMP is diagnosed.')
         END IF

         IF (latbc%buffer%lread_vn) THEN
           CALL message(routine,'VN is read from file.')
         ELSE
           CALL message(routine,'U,V are read from file.')
         END IF

       ENDIF

      ! broadcast data to prefetching and compute PE's
      ! public constant: p_comm_work_pref_compute_pe0
      CALL p_bcast(latbc%buffer%psvar,                    p_comm_work_pref_compute_pe0, p_comm_work_pref)
      CALL p_bcast(latbc%buffer%geop_ml_var,              p_comm_work_pref_compute_pe0, p_comm_work_pref)
      CALL p_bcast(latbc%buffer%hhl_var,                  p_comm_work_pref_compute_pe0, p_comm_work_pref)
      CALL p_bcast(latbc%buffer%lread_qs,                 p_comm_work_pref_compute_pe0, p_comm_work_pref)
      CALL p_bcast(latbc%buffer%lread_qr,                 p_comm_work_pref_compute_pe0, p_comm_work_pref)
      CALL p_bcast(latbc%buffer%lread_vn,                 p_comm_work_pref_compute_pe0, p_comm_work_pref)
      CALL p_bcast(latbc%buffer%lread_u_v,                p_comm_work_pref_compute_pe0, p_comm_work_pref)
      CALL p_bcast(latbc%buffer%lread_w,                  p_comm_work_pref_compute_pe0, p_comm_work_pref)

      CALL p_bcast(latbc%buffer%lread_hhl,                p_comm_work_pref_compute_pe0, p_comm_work_pref)
      CALL p_bcast(latbc%buffer%lread_theta_rho,          p_comm_work_pref_compute_pe0, p_comm_work_pref)
      CALL p_bcast(latbc%buffer%lread_ps_geop,            p_comm_work_pref_compute_pe0, p_comm_work_pref)
      CALL p_bcast(latbc%buffer%lread_pres,               p_comm_work_pref_compute_pe0, p_comm_work_pref)
      CALL p_bcast(latbc%buffer%lread_temp,               p_comm_work_pref_compute_pe0, p_comm_work_pref)
      CALL p_bcast(latbc%buffer%lconvert_omega2w,         p_comm_work_pref_compute_pe0, p_comm_work_pref)
      CALL p_bcast(latbc%buffer%lcompute_hhl_pres,        p_comm_work_pref_compute_pe0, p_comm_work_pref)
#endif

    END SUBROUTINE check_variables


    !-------------------------------------------------------------------------------------------------
    !> Replicates data needed for async prefetch on the prefetch proc.
    !  ATTENTION: The data is not completely replicated, only as far as needed for prefetch.
    !
    !  This routine has to be called by all PEs (work and prefetch)
    !
    SUBROUTINE replicate_data_on_pref_proc(patch_data, bcast_root)
      TYPE(t_patch_data),  INTENT(INOUT) :: patch_data
      INTEGER,             INTENT(IN) :: bcast_root

#ifndef NOMPI
      ! local variables
      CHARACTER(len=*), PARAMETER   :: routine = modname//"::replicate_data_on_pref_proc"
      INTEGER                       :: info_size, i, iv, nelems, nv, n, list_info, all_var, &
           &                           all_nvars, nvars, i2, ierrstat
      LOGICAL                       :: is_pref
      INTEGER, ALLOCATABLE          :: info_storage(:,:)
      TYPE(t_list_element), POINTER :: element
      TYPE(t_var_metadata)          :: info
      TYPE(t_var_list)              :: p_var_list
      ! var_list_name should have at least the length of var_list names
      CHARACTER(LEN=256)            :: var_list_name

      ! get the size - in default INTEGER words - which is needed to
      ! hold the contents of TYPE(t_var_metadata)
      info_size = SIZE(TRANSFER(info, (/ 0 /)))

      is_pref = my_process_is_pref()
      ! get the number of var lists
      IF (.NOT. is_pref) nv = nvar_lists
      CALL p_bcast(nv, bcast_root, p_comm_work_2_pref)

      IF (.NOT. is_pref) THEN
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
      IF (.NOT. is_pref) all_var = all_nvars
      CALL p_bcast(all_var, bcast_root, p_comm_work_2_pref)

      IF (all_var <= 0) RETURN

      ! allocate the array of variables
      ALLOCATE(patch_data%var_data(all_var))

      i2 = 0
      ! For each var list, get its components
      DO iv = 1, nv

         ! Send name
         IF (.NOT. is_pref) var_list_name = var_lists(iv)%p%name
         CALL p_bcast(var_list_name, bcast_root, p_comm_work_2_pref)

         IF (.NOT. is_pref) THEN
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

         IF (is_pref) THEN
            nelems = list_info
            ! Create var list
            CALL new_var_list( p_var_list, var_list_name)
         ENDIF

         ! Get the binary representation of all info members of the
         ! variables of the list and send it to the receiver.  Using
         ! the Fortran TRANSFER intrinsic may seem like a hack, but it
         ! has the advantage that it is completely independent from
         ! the actual declaration if TYPE(t_var_metadata).  Thus
         ! members may added to or removed from TYPE(t_var_metadata)
         ! without affecting the code below and we don't have an
         ! additional cross dependency between TYPE(t_var_metadata)
         ! and this module.

         ALLOCATE(info_storage(info_size, nelems), STAT=ierrstat)
         IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")

         IF (.NOT. is_pref) THEN
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

         IF (is_pref) THEN

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
         DEALLOCATE(info_storage, STAT=ierrstat)
         IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")
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
      IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")
      ALLOCATE(p_reo%own_blk(p_reo%n_own), STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")

      ! Global index of my own points
      ALLOCATE(glbidx_own(p_reo%n_own), STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")

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
      IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")
      ALLOCATE(p_reo%pe_off(0:p_n_work-1), STAT=ierrstat)
      IF (ierrstat  /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")

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
      IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")

      ! make sure all are done
      ! CALL p_barrier(comm=p_comm_work)

      CALL MPI_Allgatherv(glbidx_own, p_reo%n_own, p_int, &
           glbidx_glb, p_reo%pe_own, p_reo%pe_off, p_int, &
           p_comm_work, mpi_error)

      CALL check_mpi_error(routine, 'MPI_Allgatherv', mpi_error, .TRUE.)

      ! Get reorder_index
      ALLOCATE(p_reo%reorder_index(p_reo%n_glb), STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")

      ! Spans the complete logical domain
      ALLOCATE(reorder_index_log_dom(n_points_g), STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")
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
         WRITE (message_text, *) 'Reordering failed: n=',n, ' /= p_reo%n_glb=',p_reo%n_glb
         CALL finish(routine, message_text)
      ENDIF

      DEALLOCATE(phys_owner_mask, glbidx_own, glbidx_glb, reorder_index_log_dom, STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish(routine, "DEALLOCATE failed!")
#endif

    END SUBROUTINE set_reorder_data


    !------------------------------------------------------------------------------------------------
    !
    ! Transfers reorder data to restart PEs.
    !
    SUBROUTINE transfer_reorder_data(bcast_root, p_reo)
      INTEGER,              INTENT(IN)    :: bcast_root
      TYPE(t_reorder_data), INTENT(INOUT) :: p_reo

#ifndef NOMPI
      ! local variables
      INTEGER                             :: ierrstat
      CHARACTER(LEN=*), PARAMETER :: routine = modname//"::transfer_reorder_data"

      !transfer the global number of points, this is not known on prefetching PE
      CALL p_bcast(p_reo%n_glb, bcast_root, p_comm_work_2_pref)

      IF(my_process_is_pref())THEN

         ! on prefetch PE: n_own=0, own_ide and own_blk are not allocated
         p_reo%n_own = 0

         ! pe_own must be allocated for num_work_procs, not for p_n_work
         ALLOCATE(p_reo%pe_own(0:num_work_procs-1), STAT=ierrstat)
         IF(ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")

         ALLOCATE(p_reo%pe_off(0:num_work_procs-1), STAT=ierrstat)
         IF(ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")

         ALLOCATE(p_reo%reorder_index(p_reo%n_glb), STAT=ierrstat)
         IF(ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")
      ENDIF

      CALL p_bcast(p_reo%pe_own,        bcast_root, p_comm_work_2_pref)
      CALL p_bcast(p_reo%pe_off,        bcast_root, p_comm_work_2_pref)
      CALL p_bcast(p_reo%reorder_index, bcast_root, p_comm_work_2_pref)
#endif
    END SUBROUTINE  transfer_reorder_data


    !------------------------------------------------------------------------------------------------
    !> Initializes the remote memory window for asynchronous prefetch.
    !
    SUBROUTINE init_remote_memory_window(latbc, StrLowCasegrp)
      TYPE (t_latbc_data),        INTENT(INOUT) :: latbc
      CHARACTER(LEN=VARNAME_LEN), INTENT(IN)    :: StrLowCasegrp(:) !< grp name in lower case letter

#ifndef NOMPI
      INTEGER :: ierrstat, iv, jp, nlevs
      INTEGER (KIND=MPI_ADDRESS_KIND) :: mem_size
      LOGICAL ,ALLOCATABLE :: grp_vars_bool(:)
      CHARACTER(LEN=*), PARAMETER :: routine = modname//"::init_memory_window"

      latbc%patch_data%mem_win%mpi_win = MPI_WIN_NULL

      ! Get size and offset of the data for the input
      mem_size = 0_i8

      ALLOCATE(grp_vars_bool(latbc%buffer%ngrp_vars))
      ALLOCATE(latbc%buffer%hgrid(latbc%buffer%ngrp_vars))

      grp_vars_bool(1:latbc%buffer%ngrp_vars) = .FALSE.
      
      ! Go over all input variables
      DO iv = 1, SIZE(latbc%patch_data%var_data)
         DO jp = 1, latbc%buffer%ngrp_vars
            ! Use only the variables of time level 1 (".TL1") to determine memory sizes.
            IF((TRIM(StrLowCasegrp(jp)) == TRIM(latbc%patch_data%var_data(iv)%info%name)) .OR. &
                 & (TRIM(StrLowCasegrp(jp))//TIMELEVEL_SUFFIX//'1' == TRIM(latbc%patch_data%var_data(iv)%info%name))) THEN

               nlevs = 0
               IF(.NOT. grp_vars_bool(jp))  THEN
                  IF (latbc%patch_data%var_data(iv)%info%ndims == 2) THEN
                     nlevs = 1
                  ELSE
                     nlevs = latbc%buffer%nlev(jp) 
                  ENDIF

                  IF (nlevs == 0) CYCLE

                  SELECT CASE (latbc%patch_data%var_data(iv)%info%hgrid)

                  CASE (GRID_UNSTRUCTURED_CELL)
                     mem_size = mem_size + INT(nlevs*latbc%patch_data%cells%n_own,i8)

                     IF(my_process_is_work())THEN
                        ! allocate the buffer sizes for variables on compute processors
                        ALLOCATE(latbc%buffer%vars(jp)%buffer(nproma, nlevs, &
                          &      latbc%patch_data%nblks_c), STAT=ierrstat)
                        IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")
                     ENDIF

                     ! variable stored in cell center location
                     latbc%buffer%hgrid(jp) = latbc%patch_data%var_data(iv)%info%hgrid

                  CASE (GRID_UNSTRUCTURED_EDGE)
                     mem_size = mem_size + INT(nlevs*latbc%patch_data%edges%n_own,i8)

                     IF(my_process_is_work())THEN
                        ! allocate the buffer sizes for variables on compute processors
                        ALLOCATE(latbc%buffer%vars(jp)%buffer(nproma, nlevs, &
                          &      latbc%patch_data%nblks_e), STAT=ierrstat)
                        IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")
                     ENDIF

                     ! variable stored in edge center of a cell
                     latbc%buffer%hgrid(jp) = latbc%patch_data%var_data(iv)%info%hgrid

                  CASE DEFAULT
                     CALL finish(routine,'Unknown grid type!')
                  END SELECT

                  grp_vars_bool(jp) = .TRUE.

               ENDIF
            ENDIF
         ENDDO
      ENDDO ! vars

      IF (latbc_config%itype_latbc == LATBC_TYPE_EXT) THEN
         DO jp = 1, latbc%buffer%ngrp_vars
            IF (TRIM(latbc%buffer%mapped_name(jp)) == TRIM(latbc%buffer%geop_ml_var)) THEN
               ! Memory for GEOSP variable taken as memory equivalent to 1 level of z_ifc
               ! as the variable GEOSP doesn't exist in metadata
               mem_size = mem_size + INT(1*latbc%patch_data%cells%n_own,i8)

               IF(my_process_is_work())THEN
                  ! allocate the buffer sizes for variable 'GEOSP' on compute processors
                  ALLOCATE(latbc%buffer%vars(jp)%buffer(nproma, 1, latbc%patch_data%nblks_c), STAT=ierrstat)
                  IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")
               ENDIF

               ! variable GEOSP is stored in cell center location
               latbc%buffer%hgrid(jp) = GRID_UNSTRUCTURED_CELL
            ENDIF
         ENDDO
      ENDIF

      ! allocate amount of memory needed with MPI_Alloc_mem
      CALL allocate_mem_noncray(latbc%patch_data, MAX(mem_size,1_i8))
#endif

    END SUBROUTINE init_remote_memory_window


    !------------------------------------------------------------------------------------------------
    !> allocate amount of memory needed with MPI_Alloc_mem
    !
    !  @note Implementation for non-Cray pointers
    !
    SUBROUTINE allocate_mem_noncray(patch_data, mem_size)

      TYPE(t_patch_data),              INTENT(INOUT) :: patch_data
#ifdef NOMPI
      INTEGER,                         INTENT(IN)    :: mem_size
#else
      INTEGER (KIND=MPI_ADDRESS_KIND), INTENT(IN)    :: mem_size

      ! local variables
      CHARACTER(LEN=*), PARAMETER :: routine = modname//"::allocate_mem_noncray"
      TYPE(c_ptr)                     :: c_mem_ptr
      INTEGER                         :: mpierr, nbytes_real
      INTEGER (KIND=MPI_ADDRESS_KIND) :: mem_bytes

      ! Get the amount of bytes per REAL*4 variable (as used in MPI
      ! communication)
      CALL MPI_Type_extent(p_real_sp, nbytes_real, mpierr)

      ! For the IO PEs the amount of memory needed is 0 - allocate at least 1 word there:
      mem_bytes = mem_size*INT(nbytes_real,i8)

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

      NULLIFY(patch_data%mem_win%mem_ptr_sp)
      CALL C_F_POINTER(c_mem_ptr, patch_data%mem_win%mem_ptr_sp, (/ mem_size /) )

      ! Create memory window for communication
      patch_data%mem_win%mem_ptr_sp(:) = 0._sp
      CALL MPI_Win_create( patch_data%mem_win%mem_ptr_sp, mem_bytes, nbytes_real, MPI_INFO_NULL,&
        &                  p_comm_work_pref, patch_data%mem_win%mpi_win, mpierr )
      IF (mpierr /= 0) CALL finish(routine, "MPI error!")
#endif

    END SUBROUTINE allocate_mem_noncray

END MODULE mo_async_latbc
