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
!!                "prefetch_latbc_data" and, therein, "CALL
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

  USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_ptr, c_f_pointer

#ifndef NOMPI
    USE mpi
#endif

    ! basic modules
    USE mo_kind,                      ONLY: i8, sp
    USE mo_exception,                 ONLY: finish, message, message_text
    USE mo_mpi,                       ONLY: stop_mpi, my_process_is_pref, &
         &                                  my_process_is_mpi_test, p_real_sp, &
         &                                  p_reduce, mpi_sum, p_allgather, &
         &                                  p_allgatherv, p_bcast
    USE mo_parallel_config,           ONLY: nproma, num_prefetch_proc
    USE mo_decomposition_tools,       ONLY: t_grid_domain_decomp_info
    USE mo_model_domain,              ONLY: p_patch, t_patch
    ! MPI Communicators
    USE mo_mpi,                       ONLY: p_comm_work, p_comm_work_pref, p_comm_work_2_pref
    ! MPI Communication routines
    ! MPI Process type intrinsics
    USE mo_mpi,                       ONLY: my_process_is_work
    ! MPI Process group sizes
    USE mo_mpi,                       ONLY: num_work_procs
    ! Processor numbers
    USE mo_mpi,                       ONLY: p_pe_work, p_work_pe0, p_comm_work_pref_compute_pe0
    USE mo_time_config,               ONLY: time_config
    USE mo_reorder_info,              ONLY: t_reorder_info
    USE mo_async_latbc_types,         ONLY: t_patch_data, t_buffer, &
                                            t_latbc_data, t_mem_win
    USE mo_grid_config,               ONLY: nroot
    USE mo_async_latbc_utils,         ONLY: prefetch_latbc_data, read_init_latbc_data, async_init_latbc_data,&
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
    USE mo_limarea_config,            ONLY: latbc_config, generate_filename
    USE mo_dictionary,                ONLY: t_dictionary, dict_get, dict_init, dict_loadfile, &
         &                                  dict_finalize
    USE mo_util_string,               ONLY: add_to_list, tolower
    USE mo_util_sort,                 ONLY: quicksort
    USE mo_time_config,               ONLY: time_config
    USE mtime,                        ONLY: datetime, OPERATOR(+)
    USE mo_cdi,                       ONLY: vlistInqVarZaxis , streamOpenRead, streamInqVlist, &
         &                                  vlistNvars, zaxisInqSize, vlistInqVarName,         &
         &                                  streamClose, streamInqFiletype,                    &
         &                                  FILETYPE_NC2, FILETYPE_NC4, FILETYPE_GRB2,         &
         &                                  cdi_max_name
    USE mo_read_interface,            ONLY: nf
    USE mo_io_util,                   ONLY: read_netcdf_int_1d, t_netcdf_att_int
    USE mo_util_file,                 ONLY: util_filesize
    USE mo_util_cdi,                  ONLY: test_cdi_varID, cdiGetStringError
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

  TYPE t_var_data
    TYPE(t_var_metadata), POINTER :: info  ! Info structure for variable
  END TYPE t_var_data

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
         CALL prefetch_latbc_data(latbc, latbc_read_datetime)
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

    !--------------------------------------------------------------------------
    !
    !> Initialize data structures for prefetching boundary data
    !
    !  This routine is called after reading the namelists AND setting up
    !  the domains and variables.
    !
#ifndef NOMPI
    SUBROUTINE set_patch_data_params(patch_data, bcast_root)
      TYPE(t_patch_data), INTENT(INOUT) :: patch_data
      INTEGER, INTENT(IN) :: bcast_root

      ! local variables:
      CHARACTER(LEN=*), PARAMETER   :: routine = modname//"::set_patch_data_paraams"

      ! allocate patch data structure
      ! set number of global cells/edges/verts and patch ID

      IF (my_process_is_work()) THEN
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

    END SUBROUTINE set_patch_data_params
#endif

    ! ------------------------------------------------------------------------
    ! replicate data on prefetch proc
    ! ------------------------------------------------------------------------
#ifndef NOMPI
    SUBROUTINE set_patch_data(patch_data, cell_ro_idx, edge_ro_idx)
      TYPE(t_patch_data), INTENT(INOUT) :: patch_data
      INTEGER, ALLOCATABLE, INTENT(in) :: cell_ro_idx(:), edge_ro_idx(:)

      IF(.NOT. my_process_is_pref()) THEN

        IF (ALLOCATED(cell_ro_idx)) THEN
          CALL set_reorder_data(p_patch(1)%n_patch_cells, &
            patch_data%cell_mask, cell_ro_idx, &
            patch_data%cells)
        ELSE
          CALL set_reorder_data(p_patch(1)%n_patch_cells, &
            patch_data%cell_mask, p_patch(1)%cells%decomp_info%glb_index, &
            patch_data%cells)
        END IF

        IF (ALLOCATED(edge_ro_idx)) THEN
          CALL set_reorder_data(p_patch(1)%n_patch_edges, &
            patch_data%edge_mask, edge_ro_idx, &
            patch_data%edges)
        ELSE
          CALL set_reorder_data(p_patch(1)%n_patch_edges, &
            patch_data%edge_mask, p_patch(1)%edges%decomp_info%glb_index, &
            patch_data%edges)
        END IF

      ENDIF

      IF(.NOT. my_process_is_mpi_test()) THEN
         ! transfer reorder data to prefetch PE
         CALL transfer_reorder_data(patch_data%cells)
         CALL transfer_reorder_data(patch_data%edges)
      ENDIF

    END SUBROUTINE set_patch_data


    !> Based on a 1D list of global indices, available on the
    !> prefetching PE: Create a local LOGICAL mask on each worker PE
    !> which is TRUE, when this entry is present in the index list.
    !
    SUBROUTINE create_latbc_mask_work(n, mask, ro_idx, decomp_info)
      INTEGER, INTENT(in) :: n
      LOGICAL, INTENT(out) :: mask(n)
      INTEGER, ALLOCATABLE, INTENT(out) :: ro_idx(:)
      TYPE(t_grid_domain_decomp_info), INTENT(in) :: decomp_info
      ! local variables:
      CHARACTER(LEN=*), PARAMETER :: routine = modname//"::create_latbc_mask"
      INTEGER                   :: nr, rcnt, i, m, ierror, nfound, dummy(1)
      INTEGER, ALLOCATABLE, TARGET :: ranges(:,:)
      INTEGER, POINTER :: p_ranges(:,:)

      ! get number of ranges describing bc indices
      CALL p_bcast(nr, 0, comm=p_comm_work_2_pref)
      IF (nr > 0) THEN
        ALLOCATE(ro_idx(n), ranges(nr,3), stat=ierror)
        IF (ierror /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")
        p_ranges => ranges(:,1:2)
        CALL p_bcast(p_ranges, 0, comm=p_comm_work_2_pref)
        m = SIZE(decomp_info%glb_index)
        IF (m > 0) THEN
          ! construct partial sums so it becomes easier to compute
          ! the actual position of an index in the read-in data
          rcnt = 0
          DO i = 1, nr
            ranges(i, 3) = rcnt
            rcnt = rcnt + ranges(i, 2) - ranges(i, 1) + 1
          END DO
          CALL mask_from_range_matches(mask, ro_idx, nr, ranges, &
               decomp_info%glb_index, nfound)
          mask(m+1:n) = .FALSE.
        END IF
        dummy = p_reduce(MERGE(1, 0, nfound>0), mpi_sum, 0, comm=p_comm_work_2_pref)
      ELSE
        nfound = 0
        mask = .FALSE.
      END IF
    END SUBROUTINE create_latbc_mask_work

    SUBROUTINE mask_from_range_matches(mask, ro_idx, nr, ranges, g_idx, nfound)
      LOGICAL, INTENT(out) :: mask(*)
      INTEGER, INTENT(in) :: nr
      INTEGER, INTENT(in) :: ranges(nr,3), g_idx(:)
      INTEGER, INTENT(out) :: nfound, ro_idx(*)

      INTEGER :: n, ir, rs, re, jl, ig

      n = SIZE(g_idx)
      nfound = 0
      g_idx_loop: DO jl = 1, n
        ig = g_idx(jl)
        DO ir = 1, nr
          rs = ranges(ir, 1)
          re = ranges(ir, 2)
          IF (ig < rs) EXIT
          IF (ig <= re) THEN
            mask(jl) = .TRUE.
            nfound = nfound + 1
            ! in case of sparse latbc, global indices need to be remapped
            ro_idx(jl) = ranges(ir, 3) + ig - rs + 1
            CYCLE g_idx_loop
          END IF
        END DO
        mask(jl) = .FALSE.
        ro_idx(jl) = -1
      END DO g_idx_loop
    END SUBROUTINE mask_from_range_matches

    !> Based on a 1D list of global indices, available on the
    !> prefetching PE: Create a local LOGICAL mask on each worker PE
    !> which is TRUE, when this entry is present in the index list.
    !
    SUBROUTINE create_latbc_mask_pref(glb_indices, active_worker_count)
      INTEGER, ALLOCATABLE, INTENT(IN)    :: glb_indices(:)
      INTEGER, INTENT(out) :: active_worker_count

      CHARACTER(LEN=*), PARAMETER :: &
           routine = modname//"::create_latbc_mask_pref"
      INTEGER                   :: ierror, i, n, root_pref2work, &
           rs, re, jl, nr, dummy
      INTEGER, ALLOCATABLE      :: sorted_glb_index(:), ranges(:,:)

      n = SIZE(glb_indices)
      IF (n > 0) THEN
        ALLOCATE(sorted_glb_index(n), stat=ierror)
        IF (ierror /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")
        sorted_glb_index = glb_indices
        CALL quicksort(sorted_glb_index)
        ! count ranges
        re = sorted_glb_index(1)
        nr = 1
        DO i = 2, n
          jl = sorted_glb_index(i)
          nr = nr + MERGE(1, 0, jl /= re .AND. jl /= re + 1)
          re = jl
        END DO
        ALLOCATE(ranges(nr, 2), stat=ierror)
        IF (ierror /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")
        ! convert sorted indices into ranges
        rs = sorted_glb_index(1)
        re = sorted_glb_index(1)
        nr = 1
        DO i = 2, n
          jl = sorted_glb_index(i)
          IF (jl == re .OR. jl == re + 1) THEN
            re = jl
          ELSE
            ranges(nr, 1) = rs
            ranges(nr, 2) = re
            rs = jl
            re = jl
            nr = nr + 1
          END IF
        END DO
        ranges(nr, 1) = rs
        ranges(nr, 2) = re
      ELSE
        nr = 0
      END IF
      root_pref2work = MERGE(mpi_root, mpi_proc_null, p_pe_work == 0)
      ! bcast number of ranges describing glb_indices to clients
      IF (my_process_is_pref()) THEN
        CALL p_bcast(nr, root_pref2work, comm=p_comm_work_2_pref)
        ! bcast ranges to clients
        IF (nr > 0) THEN
         CALL p_bcast(ranges, root_pref2work, comm=p_comm_work_2_pref)
          active_worker_count = p_reduce(dummy, mpi_sum, root_pref2work, &
               comm=p_comm_work_2_pref)
        ELSE
          active_worker_count = 0
        END IF
      ENDIF
    END SUBROUTINE create_latbc_mask_pref
#endif


    !------------------------------------------------------------------------------------------------
    !> Replicates data (mainly the variable lists) needed for async prefetching
    !  on the prefetching procs.
    SUBROUTINE init_prefetch(latbc)
      TYPE(t_latbc_data), INTENT(INOUT), TARGET :: latbc

#ifndef NOMPI
      ! local variables:
      CHARACTER(LEN=*), PARAMETER :: routine = modname//"::init_prefetch"

      ! grp name in lower case letter
      CHARACTER(LEN=VARNAME_LEN), ALLOCATABLE :: StrLowCasegrp(:)

      ! Broadcast root for intercommunicator broadcasts form compute
      ! PEs to prefetching PE using p_comm_work_2_pref
      INTEGER :: bcast_root

      TYPE (t_dictionary)           :: latbc_varnames_dict
      TYPE(t_netcdf_att_int)        :: opt_att(2)            ! optional attribute values
      INTEGER                       :: ierrstat, ic, idx_c, blk_c, cell_active_ranks, edge_active_ranks
      INTEGER                       :: tlen, covered
      INTEGER                       :: fileID_latbc
      INTEGER(mpi_address_kind)     :: mem_size
      LOGICAL                       :: is_pref, is_work, is_test
      INTEGER, ALLOCATABLE :: cell_ro_idx(:), edge_ro_idx(:), var_buf_map(:)
      TYPE(t_var_data), ALLOCATABLE :: var_data(:)

      is_pref = my_process_is_pref()
      is_work = my_process_is_work()
      is_test = my_process_is_mpi_test()
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
      tlen = LEN_TRIM(latbc_config%latbc_varnames_map_file)
      IF(tlen > 0) &
         CALL dict_loadfile(latbc_varnames_dict, latbc_config%latbc_varnames_map_file(1:tlen))

      CALL set_patch_data_params(latbc%patch_data, bcast_root)

      ! --- "sparse latbc mode": read only data for boundary rows
      !
      !     this requires index information obtained from an additional
      !     grid file:
      IF (is_work) THEN
        ALLOCATE(latbc%patch_data%cell_mask(nproma, latbc%patch_data%nblks_c), &
          &      latbc%patch_data%edge_mask(nproma, latbc%patch_data%nblks_e), STAT=ierrstat)
        IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")
      END IF

      IF (latbc_config%lsparse_latbc) THEN

        CALL message(routine, "sparse LATBC read-in mode.")

        IF (is_pref .OR. is_work .AND.  p_pe_work == p_work_pe0) THEN
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

          CALL create_latbc_mask_pref(latbc%global_index%cells, &
            &                         cell_active_ranks)
          CALL create_latbc_mask_pref(latbc%global_index%edges, &
            &                         edge_active_ranks)
          ! status output
          IF (is_pref) THEN
            WRITE (0,'(3a,i0,a,i0,a)') &
              &   " ", routine, ": prefetching PE reads ", SIZE(latbc%global_index%cells), &
              &   " cells and ", SIZE(latbc%global_index%edges), " edges."
            WRITE (0,'(3a,2(i0,a))') &
              &   " ", routine, ": ", edge_active_ranks, " of ", num_work_procs, &
              &   " worker PEs are involved in the LATBC read-in."
          ENDIF
        ENDIF
        IF (is_work) THEN
          CALL create_latbc_mask_work(SIZE(latbc%patch_data%cell_mask), &
            &                         latbc%patch_data%cell_mask, &
            &                         cell_ro_idx, &
            &                         p_patch(1)%cells%decomp_info)
          CALL create_latbc_mask_work(SIZE(latbc%patch_data%edge_mask), &
            &                         latbc%patch_data%edge_mask, &
            &                         edge_ro_idx, &
            &                         p_patch(1)%edges%decomp_info)
          ! consistency check: test if all nudging points are filled by
          ! the LATBC read-in
          covered = 0
          DO ic=1,p_nh_state(1)%metrics%nudge_c_dim
            idx_c = p_nh_state(1)%metrics%nudge_c_idx(ic)
            blk_c = p_nh_state(1)%metrics%nudge_c_blk(ic)
            covered = covered &
              + MERGE(1, 0, latbc%patch_data%cell_mask(idx_c,blk_c))
          END DO
          IF (covered /= p_nh_state(1)%metrics%nudge_c_dim) THEN
            WRITE (message_text, '(2(a,i0),a)') &
              'Nudging zone width mismatch: only ', covered, ' of ', &
              p_nh_state(1)%metrics%nudge_c_dim, &
              ' nudging points are filled by the LATBC READ-in.'
            CALL finish(routine, message_text)
          END IF
        END IF
      ELSE
        IF (is_work) THEN
          latbc%patch_data%cell_mask(:,:) = .TRUE.
          latbc%patch_data%edge_mask(:,:) = .TRUE.
        END IF

        CALL message(routine, "non-sparse LATBC read-in mode.")

      END IF ! lsparse_latbc

      IF (.NOT. is_test) CALL replicate_data_on_pref_proc(var_data, bcast_root)

      ! create and transfer patch data
      CALL set_patch_data(latbc%patch_data, cell_ro_idx, edge_ro_idx)

      ! open and read file containing information of prefetch variables
      ALLOCATE(StrLowCasegrp(MAX_NUM_GRPVARS))
      IF (is_work) THEN
        CALL read_init_file(latbc, StrLowCasegrp, latbc_varnames_dict, p_patch(1), fileID_latbc)
      ELSE IF (is_pref) THEN
        CALL read_init_file(latbc, StrLowCasegrp, latbc_varnames_dict)
      ENDIF

      CALL match_var_data_to_buffer(var_buf_map, latbc%buffer, &
           var_data, StrLowCasegrp)
      DEALLOCATE(StrLowCasegrp)
      CALL infer_buffer_hgrid(latbc%buffer, var_data, var_buf_map)

      ! initialize the memory window for communication
      IF (.NOT. is_test) THEN
        IF (is_work) THEN
          CALL create_client_buffers(latbc%buffer, var_data,&
            latbc%patch_data%cells%n_own, latbc%patch_data%nblks_c, &
            latbc%patch_data%edges%n_own, latbc%patch_data%nblks_e, &
            var_buf_map, mem_size)
        ELSE IF (is_pref) THEN
          mem_size = 0_mpi_address_kind
        END IF
        ! allocate amount of memory needed with MPI_Alloc_mem
        CALL create_win(latbc%patch_data%mem_win, mem_size)
      END IF


      ! allocate input data for lateral boundary nudging
      IF (is_work) THEN
        CALL read_init_latbc_data(latbc, p_patch(1), p_int_state(1), p_nh_state(1), &
          &                       latbc%new_latbc_tlev, fileID_latbc, latbc_varnames_dict)
      ELSE IF (is_pref) THEN
        CALL async_init_latbc_data(latbc)
      ENDIF

      CALL message(routine,'Done')

      ! destroy variable name dictionaries:
      CALL dict_finalize(latbc_varnames_dict)
#endif

    END SUBROUTINE init_prefetch


    !-------------------------------------------------------------------------------------------------
    !> open files containing first variable list and analysis
    !
    SUBROUTINE read_init_file(latbc, StrLowCasegrp, latbc_varnames_dict, p_patch, fileID)
      TYPE (t_latbc_data),        INTENT(INOUT) :: latbc
      CHARACTER(LEN=VARNAME_LEN), INTENT(INOUT) :: StrLowCasegrp(:) !< grp name in lower case letter
      TYPE (t_dictionary),        INTENT(IN)    :: latbc_varnames_dict
      TYPE(t_patch), OPTIONAL,    INTENT(IN)    :: p_patch
      INTEGER,       OPTIONAL,    INTENT(OUT)   :: fileID

      CHARACTER(*), PARAMETER                   :: routine = modname//"::read_init_files"
      LOGICAL,      PARAMETER                   :: ldebug  = .FALSE.
#ifndef NOMPI
      ! local variables
      CHARACTER(LEN=VARNAME_LEN), ALLOCATABLE :: grp_vars(:), grp_vars_lc(:)
      ! dictionary which maps prefetch variable names onto
      ! GRIB2 shortnames or NetCDF var names.
      INTEGER                                 :: ierrstat, vlistID, nvars, varID, zaxisID,  &
        &                                        jp, fileID_latbc, counter, filetype, ngrp_prefetch_vars, &
        &                                        nlev_in, ncid, tlen
      INTEGER(KIND=i8)                        :: flen_latbc
      LOGICAL                                 :: l_exist, is_work
      CHARACTER(LEN=:), ALLOCATABLE           :: latbc_filename, latbc_file
      CHARACTER(LEN=CDI_MAX_NAME)             :: name, name_lc
      CHARACTER(LEN=MAX_CHAR_LENGTH)          :: cdiErrorText

      is_work = my_process_is_work()

      ! allocating buffers containing name of variables
      ALLOCATE(grp_vars(MAX_NUM_GRPVARS))

      !>Looks for variable groups ("group:xyz") and collects
      ! them to map prefetch variable names onto
      ! GRIB2 shortnames or NetCDF var names.

      ! loop over all variables and collects the variables names
      ! corresponding to the group "LATBC_PREFETCH_VARS"

      CALL collect_group(LATBC_PREFETCH_VARS, grp_vars, ngrp_prefetch_vars, loutputvars_only=.FALSE., &
        &                lremap_lonlat=.FALSE. )

      ! allocate the number of vertical levels and other fields with
      ! the same size as number of variables
      ALLOCATE(latbc%buffer%nlev(ngrp_prefetch_vars),        &
        &      latbc%buffer%mapped_name(ngrp_prefetch_vars), &
        &      latbc%buffer%internal_name(ngrp_prefetch_vars), STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")

      ! allocate the array for variable ID
      ! with same size as number of variables
      ALLOCATE(latbc%buffer%varID(ngrp_prefetch_vars), STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")

      IF (is_work .AND. p_pe_work == p_work_pe0) THEN
        ! generate file name
        latbc_filename = TRIM(generate_filename(nroot, latbc%patch_data%level, &
             &                                  time_config%tc_exp_startdate,  &
             &                                  time_config%tc_exp_startdate))
        latbc_file = TRIM(latbc_config%latbc_path)//latbc_filename
      END IF
      IF (is_work .AND.  p_pe_work == p_work_pe0) THEN

        INQUIRE (FILE=latbc_file, EXIST=l_exist)
         IF (.NOT.l_exist) THEN
            CALL finish(routine,'LATBC file not found: '//latbc_file)
         ENDIF

         ! open file
         !
         fileID_latbc = streamOpenRead(latbc_file)
         IF (fileID_latbc < 0) THEN
           CALL cdiGetStringError(fileID_latbc, cdiErrorText)
           CALL finish(routine, "File "//latbc_file//" cannot be opened: "//TRIM(cdiErrorText))
         ENDIF

         filetype = streamInqFiletype(fileID_latbc)
         IF (.NOT. ANY(filetype == [FILETYPE_NC2, FILETYPE_NC4,FILETYPE_GRB2])) THEN
           CALL finish(routine, "Unknown file type")
         END IF


        ! subroutine to read const (height level) data and to check
        ! whether some variable is specified in input file and setting
        ! flag for its further usage
        CALL check_variables(latbc%buffer, latbc%patch_data,latbc_varnames_dict, fileID_latbc)

        ! adding the variable 'GEOSP' to the list by add_to_list
        ! as the variable cannot be found in metadata variable list
        CALL add_to_list( grp_vars, ngrp_prefetch_vars, (/latbc%buffer%geop_ml_var/) , 1)

         ALLOCATE(grp_vars_lc(ngrp_prefetch_vars), stat=ierrstat)
         IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")

         ! Search name mapping for name in file
         DO jp= 1, ngrp_prefetch_vars
           grp_vars_lc(jp) = tolower(dict_get(latbc_varnames_dict, grp_vars(jp), default=grp_vars(jp)))
         ENDDO

         ! check whether the file is empty (does not work unfortunately; internal CDI error)
         flen_latbc = util_filesize(latbc_file)
         IF (flen_latbc <= 0 ) THEN
            CALL message(routine, "File "//latbc_file//" is empty")
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

            name_lc = tolower(name)
            DO jp = 1, ngrp_prefetch_vars !latbc%buffer%ngrp_vars
               IF (name_lc == grp_vars_lc(jp)) THEN
                  ! get the vertical axis ID
                  zaxisID = vlistInqVarZaxis(vlistID, varID)

                  counter = counter + 1
                  ! get the respective vertical levels for
                  ! the respective variable
                  latbc%buffer%nlev(counter) = zaxisInqSize(zaxisID)
                  !variable ID for variable to be read from file
                  latbc%buffer%varID(counter) = varID

                  ! getting the variable name in the file
                  latbc%buffer%mapped_name(counter) = name
                  tlen = LEN_TRIM(name)
                  latbc%buffer%internal_name(counter) = &
                    dict_get(latbc_varnames_dict, name(1:tlen), linverse=.TRUE., default=name(1:tlen))
                  ! getting the variable name in lower case letter
                  StrLowCasegrp(counter) = grp_vars(jp)

                  IF (ldebug) THEN
                    WRITE(0,*) '=> internal_name ',  (latbc%buffer%internal_name(counter))
                    WRITE(0,*) '=> mapped_name   ',  (latbc%buffer%mapped_name(counter))
                  END IF
               ENDIF
            ENDDO
         END DO LOOP

         DEALLOCATE(grp_vars_lc, STAT=ierrstat)
         IF (ierrstat /= SUCCESS) CALL finish(routine, "DEALLOCATE failed!")

         ! save fileID for further processing
         IF (PRESENT(fileID)) THEN
           fileID = fileID_latbc ! broadcast needed? apparently not.
         ENDIF

       END IF

      ! broadcast data to prefetching and compute PE's
      ! public constant: p_comm_work_pref_compute_pe0
      CALL p_bcast(latbc%buffer%nlev(:),           p_comm_work_pref_compute_pe0, p_comm_work_pref)
      CALL p_bcast(latbc%buffer%varID(:),          p_comm_work_pref_compute_pe0, p_comm_work_pref)
      CALL p_bcast(latbc%buffer%mapped_name(:),    p_comm_work_pref_compute_pe0, p_comm_work_pref)
      CALL p_bcast(latbc%buffer%internal_name(:),  p_comm_work_pref_compute_pe0, p_comm_work_pref)
      CALL p_bcast(StrLowCasegrp(:),               p_comm_work_pref_compute_pe0, p_comm_work_pref)
      CALL p_bcast(counter,                        p_comm_work_pref_compute_pe0, p_comm_work_pref)

      CALL p_bcast(latbc%buffer%psvar,             p_comm_work_pref_compute_pe0, p_comm_work_pref)
      CALL p_bcast(latbc%buffer%geop_ml_var,       p_comm_work_pref_compute_pe0, p_comm_work_pref)
      CALL p_bcast(latbc%buffer%hhl_var,           p_comm_work_pref_compute_pe0, p_comm_work_pref)
      CALL p_bcast(latbc%buffer%lread_qs,          p_comm_work_pref_compute_pe0, p_comm_work_pref)
      CALL p_bcast(latbc%buffer%lread_qr,          p_comm_work_pref_compute_pe0, p_comm_work_pref)
      CALL p_bcast(latbc%buffer%lread_vn,          p_comm_work_pref_compute_pe0, p_comm_work_pref)
      CALL p_bcast(latbc%buffer%lread_u_v,         p_comm_work_pref_compute_pe0, p_comm_work_pref)
      CALL p_bcast(latbc%buffer%lread_w,           p_comm_work_pref_compute_pe0, p_comm_work_pref)

      CALL p_bcast(latbc%buffer%lread_hhl,         p_comm_work_pref_compute_pe0, p_comm_work_pref)
      CALL p_bcast(latbc%buffer%lread_theta_rho,   p_comm_work_pref_compute_pe0, p_comm_work_pref)
      CALL p_bcast(latbc%buffer%lread_ps_geop,     p_comm_work_pref_compute_pe0, p_comm_work_pref)
      CALL p_bcast(latbc%buffer%lread_pres,        p_comm_work_pref_compute_pe0, p_comm_work_pref)
      CALL p_bcast(latbc%buffer%lread_temp,        p_comm_work_pref_compute_pe0, p_comm_work_pref)
      CALL p_bcast(latbc%buffer%lconvert_omega2w,  p_comm_work_pref_compute_pe0, p_comm_work_pref)
      CALL p_bcast(latbc%buffer%lcompute_hhl_pres, p_comm_work_pref_compute_pe0, p_comm_work_pref)

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
      IF (is_work) THEN
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
      IF (latbc%buffer%lcompute_hhl_pres .AND. is_work) THEN
        IF (p_pe_work == p_work_pe0) CALL nf(nf_open(latbc_file, NF_NOWRITE, ncid), routine)
        CALL latbc%latbc_data_const%vct%construct(ncid, p_work_pe0, p_comm_work)
        IF (p_pe_work == p_work_pe0) CALL nf(nf_close(ncid), routine)
      END IF

#endif
    END SUBROUTINE read_init_file


    !-------------------------------------------------------------------------------------------------
    !> Open the first of the LATBC files (start date) to determine if
    !  a variable or its alternative variable is provided as input and
    !  setting flag for its further usage.
    !
    SUBROUTINE check_variables(buffer, patch_data, latbc_dict, fileID_latbc)
      TYPE(t_buffer), INTENT(inout)  :: buffer
      TYPE(t_patch_data), INTENT(in) :: patch_data
      TYPE(t_dictionary), INTENT(in) :: latbc_dict
      INTEGER, INTENT(IN)            :: fileID_latbc

#ifndef NOMPI
      ! local variables
      CHARACTER(*), PARAMETER        :: routine = modname//"::check_variables"
      LOGICAL                        :: l_exist, lhave_ps_geop, lhave_ps, lhave_geop,  &
        &                               lhave_hhl, lhave_theta_rho, lhave_vn,          &
        &                               lhave_u, lhave_v, lhave_pres, lhave_temp


         ! Check if rain water (QR) is provided as input
         buffer%lread_qr = (test_cdi_varID(fileID_latbc, 'QR', latbc_dict) /= -1)

         ! Check if snow water (QS) is provided as input
         buffer%lread_qs = (test_cdi_varID(fileID_latbc, 'QS', latbc_dict) /= -1)


         ! --- CHECK WHICH VARIABLES ARE AVAILABLE IN THE DATA SET ---

         ! Check if vertical velocity (or OMEGA) is provided as input
         buffer%lread_w = (test_cdi_varID(fileID_latbc, 'W', latbc_dict) /= -1)

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
           buffer%hhl_var   = 'Z_IFC'
         END IF

         !
         ! Check if surface pressure (PS) or its logarithm (LNPS) is provided as input
         !
         lhave_ps = .FALSE.
         IF (test_cdi_varID(fileID_latbc, 'PS', latbc_dict) /= -1) THEN
            lhave_ps = .TRUE.
            buffer%psvar    = 'PS'
         ELSE IF (test_cdi_varID(fileID_latbc, 'LNPS', latbc_dict) /= -1) THEN
            lhave_ps = .TRUE.
            buffer%psvar    = 'LNPS'
         ENDIF

         !
         ! Check if model-level surface Geopotential is provided as GEOSP or GEOP_ML
         !
         lhave_geop = .FALSE.
         IF (test_cdi_varID(fileID_latbc, 'GEOSP', latbc_dict) /= -1) THEN
            lhave_geop  = .TRUE.
            buffer%geop_ml_var = 'GEOSP'
         ELSE IF (test_cdi_varID(fileID_latbc, 'GEOP_ML', latbc_dict) /= -1) THEN
            lhave_geop  = .TRUE.
            buffer%geop_ml_var = 'GEOP_ML'
         ELSE IF (.NOT. lhave_theta_rho .AND. .NOT. lhave_hhl) THEN
            CALL finish(routine,'Could not find model-level sfc geopotential')
         ENDIF
         lhave_ps_geop = (lhave_ps .and. lhave_geop)

         ! Check if pressure and temperature are available:
         lhave_pres = (test_cdi_varID(fileID_latbc, 'PRES', latbc_dict) /= -1)
         lhave_temp = (test_cdi_varID(fileID_latbc, 'TEMP', latbc_dict) /= -1)



         ! --- DEFINE WHICH VARIABLES SHALL BE READ FROM FILE ---

         buffer%lread_hhl         = .FALSE.
         buffer%lread_theta_rho   = .FALSE.
         buffer%lread_pres        = .FALSE.
         buffer%lread_temp        = .FALSE.
         buffer%lread_ps_geop     = .FALSE.
         buffer%lconvert_omega2w  = .FALSE.
         buffer%lcompute_hhl_pres = .FALSE.

         IF (lhave_hhl) THEN

           IF (lhave_theta_rho) THEN
             buffer%lread_hhl       = .TRUE.
             buffer%lread_theta_rho = .TRUE.
             !
           ELSE IF (lhave_pres .AND. lhave_temp) THEN
             buffer%lread_hhl       = .TRUE.
             buffer%lread_pres      = .TRUE.
             buffer%lread_temp      = .TRUE.
             !
           ELSE
             CALL finish(routine, "Non-hydrostatic LATBC data set, but neither RHO+THETA_V nor P,T provided!")
           END IF

         ELSE
           
           IF (lhave_ps_geop .AND. lhave_temp) THEN
             buffer%lread_temp        = .TRUE.
             buffer%lread_ps_geop     = .TRUE.
             buffer%lconvert_omega2w  = .TRUE.
             buffer%lcompute_hhl_pres = .TRUE.
           ELSE
             CALL finish(routine, "Hydrostatic LATBC data set, but PS,GEOP,T not provided!")
           END IF

         END IF

         buffer%lread_vn  = .FALSE.
         buffer%lread_u_v = .FALSE.
         IF (lhave_vn) THEN
           buffer%lread_vn = .TRUE.
         ELSE
           IF (lhave_u .AND. lhave_v) THEN
             buffer%lread_u_v = .TRUE.
           ELSE
             CALL finish(routine, "No VN or U&V available in LATBC data set!")
           END IF
         END IF


         !
         ! Consistency checks
         ! 
         IF (latbc_config%init_latbc_from_fg .AND. .NOT. buffer%lread_hhl) THEN
           CALL finish(routine, "Init LATBC from first guess requires BCs from non-hydrostatic model!")
         END IF


         !
         ! Write some status output:
         !
         IF (buffer%lread_theta_rho) THEN
           CALL message(routine,'Prognostic thermodynamic variables (RHO and THETA_V) are read from file.')
         ENDIF

         IF (.NOT. buffer%lread_qr) THEN
           CALL message(routine,'Rain water (QR) not available in input data.')
         ENDIF

         IF (.NOT. buffer%lread_qs) THEN
            CALL message(routine,'Snow water (QS) not available in input data.')
         ENDIF

         IF (buffer%lread_hhl) THEN
            CALL message(routine,'Input levels (HHL) are read from file.')
         ELSE
            CALL message(routine,'Input levels (HHL) are computed from sfc geopotential.')
         ENDIF

         IF (.NOT. buffer%lread_w) THEN
           CALL message(routine, "Neither W nor OMEGA provided! W is set to zero at LBCs")
         ELSE IF (buffer%lconvert_omega2w) THEN
            CALL message(routine,'Compute W from OMEGA.')
         ENDIF

         IF (buffer%lread_ps_geop) THEN
           CALL message(routine,'PS and GEOP are read from file.')
         ELSE IF (lhave_ps_geop) THEN
           CALL message(routine,'PS and GEOP are ignored.')
         END IF

         IF (buffer%lcompute_hhl_pres) THEN
           CALL message(routine,'HHL and PRES are computed based on PS and GEOP.')
         END IF

         IF (buffer%lread_pres) THEN
           CALL message(routine,'PRES is read from file.')
         ELSE
           CALL message(routine,'PRES is diagnosed.')
         END IF

         IF (buffer%lread_temp) THEN
           CALL message(routine,'TEMP is read from file.')
         ELSE
           CALL message(routine,'TEMP is diagnosed.')
         END IF

         IF (buffer%lread_vn) THEN
           CALL message(routine,'VN is read from file.')
         ELSE
           CALL message(routine,'U,V are read from file.')
         END IF


#endif

    END SUBROUTINE check_variables


    !-------------------------------------------------------------------------------------------------
    !> Replicates data needed for async prefetch on the prefetch proc.
    !  ATTENTION: The data is not completely replicated, only as far as needed for prefetch.
    !
    !  This routine has to be called by all PEs (work and prefetch)
    !
#ifndef NOMPI
    SUBROUTINE replicate_data_on_pref_proc(var_data, bcast_root)
      TYPE(t_var_data), ALLOCATABLE, INTENT(out) :: var_data(:)
      INTEGER,             INTENT(IN) :: bcast_root

      ! local variables
      CHARACTER(len=*), PARAMETER   :: routine = modname//"::replicate_data_on_pref_proc"
      INTEGER                       :: info_size, i, iv, nelems, nv, n, &
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
            DO WHILE (ASSOCIATED(element))
               !      IF(element%field%info%used_dimensions(2) == 0) CYCLE
               nvars = nvars+1
               element => element%next_list_element
            ENDDO
            all_nvars = all_nvars + nvars
         ENDDO
      ENDIF

      ! get the number of var lists
      CALL p_bcast(all_nvars, bcast_root, p_comm_work_2_pref)

      IF (all_nvars <= 0) RETURN

      ! allocate the array of variables
      ALLOCATE(var_data(all_nvars))

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
            DO WHILE (ASSOCIATED(element))
               nelems = nelems+1
               element => element%next_list_element
            ENDDO
         ENDIF

         ! Send basic info:
         CALL p_bcast(nelems, bcast_root, p_comm_work_2_pref)

         IF (is_pref) THEN
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
            DO WHILE (ASSOCIATED(element))
               i2 = i2 + 1
               var_data(i2)%info => element%field%info
               nelems = nelems+1
               info_storage(:,nelems) = TRANSFER(element%field%info, (/ 0 /))
               element => element%next_list_element
            ENDDO
         ENDIF

         ! Send binary representation of all info members

         CALL p_bcast(info_storage, bcast_root, p_comm_work_2_pref)

         IF (is_pref) THEN

            ! Insert elements into var list
           IF (nelems > 0) THEN
             ALLOCATE(p_var_list%p%first_list_element)
             element => p_var_list%p%first_list_element
             DO n = 1, nelems-1
               i2 = i2 + 1

               ! Set info structure from binary representation in info_storage
               element%field%info = TRANSFER(info_storage(:, n), info)
               var_data(i2)%info => element%field%info
               ALLOCATE(element%next_list_element)
               element => element%next_list_element
             ENDDO
             i2 = i2 + 1
             element%field%info = TRANSFER(info_storage(:, nelems), info)
             var_data(i2)%info => element%field%info
             NULLIFY(element%next_list_element)
           ELSE
             NULLIFY(p_var_list%p%first_list_element)
           END IF
         ENDIF
         DEALLOCATE(info_storage, STAT=ierrstat)
         IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")
      ENDDO

    END SUBROUTINE replicate_data_on_pref_proc


    !------------------------------------------------------------------------------------------------
    !> Sets the reorder_data for cells/edges/verts
    !  ATTENTION: This routine must only be called on work and test PE
    !             (i.e. not on prefetching PEs)
    !             The arguments don't make sense on the prefetching PE anyways
    !
    SUBROUTINE set_reorder_data(n_points, owner_mask, glb_index, ri)

      INTEGER, INTENT(IN) :: n_points        ! Local number of cells/edges/verts in logical patch
      LOGICAL, INTENT(IN) :: owner_mask(n_points) ! owner_mask for logical patch
      INTEGER, INTENT(IN) :: glb_index(:)    ! glb_index for logical patch
      TYPE(t_reorder_info), INTENT(INOUT) :: ri ! Result: reorder info

      ! local variables
      INTEGER :: i, n, ierrstat

      CHARACTER(LEN=*), PARAMETER :: routine = modname//"::set_reorder_data"

      ! Just for safety
      IF(my_process_is_pref()) CALL finish(routine, 'Must not be called on Prefetching PE')

      ! Get number of owned cells/edges/verts (without halos, physical patch only)
      ri%n_own = COUNT(owner_mask(:))

      !   WRITE(*,*) 'set_reorder_data p_pe_work ', p_pe_work , ' ri%n_own ', ri%n_own, ' n_points ', n_points

      ! Set index arrays to own cells/edges/verts
      ALLOCATE(ri%own_idx(ri%n_own), &
        &      ri%own_blk(ri%n_own), &
        &      STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")

      ! Global index of my own points
      ALLOCATE(ri%reorder_index(ri%n_own), STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")

      n = 0
      DO i = 1, n_points
         IF(owner_mask(i)) THEN
            n = n+1
            ri%own_idx(n) = idx_no(i)
            ri%own_blk(n) = blk_no(i)
            ri%reorder_index(n)  = glb_index(i)
         ENDIF
      ENDDO

      ! Get global number of points for current (physical!) patch
      ri%n_glb = -1

    END SUBROUTINE set_reorder_data


    !--------------------------------------------------------------------------
    !
    ! Transfers reorder data to restart PEs.
    !
    SUBROUTINE transfer_reorder_data(ri)
      TYPE(t_reorder_info), INTENT(INOUT) :: ri

      ! local variables
      INTEGER                             :: ierrstat, dummy(1), i, accum
      LOGICAL                             :: is_pref
      INTEGER, ALLOCATABLE                :: rcounts(:)
      CHARACTER(LEN=*), PARAMETER :: routine = modname//"::transfer_reorder_data"


      ! Gather the number of own points for every PE into ri%pe_own

      is_pref = my_process_is_pref()
      IF (is_pref) THEN

        ! on prefetch PE: n_own=0, own_ide and own_blk are not allocated
        ri%n_own = 0

        ! pe_own must be allocated for num_work_procs, not for p_n_work
        ALLOCATE(ri%pe_own(0:num_work_procs-1), &
             ri%pe_off(0:num_work_procs-1), STAT=ierrstat)
        IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")
      ENDIF

      IF (is_pref) THEN
        ! Gather the number of own points for every PE into ri%pe_own
        CALL p_allgather(dummy(1), ri%pe_own, sendcount=0, recvcount=1, &
             comm=p_comm_work_2_pref)

        ! Get offset within result array
        accum = 0
        DO i = 0, num_work_procs-1
          ri%pe_off(i) = accum
          accum = accum + ri%pe_own(i)
        ENDDO
        ri%n_glb = accum
        ALLOCATE(ri%reorder_index(accum), STAT=ierrstat)
        IF(ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")

        CALL p_allgatherv(dummy(1:0), ri%reorder_index, ri%pe_own, &
             ri%pe_off, comm=p_comm_work_2_pref)
      ELSE
        CALL p_allgather(ri%n_own, dummy, sendcount=1, recvcount=0, &
             comm=p_comm_work_2_pref)
        ALLOCATE(rcounts(num_prefetch_proc), STAT=ierrstat)
        IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")
        rcounts = 0
        CALL p_allgatherv(ri%reorder_index, dummy, rcounts, &
             rcounts, comm=p_comm_work_2_pref)
      END IF

    END SUBROUTINE  transfer_reorder_data

    SUBROUTINE match_var_data_to_buffer(map, buffer, var_data, StrLowCasegrp)
      TYPE(t_buffer), INTENT(in) :: buffer
      INTEGER, ALLOCATABLE, INTENT(out) :: map(:)
      CHARACTER(LEN=varname_len), INTENT(in) :: StrLowCasegrp(buffer%ngrp_vars)
      TYPE(t_var_data), INTENT(in) :: var_data(:)

      LOGICAL :: grp_vars_bool(buffer%ngrp_vars)
      INTEGER :: grp_tlen(buffer%ngrp_vars)
      INTEGER :: iv, jp, nlevs, ndims

      ALLOCATE(map(buffer%ngrp_vars))
      map = -1
      DO jp = 1, buffer%ngrp_vars
        grp_tlen(jp) = LEN_TRIM(StrLowCasegrp(jp))
      END DO
      grp_vars_bool = .FALSE.

      DO iv = 1, SIZE(var_data)
        ndims = var_data(iv)%info%ndims
        DO jp = 1, buffer%ngrp_vars
          ! Use only the variables of time level 1 (".TL1") to determine memory sizes.
          IF (       .NOT. grp_vars_bool(jp) &
               .AND. (StrLowCasegrp(jp) == var_data(iv)%info%name &
               &      .OR. StrLowCasegrp(jp)(1:grp_tlen(jp))//TIMELEVEL_SUFFIX//'1' == var_data(iv)%info%name)) THEN
            nlevs = MERGE(1, buffer%nlev(jp), ndims == 2)
            IF (nlevs /= 0) THEN
              map(jp) = iv
              grp_vars_bool(jp) = .TRUE.
            END IF
          END IF
        END DO
      END DO
    END SUBROUTINE match_var_data_to_buffer

    SUBROUTINE infer_buffer_hgrid(buffer, var_data, map)
      TYPE(t_buffer), INTENT(inout) :: buffer
      TYPE(t_var_data), INTENT(in) :: var_data(:)
      INTEGER, INTENT(in) :: map(buffer%ngrp_vars)

      INTEGER :: jp, iv

      ALLOCATE(buffer%hgrid(buffer%ngrp_vars))
      DO jp = 1, buffer%ngrp_vars
        iv = map(jp)
        IF (iv /= -1) THEN
          buffer%hgrid(jp) = var_data(iv)%info%hgrid
        ELSE
          buffer%hgrid(jp) = -1
        END IF
      END DO
      DO jp = 1, buffer%ngrp_vars
        IF (buffer%mapped_name(jp) == buffer%geop_ml_var) THEN
          ! variable GEOSP is stored in cell center location
          buffer%hgrid(jp) = GRID_UNSTRUCTURED_CELL
        END IF
      END DO
    END SUBROUTINE infer_buffer_hgrid
    !------------------------------------------------------------------------------------------------
    !> Initializes the remote memory window for asynchronous prefetch.
    !
    SUBROUTINE create_client_buffers(buffer, var_data, n_own_cells, nblks_c, &
         n_own_edges, nblks_e, map, mem_size)
      TYPE(t_buffer), INTENT(inout) :: buffer
      TYPE(t_var_data), INTENT(in) :: var_data(:)
      INTEGER, INTENT(in) :: n_own_cells, n_own_edges, nblks_c, nblks_e, map(:)
      INTEGER(kind=mpi_address_kind), INTENT(out) :: mem_size

      INTEGER :: ierror, iv, jp, nlevs, nblks, hgrid
      INTEGER(kind=mpi_address_kind) :: mem_size_, lvlsz
      LOGICAL :: is_work, found_unknown_grid
      CHARACTER(LEN=*), PARAMETER :: routine = modname//"::init_memory_window"

      ! Get size and offset of the data for the input
      mem_size_ = 0_mpi_address_kind

      is_work = my_process_is_work()
      found_unknown_grid = .FALSE.
      DO jp = 1, buffer%ngrp_vars
        iv = map(jp)
        IF (iv /= -1) THEN
          IF (var_data(iv)%info%ndims == 2) THEN
            nlevs = 1
          ELSE
            nlevs = buffer%nlev(jp)
          ENDIF

          IF (nlevs /= 0) THEN

            hgrid = buffer%hgrid(jp)

            SELECT CASE (hgrid)
            CASE (GRID_UNSTRUCTURED_CELL)
              ! variable stored in cell center location
              lvlsz = INT(n_own_cells, mpi_address_kind)
              nblks = nblks_c
            CASE (GRID_UNSTRUCTURED_EDGE)
              ! variable stored in edge center of a cell
              lvlsz = INT(n_own_edges, mpi_address_kind)
              nblks = nblks_e
            CASE DEFAULT
              found_unknown_grid = .TRUE.
              EXIT
            END SELECT

            mem_size_ = mem_size_ + INT(nlevs, mpi_address_kind) * lvlsz
            IF(is_work)THEN
              ! allocate the buffer sizes for variables on compute processors
              ALLOCATE(buffer%vars(jp)%buffer(nproma, nlevs, nblks), STAT=ierror)
              IF (ierror /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")
            ENDIF
          ENDIF
        ENDIF
      ENDDO ! vars
      IF (found_unknown_grid) CALL finish(routine,'Unknown grid type found!')

      DO jp = 1, buffer%ngrp_vars
        IF (buffer%mapped_name(jp) == buffer%geop_ml_var) THEN
          ! Memory for GEOSP variable taken as memory equivalent to 1 level of z_ifc
          ! as the variable GEOSP doesn't exist in metadata
          mem_size_ = mem_size_ + INT(1*n_own_cells,mpi_address_kind)

          IF(is_work)THEN
            ! allocate the buffer sizes for variable 'GEOSP' on compute processors
            ALLOCATE(buffer%vars(jp)%buffer(nproma, 1, nblks_c), STAT=ierror)
            IF (ierror /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")
          ENDIF

        ENDIF
      ENDDO
      mem_size = mem_size_

    END SUBROUTINE create_client_buffers


    !------------------------------------------------------------------------------------------------
    !> allocate amount of memory needed with MPI_Alloc_mem
    !
    !  @note Implementation for non-Cray pointers
    !
    SUBROUTINE create_win(mem_win, mem_size)

      TYPE(t_mem_win), INTENT(INOUT) :: mem_win
      INTEGER (KIND=MPI_ADDRESS_KIND), INTENT(IN)    :: mem_size

      ! local variables
      CHARACTER(LEN=*), PARAMETER :: routine = modname//"::allocate_mem_noncray"
      TYPE(c_ptr)                     :: c_mem_ptr
      INTEGER                         :: ierror, nbytes_real
      INTEGER (KIND=MPI_ADDRESS_KIND) :: mem_bytes
      REAL(sp), TARGET :: dummy(1)
      ! Get the amount of bytes per REAL*4 variable (as used in MPI
      ! communication)
      CALL MPI_Type_extent(p_real_sp, nbytes_real, ierror)

      ! For the IO PEs the amount of memory needed is 0 - allocate at least 1 word there:
      IF (mem_size > 0) THEN
        mem_bytes = mem_size*INT(nbytes_real,mpi_address_kind)

        CALL MPI_Alloc_mem(mem_bytes, MPI_INFO_NULL, c_mem_ptr, ierror)

        CALL C_F_POINTER(c_mem_ptr, mem_win%mem_ptr_sp, (/ mem_size /) )
        mem_win%mem_ptr_sp(:) = 0._sp
      ELSE
        mem_bytes = 0_mpi_address_kind
        mem_win%mem_ptr_sp => dummy
      END IF

      ! Create memory window for communication
      CALL MPI_Win_create(mem_win%mem_ptr_sp, mem_bytes, nbytes_real, &
        &                 MPI_INFO_NULL, p_comm_work_pref, mem_win%mpi_win, &
        &                 ierror )
      IF (ierror /= 0) CALL finish(routine, "MPI error!")

    END SUBROUTINE create_win
#endif

END MODULE mo_async_latbc
