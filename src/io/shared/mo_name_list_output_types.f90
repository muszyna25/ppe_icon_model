!>
!! Module containing type definitions and constants for synchronous
!! and asynchronous output.
!!
!! @author R. Johanni
!!
!! @par Revision History
!! Initial implementation  by  R. Johanni  (2011)
!! Major changes: F. Prill, DWD (2012-2013)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
!! @note: The spelling "name_list" (with underscore) is intended to make
!!        clear that this does not pertain to a FORTRAN namelist but rather
!!        to a list of names of output variables
!!
MODULE mo_name_list_output_types

  USE mo_kind,                  ONLY: wp, dp, sp
  USE mo_impl_constants,        ONLY: max_phys_dom, vname_len,                         &
    &                                 max_var_ml, max_var_pl, max_var_hl, max_var_il,  &
    &                                 MAX_TIME_LEVELS, MAX_NUM_IO_PROCS,               &
    &                                 MAX_TIME_INTERVALS, MAX_CHAR_LENGTH, MAX_NPLEVS, &
    &                                 MAX_NZLEVS, MAX_NILEVS
  USE mo_io_units,              ONLY: filename_max
  USE mo_var_metadata_types,    ONLY: t_var_metadata
  USE mo_util_uuid,             ONLY: t_uuid
  USE mo_util_string,           ONLY: tolower
  USE mo_communication,         ONLY: t_comm_gather_pattern
  USE mtime,                    ONLY: MAX_DATETIME_STR_LEN, MAX_TIMEDELTA_STR_LEN
  USE mo_output_event_types,    ONLY: t_par_output_event, MAX_EVENT_NAME_STR_LEN

  IMPLICIT NONE

  PRIVATE
  ! constants:
  PUBLIC :: l_output_phys_patch,max_z_axes
  PUBLIC :: REMAP_NONE
  PUBLIC :: REMAP_REGULAR_LATLON
  PUBLIC :: msg_io_start
  PUBLIC :: msg_io_done
  PUBLIC :: msg_io_shutdown
  PUBLIC :: IRLON, IRLAT, ILATLON
  PUBLIC :: ICELL, IEDGE, IVERT
  PUBLIC :: GRP_PREFIX
  PUBLIC :: TILE_PREFIX
  PUBLIC :: GRB2_GRID_INFO_NAME, GRB2_GRID_INFO

  ! derived data types:
  PUBLIC :: t_mem_win
  PUBLIC :: t_reorder_info
  PUBLIC :: t_grid_info
  PUBLIC :: t_patch_info
  PUBLIC :: t_patch_info_ll
  PUBLIC :: t_output_name_list
  PUBLIC :: t_rptr_5d
  PUBLIC :: t_iptr_5d
  PUBLIC :: t_var_desc
  PUBLIC :: t_fname_metadata
  PUBLIC :: t_level_selection
  PUBLIC :: t_output_file
  ! global variables
  PUBLIC :: all_events
  ! utility subroutines
  PUBLIC :: is_grid_info_var


  !------------------------------------------------------------------------------------------------
  ! CONSTANTS
  !------------------------------------------------------------------------------------------------

  ! prefix for group identifier in output namelist
  CHARACTER(len=6), PARAMETER :: GRP_PREFIX = "group:"
  ! prefix for tile-group identifier in output namelist
  CHARACTER(len=6), PARAMETER :: TILE_PREFIX = "tiles:"

  ! Tags for communication between compute PEs and I/O PEs
  INTEGER, PARAMETER :: msg_io_start    = 12345
  INTEGER, PARAMETER :: msg_io_done     = 54321
  INTEGER, PARAMETER :: msg_io_shutdown = 99999

  ! constants for better readability:
  INTEGER, PARAMETER :: REMAP_NONE            = 0
  INTEGER, PARAMETER :: REMAP_REGULAR_LATLON  = 1

  INTEGER, PARAMETER :: IRLON                 = 1
  INTEGER, PARAMETER :: IRLAT                 = 2
  INTEGER, PARAMETER :: ILATLON               = 1
  INTEGER, PARAMETER :: ICELL                 = 1
  INTEGER, PARAMETER :: IEDGE                 = 2
  INTEGER, PARAMETER :: IVERT                 = 3

  ! The following parameter decides whether physical or logical patches are output
  ! and thus whether the domain number in output name lists pertains to physical
  ! or logical patches.
  LOGICAL, PARAMETER :: l_output_phys_patch = .TRUE. !** DO NOT CHANGE - needed for GRIB output **!

  INTEGER, PARAMETER :: max_z_axes = 43

  ! Character-strings denoting the "special" GRIB2 output fields that
  ! describe the grid coordinates. These fields are ignored by most
  ! output routines and only used by "set_grid_info_grb2":
  CHARACTER(LEN=5), PARAMETER :: GRB2_GRID_INFO = "GRID:"
  CHARACTER(LEN=9), PARAMETER :: GRB2_GRID_INFO_NAME(0:3,2) = &
    &  RESHAPE( (/ "GRID:RLON", "GRID:CLON", "GRID:ELON", "GRID:VLON", &
    &              "GRID:RLAT", "GRID:CLAT", "GRID:ELAT", "GRID:VLAT" /), (/4,2/) )

  !------------------------------------------------------------------------------------------------
  ! DERIVED DATA TYPES
  !------------------------------------------------------------------------------------------------

  !> Data structure containing info on computational mesh (for output)
  !
  TYPE t_grid_info
    ! only used when copying grid info from file (grid_info_mode = GRID_INFO_BCAST):
    REAL(wp), ALLOCATABLE :: lon   (:), lat   (:)
    REAL(wp), ALLOCATABLE :: lonv(:,:), latv(:,:)

    ! only used when copying grid info from file (grid_info_mode = GRID_INFO_FILE):
    INTEGER,  ALLOCATABLE :: log_dom_index(:)
    ! Index where a point of the physical domains is in the logical domain
  END TYPE t_grid_info


  !------------------------------------------------------------------------------------------------
  ! TYPE t_reorder_info describes how local cells/edges/verts
  ! have to be reordered to get the global array.
  ! Below, "points" refers to either cells, edges or verts.
  !
  ! TODO[FP] Note that the "reorder_info" contains fields of *global*
  !          size (reorder_index). On the compute PEs these fields
  !          could be deallocated after the call to
  !          "transfer_reorder_info" in the setup phase!

  TYPE t_reorder_info
    INTEGER                    :: n_glb  ! Global number of points per physical patch
    INTEGER                    :: n_log  ! Global number of points in the associated logical patch
    INTEGER                    :: n_own  ! Number of own points (without halo, only belonging to phyiscal patch)
    ! Only set on compute PEs, set to 0 on IO PEs
    INTEGER, ALLOCATABLE       :: own_idx(:), own_blk(:)
    ! idx and blk for own points, only set on compute PEs
    INTEGER, ALLOCATABLE       :: own_dst_idx(:), own_dst_blk(:)
    ! dest idx and blk for own points, only set on sequential/test PEs
    INTEGER, ALLOCATABLE       :: pe_own(:)
    ! n_own, gathered for all compute PEs (set on all PEs)
    INTEGER, ALLOCATABLE       :: pe_off(:)
    ! offset of contributions of PEs (set on all PEs)
    INTEGER, ALLOCATABLE       :: reorder_index(:)
    ! Index how to reorder the contributions of all compute PEs
    ! into the global array (set on all PEs)

    ! grid information: geographical locations of cells, edges, and
    ! vertices which is first collected on working PE 0 - from where
    ! it will be broadcasted to the pure I/O PEs.
    TYPE (t_grid_info)         :: grid_info
  END TYPE t_reorder_info


  ! TYPE t_patch_info contains the reordering info for cells, edges and verts
  TYPE t_patch_info
    TYPE(t_reorder_info)                 :: cells
    TYPE(t_reorder_info)                 :: edges
    TYPE(t_reorder_info)                 :: verts
    INTEGER                              :: log_patch_id

    ! pointer to communication pattern for GATHER operation;
    ! corresponds to physical or logical patch, depending on
    ! "l_output_phys_patch"
    TYPE(t_comm_gather_pattern), POINTER ::  p_pat_c, p_pat_v, p_pat_e

    ! global number of points, corresponds to physical or logical
    ! patch, depending on "l_output_phys_patch"
    INTEGER                              :: nblks_glb_c, nblks_glb_v, nblks_glb_e

    ! Filename of grid file, needed only if grid information is output
    ! to NetCDF since this information is normally not read and
    ! thus not present in the patch description
    CHARACTER(LEN=filename_max)          :: grid_filename

    ! uuid of grid (provided by grid file)
    TYPE(t_uuid)                         :: grid_uuid

    ! Number of grid used (provided by grid file)
    INTEGER                              :: number_of_grid_used

    ! mode how to collect grid information (for output)
    INTEGER                              :: grid_info_mode

    ! the maximum cell connectivity reproduced from the patch
    INTEGER                              :: max_cell_connectivity
    INTEGER                              :: max_vertex_connectivity
  END TYPE t_patch_info


  !> Reordering info for regular (lon-lat) grids
  TYPE t_patch_info_ll
    TYPE(t_reorder_info)                 :: ri
    ! mode how to collect grid information (for output)
    INTEGER                              :: grid_info_mode
  END TYPE t_patch_info_ll


  TYPE t_output_name_list
    ! --------------------
    ! file name and format
    ! --------------------

    INTEGER                               :: filetype          ! One of CDI's FILETYPE_XXX constants
    CHARACTER(LEN=filename_max)           :: output_filename   ! output filename prefix
    CHARACTER(LEN=filename_max)           :: filename_format   ! output filename format (contains keywords <physdom>,<levtype> etc.)
    CHARACTER(LEN=filename_max)           :: filename_extn     ! user-specified filename extension (or "default")

    ! --------------------
    ! general settings
    ! --------------------

    INTEGER                               :: mode              ! 1 = forecast mode, 2 = climate mode
    INTEGER                               :: dom(max_phys_dom) ! domains for which this namelist is used, ending with -1
    INTEGER                               :: steps_per_file    ! Max number of output steps in one output file
    LOGICAL                               :: steps_per_file_inclfirst !< Flag. Do not count first step in files count
    CHARACTER(LEN=MAX_TIMEDELTA_STR_LEN)  :: file_interval     ! length of a file (ISO8601 duration)
    LOGICAL                               :: include_last      ! Flag whether to include the last timestep in output
    LOGICAL                               :: output_grid       ! Flag whether grid information is output (in NetCDF output)

    INTEGER                               :: taxis_tunit       ! 1 = TUNIT_SECOND, 2 = TUNIT_MINUTE, 3 TUNIT_HOUR ... (see cdi.inc)

    !> There are two alternative implementations for setting the
    !  output intervals, "output_bounds" and "output_start" /
    !  "output_end" / "output_interval". The former defines the output
    !  events relative to the simulation start (in seconds) and the
    !  latter define the output events by setting ISO8601-conforming
    !  date-time strings.
    REAL(wp)                              :: output_bounds(3*MAX_TIME_INTERVALS)

    !> Output event steps happen at regular intervals. These are given
    !  by an interval size and the time stamps for begin and end.
    CHARACTER(LEN=MAX_DATETIME_STR_LEN)   :: output_start(MAX_TIME_INTERVALS),     &
      &                                      output_end(MAX_TIME_INTERVALS),       &
      &                                      output_interval(MAX_TIME_INTERVALS)

    !> Additional days to be added to the output interval. This
    !  namelist parameter is required when the output interval has
    !  been provided as a number of seconds, e.g. which is so large
    !  that it cannot be converted into a valid ISO duration string.
    INTEGER :: additional_days(MAX_TIME_INTERVALS)

    !> ready filename prefix (=output event name)
    CHARACTER(LEN=MAX_EVENT_NAME_STR_LEN) :: ready_file

    ! --------------------
    ! variable lists
    ! --------------------

    CHARACTER(LEN=vname_len)  :: ml_varlist(max_var_ml)   ! name of model level fields
    CHARACTER(LEN=vname_len)  :: pl_varlist(max_var_pl)   ! name of pressure level fields
    CHARACTER(LEN=vname_len)  :: hl_varlist(max_var_hl)   ! name of height level fields
    CHARACTER(LEN=vname_len)  :: il_varlist(max_var_il)   ! name of isentropic level fields

    !> "stream_partitions": Split one namelist into concurrent,
    !> alternating files:
    INTEGER                   :: stream_partitions_ml, &
      &                          stream_partitions_pl, &
      &                          stream_partitions_hl, &
      &                          stream_partitions_il

    !> MPI ranks which were explicitly specified by the user:
    INTEGER                   :: pe_placement_ml(MAX_NUM_IO_PROCS), &
      &                          pe_placement_pl(MAX_NUM_IO_PROCS), &
      &                          pe_placement_hl(MAX_NUM_IO_PROCS), &
      &                          pe_placement_il(MAX_NUM_IO_PROCS)

    ! --------------------
    ! horizontal interpol.
    ! --------------------

    INTEGER  :: remap                 ! interpolate horizontally, 0: none, 1: to regular lat-lon grid, 2: to Gaussian grids, (3:...)
    LOGICAL  :: remap_internal        ! do interpolations online in the model or external (including triggering)
    INTEGER  :: lonlat_id             ! if remap=1: index of lon-lat-grid in global list "lonlat_grid_list"

    ! --------------------
    ! vertical interpol.
    ! --------------------

    CHARACTER(len=MAX_CHAR_LENGTH) :: m_levels  ! model levels (indices)

    REAL(wp) :: p_levels(MAX_NPLEVS) ! pressure levels
    REAL(wp) :: z_levels(MAX_NZLEVS) ! height levels
    REAL(wp) :: i_levels(MAX_NILEVS) ! isentropic levels

    ! -------------------------------------
    ! data operations
    ! ------------------------------------
    CHARACTER(len=MAX_CHAR_LENGTH) :: operation ! "mean"

    ! -------------------------------------
    ! Internal members, not read from input
    ! -------------------------------------

    TYPE(t_output_name_list), POINTER :: next ! Pointer to next output_name_list
  END TYPE t_output_name_list


  ! Unfortunately, Fortran does not allow arrays of pointers, so we
  ! have to define extra types
  TYPE t_rptr_5d
    REAL(wp), POINTER :: p(:,:,:,:,:)
  END TYPE t_rptr_5d

  TYPE t_sptr_5d
    REAL(sp), POINTER :: p(:,:,:,:,:)
  END TYPE t_sptr_5d

  TYPE t_iptr_5d
    INTEGER,  POINTER :: p(:,:,:,:,:)
  END TYPE t_iptr_5d


  TYPE t_var_desc
    REAL(wp), POINTER                     :: r_ptr(:,:,:,:,:)                 !< Pointer to time level independent REAL data (or NULL)
    REAL(sp), POINTER                     :: s_ptr(:,:,:,:,:)                 !< Pointer to time level independent REAL(sp) data (or NULL)
    INTEGER,  POINTER                     :: i_ptr(:,:,:,:,:)                 !< Pointer to time level independent INTEGER data (or NULL)
    TYPE(t_rptr_5d)                       :: tlev_rptr(MAX_TIME_LEVELS)       !< Pointers to time level dependent REAL data
    TYPE(t_sptr_5d)                       :: tlev_sptr(MAX_TIME_LEVELS)       !< Pointers to time level dependent REAL(sp) data
    TYPE(t_iptr_5d)                       :: tlev_iptr(MAX_TIME_LEVELS)       !< Pointers to time level dependent INTEGER data
    TYPE(t_var_metadata), POINTER         :: info_ptr                         !< Pointer to the info structure of the variable

    !> Info structure for variable: this is a modified copy of the
    !> variable's "info" data object!
    TYPE(t_var_metadata)                  :: info                             
  END TYPE t_var_desc


  !> Data structure containing variables for MPI memory window
  !
  TYPE t_mem_win
    INTEGER                               :: mpi_win                          !< MPI window for data communication
    INTEGER                               :: mpi_win_metainfo                 !< MPI window for metadata
    REAL(dp), POINTER                     :: mem_ptr_dp(:)                    !< Pointer to memory window (REAL*8)
    REAL(sp), POINTER                     :: mem_ptr_sp(:)                    !< Pointer to memory window (REAL*4)
    INTEGER,  POINTER                     :: mem_ptr_metainfo_pe0(:)          !< Pointer to variable meta-info.
  END TYPE t_mem_win


  !> Data structure containing additional meta-data for generating
  !> an output filename.
  !
  TYPE t_fname_metadata
    INTEGER                               :: steps_per_file                   !< (optional:) no. of output steps per file
    LOGICAL                               :: steps_per_file_inclfirst         !< Flag. Do not count first step in files count
    CHARACTER(LEN=MAX_TIMEDELTA_STR_LEN)  :: file_interval                    !< (optional:) length of a file (ISO8601 duration)
    INTEGER                               :: phys_patch_id                    !< ID of physical output patch
    INTEGER                               :: ilev_type                        !< level_type_ml/_pl/_hl/_il
    CHARACTER(LEN=FILENAME_MAX)           :: filename_format                  !< output filename format (contains keywords)
    CHARACTER(LEN=FILENAME_MAX)           :: filename_pref                    !< Prefix of output file name
    CHARACTER(LEN=16)                     :: extn                             !< filename extension
    INTEGER                               :: jfile_offset                     !< offset for filename numbers

    !> "stream_partitions": Split one namelist into concurrent,
    !> alternating files:
    INTEGER                               :: npartitions                      !< total no. of stream partitions
    INTEGER                               :: ifile_partition                  !< this file's partition index

  END TYPE t_fname_metadata


  TYPE t_level_selection
    !> number of selected levels
    INTEGER :: n_selected

    !> level selection as input in the form of a LOGICAL array
    !  s(1...N), where "s(i)=.TRUE." means that level "i" is selected:
    LOGICAL, ALLOCATABLE :: s(:)

    !> integer list idx(1...n_selected) containing the selected level
    !  indices.
    INTEGER, ALLOCATABLE :: global_idx(:)

    !> local index in the list of selected level indices (i.e. an
    !  integer number in the range 1...n_selected).
    INTEGER, ALLOCATABLE :: local_idx(:)
  END TYPE t_level_selection


  TYPE t_output_file
    ! The following data must be set before opening the output file:

    TYPE(t_par_output_event), POINTER     :: out_event                        !< data structure for parallel output events

    CHARACTER(LEN=filename_max)           :: filename_pref                    !< Prefix of output file name
    INTEGER                               :: output_type                      !< CDI format
    INTEGER                               :: phys_patch_id                    !< ID of physical output patch
    INTEGER                               :: log_patch_id                     !< ID of logical output patch
    INTEGER                               :: ilev_type                        !< level type: level_type_ml/level_type_pl/level_type_hl/level_type_il
    INTEGER                               :: max_vars                         !< maximum number of variables allocated
    INTEGER                               :: num_vars                         !< number of variables in use
    TYPE(t_var_desc),ALLOCATABLE          :: var_desc(:)
    TYPE(t_output_name_list), POINTER     :: name_list                        !< Pointer to corresponding output name list

    CHARACTER(LEN=vname_len), ALLOCATABLE :: name_map(:,:)                    !< mapping internal names -> names in NetCDF

    INTEGER                               :: remap                            !< Copy of remap from associated namelist
    INTEGER                               :: io_proc_id                       !< ID of process doing I/O on this file

    !> indices when one namelist has been split into concurrent,
    !> alternating files ("stream_partitions"):
    INTEGER                               :: npartitions                      !< total no. of parts
    INTEGER                               :: ifile_partition                  !< index of this file

    !> MPI rank which were explicitly specified by the user:
    INTEGER                               :: pe_placement

    ! Used for async IO only
    TYPE(t_mem_win)                       :: mem_win                          !< data structure containing variables for MPI memory window

    ! Selection of vertical levels (not necessarily present)
    TYPE (t_level_selection), POINTER     :: level_selection => NULL()        !< selection of vertical levels

    ! The following members are set during open
    INTEGER                               :: cdiFileId
    INTEGER                               :: cdiVlistId                       !< cdi vlist handler
    INTEGER                               :: cdiVlistId_orig                  !< cdi vlist handler, storing the model internal vlist id during append
    INTEGER                               :: cdiCellGridID
    INTEGER                               :: cdiSingleGridID
    INTEGER                               :: cdiVertGridID
    INTEGER                               :: cdiEdgeGridID
    INTEGER                               :: cdiLonLatGridID
    INTEGER                               :: cdiZaxisID(max_z_axes)           !< All types of possible Zaxis ID's
    INTEGER                               :: cdiTaxisID
    INTEGER                               :: cdiTaxisID_orig
    INTEGER                               :: cdiTimeIndex
    INTEGER                               :: cdiInstID                        !< output generating institute
    INTEGER                               :: cdi_grb2(3,2)                    !< geographical position: (GRID, latitude/longitude)
    LOGICAL                               :: appending = .FALSE.              !< the current file is appended (.true.), otherwise .false. 

  END TYPE t_output_file

  ! "all_events": The root I/O MPI rank "ROOT_OUTEVENT" asks all
  ! participating I/O PEs for their output event info and generates a
  ! unified output event, indicating which PE performs a write process
  ! at which step:
  TYPE(t_par_output_event), POINTER       :: all_events

CONTAINS

  !------------------------------------------------------------------------------------------------
  !> Utility routine: .TRUE. if the variable corresponds to a "grid
  !  info" variable like clon, clat, elon, elat, etc. which must be
  !  ignored by most output subroutines.
  !
  FUNCTION is_grid_info_var(varname)
    LOGICAL :: is_grid_info_var
    CHARACTER(len=*), INTENT(IN) :: varname
    ! local variables
    INTEGER :: idx

    idx = INDEX(TRIM(varname), TRIM(tolower(GRB2_GRID_INFO)))
    is_grid_info_var = (idx > 0)
  END FUNCTION is_grid_info_var

END MODULE mo_name_list_output_types
