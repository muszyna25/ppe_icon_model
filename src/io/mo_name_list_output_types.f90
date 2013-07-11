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
!! @note: The spelling "name_list" (with underscore) is intended to make
!!        clear that this does not pertain to a FORTRAN namelist but rather
!!        to a list of names of output variables
!!
MODULE mo_name_list_output_types

  USE mo_kind,                  ONLY: wp, i8
  USE mo_impl_constants,        ONLY: max_phys_dom, vname_len,                         &
    &                                 max_var_ml, max_var_pl, max_var_hl, max_var_il,  &
    &                                 MAX_TIME_LEVELS, max_bounds, max_levels
  USE mo_io_units,              ONLY: filename_max
  USE mo_var_metadata,          ONLY: t_var_metadata
  USE mo_util_uuid,             ONLY: t_uuid
  USE mo_communication,         ONLY: t_comm_pattern

  IMPLICIT NONE

  PRIVATE
  ! constants:
  PUBLIC :: l_output_phys_patch
  PUBLIC :: REMAP_NONE
  PUBLIC :: REMAP_REGULAR_LATLON
  PUBLIC :: msg_io_start
  PUBLIC :: msg_io_done
  PUBLIC :: msg_io_shutdown
  PUBLIC :: IRLON, IRLAT, ILATLON
  PUBLIC :: ICELL, IEDGE, IVERT
  PUBLIC :: sfs_name_list
  PUBLIC :: second_tos
  PUBLIC :: GRP_PREFIX
  ! derived data types:
  PUBLIC :: t_reorder_info
  PUBLIC :: t_grid_info
  PUBLIC :: t_patch_info
  PUBLIC :: t_output_name_list
  PUBLIC :: t_rptr_5d
  PUBLIC :: t_iptr_5d
  PUBLIC :: t_var_desc
  PUBLIC :: t_output_file


  !------------------------------------------------------------------------------------------------
  ! CONSTANTS
  !------------------------------------------------------------------------------------------------

  ! prefix for group identifier in output namelist
  CHARACTER(len=6), PARAMETER :: GRP_PREFIX = "group:"

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

  ! fields for which typeOfSecondFixedSurface must be re-set
  CHARACTER(LEN=12), PARAMETER :: sfs_name_list(6) =(/"z_ifc       ", "topography_c", &
    &                                                 "hbas_con    ", "htop_con    ", &
    &                                                 "hzerocl     ", "clcl        "/)
  ! typeOfSecondFixedSurface to be used
  INTEGER          , PARAMETER :: second_tos(6)    =(/101, 101, 101, 101, 101, 1/)

  ! The following parameter decides whether physical or logical patches are output
  ! and thus whether the domain number in output name lists pertains to physical
  ! or logical patches.
  LOGICAL, PARAMETER :: l_output_phys_patch = .TRUE. !** DO NOT CHANGE - needed for GRIB output **!


  !------------------------------------------------------------------------------------------------
  ! DERIVED DATA TYPES
  !------------------------------------------------------------------------------------------------

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
    INTEGER :: n_glb  ! Global number of points per physical patch
    INTEGER :: n_log  ! Global number of points in the associated logical patch
    INTEGER :: n_own  ! Number of own points (without halo, only belonging to phyiscal patch)
    ! Only set on compute PEs, set to 0 on IO PEs
    INTEGER, ALLOCATABLE :: own_idx(:), own_blk(:)
    ! idx and blk for own points, only set on compute PEs
    INTEGER, ALLOCATABLE :: own_dst_idx(:), own_dst_blk(:)
    ! dest idx and blk for own points, only set on sequential/test PEs
    INTEGER, ALLOCATABLE :: pe_own(:)
    ! n_own, gathered for all compute PEs (set on all PEs)
    INTEGER, ALLOCATABLE :: pe_off(:)
    ! offset of contributions of PEs (set on all PEs)
    INTEGER, ALLOCATABLE :: reorder_index(:)
    ! Index how to reorder the contributions of all compute PEs
    ! into the global array (set on all PEs)

    ! only used when copying grid info from file (l_grid_info_from_file = .TRUE.):
    INTEGER, ALLOCATABLE :: log_dom_index(:)
    ! Index where a point of the physical domains is in the logical domain
  END TYPE t_reorder_info


  TYPE t_grid_info
    REAL(wp), ALLOCATABLE :: lon   (:), lat   (:)
    REAL(wp), ALLOCATABLE :: lonv(:,:), latv(:,:)
  END TYPE t_grid_info


  ! TYPE t_patch_info contains the reordering info for cells, edges and verts
  TYPE t_patch_info
    TYPE(t_reorder_info) :: cells
    TYPE(t_reorder_info) :: edges
    TYPE(t_reorder_info) :: verts
    INTEGER :: log_patch_id

    ! pointer to communication pattern for GATHER operation;
    ! corresponds to physical or logical patch, depending on
    ! "l_output_phys_patch"
    TYPE(t_comm_pattern),  POINTER :: p_pat_c, p_pat_v, p_pat_e

    ! global number of points, corresponds to physical or logical
    ! patch, depending on "l_output_phys_patch"
    INTEGER :: nblks_glb_c, nblks_glb_v, nblks_glb_e

    ! grid information: geographical locations of cells, edges, and
    ! vertices which is first collected on working PE 0 - from where
    ! it will be broadcasted to the pure I/O PEs.
    TYPE (t_grid_info) :: grid_c, grid_e, grid_v

    ! Filename of grid file, needed only if grid information is output
    ! to NetCDF since this information is normally not read and
    ! thus not present in the patch description
    CHARACTER(LEN=filename_max) :: grid_filename

    ! uuid of grid (provided by grid file)
    TYPE(t_uuid) :: grid_uuid

    ! Number of grid used (provided by grid file)
    INTEGER :: number_of_grid_used
  END TYPE t_patch_info


  TYPE t_output_name_list
    ! --------------------
    ! file name and format
    ! --------------------

    INTEGER                     :: filetype          ! One of CDI's FILETYPE_XXX constants
    CHARACTER(LEN=filename_max) :: output_filename   ! output filename prefix
    CHARACTER(LEN=filename_max) :: filename_format   ! output filename format (contains keywords <physdom>,<levtype> etc.)

    ! --------------------
    ! general settings
    ! --------------------

    INTEGER          :: mode                        ! 1 = forecast mode, 2 = climate mode
    INTEGER          :: dom(max_phys_dom)           ! domains for which this namelist is used, ending with -1
    INTEGER          :: output_time_unit            ! 1 = second, 2=minute, 3=hour, 4=day, 5=month, 6=year
    INTEGER          :: steps_per_file              ! Max number of output steps in one output file
    LOGICAL          :: include_last                ! Flag whether to include the last timestep in output
    LOGICAL          :: output_grid                 ! Flag whether grid information is output (in NetCDF output)

    ! post-processing times in units defined by output_time_unit: start, end, increment:
    REAL(wp)         :: output_bounds(3,max_bounds) 

    INTEGER          :: taxis_tunit   ! 1 = TUNIT_SECOND, 2 = TUNIT_MINUTE, 3 TUNIT_HOUR ... (see cdi.inc)

    ! --------------------
    ! ready file handling
    ! --------------------

    LOGICAL                     :: lwrite_ready     ! Flag. TRUE if a "ready file" (sentinel file) should be written
    CHARACTER(LEN=filename_max) :: ready_directory  ! output directory for ready files

    ! --------------------
    ! variable lists
    ! --------------------

    CHARACTER(LEN=vname_len)  :: ml_varlist(max_var_ml)   ! name of model level fields
    CHARACTER(LEN=vname_len)  :: pl_varlist(max_var_pl)   ! name of pressure level fields
    CHARACTER(LEN=vname_len)  :: hl_varlist(max_var_hl)   ! name of height level fields
    CHARACTER(LEN=vname_len)  :: il_varlist(max_var_hl)   ! name of isentropic level fields

    ! --------------------
    ! horizontal interpol.
    ! --------------------

    INTEGER  :: remap                 ! interpolate horizontally, 0: none, 1: to regular lat-lon grid, 2: to Gaussian grids, (3:...)
    LOGICAL  :: remap_internal        ! do interpolations online in the model or external (including triggering)
    INTEGER  :: lonlat_id             ! if remap=1: index of lon-lat-grid in global list "lonlat_grid_list"

    ! --------------------
    ! vertical interpol.
    ! --------------------

    REAL(wp) :: p_levels(max_levels)  ! pressure levels [hPa]
    REAL(wp) :: h_levels(max_levels)  ! height levels
    REAL(wp) :: i_levels(max_levels)  ! isentropic levels

    ! -------------------------------------
    ! Internal members, not read from input
    ! -------------------------------------

    INTEGER  :: cur_bounds_triple     ! current output_bounds triple in use
    REAL(wp) :: next_output_time      ! next output time (in seconds simulation time)
    INTEGER  :: n_output_steps
    TYPE(t_output_name_list), POINTER :: next ! Pointer to next output_name_list
  END TYPE t_output_name_list


  ! Unfortunately, Fortran does not allow arrays of pointers, so we
  ! have to define extra types
  TYPE t_rptr_5d
    REAL(wp), POINTER :: p(:,:,:,:,:)
  END TYPE t_rptr_5d


  TYPE t_iptr_5d
    INTEGER,  POINTER :: p(:,:,:,:,:)
  END TYPE t_iptr_5d


  TYPE t_var_desc
    REAL(wp), POINTER :: r_ptr(:,:,:,:,:)         ! Pointer to time level independent REAL data (or NULL)
    INTEGER,  POINTER :: i_ptr(:,:,:,:,:)         ! Pointer to time level independent INTEGER data (or NULL)
    TYPE(t_rptr_5d) :: tlev_rptr(MAX_TIME_LEVELS) ! Pointers to time level dependent REAL data
    TYPE(t_iptr_5d) :: tlev_iptr(MAX_TIME_LEVELS) ! Pointers to time level dependent INTEGER data
    TYPE(t_var_metadata) :: info                  ! Info structure for variable
  END TYPE t_var_desc


  TYPE t_output_file
    ! The following data must be set before opening the output file:
    CHARACTER(LEN=filename_max) :: filename_pref ! Prefix of output file name
    INTEGER                     :: output_type   ! CDI format
    INTEGER                     :: phys_patch_id ! ID of physical output patch
    INTEGER                     :: log_patch_id  ! ID of logical output patch
    REAL(wp)                    :: start_time    ! start time of model domain
    REAL(wp)                    :: end_time      ! end time of model domain
    LOGICAL                     :: initialized   ! .TRUE. if vlist setup has already been called
    INTEGER                     :: ilev_type     ! level type: level_type_ml/level_type_pl/level_type_hl/level_type_il
    INTEGER                     :: max_vars      ! maximum number of variables allocated
    INTEGER                     :: num_vars      ! number of variables in use
    TYPE(t_var_desc),ALLOCATABLE :: var_desc(:)
    TYPE(t_output_name_list), POINTER :: name_list ! Pointer to corresponding output name list

    CHARACTER(LEN=vname_len), ALLOCATABLE :: name_map(:,:) ! mapping internal names -> names in NetCDF

    INTEGER                     :: remap         ! Copy of remap from associated namelist
    INTEGER                     :: io_proc_id    ! ID of process doing I/O on this file

    ! Used for async IO only
    INTEGER(i8)                 :: my_mem_win_off
    INTEGER(i8), ALLOCATABLE    :: mem_win_off(:)

    ! The following members are set during open
    CHARACTER(LEN=filename_max) :: filename           ! Actual name of output file
    CHARACTER(LEN=filename_max) :: rdy_filename       ! Actual name of ready file (if any)
    INTEGER                     :: cdiFileId
    INTEGER                     :: cdiVlistId         ! cdi vlist handler
    INTEGER                     :: cdiCellGridID
    INTEGER                     :: cdiVertGridID
    INTEGER                     :: cdiEdgeGridID
    INTEGER                     :: cdiLonLatGridID
    INTEGER                     :: cdiZaxisID(29) ! All types of possible Zaxis ID's
    INTEGER                     :: cdiTaxisID
    INTEGER                     :: cdiTimeIndex
    INTEGER                     :: cdiInstID      ! output generating institute
    INTEGER                     :: cdi_grb2(3,2)  !< geographical position: (GRID, latitude/longitude)
  END TYPE t_output_file

END MODULE mo_name_list_output_types
