!>
!! This Module is a description for listings and indices
!! used for input prefetching routine.
!!
!! @author M. Pondkule (DWD)
!!
!!
!! @par Revision History
!! Initial release by M. Pondkule, DWD (2013-11-28)
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
MODULE mo_async_latbc_types

  USE mo_kind,                  ONLY: sp
  USE mo_var_metadata_types,    ONLY: t_var_metadata, VARNAME_LEN
  USE mo_dictionary,            ONLY: DICT_MAX_STRLEN
  USE mtime,                    ONLY: event, datetime, timedelta, &
    &                                 deallocateTimedelta, deallocateEvent, deallocateDatetime
  USE mo_initicon_types,        ONLY: t_init_state, t_init_state_const
  USE mo_impl_constants,        ONLY: SUCCESS
  USE mo_exception,             ONLY: finish, message

  IMPLICIT NONE

  PRIVATE

  ! module name
  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_async_latbc_types'

  ! derived data types:
  PUBLIC :: t_latbc_data
  PUBLIC :: t_patch_data
  PUBLIC :: t_reorder_data
  PUBLIC :: t_var_data
  PUBLIC :: t_mem_win
  PUBLIC :: t_buffer
  PUBLIC :: t_size
  PUBLIC :: t_glb_indices


  !------------------------------------------------------------------------------------------------
  ! DERIVED DATA TYPES
  !------------------------------------------------------------------------------------------------


  TYPE t_reorder_data
     INTEGER              :: n_glb ! Global number of points per physical patch
     INTEGER              :: n_own ! Number of own points, only belonging to physical patch
     INTEGER, ALLOCATABLE :: reorder_index(:)
     ! Index how to reorder the contributions of all compute PEs
     ! into the global array (set on all PEs)
     ! Only set on compute PEs, set to 0 on prefetching PE
     INTEGER, ALLOCATABLE :: own_idx(:), own_blk(:)
     ! n_own, gathered for all compute PEs (set on all PEs)
     INTEGER, ALLOCATABLE :: pe_own(:)
     ! offset of contributions of PEs (set on all PEs)
     INTEGER, ALLOCATABLE :: pe_off(:)

     ! logical mask, which local points are read from input file
     LOGICAL, ALLOCATABLE :: read_mask(:,:)

     ! flag: if .TRUE., then the corresponding PE is skipped in the MPI_PUT operation
     LOGICAL              :: this_skip

     ! flag: if .TRUE., then the corresponding PE is skipped in the MPI_PUT operation
     LOGICAL, ALLOCATABLE :: pe_skip(:)

  CONTAINS
    PROCEDURE :: finalize => t_reorder_data_finalize   !< destructor
  END TYPE t_reorder_data


  !> Data structure containing variables for MPI memory window
  !
  TYPE t_mem_win
     ! Currently, we use only 1 MPI window for all input prefetching
     ! Used for async prefetch only
     INTEGER            :: mpi_win
     REAL(sp), POINTER  :: mem_ptr_sp(:) => NULL() !< Pointer to memory window (REAL*4)
  END TYPE t_mem_win


  TYPE t_size
     REAL(sp), POINTER :: buffer(:,:,:) => NULL()
  END TYPE t_size


  TYPE t_buffer
     INTEGER                                     :: ngrp_vars          ! Number of variables for prefetching
     CHARACTER(LEN=DICT_MAX_STRLEN), ALLOCATABLE :: mapped_name(:)     ! name of mapped dictionary variables
     CHARACTER(LEN=DICT_MAX_STRLEN), ALLOCATABLE :: internal_name(:)   ! corresponding internal name of variables
     CHARACTER(LEN=VARNAME_LEN),     ALLOCATABLE :: grp_vars(:)        ! name of variables for prefetching
     INTEGER,                        ALLOCATABLE :: nlev(:)            ! Size of variables for prefetching
     TYPE(t_size),                   ALLOCATABLE :: vars(:)
     INTEGER,                        ALLOCATABLE :: varID(:)           ! ID for variable to be read from file
     INTEGER,                        ALLOCATABLE :: hgrid(:)           ! CDI horizontal grid type (cell/edge)
     LOGICAL                                     :: lread_qr, lread_qs ! are qr, qs provided as input?

     LOGICAL                                     :: lread_vn           ! is vn provided as input?
     LOGICAL                                     :: lread_u_v          ! is u,v provided as input?

     ! If .FALSE., input levels are computed from sfc geopotential:
     LOGICAL                                     :: lread_hhl

     ! are prognostic thermodynamic variables (= rho and theta_v) present in the input file?
     LOGICAL                                     :: lread_theta_rho      

     ! .TRUE., if pressure is read from input
     LOGICAL                                     :: lread_pres

     ! .TRUE., if temperature is read from input
     LOGICAL                                     :: lread_temp

     ! If .TRUE., surface pressure and geopotential are available in
     ! the input file
     LOGICAL                                     :: lread_ps_geop

     ! .FALSE., if vertical component of velocity (W) is provided as input
     LOGICAL                                     :: lconvert_omega2w

     ! .TRUE., if heights are computed (hydrostatic model input):
     LOGICAL                                     :: lcompute_hhl_pres

     CHARACTER(LEN=10)                           :: psvar
     CHARACTER(LEN=10)                           :: geop_ml_var        ! model level surface geopotential
     CHARACTER(LEN=10)                           :: hhl_var

     ! input data validity DateTime
     TYPE(datetime)                              :: vDateTime          

   CONTAINS
     PROCEDURE :: finalize => t_buffer_finalize   !< destructor
  END TYPE t_buffer


  TYPE t_var_data
     TYPE(t_var_metadata) :: info  ! Info structure for variable
  END TYPE t_var_data


  ! TYPE p_patch_info contains the ordering info for cells, edges and verts
  TYPE t_patch_data
     TYPE(t_reorder_data) :: cells
     TYPE(t_reorder_data) :: edges

     TYPE(t_var_data), ALLOCATABLE :: var_data(:)

     ! used for async prefetching only
     TYPE(t_mem_win) :: mem_win  !< data structure containing variables for MPI memory window

     ! number of full and half levels
     INTEGER :: nlev
     INTEGER :: nlevp1

     INTEGER :: level     ! patch level (e.g. xx in R03Bxx)
     INTEGER :: num_vars  ! no of input prefetch variables

     ! number of cells and edges in the local patch
     INTEGER :: n_patch_cells
     INTEGER :: n_patch_edges

     ! number of cells and edges in the global patch
     INTEGER :: n_patch_cells_g
     INTEGER :: n_patch_edges_g

     ! number of blocks
     INTEGER :: nblks_c, nblks_e

   CONTAINS
     PROCEDURE :: finalize => t_patch_data_finalize   !< destructor
  END TYPE t_patch_data


  !> Sparse latbc mode: index data for boundary rows:
  !  Derived type specifying a local-to-global index mapping for
  !  extracted subgrids.
  !
  TYPE t_glb_indices
    INTEGER, ALLOCATABLE :: cells(:), edges(:)      !< (1...local) global indices for cells and edges
    INTEGER              :: n_patch_cells_g         !< total no. of global cells
    INTEGER              :: n_patch_edges_g         !< total no. of global edges
  CONTAINS
    PROCEDURE :: finalize => t_glb_indices_finalize
  END TYPE t_glb_indices


  !> Data type containing the time control, the variable buffer, and
  !> the necessary index arrays for reordering, to read in lateral
  !> boundary data into ICON.
  !
  TYPE t_latbc_data

    TYPE(datetime),  POINTER :: mtime_last_read => NULL()
    TYPE(event),     POINTER :: prefetchEvent   => NULL()
    TYPE(timedelta), POINTER :: delta_dtime     => NULL()

    ! time level indices for  latbc_data. can be 1 or 2.
    INTEGER :: new_latbc_tlev

    ! storage for time-constant height level data
    TYPE(t_init_state_const) :: latbc_data_const

    ! storage for two time-level boundary data
    TYPE(t_init_state) :: latbc_data(2)

    ! raw buffer
    TYPE(t_buffer) :: buffer

    !< indices for async latbc prefetching
    TYPE(t_patch_data) :: patch_data

    ! for sparse latbc mode: index data for boundary rows:
    TYPE(t_glb_indices) :: global_index
  CONTAINS
    PROCEDURE :: finalize => t_latbc_data_finalize
    PROCEDURE :: prev_latbc_tlev => prev_latbc_tlev

  END TYPE t_latbc_data


CONTAINS

  !> Destructor for "t_glb_indices" class.
  !
  SUBROUTINE t_glb_indices_finalize(this)
    CLASS(t_glb_indices) :: this
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = modname//'::t_glb_indices_finalize'
    INTEGER :: ierrstat=0

    !CALL message(routine, 't_glb_indices_finalize')

    IF (ALLOCATED(this%cells))  DEALLOCATE(this%cells, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "DEALLOCATE failed!")
    IF (ALLOCATED(this%edges))  DEALLOCATE(this%edges, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "DEALLOCATE failed!")
  END SUBROUTINE t_glb_indices_finalize


  SUBROUTINE t_patch_data_finalize(patch_data)
    CLASS(t_patch_data), INTENT(INOUT) :: patch_data

    !CALL message("", 't_patch_data_finalize')

    IF (ALLOCATED(patch_data%var_data))             DEALLOCATE(patch_data%var_data)
    CALL patch_data%cells%finalize()
    CALL patch_data%edges%finalize()
    ! note: we do not touch the MPI window pointer here:
    !
    ! IF (ASSOCIATED(patch_data%mem_win%mem_ptr_sp))  DEALLOCATE(patch_data%mem_win%mem_ptr_sp)
  END SUBROUTINE t_patch_data_finalize


  SUBROUTINE t_reorder_data_finalize(reorder_data)
    CLASS(t_reorder_data), INTENT(INOUT) :: reorder_data

    !CALL message("", 't_reorder_data_finalize')

    IF (ALLOCATED(reorder_data%reorder_index)) DEALLOCATE(reorder_data%reorder_index)
    IF (ALLOCATED(reorder_data%own_idx))       DEALLOCATE(reorder_data%own_idx)
    IF (ALLOCATED(reorder_data%own_blk))       DEALLOCATE(reorder_data%own_blk)
    IF (ALLOCATED(reorder_data%pe_own))        DEALLOCATE(reorder_data%pe_own)
    IF (ALLOCATED(reorder_data%pe_off))        DEALLOCATE(reorder_data%pe_off)
    IF (ALLOCATED(reorder_data%read_mask))     DEALLOCATE(reorder_data%read_mask)
    IF (ALLOCATED(reorder_data%pe_skip))       DEALLOCATE(reorder_data%pe_skip)
  END SUBROUTINE t_reorder_data_finalize


  SUBROUTINE t_buffer_finalize(buffer)
    CLASS(t_buffer), INTENT(INOUT) :: buffer
    INTEGER :: i

    !CALL message("", 't_buffer_finalize')

    IF (ALLOCATED(buffer%mapped_name))    DEALLOCATE(buffer%mapped_name)
    IF (ALLOCATED(buffer%internal_name))  DEALLOCATE(buffer%internal_name)
    IF (ALLOCATED(buffer%grp_vars))       DEALLOCATE(buffer%grp_vars)
    IF (ALLOCATED(buffer%nlev))           DEALLOCATE(buffer%nlev)
    IF (ALLOCATED(buffer%vars)) THEN
      DO i=1,SIZE(buffer%vars)
        IF (ASSOCIATED(buffer%vars(i)%buffer))  DEALLOCATE(buffer%vars(i)%buffer)      
      END DO
      DEALLOCATE(buffer%vars)
    END IF
    IF (ALLOCATED(buffer%varID))          DEALLOCATE(buffer%varID)
    IF (ALLOCATED(buffer%hgrid))          DEALLOCATE(buffer%hgrid)
  END SUBROUTINE t_buffer_finalize


  !-------------------------------------------------------------------------
  !>
  !! @par Revision History
  !! Initial version by S. Brdar, DWD (2013-06-13)
  !! Modified version by M. Pondkule, DWD (2013-04-17)
  !!
  SUBROUTINE t_latbc_data_finalize(latbc)
    CLASS(t_latbc_data), INTENT(INOUT) :: latbc

    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::t_latbc_data_finalize"
    INTEGER :: tlev

    CALL message("", 'deallocating latbc data')

    ! deallocate boundary data memory
    DO tlev = 1, 2
      CALL latbc%latbc_data(tlev)%finalize()
    END DO
    
    CALL latbc%patch_data%finalize()          ! deallocate patch data
    CALL latbc%buffer%finalize()              ! deallocate intermediate storage latbc%buffer
    CALL latbc%global_index%finalize()        ! clean up global indices data structure.

    IF (ASSOCIATED(latbc%prefetchEvent)) THEN
      CALL deallocateEvent(latbc%prefetchEvent) ! deallocate prefetch input event
    END IF

    ! deallocating date and time data structures
    IF (ASSOCIATED(latbc%mtime_last_read)) THEN
      CALL deallocateDatetime(latbc%mtime_last_read)
    END IF
    IF (ASSOCIATED(latbc%delta_dtime)) THEN
      CALL deallocateTimedelta(latbc%delta_dtime)
    END IF
  END SUBROUTINE t_latbc_data_finalize


  INTEGER FUNCTION prev_latbc_tlev(latbc)
    CLASS(t_latbc_data), INTENT(IN) :: latbc

    prev_latbc_tlev = 3 - latbc%new_latbc_tlev
  END FUNCTION prev_latbc_tlev

END MODULE mo_async_latbc_types

!------------------------------------------------------------------------------------------------
