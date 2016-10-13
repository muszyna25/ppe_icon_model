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
MODULE mo_async_latbc_types

  USE mo_kind,                  ONLY: sp
  USE mo_var_metadata_types,    ONLY: t_var_metadata, VARNAME_LEN
  USE mo_dictionary,            ONLY: DICT_MAX_STRLEN

  IMPLICIT NONE

  PRIVATE

  ! module name
  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_async_latbc_types'

  ! derived data types:
  PUBLIC :: t_patch_data
  PUBLIC :: t_reorder_data
  PUBLIC :: t_var_data
  PUBLIC :: t_mem_win
  PUBLIC :: t_v_grid
  PUBLIC :: t_buffer
  PUBLIC :: t_size

  !------------------------------------------------------------------------------------------------
  ! DERIVED DATA TYPES
  !------------------------------------------------------------------------------------------------

  TYPE t_reorder_data
     INTEGER :: n_glb  ! Global number of points per physical patch
     INTEGER :: n_own  ! Number of own points, only belonging to phyiscal patch
     INTEGER, ALLOCATABLE :: reorder_index(:)
     ! Index how to reorder the contributions of all compute PEs
     ! into the global array (set on all PEs)
     ! Only set on compute PEs, set to 0 on prefetching PE
     INTEGER, ALLOCATABLE :: own_idx(:), own_blk(:)
     ! n_own, gathered for all compute PEs (set on all PEs)
     INTEGER, ALLOCATABLE :: pe_own(:)
     ! offset of contributions of PEs (set on all PEs)                  
     INTEGER, ALLOCATABLE :: pe_off(:)
  END TYPE t_reorder_data

  !------------------------------------------------------------------------------------------------

  !> Data structure containing variables for MPI memory window
  !
  TYPE t_mem_win
     ! Currently, we use only 1 MPI window for all input prefetching
     ! Used for async prefetch only
     INTEGER            :: mpi_win
     REAL(sp), POINTER  :: mem_ptr_sp(:) !< Pointer to memory window (REAL*4)
  END TYPE t_mem_win

  !------------------------------------------------------------------------------------------------

  TYPE t_size
     REAL(sp), POINTER :: buffer(:,:,:)
  END TYPE t_size

  !------------------------------------------------------------------------------------------------

  TYPE t_buffer
     INTEGER     :: ngrp_vars    ! Number of variables for prefetching
     CHARACTER(LEN=DICT_MAX_STRLEN),ALLOCATABLE  :: mapped_name(:) ! name of mapped dictionary variables for prefetching
     CHARACTER(LEN=DICT_MAX_STRLEN),ALLOCATABLE  :: internal_name(:) ! corresponding internal name of variables
     CHARACTER(LEN=VARNAME_LEN), ALLOCATABLE :: grp_vars(:) ! name of variables for prefetching
     INTEGER, ALLOCATABLE :: nlev(:) ! Size of variables for prefetching
     TYPE(t_size), ALLOCATABLE :: vars(:) 
     INTEGER, ALLOCATABLE :: varID(:) ! variable ID for variable to be read from file
     INTEGER, ALLOCATABLE :: hgrid(:) ! CDI horizontal grid type ( element stored in either cell
                                      ! or edge location of grid )
     LOGICAL :: lread_qr, lread_qs ! are qr, qs provided as input?
     LOGICAL :: lread_vn ! is vn provided as input?
     LOGICAL :: lthd_progvars ! are prognostic thermodynamic variables (= rho and theta_v) present in the input file?
     CHARACTER(LEN=10) :: psvar
     CHARACTER(LEN=10) :: geop_ml_var ! model level surface geopotential
  END TYPE t_buffer

  TYPE(t_buffer), PUBLIC, TARGET :: latbc_buffer

  !------------------------------------------------------------------------------------------------

  TYPE t_var_data
     TYPE(t_var_metadata) :: info  ! Info structure for variable
  END TYPE t_var_data

  !------------------------------------------------------------------------------------------------
  ! TYPE t_v_grid contains the data of a vertical grid definition.
  !
  TYPE t_v_grid
     INTEGER :: type
     INTEGER :: nlevels
  END type t_v_grid

  !------------------------------------------------------------------------------------------------  
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
     INTEGER :: level
     INTEGER :: num_vars  ! no of input prefetch variables

     ! number of cells and edges in the local patch
     INTEGER :: n_patch_cells
     INTEGER :: n_patch_edges
   
     ! number of cells and edges in the global patch
     INTEGER :: n_patch_cells_g
     INTEGER :: n_patch_edges_g
    
     ! number of points, corresponds to logical patch
     INTEGER :: nblks_c, nblks_e

  END TYPE t_patch_data

END MODULE mo_async_latbc_types

  !------------------------------------------------------------------------------------------------  
