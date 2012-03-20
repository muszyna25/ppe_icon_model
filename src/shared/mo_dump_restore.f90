!>
!! This module provides routines for dumping/restoring patches
!! (including interpolation and grf state) to/from NetCDF.
!!
!! @par Revision History
!! Initial version by Rainer Johanni, Sep 2010
!!
!! @par Copyright
!! 2002-2007 by DWD and MPI-M
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
!! $Id: n/a$
!!
!-------------------------------------------------------------------------------
!
! Conventions for NetCDF output
!
! - all variables which run over nproma/nblks_[cev] are output independent of nproma
!   with the 1st dimension having the length n_patch_cells/edges/verts.
!   The remaining dimensions (if any) follow in their original order.
!
!   For example:
!
!   area(nproma,nblks_c)         is output as   area(n_patch_cells)
!   c_lin_e(nproma,2,nblks_e)    is output as   c_lin_e(n_patch_edges,2)
!   tria_north(3,nproma,nblks_e) is output as   tria_north(n_patch_edges,3)
!   lsq_qtmat_c(nproma,lsq_dim_unk,lsq_dim_c,nblks_c)
!        is output as   lsq_qtmat_c(n_patch_cells,lsq_dim_unk,lsq_dim_c)
!
! - Index pairs in the form xxx_idx/xxx_blk are joined to one linear index,
!   the variable name in the NetCDF file is xxx_index.
!
!   For example:
!
!   neighbor_idx(nproma,nblks_c) AND neighbor_blk(nproma,nblks_c)
!   is output as neighbor_index(n_patch_cells)
!
!   heli_vn_idx(nincr,nproma,nblks_e) AND heli_vn_blk(nincr,nproma,nblks_e)
!   is output as heli_vn_index(n_patch_edges,nincr)
!
! - Variables not depending on nproma/nblks are output in the shape they have,
!   idx/blk pairs are also joined and output with the suffix "index"
!
! - The variable names in the NetCDF file are chosen in a way that they
!   resemble closely the original variable name.
!
!   For example:
!   patch%cells%edge_orientation has the NetCDF name "patch.cells.edge_orientation"
!
!   Please note that a '%' is not legal in a NetCDF variable name so the
!   dot has been chosen as separator (like in the C programming language).
!
!-------------------------------------------------------------------------------
!
! How to add a new variable to NetCDF I/O
!
! - Check if all dimensions of the variable (according to the above conventions)
!   are already defined in dump_patch_state_netcdf(), if necessary add a new
!   dimension
!
! - define the variable using subroutine def_var in one of the routines
!   def_patch/def_int_state/def_grf_state
!   def_var has the following arguments:
!   - Name of the variable (freely chosen according to above conventions)
!   - NetCDF type: nf_int of nf_double, logical variables are mapped to integers
!   - up to 3 predefined dimensions
!
! - Add an I/O statement for the variable in one of the routines
!   patch_io/int_state_io/grf_state_io
!
!   For variables depending on nproma/nblocks
!     which are NOT index pairs use  bvar_io
!     which are index pairs use      bidx_io
!   For variables NOT depending on nproma/nblocks
!     which are NOT index pairs use  uvar_io
!     which are index pairs use      uidx_io
!
!   The first 2 arguments of bvar_io/bidx_io are the positions of nproma/nblocks
!   in the array dimensions (called pos_nproma/pos_nblks), e.g:
!
!   area(nproma,nblks_c)         has pos_nproma=1, pos_nblks=2
!   c_lin_e(nproma,2,nblks_e)    has pos_nproma=1, pos_nblks=3
!   tria_north(3,nproma,nblks_e) has pos_nproma=2, pos_nblks=3
!
!   Then there follows the variable name (argument 3) and the array (bvar_io)
!   or the index pair (bidx_io).
!
!   uvar_io/uidx_io take just the variable name, the array (uvar_io)
!   or the index pair (uidx_io).
!
!-------------------------------------------------------------------------------
!
! Note about the unlimited dimension:
!
! Almost every variable has the NetCDF unlimited dimension as the last dimension
! which is set to
! - the number of work PEs if l_one_file_per_patch is in effect
! - 1 if every PE writes a file of its own
!
! Although it seems to be unnecessary to use the unlimited dimension here
! since the number of work PEs a priori (and could be used as a dimension),
! it is important to do it this way because of the internal structure
! of a NetCDF file:
! Variables having only fixed dimensions are stored at a contiguous location
! which means that data belonging to one PE would be spread in pieces
! over the whole file when the number of work PEs would be a fixed dimension.
!
! Using the unlimited dimension guarantees that the data for a single PE is
! contiguous and thus can be read efficiently.
!
! The only exception to this rule are variables which are identical for every PE:
! They are only stored once (without the unlimited dimension) for l_one_file_per_patch
!
! If l_one_file_per_patch is not in effect, the unlimited dimension is not
! necessary but it is left there in order not to complicate the code.
!
!-------------------------------------------------------------------------------


MODULE mo_dump_restore

  USE mo_kind,               ONLY: wp
  USE mo_impl_constants,     ONLY: min_rlcell, max_rlcell,  &
                                   min_rledge, max_rledge, &
                                   min_rlvert, max_rlvert, &
                                   min_rlcell_int, min_rledge_int
  USE mo_exception,          ONLY: message_text, message, finish, warning
  USE mo_parallel_config, ONLY: nproma
  USE mo_run_config,         ONLY: l_one_file_per_patch, ltransport, &
     &                             num_lev, num_levp1, nshift,       &
     &                             dump_filename, dd_filename
  USE mo_dynamics_config,    ONLY: iequations
  USE mo_io_units,           ONLY: filename_max, nerr
  USE mo_model_domain,       ONLY: t_patch, p_patch_local_parent
  USE mo_grid_config,        ONLY: start_lev, n_dom, n_dom_start, lfeedback, &
                                   l_limited_area, max_childdom, dynamics_parent_grid_id, &
                                   global_cell_type
  USE mo_intp_data_strc      ! We need all from that module
  USE mo_grf_intp_data_strc  ! We need all from that module
!  USE mo_interpol_nml        ! We need all from that module
 USE mo_interpol_config      ! We need all from that module
  USE mo_gridref_config      ! We need all from that module
  USE mo_mpi,                ONLY: my_process_is_mpi_all_parallel, p_n_work, p_pe_work, &
    &                              process_mpi_io_size, my_process_is_stdio, &
    &                              my_process_is_mpi_workroot, p_int, p_comm_work, &
    &                              num_work_procs, p_barrier, get_my_mpi_work_id, p_max, p_pe
  USE mo_impl_constants_grf, ONLY: grf_bdyintp_start_c, grf_bdyintp_start_e
  USE mo_communication,      ONLY: t_comm_pattern, blk_no, idx_no, idx_1d
  USE mo_model_domimp_patches, ONLY: allocate_patch
  USE mo_intp_state,         ONLY: allocate_int_state, &
    &                              allocate_int_state_lonlat
  USE mo_grf_intp_state,     ONLY: allocate_grf_state

  USE mo_model_domimp_patches, ONLY: set_patches_grid_filename
  USE mo_lonlat_intp_config, ONLY: lonlat_intp_config
  USE mo_math_utilities,     ONLY: t_lon_lat_grid
  USE mo_util_string,        ONLY: t_keyword_list, MAX_STRING_LEN,   &
    &                              associate_keyword, with_keywords, &
    &                              int2string
  USE mo_master_nml,         ONLY: model_base_dir

  IMPLICIT NONE

  PRIVATE

  CHARACTER(LEN=*), PARAMETER :: version = '$Id$'
  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_dump_restore'

  !modules interface-------------------------------------------
  !subroutines
  PUBLIC :: dump_patch_state_netcdf
  PUBLIC :: dump_domain_decomposition
  PUBLIC :: dump_all_domain_decompositions
  PUBLIC :: restore_patches_netcdf
  PUBLIC :: restore_interpol_state_netcdf
  PUBLIC :: restore_gridref_state_netcdf

  INCLUDE 'netcdf.inc'

  INTERFACE bvar_io
    MODULE PROCEDURE bvar_io_r1
    MODULE PROCEDURE bvar_io_r2
    MODULE PROCEDURE bvar_io_r3
    MODULE PROCEDURE bvar_io_i1
    MODULE PROCEDURE bvar_io_i2
    MODULE PROCEDURE bvar_io_l1
  END INTERFACE
  INTERFACE bidx_io
    MODULE PROCEDURE bidx_io_i1
    MODULE PROCEDURE bidx_io_i2
  END INTERFACE
  INTERFACE uvar_io
    MODULE PROCEDURE uvar_io_r1
    MODULE PROCEDURE uvar_io_r2
    MODULE PROCEDURE uvar_io_i1
    MODULE PROCEDURE uvar_io_i2
  END INTERFACE
  INTERFACE uidx_io
    MODULE PROCEDURE uidx_io_i1
    MODULE PROCEDURE uidx_io_i2
  END INTERFACE
  INTERFACE check_att
    MODULE PROCEDURE check_att_int
    MODULE PROCEDURE check_att_double
    MODULE PROCEDURE check_att_int_array
    MODULE PROCEDURE check_att_double_array
  END INTERFACE

  INTEGER, PARAMETER :: grf_vec_dim_1 = 6
  INTEGER, PARAMETER :: grf_vec_dim_2 = 5

  INTEGER :: ncid
  LOGICAL :: netcdf_read = .FALSE.

  CHARACTER(LEN=filename_max) :: filename

  ! prefix for variable names (used if local parent is dumped)
  CHARACTER(LEN=64) :: prefix

  ! flag if current processor should output NetCDF defines
  LOGICAL :: output_defines

  ! number of my record within unlimited dimension
  INTEGER :: my_record

  ! Flag whether to use 1 file per patch, this is always set for dd output

  LOGICAL :: use_one_file

  ! The following dimensions must be preserved between calls
  ! thus must be global

  INTEGER :: dim_unlimited ! The unlimited dimension

  INTEGER :: dim_2, dim_3, dim_4, dim_5, dim_7, dim_6, dim_8, dim_9
  INTEGER :: dim_ncells, dim_nedges, dim_nverts
  INTEGER :: dim_nverts_per_cell, dim_nverts_per_cell_p1
  INTEGER :: dim_nedges_per_vert

  INTEGER :: dim_nrlcell, dim_nrledge, dim_nrlvert

  INTEGER :: dim_nincr
  INTEGER :: dim_lsq_dim_c_lin, dim_lsq_dim_unk_lin, dim_lsq_dim_unk2_lin
  INTEGER :: dim_lsq_dim_c_high, dim_lsq_dim_unk_high, dim_lsq_dim_unk2_high
  INTEGER :: dim_rbf_c2grad_dim
  INTEGER :: dim_rbf_vec_dim_c, dim_rbf_vec_dim_e, dim_rbf_vec_dim_v

  INTEGER :: dim_grf_vec_dim_1, dim_grf_vec_dim_2

  ! NOTE: These are either the max over all patches or same as n_patch_cells etc.
  INTEGER :: max_patch_cells, max_patch_edges, max_patch_verts

  ! Dimension IDs which pertain to global max values ...
  INTEGER :: dim_max_cells, dim_max_verts, dim_max_edges
  INTEGER :: dim_max_grf_cells, dim_max_grf_edges ! pertain to actual child

  ! ... and their actual values on the current proc
  INTEGER :: n_patch_cells, n_patch_edges, n_patch_verts
  INTEGER :: n_grf_cells, n_grf_edges ! pertain to actual child

  !-------------------------------------------------------------------------

CONTAINS
  !
  !-----------------------------------------------------------------------
  ! Help functions for NetCDF I/O
  !-----------------------------------------------------------------------
  !
  !> Checks the return value of a NetCDF function and exits in case of
  !! an error

  SUBROUTINE nf(status)

    INTEGER, INTENT(in) :: status

    IF (status /= nf_noerr) THEN
      CALL finish(modname, 'NetCDF Error: '//nf_strerror(status))
    ENDIF

  END SUBROUTINE nf

  !-----------------------------------------------------------------------
  !
  !> set_dump_restore_filename:
  !! Sets the filename for NetCDF dump/restore files
  !! from the patch_filename
  SUBROUTINE set_dump_restore_filename(patch_filename, model_base_dir)

    CHARACTER(LEN=*), INTENT(in)   :: patch_filename
    CHARACTER(len=*), INTENT(IN)   :: model_base_dir

    TYPE (t_keyword_list), POINTER :: keywords => NULL()
    CHARACTER(len=MAX_STRING_LEN)  :: proc_str

    IF(use_one_file) THEN
      proc_str = ""
    ELSE
      WRITE (proc_str,'(a,i0,a,i0,a)') "proc", p_pe_work, "of", p_n_work, "_"
    END IF

    CALL associate_keyword("<path>",     TRIM(model_base_dir),  keywords)
    CALL associate_keyword("<proc>",     TRIM(proc_str),        keywords)
    CALL associate_keyword("<gridfile>", TRIM(patch_filename),  keywords)
    filename = TRIM(with_keywords(keywords, TRIM(dump_filename)))

  END SUBROUTINE set_dump_restore_filename

  !-----------------------------------------------------------------------
  !
  !> set_dd_filename:
  !! Sets the filename for NetCDF domain decomposition files
  !! from the patch_filename
  SUBROUTINE set_dd_filename(patch_filename, model_base_dir)

    CHARACTER(LEN=*), INTENT(in)   :: patch_filename
    CHARACTER(len=*), INTENT(IN)   :: model_base_dir
    TYPE (t_keyword_list), POINTER :: keywords => NULL()

    ! Please note: Currently use_one_file is enforced in this case
    CALL associate_keyword("<path>",     TRIM(model_base_dir),  keywords)
    CALL associate_keyword("<gridfile>", TRIM(patch_filename),  keywords)
    filename = TRIM(with_keywords(keywords, TRIM(dd_filename)))

  END SUBROUTINE set_dd_filename

  !-----------------------------------------------------------------------
  !
  !> def_var:
  !! Defines a NetCDF variable with up to 3 dimensions
  !! Note that the unlimited dimension is added by default!

  SUBROUTINE def_var(var_name, var_type, dim1, dim2, dim3, add_unlim)

    CHARACTER(LEN=*), INTENT(IN) :: var_name
    INTEGER, INTENT(IN) :: var_type
    INTEGER, INTENT(IN), OPTIONAL :: dim1, dim2, dim3
    LOGICAL, INTENT(IN), OPTIONAL :: add_unlim

    INTEGER nvdims, vdims(4), varid
    LOGICAL l_add_unlim

    IF(PRESENT(add_unlim)) THEN
      l_add_unlim = add_unlim
    ELSE
      l_add_unlim = .TRUE.
    ENDIF

    nvdims = 0

    IF(PRESENT(dim1)) THEN
      nvdims = 1
      vdims(1) = dim1
      IF(PRESENT(dim2)) THEN
        nvdims = 2
        vdims(2) = dim2
        IF(PRESENT(dim3)) THEN
          nvdims = 3
          vdims(3) = dim3
        ENDIF
      ENDIF
    ENDIF

    IF(l_add_unlim) THEN
      nvdims = nvdims+1
      vdims(nvdims) = dim_unlimited
    ENDIF

    CALL nf(nf_def_var(ncid, TRIM(prefix)//var_name, var_type, nvdims, vdims, varid))

  END SUBROUTINE def_var

  !-----------------------------------------------------------------------
  !
  !> def_dim:
  !! Defines a NetCDF dimension (prefixing the name)

  SUBROUTINE def_dim(dim_name, dim_size, dim_id)

    CHARACTER(LEN=*) :: dim_name
    INTEGER, INTENT(IN) :: dim_size
    INTEGER, INTENT(OUT) :: dim_id

    CALL nf(nf_def_dim(ncid, TRIM(prefix)//dim_name, dim_size, dim_id))

  END SUBROUTINE def_dim

  !-----------------------------------------------------------------------
  !
  !> store_proc_dependent_dimensions
  !! Stores the dimension IDs of the max_... dimensions and their value on the current processor
  !! so that these dimensions can be replaced by the actual ones during IO

  SUBROUTINE store_proc_dependent_dimensions(p)

    TYPE(t_patch), INTENT(IN) :: p

    ! Get the dimension IDs for the max_... dimensions and set the actual values
    ! on the current processor for these dimensions
    ! This is only necessary if use_one_file is set, otherwise
    ! we will get an error for empty patches on the current proc

    IF(use_one_file) THEN
      CALL nf(nf_inq_dimid(ncid, TRIM(prefix)//'max_patch_cells', dim_max_cells))
      CALL nf(nf_inq_dimid(ncid, TRIM(prefix)//'max_patch_edges', dim_max_edges))
      CALL nf(nf_inq_dimid(ncid, TRIM(prefix)//'max_patch_verts', dim_max_verts))
    ELSE
      dim_max_cells = -1
      dim_max_edges = -1
      dim_max_verts = -1
    ENDIF

    n_patch_cells = p%n_patch_cells
    n_patch_edges = p%n_patch_edges
    n_patch_verts = p%n_patch_verts

    ! Reset the values which will be set only in grf_state_io
    ! but MUST have defined values during IO

    dim_max_grf_cells = -1
    dim_max_grf_edges = -1
    n_grf_cells = -1 ! for safety
    n_grf_edges = -1 ! for safety

  END SUBROUTINE store_proc_dependent_dimensions

  !-----------------------------------------------------------------------
  !
  !> get_var_info:
  !! Checks for a variable with name var_name if the type is the same
  !! as set in var_type, if the number of dimensions is as set in ndims
  !! and retreives the NetCDF variable ID (varid) and the length of
  !! every dimension (dimlen)

  SUBROUTINE get_var_info(var_name, var_type, ndims, varid, dimlen)

    CHARACTER(LEN=*), INTENT(IN) :: var_name
    INTEGER, INTENT(IN)  :: var_type, ndims
    INTEGER, INTENT(OUT) :: varid, dimlen(:)

    INTEGER i, stat, dimids(ndims+1), var_type_inq, ndims_inq, unlimdimid

    ! Get varid, check return status explicitly so we can output the variable name
    ! in case of an error (most probably a typo in var_name)

    stat = nf_inq_varid(ncid, TRIM(prefix)//var_name, varid)
    IF (stat /= nf_noerr) &
      CALL finish(modname, TRIM(prefix)//var_name//': '//nf_strerror(stat))

    ! Check if var_type is ok

    CALL nf(nf_inq_vartype(ncid, varid, var_type_inq))
    IF(var_type_inq /= var_type) &
      CALL finish(modname, TRIM(prefix)//var_name//': stored type not like expected')

    ! Check if ndims is ok

    CALL nf(nf_inq_varndims(ncid, varid, ndims_inq))
    if(ndims_inq /= ndims .AND. ndims_inq /= ndims+1) &
      CALL finish(modname, TRIM(prefix)//var_name//': stored ndims not like expected')

    ! Get dimensions

    CALL nf(nf_inq_vardimid(ncid, varid, dimids))

    DO i = 1, ndims ! Don't return the unlimited dimension
      CALL nf(nf_inq_dimlen(ncid, dimids(i), dimlen(i)))
    ENDDO

    ! Check if last dimension is the unlimited one

    IF(ndims_inq == ndims+1) THEN
      CALL nf(nf_inq_unlimdim(ncid, unlimdimid))
      IF(dimids(ndims+1) /= unlimdimid) &
        CALL finish(modname, TRIM(prefix)//var_name//': last dim not NF_UNLIMITED')
    ENDIF

    ! If the first dimension is one of the processor dependent dimensions, replace it by actual one

    IF(dimids(1) == dim_max_cells) dimlen(1) = n_patch_cells
    IF(dimids(1) == dim_max_edges) dimlen(1) = n_patch_edges
    IF(dimids(1) == dim_max_verts) dimlen(1) = n_patch_verts

    IF(dimids(1) == dim_max_grf_cells) dimlen(1) = n_grf_cells
    IF(dimids(1) == dim_max_grf_edges) dimlen(1) = n_grf_edges

  END SUBROUTINE get_var_info

  !-----------------------------------------------------------------------
  !
  !> pos_check:
  !! Checks if pos_nproma, pos_nblks is in the allowed range,
  !! checks if the stored dimensions for the variable name (in dimlen)
  !! correspond to the actual dimensions of the data array (in var_ubound),
  !! returns the position of the (optional) 2nd and 3rd dimension
  !! within the the dimensions of the data array.

  SUBROUTINE pos_check(pos_nproma, pos_nblks, var_name, var_ubound, var_start, &
    &                  dimlen, pos_dim2, pos_dim3)

    INTEGER, INTENT(IN) :: pos_nproma, pos_nblks, var_ubound(:), var_start, dimlen(:)
    CHARACTER(LEN=*), INTENT(IN) :: var_name
    INTEGER, INTENT(OUT), OPTIONAL :: pos_dim2, pos_dim3

    ! Please note that pos_nproma must always be smaller than pos_nblks!

    IF(.NOT. PRESENT(pos_dim2)) THEN
      ! Variable has 2 dimensions
      IF(pos_nproma /= 1 .OR. pos_nblks /= 2) &
        CALL finish(modname, TRIM(prefix)//var_name//': pos_nproma/pos_nblks illegal')
    ELSE
      IF(.NOT. PRESENT(pos_dim3)) THEN
        ! Variable has 3 dimensions
        IF(pos_nproma == 1 .AND. pos_nblks == 2)  THEN
          pos_dim2 = 3
        ELSEIF(pos_nproma == 1 .AND. pos_nblks == 3)  THEN
          pos_dim2 = 2
        ELSEIF(pos_nproma == 2 .AND. pos_nblks == 3)  THEN
          pos_dim2 = 1
        ELSE
          CALL finish(modname, TRIM(prefix)//var_name//': pos_nproma/pos_nblks illegal')
        ENDIF
        ! Dimension at pos_dim2 must correspond to allocated dimension in NetCDF
        IF(var_ubound(pos_dim2) /= dimlen(2)) &
          CALL finish(modname, TRIM(prefix)//var_name//': Dimension mismatch for dim2')
      ELSE
        ! Variable has 4 dimensions
        IF(pos_nproma == 1 .AND. pos_nblks == 2)  THEN
          pos_dim2 = 3; pos_dim3 = 4
        ELSEIF(pos_nproma == 1 .AND. pos_nblks == 3)  THEN
          pos_dim2 = 2; pos_dim3 = 4
        ELSEIF(pos_nproma == 1 .AND. pos_nblks == 4)  THEN
          pos_dim2 = 2; pos_dim3 = 3
        ELSEIF(pos_nproma == 3 .AND. pos_nblks == 4)  THEN
          pos_dim2 = 1; pos_dim3 = 2
        ELSE
          ! Other pos_nproma/pos_nblks combinations are possible but currently not existing!
          CALL finish(modname, TRIM(prefix)//var_name//': pos_nproma/pos_nblks illegal')
        ENDIF
        ! Dimension at pos_dim2/pos_dim3 must correspond to allocated dimension in NetCDF
        IF(var_ubound(pos_dim2) /= dimlen(2)) THEN
          CALL finish(modname, TRIM(prefix)//var_name//': Dimension mismatch for dim2')
        END IF
        IF(var_ubound(pos_dim3) /= dimlen(3)) &
          CALL finish(modname, TRIM(prefix)//var_name//': Dimension mismatch for dim3')
      ENDIF
    ENDIF

    ! Dimension at pos_nproma must be nproma

    IF(var_ubound(pos_nproma) /= nproma) &
      CALL finish(modname, TRIM(prefix)//var_name//': Dimension at pos_nproma IS NOT nproma')

    ! Dimension at pos_nblks must match dimlen(1)

    IF(var_ubound(pos_nblks) /= (dimlen(1)+var_start-2)/nproma+1) &
      CALL finish(modname, TRIM(prefix)//var_name//': Dimension at pos_nblks not like expected')

  END SUBROUTINE pos_check
  !
  !-----------------------------------------------------------------------
  !> bvar_io family:
  !! I/O of blocked variables, i.e. depending on nproma/nblks
  !! pos_nproma:    position of nproma within the indices
  !! pos_nblks:     position of nblks within the indices
  !! var_name:      name of variable
  !! var:           data array to be read/written
  !! opt_start:     optional index of the first value for I/O in var
  !-----------------------------------------------------------------------
  ! Please note:
  ! These routines work for variables with or without the unlimited dimension.
  ! If no unlimited dimension is present, the start/count values belonging
  ! to the unlimited dimension will just be ignored.
  !-----------------------------------------------------------------------
  !
  SUBROUTINE bvar_io_r1(pos_nproma, pos_nblks, var_name, var, opt_start)

    INTEGER, INTENT(IN) :: pos_nproma, pos_nblks
    CHARACTER(LEN=*), INTENT(IN) :: var_name
    REAL(wp), INTENT(INOUT) :: var(:,:)
    INTEGER, INTENT(IN), OPTIONAL :: opt_start

    INTEGER :: var_start, varid, dimlen(1)
    INTEGER :: nblks, start(2), count(2)
    REAL(wp), ALLOCATABLE :: buf(:)

    IF(PRESENT(opt_start)) THEN
      var_start = opt_start
    ELSE
      var_start = 1
    ENDIF

    CALL get_var_info(var_name, nf_double, 1, varid, dimlen)

    ! Do safety checks
    CALL pos_check(pos_nproma, pos_nblks, var_name, UBOUND(var), var_start, dimlen)

    nblks = UBOUND(var,pos_nblks)
    ALLOCATE(buf(nproma*nblks))
    buf(:) = 0._wp

    start = (/ 1, my_record /)
    count = (/ dimlen(1), 1 /)
    IF(netcdf_read) THEN
      CALL nf(nf_get_vara_double(ncid, varid, start, count, buf(var_start)))
      var(:,:) = RESHAPE(buf, (/ nproma, nblks /))
    ELSE
      buf(:) = RESHAPE(var(:,:), (/ nproma*nblks /))
      CALL nf(nf_put_vara_double(ncid, varid, start, count, buf(var_start)))
    ENDIF

    DEALLOCATE(buf)

  END SUBROUTINE bvar_io_r1
  !-----------------------------------------------------------------------
  SUBROUTINE bvar_io_r2(pos_nproma, pos_nblks, var_name, var, opt_start)

    INTEGER, INTENT(IN) :: pos_nproma, pos_nblks
    CHARACTER(LEN=*), INTENT(IN) :: var_name
    REAL(wp), INTENT(INOUT) :: var(:,:,:)
    INTEGER, INTENT(IN), OPTIONAL :: opt_start

    INTEGER :: var_start, varid, dimlen(2)
    INTEGER :: start(3), count(3), pos_dim2
    INTEGER :: nblks, i
    REAL(wp), ALLOCATABLE :: buf(:)

    IF(PRESENT(opt_start)) THEN
      var_start = opt_start
    ELSE
      var_start = 1
    ENDIF

    CALL get_var_info(var_name, nf_double, 2, varid, dimlen)

    ! Do safety checks, get position of remaining dimensions
    CALL pos_check(pos_nproma, pos_nblks, var_name, UBOUND(var), var_start, dimlen, pos_dim2)

    nblks = UBOUND(var,pos_nblks)
    ALLOCATE(buf(nproma*nblks))
    buf(:) = 0._wp

    DO i=1,UBOUND(var,pos_dim2)
      start = (/ 1, i, my_record /)
      count = (/ dimlen(1), 1, 1 /)
      IF(netcdf_read) THEN
        CALL nf(nf_get_vara_double(ncid, varid, start, count, buf(var_start)))
        IF(pos_dim2==1) var(i,:,:) = RESHAPE(buf, (/ nproma, nblks /))
        IF(pos_dim2==2) var(:,i,:) = RESHAPE(buf, (/ nproma, nblks /))
        IF(pos_dim2==3) var(:,:,i) = RESHAPE(buf, (/ nproma, nblks /))
      ELSE
        IF(pos_dim2==1) buf(:) = RESHAPE(var(i,:,:), (/ nproma*nblks /))
        IF(pos_dim2==2) buf(:) = RESHAPE(var(:,i,:), (/ nproma*nblks /))
        IF(pos_dim2==3) buf(:) = RESHAPE(var(:,:,i), (/ nproma*nblks /))
        CALL nf(nf_put_vara_double(ncid, varid, start, count, buf(var_start)))
      ENDIF
    ENDDO

    DEALLOCATE(buf)

  END SUBROUTINE bvar_io_r2
  !-----------------------------------------------------------------------
  SUBROUTINE bvar_io_r3(pos_nproma, pos_nblks, var_name, var, opt_start)

    INTEGER, INTENT(IN) :: pos_nproma, pos_nblks
    CHARACTER(LEN=*), INTENT(IN) :: var_name
    REAL(wp), INTENT(INOUT) :: var(:,:,:,:)
    INTEGER, INTENT(IN), OPTIONAL :: opt_start

    INTEGER :: var_start, varid, dimlen(3)
    INTEGER :: start(4), count(4), pos_dim2, pos_dim3
    INTEGER :: nblks, i, j
    REAL(wp), ALLOCATABLE :: buf(:)

    IF(PRESENT(opt_start)) THEN
      var_start = opt_start
    ELSE
      var_start = 1
    ENDIF

    CALL get_var_info(var_name, nf_double, 3, varid, dimlen)

    ! Do safety checks, get position of remaining dimensions
    CALL pos_check(pos_nproma, pos_nblks, var_name, UBOUND(var), var_start, &
      &            dimlen, pos_dim2, pos_dim3)

    nblks = UBOUND(var,pos_nblks)
    ALLOCATE(buf(nproma*nblks))
    buf(:) = 0._wp

    DO j=1,UBOUND(var,pos_dim3)
    DO i=1,UBOUND(var,pos_dim2)
      start = (/ 1, i, j, my_record /)
      count = (/ dimlen(1), 1, 1, 1 /)
      IF(netcdf_read) THEN
        CALL nf(nf_get_vara_double(ncid, varid, start, count, buf(var_start)))
        IF(pos_dim2==3 .AND. pos_dim3==4) var(:,:,i,j) = RESHAPE(buf, (/ nproma, nblks /))
        IF(pos_dim2==2 .AND. pos_dim3==4) var(:,i,:,j) = RESHAPE(buf, (/ nproma, nblks /))
        IF(pos_dim2==2 .AND. pos_dim3==3) var(:,i,j,:) = RESHAPE(buf, (/ nproma, nblks /))
        IF(pos_dim2==1 .AND. pos_dim3==2) var(i,j,:,:) = RESHAPE(buf, (/ nproma, nblks /))
      ELSE
        IF(pos_dim2==3 .AND. pos_dim3==4) buf(:) = RESHAPE(var(:,:,i,j), (/ nproma*nblks /))
        IF(pos_dim2==2 .AND. pos_dim3==4) buf(:) = RESHAPE(var(:,i,:,j), (/ nproma*nblks /))
        IF(pos_dim2==2 .AND. pos_dim3==3) buf(:) = RESHAPE(var(:,i,j,:), (/ nproma*nblks /))
        IF(pos_dim2==1 .AND. pos_dim3==2) buf(:) = RESHAPE(var(i,j,:,:), (/ nproma*nblks /))
        CALL nf(nf_put_vara_double(ncid, varid, start, count, buf(var_start)))
      ENDIF
    ENDDO
    ENDDO

    DEALLOCATE(buf)

  END SUBROUTINE bvar_io_r3
  !-----------------------------------------------------------------------
  SUBROUTINE bvar_io_i1(pos_nproma, pos_nblks, var_name, var, opt_start)

    INTEGER, INTENT(IN) :: pos_nproma, pos_nblks
    CHARACTER(LEN=*), INTENT(IN) :: var_name
    INTEGER, INTENT(INOUT) :: var(:,:)
    INTEGER, INTENT(IN), OPTIONAL :: opt_start

    INTEGER :: var_start, varid, dimlen(1)
    INTEGER :: nblks, start(2), count(2)
    INTEGER, ALLOCATABLE :: buf(:)

    IF(PRESENT(opt_start)) THEN
      var_start = opt_start
    ELSE
      var_start = 1
    ENDIF

    CALL get_var_info(var_name, nf_int, 1, varid, dimlen)

    ! Do safety checks
    CALL pos_check(pos_nproma, pos_nblks, var_name, UBOUND(var), var_start, dimlen)

    nblks = UBOUND(var,pos_nblks)
    ALLOCATE(buf(nproma*nblks))
    buf(:) = 0

    start = (/ 1, my_record /)
    count = (/ dimlen(1), 1 /)
    IF(netcdf_read) THEN
      CALL nf(nf_get_vara_int(ncid, varid, start, count, buf(var_start)))
      var(:,:) = RESHAPE(buf, (/ nproma, nblks /))
    ELSE
      buf(:) = RESHAPE(var(:,:), (/ nproma*nblks /))
      CALL nf(nf_put_vara_int(ncid, varid, start, count, buf(var_start)))
    ENDIF

    DEALLOCATE(buf)

  END SUBROUTINE bvar_io_i1
  !-----------------------------------------------------------------------
  SUBROUTINE bvar_io_i2(pos_nproma, pos_nblks, var_name, var, opt_start)

    INTEGER, INTENT(IN) :: pos_nproma, pos_nblks
    CHARACTER(LEN=*), INTENT(IN) :: var_name
    INTEGER, INTENT(INOUT) :: var(:,:,:)
    INTEGER, INTENT(IN), OPTIONAL :: opt_start

    INTEGER :: var_start, varid, dimlen(2)
    INTEGER :: start(3), count(3), pos_dim2
    INTEGER :: nblks, i
    INTEGER, ALLOCATABLE :: buf(:)

    IF(PRESENT(opt_start)) THEN
      var_start = opt_start
    ELSE
      var_start = 1
    ENDIF

    CALL get_var_info(var_name, nf_int, 2, varid, dimlen)

    ! Do safety checks, get position of remaining dimensions
    CALL pos_check(pos_nproma, pos_nblks, var_name, UBOUND(var), var_start, dimlen, pos_dim2)

    nblks = UBOUND(var,pos_nblks)
    ALLOCATE(buf(nproma*nblks))
    buf(:) = 0

    DO i=1,UBOUND(var,pos_dim2)
      start = (/ 1, i, my_record /)
      count = (/ dimlen(1), 1, 1 /)
      IF(netcdf_read) THEN
        CALL nf(nf_get_vara_int(ncid, varid, start, count, buf(var_start)))
        IF(pos_dim2==1) var(i,:,:) = RESHAPE(buf, (/ nproma, nblks /))
        IF(pos_dim2==2) var(:,i,:) = RESHAPE(buf, (/ nproma, nblks /))
        IF(pos_dim2==3) var(:,:,i) = RESHAPE(buf, (/ nproma, nblks /))
      ELSE
        IF(pos_dim2==1) buf(:) = RESHAPE(var(i,:,:), (/ nproma*nblks /))
        IF(pos_dim2==2) buf(:) = RESHAPE(var(:,i,:), (/ nproma*nblks /))
        IF(pos_dim2==3) buf(:) = RESHAPE(var(:,:,i), (/ nproma*nblks /))
        CALL nf(nf_put_vara_int(ncid, varid, start, count, buf(var_start)))
      ENDIF
    ENDDO

    DEALLOCATE(buf)

  END SUBROUTINE bvar_io_i2
  !-----------------------------------------------------------------------
  SUBROUTINE bvar_io_l1(pos_nproma, pos_nblks, var_name, var, opt_start)

    INTEGER, INTENT(IN) :: pos_nproma, pos_nblks
    CHARACTER(LEN=*), INTENT(IN) :: var_name
    LOGICAL, INTENT(INOUT) :: var(:,:)
    INTEGER, INTENT(IN), OPTIONAL :: opt_start

    INTEGER :: var_start, varid, dimlen(1)
    INTEGER :: nblks, start(2), count(2)
    INTEGER, ALLOCATABLE :: buf(:)

    IF(PRESENT(opt_start)) THEN
      var_start = opt_start
    ELSE
      var_start = 1
    ENDIF

    CALL get_var_info(var_name, nf_int, 1, varid, dimlen)

    ! Do safety checks
    CALL pos_check(pos_nproma, pos_nblks, var_name, UBOUND(var), var_start, dimlen)

    nblks = UBOUND(var,pos_nblks)
    ALLOCATE(buf(nproma*nblks))
    buf(:) = 0

    start = (/ 1, my_record /)
    count = (/ dimlen(1), 1 /)
    IF(netcdf_read) THEN
      CALL nf(nf_get_vara_int(ncid, varid, start, count, buf(var_start)))
      var(:,:) = RESHAPE((buf(:)/=0), (/ nproma, nblks /))
    ELSE
      buf(:) = l2i(RESHAPE(var(:,:), (/ nproma*nblks /)))
      CALL nf(nf_put_vara_int(ncid, varid, start, count, buf(var_start)))
    ENDIF

    DEALLOCATE(buf)

  CONTAINS
    ELEMENTAL INTEGER FUNCTION l2i(x) ! logical to integer
      LOGICAL, INTENT(IN) :: x
      IF(x) THEN
        l2i = 1
      ELSE
        l2i = 0
      ENDIF
    END FUNCTION l2i

  END SUBROUTINE bvar_io_l1
  !
  !-----------------------------------------------------------------------
  !> bidx_io family:
  !! I/O of blocked index variables, i.e. depending on nproma/nblks
  !! pos_nproma:    position of nproma within the indices
  !! pos_nblks:     position of nblks within the indices
  !! var_name:      name of variable
  !! var_idx:       line index of index pair to be read/written
  !! var_blk:       block index of index pair to be read/written
  !! opt_start:     optional index of the first value for I/O in var_idx/var_blk
  !-----------------------------------------------------------------------
  !
  SUBROUTINE bidx_io_i1(pos_nproma, pos_nblks, var_name, var_idx, var_blk, opt_start)

    INTEGER, INTENT(IN) :: pos_nproma, pos_nblks
    CHARACTER(LEN=*), INTENT(IN) :: var_name
    INTEGER, INTENT(INOUT) :: var_idx(:,:), var_blk(:,:)
    INTEGER, INTENT(IN), OPTIONAL :: opt_start

    INTEGER :: var_start
    INTEGER, ALLOCATABLE :: idx(:,:)

    IF(PRESENT(opt_start)) THEN
      var_start = opt_start
    ELSE
      var_start = 1
    ENDIF

    ALLOCATE(idx(UBOUND(var_idx,1),UBOUND(var_idx,2)))

    IF(netcdf_read) THEN
      call bvar_io_i1(pos_nproma, pos_nblks, var_name, idx, var_start)
      var_idx = idx_no(idx)
      var_blk = blk_no(idx)
    ELSE
      idx = idx_1d(var_idx, var_blk)
      call bvar_io_i1(pos_nproma, pos_nblks, var_name, idx, var_start)
    ENDIF

    DEALLOCATE(idx)

  END SUBROUTINE bidx_io_i1
  !-----------------------------------------------------------------------
  SUBROUTINE bidx_io_i2(pos_nproma, pos_nblks, var_name, var_idx, var_blk, opt_start)

    INTEGER, INTENT(IN) :: pos_nproma, pos_nblks
    CHARACTER(LEN=*), INTENT(IN) :: var_name
    INTEGER, INTENT(INOUT) :: var_idx(:,:,:), var_blk(:,:,:)
    INTEGER, INTENT(IN), OPTIONAL :: opt_start

    INTEGER :: var_start
    INTEGER, ALLOCATABLE :: idx(:,:,:)

    IF(PRESENT(opt_start)) THEN
      var_start = opt_start
    ELSE
      var_start = 1
    ENDIF

    ALLOCATE(idx(UBOUND(var_idx,1),UBOUND(var_idx,2),UBOUND(var_idx,3)))

    IF(netcdf_read) THEN
      call bvar_io_i2(pos_nproma, pos_nblks, var_name, idx, var_start)
      var_idx = idx_no(idx)
      var_blk = blk_no(idx)
    ELSE
      idx = idx_1d(var_idx, var_blk)
      call bvar_io_i2(pos_nproma, pos_nblks, var_name, idx, var_start)
    ENDIF

    DEALLOCATE(idx)

  END SUBROUTINE bidx_io_i2
  !
  !-----------------------------------------------------------------------
  !> uvar_io family:
  !! I/O of unblocked variables, i.e. not depending on nproma/nblks
  !! var_name:      name of variable
  !! var:           data array to be read/written
  !-----------------------------------------------------------------------
  !
  SUBROUTINE uvar_io_r1(var_name, var)

    CHARACTER(LEN=*), INTENT(IN) :: var_name
    REAL(wp), INTENT(INOUT) :: var(:)

    INTEGER :: varid, dimlen(1), start(2), count(2)

    CALL get_var_info(var_name, nf_double, 1, varid, dimlen)
    IF(ANY(UBOUND(var) > dimlen)) &
      CALL finish(modname, var_name//': Variable dimensions and NetCDF dimensions do not conform')

    ! We can output var directly since values are stored contiguous

    start = (/ 1, my_record /)
    count = (/ UBOUND(var,1), 1 /)
    IF(netcdf_read) THEN
      CALL nf(nf_get_vara_double(ncid, varid, start, count, var))
    ELSE
      CALL nf(nf_put_vara_double(ncid, varid, start, count, var))
    ENDIF

  END SUBROUTINE uvar_io_r1
  !-----------------------------------------------------------------------
  SUBROUTINE uvar_io_r2(var_name, var)

    CHARACTER(LEN=*), INTENT(IN) :: var_name
    REAL(wp), INTENT(IN) :: var(:,:)

    INTEGER :: varid, dimlen(2), start(3), count(3)

    CALL get_var_info(var_name, nf_double, 2, varid, dimlen)
    IF(ANY(UBOUND(var) > dimlen)) &
      CALL finish(modname, var_name//': Variable dimensions and NetCDF dimensions do not conform')

    ! We can output var directly since values are stored contiguous

    start = (/ 1, 1, my_record /)
    count = (/ UBOUND(var,1), UBOUND(var,2), 1 /)
    IF(netcdf_read) THEN
      CALL nf(nf_get_var_double(ncid, varid, var))
    ELSE
      CALL nf(nf_put_var_double(ncid, varid, var))
    ENDIF

  END SUBROUTINE uvar_io_r2
  !-----------------------------------------------------------------------
  SUBROUTINE uvar_io_i1(var_name, var)

    CHARACTER(LEN=*), INTENT(IN) :: var_name
    INTEGER, INTENT(IN) :: var(:)

    INTEGER :: varid, dimlen(1), start(2), count(2)

    CALL get_var_info(var_name, nf_int, 1, varid, dimlen)
    IF(ANY(UBOUND(var) > dimlen)) &
      CALL finish(modname, var_name//': Variable dimensions and NetCDF dimensions do not conform')

    ! We can output var directly since values are stored contiguous

    start = (/ 1, my_record /)
    count = (/ UBOUND(var,1), 1 /)
    IF(netcdf_read) THEN
      CALL nf(nf_get_vara_int(ncid, varid, start, count, var))
    ELSE
      CALL nf(nf_put_vara_int(ncid, varid, start, count, var))
    ENDIF

  END SUBROUTINE uvar_io_i1
  !-----------------------------------------------------------------------
  SUBROUTINE uvar_io_i2(var_name, var)

    CHARACTER(LEN=*), INTENT(IN) :: var_name
    INTEGER, INTENT(IN) :: var(:,:)

    INTEGER :: varid, dimlen(2), start(3), count(3)

    CALL get_var_info(var_name, nf_int, 2, varid, dimlen)
    IF(ANY(UBOUND(var) > dimlen)) &
      CALL finish(modname, var_name//': Variable dimensions and NetCDF dimensions do not conform')

    ! We can output var directly since values are stored contiguous

    start = (/ 1, 1, my_record /)
    count = (/ UBOUND(var,1), UBOUND(var,2), 1 /)
    IF(netcdf_read) THEN
      CALL nf(nf_get_vara_int(ncid, varid, start, count, var))
    ELSE
      CALL nf(nf_put_vara_int(ncid, varid, start, count, var))
    ENDIF

  END SUBROUTINE uvar_io_i2
  !
  !-----------------------------------------------------------------------
  !> uidx_io family:
  !! I/O of unblocked index variables, i.e. not depending on nproma/nblks
  !! var_name:      name of variable
  !! var_idx:       line index of index pair to be read/written
  !! var_blk:       block index of index pair to be read/written
  !-----------------------------------------------------------------------
  !
  SUBROUTINE uidx_io_i1(var_name, var_idx, var_blk)

    CHARACTER(LEN=*), INTENT(IN) :: var_name
    INTEGER, INTENT(INOUT) :: var_idx(:), var_blk(:)

    INTEGER, ALLOCATABLE :: idx(:)

    ALLOCATE(idx(UBOUND(var_idx,1)))

    IF(netcdf_read) THEN
      call uvar_io_i1(var_name, idx)
      var_idx = idx_no(idx)
      var_blk = blk_no(idx)
    ELSE
      idx = idx_1d(var_idx, var_blk)
      call uvar_io_i1(var_name, idx)
    ENDIF

    DEALLOCATE(idx)

  END SUBROUTINE uidx_io_i1
  !-----------------------------------------------------------------------
  SUBROUTINE uidx_io_i2(var_name, var_idx, var_blk)

    CHARACTER(LEN=*), INTENT(IN) :: var_name
    INTEGER, INTENT(INOUT) :: var_idx(:,:), var_blk(:,:)

    INTEGER, ALLOCATABLE :: idx(:,:)

    ALLOCATE(idx(UBOUND(var_idx,1),UBOUND(var_idx,2)))

    IF(netcdf_read) THEN
      call uvar_io_i2(var_name, idx)
      var_idx = idx_no(idx)
      var_blk = blk_no(idx)
    ELSE
      idx = idx_1d(var_idx, var_blk)
      call uvar_io_i2(var_name, idx)
    ENDIF

    DEALLOCATE(idx)

  END SUBROUTINE uidx_io_i2
  !
  !-----------------------------------------------------------------------
  ! Handling of communication patterns
  !-----------------------------------------------------------------------
  !
  SUBROUTINE def_comm_pat(base_name, pat)

    CHARACTER(LEN=*), INTENT(IN) :: base_name
    TYPE(t_comm_pattern), INTENT(IN) :: pat

    INTEGER :: length, dim_length

    ! Some communication patterns are not allocated at all

    IF(.NOT. ALLOCATED(pat%recv_limits)) RETURN

    ! Communication patterns are stuffed into one integer array in the NetCDF file
    ! which has the following length:

    length = 3 + 2*(p_n_work+1) + 2*pat%n_pnts + 3*pat%np_send + 3*pat%np_recv + pat%n_send

    IF(use_one_file) length = p_max(length, comm=p_comm_work)

    IF(output_defines) THEN
      CALL def_dim(base_name//'.length', length, dim_length)
      CALL def_var(base_name//'.data', nf_int, dim_length)
    ENDIF

  END SUBROUTINE def_comm_pat
  !-----------------------------------------------------------------------
  SUBROUTINE comm_pat_io(base_name, pat)

    CHARACTER(LEN=*), INTENT(IN) :: base_name
    TYPE(t_comm_pattern), INTENT(INOUT) :: pat

    INTEGER :: stat, dimid, n, length
    INTEGER, ALLOCATABLE :: pat_data(:)

    IF(netcdf_read) THEN

      ! Check if pattern length dimension is in netCDF file,
      ! otherwise this pattern will not be allocated
      stat = nf_inq_dimid (ncid, TRIM(prefix)//base_name//'.length', dimid)
      IF (stat /= nf_noerr) RETURN
      ! Get length of data
      CALL nf(nf_inq_dimlen(ncid, dimid, length))

      ! Allocate and read data
      ALLOCATE(pat_data(length))
      CALL uvar_io(base_name//'.data', pat_data)

      ! Distribute data

      pat%n_pnts  = pat_data(1)
      pat%np_recv = pat_data(2)
      pat%np_send = pat_data(3)
      n = 3

      ALLOCATE(pat%recv_limits(0:p_n_work))
      ALLOCATE(pat%send_limits(0:p_n_work))
      pat%recv_limits(0:p_n_work) = pat_data(n+1:n+p_n_work+1); n = n + p_n_work+1
      pat%send_limits(0:p_n_work) = pat_data(n+1:n+p_n_work+1); n = n + p_n_work+1

      pat%n_recv = pat%recv_limits(p_n_work)
      pat%n_send = pat%send_limits(p_n_work)

      ! Safety check
      IF(3 + 2*(p_n_work+1) + 2*pat%n_pnts + 3*pat%np_send + 3*pat%np_recv + pat%n_send > length) &
        & CALL finish(modname,'Illegal length for '//TRIM(prefix)//TRIM(base_name))

      ALLOCATE(pat%recv_src(pat%n_pnts))
      ALLOCATE(pat%recv_dst_blk(pat%n_pnts))
      ALLOCATE(pat%recv_dst_idx(pat%n_pnts))

      pat%recv_src(:) = pat_data(n+1:n+pat%n_pnts)
      n = n + pat%n_pnts
      pat%recv_dst_idx  = idx_no(pat_data(n+1:n+pat%n_pnts))
      pat%recv_dst_blk  = blk_no(pat_data(n+1:n+pat%n_pnts))
      n = n + pat%n_pnts

      ALLOCATE (pat%pelist_send(pat%np_send),   pat%pelist_recv(pat%np_recv),   &
                pat%send_startidx(pat%np_send), pat%recv_startidx(pat%np_recv), &
                pat%send_count(pat%np_send),    pat%recv_count(pat%np_recv)     )

      pat%pelist_send(:)   = pat_data(n+1:n+pat%np_send); n = n + pat%np_send
      pat%send_startidx(:) = pat_data(n+1:n+pat%np_send); n = n + pat%np_send
      pat%send_count(:)    = pat_data(n+1:n+pat%np_send); n = n + pat%np_send

      pat%pelist_recv(:)   = pat_data(n+1:n+pat%np_recv); n = n + pat%np_recv
      pat%recv_startidx(:) = pat_data(n+1:n+pat%np_recv); n = n + pat%np_recv
      pat%recv_count(:)    = pat_data(n+1:n+pat%np_recv); n = n + pat%np_recv

      ALLOCATE(pat%send_src_blk(pat%n_send))
      ALLOCATE(pat%send_src_idx(pat%n_send))

      pat%send_src_idx(:)  = idx_no(pat_data(n+1:n+pat%n_send))
      pat%send_src_blk(:)  = blk_no(pat_data(n+1:n+pat%n_send))
      n = n + pat%n_send

      DEALLOCATE(pat_data)

    ELSE

      ! Do nothing if the pattern is not allocated
      IF(.NOT. ALLOCATED(pat%recv_limits)) RETURN

      ! Communication patterns are stuffed into one integer array in the NetCDF file
      ! which has the following length:

      length = 3 + 2*(p_n_work+1) + 2*pat%n_pnts + 3*pat%np_send + 3*pat%np_recv + pat%n_send
      ALLOCATE(pat_data(length))

      pat_data(1) = pat%n_pnts
      pat_data(2) = pat%np_recv
      pat_data(3) = pat%np_send
      n = 3

      pat_data(n+1:n+p_n_work+1) = pat%recv_limits(0:p_n_work); n = n + p_n_work+1
      pat_data(n+1:n+p_n_work+1) = pat%send_limits(0:p_n_work); n = n + p_n_work+1

      pat_data(n+1:n+pat%n_pnts) = pat%recv_src(:)
      n = n + pat%n_pnts
      pat_data(n+1:n+pat%n_pnts) = idx_1d(pat%recv_dst_idx(:),pat%recv_dst_blk(:))
      n = n + pat%n_pnts

      pat_data(n+1:n+pat%np_send) = pat%pelist_send(:);   n = n + pat%np_send
      pat_data(n+1:n+pat%np_send) = pat%send_startidx(:); n = n + pat%np_send
      pat_data(n+1:n+pat%np_send) = pat%send_count(:);    n = n + pat%np_send

      pat_data(n+1:n+pat%np_recv) = pat%pelist_recv(:);   n = n + pat%np_recv
      pat_data(n+1:n+pat%np_recv) = pat%recv_startidx(:); n = n + pat%np_recv
      pat_data(n+1:n+pat%np_recv) = pat%recv_count(:);    n = n + pat%np_recv

      pat_data(n+1:n+pat%n_send)  = idx_1d(pat%send_src_idx(:),pat%send_src_blk(:))
      n = n + pat%n_send

      CALL uvar_io(base_name//'.data', pat_data)

      DEALLOCATE(pat_data)

    ENDIF

  END SUBROUTINE comm_pat_io
  !
  !-------------------------------------------------------------------------
  ! Routines for definning and I/O of patches and states
  !-------------------------------------------------------------------------
  !
  !> Determines the maximum dimensions of cells/edges/verts for NetCDF output

  SUBROUTINE set_max_patch_dims(p)

    TYPE(t_patch), INTENT(INOUT) :: p

    IF(use_one_file) THEN
      max_patch_cells = p_max(p%n_patch_cells, comm=p_comm_work)
      max_patch_edges = p_max(p%n_patch_edges, comm=p_comm_work)
      max_patch_verts = p_max(p%n_patch_verts, comm=p_comm_work)
    ELSE
      max_patch_cells = p%n_patch_cells
      max_patch_edges = p%n_patch_edges
      max_patch_verts = p%n_patch_verts
    ENDIF

  END SUBROUTINE set_max_patch_dims

  !-------------------------------------------------------------------------
  !
  !> Defines patch variables for NetCDF

  SUBROUTINE def_patch(p, lfull)

    TYPE(t_patch), INTENT(INOUT) :: p
    LOGICAL, INTENT(IN)          :: lfull ! Flag if full or basic patch is to be dumped

    INTEGER :: dim_max_childdom
    INTEGER :: dim_ncells_g, dim_nedges_g, dim_nverts_g


    ! There is nothing to do if not in define mode
    IF(.NOT. output_defines) RETURN

    ! Dimensions specific for this patch (with prefix)

    CALL def_dim('max_childdom', p%max_childdom, dim_max_childdom)

    dim_ncells = -1
    dim_nedges = -1
    dim_nverts = -1

    IF(max_patch_cells>0) &
      CALL def_dim('max_patch_cells', max_patch_cells, dim_ncells)
    IF(max_patch_edges>0) &
      CALL def_dim('max_patch_edges', max_patch_edges, dim_nedges)
    IF(max_patch_verts>0) &
      CALL def_dim('max_patch_verts', max_patch_verts, dim_nverts)

    CALL def_dim('n_patch_cells_g', p%n_patch_cells_g, dim_ncells_g)
    CALL def_dim('n_patch_edges_g', p%n_patch_edges_g, dim_nedges_g)
    CALL def_dim('n_patch_verts_g', p%n_patch_verts_g, dim_nverts_g)

    ! Variables storing the actual number of cells/edges/verts

    CALL def_var('n_patch_cells', nf_int)
    CALL def_var('n_patch_edges', nf_int)
    CALL def_var('n_patch_verts', nf_int)

    ! Variables for patch

    IF(max_patch_cells>0) THEN
      CALL def_var('patch.cells.num_edges',        nf_int   , dim_ncells)
      CALL def_var('patch.cells.parent_index',     nf_int   , dim_ncells)
      CALL def_var('patch.cells.pc_idx',           nf_int   , dim_ncells)
      CALL def_var('patch.cells.child_index',      nf_int   , dim_ncells, dim_4)
      CALL def_var('patch.cells.child_id',         nf_int   , dim_ncells)
      CALL def_var('patch.cells.neighbor_index',   nf_int   , dim_ncells, dim_nverts_per_cell)
      CALL def_var('patch.cells.edge_index',       nf_int   , dim_ncells, dim_nverts_per_cell)
      CALL def_var('patch.cells.vertex_index',     nf_int   , dim_ncells, dim_nverts_per_cell)
      CALL def_var('patch.cells.center.lon',       nf_double, dim_ncells)
      CALL def_var('patch.cells.center.lat',       nf_double, dim_ncells)
      CALL def_var('patch.cells.refin_ctrl',       nf_int   , dim_ncells)
      CALL def_var('patch.cells.decomp_domain',    nf_int   , dim_ncells)
      CALL def_var('patch.cells.glb_index',        nf_int   , dim_ncells)

     IF(lfull) THEN
      CALL def_var('patch.cells.phys_id',          nf_int   , dim_ncells)
      CALL def_var('patch.cells.edge_orientation', nf_double, dim_ncells, dim_nverts_per_cell)
      CALL def_var('patch.cells.area',             nf_double, dim_ncells)
      CALL def_var('patch.cells.f_c',              nf_double, dim_ncells)
      CALL def_var('patch.cells.owner_mask',       nf_int   , dim_ncells)
      CALL def_var('patch.cells.owner_local',      nf_int   , dim_ncells)
     ENDIF
    ENDIF
    CALL def_var  ('patch.cells.start_index',      nf_int   , dim_nrlcell, dim_max_childdom)
    CALL def_var  ('patch.cells.end_index',        nf_int   , dim_nrlcell, dim_max_childdom)
    ! patch.cells.loc_index is not saved since it is a global array and easy to calculate
    ! owner_g is stored only once in the file since it is identical on all procs
    CALL def_var  ('patch.cells.owner_g',          nf_int   , dim_ncells_g, add_unlim=.FALSE.)


    IF(max_patch_edges>0) THEN
      CALL def_var('patch.edges.parent_index',           nf_int   , dim_nedges)
      CALL def_var('patch.edges.pc_idx',                 nf_int   , dim_nedges)
      CALL def_var('patch.edges.child_index',            nf_int   , dim_nedges, dim_4)
      CALL def_var('patch.edges.child_id',               nf_int   , dim_nedges)
      CALL def_var('patch.edges.refin_ctrl',             nf_int   , dim_nedges)
      CALL def_var('patch.edges.decomp_domain',          nf_int   , dim_nedges)
      CALL def_var('patch.edges.glb_index',              nf_int   , dim_nedges)

     IF(lfull) THEN
      CALL def_var('patch.edges.phys_id',                nf_int   , dim_nedges)
      CALL def_var('patch.edges.cell_index',             nf_int   , dim_nedges, dim_2)
      CALL def_var('patch.edges.vertex_index',           nf_int   , dim_nedges, dim_4)
      CALL def_var('patch.edges.system_orientation',     nf_double, dim_nedges)
      CALL def_var('patch.edges.quad_index',             nf_int   , dim_nedges, dim_4)
      CALL def_var('patch.edges.quad_orientation',       nf_double, dim_nedges, dim_4)
      CALL def_var('patch.edges.center.lon',             nf_double, dim_nedges)
      CALL def_var('patch.edges.center.lat',             nf_double, dim_nedges)
      CALL def_var('patch.edges.primal_normal.v1',       nf_double, dim_nedges)
      CALL def_var('patch.edges.primal_normal.v2',       nf_double, dim_nedges)
      CALL def_var('patch.edges.primal_cart_normal.x1',  nf_double, dim_nedges)
      CALL def_var('patch.edges.primal_cart_normal.x2',  nf_double, dim_nedges)
      CALL def_var('patch.edges.primal_cart_normal.x3',  nf_double, dim_nedges)
      CALL def_var('patch.edges.dual_cart_normal.x1',    nf_double, dim_nedges)
      CALL def_var('patch.edges.dual_cart_normal.x2',    nf_double, dim_nedges)
      CALL def_var('patch.edges.dual_cart_normal.x3',    nf_double, dim_nedges)
      CALL def_var('patch.edges.dual_normal.v1',         nf_double, dim_nedges)
      CALL def_var('patch.edges.dual_normal.v2',         nf_double, dim_nedges)
      CALL def_var('patch.edges.primal_normal_cell.v1',  nf_double, dim_nedges, dim_2)
      CALL def_var('patch.edges.primal_normal_cell.v2',  nf_double, dim_nedges, dim_2)
      CALL def_var('patch.edges.dual_normal_cell.v1',    nf_double, dim_nedges, dim_2)
      CALL def_var('patch.edges.dual_normal_cell.v2',    nf_double, dim_nedges, dim_2)
      CALL def_var('patch.edges.primal_normal_vert.v1',  nf_double, dim_nedges, dim_4)
      CALL def_var('patch.edges.primal_normal_vert.v2',  nf_double, dim_nedges, dim_4)
      CALL def_var('patch.edges.dual_normal_vert.v1',    nf_double, dim_nedges, dim_4)
      CALL def_var('patch.edges.dual_normal_vert.v2',    nf_double, dim_nedges, dim_4)
      CALL def_var('patch.edges.primal_edge_length',     nf_double, dim_nedges)
      CALL def_var('patch.edges.inv_primal_edge_length', nf_double, dim_nedges)
      CALL def_var('patch.edges.dual_edge_length',       nf_double, dim_nedges)
      CALL def_var('patch.edges.inv_dual_edge_length',   nf_double, dim_nedges)
      CALL def_var('patch.edges.edge_vert_length',       nf_double, dim_nedges, dim_2)
      CALL def_var('patch.edges.edge_cell_length',       nf_double, dim_nedges, dim_2)
      CALL def_var('patch.edges.inv_vert_vert_length',   nf_double, dim_nedges)
      CALL def_var('patch.edges.area_edge',              nf_double, dim_nedges)
      CALL def_var('patch.edges.quad_area',              nf_double, dim_nedges)
      CALL def_var('patch.edges.f_e',                    nf_double, dim_nedges)
      CALL def_var('patch.edges.owner_mask',             nf_int   , dim_nedges)
     ENDIF
    ENDIF
    CALL def_var  ('patch.edges.start_index',            nf_int, dim_nrledge, dim_max_childdom)
    CALL def_var  ('patch.edges.end_index',              nf_int, dim_nrledge, dim_max_childdom)
    ! patch.edges.loc_index is not saved since it is a global array and easy to calculate
    ! owner_g is stored only once in the file since it is identical on all procs
    CALL def_var  ('patch.edges.owner_g', nf_int, dim_nedges_g, add_unlim=.FALSE.)


    IF(max_patch_verts>0) THEN
      CALL def_var('patch.verts.vertex.lon',       nf_double, dim_nverts)
      CALL def_var('patch.verts.vertex.lat',       nf_double, dim_nverts)
      CALL def_var('patch.verts.refin_ctrl',       nf_int   , dim_nverts)
      CALL def_var('patch.verts.decomp_domain',    nf_int   , dim_nverts)
      CALL def_var('patch.verts.glb_index',        nf_int   , dim_nverts)

     IF(lfull) THEN
      CALL def_var('patch.verts.phys_id',          nf_int   , dim_nverts)
      CALL def_var('patch.verts.neighbor_index',   nf_int   , dim_nverts, dim_nedges_per_vert)
      CALL def_var('patch.verts.cell_index',       nf_int   , dim_nverts, dim_nedges_per_vert)
      CALL def_var('patch.verts.edge_index',       nf_int   , dim_nverts, dim_nedges_per_vert)
      CALL def_var('patch.verts.edge_orientation', nf_double, dim_nverts, dim_nedges_per_vert)
      CALL def_var('patch.verts.num_edges',        nf_int   , dim_nverts)
      CALL def_var('patch.verts.dual_area',        nf_double, dim_nverts)
      CALL def_var('patch.verts.f_v',              nf_double, dim_nverts)
      CALL def_var('patch.verts.owner_mask',       nf_int   , dim_nverts)
     ENDIF
    ENDIF
    CALL def_var  ('patch.verts.start_index',      nf_int   , dim_nrlvert, dim_max_childdom)
    CALL def_var  ('patch.verts.end_index',        nf_int   , dim_nrlvert, dim_max_childdom)
    ! patch.verts.loc_index is not saved since it is a global array and easy to calculate
    ! owner_g is stored only once in the file since it is identical on all procs
    CALL def_var  ('patch.verts.owner_g',          nf_int   , dim_nverts_g, add_unlim=.FALSE.)

  END SUBROUTINE def_patch

  !-------------------------------------------------------------------------
  !
  !> Defines the communication patterns to be put into NetCDF

  SUBROUTINE def_comm_patterns(p)

    TYPE(t_patch), INTENT(IN) :: p

    CALL def_comm_pat('comm_pat_c',                   p%comm_pat_c)
    CALL def_comm_pat('comm_pat_c1',                  p%comm_pat_c1)
    CALL def_comm_pat('comm_pat_e',                   p%comm_pat_e)
    CALL def_comm_pat('comm_pat_v',                   p%comm_pat_v)
    CALL def_comm_pat('comm_pat_interpolation_c',     p%comm_pat_interpolation_c)
    CALL def_comm_pat('comm_pat_interpol_vec_grf.1',  p%comm_pat_interpol_vec_grf(1))
    CALL def_comm_pat('comm_pat_interpol_vec_grf.2',  p%comm_pat_interpol_vec_grf(2))
    CALL def_comm_pat('comm_pat_interpol_vec_grf.3',  p%comm_pat_interpol_vec_grf(3))
    CALL def_comm_pat('comm_pat_interpol_vec_grf.4',  p%comm_pat_interpol_vec_grf(4))
    CALL def_comm_pat('comm_pat_interpol_scal_grf.1', p%comm_pat_interpol_scal_grf(1))
    CALL def_comm_pat('comm_pat_interpol_scal_grf.2', p%comm_pat_interpol_scal_grf(2))
    CALL def_comm_pat('comm_pat_interpol_scal_grf.3', p%comm_pat_interpol_scal_grf(3))
    CALL def_comm_pat('comm_pat_interpol_scal_grf.4', p%comm_pat_interpol_scal_grf(4))
    CALL def_comm_pat('comm_pat_interpol_vec_ubc.1',  p%comm_pat_interpol_vec_ubc(1))
    CALL def_comm_pat('comm_pat_interpol_vec_ubc.2',  p%comm_pat_interpol_vec_ubc(2))
    CALL def_comm_pat('comm_pat_interpol_vec_ubc.3',  p%comm_pat_interpol_vec_ubc(3))
    CALL def_comm_pat('comm_pat_interpol_vec_ubc.4',  p%comm_pat_interpol_vec_ubc(4))
    CALL def_comm_pat('comm_pat_interpol_scal_ubc.1', p%comm_pat_interpol_scal_ubc(1))
    CALL def_comm_pat('comm_pat_interpol_scal_ubc.2', p%comm_pat_interpol_scal_ubc(2))
    CALL def_comm_pat('comm_pat_interpol_scal_ubc.3', p%comm_pat_interpol_scal_ubc(3))
    CALL def_comm_pat('comm_pat_interpol_scal_ubc.4', p%comm_pat_interpol_scal_ubc(4))
    CALL def_comm_pat('comm_pat_gather_c',            p%comm_pat_gather_c)
    CALL def_comm_pat('comm_pat_gather_e',            p%comm_pat_gather_e)
    CALL def_comm_pat('comm_pat_gather_v',            p%comm_pat_gather_v)

    CALL def_comm_pat('comm_pat_glb_to_loc_c',        p%comm_pat_glb_to_loc_c)
    CALL def_comm_pat('comm_pat_glb_to_loc_e',        p%comm_pat_glb_to_loc_e)
    CALL def_comm_pat('comm_pat_loc_to_glb_c_fbk',    p%comm_pat_loc_to_glb_c_fbk)
    CALL def_comm_pat('comm_pat_loc_to_glb_e_fbk',    p%comm_pat_loc_to_glb_e_fbk)

  END SUBROUTINE def_comm_patterns

  !-------------------------------------------------------------------------
  !
  !> Dumps/restores patch data to/from NetCDF file.

  SUBROUTINE patch_io(p, lfull)

    TYPE(t_patch), INTENT(INOUT) :: p
    LOGICAL, INTENT(IN)          :: lfull

    INTEGER :: varid

    CALL store_proc_dependent_dimensions(p)

    IF(.NOT.netcdf_read) THEN
      ! Output the number of cells/edges/verts on the current processor
      CALL nf(nf_inq_varid(ncid, TRIM(prefix)//'n_patch_cells', varid))
      CALL nf(nf_put_var1_int(ncid, varid, my_record, p%n_patch_cells))
      CALL nf(nf_inq_varid(ncid, TRIM(prefix)//'n_patch_edges', varid))
      CALL nf(nf_put_var1_int(ncid, varid, my_record, p%n_patch_edges))
      CALL nf(nf_inq_varid(ncid, TRIM(prefix)//'n_patch_verts', varid))
      CALL nf(nf_put_var1_int(ncid, varid, my_record, p%n_patch_verts))
    ENDIF

    IF(p%n_patch_cells>0) THEN
      CALL bvar_io(1,2,'patch.cells.num_edges',        p%cells%num_edges)
      CALL bidx_io(1,2,'patch.cells.parent_index',     p%cells%parent_idx,   p%cells%parent_blk)
      CALL bvar_io(1,2,'patch.cells.pc_idx',           p%cells%pc_idx)
      CALL bidx_io(1,2,'patch.cells.child_index',      p%cells%child_idx,    p%cells%child_blk)
      CALL bvar_io(1,2,'patch.cells.child_id',         p%cells%child_id)
      CALL bidx_io(1,2,'patch.cells.neighbor_index',   p%cells%neighbor_idx, p%cells%neighbor_blk)
      CALL bidx_io(1,2,'patch.cells.edge_index',       p%cells%edge_idx,     p%cells%edge_blk)
      CALL bidx_io(1,2,'patch.cells.vertex_index',     p%cells%vertex_idx,   p%cells%vertex_blk)
      CALL bvar_io(1,2,'patch.cells.center.lon',       p%cells%center(:,:)%lon)
      CALL bvar_io(1,2,'patch.cells.center.lat',       p%cells%center(:,:)%lat)
      CALL bvar_io(1,2,'patch.cells.refin_ctrl',       p%cells%refin_ctrl)
      CALL bvar_io(1,2,'patch.cells.decomp_domain',    p%cells%decomp_domain)
      CALL uvar_io(    'patch.cells.glb_index',        p%cells%glb_index)

     IF(lfull) THEN
      CALL bvar_io(1,2,'patch.cells.phys_id',          p%cells%phys_id)
      CALL bvar_io(1,2,'patch.cells.edge_orientation', p%cells%edge_orientation)
      CALL bvar_io(1,2,'patch.cells.area',             p%cells%area)
      CALL bvar_io(1,2,'patch.cells.f_c',              p%cells%f_c)
      CALL bvar_io(1,2,'patch.cells.owner_mask',       p%cells%owner_mask)
      CALL uvar_io(    'patch.cells.owner_local',      p%cells%owner_local)
     ENDIF
    ENDIF
    CALL uidx_io  (    'patch.cells.start_index',      p%cells%start_idx,p%cells%start_blk)
    CALL uidx_io  (    'patch.cells.end_index',        p%cells%end_idx,  p%cells%end_blk)
    IF(netcdf_read .OR. my_record==1) & ! Note: owner_g has to be stored only once
    CALL uvar_io  (    'patch.cells.owner_g',          p%cells%owner_g)


    IF(p%n_patch_edges>0) THEN
      CALL bidx_io(1,2,'patch.edges.parent_index',     p%edges%parent_idx, p%edges%parent_blk)
      CALL bvar_io(1,2,'patch.edges.pc_idx',           p%edges%pc_idx)
      CALL bidx_io(1,2,'patch.edges.child_index',      p%edges%child_idx,  p%edges%child_blk)
      CALL bvar_io(1,2,'patch.edges.child_id',         p%edges%child_id)
      CALL bvar_io(1,2,'patch.edges.refin_ctrl',             p%edges%refin_ctrl)
      CALL bvar_io(1,2,'patch.edges.decomp_domain',          p%edges%decomp_domain)
      CALL uvar_io(    'patch.edges.glb_index',              p%edges%glb_index)

     IF(lfull) THEN
      CALL bvar_io(1,2,'patch.edges.phys_id',          p%edges%phys_id)
      CALL bidx_io(1,2,'patch.edges.cell_index',       p%edges%cell_idx,   p%edges%cell_blk)
      CALL bidx_io(1,2,'patch.edges.vertex_index',     p%edges%vertex_idx, p%edges%vertex_blk)
      CALL bvar_io(1,2,'patch.edges.system_orientation',     p%edges%system_orientation)
      CALL bidx_io(1,2,'patch.edges.quad_index',             p%edges%quad_idx,   p%edges%quad_blk)
      CALL bvar_io(1,2,'patch.edges.quad_orientation',       p%edges%quad_orientation)
      CALL bvar_io(1,2,'patch.edges.center.lon',             p%edges%center(:,:)%lon)
      CALL bvar_io(1,2,'patch.edges.center.lat',             p%edges%center(:,:)%lat)
      CALL bvar_io(1,2,'patch.edges.primal_normal.v1',       p%edges%primal_normal(:,:)%v1)
      CALL bvar_io(1,2,'patch.edges.primal_normal.v2',       p%edges%primal_normal(:,:)%v2)
      CALL bvar_io(1,2,'patch.edges.primal_cart_normal.x1',  p%edges%primal_cart_normal(:,:)%x(1))
      CALL bvar_io(1,2,'patch.edges.primal_cart_normal.x2',  p%edges%primal_cart_normal(:,:)%x(2))
      CALL bvar_io(1,2,'patch.edges.primal_cart_normal.x3',  p%edges%primal_cart_normal(:,:)%x(3))
      CALL bvar_io(1,2,'patch.edges.dual_cart_normal.x1',    p%edges%dual_cart_normal(:,:)%x(1))
      CALL bvar_io(1,2,'patch.edges.dual_cart_normal.x2',    p%edges%dual_cart_normal(:,:)%x(2))
      CALL bvar_io(1,2,'patch.edges.dual_cart_normal.x3',    p%edges%dual_cart_normal(:,:)%x(3))
      CALL bvar_io(1,2,'patch.edges.dual_normal.v1',         p%edges%dual_normal(:,:)%v1)
      CALL bvar_io(1,2,'patch.edges.dual_normal.v2',         p%edges%dual_normal(:,:)%v2)
      CALL bvar_io(1,2,'patch.edges.primal_normal_cell.v1',  p%edges%primal_normal_cell(:,:,:)%v1)
      CALL bvar_io(1,2,'patch.edges.primal_normal_cell.v2',  p%edges%primal_normal_cell(:,:,:)%v2)
      CALL bvar_io(1,2,'patch.edges.dual_normal_cell.v1',    p%edges%dual_normal_cell(:,:,:)%v1)
      CALL bvar_io(1,2,'patch.edges.dual_normal_cell.v2',    p%edges%dual_normal_cell(:,:,:)%v2)
      CALL bvar_io(1,2,'patch.edges.primal_normal_vert.v1',  p%edges%primal_normal_vert(:,:,:)%v1)
      CALL bvar_io(1,2,'patch.edges.primal_normal_vert.v2',  p%edges%primal_normal_vert(:,:,:)%v2)
      CALL bvar_io(1,2,'patch.edges.dual_normal_vert.v1',    p%edges%dual_normal_vert(:,:,:)%v1)
      CALL bvar_io(1,2,'patch.edges.dual_normal_vert.v2',    p%edges%dual_normal_vert(:,:,:)%v2)
      CALL bvar_io(1,2,'patch.edges.primal_edge_length',     p%edges%primal_edge_length)
      CALL bvar_io(1,2,'patch.edges.inv_primal_edge_length', p%edges%inv_primal_edge_length)
      CALL bvar_io(1,2,'patch.edges.dual_edge_length',       p%edges%dual_edge_length)
      CALL bvar_io(1,2,'patch.edges.inv_dual_edge_length',   p%edges%inv_dual_edge_length)
      CALL bvar_io(1,2,'patch.edges.edge_vert_length',       p%edges%edge_vert_length)
      CALL bvar_io(1,2,'patch.edges.edge_cell_length',       p%edges%edge_cell_length)
      CALL bvar_io(1,2,'patch.edges.inv_vert_vert_length',   p%edges%inv_vert_vert_length)
      CALL bvar_io(1,2,'patch.edges.area_edge',              p%edges%area_edge)
      CALL bvar_io(1,2,'patch.edges.quad_area',              p%edges%quad_area)
      CALL bvar_io(1,2,'patch.edges.f_e',                    p%edges%f_e)
      CALL bvar_io(1,2,'patch.edges.owner_mask',             p%edges%owner_mask)
     ENDIF
    ENDIF
    CALL uidx_io  (    'patch.edges.start_index',            p%edges%start_idx, p%edges%start_blk)
    CALL uidx_io  (    'patch.edges.end_index',              p%edges%end_idx,   p%edges%end_blk)
    IF(netcdf_read .OR. my_record==1) & ! Note: owner_g has to be stored only once
    CALL uvar_io  (    'patch.edges.owner_g',                p%edges%owner_g)

    IF(p%n_patch_verts>0) THEN
      CALL bvar_io(1,2,'patch.verts.vertex.lon',       p%verts%vertex(:,:)%lon)
      CALL bvar_io(1,2,'patch.verts.vertex.lat',       p%verts%vertex(:,:)%lat)
      CALL bvar_io(1,2,'patch.verts.refin_ctrl',       p%verts%refin_ctrl)
      CALL bvar_io(1,2,'patch.verts.decomp_domain',    p%verts%decomp_domain)
      CALL uvar_io(    'patch.verts.glb_index',        p%verts%glb_index)

     IF(lfull) THEN
      CALL bvar_io(1,2,'patch.verts.phys_id',          p%verts%phys_id)
      CALL bidx_io(1,2,'patch.verts.neighbor_index',   p%verts%neighbor_idx, p%verts%neighbor_blk)
      CALL bidx_io(1,2,'patch.verts.cell_index',       p%verts%cell_idx, p%verts%cell_blk)
      CALL bidx_io(1,2,'patch.verts.edge_index',       p%verts%edge_idx, p%verts%edge_blk)
      CALL bvar_io(1,2,'patch.verts.edge_orientation', p%verts%edge_orientation)
      CALL bvar_io(1,2,'patch.verts.num_edges',        p%verts%num_edges)
      CALL bvar_io(1,2,'patch.verts.dual_area',        p%verts%dual_area)
      CALL bvar_io(1,2,'patch.verts.f_v',              p%verts%f_v)
      CALL bvar_io(1,2,'patch.verts.owner_mask',       p%verts%owner_mask)
     ENDIF
    ENDIF
    CALL uidx_io  (    'patch.verts.start_index',      p%verts%start_idx, p%verts%start_blk)
    CALL uidx_io  (    'patch.verts.end_index',        p%verts%end_idx,   p%verts%end_blk)
    IF(netcdf_read .OR. my_record==1) & ! Note: owner_g has to be stored only once
    CALL uvar_io  (    'patch.verts.owner_g',          p%verts%owner_g)

    IF(lfull) THEN
      CALL comm_pat_io('comm_pat_c',                   p%comm_pat_c)
      CALL comm_pat_io('comm_pat_c1',                  p%comm_pat_c1)
      CALL comm_pat_io('comm_pat_e',                   p%comm_pat_e)
      CALL comm_pat_io('comm_pat_v',                   p%comm_pat_v)
      CALL comm_pat_io('comm_pat_interpolation_c',     p%comm_pat_interpolation_c)
      CALL comm_pat_io('comm_pat_interpol_vec_grf.1',  p%comm_pat_interpol_vec_grf(1))
      CALL comm_pat_io('comm_pat_interpol_vec_grf.2',  p%comm_pat_interpol_vec_grf(2))
      CALL comm_pat_io('comm_pat_interpol_vec_grf.3',  p%comm_pat_interpol_vec_grf(3))
      CALL comm_pat_io('comm_pat_interpol_vec_grf.4',  p%comm_pat_interpol_vec_grf(4))
      CALL comm_pat_io('comm_pat_interpol_scal_grf.1', p%comm_pat_interpol_scal_grf(1))
      CALL comm_pat_io('comm_pat_interpol_scal_grf.2', p%comm_pat_interpol_scal_grf(2))
      CALL comm_pat_io('comm_pat_interpol_scal_grf.3', p%comm_pat_interpol_scal_grf(3))
      CALL comm_pat_io('comm_pat_interpol_scal_grf.4', p%comm_pat_interpol_scal_grf(4))
      CALL comm_pat_io('comm_pat_interpol_vec_ubc.1',  p%comm_pat_interpol_vec_ubc(1))
      CALL comm_pat_io('comm_pat_interpol_vec_ubc.2',  p%comm_pat_interpol_vec_ubc(2))
      CALL comm_pat_io('comm_pat_interpol_vec_ubc.3',  p%comm_pat_interpol_vec_ubc(3))
      CALL comm_pat_io('comm_pat_interpol_vec_ubc.4',  p%comm_pat_interpol_vec_ubc(4))
      CALL comm_pat_io('comm_pat_interpol_scal_ubc.1', p%comm_pat_interpol_scal_ubc(1))
      CALL comm_pat_io('comm_pat_interpol_scal_ubc.2', p%comm_pat_interpol_scal_ubc(2))
      CALL comm_pat_io('comm_pat_interpol_scal_ubc.3', p%comm_pat_interpol_scal_ubc(3))
      CALL comm_pat_io('comm_pat_interpol_scal_ubc.4', p%comm_pat_interpol_scal_ubc(4))
      CALL comm_pat_io('comm_pat_gather_c',            p%comm_pat_gather_c)
      CALL comm_pat_io('comm_pat_gather_e',            p%comm_pat_gather_e)
      CALL comm_pat_io('comm_pat_gather_v',            p%comm_pat_gather_v)

      CALL comm_pat_io('comm_pat_glb_to_loc_c',        p%comm_pat_glb_to_loc_c)
      CALL comm_pat_io('comm_pat_glb_to_loc_e',        p%comm_pat_glb_to_loc_e)
      CALL comm_pat_io('comm_pat_loc_to_glb_c_fbk',    p%comm_pat_loc_to_glb_c_fbk)
      CALL comm_pat_io('comm_pat_loc_to_glb_e_fbk',    p%comm_pat_loc_to_glb_e_fbk)
    ENDIF

  END SUBROUTINE patch_io

  !-------------------------------------------------------------------------
  !
  !> Defines interpolation state for NetCDF

  SUBROUTINE def_int_state(k_jg, p)

    INTEGER, INTENT(IN)       :: k_jg  ! domain index
    TYPE(t_patch), INTENT(IN) :: p

    ! local variables
    INTEGER :: idummy

    ! This routine does not involve any global max calls and thus we can just
    ! return if output_defines is not set

    IF(.NOT.output_defines) RETURN

    ! The definitions below are only done once for a patch (i.e. not for local parents)
    ! and thus only if prefix is blank

    IF(prefix==' ') THEN

      ! Output all variables defining interpolation as attributes (unless they are
      ! used as dimensions below)

      CALL nf(nf_put_att_int(ncid, nf_global, 'llsq_high_consv',   nf_int, 1, &
        &  MERGE(1,0,llsq_high_consv)))
      CALL nf(nf_put_att_int(ncid, nf_global, 'rbf_vec_kern_c',   nf_int, 1, rbf_vec_kern_c))
      CALL nf(nf_put_att_int(ncid, nf_global, 'rbf_vec_kern_e',   nf_int, 1, rbf_vec_kern_e))
      CALL nf(nf_put_att_int(ncid, nf_global, 'rbf_vec_kern_v',   nf_int, 1, rbf_vec_kern_v))
      CALL nf(nf_put_att_int(ncid, nf_global, 'i_cori_method',    nf_int, 1, i_cori_method))
      CALL nf(nf_put_att_int(ncid, nf_global, 'nudge_zone_width', nf_int, 1, nudge_zone_width))

      IF(p%id > 0) THEN
        CALL nf(nf_put_att_double(ncid, nf_global, 'rbf_vec_scale_c', nf_double, 1, &
          &  rbf_vec_scale_c(p%id)))
        CALL nf(nf_put_att_double(ncid, nf_global, 'rbf_vec_scale_e', nf_double, 1, &
          &  rbf_vec_scale_e(p%id)))
        CALL nf(nf_put_att_double(ncid, nf_global, 'rbf_vec_scale_v', nf_double, 1, &
          &  rbf_vec_scale_v(p%id)))
      ENDIF
      CALL nf(nf_put_att_double(ncid, nf_global, 'nudge_max_coeff', nf_double, 1, nudge_max_coeff))
      CALL nf(nf_put_att_double(ncid, nf_global, 'nudge_efold_width', nf_double, 1, &
        &  nudge_efold_width))
      CALL nf(nf_put_att_double(ncid, nf_global, 'sick_a',          nf_double, 1, sick_a))

      ! Dimensions for interpolation state

      SELECT CASE (i_cori_method)
      CASE (1,3,4)
        CALL def_dim('nincr', 14, dim_nincr)
      CASE (2)
        CALL def_dim('nincr', 10, dim_nincr)
      END SELECT

      CALL def_dim('lsq_dim_c_lin',    lsq_lin_set%dim_c, dim_lsq_dim_c_lin)
      CALL def_dim('lsq_dim_unk_lin', lsq_lin_set%dim_unk, dim_lsq_dim_unk_lin)
      idummy=(lsq_lin_set%dim_unk*lsq_lin_set%dim_unk - lsq_lin_set%dim_unk)/2
      CALL def_dim('lsq_dim_unk2_lin', idummy,     dim_lsq_dim_unk2_lin)

      CALL def_dim('lsq_dim_c_high',    lsq_high_set%dim_c, dim_lsq_dim_c_high)
      CALL def_dim('lsq_dim_unk_high', lsq_high_set%dim_unk, dim_lsq_dim_unk_high)
      idummy=(lsq_high_set%dim_unk*lsq_high_set%dim_unk - lsq_high_set%dim_unk)/2
      CALL def_dim('lsq_dim_unk2_high', idummy,   dim_lsq_dim_unk2_high)

      CALL def_dim('rbf_c2grad_dim', rbf_c2grad_dim, dim_rbf_c2grad_dim)
      CALL def_dim('rbf_vec_dim_c',  rbf_vec_dim_c,  dim_rbf_vec_dim_c)
      CALL def_dim('rbf_vec_dim_e',  rbf_vec_dim_e,  dim_rbf_vec_dim_e)
      CALL def_dim('rbf_vec_dim_v',  rbf_vec_dim_v,  dim_rbf_vec_dim_v)

    ENDIF

    IF(max_patch_cells <= 0) RETURN

    ! Variables for interpolation state

    CALL def_var('int.c_lin_e',           nf_double, dim_nedges, dim_2) ! nproma,2,nblks_e
    IF (p%cell_type == 3) THEN
    CALL def_var('int.e_bln_c_s',         nf_double, dim_ncells, dim_3) ! nproma,3,nblks_c
    CALL def_var('int.e_bln_c_u',         nf_double, dim_ncells, dim_6) ! nproma,6,nblks_c
    CALL def_var('int.e_bln_c_v',         nf_double, dim_ncells, dim_6) ! nproma,6,nblks_c
    CALL def_var('int.c_bln_avg',         nf_double, dim_ncells, dim_4) ! nproma,4,nblks_c
    CALL def_var('int.e_flx_avg',         nf_double, dim_nedges, dim_5) ! nproma,5,nblks_e
    CALL def_var('int.v_1o2_e',           nf_double, dim_nedges, dim_2) ! nproma,2,nblks_e
    ENDIF
    CALL def_var('int.e_inn_c',           nf_double, dim_ncells, dim_nverts_per_cell) ! nproma,p%cell_type,nblks_c
    IF (p%cell_type == 6) THEN
    CALL def_var('int.e_inn_v',           nf_double, dim_nverts, dim_3) ! nproma,3,nblks_v
    CALL def_var('int.e_aw_c',            nf_double, dim_ncells, dim_6) ! nproma,6,nblks_c
    CALL def_var('int.r_aw_c',            nf_double, dim_ncells, dim_6) ! nproma,6,nblks_c
    CALL def_var('int.e_aw_v',            nf_double, dim_nverts, dim_3) ! nproma,3,nblks_v
    CALL def_var('int.e_1o3_v',           nf_double, dim_nverts, dim_3) ! nproma,3,nblks_v
    CALL def_var('int.tria_aw_rhom',      nf_double, dim_nedges, dim_2) ! nproma,2,nblks_e
    ENDIF
    CALL def_var('int.verts_aw_cells',    nf_double, dim_ncells, dim_nverts_per_cell) ! nproma,p%cell_type,nblks_c
    CALL def_var('int.cells_aw_verts',    nf_double, dim_nverts, dim_nedges_per_vert) ! nproma,9-p%cell_type,nblks_v
    IF( p%cell_type == 6 ) THEN
    CALL def_var('int.tria_north',        nf_double, dim_nverts, dim_3) ! 3,nproma,nblks_v
    CALL def_var('int.tria_east',         nf_double, dim_nverts, dim_3) ! 3,nproma,nblks_v
    CALL def_var('int.hex_north',         nf_double, dim_ncells, dim_6) ! nproma,6,nblks_c
    CALL def_var('int.hex_east',          nf_double, dim_ncells, dim_6) ! nproma,6,nblks_c
    IF(i_cori_method>=3) THEN
    CALL def_var('int.quad_north',        nf_double, dim_nedges, dim_5) ! 5,nproma,nblks_e
    CALL def_var('int.quad_east',         nf_double, dim_nedges, dim_5) ! 5,nproma,nblks_e
    ENDIF
    CALL def_var('int.cno_en',            nf_double, dim_nedges, dim_2) ! nproma,2,nblks_e
    CALL def_var('int.cea_en',            nf_double, dim_nedges, dim_2) ! nproma,2,nblks_e
    CALL def_var('int.heli_coeff',        nf_double, dim_nedges, dim_nincr) ! nincr,nproma,nblks_e
    IF(i_cori_method<3) THEN
    CALL def_var('int.heli_vn_index',     nf_int,    dim_nedges, dim_nincr) ! nincr,nproma,nblks_e
    ENDIF
    ENDIF
    IF (p%cell_type == 3) THEN
    CALL def_var('int.rbf_vec_index_c',   nf_int,    dim_ncells, dim_rbf_vec_dim_c) ! rbf_vec_dim_c,nproma,nblks_c
    CALL def_var('int.rbf_vec_stencil_c', nf_int,    dim_ncells) ! nproma,nblks_c
    CALL def_var('int.rbf_vec_coeff_c',   nf_double, dim_ncells, dim_rbf_vec_dim_c, dim_2) ! rbf_vec_dim_c,2,nproma,nblks_c
    CALL def_var('int.rbf_c2grad_index',  nf_int,    dim_ncells, dim_rbf_c2grad_dim) ! rbf_c2grad_dim,nproma,nblks_c
    CALL def_var('int.rbf_c2grad_coeff',  nf_double, dim_ncells, dim_rbf_c2grad_dim, dim_2) ! rbf_c2grad_dim,2,nproma,nblks_c
    CALL def_var('int.rbf_vec_index_v',   nf_int,    dim_nverts, dim_rbf_vec_dim_v) ! rbf_vec_dim_v,nproma,nblks_v
    CALL def_var('int.rbf_vec_stencil_v', nf_int,    dim_nverts) ! nproma,nblks_v
    CALL def_var('int.rbf_vec_coeff_v',   nf_double, dim_nverts, dim_rbf_vec_dim_v, dim_2) ! rbf_vec_dim_v,2,nproma,nblks_v
    CALL def_var('int.rbf_vec_index_e',   nf_int,    dim_nedges, dim_rbf_vec_dim_e) ! rbf_vec_dim_e,nproma,nblks_e
    CALL def_var('int.rbf_vec_stencil_e', nf_int,    dim_nedges) ! nproma,nblks_e
    CALL def_var('int.rbf_vec_coeff_e',   nf_double, dim_nedges, dim_rbf_vec_dim_e) ! rbf_vec_dim_e,nproma,nblks_e
    ENDIF
    IF( ltransport .OR. iequations == 3) THEN
    CALL def_var('int.pos_on_tplane_e',   nf_double, dim_nedges, dim_8, dim_2) ! nproma,nblks_e,8,2
    CALL def_var('int.tplane_e_dotprod',  nf_double, dim_nedges, dim_4, dim_4) ! nproma,nblks_e,4,4
    CALL def_var('int.lin.lsq_dim_stencil', nf_int,   dim_ncells) ! nproma,nblks_c
    CALL def_var('int.lin.lsq_index_c',     nf_int,   dim_ncells, dim_lsq_dim_c_lin) ! nproma,nblks_c,lsq_dim_c
    CALL def_var('int.lin.lsq_weights_c',   nf_double,dim_ncells, dim_lsq_dim_c_lin) ! nproma,lsq_dim_c,nblks_c
    CALL def_var('int.lin.lsq_qtmat_c',     & ! nproma,lsq_dim_unk,lsq_dim_c,nblks_c
      &          nf_double,dim_ncells, dim_lsq_dim_unk_lin, dim_lsq_dim_c_lin)
    CALL def_var('int.lin.lsq_rmat_rdiag_c',nf_double,dim_ncells, dim_lsq_dim_unk_lin) ! nproma,lsq_dim_unk,nblks_c
    CALL def_var('int.lin.lsq_rmat_utri_c', nf_double,dim_ncells, dim_lsq_dim_unk2_lin) ! nproma,idummy,nblks_c
    CALL def_var('int.lin.lsq_pseudoinv'  , & ! nproma,lsq_dim_unk,lsq_dim_c,nblks_c
      &          nf_double,dim_ncells, dim_lsq_dim_unk_lin, dim_lsq_dim_c_lin)
    CALL def_var('int.lin.lsq_moments',     nf_double,dim_ncells, dim_lsq_dim_unk_lin) ! nproma,nblks_c,lsq_dim_unk
    CALL def_var('int.lin.lsq_moments_hat', & ! nproma,nblks_c,lsq_dim_c,lsq_dim_unk
      &          nf_double,dim_ncells, dim_lsq_dim_c_lin, dim_lsq_dim_unk_lin)

    CALL def_var('int.high.lsq_dim_stencil', nf_int,   dim_ncells) ! nproma,nblks_c
    CALL def_var('int.high.lsq_index_c',     nf_int,   dim_ncells, dim_lsq_dim_c_high) ! nproma,nblks_c,lsq_dim_c
    CALL def_var('int.high.lsq_weights_c',   nf_double,dim_ncells, dim_lsq_dim_c_high) ! nproma,lsq_dim_c,nblks_c
    CALL def_var('int.high.lsq_qtmat_c',     & ! nproma,lsq_dim_unk,lsq_dim_c,nblks_c
      &          nf_double,dim_ncells, dim_lsq_dim_unk_high, dim_lsq_dim_c_high)
    CALL def_var('int.high.lsq_rmat_rdiag_c',nf_double,dim_ncells, dim_lsq_dim_unk_high) ! nproma,lsq_dim_unk,nblks_c
    CALL def_var('int.high.lsq_rmat_utri_c', nf_double,dim_ncells, dim_lsq_dim_unk2_high) ! nproma,idummy,nblks_c
    CALL def_var('int.high.lsq_pseudoinv'  , & ! nproma,lsq_dim_unk,lsq_dim_c,nblks_c
      &          nf_double,dim_ncells, dim_lsq_dim_unk_high, dim_lsq_dim_c_high)
    CALL def_var('int.high.lsq_moments',     nf_double,dim_ncells, dim_lsq_dim_unk_high) ! nproma,nblks_c,lsq_dim_unk
    CALL def_var('int.high.lsq_moments_hat', & ! nproma,nblks_c,lsq_dim_c,lsq_dim_unk
      &          nf_double,dim_ncells, dim_lsq_dim_c_high, dim_lsq_dim_unk_high)
    END IF
    IF (p%cell_type == 3) THEN
    CALL def_var('int.geofac_qdiv',       nf_double, dim_nedges, dim_4) ! nproma,4,nblks_e
    CALL def_var('int.geofac_grdiv',      nf_double, dim_nedges, dim_5) ! nproma,5,nblks_e
    CALL def_var('int.nudgecoeff_c',      nf_double, dim_ncells) ! nproma,nblks_c
    CALL def_var('int.nudgecoeff_e',      nf_double, dim_nedges) ! nproma,nblks_e

    CALL def_var('int.gquad.qpts_tri_l.lon',nf_double, dim_ncells) ! nproma,nblks_c
    CALL def_var('int.gquad.qpts_tri_l.lat',nf_double, dim_ncells) ! nproma,nblks_c
    CALL def_var('int.gquad.qpts_tri_q.lon',nf_double, dim_ncells, dim_3) ! nproma,nblks_c,3
    CALL def_var('int.gquad.qpts_tri_q.lat',nf_double, dim_ncells, dim_3) ! nproma,nblks_c,3
    CALL def_var('int.gquad.qpts_tri_c.lon',nf_double, dim_ncells, dim_4) ! nproma,nblks_c,4
    CALL def_var('int.gquad.qpts_tri_c.lat',nf_double, dim_ncells, dim_4) ! nproma,nblks_c,4
    CALL def_var('int.gquad.weights_tri_q',nf_double, dim_3) ! 3
    CALL def_var('int.gquad.weights_tri_c',nf_double, dim_4) ! 4
    ENDIF
    CALL def_var('int.geofac_div',        nf_double, dim_ncells, dim_nverts_per_cell) ! nproma,p%cell_type,nblks_c
    CALL def_var('int.geofac_rot',        nf_double, dim_nverts, dim_nedges_per_vert) ! nproma,9-p%cell_type,nblks_v
    CALL def_var('int.geofac_n2s',        nf_double, dim_ncells, dim_nverts_per_cell_p1) ! nproma,p%cell_type+1,nblks_c
    CALL def_var('int.geofac_grg',        nf_double, dim_ncells, dim_nverts_per_cell_p1, dim_2) ! nproma,p%cell_type+1,nblks_c,2
    CALL def_var('int.cart_edge_coord',   nf_double, dim_nedges, dim_3) ! nproma,nblks_e,3
    CALL def_var('int.cart_cell_coord',   nf_double, dim_ncells, dim_3) ! nproma,nblks_c,3
    CALL def_var('int.primal_normal_ec',  nf_double, dim_ncells, dim_nverts_per_cell, dim_2) ! nproma,nblks_c,p%cell_type,2
    CALL def_var('int.edge_cell_length',  nf_double, dim_ncells, dim_nverts_per_cell) ! nproma,nblks_c,p%cell_type
    IF (p%cell_type == 6) THEN
    CALL def_var('int.dir_gradh_index1',  nf_int,    dim_nedges, dim_6) ! 6,nproma,nblks_e
    CALL def_var('int.dir_gradh_index2',  nf_int,    dim_nedges, dim_6) ! 6,nproma,nblks_e
    CALL def_var('int.dir_gradhux_c1',    nf_double, dim_nedges, dim_6) ! 6,nproma,nblks_e
    CALL def_var('int.dir_gradhux_c2',    nf_double, dim_nedges, dim_6) ! 6,nproma,nblks_e
    CALL def_var('int.strain_def_c1',    nf_double, dim_nedges, dim_6) ! 6,nproma,nblks_e
    CALL def_var('int.strain_def_c2',    nf_double, dim_nedges, dim_6) ! 6,nproma,nblks_e
    CALL def_var('int.dir_gradt_index1',  nf_int,    dim_nedges, dim_9) ! 9,nproma,nblks_e
    CALL def_var('int.dir_gradt_index2',  nf_int,    dim_nedges, dim_9) ! 9,nproma,nblks_e
    CALL def_var('int.dir_gradtxy_v1',    nf_double, dim_nedges, dim_9) ! 9,nproma,nblks_e
    CALL def_var('int.dir_gradtxy_v2',    nf_double, dim_nedges, dim_9) ! 9,nproma,nblks_e
    CALL def_var('int.dir_gradtyx_v1',    nf_double, dim_nedges, dim_9) ! 9,nproma,nblks_e
    CALL def_var('int.dir_gradtyx_v2',    nf_double, dim_nedges, dim_9) ! 9,nproma,nblks_e
    CALL def_var('int.shear_def_v1',    nf_double, dim_nedges, dim_9) ! 9,nproma,nblks_e
    CALL def_var('int.shear_def_v2',    nf_double, dim_nedges, dim_9) ! 9,nproma,nblks_e
    ENDIF

  END SUBROUTINE def_int_state

  !-------------------------------------------------------------------------
  !
  !> Defines lonlat state for NetCDF

  SUBROUTINE def_lonlat_state(k_jg)

    INTEGER, INTENT(IN)       :: k_jg  ! domain index

    TYPE (t_lon_lat_grid), POINTER :: grid ! lon-lat grid
    INTEGER :: dim_nlonlat

    grid => lonlat_intp_config(k_jg)%lonlat_grid
    ! total number of lon-lat grid points
    CALL def_dim('nlonlat', grid%total_dim, dim_nlonlat)
    ! description of lon-lat grid
    CALL nf(nf_put_att_double(ncid, nf_global, 'int.lonlat.delta',     &
      &     nf_double, 2, grid%delta(1:2)))
    CALL nf(nf_put_att_double(ncid, nf_global, 'int.lonlat.sw_corner', &
      &     nf_double, 2, grid%start_corner(1:2)))
    CALL nf(nf_put_att_double(ncid, nf_global, 'int.lonlat.poleN',     &
      &     nf_double, 2, grid%poleN(1:2)))
    CALL nf(nf_put_att_int(ncid, nf_global, 'int.lonlat.dimen',        &
      &     nf_int, 2, grid%dimen(1:2)))

    ! fields for indices and coefficients:
    ! rbf_vec_dim_c,2,nproma,nblks_lonlat
    CALL def_var('int.lonlat.rbf_vec_coeff',    &
      &          nf_double, dim_nlonlat, dim_rbf_vec_dim_c, dim_2, add_unlim=.FALSE.)
    ! rbf_c2grad_dim,2,nproma,nblks_lonlat
    CALL def_var('int.lonlat.rbf_c2grad_coeff', &
      &          nf_double, dim_nlonlat, dim_rbf_c2grad_dim, dim_2, add_unlim=.FALSE.)
    ! rbf_vec_dim_c,nproma,nblks_lonlat
    CALL def_var('int.lonlat.rbf_vec_index',    &
      &          nf_int,    dim_nlonlat, dim_rbf_vec_dim_c, add_unlim=.FALSE.)
    ! nproma,nblks_lonlat
    CALL def_var('int.lonlat.rbf_vec_stencil',  &
      &          nf_int,    dim_nlonlat, add_unlim=.FALSE.)
    ! rbf_c2grad_dim,nproma,nblks_lonlat
    CALL def_var('int.lonlat.rbf_c2grad_index', &
      &          nf_int,    dim_nlonlat, dim_rbf_c2grad_dim, add_unlim=.FALSE.)
    ! 2,nproma,nblks_lonlat
    CALL def_var('int.lonlat.rdist',            &
      &          nf_double, dim_nlonlat, dim_2, add_unlim=.FALSE.)
    ! 2,nproma,nblks_lonlat
    CALL def_var('int.lonlat.tri_idx',          &
      &          nf_int,    dim_nlonlat, dim_2, add_unlim=.FALSE.)

  END SUBROUTINE def_lonlat_state

  !-------------------------------------------------------------------------
  !
  !> Dumps/restores interpolation state data to/from NetCDF file.

  SUBROUTINE int_state_io(jg, p, pi)

    INTEGER, INTENT(IN) :: jg
    TYPE(t_patch), INTENT(IN) :: p
    TYPE(t_int_state), INTENT(INOUT) :: pi

    IF(p%n_patch_cells <= 0) RETURN

    CALL store_proc_dependent_dimensions(p)

    CALL bvar_io(1,3,'int.c_lin_e',           pi%c_lin_e          ) ! nproma,2,nblks_e
    IF (p%cell_type == 3) THEN
    CALL bvar_io(1,3,'int.e_bln_c_s',         pi%e_bln_c_s        ) ! nproma,3,nblks_c
    CALL bvar_io(1,3,'int.e_bln_c_u',         pi%e_bln_c_u        ) ! nproma,6,nblks_c
    CALL bvar_io(1,3,'int.e_bln_c_v',         pi%e_bln_c_v        ) ! nproma,6,nblks_c
    CALL bvar_io(1,3,'int.c_bln_avg',         pi%c_bln_avg        ) ! nproma,4,nblks_c
    CALL bvar_io(1,3,'int.e_flx_avg',         pi%e_flx_avg        ) ! nproma,5,nblks_e
    CALL bvar_io(1,3,'int.v_1o2_e',           pi%v_1o2_e          ) ! nproma,2,nblks_e
    ENDIF
    CALL bvar_io(1,3,'int.e_inn_c',           pi%e_inn_c          ) ! nproma,p%cell_type,nblks_c
    IF (p%cell_type == 6) THEN
    CALL bvar_io(1,3,'int.e_inn_v',           pi%e_inn_v          ) ! nproma,3,nblks_v
    CALL bvar_io(1,3,'int.e_aw_c',            pi%e_aw_c           ) ! nproma,6,nblks_c
    CALL bvar_io(1,3,'int.r_aw_c',            pi%r_aw_c           ) ! nproma,6,nblks_c
    CALL bvar_io(1,3,'int.e_aw_v',            pi%e_aw_v           ) ! nproma,3,nblks_v
    CALL bvar_io(1,3,'int.e_1o3_v',           pi%e_1o3_v          ) ! nproma,3,nblks_v
    CALL bvar_io(1,3,'int.tria_aw_rhom',      pi%tria_aw_rhom     ) ! nproma,2,nblks_e
    ENDIF
    CALL bvar_io(1,3,'int.verts_aw_cells',    pi%verts_aw_cells   ) ! nproma,p%cell_type,nblks_c
    CALL bvar_io(1,3,'int.cells_aw_verts',    pi%cells_aw_verts   ) ! nproma,9-p%cell_type,nblks_v
    IF( p%cell_type == 6 ) THEN
    CALL bvar_io(2,3,'int.tria_north',        pi%tria_north       ) ! 3,nproma,nblks_v
    CALL bvar_io(2,3,'int.tria_east',         pi%tria_east        ) ! 3,nproma,nblks_v
    CALL bvar_io(1,3,'int.hex_north',         pi%hex_north        ) ! nproma,6,nblks_c
    CALL bvar_io(1,3,'int.hex_east',          pi%hex_east         ) ! nproma,6,nblks_c
    IF (i_cori_method>=3) THEN
    CALL bvar_io(2,3,'int.quad_north',        pi%quad_north       ) ! 5,nproma,nblks_e
    CALL bvar_io(2,3,'int.quad_east',         pi%quad_east        ) ! 5,nproma,nblks_e
    ENDIF
    CALL bvar_io(1,3,'int.cno_en',            pi%cno_en           ) ! nproma,2,nblks_e
    CALL bvar_io(1,3,'int.cea_en',            pi%cea_en           ) ! nproma,2,nblks_e
    CALL bvar_io(2,3,'int.heli_coeff',        pi%heli_coeff       ) ! nincr,nproma,nblks_e
    IF (i_cori_method < 3) THEN
    CALL bidx_io(2,3,'int.heli_vn_index',     pi%heli_vn_idx, pi%heli_vn_blk) ! nincr,nproma,nblks_e
    ENDIF
    ENDIF
    IF (p%cell_type == 3) THEN
    CALL bidx_io(2,3,'int.rbf_vec_index_c',   pi%rbf_vec_idx_c, pi%rbf_vec_blk_c) ! rbf_vec_dim_c,nproma,nblks_c
    CALL bvar_io(1,2,'int.rbf_vec_stencil_c', pi%rbf_vec_stencil_c) ! nproma,nblks_c
    CALL bvar_io(3,4,'int.rbf_vec_coeff_c',   pi%rbf_vec_coeff_c  ) ! rbf_vec_dim_c,2,nproma,nblks_c
    CALL bidx_io(2,3,'int.rbf_c2grad_index',  pi%rbf_c2grad_idx, pi%rbf_c2grad_blk) ! rbf_c2grad_dim,nproma,nblks_c
    CALL bvar_io(3,4,'int.rbf_c2grad_coeff',  pi%rbf_c2grad_coeff ) ! rbf_c2grad_dim,2,nproma,nblks_c
    CALL bidx_io(2,3,'int.rbf_vec_index_v',   pi%rbf_vec_idx_v, pi%rbf_vec_blk_v) ! rbf_vec_dim_v,nproma,nblks_v
    CALL bvar_io(1,2,'int.rbf_vec_stencil_v', pi%rbf_vec_stencil_v) ! nproma,nblks_v
    CALL bvar_io(3,4,'int.rbf_vec_coeff_v',   pi%rbf_vec_coeff_v  ) ! rbf_vec_dim_v,2,nproma,nblks_v
    CALL bidx_io(2,3,'int.rbf_vec_index_e',   pi%rbf_vec_idx_e, pi%rbf_vec_blk_e) ! rbf_vec_dim_e,nproma,nblks_e
    CALL bvar_io(1,2,'int.rbf_vec_stencil_e', pi%rbf_vec_stencil_e) ! nproma,nblks_e
    CALL bvar_io(2,3,'int.rbf_vec_coeff_e',   pi%rbf_vec_coeff_e  ) ! rbf_vec_dim_e,nproma,nblks_e
    ENDIF
    IF( ltransport .OR. iequations == 3) THEN
    CALL bvar_io(1,2,'int.pos_on_tplane_e',   pi%pos_on_tplane_e  ) ! nproma,nblks_e,8,2
    CALL bvar_io(1,2,'int.tplane_e_dotprod',  pi%tplane_e_dotprod ) ! nproma,nblks_e,4,4
    CALL bvar_io(1,2,'int.lin.lsq_dim_stencil', pi%lsq_lin%lsq_dim_stencil ) ! nproma,nblks_c
    CALL bidx_io(1,2,'int.lin.lsq_index_c',     pi%lsq_lin%lsq_idx_c, pi%lsq_lin%lsq_blk_c) ! nproma,nblks_c,lsq_dim_c
    CALL bvar_io(1,3,'int.lin.lsq_weights_c',   pi%lsq_lin%lsq_weights_c   ) ! nproma,lsq_dim_c,nblks_c
    CALL bvar_io(1,4,'int.lin.lsq_qtmat_c',     pi%lsq_lin%lsq_qtmat_c     ) ! nproma,lsq_dim_unk,lsq_dim_c,nblks_c
    CALL bvar_io(1,3,'int.lin.lsq_rmat_rdiag_c',pi%lsq_lin%lsq_rmat_rdiag_c) ! nproma,lsq_dim_unk,nblks_c
    CALL bvar_io(1,3,'int.lin.lsq_rmat_utri_c', pi%lsq_lin%lsq_rmat_utri_c ) ! nproma,idummy,nblks_c
    CALL bvar_io(1,4,'int.lin.lsq_pseudoinv',   pi%lsq_lin%lsq_pseudoinv   ) ! nproma,lsq_dim_unk,lsq_dim_c,nblks_c
    CALL bvar_io(1,2,'int.lin.lsq_moments',     pi%lsq_lin%lsq_moments     ) ! nproma,nblks_c,lsq_dim_unk
    CALL bvar_io(1,2,'int.lin.lsq_moments_hat', pi%lsq_lin%lsq_moments_hat ) ! nproma,nblks_c,lsq_dim_c,lsq_dim_unk

    CALL bvar_io(1,2,'int.high.lsq_dim_stencil', pi%lsq_high%lsq_dim_stencil ) ! nproma,nblks_c
    CALL bidx_io(1,2,'int.high.lsq_index_c',     pi%lsq_high%lsq_idx_c, pi%lsq_high%lsq_blk_c) ! nproma,nblks_c,lsq_dim_c
    CALL bvar_io(1,3,'int.high.lsq_weights_c',   pi%lsq_high%lsq_weights_c   ) ! nproma,lsq_dim_c,nblks_c
    CALL bvar_io(1,4,'int.high.lsq_qtmat_c',     pi%lsq_high%lsq_qtmat_c     ) ! nproma,lsq_dim_unk,lsq_dim_c,nblks_c
    CALL bvar_io(1,3,'int.high.lsq_rmat_rdiag_c',pi%lsq_high%lsq_rmat_rdiag_c) ! nproma,lsq_dim_unk,nblks_c
    CALL bvar_io(1,3,'int.high.lsq_rmat_utri_c', pi%lsq_high%lsq_rmat_utri_c ) ! nproma,idummy,nblks_c
    CALL bvar_io(1,4,'int.high.lsq_pseudoinv',   pi%lsq_high%lsq_pseudoinv   ) ! nproma,lsq_dim_unk,lsq_dim_c,nblks_c
    CALL bvar_io(1,2,'int.high.lsq_moments',     pi%lsq_high%lsq_moments     ) ! nproma,nblks_c,lsq_dim_unk
    CALL bvar_io(1,2,'int.high.lsq_moments_hat', pi%lsq_high%lsq_moments_hat ) ! nproma,nblks_c,lsq_dim_c,lsq_dim_unk
    END IF
    IF (p%cell_type == 3) THEN
    CALL bvar_io(1,3,'int.geofac_qdiv',       pi%geofac_qdiv     ) ! nproma,4,nblks_e
    CALL bvar_io(1,3,'int.geofac_grdiv',      pi%geofac_grdiv    ) ! nproma,5,nblks_e
    CALL bvar_io(1,2,'int.nudgecoeff_c',      pi%nudgecoeff_c    ) ! nproma,nblks_c
    CALL bvar_io(1,2,'int.nudgecoeff_e',      pi%nudgecoeff_e    ) ! nproma,nblks_e

    CALL bvar_io(1,2,'int.gquad.qpts_tri_l.lon',    pi%gquad%qpts_tri_l(:,:)%lon   ) ! nproma,nblks_c
    CALL bvar_io(1,2,'int.gquad.qpts_tri_l.lat',    pi%gquad%qpts_tri_l(:,:)%lat   ) ! nproma,nblks_c
    CALL bvar_io(1,2,'int.gquad.qpts_tri_q.lon',    pi%gquad%qpts_tri_q(:,:,:)%lon ) ! nproma,nblks_c,3
    CALL bvar_io(1,2,'int.gquad.qpts_tri_q.lat',    pi%gquad%qpts_tri_q(:,:,:)%lat ) ! nproma,nblks_c,3
    CALL bvar_io(1,2,'int.gquad.qpts_tri_c.lon',    pi%gquad%qpts_tri_c(:,:,:)%lon ) ! nproma,nblks_c,4
    CALL bvar_io(1,2,'int.gquad.qpts_tri_c.lat',    pi%gquad%qpts_tri_c(:,:,:)%lat ) ! nproma,nblks_c,4
    CALL uvar_io(    'int.gquad.weights_tri_q', pi%gquad%weights_tri_q  ) ! 3
    CALL uvar_io(    'int.gquad.weights_tri_c', pi%gquad%weights_tri_c  ) ! 4
    ENDIF
    CALL bvar_io(1,3,'int.geofac_div',        pi%geofac_div      ) ! nproma,p%cell_type,nblks_c
    CALL bvar_io(1,3,'int.geofac_rot',        pi%geofac_rot      ) ! nproma,9-p%cell_type,nblks_v
    CALL bvar_io(1,3,'int.geofac_n2s',        pi%geofac_n2s      ) ! nproma,p%cell_type+1,nblks_c
    CALL bvar_io(1,3,'int.geofac_grg',        pi%geofac_grg      ) ! nproma,p%cell_type+1,nblks_c,2
    CALL bvar_io(1,2,'int.cart_edge_coord',   pi%cart_edge_coord ) ! nproma,nblks_e,3
    CALL bvar_io(1,2,'int.cart_cell_coord',   pi%cart_cell_coord ) ! nproma,nblks_c,3
    CALL bvar_io(1,2,'int.primal_normal_ec',  pi%primal_normal_ec) ! nproma,nblks_c,p%cell_type,2
    CALL bvar_io(1,2,'int.edge_cell_length',  pi%edge_cell_length) ! nproma,nblks_c,p%cell_type
    IF (p%cell_type == 6) THEN
    CALL bidx_io(2,3,'int.dir_gradh_index1',  pi%dir_gradh_i1, pi%dir_gradh_b1) ! 6,nproma,nblks_e
    CALL bidx_io(2,3,'int.dir_gradh_index2',  pi%dir_gradh_i2, pi%dir_gradh_b2) ! 6,nproma,nblks_e
    CALL bvar_io(2,3,'int.dir_gradhux_c1',    pi%dir_gradhux_c1  ) ! 6,nproma,nblks_e
    CALL bvar_io(2,3,'int.dir_gradhux_c2',    pi%dir_gradhux_c2  ) ! 6,nproma,nblks_e
    CALL bvar_io(2,3,'int.strain_def_c1',    pi%strain_def_c1  ) ! 6,nproma,nblks_e
    CALL bvar_io(2,3,'int.strain_def_c2',    pi%strain_def_c2  ) ! 6,nproma,nblks_e
    CALL bidx_io(2,3,'int.dir_gradt_index1',  pi%dir_gradt_i1, pi%dir_gradt_b1) ! 9,nproma,nblks_e
    CALL bidx_io(2,3,'int.dir_gradt_index2',  pi%dir_gradt_i2, pi%dir_gradt_b2) ! 9,nproma,nblks_e
    CALL bvar_io(2,3,'int.dir_gradtxy_v1',    pi%dir_gradtxy_v1  ) ! 9,nproma,nblks_e
    CALL bvar_io(2,3,'int.dir_gradtxy_v2',    pi%dir_gradtxy_v2  ) ! 9,nproma,nblks_e
    CALL bvar_io(2,3,'int.dir_gradtyx_v1',    pi%dir_gradtyx_v1  ) ! 9,nproma,nblks_e
    CALL bvar_io(2,3,'int.dir_gradtyx_v2',    pi%dir_gradtyx_v2  ) ! 9,nproma,nblks_e
    CALL bvar_io(2,3,'int.shear_def_v1',      pi%shear_def_v1  ) ! 9,nproma,nblks_e
    CALL bvar_io(2,3,'int.shear_def_v2',      pi%shear_def_v2  ) ! 9,nproma,nblks_e
    ENDIF

  END SUBROUTINE int_state_io

  !-------------------------------------------------------------------------
  !
  !> Dumps/restores lonlat state data to/from NetCDF file.

  SUBROUTINE lonlat_state_io(pi_lonlat)

    TYPE (t_lon_lat_intp), INTENT(INOUT) :: pi_lonlat

    ! rbf_vec_dim_c,2,nproma,nblks_lonlat
    CALL bvar_io(3,4,'int.lonlat.rbf_vec_coeff',    &
      &          pi_lonlat%rbf_vec_coeff  )
    ! rbf_c2grad_dim,2,nproma,nblks_lonlat
    CALL bvar_io(3,4,'int.lonlat.rbf_c2grad_coeff', &
      &          pi_lonlat%rbf_c2grad_coeff )
    ! rbf_vec_dim_c,nproma,nblks_lonlat
    CALL bidx_io(2,3,'int.lonlat.rbf_vec_index',    &
      &          pi_lonlat%rbf_vec_idx, pi_lonlat%rbf_vec_blk)
    ! nproma,nblks_lonlat
    CALL bvar_io(1,2,'int.lonlat.rbf_vec_stencil',  &
      &          pi_lonlat%rbf_vec_stencil)
    ! rbf_c2grad_dim,nproma,nblks_lonlat
    CALL bidx_io(2,3,'int.lonlat.rbf_c2grad_index', &
      &          pi_lonlat%rbf_c2grad_idx, pi_lonlat%rbf_c2grad_blk)
    ! 2,nproma,nblks_lonlat
    CALL bvar_io(2,3,'int.lonlat.rdist',            &
      &          pi_lonlat%rdist)
    ! 2,nproma,nblks_lonlat
    CALL bvar_io(2,3,'int.lonlat.tri_idx',          &
      &          pi_lonlat%tri_idx)

  END SUBROUTINE lonlat_state_io

  !-------------------------------------------------------------------------
  !
  !> Defines GRF state for NetCDF

  SUBROUTINE def_grf_state(p)

    TYPE(t_patch), INTENT(IN) :: p

    CHARACTER(LEN=80) ccd
    INTEGER :: jcd, js_e, je_e, js_c, je_c

    INTEGER :: max_grf_edges, max_grf_cells
    INTEGER :: dim_grf_edges, dim_grf_cells, dim_n_childdom

    ! The definitions below are only done once for a patch (i.e. not for local parents)
    ! and thus only if prefix is blank

    IF(output_defines .AND. prefix==' ') THEN
      ! Output all variables defining gridref state as attributes (unless they are
      ! used as dimensions below)
      CALL nf(nf_put_att_int(ncid, nf_global, 'rbf_vec_kern_grf_e', nf_int, 1, rbf_vec_kern_grf_e))
      CALL nf(nf_put_att_int(ncid, nf_global, 'grf_intmethod_c',    nf_int, 1, grf_intmethod_c))
      CALL nf(nf_put_att_int(ncid, nf_global, 'grf_intmethod_ct',   nf_int, 1, grf_intmethod_ct))
      CALL nf(nf_put_att_int(ncid, nf_global, 'grf_intmethod_e',    nf_int, 1, grf_intmethod_e))
      CALL nf(nf_put_att_int(ncid, nf_global, 'grf_velfbk',         nf_int, 1, grf_velfbk))
      CALL nf(nf_put_att_int(ncid, nf_global, 'grf_scalfbk',        nf_int, 1, grf_scalfbk))
      CALL nf(nf_put_att_int(ncid, nf_global, 'grf_tracfbk',        nf_int, 1, grf_tracfbk))

      IF(p%id > 0) &
        & CALL nf(nf_put_att_double(ncid, nf_global, 'rbf_scale_grf_e', nf_double, 1, &
        &         rbf_scale_grf_e(p%id)))
      CALL nf(nf_put_att_double(ncid, nf_global, 'grf_idw_exp_e12', nf_double, 1, grf_idw_exp_e12))
      CALL nf(nf_put_att_double(ncid, nf_global, 'grf_idw_exp_e34', nf_double, 1, grf_idw_exp_e34))
      CALL nf(nf_put_att_double(ncid, nf_global, 'denom_diffu_v',   nf_double, 1, denom_diffu_v))
      CALL nf(nf_put_att_double(ncid, nf_global, 'denom_diffu_t',   nf_double, 1, denom_diffu_t))

      ! Dimensions and variables for grf state

      CALL def_dim('grf_vec_dim_1', grf_vec_dim_1, dim_grf_vec_dim_1)
      CALL def_dim('grf_vec_dim_2', grf_vec_dim_2, dim_grf_vec_dim_2)
    ENDIF

    IF(p%n_childdom>0 .AND. output_defines) THEN
      CALL def_dim('n_childdom', p%n_childdom, dim_n_childdom)
      CALL def_var('grf.fbk_dom_area', nf_double, dim_n_childdom)
    ENDIF

    DO jcd = 1, p%n_childdom

      js_e = idx_1d(p%edges%start_idx(grf_bdyintp_start_e,jcd), &
        &           p%edges%start_blk(grf_bdyintp_start_e,jcd))
      je_e = idx_1d(p%edges%end_idx(min_rledge_int,jcd),p%edges%end_blk(min_rledge_int,jcd))
      js_c = idx_1d(p%cells%start_idx(grf_bdyintp_start_c,jcd), &
        &           p%cells%start_blk(grf_bdyintp_start_c,jcd))
      je_c = idx_1d(p%cells%end_idx(min_rlcell_int,jcd),p%cells%end_blk(min_rlcell_int,jcd))

      WRITE(ccd,'("ch_dom.",i0)') jcd

      IF(use_one_file) THEN
        max_grf_edges = p_max(je_e-js_e+1, comm=p_comm_work)
      ELSE
        max_grf_edges = je_e-js_e+1
      ENDIF

      IF(max_grf_edges>0 .AND. output_defines) THEN
        CALL def_dim(TRIM(ccd)//'.max_grf_edges', max_grf_edges, dim_grf_edges)

        CALL def_var(TRIM(ccd)//'.grf_vec_index_1a',   nf_int,    dim_grf_edges, dim_grf_vec_dim_1)
        CALL def_var(TRIM(ccd)//'.grf_vec_index_1b',   nf_int,    dim_grf_edges, dim_grf_vec_dim_1)
        CALL def_var(TRIM(ccd)//'.grf_vec_index_2a',   nf_int,    dim_grf_edges, dim_grf_vec_dim_2)
        CALL def_var(TRIM(ccd)//'.grf_vec_index_2b',   nf_int,    dim_grf_edges, dim_grf_vec_dim_2)

        CALL def_var(TRIM(ccd)//'.grf_vec_stencil_1a', nf_int,    dim_grf_edges)
        CALL def_var(TRIM(ccd)//'.grf_vec_stencil_1b', nf_int,    dim_grf_edges)
        CALL def_var(TRIM(ccd)//'.grf_vec_stencil_2a', nf_int,    dim_grf_edges)
        CALL def_var(TRIM(ccd)//'.grf_vec_stencil_2b', nf_int,    dim_grf_edges)

        CALL def_var(TRIM(ccd)//'.grf_vec_coeff_1a',   nf_double, dim_grf_edges, dim_grf_vec_dim_1)
        CALL def_var(TRIM(ccd)//'.grf_vec_coeff_1b',   nf_double, dim_grf_edges, dim_grf_vec_dim_1)
        CALL def_var(TRIM(ccd)//'.grf_vec_coeff_2a',   nf_double, dim_grf_edges, dim_grf_vec_dim_2)
        CALL def_var(TRIM(ccd)//'.grf_vec_coeff_2b',   nf_double, dim_grf_edges, dim_grf_vec_dim_2)

        CALL def_var(TRIM(ccd)//'.grf_dist_pe2ce',     nf_double, dim_grf_edges, dim_2)
      ENDIF

      IF(use_one_file) THEN
        max_grf_cells = p_max(je_c-js_c+1, comm=p_comm_work)
      ELSE
        max_grf_cells = je_c-js_c+1
      ENDIF

      IF(max_grf_cells>0 .AND. output_defines) THEN
        CALL def_dim(TRIM(ccd)//'.max_grf_cells', max_grf_cells, dim_grf_cells)
        CALL def_var(TRIM(ccd)//'.grf_dist_pc2cc', nf_double, dim_grf_cells, dim_4, dim_2)
      ENDIF

    ENDDO

    IF(max_patch_cells>0 .AND. output_defines) THEN
      CALL def_var('grf.fbk_wgt_c',  nf_double, dim_ncells, dim_4)
      CALL def_var('grf.fbk_wgt_ct', nf_double, dim_ncells, dim_4)
      CALL def_var('grf.fbk_wgt_e',  nf_double, dim_nedges, dim_6)
    ENDIF

  END SUBROUTINE def_grf_state

  !-------------------------------------------------------------------------
  !
  !> Dumps/restores GRF state data to/from NetCDF file.

  SUBROUTINE grf_state_io(p, pg)

    TYPE(t_patch), INTENT(IN) :: p
    TYPE(t_gridref_state), INTENT(INOUT) :: pg

    CHARACTER(LEN=80) ccd
    INTEGER :: jcd, js_e, je_e, js_c, je_c, is_e, is_c


    CALL store_proc_dependent_dimensions(p)

    IF(p%n_childdom>0) THEN
      CALL uvar_io(    'grf.fbk_dom_area',   pg%fbk_dom_area)
    ENDIF

    IF(p%n_patch_cells <= 0) RETURN

    DO jcd = 1, p%n_childdom

      is_e = p%edges%start_idx(grf_bdyintp_start_e,jcd)
      is_c = p%cells%start_idx(grf_bdyintp_start_c,jcd)

      js_e = idx_1d(p%edges%start_idx(grf_bdyintp_start_e,jcd),&
        &           p%edges%start_blk(grf_bdyintp_start_e,jcd))
      je_e = idx_1d(p%edges%end_idx(min_rledge_int,jcd),p%edges%end_blk(min_rledge_int,jcd))
      js_c = idx_1d(p%cells%start_idx(grf_bdyintp_start_c,jcd), &
        &           p%cells%start_blk(grf_bdyintp_start_c,jcd))
      je_c = idx_1d(p%cells%end_idx(min_rlcell_int,jcd),p%cells%end_blk(min_rlcell_int,jcd))

      ! Cycle if there is nothing to do.
      ! This also avoids inquiring possibly nonexisting dimensions below.
      ! IF(je_e < js_e .OR. je_c < js_c) CYCLE

      WRITE(ccd,'("ch_dom.",i0)') jcd

      ! Get the dimension IDs for the max_grf_... dimensions and set the actual values
      ! on the current processor for these dimensions
      ! This is only necessary if use_one_file is set
      IF(use_one_file) THEN
         IF (je_c >= js_c) &
              CALL nf(nf_inq_dimid(ncid,TRIM(prefix)//TRIM(ccd)//'.max_grf_cells',&
                & dim_max_grf_cells))
         IF (je_e >= js_e) &
              CALL nf(nf_inq_dimid(ncid,TRIM(prefix)//TRIM(ccd)//'.max_grf_edges',&
                & dim_max_grf_edges))
      ENDIF

      n_grf_cells = je_c-js_c+1
      n_grf_edges = je_e-js_e+1

      IF (je_e >= js_e) THEN

         CALL bidx_io(1,3,TRIM(ccd)//'.grf_vec_index_1a',   &
              &          pg%p_dom(jcd)%grf_vec_ind_1a, pg%p_dom(jcd)%grf_vec_blk_1a, is_e)
         CALL bidx_io(1,3,TRIM(ccd)//'.grf_vec_index_1b',   &
              &          pg%p_dom(jcd)%grf_vec_ind_1b, pg%p_dom(jcd)%grf_vec_blk_1b, is_e)
         CALL bidx_io(1,3,TRIM(ccd)//'.grf_vec_index_2a',   &
              &          pg%p_dom(jcd)%grf_vec_ind_2a, pg%p_dom(jcd)%grf_vec_blk_2a, is_e)
         CALL bidx_io(1,3,TRIM(ccd)//'.grf_vec_index_2b',   &
              &          pg%p_dom(jcd)%grf_vec_ind_2b, pg%p_dom(jcd)%grf_vec_blk_2b, is_e)
         
         CALL bvar_io(1,2,TRIM(ccd)//'.grf_vec_stencil_1a', &
              &          pg%p_dom(jcd)%grf_vec_stencil_1a, is_e) ! nproma_grf, isb_e:ieb_e
         CALL bvar_io(1,2,TRIM(ccd)//'.grf_vec_stencil_1b', &
              &          pg%p_dom(jcd)%grf_vec_stencil_1b, is_e) ! nproma_grf, isb_e:ieb_e
         CALL bvar_io(1,2,TRIM(ccd)//'.grf_vec_stencil_2a', &
              &          pg%p_dom(jcd)%grf_vec_stencil_2a, is_e) ! nproma_grf, isb_e:ieb_e
         CALL bvar_io(1,2,TRIM(ccd)//'.grf_vec_stencil_2b', &
              &          pg%p_dom(jcd)%grf_vec_stencil_2b, is_e) ! nproma_grf, isb_e:ieb_e

         CALL bvar_io(2,3,TRIM(ccd)//'.grf_vec_coeff_1a',   &
              &          pg%p_dom(jcd)%grf_vec_coeff_1a,   is_e) ! grf_vec_dim_1,nproma_grf,isb_e:ieb_e
         CALL bvar_io(2,3,TRIM(ccd)//'.grf_vec_coeff_1b',   &
              &          pg%p_dom(jcd)%grf_vec_coeff_1b,   is_e) ! grf_vec_dim_1,nproma_grf,isb_e:ieb_e
         CALL bvar_io(2,3,TRIM(ccd)//'.grf_vec_coeff_2a',   &
              &          pg%p_dom(jcd)%grf_vec_coeff_2a,   is_e) ! grf_vec_dim_2,nproma_grf,isb_e:ieb_e
         CALL bvar_io(2,3,TRIM(ccd)//'.grf_vec_coeff_2b',   &
              &          pg%p_dom(jcd)%grf_vec_coeff_2b,   is_e) ! grf_vec_dim_2,nproma_grf,isb_e:ieb_e

         CALL bvar_io(1,3,TRIM(ccd)//'.grf_dist_pe2ce',   &
              &          pg%p_dom(jcd)%grf_dist_pe2ce,     is_e) ! nproma_grf, 2, isb_e:ieb_e

      END IF

      IF (je_c >= js_c) THEN

         CALL bvar_io(1,4,TRIM(ccd)//'.grf_dist_pc2cc',   &
              &          pg%p_dom(jcd)%grf_dist_pc2cc,     is_c) ! nproma_grf, 4, 2, isb_c:ieb_c

      END IF

    ENDDO

    CALL bvar_io(1,2,'grf.fbk_wgt_c',  pg%fbk_wgt_c ) ! nproma,nblks_c,4
    CALL bvar_io(1,2,'grf.fbk_wgt_ct', pg%fbk_wgt_ct) ! nproma,nblks_c,4
    CALL bvar_io(1,2,'grf.fbk_wgt_e',  pg%fbk_wgt_e ) ! nproma,nblks_c,6

  END SUBROUTINE grf_state_io
  !
  !-------------------------------------------------------------------------
  ! Checking attributes and dimensions
  !-------------------------------------------------------------------------
  !
  !> check_att_int
  !! Checks if attribute "att_name" in NetCDF file has given value

  SUBROUTINE check_att_int(att_name, att_value)

    CHARACTER(LEN=*), INTENT(IN) :: att_name
    INTEGER, INTENT(IN) :: att_value

    INTEGER :: att_len, stat, val

    stat = nf_inq_attlen (ncid, nf_global, att_name, att_len)
    IF(stat /= nf_noerr) &
      & CALL finish(modname,'Attribute '//TRIM(att_name)//' not found in NetCDF file')
    IF(att_len /= 1) &
      & CALL finish(modname,'Attribute '//TRIM(att_name)//' has invalid length')

    CALL nf(nf_get_att_int(ncid, nf_global, att_name, val))
    IF(val /= att_value) THEN
      WRITE(message_text,'(a,a,i0,a,i0)') att_name,': NetCDF value = ',val,' expected: ',att_value
      CALL finish(modname,message_text)
    ENDIF

  END SUBROUTINE check_att_int

  !-------------------------------------------------------------------------
  !
  !> check_att_double
  !! Checks if attribute "att_name" in NetCDF file has given value

  SUBROUTINE check_att_double(att_name, att_value)

    CHARACTER(LEN=*), INTENT(IN) :: att_name
    REAL(wp), INTENT(IN) :: att_value

    INTEGER :: att_len, stat
    REAL(wp) :: val

    stat = nf_inq_attlen (ncid, nf_global, att_name, att_len)
    IF(stat /= nf_noerr) &
      & CALL finish(modname,'Attribute '//TRIM(att_name)//' not found in NetCDF file')
    IF(att_len /= 1) &
      & CALL finish(modname,'Attribute '//TRIM(att_name)//' has invalid length')

    CALL nf(nf_get_att_double(ncid, nf_global, att_name, val))
    IF(val /= att_value) THEN
      WRITE(message_text,'(a,a,g20.10,a,g20.10)') &
        &  att_name,': NetCDF value = ',val,' expected: ',att_value
      CALL finish(modname,message_text)
    ENDIF

  END SUBROUTINE check_att_double

  !-------------------------------------------------------------------------
  !
  !> check_att_int_array
  !! Checks if attribute "att_name" in NetCDF file has given value

  SUBROUTINE check_att_int_array(att_name, att_value)

    CHARACTER(LEN=*), INTENT(IN) :: att_name
    INTEGER, INTENT(IN) :: att_value(:)

    INTEGER :: att_len, stat
    INTEGER :: val(SIZE(att_value(:),1))

    stat = nf_inq_attlen (ncid, nf_global, att_name, att_len)
    IF(stat /= nf_noerr) &
      & CALL finish(modname,'Attribute '//TRIM(att_name)//' not found in NetCDF file')
    IF(att_len /= SIZE(att_value(:),1)) &
      & CALL finish(modname,'Attribute '//TRIM(att_name)//' has invalid length')

    CALL nf(nf_get_att_int(ncid, nf_global, att_name, val(:)))
    IF (ANY(val(:) /= att_value(:))) THEN
      WRITE(message_text,'(a,a,i0,a,i0)') att_name,': NetCDF value = ',val,' expected: ',att_value
      CALL finish(modname,message_text)
    ENDIF

  END SUBROUTINE check_att_int_array


  !-------------------------------------------------------------------------
  !
  !> check_att_double_array
  !! Checks if attribute "att_name" in NetCDF file has given value

  SUBROUTINE check_att_double_array(att_name, att_value)

    CHARACTER(LEN=*), INTENT(IN) :: att_name
    REAL(wp), INTENT(IN) :: att_value(:)

    INTEGER  :: att_len, stat
    REAL(wp) :: val(SIZE(att_value(:),1))

    stat = nf_inq_attlen (ncid, nf_global, att_name, att_len)
    IF(stat /= nf_noerr) &
      & CALL finish(modname,'Attribute '//TRIM(att_name)//' not found in NetCDF file')
    IF(att_len /= SIZE(att_value(:),1)) &
      & CALL finish(modname,'Attribute '//TRIM(att_name)//' has invalid length')

    CALL nf(nf_get_att_double(ncid, nf_global, att_name, val(:)))
    IF(ANY(val(:) /= att_value(:))) THEN
      WRITE(message_text,'(a,a,g20.10,a,g20.10)') &
        &  att_name,': NetCDF value = ',val,' expected: ',att_value
      CALL finish(modname,message_text)
    ENDIF

  END SUBROUTINE check_att_double_array
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !> restore_dim
  !! Checks if dimension "dim_name" in NetCDF file has given value
  !! Restores the given value to the netcdf value, if they differ
  SUBROUTINE restore_dim(dim_name, dim_value)

    CHARACTER(LEN=*), INTENT(IN) :: dim_name
    INTEGER, INTENT(INOUT) :: dim_value

    INTEGER :: dimid, stat, val

    stat = nf_inq_dimid (ncid, dim_name, dimid)
    IF(stat /= nf_noerr) &
      & CALL finish(modname,'Dimension '//TRIM(dim_name)//' not found in NetCDF file')

    CALL nf(nf_inq_dimlen(ncid, dimid, val))
    IF(val /= dim_value) THEN
      WRITE(message_text,'(a,a,i0,a,i0)') dim_name,': NetCDF value = ',val,' expected: ',dim_value
      CALL warning(modname,message_text)
      WRITE(message_text,'(a,i0)') " Expected value will be restored to:", val
      CALL warning(modname,message_text)
      dim_value = val
    ENDIF

  END SUBROUTINE restore_dim
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !> check_dim
  !! Checks if dimension "dim_name" in NetCDF file has given value
  SUBROUTINE check_dim(dim_name, dim_value)

    CHARACTER(LEN=*), INTENT(IN) :: dim_name
    INTEGER, INTENT(IN) :: dim_value

    INTEGER :: dimid, stat, val

    stat = nf_inq_dimid (ncid, dim_name, dimid)
    IF(stat /= nf_noerr) &
      & CALL finish(modname,'Dimension '//TRIM(dim_name)//' not found in NetCDF file')

    CALL nf(nf_inq_dimlen(ncid, dimid, val))
    IF(val /= dim_value) THEN
      WRITE(message_text,'(a,a,i0,a,i0)') dim_name,': NetCDF value = ',val,' expected: ',dim_value
      CALL finish(modname,message_text)
    ENDIF

  END SUBROUTINE check_dim

  !
  !-------------------------------------------------------------------------
  !> set_atts_and_dims:
  !! Sets common attributes and dimensions in NetCDF file

  SUBROUTINE set_atts_and_dims(p, nprocs)

    TYPE(t_patch), INTENT(IN) :: p
    INTEGER, INTENT(IN) :: nprocs

    INTEGER :: i
    CHARACTER (LEN=80) :: child_id_name, child_idl_name

    ! Output all values relevant for grid as attributes
    ! (unless they are used as dimensions below)

    CALL nf(nf_put_att_int(ncid, nf_global, 'num_work_procs', nf_int, 1, nprocs))
    CALL nf(nf_put_att_int(ncid, nf_global, 'start_lev', nf_int, 1, start_lev))
    CALL nf(nf_put_att_int(ncid, nf_global, 'n_dom', nf_int, 1, n_dom))
    IF(p%id > 0) CALL nf(nf_put_att_int(ncid, nf_global, 'lfeedback', nf_int, 1, &
      &  MERGE(1,0,lfeedback(p%id))))
    CALL nf(nf_put_att_int(ncid, nf_global, 'l_limited_area', nf_int, 1, &
      &  MERGE(1,0,l_limited_area)))

    CALL nf(nf_put_att_int(ncid, nf_global, 'patch.parent_id', nf_int, 1, p%parent_id))
    CALL nf(nf_put_att_int(ncid, nf_global, 'patch.parent_child_index', nf_int, 1, &
      &  p%parent_child_index))
    CALL nf(nf_put_att_int(ncid, nf_global, 'patch.n_childdom', nf_int, 1, p%n_childdom))
    CALL nf(nf_put_att_int(ncid, nf_global, 'patch.n_chd_total', nf_int, 1, p%n_chd_total))
    !
    DO i = 1, p%n_childdom
      WRITE(child_id_name,'(a,i0)') 'patch.child_id',i
      CALL nf(nf_put_att_int(ncid, nf_global, child_id_name, nf_int, 1, p%child_id(i)))
    ENDDO
    DO i = 1, p%n_chd_total
      WRITE(child_idl_name,'(a,i0)') 'patch.child_id_list',i
      CALL nf(nf_put_att_int(ncid, nf_global, child_idl_name, nf_int, 1, p%child_id_list(i)))
    ENDDO
    CALL nf(nf_put_att_int(ncid, nf_global, 'patch.n_proc', nf_int, 1, p%n_proc))
    CALL nf(nf_put_att_int(ncid, nf_global, 'patch.proc0', nf_int, 1, p%proc0))

    ! Set common dimensions which do not depend on the actual PE or patch.

    CALL def_dim('unlimited', NF_UNLIMITED, dim_unlimited)

    CALL def_dim('n_dim_2', 2, dim_2)
    CALL def_dim('n_dim_3', 3, dim_3)
    CALL def_dim('n_dim_4', 4, dim_4)
    CALL def_dim('n_dim_5', 5, dim_5)
    CALL def_dim('n_dim_6', 6, dim_6)
    CALL def_dim('n_dim_7', 7, dim_7)
    CALL def_dim('n_dim_8', 8, dim_8)
    CALL def_dim('n_dim_9', 9, dim_9)

    CALL def_dim('n_rlcell', max_rlcell-min_rlcell+1, dim_nrlcell)
    CALL def_dim('n_rledge', max_rledge-min_rledge+1, dim_nrledge)
    CALL def_dim('n_rlvert', max_rlvert-min_rlvert+1, dim_nrlvert)

    IF (p%cell_type == 3) THEN
      CALL def_dim('nverts_per_cell',  3, dim_nverts_per_cell)
      CALL def_dim('nverts_per_cell_plus_1',  4, dim_nverts_per_cell_p1)
      CALL def_dim('nedges_per_vert',  6, dim_nedges_per_vert)
    ELSE
      CALL def_dim('nverts_per_cell',  6, dim_nverts_per_cell)
      CALL def_dim('nverts_per_cell_plus_1',  7, dim_nverts_per_cell_p1)
      CALL def_dim('nedges_per_vert',  3, dim_nedges_per_vert)
    ENDIF

  END SUBROUTINE set_atts_and_dims

  !-------------------------------------------------------------------------
  ! Public interface routines
  !-------------------------------------------------------------------------
  !
  !> Dumps state of one patch to NetCDF, including interpolation and
  !! grid refinement variables

  SUBROUTINE dump_patch_state_netcdf(k_jg, p, pi, pg, opt_pi_lonlat)

    INTEGER, INTENT(IN)                            :: k_jg  ! domain index
    TYPE(t_patch), INTENT(INOUT)                   :: p
    TYPE(t_int_state), INTENT(INOUT)               :: pi
    TYPE(t_gridref_state), INTENT(INOUT)           :: pg
    TYPE (t_lon_lat_intp), INTENT(INOUT), OPTIONAL :: opt_pi_lonlat

    !---local variables
    INTEGER :: old_mode, i, ip
    TYPE(t_patch), POINTER :: lp ! pointer to local parent

    !-------------------------------------------------------------------------

    ! use_one_file is set according to l_one_file_per_patch
    use_one_file = l_one_file_per_patch

    netcdf_read = .FALSE. ! Set I/O routines to write mode

    CALL set_dump_restore_filename(p%grid_filename, model_base_dir)

    WRITE(message_text,'(a,a)') 'Write NetCDF file: ', TRIM(filename)
    CALL message ('', TRIM(message_text))

    ! Set pointer to local parent
    IF(my_process_is_mpi_all_parallel() .AND. p%id>n_dom_start) THEN
      lp => p_patch_local_parent(p%id)
    ELSE
      lp => NULL()
    ENDIF

    IF(use_one_file) THEN
      ! Only the first PE defines common dimensions and attributes, all others check
      output_defines = (get_my_mpi_work_id()==0)
    ELSE
      ! Every PE defines common dimensions and attributes
      output_defines = .TRUE.
    ENDIF

    prefix = ' '

    IF(output_defines) THEN
      ! Open new output file
      CALL nf(nf_set_default_format(nf_format_64bit, old_mode))
      CALL nf(nf_create(TRIM(filename), nf_clobber, ncid))
    ELSE
      ncid = -1
    ENDIF

    ! Output attributes and dimensions

    IF(output_defines) CALL set_atts_and_dims(p, num_work_procs)

    ! Emit additional dimensions and variable definitions for patch/int state/grf state

    CALL set_max_patch_dims(p)
    CALL def_patch(p, .TRUE.)
    CALL def_comm_patterns(p)
    CALL def_int_state(k_jg, p)
    IF (n_dom_start==0 .OR. n_dom > 1) CALL def_grf_state(p)

    ! The lonlat state is always only defined and output on work PE 0
    IF (k_jg>0) THEN
      IF(get_my_mpi_work_id()==0 .AND. lonlat_intp_config(k_jg)%l_enabled) &
        & CALL def_lonlat_state(k_jg)
    ENDIF

    ! Definitions for local parent

    IF(my_process_is_mpi_all_parallel() .AND. p%id>n_dom_start) THEN
      prefix = 'lp.'
      CALL set_max_patch_dims(lp)
      CALL def_patch(lp, .TRUE.)
      CALL def_comm_patterns(lp)
      CALL def_int_state(k_jg, lp)
      IF (n_dom_start==0 .OR. n_dom > 1) CALL def_grf_state(lp)
    ENDIF

    !-------------------------------------------------------------------------

    ! End of definition mode

    IF(output_defines) THEN
      CALL nf(nf_enddef(ncid))
      CALL nf(nf_close(ncid))
    ENDIF


    IF(use_one_file) THEN
      my_record = get_my_mpi_work_id()+1
    ELSE
      my_record = 1
    ENDIF

    ! Output all variables

    DO ip = 0, num_work_procs-1

      IF(use_one_file) CALL p_barrier(comm=p_comm_work)

      IF(ip /= get_my_mpi_work_id()) CYCLE

      CALL nf(nf_open(TRIM(filename), NF_WRITE, ncid))
      prefix = ' '

      CALL patch_io(p, .TRUE.)

      CALL int_state_io(k_jg, p, pi)

      IF (n_dom_start==0 .OR. n_dom > 1) CALL grf_state_io(p, pg)

      ! The lonlat state is always only defined and output on work PE 0
      IF (k_jg>0) THEN
        IF(get_my_mpi_work_id()==0 .AND. lonlat_intp_config(k_jg)%l_enabled .AND. &
           PRESENT(opt_pi_lonlat)) CALL lonlat_state_io(opt_pi_lonlat)
      ENDIF

      IF(my_process_is_mpi_all_parallel() .AND. p%id>n_dom_start) THEN ! Output for local parent
        prefix = 'lp.'
        CALL patch_io(lp, .TRUE.)
        CALL int_state_io(k_jg, lp, p_int_state_local_parent(p%id))
        IF (n_dom_start==0 .OR. n_dom > 1) &
          & CALL grf_state_io(lp, p_grf_state_local_parent(p%id))
      ENDIF

      ! Close NetCDF file
      CALL nf(nf_close(ncid))

    ENDDO

    IF(use_one_file) CALL p_barrier(comm=p_comm_work)

  END SUBROUTINE dump_patch_state_netcdf

  !-------------------------------------------------------------------------
  !
  !> Dumps domain decomposition and some additional variables which
  !! have been calculated for domain decomposition

  SUBROUTINE dump_domain_decomposition(p, lp)

    TYPE(t_patch), INTENT(INOUT) :: p                ! basic patch
    TYPE(t_patch), INTENT(INOUT), OPTIONAL :: lp     ! local parent

    !---local variables
    INTEGER :: old_mode, ip

    !-------------------------------------------------------------------------

    ! Currently, this is implemented only for one file per patch,
    ! independent of l_one_file_per_patch
    use_one_file = .TRUE.

    netcdf_read = .FALSE. ! Set I/O routines to write mode

    CALL set_dd_filename(p%grid_filename, model_base_dir)

    WRITE(message_text,'(a,a)') 'Write domain decomposition NetCDF file: ', TRIM(filename)
    CALL message ('', TRIM(message_text))

    ! Only the first PE defines common dimensions and attributes, all others check
    output_defines = (get_my_mpi_work_id()==0)

    prefix = ' '

    IF(output_defines) THEN
      ! Open new output file
      CALL nf(nf_set_default_format(nf_format_64bit, old_mode))
      CALL nf(nf_create(TRIM(filename), nf_clobber, ncid))
    ELSE
      ncid = -1
    ENDIF

    ! Output attributes and dimensions

    IF(output_defines) CALL set_atts_and_dims(p, num_work_procs)

    ! Emit additional dimensions and variable definitions for patch

    CALL set_max_patch_dims(p)
    CALL def_patch(p, .FALSE.)

    ! Definitions for local parent

    IF(PRESENT(lp)) THEN
      prefix = 'lp.'
      CALL set_max_patch_dims(lp)
      CALL def_patch(lp, .FALSE.)
    ENDIF

    !-------------------------------------------------------------------------

    ! End of definition mode

    IF(output_defines) THEN
      CALL nf(nf_enddef(ncid))
      CALL nf(nf_close(ncid))
    ENDIF

    my_record = get_my_mpi_work_id()+1

    ! Output all variables

    DO ip = 0, num_work_procs-1

      CALL p_barrier(comm=p_comm_work)

      IF(ip /= get_my_mpi_work_id()) CYCLE

      CALL nf(nf_open(TRIM(filename), NF_WRITE, ncid))
      prefix = ' '

      CALL patch_io(p, .FALSE.)

      IF(PRESENT(lp)) THEN
        prefix = 'lp.'
        CALL patch_io(lp, .FALSE.)
      ENDIF

      ! Close NetCDF file
      CALL nf(nf_close(ncid))

    ENDDO

    CALL p_barrier(comm=p_comm_work)

  END SUBROUTINE dump_domain_decomposition

  !-------------------------------------------------------------------------------------------------
  !
  !> Opens the domain decomposition dump file and dumps all patches.
  !! This routine has to be used if domain decomposition by a single proc

  SUBROUTINE dump_all_domain_decompositions(p, lp)

    TYPE(t_patch), INTENT(INOUT)           :: p(:)      ! basic patch
    TYPE(t_patch), INTENT(INOUT), OPTIONAL :: lp(:)     ! local parent

    !---local variables
    INTEGER :: old_mode, n

    !-------------------------------------------------------------------------

    ! Currently, this is implemented only for one file per patch,
    ! independent of l_one_file_per_patch
    use_one_file = .TRUE.

    netcdf_read = .FALSE. ! Set I/O routines to write mode

    CALL set_dd_filename(p(1)%grid_filename, model_base_dir)

    WRITE(message_text,'(a,a)') 'Write domain decomposition NetCDF file: ', TRIM(filename)
    CALL message ('', TRIM(message_text))

    output_defines = .TRUE.

    ! Open new output file
    CALL nf(nf_set_default_format(nf_format_64bit, old_mode))
    CALL nf(nf_create(TRIM(filename), nf_clobber, ncid))

    prefix = ' '

    ! Output attributes and dimensions

    CALL set_atts_and_dims(p(1), UBOUND(p,1))

    ! Emit additional dimensions and variable definitions for patch

    max_patch_cells = MAXVAL(p(:)%n_patch_cells)
    max_patch_edges = MAXVAL(p(:)%n_patch_edges)
    max_patch_verts = MAXVAL(p(:)%n_patch_verts)

    CALL def_patch(p(1), .FALSE.)

    ! Definitions for local parent

    IF(PRESENT(lp)) THEN
      prefix = 'lp.'
      max_patch_cells = MAXVAL(lp(:)%n_patch_cells)
      max_patch_edges = MAXVAL(lp(:)%n_patch_edges)
      max_patch_verts = MAXVAL(lp(:)%n_patch_verts)
      CALL def_patch(lp(1), .FALSE.)
    ENDIF

    ! End of definition mode

    CALL nf(nf_enddef(ncid))

    DO n = 1, UBOUND(p,1)

      my_record = n

      prefix = ' '
      CALL patch_io(p(n), .FALSE.)

      IF(PRESENT(lp)) THEN
        prefix = 'lp.'
        CALL patch_io(lp(n), .FALSE.)
      ENDIF

    ENDDO

    ! Close NetCDF file
    CALL nf(nf_close(ncid))

  END SUBROUTINE dump_all_domain_decompositions

  !-------------------------------------------------------------------------------------------------
  !
  !> Restores a single patch from NetCDF

  SUBROUTINE restore_patch_netcdf(p, lfull)

    TYPE(t_patch), INTENT(INOUT) :: p
    LOGICAL, INTENT(IN) :: lfull

    !---local variables
    INTEGER :: dimid, varid, jb, jl

    !-------------------------------------------------------------------------

    ! The patch must have the following members already set when calling this routine:
    ! - id
    ! - level
    ! - parent_id
    ! - parent_child_index
    ! - n_childdom
    ! - n_chd_total
    ! - child_id(1:n_childdom)
    ! - child_id_list(1:n_chd_total)
    ! - max_childdom
    ! - n_patch_cells_g
    ! - n_patch_edges_g
    ! - n_patch_verts_g
    ! All other members are set here, all allocatable arrays must still be unallocated.


    ! Get dimensions needed for allocating patch

    CALL nf(nf_inq_varid(ncid, TRIM(prefix)//'n_patch_cells', varid))
    CALL nf(nf_get_var1_int(ncid, varid, my_record, p%n_patch_cells))
    CALL nf(nf_inq_varid(ncid, TRIM(prefix)//'n_patch_edges', varid))
    CALL nf(nf_get_var1_int(ncid, varid, my_record, p%n_patch_edges))
    CALL nf(nf_inq_varid(ncid, TRIM(prefix)//'n_patch_verts', varid))
    CALL nf(nf_get_var1_int(ncid, varid, my_record, p%n_patch_verts))

    CALL nf(nf_inq_dimid (ncid, TRIM(prefix)//'n_patch_cells_g', dimid))
    CALL nf(nf_inq_dimlen(ncid, dimid, p%n_patch_cells_g))

    CALL nf(nf_inq_dimid (ncid, TRIM(prefix)//'n_patch_edges_g', dimid))
    CALL nf(nf_inq_dimlen(ncid, dimid, p%n_patch_edges_g))

    CALL nf(nf_inq_dimid (ncid, TRIM(prefix)//'n_patch_verts_g', dimid))
    CALL nf(nf_inq_dimlen(ncid, dimid, p%n_patch_verts_g))


    ! calculate and save values for the blocking
    !
    ! ... for the cells
    p%nblks_c       = blk_no(p%n_patch_cells)
    p%nblks_int_c   = p%nblks_c
    p%npromz_c      = p%n_patch_cells - (p%nblks_c - 1)*nproma
    p%npromz_int_c  = p%npromz_c

    ! ... for the edges
    p%nblks_e       = blk_no(p%n_patch_edges)
    p%nblks_int_e   = p%nblks_e
    p%npromz_e      = p%n_patch_edges - (p%nblks_e - 1)*nproma
    p%npromz_int_e  = p%npromz_e

    ! ... for the vertices
    p%nblks_v       = blk_no(p%n_patch_verts)
    p%nblks_int_v   = p%nblks_v
    p%npromz_v      = p%n_patch_verts - (p%nblks_v - 1)*nproma
    p%npromz_int_v  = p%npromz_v

    ! Allocate patch

    CALL allocate_patch(p)

    CALL patch_io(p, lfull)

    ! Restore local index arrays since these are not saved due to size
    CALL restore_loc_index(p%n_patch_cells, p%n_patch_cells_g, p%cells%glb_index, &
                           p%cells%loc_index)
    CALL restore_loc_index(p%n_patch_edges, p%n_patch_edges_g, p%edges%glb_index, &
                           p%edges%loc_index)
    CALL restore_loc_index(p%n_patch_verts, p%n_patch_verts_g, p%verts%glb_index, &
                           p%verts%loc_index)

  END SUBROUTINE restore_patch_netcdf

  !-------------------------------------------------------------------------
  !> Restores the local index since this is not saved

  SUBROUTINE restore_loc_index(n, n_g, glb_index, loc_index)

    INTEGER, INTENT(IN)    :: n    ! Number of local points
    INTEGER, INTENT(IN)    :: n_g  ! Number of global points
    INTEGER, INTENT(IN)    :: glb_index(:) ! Global index (for local points)
    INTEGER, INTENT(INOUT) :: loc_index(:) ! Local index to be set

    INTEGER :: i, last

    loc_index(:) = -1

    DO i = 1, n
      IF(glb_index(i)<1 .OR. glb_index(i)>UBOUND(loc_index,1)) &
        CALL finish(modname,'Got illegal global index')
      loc_index(glb_index(i)) = i
    ENDDO

    ! Set the negative values (indicating non local points) in the same way
    ! as during subdivision: To the negative number of the next value to come
    ! This shouldn't be necessary but it also doesn't hurt

    last = 0
    DO i = 1, n_g
      IF(loc_index(i)>0) THEN
        ! Increase last (unless it is an interspersed boundary point)
        IF(loc_index(i) == last+1) last = last+1
      ELSE
        loc_index(i) = -(last+1)
      ENDIF
    ENDDO

  END SUBROUTINE restore_loc_index

  !-------------------------------------------------------------------------
  !
  !> Restores all patches from NetCDF

  SUBROUTINE restore_patches_netcdf(p_patch, lfull)

    TYPE(t_patch), INTENT(inout) :: p_patch(n_dom_start:)
    LOGICAL, INTENT(IN) :: lfull ! TRUE: restore all, FALSE: restore domain decomposition

    INTEGER :: jg, jg1, jgp, i, n_chd, n_chdc, dimid
    CHARACTER (LEN=80) :: child_id_name, child_idl_name

    CALL message ('restore_patches_netcdf','start to restore patches')

    ! Currently, restoring domain decompositions (lfull==.FALSE.) is implemented only for 
    ! one file per patch, full restores are according to l_one_file_per_patch

    IF(lfull) THEN
      use_one_file = l_one_file_per_patch
    ELSE
      use_one_file = .TRUE.
    ENDIF

    CALL set_patches_grid_filename(p_patch)

    netcdf_read = .TRUE. ! Set I/O routines to read mode

    ! Set my record in input file
    IF(use_one_file) THEN
      my_record = get_my_mpi_work_id()+1
    ELSE
      my_record = 1
    ENDIF

    ! Set some basic flow control variables on the patch
    ! This is copied from mo_model_domimp_patches/import_patches !!!

    max_childdom = 0

    ! Set cell type to the value specified in the namelists for all patches
    p_patch(n_dom_start:n_dom)%cell_type = global_cell_type

    IF(n_dom_start==0) THEN
      ! The physic parent (parent of the root patch) should also be read
      p_patch(0)%id = 0
      p_patch(0)%level = start_lev-1
      p_patch(0)%parent_id = -1
      p_patch(0)%parent_child_index = 0
      p_patch(0)%n_childdom = 1
      p_patch(0)%n_chd_total = n_dom
      p_patch(0)%child_id(1) = 1
      DO jg = 1, n_dom
        p_patch(0)%child_id_list(jg) = jg
      ENDDO

      p_patch(1)%parent_child_index = 1
    ELSE
      p_patch(1)%parent_child_index = 0
    ENDIF

    DO jg = 1, n_dom

      p_patch(jg)%id = jg

      IF (jg == 1) THEN
        p_patch(jg)%level = start_lev
        p_patch(jg)%parent_id = 0
      ELSE
        ! Note: the first element of parent_id refers to jg=2
        p_patch(jg)%level = p_patch(dynamics_parent_grid_id(jg))%level + 1
        p_patch(jg)%parent_id = dynamics_parent_grid_id(jg)
      ENDIF

      n_chd = 0

      DO jg1 = jg+1, n_dom
        IF (jg == dynamics_parent_grid_id(jg1)) THEN
          n_chd = n_chd + 1
          p_patch(jg)%child_id(n_chd) = jg1
          p_patch(jg1)%parent_child_index = n_chd
        ENDIF
      ENDDO

      p_patch(jg)%n_childdom = n_chd
      max_childdom = MAX(1,max_childdom,n_chd)

      !
      ! store information about vertical levels
      !
      p_patch(jg)%nlev   = num_lev(jg)
      p_patch(jg)%nlevp1 = num_levp1(jg)

      IF (jg > 1) THEN
        IF (nshift(jg) > 0 ) THEN
          ! nshift has been modified via Namelist => use it
          p_patch(jg)%nshift = nshift(jg)
        ELSE
          ! set default value, assuming
          !- superimposed vertical levels
          !- 1 nested domain per grid level
          p_patch(jg)%nshift = num_lev(p_patch(jg)%parent_id) - num_lev(jg)
        ENDIF

        jgp = p_patch(jg)%parent_id
        p_patch(jg)%nshift_total = p_patch(jgp)%nshift_total + p_patch(jg)%nshift
      ELSE
        ! Note: the first nshift-value refers to the global domain
        p_patch(jg)%nshift = 0
        p_patch(jg)%nshift_total = 0
      ENDIF
    ENDDO

    ! Set information about total number of child domains (called recursively)
    ! and corresponding index lists

    ! Initialization
    DO jg = 1, n_dom
      p_patch(jg)%n_chd_total      = 0
      p_patch(jg)%child_id_list(:) = 0
    ENDDO

    DO jg = n_dom, 2, -1
      jg1 = p_patch(jg)%parent_id
      n_chd = p_patch(jg1)%n_chd_total
      n_chdc = p_patch(jg)%n_chd_total
      p_patch(jg1)%child_id_list(n_chd+1) = jg
      IF (n_chdc > 0) THEN
        p_patch(jg1)%child_id_list(n_chd+2:n_chd+1+n_chdc) = p_patch(jg)%child_id_list(1:n_chdc)
      ENDIF
      p_patch(jg1)%n_chd_total = n_chd+1+n_chdc
    ENDDO

    DO jg = 1, n_dom
      ! Set nshift_child in the same way as in domimp_patches
      IF (p_patch(jg)%n_childdom >= 1) THEN
        p_patch(jg)%nshift_child = p_patch(p_patch(jg)%child_id(1))%nshift
        DO jg1 = 1, p_patch(jg)%n_childdom
          IF (p_patch(p_patch(jg)%child_id(jg1))%nshift /= p_patch(jg)%nshift_child) &
          CALL finish ('mo_dump_restore:restore_patches_netcdf', &
                       'multiple nests at the same level must have the same nshift')
        ENDDO
      ELSE
        p_patch(jg)%nshift_child = 0
      ENDIF
    ENDDO

    IF (n_dom_start == 0) THEN ! reduced grid for radiation
      ! In case of n_dom_start == 0 nlev, nlevp1, nshift need to be copied from
      ! jg=1 to jg=0
      p_patch(0)%nlev   = p_patch(1)%nlev
      p_patch(0)%nlevp1 = p_patch(1)%nlevp1
      p_patch(0)%nshift = p_patch(1)%nshift
    ENDIF

    p_patch(n_dom_start:n_dom)%max_childdom =  max_childdom

    IF(my_process_is_mpi_all_parallel()) THEN

      ALLOCATE(p_patch_local_parent(n_dom_start+1:n_dom))

      DO jg = n_dom_start+1, n_dom
        jgp = p_patch(jg)%parent_id
        p_patch_local_parent(jg)%id           = p_patch(jgp)%id
        p_patch_local_parent(jg)%level        = p_patch(jgp)%level
        p_patch_local_parent(jg)%parent_id    = p_patch(jgp)%parent_id
        p_patch_local_parent(jg)%parent_child_index = p_patch(jgp)%parent_child_index
        p_patch_local_parent(jg)%n_childdom   = p_patch(jgp)%n_childdom
        p_patch_local_parent(jg)%n_chd_total  = p_patch(jgp)%n_chd_total
        p_patch_local_parent(jg)%child_id(:)  = p_patch(jgp)%child_id(:)
        p_patch_local_parent(jg)%child_id_list(:) = p_patch(jgp)%child_id_list(:)
        p_patch_local_parent(jg)%max_childdom = p_patch(jgp)%max_childdom
        p_patch_local_parent(jg)%comm   = 0 ! Not needed
        p_patch_local_parent(jg)%rank   = 0 ! Not needed
        p_patch_local_parent(jg)%n_proc = 0 ! Not needed
        p_patch_local_parent(jg)%proc0  = 0 ! Not needed
        p_patch_local_parent(jg)%nlev   = p_patch(jgp)%nlev
        p_patch_local_parent(jg)%nlevp1 = p_patch(jgp)%nlevp1
        p_patch_local_parent(jg)%nshift = 0
        p_patch_local_parent(jg)%nshift_child = 0
      ENDDO
    ENDIF

    DO jg = n_dom_start, n_dom

      IF(p_pe_work==0) WRITE(nerr,'(a,i0)') 'Restoring patch ',jg

      IF(lfull) THEN
        CALL set_dump_restore_filename(p_patch(jg)%grid_filename, model_base_dir)
      ELSE
        CALL set_dd_filename(p_patch(jg)%grid_filename, model_base_dir)
      ENDIF
!       write(0,*) "patch grid_filename:", TRIM( p_patch(jg)%grid_filename)
!       write(0,*) "dump_restore_filename:", TRIM(filename)

      prefix = ' '

      CALL nf(nf_open(TRIM(filename), NF_NOWRITE, ncid))
!       write(0,*) TRIM(filename), " is open."

      ! First check if members in patch and variables from input
      ! conform with what is stored in NetCDF file

      CALL check_att('num_work_procs', num_work_procs)
      CALL check_att('start_lev',      start_lev)
      CALL check_att('n_dom',          n_dom)
      IF(jg>0) CALL check_att('lfeedback', MERGE(1,0,lfeedback(p_patch(jg)%id)))
      CALL check_att('l_limited_area', MERGE(1,0,l_limited_area))

      CALL check_att('patch.parent_id', p_patch(jg)%parent_id)
      CALL check_att('patch.parent_child_index', p_patch(jg)%parent_child_index)
      CALL check_att('patch.n_childdom', p_patch(jg)%n_childdom)
      CALL check_att('patch.n_chd_total', p_patch(jg)%n_chd_total)

      DO i = 1, p_patch(jg)%n_childdom
        WRITE(child_id_name,'(a,i0)') 'patch.child_id',i
        CALL check_att(child_id_name, p_patch(jg)%child_id(i))
      ENDDO
      DO i = 1, p_patch(jg)%n_chd_total
        WRITE(child_idl_name,'(a,i0)') 'patch.child_id_list',i
        CALL check_att(child_idl_name, p_patch(jg)%child_id_list(i))
      ENDDO

      CALL check_dim('max_childdom', p_patch(jg)%max_childdom)

      ! nverts_per_cell must conform to p%cell_type
      CALL restore_dim('nverts_per_cell', p_patch(jg)%cell_type)

      ! Just for safety
      CALL check_dim('n_rlcell', max_rlcell-min_rlcell+1)
      CALL check_dim('n_rledge', max_rledge-min_rledge+1)
      CALL check_dim('n_rlvert', max_rlvert-min_rlvert+1)

      ! Set scalar members of patch

      CALL nf(nf_get_att_int(ncid, nf_global, 'patch.n_proc', p_patch(jg)%n_proc))
      CALL nf(nf_get_att_int(ncid, nf_global, 'patch.proc0', p_patch(jg)%proc0))

      CALL restore_patch_netcdf(p_patch(jg), lfull)

      IF(my_process_is_mpi_all_parallel() .AND. jg>n_dom_start) THEN
        prefix = 'lp.'
        CALL restore_patch_netcdf(p_patch_local_parent(jg), lfull)
      ENDIF

      CALL nf(nf_close(ncid))

    ENDDO

  END SUBROUTINE restore_patches_netcdf

  !-------------------------------------------------------------------------
  !
  !> Restores interpolation state from NetCDF for all patches

  SUBROUTINE restore_interpol_state_netcdf(p_patch, p_int_state, opt_pi_lonlat)

    TYPE(t_patch),        INTENT(INOUT) :: p_patch(n_dom_start:)
    TYPE(t_int_state),    INTENT(INOUT) :: p_int_state(n_dom_start:)
    TYPE (t_lon_lat_intp), INTENT(INOUT), OPTIONAL :: opt_pi_lonlat(n_dom_start:)

    !---local variables
    INTEGER :: jg
    TYPE (t_lon_lat_grid), POINTER :: grid ! lon-lat grid
    LOGICAL                        :: l_restore_lonlat

    !-------------------------------------------------------------------------

    ! use_one_file is set according to l_one_file_per_patch
    use_one_file = l_one_file_per_patch

    netcdf_read = .TRUE. ! Set I/O routines to read mode

    ! Set my record in input file
    IF(use_one_file) THEN
      my_record = get_my_mpi_work_id()+1
    ELSE
      my_record = 1
    ENDIF

    DO jg = n_dom_start, n_dom

      ! The lonlat state is always only defined and output on work PE 0
      IF (jg>0) THEN
        l_restore_lonlat = get_my_mpi_work_id()==0 .AND. &
          &  lonlat_intp_config(jg)%l_enabled .AND. PRESENT(opt_pi_lonlat)
      ELSE
        l_restore_lonlat = .FALSE.
      ENDIF

      IF(p_pe_work==0) WRITE(nerr,'(a,i0)') 'Restoring interpolation state ',jg

      CALL set_dump_restore_filename(p_patch(jg)%grid_filename, model_base_dir)

      CALL nf(nf_open(TRIM(filename), NF_NOWRITE, ncid))

      ! Check attributes and dimensions in NetCDF file if they conform

      CALL check_att('llsq_high_consv',   MERGE(1,0,llsq_high_consv))
      CALL check_att('rbf_vec_kern_c',   rbf_vec_kern_c)
      CALL check_att('rbf_vec_kern_e',   rbf_vec_kern_e)
      CALL check_att('rbf_vec_kern_v',   rbf_vec_kern_v)
      CALL check_att('i_cori_method',    i_cori_method)
      CALL check_att('nudge_zone_width', nudge_zone_width)

      IF(jg>0) THEN
        CALL check_att('rbf_vec_scale_c', rbf_vec_scale_c(jg))
        CALL check_att('rbf_vec_scale_e', rbf_vec_scale_e(jg))
        CALL check_att('rbf_vec_scale_v', rbf_vec_scale_v(jg))
      ENDIF
      CALL check_att('nudge_max_coeff',   nudge_max_coeff)
      CALL check_att('nudge_efold_width', nudge_efold_width)

      CALL check_dim('lsq_dim_c_lin',   lsq_lin_set%dim_c)
      CALL check_dim('lsq_dim_unk_lin', lsq_lin_set%dim_unk)

      CALL check_dim('lsq_dim_c_high',   lsq_high_set%dim_c)
      CALL check_dim('lsq_dim_unk_high', lsq_high_set%dim_unk)

      CALL check_dim('rbf_c2grad_dim', rbf_c2grad_dim)
      CALL check_dim('rbf_vec_dim_c',  rbf_vec_dim_c)
      CALL check_dim('rbf_vec_dim_e',  rbf_vec_dim_e)
      CALL check_dim('rbf_vec_dim_v',  rbf_vec_dim_v)

      prefix = ' '

      ! Allocate interpolation state
      CALL allocate_int_state(p_patch(jg), p_int_state(jg))
      ! Restore interpolation state
      CALL int_state_io(jg, p_patch(jg), p_int_state(jg))

      IF (l_restore_lonlat) THEN

        grid => lonlat_intp_config(jg)%lonlat_grid
        CALL check_att('int.lonlat.delta'     , grid%delta(1:2))
        CALL check_att('int.lonlat.sw_corner' , grid%start_corner(1:2))
        CALL check_att('int.lonlat.poleN'     , grid%poleN(1:2))
        CALL check_att('int.lonlat.dimen'     , grid%dimen(1:2))

        ! Allocate interpolation state
        CALL allocate_int_state_lonlat(jg, opt_pi_lonlat(jg))
        ! Restore interpolation state
        CALL lonlat_state_io(opt_pi_lonlat(jg))

      END IF

      IF(my_process_is_mpi_all_parallel() .AND. jg>n_dom_start) THEN
        CALL allocate_int_state(p_patch_local_parent(jg), p_int_state_local_parent(jg))
        prefix = 'lp.'
        CALL int_state_io(jg, p_patch_local_parent(jg), p_int_state_local_parent(jg))
      ENDIF

      CALL nf(nf_close(ncid))

    ENDDO

  END SUBROUTINE restore_interpol_state_netcdf

  !-------------------------------------------------------------------------
  !
  !> Restores gridref state from NetCDF for all patches

  SUBROUTINE restore_gridref_state_netcdf(p_patch, p_grf_state)

    TYPE(t_patch), INTENT(INOUT) :: p_patch(n_dom_start:)
    TYPE(t_gridref_state), INTENT(INOUT) :: p_grf_state(n_dom_start:)

    !---local variables
    INTEGER :: jg

    !-------------------------------------------------------------------------

    ! use_one_file is set according to l_one_file_per_patch
    use_one_file = l_one_file_per_patch

    netcdf_read = .TRUE. ! Set I/O routines to read mode

    ! Set my record in input file
    IF(use_one_file) THEN
      my_record = get_my_mpi_work_id()+1
    ELSE
      my_record = 1
    ENDIF

    DO jg = n_dom_start, n_dom

      IF(p_pe_work==0) WRITE(nerr,'(a,i0)') 'Restoring gridref state ',jg

      CALL set_dump_restore_filename(p_patch(jg)%grid_filename, model_base_dir)

      CALL nf(nf_open(TRIM(filename), NF_NOWRITE, ncid))

      ! Check attributes and dimensions in NetCDF file if they conform

      CALL check_att('rbf_vec_kern_grf_e', rbf_vec_kern_grf_e)
      CALL check_att('grf_intmethod_c',    grf_intmethod_c)
      CALL check_att('grf_intmethod_ct',   grf_intmethod_ct)
      CALL check_att('grf_intmethod_e',    grf_intmethod_e)
      CALL check_att('grf_velfbk',         grf_velfbk)
      CALL check_att('grf_scalfbk',        grf_scalfbk)
      CALL check_att('grf_tracfbk',        grf_tracfbk)

      IF(jg>0) CALL check_att('rbf_scale_grf_e', rbf_scale_grf_e(jg))
      CALL check_att('grf_idw_exp_e12', grf_idw_exp_e12)
      CALL check_att('grf_idw_exp_e34', grf_idw_exp_e34)
      CALL check_att('denom_diffu_v',   denom_diffu_v)
      CALL check_att('denom_diffu_t',   denom_diffu_t)

      ! Allocate interpolation state
      CALL allocate_grf_state(p_patch(jg), p_grf_state(jg))

      ! Restore grf state
      prefix = ' '
      CALL grf_state_io(p_patch(jg), p_grf_state(jg))

      IF(my_process_is_mpi_all_parallel() .AND. jg>n_dom_start) THEN
        CALL allocate_grf_state(p_patch_local_parent(jg), p_grf_state_local_parent(jg))
        prefix = 'lp.'
        CALL grf_state_io(p_patch_local_parent(jg), p_grf_state_local_parent(jg))
      ENDIF

      CALL nf(nf_close(ncid))

    ENDDO

  END SUBROUTINE restore_gridref_state_netcdf

  !-------------------------------------------------------------------------

END MODULE mo_dump_restore
