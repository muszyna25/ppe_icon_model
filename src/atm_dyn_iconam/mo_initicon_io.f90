!>
!! This module contains the I/O routines for initicon
!!
!! @author Guenther Zaengl, DWD
!!
!!
!! @par Revision History
!! First version by Guenther Zaengl, DWD (2011-07-13)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_initicon_io

  USE mo_kind,                ONLY: wp, i8
  USE mo_io_units,            ONLY: filename_max
  USE mo_parallel_config,     ONLY: nproma, p_test_run
  USE mo_run_config,          ONLY: msg_level, iqv, iqc, iqi, iqr, iqs
  USE mo_dynamics_config,     ONLY: nnow, nnow_rcf
  USE mo_model_domain,        ONLY: t_patch
  USE mo_nonhydro_types,      ONLY: t_nh_state
  USE mo_nwp_phy_types,       ONLY: t_nwp_phy_diag
  USE mo_nwp_lnd_types,       ONLY: t_lnd_state, t_lnd_prog, t_lnd_diag, t_wtr_prog
  USE mo_initicon_types,      ONLY: t_initicon_state, t_pi_atm, alb_snow_var, geop_ml_var, &
                                    ana_varnames_dict, inventory_list_fg, inventory_list_ana
  USE mo_initicon_config,     ONLY: init_mode, nlev_in,  l_sst_in, generate_filename,   &
    &                               ifs2icon_filename, dwdfg_filename, dwdana_filename, &
    &                               nml_filetype => filetype, lread_ana, lread_vn,      &
    &                               lp2cintp_incr, lp2cintp_sfcana, ltile_coldstart
  USE mo_nh_init_nest_utils,  ONLY: interpolate_increments, interpolate_sfcana
  USE mo_impl_constants,      ONLY: SUCCESS, MAX_CHAR_LENGTH, max_dom,                  &
    &                               MODE_DWDANA_INC, MODE_IAU, MODE_IAU_OLD,            &
    &                               MODE_IFSANA, MODE_COMBINED, MODE_COSMODE
  USE mo_exception,           ONLY: message, finish, message_text
  USE mo_grid_config,         ONLY: n_dom, nroot, l_limited_area
  USE mo_mpi,                 ONLY: my_process_is_stdio, p_io, p_bcast, p_comm_work,    &
    &                               p_comm_work_test, my_process_is_mpi_workroot
  USE mo_io_config,           ONLY: default_read_method
  USE mo_read_interface,      ONLY: t_stream_id, nf, openInputFile, closeFile, &
    &                               read_2d_1time, read_2d_1lev_1time, &
    &                               read_3d_1time, onCells, onEdges
  USE mo_util_cdi,            ONLY: read_cdi_2d, read_cdi_3d, t_inputParameters,        &
    &                               makeInputParameters, deleteInputParameters,         &
    &                               get_cdi_NlevRef, t_tileinfo_elt, trivial_tileinfo
  USE mo_util_string,         ONLY: one_of
  USE mo_util_file,           ONLY: util_filesize
  USE mo_ifs_coord,           ONLY: alloc_vct, init_vct, vct, vct_a, vct_b
  USE mo_lnd_nwp_config,      ONLY: nlev_soil, ntiles_total, nlev_snow, &
    &                               ntiles_water, lmulti_snow, tiles, get_tile_suffix
  USE mo_master_nml,          ONLY: model_base_dir
  USE mo_dictionary,          ONLY: dict_get, DICT_MAX_STRLEN
  USE mo_var_metadata_types,  ONLY: VARNAME_LEN
  USE mo_cdi_constants,       ONLY: FILETYPE_NC2, FILETYPE_NC4, FILETYPE_GRB2, &
    &                               streamInqVlist, streamOpenRead, cdiGetStringError
  USE mo_nwp_sfc_interp,      ONLY: smi_to_wsoil
  USE mo_util_cdi_table,      ONLY: print_cdi_summary, &
    &                               new_inventory_list, delete_inventory_list, complete_inventory_list
  USE mo_io_util,             ONLY: get_filetype
  USE mo_initicon_utils,      ONLY: initicon_inverse_post_op, allocate_extana_atm, allocate_extana_sfc
  USE mo_physical_constants,  ONLY: cpd, rd, cvd_o_rd, p0ref, vtmpc1

  IMPLICIT NONE

  INCLUDE 'netcdf.inc'

  PRIVATE

  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_io_initicon'


  ! Remark on usage of "cdiDefMissval"

  ! Inside the GRIB_API (v.1.9.18) the missing value is converted into
  ! LONG INT for a test, but the default CDI missing value is outside
  ! of the valid range for LONG INT (U. Schulzweida, bug report
  ! SUP-277). This causes a crash with INVALID OPERATION.

  ! As a workaround we can choose a different missing value in the
  ! calling subroutine (here). For the SX-9 this must lie within 53
  ! bits, because "the SX compiler generates codes using HW
  ! instructions for floating-point data instead of instructions for
  ! integers. Therefore, the operation result is not guaranteed if the
  ! value cannot be represented as integer within 53 bits."
  DOUBLE PRECISION, PARAMETER :: cdimissval = -9.E+15


  PUBLIC :: open_init_files
  PUBLIC :: close_init_files
  PUBLIC :: read_data_2d, read_data_3d
  PUBLIC :: read_extana_atm
  PUBLIC :: read_extana_sfc
  PUBLIC :: read_dwdfg_atm
  PUBLIC :: read_dwdfg_sfc
  PUBLIC :: read_dwdana_atm
  PUBLIC :: read_dwdana_sfc
  PUBLIC :: read_dwdfg_atm_ii



  CONTAINS


  !> open files containing first guess and analysis.
  !
  !  @author F. Prill, DWD
  !
  SUBROUTINE open_init_files(p_patch, fileID_fg, fileID_ana, filetype_fg, filetype_ana, &
    &                        dwdfg_file, dwdana_file)
    TYPE(t_patch),               INTENT(IN)    :: p_patch(:)
    INTEGER,                     INTENT(INOUT) :: fileID_ana(:), fileID_fg(:)     ! dim (1:n_dom)
    INTEGER,                     INTENT(INOUT) :: filetype_ana(:), filetype_fg(:) ! dim (1:n_dom)
    CHARACTER(LEN=filename_max), INTENT(INOUT) :: dwdfg_file(max_dom)             ! first guess
    CHARACTER(LEN=filename_max), INTENT(INOUT) :: dwdana_file(max_dom)            ! analysis
    ! local variables
    CHARACTER(*), PARAMETER :: routine = modname//':open_init_files'
    INTEGER :: jg, vlistID, jlev, mpi_comm
    INTEGER(KIND=i8) :: flen_fg, flen_ana                     ! filesize in bytes
    LOGICAL :: l_exist, lread_process
    CHARACTER(LEN=MAX_CHAR_LENGTH) :: cdiErrorText

    CALL cdiDefMissval(cdimissval)
    fileID_fg(:)  = -1
    fileID_ana(:) = -1
    dwdfg_file (:)=' '
    dwdana_file(:)=' '

    ! flag. if true, then this PE reads data from file and broadcasts
    lread_process = my_process_is_mpi_workroot()

    IF(lread_process) THEN
      ! first guess ("_fg")
      !
      DO jg=1,n_dom

        IF (.NOT. p_patch(jg)%ldom_active) CYCLE
        jlev = p_patch(jg)%level
        ! generate file name
        dwdfg_file(jg) = generate_filename(dwdfg_filename, model_base_dir, nroot, jlev, jg)
        INQUIRE (FILE=dwdfg_file(jg), EXIST=l_exist)
        IF (.NOT.l_exist) THEN
          CALL finish(TRIM(routine),'DWD FG file not found: '//TRIM(dwdfg_file(jg)))
        ENDIF
        IF (nml_filetype == -1) THEN
          filetype_fg(jg) = get_filetype(TRIM(dwdfg_file(jg))) ! determine filetype
        ELSE
          filetype_fg(jg) = nml_filetype
        END IF

        ! open file
        !
        WRITE (0,"(a)") " "
        WRITE (0,"(a,a)") "file inventory: ", TRIM(dwdfg_file(jg))
        fileID_fg(jg)  = streamOpenRead(TRIM(dwdfg_file(jg)))
        ! check if the file could be opened
        IF (fileID_fg(jg) < 0) THEN
          CALL cdiGetStringError(fileID_fg(jg), cdiErrorText)
          WRITE(message_text,'(4a)') 'File ', TRIM(dwdfg_file(jg)), &
               ' cannot be opened: ', TRIM(cdiErrorText)
          CALL finish(routine, TRIM(message_text))
        ENDIF

        ! check whether the file is empty (does not work unfortunately; internal CDI error)
        flen_fg = util_filesize(TRIM(dwdfg_file(jg)))
        IF (flen_fg <= 0 ) THEN
          WRITE(message_text,'(a)') 'File '//TRIM(dwdfg_file(jg))//' is empty'
          CALL message(TRIM(routine), TRIM(message_text))
          CALL finish(routine, "Arrggh!: Empty input file")
        ENDIF

        vlistID = streamInqVlist(fileID_fg(jg))

        ! create new inventory list (linked list) for first guess input
        CALL new_inventory_list(inventory_list_fg(jg))

        ! print inventory and store in linked list
        CALL print_cdi_summary(vlistID, opt_dstlist=inventory_list_fg(jg))
        CALL complete_inventory_list(filetype_fg(jg), ana_varnames_dict, inventory_list_fg(jg) )

      END DO

      ! analysis ("_ana")
      !
      IF (lread_ana) THEN
        DO jg=1,n_dom

          IF (.NOT. p_patch(jg)%ldom_active) CYCLE
          IF (lp2cintp_incr(jg) .AND. lp2cintp_sfcana(jg)) CYCLE
          jlev = p_patch(jg)%level
          ! generate file name
          dwdana_file(jg) = generate_filename(dwdana_filename, model_base_dir, nroot, jlev, jg)
          INQUIRE (FILE=dwdana_file(jg), EXIST=l_exist)
          IF (.NOT.l_exist) THEN
            CALL finish(TRIM(routine),'DWD ANA file not found: '//TRIM(dwdana_file(jg)))
          ENDIF
          IF (nml_filetype == -1) THEN
            filetype_ana(jg) = get_filetype(TRIM(dwdana_file(jg))) ! determine filetype
          ELSE
            filetype_ana(jg) = nml_filetype
          END IF

          ! open file
          !
          WRITE (0,"(a)") " "
          WRITE (0,"(a,a)") "file inventory: ", TRIM(dwdana_file(jg))
          fileID_ana(jg)  = streamOpenRead(TRIM(dwdana_file(jg)))
          ! check if the file could be opened
          IF (fileID_ana(jg) < 0) THEN
            CALL cdiGetStringError(fileID_ana(jg), cdiErrorText)
            WRITE(message_text,'(4a)') 'File ', TRIM(dwdana_file(jg)), &
                 ' cannot be opened: ', TRIM(cdiErrorText)
            CALL finish(routine, TRIM(message_text))
          ENDIF

          ! check whether the file is empty (does not work unfortunately; internal CDI error)
          flen_ana = util_filesize(TRIM(dwdana_file(jg)))
          IF (flen_ana <= 0 ) THEN
            WRITE(message_text,'(a)') 'File '//TRIM(dwdana_file(jg))//' is empty'
            CALL message(TRIM(routine), TRIM(message_text))
            CALL finish(routine, "Arrggh!: Empty analysis file")
          ENDIF

          vlistID = streamInqVlist(fileID_ana(jg))

          ! create new inventory list (linked list) for analysis input
          CALL new_inventory_list(inventory_list_ana(jg))

          ! print inventory and store in linked list
          CALL print_cdi_summary(vlistID, opt_dstlist=inventory_list_ana(jg))
          CALL complete_inventory_list(filetype_ana(jg), ana_varnames_dict, inventory_list_ana(jg) )

        END DO
      ENDIF  ! lread_ana
    END IF  ! lread_process



    IF(p_test_run) THEN
      mpi_comm = p_comm_work_test
    ELSE
      mpi_comm = p_comm_work
    ENDIF

    CALL p_bcast(filetype_fg,  p_io, mpi_comm)
    CALL p_bcast(filetype_ana, p_io, mpi_comm)
    CALL p_bcast(fileID_fg,    p_io, mpi_comm)
    CALL p_bcast(fileID_ana,   p_io, mpi_comm)

  END SUBROUTINE open_init_files


  !> close files containing first guess and analysis.
  !
  !  @author F. Prill, DWD
  !
  SUBROUTINE close_init_files(fileID_fg, fileID_ana)
    INTEGER, INTENT(INOUT) :: fileID_ana(:), fileID_fg(:) ! dim (1:n_dom)
    ! local variables
    INTEGER :: jg
    LOGICAL :: lread_process

    ! flag. if true, then this PE reads data from file and broadcasts
    lread_process = my_process_is_mpi_workroot()

    IF(lread_process) THEN
      ! first guess ("_fg")
      DO jg=1,n_dom
        IF (fileID_fg(jg) == -1) CYCLE
        CALL streamClose(fileID_fg(jg))
        CALL delete_inventory_list(inventory_list_fg(jg))
      END DO
      ! analysis ("_ana")
      IF (lread_ana) THEN
        DO jg=1,n_dom
          IF (fileID_ana(jg) == -1) CYCLE
            CALL streamClose(fileID_ana(jg))
            CALL delete_inventory_list(inventory_list_ana(jg))
        END DO
      ENDIF
    END IF
    fileID_fg(:)  = -1
    fileID_ana(:) = -1
  END SUBROUTINE close_init_files


  !-------------------------------------------------------------------------
  !> Wrapper routine for NetCDF and GRIB2 read routines, 2D case.
  !  If necessary an inverse post_op is performed for the input field
  !
  SUBROUTINE read_data_2d (parameters, filetype, varname, var_out, tileinfo, opt_checkgroup)
    TYPE(t_inputParameters), INTENT(INOUT) :: parameters
    INTEGER,                 INTENT(IN)    :: filetype       !< FILETYPE_NC2 or FILETYPE_GRB2
    CHARACTER(len=*),        INTENT(IN)    :: varname        !< var name of field to be read
    REAL(wp), POINTER,       INTENT(INOUT) :: var_out(:,:)   !< output field
    TYPE(t_tileinfo_elt),    INTENT(IN)    :: tileinfo
    CHARACTER(LEN=VARNAME_LEN), INTENT(IN), OPTIONAL :: opt_checkgroup(:) !< read only, if varname is
                                                                          !< contained in opt_checkgroup
    ! local variables
    CHARACTER(len=*), PARAMETER     :: routine ='read_data_2d'
    CHARACTER(LEN=DICT_MAX_STRLEN)  :: mapped_name
    CHARACTER(LEN=MAX_CHAR_LENGTH)  :: filetyp_read   !< filetype for log message
    LOGICAL                         :: lread          !< .FALSE.: skip reading
    CHARACTER(LEN=MAX_CHAR_LENGTH)  :: full_varname   !< varname, including tile suffix, if present

    IF (PRESENT(opt_checkgroup)) THEN
      lread = ( one_of(varname,  opt_checkgroup(:)) /= -1)
    ELSE
      lread = .TRUE.
    ENDIF
!!$      write(0,*) "lread: varname",lread, varname


    IF (lread) THEN

      SELECT CASE(filetype)
      CASE (FILETYPE_NC2, FILETYPE_NC4)
        ! Trivial name mapping for NetCDF file
        mapped_name = TRIM(varname)

        WRITE(filetyp_read,'(a)') 'NC-Read:'

      CASE (FILETYPE_GRB2)
        ! Search name mapping for name in GRIB2 file
        mapped_name = TRIM(dict_get(ana_varnames_dict, varname, default=varname))

        WRITE(filetyp_read,'(a)') 'GRB2-Read:'

      CASE DEFAULT
        CALL finish(routine, "Unknown file type")
      END SELECT



      ! Perform CDI read operation
      !
      IF(my_process_is_stdio() .AND. msg_level>10) THEN
        IF ( (tileinfo%idx == trivial_tileinfo%idx) .AND. (tileinfo%att == trivial_tileinfo%att)) THEN
          WRITE(message_text,'(a)') TRIM(filetyp_read)//' '//TRIM(mapped_name)
        ELSE
          WRITE(message_text,'(a,i2,a,i2)') TRIM(filetyp_read)//' '//TRIM(mapped_name)//', Tile:',tileinfo%idx,', Att:',tileinfo%att
        ENDIF
        CALL message(TRIM(routine), TRIM(message_text))
      ENDIF

      CALL read_cdi_2d(parameters, mapped_name, var_out, tileinfo)


      ! add tile suffix (if available) to varname
      ! This is necessary, in order to select the right post-op
      IF ((tileinfo%idx /= trivial_tileinfo%idx) .AND. (tileinfo%att /= trivial_tileinfo%att)) THEN
        full_varname = TRIM(varname//TRIM(get_tile_suffix(tileinfo%idx,tileinfo%att)))
      ELSE
        full_varname = varname
      ENDIF

      ! Perform inverse post_op on input field, if necessary
      !
      CALL initicon_inverse_post_op(full_varname, mapped_name, optvar_out2D=var_out)

    ENDIF  ! lread

  END SUBROUTINE read_data_2d



  !-------------------------------------------------------------------------
  !> Wrapper routine for NetCDF and GRIB2 read routines, 3D case.
  !
  SUBROUTINE read_data_3d(parameters, filetype, varname, nlevs, var_out, tileinfo, opt_checkgroup)
    TYPE(t_inputParameters), INTENT(INOUT) :: parameters
    INTEGER,                 INTENT(IN)    :: filetype       !< FILETYPE_NC2 or FILETYPE_GRB2
    CHARACTER(len=*),        INTENT(IN)    :: varname        !< var name of field to be read
    INTEGER,                 INTENT(IN)    :: nlevs          !< vertical levels of netcdf file
    REAL(wp), POINTER,       INTENT(INOUT) :: var_out(:,:,:) !< output field
    TYPE(t_tileinfo_elt),    INTENT(IN)    :: tileinfo
    CHARACTER(LEN=VARNAME_LEN), INTENT(IN), OPTIONAL :: opt_checkgroup(:) !< read only, if varname is
                                                                          !< contained in opt_checkgroup
    ! local variables
    CHARACTER(len=*), PARAMETER     :: routine = 'read_data_3d'
    CHARACTER(LEN=DICT_MAX_STRLEN)  :: mapped_name
    CHARACTER(LEN=MAX_CHAR_LENGTH)  :: filetyp_read   !< filetype for log message
    LOGICAL                         :: lread          !< .FALSE.: skip reading
    CHARACTER(LEN=MAX_CHAR_LENGTH)  :: full_varname   !< varname, including tile suffix, if present

    IF (PRESENT(opt_checkgroup)) THEN
      lread = ( one_of(varname,  opt_checkgroup(:)) /= -1)
    ELSE
      lread = .TRUE.
    ENDIF

    IF (lread) THEN

      SELECT CASE(filetype)
      CASE (FILETYPE_NC2, FILETYPE_NC4)
        !
        ! Trivial name mapping for NetCDF file
        mapped_name = TRIM(varname)

        WRITE(filetyp_read,'(a)') 'NC-Read:'

      CASE (FILETYPE_GRB2)
        !
        ! Search name mapping for name in GRIB2 file
        mapped_name = TRIM(dict_get(ana_varnames_dict, varname, default=varname))

        WRITE(filetyp_read,'(a)') 'GRB2-Read:'

      CASE DEFAULT
        CALL finish(routine, "Unknown file type")
      END SELECT


      ! Perform CDI read operation
      !
      IF(my_process_is_stdio() .AND. msg_level>10) THEN
        IF ( (tileinfo%idx == trivial_tileinfo%idx) .AND. (tileinfo%att == trivial_tileinfo%att)) THEN
          WRITE(message_text,'(a)') TRIM(filetyp_read)//' '//TRIM(mapped_name)
        ELSE
          WRITE(message_text,'(a,i2,a,i2)') TRIM(filetyp_read)//' '//TRIM(mapped_name)//', Tile:',tileinfo%idx,', Att:',tileinfo%att
        ENDIF
        CALL message(TRIM(routine), TRIM(message_text))
      ENDIF

      CALL read_cdi_3d (parameters, mapped_name, nlevs, var_out, tileinfo)


      ! add tile suffix (if available) to varname
      ! This is necessary, in order to select the right post-op
      IF ((tileinfo%idx /= trivial_tileinfo%idx) .AND. (tileinfo%att /= trivial_tileinfo%att)) THEN
        full_varname = TRIM(varname//TRIM(get_tile_suffix(tileinfo%idx,tileinfo%att)))
      ELSE
        full_varname = varname
      ENDIF

      ! Perform inverse post_op on input field, if necessary
      !
      ! SMI is skipped manually, since it is not contained in any of the ICON variable
      ! lists, and is thus not handled correctly by the following routine.
      IF( TRIM(mapped_name)/='smi' .AND. TRIM(mapped_name)/='SMI') &
        CALL initicon_inverse_post_op(full_varname, mapped_name, optvar_out3D=var_out)

    ENDIF  ! lread

  END SUBROUTINE read_data_3d



  !>
  !! Read horizontally interpolated external analysis (atmosphere only)
  !!
  !! Reads horizontally interpolated external analysis atmosphere data
  !! (currently IFS or COSMO) and reads in vertical coordinate table.
  !!
  !! @par Revision History
  !! Initial version by Guenther Zaengl, DWD(2011-07-14)
  !! Modification by Daniel Reinert, DWD (2012-12-18)
  !! - encapsulate reading of IFS analysis
  !!
  SUBROUTINE read_extana_atm (p_patch, initicon)

    TYPE(t_patch),          INTENT(IN)    :: p_patch(:)
    TYPE(t_initicon_state), INTENT(INOUT) :: initicon(:)

    INTEGER :: jg, jlev, jc, jk, jb, i_endidx
    LOGICAL :: l_exist

    INTEGER :: no_cells, no_levels
    INTEGER :: ncid, dimid, varid, mpi_comm
    TYPE(t_stream_id) :: stream_id
    INTEGER :: psvar_ndims, geopvar_ndims

    CHARACTER(LEN=10) :: psvar

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
      routine = 'mo_nh_initicon:read_extana_atm'

    CHARACTER(LEN=filename_max) :: ifs2icon_file(max_dom)
    LOGICAL :: lread_process
    LOGICAL :: lread_qr, lread_qs ! are qr, qs provided as input?

    !-------------------------------------------------------------------------

    ! flag. if true, then this PE reads data from file and broadcasts
    lread_process = my_process_is_mpi_workroot()

    DO jg = 1, n_dom

      jlev = p_patch(jg)%level


      ! Skip reading the atmospheric input data if a model domain
      ! is not active at initial time
      IF (.NOT. p_patch(jg)%ldom_active) CYCLE

      !
      ! generate file name
      !
      ifs2icon_file(jg) = generate_filename(ifs2icon_filename, model_base_dir, &
        &                                   nroot, jlev, jg)


      ! Read in data from IFS2ICON
      !
      IF( lread_process ) THEN

        INQUIRE (FILE=ifs2icon_file(jg), EXIST=l_exist)
        IF (.NOT.l_exist) THEN
          CALL finish(TRIM(routine),'IFS2ICON file is not found: '//TRIM(ifs2icon_file(jg)))
        ENDIF

        !
        ! open file
        !
        CALL nf(nf_open(TRIM(ifs2icon_file(jg)), NF_NOWRITE, ncid), routine)

        !
        ! get number of cells
        !
        CALL nf(nf_inq_dimid(ncid, 'ncells', dimid), routine)
        CALL nf(nf_inq_dimlen(ncid, dimid, no_cells), routine)

        !
        ! get number of vertical levels
        !
        CALL nf(nf_inq_dimid(ncid, 'lev', dimid), routine)
        CALL nf(nf_inq_dimlen(ncid, dimid, no_levels), routine)

        !
        ! check the number of cells
        !
        IF(p_patch(jg)%n_patch_cells_g /= no_cells) THEN
          CALL finish(TRIM(ROUTINE),&
          & 'Number of patch cells and cells in IFS2ICON file do not match.')
        ENDIF

        !
        ! DEFINE the number of vertical levels
        !
        IF(nlev_in /= 0)  THEN
          IF (nlev_in /= no_levels) THEN
            CALL finish(TRIM(ROUTINE), 'nlev_in has already been defined differently!')
          END IF
        ELSE
          nlev_in = no_levels
        END IF

        !
        ! Check if surface pressure (PS) or its logarithm (LNPS) is provided as input
        !
        IF (nf_inq_varid(ncid, 'PS', varid) == nf_noerr) THEN
          psvar = 'PS'
        ELSE IF (nf_inq_varid(ncid, 'LNPS', varid) == nf_noerr) THEN
          psvar = 'LNPS'
        ENDIF

        !Find out the dimension of psvar for reading purpose
        IF (nf_inq_varid(ncid, TRIM(psvar), varid) == nf_noerr) THEN
          CALL nf(nf_inq_varndims(ncid,varid,psvar_ndims), routine)
        ELSE
          CALL finish(TRIM(routine),'surface pressure var '//TRIM(psvar)//' is missing')
        ENDIF

        !
        ! Check if model-level surface Geopotential is provided as GEOSP or GEOP_ML
        !
        IF (nf_inq_varid(ncid, 'GEOSP', varid) == nf_noerr) THEN
          geop_ml_var = 'GEOSP'
        ELSE IF (nf_inq_varid(ncid, 'GEOP_ML', varid) == nf_noerr) THEN
          geop_ml_var = 'GEOP_ML'
        ELSE
          CALL finish(TRIM(routine),'Could not find model-level sfc geopotential')
        ENDIF

        !Find out the dimension of geop_ml_var for reading purpose
        IF (nf_inq_varid(ncid, TRIM(geop_ml_var), varid) == nf_noerr) THEN
          CALL nf(nf_inq_varndims(ncid,varid,geopvar_ndims), routine)
        ELSE
          CALL finish(TRIM(routine),'surface geopotential var '//TRIM(geop_ml_var)//' is missing')
        ENDIF


        !
        ! Check if rain water (QR) is provided as input
        !
        IF (nf_inq_varid(ncid, 'QR', varid) == nf_noerr) THEN
          lread_qr = .true.
        ELSE
          lread_qr = .false.
          CALL message(TRIM(routine),'Rain water (QR) not available in input data')
        ENDIF

        !
        ! Check if snow water (QS) is provided as input
        !
        IF (nf_inq_varid(ncid, 'QS', varid) == nf_noerr) THEN
          lread_qs = .true.
        ELSE
          lread_qs = .false.
          CALL message(TRIM(routine),'Snow water (QS) not available in input data')
        ENDIF

        !
        ! Check if normal velocity component (VN) is provided as input
        !
        IF (nf_inq_varid(ncid, 'VN', varid) == nf_noerr) THEN
          lread_vn = .TRUE.
        ELSE
          lread_vn = .FALSE.
        ENDIF

      ENDIF ! pe_io

      !
      ! open file
      !
      stream_id = openInputFile(ifs2icon_file(jg), p_patch(jg), &
        &                       default_read_method)

      IF(p_test_run) THEN
        mpi_comm = p_comm_work_test
      ELSE
        mpi_comm = p_comm_work
      ENDIF

      CALL p_bcast(nlev_in,       p_io, mpi_comm)
      CALL p_bcast(lread_qs,      p_io, mpi_comm)
      CALL p_bcast(lread_qr,      p_io, mpi_comm)
      CALL p_bcast(lread_vn,      p_io, mpi_comm)
      CALL p_bcast(psvar_ndims,   p_io, mpi_comm)
      CALL p_bcast(geopvar_ndims, p_io, mpi_comm)

      IF (msg_level >= 10) THEN
        WRITE(message_text,'(a)') 'surface pressure variable: '//TRIM(psvar)
        CALL message(TRIM(routine), TRIM(message_text))
        WRITE(message_text,'(a)') 'Model-level surface geopotential: '//TRIM(geop_ml_var)
        CALL message(TRIM(routine), TRIM(message_text))
        IF (.NOT. lread_vn) THEN
          WRITE(message_text,'(a)') 'No direct input of vn! vn derived from (u,v).'
          CALL message(TRIM(routine), TRIM(message_text))
        ELSE
          WRITE(message_text,'(a)') 'Direct input of vn!'
          CALL message(TRIM(routine), TRIM(message_text))
        ENDIF
      ENDIF

      ! allocate data structure
      CALL allocate_extana_atm(jg, p_patch(jg)%nblks_c, p_patch(jg)%nblks_e, initicon)

      ! start reading atmospheric fields
      !
      CALL read_3d_1time(stream_id, onCells, 'T', fill_array=initicon(jg)%atm_in%temp)

      IF (lread_vn) THEN
        CALL read_3d_1time(stream_id, onEdges, 'VN', fill_array=initicon(jg)%atm_in%vn)
      ELSE
        CALL read_3d_1time(stream_id, onCells, 'U', fill_array=initicon(jg)%atm_in%u)

        CALL read_3d_1time(stream_id, onCells, 'V', fill_array=initicon(jg)%atm_in%v)
      ENDIF


      IF (init_mode == MODE_COSMODE) THEN
        CALL read_3d_1time(stream_id, onCells, 'W', fill_array=initicon(jg)%atm_in%w_ifc)
      ELSE
        ! Note: in this case, input vertical velocity is in fact omega (Pa/s)
        CALL read_3d_1time(stream_id, onCells, 'W', fill_array=initicon(jg)%atm_in%omega)
      ENDIF

      CALL read_3d_1time(stream_id, onCells, 'QV', fill_array=initicon(jg)%atm_in%qv)
      CALL read_3d_1time(stream_id, onCells, 'QC', fill_array=initicon(jg)%atm_in%qc)
      CALL read_3d_1time(stream_id, onCells, 'QI', fill_array=initicon(jg)%atm_in%qi)

      IF (lread_qr) THEN
        CALL read_3d_1time(stream_id, onCells, 'QR', fill_array=initicon(jg)%atm_in%qr)
      ELSE
        initicon(jg)%atm_in%qr(:,:,:)=0._wp
      ENDIF

      IF (lread_qs) THEN
        CALL read_3d_1time(stream_id, onCells, 'QS', fill_array=initicon(jg)%atm_in%qs)
      ELSE
        initicon(jg)%atm_in%qs(:,:,:)=0._wp
      ENDIF

      IF (psvar_ndims==2)THEN
        CALL read_2d_1time(stream_id, onCells, TRIM(psvar), &
          &                     fill_array=initicon(jg)%atm_in%psfc)
      ELSEIF(psvar_ndims==3)THEN
        CALL read_2d_1lev_1time(stream_id, onCells, TRIM(psvar), &
          &                     fill_array=initicon(jg)%atm_in%psfc)
      ELSE
        CALL finish(TRIM(routine),'surface pressure var '//TRIM(psvar)//' dimension mismatch')
      END IF

      IF (geopvar_ndims==2)THEN
        CALL read_2d_1time(stream_id, onCells, TRIM(geop_ml_var), &
          &                     fill_array=initicon(jg)%atm_in%phi_sfc)
      ELSEIF(geopvar_ndims==3)THEN
        CALL read_2d_1lev_1time(stream_id, onCells, TRIM(geop_ml_var), &
          &                     fill_array=initicon(jg)%atm_in%phi_sfc)
      ELSE
        CALL finish(TRIM(routine),'surface geopotential var '//TRIM(geop_ml_var)//' dimension mismatch')
      END IF

      ! Allocate and read in vertical coordinate tables
      !
      ! Note that here the IFS input vertical grid is set up. This has to be distinguished
      ! from vct_a, vct_b, vct for the ICON vertical grid.
      !
      IF (init_mode == MODE_IFSANA .OR. init_mode == MODE_COMBINED) THEN

        IF (jg == 1) THEN

          ALLOCATE(vct_a(nlev_in+1), vct_b(nlev_in+1), vct(2*(nlev_in+1)))

          IF(lread_process) THEN
            CALL nf(nf_inq_varid(ncid, 'hyai', varid), routine)
            CALL nf(nf_get_var_double(ncid, varid, vct_a), routine)

            CALL nf(nf_inq_varid(ncid, 'hybi', varid), routine)
            CALL nf(nf_get_var_double(ncid, varid, vct_b), routine)
          ENDIF

          IF(p_test_run) THEN
            mpi_comm = p_comm_work_test
          ELSE
            mpi_comm = p_comm_work
          ENDIF

          CALL p_bcast(vct_a, p_io, mpi_comm)
          CALL p_bcast(vct_b, p_io, mpi_comm)


          vct(1:nlev_in+1)             = vct_a(:)
          vct(nlev_in+2:2*(nlev_in+1)) = vct_b(:)


          CALL alloc_vct(nlev_in)
          CALL init_vct(nlev_in)

        ENDIF  ! jg=1

      ELSE IF (init_mode == MODE_COSMODE) THEN ! in case of COSMO-DE initial data

        CALL read_3d_1time(stream_id, onCells, 'HHL', fill_array=initicon(jg)%atm_in%z3d_ifc)
        CALL read_3d_1time(stream_id, onCells, 'P', fill_array=initicon(jg)%atm_in%pres)

        ! Interpolate input 'z3d' and 'w' from interface levels to main levels
!$OMP PARALLEL
!$OMP DO PRIVATE (jk,jc,jb,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
        DO jb = 1,p_patch(jg)%nblks_c

          IF (jb /= p_patch(jg)%nblks_c) THEN
            i_endidx = nproma
          ELSE
            i_endidx = p_patch(jg)%npromz_c
          ENDIF

          DO jk = 1, nlev_in
            DO jc = 1, i_endidx

              initicon(jg)%atm_in%z3d(jc,jk,jb) = (initicon(jg)%atm_in%z3d_ifc(jc,jk,jb) + &
                &   initicon(jg)%atm_in%z3d_ifc(jc,jk+1,jb)) * 0.5_wp
              initicon(jg)%atm_in%w(jc,jk,jb) = (initicon(jg)%atm_in%w_ifc(jc,jk,jb) +     &
                &   initicon(jg)%atm_in%w_ifc(jc,jk+1,jb)) * 0.5_wp
            ENDDO
          ENDDO
        ENDDO
!$OMP END DO
!$OMP END PARALLEL
      ELSE

        CALL finish(TRIM(routine),'Incorrect init_mode')

      ENDIF ! init_mode = MODE_COSMODE

      ! close file
      !
      IF(lread_process) CALL nf(nf_close(ncid), routine)
      CALL closeFile(stream_id)

    ENDDO ! loop over model domains

  END SUBROUTINE read_extana_atm





  !>
  !! Read horizontally interpolated external analysis (surface only)
  !!
  !! Reads horizontally interpolated external analysis surface data
  !! Currently, only IFS data are processed here
  !!
  !! @par Revision History
  !! Initial version by Guenther Zaengl, DWD(2011-07-14)
  !! Modification by Daniel Reinert, DWD (2012-12-18)
  !! - encapsulate reading of IFS analysis
  !!
  SUBROUTINE read_extana_sfc (p_patch, initicon)

    TYPE(t_patch),          INTENT(IN)    :: p_patch(:)
    TYPE(t_initicon_state), INTENT(INOUT) :: initicon(:)

    INTEGER :: jg, jlev
    LOGICAL :: l_exist

    INTEGER :: no_cells, no_levels
    INTEGER :: ncid, dimid, varid, mpi_comm
    INTEGER :: geop_sfc_var_ndims     ! dimension of geop_sfc_var
    TYPE(t_stream_id) :: stream_id

    CHARACTER(LEN=10) :: geop_sfc_var ! surface-level surface geopotential

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
      routine = 'mo_nh_initicon:read_extana_sfc'

    CHARACTER(LEN=filename_max) :: ifs2icon_file(max_dom)
    LOGICAL :: lread_process

    !-------------------------------------------------------------------------

    ! flag. if true, then this PE reads data from file and broadcasts
    lread_process = my_process_is_mpi_workroot()

    DO jg = 1, n_dom

      jlev = p_patch(jg)%level


      ! Skip reading the atmospheric input data if a model domain
      ! is not active at initial time
      IF (.NOT. p_patch(jg)%ldom_active) CYCLE

      !
      ! generate file name
      !
      ifs2icon_file(jg) = generate_filename(ifs2icon_filename, model_base_dir, &
        &                                   nroot, jlev, jg)

      ! Read in data from IFS2ICON
      !
      IF(lread_process) THEN
        INQUIRE (FILE=ifs2icon_file(jg), EXIST=l_exist)
        IF (.NOT.l_exist) THEN
          CALL finish(TRIM(routine),'IFS2ICON file is not found: '//TRIM(ifs2icon_file(jg)))
        ENDIF

        !
        ! open file
        !
        CALL nf(nf_open(TRIM(ifs2icon_file(jg)), NF_NOWRITE, ncid), routine)

        !
        ! get number of cells
        !
        CALL nf(nf_inq_dimid(ncid, 'ncells', dimid), routine)
        CALL nf(nf_inq_dimlen(ncid, dimid, no_cells), routine)

        !
        ! get number of vertical levels
        !
        CALL nf(nf_inq_dimid(ncid, 'lev', dimid), routine)
        CALL nf(nf_inq_dimlen(ncid, dimid, no_levels), routine)

        !
        ! check the number of cells and vertical levels
        !
        IF(p_patch(jg)%n_patch_cells_g /= no_cells) THEN
          CALL finish(TRIM(ROUTINE),&
          & 'Number of patch cells and cells in IFS2ICON file do not match.')
        ENDIF

        ! Check, if the surface-level surface geopotential (GEOP_SFC) is available.
        ! If GEOP_SFC is missing, a warning will be issued and the model-level surface
        ! geopotential (GEOSP or GEOP_ML) will be used instead.
        IF (nf_inq_varid(ncid, 'GEOP_SFC', varid) == nf_noerr) THEN
          geop_sfc_var = 'GEOP_SFC'
        ELSE

          WRITE (message_text,'(a,a)')                            &
            &  'surface-level surface geopotential is missing. ', &
            &  'use model-level surface geopotential, instead.'
          CALL message(TRIM(routine),TRIM(message_text))

          ! use model level geopotential instead
          geop_sfc_var = geop_ml_var
        ENDIF
        ! inquire dimension of geop_sfc_var which can be either 2 or 3
        IF (nf_inq_varid(ncid, TRIM(geop_sfc_var), varid) == nf_noerr) THEN
          CALL nf(nf_inq_varndims(ncid,varid,geop_sfc_var_ndims), routine)
        ELSE
          CALL finish(TRIM(ROUTINE),'surface geopotential variable '//TRIM(geop_sfc_var)//' is missing')
        ENDIF

        ! Check, if the snow albedo ALB_SNOW is available.
        ! If ALB_SNOW is missing, a warning will be issued and RHO_SNOW
        ! will be used instead to determine FRESHSNOW.
        IF (nf_inq_varid(ncid, 'ALB_SNOW', varid) == nf_noerr) THEN
          WRITE (message_text,'(a,a)')                            &
            &  'snow albedo available, ', &
            &  'used to determine freshsnow.'
          alb_snow_var = 'ALB_SNOW'
        ELSE

          WRITE (message_text,'(a,a)')                            &
            &  'snow albedo is missing. ', &
            &  'use snow density value, instead.'
          CALL message(TRIM(routine),TRIM(message_text))

          alb_snow_var = 'RHO_SNOW'
        ENDIF


        ! Check, if sea surface temperature field is provided as input
        ! IF SST is missing, set l_sst_in=.FALSE.
        IF (nf_inq_varid(ncid, 'SST', varid) /= nf_noerr) THEN
          WRITE (message_text,'(a,a)')                            &
            &  'sea surface temperature not available. ', &
            &  'initialize with skin temperature, instead.'
          CALL message(TRIM(routine),TRIM(message_text))
          l_sst_in = .FALSE.     !it has to be set to FALSE
        ENDIF

      ENDIF  ! p_io

      !
      ! open file
      !
      stream_id = openInputFile(ifs2icon_file(jg), p_patch(jg), &
        &                       default_read_method)

      IF(p_test_run) THEN
        mpi_comm = p_comm_work_test
      ELSE
        mpi_comm = p_comm_work
      ENDIF

      CALL p_bcast(l_sst_in, p_io, mpi_comm)

      CALL p_bcast(alb_snow_var, p_io, mpi_comm)

      CALL p_bcast(geop_sfc_var_ndims, p_io, mpi_comm)

      ! allocate data structure
      CALL allocate_extana_sfc(jg, p_patch(jg)%nblks_c, initicon)


      ! start reading surface fields
      !
      IF (geop_sfc_var_ndims == 3) THEN
        CALL read_2d_1lev_1time(stream_id, onCells, TRIM(geop_sfc_var), &
          &                     fill_array=initicon(jg)%sfc_in%phi)
      ELSE IF (geop_sfc_var_ndims == 2) THEN
        CALL read_2d_1time(stream_id, onCells, TRIM(geop_sfc_var), &
          &                     fill_array=initicon(jg)%sfc_in%phi)
      ELSE
        CALL finish(TRIM(routine),"geop_sfc_var: Dimension mismatch")
      ENDIF

      CALL read_2d_1time(stream_id, onCells, 'SKT', &
        &                fill_array=initicon(jg)%sfc_in%tskin)

      IF ( l_sst_in) THEN
        CALL read_2d_1time(stream_id, onCells, 'SST', &
          &                fill_array=initicon(jg)%sfc_in%sst)
      ELSE
       initicon(jg)%sfc_in%sst(:,:)=0.0_wp
      END IF

      CALL read_2d_1time(stream_id, onCells, 'T_SNOW', &
        &                fill_array=initicon(jg)%sfc_in%tsnow)
      CALL read_2d_1time(stream_id, onCells, TRIM(alb_snow_var), &
        &                fill_array=initicon(jg)%sfc_in%snowalb)
      CALL read_2d_1time(stream_id, onCells, 'W_SNOW', &
        &                fill_array=initicon(jg)%sfc_in%snowweq)
      CALL read_2d_1time(stream_id, onCells, 'RHO_SNOW', &
        &                fill_array=initicon(jg)%sfc_in%snowdens)
      CALL read_2d_1time(stream_id, onCells, 'W_I', &
        &                fill_array=initicon(jg)%sfc_in%skinres)
      CALL read_2d_1time(stream_id, onCells, 'LSM', &
        &                fill_array=initicon(jg)%sfc_in%ls_mask)
      CALL read_2d_1time(stream_id, onCells, 'CI', &
        &                fill_array=initicon(jg)%sfc_in%seaice)
      CALL read_2d_1lev_1time(stream_id, onCells, 'STL1', &
        &                     fill_array=initicon(jg)%sfc_in%tsoil(:,:,1))
      CALL read_2d_1lev_1time(stream_id, onCells, 'STL2', &
        &                     fill_array=initicon(jg)%sfc_in%tsoil(:,:,2))
      CALL read_2d_1lev_1time(stream_id, onCells, 'STL3', &
        &                     fill_array=initicon(jg)%sfc_in%tsoil(:,:,3))
      CALL read_2d_1lev_1time(stream_id, onCells, 'STL4', &
        &                     fill_array=initicon(jg)%sfc_in%tsoil(:,:,4))
      CALL read_2d_1lev_1time(stream_id, onCells, 'SMIL1', &
        &                     fill_array=initicon(jg)%sfc_in%wsoil(:,:,1))
      CALL read_2d_1lev_1time(stream_id, onCells, 'SMIL2', &
        &                     fill_array=initicon(jg)%sfc_in%wsoil(:,:,2))
      CALL read_2d_1lev_1time(stream_id, onCells, 'SMIL3', &
        &                     fill_array=initicon(jg)%sfc_in%wsoil(:,:,3))
      CALL read_2d_1lev_1time(stream_id, onCells, 'SMIL4', &
        &                     fill_array=initicon(jg)%sfc_in%wsoil(:,:,4))

      ! close file
      !
      IF(lread_process) CALL nf(nf_close(ncid), routine)
      CALL closeFile(stream_id)


    ENDDO ! loop over model domains

  END SUBROUTINE read_extana_sfc



  !>
  !! Read DWD first guess (atmosphere only)
  !!
  !! Read DWD first guess (atmosphere only)
  !! First guess (FG) is read for theta_v, rho, vn, w, tke,
  !! whereas DA output is read for T, p, u, v,
  !! qv, qc, qi, qr, qs.
  !!
  !! @par Revision History
  !! Initial version by Daniel Reinert, DWD(2012-12-18)
  !! Modifications for GRIB2 : F. Prill, DWD (2013-02-19)
  !! Modifications by Daniel Reinert, DWD (2014-01-27)
  !! - split off reading of FG fields
  !!
  !!
  SUBROUTINE read_dwdfg_atm (p_patch, p_nh_state, initicon, fileID_fg, filetype_fg, dwdfg_file)

    TYPE(t_patch),          INTENT(IN)    :: p_patch(:)
    TYPE(t_nh_state),       INTENT(INOUT) :: p_nh_state(:)
    INTEGER,                INTENT(IN)    :: fileID_fg(:), filetype_fg(:)

    TYPE(t_initicon_state), INTENT(INOUT), TARGET :: initicon(:)
    CHARACTER(LEN=filename_max), INTENT(IN)       :: dwdfg_file(:)

    INTEGER :: jg
    INTEGER :: nlev, nlevp1

    INTEGER :: ngrp_vars_fg, filetype

    REAL(wp), POINTER :: my_ptr3d(:,:,:)
    TYPE(t_inputParameters) :: parameters
    CHARACTER(LEN=VARNAME_LEN), POINTER :: checkgrp(:)

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
      routine = modname//':read_dwdfg_atm'
    TYPE(t_tileinfo_elt) :: tileinfo

    !-------------------------------------------------------------------------

    tileinfo = trivial_tileinfo

    DO jg = 1, n_dom

      ! number of vertical full and half levels
      nlev   = p_patch(jg)%nlev
      nlevp1 = p_patch(jg)%nlevp1

      ! Skip reading the atmospheric input data if a model domain
      ! is not active at initial time
      IF (.NOT. p_patch(jg)%ldom_active) CYCLE

      ! save some paperwork
      ngrp_vars_fg  = initicon(jg)%ngrp_vars_fg

      !---------------------------------------!
      ! Read in DWD first guess (atmosphere)  !
      !---------------------------------------!

      IF(my_process_is_stdio() ) THEN
        CALL message (TRIM(routine), 'read atm_FG fields from '//TRIM(dwdfg_file(jg)))
      ENDIF  ! p_io

      parameters = makeInputParameters(fileID_fg(jg), p_patch(jg)%n_patch_cells_g, p_patch(jg)%comm_pat_scatter_c)

      filetype = filetype_fg(jg)
      checkgrp => initicon(jg)%grp_vars_fg(1:ngrp_vars_fg)

      ! start reading first guess (atmosphere only)
      !
      CALL read_data_3d (parameters, filetype, 'theta_v', nlev, p_nh_state(jg)%prog(nnow(jg))%theta_v, tileinfo)
      CALL read_data_3d (parameters, filetype, 'rho', nlev, p_nh_state(jg)%prog(nnow(jg))%rho, tileinfo)
      CALL read_data_3d (parameters, filetype, 'w', nlevp1, p_nh_state(jg)%prog(nnow(jg))%w, tileinfo)
      CALL read_data_3d (parameters, filetype, 'tke', nlevp1, p_nh_state(jg)%prog(nnow(jg))%tke, tileinfo)

      ! Only needed for FG-only runs; usually read from ANA
      my_ptr3d => p_nh_state(jg)%prog(nnow(jg))%tracer(:,:,:,iqv)
      CALL read_data_3d (parameters, filetype, 'qv', nlev, my_ptr3d, tileinfo, opt_checkgroup=checkgrp)

      my_ptr3d => p_nh_state(jg)%prog(nnow(jg))%tracer(:,:,:,iqc)
      CALL read_data_3d (parameters, filetype, 'qc', nlev, my_ptr3d, tileinfo, opt_checkgroup=checkgrp)

      my_ptr3d => p_nh_state(jg)%prog(nnow(jg))%tracer(:,:,:,iqi)
      CALL read_data_3d (parameters, filetype, 'qi', nlev, my_ptr3d, tileinfo, opt_checkgroup=checkgrp)

      IF ( iqr /= 0 ) THEN
        my_ptr3d => p_nh_state(jg)%prog(nnow(jg))%tracer(:,:,:,iqr)
        CALL read_data_3d (parameters, filetype, 'qr', nlev, my_ptr3d, tileinfo, opt_checkgroup=checkgrp)
      END IF

      IF ( iqs /= 0 ) THEN
        my_ptr3d => p_nh_state(jg)%prog(nnow(jg))%tracer(:,:,:,iqs)
        CALL read_data_3d (parameters, filetype, 'qs', nlev, my_ptr3d, tileinfo, opt_checkgroup=checkgrp)
      END IF

      CALL deleteInputParameters(parameters)

      !This call needs its own input parameter object because it's edge based.
      parameters = makeInputParameters(fileID_fg(jg), p_patch(jg)%n_patch_edges_g, p_patch(jg)%comm_pat_scatter_e)
      CALL read_data_3d (parameters, filetype, 'vn', nlev, p_nh_state(jg)%prog(nnow(jg))%vn, tileinfo)
      CALL deleteInputParameters(parameters)

    ENDDO ! loop over model domains

  END SUBROUTINE read_dwdfg_atm


  !>
  !! Read DWD first guess (atmosphere only) and store to initicon input state
  !! First guess (FG) is read for theta_v, rho, vn, w, tke,
  !! whereas DA output is read for T, p, u, v, 
  !! qv, qc, qi, qr, qs.
  !!
  !! @par Revision History
  !! Initial version by Daniel Reinert, DWD(2012-12-18)
  !! Modifications for GRIB2 : F. Prill, DWD (2013-02-19)
  !! Modifications by Daniel Reinert, DWD (2014-01-27)
  !! - split off reading of FG fields
  !! 
  !!
  SUBROUTINE read_dwdfg_atm_ii (p_patch, initicon, fileID_fg, filetype_fg, dwdfg_file)

    TYPE(t_patch),          INTENT(IN)    :: p_patch(:)
    INTEGER,                INTENT(IN)    :: fileID_fg(:), filetype_fg(:)

    TYPE(t_initicon_state), INTENT(INOUT), TARGET :: initicon(:)
    CHARACTER(LEN=filename_max), INTENT(IN)       :: dwdfg_file(:)

    INTEGER :: jg, jk, jc, jb, i_endidx
    INTEGER :: mpi_comm
    INTEGER :: ngrp_vars_fg, filetype

    REAL(wp) :: tempv, exner
    TYPE(t_inputParameters) :: parameters
    CHARACTER(LEN=VARNAME_LEN), POINTER :: checkgrp(:)

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
      routine = modname//':read_dwdfg_atm_ii'
    LOGICAL :: lread_process
    TYPE(t_tileinfo_elt) :: tileinfo

    !-------------------------------------------------------------------------

    tileinfo = trivial_tileinfo

    ! flag. if true, then this PE reads data from file and broadcasts
    lread_process = my_process_is_mpi_workroot()

    DO jg = 1, n_dom

      ! Skip reading the atmospheric input data if a model domain 
      ! is not active at initial time
      IF (.NOT. p_patch(jg)%ldom_active) CYCLE


      ! determine number of HALF LEVELS of generalized Z-AXIS
      ! use 'theta_v' for testing (may have chosen another one as well)
      IF(lread_process) THEN 
        IF (filetype_fg(jg) == FILETYPE_GRB2) THEN
          nlev_in = get_cdi_NlevRef(fileID_fg(jg), 'theta_v', opt_dict=ana_varnames_dict)
        ELSE
          nlev_in = get_cdi_NlevRef(fileID_fg(jg), 'theta_v')
        ENDIF
        ! number of LAYERS:
        nlev_in = nlev_in - 1
      ENDIF  ! p_io

      IF(p_test_run) THEN
        mpi_comm = p_comm_work_test
      ELSE
        mpi_comm = p_comm_work
      ENDIF
      CALL p_bcast(nlev_in, p_io, mpi_comm)

      ! allocate data structure for reading input data
      CALL allocate_extana_atm(jg, p_patch(jg)%nblks_c, p_patch(jg)%nblks_e, initicon)


      ! save some paperwork
      ngrp_vars_fg  = initicon(jg)%ngrp_vars_fg

      !---------------------------------------!
      ! Read in DWD first guess (atmosphere)  !
      !---------------------------------------!

      IF(my_process_is_stdio() ) THEN 
        CALL message (TRIM(routine), 'read atm_FG fields from '//TRIM(dwdfg_file(jg)))
      ENDIF  ! p_io

      parameters = makeInputParameters(fileID_fg(jg), p_patch(jg)%n_patch_cells_g, p_patch(jg)%comm_pat_scatter_c)

      filetype = filetype_fg(jg)
      checkgrp => initicon(jg)%grp_vars_fg(1:ngrp_vars_fg)

      ! start reading first guess (atmosphere only)
      !
      CALL read_data_3d (parameters, filetype, 'z_ifc', nlev_in+1, initicon(jg)%atm_in%z3d_ifc, tileinfo)

      CALL read_data_3d (parameters, filetype, 'theta_v', nlev_in, initicon(jg)%atm_in%theta_v, tileinfo)
      CALL read_data_3d (parameters, filetype, 'rho', nlev_in, initicon(jg)%atm_in%rho, tileinfo)
      CALL read_data_3d (parameters, filetype, 'w', nlev_in+1, initicon(jg)%atm_in%w_ifc, tileinfo)
      CALL read_data_3d (parameters, filetype, 'tke', nlev_in+1, initicon(jg)%atm_in%tke_ifc, tileinfo)

      CALL read_data_3d (parameters, filetype, 'qv', nlev_in, initicon(jg)%atm_in%qv, tileinfo)

      CALL read_data_3d (parameters, filetype, 'qc', nlev_in, initicon(jg)%atm_in%qc, tileinfo)

      CALL read_data_3d (parameters, filetype, 'qi', nlev_in, initicon(jg)%atm_in%qi, tileinfo)

      IF ( iqr /= 0 ) THEN
        CALL read_data_3d (parameters, filetype, 'qr', nlev_in, initicon(jg)%atm_in%qr, tileinfo)
      ELSE
        initicon(jg)%atm_in%qr(:,:,:) = 0._wp
      END IF

      IF ( iqs /= 0 ) THEN
        CALL read_data_3d (parameters, filetype, 'qs', nlev_in, initicon(jg)%atm_in%qs, tileinfo)
      ELSE
        initicon(jg)%atm_in%qs(:,:,:) = 0._wp
      END IF

      CALL deleteInputParameters(parameters)

      !This call needs its own input parameter object because it's edge based.
      parameters = makeInputParameters(fileID_fg(jg), p_patch(jg)%n_patch_edges_g, p_patch(jg)%comm_pat_scatter_e)
      CALL read_data_3d (parameters, filetype, 'vn', nlev_in, initicon(jg)%atm_in%vn, tileinfo)
      CALL deleteInputParameters(parameters)

      ! Interpolate half level variables from interface levels to main levels, and convert thermodynamic variables
      ! into temperature and pressure as expected by the vertical interpolation routine  
!$OMP PARALLEL
!$OMP DO PRIVATE (jk,jc,jb,i_endidx,tempv,exner) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = 1,p_patch(jg)%nblks_c
        
        IF (jb /= p_patch(jg)%nblks_c) THEN
          i_endidx = nproma
        ELSE
          i_endidx = p_patch(jg)%npromz_c
        ENDIF

        DO jk = 1, nlev_in
          DO jc = 1, i_endidx

            initicon(jg)%atm_in%z3d(jc,jk,jb) = (initicon(jg)%atm_in%z3d_ifc(jc,jk,jb) + &
              &   initicon(jg)%atm_in%z3d_ifc(jc,jk+1,jb)) * 0.5_wp
            initicon(jg)%atm_in%w(jc,jk,jb) = (initicon(jg)%atm_in%w_ifc(jc,jk,jb) +     &
              &   initicon(jg)%atm_in%w_ifc(jc,jk+1,jb)) * 0.5_wp
            initicon(jg)%atm_in%tke(jc,jk,jb) = (initicon(jg)%atm_in%tke_ifc(jc,jk,jb) + &
              &   initicon(jg)%atm_in%tke_ifc(jc,jk+1,jb)) * 0.5_wp

            exner = (initicon(jg)%atm_in%rho(jc,jk,jb)*initicon(jg)%atm_in%theta_v(jc,jk,jb)*rd/p0ref)**(1._wp/cvd_o_rd)
            tempv = initicon(jg)%atm_in%theta_v(jc,jk,jb)*exner

            initicon(jg)%atm_in%pres(jc,jk,jb) = exner**(cpd/rd)*p0ref
            initicon(jg)%atm_in%temp(jc,jk,jb) = tempv / (1._wp + vtmpc1*initicon(jg)%atm_in%qv(jc,jk,jb) - &
              (initicon(jg)%atm_in%qc(jc,jk,jb) + initicon(jg)%atm_in%qi(jc,jk,jb) +                        &
               initicon(jg)%atm_in%qr(jc,jk,jb) + initicon(jg)%atm_in%qs(jc,jk,jb)) )

          ENDDO
        ENDDO

      ENDDO
!$OMP END DO
!$OMP END PARALLEL

    ENDDO ! loop over model domains

    ! Tell the vertical interpolation routine that vn needs to be processed
    lread_vn = .TRUE.

  END SUBROUTINE read_dwdfg_atm_ii



  !>
  !! Read DA-analysis fields (atmosphere only)
  !!
  !! Depending on the initialization mode, either ful fields or increments
  !! are read (atmosphere only):
  !! MODE_DWDANA: The following full fields are read, if available
  !!              u, v, t, p, qv
  !! MODE_DWDANA_INC: The following increments are read, if available
  !!              u, v, t, p, qv
  !!
  !! @par Revision History
  !! Initial version by Daniel Reinert, DWD(2012-12-18)
  !! Modifications for GRIB2 : F. Prill, DWD (2013-02-19)
  !! Modifications by Daniel Reinert, DWD (2014-01-27)
  !! - split off reading of FG fields
  !!
  SUBROUTINE read_dwdana_atm (p_patch, p_nh_state, initicon, fileID_ana,  filetype_ana, dwdana_file)

    TYPE(t_patch),          INTENT(IN)    :: p_patch(:)
    TYPE(t_nh_state),       INTENT(INOUT) :: p_nh_state(:)
    INTEGER,                INTENT(IN)    :: fileID_ana(:), filetype_ana(:)

    TYPE(t_initicon_state), INTENT(INOUT), TARGET :: initicon(:)
    CHARACTER(LEN=filename_max), INTENT(IN)       :: dwdana_file(:)

    INTEGER :: jg, jgp
    INTEGER :: nlev

    INTEGER :: ngrp_vars_ana, filetype

    TYPE(t_pi_atm), POINTER :: my_ptr
    REAL(wp),       POINTER :: my_ptr3d(:,:,:)
    TYPE(t_inputParameters) :: parameters
    CHARACTER(LEN=VARNAME_LEN), POINTER :: checkgrp(:)

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
      routine = modname//':read_dwdana_atm'
    TYPE(t_tileinfo_elt) :: tileinfo

    !-------------------------------------------------------------------------

    tileinfo = trivial_tileinfo

    !----------------------------------------!
    ! read in DWD analysis (atmosphere)      !
    !----------------------------------------!

    DO jg = 1, n_dom

      ! Skip reading the atmospheric input data if a model domain
      ! is not active at initial time
      IF (.NOT. p_patch(jg)%ldom_active) CYCLE

      ! number of vertical full and half levels
      nlev   = p_patch(jg)%nlev

      ! save some paperwork
      ngrp_vars_ana = initicon(jg)%ngrp_vars_ana

      ! Depending on the initialization mode chosen (incremental vs. non-incremental)
      ! input fields are stored in different locations.
      IF ( ANY((/MODE_DWDANA_INC,MODE_IAU,MODE_IAU_OLD/) == init_mode) ) THEN
        IF (lp2cintp_incr(jg)) THEN
          ! Perform parent-to-child interpolation of atmospheric DA increments
          jgp = p_patch(jg)%parent_id
          CALL interpolate_increments(initicon, jgp, jg)
          CYCLE
        ELSE
          my_ptr => initicon(jg)%atm_inc
        ENDIF
      ELSE
        my_ptr => initicon(jg)%atm
      ENDIF

      IF( lread_ana .AND. my_process_is_stdio() ) THEN
        CALL message (TRIM(routine), 'read atm_ANA fields from '//TRIM(dwdana_file(jg)))
      ENDIF  ! p_io

      parameters = makeInputParameters(fileID_ana(jg), p_patch(jg)%n_patch_cells_g, p_patch(jg)%comm_pat_scatter_c)
      filetype = filetype_ana(jg)
      checkgrp => initicon(jg)%grp_vars_ana(1:ngrp_vars_ana)

      ! start reading DA output (atmosphere only)
      ! The dynamical variables temp, pres, u and v, which need further processing,
      ! are either stored in initicon(jg)%atm or initicon(jg)%atm_inc, depending on whether
      ! IAU is used or not. The moisture variables, which can be taken over directly from
      ! the Analysis, are written to the NH prognostic state
      !
      my_ptr3d => my_ptr%temp
      CALL read_data_3d (parameters, filetype, 'temp', nlev, my_ptr3d, tileinfo, opt_checkgroup=checkgrp)

      my_ptr3d => my_ptr%pres
      CALL read_data_3d (parameters, filetype, 'pres', nlev, my_ptr3d, tileinfo, opt_checkgroup=checkgrp )

      my_ptr3d => my_ptr%u
      CALL read_data_3d (parameters, filetype, 'u', nlev, my_ptr3d, tileinfo, opt_checkgroup=checkgrp )

      my_ptr3d => my_ptr%v
      CALL read_data_3d (parameters, filetype, 'v', nlev, my_ptr3d, tileinfo, opt_checkgroup=checkgrp )

      IF ( ANY((/MODE_DWDANA_INC,MODE_IAU,MODE_IAU_OLD/) == init_mode) ) THEN
        my_ptr3d => my_ptr%qv
      ELSE
        my_ptr3d => p_nh_state(jg)%prog(nnow(jg))%tracer(:,:,:,iqv)
      ENDIF

      CALL read_data_3d (parameters, filetype, 'qv', nlev, my_ptr3d, tileinfo, opt_checkgroup=checkgrp )

      ! For the time being identical to qc from FG => usually read from FG
      my_ptr3d => p_nh_state(jg)%prog(nnow(jg))%tracer(:,:,:,iqc)
      CALL read_data_3d (parameters, filetype, 'qc', nlev, my_ptr3d, tileinfo, opt_checkgroup=checkgrp )

      ! For the time being identical to qi from FG => usually read from FG
      my_ptr3d => p_nh_state(jg)%prog(nnow(jg))%tracer(:,:,:,iqi)
      CALL read_data_3d (parameters, filetype, 'qi', nlev, my_ptr3d, tileinfo, opt_checkgroup=checkgrp )

      ! For the time being identical to qr from FG => usually read from FG
      IF ( iqr /= 0 ) THEN
        my_ptr3d => p_nh_state(jg)%prog(nnow(jg))%tracer(:,:,:,iqr)
        CALL read_data_3d (parameters, filetype, 'qr', nlev, my_ptr3d, tileinfo, opt_checkgroup=checkgrp )
      END IF

      ! For the time being identical to qs from FG => usually read from FG
      IF ( iqs /= 0 ) THEN
        my_ptr3d => p_nh_state(jg)%prog(nnow(jg))%tracer(:,:,:,iqs)
        CALL read_data_3d (parameters, filetype, 'qs', nlev, my_ptr3d, tileinfo, opt_checkgroup=checkgrp )
      END IF

      CALL deleteInputParameters(parameters)

    ENDDO ! loop over model domains

  END SUBROUTINE read_dwdana_atm



  !>
  !! Read DWD first guess (land/surface only)
  !!
  !! Read DWD first guess for land/surface.
  !! First guess is read for:
  !! fr_seaice, t_ice, h_ice, t_g, qv_s, freshsnow, w_snow, w_i, t_snow,
  !! rho_snow, w_so, w_so_ice, t_so, gz0
  !!
  !!
  !! @par Revision History
  !! Initial version by Daniel Reinert, DWD(2012-12-18)
  !! Modifications by Daniel Reinert, DWD (2014-07-16)
  !! - split off reading of FG fields
  !!
  SUBROUTINE read_dwdfg_sfc (p_patch, prm_diag, p_lnd_state, initicon, fileID_fg, filetype_fg, dwdfg_file)

    TYPE(t_patch),             INTENT(IN)    :: p_patch(:)
    TYPE(t_nwp_phy_diag),      INTENT(INOUT) :: prm_diag(:)
    TYPE(t_lnd_state), TARGET, INTENT(INOUT) :: p_lnd_state(:)
    INTEGER,                   INTENT(IN)    :: fileID_fg(:), filetype_fg(:)

    TYPE(t_initicon_state), INTENT(INOUT), TARGET :: initicon(:)
    CHARACTER(LEN=filename_max), INTENT(IN)       :: dwdfg_file(:)

    INTEGER :: jg, jt, jb, jc, i_endidx, filetype

    INTEGER :: ngrp_vars_fg
    REAL(wp), POINTER :: my_ptr2d(:,:)
    REAL(wp), POINTER :: my_ptr3d(:,:,:)
    TYPE(t_lnd_prog), POINTER :: lnd_prog
    TYPE(t_lnd_diag), POINTER :: lnd_diag
    TYPE(t_wtr_prog), POINTER :: wtr_prog
    CHARACTER(LEN=VARNAME_LEN), POINTER :: checkgrp(:)
    TYPE(t_inputParameters) :: parameters

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
      routine = modname//':read_dwdfg_sfc'

    TYPE(t_tileinfo_elt) :: tileinfo

    !-------------------------------------------------------------------------

    tileinfo = trivial_tileinfo

    !----------------------------------------!
    ! read in DWD First Guess (surface)      !
    !----------------------------------------!

    DO jg = 1, n_dom

      lnd_prog => p_lnd_state(jg)%prog_lnd(nnow_rcf(jg))
      lnd_diag => p_lnd_state(jg)%diag_lnd
      wtr_prog => p_lnd_state(jg)%prog_wtr(nnow_rcf(jg))

      ! Skip reading the atmospheric input data if a model domain
      ! is not active at initial time
      IF (.NOT. p_patch(jg)%ldom_active) CYCLE

      ! save some paperwork
      ngrp_vars_fg  = initicon(jg)%ngrp_vars_fg


      IF(my_process_is_stdio() ) THEN
        CALL message (TRIM(routine), 'read sfc_FG fields from '//TRIM(dwdfg_file(jg)))
      ENDIF  ! p_io


      parameters = makeInputParameters(fileID_fg(jg), p_patch(jg)%n_patch_cells_g, p_patch(jg)%comm_pat_scatter_c)
      filetype = filetype_fg(jg)
      checkgrp => initicon(jg)%grp_vars_fg(1:ngrp_vars_fg)

      ! COSMO-DE does not provide sea ice field. In that case set fr_seaice to 0
      IF (init_mode /= MODE_COSMODE) THEN
        CALL read_data_2d(parameters, filetype, 'fr_seaice', lnd_diag%fr_seaice, tileinfo, opt_checkgroup=checkgrp)
      ELSE
!$OMP PARALLEL WORKSHARE
        lnd_diag%fr_seaice(:,:) = 0._wp
!$OMP END PARALLEL WORKSHARE
      ENDIF ! init_mode /= MODE_COSMODE


      ! sea-ice related fields
      CALL read_data_2d(parameters, filetype, 't_ice', wtr_prog%t_ice, tileinfo, opt_checkgroup=checkgrp)
      CALL read_data_2d(parameters, filetype, 'h_ice', wtr_prog%h_ice, tileinfo, opt_checkgroup=checkgrp)

      ! tile based fields
      DO jt=1, ntiles_total + ntiles_water

        IF ( .NOT. ltile_coldstart ) THEN
          tileinfo%idx = tiles(jt)%GRIB2_tile%tileIndex
          tileinfo%att = tiles(jt)%GRIB2_att%tileAttribute
        ENDIF

        my_ptr2d => lnd_prog%t_g_t(:,:,jt)
        CALL read_data_2d(parameters, filetype, 't_g', my_ptr2d, tileinfo, opt_checkgroup=checkgrp)

        my_ptr2d =>lnd_diag%qv_s_t(:,:,jt)
        CALL read_data_2d(parameters, filetype, 'qv_s', my_ptr2d, tileinfo, opt_checkgroup=checkgrp)

        ! Copy t_g_t and qv_s_t to t_g and qv_s, respectively, on limited-area domains.
        ! This is needed to avoid invalid operations in the initialization of the turbulence scheme
        IF (jt == 1 .AND. (jg > 1 .OR. l_limited_area) ) THEN
!$OMP PARALLEL DO PRIVATE(jb,jc,i_endidx)
        ! include boundary interpolation zone of nested domains and halo points
          DO jb = 1, p_patch(jg)%nblks_c
            IF (jb == p_patch(jg)%nblks_c) THEN
              i_endidx = p_patch(jg)%npromz_c
            ELSE
              i_endidx = nproma
            ENDIF

            DO jc = 1, i_endidx
              lnd_prog%t_g(jc,jb)  = lnd_prog%t_g_t(jc,jb,1)
              lnd_diag%qv_s(jc,jb) = lnd_diag%qv_s_t(jc,jb,1)
            ENDDO  ! jc

          ENDDO  ! jb
!$OMP END PARALLEL DO
        ENDIF

      END DO

      !  tile based fields
      DO jt=1, ntiles_total

        IF (.NOT. ltile_coldstart) THEN
          tileinfo%idx = tiles(jt)%GRIB2_tile%tileIndex
          tileinfo%att = tiles(jt)%GRIB2_att%tileAttribute
        ENDIF

        my_ptr2d => lnd_diag%freshsnow_t(:,:,jt)
        CALL read_data_2d(parameters, filetype, 'freshsnow', my_ptr2d, tileinfo, opt_checkgroup=checkgrp )

        my_ptr2d => lnd_diag%snowfrac_t(:,:,jt)
        CALL read_data_2d(parameters, filetype, 'snowfrac', my_ptr2d, tileinfo, opt_checkgroup=checkgrp )

        my_ptr2d => lnd_prog%w_snow_t(:,:,jt)
        CALL read_data_2d(parameters, filetype, 'w_snow', my_ptr2d, tileinfo, opt_checkgroup=checkgrp )

        my_ptr2d => lnd_prog%w_i_t(:,:,jt)
        CALL read_data_2d(parameters, filetype, 'w_i', my_ptr2d, tileinfo, opt_checkgroup=checkgrp )

        my_ptr2d => lnd_diag%h_snow_t(:,:,jt)
        CALL read_data_2d(parameters, filetype, 'h_snow', my_ptr2d, tileinfo, opt_checkgroup=checkgrp )

        my_ptr2d => lnd_prog%t_snow_t(:,:,jt)
        CALL read_data_2d(parameters, filetype, 't_snow', my_ptr2d, tileinfo, opt_checkgroup=checkgrp )

        my_ptr2d => lnd_prog%rho_snow_t(:,:,jt)
        CALL read_data_2d(parameters, filetype, 'rho_snow', my_ptr2d, tileinfo, opt_checkgroup=checkgrp )


        IF (lmulti_snow) THEN
        ! multi layer snow fields
           my_ptr3d => lnd_prog%t_snow_mult_t(:,:,:,jt)
           CALL read_data_3d (parameters, filetype, 't_snow_mult', nlev_snow+1, my_ptr3d, tileinfo, opt_checkgroup=checkgrp )

           my_ptr3d => lnd_prog%rho_snow_mult_t(:,:,:,jt)
           CALL read_data_3d (parameters, filetype, 'rho_snow_mult', nlev_snow, my_ptr3d, tileinfo, opt_checkgroup=checkgrp )

           my_ptr3d => lnd_prog%wtot_snow_t(:,:,:,jt)
           CALL read_data_3d (parameters, filetype, 'wtot_snow', nlev_snow, my_ptr3d, tileinfo, opt_checkgroup=checkgrp )

           my_ptr3d => lnd_prog%wliq_snow_t(:,:,:,jt)
           CALL read_data_3d (parameters, filetype, 'wliq_snow', nlev_snow, my_ptr3d, tileinfo, opt_checkgroup=checkgrp )

           my_ptr3d => lnd_prog%dzh_snow_t(:,:,:,jt)
           CALL read_data_3d (parameters, filetype, 'dzh_snow', nlev_snow, my_ptr3d, tileinfo, opt_checkgroup=checkgrp )
        END IF ! lmulti_snow


        ! multi layer fields
        !
        ! Note that either w_so OR smi is written to w_so_t. Which one is required depends
        ! on the initialization mode. Checking grp_vars_fg takes care of this. In case
        ! that smi is read, it is lateron converted to w_so (see smi_to_wsoil)
        my_ptr3d => lnd_prog%w_so_t(:,:,:,jt)
        CALL read_data_3d (parameters, filetype, 'w_so', nlev_soil, my_ptr3d, tileinfo, opt_checkgroup=checkgrp )
        !
        ! Note that no pointer assignment is missing here, since either W_SO or SMI is read
        ! into w_so_t. Whether SMI or W_SO must be read, is taken care of by 'checkgrp'.
        CALL read_data_3d (parameters, filetype, 'smi', nlev_soil, my_ptr3d, tileinfo, opt_checkgroup=checkgrp )


        ! Skipped in MODE_COMBINED and in MODE_COSMODE. In that case, w_so_ice
        ! is re-diagnosed in terra_multlay_init
        my_ptr3d => lnd_prog%w_so_ice_t(:,:,:,jt)
        CALL read_data_3d (parameters, filetype, 'w_so_ice', nlev_soil, my_ptr3d, tileinfo, opt_checkgroup=checkgrp )

        my_ptr3d => lnd_prog%t_so_t(:,:,:,jt)
        CALL read_data_3d (parameters, filetype, 't_so', nlev_soil+1, my_ptr3d, tileinfo, opt_checkgroup=checkgrp)

      ENDDO ! jt


      tileinfo = trivial_tileinfo

      ! Skipped in MODE_COMBINED and in MODE_COSMODE (i.e. when starting from GME soil)
      ! Instead z0 is re-initialized (see mo_nwp_phy_init)
      CALL read_data_2d(parameters, filetype, 'gz0', prm_diag(jg)%gz0, tileinfo, opt_checkgroup=checkgrp )


      ! first guess for fresh water lake fields
      !
      CALL read_data_2d(parameters, filetype, 't_mnw_lk', wtr_prog%t_mnw_lk, tileinfo, opt_checkgroup=checkgrp)
      CALL read_data_2d(parameters, filetype, 't_wml_lk', wtr_prog%t_wml_lk, tileinfo, opt_checkgroup=checkgrp)
      CALL read_data_2d(parameters, filetype, 'h_ml_lk', wtr_prog%h_ml_lk, tileinfo, opt_checkgroup=checkgrp)
      CALL read_data_2d(parameters, filetype, 't_bot_lk', wtr_prog%t_bot_lk, tileinfo, opt_checkgroup=checkgrp)
      CALL read_data_2d(parameters, filetype, 'c_t_lk', wtr_prog%c_t_lk, tileinfo, opt_checkgroup=checkgrp)
      CALL read_data_2d(parameters, filetype, 't_b1_lk', wtr_prog%t_b1_lk, tileinfo, opt_checkgroup=checkgrp)
      CALL read_data_2d(parameters, filetype, 'h_b1_lk', wtr_prog%h_b1_lk, tileinfo, opt_checkgroup=checkgrp)

      CALL deleteInputParameters(parameters)

      ! Only required, when starting from GME or COSMO soil (i.e. MODE_COMBINED or MODE_COSMODE).
      ! SMI stored in w_so_t must be converted to w_so
      IF (ANY((/MODE_COMBINED,MODE_COSMODE/) == init_mode)) THEN
        DO jt=1, ntiles_total
          CALL smi_to_wsoil(p_patch(jg), lnd_prog%w_so_t(:,:,:,jt))
        ENDDO
      END IF

    ENDDO ! loop over model domains

  END SUBROUTINE read_dwdfg_sfc



  !>
  !! Read DWD analysis (land/surface only)
  !!
  !! Read DWD analysis for land/surface.
  !! Analysis is read for:
  !! fr_seaice, t_ice, h_ice, t_g, qv_s, freshsnow, w_snow, w_i, t_snow,
  !! rho_snow, w_so, w_so_ice, t_so, gz0
  !!
  !!
  !! @par Revision History
  !! Initial version by Daniel Reinert, DWD(2012-12-18)
  !! Modifications by Daniel Reinert, DWD (2014-07-16)
  !! - split off reading of ANA fields
  !!
  SUBROUTINE read_dwdana_sfc (p_patch, p_lnd_state, initicon, fileID_ana, filetype_ana, dwdana_file)

    TYPE(t_patch),             INTENT(IN)    :: p_patch(:)
    TYPE(t_lnd_state), TARGET, INTENT(INOUT) :: p_lnd_state(:)
    INTEGER,                   INTENT(IN)    :: fileID_ana(:), filetype_ana(:)

    TYPE(t_initicon_state), INTENT(INOUT), TARGET :: initicon(:)
    CHARACTER(LEN=filename_max), INTENT(IN)       :: dwdana_file(:)

    INTEGER :: jg, jgp, jt, jb, jc, i_endidx
    INTEGER :: ngrp_vars_ana, filetype

    REAL(wp), POINTER :: my_ptr2d(:,:)
    REAL(wp), POINTER :: my_ptr3d(:,:,:)
    TYPE(t_lnd_prog), POINTER :: lnd_prog
    TYPE(t_lnd_diag), POINTER :: lnd_diag
    TYPE(t_wtr_prog), POINTER :: wtr_prog
    TYPE(t_inputParameters) :: parameters
    CHARACTER(LEN=VARNAME_LEN), POINTER :: checkgrp(:)

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
      routine = modname//':read_dwdana_sfc'
    TYPE(t_tileinfo_elt) :: tileinfo

    !-------------------------------------------------------------------------

    tileinfo = trivial_tileinfo

    !----------------------------------------!
    ! read in DWD analysis (surface)         !
    !----------------------------------------!

    DO jg = 1, n_dom

      lnd_prog => p_lnd_state(jg)%prog_lnd(nnow_rcf(jg))
      lnd_diag => p_lnd_state(jg)%diag_lnd
      wtr_prog => p_lnd_state(jg)%prog_wtr(nnow_rcf(jg))

      ! Skip reading the atmospheric input data if a model domain
      ! is not active at initial time
      IF (.NOT. p_patch(jg)%ldom_active) CYCLE

      ! save some paperwork
      ngrp_vars_ana = initicon(jg)%ngrp_vars_ana

      IF (ANY((/MODE_IAU, MODE_IAU_OLD /) == init_mode) .AND. lp2cintp_sfcana(jg)) THEN
        ! Perform parent-to-child interpolation of surface fields read from the analysis
        jgp = p_patch(jg)%parent_id
        CALL interpolate_sfcana(initicon, jgp, jg)
        CYCLE
      ENDIF

      IF(my_process_is_stdio()) THEN
        CALL message (TRIM(routine), 'read sfc_ANA fields from '//TRIM(dwdana_file(jg)))
      ENDIF   ! p_io


      ! set tile-index explicitly
      jt = 1
      parameters = makeInputParameters(fileID_ana(jg), p_patch(jg)%n_patch_cells_g, p_patch(jg)%comm_pat_scatter_c)
      filetype = filetype_ana(jg)
      checkgrp => initicon(jg)%grp_vars_ana(1:ngrp_vars_ana)

      ! sea-ice fraction
      CALL read_data_2d (parameters, filetype, 'fr_seaice', lnd_diag%fr_seaice, tileinfo, opt_checkgroup=checkgrp )
      ! sea-ice temperature
      CALL read_data_2d (parameters, filetype, 't_ice', wtr_prog%t_ice, tileinfo, opt_checkgroup=checkgrp )
      ! sea-ice height
      CALL read_data_2d (parameters, filetype, 'h_ice', wtr_prog%h_ice, tileinfo, opt_checkgroup=checkgrp )

      ! T_SO(0)
      my_ptr2d => initicon(jg)%sfc%sst(:,:)
      CALL read_data_2d (parameters, filetype, 't_so', my_ptr2d, tileinfo, opt_checkgroup=checkgrp )

      ! h_snow
      IF ( init_mode == MODE_IAU ) THEN
        my_ptr2d => initicon(jg)%sfc_inc%h_snow(:,:)
      ELSE
        my_ptr2d => lnd_diag%h_snow_t(:,:,jt)
      ENDIF
      CALL read_data_2d (parameters, filetype, 'h_snow', my_ptr2d, tileinfo, opt_checkgroup=checkgrp )

      ! w_snow
      my_ptr2d => lnd_prog%w_snow_t(:,:,jt)
      CALL read_data_2d (parameters, filetype, 'w_snow', my_ptr2d, tileinfo, opt_checkgroup=checkgrp )

      ! w_i
      my_ptr2d => lnd_prog%w_i_t(:,:,jt)
      CALL read_data_2d (parameters, filetype, 'w_i', my_ptr2d, tileinfo, opt_checkgroup=checkgrp )

      ! t_snow
      my_ptr2d => lnd_prog%t_snow_t(:,:,jt)
      CALL read_data_2d (parameters, filetype, 't_snow', my_ptr2d, tileinfo, opt_checkgroup=checkgrp )

      ! rho_snow
      my_ptr2d => lnd_prog%rho_snow_t(:,:,jt)
      CALL read_data_2d (parameters, filetype, 'rho_snow', my_ptr2d, tileinfo, opt_checkgroup=checkgrp )

      ! freshsnow
      IF ( init_mode == MODE_IAU ) THEN
        my_ptr2d => initicon(jg)%sfc_inc%freshsnow(:,:)
      ELSE
        my_ptr2d => lnd_diag%freshsnow_t(:,:,jt)
      ENDIF
      CALL read_data_2d (parameters, filetype, 'freshsnow', my_ptr2d, tileinfo, opt_checkgroup=checkgrp )

      ! w_so
      IF ( (init_mode == MODE_IAU) .OR. (init_mode == MODE_IAU_OLD) ) THEN
        my_ptr3d => initicon(jg)%sfc_inc%w_so(:,:,:)
      ELSE
        my_ptr3d => lnd_prog%w_so_t(:,:,:,jt)
      ENDIF
      CALL read_data_3d (parameters, filetype, 'w_so', nlev_soil, my_ptr3d, tileinfo, opt_checkgroup=checkgrp )

      CALL deleteInputParameters(parameters)

      ! MODE_IAU_OLD: H_SNOW, FRESHSNOW, W_SNOW and RHO_SNOW are read from analysis (full fields)
      ! Since only the first tile index is filled (see above), we fill (here) the remaining tiles 
      ! if tile approach is used. Note that for MODE_IAU this copy is skipped on purpose, 
      ! since for ltile_coldstart=.FALSE. tile information would be overwritten. 
      !
      IF (ANY((/ MODE_IAU_OLD /) == init_mode) .AND. ntiles_total>1) THEN

!$OMP PARALLEL DO PRIVATE(jb,jc,jt,i_endidx)
        DO jb = 1, p_patch(jg)%nblks_c
          IF (jb == p_patch(jg)%nblks_c) THEN
            i_endidx = p_patch(jg)%npromz_c
          ELSE
            i_endidx = nproma
          ENDIF

          DO jt = 2, ntiles_total
            DO jc = 1, i_endidx
              lnd_diag%freshsnow_t(jc,jb,jt) = lnd_diag%freshsnow_t(jc,jb,1)
              lnd_diag%h_snow_t(jc,jb,jt)    = lnd_diag%h_snow_t(jc,jb,1)
              lnd_prog%w_snow_t(jc,jb,jt)    = lnd_prog%w_snow_t(jc,jb,1)
              lnd_prog%rho_snow_t(jc,jb,jt)  = lnd_prog%rho_snow_t(jc,jb,1)
            ENDDO
          ENDDO

        ENDDO
!$OMP END PARALLEL DO
      ENDIF

    ENDDO ! loop over model domains


  END SUBROUTINE read_dwdana_sfc


END MODULE mo_initicon_io

