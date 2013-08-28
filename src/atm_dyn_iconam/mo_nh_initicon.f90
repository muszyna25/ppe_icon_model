!>
!! This module contains the I/O routines for initicon
!!
!! @author Guenther Zaengl, DWD
!!
!!
!! @par Revision History
!! First version by Guenther Zaengl, DWD (2011-07-13)
!!
!! @par Copyright
!! 2002-2010 by DWD and MPI-M
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

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_nh_initicon

  USE mo_kind,                ONLY: wp
  USE mo_io_units,            ONLY: filename_max
  USE mo_parallel_config,     ONLY: nproma, p_test_run
  USE mo_run_config,          ONLY: msg_level, nvclev, iqv, iqc, iqi, iqr, iqs
  USE mo_extpar_config,       ONLY: extpar_filename,     &
    &                               generate_extpar_filename => generate_filename
  USE mo_dynamics_config,     ONLY: nnow, nnow_rcf, nnew, nnew_rcf
  USE mo_model_domain,        ONLY: t_patch
  USE mo_nonhydro_types,      ONLY: t_nh_state, t_nh_prog, t_nh_diag
  USE mo_nwp_phy_types,       ONLY: t_nwp_phy_diag
  USE mo_nonhydrostatic_config, ONLY: kstart_moist
  USE mo_nwp_lnd_types,       ONLY: t_lnd_state
  USE mo_intp_data_strc,      ONLY: t_int_state
  USE mo_ext_data_types,      ONLY: t_external_data
  USE mo_grf_intp_data_strc,  ONLY: t_gridref_state
  USE mo_nh_initicon_types,   ONLY: t_initicon_state
  USE mo_initicon_config,     ONLY: init_mode, nlev_in, nlevsoil_in, l_sst_in,          &
    &                               ifs2icon_filename, dwdfg_filename, dwdana_filename, &
    &                               generate_filename,                                  &
    &                               nml_filetype => filetype,                           &
    &                               ana_varnames_map_file, l_ana_sfc
  USE mo_impl_constants,      ONLY: SUCCESS, MAX_CHAR_LENGTH, max_dom, MODE_DWDANA,     &
    &                               MODE_IFSANA, MODE_COMBINED, min_rlcell,             &
    &                               min_rledge, min_rledge_int, min_rlcell_int
  USE mo_impl_constants_grf,  ONLY: grf_bdywidth_c, grf_bdywidth_e
  USE mo_physical_constants,  ONLY: tf_salt, rd, cpd, cvd, p0ref, vtmpc1, grav
  USE mo_exception,           ONLY: message, finish, message_text
  USE mo_grid_config,         ONLY: n_dom, nroot
  USE mo_mpi,                 ONLY: p_pe, p_io, p_bcast, p_comm_work_test, p_comm_work
  USE mo_netcdf_read,         ONLY: read_netcdf_data, read_netcdf_data_single, nf
  USE mo_util_grib,           ONLY: read_grib_2d, read_grib_3d, get_varID
  USE mo_nh_init_utils,       ONLY: hydro_adjust, convert_thdvars, interp_uv_2_vn, init_w
  USE mo_util_phys,           ONLY: virtual_temp
  USE mo_util_string,         ONLY: tolower, difference, add_to_list, one_of
  USE mo_util_file,           ONLY: util_filesize
  USE mo_ifs_coord,           ONLY: alloc_vct, init_vct, vct, vct_a, vct_b
  USE mo_lnd_nwp_config,      ONLY: nlev_soil, ntiles_total, lmulti_snow, nlev_snow, lseaice,&
    &                               llake, isub_lake, ntiles_water
  USE mo_atm_phy_nwp_config,  ONLY: atm_phy_nwp_config
  USE mo_master_nml,          ONLY: model_base_dir
  USE mo_phyparam_soil,       ONLY: csalb_snow_min, csalb_snow_max,crhosmin_ml,crhosmax_ml
  USE mo_seaice_nwp,          ONLY: frsi_min, seaice_coldinit_nwp
  USE mo_nh_vert_interp,      ONLY: vert_interp_atm, vert_interp_sfc
  USE mo_intp_rbf,            ONLY: rbf_vec_interpol_cell
  USE mo_nh_diagnose_pres_temp,ONLY: diagnose_pres_temp
  USE mo_loopindices,         ONLY: get_indices_c, get_indices_e
  USE mo_sync,                ONLY: sync_patch_array, SYNC_E, SYNC_C
  USE mo_math_laplace,        ONLY: nabla4_vec
  USE mo_dictionary,          ONLY: t_dictionary, dict_init, dict_finalize, &
    &                               dict_loadfile, dict_get, DICT_MAX_STRLEN
  USE mo_post_op,             ONLY: perform_post_op
  USE mo_var_metadata,        ONLY: t_var_metadata, POST_OP_NONE, VARNAME_LEN
  USE mo_linked_list,         ONLY: t_list_element
  USE mo_var_list,            ONLY: get_var_name, nvar_lists, var_lists, collect_group, &
    &                               total_number_of_variables
  USE mo_var_list_element,    ONLY: level_type_ml
  USE mo_cdi_constants,       ONLY: cdiDefAdditionalKey, filetype_nc2, filetype_grb2, &
    &                               zaxisInqType, ZAXIS_HYBRID, ZAXIS_HYBRID_HALF,    &
    &                               ZAXIS_HEIGHT, vlistInqVarZaxis, vlistNvars,       &
    &                               streamInqVlist, streamOpenRead, streamInqNvars
  USE mo_nwp_sfc_interp,      ONLY: smi_to_sm_mass
  USE mo_util_cdi_table,      ONLY: print_cdi_summary
  USE mo_util_bool_table,     ONLY: init_bool_table, add_column, print_bool_table, &
    &                               t_bool_table
  USE mo_nwp_phy_state,       ONLY: prm_nwp_diag_list
  USE mo_data_flake,          ONLY: rflk_depth_bs_ref, tpl_T_f, tpl_T_r, C_T_min

  IMPLICIT NONE

  INCLUDE 'netcdf.inc'

  PRIVATE

  CHARACTER(len=*), PARAMETER :: version = '$Id$'


  TYPE(t_initicon_state), ALLOCATABLE, TARGET :: initicon(:) 

  ! NetCDF file IDs / CDI stream IDs for first guess and analysis file
  INTEGER, ALLOCATABLE :: fileID_fg(:),   fileID_ana(:)

  ! file type (NetCDF/GRB2) for first guess and analysis file
  INTEGER, ALLOCATABLE :: filetype_fg(:), filetype_ana(:)

  ! filename: first guess
  CHARACTER(LEN=filename_max) :: dwdfg_file(max_dom)

  ! filename: analysis
  CHARACTER(LEN=filename_max) :: dwdana_file(max_dom)
    

  CHARACTER(LEN=10) :: psvar 
  CHARACTER(LEN=10) :: geop_ml_var  ! model level surface geopotential
  CHARACTER(LEN=10) :: geop_sfc_var ! surface-level surface geopotential
  CHARACTER(LEN=10) :: alb_snow_var ! snow albedo


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

  ! dictionary which maps internal variable names onto
  ! GRIB2 shortnames or NetCDF var names.
  TYPE (t_dictionary) :: ana_varnames_dict


  PUBLIC :: init_icon



  CONTAINS

  !-------------
  !>
  !! SUBROUTINE init_icon
  !! Initialization routine of init_icon: Reads in either DWD or IFS analysis
  !!
  !! @par Revision History
  !! Initial version by Guenther Zaengl, DWD(2011-07-14)
  !!
  !!
  SUBROUTINE init_icon (p_patch, p_nh_state, prm_diag, p_lnd_state, p_int_state, &
    &                   p_grf_state, ext_data)

    TYPE(t_patch),          INTENT(IN)    :: p_patch(:)
    TYPE(t_nh_state),       INTENT(INOUT) :: p_nh_state(:)
    TYPE(t_nwp_phy_diag),   INTENT(INOUT) :: prm_diag(:)
    TYPE(t_lnd_state),      INTENT(INOUT) :: p_lnd_state(:)
    TYPE(t_int_state),      INTENT(IN)    :: p_int_state(:)
    TYPE(t_gridref_state),  INTENT(IN)    :: p_grf_state(:)

    TYPE(t_external_data),  INTENT(INOUT), OPTIONAL :: ext_data(:)

    INTEGER :: jg, ist

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
      routine = 'mo_nh_initicon:init_icon'



    ! Allocate initicon data type
    ALLOCATE (initicon(n_dom),                         &
      &       filetype_fg(n_dom), filetype_ana(n_dom), &
      &       fileID_fg(n_dom),   fileID_ana(n_dom),   &
      &       stat=ist)
    IF (ist /= SUCCESS)  CALL finish(TRIM(routine),'allocation for initicon failed')

    initicon(:)%atm_in%linitialized = .FALSE.
    initicon(:)%sfc_in%linitialized = .FALSE.

    ! Allocate memory for init_icon state
    CALL allocate_initicon (p_patch, initicon)


    ! Copy the topography fields and coordinate surfaces to initicon
    !
    DO jg = 1, n_dom
      initicon(jg)%topography_c(:,:) = ext_data(jg)%atm%topography_c(:,:)
      initicon(jg)%topography_v(:,:) = ext_data(jg)%atm%topography_v(:,:)

      initicon(jg)%z_ifc(:,:,:) = p_nh_state(jg)%metrics%z_ifc(:,:,:)
      initicon(jg)%z_mc (:,:,:) = p_nh_state(jg)%metrics%z_mc (:,:,:)
    ENDDO



    ! -----------------------------------------------
    ! make the CDI aware of some custom GRIB keys
    ! -----------------------------------------------

    CALL cdiDefAdditionalKey("localInformationNumber")
    CALL cdiDefAdditionalKey("localNumberOfExperiment")
    CALL cdiDefAdditionalKey("typeOfFirstFixedSurface")
    CALL cdiDefAdditionalKey("typeOfGeneratingProcess")
    CALL cdiDefAdditionalKey("backgroundProcess")


    ! -----------------------------------------------
    ! open files containing first guess and analysis
    ! and generate analysis/FG input lists
    ! -----------------------------------------------
    !
    IF ((init_mode == MODE_DWDANA) .OR.   &  ! read in DWD analysis
      & (init_mode == MODE_COMBINED)) THEN
      CALL open_init_files(p_patch, fileID_fg, fileID_ana, filetype_fg, filetype_ana, &
        &                  dwdfg_file, dwdana_file)

      ! Generate lists of fields that must be read from FG/ANA files
      !
      DO jg = 1, n_dom
        CALL create_input_groups(                                            &
          &   initicon(jg)%sfc%grp_vars_fg,  initicon(jg)%sfc%ngrp_vars_fg,  &
          &   initicon(jg)%sfc%grp_vars_ana, initicon(jg)%sfc%ngrp_vars_ana, &
          &   initicon(jg)%sfc%grp_vars_fg_default , initicon(jg)%sfc%ngrp_vars_fg_default,  &
          &   initicon(jg)%sfc%grp_vars_ana_default, initicon(jg)%sfc%ngrp_vars_ana_default, &
          &   init_mode)
      ENDDO

    END IF




    ! init ICON prognostic fields
    !
    SELECT CASE(init_mode)
    CASE(MODE_DWDANA)   ! read in DWD analysis

      CALL message(TRIM(routine),'Real-data mode: perform initialization with DWD analysis')

      ! process DWD atmosphere analysis data
      CALL process_dwdana_atm (p_patch, p_nh_state, p_int_state)

      ! process DWD land/surface analysis data
      CALL process_dwdana_sfc (p_patch, prm_diag, p_lnd_state, ext_data)

    CASE(MODE_IFSANA)   ! read in IFS analysis

      CALL message(TRIM(routine),'Real-data mode: perform initialization with IFS2ICON data')

      ! process IFS atmosphere analysis data
      CALL process_ifsana_atm (p_patch, p_nh_state, p_int_state, p_grf_state, initicon)

      ! process IFS land/surface analysis data
      CALL process_ifsana_sfc (p_patch, p_lnd_state, initicon, ext_data)

    CASE(MODE_COMBINED)

      CALL message(TRIM(routine),'Real-data mode: IFS-atm + GME-soil')

      ! process IFS atmosphere analysis data
      CALL process_ifsana_atm (p_patch, p_nh_state, p_int_state, p_grf_state, initicon)

      ! process DWD land/surface analysis
      CALL process_dwdana_sfc (p_patch, prm_diag, p_lnd_state, ext_data)

    CASE DEFAULT

      CALL finish(TRIM(routine), "Invalid operation mode!")
    END SELECT


    ! Deallocate initicon data type
    !
    CALL deallocate_initicon(initicon)
    CALL deallocate_ifs_atm (initicon)
    CALL deallocate_ifs_sfc (initicon)


    ! close first guess and analysis files
    ! 
    IF ((init_mode == MODE_DWDANA) .OR.   &
      & (init_mode == MODE_COMBINED)) THEN
      CALL close_init_files(fileID_fg, fileID_ana)
    END IF

    DEALLOCATE (initicon, filetype_fg, filetype_ana, fileID_fg, fileID_ana, stat=ist)
    IF (ist /= success) CALL finish(TRIM(routine),'deallocation for initicon failed')

    ! splitting of sea-points list into open water and sea-ice points could be placed 
    ! here, instead of nwp_phy_init/init_nwp_phy
    ! however, one needs to make sure that it is called for both restart and non-restart runs
    ! for both restart and non-restart runs. Could not be included into 
    ! mo_ext_data_state/init_index_lists due to its dependence on p_diag_lnd.
!DR    CALL init_sea_lists(p_patch, ext_data, p_diag_lnd, lseaice)

  END SUBROUTINE init_icon


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
    CHARACTER(*), PARAMETER :: routine = "mo_nh_initicon::open_init_files"
    INTEGER :: jg, vlistID, jlev, mpi_comm
    INTEGER :: flen_fg, flen_ana                     ! filesize in bytes
    LOGICAL :: l_exist

    CALL cdiDefMissval(cdimissval) 
    fileID_fg(:)  = -1
    fileID_ana(:) = -1
    
    IF(p_pe == p_io) THEN
      ! first guess ("_fg")
      DO jg=1,n_dom
        jlev = p_patch(jg)%level
        ! generate file name
        dwdfg_file(jg) = generate_filename(dwdfg_filename, model_base_dir, nroot, jlev, jg)
        INQUIRE (FILE=dwdfg_file(jg), EXIST=l_exist)
        IF (.NOT.l_exist) THEN
          CALL finish(TRIM(routine),'DWD FG file not found: '//TRIM(dwdfg_file(jg)))
        ENDIF
        filetype_fg(jg) = get_filetype(TRIM(dwdfg_file(jg))) ! determine filetype
        SELECT CASE(filetype_fg(jg))
        CASE (FILETYPE_NC2)
          CALL nf(nf_open(TRIM(dwdfg_file(jg)), NF_NOWRITE, fileID_fg(jg)), routine)
        CASE (FILETYPE_GRB2)
          WRITE (0,"(a)") " "
          WRITE (0,"(a,a)") "file inventory: ", TRIM(dwdfg_file(jg))
          fileID_fg(jg)  = streamOpenRead(TRIM(dwdfg_file(jg)))

          ! check whether the file is empty (does not work unfortunately; internal CDI error)
          flen_fg = util_filesize(TRIM(dwdfg_file(jg)))
          IF (flen_fg <= 0 ) THEN
            WRITE(message_text,'(a)') 'File '//TRIM(dwdfg_file(jg))//' is empty'
            CALL message(TRIM(routine), TRIM(message_text))
            CALL finish(routine, "Arrggh!: Empty input file")
          ENDIF

          vlistID = streamInqVlist(fileID_fg(jg))
          CALL print_cdi_summary(vlistID)
        CASE default
          CALL finish(routine, "Internal error!")
        END SELECT
      END DO
      ! analysis ("_ana")
      DO jg=1,n_dom
        jlev = p_patch(jg)%level
        ! generate file name
        dwdana_file(jg) = generate_filename(dwdana_filename, model_base_dir, nroot, jlev, jg)
        INQUIRE (FILE=dwdana_file(jg), EXIST=l_exist)
        IF (.NOT.l_exist) THEN
          CALL finish(TRIM(routine),'DWD ANA file not found: '//TRIM(dwdana_file(jg)))
        ENDIF
        filetype_ana(jg) = get_filetype(TRIM(dwdana_file(jg))) ! determine filetype
        SELECT CASE(filetype_ana(jg))
        CASE (FILETYPE_NC2)
          CALL nf(nf_open(TRIM(dwdana_file(jg)), NF_NOWRITE, fileID_ana(jg)), routine)
        CASE (FILETYPE_GRB2)
          WRITE (0,"(a)") " "
          WRITE (0,"(a,a)") "file inventory: ", TRIM(dwdana_file(jg))
          fileID_ana(jg)  = streamOpenRead(TRIM(dwdana_file(jg)))

          ! check whether the file is empty (does not work unfortunately; internal CDI error)
          flen_ana = util_filesize(TRIM(dwdana_file(jg)))
          IF (flen_ana <= 0 ) THEN
            WRITE(message_text,'(a)') 'File '//TRIM(dwdana_file(jg))//' is empty'
            CALL message(TRIM(routine), TRIM(message_text))
            CALL finish(routine, "Arrggh!: Empty analysis file")
          ENDIF

          vlistID = streamInqVlist(fileID_ana(jg))
          CALL print_cdi_summary(vlistID)
        CASE default
          CALL finish(routine, "Internal error!")
        END SELECT
      END DO
    END IF

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
    CHARACTER(*), PARAMETER :: routine = "mo_nh_initicon::close_init_files"
    INTEGER :: jg

    IF(p_pe == p_io) THEN
      ! first guess ("_fg")
      DO jg=1,n_dom
        IF (fileID_fg(jg) == -1) CYCLE
        SELECT CASE(filetype_fg(jg))
        CASE (FILETYPE_NC2)
          CALL nf(nf_close(fileID_fg(jg)), routine)
        CASE (FILETYPE_GRB2)
          CALL streamClose(fileID_fg(jg))
        CASE default
          CALL finish(routine, "Internal error!")
        END SELECT
      END DO
      ! analysis ("_ana")
      DO jg=1,n_dom
        IF (fileID_ana(jg) == -1) CYCLE
        SELECT CASE(filetype_ana(jg))
        CASE (FILETYPE_NC2)
          CALL nf(nf_close(fileID_ana(jg)), routine)
        CASE (FILETYPE_GRB2)
          CALL streamClose(fileID_ana(jg))
        CASE default
          CALL finish(routine, "Internal error!")
        END SELECT
      END DO
    END IF
    fileID_fg(:)  = -1
    fileID_ana(:) = -1
  END SUBROUTINE close_init_files


  !-------------
  !>
  !! SUBROUTINE process_dwdana_atm
  !! Initialization routine of icon:
  !! - Reads DWD first guess and DA increments (atmosphere only). 
  !!   Data are directly written to the prognostic NH state and added.
  !! - resulting fields are converted to the NH set of prognostic variables
  !!
  !! @par Revision History
  !! Initial version by Daniel Reinert, DWD(2012-12-20)
  !!
  !!
  SUBROUTINE process_dwdana_atm (p_patch, p_nh_state, p_int_state)

    TYPE(t_patch),          INTENT(IN)    :: p_patch(:)
    TYPE(t_nh_state),       INTENT(INOUT) :: p_nh_state(:)
    TYPE(t_int_state),      INTENT(IN)    :: p_int_state(:)

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
      routine = 'mo_nh_initicon:process_dwdana_atm'

!-------------------------------------------------------------------------


    ! read DWD first guess and analysis from DA for atmosphere
    ! 
    CALL read_dwdana_atm(p_patch, p_nh_state)


    ! merge first guess with DA analysis and 
    ! convert variables to the NH set of prognostic variables
    CALL create_dwdana_atm(p_patch, p_nh_state, p_int_state)

  END SUBROUTINE process_dwdana_atm



  !-------------
  !>
  !! SUBROUTINE process_dwdana_sfc
  !! Initialization routine of icon:
  !! - Reads DWD first guess (land/surface only). Data are directly
  !!   written to the prognostic model variables
  !! - reads DWD DA inrements (land/surface only). Data are written 
  !!   to intermediate initicon variables
  !! - first guess and increments are added and resulting fields are 
  !!   converted to the NH set of prognostic variables
  !!
  !! @par Revision History
  !! Initial version by Daniel Reinert, DWD(2012-12-20)
  !!
  !!
  SUBROUTINE process_dwdana_sfc (p_patch, prm_diag, p_lnd_state, ext_data)

    TYPE(t_patch),          INTENT(IN)    :: p_patch(:)
    TYPE(t_nwp_phy_diag),   INTENT(INOUT) :: prm_diag(:)
    TYPE(t_lnd_state),      INTENT(INOUT) :: p_lnd_state(:)
    TYPE(t_external_data),  INTENT(INOUT) :: ext_data(:)


    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
      routine = 'mo_nh_initicon:process_dwdana_sfc'

!-------------------------------------------------------------------------


    ! read DWD first guess and analysis for surface/land
    ! 
    CALL read_dwdana_sfc(p_patch, prm_diag, p_lnd_state)


    ! add increments to first guess and convert variables to 
    ! the NH set of prognostic variables
    CALL create_dwdana_sfc(p_patch, p_lnd_state, ext_data)

  END SUBROUTINE process_dwdana_sfc



  !-------------
  !>
  !! SUBROUTINE process_ifsana_atm
  !! Initialization routine of icon:
  !! - Reads IFS analysis data (atmosphere only)
  !! - performs vertical interpolation from intermediate IFS2ICON grid to ICON 
  !!   grid and converts variables to the NH set of prognostic variables
  !! - finally copies the results to the prognostic model variables
  !!
  !! @par Revision History
  !! Initial version by Daniel Reinert, DWD(2012-12-19)
  !!
  !!
  SUBROUTINE process_ifsana_atm (p_patch, p_nh_state, p_int_state, p_grf_state, &
    &                            initicon)

    TYPE(t_patch),          INTENT(IN)    :: p_patch(:)
    TYPE(t_nh_state),       INTENT(INOUT) :: p_nh_state(:)
    TYPE(t_int_state),      INTENT(IN)    :: p_int_state(:)
    TYPE(t_gridref_state),  INTENT(IN)    :: p_grf_state(:)

    TYPE(t_initicon_state), INTENT(INOUT) :: initicon(:)


    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
      routine = 'mo_nh_initicons:process_ifsana_atm'

!-------------------------------------------------------------------------


    ! read horizontally interpolated IFS analysis for atmosphere
    ! 
    CALL read_ifs_atm(p_patch, initicon)


    ! Perform vertical interpolation from intermediate IFS2ICON grid to ICON grid
    ! and convert variables to the NH set of prognostic variables
    !
    CALL vert_interp_atm(p_patch, p_nh_state, p_int_state, p_grf_state, initicon)

    
    ! Finally copy the results to the prognostic model variables
    !
    CALL copy_initicon2prog_atm(p_patch, initicon, p_nh_state)


  END SUBROUTINE process_ifsana_atm




  !-------------
  !>
  !! SUBROUTINE process_ifsana_sfc
  !! Initialization routine of icon:
  !! - Reads IFS analysis data (surface/land only)
  !! - performs vertical interpolation from intermediate IFS2ICON grid to ICON 
  !!   grid and converts variables to the NH set of prognostic variables
  !! - finally copies the results to the prognostic model variables
  !!
  !! @par Revision History
  !! Initial version by Daniel Reinert, DWD(2012-12-19)
  !!
  !!
  SUBROUTINE process_ifsana_sfc (p_patch, p_lnd_state, initicon, ext_data)

    TYPE(t_patch),          INTENT(IN)    :: p_patch(:)
    TYPE(t_lnd_state),      INTENT(INOUT) :: p_lnd_state(:)

    TYPE(t_initicon_state), INTENT(INOUT) :: initicon(:)
    TYPE(t_external_data),  INTENT(INOUT) :: ext_data(:)


    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
      routine = 'mo_nh_initicon:process_ifsana_sfc'

!-------------------------------------------------------------------------


    ! read horizontally interpolated IFS analysis for surface/land
    ! 
    CALL read_ifs_sfc(p_patch, initicon, ext_data)


    ! Perform vertical interpolation from intermediate IFS2ICON grid to ICON grid
    ! and convert variables to the NH set of prognostic variables
    !
    CALL vert_interp_sfc(p_patch, initicon)

    
    ! Finally copy the results to the prognostic model variables
    !
    CALL copy_initicon2prog_sfc(p_patch, initicon, p_lnd_state, ext_data)


  END SUBROUTINE process_ifsana_sfc


  !-------------------------------------------------------------------------
  !> Wrapper routine for NetCDF and GRIB2 read routines, 2D case.
  !  If necessary an inverse post_op is performed for the input field
  !
  SUBROUTINE read_data_2d (filetype, fileid, varname, glb_arr_len, loc_arr_len, &
    &                      glb_index, var_out, opt_tileidx, opt_checkgroup)

    INTEGER,          intent(IN)    :: filetype       !< FILETYPE_NC2 or FILETYPE_GRB2
    CHARACTER(len=*), INTENT(IN)    :: varname        !< var name of field to be read
    INTEGER,          INTENT(IN)    :: fileid         !< id of netcdf file or stream ID of GRIB2 file
    INTEGER,          INTENT(IN)    :: glb_arr_len    !< length of 1D field (global)
    INTEGER,          INTENT(IN)    :: loc_arr_len    !< length of 1D field (local)
    INTEGER,          INTENT(IN)    :: glb_index(:)   !< Index mapping local to global
    REAL(wp),         INTENT(INOUT) :: var_out(:,:)   !< output field
    INTEGER,          INTENT(IN), OPTIONAL :: opt_tileidx  !< tile index, encoded as "localInformationNumber"
    CHARACTER(LEN=VARNAME_LEN), INTENT(IN), OPTIONAL :: opt_checkgroup(:) !< read only, if varname is 
                                                                          !< contained in opt_checkgroup
    ! local variables
    CHARACTER(len=*), PARAMETER     :: routine = 'mo_nh_initicon:read_data_2d'
    CHARACTER(LEN=DICT_MAX_STRLEN)  :: mapped_name
    LOGICAL                         :: lread          !< .FALSE.: skip reading


    IF (PRESENT(opt_checkgroup)) THEN
      lread = ( one_of(varname,  opt_checkgroup(:)) /= -1)
!!$      write(0,*) "lread: varname",lread, varname 
    ELSE
      lread = .TRUE.
    ENDIF



    IF (lread) THEN
      ! Search name mapping for name in NetCDF/GRIB2 file
      mapped_name = TRIM(dict_get(ana_varnames_dict, varname, default=varname))

      SELECT CASE(filetype)
      CASE (FILETYPE_NC2)
        IF(p_pe == p_io .AND. msg_level>10) THEN
          WRITE(message_text,'(a)') 'NC-Read: '//TRIM(varname)
          CALL message(TRIM(routine), TRIM(message_text))
        ENDIF

        CALL read_netcdf_data (fileid, varname, glb_arr_len, loc_arr_len, &
          &                    glb_index, var_out)

      CASE (FILETYPE_GRB2)
        IF(p_pe == p_io .AND. msg_level>10) THEN
          WRITE(message_text,'(a)') 'GRB2-Read: '//TRIM(mapped_name)
          CALL message(TRIM(routine), TRIM(message_text))
        ENDIF

        CALL read_grib_2d(fileid, mapped_name, glb_arr_len, loc_arr_len, &
          &                    glb_index, var_out, opt_tileidx)

      CASE DEFAULT
        CALL finish(routine, "Unknown file type")
      END SELECT


      ! Perform inverse post_op on input field, if necessary
      !
      CALL initicon_inverse_post_op(varname, mapped_name, optvar_out2D=var_out)

    ENDIF  ! lread

  END SUBROUTINE read_data_2d



  !-------------------------------------------------------------------------
  !> Wrapper routine for NetCDF and GRIB2 read routines, 3D case.
  !
  SUBROUTINE read_data_3d (filetype, fileid, varname, glb_arr_len, loc_arr_len, &
    &                      glb_index, nlevs, var_out, opt_tileidx, opt_checkgroup)

    INTEGER,          INTENT(IN)    :: filetype       !< FILETYPE_NC2 or FILETYPE_GRB2
    CHARACTER(len=*), INTENT(IN)    :: varname        !< var name of field to be read
    INTEGER,          INTENT(IN)    :: fileid         !< id of netcdf file or stream ID of GRIB2 file
    INTEGER,          INTENT(IN)    :: nlevs          !< vertical levels of netcdf file
    INTEGER,          INTENT(IN)    :: glb_arr_len    !< length of 1D field (global)
    INTEGER,          INTENT(IN)    :: loc_arr_len    !< length of 1D field (local)
    INTEGER,          INTENT(IN)    :: glb_index(:)   !< Index mapping local to global
    REAL(wp),         INTENT(INOUT) :: var_out(:,:,:) !< output field
    INTEGER,          INTENT(IN), OPTIONAL :: opt_tileidx  !< tile index, encoded as "localInformationNumber"
    CHARACTER(LEN=VARNAME_LEN), INTENT(IN), OPTIONAL :: opt_checkgroup(:) !< read only, if varname is 
                                                                          !< contained in opt_checkgroup
    ! local variables
    CHARACTER(len=*), PARAMETER     :: routine = 'mo_nh_initicon:read_data_3d'
    CHARACTER(LEN=DICT_MAX_STRLEN)  :: mapped_name
    LOGICAL                         :: lread          !< .FALSE.: skip reading


    IF (PRESENT(opt_checkgroup)) THEN
      lread = ( one_of(varname,  opt_checkgroup(:)) /= -1)
!!$      write(0,*) "lread: varname",lread, varname 
    ELSE
      lread = .TRUE.
    ENDIF


    IF (lread) THEN
      ! Search name mapping for name in NetCDF/GRIB2 file
      mapped_name = TRIM(dict_get(ana_varnames_dict, varname, default=varname))

      SELECT CASE(filetype)
      CASE (FILETYPE_NC2)
        IF(p_pe == p_io .AND. msg_level>10 ) THEN
          WRITE(message_text,'(a)') 'NC-Read: '//TRIM(varname)
          CALL message(TRIM(routine), TRIM(message_text))
        ENDIF

        CALL read_netcdf_data_single (fileid, varname, glb_arr_len, loc_arr_len, &
          &                           glb_index, nlevs, var_out)

      CASE (FILETYPE_GRB2)
        IF(p_pe == p_io .AND. msg_level>10 ) THEN
          WRITE(message_text,'(a)') 'GRB2-Read: '//TRIM(mapped_name)
          CALL message(TRIM(routine), TRIM(message_text))
        ENDIF

        CALL read_grib_3d            (fileid, mapped_name, glb_arr_len, loc_arr_len, &
          &                           glb_index, nlevs, var_out, opt_tileidx)

      CASE DEFAULT
        CALL finish(routine, "Unknown file type")
      END SELECT

      ! Perform inverse post_op on input field, if necessary
      !    
      CALL initicon_inverse_post_op(varname, mapped_name, optvar_out3D=var_out)
    ENDIF  ! lread

  END SUBROUTINE read_data_3d



  !-------------
  !>
  !! SUBROUTINE initicon_inverse_post_op
  !! Perform inverse post_op on input field, if necessary 
  !!
  !! @par Revision History
  !! Initial version by Daniel Reinert, DWD(2013-07-05)
  !!
  !!
  SUBROUTINE initicon_inverse_post_op(varname, mapped_name, optvar_out2D, optvar_out3D)

    CHARACTER(len=*), INTENT(IN)      :: varname             !< var name of field to be read
    CHARACTER(LEN=DICT_MAX_STRLEN)    :: mapped_name         !< mapped input name
    REAL(wp), OPTIONAL, INTENT(INOUT) :: optvar_out2D(:,:)   !< 3D output field
    REAL(wp), OPTIONAL, INTENT(INOUT) :: optvar_out3D(:,:,:) !< 2D output field

    ! local variables
    INTEGER                         :: i              ! loop count
    TYPE(t_var_metadata), POINTER   :: info           ! variable metadata
    TYPE(t_list_element), POINTER   :: element
    CHARACTER(len=*), PARAMETER     :: routine = 'mo_nh_initicon:initicon_inverse_post_op'


    ! Check consistency of optional arguments
    !
    IF (PRESENT( optvar_out2D ) .AND. PRESENT( optvar_out3D )) THEN
      CALL finish(routine, 'Only one optional argument must be present')
    ELSE IF (.NOT.PRESENT( optvar_out2D ) .AND. .NOT.PRESENT( optvar_out3D )) THEN
      CALL finish(routine, 'One of 2 optional arguments must be present')
    ENDIF



    ! get metadata information for field to be read
    info => NULL()
    DO i = 1,nvar_lists
      ! loop only over model level variables
      IF (var_lists(i)%p%vlevel_type /= level_type_ml) CYCLE 

      element => NULL()
      DO
        IF(.NOT.ASSOCIATED(element)) THEN
          element => var_lists(i)%p%first_list_element
        ELSE
          element => element%next_list_element
        ENDIF
        IF(.NOT.ASSOCIATED(element)) EXIT

        ! Check for matching name (take care of suffix of
        ! time-dependent variables):
        IF (TRIM(varname) == TRIM(tolower(get_var_name(element%field)))) THEN
          info => element%field%info
          EXIT
        ENDIF
      END DO

      ! If info handle has been found, exit list loop
      IF (ASSOCIATED(info)) EXIT
    ENDDO

    IF (.NOT.ASSOCIATED(info)) THEN
      WRITE (message_text,'(a,a)') TRIM(varname), 'not found'
      CALL message('',message_text)
      CALL finish(routine, 'Varname does not match any of the ICON variable names')
    ENDIF

    ! perform post_op
    IF (info%post_op%ipost_op_type /= POST_OP_NONE) THEN
      IF(p_pe == p_io .AND. msg_level>10) THEN
        WRITE(message_text,'(a)') 'Inverse Post_op for: '//TRIM(mapped_name)
        CALL message(TRIM(routine), TRIM(message_text))
      ENDIF
      IF (PRESENT(optvar_out2D)) THEN
        CALL perform_post_op(info%post_op, optvar_out2D, opt_inverse=.TRUE.)
      ELSE IF (PRESENT(optvar_out3D)) THEN
        CALL perform_post_op(info%post_op, optvar_out3D, opt_inverse=.TRUE.)
      ENDIF
    ENDIF

  END SUBROUTINE initicon_inverse_post_op



  !-------------
  !>
  !! SUBROUTINE create_input_groups
  !! Generates groups 'in_grp_vars_fg' and 'in_grp_vars_ana', containing the names of the variables which  
  !! must be read from the FG- and ANA-File, respectively.
  !! Both groups are based in the ICON-internal output groups "in_dwd_fg_sfc" and "in_dwd_ana_sfc".
  !! In a first step it is checked, whether the ANA-File contains all members of the group 'in_grp_vars_ana'. 
  !! If a member is missing, it is removed from the group 'in_grp_vars_ana', added to the group 
  !! 'in_grp_vars_fg', and a warning is issued, however the model does not abort. In a second step it is 
  !! checked, whether the FG-File contains all members of the group 'in_grp_vars_fg'. If this is not the 
  !! case, the model aborts.
  !! At the end, a table is printed that shows which variables are part of which input group.
  !!
  !! @par Revision History
  !! Initial version by Daniel Reinert, DWD(2013-07-08)
  !!
  !!
  SUBROUTINE create_input_groups(grp_vars_fg, ngrp_vars_fg, grp_vars_ana, ngrp_vars_ana, &
    &                            grp_vars_fg_default, ngrp_vars_fg_default,              &
    &                            grp_vars_ana_default, ngrp_vars_ana_default,            &
    &                            init_mode)

    CHARACTER(LEN=VARNAME_LEN), INTENT(INOUT) :: grp_vars_fg(:)   ! vars (names) to be read from fg-file
    CHARACTER(LEN=VARNAME_LEN), INTENT(INOUT) :: grp_vars_ana(:)  ! vars (names) to be read from ana-file
    CHARACTER(LEN=VARNAME_LEN), INTENT(INOUT) :: grp_vars_fg_default(:)   ! default vars fg-file
    CHARACTER(LEN=VARNAME_LEN), INTENT(INOUT) :: grp_vars_ana_default(:)  ! default vars ana-file
    INTEGER                   , INTENT(OUT)   :: ngrp_vars_fg     ! number of fields in dwd_fg_sfc_in
    INTEGER                   , INTENT(OUT)   :: ngrp_vars_ana    ! number of fields in dwd_ana_sfc_in
    INTEGER                   , INTENT(OUT)   :: ngrp_vars_fg_default  ! default number dwd_fg_sfc_in
    INTEGER                   , INTENT(OUT)   :: ngrp_vars_ana_default ! default number dwd_ana_sfc_in
    INTEGER,                    INTENT(IN)    :: init_mode        ! initialization mode

    ! local variables
    CHARACTER(LEN=VARNAME_LEN) :: grp_vars_anafile(100)           ! ana-file inventory group
    CHARACTER(LEN=VARNAME_LEN) :: grp_vars_fgfile(100)            ! fg-file inventory group
    INTEGER :: grp_varID_anafile(100)                             ! varIDs of sfc fields in ana-file
    INTEGER :: grp_varID_fgfile (100)                             ! varIDs of sfc fields in fg-file
    INTEGER :: grp_varID_ana    (100)                             ! varIDs of sfc fields in ana-file (re-shuffled)
    INTEGER :: grp_varID_fg     (100)                             ! varIDs of sfc fields in fg-file  (re-shuffled)
    INTEGER :: nfile_totvars                                      ! total number of fields in file
    INTEGER :: ngrp_vars_anafile                                  ! number of sfc fields in ana-file
    INTEGER :: ngrp_vars_fgfile                                   ! number of sfc fields in fg-file
    CHARACTER(LEN=VARNAME_LEN) :: grp_name                        ! group name
    CHARACTER(LEN=MAX_CHAR_LENGTH) :: varname                     ! variable name
    INTEGER :: vlistID, varID, ivar, mpi_comm
    INTEGER :: zaxisID
    INTEGER :: index
!!$    INTEGER :: nlev_ana, nlev_fg                                  ! number of vertical levels (in fg/ana file)

    CHARACTER(LEN=*), PARAMETER :: routine = 'mo_nh_initicon:create_input_groups'
    TYPE(t_bool_table) :: bool_table

    IF(p_pe == p_io) THEN


      ! Collect group 'grp_vars_fg' from dwd_fg_sfc_in
      !
      grp_name ='dwd_fg_sfc_in' 
      CALL collect_group(TRIM(grp_name), grp_vars_fg, ngrp_vars_fg,    &
        &                loutputvars_only=.FALSE.,lremap_lonlat=.FALSE.)

      ! Special HACK for MODE_COMBINED:
      ! remove z0 from grp_vars_fg, if init_mode=MODE_COMBINED, since it must not be read from 
      ! GME soil input data.
      IF (init_mode == MODE_COMBINED) THEN
        CALL difference(grp_vars_fg, ngrp_vars_fg, (/'gz0'/), 1)
      ENDIF


      ! Collect group 'grp_vars_ana' from dwd_ana_sfc_in
      !
      grp_name ='dwd_ana_sfc_in' 
      CALL collect_group(TRIM(grp_name), grp_vars_ana, ngrp_vars_ana,  &
        &                loutputvars_only=.FALSE.,lremap_lonlat=.FALSE.)

      ! Lateron, grp_vars_fg and grp_vars_ana will be the groups that control the reading stuff

      ! Store default groups:
      !
      grp_vars_fg_default (1:ngrp_vars_fg ) = grp_vars_fg (1:ngrp_vars_fg )
      grp_vars_ana_default(1:ngrp_vars_ana) = grp_vars_ana(1:ngrp_vars_ana)
      ngrp_vars_fg_default  = ngrp_vars_fg
      ngrp_vars_ana_default = ngrp_vars_ana                



      !
      ! Generate file inventory lists for FG- and ANA-files
      !

      ! get ANA-file inventory list (surface fields, only)
      !
      ! get vlistID
      vlistID = streamInqVlist(fileID_ana(1))
      ! get total number of fields in analysis file
      nfile_totvars     = vlistNvars(vlistID)
      ! counter
      ngrp_vars_anafile = 0

      DO varID= 0,(nfile_totvars-1)
        ! skip atmospheric fields 
        zaxisID = vlistInqVarZaxis(vlistID,varID)
        IF (ANY((/ZAXIS_HYBRID,ZAXIS_HYBRID_HALF,ZAXIS_HEIGHT/) == zaxisInqType(zaxisID) ) ) CYCLE

        ! increase counter
        ngrp_vars_anafile = ngrp_vars_anafile + 1

        ! get variable name
        CALL vlistInqVarName(vlistID, varID, varname)

        ! In case of GRIB2 translate GRIB2 varname to netcdf varname and add to group
        IF (filetype_ana(1) == FILETYPE_GRB2) THEN
          ! GRIB2 -> Netcdf
          grp_vars_anafile(ngrp_vars_anafile) = TRIM(dict_get(ana_varnames_dict, varname, default=varname, linverse=.TRUE.))
        ELSE IF (filetype_ana(1) == FILETYPE_NC2) THEN
          grp_vars_anafile(ngrp_vars_anafile) = varname
        ENDIF
        ! store corresponding varID
        grp_varID_anafile(ngrp_vars_anafile)= varID
      ENDDO



      ! get FG-file inventory list (surface fields, only)
      !
      ! get vlistID
      vlistID = streamInqVlist(fileID_fg(1))
      ! get total number of fields in FG file
      nfile_totvars    = vlistNvars(vlistID)
      ! counter
      ngrp_vars_fgfile = 0

      DO varID= 0,(nfile_totvars-1)
        ! skip atmospheric fields 
        zaxisID = vlistInqVarZaxis(vlistID,varID)
        IF (ANY((/ZAXIS_HYBRID,ZAXIS_HYBRID_HALF,ZAXIS_HEIGHT/) == zaxisInqType(zaxisID) ) ) CYCLE

        ! increase counter
        ngrp_vars_fgfile = ngrp_vars_fgfile + 1

        ! get variable name
        CALL vlistInqVarName(vlistID, varID, varname)

        ! In case of GRIB2 translate GRIB2 varname to netcdf varname and add to group
        IF (filetype_fg(1) == FILETYPE_GRB2) THEN
          ! GRIB2 -> Netcdf
          grp_vars_fgfile(ngrp_vars_fgfile) = TRIM(dict_get(ana_varnames_dict, varname, default=varname, linverse=.TRUE.))
        ELSE IF (filetype_fg(1) == FILETYPE_NC2) THEN
          grp_vars_fgfile(ngrp_vars_fgfile) = varname
        ENDIF
        ! store corresponding varID
        grp_varID_fgfile(ngrp_vars_fgfile)= varID
      ENDDO

!!$write(0,*) "grp_vars_fg:  ", grp_vars_fg (1:ngrp_vars_fg) ,ngrp_vars_fg
!!$write(0,*) "grp_vars_ana: ", grp_vars_ana(1:ngrp_vars_ana),ngrp_vars_ana
!!$write(0,*) "grp_vars_fgfile: ", grp_vars_fgfile(1:ngrp_vars_fgfile),ngrp_vars_fgfile
!!$write(0,*) "grp_vars_anafile: ", grp_vars_anafile(1:ngrp_vars_anafile),ngrp_vars_anafile



      ! Check, whether the ANA-file inventory list contains all required analysis fields.
      ! If not, remove the corresponding file name from the group 'grp_vars_ana' and issue a warning.
      ! The missing field is added to the group 'grp_vars_fg' and thus the model tries to read 
      ! it from the FG-File.
      !
      DO ivar=1,ngrp_vars_ana_default
        index = one_of(TRIM(grp_vars_ana_default(ivar)),grp_vars_anafile(:))

        IF ( index == -1) THEN  ! variable not found
          WRITE(message_text,'(a)') 'Field '//TRIM(grp_vars_ana_default(ivar))//' not found in ANA-input file.'
          CALL message(routine, TRIM(message_text))
          WRITE(message_text,'(a)') 'Field '//TRIM(grp_vars_ana_default(ivar))//' will be read from FG-input, instead.'
          CALL message(routine, TRIM(message_text))

          ! remove missing field from analysis input-group grp_vars_ana
          CALL difference(grp_vars_ana, ngrp_vars_ana, grp_vars_ana_default(ivar:ivar), 1)

          ! add missing field to the FG-group grp_vars_fg
          CALL add_to_list(grp_vars_fg, ngrp_vars_fg, grp_vars_ana_default(ivar:ivar), 1)
        ELSE
          ! save varID
          grp_varID_ana(ivar) = grp_varID_anafile(index)
        ENDIF
      ENDDO  




!!$      ! Remove those names from grp_vars_fg which are also included in grp_vars_ana
!!$      ! Thus, only those variables are read from the first guess, which are not contained 
!!$      ! in the analysis file. By this, we avoid double reading.
!!$      !
!!$      DO ivar=1,ngrp_vars_ana
!!$
!!$        index = one_of(TRIM(grp_vars_ana(ivar)),grp_vars_fg(:))
!!$
!!$        ! If grp_vars_ana(ivar) is also a group member of grp_vars_fg(:)
!!$        IF (index /= -1) THEN ! found matching name.
!!$          !
!!$          ! Check, whether the number of vertical levels matches as well
!!$
!!$          ! get number of vertical levels for variable grp_vars_ana(ivar)
!!$          varID    = grp_varID_ana(ivar)             ! CDI variable ID with this name in analysis 
!!$          vlistID  = streamInqVlist(fileID_ana(1))
!!$          zaxisID  = vlistInqVarZaxis(vlistID, varID)
!!$          nlev_ana = zaxisInqSize(zaxisID)
!!$
!!$          ! get number of vertical levels for variable grp_vars_fg(index)
!!$          varID    = grp_varID_fg(index)             ! CDI variable ID with this name in analysis 
!!$          vlistID  = streamInqVlist(fileID_fg(1))
!!$          zaxisID  = vlistInqVarZaxis(vlistID, varID)
!!$          nlev_fg  = zaxisInqSize(zaxisID)
!!$
!!$          IF (nlev_ana == nlev_fg) THEN
!!$            ! remove grp_vars_ana(ivar) from grp_vars_fg(:)
!!$            CALL difference(grp_vars_fg, ngrp_vars_fg, grp_vars_ana(ivar:ivar), 1)
!!$          ENDIF
!!$        ENDIF          
!!$      ENDDO



      ! Check, whether the FG-file inventory group contains all fields which are needed for a 
      ! successfull model start. If not, then stop the model and issue an error. 
      DO ivar=1,ngrp_vars_fg
        index = one_of(TRIM(grp_vars_fg(ivar)),grp_vars_fgfile(:))

        IF ( index == -1) THEN   ! variable not found
          WRITE(message_text,'(a)') 'Field '//TRIM(grp_vars_fg(ivar))//' missing in FG-input file.'
          CALL finish(routine, TRIM(message_text))
        ELSE
          ! store varID
          grp_varID_fg(ivar) = grp_varID_fgfile(index)
        ENDIF
      ENDDO  


      ! Printout table
      !
      CALL message("Required surface input fields:",'Source of FG and ANA fields:')
      CALL init_bool_table(bool_table)
      CALL add_column(bool_table, "FG-SFC (default)", grp_vars_fg_default,  ngrp_vars_fg_default)
      CALL add_column(bool_table, "FG-SFC",           grp_vars_fg,          ngrp_vars_fg)
      CALL add_column(bool_table, "ANA-SFC (default)",grp_vars_ana_default, ngrp_vars_ana_default)
      CALL add_column(bool_table, "ANA-SFC",          grp_vars_ana,         ngrp_vars_ana)
      CALL print_bool_table(bool_table)
    ENDIF  ! p_pe == p_io


    IF(p_test_run) THEN
      mpi_comm = p_comm_work_test 
    ELSE
      mpi_comm = p_comm_work
    ENDIF
    CALL p_bcast(grp_vars_fg,  p_io, mpi_comm)
    CALL p_bcast(ngrp_vars_fg, p_io, mpi_comm)
    CALL p_bcast(grp_vars_ana, p_io, mpi_comm)
    CALL p_bcast(ngrp_vars_ana,p_io, mpi_comm)

  END SUBROUTINE create_input_groups



  !-------------------------------------------------------------------------
  !> @return One of CDI's FILETYPE\_XXX constants. Possible values: 2
  !          (=FILETYPE\_GRB2), 4 (=FILETYPE\_NC2)
  !
  !  The file type is determined by the setting of the "filetype"
  !  namelist parameter in "initicon_nml". If this parameter has not
  !  been set, we try to determine the file type by its extension
  !  "*.grb*" or ".nc".
  !
  FUNCTION get_filetype(filename)
    INTEGER :: get_filetype
    CHARACTER(LEN=*), INTENT(IN) :: filename
    ! local variables
    CHARACTER(len=*), PARAMETER :: routine = 'mo_nh_initicon:get_filetype'
    INTEGER :: idx
    
    get_filetype = nml_filetype
    IF (nml_filetype == -1) THEN
      idx = INDEX(tolower(filename),'.nc')
      IF (idx==0) THEN
        idx = INDEX(tolower(filename),'.grb')
        IF (idx==0) THEN
          CALL finish(routine, "File type could not be determined!")
        ELSE
          get_filetype = FILETYPE_GRB2
        END IF
      ELSE
        get_filetype = FILETYPE_NC2
      END IF
    END IF
  END FUNCTION get_filetype


  !>
  !! Read horizontally interpolated IFS analysis (atmosphere only)
  !!
  !! Reads horizontally interpolated IFS analysis atmosphere data
  !! and reads in vertical coordinate table. 
  !!
  !! @par Revision History
  !! Initial version by Guenther Zaengl, DWD(2011-07-14)
  !! Modification by Daniel Reinert, DWD (2012-12-18)
  !! - encapsulate reading of IFS analysis
  !!
  SUBROUTINE read_ifs_atm (p_patch, initicon)

    TYPE(t_patch),          INTENT(IN)    :: p_patch(:)
    TYPE(t_initicon_state), INTENT(INOUT) :: initicon(:)

    INTEGER :: jg, jlev, jk
    LOGICAL :: l_exist

    INTEGER :: no_cells, no_levels
    INTEGER :: ncid, dimid, varid, mpi_comm
    INTEGER :: ist

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
      routine = 'mo_nh_initicon:read_ifs_atm'

    CHARACTER(LEN=filename_max) :: ifs2icon_file(max_dom)

    LOGICAL :: lread_qr, lread_qs ! are qr, qs provided as input?
    LOGICAL :: lread_vn           ! is vn provided as input?

    !-------------------------------------------------------------------------


    DO jg = 1, n_dom

      jlev = p_patch(jg)%level


      ! Skip reading the atmospheric input data if a model domain 
      ! is not active at initial time
      IF (.NOT. p_patch(jg)%ldom_active) CYCLE


      ! Read in data from IFS2ICON
      !
      IF(p_pe == p_io ) THEN 

        !
        ! generate file name
        !
        ifs2icon_file(jg) = generate_filename(ifs2icon_filename, model_base_dir, &
          &                                   nroot, jlev, jg)
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
        ! Check if surface pressure (VN) is provided as input
        !
        IF (nf_inq_varid(ncid, 'VN', varid) == nf_noerr) THEN
          lread_vn = .TRUE.
        ELSE 
          lread_vn = .FALSE.
        ENDIF

      ENDIF ! pe_io

      IF(p_test_run) THEN
        mpi_comm = p_comm_work_test 
      ELSE
        mpi_comm = p_comm_work
      ENDIF

      CALL p_bcast(nlev_in,   p_io, mpi_comm)
      CALL p_bcast(lread_qs,  p_io, mpi_comm)
      CALL p_bcast(lread_qr,  p_io, mpi_comm)
      CALL p_bcast(lread_vn,  p_io, mpi_comm)


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
      CALL allocate_ifs_atm(jg, p_patch(jg)%nblks_c, initicon)

      ! start reading atmospheric fields
      !
      CALL read_netcdf_data_single (ncid, 'T', p_patch(jg)%n_patch_cells_g,           &
        &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
        &                     nlev_in,initicon(jg)%atm_in%temp)


      IF (lread_vn) THEN
        ALLOCATE(initicon(jg)%atm_in%vn(nproma,nlev_in,p_patch(jg)%nblks_e), STAT=ist)
        IF (ist /= SUCCESS) THEN
          CALL finish ( TRIM(routine), 'allocation of atm_in%vn failed')
        ENDIF
        CALL read_netcdf_data_single (ncid, 'VN', p_patch(jg)%n_patch_edges_g,          &
          &                     p_patch(jg)%n_patch_edges, p_patch(jg)%edges%glb_index, &
          &                     nlev_in,initicon(jg)%atm_in%vn)
      ELSE
        CALL read_netcdf_data_single (ncid, 'U', p_patch(jg)%n_patch_cells_g,           &
          &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
          &                     nlev_in,initicon(jg)%atm_in%u)

        CALL read_netcdf_data_single (ncid, 'V', p_patch(jg)%n_patch_cells_g,           &
          &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
          &                     nlev_in,initicon(jg)%atm_in%v)
      ENDIF


      ! Note: input vertical velocity is in fact omega (Pa/s)
      CALL read_netcdf_data_single (ncid, 'W', p_patch(jg)%n_patch_cells_g,           &
        &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
        &                     nlev_in,initicon(jg)%atm_in%omega)

      CALL read_netcdf_data_single (ncid, 'QV', p_patch(jg)%n_patch_cells_g,          &
        &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
        &                     nlev_in,initicon(jg)%atm_in%qv)

      CALL read_netcdf_data_single (ncid, 'QC', p_patch(jg)%n_patch_cells_g,          &
        &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
        &                     nlev_in,initicon(jg)%atm_in%qc)

      CALL read_netcdf_data_single (ncid, 'QI', p_patch(jg)%n_patch_cells_g,          &
        &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
        &                     nlev_in,initicon(jg)%atm_in%qi)

      IF (lread_qr) THEN
        CALL read_netcdf_data_single (ncid, 'QR', p_patch(jg)%n_patch_cells_g,          &
        &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
        &                     nlev_in,initicon(jg)%atm_in%qr)
      ELSE
        initicon(jg)%atm_in%qr(:,:,:)=0._wp
      ENDIF

      IF (lread_qs) THEN
        CALL read_netcdf_data_single (ncid, 'QS', p_patch(jg)%n_patch_cells_g,          &
        &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
        &                     nlev_in,initicon(jg)%atm_in%qs)
      ELSE
        initicon(jg)%atm_in%qs(:,:,:)=0._wp
      ENDIF

      CALL read_netcdf_data (ncid, TRIM(psvar), p_patch(jg)%n_patch_cells_g,          &
        &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
        &                     initicon(jg)%atm_in%psfc)

      CALL read_netcdf_data (ncid, TRIM(geop_ml_var), p_patch(jg)%n_patch_cells_g,    &
        &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
        &                     initicon(jg)%atm_in%phi_sfc)


      ! Allocate and read in vertical coordinate tables
      !
      IF (jg == 1) THEN

        ALLOCATE(vct_a(nlev_in+1), vct_b(nlev_in+1), vct(2*(nlev_in+1)))

        IF(p_pe == p_io) THEN
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

        nvclev = 2*(nlev_in+1)

        CALL alloc_vct(nlev_in)
        CALL init_vct(nlev_in)

        IF (msg_level >= 15) THEN
          WRITE(message_text,'(a)') 'vct table values of input data:'
          CALL message(TRIM(routine), TRIM(message_text))

          DO jk = 1, nlev_in+1
            WRITE(message_text,'(a,i4,F12.4,F12.8)') 'jk, vct_a, vct_b: ',jk, vct_a(jk), vct_b(jk)
            CALL message(TRIM(routine), TRIM(message_text))
          ENDDO
        ENDIF

      ENDIF  ! jg=1

      ! close file
      !
      IF(p_pe == p_io) CALL nf(nf_close(ncid), routine)

    ENDDO ! loop over model domains

  END SUBROUTINE read_ifs_atm





  !>
  !! Read horizontally interpolated IFS analysis (surface only)
  !!
  !! Reads horizontally interpolated IFS analysis surface data
  !!
  !! @par Revision History
  !! Initial version by Guenther Zaengl, DWD(2011-07-14)
  !! Modification by Daniel Reinert, DWD (2012-12-18)
  !! - encapsulate reading of IFS analysis
  !!
  SUBROUTINE read_ifs_sfc (p_patch, initicon, ext_data)

    TYPE(t_patch),          INTENT(IN)    :: p_patch(:)
    TYPE(t_initicon_state), INTENT(INOUT) :: initicon(:)
    TYPE(t_external_data),  INTENT(IN)    :: ext_data(:)

    INTEGER :: jg, jlev
    LOGICAL :: l_exist

    INTEGER :: no_cells, no_levels
    INTEGER :: ncid, dimid, varid, mpi_comm

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
      routine = 'mo_nh_initicon:read_ifs_sfc'

    CHARACTER(LEN=filename_max) :: ifs2icon_file(max_dom)


    !-------------------------------------------------------------------------


    DO jg = 1, n_dom

      jlev = p_patch(jg)%level


      ! Skip reading the atmospheric input data if a model domain 
      ! is not active at initial time
      IF (.NOT. p_patch(jg)%ldom_active) CYCLE


      ! Read in data from IFS2ICON
      !
      IF(p_pe == p_io ) THEN 
        !
        ! generate file name
        !
        ifs2icon_file(jg) = generate_filename(ifs2icon_filename, model_base_dir, &
          &                                   nroot, jlev, jg)
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



      IF(p_test_run) THEN
        mpi_comm = p_comm_work_test 
      ELSE
        mpi_comm = p_comm_work
      ENDIF

      CALL p_bcast(l_sst_in, p_io, mpi_comm)

      CALL p_bcast(alb_snow_var, p_io, mpi_comm)

      ! allocate data structure
      CALL allocate_ifs_sfc(jg, p_patch(jg)%nblks_c, initicon)


      ! start reading surface fields
      !
      CALL read_netcdf_data (ncid, TRIM(geop_sfc_var), p_patch(jg)%n_patch_cells_g,   &
        &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
        &                     initicon(jg)%sfc_in%phi)

      CALL read_netcdf_data (ncid, 'SKT', p_patch(jg)%n_patch_cells_g,                &
        &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
        &                     initicon(jg)%sfc_in%tskin)
      IF ( l_sst_in) THEN
       CALL read_netcdf_data (ncid, 'SST', p_patch(jg)%n_patch_cells_g,                &
         &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
         &                     initicon(jg)%sfc_in%sst)
      ELSE 
       initicon(jg)%sfc_in%sst(:,:)=0.0_wp
      END IF

      CALL read_netcdf_data (ncid, 'T_SNOW', p_patch(jg)%n_patch_cells_g,             &
        &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
        &                     initicon(jg)%sfc_in%tsnow)

      CALL read_netcdf_data (ncid,TRIM(alb_snow_var), p_patch(jg)%n_patch_cells_g,    &
        &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
        &                     initicon(jg)%sfc_in%snowalb)
 
      CALL read_netcdf_data (ncid, 'W_SNOW', p_patch(jg)%n_patch_cells_g,             &
        &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
        &                     initicon(jg)%sfc_in%snowweq)

      CALL read_netcdf_data (ncid,'RHO_SNOW', p_patch(jg)%n_patch_cells_g,            &
        &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
        &                     initicon(jg)%sfc_in%snowdens)

      CALL read_netcdf_data (ncid, 'W_I', p_patch(jg)%n_patch_cells_g,                &
        &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
        &                     initicon(jg)%sfc_in%skinres)

      CALL read_netcdf_data (ncid, 'LSM', p_patch(jg)%n_patch_cells_g,                &
        &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
        &                     initicon(jg)%sfc_in%ls_mask)

      CALL read_netcdf_data (ncid, 'CI', p_patch(jg)%n_patch_cells_g,                 &
        &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
        &                     initicon(jg)%sfc_in%seaice)

      CALL read_netcdf_data (ncid, 'STL1', p_patch(jg)%n_patch_cells_g,               &
        &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
        &                     initicon(jg)%sfc_in%tsoil(:,:,1))

      CALL read_netcdf_data (ncid, 'STL2', p_patch(jg)%n_patch_cells_g,               &
        &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
        &                     initicon(jg)%sfc_in%tsoil(:,:,2))

      CALL read_netcdf_data (ncid, 'STL3', p_patch(jg)%n_patch_cells_g,               &
        &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
        &                     initicon(jg)%sfc_in%tsoil(:,:,3))

      CALL read_netcdf_data (ncid, 'STL4', p_patch(jg)%n_patch_cells_g,               &
        &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
        &                     initicon(jg)%sfc_in%tsoil(:,:,4))

      CALL read_netcdf_data (ncid, 'SMIL1', p_patch(jg)%n_patch_cells_g,              &
        &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
        &                     initicon(jg)%sfc_in%wsoil(:,:,1))

      CALL read_netcdf_data (ncid, 'SMIL2', p_patch(jg)%n_patch_cells_g,              &
        &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
        &                     initicon(jg)%sfc_in%wsoil(:,:,2))

      CALL read_netcdf_data (ncid, 'SMIL3', p_patch(jg)%n_patch_cells_g,              &
        &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
        &                     initicon(jg)%sfc_in%wsoil(:,:,3))

      CALL read_netcdf_data (ncid, 'SMIL4', p_patch(jg)%n_patch_cells_g,              &
        &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
        &                     initicon(jg)%sfc_in%wsoil(:,:,4))



      ! close file
      !
      IF(p_pe == p_io) CALL nf(nf_close(ncid), routine)



      ! In addition, copy climatological deep-soil temperature to soil level nlev_soil
      ! These are limited to -60 deg C because less is definitely nonsense
      initicon(jg)%sfc%tsoil(:,:,nlev_soil) = MAX(213.15_wp,ext_data(jg)%atm%t_cl(:,:))

    ENDDO ! loop over model domains

  END SUBROUTINE read_ifs_sfc



  !>
  !! Read DWD first guess and analysis from DA (atmosphere only)
  !!
  !! Read DWD first guess and analysis from DA (atmosphere only)
  !! First guess (FG) is read for theta_v, rho, vn, w, tke,
  !! whereas DA output is read for T, p, u, v, 
  !! qv, qc, qi, qr, qs.
  !!
  !! @par Revision History
  !! Initial version by Daniel Reinert, DWD(2012-12-18)
  !! Modifications for GRIB2 : F. Prill, DWD (2013-02-19)
  !!
  SUBROUTINE read_dwdana_atm (p_patch, p_nh_state)

    TYPE(t_patch),          INTENT(IN)    :: p_patch(:)
    TYPE(t_nh_state),       INTENT(INOUT) :: p_nh_state(:)

    INTEGER :: jg
    INTEGER :: nlev, nlevp1

    INTEGER :: no_cells, no_cells_2, no_levels, no_levels_2
    INTEGER :: dimid

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
      routine = 'mo_nh_initicon:read_dwdana_atm'

    !-------------------------------------------------------------------------

    DO jg = 1, n_dom

      ! number of vertical full and half levels
      nlev   = p_patch(jg)%nlev
      nlevp1 = p_patch(jg)%nlevp1

      ! Skip reading the atmospheric input data if a model domain 
      ! is not active at initial time
      IF (.NOT. p_patch(jg)%ldom_active) CYCLE

      !---------------------------------------!
      ! Read in DWD first guess (atmosphere)  !
      !---------------------------------------!

      IF(p_pe == p_io ) THEN 
        CALL message (TRIM(routine), 'read atm_FG fields from '//TRIM(dwdfg_file(jg)))

        SELECT CASE(filetype_fg(jg))
        CASE (FILETYPE_NC2)
          ! read basic dimensions

          ! get number of cells
          CALL nf(nf_inq_dimid(fileID_fg(jg), 'ncells', dimid), routine)
          CALL nf(nf_inq_dimlen(fileID_fg(jg), dimid, no_cells), routine)
          CALL nf(nf_inq_dimid(fileID_fg(jg), 'ncells_2', dimid), routine)
          CALL nf(nf_inq_dimlen(fileID_fg(jg), dimid, no_cells_2), routine)
          
          ! get number of vertical levels
          CALL nf(nf_inq_dimid(fileID_fg(jg), 'lev', dimid), routine)
          CALL nf(nf_inq_dimlen(fileID_fg(jg), dimid, no_levels), routine)
          CALL nf(nf_inq_dimid(fileID_fg(jg), 'lev_2', dimid), routine)
          CALL nf(nf_inq_dimlen(fileID_fg(jg), dimid, no_levels_2), routine)

          ! -------------------------
          ! simple consistency checks
          ! -------------------------
          
          ! check the number of cells and vertical levels
          IF( p_patch(jg)%n_patch_cells_g /= no_cells .AND. &
               & p_patch(jg)%n_patch_cells_g /= no_cells_2) THEN
            CALL finish(TRIM(ROUTINE),&
                 & 'Number of patch cells and cells in DWD FG file do not match.')
          ENDIF
          
          IF(p_patch(jg)%nlev /= no_levels .AND. p_patch(jg)%nlev /= no_levels_2) THEN
            CALL finish(TRIM(ROUTINE),&
                 & 'nlev does not match the number of levels in DWD FG file.')
          ENDIF
        END SELECT
      ENDIF  ! p_io

      ! start reading first guess (atmosphere only)
      !
      CALL read_data_3d (filetype_fg(jg), fileID_fg(jg), 'theta_v', p_patch(jg)%n_patch_cells_g,     &
        &                  p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index,     &
        &                  nlev, p_nh_state(jg)%prog(nnow(jg))%theta_v)

      CALL read_data_3d (filetype_fg(jg), fileID_fg(jg), 'rho', p_patch(jg)%n_patch_cells_g,         &
        &                  p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index,     &
        &                  nlev, p_nh_state(jg)%prog(nnow(jg))%rho)

      CALL read_data_3d (filetype_fg(jg), fileID_fg(jg), 'vn', p_patch(jg)%n_patch_edges_g,          &
        &                  p_patch(jg)%n_patch_edges, p_patch(jg)%edges%glb_index,     &
        &                  nlev, p_nh_state(jg)%prog(nnow(jg))%vn)

      CALL read_data_3d (filetype_fg(jg), fileID_fg(jg), 'w', p_patch(jg)%n_patch_cells_g,           &
        &                  p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index,     &
        &                  nlevp1, p_nh_state(jg)%prog(nnow(jg))%w)

      CALL read_data_3d (filetype_fg(jg), fileID_fg(jg), 'tke', p_patch(jg)%n_patch_cells_g,         &
        &                  p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index,     &
        &                  nlevp1, p_nh_state(jg)%prog(nnow(jg))%tke)

    ENDDO ! loop over model domains


    !----------------------------------------!
    ! read in DWD analysis (atmosphere)      !
    !----------------------------------------!

    ! Analysis is read for the global domain, only.
    jg = 1

    IF(p_pe == p_io ) THEN 
      CALL message (TRIM(routine), 'read atm_ANA fields from '//TRIM(dwdana_file(jg)))

      SELECT CASE(filetype_ana(jg))
      CASE (FILETYPE_NC2)
        ! open file, read basic dimensions

        ! get number of cells
        CALL nf(nf_inq_dimid(fileID_ana(jg), 'ncells', dimid), routine)
        CALL nf(nf_inq_dimlen(fileID_ana(jg), dimid, no_cells), routine)
        ! get number of vertical levels
        CALL nf(nf_inq_dimid(fileID_ana(jg), 'lev', dimid), routine)
        CALL nf(nf_inq_dimlen(fileID_ana(jg), dimid, no_levels), routine)

        ! -------------------------
        ! simple consistency checks
        ! -------------------------
        
        ! check the number of cells and vertical levels
        IF( p_patch(jg)%n_patch_cells_g /= no_cells) THEN
          CALL finish(TRIM(ROUTINE),&
               & 'Number of patch cells and cells in DWD ANA file do not match.')
        ENDIF
        
        IF(p_patch(jg)%nlev /= no_levels) THEN
          CALL finish(TRIM(ROUTINE),&
               & 'nlev does not match the number of levels in DWD ANA file.')
        ENDIF
      END SELECT

    ENDIF  ! p_io

    ! start reading DA output (atmosphere only)
    ! The dynamical variables temp, pres, u and v, which need further processing,
    ! are stored in initicon(jg)%atm. The moisture variables, which can be taken
    ! over directly from the Analysis, are written to the NH prognostic state
    !
    CALL read_data_3d (filetype_ana(jg), fileID_ana(jg), 'temp', p_patch(jg)%n_patch_cells_g,      &
      &                  p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index,   &
      &                  nlev, initicon(jg)%atm%temp )

    CALL read_data_3d (filetype_ana(jg), fileID_ana(jg), 'pres', p_patch(jg)%n_patch_cells_g,      &
      &                  p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index,   &
      &                  nlev, initicon(jg)%atm%pres )
    
    CALL read_data_3d (filetype_ana(jg), fileID_ana(jg), 'u', p_patch(jg)%n_patch_cells_g,         &
      &                  p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index,   &
      &                  nlev, initicon(jg)%atm%u )

    CALL read_data_3d (filetype_ana(jg), fileID_ana(jg), 'v', p_patch(jg)%n_patch_cells_g,         &
      &                  p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index,   &
      &                  nlev, initicon(jg)%atm%v )

    CALL read_data_3d (filetype_ana(jg), fileID_ana(jg), 'qv', p_patch(jg)%n_patch_cells_g,        &
      &                  p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index,   &
      &                  nlev, p_nh_state(jg)%prog(nnow(jg))%tracer(:,:,:,iqv))

    ! For the time being identical to qc from FG
    CALL read_data_3d (filetype_ana(jg), fileID_ana(jg), 'qc', p_patch(jg)%n_patch_cells_g,        &
      &                  p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index,   &
      &                  nlev, p_nh_state(jg)%prog(nnow(jg))%tracer(:,:,:,iqc))

    ! For the time being identical to qi from FG
    CALL read_data_3d (filetype_ana(jg), fileID_ana(jg), 'qi', p_patch(jg)%n_patch_cells_g,        &
      &                  p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index,   &
      &                  nlev, p_nh_state(jg)%prog(nnow(jg))%tracer(:,:,:,iqi))

    ! For the time being identical to qr from FG
    CALL read_data_3d (filetype_ana(jg), fileID_ana(jg), 'qr', p_patch(jg)%n_patch_cells_g,        &
      &                  p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index,   &
      &                  nlev, p_nh_state(jg)%prog(nnow(jg))%tracer(:,:,:,iqr))

    ! For the time being identical to qs from FG
    CALL read_data_3d (filetype_ana(jg), fileID_ana(jg), 'qs', p_patch(jg)%n_patch_cells_g,        &
      &                  p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index,   &
      &                  nlev, p_nh_state(jg)%prog(nnow(jg))%tracer(:,:,:,iqs))

  END SUBROUTINE read_dwdana_atm




  !>
  !! Read DWD first guess and analysis (land/surface only)
  !!
  !! Read DWD first guess and analysis for land/surface.
  !! First guess is read for:
  !! fr_seaice, t_ice, h_ice, t_g, qv_s, freshsnow, w_snow, w_i, t_snow, 
  !! rho_snow, w_so, w_so_ice, t_so, gz0
  !!
  !! If available, 
  !!
  !! @par Revision History
  !! Initial version by Daniel Reinert, DWD(2012-12-18)
  !!
  SUBROUTINE read_dwdana_sfc (p_patch, prm_diag, p_lnd_state)

    TYPE(t_patch),          INTENT(IN)    :: p_patch(:)
    TYPE(t_nwp_phy_diag),   INTENT(INOUT) :: prm_diag(:)
    TYPE(t_lnd_state),      INTENT(INOUT) :: p_lnd_state(:)

    INTEGER :: jg, jt

    INTEGER :: no_cells, no_levels
    INTEGER :: no_depth, no_depth_2
    INTEGER :: dimid, mpi_comm
    INTEGER :: ngrp_vars_fg, ngrp_vars_ana

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
      routine = 'mo_nh_initicon:read_dwdana_sfc'

    !-------------------------------------------------------------------------


    DO jg = 1, n_dom

      ! Skip reading the atmospheric input data if a model domain 
      ! is not active at initial time
      IF (.NOT. p_patch(jg)%ldom_active) CYCLE

      ! save some paperwork
      ngrp_vars_fg  = initicon(jg)%sfc%ngrp_vars_fg
      ngrp_vars_ana = initicon(jg)%sfc%ngrp_vars_ana

      !-----------------------------------!
      ! Read in DWD first guess (surface) !
      !-----------------------------------!
      IF(p_pe == p_io ) THEN 
        CALL message (TRIM(routine), 'read sfc_FG fields from '//TRIM(dwdfg_file(jg)))

        SELECT CASE(filetype_fg(jg))
        CASE (FILETYPE_NC2)
          ! open file, read basic dimensions

          ! get number of cells
          CALL nf(nf_inq_dimid(fileID_fg(jg), 'ncells', dimid), routine)
          CALL nf(nf_inq_dimlen(fileID_fg(jg), dimid, no_cells), routine)

          CALL nf(nf_inq_dimid(fileID_fg(jg), 'depth', dimid), routine)
          CALL nf(nf_inq_dimlen(fileID_fg(jg), dimid, no_depth), routine)
          CALL nf(nf_inq_dimid(fileID_fg(jg), 'depth_2', dimid), routine)
          CALL nf(nf_inq_dimlen(fileID_fg(jg), dimid, no_depth_2), routine)
          ! -------------------------
          ! simple consistency checks
          ! -------------------------
          
          IF( p_patch(jg)%n_patch_cells_g /= no_cells ) THEN
            CALL finish(TRIM(ROUTINE),&
                 & 'Number of patch cells and cells in DWD FG file do not match.')
          ENDIF

          IF(nlev_soil /= no_depth_2 ) THEN
            CALL finish(TRIM(ROUTINE),&
                 & 'nlev_soil does not match the number of soil levels in DWD FG file.')
          ENDIF

        CASE (FILETYPE_GRB2)

          !
        CASE DEFAULT
          CALL finish(routine, "Unknown file type")
        END SELECT

      ENDIF  ! p_io

      IF(p_test_run) THEN
        mpi_comm = p_comm_work_test 
      ELSE
        mpi_comm = p_comm_work
      ENDIF


      !----------------------------------------!
      ! read in DWD First Guess (surface)      !
      !----------------------------------------!

      CALL read_data_2d (filetype_fg(jg), fileID_fg(jg), 'fr_seaice', p_patch(jg)%n_patch_cells_g, &
        &                p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index,                   &
        &                p_lnd_state(jg)%diag_lnd%fr_seaice,                                       &
        &                opt_checkgroup=initicon(jg)%sfc%grp_vars_fg(1:ngrp_vars_fg))


      ! sea-ice related fields
      CALL read_data_2d (filetype_fg(jg), fileID_fg(jg), 't_ice', p_patch(jg)%n_patch_cells_g,     &
        &                p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index,                   &
        &                p_lnd_state(jg)%prog_wtr(nnow_rcf(jg))%t_ice,                             &
        &                opt_checkgroup=initicon(jg)%sfc%grp_vars_fg(1:ngrp_vars_fg))

      CALL read_data_2d (filetype_fg(jg), fileID_fg(jg), 'h_ice', p_patch(jg)%n_patch_cells_g,   &
        &                p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index,                 &
        &                p_lnd_state(jg)%prog_wtr(nnow_rcf(jg))%h_ice,                           &
        &                opt_checkgroup=initicon(jg)%sfc%grp_vars_fg(1:ngrp_vars_fg))


      ! tile based fields
      DO jt=1, ntiles_total + ntiles_water 

        CALL read_data_2d (filetype_fg(jg), fileID_fg(jg), 't_g', p_patch(jg)%n_patch_cells_g,     &
         &                p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index,                  &
         &                p_lnd_state(jg)%prog_lnd(nnow_rcf(jg))%t_g_t(:,:,jt),                    &
         &                opt_checkgroup=initicon(jg)%sfc%grp_vars_fg(1:ngrp_vars_fg))

       ! aggregated values are not prog variables, 
       ! aggregated values are calculated in  init_nwp_phy  (PR)
       ! p_lnd_state(jg)%prog_lnd(nnow_rcf(jg))%t_g(:,:) = &
       !   &    p_lnd_state(jg)%prog_lnd(nnow_rcf(jg))%t_g_t(:,:,1)

        CALL read_data_2d (filetype_fg(jg), fileID_fg(jg), 'qv_s', p_patch(jg)%n_patch_cells_g, &
         &                p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index,              &
         &                p_lnd_state(jg)%diag_lnd%qv_s_t(:,:,jt),                             &
         &                opt_checkgroup=initicon(jg)%sfc%grp_vars_fg(1:ngrp_vars_fg))

       ! aggregated values are calculated in  init_nwp_phy  (PR)
       ! p_lnd_state(jg)%diag_lnd%qv_s(:,:) = &
       !  &    p_lnd_state(jg)%diag_lnd%qv_s_t(:,:,1)
      END DO
        !  tile based fields
      DO jt=1, ntiles_total

        CALL read_data_2d (filetype_fg(jg), fileID_fg(jg), 'freshsnow',               &
          &                p_patch(jg)%n_patch_cells_g,                               &
          &                p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index,    &
          &                p_lnd_state(jg)%diag_lnd%freshsnow_t(:,:,jt),              &
          &                opt_checkgroup=initicon(jg)%sfc%grp_vars_fg(1:ngrp_vars_fg) )

        CALL read_data_2d (filetype_fg(jg), fileID_fg(jg), 'w_snow',                  &
          &                p_patch(jg)%n_patch_cells_g,                               &
          &                p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index,    & 
          &                p_lnd_state(jg)%prog_lnd(nnow_rcf(jg))%w_snow_t(:,:,jt),   &
          &                opt_checkgroup=initicon(jg)%sfc%grp_vars_fg(1:ngrp_vars_fg) )

        CALL read_data_2d (filetype_fg(jg), fileID_fg(jg), 'w_i',                     &
          &                p_patch(jg)%n_patch_cells_g,                               &
          &                p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index,    &
          &                p_lnd_state(jg)%prog_lnd(nnow_rcf(jg))%w_i_t(:,:,jt),      &
          &                opt_checkgroup=initicon(jg)%sfc%grp_vars_fg(1:ngrp_vars_fg) )

        CALL read_data_2d (filetype_fg(jg), fileID_fg(jg), 'h_snow',                  &
          &                p_patch(jg)%n_patch_cells_g,                               &
          &                p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index,    &
          &                p_lnd_state(jg)%diag_lnd%h_snow_t(:,:,jt),                 &
          &                opt_checkgroup=initicon(jg)%sfc%grp_vars_fg(1:ngrp_vars_fg) )


        CALL read_data_2d (filetype_fg(jg), fileID_fg(jg),'t_snow',                   &
          &                p_patch(jg)%n_patch_cells_g,                               &
          &                p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index,    &
          &                p_lnd_state(jg)%prog_lnd(nnow_rcf(jg))%t_snow_t(:,:,jt),   &
          &                opt_checkgroup=initicon(jg)%sfc%grp_vars_fg(1:ngrp_vars_fg) )

        CALL read_data_2d (filetype_fg(jg), fileID_fg(jg), 'rho_snow',                &
          &                p_patch(jg)%n_patch_cells_g,                               &
          &                p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index,    &
          &                p_lnd_state(jg)%prog_lnd(nnow_rcf(jg))%rho_snow_t(:,:,jt), &
          &                opt_checkgroup=initicon(jg)%sfc%grp_vars_fg(1:ngrp_vars_fg) )



     ! multi layer fields 
        CALL read_data_3d (filetype_fg(jg), fileID_fg(jg), 'w_so',                             &
          &                p_patch(jg)%n_patch_cells_g,                                        &
          &                p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index,             &
          &                nlev_soil, p_lnd_state(jg)%prog_lnd(nnow_rcf(jg))%w_so_t(:,:,:,jt), &
          &                opt_checkgroup=initicon(jg)%sfc%grp_vars_fg(1:ngrp_vars_fg) )


        ! so far w_so_ice is re-initialized in terra_multlay_init
        CALL read_data_3d (filetype_fg(jg), fileID_fg(jg), 'w_so_ice',                    &
          &                p_patch(jg)%n_patch_cells_g,                                   &
          &                p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index,        &
          &                nlev_soil, p_lnd_state(jg)%prog_lnd(nnow_rcf(jg))%w_so_ice_t(:,:,:,jt), &
          &                opt_checkgroup=initicon(jg)%sfc%grp_vars_fg(1:ngrp_vars_fg) )

        CALL read_data_3d (filetype_fg(jg), fileID_fg(jg), 't_so',                        &
          &                p_patch(jg)%n_patch_cells_g,                                   &
          &                p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index,        &
          &                nlev_soil+1, p_lnd_state(jg)%prog_lnd(nnow_rcf(jg))%t_so_t(:,:,:,jt), &
          &                opt_checkgroup=initicon(jg)%sfc%grp_vars_fg(1:ngrp_vars_fg))

      ENDDO ! jt


      ! Skipped in MODE_COMBINED (i.e. when starting from GME soil) 
      ! Instead z0 is re-initialized (see mo_nwp_phy_init)
      CALL read_data_2d (filetype_fg(jg), fileID_fg(jg), 'gz0',                         &
        &                p_patch(jg)%n_patch_cells_g,                                   &
        &                p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index,        &
        &                prm_diag(jg)%gz0(:,:),                                         &
        &                opt_checkgroup=initicon(jg)%sfc%grp_vars_fg(1:ngrp_vars_fg) )

    ENDDO ! loop over model domains


    !----------------------------------------!
    ! read in DWD analysis (surface)         !
    !----------------------------------------!

    IF (l_ana_sfc) THEN

      ! Analysis is read for the global domain, only.
      jg = 1

      ! save some paperwork
      ngrp_vars_fg  = initicon(jg)%sfc%ngrp_vars_fg
      ngrp_vars_ana = initicon(jg)%sfc%ngrp_vars_ana


      IF(p_pe == p_io ) THEN 
        CALL message (TRIM(routine), 'read sfc_ANA fields from '//TRIM(dwdana_file(jg)))

        SELECT CASE(filetype_ana(jg))
        CASE (FILETYPE_NC2)
          ! open file, read basic dimensions

          ! get number of cells
          CALL nf(nf_inq_dimid(fileID_ana(jg), 'ncells', dimid), routine)
          CALL nf(nf_inq_dimlen(fileID_ana(jg), dimid, no_cells), routine)
          ! get number of vertical levels
          CALL nf(nf_inq_dimid(fileID_ana(jg), 'lev', dimid), routine)
          CALL nf(nf_inq_dimlen(fileID_ana(jg), dimid, no_levels), routine)

          ! -------------------------
          ! simple consistency checks
          ! -------------------------
        
          ! check the number of cells and vertical levels
          IF( p_patch(jg)%n_patch_cells_g /= no_cells) THEN
            CALL finish(TRIM(ROUTINE),&
                 & 'Number of patch cells and cells in DWD ANA file do not match.')
          ENDIF
        
          IF(p_patch(jg)%nlev /= no_levels) THEN
            CALL finish(TRIM(ROUTINE),&
                 & 'nlev does not match the number of levels in DWD ANA file.')
          ENDIF

          !------------------------------------------
          ! Check if required variables are available
          !------------------------------------------
          ! To be done for netcdf

        CASE (FILETYPE_GRB2)

          !------------------------------------------
          ! Check if required variables are available
          !------------------------------------------
          ! To be done for GRIB2

        CASE DEFAULT
          CALL finish(routine, "Unknown file type")
        END SELECT

      ENDIF   ! p_io

      ! set tile-index explicitly
      jt = 1

      ! sea-ice fraction
      CALL read_data_2d (filetype_ana(jg), fileID_ana(jg), 'fr_seaice', p_patch(jg)%n_patch_cells_g,  &
        &                p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index,                      &
        &                p_lnd_state(jg)%diag_lnd%fr_seaice,                                          &
        &                opt_checkgroup=initicon(jg)%sfc%grp_vars_ana(1:ngrp_vars_ana) )

      ! sea-ice temperature
      CALL read_data_2d (filetype_ana(jg), fileID_ana(jg), 't_ice', p_patch(jg)%n_patch_cells_g,      &
        &                p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index,                      &
        &                p_lnd_state(jg)%prog_wtr(nnow_rcf(jg))%t_ice,                                &
        &                opt_checkgroup=initicon(jg)%sfc%grp_vars_ana(1:ngrp_vars_ana) )

      ! sea-ice height
      CALL read_data_2d (filetype_ana(jg), fileID_ana(jg), 'h_ice', p_patch(jg)%n_patch_cells_g,      &
        &                p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index,                      &
        &                p_lnd_state(jg)%prog_wtr(nnow_rcf(jg))%h_ice,                                &
        &                opt_checkgroup=initicon(jg)%sfc%grp_vars_ana(1:ngrp_vars_ana) )

      ! T_SO(0)
      CALL read_data_2d (filetype_ana(jg), fileID_ana(jg), 't_so', p_patch(jg)%n_patch_cells_g,       &
        &                p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index,                      &
        &                p_lnd_state(jg)%prog_lnd(nnow_rcf(jg))%t_so_t(:,1,:,jt),                     &
        &                opt_checkgroup=initicon(jg)%sfc%grp_vars_ana(1:ngrp_vars_ana) )

      ! h_snow
      CALL read_data_2d (filetype_ana(jg), fileID_ana(jg), 'h_snow',                &
        &                p_patch(jg)%n_patch_cells_g,                               &
        &                p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index,    &
        &                p_lnd_state(jg)%diag_lnd%h_snow_t(:,:,jt),                 &
        &                opt_checkgroup=initicon(jg)%sfc%grp_vars_ana(1:ngrp_vars_ana) )

      ! w_snow
      CALL read_data_2d (filetype_ana(jg), fileID_ana(jg), 'w_snow',                &
        &                p_patch(jg)%n_patch_cells_g,                               &
        &                p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index,    &
        &                p_lnd_state(jg)%prog_lnd(nnow_rcf(jg))%w_snow_t(:,:,jt),   &
        &                opt_checkgroup=initicon(jg)%sfc%grp_vars_ana(1:ngrp_vars_ana) )

      ! w_i
      CALL read_data_2d (filetype_ana(jg), fileID_ana(jg), 'w_i',                   &
        &                p_patch(jg)%n_patch_cells_g,                               &
        &                p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index,    &
        &                p_lnd_state(jg)%prog_lnd(nnow_rcf(jg))%w_i_t(:,:,jt),      &
        &                opt_checkgroup=initicon(jg)%sfc%grp_vars_ana(1:ngrp_vars_ana) )

      ! t_snow
      CALL read_data_2d (filetype_ana(jg), fileID_ana(jg), 't_snow',                &
        &                p_patch(jg)%n_patch_cells_g,                               &
        &                p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index,    &
        &                p_lnd_state(jg)%prog_lnd(nnow_rcf(jg))%t_snow_t(:,:,jt),   &
        &                opt_checkgroup=initicon(jg)%sfc%grp_vars_ana(1:ngrp_vars_ana) )

      ! rho_snow
      CALL read_data_2d (filetype_ana(jg), fileID_ana(jg), 'rho_snow',              &
        &                p_patch(jg)%n_patch_cells_g,                               &
        &                p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index,    &
        &                p_lnd_state(jg)%prog_lnd(nnow_rcf(jg))%rho_snow_t(:,:,jt), &
        &                opt_checkgroup=initicon(jg)%sfc%grp_vars_ana(1:ngrp_vars_ana) )

      ! freshsnow
      CALL read_data_2d (filetype_ana(jg), fileID_ana(jg), 'freshsnow',             &
        &                p_patch(jg)%n_patch_cells_g,                               &
        &                p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index,    &
        &                p_lnd_state(jg)%diag_lnd%freshsnow_t(:,:,jt),              &
        &                opt_checkgroup=initicon(jg)%sfc%grp_vars_ana(1:ngrp_vars_ana) )

      ! w_so
      CALL read_data_3d (filetype_ana(jg), fileID_ana(jg), 'w_so',                           &
        &                p_patch(jg)%n_patch_cells_g,                                        &
        &                p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index,             &
        &                nlev_soil, p_lnd_state(jg)%prog_lnd(nnow_rcf(jg))%w_so_t(:,:,:,jt), &
        &                opt_checkgroup=initicon(jg)%sfc%grp_vars_ana(1:ngrp_vars_ana) )

    ENDIF  ! l_ana_sfc


    ! Only required, when starting from GME soil (so far, W_SO=SMI*1000 in GME input file)
    ! Also note, that the domain-loop is missing
    IF (init_mode == MODE_COMBINED) THEN
      jg = 1
      DO jt=1, ntiles_total
        CALL smi_to_sm_mass(p_patch(jg), p_lnd_state(jg)%prog_lnd(nnow_rcf(jg))%w_so_t(:,:,:,jt))
      ENDDO
    END IF

  END SUBROUTINE read_dwdana_sfc




  !>
  !! Analysis is created by merging the first guess with the DA output 
  !!
  !!
  !! Analysis is created by merging the first guess with the DA output 
  !! (atmosphere only).
  !! First the FG in terms of the NH prognostic set of variables
  !! is converted into p, T, u and v.
  !! Then, increments are computed as the difference between the DA output and
  !! the converted dynamical variables, and then are transformed 
  !! back to the NH prognostic set of variables and are added to the first guess.
  !!
  !! Sanity check with FG only. If the analysis is set equal to the FG, the 
  !! increments should be exactly 0. It was verified, that the nonzero values 
  !! in the increment fields are due to the GRIB packing and not due to a 
  !! coding error. I.e. ibits was increased from DATATYPE_PACK16 to 
  !! DATATYPE_PACK32 and the errors in u_incr, v_incr went down from O(10E-3) 
  !! to O(10E-7). Similarly the error in pres_incr went down from O(1) to 
  !! O(1E-1).
  !! 
  !! @par Revision History
  !! Initial version by Daniel Reinert, DWD(2012-12-18)
  !!
  SUBROUTINE create_dwdana_atm (p_patch, p_nh_state, p_int_state)

    TYPE(t_patch),    TARGET, INTENT(IN)    :: p_patch(:)
    TYPE(t_nh_state), TARGET, INTENT(INOUT) :: p_nh_state(:)
    TYPE(t_int_state),        INTENT(IN)    :: p_int_state(:)

    INTEGER :: jc,je,jk,jb,jg             ! loop indices
    INTEGER :: ist
    INTEGER :: nlev, nlevp1               ! number of vertical levels
    INTEGER :: nblks_c, nblks_e           ! number of blocks
    INTEGER :: i_nchdom
    INTEGER :: rl_start, rl_end 
    INTEGER :: i_startblk, i_endblk
    INTEGER :: i_startidx, i_endidx
    TYPE(t_nh_prog), POINTER :: p_prog_now, p_prog_now_rcf   
    TYPE(t_nh_diag), POINTER :: p_diag
    INTEGER,         POINTER :: iidx(:,:,:), iblk(:,:,:)

    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: zpres_nh, pres_incr, u_incr, v_incr, vn_incr, &
                                               nabla4_vn_incr, w_incr
    REAL(wp) :: vn_incr_smt

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
      routine = 'mo_nh_initicon:create_dwdana_atm'

    ! nondimensional diffusion coefficient for interpolated velocity increment
    REAL(wp), PARAMETER :: smtfac=0.015_wp

    !-------------------------------------------------------------------------

    ! for the time being, the generation of DWD analysis fields is implemented 
    ! for the global domain, only.
    jg = 1

    ! number of vertical levels 
    nlev      = p_patch(jg)%nlev
    nlevp1    = p_patch(jg)%nlevp1

    nblks_c   = p_patch(jg)%nblks_c
    nblks_e   = p_patch(jg)%nblks_e
    i_nchdom  = MAX(1,p_patch(jg)%n_childdom)


    ! allocate temporary arrays for nonhydrostatic pressure, DA increments and a filtering term for vn
    ! note that an explicit temperature increment is not required (see below)
    ALLOCATE(zpres_nh (nproma,nlev,nblks_c),  &
             pres_incr(nproma,nlev,nblks_c),  &
             u_incr   (nproma,nlev,nblks_c),  &
             v_incr   (nproma,nlev,nblks_c),  &
             vn_incr  (nproma,nlev,nblks_e),  &
             w_incr   (nproma,nlevp1,nblks_c),&
             nabla4_vn_incr(nproma,nlev,nblks_e), &
             STAT=ist)
    IF (ist /= SUCCESS) THEN
      CALL finish ( TRIM(routine), 'allocation of auxiliary arrays failed')
    ENDIF

    nabla4_vn_incr(:,:,:) = 0._wp

    ! define some pointers
    p_prog_now     => p_nh_state(jg)%prog(nnow(jg))
    p_prog_now_rcf => p_nh_state(jg)%prog(nnow_rcf(jg))
    p_diag         => p_nh_state(jg)%diag
    iidx           => p_patch(jg)%edges%cell_idx
    iblk           => p_patch(jg)%edges%cell_blk

    ! Recompute u and v from the first guess in order to compute the wind increment
    ! coming from the data assimilation
    CALL rbf_vec_interpol_cell(p_prog_now%vn, p_patch(jg), p_int_state(jg), p_diag%u, p_diag%v)

    ! 1) first guess in terms of rho, theta_v, qx is converted to 
    ! T, p, qx. Note, that zpres_nh is the full (nonhydrostatic) pressure field, whereas
    ! p_diag%pres is the hydrostatically integrated pressure field
    !
    ! Note that the diagnose_pres_temp routine cannot be used at this moment because
    ! the exner function still needs to be computed
    !
!$OMP PARALLEL PRIVATE(rl_start,rl_end,i_startblk,i_endblk)

    ! include boundary interpolation zone of nested domains and halo points
    rl_start = 1
    rl_end   = min_rlcell

    i_startblk = p_patch(jg)%cells%start_blk(rl_start,1)
    i_endblk   = p_patch(jg)%cells%end_blk(rl_end,i_nchdom)

!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx)
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch(jg), jb, i_startblk, i_endblk, &
        & i_startidx, i_endidx, rl_start, rl_end)

      DO jk = 1, nlev
        DO jc = i_startidx, i_endidx

          ! compute exner function
          p_prog_now%exner(jc,jk,jb) = (rd/p0ref * p_prog_now%rho(jc,jk,jb)  &
            &                        * p_prog_now%theta_v(jc,jk,jb))**(rd/cvd)

          ! compute full nonhydrostatic pressure
          zpres_nh(jc,jk,jb) = p_prog_now%exner(jc,jk,jb)**(cpd/rd) * p0ref

          ! compute virtual temperature
          p_diag%tempv(jc,jk,jb) = p_prog_now%theta_v(jc,jk,jb) &
            &                    * p_prog_now%exner(jc,jk,jb)

          ! compute temperature (currently unused - but we could use it to check the hydrostatic
          ! balance of the DA increments)
          p_diag%temp(jc,jk,jb) = p_diag%tempv(jc,jk,jb)  &
            &                   / (1._wp + vtmpc1*p_prog_now_rcf%tracer(jc,jk,jb,iqv) &
            &                   - (p_prog_now_rcf%tracer(jc,jk,jb,iqc)                &
            &                   +  p_prog_now_rcf%tracer(jc,jk,jb,iqi)                &
            &                   +  p_prog_now_rcf%tracer(jc,jk,jb,iqr)                &
            &                   +  p_prog_now_rcf%tracer(jc,jk,jb,iqs)) )

        ENDDO  ! jc
      ENDDO  ! jk

    ENDDO  ! jb
!$OMP END DO
!$OMP END PARALLEL

    ! Recompute the hydrostatically integrated pressure from the first guess
    CALL diagnose_pres_temp (p_nh_state(jg)%metrics, p_prog_now, p_prog_now_rcf, p_diag, p_patch(jg), &
      &                      opt_calc_temp=.FALSE., opt_calc_pres=.TRUE.)


    ! 2) compute DA increments and add them to the first guess
    !
!$OMP PARALLEL PRIVATE(rl_start,rl_end,i_startblk,i_endblk)

    ! include boundary interpolation zone of nested domains and halo points
    rl_start = 1
    rl_end   = min_rlcell

    i_startblk = p_patch(jg)%cells%start_blk(rl_start,1)
    i_endblk   = p_patch(jg)%cells%end_blk(rl_end,i_nchdom)

!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx)
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch(jg), jb, i_startblk, i_endblk, &
        & i_startidx, i_endidx, rl_start, rl_end)

      DO jk = 1, nlev
        DO jc = i_startidx, i_endidx

          ! pressure increment - should we verify that it is in hydrostatic balance with 
          ! the temperature increment?
          pres_incr(jc,jk,jb) = initicon(jg)%atm%pres(jc,jk,jb) - p_diag%pres(jc,jk,jb)

          ! increments for u and v - will be interpolated to edge points below
          u_incr(jc,jk,jb) = initicon(jg)%atm%u(jc,jk,jb) - p_diag%u(jc,jk,jb)
          v_incr(jc,jk,jb) = initicon(jg)%atm%v(jc,jk,jb) - p_diag%v(jc,jk,jb)

          ! add pressure increment to the nonhydrostatic pressure
          zpres_nh(jc,jk,jb) = zpres_nh(jc,jk,jb) + pres_incr(jc,jk,jb)

          ! temperature increment is not needed explicitly. Note that lateron the analysed 
          ! temperature field initicon(jg)%atm%temp, instead of the first guess 
          ! temperature field p_diag%temp is used to compute the virtual temperature 
          ! and lateron the virtual potential temperature.

        ENDDO  ! jc
      ENDDO  ! jk

    ENDDO  ! jb
!$OMP END DO


    ! include boundary interpolation zone of nested domains and the halo edges as far as possible
    rl_start = 2
    rl_end   = min_rledge_int - 2

    i_startblk = p_patch(jg)%edges%start_blk(rl_start,1)
    i_endblk   = p_patch(jg)%edges%end_blk(rl_end,i_nchdom)

!$OMP DO PRIVATE(jb,jk,je,i_startidx,i_endidx)
    DO jb = i_startblk, i_endblk

      CALL get_indices_e(p_patch(jg), jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, rl_start, rl_end)

      DO jk = 1, nlev
        DO je = i_startidx, i_endidx
          ! at cell centers the increment \vec(v_inc) is projected into the 
          ! direction of vn and then linearly interpolated to the edge midpoint 
          !
          ! should we check if the vn increments are geostrophically balanced at higher levels?
          vn_incr(je,jk,jb) = p_int_state(jg)%c_lin_e(je,1,jb)                  &
            &               *(u_incr(iidx(je,jb,1),jk,iblk(je,jb,1))            &
            &               * p_patch(jg)%edges%primal_normal_cell(je,jb,1)%v1  &
            &               + v_incr(iidx(je,jb,1),jk,iblk(je,jb,1))            &
            &               * p_patch(jg)%edges%primal_normal_cell(je,jb,1)%v2) &
            &               + p_int_state(jg)%c_lin_e(je,2,jb)                  &
            &               *(u_incr(iidx(je,jb,2),jk,iblk(je,jb,2))            & 
            &               * p_patch(jg)%edges%primal_normal_cell(je,jb,2)%v1  &
            &               + v_incr(iidx(je,jb,2),jk,iblk(je,jb,2))            &
            &               * p_patch(jg)%edges%primal_normal_cell(je,jb,2)%v2  )


        ENDDO  ! je
      ENDDO  ! jk

    ENDDO  ! jb
!$OMP ENDDO
!$OMP END PARALLEL

    !DR required to avoid crash in nabla4_vec
    CALL sync_patch_array(SYNC_E,p_patch(jg),vn_incr)

    ! Compute diffusion term 
    CALL nabla4_vec(vn_incr, p_patch(jg), p_int_state(jg), nabla4_vn_incr, opt_rlstart=5)

    ! Compute vertical wind increment consistent with the vn increment
    ! (strictly spoken, this should be done after the filtering step, but the difference is negligible)
    CALL init_w(p_patch(jg), p_int_state(jg), vn_incr, p_nh_state(jg)%metrics%z_ifc, w_incr)

!$OMP PARALLEL PRIVATE(rl_start,rl_end,i_startblk,i_endblk)

    ! include boundary interpolation zone of nested domains but no halo points (sync follows below)
    rl_start = 2
    rl_end   = min_rledge_int

    i_startblk = p_patch(jg)%edges%start_blk(rl_start,1)
    i_endblk   = p_patch(jg)%edges%end_blk(rl_end,i_nchdom)

!$OMP DO PRIVATE(jb,jk,je,i_startidx,i_endidx,vn_incr_smt)
    DO jb = i_startblk, i_endblk

      CALL get_indices_e(p_patch(jg), jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, rl_start, rl_end)

      DO jk = 1, nlev
        DO je = i_startidx, i_endidx
          ! computed filtered velocity increment 
          vn_incr_smt = vn_incr(je,jk,jb) - smtfac*nabla4_vn_incr(je,jk,jb)*p_patch(jg)%edges%area_edge(je,jb)**2

          ! add vn_incr_smt to first guess
          p_prog_now%vn(je,jk,jb) = p_prog_now%vn(je,jk,jb) + vn_incr_smt

        ENDDO  ! je
      ENDDO  ! jk

    ENDDO  ! jb
!$OMP ENDDO

    ! include boundary interpolation zone of nested domains but no halo points (sync follows below)
    rl_start = 2
    rl_end   = min_rlcell_int

    i_startblk = p_patch(jg)%cells%start_blk(rl_start,1)
    i_endblk   = p_patch(jg)%cells%end_blk(rl_end,i_nchdom)

!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx)
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch(jg), jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, rl_start, rl_end)

      DO jk = 1, nlevp1
        DO jc = i_startidx, i_endidx

          ! add w_incr to first guess
          p_prog_now%w(jc,jk,jb) = p_prog_now%w(jc,jk,jb) + w_incr(jc,jk,jb) 

        ENDDO  ! jc
      ENDDO  ! jk

    ENDDO  ! jb
!$OMP ENDDO
!$OMP END PARALLEL

    CALL sync_patch_array(SYNC_E,p_patch(jg),p_prog_now%vn)
    CALL sync_patch_array(SYNC_C,p_patch(jg),p_prog_now%w)


    ! TO DO: remove qc, where rh<90%


    ! 3) Convert analysis back to the NH set of prognostic variables
    !

    ! Compute virtual temperature
    CALL virtual_temp(p_patch(jg), initicon(jg)%atm%temp, & !in
      &               p_prog_now%tracer(:,:,:,iqv),       & !in
      &               p_prog_now%tracer(:,:,:,iqc),       & !in
      &               p_prog_now%tracer(:,:,:,iqi),       & !in
      &               p_prog_now%tracer(:,:,:,iqr),       & !in
      &               p_prog_now%tracer(:,:,:,iqs),       & !in
      &               p_diag%tempv                        ) !out  

    ! Convert thermodynamic variables into set of NH prognostic variables
    CALL convert_thdvars(p_patch(jg), zpres_nh,  & !in
      &                  p_diag%tempv,           & !in
      &                  p_prog_now%rho,         & !out
      &                  p_prog_now%exner,       & !out
      &                  p_prog_now%theta_v      ) !out


!$OMP PARALLEL
!$OMP WORKSHARE
    ! compute prognostic variable rho*theta_v
    p_prog_now%rhotheta_v(:,:,:) = p_prog_now%rho(:,:,:)   &
      &                          * p_prog_now%theta_v(:,:,:)
!$OMP END WORKSHARE
!$OMP END PARALLEL

    ! deallocate temporary arrays
    DEALLOCATE( zpres_nh, pres_incr, u_incr, v_incr, vn_incr, nabla4_vn_incr, w_incr, STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ( TRIM(routine), 'deallocation of auxiliary arrays failed' )
    ENDIF

  END SUBROUTINE create_dwdana_atm



  !-------------------------------------------------------------------------
  !>
  !! SUBROUTINE create_dwdana_sfc 
  !!
  !! Required input: patch, lnd_state
  !! Output is written on fields of NH state
  !!
  !! @par Revision History
  !! Initial version by P. Ripodas, DWD(2013-05)
  !!
  !!
  !-------------------------------------------------------------------------
  SUBROUTINE create_dwdana_sfc (p_patch, p_lnd_state, ext_data)

    TYPE(t_patch),    TARGET ,INTENT(IN)    :: p_patch(:)
    TYPE(t_lnd_state)        ,INTENT(INOUT) :: p_lnd_state(:)
    TYPE(t_external_data)    ,INTENT(INOUT) :: ext_data(:)

    INTEGER :: jg, ic, jc, jb             ! loop indices
    INTEGER :: ntlr
    INTEGER :: nblks_c   
  !-------------------------------------------------------------------------

    ! for the time being, the generation of DWD analysis fields is implemented 
    ! for the global domain, only.
    jg = 1

    nblks_c   = p_patch(jg)%nblks_c
    ntlr      = nnow_rcf(jg)
!$OMP PARALLEL 
!$OMP DO PRIVATE(jc,ic,jb) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = 1, nblks_c

        !get SST from first soil level t_so (for sea and lake points)
!CDIR NODEP,VOVERTAKE,VOB
        DO ic = 1, ext_data(jg)%atm%sp_count(jb)
           jc = ext_data(jg)%atm%idx_lst_sp(ic,jb)
           p_lnd_state(jg)%diag_lnd%t_seasfc(jc,jb) =                  & ! nproma.nlev_soil+1,nblks,ntiles_total
                                    & p_lnd_state(jg)%prog_lnd(ntlr)%t_so_t(jc,1,jb,1) 
        END DO
!CDIR NODEP,VOVERTAKE,VOB
        DO ic = 1, ext_data(jg)%atm%fp_count(jb)
          jc = ext_data(jg)%atm%idx_lst_fp(ic,jb)
          p_lnd_state(jg)%diag_lnd%t_seasfc(jc,jb) =                   &
                                    & p_lnd_state(jg)%prog_lnd(ntlr)%t_so_t(jc,1,jb,1) 
        END DO
       END DO
!$OMP END DO
!$OMP END PARALLEL


  END SUBROUTINE create_dwdana_sfc
  !-------------------------------------------------------------------------


  !>
  !! SUBROUTINE copy_initicon2prog_atm
  !! Copies atmospheric fields interpolated by init_icon to the
  !! prognostic model state variables 
  !!
  !! Required input: initicon state
  !! Output is written on fields of NH state
  !!
  !! @par Revision History
  !! Initial version by Guenther Zaengl, DWD(2011-07-28)
  !! Modification by Daniel Reinert, DWD (2012-12-19)
  !! - encapsulated surface specific part
  !!
  !!
  SUBROUTINE copy_initicon2prog_atm(p_patch, initicon, p_nh_state)

    TYPE(t_patch),          INTENT(IN) :: p_patch(:)
    TYPE(t_initicon_state), INTENT(IN) :: initicon(:)

    TYPE(t_nh_state),      INTENT(INOUT) :: p_nh_state(:)

    INTEGER :: jg, jb, jk, jc, je
    INTEGER :: nblks_c, npromz_c, nblks_e, npromz_e, nlen, nlev, nlevp1, ntl, ntlr

!$OMP PARALLEL PRIVATE(jg,nblks_c,npromz_c,nblks_e,npromz_e,nlev,nlevp1,ntl,ntlr)
    DO jg = 1, n_dom

      IF (.NOT. p_patch(jg)%ldom_active) CYCLE

      nblks_c   = p_patch(jg)%nblks_c
      npromz_c  = p_patch(jg)%npromz_c
      nblks_e   = p_patch(jg)%nblks_e
      npromz_e  = p_patch(jg)%npromz_e
      nlev      = p_patch(jg)%nlev
      nlevp1    = p_patch(jg)%nlevp1
      ntl       = nnow(jg)
      ntlr      = nnow_rcf(jg)

!$OMP DO PRIVATE(jb,jk,je,nlen) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = 1, nblks_e

        IF (jb /= nblks_e) THEN
          nlen = nproma
        ELSE
          nlen = npromz_e
        ENDIF

        ! Wind speed
        DO jk = 1, nlev
          DO je = 1, nlen
            p_nh_state(jg)%prog(ntl)%vn(je,jk,jb) = initicon(jg)%atm%vn(je,jk,jb)
          ENDDO
        ENDDO

      ENDDO
!$OMP END DO

!$OMP DO PRIVATE(jb,jk,jc,nlen) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = 1, nblks_c

        IF (jb /= nblks_c) THEN
          nlen = nproma
        ELSE
          nlen = npromz_c
        ENDIF

        ! 3D fields
        DO jk = 1, nlev
          DO jc = 1, nlen
            ! Dynamic prognostic variables on cell points
            p_nh_state(jg)%prog(ntl)%w(jc,jk,jb)       = initicon(jg)%atm%w(jc,jk,jb)
            p_nh_state(jg)%prog(ntl)%theta_v(jc,jk,jb) = initicon(jg)%atm%theta_v(jc,jk,jb)
            p_nh_state(jg)%prog(ntl)%exner(jc,jk,jb)   = initicon(jg)%atm%exner(jc,jk,jb)
            p_nh_state(jg)%prog(ntl)%rho(jc,jk,jb)     = initicon(jg)%atm%rho(jc,jk,jb)

            p_nh_state(jg)%prog(ntl)%rhotheta_v(jc,jk,jb) = &
              p_nh_state(jg)%prog(ntl)%rho(jc,jk,jb) * p_nh_state(jg)%prog(ntl)%theta_v(jc,jk,jb)

            ! Water vapor
            p_nh_state(jg)%prog(ntlr)%tracer(jc,jk,jb,iqv) = initicon(jg)%atm%qv(jc,jk,jb)
          ENDDO
        ENDDO

        ! Cloud and precipitation hydrometeors - these are supposed to be zero in the region where
        ! moisture physics is turned off
        DO jk = kstart_moist(jg), nlev
          DO jc = 1, nlen
            p_nh_state(jg)%prog(ntlr)%tracer(jc,jk,jb,iqc) = initicon(jg)%atm%qc(jc,jk,jb)
            p_nh_state(jg)%prog(ntlr)%tracer(jc,jk,jb,iqi) = initicon(jg)%atm%qi(jc,jk,jb)
            p_nh_state(jg)%prog(ntlr)%tracer(jc,jk,jb,iqr) = initicon(jg)%atm%qr(jc,jk,jb)
            p_nh_state(jg)%prog(ntlr)%tracer(jc,jk,jb,iqs) = initicon(jg)%atm%qs(jc,jk,jb)
          ENDDO
        ENDDO

        ! w at surface level
        DO jc = 1, nlen
          p_nh_state(jg)%prog(ntl)%w(jc,nlevp1,jb)      = initicon(jg)%atm%w(jc,nlevp1,jb)
          p_nh_state(jg)%prog(nnew(jg))%w(jc,nlevp1,jb) = initicon(jg)%atm%w(jc,nlevp1,jb)
        ENDDO

      ENDDO  ! jb
!$OMP END DO NOWAIT

    ENDDO  ! jg
!$OMP END PARALLEL

    ! Finally, compute exact hydrostatic adjustment for thermodynamic fields
    DO jg = 1, n_dom

      IF (.NOT. p_patch(jg)%ldom_active) CYCLE
      ntl = nnow(jg)

      CALL hydro_adjust(p_patch(jg), p_nh_state(jg)%metrics,                                  &
                        p_nh_state(jg)%prog(ntl)%rho,     p_nh_state(jg)%prog(ntl)%exner,     &
                        p_nh_state(jg)%prog(ntl)%theta_v, p_nh_state(jg)%prog(ntl)%rhotheta_v )

    ENDDO

  END SUBROUTINE copy_initicon2prog_atm




  !-------------
  !>
  !! SUBROUTINE copy_initicon2prog_sfc
  !! Copies surface fields interpolated by init_icon to the prognostic model 
  !! state variables. 
  !!
  !! Required input: initicon state
  !! Output is written on fields of land state
  !!
  !! @par Revision History
  !! Initial version by Guenther Zaengl, DWD(2011-07-28)
  !! Modification by Daniel Reinert, DWD (2012-12-19)
  !! - encapsulated surface specific part
  !! Modification by Daniel Reinert, DWD (2013-07-09)
  !! - moved sea-ice coldstart into separate subroutine
  !!
  !!
  SUBROUTINE copy_initicon2prog_sfc(p_patch, initicon, p_lnd_state, ext_data)

    TYPE(t_patch),          INTENT(IN) :: p_patch(:)
    TYPE(t_initicon_state), INTENT(IN) :: initicon(:)

    TYPE(t_lnd_state),     INTENT(INOUT) :: p_lnd_state(:)
    TYPE(t_external_data), INTENT(   IN) :: ext_data(:)

    INTEGER  :: jg, jb, jc, jt, js, jp, ic
    INTEGER  :: nblks_c, npromz_c, nlen, nlev
    REAL(wp) :: zfrice_thrhld

!$OMP PARALLEL PRIVATE(jg,nblks_c,npromz_c,nlev)
    DO jg = 1, n_dom

      IF (.NOT. p_patch(jg)%ldom_active) CYCLE

      nblks_c   = p_patch(jg)%nblks_c
      npromz_c  = p_patch(jg)%npromz_c
      nlev      = p_patch(jg)%nlev


!$OMP DO PRIVATE(jb,jc,nlen,jt,js,jp,ic,zfrice_thrhld) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = 1, nblks_c

        IF (jb /= nblks_c) THEN
          nlen = nproma
        ELSE
          nlen = npromz_c
        ENDIF


        ! ground temperature
        DO jc = 1, nlen
          p_lnd_state(jg)%prog_lnd(nnow_rcf(jg))%t_g(jc,jb) = initicon(jg)%sfc%tskin(jc,jb)
          p_lnd_state(jg)%prog_lnd(nnew_rcf(jg))%t_g(jc,jb) = initicon(jg)%sfc%tskin(jc,jb)
        ENDDO

        ! Fill also SST and sea ice fraction fields over ocean points; SST is limited to 30 deg C
        ! Note: missing values of the sea ice fraction, which may occur due to differing land-sea masks, 
        ! are indicated with -999.9; non-ocean points are filled with zero for both fields
!CDIR NODEP,VOVERTAKE,VOB
        DO ic = 1, ext_data(jg)%atm%sp_count(jb)
          jc = ext_data(jg)%atm%idx_lst_sp(ic,jb)
          IF ( l_sst_in .AND. initicon(jg)%sfc%sst(jc,jb) > 10._wp  ) THEN
            p_lnd_state(jg)%diag_lnd%t_seasfc(jc,jb) = initicon(jg)%sfc%sst(jc,jb)              
          ELSE
           p_lnd_state(jg)%diag_lnd%t_seasfc(jc,jb) = MIN(303.15_wp,initicon(jg)%sfc%tskin(jc,jb))
          ENDIF
          !
          ! In case of missing sea ice fraction values, we make use of the sea 
          ! surface temperature (tskin over ocean points). For tskin<=tf_salt, 
          ! we set the sea ice fraction to one. For tskin>tf_salt, we set it to 0.
          ! Note: tf_salt=271.45K is the salt-water freezing point
          !
          IF ( initicon(jg)%sfc%seaice(jc,jb) > -999.0_wp ) THEN
            p_lnd_state(jg)%diag_lnd%fr_seaice(jc,jb) = initicon(jg)%sfc%seaice(jc,jb) 
          ELSE    ! missing value
            IF ( initicon(jg)%sfc%tskin(jc,jb) <= tf_salt ) THEN
              p_lnd_state(jg)%diag_lnd%fr_seaice(jc,jb) = 1._wp     ! sea ice point
            ELSE
              p_lnd_state(jg)%diag_lnd%fr_seaice(jc,jb) = 0._wp     ! water point
            ENDIF
          ENDIF

          ! For fr_seaice in ]0,frsi_min[, set fr_seaice to 0
          ! For fr_seaice in ]1-frsi_min,1[, set fr_seaice to 1. This will ensure in 
          ! init_sea_lists, that sea-ice and water fractions sum up exactly to the total 
          ! sea fraction.
          IF (p_lnd_state(jg)%diag_lnd%fr_seaice(jc,jb) < frsi_min ) THEN
             p_lnd_state(jg)%diag_lnd%fr_seaice(jc,jb) = 0._wp
          ENDIF
          IF (p_lnd_state(jg)%diag_lnd%fr_seaice(jc,jb) > (1._wp-frsi_min) ) THEN
             p_lnd_state(jg)%diag_lnd%fr_seaice(jc,jb) = 1._wp
          ENDIF

        ENDDO
        ! In addition, write skin temperature to lake points, limited to 33 deg C. These will
        ! be used to initialize lake points until something more reasonable becomes available
!CDIR NODEP,VOVERTAKE,VOB
        DO ic = 1, ext_data(jg)%atm%fp_count(jb)
          jc = ext_data(jg)%atm%idx_lst_fp(ic,jb)
          p_lnd_state(jg)%diag_lnd%t_seasfc(jc,jb) = MIN(306.15_wp,initicon(jg)%sfc%tskin(jc,jb))
        ENDDO

        IF ( atm_phy_nwp_config(jg)%inwp_surface > 0 ) THEN
          DO jt = 1, ntiles_total
            DO jc = 1, nlen
              p_lnd_state(jg)%prog_lnd(nnow_rcf(jg))%t_snow_t(jc,jb,jt)           = &
                &                                                initicon(jg)%sfc%tsnow   (jc,jb)

               ! Initialize freshsnow
               ! for seapoints, freshsnow is set to 0
              IF(alb_snow_var == 'ALB_SNOW') THEN
                p_lnd_state(jg)%diag_lnd%freshsnow_t(jc,jb,jt)    =  MAX(0._wp,MIN(1._wp, &
            &                           (initicon(jg)%sfc%snowalb (jc,jb)-csalb_snow_min) &
            &                          /(csalb_snow_max-csalb_snow_min)))                 &
            &                          * REAL(NINT(ext_data(jg)%atm%fr_land(jc,jb)),wp) 
              ELSE
                p_lnd_state(jg)%diag_lnd%freshsnow_t(jc,jb,jt)    =  MAX(0._wp,MIN(1._wp, &
            &                     1._wp - ((initicon(jg)%sfc%snowalb (jc,jb)-crhosmin_ml) &
            &                    /(crhosmax_ml-crhosmin_ml))))                            &
            &                    * REAL(NINT(ext_data(jg)%atm%fr_land(jc,jb)),wp)
              END IF


              p_lnd_state(jg)%prog_lnd(nnow_rcf(jg))%w_snow_t(jc,jb,jt)           = &
                &                                                initicon(jg)%sfc%snowweq (jc,jb)
              p_lnd_state(jg)%prog_lnd(nnow_rcf(jg))%rho_snow_t(jc,jb,jt)         = &
                &                                                initicon(jg)%sfc%snowdens(jc,jb) 
              p_lnd_state(jg)%prog_lnd(nnow_rcf(jg))%w_i_t(jc,jb,jt)              = &
                &                                                initicon(jg)%sfc%skinres (jc,jb)

              p_lnd_state(jg)%prog_lnd(nnew_rcf(jg))%t_snow_t(jc,jb,jt)   = &
                &                                                initicon(jg)%sfc%tsnow   (jc,jb)
              p_lnd_state(jg)%prog_lnd(nnew_rcf(jg))%w_snow_t(jc,jb,jt)   = &
                &                                                initicon(jg)%sfc%snowweq (jc,jb)
              p_lnd_state(jg)%prog_lnd(nnew_rcf(jg))%rho_snow_t(jc,jb,jt) = &
                &                                                initicon(jg)%sfc%snowdens(jc,jb)
              p_lnd_state(jg)%prog_lnd(nnew_rcf(jg))%w_i_t(jc,jb,jt)      = &
                &                                                initicon(jg)%sfc%skinres (jc,jb)
            ENDDO
          ENDDO

          ! Multi-layer surface fields
          DO jt = 1, ntiles_total

            DO js = 0, nlev_soil
              jp = js+1 ! indexing for the ICON state field starts at 1
              DO jc = 1, nlen
                p_lnd_state(jg)%prog_lnd(nnow_rcf(jg))%t_so_t(jc,jp,jb,jt)= &
                  &                                              initicon(jg)%sfc%tsoil(jc,jb,js)
                p_lnd_state(jg)%prog_lnd(nnew_rcf(jg))%t_so_t(jc,jp,jb,jt)= &
                  &                                              initicon(jg)%sfc%tsoil(jc,jb,js)
              ENDDO
            ENDDO

            ! For soil water, no comparable layer shift exists
            DO js = 1, nlev_soil
              DO jc = 1, nlen
                p_lnd_state(jg)%prog_lnd(nnow_rcf(jg))%w_so_t(jc,js,jb,jt)= &
                  &                                              initicon(jg)%sfc%wsoil(jc,jb,js)
                p_lnd_state(jg)%prog_lnd(nnew_rcf(jg))%w_so_t(jc,js,jb,jt)= &
                  &                                              initicon(jg)%sfc%wsoil(jc,jb,js)
              ENDDO
            ENDDO

          ENDDO



          ! Coldstart for sea-ice parameterization scheme
          ! Sea-ice surface temperature is initialized with tskin from IFS.
          ! Since the seaice index list is not yet available at this stage, we loop over 
          ! all sea points and initialize points with fr_seaice>threshold. 
          ! The threshold is 0.5 without tiles and frsi_min with tiles.
          ! Note that exactly the same threshold values must be used as in init_sea_lists. 
          ! If not, you will see what you get.
          !@Pilar: This should still work out for you, since the non-sea-ice points are 
          !        now initialized during warmstart initialization 
          !        in mo_nwp_sfc_utils:nwp_surface_init
          !
          IF (lseaice) THEN
            IF ( ntiles_total == 1 ) THEN  ! no tile approach
              zfrice_thrhld = 0.5_wp
            ELSE
              zfrice_thrhld = frsi_min
            ENDIF

            CALL seaice_coldinit_nwp(nproma, zfrice_thrhld,                               &
              &         frsi    = p_lnd_state(jg)%diag_lnd%fr_seaice(:,jb),               &
              &         temp_in = initicon(jg)%sfc%tskin(:,jb),                           &
              &         tice_p  = p_lnd_state(jg)%prog_wtr(nnow_rcf(jg))%t_ice    (:,jb), &
              &         hice_p  = p_lnd_state(jg)%prog_wtr(nnow_rcf(jg))%h_ice    (:,jb), &
              &         tsnow_p = p_lnd_state(jg)%prog_wtr(nnow_rcf(jg))%t_snow_si(:,jb), &
              &         hsnow_p = p_lnd_state(jg)%prog_wtr(nnow_rcf(jg))%h_snow_si(:,jb), &
              &         tice_n  = p_lnd_state(jg)%prog_wtr(nnew_rcf(jg))%t_ice    (:,jb), &
              &         hice_n  = p_lnd_state(jg)%prog_wtr(nnew_rcf(jg))%h_ice    (:,jb), &
              &         tsnow_n = p_lnd_state(jg)%prog_wtr(nnew_rcf(jg))%t_snow_si(:,jb), &
              &         hsnow_n = p_lnd_state(jg)%prog_wtr(nnew_rcf(jg))%h_snow_si(:,jb)  )

          ENDIF  ! leseaice

          ! Cold-start initialization of the fresh-water lake model FLake.
          ! The procedure is the same as in "int2lm".
          ! Note that no lake ice is assumed at the cold start.

!DR       ! To be implemented ... CALL flake_coldstart()
          ! Make use of sfc%ls_mask in order to identify potentially problematic points, 
          ! where depth_lk>0 (lake point in ICON) but ls_mask >0.5 (land point in IFS).
!DR       ! At these points, tskin should not be used to initialize the water temperature.

          IF (llake) THEN
            ! Loop over grid boxes with lakes
            DO ic=1, ext_data(jg)%atm%fp_count(jb)
              ! Take index from the lake index list
              jc = ext_data(jg)%atm%idx_lst_fp(ic,jb)
              ! Mixed-layer depth is set equal to 10 m or to the lake depth, whichever is smaller 
              p_lnd_state(jg)%prog_wtr(nnow_rcf(jg))%h_ml_lk(jc,jb) =  & 
                &  MIN(ext_data(jg)%atm%depth_lk(jc,jb), 10._wp)
              ! Shape factor is set to its minimum value 
              p_lnd_state(jg)%prog_wtr(nnow_rcf(jg))%c_t_lk(jc,jb)  = C_T_min             
              ! Mixed-layer temperature is set to the surface temprature
              ! that should be provided for the cold-start initialization 
              ! (e.g. skin temperature can be taken as an estimate of t_wml_lk)
              ! Clearly, the mixed-layer temperature cannot be less than the fresh-water freezing point 
              p_lnd_state(jg)%prog_wtr(nnow_rcf(jg))%t_wml_lk(jc,jb) =   &
                &  MAX(initicon(jg)%sfc%tskin(jc,jb), tpl_T_f)
              ! Bottom temperature is set equal 
              ! to the mixed-layer temperarture if the water column is mixed
              ! or to the temperarture of maximum density of fresh water  
              ! if the mixed-layer depth is less than the depth to the bottom
              p_lnd_state(jg)%prog_wtr(nnow_rcf(jg))%t_bot_lk(jc,jb) = MERGE(  &
                &  p_lnd_state(jg)%prog_wtr(nnow_rcf(jg))%t_wml_lk(jc,jb), tpl_T_r,  &
                &  p_lnd_state(jg)%prog_wtr(nnow_rcf(jg))%h_ml_lk(jc,jb) >=          & 
                &  ext_data(jg)%atm%depth_lk(jc,jb) ) 
              ! Mean temperature of the water column is computed from the formula
              ! t_mnw_lk = t_wml_lk - c_t_lk * MAX(0., 1.-h_ml_lk/depth_lk) * (t_wml_lk-t_bot_lk) 
              p_lnd_state(jg)%prog_wtr(nnow_rcf(jg))%t_mnw_lk(jc,jb) =                        &
                &  p_lnd_state(jg)%prog_wtr(nnow_rcf(jg))%t_wml_lk(jc,jb) -                   & 
                &  p_lnd_state(jg)%prog_wtr(nnow_rcf(jg))%c_t_lk(jc,jb) *                     &
                &  MAX(0._wp, ( 1._wp-p_lnd_state(jg)%prog_wtr(nnow_rcf(jg))%h_ml_lk(jc,jb)/  &  
                &  ext_data(jg)%atm%depth_lk(jc,jb))) *                                       &
                &  ( p_lnd_state(jg)%prog_wtr(nnow_rcf(jg))%t_wml_lk(jc,jb) -                 &
                &    p_lnd_state(jg)%prog_wtr(nnow_rcf(jg))%t_bot_lk(jc,jb) )
              ! No ice is assumed at cold-start initilization
              p_lnd_state(jg)%prog_wtr(nnow_rcf(jg))%t_ice(jc,jb) = tpl_T_f  ! Fresh-water freezing point
              p_lnd_state(jg)%prog_wtr(nnow_rcf(jg))%h_ice(jc,jb) = 0._wp    ! Zero ice thickness
              ! Snow over lake ice is not considered explicitly
              p_lnd_state(jg)%prog_wtr(nnow_rcf(jg))%t_snow_lk(jc,jb) =  &     ! Ice surface temperature
                &  p_lnd_state(jg)%prog_wtr(nnow_rcf(jg))%t_ice(jc,jb) 
              p_lnd_state(jg)%prog_wtr(nnow_rcf(jg))%h_snow_lk(jc,jb) = 0._wp  ! Zero snow thickness
              ! Bottom sediment module of FLake is switched off  
              p_lnd_state(jg)%prog_wtr(nnow_rcf(jg))%t_b1_lk(jc,jb) = tpl_T_r            ! Reference value
              p_lnd_state(jg)%prog_wtr(nnow_rcf(jg))%h_b1_lk(jc,jb) = rflk_depth_bs_ref  ! Reference value
              ! Set lake-surface temperature equal to the mixed-layer temperature
              p_lnd_state(jg)%prog_lnd(nnow_rcf(jg))%t_g_t(jc,jb,isub_lake) =  &
                &  p_lnd_state(jg)%prog_wtr(nnow_rcf(jg))%t_wml_lk(jc,jb)
            END DO  ! End of loop over grid boxes with lakes
          ENDIF  ! llake

        ENDIF   ! inwp_surface > 0

      ENDDO
!$OMP END DO NOWAIT

    ENDDO
!$OMP END PARALLEL

  END SUBROUTINE copy_initicon2prog_sfc




  !-------------
  !>
  !! SUBROUTINE allocate_initicon
  !! Allocates the components of the initicon data type
  !!
  !! @par Revision History
  !! Initial version by Guenther Zaengl, DWD(2011-07-14)
  !!
  !!
  SUBROUTINE allocate_initicon (p_patch, initicon)

    TYPE(t_patch),          INTENT(IN)    :: p_patch(:)
    TYPE(t_initicon_state), INTENT(INOUT) :: initicon(:)

    ! Local variables: loop control and dimensions
    CHARACTER(MAX_CHAR_LENGTH), PARAMETER :: routine = "mo_nh_initicon::allocate_initicon"
    INTEGER :: jg, nlev, nlevp1, nblks_c, nblks_v, nblks_e

!-------------------------------------------------------------------------

    ! Loop over model domains
    DO jg = 1, n_dom

      nlev = p_patch(jg)%nlev
      nlevp1 = nlev + 1
      nblks_c = p_patch(jg)%nblks_c
      nblks_v = p_patch(jg)%nblks_v
      nblks_e = p_patch(jg)%nblks_e


      ! basic init_icon data
      ALLOCATE(initicon(jg)%topography_c    (nproma,nblks_c),        &
               initicon(jg)%topography_v    (nproma,nblks_v),        &
               initicon(jg)%z_ifc           (nproma,nlevp1,nblks_c), &
               initicon(jg)%z_mc            (nproma,nlev  ,nblks_c) )

      ! Allocate atmospheric output data
      ALLOCATE(initicon(jg)%atm%vn        (nproma,nlev  ,nblks_e), &
               initicon(jg)%atm%u         (nproma,nlev  ,nblks_c), &
               initicon(jg)%atm%v         (nproma,nlev  ,nblks_c), &
               initicon(jg)%atm%w         (nproma,nlevp1,nblks_c), &
               initicon(jg)%atm%temp      (nproma,nlev  ,nblks_c), &
               initicon(jg)%atm%exner     (nproma,nlev  ,nblks_c), &
               initicon(jg)%atm%pres      (nproma,nlev  ,nblks_c), &
               initicon(jg)%atm%rho       (nproma,nlev  ,nblks_c), &
               initicon(jg)%atm%theta_v   (nproma,nlev  ,nblks_c), &
               initicon(jg)%atm%qv        (nproma,nlev  ,nblks_c), &
               initicon(jg)%atm%qc        (nproma,nlev  ,nblks_c), &
               initicon(jg)%atm%qi        (nproma,nlev  ,nblks_c), &
               initicon(jg)%atm%qr        (nproma,nlev  ,nblks_c), &
               initicon(jg)%atm%qs        (nproma,nlev  ,nblks_c)  )

      ! Allocate surface output data
      ALLOCATE(initicon(jg)%sfc%tskin    (nproma,nblks_c             ), &
               initicon(jg)%sfc%sst      (nproma,nblks_c             ), &
               initicon(jg)%sfc%tsnow    (nproma,nblks_c             ), &
               initicon(jg)%sfc%snowalb  (nproma,nblks_c             ), &
               initicon(jg)%sfc%snowweq  (nproma,nblks_c             ), &
               initicon(jg)%sfc%snowdens (nproma,nblks_c             ), &
               initicon(jg)%sfc%skinres  (nproma,nblks_c             ), &
               initicon(jg)%sfc%ls_mask  (nproma,nblks_c             ), &
               initicon(jg)%sfc%seaice   (nproma,nblks_c             ), &
               initicon(jg)%sfc%tsoil    (nproma,nblks_c,0:nlev_soil ), &
               initicon(jg)%sfc%wsoil    (nproma,nblks_c,nlev_soil)     )

      ! allocate groups for list of fields that must be read during initialization
      ALLOCATE(initicon(jg)%sfc%grp_vars_fg (100)        , &
               initicon(jg)%sfc%grp_vars_ana(100)        , &
               initicon(jg)%sfc%grp_vars_fg_default (100), &
               initicon(jg)%sfc%grp_vars_ana_default(100)  )
    ENDDO ! loop over model domains

    ! read the map file into dictionary data structure:
    CALL dict_init(ana_varnames_dict, lcase_sensitive=.FALSE.)
    IF(p_pe == p_io ) THEN 
      IF(ana_varnames_map_file /= ' ') THEN
        CALL dict_loadfile(ana_varnames_dict, TRIM(ana_varnames_map_file))
      END IF
    end IF

  END SUBROUTINE allocate_initicon


  !-------------
  !>
  !! SUBROUTINE allocate_ifs_atm
  !! Allocates fields for IFS read-in
  !!
  SUBROUTINE allocate_ifs_atm (jg, nblks_c, initicon)
    INTEGER,                INTENT(IN)    :: jg, nblks_c
    TYPE(t_initicon_state), INTENT(INOUT) :: initicon(:)
    ! Local variables: loop control and dimensions
    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: routine = 'mo_nh_initicon:allocate_ifs_atm'

    IF (nlev_in == 0) THEN
      CALL finish(routine, "Number of input levels <nlev_in> not yet initialized.")
    END IF

    ! Allocate atmospheric input data
    ALLOCATE(initicon(jg)%atm_in%psfc(nproma,         nblks_c), &
      initicon(jg)%atm_in%phi_sfc(nproma,      nblks_c), &
      initicon(jg)%atm_in%pres (nproma,nlev_in,nblks_c), &
      initicon(jg)%atm_in%z3d  (nproma,nlev_in,nblks_c), &
      initicon(jg)%atm_in%temp (nproma,nlev_in,nblks_c), &
      initicon(jg)%atm_in%u    (nproma,nlev_in,nblks_c), &
      initicon(jg)%atm_in%v    (nproma,nlev_in,nblks_c), &
      initicon(jg)%atm_in%w    (nproma,nlev_in,nblks_c), &
      initicon(jg)%atm_in%omega(nproma,nlev_in,nblks_c), &
      initicon(jg)%atm_in%qv   (nproma,nlev_in,nblks_c), &
      initicon(jg)%atm_in%qc   (nproma,nlev_in,nblks_c), &
      initicon(jg)%atm_in%qi   (nproma,nlev_in,nblks_c), &
      initicon(jg)%atm_in%qr   (nproma,nlev_in,nblks_c), &
      initicon(jg)%atm_in%qs   (nproma,nlev_in,nblks_c)  )
    initicon(jg)%atm_in%linitialized = .TRUE.
  END SUBROUTINE allocate_ifs_atm


  !-------------
  !>
  !! SUBROUTINE allocate_ifs_sfc
  !! Allocates fields for IFS read-in
  !!
  SUBROUTINE allocate_ifs_sfc (jg, nblks_c, initicon)
    INTEGER,                INTENT(IN)    :: jg, nblks_c
    TYPE(t_initicon_state), INTENT(INOUT) :: initicon(:)

    ! Allocate surface input data
    ! The extra soil temperature levels are not read in; they are only used to simplify vertical interpolation
    ALLOCATE(initicon(jg)%sfc_in%phi      (nproma,nblks_c                ), &
      initicon(jg)%sfc_in%tskin    (nproma,nblks_c                ), &
      initicon(jg)%sfc_in%sst      (nproma,nblks_c                ), &
      initicon(jg)%sfc_in%tsnow    (nproma,nblks_c                ), &
      initicon(jg)%sfc_in%snowalb  (nproma,nblks_c                ), &
      initicon(jg)%sfc_in%snowweq  (nproma,nblks_c                ), &
      initicon(jg)%sfc_in%snowdens (nproma,nblks_c                ), &
      initicon(jg)%sfc_in%skinres  (nproma,nblks_c                ), &
      initicon(jg)%sfc_in%ls_mask  (nproma,nblks_c                ), &
      initicon(jg)%sfc_in%seaice   (nproma,nblks_c                ), &
      initicon(jg)%sfc_in%tsoil    (nproma,nblks_c,0:nlevsoil_in+1), &
      initicon(jg)%sfc_in%wsoil    (nproma,nblks_c,0:nlevsoil_in+1)  )
    initicon(jg)%sfc_in%linitialized = .TRUE.
  END SUBROUTINE allocate_ifs_sfc


  !-------------
  !>
  !! SUBROUTINE deallocate_initicon
  !! Deallocates the components of the initicon data type
  !!
  !! @par Revision History
  !! Initial version by Guenther Zaengl, DWD(2011-07-14)
  !!
  !!
  SUBROUTINE deallocate_initicon (initicon)

    TYPE(t_initicon_state), INTENT(INOUT) :: initicon(:)

    ! Local variables: loop control
    INTEGER :: jg

!-------------------------------------------------------------------------

    ! Loop over model domains
    DO jg = 1, n_dom

      ! basic init_icon data
      DEALLOCATE(initicon(jg)%topography_c,     &
                 initicon(jg)%topography_v,     &
                 initicon(jg)%z_ifc,            &
                 initicon(jg)%z_mc              )

      ! atmospheric output data
      DEALLOCATE(initicon(jg)%atm%vn,      &
                 initicon(jg)%atm%u,       &
                 initicon(jg)%atm%v,       &
                 initicon(jg)%atm%w,       &
                 initicon(jg)%atm%temp,    &
                 initicon(jg)%atm%exner,   &
                 initicon(jg)%atm%pres,    &  
                 initicon(jg)%atm%rho,     &
                 initicon(jg)%atm%theta_v, &
                 initicon(jg)%atm%qv,      &
                 initicon(jg)%atm%qc,      &
                 initicon(jg)%atm%qi,      &
                 initicon(jg)%atm%qr,      &
                 initicon(jg)%atm%qs       )


      ! surface output data
      DEALLOCATE(initicon(jg)%sfc%tskin,    &
                 initicon(jg)%sfc%sst,      &
                 initicon(jg)%sfc%tsnow,    &
                 initicon(jg)%sfc%snowalb,  &
                 initicon(jg)%sfc%snowweq,  &
                 initicon(jg)%sfc%snowdens, &
                 initicon(jg)%sfc%skinres,  &
                 initicon(jg)%sfc%ls_mask,  &
                 initicon(jg)%sfc%seaice,   &
                 initicon(jg)%sfc%tsoil,    &
                 initicon(jg)%sfc%wsoil     )

      ! deallocate groups for list of fields that must be read during initialization
      DEALLOCATE(initicon(jg)%sfc%grp_vars_fg  , &
                 initicon(jg)%sfc%grp_vars_ana   )
    ENDDO ! loop over model domains

    ! destroy variable name dictionaries:
    CALL dict_finalize(ana_varnames_dict)

  END SUBROUTINE deallocate_initicon


  !-------------
  !>
  !! SUBROUTINE deallocate_ifs_atm
  !! Deallocates the components of the initicon data type
  !!
  SUBROUTINE deallocate_ifs_atm (initicon)
    TYPE(t_initicon_state), INTENT(INOUT) :: initicon(:)
    ! Local variables: loop control
    INTEGER :: jg

    ! Loop over model domains
    DO jg = 1, n_dom
      IF (.NOT. initicon(jg)%atm_in%linitialized) CYCLE

      ! atmospheric input data
      DEALLOCATE(initicon(jg)%atm_in%psfc,    &
                 initicon(jg)%atm_in%phi_sfc, &
                 initicon(jg)%atm_in%pres,    &
                 initicon(jg)%atm_in%z3d,     &
                 initicon(jg)%atm_in%temp,    &
                 initicon(jg)%atm_in%u,       &
                 initicon(jg)%atm_in%v,       &
                 initicon(jg)%atm_in%w,       &
                 initicon(jg)%atm_in%omega,   &
                 initicon(jg)%atm_in%qv,      &
                 initicon(jg)%atm_in%qc,      &
                 initicon(jg)%atm_in%qi,      &
                 initicon(jg)%atm_in%qr,      &
                 initicon(jg)%atm_in%qs )
      IF (ALLOCATED(initicon(jg)%atm_in%vn)) THEN
        DEALLOCATE(initicon(jg)%atm_in%vn)
      ENDIF
      initicon(jg)%atm_in%linitialized = .FALSE.
    ENDDO ! loop over model domains
  END SUBROUTINE deallocate_ifs_atm


  !-------------
  !>
  !! SUBROUTINE deallocate_ifs_sfc
  !! Deallocates IFS fields
  !!
  SUBROUTINE deallocate_ifs_sfc (initicon)
    TYPE(t_initicon_state), INTENT(INOUT) :: initicon(:)
    ! Local variables: loop control
    INTEGER :: jg

    ! Loop over model domains
    DO jg = 1, n_dom
      IF (.NOT. initicon(jg)%sfc_in%linitialized) CYCLE
      ! surface input data
      DEALLOCATE(initicon(jg)%sfc_in%phi,      &
                 initicon(jg)%sfc_in%tskin,    &
                 initicon(jg)%sfc_in%sst,    &
                 initicon(jg)%sfc_in%tsnow,    &
                 initicon(jg)%sfc_in%snowalb,  &
                 initicon(jg)%sfc_in%snowweq,  &
                 initicon(jg)%sfc_in%snowdens, &
                 initicon(jg)%sfc_in%skinres,  &
                 initicon(jg)%sfc_in%ls_mask,  &
                 initicon(jg)%sfc_in%seaice,   &
                 initicon(jg)%sfc_in%tsoil,    &
                 initicon(jg)%sfc_in%wsoil     )
      initicon(jg)%sfc_in%linitialized = .FALSE.
    ENDDO ! loop over model domains
  END SUBROUTINE deallocate_ifs_sfc

END MODULE mo_nh_initicon

