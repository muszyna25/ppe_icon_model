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

MODULE mo_initicon

  USE mo_kind,                ONLY: wp, vp
  USE mo_parallel_config,     ONLY: nproma
  USE mo_run_config,          ONLY: iqv, iqc, iqi, iqr, iqs, iqm_max, iforcing, check_uuid_gracefully
  USE mo_dynamics_config,     ONLY: nnow, nnow_rcf
  USE mo_model_domain,        ONLY: t_patch
  USE mo_nonhydro_types,      ONLY: t_nh_state, t_nh_prog, t_nh_diag, t_nh_metrics
  USE mo_nwp_phy_types,       ONLY: t_nwp_phy_diag
  USE mo_nwp_lnd_types,       ONLY: t_lnd_state, t_lnd_prog, t_lnd_diag
  USE mo_intp_data_strc,      ONLY: t_int_state
  USE mo_ext_data_types,      ONLY: t_external_data
  USE mo_grf_intp_data_strc,  ONLY: t_gridref_state
  USE mo_initicon_types,      ONLY: t_initicon_state, ana_varnames_dict
  USE mo_initicon_config,     ONLY: init_mode, dt_iau, nlevatm_in, lvert_remap_fg, &
    &                               rho_incr_filter_wgt, lread_ana, ltile_init, &
    &                               lp2cintp_incr, lp2cintp_sfcana, ltile_coldstart, lconsistency_checks, &
    &                               niter_divdamp, niter_diffu, lanaread_tseasfc
  USE mo_nwp_tuning_config,   ONLY: max_freshsnow_inc
  USE mo_impl_constants,      ONLY: SUCCESS, MAX_CHAR_LENGTH, MODE_DWDANA,   &
    &                               MODE_IAU, MODE_IAU_OLD, MODE_IFSANA,              &
    &                               MODE_ICONVREMAP, MODE_COMBINED, MODE_COSMODE,     &
    &                               min_rlcell, INWP, min_rledge_int, grf_bdywidth_c, &
    &                               min_rlcell_int, dzsoil_icon => dzsoil
  USE mo_physical_constants,  ONLY: rd, cpd, cvd, p0ref, vtmpc1, grav, rd_o_cpd, tmelt, tf_salt
  USE mo_exception,           ONLY: message, finish
  USE mo_grid_config,         ONLY: n_dom, l_limited_area
  USE mo_nh_init_utils,       ONLY: convert_thdvars, init_w
  USE mo_nh_init_nest_utils,  ONLY: interpolate_vn_increments
  USE mo_util_phys,           ONLY: virtual_temp
  USE mo_satad,               ONLY: sat_pres_ice, spec_humi
  USE mo_lnd_nwp_config,      ONLY: nlev_soil, ntiles_total, ntiles_lnd, llake, &
    &                               isub_lake, isub_water, lsnowtile, frlnd_thrhld, &
    &                               frlake_thrhld, lprog_albsi
  USE mo_seaice_nwp,          ONLY: frsi_min
  USE mo_atm_phy_nwp_config,  ONLY: iprog_aero
  USE mo_phyparam_soil,       ONLY: cporv, cadp, crhosmaxf, crhosmin_ml, crhosmax_ml
  USE mo_nwp_soil_init,       ONLY: get_wsnow
  USE mo_nh_vert_interp,      ONLY: vert_interp_atm, vert_interp_sfc
  USE mo_intp_rbf,            ONLY: rbf_vec_interpol_cell
  USE mo_nh_diagnose_pres_temp,ONLY: diagnose_pres_temp
  USE mo_loopindices,         ONLY: get_indices_c, get_indices_e
  USE mo_sync,                ONLY: sync_patch_array, SYNC_E, SYNC_C
  USE mo_math_laplace,        ONLY: nabla2_vec, nabla4_vec
  USE mo_cdi,                 ONLY: cdiDefAdditionalKey, cdiInqMissval, FILETYPE_NC2, FILETYPE_NC4, FILETYPE_GRB2
  USE mo_flake,               ONLY: flake_coldinit
  USE mo_initicon_utils,      ONLY: fill_tile_points, init_snowtiles,             &
                                  & copy_initicon2prog_atm, copy_initicon2prog_sfc, construct_initicon, &
                                  & deallocate_initicon, deallocate_extana_atm, deallocate_extana_sfc, &
                                  & copy_fg2initicon, initVarnamesDict, printChecksums, init_aerosol
  USE mo_initicon_io,         ONLY: read_extana_atm, read_extana_sfc, fetch_dwdfg_atm, fetch_dwdana_sfc, &
                                  & process_input_dwdana_sfc, process_input_dwdana_atm, process_input_dwdfg_sfc, &
                                  & fetch_dwdfg_sfc, fetch_dwdfg_atm_ii, fetch_dwdana_atm, &
                                  & fgFilename, fgFiletype, anaFilename, anaFiletype
  USE mo_input_request_list,  ONLY: t_InputRequestList, InputRequestList_create
  USE mo_mpi,                 ONLY: my_process_is_stdio
  USE mo_input_instructions,  ONLY: t_readInstructionListPtr, readInstructionList_make, kInputSourceAna, &
                                    kInputSourceBoth, kInputSourceCold
  USE mo_util_uuid_types,     ONLY: t_uuid
  USE mo_nwp_sfc_utils,       ONLY: seaice_albedo_coldstart

  IMPLICIT NONE


  PRIVATE

  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_initicon'

  TYPE(t_initicon_state), ALLOCATABLE, TARGET :: initicon(:)

  PUBLIC :: init_icon

  DOUBLE PRECISION, PARAMETER :: kMissval = -9.E+15


  CONTAINS

  !-------------
  !>
  !! SUBROUTINE init_icon
  !! ICON initialization routine: Reads in either DWD or external (IFS/COSMO) analysis
  !!
  !! @par Revision History
  !! Initial version by Guenther Zaengl, DWD(2011-07-14)
  !!
  !!
  SUBROUTINE init_icon (p_patch,  p_int_state, p_grf_state, p_nh_state, &
    &                   ext_data, prm_diag, p_lnd_state)

    TYPE(t_patch),          INTENT(INOUT)              :: p_patch(:)
    TYPE(t_int_state),      INTENT(IN)              :: p_int_state(:)
    TYPE(t_gridref_state),  INTENT(IN)              :: p_grf_state(:)
    TYPE(t_nh_state),       INTENT(INOUT)           :: p_nh_state(:)
    TYPE(t_nwp_phy_diag),   INTENT(INOUT), OPTIONAL :: prm_diag(:)
    TYPE(t_lnd_state),      INTENT(INOUT), OPTIONAL :: p_lnd_state(:)
    TYPE(t_external_data),  INTENT(INOUT), OPTIONAL :: ext_data(:)

    CHARACTER(LEN = *), PARAMETER :: routine = modname//':init_icon'
    INTEGER :: jg, ist
    TYPE(t_readInstructionListPtr) :: inputInstructions(n_dom)

    ! Allocate initicon data type
    ALLOCATE (initicon(n_dom),                         &
      &       stat=ist)
    IF (ist /= SUCCESS)  CALL finish(TRIM(routine),'allocation for initicon failed')

    DO jg = 1, n_dom
      CALL construct_initicon(initicon(jg), p_patch(jg), ext_data(jg)%atm%topography_c, p_nh_state(jg)%metrics)
    END DO

    ! Read IN the dictionary for the variable names (IF we need it)
    CALL initVarnamesDict(ana_varnames_dict)


    ! -----------------------------------------------
    ! make the CDI aware of some custom GRIB keys
    ! -----------------------------------------------

    CALL cdiDefAdditionalKey("localInformationNumber")
    CALL cdiDefAdditionalKey("localNumberOfExperiment")
    CALL cdiDefAdditionalKey("typeOfFirstFixedSurface")
    CALL cdiDefAdditionalKey("typeOfGeneratingProcess")
    CALL cdiDefAdditionalKey("backgroundProcess")
    CALL cdiDefAdditionalKey("totalNumberOfTileAttributePairs")
    CALL cdiDefAdditionalKey("tileIndex")
    CALL cdiDefAdditionalKey("tileAttribute")



    ! -----------------------------------------------
    ! generate analysis/FG input instructions
    ! -----------------------------------------------
    DO jg = 1, n_dom
      inputInstructions(jg)%ptr => NULL()
      IF(p_patch(jg)%ldom_active) inputInstructions(jg)%ptr => readInstructionList_make(p_patch(jg), init_mode)
    END DO

    ! -----------------------------------------------
    ! READ AND process the input DATA
    ! -----------------------------------------------
    CALL print_init_mode()

    ! read and initialize ICON prognostic fields
    !
    CALL process_input_data(p_patch, inputInstructions, p_nh_state, p_int_state, p_grf_state, ext_data, prm_diag, p_lnd_state)
    CALL printChecksums(initicon, p_nh_state, p_lnd_state)

    ! Deallocate initicon data type
    !
    CALL deallocate_initicon(initicon)
    CALL deallocate_extana_atm (initicon)
    CALL deallocate_extana_sfc (initicon)

    DEALLOCATE (initicon, stat=ist)
    IF (ist /= success) CALL finish(TRIM(routine),'deallocation for initicon failed')
    DO jg = 1, n_dom
      IF(p_patch(jg)%ldom_active) THEN
        IF(my_process_is_stdio()) CALL inputInstructions(jg)%ptr%printSummary(jg)
        CALL inputInstructions(jg)%ptr%destruct()
        DEALLOCATE(inputInstructions(jg)%ptr, stat=ist)
        IF(ist /= success) CALL finish(TRIM(routine),'deallocation of an input instruction list failed')
      END IF
    END DO

    ! splitting of sea-points list into open water and sea-ice points could be placed
    ! here, instead of nwp_phy_init/init_nwp_phy
    ! however, one needs to make sure that it is called for both restart and non-restart
    ! runs. Could not be included into mo_ext_data_state/init_index_lists due to its
    ! dependence on p_diag_lnd.
!DR    CALL init_sea_lists(p_patch, ext_data, p_diag_lnd, lseaice)

  END SUBROUTINE init_icon

  ! Write an output line that informs the user of the init_mode we are using (failing the program if init_mode is invalid).
  SUBROUTINE print_init_mode()
    SELECT CASE(init_mode)
        CASE(MODE_DWDANA)   
            CALL message(modname,'MODE_DWD: perform initialization with DWD analysis')
        CASE(MODE_ICONVREMAP)   
            CALL message(modname,'MODE_VREMAP: read ICON data and perform vertical remapping')
        CASE (MODE_IAU_OLD)
            CALL message(modname,'MODE_IAU_OLD: perform initialization with incremental analysis update &
                                 &(retained for backward compatibility)')
        CASE (MODE_IAU)
            CALL message(modname,'MODE_IAU: perform initialization with incremental analysis update, including snow increments')
        CASE(MODE_IFSANA)   
            CALL message(modname,'MODE_IFS: perform initialization with IFS analysis')
        CASE(MODE_COMBINED)
            CALL message(modname,'MODE_COMBINED: IFS-atm + GME-soil')
        CASE(MODE_COSMODE)
            CALL message(modname,'MODE_COSMODE: COSMO-atm + COSMO-soil')
        CASE DEFAULT
            CALL finish(modname, "Invalid operation mode!")
    END SELECT
  END SUBROUTINE print_init_mode

  FUNCTION gridUuids(p_patch) RESULT(resultVar)
    TYPE(t_patch), INTENT(IN) :: p_patch(:)
    TYPE(t_uuid), DIMENSION(n_dom) :: resultVar

    resultVar(:) = p_patch(:)%grid_uuid
  END FUNCTION gridUuids

  ! Read the data from the first-guess file.
  SUBROUTINE read_dwdfg(p_patch, inputInstructions, p_nh_state, prm_diag, p_lnd_state)
    TYPE(t_patch), INTENT(INOUT) :: p_patch(:)
    TYPE(t_readInstructionListPtr) :: inputInstructions(n_dom)
    TYPE(t_nh_state), INTENT(INOUT) :: p_nh_state(:)
    TYPE(t_nwp_phy_diag), INTENT(INOUT), OPTIONAL :: prm_diag(:)
    TYPE(t_lnd_state), INTENT(INOUT), OPTIONAL :: p_lnd_state(:)

    CHARACTER(LEN = *), PARAMETER :: routine = modname//":read_dwdfg"
    INTEGER :: jg
    CLASS(t_InputRequestList), POINTER :: requestList

    !The input file paths & types are NOT initialized IN all modes, so we need to avoid creating InputRequestLists IN these cases.
    SELECT CASE(init_mode)
        CASE(MODE_IFSANA)    !MODE_IFSANA uses the read_extana_*() routines, which directly use NetCDF input.
            RETURN
        CASE(MODE_DWDANA, MODE_IAU_OLD, MODE_IAU, MODE_ICONVREMAP, MODE_COMBINED, MODE_COSMODE)
        CASE DEFAULT
            CALL finish(routine, "assertion failed: unknown init_mode")
    END SELECT

    ! Create a request list for all the relevant variable names.
    requestList => InputRequestList_create()
    DO jg = 1, n_dom
      IF(p_patch(jg)%ldom_active) THEN
        CALL inputInstructions(jg)%ptr%fileRequests(requestList, lIsFg = .TRUE.)
      ENDIF
    END DO

    ! Scan the input files AND distribute the relevant variables across the processes.
    DO jg = 1, n_dom
        IF(p_patch(jg)%ldom_active) THEN
            IF(my_process_is_stdio()) THEN
                CALL message (TRIM(routine), 'read atm_FG fields from '//TRIM(fgFilename(p_patch(jg))))
            ENDIF  ! p_io
            SELECT CASE(fgFiletype())
                CASE(FILETYPE_NC2, FILETYPE_NC4)
                    CALL requestList%readFile(p_patch(jg), TRIM(fgFilename(p_patch(jg))), .TRUE.)
                CASE(FILETYPE_GRB2)
                    CALL requestList%readFile(p_patch(jg), TRIM(fgFilename(p_patch(jg))), .TRUE., opt_dict = ana_varnames_dict)
                CASE DEFAULT
                    CALL finish(routine, "Unknown file TYPE")
            END SELECT
        END IF
    END DO
    IF(my_process_is_stdio()) THEN
        CALL requestList%printInventory()
        IF(lconsistency_checks) THEN
          CALL requestList%checkRuntypeAndUuids([CHARACTER(LEN=1)::], gridUuids(p_patch), lIsFg=.TRUE., &
            lHardCheckUuids=.NOT.check_uuid_gracefully)
        END IF
    END IF

    ! Fetch the input DATA from the request list.
    SELECT CASE(init_mode)
        CASE(MODE_DWDANA, MODE_IAU_OLD, MODE_IAU)
            CALL fetch_dwdfg_atm(requestList, p_patch, p_nh_state, initicon, inputInstructions)
            CALL fetch_dwdfg_sfc(requestList, p_patch, prm_diag, p_lnd_state, inputInstructions)
        CASE(MODE_ICONVREMAP)
            CALL fetch_dwdfg_atm_ii(requestList, p_patch, initicon, inputInstructions)
            CALL fetch_dwdfg_sfc(requestList, p_patch, prm_diag, p_lnd_state, inputInstructions)
        CASE(MODE_COMBINED, MODE_COSMODE)
            CALL fetch_dwdfg_sfc(requestList, p_patch, prm_diag, p_lnd_state, inputInstructions)
    END SELECT

    ! Cleanup.
    CALL requestList%destruct()
    DEALLOCATE(requestList)
  END SUBROUTINE read_dwdfg

  ! Do postprocessing of data from first-guess file.
  SUBROUTINE process_dwdfg(p_patch, p_nh_state, p_int_state, p_grf_state, ext_data, p_lnd_state, prm_diag)
    TYPE(t_patch), INTENT(IN) :: p_patch(:)
    TYPE(t_nh_state), INTENT(INOUT) :: p_nh_state(:)
    TYPE(t_int_state), INTENT(IN) :: p_int_state(:)
    TYPE(t_gridref_state), INTENT(IN) :: p_grf_state(:)
    TYPE(t_external_data), INTENT(INOUT) :: ext_data(:)
    TYPE(t_lnd_state), INTENT(INOUT), OPTIONAL :: p_lnd_state(:)
    TYPE(t_nwp_phy_diag), INTENT(INOUT), OPTIONAL :: prm_diag(:)

    CHARACTER(LEN = *), PARAMETER :: routine = modname//":process_dwdfg"

    SELECT CASE(init_mode)
        CASE(MODE_ICONVREMAP)
            CALL process_input_dwdfg_sfc (p_patch, p_lnd_state, ext_data)
        CASE(MODE_DWDANA, MODE_IAU_OLD, MODE_IAU, MODE_COMBINED, MODE_COSMODE)
            IF (lvert_remap_fg) THEN ! apply vertical remapping of FG input (requires that the number of model levels
                                     ! does not change; otherwise, init_mode = 7 must be used based on a full analysis)
                CALL copy_fg2initicon(p_patch, initicon, p_nh_state)
                CALL vert_interp_atm(p_patch, p_nh_state, p_int_state, p_grf_state, p_patch(:)%nlev, initicon, &
                &                    opt_convert_omega2w=.FALSE.)
                CALL copy_initicon2prog_atm(p_patch, initicon, p_nh_state)
            END IF
            CALL process_input_dwdfg_sfc (p_patch, p_lnd_state, ext_data)
            IF(ANY((/MODE_IAU_OLD, MODE_IAU/) == init_mode)) THEN
                ! In case of tile coldstart, fill sub-grid scale land
                ! and water points with reasonable data from
                ! neighboring grid points where possible; In case of
                ! snowtile warmstart, the index lists for snow-covered
                ! / snow-free points need to be initialized
                IF (ntiles_total > 1 .AND. ltile_init) THEN
                    CALL fill_tile_points(p_patch, p_lnd_state, ext_data, process_ana_vars=.FALSE.)
                ELSE IF (ntiles_total > 1 .AND. lsnowtile .AND. .NOT. ltile_coldstart) THEN
                    CALL init_snowtiles(p_patch, p_lnd_state, ext_data)
                END IF
            END IF
    END SELECT
    ! Init aerosol field from climatology if no first-guess data have been available
    IF (iprog_aero == 1) CALL init_aerosol(p_patch, ext_data, prm_diag)
  END SUBROUTINE process_dwdfg

  ! Read data from analysis files.
  SUBROUTINE read_dwdana(p_patch, inputInstructions, p_nh_state, p_lnd_state)
    TYPE(t_patch), INTENT(IN) :: p_patch(:)
    TYPE(t_readInstructionListPtr) :: inputInstructions(n_dom)
    TYPE(t_nh_state), INTENT(INOUT) :: p_nh_state(:)
    TYPE(t_lnd_state), INTENT(INOUT), OPTIONAL :: p_lnd_state(:)

    CHARACTER(LEN = *), PARAMETER :: routine = modname//":read_dwdana"
#ifndef __GFORTRAN__
    CHARACTER(LEN = :), ALLOCATABLE :: incrementsList(:)
#else
    CHARACTER(LEN = 9) :: incrementsList_IAU(8)
    CHARACTER(LEN = 4) :: incrementsList_IAU_OLD(6)
    CHARACTER(LEN = 1) :: incrementsList_DEFAULT(1)
#endif
    CLASS(t_InputRequestList), POINTER :: requestList
    INTEGER :: jg

    !The input file paths & types are NOT initialized IN all modes, so we need to avoid creating InputRequestLists IN these cases.
    SELECT CASE(init_mode)
        CASE(MODE_IFSANA)
            CALL read_extana_atm(p_patch, initicon)
            ! Perform vertical interpolation from intermediate
            ! IFS2ICON grid to ICON grid and convert variables to the
            ! NH set of prognostic variables
            IF (iforcing == inwp) CALL read_extana_sfc(p_patch, initicon)
            RETURN
        CASE(MODE_DWDANA, MODE_IAU_OLD, MODE_IAU, MODE_ICONVREMAP, MODE_COMBINED, MODE_COSMODE)
        CASE DEFAULT
            CALL finish(routine, "assertion failed: unknown init_mode")
    END SELECT

    ! Create a request list for all the relevant variable names.
    requestList => InputRequestList_create()
    SELECT CASE(init_mode)
        CASE(MODE_COMBINED, MODE_COSMODE)
            CALL read_extana_atm(p_patch, initicon)
    END SELECT
    DO jg = 1, n_dom
      IF(p_patch(jg)%ldom_active) THEN
        CALL inputInstructions(jg)%ptr%fileRequests(requestList, lIsFg = .FALSE.)
      ENDIF
    END DO

    ! Scan the input files AND distribute the relevant variables across the processes.
    DO jg = 1, n_dom
        IF(p_patch(jg)%ldom_active .AND. lread_ana) THEN
            IF (lp2cintp_incr(jg) .AND. lp2cintp_sfcana(jg)) CYCLE
            IF(my_process_is_stdio()) THEN
                CALL message (TRIM(routine), 'read atm_ANA fields from '//TRIM(anaFilename(p_patch(jg))))
            ENDIF  ! p_io
            SELECT CASE(anaFiletype())
                CASE(FILETYPE_NC2, FILETYPE_NC4)
                    CALL requestList%readFile(p_patch(jg), TRIM(anaFilename(p_patch(jg))), .FALSE.)
                CASE(FILETYPE_GRB2)
                    CALL requestList%readFile(p_patch(jg), TRIM(anaFilename(p_patch(jg))), .FALSE., opt_dict = ana_varnames_dict)
                CASE DEFAULT
                    CALL finish(routine, "Unknown file TYPE")
            END SELECT
        END IF
    END DO
    IF(my_process_is_stdio()) THEN
        CALL requestList%printInventory()
        IF(lconsistency_checks) THEN
! Workaround for GNU compiler (<6.0), which still does not fully support deferred length character arrays
#ifndef __GFORTRAN__
            SELECT CASE(init_mode)
                CASE(MODE_IAU)
                    incrementsList = [CHARACTER(LEN=9) :: 'u', 'v', 'pres', 'temp', 'qv', 'w_so', 'h_snow', 'freshsnow']
                CASE(MODE_IAU_OLD)
                    incrementsList = [CHARACTER(LEN=4) :: 'u', 'v', 'pres', 'temp', 'qv', 'w_so']
                CASE DEFAULT
                    incrementsList = [CHARACTER(LEN=1) :: ]
            END SELECT
            CALL requestList%checkRuntypeAndUuids(incrementsList, gridUuids(p_patch), lIsFg = .FALSE., &
              lHardCheckUuids = .NOT.check_uuid_gracefully)
#else
            SELECT CASE(init_mode)
                CASE(MODE_IAU)
                    incrementsList_IAU = (/'u        ', 'v        ', 'pres     ', 'temp     ', &
                      &                    'qv       ', 'w_so     ', 'h_snow   ', 'freshsnow'/)
                    CALL requestList%checkRuntypeAndUuids(incrementsList_IAU, gridUuids(p_patch), lIsFg = .FALSE., &
                      lHardCheckUuids = .NOT.check_uuid_gracefully)
            write(0,*) "incrementsList_IAU: ", incrementsList_IAU
                CASE(MODE_IAU_OLD)
                    incrementsList_IAU_OLD = (/'u   ', 'v   ', 'pres', 'temp', 'qv  ', 'w_so'/)
                    CALL requestList%checkRuntypeAndUuids(incrementsList_IAU_OLD, gridUuids(p_patch), lIsFg = .FALSE., &
                      lHardCheckUuids = .NOT.check_uuid_gracefully)
            write(0,*) "incrementsList_IAU_OLD: ", incrementsList_IAU_OLD
                CASE DEFAULT
                    incrementsList_DEFAULT = (/' '/)
                    CALL requestList%checkRuntypeAndUuids(incrementsList_DEFAULT, gridUuids(p_patch), lIsFg = .FALSE., &
                      lHardCheckUuids = .NOT.check_uuid_gracefully)
            END SELECT
#endif
        END IF
    END IF

    ! Fetch the input DATA from the request list.
    SELECT CASE(init_mode)
        CASE(MODE_DWDANA, MODE_IAU_OLD, MODE_IAU)
            IF(lread_ana) CALL fetch_dwdana_atm(requestList, p_patch, p_nh_state, initicon, inputInstructions)
            IF(lread_ana) CALL fetch_dwdana_sfc(requestList, p_patch, p_lnd_state, initicon, inputInstructions)
        CASE(MODE_COMBINED, MODE_COSMODE)
            IF(lread_ana) CALL fetch_dwdana_sfc(requestList, p_patch, p_lnd_state, initicon, inputInstructions)
        CASE(MODE_ICONVREMAP)
            IF(lread_ana) CALL fetch_dwdana_sfc(requestList, p_patch, p_lnd_state, initicon, inputInstructions)
    END SELECT

    ! Cleanup.
    CALL requestList%destruct()
    DEALLOCATE(requestList)
  END SUBROUTINE read_dwdana

  ! Do postprocessing of data from analysis files.
  SUBROUTINE process_dwdana(p_patch, inputInstructions, p_nh_state, p_int_state, p_grf_state, ext_data, p_lnd_state)
    TYPE(t_patch), INTENT(IN) :: p_patch(:)
    TYPE(t_readInstructionListPtr) :: inputInstructions(n_dom)
    TYPE(t_nh_state), INTENT(INOUT) :: p_nh_state(:)
    TYPE(t_int_state), INTENT(IN) :: p_int_state(:)
    TYPE(t_gridref_state), INTENT(IN) :: p_grf_state(:)
    TYPE(t_external_data), INTENT(INOUT) :: ext_data(:)
    TYPE(t_lnd_state), INTENT(INOUT), OPTIONAL :: p_lnd_state(:)

    CHARACTER(LEN = *), PARAMETER :: routine = modname//":process_dwdana"

    INTEGER :: jg, jb
    INTEGER :: i_rlstart, i_rlend, i_nchdom
    INTEGER :: i_startblk, i_endblk

    SELECT CASE(init_mode)
        CASE(MODE_DWDANA)
            ! process DWD atmosphere analysis data
            IF(lread_ana) CALL process_input_dwdana_atm(p_patch, initicon)
            ! merge first guess with DA analysis and 
            ! convert variables to the NH set of prognostic variables
            CALL create_dwdana_atm(p_patch, p_nh_state, p_int_state)
        CASE(MODE_IAU_OLD, MODE_IAU)
            ! process DWD atmosphere analysis increments
            IF(lread_ana) CALL process_input_dwdana_atm(p_patch, initicon)
            ! Compute DA increments in terms of the NH set of
            ! prognostic variables
            CALL create_dwdanainc_atm(p_patch, p_nh_state, p_int_state)
        CASE(MODE_COMBINED, MODE_IFSANA)
            ! process IFS atmosphere analysis data
            CALL vert_interp_atm(p_patch, p_nh_state, p_int_state, p_grf_state, nlevatm_in, initicon, &
            &                    opt_convert_omega2w = .TRUE.)
            ! Finally copy the results to the prognostic model
            ! variables
            CALL copy_initicon2prog_atm(p_patch, initicon, p_nh_state)
        CASE(MODE_ICONVREMAP, MODE_COSMODE)
            ! process ICON (DWD) atmosphere first-guess data (having
            ! different vertical levels than the current grid)
            CALL vert_interp_atm(p_patch, p_nh_state, p_int_state, p_grf_state, nlevatm_in, initicon, &
            &                    opt_convert_omega2w = .FALSE.)
            ! Finally copy the results to the prognostic model
            ! variables
            CALL copy_initicon2prog_atm(p_patch, initicon, p_nh_state)
    END SELECT

    SELECT CASE(init_mode)
        CASE(MODE_DWDANA, MODE_ICONVREMAP, MODE_IAU_OLD, MODE_IAU, MODE_COMBINED, MODE_COSMODE)
            ! process DWD land/surface analysis data / increments
            IF(lread_ana) CALL process_input_dwdana_sfc(p_patch, p_lnd_state, initicon)
            ! Add increments to time-shifted first guess in one go.
            ! The following CALL must not be moved after create_dwdana_sfc()!
            IF(ANY((/MODE_IAU_OLD, MODE_IAU/) == init_mode)) THEN
                CALL create_iau_sfc (p_patch, p_nh_state, p_lnd_state, ext_data)
            END IF
            ! get SST from first soil level t_so or t_seasfc
            ! perform consistency checks
            CALL create_dwdana_sfc(p_patch, p_lnd_state, ext_data, inputInstructions)
            IF (ANY((/MODE_IAU_OLD, MODE_IAU/) == init_mode) .AND. ntiles_total > 1) THEN
                ! Call neighbor-filling routine for a second time in
                ! order to ensure that fr_seaice is filled with
                ! meaningful data near coastlines if this field is
                ! read from the analysis
                CALL fill_tile_points(p_patch, p_lnd_state, ext_data, process_ana_vars=.TRUE.)
            END IF
        CASE(MODE_IFSANA)
            IF (iforcing == inwp) THEN
                ! Perform vertical interpolation from intermediate
                ! IFS2ICON grid to ICON grid and convert variables to
                ! the NH set of prognostic variables
                CALL vert_interp_sfc(p_patch, initicon)
                ! Finally copy the results to the prognostic model variables
                CALL copy_initicon2prog_sfc(p_patch, initicon, p_lnd_state, ext_data)
            END IF
    END SELECT

    SELECT CASE(init_mode)
        CASE(MODE_COMBINED,MODE_COSMODE)
            ! Cold-start initialization of the fresh-water lake model
            ! FLake. The procedure is the same as in "int2lm". Note
            ! that no lake ice is assumed at the cold start.
            IF (llake) THEN
                DO jg = 1, n_dom
                    IF (.NOT. p_patch(jg)%ldom_active) CYCLE
                    i_rlstart  = 1
                    i_rlend    = min_rlcell
                    i_nchdom   =  MAX(1,p_patch(jg)%n_childdom)
                    i_startblk = p_patch(jg)%cells%start_blk(i_rlstart,1)
                    i_endblk   = p_patch(jg)%cells%end_blk(i_rlend,i_nchdom)
!$OMP PARALLEL
!$OMP DO PRIVATE(jb)
                    DO jb = i_startblk, i_endblk
                        CALL flake_coldinit(                                        &
                            &   nflkgb      = ext_data(jg)%atm%fp_count    (jb),    &
                            &   idx_lst_fp  = ext_data(jg)%atm%idx_lst_fp(:,jb),    &
                            &   depth_lk    = ext_data(jg)%atm%depth_lk  (:,jb),    &
                            &   tskin       = p_lnd_state(jg)%prog_lnd(nnow_rcf(jg))%t_so_t(:,1,jb,1),&
                            &   t_snow_lk_p = p_lnd_state(jg)%prog_wtr(nnow_rcf(jg))%t_snow_lk(:,jb), &
                            &   h_snow_lk_p = p_lnd_state(jg)%prog_wtr(nnow_rcf(jg))%h_snow_lk(:,jb), &
                            &   t_ice_p     = p_lnd_state(jg)%prog_wtr(nnow_rcf(jg))%t_ice    (:,jb), &
                            &   h_ice_p     = p_lnd_state(jg)%prog_wtr(nnow_rcf(jg))%h_ice    (:,jb), &
                            &   t_mnw_lk_p  = p_lnd_state(jg)%prog_wtr(nnow_rcf(jg))%t_mnw_lk (:,jb), &
                            &   t_wml_lk_p  = p_lnd_state(jg)%prog_wtr(nnow_rcf(jg))%t_wml_lk (:,jb), &
                            &   t_bot_lk_p  = p_lnd_state(jg)%prog_wtr(nnow_rcf(jg))%t_bot_lk (:,jb), &
                            &   c_t_lk_p    = p_lnd_state(jg)%prog_wtr(nnow_rcf(jg))%c_t_lk   (:,jb), &
                            &   h_ml_lk_p   = p_lnd_state(jg)%prog_wtr(nnow_rcf(jg))%h_ml_lk  (:,jb), &
                            &   t_b1_lk_p   = p_lnd_state(jg)%prog_wtr(nnow_rcf(jg))%t_b1_lk  (:,jb), &
                            &   h_b1_lk_p   = p_lnd_state(jg)%prog_wtr(nnow_rcf(jg))%h_b1_lk  (:,jb), &
                            &   t_g_lk_p    = p_lnd_state(jg)%prog_lnd(nnow_rcf(jg))%t_g_t    (:,jb,isub_lake) )
                    ENDDO
!$OMP END DO
!$OMP END PARALLEL
                ENDDO
            ENDIF
    END SELECT

    !
    ! coldstart for prognostic sea-ice albedo in case that alb_si was 
    ! not found in the FG 
    !
    DO jg = 1, n_dom
      IF (.NOT. p_patch(jg)%ldom_active) CYCLE

      IF ( lprog_albsi .AND. inputInstructions(jg)%ptr%sourceOfVar('alb_si') == kInputSourceCold) THEN

        CALL seaice_albedo_coldstart(p_patch(jg), p_lnd_state(jg), ext_data(jg))

      ENDIF
    ENDDO  !jg

  END SUBROUTINE process_dwdana

  ! Reads the data from the first-guess and analysis files, and does any required processing of that input data.
  SUBROUTINE process_input_data(p_patch, inputInstructions, p_nh_state, p_int_state, p_grf_state, ext_data, prm_diag, p_lnd_state)
    TYPE(t_patch), INTENT(INOUT) :: p_patch(:)
    TYPE(t_readInstructionListPtr) :: inputInstructions(n_dom)
    TYPE(t_nh_state), INTENT(INOUT) :: p_nh_state(:)
    TYPE(t_int_state), INTENT(IN) :: p_int_state(:)
    TYPE(t_gridref_state), INTENT(IN) :: p_grf_state(:)
    TYPE(t_external_data), INTENT(INOUT) :: ext_data(:)
    TYPE(t_nwp_phy_diag), INTENT(INOUT), OPTIONAL :: prm_diag(:)
    TYPE(t_lnd_state), INTENT(INOUT), OPTIONAL :: p_lnd_state(:)

    CHARACTER(LEN = *), PARAMETER :: routine = modname//":process_input_data"

    CALL read_dwdfg(p_patch, inputInstructions, p_nh_state, prm_diag, p_lnd_state)
    CALL process_dwdfg(p_patch, p_nh_state, p_int_state, p_grf_state, ext_data, p_lnd_state, prm_diag)

    CALL read_dwdana(p_patch, inputInstructions, p_nh_state, p_lnd_state)
    ! process DWD analysis data
    CALL process_dwdana(p_patch, inputInstructions, p_nh_state, p_int_state, p_grf_state, ext_data, p_lnd_state)
  END SUBROUTINE process_input_data

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

    INTEGER :: jc,je,jk,jb,jg,jt          ! loop indices
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
    ! to sum up the water loading term
    REAL(wp), ALLOCATABLE, DIMENSION(:,:)   :: z_qsum

    REAL(wp) :: vn_incr_smt

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
      routine = modname//':create_dwdana_atm'

    ! nondimensional diffusion coefficient for interpolated velocity increment
    REAL(wp), PARAMETER :: smtfac=0.015_wp

    !-------------------------------------------------------------------------

    DO jg = 1, n_dom

      IF (.NOT. p_patch(jg)%ldom_active) CYCLE

      ! number of vertical levels
      nlev      = p_patch(jg)%nlev
      nlevp1    = p_patch(jg)%nlevp1

      nblks_c   = p_patch(jg)%nblks_c
      nblks_e   = p_patch(jg)%nblks_e
      i_nchdom  = MAX(1,p_patch(jg)%n_childdom)


      ! allocate temporary arrays for nonhydrostatic pressure, DA increments and a
      ! filtering term for vn
      ! note that an explicit temperature increment is not required (see below)
      ALLOCATE(zpres_nh (nproma,nlev,nblks_c),  &
               pres_incr(nproma,nlev,nblks_c),  &
               u_incr   (nproma,nlev,nblks_c),  &
               v_incr   (nproma,nlev,nblks_c),  &
               vn_incr  (nproma,nlev,nblks_e),  &
               w_incr   (nproma,nlevp1,nblks_c),&
               nabla4_vn_incr(nproma,nlev,nblks_e), &
               z_qsum   (nproma,nlev),          &
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


!$OMP DO PRIVATE(jb,jk,jc,jt,i_startidx,i_endidx,z_qsum)
      DO jb = i_startblk, i_endblk

        CALL get_indices_c(p_patch(jg), jb, i_startblk, i_endblk, &
          & i_startidx, i_endidx, rl_start, rl_end)

        ! Sum up the hydrometeor species for the water loading term
        z_qsum(:,:) = 0._wp
        DO jt = 2, iqm_max
          DO jk = 1, nlev
            DO jc = i_startidx, i_endidx
              z_qsum(jc,jk) = z_qsum(jc,jk) + p_prog_now_rcf%tracer(jc,jk,jb,jt)
            ENDDO
          ENDDO
        ENDDO

        DO jk = 1, nlev
          DO jc = i_startidx, i_endidx

            !******** CONSISTENCY CHECK ************
            !
            ! make sure, that due to GRIB2 roundoff errors, qv does not drop
            ! below threshhold (currently 5E-7 kg/kg)
            ! Alternative would be to increase writing precision for qv (DATATYPE_PACK24)
            ! Note: So far we are not fully convinced that the observed 'zeros' are
            ! soleyly a result of GRIB2 roundoff errors. They might also result from some
            ! numerical artifacts.
            p_prog_now_rcf%tracer(jc,jk,jb,iqv) = MAX(5.E-7_wp,                          &
             &                                       p_prog_now_rcf%tracer(jc,jk,jb,iqv))
            !******** END CONSISTENCY CHECK ********

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
              &                   / (1._wp + vtmpc1*p_prog_now_rcf%tracer(jc,jk,jb,iqv) - z_qsum(jc,jk))

          ENDDO  ! jc
        ENDDO  ! jk

      ENDDO  ! jb
!$OMP END DO
!$OMP END PARALLEL

      ! Recompute the hydrostatically integrated pressure from the first guess
      CALL diagnose_pres_temp (p_nh_state(jg)%metrics, p_prog_now, p_prog_now_rcf, p_diag, &
        &                      p_patch(jg), opt_calc_temp=.FALSE., opt_calc_pres=.TRUE.)


      IF (lread_ana) THEN

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


        ! include boundary interpolation zone of nested domains and the halo edges
        ! as far as possible
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

        ! required to avoid crash in nabla4_vec
        CALL sync_patch_array(SYNC_E,p_patch(jg),vn_incr)

        ! Compute diffusion term
        CALL nabla4_vec(vn_incr, p_patch(jg), p_int_state(jg), nabla4_vn_incr, opt_rlstart=5)

        ! Compute vertical wind increment consistent with the vn increment
        ! (strictly spoken, this should be done after the filtering step,
        ! but the difference is negligible)
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
              vn_incr_smt = vn_incr(je,jk,jb)   &
                &         - smtfac*nabla4_vn_incr(je,jk,jb)*p_patch(jg)%edges%area_edge(je,jb)**2

              ! add vn_incr_smt to first guess
              p_prog_now%vn(je,jk,jb) = p_prog_now%vn(je,jk,jb) + vn_incr_smt

            ENDDO  ! je
          ENDDO  ! jk

        ENDDO  ! jb
!$OMP ENDDO

        ! include boundary interpolation zone of nested domains but no halo points
        ! (sync follows below)
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
        IF ( iqc /= 0 .AND. iqi /= 0 .AND. iqr /= 0 .AND. iqs /= 0 ) THEN
          CALL virtual_temp(p_patch=p_patch(jg),             &
            &               temp=initicon(jg)%atm%temp,      & !in
            &               qv=p_prog_now%tracer(:,:,:,iqv), & !in
            &               qc=p_prog_now%tracer(:,:,:,iqc), & !in
            &               qi=p_prog_now%tracer(:,:,:,iqi), & !in
            &               qr=p_prog_now%tracer(:,:,:,iqr), & !in
            &               qs=p_prog_now%tracer(:,:,:,iqs), & !in
            &               temp_v=p_diag%tempv              ) !out
        ELSE IF ( iqc /= 0 .AND. iqi /= 0 ) THEN
          CALL virtual_temp(p_patch=p_patch(jg),             &
            &               temp=initicon(jg)%atm%temp,      & !in
            &               qv=p_prog_now%tracer(:,:,:,iqv), & !in
            &               qc=p_prog_now%tracer(:,:,:,iqc), & !in
            &               qi=p_prog_now%tracer(:,:,:,iqi), & !in
            &               temp_v=p_diag%tempv              ) !out
        ELSE
          CALL virtual_temp(p_patch=p_patch(jg),             &
            &               temp=initicon(jg)%atm%temp,      & !in
            &               qv=p_prog_now%tracer(:,:,:,iqv), & !in
            &               temp_v=p_diag%tempv              ) !out
        END IF

        ! Convert thermodynamic variables into set of NH prognostic variables
        CALL convert_thdvars(p_patch(jg), zpres_nh,  & !in
          &                  p_diag%tempv,           & !in
          &                  p_prog_now%rho,         & !out
          &                  p_prog_now%exner,       & !out
          &                  p_prog_now%theta_v      ) !out


      ENDIF  ! lread_ana


      ! deallocate temporary arrays
      DEALLOCATE( zpres_nh, pres_incr, u_incr, v_incr, vn_incr, nabla4_vn_incr, w_incr, z_qsum, STAT=ist )
      IF (ist /= SUCCESS) THEN
        CALL finish ( TRIM(routine), 'deallocation of auxiliary arrays failed' )
      ENDIF

    ENDDO  ! jg domain loop

  END SUBROUTINE create_dwdana_atm




  !>
  !! Compute analysis increments in terms of the NH prognostic set of variables.
  !!
  !!
  !! Compute analysis increments in terms of the NH prognostic set of variables
  !! (atmosphere only).
  !!
  !! @par Revision History
  !! Initial version by Daniel Reinert, DWD(2014-01-28)
  !!
  SUBROUTINE create_dwdanainc_atm (p_patch, p_nh_state, p_int_state)

    TYPE(t_patch),    TARGET, INTENT(IN)    :: p_patch(:)
    TYPE(t_nh_state), TARGET, INTENT(INOUT) :: p_nh_state(:)
    TYPE(t_int_state),        INTENT(IN)    :: p_int_state(:)

    INTEGER :: jc,je,jk,jb,jg,jt          ! loop indices
    INTEGER :: ist
    INTEGER :: nlev, nlevp1               ! number of vertical levels
    INTEGER :: nblks_c, nblks_e           ! number of blocks
    INTEGER :: i_nchdom
    INTEGER :: rl_start, rl_end
    INTEGER :: i_startblk, i_endblk
    INTEGER :: i_startidx, i_endidx
    TYPE(t_nh_prog), POINTER :: p_prog_now, p_prog_now_rcf
    TYPE(t_nh_diag), POINTER :: p_diag
    TYPE(t_nh_metrics), POINTER :: p_metrics
    INTEGER,         POINTER :: iidx(:,:,:), iblk(:,:,:), iqidx(:,:,:), iqblk(:,:,:)

    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: tempv_incr, nabla2_vn_incr, w_incr
    REAL(vp), ALLOCATABLE, DIMENSION(:,:,:) :: zvn_incr

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
      routine = modname//':create_dwdanainc_atm'

    ! nondimensional diffusion coefficient for interpolated velocity increment
    REAL(wp), PARAMETER :: ddfac=0.1_wp, smtfac=0.075_wp

    ! to sum up the water loading term
    REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: z_qsum

    ! For vertical filtering of mass increments
    REAL(wp), ALLOCATABLE :: rho_incr_smt(:,:), exner_ifc_incr(:,:), mass_incr_smt(:,:), mass_incr(:,:)
    REAL(wp) :: mass_incr_int(nproma), mass_incr_smt_int(nproma)
    INTEGER :: iter

    !-------------------------------------------------------------------------

    DO jg = 1, n_dom

      IF (.NOT. p_patch(jg)%ldom_active) CYCLE

      ! number of vertical levels
      nlev      = p_patch(jg)%nlev
      nlevp1    = p_patch(jg)%nlevp1

      nblks_c   = p_patch(jg)%nblks_c
      nblks_e   = p_patch(jg)%nblks_e
      i_nchdom  = MAX(1,p_patch(jg)%n_childdom)

      ! allocate temporary arrays for DA increments and a filtering term for vn
      ALLOCATE(tempv_incr(nproma,nlev,nblks_c), nabla2_vn_incr(nproma,nlev,nblks_e),   &
               rho_incr_smt(nproma,nlev), exner_ifc_incr(nproma,nlevp1),               &
               mass_incr_smt(nproma,nlev), mass_incr(nproma,nlev), z_qsum(nproma,nlev),&
               zvn_incr(nlev,nproma,nblks_e),STAT=ist)

      IF (ist /= SUCCESS) THEN
        CALL finish ( TRIM(routine), 'allocation of auxiliary arrays failed')
      ENDIF


      ! define some pointers
      p_prog_now     => p_nh_state(jg)%prog(nnow(jg))
      p_prog_now_rcf => p_nh_state(jg)%prog(nnow_rcf(jg))
      p_diag         => p_nh_state(jg)%diag
      p_metrics      => p_nh_state(jg)%metrics
      iidx           => p_patch(jg)%edges%cell_idx
      iblk           => p_patch(jg)%edges%cell_blk
      iqidx          => p_patch(jg)%edges%quad_idx
      iqblk          => p_patch(jg)%edges%quad_blk


      ! 1) Compute analysis increments in terms of the NH prognostic set of variables.
      !    Increments are computed for vn, w, exner, rho, qv. Note that a theta_v
      !    increment is not necessary.
      !    The prognostic state variables are initialized with the first guess
      !
      !
!$OMP PARALLEL PRIVATE(rl_start,rl_end,i_startblk,i_endblk)

      ! include boundary interpolation zone of nested domains and halo points
      rl_start = 1
      rl_end   = min_rlcell

      i_startblk = p_patch(jg)%cells%start_blk(rl_start,1)
      i_endblk   = p_patch(jg)%cells%end_blk(rl_end,i_nchdom)


!$OMP DO PRIVATE(jb,jk,jc,jt,i_startidx,i_endidx,rho_incr_smt,exner_ifc_incr,mass_incr_smt,mass_incr, &
!$OMP            mass_incr_int,mass_incr_smt_int,z_qsum)
      DO jb = i_startblk, i_endblk

        CALL get_indices_c(p_patch(jg), jb, i_startblk, i_endblk, &
          & i_startidx, i_endidx, rl_start, rl_end)

        ! Sum up the hydrometeor species for the water loading term
        z_qsum(:,:) = 0._wp
        DO jt = 2, iqm_max
          DO jk = 1, nlev
            DO jc = i_startidx, i_endidx
              z_qsum(jc,jk) = z_qsum(jc,jk) + p_prog_now_rcf%tracer(jc,jk,jb,jt)
            ENDDO
          ENDDO
        ENDDO

        DO jk = 1, nlev
          DO jc = i_startidx, i_endidx

            ! compute exner function (based on first guess input)
            p_prog_now%exner(jc,jk,jb) = (rd/p0ref * p_prog_now%rho(jc,jk,jb)  &
              &                        * p_prog_now%theta_v(jc,jk,jb))**(rd/cvd)

            ! compute full nonhydrostatic pressure from exner (based on first guess input)
            ! required for exner- and rho increment
            p_diag%pres(jc,jk,jb) = p0ref * (p_prog_now%exner(jc,jk,jb)**(cpd/rd))


            ! compute virtual temperature (based on first guess input)
            ! required for rho-increment
            p_diag%tempv(jc,jk,jb) = p_prog_now%theta_v(jc,jk,jb) &
              &                    * p_prog_now%exner(jc,jk,jb)

            ! compute temperature (based on first guess input)
            ! required for virtual temperature increment
            p_diag%temp(jc,jk,jb) = p_diag%tempv(jc,jk,jb)  &
              &                   / (1._wp + vtmpc1*p_prog_now_rcf%tracer(jc,jk,jb,iqv) - z_qsum(jc,jk))


            ! compute thermodynamic increments
            !
            p_diag%exner_incr(jc,jk,jb) = rd_o_cpd * p_prog_now%exner(jc,jk,jb) &
              &                  / p_diag%pres(jc,jk,jb) * initicon(jg)%atm_inc%pres(jc,jk,jb)

            tempv_incr(jc,jk,jb) = (1._wp + vtmpc1*p_prog_now_rcf%tracer(jc,jk,jb,iqv) - z_qsum(jc,jk)) &
              &                   * initicon(jg)%atm_inc%temp(jc,jk,jb)                                 &
              &                   + vtmpc1*p_diag%temp(jc,jk,jb)*initicon(jg)%atm_inc%qv(jc,jk,jb)

            p_diag%rho_incr(jc,jk,jb) = ( initicon(jg)%atm_inc%pres(jc,jk,jb) &
              &                / (rd*p_diag%tempv(jc,jk,jb)) )         &
              &                - ((p_diag%pres(jc,jk,jb) * tempv_incr(jc,jk,jb)) &
              &                / (rd*(p_diag%tempv(jc,jk,jb)**2)))

            p_diag%qv_incr(jc,jk,jb) = initicon(jg)%atm_inc%qv(jc,jk,jb)

          ENDDO  ! jc
        ENDDO  ! jk

        ! Apply vertical filtering of density increments, including a correction that ensures conservation
        ! of the vertically integrated mass increment
        !
        DO jk = 1, nlev
          DO jc = i_startidx, i_endidx
            mass_incr(jc,jk) = p_diag%rho_incr(jc,jk,jb)*p_metrics%ddqz_z_full(jc,jk,jb)
          ENDDO
        ENDDO

        ! Filtered density increment
        DO jk = 2, nlev-1
          DO jc = i_startidx, i_endidx
            rho_incr_smt(jc,jk) = (1._wp-2._wp*rho_incr_filter_wgt)*p_diag%rho_incr(jc,jk,jb) &
              + rho_incr_filter_wgt*(p_diag%rho_incr(jc,jk-1,jb)+p_diag%rho_incr(jc,jk+1,jb))
          ENDDO
        ENDDO
        DO jc = i_startidx, i_endidx
          rho_incr_smt(jc,1) = (1._wp-rho_incr_filter_wgt)*p_diag%rho_incr(jc,1,jb) &
            + rho_incr_filter_wgt*p_diag%rho_incr(jc,2,jb)
          rho_incr_smt(jc,nlev) = (1._wp-rho_incr_filter_wgt)*p_diag%rho_incr(jc,nlev,jb) &
            + rho_incr_filter_wgt*p_diag%rho_incr(jc,nlev-1,jb)
          ! correction increment for Exner pressure (zero at surface)
          exner_ifc_incr(jc,nlevp1) = 0._wp
        ENDDO

        DO jk = 1, nlev
          DO jc = i_startidx, i_endidx
            mass_incr_smt(jc,jk) = rho_incr_smt(jc,jk)*p_metrics%ddqz_z_full(jc,jk,jb)
          ENDDO
        ENDDO

        ! Vertically integrated mass increments
        DO jc = i_startidx, i_endidx
          mass_incr_int(jc)     = SUM(mass_incr(jc,1:nlev))
          mass_incr_smt_int(jc) = SUM(mass_incr_smt(jc,1:nlev))
        ENDDO

        ! Apply conservation correction for vertically integrated mass increment.
        ! Unfortunately, we have to care about possible pathological cases here
        DO jc = i_startidx, i_endidx
          IF (mass_incr_int(jc)*mass_incr_smt_int(jc) > 0._wp .AND. grav*ABS(mass_incr_int(jc)) > 0.5_wp &
              .AND. grav*ABS(mass_incr_smt_int(jc)) > 0.5_wp) THEN
            ! Multiplicative correction if mass increments are larger than 0.5 Pa and have the same sign
            DO jk = 1, nlev
              rho_incr_smt(jc,jk) = rho_incr_smt(jc,jk)                         &
               * MAX(0.98_wp,MIN(1.02_wp,mass_incr_int(jc)/mass_incr_smt_int(jc)))
            ENDDO
          ENDIF
        ENDDO

        ! integrate Exner correction increment (assuming hydrostatic balance of perturbation) and add it to the existing one
        DO jk = nlev, 1, -1
          DO jc = i_startidx, i_endidx

            exner_ifc_incr(jc,jk) = exner_ifc_incr(jc,jk+1) - rd_o_cpd*grav                   &
              * p_prog_now%exner(jc,jk,jb)/p_diag%pres(jc,jk,jb)                              &
              * (rho_incr_smt(jc,jk)-p_diag%rho_incr(jc,jk,jb))*p_metrics%ddqz_z_full(jc,jk,jb)

            p_diag%exner_incr(jc,jk,jb) = p_diag%exner_incr(jc,jk,jb) &
              + 0.5_wp*(exner_ifc_incr(jc,jk)+exner_ifc_incr(jc,jk+1))

            p_diag%rho_incr(jc,jk,jb) = rho_incr_smt(jc,jk)
          ENDDO
        ENDDO

      ENDDO  ! jb
!$OMP END DO NOWAIT



      ! 2) compute vn increments (w increment neglected)
      !
      IF (.NOT. lp2cintp_incr(jg)) THEN ! If increments are interpolated from the parent domain (see below),
                                        ! they are already provided as vn increments

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
              initicon(jg)%atm_inc%vn(je,jk,jb) = p_int_state(jg)%c_lin_e(je,1,jb)       &
                &               *(initicon(jg)%atm_inc%u(iidx(je,jb,1),jk,iblk(je,jb,1)) &
                &               * p_patch(jg)%edges%primal_normal_cell(je,jb,1)%v1       &
                &               + initicon(jg)%atm_inc%v(iidx(je,jb,1),jk,iblk(je,jb,1)) &
                &               * p_patch(jg)%edges%primal_normal_cell(je,jb,1)%v2)      &
                &               + p_int_state(jg)%c_lin_e(je,2,jb)                       &
                &               *(initicon(jg)%atm_inc%u(iidx(je,jb,2),jk,iblk(je,jb,2)) &
                &               * p_patch(jg)%edges%primal_normal_cell(je,jb,2)%v1       &
                &               + initicon(jg)%atm_inc%v(iidx(je,jb,2),jk,iblk(je,jb,2)) &
                &               * p_patch(jg)%edges%primal_normal_cell(je,jb,2)%v2  )

            ENDDO  ! je
          ENDDO  ! jk

        ENDDO  ! jb
!$OMP ENDDO
      ENDIF
!$OMP END PARALLEL

      IF (lp2cintp_incr(jg)) THEN
        ! Interpolate wind increments from parent domain (includes synchronization)
        CALL interpolate_vn_increments(initicon, p_patch(jg)%parent_id, jg)
      ELSE ! apply synchronization
        CALL sync_patch_array(SYNC_E,p_patch(jg),initicon(jg)%atm_inc%vn)
      END IF

      p_diag%vn_incr(:,:,:) = initicon(jg)%atm_inc%vn(:,:,:)

      ! Apply diffusion on wind increment
      DO iter = 1, niter_diffu
        CALL nabla2_vec(REAL(p_diag%vn_incr,wp), p_patch(jg), p_int_state(jg), nabla2_vn_incr, opt_rlstart=3)

!$OMP PARALLEL PRIVATE(rl_start,rl_end,i_startblk,i_endblk)

        ! include boundary interpolation zone of nested domains but no halo points (sync follows below)
        rl_start = 3
        rl_end   = min_rledge_int

        i_startblk = p_patch(jg)%edges%start_blk(rl_start,1)
        i_endblk   = p_patch(jg)%edges%end_blk(rl_end,i_nchdom)

!$OMP DO PRIVATE(jb,jk,je,i_startidx,i_endidx)
        DO jb = i_startblk, i_endblk

          CALL get_indices_e(p_patch(jg), jb, i_startblk, i_endblk, &
                             i_startidx, i_endidx, rl_start, rl_end)

          DO jk = 1, nlev
            DO je = i_startidx, i_endidx
              ! computed filtered velocity increment
              p_diag%vn_incr(je,jk,jb) = p_diag%vn_incr(je,jk,jb) + smtfac * &
                nabla2_vn_incr(je,jk,jb)*p_patch(jg)%edges%area_edge(je,jb)

            ENDDO  ! je
          ENDDO  ! jk

        ENDDO  ! jb
!$OMP ENDDO
!$OMP END PARALLEL
        CALL sync_patch_array(SYNC_E,p_patch(jg),p_diag%vn_incr)
      ENDDO

      ! Apply divergence damping on wind increment
      ! This is done for the global domain only because the interpolation to the nested domain(s)
      ! comes after the filtering
      IF (.NOT. lp2cintp_incr(jg)) THEN
        DO iter = 1, niter_divdamp

!$OMP PARALLEL PRIVATE(rl_start,rl_end,i_startblk,i_endblk)

          rl_start = 2
          rl_end   = min_rledge_int-2

          i_startblk = p_patch(jg)%edges%start_blk(rl_start,1)
          i_endblk   = p_patch(jg)%edges%end_blk(rl_end,i_nchdom)

!$OMP DO PRIVATE(jb,jk,je,i_startidx,i_endidx)
          DO jb = i_startblk, i_endblk

            CALL get_indices_e(p_patch(jg), jb, i_startblk, i_endblk, &
                               i_startidx, i_endidx, rl_start, rl_end)

            DO je = i_startidx, i_endidx
!DIR$ IVDEP
              DO jk = 1, nlev
                ! computed filtered velocity increment
                zvn_incr(jk,je,jb) = p_diag%vn_incr(je,jk,jb) + ddfac*p_patch(jg)%edges%area_edge(je,jb) * &
                  ( p_int_state(jg)%geofac_grdiv(je,1,jb)*p_diag%vn_incr(je,jk,jb)                         &
                  + p_int_state(jg)%geofac_grdiv(je,2,jb)*p_diag%vn_incr(iqidx(je,jb,1),jk,iqblk(je,jb,1)) &
                  + p_int_state(jg)%geofac_grdiv(je,3,jb)*p_diag%vn_incr(iqidx(je,jb,2),jk,iqblk(je,jb,2)) &
                  + p_int_state(jg)%geofac_grdiv(je,4,jb)*p_diag%vn_incr(iqidx(je,jb,3),jk,iqblk(je,jb,3)) &
                  + p_int_state(jg)%geofac_grdiv(je,5,jb)*p_diag%vn_incr(iqidx(je,jb,4),jk,iqblk(je,jb,4)) )

              ENDDO  ! je
            ENDDO  ! jk

          ENDDO  ! jb
!$OMP ENDDO

          rl_start = 3
          rl_end   = min_rledge_int

          i_startblk = p_patch(jg)%edges%start_blk(rl_start,1)
          i_endblk   = p_patch(jg)%edges%end_blk(rl_end,i_nchdom)

!$OMP DO PRIVATE(jb,jk,je,i_startidx,i_endidx)
          DO jb = i_startblk, i_endblk

            CALL get_indices_e(p_patch(jg), jb, i_startblk, i_endblk, &
                               i_startidx, i_endidx, rl_start, rl_end)

            DO je = i_startidx, i_endidx
!DIR$ IVDEP
              DO jk = 1, nlev
                ! computed filtered velocity increment
                p_diag%vn_incr(je,jk,jb) = zvn_incr(jk,je,jb) + ddfac*p_patch(jg)%edges%area_edge(je,jb) * &
                  ( p_int_state(jg)%geofac_grdiv(je,1,jb)*zvn_incr(jk,je,jb)                         &
                  + p_int_state(jg)%geofac_grdiv(je,2,jb)*zvn_incr(jk,iqidx(je,jb,1),iqblk(je,jb,1)) &
                  + p_int_state(jg)%geofac_grdiv(je,3,jb)*zvn_incr(jk,iqidx(je,jb,2),iqblk(je,jb,2)) &
                  + p_int_state(jg)%geofac_grdiv(je,4,jb)*zvn_incr(jk,iqidx(je,jb,3),iqblk(je,jb,3)) &
                  + p_int_state(jg)%geofac_grdiv(je,5,jb)*zvn_incr(jk,iqidx(je,jb,4),iqblk(je,jb,4)) )

              ENDDO  ! je
            ENDDO  ! jk

          ENDDO  ! jb
!$OMP ENDDO
!$OMP END PARALLEL

          CALL sync_patch_array(SYNC_E,p_patch(jg),p_diag%vn_incr)

        ENDDO
      ENDIF

      ! Copy filtered increment back to initicon state (this is needed to pass the filtered
      ! field to the parent-to-child interpolation for the nested domains)
      initicon(jg)%atm_inc%vn(:,:,:) = p_diag%vn_incr(:,:,:)

      ! deallocate temporary arrays
      DEALLOCATE( tempv_incr, nabla2_vn_incr, exner_ifc_incr, rho_incr_smt, mass_incr_smt, &
                  mass_incr, z_qsum, zvn_incr, STAT=ist )
      IF (ist /= SUCCESS) THEN
        CALL finish ( TRIM(routine), 'deallocation of auxiliary arrays failed' )
      ENDIF

      !
      ! If IAU-window is chosen to be zero, then analysis increments are added in one go.
      !
      IF (dt_iau == 0._wp) THEN

        ! For the special case that increments are added in one go,
        ! compute vertical wind increment consistent with the vn increment
        ! Note that here the filtered velocity increment is used.
        ALLOCATE(w_incr(nproma,nlevp1,nblks_c), STAT=ist)
        IF (ist /= SUCCESS) THEN
          CALL finish ( TRIM(routine), 'allocation of auxiliary arrays failed')
        ENDIF

        CALL sync_patch_array(SYNC_E,p_patch(jg),p_diag%vn_incr)

        CALL init_w(p_patch(jg), p_int_state(jg), REAL(p_diag%vn_incr,wp), p_nh_state(jg)%metrics%z_ifc, w_incr)

!$OMP PARALLEL PRIVATE(rl_start,rl_end,i_startblk,i_endblk)

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

              p_prog_now%exner(jc,jk,jb) = p_prog_now%exner(jc,jk,jb) + p_diag%exner_incr(jc,jk,jb)

              p_prog_now%rho(jc,jk,jb) = p_prog_now%rho(jc,jk,jb) + p_diag%rho_incr(jc,jk,jb)

              ! make sure, that due to GRIB2 roundoff errors, qv does not drop
              ! below threshhold (currently 5E-7 kg/kg)
              p_prog_now_rcf%tracer(jc,jk,jb,iqv) = MAX(5.E-7_wp,p_prog_now_rcf%tracer(jc,jk,jb,iqv)  &
                &                                 + p_diag%qv_incr(jc,jk,jb) )

              ! Remember to update theta_v
              p_prog_now%theta_v(jc,jk,jb) = (p0ref/rd) * p_prog_now%exner(jc,jk,jb)**(cvd/rd) &
                &                          / p_prog_now%rho(jc,jk,jb)

            ENDDO  ! jc
          ENDDO  ! jk

        ENDDO  ! jb
!$OMP ENDDO


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
              p_prog_now%vn(je,jk,jb) = p_prog_now%vn(je,jk,jb) + p_diag%vn_incr(je,jk,jb)
            ENDDO  ! je
          ENDDO  ! jk

        ENDDO  ! jb
!$OMP ENDDO NOWAIT


        ! include boundary interpolation zone of nested domains but no halo points
        ! (sync follows below)
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

        CALL sync_patch_array(SYNC_C,p_patch(jg),p_prog_now%w)

        ! deallocate temporary arrays
        DEALLOCATE( w_incr, STAT=ist )
        IF (ist /= SUCCESS) THEN
          CALL finish ( TRIM(routine), 'deallocation of auxiliary arrays failed' )
        ENDIF
      ENDIF  ! dt_iau = 0


      ! Recompute the hydrostatically integrated pressure from the first guess
      CALL diagnose_pres_temp (p_nh_state(jg)%metrics, p_prog_now, p_prog_now_rcf, p_diag, &
        &                      p_patch(jg), opt_calc_temp=.FALSE., opt_calc_pres=.TRUE.)

    ENDDO  ! jg domain loop

  END SUBROUTINE create_dwdanainc_atm



  !-------------------------------------------------------------------------
  !>
  !! SUBROUTINE create_iau_sfc
  !!
  !! Add increments to time-shifted first guess in one go.
  !! Increments are added for:
  !! W_SO, H_SNOW, FRESHSNW
  !!
  !! Additioanl sanity checks are performed for
  !! W_SO, H_SNOW, FRESHSNW, RHO_SNOW
  !!
  !! @par Revision History
  !! Initial version by D. Reinert, DWD (2014-07-17)
  !!
  !!
  !-------------------------------------------------------------------------
  SUBROUTINE create_iau_sfc (p_patch, p_nh_state, p_lnd_state, ext_data)

    TYPE(t_patch)             ,INTENT(IN)    :: p_patch(:)
    TYPE(t_nh_state) , TARGET ,INTENT(IN)    :: p_nh_state(:)
    TYPE(t_lnd_state), TARGET ,INTENT(INOUT) :: p_lnd_state(:)
    TYPE(t_external_data)     ,INTENT(IN)    :: ext_data(:)

    INTEGER :: jg, jb, jt, jk, jc, ic    ! loop indices
    INTEGER :: nblks_c                   ! number of blocks
    INTEGER :: rl_start, rl_end
    INTEGER :: i_startidx, i_endidx
    INTEGER :: ist

    LOGICAL :: lerr                      ! error flag

    TYPE(t_nh_diag) , POINTER :: p_diag            ! shortcut to diag state
    TYPE(t_lnd_prog), POINTER :: lnd_prog_now      ! shortcut to prognostic land state
    TYPE(t_lnd_diag), POINTER :: lnd_diag          ! shortcut to diagnostic land state

    REAL(wp) :: h_snow_t_fg(nproma,ntiles_total)   ! intermediate storage of h_snow first guess
    REAL(wp) :: snowfrac_lim

    REAL(wp), PARAMETER :: min_hsnow_inc=0.001_wp  ! minimum hsnow increment (1mm absolute value)
                                                   ! in order to avoid grib precision problems

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
      routine = modname//':create_iau_sfc'
  !-------------------------------------------------------------------------

    DO jg = 1, n_dom

      IF (.NOT. p_patch(jg)%ldom_active) CYCLE

      nblks_c   = p_patch(jg)%nblks_c
      rl_start  = 1
      rl_end    = min_rlcell

      ! save some paperwork
      p_diag       =>p_nh_state(jg)%diag
      lnd_prog_now =>p_lnd_state(jg)%prog_lnd(nnow_rcf(jg))
      lnd_diag     =>p_lnd_state(jg)%diag_lnd

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jt,jk,ic,jc,i_startidx,i_endidx,lerr,h_snow_t_fg,snowfrac_lim,ist)
      DO jb = 1, nblks_c

        ! (re)-initialize error flag
        lerr=.FALSE.

        CALL get_indices_c(p_patch(jg), jb, 1, nblks_c, &
                           i_startidx, i_endidx, rl_start, rl_end)

        ! add W_SO increment to first guess and perform some sanity checks in terms of realistic
        ! maximum/minimum values
        !
        DO jt = 1, ntiles_total

          DO jk = 1, nlev_soil
            DO ic = 1, ext_data(jg)%atm%lp_count_t(jb,jt)
              jc = ext_data(jg)%atm%idx_lst_lp_t(ic,jb,jt)

              IF (lnd_prog_now%w_so_t(jc,jk,jb,jt) <= 1.e-10_wp .AND.  &
                  cporv(ext_data(jg)%atm%soiltyp(jc,jb)) > 1.e-9_wp) THEN
                ! This should only happen for a tile coldstart; in this case,
                ! set soil water content to 50% of pore volume on newly appeared (non-dominant) land points
                lnd_prog_now%w_so_t(jc,jk,jb,jt) = 0.5_wp*cporv(ext_data(jg)%atm%soiltyp(jc,jb))*dzsoil_icon(jk)
              ELSE ! add w_so increment from SMA
                lnd_prog_now%w_so_t(jc,jk,jb,jt) = lnd_prog_now%w_so_t(jc,jk,jb,jt)  &
                   &                              + initicon(jg)%sfc_inc%w_so(jc,jk,jb)
              ENDIF

              ! Safety limits:  min=air dryness point, max=pore volume
              ist = ext_data(jg)%atm%soiltyp(jc,jb)
              SELECT CASE(ist)

                CASE (3,4,5,6,7,8) ! soil types with non-zero water content
                lnd_prog_now%w_so_t(jc,jk,jb,jt) = MIN(dzsoil_icon(jk)*cporv(ist),                                  &
                  &                                MAX(lnd_prog_now%w_so_t(jc,jk,jb,jt), dzsoil_icon(jk)*cadp(ist)) )

                CASE (9,10) ! sea water, sea ice
                ! ERROR landpoint has soiltype sea water or sea ice
                lerr = .TRUE.
              END SELECT

            ENDDO  ! ic
          ENDDO  ! jk

        ENDDO  ! jt

        IF (lerr) THEN
          CALL finish(routine, "Landpoint has invalid soiltype (sea water or sea ice)")
        ENDIF


        IF (init_mode == MODE_IAU) THEN

          ! store a copy of FG field for subsequent consistency checks
          h_snow_t_fg(:,:) = lnd_diag%h_snow_t(:,jb,:)

          ! add h_snow and freshsnow increments onto respective first guess fields
          DO jt = 1, ntiles_total

            IF (ltile_coldstart .OR. .NOT. lsnowtile) THEN
              ! Initialize snowfrac with 1 for the time being (the proper initialization follows in nwp_surface_init)
              ! This is actually needed for lsnowtile=.TRUE. because the snow cover fraction is used below in this case
              lnd_diag%snowfrac_lc_t(:,jb,jt) = 1._wp
              lnd_diag%snowfrac_t(:,jb,jt)    = 1._wp
            ENDIF

            DO ic = 1, ext_data(jg)%atm%gp_count_t(jb,jt)
              jc = ext_data(jg)%atm%idx_lst_t(ic,jb,jt)

              IF (ABS(initicon(jg)%sfc_inc%h_snow(jc,jb)) < min_hsnow_inc) THEN
                ! h_snow increment is neglected in order to avoid artefacts due to GRIB2 precision limitation
                ! minimum height: 0m; maximum height: 40m
                lnd_diag%h_snow_t   (jc,jb,jt) = MIN(40._wp,MAX(0._wp,lnd_diag%h_snow_t(jc,jb,jt)))
              ELSE
                IF (lsnowtile .AND. (jt > ntiles_lnd .OR. ltile_coldstart) ) THEN
                  ! in case of tile warmstart, add increment to snow-covered tiles only, rescaled with the snow-cover fraction
                  ! for tile coldstart, the snow increment is added in the same way as without snow tiles
                  snowfrac_lim = MAX(0.01_wp, lnd_diag%snowfrac_lc_t(jc,jb,jt))
                  lnd_diag%h_snow_t   (jc,jb,jt) = MIN(40._wp,MAX(0._wp,lnd_diag%h_snow_t(jc,jb,jt) &
                    &                            + initicon(jg)%sfc_inc%h_snow(jc,jb)/snowfrac_lim ))
                ELSE IF (lsnowtile .AND. initicon(jg)%sfc_inc%h_snow(jc,jb) > 0._wp .AND. &
                         lnd_diag%snowfrac_lc_t(jc,jb,jt) == 0._wp .AND. jt <= ntiles_lnd) THEN
                  ! if new snow is generated by the snow analysis (snowfrac_lc_t = 0 means that no corresponding
                  ! snow-covered grid point is present), the snow is added to the snow-free tile point
                  ! for the time being. Transfer to the snow-covered counterpart grid point and rescaling
                  ! is conducted after the first call to TERRA
                  lnd_diag%h_snow_t   (jc,jb,jt) = MIN(40._wp,initicon(jg)%sfc_inc%h_snow(jc,jb))
                  !
                  ! very simple initialization of snow-cover fraction for first TERRA call,
                  ! avoiding that the snow-free point disappears immediately
                  lnd_diag%snowfrac_lc_t(jc,jb,jt)  = MIN(0.9_wp,20._wp*lnd_diag%h_snow_t(jc,jb,jt))
                  lnd_diag%snowfrac_t(jc,jb,jt)     = lnd_diag%snowfrac_lc_t(jc,jb,jt)
                ELSE ! no snowtiles, snowfrac is initialized in nwp_surface_init in this case
                  ! minimum height: 0m; maximum height: 40m
                  lnd_diag%h_snow_t   (jc,jb,jt) = MIN(40._wp,MAX(0._wp,lnd_diag%h_snow_t(jc,jb,jt) &
                    &                                             + initicon(jg)%sfc_inc%h_snow(jc,jb)))
                ENDIF
              ENDIF

              ! maximum freshsnow factor: 1
              ! minimum freshsnow factor: 0
              ! two-sided limitation of freshsnow increment to +/- max_freshsnow_inc (tuning parameter)
              lnd_diag%freshsnow_t(jc,jb,jt) = MIN(1._wp,lnd_diag%freshsnow_t(jc,jb,jt)                               &
                &                            + SIGN(                                                                  &
                &                                  MIN(max_freshsnow_inc,ABS(initicon(jg)%sfc_inc%freshsnow(jc,jb))), &
                &                                  initicon(jg)%sfc_inc%freshsnow(jc,jb)                              &
                &                                  ) )

              lnd_diag%freshsnow_t(jc,jb,jt) = MAX(0._wp,lnd_diag%freshsnow_t(jc,jb,jt))


              ! adjust t_g and qv_s to new snow coming from the analysis,
              ! i.e. at new snow points, where h_snow increment > 0, but h_snow from first guess = 0
              !
              IF ( (initicon(jg)%sfc_inc%h_snow(jc,jb) >= min_hsnow_inc) .AND. (h_snow_t_fg(jc,jt) == 0._wp)) THEN
                lnd_prog_now%t_snow_t(jc,jb,jt) = MIN(tmelt,lnd_prog_now%t_snow_t(jc,jb,jt))
                lnd_prog_now%t_g_t   (jc,jb,jt) = lnd_diag%snowfrac_t(jc,jb,jt) * lnd_prog_now%t_snow_t(jc,jb,jt) &
                  &                             + (1._wp - lnd_diag%snowfrac_t(jc,jb,jt)) * lnd_prog_now%t_so_t(jc,1,jb,jt)
                lnd_diag%qv_s_t      (jc,jb,jt) = spec_humi(sat_pres_ice(lnd_prog_now%t_g_t(jc,jb,jt)),&
                  &                               p_diag%pres_sfc(jc,jb) )
              ENDIF

            ENDDO  ! ic

          ENDDO  ! jt



          ! consistency checks for rho_snow
          DO jt = 1, ntiles_total
            DO ic = 1, ext_data(jg)%atm%lp_count_t(jb,jt)
              jc = ext_data(jg)%atm%idx_lst_lp_t(ic,jb,jt)

              ! fresh snow 'missed' by the model
              IF ( (h_snow_t_fg(jc,jt) == 0._wp)                        .AND. &
                &  (initicon(jg)%sfc_inc%h_snow(jc,jb) > min_hsnow_inc) .AND. &
                &  (initicon(jg)%sfc_inc%freshsnow(jc,jb) > 0._wp) ) THEN

                lnd_prog_now%rho_snow_t(jc,jb,jt) = crhosmaxf   ! maximum density of fresh snow (150 kg/m**3)
              ENDIF

              ! old snow that is re-created by the analysis (i.e. snow melted away too fast in the model)
              IF ( (h_snow_t_fg(jc,jt) == 0._wp)                        .AND. &
                &  (initicon(jg)%sfc_inc%h_snow(jc,jb) > min_hsnow_inc) .AND. &
                &  (initicon(jg)%sfc_inc%freshsnow(jc,jb) <= 0._wp) ) THEN

                ! it is then assumed that we have 'old' snow in the model
                lnd_prog_now%rho_snow_t(jc,jb,jt) = crhosmax_ml   ! maximum density of snow (400 kg/m**3)
                lnd_diag%freshsnow_t(jc,jb,jt)    = 0._wp
              ENDIF

            ENDDO  ! ic


            ! Re-diagnose w_snow
            ! This is done in terra_multlay_init anyway, however it is safer to have consistent fields right
            ! from the beginning.
            lnd_prog_now%w_snow_t(i_startidx:i_endidx,jb,jt) = 0._wp

            CALL get_wsnow(h_snow    = lnd_diag%h_snow_t(:,jb,jt),          &
              &            rho_snow  = lnd_prog_now%rho_snow_t(:,jb,jt),    &
              &            t_snow    = lnd_prog_now%t_snow_t(:,jb,jt),      &
              &            istart    = i_startidx,                          &
              &            iend      = i_endidx,                            &
              &            soiltyp   = ext_data(jg)%atm%soiltyp_t(:,jb,jt), &
              &            w_snow    = lnd_prog_now%w_snow_t(:,jb,jt) )

          ENDDO  ! jt

        ENDIF  ! MODE_IAU

      ENDDO  ! jb
!$OMP END DO
!$OMP END PARALLEL


    ENDDO  ! jg

  END SUBROUTINE create_iau_sfc


  !-------------------------------------------------------------------------
  !>
  !! SUBROUTINE create_dwdana_sfc
  !!
  !! Required input: patch, lnd_state
  !! Output is written on fields of NH state
  !!
  !! @par Revision History
  !! Initial version by P. Ripodas, DWD(2013-05)
  !! Modification by Daniel Reinert, DWD (2013-11-20)
  !! - add consistency checks for rho_snow and w_so
  !!
  !!
  !-------------------------------------------------------------------------
  SUBROUTINE create_dwdana_sfc (p_patch, p_lnd_state, ext_data, inputInstructions)

    TYPE(t_patch),    TARGET ,INTENT(IN)    :: p_patch(:)
    TYPE(t_lnd_state)        ,INTENT(INOUT) :: p_lnd_state(:)
    TYPE(t_external_data)    ,INTENT(INOUT) :: ext_data(:)
    TYPE(t_readInstructionListPtr) :: inputInstructions(n_dom)

    INTEGER :: jg, ic, jc, jk, jb, jt, jgch, ist          ! loop indices
    INTEGER :: ntlr
    INTEGER :: nblks_c
    REAL(wp):: missval
    INTEGER :: rl_start, rl_end
    INTEGER :: i_startidx, i_endidx, i_endblk
    LOGICAL :: lp_mask(nproma)
    REAL(wp):: z_t_seasfc(nproma)              ! temporary field containing both SST 
                                               ! and lake-surface temperatures 
  !-------------------------------------------------------------------------

    ! get CDImissval
    missval = cdiInqMissval()

    DO jg = 1, n_dom

      IF (.NOT. p_patch(jg)%ldom_active) CYCLE

      nblks_c   = p_patch(jg)%nblks_c
      ntlr      = nnow_rcf(jg)

      rl_start = 1
      rl_end   = min_rlcell


      ! check, whether t_so is read from analysis
      IF (lp2cintp_sfcana(jg)) THEN
        jgch = p_patch(jg)%parent_id
      ELSE
        jgch = jg
      ENDIF

      lanaread_tseasfc(jg) = ( inputInstructions(jgch)%ptr%sourceOfVar('t_seasfc') == kInputSourceAna .OR. &
              ANY((/kInputSourceAna,kInputSourceBoth/) == inputInstructions(jgch)%ptr%sourceOfVar('t_so')) )

!$OMP PARALLEL
!$OMP DO PRIVATE(jc,ic,jk,jb,jt,i_startidx,i_endidx,lp_mask,ist,z_t_seasfc) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = 1, nblks_c

        CALL get_indices_c(p_patch(jg), jb, 1, nblks_c, &
                           i_startidx, i_endidx, rl_start, rl_end)



        IF (lanaread_tseasfc(jg)) THEN
          !
          ! SST analysis (T_SO(0) or T_SEA) was read into initicon(jg)%sfc%sst.
          ! Now copy to diag_lnd%t_seasfc for sea water points (including ice-covered ones)
          !
!CDIR NODEP,VOVERTAKE,VOB
          DO ic = 1, ext_data(jg)%atm%sp_count(jb)
             jc = ext_data(jg)%atm%idx_lst_sp(ic,jb)
             p_lnd_state(jg)%diag_lnd%t_seasfc(jc,jb) = MAX(tf_salt,initicon(jg)%sfc%sst(jc,jb))
          END DO

        ELSE  ! SST is not read from the analysis
          !
          ! get SST from first guess T_G
          !
          DO ic = 1, ext_data(jg)%atm%sp_count(jb)
            jc = ext_data(jg)%atm%idx_lst_sp(ic,jb)
            p_lnd_state(jg)%diag_lnd%t_seasfc(jc,jb) =  &
              & MAX(tf_salt, p_lnd_state(jg)%prog_lnd(ntlr)%t_g_t(jc,jb,isub_water))
            ! Ensure that t_seasfc is filled with tf_salt on completely frozen ocean points;
            IF (p_lnd_state(jg)%diag_lnd%fr_seaice(jc,jb) >= 1._wp-frsi_min) &
              p_lnd_state(jg)%diag_lnd%t_seasfc(jc,jb) = tf_salt
          END DO

        ENDIF  ! lanaread_tseasfc


        ! construct temporary field containing both SST and lake-surface temperatures
        ! which is needed for initializing T_SO at pure water points
        z_t_seasfc(:) = 0._wp
        DO ic = 1, ext_data(jg)%atm%sp_count(jb)
          jc = ext_data(jg)%atm%idx_lst_sp(ic,jb)
          z_t_seasfc(jc) = p_lnd_state(jg)%diag_lnd%t_seasfc(jc,jb)
        END DO
        IF (llake) THEN
          DO ic = 1, ext_data(jg)%atm%fp_count(jb)
            jc = ext_data(jg)%atm%idx_lst_fp(ic,jb)
            z_t_seasfc(jc) = MAX(tmelt, p_lnd_state(jg)%prog_wtr(ntlr)%t_wml_lk(jc,jb))
          END DO
        ELSE
          DO ic = 1, ext_data(jg)%atm%fp_count(jb)
            jc = ext_data(jg)%atm%idx_lst_fp(ic,jb)
            z_t_seasfc(jc) = MAX(tmelt, p_lnd_state(jg)%prog_lnd(ntlr)%t_g_t(jc,jb,isub_lake))
          END DO
        ENDIF
        !
        ! Fill T_SO with SST analysis over pure water points
        !
        ! Compute mask field for land points
        lp_mask(:) = .FALSE.
        DO ic = 1, ext_data(jg)%atm%lp_count_t(jb,1)
          jc = ext_data(jg)%atm%idx_lst_lp_t(ic,jb,1)
          lp_mask(jc) = .TRUE.
        ENDDO
        !
        DO jt = 1, ntiles_total
          DO jk = 1, nlev_soil
            DO jc = i_startidx, i_endidx
              IF (.NOT. lp_mask(jc)) THEN
                p_lnd_state(jg)%prog_lnd(ntlr)%t_so_t(jc,jk,jb,jt) = z_t_seasfc(jc)
              ENDIF
            ENDDO
          ENDDO
        ENDDO



        !***********************************!
        ! Consistency checks                !
        !***********************************!

        DO jt = 1, ntiles_total

          ! Check consistency between w_snow and rho_snow
          !
!CDIR NODEP,VOVERTAKE,VOB
          DO ic = 1, ext_data(jg)%atm%lp_count_t(jb,jt)
             jc = ext_data(jg)%atm%idx_lst_lp_t(ic,jb,jt)

             IF ( (p_lnd_state(jg)%prog_lnd(ntlr)%rho_snow_t(jc,jb,jt) < crhosmin_ml)  &
               &  .AND. ( (ext_data(jg)%atm%fr_land(jc,jb) < 0.5_wp)  .OR.             &
               &          (p_lnd_state(jg)%prog_lnd(ntlr)%w_snow_t(jc,jb,jt) >0._wp) ) &
               & )  THEN

               ! re-initialize rho_snow_t with minimum density of fresh snow (taken from TERRA)
               p_lnd_state(jg)%prog_lnd(ntlr)%rho_snow_t(jc,jb,jt) = crhosmin_ml
             ENDIF
          ENDDO  ! ic

          IF (init_mode == MODE_ICONVREMAP) THEN

            ! Constrain both rho_snow and t_snow because initial fields interpolated from a coarser grid
            ! may suffer from missing values near coasts
            DO ic = 1, ext_data(jg)%atm%lp_count_t(jb,jt)
              jc = ext_data(jg)%atm%idx_lst_lp_t(ic,jb,jt)

              p_lnd_state(jg)%prog_lnd(ntlr)%rho_snow_t(jc,jb,jt) = &
                MAX(crhosmin_ml,p_lnd_state(jg)%prog_lnd(ntlr)%rho_snow_t(jc,jb,jt))

              p_lnd_state(jg)%prog_lnd(ntlr)%t_snow_t(jc,jb,jt) = &
                MIN(p_lnd_state(jg)%prog_lnd(ntlr)%t_snow_t(jc,jb,jt), p_lnd_state(jg)%prog_lnd(ntlr)%t_g_t(jc,jb,jt))
              IF (p_lnd_state(jg)%prog_lnd(ntlr)%t_snow_t(jc,jb,jt) < tmelt-10._wp) &
                p_lnd_state(jg)%prog_lnd(ntlr)%t_snow_t(jc,jb,jt) = &
                MAX(p_lnd_state(jg)%prog_lnd(ntlr)%t_snow_t(jc,jb,jt), p_lnd_state(jg)%prog_lnd(ntlr)%t_g_t(jc,jb,jt)-10._wp)

             ENDDO
           ENDIF


          ! Catch problematic coast cases: ICON-land but GME ocean for moisture
          !
          DO jk = 1, nlev_soil
!CDIR NODEP,VOVERTAKE,VOB
            DO ic = 1, ext_data(jg)%atm%lp_count_t(jb,jt)
               jc = ext_data(jg)%atm%idx_lst_lp_t(ic,jb,jt)
               IF ((p_lnd_state(jg)%prog_lnd(nnow_rcf(jg))%w_so_t(jc,jk,jb,jt) <= 0._wp)) THEN
                  ! set dummy value (50% of pore volume)
                  p_lnd_state(jg)%prog_lnd(nnow_rcf(jg))%w_so_t(jc,jk,jb,jt) = &
                    &  0.5_wp * cporv(ext_data(jg)%atm%soiltyp_t(jc,jb,jt)) * dzsoil_icon(jk)
               ENDIF
               ! And temperature for ICON-land but COSMODE ocean
               IF ((p_lnd_state(jg)%prog_lnd(nnow_rcf(jg))%t_so_t(jc,jk,jb,jt) <= 0._wp)) THEN
                  ! set to first layer value
                  p_lnd_state(jg)%prog_lnd(nnow_rcf(jg))%t_so_t(jc,jk,jb,jt) = &
                    & p_lnd_state(jg)%prog_lnd(nnow_rcf(jg))%t_so_t(jc,1,jb,jt)
               ENDIF
            ENDDO  ! ic

            ! w_so_t, t_so_t:
            ! Search for CDI missval and replace it by meaningful value
            ! Reason: GRIB2-output fails otherwise (cumbersome values), probably due to
            ! the huge value range.
            DO jc = i_startidx, i_endidx
               IF ((p_lnd_state(jg)%prog_lnd(nnow_rcf(jg))%w_so_t(jc,jk,jb,jt) == missval)) THEN
                 p_lnd_state(jg)%prog_lnd(nnow_rcf(jg))%w_so_t(jc,jk,jb,jt) = 0._wp
               ENDIF
               IF ((p_lnd_state(jg)%prog_lnd(nnow_rcf(jg))%t_so_t(jc,jk,jb,jt) == missval)) THEN
                 p_lnd_state(jg)%prog_lnd(nnow_rcf(jg))%t_so_t(jc,jk,jb,jt) = & ! 0._wp
                   p_lnd_state(jg)%prog_lnd(nnow_rcf(jg))%t_so_t(jc,1,jb,jt)
               ENDIF

            ENDDO  ! jc

            ! Coldstart for limited-area mode still needs limitation of soil moisture to the allowed range
            IF (l_limited_area .AND. jg == 1 .AND. .NOT. lread_ana) THEN

              DO jc = i_startidx, i_endidx
                ist = ext_data(jg)%atm%soiltyp(jc,jb)
                SELECT CASE(ist)
                CASE (3,4,5,6,7,8) ! soil types with non-zero water content
                p_lnd_state(jg)%prog_lnd(nnow_rcf(jg))%w_so_t(jc,jk,jb,jt) = MIN(dzsoil_icon(jk)*cporv(ist),    &
                  &  MAX(p_lnd_state(jg)%prog_lnd(nnow_rcf(jg))%w_so_t(jc,jk,jb,jt), dzsoil_icon(jk)*cadp(ist)) )
                END SELECT
              ENDDO

            ENDIF

          ENDDO  ! jk

        ENDDO  ! jt


        ! fr_seaice, h_ice, t_ice:
        ! Search for CDI missval and replace it by meaningful value
        ! Reason: GRIB2-output fails otherwise (cumbersome values), probably due to
        ! the huge value range.
        DO jc = i_startidx, i_endidx
          IF (p_lnd_state(jg)%diag_lnd%fr_seaice(jc,jb) == missval) THEN
            p_lnd_state(jg)%diag_lnd%fr_seaice(jc,jb) = 0._wp
          ENDIF
          IF (p_lnd_state(jg)%prog_wtr(nnow_rcf(jg))%h_ice(jc,jb) == missval) THEN
            p_lnd_state(jg)%prog_wtr(nnow_rcf(jg))%h_ice(jc,jb) = 0._wp
          ENDIF
          IF (p_lnd_state(jg)%prog_wtr(nnow_rcf(jg))%t_ice(jc,jb) == missval) THEN
            p_lnd_state(jg)%prog_wtr(nnow_rcf(jg))%t_ice(jc,jb) = tf_salt
          ENDIF
        ENDDO  ! jc

      END DO  ! jb
!$OMP END DO
!$OMP END PARALLEL

      ! Fill t_seasfc on nest boundary points (needed because the turbtran initialization done in nwp_phy_init
      ! includes nest boundary points)

      IF (jg > 1) THEN

        rl_start = 1
        rl_end   = grf_bdywidth_c
        i_endblk = p_patch(jg)%cells%end_block(rl_end)

        DO jb = 1, i_endblk

          CALL get_indices_c(p_patch(jg), jb, 1, i_endblk, i_startidx, i_endidx, rl_start, rl_end)

          DO jc = i_startidx, i_endidx
            IF (ext_data(jg)%atm%fr_land(jc,jb) <= 1-frlnd_thrhld) THEN ! grid points with non-zero water fraction
              IF (lanaread_tseasfc(jg)) THEN
                p_lnd_state(jg)%diag_lnd%t_seasfc(jc,jb) = MAX(tmelt,initicon(jg)%sfc%sst(jc,jb))
                IF (ext_data(jg)%atm%fr_lake(jc,jb) < frlake_thrhld) THEN
                  p_lnd_state(jg)%prog_lnd(ntlr)%t_g_t(jc,jb,isub_water) = p_lnd_state(jg)%diag_lnd%t_seasfc(jc,jb)
                ENDIF
              ELSE IF (ext_data(jg)%atm%fr_lake(jc,jb) >= frlake_thrhld) THEN
                p_lnd_state(jg)%diag_lnd%t_seasfc(jc,jb) = MAX(tmelt, p_lnd_state(jg)%prog_lnd(ntlr)%t_g_t(jc,jb,isub_lake))
              ELSE
                p_lnd_state(jg)%diag_lnd%t_seasfc(jc,jb) = MAX(tf_salt, p_lnd_state(jg)%prog_lnd(ntlr)%t_g_t(jc,jb,isub_water))
              ENDIF
            ENDIF
          ENDDO

        ENDDO

      ENDIF

      ! This sync is needed because of the subsequent neighbor point filling
      CALL sync_patch_array(SYNC_C,p_patch(jg),p_lnd_state(jg)%diag_lnd%t_seasfc)

      ! Initialization of t_g_t(:,:,isub_water) and t_s_t(:,:,isub_water)
      ! with t_seasfc is performed in mo_nwp_sfc_utils:nwp_surface_init (nnow and nnew)

    ENDDO  ! jg domain loop

  END SUBROUTINE create_dwdana_sfc
  !-------------------------------------------------------------------------


END MODULE mo_initicon

