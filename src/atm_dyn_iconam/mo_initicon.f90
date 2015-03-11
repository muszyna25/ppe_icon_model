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

  USE mo_kind,                ONLY: wp
  USE mo_io_units,            ONLY: filename_max
  USE mo_parallel_config,     ONLY: nproma
  USE mo_run_config,          ONLY: iqv, iqc, iqi, iqr, iqs, iqm_max, iforcing
  USE mo_dynamics_config,     ONLY: nnow, nnow_rcf
  USE mo_model_domain,        ONLY: t_patch
  USE mo_nonhydro_types,      ONLY: t_nh_state, t_nh_prog, t_nh_diag, t_nh_metrics
  USE mo_nwp_phy_types,       ONLY: t_nwp_phy_diag
  USE mo_nwp_lnd_types,       ONLY: t_lnd_state, t_lnd_prog
  USE mo_intp_data_strc,      ONLY: t_int_state
  USE mo_ext_data_types,      ONLY: t_external_data
  USE mo_grf_intp_data_strc,  ONLY: t_gridref_state
  USE mo_initicon_types,      ONLY: t_initicon_state
  USE mo_initicon_config,     ONLY: init_mode, dt_iau, nlev_in, wgtfac_geobal,    &
    &                               rho_incr_filter_wgt, lread_ana,               &
    &                               lp2cintp_incr, lp2cintp_sfcana, ltile_coldstart
  USE mo_impl_constants,      ONLY: SUCCESS, MAX_CHAR_LENGTH, max_dom, MODE_DWDANA,          &
    &                               MODE_DWDANA_INC, MODE_IAU_OLD, MODE_IFSANA, MODE_ICONVREMAP, &
    &                               MODE_COMBINED, MODE_COSMODE, min_rlcell, INWP,           &
    &                               min_rledge_int, min_rlcell_int, dzsoil_icon => dzsoil
  USE mo_physical_constants,  ONLY: rd, cpd, cvd, p0ref, vtmpc1, grav, rd_o_cpd, tmelt, tf_salt
  USE mo_exception,           ONLY: message, finish
  USE mo_grid_config,         ONLY: n_dom
  USE mo_nh_init_utils,       ONLY: convert_thdvars, init_w
  USE mo_util_phys,           ONLY: virtual_temp
  USE mo_lnd_nwp_config,      ONLY: nlev_soil, ntiles_total, llake, &
    &                               isub_lake, isub_water, isub_seaice
  USE mo_phyparam_soil,       ONLY: cporv, crhosmin_ml
  USE mo_nh_vert_interp,      ONLY: vert_interp_atm, vert_interp_sfc
  USE mo_intp_rbf,            ONLY: rbf_vec_interpol_cell
  USE mo_nh_diagnose_pres_temp,ONLY: diagnose_pres_temp
  USE mo_loopindices,         ONLY: get_indices_c, get_indices_e
  USE mo_sync,                ONLY: sync_patch_array, SYNC_E, SYNC_C
  USE mo_math_laplace,        ONLY: nabla4_vec
  USE mo_math_gradients,      ONLY: grad_fd_norm
  USE mo_math_constants,      ONLY: rad2deg
  USE mo_intp_rbf,            ONLY: rbf_vec_interpol_cell
  USE mo_cdi_constants,       ONLY: cdiDefAdditionalKey, cdiInqMissval
  USE mo_flake,               ONLY: flake_coldinit
  USE mo_initicon_utils,      ONLY: create_input_groups, fill_tile_points,                        &
                                    copy_initicon2prog_atm, copy_initicon2prog_sfc, allocate_initicon, &
                                    deallocate_initicon, deallocate_extana_atm, deallocate_extana_sfc
  USE mo_initicon_io,         ONLY: open_init_files, close_init_files, read_extana_atm, read_extana_sfc, &
                                    read_dwdfg_atm, read_dwdfg_sfc, read_dwdana_atm, read_dwdana_sfc,    &
                                    read_dwdfg_atm_ii
  USE mo_util_string,         ONLY: one_of                                    

  IMPLICIT NONE


  PRIVATE

  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_initicon'

  TYPE(t_initicon_state), ALLOCATABLE, TARGET :: initicon(:) 

  ! NetCDF file IDs / CDI stream IDs for first guess and analysis file
  INTEGER, ALLOCATABLE :: fileID_fg(:),   fileID_ana(:)

  ! file type (NetCDF/GRB2) for first guess and analysis file
  INTEGER, ALLOCATABLE :: filetype_fg(:), filetype_ana(:)

  ! filename: first guess
  CHARACTER(LEN=filename_max) :: dwdfg_file(max_dom)

  ! filename: analysis
  CHARACTER(LEN=filename_max) :: dwdana_file(max_dom)

  PUBLIC :: init_icon



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

    TYPE(t_patch),          INTENT(IN)              :: p_patch(:)
    TYPE(t_int_state),      INTENT(IN)              :: p_int_state(:)
    TYPE(t_gridref_state),  INTENT(IN)              :: p_grf_state(:)
    TYPE(t_nh_state),       INTENT(INOUT)           :: p_nh_state(:)
    TYPE(t_nwp_phy_diag),   INTENT(INOUT), OPTIONAL :: prm_diag(:)
    TYPE(t_lnd_state),      INTENT(INOUT), OPTIONAL :: p_lnd_state(:)
    TYPE(t_external_data),  INTENT(INOUT), OPTIONAL :: ext_data(:)

    INTEGER :: jg, ist
    INTEGER :: jb              ! block loop index
    INTEGER :: i_startblk, i_endblk
    INTEGER :: i_rlstart, i_rlend, i_nchdom

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
      routine = modname//':init_icon'



    ! Allocate initicon data type
    ALLOCATE (initicon(n_dom),                         &
      &       filetype_fg(n_dom), filetype_ana(n_dom), &
      &       fileID_fg(n_dom),   fileID_ana(n_dom),   &
      &       stat=ist)
    IF (ist /= SUCCESS)  CALL finish(TRIM(routine),'allocation for initicon failed')

    initicon(:)%atm_in%linitialized  = .FALSE.
    initicon(:)%sfc_in%linitialized  = .FALSE.
    initicon(:)%atm%linitialized     = .FALSE.
    initicon(:)%sfc%linitialized     = .FALSE.
    initicon(:)%atm_inc%linitialized = .FALSE.
    initicon(:)%sfc_inc%linitialized = .FALSE.

    ! Allocate memory for init_icon state
    CALL allocate_initicon (p_patch, initicon)


    ! Copy the topography fields and coordinate surfaces to initicon
    !
    DO jg = 1, n_dom
      initicon(jg)%topography_c(:,:) = ext_data(jg)%atm%topography_c(:,:)

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
    IF (ANY((/MODE_DWDANA,MODE_DWDANA_INC,MODE_IAU_OLD,MODE_COMBINED,MODE_COSMODE,MODE_ICONVREMAP/) == init_mode)) THEN ! read in DWD analysis
      CALL open_init_files(p_patch, fileID_fg, fileID_ana, filetype_fg, filetype_ana, &
        &                  dwdfg_file, dwdana_file)

      ! Generate lists of fields that must be read from FG/ANA files
      !
      DO jg = 1, n_dom
        IF (.NOT. p_patch(jg)%ldom_active) CYCLE
        CALL create_input_groups(p_patch(jg),                        &
          &   initicon(jg)%grp_vars_fg,  initicon(jg)%ngrp_vars_fg,  &
          &   initicon(jg)%grp_vars_ana, initicon(jg)%ngrp_vars_ana, &
          &   initicon(jg)%grp_vars_fg_default , initicon(jg)%ngrp_vars_fg_default,  &
          &   initicon(jg)%grp_vars_ana_default, initicon(jg)%ngrp_vars_ana_default, &
          &   init_mode)
      ENDDO

    END IF




    ! init ICON prognostic fields
    !
    SELECT CASE(init_mode)
    CASE(MODE_DWDANA)   ! read in DWD analysis

      CALL message(TRIM(routine),'MODE_DWD: perform initialization with DWD analysis')

      ! process DWD atmosphere analysis data
      CALL process_dwdana_atm (p_patch, p_nh_state, p_int_state, p_grf_state)

      ! process DWD land/surface analysis data
      CALL process_dwdana_sfc (p_patch, prm_diag, p_lnd_state, ext_data)

    CASE(MODE_ICONVREMAP)   ! read in ICON prognostic variables (DWD first-guess fields) and 
                            ! perform vertical remapping

      CALL message(TRIM(routine),'MODE_VREMAP: read ICON data and perform vertical remapping')

      ! process ICON (DWD) atmosphere first-guess data (having different vertical levels than the current grid)
      CALL process_dwdana_atm (p_patch, p_nh_state, p_int_state, p_grf_state)

      ! process DWD land/surface analysis data
      CALL process_dwdana_sfc (p_patch, prm_diag, p_lnd_state, ext_data)

    CASE (MODE_DWDANA_INC)

      CALL message(TRIM(routine),'MODE_DWDANA_INC: perform initialization with '// &
        &                        ' incremental analysis update')

      ! process DWD atmosphere analysis increments
      CALL process_dwdanainc_atm (p_patch, p_nh_state, p_int_state)

      ! process DWD land/surface analysis data
      CALL process_dwdana_sfc (p_patch, prm_diag, p_lnd_state, ext_data)

    CASE (MODE_IAU_OLD)

      CALL message(TRIM(routine),'MODE_IAU_OLD: perform initialization with '// &
        &                        ' incremental analysis update')

      ! process DWD atmosphere analysis increments
      CALL process_dwdanainc_atm (p_patch, p_nh_state, p_int_state)

      ! process DWD land/surface analysis (increments)
      CALL process_dwdanainc_sfc (p_patch, prm_diag, p_lnd_state, ext_data)

    CASE(MODE_IFSANA)   ! read in IFS analysis

      CALL message(TRIM(routine),'MODE_IFS: perform initialization with IFS analysis')

      ! process IFS atmosphere analysis data
      CALL process_extana_atm (p_patch, p_nh_state, p_int_state, p_grf_state, initicon)

      IF (iforcing == inwp) THEN
        ! process IFS land/surface analysis data
        CALL process_extana_sfc (p_patch, p_lnd_state, initicon, ext_data)
      END IF

    CASE(MODE_COMBINED,MODE_COSMODE)

      IF (init_mode == MODE_COMBINED) THEN
        CALL message(TRIM(routine),'MODE_COMBINED: IFS-atm + GME-soil')
      ELSE
        CALL message(TRIM(routine),'MODE_COSMODE: COSMO-atm + COSMO-soil')
      ENDIF

      ! process IFS atmosphere analysis data
      CALL process_extana_atm (p_patch, p_nh_state, p_int_state, p_grf_state, initicon)

      ! process DWD land/surface analysis
      CALL process_dwdana_sfc (p_patch, prm_diag, p_lnd_state, ext_data)

      ! Cold-start initialization of the fresh-water lake model FLake.
      ! The procedure is the same as in "int2lm".
      ! Note that no lake ice is assumed at the cold start.
      !
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
              &   nflkgb      = ext_data(jg)%atm%fp_count    (jb),    &  ! in
              &   idx_lst_fp  = ext_data(jg)%atm%idx_lst_fp(:,jb),    &  ! in
              &   depth_lk    = ext_data(jg)%atm%depth_lk  (:,jb),    &  ! in
              ! here, a proper estimate of the lake surface temperature is required; 
              ! as neither GME nor COSMO-DE data have tiles, T_SO(0) is the best estimate
              &   tskin       = p_lnd_state(jg)%prog_lnd(nnow_rcf(jg))%t_so_t(:,1,jb,1),&  ! in
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
          ENDDO  ! jb
!$OMP END DO
!$OMP END PARALLEL
        ENDDO  ! jg
      ENDIF

    CASE DEFAULT

      CALL finish(TRIM(routine), "Invalid operation mode!")
    END SELECT


    ! Deallocate initicon data type
    !
    CALL deallocate_initicon(initicon)
    CALL deallocate_extana_atm (initicon)
    CALL deallocate_extana_sfc (initicon)

    ! close first guess and analysis files and corresponding inventory lists
    ! 
    IF (ANY((/MODE_DWDANA,MODE_DWDANA_INC,MODE_IAU_OLD,MODE_COMBINED,MODE_COSMODE/) == init_mode)) THEN
      CALL close_init_files(fileID_fg, fileID_ana)
    END IF

    DEALLOCATE (initicon, filetype_fg, filetype_ana, fileID_fg, fileID_ana, stat=ist)
    IF (ist /= success) CALL finish(TRIM(routine),'deallocation for initicon failed')

    ! splitting of sea-points list into open water and sea-ice points could be placed 
    ! here, instead of nwp_phy_init/init_nwp_phy
    ! however, one needs to make sure that it is called for both restart and non-restart
    ! runs. Could not be included into mo_ext_data_state/init_index_lists due to its 
    ! dependence on p_diag_lnd.
!DR    CALL init_sea_lists(p_patch, ext_data, p_diag_lnd, lseaice)

  END SUBROUTINE init_icon


  !-------------
  !>
  !! SUBROUTINE process_dwdana_atm
  !! Initialization routine of icon:
  !! - Reads DWD first guess and analysis (atmosphere only).
  !! - First guess is converted to T, p, qx, u, v in order to compute DA 
  !!   increments. These increments are then used to compute the analysis 
  !!   fields in terms of the NH set of prognostic variables 
  !!   (i.e. vn, rho, exner, theta_v )
  !!
  !! @par Revision History
  !! Initial version by Daniel Reinert, DWD(2012-12-20)
  !!
  !!
  SUBROUTINE process_dwdana_atm (p_patch, p_nh_state, p_int_state, p_grf_state)

    TYPE(t_patch),          INTENT(IN)    :: p_patch(:)
    TYPE(t_nh_state),       INTENT(INOUT) :: p_nh_state(:)
    TYPE(t_int_state),      INTENT(IN)    :: p_int_state(:)
    TYPE(t_gridref_state),  INTENT(IN)    :: p_grf_state(:)

!-------------------------------------------------------------------------

    IF (init_mode == MODE_ICONVREMAP) THEN

      ! read DWD first guess for atmosphere and store to initicon input state variables
      ! (input data are allowed to have a different number of model levels than the current model grid)
      CALL read_dwdfg_atm_ii (p_patch, initicon, fileID_fg, filetype_fg, dwdfg_file)

      ! Perform vertical interpolation from input ICON grid to output ICON grid
      !
      CALL vert_interp_atm(p_patch, p_nh_state, p_int_state, p_grf_state, nlev_in, initicon, &
                           opt_convert_omega2w=.FALSE.)

      ! Finally copy the results to the prognostic model variables
      !
      CALL copy_initicon2prog_atm(p_patch, initicon, p_nh_state)

    ELSE

      ! read DWD first guess for atmosphere
      CALL read_dwdfg_atm (p_patch, p_nh_state, initicon, fileID_fg, filetype_fg, dwdfg_file)

      IF(lread_ana) &   ! read DWD analysis from DA for atmosphere
        CALL read_dwdana_atm(p_patch, p_nh_state, initicon, fileID_ana, filetype_ana, dwdana_file)

      ! merge first guess with DA analysis and 
      ! convert variables to the NH set of prognostic variables
      CALL create_dwdana_atm(p_patch, p_nh_state, p_int_state)

    ENDIF

  END SUBROUTINE process_dwdana_atm



  !-------------
  !>
  !! SUBROUTINE process_dwdanainc_atm
  !! Initialization routine of icon:
  !! - Reads DWD first guess for t=T-dt_ass/2 (atmosphere only).
  !! - Reads analysis incerements for t=T (atmosphere only) in terms of 
  !!   \Delta T, \Delta u, \Delta v, \Delta p, \Delta qx
  !! - Compute analysis increments in terms of the NH set of prognostic 
  !!   variables
  !!
  !! @par Revision History
  !! Initial version by Daniel Reinert, DWD(2014-01-28)
  !!
  !!
  SUBROUTINE process_dwdanainc_atm (p_patch, p_nh_state, p_int_state)

    TYPE(t_patch),          INTENT(IN)    :: p_patch(:)
    TYPE(t_nh_state),       INTENT(INOUT) :: p_nh_state(:)
    TYPE(t_int_state),      INTENT(IN)    :: p_int_state(:)

!-------------------------------------------------------------------------


    ! read DWD first guess and analysis from DA for atmosphere
    ! 
    CALL read_dwdfg_atm (p_patch, p_nh_state, initicon, fileID_fg, filetype_fg, dwdfg_file)

    IF(lread_ana) &   
     CALL read_dwdana_atm(p_patch, p_nh_state, initicon, fileID_ana, filetype_ana, dwdana_file)


    ! Compute DA increments in terms of the NH set 
    ! of prognostic variables
    CALL create_dwdanainc_atm(p_patch, p_nh_state, p_int_state)

  END SUBROUTINE process_dwdanainc_atm



  !-------------
  !>
  !! SUBROUTINE process_dwdana_sfc
  !! Initialization routine of icon:
  !! - Reads DWD first guess (land/surface only). Data are directly
  !!   written to the prognostic model variables
  !! - reads DWD analysis (land/surface only). Data are written 
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


!!$    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
!!$      routine = modname//':process_dwdana_sfc'

!-------------------------------------------------------------------------


    ! read DWD first guess and analysis for surface/land
    ! 
    CALL read_dwdfg_sfc (p_patch, prm_diag, p_lnd_state, initicon, fileID_fg, filetype_fg, dwdfg_file)
 
    IF(lread_ana) &   
     CALL read_dwdana_sfc(p_patch, p_lnd_state, initicon, fileID_ana, filetype_ana, dwdana_file)
   
    ! get SST from first soil level t_so (for sea and lake points)
    ! perform consistency checks
    CALL create_dwdana_sfc(p_patch, p_lnd_state, ext_data)

  END SUBROUTINE process_dwdana_sfc



  !-------------
  !>
  !! SUBROUTINE process_dwdanainc_sfc
  !! Initialization routine of icon:
  !! - Reads DWD first guess (land/surface only). Data are directly
  !!   written to the prognostic model variables
  !! - reads DWD analysis fields and analysis increments (land/surface only). 
  !!   Increments are written to intermediate initicon variables and 
  !!   lateron (create_iau_sfc) added in one go.
  !!
  !! @par Revision History
  !! Initial version by Daniel Reinert, DWD(2014-07-18)
  !!
  !!
  SUBROUTINE process_dwdanainc_sfc (p_patch, prm_diag, p_lnd_state, ext_data)

    TYPE(t_patch),          INTENT(IN)    :: p_patch(:)
    TYPE(t_nwp_phy_diag),   INTENT(INOUT) :: prm_diag(:)
    TYPE(t_lnd_state),      INTENT(INOUT) :: p_lnd_state(:)
    TYPE(t_external_data),  INTENT(INOUT) :: ext_data(:)



!-------------------------------------------------------------------------


    ! read DWD first guess and analysis for surface/land
    ! 
    CALL read_dwdfg_sfc (p_patch, prm_diag, p_lnd_state, initicon, fileID_fg, filetype_fg, dwdfg_file)

    ! In case of tile coldstart, fill sub-grid scale land and water points with reasonable data
    ! from neighboring grid points where possible
    IF (ntiles_total > 1 .AND. ltile_coldstart) THEN
      CALL fill_tile_points(p_patch, p_lnd_state, ext_data)
    ENDIF

    IF(lread_ana) &   
      CALL read_dwdana_sfc(p_patch, p_lnd_state, initicon, fileID_ana, filetype_ana, dwdana_file)


    ! Add increments to time-shifted first guess in one go.
    CALL create_iau_sfc (p_patch, p_lnd_state, ext_data)

    ! get SST from first soil level t_so (for sea and lake points)
    ! perform consistency checks
    CALL create_dwdana_sfc(p_patch, p_lnd_state, ext_data)

  END SUBROUTINE process_dwdanainc_sfc



  !-------------
  !>
  !! SUBROUTINE process_extana_atm
  !! Initialization routine of icon:
  !! - Reads external analysis data (IFS or COSMO; atmosphere only)
  !! - performs vertical interpolation from intermediate IFS2ICON grid to ICON 
  !!   grid and converts variables to the NH set of prognostic variables
  !! - finally copies the results to the prognostic model variables
  !!
  !! @par Revision History
  !! Initial version by Daniel Reinert, DWD(2012-12-19)
  !!
  !!
  SUBROUTINE process_extana_atm (p_patch, p_nh_state, p_int_state, p_grf_state, &
    &                            initicon)

    TYPE(t_patch),          INTENT(IN)    :: p_patch(:)
    TYPE(t_nh_state),       INTENT(INOUT) :: p_nh_state(:)
    TYPE(t_int_state),      INTENT(IN)    :: p_int_state(:)
    TYPE(t_gridref_state),  INTENT(IN)    :: p_grf_state(:)

    TYPE(t_initicon_state), INTENT(INOUT) :: initicon(:)


!!$    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
!!$      routine = 'mo_nh_initicons:process_extana_atm'

    LOGICAL :: lomega_in

!-------------------------------------------------------------------------


    ! read horizontally interpolated external analysis for atmosphere
    ! 
    CALL read_extana_atm(p_patch, initicon)

    IF (init_mode == MODE_COSMODE) THEN
      lomega_in = .FALSE. ! in this case, w is provided as input
    ELSE
      lomega_in = .TRUE.  ! from hydrostatic models, omega is provided instead of w
    ENDIF

    ! Perform vertical interpolation from intermediate IFS2ICON grid to ICON grid
    ! and convert variables to the NH set of prognostic variables
    !
    CALL vert_interp_atm(p_patch, p_nh_state, p_int_state, p_grf_state, nlev_in, initicon, &
                         opt_convert_omega2w=lomega_in)

    
    ! Finally copy the results to the prognostic model variables
    !
    CALL copy_initicon2prog_atm(p_patch, initicon, p_nh_state)


  END SUBROUTINE process_extana_atm




  !-------------
  !>
  !! SUBROUTINE process_extana_sfc
  !! Initialization routine of icon:
  !! - Reads external analysis data (surface/land only)
  !! - performs vertical interpolation from intermediate IFS2ICON grid to ICON 
  !!   grid and converts variables to the NH set of prognostic variables
  !! - finally copies the results to the prognostic model variables
  !!
  !! @par Revision History
  !! Initial version by Daniel Reinert, DWD(2012-12-19)
  !!
  !!
  SUBROUTINE process_extana_sfc (p_patch, p_lnd_state, initicon, ext_data)

    TYPE(t_patch),          INTENT(IN)    :: p_patch(:)
    TYPE(t_lnd_state),      INTENT(INOUT) :: p_lnd_state(:)

    TYPE(t_initicon_state), INTENT(INOUT) :: initicon(:)
    TYPE(t_external_data),  INTENT(INOUT) :: ext_data(:)


!!$    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
!!$      routine = modname//':process_extana_sfc'

!-------------------------------------------------------------------------


    ! read horizontally interpolated external (currently IFS) analysis for surface/land
    ! 
    CALL read_extana_sfc(p_patch, initicon)


    ! Perform vertical interpolation from intermediate IFS2ICON grid to ICON grid
    ! and convert variables to the NH set of prognostic variables
    !
    CALL vert_interp_sfc(p_patch, initicon)

    
    ! Finally copy the results to the prognostic model variables
    !
    CALL copy_initicon2prog_sfc(p_patch, initicon, p_lnd_state, ext_data)


  END SUBROUTINE process_extana_sfc



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
    INTEGER,         POINTER :: iidx(:,:,:), iblk(:,:,:)

    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: tempv_incr, nabla4_vn_incr, w_incr

    ! Auxiliaries for artificial geostrophic balancing of pressure increments in the tropical stratosphere
    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: gradp, dpdx, dpdy
    REAL(wp) :: uincgeo, zmc, wfac

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
      routine = modname//':create_dwdanainc_atm'

    ! nondimensional diffusion coefficient for interpolated velocity increment
    REAL(wp), PARAMETER :: smtfac=0.015_wp

    ! to sum up the water loading term
    REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: z_qsum

    ! For vertical filtering of mass increments
    REAL(wp), ALLOCATABLE :: rho_incr_smt(:,:), exner_ifc_incr(:,:), mass_incr_smt(:,:), mass_incr(:,:)
    REAL(wp) :: mass_incr_int(nproma), mass_incr_smt_int(nproma)

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
      ALLOCATE(tempv_incr(nproma,nlev,nblks_c), nabla4_vn_incr(nproma,nlev,nblks_e),   &
               rho_incr_smt(nproma,nlev), exner_ifc_incr(nproma,nlevp1),               &
               mass_incr_smt(nproma,nlev), mass_incr(nproma,nlev), z_qsum(nproma,nlev),&
               gradp(nproma,nlev,nblks_e), dpdx(nproma,nlev,nblks_c),                  &
               dpdy(nproma,nlev,nblks_c), STAT=ist)
      IF (ist /= SUCCESS) THEN
        CALL finish ( TRIM(routine), 'allocation of auxiliary arrays failed')
      ENDIF

      nabla4_vn_incr(:,:,:) = 0._wp

      ! define some pointers
      p_prog_now     => p_nh_state(jg)%prog(nnow(jg))
      p_prog_now_rcf => p_nh_state(jg)%prog(nnow_rcf(jg))
      p_diag         => p_nh_state(jg)%diag
      p_metrics      => p_nh_state(jg)%metrics
      iidx           => p_patch(jg)%edges%cell_idx
      iblk           => p_patch(jg)%edges%cell_blk

      IF (wgtfac_geobal > 0._wp) THEN
        CALL grad_fd_norm(initicon(jg)%atm_inc%pres, p_patch(jg), gradp)
        CALL rbf_vec_interpol_cell(gradp, p_patch(jg), p_int_state(jg), dpdx, dpdy)
        CALL sync_patch_array(SYNC_C,p_patch(jg),dpdy) ! dpdx is unused
      ENDIF

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
!$OMP            mass_incr_int,mass_incr_smt_int,z_qsum,uincgeo,zmc,wfac)
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
              .AND. grav*ABS(mass_incr_int(jc)) > 0.5_wp) THEN
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

        ! modify zonal wind increment in order to achieve geostrophic balance with the pressure increment
        IF (wgtfac_geobal > 0._wp) THEN
          DO jk = 1, nlev
            DO jc = i_startidx, i_endidx
              uincgeo = - dpdy(jc,jk,jb)/(p_prog_now%rho(jc,jk,jb)* &
                          SIGN(MAX(1.e-5_wp,p_patch(jg)%cells%f_c(jc,jb)),p_patch(jg)%cells%f_c(jc,jb)) )
              zmc = p_metrics%z_mc(jc,jk,jb)
              IF (zmc >= 20000._wp .AND. zmc <= 40000._wp) THEN
                wfac = 1._wp
              ELSE IF (zmc >= 15000._wp .AND. zmc <= 20000._wp) THEN
                wfac = (zmc - 15000._wp) / 5000._wp
              ELSE IF (zmc >= 40000._wp .AND. zmc <= 45000._wp) THEN
                wfac = (45000._wp -  zmc) / 5000._wp
              ELSE
                wfac = 0._wp
              ENDIF
              IF (rad2deg*ABS(p_patch(jg)%cells%center(jc,jb)%lat) > 15._wp) THEN
                wfac = 0._wp
              ELSE IF (rad2deg*ABS(p_patch(jg)%cells%center(jc,jb)%lat) > 10._wp) THEN
                wfac = wfac*(15._wp-rad2deg*ABS(p_patch(jg)%cells%center(jc,jb)%lat))/5._wp
              ELSE IF (rad2deg*ABS(p_patch(jg)%cells%center(jc,jb)%lat) < 5._wp) THEN
                wfac = wfac*rad2deg*ABS(p_patch(jg)%cells%center(jc,jb)%lat)/5._wp
              ENDIF
              wfac = wgtfac_geobal*wfac
              initicon(jg)%atm_inc%u(jc,jk,jb) = wfac*uincgeo + (1._wp-wfac)*initicon(jg)%atm_inc%u(jc,jk,jb)
            ENDDO
          ENDDO
        ENDIF

      ENDDO  ! jb
!$OMP END DO NOWAIT



      ! 2) compute vn increments (w increment neglected)
      !
      IF (.NOT. lp2cintp_incr(jg)) THEN ! If increments were interpolated from the parent domain,
                                        ! they have already been converted into vn increments

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

      ! required to avoid crash in nabla4_vec (in case of interpolation from the parent domain,
      !                                        the sync has already been done)
      IF (.NOT. lp2cintp_incr(jg)) THEN
        CALL sync_patch_array(SYNC_E,p_patch(jg),initicon(jg)%atm_inc%vn)
      ENDIF

      ! Compute diffusion term 
      CALL nabla4_vec(initicon(jg)%atm_inc%vn, p_patch(jg), p_int_state(jg), nabla4_vn_incr, opt_rlstart=5)

!$OMP PARALLEL PRIVATE(rl_start,rl_end,i_startblk,i_endblk)

      ! include boundary interpolation zone of nested domains but no halo points (sync follows below)
      rl_start = 2
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
            p_diag%vn_incr(je,jk,jb) = initicon(jg)%atm_inc%vn(je,jk,jb) &
              &               - smtfac*nabla4_vn_incr(je,jk,jb)*p_patch(jg)%edges%area_edge(je,jb)**2

          ENDDO  ! je
        ENDDO  ! jk

      ENDDO  ! jb
!$OMP ENDDO
!$OMP END PARALLEL

      ! deallocate temporary arrays
      DEALLOCATE( tempv_incr, nabla4_vn_incr, exner_ifc_incr, rho_incr_smt, mass_incr_smt, &
                  mass_incr, z_qsum, STAT=ist )
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

        CALL init_w(p_patch(jg), p_int_state(jg), p_diag%vn_incr, p_nh_state(jg)%metrics%z_ifc, w_incr)

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


    ENDDO  ! jg domain loop

  END SUBROUTINE create_dwdanainc_atm



  !-------------------------------------------------------------------------
  !>
  !! SUBROUTINE create_iau_sfc 
  !!
  !! Add increments to time-shifted first guess in one go.
  !! Increments are added for: 
  !! W_SO, H_SNOW, W_SNOW, FRESHSNW
  !!
  !! Additioanl sanity checks are performed for 
  !! W_SO 
  !!
  !! @par Revision History
  !! Initial version by D. Reinert, DWD (2014-07-17)
  !!
  !!
  !-------------------------------------------------------------------------
  SUBROUTINE create_iau_sfc (p_patch, p_lnd_state, ext_data)

    TYPE(t_patch)             ,INTENT(IN)    :: p_patch(:)
    TYPE(t_lnd_state), TARGET ,INTENT(INOUT) :: p_lnd_state(:)
    TYPE(t_external_data)     ,INTENT(IN)    :: ext_data(:)

    INTEGER :: jg, jb, jt, jk, jc, ic    ! loop indices
    INTEGER :: nblks_c                   ! number of blocks
    INTEGER :: rl_start, rl_end
    INTEGER :: i_startidx, i_endidx

    LOGICAL :: lerr                      ! error flag

    TYPE(t_lnd_prog), POINTER :: lnd_prog_now      ! shortcut to prognostic land state

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
      routine = modname//':create_iau_sfc'
  !-------------------------------------------------------------------------

    DO jg = 1, n_dom

      IF (.NOT. p_patch(jg)%ldom_active) CYCLE

      nblks_c   = p_patch(jg)%nblks_c
      rl_start  = 1
      rl_end    = min_rlcell

      ! save some paperwork
      lnd_prog_now =>p_lnd_state(jg)%prog_lnd(nnow_rcf(jg))


!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jt,jk,jc,i_startidx,i_endidx,lerr)
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
              SELECT CASE(ext_data(jg)%atm%soiltyp(jc,jb))

                CASE(3)  !sand
                lnd_prog_now%w_so_t(jc,jk,jb,jt) = MIN(dzsoil_icon(jk)*0.364_wp,                                   &
                  &                                MAX(lnd_prog_now%w_so_t(jc,jk,jb,jt), dzsoil_icon(jk)*0.012_wp) &
                  &                                )
                CASE(4)  !sandyloam
                lnd_prog_now%w_so_t(jc,jk,jb,jt) = MIN(dzsoil_icon(jk)*0.445_wp,                                   &
                  &                                MAX(lnd_prog_now%w_so_t(jc,jk,jb,jt), dzsoil_icon(jk)*0.03_wp)  &
                  &                                )
                CASE(5)  !loam
                lnd_prog_now%w_so_t(jc,jk,jb,jt) = MIN(dzsoil_icon(jk)*0.455_wp,                                   &
                  &                                MAX(lnd_prog_now%w_so_t(jc,jk,jb,jt), dzsoil_icon(jk)*0.035_wp) &
                  &                                )
                CASE(6)  !clayloam
                lnd_prog_now%w_so_t(jc,jk,jb,jt) = MIN(dzsoil_icon(jk)*0.475_wp,                                   &
                  &                                MAX(lnd_prog_now%w_so_t(jc,jk,jb,jt), dzsoil_icon(jk)*0.06_wp)  &
                  &                                )
                CASE(7)  !clay
                lnd_prog_now%w_so_t(jc,jk,jb,jt) = MIN(dzsoil_icon(jk)*0.507_wp,                                   &
                  &                                MAX(lnd_prog_now%w_so_t(jc,jk,jb,jt), dzsoil_icon(jk)*0.065_wp) &
                  &                                )
                CASE(8)  !peat
                lnd_prog_now%w_so_t(jc,jk,jb,jt) = MIN(dzsoil_icon(jk)*0.863_wp,                                   &
                  &                                MAX(lnd_prog_now%w_so_t(jc,jk,jb,jt), dzsoil_icon(jk)*0.098_wp) &
                  &                                )          
                CASE(9,10)!sea water, sea ice
                ! ERROR landpoint has soiltype sea water or sea ice
                lerr = .TRUE.
              END SELECT

            ENDDO  ! ic
          ENDDO  ! jk

        ENDDO  ! jt

        IF (lerr) THEN
          CALL finish(routine, "Landpoint has invalid soiltype (sea water or sea ice)")
        ENDIF

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
  SUBROUTINE create_dwdana_sfc (p_patch, p_lnd_state, ext_data)

    TYPE(t_patch),    TARGET ,INTENT(IN)    :: p_patch(:)
    TYPE(t_lnd_state)        ,INTENT(INOUT) :: p_lnd_state(:)
    TYPE(t_external_data)    ,INTENT(INOUT) :: ext_data(:)

    INTEGER :: jg, ic, jc, jk, jb, jt, jgch          ! loop indices
    INTEGER :: ntlr
    INTEGER :: nblks_c
    REAL(wp):: missval
    INTEGER :: rl_start, rl_end 
    INTEGER :: i_startidx, i_endidx
    LOGICAL :: lanaread_tso                    ! .TRUE. T_SO(0) was read from analysis
    LOGICAL :: lp_mask(nproma)
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
      lanaread_tso = ( one_of('t_so', initicon(jgch)%grp_vars_ana(1:initicon(jgch)%ngrp_vars_ana)) /= -1)

!$OMP PARALLEL 
!$OMP DO PRIVATE(jc,ic,jk,jb,jt,i_startidx,i_endidx,lp_mask) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = 1, nblks_c

        CALL get_indices_c(p_patch(jg), jb, 1, nblks_c, &
                           i_startidx, i_endidx, rl_start, rl_end)


        IF (lanaread_tso) THEN
          !
          ! SST analysis (T_SO(0)) was read into initicon(jg)%sfc%sst.
          ! Now copy to diag_lnd%t_seasfc for water, ice and lake points
          !
!CDIR NODEP,VOVERTAKE,VOB
          DO ic = 1, ext_data(jg)%atm%sp_count(jb)
             jc = ext_data(jg)%atm%idx_lst_sp(ic,jb)
             p_lnd_state(jg)%diag_lnd%t_seasfc(jc,jb) = MAX(tf_salt,initicon(jg)%sfc%sst(jc,jb))
          END DO
!CDIR NODEP,VOVERTAKE,VOB
          DO ic = 1, ext_data(jg)%atm%fp_count(jb)
           jc = ext_data(jg)%atm%idx_lst_fp(ic,jb)
           p_lnd_state(jg)%diag_lnd%t_seasfc(jc,jb) = MAX(tmelt,initicon(jg)%sfc%sst(jc,jb))
          END DO

          ! Compute mask field for land points
          lp_mask(:) = .FALSE.
          DO ic = 1, ext_data(jg)%atm%lp_count_t(jb,1)
            jc = ext_data(jg)%atm%idx_lst_lp_t(ic,jb,1)
            lp_mask(jc) = .TRUE.
          ENDDO

          ! Fill T_SO with SST analysis over pure water points
          DO jt = 1, ntiles_total
            DO jk = 1, nlev_soil
              DO jc = i_startidx, i_endidx
                IF (.NOT. lp_mask(jc)) THEN
                  p_lnd_state(jg)%prog_lnd(ntlr)%t_so_t(jc,jk,jb,jt) = initicon(jg)%sfc%sst(jc,jb)
                ENDIF
              ENDDO
            ENDDO
          ENDDO

        ELSE  ! SST (T_SO(0)) is not read from the analysis
          !
          ! get SST from first guess T_G (for sea and lake points)
          !
!CDIR NODEP,VOVERTAKE,VOB
          DO ic = 1, ext_data(jg)%atm%sp_count(jb)
            jc = ext_data(jg)%atm%idx_lst_sp(ic,jb)
            p_lnd_state(jg)%diag_lnd%t_seasfc(jc,jb) =  &
              & MAX(tf_salt, p_lnd_state(jg)%prog_lnd(ntlr)%t_g_t(jc,jb,isub_water))
          END DO
!CDIR NODEP,VOVERTAKE,VOB
          DO ic = 1, ext_data(jg)%atm%fp_count(jb)
            jc = ext_data(jg)%atm%idx_lst_fp(ic,jb)
            p_lnd_state(jg)%diag_lnd%t_seasfc(jc,jb) = &
              & MAX(tmelt, p_lnd_state(jg)%prog_lnd(ntlr)%t_g_t(jc,jb,isub_lake))
          END DO

        ENDIF  ! lanaread_t_so



        !***********************************!
        ! Consistency checks                !
        !***********************************!

        DO jt = 1, ntiles_total

          ! Check consistency between w_snow and rho_snow
          !
!CDIR NODEP,VOVERTAKE,VOB
          DO ic = 1, ext_data(jg)%atm%lp_count_t(jb,jt)
             jc = ext_data(jg)%atm%idx_lst_lp_t(ic,jb,jt)

             IF ( (p_lnd_state(jg)%prog_lnd(nnow_rcf(jg))%rho_snow_t(jc,jb,jt) < crhosmin_ml)  &
               &  .AND. (p_lnd_state(jg)%prog_lnd(nnow_rcf(jg))%w_snow_t(jc,jb,jt) >0._wp) )  THEN

               ! re-initialize rho_snow_t with minimum density of fresh snow (taken from TERRA)
               p_lnd_state(jg)%prog_lnd(nnow_rcf(jg))%rho_snow_t(jc,jb,jt) = crhosmin_ml
             ENDIF
          ENDDO  ! ic


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


    ! Initialization of t_g_t(:,:,isub_water) with t_seasfc is performed in 
    ! mo_nwp_sfc_utils:nwp_surface_init (nnow and nnew)

    ENDDO  ! jg domain loop

  END SUBROUTINE create_dwdana_sfc
  !-------------------------------------------------------------------------


END MODULE mo_initicon

