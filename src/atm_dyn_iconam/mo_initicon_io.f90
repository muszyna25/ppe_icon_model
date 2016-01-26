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

  USE mo_kind,                ONLY: wp, dp
  USE mo_io_units,            ONLY: filename_max
  USE mo_parallel_config,     ONLY: nproma, p_test_run
  USE mo_run_config,          ONLY: msg_level, iqv, iqc, iqi, iqr, iqs
  USE mo_dynamics_config,     ONLY: nnow, nnow_rcf
  USE mo_model_domain,        ONLY: t_patch
  USE mo_nonhydro_types,      ONLY: t_nh_state, t_nh_prog
  USE mo_nwp_phy_types,       ONLY: t_nwp_phy_diag
  USE mo_nwp_lnd_types,       ONLY: t_lnd_state, t_lnd_prog, t_lnd_diag, t_wtr_prog
  USE mo_initicon_types,      ONLY: t_initicon_state, t_pi_atm, alb_snow_var, geop_ml_var
  USE mo_input_instructions,  ONLY: t_readInstructionListPtr
  USE mo_initicon_config,     ONLY: init_mode, nlevatm_in, l_sst_in, generate_filename, &
    &                               ifs2icon_filename, dwdfg_filename, dwdana_filename, &
    &                               nml_filetype => filetype, lread_vn,      &
    &                               lp2cintp_incr, lp2cintp_sfcana, ltile_coldstart,    &
    &                               lvert_remap_fg, aerosol_fg_present
  USE mo_nh_init_nest_utils,  ONLY: interpolate_increments, interpolate_sfcana
  USE mo_impl_constants,      ONLY: MAX_CHAR_LENGTH, max_dom,                           &
    &                               MODE_IAU, MODE_IAU_OLD, MODE_IFSANA, MODE_COMBINED, &
    &                               MODE_COSMODE, iss, iorg, ibc, iso4, idu
  USE mo_exception,           ONLY: message, finish, message_text
  USE mo_grid_config,         ONLY: n_dom, nroot, l_limited_area
  USE mo_mpi,                 ONLY: p_io, p_bcast, p_comm_work,    &
    &                               p_comm_work_test, my_process_is_mpi_workroot
  USE mo_io_config,           ONLY: default_read_method
  USE mo_read_interface,      ONLY: t_stream_id, nf, openInputFile, closeFile, &
    &                               read_2d_1time, read_2d_1lev_1time, &
    &                               read_3d_1time, on_cells, on_edges
  USE mo_util_cdi,            ONLY: t_inputParameters, trivial_tileId
  USE mo_ifs_coord,           ONLY: alloc_vct, init_vct, vct, vct_a, vct_b
  USE mo_lnd_nwp_config,      ONLY: ntiles_total, &
    &                               ntiles_water, lmulti_snow, tiles
  USE mo_master_config,       ONLY: getModelBaseDir
  USE mo_var_metadata_types,  ONLY: VARNAME_LEN
  USE mo_nwp_sfc_interp,      ONLY: smi_to_wsoil
  USE mo_io_util,             ONLY: get_filetype
  USE mo_initicon_utils,      ONLY: allocate_extana_atm, allocate_extana_sfc
  USE mo_physical_constants,  ONLY: cpd, rd, cvd_o_rd, p0ref, vtmpc1
  USE mo_fortran_tools,       ONLY: init
  USE mo_input_request_list,  ONLY: t_InputRequestList
  USE mo_util_string,         ONLY: int2string

  IMPLICIT NONE

  INCLUDE 'netcdf.inc'

  PRIVATE

  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_initicon_io'


  PUBLIC :: fgFilename
  PUBLIC :: fgFiletype
  PUBLIC :: anaFilename
  PUBLIC :: anaFiletype

  PUBLIC :: read_extana_atm
  PUBLIC :: read_extana_sfc

  PUBLIC :: request_dwdfg_atm
  PUBLIC :: request_dwdfg_atm_ii
  PUBLIC :: request_dwdfg_sfc
  PUBLIC :: request_dwdana_atm
  PUBLIC :: request_dwdana_sfc

  PUBLIC :: fetch_dwdfg_atm
  PUBLIC :: fetch_dwdfg_atm_ii
  PUBLIC :: fetch_dwdfg_sfc
  PUBLIC :: fetch_dwdana_atm
  PUBLIC :: fetch_dwdana_sfc

! PUBLIC :: process_input_dwdfg_atm     !This is currently not needed.
! PUBLIC :: process_input_dwdfg_atm_ii
  PUBLIC :: process_input_dwdfg_sfc
  PUBLIC :: process_input_dwdana_atm
  PUBLIC :: process_input_dwdana_sfc



  CONTAINS

  FUNCTION fgFilename(p_patch) RESULT(RESULT)
    CHARACTER(LEN = filename_max) :: RESULT
    TYPE(t_patch), INTENT(IN) :: p_patch

    RESULT = generate_filename(dwdfg_filename, getModelBaseDir(), nroot, p_patch%level, p_patch%id)
  END FUNCTION fgFilename

  FUNCTION anaFilename(p_patch) RESULT(RESULT)
    CHARACTER(LEN = filename_max) :: RESULT
    TYPE(t_patch), INTENT(IN) :: p_patch

    RESULT = generate_filename(dwdana_filename, getModelBaseDir(), nroot, p_patch%level, p_patch%id)
  END FUNCTION anaFilename

  INTEGER FUNCTION fgFiletype() RESULT(RESULT)
    IF(nml_filetype == -1) THEN
        ! get_filetype() ONLY uses the suffix, which IS already a part of the template IN dwdfg_filename.
        ! This IS why it suffices to USE the dwdfg_filename directly here without expanding it first via generate_filename().
        RESULT = get_filetype(TRIM(dwdfg_filename))
    ELSE
        RESULT = nml_filetype
    END IF
  END FUNCTION fgFiletype

  INTEGER FUNCTION anaFiletype() RESULT(RESULT)
    IF(nml_filetype == -1) THEN
        ! get_filetype() ONLY uses the suffix, which IS already a part of the template IN dwdana_filename.
        ! This IS why it suffices to USE the dwdana_filename directly here without expanding it first via generate_filename().
        RESULT = get_filetype(TRIM(dwdana_filename))
    ELSE
        RESULT = nml_filetype
    END IF
  END FUNCTION anaFiletype


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

    TYPE(t_patch), TARGET,  INTENT(IN)    :: p_patch(:)
    TYPE(t_initicon_state), INTENT(INOUT) :: initicon(:)

    INTEGER :: jg, jlev, jc, jk, jb, i_endidx
    LOGICAL :: l_exist

    INTEGER :: no_cells, no_levels, nlev_in
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
    nlev_in = 0

    DO jg = 1, n_dom

      jlev = p_patch(jg)%level

      ! Skip reading the atmospheric input data if a model domain
      ! is not active at initial time
      IF (.NOT. p_patch(jg)%ldom_active) CYCLE

      !
      ! generate file name
      !
      ifs2icon_file(jg) = generate_filename(ifs2icon_filename, getModelBaseDir(), &
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
      ! open file (NetCDF file!)
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

      nlevatm_in(:) = nlev_in

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
      CALL read_3d_1time(stream_id, on_cells, 'T', fill_array=initicon(jg)%atm_in%temp)

      IF (lread_vn) THEN
        CALL read_3d_1time(stream_id, on_edges, 'VN', fill_array=initicon(jg)%atm_in%vn)
      ELSE
        CALL read_3d_1time(stream_id, on_cells, 'U', fill_array=initicon(jg)%atm_in%u)

        CALL read_3d_1time(stream_id, on_cells, 'V', fill_array=initicon(jg)%atm_in%v)
      ENDIF


      IF (init_mode == MODE_COSMODE) THEN
        CALL read_3d_1time(stream_id, on_cells, 'W', fill_array=initicon(jg)%atm_in%w_ifc)
      ELSE
        ! Note: in this case, input vertical velocity is in fact omega (Pa/s)
        CALL read_3d_1time(stream_id, on_cells, 'W', fill_array=initicon(jg)%atm_in%omega)
      ENDIF

      CALL read_3d_1time(stream_id, on_cells, 'QV', fill_array=initicon(jg)%atm_in%qv)
      CALL read_3d_1time(stream_id, on_cells, 'QC', fill_array=initicon(jg)%atm_in%qc)
      CALL read_3d_1time(stream_id, on_cells, 'QI', fill_array=initicon(jg)%atm_in%qi)

      IF (lread_qr) THEN
        CALL read_3d_1time(stream_id, on_cells, 'QR', fill_array=initicon(jg)%atm_in%qr)
      ELSE
        initicon(jg)%atm_in%qr(:,:,:)=0._wp
      ENDIF

      IF (lread_qs) THEN
        CALL read_3d_1time(stream_id, on_cells, 'QS', fill_array=initicon(jg)%atm_in%qs)
      ELSE
        initicon(jg)%atm_in%qs(:,:,:)=0._wp
      ENDIF

      IF (psvar_ndims==2)THEN
        CALL read_2d_1time(stream_id, on_cells, TRIM(psvar), &
          &                     fill_array=initicon(jg)%atm_in%psfc)
      ELSEIF(psvar_ndims==3)THEN
        CALL read_2d_1lev_1time(stream_id, on_cells, TRIM(psvar), &
          &                     fill_array=initicon(jg)%atm_in%psfc)
      ELSE
        CALL finish(TRIM(routine),'surface pressure var '//TRIM(psvar)//' dimension mismatch')
      END IF

      IF (geopvar_ndims==2)THEN
        CALL read_2d_1time(stream_id, on_cells, TRIM(geop_ml_var), &
          &                     fill_array=initicon(jg)%atm_in%phi_sfc)
      ELSEIF(geopvar_ndims==3)THEN
        CALL read_2d_1lev_1time(stream_id, on_cells, TRIM(geop_ml_var), &
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

        CALL read_3d_1time(stream_id, on_cells, 'HHL', fill_array=initicon(jg)%atm_in%z3d_ifc)
        CALL read_3d_1time(stream_id, on_cells, 'P', fill_array=initicon(jg)%atm_in%pres)

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
      ifs2icon_file(jg) = generate_filename(ifs2icon_filename, getModelBaseDir(), &
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
        CALL read_2d_1lev_1time(stream_id, on_cells, TRIM(geop_sfc_var), &
          &                     fill_array=initicon(jg)%sfc_in%phi)
      ELSE IF (geop_sfc_var_ndims == 2) THEN
        CALL read_2d_1time(stream_id, on_cells, TRIM(geop_sfc_var), &
          &                     fill_array=initicon(jg)%sfc_in%phi)
      ELSE
        CALL finish(TRIM(routine),"geop_sfc_var: Dimension mismatch")
      ENDIF

      CALL read_2d_1time(stream_id, on_cells, 'SKT', &
        &                fill_array=initicon(jg)%sfc_in%tskin)

      IF ( l_sst_in) THEN
        CALL read_2d_1time(stream_id, on_cells, 'SST', &
          &                fill_array=initicon(jg)%sfc_in%sst)
      ELSE
       initicon(jg)%sfc_in%sst(:,:)=0.0_wp
      END IF

      CALL read_2d_1time(stream_id, on_cells, 'T_SNOW', &
        &                fill_array=initicon(jg)%sfc_in%tsnow)
      CALL read_2d_1time(stream_id, on_cells, TRIM(alb_snow_var), &
        &                fill_array=initicon(jg)%sfc_in%snowalb)
      CALL read_2d_1time(stream_id, on_cells, 'W_SNOW', &
        &                fill_array=initicon(jg)%sfc_in%snowweq)
      CALL read_2d_1time(stream_id, on_cells, 'RHO_SNOW', &
        &                fill_array=initicon(jg)%sfc_in%snowdens)
      CALL read_2d_1time(stream_id, on_cells, 'W_I', &
        &                fill_array=initicon(jg)%sfc_in%skinres)
      CALL read_2d_1time(stream_id, on_cells, 'LSM', &
        &                fill_array=initicon(jg)%sfc_in%ls_mask)
      CALL read_2d_1time(stream_id, on_cells, 'CI', &
        &                fill_array=initicon(jg)%sfc_in%seaice)
      CALL read_2d_1lev_1time(stream_id, on_cells, 'STL1', &
        &                     fill_array=initicon(jg)%sfc_in%tsoil(:,:,1))
      CALL read_2d_1lev_1time(stream_id, on_cells, 'STL2', &
        &                     fill_array=initicon(jg)%sfc_in%tsoil(:,:,2))
      CALL read_2d_1lev_1time(stream_id, on_cells, 'STL3', &
        &                     fill_array=initicon(jg)%sfc_in%tsoil(:,:,3))
      CALL read_2d_1lev_1time(stream_id, on_cells, 'STL4', &
        &                     fill_array=initicon(jg)%sfc_in%tsoil(:,:,4))
      CALL read_2d_1lev_1time(stream_id, on_cells, 'SMIL1', &
        &                     fill_array=initicon(jg)%sfc_in%wsoil(:,:,1))
      CALL read_2d_1lev_1time(stream_id, on_cells, 'SMIL2', &
        &                     fill_array=initicon(jg)%sfc_in%wsoil(:,:,2))
      CALL read_2d_1lev_1time(stream_id, on_cells, 'SMIL3', &
        &                     fill_array=initicon(jg)%sfc_in%wsoil(:,:,3))
      CALL read_2d_1lev_1time(stream_id, on_cells, 'SMIL4', &
        &                     fill_array=initicon(jg)%sfc_in%wsoil(:,:,4))

      ! close file
      !
      IF(lread_process) CALL nf(nf_close(ncid), routine)
      CALL closeFile(stream_id)


    ENDDO ! loop over model domains

  END SUBROUTINE read_extana_sfc


  !>
  !! Request the DWD first guess variables (atmosphere only)
  SUBROUTINE request_dwdfg_atm(requestList)
    CLASS(t_InputRequestList), POINTER, INTENT(INOUT) :: requestList

    CALL requestList%requestMultiple((/ 'theta_v', &
                                      & 'rho    ', &
                                      & 'w      ', &
                                      & 'tke    ', &
                                      & 'qv     ', &
                                      & 'qc     ', &
                                      & 'qi     ', &
                                      & 'vn     ' /))
    IF ( iqr /= 0 ) CALL requestList%request('qr')
    IF ( iqs /= 0 ) CALL requestList%request('qs')
    IF (lvert_remap_fg) CALL requestList%request('z_ifc')
  END SUBROUTINE request_dwdfg_atm


  !>
  !! Fetch the DWD first guess from the request list (atmosphere only)
  !! First guess (FG) is read for theta_v, rho, vn, w, tke,
  !! whereas DA output is read for T, p, u, v,
  !! qv, qc, qi, qr, qs.
  SUBROUTINE fetch_dwdfg_atm(requestList, p_patch, p_nh_state, initicon, inputInstructions)
    CLASS(t_InputRequestList), POINTER, INTENT(INOUT) :: requestList
    TYPE(t_patch), INTENT(INOUT) :: p_patch(:)
    TYPE(t_nh_state), INTENT(INOUT), TARGET :: p_nh_state(:)
    TYPE(t_initicon_state), INTENT(INOUT), TARGET :: initicon(:)
    TYPE(t_readInstructionListPtr), INTENT(INOUT) :: inputInstructions(:)

    CHARACTER(*), PARAMETER :: routine = modname//':fetch_dwdfg_atm'
    INTEGER jg
    TYPE(t_nh_prog), POINTER :: prognosticFields
    REAL(wp), POINTER :: my_ptr3d(:,:,:)

    DO jg = 1, n_dom
        IF(p_patch(jg)%ldom_active) THEN
            ! save some paperwork
            prognosticFields => p_nh_state(jg)%prog(nnow(jg))

            ! request the first guess fields (atmosphere only)
            CALL requestList%fetchRequired3d('theta_v', trivial_tileId, jg, prognosticFields%theta_v)
            CALL requestList%fetchRequired3d('rho', trivial_tileId, jg, prognosticFields%rho)
            CALL requestList%fetchRequired3d('w', trivial_tileId, jg, prognosticFields%w)
            CALL requestList%fetchRequired3d('tke', trivial_tileId, jg, prognosticFields%tke)

            ! Only needed for FG-only runs; usually read from ANA
            my_ptr3d => prognosticFields%tracer(:,:,:,iqv)
            IF(inputInstructions(jg)%ptr%wantVarFg('qv')) &
            &   CALL inputInstructions(jg)%ptr%handleErrorFg(requestList%fetch3d('qv', trivial_tileId, jg, my_ptr3d), &
            &                                                'qv', routine)

            my_ptr3d => prognosticFields%tracer(:,:,:,iqc)
            IF(inputInstructions(jg)%ptr%wantVarFg('qc')) &
            &   CALL inputInstructions(jg)%ptr%handleErrorFg(requestList%fetch3d('qc', trivial_tileId, jg, my_ptr3d), &
            &                                                'qc', routine)

            my_ptr3d => prognosticFields%tracer(:,:,:,iqi)
            IF(inputInstructions(jg)%ptr%wantVarFg('qi')) &
            &   CALL inputInstructions(jg)%ptr%handleErrorFg(requestList%fetch3d('qi', trivial_tileId, jg, my_ptr3d), &
            &                                                'qi', routine)

            IF ( iqr /= 0 ) THEN
              my_ptr3d => prognosticFields%tracer(:,:,:,iqr)
              IF(inputInstructions(jg)%ptr%wantVarFg('qr')) &
              & CALL inputInstructions(jg)%ptr%handleErrorFg(requestList%fetch3d('qr', trivial_tileId, jg, my_ptr3d), &
              &                                              'qr', routine)
            END IF

            IF ( iqs /= 0 ) THEN
              my_ptr3d => prognosticFields%tracer(:,:,:,iqs)
              IF(inputInstructions(jg)%ptr%wantVarFg('qs')) &
              & CALL inputInstructions(jg)%ptr%handleErrorFg(requestList%fetch3d('qs', trivial_tileId, jg, my_ptr3d), &
              &                                              'qs', routine)
            END IF

            IF (lvert_remap_fg) THEN
                CALL allocate_extana_atm(jg, p_patch(jg)%nblks_c, p_patch(jg)%nblks_e, initicon)
                CALL requestList%fetchRequired3D('z_ifc', trivial_tileId, jg, initicon(jg)%atm_in%z3d_ifc)
            END IF

            CALL requestList%fetchRequired3d('vn', trivial_tileId, jg, prognosticFields%vn)
        END IF
    END DO

  END SUBROUTINE fetch_dwdfg_atm

  !>
  !! Request the DWD first guess variables (atmosphere ONLY)
  SUBROUTINE request_dwdfg_atm_ii(requestList)
    CLASS(t_InputRequestList), POINTER, INTENT(INOUT) :: requestList

    CALL requestList%requestMultiple((/ 'z_ifc  ', &
                                      & 'theta_v', &
                                      & 'rho    ', &
                                      & 'w      ', &
                                      & 'tke    ', &
                                      & 'qv     ', &
                                      & 'qc     ', &
                                      & 'qi     ', &
                                      & 'vn     ' /))
    IF ( iqr /= 0 ) CALL requestList%request('qr')
    IF ( iqs /= 0 ) CALL requestList%request('qs')
  END SUBROUTINE request_dwdfg_atm_ii

  !>
  !! Fetch DWD first guess from the request list (atmosphere only) and store to initicon input state
  !! First guess (FG) is read for z_ifc, theta_v, rho, vn, w, tke,
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
  SUBROUTINE fetch_dwdfg_atm_ii(requestList, p_patch, initicon)
    CLASS(t_InputRequestList), POINTER, INTENT(INOUT) :: requestList
    TYPE(t_patch), INTENT(IN) :: p_patch(:)
    TYPE(t_initicon_state), INTENT(INOUT), TARGET :: initicon(:)

    CHARACTER(len=*), PARAMETER :: routine = modname//':fetch_dwdfg_atm_ii'
    INTEGER :: jg, jb, jk, jc, i_endidx, nlev_in
    REAL(wp) :: tempv, exner
    REAL(dp), POINTER :: levelValues(:)

    !-------------------------------------------------------------------------

    DO jg = 1, n_dom
        IF(p_patch(jg)%ldom_active) THEN  ! Skip reading the atmospheric input data if a model domain is not active at initial time

            ! determine number of HALF LEVELS of generalized Z-AXIS
            levelValues => requestList%getLevels('z_ifc', jg)
            IF(.NOT.ASSOCIATED(levelValues)) CALL finish(routine, "no DATA found for domain "//TRIM(int2string(jg))//" of &
                                                                  &required variable 'theta_v'")
            nlev_in = SIZE(levelValues, 1) - 1
            nlevatm_in(jg) = nlev_in

            CALL allocate_extana_atm(jg, p_patch(jg)%nblks_c, p_patch(jg)%nblks_e, initicon)

            ! start reading first guess (atmosphere only)
            CALL requestList%fetchRequired3d('z_ifc', trivial_tileId, jg, initicon(jg)%atm_in%z3d_ifc)
            CALL requestList%fetchRequired3d('theta_v', trivial_tileId, jg, initicon(jg)%atm_in%theta_v)
            CALL requestList%fetchRequired3d('rho', trivial_tileId, jg, initicon(jg)%atm_in%rho)
            CALL requestList%fetchRequired3d('w', trivial_tileId, jg, initicon(jg)%atm_in%w_ifc)
            CALL requestList%fetchRequired3d('tke', trivial_tileId, jg, initicon(jg)%atm_in%tke_ifc)

            CALL requestList%fetchRequired3d('qv', trivial_tileId, jg, initicon(jg)%atm_in%qv)
            CALL requestList%fetchRequired3d('qc', trivial_tileId, jg, initicon(jg)%atm_in%qc)
            CALL requestList%fetchRequired3d('qi', trivial_tileId, jg, initicon(jg)%atm_in%qi)
            IF ( iqr /= 0 ) THEN
            CALL requestList%fetchRequired3d('qr', trivial_tileId, jg, initicon(jg)%atm_in%qr)
            ELSE
            initicon(jg)%atm_in%qr(:,:,:) = 0._wp
            END IF

            IF ( iqs /= 0 ) THEN
            CALL requestList%fetchRequired3d('qs', trivial_tileId, jg, initicon(jg)%atm_in%qs)
            ELSE
            initicon(jg)%atm_in%qs(:,:,:) = 0._wp
            END IF

            CALL requestList%fetchRequired3d('vn', trivial_tileId, jg, initicon(jg)%atm_in%vn)

            ! Interpolate half level variables from interface levels to main levels, and convert thermodynamic variables
            ! into temperature and pressure as expected by the vertical interpolation routine  
!$OMP PARALLEL
!$OMP DO PRIVATE (jk,jc,jb,i_endidx,tempv,exner) ICON_OMP_DEFAULT_SCHEDULE
            DO jb = 1,p_patch(jg)%nblks_c

                IF (jb /= p_patch(jg)%nblks_c) THEN
                    i_endidx = nproma
                ELSE
                    i_endidx = p_patch(jg)%npromz_c
                END IF

                DO jk = 1, nlevatm_in(jg)
                    DO jc = 1, i_endidx

                        !XXX: This code IS a waste of resources:
                        !       1. It IS memory bound.
                        !       2. Like all other places that USE the fields z3d_ifc, w_ifc, AND tke_ifc,
                        !          it uses them purely as temporary buffers to cache the input to the interpolation.
                        !     This could be avoided by moving these calculations directly after the fetch calls above (better cache USE),
                        !     scraping the *_ifc fields, AND using a temporary buffer instead (which can be reused for the next variable, smaller memory footprint).
                        !     However, since the fields are also used IN mo_initicon_utils, mo_async_latbc_utils, AND mo_sync_latbc IN a similar pattern,
                        !     all these places would need to be rewritten to scrape the *_ifc fields.
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
                            (initicon(jg)%atm_in%qc(jc,jk,jb) + initicon(jg)%atm_in%qi(jc,jk,jb) +                      &
                             initicon(jg)%atm_in%qr(jc,jk,jb) + initicon(jg)%atm_in%qs(jc,jk,jb)) )

                    END DO
                END DO
            END DO
!$OMP END DO
!$OMP END PARALLEL
        END IF
    ENDDO

    lread_vn = .TRUE.  ! Tell the vertical interpolation routine that vn needs to be processed
  END SUBROUTINE fetch_dwdfg_atm_ii


  !>
  !! Request the DA-analysis variables (atmosphere ONLY)
  SUBROUTINE request_dwdana_atm(requestList)
    CLASS(t_InputRequestList), POINTER, INTENT(INOUT) :: requestList

    CALL requestList%requestMultiple((/ 'temp', &
                                      & 'pres', &
                                      & 'u   ', &
                                      & 'v   ' /))

    CALL requestList%requestMultiple((/ 'qv', 'qc', 'qi' /))
    IF ( iqr /= 0 ) CALL requestList%request('qr')
    IF ( iqs /= 0 ) CALL requestList%request('qs')
  END SUBROUTINE request_dwdana_atm

  !>
  !! Fetch DA-analysis DATA from the request list (atmosphere only)
  !!
  !! Depending on the initialization mode, either full fields or increments
  !! are read (atmosphere only):
  !! MODE_DWDANA: The following full fields are read, if available
  !!              u, v, t, p, qv
  !! MODE_IAO_OLD: 
  !! MODE_IAU:
  !!
  !! @par Revision History
  !! Initial version by Daniel Reinert, DWD(2012-12-18)
  !! Modifications for GRIB2 : F. Prill, DWD (2013-02-19)
  !! Modifications by Daniel Reinert, DWD (2014-01-27)
  !! - split off reading of FG fields
  !!
  SUBROUTINE fetch_dwdana_atm(requestList, p_patch, p_nh_state, initicon, inputInstructions)
    CLASS(t_InputRequestList), POINTER, INTENT(INOUT) :: requestList
    TYPE(t_patch), INTENT(IN) :: p_patch(:)
    TYPE(t_nh_state), INTENT(INOUT) :: p_nh_state(:)
    TYPE(t_initicon_state), INTENT(INOUT), TARGET :: initicon(:)
    TYPE(t_readInstructionListPtr), INTENT(INOUT) :: inputInstructions(:)

    CHARACTER(LEN = *), PARAMETER :: routine = modname//':fetch_dwdana_atm'
    INTEGER :: jg
    TYPE(t_pi_atm), POINTER :: my_ptr
    REAL(wp), POINTER :: my_ptr3d(:,:,:)

    DO jg = 1, n_dom
        IF(p_patch(jg)%ldom_active) THEN  ! Skip reading the atmospheric input data if a model domain is not active at initial time
            ! Depending on the initialization mode chosen (incremental vs. non-incremental)
            ! input fields are stored in different locations.
            IF ( ANY((/MODE_IAU,MODE_IAU_OLD/) == init_mode) ) THEN
                ! Skip this domain if its interpolated from its parent in process_input_dwdana_atm()
                IF (lp2cintp_incr(jg)) CYCLE
                my_ptr => initicon(jg)%atm_inc
            ELSE
                my_ptr => initicon(jg)%atm
            ENDIF

            ! start reading DA output (atmosphere only)
            ! The dynamical variables temp, pres, u and v, which need further processing,
            ! are either stored in initicon(jg)%atm or initicon(jg)%atm_inc, depending on whether
            ! IAU is used or not. The moisture variables, which can be taken over directly from
            ! the Analysis, are written to the NH prognostic state
            IF(inputInstructions(jg)%ptr%wantVarAna('temp')) &
            &   CALL inputInstructions(jg)%ptr%handleErrorAna(requestList%fetch3d('temp', trivial_tileId, jg, my_ptr%temp), &
            &                                                 'temp', routine)
            IF(inputInstructions(jg)%ptr%wantVarAna('pres')) &
            &   CALL inputInstructions(jg)%ptr%handleErrorAna(requestList%fetch3d('pres', trivial_tileId, jg, my_ptr%pres), &
            &                                                 'pres', routine)
            IF(inputInstructions(jg)%ptr%wantVarAna('u')) &
            &   CALL inputInstructions(jg)%ptr%handleErrorAna(requestList%fetch3d('u', trivial_tileId, jg, my_ptr%u), &
            &                                                 'u', routine)
            IF(inputInstructions(jg)%ptr%wantVarAna('v')) &
            &   CALL inputInstructions(jg)%ptr%handleErrorAna(requestList%fetch3d('v', trivial_tileId, jg, my_ptr%v), &
            &                                                 'v', routine)

            IF ( ANY((/MODE_IAU,MODE_IAU_OLD/) == init_mode) ) THEN
                IF(inputInstructions(jg)%ptr%wantVarAna('qv')) &
                &   CALL inputInstructions(jg)%ptr%handleErrorAna(requestList%fetch3d('qv', trivial_tileId, jg, my_ptr%qv), &
                &                                                 'qv', routine)
            ELSE
                my_ptr3d => p_nh_state(jg)%prog(nnow(jg))%tracer(:,:,:,iqv)
                IF(inputInstructions(jg)%ptr%wantVarAna('qv')) &
                &   CALL inputInstructions(jg)%ptr%handleErrorAna(requestList%fetch3d('qv', trivial_tileId, jg, my_ptr3d), &
                &                                                 'qv', routine)
            ENDIF

            ! For the time being, these are identical to qc, qi, qr, AND qs from FG => usually read from FG
            my_ptr3d => p_nh_state(jg)%prog(nnow(jg))%tracer(:,:,:,iqc)
            IF(inputInstructions(jg)%ptr%wantVarAna('qc')) &
            &   CALL inputInstructions(jg)%ptr%handleErrorAna(requestList%fetch3d('qc', trivial_tileId, jg, my_ptr3d), &
            &                                                 'qc', routine)
            my_ptr3d => p_nh_state(jg)%prog(nnow(jg))%tracer(:,:,:,iqi)
            IF(inputInstructions(jg)%ptr%wantVarAna('qi')) &
            &   CALL inputInstructions(jg)%ptr%handleErrorAna(requestList%fetch3d('qi', trivial_tileId, jg, my_ptr3d), &
            &                                                 'qi', routine)
            IF ( iqr /= 0 ) THEN
                my_ptr3d => p_nh_state(jg)%prog(nnow(jg))%tracer(:,:,:,iqr)
                IF(inputInstructions(jg)%ptr%wantVarAna('qr')) &
                &   CALL inputInstructions(jg)%ptr%handleErrorAna(requestList%fetch3d('qr', trivial_tileId, jg, my_ptr3d), &
                &                                                 'qr', routine)
            END IF
            IF ( iqs /= 0 ) THEN
                my_ptr3d => p_nh_state(jg)%prog(nnow(jg))%tracer(:,:,:,iqs)
                IF(inputInstructions(jg)%ptr%wantVarAna('qs')) &
                &   CALL inputInstructions(jg)%ptr%handleErrorAna(requestList%fetch3d('qs', trivial_tileId, jg, my_ptr3d), &
                &                                                 'qs', routine)
            END IF
        END IF
    ENDDO ! loop over model domains

  END SUBROUTINE fetch_dwdana_atm

  !XXX: Cannot be moved into fetch_dwdana_atm() since it needs input from fetch_dwdana_sfc()
  SUBROUTINE process_input_dwdana_atm (p_patch, initicon)
    TYPE(t_patch),          INTENT(IN)    :: p_patch(:)
    TYPE(t_initicon_state), INTENT(INOUT), TARGET :: initicon(:)

    INTEGER :: jg

    DO jg = 1, n_dom
      IF(p_patch(jg)%ldom_active) THEN
        IF (ANY((/MODE_IAU,MODE_IAU_OLD/) == init_mode)) THEN
          IF (lp2cintp_incr(jg)) THEN
            ! Perform parent-to-child interpolation of atmospheric DA increments
            CALL interpolate_increments(initicon, p_patch(jg)%parent_id, jg)
          END IF
        END IF
      END IF
    END DO
  END SUBROUTINE process_input_dwdana_atm


  !>
  !! Request the DWD first guess variables (land/surface only)
  SUBROUTINE request_dwdfg_sfc(requestList)
    CLASS(t_InputRequestList), POINTER, INTENT(INOUT) :: requestList

    CALL requestList%requestMultiple((/ 't_ice      ', &
                                      & 'h_ice      ', &
                                      & 't_g        ', &
                                      & 'qv_s       ', &
                                      & 'freshsnow  ', &
                                      & 'snowfrac_lc', &
                                      & 'w_snow     ', &
                                      & 'w_i        ', &
                                      & 'h_snow     ', &
                                      & 't_snow     ', &
                                      & 'rho_snow   ', &
                                      & 't_mnw_lk   ', &
                                      & 't_wml_lk   ', &
                                      & 'h_ml_lk    ', &
                                      & 't_bot_lk   ', &
                                      & 'c_t_lk     ', &
                                      & 't_b1_lk    ', &
                                      & 'h_b1_lk    ', &
                                      & 'aer_ss     ', &
                                      & 'aer_or     ', &
                                      & 'aer_bc     ', &
                                      & 'aer_su     ', &
                                      & 'aer_du     '/))

    IF (lmulti_snow) CALL requestList%requestMultiple((/ 't_snow_mult  ', &
                                                       & 'rho_snow_mult', &
                                                       & 'wtot_snow    ', &
                                                       & 'wliq_snow    ', &
                                                       & 'dzh_snow     ' /))    ! multi layer snow fields

    SELECT CASE(init_mode)
        CASE(MODE_COMBINED)
            CALL requestList%requestMultiple((/ 'smi      ', &
                                              & 'fr_seaice' /))
        CASE(MODE_COSMODE)
            CALL requestList%requestMultiple((/ 'smi' /))
        CASE DEFAULT
            CALL requestList%requestMultiple((/ 'w_so     ', &
                                              & 'fr_seaice' /))
    END SELECT

    CALL requestList%requestMultiple((/ 'w_so_ice', &
                                      & 't_so    ' /))
    CALL requestList%request('gz0')
  END SUBROUTINE request_dwdfg_sfc

  !>
  !! Fetch DWD first guess DATA from request list (land/surface only)
  SUBROUTINE fetch_dwdfg_sfc(requestList, p_patch, prm_diag, p_lnd_state, inputInstructions)
    CLASS(t_InputRequestList), POINTER, INTENT(INOUT) :: requestList
    TYPE(t_patch),             INTENT(IN)    :: p_patch(:)
    TYPE(t_nwp_phy_diag),      INTENT(INOUT) :: prm_diag(:)
    TYPE(t_lnd_state), TARGET, INTENT(INOUT) :: p_lnd_state(:)
    TYPE(t_readInstructionListPtr), INTENT(INOUT) :: inputInstructions(:)

    INTEGER :: jg, jt, jb, jc, i_endidx
    REAL(wp), POINTER :: my_ptr2d(:,:)
    REAL(wp), POINTER :: my_ptr3d(:,:,:)
    TYPE(t_lnd_prog), POINTER :: lnd_prog
    TYPE(t_lnd_diag), POINTER :: lnd_diag
    TYPE(t_wtr_prog), POINTER :: wtr_prog
    CHARACTER(LEN=VARNAME_LEN), POINTER :: checkgrp(:)
    TYPE(t_inputParameters) :: parameters

    CHARACTER(len=*), PARAMETER :: routine = modname//':fetch_dwdfg_sfc'

    DO jg = 1, n_dom
        IF(p_patch(jg)%ldom_active) THEN ! Skip reading the atmospheric input data if a model domain is not active at initial time
            lnd_prog => p_lnd_state(jg)%prog_lnd(nnow_rcf(jg))
            lnd_diag => p_lnd_state(jg)%diag_lnd
            wtr_prog => p_lnd_state(jg)%prog_wtr(nnow_rcf(jg))

            ! COSMO-DE does not provide sea ice field. In that case set fr_seaice to 0
            IF (init_mode /= MODE_COSMODE) THEN
                CALL fetchSurfaceWrapper('fr_seaice', jg, lnd_diag%fr_seaice)
            ELSE
!$OMP PARALLEL
            CALL init(lnd_diag%fr_seaice(:,:))
!$OMP END PARALLEL
            ENDIF ! init_mode /= MODE_COSMODE

            ! sea-ice related fields
            CALL fetchSurfaceWrapper('t_ice', jg, wtr_prog%t_ice)
            CALL fetchSurfaceWrapper('h_ice', jg, wtr_prog%h_ice)

            !These two fields are required for the processing step below, AND they are NOT initialized before this SUBROUTINE IS called, so they are fetched as required.
            !This diverges from the code that I found which READ them conditionally.
            CALL fetchRequiredTiledSurfaceWrapper('t_g', jg, ntiles_total + ntiles_water, lnd_prog%t_g_t)
            CALL fetchRequiredTiledSurfaceWrapper('qv_s', jg, ntiles_total + ntiles_water, lnd_diag%qv_s_t)

            CALL fetchTiledSurfaceWrapper('freshsnow', jg, ntiles_total, lnd_diag%freshsnow_t)
            CALL fetchTiledSurfaceWrapper('snowfrac_lc', jg, ntiles_total, lnd_diag%snowfrac_lc_t)
            CALL fetchTiledSurfaceWrapper('w_snow', jg, ntiles_total, lnd_prog%w_snow_t)
            CALL fetchTiledSurfaceWrapper('w_i', jg, ntiles_total, lnd_prog%w_i_t)
            CALL fetchTiledSurfaceWrapper('h_snow', jg, ntiles_total, lnd_diag%h_snow_t)
            CALL fetchTiledSurfaceWrapper('t_snow', jg, ntiles_total, lnd_prog%t_snow_t)
            CALL fetchTiledSurfaceWrapper('rho_snow', jg, ntiles_total, lnd_prog%rho_snow_t)

            IF (lmulti_snow) THEN
                ! multi layer snow fields
                CALL fetchTiled3dWrapper('t_snow_mult', jg, ntiles_total, lnd_prog%t_snow_mult_t)
                CALL fetchTiled3dWrapper('rho_snow_mult', jg, ntiles_total, lnd_prog%rho_snow_mult_t)
                CALL fetchTiled3dWrapper('wtot_snow', jg, ntiles_total, lnd_prog%wtot_snow_t)
                CALL fetchTiled3dWrapper('wliq_snow', jg, ntiles_total, lnd_prog%wliq_snow_t)
                CALL fetchTiled3dWrapper('dzh_snow', jg, ntiles_total, lnd_prog%dzh_snow_t)
            END IF ! lmulti_snow


            ! multi layer fields
            !
            ! Note that either w_so OR smi is written to w_so_t. Which one is required depends
            ! on the initialization mode. Checking grp_vars_fg takes care of this. In case
            ! that smi is read, it is lateron converted to w_so (see smi_to_wsoil)
            SELECT CASE(init_mode)
                CASE(MODE_COMBINED, MODE_COSMODE)
                    CALL fetchTiled3dWrapper('smi', jg, ntiles_total, lnd_prog%w_so_t)
                CASE DEFAULT
                    CALL fetchTiled3dWrapper('w_so', jg, ntiles_total, lnd_prog%w_so_t)
                    CALL fetchTiled3dWrapper('w_so_ice', jg, ntiles_total, lnd_prog%w_so_ice_t) ! w_so_ice is re-diagnosed in terra_multlay_init
            END SELECT

            CALL fetchTiled3dWrapper('t_so', jg, ntiles_total, lnd_prog%t_so_t)

            ! Skipped in MODE_COMBINED and in MODE_COSMODE (i.e. when starting from GME soil)
            ! Instead z0 is re-initialized (see mo_nwp_phy_init)
            CALL fetchSurfaceWrapper('gz0', jg, prm_diag(jg)%gz0)

            ! first guess for fresh water lake fields
            IF(ASSOCIATED(wtr_prog%t_mnw_lk)) CALL fetchSurfaceWrapper('t_mnw_lk', jg, wtr_prog%t_mnw_lk)
            IF(ASSOCIATED(wtr_prog%t_wml_lk)) CALL fetchSurfaceWrapper('t_wml_lk', jg, wtr_prog%t_wml_lk)
            IF(ASSOCIATED(wtr_prog%h_ml_lk)) CALL fetchSurfaceWrapper('h_ml_lk', jg, wtr_prog%h_ml_lk)
            IF(ASSOCIATED(wtr_prog%t_bot_lk)) CALL fetchSurfaceWrapper('t_bot_lk', jg, wtr_prog%t_bot_lk)
            IF(ASSOCIATED(wtr_prog%c_t_lk)) CALL fetchSurfaceWrapper('c_t_lk', jg, wtr_prog%c_t_lk)
            IF(ASSOCIATED(wtr_prog%t_b1_lk)) CALL fetchSurfaceWrapper('t_b1_lk', jg, wtr_prog%t_b1_lk)
            IF(ASSOCIATED(wtr_prog%h_b1_lk)) CALL fetchSurfaceWrapper('h_b1_lk', jg, wtr_prog%h_b1_lk)

            aerosol_fg_present(jg) = .TRUE.
            my_ptr2d => prm_diag(jg)%aerosol(:,iss,:)
            IF(.NOT. requestList%fetchSurface('aer_ss', trivial_tileId, jg, my_ptr2d)) aerosol_fg_present(jg) = .FALSE.
            my_ptr2d => prm_diag(jg)%aerosol(:,iorg,:)
            IF(.NOT. requestList%fetchSurface('aer_or', trivial_tileId, jg, my_ptr2d)) aerosol_fg_present(jg) = .FALSE.
            my_ptr2d => prm_diag(jg)%aerosol(:,ibc,:)
            IF(.NOT. requestList%fetchSurface('aer_bc', trivial_tileId, jg, my_ptr2d)) aerosol_fg_present(jg) = .FALSE.
            my_ptr2d => prm_diag(jg)%aerosol(:,iso4,:)
            IF(.NOT. requestList%fetchSurface('aer_su', trivial_tileId, jg, my_ptr2d)) aerosol_fg_present(jg) = .FALSE.
            my_ptr2d => prm_diag(jg)%aerosol(:,idu,:)
            IF(.NOT. requestList%fetchSurface('aer_du', trivial_tileId, jg, my_ptr2d)) aerosol_fg_present(jg) = .FALSE.
        END IF
    END DO

  CONTAINS

    ! Wrapper for requestList%fetchSurface()
    SUBROUTINE fetchSurfaceWrapper(varName, jg, field)
        CHARACTER(LEN = *), INTENT(IN) :: varName
        INTEGER, VALUE :: jg
        REAL(wp), INTENT(INOUT) :: field(:,:)

        INTEGER :: jt

        IF(inputInstructions(jg)%ptr%wantVarFg(varName)) THEN
            CALL inputInstructions(jg)%ptr%handleErrorFg(requestList%fetchSurface(varName, trivial_tileId, jg, field), &
            &                                            varName, routine)
        END IF
    END SUBROUTINE fetchSurfaceWrapper

    ! Wrapper for requestList%fetchTiledSurface() that falls back to reading copies of untiled input IF ltile_coldstart IS set.
    SUBROUTINE fetchTiledSurfaceWrapper(varName, jg, tileCount, field)
        CHARACTER(LEN = *), INTENT(IN) :: varName
        INTEGER, VALUE :: jg, tileCount
        REAL(wp), INTENT(INOUT) :: field(:,:,:)

        INTEGER :: jt

        IF(inputInstructions(jg)%ptr%wantVarFg(varName)) THEN
            IF(ltile_coldstart) THEN
                !Fake tiled input by copying the input field to all tiles.
                DO jt = 1, tileCount
                    CALL inputInstructions(jg)%ptr%handleErrorFg(requestList%fetchSurface(varName, trivial_tileId, jg, &
                    &                                                                     field(:,:,jt)), &
                    &                                            varName, routine)
                END DO
            ELSE
                !True tiled input.
                CALL inputInstructions(jg)%ptr%handleErrorFg(requestList%fetchTiledSurface(varName, jg, field), &
                &                                            varName, routine)
            END IF
        END IF
    END SUBROUTINE fetchTiledSurfaceWrapper

    ! Wrapper for requestList%fetchRequiredTiledSurface() that falls back to reading copies of untiled input IF ltile_coldstart IS set.
    SUBROUTINE fetchRequiredTiledSurfaceWrapper(varName, jg, tileCount, field)
        CHARACTER(LEN = *), INTENT(IN) :: varName
        INTEGER, VALUE :: jg, tileCount
        REAL(wp), INTENT(INOUT) :: field(:,:,:)

        INTEGER :: jt

        IF(ltile_coldstart) THEN
            !Fake tiled input by copying the input field to all tiles.
            DO jt = 1, tileCount
                CALL requestList%fetchRequiredSurface(varName, trivial_tileId, jg, field(:,:,jt))
            END DO
        ELSE
            !True tiled input.
            CALL requestList%fetchRequiredTiledSurface(varName, jg, field)
        END IF
    END SUBROUTINE fetchRequiredTiledSurfaceWrapper

    ! Wrapper for requestList%fetchTiled3d() that falls back to reading copies of untiled input IF ltile_coldstart IS set.
    SUBROUTINE fetchTiled3dWrapper(varName, jg, tileCount, field)
        CHARACTER(LEN = *), INTENT(IN) :: varName
        INTEGER, VALUE :: jg, tileCount
        REAL(wp), INTENT(INOUT) :: field(:,:,:,:)

        INTEGER :: jt

        IF(inputInstructions(jg)%ptr%wantVarFg(varName)) THEN
            IF(ltile_coldstart) THEN
                !Fake tiled input by copying the input field to all tiles.
                DO jt = 1, tileCount
                    CALL inputInstructions(jg)%ptr%handleErrorFg(requestList%fetch3d(varName, trivial_tileId, jg, &
                    &                                                                field(:,:,:,jt)), &
                    &                                            varName, routine)
                END DO
            ELSE
                !True tiled input.
                CALL inputInstructions(jg)%ptr%handleErrorFg(requestList%fetchTiled3d(varName, jg, field), &
                &                                            varName, routine)
            END IF
        END IF
    END SUBROUTINE fetchTiled3dWrapper

  END SUBROUTINE fetch_dwdfg_sfc

  !>
  !! Request the DWD analysis variables (land/surface only)
  SUBROUTINE request_dwdana_sfc(requestList)
    CLASS(t_InputRequestList), POINTER, INTENT(INOUT) :: requestList

    CALL requestList%requestMultiple((/ 'fr_seaice', &
                                      & 't_ice    ', &
                                      & 'h_ice    ', &
                                      & 't_so     ', &
                                      & 'h_snow   ', &
                                      & 'w_snow   ', &
                                      & 'w_i      ', &
                                      & 't_snow   ', &
                                      & 'rho_snow ', &
                                      & 'freshsnow', &
                                      & 'w_so     '/))
  END SUBROUTINE request_dwdana_sfc

  !>
  !! processing of DWD first-guess DATA
  SUBROUTINE process_input_dwdfg_sfc(p_patch, p_lnd_state)
    TYPE(t_patch),             INTENT(IN)    :: p_patch(:)
    TYPE(t_lnd_state), TARGET, INTENT(INOUT) :: p_lnd_state(:)

    INTEGER :: jg, jt, jb, jc, i_endidx
    TYPE(t_lnd_prog), POINTER :: lnd_prog
    TYPE(t_lnd_diag), POINTER :: lnd_diag
    CHARACTER(LEN = *), PARAMETER :: routine = modname//":process_input_dwdfg_sfc"

    DO jg = 1, n_dom
        IF(p_patch(jg)%ldom_active) THEN
            lnd_prog => p_lnd_state(jg)%prog_lnd(nnow_rcf(jg))
            lnd_diag => p_lnd_state(jg)%diag_lnd
            DO jt=1, ntiles_total + ntiles_water
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
                        END IF

                        DO jc = 1, i_endidx
                            lnd_prog%t_g(jc,jb)  = lnd_prog%t_g_t(jc,jb,1)
                            lnd_diag%qv_s(jc,jb) = lnd_diag%qv_s_t(jc,jb,1)
                        END DO  ! jc
                    END DO  ! jb
!$OMP END PARALLEL DO
                END IF
            END DO

            ! Only required, when starting from GME or COSMO soil (i.e. MODE_COMBINED or MODE_COSMODE).
            ! SMI stored in w_so_t must be converted to w_so
            IF (ANY((/MODE_COMBINED,MODE_COSMODE/) == init_mode)) THEN
                DO jt=1, ntiles_total
                    CALL smi_to_wsoil(p_patch(jg), lnd_prog%w_so_t(:,:,:,jt))
                END DO
            END IF
        END IF
    END DO
  END SUBROUTINE process_input_dwdfg_sfc


  !>
  !! Fetch DWD analysis DATA from the request list (land/surface only)
  !!
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
  SUBROUTINE fetch_dwdana_sfc(requestList, p_patch, p_lnd_state, initicon, inputInstructions)
    CLASS(t_InputRequestList), POINTER, INTENT(INOUT) :: requestList
    TYPE(t_patch),             INTENT(IN)    :: p_patch(:)
    TYPE(t_lnd_state), TARGET, INTENT(INOUT) :: p_lnd_state(:)
    TYPE(t_initicon_state), INTENT(INOUT), TARGET :: initicon(:)
    TYPE(t_readInstructionListPtr), INTENT(INOUT) :: inputInstructions(:)

    INTEGER :: jg, jt
    INTEGER :: ngrp_vars_ana, filetype
    REAL(wp), POINTER :: my_ptr2d(:,:)
    REAL(wp), POINTER :: my_ptr3d(:,:,:)
    TYPE(t_lnd_prog), POINTER :: lnd_prog
    TYPE(t_lnd_diag), POINTER :: lnd_diag
    TYPE(t_wtr_prog), POINTER :: wtr_prog

    CHARACTER(LEN = *), PARAMETER :: routine = modname//':fetch_dwdana_sfc'

    !----------------------------------------!
    ! read in DWD analysis (surface)         !
    !----------------------------------------!

    DO jg = 1, n_dom
        IF(p_patch(jg)%ldom_active) THEN
            IF (ANY((/MODE_IAU, MODE_IAU_OLD /) == init_mode) .AND. lp2cintp_sfcana(jg)) CYCLE  ! Skip this domain if it will be interpolated from its parent domain in process_input_dwdana_sfc()

            lnd_prog => p_lnd_state(jg)%prog_lnd(nnow_rcf(jg))
            lnd_diag => p_lnd_state(jg)%diag_lnd
            wtr_prog => p_lnd_state(jg)%prog_wtr(nnow_rcf(jg))

            jt = 1  ! set tile-index explicitly

            ! sea-ice fraction
            IF(inputInstructions(jg)%ptr%wantVarAna('fr_seaice')) &
            &   CALL inputInstructions(jg)%ptr%handleErrorAna(requestList%fetchSurface('fr_seaice', trivial_tileId, jg, &
            &                                                                          lnd_diag%fr_seaice), &
            &                                                 'fr_seaice', routine)
            ! sea-ice temperature
            IF(inputInstructions(jg)%ptr%wantVarAna('t_ice')) &
            &   CALL inputInstructions(jg)%ptr%handleErrorAna(requestList%fetchSurface('t_ice', trivial_tileId, jg, &
            &                                                                          wtr_prog%t_ice), &
            &                                                 't_ice', routine)
            ! sea-ice height
            IF(inputInstructions(jg)%ptr%wantVarAna('h_ice')) &
            &   CALL inputInstructions(jg)%ptr%handleErrorAna(requestList%fetchSurface('h_ice', trivial_tileId, jg, &
            &                                                                          wtr_prog%h_ice), &
            &                                                 'h_ice', routine)

            ! T_SO(0). Note that the file may contain a 3D field, of which we ONLY fetch the level at 0.0.
            my_ptr2d => initicon(jg)%sfc%sst(:,:)
            IF(inputInstructions(jg)%ptr%wantVarAna('t_so')) &
            &   CALL inputInstructions(jg)%ptr%handleErrorAna(requestList%fetch2d('t_so', 0.0_wp, trivial_tileId, jg, my_ptr2d), &
            &                                                 't_so', routine)

            ! h_snow
            IF ( init_mode == MODE_IAU ) THEN
                my_ptr2d => initicon(jg)%sfc_inc%h_snow(:,:)
            ELSE
                my_ptr2d => lnd_diag%h_snow_t(:,:,jt)
            ENDIF
            IF(inputInstructions(jg)%ptr%wantVarAna('h_snow')) &
            &   CALL inputInstructions(jg)%ptr%handleErrorAna(requestList%fetchSurface('h_snow', trivial_tileId, jg, my_ptr2d), &
            &                                                 'h_snow', routine)

            ! w_snow
            my_ptr2d => lnd_prog%w_snow_t(:,:,jt)
            IF(inputInstructions(jg)%ptr%wantVarAna('w_snow')) &
            &   CALL inputInstructions(jg)%ptr%handleErrorAna(requestList%fetchSurface('w_snow', trivial_tileId, jg, my_ptr2d), &
            &                                                 'w_snow', routine)

            ! w_i
            my_ptr2d => lnd_prog%w_i_t(:,:,jt)
            IF(inputInstructions(jg)%ptr%wantVarAna('w_i')) &
            &   CALL inputInstructions(jg)%ptr%handleErrorAna(requestList%fetchSurface('w_i', trivial_tileId, jg, my_ptr2d), &
            &                                                 'w_i', routine)

            ! t_snow
            my_ptr2d => lnd_prog%t_snow_t(:,:,jt)
            IF(inputInstructions(jg)%ptr%wantVarAna('t_snow')) &
            &   CALL inputInstructions(jg)%ptr%handleErrorAna(requestList%fetchSurface('t_snow', trivial_tileId, jg, my_ptr2d), &
            &                                                 't_snow', routine)

            ! rho_snow
            my_ptr2d => lnd_prog%rho_snow_t(:,:,jt)
            IF(inputInstructions(jg)%ptr%wantVarAna('rho_snow')) &
            &   CALL inputInstructions(jg)%ptr%handleErrorAna(requestList%fetchSurface('rho_snow', trivial_tileId, jg, my_ptr2d), &
            &                                                 'rho_snow', routine)

            ! freshsnow
            IF ( init_mode == MODE_IAU ) THEN
                my_ptr2d => initicon(jg)%sfc_inc%freshsnow(:,:)
            ELSE
                my_ptr2d => lnd_diag%freshsnow_t(:,:,jt)
            ENDIF
            IF(inputInstructions(jg)%ptr%wantVarAna('freshsnow')) &
            &   CALL inputInstructions(jg)%ptr%handleErrorAna(requestList%fetchSurface('freshsnow', trivial_tileId, jg, &
            &                                                                          my_ptr2d), &
            &                                                 'freshsnow', routine)

            ! w_so
            IF ( (init_mode == MODE_IAU) .OR. (init_mode == MODE_IAU_OLD) ) THEN
                my_ptr3d => initicon(jg)%sfc_inc%w_so(:,:,:)
            ELSE
                my_ptr3d => lnd_prog%w_so_t(:,:,:,jt)
            END IF
            IF(inputInstructions(jg)%ptr%wantVarAna('w_so')) &
            &   CALL inputInstructions(jg)%ptr%handleErrorAna(requestList%fetch3d('w_so', trivial_tileId, jg, my_ptr3d), &
            &                                                 'w_so', routine)
        END IF
    END DO ! loop over model domains

  END SUBROUTINE fetch_dwdana_sfc


  ! Fill remaining tiles for snow variables if tile approach is used
  ! Only fields that are actually read from the snow analysis are copied; note that MODE_IAU is mandatory when using tiles
  SUBROUTINE process_input_dwdana_sfc (p_patch, p_lnd_state, initicon)
    TYPE(t_patch), INTENT(IN) :: p_patch(:)
    TYPE(t_lnd_state), TARGET, INTENT(INOUT) :: p_lnd_state(:)
    TYPE(t_initicon_state), INTENT(INOUT), TARGET :: initicon(:)

    INTEGER :: jg, jt, jb, jc, i_endidx
    TYPE(t_lnd_prog), POINTER :: lnd_prog
    TYPE(t_lnd_diag), POINTER :: lnd_diag

    jt = 1
    DO jg = 1, n_dom
      IF(p_patch(jg)%ldom_active .AND. ANY((/MODE_IAU, MODE_IAU_OLD /) == init_mode)) THEN
        IF (lp2cintp_sfcana(jg)) THEN
          ! Perform parent-to-child interpolation of surface fields read from the analysis
          CALL interpolate_sfcana(initicon, p_patch(jg)%parent_id, jg)

        ELSE IF (ntiles_total>1 .AND. init_mode == MODE_IAU_OLD) THEN
          ! MODE_IAU_OLD: H_SNOW, FRESHSNOW, W_SNOW and RHO_SNOW are read from analysis (full fields)
          ! Since only the first tile index is filled (see above), we fill (here) the remaining tiles
          ! if tile approach is used. Note that for MODE_IAU this copy is skipped on purpose,
          ! since for ltile_coldstart=.FALSE. tile information would be overwritten.
          lnd_prog => p_lnd_state(jg)%prog_lnd(nnow_rcf(jg))
          lnd_diag => p_lnd_state(jg)%diag_lnd

!$OMP PARALLEL DO PRIVATE(jb,jc,jt,i_endidx)
          DO jb = 1, p_patch(jg)%nblks_c
            IF (jb == p_patch(jg)%nblks_c) THEN
              i_endidx = p_patch(jg)%npromz_c
            ELSE
              i_endidx = nproma
            END IF

            DO jt = 2, ntiles_total
              DO jc = 1, i_endidx
                lnd_diag%freshsnow_t(jc,jb,jt) = lnd_diag%freshsnow_t(jc,jb,1)
                lnd_diag%h_snow_t(jc,jb,jt)    = lnd_diag%h_snow_t(jc,jb,1)
                lnd_prog%w_snow_t(jc,jb,jt)    = lnd_prog%w_snow_t(jc,jb,1)
                lnd_prog%rho_snow_t(jc,jb,jt)  = lnd_prog%rho_snow_t(jc,jb,1)
              END DO
            END DO
          END DO
!$OMP END PARALLEL DO
        END IF
      END IF
    END DO ! loop over model domains
  END SUBROUTINE process_input_dwdana_sfc

END MODULE mo_initicon_io

