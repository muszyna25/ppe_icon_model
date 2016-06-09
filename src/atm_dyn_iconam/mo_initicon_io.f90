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
  USE mo_ext_data_types,      ONLY: t_external_data
  USE mo_nonhydro_types,      ONLY: t_nh_state, t_nh_prog
  USE mo_nwp_phy_types,       ONLY: t_nwp_phy_diag
  USE mo_nwp_lnd_types,       ONLY: t_lnd_state, t_lnd_prog, t_lnd_diag, t_wtr_prog
  USE mo_initicon_types,      ONLY: t_initicon_state, t_pi_atm, alb_snow_var, geop_ml_var
  USE mo_input_instructions,  ONLY: t_readInstructionListPtr, kInputSourceNone, kInputSourceFg, &
    &                               kInputSourceAna, kInputSourceBoth, kStateFailedFetch
  USE mo_initicon_config,     ONLY: init_mode, nlevatm_in, l_sst_in, generate_filename, &
    &                               ifs2icon_filename, dwdfg_filename, dwdana_filename, &
    &                               nml_filetype => filetype, lread_vn,      &
    &                               lp2cintp_incr, lp2cintp_sfcana, ltile_coldstart,    &
    &                               lvert_remap_fg, aerosol_fg_present
  USE mo_nh_init_nest_utils,  ONLY: interpolate_increments, interpolate_sfcana
  USE mo_impl_constants,      ONLY: MAX_CHAR_LENGTH, max_dom,                           &
    &                               MODE_IAU, MODE_IAU_OLD, MODE_IFSANA, MODE_COMBINED, &
    &                               MODE_COSMODE, iss, iorg, ibc, iso4, idu, SUCCESS
  USE mo_exception,           ONLY: message, finish, message_text
  USE mo_grid_config,         ONLY: n_dom, nroot, l_limited_area
  USE mo_mpi,                 ONLY: p_io, p_bcast, p_comm_work,    &
    &                               p_comm_work_test, my_process_is_mpi_workroot
  USE mo_io_config,           ONLY: default_read_method
  USE mo_read_interface,      ONLY: t_stream_id, nf, openInputFile, closeFile, &
    &                               read_2d_1time, read_2d_1lev_1time, &
    &                               read_3d_1time, on_cells, on_edges
  USE mo_util_cdi,            ONLY: trivial_tileId
  USE mo_ifs_coord,           ONLY: alloc_vct, init_vct, vct, vct_a, vct_b
  USE mo_lnd_nwp_config,      ONLY: ntiles_total,  l2lay_rho_snow, &
    &                               ntiles_water, lmulti_snow, lsnowtile, &
    &                               isub_lake, llake
  USE mo_master_config,       ONLY: getModelBaseDir
  USE mo_nwp_sfc_interp,      ONLY: smi_to_wsoil
  USE mo_io_util,             ONLY: get_filetype
  USE mo_initicon_utils,      ONLY: allocate_extana_atm, allocate_extana_sfc
  USE mo_physical_constants,  ONLY: cpd, rd, cvd_o_rd, p0ref, vtmpc1, tmelt
  USE mo_fortran_tools,       ONLY: init
  USE mo_input_request_list,  ONLY: t_InputRequestList
  USE mo_util_string,         ONLY: int2string
  USE mo_atm_phy_nwp_config,  ONLY: iprog_aero



  ! High level overview of how `mo_initicon` reads input data
  ! =========================================================
  !
  ! Knobs that control what is read
  ! -------------------------------
  !
  !  1. variable groups
  !
  !  2. init_mode
  !
  !  3. further flags and the `ana_varlist` field from the namelists
  !
  !
  ! Steps of reading
  ! ----------------
  !
  !  1. A list with input instructions is generated.
  !     The code for this is found in `mo_input_instructions`.
  !
  !     This is where the variable groups enter the process, all further processing depends on the variable groups only indirectly via `t_ReadInstructions`.
  !     This is also the place where the `ana_varlist` namelist parameter is evaluated.
  !
  !     The input instructions serve a dual purpose (which is a source of confusion, unfortunately):
  !     They encode from which input file a read attempt for a variable should be made,
  !     but they also collect whether such an attempt has been made, and whether that attempt was successfull.
  !
  !  2. The input instructions are used to generate a `t_InputRequestList`.
  !     This is done by the `fileRequests()` method of `t_ReadInstructionList`, and is quite straight-forward:
  !     All variables that can be read from a file are requested from that file.
  !
  !  3. The files are read, and a file inventory output is produced in the process (mo_initicon: read_dwdfg() and read_dwdana()).
  !     Note that the file inventory contains only information about variables that have been requested!
  !
  !  4. The data is fetched from the `t_InputRequestList`.
  !     This happens in `fetch_dwd...()` routines in `mo_initicon_io`.
  !
  !     This is a complicated process that takes much more conditions into account then the code generating the `t_ReadInstructions` does.
  !     Especially, this is the place where all the different namelist flags are honored, and it's the place where optional reading is handled
  !     (variables that are read if they are present in the file, but which are not necessary to start the run).
  !
  !     The `fetch_dwd...()` routines generally first ask the `t_ReadInstructions` whether they should attempt to read a variable,
  !     then they try to fetch the data from the `t_InputRequestList`, and finally inform the `t_ReadInstructions` about the result of this operation.
  !     It is then the duty of `t_ReadInstructions` to check whether a failure to read is fatal or whether it can be compensated.
  !
  !     This ask-fetch-inform triple is usually encapsulated within the wrapper functions `fetch2d()` etc. within this file.
  !     Exceptions exist, like the all-or-nothing behavior that's implemented for the aerosol variables, but they are not frequent.
  !
  !  5. The `t_ReadInstructions` is asked to output a table describing what read attempts have been made
  !     and which data is actually used for the model run.
  !
  !     Note that this is about what read attempts have *really* been made (= there was an ask-fetch-inform triple for this variable),
  !     not what read attempts *should* have been made (= `t_ReadInstructions::wantVar()` would have returned `.true.` if it had been called).
  !
  !
  ! How to add variables
  ! --------------------
  !
  ! There are two places that need to be changed to add support for a new input variable:
  !
  !  1. The variable has to be added to the relevant variable groups, so that input instructions are generated for it, and so that it's requested from the right files.
  !
  !  2. Code has to be added to one of the `fetch_dwd...()` routines in this file to retrieve the data from the file(s)
  !     and to inform the `t_ReadInstructions` about the succes/failure to do so.
  !
  !     In most cases, this is simply done by calling the right wrapper routine (`fetch2d()` and friends).
  !     If the appropriate wrapper routine is missing, just copy-paste-modify one of the existing wrappers;
  !     I have only generated the ones that were actually needed to implement the current functionality, so the list of wrappers is incomplete.

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


  TYPE :: t_fetchParams
    TYPE(t_readInstructionListPtr), ALLOCATABLE :: inputInstructions(:)
    CLASS(t_InputRequestList), POINTER :: requestList
    CHARACTER(LEN = :), ALLOCATABLE :: routine
    LOGICAL :: isFg
  END TYPE t_fetchParams

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !>
  !! Some wrapper routines to fetch DATA from an InputRequestList that first ask the inputInstructions whether a READ attempt should be made,
  !! AND proceed to tell the inputInstructions whether that attempt was successful.
  !! In the CASE of tiled input, these also take care of checking the ltile_coldstart flag, falling back to reading AND copying untiled input DATA IF it's set.

  SUBROUTINE fetch2d(params, varName, level, jg, field)
    TYPE(t_fetchParams), INTENT(INOUT) :: params
    CHARACTER(LEN = *), INTENT(IN) :: varName
    REAL(dp), VALUE :: level
    INTEGER, VALUE :: jg
    REAL(wp), INTENT(INOUT) :: field(:,:)

    LOGICAL :: fetchResult

    IF(params%inputInstructions(jg)%ptr%wantVar(varName, params%isFg)) THEN
        fetchResult = params%requestList%fetch2d(varName, level, trivial_tileId, jg, field)
        CALL params%inputInstructions(jg)%ptr%handleError(fetchResult, varName, params%routine, params%isFg)
    END IF
  END SUBROUTINE fetch2d

  SUBROUTINE fetch3d(params, varName, jg, field)
    TYPE(t_fetchParams), INTENT(INOUT) :: params
    CHARACTER(LEN = *), INTENT(IN) :: varName
    INTEGER, VALUE :: jg
    REAL(wp), INTENT(INOUT) :: field(:,:,:)

    LOGICAL :: fetchResult

    IF(params%inputInstructions(jg)%ptr%wantVar(varName, params%isFg)) THEN
        fetchResult = params%requestList%fetch3d(varName, trivial_tileId, jg, field)
        CALL params%inputInstructions(jg)%ptr%handleError(fetchResult, varName, params%routine, params%isFg)
    END IF
  END SUBROUTINE fetch3d

  FUNCTION fetchSurfaceOptional(params, varName, jg, field) RESULT(RESULT)
    TYPE(t_fetchParams), INTENT(INOUT) :: params
    CHARACTER(LEN = *), INTENT(IN) :: varName
    INTEGER, VALUE :: jg
    REAL(wp), INTENT(INOUT) :: field(:,:)
    LOGICAL :: RESULT

    IF(params%inputInstructions(jg)%ptr%wantVar(varName, params%isFg)) THEN
        RESULT = params%requestList%fetchSurface(varName, trivial_tileId, jg, field)
        CALL params%inputInstructions(jg)%ptr%optionalReadResult(RESULT, varName, params%routine, params%isFg)
    END IF
  END FUNCTION fetchSurfaceOptional

  SUBROUTINE fetchSurface(params, varName, jg, field)
    TYPE(t_fetchParams), INTENT(INOUT) :: params
    CHARACTER(LEN = *), INTENT(IN) :: varName
    INTEGER, VALUE :: jg
    REAL(wp), INTENT(INOUT) :: field(:,:)

    LOGICAL :: fetchResult

    IF(params%inputInstructions(jg)%ptr%wantVar(varName, params%isFg)) THEN
        fetchResult = params%requestList%fetchSurface(varName, trivial_tileId, jg, field)
        CALL params%inputInstructions(jg)%ptr%handleError(fetchResult, varName, params%routine, params%isFg)
    END IF
  END SUBROUTINE fetchSurface

  ! Wrapper for requestList%fetchTiledSurface() that falls back to reading copies of untiled input IF ltile_coldstart IS set.
  SUBROUTINE fetchTiledSurface(params, varName, jg, tileCount, field)
    TYPE(t_fetchParams), INTENT(INOUT) :: params
    CHARACTER(LEN = *), INTENT(IN) :: varName
    INTEGER, VALUE :: jg, tileCount
    REAL(wp), INTENT(INOUT) :: field(:,:,:)

    INTEGER :: jt
    LOGICAL :: fetchResult

    IF(params%inputInstructions(jg)%ptr%wantVar(varName, params%isFg)) THEN
        IF(ltile_coldstart) THEN
            !Fake tiled input by copying the input field to all tiles.
            fetchResult = .TRUE.
            DO jt = 1, tileCount
                fetchResult = fetchResult.AND.params%requestList%fetchSurface(varName, trivial_tileId, jg, field(:,:,jt))
            END DO
        ELSE
            !True tiled input.
            fetchResult = params%requestList%fetchTiledSurface(varName, jg, field)
        END IF
        CALL params%inputInstructions(jg)%ptr%handleError(fetchResult, varName, params%routine, params%isFg)
    END IF
  END SUBROUTINE fetchTiledSurface

  ! Wrapper for requestList%fetchTiled3d() that falls back to reading copies of untiled input IF ltile_coldstart IS set.
  SUBROUTINE fetchTiled3d(params, varName, jg, tileCount, field)
    TYPE(t_fetchParams), INTENT(INOUT) :: params
    CHARACTER(LEN = *), INTENT(IN) :: varName
    INTEGER, VALUE :: jg, tileCount
    REAL(wp), INTENT(INOUT) :: field(:,:,:,:)

    INTEGER :: jt
    LOGICAL :: fetchResult

    IF(params%inputInstructions(jg)%ptr%wantVar(varName, params%isFg)) THEN
        IF(ltile_coldstart) THEN
            !Fake tiled input by copying the input field to all tiles.
            fetchResult = .TRUE.
            DO jt = 1, tileCount
                fetchResult = fetchResult.AND.params%requestList%fetch3d(varName, trivial_tileId, jg, field(:,:,:,jt))
            END DO
        ELSE
            !True tiled input.
            fetchResult = params%requestList%fetchTiled3d(varName, jg, field)
        END IF
        CALL params%inputInstructions(jg)%ptr%handleError(fetchResult, varName, params%routine, params%isFg)
    END IF
  END SUBROUTINE fetchTiled3d

  ! Same as above for variables that are optionally read from input file.
  SUBROUTINE fetchTiled3dOptional(params, varName, jg, tileCount, field, avail)
    TYPE(t_fetchParams), INTENT(INOUT) :: params
    CHARACTER(LEN = *), INTENT(IN) :: varName
    INTEGER, VALUE :: jg, tileCount
    REAL(wp), INTENT(INOUT) :: field(:,:,:,:)
    LOGICAL, INTENT(OUT) :: avail

    INTEGER :: jt

    IF(params%inputInstructions(jg)%ptr%wantVar(varName, params%isFg)) THEN
        IF(ltile_coldstart) THEN
            !Fake tiled input by copying the input field to all tiles.
            avail = .TRUE.
            DO jt = 1, tileCount
                avail = avail.AND.params%requestList%fetch3d(varName, trivial_tileId, jg, field(:,:,:,jt))
            END DO
        ELSE
            !True tiled input.
            avail = params%requestList%fetchTiled3d(varName, jg, field)
        END IF
        CALL params%inputInstructions(jg)%ptr%optionalReadResult(avail, varName, params%routine, params%isFg)
    END IF
  END SUBROUTINE fetchTiled3dOptional

  SUBROUTINE fetchRequired3d(params, varName, jg, field)
    TYPE(t_fetchParams), INTENT(INOUT) :: params
    CHARACTER(LEN = *), INTENT(IN) :: varName
    INTEGER, VALUE :: jg
    REAL(wp), INTENT(INOUT) :: field(:,:,:)

    CALL params%requestList%fetchRequired3d(varName, trivial_tileId, jg, field)
    CALL params%inputInstructions(jg)%ptr%handleError(.TRUE., varName, params%routine, params%isFg)
  END SUBROUTINE fetchRequired3d

  ! Wrapper for requestList%fetchRequiredTiledSurface() that falls back to reading copies of untiled input IF ltile_coldstart IS set.
  SUBROUTINE fetchRequiredTiledSurface(params, varName, jg, tileCount, field)
    TYPE(t_fetchParams), INTENT(INOUT) :: params
    CHARACTER(LEN = *), INTENT(IN) :: varName
    INTEGER, VALUE :: jg, tileCount
    REAL(wp), INTENT(INOUT) :: field(:,:,:)

    INTEGER :: jt

    IF(ltile_coldstart) THEN
        !Fake tiled input by copying the input field to all tiles.
        DO jt = 1, tileCount
            CALL params%requestList%fetchRequiredSurface(varName, trivial_tileId, jg, field(:,:,jt))
        END DO
    ELSE
        !True tiled input.
        CALL params%requestList%fetchRequiredTiledSurface(varName, jg, field)
    END IF
    CALL params%inputInstructions(jg)%ptr%handleError(.TRUE., varName, params%routine, params%isFg)
  END SUBROUTINE fetchRequiredTiledSurface


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
    TYPE(t_fetchParams) :: params

    params%inputInstructions = inputInstructions
    params%requestList => requestList
    params%routine = routine
    params%isFg = .TRUE.

    DO jg = 1, n_dom
        IF(p_patch(jg)%ldom_active) THEN
            ! save some paperwork
            prognosticFields => p_nh_state(jg)%prog(nnow(jg))

            ! request the first guess fields (atmosphere only)
            CALL fetchRequired3d(params, 'theta_v', jg, prognosticFields%theta_v)
            CALL fetchRequired3d(params, 'rho', jg, prognosticFields%rho)
            CALL fetchRequired3d(params, 'w', jg, prognosticFields%w)
            CALL fetchRequired3d(params, 'tke', jg, prognosticFields%tke)

            ! Only needed for FG-only runs; usually read from ANA
            my_ptr3d => prognosticFields%tracer(:,:,:,iqv)
            CALL fetch3d(params, 'qv', jg, my_ptr3d)

            my_ptr3d => prognosticFields%tracer(:,:,:,iqc)
            CALL fetch3d(params, 'qc', jg, my_ptr3d)

            my_ptr3d => prognosticFields%tracer(:,:,:,iqi)
            CALL fetch3d(params, 'qi', jg, my_ptr3d)

            IF ( iqr /= 0 ) THEN
              my_ptr3d => prognosticFields%tracer(:,:,:,iqr)
              CALL fetch3d(params, 'qr', jg, my_ptr3d)
            END IF

            IF ( iqs /= 0 ) THEN
              my_ptr3d => prognosticFields%tracer(:,:,:,iqs)
              CALL fetch3d(params, 'qs', jg, my_ptr3d)
            END IF

            IF (lvert_remap_fg) THEN
                CALL allocate_extana_atm(jg, p_patch(jg)%nblks_c, p_patch(jg)%nblks_e, initicon)
                CALL fetchRequired3d(params, 'z_ifc', jg, initicon(jg)%atm_in%z3d_ifc)
            END IF

            CALL fetchRequired3d(params, 'vn', jg, prognosticFields%vn)
        END IF
    END DO

  END SUBROUTINE fetch_dwdfg_atm

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
  SUBROUTINE fetch_dwdfg_atm_ii(requestList, p_patch, initicon, inputInstructions)
    CLASS(t_InputRequestList), POINTER, INTENT(INOUT) :: requestList
    TYPE(t_patch), INTENT(IN) :: p_patch(:)
    TYPE(t_initicon_state), INTENT(INOUT), TARGET :: initicon(:)
    TYPE(t_readInstructionListPtr), INTENT(INOUT) :: inputInstructions(:)

    CHARACTER(len=*), PARAMETER :: routine = modname//':fetch_dwdfg_atm_ii'
    INTEGER :: jg, jb, jk, jc, i_endidx
    REAL(wp) :: tempv, exner
    REAL(dp), POINTER :: levelValues(:)
    TYPE(t_fetchParams) :: params

    params%inputInstructions = inputInstructions
    params%requestList => requestList
    params%routine = routine
    params%isFg = .TRUE.

    DO jg = 1, n_dom
        IF(p_patch(jg)%ldom_active) THEN  ! Skip reading the atmospheric input data if a model domain is not active at initial time

            ! determine number of HALF LEVELS of generalized Z-AXIS
            levelValues => requestList%getLevels('z_ifc', jg)
            IF(.NOT.ASSOCIATED(levelValues)) CALL finish(routine, "no DATA found for domain "//TRIM(int2string(jg))//" of &
                                                                  &required variable 'z_ifc'")
            nlevatm_in(jg) = SIZE(levelValues, 1) - 1

            CALL allocate_extana_atm(jg, p_patch(jg)%nblks_c, p_patch(jg)%nblks_e, initicon)

            ! start reading first guess (atmosphere only)
            CALL fetchRequired3d(params, 'z_ifc', jg, initicon(jg)%atm_in%z3d_ifc)
            CALL fetchRequired3d(params, 'theta_v', jg, initicon(jg)%atm_in%theta_v)
            CALL fetchRequired3d(params, 'rho', jg, initicon(jg)%atm_in%rho)
            CALL fetchRequired3d(params, 'w', jg, initicon(jg)%atm_in%w_ifc)
            CALL fetchRequired3d(params, 'tke', jg, initicon(jg)%atm_in%tke_ifc)

            CALL fetchRequired3d(params, 'qv', jg, initicon(jg)%atm_in%qv)
            CALL fetchRequired3d(params, 'qc', jg, initicon(jg)%atm_in%qc)
            CALL fetchRequired3d(params, 'qi', jg, initicon(jg)%atm_in%qi)
            IF ( iqr /= 0 ) THEN
            CALL fetchRequired3d(params, 'qr', jg, initicon(jg)%atm_in%qr)
            ELSE
            initicon(jg)%atm_in%qr(:,:,:) = 0._wp
            END IF

            IF ( iqs /= 0 ) THEN
            CALL fetchRequired3d(params, 'qs', jg, initicon(jg)%atm_in%qs)
            ELSE
            initicon(jg)%atm_in%qs(:,:,:) = 0._wp
            END IF

            CALL fetchRequired3d(params, 'vn', jg, initicon(jg)%atm_in%vn)

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
  !! Fetch DA-analysis DATA from the request list (atmosphere only)
  !!
  !! Depending on the initialization mode, either full fields or increments
  !! are read (atmosphere only). The following full fields are read, if available:
  !!     u, v, t, p, qv, qi, qc, qr, qs
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
    TYPE(t_fetchParams) :: params
    LOGICAL :: lHaveFg

    params%inputInstructions = inputInstructions
    params%requestList => requestList
    params%routine = routine
    params%isFg = .FALSE.

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
            CALL fetch3d(params, 'temp', jg, my_ptr%temp)
            CALL fetch3d(params, 'pres', jg, my_ptr%pres)
            CALL fetch3d(params, 'u', jg, my_ptr%u)
            CALL fetch3d(params, 'v', jg, my_ptr%v)

            IF ( ANY((/MODE_IAU,MODE_IAU_OLD/) == init_mode) ) THEN
                lHaveFg = inputInstructions(jg)%ptr%sourceOfVar('qv') == kInputSourceFg
                CALL fetch3d(params, 'qv', jg, my_ptr%qv)
                ! check whether we are using DATA from both FG AND ANA input, so that it's correctly listed IN the input source table
                IF(lHaveFg.AND.inputInstructions(jg)%ptr%sourceOfVar('qv') == kInputSourceAna) THEN
                    CALL inputInstructions(jg)%ptr%setSource('qv', kInputSourceBoth)
                END IF
            ELSE
                my_ptr3d => p_nh_state(jg)%prog(nnow(jg))%tracer(:,:,:,iqv)
                CALL fetch3d(params, 'qv', jg, my_ptr3d)
            ENDIF

            ! For the time being, these are identical to qc, qi, qr, AND qs from FG => usually read from FG
            my_ptr3d => p_nh_state(jg)%prog(nnow(jg))%tracer(:,:,:,iqc)
            CALL fetch3d(params, 'qc', jg, my_ptr3d)
            my_ptr3d => p_nh_state(jg)%prog(nnow(jg))%tracer(:,:,:,iqi)
            CALL fetch3d(params, 'qi', jg, my_ptr3d)
            IF ( iqr /= 0 ) THEN
                my_ptr3d => p_nh_state(jg)%prog(nnow(jg))%tracer(:,:,:,iqr)
                CALL fetch3d(params, 'qr', jg, my_ptr3d)
            END IF
            IF ( iqs /= 0 ) THEN
                my_ptr3d => p_nh_state(jg)%prog(nnow(jg))%tracer(:,:,:,iqs)
                CALL fetch3d(params, 'qs', jg, my_ptr3d)
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
  !! Fetch DWD first guess DATA from request list (land/surface only)
  SUBROUTINE fetch_dwdfg_sfc(requestList, p_patch, prm_diag, p_lnd_state, inputInstructions)
    CLASS(t_InputRequestList), POINTER, INTENT(INOUT) :: requestList
    TYPE(t_patch),             INTENT(IN)    :: p_patch(:)
    TYPE(t_nwp_phy_diag),      INTENT(INOUT) :: prm_diag(:)
    TYPE(t_lnd_state), TARGET, INTENT(INOUT) :: p_lnd_state(:)
    TYPE(t_readInstructionListPtr), INTENT(INOUT) :: inputInstructions(:)

    INTEGER :: jg, error
    LOGICAL :: avail
    REAL(wp), POINTER :: my_ptr2d(:,:)
    TYPE(t_lnd_prog), POINTER :: lnd_prog
    TYPE(t_lnd_diag), POINTER :: lnd_diag
    TYPE(t_wtr_prog), POINTER :: wtr_prog
    TYPE(t_fetchParams) :: params

    CHARACTER(len=*), PARAMETER :: routine = modname//':fetch_dwdfg_sfc'

    ! Workaround for the intel compiler botching implicit allocation.
    ALLOCATE(params%inputInstructions(SIZE(inputInstructions, 1)), STAT = error)
    IF(error /= SUCCESS) CALL finish(routine, "memory allocation failed")

    params%inputInstructions = inputInstructions
    params%requestList => requestList
    params%routine = routine
    params%isFg = .TRUE.

    DO jg = 1, n_dom
        IF(p_patch(jg)%ldom_active) THEN ! Skip reading the atmospheric input data if a model domain is not active at initial time
            lnd_prog => p_lnd_state(jg)%prog_lnd(nnow_rcf(jg))
            lnd_diag => p_lnd_state(jg)%diag_lnd
            wtr_prog => p_lnd_state(jg)%prog_wtr(nnow_rcf(jg))

            ! COSMO-DE does not provide sea ice field. In that case set fr_seaice to 0
            IF (init_mode /= MODE_COSMODE) THEN
                CALL fetchSurface(params, 'fr_seaice', jg, lnd_diag%fr_seaice)
            ELSE
!$OMP PARALLEL
                CALL init(lnd_diag%fr_seaice(:,:))
!$OMP END PARALLEL
            ENDIF ! init_mode /= MODE_COSMODE

            ! sea-ice related fields
            CALL fetchSurface(params, 't_ice', jg, wtr_prog%t_ice)
            CALL fetchSurface(params, 'h_ice', jg, wtr_prog%h_ice)

            !These two fields are required for the processing step below, AND they are NOT initialized before this SUBROUTINE IS called, so they are fetched as required.
            !This diverges from the code that I found which READ them conditionally.
            CALL fetchRequiredTiledSurface(params, 't_g', jg, ntiles_total + ntiles_water, lnd_prog%t_g_t)
            CALL fetchRequiredTiledSurface(params, 'qv_s', jg, ntiles_total + ntiles_water, lnd_diag%qv_s_t)

            CALL fetchTiledSurface(params, 'freshsnow', jg, ntiles_total, lnd_diag%freshsnow_t)
            IF (lsnowtile .AND. .NOT. ltile_coldstart) THEN
                CALL fetchTiledSurface(params, 'snowfrac_lc', jg, ntiles_total, lnd_diag%snowfrac_lc_t)
            ENDIF
            CALL fetchTiledSurface(params, 'w_snow', jg, ntiles_total, lnd_prog%w_snow_t)
            CALL fetchTiledSurface(params, 'w_i', jg, ntiles_total, lnd_prog%w_i_t)
            CALL fetchTiledSurface(params, 'h_snow', jg, ntiles_total, lnd_diag%h_snow_t)
            CALL fetchTiledSurface(params, 't_snow', jg, ntiles_total, lnd_prog%t_snow_t)
            CALL fetchTiledSurface(params, 'rho_snow', jg, ntiles_total, lnd_prog%rho_snow_t)

            IF (lmulti_snow) THEN
                ! multi layer snow fields
                CALL fetchTiled3d(params, 't_snow_mult', jg, ntiles_total, lnd_prog%t_snow_mult_t)
                CALL fetchTiled3d(params, 'rho_snow_mult', jg, ntiles_total, lnd_prog%rho_snow_mult_t)
                CALL fetchTiled3d(params, 'wtot_snow', jg, ntiles_total, lnd_prog%wtot_snow_t)
                CALL fetchTiled3d(params, 'wliq_snow', jg, ntiles_total, lnd_prog%wliq_snow_t)
                CALL fetchTiled3d(params, 'dzh_snow', jg, ntiles_total, lnd_prog%dzh_snow_t)
            ELSE IF (l2lay_rho_snow) THEN
                CALL fetchTiled3dOptional(params, 'rho_snow_mult', jg, ntiles_total, lnd_prog%rho_snow_mult_t, avail)
                IF(.NOT. avail) THEN
                    ! initialize top-layer snow density with average density if no input field is available
                    lnd_prog%rho_snow_mult_t(:,1,:,:) = lnd_prog%rho_snow_t(:,:,:)
                ENDIF
            END IF ! lmulti_snow


            ! multi layer fields
            !
            ! Note that either w_so OR smi is written to w_so_t. Which one is required depends
            ! on the initialization mode. Checking grp_vars_fg takes care of this. In case
            ! that smi is read, it is lateron converted to w_so (see smi_to_wsoil)
            SELECT CASE(init_mode)
                CASE(MODE_COMBINED, MODE_COSMODE)
                    CALL fetchTiled3d(params, 'smi', jg, ntiles_total, lnd_prog%w_so_t)
                CASE DEFAULT
                    CALL fetchTiled3d(params, 'w_so', jg, ntiles_total, lnd_prog%w_so_t)
                    CALL fetchTiled3d(params, 'w_so_ice', jg, ntiles_total, lnd_prog%w_so_ice_t) ! w_so_ice is re-diagnosed in terra_multlay_init
            END SELECT

            CALL fetchTiled3d(params, 't_so', jg, ntiles_total, lnd_prog%t_so_t)

            ! Skipped in MODE_COMBINED and in MODE_COSMODE (i.e. when starting from GME soil)
            ! Instead z0 is re-initialized (see mo_nwp_phy_init)
            CALL fetchSurface(params, 'gz0', jg, prm_diag(jg)%gz0)

            ! first guess for fresh water lake fields
            IF(ASSOCIATED(wtr_prog%t_mnw_lk)) CALL fetchSurface(params, 't_mnw_lk', jg, wtr_prog%t_mnw_lk)
            IF(ASSOCIATED(wtr_prog%t_wml_lk)) CALL fetchSurface(params, 't_wml_lk', jg, wtr_prog%t_wml_lk)
            IF(ASSOCIATED(wtr_prog%h_ml_lk)) CALL fetchSurface(params, 'h_ml_lk', jg, wtr_prog%h_ml_lk)
            IF(ASSOCIATED(wtr_prog%t_bot_lk)) CALL fetchSurface(params, 't_bot_lk', jg, wtr_prog%t_bot_lk)
            IF(ASSOCIATED(wtr_prog%c_t_lk)) CALL fetchSurface(params, 'c_t_lk', jg, wtr_prog%c_t_lk)
            IF(ASSOCIATED(wtr_prog%t_b1_lk)) CALL fetchSurface(params, 't_b1_lk', jg, wtr_prog%t_b1_lk)
            IF(ASSOCIATED(wtr_prog%h_b1_lk)) CALL fetchSurface(params, 'h_b1_lk', jg, wtr_prog%h_b1_lk)

            IF(iprog_aero == 1) THEN
                aerosol_fg_present(jg) = .TRUE.

                my_ptr2d => prm_diag(jg)%aerosol(:,iss,:)
                IF(.NOT.fetchSurfaceOptional(params, 'aer_ss', jg, my_ptr2d)) aerosol_fg_present(jg) = .FALSE.
                my_ptr2d => prm_diag(jg)%aerosol(:,iorg,:)
                IF(.NOT.fetchSurfaceOptional(params, 'aer_or', jg, my_ptr2d)) aerosol_fg_present(jg) = .FALSE.
                my_ptr2d => prm_diag(jg)%aerosol(:,ibc,:)
                IF(.NOT.fetchSurfaceOptional(params, 'aer_bc', jg, my_ptr2d)) aerosol_fg_present(jg) = .FALSE.
                my_ptr2d => prm_diag(jg)%aerosol(:,iso4,:)
                IF(.NOT.fetchSurfaceOptional(params, 'aer_su', jg, my_ptr2d)) aerosol_fg_present(jg) = .FALSE.
                my_ptr2d => prm_diag(jg)%aerosol(:,idu,:)
                IF(.NOT.fetchSurfaceOptional(params, 'aer_du', jg, my_ptr2d)) aerosol_fg_present(jg) = .FALSE.

                ! Reading of these variables IS purely OPTIONAL, but IF reading of one fails, we USE NONE. Inform the InputInstructions about this.
                IF(.NOT.aerosol_fg_present(jg)) THEN
                    CALL inputInstructions(jg)%ptr%setSource('aer_ss', kInputSourceNone)
                    CALL inputInstructions(jg)%ptr%setSource('aer_or', kInputSourceNone)
                    CALL inputInstructions(jg)%ptr%setSource('aer_bc', kInputSourceNone)
                    CALL inputInstructions(jg)%ptr%setSource('aer_su', kInputSourceNone)
                    CALL inputInstructions(jg)%ptr%setSource('aer_du', kInputSourceNone)
                END IF
            ELSE
                aerosol_fg_present(jg) = .FALSE.
            END IF
        END IF
    END DO

  END SUBROUTINE fetch_dwdfg_sfc


  !! processing of DWD first-guess DATA
  SUBROUTINE process_input_dwdfg_sfc(p_patch, p_lnd_state, ext_data)
    TYPE(t_patch),             INTENT(IN)    :: p_patch(:)
    TYPE(t_lnd_state), TARGET, INTENT(INOUT) :: p_lnd_state(:)
    TYPE(t_external_data)    , INTENT(INOUT) :: ext_data(:)

    INTEGER :: jg, jt, jb, jc, ic, i_endidx
    TYPE(t_lnd_prog), POINTER :: lnd_prog
    TYPE(t_wtr_prog), POINTER :: wtr_prog
    TYPE(t_lnd_diag), POINTER :: lnd_diag
    CHARACTER(LEN = *), PARAMETER :: routine = modname//":process_input_dwdfg_sfc"

    DO jg = 1, n_dom
        IF(p_patch(jg)%ldom_active) THEN
            lnd_prog => p_lnd_state(jg)%prog_lnd(nnow_rcf(jg))
            wtr_prog => p_lnd_state(jg)%prog_wtr(nnow_rcf(jg))
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


            ! Initialize t_s_t, which is not read in
            !
!$OMP PARALLEL DO PRIVATE(jb,jt,ic,jc)
            DO jb = 1, p_patch(jg)%nblks_c

              ! set t_s_t to t_so_t(1) for land tiles
              DO jt = 1, ntiles_total
                DO ic = 1, ext_data(jg)%atm%lp_count_t(jb,jt)
                  jc = ext_data(jg)%atm%idx_lst_lp_t(ic,jb,jt)
                  lnd_prog%t_s_t(jc,jb,jt) = lnd_prog%t_so_t(jc,1,jb,jt)
                ENDDO
              ENDDO  ! ntiles

              IF (llake) THEN
                ! take lake surface temperature from t_wml_lk
                DO ic = 1, ext_data(jg)%atm%fp_count(jb)
                  jc = ext_data(jg)%atm%idx_lst_fp(ic,jb)
                  lnd_prog%t_s_t(jc,jb,isub_lake) = wtr_prog%t_wml_lk(jc,jb)
                ENDDO
              ELSE
                ! without lake model, take t_g as an estimate of the lake temperature
                DO ic = 1, ext_data(jg)%atm%fp_count(jb)
                  jc = ext_data(jg)%atm%idx_lst_fp(ic,jb)
                  lnd_prog%t_s_t(jc,jb,isub_lake) = MAX(tmelt, lnd_prog%t_g_t(jc,jb,isub_lake))
                ENDDO
              ENDIF

              ! NOTE: Initialization of sea-water and sea-ice tiles 
              ! is done later in mo_nwp_sfc_utils:nwp_surface_init for 
              ! two reasons:
              ! I)  index lists for sea-ice and sea-water are 
              !     not yet available at this point.
              ! II) If an SST-analysis is read in, t_s_t(isub_water) 
              !     must be updated another time, anyway.

            ENDDO  ! jb
!$OMP END PARALLEL DO
        END IF

    END DO
  END SUBROUTINE process_input_dwdfg_sfc


  !>
  !! Fetch DWD analysis DATA from the request list (land/surface only)
  !!
  !! Analysis is read for:
  !! fr_seaice, t_ice, h_ice, t_so, h_snow, w_snow, w_i, t_snow, rho_snow, freshsnow, w_so
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

    INTEGER :: jg
    INTEGER, PARAMETER :: jt = 1
    REAL(wp), POINTER :: my_ptr2d(:,:)
    REAL(wp), POINTER :: my_ptr3d(:,:,:)
    TYPE(t_lnd_prog), POINTER :: lnd_prog
    TYPE(t_lnd_diag), POINTER :: lnd_diag
    TYPE(t_wtr_prog), POINTER :: wtr_prog
    TYPE(t_fetchParams) :: params
    LOGICAL :: lHaveFg

    CHARACTER(LEN = *), PARAMETER :: routine = modname//':fetch_dwdana_sfc'

    params%inputInstructions = inputInstructions
    params%requestList => requestList
    params%routine = routine
    params%isFg = .FALSE.

    DO jg = 1, n_dom
        IF(p_patch(jg)%ldom_active) THEN
            IF (ANY((/MODE_IAU, MODE_IAU_OLD /) == init_mode) .AND. lp2cintp_sfcana(jg)) CYCLE  ! Skip this domain if it will be interpolated from its parent domain in process_input_dwdana_sfc()

            lnd_prog => p_lnd_state(jg)%prog_lnd(nnow_rcf(jg))
            lnd_diag => p_lnd_state(jg)%diag_lnd
            wtr_prog => p_lnd_state(jg)%prog_wtr(nnow_rcf(jg))

            ! sea-ice fraction
            CALL fetchSurface(params, 'fr_seaice', jg, lnd_diag%fr_seaice)
            ! sea-ice temperature
            CALL fetchSurface(params, 't_ice', jg, wtr_prog%t_ice)
            ! sea-ice height
            CALL fetchSurface(params, 'h_ice', jg, wtr_prog%h_ice)

            ! sea-surface temperature: fetch T_SEA or, alternatively, T_SO(0)
            my_ptr2d => initicon(jg)%sfc%sst(:,:)
            CALL fetch2d(params, 't_seasfc', 0.0_wp, jg, my_ptr2d)
            !
            IF (inputInstructions(jg)%ptr%fetchStatus('t_seasfc', lIsFg=.FALSE.) == kStateFailedFetch) THEN

              ! T_SO(0). Note that the file may contain a 3D field, of which we ONLY fetch the level at 0.0.
              lHaveFg = inputInstructions(jg)%ptr%sourceOfVar('t_so') == kInputSourceFg
              CALL fetch2d(params, 't_so', 0.0_wp, jg, my_ptr2d)
              ! check whether we are using DATA from both FG AND ANA input, so that it's correctly listed IN the input source table
              IF(lHaveFg.AND.inputInstructions(jg)%ptr%sourceOfVar('t_so') == kInputSourceAna) THEN
                  CALL inputInstructions(jg)%ptr%setSource('t_so', kInputSourceBoth)
              END IF
            ENDIF


            ! h_snow
            IF ( init_mode == MODE_IAU ) THEN
                lHaveFg = inputInstructions(jg)%ptr%sourceOfVar('h_snow') == kInputSourceFg
                my_ptr2d => initicon(jg)%sfc_inc%h_snow(:,:)
                CALL fetchSurface(params, 'h_snow', jg, my_ptr2d)
                ! check whether we are using DATA from both FG AND ANA input, so that it's correctly listed IN the input source table
                IF(lHaveFg.AND.inputInstructions(jg)%ptr%sourceOfVar('h_snow') == kInputSourceAna) THEN
                    CALL inputInstructions(jg)%ptr%setSource('h_snow', kInputSourceBoth)
                END IF
            ELSE
                my_ptr2d => lnd_diag%h_snow_t(:,:,jt)
                CALL fetchSurface(params, 'h_snow', jg, my_ptr2d)
            ENDIF

            ! w_snow
            my_ptr2d => lnd_prog%w_snow_t(:,:,jt)
            CALL fetchSurface(params, 'w_snow', jg, my_ptr2d)

            ! w_i
            my_ptr2d => lnd_prog%w_i_t(:,:,jt)
            CALL fetchSurface(params, 'w_i', jg, my_ptr2d)

            ! t_snow
            my_ptr2d => lnd_prog%t_snow_t(:,:,jt)
            CALL fetchSurface(params, 't_snow', jg, my_ptr2d)

            ! rho_snow
            my_ptr2d => lnd_prog%rho_snow_t(:,:,jt)
            CALL fetchSurface(params, 'rho_snow', jg, my_ptr2d)

            ! freshsnow
            IF ( init_mode == MODE_IAU ) THEN
                lHaveFg = inputInstructions(jg)%ptr%sourceOfVar('freshsnow') == kInputSourceFg
                my_ptr2d => initicon(jg)%sfc_inc%freshsnow(:,:)
                CALL fetchSurface(params, 'freshsnow', jg, my_ptr2d)
                ! check whether we are using DATA from both FG AND ANA input, so that it's correctly listed IN the input source table
                IF(lHaveFg.AND.inputInstructions(jg)%ptr%sourceOfVar('freshsnow') == kInputSourceAna) THEN
                    CALL inputInstructions(jg)%ptr%setSource('freshsnow', kInputSourceBoth)
                END IF
            ELSE
                my_ptr2d => lnd_diag%freshsnow_t(:,:,jt)
                CALL fetchSurface(params, 'freshsnow', jg, my_ptr2d)
            ENDIF

            ! w_so
            IF ( (init_mode == MODE_IAU) .OR. (init_mode == MODE_IAU_OLD) ) THEN
                lHaveFg = inputInstructions(jg)%ptr%sourceOfVar('w_so') == kInputSourceFg
                my_ptr3d => initicon(jg)%sfc_inc%w_so(:,:,:)
                CALL fetch3d(params, 'w_so', jg, my_ptr3d)
                ! check whether we are using DATA from both FG AND ANA input, so that it's correctly listed IN the input source table
                IF(lHaveFg.AND.inputInstructions(jg)%ptr%sourceOfVar('w_so') == kInputSourceAna) THEN
                    CALL inputInstructions(jg)%ptr%setSource('w_so', kInputSourceBoth)
                END IF
            ELSE
                my_ptr3d => lnd_prog%w_so_t(:,:,:,jt)
                CALL fetch3d(params, 'w_so', jg, my_ptr3d)
            END IF
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

