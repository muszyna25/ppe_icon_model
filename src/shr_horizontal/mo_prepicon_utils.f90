!>
!! This module contains the I/O routines for prep_icon and the driver
!! routines for topography preprocessing
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
MODULE mo_prepicon_utils

  USE mo_kind
  USE mo_io_units,            ONLY: filename_max
  USE mo_parallel_config,     ONLY: nproma, p_test_run
  USE mo_run_config,          ONLY: num_lev, num_levp1, msg_level, &
    &                                dtime, nvclev
  USE mo_extpar_config,       ONLY: n_iter_smooth_topo
  USE mo_dynamics_config,     ONLY: iequations
  USE mo_nonhydrostatic_config,ONLY: ivctype
  USE mo_prepicon_nml,        ONLY: i_oper_mode, nlev_in, l_zp_out, l_w_in,&
  &                                 nlevsoil_in, l_sfc_in
  USE mo_model_domain,        ONLY: t_patch
  USE mo_impl_constants,      ONLY: SUCCESS, max_char_length, max_dom
  USE mo_exception,           ONLY: message, finish, message_text
  USE mo_grid_config,         ONLY: n_dom, nroot, start_lev, global_cell_type
  USE mo_interpolation,       ONLY: t_int_state, cells2verts_scalar
  USE mo_grf_interpolation,   ONLY: t_gridref_state
  USE mo_mpi,                 ONLY: p_pe, p_io, p_bcast, p_comm_work_test, p_comm_work
  USE mo_ext_data,            ONLY: smooth_topography, read_netcdf_data
  USE mo_atmo_control,        ONLY: p_patch
  USE mo_io_vlist,            ONLY: GATHER_C, GATHER_E, GATHER_V, num_output_vars, outvar_desc, &
                                    gather_array1, gather_array2
  USE mo_datetime,            ONLY: t_datetime, print_datetime
  USE mo_io_config,              ONLY: lkeep_in_sync
  USE mo_nh_init_utils,       ONLY: nflat, nflatlev, compute_smooth_topo, init_vert_coord, &
                                    topography_blending, topography_feedback,              &
                                    interp_uv_2_vn, init_w, convert_thdvars, virtual_temp
  USE mo_model_domain_import, ONLY: lfeedback
  USE mo_ifs_coord,           ONLY: alloc_vct, init_vct, vct, vct_a, vct_b
  USE mo_lnd_nwp_config,      ONLY: nlev_soil

  IMPLICIT NONE

  INCLUDE 'cdi.inc'
  INCLUDE 'netcdf.inc'

  PRIVATE

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

  TYPE :: t_pi_atm_in ! atmospheric input variables; surface geopotential is regarded as
                      ! atmospheric variable here because the atmospheric fields cannot be processed without it

    REAL(wp), ALLOCATABLE, DIMENSION(:,:)   :: psfc, phi_sfc
    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: temp, pres, z3d, u, v, omega, &
      &                                         w, qv, qc, qi, qr, qs

  END TYPE t_pi_atm_in

  TYPE :: t_pi_sfc_in

    REAL(wp), ALLOCATABLE, DIMENSION (:,:) :: tsnow, tskin, snowweq, snowdens, &
                                              skinres, ls_mask, seaice
    REAL(wp), ALLOCATABLE, DIMENSION (:,:,:) :: tsoil, soilwater

  END TYPE t_pi_sfc_in

  TYPE :: t_pi_atm

    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: vn, u, v, w, temp, theta_v, exner, rho, &
                                               pres, qv, qc, qi, qr, qs

  END TYPE t_pi_atm

  TYPE :: t_pi_diag

    REAL(wp), ALLOCATABLE, DIMENSION(:)     :: levels
    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: u, v, temp, pres, qv, z3d

  END TYPE t_pi_diag

  TYPE :: t_pi_sfc

    REAL(wp), ALLOCATABLE, DIMENSION (:,:) :: tsnow, tskin, snowweq, snowdens, &
                                              skinres, ls_mask, seaice
    REAL(wp), ALLOCATABLE, DIMENSION (:,:,:) :: tsoil, soilwater

  END TYPE t_pi_sfc


  TYPE :: t_prepicon_state

    REAL(wp), ALLOCATABLE, DIMENSION (:,:) :: topography_c, topography_c_smt, &
      topography_v, topography_v_smt

    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: z_ifc, z_mc

    TYPE (t_pi_atm_in) :: atm_in
    TYPE (t_pi_sfc_in) :: sfc_in
    TYPE (t_pi_atm)    :: atm
    TYPE (t_pi_diag)   :: plev
    TYPE (t_pi_diag)   :: zlev
    TYPE (t_pi_sfc)    :: sfc

  END TYPE t_prepicon_state

  TYPE(t_prepicon_state), ALLOCATABLE, TARGET :: prepicon(:) 


  INTEGER, PARAMETER :: max_outvars  = 120 ! max. number of output variables
  INTEGER, PARAMETER :: max_gridlevs = 12 ! max. number of grid levels

  ! I/O stream handler
  INTEGER, DIMENSION(max_gridlevs) ::  &
    &  streamID

  INTEGER, DIMENSION(max_gridlevs) ::  &
    &  gridCellID, gridEdgeID, gridVertexID

  INTEGER, DIMENSION(max_gridlevs) ::  &
    &  vlistID, taxisID, zaxisID_surface, zaxisID_hybrid, &
    &  zaxisID_hybrid_half, zaxisID_pres, zaxisID_hgt


  INTEGER, DIMENSION(max_outvars,max_gridlevs) ::  &
    &  varids
  ! current number of output variables, gets updated by addVar()
  INTEGER, PRIVATE :: num_varids(max_dom)

  INTEGER, SAVE :: iostep = 0
  INTEGER :: klev, nzplev

  PUBLIC :: prepicon, t_prepicon_state, t_pi_atm_in, t_pi_atm, t_pi_diag, t_pi_sfc_in, t_pi_sfc

  PUBLIC :: nzplev

  PUBLIC :: init_prepicon, init_topo_output_files, write_prepicon_output,            &
            setup_prepicon_vlist, compute_coord_fields, close_prepicon_output_files, &
            convert_variables, init_atmo_output_files, deallocate_prepicon

  CONTAINS

  !-------------
  !>
  !! SUBROUTINE init_prepicon
  !! Initialization routine of prep_icon:
  !! Reads in data and processes topography blending and feedback in the
  !! presence of nested domains
  !!
  !! @par Revision History
  !! Initial version by Guenther Zaengl, DWD(2011-07-14)
  !!
  !!
  SUBROUTINE init_prepicon (p_int_state, p_grf_state, prepicon)

    TYPE(t_int_state),     TARGET, INTENT(IN) :: p_int_state(:)
    TYPE(t_gridref_state), TARGET, INTENT(IN) :: p_grf_state(:)

    TYPE(t_prepicon_state), INTENT(INOUT) :: prepicon(:)

    INTEGER :: jg, jlev, nlev, nlevp1, nblks_c, nblks_v, nblks_e, jk
    LOGICAL :: l_exist

    INTEGER :: no_cells, no_verts, no_levels
    INTEGER :: ncid, dimid, varid, mpi_comm

    CHARACTER(len=max_char_length), PARAMETER :: &
      routine = 'mo_prepicon_utils:init_prepicon'
    CHARACTER(LEN=filename_max) :: topo_file(max_dom), ifs2icon_file(max_dom)
    CHARACTER(LEN=10) :: psvar
!-------------------------------------------------------------------------

    IF (l_zp_out) THEN ! set number of z and p levels for diagnostic output
      nzplev = 10
    ELSE
      nzplev = 1
    ENDIF

    ! Allocate memory for prep_icon state
    CALL allocate_prepicon (prepicon)


    IF (i_oper_mode >= 2 .AND. l_zp_out) THEN ! set levels for test output on pressure and height levels

      DO jg = 1, n_dom
        ! set levels - attention: ordering of the levels must be top-down as for the model levels
        DO jk = 1, nzplev
          prepicon(jg)%zlev%levels(nzplev+1-jk) = REAL(jk-1,wp)*1000._wp ! every 1000 m
        ENDDO
        ! standard pressure levels
        prepicon(jg)%plev%levels(10) = 100000._wp
        prepicon(jg)%plev%levels(9)  =  92500._wp
        prepicon(jg)%plev%levels(8)  =  85000._wp
        prepicon(jg)%plev%levels(7)  =  70000._wp
        prepicon(jg)%plev%levels(6)  =  50000._wp
        prepicon(jg)%plev%levels(5)  =  40000._wp
        prepicon(jg)%plev%levels(4)  =  30000._wp
        prepicon(jg)%plev%levels(3)  =  25000._wp
        prepicon(jg)%plev%levels(2)  =  20000._wp
        prepicon(jg)%plev%levels(1)  =  10000._wp
      ENDDO

    ENDIF

    DO jg = 1, n_dom

      jlev = p_patch(jg)%level

      IF(p_pe == p_io) THEN
        !
        ! generate file name
        !
        WRITE(topo_file(jg),'(a,i0,2(a,i2.2),a)') 'extpar_R',nroot,'B',jlev,'_DOM',jg,'.nc'

        INQUIRE (FILE=topo_file(jg), EXIST=l_exist)
        IF (.NOT.l_exist) THEN
          CALL finish(TRIM(routine),'Topography file is not found: '//TRIM(topo_file(jg)))
        ENDIF

        !
        ! open file
        !
        CALL nf(nf_open(TRIM(topo_file(jg)), NF_NOWRITE, ncid))

        !
        ! get number of cells and vertices
        !
        CALL nf(nf_inq_dimid(ncid, 'cell', dimid))
        IF (global_cell_type == 3) THEN ! triangular grid
          CALL nf(nf_inq_dimlen(ncid, dimid, no_cells))
        ELSEIF (global_cell_type == 6) THEN ! hexagonal grid
          CALL nf(nf_inq_dimlen(ncid, dimid, no_verts))
        ENDIF

        CALL nf(nf_inq_dimid(ncid, 'vertex', dimid))
        IF (global_cell_type == 3) THEN ! triangular grid
          CALL nf(nf_inq_dimlen(ncid, dimid, no_verts))
        ELSEIF (global_cell_type == 6) THEN ! hexagonal grid
          CALL nf(nf_inq_dimlen(ncid, dimid, no_cells))
        ENDIF

        !
        ! check the number of cells and verts
        !
        IF(p_patch(jg)%n_patch_cells_g /= no_cells) THEN
          CALL finish(TRIM(ROUTINE),&
          & 'Number of patch cells and cells in topography file do not match.')
        ENDIF
        IF(p_patch(jg)%n_patch_verts_g /= no_verts) THEN
          CALL finish(TRIM(ROUTINE),&
          & 'Number of patch verts and verts in topography file do not match.')
        ENDIF
      ENDIF


      !-------------------------------------------------------
      !
      ! Read topography for triangle centers and vertices
      !
      !-------------------------------------------------------

      ! triangle center

      IF (p_patch(jg)%cell_type == 3) THEN     ! triangular grid
        CALL read_netcdf_data (ncid, 'topography_c', p_patch(jg)%n_patch_cells_g,       &
          &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
          &                     prepicon(jg)%topography_c)
      ELSEIF (p_patch(jg)%cell_type == 6) THEN ! hexagonal grid
        CALL read_netcdf_data (ncid, 'topography_c', p_patch(jg)%n_patch_verts_g,       &
          &                     p_patch(jg)%n_patch_verts, p_patch(jg)%verts%glb_index, &
          &                     prepicon(jg)%topography_v)
      ENDIF

      ! triangle vertex
      IF (p_patch(jg)%cell_type == 3) THEN     ! triangular grid
        CALL read_netcdf_data (ncid, 'topography_v', p_patch(jg)%n_patch_verts_g,       &
          &                     p_patch(jg)%n_patch_verts, p_patch(jg)%verts%glb_index, &
          &                     prepicon(jg)%topography_v)
      ELSEIF (p_patch(jg)%cell_type == 6) THEN ! hexagonal grid
        CALL read_netcdf_data (ncid, 'topography_v', p_patch(jg)%n_patch_cells_g,       &
          &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
          &                     prepicon(jg)%topography_c)
      ENDIF

      ! close file
      !
      IF(p_pe == p_io) CALL nf(nf_close(ncid))

      IF(p_pe == p_io .AND. i_oper_mode >= 2) THEN ! Read in data from IFS2ICON
        !
        ! generate file name
        !
        WRITE(ifs2icon_file(jg),'(a,i0,2(a,i2.2),a)') 'ifs2icon_R',nroot,'B',jlev,'_DOM',jg,'.nc'

        INQUIRE (FILE=ifs2icon_file(jg), EXIST=l_exist)
        IF (.NOT.l_exist) THEN
          CALL finish(TRIM(routine),'IFS2ICON file is not found: '//TRIM(ifs2icon_file(jg)))
        ENDIF

        !
        ! open file
        !
        CALL nf(nf_open(TRIM(ifs2icon_file(jg)), NF_NOWRITE, ncid))

        !
        ! get number of cells
        !
        CALL nf(nf_inq_dimid(ncid, 'ncells', dimid))
        CALL nf(nf_inq_dimlen(ncid, dimid, no_cells))

        !
        ! get number of vertical levels
        !
        CALL nf(nf_inq_dimid(ncid, 'lev', dimid))
        CALL nf(nf_inq_dimlen(ncid, dimid, no_levels))

        !
        ! check the number of cells and vertical levels
        !
        IF(p_patch(jg)%n_patch_cells_g /= no_cells) THEN
          CALL finish(TRIM(ROUTINE),&
          & 'Number of patch cells and cells in IFS2ICON file do not match.')
        ENDIF

        IF(nlev_in /= no_levels) THEN
          CALL finish(TRIM(ROUTINE),&
          & 'nlev_in does not match the number of levels in IFS2ICON file.')
        ENDIF

      ENDIF

      IF (i_oper_mode >= 2) THEN

        CALL read_netcdf_data (ncid, 'T', p_patch(jg)%n_patch_cells_g,                  &
          &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
          &                     nlev_in,prepicon(jg)%atm_in%temp)

        CALL read_netcdf_data (ncid, 'U', p_patch(jg)%n_patch_cells_g,                  &
          &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
          &                     nlev_in,prepicon(jg)%atm_in%u)

        CALL read_netcdf_data (ncid, 'V', p_patch(jg)%n_patch_cells_g,                  &
          &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
          &                     nlev_in,prepicon(jg)%atm_in%v)

        IF (l_w_in) THEN ! note: input vertical velocity is in fact omega (Pa/s)
          CALL read_netcdf_data (ncid, 'W', p_patch(jg)%n_patch_cells_g,                  &
            &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
            &                     nlev_in,prepicon(jg)%atm_in%omega)
        ENDIF

        CALL read_netcdf_data (ncid, 'QV', p_patch(jg)%n_patch_cells_g,                 &
          &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
          &                     nlev_in,prepicon(jg)%atm_in%qv)

        CALL read_netcdf_data (ncid, 'QC', p_patch(jg)%n_patch_cells_g,                 &
          &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
          &                     nlev_in,prepicon(jg)%atm_in%qc)

        CALL read_netcdf_data (ncid, 'QI', p_patch(jg)%n_patch_cells_g,                 &
          &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
          &                     nlev_in,prepicon(jg)%atm_in%qi)

        CALL read_netcdf_data (ncid, 'QR', p_patch(jg)%n_patch_cells_g,                 &
          &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
          &                     nlev_in,prepicon(jg)%atm_in%qr)

        CALL read_netcdf_data (ncid, 'QS', p_patch(jg)%n_patch_cells_g,                 &
          &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
          &                     nlev_in,prepicon(jg)%atm_in%qs)

        ! Check if surface pressure (PS) or its logarithm (LNPS) is provided as input
        IF (p_pe == p_io) THEN
          IF (nf_inq_varid(ncid, 'PS', varid) == nf_noerr) THEN
            psvar = 'PS'
          ELSE IF (nf_inq_varid(ncid, 'LNPS', varid) == nf_noerr) THEN
            psvar = 'LNPS'
          ENDIF
        ENDIF

        IF (msg_level >= 10) THEN
          WRITE(message_text,'(a)') 'surface pressure variable: '//TRIM(psvar)
          CALL message('', TRIM(message_text))
        ENDIF

        CALL read_netcdf_data (ncid, TRIM(psvar), p_patch(jg)%n_patch_cells_g,          &
          &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
          &                     prepicon(jg)%atm_in%psfc)

        CALL read_netcdf_data (ncid, 'GEOSP', p_patch(jg)%n_patch_cells_g,              &
          &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
          &                     prepicon(jg)%atm_in%phi_sfc)

      ENDIF

      IF (i_oper_mode == 2) THEN ! read in additional fields when vertical interpolation 
                                 ! has already been done in IFS2ICON
        CALL read_netcdf_data (ncid, 'pres', p_patch(jg)%n_patch_cells_g,              &
          &                    p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
          &                    nlev_in,prepicon(jg)%atm_in%pres)
      ENDIF

      IF (i_oper_mode >= 2 .AND. l_sfc_in) THEN ! Read also surface data

        CALL read_netcdf_data (ncid, 'SKT', p_patch(jg)%n_patch_cells_g,                &
          &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
          &                     prepicon(jg)%sfc_in%tskin)

        CALL read_netcdf_data (ncid, 'T_SNOW', p_patch(jg)%n_patch_cells_g,             &
          &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
          &                     prepicon(jg)%sfc_in%tsnow)

        CALL read_netcdf_data (ncid, 'W_SNOW', p_patch(jg)%n_patch_cells_g,             &
          &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
          &                     prepicon(jg)%sfc_in%snowweq)

        CALL read_netcdf_data (ncid, 'RHO_SNOW', p_patch(jg)%n_patch_cells_g,           &
          &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
          &                     prepicon(jg)%sfc_in%snowdens)

        CALL read_netcdf_data (ncid, 'W_I', p_patch(jg)%n_patch_cells_g,               &
          &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
          &                     prepicon(jg)%sfc_in%skinres)

        CALL read_netcdf_data (ncid, 'LSM', p_patch(jg)%n_patch_cells_g,               &
          &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
          &                     prepicon(jg)%sfc_in%ls_mask)

        CALL read_netcdf_data (ncid, 'CI', p_patch(jg)%n_patch_cells_g,               &
          &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
          &                     prepicon(jg)%sfc_in%seaice)

        CALL read_netcdf_data (ncid, 'STL1', p_patch(jg)%n_patch_cells_g,               &
          &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
          &                     prepicon(jg)%sfc_in%tsoil(:,:,1))

        CALL read_netcdf_data (ncid, 'STL2', p_patch(jg)%n_patch_cells_g,               &
          &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
          &                     prepicon(jg)%sfc_in%tsoil(:,:,2))

        CALL read_netcdf_data (ncid, 'STL3', p_patch(jg)%n_patch_cells_g,               &
          &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
          &                     prepicon(jg)%sfc_in%tsoil(:,:,3))

        CALL read_netcdf_data (ncid, 'STL4', p_patch(jg)%n_patch_cells_g,               &
          &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
          &                     prepicon(jg)%sfc_in%tsoil(:,:,4))

        CALL read_netcdf_data (ncid, 'SWVL1', p_patch(jg)%n_patch_cells_g,              &
          &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
          &                     prepicon(jg)%sfc_in%soilwater(:,:,1))

        CALL read_netcdf_data (ncid, 'SWVL2', p_patch(jg)%n_patch_cells_g,              &
          &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
          &                     prepicon(jg)%sfc_in%soilwater(:,:,2))

        CALL read_netcdf_data (ncid, 'SWVL3', p_patch(jg)%n_patch_cells_g,              &
          &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
          &                     prepicon(jg)%sfc_in%soilwater(:,:,3))

        CALL read_netcdf_data (ncid, 'SWVL4', p_patch(jg)%n_patch_cells_g,              &
          &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
          &                     prepicon(jg)%sfc_in%soilwater(:,:,4))

      ENDIF

      IF (jg == 1 .AND. i_oper_mode == 3) THEN ! Allocate and read in vertical coordinate tables

        ALLOCATE(vct_a(nlev_in+1), vct_b(nlev_in+1), vct(2*(nlev_in+1)))

        IF(p_pe == p_io) THEN
          CALL nf(nf_inq_varid(ncid, 'hyai', varid))
          CALL nf(nf_get_var_double(ncid, varid, vct_a))

          CALL nf(nf_inq_varid(ncid, 'hybi', varid))
          CALL nf(nf_get_var_double(ncid, varid, vct_b))
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
          CALL message('', TRIM(message_text))

          DO jk = 1, nlev_in+1
            WRITE(message_text,'(a,i4,F12.4,F12.8)') 'jk, vct_a, vct_b: ',jk, vct_a(jk), vct_b(jk)
            CALL message('', TRIM(message_text))
          ENDDO
        ENDIF

      ENDIF

      ! close file
      !
      IF(p_pe == p_io .AND. i_oper_mode >= 2) CALL nf(nf_close(ncid))


      IF (n_iter_smooth_topo > 0) THEN

        ! The topography smoothing is supposed to be moved into extpar later on
        CALL smooth_topography (p_patch(jg), p_int_state(jg), prepicon(jg)%topography_c, &
                                prepicon(jg)%topography_v)

        WRITE(message_text,'(a,i4)') 'number of topography smoothing steps: ',n_iter_smooth_topo
        CALL message('', TRIM(message_text))

      ENDIF

    ENDDO ! loop over model domains

    IF (n_dom > 1) CALL topo_blending_and_fbk(p_int_state, p_grf_state, prepicon, 1)


  END SUBROUTINE init_prepicon


  RECURSIVE SUBROUTINE topo_blending_and_fbk(p_int, p_grf, prepicon, jg)

    TYPE(t_int_state),     TARGET, INTENT(IN) :: p_int(:)
    TYPE(t_gridref_state), TARGET, INTENT(IN) :: p_grf(:)

    TYPE(t_prepicon_state), INTENT(INOUT) :: prepicon(:)

    INTEGER, INTENT(IN) :: jg

    INTEGER :: jgc, jn


    ! Loop over nested domains
    DO jn = 1, p_patch(jg)%n_childdom

      jgc = p_patch(jg)%child_id(jn)


      CALL topography_blending(p_patch(jg), p_patch(jgc), p_int(jg),  &
               p_int(jgc), p_grf(jg)%p_dom(jn), jn,                   &
               prepicon(jg)%topography_c, prepicon(jgc)%topography_c, &
               prepicon(jgc)%topography_v                             )

      IF (p_patch(jgc)%n_childdom > 0) &
        CALL topo_blending_and_fbk(p_int, p_grf, prepicon, jgc)

    ENDDO

    DO jn = 1, p_patch(jg)%n_childdom

      jgc = p_patch(jg)%child_id(jn)

      IF (lfeedback(jgc)) THEN
        CALL topography_feedback(p_patch(jg), p_int(jg), p_grf(jg), jn, &
          prepicon(jg)%topography_c, prepicon(jgc)%topography_c,        &
          prepicon(jg)%topography_v                                     )
      ENDIF

    ENDDO


  END SUBROUTINE topo_blending_and_fbk

  SUBROUTINE convert_variables(p_int, prepicon)

    TYPE(t_int_state),     TARGET, INTENT(IN) :: p_int(:)

    TYPE(t_prepicon_state), INTENT(INOUT) :: prepicon(:)

    INTEGER :: jg
    REAL(wp), ALLOCATABLE :: temp_v(:,:,:)

    DO jg = 1, n_dom

      IF (i_oper_mode == 2) nlev_in = p_patch(jg)%nlev

      ALLOCATE(temp_v(nproma,SIZE(prepicon(jg)%atm_in%temp,2),p_patch(jg)%nblks_c))

      CALL virtual_temp(p_patch(jg), prepicon(jg)%atm_in%temp, prepicon(jg)%atm_in%qv, &
                        prepicon(jg)%atm_in%qc, prepicon(jg)%atm_in%qi,                &
                        prepicon(jg)%atm_in%qr, prepicon(jg)%atm_in%qs, temp_v  )

      ! Convert thermodynamic variables into set of NH prognostic variables
      CALL convert_thdvars(p_patch(jg), prepicon(jg)%atm_in%pres, temp_v, &
                           prepicon(jg)%atm%rho, prepicon(jg)%atm%exner,  &
                           prepicon(jg)%atm%theta_v                       )

      ! Convert u and v on cell points to vn at edge points
      CALL interp_uv_2_vn(p_patch(jg), p_int(jg), prepicon(jg)%atm_in%u,  &
                          prepicon(jg)%atm_in%v, prepicon(jg)%atm%vn)

      ! Initialize vertical wind field
      CALL init_w(p_patch(jg), p_int(jg), prepicon(jg)%atm%vn, &
                  prepicon(jg)%z_ifc, prepicon(jg)%atm%w)

      ! Finally, copy u, v, and moisture variables from atm_in to atm
!$OMP PARALLEL WORKSHARE
      prepicon(jg)%atm%u(:,:,:)  = prepicon(jg)%atm_in%u(:,:,:)
      prepicon(jg)%atm%v(:,:,:)  = prepicon(jg)%atm_in%v(:,:,:)
      prepicon(jg)%atm%qv(:,:,:) = prepicon(jg)%atm_in%qv(:,:,:)
      prepicon(jg)%atm%qc(:,:,:) = prepicon(jg)%atm_in%qc(:,:,:)
      prepicon(jg)%atm%qi(:,:,:) = prepicon(jg)%atm_in%qi(:,:,:)
      prepicon(jg)%atm%qr(:,:,:) = prepicon(jg)%atm_in%qr(:,:,:)
      prepicon(jg)%atm%qs(:,:,:) = prepicon(jg)%atm_in%qs(:,:,:)
!$OMP end PARALLEL WORKSHARE

      DEALLOCATE(temp_v)

    ENDDO

  END SUBROUTINE convert_variables

  SUBROUTINE compute_coord_fields(p_int, prepicon)

    TYPE(t_int_state),  INTENT(IN)       :: p_int(:)
    TYPE(t_prepicon_state), INTENT(INOUT):: prepicon(:)

    INTEGER :: jg, jgp, nblks_c, npromz_c, nblks_v, npromz_v, nlev, nlevp1
    INTEGER :: i_nchdom
    INTEGER :: nshift_total(n_dom)       !< Total shift of model top w.r.t. global domain
    LOGICAL :: l_half_lev_centr

    !------------------------------------------------------------------------

    SELECT CASE (global_cell_type)
    CASE (6)
      l_half_lev_centr = .TRUE.
      ! The HALF LEVEL where the model layer are flat, moves one layer upward.
      ! there could also be a zero there
      nflat = nflat-1
    CASE DEFAULT
      l_half_lev_centr = .FALSE.
    END SELECT

    nshift_total(1) = 0

    DO jg = 1,n_dom

      nblks_c   = p_patch(jg)%nblks_c
      npromz_c  = p_patch(jg)%npromz_c

      i_nchdom   = MAX(1,p_patch(jg)%n_childdom)

      ! number of vertical levels
      nlev   = p_patch(jg)%nlev
      nlevp1 = p_patch(jg)%nlevp1

      ! total shift of model top with respect to global domain
      IF (jg > 1) THEN
        jgp = p_patch(jg)%parent_id
        nshift_total(jg) = nshift_total(jgp) + p_patch(jg)%nshift
        nflatlev(jg)     = nflatlev(1) - nshift_total(jg)
      ENDIF
      IF (global_cell_type == 6) nflatlev(jg) = nflat
      IF (jg > 1 .AND. nshift_total(jg) > 0 .AND. nflatlev(jg) < 1) THEN
        CALL finish ('mo_prepicon_utils:compute_coord_fields', &
                     'nflat must be higher than the top of the innermost nested domain')
      ENDIF

      ! Compute smooth topography when SLEVE coordinate is used
      IF (ivctype == 2) THEN
        CALL compute_smooth_topo(p_patch(jg), p_int(jg),            &
          prepicon(jg)%topography_c, prepicon(jg)%topography_c_smt, &
          prepicon(jg)%topography_v, prepicon(jg)%topography_v_smt  )
      ENDIF

      ! Compute 3D coordinate fields for cell points (for vertex points needed
      ! only temporarily in set_nh_metrics)
      CALL init_vert_coord(prepicon(jg)%topography_c, prepicon(jg)%topography_c_smt, &
                           prepicon(jg)%z_ifc, prepicon(jg)%z_mc,                    &
                           nlev, nblks_c, npromz_c, nshift_total(jg), nflatlev(jg),  &
                           l_half_lev_centr)


    ENDDO


  END SUBROUTINE compute_coord_fields



  SUBROUTINE init_topo_output_files

    INTEGER :: jg, jlev
    INTEGER :: nlev              !< number of full levels
    CHARACTER(LEN=filename_max) :: gridtype, outputfile
    CHARACTER(LEN=filename_max) :: gridfile(max_dom)

    gridtype='icon'

    DO jg = 1, n_dom

      jlev = p_patch(jg)%level
      nlev = p_patch(jg)%nlev

      WRITE (gridfile(jg),'(a,a,i0,2(a,i2.2),a)') &
        &    TRIM(gridtype),'R',nroot,'B',jlev,'_DOM',jg,'-grid.nc'

      CALL setup_prepicon_vlist( TRIM(p_patch(jg)%grid_filename), jg )

      WRITE (outputfile,'(a,i2.2,a,i0,a,i2.2,a,i0,a)')  &
        'ICON_coord_DOM', jg, '_R', nroot, 'B', jlev, 'L', nlev, '.nc'

      IF(p_pe == p_io) CALL open_output_vlist(TRIM(outputfile), jg)

    ENDDO

  END SUBROUTINE init_topo_output_files

  SUBROUTINE init_atmo_output_files

    INTEGER :: jg, jlev
    INTEGER :: nlev              !< number of full levels
    CHARACTER(LEN=filename_max) :: gridtype, outputfile
    CHARACTER(LEN=filename_max) :: gridfile(max_dom)

    gridtype='icon'

    DO jg = 1, n_dom

      jlev = p_patch(jg)%level
      nlev = p_patch(jg)%nlev

      WRITE (gridfile(jg),'(a,a,i0,2(a,i2.2),a)') &
        &    TRIM(gridtype),'R',nroot,'B',jlev,'_DOM',jg,'-grid.nc'

      CALL setup_prepicon_vlist( TRIM(p_patch(jg)%grid_filename), jg )

      WRITE (outputfile,'(a,i2.2,a,i0,a,i2.2,a,i0,a)')  &
        'ICON_initdata_DOM', jg, '_R', nroot, 'B', jlev, 'L', nlev, '.nc'

      IF(p_pe == p_io) CALL open_output_vlist(TRIM(outputfile), jg)

    ENDDO

  END SUBROUTINE init_atmo_output_files


  !-------------
  !>
  !! SUBROUTINE allocate_prepicon
  !! Allocates the components of the prepicon data type
  !!
  !! @par Revision History
  !! Initial version by Guenther Zaengl, DWD(2011-07-14)
  !!
  !!
  SUBROUTINE allocate_prepicon (prepicon)

    TYPE(t_prepicon_state), INTENT(INOUT) :: prepicon(:)

    ! Local variables: loop control and dimensions
    INTEGER :: jg, nlev, nlevp1, nblks_c, nblks_v, nblks_e

!-------------------------------------------------------------------------

    ! Loop over model domains
    DO jg = 1, n_dom

      nlev = p_patch(jg)%nlev
      nlevp1 = nlev + 1
      nblks_c = p_patch(jg)%nblks_c
      nblks_v = p_patch(jg)%nblks_v
      nblks_e = p_patch(jg)%nblks_e

      IF (i_oper_mode == 2) nlev_in = p_patch(jg)%nlev

      ! basic prep_icon data
      ALLOCATE(prepicon(jg)%topography_c    (nproma,nblks_c),        &
               prepicon(jg)%topography_c_smt(nproma,nblks_c) ,       &
               prepicon(jg)%topography_v    (nproma,nblks_v),        &
               prepicon(jg)%topography_v_smt(nproma,nblks_v) ,       &
               prepicon(jg)%z_ifc           (nproma,nlevp1,nblks_c), &
               prepicon(jg)%z_mc            (nproma,nlev  ,nblks_c) )

      IF (i_oper_mode >= 2) THEN
        ! Allocate atmospheric input data
        ALLOCATE(prepicon(jg)%atm_in%psfc(nproma,         nblks_c), &
                 prepicon(jg)%atm_in%phi_sfc(nproma,      nblks_c), &
                 prepicon(jg)%atm_in%pres (nproma,nlev_in,nblks_c), &
                 prepicon(jg)%atm_in%z3d  (nproma,nlev_in,nblks_c), &
                 prepicon(jg)%atm_in%temp (nproma,nlev_in,nblks_c), &
                 prepicon(jg)%atm_in%u    (nproma,nlev_in,nblks_c), &
                 prepicon(jg)%atm_in%v    (nproma,nlev_in,nblks_c), &
                 prepicon(jg)%atm_in%w    (nproma,nlev_in,nblks_c), &
                 prepicon(jg)%atm_in%omega(nproma,nlev_in,nblks_c), &
                 prepicon(jg)%atm_in%qv   (nproma,nlev_in,nblks_c), &
                 prepicon(jg)%atm_in%qc   (nproma,nlev_in,nblks_c), &
                 prepicon(jg)%atm_in%qi   (nproma,nlev_in,nblks_c), &
                 prepicon(jg)%atm_in%qr   (nproma,nlev_in,nblks_c), &
                 prepicon(jg)%atm_in%qs   (nproma,nlev_in,nblks_c)  )

        ! Allocate surface input data
        ALLOCATE(prepicon(jg)%sfc_in%tskin    (nproma,nblks_c            ), &
                 prepicon(jg)%sfc_in%tsnow    (nproma,nblks_c            ), &
                 prepicon(jg)%sfc_in%snowweq  (nproma,nblks_c            ), &
                 prepicon(jg)%sfc_in%snowdens (nproma,nblks_c            ), &
                 prepicon(jg)%sfc_in%skinres  (nproma,nblks_c            ), &
                 prepicon(jg)%sfc_in%ls_mask  (nproma,nblks_c            ), &
                 prepicon(jg)%sfc_in%seaice   (nproma,nblks_c            ), &
                 prepicon(jg)%sfc_in%tsoil    (nproma,nblks_c,nlevsoil_in), &
                 prepicon(jg)%sfc_in%soilwater(nproma,nblks_c,nlevsoil_in)  )

        ! Allocate atmospheric output data
        ALLOCATE(prepicon(jg)%atm%vn        (nproma,nlev  ,nblks_e), &
                 prepicon(jg)%atm%u         (nproma,nlev  ,nblks_c), &
                 prepicon(jg)%atm%v         (nproma,nlev  ,nblks_c), &
                 prepicon(jg)%atm%w         (nproma,nlevp1,nblks_c), &
                 prepicon(jg)%atm%temp      (nproma,nlev  ,nblks_c), &
                 prepicon(jg)%atm%exner     (nproma,nlev  ,nblks_c), &
                 prepicon(jg)%atm%pres      (nproma,nlev  ,nblks_c), &
                 prepicon(jg)%atm%rho       (nproma,nlev  ,nblks_c), &
                 prepicon(jg)%atm%theta_v   (nproma,nlev  ,nblks_c), &
                 prepicon(jg)%atm%qv        (nproma,nlev  ,nblks_c), &
                 prepicon(jg)%atm%qc        (nproma,nlev  ,nblks_c), &
                 prepicon(jg)%atm%qi        (nproma,nlev  ,nblks_c), &
                 prepicon(jg)%atm%qr        (nproma,nlev  ,nblks_c), &
                 prepicon(jg)%atm%qs        (nproma,nlev  ,nblks_c)  )

        ! Allocate surface output data
        ALLOCATE(prepicon(jg)%sfc%tskin    (nproma,nblks_c          ), &
                 prepicon(jg)%sfc%tsnow    (nproma,nblks_c          ), &
                 prepicon(jg)%sfc%snowweq  (nproma,nblks_c          ), &
                 prepicon(jg)%sfc%snowdens (nproma,nblks_c          ), &
                 prepicon(jg)%sfc%skinres  (nproma,nblks_c          ), &
                 prepicon(jg)%sfc%ls_mask  (nproma,nblks_c          ), &
                 prepicon(jg)%sfc%seaice   (nproma,nblks_c          ), &
                 prepicon(jg)%sfc%tsoil    (nproma,nblks_c,nlev_soil), &
                 prepicon(jg)%sfc%soilwater(nproma,nblks_c,nlev_soil)  )

      ENDIF

      IF (i_oper_mode >= 2 .AND. l_zp_out) THEN
        ! Allocate fields for diagostic output
        ALLOCATE(prepicon(jg)%plev%levels(nzplev),                    &
                 prepicon(jg)%zlev%levels(nzplev),                    &
                 prepicon(jg)%plev%u         (nproma,nzplev,nblks_c), &
                 prepicon(jg)%plev%v         (nproma,nzplev,nblks_c), &
                 prepicon(jg)%plev%temp      (nproma,nzplev,nblks_c), &
                 prepicon(jg)%plev%pres      (nproma,nzplev,nblks_c), &
                 prepicon(jg)%plev%z3d       (nproma,nzplev,nblks_c), &
                 prepicon(jg)%plev%qv        (nproma,nzplev,nblks_c), &
                 prepicon(jg)%zlev%u         (nproma,nzplev,nblks_c), &
                 prepicon(jg)%zlev%v         (nproma,nzplev,nblks_c), &
                 prepicon(jg)%zlev%temp      (nproma,nzplev,nblks_c), &
                 prepicon(jg)%zlev%pres      (nproma,nzplev,nblks_c), &
                 prepicon(jg)%zlev%z3d       (nproma,nzplev,nblks_c), &
                 prepicon(jg)%zlev%qv        (nproma,nzplev,nblks_c)  )
      ENDIF

    ENDDO ! loop over model domains

  END SUBROUTINE allocate_prepicon

  !-------------
  !>
  !! SUBROUTINE deallocate_prepicon
  !! Deallocates the components of the prepicon data type
  !!
  !! @par Revision History
  !! Initial version by Guenther Zaengl, DWD(2011-07-14)
  !!
  !!
  SUBROUTINE deallocate_prepicon (prepicon)

    TYPE(t_prepicon_state), INTENT(INOUT) :: prepicon(:)

    ! Local variables: loop control
    INTEGER :: jg

!-------------------------------------------------------------------------

    ! Loop over model domains
    DO jg = 1, n_dom

      ! basic prep_icon data
      DEALLOCATE(prepicon(jg)%topography_c,   &
               prepicon(jg)%topography_c_smt, &
               prepicon(jg)%topography_v,     &
               prepicon(jg)%topography_v_smt, &
               prepicon(jg)%z_ifc,            &
               prepicon(jg)%z_mc              )

      IF (i_oper_mode >= 2) THEN
        ! atmospheric input data
        DEALLOCATE(prepicon(jg)%atm_in%psfc,  &
                 prepicon(jg)%atm_in%phi_sfc, &
                 prepicon(jg)%atm_in%pres,    &
                 prepicon(jg)%atm_in%z3d,     &
                 prepicon(jg)%atm_in%temp,    &
                 prepicon(jg)%atm_in%u,       &
                 prepicon(jg)%atm_in%v,       &
                 prepicon(jg)%atm_in%w,       &
                 prepicon(jg)%atm_in%omega,   &
                 prepicon(jg)%atm_in%qv,      &
                 prepicon(jg)%atm_in%qc,      &
                 prepicon(jg)%atm_in%qi,      &
                 prepicon(jg)%atm_in%qr,      &
                 prepicon(jg)%atm_in%qs )

        ! surface input data
        DEALLOCATE(prepicon(jg)%sfc_in%tskin,  &
                 prepicon(jg)%sfc_in%tsnow,    &
                 prepicon(jg)%sfc_in%snowweq,  &
                 prepicon(jg)%sfc_in%snowdens, &
                 prepicon(jg)%sfc_in%skinres,  &
                 prepicon(jg)%sfc_in%ls_mask,  &
                 prepicon(jg)%sfc_in%seaice,   &
                 prepicon(jg)%sfc_in%tsoil,    &
                 prepicon(jg)%sfc_in%soilwater )

        ! atmospheric output data
        DEALLOCATE(prepicon(jg)%atm%vn,    &
                 prepicon(jg)%atm%u,       &
                 prepicon(jg)%atm%v,       &
                 prepicon(jg)%atm%w,       &
                 prepicon(jg)%atm%temp,    &
                 prepicon(jg)%atm%exner,   &
                 prepicon(jg)%atm%pres,    &
                 prepicon(jg)%atm%rho,     &
                 prepicon(jg)%atm%theta_v, &
                 prepicon(jg)%atm%qv,      &
                 prepicon(jg)%atm%qc,      &
                 prepicon(jg)%atm%qi,      &
                 prepicon(jg)%atm%qr,      &
                 prepicon(jg)%atm%qs       )

        ! surface output data
        DEALLOCATE(prepicon(jg)%sfc%tskin,  &
                 prepicon(jg)%sfc%tsnow,    &
                 prepicon(jg)%sfc%snowweq,  &
                 prepicon(jg)%sfc%snowdens, &
                 prepicon(jg)%sfc%skinres,  &
                 prepicon(jg)%sfc%ls_mask,  &
                 prepicon(jg)%sfc%seaice,   &
                 prepicon(jg)%sfc%tsoil,    &
                 prepicon(jg)%sfc%soilwater )

      ENDIF

      IF (i_oper_mode >= 2 .AND. l_zp_out) THEN
        !  fields for diagostic output
        DEALLOCATE(prepicon(jg)%plev%levels, &
                 prepicon(jg)%zlev%levels,   &
                 prepicon(jg)%plev%u,        &
                 prepicon(jg)%plev%v,        &
                 prepicon(jg)%plev%temp,     &
                 prepicon(jg)%plev%pres,     &
                 prepicon(jg)%plev%z3d,      &
                 prepicon(jg)%plev%qv,       &
                 prepicon(jg)%zlev%u,        &
                 prepicon(jg)%zlev%v,        &
                 prepicon(jg)%zlev%temp,     &
                 prepicon(jg)%zlev%pres,     &
                 prepicon(jg)%zlev%z3d,      &
                 prepicon(jg)%zlev%qv        )
      ENDIF

    ENDDO ! loop over model domains

  END SUBROUTINE deallocate_prepicon


  SUBROUTINE setup_prepicon_vlist(grid_filename, k_jg)

    CHARACTER(len=*), INTENT(in) :: grid_filename
    INTEGER, INTENT(in) :: k_jg

    INTEGER :: ncid, dimid, varid
    INTEGER :: i_nc, i_ne, i_nv
    INTEGER :: i_ncb, i_neb, i_nvb
    INTEGER :: lnlen, ulen
    INTEGER :: nlev, nlevp1

    REAL(wp), ALLOCATABLE :: clon(:), clat(:), clonv(:), clatv(:)
    REAL(wp), ALLOCATABLE :: elon(:), elat(:), elonv(:), elatv(:)
    REAL(wp), ALLOCATABLE :: vlon(:), vlat(:), vlonv(:), vlatv(:)

    REAL(wp), ALLOCATABLE :: levels(:)

    CHARACTER(len=11) :: name
    CHARACTER(len=12) :: qname
    CHARACTER(len=10) :: dbgname
    CHARACTER(len=3)  :: cjt
    CHARACTER(LEN=1)  :: ctracer
    CHARACTER(LEN=1)  :: anextra ! number of debug fields
    CHARACTER(len=NF_MAX_NAME) :: long_name, units
    INTEGER :: i, jt
    INTEGER :: ivar
    INTEGER :: gridid, zaxisid
    INTEGER :: elemid, tableid

    CHARACTER(len=NF_MAX_NAME) :: att_txt
    INTEGER                    :: astatus

    !=========================================================================
    ! horizontal grids
    !
    CALL nf(nf_open(TRIM(grid_filename), NF_NOWRITE, ncid))
    !
    SELECT CASE (global_cell_type)
    CASE (3)
      CALL nf(nf_inq_dimid(ncid, 'cell', dimid))
    CASE (6)
      CALL nf(nf_inq_dimid(ncid, 'vertex', dimid))
    END SELECT
    CALL nf(nf_inq_dimlen(ncid, dimid, i_nc))
    !
    CALL nf(nf_inq_dimid(ncid, 'edge', dimid))
    CALL nf(nf_inq_dimlen(ncid, dimid, i_ne))
    !
    SELECT CASE (global_cell_type)
    CASE (3)
      CALL nf(nf_inq_dimid(ncid, 'vertex', dimid))
    CASE (6)
      CALL nf(nf_inq_dimid(ncid, 'cell', dimid))
    END SELECT
    CALL nf(nf_inq_dimlen(ncid, dimid, i_nv))
    !
    i_ncb = global_cell_type*i_nc
    i_neb = 4*i_ne
    i_nvb = (9-global_cell_type)*i_nv
    !
    ALLOCATE(clon(i_nc), clat(i_nc), clonv(i_ncb), clatv(i_ncb))
    ALLOCATE(elon(i_ne), elat(i_ne), elonv(i_neb), elatv(i_neb))
    ALLOCATE(vlon(i_nv), vlat(i_nv), vlonv(i_nvb), vlatv(i_nvb))
    !
    !-------------------------------------------------------------------------
    ! cell grid
    !
    gridCellID(k_jg) = gridCreate(GRID_UNSTRUCTURED, i_nc)
    CALL gridDefNvertex(gridCellID(k_jg), global_cell_type)
    !
    SELECT CASE (global_cell_type)
    CASE (3)
      name = 'clon'
    CASE (6)
      name = 'vlon'
    END SELECT
    CALL nf(nf_inq_varid(ncid, name, varid))
    CALL nf(nf_get_var_double(ncid, varid, clon))
    CALL nf(nf_get_att_text(ncid, varid, 'long_name', long_name))
    CALL nf(nf_inq_attlen(ncid, varid, 'long_name', lnlen))
    CALL nf(nf_get_att_text(ncid, varid, 'units', units))
    CALL nf(nf_inq_attlen(ncid, varid, 'units', ulen))
    !
    CALL gridDefXname(gridCellID(k_jg), name)
    CALL gridDefXvals(gridCellID(k_jg), clon)
    CALL gridDefXlongname(gridCellID(k_jg), long_name(1:lnlen))
    CALL gridDefXunits(gridCellID(k_jg), units(1:ulen))
    !
    SELECT CASE (global_cell_type)
    CASE (3)
      name = 'clat'
    CASE (6)
      name = 'vlat'
    END SELECT
    CALL nf(nf_inq_varid(ncid, name, varid))
    CALL nf(nf_get_var_double(ncid, varid, clat))
    CALL nf(nf_get_att_text(ncid, varid, 'long_name', long_name))
    CALL nf(nf_inq_attlen(ncid, varid, 'long_name', lnlen))
    CALL nf(nf_get_att_text(ncid, varid, 'units', units))
    CALL nf(nf_inq_attlen(ncid, varid, 'units', ulen))
    !
    CALL gridDefYname(gridCellID(k_jg), name)
    CALL gridDefYvals(gridCellID(k_jg), clat)
    CALL gridDefYlongname(gridCellID(k_jg), long_name(1:lnlen))
    CALL gridDefYunits(gridCellID(k_jg), units(1:ulen))
    !
    SELECT CASE (global_cell_type)
    CASE (3)
      CALL nf(nf_inq_varid(ncid, 'clon_vertices', varid))
    CASE (6)
      CALL nf(nf_inq_varid(ncid, 'vlon_vertices', varid))
    END SELECT
    CALL nf(nf_get_var_double(ncid, varid, clonv))
    !
    CALL gridDefXbounds(gridCellID(k_jg), clonv)
    !
    SELECT CASE (global_cell_type)
    CASE (3)
      CALL nf(nf_inq_varid(ncid, 'clat_vertices', varid))
    CASE (6)
      CALL nf(nf_inq_varid(ncid, 'vlat_vertices', varid))
    END SELECT
    CALL nf(nf_get_var_double(ncid, varid, clatv))
    !
    CALL gridDefYbounds(gridCellID(k_jg), clatv)
    !
    !-------------------------------------------------------------------------
    ! edge grid
    !
    gridEdgeID(k_jg) = gridCreate(GRID_UNSTRUCTURED, i_ne)
    CALL gridDefNvertex(gridEdgeID(k_jg), 4)
    !
    name = 'elon'
    CALL nf(nf_inq_varid(ncid, name, varid))
    CALL nf(nf_get_var_double(ncid, varid, elon))
    CALL nf(nf_get_att_text(ncid, varid, 'long_name', long_name))
    CALL nf(nf_inq_attlen(ncid, varid, 'long_name', lnlen))
    CALL nf(nf_get_att_text(ncid, varid, 'units', units))
    CALL nf(nf_inq_attlen(ncid, varid, 'units', ulen))
    !
    CALL gridDefXname(gridEdgeID(k_jg), name)
    CALL gridDefXvals(gridEdgeID(k_jg), elon)
    CALL gridDefXlongname(gridEdgeID(k_jg), long_name(1:lnlen))
    CALL gridDefXunits(gridEdgeID(k_jg), units(1:ulen))
    !
    name = 'elat'
    CALL nf(nf_inq_varid(ncid, name, varid))
    CALL nf(nf_get_var_double(ncid, varid, elat))
    CALL nf(nf_get_att_text(ncid, varid, 'long_name', long_name))
    CALL nf(nf_inq_attlen(ncid, varid, 'long_name', lnlen))
    CALL nf(nf_get_att_text(ncid, varid, 'units', units))
    CALL nf(nf_inq_attlen(ncid, varid, 'units', ulen))
    !
    CALL gridDefYname(gridEdgeID(k_jg), name)
    CALL gridDefYvals(gridEdgeID(k_jg), elat)
    CALL gridDefYlongname(gridEdgeID(k_jg), long_name(1:lnlen))
    CALL gridDefYunits(gridEdgeID(k_jg), units(1:ulen))
    !
    CALL nf(nf_inq_varid(ncid, 'elon_vertices', varid))
    CALL nf(nf_get_var_double(ncid, varid, elonv))
    !
    CALL gridDefXbounds(gridEdgeID(k_jg), elonv)
    !
    CALL nf(nf_inq_varid(ncid, 'elat_vertices', varid))
    CALL nf(nf_get_var_double(ncid, varid, elatv))
    !
    CALL gridDefYbounds(gridEdgeID(k_jg), elatv)
    !
    !-------------------------------------------------------------------------
    ! vertex grid
    !
    gridVertexID(k_jg) = gridCreate(GRID_UNSTRUCTURED, i_nv)
    CALL gridDefNvertex(gridVertexID(k_jg), 9-global_cell_type)
    !
    SELECT CASE (global_cell_type)
    CASE (3)
      name = 'vlon'
    CASE (6)
      name = 'clon'
    END SELECT
    CALL nf(nf_inq_varid(ncid, name, varid))
    CALL nf(nf_get_var_double(ncid, varid, vlon))
    CALL nf(nf_get_att_text(ncid, varid, 'long_name', long_name))
    CALL nf(nf_inq_attlen(ncid, varid, 'long_name', lnlen))
    CALL nf(nf_get_att_text(ncid, varid, 'units', units))
    CALL nf(nf_inq_attlen(ncid, varid, 'units', ulen))
    !
    CALL gridDefXname(gridVertexID(k_jg), name)
    CALL gridDefXvals(gridVertexID(k_jg), vlon)
    CALL gridDefXlongname(gridVertexID(k_jg), long_name(1:lnlen))
    CALL gridDefXunits(gridVertexID(k_jg), units(1:ulen))
    !
    SELECT CASE (global_cell_type)
    CASE (3)
      name = 'vlat'
    CASE (6)
      name = 'clat'
    END SELECT
    CALL nf(nf_inq_varid(ncid, name , varid))
    CALL nf(nf_get_var_double(ncid, varid, vlat))
    CALL nf(nf_get_att_text(ncid, varid, 'long_name', long_name))
    CALL nf(nf_inq_attlen(ncid, varid, 'long_name', lnlen))
    CALL nf(nf_get_att_text(ncid, varid, 'units', units))
    CALL nf(nf_inq_attlen(ncid, varid, 'units', ulen))
    !
    CALL gridDefYname(gridVertexID(k_jg), name)
    CALL gridDefYvals(gridVertexID(k_jg), vlat)
    CALL gridDefYlongname(gridVertexID(k_jg), long_name(1:lnlen))
    CALL gridDefYunits(gridVertexID(k_jg), units(1:ulen))
    !
    IF(global_cell_type==3) THEN
      CALL nf(nf_inq_varid(ncid, 'vlon_vertices', varid))
    ELSEIF(global_cell_type==6) THEN
      CALL nf(nf_inq_varid(ncid, 'clon_vertices', varid))
    ENDIF
    CALL nf(nf_get_var_double(ncid, varid, vlonv))
    !
    CALL gridDefXbounds(gridVertexID(k_jg), vlonv)
    !
    IF(global_cell_type==3) THEN
      CALL nf(nf_inq_varid(ncid, 'vlat_vertices', varid))
    ELSEIF(global_cell_type==6) THEN
      CALL nf(nf_inq_varid(ncid, 'clat_vertices', varid))
    ENDIF
    CALL nf(nf_get_var_double(ncid, varid, vlatv))
    !
    CALL gridDefYbounds(gridVertexID(k_jg), vlatv)
    !
    !-------------------------------------------------------------------------
    !
    DEALLOCATE(clon, clat, clonv, clatv)
    DEALLOCATE(elon, elat, elonv, elatv)
    DEALLOCATE(vlon, vlat, vlonv, vlatv)
    !
    !=========================================================================
    ! vertical grids
    !
    nlev   = num_lev(k_jg)
    nlevp1 = num_levp1(k_jg)

    zaxisID_surface(k_jg) = zaxisCreate(ZAXIS_SURFACE, 1)
    !
    ALLOCATE(levels(1))
    levels(1) = 0.0_wp
    CALL zaxisDefLevels(zaxisID_surface(k_jg), levels)
    DEALLOCATE(levels)
    !
    zaxisID_hybrid(k_jg)  = zaxisCreate(ZAXIS_HYBRID, nlev)
    zaxisID_hybrid_half(k_jg)  = zaxisCreate(ZAXIS_HYBRID_HALF, nlevp1)
    !
    ALLOCATE(levels(nlev))
    DO i = 1, nlev
      levels(i) = REAL(i,wp)
    END DO
    CALL zaxisDefLevels(zaxisID_hybrid(k_jg), levels)
    DEALLOCATE(levels)
    CALL zaxisDefVct(zaxisID_hybrid(k_jg), 2*nlevp1, vct(1:2*nlevp1))
    !
    ALLOCATE(levels(nlevp1))
    DO i = 1, nlevp1
      levels(i) = REAL(i,wp)
    END DO
    CALL zaxisDefLevels(zaxisID_hybrid_half(k_jg), levels)
    DEALLOCATE(levels)
    CALL zaxisDefVct(zaxisID_hybrid_half(k_jg), 2*nlevp1, vct(1:2*nlevp1))

   IF (i_oper_mode >= 2 .AND. l_zp_out) THEN
      zaxisID_pres(k_jg) = zaxisCreate(ZAXIS_PRESSURE, nzplev)
      zaxisID_hgt(k_jg)  = zaxisCreate(ZAXIS_HEIGHT, nzplev)
      ALLOCATE(levels(nzplev))
      DO i = 1, nzplev
        levels(i) = REAL(i,wp)
      END DO
      CALL zaxisDefLevels(zaxisID_pres(k_jg), levels)
      CALL zaxisDefLevels(zaxisID_hgt(k_jg), levels)
      DEALLOCATE(levels)
      CALL zaxisDefVct(zaxisID_pres(k_jg), nzplev, prepicon(k_jg)%plev%levels(1:nzplev))
      CALL zaxisDefVct(zaxisID_hgt(k_jg),  nzplev, prepicon(k_jg)%zlev%levels(1:nzplev))
    ENDIF

    !
    !=========================================================================
    ! time dimension
    !
    taxisID(k_jg) = taxisCreate(TAXIS_ABSOLUTE)
    !
    !=========================================================================
    !
    vlistID(k_jg) = vlistCreate()
    !
    !-------------------------------------------------------------------------
    ! register time axis
    !
    CALL vlistDefTaxis(vlistID(k_jg), taxisID(k_jg))
    !
    !-------------------------------------------------------------------------
    ! global attributes
    !
    !
    ! Parameters of /grid_nml/
    ! ------------------------
    CALL addGlobAttInt('nroot',nroot,vlistID(k_jg),astatus)
    CALL addGlobAttInt('start_lev',start_lev,vlistID(k_jg),astatus)
    CALL addGlobAttInt('n_dom',n_dom,vlistID(k_jg),astatus)
    !
    ! Parameters of /run_nml/
    ! -----------------------
    CALL addGlobAttInt('run_nml:global_cell_type',global_cell_type,vlistID(k_jg),astatus)
    CALL addGlobAttInt('run_nml:num_lev',num_lev(k_jg),vlistID(k_jg),astatus)

    CALL addGlobAttInt('run_nml:iequations',iequations,vlistID(k_jg),astatus)
 

    !
    ! Parameters of /nonhydrostatic_nml/
    ! ----------------------------

    IF (iequations == 3) THEN
       CALL addGlobAttInt('nonhydrostatic_nml:ivctype',ivctype,vlistID(k_jg),astatus)

    END IF


    !-------------------------------------------------------------------------
    ! register variables
    !
    varids(:,k_jg)   = 0
    ! initialize total number of varids for domain jg
    num_varids(k_jg) = 0


    CALL addVar(TimeVar('ZF3',&
    &                   'geopotential height',&
    &                   'm', 156, 128,&
    &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_hybrid(k_jg)),&
    &          k_jg)

    CALL addVar(TimeVar('ZH3',&
    &                   'half level height',&
    &                   'm', 253, 128,&
    &                    vlistID(k_jg), gridCellID(k_jg), zaxisID_hybrid_half(k_jg)),&
    &           k_jg)

    CALL addVar(TimeVar('topography_c',&
    &                   'topogaphy at cell points',&
    &                   'm', 001, 001,&
    &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
    &           k_jg)

    CALL addVar(TimeVar('topography_v',&
    &                   'topography at vertex points',&
    &                   'm', 001, 001,&
    &                   vlistID(k_jg), gridVertexID(k_jg),zaxisID_surface(k_jg)),&
    &           k_jg)

    IF (i_oper_mode >= 2) THEN

      CALL addVar(TimeVar('normal_velocity',&
      &                   'velocity normal to edge',&
      &                   'm/s', 254, 128, &
      &                   vlistID(k_jg), gridEdgeID(k_jg),zaxisID_hybrid(k_jg)),&
      &           k_jg)

      CALL addVar(TimeVar('U',&
      &                   'zonal wind',&
      &                   'm/s', 131, 128,&
      &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_hybrid(k_jg)),&
      &           k_jg)

      CALL addVar(TimeVar('V',&
      &                   'meridional wind',&
      &                   'm/s', 132, 128, &
      &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_hybrid(k_jg)),&
      &           k_jg)

      CALL addVar(TimeVar('W',&
      &                   'upward air velocity',&
      &                   'm/s', 40, 2,&
      &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_hybrid_half(k_jg)),&
      &          k_jg)

      CALL addVar(TimeVar('T', &
      &                   'temperature',&
      &                   'K',130, 128, &
      &                   vlistID(k_jg),gridCellID(k_jg),zaxisID_hybrid(k_jg)),&
      &          k_jg)

      CALL addVar(TimeVar('THETA_V', &
      &                   'virtual potential temperature',&
      &                   'K',192, 128, &
      &                   vlistID(k_jg),gridCellID(k_jg),zaxisID_hybrid(k_jg)),&
      &          k_jg)

      CALL addVar(TimeVar('EXNER',&
      &                   'Exner pressure', &
      &                   '-', 193, 128, &
      &                   vlistID(k_jg),gridCellID(k_jg),zaxisID_hybrid(k_jg)),&
      &          k_jg)

      CALL addVar(TimeVar('RHO', &
      &                   'density', &
      &                   'kg/m**3', 194, 128,&
      &                   vlistID(k_jg),gridCellID(k_jg),zaxisID_hybrid(k_jg)), &
      &           k_jg)

      CALL addVar(TimeVar('QV',&
      &                   'specific humidity',&
      &                   '(0-1)', 91, 128,&
      &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_hybrid(k_jg)),&
      &           k_jg)

      CALL addVar(TimeVar('QC',&
      &                   'cloud water',&
      &                   '(0-1)', 92, 128,&
      &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_hybrid(k_jg)),&
      &           k_jg)

      CALL addVar(TimeVar('QI',&
      &                   'cloud ice',&
      &                   '(0-1)', 93, 128,&
      &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_hybrid(k_jg)),&
      &           k_jg)

      CALL addVar(TimeVar('QR',&
      &                   'rain water',&
      &                   '(0-1)', 94, 128,&
      &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_hybrid(k_jg)),&
      &           k_jg)

      CALL addVar(TimeVar('QS',&
      &                   'snow',&
      &                   '(0-1)', 95, 128,&
      &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_hybrid(k_jg)),&
      &           k_jg)

    ENDIF

    IF (i_oper_mode >= 2 .AND. l_sfc_in) THEN

      CALL addVar(TimeVar('TSN',&
      &                   'temperature of snow layer',&
      &                   'K', 238, 128,&
      &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
      &           k_jg)

      CALL addVar(TimeVar('TSK',&
      &                   'skin temperature',&
      &                   'm', 235, 128,&
      &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
      &           k_jg)

      CALL addVar(TimeVar('SNWE',&
      &                   'snow water equivalent',&
      &                   'm', 141, 128,&
      &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
      &           k_jg)

      CALL addVar(TimeVar('SNDENS',&
      &                   'snow density',&
      &                   'kg/m**3', 33, 128,&
      &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
      &           k_jg)

      CALL addVar(TimeVar('SRC',&
      &                   'skin reservoir content',&
      &                   'm', 33, 128,&
      &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
      &           k_jg)

      CALL addVar(TimeVar('SEAICE',&
      &                   'sea-ice cover',&
      &                   '(0-1)', 33, 128,&
      &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
      &           k_jg)

      CALL addVar(TimeVar('LSM',&
      &                   'land-sea mask',&
      &                   '(0-1)', 172, 128,&
      &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
      &           k_jg)

      CALL addVar(TimeVar('TSOIL1',&
      &                   'soil temperature level 1',&
      &                   'K', 139, 128,&
      &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
      &           k_jg)

      CALL addVar(TimeVar('TSOIL2',&
      &                   'soil temperature level 2',&
      &                   'K', 139, 128,&
      &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
      &           k_jg)

      CALL addVar(TimeVar('TSOIL3',&
      &                   'soil temperature level 3',&
      &                   'K', 139, 128,&
      &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
      &           k_jg)

      CALL addVar(TimeVar('TSOIL4',&
      &                   'soil temperature level 4',&
      &                   'K', 139, 128,&
      &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
      &           k_jg)

      CALL addVar(TimeVar('TSOIL5',&
      &                   'soil temperature level 5',&
      &                   'K', 139, 128,&
      &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
      &           k_jg)

      CALL addVar(TimeVar('TSOIL6',&
      &                   'soil temperature level 6',&
      &                   'K', 139, 128,&
      &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
      &           k_jg)

      CALL addVar(TimeVar('TSOIL7',&
      &                   'soil temperature level 7',&
      &                   'K', 139, 128,&
      &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
      &           k_jg)

      CALL addVar(TimeVar('WSOIL1',&
      &                   'soil water content level 1',&
      &                   'm**3/m**3', 39, 128,&
      &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
      &           k_jg)

      CALL addVar(TimeVar('WSOIL2',&
      &                   'soil water content level 2',&
      &                   'm**3/m**3', 39, 128,&
      &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
      &           k_jg)

      CALL addVar(TimeVar('WSOIL3',&
      &                   'soil water content level 3',&
      &                   'm**3/m**3', 39, 128,&
      &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
      &           k_jg)

      CALL addVar(TimeVar('WSOIL4',&
      &                   'soil water content level 4',&
      &                   'm**3/m**3', 39, 128,&
      &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
      &           k_jg)

      CALL addVar(TimeVar('WSOIL5',&
      &                   'soil water content level 5',&
      &                   'm**3/m**3', 39, 128,&
      &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
      &           k_jg)

      CALL addVar(TimeVar('WSOIL6',&
      &                   'soil water content level 6',&
      &                   'm**3/m**3', 39, 128,&
      &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
      &           k_jg)

      CALL addVar(TimeVar('WSOIL7',&
      &                   'soil water content level 7',&
      &                   'm**3/m**3', 39, 128,&
      &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
      &           k_jg)

    ENDIF

    IF (i_oper_mode >= 2 .AND. l_zp_out) THEN
      CALL addVar(TimeVar('UZ',&
      &                   'zonal wind',&
      &                   'm/s', 131, 128,&
      &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_hgt(k_jg)),&
      &           k_jg)

      CALL addVar(TimeVar('VZ',&
      &                   'meridional wind',&
      &                   'm/s', 132, 128, &
      &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_hgt(k_jg)),&
      &           k_jg)

      CALL addVar(TimeVar('TZ', &
      &                   'temperature',&
      &                   'K',130, 128, &
      &                   vlistID(k_jg),gridCellID(k_jg),zaxisID_hgt(k_jg)),&
      &          k_jg)

      CALL addVar(TimeVar('QVZ',&
      &                   'specific humidity',&
      &                   '(0-1)', 91, 128,&
      &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_hgt(k_jg)),&
      &           k_jg)

      CALL addVar(TimeVar('PZ',&
      &                   'pressure',&
      &                   'Pa', 255, 128,&
      &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_hgt(k_jg)),&
      &          k_jg)

      CALL addVar(TimeVar('UP',&
      &                   'zonal wind',&
      &                   'm/s', 131, 128,&
      &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_pres(k_jg)),&
      &           k_jg)

      CALL addVar(TimeVar('VP',&
      &                   'meridional wind',&
      &                   'm/s', 132, 128, &
      &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_pres(k_jg)),&
      &           k_jg)

      CALL addVar(TimeVar('TP', &
      &                   'temperature',&
      &                   'K',130, 128, &
      &                   vlistID(k_jg),gridCellID(k_jg),zaxisID_pres(k_jg)),&
      &          k_jg)

      CALL addVar(TimeVar('QVP',&
      &                   'specific humidity',&
      &                   '(0-1)', 91, 128,&
      &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_pres(k_jg)),&
      &           k_jg)

      CALL addVar(TimeVar('ZP',&
      &                   'geopotential height',&
      &                   'm', 156, 128,&
      &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_pres(k_jg)),&
      &          k_jg)

    ENDIF

    CALL nf(nf_close(ncid))
    !
    !=========================================================================

    ! Create description of all output variables in vlist
    num_output_vars(k_jg) = vlistNvars(vlistID(k_jg))

    DO ivar = 1, num_output_vars(k_jg)

      gridid = vlistInqVarGrid(vlistID(k_jg), varids(ivar, k_jg))
      IF(gridid == gridCellID(k_jg)) THEN
        outvar_desc(ivar, k_jg)%type = GATHER_C
      ELSEIF(gridid == gridEdgeID(k_jg)) THEN
        outvar_desc(ivar, k_jg)%type = GATHER_E
      ELSEIF(gridid == gridVertexID(k_jg)) THEN
        outvar_desc(ivar, k_jg)%type = GATHER_V
      ELSE
        CALL finish('setup_vlist','got illegal gridid')
      ENDIF

      zaxisid = vlistInqVarZaxis(vlistID(k_jg), varids(ivar, k_jg))
      outvar_desc(ivar, k_jg)%nlev = zaxisInqSize(zaxisid)

      CALL vlistInqVarName(vlistID(k_jg), varids(ivar, k_jg), outvar_desc(ivar, k_jg)%name)

    ENDDO

  END SUBROUTINE setup_prepicon_vlist



  SUBROUTINE write_prepicon_output(datetime, z_sim_time)

    !=========================================================================
    !> Write output directly: PE 0 gathers and writes, others send
    !=========================================================================

    TYPE(t_datetime),            INTENT(in) :: datetime
    REAL(wp), OPTIONAL,          INTENT(in) :: z_sim_time(n_dom)
    INTEGER :: idate, itime
    INTEGER :: istatus
    INTEGER :: jg, ivar, n_tot

    REAL(wp), POINTER :: ptr2(:,:)
    REAL(wp), POINTER :: ptr3(:,:,:)
    LOGICAL :: reset, delete

    REAL(wp), ALLOCATABLE :: streamvar1(:), streamvar2(:,:)
    REAL(wp) :: p_sim_time
    
    idate   = cdiEncodeDate(datetime%year, datetime%month, datetime%day)
    itime   = cdiEncodeTime(datetime%hour, datetime%minute, NINT(datetime%second))

    ! Make streamvar1/streamvar2 defined everywhere

    IF(p_pe /= p_io) ALLOCATE(streamvar1(1), streamvar2(1,1))

    DO jg = 1, n_dom

      IF (PRESENT(z_sim_time)) THEN
        p_sim_time = z_sim_time(jg)
      ELSE
        p_sim_time = 1.0_wp
      ENDIF
      
      IF(p_pe == p_io) THEN
        CALL taxisDefVdate(taxisID(jg), idate)   ! YYYYMMDD
        CALL taxisDefVtime(taxisID(jg), itime)   ! HHMM

        istatus = streamDefTimestep(streamID(jg), iostep)
      ENDIF

      DO ivar = 1, num_output_vars(jg)

        CALL get_outvar_ptr_prepicon &
          & (outvar_desc(ivar,jg)%name, jg, p_sim_time, ptr2, ptr3, reset, delete)

        SELECT CASE(outvar_desc(ivar, jg)%type)
          CASE (GATHER_C)
            n_tot = p_patch(jg)%n_patch_cells_g
          CASE (GATHER_E)
            n_tot = p_patch(jg)%n_patch_edges_g
          CASE (GATHER_V)
            n_tot = p_patch(jg)%n_patch_verts_g
          CASE DEFAULT
            CALL finish('write_vlist', 'Illegal type in outvar_desc')
        END SELECT

        klev = outvar_desc(ivar, jg)%nlev

        ! Pack and output variable

        IF(ASSOCIATED(ptr2)) THEN

          IF(p_pe == p_io) ALLOCATE(streamvar1(n_tot))

          CALL gather_array1( outvar_desc(ivar, jg)%type, p_patch(jg), ptr2, &
               &                        streamvar1,outvar_desc(ivar,jg)%name )

          IF(p_pe == p_io) THEN
            CALL streamWriteVar(streamID(jg), varids(ivar,jg), streamvar1, 0)
            DEALLOCATE(streamvar1)
          ENDIF
          IF(reset) ptr2 = 0._wp
          IF(delete) DEALLOCATE(ptr2)

        ELSE

          IF(p_pe == p_io) ALLOCATE(streamvar2(n_tot, klev))
          CALL gather_array2( outvar_desc(ivar, jg)%type, p_patch(jg), ptr3,&
               &                       streamvar2,outvar_desc(ivar,jg)%name )
          IF(p_pe == p_io) THEN
            CALL streamWriteVar(streamID(jg), varids(ivar,jg), streamvar2, 0)
            DEALLOCATE(streamvar2)
          ENDIF
          IF(reset) ptr3 = 0._wp
          IF(delete) DEALLOCATE(ptr3)

        ENDIF

      ENDDO

      IF(p_pe == p_io) THEN
        IF (lkeep_in_sync) THEN
          CALL streamSync(streamID(jg))
        END IF
      END IF

    END DO

    IF(p_pe /= p_io) DEALLOCATE(streamvar1, streamvar2)


  END SUBROUTINE write_prepicon_output



  SUBROUTINE close_prepicon_output_files

    INTEGER jg

    DO jg = n_dom, 1, -1
      IF (p_pe == p_io) CALL close_output_vlist(jg)

      CALL vlistDestroy(vlistID(jg))
      CALL zaxisDestroy(zaxisID_hybrid(jg))
      CALL zaxisDestroy(zaxisID_hybrid_half(jg))
      CALL zaxisDestroy(zaxisID_surface(jg))
      CALL gridDestroy(gridVertexID(jg))
      CALL gridDestroy(gridEdgeID(jg))
      CALL gridDestroy(gridCellID(jg))

      num_output_vars(jg) = 0

    ENDDO

  END SUBROUTINE close_prepicon_output_files


  SUBROUTINE get_outvar_ptr_prepicon(varname, jg, z_sim_time,ptr2, ptr3, reset, delete)

    CHARACTER(LEN=*), INTENT(IN) :: varname
    INTEGER, INTENT(IN) :: jg
    REAL(wp), INTENT(in) :: z_sim_time
    LOGICAL, INTENT(OUT) :: reset, delete

    LOGICAL :: not_found

    REAL(wp), POINTER :: ptr2(:,:)
    REAL(wp), POINTER :: ptr3(:,:,:)
 

    ptr2 => NULL()
    ptr3 => NULL()
    reset  = .FALSE.
    delete = .FALSE.
    not_found = .FALSE.

    SELECT CASE(varname)
      CASE ('ZF3');             ptr3 => prepicon(jg)%z_mc
      CASE ('ZH3');             ptr3 => prepicon(jg)%z_ifc
      CASE ('topography_c');    ptr2 => prepicon(jg)%topography_c
      CASE ('topography_v');    ptr2 => prepicon(jg)%topography_v
      CASE ('normal_velocity'); ptr3 => prepicon(jg)%atm%vn
      CASE ('U');               ptr3 => prepicon(jg)%atm%u
      CASE ('V');               ptr3 => prepicon(jg)%atm%v
      CASE ('W');               ptr3 => prepicon(jg)%atm%w
      CASE ('T');               ptr3 => prepicon(jg)%atm%temp
      CASE ('THETA_V');         ptr3 => prepicon(jg)%atm%theta_v
      CASE ('EXNER');           ptr3 => prepicon(jg)%atm%exner
      CASE ('RHO');             ptr3 => prepicon(jg)%atm%rho
      CASE ('QV');              ptr3 => prepicon(jg)%atm%qv
      CASE ('QC');              ptr3 => prepicon(jg)%atm%qc
      CASE ('QI');              ptr3 => prepicon(jg)%atm%qi
      CASE ('QR');              ptr3 => prepicon(jg)%atm%qr
      CASE ('QS');              ptr3 => prepicon(jg)%atm%qs

      CASE ('TSN');             ptr2 => prepicon(jg)%sfc%tsnow
      CASE ('TSK');             ptr2 => prepicon(jg)%sfc%tskin
      CASE ('SNWE');            ptr2 => prepicon(jg)%sfc%snowweq
      CASE ('SNDENS');          ptr2 => prepicon(jg)%sfc%snowdens
      CASE ('SRC');             ptr2 => prepicon(jg)%sfc%skinres
      CASE ('SEAICE');          ptr2 => prepicon(jg)%sfc%seaice
      CASE ('LSM');             ptr2 => prepicon(jg)%sfc%ls_mask
      CASE ('TSOIL1');          ptr2 => prepicon(jg)%sfc%tsoil(:,:,1)
      CASE ('TSOIL2');          ptr2 => prepicon(jg)%sfc%tsoil(:,:,2)
      CASE ('TSOIL3');          ptr2 => prepicon(jg)%sfc%tsoil(:,:,3)
      CASE ('TSOIL4');          ptr2 => prepicon(jg)%sfc%tsoil(:,:,4)
      CASE ('TSOIL5');          ptr2 => prepicon(jg)%sfc%tsoil(:,:,5)
      CASE ('TSOIL6');          ptr2 => prepicon(jg)%sfc%tsoil(:,:,6)
      CASE ('TSOIL7');          ptr2 => prepicon(jg)%sfc%tsoil(:,:,7)
      CASE ('WSOIL1');          ptr2 => prepicon(jg)%sfc%soilwater(:,:,1)
      CASE ('WSOIL2');          ptr2 => prepicon(jg)%sfc%soilwater(:,:,2)
      CASE ('WSOIL3');          ptr2 => prepicon(jg)%sfc%soilwater(:,:,3)
      CASE ('WSOIL4');          ptr2 => prepicon(jg)%sfc%soilwater(:,:,4)
      CASE ('WSOIL5');          ptr2 => prepicon(jg)%sfc%soilwater(:,:,5)
      CASE ('WSOIL6');          ptr2 => prepicon(jg)%sfc%soilwater(:,:,6)
      CASE ('WSOIL7');          ptr2 => prepicon(jg)%sfc%soilwater(:,:,7)

      CASE ('UP');              ptr3 => prepicon(jg)%plev%u
      CASE ('VP');              ptr3 => prepicon(jg)%plev%v
      CASE ('TP');              ptr3 => prepicon(jg)%plev%temp
      CASE ('ZP');              ptr3 => prepicon(jg)%plev%z3d
      CASE ('QVP');             ptr3 => prepicon(jg)%plev%qv
      CASE ('UZ');              ptr3 => prepicon(jg)%zlev%u
      CASE ('VZ');              ptr3 => prepicon(jg)%zlev%v
      CASE ('TZ');              ptr3 => prepicon(jg)%zlev%temp
      CASE ('PZ');              ptr3 => prepicon(jg)%zlev%pres
      CASE ('QVZ');             ptr3 => prepicon(jg)%zlev%qv


      CASE DEFAULT;             not_found = .TRUE.
    END SELECT

    IF (not_found) CALL finish('get_outvar_ptr_prepicon', 'Unkown variable type name: '//varname)

  END SUBROUTINE get_outvar_ptr_prepicon

  SUBROUTINE open_output_vlist(vlist_filename, k_jg)

    CHARACTER(len=*), INTENT(in) :: vlist_filename
    INTEGER, INTENT(in) :: k_jg

    ! Each time a new NetCDF is created, reset "iostep" to zero
    ! (Otherwise we will get an error message from a CDI subroutine.)

    iostep = 0

    !=========================================================================
    ! open file for writing (using netCDF)
    !

    streamID(k_jg) = streamOpenWrite(TRIM(vlist_filename), FILETYPE_NC2)
    IF (streamID(k_jg) < 0) THEN
      CALL finish('setup_vlist', cdiStringError(streamID(k_jg)))
    ENDIF
    !
    CALL streamDefVlist(streamID(k_jg), vlistID(k_jg))
    !
    !=========================================================================

  END SUBROUTINE open_output_vlist

  !-------------------------------------------------------------------------------------------------

  SUBROUTINE close_output_vlist(k_jg)

    INTEGER, INTENT(in) :: k_jg

    !=========================================================================
    ! close file
    !
    CALL streamClose(streamID(k_jg))
    !
    !=========================================================================

  END SUBROUTINE close_output_vlist

  SUBROUTINE addGlobAttInt(att_name, att_int, vlist, astatus)
    CHARACTER(*), INTENT(IN)  :: att_name
    INTEGER     , INTENT(IN)  :: att_int, vlist
    INTEGER     , INTENT(OUT) :: astatus

    astatus = vlistDefAttInt(vlist,CDI_GLOBAL,TRIM(att_name),DATATYPE_INT32,1,att_int)
  END SUBROUTINE addGlobAttInt

  SUBROUTINE addGlobAttTxt(att_name, att_txt, vlist, astatus)
    CHARACTER(*), INTENT(IN)  :: att_name, att_txt
    INTEGER     , INTENT(IN)  :: vlist
    INTEGER     , INTENT(OUT) :: astatus

    astatus = vlistDefAttTxt(vlist,CDI_GLOBAL,TRIM(att_name),LEN(TRIM(att_txt)),TRIM(att_txt))
  END SUBROUTINE addGlobAttTxt

  SUBROUTINE addGlobAttTxtFromLog(att_name, boolian, vlist, astatus)
    CHARACTER(*), INTENT(IN)  :: att_name
    LOGICAL, INTENT(IN)       :: boolian
    INTEGER     , INTENT(IN)  :: vlist
    INTEGER     , INTENT(OUT) :: astatus
    CALL addGlobAttTxt(att_name, TRIM(MERGE('.true. ','.false.',boolian)), vlist, astatus)
  END SUBROUTINE addGlobAttTxtFromLog

  SUBROUTINE addGlobAttFlt(att_name, att_flt, vlist, astatus)
    CHARACTER(*), INTENT(IN)  :: att_name
    REAL(wp)    , INTENT(IN)  :: att_flt
    INTEGER     , INTENT(IN)  :: vlist
    INTEGER     , INTENT(OUT) :: astatus

    astatus = vlistDefAttFlt(vlist,CDI_GLOBAL,TRIM(att_name),DATATYPE_FLT32,1,att_flt)
  END SUBROUTINE addGlobAttFlt

  SUBROUTINE addVar(var,k_jg)
    INTEGER, INTENT(IN)    :: var,k_jg

    num_varids(k_jg)              = num_varids(k_jg) + 1
    varids(num_varids(k_jg),k_jg) = var
  END SUBROUTINE addVar

  FUNCTION TimeVar(vname,vlongname,vunit,vcode,vtable,vlist,grid,zaxis) RESULT(var)
    INTEGER                  :: var
    CHARACTER(*), INTENT(IN) :: vname, vlongname, vunit
    INTEGER, INTENT(IN)      :: vcode, vtable, vlist, grid, zaxis

    var = vlistdefvar(vlist, grid, zaxis, TIME_VARIABLE)
    CALL vlistdefvarname(vlist, var, vname)
    CALL vlistdefvarlongname (vlist, var, vlongname)
    CALL vlistdefvarunits(vlist, var, vunit)
    IF ( vcode .gt. 0 ) THEN
      CALL vlistdefvarcode(vlist, var, vcode)
    ELSE
      CALL message('WARNING:TimeVar','Prevent setting negative var code for'//TRIM(vname))
    END IF
  END FUNCTION TimeVar


  SUBROUTINE nf(status)
    INTEGER, INTENT(in) :: status

    IF (status /= nf_noerr) THEN
      CALL finish('mo_io_vlist netCDF error', nf_strerror(status))
    ENDIF

  END SUBROUTINE nf



END MODULE mo_prepicon_utils

