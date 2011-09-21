!>
!! Does I/O based on cdilib.
!!
!!
!! @par Revision History
!! Initial implementation by Luis Kornblueh (2008-04-28)
!! Outputs for NCAR workshop by A. Gassmann and J. Foerstner (2008-05-27)
!!  Modified by Marco Giorgetta, MPI-M (2009-02-26)
!!  - renamed ltracer to ltransport
!!  Modified by Marco Giorgetta, MPI-M (2009-03-18)
!!  - write variable attributes longname, units, code and table to
!!    netcdf file, for postprocessing
!!  Modified by Marco Giorgetta, MPI-M (2009-03-22)
!!  - write (selected) namelist parameters as global attributes to
!!    netcdf file, for postprocessing
!!  Modified by Almut Gassmann, MPI_M (2009-04-15)
!!  - handling of output for the nonhydrostatic model
!!  Modified by Constantin Junk, MPI-M (2010-12-06)
!!  - write namelist parameters (/echam_phy_nml/,/echam_conv_ctl/
!!    /nwp_phy_ctl/) as global attributes to netcdf file, for postprocessing
!!  Modified by Rainer Johanni (2010-12-06)
!!  - complete restructuring of the output routines to list based form
!!  Modified by Constantin Junk, MPI-M (2011-01-17)
!!  - write namelist parameters /radiation_nml/,/transport_ctl, /nonhydrostatic_ctl/
!!    /testcase_ctl/ and /nh_testcase_ctl/
!!    as global attributes to netcdf file, for postprocessing
!!  Modified by F. Prill, DWD (2011-08)
!!  - optional output of variables interpolated onto lon-lat grid.
!!
!!
!! @par Copyright
!! 2002-2008 by DWD and MPI-M
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
!!--------------------------------------------------------------------------
!!
!! Adding variables to be output:
!!
!! - add the variable definition (with name, long name etc) in setup_vlist
!! - add an entry in get_outvar_ptr_ha and/or get_outvar_ptr_nh delivering
!!   a pointer to the actual variable to be output.
!!
!! get_outvar_ptr_ha and/or get_outvar_ptr_nh return 2 additional flags:
!!
!! - reset:  If this flag is set, the variable will be reset after I/O
!!           (for some collective precipitation fields)
!!
!! - delete: If this flag is set, the variable will be deallocated after I/O,
!!           because it is not pointing to another variable but has been
!!           allocated in the routine.
!!           This is usefull to create some derived output variables
!!           on the fly in get_outvar_ptr_ha/get_outvar_ptr_nh without
!!           the need of an extra array in the diagnostics.
!!
!! This mechanism should it make easy to add/change variables without
!! having to do the same changes at several locations.
!!
!! Please note that with this mechanism setup_vlist and destruct_vlist
!! must be called from EVERY PE in order to have a common view of the
!! output variables whereas open_output_vlist and close_output_vlist
!! must be only called from the PE actually doing the output!
!!
!!--------------------------------------------------------------------------

MODULE mo_io_vlist

  !-------------------------------------------------------------------------
  !
  !    ProTeX FORTRAN source: Style 2
  !    modified for ICON project, DWD/MPI-M 2006
  !
  !-------------------------------------------------------------------------
  !
  !
  !
  USE mo_kind,                  ONLY: wp
  USE mo_exception,             ONLY: finish, message, message_text
  USE mo_datetime,              ONLY: t_datetime, print_datetime
  USE mo_impl_constants,        ONLY: max_char_length, max_dom, modelname,        &
    &                                 modelversion, icc, zml_soil,                &
    &                                 max_ntracer,                                &
    &                                 ntrac_oce, ihs_atm_temp, ihs_atm_theta,     &
    &                                 inh_atmosphere, ishallow_water,             &
    &                                 inwp, iecham,ildf_echam, ihs_ocean
  USE mo_nonhydrostatic_config, ONLY: rayleigh_coeff, damp_height, iadv_rhotheta, &
    &                                 vwind_offctr, igradp_method, exner_expol,   &
    &                                 ltheta_up_hori, ltheta_up_vert,             &
    &                                 gmres_rtol_nh, iadv_rcf, ivctype,           &
    &                                 upstr_beta, l_open_ubc, l_nest_rcf,         &
    &                                 itime_scheme_nh_atm => itime_scheme
  USE mo_ocean_nml,             ONLY: n_zlev, dzlev_m, iforc_oce,no_tracer
  USE mo_dynamics_config,       ONLY: iequations,lshallow_water,                  &
    &                                 idiv_method, divavg_cntrwgt,                &
    &                                 nold, nnow, lcoriolis
  USE mo_ha_dyn_config,         ONLY: ha_dyn_config
  USE mo_diffusion_config,      ONLY: diffusion_config
  USE mo_io_config,             ONLY: lwrite_omega, lwrite_pres, lwrite_z3,       &
    &                                 lwrite_vorticity, lwrite_divergence,        &
    &                                 lwrite_tend_phy, lwrite_radiation,          &
    &                                 lwrite_precip, lwrite_cloud, lwrite_tracer, &
    &                                 lwrite_tke,  lwrite_surface, lwrite_pzlev,  &
    &                                 lwrite_extra, inextra_2d,inextra_3d,        &
    &                                 out_filetype, out_expname,                  &
    &                                 dt_data, dt_file, lkeep_in_sync
  USE mo_nh_pzlev_config,       ONLY: nh_pzlev_config
  USE mo_parallel_config,       ONLY: nproma, p_test_run
  USE mo_extpar_config,         ONLY: itopo
  USE mo_run_config,            ONLY: num_lev, num_levp1, iforcing, lforcing,     &
    &                                 ntracer, ltransport, nsteps, dtime,         &
    &                                 ldynamics, ltestcase, lvert_nest, msg_level,&
    &                                 iqv, iqc, iqi, iqcond, ntracer_static, io3
  USE mo_grid_config,           ONLY: global_cell_type
  USE mo_echam_phy_config
  USE mo_atm_phy_nwp_config,    ONLY: atm_phy_nwp_config
  USE mo_advection_config,      ONLY: advection_config
  USE mo_echam_conv_config,     ONLY: echam_conv_config
  USE mo_lnd_nwp_config,        ONLY: nsfc_subs, nlev_snow
! USE mo_gw_hines_nml,          ONLY: lheatcal, emiss_lev, rmscon, kstar, m_min
  USE mo_vertical_coord_table,  ONLY: vct
  USE mo_model_domain_import,   ONLY: start_lev, nroot, n_dom, lfeedback, lplane, &
    & n_dom_start
  USE mo_model_domain,          ONLY: t_patch, p_patch
  USE mo_physical_constants,    ONLY: grav
  USE mo_communication,         ONLY: exchange_data, t_comm_pattern
  USE mo_mpi,                   ONLY: my_process_is_mpi_workroot, my_process_is_stdio, &
    &  my_process_is_mpi_test, my_process_is_mpi_seq, process_mpi_all_test_id,         &
    &  process_mpi_all_workroot_id, p_recv, p_send, num_work_procs
  USE mo_icoham_dyn_types,      ONLY: t_hydro_atm_prog, t_hydro_atm_diag
  USE mo_nonhydro_state,        ONLY: t_nh_prog, t_nh_diag, t_nh_diag_pz
  USE mo_oce_state,             ONLY: t_hydro_ocean_state, t_hydro_ocean_prog,       &
       &                              t_hydro_ocean_diag, t_hydro_ocean_base,        &
       &                              t_hydro_ocean_aux,                             &
       &                              v_base, set_zlev, v_ocean_state
  USE mo_oce_forcing,           ONLY: t_sfc_flx, v_sfc_flx
  USE mo_ext_data,              ONLY: t_external_ocean
  USE mo_icoham_dyn_memory,     ONLY: p_hydro_state
  USE mo_nonhydro_state,        ONLY: p_nh_state
  USE mo_nwp_lnd_state,         ONLY: p_lnd_state
  USE mo_nwp_phy_state,         ONLY: prm_diag, prm_nwp_tend !, t_nwp_phy_diag
  USE mo_nwp_lnd_state,         ONLY: t_lnd_prog, t_lnd_diag
  USE mo_echam_phy_memory,      ONLY: prm_field, prm_tend
  USE mo_icoham_sfc_indices,    ONLY: nsfc_type
  USE mo_radiation_config,      ONLY: izenith, irad_h2o,                          &
    &                                 irad_co2, irad_ch4, irad_n2o, irad_o3,      &
    &                                 irad_o2, irad_cfc11, irad_cfc12,  irad_aero
  USE mo_ha_testcases,          ONLY: ctest_name, ihs_init_type, lhs_vn_ptb,      &
    &                                 hs_vn_ptb_scale, lrh_linear_pres,           &
    &                                 rh_at_1000hpa,linit_tracer_fv
  USE mo_gw_test,               ONLY: gw_brunt_vais, gw_u0,gw_lon_deg, gw_lat_deg
  USE mo_rh_test,               ONLY: rh_wavenum, rh_init_shift_deg
  USE mo_jw_test,               ONLY: jw_uptb
  USE mo_mrw_test,              ONLY: mountctr_lon_deg, mountctr_lat_deg,         &
    &                                 mountctr_height, mount_half_width,          &
    &                                 mount_u0
  USE mo_nh_testcases,          ONLY: nh_test_name, mount_height,                 &
    &                                 torus_domain_length, nh_brunt_vais, nh_u0,  &
    &                                 nh_t0, jw_up,                               &
    &                                 rotate_axis_deg, lhs_nh_vn_ptb,             &
    &                                 hs_nh_vn_ptb_scale, qv_max, ape_sst_case
    !&                                rh_at_1000hpa,linit_tracer_fv
  USE mo_nh_mrw_exp,            ONLY: u0_mrw, mount_height_mrw,                   &
    &                                 mount_lonctr_mrw_deg, mount_latctr_mrw_deg, &
    &                                 p_int_mwbr_const, temp_i_mwbr_const,        &
    &                                 bruntvais_u_mwbr_const     
    !&                                mount_half_width
  USE mo_lonlat_intp_config,    ONLY: lonlat_intp_config
  USE mo_math_constants,        ONLY: pi
  USE mo_impl_constants,        ONLY: SUCCESS
  USE mo_math_utilities,        ONLY: t_lon_lat_grid
  USE mo_interpolation,         ONLY: t_lon_lat_intp,                               &
    &                                 rbf_interpol_lonlat, rbf_vec_interpol_lonlat, &
    &                                 verts2cells_scalar
  USE mo_intp_data_strc,        ONLY: p_int_state, p_int_state_lonlat
  USE mo_mpi,                   ONLY: p_pe
  USE mo_util_string,           ONLY: string_contains_word, toupper

  IMPLICIT NONE

  PRIVATE

  INCLUDE 'cdi.inc'
  INCLUDE 'netcdf.inc'

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

  INTEGER, PARAMETER :: max_outvars  = 250 ! max. number of output variables
  INTEGER, PARAMETER :: max_gridlevs = 12 ! max. number of grid levels
  INTEGER, PARAMETER :: max_name_len = 100

  !> level of output verbosity
  INTEGER, PARAMETER  :: dbg_level = 0

  PUBLIC :: setup_vlist, destruct_vlist,                            &
    &       open_output_vlist, close_output_vlist,                  &
    &       write_vlist, get_outvar_ptr_ha, get_outvar_ptr_nh,      &
    &       vlist_write_var, vlist_set_date_time, vlist_start_step, &
    &       de_reshape1, de_reshape2,                               &
    &       gather_array1, gather_array2, t_collected_var_ptr

  PRIVATE :: addGlobAttInt, addGlobAttTxt, addGlobAttFlt
  ! I/O stream handler
  INTEGER, DIMENSION(max_gridlevs) ::  &
    &  streamID

  INTEGER, DIMENSION(max_gridlevs) ::  &
    &  gridCellID, gridEdgeID, gridVertexID

  INTEGER, DIMENSION(max_gridlevs) ::  &
    &  vlistID, taxisID, zaxisID_surface, zaxisID_hybrid,      &
    &  zaxisID_hybrid_half, zaxisIDdepth_m, zaxisID_halfdepth, &
    &  zaxisID_depth_below_land, zaxisID_depth_below_land_p1,  &
    &  zaxisID_generic_snow, zaxisID_generic_snow_p1,          &
    &  zaxisID_pres, zaxisID_hgt

  ! NetCDF identifier to lon-lat grid
  INTEGER  :: lonlat_gridID(max_gridlevs)

  ! Flag value, indicating that no lon-lat interpolation shall be
  ! performed
  INTEGER, PARAMETER :: NO_LONLAT_INTERP = -1

  ! list of variables to be interpolated onto lon-lat grid
  ! lonlat_varids(i,j,1:2) contains var ID of lonlat variable
  ! corresponding to standard variable with var ID varids(i,j)
  !
  ! @note For scalar variables only lonlat_varids(:,:,1) is relevant.
  INTEGER  :: lonlat_varids(max_outvars,max_gridlevs,2)

  INTEGER, DIMENSION(max_outvars,max_gridlevs) ::  &
    &  varids
  ! current number of output variables, gets updated by addVar()
  INTEGER, PRIVATE :: num_varids(max_dom)

  INTEGER, SAVE :: iostep = 0

  INTEGER, PARAMETER, PUBLIC :: GATHER_C      = 1
  INTEGER, PARAMETER, PUBLIC :: GATHER_E      = 2
  INTEGER, PARAMETER, PUBLIC :: GATHER_V      = 3
  INTEGER, PARAMETER, PUBLIC :: GATHER_LONLAT = 4

  ! Descriptions of output variables

  INTEGER, PUBLIC :: num_output_vars(max_gridlevs)

  TYPE t_outvar_desc

    INTEGER ::           TYPE ! GATHER_C, GATHER_E, GATHER_V, GATHER_LONLAT
    INTEGER ::           nlev
    CHARACTER(LEN=80) :: name

  END TYPE t_outvar_desc

  !> pointer type (avoids clash of INTENT and POINTER attributes)
  TYPE t_collected_var_ptr
    REAL(wp), POINTER :: ptr(:,:,:)
  END TYPE t_collected_var_ptr

  PUBLIC :: t_outvar_desc
  INTEGER :: klev

  TYPE(t_outvar_desc), PUBLIC :: outvar_desc(max_outvars, max_gridlevs)

  ! data fields for (optional) interpolation of output variables to
  ! lon-lat grid (only allocated on IO PEs, only one field for all
  ! variables)
  TYPE t_lonlat_data
    !> interpolated variables
    REAL(wp), POINTER :: &
      &  var_3d(:,:,:), var_3d_2(:,:,:)  ! (nproma,nlev,nblks_lonlat)
  END TYPE t_lonlat_data

  TYPE(t_lonlat_data), TARGET :: var_lonlat(max_gridlevs) ! (1, ..., no. of patches)
  
CONTAINS

  !-------------------------------------------------------------------------
  !BOC
  SUBROUTINE setup_vlist(grid_filename, k_jg)

    CHARACTER(len=*), INTENT(in) :: grid_filename
    INTEGER, INTENT(in) :: k_jg

    INTEGER :: ncid, dimid, varid
    INTEGER :: i_nc, i_ne, i_nv
    INTEGER :: i_ncb, i_neb, i_nvb
    INTEGER :: lnlen, ulen, nzlevp1

    INTEGER :: nlev, nlevp1, nplev, nzlev
    INTEGER :: znlev_soil

    REAL(wp), ALLOCATABLE :: clon(:), clat(:), clonv(:), clatv(:)
    REAL(wp), ALLOCATABLE :: elon(:), elat(:), elonv(:), elatv(:)
    REAL(wp), ALLOCATABLE :: vlon(:), vlat(:), vlonv(:), vlatv(:)

    REAL(wp), ALLOCATABLE :: levels(:), levels_i(:), levels_m(:)

    CHARACTER(len=21) :: name
    CHARACTER(len=12) :: qname
    CHARACTER(len=10) :: dbgname
    CHARACTER(len=3)  :: cjt
    CHARACTER(LEN=1)  :: ctracer
    CHARACTER(len=MAX_CHAR_LENGTH) :: & !< list of tracers to initialize
      &  ctracer_list
    CHARACTER(LEN=1)  :: anextra ! number of debug fields
    CHARACTER(len=NF_MAX_NAME) :: long_name, units
    INTEGER :: i, k, jt, itracer
    INTEGER :: ivar
    INTEGER :: gridid, zaxisid
    INTEGER :: elemid, tableid

    !CHARACTER(len=NF_MAX_NAME) :: att_txt
    INTEGER                    :: astatus

    INTEGER                    :: &
      & errstat
    REAL(wp), ALLOCATABLE ::     &
      &  p_lonlat(:)        ! grid longitude pts (interpolation)
    REAL(wp) :: pi_180
    INTEGER  :: ist         ! error code

    ! ocean tracers
    CHARACTER(len=max_char_length) :: oce_tracer_names(max_ntracer),&
    &                                 oce_tracer_units(max_ntracer),&
    &                                 oce_tracer_longnames(max_ntracer)
    INTEGER                    :: oce_tracer_codes(max_ntracer)
    INTEGER                    :: oce_trace_counter
    INTEGER, PARAMETER         :: oce_max_tracer = 2

    CHARACTER(len=max_char_length) :: msg
    CHARACTER(len=max_char_length), PARAMETER :: routine = 'mo_io_vlist:setup_vlist'
    INTEGER                    :: nblks
    TYPE (t_lon_lat_grid), POINTER :: grid ! lon-lat grid

    !-------------------------------------------------------------------------

    pi_180 = ATAN(1._wp)/45._wp

    !------------------------------------------------------------------
    ! no grid refinement allowed here so far the ocean
    !------------------------------------------------------------------
    IF (iequations==ihs_ocean) THEN 
      IF (k_jg > 1 ) THEN !b
        CALL finish(TRIM(routine), ' k_jg > 1 is not allowed')
      END IF
    END IF
    ! Each time a new NetCDF is created, reset "iostep" to zero
    ! (Otherwise we will get an error message from a CDI subroutine.)

    iostep = 0

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
    
    ALLOCATE(clon(i_nc))
    
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
    DEALLOCATE(clon)


    SELECT CASE (global_cell_type)
    CASE (3)
      name = 'clat'
    CASE (6)
      name = 'vlat'
    END SELECT
    
    ALLOCATE(clat(i_nc))
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
    DEALLOCATE(clat)
    
    !
    SELECT CASE (global_cell_type)
    CASE (3)
      CALL nf(nf_inq_varid(ncid, 'clon_vertices', varid))
    CASE (6)
      CALL nf(nf_inq_varid(ncid, 'vlon_vertices', varid))
    END SELECT
    
    ALLOCATE(clonv(i_ncb))
    CALL nf(nf_get_var_double(ncid, varid, clonv))
    !
    CALL gridDefXbounds(gridCellID(k_jg), clonv)
    DEALLOCATE(clonv)

    !
    SELECT CASE (global_cell_type)
    CASE (3)
      CALL nf(nf_inq_varid(ncid, 'clat_vertices', varid))
    CASE (6)
      CALL nf(nf_inq_varid(ncid, 'vlat_vertices', varid))
    END SELECT

    ALLOCATE(clatv(i_ncb))
    CALL nf(nf_get_var_double(ncid, varid, clatv))
    !
    CALL gridDefYbounds(gridCellID(k_jg), clatv)
    DEALLOCATE(clatv)
    !
    !-------------------------------------------------------------------------
    ! edge grid
    !
    gridEdgeID(k_jg) = gridCreate(GRID_UNSTRUCTURED, i_ne)
    CALL gridDefNvertex(gridEdgeID(k_jg), 4)
    !
    name = 'elon'
    ALLOCATE(elon(i_ne))
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
    DEALLOCATE(elon)
    !
    name = 'elat'
    ALLOCATE(elat(i_ne))
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
    DEALLOCATE(elat)
    !
    ALLOCATE(elonv(i_neb))
    CALL nf(nf_inq_varid(ncid, 'elon_vertices', varid))
    CALL nf(nf_get_var_double(ncid, varid, elonv))
    !
    CALL gridDefXbounds(gridEdgeID(k_jg), elonv)
    DEALLOCATE(elonv)
    !
    ALLOCATE(elatv(i_neb))
    CALL nf(nf_inq_varid(ncid, 'elat_vertices', varid))
    CALL nf(nf_get_var_double(ncid, varid, elatv))
    !
    CALL gridDefYbounds(gridEdgeID(k_jg), elatv)
    DEALLOCATE(elatv)
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
    ALLOCATE(vlon(i_nv))
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
    DEALLOCATE(vlon)
    !
    SELECT CASE (global_cell_type)
    CASE (3)
      name = 'vlat'
    CASE (6)
      name = 'clat'
    END SELECT
    ALLOCATE(vlat(i_nv))
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
    DEALLOCATE(vlat)
    !
    IF(global_cell_type==3) THEN
      CALL nf(nf_inq_varid(ncid, 'vlon_vertices', varid))
    ELSEIF(global_cell_type==6) THEN
      CALL nf(nf_inq_varid(ncid, 'clon_vertices', varid))
    ENDIF
    ALLOCATE(vlonv(i_nvb))
    CALL nf(nf_get_var_double(ncid, varid, vlonv))
    !
    CALL gridDefXbounds(gridVertexID(k_jg), vlonv)
    DEALLOCATE(vlonv)
    !
    IF(global_cell_type==3) THEN
      CALL nf(nf_inq_varid(ncid, 'vlat_vertices', varid))
    ELSEIF(global_cell_type==6) THEN
      CALL nf(nf_inq_varid(ncid, 'clat_vertices', varid))
    ENDIF
    ALLOCATE(vlatv(i_nvb))
    CALL nf(nf_get_var_double(ncid, varid, vlatv))
    !
    CALL gridDefYbounds(gridVertexID(k_jg), vlatv)
    DEALLOCATE(vlatv)
    !
    !=========================================================================
    ! vertical grids
    !
    ! surface level
    zaxisID_surface(k_jg) = zaxisCreate(ZAXIS_SURFACE, 1)
    ALLOCATE(levels(1))
    levels(1) = 0.0_wp
    CALL zaxisDefLevels(zaxisID_surface(k_jg), levels)
    DEALLOCATE(levels)
    ! atm (pressure) height, ocean depth
    IF (iequations/=ihs_ocean) THEN ! atm 

      nlev   = num_lev(k_jg)
      nlevp1 = num_levp1(k_jg)
      ! introduce temporary variable znlev_soil, since global variable nlev_soil 
      ! is unknown to the I/O-Processor. Otherwise receive_patch_configuration in 
      ! mo_io_async complains about mismatch of levels. 
      znlev_soil = SIZE(zml_soil)-1

      zaxisID_hybrid(k_jg)      = zaxisCreate(ZAXIS_HYBRID, nlev)
      zaxisID_hybrid_half(k_jg) = zaxisCreate(ZAXIS_HYBRID_HALF, nlevp1)

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
      !
      zaxisID_depth_below_land_p1(k_jg) = zaxisCreate(ZAXIS_DEPTH_BELOW_LAND, znlev_soil+2)
      ALLOCATE(levels(znlev_soil+2))
      levels(1) = 0._wp
      DO i = 1, znlev_soil+1
      levels(i+1) = zml_soil(i)*100._wp
      END DO
      CALL zaxisDefLevels(zaxisID_depth_below_land_p1(k_jg), levels)
      DEALLOCATE(levels)


      ! Define axes for soil model
      !
      zaxisID_depth_below_land(k_jg) = zaxisCreate(ZAXIS_DEPTH_BELOW_LAND, znlev_soil+1)
      CALL zaxisDefLevels(zaxisID_depth_below_land(k_jg), zml_soil*100._wp)
      !
      zaxisID_generic_snow_p1(k_jg) = zaxisCreate(ZAXIS_GENERIC, nlev_snow+1)
      ALLOCATE(levels(nlev_snow+1))
      DO i = 1, nlev_snow+1
        levels(i) = REAL(i,wp)
      END DO
      CALL zaxisDefLevels(zaxisID_generic_snow_p1(k_jg), levels)
      DEALLOCATE(levels)
      !
      zaxisID_generic_snow(k_jg) = zaxisCreate(ZAXIS_GENERIC, nlev_snow)
      ALLOCATE(levels(nlev_snow))
      DO i = 1, nlev_snow
        levels(i) = REAL(i,wp)
      END DO
      CALL zaxisDefLevels(zaxisID_generic_snow(k_jg), levels)
      DEALLOCATE(levels)


      ! Define axes for output on p- and z-levels
      !
      IF (lwrite_pzlev) THEN
        nplev = nh_pzlev_config(k_jg)%nplev
        zaxisID_pres(k_jg) = zaxisCreate(ZAXIS_PRESSURE, nplev)
        ALLOCATE(levels(nplev))
        DO i = 1, nplev
          levels(i) = nh_pzlev_config(k_jg)%plevels(i)
        END DO
        CALL zaxisDefLevels(zaxisID_pres(k_jg), levels)
        CALL zaxisDefVct(zaxisID_pres(k_jg), nplev, levels)
        DEALLOCATE(levels)

        nzlev = nh_pzlev_config(k_jg)%nzlev
        zaxisID_hgt(k_jg)  = zaxisCreate(ZAXIS_HEIGHT, nzlev)
        ALLOCATE(levels(nzlev))
        DO i = 1, nzlev
          levels(i) = nh_pzlev_config(k_jg)%zlevels(i)
        END DO
        CALL zaxisDefLevels(zaxisID_hgt(k_jg), levels)
        CALL zaxisDefVct(zaxisID_hgt(k_jg), nzlev, levels)
        DEALLOCATE(levels)
      ENDIF

    ELSE ! oce
      zaxisIDdepth_m(k_jg)  = zaxisCreate(ZAXIS_DEPTH_BELOW_SEA, n_zlev)
      nzlevp1 = n_zlev + 1
      zaxisID_halfdepth(k_jg)  = zaxisCreate(ZAXIS_DEPTH_BELOW_SEA, nzlevp1)

      ALLOCATE(levels_i(nzlevp1))
      ALLOCATE(levels_m(n_zlev))
      CALL set_zlev(n_zlev, dzlev_m, levels_i, levels_m)
      CALL zaxisDefLevels(zaxisIDdepth_m(k_jg)   , levels_m)
      CALL zaxisDefLevels(zaxisID_halfdepth(k_jg), levels_i)
      DEALLOCATE(levels_i)
      DEALLOCATE(levels_m)
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
    CALL addGlobAtts(vlistID(k_jg),k_jg,astatus)
    IF (iequations/=ihs_ocean) THEN
      CALL addAtmAtts(vlistID(k_jg),k_jg,astatus)
    ELSE
      CALL addOceAtts(vlistID(k_jg),astatus)
    END IF

    !-------------------------------------------------------------------------
    ! create lon-lat grid data structure:

    LONLAT : IF (lonlat_intp_config(k_jg)%l_enabled) THEN
      ! id list for lon-lat interpolation
      lonlat_varids(:,k_jg,:) = NO_LONLAT_INTERP
      grid => lonlat_intp_config(k_jg)%lonlat_grid

      IF (dbg_level > 1) &
        CALL message(routine, "Allocating field for interpolated variables")
      nblks = grid%nblks
      ALLOCATE(var_lonlat(k_jg)%var_3d  (nproma,MAX(nlev,nlevp1),nblks), &
        &      var_lonlat(k_jg)%var_3d_2(nproma,MAX(nlev,nlevp1),nblks), &
        &      STAT=ist)
      IF (ist /= SUCCESS) THEN
        CALL finish (routine, 'field allocation failed')
      ENDIF

      ! define lon-lat grid
      IF (dbg_level > 1) THEN
        WRITE(message_text,*) "PE #", p_pe, "Definition of interpolation lon-lat grid (", k_jg, ")"
        CALL message(routine, message_text)
      END IF

      lonlat_gridID(k_jg) = gridCreate(GRID_LONLAT, grid%total_dim)

      CALL gridDefXsize(lonlat_gridID(k_jg), grid%dimen(1))
      CALL gridDefXname(lonlat_gridID(k_jg), 'interp_lon')
      CALL gridDefXunits(lonlat_gridID(k_jg), units(1:ulen))

      CALL gridDefYsize(lonlat_gridID(k_jg), grid%dimen(2))
      CALL gridDefYname(lonlat_gridID(k_jg), 'interp_lat')
      CALL gridDefYunits(lonlat_gridID(k_jg), units(1:ulen))

      ! define rotation info of lon-lat grid:
      CALL gridDefXpole(lonlat_gridID(k_jg), grid%poleN(1))
      CALL gridDefYpole(lonlat_gridID(k_jg), grid%poleN(2))
      ! Note: CALL gridDefAngle() not yet supported

      DO i=1,2 ! loop lon/lat
        ALLOCATE(p_lonlat(grid%dimen(i)), stat=errstat)
        IF (errstat /= SUCCESS)  CALL finish('setup_vlist','ALLOCATE failed.')      

        ! generate list of lon/lat points
        p_lonlat(1:grid%dimen(i)) =              &
          &      (/ (grid%sw_corner(i)           &
          &            +  (k-1)*grid%delta(i),   &
          &      k=1,grid%dimen(i))  /) 

        IF (i==1) THEN
          CALL gridDefXvals(lonlat_gridID(k_jg), p_lonlat)
        ELSE
          CALL gridDefYvals(lonlat_gridID(k_jg), p_lonlat)
        END IF

        DEALLOCATE(p_lonlat, stat=errstat)
        IF (errstat /= SUCCESS)  CALL finish('setup_vlist','DEALLOCATE failed.')      
      END DO
    END IF LONLAT

    !-------------------------------------------------------------------------
    ! register variables
    varids(:,k_jg)   = 0
    ! atm
    IF (iequations/=ihs_ocean) THEN
      ! initialize total number of varids for domain jg
      num_varids(k_jg) = 0

      ! get ctracer_list
      ctracer_list = advection_config(k_jg)%ctracer_list

      SELECT CASE (iforcing)
      CASE (iecham,ildf_echam)
        ! land-sea mask (1. = land, 0. = sea/lakes)
        CALL addVar(TimeVar('SLM',&
        &                   'land-sea mask',&
        &                   '', 172, 128,&
        &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
        &          k_jg)
        ! grid-box mean skin temperature
        CALL addVar(TimeVar('SKT',&
        &                   'skin temperature',&
        &                   'K', 235, 128,&
        &                    vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
        &           k_jg)
      CASE DEFAULT
      END SELECT
      ! surface pressure
      CALL addVar(TimeVar('PS',&
      &                   'surface pressure',&
      &                   'Pa', 134, 128,&
      &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
      &           k_jg)
      ! temperature
      IF (.NOT.lshallow_water) THEN
        CALL addVar(TimeVar('T', &
        &                   'temperature',&
        &                   'K', 130, 128,&
        &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_hybrid(k_jg)),&
        &           k_jg)
      ENDIF
      ! normal velocity at the edges
      CALL addVar(TimeVar('normal_velocity',&
      &                   'velocity normal to edge',&
      &                   'm/s', 254, 128, &
      &                   vlistID(k_jg), gridEdgeID(k_jg),zaxisID_hybrid(k_jg)),&
      &           k_jg)
      ! zonal wind
      CALL addVar(TimeVar('U',&
      &                   'zonal wind',&
      &                   'm/s', 131, 128,&
      &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_hybrid(k_jg)),&
      &           k_jg)
      ! meridional wind
      CALL addVar(TimeVar('V',&
      &                   'meridional wind',&
      &                   'm/s', 132, 128, &
      &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_hybrid(k_jg)),&
      &           k_jg)
      ! vertical velocity ( dp/dt )
      IF (iequations /= 3 .AND. lwrite_omega) THEN
        CALL addVar(TimeVar('OMEGA',&
        &                   'vertical velocity',&
        &                   'Pa/s', 135, 128,&
        &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_hybrid(k_jg)),&
        &          k_jg)
      END IF
      ! vertical velocity ( w ) in nonhydrostatic model (code number definition?)
      IF (iequations == 3) THEN
        CALL addVar(TimeVar('W',&
        &                   'upward air velocity',&
        &                   'm/s', 40, 2,&
        &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_hybrid_half(k_jg)),&
        &          k_jg)
      END IF
      ! pressure
      IF (lwrite_pres) THEN
        CALL addVar(TimeVar('P',&
        &                   'pressure',&
        &                   'Pa', 255, 128,&
        &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_hybrid(k_jg)),&
        &          k_jg)
      END IF
      ! geopotential height
      IF (lwrite_z3) THEN
        CALL addVar(TimeVar('ZF3',&
        &                   'geopotential height',&
        &                   'm', 156, 128,&
        &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_hybrid(k_jg)),&
        &          k_jg)
        IF (iequations == 3) THEN
          CALL addVar(TimeVar('ZH3',&
          &                   'half level height',&
          &                   'm', 253, 128,&
          &                    vlistID(k_jg), gridCellID(k_jg), zaxisID_hybrid_half(k_jg)),&
          &           k_jg)
        ENDIF
      END IF
      ! surface geopotential
      IF(iequations == 1 .OR. &
         iequations == 2 ) THEN
        CALL addVar(TimeVar('PHIS',&
        &                   'surface geopotential (orography)',&
        &                   'm**2/s**2',&
        &                   129, 128,&
        &                   vlistID(k_jg),gridCellID(k_jg),zaxisID_surface(k_jg)),&
        &           k_jg)
      ENDIF
      ! tracer fields
      IF (ntracer > 0 .AND. (iforcing == inwp .OR. ltransport)) THEN
        DO jt = 1, ntracer
          IF (lwrite_tracer(jt)) THEN
            ctracer = ctracer_list(jt:jt)
            IF (jt.EQ.1) THEN !for water vapour
             elemid=51; tableid=2   !DWD coding
            ELSEIF (jt.EQ.2) THEN !for cloud water
             elemid=31; tableid=201 !DWD coding
            ELSEIF (jt.EQ.3) THEN !for cloud ice
             elemid=33; tableid=201 !DWD coding
            ELSE !other tracers
             elemid=252; tableid=128 !default coding
            END IF
            WRITE(name,'(A1,A1)') "Q", ctracer
            CALL addVar(TimeVar(TRIM(name),TRIM(name),&
            &                   'kg/kg',elemid,tableid,&
            &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_hybrid(k_jg)),&
            &           k_jg)       
          END IF
        END DO
        IF (iforcing == inwp)THEN
         DO jt = 1, 3
          IF (lwrite_tracer(jt)) THEN
            ctracer = ctracer_list(jt:jt)
            WRITE(name,'(A2,A1)') "TQ", ctracer
            WRITE(long_name,'(A34,A1)') "vertically integrated grid-scale Q",ctracer
            CALL addVar(TimeVar(TRIM(name),TRIM(long_name),&
            &                   'kg/m**2',222,128,&
            &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
            &           k_jg)
            WRITE(name,'(A2,A1,A4)') "TQ", ctracer,"_avg"
            WRITE(long_name,'(A42,A1)') "average vertically integrated grid-scale Q",ctracer
            CALL addVar(TimeVar(TRIM(name),TRIM(long_name),&
            &                   'kg/m**2',222,128,&
            &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
            &           k_jg)
            
          END IF !lwrite_tracer(jt)
         END DO
         IF (irad_o3 == 4 .OR. irad_o3 ==6 ) THEN                  !output for O3
                        CALL addVar(TimeVar('O3','O3',&
            &                   'kg/kg',999,999,&
            &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_hybrid(k_jg)),&
            &           k_jg)  
        
         END IF
        ENDIF !iforcing == inwp



      END IF

    ! debug fields:
    ! 2d - level type "surface"
  !!$  CALL addVar(DebugVar('debug_2d_1',vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),k_jg)
  !!$  CALL addVar(DebugVar('debug_2d_2',vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),k_jg)
  !!$  CALL addVar(DebugVar('debug_2d_3',vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),k_jg)
  !!$  CALL addVar(DebugVar('debug_2d_4',vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),k_jg)
  !!$  CALL addVar(DebugVar('debug_2d_5',vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),k_jg)
  !!$  CALL addVar(DebugVar('debug_2d_6',vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),k_jg)
  !!$  CALL addVar(DebugVar('debug_2d_7',vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),k_jg)
  !!$  CALL addVar(DebugVar('debug_2d_8',vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),k_jg)
    !
    ! 3d - level type "hybrid"
  !!$  CALL addVar(DebugVar('debug_3d_1',vlistID(k_jg), gridCellID(k_jg),zaxisID_hybrid(k_jg)), k_jg)
  !!$  CALL addVar(DebugVar('debug_3d_2',vlistID(k_jg), gridCellID(k_jg),zaxisID_hybrid(k_jg)), k_jg)
  !!$  CALL addVar(DebugVar('debug_3d_3',vlistID(k_jg), gridCellID(k_jg),zaxisID_hybrid(k_jg)), k_jg)
  !!$  CALL addVar(DebugVar('debug_3d_4',vlistID(k_jg), gridCellID(k_jg),zaxisID_hybrid(k_jg)), k_jg)
  !!$  CALL addVar(DebugVar('debug_3d_5',vlistID(k_jg), gridCellID(k_jg),zaxisID_hybrid(k_jg)), k_jg)
  !!$  CALL addVar(DebugVar('debug_3d_6',vlistID(k_jg), gridCellID(k_jg),zaxisID_hybrid(k_jg)), k_jg)
  !!$  CALL addVar(DebugVar('debug_3d_7',vlistID(k_jg), gridCellID(k_jg),zaxisID_hybrid(k_jg)), k_jg)
  !!$  CALL addVar(DebugVar('debug_3d_8',vlistID(k_jg), gridCellID(k_jg),zaxisID_hybrid(k_jg)), k_jg)

      IF (lwrite_extra) THEN
        SELECT CASE (iforcing)
        CASE (iecham,ildf_echam,inwp)

        DO jt = 1, inextra_2d ! only upto 9 numbers ok!!
          WRITE(anextra,'(I1)') jt
            WRITE(dbgname,'(A9,A1)') "extra_2d_", anextra
          CALL addVar(DebugVar(TRIM(dbgname), vlistID(k_jg), gridCellID(k_jg), &
            &         zaxisID_surface(k_jg)),k_jg)
        ENDDO
        DO jt = 1, inextra_3d
          WRITE(anextra,'(I1)') jt
            WRITE(dbgname,'(A9,A1)') "extra_3d_", anextra
          CALL addVar(DebugVar(TRIM(dbgname), vlistID(k_jg), gridCellID(k_jg), &
            &         zaxisID_hybrid(k_jg)), k_jg)
        ENDDO

        END SELECT !iforcing
      ENDIF !lwrite_extra

      ! radiation
      IF(lwrite_radiation) THEN

        SELECT CASE (iforcing)
        CASE (iecham,ildf_echam,inwp)
          CALL addVar(TimeVar('cosmu0',&
          &                   'cosine of zenith angle',&
          &                   ' ', 1, 201,&
          &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
          &           k_jg)
          CALL addVar(TimeVar('flxdwswtoa',&
          &                   'downward shortwave flux at TOA',&
          &                   'W/m**2', 184, 201,&
          &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
          &           k_jg)
          IF(iforcing == iecham) THEN
            CALL addVar(TimeVar('ozone',&
            &                   'ozone mixing ratio',&
            &                   'g/g', 255, 255,&
            &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_hybrid(k_jg)),&
            &           k_jg)
          ENDIF
        END SELECT

        SELECT CASE (iforcing)
        CASE (inwp)
          CALL addVar(TimeVar('swflxsfc_avg',&
          &                   'averaged shortwave surface net flux',&
          &                   'W/m**2', 111, 2,&
          &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
          &           k_jg)
          CALL addVar(TimeVar('lwflxsfc_avg',&
          &                   'averaged longwave surface net flux',&
          &                   'W/m**2', 112, 2,&
          &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
          &           k_jg)
          CALL addVar(TimeVar('swflxtoa_avg',&
          &                   'averaged shortwave toa net flux',&
          &                   'W/m**2', 113, 2,&
          &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
          &           k_jg)
          CALL addVar(TimeVar('lwflxtoa_avg',&
          &                   'averaged longwave toa net flux',&
          &                   'W/m**2', 114, 2,&
          &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
          &           k_jg)
          CALL addVar(TimeVar('swflxsfc',&
          &                   'shortwave surface net flux',&
          &                   'W/m**2', 111, 201,&
          &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
          &           k_jg)
          CALL addVar(TimeVar('lwflxsfc',&
          &                   'longwave surface net flux',&
          &                   'W/m**2', 112, 201,&
          &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
          &           k_jg)
          CALL addVar(TimeVar('swflxtoa',&
          &                   'shortwave toa net flux',&
          &                   'W/m**2', 113, 201,&
          &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
          &           k_jg)
          CALL addVar(TimeVar('lwflxtoa',&
          &                   'longwave toa net flux',&
          &                   'W/m**2', 114, 201,&
          &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
          &           k_jg)
        END SELECT

      END IF

      ! surface precitation rates
      IF(lwrite_precip) THEN

        SELECT CASE (iforcing)
        CASE (inwp)
          CALL addVar(TimeVar('PRR_GSP',&
          &                   'grid-scale rain precipitation rate',&
          &                   'kg/s/m**2', 100, 201,&
          &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
          &           k_jg)
          CALL addVar(TimeVar('PRS_GSP',&
          &                   'grid-scale snow precipitation rate',&
          &                   'kg/s/m**2', 101, 201,&
          &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
          &           k_jg)
          CALL addVar(TimeVar('RAIN_GSP',&
          &                   'grid-scale accumulated surface rain',&
          &                   'kg/m**2', 102, 201,&
          &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
          &           k_jg)
          CALL addVar(TimeVar('SNOW_GSP',&
          &                   'grid-scale accumulated surface snow',&
          &                   'kg/m**2', 79, 2,&
          &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
          &           k_jg)
          CALL addVar(TimeVar('RAIN_CON',&
          &                   'convective accumulated surface rain',&
          &                   'kg/m**2', 113, 201,&
          &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
          &           k_jg)
          CALL addVar(TimeVar('SNOW_CON',&
          &                   'convective accumulated surface snow',&
          &                   'kg/m**2', 78, 2,&
          &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
          &           k_jg)
          CALL addVar(TimeVar('TOT_PREC',&
          &                   'grid-scale + convective accumulated surface total precipitation',&
          &                   'kg/m**2', 61, 2,&
          &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
          &           k_jg)
          CALL addVar(TimeVar('TOT_PREC_RATE_avg',&
          &                   'average grid-scale + convective surface total precipitation rate',&
          &                   'kg/m**2/s', 61, 2,&
          &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
          &           k_jg)
          CALL addVar(TimeVar('CON_PREC_RATE_avg',&
          &                   'average convective surface precipitation rate',&
          &                   'kg/m**2/s', 61, 2,&
          &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
          &           k_jg)
          CALL addVar(TimeVar('GSP_PREC_RATE_avg',&
          &                   'average grid-scale surface precipitation rate',&
          &                   'kg/m**2/s', 61, 2,&
          &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
          &           k_jg)
        CASE (iecham,ildf_echam)
          !--- aprl ---
          CALL addVar(TimeVar('APRL',&
          &           'large-scale precip amount (rain + snow) '//&
          &           'accumulated over output interval',&
          &           'kg m-2', 142, 128,&
          &           vlistID(k_jg),gridCellID(k_jg),zaxisID_surface(k_jg)),&
          &           k_jg)
          !--- aprc ---
          CALL addVar(TimeVar('APRC',&
          &           'convective precip amount (rain + snow) '//&
          &           'accumulated over output interval',&
          &           'kg m-2', 143, 128,&
          &           vlistID(k_jg),gridCellID(k_jg),zaxisID_surface(k_jg)),&
          &           k_jg)
          !--- aprs ---
          CALL addVar(TimeVar('APRS',&
          &           'snowfall amount (large scale + convective) '//&
          &           'accumulated over output interval',&
          &           'kg m-2', 144, 128,&
          &           vlistID(k_jg),gridCellID(k_jg),zaxisID_surface(k_jg)),&
          &           k_jg)
          !--- rsfl ---
          CALL addVar(TimeVar('RSFL',&
          &                   'surface rain flux due to large scale condensation',&
          &                   'kg/m**2/s', 401, 999,&
          &                   vlistID(k_jg),gridCellID(k_jg),zaxisID_surface(k_jg)),&
          &           k_jg)
          !--- rsfc ---
          CALL addVar(TimeVar('RSFC',&
          &                   'surface rain flux due to convection',&
          &                   'kg/m**2/s', 402, 999,&
          &                   vlistID(k_jg),gridCellID(k_jg),zaxisID_surface(k_jg)),&
          &           k_jg)
          !--- ssfl ---
          CALL addVar(TimeVar('SSFL',&
          &                   'surface snow flux due to large scale condensation',&
          &                   'kg/m**2/s', 403, 999,&
          &                   vlistID(k_jg),gridCellID(k_jg),zaxisID_surface(k_jg)),&
          &           k_jg)
          !--- ssfc ---
          CALL addVar(TimeVar('SSFC',&
          &                   'surface snow flux due to convection',&
          &                   'kg/m**2/s', 404, 999,&
          &                   vlistID(k_jg),gridCellID(k_jg),zaxisID_surface(k_jg)),&
          &           k_jg)
        END SELECT !iforcing
      ENDIF !lwrite_precip

      ! cloud
      IF(lwrite_cloud ) THEN
        SELECT CASE (iforcing)
        CASE (inwp)
          !--- total water vapor---
          CALL addVar(TimeVar('QV',&
          &                   'total water vapor',&
          &                   '(0-1)', 91, 128,&
          &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_hybrid(k_jg)),&
          &           k_jg)
          !--- total cloud water ---
          CALL addVar(TimeVar('QC',&
          &                   'total cloud water',&
          &                   '(0-1)', 92, 128,&
          &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_hybrid(k_jg)),&
          &           k_jg)
          !--- total cloud ice ---
          CALL addVar(TimeVar('QI',&
          &                   'total cloud ice',&
          &                   '(0-1)', 93, 128,&
          &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_hybrid(k_jg)),&
          &           k_jg)
          !--- total cloud cover ---
          CALL addVar(TimeVar('CC',&
          &                   'total cloud cover',&
          &                   '(0-1)', 94, 128,&
          &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_hybrid(k_jg)),&
          &          k_jg)

          !--- vertically integrated total water vapor---
          CALL addVar(TimeVar('TQV',&
          &                   'vertically integrated total water vapor',&
          &                   'kg/m**2', 95, 128,&
          &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
          &           k_jg)
          !--- vertically integrated total cloud water ---
          CALL addVar(TimeVar('TQC',&
          &                   'vertically integrated total cloud water',&
          &                   'kk/m**2', 96, 128,&
          &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
          &           k_jg)
          !--- vertically integrated total cloud ice ---
          CALL addVar(TimeVar('TQI',&
          &                   'vertically integrated total cloud ice',&
          &                   'kg/m**2', 97, 128,&
          &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
          &           k_jg)
          !--- cloud cover assuming Maximum-Random overlap ---
          CALL addVar(TimeVar('TCC',&
          &                   'cloud cover assuming Maximum-Random overlap',&
          &                   '(0-1)', 98, 128,&
          &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
          &           k_jg)

          !--- average over the forcast time vertically integrated total water vapor---
          CALL addVar(TimeVar('TQV_avg',&
          &                   'average over forcast vertically integrated total water vapor',&
          &                   'km/m**2', 99, 128,&
          &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
          &           k_jg)
          !--- average over the forcast time vertically integrated total cloud water ---
          CALL addVar(TimeVar('TQC_avg',&
          &                   'average over forcast vertically integrated total cloud water',&
          &                   'kg/m**2',101, 128,&
          &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
          &           k_jg)
          !--- average over the forcast time vertically integrated total cloud ice ---
          CALL addVar(TimeVar('TQI_avg',&
          &                   'average over forcast vertically integrated total cloud ice',&
          &                   'kg/m**2',102, 128,&
          &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
          &           k_jg)
          !--- average over the forecast time of the  cloud cover ---
          CALL addVar(TimeVar('TCC_avg',&
          &                   'average over the forecast time of the cloud cover',&
          &                   '(0-1)', 103, 128,&
          &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
          &           k_jg)

        CASE (iecham,ildf_echam)

          !--- aclcov ---
        CALL addVar(TimeVar('ACLCOV',&
          &                 'total cloud cover accumulated over output interval',&
          &                 's', 164, 128,&
          &                 vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),k_jg)
          !--- aclc ---
        CALL addVar(TimeVar('ACLC',&
          &                 'cloud area fraction (instantaneous)',&
          &                 '(0-1)', 162, 128,&
          &                 vlistID(k_jg), gridCellID(k_jg),zaxisID_hybrid(k_jg)),k_jg)
          !--- aclcac ---
        CALL addVar(TimeVar('ACLCAC',&
          &                 'cloud area fraction accumulated over output interval',&
          &                 's', 223, 128,&
          &                 vlistID(k_jg), gridCellID(k_jg),zaxisID_hybrid(k_jg)),k_jg)
          !--- qvi ---
        CALL addVar(TimeVar('qvi',&
          &                 'temporally and vertically integrated water vapor content',&
          &                 's kg/m**2', 230, 128,&
          &                 vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),k_jg)
          !--- xlvi ---
        CALL addVar(TimeVar('xlvi',&
          &                 'temporally and vertically integrated cloud water content',&
          &                 's kg/m**2', 231, 128,&
          &                 vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),k_jg)
          !--- xivi ---
        CALL addVar(TimeVar('xivi',&
          &                 'temporally and vertically integrated cloud ice content',&
          &                 's kg/m**2', 232, 128,&
          &                 vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),k_jg)
       !  !--- omega (for debugging) ---
       !CALL addVar(TimeVar('OMEGA_PHY',&
       !  &                 'vertical velocity in pressure coordinate',&
       !  &                 'Pa/s', 135, 128,&
       !  &                 vlistID(k_jg), gridCellID(k_jg), zaxisID_hybrid(k_jg)),k_jg)

        END SELECT !iforcing
      ENDIF !lwrite_cloud

     ! TKE
      IF(lwrite_tke ) THEN
        SELECT CASE (iforcing)
        CASE (inwp)
          !--- turbulent konetic energy---
          CALL addVar(TimeVar('TKE',&
          &                   'turbulent kinetic energy',&
          &                   'm^2/s^2', 152, 201,&
          &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_hybrid_half(k_jg)),&
          &           k_jg)
        END SELECT
      ENDIF !lwrite_tke

      IF(lwrite_surface ) THEN
        SELECT CASE (iforcing)
        CASE (inwp)
          !--- roughness length---
          CALL addVar(TimeVar ('Z0',&
               &              'roughness length',&
          &                   'm', 83, 2,&
          &                    vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
          &                    k_jg)
          !--- Temperature at Surface---
          CALL addVar(TimeVar('T_G',&
          &                   'aggregated surface temperature',&
          &                   'K', 11, 201,&
          &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
          &                   k_jg)

        IF ( atm_phy_nwp_config(k_jg)%inwp_surface == 1 ) THEN  ! TERRA
          !--- Weighted temperature at surface---
          DO jt = 1, nsfc_subs
            WRITE(cjt,'(i2)') jt
            WRITE(name,'(A,A)') "T_GT_tile_", TRIM(ADJUSTL(cjt))
            WRITE(long_name,'(A,A)') "weighted surface temperature tile ",TRIM(ADJUSTL(cjt))
            CALL addVar(TimeVar(TRIM(name),TRIM(long_name),&
            &          'K',11,201,&
            &           vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
            &           k_jg)
          ENDDO

          DO jt = 1, nsfc_subs
            WRITE(cjt,'(i2)') jt
            WRITE(name,'(A,A)') "T_SNOW_tile_", TRIM(ADJUSTL(cjt))
            WRITE(long_name,'(A,A)') "temperature of the snow-surface tile ",TRIM(ADJUSTL(cjt))
            CALL addVar(TimeVar(TRIM(name),TRIM(long_name),&
            &          'K',11,201,&
            &           vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
            &           k_jg)
          ENDDO

          DO jt = 1, nsfc_subs
            WRITE(cjt,'(i2)') jt
            WRITE(name,'(A,A)') "T_SNOW_MULT_tile_", TRIM(ADJUSTL(cjt))
            WRITE(long_name,'(A,A)') "temperature of the snow-surface tile ",TRIM(ADJUSTL(cjt))
            CALL addVar(TimeVar(TRIM(name),TRIM(long_name),&
            &          'K',11,201,&
            &           vlistID(k_jg), gridCellID(k_jg),zaxisID_generic_snow_p1(k_jg)),&
            &           k_jg)
          ENDDO

          DO jt = 1, nsfc_subs
            WRITE(cjt,'(i2)') jt
            WRITE(name,'(A,A)') "T_S_tile_", TRIM(ADJUSTL(cjt))
            WRITE(long_name,'(A,A)') "temperature of ground surface tile ",TRIM(ADJUSTL(cjt))
            CALL addVar(TimeVar(TRIM(name),TRIM(long_name),&
            &          'K',11,201,&
            &           vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
            &           k_jg)
          ENDDO

          DO jt = 1, nsfc_subs
            WRITE(cjt,'(i2)') jt
            WRITE(name,'(A,A)') "W_SNOW_tile_", TRIM(ADJUSTL(cjt))
            WRITE(long_name,'(A,A)') "water content of snow tile ",TRIM(ADJUSTL(cjt))
            CALL addVar(TimeVar(TRIM(name),TRIM(long_name),&
            &          'm H2O',11,201,&
            &           vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
            &           k_jg)
          ENDDO

          DO jt = 1, nsfc_subs
            WRITE(cjt,'(i2)') jt
            WRITE(name,'(A,A)') "RHO_SNOW_tile_", TRIM(ADJUSTL(cjt))
            WRITE(long_name,'(A,A)') "snow density tile ",TRIM(ADJUSTL(cjt))
            CALL addVar(TimeVar(TRIM(name),TRIM(long_name),&
            &          'kg/m**3',11,201,&
            &           vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
            &           k_jg)
          ENDDO

          DO jt = 1, nsfc_subs
            WRITE(cjt,'(i2)') jt
            WRITE(name,'(A,A)') "RHO_SNOW_MULT_tile_", TRIM(ADJUSTL(cjt))
            WRITE(long_name,'(A,A)') "snow density tile ",TRIM(ADJUSTL(cjt))
            CALL addVar(TimeVar(TRIM(name),TRIM(long_name),&
            &          'kg/m**3',11,201,&
            &           vlistID(k_jg), gridCellID(k_jg),zaxisID_generic_snow(k_jg)),&
            &           k_jg)
          ENDDO

          DO jt = 1, nsfc_subs
            WRITE(cjt,'(i2)') jt
            WRITE(name,'(A,A)') "W_I_tile_", TRIM(ADJUSTL(cjt))
            WRITE(long_name,'(A,A)') "water content of interception water tile",TRIM(ADJUSTL(cjt))
            CALL addVar(TimeVar(TRIM(name),TRIM(long_name),&
            &          'm H2O',11,201,&
            &           vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
            &           k_jg)

          ENDDO

          DO jt = 1, nsfc_subs
            WRITE(cjt,'(i2)') jt
            WRITE(name,'(A,A)') "T_SO_tile_", TRIM(ADJUSTL(cjt))
            WRITE(long_name,'(A,A)') "soil temperature (main level) tile",TRIM(ADJUSTL(cjt))
            CALL addVar(TimeVar(TRIM(name),TRIM(long_name),&
            &          'K',11,201,&
            &           vlistID(k_jg), gridCellID(k_jg),zaxisID_depth_below_land_p1(k_jg)),&
            &           k_jg)

          ENDDO

          DO jt = 1, nsfc_subs
            WRITE(cjt,'(i2)') jt
            WRITE(name,'(A,A)') "W_SO_tile_", TRIM(ADJUSTL(cjt))
            WRITE(long_name,'(A,A)') "total water content (ice + liquid water) tile", &
            &                        TRIM(ADJUSTL(cjt))
            CALL addVar(TimeVar(TRIM(name),TRIM(long_name),&
            &          'm H2O',11,201,&
            &           vlistID(k_jg), gridCellID(k_jg),zaxisID_depth_below_land(k_jg)),&
            &           k_jg)

          ENDDO

          DO jt = 1, nsfc_subs
            WRITE(cjt,'(i2)') jt
            WRITE(name,'(A,A)') "W_SO_ICE_tile_", TRIM(ADJUSTL(cjt))
            WRITE(long_name,'(A,A)') "ice content tile",TRIM(ADJUSTL(cjt))
            CALL addVar(TimeVar(TRIM(name),TRIM(long_name),&
            &          'm H2O',11,201,&
            &           vlistID(k_jg), gridCellID(k_jg),zaxisID_depth_below_land(k_jg)),&
            &           k_jg)

          ENDDO

          DO jt = 1, nsfc_subs
            WRITE(cjt,'(i2)') jt
            WRITE(name,'(A,A)') "WLIQ_SNOW_tile_", TRIM(ADJUSTL(cjt))
            WRITE(long_name,'(A,A)') "liquid water content in snow tile",TRIM(ADJUSTL(cjt))
            CALL addVar(TimeVar(TRIM(name),TRIM(long_name),&
            &          'm H2O',11,201,&
            &           vlistID(k_jg), gridCellID(k_jg),zaxisID_generic_snow(k_jg)),&
            &           k_jg)

          ENDDO

          DO jt = 1, nsfc_subs
            WRITE(cjt,'(i2)') jt
            WRITE(name,'(A,A)') "WTOT_SNOW_tile_", TRIM(ADJUSTL(cjt))
            WRITE(long_name,'(A,A)') "total water content in snow tile",TRIM(ADJUSTL(cjt))
            CALL addVar(TimeVar(TRIM(name),TRIM(long_name),&
            &          'm H2O',11,201,&
            &           vlistID(k_jg), gridCellID(k_jg),zaxisID_generic_snow(k_jg)),&
            &           k_jg)

          ENDDO

          DO jt = 1, nsfc_subs
            WRITE(cjt,'(i2)') jt
            WRITE(name,'(A,A)') "H_SNOW_tile_", TRIM(ADJUSTL(cjt))
            WRITE(long_name,'(A,A)') "snow height", &
            & TRIM(ADJUSTL(cjt))
            CALL addVar(TimeVar(TRIM(name),TRIM(long_name),&
            &          'm',11,201,&
            &           vlistID(k_jg), gridCellID(k_jg),zaxisID_generic_snow(k_jg)),&
            &           k_jg)

          ENDDO

          DO jt = 1, nsfc_subs
            WRITE(cjt,'(i2)') jt
            WRITE(name,'(A,A)') "DZH_SNOW_tile_", TRIM(ADJUSTL(cjt))
            WRITE(long_name,'(A,A)') "layer thickness between half levels in snow tile", &
            & TRIM(ADJUSTL(cjt))
            CALL addVar(TimeVar(TRIM(name),TRIM(long_name),&
            &          'm',11,201,&
            &           vlistID(k_jg), gridCellID(k_jg),zaxisID_generic_snow(k_jg)),&
            &           k_jg)

          ENDDO
          END IF  ! inwp_surface

          !--- Specific Humidity at Surface---
          CALL addVar(TimeVar('QV_S',&
          &                   'aggregated surface specific humidity',&
          &                   'K', 51, 201,&
          &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
          &                   k_jg)
          !--- Fluxes ....---
          CALL addVar(TimeVar('SHFL_S_avg',&
          &                   'averaged sensible heat flux at surface',&
          &                   'W/m^2', 122, 201,&  !999 == WMO here
          &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
          &                   k_jg)
          CALL addVar(TimeVar('LHFL_S_avg',&
          &                   'averaged latent heat flux at surface',&
          &                   'W/m^2', 121, 201,& !999 ==  WMO here
          &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
          &                   k_jg)
          CALL addVar(TimeVar('EVAP_RATE_avg',&
          &                   'averaged moisture flux rate (evap. rate) at surface',&
          &                   'kg/m*2/s', 121, 201,& !999 ==  WMO here
          &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
          &                   k_jg)
          CALL addVar(TimeVar('PS_s6avg',&
          &                   '6 hourly sample surface pressure average',&
          &                   'Pa', 134, 128,&
          &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
          &           k_jg)
          CALL addVar(TimeVar('T_2m',&
          &                   '2 m Temperature',&
          &                   'K', 134, 128,&
          &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
          &           k_jg)
          CALL addVar(TimeVar('T_2m_s6avg',&
          &                   '6 hourly sample 2 m Temperature average',&
          &                   'K', 134, 128,&
          &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
          &           k_jg)
          CALL addVar(TimeVar('QV_2m',&
          &                   '2 m specific humidity ',&
          &                   'Kg/kg', 134, 128,&
          &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
          &           k_jg)
          CALL addVar(TimeVar('QV_2m_s6avg',&
          &                   '6 hourly sample 2 m specific humidity average',&
          &                   'Kg/kg', 134, 128,&
          &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
          &           k_jg)
          CALL addVar(TimeVar('U_10m',&
          &                   '10 m zonal wind',&
          &                   'K', 134, 128,&
          &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
          &           k_jg)
          CALL addVar(TimeVar('V_10m',&
          &                   '10 m meridional wind',&
          &                   'K', 134, 128,&
          &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
          &           k_jg)
          CALL addVar(TimeVar('U_10m_s6avg',&
          &                   '6 hourly sample 10 m zonal wind average',&
          &                   'K', 134, 128,&
          &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
          &           k_jg)
          CALL addVar(TimeVar('V_10m_s6avg',&
          &                   '6 hourly sample 10 m meridional wind average',&
          &                   'K', 134, 128,&
          &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
          &           k_jg)


        CASE (iecham,ildf_echam)

          CALL addVar(TimeVar('evap_ac',                                         &
          &           'evaporation accumulated over output interval',            &
          &           'kg m-2', 182, 128,                                        &
          &           vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),k_jg)

          CALL addVar(TimeVar('lhflx_ac',                                        &
          &           'latent heat flux accumulated over output interval',       &
          &           'W m-2 s', 147, 128,                                       &
          &           vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),k_jg)

          CALL addVar(TimeVar('shflx_ac',                                        &
          &           'sensible heat flux accumulated over output interval',     &
          &           'W m-2 s', 146, 128,                                       &
          &           vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),k_jg)

          CALL addVar(TimeVar('u_stress_ac',                                     &
          &           'surface wind stress accumulated over output interval',    &
          &           'N m-2 s', 180, 128,                                       &
          &           vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),k_jg)

          CALL addVar(TimeVar('v_stress_ac',                                     &
          &           'surface wind stress accumulated over output interval',    &
          &           'N m-2 s', 181, 128,                                       &
          &           vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),k_jg)


          DO jt = 1,nsfc_type 

            WRITE(name,'(a,i1)') "evap_tile_",jt
            CALL addVar(TimeVar(TRIM(name), TRIM(name),'kg m-2 s-1',182,128,       &
            &           vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),k_jg)

            WRITE(name,'(a,i1)') "lhflx_tile_",jt
            CALL addVar(TimeVar(TRIM(name), TRIM(name),'W m-2',147,128,            &
            &           vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),k_jg)

            WRITE(name,'(a,i1)') "shflx_tile_",jt
            CALL addVar(TimeVar(TRIM(name), TRIM(name),'W m-2',146,128,            &
            &           vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),k_jg)

            WRITE(name,'(a,i1)') "u_stress_tile_",jt
            CALL addVar(TimeVar(TRIM(name), TRIM(name),'N m-2',180,128,            &
            &           vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),k_jg)

            WRITE(name,'(a,i1)') "v_stress_tile_",jt
            CALL addVar(TimeVar(TRIM(name), TRIM(name),'N m-2',181,128,            &
            &           vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),k_jg)

          END DO

        END SELECT !iforcing
      ENDIF !lwrite_surface

      ! Tendencies induced by physics parameterizations
      IF (lwrite_tend_phy) THEN
        SELECT CASE (iforcing)
        CASE (iecham,ildf_echam)
          ! Temperature tendencies
          CALL addVar(TimeVar('tend_temp_radsw',&
          &                   'shortwave radiative heating', &
          &                   'K/s', 101, 999,&
          &                   vlistID(k_jg),gridCellID(k_jg),zaxisID_hybrid(k_jg)),&
          &           k_jg)
          CALL addVar(TimeVar('tend_temp_radlw', &
          &                   'longwave radiative heating',&
          &                   'K/s', 102, 999,&
          &                   vlistID(k_jg),gridCellID(k_jg),zaxisID_hybrid(k_jg)),&
          &           k_jg)
          CALL addVar(TimeVar('tend_temp_cld', &
          &                   'temperature tendency caused by large scale condensation',&
          &                   'K/s', 103, 999,&
          &                   vlistID(k_jg),gridCellID(k_jg),zaxisID_hybrid(k_jg)),&
          &           k_jg)
          CALL addVar(TimeVar('tend_temp_cnv', &
          &                   'temperature tendency caused by cumulus convection',&
          &                   'K/s', 104, 999,&
          &                   vlistID(k_jg),gridCellID(k_jg),zaxisID_hybrid(k_jg)),&
          &           k_jg)
          CALL addVar(TimeVar('tend_temp_vdf', &
          &                   'temperature tendency caused by turbulent mixing',&
          &                   'K/s', 105, 999,&
          &                   vlistID(k_jg),gridCellID(k_jg),zaxisID_hybrid(k_jg)),&
          &           k_jg)
          CALL addVar(TimeVar('tend_temp_gwh', &
          &                   'temperature tendency caused by gravity wave dissipation',&
          &                   'K/s', 106, 999,&
          &                   vlistID(k_jg),gridCellID(k_jg),zaxisID_hybrid(k_jg)),&
          &           k_jg)

          ! u-wind tendency
          CALL addVar(TimeVar('tend_u_cnv',&
          &                   'zonal wind tendency caused by cumulus convection',&
          &                   'm/s2', 111, 999,&
          &                   vlistID(k_jg),gridCellID(k_jg),zaxisID_hybrid(k_jg)),&
          &           k_jg)
          CALL addVar(TimeVar('tend_u_vdf',&
          &                   'zonal wind tendency caused by turbulent mixing',&
          &                   'm/s2', 112, 999,&
          &                   vlistID(k_jg),gridCellID(k_jg),zaxisID_hybrid(k_jg)),&
          &           k_jg)
          CALL addVar(TimeVar('tend_u_gwh',&
          &                   'zonal wind tendency caused by gravity wave dissipation',&
          &                   'm/s2', 113, 999,&
          &                   vlistID(k_jg),gridCellID(k_jg),zaxisID_hybrid(k_jg)),&
          &           k_jg)

          ! v-wind tendency
          CALL addVar(TimeVar('tend_v_cnv',&
          &                   'meridional wind tendency caused by cumulus convection',&
          &                   'm/s2', 121, 999,&
          &                   vlistID(k_jg),gridCellID(k_jg),zaxisID_hybrid(k_jg)),&
          &           k_jg)
          CALL addVar(TimeVar('tend_v_vdf',&
          &                   'meridional wind tendency caused by turbulent mixing',&
          &                   'm/s2', 122, 999,&
          &                   vlistID(k_jg),gridCellID(k_jg),zaxisID_hybrid(k_jg)),&
          &           k_jg)
          CALL addVar(TimeVar('tend_v_gwh',&
          &                   'meridional wind tendency caused by gravity wave dissipation',&
          &                   'm/s2', 123, 999,&
          &                   vlistID(k_jg),gridCellID(k_jg),zaxisID_hybrid(k_jg)),&
          &           k_jg)

        ! Tracer tendencies

        IF (ntracer > 0) THEN
          DO jt = 1, iqcond
            ctracer = ctracer_list(jt:jt)
            WRITE(name,'(A6,A1,A4)') "tend_q", ctracer, "_cnv"
            CALL addVar(TimeVar(TRIM(name), &
                &               TRIM(name), &
                &               '1/s',22,128,&
                &               vlistID(k_jg), gridCellID(k_jg),zaxisID_hybrid(k_jg)),&
                &       k_jg)
          END DO

          DO jt = 1, iqcond
            ctracer = ctracer_list(jt:jt)
            WRITE(name,'(A6,A1,A4)') "tend_q", ctracer, "_vdf"
            CALL addVar(TimeVar(TRIM(name),&
            &               TRIM(name),&
            &               '1/s',23,128,&
            &               vlistID(k_jg), gridCellID(k_jg),zaxisID_hybrid(k_jg)),&
            &       k_jg)
          END DO
        ENDIF

        CASE (inwp)
          ! Temperature tendencies
          CALL addVar(TimeVar('tend_temp_radsw',&
          &                   'shortwave radiative heating', &
          &                   'K/s', 101, 999,&
          &                   vlistID(k_jg),gridCellID(k_jg),zaxisID_hybrid(k_jg)),&
          &           k_jg)
          CALL addVar(TimeVar('tend_temp_radlw', &
          &                   'longwave radiative heating',&
          &                   'K/s', 102, 999,&
          &                   vlistID(k_jg),gridCellID(k_jg),zaxisID_hybrid(k_jg)),&
          &           k_jg)
          CALL addVar(TimeVar('tend_temp_pscl', &
          &                   'temperature tendency due to grid scale microphysics',&
          &                   'K/s', 103, 999,&
          &                   vlistID(k_jg),gridCellID(k_jg),zaxisID_hybrid(k_jg)),&
          &           k_jg)
          CALL addVar(TimeVar('tend_temp_conv', &
          &                   'temperature tendency due to convection',&
          &                   'K/s', 104, 999,&
          &                   vlistID(k_jg),gridCellID(k_jg),zaxisID_hybrid(k_jg)),&
          &           k_jg)
          CALL addVar(TimeVar('tend_temp_turb', &
          &                   'temperature tendency due to turbulence',&
          &                   'K/s', 105, 999,&
          &                   vlistID(k_jg),gridCellID(k_jg),zaxisID_hybrid(k_jg)),&
          &           k_jg)
          CALL addVar(TimeVar('tend_temp_sso', &
          &                   'temperature tendency due to SSO',&
          &                   'K/s', 106, 999,&
          &                   vlistID(k_jg),gridCellID(k_jg),zaxisID_hybrid(k_jg)),&
          &           k_jg)
          CALL addVar(TimeVar('tend_temp_dyn', &
          &                   'temperature tendency due to dynamics',&
          &                   'K/s', 107, 999,&
          &                   vlistID(k_jg),gridCellID(k_jg),zaxisID_hybrid(k_jg)),&
          &           k_jg)

          ! u-wind tendency
          CALL addVar(TimeVar('tend_u_conv',&
          &                   'zonal wind tendency due to convection',&
          &                   'm/s2', 111, 999,&
          &                   vlistID(k_jg),gridCellID(k_jg),zaxisID_hybrid(k_jg)),&
          &           k_jg)
          CALL addVar(TimeVar('tend_u_turb',&
          &                   'zonal wind tendency due to turbulence',&
          &                   'm/s2', 112, 999,&
          &                   vlistID(k_jg),gridCellID(k_jg),zaxisID_hybrid(k_jg)),&
          &           k_jg)
          CALL addVar(TimeVar('tend_u_sso',&
          &                   'zonal wind tendency due to SSO',&
          &                   'm/s2', 121, 999,&
          &                   vlistID(k_jg),gridCellID(k_jg),zaxisID_hybrid(k_jg)),&
          &           k_jg)


          ! v-wind tendency
          CALL addVar(TimeVar('tend_v_conv',&
          &                   'meridional wind tendency due to convection',&
          &                   'm/s2', 122, 999,&
          &                   vlistID(k_jg),gridCellID(k_jg),zaxisID_hybrid(k_jg)),&
          &           k_jg)
          CALL addVar(TimeVar('tend_v_turb',&
          &                   'meridional wind tendency due to turbulence',&
          &                   'm/s2', 123, 999,&
          &                   vlistID(k_jg),gridCellID(k_jg),zaxisID_hybrid(k_jg)),&
          &           k_jg)
          CALL addVar(TimeVar('tend_v_sso',&
          &                   'meridional wind tendency due to SSO',&
          &                   'm/s2', 124, 999,&
          &                   vlistID(k_jg),gridCellID(k_jg),zaxisID_hybrid(k_jg)),&
          &           k_jg)

        ! Tracer tendencies

        IF (ntracer > 0) THEN
          DO jt = 1, iqcond
            ctracer = ctracer_list(jt:jt)
            WRITE(qname,'(A6,A1,A5)') "tend_q",ctracer, "_conv"
            CALL addVar(TimeVar(TRIM(qname), &
                &               TRIM(qname), &
                &               '1/s',22,128,&
                &               vlistID(k_jg), gridCellID(k_jg),zaxisID_hybrid(k_jg)),&
                &       k_jg)
          END DO

          DO jt = 1, iqcond
            ctracer = ctracer_list(jt:jt)
            WRITE(qname,'(A6,A1,A5)') "tend_q",ctracer, "_turb"
            CALL addVar(TimeVar(TRIM(qname),&
            &               TRIM(qname),&
            &               '1/s',23,128,&
            &               vlistID(k_jg), gridCellID(k_jg),zaxisID_hybrid(k_jg)),&
            &       k_jg)
          END DO

  ! KF to be implemented later
  !        DO jt = 1,  iqcond
  !          ctracer = ctracer_list(jt:jt)
  !          WRITE(qname,'(A6,A1,A5)') "tend_q", ctracer, "_pscl"
  !          CALL addVar(TimeVar(TRIM(qname),&
  !          &                   TRIM(qname),&
  !          &               '1/s',24,128,&
  !          &               vlistID(k_jg), gridCellID(k_jg),zaxisID_hybrid(k_jg)),&
  !          &       k_jg)
  !        END DO
        END IF

      END SELECT
      ENDIF !lwrite_tend_phy

      ! Vorticity
      IF (lwrite_vorticity) THEN
        CALL addVar(TimeVar('VOR',&
        &                   'vorticity',&
        &                   '1/s',138,128, &
        &                   vlistID(k_jg), gridVertexID(k_jg), zaxisID_hybrid(k_jg)),&
        &           k_jg)
      END IF

      IF (lwrite_divergence) THEN
        CALL addVar(TimeVar('DIV',&
        &                   'divergence',&
        &                   '1/s',155,128,&
        &                    vlistID(k_jg),gridCellID(k_jg),zaxisID_hybrid(k_jg)),&
        &           k_jg)
      END IF

      IF (iequations == 3) THEN
        ! virtual potential temperature
        CALL addVar(TimeVar('THETA_V', &
        &                   'virtual potential temperature',&
        &                   'K',192, 128, &
        &                   vlistID(k_jg),gridCellID(k_jg),zaxisID_hybrid(k_jg)),&
        &          k_jg)

        ! Exner pressure
        CALL addVar(TimeVar('EXNER',&
        &                   'Exner pressure', &
        &                   '-', 193, 128, &
        &                   vlistID(k_jg),gridCellID(k_jg),zaxisID_hybrid(k_jg)),&
        &          k_jg)

        ! Density
        CALL addVar(TimeVar('RHO', &
        &                   'density', &
        &                   'kg/m**3', 194, 128,&
        &                   vlistID(k_jg),gridCellID(k_jg),zaxisID_hybrid(k_jg)), &
        &           k_jg)

      ENDIF


      ! output on constant z- and/or p-levels
      !
      IF(lwrite_pzlev .AND. iequations==inh_atmosphere) THEN
        IF (nh_pzlev_config(k_jg)%lwrite_zlev) THEN
          ! zonal wind
          CALL addVar(TimeVar('U_Z',&
          &                   'zonal wind',&
          &                   'm/s', 131, 128,&
          &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_hgt(k_jg)),&
          &           k_jg)
          ! meridional wind
          CALL addVar(TimeVar('V_Z',&
          &                   'meridional wind',&
          &                   'm/s', 132, 128, &
          &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_hgt(k_jg)),&
          &           k_jg)
          CALL addVar(TimeVar('T_Z', &
          &                   'temperature',&
          &                   'K', 130, 128,&
          &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_hgt(k_jg)),&
          &           k_jg)
          CALL addVar(TimeVar('P_Z',&
          &                   'pressure',&
          &                   'Pa', 54, 128,&
          &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_hgt(k_jg)),&
          &          k_jg)


          ! prognostic tracer fields
          IF (ntracer > 0 .AND. (iforcing == inwp .OR. ltransport)) THEN
            DO jt = 1, ntracer
              IF (lwrite_tracer(jt)) THEN
                ctracer = ctracer_list(jt:jt)
                IF (jt.EQ.1) THEN !for water vapour
                  elemid=51; tableid=2   !DWD coding
                ELSEIF (jt.EQ.2) THEN !for cloud water
                 elemid=31; tableid=201 !DWD coding
                ELSEIF (jt.EQ.3) THEN !for cloud ice
                 elemid=33; tableid=201 !DWD coding
                ELSE !other tracers
                 elemid=252; tableid=128 !default coding
                END IF
                WRITE(name,'(A1,A1,A2)') "Q", ctracer,"_Z"
                CALL addVar(TimeVar(TRIM(name),TRIM(name),&
                  &                 'kg/kg',elemid,tableid,&
                  &                 vlistID(k_jg), gridCellID(k_jg),zaxisID_hgt(k_jg)),&
                  &         k_jg)          
              END IF
            END DO

            IF( lwrite_cloud ) THEN
              CALL addVar(TimeVar('QV_Z',&
              &                   'total water vapor',&
              &                   'kg/kg', 91, 128,&
              &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_hgt(k_jg)),&
              &          k_jg)
              CALL addVar(TimeVar('QC_Z',&
              &                   'total specific cloud water content',&
              &                   'kg/kg', 246, 128,&
              &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_hgt(k_jg)),&
              &          k_jg)
              CALL addVar(TimeVar('QI_Z',&
              &                   'total specific cloud ice content',&
              &                   'kg/kg', 247, 128,&
              &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_hgt(k_jg)),&
              &          k_jg)
              CALL addVar(TimeVar('CC_Z',&
              &                   'total_cloud_cover',&
              &                   '0-1', 248, 128,&
              &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_hgt(k_jg)),&
              &          k_jg)
            ENDIF
          ENDIF  ! ntracer

        ENDIF ! lwrite_zlev
        IF (nh_pzlev_config(k_jg)%lwrite_plev) THEN
          ! zonal wind
          CALL addVar(TimeVar('U_P',&
          &                   'zonal wind',&
          &                   'm/s', 131, 128,&
          &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_pres(k_jg)),&
          &           k_jg)
          ! meridional wind
          CALL addVar(TimeVar('V_P',&
          &                   'meridional wind',&
          &                   'm/s', 132, 128, &
          &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_pres(k_jg)),&
          &           k_jg)
          CALL addVar(TimeVar('T_P', &
          &                   'temperature',&
          &                   'K', 130, 128,&
          &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_pres(k_jg)),&
          &           k_jg)
          CALL addVar(TimeVar('Z', &
          &                   'geopotential',&
          &                   'm2/s2', 129, 128,&
          &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_pres(k_jg)),&
          &           k_jg)


          ! prognostic tracer fields
          IF (ntracer > 0 .AND. (iforcing == inwp .OR. ltransport)) THEN
            DO jt = 1, ntracer
              IF (lwrite_tracer(jt)) THEN
                ctracer = ctracer_list(jt:jt)
                IF (jt.EQ.1) THEN !for water vapour
                  elemid=51; tableid=2   !DWD coding
                ELSEIF (jt.EQ.2) THEN !for cloud water
                 elemid=31; tableid=201 !DWD coding
                ELSEIF (jt.EQ.3) THEN !for cloud ice
                 elemid=33; tableid=201 !DWD coding
                ELSE !other tracers
                 elemid=252; tableid=128 !default coding
                END IF
                WRITE(name,'(A1,A1,A2)') "Q", ctracer,"_P"
                CALL addVar(TimeVar(TRIM(name),TRIM(name),&
                  &                 'kg/kg',elemid,tableid,&
                  &                 vlistID(k_jg), gridCellID(k_jg),zaxisID_pres(k_jg)),&
                  &         k_jg)          
              END IF
            END DO

            IF( lwrite_cloud ) THEN
              CALL addVar(TimeVar('QV_P',&
              &                   'total water vapor',&
              &                   'kg/kg', 91, 128,&
              &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_pres(k_jg)),&
              &          k_jg)
              CALL addVar(TimeVar('QC_P',&
              &                   'total specific cloud water content',&
              &                   'kg/kg', 246, 128,&
              &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_pres(k_jg)),&
              &          k_jg)
              CALL addVar(TimeVar('QI_P',&
              &                   'total specific cloud ice content',&
              &                   'kg/kg', 247, 128,&
              &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_pres(k_jg)),&
              &          k_jg)
              CALL addVar(TimeVar('CC_P',&
              &                   'total cloud cover',&
              &                   '0-1', 248, 128,&
              &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_pres(k_jg)),&
              &          k_jg)
            ENDIF

          ENDIF  ! ntracer


        ENDIF ! lwrite_plev
      ENDIF !lwrite_pzlev

    ELSE 
      ! ocean
      ! 3-dim lsm-masks
      CALL addVar(ConstVar('wet_c',&
      &                    '3d lsm on cells',&
      &                    '', 1, 128,&
      &                     vlistID(k_jg),&
      &                     gridCellID(k_jg), &
      &                     zaxisIDdepth_m(k_jg)),&
      &           k_jg)
      CALL addVar(ConstVar('wet_e',&
      &                    '3d lsm on edges',&
      &                    '', 1, 128,&
      &                     vlistID(k_jg),&
      &                     gridEdgeID(k_jg), &
      &                     zaxisIDdepth_m(k_jg)),&
      &           k_jg)
      CALL addVar(TimeVar('ELEV',&
      &                   'surface elevation at cell center',&
      &                   'm', 1, 128,&
      &                    vlistID(k_jg),&
      &                    gridCellID(k_jg), &
      &                    zaxisID_surface(k_jg)),&
      &           k_jg)
    !  CALL addVar(TimeVar('forc-u',&
    !  &                   'u-forcing component at centers',&
    !  &                   'N/m2',13,128,&
    !  &                   vlistID(k_jg),&
    !  &                   gridCellID(k_jg),&
    !  &                   zaxisID_surface(k_jg)),&
    !  &           k_jg)
    !  CALL addVar(TimeVar('forc-v',&
    !  &                   'v-forcing component at centers',&
    !  &                   'N/m2',14,128,&
    !  &                   vlistID(k_jg),&
    !  &                   gridCellID(k_jg),&
    !  &                   zaxisID_surface(k_jg)),&
    !  &           k_jg)
    !  CALL addVar(TimeVar('horz-adv',&
    !  &                   'nonlin Cor ',&
    !  &                   'm/s',2,128,&
    !  &                   vlistID(k_jg),&
    !  &                   gridEdgeID(k_jg), &
    !  &                   zaxisIDdepth_m(k_jg)),&
    !  &           k_jg)
    !  
    !  CALL addVar(TimeVar('Ekin-grad',&
    !  &                   'gradient Ekin ',&
    !  &                   'm/s',2,128,&
    !  &                   vlistID(k_jg),&
    !  &                   gridEdgeID(k_jg), &
    !  &                   zaxisIDdepth_m(k_jg)),&
    !  &           k_jg)
    !  CALL addVar(TimeVar('flux-u',&
    !  &                   'sfc_flux u at centers',&
    !  &                   'N/m2',13,128,&
    !  &                   vlistID(k_jg),&
    !  &                   gridCellID(k_jg),&
    !  &                   zaxisID_surface(k_jg)),&
    !  &           k_jg)
    !  CALL addVar(TimeVar('flux-v',&
    !  &                   'sfc_flux v at centers',&
    !  &                   'N/m2',14,128,&
    !  &                   vlistID(k_jg),&
    !  &                   gridCellID(k_jg),&
    !  &                   zaxisID_surface(k_jg)),&
    !  &           k_jg)
    !  
    !  CALL addVar(TimeVar('flux-VN',&
    !  &                   'sfc-flux at edge',&
    !  &                   'm/s',2,128,&
    !  &                   vlistID(k_jg),&
    !  &                   gridEdgeID(k_jg), &
    !  &                   zaxisID_surface(k_jg)),&
    !  &           k_jg)
    !  
    !  CALL addVar(TimeVar('VN',&
    !  &                   'velocity normal to edge',&
    !  &                   'm/s',2,128,&
    !  &                   vlistID(k_jg),&
    !  &                   gridEdgeID(k_jg), &
    !  &                   zaxisIDdepth_m(k_jg)),&
    !  &           k_jg)
      CALL addVar(TimeVar('VORT',&
      &                   'vorticity at vertices',&
      &                   '1/s',3,128,&
      &                   vlistID(k_jg),&
      &                   gridVertexID(k_jg),&
      &                   zaxisIDdepth_m(k_jg)),&
      &          k_jg)
      CALL addVar(TimeVar('u-veloc',&
      &                   'u-velocity component at centers',&
      &                   'm/s',4,128,&
      &                   vlistID(k_jg),&
      &                   gridCellID(k_jg),&
      &                   zaxisIDdepth_m(k_jg)),&
      &           k_jg)
      CALL addVar(TimeVar('v-veloc',&
      &                   'v-velocity component at centers',&
      &                   'm/s',5,128,&
      &                   vlistID(k_jg), &
      &                   gridCellID(k_jg), &
      &                   zaxisIDdepth_m(k_jg)),&
      &           k_jg)
      CALL addVar(TimeVar('W',&
      &                   'vertical velocity at cells',&
      &                   'm/s', 6, 128,&
      &                   vlistid(k_jg),&
      &                   gridcellid(k_jg),&
      &                   zaxisid_halfdepth(k_jg)),&
      &           k_jg)
 !     CALL addVar(TimeVar('press-grad',&
 !     &                   'pressure-gradient at edges',&
 !     &                   'm/s',5,128,&
 !     &                   vlistID(k_jg), &
 !     &                   gridEdgeID(k_jg), &
 !     &                   zaxisIDdepth_m(k_jg)),&
 !     &                   k_jg)
 !     CALL addVar(TimeVar('rho',&
 !     &                   'density cells',&
 !     &                   'kg/m**3', 6, 128,&
 !     &                   vlistid(k_jg),&
 !     &                   gridcellid(k_jg),&
 !     &                   zaxisIDdepth_m(k_jg)),&
 !     &           k_jg)



      CALL nf(nf_close(ncid))

      ! tracer fields
      !------------------------------------------------------------------
      ! set tracer names
      !------------------------------------------------------------------
      oce_tracer_names(1)     = 'T'
      oce_tracer_names(2)     = 'S'
      oce_tracer_longnames(1) = 'potential temperature'
      oce_tracer_longnames(2) = 'salinity'
      oce_tracer_units(1)     = 'deg C'
      oce_tracer_units(2)     = 'psu'
      oce_tracer_codes(1)     = 200
      oce_tracer_codes(2)     = 201
      IF (no_tracer > oce_max_tracer) THEN
        CALL message(TRIM(routine), 'no_tracer is larger than oce_max_tracer -> limitted to 2')
        oce_trace_counter = oce_max_tracer
      ELSE
        oce_trace_counter = no_tracer
      END IF
      DO itracer = 1, oce_trace_counter
        msg='Create tracer: '//TRIM(oce_tracer_names(itracer))
        CALL message(TRIM(routine), TRIM(msg))
        CALL addVar(TimeVar(TRIM(oce_tracer_names(itracer)),TRIM(oce_tracer_longnames(itracer)),&
        &                   oce_tracer_units(itracer),oce_tracer_codes(itracer),128,        &
        &                   vlistID(k_jg), gridCellID(k_jg),zaxisIDdepth_m(k_jg)),  &
        &           k_jg)
      END DO

    END IF  ! ocean

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
      ELSEIF(gridid == lonlat_gridID(k_jg)) THEN
        outvar_desc(ivar, k_jg)%type = GATHER_LONLAT
      ELSE
        CALL finish(routine,'got illegal gridid')
      ENDIF

      zaxisid = vlistInqVarZaxis(vlistID(k_jg), varids(ivar, k_jg))
      outvar_desc(ivar, k_jg)%nlev = zaxisInqSize(zaxisid)

      CALL vlistInqVarName(vlistID(k_jg), varids(ivar, k_jg), outvar_desc(ivar, k_jg)%name)

    ENDDO

  END SUBROUTINE setup_vlist


  !-------------------------------------------------------------------------------------------------
  SUBROUTINE destruct_vlist(k_jg)

    INTEGER, INTENT(in) :: k_jg
    INTEGER             :: ist ! error status
    !=========================================================================
    ! In parallel mode only 1 PE is writing the output

    IF(.NOT. my_process_is_stdio()) RETURN

    !=========================================================================
    ! cleanup memory
    !
    CALL vlistDestroy(vlistID(k_jg))
    !    CALL taxisDestroy(taxisID(k_jg))
    CALL zaxisDestroy(zaxisID_surface(k_jg))
    CALL gridDestroy(gridVertexID(k_jg))
    CALL gridDestroy(gridEdgeID(k_jg))
    CALL gridDestroy(gridCellID(k_jg))
    IF (iequations/=ihs_ocean) THEN ! atm 
      CALL zaxisDestroy(zaxisID_hybrid(k_jg))
      CALL zaxisDestroy(zaxisID_hybrid_half(k_jg))
    ELSE
      CALL zaxisDestroy(zaxisIDdepth_m(k_jg))
      CALL zaxisDestroy(zaxisID_halfdepth(k_jg))
    END IF

    LONLAT : IF (lonlat_intp_config(k_jg)%l_enabled) THEN
      CALL gridDestroy(lonlat_gridID(k_jg))
      DEALLOCATE(var_lonlat(k_jg)%var_3d, var_lonlat(k_jg)%var_3d_2, &
        &        STAT=ist)
      IF (ist /= SUCCESS) THEN
        CALL finish ('mo_io_vlist:destruct_vlist', 'deallocation failed')
      ENDIF
    END IF LONLAT

    num_output_vars(k_jg) = 0
    !
    !=========================================================================

  END SUBROUTINE destruct_vlist


  !-------------------------------------------------------------------------------------------------
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


  !-------------------------------------------------------------------------
  !> Write vlist variable.
  !! This is a wrapper in order to keep streamID and varids private.
  !! Called by asynchronous IO routines.
  !!
  !! Important note: 
  !! Lon-lat interpolation of variables is not yet functional in ICON!
  !! In module "mo_atmo_model", the vertical coordinate table "vct"
  !! is not initialized by pure IO proc when called in this
  !! order. 

  SUBROUTINE vlist_write_var(ivar, k_jg, var, nlev, typ)
    INTEGER,       INTENT(in)  :: ivar, k_jg
    REAL(wp),      INTENT(in)  :: var(*)
    INTEGER,       INTENT(in)  :: nlev ! number of levels
    INTEGER,       INTENT(in)  :: typ  ! variable type

    ! local variables
    CHARACTER(*), PARAMETER :: routine = TRIM("mo_io_vlist:vlist_write_var")
    LOGICAL            :: l_interpolate_lonlat
    INTEGER            :: ierrstat, dim1, dim3, n_tot
    TYPE(t_collected_var_ptr) :: var_3d

    CALL streamWriteVar(streamID(k_jg), varids(ivar, k_jg), var, 0)

    ! Flag. True, if this variable is due to interpolation onto
    ! lon-lat grid
    l_interpolate_lonlat = (lonlat_intp_config(k_jg)%l_enabled) .AND. &
      &                    (lonlat_varids(ivar,k_jg,1) /= NO_LONLAT_INTERP)
    l_interpolate_lonlat = .FALSE. ! see comment above

    ! perform interpolation onto lon-lat grid
    IF ((l_interpolate_lonlat) .AND. (my_process_is_stdio())) THEN
      
      ! reshape 1D input variable to original 2D/3D array
      dim1 = nproma
      IF(typ == GATHER_C) THEN
        dim3 = (p_patch(k_jg)%n_patch_cells_g-1)/nproma+1
      ELSE IF(typ == GATHER_E) THEN
        dim3 = (p_patch(k_jg)%n_patch_edges_g-1)/nproma+1
      ELSE IF(typ == GATHER_V) THEN
        dim3 = (p_patch(k_jg)%n_patch_verts_g-1)/nproma+1
      ELSE
        CALL finish(routine,'Illegal type parameter')
      END IF

      ALLOCATE(var_3d%ptr(dim1, nlev, dim3), STAT=ierrstat)
      IF (ierrstat /= SUCCESS) THEN
        CALL finish (routine, 'allocation failed')
      ENDIF      
      n_tot = dim1*nlev*dim3
      var_3d%ptr = RESHAPE(var(1:n_tot), (/ dim1, nlev, dim3 /))

      CALL postprocess_lonlat(ivar, typ, var_3d, k_jg, nlev)
      
      DEALLOCATE(var_3d%ptr, STAT=ierrstat)
      IF (ierrstat /= SUCCESS) THEN
        CALL finish (routine, 'deallocation failed')
      ENDIF      
      
    END IF ! (l_interpolate_lonlat)

  END SUBROUTINE vlist_write_var


  !-------------------------------------------------------------------------
  !> Set date and time in taxis
  !! This is a wrapper to keep variables private

  SUBROUTINE vlist_set_date_time(k_jg, idate, itime)
    INTEGER, INTENT(in) :: k_jg, idate, itime

    CALL taxisDefVdate(taxisID(k_jg), idate)   ! YYYYMMDD
    CALL taxisDefVtime(taxisID(k_jg), itime)   ! HHMM

  END SUBROUTINE vlist_set_date_time


  !-------------------------------------------------------------------------
  !> Set timestep
  !! This is a wrapper to keep variables private

  SUBROUTINE vlist_start_step(k_jg, istep)
    INTEGER, INTENT(in) :: k_jg, istep
    INTEGER :: istatus

    istatus = streamDefTimestep(streamID(k_jg), istep)

  END SUBROUTINE vlist_start_step


  !-------------------------------------------------------------------------------------------------
  !>
  !! Get pointer to an output variable given its name
  !! If reset is set on exit, the variable has to be set to 0 after I/O
  !! If delete is set after I/O, the pointer has to be deallocated
  !! (because it is not pointing to another variable but has been allocated here).

  SUBROUTINE get_outvar_ptr_ha(varname, jg, ptr2, ptr3, reset, delete)

    CHARACTER(LEN=*), INTENT(IN) :: varname
    INTEGER, INTENT(IN) :: jg

    REAL(wp), POINTER :: ptr2(:,:)
    REAL(wp), POINTER :: ptr3(:,:,:)
    LOGICAL, INTENT(OUT) :: reset, delete

    INTEGER :: jt
    INTEGER :: nlevp1
    LOGICAL :: not_found
    CHARACTER(LEN=1) :: ctracer, ctile
    CHARACTER(len=MAX_CHAR_LENGTH) :: & !< list of tracers to initialize
      &  ctracer_list

    TYPE(t_hydro_atm_prog), POINTER :: p_prog
    TYPE(t_hydro_atm_diag), POINTER :: p_diag

    p_prog => p_hydro_state(jg)%prog_out
    p_diag => p_hydro_state(jg)%diag_out

    ptr2 => NULL()
    ptr3 => NULL()
    reset  = .FALSE.
    delete = .FALSE.
    not_found = .FALSE.

    ! get ctracer_list
    ctracer_list = advection_config(jg)%ctracer_list

    nlevp1 = num_levp1(jg)

    SELECT CASE(varname)
      CASE ('SLM');             ptr2 => prm_field(jg)%lsmask
      CASE ('SKT');             ptr2 => prm_field(jg)%tsfc(:,:)
      CASE ('SKT_TILE');        ptr3 => prm_field(jg)%tsfc_tile(:,:,:)
      CASE ('PS');              ptr2 => p_prog%pres_sfc
      CASE ('T');               ptr3 => p_prog%temp
      CASE ('normal_velocity'); ptr3 => p_prog%vn
      CASE ('U');               ptr3 => p_diag%u
      CASE ('V');               ptr3 => p_diag%v
      CASE ('OMEGA');           ptr3 => p_diag%wpres_mc
      CASE ('P');               ptr3 => p_diag%pres_mc
      CASE ('ZF3');             ptr3 => dup3(p_diag%geo_mc/grav); delete = .TRUE.
      CASE ('ozone');           ptr3 => prm_field(jg)%o3
      CASE ('PHIS');            ptr2 => p_diag%geo_ic(:,nlevp1,:)
      CASE ('cosmu0');          ptr2 => prm_field(jg)%cosmu0
      CASE ('flxdwswtoa');      ptr2 => prm_field(jg)%flxdwswtoa
      CASE ('APRL');            ptr2 => prm_field(jg)%aprl(:,:);   reset = .TRUE.
      CASE ('APRC');            ptr2 => prm_field(jg)%aprc(:,:);   reset = .TRUE.
      CASE ('APRS');            ptr2 => prm_field(jg)%aprs(:,:);   reset = .TRUE.
      CASE ('RSFL');            ptr2 => prm_field(jg)%rsfl(:,:)
      CASE ('RSFC');            ptr2 => prm_field(jg)%rsfc(:,:)
      CASE ('SSFL');            ptr2 => prm_field(jg)%ssfl(:,:)
      CASE ('SSFC');            ptr2 => prm_field(jg)%ssfc(:,:)
      CASE ('ACLC');            ptr3 => prm_field(jg)%aclc
      CASE ('ACLCAC');          ptr3 => prm_field(jg)%aclcac;      reset = .TRUE.
      CASE ('ACLCOV');          ptr2 => prm_field(jg)%aclcov(:,:); reset = .TRUE.
      CASE ('qvi');             ptr2 => prm_field(jg)%qvi (:,:);   reset = .TRUE.
      CASE ('xlvi');            ptr2 => prm_field(jg)%xlvi(:,:);   reset = .TRUE.
      CASE ('xivi');            ptr2 => prm_field(jg)%xivi(:,:);   reset = .TRUE.

      CASE ('evap_ac');         ptr2 => prm_field(jg)%    evap_ac(:,:); reset = .TRUE.
      CASE ('lhflx_ac');        ptr2 => prm_field(jg)%   lhflx_ac(:,:); reset = .TRUE.
      CASE ('shflx_ac');        ptr2 => prm_field(jg)%   shflx_ac(:,:); reset = .TRUE.
      CASE ('u_stress_ac');     ptr2 => prm_field(jg)%u_stress_ac(:,:); reset = .TRUE.
      CASE ('v_stress_ac');     ptr2 => prm_field(jg)%v_stress_ac(:,:); reset = .TRUE.

      CASE ('OMEGA_PHY');       ptr3 => prm_field(jg)%omega
      CASE ('tend_temp_radsw'); ptr3 => prm_tend(jg)%temp_radsw
      CASE ('tend_temp_radlw'); ptr3 => prm_tend(jg)%temp_radlw
      CASE ('tend_temp_cld');   ptr3 => prm_tend(jg)%temp_cld
      CASE ('tend_temp_cnv');   ptr3 => prm_tend(jg)%temp_cnv
      CASE ('tend_temp_vdf');   ptr3 => prm_tend(jg)%temp_vdf
      CASE ('tend_temp_gwh');   ptr3 => prm_tend(jg)%temp_gwh
      CASE ('tend_u_cnv');      ptr3 => prm_tend(jg)%u_cnv
      CASE ('tend_u_vdf');      ptr3 => prm_tend(jg)%u_vdf
      CASE ('tend_u_gwh');      ptr3 => prm_tend(jg)%u_gwh
      CASE ('tend_v_cnv');      ptr3 => prm_tend(jg)%v_cnv
      CASE ('tend_v_vdf');      ptr3 => prm_tend(jg)%v_vdf
      CASE ('tend_v_gwh');      ptr3 => prm_tend(jg)%v_gwh
      CASE ('VOR');             ptr3 => p_diag%rel_vort
      CASE ('DIV');             ptr3 => p_diag%div
      !
!!$      CASE ('debug_2d_1');      ptr2 => prm_field(jg)%debug_2d_1
!!$      CASE ('debug_2d_2');      ptr2 => prm_field(jg)%debug_2d_2
!!$      CASE ('debug_2d_3');      ptr2 => prm_field(jg)%debug_2d_3
!!$      CASE ('debug_2d_4');      ptr2 => prm_field(jg)%debug_2d_4
!!$      CASE ('debug_2d_5');      ptr2 => prm_field(jg)%debug_2d_5
!!$      CASE ('debug_2d_6');      ptr2 => prm_field(jg)%debug_2d_6
!!$      CASE ('debug_2d_7');      ptr2 => prm_field(jg)%debug_2d_7
!!$      CASE ('debug_2d_8');      ptr2 => prm_field(jg)%debug_2d_8
      !
!!$      CASE ('debug_3d_1');      ptr3 => prm_field(jg)%debug_3d_1
!!$      CASE ('debug_3d_2');      ptr3 => prm_field(jg)%debug_3d_2
!!$      CASE ('debug_3d_3');      ptr3 => prm_field(jg)%debug_3d_3
!!$      CASE ('debug_3d_4');      ptr3 => prm_field(jg)%debug_3d_4
!!$      CASE ('debug_3d_5');      ptr3 => prm_field(jg)%debug_3d_5
!!$      CASE ('debug_3d_6');      ptr3 => prm_field(jg)%debug_3d_6
!!$      CASE ('debug_3d_7');      ptr3 => prm_field(jg)%debug_3d_7
!!$      CASE ('debug_3d_8');      ptr3 => prm_field(jg)%debug_3d_8
      !
      CASE DEFAULT;             not_found = .TRUE.
    END SELECT

    ! If not found in the list above, check for tracers, tracer tendencies, or 
    ! tile-specific quantities

    IF(not_found) THEN

      DO jt = 1, ntracer
        ctracer = ctracer_list(jt:jt)
        IF(varname == 'Q'//ctracer) THEN
          ptr3 => p_prog%tracer(:,:,:,jt)
          RETURN
        ENDIF
      ENDDO

      DO jt = 1, iqcond
        ctracer = ctracer_list(jt:jt)
        IF(varname == 'tend_q'//ctracer//'_cnv') THEN
          ptr3 => prm_tend(jg)%q_cnv(:,:,:,jt)
          RETURN
        ENDIF
        IF(varname == 'tend_q'//ctracer//'_vdf') THEN
          ptr3 => prm_tend(jg)%q_vdf(:,:,:,jt)
          RETURN
        ENDIF
      ENDDO

      DO jt = 1, nsfc_type 
        WRITE(ctile,'(i1)') jt

        IF(varname == 'evap_tile_'//ctile) THEN
          ptr2 => prm_field(jg)%evap_tile(:,:,jt)
          RETURN
        ENDIF

        IF(varname == 'lhflx_tile_'//ctile) THEN
          ptr2 => prm_field(jg)%lhflx_tile(:,:,jt)
          RETURN
        ENDIF

        IF(varname == 'shflx_tile_'//ctile) THEN
          ptr2 => prm_field(jg)%shflx_tile(:,:,jt)
          RETURN
        ENDIF

        IF(varname == 'u_stress_tile_'//ctile) THEN
          ptr2 => prm_field(jg)%u_stress_tile(:,:,jt)
          RETURN
        ENDIF

        IF(varname == 'v_stress_tile_'//ctile) THEN
          ptr2 => prm_field(jg)%v_stress_tile(:,:,jt)
          RETURN
        ENDIF

      ENDDO

      ! If we are here, the varname was definitly not found
      CALL finish('get_outvar_ptr_ha', 'Unkown variable type name: '//varname)
    ENDIF

!  CONTAINS


  END SUBROUTINE get_outvar_ptr_ha


  !-------------------------------------------------------------------------------------------------
  !>
  !! Get pointer to an output variable given its name (for NH state)
  !! If reset is set on exit, the variable has to be set to 0 after I/O
  !! If delete is set after I/O, the pointer has to be deallocated
  !! (because it is not pointing to another variable but has been allocated here).

  SUBROUTINE get_outvar_ptr_nh(varname, jg, z_sim_time,ptr2, ptr3, reset, delete)

    CHARACTER(LEN=*), INTENT(IN) :: varname
    INTEGER, INTENT(IN) :: jg
    REAL(wp), INTENT(in) :: z_sim_time
    LOGICAL, INTENT(OUT) :: reset, delete

    INTEGER :: jt
    LOGICAL :: not_found
    CHARACTER(LEN=1) :: ctracer
    CHARACTER(len=MAX_CHAR_LENGTH) :: & !< list of tracers to initialize
      &  ctracer_list

    CHARACTER(LEN=1) :: anextra
    CHARACTER(LEN=21):: name
    CHARACTER(LEN=3) :: cjt

    REAL(wp), POINTER :: ptr2(:,:)
    REAL(wp), POINTER :: ptr3(:,:,:)

    TYPE(t_nh_prog),    POINTER :: p_prog
    TYPE(t_nh_diag),    POINTER :: p_diag
    TYPE(t_nh_diag_pz), POINTER :: p_diag_p
    TYPE(t_nh_diag_pz), POINTER :: p_diag_z

 !   TYPE(t_nwp_phy_diag) :: prm_diag
 !   TYPE(t_nwp_phy_tend) :: prm_nwp_tend

    TYPE(t_lnd_prog), POINTER :: p_prog_lnd
    TYPE(t_lnd_diag), POINTER :: p_diag_lnd

    p_prog => p_nh_state(jg)%prog(nnow(jg))
    p_diag => p_nh_state(jg)%diag

    IF (lwrite_pzlev) THEN
      p_diag_p => p_nh_state(jg)%diag_p
      p_diag_z => p_nh_state(jg)%diag_z
    ENDIF

    IF (iforcing==inwp) THEN
     !p_prm_diag     => prm_diag(jg)
     !p_prm_nwp_tend => prm_nwp_tend(jg)
      p_prog_lnd     => p_lnd_state(jg)%prog_lnd(nnow(jg))
      p_diag_lnd     => p_lnd_state(jg)%diag_lnd
    ENDIF

    ptr2 => NULL()
    ptr3 => NULL()
    reset  = .FALSE.
    delete = .FALSE.
    not_found = .FALSE.

    SELECT CASE(varname)
      CASE ('PS');              ptr2 => p_diag%pres_sfc
      CASE ('PS_s6avg');        ptr2 => p_diag%pres_sfc_s6avg
      CASE ('T');               ptr3 => p_diag%temp
      CASE ('T_2m' );           ptr2 => prm_diag(jg)%t_2m
      CASE ('T_2m_s6avg' );     ptr2 => prm_diag(jg)%t_2m_s6avg
      CASE ('QV_2m' );           ptr2 => prm_diag(jg)%qv_2m
      CASE ('QV_2m_s6avg' );     ptr2 => prm_diag(jg)%qv_2m_s6avg
      CASE ('normal_velocity'); ptr3 => p_prog%vn
      CASE ('U');               ptr3 => p_diag%u
      CASE ('V');               ptr3 => p_diag%v
      CASE ('W');               ptr3 => p_prog%w
      CASE ('U_10m');           ptr2 => prm_diag(jg)%u_10m
      CASE ('V_10m');           ptr2 => prm_diag(jg)%v_10m
      CASE ('U_10m_s6avg');     ptr2 => prm_diag(jg)%u_10m_s6avg
      CASE ('V_10m_s6avg');     ptr2 => prm_diag(jg)%v_10m_s6avg
      CASE ('P');               ptr3 => p_diag%pres
      CASE ('QV');              ptr3 => prm_diag(jg)%tot_cld(:,:,:,iqv)
      CASE ('QC');              ptr3 => prm_diag(jg)%tot_cld(:,:,:,iqc)
      CASE ('QI');              ptr3 => prm_diag(jg)%tot_cld(:,:,:,iqi)
      CASE ('CC');              ptr3 => prm_diag(jg)%tot_cld(:,:,:,icc)
      CASE ('TQV');             ptr2 => prm_diag(jg)%tot_cld_vi(:,:,iqv)
      CASE ('TQC');             ptr2 => prm_diag(jg)%tot_cld_vi(:,:,iqc)
      CASE ('TQI');             ptr2 => prm_diag(jg)%tot_cld_vi(:,:,iqi)
      CASE ('TCC');             ptr2 => prm_diag(jg)%tot_cld_vi(:,:,icc)
      CASE ('TQV_avg');           ptr2 => prm_diag(jg)%tot_cld_vi_avg(:,:,iqv)
      CASE ('TQC_avg');           ptr2 => prm_diag(jg)%tot_cld_vi_avg(:,:,iqc)
      CASE ('TQI_avg');           ptr2 => prm_diag(jg)%tot_cld_vi_avg(:,:,iqi)
      CASE ('TCC_avg');           ptr2 => prm_diag(jg)%tot_cld_vi_avg(:,:,icc)
      CASE ('ZF3');             ptr3 => p_nh_state(jg)%metrics%z_mc
      CASE ('ZH3');             ptr3 => p_nh_state(jg)%metrics%z_ifc
      CASE ('PRR_GSP');         ptr2 => prm_diag(jg)%tracer_rate(:,:,1)
      CASE ('PRS_GSP');         ptr2 => prm_diag(jg)%tracer_rate(:,:,2)
      CASE ('RAIN_GSP');        ptr2 => prm_diag(jg)%rain_gsp(:,:)
      CASE ('SNOW_GSP');        ptr2 => prm_diag(jg)%snow_gsp(:,:)
      CASE ('RAIN_CON');        ptr2 => prm_diag(jg)%rain_con(:,:)
      CASE ('SNOW_CON');        ptr2 => prm_diag(jg)%snow_con(:,:)
      CASE ('TOT_PREC');        ptr2 => prm_diag(jg)%tot_prec(:,:)
      CASE ('TOT_PREC_RATE_avg'); ptr2 => prm_diag(jg)%tot_prec_rate_avg(:,:)
      CASE ('CON_PREC_RATE_avg'); ptr2 => prm_diag(jg)%con_prec_rate_avg(:,:)
      CASE ('GSP_PREC_RATE_avg'); ptr2 => prm_diag(jg)%gsp_prec_rate_avg(:,:)
      CASE ('cosmu0');          ptr2 => prm_diag(jg)%cosmu0(:,:)
      CASE ('flxdwswtoa');      ptr2 => prm_diag(jg)%flxdwswtoa(:,:)
      CASE ('swflxsfc');        ptr2 => prm_diag(jg)%swflxsfc(:,:)
      CASE ('lwflxsfc');        ptr2 => prm_diag(jg)%lwflxsfc(:,:)
      CASE ('swflxtoa');        ptr2 => prm_diag(jg)%swflxtoa(:,:)
      CASE ('lwflxtoa');        ptr2 => prm_diag(jg)%lwflxall(:,1,:)
      CASE ('swflxsfc_avg');    ptr2 => prm_diag(jg)%swflxsfc_avg(:,:)
      CASE ('lwflxsfc_avg');    ptr2 => prm_diag(jg)%lwflxsfc_avg(:,:)       
      CASE ('swflxtoa_avg');    ptr2 => prm_diag(jg)%swflxtoa_avg(:,:)
      CASE ('lwflxtoa_avg');    ptr2 => prm_diag(jg)%lwflxtoa_avg(:,:)
      CASE ('T_G');             ptr2 => p_prog_lnd%t_g
      CASE ('QV_S');            ptr2 => p_diag_lnd%qv_s
      CASE ('SHFL_S_avg');      ptr2 => prm_diag(jg)%shfl_s_avg 
      CASE ('LHFL_S_avg');      ptr2 => prm_diag(jg)%lhfl_s_avg
      CASE ('EVAP_RATE_avg');   ptr2 => prm_diag(jg)%qhfl_s_avg     
      CASE ('VOR');             ptr3 => p_diag%omega_z
      CASE ('DIV');             ptr3 => p_diag%div
      CASE ('THETA_V');         ptr3 => p_prog%theta_v
      CASE ('EXNER');           ptr3 => p_prog%exner
      CASE ('RHO');             ptr3 => p_prog%rho
      CASE ('TKE')
        IF (atm_phy_nwp_config(jg)%inwp_turb.EQ.1) THEN
                                ptr3 => dup3(0.5_wp*SQRT(p_prog%tke(:,:,:))); delete = .TRUE.
           ELSE
                                ptr3 => p_prog%tke
           ENDIF
      CASE ('Z0')
        IF (atm_phy_nwp_config(jg)%inwp_turb.EQ.1) THEN
                                ptr2 => dup2(prm_diag(jg)%gz0(:,:)/grav); delete = .TRUE.
                              ELSE
                                ptr2 => prm_diag(jg)%z0m(:,:)
        END IF
      CASE ('tend_temp_radsw'); ptr3 => prm_nwp_tend(jg)%ddt_temp_radsw(:,:,:)
      CASE ('tend_temp_radlw'); ptr3 => prm_nwp_tend(jg)%ddt_temp_radlw(:,:,:)
      CASE ('tend_temp_pscl');  ptr3 => prm_nwp_tend(jg)%ddt_temp_pscl (:,:,:)
      CASE ('tend_temp_conv');  ptr3 => prm_nwp_tend(jg)%ddt_temp_pconv(:,:,:)
      CASE ('tend_temp_turb');  ptr3 => prm_nwp_tend(jg)%ddt_temp_turb (:,:,:)
      CASE ('tend_temp_sso');   ptr3 => prm_nwp_tend(jg)%ddt_temp_sso  (:,:,:)
      CASE ('tend_temp_dyn');   ptr3 => p_diag%ddt_temp_dyn(:,:,:)
      CASE ('tend_u_conv');     ptr3 => prm_nwp_tend(jg)%ddt_u_pconv   (:,:,:)
      CASE ('tend_u_turb');     ptr3 => prm_nwp_tend(jg)%ddt_u_turb    (:,:,:)
      CASE ('tend_u_sso');      ptr3 => prm_nwp_tend(jg)%ddt_u_turb    (:,:,:)
      CASE ('tend_v_conv');     ptr3 => prm_nwp_tend(jg)%ddt_v_pconv   (:,:,:)
      CASE ('tend_v_turb');     ptr3 => prm_nwp_tend(jg)%ddt_v_turb    (:,:,:)
      CASE ('tend_v_sso');      ptr3 => prm_nwp_tend(jg)%ddt_v_turb    (:,:,:)
      CASE ('U_Z');             ptr3 => p_diag_z%u     (:,:,:)
      CASE ('V_Z');             ptr3 => p_diag_z%v     (:,:,:)
      CASE ('T_Z');             ptr3 => p_diag_z%temp  (:,:,:)
      CASE ('P_Z');             ptr3 => p_diag_z%pres  (:,:,:)
      CASE ('QV_Z');            ptr3 => p_diag_z%tot_cld(:,:,:,iqv)
      CASE ('QC_Z');            ptr3 => p_diag_z%tot_cld(:,:,:,iqc)
      CASE ('QI_Z');            ptr3 => p_diag_z%tot_cld(:,:,:,iqi)
      CASE ('CC_Z');            ptr3 => p_diag_z%tot_cld(:,:,:,icc)
      CASE ('U_P');             ptr3 => p_diag_p%u     (:,:,:)
      CASE ('V_P');             ptr3 => p_diag_p%v     (:,:,:)
      CASE ('T_P');             ptr3 => p_diag_p%temp  (:,:,:)
      CASE ('Z');               ptr3 => dup3(grav*p_diag_p%geopot(:,:,:))
      CASE ('QV_P');            ptr3 => p_diag_p%tot_cld(:,:,:,iqv)
      CASE ('QC_P');            ptr3 => p_diag_p%tot_cld(:,:,:,iqc)
      CASE ('QI_P');            ptr3 => p_diag_p%tot_cld(:,:,:,iqi)
      CASE ('CC_P');            ptr3 => p_diag_p%tot_cld(:,:,:,icc)
      CASE DEFAULT;             not_found = .TRUE.
    END SELECT

    ! If not found in the list above, check for tracers


    ! get ctracer_list
    ctracer_list = advection_config(jg)%ctracer_list

    IF(not_found) THEN
      DO jt = 1, ntracer ! all tracer
        ctracer = ctracer_list(jt:jt)
        IF(varname == 'Q'//ctracer) THEN
          ptr3 => p_prog%tracer(:,:,:,jt)
          RETURN
        ENDIF
        IF(varname == 'Q'//ctracer//'_P') THEN
          ptr3 => p_diag_p%tracer(:,:,:,jt)
          RETURN
        ENDIF
        IF(varname == 'Q'//ctracer//'_Z') THEN
          ptr3 => p_diag_z%tracer(:,:,:,jt)
          RETURN
        ENDIF
      ENDDO

      DO jt = 1, 3
        ctracer = ctracer_list(jt:jt)
        IF(varname == 'TQ'//ctracer) THEN
          ptr2 => p_diag%tracer_vi(:,:,jt)
          RETURN
        ENDIF
        IF(varname == 'TQ'//ctracer//'_avg') THEN
          ptr2 => p_diag%tracer_vi_avg(:,:,jt)
          RETURN
        ENDIF
      ENDDO
      DO jt = 1, iqcond ! only hydromet-tendencies
        ctracer = ctracer_list(jt:jt)
        IF(varname == 'tend_q'//ctracer//'_conv') THEN
          ptr3 => prm_nwp_tend(jg)%ddt_tracer_pconv(:,:,:,jt)
          RETURN
        ENDIF
        IF(varname == 'tend_q'//ctracer//'_turb') THEN
          ptr3 => prm_nwp_tend(jg)%ddt_tracer_turb(:,:,:,jt)
          RETURN
        ENDIF
!        IF(varname == 'tend_q'//ctracer//'_pscl') THEN
!          ptr3 => prm_nwp_tend(jg)%ddt_tracer_gscp(:,:,:,jt)
!          RETURN
!        ENDIF
      ENDDO

      IF (irad_o3 == 4 .OR. irad_o3 == 6) THEN
        IF(varname == 'O3') THEN
          ptr3 => p_prog%tracer(:,:,:,io3)
          RETURN
        ENDIF
      END IF

      DO jt = 1, inextra_2d ! 2d debug variables
        WRITE(anextra,'(I1)') jt
        IF(varname == 'extra_2d_'//anextra) THEN
          ptr2 => p_diag%extra_2d(:,:,jt)
          RETURN
        ENDIF
      ENDDO

      DO jt = 1, inextra_3d ! 3d debugging
        WRITE(anextra,'(I1)') jt
        IF(varname == 'extra_3d_'//anextra) THEN
          ptr3 => p_diag%extra_3d(:,:,:,jt)
          RETURN
        ENDIF
      ENDDO

      ! check for land variables
      !
      IF ( atm_phy_nwp_config(jg)%inwp_surface == 1 ) THEN  ! TERRA

      DO jt = 1, nsfc_subs
        WRITE(cjt, '(i2)') jt
        WRITE(name,'(A,A)') "T_GT_tile_", TRIM(ADJUSTL(cjt))
        IF(varname == TRIM(name)) THEN
          ptr2 => p_prog_lnd%t_gt(:,:,jt)
          RETURN
        ENDIF

        WRITE(name,'(A,A)') "T_SNOW_tile_", TRIM(ADJUSTL(cjt))
        IF(varname == TRIM(name)) THEN
          ptr2 => p_prog_lnd%t_snow(:,:,jt)
          RETURN
        ENDIF

        WRITE(name,'(A,A)') "T_SNOW_MULT_tile_", TRIM(ADJUSTL(cjt))
        IF(varname == TRIM(name)) THEN
          ptr3 => p_prog_lnd%t_snow_mult(:,:,:,jt)
          RETURN
        ENDIF

        WRITE(name,'(A,A)') "T_S_tile_", TRIM(ADJUSTL(cjt))
        IF(varname == TRIM(name)) THEN
          ptr2 => p_prog_lnd%t_s(:,:,jt)
          RETURN
        ENDIF

        WRITE(name,'(A,A)') "W_SNOW_tile_", TRIM(ADJUSTL(cjt))
        IF(varname == TRIM(name)) THEN
          ptr2 => p_prog_lnd%w_snow(:,:,jt)
          RETURN
        ENDIF

        WRITE(name,'(A,A)') "RHO_SNOW_tile_", TRIM(ADJUSTL(cjt))
        IF(varname == TRIM(name)) THEN
          ptr2 => p_prog_lnd%rho_snow(:,:,jt)
          RETURN
        ENDIF

        WRITE(name,'(A,A)') "RHO_SNOW_MULT_tile_", TRIM(ADJUSTL(cjt))
        IF(varname == TRIM(name)) THEN
          ptr3 => p_prog_lnd%rho_snow_mult(:,:,:,jt)
          RETURN
        ENDIF

        WRITE(name,'(A,A)') "W_I_tile_", TRIM(ADJUSTL(cjt))
        IF(varname == TRIM(name)) THEN
          ptr2 => p_prog_lnd%w_i(:,:,jt)
          RETURN
        ENDIF

        WRITE(name,'(A,A)') "T_SO_tile_", TRIM(ADJUSTL(cjt))
        IF(varname == TRIM(name)) THEN
          ptr3 => p_prog_lnd%t_so(:,:,:,jt)
          RETURN
        ENDIF

        WRITE(name,'(A,A)') "W_SO_tile_", TRIM(ADJUSTL(cjt))
        IF(varname == TRIM(name)) THEN
          ptr3 => p_prog_lnd%w_so(:,:,:,jt)
          RETURN
        ENDIF

        WRITE(name,'(A,A)') "W_SO_ICE_tile_", TRIM(ADJUSTL(cjt))
        IF(varname == TRIM(name)) THEN
          ptr3 => p_prog_lnd%w_so_ice(:,:,:,jt)
          RETURN
        ENDIF

        WRITE(name,'(A,A)') "WLIQ_SNOW_tile_", TRIM(ADJUSTL(cjt))
        IF(varname == TRIM(name)) THEN
          ptr3 => p_prog_lnd%wliq_snow(:,:,:,jt)
          RETURN
        ENDIF

        WRITE(name,'(A,A)') "WTOT_SNOW_tile_", TRIM(ADJUSTL(cjt))
        IF(varname == TRIM(name)) THEN
          ptr3 => p_prog_lnd%wtot_snow(:,:,:,jt)
          RETURN
        ENDIF

        WRITE(name,'(A,A)') "H_SNOW_tile_", TRIM(ADJUSTL(cjt))
        IF(varname == TRIM(name)) THEN
          ptr2 => p_diag_lnd%h_snow(:,:,jt)
          RETURN
        ENDIF

        WRITE(name,'(A,A)') "DZH_SNOW_tile_", TRIM(ADJUSTL(cjt))
        IF(varname == TRIM(name)) THEN
          ptr3 => p_prog_lnd%dzh_snow(:,:,:,jt)
          RETURN
        ENDIF
      ENDDO

      ENDIF

      ! If we are here, the varname was definitly not found
      CALL finish('get_outvar_ptr_nh', 'Unknown variable type name: '//varname)
    ENDIF

  END SUBROUTINE get_outvar_ptr_nh


  !-------------------------------------------------------------------------------------------------
  SUBROUTINE get_outvar_ptr_oce(varname,jg,ptr2d,ptr3d,reset,delete)

    CHARACTER(LEN=*), INTENT(IN) :: varname
    INTEGER, INTENT(IN) :: jg
    REAL(wp), POINTER :: ptr2d(:,:)
    REAL(wp), POINTER :: ptr3d(:,:,:)


    LOGICAL, INTENT(OUT) :: reset, delete

    LOGICAL :: not_found

    TYPE(t_hydro_ocean_prog), POINTER  :: p_prog
    TYPE(t_hydro_ocean_diag), POINTER  :: p_diag
    TYPE(t_hydro_ocean_aux),  POINTER  :: p_aux
    TYPE(t_sfc_flx),          POINTER  :: forcing

    ! pointer to components of state variable
    p_prog  => v_ocean_state(jg)%p_prog(nold(jg))
    p_diag  => v_ocean_state(jg)%p_diag
    p_aux   => v_ocean_state(jg)%p_aux
    forcing => v_sfc_flx

    ptr2d     => NULL()
    ptr3d     => NULL()
    reset     = .FALSE.
    delete    = .FALSE.
    not_found = .FALSE.

    SELECT CASE(varname)
      CASE ('wet_c');        ptr3d => v_base%wet_c
      CASE ('wet_e');        ptr3d => v_base%wet_e
      CASE ('ELEV');         ptr2d => p_prog%h
      CASE ('forc-u');       ptr2d => forcing%forc_wind_u
      CASE ('forc-v');       ptr2d => forcing%forc_wind_v
      CASE('horz-adv');      ptr3d => p_diag%veloc_adv_horz
      CASE('Ekin-grad');     ptr3d => p_diag%grad
      CASE ('flux-u');       ptr2d => p_aux%bc_top_u
      CASE ('flux-v');       ptr2d => p_aux%bc_top_v
      CASE ('flux-VN');      ptr2d => p_aux%bc_top_vn
      CASE ('VN');           ptr3d => p_prog%vn
      CASE ('T');            ptr3d => p_prog%tracer(:,:,:,1)
      CASE ('S');            ptr3d => p_prog%tracer(:,:,:,2)
      CASE ('VORT');         ptr3d => p_diag%vort
      CASE ('u-veloc');      ptr3d => p_diag%u
      CASE ('v-veloc');      ptr3d => p_diag%v
      CASE ('W');            ptr3d => p_diag%w
      !CASE('press-grad');    ptr3d => p_diag%press_grad
      !CASE ('rho');          ptr3d => p_diag%rho

      CASE DEFAULT;    not_found = .TRUE.
    END SELECT
    
    IF(not_found) THEN
      CALL finish('get_outvar_ptr_oce', 'Unknown variable type name: '//varname)
    END IF
    !   DO jt = 1, ntracer ! all tracer
    !     ctracer = ctracer_list(jt:jt)
    !     IF(varname == 'Q'//ctracer) THEN
    !       ptr3 => p_prog%tracer(:,:,:,jt)
    !       RETURN
    !     ENDIF
    !   ENDDO
    !   DO jt = 1, iqcond ! only hydromet-tendencies
    !     ctracer = ctracer_list(jt:jt)
    !     IF(varname == 'tend_q'//ctracer//'_conv') THEN
    !       ptr3 => prm_nwp_tend(jg)%ddt_tracer_pconv(:,:,:,jt)
    !       RETURN
    !     ENDIF
    !     IF(varname == 'tend_q'//ctracer//'_turb') THEN
    !       ptr3 => prm_nwp_tend(jg)%ddt_tracer_turb(:,:,:,jt)
    !       RETURN
    !     ENDIF
!   !      IF(varname == 'tend_q'//ctracer//'_pscl') THEN
!   !        ptr3 => prm_nwp_tend(jg)%ddt_tracer_gscp(:,:,:,jt)
!   !        RETURN
!   !      ENDIF
    !   ENDDO

    !   DO jt = 1, inextra_2d ! 2d debug variables
    !     WRITE(anextra,'(I1)') jt
    !     IF(varname == 'extra_2d_'//anextra) THEN
    !       ptr2 => p_diag%extra_2d(:,:,jt)
    !       RETURN
    !     ENDIF
    !   ENDDO

    !   DO jt = 1, inextra_3d ! 3d debugging
    !     WRITE(anextra,'(I1)') jt
    !     IF(varname == 'extra_3d_'//anextra) THEN
    !       ptr3 => p_diag%extra_3d(:,:,:,jt)
    !       RETURN
    !     ENDIF
    !   ENDDO
    !  ENDIF
    ! If we are here, the varname was definitly not found

  END SUBROUTINE get_outvar_ptr_oce


  !-------------------------------------------------------------------------------------------------
  SUBROUTINE write_vlist (datetime, z_sim_time)

    !=========================================================================
    !> Write output directly: PE 0 gathers and writes, others send
    !=========================================================================

    TYPE(t_datetime),            INTENT(in) :: datetime
    REAL(wp), OPTIONAL,          INTENT(in) :: z_sim_time(n_dom)
    INTEGER :: idate, itime
    INTEGER :: istatus, ierrstat
    INTEGER :: jg, ivar, n_tot

    REAL(wp), POINTER :: ptr2(:,:)
    REAL(wp), POINTER :: ptr3(:,:,:)
    LOGICAL :: reset, delete

    REAL(wp), ALLOCATABLE :: streamvar1(:), streamvar2(:,:)
    ! complete 3d field (nproma,nlev,nblks)
    TYPE(t_collected_var_ptr) :: collected_var_3d
    REAL(wp) :: p_sim_time
    INTEGER  :: nlev
    LOGICAL l_interpolate_lonlat
    
    idate   = cdiEncodeDate(datetime%year, datetime%month, datetime%day)
    itime   = cdiEncodeTime(datetime%hour, datetime%minute, NINT(datetime%second))

    ! Make streamvar1/streamvar2 defined everywhere

        IF(.NOT. my_process_is_stdio()) ALLOCATE(streamvar1(1), streamvar2(1,1))

    DO jg = 1, n_dom

      IF (PRESENT(z_sim_time)) THEN
        p_sim_time = z_sim_time(jg)
      ELSE
        p_sim_time = 1.0_wp
      ENDIF
      
      IF(my_process_is_stdio()) THEN
        CALL taxisDefVdate(taxisID(jg), idate)   ! YYYYMMDD
        CALL taxisDefVtime(taxisID(jg), itime)   ! HHMM

        istatus = streamDefTimestep(streamID(jg), iostep)
      ENDIF

      DO ivar = 1, num_output_vars(jg)

        ! lon-lat interpolated variables are processed together with
        ! their source variables:
        IF (outvar_desc(ivar, jg)%type == GATHER_LONLAT) CYCLE

        ! Flag. True, if this variable is due to interpolation onto
        ! lon-lat grid
        l_interpolate_lonlat = (lonlat_intp_config(jg)%l_enabled) .AND. &
          &                    (lonlat_varids(ivar,jg,1) /= NO_LONLAT_INTERP)

        ! Get a pointer to the variable
        SELECT CASE (iequations)
          CASE (ishallow_water, ihs_atm_temp, ihs_atm_theta)
            CALL get_outvar_ptr_ha(outvar_desc(ivar,jg)%name, jg, ptr2, ptr3, reset, delete)
          CASE (inh_atmosphere)
            CALL get_outvar_ptr_nh &
              & (outvar_desc(ivar,jg)%name, jg, p_sim_time, ptr2, ptr3, reset, delete)
          CASE (ihs_ocean)
            CALL get_outvar_ptr_oce(outvar_desc(ivar,jg)%name, jg, ptr2, ptr3,reset, delete)
          CASE DEFAULT
            CALL finish('write_vlist','Unsupported value of iequations')
        END SELECT

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
        IF_2D_3D : IF(ASSOCIATED(ptr2)) THEN

          IF(my_process_is_stdio()) ALLOCATE(streamvar1(n_tot))

          CALL gather_array1( outvar_desc(ivar, jg)%type, p_patch(jg), ptr2, &
            &                 streamvar1,outvar_desc(ivar,jg)%name, collected_var_3d )

          IF(my_process_is_stdio()) THEN
            CALL streamWriteVar(streamID(jg), varids(ivar,jg), streamvar1, 0)
            DEALLOCATE(streamvar1)
          ENDIF
          IF(reset) ptr2 = 0._wp
          IF(delete) DEALLOCATE(ptr2)

          nlev = 1

        ELSE
          IF(my_process_is_stdio()) ALLOCATE(streamvar2(n_tot, klev))

          CALL gather_array2( outvar_desc(ivar, jg)%type, p_patch(jg), ptr3,  &
            &                 streamvar2,outvar_desc(ivar,jg)%name, collected_var_3d )

          IF(my_process_is_stdio()) THEN
            CALL streamWriteVar(streamID(jg), varids(ivar,jg), streamvar2, 0)
            DEALLOCATE(streamvar2)
          ENDIF

          nlev = SIZE(ptr3(:,:,:), 2)
          IF(reset) ptr3 = 0._wp
          IF(delete) DEALLOCATE(ptr3)

        END IF IF_2D_3D

        ! perform interpolation onto lon-lat grid
        IF ((l_interpolate_lonlat) .AND. (my_process_is_stdio())) THEN
          CALL postprocess_lonlat(ivar, outvar_desc(ivar, jg)%type, &
            &                     collected_var_3d, jg, nlev)
        END IF ! (l_interpolate_lonlat)

        ! clean up: collected 3d field on triangular grid no longer needed
        DEALLOCATE(collected_var_3d%ptr, STAT=ierrstat)
        IF (ierrstat /= SUCCESS) THEN
          CALL finish ('mo_io_vlist/write_vlist', 'deallocation failed')
        ENDIF

      ENDDO ! var loop

      IF(my_process_is_stdio()) THEN
        IF (lkeep_in_sync) THEN
          CALL streamSync(streamID(jg))
        END IF
      END IF

    END DO ! loop over patches

    IF(.NOT. my_process_is_stdio()) DEALLOCATE(streamvar1, streamvar2)

    iostep = iostep+1

  END SUBROUTINE write_vlist


  !-------------------------------------------------------------------------
  SUBROUTINE nf(status)
    INTEGER, INTENT(in) :: status

    IF (status /= nf_noerr) THEN
      CALL finish('mo_io_vlist netCDF error', nf_strerror(status))
    ENDIF

  END SUBROUTINE nf


  !-------------------------------------------------------------------------
  !
  !
  !>
  !! For output, the fields have to be provided in the form:.
  !!
  !! For output, the fields have to be provided in the form:
  !! field(vector_length,nlevs).
  !! Therefore, the reshaping performed for optimization should be reversed.
  !! This is the version for a horizontal field.
  !!
  !! The optional parameter "out_field_2d" allows the calling procedure
  !! to work with the collected 2d variable (not only with the RESHAPEd
  !! variable "out_field"). However, in this case the caller is in
  !! charge of deallocating this field.
  !!
  !! @par Revision History
  !! Initial release by Almut Gassmann, MPI-M, (2008-11-19)
  !!
  SUBROUTINE gather_array1(typ,p_patch,in_field,out_field,opt_name,opt_out_field_2d)
    !
    INTEGER,                   INTENT(in)    :: typ
    TYPE(t_patch), TARGET,     INTENT(in)    :: p_patch
    ! 2dimensional input field (nblks,nproma)
    REAL(wp),                  INTENT(in)    :: in_field(:,:)
    ! one vector line version of the input
    REAL(wp),                  INTENT(inout) :: out_field(:)
    CHARACTER(LEN=*), OPTIONAL,INTENT(IN)    :: opt_name
    ! quasi-2d output (nproma,1,nblks)
    TYPE(t_collected_var_ptr), INTENT(OUT), OPTIONAL :: opt_out_field_2d

    REAL(wp), ALLOCATABLE :: out_field2(:,:)
    INTEGER               :: ierrstat

    !-----------------------------------------------------------------------

    IF (my_process_is_stdio()) THEN
      ALLOCATE(out_field2(UBOUND(out_field,1),1))
    ELSE
      ALLOCATE(out_field2(0,0))
    ENDIF
    IF (PRESENT(opt_name)) THEN
      CALL gather_array2(typ,p_patch,&
                         RESHAPE(in_field,(/UBOUND(in_field,1),1,UBOUND(in_field,2)/)), &
                         out_field2, opt_name, opt_out_field_2d)
    ELSE
      CALL gather_array2(typ,p_patch,&
                         RESHAPE(in_field,(/UBOUND(in_field,1),1,UBOUND(in_field,2)/)), &
                         out_field2, opt_out_field_3d=opt_out_field_2d)
    ENDIF

    IF(my_process_is_stdio()) THEN
      out_field(:) = out_field2(:,1)
    ENDIF
    ! if required, return 2D field as (nproma, nblks) array
    IF (.NOT. PRESENT(opt_out_field_2d)) THEN
      DEALLOCATE(opt_out_field_2d%ptr, STAT=ierrstat)
      IF (ierrstat /= SUCCESS) THEN
        CALL finish ('mo_io_vlist/gather_array1', 'deallocation failed')
      ENDIF
    END IF
    DEALLOCATE(out_field2)

  END SUBROUTINE gather_array1


  !-------------------------------------------------------------------------
  !
  !
  !>
  !! For output, the fields have to be provided in the form:
  !!   field(vector_length,nlevs).
  !! Therefore, the reshaping performed for optimization should be reversed.
  !! This is the version for a 3-d field
  !!
  !! The optional parameter "out_field_3d" allows the calling procedure
  !! to work with the collected 3d variable (not only with the RESHAPEd
  !! variable "out_field"). However, in this case the caller is in
  !! charge of deallocating this field.
  !!
  !! @par Revision History
  !! Initial release by Almut Gassmann, MPI-M, (2008-11-19)
  !!
  SUBROUTINE gather_array2(typ,p_patch,in_field,out_field, opt_name, opt_out_field_3d)

    INTEGER,     INTENT(in)                  :: typ
    TYPE(t_patch), INTENT(in), TARGET        :: p_patch
    ! 3d input (nproma,nlev,nblks)
    REAL(wp),    INTENT(in)                  :: in_field(:,:,:)
    ! 2d output (length,nlev)
    REAL(wp),    INTENT(inout)               :: out_field(:,:)
    CHARACTER(LEN=*), OPTIONAL,INTENT(IN)    :: opt_name
    ! 3d output (nproma,nlev,nblks)
    TYPE(t_collected_var_ptr), INTENT(OUT), OPTIONAL :: opt_out_field_3d

    ! local variables
    REAL(wp), POINTER             :: tmp_field(:,:,:)
    INTEGER                       :: nblks, npromz, jb, jl, jk, jend, &
      &                              dim1, dim3, ierrstat
    INTEGER                       :: isize_out, isize_lev               ! array size of output
    TYPE(t_comm_pattern), POINTER :: p_comm_pat

    !-----------------------------------------------------------------------

    IF(UBOUND(in_field,1) /= nproma) THEN
      CALL finish('mo_io_vlist/gather_array2','Illegal 1st array dimension')
    ENDIF
!     IF(p_io/=p_test_pe .AND. p_io/=p_work_pe0) THEN ! Safety check only
!       CALL finish('mo_io_vlist/gather_array2','Illegal I/O PE number for this routine')
!     ENDIF

    IF(typ == GATHER_C) THEN

      IF(UBOUND(in_field,3) /= p_patch%nblks_c) &
        CALL finish('mo_io_vlist/gather_array2','Illegal 3rd array dimension')
      dim1 = nproma
      dim3 = (p_patch%n_patch_cells_g-1)/nproma+1

      p_comm_pat => p_patch%comm_pat_gather_c
      nblks      =  p_patch%nblks_c
      npromz     =  p_patch%npromz_c

    ELSE IF(typ == GATHER_E) THEN

      IF(UBOUND(in_field,3) /= p_patch%nblks_e) &
        CALL finish('mo_io_vlist/gather_array2','Illegal 3rd array dimension')
      dim1 = nproma
      dim3 = (p_patch%n_patch_edges_g-1)/nproma+1

      p_comm_pat => p_patch%comm_pat_gather_e
      nblks      =  p_patch%nblks_e
      npromz     =  p_patch%npromz_e

    ELSE IF(typ == GATHER_V) THEN

      IF(UBOUND(in_field,3) /= p_patch%nblks_v) &
        CALL finish('mo_io_vlist/gather_array2','Illegal 3rd array dimension')
      dim1 = nproma
      dim3 = (p_patch%n_patch_verts_g-1)/nproma+1

      p_comm_pat => p_patch%comm_pat_gather_v
      nblks      =  p_patch%nblks_v
      npromz     =  p_patch%npromz_v

    ELSE

      CALL finish('mo_io_vlist/gather_array2','Illegal type parameter')

      ! To get rid of compiler warnings (by gcc) about variables which may be used uninitialized,
      ! define these varaibles also here. They are not used since the "finish" above stops
      ! the model integration.
      p_comm_pat => p_patch%comm_pat_gather_c
      nblks      =  p_patch%nblks_c
      npromz     =  p_patch%npromz_c

    ENDIF

    IF (PRESENT(opt_out_field_3d)) THEN
      ALLOCATE(opt_out_field_3d%ptr(dim1, UBOUND(in_field,2), dim3), STAT=ierrstat)
      IF (ierrstat /= SUCCESS) THEN
        CALL finish ('mo_io_vlist/gather_array2', 'allocation failed')
      ENDIF
      tmp_field => opt_out_field_3d%ptr
    ELSE
      ALLOCATE(tmp_field(dim1, UBOUND(in_field,2), dim3))
    END IF
    tmp_field(:,:,:)=0.0_wp

    IF(p_test_run) THEN
      IF(.NOT. my_process_is_mpi_test()) THEN
        ! Gather all data on process_mpi_all_workroot_id and send it to process_mpi_test_id for verification
        CALL exchange_data(p_comm_pat, RECV=tmp_field, SEND=in_field)
        IF(my_process_is_mpi_workroot()) CALL p_send(tmp_field, process_mpi_all_test_id, 1)
      ELSE
        ! Receive result from parallel worker PEs and check for correctness
        CALL p_recv(tmp_field, process_mpi_all_workroot_id, 1)
        DO jb = 1, nblks
          jend = nproma
          IF(jb==nblks) jend = npromz
          DO jl = 1, jend
            IF(ANY(tmp_field(jl,:,jb) /= in_field(jl,:,jb))) THEN
                IF (PRESENT(opt_name)) THEN
                  WRITE(0,*)'Error ',opt_name,jl,jb !,tmp_field(jl,:,jb),in_field(jl,:,jb)
                ELSE
                  WRITE(0,*)'Error ',jl,jb !,tmp_field(jl,:,jb),in_field(jl,:,jb)
               ENDIF
              CALL message('mo_io_vlist/gather_array2','Sync error test PE/worker PEs')
            ENDIF
          ENDDO
        ENDDO
      ENDIF
    ELSE
      IF(my_process_is_mpi_seq()) THEN
        ! We are running on 1 PE. Thus, just copy in_field.
        DO jb= 1, nblks
          jend = nproma
          IF(jb==nblks) jend = npromz
          DO jl = 1, jend
            tmp_field(jl,:,jb) = in_field(jl,:,jb)
          ENDDO
        ENDDO
      ELSE
        ! Gather all data on process_mpi_all_workroot_id
        CALL exchange_data(p_comm_pat, RECV=tmp_field, SEND=in_field)
      ENDIF
    ENDIF

    IF(my_process_is_stdio()) THEN
      isize_out = SIZE(out_field,1)
      isize_lev = SIZE(in_field,2)

      DO jk = 1, isize_lev
        out_field(:,jk) = RESHAPE(tmp_field(:,jk,:),(/isize_out/))
      ENDDO
    ENDIF

    IF (.NOT. PRESENT(opt_out_field_3d)) &
      DEALLOCATE(tmp_field)

  END SUBROUTINE gather_array2


  !-------------------------------------------------------------------------
  !-------------------------------------------------------------------------
  !
  ! Below are the original routines de_reshape1/de_reshape2
  ! These are still used from other modules
  !
  !-------------------------------------------------------------------------
  !-------------------------------------------------------------------------
  !
  !
  !>
  !! For output, the fields have to be provided in the form:.
  !!
  !! For output, the fields have to be provided in the form:
  !! field(vector_length,nlevs).
  !! Therefore, the reshaping performed for optimization should be reversed.
  !! This is the version for a horizontal field.
  !!
  !! @par Revision History
  !! Initial release by Almut Gassmann, MPI-M, (2008-11-19)
  !!
  SUBROUTINE de_reshape1(in_field,out_field)
    !

    REAL(wp), INTENT(in)    :: in_field(:,:)    ! 2dimensional input field (nblks,nproma)

    REAL(wp), INTENT(inout) :: out_field(:)   ! one vector line version of the input

    LOGICAL, ALLOCATABLE    :: lmask(:)
    INTEGER :: idiscrep, isize_in, isize_out ! array sizes of input and output,
    ! and their discrepancy

    !-----------------------------------------------------------------------

    isize_in  = SIZE(in_field)
    isize_out = SIZE(out_field)
    idiscrep = isize_in-isize_out

    IF(idiscrep /= 0 )THEN
      ALLOCATE (lmask(isize_in))
      lmask(1:isize_out) = .TRUE.
      lmask(isize_out+1:isize_in) = .FALSE.
      out_field = PACK(RESHAPE(in_field,(/isize_in/)),lmask)
      DEALLOCATE (lmask)
    ELSE
      out_field = RESHAPE( in_field,(/isize_out/))
    ENDIF

  END SUBROUTINE de_reshape1


  !-------------------------------------------------------------------------
  !
  !
  !>
  !! For output, the fields have to be provided in the form:.
  !!
  !! For output, the fields have to be provided in the form:
  !! field(vector_length,nlevs).
  !! Therefore, the reshaping performed for optimization should be reversed.
  !! This is the version for a 3-d field
  !!
  !! @par Revision History
  !! Initial release by Almut Gassmann, MPI-M, (2008-11-19)
  !!
  SUBROUTINE de_reshape2(in_field,out_field)
    !

    REAL(wp), INTENT(in)    :: in_field(:,:,:)  ! 3d input (nproma,nlev,nblks)

    REAL(wp), INTENT(inout) :: out_field(:,:)   ! 2d output (length,nlev)

    LOGICAL, ALLOCATABLE    :: lmask(:)
    INTEGER :: idiscrep, isize_in, isize_out, isize_lev, jk
    ! array sizes of input and output, and their discrepancy

    !-----------------------------------------------------------------------

    isize_in  = SIZE(in_field,1)*SIZE(in_field,3)
    isize_out = SIZE(out_field,1)
    isize_lev = SIZE(in_field,2)
    idiscrep = isize_in-isize_out

    IF (idiscrep /= 0 )THEN
      ALLOCATE (lmask(isize_in))
      lmask(1:isize_out) = .TRUE.
      lmask(isize_out+1:isize_in) = .FALSE.
    ENDIF

    DO jk = 1, isize_lev
      IF (idiscrep /= 0 )THEN
        out_field(:,jk) = PACK(RESHAPE(in_field(:,jk,:),(/isize_in/)),lmask)
      ELSE
        out_field(:,jk) =      RESHAPE(in_field(:,jk,:),(/isize_out/))
      ENDIF
    ENDDO

    IF (idiscrep /= 0 )THEN
      DEALLOCATE (lmask)
    ENDIF

  END SUBROUTINE de_reshape2


  !-------------------------------------------------------------------------------------------------
  SUBROUTINE addGlobAttInt(att_name, att_int, vlist, astatus)
    CHARACTER(*), INTENT(IN)  :: att_name
    INTEGER     , INTENT(IN)  :: att_int, vlist
    INTEGER     , INTENT(OUT) :: astatus

    astatus = vlistDefAttInt(vlist,CDI_GLOBAL,TRIM(att_name),DATATYPE_INT32,1,att_int)
  END SUBROUTINE addGlobAttInt


  !-------------------------------------------------------------------------------------------------
  SUBROUTINE addGlobAttTxt(att_name, att_txt, vlist, astatus)
    CHARACTER(*), INTENT(IN)  :: att_name, att_txt
    INTEGER     , INTENT(IN)  :: vlist
    INTEGER     , INTENT(OUT) :: astatus

    astatus = vlistDefAttTxt(vlist,CDI_GLOBAL,TRIM(att_name),LEN(TRIM(att_txt)),TRIM(att_txt))
  END SUBROUTINE addGlobAttTxt


  !-------------------------------------------------------------------------------------------------
  SUBROUTINE addGlobAttTxtFromLog(att_name, boolian, vlist, astatus)
    CHARACTER(*), INTENT(IN)  :: att_name
    LOGICAL, INTENT(IN)       :: boolian
    INTEGER     , INTENT(IN)  :: vlist
    INTEGER     , INTENT(OUT) :: astatus
    CALL addGlobAttTxt(att_name, TRIM(MERGE('.true. ','.false.',boolian)), vlist, astatus)
  END SUBROUTINE addGlobAttTxtFromLog


  !-------------------------------------------------------------------------------------------------
  SUBROUTINE addGlobAttFlt(att_name, att_flt, vlist, astatus)
    CHARACTER(*), INTENT(IN)  :: att_name
    REAL(wp)    , INTENT(IN)  :: att_flt
    INTEGER     , INTENT(IN)  :: vlist
    INTEGER     , INTENT(OUT) :: astatus

    astatus = vlistDefAttFlt(vlist,CDI_GLOBAL,TRIM(att_name),DATATYPE_FLT32,1,att_flt)
  END SUBROUTINE addGlobAttFlt


  !-------------------------------------------------------------------------------------------------
  SUBROUTINE addVar(var,k_jg)
    INTEGER, INTENT(IN)    :: var, k_jg
    ! local variables
    CHARACTER(len=max_name_len) :: zname
    INTEGER                     :: ivar_idx

    ! get variable name
    CALL vlistInqVarName(vlistID(k_jg), var, zname);

    ! check if this variable is due to lon-lat interpolation
    LONLAT : IF (.NOT. lonlat_intp_config(k_jg)%l_enabled) THEN 
      num_varids(k_jg)      = num_varids(k_jg) + 1
      ivar_idx              = num_varids(k_jg)
      varids(ivar_idx,k_jg) = var
    ELSE
      IF ((toupper(lonlat_intp_config(k_jg)%zlist) == 'ALL') .OR.  &
        & (string_contains_word(zname,                             &
        &     lonlat_intp_config(k_jg)%zlist,                      &
        &     lonlat_intp_config(k_jg)%n_list,                     &
        &     lonlat_intp_config(k_jg)%pos_list,                   &
        &     lonlat_intp_config(k_jg)%ilen_list))) THEN
        CALL addVar_lonlat(var, k_jg)
      ELSE
        num_varids(k_jg)      = num_varids(k_jg) + 1
        ivar_idx              = num_varids(k_jg)
        varids(ivar_idx,k_jg) = var
      END IF
    END IF LONLAT

  END SUBROUTINE addVar


  !-------------------------------------------------------------------------------------------------
  !> defines a new output variable which is then also interpolated to
  !> a lon-lat grid
  SUBROUTINE addVar_lonlat(varID, k_jg)
    INTEGER, INTENT(IN)     :: varID   ! variable ID
    INTEGER, INTENT(IN)     :: k_jg    ! patch index
    ! local variables
    CHARACTER(*), PARAMETER :: routine = TRIM("mo_io_vlist:addVar_lonlat")
    INTEGER                 :: igridID, ivar_idx

    ! add standard variable
    num_varids(k_jg)       = num_varids(k_jg) + 1
    ivar_idx               = num_varids(k_jg)
    varids(ivar_idx,k_jg)  = varID
    
    ! determine the grid on which the variable is defined:
    igridID = vlistInqVarGrid(vlistID(k_jg), varID)

    ! depending on the grid on which the variable is defined, add one
    ! or two (gradients) variables for lon-lat grid
    IF (igridID == gridCellID(k_jg)) THEN

      num_varids(k_jg)              = num_varids(k_jg) + 1
      varids(num_varids(k_jg),k_jg) = var_clone(vlistID(k_jg), varID, &
        &                             lonlat_gridID(k_jg), "_LONLAT")
      ! store association between variables
      lonlat_varids(ivar_idx,k_jg,1) = varids(num_varids(k_jg),k_jg)

    ELSE IF (igridID == gridEdgeID(k_jg)) THEN
      
      ! add x-component:
      num_varids(k_jg)               = num_varids(k_jg) + 1
      varids(num_varids(k_jg),k_jg)  = var_clone(vlistID(k_jg), varID, &
        &                             lonlat_gridID(k_jg), "_LONLAT_X")
      ! store association between variables
      lonlat_varids(ivar_idx,k_jg,1) = varids(num_varids(k_jg),k_jg)

      ! add y-component:
      num_varids(k_jg)               = num_varids(k_jg) + 1
      varids(num_varids(k_jg),k_jg)  = var_clone(vlistID(k_jg), varID, &
        &                             lonlat_gridID(k_jg), "_LONLAT_Y")
      ! store association between variables
      lonlat_varids(ivar_idx,k_jg,2) = varids(num_varids(k_jg),k_jg)

    END IF

  END SUBROUTINE addVar_lonlat


  !-------------------------------------------------------------------------------------------------
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


  !-------------------------------------------------------------------------------------------------
  FUNCTION ConstVar(vname,vlongname,vunit,vcode,vtable,vlist,grid,zaxis) RESULT(var)
    INTEGER                  :: var
    CHARACTER(*), INTENT(IN) :: vname, vlongname, vunit
    INTEGER, INTENT(IN)      :: vcode, vtable, vlist, grid, zaxis

    var = vlistdefvar(vlist, grid, zaxis, TIME_CONSTANT)
    CALL vlistdefvarname(vlist, var, vname)
    CALL vlistdefvarlongname (vlist, var, vlongname)
    CALL vlistdefvarunits(vlist, var, vunit)
    IF ( vcode .gt. 0 ) THEN
      CALL vlistdefvarcode(vlist, var, vcode)
    ELSE
      CALL message('WARNING:TimeVar','Prevent setting negative var code for'//TRIM(vname))
    END IF
  END FUNCTION ConstVar


  !-------------------------------------------------------------------------------------------------
  FUNCTION DebugVar(vname,vlist,grid,zaxis) RESULT(var)
    INTEGER                  :: var
    CHARACTER(*), INTENT(IN) :: vname
    INTEGER, INTENT(IN)      :: vlist, grid, zaxis

    var  = vlistdefvar(vlist, grid, zaxis, TIME_VARIABLE)
    CALL vlistdefvarname(vlist, var, vname)

  END FUNCTION DebugVar


  !-------------------------------------------------------------------------------------------------
  !> Duplicates a given "ConstVar", "TimeVar" or "DebugVar"
  FUNCTION var_clone(vlistID, varID, new_gridID, suffix) RESULT(var)

    INTEGER, INTENT(IN)      :: vlistID, varID, new_gridID
    CHARACTER(*), INTENT(IN) :: suffix
    INTEGER                  :: var
    ! local variables
    INTEGER                      :: gridID, zaxisID, timeID, vcode
    CHARACTER(len=max_name_len)  :: vname, vlongname, vunits

    ! inquire about vlist, zaxis, TIME_VARIABLE/TIME_CONSTANT, vname,
    ! vlongname, vunit, vcode
    CALL vlistInqVar(vlistID, varID, gridID, zaxisID, timeID)
    CALL vlistInqVarName(vlistID, varID, vname)
    CALL vlistInqVarLongname(vlistID, varID, vlongname)
    CALL vlistInqVarUnits(vlistID, varID, vunits);
    vcode = vlistInqVarCode(vlistID, varID);

    var = vlistdefvar(vlistID, new_gridID, zaxisID, timeID)
    CALL vlistdefvarname(vlistID, var, TRIM(TRIM(vname)//suffix))

    CALL vlistdefvarlongname (vlistID, var, vlongname)
    CALL vlistdefvarunits(vlistID, var, vunits)
    IF ( vcode .GT. 0 ) THEN
      CALL vlistdefvarcode(vlistID, var, vcode)
    ELSE
      CALL message('WARNING:var_clone','Prevent setting negative var code for'//TRIM(vname))
    END IF

  END FUNCTION var_clone


  !-------------------------------------------------------------------------------------------------
  ! dup3: duplicates an array and returns a pointer to the duplicate
  FUNCTION dup3(arr)
    REAL(wp), INTENT(IN) :: arr(:,:,:)
    REAL(wp), POINTER :: dup3(:,:,:)

    ALLOCATE(dup3(size(arr,1),size(arr,2),size(arr,3)))

    dup3 = arr

  END FUNCTION dup3


  !-------------------------------------------------------------------------------------------------
  ! dup2: duplicates an array and returns a pointer to the duplicate
  FUNCTION dup2(arr)
    REAL(wp), INTENT(IN) :: arr(:,:)
    REAL(wp), POINTER :: dup2(:,:)

    ALLOCATE(dup2(size(arr,1),size(arr,2)))

    dup2 = arr

  END FUNCTION dup2


  !-------------------------------------------------------------------------------------------------
  SUBROUTINE addGlobAtts(vlist,k_jg,astatus)
    INTEGER, INTENT(IN) :: vlist, k_jg
    INTEGER             :: astatus
    ! Model name and version
    ! ----------------------
    CALL addGlobAttTxt('model-version',TRIM(modelname)//'-'//TRIM(modelversion),&
    &                  vlist,astatus)
    CALL addGlobAttInt('mpi:num_work_procs',num_work_procs,vlist,astatus)
    !
    ! Parameters of /grid_nml/
    ! ------------------------
    CALL addGlobAttInt('grid_nml:nroot',nroot,vlist,astatus)
    CALL addGlobAttInt('grid_nml:start_lev',start_lev,vlist,astatus)
    CALL addGlobAttInt('grid_nml:n_dom',n_dom,vlist,astatus)
    CALL addGlobAttTxtFromLog('grid_nml:lfeedback', lfeedback(k_jg),vlist,astatus)
    CALL addGlobAttTxtFromLog('grid_nml:lplane',lplane,vlist,astatus)
    CALL addGlobAttTxtFromLog('grid_nml:lvert_nest',lvert_nest,vlist,astatus)
    !
    ! Parameters of /run_nml/
    ! -----------------------
    CALL addGlobAttInt('run_nml:global_cell_type',global_cell_type,vlist,astatus)
    CALL addGlobAttInt('run_nml:num_lev',num_lev(k_jg),vlist,astatus)
    CALL addGlobAttFlt('run_nml:dtime',dtime,vlist,astatus)
    CALL addGlobAttInt('run_nml:iequations',iequations,vlist,astatus)
    CALL addGlobAttInt('run_nml:ntracer',ntracer,vlist,astatus)
    CALL addGlobAttTxtFromLog('run_nml:ldynamics',ldynamics,vlist,astatus)
    CALL addGlobAttTxtFromLog('run_nml:ltransport',ltransport,vlist,astatus)
    CALL addGlobAttInt('run_nml:iforcing',iforcing,vlist,astatus)
    CALL addGlobAttTxtFromLog('run_nml:ltestcase',ltestcase,vlist,astatus)
    CALL addGlobAttInt('run_nml:nproma',nproma,vlist,astatus)
    CALL addGlobAttInt('run_nml:nsteps',nsteps,vlist,astatus)
    CALL addGlobAttInt('run_nml:itopo',itopo,vlist,astatus)
    CALL addGlobAttInt('run_nml:msg_level',msg_level,vlist,astatus)
    !
    ! Parameters of /io_nml/
    ! ----------------------
    CALL addGlobAttTxt('io_nml:out_expname',TRIM(out_expname),vlist,astatus)
    CALL addGlobAttFlt('io_nml:dt_data',dt_data,vlist,astatus)
    CALL addGlobAttFlt('io_nml:dt_file',dt_file,vlist,astatus)
    CALL addGlobAttInt('io_nml:out_filetype',out_filetype,vlist,astatus)
    CALL addGlobAttTxtFromLog('io_nml:lwrite_vorticity',lwrite_vorticity,vlist,astatus)
    CALL addGlobAttTxtFromLog('io_nml:lwrite_pres',lwrite_pres,vlist,astatus)
    CALL addGlobAttTxtFromLog('io_nml:lwrite_z3',lwrite_z3,vlist,astatus)
    CALL addGlobAttTxtFromLog('io_nml:lwrite_omega',lwrite_omega,vlist,astatus)
    !
    ! Parameters of /dynamics_nml/
    ! ----------------------------
    CALL addGlobAttInt('dynamics_nml:idiv_method',idiv_method,vlist,astatus)
    CALL addGlobAttFlt('dynamics_nml:divavg_cntrwgt',divavg_cntrwgt,vlist,astatus)
    CALL addGlobAttTxtFromLog('dynamics_nml:lcoriolis',lcoriolis, &
                              vlist,astatus)

  END SUBROUTINE addGlobAtts


  !-------------------------------------------------------------------------------------------------
  SUBROUTINE addAtmAtts(vlist,k_jg,astatus)
    INTEGER, INTENT(IN) :: vlist,k_jg
    INTEGER             :: astatus

    INTEGER             :: jt
    CHARACTER(len=3)    :: cjt
    !----------------------------
    ! namelist/ha_dyn_nml/
    !----------------------------
     IF (iequations == 1 .OR. iequations == 2) THEN

       CALL addGlobAttInt('ha_dyn_nml:itime_scheme', &
            ha_dyn_config%itime_scheme,vlist,astatus)
       CALL addGlobAttTxtFromLog('ha_dyn_nml:lsi_3d',&
            ha_dyn_config%lsi_3d,vlist,astatus)
       CALL addGlobAttInt('ha_dyn_nml:ileapfrog_startup', &
            ha_dyn_config%ileapfrog_startup,vlist,astatus)
       CALL addGlobAttFlt('ha_dyn_nml:asselin_coeff',&
            ha_dyn_config%asselin_coeff,vlist,astatus)
       CALL addGlobAttFlt('ha_dyn_nml:si_2tls',&
            ha_dyn_config%si_2tls,vlist,astatus)
       CALL addGlobAttFlt('ha_dyn_nml:si_coeff',&
            ha_dyn_config%si_coeff,vlist,astatus)
       CALL addGlobAttFlt('ha_dyn_nml:si_offctr',&
            ha_dyn_config%si_offctr,vlist,astatus)
       CALL addGlobAttFlt('ha_dyn_nml:si_rtol',&
            ha_dyn_config%si_rtol,vlist,astatus)
       CALL addGlobAttFlt('ha_dyn_nml:si_cmin',&
            ha_dyn_config%si_cmin,vlist,astatus)
       CALL addGlobAttInt('ha_dyn_nml:si_expl_scheme',&
            ha_dyn_config%si_expl_scheme,vlist,astatus)
       CALL addGlobAttTxtFromLog('ha_dyn_nml:ldry_dycore',&
            ha_dyn_config%ldry_dycore,vlist,astatus)
       CALL addGlobAttTxtFromLog('ha_dyn_nml:lref_temp',&
            ha_dyn_config%lref_temp,vlist,astatus)
     ENDIF


     !
     ! Parameters of /nonhydrostatic_nml/
     ! ----------------------------
     IF (iequations == 3) THEN
        CALL addGlobAttInt('nonhydrostatic_nml:itime_scheme',&
                          & itime_scheme_nh_atm,vlist,astatus)
        CALL addGlobAttInt('nonhydrostatic_nml:ivctype',ivctype,vlist,astatus)
        CALL addGlobAttInt('nonhydrostatic_nml:iadv_rcf',iadv_rcf,vlist,astatus)

        IF (global_cell_type == 3) THEN
          CALL addGlobAttFlt('nonhydrostatic_nml:rayleigh_coeff',   &
                  &         rayleigh_coeff(k_jg),vlist,astatus)
          CALL addGlobAttFlt('nonhydrostatic_nml:damp_height',      &
                  &         damp_height(k_jg),vlist,astatus)
          CALL addGlobAttFlt('nonhydrostatic_nml:vwind_offctr',     &
                  &         vwind_offctr,vlist,astatus)
          CALL addGlobAttTxtFromLog('nonhydrostatic_nml:l_nest_rcf',&
                  &         l_nest_rcf,vlist,astatus)
          CALL addGlobAttInt('nonhydrostatic_nml:iadv_rhotheta',    &
                  &         iadv_rhotheta,vlist,astatus)
          CALL addGlobAttInt('nonhydrostatic_nml:igradp_method',    &
                  &         igradp_method,vlist,astatus)
          CALL addGlobAttFlt('nonhydrostatic_nml:exner_expol',      &
                  &         exner_expol,vlist,astatus)
          CALL addGlobAttTxtFromLog('nonhydrostatic_nml:l_open_ubc',&
                  &         l_open_ubc,vlist,astatus)
        ELSEIF (global_cell_type == 6) THEN
          CALL addGlobAttTxtFromLog('nonhydrostatic_nml:ltheta_up_hori', &
                  &         ltheta_up_hori,vlist,astatus)
          CALL addGlobAttTxtFromLog('nonhydrostatic_nml:ltheta_up_vert',&
                  &         ltheta_up_vert,vlist,astatus)
          CALL addGlobAttFlt('nonhydrostatic_nml:gmres_rtol_nh',      &
                  &         gmres_rtol_nh,vlist,astatus)
          CALL addGlobAttFlt('nonhydrostatic_nml:upstr_beta',      &
                  &         upstr_beta,vlist,astatus)

        ENDIF
     END IF
     !
     ! Parameters of /diffusion_nml/
     ! -----------------------------
     CALL addGlobAttInt('diffusion_nml:hdiff_order',            &
       &  diffusion_config(k_jg)%hdiff_order,vlist,astatus)
     CALL addGlobAttFlt('diffusion_nml:hdiff_efdt_ratio',       &
       &  diffusion_config(k_jg)%hdiff_efdt_ratio,vlist,astatus)
     CALL addGlobAttFlt('diffusion_nml:hdiff_multfac',          &
       &  diffusion_config(k_jg)%hdiff_multfac,vlist,astatus)
     CALL addGlobAttFlt('diffusion_nml:hdiff_tv_ratio',         &
       &  diffusion_config(k_jg)%hdiff_tv_ratio,vlist,astatus)
     CALL addGlobAttTxtFromLog('diffusion_nml:lhdiff_temp',     &
       &  diffusion_config(k_jg)%lhdiff_temp,vlist,astatus)
     CALL addGlobAttTxtFromLog('diffusion_nml:lhdiff_vn',       &
       &  diffusion_config(k_jg)%lhdiff_vn,vlist,astatus)
     CALL addGlobAttFlt('diffusion_nml:hdiff_smag_fac',         &
       &  diffusion_config(k_jg)%hdiff_smag_fac,vlist,astatus)
     !
     ! Parameters of /transport_nml/
     ! -----------------------------
     IF (ltransport) THEN
        CALL addGlobAttInt('transport_nml:ihadv_tracer',       &
          &  advection_config(k_jg)%ihadv_tracer(max_ntracer),vlist,astatus)
        CALL addGlobAttInt('transport_nml:ivadv_tracer',       &
          &  advection_config(k_jg)%ivadv_tracer(max_ntracer),vlist,astatus)
        CALL addGlobAttTxtFromLog('transport_nml:lvadv_tracer',&
          &  advection_config(k_jg)%lvadv_tracer,vlist,astatus)
        CALL addGlobAttTxtFromLog('transport_nml:lstrang',     &
          &  advection_config(k_jg)%lstrang,vlist,astatus)
        CALL addGlobAttInt('transport_nml:itype_vlimit',       &
          &  advection_config(k_jg)%itype_vlimit(max_ntracer),vlist,astatus)
        CALL addGlobAttInt('transport_nml:itype_hlimit',       &
          &  advection_config(k_jg)%itype_hlimit(max_ntracer),vlist,astatus)
        CALL addGlobAttInt('transport_nml:iord_backtraj',      &
          &  advection_config(k_jg)%iord_backtraj,vlist,astatus)
        CALL addGlobAttInt('transport_nml:igrad_c_miura',      &
          &  advection_config(k_jg)%igrad_c_miura,vlist,astatus)
        CALL addGlobAttTxtFromLog('transport_nml:lclip_tracer',&
          &  advection_config(k_jg)%lclip_tracer,vlist,astatus)
        CALL addGlobAttInt('transport_nml:iadv_slev',          &
          &  advection_config(k_jg)%iadv_slev(max_ntracer),vlist,astatus)
        CALL addGlobAttFlt('transport_nml:upstr_beta_adv',     &
          &  advection_config(k_jg)%upstr_beta_adv,vlist,astatus)
     END IF
     !
     IF (ntracer > 0) THEN
       DO jt=1,ntracer
         WRITE(cjt,'(i3.3)') jt
         CALL addGlobAttTxtFromLog('io_nml:lwrite_tracer(' // TRIM(ADJUSTL(cjt)) // ')', &
         &                  lwrite_tracer(jt),vlist,astatus)
       END DO
     END IF
     CALL addGlobAttTxtFromLog('io_nml:lwrite_precip',lwrite_precip,vlist,astatus)
     CALL addGlobAttTxtFromLog('io_nml:lwrite_cloud',lwrite_cloud,vlist,astatus)
     CALL addGlobAttTxtFromLog('io_nml:lwrite_divergence',lwrite_divergence,vlist,astatus)


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! Parameters of physics parameterizations
  !------------------------------------------
     IF (lforcing) THEN
       SELECT CASE (iforcing)

       CASE (iecham,ildf_echam)
         !
         !!! Parameters of /echam_phy_nml/
         !--------------------------------
         CALL addGlobAttTxtFromLog('echam_phy_nml:lrad',get_lrad(),vlist,astatus)
         CALL addGlobAttTxtFromLog('echam_phy_nml:lvdiff',get_lvdiff(),vlist,astatus)
         CALL addGlobAttTxtFromLog('echam_phy_nml:lconv',get_lconv(),vlist,astatus)
         CALL addGlobAttTxtFromLog('echam_phy_nml:lcond',get_lcond(),vlist,astatus)
         CALL addGlobAttTxtFromLog('echam_phy_nml:lcover',get_lcover(),vlist,astatus)
         CALL addGlobAttTxtFromLog('echam_phy_nml:lgw_hines',get_lgw_hines(),vlist,astatus)
         !
         !!! Parameters of /echam_conv_nml/
         !---------------------------------
         IF ( get_lconv() ) THEN
            CALL addGlobAttTxtFromLog('echam_conv_nml:lmfpen', &
                 echam_conv_config%lmfpen,vlist,astatus)
            CALL addGlobAttTxtFromLog('echam_conv_nml:lmfmid', &
                 echam_conv_config%lmfmid,vlist,astatus)
           !CALL addGlobAttTxtFromLog('echam_conv_nml:lmfscv', &
           !     echam_conv_config%lmfscv,vlist,astatus)
            CALL addGlobAttTxtFromLog('echam_conv_nml:lmfdd',  &
                 echam_conv_config%lmfdd,vlist,astatus)
            CALL addGlobAttTxtFromLog('echam_conv_nml:lmfdudv',&
                 echam_conv_config%lmfdudv,vlist,astatus)
            CALL addGlobAttInt('echam_conv_nml:iconv', &
                 echam_conv_config%iconv,vlist,astatus)
            CALL addGlobAttFlt('echam_conv_nml:cmftau', &
                 echam_conv_config%cmftau,vlist,astatus)
            CALL addGlobAttFlt('echam_conv_nml:cmfctop', &
                 echam_conv_config%cmfctop,vlist,astatus)
            CALL addGlobAttFlt('echam_conv_nml:cprcon', &
                 echam_conv_config%cprcon,vlist,astatus)
            CALL addGlobAttFlt('echam_conv_nml:cminbuoy',&
                 echam_conv_config%cminbuoy,vlist,astatus)
            CALL addGlobAttFlt('echam_conv_nml:entrpen', &
                 echam_conv_config%entrpen,vlist,astatus)
            CALL addGlobAttFlt('echam_conv_nml:dlev',    &
                 echam_conv_config%dlev,vlist,astatus)
         END IF
      !
!!$        !!! Parameters of /gw_hines_nml/
!!$        !-------------------------------
!!$        IF (lgw_hines) THEN
!!$          CALL addGlobAttTxtFromLog('gw_hines_nml:lheatcal',lheatcal,vlist,astatus)
!!$          CALL addGlobAttInt('gw_hines_nml:emiss_lev',emiss_lev,vlist,astatus)
!!$          CALL addGlobAttFlt('gw_hines_nml:rmscon',rmscon,vlist,astatus)
!!$          CALL addGlobAttFlt('gw_hines_nml:kstar',kstar,vlist,astatus)
!!$          CALL addGlobAttFlt('gw_hines_nml:m_min',m_min,vlist,astatus)
!!$        END IF

       CASE (inwp)
         !!! Parameters of /nwp_phy_nml/
         !-------------------------------
        CALL addGlobAttInt('nwp_phy_nml:inwp_gscp',atm_phy_nwp_config(k_jg)%inwp_gscp,&
          &vlist,astatus)
        CALL addGlobAttInt('nwp_phy_nml:inwp_convection',atm_phy_nwp_config(k_jg)%inwp_convection,&
          &vlist,astatus)
        CALL addGlobAttInt('nwp_phy_nml:inwp_cldcover',atm_phy_nwp_config(k_jg)%inwp_cldcover,&
          &vlist,astatus)
        CALL addGlobAttInt('nwp_phy_nml:inwp_radiation',atm_phy_nwp_config(k_jg)%inwp_radiation,&
          &vlist,astatus)
        CALL addGlobAttInt('nwp_phy_nml:inwp_satad',atm_phy_nwp_config(k_jg)%inwp_satad,&
          &vlist,astatus)
        CALL addGlobAttInt('nwp_phy_nml:inwp_turb',atm_phy_nwp_config(k_jg)%inwp_turb,&
          &vlist,astatus)
        CALL addGlobAttInt('nwp_phy_nml:inwp_sso',atm_phy_nwp_config(k_jg)%inwp_sso,&
          &vlist,astatus)
        CALL addGlobAttFlt('nwp_phy_nml:dt_conv',atm_phy_nwp_config(k_jg)%dt_conv ,&
          &vlist,astatus)
        CALL addGlobAttFlt('nwp_phy_nml:dt_rad',atm_phy_nwp_config(k_jg)%dt_rad  ,&
          &vlist,astatus)
        CALL addGlobAttFlt('nwp_phy_nml:dt_sso',atm_phy_nwp_config(k_jg)%dt_sso ,&
          &vlist,astatus)
       END SELECT
     END IF

     !
     ! Parameters of /radiation_nml/
     ! -----------------------------
     CALL addGlobAttFlt('radiation_nml:dt_rad',atm_phy_nwp_config(k_jg)%dt_rad,&
       &               vlist,astatus)
     CALL addGlobAttInt('radiation_nml:izenith',izenith,vlist,astatus)
     CALL addGlobAttInt('radiation_nml:irad_h2o',irad_h2o,vlist,astatus)
     CALL addGlobAttInt('radiation_nml:irad_co2',irad_co2,vlist,astatus)
     CALL addGlobAttInt('radiation_nml:irad_ch4',irad_ch4,vlist,astatus)
     CALL addGlobAttInt('radiation_nml:irad_n2o',irad_n2o,vlist,astatus)
     CALL addGlobAttInt('radiation_nml:irad_o3',irad_o3,vlist,astatus)
     CALL addGlobAttInt('radiation_nml:irad_o2',irad_o2,vlist,astatus)
     CALL addGlobAttInt('radiation_nml:irad_cfc11',irad_cfc11,vlist,astatus)
     CALL addGlobAttInt('radiation_nml:irad_cfc12',irad_cfc12,vlist,astatus)
     CALL addGlobAttInt('radiation_nml:irad_aero',irad_aero,vlist,astatus)
  !
  ! Parameters of testcases (hydrostatic and non-hydrostatic
  ! --------------------------------------------------------
     IF (ltestcase) THEN
       !
       ! Parameters of /testcase_nml/
       ! -----------------------------
       IF(iequations == 0 .OR. &
          iequations == 1 .OR. &
          iequations == 2 ) THEN
       !
          CALL addGlobAttTxt('testcase_nml:ctest_name',TRIM(ctest_name),vlist,astatus)
          !
          IF ( TRIM(ctest_name) == 'GW') THEN
             CALL addGlobAttFlt('testcase_nml:gw_brunt_vais',gw_brunt_vais,vlist,astatus)
             CALL addGlobAttFlt('testcase_nml:gw_u0',gw_u0,vlist,astatus)
             CALL addGlobAttFlt('testcase_nml:gw_lon_deg',gw_lon_deg,vlist,astatus)
             CALL addGlobAttFlt('testcase_nml:gw_lat_deg',gw_lat_deg,vlist,astatus)
          !
          ELSEIF ( TRIM(ctest_name) == 'JWw' ) THEN
             CALL addGlobAttFlt('testcase_nml:jw_uptb',jw_uptb,vlist,astatus)
          !
          ELSEIF ( TRIM(ctest_name) == 'MRW' .OR. TRIM(ctest_name) == 'MRW2'  ) THEN
             CALL addGlobAttFlt('testcase_nml:mountctr_lon_deg',  &
                       &        mountctr_lon_deg,vlist,astatus)
             CALL addGlobAttFlt('testcase_nml:mountctr_lat_deg',  &
                       &        mountctr_lat_deg,vlist,astatus)
             CALL addGlobAttFlt('testcase_nml:mountctr_height',   &
                       &        mountctr_height,vlist,astatus)
             CALL addGlobAttFlt('testcase_nml:mount_half_width',  &
                       &        mount_half_width,vlist,astatus)
             CALL addGlobAttFlt('testcase_nml:mount_u0',          &
                       &        mount_u0,vlist,astatus)
          !
          ELSEIF ( TRIM(ctest_name) == 'RH' ) THEN
             CALL addGlobAttInt('testcase_nml:rh_wavenum',rh_wavenum,vlist,astatus)
             CALL addGlobAttFlt('testcase_nml:rh_init_shift_deg', &
                       &        rh_init_shift_deg,vlist,astatus)
          !
          ELSEIF ( TRIM(ctest_name) == 'HS' ) THEN
             CALL addGlobAttInt('testcase_nml:ihs_init_type',ihs_init_type,vlist,astatus)
             CALL addGlobAttTxtFromLog('testcase_nml:lhs_vn_ptb',lhs_vn_ptb,vlist,astatus)
             CALL addGlobAttFlt('testcase_nml:hs_vn_ptb_scale',   &
                       &        hs_vn_ptb_scale,vlist,astatus)
          !
          ELSEIF ( TRIM(ctest_name) == 'JWw-Moist' .OR. TRIM(ctest_name) == 'APE' &
                   .OR. TRIM(ctest_name) == 'LDF-Moist' ) THEN
             CALL addGlobAttTxtFromLog('testcase_nml:lrh_linear_pres',&
                       &        lrh_linear_pres,vlist,astatus)
             CALL addGlobAttFlt('testcase_nml:rh_at_1000hpa',         &
                       &        rh_at_1000hpa,vlist,astatus)
          !
          ELSEIF ( TRIM(ctest_name) == 'PA' ) THEN
             CALL addGlobAttTxtFromLog('testcase_nml:linit_tracer_fv',&
                       &        linit_tracer_fv,vlist,astatus)
          !
          END IF

       !
       ! Parameters of /nh_testcase_nml/
       ! -----------------------------
       ELSEIF (iequations == 3) THEN
       !
          CALL addGlobAttTxt('nh_testcase_nml:nh_test_name', &
                   &         TRIM(nh_test_name),vlist,astatus)
          !
          IF ( TRIM(nh_test_name) == 'jabw') THEN
             CALL addGlobAttFlt('nh_testcase_nml:jw_up',jw_up,vlist,astatus)
          !
          ELSEIF ( TRIM(nh_test_name) == 'mrw2_nh') THEN
             CALL addGlobAttFlt('nh_testcase_nml:u0_mrw',u0_mrw,vlist,astatus)
          !
          ELSEIF ( TRIM(nh_test_name) == 'mrw2_nh' .OR. TRIM(nh_test_name) == 'mwbr_const') THEN
             CALL addGlobAttFlt('nh_testcase_nml:mount_height_mrw',     &
                  &             mount_height_mrw,vlist,astatus)
             !CALL addGlobAttFlt('nh_testcase_nml:mount_half_width',     &
             !     &             mount_half_width,vlist,astatus)
             CALL addGlobAttFlt('nh_testcase_nml:mount_lonctr_mrw_deg', &
                  &             mount_lonctr_mrw_deg,vlist,astatus)
             CALL addGlobAttFlt('nh_testcase_nml:mount_latctr_mrw_deg', &
                  &             mount_latctr_mrw_deg,vlist,astatus)
          !
          ELSEIF ( TRIM(nh_test_name) == 'mwbr_const') THEN
             CALL addGlobAttFlt('nh_testcase_nml:u0_mrw',                &
                  &             u0_mrw,vlist,astatus)
             CALL addGlobAttFlt('nh_testcase_nml:temp_i_mwbr_const',     &
                  &             temp_i_mwbr_const,vlist,astatus)
             CALL addGlobAttFlt('nh_testcase_nml:p_int_mwbr_const',      &
                  &             p_int_mwbr_const,vlist,astatus)
             CALL addGlobAttFlt('nh_testcase_nml:bruntvais_u_mwbr_const',&
                  &             bruntvais_u_mwbr_const,vlist,astatus)
          !
          ELSEIF ( TRIM(nh_test_name) == 'bell') THEN
             CALL addGlobAttFlt('nh_testcase_nml:mount_height',          &
                  &             mount_height,vlist,astatus)
             CALL addGlobAttFlt('nh_testcase_nml:torus_domain_length',   &
                  &             torus_domain_length,vlist,astatus)
             CALL addGlobAttFlt('nh_testcase_nml:nh_brunt_vais',         &
                  &             nh_brunt_vais,vlist,astatus)
             CALL addGlobAttFlt('nh_testcase_nml:nh_u0',nh_u0,vlist,astatus)
             CALL addGlobAttFlt('nh_testcase_nml:nh_t0',nh_t0,vlist,astatus)
          !
          ELSEIF ( TRIM(nh_test_name) == 'PA') THEN
             CALL addGlobAttFlt('nh_testcase_nml:rotate_axis_deg',       &
                  &             rotate_axis_deg,vlist,astatus)
             !CALL addGlobAttTxtFromLog('nh_testcase_nml:linit_tracer_fv',&
             !    &                     linit_tracer_fv,vlist,astatus)
          !
          ELSEIF ( TRIM(nh_test_name) == 'HS_nh') THEN
             CALL addGlobAttTxtFromLog('nh_testcase_nml:lhs_nh_vn_ptb',  &
                 &                     lhs_nh_vn_ptb,vlist,astatus)
             CALL addGlobAttFlt('nh_testcase_nml:hs_nh_vn_ptb_scale',    &
                 &                     hs_nh_vn_ptb_scale,vlist,astatus)
          !
          ELSEIF ( TRIM(nh_test_name) == 'jabw' .OR. TRIM(nh_test_name) == 'mrw') THEN
             CALL addGlobAttFlt('nh_testcase_nml:qv_max',qv_max,vlist,astatus)
             !CALL addGlobAttFlt('nh_testcase_nml:rh_at_1000hpa',         &
             !    &              rh_at_1000hpa,vlist,astatus)
          !
          ELSEIF ( TRIM(nh_test_name) == 'APE_nh') THEN
             CALL addGlobAttTxt('nh_testcase_nml:ape_sst_case',          &
                 &              TRIM(ape_sst_case),vlist,astatus)
          !
          END IF
       END IF
     END IF
  END SUBROUTINE addAtmAtts


  !-------------------------------------------------------------------------------------------------
  SUBROUTINE addOceAtts(vlist,istatus)
    INTEGER, INTENT(IN) :: vlist
    INTEGER             :: istatus
  END SUBROUTINE addOceAtts


  !-------------------------------------------------------------------------------------------------
  !> Interpolate variable onto associated lon-lat grid variable
  SUBROUTINE postprocess_lonlat(ivar, itype, var_3d, k_jg, nlev)

    INTEGER,               INTENT(IN) :: ivar           ! source variable ID
    INTEGER,               INTENT(IN) :: itype          ! variable type (cell/edge/vertex based)
    TYPE(t_collected_var_ptr), INTENT(IN) :: var_3d     ! source variable (nproma,nlev,nblks)
    INTEGER,               INTENT(IN) :: k_jg           ! patch index
    INTEGER,               INTENT(IN) :: nlev           ! number of levels
    ! local variables
    CHARACTER(*), PARAMETER :: routine = TRIM("mo_io_vlist:postprocess_lonlat")
    INTEGER                 :: ierrstat, jk, n_tot
    REAL(wp), ALLOCATABLE   :: streamvar(:,:), tmp_var(:,:,:)
    INTEGER                 :: nblks_lonlat, npromz_lonlat
    TYPE (t_lon_lat_grid), POINTER :: grid ! lon-lat grid

    ! ------------------------------------------------

    ! perform RBF interpolation of field, depending on the grid on
    ! which the variable is defined.
    IF (dbg_level > 1) THEN
      WRITE(message_text,*) "PE #", p_pe, ": processing patch ", k_jg
      CALL message(routine, message_text)
    END IF

    ! "streamvar": temporary variable
    grid => lonlat_intp_config(k_jg)%lonlat_grid
    n_tot = grid%total_dim
    ALLOCATE(streamvar(n_tot,nlev), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) THEN
      CALL finish (routine, 'allocation failed')
    ENDIF

    ! for cell-based variables: interpolate gradients (finite differences) and reconstruct
    IF (itype == GATHER_C) THEN
      CALL rbf_interpol_lonlat( k_jg, p_patch(k_jg), p_int_state_lonlat(k_jg), &
        &                       var_3d%ptr(:,:,:), var_lonlat(k_jg)%var_3d(:,:,:), nlev )

    ELSEIF (itype == GATHER_V) THEN

      ! vertex-based variables (vorticity "VOR") are treated in a
      ! special way: They are first interpolated onto the cell centers
      ! and afterwards treated as scalar variables of type "GATHER_C"

      ! Note: The whole process is overly expensive! For the case that
      ! there are other vertex-based variables for output, this
      ! routine must be optimized!

      ! this requires a temporary variable:
      ALLOCATE(tmp_var(nproma, nlev, p_patch(k_jg)%nblks_c), STAT=ierrstat)
      IF (ierrstat /= SUCCESS) THEN
        CALL finish (routine, 'allocation failed')
      ENDIF
      CALL verts2cells_scalar( var_3d%ptr(:,:,:), p_patch(k_jg), &
        &                      p_int_state(k_jg)%verts_aw_cells, &
        &                      tmp_var(:,:,:), 1, nlev )
      CALL rbf_interpol_lonlat( k_jg, p_patch(k_jg), p_int_state_lonlat(k_jg), &
        &                       tmp_var(:,:,:), var_lonlat(k_jg)%var_3d(:,:,:), nlev )
      ! clean up:
      DEALLOCATE(tmp_var, STAT=ierrstat)
      IF (ierrstat /= SUCCESS) THEN
        CALL finish (routine, 'deallocation failed')
      ENDIF

    ELSEIF (itype == GATHER_E) THEN
      
      ! for edge-based variables: simple interpolation
      nblks_lonlat  = grid%nblks
      npromz_lonlat = grid%npromz
      CALL rbf_vec_interpol_lonlat( var_3d%ptr(:,:,:), p_patch(k_jg),     &
        &                           p_int_state_lonlat(k_jg),         &
        &                           var_lonlat(k_jg)%var_3d(:,:,:),   &
        &                           var_lonlat(k_jg)%var_3d_2(:,:,:), &
        &                           nblks_lonlat, npromz_lonlat,      &
        &                           1, nlev )
      
      ! reshape lon-lat variable (nproma,nlev,nblks) -> (ntot,nlev)       
      DO jk = 1, nlev
        streamvar(:,jk) = RESHAPE(var_lonlat(k_jg)%var_3d_2(:,jk,:), (/n_tot/))
      ENDDO
      
      ! write Y-component of the interpolated variable
      IF (dbg_level > 1) THEN
        WRITE(message_text,*) "PE #", p_pe, ": write  interpolated 2D variable"
        CALL message(routine, message_text)
      END IF
      CALL streamWriteVar(streamID(k_jg), lonlat_varids(ivar,k_jg,2), streamvar, 0)
        
    END IF

    ! reshape lon-lat variable (nproma,nlev,nblks) -> (ntot,nlev)       
    DO jk = 1, nlev
      streamvar(:,jk) = RESHAPE(var_lonlat(k_jg)%var_3d(:,jk,:), (/n_tot/))
    ENDDO

    ! actually write the interpolated variable (scalar/X-component)
    IF (dbg_level > 1) THEN
      WRITE(message_text,*) "PE #", p_pe, ": write  interpolated 2D variable"
      CALL message(routine, message_text)
    END IF
    CALL streamWriteVar(streamID(k_jg), lonlat_varids(ivar,k_jg,1), streamvar, 0)

    ! clean up:
    DEALLOCATE(streamvar, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) THEN
      CALL finish (routine, 'deallocation failed')
    ENDIF

  END SUBROUTINE postprocess_lonlat

END MODULE mo_io_vlist
