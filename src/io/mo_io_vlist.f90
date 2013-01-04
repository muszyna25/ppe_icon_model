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
!!  Modification by daniel Reinert, DWD (2012-03-22)
!! - some IO-routines, which might be of future use, have been moved to
!!   mo_io_util, since this module will be removed at some time.
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
  USE mo_ocean_nml,             ONLY: n_zlev, dzlev_m, iforc_oce, no_tracer,      &
    &                                 temperature_relaxation, i_sea_ice,          &
    &                                 irelax_2d_S !, i_apply_bulk
  USE mo_dynamics_config,       ONLY: iequations,lshallow_water,                  &
    &                                 idiv_method, divavg_cntrwgt,                &
    &                                 nold, nnow, nnow_rcf, lcoriolis
  USE mo_ha_dyn_config,         ONLY: ha_dyn_config
  USE mo_diffusion_config,      ONLY: diffusion_config
  USE mo_io_config,             ONLY: lwrite_omega, lwrite_pres, lwrite_z3,       &
    &                                 lwrite_vorticity, lwrite_divergence,        &
    &                                 lwrite_tend_phy, lwrite_radiation,          &
    &                                 lwrite_precip, lwrite_cloud, lwrite_tracer, &
    &                                 lwrite_tke,  lwrite_surface,                &
    &                                 lwrite_dblprec, lwrite_oce_timestepping,    &
    &                                 lwrite_extra, lwrite_decomposition,         &
    &                                 inextra_2d,inextra_3d,        &
    &                                 out_filetype, out_expname,                  &
    &                                 dt_data, dt_file, lkeep_in_sync,            &
    &                                 lflux_avg
  USE mo_io_util,               ONLY: gather_array1, gather_array2, outvar_desc,  &
    &                                 t_outvar_desc, GATHER_C, GATHER_E, GATHER_V,&
    &                                 max_outvars, max_gridlevs,   &
    &                                 t_collected_var_ptr, num_output_vars
  USE mo_nh_pzlev_config,       ONLY: nh_pzlev_config
  USE mo_parallel_config,       ONLY: nproma
  USE mo_extpar_config,         ONLY: itopo
  USE mo_run_config,            ONLY: num_lev, num_levp1, iforcing, lforcing,     &
    &                                 ntracer, ltransport, nsteps, dtime,         &
    &                                 dtime_adv, ldynamics, ltestcase,            &
    &                                 lvert_nest, msg_level, iqv, iqc, iqi,       &
    &                                 nqtendphy
  USE mo_grid_config,           ONLY: global_cell_type
  USE mo_echam_phy_config
  USE mo_atm_phy_nwp_config,    ONLY: atm_phy_nwp_config
  USE mo_advection_config,      ONLY: advection_config
  USE mo_echam_conv_config,     ONLY: echam_conv_config
  USE mo_lnd_nwp_config,        ONLY: ntiles_total, nlev_snow
! USE mo_gw_hines_nml,          ONLY: lheatcal, emiss_lev, rmscon, kstar, m_min
  USE mo_vertical_coord_table,  ONLY: vct
  USE mo_grid_config,           ONLY: start_lev, nroot, n_dom, lfeedback, lplane, &
    &                                 n_dom_start
  USE mo_model_domain,          ONLY: t_patch, p_patch
  USE mo_physical_constants,    ONLY: grav
  USE mo_mpi,                   ONLY: my_process_is_stdio, p_recv, p_send, &
    &                                 num_work_procs, get_my_mpi_all_id
  USE mo_icoham_dyn_types,      ONLY: t_hydro_atm_prog, t_hydro_atm_diag
  USE mo_nonhydro_types,        ONLY: t_nh_prog, t_nh_diag
  USE mo_opt_diagnostics,       ONLY: t_nh_diag_pz
  USE mo_oce_state,             ONLY: t_hydro_ocean_state, t_hydro_ocean_prog,       &
       &                              t_hydro_ocean_diag, t_hydro_ocean_base,        &
       &                              t_hydro_ocean_aux,                             &
       &                              v_base, set_zlev, v_ocean_state
!  USE mo_oce_forcing,           ONLY: t_sfc_flx, v_sfc_flx
  ! #
  USE mo_sea_ice_types,                ONLY: t_sfc_flx, v_sfc_flx, t_sea_ice, v_sea_ice
  USE mo_icoham_dyn_memory,     ONLY: p_hydro_state
  USE mo_nonhydro_state,        ONLY: p_nh_state
  USE mo_nwp_lnd_types,         ONLY: t_lnd_prog, t_lnd_diag
  USE mo_nwp_lnd_state,         ONLY: p_lnd_state
  USE mo_nwp_phy_state,         ONLY: prm_diag, prm_nwp_tend !, t_nwp_phy_diag
  USE mo_ext_data_state,        ONLY: ext_data
  USE mo_echam_phy_memory,      ONLY: prm_field, prm_tend
  USE mo_icoham_sfc_indices,    ONLY: nsfc_type, iice
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
    &                                 hs_nh_vn_ptb_scale, qv_max,                 &
    &                                 ape_sst_case, ape_sst_val
    !&                                rh_at_1000hpa,linit_tracer_fv
  USE mo_nh_mrw_exp,            ONLY: u0_mrw, mount_height_mrw,                   &
    &                                 mount_lonctr_mrw_deg, mount_latctr_mrw_deg, &
    &                                 p_int_mwbr_const, temp_i_mwbr_const,        &
    &                                 bruntvais_u_mwbr_const     
    !&                                mount_half_width
  USE mo_math_constants,        ONLY: pi
  USE mo_impl_constants,        ONLY: SUCCESS
  USE mo_intp,                  ONLY: verts2cells_scalar
  USE mo_intp_data_strc,        ONLY: p_int_state
  USE mo_mpi,                   ONLY: p_pe
  USE mo_util_string,           ONLY: string_contains_word, toupper
  USE mo_oce_physics,           ONLY: t_ho_params, v_params
  USE mo_linked_list,           ONLY: t_list_element
  USE mo_var_list,              ONLY: nvar_lists, var_lists
  IMPLICIT NONE

  PRIVATE

  INCLUDE 'cdi.inc'
  INCLUDE 'netcdf.inc'

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

  INTEGER, PARAMETER :: max_name_len = 100

  !> level of output verbosity
  INTEGER, PARAMETER  :: dbg_level = 0

  PUBLIC :: setup_vlist, destruct_vlist,                            &
    &       open_output_vlist, close_output_vlist,                  &
    &       write_vlist, get_outvar_ptr_ha, get_outvar_ptr_nh,      &
    &       get_outvar_ptr_oce,                                     &
    &       vlist_write_var, vlist_set_date_time, vlist_start_step, &
    &       de_reshape1, de_reshape2,                               &
    &       addGlobAtts, addAtmAtts, addOceAtts, translate_vars

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
    &  zaxisID_generic_ice,                                    &
    &  zaxisID_pres, zaxisID_hgt

  INTEGER, DIMENSION(max_outvars,max_gridlevs) ::  &
    &  varids
  ! current number of output variables, gets updated by addVar()
  INTEGER, PRIVATE :: num_varids(max_dom)

  INTEGER, SAVE :: iostep = 0


  PUBLIC :: t_outvar_desc
  INTEGER :: klev

  REAL(wp),POINTER,DIMENSION(:,:) :: cell_owner => NULL()
CONTAINS

  !-------------------------------------------------------------------------
  !BOC
  SUBROUTINE setup_vlist(grid_filename, k_jg, l_do_io)

    CHARACTER(len=*), INTENT(in) :: grid_filename
    INTEGER, INTENT(in) :: k_jg
    LOGICAL, INTENT(in) :: l_do_io

    INTEGER :: ncid, dimid, varid
    INTEGER :: i_nc, i_ne, i_nv
    INTEGER :: i_ncb, i_neb, i_nvb
    INTEGER :: lnlen, ulen, nzlevp1

    INTEGER :: nlev, nlevp1
    INTEGER :: znlev_soil

    REAL(wp), ALLOCATABLE :: clon(:), clat(:), clonv(:), clatv(:)
    REAL(wp), ALLOCATABLE :: elon(:), elat(:), elonv(:), elatv(:)
    REAL(wp), ALLOCATABLE :: vlon(:), vlat(:), vlonv(:), vlatv(:)

    REAL(wp), ALLOCATABLE :: levels(:), levels_i(:), levels_m(:)

    CHARACTER(len=21) :: name
    CHARACTER(len=12) :: qname
    CHARACTER(len=10) :: dbgname
    CHARACTER(len=3)  :: cjt
    CHARACTER(len=4)  :: sufix
    CHARACTER(len=3)  :: prefix
    CHARACTER(len=8)  :: meaning
    CHARACTER(len=10) :: varunits
    CHARACTER(LEN=1)  :: ctracer
    CHARACTER(len=MAX_CHAR_LENGTH) :: & !< list of tracers to initialize
      &  ctracer_list
    CHARACTER(LEN=1)  :: anextra ! number of debug fields
    CHARACTER(len=NF_MAX_NAME) :: long_name, units
    INTEGER :: i, jt, itracer
    INTEGER :: ivar
    INTEGER :: gridid, zaxisid
    INTEGER :: elemid, elemid2,tableid

    !CHARACTER(len=NF_MAX_NAME) :: att_txt
    INTEGER                    :: astatus

    REAL(wp) :: pi_180

    ! ocean tracers
    CHARACTER(len=max_char_length) :: oce_tracer_names(max_ntracer),&
    &                                 oce_tracer_units(max_ntracer),&
    &                                 oce_tracer_longnames(max_ntracer)
    INTEGER                    :: oce_tracer_codes(max_ntracer)
    INTEGER                    :: oce_trace_counter
    INTEGER, PARAMETER         :: oce_max_tracer = 2

    CHARACTER(len=max_char_length) :: msg
    CHARACTER(len=max_char_length), PARAMETER :: routine = 'mo_io_vlist:setup_vlist'

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
    IF(l_do_io) THEN
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
    ELSE
      i_nc = p_patch(k_jg)%n_patch_cells_g
      i_ne = p_patch(k_jg)%n_patch_edges_g
      i_nv = p_patch(k_jg)%n_patch_verts_g
    ENDIF

    !
    i_ncb = global_cell_type*i_nc
    i_neb = 4*i_ne
    i_nvb = (9-global_cell_type)*i_nv
    !
    !-------------------------------------------------------------------------
    ! For parallel runs, setting up the vlist serves two purposes:
    ! - on processes which actually do I/O, the vlist must be set up
    !   for doing the CDI output later on
    ! - on all other processes, it is just used to gather the settings
    !   in outvar_desc at the end of this routine, i.e. all the calls
    !   below just serve for setting up a list of variables to be output
    !
    ! Since some calls during vlist setup (gridDefX/Yvals, gridDefX/Ybounds)
    ! reserve a considerable amount of memory and need a considerable
    ! amount of I/O, they are left away on tasks which actually don't do I/O.
    !-------------------------------------------------------------------------
    ! cell grid
    !
    gridCellID(k_jg) = gridCreate(GRID_UNSTRUCTURED, i_nc)
    CALL gridDefNvertex(gridCellID(k_jg), global_cell_type)
    !
    IF(l_do_io) THEN
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
    ENDIF
    !
    !-------------------------------------------------------------------------
    ! edge grid
    !
    gridEdgeID(k_jg) = gridCreate(GRID_UNSTRUCTURED, i_ne)
    CALL gridDefNvertex(gridEdgeID(k_jg), 4)
    !
    IF(l_do_io) THEN
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
    ENDIF
    !
    !-------------------------------------------------------------------------
    ! vertex grid
    !
    gridVertexID(k_jg) = gridCreate(GRID_UNSTRUCTURED, i_nv)
    CALL gridDefNvertex(gridVertexID(k_jg), 9-global_cell_type)
    !
    IF(l_do_io) THEN
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

    !-------------------------------------------------------------------------

    ! Close NetCDF file, it is not needed any more
      CALL nf(nf_close(ncid))
    ENDIF
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
      znlev_soil = SIZE(zml_soil)

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
      zaxisID_depth_below_land_p1(k_jg) = zaxisCreate(ZAXIS_DEPTH_BELOW_LAND, znlev_soil+1)
      ALLOCATE(levels(znlev_soil+1))
      levels(1) = 0._wp
      DO i = 1, znlev_soil
      levels(i+1) = zml_soil(i)*100._wp
      END DO
      CALL zaxisDefLevels(zaxisID_depth_below_land_p1(k_jg), levels)
      DEALLOCATE(levels)


      ! Define axes for soil model
      !
      zaxisID_depth_below_land(k_jg) = zaxisCreate(ZAXIS_DEPTH_BELOW_LAND, znlev_soil)
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

    ELSE ! oce
      zaxisIDdepth_m(k_jg)  = zaxisCreate(ZAXIS_DEPTH_BELOW_SEA, n_zlev)
      nzlevp1 = n_zlev + 1
      zaxisID_halfdepth(k_jg)  = zaxisCreate(ZAXIS_DEPTH_BELOW_SEA, nzlevp1)

      ALLOCATE(levels_i(nzlevp1))
      ALLOCATE(levels_m(n_zlev))
      CALL set_zlev(levels_i, levels_m)
      CALL zaxisDefLevels(zaxisIDdepth_m(k_jg)   , levels_m)
      CALL zaxisDefLevels(zaxisID_halfdepth(k_jg), levels_i)
      DEALLOCATE(levels_i)
      DEALLOCATE(levels_m)
    ENDIF
    ! The ice axis must exist for atmosphere and ocean
    zaxisID_generic_ice(k_jg) = zaxisCreate(ZAXIS_GENERIC, 1) !TOOE v_ice%kice
    !
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
    ! register variables
    varids(:,k_jg)   = 0
    ! initialize total number of varids for domain jg
    num_varids(k_jg) = 0
    ! atm
    IF (iequations/=ihs_ocean) THEN

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
        &                   'Pa', 54, 128,&
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
      IF(iequations == 1 .OR. iequations == 2 ) THEN
        CALL addVar(TimeVar('PHIS',&
        &                   'surface geopotential (orography)',&
        &                   'm**2/s**2',&
        &                   129, 128,&
        &                   vlistID(k_jg),gridCellID(k_jg),zaxisID_surface(k_jg)),&
        &           k_jg)
      END IF
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
            IF (jt.EQ.1) THEN !for water vapour
             elemid=80; tableid=128   
             elemid2=81
            ELSEIF (jt.EQ.2) THEN !for cloud water
             elemid=82; tableid=128
             elemid2=83 
            ELSEIF (jt.EQ.3) THEN !for cloud ice
             elemid=84; tableid=128
             elemid2=85 
            ELSE !other tracers
             elemid=86; tableid=128 !default coding
             elemid2=87
            END IF
            WRITE(name,'(A2,A1)') "TQ", ctracer
            WRITE(long_name,'(A34,A1)') "vertically integrated grid-scale Q",ctracer
            CALL addVar(TimeVar(TRIM(name),TRIM(long_name),&
            &                   'kg/m**2',elemid,tableid,&
            &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
            &           k_jg)
            WRITE(name,'(A2,A1,A4)') "TQ", ctracer,"_avg"
            WRITE(long_name,'(A42,A1)') "average vertically integrated grid-scale Q",ctracer
            CALL addVar(TimeVar(TRIM(name),TRIM(long_name),&
            &                   'kg/m**2',elemid2,tableid,&
            &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
            &           k_jg)

          END IF !lwrite_tracer(jt)
         END DO
         IF ( irad_o3 == 4 .OR. irad_o3 == 6 .OR. irad_o3 == 7 ) THEN     !output for O3
                        CALL addVar(TimeVar('O3','O3',&
            &                   'kg/kg',203,128,&
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
              &                   'g/g', 203, 128,&
              &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_hybrid(k_jg)),&
              &           k_jg)

            !KF this variable should only go into the output as an averaged one 
            WRITE(name,'(A14)') "dlwfsfc_dT_avg"
            WRITE(long_name,'(A8,A27)') "averaged", " longwave surface net flux T-tend"
            CALL addVar(TimeVar(TRIM(name),&
              &                   TRIM(long_name),&
              &                   'W/m**2/K', 112, 128,&
              &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),k_jg)


            IF (lflux_avg ) THEN
              sufix = "_avg"
              meaning = "averaged"
            ELSE 
              sufix = ""
              meaning = "instant."     
            END IF

            WRITE(name,'(A8,A4)') "swflxsfc", sufix
            WRITE(long_name,'(A8,A27)') meaning, " shortwave surface net flux"
            CALL addVar(TimeVar(TRIM(name),&
              &                   TRIM(long_name),&
              &                   'W/m**2', 111, 128,&
              &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),k_jg)

            WRITE(name,'(A8,A4)') "lwflxsfc", sufix
            WRITE(long_name,'(A8,A27)') meaning, " longwave  surface net flux"
            CALL addVar(TimeVar(TRIM(name),&
              &                   TRIM(long_name),&
              &                   'W/m**2', 112, 128,&
              &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)), k_jg)


            WRITE(name,'(A8,A4)') "swflxtoa", sufix
            WRITE(long_name,'(A8,A23)') meaning, " shortwave toa net flux"
            CALL addVar(TimeVar(TRIM(name),&
              &                   TRIM(long_name),&
              &                   'W/m**2', 113, 128,&
              &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)), k_jg)

            WRITE(name,'(A8,A4)') "lwflxtoa", sufix
            WRITE(long_name,'(A8,A23)') meaning, " longwave  toa net flux"
            CALL addVar(TimeVar(TRIM(name),&
              &                   TRIM(long_name),&
              &                   'W/m**2', 114, 128,&
              &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)), k_jg)
          ENDIF
        END SELECT

        SELECT CASE (iforcing)
        CASE (inwp)
          IF (lflux_avg ) THEN
            prefix = "A"
            meaning = "averaged"
            varunits = "W/m**2"
          ELSE 
            prefix = "ACC"
            meaning = "accumul." 
            varunits = "J/m**2"    
          END IF


          WRITE(name,'(A,A5)') TRIM(prefix),"SOB_S"
          WRITE(long_name,'(A8,A27)') meaning, " shortwave surface net flux"
          CALL addVar(TimeVar(TRIM(name),&
            &                   TRIM(long_name),&
            &                   TRIM(varunits), 111, 128,&
            &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
            &           k_jg)

          WRITE(name,'(A,A5)') TRIM(prefix),"THB_S"
          WRITE(long_name,'(A8,A27)') meaning, " longwave  surface net flux"
          CALL addVar(TimeVar(TRIM(name),&
            &                   TRIM(long_name),&
            &                   TRIM(varunits), 112, 128,&
            &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
            &           k_jg)
          WRITE(name,'(A,A5)') TRIM(prefix),"SOB_T"
          WRITE(long_name,'(A8,A23)') meaning, " shortwave toa net flux"
          CALL addVar(TimeVar(TRIM(name),&
            &                   TRIM(long_name),&
            &                   TRIM(varunits), 113, 128,&
            &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
            &           k_jg)
          WRITE(name,'(A,A5)') TRIM(prefix),"THB_T"
          WRITE(long_name,'(A8,A23)') meaning, " longwave  toa net flux"
          CALL addVar(TimeVar(TRIM(name),&
            &                   TRIM(long_name),&
            &                   TRIM(varunits), 114, 128,&
            &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
            &           k_jg)

        END SELECT

        SELECT CASE (iforcing)
        CASE (inwp)
          CALL addVar(TimeVar('SOB_S',&
            &                   'shortwave surface net flux',&
            &                   'W/m**2', 176, 128,&
            &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
            &           k_jg)
          CALL addVar(TimeVar('THB_S',&
            &                   'longwave surface net flux',&
            &                   'W/m**2', 177, 128,&
            &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
            &           k_jg)
          CALL addVar(TimeVar('SOB_T',&
            &                   'shortwave toa net flux',&
            &                   'W/m**2', 178, 128,&
            &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
            &           k_jg)
          CALL addVar(TimeVar('THB_T',&
            &                   'longwave toa net flux',&
            &                   'W/m**2', 179, 128,&
            &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
            &           k_jg)
          CALL addVar(TimeVar('ALB_RAD',&
            &                   'diffuse solar surface albedo',&
            &                   ' ', 243, 128,&
            &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
            &           k_jg)
        END SELECT

      END IF !lwrite_radiation

      ! surface precitation rates
      IF(lwrite_precip) THEN

        SELECT CASE (iforcing)
        CASE (inwp)
          CALL addVar(TimeVar('PRR_GSP',&
          &                   'grid-scale rain precipitation rate',&
          &                   'kg/s/m**2', 104, 128,&
          &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
          &           k_jg)
          CALL addVar(TimeVar('PRS_GSP',&
          &                   'grid-scale snow precipitation rate',&
          &                   'kg/s/m**2', 105, 128,&
          &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
          &           k_jg)
          CALL addVar(TimeVar('RAIN_GSP',&
          &                   'grid-scale accumulated surface rain',&
          &                   'kg/m**2', 106, 128,&
          &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
          &           k_jg)
          CALL addVar(TimeVar('SNOW_GSP',&
          &                   'grid-scale accumulated surface snow',&
          &                   'kg/m**2', 107, 128,&
          &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
          &           k_jg)
          CALL addVar(TimeVar('RAIN_CON',&
          &                   'convective accumulated surface rain',&
          &                   'kg/m**2', 108, 128,&
          &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
          &           k_jg)
          CALL addVar(TimeVar('SNOW_CON',&
          &                   'convective accumulated surface snow',&
          &                   'kg/m**2', 109, 128,&
          &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
          &           k_jg)
          CALL addVar(TimeVar('TOT_PREC',&
          &                   'grid-scale + convective accumulated surface total precipitation',&
          &                   'kg/m**2', 228, 128,&
          &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
          &           k_jg)
          CALL addVar(TimeVar('TOT_PREC_RATE_avg',&
          &                   'average grid-scale + convective surface total precipitation rate',&
          &                   'kg/m**2/s', 110, 128,&
          &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
          &           k_jg)
          CALL addVar(TimeVar('CON_PREC_RATE_avg',&
          &                   'average convective surface precipitation rate',&
          &                   'kg/m**2/s', 143, 128,&
          &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
          &           k_jg)
          CALL addVar(TimeVar('GSP_PREC_RATE_avg',&
          &                   'average grid-scale surface precipitation rate',&
          &                   'kg/m**2/s', 50, 128,&
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

!!$ TR JSBACH output for testing
      IF (echam_phy_config%ljsbach) THEN

      CALL addVar(TimeVar('surface_temperature',&
      &                   'temperature of land surface',&
      &                   'K', 21, 128,&
      &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
      &           k_jg)
      CALL addVar(TimeVar('surface_temperature_old',&
      &                   'temperature of land surface previous time step',&
      &                   'K', 22, 128,&
      &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
      &           k_jg)

      CALL addVar(TimeVar('surface_temperature_rad',&
      &                   'radiative temperature',&
      &                   'K', 23, 128,&
      &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
      &           k_jg)

      CALL addVar(TimeVar('surface_temperature_eff',&
      &                   'effective temperature',&
      &                   'K', 24, 128,&
      &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
      &           k_jg)

      CALL addVar(TimeVar('c_soil_temperature1',&
      &                   'soil temperature parameter c layer 1',&
      &                   ' ', 25, 128,&
      &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
      &           k_jg)
      CALL addVar(TimeVar('c_soil_temperature2',&
      &                   'soil temperature parameter c layer 2',&
      &                   ' ', 26, 128,&
      &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
      &           k_jg)
      CALL addVar(TimeVar('c_soil_temperature3',&
      &                   'soil temperature parameter c layer 3',&
      &                   ' ', 27, 128,&
      &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
      &           k_jg)
      CALL addVar(TimeVar('c_soil_temperature4',&
      &                   'soil temperature parameter c layer 4',&
      &                   ' ', 28, 128,&
      &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
      &           k_jg)
      CALL addVar(TimeVar('c_soil_temperature5',&
      &                   'soil temperature parameter c layer 5',&
      &                   ' ', 29, 128,&
      &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
      &           k_jg)
      CALL addVar(TimeVar('d_soil_temperature1',&
      &                   'soil temperature parameter d layer 1',&
      &                   ' ', 30, 128,&
      &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
      &           k_jg)
      CALL addVar(TimeVar('d_soil_temperature2',&
      &                   'soil temperature parameter d layer 2',&
      &                   ' ', 31, 128,&
      &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
      &           k_jg)
      CALL addVar(TimeVar('d_soil_temperature3',&
      &                   'soil temperature parameter d layer 3',&
      &                   ' ', 32, 128,&
      &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
      &           k_jg)
      CALL addVar(TimeVar('d_soil_temperature4',&
      &                   'soil temperature parameter d layer 4',&
      &                   ' ', 33, 128,&
      &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
      &           k_jg)
      CALL addVar(TimeVar('d_soil_temperature5',&
      &                   'soil temperature parameter d layer 5',&
      &                   ' ', 34, 128,&
      &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
      &           k_jg)
      CALL addVar(TimeVar('soil_temperature1',&
      &                   'soil temperature layer 1',&
      &                   ' ', 35, 128,&
      &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
      &           k_jg)
      CALL addVar(TimeVar('soil_temperature2',&
      &                   'soil temperature layer 2',&
      &                   ' ', 36, 128,&
      &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
      &           k_jg)
      CALL addVar(TimeVar('soil_temperature3',&
      &                   'soil temperature layer 3',&
      &                   ' ', 37, 128,&
      &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
      &           k_jg)
      CALL addVar(TimeVar('soil_temperature4',&
      &                   'soil temperature layer 4',&
      &                   ' ', 38, 128,&
      &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
      &           k_jg)
      CALL addVar(TimeVar('soil_temperature5',&
      &                   'soil temperature layer 5',&
      &                   ' ', 39, 128,&
      &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
      &           k_jg)
      CALL addVar(TimeVar('heat_capacity',&
      &                   'heat capacity of soil layer 1',&
      &                   ' ', 40, 128,&
      &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
      &           k_jg)
      CALL addVar(TimeVar('ground_heat_flux',&
      &                   'ground heat flux',&
      &                   'Wm-2', 41, 128,&
      &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
      &           k_jg)
      CALL addVar(TimeVar('swnet',&
      &                   'swnet',&
      &                   'Wm-2', 42, 128,&
      &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
      &           k_jg)
      CALL addVar(TimeVar('time_steps_soil',&
      &                   'time_steps_soil',&
      &                   '', 43, 128,&
      &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
      &           k_jg)
     CALL addVar(TimeVar('moisture1',&
      &                   'soil moisture layer 1',&
      &                   'm', 44, 128,&
      &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
      &           k_jg)
     CALL addVar(TimeVar('moisture2',&
      &                   'soil moisture layer 2',&
      &                   'm', 45, 128,&
      &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
      &           k_jg)
     CALL addVar(TimeVar('moisture3',&
      &                   'soil moisture layer 3',&
      &                   'm', 46, 128,&
      &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
      &           k_jg)
     CALL addVar(TimeVar('moisture4',&
      &                   'soil moisture layer 4',&
      &                   'm', 47, 128,&
      &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
      &           k_jg)
     CALL addVar(TimeVar('moisture5',&
      &                   'soil moisture layer 5',&
      &                   'm', 48, 128,&
      &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
      &           k_jg)
     CALL addVar(TimeVar('moisture_all',&
      &                   'sum of soil moisture in all layers',&
      &                   'm', 49, 128,&
      &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
      &           k_jg)
     CALL addVar(TimeVar('skin_reservoir',&
      &                   'skin reservoir',&
      &                   'm', 51, 128,&
      &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
      &           k_jg)
     CALL addVar(TimeVar('snow',&
      &                   'snow depth',&
      &                   'm', 52, 128,&
      &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
      &           k_jg)
     CALL addVar(TimeVar('albvisdir',&
      &                   'albedo visible',&
      &                   '', 53, 128,&
      &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
      &           k_jg)
     CALL addVar(TimeVar('albnirdir',&
      &                   'albedo NIR',&
      &                   '', 54, 128,&
      &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
      &           k_jg)
      END IF ! ljsbach

        END SELECT !iforcing
      ENDIF !lwrite_precip

      ! cloud
      IF(lwrite_cloud ) THEN
        SELECT CASE (iforcing)
        CASE (inwp)
          !--- total water vapor---
          CALL addVar(TimeVar('QV',&
          &                   'total water vapor',&
          &                   '(0-1)', 133, 128,&
          &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_hybrid(k_jg)),&
          &           k_jg)
          !--- total cloud water ---
          CALL addVar(TimeVar('QC',&
          &                   'total cloud water',&
          &                   '(0-1)', 246, 128,&
          &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_hybrid(k_jg)),&
          &           k_jg)
          !--- total cloud ice ---
          CALL addVar(TimeVar('QI',&
          &                   'total cloud ice',&
          &                   '(0-1)', 247, 128,&
          &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_hybrid(k_jg)),&
          &           k_jg)
          !--- total cloud cover ---
          CALL addVar(TimeVar('CC',&
          &                   'total cloud cover',&
          &                   '(0-1)', 248, 128,&
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
          &                   'kk/m**2', 78, 128,&
          &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
          &           k_jg)
          !--- vertically integrated total cloud ice ---
          CALL addVar(TimeVar('TQI',&
          &                   'vertically integrated total cloud ice',&
          &                   'kg/m**2', 79, 128,&
          &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
          &           k_jg)
          !--- cloud cover assuming Maximum-Random overlap ---
          CALL addVar(TimeVar('TCC',&
          &                   'cloud cover assuming Maximum-Random overlap',&
          &                   '(0-1)', 164, 128,&
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
          !--- turbulent kinetic energy---
          CALL addVar(TimeVar('TKE',&
          &                   'turbulent kinetic energy',&
          &                   'm^2/s^2', 152, 201,&
          &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_hybrid_half(k_jg)),&
          &           k_jg)
!DR to be implemented
!!$          !--- TKE-tendency ---
!!$          CALL addVar(TimeVar('tend_tke',&
!!$          &                   'turbulent kinetic energy tendency',&
!!$          &                   'm2/s3', 125, 999,&
!!$          &                   vlistID(k_jg),gridCellID(k_jg),zaxisID_hybrid(k_jg)),&
!!$          &           k_jg)
          !--- turbulent transfer coefficients for momentum ---
          CALL addVar(TimeVar('TCM',&
          &                   'turbulent transfer coefficients for momentum',&
          &                   '', 170, 201,&
          &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
          &           k_jg)
          !--- turbulent transfer coefficients for heat ---
          CALL addVar(TimeVar('TCH',&
          &                   'turbulent transfer coefficients for heat ',&
          &                   '', 171, 201,&
          &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
          &           k_jg)
        END SELECT
      ENDIF !lwrite_tke

      IF(lwrite_surface ) THEN
        SELECT CASE (iforcing)
        CASE (inwp)
          !--- roughness length---
          CALL addVar(TimeVar ('Z0',&
               &              'roughness length',&
          &                   'm', 173, 128,&
          &                    vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
          &                    k_jg)
          !--- Temperature at Surface---
          CALL addVar(TimeVar('T_G',&
          &                   'aggregated surface temperature',&
          &                   'K', 235, 128,&
          &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
          &                   k_jg)

        IF ( atm_phy_nwp_config(k_jg)%inwp_surface == 1 ) THEN  ! TERRA
          !--- Weighted temperature at surface---
          DO jt = 1, ntiles_total
            WRITE(cjt,'(i2)') jt
            WRITE(name,'(A,A)') "T_GT_tile_", TRIM(ADJUSTL(cjt))
            WRITE(long_name,'(A,A)') "weighted surface temperature tile ",TRIM(ADJUSTL(cjt))
            CALL addVar(TimeVar(TRIM(name),TRIM(long_name),&
            &          'K',11,201,&
            &           vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
            &           k_jg)
          ENDDO

          DO jt = 1, ntiles_total
            WRITE(cjt,'(i2)') jt
            WRITE(name,'(A,A)') "T_SNOW_tile_", TRIM(ADJUSTL(cjt))
            WRITE(long_name,'(A,A)') "temperature of the snow-surface tile ",TRIM(ADJUSTL(cjt))
            CALL addVar(TimeVar(TRIM(name),TRIM(long_name),&
            &          'K',11,201,&
            &           vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
            &           k_jg)
          ENDDO

          DO jt = 1, ntiles_total
            WRITE(cjt,'(i2)') jt
            WRITE(name,'(A,A)') "T_SNOW_MULT_tile_", TRIM(ADJUSTL(cjt))
            WRITE(long_name,'(A,A)') "temperature of the snow-surface tile ",TRIM(ADJUSTL(cjt))
            CALL addVar(TimeVar(TRIM(name),TRIM(long_name),&
            &          'K',11,201,&
            &           vlistID(k_jg), gridCellID(k_jg),zaxisID_generic_snow_p1(k_jg)),&
            &           k_jg)
          ENDDO

          DO jt = 1, ntiles_total
            WRITE(cjt,'(i2)') jt
            WRITE(name,'(A,A)') "T_S_tile_", TRIM(ADJUSTL(cjt))
            WRITE(long_name,'(A,A)') "temperature of ground surface tile ",TRIM(ADJUSTL(cjt))
            CALL addVar(TimeVar(TRIM(name),TRIM(long_name),&
            &          'K',11,201,&
            &           vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
            &           k_jg)
          ENDDO

          DO jt = 1, ntiles_total
            WRITE(cjt,'(i2)') jt
            WRITE(name,'(A,A)') "W_SNOW_tile_", TRIM(ADJUSTL(cjt))
            WRITE(long_name,'(A,A)') "water content of snow tile ",TRIM(ADJUSTL(cjt))
            CALL addVar(TimeVar(TRIM(name),TRIM(long_name),&
            &          'm H2O',11,201,&
            &           vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
            &           k_jg)
          ENDDO

          DO jt = 1, ntiles_total
            WRITE(cjt,'(i2)') jt
            WRITE(name,'(A,A)') "RHO_SNOW_tile_", TRIM(ADJUSTL(cjt))
            WRITE(long_name,'(A,A)') "snow density tile ",TRIM(ADJUSTL(cjt))
            CALL addVar(TimeVar(TRIM(name),TRIM(long_name),&
            &          'kg/m**3',11,201,&
            &           vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
            &           k_jg)
          ENDDO

          DO jt = 1, ntiles_total
            WRITE(cjt,'(i2)') jt
            WRITE(name,'(A,A)') "RHO_SNOW_MULT_tile_", TRIM(ADJUSTL(cjt))
            WRITE(long_name,'(A,A)') "snow density tile ",TRIM(ADJUSTL(cjt))
            CALL addVar(TimeVar(TRIM(name),TRIM(long_name),&
            &          'kg/m**3',11,201,&
            &           vlistID(k_jg), gridCellID(k_jg),zaxisID_generic_snow(k_jg)),&
            &           k_jg)
          ENDDO

          DO jt = 1, ntiles_total
            WRITE(cjt,'(i2)') jt
            WRITE(name,'(A,A)') "W_I_tile_", TRIM(ADJUSTL(cjt))
            WRITE(long_name,'(A,A)') "water content of interception water tile",TRIM(ADJUSTL(cjt))
            CALL addVar(TimeVar(TRIM(name),TRIM(long_name),&
            &          'm H2O',11,201,&
            &           vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
            &           k_jg)

          ENDDO

          DO jt = 1, ntiles_total
            WRITE(cjt,'(i2)') jt
            WRITE(name,'(A,A)') "T_SO_tile_", TRIM(ADJUSTL(cjt))
            WRITE(long_name,'(A,A)') "soil temperature (main level) tile",TRIM(ADJUSTL(cjt))
            CALL addVar(TimeVar(TRIM(name),TRIM(long_name),&
            &          'K',11,201,&
            &           vlistID(k_jg), gridCellID(k_jg),zaxisID_depth_below_land_p1(k_jg)),&
            &           k_jg)

          ENDDO

          DO jt = 1, ntiles_total
            WRITE(cjt,'(i2)') jt
            WRITE(name,'(A,A)') "W_SO_tile_", TRIM(ADJUSTL(cjt))
            WRITE(long_name,'(A,A)') "total water content (ice + liquid water) tile", &
            &                        TRIM(ADJUSTL(cjt))
            CALL addVar(TimeVar(TRIM(name),TRIM(long_name),&
            &          'm H2O',11,201,&
            &           vlistID(k_jg), gridCellID(k_jg),zaxisID_depth_below_land(k_jg)),&
            &           k_jg)

          ENDDO

          DO jt = 1, ntiles_total
            WRITE(cjt,'(i2)') jt
            WRITE(name,'(A,A)') "W_SO_ICE_tile_", TRIM(ADJUSTL(cjt))
            WRITE(long_name,'(A,A)') "ice content tile",TRIM(ADJUSTL(cjt))
            CALL addVar(TimeVar(TRIM(name),TRIM(long_name),&
            &          'm H2O',11,201,&
            &           vlistID(k_jg), gridCellID(k_jg),zaxisID_depth_below_land(k_jg)),&
            &           k_jg)

          ENDDO

          DO jt = 1, ntiles_total
            WRITE(cjt,'(i2)') jt
            WRITE(name,'(A,A)') "WLIQ_SNOW_tile_", TRIM(ADJUSTL(cjt))
            WRITE(long_name,'(A,A)') "liquid water content in snow tile",TRIM(ADJUSTL(cjt))
            CALL addVar(TimeVar(TRIM(name),TRIM(long_name),&
            &          'm H2O',11,201,&
            &           vlistID(k_jg), gridCellID(k_jg),zaxisID_generic_snow(k_jg)),&
            &           k_jg)

          ENDDO

          DO jt = 1, ntiles_total
            WRITE(cjt,'(i2)') jt
            WRITE(name,'(A,A)') "WTOT_SNOW_tile_", TRIM(ADJUSTL(cjt))
            WRITE(long_name,'(A,A)') "total water content in snow tile",TRIM(ADJUSTL(cjt))
            CALL addVar(TimeVar(TRIM(name),TRIM(long_name),&
            &          'm H2O',11,201,&
            &           vlistID(k_jg), gridCellID(k_jg),zaxisID_generic_snow(k_jg)),&
            &           k_jg)

          ENDDO

          DO jt = 1, ntiles_total
            WRITE(cjt,'(i2)') jt
            WRITE(name,'(A,A)') "H_SNOW_tile_", TRIM(ADJUSTL(cjt))
            WRITE(long_name,'(A,A)') "snow height", &
            & TRIM(ADJUSTL(cjt))
            CALL addVar(TimeVar(TRIM(name),TRIM(long_name),&
            &          'm',11,201,&
            &           vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
            &           k_jg)

          ENDDO

          DO jt = 1, ntiles_total
            WRITE(cjt,'(i2)') jt
            WRITE(name,'(A,A)') "FRESHSNOW_tile_", TRIM(ADJUSTL(cjt))
            WRITE(long_name,'(A,A)') "age of snow indicator (top layer)", &
            & TRIM(ADJUSTL(cjt))
            CALL addVar(TimeVar(TRIM(name),TRIM(long_name),&
            &          'm',11,201,&
            &           vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
            &           k_jg)

          ENDDO

          DO jt = 1, ntiles_total
            WRITE(cjt,'(i2)') jt
            WRITE(name,'(A,A)') "SNOWFRAC_tile_", TRIM(ADJUSTL(cjt))
            WRITE(long_name,'(A,A)') "snow-cover fraction", &
            & TRIM(ADJUSTL(cjt))
            CALL addVar(TimeVar(TRIM(name),TRIM(long_name),&
            &          'm',11,201,&
            &           vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
            &           k_jg)

          ENDDO

          DO jt = 1, ntiles_total
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
          &                   'K', 233, 128,&
          &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
          &                   k_jg)
          !--- Fluxes ....---
          IF (lflux_avg ) THEN
            prefix = "A"
            meaning = "averaged"
            varunits= "W/m**2"
          ELSE
            prefix = "ACC"
            meaning = "accumul." 
            varunits= "J/m**2"    
          END IF

          WRITE(name,'(A,A6)') TRIM(prefix),"SHFL_S"
          WRITE(long_name,'(A8,A30)') meaning, " sensible heat flux at surface"
          CALL addVar(TimeVar(TRIM(name),&
          &                   TRIM(long_name),&
          &                   TRIM(varunits), 255, 201,&  !999 == WMO here
          &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
          &                   k_jg)
          WRITE(name,'(A,A6)') TRIM(prefix),"LHFL_S"
          WRITE(long_name,'(A8,A30)') meaning, " latent   heat flux at surface"
          CALL addVar(TimeVar(TRIM(name),&
          &                   TRIM(long_name),&
          &                   TRIM(varunits), 91, 128,& !999 ==  WMO here
          &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
          &                   k_jg)
          CALL addVar(TimeVar('SHFL_S',&
          &                   'sensible heat flux at surface',&
          &                   'W/m**2', 92, 128,& 
          &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
          &                   k_jg)
          CALL addVar(TimeVar('LHFL_S',&
          &                   'latent heat flux at surface',&
          &                   'W/m**2', 147, 255,& 
          &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
          &                   k_jg)
          CALL addVar(TimeVar('EVAP_RATE_avg',&
          &                   'averaged moisture flux rate (evap. rate) at surface',&
          &                   'kg/m**2/s', 100, 128,& 
          &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
          &                   k_jg)
          CALL addVar(TimeVar('PS_s6avg',&
          &                   '6 hourly sample surface pressure average',&
          &                   'Pa', 117, 128,&
          &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
          &           k_jg)
          CALL addVar(TimeVar('T_2M',&
          &                   '2 m Temperature',&
          &                   'K', 167, 128,&
          &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
          &           k_jg)
          CALL addVar(TimeVar('T_2M_s6avg',&
          &                   '6 hourly sample 2 m Temperature average',&
          &                   'K', 120, 128,&
          &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
          &           k_jg)
          CALL addVar(TimeVar('QV_2M',&
          &                   '2 m specific humidity ',&
          &                   'Kg/kg', 115, 128,&
          &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
          &           k_jg)
          CALL addVar(TimeVar('QV_2M_s6avg',&
          &                   '6 hourly sample 2 m specific humidity average',&
          &                   'Kg/kg', 116, 128,&
          &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
          &           k_jg)
          CALL addVar(TimeVar('U_10M',&
          &                   '10 m zonal wind',&
          &                   'K', 165, 128,&
          &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
          &           k_jg)
          CALL addVar(TimeVar('V_10M',&
          &                   '10 m meridional wind',&
          &                   'K', 166, 128,&
          &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
          &           k_jg)
          CALL addVar(TimeVar('U_10M_s6avg',&
          &                   '6 hourly sample 10 m zonal wind average',&
          &                   'K', 118, 128,&
          &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
          &           k_jg)
          CALL addVar(TimeVar('V_10M_s6avg',&
          &                   '6 hourly sample 10 m meridional wind average',&
          &                   'K', 119, 128,&
          &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),&
          &           k_jg)


        CASE (iecham,ildf_echam)

          CALL addVar(TimeVar('TOTPREC_AVG',&
          &                   'average over output total preipitation flux',&
          &                   'kg/m**2/s', 110, 128,&
          &                   vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),k_jg)
          CALL addVar(TimeVar('evap_avg',                                         &
          &           'evaporation accumulated over output interval',            &
          &           'kg m-2', 182, 128,                                        &
          &           vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),k_jg)

          CALL addVar(TimeVar('lhflx_avg',                                        &
          &           'latent heat flux accumulated over output interval',       &
          &           'W m-2 s', 147, 128,                                       &
          &           vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),k_jg)

          CALL addVar(TimeVar('shflx_avg',                                        &
          &           'sensible heat flux accumulated over output interval',     &
          &           'W m-2 s', 146, 128,                                       &
          &           vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),k_jg)

          CALL addVar(TimeVar('u_stress_avg',                                     &
          &           'surface wind stress accumulated over output interval',    &
          &           'N m-2 s', 180, 128,                                       &
          &           vlistID(k_jg), gridCellID(k_jg),zaxisID_surface(k_jg)),k_jg)

          CALL addVar(TimeVar('v_stress_avg',                                     &
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

            WRITE(name,'(a,i1)') "dshflx_dT_avg_tile",jt
            CALL addVar(TimeVar(TRIM(name), TRIM(name),'W m-2 K-1',255,255,        &
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
          DO jt = 1, ntracer
            ctracer = ctracer_list(jt:jt)
            WRITE(name,'(A6,A1,A4)') "tend_q", ctracer, "_cnv"
            CALL addVar(TimeVar(TRIM(name), &
                &               TRIM(name), &
                &               '1/s',22,128,&
                &               vlistID(k_jg), gridCellID(k_jg),zaxisID_hybrid(k_jg)),&
                &       k_jg)
          END DO

          DO jt = 1, ntracer
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
          CALL addVar(TimeVar('tend_temp_drag', &
          &                   'temperature tendency due to SSO, GWD and Rayleigh friction',&
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
          &                   'm/s2', 113, 999,&
          &                   vlistID(k_jg),gridCellID(k_jg),zaxisID_hybrid(k_jg)),&
          &           k_jg)
          CALL addVar(TimeVar('tend_u_gwd',&
          &                   'zonal wind tendency due to non-orographic GWD',&
          &                   'm/s2', 114, 999,&
          &                   vlistID(k_jg),gridCellID(k_jg),zaxisID_hybrid(k_jg)),&
          &           k_jg)
          CALL addVar(TimeVar('tend_u_raylfric',&
          &                   'zonal wind tendency due to Rayleigh friction',&
          &                   'm/s2', 115, 999,&
          &                   vlistID(k_jg),gridCellID(k_jg),zaxisID_hybrid(k_jg)),&
          &           k_jg)

          ! v-wind tendency
          CALL addVar(TimeVar('tend_v_conv',&
          &                   'meridional wind tendency due to convection',&
          &                   'm/s2', 121, 999,&
          &                   vlistID(k_jg),gridCellID(k_jg),zaxisID_hybrid(k_jg)),&
          &           k_jg)
          CALL addVar(TimeVar('tend_v_turb',&
          &                   'meridional wind tendency due to turbulence',&
          &                   'm/s2', 122, 999,&
          &                   vlistID(k_jg),gridCellID(k_jg),zaxisID_hybrid(k_jg)),&
          &           k_jg)
          CALL addVar(TimeVar('tend_v_sso',&
          &                   'meridional wind tendency due to SSO',&
          &                   'm/s2', 123, 999,&
          &                   vlistID(k_jg),gridCellID(k_jg),zaxisID_hybrid(k_jg)),&
          &           k_jg)
          CALL addVar(TimeVar('tend_v_gwd',&
          &                   'meridional wind tendency due to non-orographic GWD',&
          &                   'm/s2', 124, 999,&
          &                   vlistID(k_jg),gridCellID(k_jg),zaxisID_hybrid(k_jg)),&
          &           k_jg)
          CALL addVar(TimeVar('tend_v_raylfric',&
          &                   'meridional wind tendency due to Rayleigh friction',&
          &                   'm/s2', 125, 999,&
          &                   vlistID(k_jg),gridCellID(k_jg),zaxisID_hybrid(k_jg)),&
          &           k_jg)

        ! Tracer tendencies
        IF (ntracer > 0) THEN
          DO jt = 1, nqtendphy
            ctracer = ctracer_list(jt:jt)
            WRITE(qname,'(A6,A1,A5)') "tend_q",ctracer, "_conv"
            CALL addVar(TimeVar(TRIM(qname), &
                &               TRIM(qname), &
                &               '1/s',22,128,&
                &               vlistID(k_jg), gridCellID(k_jg),zaxisID_hybrid(k_jg)),&
                &       k_jg)
          END DO

          DO jt = 1, nqtendphy
            ctracer = ctracer_list(jt:jt)
            WRITE(qname,'(A6,A1,A5)') "tend_q",ctracer, "_turb"
            CALL addVar(TimeVar(TRIM(qname),&
            &               TRIM(qname),&
            &               '1/s',23,128,&
            &               vlistID(k_jg), gridCellID(k_jg),zaxisID_hybrid(k_jg)),&
            &       k_jg)
          END DO

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

      ! Sea ice - mostly for debuging purposes
      IF ( iice <= nsfc_type ) THEN
        CALL addVar(TimeVar('ice_Tsurf','surface temperature of snow/ice','C',&
        &         100,128,                    &
        &         vlistID(k_jg),gridCellID(k_jg),zaxisID_generic_ice(k_jg)),k_jg)
        CALL addVar(TimeVar('ice_T1','temperature of the upper ice layer','C',&
        &         100,128,                    &
        &         vlistID(k_jg),gridCellID(k_jg),zaxisID_generic_ice(k_jg)),k_jg)
        CALL addVar(TimeVar('ice_T2','temperature of the lower ice layer','C',&
        &         100,128,                    &
        &         vlistID(k_jg),gridCellID(k_jg),zaxisID_generic_ice(k_jg)),k_jg)
       CALL addVar(TimeVar('ice_hi','ice thickness','m',&
        &         100,128,                    &
        &         vlistID(k_jg),gridCellID(k_jg),zaxisID_generic_ice(k_jg)),k_jg)
        CALL addVar(TimeVar('ice_conc','ice concentration in each ice class','',&
        &         100,128,                    &
        &         vlistID(k_jg),gridCellID(k_jg),zaxisID_generic_ice(k_jg)),k_jg)
      ENDIF

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
      CALL addVar(ConstVar('rbasin_c',&
      &                    '2d basin ID on cells',&
      &                    '', 1, 128,&
      &                     vlistID(k_jg),&
      &                     gridCellID(k_jg), &
      &                     zaxisID_surface(k_jg)),&
      &           k_jg)
      CALL addVar(ConstVar('rregio_c',&
      &                    '2d region ID on cells',&
      &                    '', 1, 128,&
      &                     vlistID(k_jg),&
      &                     gridCellID(k_jg), &
      &                     zaxisID_surface(k_jg)),&
      &           k_jg)
   !  CALL addVar(ConstVar('dolic_c',&
   !  &                    'deepest ocean layer on cells',&
   !  &                    '', 1, 128,&
   !  &                     vlistID(k_jg),&
   !  &                     gridCellID(k_jg), &
   !  &                     zaxisID_surface(k_jg)),&
   !  &           k_jg)
   !  CALL addVar(ConstVar('dolic_e',&
   !  &                    'deepest ocean layer on edges',&
   !  &                    '', 1, 128,&
   !  &                     vlistID(k_jg),&
   !  &                     gridEdgeID(k_jg), &
   !  &                     zaxisID_surface(k_jg)),&
   !  &           k_jg)
      CALL addVar(TimeVar('ELEV',&
      &                   'surface elevation at cell center',&
      &                   'm', 1, 128,&
      &                    vlistID(k_jg),&
      &                    gridCellID(k_jg), &
      &                    zaxisID_surface(k_jg)),&
      &           k_jg)
      CALL addVar(TimeVar('Ekin',&
      &                   'kinetic energy at centers',&
      &                   'm/s',4,128,&
      &                   vlistID(k_jg),&
      &                   gridCellID(k_jg),&
      &                   zaxisIDdepth_m(k_jg)),&
      &           k_jg)
  ! IF (iforc_oce > 10) THEN
  !   CALL addVar(TimeVar('forc_u',&
  !   &                   'u-forcing component at centers',&
  !   &                   'N/m2',13,128,&
  !   &                   vlistID(k_jg),&
  !   &                   gridCellID(k_jg),&
  !   &                   zaxisID_surface(k_jg)),&
  !   &           k_jg)
  !   CALL addVar(TimeVar('forc_v',&
  !   &                   'v-forcing component at centers',&
  !   &                   'N/m2',14,128,&
  !   &                   vlistID(k_jg),&
  !   &                   gridCellID(k_jg),&
  !   &                   zaxisID_surface(k_jg)),&
  !   &           k_jg)
  ! END IF
    IF (temperature_relaxation >= 1 ) THEN
      CALL addVar(TimeVar('forc_tdata',&
      &                   'temperature relaxation data',&
      &                   'K',15,128,&
      &                   vlistID(k_jg),&
      &                   gridCellID(k_jg),&
      &                   zaxisID_surface(k_jg)),&
      &           k_jg)
    END IF 
    IF (temperature_relaxation /= 0 ) THEN
      CALL addVar(TimeVar('forc_t',&
      &                   'temperature relaxation flux',&
      &                   'K*m/s',15,128,&
      &                   vlistID(k_jg),&
      &                   gridCellID(k_jg),&
      &                   zaxisID_surface(k_jg)),&
      &           k_jg)
    END IF 
    !IF (temperature_relaxation /= 0 .OR. i_sea_ice >= 1 .OR. i_apply_bulk==1 ) THEN
    IF (iforc_oce > 11) THEN   !  e.g. OMIP forcing or coupled
      CALL addVar(TimeVar('forc_hflx',&
      &                   'net surface heat flux',&
      &                   'W/m2',16,128,&
      &                   vlistID(k_jg),&
      &                   gridCellID(k_jg),&
      &                   zaxisID_surface(k_jg)),&
      &           k_jg)
      CALL addVar(TimeVar('forc_swflx',&
      &                   'short wave surface heat flux',&
      &                   'W/m2',16,128,&
      &                   vlistID(k_jg),&
      &                   gridCellID(k_jg),&
      &                   zaxisID_surface(k_jg)),&
      &           k_jg)
      CALL addVar(TimeVar('forc_lwflx',&
      &                   'long wave surface heat flux',&
      &                   'W/m2',16,128,&
      &                   vlistID(k_jg),&
      &                   gridCellID(k_jg),&
      &                   zaxisID_surface(k_jg)),&
      &           k_jg)
      CALL addVar(TimeVar('forc_ssflx',&
      &                   'sensible surface heat flux',&
      &                   'W/m2',16,128,&
      &                   vlistID(k_jg),&
      &                   gridCellID(k_jg),&
      &                   zaxisID_surface(k_jg)),&
      &           k_jg)
      CALL addVar(TimeVar('forc_slflx',&
      &                   'latent surface heat flux',&
      &                   'W/m2',16,128,&
      &                   vlistID(k_jg),&
      &                   gridCellID(k_jg),&
      &                   zaxisID_surface(k_jg)),&
      &           k_jg)
    END IF 
  ! IF (irelax_2d_S >= 1 ) THEN
  !   CALL addVar(TimeVar('forc_sdata',&
  !   &                   'salinity relaxation data',&
  !   &                   'psu',15,128,&
  !   &                   vlistID(k_jg),&
  !   &                   gridCellID(k_jg),&
  !   &                   zaxisID_surface(k_jg)),&
  !   &           k_jg)
  ! END IF 
  ! IF (irelax_2d_S /= 0 ) THEN
  !   CALL addVar(TimeVar('forc_s',&
  !   &                   'salinity relaxation flux at centers',&
  !   &                   'psu*m/s',15,128,&
  !   &                   vlistID(k_jg),&
  !   &                   gridCellID(k_jg),&
  !   &                   zaxisID_surface(k_jg)),&
  !   &           k_jg)
  ! END IF 
    IF (irelax_2d_S /= 0 ) THEN
      CALL addVar(TimeVar('forc_fwfx',&
      &                   'diagnosed net freshwater flux',&
!     &                   'm/s',16,128,&
      &                   'm/month',16,128,&
      &                   vlistID(k_jg),&
      &                   gridCellID(k_jg),&
      &                   zaxisID_surface(k_jg)),&
      &           k_jg)
    END IF 
   !  CALL addVar(TimeVar('horz_adv',&
   !  &                   'nonlin Cor ',&
   !  &                   'm/s',2,128,&
   !  &                   vlistID(k_jg),&
   !  &                   gridEdgeID(k_jg), &
   !  &                   zaxisIDdepth_m(k_jg)),&
   !  &           k_jg)

   !  CALL addVar(TimeVar('Ekin_grad',&
   !  &                   'gradient Ekin ',&
   !  &                   'm/s',2,128,&
   !  &                   vlistID(k_jg),&
   !  &                   gridEdgeID(k_jg), &
   !  &                   zaxisIDdepth_m(k_jg)),&
   !  &           k_jg)
   !  CALL addVar(TimeVar('flux_u',&
   !  &                   'sfc_flux u at centers',&
   !  &                   'N/m2',13,128,&
   !  &                   vlistID(k_jg),&
   !  &                   gridCellID(k_jg),&
   !  &                   zaxisID_surface(k_jg)),&
   !  &           k_jg)
   !  CALL addVar(TimeVar('flux_v',&
   !  &                   'sfc_flux v at centers',&
   !  &                   'N/m2',14,128,&
   !  &                   vlistID(k_jg),&
   !  &                   gridCellID(k_jg),&
   !  &                   zaxisID_surface(k_jg)),&
   !  &           k_jg)

   !  CALL addVar(TimeVar('flux_VN',&
   !  &                   'sfc-flux at edge',&
   !  &                   'm/s',2,128,&
   !  &                   vlistID(k_jg),&
   !  &                   gridEdgeID(k_jg), &
   !  &                   zaxisID_surface(k_jg)),&
   !  &           k_jg)
      CALL addVar(TimeVar('VN',&
      &                   'velocity normal to edge',&
      &                   'm/s',2,128,&
      &                   vlistID(k_jg),&
      &                   gridEdgeID(k_jg), &
      &                   zaxisIDdepth_m(k_jg)),&
      &           k_jg)
      CALL addVar(TimeVar('VORT',&
      &                   'vorticity at vertices',&
      &                   '1/s',3,128,&
      &                   vlistID(k_jg),&
      &                   gridVertexID(k_jg),&
      &                   zaxisIDdepth_m(k_jg)),&
      &          k_jg)
      CALL addVar(TimeVar('u',&
      &                   'u-velocity component at centers',&
      &                   'm/s',4,128,&
      &                   vlistID(k_jg),&
      &                   gridCellID(k_jg),&
      &                   zaxisIDdepth_m(k_jg)),&
      &           k_jg)
      CALL addVar(TimeVar('v',&
      &                   'v-velocity component at centers',&
      &                   'm/s',5,128,&
      &                   vlistID(k_jg), &
      &                   gridCellID(k_jg), &
      &                   zaxisIDdepth_m(k_jg)),&
      &           k_jg)
    ! CALL addVar(TimeVar('Vert_Veloc_Adv',&
    ! &                   'vertical-velocity advection at edges',&
    ! &                   'm/s',5,128,&
    ! &                   vlistID(k_jg), &
    ! &                   gridEdgeID(k_jg), &
    ! &                   zaxisIDdepth_m(k_jg)),&
    ! &           k_jg)
      CALL addVar(TimeVar('u_vint',&
      &                   'barotropic zonal velocity at centers',&
      &                   'm*m/s',4,128,&
      &                   vlistID(k_jg),&
      &                   gridCellID(k_jg),&
      &                   zaxisID_surface(k_jg)),&
      &           k_jg)
      CALL addVar(TimeVar('W',&
      &                   'vertical velocity at cells',&
      &                   'm/s', 6, 128,&
      &                   vlistid(k_jg),&
      &                   gridCellID(k_jg),&
      &                   zaxisID_halfdepth(k_jg)),&
      &           k_jg)      
      CALL addVar(TimeVar('Vert_Mixing_V',&
      &                   'vertical mixing coeff veloc',&
      &                   'm^2/s', 6, 128,&
      &                   vlistid(k_jg),&
      &                   gridEdgeID(k_jg),&
      &                   zaxisID_halfdepth(k_jg)),&
      &           k_jg)
      CALL addVar(TimeVar('Horz_Mixing_V',&
      &                   'horizontal mixing coeff veloc',&
      &                   'm^2/s', 6, 128,&
      &                   vlistid(k_jg),&
      &                   gridEdgeID(k_jg),&
      &                   zaxisIDdepth_m(k_jg)),&
      &           k_jg)
      IF (no_tracer > 0) THEN
        CALL addVar(TimeVar('Vert_Mixing_T',&
        &                   'vertical mixing coeff temp',&
        &                   'm^2/s', 6, 128,&
        &                   vlistid(k_jg),&
        &                   gridCellID(k_jg),&
        &                   zaxisID_halfdepth(k_jg)),&
        &           k_jg)
      END IF
      CALL addVar(TimeVar('press_grad',&
      &                   'pressure-gradient at edges',&
      &                   'm/s',5,128,&
      &                   vlistID(k_jg), &
      &                   gridEdgeID(k_jg), &
      &                   zaxisIDdepth_m(k_jg)),&
      &                   k_jg)
      CALL addVar(TimeVar('rho',&
      &                   'density cells',&
      &                   'kg/m**3', 6, 128,&
      &                   vlistid(k_jg),&
      &                   gridCellID(k_jg),&
      &                   zaxisIDdepth_m(k_jg)),&
      &           k_jg)



   ! sea ice
     IF (i_sea_ice >= 1 ) THEN
!       CALL addVar(TimeVar('p_ice_isice','','',&
!       &         100,128,                    &
!       &         vlistID(k_jg),gridCellID(k_jg),zaxisID_generic_ice(k_jg)),k_jg)
!      CALL addVar(TimeVar('p_ice_alb','albedo of the snow/ice system','',&
!      &         100,128,                    &
!      &         vlistID(k_jg),gridCellID(k_jg),zaxisID_generic_ice(k_jg)),k_jg)
       CALL addVar(TimeVar('p_ice_Tsurf','surface temperature of snow/ice','C',&
       &         100,128,                    &
       &         vlistID(k_jg),gridCellID(k_jg),zaxisID_generic_ice(k_jg)),k_jg)
!      CALL addVar(TimeVar('p_ice_T1','temperature of the upper ice layer','C',&
!      &         100,128,                    &
!      &         vlistID(k_jg),gridCellID(k_jg),zaxisID_generic_ice(k_jg)),k_jg)
!      CALL addVar(TimeVar('p_ice_T2','temperature of the lower ice layer','C',&
!      &         100,128,                    &
!      &         vlistID(k_jg),gridCellID(k_jg),zaxisID_generic_ice(k_jg)),k_jg)
!       CALL addVar(TimeVar('p_ice_E1','energy content of the upper ice layer','Jm/kg',&
!       &         100,128,                    &
!       &         vlistID(k_jg),gridCellID(k_jg),zaxisID_generic_ice(k_jg)),k_jg)
!       CALL addVar(TimeVar('p_ice_E2','energy content of the second ice layer','Jm/kg',&
!       &         100,128,                    &
!       &         vlistID(k_jg),gridCellID(k_jg),zaxisID_generic_ice(k_jg)),k_jg)
       CALL addVar(TimeVar('p_ice_hi','ice thickness','m',&
       &         100,128,                    &
       &         vlistID(k_jg),gridCellID(k_jg),zaxisID_generic_ice(k_jg)),k_jg)
!      CALL addVar(TimeVar('p_ice_hs','snow thickness','m',&
!      &         100,128,                    &
!      &         vlistID(k_jg),gridCellID(k_jg),zaxisID_generic_ice(k_jg)),k_jg)
!       CALL addVar(TimeVar('p_ice_Qtop','energy flux available for surface melting','W/m2',&
!       &         100,128,                    &
!       &         vlistID(k_jg),gridCellID(k_jg),zaxisID_generic_ice(k_jg)),k_jg)
!       CALL addVar(TimeVar('p_ice_Qbot','energy flux at ice-ocean interface','W/m2',&
!       &         100,128,                    &
!       &         vlistID(k_jg),gridCellID(k_jg),zaxisID_generic_ice(k_jg)),k_jg)
!       CALL addVar(TimeVar('p_ice_heatocei','energy to ocean when all ice is melted','J',&
!       &         100,128,                    &
!       &         vlistID(k_jg),gridCellID(k_jg),zaxisID_generic_ice(k_jg)),k_jg)
!       CALL addVar(TimeVar('p_ice_snow_to_ice','amount of snow that is transformed to ice','m',&
!       &         100,128,                    &
!       &         vlistID(k_jg),gridCellID(k_jg),zaxisID_generic_ice(k_jg)),k_jg)
!       CALL addVar(TimeVar('p_ice_surfmelt','surface melt water running into ocean','m',&
!       &         100,128,                    &
!       &         vlistID(k_jg),gridCellID(k_jg),zaxisID_generic_ice(k_jg)),k_jg)
!       CALL addVar(TimeVar('p_ice_surfmeltT','mean temperature of surface melt water','C',&
!       &         100,128,                    &
!       &         vlistID(k_jg),gridCellID(k_jg),zaxisID_generic_ice(k_jg)),k_jg)
!       CALL addVar(TimeVar('p_ice_evapwi','amount of evaporated water if no ice is left','kg/m2',&
!       &         100,128,                    &
!       &         vlistID(k_jg),gridCellID(k_jg),zaxisID_generic_ice(k_jg)),k_jg)
!       CALL addVar(TimeVar('p_ice_conc','ice concentration in each ice class','',&
!       &         100,128,                    &
!       &         vlistID(k_jg),gridCellID(k_jg),zaxisID_generic_ice(k_jg)),k_jg)
     ! 2D
!       CALL addVar(TimeVar('p_ice_u','zonal velocity','m/s',&
!       &         100,128,                    &
!       &         vlistID(k_jg),gridCellID(k_jg),zaxisID_surface(k_jg)),k_jg)
!       CALL addVar(TimeVar('p_ice_v','meridional velocity','m/s',&
!       &         100,128,                    &
!       &         vlistID(k_jg),gridCellID(k_jg),zaxisID_surface(k_jg)),k_jg)
       CALL addVar(TimeVar('p_ice_concSum','total ice concentration','',&
       &         100,128,                    &
       &         vlistID(k_jg),gridCellID(k_jg),zaxisID_surface(k_jg)),k_jg)
!       CALL addVar(TimeVar('p_ice_newice','new-ice growth in open water','m',&
!       &         100,128,                    &
!       &         vlistID(k_jg),gridCellID(k_jg),zaxisID_surface(k_jg)),k_jg)
!       CALL addVar(TimeVar('p_ice_zUnderIce','water in upper ocean grid cell below ice','m',&
!       &         100,128,                    &
!       &         vlistID(k_jg),gridCellID(k_jg),zaxisID_surface(k_jg)),k_jg)
  !    CALL addVar(TimeVar('p_ice_hi_lim','','',&
  !    &         100,128,                    &
  !    &         vlistID(k_jg),gridCellID(k_jg),zaxisID_surface(k_jg)),k_jg)
      ENDIF

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

        IF ( lwrite_oce_timestepping ) THEN
          CALL addVar(TimeVar('g_n_c_v_'//TRIM(oce_tracer_names(itracer)),&
            &'g_n_c_v for '//TRIM(oce_tracer_names(itracer)),'',255,128,&
            &vlistID(k_jg),gridCellID(k_jg),zaxisIDdepth_m(k_jg)),k_jg)
          CALL addVar(TimeVar('g_n_c_h_'//TRIM(oce_tracer_names(itracer)),&
            &'g_n_c_h for '//TRIM(oce_tracer_names(itracer)),'',255,128,&
            &vlistID(k_jg),gridCellID(k_jg),zaxisIDdepth_m(k_jg)),k_jg)
          CALL addVar(TimeVar('g_nm1_c_v_'//TRIM(oce_tracer_names(itracer)),&
            &'g_nm1_c_v for '//TRIM(oce_tracer_names(itracer)),'',255,128,&
            &vlistID(k_jg),gridCellID(k_jg),zaxisIDdepth_m(k_jg)),k_jg)
          CALL addVar(TimeVar('g_nm1_c_h_'//TRIM(oce_tracer_names(itracer)),&
            &'g_nm1_c_h for '//TRIM(oce_tracer_names(itracer)),'',255,128,&
            &vlistID(k_jg),gridCellID(k_jg),zaxisIDdepth_m(k_jg)),k_jg)
          CALL addVar(TimeVar('g_nimd_c_v_'//TRIM(oce_tracer_names(itracer)),&
            &'g_nimd_c_v for '//TRIM(oce_tracer_names(itracer)),'',255,128,&
            &vlistID(k_jg),gridCellID(k_jg),zaxisIDdepth_m(k_jg)),k_jg)
          CALL addVar(TimeVar('g_nimd_c_h_'//TRIM(oce_tracer_names(itracer)),&
            &'g_nimd_c_h for '//TRIM(oce_tracer_names(itracer)),'',255,128,&
            &vlistID(k_jg),gridCellID(k_jg),zaxisIDdepth_m(k_jg)),k_jg)
        END IF
      END DO
      IF ( lwrite_oce_timestepping ) THEN
        CALL addVar(TimeVar('g_n','g_n','',255,128,&
          &vlistID(k_jg),gridEdgeID(k_jg),zaxisIDdepth_m(k_jg)),k_jg)
        CALL addVar(TimeVar('g_nm1','g_m1n','',255,128,&
          &vlistID(k_jg),gridEdgeID(k_jg),zaxisIDdepth_m(k_jg)),k_jg)
        CALL addVar(TimeVar('g_nimd','g_nimd','',255,128,&
          &vlistID(k_jg),gridEdgeID(k_jg),zaxisIDdepth_m(k_jg)),k_jg)
      END IF

    END IF  ! ocean

    IF (lwrite_decomposition) THEN
      CALL addVar(ConstVar('cell_owner',&
      &                    'MPI-ownership of cells',&
      &                    '', 1, 128,&
      &                     vlistID(k_jg),&
      &                     gridCellID(k_jg), &
      &                     zaxisID_surface(k_jg)),k_jg)
!     CALL addVar(ConstVar('edge_owner',&
!     &                    'MPI-ownership of edges',&
!     &                    '', 1, 128,&
!     &                     vlistID(k_jg),&
!     &                     gridEdgeID(k_jg), &
!     &                     zaxisID_surface(k_jg)),k_jg)
!     CALL addVar(ConstVar('vertex_owner',&
!     &                    'MPI-ownership of vertices',&
!     &                    '', 1, 128,&
!     &                     vlistID(k_jg),&
!     &                     gridVertexID(k_jg), &
!     &                     zaxisID_surface(k_jg)),k_jg)
      ALLOCATE(cell_owner(nproma, p_patch(k_jg)%nblks_c))
!     ALLOCATE(edge_owner(nproma, p_patch(1)%nblks_e))
!     ALLOCATE(vertex_owner(nproma, p_patch(1)%nblks_v))
      cell_owner = get_my_mpi_all_id()
    ENDIF

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
        CALL finish(routine,'got illegal gridid')
      ENDIF

      zaxisid = vlistInqVarZaxis(vlistID(k_jg), varids(ivar, k_jg))
      outvar_desc(ivar, k_jg)%nlev = zaxisInqSize(zaxisid)

      CALL vlistInqVarName(vlistID(k_jg), varids(ivar, k_jg), outvar_desc(ivar, k_jg)%name)

      ! write double precision
      IF ( lwrite_dblprec ) CALL vlistDefVarDatatype(vlistID(k_jg),&
        &                                            varids(ivar,k_jg),&
        &                                            DATATYPE_FLT64)

    ENDDO

  END SUBROUTINE setup_vlist


  !-------------------------------------------------------------------------------------------------
  SUBROUTINE destruct_vlist(k_jg)

    INTEGER, INTENT(in) :: k_jg

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

  SUBROUTINE vlist_write_var(ivar, k_jg, var)
    INTEGER,       INTENT(in)  :: ivar, k_jg
    REAL(wp),      INTENT(in)  :: var(*)

    ! local variables
    CHARACTER(*), PARAMETER :: routine = TRIM("mo_io_vlist:vlist_write_var")
    INTEGER            :: n_tot

    CALL streamWriteVar(streamID(k_jg), varids(ivar, k_jg), var, 0)

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

!      WRITE(0,*)'varname=',varname
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
!!$ TR: JSBACH testing
      CASE ('surface_temperature');  ptr2 => prm_field(jg)%surface_temperature
      CASE ('surface_temperature_old');  ptr2 => prm_field(jg)%surface_temperature_old
      CASE ('surface_temperature_rad');  ptr2 => prm_field(jg)%surface_temperature_rad
      CASE ('surface_temperature_eff');  ptr2 => prm_field(jg)%surface_temperature_eff

      CASE ('c_soil_temperature1');  ptr2 => prm_field(jg)%c_soil_temperature1
      CASE ('c_soil_temperature2');  ptr2 => prm_field(jg)%c_soil_temperature2
      CASE ('c_soil_temperature3');  ptr2 => prm_field(jg)%c_soil_temperature3
      CASE ('c_soil_temperature4');  ptr2 => prm_field(jg)%c_soil_temperature4
      CASE ('c_soil_temperature5');  ptr2 => prm_field(jg)%c_soil_temperature5

      CASE ('d_soil_temperature1');  ptr2 => prm_field(jg)%d_soil_temperature1
      CASE ('d_soil_temperature2');  ptr2 => prm_field(jg)%d_soil_temperature2
      CASE ('d_soil_temperature3');  ptr2 => prm_field(jg)%d_soil_temperature3
      CASE ('d_soil_temperature4');  ptr2 => prm_field(jg)%d_soil_temperature4
      CASE ('d_soil_temperature5');  ptr2 => prm_field(jg)%d_soil_temperature5

      CASE ('soil_temperature1');  ptr2 => prm_field(jg)%soil_temperature1
      CASE ('soil_temperature2');  ptr2 => prm_field(jg)%soil_temperature2
      CASE ('soil_temperature3');  ptr2 => prm_field(jg)%soil_temperature3
      CASE ('soil_temperature4');  ptr2 => prm_field(jg)%soil_temperature4
      CASE ('soil_temperature5');  ptr2 => prm_field(jg)%soil_temperature5

      CASE ('heat_capacity');  ptr2 => prm_field(jg)%heat_capacity
      CASE ('ground_heat_flux');  ptr2 => prm_field(jg)%ground_heat_flux
      CASE ('swnet');  ptr2 => prm_field(jg)%swnet
      CASE ('time_steps_soil');  ptr2 => prm_field(jg)%time_steps_soil
      CASE ('moisture1');   ptr2 => prm_field(jg)%moisture1
      CASE ('moisture2');   ptr2 => prm_field(jg)%moisture2
      CASE ('moisture3');   ptr2 => prm_field(jg)%moisture3
      CASE ('moisture4');   ptr2 => prm_field(jg)%moisture4
      CASE ('moisture5');   ptr2 => prm_field(jg)%moisture5
      CASE ('moisture_all');   ptr2 => prm_field(jg)%moisture_all
      CASE ('skin_reservoir');   ptr2 => prm_field(jg)%skin_reservoir
      CASE ('snow');   ptr2 => prm_field(jg)%snow
      CASE ('albvisdir');   ptr2 => prm_field(jg)%albvisdir
      CASE ('albnirdir');   ptr2 => prm_field(jg)%albnirdir
        !KF  the reset command can only be used for 'plain' fields
      CASE ('swflxsfc_avg')
                                ptr2 => dup2(prm_field(jg)% swflxsfc_avg(:,:)/dt_data)
                                             prm_field(jg)% swflxsfc_avg(:,:)=0.0_wp
      CASE ('lwflxsfc_avg') 
                               ptr2 => dup2(prm_field(jg)% lwflxsfc_avg(:,:)/dt_data)
                                             prm_field(jg)% lwflxsfc_avg(:,:)=0.0_wp
      CASE ('dlwfsfc_dT_avg') 
                               ptr2 => dup2(prm_field(jg)%dlwflxsfc_dT_avg(:,:)/dt_data)
                                             prm_field(jg)%dlwflxsfc_dT_avg(:,:)=0.0_wp
      CASE ('swflxtoa_avg')
                               ptr2 => dup2(prm_field(jg)% swflxtoa_avg(:,:)/dt_data)
                                            prm_field(jg)% swflxtoa_avg(:,:) = 0.0_wp
      CASE ('lwflxtoa_avg')
                               ptr2 => dup2(prm_field(jg)% lwflxtoa_avg(:,:)/dt_data)
                                            prm_field(jg)% lwflxtoa_avg(:,:) = 0.0_wp
      CASE ('TOTPREC_AVG')
                               ptr2 => dup2(prm_field(jg)% totprec_avg(:,:)/dt_data)
                                       prm_field(jg)% totprec_avg(:,:) = 0.0_wp
      CASE ('evap_avg')
                               ptr2 => dup2(prm_field(jg)%    evap_avg(:,:)/dt_data)
                               prm_field(jg)%    evap_avg(:,:) = 0.0_wp
      CASE ('lhflx_avg')
                               ptr2 => dup2(prm_field(jg)%   lhflx_avg(:,:)/dt_data)
                                            prm_field(jg)%   lhflx_avg(:,:) = 0.0_wp
      CASE ('shflx_avg')
                              ptr2 => dup2(prm_field(jg)%   shflx_avg(:,:)/dt_data)
                                           prm_field(jg)%   shflx_avg(:,:) = 0.0_wp
      CASE ('u_stress_avg') 
                               ptr2 => dup2(prm_field(jg)%u_stress_avg(:,:)/dt_data)
                                            prm_field(jg)%u_stress_avg(:,:) = 0.0_wp
      CASE ('v_stress_avg')  
                               ptr2 => dup2(prm_field(jg)%v_stress_avg(:,:)/dt_data)
                               prm_field(jg)%v_stress_avg(:,:) = 0.0_wp

      CASE ('swflxsfc');    ptr2 => prm_field(jg)% swflxsfc(:,:)
      CASE ('lwflxsfc');    ptr2 => prm_field(jg)% lwflxsfc(:,:)
      CASE ('swflxtoa');    ptr2 => prm_field(jg)% swflxtoa(:,:)
      CASE ('lwflxtoa');    ptr2 => prm_field(jg)% lwflxtoa(:,:)

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
      CASE ('cell_owner');      ptr2 => cell_owner
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

      ! Sea ice - mostly for debuging purposes
      CASE ('ice_Tsurf');       ptr3 => prm_field(jg)%Tsurf
      CASE ('ice_T1');          ptr3 => prm_field(jg)%T1
      CASE ('ice_T2');          ptr3 => prm_field(jg)%T2
      CASE ('ice_hi');          ptr3 => prm_field(jg)%hi
      CASE ('ice_conc');        ptr3 => prm_field(jg)%conc

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

      DO jt = 1, ntracer
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

        IF(varname == 'dshflx_dT_avg_tile'//ctile) THEN
          ptr2 => dup2(prm_field(jg)%dshflx_dT_avg_tile(:,:,jt)/dt_data)
          prm_field(jg)%dshflx_dT_avg_tile(:,:,jt)= 0.0_wp
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

  SUBROUTINE get_outvar_ptr_nh(varname, jg, ptr2, ptr3, reset, delete)

    CHARACTER(LEN=*), INTENT(IN) :: varname
    INTEGER, INTENT(IN) :: jg
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

    TYPE(t_nh_prog),    POINTER :: p_prog, p_prog_rcf
    TYPE(t_nh_diag),    POINTER :: p_diag

    TYPE(t_lnd_prog),   POINTER :: p_prog_lnd
    TYPE(t_lnd_diag),   POINTER :: p_diag_lnd

    p_prog     => p_nh_state(jg)%prog(nnow(jg))
    p_prog_rcf => p_nh_state(jg)%prog(nnow_rcf(jg))
    p_diag     => p_nh_state(jg)%diag

    IF (iforcing==inwp) THEN
      p_prog_lnd     => p_lnd_state(jg)%prog_lnd(nnow_rcf(jg))
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
      CASE ('T_2M' );           ptr2 => prm_diag(jg)%t_2m
      CASE ('T_2M_s6avg' );     ptr2 => prm_diag(jg)%t_2m_s6avg
      CASE ('QV_2M' );           ptr2 => prm_diag(jg)%qv_2m
      CASE ('QV_2M_s6avg' );     ptr2 => prm_diag(jg)%qv_2m_s6avg
      CASE ('normal_velocity'); ptr3 => p_prog%vn
      CASE ('U');               ptr3 => p_diag%u
      CASE ('V');               ptr3 => p_diag%v
      CASE ('W');               ptr3 => p_prog%w
      CASE ('U_10M');           ptr2 => prm_diag(jg)%u_10m
      CASE ('V_10M');           ptr2 => prm_diag(jg)%v_10m
      CASE ('U_10M_s6avg');     ptr2 => prm_diag(jg)%u_10m_s6avg
      CASE ('V_10M_s6avg');     ptr2 => prm_diag(jg)%v_10m_s6avg
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
      CASE ('PRR_GSP');         ptr2 => prm_diag(jg)%rain_gsp_rate(:,:)
      CASE ('PRS_GSP');         ptr2 => prm_diag(jg)%snow_gsp_rate(:,:)
      CASE ('RAIN_GSP');        ptr2 => prm_diag(jg)%rain_gsp(:,:)
      CASE ('SNOW_GSP');        ptr2 => prm_diag(jg)%snow_gsp(:,:)
      CASE ('RAIN_CON');        ptr2 => prm_diag(jg)%rain_con(:,:)
      CASE ('SNOW_CON');        ptr2 => prm_diag(jg)%snow_con(:,:)
      CASE ('TOT_PREC');        ptr2 => prm_diag(jg)%tot_prec(:,:)
      CASE ('TOT_PREC_RATE_avg'); ptr2 => prm_diag(jg)%tot_prec_rate_avg(:,:)
      CASE ('CON_PREC_RATE_avg'); ptr2 => prm_diag(jg)%con_prec_rate_avg(:,:)
      CASE ('GSP_PREC_RATE_avg'); ptr2 => prm_diag(jg)%gsp_prec_rate_avg(:,:)
      CASE ('ALB_RAD');         ptr2 => prm_diag(jg)%albvisdif(:,:)
      CASE ('cosmu0');          ptr2 => prm_diag(jg)%cosmu0(:,:)
      CASE ('flxdwswtoa');      ptr2 => prm_diag(jg)%flxdwswtoa(:,:)
      CASE ('SOB_S');           ptr2 => prm_diag(jg)%swflxsfc(:,:)
      CASE ('THB_S');           ptr2 => prm_diag(jg)%lwflxsfc(:,:)
      CASE ('SOB_T');           ptr2 => prm_diag(jg)%swflxtoa(:,:)
      CASE ('THB_T');           ptr2 => prm_diag(jg)%lwflxall(:,1,:)
      CASE ('ASOB_S');          ptr2 => prm_diag(jg)%swflxsfc_a(:,:)
      CASE ('ATHB_S');          ptr2 => prm_diag(jg)%lwflxsfc_a(:,:)       
      CASE ('ASOB_T');          ptr2 => prm_diag(jg)%swflxtoa_a(:,:)
      CASE ('ATHB_T');          ptr2 => prm_diag(jg)%lwflxtoa_a(:,:)
      CASE ('ACCSOB_S');        ptr2 => prm_diag(jg)%swflxsfc_a(:,:)
      CASE ('ACCTHB_S');        ptr2 => prm_diag(jg)%lwflxsfc_a(:,:)       
      CASE ('ACCSOB_T');        ptr2 => prm_diag(jg)%swflxtoa_a(:,:)
      CASE ('ACCTHB_T');        ptr2 => prm_diag(jg)%lwflxtoa_a(:,:)
      CASE ('T_G');             ptr2 => p_prog_lnd%t_g
      CASE ('QV_S');            ptr2 => p_diag_lnd%qv_s
      CASE ('ASHFL_S');         ptr2 => prm_diag(jg)%shfl_s_a 
      CASE ('ALHFL_S');         ptr2 => prm_diag(jg)%lhfl_s_a
      CASE ('ACCSHFL_S');       ptr2 => prm_diag(jg)%shfl_s_a 
      CASE ('ACCLHFL_S');       ptr2 => prm_diag(jg)%lhfl_s_a
      CASE ('SHFL_S')
       IF   (atm_phy_nwp_config(jg)%inwp_turb.EQ.1) THEN  
         ptr2 => dup2(-1.*prm_diag(jg)%shfl_s(:,:)); delete = .TRUE.
       ELSEIF  (atm_phy_nwp_config(jg)%inwp_turb.EQ.4) THEN 
         ptr2 => prm_diag(jg)%shfl_s
       ELSE
         ptr2 => prm_diag(jg)%shfl_s
       ENDIF
      CASE ('LHFL_S')
       IF   (atm_phy_nwp_config(jg)%inwp_turb.EQ.1) THEN  
         ptr2 =>  dup2(-1.*prm_diag(jg)%lhfl_s(:,:)); delete = .TRUE.
       ELSEIF  (atm_phy_nwp_config(jg)%inwp_turb.EQ.4) THEN 
         ptr2 => prm_diag(jg)%lhfl_s
       ELSE
         ptr2 => prm_diag(jg)%lhfl_s
       ENDIF
      CASE ('EVAP_RATE_avg');   ptr2 => prm_diag(jg)%qhfl_s_avg     
      CASE ('VOR');             ptr3 => p_diag%omega_z
      CASE ('DIV');             ptr3 => p_diag%div
      CASE ('THETA_V');         ptr3 => p_prog%theta_v
      CASE ('EXNER');           ptr3 => p_prog%exner
      CASE ('RHO');             ptr3 => p_prog%rho
      CASE ('TKE');             ptr3 => p_prog_rcf%tke
      CASE ('TCM');             ptr2 => prm_diag(jg)%tcm
      CASE ('TCH');             ptr2 => prm_diag(jg)%tch
      CASE ('Z0')
        IF (atm_phy_nwp_config(jg)%inwp_turb.EQ.1 .OR.  &
         &  atm_phy_nwp_config(jg)%inwp_turb.EQ.2) THEN
                                ptr2 => dup2(prm_diag(jg)%gz0(:,:)/grav); delete = .TRUE.
                              ELSE
                                ptr2 => prm_diag(jg)%z0m(:,:)
        END IF
      CASE ('tend_temp_radsw'); ptr3 => prm_nwp_tend(jg)%ddt_temp_radsw(:,:,:)
      CASE ('tend_temp_radlw'); ptr3 => prm_nwp_tend(jg)%ddt_temp_radlw(:,:,:)
      CASE ('tend_temp_conv');  ptr3 => prm_nwp_tend(jg)%ddt_temp_pconv(:,:,:)
      CASE ('tend_temp_turb');  ptr3 => prm_nwp_tend(jg)%ddt_temp_turb (:,:,:)
      CASE ('tend_temp_drag');  ptr3 => prm_nwp_tend(jg)%ddt_temp_drag (:,:,:)
      CASE ('tend_temp_dyn');   ptr3 => p_diag%ddt_temp_dyn(:,:,:)
      CASE ('tend_u_conv');     ptr3 => prm_nwp_tend(jg)%ddt_u_pconv   (:,:,:)
      CASE ('tend_u_turb');     ptr3 => prm_nwp_tend(jg)%ddt_u_turb    (:,:,:)
      CASE ('tend_u_sso');      ptr3 => prm_nwp_tend(jg)%ddt_u_sso     (:,:,:)
      CASE ('tend_u_gwd');      ptr3 => prm_nwp_tend(jg)%ddt_u_gwd     (:,:,:)
      CASE ('tend_u_raylfric'); ptr3 => prm_nwp_tend(jg)%ddt_u_raylfric(:,:,:)
      CASE ('tend_v_conv');     ptr3 => prm_nwp_tend(jg)%ddt_v_pconv   (:,:,:)
      CASE ('tend_v_turb');     ptr3 => prm_nwp_tend(jg)%ddt_v_turb    (:,:,:)
      CASE ('tend_v_sso');      ptr3 => prm_nwp_tend(jg)%ddt_v_sso     (:,:,:)
      CASE ('tend_v_gwd');      ptr3 => prm_nwp_tend(jg)%ddt_v_gwd     (:,:,:)
      CASE ('tend_v_raylfric'); ptr3 => prm_nwp_tend(jg)%ddt_v_raylfric(:,:,:)
      CASE ('cell_owner');      ptr2 => cell_owner
      CASE DEFAULT;             not_found = .TRUE.
    END SELECT

    ! If not found in the list above, check for tracers


    ! get ctracer_list
    ctracer_list = advection_config(jg)%ctracer_list

    IF(not_found) THEN
      DO jt = 1, ntracer ! all tracer
        ctracer = ctracer_list(jt:jt)
        IF(varname == 'Q'//ctracer) THEN
          ptr3 => p_prog_rcf%tracer(:,:,:,jt)
          RETURN
        ENDIF
      ENDDO

      DO jt = 1, 5
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
      DO jt = 1, nqtendphy ! only hydromet-tendencies (qv,qc,qi)
        ctracer = ctracer_list(jt:jt)
        IF(varname == 'tend_q'//ctracer//'_conv') THEN
          ptr3 => prm_nwp_tend(jg)%ddt_tracer_pconv(:,:,:,jt)
          RETURN
        ENDIF
        IF(varname == 'tend_q'//ctracer//'_turb') THEN
          ptr3 => prm_nwp_tend(jg)%ddt_tracer_turb(:,:,:,jt)
          RETURN
        ENDIF
      ENDDO

      IF ( irad_o3 == 4 .OR. irad_o3 == 6 .OR. irad_o3 == 7 ) THEN
        IF(varname == 'O3') THEN
          ptr3 => ext_data(jg)%atm%o3(:,:,:)
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

      DO jt = 1, ntiles_total
        WRITE(cjt, '(i2)') jt
        WRITE(name,'(A,A)') "T_GT_tile_", TRIM(ADJUSTL(cjt))
        IF(varname == TRIM(name)) THEN
          ptr2 => p_prog_lnd%t_g_t(:,:,jt)
          RETURN
        ENDIF

        WRITE(name,'(A,A)') "T_SNOW_tile_", TRIM(ADJUSTL(cjt))
        IF(varname == TRIM(name)) THEN
          ptr2 => p_prog_lnd%t_snow_t(:,:,jt)
          RETURN
        ENDIF

        WRITE(name,'(A,A)') "T_SNOW_MULT_tile_", TRIM(ADJUSTL(cjt))
        IF(varname == TRIM(name)) THEN
          ptr3 => p_prog_lnd%t_snow_mult_t(:,:,:,jt)
          RETURN
        ENDIF

        WRITE(name,'(A,A)') "T_S_tile_", TRIM(ADJUSTL(cjt))
        IF(varname == TRIM(name)) THEN
          ptr2 => p_prog_lnd%t_s_t(:,:,jt)
          RETURN
        ENDIF

        WRITE(name,'(A,A)') "W_SNOW_tile_", TRIM(ADJUSTL(cjt))
        IF(varname == TRIM(name)) THEN
          ptr2 => p_prog_lnd%w_snow_t(:,:,jt)
          RETURN
        ENDIF

        WRITE(name,'(A,A)') "RHO_SNOW_tile_", TRIM(ADJUSTL(cjt))
        IF(varname == TRIM(name)) THEN
          ptr2 => p_prog_lnd%rho_snow_t(:,:,jt)
          RETURN
        ENDIF

        WRITE(name,'(A,A)') "RHO_SNOW_MULT_tile_", TRIM(ADJUSTL(cjt))
        IF(varname == TRIM(name)) THEN
          ptr3 => p_prog_lnd%rho_snow_mult_t(:,:,:,jt)
          RETURN
        ENDIF

        WRITE(name,'(A,A)') "W_I_tile_", TRIM(ADJUSTL(cjt))
        IF(varname == TRIM(name)) THEN
          ptr2 => p_prog_lnd%w_i_t(:,:,jt)
          RETURN
        ENDIF

        WRITE(name,'(A,A)') "T_SO_tile_", TRIM(ADJUSTL(cjt))
        IF(varname == TRIM(name)) THEN
          ptr3 => p_prog_lnd%t_so_t(:,:,:,jt)
          RETURN
        ENDIF

        WRITE(name,'(A,A)') "W_SO_tile_", TRIM(ADJUSTL(cjt))
        IF(varname == TRIM(name)) THEN
          ptr3 => p_prog_lnd%w_so_t(:,:,:,jt)
          RETURN
        ENDIF

        WRITE(name,'(A,A)') "W_SO_ICE_tile_", TRIM(ADJUSTL(cjt))
        IF(varname == TRIM(name)) THEN
          ptr3 => p_prog_lnd%w_so_ice_t(:,:,:,jt)
          RETURN
        ENDIF

        WRITE(name,'(A,A)') "WLIQ_SNOW_tile_", TRIM(ADJUSTL(cjt))
        IF(varname == TRIM(name)) THEN
          ptr3 => p_prog_lnd%wliq_snow_t(:,:,:,jt)
          RETURN
        ENDIF

        WRITE(name,'(A,A)') "WTOT_SNOW_tile_", TRIM(ADJUSTL(cjt))
        IF(varname == TRIM(name)) THEN
          ptr3 => p_prog_lnd%wtot_snow_t(:,:,:,jt)
          RETURN
        ENDIF

        WRITE(name,'(A,A)') "H_SNOW_tile_", TRIM(ADJUSTL(cjt))
        IF(varname == TRIM(name)) THEN
          ptr2 => p_diag_lnd%h_snow_t(:,:,jt)
          RETURN
        ENDIF

        WRITE(name,'(A,A)') "FRESHSNOW_tile_", TRIM(ADJUSTL(cjt))
        IF(varname == TRIM(name)) THEN
          ptr2 => p_diag_lnd%freshsnow_t(:,:,jt)
          RETURN
        ENDIF

        WRITE(name,'(A,A)') "SNOWFRAC_tile_", TRIM(ADJUSTL(cjt))
        IF(varname == TRIM(name)) THEN
          ptr2 => p_diag_lnd%snowfrac_t(:,:,jt)
          RETURN
        ENDIF

        WRITE(name,'(A,A)') "DZH_SNOW_tile_", TRIM(ADJUSTL(cjt))
        IF(varname == TRIM(name)) THEN
          ptr3 => p_prog_lnd%dzh_snow_t(:,:,:,jt)
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

    TYPE(t_hydro_ocean_prog), POINTER :: p_prog
    TYPE(t_hydro_ocean_diag), POINTER :: p_diag
    TYPE(t_hydro_ocean_aux),  POINTER :: p_aux
    TYPE(t_sfc_flx),          POINTER :: forcing
    TYPE(t_ho_params),        POINTER :: p_params
    TYPE(t_sea_ice),          POINTER :: p_ice
    REAL(wp),                 POINTER :: r_isice(:,:,:)
    REAL(wp),                 POINTER :: r_dolic_c(:,:), r_dolic_e(:,:)
    INTEGER                           :: s_isice(3), shp_dolic(2)

    p_prog  => v_ocean_state(jg)%p_prog(nold(jg))
    p_diag  => v_ocean_state(jg)%p_diag
    p_aux   => v_ocean_state(jg)%p_aux
    forcing => v_sfc_flx 
    p_params=> v_params
    p_ice   => v_sea_ice

    ptr2d     => NULL()
    ptr3d     => NULL()
    reset     = .FALSE.
    delete    = .FALSE.
    not_found = .FALSE.

    SELECT CASE(varname)
      CASE ('wet_c');        ptr3d => v_base%wet_c
      CASE ('wet_e');        ptr3d => v_base%wet_e
      CASE ('rbasin_c');     ptr2d => v_base%rbasin_c(:,:)
      CASE ('rregio_c');     ptr2d => v_base%rregio_c(:,:)
      CASE ('dolic_c')
        shp_dolic = SHAPE(v_base%dolic_c)
        ALLOCATE(r_dolic_c(shp_dolic(1),shp_dolic(2)))
        r_dolic_c(:,:) = REAL(v_base%dolic_c(:,:))
        ptr2d => r_dolic_c
      CASE ('dolic_e')
        shp_dolic = SHAPE(v_base%dolic_e)
        ALLOCATE(r_dolic_e(shp_dolic(1),shp_dolic(2)))
        r_dolic_e(:,:) = REAL(v_base%dolic_e(:,:))
        ptr2d => r_dolic_e
      CASE ('ELEV');         ptr2d => p_prog%h
      CASE ('forc_u');       ptr2d => forcing%forc_wind_u
      CASE ('forc_v');       ptr2d => forcing%forc_wind_v
      CASE ('forc_tdata');   ptr2d => forcing%forc_tracer_relax(:,:,1)
      CASE ('forc_t');       ptr2d => forcing%forc_tracer(:,:,1)
      CASE ('forc_hflx');    ptr2d => forcing%forc_hflx(:,:)
      CASE ('forc_swflx');   ptr2d => forcing%forc_swflx(:,:)
      CASE ('forc_lwflx');   ptr2d => forcing%forc_lwflx(:,:)
      CASE ('forc_ssflx');   ptr2d => forcing%forc_ssflx(:,:)
      CASE ('forc_slflx');   ptr2d => forcing%forc_slflx(:,:)
      CASE ('forc_sdata');   ptr2d => forcing%forc_tracer_relax(:,:,2)
      CASE ('forc_s');       ptr2d => forcing%forc_tracer(:,:,2)
      CASE ('forc_fwfx');    ptr2d => forcing%forc_fwfx(:,:)
      CASE('horz_adv');      ptr3d => p_diag%veloc_adv_horz
      CASE('Ekin_grad');     ptr3d => p_diag%grad
      CASE('Ekin');          ptr3d => p_diag%kin
      CASE ('flux_u');       ptr2d => p_aux%bc_top_u
      CASE ('flux_v');       ptr2d => p_aux%bc_top_v
      CASE ('flux-VN');      ptr2d => p_aux%bc_top_vn
      CASE ('VN');           ptr3d => p_prog%vn
      CASE ('T');            ptr3d => p_prog%tracer(:,:,:,1)
      CASE ('S');            ptr3d => p_prog%tracer(:,:,:,2)
      CASE ('g_n');          ptr3d => p_aux%g_n(:,:,:)
      CASE ('g_nm1');        ptr3d => p_aux%g_nm1(:,:,:)
      CASE ('g_nimd');       ptr3d => p_aux%g_nimd(:,:,:)
      CASE ('g_n_c_h_T');    ptr3d => p_aux%g_n_c_h(:,:,:,1)
      CASE ('g_n_c_h_S');    ptr3d => p_aux%g_n_c_h(:,:,:,2)
      CASE ('g_n_c_v_T');    ptr3d => p_aux%g_n_c_v(:,:,:,1)
      CASE ('g_n_c_v_S');    ptr3d => p_aux%g_n_c_v(:,:,:,2)
      CASE ('g_nm1_c_h_T');  ptr3d => p_aux%g_nm1_c_h(:,:,:,1)
      CASE ('g_nm1_c_h_S');  ptr3d => p_aux%g_nm1_c_h(:,:,:,2)
      CASE ('g_nm1_c_v_T');  ptr3d => p_aux%g_nm1_c_h(:,:,:,1)
      CASE ('g_nm1_c_v_S');  ptr3d => p_aux%g_nm1_c_h(:,:,:,2)
      CASE ('g_nimd_c_h_T'); ptr3d => p_aux%g_nimd_c_v(:,:,:,1)
      CASE ('g_nimd_c_h_S'); ptr3d => p_aux%g_nimd_c_v(:,:,:,2)
      CASE ('g_nimd_c_v_T'); ptr3d => p_aux%g_nimd_c_v(:,:,:,1)
      CASE ('g_nimd_c_v_S'); ptr3d => p_aux%g_nimd_c_v(:,:,:,2)
      CASE ('VORT');         ptr3d => p_diag%vort
      CASE ('u');            ptr3d => p_diag%u
      CASE ('v');            ptr3d => p_diag%v
      CASE ('W');            ptr3d => p_diag%w
      CASE ('u_vint');       ptr2d => p_diag%u_vint
      CASE('Vert_Veloc_Adv');ptr3d => p_diag%veloc_adv_vert
      CASE('press_grad');    ptr3d => p_diag%press_grad
      CASE ('rho');          ptr3d => p_diag%rho   
      CASE('Vert_Mixing_V'); ptr3d => p_params%A_veloc_v
        CASE('Vert_Mixing_T')
          IF (no_tracer > 0) ptr3d => p_params%A_tracer_v(:,:,:,1)
      CASE('Horz_Mixing_V'); ptr3d => p_params%K_veloc_h
      CASE ('cell_owner');   ptr2d => cell_owner
      ! sea ice variables
      CASE('p_ice_isice')
        s_isice = SHAPE(p_ice%isice)
        ALLOCATE(r_isice(s_isice(1),s_isice(2),s_isice(3)))
        WHERE(p_ice%isice)
          r_isice = 1.0_wp
        ELSEWHERE
          r_isice = 0.0_wp
        ENDWHERE
        ptr3d => r_isice
      CASE('p_ice_alb');         ptr3d => p_ice%alb
      CASE('p_ice_Tsurf');       ptr3d => p_ice%Tsurf
      CASE('p_ice_T1');          ptr3d => p_ice%T1
      CASE('p_ice_T2');          ptr3d => p_ice%T2
      CASE('p_ice_E1');          ptr3d => p_ice%E1
      CASE('p_ice_E2');          ptr3d => p_ice%E2
      CASE('p_ice_hi');          ptr3d => p_ice%hi
      CASE('p_ice_hs');          ptr3d => p_ice%hs
      CASE('p_ice_Qtop');        ptr3d => p_ice%Qtop
      CASE('p_ice_Qbot');        ptr3d => p_ice%Qbot
      CASE('p_ice_heatocei');    ptr3d => p_ice%heatocei
      CASE('p_ice_snow_to_ice'); ptr3d => p_ice%snow_to_ice
      CASE('p_ice_surfmelt');    ptr3d => p_ice%surfmelt
      CASE('p_ice_surfmeltT');   ptr3d => p_ice%surfmeltT
      CASE('p_ice_evapwi');      ptr3d => p_ice%evapwi
      CASE('p_ice_conc');        ptr3d => p_ice%conc
      CASE('p_ice_u');           ptr2d => p_ice%u
      CASE('p_ice_v');           ptr2d => p_ice%v
      CASE('p_ice_concSum');     ptr2d => p_ice%concSum
      CASE('p_ice_newice');      ptr2d => p_ice%newice
      CASE('p_ice_zUnderIce');   ptr2d => p_ice%zUnderIce

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
    !   DO jt = 1, nqtendphy ! only hydromet-tendencies (qv,qc,qi)
    !     ctracer = ctracer_list(jt:jt)
    !     IF(varname == 'tend_q'//ctracer//'_conv') THEN
    !       ptr3 => prm_nwp_tend(jg)%ddt_tracer_pconv(:,:,:,jt)
    !       RETURN
    !     ENDIF
    !     IF(varname == 'tend_q'//ctracer//'_turb') THEN
    !       ptr3 => prm_nwp_tend(jg)%ddt_tracer_turb(:,:,:,jt)
    !       RETURN
    !     ENDIF
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

        ! Get a pointer to the variable
        SELECT CASE (iequations)
          CASE (ishallow_water, ihs_atm_temp, ihs_atm_theta)
            CALL get_outvar_ptr_ha(outvar_desc(ivar,jg)%name, jg, ptr2, ptr3, reset, delete)
          CASE (inh_atmosphere)
            CALL get_outvar_ptr_nh &
              & (outvar_desc(ivar,jg)%name, jg, ptr2, ptr3, reset, delete)
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

    num_varids(k_jg)      = num_varids(k_jg) + 1
    ivar_idx              = num_varids(k_jg)
    varids(ivar_idx,k_jg) = var

  END SUBROUTINE addVar


  !>
  !! Utility function: Translates NetCDF name to internal variable name.
  !!
  !! In the "traditional" vlist output mode, variable names appearing
  !! in the NetCDF output were defined by addVar commands. For
  !! example, the temperature was named "T". In contrast to this, the
  !! new output mode (MODULE mo_name_list_output) uses the
  !! varlistelement%field%info%name string as output name, which is
  !! set by add_var. 
  !! A name dictionary has been implemented by Rainer Johanni in
  !! r7514. This subroutine serves the purpose to create such a
  !! translation table automatically.
  !!
  !! @par Revision History
  !! Initial implementation by F. Prill, DWD (2011-12-22)
  !!
  SUBROUTINE translate_vars(k_jg)
    INTEGER, INTENT(IN)           :: k_jg
    ! local variables
    TYPE(t_list_element), POINTER :: list_element
    INTEGER                       :: i, nindex, idx, ivar
    REAL(wp),             POINTER :: ptr2(:,:)
    REAL(wp),             POINTER :: ptr3(:,:,:)
    LOGICAL                       :: reset, delete

    ! loop over io_vlist variables:
    DO ivar = 1, num_output_vars(k_jg)
      ! get pointer to variable
      CALL get_outvar_ptr_nh &
        & (outvar_desc(ivar,k_jg)%name, k_jg, ptr2, ptr3, reset, delete)
      ! loop over "add_var" variable list and compare to pointer
      DO i = 1, nvar_lists
        list_element => var_lists(i)%p%first_list_element
        DO WHILE (ASSOCIATED(list_element))
          IF (list_element%field%info%lcontained) THEN 
            nindex = list_element%field%info%ncontained
          ELSE
            nindex = 1
          ENDIF
          idx = INDEX(list_element%field%info%name,'.TL')
          
          IF (ASSOCIATED(ptr2,list_element%field%r_ptr(:,:,nindex,1,1)) .OR.  &
            & ASSOCIATED(ptr3,list_element%field%r_ptr(:,:,:,nindex,1))) THEN
            
            IF(idx == 0) THEN
              WRITE (0,*) TRIM(list_element%field%info%name), "  ", &
                & outvar_desc(ivar,k_jg)%name
            ELSE
              WRITE (0,*) TRIM(list_element%field%info%name(1:idx-1)), "  ", &
                & outvar_desc(ivar,k_jg)%name
            END IF
          END IF
          list_element => list_element%next_list_element
        END DO ! while
      END DO !i
    END DO ! ivar

  END SUBROUTINE translate_vars


  !-------------------------------------------------------------------------------------------------
  FUNCTION TimeVar(vname,vlongname,vunit,vcode,vtable,vlist,grid,zaxis) RESULT(var)
    INTEGER                  :: var
    CHARACTER(*), INTENT(IN) :: vname, vlongname, vunit
    INTEGER, INTENT(IN)      :: vcode, vtable, vlist, grid, zaxis

    INTEGER :: vartype
    
    vartype = TIME_VARIABLE
    
    var = vlistdefvar(vlist, grid, zaxis, vartype)
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
    CALL addGlobAttFlt('run_nml:dtime_adv',dtime_adv,vlist,astatus)
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
             IF (TRIM(ape_sst_case) == 'ape_sst_const') THEN
               CALL addGlobAttFlt('nh_testcase_nml:ape_sst_val',           &
                 &              ape_sst_val,vlist,astatus)
             END IF
          !
          ELSEIF ( TRIM(nh_test_name) == 'dcmip_tc_52') THEN
             CALL addGlobAttTxt('nh_testcase_nml:ape_sst_case',          &
                 &              TRIM(ape_sst_case),vlist,astatus)
             CALL addGlobAttFlt('nh_testcase_nml:ape_sst_val',           &
                 &              ape_sst_val,vlist,astatus)
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

END MODULE mo_io_vlist
