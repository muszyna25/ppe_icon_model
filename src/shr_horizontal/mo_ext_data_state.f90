!>
!! Allocation/deallocation of external parameter state
!!
!! This module contains routines for setting up the external data state.
!!
!! @author Daniel Reinert, DWD
!! @author Hermann Asensio, DWD
!!
!!
!! @par Revision History
!! Initial revision by Daniel Reinert, DWD (2010-07-12)
!! Modification by Hermann Asensio, DWD (2010-07-16)
!!  - add miscellaneous variables for external parameters
!! Modification by Daniel Reinert, DWD (2011-05-03)
!! - Memory allocation method changed from explicit allocation to Luis'
!!   infrastructure
!! Modification by Daniel Reinert, DWD (2012-02-23)
!! - Routine smooth_topography moved to a new module named mo_smooth_topo
!! Modification by Daniel Reinert, DWD (2012-03-22)
!! - Type declaration moved to new module mo_ext_data_types
!! Modification by Daniel Reinert, DWD (2015-06-11)
!! - reading and other stuff moved to new module mo_ext_data_init
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

MODULE mo_ext_data_state

  USE mo_kind,               ONLY: wp
  USE mo_impl_constants,     ONLY: inwp, iecham, MODIS, ildf_echam,                &
    &                              ihs_atm_temp, ihs_atm_theta, io3_clim, io3_ape, &
    &                              HINTP_TYPE_LONLAT_NNB, MAX_CHAR_LENGTH,         &
    &                              SSTICE_ANA, SSTICE_ANA_CLINC, SSTICE_CLIM,      &
    &                              SSTICE_AVG_MONTHLY, SSTICE_AVG_DAILY
  USE mo_cdi_constants,      ONLY: GRID_UNSTRUCTURED_CELL, GRID_CELL
  USE mo_exception,          ONLY: message, finish
  USE mo_model_domain,       ONLY: t_patch
  USE mo_ext_data_types,     ONLY: t_external_data, t_external_atmos_td, &
    &                              t_external_atmos
  USE mo_linked_list,        ONLY: t_var_list
  USE mo_var_groups,         ONLY: groups
  USE mo_var_metadata_types, ONLY: POST_OP_SCALE, POST_OP_LUC, CLASS_TILE
  USE mo_var_metadata,       ONLY: post_op, create_hor_interp_metadata
  USE mo_var_list,           ONLY: new_var_list, delete_var_list, add_var, add_ref, &
    &                              default_var_list_settings
  USE mo_cf_convention,      ONLY: t_cf_var
  USE mo_grib2,              ONLY: t_grib2_var, grib2_var, t_grib2_int_key, &
    &                              OPERATOR(+)
  USE mo_parallel_config,    ONLY: nproma
  USE mo_io_config,          ONLY: lnetcdf_flt64_output
  USE mo_grid_config,        ONLY: n_dom
  USE mo_run_config,         ONLY: iforcing
  USE mo_dynamics_config,    ONLY: iequations
  USE mo_lnd_nwp_config,     ONLY: ntiles_total, ntiles_water, llake, &
    &                              sstice_mode
  USE mo_radiation_config,   ONLY: irad_o3, albedo_type
  USE mo_extpar_config,      ONLY: i_lctype, nclass_lu, nmonths_ext, itype_vegetation_cycle
  USE mo_cdi,                ONLY: DATATYPE_PACK16, DATATYPE_FLT32, DATATYPE_FLT64, &
    &                              TSTEP_CONSTANT, TSTEP_MAX, TSTEP_AVG,            &
    &                              GRID_UNSTRUCTURED
  USE mo_zaxis_type,         ONLY: ZA_REFERENCE, ZA_LAKE_BOTTOM, ZA_SURFACE, &
    &                              ZA_HEIGHT_2M, ZA_PRESSURE

  IMPLICIT NONE


  PRIVATE

  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_ext_data_state'

  ! necessary information when reading ozone from file
  CHARACTER(len=6)  :: levelname
  CHARACTER(len=6)  :: cellname
  CHARACTER(len=5)  :: o3name
  CHARACTER(len=20) :: o3unit
  !
  INTEGER :: nlev_o3
  INTEGER :: nmonths

  ! variables
  PUBLIC :: nmonths
  PUBLIC :: nlev_o3
  PUBLIC :: levelname
  PUBLIC :: cellname
  PUBLIC :: o3name
  PUBLIC :: o3unit

  ! state
  PUBLIC :: ext_data

  ! subroutines
  PUBLIC :: construct_ext_data
  PUBLIC :: destruct_ext_data


  TYPE(t_external_data),TARGET, ALLOCATABLE :: &
    &  ext_data(:)  ! n_dom

!-------------------------------------------------------------------------

CONTAINS



  !-------------------------------------------------------------------------
  !>
  !! Top-level procedure for building external data structure
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2010-07-12)
  !!
  SUBROUTINE construct_ext_data (p_patch, ext_data)

    TYPE(t_patch),          INTENT(IN)    :: p_patch(:)
    TYPE(t_external_data),  INTENT(INOUT) :: ext_data(:)

    INTEGER :: jg

    CHARACTER(len=MAX_CHAR_LENGTH) :: listname

    CHARACTER(len=max_char_length), PARAMETER :: &
      routine = modname//':construct_ext_data'

!-------------------------------------------------------------------------


    CALL message (TRIM(routine), 'Construction of data structure for ' // &
      &                          'external data started')

    ! Build external data list for constant-in-time fields for the atm model
    DO jg = 1, n_dom
      WRITE(listname,'(a,i2.2)') 'ext_data_atm_D',jg
      CALL new_ext_data_atm_list(p_patch(jg), ext_data(jg)%atm,       &
        &                        ext_data(jg)%atm_list, TRIM(listname))
    END DO

    ! Build external data list for time-dependent fields
    IF (iforcing > 1 ) THEN ! further distinction is made inside
      DO jg = 1, n_dom
        WRITE(listname,'(a,i2.2)') 'ext_data_atm_td_D',jg
        CALL new_ext_data_atm_td_list(p_patch(jg), ext_data(jg)%atm_td,       &
          &                           ext_data(jg)%atm_td_list, TRIM(listname))
      END DO
    END IF

    CALL message (TRIM(routine), 'Construction of data structure for ' // &
      &                          'external data finished')

  END SUBROUTINE construct_ext_data

  !-------------------------------------------------------------------------
  !
  !
  !>
  !! Allocation of atmospheric external data structure
  !!
  !! Allocation of atmospheric external data structure (constant in time
  !! elements).
  !!
  !! Initialization of elements with zero.
  !!
  !! @par Revision History
  !! Initial release by Daniel Reinert (2011-05-03)
  !! Statements that assign initial value added by Hui Wan (MPI-M, 2011-05-30)
  !!
  SUBROUTINE new_ext_data_atm_list ( p_patch, p_ext_atm, p_ext_atm_list, &
    &                                listname)
!
    TYPE(t_patch), TARGET , INTENT(IN)   :: & !< current patch
      &  p_patch

    TYPE(t_external_atmos), INTENT(INOUT):: & !< current external data structure
      &  p_ext_atm

    TYPE(t_var_list)      , INTENT(INOUT):: p_ext_atm_list !< current external data list

    CHARACTER(len=*)      , INTENT(IN)   :: & !< list name
      &  listname

    TYPE(t_cf_var)    :: cf_desc, new_cf_desc
    TYPE(t_grib2_var) :: grib2_desc

    INTEGER :: jg

    INTEGER :: nlev          !< number of vertical levels

    INTEGER :: nblks_c       !< number of cell blocks to allocate

    INTEGER :: shape2d_c(2)
    INTEGER :: shape3d_c(3)
    INTEGER :: shape3d_sfc(3), shape3d_nt(3), shape3d_ntw(3)

    INTEGER :: ibits         !< "entropy" of horizontal slice
    INTEGER :: datatype_flt  !< floating point accuracy in NetCDF output

    INTEGER          :: jsfc
    CHARACTER(LEN=2) :: csfc

    !--------------------------------------------------------------

    !determine size of arrays
    nblks_c = p_patch%nblks_c

    ! get patch ID
    jg = p_patch%id
    ibits = DATATYPE_PACK16   ! "entropy" of horizontal slice

    ! number of vertical levels
    nlev = p_patch%nlev

    IF ( lnetcdf_flt64_output ) THEN
      datatype_flt = DATATYPE_FLT64
    ELSE
      datatype_flt = DATATYPE_FLT32
    ENDIF

    ! predefined array shapes
    shape2d_c  = (/ nproma, nblks_c /)
    shape3d_c  = (/ nproma, nlev, nblks_c       /)
    shape3d_sfc= (/ nproma, nblks_c, nclass_lu(jg) /)
    shape3d_nt = (/ nproma, nblks_c, ntiles_total     /)
    shape3d_ntw = (/ nproma, nblks_c, ntiles_total + ntiles_water /)


    !------------------------------
    ! Ensure that all pointers have a defined association status.
    !------------------------------
    NULLIFY(p_ext_atm%topography_c,    &
      &     p_ext_atm%grad_topo,       &
      &     p_ext_atm%topo_t2mclim,    &
      &     p_ext_atm%fis,             &
      &     p_ext_atm%o3,              &
      &     p_ext_atm%llsm_atm_c,      &
      &     p_ext_atm%llake_c,         &
      &     p_ext_atm%fr_land,         &
      &     p_ext_atm%fr_glac,         &
      &     p_ext_atm%fr_land_smt,     &
      &     p_ext_atm%fr_glac_smt,     &
      &     p_ext_atm%z0,              &
      &     p_ext_atm%fr_lake,         &
      &     p_ext_atm%depth_lk,        &
      &     p_ext_atm%fetch_lk,        &
      &     p_ext_atm%dp_bs_lk,        &
      &     p_ext_atm%t_bs_lk,         &
      &     p_ext_atm%gamso_lk,        &
      &     p_ext_atm%sso_stdh,        &
      &     p_ext_atm%sso_stdh_raw,    &
      &     p_ext_atm%sso_gamma,       &
      &     p_ext_atm%sso_theta,       &
      &     p_ext_atm%sso_sigma,       &
      &     p_ext_atm%plcov_mx,        &
      &     p_ext_atm%plcov,           &
      &     p_ext_atm%plcov_t,         &
      &     p_ext_atm%lai_mx,          &
      &     p_ext_atm%lai,             &
      &     p_ext_atm%sai,             &
      &     p_ext_atm%sai_t,           &
      &     p_ext_atm%tai,             &
      &     p_ext_atm%tai_t,           &
      &     p_ext_atm%laifac_t,        &
      &     p_ext_atm%eai,             &
      &     p_ext_atm%eai_t,           &
      &     p_ext_atm%rootdp,          &
      &     p_ext_atm%rootdp_t,        &
      &     p_ext_atm%for_e,           &
      &     p_ext_atm%for_d,           &
      &     p_ext_atm%rsmin,           &
      &     p_ext_atm%rsmin2d_t,       &
      &     p_ext_atm%ndvi_max,        &
      &     p_ext_atm%ndviratio,       &
      &     p_ext_atm%idx_lst_lp,      &
      &     p_ext_atm%idx_lst_sp,      &
      &     p_ext_atm%idx_lst_fp,      &
      &     p_ext_atm%idx_lst_lp_t,    &
      &     p_ext_atm%idx_lst_t,       &
      &     p_ext_atm%idx_lst_spw,     &
      &     p_ext_atm%idx_lst_spi,     &
      &     p_ext_atm%snowtile_flag_t, &
      &     p_ext_atm%lc_class_t,      &
      &     p_ext_atm%lc_frac_t,       &
      &     p_ext_atm%frac_t,          &
      &     p_ext_atm%inv_frland_from_tiles, &
      &     p_ext_atm%soiltyp,         &
      &     p_ext_atm%soiltyp_t,       &
      &     p_ext_atm%t_cl,            &
      &     p_ext_atm%emis_rad,        &
      &     p_ext_atm%lu_class_fraction, &
      &     p_ext_atm%alb_dif,         &
      &     p_ext_atm%albuv_dif,       &
      &     p_ext_atm%albni_dif,       &
      &     p_ext_atm%lsm_ctr_c,       &
      &     p_ext_atm%elevation_c,     &
      &     p_ext_atm%emis_rad         )


    !
    ! Register a field list and apply default settings
    !
    CALL new_var_list( p_ext_atm_list, TRIM(listname), patch_id=p_patch%id )
    CALL default_var_list_settings( p_ext_atm_list,            &
                                  & lrestart=.FALSE.  )


    ! topography height at cell center
    !
    ! topography_c  p_ext_atm%topography_c(nproma,nblks_c)
    cf_desc    = t_cf_var('surface_height', 'm', &
      &                   'geometric height of the earths surface above sea level', datatype_flt)
    grib2_desc = grib2_var( 0, 3, 6, ibits, GRID_UNSTRUCTURED, GRID_CELL)  &
      &           + t_grib2_int_key("typeOfSecondFixedSurface", 101)
    CALL add_var( p_ext_atm_list, 'topography_c', p_ext_atm%topography_c,  &
      &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,             &
      &           grib2_desc, ldims=shape2d_c, loutput=.TRUE.,             &
      &           isteptype=TSTEP_CONSTANT )


    ! gradient of topography height at cell center
    !
    ! grad_topo     p_ext_atm%grad_topo(2,nproma,nblks_c)
    cf_desc    = t_cf_var('grad_surface_height', 'm m-1', &
      &                   'gradient of geometric height of the earths surface above sea level', datatype_flt)
    grib2_desc = grib2_var( 255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_ext_atm_list, 'grad_topo', p_ext_atm%grad_topo,        &
      &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,             &
      &           grib2_desc, ldims=(/2,nproma,nblks_c/), loutput=.FALSE.,  &
      &           isteptype=TSTEP_CONSTANT )


    IF (itype_vegetation_cycle > 1) THEN
      ! interpolated topographic height for T2M climatology data
      !
      ! topo_t2mclim     p_ext_atm%topo_t2mclim(nproma,nblks_c)
      cf_desc    = t_cf_var('surface_height_of_T2M_climatology', 'm', &
        &                   'interpolated topographic height for T2M climatology data', datatype_flt)
      grib2_desc = grib2_var( 255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'topo_t2mclim', p_ext_atm%topo_t2mclim,  &
        &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,             &
        &           grib2_desc, ldims=shape2d_c, loutput=.TRUE.,             &
        &           isteptype=TSTEP_CONSTANT )
    ENDIF

    ! geopotential (s)
    !
    ! fis          p_ext_atm%fis(nproma,nblks_c)
    cf_desc    = t_cf_var('Geopotential_(s)', 'm2 s-2', &
      &                   'Geopotential (s)', datatype_flt)
    grib2_desc = grib2_var( 0, 3, 4, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_ext_atm_list, 'fis', p_ext_atm%fis,           &
      &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,    &
      &           grib2_desc, ldims=shape2d_c, loutput=.TRUE.,    &
      &           isteptype=TSTEP_CONSTANT )


    IF ( iforcing == inwp ) THEN

      ! ozone mixing ratio
      !
      ! o3            p_ext_atm%o3(nproma,nlev,nblks_c)
      cf_desc    = t_cf_var('ozone mixing ratio', 'kg kg-1', &
        &                   'ozone mixing ratio', datatype_flt)
      grib2_desc = grib2_var( 0, 14, 1, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'o3', p_ext_atm%o3,                      &
        &           GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc,              &
        &           grib2_desc, ldims=shape3d_c, loutput=.TRUE. )

      ! external parameter for NWP forcing

      ! land sea mask for cells (LOGICAL)
      ! Note: Here "loutput" is set to .FALSE. since the output
      !       scheme operates on REAL model variables only and
      !       throws an error on this.
      !
      ! llsm_atm_c    p_ext_atm%llsm_atm_c(nproma,nblks_c)
      cf_desc    = t_cf_var('land_sea_mask_(cell)', '-', &
        &                   'land sea mask (cell)', datatype_flt)
      grib2_desc = grib2_var( 2, 0, 0, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'llsm_atm_c', p_ext_atm%llsm_atm_c, &
        &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,        &
        &           grib2_desc, ldims=shape2d_c, loutput=.FALSE.,       &
        &           isteptype=TSTEP_CONSTANT )

      ! llake_c    p_ext_atm%llake_c(nproma,nblks_c)
      cf_desc    = t_cf_var('lake_mask_(cell)', '-', &
        &                   'lake mask (cell)', datatype_flt)
      grib2_desc = grib2_var( 2, 0, 0, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'llake_c', p_ext_atm%llake_c,       &
        &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,        &
        &           grib2_desc, ldims=shape2d_c, loutput=.FALSE.,       &
        &           isteptype=TSTEP_CONSTANT )


      ! land fraction
      !
      ! fr_land      p_ext_atm%fr_land(nproma,nblks_c)
      cf_desc    = t_cf_var('land_area_fraction', '-', 'Fraction land', datatype_flt)
      grib2_desc = grib2_var( 2, 0, 0, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'fr_land', p_ext_atm%fr_land,   &
        &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,    &
        &           grib2_desc, ldims=shape2d_c, loutput=.TRUE.,    &
        &           isteptype=TSTEP_CONSTANT,                       &
        &           in_group=groups("dwd_fg_sfc_vars","ICON_INI_OUT") )


      ! glacier fraction
      !
      ! fr_glac      p_ext_atm%fr_glac(nproma,nblks_c)
      cf_desc    = t_cf_var('glacier_area_fraction', '-', 'Fraction glacier', datatype_flt)
      grib2_desc = grib2_var( 2, 0, 192, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'fr_glac', p_ext_atm%fr_glac,   &
        &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,    &
        &           grib2_desc, ldims=shape2d_c, loutput=.TRUE. )


      ! maybe the next one (fr_land_smt)
      ! should be moved into corresponding if block

      ! land fraction (smoothed)
      !
      ! fr_land_smt  p_ext_atm%fr_land_smt(nproma,nblks_c)
      cf_desc    = t_cf_var('land_area_fraction_(smoothed)', '-', &
        &                   'land area fraction (smoothed)', datatype_flt)
      grib2_desc = grib2_var( 2, 0, 0, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'fr_land_smt', p_ext_atm%fr_land_smt, &
        &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,          &
        &           grib2_desc, ldims=shape2d_c, loutput=.FALSE.,         &
        &           isteptype=TSTEP_CONSTANT )


      ! glacier area fraction (smoothed)
      !
      ! fr_glac_smt  p_ext_atm%fr_glac_smt(nproma,nblks_c)
      cf_desc    = t_cf_var('glacier_area_fraction_(smoothed)', '-', &
        &                   'glacier area fraction (smoothed)', datatype_flt)
      grib2_desc = grib2_var( 2, 0, 192, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'fr_glac_smt', p_ext_atm%fr_glac_smt, &
        &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,          &
        &           grib2_desc, ldims=shape2d_c, loutput=.FALSE. )

      ! roughness length
      !
      ! z0           p_ext_atm%z0(nproma,nblks_c)
      cf_desc    = t_cf_var('roughtness_length', 'm', 'roughtness length', datatype_flt)
      grib2_desc = grib2_var( 2, 0, 1, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'z0', p_ext_atm%z0,             &
        &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,    &
        &           grib2_desc, ldims=shape2d_c, loutput=.FALSE. )

      !
      ! fr_lake and lake depth are needed, even if the lake model is switched off
      !

      ! fraction lake
      !
      ! fr_lake      p_ext_atm%fr_lake(nproma,nblks_c)
      cf_desc    = t_cf_var('fraction_lake', '-', 'fraction lake', datatype_flt)
      grib2_desc = grib2_var( 1, 2, 2, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'fr_lake', p_ext_atm%fr_lake,   &
        &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,    &
        &           grib2_desc, ldims=shape2d_c, loutput=.TRUE.,   &
        &           isteptype=TSTEP_CONSTANT )


      ! lake depth
      !
      ! depth_lk     p_ext_atm%depth_lk(nproma,nblks_c)
      cf_desc    = t_cf_var('lake_depth', 'm', 'lake depth', datatype_flt)
      grib2_desc = grib2_var( 1, 2, 0, ibits, GRID_UNSTRUCTURED, GRID_CELL)  &
        &           + t_grib2_int_key("typeOfFirstFixedSurface", 1)
      CALL add_var( p_ext_atm_list, 'depth_lk', p_ext_atm%depth_lk, &
        &           GRID_UNSTRUCTURED_CELL, ZA_LAKE_BOTTOM, cf_desc,&
        &           grib2_desc, ldims=shape2d_c, loutput=.TRUE.,    &
        &           isteptype=TSTEP_CONSTANT )

      IF (llake) THEN

        ! fetch_lk     p_ext_atm%fetch_lk(nproma,nblks_c)
        cf_desc    = t_cf_var('fetch_lk', 'm', 'wind fetch over lake', datatype_flt)
        grib2_desc = grib2_var( 0, 2, 33, ibits, GRID_UNSTRUCTURED, GRID_CELL)
        CALL add_var( p_ext_atm_list, 'fetch_lk', p_ext_atm%fetch_lk, &
          &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,    &
          &           grib2_desc, ldims=shape2d_c, loutput=.TRUE.,    &
          &           isteptype=TSTEP_CONSTANT )


        ! dp_bs_lk     p_ext_atm%dp_bs_lk(nproma,nblks_c)
        cf_desc    = t_cf_var('dp_bs_lk', 'm', &
          &          'depth of thermally active layer of bot. sediments.', datatype_flt)
        grib2_desc = grib2_var( 1, 2, 3, ibits, GRID_UNSTRUCTURED, GRID_CELL)
        CALL add_var( p_ext_atm_list, 'dp_bs_lk', p_ext_atm%dp_bs_lk, &
          &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,    &
          &           grib2_desc, ldims=shape2d_c, loutput=.TRUE.,    &
          &           isteptype=TSTEP_CONSTANT )


        ! t_bs_lk     p_ext_atm%t_bs_lk(nproma,nblks_c)
        cf_desc    = t_cf_var('t_bs_lk', 'm', &
          &          'clim. temp. at bottom of thermally active layer of sediments', &
          &          datatype_flt)
        grib2_desc = grib2_var( 1, 2, 4, ibits, GRID_UNSTRUCTURED, GRID_CELL)
        CALL add_var( p_ext_atm_list, 't_bs_lk', p_ext_atm%t_bs_lk,   &
          &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,    &
          &           grib2_desc, ldims=shape2d_c, loutput=.TRUE.,    &
          &           isteptype=TSTEP_CONSTANT )


        ! gamso_lk     p_ext_atm%gamso_lk(nproma,nblks_c)
        cf_desc    = t_cf_var('gamso_lk', 'm', &
          &          'attenuation coefficient of lake water with respect to sol. rad.', &
          &          datatype_flt)
        grib2_desc = grib2_var( 1, 2, 11, ibits, GRID_UNSTRUCTURED, GRID_CELL)
        CALL add_var( p_ext_atm_list, 'gamso_lk', p_ext_atm%gamso_lk, &
          &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,    &
          &           grib2_desc, ldims=shape2d_c, loutput=.TRUE.,    &
          &           isteptype=TSTEP_CONSTANT )

      ENDIF



      !--------------------------------
      ! sub-gridscale orography
      !--------------------------------

      ! Standard deviation of sub-grid scale orography
      !
      ! sso_stdh     p_ext_atm%sso_stdh(nproma,nblks_c)
      cf_desc    = t_cf_var('standard_deviation_of_height', 'm',    &
        &                   'Standard deviation of sub-grid scale orography', datatype_flt)
      grib2_desc = grib2_var( 0, 3, 20, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'sso_stdh', p_ext_atm%sso_stdh, &
        &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,    &
        &           grib2_desc, ldims=shape2d_c, loutput=.TRUE.,    &
        &           isteptype=TSTEP_CONSTANT )


      ! field derived from sso_stdh used for pat_len in turbulence scheme
      ! for the time being, it is the same as sso_stdh except for not being adjusted to orography smoothing
      !
      ! sso_stdh_raw     p_ext_atm%sso_stdh_raw(nproma,nblks_c)
      cf_desc    = t_cf_var('standard_deviation_of_height', 'm',    &
        &                   'Standard deviation of sub-grid scale orography', datatype_flt)
      grib2_desc = grib2_var( 255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'sso_stdh_raw', p_ext_atm%sso_stdh_raw, &
        &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,    &
        &           grib2_desc, ldims=shape2d_c, loutput=.TRUE.,    &
        &           isteptype=TSTEP_CONSTANT )


      ! Anisotropy of sub-gridscale orography
      !
      ! sso_gamma    p_ext_atm%sso_gamma(nproma,nblks_c)
      cf_desc    = t_cf_var('anisotropy_factor', '-',&
        &                   'Anisotropy of sub-gridscale orography', datatype_flt)
      grib2_desc = grib2_var( 0, 3, 24, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'sso_gamma', p_ext_atm%sso_gamma, &
        &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,      &
        &           grib2_desc, ldims=shape2d_c, loutput=.TRUE.,      &
        &           isteptype=TSTEP_CONSTANT )



      ! Angle of sub-gridscale orography
      !
      ! sso_theta    p_ext_atm%sso_theta(nproma,nblks_c)
      cf_desc    = t_cf_var('angle_of_principal_axis', 'radians',&
        &                   'Angle of sub-gridscale orography', datatype_flt)
      grib2_desc = grib2_var( 0, 3, 21, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'sso_theta', p_ext_atm%sso_theta, &
        &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,      &
        &           grib2_desc, ldims=shape2d_c, loutput=.TRUE.,      &
        &           isteptype=TSTEP_CONSTANT )



      ! Slope of sub-gridscale orography
      !
      ! sso_sigma    p_ext_atm%sso_sigma(nproma,nblks_c)
      cf_desc    = t_cf_var('slope_of_terrain', '-',&
        &                   'Slope of sub-gridscale orography', datatype_flt)
      grib2_desc = grib2_var( 0, 3, 22, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'sso_sigma', p_ext_atm%sso_sigma, &
        &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,      &
        &           grib2_desc, ldims=shape2d_c, loutput=.TRUE.,      &
        &           isteptype=TSTEP_CONSTANT )





      !--------------------------------
      ! vegetation parameters
      !--------------------------------

      ! Plant covering degree in the vegetation phase
      !
      ! plcov_mx     p_ext_atm%plcov_mx(nproma,nblks_c)
      cf_desc    = t_cf_var('vegetation_area_fraction_vegetation_period', '-',&
        &                   'Plant covering degree in the vegetation phase', datatype_flt)
      grib2_desc = grib2_var( 2, 0, 4, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'plcov_mx', p_ext_atm%plcov_mx, &
        &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,    &
        &           grib2_desc, ldims=shape2d_c, loutput=.FALSE.,   &
        &           isteptype=TSTEP_MAX )

      ! plcov     p_ext_atm%plcov(nproma,nblks_c)
      cf_desc    = t_cf_var('vegetation_area_fraction_vegetation_period', '-',&
        &                   'Plant covering degree in the vegetation phase', datatype_flt)
      new_cf_desc= t_cf_var('vegetation_area_fraction_vegetation_period', '%',&
        &                   'Plant covering degree in the vegetation phase', datatype_flt)
      grib2_desc = grib2_var( 2, 0, 4, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'plcov', p_ext_atm%plcov,       &
        &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,    &
        &           grib2_desc, ldims=shape2d_c, loutput=.TRUE.,    &
        &           isteptype=TSTEP_CONSTANT,                       &
        &           post_op=post_op(POST_OP_SCALE, arg1=100._wp,    &
        &                 new_cf=new_cf_desc) )


      ! plcov_t     p_ext_atm%plcov_t(nproma,nblks_c,ntiles_total)
      cf_desc    = t_cf_var('vegetation_area_fraction_vegetation_period', '-',&
        &                   'Plant covering degree in the vegetation phase', datatype_flt)
      grib2_desc = grib2_var( 2, 0, 4, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'plcov_t', p_ext_atm%plcov_t,    &
        &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,     &
        &           grib2_desc, ldims=shape3d_nt, lcontainer=.TRUE., &
        &           loutput=.FALSE. )

      ALLOCATE(p_ext_atm%plcov_t_ptr(ntiles_total))
      DO jsfc = 1,ntiles_total
        WRITE(csfc,'(i2)') jsfc
        CALL add_ref( p_ext_atm_list, 'plcov_t',                         &
               & 'plcov_t_'//ADJUSTL(TRIM(csfc)),                        &
               & p_ext_atm%plcov_t_ptr(jsfc)%p_2d,                       &
               & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                     &
               & t_cf_var('plcov_t_'//csfc, '', '', datatype_flt),     &
               & grib2_desc,                                             &
               & ldims=shape2d_c, loutput=.TRUE.)
      ENDDO



      ! Max Leaf area index
      !
      ! lai_mx       p_ext_atm%lai_mx(nproma,nblks_c)
      cf_desc    = t_cf_var('leaf_area_index_vegetation_period', '-',&
        &                   'Leaf Area Index Maximum', datatype_flt)
      grib2_desc = grib2_var( 2, 0, 28, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'lai_mx', p_ext_atm%lai_mx,     &
        &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,    &
        &           grib2_desc, ldims=shape2d_c, loutput=.FALSE.,   &
        &           isteptype=TSTEP_MAX )

      ! Leaf area index (aggregated)
      !
      ! lai       p_ext_atm%lai(nproma,nblks_c)
      cf_desc    = t_cf_var('leaf_area_index_vegetation_period', '-',&
        &                   'Leaf Area Index', datatype_flt)
      grib2_desc = grib2_var( 2, 0, 28, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'lai', p_ext_atm%lai,           &
        &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,    &
        &           grib2_desc, ldims=shape2d_c, loutput=.TRUE.,   &
        &           isteptype=TSTEP_CONSTANT )

      ! Surface area index (aggregated)
      !
      ! sai        p_ext_atm%sai(nproma,nblks_c)
      cf_desc    = t_cf_var('sai', ' ','surface area index', datatype_flt)
      grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'sai', p_ext_atm%sai,            &
        &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,     &
        &           grib2_desc, ldims=shape2d_c, loutput=.TRUE.)

      ! Surface area index
      !
      ! sai_t       p_ext_atm%sai_t(nproma,nblks_c,ntiles_total+ntiles_water)
      cf_desc    = t_cf_var('surface_area_index_vegetation_period', '-',&
        &                   'Surface Area Index', datatype_flt)
      grib2_desc = grib2_var( 255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'sai_t', p_ext_atm%sai_t,     &
        &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,  &
        &           grib2_desc, ldims=shape3d_ntw, loutput=.FALSE. )

      ! Transpiration area index (aggregated)
      !
      ! tai         p_ext_atm%tai(nproma,nblks_c)
      cf_desc    = t_cf_var('tai', ' ','transpiration area index', datatype_flt)
      grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'tai', p_ext_atm%tai,         &
        &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,  &
        &           grib2_desc, ldims=shape2d_c, loutput=.TRUE.)

      ! Transpiration area index
      !
      ! tai_t       p_ext_atm%tai_t(nproma,nblks_c,ntiles_total)
      cf_desc    = t_cf_var('transpiration_area_index_vegetation_period', '-',&
        &                   'Transpiration Area Index', datatype_flt)
      grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'tai_t', p_ext_atm%tai_t,     &
        &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,  &
        &           grib2_desc, ldims=shape3d_nt, loutput=.FALSE. )

      ! ratio between current LAI and laimax
      !
      ! laifac_t       p_ext_atm%laifac_t(nproma,nblks_c,ntiles_total)
      cf_desc    = t_cf_var('lai_ratio', '-',&
        &                   'ratio between current LAI and laimax', datatype_flt)
      grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'laifac_t', p_ext_atm%laifac_t,&
        &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,  &
        &           grib2_desc, ldims=shape3d_nt, loutput=.FALSE. )


      ! Evaporative area index (aggregated)
      !
      ! eai        p_ext_atm%eai(nproma,nblks_c)
      cf_desc    = t_cf_var('eai', ' ','(evaporative) earth area index', datatype_flt)
      grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'eai', p_ext_atm%eai,         &
        &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,  &
        &           grib2_desc, ldims=shape2d_c, loutput=.TRUE.)

      ! Evaporative area index
      !
      ! eai_t       p_ext_atm%eai_t(nproma,nblks_c,ntiles_total)
      cf_desc    = t_cf_var('evaporative_surface_area_index_vegetation_period', '-',&
        &                   'Earth Area (evaporative surface area) Index', datatype_flt)
      grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'eai_t', p_ext_atm%eai_t,     &
        &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,  &
        &           grib2_desc, ldims=shape3d_nt, loutput=.FALSE. )


      ! root depth of vegetation
      !
      ! rootdp      p_ext_atm%rootdp(nproma,nblks_c)
      cf_desc    = t_cf_var('root_depth_of_vegetation', 'm',&
        &                   'root depth of vegetation', datatype_flt)
      grib2_desc = grib2_var( 2, 0, 32, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'rootdp', p_ext_atm%rootdp,     &
        &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,    &
        &           grib2_desc, ldims=shape2d_c, loutput=.TRUE.,    &
        &           isteptype=TSTEP_CONSTANT )

      ! rootdp_t      p_ext_atm%rootdp_t(nproma,nblks_c,ntiles_total)
      cf_desc    = t_cf_var('root_depth_of_vegetation', 'm',&
        &                   'root depth of vegetation', datatype_flt)
      grib2_desc = grib2_var( 2, 0, 32, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'rootdp_t', p_ext_atm%rootdp_t, &
        &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,    &
        &           grib2_desc, ldims=shape3d_nt, loutput=.FALSE. )

      ! evergreen forest
      !
      ! for_e        p_ext_atm%for_e(nproma,nblks_c)
      cf_desc    = t_cf_var('fraction_of_evergreen_forest_cover', '-',&
        &                   'Fraction of evergreen forest', datatype_flt)
      grib2_desc = grib2_var( 2, 0, 29, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'for_e', p_ext_atm%for_e,       &
        &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,    &
        &           grib2_desc, ldims=shape2d_c, loutput=.FALSE. )



      ! deciduous forest
      !
      ! for_d     p_ext_atm%for_d(nproma,nblks_c)
      cf_desc    = t_cf_var('fraction_of_deciduous_forest_cover', '-',&
        &                   'Fraction of deciduous forest', datatype_flt)
      grib2_desc = grib2_var( 2, 0, 30, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'for_d', p_ext_atm%for_d,       &
        &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,    &
        &           grib2_desc, ldims=shape2d_c )


      ! Minimal stomata resistence
      !
      ! rsmin        p_ext_atm%rsmin(nproma,nblks_c)
      cf_desc    = t_cf_var('RSMIN', 's m-1', 'Minimal stomata resistence', datatype_flt)
      grib2_desc = grib2_var( 2, 0, 16, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'rsmin', p_ext_atm%rsmin,       &
        &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,    &
        &           grib2_desc, ldims=shape2d_c, loutput=.TRUE.,    &
        &           isteptype=TSTEP_CONSTANT )

      ! rsmin2d_t        p_ext_atm%rsmin2d_t(nproma,nblks_c,ntiles_total)
      cf_desc    = t_cf_var('RSMIN', 's m-1', 'Minimal stomata resistence', datatype_flt)
      grib2_desc = grib2_var( 2, 0, 16, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'rsmin2d_t', p_ext_atm%rsmin2d_t,       &
        &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,            &
        &           grib2_desc, ldims=shape3d_nt, loutput=.FALSE. )

      ! NDVI yearly maximum
      !
      ! ndvi_max        p_ext_atm%ndvi_max(nproma,nblks_c)
      cf_desc    = t_cf_var('normalized_difference_vegetation_index', '-', &
        &                   'NDVI yearly maximum', datatype_flt)
      grib2_desc = grib2_var( 2, 0, 31, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'ndvi_max', p_ext_atm%ndvi_max, &
        &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,    &
        &           grib2_desc, ldims=shape2d_c, loutput=.TRUE.  )

      ! proportion of actual value/maximum NDVI (at ini_datetime)
      !
      ! ndviratio        p_ext_atm%ndviratio(nproma,nblks_c)
      cf_desc    = t_cf_var('normalized_difference_vegetation_index', '-',     &
        &                   '(monthly) proportion of actual value/maximum ' // &
        &                   'NDVI (at init time)', datatype_flt)
      grib2_desc = grib2_var( 2, 0, 192, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'ndviratio', p_ext_atm%ndviratio, &
        &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,      &
        &           grib2_desc, ldims=shape2d_c, loutput=.TRUE.  )

      ! Control fields for tile approach
      ! idx_lst_lp          p_ext_atm%idx_lst_lp(nproma,nblks_c)
      cf_desc    = t_cf_var('land point index list', '-', &
        &                   'land point index list', datatype_flt)
      grib2_desc = grib2_var( 255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'idx_lst_lp', p_ext_atm%idx_lst_lp, &
        &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,        &
        &           grib2_desc, ldims=shape2d_c, loutput=.FALSE. )

      ! idx_lst_sp          p_ext_atm%idx_lst_sp(nproma,nblks_c)
      cf_desc    = t_cf_var('sea point index list', '-', &
        &                   'sea point index list', datatype_flt)
      grib2_desc = grib2_var( 255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'idx_lst_sp', p_ext_atm%idx_lst_sp, &
        &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,        &
        &           grib2_desc, ldims=shape2d_c, loutput=.FALSE. )

      ! idx_lst_fp          p_ext_atm%idx_lst_sp(nproma,nblks_c)
      cf_desc    = t_cf_var('lake point index list', '-', &
        &                   'lake point index list', datatype_flt)
      grib2_desc = grib2_var( 255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'idx_lst_fp', p_ext_atm%idx_lst_fp, &
        &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,        &
        &           grib2_desc, ldims=shape2d_c, loutput=.FALSE. )

      ! idx_lst_lp_t        p_ext_atm%idx_lst_lp_t(nproma,nblks_c,ntiles_total)
      cf_desc    = t_cf_var('static land tile point index list', '-', &
        &                   'static land tile point index list', datatype_flt)
      grib2_desc = grib2_var( 255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'idx_lst_lp_t', p_ext_atm%idx_lst_lp_t, &
        &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,            &
        &           grib2_desc, ldims=shape3d_nt, loutput=.FALSE. )

      ! idx_lst_t        p_ext_atm%idx_lst_t(nproma,nblks_c,ntiles_total)
      cf_desc    = t_cf_var('dynamic land tile point index list', '-', &
        &                   'dynamic land tile point index list', datatype_flt)
      grib2_desc = grib2_var( 255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'idx_lst_t', p_ext_atm%idx_lst_t, &
        &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,      &
        &           grib2_desc, ldims=shape3d_nt, loutput=.FALSE. )


      ! idx_lst_spw      p_ext_atm%idx_lst_spw(nproma,nblks_c)
      cf_desc    = t_cf_var('sea water point index list', '-', &
        &                   'sea water point index list', datatype_flt)
      grib2_desc = grib2_var( 255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'idx_lst_spw', p_ext_atm%idx_lst_spw, &
        &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,          &
        &           grib2_desc, ldims=shape2d_c, loutput=.FALSE. )

      ! idx_lst_spi      p_ext_atm%idx_lst_spi(nproma,nblks_c)
      cf_desc    = t_cf_var('sea ice point index list', '-', &
        &                   'sea ice point index list', datatype_flt)
      grib2_desc = grib2_var( 255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'idx_lst_spi', p_ext_atm%idx_lst_spi, &
        &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,          &
        &           grib2_desc, ldims=shape2d_c, loutput=.FALSE. )


      ! snowtile_flag_t   p_ext_atm%snowtile_flag_t(nproma,nblks_c,ntiles_total)
      ! -1: no separation between snow tile and snow-free tile
      !  0: inactive
      !  1: active
      !  2: newly activated; initialization from corresponding tile required
      cf_desc    = t_cf_var('flag of activity', '-', &
        &                   'flag of activity', datatype_flt)
      grib2_desc = grib2_var( 255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'snowtile_flag_t', p_ext_atm%snowtile_flag_t, &
        &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,                  &
        &           grib2_desc, ldims=shape3d_nt, loutput=.FALSE. )

      ! not sure if these dimensions are supported by add_var...
      ALLOCATE(p_ext_atm%lp_count(nblks_c), p_ext_atm%gp_count_t(nblks_c,ntiles_total), &
               p_ext_atm%lp_count_t(nblks_c,ntiles_total) )
      ALLOCATE(p_ext_atm%sp_count (nblks_c),p_ext_atm%fp_count (nblks_c))

      ! allocate grid point counts per block for dynamic ocean ice/water point
      ! index lists
      ALLOCATE(p_ext_atm%spw_count(nblks_c),p_ext_atm%spi_count(nblks_c))



      ! lc_class_t        p_ext_atm%lc_class_t(nproma,nblks_c,ntiles_total+ntiles_water)
      cf_desc    = t_cf_var('tile point land cover class', '-', &
        &                   'tile point land cover class', datatype_flt)
      grib2_desc = grib2_var( 2, 0, 35, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'lc_class_t', p_ext_atm%lc_class_t, &
        &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,        &
        &           grib2_desc, ldims=shape3d_ntw,                      &
        &           loutput=.FALSE., lcontainer=.TRUE. )

      ! fill the separate variables belonging to the container lc_class_t
      ALLOCATE(p_ext_atm%lc_class_t_ptr(ntiles_total+ntiles_water))
      DO jsfc = 1,ntiles_total + ntiles_water
      WRITE(csfc,'(i2)') jsfc
      CALL add_ref( p_ext_atm_list, 'lc_class_t', 'lc_class_t_'//TRIM(ADJUSTL(csfc)),  &
        &           p_ext_atm%lc_class_t_ptr(jsfc)%p_2d,                               &
        &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                &
        &           t_cf_var('lc_class_t_'//csfc, '-', '', datatype_flt),            &
        &           grib2_desc,                                                        &
        &           hor_interp=create_hor_interp_metadata(hor_intp_type=HINTP_TYPE_LONLAT_NNB),&
        &           var_class=CLASS_TILE,                                              &
        &           ldims=shape2d_c, loutput=.TRUE.,                                   &
        &           post_op=post_op(POST_OP_LUC, new_cf=cf_desc, arg1=i_lctype(jg)) )
      ENDDO



      ! lc_frac_t        p_ext_atm%lc_frac_t(nproma,nblks_c,ntiles_total+ntiles_water)
      cf_desc    = t_cf_var('lc_frac_t', '-', &
        &                   'tile point land cover fraction list', datatype_flt)
      grib2_desc = grib2_var( 255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'lc_frac_t', p_ext_atm%lc_frac_t, &
        &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,      &
        &           grib2_desc, ldims=shape3d_ntw, loutput=.FALSE. )


      ! frac_t        p_ext_atm%frac_t(nproma,nblks_c,ntiles_total+ntiles_water)
      cf_desc    = t_cf_var('frac_t', '-', &
        &                   'tile point area fraction list', datatype_flt)
      grib2_desc = grib2_var( 2, 0, 36, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'frac_t', p_ext_atm%frac_t,   &
        &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,  &
        &           grib2_desc, ldims=shape3d_ntw, loutput=.FALSE., lcontainer=.TRUE.)

      ! fill the separate variables belonging to the container frac_t
      ALLOCATE(p_ext_atm%frac_t_ptr(ntiles_total+ntiles_water))
      DO jsfc = 1,ntiles_total + ntiles_water
      WRITE(csfc,'(i2)') jsfc
      CALL add_ref( p_ext_atm_list, 'frac_t', 'frac_t_'//TRIM(ADJUSTL(csfc)),  &
        &           p_ext_atm%frac_t_ptr(jsfc)%p_2d,                           &
        &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                        &
        &           t_cf_var('frac_t_'//csfc, '-', '', datatype_flt),        &
        &           grib2_desc,                                                &
        &           var_class=CLASS_TILE,                                      &
        &           ldims=shape2d_c, loutput=.TRUE. )
      ENDDO


      ! inv_frland_from_tiles      p_ext_atm%inv_frland_from_tiles(nproma,nblks_c)
      cf_desc    = t_cf_var('inv_frland_from_tiles', '-', &
        &                   'inverse of fr_land derived from land tiles', datatype_flt)
      grib2_desc = grib2_var( 255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'inv_frland_from_tiles', p_ext_atm%inv_frland_from_tiles,&
        &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,                          &
        &           grib2_desc, ldims=shape2d_c, loutput=.FALSE.)


      ! Storage for table values - not sure if these dimensions are supported by add_var
      ALLOCATE(p_ext_atm%z0_lcc(nclass_lu(jg)),         & ! Land-cover related roughness length
                p_ext_atm%z0_lcc_min(nclass_lu(jg)),    & ! Minimum land-cover related roughness length
                p_ext_atm%plcovmax_lcc(nclass_lu(jg)),  & ! Maximum plant cover fraction for each land-cover class
                p_ext_atm%laimax_lcc(nclass_lu(jg)),    & ! Maximum leaf area index for each land-cover class
                p_ext_atm%rootdmax_lcc(nclass_lu(jg)),  & ! Maximum root depth each land-cover class
                p_ext_atm%stomresmin_lcc(nclass_lu(jg)),& ! Minimum stomata resistance for each land-cover class
                p_ext_atm%snowalb_lcc(nclass_lu(jg)),   & ! Albedo in case of snow cover for each land-cover class
                p_ext_atm%snowtile_lcc(nclass_lu(jg))   ) ! Specification of snow tiles for land-cover class


      !--------------------------------
      ! soil parameters
      !--------------------------------

      ! soil type
      !
      ! soiltyp      p_ext_atm%soiltyp(nproma,nblks_c)
      cf_desc    = t_cf_var('soil_type', '-','soil type', datatype_flt)
      grib2_desc = grib2_var( 2, 3, 196, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'soiltyp', p_ext_atm%soiltyp,   &
        &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,    &
        &           grib2_desc, ldims=shape2d_c, loutput=.TRUE.,    &
        &           hor_interp=create_hor_interp_metadata(          &
        &               hor_intp_type=HINTP_TYPE_LONLAT_NNB ),      &
        &           isteptype=TSTEP_CONSTANT )

      ! soiltyp_t      p_ext_atm%soiltyp_t(nproma,nblks_c,ntiles_total)
      cf_desc    = t_cf_var('soil_type', '-','soil type', datatype_flt)
      grib2_desc = grib2_var( 2, 3, 196, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'soiltyp_t', p_ext_atm%soiltyp_t,   &
        &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,        &
        &           grib2_desc, ldims=shape3d_nt, loutput=.FALSE. )


      ! Climat. temperature
      ! Climat. temperature 2m above ground. However, this temperature is used
      ! to initialize the climatological layer of the soil model (lowermost layer)
      !
      ! t_cl         p_ext_atm%t_cl(nproma,nblks_c)
      cf_desc    = t_cf_var('soil_temperature', 'K',                  &
        &                   'CRU near surface temperature climatology', datatype_flt)
      grib2_desc = grib2_var( 0, 0, 0, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_ext_atm_list, 't_cl', p_ext_atm%t_cl,           &
        &           GRID_UNSTRUCTURED_CELL, ZA_HEIGHT_2M, cf_desc,    &
        &           grib2_desc, ldims=shape2d_c, loutput=.TRUE.,      &
        &           isteptype=TSTEP_CONSTANT )

      IF (itype_vegetation_cycle > 1) THEN
        ! t2m_clim         p_ext_atm%t2m_clim(nproma,nblks_c)
        cf_desc    = t_cf_var('2m_temperature', 'K',                  &
          &                   'T2M interpolated from monthly climatology', datatype_flt)
        grib2_desc = grib2_var( 0, 0, 0, ibits, GRID_UNSTRUCTURED, GRID_CELL)
        CALL add_var( p_ext_atm_list, 't2m_clim', p_ext_atm%t2m_clim,   &
          &           GRID_UNSTRUCTURED_CELL, ZA_HEIGHT_2M, cf_desc,    &
          &           grib2_desc, ldims=shape2d_c, loutput=.TRUE.,      &
          &           isteptype=TSTEP_CONSTANT )

        ! t2m_clim_hc         p_ext_atm%t2m_clim_hc(nproma,nblks_c)
        cf_desc    = t_cf_var('Height-corrected 2m_temperature', 'K',                  &
          &                   'Height-corrected T2M interpolated from monthly climatology', datatype_flt)
        grib2_desc = grib2_var( 0, 0, 0, ibits, GRID_UNSTRUCTURED, GRID_CELL)
        CALL add_var( p_ext_atm_list, 't2m_clim_hc', p_ext_atm%t2m_clim_hc,   &
          &           GRID_UNSTRUCTURED_CELL, ZA_HEIGHT_2M, cf_desc,    &
          &           grib2_desc, ldims=shape2d_c, loutput=.TRUE.,      &
          &           isteptype=TSTEP_CONSTANT )

        ! t2m_climgrad         p_ext_atm%t2m_climgrad(nproma,nblks_c)
        cf_desc    = t_cf_var('2m_temperature_gradient', 'K/month',      &
          &                   'climatology T2M gradient', datatype_flt)
        grib2_desc = grib2_var( 0, 0, 0, ibits, GRID_UNSTRUCTURED, GRID_CELL)
        CALL add_var( p_ext_atm_list, 't2m_climgrad', p_ext_atm%t2m_climgrad, &
          &           GRID_UNSTRUCTURED_CELL, ZA_HEIGHT_2M, cf_desc,    &
          &           grib2_desc, ldims=shape2d_c, loutput=.TRUE.,      &
          &           isteptype=TSTEP_CONSTANT )

      ENDIF

      ! longwave surface emissivity
      !
      ! emis_rad     p_ext_atm%emis_rad(nproma,nblks_c)
      cf_desc    = t_cf_var('emis_rad', '-', 'longwave surface emissivity', datatype_flt)
      grib2_desc = grib2_var( 2, 3, 199, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'emis_rad', p_ext_atm%emis_rad, &
        &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,    &
        &           grib2_desc, ldims=shape2d_c, loutput=.TRUE.,    &
        &           isteptype=TSTEP_CONSTANT )


      ! landuse class fraction
      !
      ! lu_class_fraction    p_ext_atm%lu_class_fraction(nproma,nblks_c,nclass_lu)
      cf_desc    = t_cf_var('lu_class_fraction', '-', 'landuse class fraction', datatype_flt)
      grib2_desc = grib2_var( 2, 0, 36, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'lu_class_fraction', p_ext_atm%lu_class_fraction, &
        &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,    &
        &           grib2_desc, ldims=shape3d_sfc, loutput=.FALSE. )


      !--------------------------------
      ! If MODIS albedo is used
      !--------------------------------
      IF ( albedo_type == MODIS) THEN

        ! Shortwave broadband albedo for diffuse radiation (0.3 - 5.0 micron), snow-free
        !
        ! alb_dif    p_ext_atm%alb_dif(nproma,nblks_c,ntimes)
        cf_desc    = t_cf_var('Shortwave_albedo_diffuse', '-', &
          &                   'Shortwave albedo for diffuse radiation', datatype_flt)
        grib2_desc = grib2_var(0, 19, 18, ibits, GRID_UNSTRUCTURED, GRID_CELL)
        CALL add_var( p_ext_atm_list, 'alb_dif', p_ext_atm%alb_dif,               &
          &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,    &
          &           ldims=shape2d_c, loutput=.TRUE.                            )

        ! UV visible broadband albedo for diffuse radiation (0.3 - 0.7 micron)
        !
        ! albuv_dif    p_ext_atm%albuv_dif(nproma,nblks_c,ntimes)
        cf_desc    = t_cf_var('UV_visible_albedo_diffuse', '-', &
          &                   'UV visible albedo for diffuse radiation', datatype_flt)
        grib2_desc = grib2_var(0, 19, 222, ibits, GRID_UNSTRUCTURED, GRID_CELL)
        CALL add_var( p_ext_atm_list, 'albuv_dif', p_ext_atm%albuv_dif,           &
          &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,    &
          &           ldims=shape2d_c, loutput=.TRUE.                             )

        ! Near IR broadband albedo for diffuse radiation (0.7 - 5.0 micron)
        !
        ! albni_dif    p_ext_atm%albni_dif(nproma,nblks_c,ntimes)
        cf_desc    = t_cf_var('Near_IR_albedo_diffuse', '-', &
          &                   'Near IR albedo for diffuse radiation', datatype_flt)
        grib2_desc = grib2_var(0, 19, 223, ibits, GRID_UNSTRUCTURED, GRID_CELL)
        CALL add_var( p_ext_atm_list, 'albni_dif', p_ext_atm%albni_dif,           &
          &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,    &
          &           ldims=shape2d_c, loutput=.TRUE.                             )

      END IF  ! albedo_type

    END IF ! iforcing


    IF ( iforcing == iecham .OR. iforcing == ildf_echam ) THEN

      ! atmosphere land-sea-mask at surface on cell centers
      ! lsm_ctr_c  p_ext_atm%lsm_ctr_c(nproma,nblks_c)
      cf_desc    = t_cf_var('Atmosphere model land-sea-mask at cell center', '-2/-1/1/2', &
        &                   'Atmosphere model land-sea-mask', datatype_flt)
      grib2_desc = grib2_var( 192, 140, 219, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'lsm_ctr_c', p_ext_atm%lsm_ctr_c,        &
        &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,             &
        grib2_desc, ldims=shape2d_c )

      IF (iequations == ihs_atm_temp .OR. iequations == ihs_atm_theta ) THEN
        ! elevation p_ext_atm%elevation_c(nproma,nblks_c)
        cf_desc    = t_cf_var('elevation at cell center', 'm', &
          &                     'elevation', datatype_flt)
        grib2_desc = grib2_var( 192, 140, 219, ibits, GRID_UNSTRUCTURED, GRID_CELL)
        CALL add_var( p_ext_atm_list, 'elevation_c', p_ext_atm%elevation_c,        &
          &             GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,          &
          grib2_desc, ldims=shape2d_c )
      END IF

      ! longwave surface emissivity
      !
      ! emis_rad     p_ext_atm%emis_rad(nproma,nblks_c)
      cf_desc    = t_cf_var('emis_rad', '-', 'longwave surface emissivity', datatype_flt)
      grib2_desc = grib2_var( 2, 3, 199, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'emis_rad', p_ext_atm%emis_rad, &
        &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,    &
        &           grib2_desc, ldims=shape2d_c, loutput=.FALSE. )

      ! HDmodel land-sea-mask at surface on cell centers
      ! lsm_hd_c   p_ext_atm%lsm_hd_c(nproma,nblks_c)
      cf_desc    = t_cf_var('HD model land-sea-mask at cell center', '-2/-1/1/2', &
        &                   'HD model land-sea-mask', datatype_flt)
      grib2_desc = grib2_var( 192, 140, 219, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'lsm_hd_c', p_ext_atm%lsm_hd_c,          &
        &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,             &
        grib2_desc, ldims=shape2d_c )

    END IF

  END SUBROUTINE new_ext_data_atm_list

  !-------------------------------------------------------------------------
  !
  !
  !>
  !! Allocation of atmospheric external data structure (time dependent)
  !!
  !! Allocation of atmospheric external data structure (time dependent
  !! elements).
  !!
  !! Initialization of elements with zero.
  !!
  !! @par Revision History
  !! Initial release by Daniel Reinert (2011-05-03)
  !!
  SUBROUTINE new_ext_data_atm_td_list ( p_patch, p_ext_atm_td, &
    &                               p_ext_atm_td_list, listname)
!
    TYPE(t_patch), TARGET    , INTENT(IN)   :: & !< current patch
      &  p_patch

    TYPE(t_external_atmos_td), INTENT(INOUT):: & !< current external data structure
      &  p_ext_atm_td

    TYPE(t_var_list)         , INTENT(INOUT):: & !< current external data list
      &  p_ext_atm_td_list

    CHARACTER(len=*)         , INTENT(IN)   :: & !< list name
      &  listname

    TYPE(t_cf_var)    :: cf_desc
    TYPE(t_grib2_var) :: grib2_desc

    INTEGER :: nblks_c      !< number of cell blocks to allocate
    INTEGER :: jg           !< patch ID

    INTEGER :: shape3d_c(3)
    INTEGER :: shape4d_c(4)
    INTEGER :: shape3d_sstice(3)

    INTEGER :: ibits         !< "entropy" of horizontal slice
    INTEGER :: datatype_flt  !< floating point accuracy in NetCDF output

    CHARACTER(len=max_char_length), PARAMETER :: &
      routine = modname//':new_ext_data_atm_td_list'
    !--------------------------------------------------------------

    !determine size of arrays
    nblks_c = p_patch%nblks_c

    ! get patch ID
    jg = p_patch%id

    ibits  = 16   ! "entropy" of horizontal slice

    IF ( lnetcdf_flt64_output ) THEN
      datatype_flt = DATATYPE_FLT64
    ELSE
      datatype_flt = DATATYPE_FLT32
    ENDIF

    ! predefined array shapes
    shape3d_c   = (/ nproma, nblks_c, nmonths_ext(jg)  /)
    shape4d_c   = (/ nproma, nlev_o3, nblks_c, nmonths /)

    IF (iforcing == inwp) THEN
      SELECT CASE (sstice_mode)
        CASE(SSTICE_ANA)
          ! nothing to do
        CASE(SSTICE_ANA_CLINC, SSTICE_CLIM)
          shape3d_sstice = (/ nproma, nblks_c, 12 /)
        CASE(SSTICE_AVG_MONTHLY)
          shape3d_sstice = (/ nproma, nblks_c,  2 /)
        CASE(SSTICE_AVG_DAILY)
          CALL finish (TRIM(routine), 'sstice_mode=4  not implemented!')
        CASE DEFAULT
          CALL finish (TRIM(routine), 'sstice_mode not valid!')
      END SELECT
    END IF

    !
    ! Register a field list and apply default settings
    !
    CALL new_var_list( p_ext_atm_td_list, TRIM(listname), patch_id=jg )
    CALL default_var_list_settings( p_ext_atm_td_list,         &
                                  & lrestart=.FALSE.,          &
                                  & loutput=.TRUE.  )


    !--------------------------------
    ! radiation parameters
    !--------------------------------


    IF (iforcing == inwp) THEN

    ! ozone on pressure levels
    ! ATTENTION: a GRIB2 number will go to
    ! the ozone mass mixing ratio...
    !
    IF ( irad_o3 == io3_clim .OR. irad_o3 == io3_ape ) THEN

      CALL message(routine, 'generate ext ozone field')

      ! o3  main height level from read-in file
      cf_desc    = t_cf_var('O3_zf', 'm',   &
        &                   'ozone geometric height level', datatype_flt)
      grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_ext_atm_td_list, 'O3_zf', p_ext_atm_td%zf,  &
        &           GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc,   &
        &           grib2_desc, ldims=(/nlev_o3/), loutput=.FALSE.  )

      ! o3  main pressure level from read-in file
      cf_desc    = t_cf_var('O3_pf', 'Pa',   &
        &                   'ozone main pressure level', datatype_flt)
      grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_ext_atm_td_list, 'O3_pf', p_ext_atm_td%pfoz, &
        &           GRID_UNSTRUCTURED_CELL, ZA_PRESSURE, cf_desc,  &
        &           grib2_desc, ldims=(/nlev_o3/), loutput=.FALSE.  )

      ! o3  intermediate pressure level
      cf_desc    = t_cf_var('O3_ph', 'Pa',   &
        &                   'ozone intermediate pressure level', datatype_flt)
      grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_ext_atm_td_list, 'O3_ph', p_ext_atm_td%phoz, &
        &           GRID_UNSTRUCTURED_CELL, ZA_PRESSURE, cf_desc,  &
        &           grib2_desc, ldims=(/nlev_o3+1/), loutput=.FALSE.  )

      ! o3       p_ext_atm_td%o3(nproma,nlev_o3,nblks_c,nmonths)
      cf_desc    = t_cf_var('O3', TRIM(o3unit),   &
        &                   'mole_fraction_of_ozone_in_air', datatype_flt)
      grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_ext_atm_td_list, 'O3', p_ext_atm_td%O3, &
        &           GRID_UNSTRUCTURED_CELL, ZA_PRESSURE, cf_desc, &
        &           grib2_desc, ldims=shape4d_c, loutput=.FALSE.  )

    END IF ! irad_o3

    ! Black carbon aerosol
    !
    ! aer_bc       p_ext_atm_td%aer_bc(nproma,nblks_c,ntimes)
    cf_desc    = t_cf_var('aerosol optical thickness of black carbon', '-',   &
      &                   'atmosphere_absorption_optical_thickness_due_to_' //&
      &                   'black_carbon_ambient_aerosol', datatype_flt)
    grib2_desc = grib2_var( 0, 20, 102, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_ext_atm_td_list, 'aer_bc', p_ext_atm_td%aer_bc, &
      &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,      &
      &           grib2_desc, ldims=shape3d_c, loutput=.FALSE.,     &
      &           isteptype=TSTEP_AVG )  ! Meta info constituentType missing


    ! Dust aerosol
    !
    ! aer_dust     p_ext_atm_td%aer_dust(nproma,nblks_c,ntimes)
    cf_desc    = t_cf_var('aot_dust', '-', &
      &                   'atmosphere absorption optical thickness due '//  &
      &                   'to dust ambient aerosol', datatype_flt)
    grib2_desc = grib2_var( 0, 20, 102, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_ext_atm_td_list, 'aer_dust', p_ext_atm_td%aer_dust, &
      &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
      &           ldims=shape3d_c, loutput=.FALSE.,                        &
      &           isteptype=TSTEP_AVG )  ! Meta info constituentType missing


    ! Organic aerosol
    !
    ! aer_org      p_ext_atm_td%aer_org(nproma,nblks_c,ntimes)
    cf_desc    = t_cf_var('aot_org', '-', &
      &                   'atmosphere absorption optical thickness due '//  &
      &                   'to particulate organic matter ambient aerosol', datatype_flt)
    grib2_desc = grib2_var( 0, 20, 102, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_ext_atm_td_list, 'aer_org', p_ext_atm_td%aer_org,     &
      &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,&
      &           ldims=shape3d_c, loutput=.FALSE.,                       &
      &           isteptype=TSTEP_AVG )  ! Meta info constituentType missing


    ! Sulfate aerosol
    !
    ! aer_so4      p_ext_atm_td%aer_so4(nproma,nblks_c,ntimes)
    cf_desc    = t_cf_var('aot_so4', '-', &
      &                   'atmosphere absorption optical thickness due '//  &
      &                   'to sulfate_ambient_aerosol', datatype_flt)
    grib2_desc = grib2_var( 0, 20, 102, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_ext_atm_td_list, 'aer_so4', p_ext_atm_td%aer_so4, &
      &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,&
      &           ldims=shape3d_c, loutput=.FALSE.,                       &
      &           isteptype=TSTEP_AVG )  ! Meta info constituentType missing


    ! Seasalt aerosol
    !
    ! aer_ss       p_ext_atm_td%aer_ss(nproma,nblks_c,ntimes)
    cf_desc    = t_cf_var('aot_ss', '-', &
      &                   'atmosphere absorption optical thickness due '//  &
      &                   'to seasalt_ambient_aerosol', datatype_flt)
    grib2_desc = grib2_var( 0, 20, 102, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_ext_atm_td_list, 'aer_ss', p_ext_atm_td%aer_ss, &
      &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,&
      &           ldims=shape3d_c, loutput=.FALSE.,                       &
      &           isteptype=TSTEP_AVG )  ! Meta info constituentType missing


    !--------------------------------
    ! vegetation parameters
    !--------------------------------

    ! (monthly) proportion of actual value/maximum NDVI
    !
    ! ndvi_mrat     p_ext_atm_td%ndvi_mrat(nproma,nblks_c,ntimes)
    cf_desc    = t_cf_var('normalized_difference_vegetation_index', '-', &
      &                   '(monthly) proportion of actual value/maximum ' // &
      &                   'normalized differential vegetation index', datatype_flt)
    grib2_desc = grib2_var( 2, 0, 192, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_ext_atm_td_list, 'ndvi_mrat', p_ext_atm_td%ndvi_mrat,  &
      &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
      &           ldims=shape3d_c, loutput=.FALSE.,                         &
      &           isteptype=TSTEP_AVG )



    !--------------------------------
    ! If MODIS albedo is used
    !--------------------------------
    IF ( albedo_type == MODIS) THEN

      ! (monthly)  Shortwave broadband albedo for diffuse radiation (0.3 - 5.0 micron), snow-free
      !
      ! alb_dif    p_ext_atm_td%alb_dif(nproma,nblks_c,ntimes)
      cf_desc    = t_cf_var('Shortwave_albedo_diffuse', '-', &
        &                   'Shortwave albedo for diffuse radiation', datatype_flt)
      grib2_desc = grib2_var(0, 19, 18, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_ext_atm_td_list, 'alb_dif', p_ext_atm_td%alb_dif,         &
        &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,    &
        &           ldims=shape3d_c, loutput=.FALSE.,                           &
        &           isteptype=TSTEP_AVG )

      ! (monthly)  UV visible broadband albedo for diffuse radiation (0.3 - 0.7 micron)
      !
      ! albuv_dif    p_ext_atm_td%albuv_dif(nproma,nblks_c,ntimes)
      cf_desc    = t_cf_var('UV_visible_albedo_diffuse', '-', &
        &                   'UV visible albedo for diffuse radiation', datatype_flt)
      grib2_desc = grib2_var(0, 19, 222, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_ext_atm_td_list, 'albuv_dif', p_ext_atm_td%albuv_dif,     &
        &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,    &
        &           ldims=shape3d_c, loutput=.FALSE.,                           &
        &           isteptype=TSTEP_AVG )

      ! (monthly)  Near IR broadband albedo for diffuse radiation (0.7 - 5.0 micron)
      !
      ! albni_dif    p_ext_atm_td%albni_dif(nproma,nblks_c,ntimes)
      cf_desc    = t_cf_var('Near_IR_albedo_diffuse', '-', &
        &                   'Near IR albedo for diffuse radiation', datatype_flt)
      grib2_desc = grib2_var(0, 19, 223, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_ext_atm_td_list, 'albni_dif', p_ext_atm_td%albni_dif,     &
        &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,    &
        &           ldims=shape3d_c, loutput=.FALSE.,                           &
        &           isteptype=TSTEP_AVG )

    ENDIF  ! albedo_type


    IF (itype_vegetation_cycle > 1) THEN
      ! t2m_m     p_ext_atm_td%t2m_m(nproma,nblks_c,ntimes)
      cf_desc    = t_cf_var('t2m_m', 'K', &
        &                   '(monthly) 2-metre temperature ', datatype_flt)
      grib2_desc = grib2_var(0, 0, 0, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_ext_atm_td_list, 't2m_m', p_ext_atm_td%t2m_m, &
        &           GRID_UNSTRUCTURED_CELL, ZA_HEIGHT_2M, cf_desc, grib2_desc,&
        &           ldims=shape3d_c, loutput=.FALSE. )
    ENDIF

    !--------------------------------
    !SST and sea ice fraction
    !--------------------------------
    SELECT CASE (sstice_mode)
      CASE (SSTICE_ANA_CLINC)  ! SST is read from analysis and is updated by climatological increments 
                               ! on a daily basis. Therefore, sst_m is required to store the monthly fields
        ! sst_m     p_ext_atm_td%sst_m(nproma,nblks_c,ntimes)
        cf_desc    = t_cf_var('sst_m', 'K', &
          &                   '(monthly) sea surface temperature '  &
          &                   , datatype_flt)
        grib2_desc = grib2_var(10 ,3 ,0, ibits, GRID_UNSTRUCTURED, GRID_CELL)
        CALL add_var( p_ext_atm_td_list, 'sst_m', p_ext_atm_td%sst_m, &
          &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,&
          &           ldims=shape3d_sstice, loutput=.FALSE. )
        !
      CASE (SSTICE_CLIM,SSTICE_AVG_MONTHLY,SSTICE_AVG_DAILY)
        !
        ! sst_m     p_ext_atm_td%sst_m(nproma,nblks_c,ntimes)
        cf_desc    = t_cf_var('sst_m', 'K', &
          &                   '(monthly) sea surface temperature '  &
          &                   , datatype_flt)
        grib2_desc = grib2_var(10 ,3 ,0, ibits, GRID_UNSTRUCTURED, GRID_CELL)
        CALL add_var( p_ext_atm_td_list, 'sst_m', p_ext_atm_td%sst_m, &
          &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,&
          &           ldims=shape3d_sstice, loutput=.FALSE. )
        !
        ! fr_ice_m     p_ext_atm_td%fr_ice_m(nproma,nblks_c,ntimes)
        cf_desc    = t_cf_var('fr_ice_m', '(0-1)', &
          &                   '(monthly) sea ice fraction '  &
          &                   , datatype_flt)
        grib2_desc = grib2_var( 192,128 ,31 , ibits, GRID_UNSTRUCTURED, GRID_CELL)
        CALL add_var( p_ext_atm_td_list, 'fr_ice_m', p_ext_atm_td%fr_ice_m, &
          &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,&
          &           ldims=shape3d_sstice, loutput=.FALSE. )
        !
      CASE default
        ! do nothing
        !
    END SELECT

    ENDIF ! inwp

  END SUBROUTINE new_ext_data_atm_td_list
  !-------------------------------------------------------------------------



  !-------------------------------------------------------------------------
  !>
  !! Destruct external data data structure and lists
  !!
  !! Destruct external data data structure and lists
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2011-05-04)
  !!
  SUBROUTINE destruct_ext_data

    INTEGER :: jg
    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
      routine = modname//':destruct_ext_data'
    !-------------------------------------------------------------------------

    CALL message (TRIM(routine), 'Destruction of data structure for ' // &
      &                          'external data started')

    DO jg = 1,n_dom
      ! Delete list of constant in time atmospheric elements
      CALL delete_var_list( ext_data(jg)%atm_list )
    ENDDO

    IF (iforcing > 1 ) THEN
    DO jg = 1,n_dom
      ! Delete list of time-dependent atmospheric elements
      CALL delete_var_list( ext_data(jg)%atm_td_list )
    ENDDO
    END IF

    CALL message (TRIM(routine), 'Destruction of data structure for ' // &
      &                          'external data finished')

  END SUBROUTINE destruct_ext_data
  !-------------------------------------------------------------------------

END MODULE mo_ext_data_state
