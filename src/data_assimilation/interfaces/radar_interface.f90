!+ Source module for the radar forward operator of the COSMO-model
!------------------------------------------------------------------------------
!NEC$ options "-finline-max-depth=3 -finline-max-function-size=1000"
MODULE radar_interface

  !------------------------------------------------------------------------------
  !
  ! Description:
  !   This module provides some routines which are necessary for the coupling
  !   of the radar forward operator to the ICON model.
  !
  ! Method:
  !   See subroutines below
  !
  ! Current Code Owner: DWD, Ulrich Blahak
  !  phone:  +49  69  8062 2393
  !  fax:    +49  69  8062 3721
  !  email:  ulrich.blahak@dwd.de
  !
  ! History:
  ! Version    Date       Name
  ! ---------- ---------- ----
  ! 1.1        2013/12/11 Ulrich Blahak
  !  Initial release
  ! 1.2        2017/06/06 Alberto de Lozar
  !  Warm bubbles added
  !
  ! Code Description:
  ! Language: Fortran 90.
  ! Software Standards: "European Standards for Writing and
  ! Documenting Exchangeable Fortran 90 Code".
  !==============================================================================
  !

  !------------------------------------------------------------------------------
  !
  ! Declarations:
  !
  ! Modules used:

  !==============================================================================
  !
  ! ICON modules:
  
  USE mo_kind,                  ONLY: sp, dp, wp
  USE mo_mpi,                   ONLY: my_process_is_mpi_workroot
  USE mo_exception,             ONLY: finish
  USE mo_grid_config,           ONLY: start_time, end_time ! in seconds since experiment start for each domain
  USE mo_parallel_config,       ONLY: nproma, blk_no, idx_no, idx_1d
  USE mo_nonhydro_state,        ONLY: p_nh_state
  USE mo_nwp_phy_state,         ONLY: prm_diag
  USE mo_model_domain,          ONLY: p_patch, t_patch
  USE mo_atm_phy_nwp_config,   ONLY: atm_phy_nwp_config
  USE mo_impl_constants,        ONLY: min_rlcell, min_rlcell_int, max_dom
  USE mo_impl_constants_grf,    ONLY: grf_bdywidth_c
  USE mo_loopindices,           ONLY: get_indices_c
  USE mo_time_config,           ONLY: time_config
  USE mo_util_mtime,            ONLY: getElapsedSimTimeInSeconds
  USE mo_run_config,            ONLY: dtime, nsteps, ntracer, ltimer, &
       &                              iqc,  iqi,  iqr,  iqs,  iqg,  iqh,  iqv, &
       &                              iqnc, iqni, iqnr, iqns, iqng, iqnh
!  USE mo_run_config,            ONLY: iqgl, qghl   ! later, if inwp_gscp=7 will be available
  USE mo_gribout_config,        ONLY: gribout_config
  USE mtime,                    ONLY: datetimeToString
  USE mo_vertical_coord_table,  ONLY: vct_a
  USE mo_physical_constants,    ONLY: rv           , & !> gas constant for water vapour
       &                              earth_radius , & !> [m]    average radius
       &                              rhoh2o       , & !> density of liquid water (kg/m^3)
       &                              rhoice       , & !> density of pure ice (kg/m^3)
       &                              K_w_0, K_i_0, Tmelt
  USE mo_intp_rbf,              ONLY: rbf_vec_interpol_cell
  USE mo_nh_diagnose_pres_temp, ONLY: diagnose_pres_temp
  USE mo_intp_data_strc,        ONLY: p_int_state
  USE mo_communication,         ONLY: exchange_data
  USE mo_2mom_mcrph_main,       ONLY: init_2mom_scheme,      &
       &                              rain_coeffs_2mom => rain_coeffs  ! contains the parameters for the mue-Dm-relation
  USE mo_2mom_mcrph_types,      ONLY: particle, particle_frozen, particle_lwf

  USE mo_timer,                 ONLY: timer_start, timer_stop, &
       &                              timer_radar_ini      , &
       &                              timer_radar_prep_compute, &
       &                              timer_radar_composites,&
       &                              timer_radar_bubbles  , &
       &                              timer_radar_compgrid , &
       &                              timer_radar_comm     , &
       &                              timer_radar_ongeom   , &
       &                              timer_radar_comppolar, &
       &                              timer_radar_out      , & 
       &                              timer_radar_barrier

  USE mo_opt_nwp_diagnostics,   ONLY: compute_field_dbz_1mom, compute_field_dbz_2mom
  USE gscp_data,                ONLY: cloud_num

!!$ There are other parameters available for the 1mom-scheme, but these are
!!$  not yet coupled explicitly to the EMVORADO 1mom reflectivity routines (at the moment
!!$  hardcoded)
!  USE gscp_data,                ONLY: cloud_num, isnow_n0temp, zami, mu_rain, zams, zams_gr, zbms, &
!    &                                 znimax_Thom, ztmix, zthn, mma, mmb, zcnue
  
  !==============================================================================
  !
  ! EMVORADO modules:

#ifdef HAVE_RADARFWO
  
  USE radar_data, ONLY : &
       time_mod_sec,    & ! model forecast time in [s] since model start
       ie_fwo, je_fwo, ke_fwo, &
       miss_threshold, miss_value, missval_int, &
       cmaxlen,         &
       num_compute_fwo,     & ! number of compute PEs
       num_radar,       & ! total number of radar PEs (num_compute + num_radario)
       num_radario,     & ! total number of radar-IO PEs
       radario_master_dom,& ! root-PEs of the radario group for each active radar domain (in the radar_dom comm., not radario_dom-comm.!!!)
       my_cart_id_fwo,  & ! rank of this PE (=subdomain) in the cartesian communicator
       my_radar_id,     & ! rank of this PE in the radar communicator (cart+radario)
       icomm_cart_fwo,  & ! communicator for the virtual cartesian topology
       icomm_radar,     & ! communicator for the group of radar-IO PEs (all domains) + compute PEs
       icomm_radar_dom, & ! communicator for the group of radar-IO PEs of each domain + compute PEs
       lcompute_pe_fwo, & ! indicates whether this is a compute PE or not
       radar_meta_type, & ! TYPE to hold the radar meta informations
       radar_grid_type, & ! TYPE to hold the pointers to the data on the aux azimutal slice grid
       i_fwo_prep_compute,& ! Timing flag
       i_fwo_bubbles   ,& ! Timing flag
       i_fwo_composites,& ! Timing flag
       i_fwo_ini,       &  ! Timing flag for the initialization of the forward operator
       i_fwo_compgrid,  &  ! Timing flag for computations on the model grid
       i_fwo_comm,      &  ! Timing flag for MPI-communications
       i_fwo_ongeom,    &  ! Timing flag for the ray tracing in online beam propagation
       i_fwo_comppolar, &  ! Timing flag for interpolation of reflectivity and radial wind
                             !  from model grid points to the radar bins/auxiliary azi slice grid
       i_fwo_out,       &  ! Timing flag for output (collecting simulated data on one PE per station, sorting, ASCII-output,
                             !  reading obs data, producing feedback files)
       i_fwo_barrier,   &  ! Timing flag for barrier waiting in MPI-communications (measure for load imbalance)
       r_earth_dp,         & ! mean radius of the earth
       rho_w_model, rho_ice_model, K_w_model, K_ice_model, t0_melt_model, itype_gscp_model,        &
       pi,             &
       degrad,         & ! factor for transforming degree to rad
       raddeg,         & ! factor for transforming rad to degree
       fdbk_meta_type, &
       composite_meta_type, &
       list_domains_for_radar, ndoms_max, &
       bubble_list_type, bubble_list, nautobubbles_max

  USE radar_data_namelist, ONLY :  &
       ldebug_radsim,              &
       htop,                       & ! Maximum height MSL for radar computations ( <= model domain top height )
       supob_cart_resolution,      & ! Grid spacing of the superobservation grid in m
       bubble_type, bubble_heatingrate, bubble_timespan, bubble_dT, &
       bubble_centz, bubble_radx, bubble_rady, bubble_radz, bubble_rotangle, &
       bubble_dT_noise, bubble_holdrhconst, bubble_addnoise_T

  ! .. Local domain bounds for computations of grid point values on the model grid, 
  !      e.g., radar_sb_mie_vec(), ..., in radar_mie_iface_cosmo.f90:
  USE radar_data_mie, ONLY : &
       ilow_modelgrid, iup_modelgrid, jlow_modelgrid, jup_modelgrid, klow_modelgrid, kup_modelgrid, &
       rain, cloud, snow, ice, graupel, hail, rain_coeffs, &
       Tmax_i_modelgrid, Tmax_s_modelgrid, Tmax_g_modelgrid, Tmax_h_modelgrid, &
       lgsp_fwo, itype_gscp_fwo, pi6 => pi6_dp

  USE radar_parallel_utilities, ONLY :  &
       global_values_radar, distribute_values_radar, distribute_path_radar

       
  USE radar_utilities, ONLY : &
       smoother, geo_heading, geo_dist, get_range2, polar2geo_xy, &
       el_loc_43, &
       ind2sub2D, sub2ind2D, &
       ind2sub3D, sub2ind3D, &
       ind2sub4D, sub2ind4D, &
       ind2sub5D, sub2ind5D, &
       phirot2phi,           &
       rlarot2rla,           &
       phi2phirot,           &
       rla2rlarot,           &
       get_utc_date,         &
       init_vari

#endif

!==============================================================================

#ifndef NOMPI
  USE mpi
#endif
  
!==============================================================================

  IMPLICIT NONE

!==============================================================================

#ifdef NOMPI
  INCLUDE "nompi_mpif.h"
#endif

!==============================================================================

! default private
  PRIVATE

!==============================================================================

#ifndef HAVE_RADARFWO

!**************************************************************************************************
!**************************************************************************************************
!
! Dummies for global EMVORADO variables when not using the EMVORADO modules,
!  to enable ICON compilation without EMVORADO.
! This compilation mode should be possible for the sole purpose to make sure that needed ICON
!  functionality for EMVORADO cannot be accidentially deleted by
!  other model developers.
!
  
  INTEGER, PARAMETER :: ndoms_max = 5
  INTEGER, PARAMETER :: cmaxlen   = 30
  INTEGER, PARAMETER :: my_radar_id = 0
  INTEGER, PARAMETER :: my_cart_id_fwo = 0
  INTEGER, PARAMETER :: num_compute_fwo = 0
  LOGICAL            :: lcompute_pe_fwo
  REAL(kind=dp)      :: r_earth_dp
  REAL(kind=dp)      :: raddeg, degrad
  INTEGER            :: i_fwo_prep_compute, &
                        i_fwo_bubbles      , &
                        i_fwo_composites   , &
                        i_fwo_ini          , &
                        i_fwo_compgrid     , &
                        i_fwo_comm         , &
                        i_fwo_ongeom       , &
                        i_fwo_comppolar    , &
                        i_fwo_out          , &
                        i_fwo_barrier
  INTEGER            :: ie_fwo, je_fwo, ke_fwo, itype_gscp_fwo
  REAL(kind=wp)      :: K_ice_model, K_w_model, rho_ice_model, rho_w_model, t0_melt_model
  logical            :: lgsp_fwo

! End of dummies
!
!**************************************************************************************************
!**************************************************************************************************

#endif

!==============================================================================

#ifdef HAVE_RADARFWO
  REAL(kind=dp) :: &
       dom_center_lon(ndoms_max) = HUGE(1.0_dp), & ! Centers of radar-active domains
       dom_center_lat(ndoms_max) = HUGE(1.0_dp)    ! Centers of radar-active domains
#endif
  
  REAL(kind=wp) :: &
       pi_wp=4.0_wp*ATAN(1.0_wp)                   ! pi in ICON working precision

  REAL(kind=dp) :: &
       pi_dp=4.0_dp*ATAN(1.0_dp)                   ! pi in double precision

#ifdef HAVE_RADARFWO
  ! meta data for NetCDF feedback files ("fof") for each radar-active domain:
  TYPE(fdbk_meta_type) ::   fdbk_meta_container(ndoms_max)
#endif

!==============================================================================

  ! Type for the auxiliary grid to hold the nearest ICON cell center index information:
  TYPE t_cindex_grid

    INTEGER       :: nlon_tot     ! no. grid points for aux. grid in bounding box (BB) of total radar-covered region
    INTEGER       :: nlat_tot
    INTEGER       :: nlon         ! no. grid points for aux. grid in BB of radar-covered region on local PE
    INTEGER       :: nlat
    REAL(kind=dp) :: pollon       ! rotated pole of both total and local grids for BBs
    REAL(kind=dp) :: pollat
    REAL(kind=dp) :: polgam
    REAL(kind=dp) :: dlon         ! grid spacing of both total and local grids for BBs
    REAL(kind=dp) :: dlat
    REAL(kind=dp) :: startlon_tot ! lower left corner of grid for total BB
    REAL(kind=dp) :: startlat_tot
    REAL(kind=dp) :: startlon     ! lower left corner of grid for BB on local PE
    REAL(kind=dp) :: startlat
    REAL(kind=dp) :: r_earth      ! Earth's radius of grids for both global and local BBs

    INTEGER,       ALLOCATABLE, DIMENSION(:,:) :: cind      ! index of nearest ICON cell centroid in aux. grid on local PE
    INTEGER,       ALLOCATABLE, DIMENSION(:,:) :: cind_glob ! global index of nearest ICON cell centroid in aux. grid on local PE
    INTEGER,       ALLOCATABLE, DIMENSION(:,:) :: idx       ! nproma-index of nearest ICON cell centroid in aux. grid on local PE
    INTEGER,       ALLOCATABLE, DIMENSION(:,:) :: blk       ! nblk-index of nearest ICON cell centroid in aux. grid on local PE
    REAL(kind=dp), ALLOCATABLE, DIMENSION(:,:) :: dist      ! distance to nearest ICON cell centroid

  END TYPE t_cindex_grid

  TYPE(t_cindex_grid) :: cindex_grid(ndoms_max)

  ! Setup of the vectors for storing model index informations and interpolation weights,
  !  which is model specific. ICON with nearest neighbour only needs one weight, which is for the vertical interpolation:
  INTEGER, PARAMETER :: ninterp_weights_model = 1

!==============================================================================

  ! maximum allowed number of model domains in an ICON run:
  INTEGER, PARAMETER :: ndoms_max_model = max_dom

  ! Dynamical variables and other model variables:
  REAL(KIND=wp), POINTER :: &
    hhl      (:,:,:)   => NULL(),   & ! hhl half levels
    hfl      (:,:,:)   => NULL(),   & ! hfl full levels
    u        (:,:,:)   => NULL(),   & ! U at nnew
    v        (:,:,:)   => NULL(),   & ! V at nnew
    w        (:,:,:)   => NULL(),   & ! W at nnew
    t        (:,:,:)   => NULL(),   & ! T at nnew
    p        (:,:,:)   => NULL(),   & ! P at nnew
    rho      (:,:,:)   => NULL(),   & ! RHO at nnew
    rho0     (:,:,:)   => NULL()      ! RHO0

  ! Microphysics tracers of EMVORADO, which are simply retrieved as pointers 
  !  to COSMO tracers and are of same precision as in the COSMO-model:
  REAL(KIND=wp), POINTER :: &
    qv       (:,:,:)   => NULL(),   & ! QV
    qc       (:,:,:)   => NULL(),   & ! QC
    qi       (:,:,:)   => NULL(),   & ! QI
    qr       (:,:,:)   => NULL(),   & ! QR
    qs       (:,:,:)   => NULL(),   & ! QS
    qg       (:,:,:)   => NULL(),   & ! QG
    qnc_s    (:,:)     => NULL()      ! QNC_S (near surface value)

  ! dummy memory for pointing some undefined hydrometeors to,
  !  either in init stage, when p_nh_state is not yet allocated,
  !  or for qgl and qhl in case of 2-moment scheme without liquid water fraction:

  TYPE t_dummy
    REAL(KIND=wp), POINTER :: dummy0(:,:,:) => NULL()
    REAL(KIND=wp), POINTER :: qnc_s (:,:)   => NULL()
  END TYPE t_dummy

  TYPE(t_dummy), DIMENSION(ndoms_max_model) :: dum_dom
  
  ! Other fields for the operator, which should have the
  !  dimension of the model fields: Will be allocated in alloc_aux_model_variables():
  REAL(KIND=dp), ALLOCATABLE :: &
       rlat     (:,:)    ,   & ! rlat
       rlon     (:,:)    ,   & ! rlon
       z_radar(:,:,:), z_ext(:,:,:), vt_radar(:,:,:)

  REAL(KIND=wp), ALLOCATABLE :: &
       vapor_pres(:,:,:)

  LOGICAL :: lalloc_qi, lalloc_qs, lalloc_qg, lalloc_qh

  REAL(KIND=wp), POINTER :: &
    qh        (:,:,:)  => NULL(),   & ! QH
    qnc       (:,:,:)  => NULL(),   & ! NCCLOUD
    qni       (:,:,:)  => NULL(),   & ! NCICE
    qnr       (:,:,:)  => NULL(),   & ! NCRRAIN
    qns       (:,:,:)  => NULL(),   & ! NCSNOW
    qng       (:,:,:)  => NULL(),   & ! NCGRAUPEL
    qnh       (:,:,:)  => NULL(),   & ! NCHAIL
    qgl       (:,:,:)  => NULL(),   & ! QGL
    qhl       (:,:,:)  => NULL()      ! QHL

  ! geometric information for hydrometeor testpattern (ltestpattern_hydrometeors=.TRUE.):
  TYPE(t_patch), POINTER :: p_patch_for_testpattern

  !==============================================================================
  !
  ! Provisional local bubble parameter vectors. Have to be replaced by ICON
  !  parameters once the bubble generator is implemented there!

#ifdef HAVE_RADARFWO
  
  LOGICAL ::    ltempdist(nautobubbles_max), &  ! Switch(es) (up to 50) to release temperature disturbances  (LOGICAL)
       ladd_bubblenoise_t(nautobubbles_max), &  ! Switch(es) (up to 50) to overlay temperature disturbances (LOGICAL)
       lbub_rhconst(nautobubbles_max)           ! For each bubble, specify const. rel. hum. during temperature disturbance (LOGICAL)

  INTEGER :: ntstep_bubble(nautobubbles_max), & ! Beginning of bubble (in time steps)
       bub_zeitzaehler(nautobubbles_max)        ! Counter of bubbles. it must be set to zero before each start

  CHARACTER(len=12) :: ctype_tempdist(nautobubbles_max)   ! Type of temperature disturbance(s)

  REAL(kind=wp), DIMENSION(nautobubbles_max) :: &
       htempdist,                        &  ! Time for release of temperature disturbances [h]
       bub_centlon,                      &  ! Center (i) of temperature disturbances [grid points, real]
       bub_centlat,                      &  ! Center (j) of temperature disturbances [grid points, real]
       bub_centz,                        &  ! Center (Z) of temperature disturbances [m]
       bub_timespan,                     &  ! Duration of release of temperature disturbance(s) [# time steps]
       bub_radx,                         &  ! Length scale (X) of temperature disturbances [m]
       bub_rady,                         &  ! Length scale (Y) of temperature disturbances [m]
       bub_radz,                         &  ! Length scale (Z) of temperature disturbances [m]
       bub_rotangle,                     &  ! Rotation angle of main axes of temperature disturbances [degrees]
       bub_heatingrate,                  &  ! Constant heating rate [K*s**-1] for cos-hrd bubble
       bub_dT,                           &  ! Temperature increment of cos bubble [K]
       bub_dT_bubblenoise                   ! Rel. amplitude of the noise on the temperature disturbance(s) [K]

#endif

!==============================================================================

  PUBLIC abort_run, setup_auxgrid_for_cellindex, diagnose_and_store_uv_pres_temp, &
         get_model_variables, get_model_hydrometeors, get_model_config_for_radar, &
         setup_runtime_timings, get_runtime_timings
  
#ifdef HAVE_RADARFWO
  PUBLIC get_model_time_sec, get_model_inputdir, get_model_outputdir,               &
         get_datetime_ini, get_datetime_act, get_loc_domain,                        &
         bottomlevel, bottomlevel_stag, toplevel, levelincr, get_model_top_height,  &
         get_fdbk_metadata, get_composite_metadata,                                 &
         set_testpattern_hydrometeors_mg, initialize_tmax_1mom_vec_par, finalize_tmax,   &
         initialize_tmax_2mom_vec_par,                                                &
         get_obstime_ind_of_currtime, get_obs_time_tolerance,                       &
         check_obstime_within_forecast, check_obstime_within_modelrun,              &
         check_if_currtime_is_obstime, get_domain_starttime_in_sec,                 &
         alloc_aux_model_variables, dealloc_aux_model_variables,                    &
         it_is_time_for_radar, num_regular_obstimes,                                &
         it_is_time_for_bubblecheck,                                                &
         grid_length_model, one_level_up, one_level_down, ndoms_max_model,          &
         get_dbz3dlin_with_model_method_1mom, get_dbz3dlin_with_model_method_2mom

  PUBLIC t, rho, qv, qc, qi, qr, qs ,qg, qh,            &
         qnc, qni, qnr, qns ,qng, qnh, qgl, qhl, qnc_s, &
         rlat, rlon, hhl, rho0, u, v, w, p,             &
         lalloc_qi, lalloc_qs, lalloc_qg, lalloc_qh,    &
         z_radar, z_ext, vt_radar, hfl, vapor_pres,     &
         pi_wp

  PUBLIC geo2model_cellindex, geo2model_coord_domaincheck, get_lonlat_domain_center, &
         setup_model2radarbins_vec, calc_vert_weight_vec, setup_model2azislices_vec, &
         interp_model2radarbins_scalar, interp2d_model2radarbins_scalar,             &
         interp_model2radarbins_vr,                                                  &
         interp_model2azislices_scalar, interp_model2azislices_vr, calc_vert_weight, &
         interp2D_model2geo_horiz_scalar, interp3D_model2geo_scalar,                 &
         get_domaincenter_global,                                                    &
         get_rotlatlon_domain_for_superobing, set_fdbk_metadata

#ifdef GRIBAPI
  PUBLIC grib2_add_modelspec_info
#endif

  PUBLIC trigger_warm_bubbles
#endif

!==============================================================================
! Interface Blocks for overloaded procedures:
!==============================================================================

#ifdef HAVE_RADARFWO
  
INTERFACE geo2model_cellindex
  MODULE PROCEDURE            &
       geo2model_cellindex_int,           &
       geo2model_cellindex_real
END INTERFACE geo2model_cellindex

INTERFACE check_obstime_within_forecast
  MODULE PROCEDURE    &
       check_obstime_within_fcst_scal, &
       check_obstime_within_fcst_vec
END INTERFACE

INTERFACE check_obstime_within_modelrun
  MODULE PROCEDURE    &
       check_obstime_within_run_scal, &
       check_obstime_within_run_vec
END INTERFACE

#endif

!==============================================================================
!==============================================================================

CONTAINS

!==============================================================================
!==============================================================================



!**************************************************************************************************
!**************************************************************************************************
!
! THE FOLLOWING ROUTINES ARE ALSO COMPILED WHEN NOT USING EMVORADO,
!  IN ORDER TO MAKE SURE THAT ICON VARIABLES AND FUNCTIONALITY NEEDED FOR EMVORADO
!  CANNOT ACCIDENTIALLY BE DELETED FROM THE ICON CODE BY OTHER DEVELOPERS.
!
!**************************************************************************************************
!**************************************************************************************************

  !==============================================================================
  !
  ! Emulate COSMOs model_abort for EMVORADO, for mo_fdbk.f90 and mo_fdbk_io.f90:
  !
  !==============================================================================

  SUBROUTINE abort_run (my_id, ierrorcode, errorstring, routine, mpi_error)
    
    ! Parameter list:
    
    INTEGER ,INTENT(IN)   ::                                     &
         my_id,        & ! id of this processor
         ierrorcode      ! self-defined integer code of the error detected
    
    CHARACTER (LEN=*), INTENT(IN) ::     &
         errorstring   ! self-defined error message
    CHARACTER (LEN=*), INTENT(IN) ::     &
         routine       ! calling routine
    INTEGER , OPTIONAL, INTENT(IN)   ::                          &
         mpi_error   ! error-code of the message passing library
    
    CHARACTER(len=250) :: loc_name
    CHARACTER(len=250) :: loc_err
    INTEGER            :: nzlen, nzerrcode

    loc_name(:) = ' '
    WRITE(loc_name, '(a,i4.4)') ' on proc ', my_id
    loc_name = TRIM(routine)//TRIM(loc_name)

    loc_err(:) = ' '
    IF (PRESENT(mpi_error)) THEN
      ! this is parallel mode
#ifndef NOMPI
      CALL MPI_ERROR_STRING (mpi_error, loc_err, nzlen, nzerrcode)
#endif
      loc_err = TRIM(errorstring)//' : '//TRIM(loc_err)
    ELSE
      loc_err = TRIM(errorstring)
    END IF

    IF (my_id == 0 .OR. my_process_is_mpi_workroot()) THEN
      WRITE (*,*) REPEAT('*', 80)
      WRITE (*,'(a,":",/,a,/,a,i7)') TRIM(loc_name), TRIM(loc_err), 'error code = ', ierrorcode
      WRITE (*,*) REPEAT('*', 80)
    END IF

    CALL finish(TRIM(loc_name), TRIM(loc_err), ierrorcode)

  END SUBROUTINE abort_run

  
  !============================================================================
  ! 
  ! Subroutine for getting the necessary configuration variables of the model.
  ! Has to be called during model initialization stage.
  ! This has to be done also for luse_radarfwo = .FALSE., because it also
  ! concerns the "normal" DBZ gridpoint output (calc_dbz_vec()), which can be
  ! requested in the model GRIBOUT namelists even without using the "rest"
  ! of the radar forward operator. 
  ! The routine has to be called BEFORE the GRIBOUT-namelist(s) is/are read,
  ! but after the input of the basic domain definition and MPI parameters.
  !
  ! Method: Initialize parameters in radar_data.f90, which have to do with the model
  !  configuration (cartesian grid, PEs, MPI), with the respective parameters of the COSMO-model.
  !
  ! If called more than once, most actions are done only during in the first call!
  !
  ! Exceptions are the timing flags, which are initialized in every call, because
  !  in COSMO, the first call is too early for them to be properly initialized.
  !
  !============================================================================
  
  SUBROUTINE get_model_config_for_radar ( idom )

    IMPLICIT NONE

    INTEGER, INTENT(in)    :: idom

    CHARACTER(len=80)      :: yzroutine
    CHARACTER(len=cmaxlen) :: yerrmsg

    yzroutine(:) = ' '
    yzroutine    = 'get_model_config_for_radar'

    ! Set up domain decomposition of COSMO:
    IF (lcompute_pe_fwo) THEN
      ie_fwo = nproma
      je_fwo = p_patch(idom)%nlev
      ke_fwo = p_patch(idom)%nblks_c
    ELSE
      ! not needed on pure io PEs
      ie_fwo = -1
      je_fwo = -1
      ke_fwo = -1
    END IF
           
    ! Set up parameters for the reflectivity calculation on the model grid
    !  which are related to the cloud microphysics parameterization or output level index:
    lgsp_fwo       = ( atm_phy_nwp_config(idom)%inwp_gscp > 0 )

    itype_gscp_model = atm_phy_nwp_config(idom)%inwp_gscp
    
    SELECT CASE(itype_gscp_model)
    CASE (0)
      ! No microphysics in this case, we just need a dummy value here:
      itype_gscp_fwo = 3
    CASE (1)
      itype_gscp_fwo = 3
    CASE (2)
      itype_gscp_fwo = 4
    CASE (4,5,6,7)
      ! 2-moment scheme: subtype not important, just a number >= 2000 for EMVORADO
      itype_gscp_fwo = 2001
    CASE default
      yerrmsg(:) = ' '
      WRITE (yerrmsg,'(a,i3)') 'Error inwp_gscp: scheme not implemented in EMVORADO: inwp_gscp = ', &
           itype_gscp_model
      CALL abort_run(my_radar_id, 79432, yerrmsg, yzroutine)
    END SELECT
    
    ! Set up the timer IDs. This extra call is necessary in case of stand-alone
    !  use in calc_dbz_vec(). Otherwise it is called at the beginning of organize_radar('init').
    CALL setup_runtime_timings ()

    ! Some necessary constants retrieved from the model:
    r_earth_dp    = earth_radius ! radius of spherical earth
    rho_w_model   = rhoh2o       ! dens. of pure water (in wp)
    rho_ice_model = rhoice       ! dens. of pure ice (in wp)
    K_w_model     = K_w_0        ! dielectric constant for water in the microwave region (in wp)
    K_ice_model   = K_i_0        ! dielectric constant for ice   in the microwave region (in wp)
    t0_melt_model = Tmelt        ! melting temperature (in wp)
    
  END SUBROUTINE get_model_config_for_radar

  
  !============================================================================
  ! 
  ! Subroutine for getting the microphysics tracers.
  ! 
  !============================================================================
  
  SUBROUTINE get_model_hydrometeors (ntlev, idom)

    !--------------------------------------------------------------------------
    !
    ! Description:
    !   This subroutine retrieves the globally needed microphysics tracers
    !   out of the tracer structure of COSMO and returns pointers to the 3D fields
    !
    !--------------------------------------------------------------------------

    IMPLICIT NONE

    ! Parameters
    INTEGER, INTENT(in) :: ntlev, idom

    ! Local variables:
    CHARACTER(LEN=25)       :: yzroutine = 'get_model_hydrometeors'
    CHARACTER(LEN=cmaxlen)  :: yerrmsg
    INTEGER                 :: ierror, ni, nk, nkp1, nj

    ierror  = 0
    yerrmsg = ' '

    ! field dimensions: (Note: here nk is the number of height levels, but in other
    !  parts of emvorado this might be nj instead. nk might there be the number
    !  of blocks.
    ni   = nproma
    nk   = p_patch(idom)%nlev
    nkp1 = p_patch(idom)%nlevp1
    nj   = p_patch(idom)%nblks_c

    ! Clean up pointers before new association:
    IF (ASSOCIATED(qv )) NULLIFY(qv )
    IF (ASSOCIATED(qc )) NULLIFY(qc )
    IF (ASSOCIATED(qr )) NULLIFY(qr )
    IF (ASSOCIATED(qi )) NULLIFY(qi )
    IF (ASSOCIATED(qs )) NULLIFY(qs )
    IF (ASSOCIATED(qg )) NULLIFY(qg )
    IF (ASSOCIATED(qh )) NULLIFY(qh )
    IF (ASSOCIATED(qnc)) NULLIFY(qnc)
    IF (ASSOCIATED(qnr)) NULLIFY(qnr)
    IF (ASSOCIATED(qni)) NULLIFY(qni)
    IF (ASSOCIATED(qns)) NULLIFY(qns)
    IF (ASSOCIATED(qng)) NULLIFY(qng)
    IF (ASSOCIATED(qnh)) NULLIFY(qnh)
    IF (ASSOCIATED(qgl)) NULLIFY(qgl)
    IF (ASSOCIATED(qhl)) NULLIFY(qhl)
    IF (ASSOCIATED(qnc_s)) NULLIFY(qnc_s)
    
    IF (ALLOCATED(p_nh_state)) THEN  
      
      qv  => p_nh_state(idom)%prog(ntlev)%tracer(:,:,:,iqv)
      qc  => p_nh_state(idom)%prog(ntlev)%tracer(:,:,:,iqc)
      qr  => p_nh_state(idom)%prog(ntlev)%tracer(:,:,:,iqr)
      qi  => p_nh_state(idom)%prog(ntlev)%tracer(:,:,:,iqi)
      lalloc_qi = .TRUE.
      qs  => p_nh_state(idom)%prog(ntlev)%tracer(:,:,:,iqs)
      lalloc_qs = .TRUE.
      IF (iqg > 0 .AND. iqg <= ntracer) THEN
        qg  => p_nh_state(idom)%prog(ntlev)%tracer(:,:,:,iqg)
        lalloc_qg = .TRUE.
      ELSE
        lalloc_qg = .FALSE.
      END IF
      
      IF (iqh > 0 .AND. iqh <= ntracer) THEN
        qh  => p_nh_state(idom)%prog(ntlev)%tracer(:,:,:,iqh)
        lalloc_qh = .TRUE.
      ELSE
        lalloc_qh = .FALSE.
      END IF
      IF (iqnc > 0 .AND. iqnc <= ntracer) THEN
        qnc => p_nh_state(idom)%prog(ntlev)%tracer(:,:,:,iqnc)
      END IF
      IF (iqnr > 0 .AND. iqnr <= ntracer) THEN
        qnr => p_nh_state(idom)%prog(ntlev)%tracer(:,:,:,iqnr)
      END IF
      IF (iqni > 0 .AND. iqni <= ntracer) THEN
        qni => p_nh_state(idom)%prog(ntlev)%tracer(:,:,:,iqni)
      END IF
      IF (iqns > 0 .AND. iqns <= ntracer) THEN
        qns => p_nh_state(idom)%prog(ntlev)%tracer(:,:,:,iqns)
      END IF
      IF (iqng > 0 .AND. iqng <= ntracer) THEN
        qng => p_nh_state(idom)%prog(ntlev)%tracer(:,:,:,iqng)
      END IF
      IF (iqnh > 0 .AND. iqnh <= ntracer) THEN
        qnh => p_nh_state(idom)%prog(ntlev)%tracer(:,:,:,iqnh)
      END IF
!      IF (iqgl > 0 .AND. iqgl <= ntracer) THEN
!        ! 2mom scheme with liquid water fraction of graupel qgl:
!        qgl => p_nh_state(idom)%prog(ntlev)%tracer(:,:,:,iqgl)
!      ELSE IF (iqng > 0 .AND. iqng <= ntracer) THEN
        ! 2mom scheme without lwf of graupel. qgl needs to be allocated
        ! with 0.0:
        IF (.NOT. ASSOCIATED(dum_dom(idom)%dummy0)) THEN
          ALLOCATE(dum_dom(idom)%dummy0(ni,nkp1,nj))
          dum_dom(idom)%dummy0 = 0.0_wp
        END IF
        qgl => dum_dom(idom)%dummy0(:,1:nk,:)
!      END IF
!      IF (iqhl > 0 .AND. iqhl <= ntracer) THEN
!        ! 2mom scheme with liquid water fraction of hail qhl:
!        qhl => p_nh_state(idom)%prog(ntlev)%tracer(:,:,:,iqhl)
!      ELSE IF (iqnh > 0 .AND. iqnh <= ntracer) THEN
        ! 2mom scheme without lwf of hail. qhl needs to be allocated
        ! with 0.0:
        IF (.NOT. ASSOCIATED(dum_dom(idom)%dummy0)) THEN
          ALLOCATE(dum_dom(idom)%dummy0(ni,nkp1,nj))
          dum_dom(idom)%dummy0 = 0.0_wp
        END IF
        qhl => dum_dom(idom)%dummy0(:,1:nk,:)
!      END IF
      
      IF (atm_phy_nwp_config(idom)%icpl_aero_gscp == 2) THEN
        ! Not yet implemented in microphysics! We give a dummy value here.
        IF (.NOT. ASSOCIATED(dum_dom(idom)%qnc_s)) THEN
          ALLOCATE(dum_dom(idom)%qnc_s(ni,nj)) 
          dum_dom(idom)%qnc_s(:,:) = cloud_num               ! 1/kg
        END IF
        qnc_s => dum_dom(idom)%qnc_s(:,:)
      ELSE IF (atm_phy_nwp_config(idom)%icpl_aero_gscp == 1) THEN
        qnc_s => prm_diag(idom)%cloud_num(:,:) ! neglect difference of 1/m^3 and 1/kg for this near-surface value
      ELSE
        IF (.NOT. ASSOCIATED(dum_dom(idom)%qnc_s)) THEN
          ALLOCATE(dum_dom(idom)%qnc_s(ni,nj)) 
          dum_dom(idom)%qnc_s(:,:) = cloud_num               ! 1/kg
        END IF
        qnc_s => dum_dom(idom)%qnc_s(:,:)
      END IF

      ! Note: computations include halos and boundaries for consistency to EMVORADO internal computations
      
      IF (itype_gscp_fwo >= 100) THEN
        CALL init_2mom_types ()
      ELSE
        CALL init_1mom_types (itype_gscp_fwo)
      END IF
      IF (ldebug_radsim .AND. my_radar_id == 0) THEN
        WRITE (*,'(a)') 'INFO EMVORADO '//TRIM(yzroutine)//':'
        WRITE (*,'(T10,a,T30,a)') 'Cloud type:'  , TRIM(cloud%name)
        WRITE (*,'(T10,a,T30,a)') 'Rain type:'   , TRIM(rain%name)
        WRITE (*,'(T10,a,T30,a)') 'Ice type:'    , TRIM(ice%name)
        WRITE (*,'(T10,a,T30,a)') 'Snow type:'   , TRIM(snow%name)
        WRITE (*,'(T10,a,T30,a)') 'Graupel type:', TRIM(graupel%name)
      END IF

    ELSE

      ! This is necessary, because this routine is called also at a very early
      !  initialization state in ICON by another emvorado routine, and this should
      !  not lead to a crash.
      
      IF (.NOT. ASSOCIATED(dum_dom(idom)%dummy0)) THEN
        ALLOCATE(dum_dom(idom)%dummy0(ni,nkp1,nj))
        dum_dom(idom)%dummy0 = 0.0_wp
      END IF
      qv  => dum_dom(idom)%dummy0(:,1:nk,:)
      qc  => dum_dom(idom)%dummy0(:,1:nk,:)
      qr  => dum_dom(idom)%dummy0(:,1:nk,:)
      qi  => dum_dom(idom)%dummy0(:,1:nk,:)
      lalloc_qi = .TRUE.
      qs  => dum_dom(idom)%dummy0(:,1:nk,:)
      IF (iqg > 0 .AND. iqg <= ntracer) THEN
        qg  => dum_dom(idom)%dummy0(:,1:nk,:)
        lalloc_qg = .TRUE.
      ELSE
        lalloc_qg = .FALSE.
      END IF
      IF (iqh > 0 .AND. iqh <= ntracer) THEN
        qh  => dum_dom(idom)%dummy0(:,1:nk,:)
        lalloc_qh = .TRUE.
      ELSE
        lalloc_qh = .FALSE.
      END IF
      IF (iqnc > 0 .AND. iqnc <= ntracer) THEN
        qnc => dum_dom(idom)%dummy0(:,1:nk,:)
      END IF
      IF (iqnr > 0 .AND. iqnr <= ntracer) THEN
        qnr => dum_dom(idom)%dummy0(:,1:nk,:)
      END IF
      IF (iqni > 0 .AND. iqni <= ntracer) THEN
        qni => dum_dom(idom)%dummy0(:,1:nk,:)
      END IF
      IF (iqns > 0 .AND. iqns <= ntracer) THEN
        qns => dum_dom(idom)%dummy0(:,1:nk,:)
      END IF
      IF (iqng > 0 .AND. iqng <= ntracer) THEN
        qng => dum_dom(idom)%dummy0(:,1:nk,:)
        qgl => dum_dom(idom)%dummy0(:,1:nk,:)
      END IF
      IF (iqnh > 0 .AND. iqnh <= ntracer) THEN
        qnh => dum_dom(idom)%dummy0(:,1:nk,:)
        qhl => dum_dom(idom)%dummy0(:,1:nk,:)
      END IF
      qnc_s => dum_dom(idom)%dummy0(:,1,:)
      
    END IF

  END SUBROUTINE get_model_hydrometeors

  !------------------------------------------------------------------------------

  !============================================================================
  ! 
  ! Subroutine for initializing derived types for the hydrometeor parameters
  !  of the 2-moment scheme
  !
  !============================================================================

  SUBROUTINE init_2mom_types ()
    IMPLICIT NONE

    TYPE(particle)        :: cloud_2mom, rain_2mom
    TYPE(particle_frozen) :: ice_2mom, snow_2mom, graupel_2mom, hail_2mom
!!$ LWF scheme not yet implemented:    CLASS(particle_lwf)    :: graupel_lwf_2mom, hail_lwf_2mom

    CALL init_2mom_scheme(cloud_2mom,rain_2mom,ice_2mom,snow_2mom,graupel_2mom,hail_2mom)

#ifdef HAVE_RADARFWO
    
    !------------------------------------------------------------------------------
    !
    ! Changed EMVORADO-side MGD parameters mu/nu from Seifert mass-based notation used in
    ! COSMO/ICON 2-mom scheme to general diameter-based notation:
    !
    ! (Seifert) N(m) = N0 * m^nu_x * exp(-lam*m^mu_x)  -->
    ! (general) N(D) = N0 * D^mu_D * exp(-lam*D^nu_D)
    !
    ! (NOTE: in ICON/COSMO 2-mom, they remain to be in Seifert notation)
    !
    ! Conversion of mass-size relation parameters a_geo and b_geo from mass-space
    ! (as used in COSMO/ICON 2-mom) to D-space:
    !
    ! (mass-space) D = a_m * m^(b_m) -->
    ! (D-space)    m = a_D * D^(b_D)
    !
    ! related by:
    ! a_D = (1/a_m)^(1/b_m)
    ! b_D = 1/b_m
    !
    !
    ! NOTE: Remaining a&b parameters a/b_vel (mass-fallspeed parameters) and
    !       a/b_ven are still in the original (mass-space) SB notation:
    !
    ! v_T = a_v * m^(b_v)
    !
    ! The switch to D-notation might be done in the future.
    !
    !------------------------------------------------------------------------------

    cloud%name  = cloud_2mom%name
    ! turn around MGD parameters & convert from mu/nu_m to mu/nu_D
    cloud%mu    = (1.0d0/cloud_2mom%b_geo)*(cloud_2mom%nu+1.0d0)-1.0d0
    cloud%nu    = (1.0d0/cloud_2mom%b_geo)*cloud_2mom%mu
    cloud%x_max = cloud_2mom%x_max
    cloud%x_min = cloud_2mom%x_min
    ! convert a/b_geo_m to a/b_geo_D
    cloud%a_geo = (1.0d0/cloud_2mom%a_geo)**(1.0d0/cloud_2mom%b_geo)
    cloud%b_geo = (1.0d0/cloud_2mom%b_geo)
    ! leave those as is for now
    cloud%a_vel = cloud_2mom%a_vel
    cloud%b_vel = cloud_2mom%b_vel
    cloud%a_ven = cloud_2mom%a_ven
    cloud%b_ven = cloud_2mom%b_ven
    ! in 2mom, n0 just a dummy (but needs to be positive)
    cloud%n0_const = 1.0d0

    ice%name  = ice_2mom%name
    ! turn around MGD parameters & convert from mu/nu_m to mu/nu_D
    ice%mu    = (1.0d0/ice_2mom%b_geo)*(ice_2mom%nu+1.0d0)-1.0d0
    ice%nu    = (1.0d0/ice_2mom%b_geo)*ice_2mom%mu
    ice%x_max = ice_2mom%x_max
    ice%x_min = ice_2mom%x_min
    ! convert a/b_geo_m to a/b_geo_D
    ice%a_geo = (1.0d0/ice_2mom%a_geo)**(1.0d0/ice_2mom%b_geo)
    ice%b_geo = (1.0d0/ice_2mom%b_geo)
    ! leave those as is for now
    ice%a_vel = ice_2mom%a_vel
    ice%b_vel = ice_2mom%b_vel
    ice%a_ven = ice_2mom%a_ven
    ice%b_ven = ice_2mom%b_ven
    ! in 2mom, n0 just a dummy (but needs to be positive)
    ice%n0_const = 1.0d0

    rain%name  = rain_2mom%name
    ! turn around MGD parameters & convert from mu/nu_m to mu/nu_D
    rain%mu    = (1.0d0/rain_2mom%b_geo)*(rain_2mom%nu+1.0d0)-1.0d0
    rain%nu    = (1.0d0/rain_2mom%b_geo)*rain_2mom%mu
    rain%x_max = rain_2mom%x_max
    rain%x_min = rain_2mom%x_min
    ! convert a/b_geo_m to a/b_geo_D
    rain%a_geo = (1.0d0/rain_2mom%a_geo)**(1.0d0/rain_2mom%b_geo)
    rain%b_geo = (1.0d0/rain_2mom%b_geo)
    ! leave those as is for now
    rain%a_vel = rain_2mom%a_vel
    rain%b_vel = rain_2mom%b_vel
    rain%a_ven = rain_2mom%a_ven
    rain%b_ven = rain_2mom%b_ven
    ! in 2mom, n0 just a dummy (but needs to be positive)
    rain%n0_const = 1.0d0

    snow%name  = snow_2mom%name
    ! turn around MGD parameters & convert from mu/nu_m to mu/nu_D
    snow%mu    = (1.0d0/snow_2mom%b_geo)*(snow_2mom%nu+1.0d0)-1.0d0
    snow%nu    = (1.0d0/snow_2mom%b_geo)*snow_2mom%mu
    snow%x_max = snow_2mom%x_max
    snow%x_min = snow_2mom%x_min
    ! convert a/b_geo_m to a/b_geo_D
    snow%a_geo = (1.0d0/snow_2mom%a_geo)**(1.0d0/snow_2mom%b_geo)
    snow%b_geo = (1.0d0/snow_2mom%b_geo)
    ! leave those as is for now
    snow%a_vel = snow_2mom%a_vel
    snow%b_vel = snow_2mom%b_vel
    snow%a_ven = snow_2mom%a_ven
    snow%b_ven = snow_2mom%b_ven
    ! in 2mom, n0 just a dummy (but needs to be positive)
    snow%n0_const = 1.0d0

    graupel%name  = graupel_2mom%name
    ! turn around MGD parameters & convert from mu/nu_m to mu/nu_D
    graupel%mu    = (1.0d0/graupel_2mom%b_geo)*(graupel_2mom%nu+1.0d0)-1.0d0
    graupel%nu    = (1.0d0/graupel_2mom%b_geo)*graupel_2mom%mu
    graupel%x_max = graupel_2mom%x_max
    graupel%x_min = graupel_2mom%x_min
    ! convert a/b_geo_m to a/b_geo_D
    graupel%a_geo = (1.0d0/graupel_2mom%a_geo)**(1.0d0/graupel_2mom%b_geo)
    graupel%b_geo = (1.0d0/graupel_2mom%b_geo)
    ! leave those as is for now
    graupel%a_vel = graupel_2mom%a_vel
    graupel%b_vel = graupel_2mom%b_vel
    graupel%a_ven = graupel_2mom%a_ven
    graupel%b_ven = graupel_2mom%b_ven
    ! in 2mom, n0 just a dummy (but needs to be positive)
    graupel%n0_const = 1.0d0

    hail%name  = hail_2mom%name
    ! turn around MGD parameters & convert from mu/nu_m to mu/nu_D
    hail%mu    = (1.0d0/hail_2mom%b_geo)*(hail_2mom%nu+1.0d0)-1.0d0
    hail%nu    = (1.0d0/hail_2mom%b_geo)*hail_2mom%mu
    hail%x_max = hail_2mom%x_max
    hail%x_min = hail_2mom%x_min
    ! convert a/b_geo_m to a/b_geo_D
    hail%a_geo = (1.0d0/hail_2mom%a_geo)**(1.0d0/hail_2mom%b_geo)
    hail%b_geo = (1.0d0/hail_2mom%b_geo)
    ! leave those as is for now
    hail%a_vel = hail_2mom%a_vel
    hail%b_vel = hail_2mom%b_vel
    hail%a_ven = hail_2mom%a_ven
    hail%b_ven = hail_2mom%b_ven
    ! in 2mom, n0 just a dummy (but needs to be positive)
    hail%n0_const = 1.0d0

    ! No need to rename these; they DO actually refer to D-space mu
    ! (follows Seifert(2008), which applies the general MGD notation).
    ! see also FUNCTION mu_d_relation_seifert in radar_mie_utils.f90.
    rain_coeffs%cmu0 = rain_coeffs_2mom%cmu0
    rain_coeffs%cmu1 = rain_coeffs_2mom%cmu1
    rain_coeffs%cmu2 = rain_coeffs_2mom%cmu2
    rain_coeffs%cmu3 = rain_coeffs_2mom%cmu3
    rain_coeffs%cmu4 = rain_coeffs_2mom%cmu4
    rain_coeffs%cmu5 = rain_coeffs_2mom%cmu5

#endif

  END SUBROUTINE init_2mom_types

!----------------------------------------------------------------------------------------
! Particle definitions for 1-mom cases.
! Settings taken from [vt]hydroparams_1mom subroutine(s).
!
! NOTE:
! - mu and nu for general DSD notation: N(D) = N0 * D^mu * exp(-lam * D^nu)
! - mu, nu, a/b_geo in diameter (D) space (not in mass (x)!)
! - but: a/b_vel in mass (x) space
! - unknown/irrelevant parameters set to fill value

  SUBROUTINE init_1mom_types(itype_gscp)
    IMPLICIT NONE

    INTEGER, INTENT(in) :: itype_gscp

    cloud%name  = 'cloud1mom'            !.name...Bezeichnung der Partikelklasse
    cloud%mu    = 3.0000d0               !.mu.....Breiteparameter der Verteil.
    cloud%nu    = 3.0000d0               !.nu.....Exp.-parameter der Verteil.
    cloud%n0_const = 1.0d0                !.n0.....Scaling parameter of distribution (only set if constant)
    cloud%a_geo = rho_w_model*pi6       !.a_geo..Koeff. Geometrie
    cloud%b_geo = 3.0000d0               !.b_geo..Koeff. Geometrie
    ! here: x_min == x_max == x(D_c=20um)
    ! medial mass of cloud droplets; monodisperse distribution
    cloud%x_max = cloud%a_geo*20d-6**cloud%b_geo
                                        !.x_max..maximale Teilchenmasse
    cloud%x_min = cloud%x_max            !.x_min..minimale Teilchenmasse
    ! Parameters for terminal velocity: vt = a_velD * D^b_velD = a_velx * x^b_velx
    !  (Stokes law, T=273.15 K, dyn. visc. eta after Sutherland-formula, Archimedes buoyancy neglected)
    !   v = g*rho_w/(18*eta) * D^2
    ! first, set them as D-space parameters (as were given in vthydroparams_1mom)
    cloud%a_vel = 3.17204d7             !.a_vel..Koeff. Fallgesetz
    cloud%b_vel = 2.00000d0             !.b_vel..Koeff. Fallgesetz
    ! now convert to x-space (more practical to use in code)
    cloud%b_vel = cloud%b_vel/cloud%b_geo
    cloud%a_vel = cloud%a_vel/cloud%a_geo**cloud%b_vel
    cloud%a_ven = miss_value             !.a_ven..Koeff. Ventilation
    cloud%b_ven = miss_value             !.b_ven..Koeff. Ventilation

    IF (itype_gscp < 4) THEN
      rain%name  = 'rain1mom_gscp.lt.4' !.name...Bezeichnung der Partikelklasse
      rain%mu    = 0.0000d0            !.mu.....Breiteparameter der Verteil.
      rain%n0_const = 8.0d6 * EXP(3.2d0*rain%mu) * (1d-2)**(-rain%mu) !* rain_n0_factor
                                        !.n0.....Scaling parameter of distribution (only set if constant)
    ELSE
      rain%name  = 'rain1mom_gscp.ge.4' !.name...Bezeichnung der Partikelklasse
      rain%mu    = 0.5000d0             !.mu.....Breiteparameter der Verteil.
      rain%n0_const = 8d6 * EXP(3.2d0*rain%mu) * (1d-2)**(-rain%mu) !* rain_n0_factor
                                        !.n0.....Scaling parameter of distribution (only set if constant)
    END IF
    rain%nu    = 1.0000d0               !.nu.....Exp.-parameter der Verteil.
    rain%x_max = miss_value             !.x_max..maximale Teilchenmasse
    rain%x_min = miss_value             !.x_min..minimale Teilchenmasse
    rain%a_geo = rho_w_model*pi6        !.a_geo..Koeff. Geometrie
    rain%b_geo = 3.0000d0               !.b_geo..Koeff. Geometrie
    ! Parameters for terminal velocity: vt = a_velD * D^b_velD = a_velx * x^b_velx
    ! first, set them as D-space parameters (as were given in vthydroparams_1mom)
    rain%a_vel = 130.00d0               !.a_vel..Koeff. Fallgesetz
    rain%b_vel = 0.5000d0  !0.25d0      !.b_vel..Koeff. Fallgesetz
    ! now convert to x-space (more practical to use in code)
    rain%b_vel = rain%b_vel/rain%b_geo
    rain%a_vel = rain%a_vel/rain%a_geo**rain%b_vel
    rain%a_ven = miss_value             !.a_ven..Koeff. Ventilation
    rain%b_ven = miss_value             !.b_ven..Koeff. Ventilation

    ice%name  = 'ice1mom'               !.name...Bezeichnung der Partikelklasse
    ice%mu    = 1.0000d0                !.mu.....Breiteparameter der Verteil.
    ice%nu    = 1.0000d0                !.nu.....Exp.-parameter der Verteil.
    ice%n0_const = 1.0d0                  !.n0.....Scaling parameter of distribution (only set if constant)
    ice%a_geo = 130.00d0                !.a_geo..Koeff. Geometrie
    ice%b_geo = 3.0000d0                !.b_geo..Koeff. Geometrie
    ! here: x_min == x_max == x(D_c=100um)
    ! medial mass of cloud ice crystals; monodisperse distribution
    ice%x_max = ice%a_geo*100d-6**ice%b_geo
                                        !.x_max..maximale Teilchenmasse
    ice%x_min = ice%x_max               !.x_min..minimale Teilchenmasse
    ! Parameters for terminal velocity: vt = a_velD * D^b_velD = a_velx * x^b_velx
    ! first, set them as D-space parameters (as were given in vthydroparams_1mom)
    ice%a_vel = 4.9d0 * 0.75d0  ! the 0.75 is Uli's tuning of the snow relation in the COSMO-Docs
                                        !.a_vel..Koeff. Fallgesetz
    ice%b_vel = 0.2500d0                !.b_vel..Koeff. Fallgesetz
    ! now convert to x-space (more practical to use in code)
    ice%b_vel = ice%b_vel/ice%b_geo
    ice%a_vel = ice%a_vel/ice%a_geo**ice%b_vel
    ice%a_ven = miss_value              !.a_ven..Koeff. Ventilation
    ice%b_ven = miss_value              !.b_ven..Koeff. Ventilation

    IF (itype_gscp == 4) THEN
      snow%name  = 'snow1mom_gscp.eq.4' !.name...Bezeichnung der Partikelklasse
      snow%a_geo = 0.0380d0             !.a_geo..Koeff. Geometrie
    ELSE
      snow%name  = 'snow1mom_gscp.ne.4' !.name...Bezeichnung der Partikelklasse
      snow%a_geo = 0.0690d0             !.a_geo..Koeff. Geometrie
    END IF
    snow%b_geo = 2.0000d0               !.b_geo..Koeff. Geometrie
    snow%mu    = 0.0000d0               !.mu.....Breiteparameter der Verteil.
    snow%nu    = 1.0000d0               !.nu.....Exp.-parameter der Verteil.
    snow%n0_const = 1.0d0                 !.n0.....Scaling parameter of distribution (only set if constant)
    ! x_min/max set from max limits for particle type in 2mom case
    snow%x_max = 2.00d-05               !.x_max..maximale Teilchenmasse
    snow%x_min = 1.00d-12               !.x_min..minimale Teilchenmasse
    ! Parameters for terminal velocity: vt = a_velD * D^b_velD = a_velx * x^b_velx
    ! first, set them as D-space parameters (as were given in vthydroparams_1mom)
    snow%a_vel = 25.000d0               !.a_vel..Koeff. Fallgesetz
    snow%b_vel = 0.5000d0               !.b_vel..Koeff. Fallgesetz
    ! now convert to x-space (more practical to use in code)
    snow%b_vel = snow%b_vel/snow%b_geo
    snow%a_vel = snow%a_vel/snow%a_geo**snow%b_vel
    snow%a_ven = miss_value             !.a_ven..Koeff. Ventilation
    snow%b_ven = miss_value             !.b_ven..Koeff. Ventilation

    graupel%name  = 'graupel1mom'       !.name...Bezeichnung der Partikelklasse
    graupel%mu    = 0.0000d0            !.mu.....Breiteparameter der Verteil.
    graupel%nu    = 1.0000d0            !.nu.....Exp.-parameter der Verteil.
    graupel%n0_const = 4.0d6 ! * 0.1  ! factor 0.1 seems to be too low, maybe 0.2?
                                        !.n0.....Scaling parameter of distribution (only set if constant)
    ! x_min/max set from max limits for particle type in 2mom case
    graupel%x_max = 5.00d-04            !.x_max..maximale Teilchenmasse
    !graupel%x_min = 1.00d-09            !.x_min..minimale Teilchenmasse
    ! adapted to avoid Dref4fmelt-calc in zradar_wetgr_mie_vec be noticably low-cut
    graupel%x_min = 1.00d-10            !.x_min..minimale Teilchenmasse
    graupel%a_geo = 169.60d0            !.a_geo..Koeff. Geometrie
    graupel%b_geo = 3.1000d0            !.b_geo..Koeff. Geometrie
    ! Parameters for terminal velocity: vt = a_velD * D^b_velD = a_velx * x^b_velx
    ! first, set them as D-space parameters (as were given in vthydroparams_1mom)
    graupel%a_vel = 442.00d0            !.a_vel..Koeff. Fallgesetz
    graupel%b_vel = 0.8900d0            !.b_vel..Koeff. Fallgesetz
    ! now convert to x-space (more practical to use in code)
    graupel%b_vel = graupel%b_vel/graupel%b_geo
    graupel%a_vel = graupel%a_vel/graupel%a_geo**graupel%b_vel
    graupel%a_ven = miss_value          !.a_ven..Koeff. Ventilation
    graupel%b_ven = miss_value          !.b_ven..Koeff. Ventilation

    ! no hail in 1mom, but to avoid issues with routines used by 1- and 2-mom, define a rudimentary set.
    hail%name  = 'hail1mom'             !.name...Bezeichnung der Partikelklasse
    hail%mu    = miss_value             !.mu.....Breiteparameter der Verteil.
    hail%nu    = miss_value             !.nu.....Exp.-parameter der Verteil.
    hail%n0_const = 1.0d0                 !.n0.....Scaling parameter of distribution (only set if constant)
    hail%x_max = miss_value             !.x_max..maximale Teilchenmasse
    hail%x_min = miss_value             !.x_min..minimale Teilchenmasse
    hail%a_geo = miss_value             !.a_geo..Koeff. Geometrie
    hail%b_geo = miss_value             !.b_geo..Koeff. Geometrie
    hail%a_vel = miss_value             !.a_vel..Koeff. Fallgesetz
    hail%b_vel = miss_value             !.b_vel..Koeff. Fallgesetz
    hail%a_ven = miss_value             !.a_ven..Koeff. Ventilation
    hail%b_ven = miss_value             !.b_ven..Koeff. Ventilation

  END SUBROUTINE init_1mom_types

  !============================================================================
  ! 
  ! Subroutine for getting  the microphysics tracers.
  ! 
  !============================================================================
  
  SUBROUTINE get_model_variables (ntlev_dyn, ntlev_qx, idom)

    !--------------------------------------------------------------------------
    !
    ! Description:
    !   This subroutine connects the dynamic model variables with their
    !    associated pointers in the forward operator and retrieves some
    !    model constants.
    !
    !--------------------------------------------------------------------------

    IMPLICIT NONE

    ! Parameters
    INTEGER, INTENT(in) :: ntlev_dyn, ntlev_qx, idom

    ! Local variables:
    CHARACTER(LEN=25)       :: yzroutine = 'get_model_variables'
    CHARACTER(LEN=cmaxlen)  :: yerrmsg
    INTEGER                 :: ierror, kbub, ncomp, ni, nk, nkp1, nj, k, nk_outer

    ! For checking if the last call was on the same domain at same time and time levels:
    INTEGER, SAVE           :: ntlev_dyn_lastcall = -HUGE(1), &
                               ntlev_qx_lastcall  = -HUGE(1), &
                               idom_lastcall      = -HUGE(1)
    REAL(kind=dp), SAVE     :: time_mod_lastcall  = -HUGE(1.0_dp)

    ierror  = 0
    yerrmsg = ' '

    ni   = nproma
    nk   = p_patch(idom)%nlev
    nkp1 = p_patch(idom)%nlevp1
    nj   = p_patch(idom)%nblks_c

    ! Clean up pointers before new association:
    IF (ASSOCIATED(hhl )) NULLIFY(hhl )
    IF (ASSOCIATED(hfl )) NULLIFY(hfl )
    IF (ASSOCIATED(u   )) NULLIFY(u   )
    IF (ASSOCIATED(v   )) NULLIFY(v   )
    IF (ASSOCIATED(w   )) NULLIFY(w   )
    IF (ASSOCIATED(t   )) NULLIFY(t   )
    IF (ASSOCIATED(p   )) NULLIFY(p   )
    IF (ASSOCIATED(rho )) NULLIFY(rho )
    IF (ASSOCIATED(rho0)) NULLIFY(rho0)
    
    IF (ALLOCATED(p_nh_state)) THEN

      ! Diagnostic variables might have to be updated first:
      IF (p_patch(idom)%ldom_active) THEN
        IF (get_model_time_sec() /= time_mod_lastcall .OR. idom_lastcall /= idom .OR. &
             ntlev_dyn_lastcall /= ntlev_dyn .OR. ntlev_qx_lastcall /= ntlev_qx) THEN
          CALL diagnose_and_store_uv_pres_temp (ntlev_dyn, ntlev_qx, idom)
          time_mod_lastcall  = get_model_time_sec()
          idom_lastcall      = idom
          ntlev_dyn_lastcall = ntlev_dyn
          ntlev_qx_lastcall  = ntlev_qx
        END IF
      END IF

      ! Link the ICON model prognostic variables and other fields:
      hhl  => p_nh_state(idom)%metrics%z_ifc(:,:,:)
      hfl  => p_nh_state(idom)%metrics%z_mc(:,:,:)
      u    => p_nh_state(idom)%diag%u(:,:,:)
      v    => p_nh_state(idom)%diag%v(:,:,:)
      w    => p_nh_state(idom)%prog(ntlev_dyn)%w(:,:,:)
      t    => p_nh_state(idom)%diag%temp(:,:,:)
      p    => p_nh_state(idom)%diag%pres(:,:,:)
      rho  => p_nh_state(idom)%prog(ntlev_dyn)%rho(:,:,:)
      rho0 => rho(:,:,:)

    ELSE

      IF (.NOT. ASSOCIATED(dum_dom(idom)%dummy0)) ALLOCATE(dum_dom(idom)%dummy0(ni,nkp1,nj))
      dum_dom(idom)%dummy0 = 0.0_wp
      hhl  => dum_dom(idom)%dummy0(:,1:nkp1,:)
      hfl  => dum_dom(idom)%dummy0(:,1:nk,:)
      u    => dum_dom(idom)%dummy0(:,1:nk,:)
      v    => dum_dom(idom)%dummy0(:,1:nk,:)
      w    => dum_dom(idom)%dummy0(:,1:nkp1,:)
      t    => dum_dom(idom)%dummy0(:,1:nk,:)
      p    => dum_dom(idom)%dummy0(:,1:nk,:)
      rho  => dum_dom(idom)%dummy0(:,1:nk,:)
      rho0 => dum_dom(idom)%dummy0(:,1:nk,:)

    END IF

    ! For ltestpattern_hydrometeors=.TRUE.:
    p_patch_for_testpattern => p_patch(idom)

  END SUBROUTINE get_model_variables

  !------------------------------------------------------------------------------

  SUBROUTINE diagnose_and_store_uv_pres_temp (ntlev_dyn, ntlev_qx, jg)

    INTEGER, INTENT(in) :: ntlev_dyn, ntlev_qx, jg

    CALL rbf_vec_interpol_cell ( p_nh_state(jg)%prog(ntlev_dyn)%vn,               &
         &                       p_patch(jg), p_int_state(jg),                    &
         &                       p_nh_state(jg)%diag%u, p_nh_state(jg)%diag%v )

    CALL diagnose_pres_temp ( p_nh_state(jg)%metrics, p_nh_state(jg)%prog(ntlev_dyn), &
         &                    p_nh_state(jg)%prog(ntlev_qx),                      &
         &                    p_nh_state(jg)%diag, p_patch(jg),                   &
         &                    opt_calc_temp=.TRUE.,                               &
         &                    opt_calc_pres=.TRUE. )

  END SUBROUTINE diagnose_and_store_uv_pres_temp

  !============================================================================

  !============================================================================
  !============================================================================

  !============================================================================
  !
  ! Subroutine for computing the linear reflectivity on any 3D grid by the same
  ! intrinsic method of the hosting model (ICON in this case).
  !
  ! Is called in calc_dbz_vec_modelgrid() and calc_dbz_vec_generic()
  !   in radar_mie_iface_cosmo.f90, if dbz%itype_refl = 4
  !
  ! Versions for the 1mom- and 2mom microphysics schemes of ICON
  !
  ! Before calling these routines, the following has to be called:
  !   CALL get_model_config_for_radar (idom)
  !
  !============================================================================
  
  SUBROUTINE get_dbz3dlin_with_model_method_1mom (ni, nj, nk,           &
       itype_gscp_model_in, rho_water, rho_ice, K_water, K_ice, t0melt, &
       t_in, rho_in, qc_in, qr_in, &
       qi_in, qs_in, qnc_s, z_radar, qg_in )

    INTEGER , INTENT(in)                    :: ni, nj, nk
    INTEGER , INTENT(in)                    :: itype_gscp_model_in
    REAL(wp), INTENT(in)                    :: rho_water, rho_ice, K_water, K_ice, t0melt
    REAL(wp), INTENT(in), DIMENSION(:,:,:)  :: t_in, rho_in, qc_in, qr_in, &
                                               qi_in, qs_in
    REAL(wp), INTENT(in), DIMENSION(:,:)    :: qnc_s
    REAL(wp), INTENT(in), OPTIONAL          :: qg_in(:,:,:)
    REAL(wp), INTENT(inout)                 :: z_radar(:,:,:)

    CHARACTER(len=cmaxlen)                  :: yerrmsg
    CHARACTER(len=80)                       :: yzroutine

    yzroutine(:) = ' '
    yzroutine    = 'get_dbz3dlin_with_model_method_1mom'

    SELECT CASE ( itype_gscp_model_in )
      
    CASE ( 1 )

      ! Note: computations include halos and boundaries for consistency to EMVORADO internal computations
      CALL compute_field_dbz_1mom( npr       = ni,                   &
                                   nlev      = nj,                   &
                                   nblks     = nk,                   &
                                   startblk  = 1,                    &
                                   endblk    = nk,                   &
                                   jk_start  = 1,                    &
                                   startidx1 = 1,                    &
                                   endidx2   = ni,                   &
                                   lmessage_light = (ldebug_radsim .AND. my_radar_id == 0), &
                                   lmessage_full  = (ldebug_radsim), &
                                   my_id_for_message = my_radar_id,  &
                                   rho_w     = rho_water,            &
                                   rho_ice   = rho_ice,              &
                                   K_w       = K_water,              &
                                   K_ice     = K_ice,                &
                                   T_melt    = t0melt,               &
                                   igscp     = itype_gscp_model_in,  &
                                   q_crit_radar = 1e-8_wp,           &
                                   T         = t(:,:,:),             &
                                   rho       = rho(:,:,:),           &
                                   q_cloud   = qc_in(:,:,:),         &
                                   q_ice     = qi_in(:,:,:),         &
                                   q_rain    = qr_in(:,:,:),         &
                                   q_snow    = qs_in(:,:,:),         &
                                   n_cloud_s = qnc_s(:,:),           &
                                   z_radar   = z_radar(:,:,:)        )

    CASE ( 2 )

      ! Note: computations include halos and boundaries for consistency to EMVORADO internal computations
      CALL compute_field_dbz_1mom( npr       = ni,                   &
                                   nlev      = nj,                   &
                                   nblks     = nk,                   &
                                   startblk  = 1,                    &
                                   endblk    = nk,                   &
                                   jk_start  = 1,                    &
                                   startidx1 = 1,                    &
                                   endidx2   = ni,                   &
                                   lmessage_light = (ldebug_radsim .AND. my_radar_id == 0), &
                                   lmessage_full  = (ldebug_radsim), &
                                   my_id_for_message = my_radar_id,  &
                                   rho_w     = rho_water,            &
                                   rho_ice   = rho_ice,              &
                                   K_w       = K_water,              &
                                   K_ice     = K_ice,                &
                                   T_melt    = t0melt,               &
                                   igscp     = itype_gscp_model_in,  &
                                   q_crit_radar = 1e-8_wp,           &
                                   T         = t(:,:,:),             &
                                   rho       = rho(:,:,:),           &
                                   q_cloud   = qc_in(:,:,:),         &
                                   q_ice     = qi_in(:,:,:),         &
                                   q_rain    = qr_in(:,:,:),         &
                                   q_snow    = qs_in(:,:,:),         &
                                   q_graupel = qg_in(:,:,:),         &
                                   n_cloud_s = qnc_s(:,:),           &
                                   z_radar   = z_radar(:,:,:)        )
            
    CASE DEFAULT
      
      yerrmsg(:) = ' '
      WRITE (yerrmsg,'(a,i3)') 'Error itype_gscp_model_in: scheme not implemented in EMVORADO: itype_gscp_model_in = ', &
           itype_gscp_model_in
      CALL abort_run(my_radar_id, 35621, yerrmsg, yzroutine)

    END SELECT

    

  END SUBROUTINE get_dbz3dlin_with_model_method_1mom

  !============================================================================

  SUBROUTINE get_dbz3dlin_with_model_method_2mom (ni, nj, nk,           &
       itype_gscp_model_in, rho_water, rho_ice, K_water, K_ice, t0melt, &
       t_in, rho_in, qc_in, qr_in, qi_in, qs_in, qg_in, qh_in,          &
       qnc_in, qnr_in, qni_in, qns_in, qng_in, qnh_in, qgl_in, qhl_in,  &
       z_radar )

!!$ qgl, qhl implementieren
    
    INTEGER , INTENT(in)                    :: ni, nj, nk
    INTEGER , INTENT(in)                    :: itype_gscp_model_in
    REAL(wp), INTENT(in)                    :: rho_water, rho_ice, K_water, K_ice, t0melt
    REAL(wp), INTENT(in), DIMENSION(:,:,:)  :: t_in, rho_in, &
                                               qc_in,  qr_in,  qi_in,  qs_in,  qg_in,  qh_in, &
                                               qnc_in, qnr_in, qni_in, qns_in, qng_in, qnh_in, &
                                               qgl_in, qhl_in
    REAL(wp), INTENT(inout)                 :: z_radar(:,:,:)

    CHARACTER(len=cmaxlen)                  :: yerrmsg
    CHARACTER(len=80)                       :: yzroutine

    yzroutine(:) = ' '
    yzroutine    = 'get_dbz3dlin_with_model_method_2mom'


    SELECT CASE ( itype_gscp_model_in )

    CASE ( 4, 5, 6 )

      ! Note: computations include halos and boundaries for consistency to EMVORADO internal computations
      CALL compute_field_dbz_2mom( npr       = ni,                   &
                                   nlev      = nj,                   &
                                   nblks     = nk,                   &
                                   startblk  = 1,                    &
                                   endblk    = nk,                   &
                                   jk_start  = 1,                    &
                                   startidx1 = 1,                    &
                                   endidx2   = ni,                   &
                                   lmessage_light = (ldebug_radsim .AND. my_radar_id == 0), &
                                   lmessage_full  = (ldebug_radsim), &
                                   my_id_for_message = my_radar_id,  &
                                   rho_w     = rho_water,            &
                                   rho_ice   = rho_ice,              &
                                   K_w       = K_water,              &
                                   K_ice     = K_ice,                &
                                   T_melt    = t0melt,               &
                                   q_crit_radar = 1e-8_wp,           &
                                   T         = t(:,:,:),             &
                                   rho       = rho(:,:,:),           &
                                   q_cloud   = qc_in(:,:,:),         &
                                   q_ice     = qi_in(:,:,:),         &
                                   q_rain    = qr_in(:,:,:),         &
                                   q_snow    = qs_in(:,:,:),         &
                                   q_graupel = qg_in(:,:,:),         &
                                   q_hail    = qh_in(:,:,:),         &
                                   n_cloud   = qnc_in(:,:,:),        &
                                   n_ice     = qni_in(:,:,:),        &
                                   n_rain    = qnr_in(:,:,:),        &
                                   n_snow    = qns_in(:,:,:),        &
                                   n_graupel = qng_in(:,:,:),        &
                                   n_hail    = qnh_in(:,:,:),        &
                                   z_radar   = z_radar(:,:,:)        )
      
    CASE ( 7 )

      ! Note: computations include halos and boundaries for consistency to EMVORADO internal computations
      CALL compute_field_dbz_2mom( npr       = ni,                   &
                                   nlev      = nj,                   &
                                   nblks     = nk,                   &
                                   startblk  = 1,                    &
                                   endblk    = nk,                   &
                                   jk_start  = 1,                    &
                                   startidx1 = 1,                    &
                                   endidx2   = ni,                   &
                                   lmessage_light = (ldebug_radsim .AND. my_radar_id == 0), &
                                   lmessage_full  = (ldebug_radsim), &
                                   my_id_for_message = my_radar_id,  &
                                   rho_w     = rho_water,            &
                                   rho_ice   = rho_ice,              &
                                   K_w       = K_water,              &
                                   K_ice     = K_ice,                &
                                   T_melt    = t0melt,               &
                                   q_crit_radar = 1e-8_wp,           &
                                   T         = t(:,:,:),             &
                                   rho       = rho(:,:,:),           &
                                   q_cloud   = qc_in(:,:,:),         &
                                   q_ice     = qi_in(:,:,:),         &
                                   q_rain    = qr_in(:,:,:),         &
                                   q_snow    = qs_in(:,:,:),         &
                                   q_graupel = qg_in(:,:,:),         &
                                   q_hail    = qh_in(:,:,:),         &
                                   n_cloud   = qnc_in(:,:,:),        &
                                   n_ice     = qni_in(:,:,:),        &
                                   n_rain    = qnr_in(:,:,:),        &
                                   n_snow    = qns_in(:,:,:),        &
                                   n_graupel = qng_in(:,:,:),        &
                                   n_hail    = qnh_in(:,:,:),        &
                                   ql_graupel= qgl_in(:,:,:),        &
                                   ql_hail   = qhl_in(:,:,:),        &
                                   z_radar   = z_radar(:,:,:)        )
      
    CASE DEFAULT
      
      yerrmsg(:) = ' '
      WRITE (yerrmsg,'(a,i3)') 'Error itype_gscp_model_in: scheme not implemented in EMVORADO: itype_gscp_model_in = ', &
           itype_gscp_model_in
      CALL abort_run(my_radar_id, 35621, yerrmsg, yzroutine)

    END SELECT

    

  END SUBROUTINE get_dbz3dlin_with_model_method_2mom

  !==============================================================================
  !=================================================================================

  !============================================================================

  !============================================================================
  ! 
  ! Setup of the auxiliary rotated lat/lon grid to find the index of
  !  nearest ICON cells for given geographic lat/lon coordinates
  ! 
  !============================================================================
  
  SUBROUTINE setup_auxgrid_for_cellindex (idom, nradsta, &
#ifdef HAVE_RADARFWO
                                          rsm, &
#endif
                                          ldebug)

    ! INPUT variables:
    !-----------------

    INTEGER, INTENT(in)               :: idom, nradsta
#ifdef HAVE_RADARFWO
    TYPE(radar_meta_type), INTENT(in) :: rsm(:)
#endif
    LOGICAL, INTENT(in)               :: ldebug

    ! Local variables:
    !-----------------

    INTEGER, PARAMETER :: oversamp_fac = 7

    LOGICAL           :: not_owned
    INTEGER           :: idom_fwo, i, ii, iidx, k, kk, iblk, i_startblk, i_endblk, is, ie, in, jn, &
         &               inear(nproma,p_patch(idom)%nblks_c), jnear(nproma,p_patch(idom)%nblks_c), &
         &               i_idx, i_blk, ni_over, nj_over, num_neigh, &
         &               global_idx_cell(nproma,p_patch(idom)%nblks_c)

    REAL(KIND=dp)     :: rlon, rlat, clon_min, clon_max, clat_min, clat_max, &
         &               dist, distlat, distlon, tmplon, tmplat, edge_length, &
         &               dlon_edge, dlat_edge, maxx_edge_length

    REAL(KIND=dp), DIMENSION(nproma,p_patch(idom)%nblks_c)       :: &
         &               max_edge_length

    REAL(KIND=dp), DIMENSION(nproma,p_patch(idom)%nblks_c,3)     :: &
         &               lon_vertex, lat_vertex

    REAL(KIND=dp), DIMENSION(nproma,p_patch(idom)%nblks_c)       :: &
         &               clon, clat

    INTEGER, ALLOCATABLE, DIMENSION(:,:) :: counter

    CHARACTER(len=32)  :: yzroutine

    CHARACTER(len=cmaxlen) :: testfile


    ! Code:
    !------

    yzroutine(:) = ' '
    yzroutine    = 'setup_auxgrid_for_cellindex'

    IF (ldebug .OR. my_radar_id == 0) WRITE(*,*) TRIM(yzroutine)//' on proc ', my_radar_id


    ! Internal domain index in EMVORADO of the idom_model'th domain:
#ifdef HAVE_RADARFWO
    idom_fwo = list_domains_for_radar(idom)
#endif

    !====================================================================================
    !
    ! 1) Find the lokal bounding box of the radar-covered region in terms
    !    of a rotated lat/lon grid centered around the center of the radar-covered region.
    !    First, compute the global bounding box, then the lokal box for the specific
    !    processor domain:

    ! 1a) Find the center of the region:
    !     ------------------------------
    clon_min = HUGE(1.0_dp)
    clon_max = -HUGE(1.0_dp)
    clat_min = HUGE(1.0_dp)
    clat_max = -HUGE(1.0_dp)
#ifdef HAVE_RADARFWO
    DO i=1, nradsta
      clon_min = MIN(clon_min, rsm(i)%lon)
      clon_max = MAX(clon_max, rsm(i)%lon)
      clat_min = MIN(clat_min, rsm(i)%lat)
      clat_max = MAX(clat_max, rsm(i)%lat)
    END DO
#endif

    cindex_grid(idom_fwo)%r_earth  = r_earth_dp

    ! Radars spread across the dateline? We assume radar station spreading over areas smaller than half of the world.
    ! Center point of all the radar stations should be the reference point for the auxiliary grid:
    IF (clon_max-clon_min > 180.0) THEN
      cindex_grid(idom_fwo)%pollon = (clon_min+360.0_dp + clon_max) * 0.5_dp
     ELSE
      cindex_grid(idom_fwo)%pollon = (clon_max          + clon_min) * 0.5_dp
    END IF
    cindex_grid(idom_fwo)%pollat = 0.5_dp * (clat_max + clat_min)

    IF (cindex_grid(idom_fwo)%pollat >= 0.0) THEN
      ! Center is on the Northern hemisphere, so rotated pole is computed in the following way:
      cindex_grid(idom_fwo)%pollon = MODULO(cindex_grid(idom_fwo)%pollon, 360.0_dp) - 180.0_dp
      cindex_grid(idom_fwo)%pollat = 90.0_dp - cindex_grid(idom_fwo)%pollat
      cindex_grid(idom_fwo)%polgam = 0.0_dp
    ELSE
      ! Center is on the Southern hemisphere, so rotated pole is computed in the following way:
      cindex_grid(idom_fwo)%pollon = MODULO(cindex_grid(idom_fwo)%pollon, 360.0_dp)
      cindex_grid(idom_fwo)%pollat = -cindex_grid(idom_fwo)%pollat
      cindex_grid(idom_fwo)%polgam = 0.0_dp
    END IF

    ! 1b) Find the total bounding box of the radar-covered region:
    !     ---------------------------------------------------------
    clon_min = HUGE(1.0_dp)
    clon_max = -HUGE(1.0_dp)
    clat_min = HUGE(1.0_dp)
    clat_max = -HUGE(1.0_dp)
    DO i=1, nradsta

#ifdef HAVE_RADARFWO
      ! Radar station coordinates in rotated grid:
      CALL geo2rotll_coord (rsm(i)%lon, rsm(i)%lat, &
                            cindex_grid(idom_fwo)%pollon, cindex_grid(idom_fwo)%pollat, cindex_grid(idom_fwo)%polgam, &
                            rlon, rlat)

      ! Radar range in terms of longitudial / latitudinal differences:
      distlat = 1.1_dp * rsm(i)%ra_inc*rsm(i)%nra / cindex_grid(idom_fwo)%r_earth * raddeg
      distlon = 1.1_dp * rsm(i)%ra_inc*rsm(i)%nra / ( cindex_grid(idom_fwo)%r_earth * COS((ABS(rlat)+distlat)*degrad) ) * raddeg
#endif

      clon_min = MIN(clon_min, rlon - distlon)
      clon_max = MAX(clon_max, rlon + distlon)
      clat_min = MIN(clat_min, rlat - distlat)
      clat_max = MAX(clat_max, rlat + distlat)

    END DO

    IF (ldebug .AND. my_cart_id_fwo == 0) THEN
      WRITE (*,'(a,4(1x,f12.5))') 'INFO '//TRIM(yzroutine)//': Total (rotated lat/lon) bounding box of radar-covered area ', &
           clon_min, clon_max, clat_min, clat_max
    END IF

    cindex_grid(idom_fwo)%startlon_tot = clon_min
    cindex_grid(idom_fwo)%startlat_tot = clat_min
#ifdef HAVE_RADARFWO
    cindex_grid(idom_fwo)%dlon = grid_length_model(idom) / (cindex_grid(idom_fwo)%r_earth * oversamp_fac) * raddeg
#endif
    cindex_grid(idom_fwo)%dlat = cindex_grid(idom_fwo)%dlon
    cindex_grid(idom_fwo)%nlon_tot = FLOOR((clon_max-clon_min)/cindex_grid(idom_fwo)%dlon) + 1
    cindex_grid(idom_fwo)%nlat_tot = FLOOR((clat_max-clat_min)/cindex_grid(idom_fwo)%dlat) + 1

    IF (cindex_grid(idom_fwo)%nlon_tot <= 1 .OR. cindex_grid(idom_fwo)%nlat_tot <= 1) THEN
      WRITE (*,'(a,i3)') 'WARNING '//TRIM(yzroutine)//': Radar-covered area smaller than a grid triangle in domain ', idom
    END IF

    ! 1b) Store and exchange the grid cell center lon/lats:
    !     -------------------------------------------------

    i_startblk = p_patch(idom) % cells % start_block(grf_bdywidth_c+1)
    i_endblk   = p_patch(idom) % cells % end_block(min_rlcell_int)       ! excluding halo cells

    clon(:,:) = -999.99_dp
    clat(:,:) = -999.99_dp
!$OMP PARALLEL PRIVATE(i,k,is,ie)
!$OMP DO
    DO k = i_startblk, i_endblk

      CALL get_indices_c(p_patch(idom), k, i_startblk, i_endblk, is, ie, grf_bdywidth_c+1, min_rlcell_int)

      DO i = is, ie

        ! local grid cell lon/lat:
        clon(i,k) = p_patch(idom) % cells % center(i,k) % lon * raddeg
        clat(i,k) = p_patch(idom) % cells % center(i,k) % lat * raddeg

      END DO
    END DO
!$OMP END DO
!$OMP END PARALLEL

    IF (num_compute_fwo > 1) THEN
      CALL exchange_data(p_patch(idom)%comm_pat_c, clon)
      CALL exchange_data(p_patch(idom)%comm_pat_c, clat)
    END IF


    ! 1c) Find the local bounding box on the PE subdomain including the halo cells for later double-assignment search:
    !     ------------------------------------------------------------------------------------------------------------

    i_startblk = p_patch(idom) % cells % start_block(grf_bdywidth_c+1)
    i_endblk   = p_patch(idom) % cells % end_block(min_rlcell)       ! including halo cells

    clon_min = HUGE(1.0_dp)   
    clon_max = -HUGE(1.0_dp)
    clat_min = HUGE(1.0_dp)
    clat_max = -HUGE(1.0_dp)
!$OMP PARALLEL PRIVATE(dlon_edge,dlat_edge)
    dlon_edge = -HUGE(1.0_dp)
#ifdef HAVE_RADARFWO
    dlat_edge = 2.0_dp * grid_length_model(idom) / ( cindex_grid(idom_fwo)%r_earth) * raddeg
#endif
!$OMP DO PRIVATE(i,k,is,ie,rlon,rlat) REDUCTION(min:clon_min,clat_min), &
!$OMP&   REDUCTION(max:clon_max,clat_max)
    DO k = i_startblk, i_endblk

      CALL get_indices_c(p_patch(idom), k, i_startblk, i_endblk, is, ie, grf_bdywidth_c+1, min_rlcell)

      DO i = is, ie

        ! The next IF is needed because not all halo cells are exchanged by exchange_data() above.
        ! I think that halo cells on the outer boundary are part of the loop, but not part of the inner
        !  boundaries which have been exchanged.
        IF (clon(i,k) > -900.0_dp) THEN

          ! ... converted to rotated grid rlon/rlat:
          CALL geo2rotll_coord (clon(i,k), clat(i,k), &
                                cindex_grid(idom_fwo)%pollon, cindex_grid(idom_fwo)%pollat, cindex_grid(idom_fwo)%polgam, &
                                rlon, rlat)

          ! ... min/max bounding box + some frame representing approx. 2 triangles for safety:
#ifdef HAVE_RADARFWO
          dlon_edge = 2.0_dp * grid_length_model(idom) / ( cindex_grid(idom_fwo)%r_earth * COS(rlat*degrad) ) * raddeg
#endif
          clon_min = MIN(clon_min, rlon - dlon_edge)
          clon_max = MAX(clon_max, rlon + dlon_edge)
          clat_min = MIN(clat_min, rlat - dlat_edge)
          clat_max = MAX(clat_max, rlat + dlat_edge)

        END IF

      END DO
    END DO
!$OMP END DO
!$OMP END PARALLEL
    
    IF (clon_min < cindex_grid(idom_fwo)%startlon_tot+(cindex_grid(idom_fwo)%nlon_tot-1)*cindex_grid(idom_fwo)%dlon .OR. &
        clon_max > cindex_grid(idom_fwo)%startlon_tot ) THEN
      ! lower left corner or upper right corner inside of the radar-covered area:
      cindex_grid(idom_fwo)%startlon = cindex_grid(idom_fwo)%startlon_tot + &
           FLOOR((clon_min-cindex_grid(idom_fwo)%startlon_tot)/cindex_grid(idom_fwo)%dlon) * cindex_grid(idom_fwo)%dlon
      cindex_grid(idom_fwo)%nlon = CEILING((clon_max-cindex_grid(idom_fwo)%startlon) / cindex_grid(idom_fwo)%dlon)
    ELSE
      cindex_grid(idom_fwo)%nlon = 0
    END IF

    IF (clat_min < cindex_grid(idom_fwo)%startlat_tot+(cindex_grid(idom_fwo)%nlat_tot-1)*cindex_grid(idom_fwo)%dlat .OR. &
        clat_max > cindex_grid(idom_fwo)%startlat_tot ) THEN
      ! lower left corner or upper right corner inside of the radar-covered area:
      cindex_grid(idom_fwo)%startlat = cindex_grid(idom_fwo)%startlat_tot + &
           FLOOR((clat_min-cindex_grid(idom_fwo)%startlat_tot)/cindex_grid(idom_fwo)%dlat) * cindex_grid(idom_fwo)%dlat
      cindex_grid(idom_fwo)%nlat = CEILING((clat_max-cindex_grid(idom_fwo)%startlat) / cindex_grid(idom_fwo)%dlat)
    ELSE
      cindex_grid(idom_fwo)%nlat = 0
    END IF

    IF (ldebug) THEN
      IF (cindex_grid(idom_fwo)%nlon == 0 .OR. cindex_grid(idom_fwo)%nlat == 0) THEN
        WRITE (*,'(a,i5)') 'INFO '//TRIM(yzroutine)//': No radar-area on work PE ', my_cart_id_fwo
      ELSE
        WRITE (*,'(a,i5,2(/,i5,1x,a,f12.5,T36,a,i5,T53,a,f12.7,T79,a,f12.5))') 'INFO '//TRIM(yzroutine)//&
             ': Bounding box of radar-covered area on work PE ', my_cart_id_fwo, &
             my_cart_id_fwo, 'startlon = ',cindex_grid(idom_fwo)%startlon, &
             'nlon = ',cindex_grid(idom_fwo)%nlon, &
             'dlon = ',cindex_grid(idom_fwo)%dlon, &
             'endlon = ',cindex_grid(idom_fwo)%startlon + (cindex_grid(idom_fwo)%nlon-1)*cindex_grid(idom_fwo)%dlon, &
             my_cart_id_fwo, 'startlat = ',cindex_grid(idom_fwo)%startlat, &
             'nlat = ',cindex_grid(idom_fwo)%nlat, &
             'dlat = ',cindex_grid(idom_fwo)%dlat, &
             'endlat = ',cindex_grid(idom_fwo)%startlat + (cindex_grid(idom_fwo)%nlat-1)*cindex_grid(idom_fwo)%dlat
      END IF
    END IF


    ! 1d) Store and exchange the global cell indices for later use:
    !     ---------------------------------------------------------

    i_startblk = p_patch(idom) % cells % start_block(grf_bdywidth_c+1)
    i_endblk   = p_patch(idom) % cells % end_block(min_rlcell_int)  ! excluding halo cells

    global_idx_cell(:,:) = -HUGE(1) 
!$OMP PARALLEL
!$OMP DO PRIVATE(i,k,is,ie)
    DO k = i_startblk, i_endblk

      CALL get_indices_c(p_patch(idom), k, i_startblk, i_endblk, is, ie, grf_bdywidth_c+1, min_rlcell_int)

      DO i = is, ie

        ! Find and store global index of cell, also for later re-use:
        ii = idx_1d(i,k)
        global_idx_cell(i,k) = p_patch(idom) % cells % decomp_info % glb_index(ii)

      END DO
    END DO
!$OMP END DO
!$OMP END PARALLEL

    IF (num_compute_fwo > 1) THEN
      CALL exchange_data(p_patch(idom)%comm_pat_c, global_idx_cell)
    END IF


    ! 1e) Store and exchange the cell vertex coordinates and max(edge_length) for later use:
    !     ----------------------------------------------------------------------------------

    i_startblk = p_patch(idom) % cells % start_block(grf_bdywidth_c+1)
    i_endblk   = p_patch(idom) % cells % end_block(min_rlcell_int)  ! exluding the halo cells

    max_edge_length(:,:) = -HUGE(1.0_dp)
    lon_vertex(:,:,:)    = -999.99_dp
    lat_vertex(:,:,:)    = -999.99_dp
!$OMP PARALLEL
!$OMP DO PRIVATE(i,k,is,ie,i_idx,i_blk,edge_length,rlon,rlat)
    DO k = i_startblk, i_endblk

      CALL get_indices_c(p_patch(idom), k, i_startblk, i_endblk, is, ie, grf_bdywidth_c+1, min_rlcell_int)

      DO ii = 1, 3
        DO i = is, ie

          i_idx       = p_patch(idom) % cells % edge_idx(i,k,ii)
          i_blk       = p_patch(idom) % cells % edge_blk(i,k,ii)
          edge_length = p_patch(idom) % edges % primal_edge_length(i_idx,i_blk)

          max_edge_length(i,k) = MAX(max_edge_length(i,k), edge_length)

          i_idx              = p_patch(idom) % cells % vertex_idx(i,k,ii)
          i_blk              = p_patch(idom) % cells % vertex_blk(i,k,ii)
          lon_vertex(i,k,ii) = p_patch(idom) % verts % vertex(i_idx,i_blk) % lon * raddeg
          lat_vertex(i,k,ii) = p_patch(idom) % verts % vertex(i_idx,i_blk) % lat * raddeg

          ! ... converted to rotated grid rlon/rlat:
          CALL geo2rotll_coord (lon_vertex(i,k,ii), lat_vertex(i,k,ii), &
               cindex_grid(idom_fwo)%pollon, cindex_grid(idom_fwo)%pollat, cindex_grid(idom_fwo)%polgam, &
               rlon, rlat)

          lon_vertex(i,k,ii) = rlon
          lat_vertex(i,k,ii) = rlat

        END DO
      END DO
        
    END DO
!$OMP END DO
!$OMP END PARALLEL
   
    IF (num_compute_fwo > 1) THEN
      CALL exchange_data(p_patch(idom)%comm_pat_c, max_edge_length(:,:))
      DO ii=1, 3
        CALL exchange_data(p_patch(idom)%comm_pat_c, lon_vertex(:,:,ii))
        CALL exchange_data(p_patch(idom)%comm_pat_c, lat_vertex(:,:,ii))
      END DO
    END IF


    IF (cindex_grid(idom_fwo)%nlon > 0 .AND. cindex_grid(idom_fwo)%nlat > 0) THEN

      !====================================================================================
      !
      ! 2) Allocate arrays for storing the index and the distance of the triangle
      !    with the nearest ICON cell center to each grid point of the rotated lat/lon aux.
      !    grid:
      !

      ALLOCATE( cindex_grid(idom_fwo)%cind ( cindex_grid(idom_fwo)%nlon, cindex_grid(idom_fwo)%nlat ) )
      
      ALLOCATE( cindex_grid(idom_fwo)%cind_glob ( cindex_grid(idom_fwo)%nlon, cindex_grid(idom_fwo)%nlat ) )
      
      ALLOCATE( cindex_grid(idom_fwo)%idx  ( cindex_grid(idom_fwo)%nlon, cindex_grid(idom_fwo)%nlat ) )

      ALLOCATE( cindex_grid(idom_fwo)%blk  ( cindex_grid(idom_fwo)%nlon, cindex_grid(idom_fwo)%nlat ) )

      ALLOCATE( cindex_grid(idom_fwo)%dist ( cindex_grid(idom_fwo)%nlon, cindex_grid(idom_fwo)%nlat ) )

      ALLOCATE( counter ( cindex_grid(idom_fwo)%nlon, cindex_grid(idom_fwo)%nlat ) )

      CALL init_vari(cindex_grid(idom_fwo)%cind(:,:)      , -HUGE(1) )
      CALL init_vari(cindex_grid(idom_fwo)%cind_glob(:,:) , -HUGE(1) )
      CALL init_vari(cindex_grid(idom_fwo)%idx(:,:)       , -HUGE(1) )
      CALL init_vari(cindex_grid(idom_fwo)%blk(:,:)       , -HUGE(1) )
      CALL init_vari(cindex_grid(idom_fwo)%dist(:,:)      ,  HUGE(1.0_dp) )
      CALL init_vari(counter(:,:)                         ,  0 )

      !====================================================================================
      !
      ! 3) Determine and store index and distance of nearest ICON cell for interior cells.
      !    In case of double assignments (lon/lat point is exactly on an edge), take
      !     the cell with the larger global grid index. This should lead to reproducible
      !     results.
      !

      ! 3a) Find nearest points to the ICON interior cells + halo cells:
      !     ------------------------------------------------------------

      i_startblk = p_patch(idom) % cells % start_block(grf_bdywidth_c+1)
      i_endblk   = p_patch(idom) % cells % end_block(min_rlcell)  ! including halo cells

      inear(:,:) = -HUGE(1)
      jnear(:,:) = -HUGE(1)
!$OMP PARALLEL
!$OMP DO PRIVATE(i,k,is,ie,rlon,rlat,tmplon,tmplat)
      DO k = i_startblk, i_endblk

        CALL get_indices_c(p_patch(idom), k, i_startblk, i_endblk, is, ie, grf_bdywidth_c+1, min_rlcell)
!NEC$ ivdep
        DO i = is, ie

          ! The next IF is needed because not all halo cells are exchanged by exchange_data() above.
          ! I think that halo cells on the outer boundary are part of the loop, but not part of the inner
          !  boundaries which have been exchanged.
          IF (clon(i,k) > -900.0_dp) THEN

            ! ... converted to rotated grid rlon/rlat:
            CALL geo2rotll_coord (clon(i,k), clat(i,k), &
                 cindex_grid(idom_fwo)%pollon, cindex_grid(idom_fwo)%pollat, cindex_grid(idom_fwo)%polgam, &
                 rlon, rlat)

            inear(i,k) = NINT( ( rlon - cindex_grid(idom_fwo)%startlon ) / cindex_grid(idom_fwo)%dlon ) + 1
            jnear(i,k) = NINT( ( rlat - cindex_grid(idom_fwo)%startlat ) / cindex_grid(idom_fwo)%dlat ) + 1

            ! Compute cindex_grid(idom_fwo)%cind(inear,jnear):
            cindex_grid(idom_fwo)%idx(inear(i,k),jnear(i,k)) = i
            cindex_grid(idom_fwo)%blk(inear(i,k),jnear(i,k)) = k
            CALL sub2ind2D(i, k, nproma, cindex_grid(idom_fwo)%cind(inear(i,k),jnear(i,k)) )
            counter(inear(i,k),jnear(i,k))   = counter(inear(i,k),jnear(i,k)) + 1
            
            ! Compute cindex_grid(idom_fwo)%dist(inear,jnear):
            tmplon = cindex_grid(idom_fwo)%startlon + (inear(i,k)-1)*cindex_grid(idom_fwo)%dlon
            tmplat = cindex_grid(idom_fwo)%startlat + (jnear(i,k)-1)*cindex_grid(idom_fwo)%dlat
#ifdef HAVE_RADARFWO
            cindex_grid(idom_fwo)%dist(inear(i,k),jnear(i,k)) = geo_dist(rlon, rlat, tmplon, tmplat, &
                                                                         cindex_grid(idom_fwo)%r_earth, 0.0_dp)
#endif
            ! Store global index of nearest neighbour cell:
            cindex_grid(idom_fwo)%cind_glob(inear(i,k),jnear(i,k)) = global_idx_cell(i,k)

          END IF

        END DO
      END DO
!$OMP END DO
!$OMP END PARALLEL
      
    ELSE
      
      inear(:,:) = -HUGE(1)
      jnear(:,:) = -HUGE(1)

    END IF

    ! 3b) Fill space between nearest points to ICON interior+halo cells:
    !     --------------------------------------------------------------

    IF (cindex_grid(idom_fwo)%nlon > 0 .AND. cindex_grid(idom_fwo)%nlat > 0) THEN

      i_startblk = p_patch(idom) % cells % start_block(grf_bdywidth_c+1)
      i_endblk   = p_patch(idom) % cells % end_block(min_rlcell)    ! including halo cells

!$OMP PARALLEL
!$OMP DO PRIVATE(i,k,is,ie,rlon,rlat,ii,kk,ni_over,nj_over,in,jn,tmplon,tmplat,dist)
      DO k = i_startblk, i_endblk

        CALL get_indices_c(p_patch(idom), k, i_startblk, i_endblk, is, ie, grf_bdywidth_c+1, min_rlcell)
        
        DO i = is, ie

          ! The next IF is needed because not all halo cells are exchanged by exchange_data() above.
          ! I think that halo cells on the outer boundary are part of the loop, but not part of the inner
          !  boundaries which have been exchanged.
          IF (inear(i,k) > -HUGE(1)) THEN

            ii   = cindex_grid(idom_fwo)%cind(inear(i,k),jnear(i,k))

            rlon = cindex_grid(idom_fwo)%startlon + (inear(i,k)-1) * cindex_grid(idom_fwo)%dlon
            rlat = cindex_grid(idom_fwo)%startlat + (jnear(i,k)-1) * cindex_grid(idom_fwo)%dlat

            ! determine search region in terms of auxiliary oversampling grid points:
            ni_over = FLOOR( max_edge_length(i,k) / ( cindex_grid(idom_fwo)%r_earth*cindex_grid(idom_fwo)%dlon*degrad * &
                 COS( (rlat+0.5_dp*raddeg*max_edge_length(i,k)/cindex_grid(idom_fwo)%r_earth) * degrad ) ) ) + 5
            ni_over = ni_over + MOD(ni_over+1, 2)  ! makes ni_over an odd integer
            nj_over = FLOOR( max_edge_length(i,k) / ( cindex_grid(idom_fwo)%r_earth*cindex_grid(idom_fwo)%dlat*degrad ) ) + 5
            nj_over = nj_over + MOD(nj_over+1, 2)  ! makes nj_over an odd integer

            DO jn = jnear(i,k)-nj_over/2 , jnear(i,k)+nj_over/2
              DO in = inear(i,k)-ni_over/2 , inear(i,k)+ni_over/2
            
                IF ( ( in /= inear(i,k) .AND. in >= 1 .AND. in <= cindex_grid(idom_fwo)%nlon ) .OR. &
                     ( jn /= jnear(i,k) .AND. jn >= 1 .AND. jn <= cindex_grid(idom_fwo)%nlat ) ) THEN
              
                  tmplon = cindex_grid(idom_fwo)%startlon + (in-1) * cindex_grid(idom_fwo)%dlon
                  tmplat = cindex_grid(idom_fwo)%startlat + (jn-1) * cindex_grid(idom_fwo)%dlat
#ifdef HAVE_RADARFWO
                  dist = geo_dist(rlon, rlat, tmplon, tmplat, cindex_grid(idom_fwo)%r_earth, 0.0_dp)
                  
                  ! Check if auxiliary point is inside the triangle. If yes, check by inspection
                  !  of the stored global nearest cell index, if it has already
                  !  been assigned to another cell. If yes, the cell with the larger global
                  !  index shall win:
                  IF ( is_inside_triangle(tmplon, tmplat, lon_vertex(i,k,1:3), lat_vertex(i,k,1:3)) ) THEN
                    kk = cindex_grid(idom_fwo)%cind_glob(in,jn)
                    IF ( global_idx_cell(i,k) > kk ) THEN
                      cindex_grid(idom_fwo)%cind(in,jn)      = ii
                      cindex_grid(idom_fwo)%idx (in,jn)      = i
                      cindex_grid(idom_fwo)%blk (in,jn)      = k
                      cindex_grid(idom_fwo)%dist(in,jn)      = dist
                      cindex_grid(idom_fwo)%cind_glob(in,jn) = global_idx_cell(i,k)
                      counter(in,jn)                         = counter(in,jn) + 1
                    END IF
                  END IF
#endif
                
                END IF
              
              END DO
            END DO
          
          END IF

        END DO
      END DO
!$OMP END DO
!$OMP END PARALLEL
      
    END IF


    IF (cindex_grid(idom_fwo)%nlon > 0 .AND. cindex_grid(idom_fwo)%nlat > 0 .AND. ldebug) THEN

!===============================================================================
!    Debug file output (commented out, but left in the code for eventual later use):
!
      testfile(:) = ' '
      WRITE (testfile, '("testdomain_",i4.4,".dat")') my_cart_id_fwo
      OPEN(350+my_cart_id_fwo, file=TRIM(testfile), status='replace', form='formatted')

      testfile(:) = ' '
      WRITE (testfile, '(i4)') MIN(cindex_grid(idom_fwo)%nlon, 200)
      testfile = ADJUSTL(testfile)
      DO k = 1, cindex_grid(idom_fwo)%nlat
        WRITE (350+my_cart_id_fwo, '('//TRIM(testfile)//'(i6))') &
             (cindex_grid(idom_fwo)%cind(i,k), i=1,MIN(cindex_grid(idom_fwo)%nlon, 200))
      END DO
      CLOSE(350+my_cart_id_fwo)


      testfile(:) = ' '
      WRITE (testfile, '("testdist_",i4.4,".dat")') my_cart_id_fwo
      OPEN(350+my_cart_id_fwo, file=TRIM(testfile), status='replace', form='formatted')
      
      testfile(:) = ' '
      WRITE (testfile, '(i4)') MIN(cindex_grid(idom_fwo)%nlon, 200)
      testfile = ADJUSTL(testfile)
      DO k = 1, cindex_grid(idom_fwo)%nlat
        WRITE (350+my_cart_id_fwo, '('//TRIM(testfile)//'(f8.0))') &
             (cindex_grid(idom_fwo)%dist(i,k)*1e-3_dp, i=1,MIN(cindex_grid(idom_fwo)%nlon,200))
      END DO
      CLOSE(350+my_cart_id_fwo)


      testfile(:) = ' '
      WRITE (testfile, '("testcounter_",i4.4,".dat")') my_cart_id_fwo
      OPEN(350+my_cart_id_fwo, file=TRIM(testfile), status='replace', form='formatted')
      
      testfile(:) = ' '
      WRITE (testfile, '(i4)') MIN(cindex_grid(idom_fwo)%nlon, HUGE(1))
      testfile = ADJUSTL(testfile)
      DO k = 1, cindex_grid(idom_fwo)%nlat
        WRITE (350+my_cart_id_fwo, '('//TRIM(testfile)//'(i3))') &
             (counter(i,k), i=1,MIN(cindex_grid(idom_fwo)%nlon,HUGE(1)))
      END DO
      CLOSE(350+my_cart_id_fwo)
!
! END debug output
!===============================================================================

    
      ! 3d) Search for multiple assignments of ICON cells to aux. grid points.
      !     If such points are found a warning is issued:
      !     ------------------------------------------------------------------

!$OMP PARALLEL
!$OMP DO PRIVATE(i,k)
      DO k = 2, cindex_grid(idom_fwo)%nlat - 1
        DO i = 2, cindex_grid(idom_fwo)%nlon - 1
        
          IF (counter(i,k) > 1) THEN
            ! This should not have happened:
            WRITE (*,'(a,i5,a,i5,a,i5,a,i5,a)') 'WARNING '//TRIM(yzroutine)//' work PE ', my_cart_id_fwo, &
                 ' : aux. rot-lat/lon point (', i, ',', k, &
                 ') had more than one (n=',counter(i,k),') assignments to nearest neighbouring ICON cells!'
          END IF
          
        END DO
      END DO
!$OMP END DO
!$OMP END PARALLEL

      ! 3e) Search for "holes" in the aux. grid, which are cells with no assignment and
      !     surrounded by at least 7 assigned neighbours. If such holes are found a warning
      !     is issued:
      !     -------------------------------------------------------------------------------

!$OMP PARALLEL
!$OMP DO PRIVATE(i,k,num_neigh,tmplon,tmplat,in,jn,iidx,iblk)
      DO k = 2, cindex_grid(idom_fwo)%nlat - 1
        DO i = 2, cindex_grid(idom_fwo)%nlon - 1
        
          num_neigh = SUM( MAX( MIN( counter(i-1:i+1,k-1:k+1), 1), 0) )

          IF (counter(i,k) == 0 .AND. num_neigh >= 7) THEN
            ! This is a hole:
            WRITE (*,'(a,i5,a,i5,a,i5,a)') 'WARNING '//TRIM(yzroutine)//' work PE ', my_cart_id_fwo, &
                 ' : aux. rot-lat/lon point (', i, ',', k, &
                 ') is an isolated hole in the auxiliary grid! No nearest ICON cell has been assigned!'
            
            tmplon = cindex_grid(idom_fwo)%startlon + (i-1) * cindex_grid(idom_fwo)%dlon
            tmplat = cindex_grid(idom_fwo)%startlat + (k-1) * cindex_grid(idom_fwo)%dlat
            DO jn=k-1,k+1
              DO in=i-1,i+1
                IF (in /= i .OR. jn /= k) THEN
                  iidx = cindex_grid(idom_fwo)%idx(in,jn)
                  iblk = cindex_grid(idom_fwo)%blk(in,jn)
                  IF (iidx > -HUGE(1) .AND. iblk > -HUGE(1)) THEN
                    PRINT '(a,2i5,4(1x,es21.14),"   ",4(1x,es21.14))', '  DEBUG  ', &
                         i, k, tmplon, lon_vertex(iidx,iblk,1:3), tmplat, lat_vertex(iidx,iblk,1:3)
                  ELSE
                    PRINT '(a,2i5,4(1x,es21.14),"   ",4(1x,es21.14))', '  DEBUG  ', &
                         i, k, tmplon, -999.99_dp, -999.99_dp, -999.99_dp, tmplat, -999.99_dp, -999.99_dp, -999.99_dp
                  END IF
                END IF
              END DO
            END DO
            
          END IF
          
        END DO
      END DO
!$OMP END DO
!$OMP END PARALLEL

      ! Clean up memory
      DEALLOCATE (counter)

    END IF  ! PE domain overlaps with radar-covered area

    IF (ldebug) WRITE(*,*) 'Done with '//TRIM(yzroutine)//' on proc ', my_radar_id

  END SUBROUTINE setup_auxgrid_for_cellindex

  !------------------------------------------------------------------------------

  FUNCTION is_inside_triangle(lon, lat, lon_verts, lat_verts) RESULT (inside)

    REAL(kind=dp), INTENT(in) :: lon, lat, lon_verts(3), lat_verts(3)
    LOGICAL                   :: inside, is_on_left(3)
    REAL(kind=dp)             :: xr(3), yr(3), xt(3), yt(3)
    
    ! direction vectors of the 3 lines:
    xr(1) = lon_verts(2)-lon_verts(1)
    xr(2) = lon_verts(3)-lon_verts(2)
    xr(3) = lon_verts(1)-lon_verts(3)
    yr(1) = lat_verts(2)-lat_verts(1)
    yr(2) = lat_verts(3)-lat_verts(2)
    yr(3) = lat_verts(1)-lat_verts(3)

    ! translated coordinates of point (lon,lat), so that beginnings of the 3 lines are the new coordinate origins:
    xt(1) = lon-lon_verts(1)
    xt(2) = lon-lon_verts(2)
    xt(3) = lon-lon_verts(3)
    yt(1) = lat-lat_verts(1)
    yt(2) = lat-lat_verts(2)
    yt(3) = lat-lat_verts(3)

    !     /xr xt\       means that the volume of the parallelogram spanned by
    ! det |     | >= 0  the two vectors r and t is mathematically positively oriented,
    !     \yr yt/       i.e., the end point of t is to the "left" of r.
    !                  
    is_on_left(1) = xr(1)*yt(1)-xt(1)*yr(1) >= 0.0_dp
    is_on_left(2) = xr(2)*yt(2)-xt(2)*yr(2) >= 0.0_dp
    is_on_left(3) = xr(3)*yt(3)-xt(3)*yr(3) >= 0.0_dp

    inside = ALL(is_on_left)

  END FUNCTION is_inside_triangle

  !============================================================================

  SUBROUTINE setup_runtime_timings ()

    IMPLICIT NONE

    ! Set up runtime timings (will appear in the COSMO YUTIMING output):
    i_fwo_prep_compute = timer_radar_prep_compute
    i_fwo_bubbles      = timer_radar_bubbles
    i_fwo_composites   = timer_radar_composites
    i_fwo_ini          = timer_radar_ini
    i_fwo_compgrid     = timer_radar_compgrid
    i_fwo_comm         = timer_radar_comm
    i_fwo_ongeom       = timer_radar_ongeom
    i_fwo_comppolar    = timer_radar_comppolar
    i_fwo_out          = timer_radar_out
    i_fwo_barrier      = timer_radar_barrier  

  END SUBROUTINE setup_runtime_timings

  ! call get_runtime_timings ( iflag0, lstart=.TRUE.) at the start of organize_radar() to start timer iflag0
  ! call get_runtime_timings ( iflag1...N ) for all calls except the last to start timers iflag1...N
  ! call get_runtime_timings ( iflagN, lstop=.TRUE.) last timer call repeated for iflag
  !  with the stop signal, to finalize the timers

  SUBROUTINE get_runtime_timings ( iflag, lstart, lstop )

    IMPLICIT NONE

    INTEGER, INTENT(in)           :: iflag
    LOGICAL, INTENT(in), OPTIONAL :: lstart, lstop

    LOGICAL :: zlstart, zlstop
    INTEGER :: lasttimer = -999

    IF (ltimer) THEN

      IF (PRESENT(lstart)) THEN
        zlstart = lstart
      ELSE
        zlstart = .FALSE.
      END IF

      IF (PRESENT(lstop)) THEN
        zlstop = lstop
      ELSE
        zlstop = .FALSE.
      END IF

      IF (lasttimer == -999 .AND. .NOT.zlstart) THEN
        WRITE (*,*) 'ERROR get_runtime_timings() radar: wrong calling sequence, forgot call with lstart=.TRUE.'
        STOP
      END IF

      IF (zlstart) THEN
        CALL timer_start (iflag)
        lasttimer = iflag
      ELSE IF (zlstop) THEN
        CALL timer_stop (lasttimer)
        lasttimer = -999
      ELSE
        CALL timer_stop (lasttimer)
        CALL timer_start(iflag)
        lasttimer = iflag      
      END IF

    END IF

  END SUBROUTINE get_runtime_timings

!============================================================================
!============================================================================



  
!**************************************************************************************************
!**************************************************************************************************
!
! THE FOLLOWING ROUTINES CAN ONLY BE COMPILED WHEN USING EMVORADO
!  BECAUSE OF DEPENDENCIES TO THE EMVORADO CORE MODULES.
!
!**************************************************************************************************
!**************************************************************************************************

!============================================================================
!============================================================================

#ifdef HAVE_RADARFWO

  !============================================================================
  ! 
  ! Function for getting the lower bound of the typical model grid distance
  !  in units of length (meters).
  !
  !============================================================================

  PURE FUNCTION grid_length_model (idom) RESULT (len)

    IMPLICIT NONE

    INTEGER, INTENT(in) :: idom
    REAL(KIND=dp) :: len

    len = p_patch(idom) % geometry_info % mean_characteristic_length

  END FUNCTION grid_length_model

  PURE FUNCTION grid_length_model_x (idom) RESULT (len)

    IMPLICIT NONE

    INTEGER, INTENT(in) :: idom
    REAL(KIND=dp) :: len

    len = grid_length_model(idom)

  END FUNCTION grid_length_model_x

  PURE FUNCTION grid_length_model_y (idom) RESULT (len)

    IMPLICIT NONE

    INTEGER, INTENT(in) :: idom
    REAL(KIND=dp) :: len

    len = grid_length_model(idom)

  END FUNCTION grid_length_model_y

  !------------------------------------------------------------------------------

!================================================================
!
! Domain start and end indices for computations
!
! This is to have the possibility to do computations only
! at those grid points which are absolutely necessary to
! save time.
!
! Has to be called at first before any computations are done
!  (usually at the beginning of the interface routines
!   radar_sb_mie_vec(), ...)
!
!================================================================

  SUBROUTINE get_loc_domain(idom,zlinc_boundlines)

    IMPLICIT NONE

    INTEGER, INTENT(in) :: idom
    LOGICAL, INTENT(in) :: zlinc_boundlines
    
    IF (zlinc_boundlines) THEN
      ! .. Full domain of this PE, i.e., including boundary lines:
      ilow_modelgrid = 1
      iup_modelgrid  = ie_fwo
      jlow_modelgrid = 1
      jup_modelgrid  = je_fwo
!      klow_modelgrid = 1
!      kup_modelgrid  = ke_fwo
! Use reduced domain in any case for ICON!
      klow_modelgrid = p_patch(idom) % cells % start_block(grf_bdywidth_c+1)
      kup_modelgrid  = p_patch(idom) % cells % end_block(min_rlcell_int)
    ELSE
      ! .. reduced domain of this PE, i.e., excluding interior
      !    boundary lines but including boudary halo at the outer
      !    model domain boundaries:
      ilow_modelgrid = 1
      iup_modelgrid  = ie_fwo   ! the few halo/boundary points in the first and last block are not explicitly excluded here for simplicity
      jlow_modelgrid = 1
      jup_modelgrid  = je_fwo
      klow_modelgrid = p_patch(idom) % cells % start_block(grf_bdywidth_c+1)
      kup_modelgrid  = p_patch(idom) % cells % end_block(min_rlcell_int)
    END IF


  END SUBROUTINE get_loc_domain


  !------------------------------------------------------------------------------

  !============================================================================
  ! 
  ! Function for getting the initial date of the model run (String YYYYMMDDhhmmss).
  !
  !============================================================================

  FUNCTION get_datetime_ini () RESULT (d)

    IMPLICIT NONE

    CHARACTER(len=14) :: d
    CHARACTER(len=32) :: tmp_d

    tmp_d(:) = ' '
    CALL datetimeToString( time_config % tc_exp_startdate, tmp_d )
    d(:) = '0'
    d = tmp_d(1:4)//tmp_d(6:7)//tmp_d(9:10)//tmp_d(12:13)//tmp_d(15:16)//tmp_d(18:19)

  END FUNCTION get_datetime_ini

  !------------------------------------------------------------------------------

  !============================================================================
  ! 
  ! Function for getting the actual date of the model run (String YYYYMMDDhhmmss).
  !
  !============================================================================

  FUNCTION get_datetime_act (l_round_to_minute) RESULT (actdate)

    IMPLICIT NONE

    LOGICAL, OPTIONAL       :: l_round_to_minute

    CHARACTER(len=14)       :: actdate

    INTEGER                 :: nactday, zitype_calendar, ntstep
    CHARACTER (LEN=28)      :: dum1
    REAL (KIND=dp)          :: dum2, dtloc
    CHARACTER(len=14)       :: ydate_ini
    LOGICAL                 :: l_round

    ydate_ini = get_datetime_ini ()

    IF (PRESENT(l_round_to_minute)) THEN
      l_round = l_round_to_minute
    ELSE
      l_round = .FALSE.
    END IF
    IF (l_round) THEN
      dtloc = 60.0_dp
    ELSE
      dtloc = dtime
    END IF
    ntstep = NINT(get_model_time_sec () / dtloc)

    zitype_calendar = 0       !  = 0: gregorian calendar (default)
                              !    (but this needs a bug fix in get_utc_date,
                              !    because up to now we only have julian calendar)
                              !  = 1: every year has 360 days
                              !  = 2: every year has 365 days
    actdate(:) = '0'
    IF (ntstep > 0) THEN
      CALL get_utc_date (ntstep, ydate_ini, dtloc, zitype_calendar, actdate, dum1, nactday, dum2)
    ELSE
      actdate = get_datetime_ini ()
    END IF

  END FUNCTION get_datetime_act

  !------------------------------------------------------------------------------

  !============================================================================
  ! 
  ! Subroutine for getting the model time in seconds since model start.
  !
  !============================================================================

  PURE FUNCTION get_model_time_sec () RESULT(time_mod)

    IMPLICIT NONE

    REAL(KIND=dp) :: time_mod

    time_mod = REAL(time_mod_sec, kind=dp)   ! [s] since tc_exp_startdate

  END FUNCTION get_model_time_sec

  !------------------------------------------------------------------------------

  !============================================================================
  ! 
  ! Function for checking if any of the observation times of a radar
  !  falls within the current model time step. This function is used to trigger
  !  the forward simulation for a single radar station.
  !
  !============================================================================


  FUNCTION it_is_time_for_radar ( obs_times ) RESULT (it_is_time)

    IMPLICIT NONE

    REAL(KIND=dp), INTENT(in) :: obs_times(:)
    LOGICAL                       :: it_is_time

    REAL(KIND=dp) :: time_mod

    time_mod = get_model_time_sec ()  ! [s] since tc_exp_startdate

    it_is_time = ANY( obs_times > -900.0_dp  .AND. &
         ( obs_times - 0.5_dp*dtime < time_mod .AND. time_mod <= obs_times + 0.5_dp*dtime ) )

  END FUNCTION it_is_time_for_radar

  !------------------------------------------------------------------------------

  !============================================================================
  !
  ! Function for checking if a given model time is one of the bubble
  !  checking times for the automatic bubble generator (ldo_bubbles=.TRUE.).
  ! Bubble checkting times are characterized by a regular time interval
  !  dt_bubblecheck [s] and a time offset t_offset_bubblecheck [s] (relative
  !  to model start) for the first checking interval.
  ! If t_offset_bubblecheck is 0.0, the model starttime itself is omitted.
  !
  !============================================================================

  FUNCTION it_is_time_for_bubblecheck (time_mod, dt_bubblecheck, t_offset_bubblecheck) &
       RESULT (it_is_time)

    IMPLICIT NONE

    REAL(kind=dp), INTENT(in) :: time_mod             ! actual time in seconds since model start
    REAL(kind=dp), INTENT(in) :: dt_bubblecheck       ! time interval between checks
    REAL(kind=dp), INTENT(in) :: t_offset_bubblecheck ! time offset for the check intervals relative to model start (=0.0 s)

    LOGICAL                   :: it_is_time

    INTEGER                   :: n
    REAL(kind=dp)             :: tol, checktime

    tol = 0.5_dp * dtime
    checktime = time_mod - t_offset_bubblecheck
    n = NINT(checktime/dt_bubblecheck)

    ! Only times after the model start are possible candidates for bubble checking:
    IF (n > 0 .AND. time_mod > dt_bubblecheck-tol) THEN
      it_is_time = (n*dt_bubblecheck-tol < checktime .AND. checktime <= n*dt_bubblecheck+tol)
    ELSE
      it_is_time = .FALSE.
    END IF

  END FUNCTION it_is_time_for_bubblecheck

  !------------------------------------------------------------------------------

  !============================================================================
  ! 
  ! Function for computing the index of the nearest obs_time to the current
  !  model time. The obs times are given in an input vector. If none of the
  !  obs times falls within the current model time step +/- 0.5*dt, no
  !  index can be found and -999 is returned.
  !
  !============================================================================

  FUNCTION get_obstime_ind_of_currtime ( obs_times ) RESULT (i_time)

    IMPLICIT NONE

    REAL(KIND=dp), INTENT(in) :: obs_times(:)
    INTEGER       :: i_time

    REAL(KIND=dp)             :: time_mod
    INTEGER       :: n

    time_mod = get_model_time_sec ()  ! [s] since tc_exp_startdate

    i_time = -999
    DO n = 1, UBOUND(obs_times,1)
      IF ((obs_times(n) - 0.5_dp*dtime < time_mod .AND. time_mod <= obs_times(n) + 0.5_dp*dtime)) THEN
        i_time = n
        EXIT
      END IF
    END DO
    
  END FUNCTION get_obstime_ind_of_currtime

  !------------------------------------------------------------------------------

  !============================================================================
  ! 
  ! Function for checking if one specific observation time of a radar
  !  falls within the current model time step.
  !
  ! Method: nearest neighbour in time
  !
  !============================================================================

  ELEMENTAL FUNCTION check_if_currtime_is_obstime ( obs_time ) RESULT (it_is_time)

    IMPLICIT NONE

    REAL(KIND=dp), INTENT(in) :: obs_time
    LOGICAL                   :: it_is_time

    REAL(KIND=dp) :: time_mod

    time_mod = get_model_time_sec ()  ! [s] since tc_exp_startdate

    it_is_time = ( obs_time > -900.0_dp .AND. ABS(time_mod - obs_time) <= 0.5_dp*dtime )

  END FUNCTION check_if_currtime_is_obstime

  !------------------------------------------------------------------------------

  FUNCTION get_domain_starttime_in_sec (idom) RESULT (time)

    INTEGER, INTENT(in) :: idom
    REAL(KIND=dp) :: time

    ! Start time of the domain in experiment (in seconds from tc_exp_startdate)
    time = MAX(start_time(idom), 0.0_wp)   ! [s] since tc_exp_startdate

  END FUNCTION get_domain_starttime_in_sec

  !------------------------------------------------------------------------------

  FUNCTION get_domain_runstarttime_in_sec (idom) RESULT (time)

    INTEGER, INTENT(in) :: idom
    REAL(KIND=dp) :: time
    REAL(KIND=dp) :: current_run_start

    ! Start time of the domain in the current run (in seconds from tc_exp_startdate)
    current_run_start = REAL( getElapsedSimTimeInSeconds(time_config%tc_startdate), KIND=dp)
    time = MAX( MAX(start_time(idom), 0.0_wp), current_run_start )   ! [s] since tc_exp_startdate

  END FUNCTION get_domain_runstarttime_in_sec

  !------------------------------------------------------------------------------

  FUNCTION get_domain_endtime_in_sec (idom) RESULT (time)

    INTEGER, INTENT(in) :: idom
    REAL(KIND=dp) :: time
    REAL(KIND=dp) :: current_exp_end

    ! End time of the domain in experiment (in seconds from tc_exp_startdate)
    current_exp_end = REAL( getElapsedSimTimeInSeconds(time_config%tc_exp_stopdate), KIND=dp)
    time = MIN ( end_time(idom), current_exp_end)   ! [s] since tc_exp_startdate

  END FUNCTION get_domain_endtime_in_sec

  !------------------------------------------------------------------------------

  FUNCTION get_domain_runendtime_in_sec (idom) RESULT (time)

    INTEGER, INTENT(in) :: idom
    REAL(KIND=dp) :: time
    REAL(KIND=dp) :: current_run_end

    ! End time of the domain in the current run (in seconds from tc_exp_startdate)
    current_run_end = REAL( getElapsedSimTimeInSeconds(time_config%tc_stopdate), KIND=dp)
    time = MIN( end_time(idom), current_run_end )   ! [s] since tc_exp_startdate

  END FUNCTION get_domain_runendtime_in_sec

  !============================================================================
  ! 
  ! Function for checking if one specific observation time of a radar
  !  falls within the total forecast time span of the model.
  !
  ! Method: nearest neighbour in time
  !
  !============================================================================

  FUNCTION check_obstime_within_fcst_scal ( idom, obs_time ) RESULT (within)

    IMPLICIT NONE

    INTEGER, INTENT(in) :: idom
    REAL(KIND=dp),           INTENT(in) :: obs_time
    LOGICAL                             :: within

    within = ( obs_time >= get_domain_starttime_in_sec(idom)-0.5_dp*dtime .AND. &
         obs_time <= get_domain_endtime_in_sec(idom)+0.5_dp*dtime )

  END FUNCTION check_obstime_within_fcst_scal

  !------------------------------------------------------------------------------

  FUNCTION check_obstime_within_fcst_vec ( idom, obs_time ) RESULT (within)

    IMPLICIT NONE

    INTEGER, INTENT(in) :: idom
    REAL(KIND=dp), INTENT(in)           :: obs_time(:)
    LOGICAL                             :: within(1:size(obs_time))

    within = ( obs_time >= get_domain_starttime_in_sec(idom)-0.5_dp*dtime .AND. &
         obs_time <= get_domain_endtime_in_sec(idom)+0.5_dp*dtime )

  END FUNCTION check_obstime_within_fcst_vec

  !------------------------------------------------------------------------------

  !============================================================================
  ! 
  ! Function for checking if one specific observation time of a radar
  !  falls within the time span of the actual model run (for restart runs
  !  this is less than the forecast time span!).
  !
  ! Method: nearest neighbour in time
  !
  !============================================================================

  FUNCTION check_obstime_within_run_scal ( idom, obs_time ) RESULT (within)

    IMPLICIT NONE

    INTEGER, INTENT(in) :: idom
    REAL(KIND=dp),           INTENT(in) :: obs_time
    LOGICAL                             :: within

    within = ( obs_time >= get_domain_runstarttime_in_sec(idom)-0.5_dp*dtime .AND. &
         obs_time <= get_domain_runendtime_in_sec(idom)+0.5_dp*dtime )

  END FUNCTION check_obstime_within_run_scal

  !------------------------------------------------------------------------------

  FUNCTION check_obstime_within_run_vec ( idom, obs_time ) RESULT (within)

    IMPLICIT NONE

    INTEGER, INTENT(in) :: idom
    REAL(KIND=dp), INTENT(in)           :: obs_time(:)
    LOGICAL                             :: within(1:size(obs_time))

    within = ( obs_time >= get_domain_runstarttime_in_sec(idom)-0.5_dp*dtime .AND. &
         obs_time <= get_domain_runendtime_in_sec(idom)+0.5_dp*dtime )

  END FUNCTION check_obstime_within_run_vec

  !------------------------------------------------------------------------------

  !============================================================================
  ! 
  ! Function for getting the needed tolerance (+/-) to match observation times
  !  against the actual (discrete) model time.
  !
  ! Method: for COSMO, this is just half of the model timestep dt.
  !
  !============================================================================

  FUNCTION get_obs_time_tolerance( idom ) RESULT (tol)

    IMPLICIT NONE

    INTEGER, INTENT(in) :: idom
    REAL(KIND=dp) :: tol

    tol = 0.5_dp * dtime

  END FUNCTION get_obs_time_tolerance

  !============================================================================
  ! 
  ! Function for computing the number of regularly spaced obs times within
  !  the model run for a specific domain, based on the model forecast time range and the
  !  regular observation time interval dt_obs. Up to now, these intervals are assumed
  !  to start at t=0.0, not at the true start time of the domain.
  !
  !============================================================================

  FUNCTION num_regular_obstimes ( idom, dt_obs ) RESULT (n_times)

    IMPLICIT NONE

    INTEGER, INTENT(in) :: idom
    REAL(KIND=dp),           INTENT(in) :: dt_obs
    INTEGER       :: n_times

    REAL(KIND=dp) :: time_sec

    time_sec = get_domain_endtime_in_sec(idom) - get_domain_starttime_in_sec(idom)
    n_times = INT(time_sec/dt_obs) + 1

  END FUNCTION num_regular_obstimes


  !------------------------------------------------------------------------------

  !============================================================================
  ! 
  ! Subroutine for getting the model input directory. This will be the
  ! default input directory of radar operator input observation files, if the namelist
  ! variable "ydirradarin" is not specified.
  !
  !============================================================================

  FUNCTION get_model_inputdir ( ) RESULT (ydirin_model)

    IMPLICIT NONE

    CHARACTER(len=cmaxlen) :: ydirin_model

    ydirin_model(:) = ' '

  END FUNCTION get_model_inputdir
  
  !============================================================================
  ! 
  ! Subroutine for getting the model output directory. This will be the
  ! default output directory of radar operator output files, if the namelist
  ! variable "ydirradarout" is not specified.
  !
  !============================================================================

  FUNCTION get_model_outputdir ( ) RESULT (ydirout_model)

    IMPLICIT NONE

    CHARACTER(len=cmaxlen) :: ydirout_model

    ydirout_model(:) = ' '

  END FUNCTION get_model_outputdir
  
  !------------------------------------------------------------------------------

  !=================================================================================
  !=================================================================================

#ifdef NUDGING
  SUBROUTINE set_fdbk_metadata (idom_model, ldebug)

    INTEGER, INTENT(in) :: idom_model
    LOGICAL, INTENT(in) :: ldebug

    INTEGER :: idom, ierr

    CHARACTER(len=cmaxlen) ::  &
         yncglob_institution, & ! originating center name
         yncglob_source         ! program name and version

    REAL(kind=dp) :: &
         hversta, & ! start of verification period in 'model integ. hours'
         hverend    ! end of verification period in 'model integr. hours'

    idom = list_domains_for_radar(idom_model)

    IF (my_radar_id == 0) THEN

      hversta = get_domain_starttime_in_sec(idom_model) / 3600.0_dp + 0.001_dp
      hverend = get_domain_endtime_in_sec(idom_model) / 3600.0_dp

      yncglob_source(:)   = ' '
      yncglob_source      = 'ICON ILAM'
      yncglob_institution(:) = ' '
      yncglob_institution    = 'DWD'

      ! This is one of the workers, so p_patch is defined:
      fdbk_meta_container(idom)%ie_tot  = p_patch(idom_model)%n_patch_cells_g
      fdbk_meta_container(idom)%je_tot  = 1
      fdbk_meta_container(idom)%ke_tot  = p_patch(idom_model)%nlev
      fdbk_meta_container(idom)%pollon  = -180.0_dp
      fdbk_meta_container(idom)%pollat  = 90.0_dp
      fdbk_meta_container(idom)%polgam  = 0.0_dp
      fdbk_meta_container(idom)%dlon    = 0.0_dp
      fdbk_meta_container(idom)%dlat    = 0.0_dp
      fdbk_meta_container(idom)%startlon_tot      = 0.0_dp
      fdbk_meta_container(idom)%startlat_tot      = 0.0_dp
      fdbk_meta_container(idom)%iveri_ens_member  = gribout_config(idom_model)%perturbationNumber
      fdbk_meta_container(idom)%nvers             = gribout_config(idom_model)%localNumberOfExperiment
      fdbk_meta_container(idom)%yglatt_institution(:) = ' '
      fdbk_meta_container(idom)%yglatt_institution    = TRIM(yncglob_institution)
      fdbk_meta_container(idom)%yglatt_source(:) = ' '
      fdbk_meta_container(idom)%yglatt_source    = TRIM(yncglob_source)
      fdbk_meta_container(idom)%hversta  = hversta
      fdbk_meta_container(idom)%hverend  = hverend
    END IF

    IF (num_radar > 1 .AND. num_radario > 0) THEN
      CALL distribute_values_radar (fdbk_meta_container(idom)%ie_tot            , 1, 0, icomm_radar, ierr)
      CALL distribute_values_radar (fdbk_meta_container(idom)%je_tot            , 1, 0, icomm_radar, ierr)
      CALL distribute_values_radar (fdbk_meta_container(idom)%ke_tot            , 1, 0, icomm_radar, ierr)
      CALL distribute_values_radar (fdbk_meta_container(idom)%pollon            , 1, 0, icomm_radar, ierr)
      CALL distribute_values_radar (fdbk_meta_container(idom)%pollat            , 1, 0, icomm_radar, ierr)
      CALL distribute_values_radar (fdbk_meta_container(idom)%polgam            , 1, 0, icomm_radar, ierr)
      CALL distribute_values_radar (fdbk_meta_container(idom)%dlon              , 1, 0, icomm_radar, ierr)
      CALL distribute_values_radar (fdbk_meta_container(idom)%dlat              , 1, 0, icomm_radar, ierr)
      CALL distribute_values_radar (fdbk_meta_container(idom)%startlon_tot      , 1, 0, icomm_radar, ierr)
      CALL distribute_values_radar (fdbk_meta_container(idom)%startlat_tot      , 1, 0, icomm_radar, ierr)
      CALL distribute_values_radar (fdbk_meta_container(idom)%iveri_ens_member  , 1, 0, icomm_radar, ierr)
      CALL distribute_values_radar (fdbk_meta_container(idom)%nvers             , 1, 0, icomm_radar, ierr)
      CALL distribute_values_radar (fdbk_meta_container(idom)%hversta           , 1, 0, icomm_radar, ierr)
      CALL distribute_values_radar (fdbk_meta_container(idom)%hverend           , 1, 0, icomm_radar, ierr)
      CALL distribute_path_radar   (fdbk_meta_container(idom)%yglatt_institution, icomm_radar)
      CALL distribute_path_radar   (fdbk_meta_container(idom)%yglatt_source     , icomm_radar)
    END IF

  END SUBROUTINE set_fdbk_metadata

  SUBROUTINE get_fdbk_metadata ( idom_model, fdbk_meta_data )

    IMPLICIT NONE

    INTEGER, INTENT(in) :: idom_model

    TYPE(fdbk_meta_type), INTENT(out)  :: fdbk_meta_data

    fdbk_meta_data = fdbk_meta_container(list_domains_for_radar(idom_model))

  END SUBROUTINE get_fdbk_metadata
#endif

  !=================================================================================

  ! Meta data for the composite lat/lon grid. At the moment, this has to
  !  be equal to the COSMO grid, because we want to use the standard output facilities
  !  of COSMO for grib output. And these expect the output fields to be
  !  distributed among the compute nodes.
  ! However, the code is written in a way that it is possible and easy
  !  to separate both grids in the future.

  ! For ICON, this is used to define the default before namelist reading
  SUBROUTINE get_composite_metadata ( idom, composite_meta_data )

    IMPLICIT NONE

    ! Parameters
    INTEGER, INTENT(in) :: idom  ! At the moment only a dummy, but in the future the composite grid might
                                 !  depend on the domain

    TYPE(composite_meta_type), INTENT(out)  :: composite_meta_data

    ! .. The following COSMO-DE grid is the default, which can be changed
    !     in the RADARSIM_PARAMS namelist:
    composite_meta_data%ni       = 421
    composite_meta_data%nj       = 461
    composite_meta_data%pollon   = -170.0_dp
    composite_meta_data%pollat   = 40.0_dp
    composite_meta_data%polgam   = 0.0_dp
    composite_meta_data%dlon     = 0.025_dp
    composite_meta_data%dlat     = 0.025_dp
    composite_meta_data%startlon = -5.0_dp
    composite_meta_data%startlat = -5.0_dp
    composite_meta_data%r_earth  = r_earth_dp ! Not changeable via namelist!

  END SUBROUTINE get_composite_metadata

  !=================================================================================

  SUBROUTINE alloc_aux_model_variables (idom, lonline)

    !--------------------------------------------------------------------------
    !
    ! Description:
    !   This subroutine allocates auxiliary fields which are needed locally in the
    !   forward operator and which have to have the same dimensions as the
    !   model variables
    !
    !--------------------------------------------------------------------------

    IMPLICIT NONE

    INTEGER, INTENT(in) :: idom
    LOGICAL, INTENT(in) :: lonline
    INTEGER             :: ni, nj, nk, i, k

    ni = nproma
    nj = p_patch(idom)%nlev
    nk = p_patch(idom)%nblks_c

    IF (.NOT.ALLOCATED(z_radar))    ALLOCATE(z_radar(ni,nj,nk))
    IF (.NOT.ALLOCATED(z_ext))      ALLOCATE(z_ext(ni,nj,nk))
    IF (.NOT.ALLOCATED(vt_radar))   ALLOCATE(vt_radar(ni,nj,nk))
    IF (.NOT.ALLOCATED(rlat))       ALLOCATE(rlat(ni,nk))
    IF (.NOT.ALLOCATED(rlon))       ALLOCATE(rlon(ni,nk))

    DO k=1, nk
      DO i=1, ni
        rlat(i,k) = p_patch(idom) % cells % center(i,k) % lat * raddeg
        rlon(i,k) = p_patch(idom) % cells % center(i,k) % lon * raddeg
      END DO
    END DO

    IF (lonline) THEN
      IF (.NOT.ALLOCATED(vapor_pres)) ALLOCATE(vapor_pres(ni,nj,nk))
      vapor_pres = rho * qv * rv * t
    END IF

  END SUBROUTINE alloc_aux_model_variables

  SUBROUTINE dealloc_aux_model_variables ()

    !--------------------------------------------------------------------------
    !
    ! Description:
    !   This subroutine allocates auxiliary fields which are needed locally in the
    !   forward operator and which have to have the same dimensions as the
    !   model variables
    !
    !--------------------------------------------------------------------------

    IMPLICIT NONE

    IF (ALLOCATED(z_radar))    DEALLOCATE(z_radar)
    IF (ALLOCATED(z_ext))      DEALLOCATE(z_ext)
    IF (ALLOCATED(vt_radar))   DEALLOCATE(vt_radar)
    IF (ALLOCATED(rlat))       DEALLOCATE(rlat)
    IF (ALLOCATED(rlon))       DEALLOCATE(rlon)
    IF (ALLOCATED(vapor_pres)) DEALLOCATE(vapor_pres)

  END SUBROUTINE dealloc_aux_model_variables

  !=================================================================================

  SUBROUTINE get_model_top_height (idom_model, htop_model)

    !--------------------------------------------------------------------------
    !
    ! Description:
    !   This subroutine determines the top height MSL up to which radar
    !   computations should be performed. This is usually the top height
    !   of the model domain in the vertical. This subroutine has to be
    !   called on all compute- and asynchroneous-output-PEs. The latter
    !   require an MPI communication from the compute- to the output-PEs.
    !
    !--------------------------------------------------------------------------

    IMPLICIT NONE

    ! Parameters
    INTEGER, INTENT(in)        :: idom_model
    REAL(kind=dp), INTENT(out) :: htop_model
    
    CHARACTER (LEN=cmaxlen)    :: yzerrmsg
    INTEGER                    :: mpierror

    IF (lcompute_pe_fwo) THEN
      htop_model = vct_a(1+p_patch(idom_model)%nshift_total)
    ELSE
      htop_model = -HUGE(1.0_dp)
    END IF
    IF (num_radar > 1 .AND. num_radario > 0) THEN
      ! Asynchroneous output on separate IO-PEs, therefore communicate htop
      !  to the output PEs using a collective MAX operation:
      yzerrmsg(:) = ' '
      CALL global_values_radar(htop_model, 'MAX', icomm_radar, -1, yzerrmsg, mpierror)
    END IF

  END SUBROUTINE get_model_top_height

  !=================================================================================
  !=================================================================================

  !------------------------------------------------------------------------------

  !=======================================================================================
  !
  ! Functions to initialize the Tmax-parameter (2D (i,j) on model grid) for the
  !  the degree-of-melting parameterization. Tmax is the maximum value in each grid column,
  !  without any neighbourhood search.
  ! Needed by "calc_dbz_vec()" and "calc_fallspeed_vec()" from module "radar_mie_iface_cosmo":
  !
  !=======================================================================================

  SUBROUTINE initialize_tmax_2mom_vec_par(neigh)

    IMPLICIT NONE
    REAL(KIND=dp), INTENT(in) :: neigh

    INTEGER :: iu, ju, i, j, k, im, jm, ierr
    REAL(KIND=dp) :: qithresh, qsthresh, qgthresh, qhthresh

    ! Q - Thresholds for estimation Tmax for degree-of-melting-parameterization
    qithresh = 1d-8
    qsthresh = 1d-8
    qgthresh = 1d-8

    ALLOCATE(Tmax_i_modelgrid(ie_fwo,ke_fwo),Tmax_s_modelgrid(ie_fwo,ke_fwo),Tmax_g_modelgrid(ie_fwo,ke_fwo))

    ! Initial (minimal) values of Tmax_x: 
    CALL init_vari(Tmax_i_modelgrid, 3.0d0 + 273.16)
    CALL init_vari(Tmax_s_modelgrid, 3.0d0 + 273.16)
    CALL init_vari(Tmax_g_modelgrid, 3.0d0 + 273.16)

    ! Tmax is the highest temperature of a point within a neighbourhood (neigh) of the point im, jm, 
    ! where liquid-water-content is still above a small threshold (qthresh).
    ! Due to this, there is a strong dependency on the parameterization of melting applied in cloud physic package.

    ! First, compute max. Temperature where hydrometeors are present within each vertical column:
    DO jm=1,je_fwo
!!$omp parallel do collapse(2) private(k,im)
!$omp parallel do private(k,im)
      DO k=1,ke_fwo
        DO im=1,ie_fwo
          
          IF (qi(im,jm,k) > qithresh .AND. qni(im,jm,k) > 1.0d0) THEN
            Tmax_i_modelgrid(im,k) = MAX(Tmax_i_modelgrid(im,k), t(im,jm,k))
          END IF
          IF (qs(im,jm,k) > qsthresh .AND. qns(im,jm,k) > 1.0d0) THEN
            Tmax_s_modelgrid(im,k) = MAX(Tmax_s_modelgrid(im,k), t(im,jm,k))
          END IF
          IF (qg(im,jm,k) > qgthresh .AND. qng(im,jm,k) > 1.0d-3) THEN
            Tmax_g_modelgrid(im,k) = MAX(Tmax_g_modelgrid(im,k), t(im,jm,k))
          END IF

        END DO
      END DO
!$omp end parallel do
    END DO
    
    ! Maximum acceptable values of Tmax_x_modelgrid: 
!$omp parallel
!$omp workshare
    Tmax_i_modelgrid = MIN(Tmax_i_modelgrid , 10.0d0 + 273.16)
    Tmax_s_modelgrid = MIN(Tmax_s_modelgrid , 10.0d0 + 273.16)
    Tmax_g_modelgrid = MIN(Tmax_g_modelgrid , 15.0d0 + 273.16)
!$omp end workshare
!$omp end parallel
   
    IF (lalloc_qh) THEN
      
      qhthresh = 1d-8
      
      ALLOCATE(Tmax_h_modelgrid(ie_fwo,ke_fwo))
      ! Initial (minimal) value of Tmax_h_modelgrid: 
      CALL init_vari(Tmax_h_modelgrid, 5.0d0 + 273.16)

      DO jm=1,je_fwo
!!$omp parallel do collapse(2) private(im,k)
!$omp parallel do private(im,k)
        DO k=1,ke_fwo
          DO im=1,ie_fwo
            
            IF (qh(im,jm,k) > qhthresh .AND. qnh(im,jm,k) > 1.0d-3) THEN
              Tmax_h_modelgrid(im,k) = MAX(Tmax_h_modelgrid(im,k), t(im,jm,k))
            END IF
            
          END DO
        END DO
!$omp end parallel do
      END DO
      
      ! Maximum acceptable value of Tmax_h_modelgrid: 
!$omp parallel
!$omp workshare
      Tmax_h_modelgrid = MIN(Tmax_h_modelgrid , 30.0d0 + 273.16)
!$omp end workshare
!$omp end parallel

    END IF


    RETURN
  END SUBROUTINE initialize_tmax_2mom_vec_par

  SUBROUTINE initialize_tmax_1mom_vec_par(neigh)

    IMPLICIT NONE
    REAL(KIND=dp), INTENT(in) :: neigh

    INTEGER :: iu, ju, i, j, k, im, jm, ierr
    REAL(KIND=dp) :: qithresh, qsthresh, qgthresh, qhthresh

    !!!!******** neigh is ignored here ********!!!
    

    ! Tmax is the highest temperature of a point within a neighbourhood (neigh) of the point im, jm, 
    ! where liquid-water-content is still above a small threshold (qthresh).
    ! Due to this, there is a strong dependency on the parameterization of melting applied in cloud physic package.

    ! First, compute max. Temperature where hydrometeors are present within each vertical column:
    IF (lalloc_qi) THEN

      ALLOCATE(Tmax_i_modelgrid(ie_fwo,ke_fwo))

      ! Q - Thresholds for estimation Tmax for degree-of-melting-parameterization
      qithresh = 1d-8

      ! Initial (minimal) value of Tmax_g: 
      CALL init_vari(Tmax_i_modelgrid, 3.0d0 + 273.16)

      DO jm=1,je_fwo
!!$omp parallel do collapse(2) private(im,k)
!$omp parallel do private(im,k)
        DO k=1,ke_fwo
          DO im=1,ie_fwo
          
            IF (qi(im,jm,k) >= qithresh) THEN
              Tmax_i_modelgrid(im,k) = MAX(Tmax_i_modelgrid(im,k), t(im,jm,k))
            END IF

          END DO
        END DO
!$omp end parallel do
      END DO

      ! Maximum acceptable value of Tmax_i_modelgrid: 
!$omp parallel
!$omp workshare
      Tmax_i_modelgrid = MIN(Tmax_i_modelgrid , 10.0d0 + 273.16)
!$omp end workshare
!$omp end parallel

    END IF

    IF (lalloc_qs) THEN

      ALLOCATE(Tmax_s_modelgrid(ie_fwo,ke_fwo))
      
      ! Q - Thresholds for estimation Tmax for degree-of-melting-parameterization
      qsthresh = 1d-8

      ! Initial (minimal) value of Tmax_s: 
      CALL init_vari(Tmax_s_modelgrid, 3.0d0 + 273.16)

      DO jm=1,je_fwo
!!$omp parallel do collapse(2) private(k,im)
!$omp parallel do  private(k,im)
        DO k=1,ke_fwo
          DO im=1,ie_fwo
          
            IF (qs(im,jm,k) >= qsthresh) THEN
              Tmax_s_modelgrid(im,k) = MAX(Tmax_s_modelgrid(im,k), t(im,jm,k))
            END IF

          END DO
        END DO
!$omp end parallel do
      END DO
      
      ! Maximum acceptable value of Tmax_s_modelgrid: 
!$omp parallel
!$omp workshare
      Tmax_s_modelgrid = MIN(Tmax_s_modelgrid , 10.0d0 + 273.16)
!$omp end workshare
!$omp end parallel

    END IF

    IF (lalloc_qg) THEN
      
      ALLOCATE(Tmax_g_modelgrid(ie_fwo,ke_fwo))

      ! Q - Thresholds for estimation Tmax for degree-of-melting-parameterization
      qgthresh = 1d-8

      ! Initial (minimal) value of Tmax_g: 
      CALL init_vari(Tmax_g_modelgrid, 3.0d0 + 273.16)

      DO jm=1,je_fwo
!!$omp parallel do collapse(2) private(im,k)
!$omp parallel do private(im,k)
        DO k=1,ke_fwo
           DO im=1,ie_fwo
          
            IF (qg(im,jm,k) >= qgthresh) THEN
              Tmax_g_modelgrid(im,k) = MAX(Tmax_g_modelgrid(im,k), t(im,jm,k))
            END IF
            
          END DO
        END DO
!$omp end parallel do
      END DO
      
      ! Maximum acceptable value of Tmax_g_modelgrid: 
!$omp parallel
!$omp workshare
      Tmax_g_modelgrid = MIN(Tmax_g_modelgrid , 15.0d0 + 273.16)
!$omp end workshare
!$omp end parallel
      
    END IF


    RETURN
  END SUBROUTINE initialize_tmax_1mom_vec_par


  ! Clean-up of Tmax_XX:
  SUBROUTINE finalize_tmax()

    IMPLICIT NONE

    IF (ALLOCATED(Tmax_i_modelgrid)) DEALLOCATE(Tmax_i_modelgrid)
    IF (ALLOCATED(Tmax_s_modelgrid)) DEALLOCATE(Tmax_s_modelgrid)
    IF (ALLOCATED(Tmax_g_modelgrid)) DEALLOCATE(Tmax_g_modelgrid)
    IF (ALLOCATED(Tmax_h_modelgrid)) DEALLOCATE(Tmax_h_modelgrid)

    RETURN
  END SUBROUTINE finalize_tmax

  !============================================================================
  !============================================================================

  !============================================================================
  ! 
  ! Subroutine for specifying an artificial test pattern of hydrometeors
  !  for testing the forward operator during one timestep. Is used in the
  !  Testsuite if ltestpattern_hydrometeors=.TRUE.
  ! 
  ! This routine is for a testpattern on the model grid. It has to be called
  ! AFTER get_model_hydrometeors() and get_model_variables()!!!
  !
  !============================================================================
  
  SUBROUTINE set_testpattern_hydrometeors_mg

    IMPLICIT NONE

    INTEGER :: i, j, k, i_startblk, i_endblk, is, ie, ierr
    REAL(kind=wp), ALLOCATABLE, SAVE :: spatial_modulation(:,:,:), spatial_modulation_2d(:,:), &
         lon_mat(:,:), lat_mat(:,:), t_start(:,:,:)
    REAL(kind=wp) :: lon_min, lon_max, lat_min, lat_max, lon_span, lat_span, &
                     lon_fold, lat_fold
    CHARACTER(len=cmaxlen)  :: yerrmsg

    WRITE (*,'(a,i4,a,f10.0,a)') '*** WARNING: Experimental run with pre-specified test pattern of '//&
         'hydrometeor contents in radar_interface.f90 on my_cart_id_fwo = ', my_cart_id_fwo, &
         ' at time ', get_model_time_sec(), ' ! ***'

    IF (.NOT.ALLOCATED(spatial_modulation)) THEN

      ALLOCATE(spatial_modulation(ie_fwo,je_fwo,ke_fwo)) ! 3D
      ALLOCATE(spatial_modulation_2d(ie_fwo,ke_fwo))    ! 2D horiztonal

      ALLOCATE(lon_mat(ie_fwo,ke_fwo), lat_mat(ie_fwo,ke_fwo))

      ! .. 4x3 checkerboard pattern on total domain, provided the domain does not cross the date line:

      lon_mat = -HUGE(1.0_dp)
      lat_mat = -HUGE(1.0_dp)

      i_startblk = p_patch_for_testpattern % cells % start_block(grf_bdywidth_c+1)
      i_endblk   = p_patch_for_testpattern % cells % end_block(min_rlcell_int)

      DO k = i_startblk, i_endblk
        CALL get_indices_c(p_patch_for_testpattern, k, i_startblk, i_endblk, is, ie, grf_bdywidth_c+1, min_rlcell_int)
        DO i = is, ie
          lon_mat(i,k) = p_patch_for_testpattern % cells % center(i,k) % lon * raddeg
          lat_mat(i,k) = p_patch_for_testpattern % cells % center(i,k) % lat * raddeg
        END DO
      END DO

      lon_min = MINVAL(lon_mat, mask=lon_mat>-900.0_wp)
      lon_max = MAXVAL(lon_mat, mask=lon_mat>-900.0_wp)
      lat_min = MINVAL(lat_mat, mask=lat_mat>-900.0_wp)
      lat_max = MAXVAL(lat_mat, mask=lat_mat>-900.0_wp)

      ! Get the global min/max of lat/lon:
      IF (num_compute_fwo > 1) THEN
        CALL global_values_radar(lon_min, 'MIN', icomm_cart_fwo, -1, yerrmsg, ierr)
        CALL global_values_radar(lon_max, 'MAX', icomm_cart_fwo, -1, yerrmsg, ierr)
        CALL global_values_radar(lat_min, 'MIN', icomm_cart_fwo, -1, yerrmsg, ierr)
        CALL global_values_radar(lat_max, 'MAX', icomm_cart_fwo, -1, yerrmsg, ierr)
      END IF

      lon_span = (lon_max - lon_min) / 4.0_wp
      lat_span = (lat_max - lat_min) / 3.0_wp

      spatial_modulation_2d(:,:) = 0.0_wp
      DO k=1, ke_fwo
        DO i=1, ie_fwo
          IF (lon_mat(i,k) > -HUGE(1.0_dp)+1e-20_dp .AND. lat_mat(i,k) > -HUGE(1.0_dp)+1e-20_dp) THEN
            lon_fold = MODULO(lon_mat(i,k)-lon_min, lon_span)
            lat_fold = MODULO(lat_mat(i,k)-lat_min, lat_span)
            spatial_modulation_2d(i,k) = 1.0 + 0.25 * ( &
                 (lon_fold-0.5_wp*lon_span)/(0.5_wp*lon_span) + &
                 (lat_fold-0.5_wp*lat_span)/(0.5_wp*lat_span)   &
                 )
          END IF
        END DO
      END DO
      
      ! Distribute to vertical levels (j):
      DO j=1, je_fwo
        spatial_modulation(:,j,:) = spatial_modulation_2d(:,:)
      END DO

      ALLOCATE(t_start(ie_fwo,je_fwo,ke_fwo)) ! 3D
      t_start(:,:,:) = t(:,:,:)

    END IF

    ! Vertical profiles for model variables:
    qc(:,1:je_fwo-10,:) = 0.5e-3_wp * spatial_modulation(:,1:je_fwo-10,:)
    qc(:,je_fwo-9:je_fwo,:) = 0.0e-3_wp * spatial_modulation(:,je_fwo-9:je_fwo,:)
    IF (lalloc_qi) qi = 0.1e-3_wp * spatial_modulation
    IF (ASSOCIATED(qr)) THEN
      qr = 1.5e-3_wp * spatial_modulation
      WHERE (t_start < 248.16) qr = 0.0_wp
    END IF
    IF (lalloc_qs) THEN
      qs = 1.0e-3_wp * spatial_modulation
      WHERE (t_start > 282.16) qs = 0.0_wp
    END IF
    IF (lalloc_qg) THEN
      qg = 2.0e-3_wp * spatial_modulation
      WHERE (t_start > 284.16) qg = 0.0_wp
      WHERE (t_start < 238.16) qg = 0.0_wp
    END IF

    IF (lalloc_qh) THEN
      qh = 4.0e-3_wp * spatial_modulation
      WHERE (t_start > 286.16) qh = 0.0_wp
      WHERE (t_start < 243.16) qh = 0.0_wp
    END IF
    IF (ASSOCIATED(qnc)) THEN
      WHERE (rho0 > 1e-7_wp) qnc = 5e7 / rho0 * spatial_modulation
    END IF
    IF (ASSOCIATED(qni)) THEN 
      WHERE (rho0 > 1e-7_wp) qni = 5e4 / rho0 * spatial_modulation
    END IF
    IF (ASSOCIATED(qnr)) THEN
      WHERE (rho0 > 1e-7_wp) qnr = 1e3 / rho0 * spatial_modulation
      WHERE (qr < 1e-7_wp) qnr = 0.0_wp
    END IF
    IF (ASSOCIATED(qns)) THEN
      WHERE (rho0 > 1e-7_wp) qns = 1e4 / rho0 * spatial_modulation
      WHERE (qs < 1e-7_wp) qns = 0.0_wp
    END IF
    IF (ASSOCIATED(qng)) THEN
      qng = 1e4_wp * spatial_modulation
      WHERE (qg < 1e-7_wp) qng = 0.0_wp
    END IF
    IF (ASSOCIATED(qnh)) THEN
      qnh = 0.5e3_wp * spatial_modulation
      WHERE (qh < 1e-7_wp) qnh = 0.0_wp
    END IF

  END SUBROUTINE set_testpattern_hydrometeors_mg

!================================================================================

  !============================================================================
  ! 
  ! Functions to return the index limits of the vertical grid, in order
  !  to provide a flexible means of changing the level ordering when coupling
  !  to another model
  ! 
  !============================================================================

  PURE FUNCTION bottomlevel () RESULT (kbot)
    INTEGER :: kbot
    kbot = je_fwo
  END FUNCTION bottomlevel

  PURE FUNCTION bottomlevel_stag () RESULT (kbot)
    INTEGER :: kbot
    kbot = je_fwo + 1
  END FUNCTION bottomlevel_stag

  PURE FUNCTION toplevel () RESULT (ktop)
    INTEGER :: ktop
    ktop = 1
  END FUNCTION toplevel

  PURE FUNCTION levelincr () RESULT (kincr)
    INTEGER :: kincr
    kincr = -1
  END FUNCTION levelincr


  !============================================================================
  ! 
  ! Compute the level index below the level k in a 3D model field at mass positions:
  ! 
  !============================================================================

  ELEMENTAL FUNCTION one_level_down (k) RESULT (ko)

    IMPLICIT NONE

    INTEGER, INTENT(in) :: k
    INTEGER             :: ko

    ko = MIN(k+1, je_fwo)

  END FUNCTION one_level_down

  !============================================================================
  ! 
  ! Compute the level index above the level k in a 3D model field:
  ! 
  !============================================================================

  ELEMENTAL FUNCTION one_level_up (k) RESULT (ku)

    IMPLICIT NONE

    INTEGER, INTENT(in) :: k
    INTEGER             :: ku

    ku = MAX(k-1, 1)

  END FUNCTION one_level_up

!================================================================================

  !------------------------------------------------------------------------------

  !============================================================================
  ! 
  ! Subroutines related to interpolations from the model grid to the radar bins
  ! 
  !============================================================================


  !============================================================================
  ! 
  ! Find the center of the global domain by inspecting all local domains.
  ! Then do collective communication over all workers+radarIO-PEs.
  ! Has to be called on all radar PEs, including the radarIO-PEs.
  ! 
  !============================================================================
  
  SUBROUTINE get_domaincenter_global (idom_model, ldebug)

    INTEGER, INTENT(in) :: idom_model
    LOGICAL, INTENT(in) :: ldebug

    INTEGER             :: idom_fwo, i, k, i_startblk, i_endblk, is, ie, mpierror
    REAL(KIND=dp)       :: clon, clat, clon_min, clon_max, clat_min, clat_max, tmp
    CHARACTER(len=cmaxlen)  :: yzerrmsg

    ! 1) Find the local bounding box on the PE subdomain including the halo cells for later double-assignment search:
    !    ------------------------------------------------------------------------------------------------------------

    IF (ldebug) WRITE(*,*) 'get_domaincenter_global on proc ', my_radar_id

    clon_min = HUGE(1.0_dp)
    clon_max = -HUGE(1.0_dp)
    clat_min = HUGE(1.0_dp)
    clat_max = -HUGE(1.0_dp)
    IF ( lcompute_pe_fwo) THEN

      i_startblk = p_patch(idom_model) % cells % start_block(grf_bdywidth_c+1)
      i_endblk   = p_patch(idom_model) % cells % end_block(min_rlcell_int)

      DO k = i_startblk, i_endblk

        CALL get_indices_c(p_patch(idom_model), k, i_startblk, i_endblk, is, ie, grf_bdywidth_c+1, min_rlcell_int)

        DO i = is, ie

          ! local grid cell lon/lat:
          clon = p_patch(idom_model) % cells % center(i,k) % lon * raddeg
          clat = p_patch(idom_model) % cells % center(i,k) % lat * raddeg
          
          ! ... min/max bounding box:
          clon_min = MIN(clon_min, clon)
          clon_max = MAX(clon_max, clon)
          clat_min = MIN(clat_min, clat)
          clat_max = MAX(clat_max, clat)

        END DO
      END DO

      ! Estimate if local domain may be across the date line. If yes,
      !  make the min the new max by adding 360 degrees:
      IF (clon_max-clon_min > 180.0_dp) THEN
        tmp = clon_max
        clon_max = clon_min + 360.0_dp
        clon_min = tmp
      END IF
    END IF

    ! In parallel runs, find global min/max of local BBs and distribute to all radar PEs,
    !  including the asynchroneous radar IO PEs:
    IF (num_radar > 1) THEN
      yzerrmsg(:) = ' '
      CALL global_values_radar(clon_min, 'MIN', icomm_radar, -1, yzerrmsg, mpierror)
      CALL global_values_radar(clon_max, 'MAX', icomm_radar, -1, yzerrmsg, mpierror)
      CALL global_values_radar(clat_min, 'MIN', icomm_radar, -1, yzerrmsg, mpierror)
      CALL global_values_radar(clat_max, 'MAX', icomm_radar, -1, yzerrmsg, mpierror)
    END IF
    
    ! Estimate if the total domain may be across the date line. If yes,
    !  make the min the new max by adding 360 degrees:
    IF (clon_max-clon_min > 180.0_dp) THEN
      tmp = clon_max
      clon_max = clon_min + 360.0_dp
      clon_min = tmp
    END IF

    ! The domain center is defined as the mean value of the global bounding box margins
    !  and is stored on internal fields of radar_interface.f90:
    idom_fwo =  list_domains_for_radar(idom_model)
    dom_center_lon(idom_fwo) = MODULO(0.5_dp*(clon_max+clon_min), 360.0_dp)
    dom_center_lat(idom_fwo) =        0.5_dp*(clat_max+clat_min)

    IF (ldebug) THEN
      WRITE (*, '(a,i3,a,/,a,f12.5,/,a,f12.5)') 'INFO get_domaincenter_global(): '// &
           'center of radar-active domain ', idom_model, ':', &
           '  lon_center = ', dom_center_lon(idom_fwo), &
           '  lat_center = ', dom_center_lat(idom_fwo)
    END IF

    IF (ldebug) WRITE(*,*) 'Done with get_domaincenter_global on proc ', my_radar_id

  END SUBROUTINE get_domaincenter_global

  !============================================================================
  !============================================================================
  
  SUBROUTINE setup_model2radarbins_vec (idom_model, rs_meta, nobsmax, lon_r, lat_r, alt_r, &
       nobs, nobs_above_surface, ind_intp, w_intp, hl_loc)

    ! INPUT variables:
    !-----------------

    INTEGER, INTENT(in)                  :: idom_model
    TYPE(radar_meta_type), INTENT(in)    :: rs_meta

    INTEGER, INTENT(in)       :: nobsmax   ! number of data in rlon, rlat, ...
    REAL(KIND=dp), INTENT(in) :: lon_r(nobsmax), lat_r(nobsmax), alt_r(nobsmax)

    ! OUTPUT variables:
    !-----------------

    INTEGER, INTENT(out) :: nobs, nobs_above_surface  ! actual numbers of obs on this PE
    INTEGER, POINTER     :: ind_intp(:,:)
    REAL(KIND=dp), POINTER ::   &
         w_intp(:,:) ,& ! array of horizontal interpolation weights for each observation (first dimension)
                        ! and each spatial direction i,j,k (second dimension) 
         hl_loc(:)  ! radar ray heights according to 4/3 earth model

    ! Local variables:
    !-----------------

    INTEGER        :: irp, &
         nrp, irploc, &
         nobs_above_sfc, nobs_below_sfc

    REAL    (KIND=dp), ALLOCATABLE :: &
         alt_rtmp(:) ,& ! array of altitudes above sea level for temporary radar point (nrp)
         w_intptmp(:,:) ! array of horizontal interpolation weights for each temporary radar point

    INTEGER , ALLOCATABLE ::&
         idx(:), blk(:),  &
         irp_tmp(:) , &
         ind_intptmp(:,:)   ! the second dimension consists of:
                            ! 1:     the continuous number of the radar bin
                            ! 2,3,4: the model indices in the 3 coordinate directions
 
    ALLOCATE(irp_tmp(nobsmax)) 
    ALLOCATE(idx(nobsmax), blk(nobsmax))
    CALL init_vari(irp_tmp(:), -1      )
    CALL init_vari(idx(:)    , -HUGE(1))
    CALL init_vari(blk(:)    , -HUGE(1))

    ! calculate idx and blk indices of enclosing ICON grid cell:
!$omp parallel do
    DO irp = 1, nobsmax

      CALL geo2cell_index2d (cindex_grid(list_domains_for_radar(idom_model)), p_patch(idom_model), &
                             lon_r(irp), lat_r(irp), idx(irp), blk(irp))

    ENDDO
!$omp end parallel do

    ! To gain the indices i,j of model grids southwest to radar points and 
    ! check if the radar point is in the local domain, based on lon i and lat j of radar points
    nrp = 0

    DO irp = 1, nobsmax

      IF ( idx(irp) > -HUGE(1) .AND. blk(irp) > -HUGE(1) .AND. &
           (alt_r(irp)  <= htop) ) THEN
        
        nrp = nrp + 1
        
        ! store the continuously running index of the point w.r.t. azi,range,ele, which has been
        !  set up by a call to sub2ind3D() in rad2geo_const_vec():
! This hinders OMP parallelization!        
        irp_tmp(nrp) = irp
        
      END IF
      
    END DO

    ALLOCATE(ind_intptmp(nrp,4))                    ! 4 = number of model coordinates + 1 (for cont. 1D-index of radar coords)
    ALLOCATE(w_intptmp(nrp,ninterp_weights_model))  ! 1 = number of needed interpolation weights for ICON
    ALLOCATE(alt_rtmp(nrp))

    CALL init_vari(w_intptmp, -1.0_dp)
    CALL init_vari(ind_intptmp, -1)
    CALL init_vari(alt_rtmp, -1.0_dp)


!CDIR NODEP,VOVERTAKE,VOB
!$omp parallel do private(irploc)
    DO irp = 1, nrp

      irploc = irp_tmp(irp)

      ! store continuously running index over radar coordinates and corresponging horiz. indices idx and blk
      !  of enclosing ICON cell for later local sorting:
      ind_intptmp(irp,1) = irploc
      ind_intptmp(irp,2) = idx(irploc)
      ind_intptmp(irp,4) = blk(irploc)
      
      ! compute interpolation weights in i and j directions (THIS IS MODEL SPECIFIC!)
      !!$ --> Nothing to be done for nearest neighbour in ICON!

      alt_rtmp(irp) = alt_r(irploc)
      
    ENDDO
!$omp end parallel do
    
    DEALLOCATE(irp_tmp)

    ! To gain the indices k of model grids above radar points (THIS IS MODEL SPECIFIC)
    IF (nrp > 0) THEN
      CALL calc_vert_weight_vec(nrp,alt_rtmp(1:nrp),ind_intptmp(1:nrp,2),ind_intptmp(1:nrp,4),&
                                ind_intptmp(1:nrp,3),w_intptmp(1:nrp,1))
    END IF

    ! Obtain the final set of observable radar points above the surface in the local domain

    ! Allocate and fill final arrays of interpolation coefficients and indices of current radar station
    !  (for safety: deallocate first, if associated)
    IF (ASSOCIATED(w_intp))   DEALLOCATE (w_intp)
    IF (ASSOCIATED(ind_intp)) DEALLOCATE (ind_intp)
    IF (ASSOCIATED(hl_loc))   DEALLOCATE (hl_loc)
    ALLOCATE(w_intp(nrp, ninterp_weights_model))
    ALLOCATE(ind_intp(nrp, 2)) ! index 1: combined model index interpolation; index 2: combined radar index
    ALLOCATE(hl_loc(nrp))

    ! First, sort radar points above surface into the weight fields:
    nobs_above_sfc = 0 
!CDIR VOVERTAKE,VOB
    DO irp = 1, nrp

      ! Check if the radar point is in the local domain, based on height k (ind_intptmp(irp,3) <= 0 means below the surface)
      IF (ind_intptmp(irp,3) > 0) THEN

        nobs_above_sfc = nobs_above_sfc + 1

! This hinders OMP parallelization!        
        w_intp(nobs_above_sfc,1) = w_intptmp(irp,1)
        
        CALL sub2ind3D(ind_intptmp(irp,2), ind_intptmp(irp,3), ind_intptmp(irp,4), &
                       nproma, je_fwo, ind_intp(nobs_above_sfc,1))

        ind_intp(nobs_above_sfc,2) = ind_intptmp(irp,1)

        hl_loc(nobs_above_sfc) = alt_rtmp(irp)

      ENDIF

    ENDDO

    ! Number of obs above surface:
    nobs_above_surface = nobs_above_sfc
    ! Total number of obs (including points below surface):
    nobs               = nrp

    ! Then, keep track of the positions of the points below surface:
    !   (Here, later the missing value -988.88 will be set for radial wind and reflectivity)
    nobs_below_sfc = nobs_above_sfc
!CDIR VOVERTAKE,VOB
    DO irp = 1, nrp
      
      IF (ind_intptmp(irp,3) == -1) THEN
        nobs_below_sfc = nobs_below_sfc + 1
        w_intp(nobs_below_sfc,1)   = -999.99_dp
        ind_intp(nobs_below_sfc,1) = -999
        ind_intp(nobs_below_sfc,2) = ind_intptmp(irp,1)
        hl_loc(nobs_below_sfc)     = alt_rtmp(irp)          
      END IF

    END DO

    ! deallocate auxiliary fields
    DEALLOCATE(w_intptmp)
    DEALLOCATE(ind_intptmp)
    DEALLOCATE(alt_rtmp)

  END SUBROUTINE setup_model2radarbins_vec

!==============================================================================

  SUBROUTINE setup_model2azislices_vec (idom_model, rs_grid, lon_g, lat_g)


    ! INPUT variables:
    !-----------------

    INTEGER, INTENT(in)       :: idom_model
    REAL(KIND=dp), INTENT(in) :: lon_g(:,:), lat_g(:,:)

    ! INOUT variables:
    !-----------------

    TYPE(radar_grid_type), INTENT(inout) :: rs_grid

    ! Local variables:
    !-----------------

    INTEGER        :: ngrd, ngrdmax, i, j, k, m, n, o, offset_i, offset_j, idx, blk

    REAL    (KIND=dp)          :: &
         wi,& ! interpolation weight in i-direction
         wj,& ! interpolation weight in j-direction
         rlon_g, & ! rotated longitude for a auxiliary grid
         rlat_g    ! rotated latitude for a auxiliary grid

    REAL    (KIND=dp),  ALLOCATABLE :: &
         hl(:)     ! array of geographical heights for each auxiliary grid  

    INTEGER , ALLOCATABLE     :: &
         ind_intptmp(:,:)  ! array of indices for each auxiliary grid (first dimension)
                           ! the second dimension consists of:
                           ! 1:     the continuous index of the model grid cell associated with the observation
                           ! 2:     the continuous index representing aux grid points in
                           !        azimuthal, arc length and vertical direction (m,n,k)


    ! allocate and initialize local aux arrays with maximum possible vector length:
    ngrdmax = (rs_grid%nal+1) * (rs_grid%naz_nbl) *  je_fwo

    ALLOCATE(ind_intptmp(ngrdmax,2))
    ind_intptmp = -1

    ! loop over all grid points
    ngrd = 0     

    ALLOCATE(hl(ngrdmax))
    CALL init_vari(hl, -999.99_dp)

    DO n = 1, rs_grid%nal+1      ! loop over arc length
      DO m = 1, rs_grid%naz_nbl  ! loop over azimuths

        ! calculate index idx and blk of ICON cell that contains the radar point
        CALL geo2cell_index2d (cindex_grid(list_domains_for_radar(idom_model)), p_patch(idom_model), &
                               lon_g(m,n), lat_g(m,n), idx, blk)

        ! find points within the local computational processor domain:
        IF ( idx > -HUGE(1) .AND. blk > -HUGE(1) ) THEN

          DO k = je_fwo,1,-1         ! loop over heights

            ngrd = ngrd + 1

! This hinders OMP parallelization!
            hl(ngrd) = hfl(idx, k, blk)

            ! continuous index representing the upper ICON cell idx,k,blk
            CALL sub2ind3D(idx, k, blk, nproma, je_fwo, ind_intptmp(ngrd,1))

            ! continuous index representing aux grid point azi, arc dist, height
            CALL sub2ind3D(m, n, k, rs_grid%naz_nbl, rs_grid%nal+1, ind_intptmp(ngrd,2))

          END DO

        END IF

      END DO
    END DO

    ! allocate arrays of interpolation coefficients and indices of current radar station
    !  (for safety: deallocate first, if associated)
    IF (ASSOCIATED(rs_grid%ind_intp)) DEALLOCATE (rs_grid%ind_intp)
    IF (ASSOCIATED(rs_grid%hl_grd))   DEALLOCATE (rs_grid%hl_grd)

    ALLOCATE(rs_grid%w_intp(1,2))   ! just a dummy in case of ICON
    ALLOCATE(rs_grid%ind_intp(ngrd,2))
    ALLOCATE(rs_grid%hl_grd(ngrd))

    rs_grid%ngrd = ngrd

    IF (ngrd > 0) THEN

      ! copy values from auxiliary variables
      rs_grid%w_intp   = 1.0_dp    ! just a dummy in case of ICON
!$omp parallel workshare
      rs_grid%ind_intp = ind_intptmp(1:ngrd,:)
      rs_grid%hl_grd   = hl(1:ngrd)
!$omp end parallel workshare

    END IF

    DEALLOCATE(ind_intptmp)

  END SUBROUTINE setup_model2azislices_vec

!==============================================================================

  SUBROUTINE calc_vert_weight(alt,hhl,hfl,ke,kout,wk)

    !------------------------------------------------------------------------------
    !
    ! Description: Calculation of vertical interpolation weight and level index
    !              above the radar point for the interpolation of the model grids
    !              to the radar points. This routine assumes that all input radar
    !              points are horizontally within the model domain.
    !
    !------------------------------------------------------------------------------
    ! 
    ! Subroutine / Function arguments
    !------------------------------------------------------------------------------
    ! Parameter list:
    INTEGER , INTENT (IN)     ::        &
         ke               ! vertical dimension of model levels

    REAL (KIND=dp), INTENT (IN)           ::        &
         alt,           & ! height a.s.l of radar point
         hhl(je_fwo+1), & ! model half level height of the enclosing ICON grid cell
         hfl(je_fwo)      ! model full level height of the enclosing ICON grid cell

    INTEGER , INTENT (OUT)    ::        &
         kout             ! lowest level index above radar point
                          ! -1 if radar point below topography
                          ! ke if radar point above topo and below first model full level

    REAL (KIND=dp), INTENT (OUT)          ::        &
         wk               ! vertical interpolation weight
                          ! -1.0 if radar point below topography
                          !  0.0 if radar point above topo and below first model full level

    ! Local scalars:
    INTEGER  :: k
    REAL    (KIND=dp)    :: hfl_up, hfl_low, hsurf_avg

    ! Local arrays:

    !- End of header
    !==============================================================================


    !initialize vertical interpolation weight and level index
    wk = -1.0_dp
    kout = -1

    ! find vertical index below height of radar point

    hfl_up = hfl(je_fwo)

    ! if radar point is below surface return undefined
    hsurf_avg = hhl(je_fwo+1)

    IF (alt < hsurf_avg) RETURN

    ! if radar point is below lowest full level
    IF (alt < hfl_up) THEN
      wk = 0.0_dp
      kout = je_fwo
      RETURN 
    END IF

    ! loop over levels, start at the surface and proceed upwards
    DO k = je_fwo,1,-1

      hfl_low = hfl_up
      hfl_up  = hfl(k)

      IF (hfl_up > alt) THEN
        ! radar point found:
        wk = (alt-hfl_up)/(hfl_low-hfl_up)
        kout = k
        EXIT
      END IF

    END DO

  END SUBROUTINE calc_vert_weight

  SUBROUTINE calc_vert_weight_vec(np,alt,idxin,blkin,kout,wk)

    !------------------------------------------------------------------------------
    !
    ! Description: Calculation of vertical interpolation weight and level index
    !              above the radar point for the interpolation of the model grids
    !              to the radar points. This routine assumes that all input radar
    !              points are horizontally within the model domain.
    !
    !------------------------------------------------------------------------------
    ! 
    ! Subroutine / Function arguments
    !------------------------------------------------------------------------------
    ! Parameter list:


    INTEGER , INTENT (IN)     ::   &
         np,                       &             
         idxin(np),                &
         blkin(np)            

    REAL (KIND=dp), INTENT (IN)           ::        &
         alt(np)              ! height a.s.l of radar point

    INTEGER , INTENT (OUT)    ::        &
         kout(np)             ! lowest level index above radar point
                              ! -1 if radar point below topography
                              ! ke if radar point above topo and below first model full level

    REAL (KIND=dp), INTENT (OUT)          ::        &
         wk(np)               ! vertical interpolation weight
                              ! -1.0 if radar point below topography
                              !  0.0 if radar point above topo and below first model full level

    INTEGER  :: k, irp, flag(np) ! use flag to check if the radar point has already been done
    REAL    (KIND=dp)    :: hfl_low, hsurf_avg

    ! Local arrays:
    REAL    (KIND=dp)    :: hfl_up(np) ! full level heights of the 4 grid points
                                                                   ! surrounding the radar point

    !- End of header
    !==============================================================================

    !initialize vertical interpolation weight and level index
    wk = -1.0_dp
    kout = -1
    flag = 0

    ! Check if the radar point is below topography.
    ! If yes, then flag = 1 (wk,kout = default value); Elseif below the first main level, then flag = 1, wk = 0, kout = ke.
!$omp parallel private (hsurf_avg)
!$omp do 
    DO irp = 1, np

      hfl_up(irp) = hfl(idxin(irp), je_fwo, blkin(irp))

      hsurf_avg = hhl(idxin(irp), je_fwo+1, blkin(irp))

      IF (alt(irp) < hsurf_avg) THEN

        flag(irp) = 1

      ELSEIF (hfl_up(irp) > alt(irp)) THEN

        wk(irp) = 0.0_dp
        kout(irp) = je_fwo
        flag(irp) = 1

      ENDIF

    ENDDO
!$omp end do
!$omp end parallel

    ! Loop over k from second lowest level to the top
    ! Determine the k indice above the radar point
    DO k = je_fwo-1,1,-1

!$omp parallel private (hfl_low)
!$omp do 
      DO irp =1, np

        IF (flag(irp) == 0) THEN

          hfl_low = hfl_up(irp)

          hfl_up(irp) = hfl(idxin(irp), k, blkin(irp))

          IF (hfl_up(irp) > alt(irp)) THEN

            ! radar point found:
            wk(irp) = (alt(irp)-hfl_up(irp))/(hfl_low-hfl_up(irp))

            kout(irp) = k

            flag(irp) = 1

          ENDIF

        ENDIF

      ENDDO
!$omp end do
!$omp end parallel

    ENDDO

  END SUBROUTINE calc_vert_weight_vec

!==============================================================================

  SUBROUTINE geo2cell_index2d (cidx_grid, p_patch_dom, lon_geo, lat_geo, idx, blk)

    IMPLICIT NONE
    TYPE(t_cindex_grid), INTENT(in) :: cidx_grid
    TYPE(t_patch)      , INTENT(in) :: p_patch_dom
    REAL(KIND=dp), INTENT(in)       :: lon_geo, lat_geo  ! geogr. lon/lat
    INTEGER, INTENT(out)            :: idx, blk

    REAL(KIND=dp)                   :: lon_rot, lat_rot  ! rotated lon/lat
    REAL(KIND=dp)                   :: tmplon, tmplat    ! geogr. lon/lat
    REAL(KIND=dp)                   :: dist, disttmp
    INTEGER                         :: i, k, ii,kk, idxtmp, blktmp

    IF ( ALLOCATED(cidx_grid%cind) ) THEN

      CALL geo2rotll_coord (lon_geo, lat_geo, &
                            cidx_grid%pollon, cidx_grid%pollat, cidx_grid%polgam, &
                            lon_rot, lat_rot)

      i = FLOOR( (lon_rot-cidx_grid%startlon) / cidx_grid%dlon ) + 1
      k = FLOOR( (lat_rot-cidx_grid%startlat) / cidx_grid%dlat ) + 1
      
      IF ( i > 0 .AND. i <= cidx_grid%nlon .AND. &
           k > 0 .AND. k <= cidx_grid%nlat ) THEN

        idx = cidx_grid%idx(i,k)
        blk = cidx_grid%blk(i,k)

        IF (idx > -HUGE(1)) THEN
          ! search in the 8 neighbours if any of their nearest cells might be closer:
          tmplon = p_patch_dom % cells % center(idx,blk) % lon * raddeg
          tmplat = p_patch_dom % cells % center(idx,blk) % lon * raddeg
          dist = geo_dist(lon_geo, lat_geo, tmplon, tmplat, cidx_grid%r_earth, 0.0_dp)
          DO kk = MAX(k-1,1), MIN(k+1,cidx_grid%nlat)
            DO ii = MAX(i-1,1), MIN(i+1,cidx_grid%nlon)
              idxtmp = cidx_grid%idx(ii,kk)
              blktmp = cidx_grid%blk(ii,kk)
              IF ( idxtmp  > -HUGE(1) ) THEN
                tmplon = p_patch_dom % cells % center(idxtmp,blktmp) % lon * raddeg
                tmplat = p_patch_dom % cells % center(idxtmp,blktmp) % lon * raddeg
                disttmp = geo_dist(lon_geo, lat_geo, tmplon, tmplat, cidx_grid%r_earth, 0.0_dp)
                IF ( disttmp < dist ) THEN
                  dist = disttmp
                  idx  = idxtmp
                  blk  = blktmp
                END IF
              END IF
            END DO
          END DO
          
          ! if the indentified nearest cell is not among the interior cells,
          !  it will be found on another PE and is flagged as not found here:
          IF (.NOT. p_patch_dom % cells % decomp_info % owner_mask(idx,blk)) THEN
            idx = -HUGE(1)
            blk = -HUGE(1)
          END IF

        END IF

      ELSE
        idx = -HUGE(1)
        blk = -HUGE(1)
      END IF

    ELSE

      ! Not on a worker PE or no overlap with radar-covered area or setup_auxgrid_for_cellindex() has not yet been called, so we
      !  cannot give a correct index:

      idx = -HUGE(1)
      blk = -HUGE(1)

    END IF

  END SUBROUTINE geo2cell_index2d

  FUNCTION geo2cell_index1d (cidx_grid, p_patch_dom, lon_geo, lat_geo) RESULT (index1d)

    IMPLICIT NONE
    TYPE(t_cindex_grid), INTENT(in) :: cidx_grid
    TYPE(t_patch)      , INTENT(in) :: p_patch_dom
    REAL(KIND=dp), INTENT(in)       :: lon_geo, lat_geo
    INTEGER :: index1d

    REAL(KIND=dp)                   :: lon_rot, lat_rot    ! rotated lon/lat
    INTEGER                         :: i, k, idx, blk

    IF ( ALLOCATED(cidx_grid%cind) ) THEN

      CALL geo2rotll_coord (lon_geo, lat_geo, &
                            cidx_grid%pollon, cidx_grid%pollat, cidx_grid%polgam, &
                            lon_rot, lat_rot)

      CALL geo2cell_index2d (cidx_grid, p_patch_dom, lon_geo, lat_geo, idx, blk)

      IF (idx > -HUGE(1)) THEN
        i = FLOOR( (lon_rot-cidx_grid%startlon) / cidx_grid%dlon ) + 1
        k = FLOOR( (lat_rot-cidx_grid%startlat) / cidx_grid%dlat ) + 1
        index1d = cidx_grid%cind(i,k)
      ELSE
        index1d = -HUGE(1)
      END IF

    ELSE

      ! Not on a worker PE or setup_auxgrid_for_cellindex() has not yet been called, so we
      !  cannot give a correct index:

      index1d = -HUGE(1)     

    END IF

  END FUNCTION geo2cell_index1d

  FUNCTION geo2cell_index1d_glob (cidx_grid,  p_patch_dom, lon_geo, lat_geo) RESULT (index1d)

    IMPLICIT NONE
    TYPE(t_cindex_grid), INTENT(in) :: cidx_grid
    TYPE(t_patch)      , INTENT(in) :: p_patch_dom
    REAL(KIND=dp), INTENT(in)       :: lon_geo, lat_geo
    INTEGER                         :: index1d

    REAL(KIND=dp)                   :: lon_rot, lat_rot    ! rotated lon/lat
    INTEGER                         :: i, k, idx, blk

    IF ( ALLOCATED(cidx_grid%cind_glob) ) THEN

      CALL geo2rotll_coord (lon_geo, lat_geo, &
                            cidx_grid%pollon, cidx_grid%pollat, cidx_grid%polgam, &
                            lon_rot, lat_rot)

      CALL geo2cell_index2d (cidx_grid, p_patch_dom, lon_geo, lat_geo, idx, blk)

      IF (idx > -HUGE(1)) THEN
        i = FLOOR( (lon_rot-cidx_grid%startlon) / cidx_grid%dlon) + 1
        k = FLOOR( (lat_rot-cidx_grid%startlat) / cidx_grid%dlat) + 1
        index1d = cidx_grid%cind_glob(i,k)
      ELSE
        index1d = -HUGE(1)
      END IF

    ELSE

      ! Not on a worker PE or setup_auxgrid_for_cellindex() has not yet been called, so we
      !  cannot give a correct index:

      index1d = -HUGE(1)

    END IF

  END FUNCTION geo2cell_index1d_glob

  SUBROUTINE rotll2geo_coord (lon_rot, lat_rot, pollon, pollat, polgam, lon_geo, lat_geo)
    
    IMPLICIT NONE

    REAL(KIND=dp), INTENT(in)       :: lon_rot, lat_rot    ! rotated lon/lat
    REAL(KIND=dp), INTENT(in)       :: pollon, pollat, polgam
    REAL(KIND=dp), INTENT(out)      :: lon_geo, lat_geo    ! geogr. lon/lat

    lon_geo = rlarot2rla(lat_rot,lon_rot,pollat,pollon,polgam)
    lat_geo = phirot2phi(lat_rot,lon_rot,pollat,pollon,polgam)


  END SUBROUTINE rotll2geo_coord

  SUBROUTINE geo2rotll_coord (lon_geo, lat_geo, pollon, pollat, polgam, lon_rot, lat_rot)
    
    IMPLICIT NONE

    REAL(KIND=dp), INTENT(in)       :: lon_geo, lat_geo    ! geogr. lon/lat
    REAL(KIND=dp), INTENT(in)       :: pollon, pollat, polgam
    REAL(KIND=dp), INTENT(out)      :: lon_rot, lat_rot    ! rotated lon/lat

    
    lon_rot = rla2rlarot(lat_geo,lon_geo,pollat,pollon,polgam) ! rotated longitude
    lat_rot = phi2phirot(lat_geo,lon_geo,pollat,pollon)        ! rotated latitude


  END SUBROUTINE geo2rotll_coord

  SUBROUTINE geo2model_coord_domaincheck (idom_model, lon_geo, lat_geo, x_model, y_model, &
                                 is_inside, rlon_min, rlon_max, rlat_min, rlat_max)
    
    ! To be called by all radar PEs AFTER setup_auxgrid_for_cellindex() !!!

    IMPLICIT NONE

    INTEGER, INTENT(in)             :: idom_model
    REAL(KIND=dp), INTENT(in)       :: lon_geo, lat_geo    ! geogr. lon/lat
    REAL(KIND=dp), INTENT(out)      :: x_model, y_model    ! x and y-coord of the model (for COSMO: rotated lon/lat)
                                                           !  (for ICON: this would be just geogr. lon/lat)
    LOGICAL,           INTENT(out)  :: is_inside           ! returns .TRUE. if the point is inside the model domain

    REAL(KIND=dp), INTENT(out), OPTIONAL  :: rlon_min, rlon_max, rlat_min, rlat_max
    INTEGER                         :: mpierror
    LOGICAL                         :: is_inside_loc
    CHARACTER(len=cmaxlen)          :: yzerrmsg

    is_inside_loc = is_inside_loc_domain (idom_model, lon_geo, lat_geo)

    IF (num_radar > 1) THEN
      CALL global_values_radar(is_inside_loc, 'OR', icomm_radar, -1, yzerrmsg, mpierror)
    END IF
    is_inside = is_inside_loc

    x_model = lon_geo
    y_model = lat_geo

    ! Domain corners of computational domain (excluding the nboundlines) in rotated lon/lat coordinates:
    IF (PRESENT(rlon_min)) rlon_min = -999.99_dp
    IF (PRESENT(rlon_max)) rlon_max = -999.99_dp
    IF (PRESENT(rlat_min)) rlat_min = -999.99_dp
    IF (PRESENT(rlat_max)) rlat_max = -999.99_dp

  END SUBROUTINE geo2model_coord_domaincheck

  SUBROUTINE geo2model_cellindex_real (idom_model, lon_geo, lat_geo, i_model, k_model)

    ! To be called AFTER setup_auxgrid_for_cellindex() !!!

    IMPLICIT NONE

    INTEGER, INTENT(in)        :: idom_model
    REAL(KIND=dp), INTENT(in)  :: lon_geo, lat_geo    ! geogr. lon/lat
    REAL(KIND=dp), INTENT(out) :: i_model, k_model    ! horizontal index of model grid cell
                                                            !  which contains lon_geo/lat_geo
    i_model = REAL(geo2cell_index1d_glob (cindex_grid(list_domains_for_radar(idom_model)), p_patch(idom_model), &
         lon_geo, lat_geo), kind=dp)
    k_model = 1.0_dp

  END SUBROUTINE geo2model_cellindex_real

  SUBROUTINE geo2model_cellindex_int (idom_model, lon_geo, lat_geo, i_model, k_model)

    ! To be called AFTER setup_auxgrid_for_cellindex() !!!

    IMPLICIT NONE
    
    INTEGER, INTENT(in)        :: idom_model
    REAL(KIND=dp), INTENT(in)  :: lon_geo, lat_geo    ! geogr. lon/lat
    INTEGER, INTENT(out)       :: i_model, k_model    ! horizontal index (nproma,nblk) of model grid cell
                                                      !  which contains lon_geo/lat_geo
    i_model = geo2cell_index1d_glob (cindex_grid(list_domains_for_radar(idom_model)), p_patch(idom_model), &
         lon_geo, lat_geo)
    k_model = 1

  END SUBROUTINE geo2model_cellindex_int

  SUBROUTINE model_cellindex2geo (idom_model, i_model, k_model, lon_geo, lat_geo)

    IMPLICIT NONE

    INTEGER, INTENT(in)        :: idom_model
    INTEGER, INTENT(in)        :: i_model, k_model  ! horizontal index (nproma,nblk) of model grid cell
    REAL(KIND=dp), INTENT(out) :: lon_geo, lat_geo  ! geogr. lon/lat at the cell center

  
    lon_geo = p_patch(idom_model) % cells % center(i_model,k_model) % lon * raddeg
    lat_geo = p_patch(idom_model) % cells % center(i_model,k_model) % lat * raddeg

  END SUBROUTINE model_cellindex2geo

  FUNCTION is_inside_loc_domain (idom_model, lon_geo, lat_geo) RESULT (is_inside)
    
    ! To be called AFTER setup_auxgrid_for_cellindex() !!!

    IMPLICIT NONE

    INTEGER, INTENT(in)        :: idom_model          ! ICON model domain
    REAL(KIND=dp), INTENT(in)  :: lon_geo, lat_geo    ! geogr. lon/lat
    LOGICAL                    :: is_inside           ! returns .TRUE. if the point is inside the local processor domain

    is_inside = ( geo2cell_index1d (cindex_grid(list_domains_for_radar(idom_model)), p_patch(idom_model), &
                                                                                 lon_geo, lat_geo) > -HUGE(1) )
    
  END FUNCTION is_inside_loc_domain

  SUBROUTINE get_lonlat_domain_center (idom_model, lon_center, lat_center)

    IMPLICIT NONE

    INTEGER, INTENT(in) :: idom_model
    REAL(KIND=dp), INTENT(out) :: lon_center, lat_center  ! geogr. lon/lat of total domain center
    
    ! Is called BEFORE setup_auxgrid_for_cellindex(), so we cannot use the radar-covered domain
    !  to compute it's center. We simply have no idea on the domain center here,
    !  so we specify the middle of Germany as center of ICON's domain 1:

    lon_center = dom_center_lon(list_domains_for_radar(idom_model))
    lat_center = dom_center_lat(list_domains_for_radar(idom_model))
    
  END SUBROUTINE get_lonlat_domain_center

  SUBROUTINE get_rotlatlon_domain_for_superobing (idom_model,pollon_s,pollat_s,polgam_s,startlon_s,startlat_s,endlon_s,endlat_s)

    IMPLICIT NONE

    INTEGER, INTENT(in)        :: idom_model
    REAL(KIND=dp), INTENT(out) :: pollon_s,pollat_s,polgam_s,startlon_s,startlat_s,endlon_s,endlat_s

    INTEGER                    :: idom_r, ierr
    REAL(KIND=dp)              :: cart_width_deg

    ! Define the domain in terms of a rectangular rotated lat/lon grid
    !  for superobing. In case of ICON, this is the radar-covered area:

    idom_r = list_domains_for_radar(idom_model)

    IF (my_radar_id == 0) THEN
      ! The bounding box of the radar-covered area is the grid for superobing:
      pollon_s   = cindex_grid(idom_r) % pollon
      pollat_s   = cindex_grid(idom_r) % pollat
      polgam_s   = cindex_grid(idom_r) % polgam
      startlon_s = cindex_grid(idom_r) % startlon_tot
      startlat_s = cindex_grid(idom_r) % startlat_tot
      endlon_s   = startlon_s + (cindex_grid(idom_r)%nlon_tot-1) * cindex_grid(idom_r)%dlon
      endlat_s   = startlat_s + (cindex_grid(idom_r)%nlat_tot-1) * cindex_grid(idom_r)%dlat
      ! Extend the bounding box of the radar-covered area by +/- one cart_resolution:
      cart_width_deg = supob_cart_resolution/cindex_grid(idom_r)%r_earth*raddeg
      startlon_s = startlon_s - cart_width_deg
      startlat_s = startlat_s - cart_width_deg
      endlon_s   = endlon_s   + cart_width_deg
      endlat_s   = endlat_s   + cart_width_deg
    END IF

    IF (num_radar > 1 .AND. num_radario > 0) THEN
      CALL distribute_values_radar (pollon_s  , 1, 0, icomm_radar, ierr)
      CALL distribute_values_radar (pollat_s  , 1, 0, icomm_radar, ierr)
      CALL distribute_values_radar (polgam_s  , 1, 0, icomm_radar, ierr)
      CALL distribute_values_radar (startlon_s, 1, 0, icomm_radar, ierr)
      CALL distribute_values_radar (startlat_s, 1, 0, icomm_radar, ierr)
      CALL distribute_values_radar (endlon_s  , 1, 0, icomm_radar, ierr)
      CALL distribute_values_radar (endlat_s  , 1, 0, icomm_radar, ierr)
    END IF

  END SUBROUTINE get_rotlatlon_domain_for_superobing

  SUBROUTINE interp2D_model2geo_horiz_scalar (idom_model, field3D_model, k_int, lon_geo, lat_geo, dat2D_geo, found)
    
    IMPLICIT NONE

    !------------------------------------------------------------------------------
    !
    ! Description: Horizontal interpolation of values from model grid (scalar positions) to a
    !              geographic location at a specific height index.
    !
    ! If the lon/lat is not within the local processor domain, a missing value of 
    ! -9999.99 is returned and the flag "found" is set to .FALSE.
    !
    !------------------------------------------------------------------------------
    ! 
    ! Subroutine / Function arguments
    !------------------------------------------------------------------------------
    !
    ! Parameter list:

    INTEGER, INTENT(in)             :: idom_model
    REAL(KIND=dp), INTENT(in)       :: field3D_model(:,:,:) ! ie,je,ke
    INTEGER, INTENT(in)             :: k_int                ! level index of the height layer for horiz. interpolation
    REAL(KIND=dp), INTENT(in)       :: lon_geo, lat_geo     ! interpolation point in geogr. lon/lat
    REAL(KIND=dp), INTENT(out)      :: dat2D_geo            ! interpolated value at location lon/lat
    LOGICAL      , INTENT(out)      :: found

    !------------------------------------------------------------------------------
    ! Local variables:
    INTEGER :: idx, blk, index1d

    dat2D_geo = -9999.99

    ! get index idx, blk of nearest ICON grid cell to the geographic point:
    CALL geo2cell_index2d ( cindex_grid(list_domains_for_radar(idom_model)), p_patch(idom_model), &
         lon_geo, lat_geo, idx, blk )
    found = ( idx > -HUGE(1) )

    IF (found) THEN
      dat2D_geo  = field3D_model(idx,k_int,blk)
    END IF

  END SUBROUTINE interp2D_model2geo_horiz_scalar

  SUBROUTINE interp3D_model2geo_scalar (idom_model, field3D_model, lon_geo, lat_geo, height_msl, dat3D_geo, k, found)
    
    IMPLICIT NONE

    !------------------------------------------------------------------------------
    !
    ! Description: Interpolation of values from model grid (scalar positions) to a
    !              geographic location and height.
    !
    ! If the lon/lat is not within the local model domain or below the surface or
    !  above model top, a missing value of -999.99 is returned and the flag "found"
    !  is set to .FALSE.
    !
    !------------------------------------------------------------------------------
    ! 
    ! Subroutine / Function arguments
    !------------------------------------------------------------------------------
    !
    ! Parameter list:

    INTEGER, INTENT(in)             :: idom_model
    REAL(KIND=dp), INTENT(in)       :: field3D_model(:,:,:)         ! ie,je,ke
    REAL(KIND=dp), INTENT(in)       :: lon_geo, lat_geo, height_msl ! interpolation point in geogr. lon/lat/height
    REAL(KIND=dp), INTENT(out)      :: dat3D_geo                    ! interpolated value at location lon/lat
    INTEGER, INTENT(out)            :: k                            ! index of the model level just above dat3D_geo
    LOGICAL, INTENT(out)            :: found

    !------------------------------------------------------------------------------
    ! Local variables:
    REAL(KIND=dp) :: f1, f2, wk
    INTEGER :: idx, blk, ko


    dat3D_geo = -9999.99

    ! get index idx, blk of nearest ICON grid cell to the geographic point:
    CALL geo2cell_index2d ( cindex_grid(list_domains_for_radar(idom_model)), p_patch(idom_model), &
         lon_geo, lat_geo, idx, blk )

    ! Interpolate the value at location lon_geo/lat_geo:


    IF ( idx > -HUGE(1) ) THEN

      ! The radar station is in this processors sub-domain
      ! and the interpolated surface height of its position
      ! can be distributed to all other nodes:

      ! calculate vertical index k above radar station and vertical interpolation weight wk
      CALL calc_vert_weight(height_msl,hhl(idx,:,blk),hfl(idx,:,blk),je_fwo,k,wk)

      ! If height_msl below the lowest main level but above orography, k = ke and wk = 0.0 is returned.
      ! This means that the vertical interpolation below does a constant extrapolation towards the ground.

      IF (k >= 1) THEN

        ko = MIN(k+1,je_fwo)

        f1 = field3D_model(idx, k,  blk) 
        f2 = field3D_model(idx, ko, blk)

        dat3D_geo = f1*(1.0_dp-wk) + f2*wk
        
        found = .TRUE.

      ELSE

        found = .FALSE.

      END IF

    ELSE

      found = .FALSE.

    END IF


  END SUBROUTINE interp3D_model2geo_scalar

  SUBROUTINE interp_model2radarbins_scalar (dat_model, nobs_above_sfc, &
       ind_intp, w_intp, dat_radarbins)

    IMPLICIT NONE

    !------------------------------------------------------------------------------
    !
    ! Description: Interpolation from model grid (scalar positions) to radar bins, 
    !              based on the indices and weights stored in rs_data-structure.
    !              Allocate and store data in the vector "dat_radarbins".
    !
    !------------------------------------------------------------------------------
    ! 
    ! Subroutine / Function arguments
    !------------------------------------------------------------------------------
    !
    ! Parameter list:

    INTEGER, INTENT(in)             :: nobs_above_sfc    ! from rs_data-structure
    REAL(KIND=dp), INTENT(in)       :: dat_model(:,:,:)  ! ie,je,ke
    INTEGER, INTENT(in) :: ind_intp(:,:)     ! dim = (nobs,2), interp. indices from rs_data-structure
    REAL(KIND=dp), INTENT(in)       :: w_intp(:,:)       ! dim = (nobs,3), interp. weights from rs_data-structure
    REAL(KIND=dp), INTENT(inout)    :: dat_radarbins(:)  ! nobs

    !------------------------------------------------------------------------------
    ! Local variables:

    INTEGER          :: iobs, nk, k, ko, idx, blk
    REAL(KIND=dp)    :: wk,              & ! interpolation weight in k-direction
                        f1, f2

    ! Loop over observation points above SFC:
!CDIR NODEP,VOVERTAKE,VOB
!$omp parallel private (nk,idx,blk,k,wk,ko,f1,f2)
!$omp do 
    DO iobs = 1, nobs_above_sfc

      ! for each radar point determine the continuous number nk of the grid cell
      ! in the ICON column above the grid points
      nk = ind_intp(iobs,1)

      ! calculate vertical index k above radar station and vertical interpolation weight wk
      CALL ind2sub3D(nk, nproma, je_fwo, idx, k, blk)

      ! determine interpolation weights
      wk = w_intp(iobs,1)

      ! interpolate Z (mm^6/m^3) of the 8 surrounding mass points trilinearly to the radar point
      ko = MIN(k+1,je_fwo)

      ! If height of interpolation point is below the lowest main level but above orography, wk = 0.0 and k = ke.
      ! This means that the vertical interpolation below does a constant extrapolation towards the ground.

      f1 = dat_model(idx ,k  ,blk)
      f2 = dat_model(idx ,ko ,blk)

      dat_radarbins(iobs) = f1*(1.0_dp-wk) + f2*wk

    END DO    ! loop over radar points    
!$omp end do
!$omp end parallel

  END SUBROUTINE interp_model2radarbins_scalar

  SUBROUTINE interp2d_model2radarbins_scalar (dat_model, nobs_above_sfc, &
       ind_intp, w_intp, dat_radarbins)

    IMPLICIT NONE

    !------------------------------------------------------------------------------
    !
    ! Description: Interpolation from model grid (scalar positions) to radar bins, 
    !              based on the indices and weights stored in rs_data-structure.
    !              Allocate and store data in the vector "dat_radarbins".
    !              This routine is specifically designed for variables which are
    !              horizontal 2D and provide one value for each vertical column,
    !              i.e., the Tmax_XX_modelgrid fields for the reflectivity calculation.
    !
    !------------------------------------------------------------------------------
    ! 
    ! Subroutine / Function arguments
    !------------------------------------------------------------------------------
    !
    ! Parameter list:

    INTEGER, INTENT(in)             :: nobs_above_sfc    ! from rs_data-structure
    REAL(KIND=dp), INTENT(in)       :: dat_model(:,:)    ! idx,blk
    INTEGER, INTENT(in) :: ind_intp(:,:)     ! dim = (nobs,2), interp. indices from rs_data-structure
    REAL(KIND=dp), INTENT(in)       :: w_intp(:,:)       ! dim = (nobs,3), interp. weights from rs_data-structure
    REAL(KIND=dp), INTENT(inout)    :: dat_radarbins(:)  ! nobs

    !------------------------------------------------------------------------------
    ! Local variables:

    INTEGER          :: iobs, nk, k, idx, blk
    REAL(KIND=dp)    :: f1, f2

    ! Loop over observation points above SFC:
!CDIR NODEP,VOVERTAKE,VOB
!$omp parallel private (nk,idx,blk,k)
!$omp do 
    DO iobs = 1, nobs_above_sfc

      ! for each radar point determine the continuous number nk of the grid cell
      ! in the ICON column above the grid points
      nk = ind_intp(iobs,1)

      ! calculate vertical index k above radar station and vertical interpolation weight wk
      CALL ind2sub3D(nk, nproma, je_fwo, idx, k, blk)

      dat_radarbins(iobs) = dat_model(idx,blk)

    END DO    ! loop over radar points    
!$omp end do
!$omp end parallel

  END SUBROUTINE interp2d_model2radarbins_scalar

!==============================================================================

  SUBROUTINE interp_model2radarbins_vr (u, v, w, nobs_above_sfc, &
       ind_intp, w_intp, u_rp, v_rp, w_rp)

    IMPLICIT NONE

    !------------------------------------------------------------------------------
    !
    ! Description: Interpolation of u,v,w from model grid (u, v on cell centers, w staggered)
    !              first to mass points and then to radar bins, 
    !              based on the indices and weights stored in rs_data-structure.
    !              In ICON, u and v are the geographic wind components and have to be diagnosed
    !              from the cell edge normal components earlier.
    !              w is the vertical wind contribution without the hydrometeor terminal fall speed.
    !              The latter is added later in the operator during the output stage.
    !
    !------------------------------------------------------------------------------
    ! 
    ! Subroutine / Function arguments
    !------------------------------------------------------------------------------
    !
    ! Parameter list:

    INTEGER, INTENT(in)           :: nobs_above_sfc   ! from rs_data-structure
    REAL(KIND=dp), INTENT(in)     :: u(:,:,:), &      ! idx,ke,blk
                                     v(:,:,:), &
                                     w(:,:,:)
    INTEGER, INTENT(in)           :: ind_intp(:,:)  ! dim = (nobs,2), interp. indices from rs_data-structure
    REAL(KIND=dp), INTENT(in)     :: w_intp(:,:)    ! dim = (nobs,1), interp. weights from rs_data-structure
    REAL(KIND=dp), INTENT(inout)  :: u_rp(:), &     ! u at radar points, dim = (at least nobs_above_sfc)
                                     v_rp(:), &     ! v at radar points
                                     w_rp(:)        ! w at radar points

    !------------------------------------------------------------------------------
    ! Local variables:

    INTEGER        :: iobs, nk, idx, blk, k, ku, ko
    REAL(KIND=dp)  :: wk,                      & ! interpolation weight in k-direction
                      uu1, uu2, vv1, vv2, ww1, ww2

    ! Loop over observation points above SFC:
!CDIR NODEP,VOVERTAKE,VOB
!$omp parallel private (nk,idx,blk,k,wk,ko,uu1,uu2,vv1,vv2,ku,ww1,ww2)
!$omp do 
    DO iobs = 1, nobs_above_sfc

      ! for each radar point determine the continuous number nk of the grid cell
      ! in the ICON column above the grid points
      nk = ind_intp(iobs,1)

      ! calculate vertical index k above radar station and vertical interpolation weight wk
      CALL ind2sub3D(nk, nproma, je_fwo, idx, k, blk)

      ! determine interpolation weights
      wk = w_intp(iobs,1)

      ! If height of interpolation point is below the lowest main level but above orography, wk = 0.0 and k = ke.
      ! This means that the vertical interpolation below does a constant extrapolation towards the ground.

      ko = MIN(k+1,je_fwo)
      
      uu1 = u(idx, k,  blk) 
      uu2 = u(idx, ko, blk)

      u_rp(iobs) = uu1*(1.0_dp-wk) + uu2*wk

      vv1 = v(idx, k,  blk) 
      vv2 = v(idx, ko, blk)

      v_rp(iobs) = vv1*(1.0_dp-wk) + vv2*wk

      ko = MIN(k+1,je_fwo+1)
      ku = MIN(k+2,je_fwo+1)

      ww1 = 0.5_dp*(w(idx, ko, blk) + w(idx, ku,blk))
      ww2 = 0.5_dp*(w(idx, ku, blk) + w(idx, k ,blk))

      w_rp(iobs) = ww1*(1.0_dp-wk) + ww2*wk

    END DO    ! loop over radar points    
!$omp end do
!$omp end parallel

  END SUBROUTINE interp_model2radarbins_vr

!==============================================================================

  SUBROUTINE interp_model2azislices_scalar (dat_model, ngrd, &
       ind_intp, w_intp, dat_grd, at_k_upper, at_k_lower)

    IMPLICIT NONE

    !------------------------------------------------------------------------------
    !
    ! Description: Horizontal interpolation from model grid (scalar positions) to
    !              the auxiliary azimutal-slices grid, which is needed for
    !              online beam propagation (ray tracing).
    !              Based on the indices and weights stored in rs_grid-structure.
    !
    !------------------------------------------------------------------------------
    ! 
    ! Subroutine / Function arguments
    !------------------------------------------------------------------------------
    !
    ! Parameter list:

    INTEGER, INTENT(in)             :: ngrd              ! from rs_grid-structure
    REAL(KIND=dp), INTENT(in)       :: dat_model(:,:,:)  ! ie,je,ke
    INTEGER, INTENT(in) :: ind_intp(:,:)     ! dim = (nobs,2), interp. indices from rs_grid-structure
    REAL(KIND=dp), INTENT(in)       :: w_intp(:,:)       ! dim = (nobs,2), interp. weights from rs_grid-structure (dummy)
    REAL(KIND=dp), INTENT(inout)    :: dat_grd(:)        ! dim = (ngrd)

    LOGICAL, OPTIONAL, INTENT(in)   :: at_k_upper, &     ! do the interpolation for one height level up
                                       at_k_lower        ! do the interpolation for one height level down

    !------------------------------------------------------------------------------
    ! Local variables:

    INTEGER :: igrd, nk, idx, blk, ku, k, ko
    LOGICAL                 :: at_k_upper_loc, at_k_lower_loc

    ! Catch optional arguments:
    IF (PRESENT(at_k_upper)) THEN
      at_k_upper_loc = at_k_upper
    ELSE
      at_k_upper_loc = .FALSE.
    END IF

    IF (PRESENT(at_k_lower)) THEN
      at_k_lower_loc = at_k_lower
    ELSE
      at_k_lower_loc = .FALSE.
    END IF

    IF (at_k_upper_loc .AND. at_k_lower_loc) THEN
      WRITE (*,*) 'ERROR in call to interp_model2azislices_scalar: at_k_upper and at_k_lower cannot be both .TRUE.! STOP!'
      STOP
    END IF

    !------------------------------------------------------------------------------
    ! loop over observation points

    IF ( at_k_upper_loc ) THEN

!CDIR NODEP,VOVERTAKE,VOB
!$omp parallel private (nk,idx,blk,ku)
!$omp do 
      DO igrd = 1, ngrd

        ! for each radar point determine the continuous number nk of the grid cell
        ! in the ICON column above the grid points
        nk = ind_intp(igrd,1)
        
        ! calculate vertical index k above radar station and vertical interpolation weight wk
        CALL ind2sub3D(nk, nproma, je_fwo, idx, k, blk)
        
        ku = MAX(1,k-1)  ! one layer above k

        dat_grd(igrd) =  dat_model(idx ,ku, blk )
                
      END DO    ! loop over radar points
!$omp end do
!$omp end parallel

    ELSE IF ( at_k_lower_loc ) THEN

!$omp parallel private (nk,idx,blk,ko)
!$omp do 
      DO igrd = 1, ngrd

        ! for each radar point determine the continuous number nk of the grid cell
        ! in the ICON column above the grid points
        nk = ind_intp(igrd,1)
        
        ! calculate vertical index k above radar station and vertical interpolation weight wk
        CALL ind2sub3D(nk, nproma, je_fwo, idx, k, blk)
        
        ko = MIN(k+1,je_fwo)  ! one layer below k
        
        dat_grd(igrd) =  dat_model(idx ,ko, blk )
                
      END DO    ! loop over radar points
!$omp end do
!$omp end parallel

    ELSE

!$omp parallel private (nk,idx,blk,k)
!$omp do 
      DO igrd = 1, ngrd

        ! for each radar point determine the continuous number nk of the grid cell
        ! in the ICON column above the grid points
        nk = ind_intp(igrd,1)
        
        ! calculate vertical index k above radar station and vertical interpolation weight wk
        CALL ind2sub3D(nk, nproma, je_fwo, idx, k, blk)
        
        dat_grd(igrd) =  dat_model(idx ,k, blk )
                
      END DO    ! loop over radar points
!$omp end do
!$omp end parallel

    END IF

  END SUBROUTINE interp_model2azislices_scalar

!==============================================================================

  SUBROUTINE interp_model2azislices_vr (u, v, w, ngrd, &
       ind_intp, w_intp, u_grd, v_grd, w_grd)

    IMPLICIT NONE

    !------------------------------------------------------------------------------
    !
    ! Description: Mapping of u, v, w from model grid (u/v at scalar, w at staggered positions)
    !              frist to scalar points and then to the auxiliary azimutal-slices grid,
    !              which is needed for online beam propagation (ray tracing).
    !              Based on the indices and weights stored in rs_grid-structure.
    !
    !------------------------------------------------------------------------------
    ! 
    ! Subroutine / Function arguments
    !------------------------------------------------------------------------------
    !
    ! Parameter list:

    INTEGER, INTENT(in)           :: ngrd             ! from rs_grid-structure
    REAL(KIND=dp), INTENT(in)     :: u(:,:,:), &      ! ie,je,ke
                                     v(:,:,:), &
                                     w(:,:,:)
    INTEGER, INTENT(in) :: ind_intp(:,:) ! dim = (nobs,2), interp. indices from rs_grid-structure
    REAL(KIND=dp), INTENT(in)     :: w_intp(:,:)   ! dim = (nobs,2), horiz. interp. weights, here just DUMMYies
    REAL(KIND=dp), INTENT(inout)  :: u_grd(:), &   ! u at aux grid points, dim = (ngrd)
                                     v_grd(:), &   ! v at aux grid points
                                     w_grd(:)      ! w at aux grid points

    !------------------------------------------------------------------------------
    ! Local variables:

    INTEGER         :: igrd, nk, idx, blk, k, ko
    
    ! loop over observation points
       
!CDIR NODEP,VOVERTAKE,VOB
!$omp parallel private (nk,idx,blk,k,ko)
!$omp do 
    DO igrd = 1, ngrd

      ! for each radar point determine the continuous number nk of the grid cell
      ! in the ICON column above the grid points
      nk = ind_intp(igrd,1)
          
      ! calculate vertical index k above radar station and vertical interpolation weight wk
      CALL ind2sub3D(nk, nproma, je_fwo, idx, k, blk)

      u_grd(igrd) = u(idx, k, blk)
      v_grd(igrd) = v(idx, k, blk)

      ko = MIN(k+1,je_fwo+1)

      w_grd(igrd) = 0.5_dp*(w(idx ,ko, blk) + w(idx ,k ,blk ))

    END DO    ! loop over radar points
!$omp end do
!$omp end parallel

  END SUBROUTINE interp_model2azislices_vr

!==============================================================================
!================================================================================


!==============================================================================
!
! Subroutine to pass grib2 model specific header information to EMVORADO's own 
!  dbz-composite output. This will only work correctly if the respective grib2 namelist
!  parameters are specified in the /gribout_nml/ namelist(s) for the actual model domain.
!
! If not, we will assume center "78" (DWD) below so that grib2-writing will not crash
!  but the header informations about experiment ID and ensemble meta data
!  will be most likely wrong.
!
!==============================================================================

#ifdef GRIBAPI
  SUBROUTINE grib2_add_modelspec_info(igribid, idom_model, error)

    USE grib_api,        ONLY : grib_set, GRIB_SUCCESS, grib_get_error_string

    IMPLICIT NONE

    INTEGER, INTENT(in)  :: igribid
    INTEGER, INTENT(in)  :: idom_model
    INTEGER, INTENT(out) :: error

    INTEGER              :: ierr, iprdeftnr, ncenter, nsubcenter
    LOGICAL              :: leps

    CHARACTER(len=*), PARAMETER :: yzroutine = 'grib2_add_modelspec_info'

    error = 0

    leps = (gribout_config(idom_model)%numberOfForecastsInEnsemble > -1)

    IF (gribout_config(idom_model)%generatingCenter == -1) THEN
      ! grib2 configuration informations have not been properly specified
      ! in the /gribout_nml/ namelist(s) for this domain. We set "DWD" as
      ! the provisional generating center, but use the ICON defaults for
      ! the rest of the parameters.
      ncenter = 78
      nsubcenter = 255
      WRITE (*,'(a,/,a,i3)') 'WARNING '//TRIM(yzroutine)//': grib2 header infos on generating model and local use section', &
           'and ensemble members are wrong. Please define these informations in the /gribout_nml/ namelist for domain ', &
           idom_model
    ELSE
      ncenter    = gribout_config(idom_model)%generatingCenter
      nsubcenter = gribout_config(idom_model)%generatingSubcenter
    END IF

    ! local section has to be deleted first before changing the centre 
    CALL grib_set (igribid, 'grib2LocalSectionPresent',        0, ierr)   ! delete local section
    CALL check_grib_err(ierr, 'grib2LocalSectionPresent')
    error = error + ierr
    
    CALL grib_set (igribid, 'centre', ncenter      , ierr)   ! originating centre
    CALL check_grib_err(ierr, 'centre')
    error = error + ierr
    CALL grib_set (igribid, 'subCentre', nsubcenter, ierr)   ! originating subcentre
    CALL check_grib_err(ierr, 'subCentre')
    error = error + ierr


    ! Type of processed data
    ! this is set to "control and perturbed forecast products"
    CALL grib_set (igribid,'typeOfProcessedData', gribout_config(idom_model)%typeOfProcessedData, ierr)
    CALL check_grib_err(ierr,'typeOfProcessedData')
    error = error + ierr

   
    ! Local Use Section
    ! -----------------

    ! For every COSMO application we require a local use section
    CALL grib_set (igribid, 'grib2LocalSectionPresent', 1, ierr)
    CALL check_grib_err(ierr, 'grib2LocalSectionPresent')
    error = error + ierr


    CALL grib_set (igribid, 'localDefinitionNumber', gribout_config(idom_model)%localDefinitionNumber, ierr)
    CALL check_grib_err(ierr, 'localDefinitionNumber')
    error = error + ierr
    CALL grib_set (igribid, 'productionStatusOfProcessedData', gribout_config(idom_model)%productionStatusOfProcessedData, ierr)
    CALL check_grib_err(ierr, 'productionStatusOfProcessedData')
    error = error + ierr
    CALL grib_set (igribid, 'localNumberOfExperiment', gribout_config(idom_model)%localNumberOfExperiment, ierr)
    CALL check_grib_err(ierr, 'localNumberOfExperiment')
    error = error + ierr
    CALL grib_set (igribid,'generatingProcessIdentifier', gribout_config(idom_model)%generatingProcessIdentifier, ierr)
    CALL check_grib_err(ierr,'generatingProcessIdentifier')
    error = error + ierr
    CALL grib_set (igribid,'backgroundProcess', gribout_config(idom_model)%backgroundProcess, ierr)
    CALL check_grib_err(ierr,'backgroundProcess')
    error = error + ierr
    CALL grib_set (igribid,'typeOfGeneratingProcess', gribout_config(idom_model)%typeOfGeneratingProcess, ierr)
    CALL check_grib_err(ierr,'typeOfGeneratingProcess')
    error = error + ierr

    IF (leps) THEN
      iprdeftnr =  1     ! individual ensemble forecast
    ELSE
      iprdeftnr =  0     ! Standard products
    ENDIF
    CALL grib_set (igribid, 'productDefinitionTemplateNumber', iprdeftnr, ierr)
    CALL check_grib_err(ierr, 'productDefinitionTemplateNumber')
    error = error + ierr

    IF (leps) THEN
      CALL grib_set (igribid, 'localTypeOfEnsembleForecast', gribout_config(idom_model)%localTypeOfEnsembleForecast, ierr)
      CALL check_grib_err(ierr, 'localTypeOfEnsembleForecast')
      error = error + ierr
      CALL grib_set (igribid,'typeOfEnsembleForecast', gribout_config(idom_model)%typeOfEnsembleForecast, ierr)
      CALL check_grib_err(ierr,'typeOfEnsembleForecast')
      error = error + ierr
      CALL grib_set (igribid,'perturbationNumber', gribout_config(idom_model)%perturbationNumber, ierr) ! Ensemble Forecast
      CALL check_grib_err(ierr,'perturbationNumber')
      error = error + ierr
      CALL grib_set (igribid,'numberOfForecastsInEnsemble', gribout_config(idom_model)%numberOfForecastsInEnsemble, ierr) ! Ensemble Forecast
      CALL check_grib_err(ierr,'numberOfForecastsInEnsemble')
      error = error + ierr
    END IF

  CONTAINS

    SUBROUTINE check_grib_err (griberr, cident)

      INTEGER         , INTENT(in) :: griberr
      CHARACTER(len=*), INTENT(in) :: cident

      CHARACTER(len=250)           :: griberrmsg

      IF (griberr /= GRIB_SUCCESS) THEN
        griberrmsg(:) = ' '
        CALL grib_get_error_string(griberr, griberrmsg)
        WRITE (*,'(a)') 'ERROR grib_set('//TRIM(cident)//') in '//TRIM(yzroutine)//': '//TRIM(griberrmsg)
      END IF

    END SUBROUTINE check_grib_err

  END SUBROUTINE grib2_add_modelspec_info
#endif

!================================================================================
!================================================================================

  !======================================================================================
  !
  ! Trigger routine for automatic warm bubbles (ldo_bubbles = .TRUE.)
  !
  ! Pass the information about missing cells, which have been detected previously by
  !  SR detect_missing_cells() from radar_bubblegen.f90, to the ICON model
  !  by passing the relevant informations (shape, amplitude, location, time)
  !  to respective data vectors used in ICON.
  !
  ! This routine handles both synchroneous and asynchroneous radar IO:
  !
  ! - in the synchroneous case, it is called right after detect_missing_cells()
  !   in organize_radar().
  ! - in the asynchroneous case, it is called at a later model time on the
  !   the workers (determined by namelist parameter "t_offset_bubble_trigger_async"),
  !   so that the workers can continue with model integration while the
  !   asynchroneous IO PEs construct the composites and detect the missing cells.
  !
  ! In any case, this routine must be called on all PEs belonging to the
  !  icomm_radar_dom(idom_model) group of PEs, i.e., if lradar_pe_dom(idom_in) = .TRUE.
  !
  !======================================================================================


  SUBROUTINE trigger_warm_bubbles (idom_model, time_mod_sec, dt_advect, zlow_meanwind_for_advect, &
                                   zup_meanwind_for_advect)

    INTEGER, intent(in)        :: idom_model   ! No. of the model domain in the hosting model
    REAL (kind=dp), INTENT(in) :: time_mod_sec ! model time in seconds since model start

    REAL (kind=dp), INTENT(in) :: dt_advect                 ! Time scale for downstream advection of automatic bubbles [seconds]
    REAL (kind=dp), INTENT(in) :: zlow_meanwind_for_advect  ! The lower bound of averaging height interval for bubble advection speed [meters AMSL]
    REAL (kind=dp), INTENT(in) :: zup_meanwind_for_advect   ! The upper bound of averaging height interval for bubble advection speed [meters AMSL]

    LOGICAL, SAVE :: first_call = .TRUE.

    INTEGER       :: zntstep_bub, i, k, ierror
    LOGICAL       :: is_inside
    REAL(kind=dp) :: rlon_model, rlat_model
    CHARACTER(len=cmaxlen) :: yerrmsg
    CHARACTER(len=32) :: yzroutine

    yzroutine(:) = ' '
    yzroutine    = 'trigger_warm_bubbles'
    ierror = 0

    IF (ldebug_radsim) WRITE (*,*) TRIM(yzroutine), ' on proc ', my_radar_id
    
    !----------------------------------------------------------------------------------------
    ! .. In the first call to this routine, gather common constant parameters for 
    !    automatic bubbles from NAMELIST parameters in /RADARSIM_PARAMS/ and store 
    !    them in the respectivedata structures of src_artifdata.f90 for later triggering:

    IF (first_call) THEN
      first_call = .FALSE.

      ctype_tempdist(:)      =  bubble_type             ! Type of perturbation 'cos-hrd', 'cos-instant'
      bub_centz(:)           =  bubble_centz            ! Center (Z) of temperature disturbances [m]
      bub_radx(:)            =  bubble_radx             ! Length scale (X) of temperature disturbances [m]
      bub_rady(:)            =  bubble_rady             ! Length scale (Y) of temperature disturbances [m]
      bub_radz(:)            =  bubble_radz             ! Length scale (Z) of temperature disturbances [m]
      bub_rotangle(:)        =  bubble_rotangle         ! Rotation angle of main axes of temperature disturbances [degrees]
      bub_heatingrate(:)     =  bubble_heatingrate      ! Constant heating rate [K*s**-1] for cos-hrd bubble
      bub_dT(:)              =  bubble_dT               ! Temperature increment [K] for cos bubble
      ladd_bubblenoise_t(:)  =  bubble_addnoise_T       ! Switch(es) (up to 50) to overlay temperature disturbances (LOGICAL)
      bub_dT_bubblenoise(:)  =  bubble_dT_noise         ! In case of ladd_bubblenoise_t=.true., relative noise level, such that
      ! dT_bubble = dT_bubble * (1 + bub_dT_bubblenoise * random_noise[-1,1])
      lbub_rhconst(:)        =  bubble_holdrhconst      ! whether or not to keep the relative humidity constant during heating
      IF (TRIM(bubble_type) == 'cos-instant') THEN
        bub_timespan(:)      =  1
      ELSE
        bub_timespan(:)      =  CEILING(bubble_timespan/dtime) ! timespan of the heating in number of time steps
      END IF

    END IF

    !----------------------------------------------------------------------------------------
    ! .. Distribute the number of bubbles in the bubble_list to all PEs, to
    !    see if there have been any new bubbles to trigger.
    !
    !*********** THIS IS THE SYNCHRONIZATION POINT WITH THE WORKER PEs! *********************

    IF (num_radar > 1) THEN
      CALL distribute_values_radar(bubble_list%nbubbles, 1, radario_master_dom(idom_model), &
                                   icomm_radar_dom(idom_model), ierror)
    END IF

    !----------------------------------------------------------------------------------------
    ! .. If detect_missing_cells() was active recently, there might be missing cells in the bubble_list.
    !    In this case, transfer them to the respective data structures of src_artifdata.f90 for later triggering:

    IF (bubble_list%nbubbles > 0) THEN
     
      !----------------------------------------------------------------------------------------
      ! .. Initialization of main switch list for bubbles, to deactivate all "old" automatic bubbles:

      ltempdist(1:nautobubbles_max) = .FALSE.

      !----------------------------------------------------------------------------------------
      ! .. If there have been new missing cells detected, copy as much as possible of them to
      !    the ICON namelist vectors to trigger bubbles in the next timestep:

      IF (num_radar > 1) THEN
        CALL distribute_values_radar( bubble_list%nbubbles_reject, 1,             radario_master_dom(idom_model), &
                                      icomm_radar_dom(idom_model), ierror)
        CALL distribute_values_radar( bubble_list%bub_centlon,  nautobubbles_max, radario_master_dom(idom_model), &
                                      icomm_radar_dom(idom_model), ierror)
        CALL distribute_values_radar( bubble_list%bub_centlat,  nautobubbles_max, radario_master_dom(idom_model), &
                                      icomm_radar_dom(idom_model), ierror)
        CALL distribute_values_radar( bubble_list%bub_timestamp,nautobubbles_max, radario_master_dom(idom_model), &
                                      icomm_radar_dom(idom_model), ierror)
      END IF

      IF (lcompute_pe_fwo) THEN
      
        !----------------------------------------------------------------------------------------
        ! .. Advect bubbles (update-in-place and redistribute among icomm_radar) by calling the model specific advection routine:
        !    Takes into account also the time difference between the actual model time and the timestamp of the bubble,
        !    so that this routine can also be used for the asynchroneous IO mode, when the call to trigger_warm_bubbles
        !    is actually at the next radar output time step, not the actual time step:

        CALL advect_bubbles ( idom_model, time_mod_sec, bubble_list, dt_advect, &
                              zlow_meanwind_for_advect, zup_meanwind_for_advect)
       
        ! .. Define the new bubble positions and times from the list of newly detected missing cells:
        k = 0
        DO i=1, bubble_list%nbubbles

          ! Check if bubble is inside modeldomain:
          CALL geo2model_coord_domaincheck(idom_model, bubble_list%bub_centlon(i), bubble_list%bub_centlat(i), &
                                           rlon_model, rlat_model, is_inside)
          IF (is_inside) THEN
            k = k + 1
            bub_centlon(k)     = bubble_list%bub_centlon(i)
            bub_centlat(k)     = bubble_list%bub_centlat(i)
            zntstep_bub        = NINT( time_mod_sec / dtime )
            ntstep_bubble(k)   = zntstep_bub + 1         ! This really denotes the bubble time
            htempdist(k)       = (ntstep_bubble(k) * dtime) / 3600.0_dp ! This is "decoration"
            ltempdist(k)       = .TRUE.                  ! This activates the bubble
            bub_zeitzaehler(k) = 0                       ! This resets the time counter for the bubble

            IF (my_cart_id_fwo == 0) THEN
              WRITE(*,'(a,f0.4,a,f0.4,a,f0.1,a,f0.1,a)') 'INFO '//TRIM(yzroutine)//': Bubble of type '//TRIM(bubble_type)//  &
                   ' triggered at lon = ', bubble_list%bub_centlon(i), ' and lat = ', bubble_list%bub_centlat(i), &
                   ' with timestamp = ', bubble_list%bub_timestamp(i), &
                   ' s triggered at t = ', time_mod_sec, ' s'
            END IF
          ELSE
            
            bubble_list%nbubbles_reject = bubble_list%nbubbles_reject + 1
            IF (my_cart_id_fwo == 0) THEN
              WRITE(*,'(a,f0.4,a,f0.4,a,f0.1,a)') 'WARNING'//TRIM(yzroutine)//': Bubble location after advection '// &
                   ' outside of model domain at lon = ', bubble_list%bub_centlon(i), ' and lat = ', bubble_list%bub_centlat(i), &
                   ' with timestamp = ', bubble_list%bub_timestamp(i), ' s'
            END IF
          END IF
          
        END DO
        bubble_list%nbubbles = k

        IF (my_cart_id_fwo == 0) THEN
          WRITE(*,'(a,i3,a,f0.1,a)') 'INFO SUMMARY '//TRIM(yzroutine)//': ', &
               bubble_list%nbubbles, ' automatic bubble(s) triggered at time = ', time_mod_sec, ' s'
        END IF
        
      END IF

      !----------------------------------------------------------------------------------------
      ! .. This round of bubbles has been enforced in the model by setting their respective
      !    model parameters. The bubble_list now has to be reset to neutral values, until the
      !    next detection round in detect_missing_cells() from radar_bubblegen.f90 gives new bubbles:

      bubble_list%nbubbles         = 0
      bubble_list%nbubbles_reject  = 0
      bubble_list%bub_centlon(:)   = -HUGE(1.0_dp)
      bubble_list%bub_centlat(:)   = -HUGE(1.0_dp)
      bubble_list%bub_timestamp(:) = -HUGE(1.0_dp)

    END IF

    IF (ldebug_radsim) WRITE (*,*)  'done with '//TRIM(yzroutine)//' on proc ', my_radar_id

  END SUBROUTINE trigger_warm_bubbles
  
  !================================================================================

  !================================================================================
  !
  ! Advect bubbles downstream for a given advection time step dt_advect plus
  !  a time difference which is given by the actual model time time_mod_sec
  !  and the bubble timestamp.
  !
  ! Bubble parameters in the INOUT list bubble_list are updated in place.
  !
  ! If advected positions are outside model domain or too close to the boundaries,
  !  bubbles are removed from the list.
  !
  !================================================================================
  
  SUBROUTINE advect_bubbles ( idom_model, time_mod_sec, bubble_list, dt_advect, &
                              zlow_meanwind_for_advect, zup_meanwind_for_advect)

    INTEGER, intent(in)        :: idom_model   ! No. of the model domain in the hosting model (dummy for COSMO)
    REAL (kind=dp), INTENT(in) :: time_mod_sec ! model time in seconds since model start
    TYPE(bubble_list_type), INTENT(inout) :: bubble_list
    
    REAL (kind=dp), INTENT(in) :: dt_advect                 ! Time scale for downstream advection of automatic bubbles [seconds]
    REAL (kind=dp), INTENT(in) :: zlow_meanwind_for_advect  ! The lower bound of averaging height interval for bubble advection speed [meters AMSL]
    REAL (kind=dp), INTENT(in) :: zup_meanwind_for_advect   ! The upper bound of averaging height interval for bubble advection speed [meters AMSL]

    
    INTEGER                    :: i, k, ierror
    REAL (kind=dp), DIMENSION(nautobubbles_max)  :: &
         lon_model, lat_model, bub_centlon_orig, bub_centlat_orig
    REAL (kind=dp)             :: dt_advect_total, dummylon, dummylat, dist_x, dist_y
    REAL (kind=dp)             :: u_bub, v_bub, dz_sum, hhl_u, hhl_o, zdz, zml, utmp, vtmp
    LOGICAL                    :: is_inside, found
    CHARACTER(LEN=50)          :: yzroutine
    CHARACTER(len=cmaxlen)     :: yerrmsg

    
    yzroutine(:) = ' '
    yzroutine    = 'advect_bubbles'
    yerrmsg(:)   = ' '

    IF (ldebug_radsim) WRITE (*,*) TRIM(yzroutine), ' on proc ', my_radar_id

    IF (lcompute_pe_fwo) THEN

      IF (bubble_list%nbubbles > 0) THEN

        DO i=1, bubble_list%nbubbles
          
          IF (ABS(dt_advect) >= 1e-4_wp .OR. time_mod_sec > bubble_list%bub_timestamp(i)) THEN

            dt_advect_total = dt_advect + MAX(time_mod_sec-bubble_list%bub_timestamp(i), 0.0_dp)

            bub_centlon_orig(i) = bubble_list%bub_centlon(i)
            bub_centlat_orig(i) = bubble_list%bub_centlat(i)
             
            IF (ldebug_radsim) THEN
              IF (time_mod_sec > bubble_list%bub_timestamp(i) .AND. my_cart_id_fwo == 0) THEN
                WRITE (*,'(a,f0.4,a,f0.4,a,f0.1,a,f0.1,a)') &
                     'INFO '//TRIM(yzroutine)//': automatic bubble at position lon = ', &
                     bub_centlon_orig(i),' lat = ',bub_centlat_orig(i), &
                     'from time ', bubble_list%bub_timestamp(i),' s will be advected to the actual model time ', &
                     time_mod_sec,' s!'
              END IF
            END IF
            
            IF ( is_inside_loc_domain (idom_model, bub_centlon_orig(i), bub_centlat_orig(i)) ) THEN

              ! .. Bubble has been found on this PE, so compute advection

              ! .. Advection parameters:
              !     Average horizontal wind components from 3000 to 6000 m, integral average using the trapezoidal rule:
              u_bub  = 0.0_wp
              v_bub  = 0.0_wp
              dz_sum = 0.0_wp
              DO k = je_fwo, 1, -1

                CALL interp2D_model2geo_horiz_scalar(idom_model, hhl, k  , bub_centlon_orig(i), bub_centlat_orig(i), hhl_o, found)
                CALL interp2D_model2geo_horiz_scalar(idom_model, hhl, k+1, bub_centlon_orig(i), bub_centlat_orig(i), hhl_u, found)

                zml = 0.5_wp * (hhl_u + hhl_o)
                IF (zml >= zlow_meanwind_for_advect .AND. zml <= zup_meanwind_for_advect) THEN
                  zdz =   MIN(hhl_o  , zup_meanwind_for_advect ) - MAX(hhl_u, zlow_meanwind_for_advect)
                  CALL interp2D_model2geo_horiz_scalar(idom_model, u, k, bub_centlon_orig(i), bub_centlat_orig(i), utmp, found)
                  CALL interp2D_model2geo_horiz_scalar(idom_model, v, k, bub_centlon_orig(i), bub_centlat_orig(i), vtmp, found)
                  u_bub = u_bub + utmp
                  v_bub = v_bub + vtmp
                  dz_sum = dz_sum + zdz
                END IF
              END DO
              IF (dz_sum > 0.0_wp) THEN
                u_bub = u_bub / dz_sum
                v_bub = v_bub / dz_sum
              END IF

              ! .. Advection step:
              dist_x = u_bub*dt_advect_total
              dist_y = v_bub*dt_advect_total
              ! .. Compute the geographic location which is located at an east-west-distance of dist_x [m]
              !    and a south-north-distance dist_y [m] from the original bubble location.
              !    This computation assumes that x and y are the pseudo-kartesian coordinates of an azimutal
              !    equidistant map projection centered around the original bubble location:
              CALL polar2geo_xy (bub_centlon_orig(i), bub_centlat_orig(i), r_earth_dp, dist_x, dist_y, lon_model(i), lat_model(i))

            ELSE

              ! .. Bubble is not on this PE, so set indices to a dummy value which will not "win" in the
              !     global_values_radar('MAX') calls below:

              lon_model(i) = -HUGE(1.0_dp)
              lat_model(i) = -HUGE(1.0_dp)

            END IF
            
          END IF

        END DO
         
        ! Distribute the corrected positions to all compute PEs. We use a global MAX operation, because
        !  only on the PE(s) where the bubble was found the indices are > -999.99:
        IF (num_compute_fwo > 1) THEN
          CALL global_values_radar(lon_model, bubble_list%nbubbles, 'MAX', icomm_cart_fwo, -1, yerrmsg, ierror)
          CALL global_values_radar(lat_model, bubble_list%nbubbles, 'MAX', icomm_cart_fwo, -1, yerrmsg, ierror)
        END IF

        ! Compute and check new positions:
        k = 0
        DO i = 1, bubble_list%nbubbles
          IF (my_cart_id_fwo == 0) THEN
            IF (lon_model(i) < -0.9_dp*HUGE(1.0_dp) .OR. lat_model(i) < -0.9_dp*HUGE(1.0_dp)) THEN
              WRITE (*,'(a,f0.4,a,f0.4,a)') 'WARNING '//TRIM(yzroutine)//': Automatic bubble position (lon = ', &
                   bub_centlon_orig(i),' lat = ',bub_centlat_orig(i), &
                   ' ) not found on any PE during downstream advection! Should not have happened!'
            END IF
          END IF
          ! Remove advected bubbles which are too close to the domain boundaries:
!!$          IF ( WHAT IS THE CORRECT CHECK ???) THEN
!!$            IF (my_cart_id_fwo == 0) THEN
!!$              WRITE (*,'(a,f0.4,a,f0.4,a)') 'WARNING '//TRIM(yzroutine)// &
!!$                   ': Advected automatic bubble too close to boundary, removed ( lon = ', &
!!$                   lat_model(i), ' lat = ', lat_model(i), ' )!'
!!$            END IF
!!$            lon_model(i) = -HUGE(1.0_dp)
!!$            lat_model(i) = -HUGE(1.0_dp)
!!$          END IF
          ! Remove advected bubbles which are outside model domain:
          CALL geo2model_coord_domaincheck(idom_model, lon_model(i), lat_model(i), &
                                           dummylon, dummylat, is_inside)
          IF ( .NOT. is_inside ) THEN
            IF (my_cart_id_fwo == 0) THEN
              WRITE (*,'(a,f0.4,a,f0.4,a)') 'WARNING '//TRIM(yzroutine)// &
                   ': Advected automatic bubble outside model domain, removed ( lon = ', &
                   lon_model(i), ' lat = ', lat_model(i), ' )!'
            END IF
            lon_model(i) = -HUGE(1.0_dp)
            lat_model(i) = -HUGE(1.0_dp)
          END IF
          IF (lon_model(i) >= -0.9_dp*HUGE(1.0_dp) .AND. lat_model(i) >= -0.9_dp*HUGE(1.0_dp)) THEN
            ! Advected bubble is inside model domain, so update position:
            k = k + 1
            bubble_list%bub_centlon(k) = lon_model(i)
            bubble_list%bub_centlat(k) = lat_model(i)
            IF (my_cart_id_fwo == 0 .AND. ldebug_radsim) THEN
              WRITE (*,'(a,f0.4,a,f0.4,a,f0.4,a,f0.4,a,f0.1,a)') 'INFO '//TRIM(yzroutine)// &
                   ': bubble advected from position ( lon = ',bub_centlon_orig(i),' lat = ',bub_centlat_orig(i), &
                   ' ) to ( lon = ', bubble_list%bub_centlon(k), &
                          ' lat = ', bubble_list%bub_centlat(k),' ) at step = ',time_mod_sec, ' s'
            END IF
          END IF
        END DO

        ! Clean the bubble list for the rejected bubbles after advection:
        bubble_list%bub_centlon(k+1:bubble_list%nbubbles)   = -HUGE(1.0_dp)
        bubble_list%bub_centlat(k+1:bubble_list%nbubbles)   = -HUGE(1.0_dp)
        bubble_list%bub_timestamp(k+1:bubble_list%nbubbles) = -HUGE(1.0_dp)
        bubble_list%nbubbles_reject = bubble_list%nbubbles_reject + bubble_list%nbubbles - k
        bubble_list%nbubbles        = k
        
     END IF  ! nbubbles > 0

    END IF  ! lcompute_pe_fwo

    IF (ldebug_radsim) WRITE (*,*)  'done with '//TRIM(yzroutine)//' on proc ', my_radar_id

  END SUBROUTINE advect_bubbles

!================================================================================

#endif               /* HAVE_RADARFWO */

END MODULE radar_interface
