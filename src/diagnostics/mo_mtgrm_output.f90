!>
!! Data structures and subroutines for meteogram output.
!!
!! The sampling intervals for meteogram data are independent 
!! from global output steps. Values are buffered in memory until
!! the next field output is invoked.
!! Before each write operation, data is gathered from all working
!! PEs by the IO PE and written to a NetCDF file.
!!
!!
!! How to add new variables for sampling:
!! --------------------------------------
!! In SR "meteogram_setup_variables", insert
!!
!! a) for volume variables:
!!      CALL add_atmo_var(VAR_GROUP_ATMO_ML, "myvarname", "myvarunit", "long name", jg, <state_var>)
!!      where "VAR_GROUP_ATMO_ML" (or "VAR_GROUP_ATMO_HL" or "VAR_GROUP_SOIL_HL", ...)
!!      determine the level heights of this variable.
!!
!! b) for surface variables:
!!      CALL add_sfc_var(VAR_GROUP_SFC, "myvarname", "myvarunit", "long name", jg, <state_var>)
!!
!!
!! Known limitations:
!! ------------------
!! - So far, only NetCDF file format is supported, which is written by
!!   IO PE #0.
!! - ASCII output data (similar to COSMO) must be generated in a
!!   post-processing step.
!! - In case of an application crash, latest meteogram data, which has
!!   not yet been written to hard disk, may be lost.
!!
!! @author F. Prill, DWD
!!
!! @par Revision History
!! Initial implementation  by  F. Prill, DWD (2011-08-22)
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
!!
!! TODO[FP] : use the same GNAT data structure as for the RBF
!!            coefficient computation!
!! TODO[FP] : extend list of meteogram variables.
!! TODO[FP] : use cdi functionality instead of direct NetCDF access.
!! TODO[FP] : MPI communication of height levels and header info is 
!!            necessary only once at the beginning.

MODULE mo_meteogram_output

  USE mo_kind,                  ONLY: wp
  USE mo_datetime,              ONLY: t_datetime, iso8601
  USE mo_exception,             ONLY: message, message_text, finish
  USE mo_mpi,                   ONLY: p_n_work, my_process_is_stdio,      &
    &                                 get_my_mpi_all_id, p_wait,          &
    &                                 p_send_packed, p_irecv_packed,      &
    &                                 p_pack_int,    p_pack_real,         &
    &                                 p_pack_int_1d, p_pack_real_1d,      &
    &                                 p_pack_real_2d,                     &
    &                                 p_unpack_int,    p_unpack_real,     &
    &                                 p_unpack_int_1d, p_unpack_real_1d,  &
    &                                 p_unpack_real_2d,                   &
    &                                 get_my_mpi_all_id,                  &
    &                                 get_mpi_all_workroot_id,            &
    &                                 p_real_dp_byte,                     &
    &                                 MPI_ANY_SOURCE
  USE mo_model_domain,          ONLY: t_patch
  USE mo_parallel_config,       ONLY: nproma
  USE mo_impl_constants,        ONLY: inwp, max_dom, SUCCESS, zml_soil
  USE mo_communication,         ONLY: idx_1d, blk_no, idx_no
  USE mo_ext_data,              ONLY: t_external_data
  USE mo_nonhydro_state,        ONLY: t_nh_state, t_nh_prog, t_nh_diag,   &
    &                                 t_nh_metrics
  USE mo_nwp_phy_state,         ONLY: t_nwp_phy_diag
  USE mo_nwp_lnd_state,         ONLY: t_lnd_state, t_lnd_prog, t_lnd_diag
  USE mo_cf_convention,         ONLY: t_cf_var, t_cf_global
  ! TODO[FP] : When using an already built GNAT, not all of the
  ! following USEs will be necessary:
  USE mo_gnat_gridsearch,       ONLY: gnat_init_grid, gnat_destroy, gnat_tree,&
    &                                 gnat_query_containing_triangles,        &
    &                                 gnat_merge_distributed_queries, gk
  USE mo_dynamics_config,       ONLY: nnow
  USE mo_io_config,             ONLY: lwrite_extra, inextra_2d, inextra_3d
  USE mo_run_config,            ONLY: iqv, iqc, iqi
  USE mo_util_string,           ONLY: int2string
  USE mo_meteogram_config,      ONLY: t_meteogram_output_config, t_station_list, &
    &                                 FTYPE_NETCDF, MAX_NAME_LENGTH, MAX_NUM_STATIONS
  
  IMPLICIT NONE
  
  PRIVATE
  CHARACTER(LEN=*), PARAMETER :: modname     = 'mo_meteogram_output'
  INTEGER,          PARAMETER :: dbg_level   = 0

  INCLUDE 'netcdf.inc'

  ! IO routines.
  ! called collectively, though non-IO PEs are occupied
  ! only for the case of distributed write mode.
  PUBLIC ::  meteogram_init  
  PUBLIC ::  meteogram_is_sample_step
  PUBLIC ::  meteogram_sample_vars
  PUBLIC ::  meteogram_finalize  
  PUBLIC ::  meteogram_flush_file

  INTEGER, PARAMETER :: MAX_TIME_STAMPS      =  100  !< max. number of time stamps
  INTEGER, PARAMETER :: MAX_NVARS            =   30  !< max. number of sampled 3d vars
  INTEGER, PARAMETER :: MAX_NSFCVARS         =   10  !< max. number of sampled surface vars
  INTEGER, PARAMETER :: MAX_DESCR_LENGTH     =  128  !< length of info strings (see cf_convention)

  ! arbitrarily chosen value for buffer size (somewhat large for safety reasons)
  INTEGER, PARAMETER :: MAX_HEADER_SIZE      =  128  !< *p_real_dp_byte

  !! Groups of variables; using this we can distinguish different height axes
  INTEGER, PARAMETER :: VAR_GROUP_ATMO_ML    =    1  !< variables defined on model levels
  INTEGER, PARAMETER :: VAR_GROUP_ATMO_HL    =    2  !< variables defined on half levels
  INTEGER, PARAMETER :: VAR_GROUP_SURFACE    =    3  !< surface variables
  INTEGER, PARAMETER :: VAR_GROUP_SOIL_ML    =    4  !< variables defined on soil half levels
  INTEGER, PARAMETER :: VAR_GROUP_SOIL_MLp2  =    5  !< height levels [0m, soil half levels, -14.58m]

  !>
  !! Generic interface for adding atmospheric vars to list (required
  !! to cope with 4d vars, e.g. with tiles):
  INTERFACE add_atmo_var
    module procedure add_atmo_var_3d
    module procedure add_atmo_var_4d
  END INTERFACE

  !>
  !! Generic interface for adding surface vars to list (required
  !! to cope with 3d vars, e.g. with tiles):
  INTERFACE add_sfc_var
    module procedure add_sfc_var_2d
    module procedure add_sfc_var_3d
  END INTERFACE
  
  !>
  !! Storage for information on a single variable
  !!
  TYPE t_var_info
    TYPE(t_cf_var)        :: cf                   !< variable name, unit
    INTEGER               :: igroup_id            !< variable group (surface vars, soil temperatures, ...)
    INTEGER               :: nlevs                !< number of levels for this variable
    INTEGER,  POINTER     :: levels(:)            !< level indices (1:nlevs)
    REAL(wp), POINTER     :: p_source_3d(:,:,:)   !< pointer to source array  (nproma, nlev, nblk)
  END TYPE t_var_info

  !>
  !! Storage for information on a single surface variable
  !!
  TYPE t_sfc_var_info
    TYPE(t_cf_var)                 :: cf              !< variable name, unit
    INTEGER                        :: igroup_id       !< variable group (surface vars, soil temperatures, ...)
    REAL(wp), POINTER              :: p_source(:,:)   !< pointer to source array
  END TYPE t_sfc_var_info

  !>
  !! Value buffer for a single variable of a station.
  !!
  TYPE t_var_buffer
    REAL(wp), POINTER              :: heights(:)     !< level heights
    REAL(wp), POINTER              :: values(:,:)    !< sampled data for different levels (1:nlevs,time)
  END TYPE t_var_buffer

  !>
  !! Value buffer for a single surface variable of a station.
  !!
  TYPE t_sfc_var_buffer
    REAL(wp), POINTER              :: values(:)      !< sampled data (1:time)
  END TYPE t_sfc_var_buffer

  !>
  !! Data structure containing time slice info.
  !!
  TYPE t_time_stamp
    INTEGER           :: istep !< iteration step of model
    CHARACTER(len=16) :: zdate !< date and time of point sample (iso8601)
  END TYPE t_time_stamp

  !>
  !! Data structure containing meteogram data and meta info for a
  !! single station.
  !!
  !! Apart from header info, data structures of this type buffer point
  !! values for a meteogram between file I/O.
  !! The time slices and variables where the values are collected are
  !! defined outside of this data structure in a record of type
  !! t_meteogram_data.
  !!
  !! Note: This info is different for different patches.
  !!
  TYPE t_meteogram_station
    ! Meteogram header (information on location, ...)
    INTEGER                       :: station_idx(2)   !< (idx,block) of station specification
    INTEGER                       :: tri_idx(2)       !< triangle index (global idx,block)
    INTEGER                       :: tri_idx_local(2) !< triangle index (idx,block)
    INTEGER                       :: owner            !< proc ID where station is located.    
    REAL(wp)                      :: hsurf            !< surface height
    REAL(wp)                      :: frland           !< fraction of land
    REAL(wp)                      :: fc               !< Coriolis parameter
    INTEGER                       :: soiltype         !< soil type

    ! Buffer for currently stored meteogram values.
    TYPE(t_var_buffer),     POINTER :: var(:)         !< sampled data (1:nvars)
    TYPE(t_sfc_var_buffer), POINTER :: sfc_var(:)     !< sampled data (1:nsfcvars)
  END TYPE t_meteogram_station

  !>
  !! Storage for information on the set of collected variables for
  !! several stations.
  !!
  TYPE t_meteogram_data
    ! variable info:
    INTEGER                        :: nvars, nsfcvars  !< number of sampled variables and surface variables
    INTEGER                        :: max_nlevs        !< maximum no. of levels for variables
    TYPE(t_var_info),     POINTER  :: var_info(:)      !< info for each variable (1:nvars)
    TYPE(t_sfc_var_info), POINTER  :: sfc_var_info(:)  !< info for each surface variable (1:nsfcvars)
    ! time stamp info:
    INTEGER                        :: icurrent         !< current time stamp index
    TYPE(t_time_stamp), POINTER    :: time_stamp(:)    !< info on sample times
    ! value buffers:
    TYPE(t_meteogram_station), POINTER :: station(:,:) !< meteogram data and meta info for each station (idx,blk).
    INTEGER                        :: nstations, nblks, npromz
  END TYPE t_meteogram_data
  
  !>
  !! Data structure specifying output file for meteogram data.
  !!
  TYPE t_meteogram_file
    INTEGER                        :: ftype   !< file type (NetCDF, ...)
    CHARACTER(len=MAX_NAME_LENGTH) :: zname   !< file name string
    INTEGER                        :: file_id !< meteogram file ID
    TYPE(t_cf_global)              :: cf      !< meta info
  END TYPE t_meteogram_file

  !>
  !! Data structure specifying NetCDF IDs
  !!
  TYPE t_ncid
    INTEGER  :: nstations, nvars, charid, station_name, station_lat, station_lon,         &
      &         station_idx, station_blk, station_hsurf, station_frland, station_fc,      &
      &         station_soiltype, nsfcvars, var_name, var_unit, sfcvar_name, sfcvar_unit, &
      &         var_group_id, sfcvar_group_id, var_nlevs, max_nlevs, var_levels, timeid,  &
      &         time_step, dateid, var_values, sfcvar_values, var_heights, var_longname,  &
      &         sfcvar_longname
  END TYPE t_ncid

  !>
  !! Data structure containing internal indices for variables
  !!
  TYPE t_var
    INTEGER :: no_atmo_vars       !< number of atmo variables declared so far
    INTEGER :: no_sfc_vars        !< number of surface variables declared so far
  END TYPE t_var

  ! -------------------------------------------------------------------------------------------

  !! -- module data: --
  TYPE(t_meteogram_data), SAVE, TARGET  :: meteogram_local_data(1:max_dom)  !< meteogram data local to this PE
  TYPE(t_meteogram_data), SAVE, TARGET  :: meteogram_global_data(1:max_dom) !< collected buffers (on IO PE)
  TYPE(t_meteogram_file), SAVE, TARGET  :: meteogram_file_info(1:max_dom)   !< meteogram file handle etc.

  TYPE(t_ncid),       SAVE, TARGET      :: ncid_list(1:max_dom)             !< NetCDF dimension IDs
  TYPE(t_var),        SAVE, TARGET      :: var_list(1:max_dom)              !< internal indices of variables
  
  !! -- data for distributed meteogram sampling (MPI) --
  CHARACTER,          SAVE, ALLOCATABLE :: msg_buffer(:,:)                  !< MPI buffer for station data
  INTEGER,            SAVE              :: max_buf_size                     !< max buffer size for MPI messages
  LOGICAL,            SAVE              :: l_is_sender, l_is_receiver
  INTEGER,            SAVE              :: io_rank                          !< rank of PE which gathers data
  INTEGER,            SAVE              :: owner(MAX_NUM_STATIONS)          !< rank of sender PE for each station

CONTAINS

  !>
  !! Set up list of variables for sampling.
  !!
  !! Note: "extra2d" and "extra3d" variable fields are defined below.
  !!
  !! @par Revision History
  !! Initial implementation  by  F. Prill, DWD (2011-11-09)
  !!
  SUBROUTINE meteogram_setup_variables(ext_data, p_nh_state, prm_diag, p_lnd_state, jg)
    ! atmosphere external data
    TYPE(t_external_data),               INTENT(IN) :: ext_data
    ! nonhydrostatic state
    TYPE(t_nh_state), TARGET,            INTENT(IN)  :: p_nh_state
    ! physical model state and other auxiliary variables
    TYPE(t_nwp_phy_diag),                INTENT(IN)  :: prm_diag
    ! model state for the NWP land physics
    TYPE(t_lnd_state), TARGET,           INTENT(IN)  :: p_lnd_state
    ! patch index
    INTEGER,                             INTENT(IN) :: jg

    TYPE(t_nh_prog)          , POINTER :: prog
    TYPE(t_nh_diag)          , POINTER :: diag
    TYPE(t_nh_metrics)       , POINTER :: metrics
    TYPE(t_lnd_prog)         , POINTER :: p_lnd_prog
    TYPE(t_lnd_diag)         , POINTER :: p_lnd_diag

    diag       => p_nh_state%diag
    prog       => p_nh_state%prog(nnow(jg))
    metrics    => p_nh_state%metrics
    p_lnd_prog => p_lnd_state%prog_lnd(nnow(jg))
    p_lnd_diag => p_lnd_state%diag_lnd

    ! -- atmosphere
    CALL add_atmo_var(VAR_GROUP_ATMO_ML, "P", "Pa", "Pressure", jg, diag%pres(:,:,:))
    CALL add_atmo_var(VAR_GROUP_ATMO_ML, "T", "K", "Temperature", jg, diag%temp(:,:,:))
    CALL add_atmo_var(VAR_GROUP_ATMO_ML, "QV", "kg kg-1", "specific humidity", jg, &
      &               prog%tracer_ptr(iqv)%p_3d(:,:,:))
    CALL add_atmo_var(VAR_GROUP_ATMO_ML, "QC", "kg kg-1", "specific cloud water content", &
      &               jg, prog%tracer_ptr(iqc)%p_3d(:,:,:))
    CALL add_atmo_var(VAR_GROUP_ATMO_ML, "QI", "kg kg-1", "specific cloud ice content", &
      &               jg, prog%tracer_ptr(iqi)%p_3d(:,:,:))
    CALL add_atmo_var(VAR_GROUP_ATMO_ML, "RHO", "kg/m^3", "Density", jg, prog%rho(:,:,:))
    CALL add_atmo_var(VAR_GROUP_ATMO_ML, "PEXNER", "-", "Exner pressure", jg, prog%exner(:,:,:))
    CALL add_atmo_var(VAR_GROUP_ATMO_ML, "THETAV", "K", "virtual potential temperature", &
      &               jg, prog%theta_v(:,:,:))
    CALL add_atmo_var(VAR_GROUP_ATMO_ML, "RHOTHETAV", "K*kg/m^3", "rho*theta_v", jg, &
      &               prog%rhotheta_v(:,:,:))
    CALL add_atmo_var(VAR_GROUP_ATMO_ML, "U", "m/s", "zonal wind", jg, diag%u(:,:,:))
    CALL add_atmo_var(VAR_GROUP_ATMO_ML, "V", "m/s", "meridional wind", jg, diag%v(:,:,:))
    CALL add_atmo_var(VAR_GROUP_ATMO_HL, "W", "m/s", "orthogonal vertical wind", jg, prog%w(:,:,:))

    ! -- soil
    CALL add_atmo_var(VAR_GROUP_SOIL_MLp2, "t_so", "K", "soil temperature", jg, &
      &               p_lnd_prog%t_so(:,:,:,:))
    CALL add_atmo_var(VAR_GROUP_SOIL_ML, "w_so", "m H2O", &
      &               "total water content (ice + liquid water)", &
      &               jg, p_lnd_prog%w_so(:,:,:,:))
    CALL add_atmo_var(VAR_GROUP_SOIL_ML, "w_so_ice", "m H2O", &
      &               "ice content", jg, p_lnd_prog%w_so_ice(:,:,:,:))

    ! -- surface
    CALL add_sfc_var(VAR_GROUP_SURFACE,  "P_SFC", "Pa", "surface pressure", jg, diag%pres_sfc(:,:))
    CALL add_sfc_var(VAR_GROUP_SURFACE,  "PL_Cov", "-", "ground fraction covered by plants", &
      &              jg, ext_data%atm%plcov_mx(:,:))
    CALL add_sfc_var(VAR_GROUP_SURFACE,  "LA_Ind", "-", "leaf area index (vegetation period)", &
      &              jg, ext_data%atm%lai_mx(:,:))
    CALL add_sfc_var(VAR_GROUP_SURFACE,  "RO_Dept", "-", "root depth", jg, &
      &              ext_data%atm%rootdp(:,:))
    CALL add_sfc_var(VAR_GROUP_SURFACE,  "RL*g", "m", "roughness length", jg, prm_diag%gz0(:,:))
    CALL add_sfc_var(VAR_GROUP_SURFACE,  "qv_s", "kg/kg", "specific humidity at the surface", &
      &              jg, p_lnd_diag%qv_s(:,:))
    CALL add_sfc_var(VAR_GROUP_SURFACE,  "w_i", "m H2O", "water content of interception water", &
      &              jg, p_lnd_prog%w_i(:,:,:))
    CALL add_sfc_var(VAR_GROUP_SURFACE,  "w_snow", "m H2O", "water content of snow", &
      &              jg, p_lnd_prog%w_snow(:,:,:))

  END SUBROUTINE meteogram_setup_variables


  !>
  !! Initialize meteogram data buffer, allocating storage.
  !! This is a collective operation.
  !!
  !! @par Revision History
  !! Initial implementation  by  F. Prill, DWD (2011-08-22)
  !!
  SUBROUTINE meteogram_init(meteogram_output_config, ptr_patch, ext_data, p_nh_state, &
    &                       prm_diag, p_lnd_state, iforcing, jg)
    ! station data from namelist
    TYPE(t_meteogram_output_config), TARGET, INTENT(IN) :: meteogram_output_config
    ! data structure containing grid info:
    TYPE(t_patch),                       INTENT(IN) :: ptr_patch
    ! atmosphere external data
    TYPE(t_external_data),               INTENT(IN) :: ext_data
    ! nonhydrostatic state
    TYPE(t_nh_state), TARGET,            INTENT(IN)  :: p_nh_state
    ! physical model state and other auxiliary variables
    TYPE(t_nwp_phy_diag),                INTENT(IN)  :: prm_diag
    ! model state for the NWP land physics
    TYPE(t_lnd_state), TARGET,           INTENT(IN)  :: p_lnd_state
    ! parameterized forcing (right hand side) of dynamics, affects
    ! topography specification, see "mo_extpar_config"
    INTEGER,                             INTENT(IN) :: iforcing
    ! patch index
    INTEGER,                             INTENT(IN) :: jg

    ! local variables:
    CHARACTER(*), PARAMETER :: routine = TRIM("mo_meteogram_output:meteogram_init")

    INTEGER      :: ithis_nlocal_pts, nblks, npromz,  &
      &             nstations, ierrstat,              &
      &             jb, jc, glb_index, i_startidx,    &
      &             i_endidx, jc_station, jb_station, &
      &             istation, my_id, nvars, nlevs,    &
      &             ivar, nsfcvars, iextra,           &
      &             nblks_global, npromz_global, ilev
    REAL(gk)     :: in_points(nproma,meteogram_output_config%nblks,2) !< geographical locations
    REAL(gk)     :: min_dist(nproma,meteogram_output_config%nblks)    !< minimal distance
    ! list of triangles containing lon-lat grid points (first dim: index and block)
    INTEGER      :: tri_idx(2,nproma,meteogram_output_config%nblks)
    INTEGER      :: nlocal_pts(p_n_work)
    LOGICAL      :: l_is_io_pe
    REAL(wp)     :: pi_180
    INTEGER      :: max_var_size, max_sfcvar_size
    REAL(wp), ALLOCATABLE              :: hlevels(:)
    TYPE(t_meteogram_data)   , POINTER :: meteogram_data
    TYPE(t_meteogram_station), POINTER :: p_station
    TYPE(t_nh_diag)          , POINTER :: diag
    TYPE(t_nh_metrics)       , POINTER :: metrics
    TYPE(t_cf_global)        , POINTER :: cf  !< meta info

    pi_180 = ATAN(1._wp)/45._wp

    meteogram_data => meteogram_local_data(jg)

    var_list(jg)%no_atmo_vars = 0
    var_list(jg)%no_sfc_vars  = 0

    diag       => p_nh_state%diag
    metrics    => p_nh_state%metrics

    ! set meta data (appears in NetCDF output file)
    cf => meteogram_file_info(jg)%cf
    cf%title       = 'ICON Meteogram File'
    cf%institution = 'Max Planck Institute for Meteorology/Deutscher Wetterdienst'
    cf%source      = 'icon-dev'
    cf%history     = ''

    ! ------------------------------------------------------------
    ! Distribute stations, determine number of stations located on
    ! this PE:
    ! ------------------------------------------------------------

    ! build an array of geographical coordinates from station list:
    ! in_points(...)
    nstations = meteogram_output_config%nstations
    nblks     = meteogram_output_config%nblks
    npromz    = meteogram_output_config%npromz   

    DO jb=1,nblks
      i_startidx = 1
      i_endidx   = nproma
      IF (jb == nblks) i_endidx = npromz

      DO jc=i_startidx,i_endidx    
        in_points(jc,jb,:) = &
          &  (/ meteogram_output_config%station_list(jc,jb)%location%lon, &
          &     meteogram_output_config%station_list(jc,jb)%location%lat  /) * pi_180
      END DO
    END DO

    ! build GNAT data structure
    CALL gnat_init_grid(ptr_patch)
    ! perform proximity query
    CALL gnat_query_containing_triangles(ptr_patch, gnat_tree, in_points(:,:,:),    &
      &                                  nproma, nblks, npromz,                     &
      &                                  tri_idx(:,:,:), min_dist(:,:))
    CALL gnat_merge_distributed_queries(ptr_patch, nstations, nproma, nblks, min_dist,  &
      &                                 tri_idx(:,:,:), in_points(:,:,:),               &
      &                                 nlocal_pts(:), owner(:), ithis_nlocal_pts)
    nblks    = ithis_nlocal_pts/nproma + 1
    npromz   = ithis_nlocal_pts - (nblks-1)*nproma
    meteogram_data%nstations = ithis_nlocal_pts
    meteogram_data%nblks     = nblks
    meteogram_data%npromz    = npromz
    ! clean up
    CALL gnat_destroy()

    ! ------------------------------------------------------------
    ! Initialize local data structure, fill header
    ! ------------------------------------------------------------

    meteogram_data%icurrent  = 0 ! reset current sample index
    meteogram_data%max_nlevs = 1
    ALLOCATE(meteogram_data%time_stamp(MAX_TIME_STAMPS), &
      &      meteogram_data%sfc_var_info(MAX_NSFCVARS),  &
      &      meteogram_data%var_info(MAX_NVARS),         &
      &      stat=ierrstat)
    IF (ierrstat /= SUCCESS) THEN
      CALL finish (routine, 'ALLOCATE of meteogram data structures failed')
    ENDIF

    ! set up list of variables:
    CALL meteogram_setup_variables(ext_data, p_nh_state, prm_diag, p_lnd_state, jg)
    
    IF (lwrite_extra) THEN
      ! Variable: Extra 2D
      DO iextra=1,inextra_2d
        CALL add_sfc_var(VAR_GROUP_SURFACE, "EXTRA2D"//int2string(iextra), "", &
          &              "-", jg, diag%extra_2d(:,:, iextra))
      END DO
      ! Variable: Extra 3D
      DO iextra=1,inextra_3d
        CALL add_atmo_var(VAR_GROUP_ATMO_ML, "EXTRA3D"//int2string(iextra), "", &
          &               "-", jg, diag%extra_3d(:,:,:,iextra))
      END DO
    END IF

    meteogram_data%nsfcvars  = var_list(jg)%no_sfc_vars
    meteogram_data%nvars     = var_list(jg)%no_atmo_vars
    nsfcvars = meteogram_data%nsfcvars
    nvars    = meteogram_data%nvars

    ! set up list of stations:
    ALLOCATE(meteogram_data%station(nproma, nblks), stat=ierrstat)
    IF (ierrstat /= SUCCESS) THEN
      CALL finish (routine, 'ALLOCATE of meteogram data structures failed')
    ENDIF

    my_id = get_my_mpi_all_id()
    istation = 0
    DO jb=1,nblks
      i_startidx = 1
      i_endidx   = nproma
      IF (jb == nblks) i_endidx = npromz

      DO jc=i_startidx,i_endidx
        ! find corresponding entry in station list:
        DO
          istation = istation + 1
          IF (owner(istation) == my_id) EXIT
        END DO
        jb_station = istation/nproma + 1
        jc_station = istation - (jb_station-1)*nproma
        meteogram_data%station(jc,jb)%station_idx = (/ jc_station, jb_station /)
        ! set owner ID:
        meteogram_data%station(jc,jb)%owner = my_id
        ! set local triangle index, block:
        meteogram_data%station(jc,jb)%tri_idx_local(1:2) = tri_idx(1:2,jc,jb)
        ! translate local index to global index:
        glb_index = ptr_patch%cells%glb_index(idx_1d(tri_idx(1,jc,jb), &
          &                                          tri_idx(2,jc,jb)))
        meteogram_data%station(jc,jb)%tri_idx(1:2) =  &
          &  (/ idx_no(glb_index), blk_no(glb_index) /)
        ! set Coriolis parameter for station
        meteogram_data%station(jc,jb)%fc           =  &
          &  ptr_patch%cells%f_c(tri_idx(1,jc,jb), tri_idx(2,jc,jb))
        ! set station information on height, soil type etc.:
        SELECT CASE ( iforcing )
        CASE ( inwp ) ! NWP physics
          meteogram_data%station(jc,jb)%hsurf    =  &
            &  ext_data%atm%topography_c(tri_idx(1,jc,jb), tri_idx(2,jc,jb))
          meteogram_data%station(jc,jb)%frland   =  &
            &  ext_data%atm%fr_land(tri_idx(1,jc,jb), tri_idx(2,jc,jb))
          meteogram_data%station(jc,jb)%soiltype =  &
            &  ext_data%atm%soiltyp(tri_idx(1,jc,jb), tri_idx(2,jc,jb))
        CASE DEFAULT
          meteogram_data%station(jc,jb)%hsurf    =  0._wp
          meteogram_data%station(jc,jb)%frland   =  0._wp
          meteogram_data%station(jc,jb)%soiltype =  0
        END SELECT
        ! initialize value buffer and set level heights:
        ALLOCATE(meteogram_data%station(jc,jb)%var(nvars),    &
          &      hlevels(meteogram_local_data(jg)%max_nlevs), &
          &      stat=ierrstat)
        IF (ierrstat /= SUCCESS) THEN
          CALL finish (routine, 'ALLOCATE of meteogram data structures failed')
        ENDIF
        DO ivar=1,nvars
          nlevs = meteogram_data%var_info(ivar)%nlevs
          ALLOCATE(meteogram_data%station(jc,jb)%var(ivar)%values(nlevs, MAX_TIME_STAMPS), &
            &      meteogram_data%station(jc,jb)%var(ivar)%heights(nlevs),                 &
            &      stat=ierrstat)
          IF (ierrstat /= SUCCESS) THEN
            CALL finish (routine, 'ALLOCATE of meteogram data structures failed')
          ENDIF
          ! initialize level heights:
          SELECT CASE(meteogram_data%var_info(ivar)%igroup_id)
          CASE(VAR_GROUP_ATMO_ML)
            ! model level heights
            hlevels(1:nlevs) = metrics%z_mc(tri_idx(1,jc,jb), 1:nlevs, tri_idx(2,jc,jb))
          CASE(VAR_GROUP_ATMO_HL)
            ! half level heights
            hlevels(1:nlevs) = metrics%z_ifc(tri_idx(1,jc,jb), 1:nlevs, tri_idx(2,jc,jb))
          CASE(VAR_GROUP_SOIL_ML)
            ! soil half level heights
            hlevels(1:nlevs) = zml_soil(1:nlevs)
          CASE(VAR_GROUP_SOIL_MLp2)
            ! soil half level heights PLUS surface level
            hlevels(1)       = 0._wp
            hlevels(2:nlevs) = zml_soil(1:(nlevs-1))
          CASE DEFAULT
            CALL finish (routine, 'Invalid group ID.')
          END SELECT
          DO ilev = 1,nlevs
            meteogram_data%station(jc,jb)%var(ivar)%heights(ilev) =    &
              &  hlevels(meteogram_data%var_info(ivar)%levels(ilev))
          END DO
        END DO
        DEALLOCATE(hlevels, stat=ierrstat)
        IF (ierrstat /= SUCCESS) THEN
          CALL finish (routine, 'DEALLOCATE of helper variable failed.')
        ENDIF
        ! initialize value buffer for surface variables:
        ALLOCATE(meteogram_data%station(jc,jb)%sfc_var(nsfcvars), &
          &      stat=ierrstat)
        IF (ierrstat /= SUCCESS) THEN
          CALL finish (routine, 'ALLOCATE of meteogram data structures failed')
        ENDIF
        DO ivar=1,nsfcvars
          ALLOCATE(meteogram_data%station(jc,jb)%sfc_var(ivar)%values(MAX_TIME_STAMPS), &
            &      stat=ierrstat)
          IF (ierrstat /= SUCCESS) THEN
            CALL finish (routine, 'ALLOCATE of meteogram data structures failed')
          ENDIF
        END DO
      END DO
    END DO

    ! ------------------------------------------------------------
    ! If this is the IO PE: initialize global data structure
    ! ------------------------------------------------------------

    IO_PE : IF (my_process_is_stdio() .AND. .NOT. meteogram_output_config%ldistributed) THEN

      meteogram_global_data(jg)%nstations =  meteogram_output_config%nstations
      meteogram_global_data(jg)%nblks     =  meteogram_output_config%nblks    
      meteogram_global_data(jg)%npromz    =  meteogram_output_config%npromz   

      ! Note: variable info is not duplicated
      meteogram_global_data(jg)%nvars     =  meteogram_local_data(jg)%nvars
      meteogram_global_data(jg)%nsfcvars  =  meteogram_local_data(jg)%nsfcvars
      meteogram_global_data(jg)%max_nlevs =  meteogram_local_data(jg)%max_nlevs
      meteogram_global_data(jg)%var_info      =>  meteogram_local_data(jg)%var_info
      meteogram_global_data(jg)%sfc_var_info  =>  meteogram_local_data(jg)%sfc_var_info
      meteogram_global_data(jg)%time_stamp    =>  meteogram_local_data(jg)%time_stamp

      nblks_global  = meteogram_global_data(jg)%nblks
      npromz_global = meteogram_global_data(jg)%npromz
      ALLOCATE(meteogram_global_data(jg)%station(nproma, nblks_global), stat=ierrstat)
      IF (ierrstat /= SUCCESS) THEN
        CALL finish (routine, 'ALLOCATE of meteogram data structures failed')
      ENDIF

      DO jb=1,nblks_global
        i_startidx = 1
        i_endidx   = nproma
        IF (jb == nblks_global) i_endidx = npromz_global
        
        DO jc=i_startidx,i_endidx
          p_station => meteogram_global_data(jg)%station(jc,jb)

          ALLOCATE(p_station%var(nvars), stat=ierrstat)
          IF (ierrstat /= SUCCESS) &
            CALL finish (routine, 'ALLOCATE of meteogram data structures failed')
          
          DO ivar=1,nvars
            nlevs = meteogram_data%var_info(ivar)%nlevs
            ALLOCATE(p_station%var(ivar)%values(nlevs, MAX_TIME_STAMPS), &
              &      p_station%var(ivar)%heights(nlevs),                 &
              &      stat=ierrstat)
            IF (ierrstat /= SUCCESS) &
              CALL finish (routine, 'ALLOCATE of meteogram data structures failed')
          END DO
          ALLOCATE(p_station%sfc_var(nsfcvars), stat=ierrstat)
          IF (ierrstat /= SUCCESS) THEN
            CALL finish (routine, 'ALLOCATE of meteogram data structures failed')
          ENDIF
          DO ivar=1,nsfcvars
            ALLOCATE(p_station%sfc_var(ivar)%values(MAX_TIME_STAMPS), stat=ierrstat)
            IF (ierrstat /= SUCCESS) THEN
              CALL finish (routine, 'ALLOCATE of meteogram data structures failed')
            ENDIF
          END DO

        END DO
      END DO

    END IF IO_PE

    ! ------------------------------------------------------------
    ! initialize MPI buffer
    ! ------------------------------------------------------------
    
    IF (.NOT. meteogram_output_config%ldistributed) THEN
      io_rank       = get_mpi_all_workroot_id()
      l_is_io_pe    = (get_my_mpi_all_id() == io_rank)
      l_is_sender   = .NOT. l_is_io_pe
      l_is_receiver = l_is_io_pe

      ! compute maximum buffer size for MPI messages:
      ! (max_var_size: contains also height levels)
      max_var_size    = (MAX_TIME_STAMPS+1)*p_real_dp_byte*meteogram_data%max_nlevs
      max_sfcvar_size = MAX_TIME_STAMPS*p_real_dp_byte
      max_buf_size    = MAX_HEADER_SIZE*p_real_dp_byte          &
        &               + meteogram_data%nvars*max_var_size         &
        &               + meteogram_data%nsfcvars*max_sfcvar_size 

      ! allocate buffer:
      IF (l_is_receiver) THEN
        ALLOCATE(msg_buffer(max_buf_size, MAX_NUM_STATIONS), stat=ierrstat)  
        IF (ierrstat /= SUCCESS) THEN
          CALL finish (routine, 'ALLOCATE of meteogram message buffer failed')
        ENDIF
      ELSE
        ! allocate buffer:
        ALLOCATE(msg_buffer(max_buf_size, 1), stat=ierrstat)  
        IF (ierrstat /= SUCCESS) THEN
          CALL finish (routine, 'ALLOCATE of meteogram message buffer failed')
        ENDIF
      END IF

    END IF

    ! ------------------------------------------------------------
    ! If this is the IO PE: open NetCDF file
    ! ------------------------------------------------------------

    CALL meteogram_open_file(meteogram_output_config, jg)

  END SUBROUTINE meteogram_init


  !>
  !! @return .TRUE. if meteogram data will be recorded for this step.
  !!
  !! @par Revision History
  !! Initial implementation  by  F. Prill, DWD (2011-08-22)
  !!
  FUNCTION meteogram_is_sample_step(meteogram_output_config, cur_step)
    LOGICAL :: meteogram_is_sample_step
    ! station data from namelist
    TYPE(t_meteogram_output_config), TARGET, INTENT(IN) :: meteogram_output_config
    INTEGER,          INTENT(IN)  :: cur_step     !< current model iteration step

    meteogram_is_sample_step = &
      &  meteogram_output_config%lenabled               .AND. &
      &  (cur_step >= meteogram_output_config%n0_mtgrm) .AND. &
      &  (MOD((cur_step - meteogram_output_config%n0_mtgrm),  &
      &       meteogram_output_config%ninc_mtgrm) == 0)
    
  END FUNCTION meteogram_is_sample_step

  
  !>
  !! Adds values for current model time to buffer.
  !! For gathered NetCDF output, this is a collective operation,
  !! otherwise this is a local operation.
  !!
  !! @par Revision History
  !! Initial implementation  by  F. Prill, DWD (2011-08-22)
  !!
  SUBROUTINE meteogram_sample_vars(jg, cur_step, cur_datetime, ierr)
    INTEGER,          INTENT(IN)  :: jg           !< patch index
    INTEGER,          INTENT(IN)  :: cur_step     !< current model iteration step
    TYPE(t_datetime), INTENT(IN)  :: cur_datetime !< date and time of point sample
    INTEGER,          INTENT(OUT) :: ierr         !< error code (e.g. buffer overflow)
    ! local variables
    CHARACTER(*), PARAMETER :: routine = TRIM("mo_meteogram_output:meteogram_sample_vars")
    INTEGER :: jb, jc, i_startidx, i_endidx, ilev, ivar,  &
      &        i_tstep, iidx, iblk
    TYPE(t_meteogram_data), POINTER :: meteogram_data

    meteogram_data => meteogram_local_data(jg)
    ierr = 0

    IF (dbg_level > 0) THEN
      WRITE(message_text,*) "Sampling at step=", cur_step
      CALL message(routine, TRIM(message_text))
    END IF

    ! increase time step counter
    meteogram_data%icurrent = meteogram_data%icurrent + 1
    i_tstep = meteogram_data%icurrent
    IF (i_tstep > MAX_TIME_STAMPS) THEN
      ! buffer full
      ierr = -1
      RETURN
    END IF

    meteogram_data%time_stamp(i_tstep)%istep = cur_step
    meteogram_data%time_stamp(i_tstep)%zdate = iso8601(cur_datetime)

    ! fill time step with values
    DO jb=1,meteogram_data%nblks
      i_startidx = 1
      i_endidx   = nproma
      IF (jb == meteogram_data%nblks) i_endidx = meteogram_data%npromz

      DO jc=i_startidx,i_endidx
        iidx  = meteogram_data%station(jc,jb)%tri_idx_local(1)
        iblk  = meteogram_data%station(jc,jb)%tri_idx_local(2)

        DO ivar=1,var_list(jg)%no_atmo_vars
          IF (ASSOCIATED(meteogram_data%var_info(ivar)%p_source_3d)) THEN
            DO ilev=1,meteogram_data%var_info(ivar)%nlevs
              meteogram_data%station(jc,jb)%var(ivar)%values(ilev, i_tstep) =    &
                &  meteogram_data%var_info(ivar)%p_source_3d(                    &
                &       iidx, meteogram_data%var_info(ivar)%levels(ilev), iblk )
            END DO
          ELSE 
            CALL finish (routine, 'Source array not associated!')
          END IF
        END DO

        DO ivar=1,var_list(jg)%no_sfc_vars
          meteogram_data%station(jc,jb)%sfc_var(ivar)%values(i_tstep) =  &
            &  meteogram_data%sfc_var_info(ivar)%p_source( iidx, iblk )
        END DO

      END DO
    END DO

  END SUBROUTINE meteogram_sample_vars


  !>
  !! Destroy meteogram data structure.
  !! For gathered NetCDF output, this is a collective operation,
  !! otherwise this is a local operation.
  !!
  !! @par Revision History
  !! Initial implementation  by  F. Prill, DWD (2011-08-22)
  !!
  SUBROUTINE meteogram_finalize(meteogram_output_config, jg)
    ! station data from namelist
    TYPE(t_meteogram_output_config), TARGET, INTENT(IN) :: meteogram_output_config
    INTEGER, INTENT(IN)  :: jg    !< patch index
    ! local variables:
    CHARACTER(*), PARAMETER     :: routine = TRIM("mo_meteogram_output:meteogram_finalize")
    INTEGER                     :: ierrstat, jb, jc, i_startidx, i_endidx, &
      &                            nvars, nsfcvars, ivar
    TYPE(t_meteogram_data), POINTER :: meteogram_data

    meteogram_data => meteogram_local_data(jg)

    ! ------------------------------------------------------------
    ! If this is the IO PE: close NetCDF file
    ! ------------------------------------------------------------

    CALL meteogram_close_file(meteogram_output_config, jg)

    DEALLOCATE(meteogram_data%time_stamp, stat=ierrstat)
    IF (ierrstat /= SUCCESS) THEN
      CALL finish (routine, 'DEALLOCATE of meteogram data structures failed')
    ENDIF

    nvars    = meteogram_data%nvars
    nsfcvars = meteogram_data%nsfcvars
    DO ivar=1,nvars
      IF (ASSOCIATED(meteogram_data%var_info(ivar)%levels)) THEN
        DEALLOCATE(meteogram_data%var_info(ivar)%levels, stat=ierrstat)
        IF (ierrstat /= SUCCESS) THEN
          CALL finish (routine, 'DEALLOCATE of meteogram data structures failed')
        ENDIF
      END IF
    END DO
    DEALLOCATE(meteogram_data%var_info, stat=ierrstat)
    IF (ierrstat /= SUCCESS) THEN
      CALL finish (routine, 'DEALLOCATE of meteogram data structures failed')
    ENDIF
    DEALLOCATE(meteogram_data%sfc_var_info, stat=ierrstat)
    IF (ierrstat /= SUCCESS) THEN
      CALL finish (routine, 'DEALLOCATE of meteogram data structures failed')
    ENDIF
    
    DO jb=1,meteogram_data%nblks
      i_startidx = 1
      i_endidx   = nproma
      IF (jb == meteogram_data%nblks) &
        &  i_endidx = meteogram_data%npromz

      DO jc=i_startidx,i_endidx
        DO ivar=1,nvars
          DEALLOCATE(meteogram_data%station(jc,jb)%var(ivar)%values,  &
            &        meteogram_data%station(jc,jb)%var(ivar)%heights, &
            &        stat=ierrstat)
          IF (ierrstat /= SUCCESS) THEN
            CALL finish (routine, 'DEALLOCATE of meteogram data structures failed')
          ENDIF
        END DO
        DEALLOCATE(meteogram_data%station(jc,jb)%var, stat=ierrstat)
        IF (ierrstat /= SUCCESS) THEN
          CALL finish (routine, 'DEALLOCATE of meteogram data structures failed')
        ENDIF
        DO ivar=1,nsfcvars
          DEALLOCATE(meteogram_data%station(jc,jb)%sfc_var(ivar)%values, &
            &        stat=ierrstat)
          IF (ierrstat /= SUCCESS) THEN
            CALL finish (routine, 'DEALLOCATE of meteogram data structures failed')
          ENDIF
        END DO
        DEALLOCATE(meteogram_data%station(jc,jb)%sfc_var, stat=ierrstat)
        IF (ierrstat /= SUCCESS) THEN
          CALL finish (routine, 'DEALLOCATE of meteogram data structures failed')
        ENDIF

      END DO
    END DO
    DEALLOCATE(meteogram_data%station, stat=ierrstat)
    IF (ierrstat /= SUCCESS) THEN
      CALL finish (routine, 'DEALLOCATE of meteogram data structures failed')
    ENDIF
    meteogram_data%nvars     = 0
    meteogram_data%nsfcvars  = 0
    meteogram_data%nstations = 0

    ! deallocate MPI buffer
    IF (.NOT. meteogram_output_config%ldistributed) THEN
      DEALLOCATE(msg_buffer, stat=ierrstat)
      IF (ierrstat /= SUCCESS) THEN
        CALL finish (routine, 'DEALLOCATE of MPI message buffer failed')
      ENDIF
    END IF

    ! deallocate global meteogram data

    IO_PE : IF (my_process_is_stdio() .AND. .NOT. meteogram_output_config%ldistributed) THEN
    
      DO jb=1,meteogram_global_data(jg)%nblks
        i_startidx = 1
        i_endidx   = nproma
        IF (jb == meteogram_global_data(jg)%nblks) i_endidx = meteogram_global_data(jg)%npromz
        
        DO jc=i_startidx,i_endidx
          DO ivar=1,nvars
            DEALLOCATE(meteogram_global_data(jg)%station(jc,jb)%var(ivar)%values,  &
              &        meteogram_global_data(jg)%station(jc,jb)%var(ivar)%heights, &
              &        stat=ierrstat)
            IF (ierrstat /= SUCCESS) &
              CALL finish (routine, 'ALLOCATE of meteogram data structures failed')
          END DO
          DEALLOCATE(meteogram_global_data(jg)%station(jc,jb)%var, stat=ierrstat)
          IF (ierrstat /= SUCCESS) &
            CALL finish (routine, 'ALLOCATE of meteogram data structures failed')

          DO ivar=1,nsfcvars
            DEALLOCATE(meteogram_global_data(jg)%station(jc,jb)%sfc_var(ivar)%values, &
              &        stat=ierrstat)
            IF (ierrstat /= SUCCESS) &
              CALL finish (routine, 'ALLOCATE of meteogram data structures failed')
          END DO
          DEALLOCATE(meteogram_global_data(jg)%station(jc,jb)%sfc_var, stat=ierrstat)
          IF (ierrstat /= SUCCESS) &
            CALL finish (routine, 'ALLOCATE of meteogram data structures failed')

        END DO
      END DO
      DEALLOCATE(meteogram_global_data(jg)%station, stat=ierrstat)
      IF (ierrstat /= SUCCESS) &
        CALL finish (routine, 'ALLOCATE of meteogram data structures failed')
      
    END IF IO_PE

  END SUBROUTINE meteogram_finalize


  !>
  !! IO PE gathers all buffer information from working PEs and copies
  !! the contents to the global meteogram buffer.
  !! Afterwards, all working PEs flush their local buffers.
  !!
  !! For gathered NetCDF output, this is a collective operation.
  !!
  !! @par Revision History
  !! Initial implementation  by  F. Prill, DWD (2011-09-10)
  !!
  SUBROUTINE meteogram_collect_buffers(meteogram_output_config, jg, io_rank)
    ! station data from namelist
    TYPE(t_meteogram_output_config), TARGET, INTENT(IN) :: meteogram_output_config
    INTEGER, INTENT(IN)  :: jg       !< patch index
    INTEGER, INTENT(IN)  :: io_rank  !< MPI rank of root process

#ifndef NOMPI
    ! local variables
    INTEGER     :: station_idx(2), position, icurrent,   &
      &            jb, jc, i_startidx, i_endidx,         &
      &            istation, ivar, nlevs
    INTEGER     :: myrank
    TYPE(t_meteogram_data),    POINTER :: meteogram_data
    TYPE(t_meteogram_station), POINTER :: p_station

    meteogram_data => meteogram_local_data(jg)

    ! get global rank and MPI communicator
    myrank = get_my_mpi_all_id()

    ! global time stamp index
    ! Note: We assume that this value is identical for all PEs
    icurrent = meteogram_data%icurrent

    ! -- RECEIVER CODE --
    RECEIVER : IF (l_is_receiver) THEN
      ! launch MPI message requests for station data on foreign PEs
      DO istation=1,meteogram_output_config%nstations
        IF (owner(istation) /= myrank) THEN
          CALL p_irecv_packed(msg_buffer(:,istation), MPI_ANY_SOURCE, istation, max_buf_size)
        END IF
      END DO

      ! wait for messages to arrive:
      CALL p_wait()

      ! unpack received messages:
      jc = 0
      jb = 1
      DO istation=1,meteogram_output_config%nstations
        IF (owner(istation) /= myrank) THEN
          position = 0
          CALL p_unpack_int_1d(msg_buffer(:,istation),max_buf_size, position, station_idx(:),2)
          p_station => meteogram_global_data(jg)%station(station_idx(1), station_idx(2))
          p_station%station_idx(1:2) = station_idx(1:2)
          ! unpack header information
          CALL p_unpack_int_1d(msg_buffer(:,istation),max_buf_size, position, &
            &                  p_station%tri_idx(:),2)
          CALL p_unpack_int_1d(msg_buffer(:,istation),max_buf_size, position, &
            &                  p_station%tri_idx_local(:),2)
          CALL p_unpack_int(msg_buffer(:,istation),max_buf_size, position, p_station%owner)
          CALL p_unpack_real(msg_buffer(:,istation),max_buf_size, position, p_station%hsurf)
          CALL p_unpack_real(msg_buffer(:,istation),max_buf_size, position, p_station%frland)
          CALL p_unpack_real(msg_buffer(:,istation),max_buf_size, position, p_station%fc)
          CALL p_unpack_int(msg_buffer(:,istation),max_buf_size, position, p_station%soiltype)

          ! unpack heights and meteogram data:
          DO ivar=1,meteogram_data%nvars
            nlevs = meteogram_data%var_info(ivar)%nlevs
            CALL p_unpack_real_1d(msg_buffer(:,istation),max_buf_size, position, &
              &                   p_station%var(ivar)%heights(:), nlevs)
            CALL p_unpack_real_2d(msg_buffer(:,istation),max_buf_size, position, &
              &                   p_station%var(ivar)%values(:,:), nlevs*icurrent)
          END DO

          DO ivar=1,meteogram_data%nsfcvars
            CALL p_unpack_real_1d(msg_buffer(:,istation),max_buf_size, position, &
              &                   p_station%sfc_var(ivar)%values(:), icurrent)
          END DO
        ELSE
          ! this PE is both sender and receiver - direct copy:
          jc = jc + 1
          IF (jc > nproma) THEN
            jc = 1
            jb = jb + 1
          END IF
          station_idx(1:2) = meteogram_data%station(jc,jb)%station_idx(1:2)

          p_station => meteogram_global_data(jg)%station(station_idx(1),station_idx(2))
          p_station%station_idx(1:2)   = meteogram_data%station(jc,jb)%station_idx(1:2)
          p_station%tri_idx(1:2)       = meteogram_data%station(jc,jb)%tri_idx(1:2)
          p_station%tri_idx_local(1:2) = meteogram_data%station(jc,jb)%tri_idx_local(1:2)
          p_station%owner              = meteogram_data%station(jc,jb)%owner
          p_station%hsurf              = meteogram_data%station(jc,jb)%hsurf
          p_station%frland             = meteogram_data%station(jc,jb)%frland
          p_station%fc                 = meteogram_data%station(jc,jb)%fc
          p_station%soiltype           = meteogram_data%station(jc,jb)%soiltype
          ! copy heights and meteogram data
          DO ivar=1,meteogram_data%nvars
            nlevs = meteogram_data%var_info(ivar)%nlevs
            p_station%var(ivar)%heights(1:nlevs) =  &
              &  meteogram_data%station(jc,jb)%var(ivar)%heights(1:nlevs)
            p_station%var(ivar)%values(1:nlevs, 1:icurrent) =  &
              &  meteogram_data%station(jc,jb)%var(ivar)%values(1:nlevs, 1:icurrent)
          END DO
          DO ivar=1,meteogram_data%nsfcvars
            p_station%sfc_var(ivar)%values(1:icurrent) =  &
              &  meteogram_data%station(jc,jb)%sfc_var(ivar)%values(1:icurrent)
          END DO
        END IF
      END DO
    END IF RECEIVER

    ! -- SENDER CODE --
    SENDER : IF (l_is_sender) THEN
      ! pack station into buffer; send it
      DO jb=1,meteogram_data%nblks
        i_startidx = 1
        i_endidx   = nproma
        IF (jb == meteogram_data%nblks) i_endidx = meteogram_data%npromz
        DO jc=i_startidx,i_endidx

          ! Meteogram header (information on location, ...)
          p_station => meteogram_data%station(jc,jb)
          position = 0
          CALL p_pack_int_1d(p_station%station_idx(:), 2, msg_buffer(:,1), max_buf_size, position)
          CALL p_pack_int_1d(p_station%tri_idx(:), 2, msg_buffer(:,1), max_buf_size, position)
          CALL p_pack_int_1d(p_station%tri_idx_local(:), 2, msg_buffer(:,1), &
            &                max_buf_size, position)
          CALL p_pack_int (p_station%owner, msg_buffer(:,1), max_buf_size, position)
          CALL p_pack_real(p_station%hsurf, msg_buffer(:,1), max_buf_size, position)
          CALL p_pack_real(p_station%frland, msg_buffer(:,1), max_buf_size, position)
          CALL p_pack_real(p_station%fc, msg_buffer(:,1), max_buf_size, position)
          CALL p_pack_int (p_station%soiltype, msg_buffer(:,1), max_buf_size, position)
          ! pack heights and meteogram data:
          DO ivar=1,meteogram_data%nvars
            nlevs = meteogram_data%var_info(ivar)%nlevs
            CALL p_pack_real_1d(p_station%var(ivar)%heights(:), nlevs,     &
              &                 msg_buffer(:,1), max_buf_size, position)
            CALL p_pack_real_2d(p_station%var(ivar)%values(:,:), nlevs*icurrent, &
              &                 msg_buffer(:,1), max_buf_size, position)
          END DO
          DO ivar=1,meteogram_data%nsfcvars
            CALL p_pack_real_1d(p_station%sfc_var(ivar)%values(:), icurrent,     &
              &                 msg_buffer(:,1), max_buf_size, position)
          END DO

          ! (blocking) send of packed station data to IO PE:
          istation = nproma*(p_station%station_idx(2) - 1) + p_station%station_idx(1)
          CALL p_send_packed(msg_buffer(:,1), io_rank, istation, position)
        END DO
      END DO

      ! reset buffer on sender side
      IF (.NOT. l_is_receiver) &
        meteogram_data%icurrent = 0

    END IF SENDER

#endif
  END SUBROUTINE meteogram_collect_buffers


  !>
  !! The IO PE creates and opens a disk file for output.
  !! For gathered NetCDF output, this is a collective operation,
  !! otherwise this is a local operation for the IO PE.
  !!
  !! @par Revision History
  !! Initial implementation  by  F. Prill, DWD (2011-08-22)
  !!
  SUBROUTINE meteogram_open_file(meteogram_output_config, jg)
    ! station data from namelist
    TYPE(t_meteogram_output_config), TARGET, INTENT(IN) :: meteogram_output_config
    ! patch index
    INTEGER,                             INTENT(IN) :: jg
    ! local variables:
    CHARACTER(*), PARAMETER :: routine = TRIM("mo_meteogram_output:meteogram_open_file")
    INTEGER                     :: jb, jc, i_startidx, i_endidx, old_mode, ncfile, &
      &                            istation, ivar, nvars, nsfcvars, nlevs
    TYPE(t_ncid),       POINTER :: ncid
    INTEGER                     :: station_name_dims(2), var_name_dims(2), &
      &                            var_level_dims(2), time_string_dims(2), &
      &                            var_dims(4),  sfcvar_dims(3),           &
      &                            height_level_dims(3),                   &
      &                            istart2(2), icount2(2)
    TYPE(t_station_list)  , POINTER :: this_station
    TYPE(t_meteogram_data), POINTER :: meteogram_data
    TYPE(t_cf_global)     , POINTER :: cf  !< meta info

    IF (meteogram_output_config%ftype /= FTYPE_NETCDF) THEN
      CALL finish(routine, "Output format not yet implemented.")
    END IF

    ! In "non-distributed" mode, station data is gathered by PE #0
    ! which writes a single file.
    ! Note that info on variables is not copied to the global data set
    ! (we use the local meteogram_data there).

    IF (.NOT. meteogram_output_config%ldistributed) THEN
      CALL meteogram_collect_buffers(meteogram_output_config, jg, io_rank)
      meteogram_data => meteogram_global_data(jg)
    ELSE
      meteogram_data => meteogram_local_data(jg)
      ! skip routine, if this PE has nothing to do...
      IF (meteogram_local_data(jg)%nstations == 0) RETURN
    END IF

    ! skip routine, if this PE has nothing to do...
    IF  (.NOT. meteogram_output_config%ldistributed .AND.  &
      & (.NOT. l_is_receiver)) RETURN

    ncid => ncid_list(jg)
    nvars    = meteogram_data%nvars
    nsfcvars = meteogram_data%nsfcvars

    ! create a file name for this PE:
    CALL meteogram_create_filename(meteogram_output_config, jg)

    ! create NetCDF file:
    CALL nf(nf_set_default_format(nf_format_64bit, old_mode))
    CALL nf(nf_create(TRIM(meteogram_file_info(jg)%zname), nf_clobber, &
      &               meteogram_file_info(jg)%file_id))
    ncfile = meteogram_file_info(jg)%file_id
    CALL nf(nf_set_fill(ncfile, nf_nofill, old_mode))

    cf => meteogram_file_info(jg)%cf
    CALL nf(nf_put_att_text(ncfile, NF_GLOBAL, 'title',       &
      &                     LEN_TRIM(cf%title),       TRIM(cf%title)))
    CALL nf(nf_put_att_text(ncfile, NF_GLOBAL, 'history',     &
      &                     LEN_TRIM(cf%history),     TRIM(cf%history)))
    CALL nf(nf_put_att_text(ncfile, NF_GLOBAL, 'institution', &
      &                     LEN_TRIM(cf%institution), TRIM(cf%institution)))
    CALL nf(nf_put_att_text(ncfile, NF_GLOBAL, 'source',      &
      &                     LEN_TRIM(cf%source),      TRIM(cf%source)))
    CALL nf(nf_put_att_text(ncfile, NF_GLOBAL, 'comment',     &
      &                     LEN_TRIM(cf%comment),     TRIM(cf%comment)))
    CALL nf(nf_put_att_text(ncfile, NF_GLOBAL, 'references',  &
      &                     LEN_TRIM(cf%references),  TRIM(cf%references)))

    ! for the definition of a character-string variable define
    ! character-position dimension for strings
    CALL nf(nf_def_dim(ncfile, "stringlen", MAX_DESCR_LENGTH, ncid%charid))
    ! station header:
    CALL nf(nf_def_dim(ncfile, 'nstations',  meteogram_data%nstations, ncid%nstations))
    ! write variables:
    CALL nf(nf_def_dim(ncfile, 'nvars',      meteogram_data%nvars,     ncid%nvars))
    IF (meteogram_data%nsfcvars > 0) THEN
      CALL nf(nf_def_dim(ncfile, 'nsfcvars', meteogram_data%nsfcvars,  ncid%nsfcvars))
    END IF
    CALL nf(nf_def_dim(ncfile, 'max_nlevs',  meteogram_data%max_nlevs, ncid%max_nlevs))
    ! create time dimension:
    CALL nf(nf_def_dim(ncfile, 'time', NF_UNLIMITED, ncid%timeid))
    
    ! create station variables:
    station_name_dims = (/ ncid%charid, ncid%nstations /)
    CALL nf(nf_def_var(ncfile, "station_name", NF_CHAR, 2, station_name_dims(:), &
      &                ncid%station_name))
    CALL nf_add_descr("Station name (character string)", ncfile, ncid%station_name)
    CALL nf(nf_def_var(ncfile, "station_lon", NF_DOUBLE, 1, ncid%nstations, &
      &                ncid%station_lon))
    CALL nf_add_descr("Longitude of meteogram station", ncfile, ncid%station_lon)
    CALL nf(nf_def_var(ncfile, "station_lat", NF_DOUBLE, 1, ncid%nstations, &
      &                ncid%station_lat))
    CALL nf_add_descr("Latitude of meteogram station", ncfile, ncid%station_lat)
    CALL nf(nf_def_var(ncfile, "station_idx", NF_INT, 1, ncid%nstations, &
      &                ncid%station_idx))
    CALL nf_add_descr("Global triangle adjacent to meteogram station (index)", &
      &               ncfile, ncid%station_idx)
    CALL nf(nf_def_var(ncfile, "station_blk", NF_INT, 1, ncid%nstations, &
      &                ncid%station_blk))
    CALL nf_add_descr("Global triangle adjacent to meteogram station (block)", &
      &               ncfile, ncid%station_blk)
    CALL nf(nf_def_var(ncfile, "station_hsurf", NF_DOUBLE, 1, ncid%nstations, &
      &                ncid%station_hsurf))
    CALL nf_add_descr("Meteogram station surface height", ncfile, ncid%station_hsurf)
    CALL nf(nf_def_var(ncfile, "station_frland", NF_DOUBLE, 1, ncid%nstations, &
      &                ncid%station_frland))
    CALL nf_add_descr("Meteogram station land fraction", ncfile, ncid%station_frland)
    CALL nf(nf_def_var(ncfile, "station_fc", NF_DOUBLE, 1, ncid%nstations, &
      &                ncid%station_fc))
    CALL nf_add_descr("Meteogram station Coriolis parameter", ncfile, ncid%station_fc)
    CALL nf(nf_def_var(ncfile, "station_soiltype", NF_INT, 1, ncid%nstations, &
      &                ncid%station_soiltype))
    CALL nf_add_descr("Meteogram station soil type", ncfile, ncid%station_soiltype)

    ! create variable info fields:
    ! volume variables
    var_name_dims = (/ ncid%charid, ncid%nvars /)
    CALL nf(nf_def_var(ncfile, "var_name", NF_CHAR, 2, var_name_dims(:), &
      &                ncid%var_name))
    CALL nf_add_descr("Variable name (character string)", ncfile, ncid%var_name)
    CALL nf(nf_def_var(ncfile, "var_long_name", NF_CHAR, 2, var_name_dims(:), &
      &                ncid%var_longname))
    CALL nf_add_descr("Variable name (long, character string)", ncfile, ncid%var_longname)
    CALL nf(nf_def_var(ncfile, "var_unit", NF_CHAR, 2, var_name_dims(:), &
      &                ncid%var_unit))
    CALL nf_add_descr("Variable unit (character string)", ncfile, ncid%var_unit)
    CALL nf(nf_def_var(ncfile, "var_group_id", NF_INT, 1, ncid%nvars, &
      &                ncid%var_group_id))
    CALL nf_add_descr("Variable group ID", ncfile, ncid%var_group_id)
    CALL nf(nf_def_var(ncfile, "var_nlevs", NF_INT, 1, ncid%nvars, &
      &                ncid%var_nlevs))
    CALL nf_add_descr("No. of levels for volume variable", ncfile, ncid%var_nlevs)
    var_level_dims = (/ ncid%max_nlevs, ncid%nvars /)
    CALL nf(nf_def_var(ncfile, "var_levels", NF_DOUBLE, 2, var_level_dims(:), &
      &                ncid%var_levels))
    CALL nf_add_descr("Volume variable levels (indices)", ncfile, ncid%var_levels)

    ! surface variables:
    IF (meteogram_data%nsfcvars > 0) THEN
      var_name_dims = (/ ncid%charid, ncid%nsfcvars /)
      CALL nf(nf_def_var(ncfile, "sfcvar_name", NF_CHAR, 2, var_name_dims(:), &
        &                ncid%sfcvar_name))
      CALL nf_add_descr("Surface variable name (character string)", ncfile, ncid%sfcvar_name)
      CALL nf(nf_def_var(ncfile, "sfcvar_long_name", NF_CHAR, 2, var_name_dims(:), &
        &                ncid%sfcvar_longname))
      CALL nf_add_descr("Surface variable name (long, character string)", &
        &               ncfile, ncid%sfcvar_longname)
      CALL nf(nf_def_var(ncfile, "sfcvar_unit", NF_CHAR, 2, var_name_dims(:), &
        &                ncid%sfcvar_unit))
      CALL nf_add_descr("Surface variable unit (character string)", ncfile, ncid%sfcvar_unit)
      CALL nf(nf_def_var(ncfile, "sfcvar_group_id", NF_INT, 1, ncid%nsfcvars, &
        &                ncid%sfcvar_group_id))
      CALL nf_add_descr("Surface variable group ID", ncfile, ncid%sfcvar_group_id)
    END IF

    ! create variables for time slice info:
    CALL nf(nf_def_var(ncfile, "time_step", NF_INT, 1, ncid%timeid, &
      &                ncid%time_step))
    CALL nf_add_descr("Time step indices", ncfile, ncid%time_step)
    time_string_dims = (/ ncid%charid, ncid%timeid /)
    CALL nf(nf_def_var(ncfile, "date", NF_CHAR, 2, time_string_dims(:), &
      &                ncid%dateid))
    CALL nf_add_descr("Sample dates (character string)", ncfile, ncid%dateid)

    ! height levels
    height_level_dims = (/ ncid%nstations, ncid%nvars, ncid%max_nlevs /)
    CALL nf(nf_def_var(ncfile, "heights", NF_DOUBLE, 3, height_level_dims(:), &
      &                ncid%var_heights))
    CALL nf_add_descr("level heights for volume variables", ncfile, ncid%var_heights)

    ! add value buffer for volume variables:
    var_dims = (/ ncid%nstations, ncid%nvars, ncid%max_nlevs, ncid%timeid /)
    CALL nf(nf_def_var(ncfile, "values", NF_DOUBLE, 4, var_dims(:), &
      &                ncid%var_values))
    CALL nf_add_descr("value buffer for volume variables", ncfile, ncid%var_values)
    ! add value buffer for surface variables:
    IF (meteogram_data%nsfcvars > 0) THEN
      sfcvar_dims = (/ ncid%nstations, ncid%nsfcvars, ncid%timeid /)
      CALL nf(nf_def_var(ncfile, "sfcvalues", NF_DOUBLE, 3, sfcvar_dims(:), &
        &                ncid%sfcvar_values))
      CALL nf_add_descr("value buffer for surface variables", ncfile, ncid%sfcvar_values)
    END IF

    ! ----------------------
    ! End of definition mode
    CALL nf(nf_enddef(ncfile))

    DO ivar=1,nvars
      CALL nf(nf_put_vara_text(ncfile, ncid%var_name, (/ 1, ivar /), &
        &        (/ LEN(TRIM(meteogram_data%var_info(ivar)%cf%standard_name)), 1 /), &
        &        TRIM(meteogram_data%var_info(ivar)%cf%standard_name)))
      CALL nf(nf_put_vara_text(ncfile, ncid%var_longname, (/ 1, ivar /), &
        &        (/ LEN(TRIM(meteogram_data%var_info(ivar)%cf%long_name)), 1 /), &
        &        TRIM(meteogram_data%var_info(ivar)%cf%long_name)))
      CALL nf(nf_put_vara_text(ncfile, ncid%var_unit, (/ 1, ivar /), &
        &        (/ LEN(TRIM(meteogram_data%var_info(ivar)%cf%units)), 1 /), &
        &        TRIM(meteogram_data%var_info(ivar)%cf%units)))
      CALL nf(nf_put_vara_int(ncfile, ncid%var_group_id, ivar, 1, &
        &        meteogram_data%var_info(ivar)%igroup_id))
      CALL nf(nf_put_vara_int(ncfile, ncid%var_nlevs, ivar, 1, &
        &        meteogram_data%var_info(ivar)%nlevs))
      istart2 = (/ 1, ivar /)
      icount2 = (/ meteogram_data%var_info(ivar)%nlevs, 1 /)
      CALL nf(nf_put_vara_int(ncfile, ncid%var_levels, istart2, icount2, &
        &                     meteogram_data%var_info(ivar)%levels(:)))
    END DO

    DO ivar=1,nsfcvars
      CALL nf(nf_put_vara_text(ncfile, ncid%sfcvar_name, (/ 1, ivar /), &
        &        (/ LEN(TRIM(meteogram_data%sfc_var_info(ivar)%cf%standard_name)), 1 /), &
        &        TRIM(meteogram_data%sfc_var_info(ivar)%cf%standard_name)))
      CALL nf(nf_put_vara_text(ncfile, ncid%sfcvar_longname, (/ 1, ivar /), &
        &        (/ LEN(TRIM(meteogram_data%sfc_var_info(ivar)%cf%long_name)), 1 /), &
        &        TRIM(meteogram_data%sfc_var_info(ivar)%cf%long_name)))
      CALL nf(nf_put_vara_text(ncfile, ncid%sfcvar_unit, (/ 1, ivar /), &
        &        (/ LEN(TRIM(meteogram_data%sfc_var_info(ivar)%cf%units)), 1 /), &
        &        TRIM(meteogram_data%sfc_var_info(ivar)%cf%units)))
      CALL nf(nf_put_vara_int(ncfile, ncid%sfcvar_group_id, ivar, 1, &
        &        meteogram_data%sfc_var_info(ivar)%igroup_id))
    END DO

    istation = 1
    DO jb=1,meteogram_data%nblks
      i_startidx = 1
      i_endidx   = nproma
      IF (jb == meteogram_data%nblks) i_endidx = meteogram_data%npromz

      DO jc=i_startidx,i_endidx
        this_station => meteogram_output_config%station_list(           &
          &               meteogram_data%station(jc,jb)%station_idx(1), &
          &               meteogram_data%station(jc,jb)%station_idx(2))
        CALL nf(nf_put_vara_text(ncfile, ncid%station_name, (/ 1, istation /), &
          &                      (/ LEN(TRIM(this_station%zname)), 1 /), &
          &                      TRIM(this_station%zname)))
        CALL nf(nf_put_vara_double(ncfile, ncid%station_lon, istation, 1, &
          &                        this_station%location%lon))
        CALL nf(nf_put_vara_double(ncfile, ncid%station_lat, istation, 1, &
          &                        this_station%location%lat))
        CALL nf(nf_put_vara_int(ncfile, ncid%station_idx, istation, 1, &
          &                     meteogram_data%station(jc,jb)%tri_idx(1)))
        CALL nf(nf_put_vara_int(ncfile, ncid%station_blk, istation, 1, &
          &                     meteogram_data%station(jc,jb)%tri_idx(2)))
        CALL nf(nf_put_vara_double(ncfile, ncid%station_hsurf, istation, 1, &
          &                        meteogram_data%station(jc,jb)%hsurf))
        CALL nf(nf_put_vara_double(ncfile, ncid%station_frland, istation, 1, &
          &                        meteogram_data%station(jc,jb)%frland))
        CALL nf(nf_put_vara_double(ncfile, ncid%station_fc, istation, 1, &
          &                        meteogram_data%station(jc,jb)%fc))
        CALL nf(nf_put_vara_int(ncfile, ncid%station_soiltype, istation, 1, &
          &                     meteogram_data%station(jc,jb)%soiltype))

        ! model level heights
        DO ivar=1,nvars
          nlevs = meteogram_data%var_info(ivar)%nlevs
          CALL nf(nf_put_vara_double(ncfile, ncid%var_heights,     &
            &                        (/ istation, ivar,     1 /),  &
            &                        (/        1,    1, nlevs /),  &
            &                        meteogram_data%station(jc,jb)%var(ivar)%heights(1:nlevs)))
        END DO

        istation = istation + 1 

      END DO
    END DO

  END SUBROUTINE meteogram_open_file


  !>
  !! The IO PE writes the global meteogram buffer to the output
  !! file. Afterwards, the global meteogram buffer is cleared.
  !!
  !! For gathered NetCDF output, this is a collective operation,
  !! otherwise this is a local operation for the IO PE.
  !!
  !! @par Revision History
  !! Initial implementation  by  F. Prill, DWD (2011-08-22)
  !!
  SUBROUTINE meteogram_flush_file(meteogram_output_config, jg)
    ! station data from namelist
    TYPE(t_meteogram_output_config), TARGET, INTENT(IN) :: meteogram_output_config
    INTEGER, INTENT(IN)         :: jg       !< patch index
    ! local variables:
    CHARACTER(*), PARAMETER :: routine = TRIM("mo_meteogram_output:meteogram_flush_file")
    INTEGER                     :: ncfile,  totaltime, itime, istation, ivar, &
      &                            jb, jc, i_startidx, i_endidx, nlevs,       &
      &                            nvars, nsfcvars
    TYPE(t_meteogram_data), POINTER :: meteogram_data
    TYPE(t_ncid)          , POINTER :: ncid
    INTEGER                         :: istart4(4), icount4(4)
  
    ncid => ncid_list(jg)
    ncfile = meteogram_file_info(jg)%file_id

    ! In "non-distributed" mode, station data is gathered by PE #0
    ! which writes a single file:
    IF (meteogram_output_config%ldistributed) THEN
      meteogram_data => meteogram_local_data(jg)
    ELSE
      CALL meteogram_collect_buffers(meteogram_output_config, jg, io_rank)
      meteogram_data => meteogram_global_data(jg)
    END IF

    ! skip routine, if this PE has nothing to do...
    IF (meteogram_data%nstations == 0) RETURN
    IF  (.NOT. meteogram_output_config%ldistributed .AND.  &
      & (.NOT. l_is_receiver)) RETURN

    nvars    = meteogram_data%nvars
    nsfcvars = meteogram_data%nsfcvars

    IF (dbg_level > 0) THEN
      WRITE(message_text,*) "Meteogram"
      CALL message(routine, TRIM(message_text))    
    END IF

    ! inquire about current number of records in file:
    CALL nf(nf_inq_dimlen(ncfile, ncid%timeid, totaltime))

    ! write time stamp info:
    DO itime=1,meteogram_local_data(jg)%icurrent

      CALL nf(nf_put_vara_text(ncfile, ncid%dateid, (/ 1, totaltime+itime /), &
        &                      (/ LEN(TRIM(meteogram_data%time_stamp(itime)%zdate)), 1 /), &
        &                      TRIM(meteogram_data%time_stamp(itime)%zdate)))
      CALL nf(nf_put_vara_int(ncfile, ncid%time_step, totaltime+itime, 1, &
        &                     meteogram_data%time_stamp(itime)%istep))

      ! write meteogram buffer:
      istation = 1
      DO jb=1,meteogram_data%nblks
        i_startidx = 1
        i_endidx   = nproma
        IF (jb == meteogram_data%nblks)  &
          &  i_endidx = meteogram_data%npromz

        DO jc=i_startidx,i_endidx

          ! volume variables:
          DO ivar=1,nvars
            nlevs = meteogram_data%var_info(ivar)%nlevs
            istart4 = (/ istation, ivar, 1, totaltime+itime /)
            icount4 = (/ 1, 1, nlevs, 1 /)
            CALL nf(nf_put_vara_double(ncfile, ncid%var_values,                   &
              &                        istart4, icount4,                          &
              &                        meteogram_data%station(jc,jb)%var(ivar)%values(1:nlevs, &
              &                        itime)))
          END DO
          ! surface variables:
          DO ivar=1,nsfcvars
            CALL nf(nf_put_vara_double(ncfile, ncid%sfcvar_values,                &
              &                        (/ istation, ivar, totaltime+itime /),     &
              &                        (/ 1, 1, 1 /),                             &
              &                        meteogram_data%station(jc,jb)%sfc_var(ivar)%values(itime)))
          END DO

          istation = istation + 1
        END DO
      END DO
    END DO

    ! finally, reset buffer counter for new data
    meteogram_local_data(jg)%icurrent = 0

  END SUBROUTINE meteogram_flush_file


  !>
  !! The IO PE closes the meteogram output file.
  !! For gathered NetCDF output, this is a collective operation,
  !! otherwise this is a local operation for the IO PE.
  !!
  !! @par Revision History
  !! Initial implementation  by  F. Prill, DWD (2011-08-22)
  !!
  SUBROUTINE meteogram_close_file(meteogram_output_config, jg)
    ! station data from namelist
    TYPE(t_meteogram_output_config), TARGET, INTENT(IN) :: meteogram_output_config
    INTEGER, INTENT(IN)  :: jg    !< patch index

    ! write remaining buffers:
    CALL meteogram_flush_file(meteogram_output_config, jg)

    ! Close NetCDF file
    ! skip routine, if this PE has nothing to do...
    IF (meteogram_output_config%ldistributed) THEN
      IF (meteogram_local_data(jg)%nstations > 0) &
        CALL nf(nf_close(meteogram_file_info(jg)%file_id))
    ELSE
      IF (io_rank == get_my_mpi_all_id()) &
        CALL nf(nf_close(meteogram_file_info(jg)%file_id))
    END IF
  END SUBROUTINE meteogram_close_file


  !>
  !! @return file name of meteogram file.
  !! This is a local operation.
  !!
  !! @par Revision History
  !! Initial implementation  by  F. Prill, DWD (2011-08-22)
  !!
  SUBROUTINE meteogram_create_filename (meteogram_output_config, jg)
    
    ! station data from namelist
    TYPE(t_meteogram_output_config), TARGET, INTENT(IN) :: meteogram_output_config
    ! patch index
    INTEGER,                             INTENT(IN) :: jg
    ! Local variables
    INTEGER :: my_id

    my_id = get_my_mpi_all_id()

    SELECT CASE (meteogram_output_config%ftype)
    CASE (FTYPE_NETCDF)
      IF (meteogram_output_config%ldistributed) THEN
        WRITE (meteogram_file_info(jg)%zname,'(a,i3.3,a,i3.3,a)') "PE", my_id, "_patch", jg, ".nc"
      ELSE
        WRITE (meteogram_file_info(jg)%zname,'(a,i3.3,a)') "patch", jg, ".nc"
      END IF
    END SELECT
    meteogram_file_info(jg)%zname = &
      &  TRIM(meteogram_output_config%zprefix)//TRIM(meteogram_file_info(jg)%zname)
  END SUBROUTINE meteogram_create_filename


  !>
  !!  Help functions for NetCDF I/O
  !!
  !!  Checks the return value of a NetCDF function and exits in case of
  !!  an error
  SUBROUTINE nf(status)

    INTEGER, INTENT(in) :: status !< NetCDF error code

    IF (status /= nf_noerr) THEN
      CALL finish(modname, 'NetCDF Error: '//nf_strerror(status))
    ENDIF

  END SUBROUTINE nf


  !>
  !!  Help functions for NetCDF I/O
  !!
  !!  Adds a string attribute containing variable description.
  SUBROUTINE nf_add_descr(description_str, ncfile, var_id)

    CHARACTER(LEN=*), INTENT(in) :: description_str
    INTEGER         , INTENT(in) :: ncfile, var_id
    CHARACTER(LEN=*), PARAMETER  :: descr_label = "description"

    CALL nf(nf_put_att_text(ncfile, var_id, descr_label, &
      &         LEN(TRIM(description_str)), TRIM(description_str)))

  END SUBROUTINE nf_add_descr


  !>
  !!  Utility function (3d formulation of generic interface).
  !!
  !!  Registers a new atmospheric (volume) variable.
  SUBROUTINE add_atmo_var_3d(igroup_id, zname, zunit, zlong_name, jg, source)
    CHARACTER(LEN=*),  INTENT(IN)    :: zname, zunit, zlong_name
    INTEGER,           INTENT(IN)    :: igroup_id, jg
    REAL(wp), TARGET,  INTENT(IN)    :: source(:,:,:)   !< source array
    ! Local variables
    CHARACTER(*), PARAMETER :: routine = TRIM("mo_meteogram_output:add_atmo_var")
    TYPE(t_meteogram_data), POINTER :: meteogram_data
    INTEGER                         :: ierrstat, nlev, ilev, ivar

    nlev = SIZE(source, 2) ! get level no from array dimensions
    meteogram_data => meteogram_local_data(jg)

    ! create new variable index
    var_list(jg)%no_atmo_vars = var_list(jg)%no_atmo_vars + 1
    ivar = var_list(jg)%no_atmo_vars

    ! create meteogram data structure
    meteogram_data%var_info(ivar)%cf%standard_name = TRIM(zname)
    meteogram_data%var_info(ivar)%cf%long_name     = TRIM(zlong_name)
    meteogram_data%var_info(ivar)%cf%units         = TRIM(zunit)
    meteogram_data%var_info(ivar)%igroup_id        = igroup_id
    meteogram_data%var_info(ivar)%nlevs            = nlev
    meteogram_data%var_info(ivar)%p_source_3d => source
    meteogram_data%max_nlevs = MAX(meteogram_data%max_nlevs,  meteogram_data%var_info(ivar)%nlevs)
    ALLOCATE(meteogram_data%var_info(ivar)%levels(nlev), stat=ierrstat)
    IF (ierrstat /= SUCCESS) THEN
      CALL finish (routine, 'ALLOCATE of meteogram data structures failed')
    ENDIF
    meteogram_data%var_info(ivar)%levels = (/ (ilev, ilev=1,nlev) /)
  END SUBROUTINE add_atmo_var_3d


  !>
  !!  Utility function (4d formulation of generic interface).
  !!  Adds the 4d var as separate 3d var slices.
  !!
  !!  Registers a new atmospheric (volume) variable.
  SUBROUTINE add_atmo_var_4d(igroup_id, zname, zunit, zlong_name, jg, source)
    CHARACTER(LEN=*),  INTENT(IN)    :: zname, zunit, zlong_name
    INTEGER,           INTENT(IN)    :: igroup_id, jg
    REAL(wp), TARGET,  INTENT(IN)    :: source(:,:,:,:)   !< source array
    ! Local variables
    CHARACTER(*), PARAMETER :: routine = TRIM("mo_meteogram_output:add_atmo_var")
    TYPE(t_meteogram_data), POINTER :: meteogram_data
    INTEGER                         :: ierrstat, nlev, ilev, ivar, nidx, isource_idx

    nlev = SIZE(source, 2) ! get level no from array dimensions
    nidx = SIZE(source, 4) ! get number of 4d var indices (e.g. tile number)
    meteogram_data => meteogram_local_data(jg)

    ! add the 4d var as separate 3d var slices
    DO isource_idx=1,nidx
      ! create new variable index
      var_list(jg)%no_atmo_vars = var_list(jg)%no_atmo_vars + 1
      ivar = var_list(jg)%no_atmo_vars
      
      ! create meteogram data structure
      meteogram_data%var_info(ivar)%cf%standard_name = TRIM(zname)//"_"//int2string(isource_idx)
      meteogram_data%var_info(ivar)%cf%long_name     = TRIM(zlong_name)
      meteogram_data%var_info(ivar)%cf%units         = TRIM(zunit)
      meteogram_data%var_info(ivar)%igroup_id        = igroup_id
      meteogram_data%var_info(ivar)%nlevs            = nlev
      meteogram_data%var_info(ivar)%p_source_3d => source(:,:,:,isource_idx)
      
      meteogram_data%max_nlevs = MAX(meteogram_data%max_nlevs,  &
        &                            meteogram_data%var_info(ivar)%nlevs)
      ALLOCATE(meteogram_data%var_info(ivar)%levels(nlev), stat=ierrstat)
      IF (ierrstat /= SUCCESS) THEN
        CALL finish (routine, 'ALLOCATE of meteogram data structures failed')
      ENDIF
      meteogram_data%var_info(ivar)%levels = (/ (ilev, ilev=1,nlev) /)
    END DO
  END SUBROUTINE add_atmo_var_4d


  !>
  !!  Utility function.
  !!
  !!  Registers a new surface variable.
  SUBROUTINE add_sfc_var_2d(igroup_id, zname, zunit, zlong_name, jg, source)
    CHARACTER(LEN=*),  INTENT(IN)    :: zname, zunit, zlong_name
    INTEGER,           INTENT(IN)    :: igroup_id, jg
    REAL(wp), TARGET,  INTENT(IN)    :: source(:,:)   !< source array
    ! Local variables
    TYPE(t_meteogram_data), POINTER :: meteogram_data
    INTEGER                 :: ivar
 
    meteogram_data => meteogram_local_data(jg)
    ! create new variable index
    var_list(jg)%no_sfc_vars = var_list(jg)%no_sfc_vars + 1
    ivar = var_list(jg)%no_sfc_vars
    ! create meteogram data structure
    meteogram_data%sfc_var_info(ivar)%cf%standard_name = TRIM(zname)
    meteogram_data%sfc_var_info(ivar)%cf%long_name     = TRIM(zlong_name)
    meteogram_data%sfc_var_info(ivar)%cf%units         = TRIM(zunit)
    meteogram_data%sfc_var_info(ivar)%igroup_id        = igroup_id
    meteogram_data%sfc_var_info(ivar)%p_source  => source
  END SUBROUTINE add_sfc_var_2d


  !>
  !!  Utility function.
  !!
  !!  Registers a new surface variable.
  SUBROUTINE add_sfc_var_3d(igroup_id, zname, zunit, zlong_name, jg, source)
    CHARACTER(LEN=*),  INTENT(IN)    :: zname, zunit, zlong_name
    INTEGER,           INTENT(IN)    :: igroup_id, jg
    REAL(wp), TARGET,  INTENT(IN)    :: source(:,:,:)   !< source array
    ! Local variables
    TYPE(t_meteogram_data), POINTER :: meteogram_data
    INTEGER                 :: ivar, nidx, isource_idx
 
    meteogram_data => meteogram_local_data(jg)
    nidx = SIZE(source, 3) ! get number of 4d var indices (e.g. tile number)

    ! add the 3d var as separate 2d var slices
    DO isource_idx=1,nidx
      ! create new variable index
      var_list(jg)%no_sfc_vars = var_list(jg)%no_sfc_vars + 1
      ivar = var_list(jg)%no_sfc_vars
      ! create meteogram data structure
      meteogram_data%sfc_var_info(ivar)%cf%standard_name = &
        &   TRIM(zname)//"_"//int2string(isource_idx)
      meteogram_data%sfc_var_info(ivar)%cf%long_name     = TRIM(zlong_name)
      meteogram_data%sfc_var_info(ivar)%cf%units         = TRIM(zunit)
      meteogram_data%sfc_var_info(ivar)%igroup_id        = igroup_id
      meteogram_data%sfc_var_info(ivar)%p_source  => source(:,:,isource_idx)
    END DO
  END SUBROUTINE add_sfc_var_3d

END MODULE mo_meteogram_output
