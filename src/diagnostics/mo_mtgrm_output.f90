!>
!! Data structures and subroutines for meteogram output.
!!
!! The sampling intervals for meteogram data are independent 
!! from global output steps. Values are buffered in memory until
!! the next field output is invoked.
!! Before each write operation, data is gathered from all working
!! PEs by the IO PE and written to a NetCDF file.
!!
!! Known limitations:
!! - So far, only NetCDF file format is supported, which is written by
!!   IO PE #0.
!! - ASCII output data (similar to COSMO) must be generated in a
!!   post-processing step.
!! - In case of an application crash, latest meteogram data, which has
!!   not yet been written to hard disk, may be lost.
!!
!! How to add new variables for sampling:
!! --------------------------------------
!!
!! a) Volume variables:
!!    - Add a new index variable "MYVAR" in data structure "TYPE t_var"
!!    - In SR "mtgrm_init": 
!!      Increase total number of sampling variables
!!        "mtgrm_data%nvars"
!!      Insert
!!        CALL add_atmo_var(VAR%MYVAR, VAR_GROUP_ATMO, "myvarname", "myvarunit", jg, nlev)
!!    - In SR "mtgrm_sample_vars": Insert sampling of values, for example
!!        DO ilev=1,mtgrm_data%var_info(VAR%MYVAR)%nlevs
!!          mtgrm_data%station(jc,jb)%var(VAR%MYVAR)%values(ilev, i_tstep) =  &
!!          &  prog%<myvar>(iidx, mtgrm_data%var_info(VAR%MYVAR)%levels(ilev), iblk)
!!        END DO
!!
!! b) Surface variables:
!!    - Add a new index variable "MYVAR" in data structure "TYPE t_var"
!!    - In SR "mtgrm_init": 
!!      Increase total number of sampling surface variables
!!        "mtgrm_data%nsfcvars"
!!      Insert
!!        CALL add_sfc_var(VAR%MYVAR, VAR_GROUP_SFC, "myvarname", "myvarunit", jg)
!!    - In SR "mtgrm_sample_vars": Insert sampling of values, for example
!!          mtgrm_data%station(jc,jb)%sfcvar(VAR%MYVAR)%values(i_tstep) =  &
!!          &  prog%<myvar>(iidx, iblk)
!!
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
!! TODO[FP] : Extend list of meteogram variables
!! TODO[FP] : Subroutine "mtgrm_collect_buffers" contains
!!            MPI calls which should be moved to MODULE mo_mpi.

MODULE mo_mtgrm_output

  USE mo_kind,                  ONLY: wp
  USE mo_datetime,              ONLY: t_datetime, iso8601
  USE mo_exception,             ONLY: message, message_text, finish
  USE mo_mpi,                   ONLY: p_n_work, my_process_is_stdio, &
    &                                 get_my_mpi_all_id,             &
    ! TODO[FP] : When encapsulating MPI_PACK calls into mo_mpi,
    !            the following USEs are no longer required
    &                                 get_my_mpi_all_id, &
    &                                 get_my_mpi_communicator, p_real_dp, &
    &                                 p_int, get_mpi_all_workroot_id,     &
    &                                 p_real_dp_byte
#ifndef NOMPI
  USE mo_mpi,                   ONLY: MPI_STATUS_SIZE
#endif
  USE mo_model_domain,          ONLY: t_patch
  USE mo_parallel_config,       ONLY: nproma
  USE mo_impl_constants,        ONLY: inwp, max_dom, SUCCESS
  USE mo_communication,         ONLY: idx_1d, blk_no, idx_no
  USE mo_ext_data,              ONLY: t_external_data
  USE mo_nonhydro_state,        ONLY: t_nh_state, t_nh_prog, t_nh_diag
  ! TODO[FP] : When using an already built GNAT, not all of the
  ! following USEs will be necessary:
  USE mo_gnat_gridsearch,       ONLY: gnat_init_grid, gnat_destroy, gnat_tree,&
    &                                 gnat_query_containing_triangles,        &
    &                                 gnat_merge_distributed_queries, gk
  USE mo_dynamics_config,       ONLY: nnow
  USE mo_io_config,             ONLY: lwrite_extra, inextra_2d, inextra_3d
  USE mo_util_string,           ONLY: int2string
  USE mo_mtgrm_config,          ONLY: t_mtgrm_output_config, t_station_list, &
    &                                 FTYPE_NETCDF, MAX_NAME_LENGTH, MAX_NUM_STATIONS
  ! TODO[FP] : When encapsulating MPI_PACK calls into mo_mpi,
  !            the following USE is no longer required:
#ifndef NOMPI
  USE MPI
#endif
  
  IMPLICIT NONE
  
  PRIVATE
  CHARACTER(LEN=*), PARAMETER :: modname     = 'mo_mtgrm_output'
  INTEGER,          PARAMETER :: dbg_level   = 0

  INCLUDE 'netcdf.inc'

  ! IO routines.
  ! called collectively, though non-IO PEs are occupied
  ! only for the case of distributed write mode.
  PUBLIC ::  mtgrm_init  
  PUBLIC ::  mtgrm_is_sample_step
  PUBLIC ::  mtgrm_sample_vars
  PUBLIC ::  mtgrm_finalize  
  PUBLIC ::  mtgrm_flush_file

  INTEGER, PARAMETER :: MAX_TIME_STAMPS      =  100  !<  max. number of time stamps
  INTEGER, PARAMETER :: MAX_EXTRA2D          =    4  !<  max. number of extra 2D variables (inextra_2d)
  INTEGER, PARAMETER :: MAX_EXTRA3D          =    4  !<  max. number of extra 3D variables (inextra_3d)

  ! arbitrarily chosen value for buffer size (somewhat large for safety reasons)
  INTEGER, PARAMETER :: MAX_HEADER_SIZE      =  128 ! *p_real_dp_byte

  INTEGER, PARAMETER :: VAR_GROUP_ATMO       =    1  
  INTEGER, PARAMETER :: VAR_GROUP_SURFACE    =    2  


  !>
  !! Storage for information on a single variable
  !!
  TYPE t_var_info
    CHARACTER(len=MAX_NAME_LENGTH) :: zname, zunit  !< variable name, unit
    INTEGER                        :: igroup_id     !< variable group (surface vars, soil temperatures, ...)
    INTEGER                        :: nlevs         !< number of levels for this variable
    INTEGER, POINTER               :: levels(:)     !< level indices (1:nlevs)
  END TYPE t_var_info

  !>
  !! Storage for information on a single surface variable
  !!
  TYPE t_sfc_var_info
    CHARACTER(len=MAX_NAME_LENGTH) :: zname, zunit  !< variable name, unit
    INTEGER                        :: igroup_id     !< variable group (surface vars, soil temperatures, ...)
  END TYPE t_sfc_var_info

  !>
  !! Value buffer for a single variable of a station.
  !!
  TYPE t_var_buffer
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
  !! t_mtgrm_data.
  !!
  !! Note: This info is different for different patches.
  !!
  TYPE t_mtgrm_station
    ! Meteogram header (information on location, ...)
    INTEGER                       :: station_idx(2)   !< (idx,block) of station specification
    INTEGER                       :: tri_idx(2)       !< triangle index (global idx,block)
    INTEGER                       :: tri_idx_local(2) !< triangle index (idx,block)
    INTEGER                       :: owner      !< proc ID where station is located.    
    REAL(wp)                      :: hsurf      !< surface height
    REAL(wp)                      :: frland     !< fraction of land
    REAL(wp)                      :: fc         !< Coriolis parameter
    INTEGER                       :: soiltype   !< soil type

    ! Buffer for currently stored meteogram values.
    TYPE(t_var_buffer),     POINTER :: var(:)       !< sampled data (1:nvars)
    TYPE(t_sfc_var_buffer), POINTER :: sfc_var(:)   !< sampled data (1:nsfcvars)
  END TYPE t_mtgrm_station

  !>
  !! Storage for information on the set of collected variables for
  !! several stations.
  !!
  TYPE t_mtgrm_data
    ! variable info:
    INTEGER                        :: nvars, nsfcvars !< number of sampled variables and surface variables
    INTEGER                        :: max_nlevs       !< maximum no. of levels for variables
    TYPE(t_var_info),     POINTER  :: var_info(:)     !< info for each variable (1:nvars)
    TYPE(t_sfc_var_info), POINTER  :: sfc_var_info(:) !< info for each surface variable (1:nsfcvars)
    ! time stamp info:
    INTEGER                        :: icurrent        !< current time stamp index
    TYPE(t_time_stamp), POINTER    :: time_stamp(:)   !< info on sample times
    ! value buffers:
    TYPE(t_mtgrm_station), POINTER :: station(:,:) !< meteogram data and meta info for each station (idx,blk).
    INTEGER                        :: nstations, nblks, npromz
  END TYPE t_mtgrm_data
  
  !>
  !! Data structure specifying output file for meteogram data.
  !!
  TYPE t_mtgrm_file
    INTEGER                        :: ftype   !< file type (NetCDF, ...)
    CHARACTER(len=MAX_NAME_LENGTH) :: zname   !< file name string
    INTEGER                        :: file_id !< meteogram file ID
  END TYPE t_mtgrm_file

  !>
  !! Data structure specifying NetCDF IDs
  !!
  TYPE t_ncid
    INTEGER  :: nstations, nvars, charid, station_name, station_lat, station_lon,         &
      &         station_idx, station_blk, station_hsurf, station_frland, station_fc,      &
      &         station_soiltype, nsfcvars, var_name, var_unit, sfcvar_name, sfcvar_unit, &
      &         var_group_id, sfcvar_group_id, var_nlevs, max_nlevs, var_levels, timeid,  &
      &         time_step, dateid, var_values, sfcvar_values
  END TYPE t_ncid

  !>
  !! Data structure containing internal indices for variables
  !!
  TYPE t_var
    INTEGER :: no_atmo_vars       !< number of atmo variables declared so far
    INTEGER :: no_sfc_vars        !< number of surface variables declared so far

    ! The following variables are declared for convenience. Their values
    ! are set when defining the meteogram variables. Elsewhere they
    ! should be treated as constants.
    INTEGER :: P, T, RHO, EXNER, THETAV, RHOTHETAV, &
      &        U, V, W, EXTRA2D(MAX_EXTRA2D), EXTRA3D(MAX_EXTRA3D), &
      &        P_SFC
  END TYPE t_var

  ! -------------------------------------------------------------------------------------------

  !! -- module data: --
  TYPE(t_mtgrm_data), SAVE, TARGET      :: mtgrm_local_data(1:max_dom)   !< meteogram data local to this PE
  TYPE(t_mtgrm_data), SAVE, TARGET      :: mtgrm_global_data(1:max_dom)  !< collected buffers (on IO PE)
  TYPE(t_mtgrm_file), SAVE              :: mtgrm_file_info(1:max_dom)    !< meteogram file handle etc.

  TYPE(t_ncid),       SAVE, TARGET      :: ncid_list(1:max_dom)          !< NetCDF dimension IDs
  TYPE(t_var),        SAVE, TARGET      :: var_list(1:max_dom)           !< internal indices of variables
  
  !! -- data for distributed meteogram sampling (MPI) --
  CHARACTER,          SAVE, ALLOCATABLE :: msg_buffer(:,:)               !< MPI buffer for station data
  INTEGER,            SAVE              :: max_buf_size                  !< max buffer size for MPI messages
  LOGICAL,            SAVE              :: l_is_sender, l_is_receiver
  INTEGER,            SAVE              :: io_rank                       !< rank of PE which gathers data
  INTEGER,            SAVE              :: owner(MAX_NUM_STATIONS)       !< rank of sender PE for each station

CONTAINS

  !>
  !! Initialize meteogram data buffer, allocating storage.
  !! This is a collective operation.
  !!
  !! @par Revision History
  !! Initial implementation  by  F. Prill, DWD (2011-08-22)
  !!
  SUBROUTINE mtgrm_init(mtgrm_output_config, ptr_patch, ext_data, iforcing, jg)
    ! station data from namelist
    TYPE(t_mtgrm_output_config), TARGET, INTENT(IN) :: mtgrm_output_config
    ! data structure containing grid info:
    TYPE(t_patch),                       INTENT(IN) :: ptr_patch
    ! atmosphere external data
    TYPE(t_external_data),               INTENT(IN) :: ext_data
    ! parameterized forcing (right hand side) of dynamics, affects
    ! topography specification, see "mo_extpar_config"
    INTEGER,                             INTENT(IN) :: iforcing
    ! patch index
    INTEGER,                             INTENT(IN) :: jg

    ! local variables:
    CHARACTER(*), PARAMETER :: routine = TRIM("mo_mtgrm_output:mtgrm_init")

    INTEGER      :: ithis_nlocal_pts, nblks, npromz,  &
      &             nstations, ierrstat,              &
      &             jb, jc, glb_index, i_startidx,    &
      &             i_endidx, jc_station, jb_station, &
      &             istation, my_id, nvars, nlevs,    &
      &             ivar, nsfcvars, iextra,           &
      &             nblks_global, npromz_global
    REAL(gk)     :: in_points(nproma,mtgrm_output_config%nblks,2) !< geographical locations
    REAL(gk)     :: min_dist(nproma,mtgrm_output_config%nblks)    !< minimal distance
    ! list of triangles containing lon-lat grid points (first dim: index and block)
    INTEGER      :: tri_idx(2,nproma,mtgrm_output_config%nblks)
    INTEGER      :: nlocal_pts(p_n_work)
    LOGICAL      :: l_is_io_pe
    REAL(wp)     :: pi_180
    INTEGER      :: max_var_size, max_sfcvar_size
    TYPE(t_mtgrm_data)   , POINTER :: mtgrm_data
    TYPE(t_var)          , POINTER :: VAR
    TYPE(t_mtgrm_station), POINTER :: p_station

    pi_180 = ATAN(1._wp)/45._wp

    mtgrm_data => mtgrm_local_data(jg)

    VAR => var_list(jg)    
    VAR%no_atmo_vars = 0
    VAR%no_sfc_vars  = 0

    ! ------------------------------------------------------------
    ! Distribute stations, determine number of stations located on
    ! this PE:
    ! ------------------------------------------------------------

    ! build an array of geographical coordinates from station list:
    ! in_points(...)
    nstations = mtgrm_output_config%nstations
    nblks     = mtgrm_output_config%nblks
    npromz    = mtgrm_output_config%npromz   

    DO jb=1,nblks
      i_startidx = 1
      i_endidx   = nproma
      IF (jb == nblks) i_endidx = npromz

      DO jc=i_startidx,i_endidx    
        in_points(jc,jb,:) = (/ mtgrm_output_config%station_list(jc,jb)%location%lon, &
          &                     mtgrm_output_config%station_list(jc,jb)%location%lat  /) * pi_180
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
    mtgrm_data%nstations = ithis_nlocal_pts
    mtgrm_data%nblks     = nblks
    mtgrm_data%npromz    = npromz
    ! clean up
    CALL gnat_destroy()

    ! ------------------------------------------------------------
    ! Initialize local data structure, fill header
    ! ------------------------------------------------------------

    mtgrm_data%icurrent = 0 ! reset current sample index
    ALLOCATE(mtgrm_data%time_stamp(MAX_TIME_STAMPS), stat=ierrstat)
    IF (ierrstat /= SUCCESS) THEN
      CALL finish (routine, 'ALLOCATE of meteogram data structures failed')
    ENDIF

    ! set up list of variables:
    ! Note: adjust these values when adding new variables:
    mtgrm_data%nsfcvars  = 1
    mtgrm_data%nvars     = 9

    IF (lwrite_extra)  &
      mtgrm_data%nsfcvars = mtgrm_data%nsfcvars + inextra_2d
    nsfcvars = mtgrm_data%nsfcvars
    ALLOCATE(mtgrm_data%sfc_var_info(nsfcvars), stat=ierrstat)
    IF (ierrstat /= SUCCESS) THEN
      CALL finish (routine, 'ALLOCATE of meteogram data structures failed')
    ENDIF

    IF (lwrite_extra)  &
      mtgrm_data%nvars = mtgrm_data%nvars + inextra_3d
    mtgrm_data%max_nlevs = 1
    nvars = mtgrm_data%nvars
    ALLOCATE(mtgrm_data%var_info(nvars), stat=ierrstat)
    IF (ierrstat /= SUCCESS) THEN
      CALL finish (routine, 'ALLOCATE of meteogram data structures failed')
    ENDIF

    ! Variable: Pressure
    CALL add_atmo_var(VAR%P,         VAR_GROUP_ATMO, "P",         "Pa",       jg, ptr_patch%nlev)
    ! Variable: Temperature
    CALL add_atmo_var(VAR%T,         VAR_GROUP_ATMO, "T",         "K",        jg, ptr_patch%nlev)
    ! Variable: Density
    CALL add_atmo_var(VAR%RHO,       VAR_GROUP_ATMO, "RHO",       "kg/m^3",   jg, ptr_patch%nlev)
    ! Variable: Exner pressure
    CALL add_atmo_var(VAR%EXNER,     VAR_GROUP_ATMO, "PEXNER",    "-",        jg, ptr_patch%nlev)
    ! Variable: virtual potential temperature
    CALL add_atmo_var(VAR%THETAV,    VAR_GROUP_ATMO, "THETAV",    "K",        jg, ptr_patch%nlev)
    ! Variable: rho*theta_v
    CALL add_atmo_var(VAR%RHOTHETAV, VAR_GROUP_ATMO, "RHOTHETAV", "K*kg/m^3", jg, ptr_patch%nlev)
    ! Variable: zonal wind
    CALL add_atmo_var(VAR%U,         VAR_GROUP_ATMO, "U",         "m/s",      jg, ptr_patch%nlev)
    ! Variable: meridional wind
    CALL add_atmo_var(VAR%V,         VAR_GROUP_ATMO, "V",         "m/s",      jg, ptr_patch%nlev)
    ! Variable: orthogonal vertical wind
    CALL add_atmo_var(VAR%W,         VAR_GROUP_ATMO, "W",         "m/s",      jg, ptr_patch%nlevp1)

    ! Surface pressure
    CALL add_sfc_var(VAR%P_SFC, VAR_GROUP_SURFACE, "P_SFC", "Pa", jg)

    IF (lwrite_extra) THEN
      ! Variable: Extra 2D
      DO iextra=1,inextra_2d
        CALL add_sfc_var(VAR%EXTRA2D(iextra), VAR_GROUP_SURFACE, "EXTRA2D"//int2string(iextra), &
          &              "-", jg)
      END DO
      ! Variable: Extra 3D
      DO iextra=1,inextra_3d
        CALL add_atmo_var(VAR%EXTRA3D(iextra), VAR_GROUP_ATMO, "EXTRA3D"//int2string(iextra), &
          &               "-", jg, ptr_patch%nlev)
      END DO
    END IF

    ! set up list of stations:
    ALLOCATE(mtgrm_data%station(nproma, nblks), stat=ierrstat)
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
        mtgrm_data%station(jc,jb)%station_idx = (/ jc_station, jb_station /)
        ! set owner ID:
        mtgrm_data%station(jc,jb)%owner = my_id
        ! set local triangle index, block:
        mtgrm_data%station(jc,jb)%tri_idx_local(1:2) = tri_idx(1:2,jc,jb)
        ! translate local index to global index:
        glb_index = ptr_patch%cells%glb_index(idx_1d(tri_idx(1,jc,jb), &
          &                                          tri_idx(2,jc,jb)))
        mtgrm_data%station(jc,jb)%tri_idx(1:2) =  &
          &  (/ idx_no(glb_index), blk_no(glb_index) /)
        ! set Coriolis parameter for station
        mtgrm_data%station(jc,jb)%fc       =  &
          &  ptr_patch%cells%f_c(tri_idx(1,jc,jb), tri_idx(2,jc,jb))
        ! set station information on height, soil type etc.:
        SELECT CASE ( iforcing )
        CASE ( inwp ) ! NWP physics
          mtgrm_data%station(jc,jb)%hsurf    =  &
            &  ext_data%atm%topography_c(tri_idx(1,jc,jb), tri_idx(2,jc,jb))
          mtgrm_data%station(jc,jb)%frland   =  &
            &  ext_data%atm%fr_land(tri_idx(1,jc,jb), tri_idx(2,jc,jb))
          mtgrm_data%station(jc,jb)%soiltype =  &
            &  ext_data%atm%soiltyp(tri_idx(1,jc,jb), tri_idx(2,jc,jb))
        CASE DEFAULT
          mtgrm_data%station(jc,jb)%hsurf    =  0._wp
          mtgrm_data%station(jc,jb)%frland   =  0._wp
          mtgrm_data%station(jc,jb)%soiltype =  0
        END SELECT
        ! initialize value buffer:
        ALLOCATE(mtgrm_data%station(jc,jb)%var(nvars), &
          &      stat=ierrstat)
        IF (ierrstat /= SUCCESS) THEN
          CALL finish (routine, 'ALLOCATE of meteogram data structures failed')
        ENDIF
        DO ivar=1,nvars
          nlevs = mtgrm_data%var_info(ivar)%nlevs
          ALLOCATE(mtgrm_data%station(jc,jb)%var(ivar)%values(nlevs, MAX_TIME_STAMPS), &
            &      stat=ierrstat)
          IF (ierrstat /= SUCCESS) THEN
            CALL finish (routine, 'ALLOCATE of meteogram data structures failed')
          ENDIF
        END DO
        ! initialize value buffer for surface variables:
        ALLOCATE(mtgrm_data%station(jc,jb)%sfc_var(nsfcvars), &
          &      stat=ierrstat)
        IF (ierrstat /= SUCCESS) THEN
          CALL finish (routine, 'ALLOCATE of meteogram data structures failed')
        ENDIF
        DO ivar=1,nsfcvars
          ALLOCATE(mtgrm_data%station(jc,jb)%sfc_var(ivar)%values(MAX_TIME_STAMPS), &
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

    IO_PE : IF (my_process_is_stdio() .AND. .NOT. mtgrm_output_config%ldistributed) THEN

      mtgrm_global_data(jg)%nstations =  mtgrm_output_config%nstations
      mtgrm_global_data(jg)%nblks     =  mtgrm_output_config%nblks    
      mtgrm_global_data(jg)%npromz    =  mtgrm_output_config%npromz   

      ! Note: variable info is not duplicated
      mtgrm_global_data(jg)%nvars     =  mtgrm_local_data(jg)%nvars
      mtgrm_global_data(jg)%nsfcvars  =  mtgrm_local_data(jg)%nsfcvars
      mtgrm_global_data(jg)%max_nlevs =  mtgrm_local_data(jg)%max_nlevs
      mtgrm_global_data(jg)%var_info      =>  mtgrm_local_data(jg)%var_info
      mtgrm_global_data(jg)%sfc_var_info  =>  mtgrm_local_data(jg)%sfc_var_info
      mtgrm_global_data(jg)%time_stamp    =>  mtgrm_local_data(jg)%time_stamp

      nblks_global  = mtgrm_global_data(jg)%nblks
      npromz_global = mtgrm_global_data(jg)%npromz
      ALLOCATE(mtgrm_global_data(jg)%station(nproma, nblks_global), stat=ierrstat)
      IF (ierrstat /= SUCCESS) THEN
        CALL finish (routine, 'ALLOCATE of meteogram data structures failed')
      ENDIF

      DO jb=1,nblks_global
        i_startidx = 1
        i_endidx   = nproma
        IF (jb == nblks_global) i_endidx = npromz_global
        
        DO jc=i_startidx,i_endidx
          p_station => mtgrm_global_data(jg)%station(jc,jb)

          ALLOCATE(p_station%var(nvars), stat=ierrstat)
          IF (ierrstat /= SUCCESS) &
            CALL finish (routine, 'ALLOCATE of meteogram data structures failed')
          
          DO ivar=1,nvars
            nlevs = mtgrm_data%var_info(ivar)%nlevs
            ALLOCATE(p_station%var(ivar)%values(nlevs, MAX_TIME_STAMPS), stat=ierrstat)
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
    
    IF (.NOT. mtgrm_output_config%ldistributed) THEN
      io_rank       = get_mpi_all_workroot_id()
      l_is_io_pe    = (get_my_mpi_all_id() == io_rank)
      l_is_sender   = .NOT. l_is_io_pe
      l_is_receiver = l_is_io_pe

      ! compute maximum buffer size for MPI messages:
      max_var_size    = MAX_TIME_STAMPS*p_real_dp_byte*mtgrm_data%max_nlevs
      max_sfcvar_size = MAX_TIME_STAMPS*p_real_dp_byte
      max_buf_size    = MAX_HEADER_SIZE*p_real_dp_byte          &
        &               + mtgrm_data%nvars*max_var_size         &
        &               + mtgrm_data%nsfcvars*max_sfcvar_size 

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

    CALL mtgrm_open_file(mtgrm_output_config, jg)

  END SUBROUTINE mtgrm_init


  !>
  !! @return .TRUE. if meteogram data will be recorded for this step.
  !!
  !! @par Revision History
  !! Initial implementation  by  F. Prill, DWD (2011-08-22)
  !!
  FUNCTION mtgrm_is_sample_step(mtgrm_output_config, cur_step)
    LOGICAL :: mtgrm_is_sample_step
    ! station data from namelist
    TYPE(t_mtgrm_output_config), TARGET, INTENT(IN) :: mtgrm_output_config
    INTEGER,          INTENT(IN)  :: cur_step     !< current model iteration step

    mtgrm_is_sample_step = &
      &  mtgrm_output_config%lenabled               .AND. &
      &  (cur_step >= mtgrm_output_config%n0_mtgrm) .AND. &
      &  (MOD((cur_step - mtgrm_output_config%n0_mtgrm),  &
      &       mtgrm_output_config%ninc_mtgrm) == 0)
    
  END FUNCTION mtgrm_is_sample_step

  
  !>
  !! Adds values for current model time to buffer.
  !! For gathered NetCDF output, this is a collective operation,
  !! otherwise this is a local operation.
  !!
  !! @par Revision History
  !! Initial implementation  by  F. Prill, DWD (2011-08-22)
  !!
  SUBROUTINE mtgrm_sample_vars(p_nh_state, jg, cur_step, cur_datetime, ierr)
    TYPE(t_nh_state), TARGET, INTENT(IN)  :: p_nh_state
    INTEGER,          INTENT(IN)  :: jg           !< patch index
    INTEGER,          INTENT(IN)  :: cur_step     !< current model iteration step
    TYPE(t_datetime), INTENT(IN)  :: cur_datetime !< date and time of point sample
    INTEGER,          INTENT(OUT) :: ierr         !< error code (e.g. buffer overflow)
    ! local variables
    CHARACTER(*), PARAMETER :: routine = TRIM("mo_mtgrm_output:mtgrm_sample_vars")
    INTEGER :: jb, jc, i_startidx, i_endidx, ilev, &
      &        i_tstep, iidx, iblk, iextra, var_extra
    TYPE(t_mtgrm_data), POINTER :: mtgrm_data
    TYPE(t_nh_prog), POINTER :: prog
    TYPE(t_nh_diag), POINTER :: diag
    TYPE(t_var), POINTER :: VAR

    mtgrm_data => mtgrm_local_data(jg)
    ierr = 0

    VAR => var_list(jg) 
    diag => p_nh_state%diag
    prog => p_nh_state%prog(nnow(jg))

    IF (dbg_level > 0) THEN
      WRITE(message_text,*) "Sampling at step=", cur_step
      CALL message(routine, TRIM(message_text))
    END IF

    ! increase time step counter
    mtgrm_data%icurrent = mtgrm_data%icurrent + 1
    i_tstep = mtgrm_data%icurrent
    IF (i_tstep > MAX_TIME_STAMPS) THEN
      ! buffer full
      ierr = -1
      RETURN
    END IF

    mtgrm_data%time_stamp(i_tstep)%istep = cur_step
    mtgrm_data%time_stamp(i_tstep)%zdate = iso8601(cur_datetime)

    ! fill time step with values
    DO jb=1,mtgrm_data%nblks
      i_startidx = 1
      i_endidx   = nproma
      IF (jb == mtgrm_data%nblks) i_endidx = mtgrm_data%npromz

      DO jc=i_startidx,i_endidx
        iidx  = mtgrm_data%station(jc,jb)%tri_idx_local(1)
        iblk  = mtgrm_data%station(jc,jb)%tri_idx_local(2)

        DO ilev=1,mtgrm_data%var_info(VAR%P)%nlevs
          mtgrm_data%station(jc,jb)%var(VAR%P)%values(ilev, i_tstep) =  &
            &  diag%pres(iidx, mtgrm_data%var_info(VAR%P)%levels(ilev), iblk)
        END DO
        DO ilev=1,mtgrm_data%var_info(VAR%T)%nlevs
          mtgrm_data%station(jc,jb)%var(VAR%T)%values(ilev, i_tstep) =  &
            &  diag%temp(iidx, mtgrm_data%var_info(VAR%T)%levels(ilev), iblk)
        END DO
        DO ilev=1,mtgrm_data%var_info(VAR%RHO)%nlevs
          mtgrm_data%station(jc,jb)%var(VAR%RHO)%values(ilev, i_tstep) =  &
            &  prog%rho(iidx, mtgrm_data%var_info(VAR%RHO)%levels(ilev), iblk)
        END DO
        DO ilev=1,mtgrm_data%var_info(VAR%EXNER)%nlevs
          mtgrm_data%station(jc,jb)%var(VAR%EXNER)%values(ilev, i_tstep) =  &
          &  prog%exner(iidx, mtgrm_data%var_info(VAR%EXNER)%levels(ilev), iblk)
        END DO
        DO ilev=1,mtgrm_data%var_info(VAR%THETAV)%nlevs
          mtgrm_data%station(jc,jb)%var(VAR%THETAV)%values(ilev, i_tstep) =  &
          &  prog%theta_v(iidx, mtgrm_data%var_info(VAR%THETAV)%levels(ilev), iblk)
        END DO
        DO ilev=1,mtgrm_data%var_info(VAR%RHOTHETAV)%nlevs
          mtgrm_data%station(jc,jb)%var(VAR%RHOTHETAV)%values(ilev, i_tstep) =  &
          &  prog%rhotheta_v(iidx, mtgrm_data%var_info(VAR%RHOTHETAV)%levels(ilev), iblk)
        END DO
        DO ilev=1,mtgrm_data%var_info(VAR%U)%nlevs
          mtgrm_data%station(jc,jb)%var(VAR%U)%values(ilev, i_tstep) =  &
          &  diag%u(iidx, mtgrm_data%var_info(VAR%U)%levels(ilev), iblk)
        END DO
        DO ilev=1,mtgrm_data%var_info(VAR%V)%nlevs
          mtgrm_data%station(jc,jb)%var(VAR%V)%values(ilev, i_tstep) =  &
          &  diag%v(iidx, mtgrm_data%var_info(VAR%V)%levels(ilev), iblk)
        END DO
        DO ilev=1,mtgrm_data%var_info(VAR%W)%nlevs
          mtgrm_data%station(jc,jb)%var(VAR%W)%values(ilev, i_tstep) =  &
          &  prog%w(iidx, mtgrm_data%var_info(VAR%W)%levels(ilev), iblk)
        END DO

        mtgrm_data%station(jc,jb)%sfc_var(VAR%P_SFC)%values(i_tstep) =  &
          &  diag%pres_sfc(iidx, iblk)

        IF (lwrite_extra) THEN
          ! Variable: Extra 2D
          DO iextra=1,inextra_2d
            var_extra = VAR%EXTRA2D(iextra)
            mtgrm_data%station(jc,jb)%sfc_var(var_extra)%values(i_tstep) =  &
              &  diag%extra_2d(iidx, iblk, iextra)
          END DO
          ! Variable: Extra 3D
          DO iextra=1,inextra_3d
            var_extra = VAR%EXTRA3D(iextra)
            DO ilev=1,mtgrm_data%var_info(var_extra)%nlevs
              mtgrm_data%station(jc,jb)%var(var_extra)%values(ilev, i_tstep) =  &
                &  diag%extra_3d(iidx, mtgrm_data%var_info(var_extra)%levels(ilev), iblk, iextra)
            END DO
          END DO
        END IF
      END DO
    END DO

  END SUBROUTINE mtgrm_sample_vars


  !>
  !! Destroy meteogram data structure.
  !! For gathered NetCDF output, this is a collective operation,
  !! otherwise this is a local operation.
  !!
  !! @par Revision History
  !! Initial implementation  by  F. Prill, DWD (2011-08-22)
  !!
  SUBROUTINE mtgrm_finalize(mtgrm_output_config, jg)
    ! station data from namelist
    TYPE(t_mtgrm_output_config), TARGET, INTENT(IN) :: mtgrm_output_config
    INTEGER, INTENT(IN)  :: jg    !< patch index
    ! local variables:
    CHARACTER(*), PARAMETER     :: routine = TRIM("mo_mtgrm_output:mtgrm_finalize")
    INTEGER                     :: ierrstat, jb, jc, i_startidx, i_endidx, &
      &                            nvars, nsfcvars, ivar
    TYPE(t_mtgrm_data), POINTER :: mtgrm_data

    mtgrm_data => mtgrm_local_data(jg)

    ! ------------------------------------------------------------
    ! If this is the IO PE: close NetCDF file
    ! ------------------------------------------------------------

    CALL mtgrm_close_file(mtgrm_output_config, jg)

    DEALLOCATE(mtgrm_data%time_stamp, stat=ierrstat)
    IF (ierrstat /= SUCCESS) THEN
      CALL finish (routine, 'DEALLOCATE of meteogram data structures failed')
    ENDIF

    nvars    = mtgrm_data%nvars
    nsfcvars = mtgrm_data%nsfcvars
    DO ivar=1,nvars
      IF (ASSOCIATED(mtgrm_data%var_info(ivar)%levels)) THEN
        DEALLOCATE(mtgrm_data%var_info(ivar)%levels, stat=ierrstat)
        IF (ierrstat /= SUCCESS) THEN
          CALL finish (routine, 'DEALLOCATE of meteogram data structures failed')
        ENDIF
      END IF
    END DO
    DEALLOCATE(mtgrm_data%var_info, stat=ierrstat)
    IF (ierrstat /= SUCCESS) THEN
      CALL finish (routine, 'DEALLOCATE of meteogram data structures failed')
    ENDIF
    DEALLOCATE(mtgrm_data%sfc_var_info, stat=ierrstat)
    IF (ierrstat /= SUCCESS) THEN
      CALL finish (routine, 'DEALLOCATE of meteogram data structures failed')
    ENDIF
    
    DO jb=1,mtgrm_data%nblks
      i_startidx = 1
      i_endidx   = nproma
      IF (jb == mtgrm_data%nblks) &
        &  i_endidx = mtgrm_data%npromz

      DO jc=i_startidx,i_endidx
        DO ivar=1,nvars
          DEALLOCATE(mtgrm_data%station(jc,jb)%var(ivar)%values, &
            &        stat=ierrstat)
          IF (ierrstat /= SUCCESS) THEN
            CALL finish (routine, 'DEALLOCATE of meteogram data structures failed')
          ENDIF
        END DO
        DEALLOCATE(mtgrm_data%station(jc,jb)%var, stat=ierrstat)
        IF (ierrstat /= SUCCESS) THEN
          CALL finish (routine, 'DEALLOCATE of meteogram data structures failed')
        ENDIF
        DO ivar=1,nsfcvars
          DEALLOCATE(mtgrm_data%station(jc,jb)%sfc_var(ivar)%values, &
            &        stat=ierrstat)
          IF (ierrstat /= SUCCESS) THEN
            CALL finish (routine, 'DEALLOCATE of meteogram data structures failed')
          ENDIF
        END DO
        DEALLOCATE(mtgrm_data%station(jc,jb)%sfc_var, stat=ierrstat)
        IF (ierrstat /= SUCCESS) THEN
          CALL finish (routine, 'DEALLOCATE of meteogram data structures failed')
        ENDIF

      END DO
    END DO
    DEALLOCATE(mtgrm_data%station, stat=ierrstat)
    IF (ierrstat /= SUCCESS) THEN
      CALL finish (routine, 'DEALLOCATE of meteogram data structures failed')
    ENDIF
    mtgrm_data%nvars     = 0
    mtgrm_data%nsfcvars  = 0
    mtgrm_data%nstations = 0

    ! deallocate MPI buffer
    IF (.NOT. mtgrm_output_config%ldistributed) THEN
      DEALLOCATE(msg_buffer, stat=ierrstat)
      IF (ierrstat /= SUCCESS) THEN
        CALL finish (routine, 'DEALLOCATE of MPI message buffer failed')
      ENDIF
    END IF

    ! deallocate global meteogram data

    IO_PE : IF (my_process_is_stdio() .AND. .NOT. mtgrm_output_config%ldistributed) THEN
    
      DO jb=1,mtgrm_global_data(jg)%nblks
        i_startidx = 1
        i_endidx   = nproma
        IF (jb == mtgrm_global_data(jg)%nblks) i_endidx = mtgrm_global_data(jg)%npromz
        
        DO jc=i_startidx,i_endidx
          DO ivar=1,nvars
            DEALLOCATE(mtgrm_global_data(jg)%station(jc,jb)%var(ivar)%values, stat=ierrstat)
            IF (ierrstat /= SUCCESS) &
              CALL finish (routine, 'ALLOCATE of meteogram data structures failed')
          END DO
          DEALLOCATE(mtgrm_global_data(jg)%station(jc,jb)%var, stat=ierrstat)
          IF (ierrstat /= SUCCESS) &
            CALL finish (routine, 'ALLOCATE of meteogram data structures failed')

          DO ivar=1,nsfcvars
            DEALLOCATE(mtgrm_global_data(jg)%station(jc,jb)%sfc_var(ivar)%values, stat=ierrstat)
            IF (ierrstat /= SUCCESS) &
              CALL finish (routine, 'ALLOCATE of meteogram data structures failed')
          END DO
          DEALLOCATE(mtgrm_global_data(jg)%station(jc,jb)%sfc_var, stat=ierrstat)
          IF (ierrstat /= SUCCESS) &
            CALL finish (routine, 'ALLOCATE of meteogram data structures failed')

        END DO
      END DO
      DEALLOCATE(mtgrm_global_data(jg)%station, stat=ierrstat)
      IF (ierrstat /= SUCCESS) &
        CALL finish (routine, 'ALLOCATE of meteogram data structures failed')
      
    END IF IO_PE

  END SUBROUTINE mtgrm_finalize


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
  SUBROUTINE mtgrm_collect_buffers(mtgrm_output_config, jg, io_rank)
    ! station data from namelist
    TYPE(t_mtgrm_output_config), TARGET, INTENT(IN) :: mtgrm_output_config
    INTEGER, INTENT(IN)  :: jg       !< patch index
    INTEGER, INTENT(IN)  :: io_rank  !< MPI rank of root process

#ifndef NOMPI
    ! local variables
    CHARACTER(*), PARAMETER :: routine = TRIM("mo_mtgrm_output:mtgrm_collect_buffers")
    INTEGER     :: station_idx(2), position, icurrent,   &
      &            jb, jc, i_startidx, i_endidx, ierr,   &
      &            istation, ivar, nlevs
    INTEGER     :: irecv_req(MAX_NUM_STATIONS)
    INTEGER     :: irecv_status(MPI_STATUS_SIZE, MAX_NUM_STATIONS)
    INTEGER     :: myrank, comm
    TYPE(t_mtgrm_data),    POINTER :: mtgrm_data
    TYPE(t_mtgrm_station), POINTER :: p_station

    mtgrm_data => mtgrm_local_data(jg)
    ! skip routine, if there is nothing to do...
    IF (mtgrm_global_data(jg)%nstations == 0) RETURN

    ! get global rank and MPI communicator
    myrank = get_my_mpi_all_id()
    comm   = get_my_mpi_communicator()

    ! global time stamp index
    ! Note: We assume that this value is identical for all PEs
    icurrent = mtgrm_data%icurrent

    ! -- RECEIVER CODE --
    RECEIVER : IF (l_is_receiver) THEN
      ! launch MPI message requests for station data on foreign PEs
      irecv_req(:) = MPI_REQUEST_NULL
      DO istation=1,mtgrm_output_config%nstations
        IF (owner(istation) /= myrank) THEN
          CALL MPI_IRECV(msg_buffer(:,istation), max_buf_size, MPI_PACKED, MPI_ANY_SOURCE, &
            &            istation, comm, irecv_req(istation), ierr)
          IF (ierr /= MPI_SUCCESS) CALL finish (routine, 'MPI function call failed')
        END IF
      END DO

      ! wait for messages to arrive:
      CALL MPI_WAITALL(mtgrm_output_config%nstations, irecv_req(:), irecv_status(:,:), ierr)
      IF (ierr /= MPI_SUCCESS) CALL finish (routine, 'MPI function call failed')

      ! unpack received messages:
      jc = 0
      jb = 1
      DO istation=1,mtgrm_output_config%nstations
        IF (owner(istation) /= myrank) THEN
          position = 0
          CALL MPI_UNPACK(msg_buffer(:,istation), max_buf_size, position, station_idx(:), &
            &             2, p_int, comm, ierr)
          IF (ierr /= MPI_SUCCESS) CALL finish (routine, 'MPI function call failed')
          p_station => mtgrm_global_data(jg)%station(station_idx(1), station_idx(2))
          p_station%station_idx(1:2) = station_idx(1:2)

          ! unpack header information
          CALL MPI_UNPACK(msg_buffer(:,istation), max_buf_size, position, p_station%tri_idx(:), &
            &             2, p_int, comm, ierr)
          IF (ierr /= MPI_SUCCESS) CALL finish (routine, 'MPI function call failed')
          CALL MPI_UNPACK(msg_buffer(:,istation), max_buf_size, position, &
            &             p_station%tri_idx_local(:), &
            &             2, p_int, comm, ierr)
          IF (ierr /= MPI_SUCCESS) CALL finish (routine, 'MPI function call failed')
          CALL MPI_UNPACK(msg_buffer(:,istation), max_buf_size, position, p_station%owner, &
            &             1, p_int, comm, ierr)
          IF (ierr /= MPI_SUCCESS) CALL finish (routine, 'MPI function call failed')
          CALL MPI_UNPACK(msg_buffer(:,istation), max_buf_size, position, p_station%hsurf,      &
            &             1, p_real_dp, comm, ierr)
          IF (ierr /= MPI_SUCCESS) CALL finish (routine, 'MPI function call failed')
          CALL MPI_UNPACK(msg_buffer(:,istation), max_buf_size, position, p_station%frland,     &
            &             1, p_real_dp, comm, ierr)
          IF (ierr /= MPI_SUCCESS) CALL finish (routine, 'MPI function call failed')
          CALL MPI_UNPACK(msg_buffer(:,istation), max_buf_size, position, p_station%fc,         &
            &             1, p_real_dp, comm, ierr)
          IF (ierr /= MPI_SUCCESS) CALL finish (routine, 'MPI function call failed')
          CALL MPI_UNPACK(msg_buffer(:,istation), max_buf_size, position, p_station%soiltype,   &
            &             1, p_int, comm, ierr)
          IF (ierr /= MPI_SUCCESS) CALL finish (routine, 'MPI function call failed')

          ! unpack meteogram data:
          DO ivar=1,mtgrm_data%nvars
            nlevs = mtgrm_data%var_info(ivar)%nlevs
            CALL MPI_UNPACK(msg_buffer(:,istation), max_buf_size, position, &
              &             p_station%var(ivar)%values(:,:),                &
              &             nlevs*icurrent, p_real_dp, comm, ierr)
            IF (ierr /= MPI_SUCCESS) CALL finish (routine, 'MPI function call failed')
          END DO

          DO ivar=1,mtgrm_data%nsfcvars
            CALL MPI_UNPACK(msg_buffer(:,istation), max_buf_size, position, &
              &             p_station%sfc_var(ivar)%values(:),              &
              &             icurrent, p_real_dp, comm, ierr)
            IF (ierr /= MPI_SUCCESS) CALL finish (routine, 'MPI function call failed')
          END DO
        ELSE
          ! this PE is both sender and receiver - direct copy:
          jc = jc + 1
          IF (jc > nproma) THEN
            jc = 1
            jb = jb + 1
          END IF
          station_idx(1:2) = mtgrm_data%station(jc,jb)%station_idx(1:2)

          p_station => mtgrm_global_data(jg)%station(station_idx(1),station_idx(2))
          p_station%station_idx(1:2)   = mtgrm_data%station(jc,jb)%station_idx(1:2)
          p_station%tri_idx(1:2)       = mtgrm_data%station(jc,jb)%tri_idx(1:2)
          p_station%tri_idx_local(1:2) = mtgrm_data%station(jc,jb)%tri_idx_local(1:2)
          p_station%owner              = mtgrm_data%station(jc,jb)%owner
          p_station%hsurf              = mtgrm_data%station(jc,jb)%hsurf
          p_station%frland             = mtgrm_data%station(jc,jb)%frland
          p_station%fc                 = mtgrm_data%station(jc,jb)%fc
          p_station%soiltype           = mtgrm_data%station(jc,jb)%soiltype
          ! copy meteogram data
          DO ivar=1,mtgrm_data%nvars
            nlevs = mtgrm_data%var_info(ivar)%nlevs
            p_station%var(ivar)%values(1:nlevs, 1:icurrent) =  &
              &  mtgrm_data%station(jc,jb)%var(ivar)%values(1:nlevs, 1:icurrent)
          END DO
          DO ivar=1,mtgrm_data%nsfcvars
            p_station%sfc_var(ivar)%values(1:icurrent) =  &
              &  mtgrm_data%station(jc,jb)%sfc_var(ivar)%values(1:icurrent)
          END DO
        END IF
      END DO
    END IF RECEIVER

    ! -- SENDER CODE --
    SENDER : IF (l_is_sender) THEN
      ! pack station into buffer; send it
      DO jb=1,mtgrm_data%nblks
        i_startidx = 1
        i_endidx   = nproma
        IF (jb == mtgrm_data%nblks) i_endidx = mtgrm_data%npromz
        DO jc=i_startidx,i_endidx

          ! Meteogram header (information on location, ...)
          p_station => mtgrm_data%station(jc,jb)
          position = 0
          CALL MPI_PACK(p_station%station_idx(:), 2, p_int,     msg_buffer(:,1), max_buf_size, &
            &           position, comm, ierr)
          IF (ierr /= MPI_SUCCESS) CALL finish (routine, 'MPI function call failed')
          CALL MPI_PACK(p_station%tri_idx(:),     2, p_int,     msg_buffer(:,1), max_buf_size, &
            &           position, comm, ierr)
          IF (ierr /= MPI_SUCCESS) CALL finish (routine, 'MPI function call failed')
          CALL MPI_PACK(p_station%tri_idx_local(:),2,p_int,     msg_buffer(:,1), max_buf_size, &
            &           position, comm, ierr)
          IF (ierr /= MPI_SUCCESS) CALL finish (routine, 'MPI function call failed')
          CALL MPI_PACK(p_station%owner,          1 ,p_int,     msg_buffer(:,1), max_buf_size, &
            &           position, comm, ierr)
          IF (ierr /= MPI_SUCCESS) CALL finish (routine, 'MPI function call failed')
          CALL MPI_PACK(p_station%hsurf,          1, p_real_dp, msg_buffer(:,1), max_buf_size, &
            &           position, comm, ierr)
          IF (ierr /= MPI_SUCCESS) CALL finish (routine, 'MPI function call failed')
          CALL MPI_PACK(p_station%frland,         1, p_real_dp, msg_buffer(:,1), max_buf_size, &
            &           position, comm, ierr)
          IF (ierr /= MPI_SUCCESS) CALL finish (routine, 'MPI function call failed')
          CALL MPI_PACK(p_station%fc,             1, p_real_dp, msg_buffer(:,1), max_buf_size, &
            &           position, comm, ierr)
          IF (ierr /= MPI_SUCCESS) CALL finish (routine, 'MPI function call failed')
          CALL MPI_PACK(p_station%soiltype,       1, p_int,     msg_buffer(:,1), max_buf_size, &
            &           position, comm, ierr)
          IF (ierr /= MPI_SUCCESS) CALL finish (routine, 'MPI function call failed')
          ! pack meteogram data:
          DO ivar=1,mtgrm_data%nvars
            nlevs = mtgrm_data%var_info(ivar)%nlevs
            CALL MPI_PACK(p_station%var(ivar)%values(:,:), nlevs*icurrent, &
              &           p_real_dp, msg_buffer(:,1),                      &
              &           max_buf_size, position, comm, ierr)
            IF (ierr /= MPI_SUCCESS) CALL finish (routine, 'MPI function call failed')
          END DO
          DO ivar=1,mtgrm_data%nsfcvars
            CALL MPI_PACK(p_station%sfc_var(ivar)%values(:), icurrent,     &
              &           p_real_dp, msg_buffer(:,1),                      &
              &           max_buf_size, position, comm, ierr)
            IF (ierr /= MPI_SUCCESS) CALL finish (routine, 'MPI function call failed')
          END DO

          ! (blocking) send of packed station data to IO PE:
          istation = nproma*(p_station%station_idx(2) - 1) + p_station%station_idx(1)
          CALL MPI_SEND(msg_buffer, position, MPI_PACKED, io_rank, istation, comm, ierr);
          IF (ierr /= MPI_SUCCESS) CALL finish (routine, 'MPI function call failed')
        END DO
      END DO

      ! reset buffer on sender side
      IF (.NOT. l_is_receiver) &
        mtgrm_data%icurrent = 0

    END IF SENDER

#endif
  END SUBROUTINE mtgrm_collect_buffers


  !>
  !! The IO PE creates and opens a disk file for output.
  !! For gathered NetCDF output, this is a collective operation,
  !! otherwise this is a local operation for the IO PE.
  !!
  !! @par Revision History
  !! Initial implementation  by  F. Prill, DWD (2011-08-22)
  !!
  SUBROUTINE mtgrm_open_file(mtgrm_output_config, jg)
    ! station data from namelist
    TYPE(t_mtgrm_output_config), TARGET, INTENT(IN) :: mtgrm_output_config
    ! patch index
    INTEGER,                             INTENT(IN) :: jg
    ! local variables:
    CHARACTER(*), PARAMETER :: routine = TRIM("mo_mtgrm_output:mtgrm_open_file")
    INTEGER                     :: jb, jc, i_startidx, i_endidx, old_mode, ncfile, &
      &                            istation, ivar, nvars, nsfcvars
    TYPE(t_ncid),       POINTER :: ncid
    INTEGER                     :: station_name_dims(2), var_name_dims(2), &
      &                            var_level_dims(2), time_string_dims(2), &
      &                            var_dims(4),  sfcvar_dims(3),           &
      &                            istart2(2), icount2(2)
    CHARACTER(len=MAX_NAME_LENGTH) :: description_str
    TYPE(t_station_list), POINTER  :: this_station
    TYPE(t_mtgrm_data), POINTER    :: mtgrm_data

    IF (mtgrm_output_config%ftype /= FTYPE_NETCDF) THEN
      CALL finish(routine, "Output format not yet implemented.")
    END IF

    ! In "non-distributed" mode, station data is gathered by PE #0
    ! which writes a single file.
    ! Note that info on variables is not copied to the global data set
    ! (we use the local mtgrm_data there).

    IF (.NOT. mtgrm_output_config%ldistributed) THEN
      CALL mtgrm_collect_buffers(mtgrm_output_config, jg, io_rank)
      mtgrm_data => mtgrm_global_data(jg)
    ELSE
      mtgrm_data => mtgrm_local_data(jg)
      ! skip routine, if this PE has nothing to do...
      IF (mtgrm_local_data(jg)%nstations == 0) RETURN
    END IF

    ! skip routine, if this PE has nothing to do...
    IF  (.NOT. mtgrm_output_config%ldistributed .AND.  &
      & (.NOT. l_is_receiver)) RETURN

    ncid => ncid_list(jg)
    nvars    = mtgrm_data%nvars
    nsfcvars = mtgrm_data%nsfcvars

    ! create a file name for this PE:
    CALL mtgrm_create_filename(mtgrm_output_config, jg)

    ! create NetCDF file:
    CALL nf(nf_set_default_format(nf_format_64bit, old_mode))
    CALL nf(nf_create(TRIM(mtgrm_file_info(jg)%zname), nf_clobber, &
      &               mtgrm_file_info(jg)%file_id))
    ncfile = mtgrm_file_info(jg)%file_id
    CALL nf(nf_set_fill(ncfile, nf_nofill, old_mode))

    description_str = "ICON Meteogram File"
    CALL nf(nf_put_att_text(ncfile, NF_GLOBAL, "info", &
      &         LEN(TRIM(description_str)), TRIM(description_str)))

    ! for the definition of a character-string variable define
    ! character-position dimension for strings of max length 40
    CALL nf(nf_def_dim(ncfile, "charid", MAX_NAME_LENGTH, ncid%charid))
    ! station header:
    CALL nf(nf_def_dim(ncfile, 'nstations',  mtgrm_data%nstations,   ncid%nstations))
    ! write variables:
    CALL nf(nf_def_dim(ncfile, 'nvars',      mtgrm_data%nvars, ncid%nvars))
    IF (mtgrm_data%nsfcvars > 0) &
      CALL nf(nf_def_dim(ncfile, 'nsfcvars', mtgrm_data%nsfcvars,  ncid%nsfcvars))
    CALL nf(nf_def_dim(ncfile, 'max_nlevs',  mtgrm_data%max_nlevs, ncid%max_nlevs))
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
    IF (mtgrm_data%nsfcvars > 0) THEN
      var_name_dims = (/ ncid%charid, ncid%nsfcvars /)
      CALL nf(nf_def_var(ncfile, "sfcvar_name", NF_CHAR, 2, var_name_dims(:), &
        &                ncid%sfcvar_name))
      CALL nf_add_descr("Surface variable name (character string)", ncfile, ncid%sfcvar_name)
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

    ! add value buffer for volume variables:
    var_dims = (/ ncid%nstations, ncid%nvars, ncid%max_nlevs, ncid%timeid /)
    CALL nf(nf_def_var(ncfile, "values", NF_DOUBLE, 4, var_dims(:), &
      &                ncid%var_values))
    CALL nf_add_descr("value buffer for volume variables", ncfile, ncid%var_values)
    ! add value buffer for surface variables:
    IF (mtgrm_data%nsfcvars > 0) THEN
      sfcvar_dims = (/ ncid%nstations, ncid%nvars, ncid%timeid /)
      CALL nf(nf_def_var(ncfile, "sfcvalues", NF_DOUBLE, 3, sfcvar_dims(:), &
        &                ncid%sfcvar_values))
      CALL nf_add_descr("value buffer for surface variables", ncfile, ncid%sfcvar_values)
    END IF

    ! ----------------------
    ! End of definition mode
    CALL nf(nf_enddef(ncfile))

    DO ivar=1,nvars
      CALL nf(nf_put_vara_text(ncfile, ncid%var_name, (/ 1, ivar /), &
        &                      (/ LEN(TRIM(mtgrm_data%var_info(ivar)%zname)), 1 /), &
        &                      TRIM(mtgrm_data%var_info(ivar)%zname)))
      CALL nf(nf_put_vara_text(ncfile, ncid%var_unit, (/ 1, ivar /), &
        &                      (/ LEN(TRIM(mtgrm_data%var_info(ivar)%zunit)), 1 /), &
        &                      TRIM(mtgrm_data%var_info(ivar)%zunit)))
      CALL nf(nf_put_vara_int(ncfile, ncid%var_group_id, ivar, 1, &
        &                     mtgrm_data%var_info(ivar)%igroup_id))
      CALL nf(nf_put_vara_int(ncfile, ncid%var_nlevs, ivar, 1, &
        &                     mtgrm_data%var_info(ivar)%nlevs))
      istart2 = (/ 1, ivar /)
      icount2 = (/ mtgrm_data%var_info(ivar)%nlevs, 1 /)
      CALL nf(nf_put_vara_int(ncfile, ncid%var_levels, istart2, icount2, &
        &                     mtgrm_data%var_info(ivar)%levels(:)))
    END DO

    DO ivar=1,nsfcvars
      CALL nf(nf_put_vara_text(ncfile, ncid%sfcvar_name, (/ 1, ivar /), &
        &                      (/ LEN(TRIM(mtgrm_data%sfc_var_info(ivar)%zname)), 1 /), &
        &                      TRIM(mtgrm_data%sfc_var_info(ivar)%zname)))
      CALL nf(nf_put_vara_text(ncfile, ncid%sfcvar_unit, (/ 1, ivar /), &
        &                      (/ LEN(TRIM(mtgrm_data%sfc_var_info(ivar)%zunit)), 1 /), &
        &                      TRIM(mtgrm_data%sfc_var_info(ivar)%zunit)))
      CALL nf(nf_put_vara_int(ncfile, ncid%sfcvar_group_id, ivar, 1, &
        &                     mtgrm_data%sfc_var_info(ivar)%igroup_id))
    END DO

    istation = 1
    DO jb=1,mtgrm_data%nblks
      i_startidx = 1
      i_endidx   = nproma
      IF (jb == mtgrm_data%nblks) i_endidx = mtgrm_data%npromz

      DO jc=i_startidx,i_endidx
        this_station => mtgrm_output_config%station_list(           &
          &               mtgrm_data%station(jc,jb)%station_idx(1), &
          &               mtgrm_data%station(jc,jb)%station_idx(2))
        CALL nf(nf_put_vara_text(ncfile, ncid%station_name, (/ 1, istation /), &
          &                      (/ LEN(TRIM(this_station%zname)), 1 /), &
          &                      TRIM(this_station%zname)))
        CALL nf(nf_put_vara_double(ncfile, ncid%station_lon, istation, 1, &
          &                        this_station%location%lon))
        CALL nf(nf_put_vara_double(ncfile, ncid%station_lat, istation, 1, &
          &                        this_station%location%lat))
        CALL nf(nf_put_vara_int(ncfile, ncid%station_idx, istation, 1, &
          &                     mtgrm_data%station(jc,jb)%tri_idx(1)))
        CALL nf(nf_put_vara_int(ncfile, ncid%station_blk, istation, 1, &
          &                     mtgrm_data%station(jc,jb)%tri_idx(2)))
        CALL nf(nf_put_vara_double(ncfile, ncid%station_hsurf, istation, 1, &
          &                        mtgrm_data%station(jc,jb)%hsurf))
        CALL nf(nf_put_vara_double(ncfile, ncid%station_frland, istation, 1, &
          &                        mtgrm_data%station(jc,jb)%frland))
        CALL nf(nf_put_vara_double(ncfile, ncid%station_fc, istation, 1, &
          &                        mtgrm_data%station(jc,jb)%fc))
        CALL nf(nf_put_vara_int(ncfile, ncid%station_soiltype, istation, 1, &
          &                     mtgrm_data%station(jc,jb)%soiltype))

        istation = istation + 1 

      END DO
    END DO

  END SUBROUTINE mtgrm_open_file


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
  SUBROUTINE mtgrm_flush_file(mtgrm_output_config, jg)
    ! station data from namelist
    TYPE(t_mtgrm_output_config), TARGET, INTENT(IN) :: mtgrm_output_config
    INTEGER, INTENT(IN)         :: jg       !< patch index
    ! local variables:
    CHARACTER(*), PARAMETER :: routine = TRIM("mo_mtgrm_output:mtgrm_flush_file")
    INTEGER                     :: ncfile,  totaltime, itime, istation, ivar, &
      &                            jb, jc, i_startidx, i_endidx, nlevs,       &
      &                            nvars, nsfcvars
    TYPE(t_mtgrm_data), POINTER :: mtgrm_data
    TYPE(t_ncid),       POINTER :: ncid
    INTEGER                     :: istart4(4), icount4(4)
  
    ncid => ncid_list(jg)
    ncfile = mtgrm_file_info(jg)%file_id

    ! In "non-distributed" mode, station data is gathered by PE #0
    ! which writes a single file:
    IF (mtgrm_output_config%ldistributed) THEN
      mtgrm_data => mtgrm_local_data(jg)
    ELSE
      CALL mtgrm_collect_buffers(mtgrm_output_config, jg, io_rank)
      mtgrm_data => mtgrm_global_data(jg)
    END IF

    ! skip routine, if this PE has nothing to do...
    IF (mtgrm_data%nstations == 0) RETURN
    IF  (.NOT. mtgrm_output_config%ldistributed .AND.  &
      & (.NOT. l_is_receiver)) RETURN

    nvars    = mtgrm_data%nvars
    nsfcvars = mtgrm_data%nsfcvars

    IF (dbg_level > 0) THEN
      WRITE(message_text,*) "Meteogram"
      CALL message(routine, TRIM(message_text))    
    END IF

    ! inquire about current number of records in file:
    CALL nf(nf_inq_dimlen(ncfile, ncid%timeid, totaltime))

    ! write time stamp info:
    DO itime=1,mtgrm_local_data(jg)%icurrent

      CALL nf(nf_put_vara_text(ncfile, ncid%dateid, (/ 1, totaltime+itime /), &
        &                      (/ LEN(TRIM(mtgrm_data%time_stamp(itime)%zdate)), 1 /), &
        &                      TRIM(mtgrm_data%time_stamp(itime)%zdate)))
      CALL nf(nf_put_vara_int(ncfile, ncid%time_step, totaltime+itime, 1, &
        &                     mtgrm_data%time_stamp(itime)%istep))

      ! write meteogram buffer:
      istation = 1
      DO jb=1,mtgrm_data%nblks
        i_startidx = 1
        i_endidx   = nproma
        IF (jb == mtgrm_data%nblks)  &
          &  i_endidx = mtgrm_data%npromz

        DO jc=i_startidx,i_endidx

          ! volume variables:
          DO ivar=1,nvars
            nlevs = mtgrm_data%var_info(ivar)%nlevs
            istart4 = (/ istation, ivar, 1, totaltime+itime /)
            icount4 = (/ 1, 1, nlevs, 1 /)
            CALL nf(nf_put_vara_double(ncfile, ncid%var_values,                   &
              &                        istart4, icount4,                          &
              &                        mtgrm_data%station(jc,jb)%var(ivar)%values(1:nlevs, itime)))
          END DO
          ! surface variables:
          DO ivar=1,nsfcvars
            CALL nf(nf_put_vara_double(ncfile, ncid%sfcvar_values,                &
              &                        (/ istation, ivar, totaltime+itime /),     &
              &                        (/ 1, 1, 1 /),                             &
              &                        mtgrm_data%station(jc,jb)%sfc_var(ivar)%values(itime)))
          END DO

          istation = istation + 1
        END DO
      END DO
    END DO

    ! finally, reset buffer counter for new data
    mtgrm_local_data(jg)%icurrent = 0

  END SUBROUTINE mtgrm_flush_file


  !>
  !! The IO PE closes the meteogram output file.
  !! For gathered NetCDF output, this is a collective operation,
  !! otherwise this is a local operation for the IO PE.
  !!
  !! @par Revision History
  !! Initial implementation  by  F. Prill, DWD (2011-08-22)
  !!
  SUBROUTINE mtgrm_close_file(mtgrm_output_config, jg)
    ! station data from namelist
    TYPE(t_mtgrm_output_config), TARGET, INTENT(IN) :: mtgrm_output_config
    INTEGER, INTENT(IN)  :: jg    !< patch index

    ! write remaining buffers:
    CALL mtgrm_flush_file(mtgrm_output_config, jg)

    ! skip routine, if this PE has nothing to do...
    IF  (mtgrm_output_config%ldistributed .OR.  &
      & (io_rank == get_my_mpi_all_id())) THEN
      ! Close NetCDF file
      CALL nf(nf_close(mtgrm_file_info(jg)%file_id))
    END IF
  END SUBROUTINE mtgrm_close_file


  !>
  !! @return file name of meteogram file.
  !! This is a local operation.
  !!
  !! @par Revision History
  !! Initial implementation  by  F. Prill, DWD (2011-08-22)
  !!
  SUBROUTINE mtgrm_create_filename (mtgrm_output_config, jg)
    
    ! station data from namelist
    TYPE(t_mtgrm_output_config), TARGET, INTENT(IN) :: mtgrm_output_config
    ! patch index
    INTEGER,                             INTENT(IN) :: jg
    ! Local variables
    INTEGER :: my_id

    my_id = get_my_mpi_all_id()

    SELECT CASE (mtgrm_output_config%ftype)
    CASE (FTYPE_NETCDF)
      IF (mtgrm_output_config%ldistributed) THEN
        WRITE (mtgrm_file_info(jg)%zname,'(a,i3.3,a,i3.3,a)') "PE", my_id, "_patch", jg, ".nc"
      ELSE
        WRITE (mtgrm_file_info(jg)%zname,'(a,i3.3,a)') "patch", jg, ".nc"
      END IF
    END SELECT
    mtgrm_file_info(jg)%zname = TRIM(mtgrm_output_config%zprefix)//TRIM(mtgrm_file_info(jg)%zname)
  END SUBROUTINE mtgrm_create_filename


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
  !!  Utility function.
  !!
  !!  Registers a new atmospheric (volume) variable.
  SUBROUTINE add_atmo_var(ivar, igroup_id, zname, zunit, jg, nlev)
    CHARACTER(LEN=*), INTENT(IN)    :: zname, zunit
    INTEGER,          INTENT(INOUT) :: ivar
    INTEGER,          INTENT(IN)    :: igroup_id, nlev, jg
    ! Local variables
    CHARACTER(*), PARAMETER :: routine = TRIM("mo_mtgrm_output:add_atmo_var")
    TYPE(t_mtgrm_data), POINTER     :: mtgrm_data
    INTEGER                         :: ierrstat, ilev

    mtgrm_data => mtgrm_local_data(jg)

    ! create new variable index
    var_list(jg)%no_atmo_vars = var_list(jg)%no_atmo_vars + 1
    ivar = var_list(jg)%no_atmo_vars
    IF (ivar > mtgrm_data%nvars) THEN
      CALL finish(routine, 'Number of sampling variables exceeds preset value "nvars"')
    END IF

    ! create meteogram data structure
    mtgrm_data%var_info(ivar)%zname     = TRIM(zname)
    mtgrm_data%var_info(ivar)%zunit     = TRIM(zunit)
    mtgrm_data%var_info(ivar)%igroup_id = igroup_id
    mtgrm_data%var_info(ivar)%nlevs     = nlev
    mtgrm_data%max_nlevs = MAX(mtgrm_data%max_nlevs,  mtgrm_data%var_info(ivar)%nlevs)
    ALLOCATE(mtgrm_data%var_info(ivar)%levels(nlev), stat=ierrstat)
    IF (ierrstat /= SUCCESS) THEN
      CALL finish (routine, 'ALLOCATE of meteogram data structures failed')
    ENDIF
    mtgrm_data%var_info(ivar)%levels = (/ (ilev, ilev=1,nlev) /)
  END SUBROUTINE add_atmo_var


  !>
  !!  Utility function.
  !!
  !!  Registers a new surface variable.
  SUBROUTINE add_sfc_var(ivar, igroup_id, zname, zunit, jg)
    CHARACTER(LEN=*), INTENT(IN)    :: zname, zunit
    INTEGER,          INTENT(INOUT) :: ivar
    INTEGER,          INTENT(IN)    :: igroup_id, jg
    ! Local variables
    CHARACTER(*), PARAMETER :: routine = TRIM("mo_mtgrm_output:add_sfc_var")
    TYPE(t_mtgrm_data), POINTER     :: mtgrm_data
 
    mtgrm_data => mtgrm_local_data(jg)

    ! create new variable index
    var_list(jg)%no_sfc_vars = var_list(jg)%no_sfc_vars + 1
    ivar = var_list(jg)%no_sfc_vars
    IF (ivar > mtgrm_data%nsfcvars) THEN
      CALL finish(routine, 'Number of sampling variables exceeds preset value "nsfcvars"')
    END IF
    ! create meteogram data structure
    mtgrm_data%sfc_var_info(ivar)%zname     = TRIM(zname)
    mtgrm_data%sfc_var_info(ivar)%zunit     = TRIM(zunit)
    mtgrm_data%sfc_var_info(ivar)%igroup_id = igroup_id
  END SUBROUTINE add_sfc_var

END MODULE mo_mtgrm_output

