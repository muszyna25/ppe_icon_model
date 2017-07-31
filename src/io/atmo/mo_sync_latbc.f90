!>
!! This module contains the I/O routines for lateral boundary nudging
!!
!! @author S. Brdar (DWD)
!!
!!
!! @par Revision History
!! Initial release by S. Brdar, DWD (2013-06-13)
!! Allow boundary data from the ICON output by S. Brdar, DWD (2013-07-19)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
!!
!!

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_sync_latbc
!-------------------------------------------------------------------------
!
!    ProTeX FORTRAN source: Style 2
!    modified for ICON project, DWD/MPI-M 2006
!
!-------------------------------------------------------------------------
!
!

  USE mo_kind,                ONLY: wp
  USE mo_parallel_config,     ONLY: nproma, p_test_run
  USE mo_time_config,         ONLY: time_config
  USE mo_model_domain,        ONLY: t_patch
  USE mo_grid_config,         ONLY: nroot
  USE mo_exception,           ONLY: message, message_text, finish
  USE mo_impl_constants,      ONLY: MAX_CHAR_LENGTH, MODE_COSMO, MODE_COMBINED, MODE_IFSANA, SUCCESS
  USE mo_impl_constants_grf,  ONLY: grf_bdywidth_c, grf_bdywidth_e
  USE mo_mpi,                 ONLY: p_io, p_bcast, my_process_is_stdio,       &
                                    p_comm_work_test, p_comm_work
  USE mo_io_units,            ONLY: filename_max
  USE mo_nonhydro_types,      ONLY: t_nh_state
  USE mo_intp_data_strc,      ONLY: t_int_state
  USE mo_nh_vert_interp,      ONLY: vert_interp
  USE mo_util_phys,           ONLY: virtual_temp
  USE mo_util_string,         ONLY: int2string
  USE mo_nh_init_utils,       ONLY: interp_uv_2_vn, convert_thdvars, convert_omega2w, &
    &                               compute_input_pressure_and_height
  USE mo_io_config,           ONLY: default_read_method
  USE mo_read_interface,      ONLY: nf, openInputFile, closeFile, read_3D, read_3D_1time,&
                                    read_2D_1lev_1time, read_2D_1time, t_stream_id, on_cells, &
                                    on_edges
  USE mo_sync,                ONLY: SYNC_E, SYNC_C, sync_patch_array
  USE mo_initicon_types,      ONLY: t_init_state, t_init_state_const
  USE mo_loopindices,         ONLY: get_indices_c, get_indices_e
  USE mo_limarea_config,      ONLY: latbc_config, generate_filename, LATBC_TYPE_EXT
  USE mo_ext_data_types,      ONLY: t_external_data
  USE mo_run_config,          ONLY: iqv, iqc, iqi, iqr, iqs, ltransport
  USE mo_initicon_config,     ONLY: init_mode
  USE mtime,                  ONLY: datetime, newDatetime, deallocateDatetime,              &
    &                               divideDatetimeDifferenceInSeconds, deallocateTimedelta, &
    &                               divisionquotienttimespan,                               &
    &                               timedelta, newTimedelta, OPERATOR(*), OPERATOR(+),      &
    &                               getTotalMilliSecondsTimeDelta, OPERATOR(>=)

  IMPLICIT NONE

  ! required for reading netcdf files
  INCLUDE 'netcdf.inc'

  PRIVATE


  LOGICAL                :: lread_qr, lread_qs  ! are qr, qs provided as input?
  LOGICAL                :: lread_vn            ! is vn provided as input?
  CHARACTER(LEN=10)      :: psvar
  CHARACTER(LEN=10)      :: geop_ml_var         ! model level surface geopotential
  INTEGER                :: latbc_fileid, &
                            read_latbc_tlev, &  ! time level indices for  p_latbc_data. can be 0 or 1.
                            last_latbc_tlev     ! last_ext_tlev is the last written time level index

  ! ---------------------------------------------------------------------------------
  TYPE(datetime), POINTER :: last_latbc_mtime    ! last read time step (mtime)
  TYPE(t_init_state)      :: p_latbc_data(2)     ! storage for two time-level boundary data

  ! storage for time-constant height level data
  TYPE(t_init_state_const), TARGET :: p_latbc_data_const


  PUBLIC :: prepare_latbc_data, read_latbc_data, deallocate_latbc_data,   &
    &       p_latbc_data, latbc_fileid, read_latbc_tlev, last_latbc_tlev, &
    &       update_lin_interc

  !> module name
  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_sync_latbc'

  CONTAINS


    SUBROUTINE allocate_latbc_data(p_patch, p_nh_state, ext_data, nlev_in)
      TYPE(t_patch),          INTENT(IN)   :: p_patch
      TYPE(t_nh_state),       INTENT(IN)   :: p_nh_state  !< nonhydrostatic state on the global domain
      TYPE(t_external_data),  INTENT(IN)   :: ext_data    !< external data on the global domain
      INTEGER, INTENT(IN) :: nlev_in !< number of vertical levels in the boundary data

      CHARACTER(LEN=*), PARAMETER :: routine = modname//"::allocate_latbc_data"
      INTEGER :: tlev, nlev, nlevp1, nblks_c, nblks_e, ierrstat
      LOGICAL :: l_initialized

      ! --- skip this routine if it has already been called ----
      l_initialized = p_latbc_data_const%linitialized
      IF (l_initialized)  RETURN
      ! --------------------------------------------------------

      nlev    = p_patch%nlev
      nlevp1  = p_patch%nlevp1
      nblks_c = p_patch%nblks_c
      nblks_e = p_patch%nblks_e

      ! Basic icon_remap data

      p_latbc_data_const%linitialized = .TRUE.
      ALLOCATE(p_latbc_data_const%z_mc_in(nproma,nlev_in,nblks_c), STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")        

      ! topography and metrics are time independent
      p_latbc_data_const%topography_c => ext_data%atm%topography_c
      p_latbc_data_const%z_ifc        => p_nh_state%metrics%z_ifc
      p_latbc_data_const%z_mc         => p_nh_state%metrics%z_mc

      DO tlev = 1, 2
        NULLIFY(p_latbc_data(tlev)%atm_in%tke)

        p_latbc_data(tlev)%atm_in%nlev = nlev_in
       
        ! Allocate atmospheric input data
        ALLOCATE(&
          p_latbc_data(tlev)%atm_in%pres (nproma,nlev_in,nblks_c), &
          p_latbc_data(tlev)%atm_in%temp (nproma,nlev_in,nblks_c), &
          p_latbc_data(tlev)%atm_in%u    (nproma,nlev_in,nblks_c), &
          p_latbc_data(tlev)%atm_in%v    (nproma,nlev_in,nblks_c), &
          p_latbc_data(tlev)%atm_in%w    (nproma,nlev_in,nblks_c), &
          p_latbc_data(tlev)%atm_in%qv   (nproma,nlev_in,nblks_c), &
          p_latbc_data(tlev)%atm_in%qc   (nproma,nlev_in,nblks_c), &
          p_latbc_data(tlev)%atm_in%qi   (nproma,nlev_in,nblks_c), &
          p_latbc_data(tlev)%atm_in%qr   (nproma,nlev_in,nblks_c), &
          p_latbc_data(tlev)%atm_in%qs   (nproma,nlev_in,nblks_c)  )
      
        ! Allocate atmospheric output data
        ALLOCATE(&
          p_latbc_data(tlev)%atm%temp     (nproma,nlev  ,nblks_c), &
          p_latbc_data(tlev)%atm%exner    (nproma,nlev  ,nblks_c), &
          p_latbc_data(tlev)%atm%pres     (nproma,nlev  ,nblks_c) )

        ALLOCATE(p_latbc_data(tlev)%atm%vn       (nproma,nlev  ,nblks_e), &
          p_latbc_data(tlev)%atm%u        (nproma,nlev  ,nblks_c), &
          p_latbc_data(tlev)%atm%v        (nproma,nlev  ,nblks_c), &
          p_latbc_data(tlev)%atm%w        (nproma,nlevp1,nblks_c))

        ALLOCATE(&
          p_latbc_data(tlev)%atm%rho      (nproma,nlev  ,nblks_c), &
          p_latbc_data(tlev)%atm%theta_v  (nproma,nlev  ,nblks_c), &
          p_latbc_data(tlev)%atm%qv       (nproma,nlev  ,nblks_c), &
          p_latbc_data(tlev)%atm%qc       (nproma,nlev  ,nblks_c), &
          p_latbc_data(tlev)%atm%qi       (nproma,nlev  ,nblks_c), &
          p_latbc_data(tlev)%atm%qr       (nproma,nlev  ,nblks_c), &
          p_latbc_data(tlev)%atm%qs       (nproma,nlev  ,nblks_c)  )
        
        ! allocate anyway (sometimes not needed)
        ALLOCATE(p_latbc_data(tlev)%atm_in%vn(nproma,nlev_in,p_patch%nblks_e))
              
        p_latbc_data(tlev)%const => p_latbc_data_const

      END DO
     
    END SUBROUTINE allocate_latbc_data

  !-------------------------------------------------------------------------
  !>
  !! @par Revision History
  !! Initial version by S. Brdar, DWD (2013-06-13)
  !!
  SUBROUTINE prepare_latbc_data(p_patch, p_int_state, p_nh_state, ext_data)
    TYPE(t_patch),          INTENT(IN)   :: p_patch
    TYPE(t_int_state),      INTENT(IN)   :: p_int_state
    TYPE(t_nh_state),       INTENT(INOUT):: p_nh_state  !< nonhydrostatic state on the global domain
    TYPE(t_external_data),  INTENT(IN)   :: ext_data    !< external data on the global domain

    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::allocate_latbc_data"
    INTEGER                              :: errno, step
    TYPE(divisionquotienttimespan)       :: tq     
    TYPE(timedelta), POINTER             :: td

    last_latbc_mtime    => newDatetime(time_config%tc_exp_startdate, errno)
    IF (errno /= 0)  CALL finish(routine, "Error in initialization of last_latbc_mtime.")

    ! last reading-in time is the current time
    ! mtime calculation: 
    !
    td => newTimedelta(latbc_config%dtime_latbc_mtime, errno)
    IF (errno /= 0)  CALL finish(routine, "Error in initialization of dtime_latbc time delta.")
    !   last_latbc_datetime := tc_startdate + dtime_latbc * FLOOR( (tc_current_date - tc_startdate)/dtime_latbc )

    CALL divideDatetimeDifferenceInSeconds(time_config%tc_current_date, time_config%tc_exp_startdate, &
      &                                    latbc_config%dtime_latbc_mtime, tq)

    step = INT(tq%quotient)
    td   = td * step
    last_latbc_mtime = last_latbc_mtime + td
    CALL deallocateTimedelta(td)

    ! prepare read/last indices
    read_latbc_tlev = 1   ! read in the first time-level slot
    last_latbc_tlev = 2

    ! read first two time steps
    CALL read_latbc_data( p_patch, p_nh_state, ext_data, p_int_state,           &
      &                   time_config%tc_current_date, lopt_check_read=.FALSE., &
      &                   lopt_time_incr=.FALSE.                              )
    CALL read_latbc_data( p_patch, p_nh_state, ext_data, p_int_state,           &
      &                   time_config%tc_current_date, lopt_check_read=.FALSE., &
      &                   lopt_time_incr=.TRUE.                               )

  END SUBROUTINE prepare_latbc_data
  !-------------------------------------------------------------------------



  !-------------------------------------------------------------------------
  !>
  !! Read horizontally interpolated atmospheric boundary data
  !!
  !! The subroutine reads atmospheric boundary data and projects on
  !! the ICON global grid
  !!
  !! @par Revision History
  !! Initial version by S. Brdar, DWD (2013-07-19)
  !!
  SUBROUTINE read_latbc_data( p_patch, p_nh_state, ext_data, p_int, mtime_date, &
    &                         lopt_check_read, lopt_time_incr )
    TYPE(t_patch),          INTENT(IN)    :: p_patch
    TYPE(t_nh_state),       INTENT(INOUT) :: p_nh_state  !< nonhydrostatic state on the global domain
    TYPE(t_external_data),  INTENT(IN)    :: ext_data    !< external data on the global domain
    TYPE(t_int_state),      INTENT(IN)    :: p_int
    TYPE(datetime),         pointer       :: mtime_date       !< current time (mtime format)
    LOGICAL,      INTENT(IN), OPTIONAL    :: lopt_check_read
    LOGICAL,      INTENT(IN), OPTIONAL    :: lopt_time_incr  !< increment latbc_datetime

    LOGICAL                               :: lcheck_read
    LOGICAL                               :: ltime_incr

    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::read_latbc_data"

    IF (PRESENT(lopt_check_read)) THEN
      lcheck_read = lopt_check_read
    ELSE
      lcheck_read = .TRUE.
    ENDIF

    IF (PRESENT(lopt_time_incr)) THEN
      ltime_incr = lopt_time_incr
    ELSE
      ltime_incr = .TRUE.
    ENDIF

    !
    ! do we need to read boundary data
    !
    IF (lcheck_read .AND. (last_latbc_mtime >= mtime_date)) RETURN

    ! Prepare the last_latbc_datetime for the next time level
    IF (ltime_incr) THEN
      last_latbc_mtime = last_latbc_mtime + latbc_config%dtime_latbc_mtime
    ENDIF

    ! Adjust read/last indices
    !
    ! New boundary data time-level is always read in p_latbc_data(read_latbc_tlev),
    ! whereas p_latbc_data(last_latbc_tlev) always holds the last read boundary data
    !
    read_latbc_tlev = last_latbc_tlev
    last_latbc_tlev = 3 - read_latbc_tlev

    !
    ! start reading boundary data
    !
    IF (latbc_config%itype_latbc == LATBC_TYPE_EXT) THEN
      CALL read_latbc_ifs_data(  p_patch, p_nh_state, ext_data, p_int )
    ELSE
      CALL read_latbc_icon_data( p_patch, p_nh_state, ext_data, p_int )
    ENDIF

    ! Compute tendencies for nest boundary update
    IF(ltime_incr) &
      CALL compute_boundary_tendencies(p_patch, p_nh_state)

  END SUBROUTINE read_latbc_data
  !-------------------------------------------------------------------------



  !-------------------------------------------------------------------------
  !>
  !! Read horizontally interpolated atmospheric ICON output
  !!
  !! @par Revision History
  !! Initial version by S. Brdar, DWD (2013-07-25)
  !!
  SUBROUTINE read_latbc_icon_data(p_patch, p_nh_state, ext_data, p_int)
    TYPE(t_patch), TARGET,  INTENT(IN)  :: p_patch
    TYPE(t_nh_state),       INTENT(IN)  :: p_nh_state  !< nonhydrostatic state on the global domain
    TYPE(t_external_data),  INTENT(IN)  :: ext_data    !< external data on the global domain
    TYPE(t_int_state),      INTENT(IN)  :: p_int

    ! local variables
    INTEGER                             :: dimid, no_cells, &
                                           latbc_ncid, no_levels
    TYPE(t_stream_id)                   :: latbc_stream_id
    LOGICAL                             :: l_exist
    REAL(wp)                            :: temp_v(nproma,p_patch%nlev,p_patch%nblks_c)
    INTEGER                             :: tlev, mpi_comm

    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::read_latbc_icon_data"
    CHARACTER(LEN=filename_max)         :: latbc_filename, latbc_full_filename

    ! NOTE: this subroutine is called only for
    !       (latbc_config%itype_latbc /= LATBC_TYPE_EXT). It is
    !       therefore implicitly assumed that no. of levels in input
    !       file matches the no. of levels in the current ICON run.

    tlev    = read_latbc_tlev

    IF(p_test_run) THEN
      mpi_comm = p_comm_work_test
    ELSE
      mpi_comm = p_comm_work
    ENDIF

    latbc_filename = generate_filename(nroot, p_patch%level, last_latbc_mtime)

    latbc_full_filename = TRIM(latbc_config%latbc_path)//TRIM(latbc_filename)

    IF (my_process_is_stdio()) THEN
      WRITE(message_text,'(a,a)') 'reading from ', TRIM(latbc_filename)
      CALL message(TRIM(routine), message_text)
      INQUIRE (FILE=TRIM(ADJUSTL(latbc_full_filename)), EXIST=l_exist)
      IF (.NOT. l_exist) THEN
        WRITE (message_text,'(a,a)') 'file not found: ', TRIM(latbc_filename)
        CALL finish(TRIM(routine), message_text)
      ENDIF

      !
      ! open file
      !
      CALL nf(nf_open(TRIM(ADJUSTL(latbc_full_filename)), NF_NOWRITE, latbc_ncid), routine)

      !
      ! get number of cells
      !
      CALL nf(nf_inq_dimid(latbc_ncid, 'ncells', dimid), routine)
      CALL nf(nf_inq_dimlen(latbc_ncid, dimid, no_cells), routine)

      !
      ! get number of vertical levels
      !
      CALL nf(nf_inq_dimid(latbc_ncid, 'height', dimid), routine)
      CALL nf(nf_inq_dimlen(latbc_ncid, dimid, no_levels), routine)

      ! consistency check
      IF (p_patch%nlev /= no_levels) THEN
        CALL finish(routine, "Number of levels does not match: "//int2string(no_levels, '(i0)')//&
          &" /= "//int2string(p_patch%nlev, '(i0)')//"!")
      END IF

      !
      ! check the number of cells and vertical levels
      !
      IF(p_patch%n_patch_cells_g /= no_cells) THEN
        WRITE(message_text,*) 'n_patch_cells_g does not match', &
          &   p_patch%n_patch_cells_g, ', ', no_cells
        CALL finish(TRIM(routine),TRIM(message_text))
      ENDIF

    END IF ! my_process_is_stdio()

    ! broadcast the no. of vertical input levels and allocate the data
    ! structure (if necessary):
    CALL p_bcast(no_levels, p_io, mpi_comm)
    CALL allocate_latbc_data(p_patch, p_nh_state, ext_data, no_levels)

    latbc_stream_id = openInputFile(latbc_full_filename, p_patch, &
      &                             default_read_method)
    !
    ! read prognostic 3d fields
    !
    CALL read_3D(stream_id=latbc_stream_id, location=on_cells, &
      &          variable_name='temp', fill_array=p_latbc_data(tlev)%atm%temp)
    CALL read_3D(stream_id=latbc_stream_id, location=on_cells, &
      &          variable_name='u', fill_array=p_latbc_data(tlev)%atm_in%u)
    CALL read_3D(stream_id=latbc_stream_id, location=on_cells, &
      &          variable_name='v', fill_array=p_latbc_data(tlev)%atm_in%v)
    CALL read_3D(stream_id=latbc_stream_id, location=on_cells, &
      &          variable_name='w', fill_array=p_latbc_data(tlev)%atm%w)
    CALL read_3D(stream_id=latbc_stream_id, location=on_cells, &
      &          variable_name='pres', fill_array=p_latbc_data(tlev)%atm%pres)
    CALL read_3D(stream_id=latbc_stream_id, location=on_cells, &
      &          variable_name='qv', fill_array=p_latbc_data(tlev)%atm%qv)
    CALL read_3D(stream_id=latbc_stream_id, location=on_cells, &
      &          variable_name='qc', fill_array=p_latbc_data(tlev)%atm%qc)
    CALL read_3D(stream_id=latbc_stream_id, location=on_cells, &
      &          variable_name='qi', fill_array=p_latbc_data(tlev)%atm%qi)
    CALL read_3D(stream_id=latbc_stream_id, location=on_cells, &
      &          variable_name='qr', fill_array=p_latbc_data(tlev)%atm%qr)
    CALL read_3D(stream_id=latbc_stream_id, location=on_cells, &
      &          variable_name='qs', fill_array=p_latbc_data(tlev)%atm%qs)
    !
    ! Convert u and v on cell points to vn at edge points
    !
    CALL interp_uv_2_vn( p_patch, p_int, p_latbc_data(tlev)%atm_in%u,                   &
      &                  p_latbc_data(tlev)%atm_in%v, p_latbc_data(tlev)%atm%vn )

    !
    ! Compute virtual temperature
    !
    CALL virtual_temp( p_patch, p_latbc_data(tlev)%atm%temp, p_latbc_data(tlev)%atm%qv, &
      &                p_latbc_data(tlev)%atm%qc, p_latbc_data(tlev)%atm%qi,            &
      &                p_latbc_data(tlev)%atm%qr, p_latbc_data(tlev)%atm%qs,            &
      &                temp_v=temp_v )

    !
    ! Compute NH prognostic thermodynamical variables
    !
    CALL convert_thdvars( p_patch, p_latbc_data(tlev)%atm%pres, temp_v,                 &
      &                   p_latbc_data(tlev)%atm%rho,                                   &
      &                   p_latbc_data(tlev)%atm%exner,                                 &
      &                   p_latbc_data(tlev)%atm%theta_v )

    CALL sync_patch_array(SYNC_E, p_patch, p_latbc_data(tlev)%atm%vn)
    CALL sync_patch_array(SYNC_C, p_patch, p_latbc_data(tlev)%atm%w)
    CALL sync_patch_array(SYNC_C, p_patch, p_latbc_data(tlev)%atm%theta_v)
    CALL sync_patch_array(SYNC_C, p_patch, p_latbc_data(tlev)%atm%rho)

    !
    ! close file
    !
    IF (my_process_is_stdio()) THEN
      WRITE(message_text,'(a,a)') 'closing file ', TRIM(latbc_filename)
      CALL message(TRIM(routine), message_text)
      CALL nf(nf_close(latbc_ncid), routine)
    END IF
    CALL closeFile(latbc_stream_id)

  END SUBROUTINE read_latbc_icon_data
  !-------------------------------------------------------------------------



  !-------------------------------------------------------------------------
  !>
  !! Read horizontally interpolated atmospheric IFS analysis
  !!
  !! This subroutine is a simplification of mo_nh_initicons:process_ifsana_atm.
  !! The following steps are performed:
  !! - read atmospheric IFS analysis data,
  !! - interpolate vertically from intermediate IFS2ICON grid to ICON grid
  !!   and compute the prognostic NH variable set,
  !!  !! @par Revision History
  !! Initial version by S. Brdar, DWD (2013-06-13)
  !!
  SUBROUTINE read_latbc_ifs_data(p_patch, p_nh_state, ext_data, p_int)
    TYPE(t_patch),          INTENT(IN)  :: p_patch
    TYPE(t_nh_state),       INTENT(IN)  :: p_nh_state  !< nonhydrostatic state on the global domain
    TYPE(t_external_data),  INTENT(IN)  :: ext_data    !< external data on the global domain
    TYPE(t_int_state),      INTENT(IN)  :: p_int

    ! local variables
    INTEGER                             :: mpi_comm, dimid, no_cells, &
                                           latbc_ncid, no_levels, varid
    TYPE(t_stream_id)                   :: latbc_stream_id
    LOGICAL                             :: l_exist, lconvert_omega2w
    INTEGER                             :: jc, jk, jb, i_endidx

    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::read_latbc_ifs_data"
    CHARACTER(LEN=filename_max)         :: latbc_filename, latbc_full_filename
    INTEGER                             :: tlev, nlev, ierrstat, nblks_c
    REAL(wp), ALLOCATABLE               :: psfc(:,:), phi_sfc(:,:), z_ifc_in(:,:,:), &
      &                                    w_ifc(:,:,:), omega(:,:,:)

    tlev      = read_latbc_tlev
    nblks_c   = p_patch%nblks_c

    IF(p_test_run) THEN
      mpi_comm = p_comm_work_test
    ELSE
      mpi_comm = p_comm_work
    ENDIF

    latbc_filename = generate_filename(nroot, p_patch%level, last_latbc_mtime)

    latbc_full_filename = TRIM(latbc_config%latbc_path)//TRIM(latbc_filename)

    IF (my_process_is_stdio()) THEN
      WRITE(message_text,'(a,a)') 'reading boundary data: ', TRIM(latbc_filename)
      CALL message(TRIM(routine), message_text)
      INQUIRE (FILE=TRIM(ADJUSTL(latbc_full_filename)), EXIST=l_exist)
      IF (.NOT. l_exist) THEN
        WRITE (message_text,'(a,a)') 'file not found:', TRIM(latbc_filename)
        CALL finish(TRIM(routine), message_text)
      ENDIF

      !
      ! open file
      !
      CALL nf(nf_open(TRIM(ADJUSTL(latbc_full_filename)), NF_NOWRITE, latbc_ncid), routine)

      !
      ! get number of cells
      !
      CALL nf(nf_inq_dimid(latbc_ncid, 'ncells', dimid), routine)
      CALL nf(nf_inq_dimlen(latbc_ncid, dimid, no_cells), routine)

      !
      ! get number of vertical levels
      !
      CALL nf(nf_inq_dimid(latbc_ncid, 'lev', dimid), routine)
      CALL nf(nf_inq_dimlen(latbc_ncid, dimid, no_levels), routine)

      !
      ! check the number of cells
      !
      IF(p_patch%n_patch_cells_g /= no_cells) THEN
        WRITE(message_text,*) 'n_patch_cells_g does not match', &
          &   p_patch%n_patch_cells_g, ', ', no_cells
        CALL finish(TRIM(routine),TRIM(message_text))
      ENDIF

      !
      ! Check if surface pressure (PS) or its logarithm (LNPS) is provided as input
      !
      IF (nf_inq_varid(latbc_ncid, 'PS', varid) == nf_noerr) THEN
        psvar = 'PS'
      ELSE IF (nf_inq_varid(latbc_ncid, 'LNPS', varid) == nf_noerr) THEN
        psvar = 'LNPS'
      ENDIF

      !
      ! Check if model-level surface Geopotential is provided as GEOSP or GEOP_ML
      !
      IF (nf_inq_varid(latbc_ncid, 'GEOSP', varid) == nf_noerr) THEN
        geop_ml_var = 'GEOSP'
      ELSE IF (nf_inq_varid(latbc_ncid, 'GEOP_ML', varid) == nf_noerr) THEN
        geop_ml_var = 'GEOP_ML'
      ELSE
        CALL finish(TRIM(routine),'Could not find model-level sfc geopotential')
      ENDIF

      !
      ! Check if rain water (QR) is provided as input
      !
      IF (nf_inq_varid(latbc_ncid, 'QR', varid) == nf_noerr) THEN
        lread_qr = .true.
      ELSE
        lread_qr = .false.
        CALL message(TRIM(routine),'Rain water (QR) not available in input data')
      ENDIF

      !
      ! Check if snow water (QS) is provided as input
      !
      IF (nf_inq_varid(latbc_ncid, 'QS', varid) == nf_noerr) THEN
        lread_qs = .true.
      ELSE
        lread_qs = .false.
        CALL message(TRIM(routine),'Snow water (QS) not available in input data')
      ENDIF

      !
      ! Check if normal wind (VN) is provided as input
      !
      IF (nf_inq_varid(latbc_ncid, 'VN', varid) == nf_noerr) THEN
        lread_vn = .TRUE.
      ELSE
        lread_vn = .FALSE.
      ENDIF

    END IF ! my_process_is_stdio()

    CALL p_bcast(lread_qs, p_io, mpi_comm)
    CALL p_bcast(lread_qr, p_io, mpi_comm)
    CALL p_bcast(lread_vn, p_io, mpi_comm)

    ! broadcast the no. of vertical input levels and allocate the data
    ! structure (if necessary):
    CALL p_bcast(no_levels, p_io, mpi_comm)
    CALL allocate_latbc_data(p_patch, p_nh_state, ext_data, no_levels)

    latbc_stream_id = openInputFile(latbc_full_filename, p_patch, &
      &                             default_read_method)

    !
    ! read IFS data
    !
    CALL read_3D_1time(stream_id=latbc_stream_id, location=on_cells, &
      &          variable_name='T', fill_array=p_latbc_data(tlev)%atm_in%temp)

    IF (lread_vn) THEN
      CALL read_3D_1time(stream_id=latbc_stream_id, location=on_edges, &
        &          variable_name='VN', fill_array= p_latbc_data(tlev)%atm_in%vn)
    ELSE
      CALL read_3D_1time(stream_id=latbc_stream_id, location=on_cells, &
        &          variable_name='U', fill_array=p_latbc_data(tlev)%atm_in%u)
      CALL read_3D_1time(stream_id=latbc_stream_id, location=on_cells, &
        &          variable_name='V', fill_array=p_latbc_data(tlev)%atm_in%v)
    ENDIF

    IF (init_mode /= MODE_COSMO) THEN
      lconvert_omega2w = .TRUE.

      ! allocate temporary array:
      nlev = SIZE(p_latbc_data_const%z_mc_in,2)
      ALLOCATE(omega(nproma, nlev, nblks_c), STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")

      CALL read_3D_1time(stream_id=latbc_stream_id, location=on_cells, &
        &          variable_name='W', fill_array=omega)
    ELSE
      lconvert_omega2w = .FALSE.

      ! allocate temporary array:
      nlev = SIZE(p_latbc_data_const%z_mc_in,2)
      ALLOCATE(w_ifc(nproma,    (nlev+1), nblks_c), STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")

      CALL read_3D_1time(stream_id=latbc_stream_id, location=on_cells, &
        &          variable_name='W', fill_array=w_ifc)
    ENDIF

    IF (init_mode == MODE_COSMO) THEN
      ! allocate temporary array:
      nlev = SIZE(p_latbc_data_const%z_mc_in,2)
      ALLOCATE(z_ifc_in(nproma, (nlev+1), nblks_c), STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")

      IF (nlev /= no_levels) THEN
        WRITE (message_text, *) "Non-matching no. of levels: nlev=", nlev, "; no_levels=",no_levels
        CALL finish(routine, message_text)
      END IF

      CALL read_3D_1time(stream_id=latbc_stream_id, location=on_cells, &
        &          variable_name='HHL', fill_array=z_ifc_in)

      ! Interpolate input 'z3d' and 'w' from the interface levels to the main levels
      !
!$OMP PARALLEL
!$OMP DO PRIVATE (jk,jc,jb,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = 1,p_patch%nblks_c

        IF (jb /= p_patch%nblks_c) THEN
          i_endidx = nproma
        ELSE
          i_endidx = p_patch%npromz_c
        ENDIF

        DO jk = 1, no_levels
          DO jc = 1, i_endidx

        ! Note: In future, we want to z3d from boundary data.
        !
        p_latbc_data_const%z_mc_in(jc,jk,jb)  = (z_ifc_in(jc,jk,jb) + z_ifc_in(jc,jk+1,jb)) * 0.5_wp
        p_latbc_data(tlev)%atm_in%w(jc,jk,jb) = (w_ifc(jc,jk,jb)    + w_ifc(jc,jk+1,jb)) * 0.5_wp

      ENDDO
    ENDDO
  ENDDO
!$OMP END DO
!$OMP END PARALLEL

    ENDIF

    CALL read_3D_1time(stream_id=latbc_stream_id, location=on_cells, &
      &          variable_name='QV', fill_array=p_latbc_data(tlev)%atm_in%qv)
    CALL read_3D_1time(stream_id=latbc_stream_id, location=on_cells, &
      &          variable_name='QC', fill_array=p_latbc_data(tlev)%atm_in%qc)
    CALL read_3D_1time(stream_id=latbc_stream_id, location=on_cells, &
      &          variable_name='QI', fill_array=p_latbc_data(tlev)%atm_in%qi)
    IF (lread_qr) THEN
      CALL read_3D_1time(stream_id=latbc_stream_id, location=on_cells, &
        &          variable_name='QR', fill_array=p_latbc_data(tlev)%atm_in%qr)
    ELSE
      p_latbc_data(tlev)%atm_in%qr(:,:,:)=0._wp
    ENDIF

    IF (lread_qs) THEN
      CALL read_3D_1time(stream_id=latbc_stream_id, location=on_cells, &
        &          variable_name='QS', fill_array=p_latbc_data(tlev)%atm_in%qs)
    ELSE
      p_latbc_data(tlev)%atm_in%qs(:,:,:)=0._wp
    ENDIF

    ! allocate local temporary arrays:
    ALLOCATE(psfc(nproma, p_patch%nblks_c), phi_sfc(nproma, p_patch%nblks_c))

    IF (init_mode == MODE_COSMO) THEN
      CALL read_2D_1time(stream_id=latbc_stream_id, location=on_cells, &
       &                     variable_name=TRIM(psvar), fill_array=psfc)
      CALL read_2D_1time(stream_id=latbc_stream_id, location=on_cells, &
        &                     variable_name=TRIM(geop_ml_var), fill_array=phi_sfc)
    ELSE
      CALL read_2D_1lev_1time(stream_id=latbc_stream_id, location=on_cells, &
        &                     variable_name=TRIM(psvar), fill_array=psfc)
      CALL read_2D_1lev_1time(stream_id=latbc_stream_id, location=on_cells, &
        &                     variable_name=TRIM(geop_ml_var), fill_array=phi_sfc)
    END IF

    ! compute pressure and height of input data, using the IFS routines
    IF ((init_mode == MODE_IFSANA) .OR. (init_mode == MODE_COMBINED)) THEN ! i.e. atmospheric data from IFS
      CALL compute_input_pressure_and_height(p_patch, psfc, phi_sfc, p_latbc_data(tlev))
    END IF

    IF (init_mode == MODE_COSMO) THEN
      CALL read_3D_1time(stream_id=latbc_stream_id, location=on_cells, &
        &          variable_name='P', fill_array=p_latbc_data(tlev)%atm_in%pres)
    ENDIF

    !
    ! close file
    !
    IF (my_process_is_stdio()) THEN
      WRITE(message_text,'(a,a)') 'closing file ', TRIM(latbc_filename)
      CALL message(TRIM(routine), message_text)
      CALL nf(nf_close(latbc_ncid), routine)
    END IF
    CALL closeFile(latbc_stream_id)

    IF (lconvert_omega2w) THEN
      CALL convert_omega2w(omega, &
        &                  p_latbc_data(tlev)%atm_in%w,     &
        &                  p_latbc_data(tlev)%atm_in%pres,  &
        &                  p_latbc_data(tlev)%atm_in%temp,  &
        &                  p_patch%nblks_c, p_patch%npromz_c,   &
        &                  p_latbc_data(tlev)%atm_in%nlev )
    END IF

    ! cleanup
    IF (ALLOCATED(omega))    DEALLOCATE(omega)
    IF (ALLOCATED(psfc))     DEALLOCATE(psfc)
    IF (ALLOCATED(phi_sfc))  DEALLOCATE(phi_sfc)
    IF (ALLOCATED(z_ifc_in)) DEALLOCATE(z_ifc_in)
    IF (ALLOCATED(w_ifc))    DEALLOCATE(w_ifc)

    !
    ! perform vertical interpolation of horizonally interpolated analysis data
    !
    CALL vert_interp(p_patch, p_int, p_nh_state%metrics, p_latbc_data(tlev), opt_use_vn=lread_vn)

  END SUBROUTINE read_latbc_ifs_data
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !>
  !! @par Revision History
  !! Initial version by S. Brdar, DWD (2013-06-13)
  !!
  SUBROUTINE deallocate_latbc_data()
    INTEGER             :: tlev
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::deallocate_latbc_data"

    WRITE(message_text,'(a,a)') 'deallocating latbc data'
    CALL message(TRIM(routine), message_text)

    !
    ! deallocate boundary data memory
    !
    DO tlev = 1, 2
      CALL p_latbc_data(tlev)%finalize()
    END DO

    CALL p_latbc_data_const%finalize()

    ! clean up mtime calendar objects
    CALL deallocateDatetime(last_latbc_mtime)

  END SUBROUTINE deallocate_latbc_data
  !-------------------------------------------------------------------------



  !-------------------------------------------------------------------------
  !>
  !! @par Revision History
  !! Initial version by G. Zaengl, DWD (2013-10-22)
  !!
  SUBROUTINE compute_boundary_tendencies ( p_patch, p_nh )
    TYPE(t_patch),    INTENT(IN)    :: p_patch
    TYPE(t_nh_state), INTENT(INOUT) :: p_nh

    ! Local variables
    INTEGER                         :: i_startblk, i_endblk, &
      &                                i_startidx, i_endidx
    INTEGER                         :: jc, jk, jb, je
    INTEGER                         :: nlev, nlevp1
    REAL(wp)                        :: rdt

    nlev = p_patch%nlev
    nlevp1 = p_patch%nlevp1

    ! Inverse value of boundary update frequency
    rdt = 1._wp/latbc_config%dtime_latbc

!$OMP PARALLEL PRIVATE(i_startblk,i_endblk)

    ! a) Boundary tendency of horizontal velocity
    i_startblk = p_patch%edges%start_blk(1,1)
    i_endblk   = p_patch%edges%end_blk(grf_bdywidth_e,1)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk,je) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, 1, grf_bdywidth_e)

      DO jk = 1, nlev
        DO je = i_startidx, i_endidx
          p_nh%diag%grf_tend_vn(je,jk,jb) = rdt * (            &
            &   p_latbc_data(read_latbc_tlev)%atm%vn(je,jk,jb) &
            & - p_latbc_data(last_latbc_tlev)%atm%vn(je,jk,jb) )
        ENDDO
      ENDDO
    ENDDO
!$OMP END DO

    ! b) Boundary tendencies of variables at cell centers
    i_startblk = p_patch%cells%start_blk(1,1)
    i_endblk   = p_patch%cells%end_blk(grf_bdywidth_c,1)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk,jc) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, 1, grf_bdywidth_c)

      DO jk = 1, nlev
        DO jc = i_startidx, i_endidx

          p_nh%diag%grf_tend_rho(jc,jk,jb) = rdt * (            &
            &   p_latbc_data(read_latbc_tlev)%atm%rho(jc,jk,jb) &
            & - p_latbc_data(last_latbc_tlev)%atm%rho(jc,jk,jb) )

          p_nh%diag%grf_tend_thv(jc,jk,jb) = rdt * (                &
            &   p_latbc_data(read_latbc_tlev)%atm%theta_v(jc,jk,jb) &
            & - p_latbc_data(last_latbc_tlev)%atm%theta_v(jc,jk,jb) )

          p_nh%diag%grf_tend_w(jc,jk,jb) = rdt * (            &
            &   p_latbc_data(read_latbc_tlev)%atm%w(jc,jk,jb) &
            & - p_latbc_data(last_latbc_tlev)%atm%w(jc,jk,jb) )

        ENDDO
      ENDDO

      DO jc = i_startidx, i_endidx
        p_nh%diag%grf_tend_w(jc,nlevp1,jb) = rdt * (            &
          &   p_latbc_data(read_latbc_tlev)%atm%w(jc,nlevp1,jb) &
          & - p_latbc_data(last_latbc_tlev)%atm%w(jc,nlevp1,jb) )
      ENDDO

      IF (ltransport) THEN
        DO jk = 1, nlev
          DO jc = i_startidx, i_endidx
            p_nh%diag%grf_tend_tracer(jc,jk,jb,iqv) =  rdt * (   &
              &   p_latbc_data(read_latbc_tlev)%atm%qv(jc,jk,jb) &
              & - p_latbc_data(last_latbc_tlev)%atm%qv(jc,jk,jb) )

            p_nh%diag%grf_tend_tracer(jc,jk,jb,iqc) =  rdt * (   &
              &   p_latbc_data(read_latbc_tlev)%atm%qc(jc,jk,jb) &
              & - p_latbc_data(last_latbc_tlev)%atm%qc(jc,jk,jb) )

            p_nh%diag%grf_tend_tracer(jc,jk,jb,iqi) =  rdt * (   &
              &   p_latbc_data(read_latbc_tlev)%atm%qi(jc,jk,jb) &
              & - p_latbc_data(last_latbc_tlev)%atm%qi(jc,jk,jb) )

            p_nh%diag%grf_tend_tracer(jc,jk,jb,iqr) =  rdt * (   &
              &   p_latbc_data(read_latbc_tlev)%atm%qr(jc,jk,jb) &
              & - p_latbc_data(last_latbc_tlev)%atm%qr(jc,jk,jb) )

            p_nh%diag%grf_tend_tracer(jc,jk,jb,iqs) =  rdt * (   &
              &   p_latbc_data(read_latbc_tlev)%atm%qs(jc,jk,jb) &
              & - p_latbc_data(last_latbc_tlev)%atm%qs(jc,jk,jb) )

          ENDDO
        ENDDO
      ENDIF

    ENDDO
!$OMP END DO
!$OMP END PARALLEL


  END SUBROUTINE compute_boundary_tendencies
  !-------------------------------------------------------------------------



  !-------------------------------------------------------------------------
  !>
  !  Update linear interpolation coefficients for a given time stamp
  !
  !! @par Revision History
  !! Initial version by S. Brdar, DWD (2013-08-02)
  !!
  SUBROUTINE update_lin_interc( mtime_date )
    TYPE(datetime),     POINTER    :: mtime_date
    ! local variables
    CHARACTER(MAX_CHAR_LENGTH), PARAMETER :: routine = modname//"::update_lin_interc"
    TYPE(divisionquotienttimespan)       :: tq     
    REAL(wp)                             :: dtime_latbc_in_ms
#ifdef _MTIME_DEBUG
    REAL(wp)                             :: d1
#endif

    ! calendar calculation: mtime data structure
    !
    ! lc1 := (last_latbc_datetime - this_datetime) / dtime_latbc
    ! lc2 := 1 - lc1
    !
    ! We have 0 <= lc1 <= 1, because "last_latbc_datetime" (despite
    ! its name) denotes the *next* lateral bc read event.
    !
    dtime_latbc_in_ms = getTotalMilliSecondsTimeDelta(latbc_config%dtime_latbc_mtime, mtime_date)
    CALL divideDatetimeDifferenceInSeconds(last_latbc_mtime, mtime_date, latbc_config%dtime_latbc_mtime, tq)
    latbc_config%lc1 = REAL(tq%remainder_in_ms,wp)/dtime_latbc_in_ms
    latbc_config%lc2 = 1._wp - latbc_config%lc1

  END SUBROUTINE update_lin_interc
  !-------------------------------------------------------------------------


END MODULE mo_sync_latbc
