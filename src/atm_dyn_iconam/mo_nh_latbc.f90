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
!! @par Copyright
!! 2002-2013 by DWD and MPI-M
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

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_nh_latbc
!-------------------------------------------------------------------------
!
!    ProTeX FORTRAN source: Style 2
!    modified for ICON project, DWD/MPI-M 2006
!
!-------------------------------------------------------------------------
!
!

  USE mo_kind,                ONLY: wp, i8
  USE mo_parallel_config,     ONLY: nproma, p_test_run
  USE mo_model_domain,        ONLY: t_patch
  USE mo_grid_config,         ONLY: nroot
  USE mo_exception,           ONLY: message, message_text, finish
  USE mo_impl_constants,      ONLY: SUCCESS, MAX_CHAR_LENGTH
  USE mo_mpi,                 ONLY: p_io, p_bcast, my_process_is_stdio,       &
                                    p_comm_work_test, p_comm_work
  USE mo_io_units,            ONLY: filename_max
  USE mo_util_string,         ONLY: MAX_STRING_LEN
  USE mo_nonhydro_types,      ONLY: t_nh_state
  USE mo_intp_data_strc,      ONLY: t_int_state
  USE mo_nh_vert_interp,      ONLY: vert_interp
  USE mo_util_phys,           ONLY: virtual_temp
  USE mo_nh_init_utils,       ONLY: interp_uv_2_vn, convert_thdvars
  USE mo_sync,                ONLY: sync_patch_array, SYNC_E, SYNC_C
#ifdef NOMPI
  USE mo_mpi,                 ONLY: my_process_is_mpi_all_seq
#endif
  USE mo_netcdf_read,         ONLY: nf, read_netcdf_data_single,        &
                                    read_netcdf_data
  USE mo_nh_initicon_types,   ONLY: t_initicon_state
  USE mo_datetime,            ONLY: t_datetime
  USE mo_time_config,         ONLY: time_config
  USE mo_initicon_config,     ONLY: nlev_in
  USE mo_limarea_config,      ONLY: latbc_config, generate_filename
  USE mo_ext_data_types,      ONLY: t_external_data
  
  IMPLICIT NONE
  
  ! required for reading netcdf files
  INCLUDE 'netcdf.inc'

  PRIVATE

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

  LOGICAL                :: lread_qr, lread_qs  ! are qr, qs provided as input?
  LOGICAL                :: lread_vn            ! is vn provided as input?
  CHARACTER(LEN=10)      :: psvar
  CHARACTER(LEN=10)      :: geop_ml_var         ! model level surface geopotential
  INTEGER                :: latbc_fileid, &
                            read_latbc_tlev, &  ! time level indices for  p_latbc_data. can be 0 or 1.
                            last_latbc_tlev     ! last_ext_tlev is the last written time level index
  TYPE(t_datetime)       :: last_latbc_datetime ! last read time step
  TYPE(t_initicon_state) :: p_latbc_data(2)     ! storage for two time-level boundary data

  PUBLIC :: allocate_latbc_data, read_latbc_data, &
            deallocate_latbc_data, p_latbc_data, &
            latbc_fileid, read_latbc_tlev, last_latbc_tlev, last_latbc_datetime

  CONTAINS

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !>
  !! @par Revision History
  !! Initial version by S. Brdar, DWD (2013-06-13)
  !!
  SUBROUTINE allocate_latbc_data(patch, p_nh_state, ext_data)
    TYPE(t_patch),          INTENT(IN)  :: patch
    TYPE(t_nh_state),       INTENT(IN)  :: p_nh_state  !< nonhydrostatic state on the global domain
    TYPE(t_external_data),  INTENT(IN)  :: ext_data    !< external data on the global domain

    ! local variables
    INTEGER       :: mpi_comm
    INTEGER       :: tlev
    INTEGER       :: nlev, nlevp1, nblks_c, nblks_v, nblks_e
    CHARACTER(MAX_CHAR_LENGTH), PARAMETER :: routine = &
      "mo_nh_latbc::allocate_latbc_data"

    CALL message(TRIM(routine),'start')
    IF (nlev_in == 0) THEN
      CALL finish(TRIM(routine), "Number of input levels <nlev_in> not yet initialized.")
    END IF
     
    IF(p_test_run) THEN
      mpi_comm = p_comm_work_test
    ELSE
      mpi_comm = p_comm_work
    ENDIF

    last_latbc_datetime = time_config%ini_datetime

    nlev = patch%nlev
    nlevp1 = nlev + 1
    nblks_c = patch%nblks_c
    nblks_v = patch%nblks_v
    nblks_e = patch%nblks_e

    DO tlev = 1, 2
      ! Basic icon_remap data
      ALLOCATE(p_latbc_data(tlev)%topography_c(nproma,nblks_c),         &
               p_latbc_data(tlev)%topography_v(nproma,nblks_v),         &
               p_latbc_data(tlev)%z_ifc       (nproma,nlevp1,nblks_c),  &
               p_latbc_data(tlev)%z_mc        (nproma,nlev  ,nblks_c)   )

      ! Allocate atmospheric input data
      ALLOCATE(p_latbc_data(tlev)%atm_in%psfc(nproma,         nblks_c), &
               p_latbc_data(tlev)%atm_in%phi_sfc(nproma,      nblks_c), &
               p_latbc_data(tlev)%atm_in%pres (nproma,nlev_in,nblks_c), &
               p_latbc_data(tlev)%atm_in%z3d  (nproma,nlev_in,nblks_c), &
               p_latbc_data(tlev)%atm_in%temp (nproma,nlev_in,nblks_c), &
               p_latbc_data(tlev)%atm_in%u    (nproma,nlev_in,nblks_c), &
               p_latbc_data(tlev)%atm_in%v    (nproma,nlev_in,nblks_c), &
               p_latbc_data(tlev)%atm_in%w    (nproma,nlev_in,nblks_c), &
               p_latbc_data(tlev)%atm_in%omega(nproma,nlev_in,nblks_c), &
               p_latbc_data(tlev)%atm_in%qv   (nproma,nlev_in,nblks_c), &
               p_latbc_data(tlev)%atm_in%qc   (nproma,nlev_in,nblks_c), &
               p_latbc_data(tlev)%atm_in%qi   (nproma,nlev_in,nblks_c), &
               p_latbc_data(tlev)%atm_in%qr   (nproma,nlev_in,nblks_c), &
               p_latbc_data(tlev)%atm_in%qs   (nproma,nlev_in,nblks_c)  )

      ! Allocate atmospheric output data
      ALLOCATE(p_latbc_data(tlev)%atm%vn       (nproma,nlev  ,nblks_e), &
               p_latbc_data(tlev)%atm%u        (nproma,nlev  ,nblks_c), &
               p_latbc_data(tlev)%atm%v        (nproma,nlev  ,nblks_c), &
               p_latbc_data(tlev)%atm%w        (nproma,nlevp1,nblks_c), &
               p_latbc_data(tlev)%atm%temp     (nproma,nlev  ,nblks_c), &
               p_latbc_data(tlev)%atm%exner    (nproma,nlev  ,nblks_c), &
               p_latbc_data(tlev)%atm%pres     (nproma,nlev  ,nblks_c), &
               p_latbc_data(tlev)%atm%rho      (nproma,nlev  ,nblks_c), &
               p_latbc_data(tlev)%atm%theta_v  (nproma,nlev  ,nblks_c), &
               p_latbc_data(tlev)%atm%qv       (nproma,nlev  ,nblks_c), &
               p_latbc_data(tlev)%atm%qc       (nproma,nlev  ,nblks_c), &
               p_latbc_data(tlev)%atm%qi       (nproma,nlev  ,nblks_c), &
               p_latbc_data(tlev)%atm%qr       (nproma,nlev  ,nblks_c), &
               p_latbc_data(tlev)%atm%qs       (nproma,nlev  ,nblks_c)  )
               
      ! allocate anyway (sometimes not needed)
      ALLOCATE(p_latbc_data(tlev)%atm_in%vn(nproma,nlev_in,patch%nblks_e))

      ! topography and metrics are time independent
!$OMP PARALLEL
!$OMP WORKSHARE
      p_latbc_data(tlev)%topography_c(:,:) = ext_data%atm%topography_c(:,:)
      p_latbc_data(tlev)%topography_v(:,:) = ext_data%atm%topography_v(:,:)
      p_latbc_data(tlev)%z_ifc(:,:,:) = p_nh_state%metrics%z_ifc(:,:,:)
      p_latbc_data(tlev)%z_mc (:,:,:) = p_nh_state%metrics%z_mc (:,:,:) 
!$OMP END WORKSHARE
!$OMP END PARALLEL

    END DO

    CALL message(TRIM(routine),'done')

  END SUBROUTINE allocate_latbc_data
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
  SUBROUTINE read_latbc_data(p_patch, p_nh_state, p_int, ext_data, latbc_datetime, tlev)
    TYPE(t_patch), TARGET,  INTENT(IN)  :: p_patch(:)
    TYPE(t_nh_state),       INTENT(IN)  :: p_nh_state  !< nonhydrostatic state on the global domain
    TYPE(t_int_state),      INTENT(IN)  :: p_int
    TYPE(t_external_data),  INTENT(IN)  :: ext_data    !< external data on the global domain
    TYPE(t_datetime),       INTENT(IN)  :: latbc_datetime
    INTEGER,                INTENT(IN)  :: tlev        !< latbc data slot id used for storing time
                                                       !  level boundary data (values are 1 or 2)
    CHARACTER(MAX_CHAR_LENGTH), PARAMETER :: routine = "mo_nh_latbc::read_latbc_data"
                                                       
    SELECT CASE (latbc_config%itype_latbc)
    CASE(1)
      CALL message(TRIM(routine), 'IFS boundary data')
      CALL read_latbc_ifs_data(  p_patch, p_nh_state, p_int, ext_data, latbc_datetime, tlev)
    CASE(2)
      CALL message(TRIM(routine), 'ICON output boundary data')
      CALL read_latbc_icon_data( p_patch, p_nh_state, p_int, ext_data, latbc_datetime, tlev)
    END SELECT
  END SUBROUTINE read_latbc_data
  !-------------------------------------------------------------------------


  
  !-------------------------------------------------------------------------
  !>
  !! Read horizontally interpolated atmospheric ICON output
  !!
  !! @par Revision History
  !! Initial version by S. Brdar, DWD (2013-07-25)
  !!
  SUBROUTINE read_latbc_icon_data(p_patch, p_nh_state, p_int, ext_data, latbc_datetime, tlev)
    TYPE(t_patch), TARGET,  INTENT(IN)  :: p_patch(:)
    TYPE(t_nh_state),       INTENT(IN)  :: p_nh_state  !< nonhydrostatic state on the global domain
    TYPE(t_int_state),      INTENT(IN)  :: p_int
    TYPE(t_external_data),  INTENT(IN)  :: ext_data    !< external data on the global domain
    TYPE(t_datetime),       INTENT(IN)  :: latbc_datetime
    INTEGER,                INTENT(IN)  :: tlev        !< latbc data slot id used for storing time
                                                       !  level boundary data (values are 1 or 2)

    ! local variables
    TYPE(t_patch), POINTER              :: patch
    INTEGER                             :: rl_start, rl_end
    INTEGER                             :: jc,je,jk,jb                ! loop indices
    INTEGER                             :: nblks_c, nblks_e           ! number of blocks
    INTEGER                             :: i_startblk, i_endblk
    INTEGER                             :: i_startidx, i_endidx
    INTEGER                             :: mpi_comm, ist, dimid, no_cells, &
                                           latbc_fileid, no_levels
    LOGICAL                             :: l_exist
    INTEGER, POINTER                    :: iidx(:,:,:), iblk(:,:,:)
    REAL(wp)                            :: temp_v(nproma,p_patch(1)%nlev,p_patch(1)%nblks_c), &
      &                                    vn, w, rho, theta_v

    CHARACTER(MAX_CHAR_LENGTH), PARAMETER :: routine = "mo_nh_latbc::read_latbc_data"
    CHARACTER(LEN=filename_max)           :: latbc_full_filename

    patch          => p_patch(1)                                        
    iidx           => patch%edges%cell_idx
    iblk           => patch%edges%cell_blk

    IF (nlev_in == 0) THEN
      CALL finish(TRIM(routine), "Number of input levels <nlev_in> not yet initialized.")
    END IF
      
    IF(p_test_run) THEN
      mpi_comm = p_comm_work_test
    ELSE
      mpi_comm = p_comm_work
    ENDIF

    IF (my_process_is_stdio()) THEN
      latbc_full_filename = generate_filename(nroot, patch%level, latbc_datetime)

      WRITE(message_text,'(a,a)') 'reading from ', TRIM(latbc_full_filename)
      CALL message(TRIM(routine), message_text)
      INQUIRE (FILE=TRIM(ADJUSTL(latbc_full_filename)), EXIST=l_exist)
      IF (.NOT. l_exist) THEN
        WRITE (message_text,'(a,a)') 'file not found:', TRIM(latbc_full_filename)
        CALL finish(TRIM(routine), message_text)
      ENDIF

      !
      ! open file
      !
      CALL nf(nf_open(TRIM(ADJUSTL(latbc_full_filename)), NF_NOWRITE, latbc_fileid), routine)

      !
      ! get number of cells
      !
      CALL nf(nf_inq_dimid(latbc_fileid, 'cell', dimid), routine)
      CALL nf(nf_inq_dimlen(latbc_fileid, dimid, no_cells), routine)

      !
      ! get number of vertical levels
      !
      CALL nf(nf_inq_dimid(latbc_fileid, 'lev', dimid), routine)
      CALL nf(nf_inq_dimlen(latbc_fileid, dimid, no_levels), routine)

      !
      ! check the number of cells and vertical levels
      !
      IF(patch%n_patch_cells_g /= no_cells) THEN
        CALL finish(TRIM(routine),&
        & 'n_patch_cells_g and cells in IFS2ICON file do not match.')
      ENDIF

      WRITE(*,'(a,2i5)') 'nlev_in, no_levels', nlev_in, no_levels
      IF(nlev_in /= no_levels+1) THEN
        CALL finish(TRIM(routine),&
        & 'nlev_in does not match the number of levels in IFS2ICON file.')
      ENDIF

    END IF ! my_process_is_stdio()
      
    !
    ! read prognostic 3d fields
    !
    CALL read_netcdf_data_single( latbc_fileid, 'temp', patch%n_patch_cells_g,          &
      &                           patch%n_patch_cells,  patch%cells%glb_index,          &
      &                           patch%nlev,           p_latbc_data(tlev)%atm%temp     )

    CALL read_netcdf_data_single( latbc_fileid, 'u',    patch%n_patch_cells_g,          &
      &                           patch%n_patch_cells,  patch%cells%glb_index,          &
      &                           patch%nlev,           p_latbc_data(tlev)%atm_in%u     )

    CALL read_netcdf_data_single( latbc_fileid, 'v',    patch%n_patch_cells_g,          &
      &                           patch%n_patch_cells,  patch%cells%glb_index,          &
      &                           patch%nlev,           p_latbc_data(tlev)%atm_in%v     )

    CALL read_netcdf_data_single( latbc_fileid, 'w',    patch%n_patch_cells_g,          &
      &                           patch%n_patch_cells,  patch%cells%glb_index,          &
      &                           patch%nlevp1,         p_latbc_data(tlev)%atm%w        )

    CALL read_netcdf_data_single( latbc_fileid, 'qv',   patch%n_patch_cells_g,          &
      &                           patch%n_patch_cells,  patch%cells%glb_index,          &
      &                           patch%nlev,           p_latbc_data(tlev)%atm%qv       )

    CALL read_netcdf_data_single( latbc_fileid, 'qc',   patch%n_patch_cells_g,          &
      &                           patch%n_patch_cells,  patch%cells%glb_index,          &
      &                           patch%nlev,           p_latbc_data(tlev)%atm%qc       )

    CALL read_netcdf_data_single( latbc_fileid, 'qi',   patch%n_patch_cells_g,          &
      &                           patch%n_patch_cells,  patch%cells%glb_index,          &
      &                           patch%nlev,           p_latbc_data(tlev)%atm%qi       )

    CALL read_netcdf_data_single( latbc_fileid, 'qr',   patch%n_patch_cells_g,          &
      &                           patch%n_patch_cells,  patch%cells%glb_index,          &
      &                           patch%nlev,           p_latbc_data(tlev)%atm%qr       )

    CALL read_netcdf_data_single( latbc_fileid, 'qs',   patch%n_patch_cells_g,          &
      &                           patch%n_patch_cells,  patch%cells%glb_index,          &
      &                           patch%nlev,           p_latbc_data(tlev)%atm%qs       )

    CALL read_netcdf_data_single( latbc_fileid, 'pres', patch%n_patch_cells_g,          &
      &                           patch%n_patch_cells,  patch%cells%glb_index,          &
      &                           patch%nlev,           p_latbc_data(tlev)%atm%pres     )

    !
    ! Convert u and v on cell points to vn at edge points
    !
    CALL interp_uv_2_vn( patch, p_int, p_latbc_data(tlev)%atm_in%u,                     &
      &                  p_latbc_data(tlev)%atm_in%v, p_latbc_data(tlev)%atm%vn )

    !
    ! Compute virtual temperature
    !
    CALL virtual_temp( patch, p_latbc_data(tlev)%atm%temp, p_latbc_data(tlev)%atm%qv,   &
      &                p_latbc_data(tlev)%atm%qc, p_latbc_data(tlev)%atm%qi,            &
      &                p_latbc_data(tlev)%atm%qr, p_latbc_data(tlev)%atm%qs,            &
      &                temp_v )

    !
    ! Compute NH prognostic thermodynamical variables 
    !
    CALL convert_thdvars( patch, p_latbc_data(tlev)%atm%pres, temp_v,                   &
      &                   p_latbc_data(tlev)%atm%rho,                                   &
      &                   p_latbc_data(tlev)%atm%exner,                                 &
      &                   p_latbc_data(tlev)%atm%theta_v )

    CALL sync_patch_array(SYNC_E, patch, p_latbc_data(tlev)%atm%vn)
    CALL sync_patch_array(SYNC_C, patch, p_latbc_data(tlev)%atm%w)
    CALL sync_patch_array(SYNC_C, patch, p_latbc_data(tlev)%atm%theta_v)
    CALL sync_patch_array(SYNC_C, patch, p_latbc_data(tlev)%atm%rho)

    !
    ! close file
    !
    IF (my_process_is_stdio()) THEN
      WRITE(message_text,'(a,a)') 'closing file', TRIM(latbc_full_filename)
      CALL message(TRIM(routine), message_text)
      CALL nf(nf_close(latbc_fileid), routine)
    END IF

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
  !!
  !! @par Revision History
  !! Initial version by S. Brdar, DWD (2013-06-13)
  !!
  SUBROUTINE read_latbc_ifs_data(p_patch, p_nh_state, p_int, ext_data, latbc_datetime, tlev)
    TYPE(t_patch), TARGET,  INTENT(IN)  :: p_patch(:)
    TYPE(t_nh_state),       INTENT(IN)  :: p_nh_state  !< nonhydrostatic state on the global domain
    TYPE(t_int_state),      INTENT(IN)  :: p_int
    TYPE(t_external_data),  INTENT(IN)  :: ext_data    !< external data on the global domain
    TYPE(t_datetime),       INTENT(IN)  :: latbc_datetime
    INTEGER,                INTENT(IN)  :: tlev        !< latbc data slot id used for storing time
                                                       !  level boundary data (values are 1 or 2)

    ! local variables
    TYPE(t_patch), POINTER              :: patch
    INTEGER                             :: mpi_comm, ist, dimid, no_cells, &
                                            latbc_fileid, no_levels, varid
    LOGICAL                             :: l_exist

    CHARACTER(MAX_CHAR_LENGTH), PARAMETER :: routine = "mo_nh_latbc::read_latbc_data"
    CHARACTER(LEN=filename_max)           :: latbc_full_filename
                                            
    patch => p_patch(1)                                        
    IF (nlev_in == 0) THEN
      CALL finish(TRIM(routine), "Number of input levels <nlev_in> not yet initialized.")
    END IF
      
    IF(p_test_run) THEN
      mpi_comm = p_comm_work_test
    ELSE
      mpi_comm = p_comm_work
    ENDIF

    IF (my_process_is_stdio()) THEN
      latbc_full_filename = generate_filename(nroot, patch%level, latbc_datetime)

      WRITE(message_text,'(a,a)') 'reading from', TRIM(latbc_full_filename)
      CALL message(TRIM(routine), message_text)
      INQUIRE (FILE=TRIM(ADJUSTL(latbc_full_filename)), EXIST=l_exist)
      IF (.NOT. l_exist) THEN
        WRITE (message_text,'(a,a)') 'file not found:', TRIM(latbc_full_filename)
        CALL finish(TRIM(routine), message_text)
      ENDIF

      !
      ! open file
      !
      CALL nf(nf_open(TRIM(ADJUSTL(latbc_full_filename)), NF_NOWRITE, latbc_fileid), routine)

      !
      ! get number of cells
      !
      CALL nf(nf_inq_dimid(latbc_fileid, 'ncells', dimid), routine)
      CALL nf(nf_inq_dimlen(latbc_fileid, dimid, no_cells), routine)

      !
      ! get number of vertical levels
      !
      CALL nf(nf_inq_dimid(latbc_fileid, 'lev', dimid), routine)
      CALL nf(nf_inq_dimlen(latbc_fileid, dimid, no_levels), routine)

      !
      ! check the number of cells and vertical levels
      !
      IF(patch%n_patch_cells_g /= no_cells) THEN
        CALL finish(TRIM(routine),&
        & 'n_patch_cells_g and cells in IFS2ICON file do not match.')
      ENDIF

      IF(nlev_in /= no_levels) THEN
        CALL finish(TRIM(routine),&
        & 'nlev_in does not match the number of levels in IFS2ICON file.')
      ENDIF

      !
      ! Check if surface pressure (PS) or its logarithm (LNPS) is provided as input
      !
      IF (nf_inq_varid(latbc_fileid, 'PS', varid) == nf_noerr) THEN
        psvar = 'PS'
      ELSE IF (nf_inq_varid(latbc_fileid, 'LNPS', varid) == nf_noerr) THEN
        psvar = 'LNPS'
      ENDIF

      !
      ! Check if model-level surface Geopotential is provided as GEOSP or GEOP_ML
      !
      IF (nf_inq_varid(latbc_fileid, 'GEOSP', varid) == nf_noerr) THEN
        geop_ml_var = 'GEOSP'
      ELSE IF (nf_inq_varid(latbc_fileid, 'GEOP_ML', varid) == nf_noerr) THEN
        geop_ml_var = 'GEOP_ML'
      ELSE
        CALL finish(TRIM(routine),'Could not find model-level sfc geopotential')
      ENDIF

      !
      ! Check if rain water (QR) is provided as input
      !
      IF (nf_inq_varid(latbc_fileid, 'QR', varid) == nf_noerr) THEN
        lread_qr = .true.
      ELSE
        lread_qr = .false.
        CALL message(TRIM(routine),'Rain water (QR) not available in input data')
      ENDIF

      !
      ! Check if snow water (QS) is provided as input
      !
      IF (nf_inq_varid(latbc_fileid, 'QS', varid) == nf_noerr) THEN
        lread_qs = .true.
      ELSE
        lread_qs = .false.
        CALL message(TRIM(routine),'Snow water (QS) not available in input data')
      ENDIF
      
      !
      ! Check if surface pressure (VN) is provided as input
      !
      IF (nf_inq_varid(latbc_fileid, 'VN', varid) == nf_noerr) THEN
        lread_vn = .TRUE.
      ELSE
        lread_vn = .FALSE.
      ENDIF

    END IF ! my_process_is_stdio()


    CALL p_bcast(lread_qs, p_io, mpi_comm)
    CALL p_bcast(lread_qr, p_io, mpi_comm)
    CALL p_bcast(lread_vn, p_io, mpi_comm)

    
    !
    ! read correct prognostic 3d fields from IFS2ICON file
    !
    CALL read_netcdf_data_single( latbc_fileid, 'T', patch%n_patch_cells_g,           &
                                  patch%n_patch_cells, patch%cells%glb_index,         &
                                  nlev_in, p_latbc_data(tlev)%atm_in%temp )

    IF (lread_vn) THEN
      CALL read_netcdf_data_single( latbc_fileid, 'VN', patch%n_patch_edges_g,        &
        &                     patch%n_patch_edges, patch%edges%glb_index,             &
        &                     nlev_in, p_latbc_data(tlev)%atm_in%vn )
    ELSE
      CALL read_netcdf_data_single( latbc_fileid, 'U', patch%n_patch_cells_g,         &
        &                     patch%n_patch_cells, patch%cells%glb_index,             &
        &                     nlev_in, p_latbc_data(tlev)%atm_in%u )

      CALL read_netcdf_data_single( latbc_fileid, 'V', patch%n_patch_cells_g,         &
        &                     patch%n_patch_cells, patch%cells%glb_index,             &
        &                     nlev_in, p_latbc_data(tlev)%atm_in%v )
    ENDIF

    CALL read_netcdf_data_single( latbc_fileid, 'W', patch%n_patch_cells_g,           &
        &                     patch%n_patch_cells, patch%cells%glb_index,             &
        &                     nlev_in, p_latbc_data(tlev)%atm_in%omega )

    CALL read_netcdf_data_single( latbc_fileid, 'QV', patch%n_patch_cells_g,          &
        &                     patch%n_patch_cells, patch%cells%glb_index,             &
        &                     nlev_in, p_latbc_data(tlev)%atm_in%qv )

    CALL read_netcdf_data_single( latbc_fileid, 'QC', patch%n_patch_cells_g,          &
        &                     patch%n_patch_cells, patch%cells%glb_index,             &
        &                     nlev_in, p_latbc_data(tlev)%atm_in%qc )

    CALL read_netcdf_data_single( latbc_fileid, 'QI', patch%n_patch_cells_g,          &
        &                     patch%n_patch_cells, patch%cells%glb_index,             &
        &                     nlev_in, p_latbc_data(tlev)%atm_in%qi )

    IF (lread_qr) THEN
      CALL read_netcdf_data_single( latbc_fileid, 'QR', patch%n_patch_cells_g,        &
        &                     patch%n_patch_cells, patch%cells%glb_index,             &
        &                     nlev_in, p_latbc_data(tlev)%atm_in%qr )
    ELSE
      p_latbc_data(tlev)%atm_in%qr(:,:,:)=0._wp
    ENDIF

    IF (lread_qs) THEN
      CALL read_netcdf_data_single( latbc_fileid, 'QS', patch%n_patch_cells_g,        &
        &                     patch%n_patch_cells, patch%cells%glb_index,             &
        &                     nlev_in, p_latbc_data(tlev)%atm_in%qs )
    ELSE
      p_latbc_data(tlev)%atm_in%qs(:,:,:)=0._wp
    ENDIF

    CALL read_netcdf_data( latbc_fileid, TRIM(psvar), patch%n_patch_cells_g,          &
        &                     patch%n_patch_cells, patch%cells%glb_index,             &
        &                     p_latbc_data(tlev)%atm_in%psfc )

    CALL read_netcdf_data( latbc_fileid, TRIM(geop_ml_var), patch%n_patch_cells_g,    &
        &                     patch%n_patch_cells, patch%cells%glb_index,             &
        &                     p_latbc_data(tlev)%atm_in%phi_sfc )
      
    !
    ! close file
    !
    IF (my_process_is_stdio()) THEN
      WRITE(message_text,'(a,a)') 'closing file', TRIM(latbc_full_filename)
      CALL message(TRIM(routine), message_text)
      CALL nf(nf_close(latbc_fileid), routine)
    END IF

    !
    ! perform vertical interpolation of horizonally interpolated analysis data
    !
    CALL vert_interp(patch, p_int, p_nh_state%metrics, p_latbc_data(tlev))

  END SUBROUTINE read_latbc_ifs_data
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !>
  !! @par Revision History
  !! Initial version by S. Brdar, DWD (2013-06-13)
  !!
  SUBROUTINE deallocate_latbc_data(patch)
    TYPE(t_patch), INTENT(INOUT) :: patch

    ! local variables
    INTEGER             :: mpi_comm
    INTEGER             :: latbc_fileid, tlev
    INTEGER             :: dimid, no_cells, no_levels
    LOGICAL             :: l_exist, l_all_prog_vars
    INTEGER             :: nlev, nlevp1, nblks_c, nblks_v, nblks_e
    CHARACTER(MAX_CHAR_LENGTH), PARAMETER :: routine = &
      "mo_nh_latbc::deallocate_latbc_data"
     
    WRITE(message_text,'(a,a)') 'deallocating latbc data'
    CALL message(TRIM(routine), message_text)
      
    IF(p_test_run) THEN
      mpi_comm = p_comm_work_test
    ELSE
      mpi_comm = p_comm_work
    ENDIF

    nlev = patch%nlev
    nlevp1 = nlev + 1
    nblks_c = patch%nblks_c
    nblks_v = patch%nblks_v
    nblks_e = patch%nblks_e
    
    ! 
    ! deallocate boundary data memory
    !
    DO tlev = 1, 2
      DEALLOCATE(p_latbc_data(tlev)%atm_in%psfc, &
                 p_latbc_data(tlev)%atm_in%phi_sfc, &
                 p_latbc_data(tlev)%atm_in%pres, &
                 p_latbc_data(tlev)%atm_in%temp, &
                 p_latbc_data(tlev)%atm_in%z3d, &
                 p_latbc_data(tlev)%atm_in%u, &
                 p_latbc_data(tlev)%atm_in%v, &
                 p_latbc_data(tlev)%atm_in%w, &
                 p_latbc_data(tlev)%atm_in%omega, &
                 p_latbc_data(tlev)%atm_in%qv, &
                 p_latbc_data(tlev)%atm_in%qc, &
                 p_latbc_data(tlev)%atm_in%qi, &
                 p_latbc_data(tlev)%atm_in%qr, &
                 p_latbc_data(tlev)%atm_in%qs )

      ! Allocate atmospheric output data
      DEALLOCATE(p_latbc_data(tlev)%atm%vn, &
                 p_latbc_data(tlev)%atm%u, &
                 p_latbc_data(tlev)%atm%v, &
                 p_latbc_data(tlev)%atm%w, &
                 p_latbc_data(tlev)%atm%temp, &
                 p_latbc_data(tlev)%atm%exner, &
                 p_latbc_data(tlev)%atm%pres, &
                 p_latbc_data(tlev)%atm%rho, &
                 p_latbc_data(tlev)%atm%theta_v, &
                 p_latbc_data(tlev)%atm%qv, &
                 p_latbc_data(tlev)%atm%qc, &
                 p_latbc_data(tlev)%atm%qi, &
                 p_latbc_data(tlev)%atm%qr, &
                 p_latbc_data(tlev)%atm%qs )

      DEALLOCATE(p_latbc_data(tlev)%atm_in%vn)
    END DO

  END SUBROUTINE deallocate_latbc_data
  !-------------------------------------------------------------------------


END MODULE mo_nh_latbc
