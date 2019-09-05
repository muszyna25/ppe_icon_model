MODULE mo_reader_sst_sic

  USE mo_kind,                    ONLY: dp, wp
  USE mo_parallel_config,         ONLY: get_nproma
  USE mo_exception,               ONLY: finish
  USE mo_reader_abstract,         ONLY: t_abstract_reader
  USE mo_util_mtime,              ONLY: t_datetime_ptr, mtime_convert_netcdf_units
  USE mo_io_units,                ONLY: FILENAME_MAX
  USE mo_model_domain,            ONLY: t_patch
  USE mo_read_interface,          ONLY: nf
  USE mtime,                      ONLY: newdatetime, datetime, deallocateDatetime, &
       &                                OPERATOR(*), OPERATOR(+), &
       &                                datetimetostring, max_datetime_str_len, timedelta,  &
       &                                newTimeDelta, deallocateTimedelta
  USE mo_mpi,                     ONLY: my_process_is_stdio, my_process_is_mpi_workroot, &
       &                                process_mpi_root_id, p_comm_work, p_bcast, p_pe_work
  USE mo_read_netcdf_distributed, ONLY: distrib_nf_open, distrib_read, distrib_nf_close, &
       &                                idx_lvl_blk

  IMPLICIT NONE

  PRIVATE

  INCLUDE 'netcdf.inc'

  PUBLIC :: t_sst_sic_reader

  TYPE, EXTENDS(t_abstract_reader) :: t_sst_sic_reader

    TYPE(t_patch), POINTER      :: p_patch => NULL()
    CHARACTER(len=NF_MAX_NAME)  :: varnames(2)
    CHARACTER(len=FILENAME_MAX) :: filename
    INTEGER                     :: fileid, dist_fileid
    LOGICAL                     :: lopened = .FALSE.
    
  CONTAINS

    PROCEDURE :: init            => sst_sic_init_reader
    PROCEDURE :: get_one_timelev => sst_sic_get_one_timelevel
    PROCEDURE :: get_times       => sst_sic_get_times
    PROCEDURE :: deinit          => sst_sic_deinit_reader

    PROCEDURE :: get_nblks       => sst_sic_get_nblks
    PROCEDURE :: get_npromz      => sst_sic_get_npromz

  END TYPE t_sst_sic_reader

  CHARACTER(len=*), PARAMETER :: modname = 'mo_reader_sst_sic'

CONTAINS

  SUBROUTINE sst_sic_init_reader(this, p_patch, filename)
    CLASS(t_sst_sic_reader),    INTENT(inout) :: this
    TYPE(t_patch),      TARGET, INTENT(in   ) :: p_patch
    CHARACTER(len=*),           INTENT(in   ) :: filename

    CHARACTER(len=*), PARAMETER :: routine = 'sst_sic_init_reader'

    this%filename = TRIM(filename)
    this%varnames(1) = "SST"
    this%varnames(2) = "SIC"

    this%p_patch => p_patch
    
    IF (.NOT. this%lopened) THEN
      IF (my_process_is_mpi_workroot()) THEN
        CALL nf(nf_open(this%filename, nf_nowrite, this%fileid), routine)
      ENDIF
      this%dist_fileid = distrib_nf_open(TRIM(this%filename))
      this%lopened = .TRUE.
    ENDIF
  END SUBROUTINE sst_sic_init_reader

  SUBROUTINE sst_sic_get_times (this, times)

    CLASS(t_sst_sic_reader),           INTENT(inout) :: this
    TYPE(t_datetime_ptr), ALLOCATABLE, INTENT(  out) :: times(:)

    INTEGER                    :: tvid, tdid
    CHARACTER(len=NF_MAX_NAME) :: unit_att
    REAL(wp), ALLOCATABLE      :: times_read(:)
    INTEGER                    :: ntimes
    TYPE(datetime),  POINTER   :: startdatetime
    TYPE(timedelta), POINTER   :: timeAxisUnit
    INTEGER                    :: i
    CHARACTER(len=*), PARAMETER :: routine = 'sst_sic_get_times'

    IF (my_process_is_mpi_workroot()) THEN

      CALL nf(nf_inq_varid(this%fileid, "time", tvid), routine)
      CALL nf(nf_inq_dimid(this%fileid, "time", tdid), routine)
      CALL nf(nf_inq_dimlen(this%fileid, tdid, ntimes), routine)

      ALLOCATE(times_read(ntimes))

      CALL nf(nf_get_var_double(this%fileid, tvid, times_read), routine)
      CALL nf(nf_get_att_text(this%fileid,   tvid, "units", unit_att), routine)

    ENDIF

    CALL p_bcast(ntimes, process_mpi_root_id, p_comm_work)
    IF (.NOT. ALLOCATED(times_read)) THEN
      ALLOCATE(times_read(ntimes))
    ENDIF
    CALL p_bcast(times_read, process_mpi_root_id, p_comm_work)
    CALL p_bcast(unit_att,   process_mpi_root_id, p_comm_work)

    CALL mtime_convert_netcdf_units(unit_att, startdatetime, timeAxisUnit)

    ALLOCATE(times(ntimes))
    DO i = 1, ntimes
      times(i)%ptr => newdatetime("0000-01-01T00:00:00.000")
      IF (REAL(INT(times_read(i)),wp) /= times_read(i)) THEN
        CALL finish(routine, "Timestamp cannot be cast to INT safely. Check your input data.")
      ENDIF
      times(i)%ptr =  startdatetime + INT(times_read(i)) * timeAxisUnit
    ENDDO

    CALL deallocateDatetime(startdatetime)
    CALL deallocateTimedelta(timeAxisUnit)

  END SUBROUTINE sst_sic_get_times

  SUBROUTINE sst_sic_get_one_timelevel(this, timelevel, varname, dat)
    CLASS(t_sst_sic_reader),   INTENT(inout) :: this
    INTEGER,               INTENT(in   ) :: timelevel
    CHARACTER(len=*),      INTENT(in   ) :: varname
    REAL(dp), ALLOCATABLE, INTENT(  out) :: dat(:,:,:,:)

    ALLOCATE(dat(get_nproma(), 1, this%p_patch%nblks_c, 1))
    dat(:,:,:,:) = -1.0_dp
    IF (.NOT. this%lopened) THEN
      CALL finish(modname, '6 hourly SST/Seaice file not open!') 
    END IF
    CALL distrib_read(this%dist_fileid, varname, dat(:,:,:,1), &
         &            1, idx_lvl_blk, this%p_patch%cells%dist_io_data, timelevel, timelevel)
    
    CALL sst_sic_replace_missval(this, dat, -1.0_wp)
  
  END SUBROUTINE sst_sic_get_one_timelevel

  SUBROUTINE sst_sic_replace_missval (this, dat, new_missval)
    CLASS(t_sst_sic_reader), INTENT(inout) :: this
    REAL(wp),                INTENT(inout) :: dat(:,:,:,:)
    REAL(wp),                INTENT(in   ) :: new_missval

    WHERE (dat < -1e10_wp)
      dat = new_missval
    END WHERE 
  END SUBROUTINE sst_sic_replace_missval

  FUNCTION sst_sic_get_nblks (this) RESULT(nblks)
    CLASS(t_sst_sic_reader), INTENT(in   ) :: this
    INTEGER                                :: nblks
    nblks = this%p_patch%nblks_c
  END FUNCTION sst_sic_get_nblks

  FUNCTION sst_sic_get_npromz (this) RESULT(npromz)
    CLASS(t_sst_sic_reader), INTENT(in   ) :: this
    INTEGER                                :: npromz
    npromz = this%p_patch%npromz_c
  END FUNCTION sst_sic_get_npromz

  SUBROUTINE sst_sic_deinit_reader(this)
    CLASS(t_sst_sic_reader), INTENT(inout) :: this
    CHARACTER(len=*), PARAMETER :: routine = 'sst_sic_deinit_reader'
    IF (ASSOCIATED(this%p_patch)) NULLIFY(this%p_patch)
    IF (this%lopened) THEN
      IF (my_process_is_mpi_workroot()) THEN
        CALL nf(nf_close(this%fileid), routine)
      END IF
      CALL distrib_nf_close(this%dist_fileid)
    END IF
  END SUBROUTINE sst_sic_deinit_reader

END MODULE mo_reader_sst_sic
