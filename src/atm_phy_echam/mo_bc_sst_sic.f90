!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!>
!! Preliminary read and time interpolation routine for monthly SST and sea ice data
!! 
!! This is  a clone of the respective ECHAM routine
!!
!! U. Schlese, DKRZ,  May 1993, original version
!! U. Schulzweida, MPI, May 1999, netCDF version
!! L. Kornblueh, MPI, November 2001, cleanup for parallel environment
!! U. Schulzweida, MPI, May 2002, blocking (nproma)
!! L. Kornblueh, MPI, February 2013, adapted as temporary reader in ICON using cdi
!!
!! TODO: ctfreez in echam = 271.38, this is 271.45 K
!
#include "icon_contiguous_defines.h"
MODULE mo_bc_sst_sic
  
  USE mo_kind,               ONLY: dp, i8
  USE mo_exception,          ONLY: finish, message, message_text
#ifdef _OPENACC
  USE mo_exception,          ONLY: warning
#endif
  USE mo_mpi,                ONLY: my_process_is_mpi_workroot, p_bcast, &
    &                              process_mpi_root_id, p_comm_work
  USE mo_model_domain,       ONLY: t_patch
  USE mo_grid_config,        ONLY: n_dom
  USE mo_parallel_config,    ONLY: nproma
  USE mo_physical_constants, ONLY: tf_salt !, tmelt
  USE mo_impl_constants,     ONLY: MAX_CHAR_LENGTH, max_dom
  USE mo_cdi,                ONLY: streamOpenRead, streamInqVlist, streamClose, &
    & vlistInqTaxis, streamInqTimestep, taxisInqVdate, streamReadVarSlice
  USE mo_util_cdi,           ONLY: cdiGetStringError
  USE mo_bcs_time_interpolation, ONLY: t_time_interpolation_weights

  IMPLICIT NONE

  PRIVATE

  TYPE t_ext_sea
    REAL(dp), CONTIGUOUS_POINTER :: sst(:,:,:) => NULL()
    REAL(dp), CONTIGUOUS_POINTER :: sic(:,:,:) => NULL()
  END TYPE t_ext_sea

  TYPE(t_ext_sea), TARGET :: ext_sea(max_dom)

  PUBLIC :: read_bc_sst_sic
  PUBLIC :: bc_sst_sic_time_interpolation
  PUBLIC :: get_current_bc_sst_sic_year

  INTEGER(i8), SAVE :: current_year = -1

CONTAINS
  
  SUBROUTINE read_bc_sst_sic(year, p_patch)
    INTEGER(i8),   INTENT(IN) :: year
    TYPE(t_patch), INTENT(IN) :: p_patch
    INTEGER :: jg
    CHARACTER(len=16) :: fn

    jg = p_patch%id

    IF (n_dom > 1) THEN
      WRITE(fn, '(a,i2.2)') 'bc_sst_DOM', jg, '.nc'
    ELSE
      fn = 'bc_sst.nc'
    END IF
    IF (my_process_is_mpi_workroot()) THEN
      WRITE(message_text,'(3a,i0)') 'Read SST from ', TRIM(fn), ' for ', year
      CALL message('',message_text)
    ENDIF
    IF (.NOT.ASSOCIATED(ext_sea(jg)%sst)) THEN
      ALLOCATE (ext_sea(jg)%sst(nproma, p_patch%nblks_c, 0:13))
      !$ACC ENTER DATA PCREATE( ext_sea(jg)%sst)
    ENDIF
    CALL read_sst_sic_data(p_patch, ext_sea(jg)%sst, TRIM(fn), year)

    IF (n_dom > 1) THEN
      WRITE(fn, '(a,i2.2)') 'bc_sic_DOM', jg, '.nc'
    ELSE
      fn = 'bc_sic.nc'
    END IF
    IF (my_process_is_mpi_workroot()) THEN
      WRITE(message_text,'(3a,i0)') 'Read sea ice from ',TRIM(fn),' for ',year
      CALL message('',message_text)
    ENDIF
    IF (.NOT.ASSOCIATED(ext_sea(jg)%sic)) THEN
      ALLOCATE (ext_sea(jg)%sic(nproma, p_patch%nblks_c, 0:13))
      !$ACC ENTER DATA PCREATE( ext_sea(jg)%sic)
    ENDIF
    CALL read_sst_sic_data(p_patch, ext_sea(jg)%sic, TRIM(fn), year)
    
    IF (jg==n_dom) current_year = year

#ifdef _OPENACC
    CALL warning("GPU:read_bc_sst_sic", "GPU device synchronization")
#endif
    !$ACC UPDATE DEVICE( ext_sea(jg)%sst, ext_sea(jg)%sic )

  END SUBROUTINE read_bc_sst_sic
  
  SUBROUTINE read_sst_sic_data(p_patch, dst, fn, y)
!TODO: switch to reading via mo_read_netcdf_distributed?
    TYPE(t_patch), INTENT(in) :: p_patch
    REAL(dp), CONTIGUOUS_ARGUMENT(INOUT) :: dst(:,:,0:)
    CHARACTER(len=*), INTENT(IN) :: fn
    INTEGER(i8), INTENT(in) :: y
    REAL(dp), ALLOCATABLE :: zin(:)
    REAL(dp) :: dummy(0)
    INTEGER :: vlID, taxID, tsID, ts_idx, strID, nmiss, vd, vy, vm, ts_found
    LOGICAL :: found_last_ts, lexist
    CHARACTER(LEN=MAX_CHAR_LENGTH) :: cdiErrorText

    IF (my_process_is_mpi_workroot()) THEN
      INQUIRE (file=fn, exist=lexist)
      IF (.NOT.lexist) THEN
        WRITE (message_text, '(3a)') 'Could not open file ', fn, ': run terminated.'
        CALL finish ('mo_bc_sst_sic:read_sst_sic_data', message_text)
      ENDIF
      strID = streamOpenRead(fn)
      IF ( strID < 0 ) THEN
        CALL cdiGetStringError(strID, cdiErrorText)
        WRITE (message_text, '(4a)') 'Could not open file ', fn, ': ', cdiErrorText
        CALL finish('mo_bc_sst_sic:read_sst_sic_data', message_text)
      END IF
      vlID = streamInqVlist(strID)
      taxID = vlistInqTaxis(vlID)
      tsID = 0
      found_last_ts = .FALSE.
      ts_found = 0
      ALLOCATE(zin(p_patch%n_patch_cells_g))
      DO WHILE (.NOT. found_last_ts)
        IF (streamInqTimestep(strID, tsID) == 0) EXIT
        vd = taxisInqVdate(taxID)
        vy = vd/10000
        vm = (vd/100)-vy*100
        ts_idx = -1
        IF (INT(vy,i8) == y-1_i8 .AND. vm == 12) THEN
          ts_idx = 0
          ts_found = ts_found + 1
        ELSE IF (INT(vy,i8) == y) THEN
          ts_idx = vm
          ts_found = ts_found + 1
        ELSE IF (INT(vy,i8) == y+1_i8 .AND. vm == 1) THEN
          ts_idx = 13
          ts_found = ts_found + 1
          found_last_ts = .TRUE.
        END IF
        IF (ts_idx /= -1) THEN
          CALL streamReadVarSlice(strID, 0, 0, zin, nmiss)
          CALL p_bcast(ts_idx, process_mpi_root_id, p_comm_work)
          dst(:,SIZE(dst,2),ts_idx) = 0._dp
          CALL p_patch%comm_pat_scatter_c%distribute(zin, dst(:,:,ts_idx), .FALSE.)
        ENDIF
        tsID = tsID+1
      END DO
      IF (ts_found < 14) THEN
          CALL finish ('mo_bc_sst_sic:read_sst_sic_data', &
            & 'could not read required data from input file')
      END IF
      ts_idx = -1
      CALL p_bcast(ts_idx, process_mpi_root_id, p_comm_work)
      DEALLOCATE(zin)
      CALL streamClose(strID)
    ELSE
      ts_idx = 0
      DO
        CALL p_bcast(ts_idx, process_mpi_root_id, p_comm_work)
        IF(ts_idx .EQ. -1) EXIT
        dst(:,SIZE(dst,2),ts_idx) = 0._dp
        CALL p_patch%comm_pat_scatter_c%distribute(dummy, dst(:,:,ts_idx), .FALSE.)
      END DO
    END IF
  END SUBROUTINE read_sst_sic_data

  SUBROUTINE bc_sst_sic_time_interpolation(tiw, tsw, seaice, siced, p_patch, mask, l_init, lopenacc)
    
    TYPE( t_time_interpolation_weights), INTENT(in) :: tiw
    REAL(dp)       , INTENT(inout) :: tsw(:,:) 
    REAL(dp)       , INTENT(out) :: seaice(:,:) 
    REAL(dp)       , INTENT(out) :: siced(:,:) 
    TYPE(t_patch)  , INTENT(in)  :: p_patch
    LOGICAL        , INTENT(in)  :: mask(:,:)  !< logical mask, indicating where to apply tsw and sea ice/depth
    LOGICAL        , INTENT(in)  :: l_init     !< switch for first call at initialization
    LOGICAL, OPTIONAL, INTENT(in) :: lopenacc  ! Flag to run on GPU
    ! If l_init=.FALSE., tsw is only computed where mask==.TRUE (this is used
    ! at all time steps in the time loop). If l_init=.TRUE., tsw is computed 
    ! everywhere (this is used during initialization in order to initialize ts_tile(:,:,iwtr)
    ! also over land/lakes). 
    ! Note that lakes and ocean/sea ice are mutually exclusive, i.e. a cell cannot
    ! contain both lake and ocean/sea ice.

    REAL(dp) :: zts(SIZE(tsw,1),SIZE(tsw,2))
    REAL(dp) :: zic(SIZE(tsw,1),SIZE(tsw,2))
    REAL(dp) :: ztsw(SIZE(tsw,1),SIZE(tsw,2))
    REAL(dp), CONTIGUOUS_POINTER :: sst(:,:,:), sic(:,:,:)

    INTEGER  :: jc, jb, jg
#ifdef _OPENACC
    LOGICAL  :: lzopenacc

    IF (PRESENT(lopenacc)) THEN
      lzopenacc = lopenacc
    ELSE
      lzopenacc = .FALSE.
    ENDIF
#endif

    jg = p_patch%id

    sst => ext_sea(jg)%sst
    sic => ext_sea(jg)%sic
    !$ACC DATA PRESENT( tsw, seaice, siced, mask, sst, sic, p_patch%cells%center ) &
    !$ACC       CREATE( zts, zic, ztsw )                                                            &
    !$ACC           IF( lzopenacc )

    !$ACC PARALLEL DEFAULT(PRESENT) IF( lzopenacc )
    !$ACC LOOP SEQ
!$omp parallel private(jb,jc)
!$omp do
    DO jb = 1, SIZE(zts,2)
      !$ACC LOOP GANG VECTOR
      DO jc = 1, SIZE(zts,1)
        zts(jc,jb) = tiw%weight1 * sst(jc,jb,tiw%month1_index) + tiw%weight2 * sst(jc,jb,tiw%month2_index)
        zic(jc,jb) = tiw%weight1 * sic(jc,jb,tiw%month1_index) + tiw%weight2 * sic(jc,jb,tiw%month2_index)
      END DO
    END DO
!$omp end do nowait
    !$ACC END PARALLEL

    !TODO: missing siced needs to be added

    !$ACC PARALLEL DEFAULT(PRESENT) IF( lzopenacc )
    !$ACC LOOP SEQ
!$omp do
    DO jb = 1, SIZE(tsw,2)
      !$ACC LOOP GANG VECTOR
      DO jc = 1, SIZE(tsw,1)
        ! assuming input data is in percent
        seaice(jc,jb) = zic(jc,jb) * MERGE(0.01_dp, 0.0_dp, mask(jc,jb))
      END DO
    END DO
!$omp end do nowait
    !$ACC END PARALLEL
    !$ACC PARALLEL DEFAULT(PRESENT) IF( lzopenacc )
    !$ACC LOOP SEQ
!$omp do
    DO jb = 1, SIZE(tsw,2)
      !$ACC LOOP GANG VECTOR
      DO jc = 1, SIZE(tsw,1)
        seaice(jc,jb) = MERGE(0.99_dp, seaice(jc,jb), seaice(jc,jb) > 0.99_dp)
        seaice(jc,jb) = MERGE(0.0_dp, seaice(jc,jb), seaice(jc,jb) <= 0.01_dp)

        ztsw(jc,jb) = MERGE(tf_salt, MAX(zts(jc,jb), tf_salt), seaice(jc,jb) > 0.0_dp) 
      END DO
    END DO
!$omp end do nowait
    !$ACC END PARALLEL

    IF (l_init) THEN
      !$ACC PARALLEL DEFAULT(PRESENT) IF( lzopenacc )
      !$ACC LOOP SEQ
!$omp do
      DO jb = 1, SIZE(tsw,2)
        !$ACC LOOP GANG VECTOR
        DO jc = 1, SIZE(tsw,1)
          tsw(jc,jb) = ztsw(jc,jb)
        END DO
      END DO
!$omp end do nowait
      !$ACC END PARALLEL
    ELSE
      !$ACC PARALLEL DEFAULT(PRESENT) IF( lzopenacc )
      !$ACC LOOP SEQ
!$omp do
      DO jb = 1, SIZE(tsw,2)
        !$ACC LOOP GANG VECTOR
        DO jc = 1, SIZE(tsw,1)
          IF (mask(jc,jb)) THEN
            tsw(jc,jb) = ztsw(jc,jb)
          END IF
        END DO
      END DO
!$omp end do nowait
      !$ACC END PARALLEL
    END IF

    !$ACC PARALLEL DEFAULT(PRESENT) IF( lzopenacc )
    !$ACC LOOP SEQ
!$omp do
    DO jb = 1, SIZE(tsw,2)
      !$ACC LOOP GANG VECTOR
      DO jc = 1, SIZE(tsw,1)
        IF (seaice(jc,jb) > 0.0_dp) THEN
          siced(jc,jb) = MERGE(2._dp, 1._dp, p_patch%cells%center(jc,jb)%lat > 0.0_dp)
        ELSE
          siced(jc,jb) = 0._dp
        END IF
      END DO
    END DO
!$omp end do nowait
!$omp end parallel
    !$ACC END PARALLEL

    !CALL message('','Interpolated sea surface temperature and sea ice cover.')

    !$ACC END DATA

  END SUBROUTINE bc_sst_sic_time_interpolation

  FUNCTION get_current_bc_sst_sic_year() RESULT(this_year)
    INTEGER(i8) :: this_year
    this_year = current_year
  END FUNCTION get_current_bc_sst_sic_year

END MODULE mo_bc_sst_sic
