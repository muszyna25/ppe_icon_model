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
MODULE mo_bc_sst_sic
  
  USE mo_kind,               ONLY: dp, i8
  USE mo_exception,          ONLY: finish, message, message_text
  USE mo_mpi,                ONLY: my_process_is_mpi_workroot
  USE mo_scatter,            ONLY: scatter_time_array
  USE mo_model_domain,       ONLY: t_patch
  USE mo_parallel_config,    ONLY: nproma
  USE mo_physical_constants, ONLY: tf_salt !, tmelt
  USE mo_impl_constants,     ONLY: MAX_CHAR_LENGTH
  USE mo_cdi,                ONLY: streamOpenRead, streamInqVlist, gridInqSize,      &
    &                              vlistInqTaxis, streamInqTimestep, taxisInqVdate,  &
    &                              vlistInqVarGrid, streamClose, streamReadVarslice
  USE mo_util_cdi,           ONLY: cdiGetStringError
  USE mo_bcs_time_interpolation, ONLY: t_time_interpolation_weights

  IMPLICIT NONE

  PRIVATE
  
  REAL(dp), POINTER :: sst(:,:,:) => NULL()
  REAL(dp), POINTER :: sic(:,:,:) => NULL()
  
  CHARACTER(len=*), PARAMETER :: sst_fn = 'bc_sst.nc'
  CHARACTER(len=*), PARAMETER :: sic_fn = 'bc_sic.nc'

  PUBLIC :: read_bc_sst_sic
  PUBLIC :: bc_sst_sic_time_interpolation
  PUBLIC :: get_current_bc_sst_sic_year

  INTEGER(i8), SAVE :: current_year = -1

CONTAINS
  
  SUBROUTINE read_bc_sst_sic(year, p_patch)

    INTEGER(i8),   INTENT(in) :: year
    TYPE(t_patch), INTENT(in) :: p_patch
   
    REAL(dp), POINTER :: zin(:,:) => NULL()
    
    LOGICAL :: lexist
    
    ! global
    IF (.NOT. ASSOCIATED(zin)) ALLOCATE(zin(p_patch%n_patch_cells_g, 0:13))

    IF (my_process_is_mpi_workroot()) THEN
   
      WRITE(message_text,'(a,a,a,i0)') &
           'Read SST from ', sst_fn, ' for ', year
      CALL message('',message_text)
      
      INQUIRE (file=sst_fn, exist=lexist)
      IF (lexist) THEN
        CALL read_sst_sic_data(sst_fn, year, zin)
      ELSE
        WRITE (message_text,*) 'Could not open file ',sst_fn
        CALL message('',message_text)
        CALL finish ('mo_bc_sst_sic:read_bc_sst_sic', 'run terminated.')
      ENDIF

    ENDIF

    ! local
    IF (.NOT. ASSOCIATED(sst)) ALLOCATE (sst(nproma, p_patch%nblks_c, 0:13))
    CALL scatter_time_array(zin, sst, p_patch%cells%decomp_info%glb_index)

    IF (my_process_is_mpi_workroot()) THEN

      WRITE(message_text,'(a,a,a,i0)') &
           'Read sea ice from ', sic_fn, ' for ', year
      CALL message('',message_text)
      
      INQUIRE (file=sic_fn, exist=lexist)
      IF (lexist) THEN
        CALL read_sst_sic_data(sic_fn, year, zin)
      ELSE
        WRITE (message_text,*) 'Could not open file ', sic_fn
        CALL message('',message_text)
        CALL finish ('mo_bc_sst_sic:read_bc_sst_sic', 'run terminated.')
      ENDIF
      
    ENDIF

    ! local
    IF (.NOT. ASSOCIATED(sic)) ALLOCATE (sic(nproma, p_patch%nblks_c, 0:13))
    CALL scatter_time_array(zin, sic, p_patch%cells%decomp_info%glb_index)
    
    IF (ASSOCIATED(zin)) DEALLOCATE(zin)
    
    current_year = year

  END SUBROUTINE read_bc_sst_sic
  
  SUBROUTINE read_sst_sic_data(fn, y, zin)
    
    CHARACTER(len=*), INTENT(in) :: fn
    INTEGER(i8), INTENT(in) :: y
    REAL(dp), POINTER :: zin(:,:)
    
    INTEGER :: ngridsize
    INTEGER(i8) :: ym1, yp1
    
    INTEGER :: taxisID
    INTEGER :: vlistID, varID, streamID, tsID
    INTEGER :: nmiss, status, vdate, vyear, vmonth

    REAL(dp), ALLOCATABLE :: buffer(:)

    CHARACTER(LEN=MAX_CHAR_LENGTH) :: cdiErrorText

    ym1 = y-1
    yp1 = y+1

    streamID = streamOpenRead(fn)
    IF ( streamID < 0 ) THEN
      CALL cdiGetStringError(streamID, cdiErrorText)
      WRITE(message_text,*) TRIM(cdiErrorText)
      CALL finish('mo_bc_sst_sic:read_sst_sic_data', message_text)
    END IF
    
    vlistID = streamInqVlist(streamID)
    varID = 0
    ngridsize = gridInqSize(vlistInqVarGrid(vlistID, varID))

    ALLOCATE(buffer(ngridsize))

    taxisID = vlistInqTaxis(vlistID)
    tsID = 0

    DO
      status = streamInqTimestep(streamID, tsID)
      IF ( status == 0 ) EXIT
      vdate = taxisInqVdate(taxisID)
      vyear = vdate/10000
      vmonth = (vdate/100)-vyear*100
      IF (INT(vyear,i8) == ym1 .AND. vmonth == 12) THEN
        CALL streamReadVarslice(streamID, varID, 0, buffer, nmiss)
        zin(:,0) = buffer(:)
      ELSE IF (INT(vyear,i8) == y) THEN
        CALL streamReadVarslice(streamID, varID, 0, buffer, nmiss)
        zin(:,vmonth) = buffer(:)
      ELSE IF (INT(vyear,i8) == yp1 .AND. vmonth == 1) THEN
        CALL streamReadVarslice(streamID, varID, 0, buffer, nmiss)
        zin(:,13) = buffer(:)
        EXIT
      ENDIF
      tsID = tsID+1
    END DO

    CALL streamClose(streamID)
    
    DEALLOCATE(buffer)

  END SUBROUTINE read_sst_sic_data

  SUBROUTINE bc_sst_sic_time_interpolation(tiw, slf, tsw, seaice, siced, p_patch)
    TYPE( t_time_interpolation_weights), INTENT(in) :: tiw
    REAL(dp)       , INTENT(in)  :: slf(:,:) 
    REAL(dp)       , INTENT(out) :: tsw(:,:) 
    REAL(dp)       , INTENT(out) :: seaice(:,:) 
    REAL(dp)       , INTENT(out) :: siced(:,:) 
    TYPE(t_patch)  , INTENT(in)  :: p_patch

    REAL(dp) :: zts(SIZE(tsw,1),SIZE(tsw,2))
    REAL(dp) :: zic(SIZE(tsw,1),SIZE(tsw,2))

    zts(:,:) = tiw%weight1 * sst(:,:,tiw%month1_index) + tiw%weight2 * sst(:,:,tiw%month2_index)
    zic(:,:) = tiw%weight1 * sic(:,:,tiw%month1_index) + tiw%weight2 * sic(:,:,tiw%month2_index)

    !TODO: missing siced needs to be added

    WHERE (slf(:,:) < 0.5_dp)
      seaice(:,:) = zic(:,:)*0.01_dp               ! assuming input data is in percent
      ! seaice(:,:) = MAX(0.0_dp, MIN(0.99_dp, zic(:,:)))
      seaice(:,:) = MERGE(0.99_dp, seaice(:,:), seaice(:,:) > 0.99_dp)
      ! IF (seaice(:,:) <= 0.01_dp) seaice(:,:) = 0.0_dp
      seaice(:,:) = MERGE(0.0_dp, seaice(:,:), seaice(:,:) <= 0.01_dp)
      ! IF (seaice(:,:) > 0.0_dp) THEN           ! ice
      !   tsw(:,:)=tf_salt              
      ! ELSE                                     ! water
      !   tsw(:,:)=MAX(zts(:,:), tf_salt)
      ! END IF
      tsw(:,:) = MERGE(tf_salt, MAX(zts(:,:), tf_salt), seaice(:,:) > 0.0_dp) 
      WHERE (p_patch%cells%center(:,:)%lat > 0.0_dp)
         siced(:,:) = MERGE(2.0_dp, 0.0_dp, seaice(:,:) > 0.0_dp) 
      ELSEWHERE
         siced(:,:) = MERGE(1.0_dp, 0.0_dp, seaice(:,:) > 0.0_dp)
      ENDWHERE
    ELSEWHERE                                  ! land
      seaice(:,:) = 0.0_dp
      siced(:,:)  = 0.0_dp
      !TODO: check tsw/i/l sequence,dummy setting to some reasonable value for land and ice
      tsw(:,:) = zts(:,:)
    ENDWHERE

    !CALL message('','Interpolated sea surface temperature and sea ice cover.')

  END SUBROUTINE bc_sst_sic_time_interpolation

  FUNCTION get_current_bc_sst_sic_year() RESULT(this_year)
    INTEGER(i8) :: this_year
    this_year = current_year
  END FUNCTION get_current_bc_sst_sic_year

END MODULE mo_bc_sst_sic
