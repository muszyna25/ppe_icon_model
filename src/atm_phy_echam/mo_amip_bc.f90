!>
!! Preliminary read and time interpolation routine for AMIP SST and sea ice data
!! 
!! This is  a clone of the respective ECHAM routine
!!
!! U. Schlese, DKRZ,  May 1993, original version
!! U. Schulzweida, MPI, May 1999, netCDF version
!! L. Kornblueh, MPI, November 2001, cleanup for parallel environment
!! U. Schulzweida, MPI, May 2002, blocking (nproma)
!! L. Kornblueh, MPI, February 2013, adapted as temporary reader in ICON using cdi
!!
MODULE mo_amip_bc
  
  USE mo_kind,               ONLY: dp
  USE mo_exception,          ONLY: finish, message, message_text
  USE mo_mpi,                ONLY: my_process_is_io, p_io, p_comm_work, p_bcast
  USE mo_model_domain,       ONLY: t_patch
  USE mo_communication,      ONLY: idx_no, blk_no
  USE mo_parallel_config,    ONLY: nproma
  USE mo_datetime,           ONLY: t_datetime, add_time, date_to_time, idaylen, rdaylen
  USE mo_run_config,         ONLY: dtime
  !TODO: ctfreez in echam = 271.38, this is 271.45 K
  USE mo_physical_constants, ONLY: tmelt, tf_salt 

  IMPLICIT NONE

  PRIVATE

  INCLUDE 'cdi.inc'
  
  REAL(dp), ALLOCATABLE :: sst(:,:,:)
  REAL(dp), ALLOCATABLE :: sic(:,:,:)
  
  CHARACTER(len=*), PARAMETER :: sst_fn = 'bc_sst.nc'
  CHARACTER(len=*), PARAMETER :: sic_fn = 'bc_sic.nc'

  ! weighting factors and indices for time interpolation

  REAL(dp):: wgt1, wgt2
  INTEGER :: nmw1, nmw2

  REAL(dp):: wgtd1, wgtd2
  INTEGER :: ndw1, ndw2

  PUBLIC :: read_amip_bc
  PUBLIC :: amip_time_weights
  PUBLIC :: amip_time_interpolation

CONTAINS
  
  SUBROUTINE read_amip_bc(year, p_patch)

    INTEGER,       INTENT(in) :: year
    TYPE(t_patch), INTENT(in) :: p_patch
   
    REAL(dp), POINTER :: zin(:,:) => NULL()
    
    LOGICAL       :: lexist
    
    INTEGER :: j, k, jl, jb
    
    IF (my_process_is_io()) THEN
    
      WRITE(message_text,'(a,a,a,i0)') &
           'Read AMIP SST from ', sst_fn, ' for ', year
      CALL message('',message_text)
      
      INQUIRE (file=sst_fn, exist=lexist)
      IF (lexist) THEN
        CALL read_amip_data(sst_fn, year, zin)
      ELSE
        WRITE (message_text,*) 'Could not open file ',sst_fn
        CALL message('',message_text)
        CALL finish ('mo_amip_bc:read_amip_bc', 'run terminated.')
      ENDIF

    ENDIF

    ! local
    IF (.NOT. ALLOCATED(sst)) ALLOCATE (sst(nproma, p_patch%nblks_c, 0:13))      
    ! global
    IF (.NOT. ASSOCIATED(zin)) ALLOCATE(zin(p_patch%n_patch_cells_g, 0:13))
    
    CALL p_bcast(zin(:,:), p_io, p_comm_work)
    DO k = 0, 13
      DO j = 1, SIZE(p_patch%cells%glb_index)
        jb = blk_no(j)
        jl = idx_no(j)
        sst(jl,jb,k) = zin(p_patch%cells%glb_index(j),k)
      ENDDO
    ENDDO

    IF (my_process_is_io()) THEN

      WRITE(message_text,'(a,a,a,i0)') &
           'Read AMIP sea ice from ', sic_fn, ' for ', year
      CALL message('',message_text)
      
      INQUIRE (file=sic_fn, exist=lexist)
      IF (lexist) THEN
        CALL read_amip_data(sic_fn, year, zin)
      ELSE
        WRITE (message_text,*) 'Could not open file ', sic_fn
        CALL message('',message_text)
        CALL finish ('mo_amip_bc:read_amip_bc', 'run terminated.')
      ENDIF
      
    ENDIF

    ! local
    IF (.NOT. ALLOCATED(sic)) ALLOCATE (sic(nproma, p_patch%nblks_c, 0:13))      
    ! global part done before

    CALL p_bcast(zin(:,:), p_io, p_comm_work)
    DO k = 0, 13
      DO j = 1, SIZE(p_patch%cells%glb_index)
        jb = blk_no(j)
        jl = idx_no(j)
        sic(jl,jb,k) = zin(p_patch%cells%glb_index(j),k)
      ENDDO
    ENDDO
    
    IF (ASSOCIATED(zin)) DEALLOCATE(zin)
    
  END SUBROUTINE read_amip_bc
  
  SUBROUTINE read_amip_data(fn, y, zin)
    
    CHARACTER(len=*), INTENT(in) :: fn
    INTEGER, INTENT(in) :: y
    REAL(dp), POINTER, INTENT(inout) :: zin(:,:)
    
    INTEGER :: ngridsize
    INTEGER :: ym1, yp1
    
    INTEGER :: taxisID
    INTEGER :: vlistID, varID, streamID, tsID
    INTEGER :: nmiss, status, vdate, vyear, vmonth
    
    ym1 = y-1
    yp1 = y+1

    streamID = streamOpenRead(fn)
    IF ( streamID < 0 ) THEN
      WRITE(message_text,*) cdiStringError(streamID)
      CALL finish('mo_amip_bc:read_amip_data', message_text)
    END IF
    
    vlistID = streamInqVlist(streamID)
    varID = 0
    ngridsize = gridInqSize(vlistInqVarGrid(vlistID, varID))
    ALLOCATE(zin(ngridsize,0:13))
    
    taxisID = vlistInqTaxis(vlistID)
    tsID = 0
    
    DO
      status = streamInqTimestep(streamID, tsID)
      IF ( status == 0 ) EXIT
      vdate = taxisInqVdate(taxisID)
      vyear = vdate/10000
      vmonth = (vdate/100)-vyear*100
      IF (vyear == ym1 .AND. vmonth == 12) THEN
        CALL streamReadVar(streamID, varID, zin(:,0), nmiss)
      ELSE IF (vyear == y) THEN
        CALL streamReadVar(streamID, varID, zin(:,vmonth), nmiss)
      ELSE IF (vyear == yp1 .AND. vmonth == 1) THEN
        CALL streamReadVar(streamID, varID, zin(:,13), nmiss)
        EXIT
      ENDIF
      tsID = tsID+1
    END DO
    
    CALL streamClose(streamID)
    
  END SUBROUTINE read_amip_data

  SUBROUTINE amip_time_weights(current_date)

    TYPE(t_datetime), INTENT(in) :: current_date 

    ! calculates weighting factores for AMIP sst and sea ice

    TYPE(t_datetime) :: next_date

    TYPE(t_datetime) :: date_monm1, date_mon, date_monp1
    INTEGER   :: yr, mo, dy, hr, mn, se
    INTEGER   :: isec
    INTEGER   :: imp1, imm1, imlenm1, imlen, imlenp1
    REAL (dp) :: zsec, zdayl
    REAL (dp) :: zmohlf, zmohlfp1, zmohlfm1
    REAL (dp) :: zdh, zdhp1, zdhm1
    
    ! time of next timestep and split
    !
    next_date = current_date
    CALL add_time(dtime,0,0,0,next_date)
    CALL date_to_time(next_date)    

    yr = next_date%year
    mo = next_date%month 
    dy = next_date%day   
    hr = next_date%hour  
    mn = next_date%minute
    se = INT(next_date%second)

    ! month index for AMIP data  (0..13)
    imp1 = mo+1
    imm1 = mo-1
      
    ! determine length of months and position within current month

    date_monm1%calendar = next_date%calendar
    IF (imm1 ==  0) THEN
      date_monm1%year = yr-1;  date_monm1%month = 12;   date_monm1%day = 1;
    ELSE
      date_monm1%year = yr;    date_monm1%month = imm1; date_monm1%day = 1;
    ENDIF
    date_monm1%hour = 0;   date_monm1%minute = 0; date_monm1%second   = 0;
    CALL date_to_time(date_monm1)

    date_monp1%calendar = next_date%calendar
    IF (imp1 == 13) THEN
      date_monp1%year = yr+1;  date_monp1%month = 1;    date_monp1%day = 1;
    ELSE
      date_monp1%year = yr;    date_monp1%month = imp1; date_monp1%day = 1;
    ENDIF
    date_monp1%hour     = 0;   date_monp1%minute   = 0;    date_monp1%second   = 0;
    CALL date_to_time(date_monp1)

    imlenm1 = date_monm1%monlen
    imlen   = date_mon%monlen
    imlenp1 = date_monp1%monlen
      
    zdayl    = rdaylen
    zmohlfm1 = imlenm1*zdayl*0.5_dp
    zmohlf   = imlen  *zdayl*0.5_dp
    zmohlfp1 = imlenp1*zdayl*0.5_dp
    
    ! weighting factors for first/second half of month
      
    nmw1   = mo
      
    ! seconds in the present month
    isec = (dy-1) * idaylen + INT(next_date%second)
    zsec = REAL(isec,dp)
      
    IF(zsec <= zmohlf) THEN                     ! first part of month
      wgt1   = (zmohlfm1+zsec)/(zmohlfm1+zmohlf)
      wgt2   = 1.0_dp-wgt1
      nmw2   = imm1
    ELSE                                        ! second part of month
      wgt2   = (zsec-zmohlf)/(zmohlf+zmohlfp1)
      wgt1   = 1.0_dp-wgt2
      nmw2   = imp1
    ENDIF
      
    ! weighting factors for first/second half of day
      
    ndw1   = 2
      
    zsec = REAL(next_date%second, dp)
    zdh   = 12.0_dp*3600.0_dp
    zdhm1 = zdh
    zdhp1 = zdh
    IF( zsec <= zdh ) THEN                     ! first part of day
      wgtd1  = (zdhm1+zsec)/(zdhm1+zdh)
      wgtd2  = 1.0_dp-wgtd1
      ndw2   = 1
    ELSE                                       ! second part of day
      wgtd2  = (zsec-zdh)/(zdh+zdhp1)
      wgtd1  = 1.0_dp-wgtd2
      ndw2   = 3
    ENDIF
    
  END SUBROUTINE amip_time_weights

  SUBROUTINE amip_time_interpolation(seaice, tsw, slf)
    REAL(dp), INTENT(out) :: seaice(:,:) 
    REAL(dp), INTENT(out) :: tsw(:,:) 
    REAL(dp), INTENT(in) :: slf(:,:) 

    REAL(dp) :: zts(SIZE(tsw,1),SIZE(tsw,2))
    REAL(dp) :: zic(SIZE(tsw,1),SIZE(tsw,2))

    zts(:,:) = wgt1 * sst(:,:,nmw1) + wgt2 * sst(:,:,nmw2)
    zic(:,:) = wgt1 * sic(:,:,nmw1) + wgt2 * sic(:,:,nmw2)

    !TODO: missing siced needs to be added

    WHERE (slf(:,:) < 1.0_dp)
      zic(:,:) = zic(:,:)*0.01_dp               ! assuming input data is in percent
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
    ELSEWHERE                                  ! land
      seaice(:,:) = 0.0_dp
      !TODO check tsw/i/l sequence,dummy setting to some reasonable value for land and ice
      tsw(:,:) = tmelt
    ENDWHERE

  END SUBROUTINE amip_time_interpolation

END MODULE mo_amip_bc
