SUBROUTINE lcariolle_init(                            &
         & open_file,        close_file,              &
         & read_3d_var,      read_1d_var,             &
         & NCX,              nlev                     )
USE mo_lcariolle_kind,         ONLY: wp,wi
USE mo_lcariolle_types,        ONLY: nlatx,nlevx,nmonthx,pvi,avi
IMPLICIT NONE
INTEGER, EXTERNAL          :: open_file
EXTERNAL close_file, read_3d_var, read_1d_var
INTEGER, INTENT(in)        :: NCX, nlev 
INTEGER                 :: file_id
CHARACTER(LEN=17)            :: fname='cariolle_coeff.nc'
REAL(wp)                     :: a_3d(nmonthx,nlevx,nlatx)
REAL(wp)                     :: deg2rad

file_id=open_file(fname)
CALL read_1d_var(file_id, 'lat',nlatx,pvi%rlat(1:nlatx))
CALL read_1d_var(file_id, 'plev',nlevx,pvi%plev)
pvi%plev=pvi%plev*100._wp
CALL read_3d_var(file_id, 'a1',nmonthx,nlevx,nlatx,a_3d)
pvi%a1(1:nlatx,:,1:12)=RESHAPE(a_3d,(/nlatx,nlevx,nmonthx/),ORDER=(/3,2,1/))
CALL read_3d_var(file_id, 'a2',nmonthx,nlevx,nlatx,a_3d)
pvi%a2(1:nlatx,:,1:12)=RESHAPE(a_3d,(/nlatx,nlevx,nmonthx/),ORDER=(/3,2,1/))
CALL read_3d_var(file_id, 'a3',nmonthx,nlevx,nlatx,a_3d)
pvi%a3(1:nlatx,:,1:12)=RESHAPE(a_3d,(/nlatx,nlevx,nmonthx/),ORDER=(/3,2,1/))
CALL read_3d_var(file_id, 'a4',nmonthx,nlevx,nlatx,a_3d)
pvi%a4(1:nlatx,:,1:12)=RESHAPE(a_3d,(/nlatx,nlevx,nmonthx/),ORDER=(/3,2,1/))
CALL read_3d_var(file_id, 'a5',nmonthx,nlevx,nlatx,a_3d)
pvi%a5(1:nlatx,:,1:12)=RESHAPE(a_3d,(/nlatx,nlevx,nmonthx/),ORDER=(/3,2,1/))
CALL read_3d_var(file_id, 'a6',nmonthx,nlevx,nlatx,a_3d)
pvi%a6(1:nlatx,:,1:12)=RESHAPE(a_3d,(/nlatx,nlevx,nmonthx/),ORDER=(/3,2,1/))
CALL read_3d_var(file_id, 'a7',nmonthx,nlevx,nlatx,a_3d)
pvi%a7(1:nlatx,:,1:12)=RESHAPE(a_3d,(/nlatx,nlevx,nmonthx/),ORDER=(/3,2,1/))
CALL read_3d_var(file_id, 'a8',nmonthx,nlevx,nlatx,a_3d)
pvi%a8(1:nlatx,:,1:12)=RESHAPE(a_3d,(/nlatx,nlevx,nmonthx/),ORDER=(/3,2,1/))
CALL close_file(file_id)
IF (pvi%rlat(1)<0._wp) THEN
  pvi%l_lat_sn=.TRUE.
  pvi%rlat(0)=-90._wp
  pvi%rlat(nlatx+1)=90._wp
ELSE
  pvi%l_lat_sn=.FALSE.
  pvi%rlat(0)=90._wp
  pvi%rlat(nlatx+1)=-90._wp
END IF
deg2rad=acos(-1._wp)/180._wp
pvi%rlat=deg2rad*pvi%rlat
pvi%delta_lat=ABS(pvi%rlat(1)-pvi%rlat(2))
! complete the fields for easier interpolation routines
pvi%a1(1:nlatx,:,0)=pvi%a1(1:nlatx,:,12)
pvi%a2(1:nlatx,:,0)=pvi%a2(1:nlatx,:,12)
pvi%a3(1:nlatx,:,0)=pvi%a3(1:nlatx,:,12)
pvi%a4(1:nlatx,:,0)=pvi%a4(1:nlatx,:,12)
pvi%a5(1:nlatx,:,0)=pvi%a5(1:nlatx,:,12)
pvi%a6(1:nlatx,:,0)=pvi%a6(1:nlatx,:,12)
pvi%a7(1:nlatx,:,0)=pvi%a7(1:nlatx,:,12)
pvi%a8(1:nlatx,:,0)=pvi%a8(1:nlatx,:,12)
pvi%a1(1:nlatx,:,13)=pvi%a1(1:nlatx,:,1)
pvi%a2(1:nlatx,:,13)=pvi%a2(1:nlatx,:,1)
pvi%a3(1:nlatx,:,13)=pvi%a3(1:nlatx,:,1)
pvi%a4(1:nlatx,:,13)=pvi%a4(1:nlatx,:,1)
pvi%a5(1:nlatx,:,13)=pvi%a5(1:nlatx,:,1)
pvi%a6(1:nlatx,:,13)=pvi%a6(1:nlatx,:,1)
pvi%a7(1:nlatx,:,13)=pvi%a7(1:nlatx,:,1)
pvi%a8(1:nlatx,:,13)=pvi%a8(1:nlatx,:,1)
pvi%a1(0,:,:)      =pvi%a1(1,:,:)
pvi%a1(nlatx+1,:,:)=pvi%a1(nlatx,:,:)
pvi%a2(0,:,:)      =pvi%a2(1,:,:)
pvi%a2(nlatx+1,:,:)=pvi%a2(nlatx,:,:)
pvi%a3(0,:,:)      =pvi%a3(1,:,:)
pvi%a3(nlatx+1,:,:)=pvi%a3(nlatx,:,:)
pvi%a4(0,:,:)      =pvi%a4(1,:,:)
pvi%a4(nlatx+1,:,:)=pvi%a4(nlatx,:,:)
pvi%a5(0,:,:)      =pvi%a5(1,:,:)
pvi%a5(nlatx+1,:,:)=pvi%a5(nlatx,:,:)
pvi%a6(0,:,:)      =pvi%a6(1,:,:)
pvi%a6(nlatx+1,:,:)=pvi%a6(nlatx,:,:)
pvi%a7(0,:,:)      =pvi%a7(1,:,:)
pvi%a7(nlatx+1,:,:)=pvi%a7(nlatx,:,:)
pvi%a8(0,:,:)      =pvi%a8(1,:,:)
pvi%a8(nlatx+1,:,:)=pvi%a8(nlatx,:,:)
ALLOCATE(avi%tmprt(NCX,nlev))
ALLOCATE(avi%vmr2molm2(NCX,nlev))
ALLOCATE(avi%pres(NCX,nlev))
ALLOCATE(avi%o3_vmr(NCX,nlev))
ALLOCATE(avi%cell_center_lat(NCX))
END SUBROUTINE lcariolle_init
