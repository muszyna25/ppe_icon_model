SUBROUTINE cariolle_init( open_file,        close_file,              &
                        & read_3d_var,      read_1d_var,             &
                        & NCX,              nlev                     )
!USE mo_cariolle_kind,         ONLY: wp,wi
USE mo_kind,                  ONLY: wp
USE mo_cariolle_types,        ONLY: nlatx,nlevx,nmonthx,pvi,avi
IMPLICIT NONE
INTEGER, EXTERNAL          :: open_file
EXTERNAL close_file, read_3d_var, read_1d_var
INTEGER, INTENT(in)        :: NCX, nlev 
INTEGER                 :: file_id
CHARACTER(LEN=17)            :: fname='cariolle_coeff.nc'
REAL(wp)                     :: a_3d(nmonthx,nlevx,nlatx)
REAL(wp)                     :: deg2rad

file_id=open_file(fname)
CALL read_1d_var(file_id, 'lat',nlatx,pvi%rlat)
CALL read_1d_var(file_id, 'plev',nlevx,pvi%plev)
pvi%plev=pvi%plev*100._wp
CALL read_3d_var(file_id, 'a1',nmonthx,nlevx,nlatx,a_3d)
pvi%a1=RESHAPE(a_3d,(/nlatx,nlevx,nmonthx/),ORDER=(/3,2,1/))
CALL read_3d_var(file_id, 'a2',nmonthx,nlevx,nlatx,a_3d)
pvi%a2=RESHAPE(a_3d,(/nlatx,nlevx,nmonthx/),ORDER=(/3,2,1/))
CALL read_3d_var(file_id, 'a3',nmonthx,nlevx,nlatx,a_3d)
pvi%a3=RESHAPE(a_3d,(/nlatx,nlevx,nmonthx/),ORDER=(/3,2,1/))
CALL read_3d_var(file_id, 'a4',nmonthx,nlevx,nlatx,a_3d)
pvi%a4=RESHAPE(a_3d,(/nlatx,nlevx,nmonthx/),ORDER=(/3,2,1/))
CALL read_3d_var(file_id, 'a5',nmonthx,nlevx,nlatx,a_3d)
pvi%a5=RESHAPE(a_3d,(/nlatx,nlevx,nmonthx/),ORDER=(/3,2,1/))
CALL read_3d_var(file_id, 'a6',nmonthx,nlevx,nlatx,a_3d)
pvi%a6=RESHAPE(a_3d,(/nlatx,nlevx,nmonthx/),ORDER=(/3,2,1/))
CALL read_3d_var(file_id, 'a7',nmonthx,nlevx,nlatx,a_3d)
pvi%a7=RESHAPE(a_3d,(/nlatx,nlevx,nmonthx/),ORDER=(/3,2,1/))
CALL read_3d_var(file_id, 'a8',nmonthx,nlevx,nlatx,a_3d)
pvi%a8=RESHAPE(a_3d,(/nlatx,nlevx,nmonthx/),ORDER=(/3,2,1/))
CALL close_file(file_id)
deg2rad=acos(-1._wp)/180._wp
pvi%rlat=deg2rad*pvi%rlat
pvi%lat_shift=ABS(pvi%rlat(1))
pvi%delta_lat=ABS(pvi%rlat(1)-pvi%rlat(2))
ALLOCATE(avi%tmprt(NCX,nlev))
ALLOCATE(avi%vmr2molm2(NCX,nlev))
ALLOCATE(avi%o3_vmr(NCX,nlev))
ALLOCATE(avi%cell_center_lat(NCX))
END SUBROUTINE cariolle_init
