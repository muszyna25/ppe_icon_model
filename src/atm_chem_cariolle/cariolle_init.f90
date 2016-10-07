SUBROUTINE cariolle_init(open_file,close_file,read_3d_var,read_1d_var)
!USE mo_cariolle_kind,         ONLY: wp,wi
USE mo_kind,                  ONLY: wp
USE mo_cariolle_types,        ONLY: nlatx,nlevx,nmonthx,pvi
IMPLICIT NONE
INTEGER, EXTERNAL        :: open_file
EXTERNAL close_file, read_3d_var, read_1d_var
INTEGER                 :: file_id
CHARACTER(LEN=17)            :: fname='cariolle_coeff.nc'
REAL(wp)                     :: a_3d(nmonthx,nlevx,nlatx)

file_id=open_file(fname)
write(*,*) 'file_id=',file_id
CALL read_1d_var(file_id, 'lat',nlatx,pvi%rlat)
CALL read_1d_var(file_id, 'plev',nlevx,pvi%plev)
pvi%plev=pvi%plev*100._wp
write(*,*) 'lat=',pvi%rlat
write(*,*) 'plev=',pvi%plev
write(*,*) 'SIZE(a)=',SIZE(a_3d,1),SIZE(a_3d,2),SIZE(a_3d,3)
CALL read_3d_var(file_id, 'a1',nmonthx,nlevx,nlatx,a_3d)
pvi%a1=RESHAPE(a_3d,(/nlatx,nlevx,nmonthx/),ORDER=(/3,2,1/))
 write(*,*) 'SIZE(pvi%a1)=',SIZE(pvi%a1,1),SIZE(pvi%a1,2),SIZE(pvi%a1,3),pvi%a1(1,2,3)
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

END SUBROUTINE cariolle_init
