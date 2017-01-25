!>
!! @brief This subroutine initializes the data structures for the
!! Cariolle interactive ozone scheme Cariolle et al.:
!! Atmos. Chem. Phys. 7, 2183 (2007), and reads all necessary coefficients.
!! documentation: cr2016_10_22_rjs
!!
!! @author Sebastian Rast, MPI-M
!!
!! @par Revision History
!!  Original version Sebastian Rast (2016)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
SUBROUTINE lcariolle_init(                            &
         & open_file,        close_file,              &
         & read_3d_var,      read_1d_var,             &
         & get_constants,    NCX,                     &
         & nlev                                       )
USE mo_lcariolle_kind,         ONLY: wp,wi
USE mo_lcariolle_types,        ONLY: &
     & nlatx,nlevx,nmonthx, & !< number of latitudes, levels, months in climatology
     & pvi,avi                !< derived types: pvi (inside Cariolle),
                              !< avi (variables passed from host model to this submodel)
IMPLICIT NONE
INTEGER, EXTERNAL            :: open_file
EXTERNAL close_file, read_3d_var, read_1d_var !< reading from netcdf files
REAL(wp),EXTERNAL            :: get_constants !< defines physical constants
INTEGER, INTENT(in)          :: NCX,        & !< number of columns as in calling subprograms
                              & nlev          !< number of levels
INTEGER                      :: file_id
CHARACTER(LEN=17)            :: fname='cariolle_coeff.nc' !< standard input file name for Cariolle
                                                          !< coefficients
REAL(wp)                     :: a_3d(nmonthx,nlevx,nlatx) !< temporary read array
REAL(wp)                     :: deg2rad                   !< conversion factor degrees to radiant

! read latitudes, pressure levels and coefficients A_1,...,A_8 for Cariolle scheme
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
! decide whether latitudes are from S->N and convert to radiant
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
! complete the fields wrt to months for easier interpolation routines
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
! complete the fields wrt to latitude for easier interpolation routines
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
! get physical constants
pvi%avogadro=get_constants('avogadro')
! allocate fields for variables passed from host model to submodel
ALLOCATE(avi%tmprt(NCX,nlev))
ALLOCATE(avi%vmr2molm2(NCX,nlev))
ALLOCATE(avi%pres(NCX,nlev))
ALLOCATE(avi%o3_vmr(NCX,nlev))
ALLOCATE(avi%cell_center_lat(NCX))
ALLOCATE(avi%lday(NCX))
END SUBROUTINE lcariolle_init
