#!/bin/ksh

#  Calculate barotropic stream function from icon ocean model output
#
#  Author: Stephan Lorenz, MPIfMet, 06/2012
#
#  Method:
#   - vertical interpolation is done in the model, variable u_vint
#   - interpolation done by cdo remapnn, nearest neighbor interpolation
#   - convert to SERVICE format, calculate meridional integral incl. using
#     land-sea-mask, fortran program psiread
#   - convert badk to NETCDF format using grid description file
#   - simple plot using NCL with shell-script nclsh
#          

set -e

### input parameter
#nlon=72
#nlat=36
# 1x1 deg resolution
nlon=360
nlat=180
resol=r${nlon}x${nlat}
avgfile=avg.5y.xom.bliz.r9368.56
inpfile=uvint.${resol}.avg56
outfile=psi.${resol}.avg56
plotfile=nclpsi.avg56
### input parameter end

# select wet_c (lsm) and u_vint (vertical integral of u)
if [ ! -s $inpfile.nc ]; then
  cdo remapnn,$resol -selvar,u_vint,wet_c $avgfile.nc $inpfile.nc
fi

# store grid description, convert to service format
cdo griddes $inpfile.nc > griddes.$resol
cdo -f srv copy $inpfile.nc $inpfile.srv

cat >psiread.f90 <<EOF
!-------------------------------------------------------------------------  
!
!
!!  Calculation of horizontal stream function
!
!>
!!
!! @par Revision History
!! Developed  by  Stephan Lorenz, MPI-M (2012).
!!  based on code from MPIOM
!
! TODO: implement variable output dimension (1 deg resolution) and smoothing extent
!! 
PROGRAM psiread

IMPLICIT NONE

  INTEGER, PARAMETER ::  rho_ref = 1025.022            ! reference density

! INTEGER, PARAMETER ::  nlat = 36                     ! meridional dimension of regular grid
! INTEGER, PARAMETER ::  nlon = 72                     ! zonal dimension of regular grid
  INTEGER, PARAMETER ::  nlat = $nlat                  ! meridional dimension of regular grid
  INTEGER, PARAMETER ::  nlon = $nlon                  ! zonal dimension of regular grid

  INTEGER, PARAMETER ::  nlev = 20                     ! vertical dimension of experiment

  ! smoothing area is 2*jsmth-1 lat/lon areas of 1 deg
  INTEGER, PARAMETER ::  jsmth = 3                  
  INTEGER            :: jb, jc, jk, i_startidx, i_endidx
  INTEGER            :: jlat, jlon, jx, jy
  INTEGER            :: isrv(8)


  REAL               :: z_lat_dist, erad, pi
  REAL               :: z_uint_reg(nlon,nlat)     ! vertical integral on regular grid
  REAL               :: psi_reg(nlon,nlat)        ! horizontal stream function
  REAL               :: wet_c(nlon,nlev,nlat)     ! slm

  !CHARACTER(len=max_char_length), PARAMETER :: routine = ('mo_oce_diagnostics:calc_psi')

  !-----------------------------------------------------------------------

  psi_reg(:,:)    = 0.0
  z_uint_reg(:,:) = 0.0

  ! test calculation - Pacific; note that first row is at Antarctica
  ! latitude  ~  50S = 40 north of SP
  ! longitude ~ 270E = 90W
  jy = 1 + 40*nlat/180
  jx = 1 + 270*nlat/360


  ! (1) barotropic system - done in ICON ocean model:
  !     vertical integration of zonal velocity times vertical layer thickness [m/s*m]
 ! u_vint(:,:)     = 0.0_wp
 ! DO jb = all_cells%start_block, all_cells%end_block
 !   CALL get_index_range(all_cells, jb, i_startidx, i_endidx)
 !   DO jk = 1, n_zlev
 !     DO jc = i_startidx, i_endidx
 !       delta_z = v_base%del_zlev_m(jk)
 !       IF (jk == 1) delta_z = v_base%del_zlev_m(jk) + h(jc,jb)
 !       u_vint(jc,jb) = u_vint(jc,jb) - u(jc,jk,jb)*delta_z*v_base%wet_c(jc,jk,jb)
 !     END DO
 !   END DO
 ! END DO

  ! (2) read barotropic system: first wet_c
  open (11,file="$inpfile.srv", form='unformatted')
  open (80,file="$outfile.srv", form='unformatted')

  !open (11,file='uvint.r72x36.avg45.srv', form='unformatted')
  !open (80,file=  'psi.r72x36.avg45.srv', form='unformatted')

  do jk=1,nlev
    read (11) isrv
    write (*,*) isrv
    read (11) wet_c(:,jk,:)
  enddo

  read (11) isrv
  write (*,*) isrv
  read (11) z_uint_reg(:,:)
  write(*,*) 'jx=',jx,' jy=',jy,' read uvint=',z_uint_reg(jx,jy)


  ! (3) calculate meridional integral on regular grid starting from south pole:

  DO jlat = nlat-1, 1, -1
    z_uint_reg(:,jlat) = z_uint_reg(:,jlat) + z_uint_reg(:,jlat+1)
  END DO
  write(*,*) 'jx=',jx,' jy=',jy,' int. uvint=',z_uint_reg(jx,jy)

  ! (4) calculate stream function: scale with length of meridional resolution:

  erad = 6.371229e6                 !  earth's radius [m]
  pi   = 3.141592653
  z_lat_dist = pi/real(nlat)*erad   !  z_lat_dist = dlat* pi*R/180 ; dlat=180/nlat

  !psi_reg(:,:) = z_uint_reg(:,:) * z_lat_dist * rho_ref * wet_c(:,1,:) * 1.0e-9 ! e+9 [kg/s]
  psi_reg(:,:) = z_uint_reg(:,:) * z_lat_dist * wet_c(:,1,:) * 1.0e-6           ! e+6 [m3/s]

  write(*,*) 'jx=',jx,' jy=',jy,' wet_c     =',wet_c(jx,1,jy)
  write(*,*) 'jx=',jx,' jy=',jy,' psi_reg   =',psi_reg(jx,jy)
  write(*,*) 'write global PSI at idate:', isrv(3)

  write(80) (isrv(jb),jb=1,8)
  write(80) ((psi_reg(jlon,jlat),jlon=1,nlon),jlat=1,nlat)

! write(82,*) (isrv(jb),jb=1,8)
! do jlat=1,nlat
!     write(82,*) 'jlat=',jlat
!     write(82,'(1p10e12.3)') (psi_reg(jlon,jlat),jlon=1,nlon)
!      write(82,'(1p10e12.3)') (wet_c(jlon,1,jlat),jlon=1,nlon)
! enddo
! write(83,*) (isrv(jb),jb=1,8), 'WET(1)'
! write(83,'(72i1)') ((int(wet_c(jlon,1,jlat)+0.1),jlon=1,nlon),jlat=1,nlat)
! write(83,*) (isrv(jb),jb=1,8), 'WET(16)'
! write(83,'(72i1)') ((int(wet_c(jlon,16,jlat)+0.1),jlon=1,nlon),jlat=1,nlat)

END PROGRAM psiread
!-------------------------------------------------------------------------  
EOF

gfortran -o psiread.x psiread.f90 
./psiread.x

# convert back to netcdf
cdo -f nc -g griddes.$resol chvar,var4,psi $outfile.srv $outfile.nc

# plot with nclsh:
nclsh /pool/data/ICON/tools/icon_plot.ncl -altLibDir=/pool/data/ICON/tools -iFile=$outfile.nc -oFile=$plotfile -varName=psi -timeStep=0 -oType=png

exit


