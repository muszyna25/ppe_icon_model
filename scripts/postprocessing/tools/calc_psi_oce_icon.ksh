#!/bin/ksh

#  Calculate barotropic stream function from icon ocean model output
#
#  Author: Stephan Lorenz, MPIfMet, 06/2012
#
#  Input : averaged icon-ocean standard output file including variables u_vint and wet_c 
#  Output: interpolated input/output/plot-files named u_vint/psi.r360x180/nclpsi
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
# 1x1 deg resolution
nlon=360
nlat=180

# direct assignment of filenames
#avgfile=avg.5y.xom.bliz.r9368.45
#filestr=avg45.2164ym
#avgfile=timmean_2420-2470
#filestr=2420-2470ym

# via parameter
avgfile=${1}
filestr=${2}
weightfile=${3}

echo "Input file is '$avgfile'"
echo "Tag is '$filestr'"
echo "Weightsfile is '$weightfile'"
### input parameter end

# file names
resol=r${nlon}x${nlat}
inpfile=uvint.${resol}.$filestr
outfile=psi.${resol}.$filestr
plotfile=nclpsi.$filestr

echo " working on files: input: $inpfile ; output: $outfile ; plot: $plotfile ; resolution: $resol"

# select wet_c (lsm) and u_vint (vertical integral of u)
if [ ! -s $inpfile.srv ]; then
  if [ -z "$weightfile" ]; then
    echo "weightfile is not given! Use remapnn ..."
    cdo -P 8 remapnn,$resol -selvar,u_vint,wet_c $avgfile $inpfile.nc
  else
    echo "weightfile is given:'$weightfile'. Use remap ..."
    cdo remap,$resol,$weightfile -selvar,u_vint,wet_c, $avgfile $inpfile.nc
  fi

  # store grid description, convert to service format
  echo "cdo -v griddes $inpfile.nc > griddes.$resol"
  cdo griddes $inpfile.nc > griddes.$resol

  echo "cdo -v -f srv copy $inpfile.nc $inpfile.srv"
  cdo -f srv copy $inpfile.nc $inpfile.srv
  rm $inpfile.nc
else
  echo "Use input for the psi-computation:'$inpfile.srv'"
fi

cat > scr-psiread.f90 <<EOF
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


  do jk=1,nlev
    read (11) isrv
    read (11) wet_c(:,jk,:)

    IF (jk == 1) THEN
      write(80) (isrv(jb),jb=1,8)
      write(80) ((wet_c(jlon,jk,jlat),jlon=1,nlon),jlat=1,nlat)
    END IF
  enddo

  read (11) isrv
  write (*,*) isrv
  read (11) z_uint_reg(:,:)
! write(*,*) 'jx=',jx,' jy=',jy,' read uvint=',z_uint_reg(jx,jy)


  ! (3) calculate meridional integral on regular grid starting from south pole:

  DO jlat = nlat-1, 1, -1
    z_uint_reg(:,jlat) = z_uint_reg(:,jlat) + z_uint_reg(:,jlat+1)
  END DO
  ! DO jlat = 2, nlat
  !   z_uint_reg(:,jlat) = z_uint_reg(:,jlat) + z_uint_reg(:,jlat-1)
  ! END DO
  write(*,*) 'jx=',jx,' jy=',jy,' int. uvint=',z_uint_reg(jx,jy)

  ! (4) calculate stream function: scale with length of meridional resolution:

  erad = 6.371229e6                 !  earth's radius [m]
  pi   = 3.141592653
  z_lat_dist = pi/real(nlat)*erad   !  z_lat_dist = dlat* pi*R/180 ; dlat=180/nlat

  !psi_reg(:,:) = z_uint_reg(:,:) * z_lat_dist * rho_ref * wet_c(:,1,:) * 1.0e-9 ! e+9 [kg/s]
  psi_reg(:,:) = z_uint_reg(:,:) * z_lat_dist * wet_c(:,1,:) * 1.0e-6           ! e+6 [m3/s]


  write(80) (isrv(jb),jb=1,8)
  write(80) ((psi_reg(jlon,jlat),jlon=1,nlon),jlat=1,nlat)


END PROGRAM psiread
!-------------------------------------------------------------------------  
EOF

gfortran -o scr-psiread.x scr-psiread.f90
./scr-psiread.x
rm scr-psiread.*

# convert back to netcdf
cdo -f nc -g griddes.$resol chvar,var4,psi -chvar,var1,wet_c $outfile.srv $outfile.nc
rm $outfile.srv

# plot with nclsh:
nclsh $ICONPLOT \
  -iFile=$outfile.nc -oFile=$plotfile -varName=psi -timeStep=0 -oType=ps \
  -selMode=manual -minVar=-60 -maxVar=-10 -numLevs=25 -bStrg=' ' -maskName=wet_c \
  -withLineLabels \
 -mapLLC=-77,20 -mapURC=-70,35 #north-atlantic gyre
#  -mapLLC=120,10 -mapURC=160,40 # japanise WBC
# -plotLevs=-150,-100,-75,-50,-30,-20,-15,-10,-5,0,5,10,15,20,30,50,75,100,150 -withLineLabels
# -plotLevs=-150,-100,-75,-50,-30,-20,-15,-10,-5,0,5,10,15,20,30,50,75,100,150 -withLineLabels
# -maxView \
# -selMode=manual -minVar=-250 -maxVar=250 -numLevs=20
# -selMode=manual -minVar=-250 -maxVar=200 -numLevs=15
# -selMode=manual -minVar=-300 -maxVar=150 -numLevs=15
# -selMode=manual -minVar=-240 -maxVar=210 -numLevs=15
# -maskName=wet_c -selMode=manual -minVar=-240 -maxVar=210 -numLevs=15 \
# -maskName=wet_c -selMode=manual -minVar=-100 -maxVar=100 -numLevs=20 \

exit

# -selMode=halflog -scaleLimit=1

