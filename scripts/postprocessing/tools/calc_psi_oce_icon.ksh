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
#   - horizontal interpolation done by cdo remapcon, conservative remapping
#   - convert to SERVICE format, calculate meridional integral incl. using
#     land-sea-mask, fortran program psiread
#   - convert badk to NETCDF format using grid description file
#   - simple plot using NCL with shell-script nclsh
#          

set -e

### input parameter

# Resolution:
# cdo remapnn and remapdis generate shaky isolines
# cdo remapcon generates missing values at Poles if resolution of regular target grid is too fine
#   R2B04-r180x90 and R2B05-r360x180 match well

# direct assignment of filenames
#avgfile=avg.5y.xom.bliz.r9368.45
#filestr=avg45.2164ym
#avgfile=timmean_2420-2470
#filestr=2420-2470ym
#resol=r${nlon}x${nlat}
#resol='r180x90'

# via parameter
avgfile=${1}
filestr=${2}
resol=${3:-r180x90}
weightfile=${4}

echo "PSI: Input file is '$avgfile'"
echo "PSI: Tag is '$filestr'"
echo "PSI: Resolution is '$resol'"
echo "PSI: Weightsfile is '$weightfile'"
### input parameter end

# file names
inpfile=uvint.${resol}.$filestr
outfile=psi.${resol}.$filestr
plotfile=nclpsi.$filestr

rlon=${resol%x*}
nlon=${rlon#r}
nlat=${resol#*x}

echo "PSI: Resolution used is $resol"
echo "PSI: Working on files: input: $avgfile ; output: $outfile ; plot: $plotfile"

# select wet_c (lsm at surface) and u_vint (vertical integral of u)
#  - sellevidx: select first level index at surface
if [ ! -s $inpfile.srv ]; then
  if [ -z "$weightfile" ]; then
    echo "PSI: weightfile is not given! Use remapcon for u_vint and remapnn for wet_c ..."
    #cdo -P 8 remapdis,$resol -selvar,u_vint $avgfile xsrc.${resol}_uint.nc
    cdo -P 8 remapcon,$resol -selvar,u_vint $avgfile xsrc.${resol}_uint.nc
    cdo -P 8 remapnn,$resol  -selvar,wet_c -sellevidx,1 $avgfile xsrc.${resol}_wetc.nc
    #cdo merge xsrc.${resol}_wetc.nc xsrc.${resol}_uint.nc $inpfile.nc
    cdo merge xsrc.${resol}_uint.nc xsrc.${resol}_wetc.nc $inpfile.nc
    # Memory fault in r2b5?
    #cdo -P 8 merge -remapnn,$resol  -selvar,wet_c -sellevidx,1 $avgfile -selvar,u_vint $avgfile $inpfile.nc
    rm xsrc.${resol}_wetc.nc xsrc.${resol}_uint.nc
  else
    echo "PSI: weightfile is given:'$weightfile'. Use remap ..."
    cdo -P 8 remap,$resol,$weightfile -selvar,u_vint,wet_c, $avgfile $inpfile.nc
  fi

  # store grid description, convert to service format
  echo "PSI: cdo -v griddes $inpfile.nc > griddes.$resol"
  cdo griddes $inpfile.nc > griddes.$resol

  echo "PSI: cdo -v -f srv copy $inpfile.nc $inpfile.srv"
  cdo -f srv copy $inpfile.nc $inpfile.srv
  rm $inpfile.nc
else
  echo "PSI: Use input for the psi-computation:'$inpfile.srv'"
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
!   ignore vertical dimension
!
! TODO: diffuse from ocean to land, cut land points
! TODO: implement variable output dimension (1 deg resolution) and smoothing extent
!! 
PROGRAM psiread

IMPLICIT NONE

INTEGER, PARAMETER ::  rho_ref = 1025.022            ! reference density

INTEGER, PARAMETER ::  nlat = $nlat                  ! meridional dimension of regular grid
INTEGER, PARAMETER ::  nlon = $nlon                  ! zonal dimension of regular grid

! smoothing area is 2*jsmth-1 lat/lon areas of 1 deg
INTEGER, PARAMETER ::  jsmth = 3                  
INTEGER            :: jb, jc, i_startidx, i_endidx
INTEGER            :: jlat, jlon, jx, jy
INTEGER            :: isrv(8), isrvu(8)


REAL               :: z_lat_dist, erad, pi
REAL               :: z_uint_reg(nlon,nlat)     ! vertical integral on regular grid
REAL               :: psi_reg(nlon,nlat)        ! horizontal stream function
REAL               :: wet_c(nlon,nlat)          ! slm

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

! (2) read barotropic system: after cdo merge above, now read first uint then wet_c

  open (11,file="$inpfile.srv", form='unformatted')
  open (80,file="$outfile.srv", form='unformatted')

  read (11) isrvu
  write (*,*) isrvu
  read (11) z_uint_reg(:,:)
! write(*,*) 'jx=',jx,' jy=',jy,' read uvint=',z_uint_reg(jx,jy)

  read (11) isrv
  read (11) wet_c(:,:)

  write(80) (isrv(jb),jb=1,8)
  write(80) ((wet_c(jlon,jlat),jlon=1,nlon),jlat=1,nlat)


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
  !psi_reg(:,:) = z_uint_reg(:,:) * z_lat_dist * wet_c(:,1,:) * 1.0e-6           ! e+6 [m3/s]
  psi_reg(:,:) = z_uint_reg(:,:) * z_lat_dist * wet_c(:,:) * 1.0e-6


  write(80) (isrvu(jb),jb=1,8)
  write(80) ((psi_reg(jlon,jlat),jlon=1,nlon),jlat=1,nlat)


END PROGRAM psiread
!-------------------------------------------------------------------------  
EOF

echo "PSI: compile and run program scr-psiread.x"
gfortran -o scr-psiread.x scr-psiread.f90
./scr-psiread.x
rm scr-psiread.*

# convert back to netcdf
echo "PSI: cdo -f nc -g griddes.$resol chvar,var4,psi -chvar,var1,wet_c $outfile.srv $outfile.nc"
cdo -f nc -g griddes.$resol chvar,var4,psi -chvar,var1,wet_c $outfile.srv $outfile.nc
rm $outfile.srv

# plot with nclsh:

echo "PSI: plot using icon_plot.ncl:"
nclsh /pool/data/ICON/tools/icon_plot.ncl -altLibDir=/pool/data/ICON/tools \
  -iFile=$outfile.nc -oFile=$plotfile -varName=psi -timeStep=0 -oType=eps \
  -maskName=wet_c -selMode=manual -minVar=-150 -maxVar=150 -numLevs=15 \
  -plotLevs=-150,-100,-75,-50,-30,-20,-10,-5,0,5,10,20,30,50,75,100,150 -withLineLabels

# nclsh $ICONPLOT \
#   -iFile=$outfile.nc -oFile=$plotfile -varName=psi -timeStep=0 -oType=ps \
#   -selMode=manual -minVar=-60 -maxVar=-10 -numLevs=25 -bStrg=' ' -maskName=wet_c \
#   -withLineLabels \

# -mapLLC=-77,20 -mapURC=-70,35 #north-atlantic gyre
#  -mapLLC=120,10 -mapURC=160,40 # japanise WBC
# -plotLevs=-150,-100,-75,-50,-30,-20,-15,-10,-5,0,5,10,15,20,30,50,75,100,150 -withLineLabels
# -plotLevs=-150,-100,-75,-50,-30,-20,-10,-5,0,5,10,20,30,50,75,100,150 -withLineLabels
# -maxView \

# -selMode=manual -minVar=-240 -maxVar=210 -numLevs=15
# -maskName=wet_c -selMode=manual -minVar=-240 -maxVar=210 -numLevs=15 \
# -maskName=wet_c -selMode=manual -minVar=-100 -maxVar=100 -numLevs=20 \

exit

# -selMode=halflog -scaleLimit=1

