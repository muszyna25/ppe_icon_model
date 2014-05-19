#!/usr/bin/env bash

set -ex

# calculates sea ice thickness and snow times ice concentration for NH and SH
# plots thickness and concentration for March/September using nclsh
# arranges 4 plots for SH and NH together using convert
#  - update path etc.
#  - SLO 29.04.2014

expfile=$1    # long filename
suffix=$2     # nh or sh
selyear=2023
echo $selyear
#datfile=dat.r16207.def.icets.$suffix.${selyear}y.nc
expname=$(basename $(pwd))
#brfname=${expname##.}
#brfname=r16289.leadcl025
#brfname=r16427.icdyn.lc25
#brfname=r16462.iceth.lc25
brfname=r17044.def
#datfile=dat.$expname.icets.$suffix.${selyear}y.nc
datfile=dat.$brfname.icets.$suffix.${selyear}y.nc
cal=y      # calculate nh or sh ice + ts variables
ncsl=y     # plot using nclsh
conv=y     # convert plots 4 in 1

vi=hixc    # thickness times concentration
vi=hi      # thickness only


nclsh="nclsh /scratch/mpi/CC/mh0287/users/m211032/Icon/icon-dev.tst/scripts/postprocessing/tools/icon_plot.ncl -altLibDir=/scratch/mpi/CC/mh0287/users/m211032/Icon/icon-dev.tst/scripts/postprocessing/tools"
nclsh="nclsh /pool/data/ICON/tools/icon_plot.ncl -altLibDir=/pool/data/ICON/tools"

if [[ $cal = y ]]; then

  if [[ "$suffix" = "nh" ]]; then
   cdo -sellonlatbox,0,360,50,90   \
       -selvar,wet_c,t_acc,s_acc,hi_acc,hs_acc,conc_acc,ice_u_acc,ice_v_acc \
       -selyear,$selyear $expfile xscr.icets.nc
  fi

  if [[ "$suffix" = "sh" ]]; then
   cdo -sellonlatbox,0,360,-90,-50 \
       -selvar,wet_c,t_acc,s_acc,hi_acc,hs_acc,conc_acc,ice_u_acc,ice_v_acc \
       -selyear,$selyear $expfile xscr.icets.nc
  fi

  cdo chname,hi_acc,hixc_acc -mul -selvar,hi_acc xscr.icets.nc -selvar,conc_acc xscr.icets.nc xscr.hiXconc.nc
  cdo chname,hs_acc,hsxc_acc -mul -selvar,hs_acc xscr.icets.nc -selvar,conc_acc xscr.icets.nc xscr.hsXconc.nc
  rm -f xscr.hihsXconc.nc
  cdo merge xscr.hiXconc.nc xscr.hsXconc.nc xscr.hihsXconc.nc
  rm -f $datfile
  cdo merge xscr.icets.nc xscr.hihsXconc.nc $datfile
  rm -f xscr.*nc
  
fi  # calculate
  
# plot mit icon_plot.ncl:

if [[ $ncsl = y ]]; then

  ts=2       # timestep
  vari=${vi}_acc
  iceplot=pl.$brfname.$suffix.${selyear}y.$vi.$ts
  if [[ "$suffix" = "nh" ]]; then mapt=NHps; fi
  if [[ "$suffix" = "sh" ]]; then mapt=SHps; fi
  
  $nclsh \
     -iFile=$datfile -oFile=$iceplot -timeStep=$ts -mapType=$mapt -varName=$vari \
     -plotLevs=0.01,0.1,0.3,0.5,0.7,1,1.5,2,3,4,5 \
     -colormap=WhiteBlueGreenYellowRed \
  #  -showGrid \
  
  convert -density 90 $iceplot.eps $iceplot.png
  
  ts=8
  iceplot=pl.$brfname.$suffix.${selyear}y.$vi.$ts

  $nclsh \
     -iFile=$datfile -oFile=$iceplot -timeStep=$ts -mapType=$mapt -varName=$vari \
     -plotLevs=0.01,0.1,0.3,0.5,0.7,1,1.5,2,3,4,5 \
     -colormap=WhiteBlueGreenYellowRed \
  #  -showGrid \

  convert -density 90 $iceplot.eps $iceplot.png
  
  ts=2
  vc=conc
  varc=${vc}_acc
  iceplot=pl.$brfname.$suffix.${selyear}y.$vc.$ts
  
  $nclsh \
     -iFile=$datfile -oFile=$iceplot -timeStep=$ts -mapType=$mapt -varName=$varc \
     -plotLevs=0.01,0.05,0.1,0.3,0.5,0.7,0.8,0.9,0.99 \
     -colormap=WhiteBlueGreenYellowRed \
  #  -showGrid \
  
  convert -density 90 $iceplot.eps $iceplot.png
  
  ts=8
  iceplot=pl.$brfname.$suffix.${selyear}y.$vc.$ts
  
  $nclsh \
     -iFile=$datfile -oFile=$iceplot -timeStep=$ts -mapType=$mapt -varName=$varc \
     -plotLevs=0.01,0.05,0.1,0.3,0.5,0.7,0.8,0.9,0.99 \
     -colormap=WhiteBlueGreenYellowRed \
  #  -showGrid \
  
  convert -density 90 $iceplot.eps $iceplot.png

fi  # ncsl

if [[ $conv = y ]]; then

  x1dat=xscr1.png
  x2dat=xscr2.png
  splnam=pl.$brfname.$suffix.${selyear}y
  quadnhigh=quadpl.high.$brfname.$suffix.${selyear}y.png
  quadnlow=quadpl.low.$brfname.$suffix.${selyear}y.png
  
  convert +append $splnam.${vc}.2.png $splnam.${vc}.8.png $x1dat
  convert +append $splnam.${vi}.2.png $splnam.${vi}.8.png $x2dat
  convert -geometry 1200x1200 -append $x1dat $x2dat $quadnhigh
  convert -geometry  700x700  -append $x1dat $x2dat $quadnlow

fi  # conv
  
exit
  

  
  
  
  
  
  
  
