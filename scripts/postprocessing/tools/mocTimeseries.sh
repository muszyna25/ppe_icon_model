#!/usr/bin/env bash

set -x

# file patterns:
#  - 'MOC*Z' : if more than one file, files are catted to one annual mean file 'moc*'
#  - moc*    : one annual mean file to be plotted
mocPattern='MOC*Z'
if [[ ! -z "$1" ]];then
  mocPattern=$1
fi

mocVar=var778 #AMOC
if [[ ! -z "$2" ]]; then
  mocVar=$2
fi

# calculate the moc timeseries in 1000 depth
# cat, yearmean, netcdf of MOC to moc.ym.nc
files=$(echo $mocPattern | wc -w)
if [[ $files -gt 1 ]]; then
  cdo -f nc -r -yearmean -cat "${mocPattern}"  moc.ym.nc
  mocPattern=moc.ym.nc
fi

# "moc.dat" - plot only:
if [[ "$mocPattern" != "moc.dat" ]]; then

  #  select AMOC at 1000m depth and write ascii data for gnuplot
  if [[ "20" = $(cdo nlevel -selvar,$mocVar $(ls -1 $mocPattern | head -n 1)) ]]; then
    cdo outputkey,date,value -mulc,1.e-9 -fldmean -selname,$mocVar -sellonlatbox,0,1,40,60 \
         -sellevel,1000 -setgrid,r1x180 $mocPattern > moc.dat
  else
    cdo outputkey,date,value -mulc,1.e-9 -fldmean -selname,$mocVar -sellonlatbox,0,1,40,60 \
         -intlevel,1000 -sellevel,900/1100 -setgrid,r1x180 $mocPattern > moc.dat
  fi
fi

gnuplot -p <<EOF
#set terminal postscript color
#set output 'moc-tims.ps'
#set terminal png size 1000,600
#set output 'moc-tims.png'

set size 0.98,0.95
set origin 0.02,0.02

set style data lines
set grid
set title font "Helvetica-Bold,16"
set title  "N-Atlantic Meridional Overturning Circulation (MOC)"
#set title  "N-Atlantic Meridional Overturning Circulation (MOC)" font "Helvetica-Bold,17"

set timefmt x "%Y-%m-%d"
set xdata time
set format x "%Y"
#set xtics 20
#set mxtics 5
set xlabel "simulation years" font "Times-Bold,14

#set yrange [8.1:27.8]
set ytics 4
set mytics 4   # minor tickmarks
set ylabel "[SV]" font "Times-Bold,14

set key left
#plot 'moc.dat' using 1:2 w l t
plot 'moc.dat' using 1:2 w l t 'AMOC, 1000m, 40N-60N' lw 3 lt -1
EOF

#gwenview moc-tims.png
