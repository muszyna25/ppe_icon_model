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
# cat, yearmean, mulc, netcdf of MOC to scr_moc.ym.nc
files=$(echo $mocPattern | wc -w)
if [[ $files -eq 1 ]]; then
  ln -s $mocPattern scr_moc.ym.nc
else
  cdo -f nc -r -mulc,1.e-9 -yearmean -cat "${mocPattern}"  scr_moc.ym.nc
fi
#  select AMOC at 1000m depth and write ascii data for gnuplot
#if [[ "20" = $(cdo nlevel $(ls -1 $mocPattern | head -n 1)) ]]; then
cdo outputkey,date,value -fldmean -selname,$mocVar -sellonlatbox,0,1,40,60 -sellevel,1000 -setgrid,r1x180 scr_moc.ym.nc > moc.dat
#else
#cdo outputkey,date,value -fldmean -selname,$mocVar -sellonlatbox,0,1,40,60 -intlevel,1000 -sellevel,900/1100 -setgrid,r1x180 scr_moc.ym.nc > moc.dat
#fi
rm scr_moc.ym.nc

gnuplot -p <<EOF
#set term postscript color
#set output 'moc-tims.ps'
#set term png
#set output 'moc-tims.png'

set style data lines
set grid
set title font "Helvetica-Bold,18"
set title  "N-Atlantic Meridional Overturning Circulation (MOC)"
#set title  "N-Atlantic Meridional Overturning Circulation (MOC)" font "Helvetica-Bold,17"
set font "Times-Roman,15"
set xlabel "[simulation years]"
set key left

#set yrange [8.1:27.8]
set ytics 2
set mytics 4

set size 0.95,0.95
set origin 0.05,0.05
set timefmt x "%Y-%m-%d"
set xdata time
set format x "%Y"
#set xtics 20
#set mxtics 5
#plot 'moc.dat' using 1:2 w l t
plot 'moc.dat' using 1:2 w l t 'AMOC, 1000m, 40N-60N' lw 3 lt -1
EOF
