#!/usr/bin/env bash

mocPattern='MOC*Z'
if [[ ! -z "$1" ]];then
  mocPattern=$1
fi
mocVar=var778 #AMOC
if [[ ! -z "$2" ]]; then
  mocVar=$2
fi

# calculate the moc timeseries in 1000 depth
cdo -f nc -r -outputkey,date,value -fldmean -selname,$mocVar -mulc,1.e-9 -sellonlatbox,0,1,40,60 -sellevel,1000 -yearmean -setgrid,r1x180 -cat "${mocPattern}"  > moc.dat

gnuplot -p <<EOF
set timefmt x "%Y-%m-%d"
set xdata time
set grid
set format x "%Y"
plot 't.dat' using 1:2 w l t 'AMOC, 1000m, 50N (interpol. from 40N-60N)'
EOF
