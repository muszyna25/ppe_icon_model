#!/bin/bash

set -x

cp iconplot.zonal.parallel.s iconplot.zonal.parallel.list \
   iconplot.zonal-vert.ncl iconplot.zonal-time.ncl \
   iconplot.maps.ncl $TMPDIR
cd $TMPDIR
pwd

/e/uhome/for1han/bin/pshell -p7 -f iconplot.zonal.parallel.list
iconplot.zonal.parallel.s
ncl iconplot.maps.ncl
