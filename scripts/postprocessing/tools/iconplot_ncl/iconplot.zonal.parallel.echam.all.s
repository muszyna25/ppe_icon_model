#!/bin/bash

set -x

cp iconplot.zonal.parallel.echam.s iconplot.zonal.parallel.echam.list \
   iconplot.zonal-vert.echam.ncl iconplot.zonal-time.echam.ncl \
   iconplot.maps.echam.ncl $TMPDIR
cd $TMPDIR
pwd

/e/uhome/for1han/bin/pshell -p7 -f iconplot.zonal.parallel.echam.list
iconplot.zonal.parallel.echam.s
ncl iconplot.maps.echam.ncl
