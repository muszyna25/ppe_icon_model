#! /bin/bash

INFILE=$1
OUTFILE=$2

cdo -f nc -remap,/pool/data/ICON/grids/private/r2b4_amip/t63grid.nc,/pool/data/ICON/grids/private/r2b4_amip/remapcon_r2b4_amip_cell_t63grid.nc $INFILE $OUTFILE
