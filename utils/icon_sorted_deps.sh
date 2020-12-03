#!/bin/sh

MY_DIR=`dirname "$0"`

"${MY_DIR}/mkhelper/deplist.py" -t icon -f - <<_EOF
icon: cuda cub mpi netcdf-fortran rte-rrtmgp ecrad rttov sct yaxt cdi serialbox2 mtime blas lapack yac tixi self eccodes hdf5 zlib
netcdf-fortran: netcdf
netcdf: hdf5 zlib
cdi: eccodes netcdf aec mpi yaxt
yac: lapack mtime xml2 netcdf mpi
lapack: blas
sct: hdf5 mpi
hdf5: mpi aec zlib
yaxt: mpi
tixi: xml2
rttov: netcdf-fortran hdf5
ecrad: netcdf-fortran
rte-rrtmgp: netcdf-fortran
serialbox2: netcdf stdc++
eccodes: aec
cub: cuda stdc++
cuda: stdc++
xml2: zlib
_EOF

