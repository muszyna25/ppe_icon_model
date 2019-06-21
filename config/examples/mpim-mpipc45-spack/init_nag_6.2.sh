BUILD_ENV='export NAG_KUSARI_FILE="license.mpimet.mpg.de:";'

PLAIN_FC='/sw/stretch-x64/nag/nag-6.2/bin/nagfor'
PLAIN_CC='/usr/bin/gcc'

SW_ROOT='/scratch/local1/icon-bootstrap/opt/nag-6.2'

MPICH_ROOT="${SW_ROOT}/mpich-3.2.1-w7ksmw2"
MPICH_FC="${MPICH_ROOT}/bin/mpif90"
MPICH_CC="${MPICH_ROOT}/bin/mpicc"
MPICH_FCFLAGS="-I${MPICH_ROOT}/include"
MPICH_CPPFLAGS="-I${MPICH_ROOT}/include"
MPICH_LDFLAGS="-L${MPICH_ROOT}/lib"
MPICH_LIBS='-lmpifort -lmpi'
MPICH_LAUNCH="${MPICH_ROOT}/bin/mpirun"

HDF5_ROOT="${SW_ROOT}/hdf5-1.10.3-aehnkhs"
HDF5_CPPFLAGS="-I${HDF5_ROOT}/include"
HDF5_LDFLAGS="-L${HDF5_ROOT}/lib"
HDF5_LIBS='-lhdf5'

NETCDF_ROOT="${SW_ROOT}/netcdf-4.6.1-yalf6ng"
NETCDF_CPPFLAGS="-I${NETCDF_ROOT}/include"
NETCDF_LDFLAGS="-L${NETCDF_ROOT}/lib"
NETCDF_LIBS='-lnetcdf'

NETCDFF_ROOT="${SW_ROOT}/netcdf-fortran-4.4.4-63bnpcj"
NETCDFF_FCFLAGS="-I${NETCDFF_ROOT}/include"
NETCDFF_LDFLAGS="-L${NETCDFF_ROOT}/lib"
NETCDFF_LIBS='-lnetcdff'

RTE_RRTMGP_ROOT="${SW_ROOT}/rte-rrtmgp-develop-yu5t3nd"
RTE_RRTMGP_FCFLAGS="-I${RTE_RRTMGP_ROOT}/include"
RTE_RRTMGP_LDFLAGS="-L${RTE_RRTMGP_ROOT}/lib"
RTE_RRTMGP_LIBS='-lrrtmgp -lrte'

ECRAD_ROOT="${SW_ROOT}/ecrad-1.1.0-mc3oq5g"
ECRAD_FCFLAGS="-I${ECRAD_ROOT}/include"
ECRAD_LDFLAGS="-L${ECRAD_ROOT}/lib"
ECRAD_LIBS='-lradiation -lifsrrtm -lutilities -lifsaux'

SCT_ROOT="${SW_ROOT}/libsct-develop-4srm5jo"
SCT_FCFLAGS="-I${SCT_ROOT}/include"
SCT_LDFLAGS="-L${SCT_ROOT}/lib"
SCT_LIBS='-lsct'

YAXT_ROOT="${SW_ROOT}/yaxt-develop-bstwu2r"
YAXT_FCFLAGS="-I${YAXT_ROOT}/include"
YAXT_LDFLAGS="-L${YAXT_ROOT}/lib"
YAXT_LIBS='-lyaxt'

GRIBAPI_ROOT="${SW_ROOT}/grib-api-1.24.0-i3glwdl"
GRIBAPI_CPPFLAGS="-I${GRIBAPI_ROOT}/include"
GRIBAPI_LDFLAGS="-L${GRIBAPI_ROOT}/lib"
GRIBAPI_LIBS='-lgrib_api'

CDI_ROOT="${SW_ROOT}/libcdi-1.8.2-j6ugmwh"
CDI_FCFLAGS="-I${CDI_ROOT}/include"
CDI_LDFLAGS="-L${CDI_ROOT}/lib"
CDI_LIBS='-lcdi_f2003 -lcdi'

MTIME_ROOT="${SW_ROOT}/libmtime-1.0.8-p1-psgf6sn"
MTIME_FCFLAGS="-I${MTIME_ROOT}/include"
MTIME_CPPFLAGS="-I${MTIME_ROOT}/include"
MTIME_LDFLAGS="-L${MTIME_ROOT}/lib"
MTIME_LIBS='-lmtime'

BLAS_ROOT="${SW_ROOT}/netlib-lapack-3.8.0-footyva"
BLAS_LDFLAGS="-L${BLAS_ROOT}/lib"
BLAS_LIBS='-lblas'

LAPACK_ROOT="${SW_ROOT}/netlib-lapack-3.8.0-footyva"
LAPACK_LDFLAGS="-L${LAPACK_ROOT}/lib"
LAPACK_LIBS='-llapack'

XML2_ROOT="${SW_ROOT}/libxml2-2.9.8-chuop5e"
XML2_CPPFLAGS="-I${XML2_ROOT}/include/libxml2"
XML2_LDFLAGS="-L${XML2_ROOT}/lib"
XML2_LIBS='-lxml2'

YAC_ROOT="${SW_ROOT}/yac-1.5.2-rc-pstruto"
YAC_FCFLAGS="-I${YAC_ROOT}/include"
YAC_LDFLAGS="-L${YAC_ROOT}/lib"
YAC_LIBS='-lyac'

SELF_ROOT="${SW_ROOT}/libself-0.2-iefxuap"
SELF_FCFLAGS="-I${SELF_ROOT}/include"
SELF_LDFLAGS="-L${SELF_ROOT}/lib"
SELF_LIBS='-lself'

