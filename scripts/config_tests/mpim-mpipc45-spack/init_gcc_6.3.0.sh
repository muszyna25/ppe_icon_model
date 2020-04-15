BUILD_ENV=

PLAIN_FC='/usr/bin/gfortran'
PLAIN_CC='/usr/bin/gcc'

SW_ROOT='/scratch/local1/icon-bootstrap/opt/gcc-6.3.0'

MPICH_ROOT="${SW_ROOT}/mpich-3.2.1-7rbrrmk"
MPICH_FC="${MPICH_ROOT}/bin/mpif90"
MPICH_CC="${MPICH_ROOT}/bin/mpicc"
MPICH_FCFLAGS="-I${MPICH_ROOT}/include"
MPICH_CPPFLAGS="-I${MPICH_ROOT}/include"
MPICH_LDFLAGS="-L${MPICH_ROOT}/lib"
MPICH_LIBS='-lmpifort -lmpi'
MPICH_LAUNCH="${MPICH_ROOT}/bin/mpirun"

HDF5_ROOT="${SW_ROOT}/hdf5-1.10.3-45iykje"
HDF5_CPPFLAGS="-I${HDF5_ROOT}/include"
HDF5_LDFLAGS="-L${HDF5_ROOT}/lib"
HDF5_LIBS='-lhdf5'

HDF5_PARALLEL_ROOT="${SW_ROOT}/hdf5-1.10.3-yjykbdp"
HDF5_PARALLEL_CPPFLAGS="-I${HDF5_PARALLEL_ROOT}/include"
HDF5_PARALLEL_LDFLAGS="-L${HDF5_PARALLEL_ROOT}/lib"
HDF5_PARALLEL_LIBS='-lhdf5'

NETCDF_ROOT="${SW_ROOT}/netcdf-4.6.1-k4x5yne"
NETCDF_CPPFLAGS="-I${NETCDF_ROOT}/include"
NETCDF_LDFLAGS="-L${NETCDF_ROOT}/lib"
NETCDF_LIBS='-lnetcdf'

NETCDF_PARALLEL_ROOT="${SW_ROOT}/netcdf-4.6.1-ufgmniu"
NETCDF_PARALLEL_CPPFLAGS="-I${NETCDF_PARALLEL_ROOT}/include"
NETCDF_PARALLEL_LDFLAGS="-L${NETCDF_PARALLEL_ROOT}/lib"
NETCDF_PARALLEL_LIBS='-lnetcdf'

NETCDFF_ROOT="${SW_ROOT}/netcdf-fortran-4.4.4-aaix567"
NETCDFF_FCFLAGS="-I${NETCDFF_ROOT}/include"
NETCDFF_LDFLAGS="-L${NETCDFF_ROOT}/lib"
NETCDFF_LIBS='-lnetcdff'

NETCDFF_PARALLEL_ROOT="${SW_ROOT}/netcdf-fortran-4.4.4-hc2po4m"
NETCDFF_PARALLEL_FCFLAGS="-I${NETCDFF_PARALLEL_ROOT}/include"
NETCDFF_PARALLEL_LDFLAGS="-L${NETCDFF_PARALLEL_ROOT}/lib"
NETCDFF_PARALLEL_LIBS='-lnetcdff'

RTE_RRTMGP_ROOT="${SW_ROOT}/rte-rrtmgp-develop-obwldg6"
RTE_RRTMGP_FCFLAGS="-I${RTE_RRTMGP_ROOT}/include"
RTE_RRTMGP_LDFLAGS="-L${RTE_RRTMGP_ROOT}/lib"
RTE_RRTMGP_LIBS='-lrrtmgp -lrte'

ECRAD_ROOT="${SW_ROOT}/ecrad-1.1.0-bjpih4b"
ECRAD_FCFLAGS="-I${ECRAD_ROOT}/include"
ECRAD_LDFLAGS="-L${ECRAD_ROOT}/lib"
ECRAD_LIBS='-lradiation -lifsrrtm -lutilities -lifsaux'

SCT_ROOT="${SW_ROOT}/libsct-develop-b7cvy6w"
SCT_FCFLAGS="-I${SCT_ROOT}/include"
SCT_LDFLAGS="-L${SCT_ROOT}/lib"
SCT_LIBS='-lsct'

YAXT_ROOT="${SW_ROOT}/yaxt-develop-j74yauy"
YAXT_FCFLAGS="-I${YAXT_ROOT}/include"
YAXT_LDFLAGS="-L${YAXT_ROOT}/lib"
YAXT_LIBS='-lyaxt'

ECCODES_ROOT="${SW_ROOT}/eccodes-2.5.0-g6ta3it"
ECCODES_CPPFLAGS="-I${ECCODES_ROOT}/include"
ECCODES_LDFLAGS="-L${ECCODES_ROOT}/lib"
ECCODES_LIBS='-leccodes'

GRIBAPI_ROOT="${SW_ROOT}/grib-api-1.24.0-35xwc3t"
GRIBAPI_CPPFLAGS="-I${GRIBAPI_ROOT}/include"
GRIBAPI_LDFLAGS="-L${GRIBAPI_ROOT}/lib"
GRIBAPI_LIBS='-lgrib_api'

CDI_ROOT="${SW_ROOT}/libcdi-1.8.2-f52z2uw"
CDI_FCFLAGS="-I${CDI_ROOT}/include"
CDI_LDFLAGS="-L${CDI_ROOT}/lib"
CDI_LIBS='-lcdi_f2003 -lcdi'

CDI_PIO_ROOT="${SW_ROOT}/libcdi-pio-develop-bx5jzjx"
CDI_PIO_FCFLAGS="-I${CDI_PIO_ROOT}/include"
CDI_PIO_LDFLAGS="-L${CDI_PIO_ROOT}/lib"
CDI_PIO_LIBS='-lcdi_f2003 -lcdipio -lcdi'

SERIALBOX2_ROOT="${SW_ROOT}/serialbox2-2.5.3-qhmkdrw"
SERIALBOX2_FCFLAGS="-I${SERIALBOX2_ROOT}/include"
SERIALBOX2_LDFLAGS="-L${SERIALBOX2_ROOT}/lib"
SERIALBOX2_LIBS='-lSerialboxFortranShared'

MTIME_ROOT="${SW_ROOT}/libmtime-1.0.8-p1-zwn6cyh"
MTIME_FCFLAGS="-I${MTIME_ROOT}/include"
MTIME_CPPFLAGS="-I${MTIME_ROOT}/include"
MTIME_LDFLAGS="-L${MTIME_ROOT}/lib"
MTIME_LIBS='-lmtime'

BLAS_ROOT="${SW_ROOT}/netlib-lapack-3.8.0-aayvlps"
BLAS_LDFLAGS="-L${BLAS_ROOT}/lib"
BLAS_LIBS='-lblas'

LAPACK_ROOT="${SW_ROOT}/netlib-lapack-3.8.0-aayvlps"
LAPACK_LDFLAGS="-L${LAPACK_ROOT}/lib"
LAPACK_LIBS='-llapack'

XML2_ROOT="${SW_ROOT}/libxml2-2.9.8-blvm4bl"
XML2_CPPFLAGS="-I${XML2_ROOT}/include/libxml2"
XML2_LDFLAGS="-L${XML2_ROOT}/lib"
XML2_LIBS='-lxml2'

YAC_ROOT="${SW_ROOT}/yac-1.5.2-rc-cg6mmw7"
YAC_FCFLAGS="-I${YAC_ROOT}/include"
YAC_LDFLAGS="-L${YAC_ROOT}/lib"
YAC_LIBS='-lyac'

TIXI_ROOT="${SW_ROOT}/icon-tixi-develop-oxmq44j"
TIXI_FCFLAGS="-I${TIXI_ROOT}/include"
TIXI_LDFLAGS="-L${TIXI_ROOT}/lib"
TIXI_LIBS='-licon_tixi'

SELF_ROOT="${SW_ROOT}/libself-0.2-fq3unll"
SELF_FCFLAGS="-I${SELF_ROOT}/include"
SELF_LDFLAGS="-L${SELF_ROOT}/lib"
SELF_LIBS='-lself'

CLAW_ROOT='/scratch/local1/claw/install'
CLAW="${CLAW_ROOT}/bin/clawfc"