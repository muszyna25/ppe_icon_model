#!/bin/bash

set -eu
#module purge 
#module load PrgEnv-gnu
#module load cray-hdf5
#module load cray-netcdf
#module load libxml2

SCRIPT_DIR=$(cd "$(dirname "$0")"; pwd)
ICON_DIR=$(cd "${SCRIPT_DIR}/../.."; pwd)

MPICH_ROOT='/opt/cray/pe/mpich/8.1.4/ofi/gnu/9.1'
MPICH_LIBS='-lmpifort -lmpi'

HDF5_ROOT=$HDF5_DIR
HDF5_LIBS='-lhdf5'

NETCDF_ROOT='/opt/cray/pe/netcdf-hdf5parallel/default/gnu/9.1'
NETCDF_LIBS='-lnetcdf'

NETCDFF_ROOT='/opt/cray/pe/netcdf-hdf5parallel/default/gnu/9.1'
NETCDFF_LIBS='-lnetcdff'

BLAS_LAPACK_ROOT='/opt/cray/pe/libsci/21.04.1.1/CRAY/9.0/x86_64/'
BLAS_LAPACK_LIBS='-lsci_cray_mpi -lsci_cray'

XML2_ROOT=$LIBXML2_DIR
XML2_LIBS='-lxml2'

###########################
MPI_LAUNCH='/usr/bin/srun'

BUILD_ENV="export LD_LIBRARY_PATH=\"${HDF5_ROOT}/lib:${NETCDF_ROOT}/lib:${NETCDFF_ROOT}/lib:\${LD_LIBRARY_PATH}\";"

CC='cc'
CFLAGS='-g -march=native -mpc64'
ICON_CFLAGS='-O2'
CPPFLAGS="-I${MPICH_ROOT}/include -I${HDF5_ROOT}/include -I${NETCDF_ROOT}/include -I${XML2_ROOT}/include/libxml2"

FC='ftn'
FCFLAGS="-I${MPICH_ROOT}/include -I${NETCDFF_ROOT}/include -std=f2008 -fmodule-private -fimplicit-none -fmax-identifier-length=63 -Wall -Wcharacter-truncation -Wconversion -Wunderflow -Wunused-parameter -Wno-surprising -fall-intrinsics -g -march=native -mpc64"
ICON_FCFLAGS='-fallow-argument-mismatch -fbounds-check -fstack-protector-all -finit-real=nan -finit-integer=-2147483648 -finit-character=127 -O2'

LDFLAGS="-L${MPICH_ROOT}/lib -L${XML2_ROOT}/lib -L${HDF5_ROOT}/lib -L${NETCDF_ROOT}/lib -L${NETCDFF_ROOT}/lib -L${BLAS_LAPACK_ROOT}/lib"
LIBS="-Wl,--as-needed ${XML2_LIBS} ${NETCDFF_LIBS} ${NETCDF_LIBS} ${HDF5_LIBS} ${BLAS_LAPACK_LIBS}"

EXTRA_CONFIG_ARGS='--enable-yaxt --enable-rte-rrtmgp --enable-gpu=no --disable-silent-rules'

"${ICON_DIR}/configure" \
BUILD_ENV="${BUILD_ENV}" \
ARFLAGS="crv" \
CC="${CC}" \
CFLAGS="${CFLAGS}" \
ICON_CFLAGS="${ICON_CFLAGS}" \
CPPFLAGS="${CPPFLAGS}" \
FC="${FC}" \
FCFLAGS="${FCFLAGS}" \
ICON_FCFLAGS="${ICON_FCFLAGS}" \
LDFLAGS="${LDFLAGS}" \
LIBS="${LIBS}" \
MPI_LAUNCH="${MPI_LAUNCH}" \
${EXTRA_CONFIG_ARGS} \
"$@"

for arg in "$@"; do
  case $arg in
    -help | --help | --hel | --he | -h | -help=r* | --help=r* | --hel=r* | --he=r* | -hr* | -help=s* | --help=s* | --hel=s* | --he=s* | -hs*)
      test -n "${EXTRA_CONFIG_ARGS}" && echo '' && echo "This wrapper script ('$0') calls the configure script with the following extra arguments, which might override the default values listed above: ${EXTRA_CONFIG_ARGS}"
      exit 0 ;;
  esac
done

# Generate input file for runscript generation:
./config.status --file=run/set-up.info
echo 'use_mpi_root="openmpi"' >> run/set-up.info

# Copy runscript-related files when building out-of-source:
if test $(pwd) != $(cd "${ICON_DIR}"; pwd); then
  echo "Copying runscript input files from the source directory..."
  rsync -uavz ${ICON_DIR}/run . --exclude='*in' --exclude='.*'
  rsync -uavz ${ICON_DIR}/externals . --exclude='.git' --exclude='*.f90' --exclude='*.F90' --exclude='*.c' --exclude='*.h' --exclude='*.Po' --exclude='tests' --exclude='rrtmgp*.nc' --exclude='*.mod' --exclude='*.o'
  rsync -uavz ${ICON_DIR}/make_runscripts .
  ln -sf ${ICON_DIR}/data
  ln -sf ${ICON_DIR}/vertical_coord_tables
fi

