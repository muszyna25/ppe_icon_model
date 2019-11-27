#!/bin/bash

set -eu

MY_DIR=$(cd "$(dirname "$0")"; pwd)
ICON_DIR="${MY_DIR}/../../.."

. "${MY_DIR}/init_intel_16.0.2.sh"

FC=$MPICH_FC
CC=$PLAIN_CC

FCFLAGS="${NETCDFF_FCFLAGS} ${CDI_FCFLAGS} ${MTIME_FCFLAGS} ${YAC_FCFLAGS} ${SELF_FCFLAGS}"
ICON_FCFLAGS='-O2 -g'
CFLAGS='-gcc-name=/usr/bin/gcc'
ICON_CFLAGS='-O2 -g'
CPPFLAGS=
LDFLAGS="${NETCDF_LDFLAGS} ${NETCDFF_LDFLAGS} ${CDI_LDFLAGS} ${MTIME_LDFLAGS} ${BLAS_LDFLAGS} ${LAPACK_LDFLAGS} ${XML2_LDFLAGS} ${YAC_LDFLAGS} ${SELF_LDFLAGS}"
LIBS="${SELF_LIBS} ${YAC_LIBS} ${XML2_LIBS} ${LAPACK_LIBS} ${BLAS_LIBS} ${MTIME_LIBS} ${CDI_LIBS} ${NETCDFF_LIBS} ${NETCDF_LIBS}"

# Set LD_LIBRARY_PATH inside BUILD_ENV to enable the location of dynamic
# libraries at the configure time and when running 'make check':
BUILD_ENV="$BUILD_ENV LD_LIBRARY_PATH=\"`echo " ${LDFLAGS}" | sed -e 's/[ ][ ]*/ /g;s/-L[ ]/-L/g;s/[ ]\(-[^L]\|[^-]\)[^ ]*//g;s/[ ]*-L/:/g;s/^://'`:\$LD_LIBRARY_PATH\"; export LD_LIBRARY_PATH;"

"${ICON_DIR}/configure" \
BUILD_ENV="$BUILD_ENV" \
CC="$CC" \
CFLAGS="$CFLAGS" \
CPPFLAGS="$CPPFLAGS" \
FC="$FC" \
FCFLAGS="$FCFLAGS" \
ICON_CFLAGS="$ICON_CFLAGS" \
ICON_FCFLAGS="$ICON_FCFLAGS" \
LDFLAGS="$LDFLAGS" \
LIBS="$LIBS" \
--enable-explicit-fpp \
--with-external-cdi \
--with-external-mtime \
--with-external-self \
--with-external-yac \
"$@"

