#!/bin/ksh

#  compdev_mpipc.ksh
#  Simple script to quickly compile and link files to update solely the executable
#  "control_model" without compiling sources other than specified

#  Author: Stephan Lorenz - 2010/02
#  Modified - 2010/11

#   - much faster than "make" in case of simple development of single routines
#   - runs after successful configure/make to update control_model
#   - compiler (gcc/nag/intel/pgi) and several files can be chosen as input
#     where the path down of src must be specified
#   - without arguments only the linker is invoked

# Compiles a module file for control_model
#  - caution: new module binaries are written to src/module.o, i.e. not to modules/module.mod 
#  - caution: compile as much as necessary since no dependicies as in make are known
# Revisions:
#  - now working directly in icon-dev/src - compiling only two sources (2010-06-07)
#  - now reading objects from Makefile (2010-08-02)
#  - now compiling sources from input arguments (2010-11-18)
#  - now running on all system via $arch (2010-11-18)
#  - also on multiple builds (2011-02-09)

set -e
# you must start this script from below the trunk/branch - e.g. from icon-dev/run
# in case of multiple builds start below the specific dir - e.g. icon-dev/nag_nMnO/run
dir=$(pwd -P)
basedir=${dir%/*}
# determine architecture
arch=`ls ${basedir}/build`

# go to basedir
cd ..
# get compiler from Makefile
fccomp=$(grep FC Makefile)
COMP=${fccomp##* }

# go to src directory:
cd src

# objects to link are the same for all compilers but may change with revisions:

# get objects from Makefile
OBS=$(sed -n '/OBJS =/,/^$/p' ../build/$arch/src/Makefile | sed 's/OBJS =//' | sed s'/\\//')
objects="-o ../bin/control_model control_model.o $OBS"
#echo objects are: $objects

if [ $COMP == "ifort" ] ; then

  # compiler and options intel:
  compopt="ifort -I../include -I/sw/etch-ia32/netcdf-3.6.3/include    -I/sw/etch-ia32/mpich2-1.2.1-intel11/include   -module ../module -I../module -O3 -msse2 -mieee-fp -fpe0 -pc64 -fpp -traceback -DNOMPI -c ../../../src"

  # loader with objects and options intel:
  loadobj="ifort -I../module -O3 -msse2 -mieee-fp -fpe0 -pc64 -fpp -traceback $objects -L/sw/etch-ia32/mpich2-1.2.1-intel11/lib -lmpichf90 -lmpich -lpthread -lrt"

elif [ $COMP == "pgi" ] ; then

  # compiler and options pgi:
  compopt="pgf95 -I../include -I/sw/etch-ia32/netcdf-3.6.3/include -I/sw/etch-ia32/mpich2-1.2.1-pgi9/include -module ../module -I../module -fast -fastsse -g -Mpreprocess -Mchkstk  -Mrecursive -c ../../../src"

  # loader with objects and options pgi:
  loadobj="pgf95 -I../module -fast -fastsse -g -Mpreprocess -Mchkstk  -Mrecursive $objects
  /zmaw/sw/sw/etch-ia32/gcc-4.3.3/bin/../lib/gcc/$arch/4.3.3/libgcc.a    -L/sw/etch-ia32/mpich2-1.2.1-pgi9/lib -lmpichf90 -lmpich -lpthread -lrt"


elif [ $COMP == "nagfor" ] ; then

  # compiler and options nag:
  compopt="nagfor -I../include -I/sw/etch-ia32/netcdf-3.6.3/include    -I/sw/etch-ia32/mpich2-1.2.1-nag52/include   -mdir ../module -I../module -float-store -gline -g -maxcontin=99 -fpp -f2003 -nan -C=all -wmismatch=mpi_send,mpi_isend,mpi_recv,mpi_irecv,mpi_bcast,mpi_allreduce,nf_put_att_double,nf_def_var,nf_get_att_int -w=uep  -c ../../../src"

  # loader with objects and options nag:
  loadobj="nagfor -I../module -float-store -gline -g -maxcontin=99 -fpp -f2003 -C=all
  -wmismatch=mpi_send,mpi_isend,mpi_recv,mpi_irecv,mpi_bcast,mpi_allreduce,nf_put_att_double,nf_def_var,nf_get_att_int
  -w=uep $objects -L../lib -lsupport  -L../lib -llapack -lblas -L/sw/etch-ia32/netcdf-3.6.3/lib
  -lnetcdf /zmaw/sw/sw/etch-ia32/gcc-4.3.3/bin/../lib/gcc/$arch/4.3.3/libgcc.a    -L/sw/etch-ia32/mpich2-1.2.1-nag52/lib -lmpichf90 -lmpich -lpthread -lrt"



elif [ $COMP == "gfortran" ] ; then

  # compiler and options gcc:
  compopt="gfortran -I../include -I/sw/etch-ia32/netcdf-3.6.3/include    -I/sw/etch-ia32/mpich2-1.2.1-gcc43/include   -J../module -I../module -march=native -O -ffast-math -xf95-cpp-input -std=f2003 -fmodule-private -fimplicit-none -fmax-identifier-length=31 -ffree-line-length-99 -Wall -Wcharacter-truncation -Wconversion -Wunderflow -Wunused-parameter -fbacktrace -finit-real=nan -c ../../../src"


  # loader with objects and options gcc:
  loadobj="gfortran -I../module -march=native -O -ffast-math  -std=f2003 -fmodule-private
  -fimplicit-none -fmax-identifier-length=31 -ffree-line-length-99 -Wall -Wcharacter-truncation
  -Wconversion -Wunderflow -Wunused-parameter -fbacktrace -finit-real=nan $objects -L../lib
  -lsupport  -L../lib -llapack -lblas -L/sw/etch-ia32/netcdf-3.6.3/lib -lnetcdf
  /zmaw/sw/sw/etch-ia32/gcc-4.3.3/bin/../lib/gcc/$arch/4.3.3/libgcc.a    -L/sw/etch-ia32/mpich2-1.2.1-gcc43/lib -lmpichf90 -lmpich -lpthread -lrt"

fi

# change to build directory - starting from below the trunk/branch - e.g. from run
cd ../build/$arch/src

# do compiling and loading
#  - files to compile in input
for i do
  echo " *** *** *** compiler $COMP compiles file $i"
  $compopt/$i
  shift
done

echo " *** *** *** linking objects to control_model"
$loadobj

