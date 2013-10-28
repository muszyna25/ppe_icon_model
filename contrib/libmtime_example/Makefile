BUILDDIR=../../build/x86_64-unknown-linux-gnu/
LIBS+=-L${BUILDDIR}/lib/ -L${BUILDDIR}/externals/mtime/src -lmtime
INCLUDE+=-I${BUILDDIR}/module
all:
	gfortran example.f90 ${LIBS} ${INCLUDE} -o example -g 
