export

SHELL = /bin/sh

ARCH = x86_64
OS   = darwin11.4.2

HOST = wanglung
SITE = zmaw.de

srcdir = .


prefix = .
exec_prefix = build/x86_64-apple-darwin11.4.2

bindir = ${exec_prefix}/bin
libdir = ${exec_prefix}/lib

NETCDFROOT     = /opt/local
NETCDFLIBPATH  = /opt/local/lib
NETCDF_LIB     = -L$(NETCDFLIBPATH) -lnetcdff -lnetcdf
NETCDF_INCLUDE = -I$(NETCDFROOT)/include

HDF5ROOT       = /opt/local
HDF5_LIB       = -L$(HDF5ROOT)/lib -lhdf5_hl -lhdf5
HDF5_INCLUDE   = -I$(HDF5ROOT)/include

SZIPROOT       = /opt/local
SZIP_LIB       = -L$(SZIPROOT)/lib -lsz
SZIP_INCLUDE   = -I$(SZIPROOT)/include

ZLIBROOT       = /opt/local
ZLIB_LIB       = -L$(ZLIBROOT)/lib -lz
ZLIB_INCLUDE   = -I$(ZLIBROOT)/include

MPIROOT        = /opt/local
MPI_LIB        = -L$(MPIROOT)/lib -lmpi_f90 -lmpi_f77 -lmpi -lopen-rte -lopen-pal -lutil
MPI_INCLUDE    = -I$(MPIROOT)/lib -I$(MPIROOT)/include/openmpi/

LAPACKROOT     = 
LAPACK_LIB_PATH= 
LAPACK_LIB     = -L../lib -llapack -lblas

PROFILE_LIB     = 
PROFILE_INCLUDE = 

METIS_LIB      = 

LIBS           = -L../lib -lsupport -lcurl $(LAPACK_LIB) $(NETCDF_LIB) $(HDF5_LIB) $(SZIP_LIB) $(ZLIB_LIB) $(MPI_LIB) $(METIS_LIB) $(PROFILE_LIB)
INCLUDE        = -I../include -I../../../src/include $(NETCDF_INCLUDE) $(HDF5_INCLUDE) $(SZIP_INCLUDE) $(ZLIB_INCLUDE) $(MPI_INCLUDE) $(PROFILE_INCLUDE)
INCLUDES       = $(INCLUDE)

AS             = as

CC             = gcc
CFLAGS         = $(INCLUDE) -march=core2 -O -DpgiFortran -DHAVE_LIBNETCDF -DHAVE_CF_INTERFACE -std=gnu99 -g -D__ICON__
FC             = gfortran
FFLAGS         = $(INCLUDES) -J../module -I../module -march=core2 -O -ffast-math -xf95-cpp-input -std=f2003 -fmodule-private -fimplicit-none -fmax-identifier-length=63 -ffree-line-length-132 -Wall -Wcharacter-truncation -Wconversion -Wunderflow -Wunused-parameter -g -fbacktrace -fbounds-check -finit-real=nan -finit-integer=-2147483648 -finit-character=127 -D__ICON__ 
F77            = gfortran
F77FLAGS       = -march=core2 -O -ffast-math

AR             = ar
ARFLAGS        = crv

LDFLAGS        = -I../module -march=core2 -O -ffast-math  -std=f2003 -fmodule-private -fimplicit-none -fmax-identifier-length=63 -ffree-line-length-132 -Wall -Wcharacter-truncation -Wconversion -Wunderflow -Wunused-parameter -g -fbacktrace -fbounds-check -finit-real=nan -finit-integer=-2147483648 -finit-character=127 -D__ICON__

SRCDIRS        =  blas lapack support src
OBJDIRS        =  build/x86_64-apple-darwin11.4.2/blas build/x86_64-apple-darwin11.4.2/lapack build/x86_64-apple-darwin11.4.2/support build/x86_64-apple-darwin11.4.2/src

.PHONY: doc

all:
	@for DIR in $(OBJDIRS) ;\
	  do \
	    back=`pwd`; \
	    cd $$DIR ;\
	    $(MAKE) ; status=$$? ; \
	    if [ $$status != 0 ] ; then \
	      echo "Exit status from make was $$status" ; exit $$status ; \
	    fi ; \
	    cd $$back ; \
	  done 
control:
	@for DIR in $(OBJDIRS) ;\
	  do LASTDIR=$$DIR ;\
	done ;\
	back=`pwd` ;\
	cd $$LASTDIR ;\
	$(MAKE) control_model  ;\
	cd $$back

one:
	@for DIR in $(OBJDIRS) ;\
	  do LASTDIR=$$DIR ;\
	done ;\
	back=`pwd` ;\
	cd $$LASTDIR ;\
	$(MAKE) $(name)  ;\
	cd $$back


install:
	@for DIR in $(OBJDIRS) ;\
	  do \
	  (cd $$DIR ;\
	  $(MAKE) install ; if [ $$? != 0 ] ; then \
	        echo "Exit status from make was $$?" ; exit 1 ; fi ;) ;\
	  done

clean:
	@for DIR in $(OBJDIRS) ;\
	  do \
	  (cd $$DIR ;\
	  $(MAKE) clean ; if [ $$? != 0 ] ; then \
	        echo "Exit status from make was $$?" ; exit 1 ; fi ;) ;\
	  done
	-rm -f ${exec_prefix}/bin/*
	-rm -f ${exec_prefix}/lib/*.a
	-rm -f ${exec_prefix}/module/*.mod  
	-rm -f ${exec_prefix}/src/*.o   
	-rm -rf html/[a-z]*

distclean:
	-rm -rf build
	-rm Makefile
	-rm build_command
	-rm config.log
	-rm config.status
	-rm config/config.h
	-rm config/mh-config-use
	-rm config/set-up.info
	-rm -rf doc/html
	-rm -rf doc/latex
	-rm -rf html/[a-z]*

doc:
	doxygen doc/resources/doxyfile_icon_html
	@echo 
	@echo "Start of HTML documentation: doc/html/index.html"
	@echo 

pdf: 
	doxygen doc/resources/doxyfile_icon_pdf

index:
	-rm -rf html/[a-z]*
	scripts/f2html_scripts/f2html.pl -f scripts/f2html_scripts/fgenrc -d html $(SRCDIRS)

checkstyle:
	scripts/codestyle_scripts/process_src -v
