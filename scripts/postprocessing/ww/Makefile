# Makefile for ww_main
#
PROGRAM = ww_main

VPATH = ../../../src/atm_phy_nwp

FC = ifort

FFDEBUG = # -debug -ftrapuv -warn -fp-stack-check -fstack-protector-all -fpe-all=0 -check arg_temp_created,format,output_conversion,uninit,pointers,stack,bounds
FFLAGS  = -I. ${GRIB_API_INCLUDE} -O3 -fpp1 -traceback -fpe0 -fp-speculation safe -fp-model source -sox $(FFDEBUG) -list -Winline -DONLYWW

LD  = $(FC)

LDFLAGS  = ${GRIB_API_LIB_FTN} # -Wl,"--print-map"

OBJS  = mo_wwonly.o mo_nwp_ww.o ww_main.o

.SUFFIXES: .f90

%.o: %.f90
	$(FC) $(FFLAGS) -c $<

$(PROGRAM): $(OBJS)
	$(LD) $(OBJS) $(LDFLAGS) -o $(PROGRAM)

mo_wwonly.o: mo_wwonly.f90

mo_nwp_ww.o: mo_nwp_ww.f90 mo_wwonly.o

ww_main.o: ww_main.f90 mo_nwp_ww.o mo_wwonly.o

clean:
	 rm $(OBJS) mo_nwp_ww.mod mo_wwonly.mod
