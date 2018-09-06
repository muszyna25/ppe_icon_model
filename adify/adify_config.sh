#!/bin/bash
#
#Define your AD system behaviour here:
#
#  set path to your ECCODES installation
#
#ECCODES=${HOME}/local/gf6-eccodes
ECCODES=/usr
#
#  specify other configure options here
#
CFGOPTS="--disable-atmo --disable-psrad --disable-jsbach --disable-serialization --disable-testbed --without-yac"
#
#CFGPREFIX="FC=gfortran-6 CC=gcc-6 "
#
# directory with AD includes (relative to repo base)
ADINCS=adify/include
#
##
##
##
##
##
##########################################################################
##                                                                      ##
##  DO **NOT** CHANGE ANYTHING BELOW THIS LINE  !!!!!!!!!!!!!!!!!!!     ##
##                                                                      ## 
##########################################################################
##
#MODELIST="A1S T1S PSV CPTEST ADJ_IO_TEST"
MODELIST="PSV CPTEST ADJ_IO_TEST"
#
##
BASEDIR=`dirname $0`
DATADIR=${BASEDIR}/tmpdata
SEDDIR=${BASEDIR}/sed
FIXBLDORDR="${BASEDIR}/adify_fix_build_order.sh"
UDEPDIR=config/usedep
TDIR=${BASEDIR}/__tmp__
#
function usage() {
    ML=`echo ${MODELIST} | tr ' ' '|'`
    echo "  Usage :: $0 [-s] [-d] [-w] [-r] [-h] [--usedep_clean] {${ML}}"
    echo "    Options :"
    echo "        -d    debug mode in FFLAGS and LDFLAGS     ( --DEBUG  )"
    echo "        -s    serial mode ( -DNOMPI )              ( --SERIAL )"
    echo "        -w    suppress warning messages            ( --NOWARN )" 
    echo "        -r    rebuild dependency info              ( --REFIX )"
    echo "        -n    do not clean before configuration    ( --NOCLN )"
    echo ""
    echo "        --usedep_clean  run 'make clean' in USEDEP directory."
    echo "                         'ragel'  needs to be installed.."
    echo ""
    echo "        -h    print help message "
    exit
}
function help() {
    usage
}
#
NOCLN=0
DEPRFX=0
NOWARN=0
NOMPI=1
DBG=0
UDEPCLN=0
MODE=""
#
UNKWON=0
while [ "${1#-}" != "$1" ] ; do
    OPT=`echo $1 | tr 'a-z' 'A-Z'`
    case ${OPT} in
	-D | --DEBUG )
	    DBG=1
	    ;;
	-S | --SERIAL )
	    NOMPI=1
	    ;;
	-W | --NOWARN )
	    NOWARN=1
	    ;;
	-R | --REFIX )
	    DEPRFX=1
	    ;;
	-N | --NOCLN )
	    NOCLN=1
	    ;;
	--USEDEP_CLEAN )
	    UDEPCLN=1
	    ;;
	-H | --HELP )
	    help
	    ;;
	*)
	    echo "Error: Unknown option : $1"
	    UNKWON=1
	    ;;
    esac
    shift
done
if [ 1 == "${UNKWON}" ] ; then exit; fi
#
MODE=$1  
if [ -z "${MODE}" ] ; then
    echo "Error :: ADmode missing."
    usage
    exit
fi
#
#
MODE=`echo ${MODE} | tr 'a-z' 'A-Z'`
UNKNWN=1
for i in ${MODELIST} ; do
    if [ "$i" == "${MODE}" ] ; then UNKNWN=0; fi
done
if [ 0 -lt "${UNKNWN}" ] ; then
    echo "Error: Unknown mode :: ${MODE}"
    usage 
fi
##
##
##
case ${MODE} in
    "PSV" )
	MODEDEF=" -D__COMPAD_${MODE}__ "
	;;
    "A1S" | "T1S" | "CPTEST" |  "ADJ_IO_TEST")
	MODEDEF=" -D__COMPAD_${MODE}__ -D__COMPAD_DECLARATIONS__ "
	;;
    *)
	echo "Error: Unknown mode :: ${MODE}"
	usage
esac
#
#  AD mode dependend special settings
#
case ${MODE} in
    "A1S" )
	MODEDEF="${MODEDEF} -D__COMPAD_ADJLOOP__ -D__NO_SEAICE_ADJOINTS__ -I../../../dco_fortran/gfortran/a1s/include -I../../../${ADINCS} "
	;;
    "T1S" )
	MODEDEF="${MODEDEF} -I../../../dco_fortran/gfortran/t1s/include "
	;;
    "CPTEST" )
	MODEDEF="${MODEDEF} -D__COMPAD_ADJLOOP__-D__NO_SEAICE_ADJOINTS__ -I../../../${ADINCS} "
##	MODEDEF="${MODEDEF}"
	;;
    "ADJ_IO_TEST" )
	MODEDEF="${MODEDEF} -D__COMPAD_ADJLOOP__-D__NO_SEAICE_ADJOINTS__ -I../../../${ADINCS} "
##	MODEDEF="${MODEDEF}"
	;;
esac
##
##
##
#
if [ -n "${NOMPI}" ] ; then
    if [ ${NOMPI} -gt 0 ] ; then
	CFGOPTS="${CFGOPTS}  --without-mpi"
    fi
else
    NOMPI=0
fi
#
mkdir -p ${DATADIR}
#
#########################################
##
##  Run original icon configure
##
#########################################
#
#./configure --enable-eccodes=${HOME}/local/gf6-eccodes --disable-atmo --disable-jsbach --disable-serialization --disable-testbed --without-yac
#./configure --enable-eccodes=/usr --disable-atmo --disable-jsbach --disable-serialization --disable-testbed --without-yac --without-mpi
#
for i in ${CFGPREFIX} ; do echo "export $i" ; export $i ; done 
echo "./configure --enable-eccodes=${ECCODES} ${CFGOPTS}"
./configure --enable-eccodes=${ECCODES} ${CFGOPTS}
for i in ${CFGPREFIX} ; do V=`echo $i | sed 's/=.*$//' `; unset $V ; done 
#
#
#
#########################################
##
##  Patch generated Makefile in main dir
##    Part 1
##
#########################################
#
#  patch generated Makefile in  main directory
#
if [ -f Makefile.generated_by_config ] ; then
    if [ config.log -nt Makefile.generated_by_config ] ; then
	cp Makefile Makefile.generated_by_config
    fi
else
    cp Makefile Makefile.generated_by_config
fi
#
#  remove not wanted libraries
#
sed -i 's/\(LAPACK_LIB[ \t]*= \).*$/\1/' Makefile
sed -i 's/ -ltixi //' Makefile
sed -i 's/\$(LAPACK_LIB) //' Makefile
###sed -i '/^SRCDIRS[ \t]*=/{s,externals/tixi/src ,,;s/lapack //;s/blas //;}' Makefile
sed -i '/^SRCDIRS[ \t]*=/{s/[^ \t]*tixi[^ \t]* //;s/[^ \t]*lapack[^ \t]* //;s/[^ \t]*blas[^ \t]* //;s/[^ \t]*yac[^ \t]* //;}' Makefile
sed -i '/^OBJDIRS[ \t]*=/{s/[^ \t]*tixi[^ \t]*// ;s/[^ \t]*lapack[^ \t]* //;s/[^ \t]*blas[^ \t]* //;s/[^ \t]*yac[^ \t]* //;}' Makefile
sed -i 's/ -DHAVE_FC_ATTRIBUTE_CONTIGUOUS//'  Makefile
#
# create new PREParation ICONBLD variable and target
#
PAT=`grep 'PREPOBJDIRS' Makefile`
if [ -z "${PAT}" ] ; then
    sed -i '/OBJDIRS\s*=/{h;s/^/PREP/;s,\s*\([-_a-Z0-9/]*/src\)\s*$,\nICONBLDDIR = \1,;H;g;}' Makefile
#    sed -i '/^\s*all:/,/^\s*done/{H;/done/{g;p;s/all:/PREP:/;P;s/\(OBJDIR\)/PREP\1/;p}}' Makefile
    sed -in '/^\s*all:/,/^\s*done/{H;/done/{p;g;s/all:/PREP:/;P;s/\(OBJDIR\)/PREP\1/;}}' Makefile

fi
#
# create new clean and o(irignal)clean
#
PAT=`grep 'oclean:' Makefile`
if [ -z "${PAT}" ] ; then
    sed -i 's/clean:/oclean:/' Makefile
    sed -i '/^oclean:/{
    	i \
clean:\
	@cd ${ICONBLDDIR} && rm -f *.o libicon.a ; \
	@for i in `find src -name "*.f90"` ; do n=build/x86_64-unknown-linux-gnu/module/`basename $$i .f90`.mod; if [ -f $$n ] ; then rm -f $$n; fi; done \
	@cd ${ICONBLDDIR}/../bin && rm -f icon test_insert_dimension ; \

    	a \
	rm -f Makefile.generated_by_config
	}' Makefile
fi
#
# remove warnings if requested
#

if [ -n "${NOWARN}" ] ; then
    if [ ${NOWARN} -gt 0 ] ; then
	sed -i '/FFLAGS\s*=/s/\s-W[-a-Z_0-9]*//g' Makefile
    fi
fi
#
#
#########################################
##
##  Extract ICON build directory
##
#########################################
#
ICONBLDDIR=`sed -n '/^\s*ICONBLDDIR\s*=\s*/{s/^\s*ICONBLDDIR\s*=\s*//;s/\s*$//;p}'  Makefile`
echo "ICONBLDDIR  :: ${ICONBLDDIR}"
echo ${ICONBLDDIR} > ${DATADIR}/iconblddir
#
#########################################
##
##  Build all supporting libraries
##
#########################################
#
PPWD=$PWD
# remove mo_cdi.o to trigger AD mode dependent recompilation
rm -f ${ICONBLDDIR}/../support/mo_cdi.o
# finsh preparation
make PREP
cd ${ICONBLDDIR}
../../../config/pvcs.pl --srcdir ../../..
cd $PPWD
#
#########################################
##
##  Patch generated Makefile in main dir
##    Part 2
##
#########################################
#
#  transform compiler and linker flags

#
sed -i "/FFLAGS\s*=/{s/\s*-D__A1S__//; s/\s*-D__T1S__//; s/\s*-D__PSV__//; s/\s*-D__CPTEST__//; s/\s*-D__ADJ_IO_TEST__//;s:$: ${MODEDEF}:;}" Makefile
sed -i "/FFLAGS\s*=/{s/\s*-O[0-4]*\s/ /;s/\s-g\s/ /;}" Makefile
sed -i "/LDFLAGS\s*=/{s/\s*-O[0-4]*\s/ /;s/\s-g\s/ /;}" Makefile
if [ 0 -lt "${DBG}" ] ; then
	    echo "debug mode"
    sed -i "/FFLAGS\s*=/{s/$/ -O0 -g /;}" Makefile
    sed -i "/LDFLAGS\s*=/{s/$/ -O0 -g /;}" Makefile
else
    sed -i "/FFLAGS\s*=/{s/$/ -O2 /;}" Makefile
    sed -i "/LDFLAGS\s*=/{s/$/ -O3 /;}" Makefile
fi
#
#########################################
##
##  build ocean dependency
##
#########################################
#
WD=${PWD}
cd ${UDEPDIR}
EDEFS=""
if [ ${NOMPI} -gt 0 ] ; then
    EDEFS="${EDEFS} NOMPI"
fi
if [ ${UDEPCLN} -gt 0 ] ; then
    make EXTRADEF="${EDEFS}" clean
fi
make EXTRADEF="${EDEFS}" ocean.dep ocean.files
cd ${WD}
unset WD
#
#########################################
##
##  patch generated Makefile in build directory 
##
#########################################
#
DIR=${ICONBLDDIR}
if [ -f ${DIR}/Makefile.generated_by_config ] ; then
    if [ config.log -nt ${DIR}/Makefile.generated_by_config ] ; then
	cp ${DIR}/Makefile ${DIR}/Makefile.generated_by_config
    fi
else
    cp ${DIR}/Makefile ${DIR}/Makefile.generated_by_config
fi
#
CLN=`grep -A2 "clean:" ${DIR}/Makefile`
sed -f ${SEDDIR}/adify_config.sed -i ${DIR}/Makefile
echo -e "${CLN}\n\trm -f Makefile.generated_by_config\n" >> ${DIR}/Makefile
#
# append dependency information to Makefile
#  remove possible error messages in ocean.dep (no idea why they are sent to stdio ...)
#
cat ${UDEPDIR}/ocean.dep | sed '/^\s*Error\s*(.*):/d;/^\s*Error\s*somewhere\s*around/,+1d' >> ${DIR}/Makefile
#  cat ${UDEPDIR}/ocean.dep >> ${DIR}/Makefile
#
#
# replace OCEFILES_DUMMY string with content of ${OCEFILES} in Makefile
#   remove possible error messages in ocean.files (no idea why they are sent to stdio ...)
#
sed -n '1,/OCEFILES_DUMMY/{/OCEFILES_DUMMY/d;p}' ${DIR}/Makefile > ${DIR}/Makefile.p1
#sed    's,^.*/,,;$!s/$/ \\/;' ${UDEPDIR}/ocean.files > ${DIR}/Makefile.p2 
sed    '/^\s*Error\s*(.*):/d;/^\s*Error\s*somewhere\s*around/,+1d;s,^.*/,,;$!s/$/ \\/;' ${UDEPDIR}/ocean.files > ${DIR}/Makefile.p2 
sed -n '/OCEFILES_DUMMY/,${/OCEFILES_DUMMY/d;p}' ${DIR}/Makefile > ${DIR}/Makefile.p3
#
cat ${DIR}/Makefile.p1 >  ${DIR}/Makefile
cat ${DIR}/Makefile.p2 >>  ${DIR}/Makefile
cat ${DIR}/Makefile.p3 >>  ${DIR}/Makefile
#
#
#
sed -i 's/^\(\t$(AR)\)\s*\($(ARFLAGS)\)\s*\($@\) $(OBJS)/\1 \2 \3 *.o ; \1 -d \3 version.o/' ${DIR}/Makefile
#
#########################################
##
##  clean up
##
#########################################
#
DIR=${ICONBLDDIR}
rm ${DIR}/Makefile.p1 ${DIR}/Makefile.p2 ${DIR}/Makefile.p3
#
#
#########################################
##
##  build and use source list with fixed depency
##
#########################################
#
if [ -n "${DEPRFX}" ] ; then
    if [ ${DEPRFX} -gt 0 ] ; then
	rm -f ${DATADIR}/sources.lst ${DATADIR}/files.lst
    fi
fi
#
if [ ! -f "${DATADIR}/sources.lst" -o ! -f  "${DATADIR}/files.lst" ] ; then
   echo "$0 :: create adify_sources.lst and adify_files.lst"
   echo "$0 :: calling now ${FIXBLDORDR}  .. this will take some time .. "
   bash ${FIXBLDORDR}
else
   echo "$0 :: Reusing existing adify_sources.lst"
fi 
#
# replace source file list in Makefile of build directory with place holder for new file list
#
if [  -f "${DATADIR}/sources.lst" ] ; then
   echo "$0 :: correcting build order in  ${ICONBLDDIR}/Makefile"
   sed -i '/\s*SRCS\s*=\s*\\s*$/,/^\s*$/{//!d;s/^\(\s*SRCS\s*=\).*$/\1 \\/};/\s*SRCS\s*=\s*\\s*$/r '"${BASEDIR}"'/adify_sources.lst'  ${ICONBLDDIR}/Makefile
   ## sed -i '/\s*SRCS\s*=\s*\\s*$/,/^\s*$/{//!d;s/^\(\s*SRCS\s*=\).*$/\1 \\/};/\s*SRCS\s*=\s*\\s*$/r ADIFY_OES/__tmp__/02files.sources'  ${ICONBLDDIR}/Makefile
else
    echo ""
    echo ""
    echo "$0 :: Warning: Cannot correct build order in  ${ICONBLDDIR}/Makefile : File ${DATADIR}/sources.lst does not exist."
    echo ""
    echo ""
fi
#
# create pure file list (without tailing backslash)
#
sed 's/\\//' ${DATADIR}/sources.lst >   ${DATADIR}/fnames.lst
#
#
#########################################
##
##  add fixed build order dependencies
##
#########################################
#
sed -i '/#*\s*[Ff]ixed\s*build\s*order\s*dependencies/,/#*\s*[Ee][Nn][Dd]\s*[Ff]ixed\s*build\s*order\s*dependencies/d' ${ICONBLDDIR}/Makefile
/bin/cat <<EOF1 >> ${ICONBLDDIR}/Makefile
## Fixed build order dependencies
##
EOF1
ll=''; for l in `cat ${DATADIR}/fnames.lst ` ; do l=`basename $l .f90`; if [ -n "${ll}" ]  ; then echo "${l}.o : ${ll}.o" >> ${ICONBLDDIR}/Makefile ; fi ; ll=$l ; done

/bin/cat <<EOF2 >> ${ICONBLDDIR}/Makefile
##
## END Fixed build order dependencies
#
EOF2

#
#
#########################################
##
##  add AD dependencies
##
#########################################
#
sed -i '/#*\s*[Aa][Dd]\s*include\s*dependencies/,/#*\s*[Ee][Nn][Dd]\s*[Aa][Dd]\s*include\s*dependencies/d' ${ICONBLDDIR}/Makefile
/bin/cat <<EOF3 >> ${ICONBLDDIR}/Makefile
## AD include dependencies
##
`grep -r 'adify_' src | sed -n '/:\s*#\s*include\s*/{s/"//g;s/\.f90/\.o/g;s,src.*/,,;s,#\s*include\s*,  ../../../adify/include/,g;p}' `
##
## End AD include dependencies
EOF3

#
#
#########################################
##
##  extract compilation command
##
#########################################
#
#
mkdir -p ${TDIR}
mkdir -p ${DATADIR}
#
echo "$0 :: start prepartion .."
if [ -n "${NOCLN}" ] ; then
    if [ ${NOCLN} -lt 1 ] ; then
	make clean
    fi
fi
make --no-print-directory -n > ${TDIR}/00make.lst
sed -n 's/^\([^\s]*\s*.*\s\)\([a-zA-Z0-9_\.\/]*\.f90\)$/\2/p' ${TDIR}/00make.lst > ${TDIR}/00files.lst
sed -n 's/^\([^\s]*\s*.*\s\)\([a-zA-Z0-9_\.\/]*\.f90\)$/\1/p' ${TDIR}/00make.lst | sort | uniq > ${TDIR}/00compile.cmd
if [ 1 -ne `cat ${TDIR}/00compile.cmd | wc -l ` ] ; then
    echo "ERROR: compilation command not uniq :"
    cat  ${TDIR}/00compile.cmd
    exit 11
fi
cp ${TDIR}/00compile.cmd ${DATADIR}/compile.cmd
#
sed -n '/gcc-ar /,$p'  ${TDIR}/00make.lst > ${TDIR}/00link.cmd
cp ${TDIR}/00link.cmd ${DATADIR}/link.cmd
#
#
#
rm -rf ${TDIR}
