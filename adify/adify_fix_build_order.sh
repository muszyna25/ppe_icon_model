#!/bin/bash
PPWD=${PWD}
BASEDIR=`dirname $0`
DATADIR=${BASEDIR}/tmpdata
VRBS=0
#
tell() {
    if [ 0 -lt "${VRBS}" ] ; then
	MSG=$1
	SEP=$2
	CMD=$3
	WMSG=$4
	echo -e "$MSG"
	if [ -n "${CMD}" ] ; then
	    ${CMD}
	fi
	if [ -n "${WMSG}" ] ; then
	    read -p "${WMSG}" AA
	fi
	if [ -n "${SEP}" ] ; then
	    echo -e ${SEP}
	else
	    echo "--------------------------"
	fi
    fi
}
#
TDIR=`dirname $0`/__tmp__
mkdir -p ${TDIR}
#
echo "$0 :: start prepartion .."
make clean
make -n > ${TDIR}/00make.lst
sed -n 's/^\([^\s]*\s*.*\s\)\([a-zA-Z0-9_\.\/]*\.f90\)$/\2/p' ${TDIR}/00make.lst > ${TDIR}/00files.lst
sed -n 's/^\([^\s]*\s*.*\s\)\([a-zA-Z0-9_\.\/]*\.f90\)$/\1/p' ${TDIR}/00make.lst | sort | uniq > ${TDIR}/00compile.cmd
if [ 1 -ne `cat ${TDIR}/00compile.cmd | wc -l ` ] ; then
    echo "ERROR: compilation command not uniq :"
    cat  ${TDIR}/00compile.cmd
    exit 11
fi
#tell00 "Compilation command  :"  "cat ${TDIR}/00compile.cmd"

CMPL=`cat ${TDIR}/00compile.cmd`
rm -f ${TDIR}/01files.lst


compile() {
    FLS=$1
    FILE=$2
    CSTOP=0
#    echo -e '\n\n\n'
#    echo "FLS  :: ${FLS}"
#    echo "FILE :: ${FILE}"
    while [ 1 -ne "${CSTOP}" ] ; do
	tell " Compiling now ${FILE} .."
	OFILE=`basename ${FILE} .f90`.o
	make --no-print-directory name=${OFILE} one 1>___00.log 2>&1
	RES=`grep "Fatal Error:" ___00.log`
##	tell "   compile ${FILE} RES00 :: ${RES}" "::::::::::::::::::::"
##	tell "   compile log file : " ";;;;;;;;;;;;;;;;;;"  "cat ___00.log"
	if [ -n "${RES}" ] ; then
	    RES=13
	else
	    RES=0
	fi
	tell "   compile ${FILE} RES  :: ${RES}" "::::::::::::::::::::"

	if [ 0 -ne ${RES} ] ; then
	    MMOD=`grep "Fatal Error: Can't open module file ‘"   ___00.log`
	    if [ -n "${MMOD}" ] ; then
		MMOD=`sed -n "s/\s*Fatal\s*Error:\s*Can't\s*open\s*module\s*file\s*\‘//p"  ___00.log | sed 's/\.mod.*$//' | sed -n '1p'`
		tell " Error: Missing module : ${MMOD}"
		compile "${FLS}" ${MMOD}
	    else
		tell "Unexpected ERROR  during compilation of ${MMOD} :" "!!!!!!!!!!!!!!!!" "cat ___00.log" 
		exit 12
	    fi
	fi
	CSTOP=1
    done
    LFILES=`sed -n '/^gfortran.*\.f90\s*$/p' ___00.log`
#    tell "     Returning fromn successfull sub - compilation :: " "===================="  "echo ${LFILES}" " PRESS "
    sed -n '/^gfortran.*\.f90\s*$/p' ___00.log  >>  ${TDIR}/01files.lst
#    sed -n '/^gfortran.*\.f90\s*$/p' ___00.log
    sed -n 's,^gfortran.*-c.*/,,p'  ___00.log | sed 's/^/    Missing module dependency found ::: /'
}


echo "$0 :: start passive compilation to detect missing module dependencies .."
STOP=0
while [ 1 -ne "${STOP}" ] ; do
    make --no-print-directory 1>___00.log 2>&1
#    make --no-print-directory 2>&1 | tee ___00.log
#    tell "\n\n\n  Back from MAIN MAKE with ${RES} ... \n\n\n"
    RES=`grep "Fatal Error:" ___00.log`
##   tell "   compile ${FILE} RES00 :: ${RES}" "::::::::::::::::::::"
##   tell "   compile log file : " ";;;;;;;;;;;;;;;;;;"  "cat ___00.log"
    if [ -n "${RES}" ] ; then
	RES=13
    else
	RES=0
    fi
    if [ 0 -ne ${RES} ] ; then
	sed -n '/^gfortran.*\.f90\s*$/p' ___00.log | sed '$d'  >> ${TDIR}/01files.lst
	##sed -n '/^gfortran.*\.f90\s*$/p' ___00.log | sed '$d'
	MMOD=`grep "Fatal Error: Can't open module file ‘"   ___00.log`
	if [ -n "${MMOD}" ] ; then
	    MMOD=`sed -n "s/\s*Fatal\s*Error:\s*Can't\s*open\s*module\s*file\s*\‘//p"  ___00.log | sed 's/\.mod.*$//' | sed -n '1p'`

	    tell " Error: Missing module : ${MMOD}"
	    
	    compile "${PPWD}/01files.lst" ${MMOD}

	    ##tell "" "" ""  "                          BACK in MAIN   ${STOP} ...."
	    
	    ##exit 12
	else
	    echo "ERROR  during compilation :"
	    cat ___00.log
	    exit 12
	fi
    else
	sed -n '/^gfortran.*\.f90\s*$/p' ___00.log  >> ${TDIR}/01files.lst
	##sed -n '/^gfortran.*\.f90\s*$/p' ___00.log
	STOP=1
    fi
done

sed 's/^gfortran.*-c//'  ${TDIR}/01files.lst | uniq > ${TDIR}/02files.lst
sed 's,^\(.*/\),,;s,$, \\,;$s/\\$//' ${TDIR}/02files.lst > ${TDIR}/02files.sources
sed 's,^\(.*/\)\(.*\)\.\(.*\),\2.o: \1\2\.\3,' ${TDIR}/02files.lst > ${TDIR}/02files.deps
#
# copy stuff that has to be kept for adify_config.sh
#
mkdir -p ${DATADIR}
cp ${BASEDIR}/__tmp__/02files.sources  ${DATADIR}/sources.lst
sed 's,^[ \t\./]*,,' ${BASEDIR}/__tmp__/02files.lst > ${DATADIR}/files.lst
#
#  cleaning up
#
rm -rf ${TDIR}
