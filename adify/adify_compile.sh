#!/bin/bash
#
#SRC=../../../src/
BASEDIR=`dirname $0`
DATADIR=${BASEDIR}/tmpdata
if [ ! \( -f ${DATADIR}/iconblddir -a -f ${DATADIR}/compile.cmd -a -f ${DATADIR}/link.cmd -a -f ${DATADIR}/files.lst \)  ] ; then
    echo -e "\n\n  Error: Run adify/adify_config.sh first."
    echo -e "         To see available options Run: bash adify/adify_config.sh --help\n\n"
    exit 3
fi
DEST=`cat ${DATADIR}/iconblddir`
#
FCCMD=`cat ${DATADIR}/compile.cmd`
LKCMD=`cat ${DATADIR}/link.cmd`
#
function show() {
    WHT=$1
    while [ -n "${WHT}" ] ; do
	echo -e "${WHT}  :: \n${!WHT}\n==================="
	shift 1
	WHT=$1
    done
}
#show BASEDIR DEST FCCMD LKCMD
#
SRC=`echo  ${DEST} | sed 's/[-a-zA-Z0-9_]*/\.\./g' `
show SRC
#
PPWD=${PWD}
#
cd $DEST
for i in `cat ${PPWD}/${DATADIR}/files.lst` ; do
    echo "  Working on $i"
    ON=`basename $i .f90`.o
    if [ ! -f ${ON} ] ; then
        ${FCCMD} ${SRC}/${i}
	RES=$?
	echo "RES ::: ${RES}"
	if [ 0 -lt ${RES} ] ; then break; fi
    fi
done
bash ${PPWD}/${DATADIR}/link.cmd
cd ${PPWD}
