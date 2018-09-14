#!/bin/bash
PNAM="OCE-TimeStep-"
PBAS="CPTEST_"
PBGN="${PBAS}BEGIN"
PEND="${PBAS}END"
PDEL=" mo_name_list_output::write_name_list_output mo_ocean_output:output_ocean: "
#
DIFF=0
if [ "-d" == "$1" -o "-D" == "$1" ] ; then DIFF=1; shift; fi
F=$1
shift
if [ "-d" == "$1" -o "-D" == "$1" ] ; then DIFF=1; shift; fi
DPRG=${1:-diff}
shift
if [ "-d" == "$1" -o "-D" == "$1" ] ; then DIFF=1; shift; fi
#
if [ -z "${F}" ] ; then
    echo "Usage:  $0 [-d] logfile [diffprog] "
    exit
fi
if [ ! -f ${F} ] ; then
    echo "ERROR:  logfile '$F' not found."
    exit
fi
#
D=CPTS_`echo ${F} | sed 's/\..*$//'`
ANS="Y"
if [ -d ${D} ] ; then
    read -p "WARNING: Directory ${D} exists. Remove ?? (Y/N) : "  ANS
    if [ "y" == "${ANS}" -o "Y" == "${ANS}" ] ; then
	rm -rf ${D}
    else
	echo "ERROR: Directory ${D} exists."
	exit
    fi
fi
echo "Create directory ${D}"
mkdir -p ${D}
MTS=72
##sed -n '/CPTEST_BEGIN OCE-TimeStep-/,/CPTEST_END OCE-TimeStep-/p' $F
for i in `sed -n "s/^[ \t]*${PBGN}[ \t]*\(${PNAM}[0-9-]*\)/\1/p" $F` ; do
    sed -n "/^[ \t]*${PBGN}[ \t]*$i\$/,/^[ \t]*${PEND}[ \t]*$i\$/{/${PBAS}/d;p}" $F > $D/$i
    TS=`echo $i | sed 's/^.*Step-//;s/-.*$//;'`
    if [ "${MTS}" -lt "${TS}" ] ; then MTS=${TS}; echo "new max :: ${MTS}"; fi
done
##
for i in ${PDEL} ; do
    sed -i "/$i/d" $D/*
done
#
for i in `seq ${MTS}` ; do
    FN1=$D/${PNAM}${i}-1
    if [ -f ${FN1} ] ; then
	echo "Checking ${FN1} ..."
	let MTS1=${i}+1
	for j in `seq 2 ${MTS1}` ; do
	    FN2=$D/${PNAM}${i}-${j}
	    if [ -f ${FN2} ] ; then
		if [ 1 == "${DIFF}" ] ; then
		    cmp -s $FN1  $FN2
		    if [ $? -gt 0 ] ; then
			echo -e "-------  $FN1  $FN2  -------------"
			${DPRG}  $FN1  $FN2
			echo -e "===================================\n"
		    fi
		else
		    cmp $FN1  $FN2
		fi
	    fi
	done
    fi
done
