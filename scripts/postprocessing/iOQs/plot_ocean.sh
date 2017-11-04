#! /usr/bin/env bash
set -e

set -x # DEBUG

export DATADIR=$DATADIR
export WORKDIR=$WORKDIR
export EXPID=$EXPID
export CDO=${CDO:-cdo -s -P 4}
export Y1=$Y1
export Y2=$Y2
export CHUNK=${CHUNK:-1}
export GRID=$GRID
export LEV=$LEV
export POOL=$POOL
export BINDIR=$(dirname $0)

run_bg () {(
    trap 'status=$?; [ $status != 0 ] && echo $1 $status >> status' EXIT
    time "$@"
) &}

rm -f status

run_bg $BINDIR/OQs_icon_surf.sh
run_bg $BINDIR/OQs_icon_subsurf.sh
run_bg $BINDIR/OQs_icon_moc.sh

wait
cd ${EXPID}_${Y1}-${Y2}
cat *_moc*.lst_ > ${EXPID}.lst 
cat *_surf_*.lst_ >> ${EXPID}.lst 
cat *_surf39_*.lst_ >> ${EXPID}.lst 
cat *_subsurf*.lst_ >> ${EXPID}.lst 


$BINDIR/create_plot_browser -t $EXPID ${EXPID}.lst > index.html


if [ -e status ]
then
    echo "Sorry: errors during plotting: $(<status)" >&2
    exit 1
fi
