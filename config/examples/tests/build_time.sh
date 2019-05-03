#!/bin/bash

set -eu
set -o pipefail

curr_dir=$(pwd)
root_dir=$(cd "$(dirname "$0")"; pwd)
work_dir=$(basename -s '.sh' "$0")
rm -rf "${work_dir}" && mkdir -p "${work_dir}" && cd "${work_dir}"
work_dir=$(pwd)

cfg='mpim-mpipc45-spack/gcc.external'
time_cmd=$(which time)

"${root_dir}/../${cfg}"
make -j8 depend
for i in {1..5}; do
  make mostlyclean
  "$time_cmd" -f '%e' -o time.txt -a make -j8
  "$time_cmd" -f '%e' -o time_remake.txt -a make -j8
done
echo
echo '********************'
echo
echo 'Build time:'
awk '{s+=$1; ssq+=$1^2}END{print "Total runs:", NR, "\nAverage wall time:", s/NR, "seconds\nStandard deviation:", sqrt((ssq-s^2/NR)/(NR-1))}' time.txt
echo
echo 'Remake time:'
awk '{s+=$1; ssq+=$1^2}END{print "Total runs:", NR, "\nAverage wall time:", s/NR, "seconds\nStandard deviation:", sqrt((ssq-s^2/NR)/(NR-1))}' time_remake.txt

cd "${curr_dir}"
