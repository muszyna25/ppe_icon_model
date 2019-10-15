#!/bin/bash
set -eu
set -o pipefail

root_dir=$(cd "$(dirname "$0")"; pwd)
work_dir="${root_dir}/$(basename "$0" '.sh')"
rm -rf "${work_dir}" && mkdir -p "${work_dir}" && cd "${work_dir}"

###############################################################################

config='mpim-mpipc45-spack/gcc.external'

time_cmd=$(which time)

eval "${root_dir}/../${config}"
make -j8 depend

###############################################################################

for i in {1..5}; do
  make mostlyclean
  "${time_cmd}" -f '%e' -o time.txt -a make -j8
  "${time_cmd}" -f '%e' -o time_remake.txt -a make -j8
done

###############################################################################

echo
echo '********************'
echo
echo 'Build time:'
awk '{s+=$1; ssq+=$1^2}END{print "Total runs:", NR, "\nAverage wall time:", s/NR, "seconds\nStandard deviation:", sqrt((ssq-s^2/NR)/(NR-1))}' time.txt
echo
echo 'Remake time:'
awk '{s+=$1; ssq+=$1^2}END{print "Total runs:", NR, "\nAverage wall time:", s/NR, "seconds\nStandard deviation:", sqrt((ssq-s^2/NR)/(NR-1))}' time_remake.txt
