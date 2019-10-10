#!/bin/bash
set -eu
set -o pipefail

root_dir=$(cd "$(dirname "$0")"; pwd)
work_dir="${root_dir}/$(basename "$0" '.sh')"
rm -rf "${work_dir}" && mkdir -p "${work_dir}" && cd "${work_dir}"

###############################################################################

config='mpim-mpipc45-spack/gcc.no_jsbach'

time_cmd=$(which time)

eval "${root_dir}/../${config}"

###############################################################################

for i in {1..5}; do
  "$time_cmd" -f '%e' -o time.txt -a make "$@" depend
  touch icon.mk
done

###############################################################################

echo
echo '********************'
echo
awk '{s+=$1; ssq+=$1^2}END{print "Total runs:", NR, "\nAverage wall time:", s/NR, "seconds\nStandard deviation:", sqrt((ssq-s^2/NR)/(NR-1))}' time.txt
