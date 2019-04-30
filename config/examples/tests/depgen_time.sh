#!/bin/bash

set -eu
set -o pipefail

curr_dir=$(pwd)
root_dir=$(cd "$(dirname "$0")"; pwd)
work_dir=$(basename -s '.sh' "$0")
rm -rf "${work_dir}" && mkdir -p "${work_dir}" && cd "${work_dir}"
work_dir=$(pwd)

cfg='gcc.external'
time_cmd=$(which time)

"${root_dir}/../${cfg}"
for i in {1..10}; do
  "$time_cmd" -f '%e' -o time.txt -a make -j8 depend
  rm -rf ./deps
done
echo
echo '********************'
echo
awk '{s+=$1; c++}END{print "Total runs:", c, "\nAverage wall time:", s/c, "seconds"}' time.txt

cd "${curr_dir}"
