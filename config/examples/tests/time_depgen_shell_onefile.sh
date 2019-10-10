#!/bin/bash
set -eu
set -o pipefail

root_dir=$(cd "$(dirname "$0")"; pwd)
work_dir="${root_dir}/$(basename "$0" '.sh')"
rm -rf "${work_dir}" && mkdir -p "${work_dir}" && cd "${work_dir}"

###############################################################################

icon_dir="${root_dir}/../../.."

python='python'
depgen="${icon_dir}/utils/mkhelper/depgen.py"
depgen_args='--pp-enable --pp-eval-expr --fc-enable'
max_file_count=1000

time_cmd=$(which time)

set +e
files=$(find "${icon_dir}/src" -name '*.f90' | sort | head -${max_file_count})
set -e
echo "Total number of files: $(echo "$files" | wc -l)"

###############################################################################

for i in {1..5}; do
  cat $files | "${time_cmd}" -f '%e' -o time.txt -a "${python}" "${depgen}" --obj-name all.o --src-name all.f90 -o all.f90.d --src-roots="${icon_dir}" ${depgen_args} -- -Jmods -I"${icon_dir}/src/include"
done

###############################################################################

echo
echo '********************'
echo
awk '{s+=$1; ssq+=$1^2}END{print "Total runs:", NR, "\nAverage wall time:", s/NR, "seconds\nStandard deviation:", sqrt((ssq-s^2/NR)/(NR-1))}' time.txt
