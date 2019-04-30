#!/bin/bash

set -eu
set -o pipefail

curr_dir=$(pwd)
root_dir=$(cd "$(dirname "$0")"; pwd)
work_dir=$(basename -s '.sh' "$0")
rm -rf "${work_dir}" && mkdir -p "${work_dir}" && cd "${work_dir}"
work_dir=$(pwd)

if test $# -eq 0; then
  vendors='gcc intel nag'
else
  vendors=$@
fi
modes='bundled external'

for vendor in ${vendors}; do
  mkdir ${vendor} && cd ${vendor}
  for mode in ${modes}; do
    mkdir ${mode} && cd ${mode}
    "${root_dir}/../${vendor}.${mode}" 2>&1 | tee "${work_dir}/${vendor}.${mode}.log"
    make -j8 2>&1 | tee -a "${work_dir}/${vendor}.${mode}.log"
    cd ..
  done
  cd ..
done

for vendor in ${vendors}; do
  cd ${vendor}
  for mode in ${modes}; do
    cd ${mode}
    echo "Running make for the second time in $(pwd)..." | tee -a "${work_dir}/${vendor}.${mode}.log"
    make 2>&1 | tee -a "${work_dir}/${vendor}.${mode}.log"
    cd ..
  done
  cd ..
done

cd "${curr_dir}"
