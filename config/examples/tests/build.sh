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

case " $vendors " in
  *\ nag\ *) export NAG_KUSARI_FILE='license.mpimet.mpg.de:';;
esac

for vendor in ${vendors}; do
  vendor_tests=$(find "${root_dir}/../" -name ${vendor}'.*' -type f -executable | sort)
  for vendor_test in ${vendor_tests}; do
    vendor_test_dir=$(basename "${vendor_test}")
    mkdir "${vendor_test_dir}"
    cd "${vendor_test_dir}"
    echo "*************************"
    pwd
    echo "*************************"
    log_file="${work_dir}/$(basename "${vendor_test}").log"
    "${vendor_test}" 2>&1 | tee "${log_file}"
    make -j8 2>&1 | tee -a "${log_file}"
    cd ..
  done
done

for vendor in ${vendors}; do
  vendor_tests=$(find "${root_dir}/../" -name ${vendor}'.*' -type f -executable | sort)
  for vendor_test in ${vendor_tests}; do
    vendor_test_dir=$(basename "${vendor_test}")
    cd "${vendor_test_dir}"
    echo "*************************"
    pwd
    echo "*************************"
    log_file="${work_dir}/$(basename "${vendor_test}").log"
    echo "Running make for the second time in $(pwd)..." | tee -a "${log_file}"
    make 2>&1 | tee -a "${log_file}"
    cd ..
  done
done

cd "${curr_dir}"
