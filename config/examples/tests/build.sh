#!/bin/bash

set -eu
set -o pipefail

root_dir=$(cd "$(dirname "$0")"; pwd)
work_dir="${root_dir}/$(basename -s '.sh' "$0")"
rm -rf "${work_dir}" && mkdir -p "${work_dir}" && cd "${work_dir}"

if test $# -eq 0; then
  vendors='gcc intel nag'
else
  vendors=$@
fi

test_suite='mpim-mpipc45-spack'

art_repo='/scratch/local1/icon-extra-externals/art'
art_dir="$root_dir/../../../src/art"
if test ! -d "$art_dir" && test -d "$art_repo"; then
  git clone "$art_repo" "$art_dir"
fi

ecrad_repo='/scratch/local1/icon-extra-externals/ecrad'
ecrad_dir="$root_dir/../../../externals/ecrad"
if test ! -d "$ecrad_dir" && test -d "$ecrad_repo"; then
  git clone "$ecrad_repo" "$ecrad_dir"
fi

for vendor in ${vendors}; do
  vendor_tests=$(find "${root_dir}/../${test_suite}" -name ${vendor}'.*' -type f -executable | sort)
  for vendor_test in ${vendor_tests}; do
    vendor_test_dir=$(basename "${vendor_test}")
    mkdir "${vendor_test_dir}"
    cd "${vendor_test_dir}"
    echo "*************************"
    pwd
    echo "*************************"
    log_file="${work_dir}/$(basename "${vendor_test}").log"
    "${vendor_test}" 2>&1 | tee "${log_file}"
    echo "Running 'make -j8' in $(pwd)" | tee -a "${log_file}"
    make -j8 2>&1 | tee -a "${log_file}"
    echo "Running make for the second time in $(pwd)..." | tee -a "${log_file}"
    make 2>&1 | tee -a "${log_file}" | tee make_2.log
    if grep '^\(  GEN      version\.c\|perl .*pvcs.pl\)' make_2.log 2>&1 >/dev/null; then :
      else echo -e "\nMake did not try to regenerate 'version.c'" && exit 1; fi
    if grep ' [ ]*\(DEPGEN\|FC\|CC\|CCLD\|FCLD\|-o icon\)' make_2.log 2>&1 >/dev/null; then
      echo -e "\nThe second call of make did some compilation/linking work" && exit 1; fi
    echo "Running make for the third time in $(pwd)..." | tee -a "${log_file}"
    make 2>&1 | tee -a "${log_file}" | tee make_3.log
    if diff make_2.log make_3.log 2>&1 >/dev/null; then :
      else echo -e "\nSecond and third calls of make gave different output" && exit 1; fi
    rm make_2.log  make_3.log
    echo "Running 'make -j8 distclean' in $(pwd)..." | tee -a "${log_file}"
    make -j8 distclean 2>&1 | tee -a "${log_file}"
    if find -mindepth 1 -print -quit 2>/dev/null | grep -q . 2>&1 >/dev/null; then
      echo -e "\n'make -j8 distclean' did not delete all files" && exit 1; fi
    cd ..
    find "${vendor_test_dir}" -type d -empty -delete
  done
done

