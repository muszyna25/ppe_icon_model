#!/bin/bash

set -eu
set -o pipefail

root_dir=$(cd "$(dirname "$0")"; pwd)
work_dir="${root_dir}/$(basename "$0" '.sh')"
rm -rf "${work_dir}" && mkdir -p "${work_dir}" && cd "${work_dir}"

if test $# -eq 0; then
  vendors='gcc intel nag cray'
else
  vendors=$@
fi

skip_tests=no
test_suite=
xfail_test_names=
failed_test_names=

case $(uname -s) in
  Darwin)
    test_suite='darwin-macports' ;;
  *)
    case $(host $(uname -n)) in
      mpipc45.mpimet.mpg.de\ *)
        test_suite='mpim-mpipc45-spack'
        xfail_test_names='gcc.cc_nagfor nag.cc_nagfor nag.disable_rpath' ;;
      mlogin*.hpc.dkrz.de\ *)
        test_suite='dkrz-mistral' ;;
      daint*.login.cscs.ch\ *)
        test_suite='cscs-daint'
        xfail_test_names='cray.bundled cray.bundled_sct cray.bundled_serialize' ;;
      xce*.dwd.de\ *)
        test_suite='dwd-xce'
        xfail_test_names='cray.bundled_dynamic cray.bundled_static' ;;
      *)
        echo "Unknown test environment" && exit 1 ;;
    esac ;;
esac

art_repo='git@gitlab.dkrz.de:m300488/art.git'
art_dir="${root_dir}/../../../externals/art"
if test ! -d "${art_dir}"; then
  git clone "${art_repo}" "${art_dir}" || echo "WARNING: failed to clone ICON-ART repository."
elif test -d "${art_dir}/.git"; then
  cd "${art_dir}" && git pull && cd - >/dev/null 2>&1 || echo "WARNING: failed to update ICON-ART repository."
else
  echo "The directory "${art_dir}" exists but it is not an ICON-ART repo." && exit 1
fi

for vendor in ${vendors}; do
  vendor_tests=$(find "${root_dir}/../${test_suite}" -name ${vendor}'.*' -type f | sort)

  for vendor_test in ${vendor_tests}; do
    test -x "${vendor_test}" || continue
    vendor_test_name=$(basename "${vendor_test}")
    mkdir "${vendor_test_name}"
    cd "${vendor_test_name}"
    echo "*************************"
    pwd
    echo "*************************"
    log_file="${work_dir}/$(basename "${vendor_test}").log"

    "${vendor_test}" 2>&1 | tee "${log_file}"

    echo "Running 'make -j8' in $(pwd)" | tee -a "${log_file}"
    make -j8 2>&1 | tee -a "${log_file}"
    if ldd ./icon | grep -q 'not found' >/dev/null 2>&1; then
      echo -e "\nThe executable of ICON ('./icon') does not contain all required RPATHS" && exit 1; fi

    echo "Running 'make' for the second time in $(pwd)..." | tee -a "${log_file}"
    make 2>&1 | tee -a "${log_file}" | tee make_2.log
    if grep '^\(  GEN      version\.c\|perl .*pvcs.pl\)' make_2.log >/dev/null 2>&1; then :
      else echo -e "\nMake did not try to regenerate 'version.c'" && exit 1; fi
    if grep ' [ ]*\(DEPGEN\|FC\|CC\|CCLD\|FCLD\|-o icon\)' make_2.log >/dev/null 2>&1; then
      echo -e "\nThe second call of make did some compilation/linking work" && exit 1; fi

    echo "Running 'make' for the third time in $(pwd)..." | tee -a "${log_file}"
    make 2>&1 | tee -a "${log_file}" | tee make_3.log
    if diff make_2.log make_3.log >/dev/null 2>&1; then :
      else echo -e "\nSecond and third calls of make gave different output" && exit 1; fi
    rm make_2.log  make_3.log

    if test "x${skip_tests}" = xno; then
      case " ${xfail_test_names} " in
        *" ${vendor_test_name} "*)
          echo "Running 'make -j8 check' in $(pwd) (which is expected to fail)..." | tee -a "${log_file}"
          if make -j8 check 2>&1 | tee -a "${log_file}"; then
            echo -e "\n'make -j8 check' is expected to fail but did not" && exit 1
          else
            if test x"${failed_test_names}" = x; then
              failed_test_names=${vendor_test_name}
            else
              failed_test_names=" ${vendor_test_name}"
            fi
          fi ;;
        *)
          echo "Running 'make -j8 check' in $(pwd)..." | tee -a "${log_file}"
          make -j8 check 2>&1 | tee -a "${log_file}"
      esac
    else echo "Skipping 'make check' in $(pwd)..."; fi

    echo "Running 'make -j8 distclean' in $(pwd)..." | tee -a "${log_file}"
    make -j8 distclean 2>&1 | tee -a "${log_file}"
    if find . -mindepth 1 -print -quit 2>/dev/null | grep -q . >/dev/null 2>&1; then
      echo -e "\n'make -j8 distclean' did not delete all files" && exit 1; fi
    cd ..
    find "${vendor_test_name}" -type d -empty -delete
  done
done

if test "x${skip_tests}" = xno; then
  if test "x${failed_test_names}" = x; then
    echo -e "\nAll tests passed."
  else
    echo -e "\nThe building tests passed, 'make check' failed for \"${failed_test_names}\" as expected."
  fi
else
  echo -e "\nAll tests passed ('make check' skipped)."
fi

