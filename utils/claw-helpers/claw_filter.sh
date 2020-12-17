#!/bin/bash

###################################################################
# Print the help for the CLAW Filter
# Arguments:
#   $1: program name
###################################################################
function claw_filter::print_help() {
  cat <<EOF
usage: $1 <OPTIONS> <INPUTFILE> ...

Checks whether Fortran files require CLAW code transformation.

Names of the files that require the code transformation are printed to the
standard output stream. If none of the input files require the transformation,
the program exits with a non-zero code.

CLAW Filter options:
   -h,--help                  : print usage.
   -f,--force                 : print all input filenames and exit.
EOF
}

###################################################################
# Process input parameters and set variables accordingly
# Global:
#   f_files, force_translation
###################################################################
function claw_filter::set_parameters() {
  while [[ -n "$1" ]]; do
    case "$1" in
    *.f90 | *.f | *.F90 | *.F)
      f_files+=("$1")
      ;;
    -h | --help)
      claw_filter::print_help "$(basename "$0")"
      exit 0
      ;;
    -f | --force) force_translation=true ;;
    esac
    shift
  done
}

force_translation=false
f_files_transformation=()
f_files=()

claw_filter::set_parameters "${@+"$@"}"

readonly force_translation

for input_file in "${f_files[@]}"; do
  if ! [[ ${force_translation} == true ]]; then
    num_directives="$(grep --count --ignore-case "!\$claw" "${input_file}")"
    num_omp_compile_guard="$(grep --count --ignore-case "!\$omp claw-guard" "${input_file}")"
    num_acc_compile_guard="$(grep --count --ignore-case "!\$acc claw-guard" "${input_file}")"
    if [[ "${num_directives}" == "0" ]] &&
      [[ "${num_omp_compile_guard}" == "0" ]] &&
      [[ "${num_acc_compile_guard}" == "0" ]]; then
      : # The file does no contain $claw.
    else
      f_files_transformation+=("${input_file}")
    fi
  else
    f_files_transformation+=("${input_file}")
  fi
done

[[ ${#f_files_transformation[@]} -eq 0 ]] && exit 1

for input_file in "${f_files_transformation[@]}"; do
  echo "${input_file}"
done
