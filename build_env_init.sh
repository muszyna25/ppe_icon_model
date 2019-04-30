# This file contains some helper shell functions, which help to initialize the
# building environment. It is supposed to be sourced inside the script passed
# as BUILD_ENV parameter of the configure script of ICON (see example below).
# It is recommended to keep this script as portable as possible:
# https://www.gnu.org/software/autoconf/manual/autoconf-2.69/html_node/Portable-Shell.html
#
# Example:
# ICON_DIR=`pwd`
# ${ICON_DIR}/configure BUILD_ENV='. ${ICON_DIR}/build_env_init.sh; switch_for_module PrgEnv-cray/6.0.4;'

function switch_for_module {
  packageName=`echo $1 | cut -d/ -f1`
  # same name:
  module unload `module -t list 2>&1 | sed -n '/'$packageName'/p'`
  # conflicts:
  module unload `module show $packageName 2>&1 | sed -n 's/^conflict[ \t][ \t]*\([^ \t][^ \t]*\)/\1/p'`
  module load $1
}
