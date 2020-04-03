#!/bin/sh

fn_error()
{
  echo "ERROR: $2" >&2 && exit $1
}

fn_warn()
{
  echo "WARNING: $1" >&2
}

fn_rsync()
{
  cmd="rsync -uavz $1 $2 $3"
  echo "$cmd"
  eval "$cmd" || fn_error 4 "failed to copy '$1' to '$2'"
}

fn_copy_file()
{
  cmd="$mkdir_p \"$2\" && cp \"$1\" \"$2\""
  echo "$cmd"
  eval "$cmd" || fn_error 4 "failed to copy '$1' to '$2'"
}

# Default values:
builddir=`pwd`
force=no

prev=
for option; do
  if test -n "$prev"; then
    eval $prev=\$option
    prev=
    continue
  fi

  optarg=
  case $option in
    *=?*) optarg=`expr "X$option" : '[^=]*=\(.*\)'` ;;
  esac

  case $option in
    --help | -h)
      cat <<_EOF
Usage: $0 [OPTION]...

Generates and transfers files required for running ICON from source and build
directories to PREFIX directory. Paths to source and PREFIX directories are read
from file 'config.status' generated by the configure script of ICON.

A possible usage scenario is to build ICON on a fast local disk partition of a
login node and then transfer data required for running the model to a shared
partition available on compute nodes.

Defaults for the options are specified in brackets.

Options:
  -h, --help       display this help and exit
  -f, --force      disable checks preventing unexpected results
  --builddir=DIR   look for 'config.status' in DIR [`pwd`]

_EOF
      exit 0 ;;
    --builddir)
      prev=builddir ;;
    --builddir=*)
      builddir=$optarg ;;
    --force | -f)
      force=yes ;;
    *)
      fn_error 1 "unrecognized option: '$option': try '$0 --help' for more information" ;;
  esac
done

test -z $builddir && fn_error 1 "expected a non-empty value for --builddir"
test -f "$builddir/config.status" || fn_error 2 "cannot find config.status in $builddir: try '$0 --help' for more information"
test -f "$builddir/bin/icon" || fn_error 2 "cannot find 'bin/icon' in '$builddir': you need to build first: 'cd \"$builddir\" && make'"

eval "`echo "top_srcdir='@abs_top_srcdir@'; prefix='@prefix@'; mkdir_p='@MKDIR_P@'" | "$builddir/config.status" --file=-`"

test x"$force" = xno && test x"$prefix" = 'x/usr/local' && fn_error 3 \
"prefix equals to the default value '/usr/local': rerun the configure script
specifing an installation prefix other than '/usr/local' using '--prefix' or
rerun this script with an option '--force' to disable this check"

# We do not run 'make install' in $builddir because the user might have
# specified '--exec-prefix' or '--bindir' configure options, which is valid but
# potentially leads to the executables being installed in a directory other than
# $prefix/bin. Therefore, we simply copy '$builddir/bin/icon' to '$prefix/bin'.

fn_copy_file "$builddir/bin/icon" "$prefix/bin/"

if test -f "$builddir/run/set-up.info"; then
  fn_copy_file "$builddir/run/set-up.info" "$prefix/run/"
else
  fn_warn "cannot find 'run/set-up.info' in '$builddir': trying to generate it..."
  cmd="\"$builddir/run/collect.set-up.info\" \"$prefix/run/set-up.info\""
  echo "$cmd"
  eval "$cmd" || fn_error 4 "failed to generate '$prefix/run/set-up.info'"
fi

fn_rsync "$top_srcdir/run" "$prefix/" "--exclude='*.in' --exclude='.*'"
fn_rsync "$top_srcdir/vertical_coord_tables" "$prefix/"
fn_rsync "$top_srcdir/externals" "$prefix/" "--exclude='.git' --exclude='*.f90' --exclude='*.F90' --exclude='*.c' --exclude='*.h' --exclude='*.Po' --exclude='tests' --exclude='rrtmgp*.nc' --exclude='*.mod' --exclude='*.o'"
fn_rsync "$top_srcdir/data" "$prefix/"

fn_rsync "$top_srcdir/make_runscripts" "$prefix/"

cat <<_EOF

********************************************************************************
All runtime files have been successfully transferred to '$prefix'
You can proceed with the following command:

cd "$prefix" && ./make_runscripts

_EOF
