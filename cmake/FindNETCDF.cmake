# Find the native NETCDF includes and library
#
#  NETCDF_INCLUDE_DIRS - where to find netcdf include, etc.
#  NETCDF_LIBRARIES    - List of libraries when using netcdf.
#  NETCDF_FOUND        - True if netcdf libraries are found.
#
#------------------------------------------------------------------------------
#

set (_netcdf_root "$ENV{NETCDFROOT}")
message (STATUS "NETCDFROOT: ${_netcdf_root}")

# first find the netCDF configuration information program nc-config

find_program (_netcdf_config_executable
  NAMES nc-config
  HINTS ENV NETCDFROOT
  PATH_SUFFIXES bin
  DOC "netcdf configuration executable"
  NO_DEFAULT_PATH
)
mark_as_advanced (_netcdf_config_executable)

# and/or ncdump to catch the version of old netcdf installations

find_program (_netcdf_dump_executable
  NAMES ncdump
  HINTS ENV NETCDFROOT
  PATH_SUFFIXES bin
  DOC "netcdf dump executable"
  NO_DEFAULT_PATH
)
mark_as_advanced (_netcdf_dump_executable)

if (_netcdf_config_executable MATCHES "NOTFOUND")

  # this is netcdf version 3, get full version number first

  execute_process (COMMAND ${_netcdf_dump_executable} ERROR_VARIABLE _ncdump_version)
  string (REGEX REPLACE ".*\"([.0-9]+)\".*" "\\1" _netcdf_version "${_ncdump_version}")

  # check for C header and C library

else()

  execute_process (COMMAND ${_netcdf_config_executable} --version OUTPUT_VARIABLE _netcdf_version) 
  string (REGEX MATCHALL "([.0-9]+)" _netcdf_version "${_netcdf_version}")
  string (REGEX REPLACE ".*\\.([0-9]+).*" "\\1" _netcdf_minor_version "${_netcdf_version}")
  message (STATUS "netCDF minor version: ${_netcdf_vminor_ersion}")

endif()
message (STATUS "netCDF version: ${_netcdf_version}")

find_path(NETCDF_INCLUDE_DIR netcdf.h 
  HINTS ${_netcdf_root}/include)

find_library(
  NETCDF_Fortran_LIBRARY 
  NAMES netcdff
  HINTS ${_netcdf_root}/lib
)

find_library(
  NETCDF_C_LIBRARY 
  NAMES netcdf
  HINTS ${_netcdf_root}/lib
)

mark_as_advanced(NETCDF_LIBRARY NETCDF_INCLUDE_DIR)

message (STATUS "netCDF include directory: ${NETCDF_INCLUDE_DIR}")
message (STATUS "netCDF libraries        : ${NETCDF_Fortran_LIBRARY} ${NETCDF_C_LIBRARY}")

# Per-recommendation
set (NETCDF_INCLUDE_DIRS "${NETCDF_INCLUDE_DIR}")
set (NETCDF_LIBRARIES    "${NETCDF_Fortran_LIBRARY} ${NETCDF_C_LIBRARY}")

# handle the QUIETLY and REQUIRED arguments and set NETCDF_FOUND to TRUE if 
# all listed variables are TRUE
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(NETCDF DEFAULT_MSG _netcdf_version)