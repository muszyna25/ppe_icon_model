# Find the native NETCDF includes and library
#
#  NETCDF_INCLUDE_DIRS - where to find netcdf include, etc.
#  NETCDF_LIBRARIES    - List of libraries when using netcdf.
#  NETCDF_FOUND        - True if netcdf libraries are found.
#
#------------------------------------------------------------------------------
#
# first find the netCDF configuration information program nc-config
#

set (_netcdf_root "$ENV{NETCDF_ROOT}")
message (STATUS "NETCDF_ROOT: ${_netcdf_root}")

find_program (_netcdf_config_executable
  NAMES nc-config
  HINTS ENV NETCDF_ROOT
  PATH_SUFFIXES bin
  DOC "netcdf configuration executable"
  NO_DEFAULT_PATH
)
mark_as_advanced (_netcdf_config_executable)
message (STATUS "nc-config: ${_netcdf_config_executable}")

find_program (_netcdf_dump_executable
  NAMES ncdump
  HINTS ENV NETCDF_ROOT
  PATH_SUFFIXES bin
  DOC "netcdf dump executable"
  NO_DEFAULT_PATH
)
mark_as_advanced (_netcdf_dump_executable)
message (STATUS "ncdump   : ${_netcdf_dump_executable}")

execute_process (COMMAND ${_netcdf_config_executable} --version OUTPUT_VARIABLE _netcdf_version) 
string (REGEX MATCHALL "([.0-9]+)" _netcdf_version "${_netcdf_version}")

if (NOT _netcdf_version MATCHES "^4.*")
#  execute_process (COMMAND ${_netcdf_dump_executable OUTPUT_VARIABLE _ncdump}
  set (_netcdf_version "3")
  message (STATUS "netCDF version -> : ${_ncdump}")
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