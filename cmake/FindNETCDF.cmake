# Find the native NETCDF includes and library
#
#  NETCDF_INCLUDE_DIRS - where to find netcdf include, etc.
#  NETCDF_LIBRARIES    - List of libraries when using netcdf.
#  NETCDF_FOUND        - True if netcdf libraries are found.
#
#------------------------------------------------------------------------------
#
function (_reduce_list _list _prefix )
  set (_result)
  foreach (_member ${ARGN})
    if (_member MATCHES "(.*)${_prefix}/(.*)")
      string (REGEX REPLACE "${_prefix}" "" _member ${_member})
      list (APPEND _result ${_member})
    endif()
  endforeach()
  set (${_list} ${_result} PARENT_SCOPE)
endfunction(_reduce_list)
#
#------------------------------------------------------------------------------
#

set (_netcdf_root "$ENV{NETCDFROOT}")

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

  execute_process (COMMAND ${_netcdf_dump_executable} 
    ERROR_VARIABLE _ncdump_version)
  string (REGEX REPLACE ".*\"([.0-9]+)\".*" "\\1" _netcdf_version "${_ncdump_version}")

  # check for C header and C library

else()

  execute_process (COMMAND ${_netcdf_config_executable} --version 
    OUTPUT_VARIABLE _netcdf_version) 
  string (REGEX MATCHALL "([.0-9]+)" _netcdf_version "${_netcdf_version}")
  string (REGEX REPLACE "[0-9]+\\.([0-9]+).*" "\\1" _netcdf_minor_version "${_netcdf_version}")

endif()

# catch C first, all versions have the same storage pattern
# do we have hdf5 based netcdf support build in

set (_netcdf_nc4_support "no") 
if (_netcdf_version VERSION_GREATER 3)
  execute_process (COMMAND ${_netcdf_config_executable} --has-nc4 
    OUTPUT_VARIABLE _netcdf_nc4_support 
    OUTPUT_STRIP_TRAILING_WHITESPACE) 
endif()

find_path(NETCDF_C_INCLUDE_DIR netcdf.h 
  HINTS ${_netcdf_root}/include)

find_library(
  NETCDF_C_LIBRARY 
  NAMES netcdf
  HINTS ${_netcdf_root}/lib
)

mark_as_advanced(NETCDF_C_LIBRARY NETCDF_C_INCLUDE_DIR)

message (STATUS "netCDF C include directory       : ${NETCDF_C_INCLUDE_DIR}")
message (STATUS "netCDF C libraries               : ${NETCDF_C_LIBRARY}")

if (_netcdf_nc4_support)

  message (STATUS "netCDF installation supports nc4 format (hdf5 based)")

  # need to find the required libraries, information can only be retrieved 
  # from --cflags option of nc-config 
  
  execute_process (COMMAND ${_netcdf_config_executable} --cflags 
    OUTPUT_VARIABLE _netcdf_cflags
    OUTPUT_STRIP_TRAILING_WHITESPACE) 

  string (REGEX MATCHALL "-I([^\ ]+\ |[^\ ]+$)" _cflags_list "${_netcdf_cflags}")  
  _reduce_list (_all_includes_list "-I" ${_cflags_list})
  foreach (_include ${_all_includes_list})
    string (TOLOWER  "${_include}" _include_lc)
    if (_include_lc MATCHES "hdf")
      string (REGEX REPLACE "/include" "" _hdf5_root "${_include}")
      message (STATUS "HDF5 root directory: ${_hdf5_root}") 
    endif()
  endforeach()

endif()


# catch Fortran with different storage patterns
#
# 1. netcdf 3 has only one way we know of: it is included in the 
#    C library (interface via cfortran.h) and an include file must 
#    be available. Need to check if include file exists. If yes
#    C paths are sufficient for proper linking.
#
# 2. netcdf 4 before netcdf 4.2
#
# 3. netcdf 4.2 and newer
#

if (_netcdf_version VERSION_LESS 4.0)

  find_path(NETCDF_Fortran_INCLUDE_DIR netcdf.inc
    HINTS ${_netcdf_root}/include)

  set (NETCDF_Fortran_MODULE_DIR "n/a")

  find_library(
    NETCDF_Fortran_LIBRARY 
    NAMES netcdf
    HINTS ${_netcdf_root}/lib)

elseif (_netcdf_version VERSION_LESS 4.2)

  message (STATUS "using netcdf 4.0 search pattern")

  find_path(NETCDF_Fortran_INCLUDE_DIR netcdf.inc 
    HINTS ${_netcdf_root}/include)

  find_path(NETCDF_Fortran_MODULE_DIR netcdf.mod 
    HINTS ${_netcdf_root}/include)

  find_library(
    NETCDF_Fortran_LIBRARY 
    NAMES netcdff
    HINTS ${_netcdf_root}/lib
)

else()

  message (STATUS "using netcdf 4.2 search pattern")

  set (_netcdf_root "$ENV{NETCDFFROOT}")
  message (STATUS "NETCDFFROOT: ${_netcdf_root}")

endif()

message (STATUS "netCDF Fortran include directory : ${NETCDF_Fortran_INCLUDE_DIR}")
message (STATUS "netCDF Fortran module directory  : ${NETCDF_Fortran_MODULE_DIR}")
message (STATUS "netCDF Fortran libraries         : ${NETCDF_Fortran_LIBRARY}")

# Per-recommendation
set (NETCDF_INCLUDE_DIRS "${NETCDF_INCLUDE_DIR}")
set (NETCDF_LIBRARIES    "${NETCDF_Fortran_LIBRARY} ${NETCDF_C_LIBRARY}")

# handle the QUIETLY and REQUIRED arguments and set NETCDF_FOUND to TRUE if 
# all listed variables are TRUE
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(NETCDF DEFAULT_MSG _netcdf_version)