# Find HDF5 includes and library
#
#  HDF5_INCLUDE_DIRS - where to find netcdf include, etc.
#  HDF5_LIBRARIES    - List of libraries when using netcdf.
#  HDF5_FOUND        - True if netcdf libraries are found.
#
#------------------------------------------------------------------------------
#
message (STATUS "Looking for hdf5 library")

if (NETCDF_REQUIRES_HDF5)
  set (_hdf5_root "${NETCDF_REQUIRED_HDF5_ROOT}")
else ()
  set (_hdf5_root "$ENV{HDF5ROOT}")
endif ()

find_program (_hdf5_c_compiler_executable
  NAMES h5cc
  HINTS ${_hdf5_root}
  PATH_SUFFIXES bin
  NO_DEFAULT_PATH)

execute_process (COMMAND ${_hdf5_c_compiler_executable} -show
  OUTPUT_VARIABLE _hdf5_c_compile_line
  OUTPUT_STRIP_TRAILING_WHITESPACE)

if (_hdf5_c_compile_line)

  string (REGEX MATCHALL "-I([^\" ]+)" include_path_flags "${_hdf5_c_compile_line}")
  foreach (ipath ${include_path_flags})
    string (REGEX REPLACE "^-I" "" ipath ${ipath})
    string (REGEX REPLACE "//" "/" ipath ${ipath})
    list (APPEND include_paths ${ipath})
  endforeach ()

  string (REGEX MATCHALL "-L([^\" ]+|\"[^\"]+\")" library_path_flags "${_hdf5_c_compile_line}")
  foreach (lpath ${library_path_flags})
    string (REGEX REPLACE "^-L" "" lpath ${lpath})
    string (REGEX REPLACE "//" "/" lpath ${lpath})
    list (APPEND library_paths ${lpath})
  endforeach ()
  
  string (REGEX MATCHALL "[, ]-l([^\", ]+)" library_name_flags "${_hdf5_c_compile_line}")
  foreach (lib ${library_name_flags})
    string (REGEX REPLACE "^[, ]-l" "" lib ${lib})
    list (APPEND libraries ${lib})
  endforeach ()
  
endif ()

foreach (dir ${include_paths})
  list (APPEND HDF5_C_INCLUDE_DIRS "${dir}/include")
endforeach ()

find_path (HDF5_C_INCLUDE_DIR hdf5.h
  HINTS ${_hdf5_root}/include)

mark_as_advanced (HDF5_C_INCLUDE_DIR)

list (APPEND HDF5_INCLUDE_DIRS ${HDF5_C_INCLUDE_DIR})

set (_hdf5_libraries hdf5_hl hdf5)
foreach (_lib ${_hdf5_libraries})
  find_library(
  _hdf5_c_library_${_lib} 
  NAMES ${_lib}
  HINTS ${_hdf5_root}/lib
  NO_DEFAULT_PATH)
  list (APPEND HDF5_C_LIBRARIES "${_hdf5_c_library_${_lib}}")
endforeach ()

foreach (_lib ${libraries})
  find_library(
    found_${_lib} 
    NAMES ${_lib}
    HINTS ${library_paths})
  list (APPEND _additional_libraries "${found_${_lib}}")
endforeach ()

foreach (_lib ${_additional_libraries})
  list (APPEND HDF5_C_LIBRARIES "${_lib}")
endforeach ()

mark_as_advanced (HDF5_C_LIBRARIES)

# If the HDF5 include directory was found, open H5pubconf.h to determine if
# HDF5 was compiled with parallel IO support

set (HDF5_IS_PARALLEL FALSE)
foreach (_dir IN LISTS HDF5_INCLUDE_DIRS)
  if (EXISTS "${_dir}/H5pubconf.h")
    file (STRINGS "${_dir}/H5pubconf.h"
      HDF5_HAVE_PARALLEL_DEFINE
      REGEX "HAVE_PARALLEL 1")
    if (HDF5_HAVE_PARALLEL_DEFINE)
      set (HDF5_IS_PARALLEL TRUE)
    endif ()
  endif ()
endforeach ()
set (HDF5_IS_PARALLEL ${HDF5_IS_PARALLEL} CACHE BOOL
  "HDF5 library compiled with parallel IO support")
mark_as_advanced (HDF5_IS_PARALLEL)

message ("   HDF5 C include directory : ${HDF5_C_INCLUDE_DIR}")
message ("   HDF5 C library directory : ${HDF5_C_LIBRARY_DIR}")
message ("   HDF5 C libraries         :")
foreach (_lib ${HDF5_C_LIBRARIES})
  message ("   ${_lib}")
endforeach () 

message ("   HDF5 parallel support    : ${HDF5_IS_PARALLEL}")

set (HDF5_INCLUDE_DIRS ${HDF5_C_INCLUDE_DIR})
set (HDF5_LIBRARY_DIRS ${HDF5_C_LIBRARY_DIR})
set (HDF5_LIBRARIES    ${HDF5_C_LIBRARIES})

include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (HDF5 DEFAULT_MSG _hdf5_root)

mark_as_advanced (HDF5_INCLUDE_DIRS HDF5_LIBRARIES HDF5_LIBRARY_DIRS)
