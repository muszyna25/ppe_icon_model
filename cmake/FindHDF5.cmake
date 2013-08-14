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
else()
  set (_hdf5_root "$ENV{HDF5ROOT}")
endif()

find_program (_hdf5_c_compiler_executable
  NAMES h5cc
  HINTS ${_hdf5_root}
  PATH_SUFFIXES bin
  DOC "hdf5 compiler wrapper, used to retrieve hdf5 compile flags only."
  NO_DEFAULT_PATH)
mark_as_advanced (_hdf5_c_compiler_executable)

execute_process (COMMAND ${_hdf5_c_compiler_executable} -show
  OUTPUT_VARIABLE _hdf5_c_compile_line
  OUTPUT_STRIP_TRAILING_WHITESPACE)

set (HDF5_C_LIBRARY_NAMES_INIT hdf5_hl hdf5)

if (_hdf5_c_compile_line)

  string (REGEX MATCHALL "-I([^\" ]+)" include_path_flags "${_hdf5_c_compile_line}")
  foreach (ipath ${include_path_flags})
    string (REGEX REPLACE "^-I" "" ipath ${ipath})
    string (REGEX REPLACE "//" "/" ipath ${ipath})
    list (APPEND include_paths ${ipath})
  endforeach()
  
  string (REGEX MATCHALL "-D[^ ]*" definition_flags "${_hdf5_c_compile_line}" )
  foreach (def ${definition_flags})
    list (APPEND definitions ${def})
  endforeach()
  
  string( REGEX MATCHALL "-L([^\" ]+|\"[^\"]+\")" library_path_flags "${_hdf5_c_compile_line}")
  foreach (lpath ${library_path_flags})
    string (REGEX REPLACE "^-L" "" lpath ${lpath})
    string (REGEX REPLACE "//" "/" lpath ${lpath})
    list (APPEND library_paths ${lpath})
  endforeach()
  
  string (REGEX MATCHALL "[, ]-l([^\", ]+)" library_name_flags "${_hdf5_c_compile_line}")
  foreach (lib ${library_name_flags})
    string (REGEX REPLACE "^[, ]-l" "" lib ${lib})
    list (APPEND libraries ${lib})
  endforeach()
  
endif()

foreach (dir ${include_paths})
  list (APPEND HDF5_C_INCLUDE_FLAGS "${dir}/include")
endforeach()

find_path( HDF5_C_INCLUDE_DIR hdf5.h
  HINTS ${HDF5_C_INCLUDE_FLAGS} ENV HDF5_ROOT
  PATH_SUFFIXES include
)
mark_as_advanced( HDF5_C_INCLUDE_DIR )

list( APPEND HDF5_INCLUDE_DIRS ${HDF5_C_INCLUDE_DIR} )

set( HDF5_C_LIBRARY_NAMES ${HDF5_C_LIBRARY_NAMES_INIT} ${HDF5_C_LIBRARY_NAMES} )

# find the HDF5 libraries
foreach( LIB ${HDF5_C_LIBRARY_NAMES} )
  if( UNIX AND HDF5_USE_STATIC_LIBRARIES )
    # According to bug 1643 on the CMake bug tracker, this is the
    # preferred method for searching for a static library.
    # See http://www.cmake.org/Bug/view.php?id=1643.  We search
    # first for the full static library name, but fall back to a
    # generic search on the name if the static search fails.
    set( THIS_LIBRARY_SEARCH_DEBUG lib${LIB}d.a ${LIB}d )
    set( THIS_LIBRARY_SEARCH_RELEASE lib${LIB}.a ${LIB} )
  else()
    set( THIS_LIBRARY_SEARCH_DEBUG ${LIB}d )
    set( THIS_LIBRARY_SEARCH_RELEASE ${LIB} )
  endif()
  find_library( HDF5_${LIB}_LIBRARY_DEBUG
    NAMES ${THIS_LIBRARY_SEARCH_DEBUG}
    HINTS ${HDF5_C_LIBRARY_DIRS}
    ENV HDF5_ROOT
    PATH_SUFFIXES lib Lib )
  find_library( HDF5_${LIB}_LIBRARY_RELEASE
    NAMES ${THIS_LIBRARY_SEARCH_RELEASE}
    HINTS ${HDF5_C_LIBRARY_DIRS}
    ENV HDF5_ROOT
    PATH_SUFFIXES lib Lib )
  
  include (SelectLibraryConfigurations)
  select_library_configurations( HDF5_${LIB} )

  list( APPEND HDF5_C_LIBRARIES_RELEASE
    ${HDF5_${LIB}_LIBRARY_RELEASE} )
endforeach()
list( APPEND HDF5_LIBRARY_DIRS ${HDF5_C_LIBRARY_DIRS} )

# We may have picked up some duplicates in various lists during the above
# process for the language bindings (both the C and C++ bindings depend on
# libz for example).  Remove the duplicates.
if( HDF5_INCLUDE_DIRS )
  list( REMOVE_DUPLICATES HDF5_INCLUDE_DIRS )
endif()
if( HDF5_LIBRARIES_DEBUG )
  list( REMOVE_DUPLICATES HDF5_LIBRARIES_DEBUG )
endif()
if( HDF5_LIBRARIES_RELEASE )
  list( REMOVE_DUPLICATES HDF5_LIBRARIES_RELEASE )
endif()
if( HDF5_LIBRARY_DIRS )
  list( REMOVE_DUPLICATES HDF5_LIBRARY_DIRS )
endif()

# Construct the complete list of HDF5 libraries with debug and optimized
# variants when the generator supports them.
if( CMAKE_CONFIGURATION_TYPES OR CMAKE_BUILD_TYPE )
  set( HDF5_LIBRARIES
    debug ${HDF5_LIBRARIES_DEBUG}
    optimized ${HDF5_LIBRARIES_RELEASE} )
else()
  set( HDF5_LIBRARIES ${HDF5_LIBRARIES_RELEASE} )
endif()

# If the HDF5 include directory was found, open H5pubconf.h to determine if
# HDF5 was compiled with parallel IO support
set( HDF5_IS_PARALLEL FALSE )
foreach( _dir IN LISTS HDF5_INCLUDE_DIRS )
  if( EXISTS "${_dir}/H5pubconf.h" )
    file( STRINGS "${_dir}/H5pubconf.h"
      HDF5_HAVE_PARALLEL_DEFINE
      REGEX "HAVE_PARALLEL 1" )
    if( HDF5_HAVE_PARALLEL_DEFINE )
      set( HDF5_IS_PARALLEL TRUE )
    endif()
  endif()
endforeach()
set( HDF5_IS_PARALLEL ${HDF5_IS_PARALLEL} CACHE BOOL
  "HDF5 library compiled with parallel IO support" )
mark_as_advanced( HDF5_IS_PARALLEL )

include (FindPackageHandleStandardArgs)
find_package_handle_standard_args(HDF5 DEFAULT_MSG HDF5_LIBRARIES HDF5_INCLUDE_DIRS)

mark_as_advanced(HDF5_INCLUDE_DIRS HDF5_LIBRARIES HDF5_DEFINTIONS HDF5_LIBRARY_DIRS HDF5_C_COMPILER_EXECUTABLE)
