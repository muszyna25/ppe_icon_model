# Find HDF5 includes and library
#
#  HDF5_INCLUDE_DIRS - where to find netcdf include, etc.
#  HDF5_LIBRARIES    - List of libraries when using netcdf.
#  HDF5_FOUND        - True if netcdf libraries are found.
#
#------------------------------------------------------------------------------
#
#include (SelectLibraryConfigurations)
# Parse a compile line for definitions, includes, library paths, and libraries.
#
macro (_hdf5_parse_compile_line
    compile_line
    include_paths
    definitions
    library_paths
    libraries)

  # Match the include paths
  string (REGEX MATCHALL "-I([^\" ]+)" include_path_flags "${${compile_line}}")
  foreach (ipath ${include_path_flags})
    string (REGEX REPLACE "^-I" "" ipath ${ipath})
    string (REGEX REPLACE "//" "/" ipath ${ipath})
    list (APPEND ${include_paths} ${ipath})
  endforeach()
  
  # Match the definitions
  string (REGEX MATCHALL "-D[^ ]*" definition_flags "${${compile_line}}" )
  foreach (def ${definition_flags})
    list (APPEND ${definitions} ${def})
  endforeach()
  
  # Match the library paths
  string( REGEX MATCHALL "-L([^\" ]+|\"[^\"]+\")" library_path_flags "${${compile_line}}")
   
  foreach (lpath ${library_path_flags})
    string (REGEX REPLACE "^-L" "" lpath ${lpath})
    string (REGEX REPLACE "//" "/" lpath ${lpath})
    list (APPEND ${library_paths} ${lpath})
  endforeach()
  
  # now search for the library names specified in the compile line (match -l...)
  # match only -l's preceded by a space or comma
  # this is to exclude directory names like xxx-linux/
  string (REGEX MATCHALL "[, ]-l([^\", ]+)" library_name_flags "${${compile_line}}" )
  # strip the -l from all of the library flags and add to the search list
  foreach (lib ${library_name_flags})
    string (REGEX REPLACE "^[, ]-l" "" lib ${lib})
    list (APPEND ${libraries} ${lib})
  endforeach()

endmacro()
#
#------------------------------------------------------------------------------
#
message (STATUS "Looking for hdf5 library")

set (_hdf5_root "$ENV{HDF5ROOT}")

find_program (_hdf5_c_compiler_executable
  NAMES h5cc
  HINTS ENV HDF5ROOT
  PATH_SUFFIXES bin
  DOC "hdf5 compiler wrapper, used to retrieve hdf5 compile flags only."
  NO_DEFAULT_PATH)
mark_as_advanced (_hdf5_c_compiler_executable)


if (HDF5_INCLUDE_DIRS AND HDF5_LIBRARIES)

  # Do nothing: we already have HDF5_INCLUDE_PATH and HDF5_LIBRARIES in the
  # cache, it would be a shame to override them

else()

  _HDF5_invoke_compiler( C HDF5_C_COMPILE_LINE HDF5_C_RETURN_VALUE )
  _HDF5_invoke_compiler( CXX HDF5_CXX_COMPILE_LINE HDF5_CXX_RETURN_VALUE )
  
  if( NOT HDF5_FIND_COMPONENTS )

    set( HDF5_LANGUAGE_BINDINGS "C" )

  else()

    # add the extra specified components, ensuring that they are valid.
    foreach( component ${HDF5_FIND_COMPONENTS} )
      list( FIND HDF5_VALID_COMPONENTS ${component} component_location )
      if( ${component_location} EQUAL -1 )
        message( FATAL_ERROR 
          "\"${component}\" is not a valid HDF5 component." )
      else()
        list( APPEND HDF5_LANGUAGE_BINDINGS ${component} )
      endif()
    endforeach()
  endif()
  
  # seed the initial lists of libraries to find with items we know we need
  set( HDF5_C_LIBRARY_NAMES_INIT hdf5_hl hdf5 )
  
  foreach( LANGUAGE ${HDF5_LANGUAGE_BINDINGS} )
    if( HDF5_${LANGUAGE}_COMPILE_LINE )
      _HDF5_parse_compile_line( HDF5_${LANGUAGE}_COMPILE_LINE
        HDF5_${LANGUAGE}_INCLUDE_FLAGS
        HDF5_${LANGUAGE}_DEFINITIONS
        HDF5_${LANGUAGE}_LIBRARY_DIRS
        HDF5_${LANGUAGE}_LIBRARY_NAMES
        )
      
      # take a guess that the includes may be in the 'include' sibling directory
      # of a library directory.
      foreach( dir ${HDF5_${LANGUAGE}_LIBRARY_DIRS} )
        list( APPEND HDF5_${LANGUAGE}_INCLUDE_FLAGS ${dir}/../include )
      endforeach()
    endif()
    
    # set the definitions for the language bindings.
    list( APPEND HDF5_DEFINITIONS ${HDF5_${LANGUAGE}_DEFINITIONS} )
    
    # find the HDF5 include directories
    find_path( HDF5_${LANGUAGE}_INCLUDE_DIR hdf5.h
      HINTS
      ${HDF5_${LANGUAGE}_INCLUDE_FLAGS}
      ENV
      HDF5_ROOT
      PATHS
      $ENV{HOME}/.local/include
      PATH_SUFFIXES
      include
      Include
      )
    mark_as_advanced( HDF5_${LANGUAGE}_INCLUDE_DIR )
    list( APPEND HDF5_INCLUDE_DIRS ${HDF5_${LANGUAGE}_INCLUDE_DIR} )
    
    set( HDF5_${LANGUAGE}_LIBRARY_NAMES
      ${HDF5_${LANGUAGE}_LIBRARY_NAMES_INIT}
      ${HDF5_${LANGUAGE}_LIBRARY_NAMES} )
    
    # find the HDF5 libraries
    foreach( LIB ${HDF5_${LANGUAGE}_LIBRARY_NAMES} )
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
        HINTS ${HDF5_${LANGUAGE}_LIBRARY_DIRS}
        ENV HDF5_ROOT
        PATH_SUFFIXES lib Lib )
      find_library( HDF5_${LIB}_LIBRARY_RELEASE
        NAMES ${THIS_LIBRARY_SEARCH_RELEASE}
        HINTS ${HDF5_${LANGUAGE}_LIBRARY_DIRS}
        ENV HDF5_ROOT
        PATH_SUFFIXES lib Lib )
      select_library_configurations( HDF5_${LIB} )
      # even though we adjusted the individual library names in
      # select_library_configurations, we still need to distinguish
      # between debug and release variants because HDF5_LIBRARIES will
      # need to specify different lists for debug and optimized builds.
      # We can't just use the HDF5_${LIB}_LIBRARY variable (which was set
      # up by the selection macro above) because it may specify debug and
      # optimized variants for a particular library, but a list of
      # libraries is allowed to specify debug and optimized only once.
      list( APPEND HDF5_${LANGUAGE}_LIBRARIES_DEBUG
        ${HDF5_${LIB}_LIBRARY_DEBUG} )
      list( APPEND HDF5_${LANGUAGE}_LIBRARIES_RELEASE
        ${HDF5_${LIB}_LIBRARY_RELEASE} )
    endforeach()
    list( APPEND HDF5_LIBRARY_DIRS ${HDF5_${LANGUAGE}_LIBRARY_DIRS} )
    
    # Append the libraries for this language binding to the list of all
    # required libraries.
    list( APPEND HDF5_LIBRARIES_DEBUG
      ${HDF5_${LANGUAGE}_LIBRARIES_DEBUG} )
    list( APPEND HDF5_LIBRARIES_RELEASE
      ${HDF5_${LANGUAGE}_LIBRARIES_RELEASE} )
  endforeach()
  
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
  
endif()


include (FindPackageHandleStandardArgs)
find_package_handle_standard_args(HDF5 DEFAULT_MSG HDF5_LIBRARIES HDF5_INCLUDE_DIRS)

mark_as_advanced(HDF5_INCLUDE_DIRS HDF5_LIBRARIES HDF5_DEFINTIONS HDF5_LIBRARY_DIRS HDF5_C_COMPILER_EXECUTABLE)
