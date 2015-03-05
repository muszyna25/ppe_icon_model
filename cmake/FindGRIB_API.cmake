find_package(PkgConfig QUIET)

# C interface

pkg_check_modules(MY_GRIB_API QUIET grib_api)

set(GRIB_API_DEFINITIONS ${MY_GRIB_API_CFLAGS_OTHER})

find_path(GRIB_API_INCLUDE_DIR NAMES grib_api.h
  HINTS
  ${MY_GRIB_API_INCLUDE_DIRS})

find_library(GRIB_API_LIBRARIES NAMES grib_api
  HINTS
  ${MY_GRIB_API_LIBRARY_DIRS})

if(MY_GRIB_API_VERSION)
    set(GRIB_API_VERSION ${MY_GRIB_API_VERSION})
else()
    set(GRIB_API_VERSION "")
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(grib_api
                                  REQUIRED_VARS GRIB_API_LIBRARIES GRIB_API_INCLUDE_DIR
                                  VERSION_VAR GRIB_API_VERSION_STRING)

mark_as_advanced(GRIB_API_INCLUDE_DIR GRIB_API_LIBRARIES)

# Fortran interface

pkg_check_modules(MY_GRIB_API_F90 QUIET grib_api_f90)

set(GRIB_API_F90_DEFINITIONS ${MY_GRIB_API_F90_CFLAGS_OTHER})

find_path(GRIB_API_F90_INCLUDE_DIR NAMES grib_api_f77.h
  HINTS
  ${MY_GRIB_API_F90_INCLUDE_DIRS})

find_library(GRIB_API_F90_LIBRARIES NAMES grib_api_f90
  HINTS
  ${MY_GRIB_API_F90_LIBRARY_DIRS})

if(MY_GRIB_API_F90_VERSION)
    set(GRIB_API_F90_VERSION ${MY_GRIB_API_F90_VERSION})
else()
    set(GRIB_API_F90_VERSION "")
endif()

find_package_handle_standard_args(grib_api_f90
                                  REQUIRED_VARS GRIB_API_F90_LIBRARIES GRIB_API_F90_INCLUDE_DIR
                                  VERSION_VAR GRIB_API_F90_VERSION_STRING)
