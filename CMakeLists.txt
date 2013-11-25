#_________________________________________________________________________________________________
#
cmake_minimum_required (VERSION 2.8 FATAL_ERROR)
#_________________________________________________________________________________________________
#
project (icon C Fortran)
add_definitions (-D__ICON__)
add_definitions (-D__JSBACH__)
#_________________________________________________________________________________________________
#
message (STATUS "Selected system                      : ${CMAKE_SYSTEM_NAME}")
message (STATUS "Selected generator                   : ${CMAKE_GENERATOR}")
message (STATUS "User specified Fortran compiler      : $ENV{FC}")
message (STATUS "Selected Fortran Compiler            : ${CMAKE_Fortran_COMPILER_ID}")
#_________________________________________________________________________________________________
#
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
#_________________________________________________________________________________________________
#
set (C_C99_FLAG "-std=gnu99")
set (CMAKE_C_FLAGS  "${C_C99_FLAG} $ENV{CFLAGS}")

set (Fortran_PREPROCESS_FLAG "-cpp")
set (CMAKE_Fortran_FLAGS  "${Fortran_PREPROCESS_FLAG} $ENV{FCFLAGS}")

set (CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
#_________________________________________________________________________________________________
#
include_directories (
    ${PROJECT_SOURCE_DIR}/src/include 
    ${PROJECT_SOURCE_DIR}/include 
    ${PROJECT_SOURCE_DIR}/support 
    ${PROJECT_SOURCE_DIR}/externals/mtime/include
    ${PROJECT_BINARY_DIR}/lib
)
#_________________________________________________________________________________________________
include(CheckIncludeFiles)

check_include_files(execinfo.h HAVE_EXECINFO_H)
check_include_files(ucontext.h HAVE_UCONTEXT_H)

configure_file(${PROJECT_SOURCE_DIR}/config/config.h.in ${PROJECT_SOURCE_DIR}/include/config.h)
#_________________________________________________________________________________________________
#
#include(CheckLibraryExists)
#
#check_function_exists(clock_gettime HAVE_CLOCK_GETTIME)
#
#if (HAVE_CLOCK_GETTIME)
#   set(RT_LIBRARY "")
#else()
#   check_library_exists(rt clock_gettime "" HAVE_LIBRT)
#   if (HAVE_LIBRT)
#     set(RT_LIBRARY "-lrt")
#   endif ()
#endif ()
#_________________________________________________________________________________________________
#
find_package (MPI)
include_directories (${MPI_Fortran_INCLUDE_PATH})

message (STATUS "MPI C include path                   : ${MPI_C_INCLUDE_PATH}")
message (STATUS "MPI C libraries                      : ${MPI_C_LIBRARIES}") 
message (STATUS "MPI Fortran include path             : ${MPI_Fortran_INCLUDE_PATH}")
message (STATUS "MPI Fortran libraries                : ${MPI_Fortran_LIBRARIES}") 

message (STATUS "MPI mpiexec                          : ${MPIEXEC}") 
message (STATUS "MPI flag for MPI task count          : ${MPIEXEC_NUMPROC_FLAG}") 

find_package (NETCDF)
include_directories (${NETCDF_INCLUDE_DIRS})

if (${NETCDF_REQUIRES_HDF5} STREQUAL "yes")
  find_package(HDF5)
endif()

add_definitions (-DHAVE_LIBNETCDF -DHAVE_CF_INTERFACE)
#_________________________________________________________________________________________________
#
add_subdirectory (blas)
add_subdirectory (lapack)
add_subdirectory (support)
add_subdirectory (externals/mtime/src)
#_________________________________________________________________________________________________
#
add_subdirectory (src)
#_________________________________________________________________________________________________



