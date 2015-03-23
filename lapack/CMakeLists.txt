file (GLOB lapackfiles "*.f")

#set_source_files_properties (${lapackfiles} PROPERTIES COMPILE_FLAGS "-std=legacy")

set (CMAKE_Fortran_FLAGS "$ENV{F77FLAGS}")

add_library (lapack ${lapackfiles})

set_target_properties (lapack PROPERTIES LINK_FLAGS "${COMMON_LINK_FLAGS}")
