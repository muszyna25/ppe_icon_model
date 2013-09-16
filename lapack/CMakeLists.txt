file (GLOB lapackfiles "*.f")

set_source_files_properties (${lapackfiles} PROPERTIES COMPILE_FLAGS "$ENV{FFLAGS}")

add_library (lapack ${lapackfiles})

set_target_properties (lapack PROPERTIES LINK_FLAGS "${COMMON_LINK_FLAGS}")
