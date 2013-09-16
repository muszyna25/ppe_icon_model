file (GLOB blasfiles "*.f")

set_source_files_properties (${blasfiles} PROPERTIES COMPILE_FLAGS "$ENV{FFLAGS}")

add_library (blas ${blasfiles})

set_target_properties (blas PROPERTIES LINK_FLAGS "${COMMON_LINK_FLAGS}")
