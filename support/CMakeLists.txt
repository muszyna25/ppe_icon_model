file (GLOB supportfiles "*.c")

set_source_files_properties (${supportfiles} PROPERTIES COMPILE_FLAGS "$ENV{FFLAGS}")

add_library (support ${supportfiles})

set_target_properties (support PROPERTIES LINK_FLAGS "${COMMON_LINK_FLAGS}")
