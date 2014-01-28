file (GLOB supportfiles "*.c")

add_library (support ${supportfiles})

set_target_properties (support PROPERTIES LINK_FLAGS "${COMMON_LINK_FLAGS}")
