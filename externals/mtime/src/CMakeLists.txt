file (GLOB mtimeCFiles "*.[c]")

add_library (mtime ${mtimeCFiles})

set_target_properties (mtime PROPERTIES LINK_FLAGS "${COMMON_LINK_FLAGS}")
