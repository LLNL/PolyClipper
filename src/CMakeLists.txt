#-------------------------------------------------------------------------------
# Build PolyClipper packages 
#-------------------------------------------------------------------------------
include(${CMAKE_MODULE_PATH}/PolyClipperAddLibs.cmake)
include_directories(.)

# Build the PolyClipper library
set(PolyClipper_headers polyclipper.hh)
set(PolyClipper_sources polyclipper2d.cc polyclipper3d.cc)
polyclipper_add_cxx_library(PolyClipper)
install(TARGETS              PolyClipper
        EXPORT               ${PROJECT_NAME}-targets
        DESTINATION          lib
        INCLUDES DESTINATION include)
install(FILES       polyclipper.hh
        DESTINATION include/PolyClipper)

# # Are we building the Python bindings?
# if(NOT ENABLE_CXXONLY)
#   add_subdirectory(src/Pybind11Wraps)
# endif()

install(EXPORT ${PROJECT_NAME}-targets DESTINATION lib/cmake)
