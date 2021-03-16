#-------------------------------------------------------------------------------
# Build PolyClipper packages 
#-------------------------------------------------------------------------------
include(${CMAKE_MODULE_PATH}/PolyClipperAddLibs.cmake)
include_directories(.)

# PolyClipper itself is entirely C++ headers
add_library(PolyClipper INTERFACE)
target_include_directories(
  PolyClipper
  INTERFACE $<BUILD_INTERFACE:${${PROJECT_NAME}_SOURCE_DIR}/src>
            $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
  )

set(PolyClipper_headers 
    polyclipper2d.hh
    polyclipper2dImpl.hh
    polyclipper3d.hh
    polyclipper3dImpl.hh
    polyclipper_adapter.hh
    polyclipper_plane.hh
    polyclipper_serialize.hh
    polyclipper_serializeImpl.hh
    polyclipper_utilities.hh
    polyclipper_vector2d.hh
    polyclipper_vector3d.hh)
# set(PolyClipper_sources polyclipper2d.cc polyclipper3d.cc)
# polyclipper_add_cxx_library(PolyClipper)
# install(TARGETS              PolyClipper
#         EXPORT               ${PROJECT_NAME}-targets
#         DESTINATION          lib
#         INCLUDES DESTINATION include)
install(FILES       ${PolyClipper_headers}
        DESTINATION include)

# Are we building the Python bindings?
if(NOT ENABLE_CXXONLY)
  add_subdirectory(Pybind11Wraps)
endif()

# install(EXPORT ${PROJECT_NAME}-targets DESTINATION lib/cmake)
