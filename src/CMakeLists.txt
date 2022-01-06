#-------------------------------------------------------------------------------
# Build PolyClipper packages 
#-------------------------------------------------------------------------------
include_directories(.)

# PolyClipper itself is entirely C++ headers
add_library(PolyClipperAPI INTERFACE)
target_include_directories(
  PolyClipperAPI
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
install(FILES       ${PolyClipper_headers}
        DESTINATION include)

# Are we building the Python bindings?
if(NOT ENABLE_CXXONLY)
  add_subdirectory(PYB11)
endif()
