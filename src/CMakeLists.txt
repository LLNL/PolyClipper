#-------------------------------------------------------------------------------
# Build PolyClipper packages 
#-------------------------------------------------------------------------------
include_directories(.)

# PolyClipper itself is entirely C++ headers
add_library(PolyClipperAPI INTERFACE)
target_include_directories(
  PolyClipperAPI
  INTERFACE $<BUILD_INTERFACE:${${PROJECT_NAME}_SOURCE_DIR}/src>
            $<INSTALL_INTERFACE:${POLYCLIPPER_INSTALL_DIR}>
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
        DESTINATION ${POLYCLIPPER_INSTALL_DIR})

# Are we building the Python bindings?
if(${POLYCLIPPER_MODULE_GEN})
  add_subdirectory(PYB11)
endif()
