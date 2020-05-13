#-------------------------------------------------------------------------------
# polyclipper_add_cxx_library package_name
#-------------------------------------------------------------------------------
function(polyclipper_add_cxx_library package_name)
  if(ENABLE_STATIC_CXXONLY)
    blt_add_library(NAME        ${package_name}
                    HEADERS     ${${package_name}_headers}
                    SOURCES     ${${package_name}_sources}
                    DEPENDS_ON  -Wl,--start-group ${polyclipper_blt_depends} -Wl,--end-group
                    SHARED      FALSE
                    )
  else()
    blt_add_library(NAME        ${package_name}
                    HEADERS     ${${package_name}_headers}
                    SOURCES     ${${package_name}_sources}
                    DEPENDS_ON  -Wl,--start-group  ${polyclipper_blt_depends} -Wl,--end-group
                    SHARED      TRUE
                    )
  endif()

  # Are there any additional depends?
  if(polyclipper_depends)
    add_dependencies(Polyclipper_${package_name} ${polyclipper_depends})
  endif()

  # Set the r-path of the C++ lib such that it is independent of the build dir when installed
  set_target_properties(${package_name} PROPERTIES
                        INSTALL_RPATH   ${POLYCLIPPER_INSTALL_DIR}/lib
                        )
endfunction()

#----------------------------------------------------------------------
# polyclipper_add_pybind11_library package_name
#----------------------------------------------------------------------
function(polyclipper_add_pybind11_library package_name)
  include(${CMAKE_MODULE_PATH}/PYB11Generator.cmake)

  set(PYB11_MODULE_NAME ${package_name})
  PYB11_GENERATE_BINDINGS()

  set(MODULE_NAME ${PYB11_MODULE_NAME})
  set(GENERATED_SOURCE ${PYB11_GENERATED_SOURCE})

  blt_add_library(
    NAME         ${MODULE_NAME}
    SOURCES      ${GENERATED_SOURCE} ${${package_name}_ADDITIONAL_SOURCES}
    DEPENDS_ON   -Wl,--start-group ${POLYCLIPPER_PYTHON_DEPENDS} -Wl,--end-group ${${package_name}_ADDITIONAL_DEPENDS} ${polyclipper_blt_depends}
    INCLUDES     ${${package_name}_ADDITIONAL_INCLUDES}
    OUTPUT_NAME  ${PYB11_MODULE_NAME}
    CLEAR_PREFIX TRUE
    SHARED       TRUE
    )

  install(
    TARGETS ${MODULE_NAME}
    DESTINATION ${POLYCLIPPER_PYTHON_INSTALL}
    )

  # Are there any additional depends?
  if (polyclipper_py_depends OR polyclipper_depends)
    add_dependencies(${MODULE_NAME} ${polyclipper_py_depends} ${polyclipper_depends})
  endif()

  # Set the r-path of the C++ lib such that it is independent of the build dir when installed
  set_target_properties(${MODULE_NAME}  PROPERTIES
                        INSTALL_RPATH   ${POLYCLIPPER_INSTALL_DIR}/lib
                        )
endfunction()
