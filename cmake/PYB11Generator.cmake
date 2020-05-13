#-----------------------------------------------------------------------------------
# PYB11_GENERATE_BINDINGS
#     - Generates the Python bindings for each module in the list
#     - Generates python stamp files for listing python dependency file to help
#       detecting changes in the pyb11 python files at build time
#
# Variables that must be set before calling PYB11_GENERATE_BINDINGS:
#   PYB11_MODULE_NAME
#     - Pyb11 module to be generated
#   PYTHON_EXE
#     - Python executable
#   <PYB11_MODULE_NAME>_DEPENDS
#     - Any target dependencies that must be built before generating the module
#
# To get the names of the generated source
# use: ${PYB11_GENERATED_SOURCE}
#-----------------------------------------------------------------------------------

macro(PYB11_GENERATE_BINDINGS PYB11_MODULE_NAME)
  set(PYB11_SOURCE "${PYB11_MODULE_NAME}MOD.py")
  set(PYB11_GENERATED_SOURCE "${PYB11_MODULE_NAME}MOD.cc")

  message("**** ${PYB11_MODULE_NAME}")
  message("**** ${PYB11_SOURCE}")
  message("**** ${PYB11_GENERATED_SOURCE}")

  # List directories in which spheral .py files can be found.
  set(PYTHON_ENV 
      "${CMAKE_MODULE_PATH}/PYB11Generator:"
      "${PROJECT_SOURCE_DIR}/src/Pybind11Wraps:"
      "${PROJECT_SOURCE_DIR}/src/Pybind11Wraps/${PYB11_MODULE_NAME}:")

  # Format list into a one line shell friendly format
  STRING(REPLACE ";" "<->" PYTHON_ENV_STR ${PYTHON_ENV})
  message("**** ${PYTHON_ENV_STR}")

  # Generating python stamp files to detect changes in PYB11_SOURCE and
  # its included modules
  if(EXISTS ${PYTHON_EXE})
    # Python must exist to generate at config time
    if(NOT EXISTS "${CMAKE_CURRENT_BINARY_DIR}/${PYB11_MODULE_NAME}_stamp.cmake")
      # Generate stamp files at config time
      execute_process(COMMAND env PYTHONPATH=\"${PYTHON_ENV_STR}\"
                      ${PYTHON_EXE} ${CMAKE_MODULE_PATH}/moduleCheck.py 
                      ${PYB11_MODULE_NAME}
                      ${CMAKE_CURRENT_SOURCE_DIR}/${PYB11_SOURCE}
                      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                      )
    endif()

    # Include list of dependent python files
    include(${CMAKE_CURRENT_BINARY_DIR}/${PYB11_MODULE_NAME}_stamp.cmake)
  endif()

  # Always regenerate the stamp files at build time. Any change in the stamp file
  # will trigger a rebuild of the target pyb11 module
  add_custom_target(${PYB11_MODULE_NAME}_stamp ALL
                    COMMAND env PYTHONPATH=\"${PYTHON_ENV_STR}\"
                    ${PYTHON_EXE} ${CMAKE_MODULE_PATH}/moduleCheck.py
                    ${PYB11_MODULE_NAME}
                    ${CMAKE_CURRENT_SOURCE_DIR}/${PYB11_SOURCE}
                    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                    #DEPENDS python-install
                    )

  # Generate the actual pyb11 module cpp source file
  add_custom_command(OUTPUT ${PYB11_GENERATED_SOURCE}
                     COMMAND env PYTHONPATH=\"${PYTHON_ENV_STR}\"
                     ${PYTHON_EXE} -c
                     'import os \; print(\"${PYB11_MODULE_NAME} in directory \", os.getcwd()) \;
                     from PYB11Generator import * \;
                     import ${PYB11_MODULE_NAME}MOD as ${PYB11_MODULE_NAME}\; 
                     PYB11generateModule(${PYB11_MODULE_NAME}) '
                     DEPENDS ${PYB11_MODULE_NAME}_stamp ${${PYB11_MODULE_NAME}_DEPENDS} ${PYB11_SOURCE}
                     )

endmacro()
