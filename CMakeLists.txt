# CMakeLists to build the PolyClipper library.
cmake_minimum_required(VERSION 3.14)
project(PolyClipper LANGUAGES CXX)

#-------------------------------------------------------------------------------
# Configure CMake
#-------------------------------------------------------------------------------
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -w")
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
set(CMAKE_EXPORT_COMPILE_COMMANDS On)

#-------------------------------------------------------------------------------
# Options
#-------------------------------------------------------------------------------
option(ENABLE_CXXONLY "Do not build Python language bindings" OFF)
option(ENABLE_DOCS "enable the sphinx documentation" OFF)
option(ENABLE_BOUNDCHECKING "enable Gnu compiler bound checking on memory" OFF)

#-------------------------------------------------------------------------------
# Install
#-------------------------------------------------------------------------------
set(POLYCLIPPER_INSTALL_DIR ${CMAKE_INSTALL_PREFIX})
message("-- install path: ${CMAKE_INSTALL_PREFIX}")

#-------------------------------------------------------------------------------
# Include standard build system logic and options / definitions
#-------------------------------------------------------------------------------
if (ENABLE_BOUNDCHECKING)
  add_definitions(-D_GLIBCXX_DEBUG=1)
endif()

#-------------------------------------------------------------------------------
# Are we building the python package?
# PolyClipper doesn't actually build any compiled objects in the C++ only case,
# but you can still install the C++ headers.
#-------------------------------------------------------------------------------
if (ENABLE_CXXONLY)
  add_definitions(-DCXXONLY=1)

else()
  if (NOT PYTHON_EXE)
    find_package(Python3 COMPONENTS Interpreter Development)
    set(PYTHON_EXE "${Python3_EXECUTABLE}")
  endif()
  message("-- PYTHON_EXE: ${PYTHON_EXE}")

  # We need where to find PYB11Generator
  if (NOT PYB11GEN_PATH)
    set(PYB11GEN_PATH "${CMAKE_SOURCE_DIR}/extern/PYB11Generator")
  endif()
  set (PYB11GENERATOR_ROOT_DIR ${PYB11GEN_PATH})  # PYB11Generator requires PYB11GENERATOR_ROOT_DIR
  add_subdirectory(${PYB11GEN_PATH}/extern/pybind11)
  include(${PYB11GEN_PATH}/cmake/PYB11Generator.cmake)
  message("-- PYB11GEN_PATH: ${PYB11GEN_PATH}")

  # Install python module path
  if (NOT POLYCLIPPER_PYTHON_INSTALL)
    set(POLYCLIPPER_PYTHON_INSTALL ${CMAKE_INSTALL_PREFIX}/lib/python${Python3_VERSION_MAJOR}.${Python3_VERSION_MINOR}/site-packages/PolyClipper)
  endif()
  message("-- Python module install path ${POLYCLIPPER_PYTHON_INSTALL}")
endif()

set(POLYCLIPPER_PYTHON_DEPENDS )

#-------------------------------------------------------------------------------
# Prepare to build the src
#-------------------------------------------------------------------------------
add_subdirectory(src)

#-------------------------------------------------------------------------------
# Add the documentation
#-------------------------------------------------------------------------------
if (ENABLE_DOCS)
  add_subdirectory(docs)
endif()

#-------------------------------------------------------------------------------
# Add any tests
#-------------------------------------------------------------------------------
add_subdirectory(test/test_array_vector)
