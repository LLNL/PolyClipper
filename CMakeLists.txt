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
option(POLYCLIPPER_ENABLE_CXXONLY "Only install the C++ headers" OFF)
option(POLYCLIPPER_MODULE_GEN "Generate the PolyClipper modules" ON)
option(POLYCLIPPER_ENABLE_DOCS "enable the sphinx documentation" OFF)
option(ENABLE_BOUNDCHECKING "enable Gnu compiler bound checking on memory" OFF)

if(ENABLE_CXXONLY)
  set(POLYCLIPPER_ENABLE_CXXONLY ${ENABLE_CXXONLY})
endif()
#-------------------------------------------------------------------------------
# Install
#-------------------------------------------------------------------------------
if(NOT POLYCLIPPER_INSTALL_DIR)
  set(POLYCLIPPER_INSTALL_DIR include)
endif()
message("-- PolyClipper install path: ${CMAKE_INSTALL_PREFIX}/${POLYCLIPPER_INSTALL_DIR}")

#-------------------------------------------------------------------------------
# Include standard build system logic and options / definitions
#-------------------------------------------------------------------------------
if (${ENABLE_BOUNDCHECKING})
  add_definitions(-D_GLIBCXX_DEBUG=1)
endif()

#-------------------------------------------------------------------------------
# Are we building the python package?
# PolyClipper doesn't actually build any compiled objects in the C++ only case,
# but you can still install the C++ headers.
#-------------------------------------------------------------------------------
if (${POLYCLIPPER_ENABLE_CXXONLY})
  add_definitions(-DCXXONLY=1)
  set(POLYCLIPPER_MODULE_GEN OFF)
else()
  if (NOT PYTHON_EXE)
    find_package(Python3 COMPONENTS Interpreter Development)
    set(PYTHON_EXE "${Python3_EXECUTABLE}")
    message("-- PYTHON_EXE: ${PYTHON_EXE}")
  endif()

  # We need where to find PYB11Generator
  if (NOT PYB11GENERATOR_ROOT_DIR)
    set(PYB11GENERATOR_ROOT_DIR "${CMAKE_SOURCE_DIR}/extern/PYB11Generator")
    message("-- PYB11GENERATOR_ROOT_DIR: ${PYB11GENERATOR_ROOT_DIR}")
  endif()
  # Previous variable for PYB11 directory path
  set(PYB11GEN_PATH ${PYB11GENERATOR_ROOT_DIR})
  include(${PYB11GENERATOR_ROOT_DIR}/cmake/PYB11Generator.cmake)

  # Install python module path
  if (NOT POLYCLIPPER_PYTHON_INSTALL)
    #set(POLYCLIPPER_PYTHON_INSTALL ${CMAKE_INSTALL_PREFIX}/lib/python${Python3_VERSION_MAJOR}.${Python3_VERSION_MINOR}/site-packages/PolyClipper)
    set(POLYCLIPPER_PYTHON_INSTALL ${CMAKE_INSTALL_PREFIX}/lib/python/site-packages/PolyClipper)
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
if (${POLYCLIPPER_ENABLE_DOCS})
  add_subdirectory(docs)
endif()

#-------------------------------------------------------------------------------
# Add any tests
#-------------------------------------------------------------------------------
add_subdirectory(test/test_array_vector)
