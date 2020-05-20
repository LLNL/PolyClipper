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
# Configure and Include blt
#-------------------------------------------------------------------------------
set(ENABLE_OPENMP ON CACHE BOOL "")

if (NOT POLYCLIPPER_BLT_DIR) 
  set (POLYCLIPPER_BLT_REL_DIR "${PROJECT_SOURCE_DIR}/cmake/blt" CACHE PATH "")
  get_filename_component(POLYCLIPPER_BLT_DIR "${POLYCLIPPER_BLT_REL_DIR}" ABSOLUTE)
endif()

if (NOT EXISTS "${POLYCLIPPER_BLT_DIR}/SetupBLT.cmake")
    message(FATAL_ERROR 
            "${POLYCLIPPER_BLT_DIR} is not present.\n"
            "call cmake with -DPOLYCLIPPER_BLT_DIR=/your/installation/of/blt\n")
endif()

include(${POLYCLIPPER_BLT_DIR}/SetupBLT.cmake)

#-------------------------------------------------------------------------------
# Install
#-------------------------------------------------------------------------------
set(POLYCLIPPER_INSTALL_DIR "" CACHE STRING "Directory to install PolyClipper.")
set(CMAKE_INSTALL_PREFIX ${CMAKE_BINARY_DIR}/PolyClipper)
if (POLYCLIPPER_INSTALL_DIR)
  set(CMAKE_INSTALL_PREFIX ${POLYCLIPPER_INSTALL_DIR})
endif()
message("-- install path ${CMAKE_INSTALL_PREFIX}")

#-------------------------------------------------------------------------------
# Include standard build system logic and options / definitions
#-------------------------------------------------------------------------------
option(ENABLE_STATIC_CXXONLY "enable C++ only build with static libraries" OFF)
option(ENABLE_CXXONLY "enable C++ only build without python bindings" OFF)

if (ENABLE_STATIC_CXXONLY)
  set(ENABLE_CXXONLY ON)
endif()

if (ENABLE_OPENMP)
  list(APPEND polyclipper_blt_depends openmp)
endif()

if (ENABLE_BOUNDCHECKING)
  add_definitions(-D_GLIBCXX_DEBUG=1)
endif()

#-------------------------------------------------------------------------------
# Are we building the python package?
#-------------------------------------------------------------------------------
if (ENABLE_CXXONLY)
  add_definitions(-DCXXONLY=1)

else()
  if (NOT PYTHON_EXE)
    find_package(Python2 COMPONENTS Interpreter Development)
    set(PYTHON_EXE "${Python2_EXECUTABLE}")
  endif()
  message("-- PYTHON_EXE ${PYTHON_EXE}")

  # We need pybind11 paths to compile
  if (LOOKUP_PYBIND11_INCLUDE_PATH)
    execute_process(COMMAND         ${PYTHON_EXE} -m pybind11 --includes
                    OUTPUT_VARIABLE PYBIND11_INCLUDE_PATH)
    string(REPLACE "-I" "" PYBIND11_INCLUDE_PATH "${PYBIND11_INCLUDE_PATH}")
  elseif (NOT PYBIND11_INCLUDE_PATH)
    set(PYBIND11_INCLUDE_PATH "${CMAKE_MODULE_PATH}/pybind11/include ${Python2_INCLUDE_DIRS}")
  endif()
  set(PYBIND11_INCLUDE_PATHS ${PYBIND11_INCLUDE_PATH})
  separate_arguments(PYBIND11_INCLUDE_PATHS)
  include_directories(${PYBIND11_INCLUDE_PATHS})
  message("-- PYBIND11_INCLUDE_PATH ${PYBIND11_INCLUDE_PATH}")

  # Also PYB11Generator
  if (NOT PYB11GEN_PATH)
    set(PYB11GEN_PATH "${CMAKE_MODULE_PATH}/PYB11Generator")
  endif()
  message("-- PYB11GEN_PATH ${PYB11GEN_PATH}")

  # Install python module path
  if (NOT POLYCLIPPER_PYTHON_INSTALL)
    set(POLYCLIPPER_PYTHON_INSTALL "${POLYCLIPPER_INSTALL_DIR}/lib/python/site-packages/PolyClipper")
  endif()
  #set(POLYCLIPPER_PYTHON_INSTALL "" CACHE STRING "Directory to install PolyClipper Python module.")
  message("-- Python module install path ${POLYCLIPPER_PYTHON_INSTALL}")
endif()

set(POLYCLIPPER_PYTHON_DEPENDS PolyClipper)

#-------------------------------------------------------------------------------
# Prepare to build the src
#-------------------------------------------------------------------------------
add_subdirectory(src)

#-------------------------------------------------------------------------------
# Add the documentation
#-------------------------------------------------------------------------------
add_subdirectory(docs)
