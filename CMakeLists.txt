# CMakeLists to build the PolyClipper library.
cmake_minimum_required(VERSION 3.14)
project(PolyClipper LANGUAGES CXX)

################################
# Configure CMake
################################
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -w")
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
set(CMAKE_EXPORT_COMPILE_COMMANDS On)

################################
# Configure and Include blt
################################
set(ENABLE_OPENMP ON CACHE BOOL "")

if(NOT POLYCLIPPER_BLT_DIR) 
  set (POLYCLIPPER_BLT_REL_DIR "${PROJECT_SOURCE_DIR}/cmake/blt" CACHE PATH "")
  get_filename_component(POLYCLIPPER_BLT_DIR "${POLYCLIPPER_BLT_REL_DIR}" ABSOLUTE)
endif()

if (NOT EXISTS "${POLYCLIPPER_BLT_DIR}/SetupBLT.cmake")
    message(FATAL_ERROR 
            "${POLYCLIPPER_BLT_DIR} is not present.\n"
            "call cmake with -DPOLYCLIPPER_BLT_DIR=/your/installation/of/blt\n")
endif()

include(${POLYCLIPPER_BLT_DIR}/SetupBLT.cmake)

################################
# Install
################################
set(POLYCLIPPER_INSTALL_DIR "" CACHE STRING "Directory to install Polyclipper.")

################################
# Include standard build system logic and options / definitions
################################
option(ENABLE_CXXONLY "enable C++ only build without python bindings" OFF)

if(ENABLE_CXXONLY)
  set(CMAKE_INSTALL_PREFIX ${CMAKE_BINARY_DIR}/PolyClipper)
  if(POLYCLIPPER_INSTALL_DIR)
    set(CMAKE_INSTALL_PREFIX ${POLYCLIPPER_INSTALL_DIR})
  endif()
else()
  set(CMAKE_INSTALL_PREFIX ${PYTHON_SITE_PACKAGE_DIR}/Polyclipper)
endif()

if(ENABLE_OPENMP)
  list(APPEND polyclipper_blt_depends openmp)
endif()

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    add_definitions("-DDEBUG=1")
else()
    add_definitions("-DDEBUG=0")
endif()

if (ENABLE_BOUNDCHECKING)
  add_definitions(-D_GLIBCXX_DEBUG=1)
endif()

if(ENABLE_CXXONLY)
  add_definitions(-DCXXONLY=1)
endif()

set(POLYCLIPPER_PYTHON_DEPENDS )
set(POLYCLIPPER_PYTHON_INSTALL ${PROJECT_BINARY_DIR}/lib)

#-------------------------------------------------------------------------------
# Prepare to build the src
#-------------------------------------------------------------------------------
add_subdirectory(src)
