cmake_minimum_required (VERSION 2.6)
project (Super4PCS)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
set(SRC_DIR ${PROJECT_SOURCE_DIR}/src/super4pcs/)

################################################################################
## Configuration file
################################################################################
set (Super4PCS_VERSION_MAJOR 0)
set (Super4PCS_VERSION_MINOR 1)

# configure a header file to pass some of the CMake settings
# to the source code
#configure_file (
#  "${PROJECT_SOURCE_DIR}/super4pcs_test.cc"
#  )

# Set option to enable timing recording
OPTION (SUPER4PCS_COMPILE_TESTS "Enable testing" TRUE)
OPTION (IO_USE_OPENCV "Use OpenCV for texture loading" TRUE)
OPTION (ENABLE_TIMING "Enable computation time recording" FALSE)
OPTION (SUPER4PCS_USE_CHEALPIX "Use Chealpix for orientation filtering (deprecated)" FALSE)
OPTION (DL_DATASETS "Download demo datasets and associated run scripts" FALSE)


# add the binary tree to the search path for include files
include_directories("${PROJECT_BINARY_DIR}")

################################################################################
## GUARDS
################################################################################

# guard against in-source builds (source: Eigen)
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  message(FATAL_ERROR "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there. You may need to remove CMakeCache.txt. ")
endif()

# guard against bad build-type strings (source: Eigen)
if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()
string(TOLOWER "${CMAKE_BUILD_TYPE}" cmake_build_type_tolower)
if(    NOT cmake_build_type_tolower STREQUAL "debug"
   AND NOT cmake_build_type_tolower STREQUAL "release"
   AND NOT cmake_build_type_tolower STREQUAL "relwithdebinfo")
  message(FATAL_ERROR "Unknown build type \"${CMAKE_BUILD_TYPE}\". Allowed values are Debug, Release, RelWithDebInfo (case-insensitive).")
endif()

################################################################################
## Compiler configuration
################################################################################
if( NOT cmake_build_type_tolower STREQUAL "release" )
  add_definitions(-DDEBUG)
endif()

if (MSVC)
  if (MSVC_VERSION LESS 1900)
        message(FATAL_ERROR "Requires Microsoft Visual Studio Compiler version 14.0 or above.")
  endif()
endif()

# Add the c++11 flag, whatever it is
# This test is by-passed for MSVC compiler, which support C++11 features natively
# without requiring a flag (after version 1800)
include(CheckCXXCompilerFlag)
if (NOT MSVC)
  CHECK_CXX_COMPILER_FLAG(-std=c++11 COMPILER_SUPPORTS_CXX11)
  CHECK_CXX_COMPILER_FLAG(-std=c++0x COMPILER_SUPPORTS_CXX0X)
  if(COMPILER_SUPPORTS_CXX11)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
  elseif(COMPILER_SUPPORTS_CXX0X)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
  else()
        message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
  endif()
endif()

################################################################################
## Apply user-defined configuration
################################################################################
if (ENABLE_TIMING)
    ADD_DEFINITIONS( -DTEST_GLOBAL_TIMINGS )
    message (STATUS "OPT: Enable timing recording")
endif (ENABLE_TIMING)

if (DL_DATASETS)
    include(ExternalProject)

    ExternalProject_Add(dl-dataset-demo1
        SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/datasets/demo1
        URL https://www.irit.fr/~Nicolas.Mellado/dl/datasets/point-clouds/Super4PCS-dataset-demo1.zip
        URL_MD5 68796d9b4efa7edee3ed6f55b92cf7ad
        CONFIGURE_COMMAND ""
        BUILD_COMMAND ""
        INSTALL_COMMAND ""
    )

add_custom_target(dl-datasets DEPENDS dl-dataset-demo1)

endif (DL_DATASETS)

################################################################################
## Dependencies
################################################################################

add_subdirectory(${PROJECT_SOURCE_DIR}/3rdparty)
if(SUPER4PCS_USE_CHEALPIX)
    add_definitions(-DSUPER4PCS_USE_CHEALPIX)
    include_directories(${Chealpix_INCLUDE_DIR})
endif(SUPER4PCS_USE_CHEALPIX)

## Eigen: automatic configuration:
##  1. if EIGEN3_INCLUDE_DIR is set, use it directly
##  2. if any, and version >= 3.3.x, use system version
##  3. otherwise, download (if required) and use submodule
if( NOT EIGEN3_INCLUDE_DIR )
    find_package(Eigen3 QUIET)
    if( (NOT Eigen3_FOUND) OR EIGEN_VERSION_NUMBER VERSION_LESS 3.3 )

        set(EIGEN3_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/3rdparty/Eigen")

        if( NOT EXISTS ${EIGEN3_INCLUDE_DIR}/signature_of_eigen3_matrix_library )
            execute_process(COMMAND git submodule update --init -- ${EIGEN3_INCLUDE_DIR}
                    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
        endif( NOT EXISTS ${EIGEN3_INCLUDE_DIR}/signature_of_eigen3_matrix_library )
    endif()
endif( NOT EIGEN3_INCLUDE_DIR )
include_directories(${EIGEN3_INCLUDE_DIR})
message(STATUS "-- Eigen3 root path: ${EIGEN3_INCLUDE_DIR}")

################################################################################
## Tests
################################################################################

if(SUPER4PCS_COMPILE_TESTS)
        include(cmake/ConfigureTesting.cmake)
        add_subdirectory(tests EXCLUDE_FROM_ALL)
endif(SUPER4PCS_COMPILE_TESTS)

################################################################################
## Project files
################################################################################

include_directories(${SRC_DIR})
set(Super4PCS_SRC
    ${SRC_DIR}/algorithms/match4pcsBase.cc
    ${SRC_DIR}/algorithms/4pcs.cc
    ${SRC_DIR}/algorithms/super4pcs.cc
)
set(Super4PCS_INCLUDE
    ${SRC_DIR}/sampling.h
    ${SRC_DIR}/shared4pcs.h
    ${SRC_DIR}/pairExtraction.h
    ${SRC_DIR}/algorithms/match4pcsBase.h
    ${SRC_DIR}/algorithms/4pcs.h
    ${SRC_DIR}/algorithms/super4pcs.h
    ${SRC_DIR}/pairCreationFunctor.h
)

add_subdirectory(${SRC_DIR}/io)
add_subdirectory(${SRC_DIR}/utils)
add_subdirectory(${SRC_DIR}/accelerators)

add_library(super4pcs_algo  ${Super4PCS_SRC} ${Super4PCS_INCLUDE})

# add the executable
add_executable(${PROJECT_NAME} ${SRC_DIR}/super4pcs_test.cc)

TARGET_LINK_LIBRARIES(super4pcs_algo super4pcs_io super4pcs_accel)
TARGET_LINK_LIBRARIES(${PROJECT_NAME} super4pcs_algo super4pcs_utils)
