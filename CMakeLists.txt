cmake_minimum_required (VERSION 2.6)
project (Super4PCS)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

set(SRC_DIR ${PROJECT_SOURCE_DIR}/super4pcs)
 
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
set (ENABLE_TIMING 0)

 
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


################################################################################
## Apply user-defined configuration
################################################################################
if (${ENABLE_TIMING})
    ADD_DEFINITIONS( -DTEST_GLOBAL_TIMINGS )
    message ("OPT: Enable timing recording")
endif (${ENABLE_TIMING})



################################################################################
## Dependencies
################################################################################

add_subdirectory(${PROJECT_SOURCE_DIR}/3rdparty/chealpix EXCLUDE_FROM_ALL)
include_directories(${Chealpix_INCLUDE_DIR})

# switch from C to CXX11
ADD_DEFINITIONS( -std=c++11 )

find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR}) 

find_package( OpenCV REQUIRED core highgui)

# Add package imgcodec required for cv::imread with OpenCV 3.0
if (${OpenCV_VERSION} GREATER "3.0.0")
    find_package( OpenCV REQUIRED imgcodecs )
endif (${OpenCV_VERSION} GREATER "3.0.0")

include_directories( ${OpenCV_INCLUDE_DIRS} )


find_path(ANN_DIR ann_1.1.2)
IF(DEFINED ANN_DIR)
  set(ANN_INCLUDE_DIR "${ANN_DIR}/include")
  include_directories( ${ANN_INCLUDE_DIR} )
  
  set(ANN_LIBS "${ANN_DIR}/lib/libANN.a")
ENDIF(DEFINED ANN_DIR)


################################################################################
## Tests
################################################################################

include(cmake/ConfigureTesting.cmake)
add_subdirectory(tests EXCLUDE_FROM_ALL)

################################################################################
## Project files
################################################################################

include_directories(${SRC_DIR})
set(Super4PCS_SRC
    ${SRC_DIR}/4pcs.cc
    ${SRC_DIR}/super4pcs.cc
)
set(Super4PCS_INCLUDE
    ${SRC_DIR}/shared4pcs.h
    ${SRC_DIR}/4pcs.h
    ${SRC_DIR}/pairCreationFunctor.h
)

add_subdirectory(${SRC_DIR}/io)
add_subdirectory(${SRC_DIR}/accelerators)

add_library(super4pcs_algo  ${Super4PCS_SRC} ${Super4PCS_INCLUDE})
 
# add the executable
add_executable(${PROJECT_NAME} ${SRC_DIR}/super4pcs_test.cc)

TARGET_LINK_LIBRARIES(super4pcs_algo super4pcs_io super4pcs_accel)
TARGET_LINK_LIBRARIES(${PROJECT_NAME} super4pcs_algo)
TARGET_LINK_LIBRARIES(${PROJECT_NAME} ${OpenCV_LIBS} )
TARGET_LINK_LIBRARIES(${PROJECT_NAME} ${ANN_LIBS} )
TARGET_LINK_LIBRARIES(${PROJECT_NAME} ${Chealpix_LIBS} )
