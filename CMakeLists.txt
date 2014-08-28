cmake_minimum_required (VERSION 2.6)
project (Super4PCS)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
 
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

 
# add the binary tree to the search path for include files
include_directories("${PROJECT_BINARY_DIR}")
# The version number.

################################################################################
## Compiler configuration 
################################################################################
ADD_DEFINITIONS(
    -std=c++11
)


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
## Dependencies
################################################################################

find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR}) 

find_package( OpenCV REQUIRED core highgui)
include_directories( ${OpenCV_INCLUDE_DIRS} )


IF(DEFINED ANN_DIR)
  set(ANN_INCLUDE_DIR "${ANN_DIR}/include")
  include_directories( ${ANN_INCLUDE_DIR} )
  
  set(ANN_LIBS "${ANN_DIR}/lib/libANN.a")
ENDIF(DEFINED ANN_DIR)

add_subdirectory(${PROJECT_SOURCE_DIR}/dep/chealpix)
include_directories(${Chealpix_INCLUDE_DIR})

message ("${Chealpix_LIBS}")

################################################################################
## Project files
################################################################################

include_directories(${CMAKE_SOURCE_DIR})
set(Super4PCS_SRC
    io/io.cc
    4pcs.cc
    super4pcs.cc
    accelerators/normalHealSet.cpp
    super4pcs_test.cc
)
 
# add the executable
add_executable(${PROJECT_NAME} ${Super4PCS_SRC})

#add_dependencies(${PROJECT_NAME} ${Chealpix_LIBS})

TARGET_LINK_LIBRARIES(${PROJECT_NAME} ${OpenCV_LIBS} )
TARGET_LINK_LIBRARIES(${PROJECT_NAME} ${ANN_LIBS} )
TARGET_LINK_LIBRARIES(${PROJECT_NAME} ${Chealpix_LIBS} )
