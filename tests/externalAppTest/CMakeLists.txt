project(Super4PCS-externalAppTest)
cmake_minimum_required(VERSION 3.3)

set (CMAKE_CXX_STANDARD 11)

find_package(Super4PCS REQUIRED)
include_directories(${Super4PCS_INCLUDE_DIR})
link_directories(${Super4PCS_LIB_DIR})

add_executable(${PROJECT_NAME} main.cpp)
target_link_libraries(${PROJECT_NAME} ${Super4PCS_LIBRARIES})

message (STATUS "Super4PCS_INCLUDE_DIR: ${Super4PCS_INCLUDE_DIR}")
message (STATUS "Super4PCS_LIB_DIR    : ${Super4PCS_LIB_DIR}")
message (STATUS "Super4PCS_LIBRARIES  : ${Super4PCS_LIBRARIES}")


# we also need eigen
find_package( Eigen3 REQUIRED )
include_directories( ${EIGEN3_INCLUDE_DIR} )

install( TARGETS ${PROJECT_NAME} RUNTIME DESTINATION .  )

