set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(accel_ROOT "${SRC_DIR}/accelerators")

find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

#############################################
set(accel_SRC
    ${accel_ROOT}/normalHealSet.cpp
    PARENT_SCOPE)
set(accel_INCLUDE 
    ${accel_ROOT}/intersection.h
    ${accel_ROOT}/node.h
    ${accel_ROOT}/normalHealSet.h
    ${accel_ROOT}/normalset.h
    ${accel_ROOT}/normalset.hpp
    ${accel_ROOT}/primitives.h
    ${accel_ROOT}/utils.h 
    PARENT_SCOPE)


