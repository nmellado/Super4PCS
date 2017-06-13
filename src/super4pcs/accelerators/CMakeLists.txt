set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(accel_ROOT "${SRC_DIR}/accelerators")

find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR}  ${Chealpix_INCLUDE_DIR})

#############################################
set(accel_SRC
    ${accel_ROOT}/normalHealSet.cpp)
set(accel_INCLUDE 
    ${accel_ROOT}/kdtree.h
    ${accel_ROOT}/pairExtraction/bruteForceFunctor.h
    ${accel_ROOT}/pairExtraction/intersectionFunctor.h
    ${accel_ROOT}/pairExtraction/intersectionNode.h
    ${accel_ROOT}/pairExtraction/intersectionPrimitive.h
    ${accel_ROOT}/normalHealSet.h
    ${accel_ROOT}/normalset.h
    ${accel_ROOT}/normalset.hpp
    ${accel_ROOT}/bbox.h
    ${accel_ROOT}/utils.h)

add_library(super4pcs_accel ${accel_SRC} ${accel_INCLUDE})
TARGET_LINK_LIBRARIES(super4pcs_accel ${Chealpix_LIBS} )
add_dependencies(super4pcs_accel  chealpix)

