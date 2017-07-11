set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(accel_ROOT "${SRC_DIR}/accelerators")

find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

#############################################

set(accel_INCLUDE
    ${accel_ROOT}/kdtree.h
    ${accel_ROOT}/pairExtraction/bruteForceFunctor.h
    ${accel_ROOT}/pairExtraction/intersectionFunctor.h
    ${accel_ROOT}/pairExtraction/intersectionNode.h
    ${accel_ROOT}/pairExtraction/intersectionPrimitive.h
    ${accel_ROOT}/normalset.h
    ${accel_ROOT}/normalset.hpp
    ${accel_ROOT}/bbox.h
    ${accel_ROOT}/utils.h)

if(SUPER4PCS_USE_CHEALPIX)
    include_directories(${Chealpix_INCLUDE_DIR})
    set(chealpix_files
        ${accel_ROOT}/normalHealSet.cpp
        ${accel_ROOT}/normalHealSet.h)

    add_library(super4pcs_accel ${accel_INCLUDE} ${chealpix_files})
    target_link_libraries(super4pcs_accel ${Chealpix_LIBS} )
    add_dependencies(super4pcs_accel  chealpix)

    add_definitions(-DUSE_OPENCV)
else ()
    add_library(super4pcs_accel INTERFACE)
    target_include_directories(super4pcs_accel INTERFACE ${accel_ROOT})

    if(CMAKE_VERSION VERSION_GREATER 3.1)
        target_sources(super4pcs_accel INTERFACE ${accel_INCLUDE})
    endif(CMAKE_VERSION VERSION_GREATER 3.1)
endif(SUPER4PCS_USE_CHEALPIX)
