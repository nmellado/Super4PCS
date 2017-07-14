set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(accel_ROOT "${SRC_DIR}/accelerators")

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
    if(CMAKE_VERSION VERSION_GREATER 3.1)
        add_library(super4pcs_accel INTERFACE)
        target_include_directories(super4pcs_accel INTERFACE ${accel_ROOT})
        target_sources(super4pcs_accel INTERFACE ${accel_INCLUDE})
    else ()
        add_library(super4pcs_accel STATIC ${accel_INCLUDE})
        include_directories(${accel_ROOT})
        set_target_properties(super4pcs_accel PROPERTIES LINKER_LANGUAGE CXX)
    endif(CMAKE_VERSION VERSION_GREATER 3.1)
endif(SUPER4PCS_USE_CHEALPIX)
