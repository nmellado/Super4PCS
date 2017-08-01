set(accel_ROOT "${SRC_DIR}/super4pcs/accelerators")

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
        target_sources(super4pcs_accel INTERFACE ${accel_INCLUDE})
    else ()
        add_library(super4pcs_accel STATIC ${accel_INCLUDE})
        set_target_properties(super4pcs_accel PROPERTIES LINKER_LANGUAGE CXX)
    endif(CMAKE_VERSION VERSION_GREATER 3.1)
endif(SUPER4PCS_USE_CHEALPIX)


install( TARGETS super4pcs_accel
    RUNTIME DESTINATION ${INSTALL_DIR}/bin
    LIBRARY DESTINATION ${INSTALL_DIR}/lib
    ARCHIVE DESTINATION ${INSTALL_DIR}/lib )
