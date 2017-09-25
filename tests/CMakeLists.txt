project (Super4PCS-Tests)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

include_directories(${SRC_DIR})

#############################################
## start by downloading standford datasets, that will be used
## to test the approach
include(ExternalProject)
ExternalProject_Add(dataset-standford-armadillo
    SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/datasets/armadillo
    URL http://graphics.stanford.edu/pub/3Dscanrep/armadillo/Armadillo_scans.tar.gz
    URL_MD5 437eafa8d3d07067e17a42c727a0c28d
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
)
ExternalProject_Add(dataset-standford-bunny
    SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/datasets/bunny
    URL http://graphics.stanford.edu/pub/3Dscanrep/bunny.tar.gz
    URL_MD5 232bb90cac5f2ac77140e09a59bcffe4
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
)

#############################################
## testing files
set (testing_SRCS
    testing.h
)

#############################################
## externalAppTestD

include(ExternalProject)
ExternalProject_Add(
        externalAppTestTarget
        # where the source will live
        SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/externalAppTest"

        CMAKE_ARGS
            -DCMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_BINARY_DIR}
            -DCMAKE_MODULE_PATH=${CMAKE_MODULE_PATH}
            -DCMAKE_PREFIX_PATH=${CMAKE_INSTALL_PREFIX}/lib/cmake
            -DCMAKE_INSTALL_RPATH=${CMAKE_INSTALL_PREFIX}/lib
            -DEIGEN3_INCLUDE_DIR=${EIGEN3_INCLUDE_DIR}
            -DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}

        # override default behaviours
        UPDATE_COMMAND ""
        EXCLUDE_FROM_ALL
)

add_dependencies(externalAppTestTarget ForceSuper4PCSInstallation)
add_dependencies(buildtests externalAppTestTarget)
add_test(NAME externalAppTest
         COMMAND Super4PCS-externalAppTest)

#############################################
## pair extraction
set(pair_extraction_SRCS
    pair_extraction.cc
)
add_executable(pair_extraction ${pair_extraction_SRCS} ${testing_SRCS})
add_dependencies(buildtests pair_extraction)
add_test(NAME pair_extraction
         #CONFIGURATIONS Release
         COMMAND pair_extraction)
target_link_libraries(pair_extraction super4pcs_accel super4pcs_algo super4pcs_utils)
if(SUPER4PCS_USE_CHEALPIX)
    target_link_libraries(pair_extraction ${Chealpix_LIBS} )
endif(SUPER4PCS_USE_CHEALPIX)

#############################################
## quad extraction
#set(quad_extraction_SRCS
#    quad_extraction.cc
#)
#add_executable(quad_extraction ${quad_extraction_SRCS} ${testing_SRCS})
#add_dependencies(buildtests quad_extraction)
#add_test(NAME quad_extraction
#         #CONFIGURATIONS Release
#         COMMAND quad_extraction)
#target_link_libraries(quad_extraction super4pcs_accel super4pcs_algo)
#if(SUPER4PCS_USE_CHEALPIX)
#    target_link_libraries(quad_extraction ${Chealpix_LIBS} )
#endif(SUPER4PCS_USE_CHEALPIX)

#############################################
## matching
set(matching_SRCS
    matching.cc
)

find_package(Boost 1.54 COMPONENTS system filesystem QUIET)
if (Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
    add_executable(matching ${matching_SRCS} ${testing_SRCS})
    add_dependencies(buildtests matching)
    set (matching_subtests 0 1)
    foreach(subtest_id ${matching_subtests})
        add_test(NAME "matching_${subtest_id}"
                 #CONFIGURATIONS Release
                 COMMAND matching ${subtest_id})
        set_tests_properties ( "matching_${subtest_id}" PROPERTIES TIMEOUT 600) # in seconds
    endforeach()
    target_link_libraries(matching super4pcs_io super4pcs_accel super4pcs_algo super4pcs_utils)
    target_link_libraries(matching ${Boost_FILESYSTEM_LIBRARY} ${Boost_SYSTEM_LIBRARY})
    add_dependencies(matching dataset-standford-armadillo dataset-standford-bunny)
    if(SUPER4PCS_USE_CHEALPIX)
        target_link_libraries(matching ${Chealpix_LIBS} )
    endif(SUPER4PCS_USE_CHEALPIX)
endif (Boost_FOUND)
