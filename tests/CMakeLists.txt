set(CMAKE_INCLUDE_CURRENT_DIR ON)

include_directories(${SRC_DIR})
include_directories(${SRC_DIR}/accelerators)
include_directories(${SRC_DIR}/io)

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
## pair extraction
set(pair_extraction_SRCS
    pair_extraction.cc
)
add_executable(pair_extraction ${pair_extraction_SRCS} ${testing_SRCS})
add_dependencies(buildtests pair_extraction)
add_test(NAME pair_extraction
         #CONFIGURATIONS Release
         COMMAND pair_extraction)


#############################################
## quad extraction
set(quad_extraction_SRCS
    quad_extraction.cc
)
add_executable(quad_extraction ${quad_extraction_SRCS} ${testing_SRCS})
add_dependencies(buildtests quad_extraction)
add_test(NAME quad_extraction
         #CONFIGURATIONS Release
         COMMAND quad_extraction)

#############################################
## matching
set(matching_SRCS
    matching.cc
)

find_package(Boost 1.46.0 COMPONENTS system filesystem REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

add_executable(matching ${matching_SRCS} ${testing_SRCS})
add_dependencies(buildtests matching dataset-standford-armadillo dataset-standford-bunny)

set (matching_subtests 0 1)
foreach(subtest_id ${matching_subtests})
    add_test(NAME "matching_${subtest_id}"
             #CONFIGURATIONS Release
             COMMAND matching ${subtest_id})
    set_tests_properties ( "matching_${subtest_id}" PROPERTIES TIMEOUT 600) # in seconds
endforeach()

target_link_libraries(quad_extraction ${Chealpix_LIBS} super4pcs_accel super4pcs_algo)
target_link_libraries(pair_extraction super4pcs_accel super4pcs_algo super4pcs_utils)
TARGET_LINK_LIBRARIES(pair_extraction ${Chealpix_LIBS} )
target_link_libraries(matching super4pcs_io super4pcs_accel super4pcs_algo super4pcs_utils)
target_link_libraries(matching ${Boost_FILESYSTEM_LIBRARY} ${Boost_SYSTEM_LIBRARY})
TARGET_LINK_LIBRARIES(matching ${Chealpix_LIBS} )
