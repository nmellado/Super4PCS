set(CMAKE_INCLUDE_CURRENT_DIR ON)

find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

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
## pair extraction
set(pair_extraction_SRCS
    pair_extraction.cc
)
add_executable(pair_extraction ${pair_extraction_SRCS})
add_dependencies(buildtests pair_extraction)
add_test(pair_extraction pair_extraction)

#############################################
## matching
set(matching_SRCS
    matching.cc
)

find_package(Boost 1.55.0 COMPONENTS system filesystem REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

add_executable(matching ${matching_SRCS} ${io_SRC})
add_dependencies(buildtests matching dataset-standford-armadillo dataset-standford-bunny)
add_test(matching matching)

target_link_libraries(matching super4pcs_io super4pcs_accel super4pcs_algo)
target_link_libraries(matching ${Boost_FILESYSTEM_LIBRARY} ${Boost_SYSTEM_LIBRARY})
TARGET_LINK_LIBRARIES(matching ${ANN_LIBS} )
TARGET_LINK_LIBRARIES(matching ${Chealpix_LIBS} )
