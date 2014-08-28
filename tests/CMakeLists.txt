set(CMAKE_INCLUDE_CURRENT_DIR ON)

find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

include_directories(${CMAKE_SOURCE_DIR}/${CMAKE_PROJECT_NAME})
include_directories(${CMAKE_SOURCE_DIR}/${CMAKE_PROJECT_NAME}/accelerators)

#############################################
set(pair_extraction_SRCS 
    pair_extraction.cc 
)
add_executable(pair_extraction ${pair_extraction_SRCS})
add_dependencies(buildtests pair_extraction)
add_test(pair_extraction pair_extraction)
