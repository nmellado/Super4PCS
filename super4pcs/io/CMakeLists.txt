set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(io_ROOT "${SRC_DIR}/io")

find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

#############################################
set(io_SRC
    ${io_ROOT}/io.cc
  PARENT_SCOPE)
set(io_INCLUDE
    ${io_ROOT}/io.h 
    ${io_ROOT}/io_ply.h
  PARENT_SCOPE)


