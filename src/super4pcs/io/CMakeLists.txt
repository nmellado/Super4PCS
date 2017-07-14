set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(io_ROOT "${SRC_DIR}/io")

#############################################
set(io_SRC
    ${io_ROOT}/io.cc)
set(io_INCLUDE
    ${io_ROOT}/io.h
    ${io_ROOT}/io_ply.h)

# use library target, to make accessible to test routines
add_library(super4pcs_io ${io_SRC} ${io_INCLUDE})


if(IO_USE_OPENCV)
    find_package( OpenCV REQUIRED core highgui)
    if (OpenCV_FOUND)
        message("[IO] texture loading with OpenCV enabled")
        add_definitions(-DUSE_OPENCV)
        include_directories( ${OpenCV_INCLUDE_DIRS} )
        target_link_libraries(super4pcs_io ${OpenCV_LIBS} )
    endif (OpenCV_FOUND)
endif(IO_USE_OPENCV)


