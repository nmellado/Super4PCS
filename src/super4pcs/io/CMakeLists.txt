project (Super4PCS-IOLib)

set(io_ROOT "${SRC_DIR}/super4pcs/io")

file(GLOB_RECURSE io_sources ${io_ROOT}/*.cc)
file(GLOB_RECURSE io_headers ${io_ROOT}/*.h ${io_ROOT}/*.hpp)

add_library(super4pcs_io ${io_sources} ${io_headers})

if(IO_USE_OPENCV)
    find_package( OpenCV REQUIRED)
    if (OpenCV_FOUND)
        message(STATUS "[Deps] Use OpenCV ${OpenCV_VERSION} for image loading")
        add_definitions(-DUSE_OPENCV)
        include_directories( ${OpenCV_INCLUDE_DIRS} )
        target_link_libraries(super4pcs_io ${OpenCV_LIBS} )
    else ()
        message(STATUS "[Deps] OpenCV not found")
    endif (OpenCV_FOUND)
endif(IO_USE_OPENCV)

set(targets_export_name "${PROJECT_NAME}Targets")
install( TARGETS super4pcs_io
    EXPORT "${targets_export_name}"
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    INCLUDES DESTINATION include  )
