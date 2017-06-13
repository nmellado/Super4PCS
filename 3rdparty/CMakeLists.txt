
# EXTERNALS
# have ExternalProject available
include(ExternalProject)

ExternalProject_Add(
        cfitsio
        # where the source will live
        SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/cfitsio"
        INSTALL_DIR "${CMAKE_BINARY_DIR}/3rdparty/cfitsio"

        CMAKE_ARGS
            -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
            -DBUILD_SHARED_LIBS=True
            -DCMAKE_C_FLAGS=-lm

        # override default behaviours
        UPDATE_COMMAND ""
)


add_subdirectory(${PROJECT_SOURCE_DIR}/3rdparty/chealpix EXCLUDE_FROM_ALL)
include_directories(${Chealpix_INCLUDE_DIR})
