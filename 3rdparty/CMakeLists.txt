

if (SUPER4PCS_USE_CHEALPIX)

# EXTERNALS
# have ExternalProject available
include(ExternalProject)
message("[Super4PCS] Use Chealpix")

ExternalProject_Add(
        cfitsio
        # where the source will live
        SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/cfitsio"
        INSTALL_DIR "${CMAKE_BINARY_DIR}/3rdparty/cfitsio"

        CMAKE_ARGS
            -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>

        # override default behaviours
        UPDATE_COMMAND ""
)


add_subdirectory(${PROJECT_SOURCE_DIR}/3rdparty/chealpix EXCLUDE_FROM_ALL)
endif (SUPER4PCS_USE_CHEALPIX)
