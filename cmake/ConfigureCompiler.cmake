if( NOT cmake_build_type_tolower STREQUAL "release" )
  add_definitions(-DDEBUG)
endif()

if (MSVC)
  if (MSVC_VERSION LESS 1900)
        message(FATAL_ERROR "Requires Microsoft Visual Studio Compiler version 14.0 or above.")
  endif()

  # remove exceptions from default args
  add_definitions(-D_HAS_EXCEPTIONS=0)
  # disable secure CRT warnings
  add_definitions(-D_CRT_SECURE_NO_WARNINGS)
  add_definitions(-D_SCL_SECURE_NO_WARNINGS)
endif()

set (CMAKE_CXX_STANDARD 11)

if (NOT CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  if( cmake_build_type_tolower STREQUAL "release" )
    find_package(OpenMP)

    if(OPENMP_FOUND)
      message(STATUS "Enable OpenMP")
      add_definitions("-DSUPER4PCS_USE_OPENMP")
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    endif(OPENMP_FOUND)
  endif()
endif()

find_package(Meshlab QUIET)

if(MESHLAB_FOUND)
    # To ease use in shared libraries, even compiling statics Sup4pcs libs
    set(CMAKE_POSITION_INDEPENDENT_CODE ON)
    message(STATUS "Enable position independent code")
endif(MESHLAB_FOUND)
