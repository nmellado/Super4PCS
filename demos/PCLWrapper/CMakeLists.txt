project(Super4PCS-PCLWrapper)

find_package( PCL 1.8 REQUIRED QUIET )

set(Demo_PCL_SRC
    ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp
)
file(GLOB_RECURSE Demo_PCL_headers
                  ${CMAKE_CURRENT_SOURCE_DIR}/*.h
                  ${CMAKE_CURRENT_SOURCE_DIR}/*.hpp)

include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )
add_executable(${PROJECT_NAME} ${Demo_PCL_SRC} ${Demo_PCL_headers})
target_link_libraries(${PROJECT_NAME} super4pcs_algo super4pcs_utils)

if( PCL_FOUND )
  # [C/C++]>[General]>[Additional Include Directories]
  include_directories( ${PCL_INCLUDE_DIRS} )

  # [C/C++]>[Preprocessor]>[Preprocessor Definitions]
  add_definitions( ${PCL_DEFINITIONS} )

  # For Use Not PreCompiled Features
  #add_definitions( -DPCL_NO_PRECOMPILE )

  # [Linker]>[General]>[Additional Library Directories]
  link_directories( ${PCL_LIBRARY_DIRS} )

  # [Linker]>[Input]>[Additional Dependencies]
  target_link_libraries( ${PROJECT_NAME} ${PCL_LIBRARIES} )
endif()


install( TARGETS ${PROJECT_NAME}
         RUNTIME DESTINATION bin)
install( DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/
         DESTINATION include
         FILES_MATCHING REGEX ".*(h|hpp)$")
