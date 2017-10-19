project(Super4PCS-CGALWrapper)

find_package( CGAL QUIET )
file(GLOB_RECURSE Demo_CGAL_headers
                  ${CMAKE_CURRENT_SOURCE_DIR}/*.h
                  ${CMAKE_CURRENT_SOURCE_DIR}/*.hpp)


if( CGAL_FOUND AND CGAL_VERSION AND ${CGAL_VERSION} VERSION_EQUAL 4.11 )
  message ("CGAL ${CGAL_VERSION} found. Compiling wrapper demo" )

  set(Demo_CGAL_SRC ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp )

  include( ${CGAL_USE_FILE} )

  include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )
  add_executable(${PROJECT_NAME} ${Demo_CGAL_SRC} ${Demo_CGAL_headers})
  target_link_libraries(${PROJECT_NAME} super4pcs_algo super4pcs_utils)

  install( TARGETS ${PROJECT_NAME} RUNTIME DESTINATION bin)

endif()

install( DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/
         DESTINATION include
         FILES_MATCHING REGEX ".*(h|hpp)$")

