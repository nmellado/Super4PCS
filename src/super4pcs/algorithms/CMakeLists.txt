project (Super4PCS-AlgoLib)

set(algo_ROOT "${SRC_DIR}/super4pcs/algorithms")

file(GLOB_RECURSE algo_sources ${algo_ROOT}/*.cc)
file(GLOB_RECURSE algo_headers ${algo_ROOT}/*.h ${algo_ROOT}/*.hpp)

add_library(super4pcs_algo ${algo_headers} ${algo_sources})
# we need that since super4pcs_accel is an header only library
if(SUPER4PCS_USE_CHEALPIX)
    target_link_libraries(super4pcs_algo ${Chealpix_LIBS} )
endif(SUPER4PCS_USE_CHEALPIX)
target_link_libraries(super4pcs_algo super4pcs_io super4pcs_accel)

set(targets_export_name "${PROJECT_NAME}Targets")
install( TARGETS super4pcs_algo
    EXPORT "${targets_export_name}"
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    INCLUDES DESTINATION include  )
