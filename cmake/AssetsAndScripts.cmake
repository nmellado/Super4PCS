
# SCRIPTS
file (GLOB Scripts_SRC ${SCRIPTS_DIR}/*.sh
                       ${SCRIPTS_DIR}/*.bat)

install( FILES ${Scripts_SRC}
         DESTINATION ${INSTALL_DIR}/scripts
         PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ)
add_custom_target(super4pcs_scripts_IDE SOURCES ${Scripts_SRC})


# ASSETS
install( DIRECTORY ${ASSETS_DIR}/ DESTINATION ${INSTALL_DIR}/assets )
