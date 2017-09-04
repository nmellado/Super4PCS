project(Meshlab-GlobalRegistrationPlugin)

find_package(Meshlab QUIET)

if(MESHLAB_FOUND)
    message(STATUS "Preparing Meshlab plugin")

    set(PLUGIN_NAME "filter_globalregistration")
    set(PLUGIN_SRC_DIR              "${CMAKE_CURRENT_SOURCE_DIR}/${PLUGIN_NAME}")
    set(PLUGIN_PRO_INC              "${CMAKE_CURRENT_SOURCE_DIR}/${PLUGIN_NAME}.pro.in")
    set(PLUGIN_RDM_INC              "${CMAKE_CURRENT_SOURCE_DIR}/readme.md.in")
    set(PLUGIN_INSTALL_DIR          "${MESHLAB_PLUGINS_SOURCE_DIR}/${PLUGIN_NAME}")
    set(PLUGIN_RELATIVE_INSTALL_DIR "${MESHLAB_RELATIVE_PLUGINS_SOURCE_DIR}/${PLUGIN_NAME}")
    set(PLUGIN_PRO_GEN              "${PLUGIN_INSTALL_DIR}/${PLUGIN_NAME}.pro")
    set(PLUGIN_RDM_GEN              "${PLUGIN_INSTALL_DIR}/readme.md")

    configure_file(${PLUGIN_PRO_INC} ${PLUGIN_PRO_GEN} @ONLY)
    configure_file(${PLUGIN_RDM_INC} ${PLUGIN_RDM_GEN} @ONLY)

    add_custom_target(internal_install_super4pcs_to_meshlab_externals
       COMMAND ${CMAKE_COMMAND} -E copy_directory "${Super4PCS_INSTALL_DIR}/include/" "${MESHLAB_EXTERNAL_INCLUDE_DIR}"
       COMMAND ${CMAKE_COMMAND} -E copy_directory "${Super4PCS_INSTALL_DIR}/lib/"     "${MESHLAB_EXTERNAL_LIBRARIES_DIR}"
       COMMENT "Install Super4PCS library to Meshlab externals (${MESHLAB_EXTERNAL_INCLUDE_DIR})" VERBATIM
       DEPENDS ForceSuper4PCSInstallation)

    add_custom_target(internal_install_meshlab_plugin_source
       COMMAND ${CMAKE_COMMAND} -E copy_directory "${PLUGIN_SRC_DIR}" "${PLUGIN_INSTALL_DIR}"
       COMMENT "Install plugin sources at ${PLUGIN_INSTALL_DIR}" VERBATIM
       )

    # Main target visible in the project
    add_custom_target(Meshlab_install_plugin
       COMMENT "Performing Meshlab plugin installation (manual compilation is still requires. See ${PLUGIN_INSTALL_DIR}/readme.md)"
       DEPENDS internal_install_super4pcs_to_meshlab_externals internal_install_meshlab_plugin_source)
else()
     message(WARNING "Skipping Meshlab plugin preparation - Meshlab not found")
endif(MESHLAB_FOUND)


