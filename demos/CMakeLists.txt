project (Super4PCS-Demos)

add_subdirectory(Super4PCS)
add_subdirectory(MeshlabPlugin)


find_package( PCL 1.8 QUIET )

if( PCL_FOUND )
    add_subdirectory(PCLWrapper)
endif( PCL_FOUND )
