PROJECT(GLTOOL)

add_definitions(-DFONT_PATH="${FUROO_SOURCE_DIR}/assets/arial.ttf")

set(CMAKE_EXPORT_COMPILE_COMMANDS 1)

#if (WIN32)
    #   set(PONOS_INCLUDES ${PONOS_SOURCE_DIR} ${AERGIA_SOURCE_DIR} ${GLAD_INCLUDES} ${NANOGUI_INCLUDES} ${EIGEN_INCLUDES} ${NANOVG_INCLUDES} PARENT_SCOPE)
    #   set(PONOS_EXTRA_LIBS ${NANOGUI_LIBS} ${GLAD_LIBS} PARENT_SCOPE)
#    set(PONOS_LIBS
#            "C:/Projects/ponos/build/ponos/ponosd.lib"
#            "C:/Projects/ponos/build/aergia/aergiad.lib"
#            "C:/Projects/ponos/build/external_build/nanogui/nanoguid.lib"
#            "C:/Projects/ponos/build/external_build/glad/gladd.lib"
#            )
#    set(PONOS_INCLUDES
#            "${FUROO_SOURCE_DIR}/gltool/external/ponos/aergia"
#            "${FUROO_SOURCE_DIR}/gltool/external/ponos/ponos"
#            "C:/Projects/ponos/build/external_build/glad/include"
#            "${FUROO_SOURCE_DIR}/gltool/external/ponos/external/nanogui/include"
#            "${FUROO_SOURCE_DIR}/gltool/external/ponos/external/nanogui/ext/nanovg/src"
#            )
#else (WIN32)
    add_subdirectory("${FUROO_SOURCE_DIR}/gltool/external/ponos" "external_build/ponos")
    set(PONOS_LIBS ponos aergia ${PONOS_LIBRARIES})
    set(PONOS_INCLUDES ${PONOS_INCLUDES})
    # set(PONOS_LIBS aergia ponos ${PONOS_EXTRA_LIBS})
#endif (WIN32)


#add_executable(pic2 pic2.cpp
#        boundary_handler2_model.h
#        boundary_handler2_model.cpp
#        simulation_regular_domain_model.h
#        simulation_regular_domain_model.cpp particle_system2_model.h particle_system2_model.cpp ponos_wrapper.h ponos_wrapper.cpp SearchRenderer.cpp SearchRenderer.h)
#target_link_libraries(pic2 furoo ponos aergia ${PONOS_EXTRA_LIBS})
#target_include_directories(pic2 PUBLIC ${FUROO_SOURCE_DIR}/src ${PONOS_INCLUDES})
#add_executable(pic2qt pic2qt.cpp
#        boundary_handler2_model.h
#        boundary_handler2_model.cpp
#        particle_system2_model.h particle_system2_model.cpp ponos_wrapper.h ponos_wrapper.cpp simulation_quad_tree_domain_model.h simulation_quad_tree_domain_model.cpp SearchRenderer.cpp SearchRenderer.h)
#target_link_libraries(pic2qt furoo ${PONOS_LIBRARIES})
set(GLTOOL_HEADERS
        ${GLTOOL_SOURCE_DIR}/external/lodepng.h
        # ${GLTOOL_SOURCE_DIR}/boundary_handler2_model.h
        ${GLTOOL_SOURCE_DIR}/cell_graph2_model.h
        # ${GLTOOL_SOURCE_DIR}/cell_centered_graph2_model.h
        #${GLTOOL_SOURCE_DIR}/simulation_regular_domain_model.h
        ${GLTOOL_SOURCE_DIR}/mesh_model.h
        ${GLTOOL_SOURCE_DIR}/particle_system_model.h
        ${GLTOOL_SOURCE_DIR}/particle_system_model.inl
        #${GLTOOL_SOURCE_DIR}/particle_system2_model.h
        ${GLTOOL_SOURCE_DIR}/ponos_wrapper.h
        ${GLTOOL_SOURCE_DIR}/solid_model.h
        ${GLTOOL_SOURCE_DIR}/graphic_element_set.h
        #${GLTOOL_SOURCE_DIR}/simulation_quad_tree_domain_model.h
        #${GLTOOL_SOURCE_DIR}/SearchRenderer.h
        )

set(GLTOOL_SOURCES
        ${GLTOOL_SOURCE_DIR}/cell_graph2_model.cpp
        ${GLTOOL_SOURCE_DIR}/graphic_element_set.cpp
        #${GLTOOL_SOURCE_DIR}/cell_centered_graph2_model.cpp
        ${GLTOOL_SOURCE_DIR}/external/lodepng.cpp
        #${GLTOOL_SOURCE_DIR}/boundary_handler2_model.cpp
        ${GLTOOL_SOURCE_DIR}/ponos_wrapper.cpp
        #${GLTOOL_SOURCE_DIR}/particle_system2_model.cpp
        #${GLTOOL_SOURCE_DIR}/simulation_quad_tree_domain_model.cpp
        #${GLTOOL_SOURCE_DIR}/SearchRenderer.cpp
        )

 add_subdirectory(tests)
add_subdirectory(debug)
