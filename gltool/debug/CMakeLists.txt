include_directories(${GLTOOL_SOURCE_DIR})
include_directories("${GLTOOL_SOURCE_DIR}/debug")
SET(TEST_SOURCES
        test_utils.cpp
        )
SET(TEST_HEADERS
        test_utils.h
        )
#add_executable(full_cell_centered_graph2_pic_debug full_cell_centered_graph2_pic.cpp ${GLTOOL_HEADERS}
#        ${GLTOOL_SOURCES} ${TEST_SOURCES} ${TEST_HEADERS})
#target_link_libraries(full_cell_centered_graph2_pic_debug furoo ${PONOS_LIBRARIES})

#add_executable(visualize_simulation_step visualize_simulation_step.cpp ${GLTOOL_HEADERS}
#        ${GLTOOL_SOURCES} ${TEST_SOURCES} ${TEST_HEADERS})
#target_link_libraries(visualize_simulation_step furoo ${PONOS_LIBRARIES})

#add_executable(flip_debugger flip_debugger.cpp ui.h ui.inl
#${GLTOOL_HEADERS} ${GLTOOL_SOURCES}
#        )

#target_link_libraries(flip_debugger furoo ponos aergia ${PONOS_EXTRA_LIBS})
#target_include_directories(flip_debugger PUBLIC ${FUROO_SOURCE_DIR}/src ${PONOS_INCLUDES})

set(CMAKE_CXX_STANDARD 17)
#set(CMAKE_CXX_FLAGS_DEBUG "-std=c++17 -lstdc++fs")

add_executable(furoo_debugger main.cpp ui.h ui.inl ${GLTOOL_HEADERS} ${GLTOOL_SOURCES})
target_link_libraries(furoo_debugger furoo ${PONOS_LIBS} -lstdc++fs ${PONOS_EXTRA_LIBS})
target_include_directories(furoo_debugger PUBLIC ${FUROO_SOURCE_DIR}/src ${PONOS_INCLUDES} ${GLTOOL_HEADERS})
target_compile_options(furoo_debugger PUBLIC -std=c++17)

add_executable(pic_debugger pic_debugger.cpp ui.h ui.inl ${GLTOOL_HEADERS} ${GLTOOL_SOURCES})
target_link_libraries(pic_debugger furoo ${PONOS_LIBS} -lstdc++fs ${PONOS_EXTRA_LIBS})
target_include_directories(pic_debugger PUBLIC ${FUROO_SOURCE_DIR}/src ${PONOS_INCLUDES} ${GLTOOL_HEADERS})
target_compile_options(pic_debugger PUBLIC -std=c++17)

add_executable(pic_debugger3 pic_debugger3.cpp ui.h ui.inl ${GLTOOL_HEADERS} ${GLTOOL_SOURCES})
target_link_libraries(pic_debugger3 furoo ${PONOS_LIBS} -lstdc++fs ${PONOS_EXTRA_LIBS})
target_include_directories(pic_debugger3 PUBLIC ${FUROO_SOURCE_DIR}/src ${PONOS_INCLUDES} ${GLTOOL_HEADERS})
target_compile_options(pic_debugger3 PUBLIC -std=c++17)

add_executable(custom_simulation_debugger custom_simulation_debugger.cpp ui.h ui.inl ${GLTOOL_HEADERS} ${GLTOOL_SOURCES})
target_link_libraries(custom_simulation_debugger furoo ${PONOS_LIBS} -lstdc++fs ${PONOS_EXTRA_LIBS})
target_include_directories(custom_simulation_debugger PUBLIC ${FUROO_SOURCE_DIR}/src ${PONOS_INCLUDES} ${GLTOOL_HEADERS})
target_compile_options(custom_simulation_debugger PUBLIC -std=c++17)

add_executable(stencil_debugger stencil_debugger.cpp ui.h ui.inl ${GLTOOL_HEADERS} ${GLTOOL_SOURCES})
target_link_libraries(stencil_debugger furoo ${PONOS_LIBS} -lstdc++fs ${PONOS_EXTRA_LIBS})
target_include_directories(stencil_debugger PUBLIC ${FUROO_SOURCE_DIR}/src ${PONOS_INCLUDES} ${GLTOOL_HEADERS})
target_compile_options(stencil_debugger PUBLIC -std=c++17)

add_executable(obj_test obj_test.cpp ui.h ui.inl ${GLTOOL_HEADERS} ${GLTOOL_SOURCES})
target_link_libraries(obj_test furoo ${PONOS_LIBS} -lstdc++fs ${PONOS_EXTRA_LIBS})
target_include_directories(obj_test PUBLIC ${FUROO_SOURCE_DIR}/src ${PONOS_INCLUDES} ${GLTOOL_HEADERS})
target_compile_options(obj_test PUBLIC -std=c++17)