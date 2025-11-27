include_directories(${GLTOOL_SOURCE_DIR} ${FUROO_SOURCE_DIR}/src ${PONOS_INCLUDES})
SET(TEST_SOURCES
        test_utils.cpp
        )
SET(TEST_HEADERS
        test_utils.h
        )
#add_executable(advection_test advection_test.cpp ${GLTOOL_HEADERS} ${GLTOOL_SOURCES} ${TEST_SOURCES} ${TEST_HEADERS})
#target_link_libraries(advection_test furoo ${PONOS_LIBRARIES})
#add_executable(particle_grid_test particle_grid_test.cpp ${GLTOOL_HEADERS} ${GLTOOL_SOURCES} ${TEST_SOURCES} ${TEST_HEADERS})
#target_link_libraries(particle_grid_test furoo ${PONOS_LIBRARIES})
#add_executable(pic_qt_test pic_qt_test.cpp ${GLTOOL_HEADERS} ${GLTOOL_SOURCES} ${TEST_SOURCES} ${TEST_HEADERS})
#target_link_libraries(pic_qt_test furoo ${PONOS_LIBRARIES})
#add_executable(pic_adap_test pic_adap_test.cpp ${GLTOOL_HEADERS} ${GLTOOL_SOURCES} ${TEST_SOURCES} ${TEST_HEADERS})
#target_link_libraries(pic_adap_test furoo ${PONOS_LIBRARIES})
#add_executable(simulation2 simulation2_test.cpp ${GLTOOL_HEADERS}
#        ${GLTOOL_SOURCES} ${TEST_SOURCES} ${TEST_HEADERS})
#target_link_libraries(simulation2 furoo ${PONOS_LIBRARIES})
add_executable(cell_graph2_test cell_graph2_test.cpp ${GLTOOL_HEADERS}
        ${GLTOOL_SOURCES} ${TEST_SOURCES} ${TEST_HEADERS})
#target_include_directories(cell_graph2_test PUBLIC ${FUROO_SOURCE_DIR}/src ${PONOS_INCLUDES})
target_link_libraries(cell_graph2_test furoo ${PONOS_LIBS})
add_executable(cell_graph3_test cell_graph3_test.cpp ${GLTOOL_HEADERS}
        ${GLTOOL_SOURCES} ${TEST_SOURCES} ${TEST_HEADERS})
target_link_libraries(cell_graph3_test furoo ${PONOS_LIBS})
#add_executable(ping_pong_test ping_pong_test.cpp ${GLTOOL_HEADERS}
#        ${GLTOOL_SOURCES} ${TEST_SOURCES} ${TEST_HEADERS})
#target_link_libraries(ping_pong_test furoo ${PONOS_LIBRARIES})
#add_executable(full_cell_graph2_pic full_cell_graph2_pic.cpp ${GLTOOL_HEADERS}
#        ${GLTOOL_SOURCES} ${TEST_SOURCES} ${TEST_HEADERS})
#target_link_libraries(full_cell_graph2_pic furoo ${PONOS_LIBRARIES})
#add_executable(full_cell_centered_graph2_pic full_cell_centered_graph2_pic.cpp ${GLTOOL_HEADERS}
#        ${GLTOOL_SOURCES} ${TEST_SOURCES} ${TEST_HEADERS})
#target_link_libraries(full_cell_centered_graph2_pic furoo ${PONOS_LIBRARIES})


add_executable(transfer_face_vertex_test transfer_face_vertex_test.cpp ${GLTOOL_SOURCES} ${GLTOOL_HEADERS})
target_link_libraries(transfer_face_vertex_test furoo ponos aergia ${PONOS_EXTRA_LIBS})
target_include_directories(transfer_face_vertex_test PUBLIC ${FUROO_SOURCE_DIR}/src ${PONOS_INCLUDES})
add_executable(bilinear_cellgraph bilinear_cellgraph.cpp ${GLTOOL_SOURCES} ${GLTOOL_HEADERS})
target_link_libraries(bilinear_cellgraph furoo ponos aergia ${PONOS_EXTRA_LIBS})
target_include_directories(bilinear_cellgraph PUBLIC ${FUROO_SOURCE_DIR}/src ${PONOS_INCLUDES})
add_executable(custom_fields_cellgraph2 custom_fields_cellgraph2.cpp ${GLTOOL_SOURCES} ${GLTOOL_HEADERS})
target_link_libraries(custom_fields_cellgraph2 furoo ponos aergia ${PONOS_EXTRA_LIBS})
target_include_directories(custom_fields_cellgraph2 PUBLIC ${FUROO_SOURCE_DIR}/src ${PONOS_INCLUDES})
