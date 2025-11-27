#include_directories(${GLTOOL_SOURCE_DIR})


#add_executable(full_cell_centered_graph2_pic_cluster full_cell_centered_graph2_pic.cpp)
#target_link_libraries(full_cell_centered_graph2_pic_cluster furoo)

add_executable(lid_driven_cavity lid_drive_cavity.cpp)
target_link_libraries(lid_driven_cavity furoo)

add_executable(poisson poisson.cpp)
target_link_libraries(poisson furoo)

add_executable(pic_custom pic_custom.cpp)
target_link_libraries(pic_custom furoo)

add_executable(poissonNeumann poissonNeumann.cpp)
target_link_libraries(poissonNeumann furoo)

add_executable(classifyParticles classifyParticles.cpp)
target_link_libraries(classifyParticles furoo)
add_executable(stats stats.cpp)
target_link_libraries(stats furoo)
add_executable(cellGraphToQuadMesh cellGraphToQuadMesh.cpp)
target_link_libraries(cellGraphToQuadMesh furoo)
