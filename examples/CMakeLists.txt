set(target furoo_examples)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

file(GLOB SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)

include_directories("${FUROO_SOURCE_DIR}/src")
target_link_libraries(furoo)
add_executable(${target} ${SOURCES})

if(WIN32)
target_link_libraries(${target} furoo gtest)
else(WIN32)
target_link_libraries(${target} furoo gtest gcov pthread)
endif(WIN32)
