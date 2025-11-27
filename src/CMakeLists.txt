file(GLOB HEADERS
        ${CMAKE_CURRENT_SOURCE_DIR}/*.h
        ${CMAKE_CURRENT_SOURCE_DIR}/animation/*.h
        ${CMAKE_CURRENT_SOURCE_DIR}/algorithms/*.h
        ${CMAKE_CURRENT_SOURCE_DIR}/apps/*.h
        ${CMAKE_CURRENT_SOURCE_DIR}/apps/*.inl
        ${CMAKE_CURRENT_SOURCE_DIR}/common/*.h
        ${CMAKE_CURRENT_SOURCE_DIR}/common/*.inl
        ${CMAKE_CURRENT_SOURCE_DIR}/geometry/*.h
        ${CMAKE_CURRENT_SOURCE_DIR}/geometry/*.inl
        ${CMAKE_CURRENT_SOURCE_DIR}/blas/*.h
        ${CMAKE_CURRENT_SOURCE_DIR}/blas/*.inl
        ${CMAKE_CURRENT_SOURCE_DIR}/physics/*.h
        ${CMAKE_CURRENT_SOURCE_DIR}/physics/*.inl
        ${CMAKE_CURRENT_SOURCE_DIR}/structures/*.h
        ${CMAKE_CURRENT_SOURCE_DIR}/structures/*.inl
        )

file(GLOB SOURCES
        ${CMAKE_CURRENT_SOURCE_DIR}/animation/*.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/algorithms/*.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/apps/pic_custom.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/geometry/*.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/blas/*.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/common/*.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/physics/*.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/structures/*.cpp
        )

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
add_library(furoo STATIC ${SOURCES} ${HEADERS})
target_link_libraries(furoo PUBLIC OpenMP::OpenMP_CXX)
set_target_properties(furoo PROPERTIES
        OUTPUT_NAME "furoo"
        FOLDER "FUROO")
