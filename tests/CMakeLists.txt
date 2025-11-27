set(target furoo_unit_tests)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

file(GLOB SOURCES
        # ${CMAKE_CURRENT_SOURCE_DIR}/transform_test.cpp
        # ${CMAKE_CURRENT_SOURCE_DIR}/distribution_test.cpp
        # ${CMAKE_CURRENT_SOURCE_DIR}/rbf_test.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/point_z_grid_test.cpp
        )

include_directories("${FUROO_SOURCE_DIR}/src")
target_link_libraries(furoo)
add_executable(${target} ${SOURCES} main.cpp)

if (WIN32)
    target_link_libraries(${target} furoo gtest)
elseif (APPLE)
    target_link_libraries(${target} furoo gtest pthread)
else ()
    target_link_libraries(${target} furoo gtest gcov pthread)
endif (WIN32)

add_custom_target(tests
        COMMAND furoo_unit_tests DEPENDS furoo_unit_tests
        COMMAND lcov --directory .. --capture -output-file cov.info
        COMMAND lcov --remove cov.info '/usr/*' --output-file cov.info
        COMMAND lcov --remove cov.info '*/external/*' --output-file cov.info.cleaned
        COMMAND genhtml cov.info.cleaned --output-directory out
        )

add_subdirectory(simulations)
