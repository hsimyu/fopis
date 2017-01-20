set(GRID_TEST_SRCS
    "grid_test.cpp"
)

add_executable(grid_test ${SRCS} ${TEST_FIXTURE_SRCS} ${GRID_TEST_SRCS})

target_link_libraries(grid_test gtest gtest_main ${CXX_REQUIRED_LIBS})

add_test(
    NAME grid
    COMMAND $<TARGET_FILE:grid_test>
)
