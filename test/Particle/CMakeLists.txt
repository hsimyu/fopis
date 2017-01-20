set(PARTICLE_TEST_SRCS
    "particle_test.cpp"
)

add_executable(particle_test ${SRCS} ${TEST_FIXTURE_SRCS} ${PARTICLE_TEST_SRCS})

target_link_libraries(particle_test gtest gtest_main ${CXX_REQUIRED_LIBS})

add_test(
    NAME particle
    COMMAND $<TARGET_FILE:particle_test>
)
