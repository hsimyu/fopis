set(SRCS
    "particle_test.cpp"
)

add_executable(particle_test ${SRCS})

target_link_libraries(particle_test
    gtest
    gtest_main
    )

add_test(
    NAME particle
    COMMAND $<TARGET_FILE:particle_test>
)
