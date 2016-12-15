set(SRCS
    "particle_test.cpp"
)

add_executable(particle_test ${SRCS})

target_link_libraries(particle_test
    gtest
    gtest_main
    ${MKL_INTEL_LP64}
    ${MKL_CORE}
    )

add_test(
    NAME particle
    COMMAND $<TARGET_FILE:particle_test>
)
