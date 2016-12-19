set(FIELD_TEST_SRCS
    "field_test.cpp"
)

add_executable(field_test ${SRCS} ${FIELD_TEST_SRCS})

target_link_libraries(field_test
    gtest
    gtest_main
    ${MKL_INTEL_LP64}
    ${MKL_BLAS95}
    ${MKL_LAPACK}
    ${MKL_SEQ}
    ${MKL_CORE}
    )

add_test(
    NAME field
    COMMAND $<TARGET_FILE:field_test>
)
