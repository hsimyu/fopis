set(FIELD_TEST_SRCS
    "field_test.cpp"
)

add_executable(field_test ${SRCS} ${FIELD_TEST_SRCS})

target_link_libraries(field_test gtest gtest_main ${CXX_REQUIRED_LIBS})

add_test(
    NAME field
    COMMAND $<TARGET_FILE:field_test>
)
