set(TEST_SRCS
    "Particle"
    "Field"
    )

add_definitions(-DBUILD_TEST)
subdirs(${TEST_SRCS})
