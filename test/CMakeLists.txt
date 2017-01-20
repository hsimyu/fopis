file(GLOB TEST_FIXTURE_SRCS *.cpp)
message(STATUS "TEST_FIXTURE_SRCS = ${TEST_FIXTURE_SRCS}")
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

set(TEST_SRCS
    "Particle"
    "Field"
    "Grid"
)

add_definitions(-DBUILD_TEST)
subdirs(${TEST_SRCS})
