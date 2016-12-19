cmake_minimum_required(VERSION 2.8)
project(TDPIC)

# options
option(BUILD_FOR_HOST "Build for the login node" ON)
option(COPY_SAMPLE_INPUT "Copy input.json to build dir" OFF)
option(DEBUG          "Build in Debug Mode"      OFF)
option(WITH_TEST      "Build with Test"          ON)

# add submodules
include_directories(${CMAKE_SOURCE_DIR}/picojson)
include_directories(${CMAKE_SOURCE_DIR}/include)

# check host
message(STATUS "CMAKE_SYSTEM_NAME = ${CMAKE_SYSTEM_NAME}")
if(${CMAKE_SYSTEM_NAME} STREQUAL "Windows")
    set(HOME "$ENV{USERPROFILE}")
    set(CMAKE_SUPPRESS_REGENERATION true)
elseif(${CMAKE_SYSTEM_NAME} STREQUAL "Darwin")
    set(HOME "$ENV{HOME}")
    find_package( Boost REQUIRED )
    include_directories( ${Boost_INCLUDE_DIRS} )

    # Doxygen
    FIND_PACKAGE(Doxygen)
    IF(DOXYGEN_FOUND)
        SET(DOXYGEN_CONF_FILE "doxygen.conf")
        SET(DOXYGEN_CONF_PATH ${CMAKE_CURRENT_BINARY_DIR}/${DOXYGEN_CONF_FILE})
        SET(DOXYGEN_SOURCE_DIR ${CMAKE_SOURCE_DIR})
        SET(DOXYGEN_TARGET "doc" )

        CONFIGURE_FILE(${DOXYGEN_CONF_FILE}.in ${DOXYGEN_CONF_PATH})

        ADD_CUSTOM_TARGET(${DOXYGEN_TARGET}
            ${DOXYGEN_EXECUTABLE} ${DOXYGEN_CONF_PATH}
            DEPENDS  ${DOXYGEN_CONF_PATH})
    ELSE(DOXYGEN_FOUND)
        MESSAGE (WARNING "Doxygen binary couldn't be found")
    ENDIF(DOXYGEN_FOUND)
elseif(${CMAKE_SYSTEM_NAME} STREQUAL "Linux")
    set(HOME "$ENV{HOME}")
endif()
message(STATUS "HOME              = ${HOME}")

# if(${CMAKE_SYSTEM_NAME} STREQUAL "CrayLinuxEnvironment")
if(${CMAKE_SYSTEM_NAME} STREQUAL "Linux")

    #set(CRAYPE_LINK_TYPE "dynamic")

    if(${DEBUG})
        set(CMAKE_BUILD_TYPE Debug)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0 -g -Wall -mcmodel=medium -fpic -dynamic -shared-intel -mkl")
        add_definitions(-DDEBUG)
    else()
        set(CMAKE_BUILD_TYPE Release)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -mcmodel=medium -fpic -dynamic -shared-intel -mkl")
    endif()

    if(BUILD_FOR_HOST)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -xHost")
    endif()

else()
    # Flags for Local build
    if(${DEBUG})
        set(CMAKE_BUILD_TYPE Debug)
        set(CMAKE_CXX_FLAGS "-O0 -g -Wall -m64")
        add_definitions(-DDEBUG)
    else()
        set(CMAKE_BUILD_TYPE Release)
        set(CMAKE_CXX_FLAGS "-O2 -m64")
    endif()

    # options for mkl
    set(MKLROOT "/opt/intel/mkl")
    set(MKL_INCLUDE_DIRS "${MKLROOT}/include")
    set(MKL_LIBRARIES "${MKLROOT}/lib")
    # add_definitions(-DMKL_LP64)
    link_directories(${MKL_LIBRARIES})
    include_directories(${MKL_INCLUDE_DIRS})
    find_library(MKL_CORE mkl_core HINTS ${MKL_LIBRARIES})
    find_library(MKL_INTEL_LP64 mkl_intel_lp64 HINTS ${MKL_LIBRARIES})
    find_library(MKL_INTEL_THREAD mkl_intel_thread HINTS ${MKL_LIBRARIES})
    find_library(MKL_BLAS95 mkl_blas95_lp64 HINTS ${MKL_LIBRARIES})
    find_library(MKL_SEQ mkl_sequential HINTS ${MKL_LIBRARIES})
    find_library(MKL_LAPACK mkl_lapack95_lp64 HINTS ${MKL_LIBRARIES})
endif()

# -- check C++11 --
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)

if(COMPILER_SUPPORTS_CXX11)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
elseif(COMPILER_SUPPORTS_CXX0X)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
else()
    message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support.")
endif()
# -----------------

# add include
INCLUDE_DIRECTORIES("${HOME}/local/include")
file(GLOB SRCS src/*.cpp)
message("SOURCES: ${SRCS}")

add_executable(tdpic ${SRCS})
target_link_libraries(tdpic ${Boost_LIBRARIES})

if(NOT (${CMAKE_SYSTEM_NAME} STREQUAL "Linux") )
    target_link_libraries(tdpic ${MKL_INTEL_LP64} ${MKL_BLAS95} ${MKL_LAPACK} ${MKL_SEQ} ${MKL_CORE})
endif()

if(${COPY_SAMPLE_INPUT})
    add_custom_command(
        OUTPUT copy_json
        COMMAND cp "${CMAKE_CURRENT_SOURCE_DIR}/src/input.json.sample" "${CMAKE_BINARY_DIR}/input.json"
    )

    add_custom_target(copy_sample_input DEPENDS copy_json)
    add_dependencies(tdpic copy_sample_input)
endif()

# compiler info
message(STATUS "CMAKE_CXX_COMPILER = ${CMAKE_CXX_COMPILER}")
message(STATUS "CMAKE_CXX_FLAGS = ${CMAKE_CXX_FLAGS}")

if(${WITH_TEST})
    enable_testing()
    add_subdirectory(${CMAKE_SOURCE_DIR}/googletest)
    include_directories(${CMAKE_SOURCE_DIR}/googletest/googletest/include)
    add_subdirectory(test)
endif()
