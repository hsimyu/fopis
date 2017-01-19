cmake_minimum_required(VERSION 2.8)
project(TDPIC)

# options
option(BUILD_FOR_HOST    "Build for the login node"      ON)
option(COPY_SAMPLE_INPUT "Copy input.json to build dir"  ON)
option(DEBUG             "Build in Debug Mode"           ON)

set(CXX_REQUIRED_LIBS "")

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
    find_package( Boost REQUIRED COMPONENTS filesystem system)
    include_directories( ${Boost_INCLUDE_DIRS} )
    set(CXX_REQUIRED_LIBS ${CXX_REQUIRED_LIBS} ${Boost_LIBRARIES})
    message(STATUS "Boost: ${Boost_LIBRARIES}")

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

set(CXX_WARNING_FLAGS "-Wimplicit -Wreturn-type -Wswitch -Wcomment -Wformat=2 -Wuninitialized -Wcast-qual -Wcast-align -Wwrite-strings -Wpointer-arith -Winit-self")

# if(${CMAKE_SYSTEM_NAME} STREQUAL "CrayLinuxEnvironment")
if(${CMAKE_SYSTEM_NAME} STREQUAL "Linux")

    #set(CRAYPE_LINK_TYPE "dynamic")
    set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g ${CXX_WARNING_FLAGS} -mcmodel=medium -fpic -dynamic -shared-intel -mkl")
    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -mcmodel=medium -fpic -dynamic -shared-intel -mkl")

    if(${DEBUG})
        set(CMAKE_BUILD_TYPE Debug)
        add_definitions(-DDEBUG)
    else()
        set(CMAKE_BUILD_TYPE Release)
    endif()

    if(BUILD_FOR_HOST)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -xHost")
    endif()

else()
    # Flags for Local build

    # OpenMPI
    set(CMAKE_C_COMPILER mpicc)
    set(CMAKE_CXX_COMPILER mpicxx)
    set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g ${CXX_WARNING_FLAGS} -m64")
    set(CMAKE_CXX_FLAGS_RELEASE "-O2 -m64")

    if(${DEBUG})
        set(CMAKE_BUILD_TYPE Debug)
        add_definitions(-DDEBUG)
    else()
        set(CMAKE_BUILD_TYPE Release)
    endif()

    # options for mkl
    set(MKLROOT "/opt/intel/mkl")
    set(MKL_INCLUDE_DIRS "${MKLROOT}/include")
    set(MKL_LIBRARIES "${MKLROOT}/lib")
    # add_definitions(-DMKL_ILP64)
    link_directories(${MKL_LIBRARIES})
    include_directories(${MKL_INCLUDE_DIRS})
    find_library(MKL_CORE mkl_core HINTS ${MKL_LIBRARIES})
    find_library(MKL_INTEL_LP64 mkl_intel_lp64 HINTS ${MKL_LIBRARIES})
    find_library(MKL_INTEL_THREAD mkl_intel_thread HINTS ${MKL_LIBRARIES})
    find_library(MKL_BLAS95 mkl_blas95_lp64 HINTS ${MKL_LIBRARIES})
    find_library(MKL_SEQ mkl_sequential HINTS ${MKL_LIBRARIES})
    find_library(MKL_LAPACK mkl_lapack95_lp64 HINTS ${MKL_LIBRARIES})

    set(CXX_REQUIRED_LIBS ${CXX_REQUIRED_LIBS} ${MKL_INTEL_LP64} ${MKL_BLAS95} ${MKL_LAPACK} ${MKL_SEQ} ${MKL_CORE})

    # options for hdf5
    set(HOME_ROOT "${HOME}/local/")
    set(HOME_LIBRARIES "${HOME_ROOT}/lib")
    find_library(LIB_SILO silo HINTS ${HOME_LIBRARIES})

    set(CXX_REQUIRED_LIBS ${CXX_REQUIRED_LIBS} ${LIB_SILO})
    message(STATUS "CXX_REQUIRED_LIBS: ${CXX_REQUIRED_LIBS}")
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

add_executable(tdpic ${SRCS})
target_link_libraries(tdpic ${CXX_REQUIRED_LIBS})

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
if(${DEBUG})
    message(STATUS "CMAKE_CXX_FLAGS = ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG}")
else()
    message(STATUS "CMAKE_CXX_FLAGS = ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE}")
endif()
