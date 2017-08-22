cmake_minimum_required(VERSION 2.8)
project(TDPIC)

# options
option(BUILD_FOR_HOST    "Build for the login node."            ON)
option(COPY_DEBUG_MPIRUN_SCRIPT   "Copy gdb_mpirun.sh to build dir." ON)
option(COPY_SAMPLE_INPUT "Copy input.json to build dir."        ON)
option(DEBUG             "Build in Debug Mode."                 ON)
option(CHARGE_CONSERVATION "Check charge conservation on each step." ON)
option(WITH_TEST         "Build Test code."                     ON)

set(CXX_REQUIRED_LIBS "")
execute_process (COMMAND date +%Y%m%d OUTPUT_VARIABLE CMAKE_DATE OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process (COMMAND git rev-parse --short HEAD OUTPUT_VARIABLE CMAKE_REVISION OUTPUT_STRIP_TRAILING_WHITESPACE)

set(TDPIC_VERSION 0.1.2)
set(TDPIC_DATE ${CMAKE_DATE})
set(TDPIC_REVISION ${CMAKE_REVISION})

configure_file(${CMAKE_SOURCE_DIR}/include/tdpic_configure.h.in ${CMAKE_SOURCE_DIR}/include/tdpic_configure.h @ONLY)

# add submodules
include_directories(${CMAKE_SOURCE_DIR}/picojson)
include_directories(${CMAKE_SOURCE_DIR}/SimpleXdmf/include)
include_directories(${CMAKE_SOURCE_DIR}/HighFive/include)
include_directories(${CMAKE_SOURCE_DIR}/include)

# check host
message(STATUS "CMAKE_SYSTEM_NAME = ${CMAKE_SYSTEM_NAME}")
if(${CMAKE_SYSTEM_NAME} STREQUAL "Windows")
    set(HOME "$ENV{USERPROFILE}")
    set(CMAKE_SUPPRESS_REGENERATION true)
    set(BOOST_USE_STATIC_LIBS ON)

    find_package( Boost REQUIRED COMPONENTS filesystem system)
    include_directories( ${Boost_INCLUDE_DIRS} )
    set(CXX_REQUIRED_LIBS ${CXX_REQUIRED_LIBS} ${Boost_LIBRARIES})
    message(STATUS "Boost: ${Boost_LIBRARIES}")

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
            ${DOXYGEN_EXECUTABLE} ${DOXYGEN_CONF_PATH} && open ${CMAKE_BINARY_DIR}/doc/html/index.html
            DEPENDS  ${DOXYGEN_CONF_PATH})
    ELSE(DOXYGEN_FOUND)
        MESSAGE (WARNING "Doxygen binary couldn't be found")
    ENDIF(DOXYGEN_FOUND)
elseif(${CMAKE_SYSTEM_NAME} STREQUAL "Linux")
    set(HOME "$ENV{HOME}")
    find_package( Boost REQUIRED COMPONENTS filesystem system)
    include_directories( ${Boost_INCLUDE_DIRS} )
    set(CXX_REQUIRED_LIBS ${CXX_REQUIRED_LIBS} ${Boost_LIBRARIES})
    message(STATUS "Boost: ${Boost_LIBRARIES}")
endif()
message(STATUS "HOME              = ${HOME}")

set(CXX_WARNING_FLAGS "-Wimplicit -Wreturn-type -Wswitch -Wcomment -Wformat=2 -Wuninitialized -Wcast-qual -Wcast-align -Wwrite-strings -Wpointer-arith -Winit-self -Wno-unused-local-typedefs")

if(${CMAKE_SYSTEM_NAME} STREQUAL "Linux")
    #set(CRAYPE_LINK_TYPE "dynamic")
    set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g3 -C -traceback ${CXX_WARNING_FLAGS} -mcmodel=medium -fpic -dynamic -shared-intel -mkl")
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

    message(STATUS "CXX_REQUIRED_LIBS: ${CXX_REQUIRED_LIBS}")
elseif(${CMAKE_SYSTEM_NAME} STREQUAL "Windows")

    # MS MPI
    # set(CMAKE_C_COMPILER mpicc)
    # set(CMAKE_CXX_COMPILER mpicxx)
    set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g3 ${CXX_WARNING_FLAGS} -m64")
    set(CMAKE_CXX_FLAGS_RELEASE "-O2 -m64")

    # MS MPI Dirs
    set(MSMPI_INCDIR "$ENV{MSMPI_INC}")
    set(MSMPI_LIBDIR "$ENV{MSMPI_LIB64}")
    link_directories(${MSMPI_LIBDIR})
    include_directories(${MSMPI_INCDIR})
    find_library(LIB_MSMPI msmpi HINTS ${MSMPI_LIBDIR})
    set(CXX_REQUIRED_LIBS ${CXX_REQUIRED_LIBS} ${LIB_MSMPI})

    if(${DEBUG})
        set(CMAKE_BUILD_TYPE Debug)
        add_definitions(-DDEBUG)
    else()
        set(CMAKE_BUILD_TYPE Release)
    endif()

    message(STATUS "CXX_REQUIRED_LIBS: ${CXX_REQUIRED_LIBS}")
else()
    # Flags for Darwin build

    # OpenMPI
    set(CMAKE_C_COMPILER /usr/local/bin/mpicc)
    set(CMAKE_CXX_COMPILER /usr/local/bin/mpicxx)
    set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g3 ${CXX_WARNING_FLAGS} -m64")
    set(CMAKE_CXX_FLAGS_RELEASE "-O2 -m64")

    if(${DEBUG})
        set(CMAKE_BUILD_TYPE Debug)
        add_definitions(-DDEBUG)
    else()
        set(CMAKE_BUILD_TYPE Release)
    endif()

    # MPI LIBRARIES
    find_library(MPI_LIBS mpi HINTS /usr/local/lib)
    set(CXX_REQUIRED_LIBS ${CXX_REQUIRED_LIBS} ${MPI_LIBS})

    # options for silo
    set(HOME_ROOT "${HOME}/local/")
    set(HOME_LIBRARIES "${HOME_ROOT}/lib")
    find_library(LIB_HDF5 hdf5 HINTS ${HOME_LIBRARIES})

    set(CXX_REQUIRED_LIBS ${CXX_REQUIRED_LIBS} ${LIB_Z} ${LIB_HDF5})
    message(STATUS "CXX_REQUIRED_LIBS: ${CXX_REQUIRED_LIBS}")
endif()

# enable computation validation
if(${CHARGE_CONSERVATION})
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DCHARGE_CONSERVATION")
endif()

# -- check C++11 --
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++14" COMPILER_SUPPORTS_CXX14)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)

if(COMPILER_SUPPORTS_CXX14)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
elseif(COMPILER_SUPPORTS_CXX11)
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

if(${COPY_DEBUG_MPIRUN_SCRIPT})
    add_custom_command(
        OUTPUT copy_debug_mpirun_command
        COMMAND cp "${CMAKE_CURRENT_SOURCE_DIR}/scripts/lldb_mpirun.sh" "${CMAKE_BINARY_DIR}/lldb_mpirun.sh"
    )

    add_custom_target(copy_debug_mpirun DEPENDS copy_debug_mpirun_command)
    add_dependencies(tdpic copy_debug_mpirun)
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
if(${DEBUG})
    message(STATUS "CMAKE_CXX_FLAGS = ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG}")
else()
    message(STATUS "CMAKE_CXX_FLAGS = ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE}")
endif()

if(${WITH_TEST})
    enable_testing()
    add_subdirectory(${CMAKE_SOURCE_DIR}/googletest)
    include_directories(${CMAKE_SOURCE_DIR}/googletest/googletest/include)
    add_subdirectory(test)
endif()
