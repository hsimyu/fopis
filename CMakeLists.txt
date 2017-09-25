cmake_minimum_required(VERSION 2.8)
project(TDPIC)

# options
option(COPY_SAMPLE_INPUT "Copy input.json to build dir." ON)
option(DEBUG "Build in Debug Mode." ON)
option(CHARGE_CONSERVATION "Check charge conservation on each step." ON)
option(WITH_TEST "Build Test code." OFF)
set(LOCAL_LIBRARYDIR "NOTFOUND" CACHE PATH "Specified local library path to include and link.")

execute_process (COMMAND date +%Y%m%d OUTPUT_VARIABLE CMAKE_DATE OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process (COMMAND git rev-parse --short HEAD OUTPUT_VARIABLE CMAKE_REVISION OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process (COMMAND git describe --abbrev=0 OUTPUT_VARIABLE CMAKE_TAG OUTPUT_STRIP_TRAILING_WHITESPACE)

set(TDPIC_VERSION ${CMAKE_TAG})
set(TDPIC_DATE ${CMAKE_DATE})
set(TDPIC_REVISION ${CMAKE_REVISION})

configure_file(${CMAKE_SOURCE_DIR}/include/tdpic_configure.h.in ${CMAKE_SOURCE_DIR}/include/tdpic_configure.h @ONLY)

set(CXX_REQUIRED_LIBS "")
# add local library
if(LOCAL_LIBRARYDIR)
    message(STATUS "LOCAL_LIBRARYDIR = ${LOCAL_LIBRARYDIR}")
    link_directories(${LOCAL_LIBRARYDIR}/lib)
    include_directories(${LOCAL_LIBRARYDIR}/include)
    set(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH} ${LOCAL_LIBRARYDIR}/cmake ${LOCAL_LIBRARYDIR}/share/cmake)
    message(STATUS "CMAKE_PREFIX_PATH = ${CMAKE_PREFIX_PATH}")
endif()

# add include and submodules
include_directories(${CMAKE_SOURCE_DIR}/picojson)
include_directories(${CMAKE_SOURCE_DIR}/SimpleVTK/include)
include_directories(${CMAKE_SOURCE_DIR}/HighFive/include)
include_directories(${CMAKE_SOURCE_DIR}/include)

# check host
message(STATUS "CMAKE_SYSTEM_NAME = ${CMAKE_SYSTEM_NAME}")
if(${CMAKE_SYSTEM_NAME} STREQUAL "Windows")
    set(HOME "$ENV{USERPROFILE}")
elseif(${CMAKE_SYSTEM_NAME} STREQUAL "Darwin")
    set(HOME "$ENV{HOME}")
elseif(${CMAKE_SYSTEM_NAME} STREQUAL "Linux" OR ${CMAKE_SYSTEM_NAME} STREQUAL "CrayLinuxEnvironment")
    set(HOME "$ENV{HOME}")
endif()
message(STATUS "HOME PATH = ${HOME}")

# Boost
set(Boost_USE_STATIC_LIBS ON)
find_package(Boost REQUIRED COMPONENTS filesystem system)
link_directories( ${Boost_LIBRARY_DIRS} )
include_directories( ${Boost_INCLUDE_DIRS} )
set(CXX_REQUIRED_LIBS ${CXX_REQUIRED_LIBS} ${Boost_LIBRARIES})
message(STATUS "Boost LIBRARY DIR: ${Boost_LIBRARY_DIRS}")
message(STATUS "Boost INCLUDE DIR: ${Boost_INCLUDE_DIRS}")
message(STATUS "Boost LIBRARIES  : ${Boost_LIBRARIES}")

# HDF5 and zlib
set (LIB_TYPE STATIC)
string(TOLOWER ${LIB_TYPE} SEARCH_TYPE)

if(${CMAKE_SYSTEM_NAME} STREQUAL "Linux" OR ${CMAKE_SYSTEM_NAME} STREQUAL "Darwin")
    set(HDF5_LIBRARY "${LOCAL_LIBRARYDIR}/lib/libhdf5.a;${LOCAL_LIBRARYDIR}/lib/libz.a")
    set(CXX_REQUIRED_LIBS ${CXX_REQUIRED_LIBS} ${HDF5_LIBRARY})
elseif(${CMAKE_SYSTEM_NAME} STREQUAL "CrayLinuxEnvironment")
    link_directories( "$ENV{HDF5_DIR}/lib" )
    include_directories( "$ENV{HDF5_INCLUDE_OPTS}" )
    find_library(HDF5_LIBRARY NAMES hdf5 PATHS "$ENV{HDF5_DIR}/lib")
    find_library(ZLIB_LIBRARY NAMES z)
    set(CXX_REQUIRED_LIBS ${CXX_REQUIRED_LIBS} ${HDF5_LIBRARY} ${ZLIB_LIBRARY})
else()
    find_package (HDF5 NAMES hdf5 REQUIRED COMPONENTS C ${SEARCH_TYPE} HINTS ${LOCAL_LIBRARYDIR}/lib)
    link_directories(${HDF5_LIBRARY_DIRS}) # for zlib inclusion
    set(CXX_REQUIRED_LIBS ${CXX_REQUIRED_LIBS} ${HDF5_C_${LIB_TYPE}_LIBRARY})
endif()

# Doxygen (optional)
find_package(Doxygen)
IF(DOXYGEN_FOUND)
    SET(DOXYGEN_CONF_FILE "doxygen.conf")
    SET(DOXYGEN_CONF_PATH ${CMAKE_CURRENT_BINARY_DIR}/${DOXYGEN_CONF_FILE})
    SET(DOXYGEN_SOURCE_DIR ${CMAKE_SOURCE_DIR})
    SET(DOXYGEN_TARGET "doc" )

    CONFIGURE_FILE(${DOXYGEN_CONF_FILE}.in ${DOXYGEN_CONF_PATH})

    ADD_CUSTOM_TARGET(${DOXYGEN_TARGET}
        ${DOXYGEN_EXECUTABLE} ${DOXYGEN_CONF_PATH} && open ${CMAKE_BINARY_DIR}/doc/html/index.html
        DEPENDS  ${DOXYGEN_CONF_PATH})
ENDIF(DOXYGEN_FOUND)

if(${CMAKE_SYSTEM_NAME} STREQUAL "Linux" OR ${CMAKE_SYSTEM_NAME} STREQUAL "CrayLinuxEnvironment")
    set(CXX_WARNING_FLAGS "-Wimplicit -Wreturn-type -Wswitch -Wcomment -Wformat=2 -Wuninitialized -Wwrite-strings -Wpointer-arith -Winit-self")
    set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g3 -C -traceback ${CXX_WARNING_FLAGS} -mcmodel=medium -fpic -dynamic -shared-intel -mkl -qopenmp")
    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -mcmodel=medium -fpic -dynamic -shared-intel -mkl -qopenmp")
elseif(${CMAKE_SYSTEM_NAME} STREQUAL "Windows")
    set(CXX_WARNING_FLAGS "/W4")

    # MS MPI
    set(CMAKE_CXX_FLAGS_DEBUG "/MDd /GL- /Zi ${CXX_WARNING_FLAGS}")
    set(CMAKE_CXX_FLAGS_RELEASE "/MD /O2")

    # MS MPI Dirs
    set(MSMPI_INCDIR "$ENV{MSMPI_INC}")
    set(MSMPI_LIBDIR "$ENV{MSMPI_LIB64}")
    link_directories(${MSMPI_LIBDIR})
    include_directories(${MSMPI_INCDIR})
    find_library(LIB_MSMPI msmpi HINTS ${MSMPI_LIBDIR})
    set(CXX_REQUIRED_LIBS ${CXX_REQUIRED_LIBS} ${LIB_MSMPI})

    FILE(WRITE "${CMAKE_CURRENT_BINARY_DIR}/mpirun_tdpic.vcxproj.user"
        "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n"
        "<Project ToolsVersion=\"15.0\">\n"
        "<PropertyGroup Condition=\"'$(Configuration)|$(Platform)'=='Debug|x64'\">\n"
        "<LocalDebuggerCommand>$ENV{MSMPI_BIN}mpiexec.exe</LocalDebuggerCommand>\n"
        "<LocalDebuggerCommandArguments>-np 4 $(ProjectName).$(OutputType)</LocalDebuggerCommandArguments>\n"
        "<DebuggerFlavor>WindowsLocalDebugger</DebuggerFlavor>\n"
        "<LocalDebuggerWorkingDirectory>$(ProjectDir)$(Configuration)</LocalDebuggerWorkingDirectory>\n"
        "</PropertyGroup>\n"
        "</Project>\n")
else()
    # Flags for Darwin build
    set(CXX_WARNING_FLAGS "-Wimplicit -Wreturn-type -Wswitch -Wcomment -Wformat=2 -Wuninitialized -Wcast-align -Wwrite-strings -Wpointer-arith -Winit-self -Wno-unused-local-typedefs")

    # OpenMPI compilers
    set(CMAKE_CXX_COMPILER /usr/local/opt/llvm/bin/clang++)
    link_directories(/usr/local/opt/llvm/lib)

    # CXX Flags
    set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g3 ${CXX_WARNING_FLAGS} -m64 -fopenmp")
    set(CMAKE_CXX_FLAGS_RELEASE "-O2 -m64 -fopenmp")

    # MPI LIBRARIES
    find_library(MPI_LIBS mpi HINTS /usr/local/lib)
    set(CXX_REQUIRED_LIBS ${CXX_REQUIRED_LIBS} ${MPI_LIBS})
endif()
message(STATUS "CXX_REQUIRED_LIBS: ${CXX_REQUIRED_LIBS}")

if(${DEBUG})
    set(CMAKE_BUILD_TYPE Debug)
    add_definitions(-DDEBUG)
else()
    set(CMAKE_BUILD_TYPE Release)
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
    message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++14 support.")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
elseif(COMPILER_SUPPORTS_CXX0X)
    message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support.")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
else()
    message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++0x support.")
endif()
# -----------------
file(GLOB_RECURSE SRCS 
    ${PROJECT_SOURCE_DIR}/src/*.cpp
)
message(STATUS "Current source files: ${SRCS}")

add_executable(tdpic ${SRCS})
target_link_libraries(tdpic ${CXX_REQUIRED_LIBS})

install(
    TARGETS tdpic
    DESTINATION .
)

install(
    DIRECTORY input_samples
    DESTINATION .
)

add_executable(mpirun_tdpic EXCLUDE_FROM_ALL ${SRCS})
target_link_libraries(mpirun_tdpic ${CXX_REQUIRED_LIBS})

if(${COPY_SAMPLE_INPUT})
    if(${CMAKE_SYSTEM_NAME} STREQUAL "Windows")
        # add_custom_command(
        #     OUTPUT copy_json
        #     COMMAND copy "\"${CMAKE_CURRENT_SOURCE_DIR}/input_samples/input.json\"" "\"${CMAKE_BINARY_DIR}/input.json\""
        # )
        # add_custom_command(
        #     OUTPUT copy_objects_dir
        #     COMMAND copy "${CMAKE_CURRENT_SOURCE_DIR}\\input_samples\\objects" "${CMAKE_BINARY_DIR}\\"
        # )
        # add_custom_target(copy_sample_input DEPENDS copy_json)
        # add_dependencies(tdpic copy_sample_input)
    else()
        add_custom_command(
            OUTPUT copy_json
            COMMAND cp "${CMAKE_CURRENT_SOURCE_DIR}/input_samples/input.json" "${CMAKE_BINARY_DIR}/input.json"
        )

        add_custom_command(
            OUTPUT copy_objects_dir
            COMMAND cp -r "${CMAKE_CURRENT_SOURCE_DIR}/input_samples/objects" "${CMAKE_BINARY_DIR}/"
        )
        add_custom_target(copy_sample_input DEPENDS copy_json copy_objects_dir)
        add_dependencies(tdpic copy_sample_input)
    endif()
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
