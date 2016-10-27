cmake_minimum_required(VERSION 2.8)
project(TDPIC)

# options
option(BUILD_FOR_HOST "Build for the login node" ON)
option(DEBUG          "Build in Debug Mode"      OFF)

# check host
message(STATUS "CMAKE_CXX_FLAGS_DEBUG = ${CMAKE_CXX_FLAGS_DEBUG}")
message(STATUS "CMAKE_SYSTEM_NAME = ${CMAKE_SYSTEM_NAME}")
if(${CMAKE_SYSTEM_NAME} STREQUAL "Windows")
    set(HOME "$ENV{USERPROFILE}")
    set(CMAKE_SUPPRESS_REGENERATION true)
elseif(${CMAKE_SYSTEM_NAME} STREQUAL "Darwin")
    set(HOME "$ENV{HOME}")
    find_package( Boost REQUIRED )
    include_directories( ${Boost_INCLUDE_DIRS} )
elseif(${CMAKE_SYSTEM_NAME} STREQUAL "Linux")
    set(HOME "$ENV{HOME}")
endif()
message(STATUS "HOME              = ${HOME}")

# if(${CMAKE_SYSTEM_NAME} STREQUAL "CrayLinuxEnvironment")
if(${CMAKE_SYSTEM_NAME} STREQUAL "Linux")
    #set(CRAYPE_LINK_TYPE "dynamic")
    set(CMAKE_CXX_COMPILER "CC")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -mcmodel=medium -fpic -dynamic -shared-intel -mkl")
    if(BUILD_FOR_HOST)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -xHost")
    endif()
else()
    # Flags for Local build
    if(${DEBUG})
        set(CMAKE_BUILD_TYPE Debug)
        set(CMAKE_CXX_FLAGS "-O0 -g -Wall")
    else()
        set(CMAKE_BUILD_TYPE Release)
        set(CMAKE_CXX_FLAGS "-O2")
    endif()
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

set(SRCS
    "src/main.cpp"
    "src/particle.cpp"
    "src/position.cpp"
    "src/velocity.cpp"
    "src/field.cpp"
    "src/initialize.cpp"
    "src/utils.cpp"
)

add_executable(tdpic ${SRCS})
target_link_libraries(tdpic ${Boost_LIBRARIES})
