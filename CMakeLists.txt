cmake_minimum_required(VERSION 2.8)
project(2DPIC)
#set(CMAKE_SYSTEM_NAME CrayLinuxEnvironment)
#set(CRAYPE_LINK_TYPE "dynamic")

# options
option(BUILD_FOR_HOST "Build for the login node" ON)

# check host
message("CMAKE_SYSTEM_NAME = ${CMAKE_SYSTEM_NAME}")
if(${CMAKE_SYSTEM_NAME} STREQUAL "Windows")
    set(HOME "$ENV{USERPROFILE}")
    set(CMAKE_SUPPRESS_REGENERATION true)
endif()
message("HOME              = ${HOME}")

if(${CMAKE_SYSTEM_NAME} STREQUAL "CrayLinuxEnvironment")
    #set(CMAKE_C_COMPILER "cc")
    #set(CMAKE_C_FLAGS "-O2 -mcmodel=medium -fpic -dynamic -shared-intel -mkl")
    set(CMAKE_CXX_COMPILER "CC")
    set(CMAKE_CXX_FLAGS "-O2 -mcmodel=medium -fpic -dynamic -shared-intel -mkl")
    #set(CMAKE_Fortran_COMPILER "ftn")
    #set(CMAKE_Fortran_FLAGS "-O2 -mcmodel=medium -fpic -dynamic -shared-intel -mkl")

    if(BUILD_FOR_HOST)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -xHost")
    endif()
else()
    # Flags for Local build
    set(CMAKE_CXX_FLAGS "-O2")
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

# Eigen3
#INCLUDE_DIRECTORIES("${HOME}/local/include/eigen3")

# Boost
INCLUDE_DIRECTORIES("${HOME}/local/include")

set(SRCS
    "src/main.cpp"
    "src/2dpic.cpp"
    "src/initialize.cpp"
    "src/utils.cpp"
)

add_executable(2dpic ${SRCS})
target_link_libraries(2dpic ${Boost_LIBRARIES})
