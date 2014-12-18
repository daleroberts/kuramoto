cmake_minimum_required(VERSION 2.8.4)
project(kuramoto)

find_package(PkgConfig)
pkg_check_modules(EIGEN3 REQUIRED eigen3)
include_directories(${EIGEN3_INCLUDE_DIRS})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DDEBUG")

set(SOURCE_FILES kuramoto.cc graph.cc)
set(HEADER_FILES graph.h utils.h variates.h)
source_group("Header Files" FILES ${HEADER_FILES})

add_executable(kuramoto ${SOURCE_FILES})
