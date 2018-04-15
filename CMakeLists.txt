cmake_minimum_required (VERSION 3.9)
project(mesh)
set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_CXX_STANDARD 11)

include_directories(headers)

file(GLOB SOURCES "src/*.cpp")
add_executable(mesh ${SOURCES})
