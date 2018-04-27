cmake_minimum_required (VERSION 3.9)
project(main)
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_STANDARD 11)

include_directories(headers)

file(GLOB SOURCES "src/*.cpp")
add_executable(main ${SOURCES})
