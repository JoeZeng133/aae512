project(aae512)
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_STANDARD 14)

#However, the file(GLOB...) allows for wildcard additions:
file(GLOB SOURCES "./hw1.cpp")

add_executable(aae512 ${SOURCES})
