project(aae512)
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_STANDARD 14)

find_package( Armadillo REQUIRED )
include_directories(${ARMADILLO_INCLUDE_DIRS})

#However, the file(GLOB...) allows for wildcard additions:
file(GLOB SOURCES "./hw2.cpp")
add_executable(aae512 ${SOURCES})
target_link_libraries(aae512 ${ARMADILLO_LIBRARIES})
