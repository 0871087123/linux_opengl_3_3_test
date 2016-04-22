cmake_minimum_required(VERSION 2.7)
find_package(PkgConfig REQUIRED)
pkg_search_module(GLFW REQUIRED glfw3)
pkg_search_module(GLEW REQUIRED glew)
pkg_search_module(ASSIMP REQUIRED assimp)
 include_directories(${GLFW_INCLUDE_DIRS})
 include_directories(${GLEW_INCLUDE_DIRS})
include_directories(${ASSIMP_INCLUDE_DIRS})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g --std=c++11")
link_libraries(GL ${GLEW_STATIC_LIBRARIES}  ${GLFW_STATIC_LIBRARIES} dl X11 ${ASSIMP_STATIC_LIBRARIES})

add_executable(gltest main.cpp)
