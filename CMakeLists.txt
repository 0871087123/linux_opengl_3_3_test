cmake_minimum_required(VERSION 2.7)
find_package(PkgConfig REQUIRED)
pkg_search_module(GLFW REQUIRED glfw3)
pkg_search_module(GLEW REQUIRED glew)
 include_directories(${GLFW_INCLUDE_DIRS})
 include_directories(${GLEW_INCLUDE_DIRS})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g --std=c++11")
link_libraries(GL ${GLEW_STATIC_LIBRARIES}  ${GLFW_STATIC_LIBRARIES} dl X11 )

add_executable(gltest main.cpp)
