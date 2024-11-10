cmake_minimum_required(VERSION 3.10)

project(tri-ben)

set (CMAKE_CXX_STANDARD 14)

include_directories(include)

file(GLOB SOURCES "src/*.cpp")

set(CMAKE_CXX_FLAGS " -Wextra -g -O3")

add_executable(tri-ben.exe ${SOURCES})

target_link_libraries(tri-ben.exe -lsfml-graphics -lsfml-window -lsfml-system)