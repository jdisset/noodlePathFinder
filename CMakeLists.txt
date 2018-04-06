cmake_minimum_required(VERSION 3.1)
project(NoodlePathFinder CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "-O3 -g -Wall -Wextra -pedantic ")

add_executable(noodlepath src/utils.hpp src/vec.hpp src/main.cpp)

add_executable(launchtest src/utils.hpp src/vec.hpp ../test/main.cpp)
