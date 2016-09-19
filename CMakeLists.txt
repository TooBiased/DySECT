cmake_minimum_required (VERSION 2.6)
project (space-grow-table)

set (CMAKE_CXX_FLAGS "-std=c++14 -msse4.2 -Wall -Wextra -O3")

include_directories(.)

add_executable(in source/in_test.cpp)
