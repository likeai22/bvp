cmake_minimum_required(VERSION 3.12)
project(bvp)

set(CMAKE_CXX_STANDARD 14)

add_executable(bvp main.cpp BVP.cpp BVP.h functions_bvp.cpp functions_bvp.h test.cpp test.h)