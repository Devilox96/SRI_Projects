cmake_minimum_required(VERSION 3.10)
project(ShallowWater)

set(CMAKE_CXX_STANDARD 11)

add_subdirectory(Include/dMath)

add_executable(ShallowWater
        main.cpp
        Include/Calculation/Calculation.cpp)

target_link_libraries(ShallowWater dMathLib)