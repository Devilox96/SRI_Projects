cmake_minimum_required(VERSION 3.13)
project(ShallowWaterNew)

set(CMAKE_CXX_STANDARD 14)

add_subdirectory(Libs/dMath/Core)
add_subdirectory(Libs/dMath/NumerCalc)

add_subdirectory(Solver)

add_executable(ShallowWaterNew main.cpp Solver/Solver.cpp Solver/Solver.h)
target_link_libraries(ShallowWaterNew
        dMathCoreLib
        dMathNumerCalcLib
        SolverLib)