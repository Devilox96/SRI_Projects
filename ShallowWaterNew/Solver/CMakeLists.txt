project(SolverLib)

set(Source
        Solver.cpp
        RossbySolver.cpp)

add_library(SolverLib ${Source})
target_link_libraries(SolverLib
        dMathCoreLib
        dMathNumerCalcLib)