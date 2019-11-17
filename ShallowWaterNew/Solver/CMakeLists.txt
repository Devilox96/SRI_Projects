project(SolverLib)

set(Source
        Solver.cpp
        )

add_library(SolverLib ${Source})
target_link_libraries(SolverLib
        dMathCoreLib
        dMathNumerCalcLib)