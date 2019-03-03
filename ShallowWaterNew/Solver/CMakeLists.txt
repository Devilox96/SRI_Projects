project(DolverLib)

set(Source
        Solver.cpp
        OutputDataFormat.cpp)

add_library(SolverLib ${Source})
target_link_libraries(SolverLib
        dMathCoreLib
        dMathNumerCalcLib)