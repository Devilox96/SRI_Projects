project(dMathLib)

add_subdirectory(Core)
add_subdirectory(NumerCalc)

set(Source
        dMath.h)

set(Libs
        dMathCoreLib
        dMathNumerCalcLib)

add_library(dMathLib ${Source})
target_link_libraries(dMathLib ${Libs})
