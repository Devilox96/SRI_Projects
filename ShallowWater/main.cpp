#include <iostream>
#include <vector>
//-----------------------------//
#include "Include/dMath/dMath.h"
#include "Include/Calculation/Calculation.h"
//-----------------------------//
int main() {
//    Calculation Method(50, 50);

    dVectorND <int> Test1({1, 2, 3, 4, 5});
    dVectorND <int> Test2({1, 2, 3, 4, 5});
    dVectorND <int> Test3({1, 6, 3, 4, 5});

    std::cout << (Test1 == Test2) << std::endl;
    std::cout << (Test1 == Test3) << std::endl;

    Test1 = Test3;

    std::cout << (Test1 == Test2) << std::endl;
    std::cout << (Test1 == Test3) << std::endl;

    return 0;
}