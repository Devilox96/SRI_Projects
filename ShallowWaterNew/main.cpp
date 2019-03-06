#include <iostream>
//-----------------------------//
#include "Solver/Solver.h"
#include "Solver/OutputDataFormat.h"
//-----------------------------//
int main() {
    auto TestSolver = new Solver;
//    auto DataFormat = new OutputDataFormat <double> ("Output.dat");

    for (int i = 0; i < 50; i++) {
        TestSolver -> Calc(10);
        TestSolver -> SaveData("Output_" + std::to_string(i) + ".dat");
        std::cout << "Step: " << i * 10 << std::endl;
    }

    return 0;
}