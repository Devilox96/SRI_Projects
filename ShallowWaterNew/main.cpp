#include <iostream>
//-----------------------------//
#include "Solver/Solver.h"
#include "Solver/OutputDataFormat.h"
//-----------------------------//
int main() {
//    auto TestSolver = new Solver;
////    auto DataFormat = new OutputDataFormat <double> ("Output.dat");
//
//    for (int i = 0; i < 50; i++) {
//        TestSolver -> Calc(10);
//        TestSolver -> SaveData("Output_" + std::to_string(i) + ".dat");
//        std::cout << "Step: " << i * 10 << std::endl;
//    }

//    std::vector <std::vector <dVectorND <double>>> DataArray;

    auto MHDSolver = new dRichtmyerSolver2D(0.0005, 0.1, 0.1);
    MHDSolver -> SetInitialData(200, 200);

    for (int i = 0; i < 3000; i++) {
        MHDSolver -> Solve();

        if (i % 10 == 0) {
            MHDSolver -> SaveData("Output_" + std::to_string(i / 10) + ".dat");
            std::cout << "Step: " << i << std::endl;

            std::cout << "Full energy: " << MHDSolver -> GetFullEnergy() << std::endl;
        }
    }

    return 0;
}