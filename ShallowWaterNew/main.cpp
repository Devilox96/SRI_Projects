#include <iostream>
//-----------------------------//
//#include "Solver/Solver.h"
//#include "Solver/RossbySolver.h"
#include "Solver/TestSolver.h"
//-----------------------------//
int main() {
//    auto MHDSolver = new dRichtmyerSolver2D(0.00005, 0.1, 0.1);
////    auto MHDSolver = new dLaxFriedrichsSolver2D(0.00005, 0.1, 0.1);
//    MHDSolver -> SetSavePath("./Data/");
//    MHDSolver -> SetInitialData(200, 200);
//    double InitialFullEnergy = MHDSolver -> GetFullEnergy();
//    double Time = 0.0;
//
//    for (int i = 0; i < 1000; i++) {
//        double NewTime = (MHDSolver -> getStepY()) / (MHDSolver -> getMaxSpeed() + 0.343) * 0.5;
//
//        MHDSolver -> setTimeStep(NewTime);
//        MHDSolver->solveGridRichtmyerZwas();
////        MHDSolver->solveGridRichtmyer();
////        MHDSolver->solveGrid();
////        std::cout << NewTime << std::endl;
//        Time += MHDSolver -> getStepTime();
//
//        if (i % 50 == 0) {
//            MHDSolver -> AppendData();
//            std::cout << "Step: " << i << "\tTime:" << Time << std::endl;
//
//
//            double FullEnergy = MHDSolver -> GetFullEnergy();
//
//            std::cout << "Full energy: " << FullEnergy << std::endl;
//
//            if (FullEnergy > InitialFullEnergy * 1.1) {
//                break;
//            }
//        }
//    }
//
//    MHDSolver -> SaveData();







//    auto Solver = new RossbySolver(0.00005, 111e+03, 111e+03);
//
//    Solver -> setSavePath("./Data/");
//    Solver -> initGrid(180, 88);
//    Solver -> solve();



    auto Solver = new TestSolver;
    Solver -> setSavePath("./Data/");
    Solver -> openFiles();
    Solver -> solveCustom();

//    system("python3.6 Plotting.py 200 200 -100 100 -100 100");

    return 0;
}