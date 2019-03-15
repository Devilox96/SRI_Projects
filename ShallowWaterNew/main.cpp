#include <iostream>
//-----------------------------//
#include "Solver/Solver.h"
#include "Solver/OutputDataFormat.h"
//-----------------------------//
int main() {
    auto MHDSolver = new dRichtmyerSolver2D(0.00005, 0.1, 0.1);
    MHDSolver -> SetSavePath("./Data/");
    MHDSolver -> SetInitialData(200, 200);
    double InitialFullEnergy = MHDSolver -> GetFullEnergy();

    for (int i = 0; i < 20000; i++) {
        MHDSolver -> SetTimeStep(0.2 * (MHDSolver -> GetXStep()) / MHDSolver -> GetMaxAmplitude());
        MHDSolver -> Solve();
        std::cout << 0.2 * (MHDSolver -> GetXStep()) / MHDSolver -> GetMaxAmplitude() << std::endl;

        if (i % 10 == 0) {
            MHDSolver -> AppendData();
            std::cout << "Step: " << i << std::endl;


            double FullEnergy = MHDSolver -> GetFullEnergy();

            std::cout << "Full energy: " << FullEnergy << std::endl;

            if (FullEnergy > InitialFullEnergy * 1.001) {
                break;
            }
        }
    }

    MHDSolver -> SaveData();

    system("python3.6 Plotting.py 200 200 -100 100 -100 100");

    return 0;
}