#include <iostream>
//-----------------------------//
#include "Solver/Solver.h"
//-----------------------------//
int main() {
    auto MHDSolver = new dRichtmyerSolver2D(0.00005, 0.1, 0.1);
    MHDSolver -> SetSavePath("./Data/");
    MHDSolver -> SetInitialData(200, 200);
    double InitialFullEnergy = MHDSolver -> GetFullEnergy();
    double Time = 0.0;

    for (int i = 0; i < 20000; i++) {
        MHDSolver -> SetTimeStep(0.3 * (MHDSolver -> GetXStep()) / MHDSolver -> GetMaxAmplitude());
        MHDSolver -> Solve();
        std::cout << 0.3 * (MHDSolver -> GetXStep()) / MHDSolver -> GetMaxAmplitude() << std::endl;
        Time += MHDSolver -> GetTimeStep();

        if (i % 10 == 0) {
            MHDSolver -> AppendData();
            std::cout << "Step: " << i << "\tTime:" << Time << std::endl;


            double FullEnergy = MHDSolver -> GetFullEnergy();

            std::cout << "Full energy: " << FullEnergy << std::endl;

            if (FullEnergy > InitialFullEnergy * 1.1) {
                break;
            }
        }
    }

    MHDSolver -> SaveData();

    system("python3.6 Plotting.py 200 200 -100 100 -100 100");

    return 0;
}