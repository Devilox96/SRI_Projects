#include "Solver.h"
//-----------------------------//
Solver::Solver() {
    Test = new dRichtmyerSolver2D(0.0005, 0.1, 0.1);

    InitGrid(10.0, 0.0, 0.0, 200);

    Grid[100][100][0] = 12.0;
}
Solver::~Solver() {
    delete Test;
}

void Solver :: Calc(unsigned int StepsP) {
    std::ofstream EnergyFile;
    EnergyFile.open("Energy.dat", std::ios::app);


    unsigned int PointNum = 200;

    std :: vector <std :: vector <dVectorND <double>>> TempGrid = Grid;

    for (unsigned int iter = 0; iter < StepsP; iter++) {
        double Energy = 0.0;

        for (unsigned int i = 0; i < PointNum; i++) {
            for (unsigned int j = 0; j < PointNum; j++) {
                dVectorND <double> TempH = Test -> Solve(Grid, i, j, 200, 200);

                TempGrid[i][j] = TempH;

                Energy += ( pow((Grid[i][j][0] - TempGrid[i][j][0]) / Test -> GetTimeStep(), 2.0) +
                            pow(TempGrid[i][j][1] / TempGrid[i][j][0], 2.0) +
                            pow(TempGrid[i][j][2] / TempGrid[i][j][0], 2.0));
            }
        }

        for (unsigned int i = 0; i < PointNum; i++) {
            Energy += ( pow((Grid[i][0][0] - TempGrid[i][0][0]) / Test -> GetTimeStep(), 1.0) +
                        pow(TempGrid[i][0][1] / TempGrid[i][0][0], 2.0) +
                        pow(TempGrid[i][0][2] / TempGrid[i][0][0], 2.0));
        }
        for (unsigned int i = 1; i < PointNum - 1; i++) {
            Energy += ( pow((Grid[0][i][0] - TempGrid[0][i][0] / Test -> GetTimeStep()), 2.0) +
                        pow(TempGrid[0][i][1] / TempGrid[0][i][0], 2.0) +
                        pow(TempGrid[0][i][2] / TempGrid[0][i][0], 2.0));
        }

        std::cout << Energy << " : " << iter << std::endl;
        EnergyFile << Energy << std::endl;

        Grid = TempGrid;
    }

    EnergyFile.close();
}

void Solver::InitGrid(double ExcitationP, double VXP, double VYP, int PointsNumP) {
    for (int i = 0; i < PointsNumP; i++) {
        Grid.emplace_back(std :: vector <dVectorND <double>>());

        for (int j = 0; j < PointsNumP; j++) {
            Grid.back().emplace_back(dVectorND <double>({ExcitationP, VXP * ExcitationP, VYP * ExcitationP, -1, 1, 1}));
        }
    }
}