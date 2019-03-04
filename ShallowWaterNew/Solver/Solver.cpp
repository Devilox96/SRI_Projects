#include "Solver.h"
//-----------------------------//
Solver::Solver() {
    Test = new dRichtmyerSolver(0.001, 0.1, 0.1);

    InitGrid(10.0, 0.0, 0.0, 200);

    Grid[100][100][0] = 12.0;
}
Solver::~Solver() {
    delete Test;
}

void Solver :: Calc(unsigned int StepsP) {
    auto Start = std::chrono::system_clock::now();
    //----------//
    std::ofstream EnergyFile;
    EnergyFile.open("Energy.dat", std::ios::app);


    unsigned int PointNum = 200;

//    InitGrid(10.0, 0.0, 0.0, PointNum);
//
//    Grid[100][100][0] = 12.0;

    std :: vector <std :: vector <dVectorND <double>>> TempGrid = Grid;

    for (unsigned int iter = 0; iter < StepsP; iter++) {
        double Energy = 0.0;

        for (unsigned int i = 1; i < PointNum - 1; i++) {
            for (unsigned int j = 1; j < PointNum - 1; j++) {
                dVectorND <double> TempH = Test -> Solve2D(Grid[i][j], Grid[i - 1][j], Grid[i + 1][j], Grid[i][j - 1], Grid[i][j + 1]);

                TempGrid[i][j] = TempH;

                Energy += ( pow((Grid[i][j][0] - TempGrid[i][j][0]) / Test -> GetTimeStep(), 2.0) +
                            pow(TempGrid[i][j][1] / TempGrid[i][j][0], 2.0) +
                            pow(TempGrid[i][j][2] / TempGrid[i][j][0], 2.0));
            }
        }

        for (unsigned int i = 0; i < PointNum; i++) {
            TempGrid[i][0] = TempGrid[i][1];
            TempGrid[i][PointNum - 1] = TempGrid[i][PointNum - 2];
            TempGrid[0][i] = TempGrid[1][i];
            TempGrid[PointNum - 1][i] = TempGrid[PointNum - 2][i];
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

    //----------//
    auto Stop = std::chrono::system_clock::now();

    std::chrono::duration <double> Time = Stop - Start;
    std::cout << Time.count() << std::endl;
}

void Solver::InitGrid(double ExcitationP, double VXP, double VYP, int PointsNumP) {
    for (int i = 0; i < PointsNumP; i++) {
        Grid.emplace_back(std :: vector <dVectorND <double>>());

        for (int j = 0; j < PointsNumP; j++) {
            Grid.back().emplace_back(dVectorND <double>({ExcitationP, VXP * ExcitationP, VYP * ExcitationP, -1, 1, 1}));
        }
    }
}