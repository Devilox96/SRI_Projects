#ifndef SOLVER_H
#define SOLVER_H
//-----------------------------//
#include <vector>
#include <chrono>
#include <fstream>
//-----------------------------//
#include "../Libs/dMath/Core/dVectorND.h"
#include "../Libs/dMath/NumerCalc/dRichtmyer.h"
//-----------------------------//
class Solver {
public:
    Solver();
    ~Solver();

    void Calc(unsigned int StepsP);

    std :: vector <std :: vector <dVectorND <double>>> Grid;
    dRichtmyerSolver* Test;

    void SaveData(const std::string& PathP) {
        std::ofstream OutputL;
        OutputL.open(PathP);

        for (int i = 0; i < 200; i++) {
            for (int j = 0; j < 200; j++) {
                OutputL << Grid[i][j][0] << "\t";
            }

            OutputL << std::endl;
        }

        OutputL.close();
    }
private:
    void InitGrid(double ExcitationP, double VXP, double VYP, int PointsNumP);
};
//-----------------------------//
#endif
