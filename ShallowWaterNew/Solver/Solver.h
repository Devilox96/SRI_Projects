#ifndef SOLVER_H
#define SOLVER_H
//-----------------------------//
#include <vector>
#include <chrono>
#include <fstream>
//-----------------------------//
#include "../Libs/dMath/Core/dVectorND.h"
#include "../Libs/dMath/NumerCalc/dRichtmyer2D.h"
//-----------------------------//
class dRichtmyerSolver2D : public dRichtmyer2D <dVectorND <double>> {
public:
    dRichtmyerSolver2D() = default;
    dRichtmyerSolver2D(double TimeStepP, double xStepP, double yStepP) {
        TimeStep = TimeStepP;
        xStep = xStepP;
        yStep = yStepP;

        Homogeneous = true;
    }
    ~dRichtmyerSolver2D() = default;
private:
    const double g = 9.81;
    const double B_0 = 10.0;
    const double f_0 = 10.0;

    //----------//

    dVectorND <double> xFunc(const dVectorND <double>& U) override {
        return dVectorND <double> ({U[1],
                                    (pow(U[1], 2.0) - pow(U[3], 2.0)) / U[0] + 0.5 * g * pow(U[0], 2.0),
                                    (U[1] * U[2] - U[3] * U[4]) / U[0],
                                    0,
                                    (U[1] * U[4] - U[2] * U[3]) / U[0],
                                    B_0 * U[1] / U[0]});
    }
    dVectorND <double> yFunc(const dVectorND <double>& U) override {
        return dVectorND <double> ({U[2],
                                    (U[1] * U[2] - U[3] * U[4]) / U[0],
                                    (pow(U[2], 2.0) - pow(U[4], 2.0)) / U[0] + 0.5 * g * pow(U[0], 2.0),
                                    (U[2] * U[3] - U[1] * U[4]) / U[0],
                                    0,
                                    B_0 * U[2] / U[0]});
    }
    dVectorND <double> AbsValFunc(const dVectorND <double>& U) override {
        return dVectorND <double> ({0,
                                    B_0 * U[3] / U[0] - U[2] * f_0,
                                    B_0 * U[4] / U[0] - U[1] * f_0,
                                    -B_0 * U[1] / U[0],
                                    -B_0 * U[2] / U[0],
                                    0});
    }
};
//-----------------------------//
class Solver {
public:
    Solver();
    ~Solver();

    void Calc(unsigned int StepsP);

    std :: vector <std :: vector <dVectorND <double>>> Grid;
    dRichtmyerSolver2D* Test;

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
