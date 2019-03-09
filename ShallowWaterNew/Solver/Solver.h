#ifndef SOLVER_H
#define SOLVER_H
//-----------------------------//
#include <vector>
#include <chrono>
#include <fstream>
#include <cmath>
//-----------------------------//
#include "../Libs/dMath/Core/dVectorND.h"
#include "../Libs/dMath/NumerCalc/dRichtmyer2D.h"
#include "OutputDataFormat.h"
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

    void SetInitialData(unsigned long xSizeP, unsigned long ySizeP) {
        xSize = xSizeP;
        ySize = ySizeP;

        DataFirst.resize(xSizeP);
        DataSecond.resize(xSizeP);

        for (unsigned long i = 0; i < xSizeP; i++) {
            for (unsigned long j = 0; j < ySizeP; j++) {
                DataFirst[i].emplace_back(dVectorND <double> ({10.0, 0.0, 0.0, 0.0, 0.0, 0.0}));
                DataSecond[i].emplace_back(dVectorND <double> ({10.0, 0.0, 0.0, 0.0, 0.0, 0.0}));
            }
        }

        DataFirst[xSizeP / 2][ySizeP / 2][0] = 12.0;
    }
    double GetFullEnergy() {
        double EnergyL = 0.0;

        for (const auto& LineI : (*CurrentData)) {
            for (const auto& ValueI : LineI) {
                EnergyL += (g * ValueI[0] + pow(ValueI[1] / ValueI[0], 2.0) + pow(ValueI[2] / ValueI[0], 2.0));
            }
        }

        return EnergyL;
    }
    void SaveData(const std::string& PathP) {
        std::ofstream OutputL;
        OutputL.open(PathP);

        for (int i = 0; i < xSize; i++) {
            for (int j = 0; j < ySize; j++) {
                OutputL << (*CurrentData)[i][j][0] << "\t";
            }

            OutputL << std::endl;
        }

        OutputL.close();
    }
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

    //----------//

    OutputDataFormat <double>* AmplitudeSaver;
    OutputDataFormat <double>* xVelocitySaver;
    OutputDataFormat <double>* yVelocitySaver;
    OutputDataFormat <double>* xFieldSaver;
    OutputDataFormat <double>* yFieldSaver;
};
//-----------------------------//
#endif
