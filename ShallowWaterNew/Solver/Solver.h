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
//-----------------------------//
class dRichtmyerSolver2D : public dRichtmyer2D <dVectorND <double>> {
public:
    dRichtmyerSolver2D(double TimeStepP, double xStepP, double yStepP) {
        TimeStep = TimeStepP;
        xStep = xStepP;
        yStep = yStepP;

        Homogeneous = false;
    }
    ~dRichtmyerSolver2D() = default;

    void SetInitialData(unsigned long xSizeP, unsigned long ySizeP) {
        xSize = xSizeP;
        ySize = ySizeP;

        DataFirst.resize(xSizeP);
        DataSecond.resize(xSizeP);

        for (unsigned long i = 0; i < xSizeP; i++) {
            for (unsigned long j = 0; j < ySizeP; j++) {
                DataFirst[i].emplace_back(dVectorND <double> ({10.0, 0.0, 0.0, 0.0, 0.0}));
                DataSecond[i].emplace_back(dVectorND <double> ({10.0, 0.0, 0.0, 0.0, 0.0}));
            }
        }

        SetExcitation();

        for (unsigned long i = 0; i < ySizeP; i++) {
            Gradient.emplace_back(0.1 / ySizeP * i);
        }

        AmplitudeFile.open(SavePath + "Amplitude.dat");
        xVelocityFile.open(SavePath + "xVelocity.dat");
        yVelocityFile.open(SavePath + "yVelocity.dat");
        xFieldFile.open(SavePath + "xField.dat");
        yFieldFile.open(SavePath + "yField.dat");
    }
    void SetSavePath(const std::string& PathP) {
        SavePath = PathP;
    }
    void AppendData() {
        for (int j = 0; j < ySize; j++) {
            for (int i = 0; i < xSize; i++) {
                AmplitudeFile << (*CurrentData)[i][j][0] << "\t";
                xVelocityFile << (*CurrentData)[i][j][1] << "\t";
                yVelocityFile << (*CurrentData)[i][j][2] << "\t";
                xFieldFile << (*CurrentData)[i][j][3] << "\t";
                yFieldFile << (*CurrentData)[i][j][4] << "\t";
            }

            AmplitudeFile << std::endl;
            xVelocityFile << std::endl;
            yVelocityFile << std::endl;
            xFieldFile << std::endl;
            yFieldFile << std::endl;
        }
    }

    double GetFullEnergy() {
        double EnergyL = 0.0;

        for (const auto& LineI : (*CurrentData)) {
            for (const auto& ValueI : LineI) {
                EnergyL += (g * (ValueI[0] - 10.0) +
                            pow(ValueI[1] / ValueI[0], 2.0) +
                            pow(ValueI[2] / ValueI[0], 2.0) +
                            pow(ValueI[3] / ValueI[0], 2.0) +
                            pow(ValueI[4] / ValueI[0], 2.0));
            }
        }

        return EnergyL;
    }
    double GetMaxAmplitude() {
        double AbsMaxL = 0.0;

        for (unsigned long i = 0; i < xSize; i++) {
            for (unsigned long j = 0; j < ySize; j++) {
                double Value =  fabs((*CurrentData)[i][j][0]) +
                                fabs((*CurrentData)[i][j][1]) +
                                fabs((*CurrentData)[i][j][2]) +
                                fabs((*CurrentData)[i][j][3]) +
                                fabs((*CurrentData)[i][j][4]);

                if (AbsMaxL < Value) {
                    AbsMaxL = Value;
                }
            }
        }

        return AbsMaxL;
    }

    void SaveData() {
        AmplitudeFile.close();
        xVelocityFile.close();
        yVelocityFile.close();
        xFieldFile.close();
        yFieldFile.close();
    }
private:
    const double g = 9.81;
    const double B_0 = 0.5;
    const double f_0 = 0.1;

    std::vector <double> Gradient;

    //----------//

    dVectorND <double> xFunc(const dVectorND <double>& U) override {
        return dVectorND <double> ({U[1],
                                    (pow(U[1], 2.0) - pow(U[3], 2.0)) / U[0] + 0.5 * g * pow(U[0], 2.0),
                                    (U[1] * U[2] - U[3] * U[4]) / U[0],
                                    0,
                                    (U[1] * U[4] - U[2] * U[3]) / U[0]});
    }
    dVectorND <double> yFunc(const dVectorND <double>& U) override {
        return dVectorND <double> ({U[2],
                                    (U[1] * U[2] - U[3] * U[4]) / U[0],
                                    (pow(U[2], 2.0) - pow(U[4], 2.0)) / U[0] + 0.5 * g * pow(U[0], 2.0),
                                    (U[2] * U[3] - U[1] * U[4]) / U[0],
                                    0});
    }
    dVectorND <double> AbsValFunc(int xPosP, int yPosP) override {
        return dVectorND <double> ({0,
                                    B_0 * (*CurrentData)[xPosP][yPosP][3] / (*CurrentData)[xPosP][yPosP][0] - (*CurrentData)[xPosP][yPosP][2] * (f_0 + Gradient[yPosP]),
                                    B_0 * (*CurrentData)[xPosP][yPosP][4] / (*CurrentData)[xPosP][yPosP][0] - (*CurrentData)[xPosP][yPosP][1] * (f_0 + Gradient[yPosP]),
                                    -B_0 * (*CurrentData)[xPosP][yPosP][1] / (*CurrentData)[xPosP][yPosP][0],
                                    -B_0 * (*CurrentData)[xPosP][yPosP][2] / (*CurrentData)[xPosP][yPosP][0]});
    }
    dVectorND <double> Viscosity(int xPosP, int yPosP) override {
        long xIndex_plus_1;
        long xIndex_minus_1;
        long yIndex_plus_1;
        long yIndex_minus_1;

        xIndex_plus_1 = (xPosP + 1 == xSize ? 0 : xPosP + 1);
        xIndex_minus_1 = (xPosP - 1 < 0 ? xSize - 1 : xPosP - 1);
        yIndex_plus_1 = (yPosP + 1 == ySize ? 0 : yPosP + 1);
        yIndex_minus_1 = (yPosP - 1 < 0 ? ySize - 1 : yPosP - 1);

        double v_x_xx = ((*CurrentData)[xIndex_minus_1][yPosP][1] / (*CurrentData)[xIndex_minus_1][yPosP][0] +
                         (*CurrentData)[xPosP][yPosP][1] / (*CurrentData)[xPosP][yPosP][0] * 2 +
                         (*CurrentData)[xIndex_plus_1][yPosP][1] / (*CurrentData)[xIndex_plus_1][yPosP][0]) /
                        (xStep * xStep);
        double v_x_yy = ((*CurrentData)[xPosP][yIndex_minus_1][1] / (*CurrentData)[xPosP][yIndex_minus_1][0] +
                         (*CurrentData)[xPosP][yPosP][1] / (*CurrentData)[xPosP][yPosP][0] * 2 +
                         (*CurrentData)[xPosP][yIndex_plus_1][1] / (*CurrentData)[xPosP][yIndex_plus_1][0]) /
                        (yStep * yStep);
        double v_y_xx = ((*CurrentData)[xIndex_minus_1][yPosP][2] / (*CurrentData)[xIndex_minus_1][yPosP][0] +
                         (*CurrentData)[xPosP][yPosP][2] / (*CurrentData)[xPosP][yPosP][0] * 2 +
                         (*CurrentData)[xIndex_plus_1][yPosP][2] / (*CurrentData)[xIndex_plus_1][yPosP][0]) /
                        (xStep * xStep);
        double v_y_yy = ((*CurrentData)[xPosP][yIndex_minus_1][2] / (*CurrentData)[xPosP][yIndex_minus_1][0] +
                         (*CurrentData)[xPosP][yPosP][2] / (*CurrentData)[xPosP][yPosP][0] * 2 +
                         (*CurrentData)[xPosP][yIndex_plus_1][2] / (*CurrentData)[xPosP][yIndex_plus_1][0]) /
                        (yStep * yStep);
        double B_x_xx = ((*CurrentData)[xIndex_minus_1][yPosP][3] / (*CurrentData)[xIndex_minus_1][yPosP][0] +
                         (*CurrentData)[xPosP][yPosP][3] / (*CurrentData)[xPosP][yPosP][0] * 2 +
                         (*CurrentData)[xIndex_plus_1][yPosP][3] / (*CurrentData)[xIndex_plus_1][yPosP][0]) /
                        (xStep * xStep);
        double B_x_yy = ((*CurrentData)[xPosP][yIndex_minus_1][3] / (*CurrentData)[xPosP][yIndex_minus_1][0] +
                         (*CurrentData)[xPosP][yPosP][3] / (*CurrentData)[xPosP][yPosP][0] * 2 +
                         (*CurrentData)[xPosP][yIndex_plus_1][3] / (*CurrentData)[xPosP][yIndex_plus_1][0]) /
                        (yStep * yStep);
        double B_y_xx = ((*CurrentData)[xIndex_minus_1][yPosP][4] / (*CurrentData)[xIndex_minus_1][yPosP][0] +
                         (*CurrentData)[xPosP][yPosP][4] / (*CurrentData)[xPosP][yPosP][0] * 2 +
                         (*CurrentData)[xIndex_plus_1][yPosP][4] / (*CurrentData)[xIndex_plus_1][yPosP][0]) /
                        (xStep * xStep);
        double B_y_yy = ((*CurrentData)[xPosP][yIndex_minus_1][4] / (*CurrentData)[xPosP][yIndex_minus_1][0] +
                         (*CurrentData)[xPosP][yPosP][4] / (*CurrentData)[xPosP][yPosP][0] * 2 +
                         (*CurrentData)[xPosP][yIndex_plus_1][4] / (*CurrentData)[xPosP][yIndex_plus_1][0]) /
                        (yStep * yStep);

        return dVectorND <double> ({0,
                                    (*CurrentData)[xPosP][yPosP][0] * (v_x_xx + v_x_yy),
                                    (*CurrentData)[xPosP][yPosP][0] * (v_y_xx + v_y_yy),
                                    (*CurrentData)[xPosP][yPosP][0] * (B_x_xx + B_x_yy),
                                    (*CurrentData)[xPosP][yPosP][0] * (B_y_xx + B_y_yy)});
    }

    //----------//

    std::ofstream AmplitudeFile;
    std::ofstream xVelocityFile;
    std::ofstream yVelocityFile;
    std::ofstream xFieldFile;
    std::ofstream yFieldFile;

    std::string SavePath = "./";

    //----------//

    void SetPlotParameters(long xMinP, long xMaxP, long yMinP, long yMaxP) {
        AmplitudeFile   << xSize        << "\t"     << ySize    << "\t"
                        << xMinP        << "\t"     << xMaxP    << "\t"
                        << yMinP        << "\t"     << yMaxP    << "\t"
                        << "\n";
        xVelocityFile   << xSize        << "\t"     << ySize    << "\t"
                        << xMinP        << "\t"     << xMaxP    << "\t"
                        << yMinP        << "\t"     << yMaxP    << "\t"
                        << "\n";
        yVelocityFile   << xSize        << "\t"     << ySize    << "\t"
                        << xMinP        << "\t"     << xMaxP    << "\t"
                        << yMinP        << "\t"     << yMaxP    << "\t"
                        << "\n";
        xFieldFile      << xSize        << "\t"     << ySize    << "\t"
                        << xMinP        << "\t"     << xMaxP    << "\t"
                        << yMinP        << "\t"     << yMaxP    << "\t"
                        << "\n";
        yFieldFile      << xSize        << "\t"     << ySize    << "\t"
                        << xMinP        << "\t"     << xMaxP    << "\t"
                        << yMinP        << "\t"     << yMaxP    << "\t"
                        << "\n";
    }
    void SetExcitation() {
        double Fraction = 1.0 / (2.0 * M_PI * 20.0 * xStep * 20.0 * yStep);

        for (int i = -100; i < 100; i++) {
            for (int j = -100; j < 100; j++) {
                (*CurrentData)[xSize / 2 + i][ySize / 2 + j][0] = 10.0 + Fraction * exp(-0.5 * (pow(i / 20.0, 2.0) + pow(j / 20.0, 2.0)));
            }
        }
//        for (int i = 0; i < 200; i++) {
//            for (int j = 0; j < 200; j++) {
//                (*CurrentData)[i][j][0] = 12 - 2.0 / 200 * i;
//            }
//        }
    }
};
//-----------------------------//
#endif
