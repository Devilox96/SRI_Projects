#ifndef SOLVER_H
#define SOLVER_H
//-----------------------------//
#include <vector>
#include <chrono>
#include <fstream>
#include <cmath>
//-----------------------------//
#include "../Libs/dMath/Core/dVector.h"
#include "../Libs/dMath/NumerCalc/dRichtmyer2D.h"
#include "../Libs/dMath/NumerCalc/dLaxFriedrichs2D.h"
//-----------------------------//
class dRichtmyerSolver2D : public dRichtmyer2D <dVector <double, 5>> {
public:
    dRichtmyerSolver2D(double TimeStepP, double xStepP, double yStepP) {
        mStepTime = TimeStepP;
        mStepX = xStepP;
        mStepY = yStepP;

//        Homogeneous = false;
    }
    ~dRichtmyerSolver2D() = default;

    void SetInitialData(unsigned long xSizeP, unsigned long ySizeP) {
        xSize = xSizeP;
        ySize = ySizeP;

        DataFirst.resize(xSizeP);
        DataSecond.resize(xSizeP);

        for (unsigned long i = 0; i < xSizeP; i++) {
            for (unsigned long j = 0; j < ySizeP; j++) {
                DataFirst[i].emplace_back(dVector <double, 5> (10.0, 0.0, 0.0, 0.0, 0.0));
                DataSecond[i].emplace_back(dVector <double, 5> (10.0, 0.0, 0.0, 0.0, 0.0));
            }
        }

        SetExcitation();
//        setSlope();

        for (unsigned long i = 0; i < ySizeP; i++) {
//            Gradient.emplace_back(0.1 / ySizeP * i);
            Gradient.emplace_back(14.58e-05 * i);
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
    double getMaxSpeed() {
        double VelMax = 0.0;

        for (unsigned long i = 0; i < xSize; i++) {
            for (unsigned long j = 0; j < ySize; j++) {
                double CurVal = 0.0;

                if ((*CurrentData)[i][j][0] != 0) {
                    CurVal = sqrt((*CurrentData)[i][j][1] * (*CurrentData)[i][j][1] +
                                  (*CurrentData)[i][j][2] * (*CurrentData)[i][j][2]) / (*CurrentData)[i][j][0];
                }

                if (VelMax < CurVal) {
                    VelMax = CurVal;
                }
            }
        }

        return VelMax;
    }

    void SaveData() {
        AmplitudeFile.close();
        xVelocityFile.close();
        yVelocityFile.close();
        xFieldFile.close();
        yFieldFile.close();
    }

    double getStepTime() {
        return mStepTime;
    }

    double getStepX() {
        return mStepX;
    }
    double getStepY() {
        return mStepY;
    }

    //----------//

    void solveGridRichtmyer() {
        long xIndex_plus_1;
        long xIndex_plus_2;
        long xIndex_minus_1;
        long xIndex_minus_2;
        long yIndex_plus_1;
        long yIndex_plus_2;
        long yIndex_minus_1;
        long yIndex_minus_2;

        for (int i = 0; i < xSize; i++) {
            for (int j = 0; j < ySize; j++) {
                xIndex_plus_1 = (i + 1 == xSize ? 0 : i + 1);
                xIndex_plus_2 = (i + 2 >= xSize ? i + 2 - xSize : i + 2);
                xIndex_minus_1 = (i - 1 < 0 ? xSize - 1 : i - 1);
                xIndex_minus_2 = (i - 2 < 0 ? xSize + i - 2 : i - 2);
                yIndex_plus_1 = (j + 1 == ySize ? 0 : j + 1);
                yIndex_plus_2 = (j + 2 >= ySize ? j + 2 - ySize : j + 2);
                yIndex_minus_1 = (j - 1 < 0 ? ySize - 1 : j - 1);
                yIndex_minus_2 = (j - 2 < 0 ? ySize + j - 2 : j - 2);

                (*TempData)[i][j] = solve((*CurrentData)[i][j],
                                          (*CurrentData)[xIndex_minus_2][j],
                                          (*CurrentData)[xIndex_plus_2][j],
                                          (*CurrentData)[i][yIndex_minus_2],
                                          (*CurrentData)[i][yIndex_plus_2],
                                          (*CurrentData)[xIndex_plus_1][yIndex_plus_1],
                                          (*CurrentData)[xIndex_minus_1][yIndex_minus_1],
                                          (*CurrentData)[xIndex_plus_1][yIndex_minus_1],
                                          (*CurrentData)[xIndex_minus_1][yIndex_plus_1]) -
                        viscosity(i, j) * 0.00007;
            }
        }

        std::swap(CurrentData, TempData);
    }
    void solveGridRichtmyerZwas() {
        long xIndex_plus_1;
        long xIndex_minus_1;
        long yIndex_plus_1;
        long yIndex_minus_1;

        dVector <double, 5> Extra(0, 0, 0, 0, 0);

        for (int i = 0; i < xSize; i++) {
            for (int j = 0; j < ySize; j++) {
                xIndex_plus_1 = (i + 1 == xSize ? 0 : i + 1);
                xIndex_minus_1 = (i - 1 < 0 ? xSize - 1 : i - 1);
                yIndex_plus_1 = (j + 1 == ySize ? 0 : j + 1);
                yIndex_minus_1 = (j - 1 < 0 ? ySize - 1 : j - 1);

                Extra = nonhomogenPart(i, j) - viscosity(i, j) * 0.000025;

                (*TempData)[i][j] = solveZwas((*CurrentData)[i][j],
                                              (*CurrentData)[xIndex_minus_1][j],
                                              (*CurrentData)[xIndex_plus_1][j],
                                              (*CurrentData)[i][yIndex_minus_1],
                                              (*CurrentData)[i][yIndex_plus_1],
                                              (*CurrentData)[xIndex_plus_1][yIndex_plus_1],
                                              (*CurrentData)[xIndex_minus_1][yIndex_minus_1],
                                              (*CurrentData)[xIndex_plus_1][yIndex_minus_1],
                                              (*CurrentData)[xIndex_minus_1][yIndex_plus_1],
                                              Extra);
            }
        }

        std::swap(CurrentData, TempData);
    }
private:
    const double g = 9.81e-03;
    const double B_0 = 65.0e-02;
    const double f_0 = 14.58e-05;

    unsigned long xSize = 1;
    unsigned long ySize = 1;

    std::vector <double> Gradient;

    std::vector <std::vector <dVector <double, 5>>> DataFirst;
    std::vector <std::vector <dVector <double, 5>>> DataSecond;

    std::vector <std::vector <dVector <double, 5>>>* CurrentData = &DataFirst;
    std::vector <std::vector <dVector <double, 5>>>* TempData = &DataSecond;

    //----------//

    dVector <double, 5> funcX(const dVector <double, 5>& U) override {
        return dVector <double, 5> (U[1],
                                    (pow(U[1], 2.0) - pow(U[3], 2.0)) / U[0] + 0.5 * g * pow(U[0], 2.0),
                                    (U[1] * U[2] - U[3] * U[4]) / U[0],
                                    0,
                                    (U[1] * U[4] - U[2] * U[3]) / U[0]);
    }
    dVector <double, 5> funcY(const dVector <double, 5>& U) override {
        return dVector <double, 5> (U[2],
                                    (U[1] * U[2] - U[3] * U[4]) / U[0],
                                    (pow(U[2], 2.0) - pow(U[4], 2.0)) / U[0] + 0.5 * g * pow(U[0], 2.0),
                                    (U[2] * U[3] - U[1] * U[4]) / U[0],
                                    0);
    }
    dVector <double, 5> nonhomogenPart(int xPosP, int yPosP) {
        return dVector <double, 5> (
                0,
                B_0 * (*CurrentData)[xPosP][yPosP][3] / (*CurrentData)[xPosP][yPosP][0] - (*CurrentData)[xPosP][yPosP][2] * (f_0 + Gradient[yPosP]),
                B_0 * (*CurrentData)[xPosP][yPosP][4] / (*CurrentData)[xPosP][yPosP][0] - (*CurrentData)[xPosP][yPosP][1] * (f_0 + Gradient[yPosP]),
                -B_0 * (*CurrentData)[xPosP][yPosP][1] / (*CurrentData)[xPosP][yPosP][0],
                -B_0 * (*CurrentData)[xPosP][yPosP][2] / (*CurrentData)[xPosP][yPosP][0]);
    }
    dVector <double, 5> viscosity(int xPosP, int yPosP) {
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
                        (mStepX * mStepX);
        double v_x_yy = ((*CurrentData)[xPosP][yIndex_minus_1][1] / (*CurrentData)[xPosP][yIndex_minus_1][0] +
                         (*CurrentData)[xPosP][yPosP][1] / (*CurrentData)[xPosP][yPosP][0] * 2 +
                         (*CurrentData)[xPosP][yIndex_plus_1][1] / (*CurrentData)[xPosP][yIndex_plus_1][0]) /
                        (mStepY * mStepY);
        double v_y_xx = ((*CurrentData)[xIndex_minus_1][yPosP][2] / (*CurrentData)[xIndex_minus_1][yPosP][0] +
                         (*CurrentData)[xPosP][yPosP][2] / (*CurrentData)[xPosP][yPosP][0] * 2 +
                         (*CurrentData)[xIndex_plus_1][yPosP][2] / (*CurrentData)[xIndex_plus_1][yPosP][0]) /
                        (mStepX * mStepX);
        double v_y_yy = ((*CurrentData)[xPosP][yIndex_minus_1][2] / (*CurrentData)[xPosP][yIndex_minus_1][0] +
                         (*CurrentData)[xPosP][yPosP][2] / (*CurrentData)[xPosP][yPosP][0] * 2 +
                         (*CurrentData)[xPosP][yIndex_plus_1][2] / (*CurrentData)[xPosP][yIndex_plus_1][0]) /
                        (mStepY * mStepY);
//        double B_x_xx = ((*CurrentData)[xIndex_minus_1][yPosP][3] / (*CurrentData)[xIndex_minus_1][yPosP][0] +
//                         (*CurrentData)[xPosP][yPosP][3] / (*CurrentData)[xPosP][yPosP][0] * 2 +
//                         (*CurrentData)[xIndex_plus_1][yPosP][3] / (*CurrentData)[xIndex_plus_1][yPosP][0]) /
//                        (mStepX * mStepX);
//        double B_x_yy = ((*CurrentData)[xPosP][yIndex_minus_1][3] / (*CurrentData)[xPosP][yIndex_minus_1][0] +
//                         (*CurrentData)[xPosP][yPosP][3] / (*CurrentData)[xPosP][yPosP][0] * 2 +
//                         (*CurrentData)[xPosP][yIndex_plus_1][3] / (*CurrentData)[xPosP][yIndex_plus_1][0]) /
//                        (mStepY * mStepY);
//        double B_y_xx = ((*CurrentData)[xIndex_minus_1][yPosP][4] / (*CurrentData)[xIndex_minus_1][yPosP][0] +
//                         (*CurrentData)[xPosP][yPosP][4] / (*CurrentData)[xPosP][yPosP][0] * 2 +
//                         (*CurrentData)[xIndex_plus_1][yPosP][4] / (*CurrentData)[xIndex_plus_1][yPosP][0]) /
//                        (mStepX * mStepX);
//        double B_y_yy = ((*CurrentData)[xPosP][yIndex_minus_1][4] / (*CurrentData)[xPosP][yIndex_minus_1][0] +
//                         (*CurrentData)[xPosP][yPosP][4] / (*CurrentData)[xPosP][yPosP][0] * 2 +
//                         (*CurrentData)[xPosP][yIndex_plus_1][4] / (*CurrentData)[xPosP][yIndex_plus_1][0]) /
//                        (mStepY * mStepY);

//        return dVector <double, 5> (0,
//                                    (*CurrentData)[xPosP][yPosP][0] * (v_x_xx + v_x_yy),
//                                    (*CurrentData)[xPosP][yPosP][0] * (v_y_xx + v_y_yy),
//                                    (*CurrentData)[xPosP][yPosP][0] * (B_x_xx + B_x_yy),
//                                    (*CurrentData)[xPosP][yPosP][0] * (B_y_xx + B_y_yy));
        return dVector <double, 5> (0,
                                    (*CurrentData)[xPosP][yPosP][0] * (v_x_xx + v_x_yy),
                                    (*CurrentData)[xPosP][yPosP][0] * (v_y_xx + v_y_yy),
                                    0,
                                    0);

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
        double Fraction = 1.0 / (2.0 * M_PI * 20.0 * mStepX * 20.0 * mStepY);

        for (int i = -100; i < 100; i++) {
            for (int j = -100; j < 100; j++) {
                (*CurrentData)[xSize / 2 + i][ySize / 2 + j][0] = 10.0 + Fraction * exp(-0.5 * (pow(i / 20.0, 2.0) + pow(j / 20.0, 2.0))) * 10.0;
            }
        }
//        for (int i = 0; i < 200; i++) {
//            for (int j = 0; j < 200; j++) {
//                (*CurrentData)[i][j][0] = 12 - 2.0 / 200 * i;
//            }
//        }
    }
    void setSlope() {
        for (int i = 0; i < ySize; i++) {
            for (int j = 0; j < xSize; j++) {
                (*CurrentData)[j][i][0] = 10.0 + double(i) / ySize * 0.01;
//                (*CurrentData)[j][i][1] = (*CurrentData)[j][i][0] * 0.01 * double(i) / ySize;
            }
        }
    }
};
//-----------------------------//
class dLaxFriedrichsSolver2D : public dLaxFriedrichs2D <dVector <double, 5>> {
public:
    dLaxFriedrichsSolver2D(double TimeStepP, double xStepP, double yStepP) {
        mStepTime = TimeStepP;
        mStepX = xStepP;
        mStepY = yStepP;

//        Homogeneous = false;
    }
    ~dLaxFriedrichsSolver2D() = default;

    void SetInitialData(unsigned long xSizeP, unsigned long ySizeP) {
        xSize = xSizeP;
        ySize = ySizeP;

        DataFirst.resize(xSizeP);
        DataSecond.resize(xSizeP);

        for (unsigned long i = 0; i < xSizeP; i++) {
            for (unsigned long j = 0; j < ySizeP; j++) {
                DataFirst[i].emplace_back(dVector <double, 5> (10.0, 0.0, 0.0, 0.0, 0.0));
                DataSecond[i].emplace_back(dVector <double, 5> (10.0, 0.0, 0.0, 0.0, 0.0));
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

    double getMaxSpeed() {
        double VelMax = 0.0;

        for (unsigned long i = 0; i < xSize; i++) {
            for (unsigned long j = 0; j < ySize; j++) {
                double CurVal = 0.0;

                if ((*CurrentData)[i][j][0] != 0) {
                    CurVal = sqrt((*CurrentData)[i][j][1] * (*CurrentData)[i][j][1] +
                                         (*CurrentData)[i][j][2] * (*CurrentData)[i][j][2]) / (*CurrentData)[i][j][0];
                }

                if (VelMax < CurVal) {
                    VelMax = CurVal;
                }
            }
        }

        return VelMax;
    }

    void SaveData() {
        AmplitudeFile.close();
        xVelocityFile.close();
        yVelocityFile.close();
        xFieldFile.close();
        yFieldFile.close();
    }

    double getStepTime() {
        return mStepTime;
    }

    double getStepX() {
        return mStepX;
    }
    double getStepY() {
        return mStepY;
    }

    //----------//

    void solveGrid() {
        long xIndex_plus_1;
        long xIndex_minus_1;
        long yIndex_plus_1;
        long yIndex_minus_1;

        for (int i = 0; i < xSize; i++) {
            for (int j = 0; j < ySize; j++) {
                xIndex_plus_1 = (i + 1 == xSize ? 0 : i + 1);
                xIndex_minus_1 = (i - 1 < 0 ? xSize - 1 : i - 1);
                yIndex_plus_1 = (j + 1 == ySize ? 0 : j + 1);
                yIndex_minus_1 = (j - 1 < 0 ? ySize - 1 : j - 1);

                (*TempData)[i][j] = solve((*CurrentData)[xIndex_minus_1][j],
                                              (*CurrentData)[xIndex_plus_1][j],
                                              (*CurrentData)[i][yIndex_minus_1],
                                              (*CurrentData)[i][yIndex_plus_1]);
            }
        }

        std::swap(CurrentData, TempData);
    }
private:
    const double g = 9.81e-03;
    const double B_0 = 0.5;
    const double f_0 = 0.1;

    unsigned long xSize = 1;
    unsigned long ySize = 1;

    std::vector <double> Gradient;

    std::vector <std::vector <dVector <double, 5>>> DataFirst;
    std::vector <std::vector <dVector <double, 5>>> DataSecond;

    std::vector <std::vector <dVector <double, 5>>>* CurrentData = &DataFirst;
    std::vector <std::vector <dVector <double, 5>>>* TempData = &DataSecond;

    //----------//

    dVector <double, 5> funcX(const dVector <double, 5>& U) override {
        return dVector <double, 5> (U[1],
                                    (pow(U[1], 2.0) - pow(U[3], 2.0)) / U[0] + 0.5 * g * pow(U[0], 2.0),
                                    (U[1] * U[2] - U[3] * U[4]) / U[0],
                                    0,
                                    (U[1] * U[4] - U[2] * U[3]) / U[0]);
    }
    dVector <double, 5> funcY(const dVector <double, 5>& U) override {
        return dVector <double, 5> (U[2],
                                    (U[1] * U[2] - U[3] * U[4]) / U[0],
                                    (pow(U[2], 2.0) - pow(U[4], 2.0)) / U[0] + 0.5 * g * pow(U[0], 2.0),
                                    (U[2] * U[3] - U[1] * U[4]) / U[0],
                                    0);
    }
//    dVectorND <double> AbsValFunc(int xPosP, int yPosP) override {
//        return dVectorND <double> ({0,
//                                    B_0 * (*CurrentData)[xPosP][yPosP][3] / (*CurrentData)[xPosP][yPosP][0] - (*CurrentData)[xPosP][yPosP][2] * (f_0 + Gradient[yPosP]),
//                                    B_0 * (*CurrentData)[xPosP][yPosP][4] / (*CurrentData)[xPosP][yPosP][0] - (*CurrentData)[xPosP][yPosP][1] * (f_0 + Gradient[yPosP]),
//                                    -B_0 * (*CurrentData)[xPosP][yPosP][1] / (*CurrentData)[xPosP][yPosP][0],
//                                    -B_0 * (*CurrentData)[xPosP][yPosP][2] / (*CurrentData)[xPosP][yPosP][0]});
//    }
//    dVectorND <double> viscosity(int xPosP, int yPosP) override {
//        long xIndex_plus_1;
//        long xIndex_minus_1;
//        long yIndex_plus_1;
//        long yIndex_minus_1;
//
//        xIndex_plus_1 = (xPosP + 1 == xSize ? 0 : xPosP + 1);
//        xIndex_minus_1 = (xPosP - 1 < 0 ? xSize - 1 : xPosP - 1);
//        yIndex_plus_1 = (yPosP + 1 == ySize ? 0 : yPosP + 1);
//        yIndex_minus_1 = (yPosP - 1 < 0 ? ySize - 1 : yPosP - 1);
//
//        double v_x_xx = ((*CurrentData)[xIndex_minus_1][yPosP][1] / (*CurrentData)[xIndex_minus_1][yPosP][0] +
//                         (*CurrentData)[xPosP][yPosP][1] / (*CurrentData)[xPosP][yPosP][0] * 2 +
//                         (*CurrentData)[xIndex_plus_1][yPosP][1] / (*CurrentData)[xIndex_plus_1][yPosP][0]) /
//                        (xStep * xStep);
//        double v_x_yy = ((*CurrentData)[xPosP][yIndex_minus_1][1] / (*CurrentData)[xPosP][yIndex_minus_1][0] +
//                         (*CurrentData)[xPosP][yPosP][1] / (*CurrentData)[xPosP][yPosP][0] * 2 +
//                         (*CurrentData)[xPosP][yIndex_plus_1][1] / (*CurrentData)[xPosP][yIndex_plus_1][0]) /
//                        (yStep * yStep);
//        double v_y_xx = ((*CurrentData)[xIndex_minus_1][yPosP][2] / (*CurrentData)[xIndex_minus_1][yPosP][0] +
//                         (*CurrentData)[xPosP][yPosP][2] / (*CurrentData)[xPosP][yPosP][0] * 2 +
//                         (*CurrentData)[xIndex_plus_1][yPosP][2] / (*CurrentData)[xIndex_plus_1][yPosP][0]) /
//                        (xStep * xStep);
//        double v_y_yy = ((*CurrentData)[xPosP][yIndex_minus_1][2] / (*CurrentData)[xPosP][yIndex_minus_1][0] +
//                         (*CurrentData)[xPosP][yPosP][2] / (*CurrentData)[xPosP][yPosP][0] * 2 +
//                         (*CurrentData)[xPosP][yIndex_plus_1][2] / (*CurrentData)[xPosP][yIndex_plus_1][0]) /
//                        (yStep * yStep);
//        double B_x_xx = ((*CurrentData)[xIndex_minus_1][yPosP][3] / (*CurrentData)[xIndex_minus_1][yPosP][0] +
//                         (*CurrentData)[xPosP][yPosP][3] / (*CurrentData)[xPosP][yPosP][0] * 2 +
//                         (*CurrentData)[xIndex_plus_1][yPosP][3] / (*CurrentData)[xIndex_plus_1][yPosP][0]) /
//                        (xStep * xStep);
//        double B_x_yy = ((*CurrentData)[xPosP][yIndex_minus_1][3] / (*CurrentData)[xPosP][yIndex_minus_1][0] +
//                         (*CurrentData)[xPosP][yPosP][3] / (*CurrentData)[xPosP][yPosP][0] * 2 +
//                         (*CurrentData)[xPosP][yIndex_plus_1][3] / (*CurrentData)[xPosP][yIndex_plus_1][0]) /
//                        (yStep * yStep);
//        double B_y_xx = ((*CurrentData)[xIndex_minus_1][yPosP][4] / (*CurrentData)[xIndex_minus_1][yPosP][0] +
//                         (*CurrentData)[xPosP][yPosP][4] / (*CurrentData)[xPosP][yPosP][0] * 2 +
//                         (*CurrentData)[xIndex_plus_1][yPosP][4] / (*CurrentData)[xIndex_plus_1][yPosP][0]) /
//                        (xStep * xStep);
//        double B_y_yy = ((*CurrentData)[xPosP][yIndex_minus_1][4] / (*CurrentData)[xPosP][yIndex_minus_1][0] +
//                         (*CurrentData)[xPosP][yPosP][4] / (*CurrentData)[xPosP][yPosP][0] * 2 +
//                         (*CurrentData)[xPosP][yIndex_plus_1][4] / (*CurrentData)[xPosP][yIndex_plus_1][0]) /
//                        (yStep * yStep);
//
//        return dVectorND <double> ({0,
//                                    (*CurrentData)[xPosP][yPosP][0] * (v_x_xx + v_x_yy),
//                                    (*CurrentData)[xPosP][yPosP][0] * (v_y_xx + v_y_yy),
//                                    (*CurrentData)[xPosP][yPosP][0] * (B_x_xx + B_x_yy),
//                                    (*CurrentData)[xPosP][yPosP][0] * (B_y_xx + B_y_yy)});
//    }

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
        double Fraction = 1.0 / (2.0 * M_PI * 20.0 * mStepX * 20.0 * mStepY);

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
