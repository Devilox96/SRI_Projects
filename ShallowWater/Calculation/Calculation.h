#ifndef CALCULATION_H
#define CALCULATION_H
//-----------------------------//
#include <zconf.h>
#include <iostream>
#include <vector>
#include <cmath>
//-----------------------------//
#include <QObject>
#include <QVector>
//-----------------------------//
#include "Libs/dMath/Core/dVectors.h"
#include "Libs/dMath/NumerCalc/dLaxFriedrichs.h"
//-----------------------------//
class Calculation : public QObject {
    Q_OBJECT
public:
    Calculation(uint xGridNumP, uint yGridNumP);
    ~Calculation() override = default;

    void GetSolution(uint xPosP, uint yPosP, dVectorND <double>& TargetVectorP);
    void Solve(uint yPosP);
    void DisplayData();
private:
    double g = 1;
    double B_0 = 0;
    double f_0 = 0;
    double beta = 0;

    uint xGridNum;
    uint yGridNum;

    bool ActiveGrid = false;

    std::vector <std::vector <dVectorND <double>>> Grid[2];

    std::vector <double> f;

    //----------//

    dLaxFriedrichs* Solver;

    //----------//

    void InitGrid();

    //----------//

    dVectorND <double> CalcX(const dVectorND <double>& SolutionP);
    dVectorND <double> CalcY(const dVectorND <double>& SolutionP);
    dVectorND <double> CalcF(const dVectorND <double>& SolutionP);
    dVectorND <double> CalcR(uint yPosP, const dVectorND <double>& SolutionP);
signals:
    void SendDisplayDataSignal(const QVector <QVector <double>>& DataP, double xGridStepP, double yGridStepP);
};
//-----------------------------//
#endif
