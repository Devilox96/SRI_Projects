#ifndef CALCULATION_H
#define CALCULATION_H
//-----------------------------//
#include <zconf.h>
#include <iostream>
#include <vector>
#include <cmath>
//-----------------------------//
#include "../dMath/Core/dVectors.h"
//-----------------------------//
class Calculation {
public:
    Calculation(uint xGridNumP, uint yGridNumP);
    ~Calculation() = default;

    void GetSolution(uint xPosP, uint yPosP, dVectorND <double>& TargetVectorP);
private:
    double g = 1;
    double B_0 = 0;
    double f_0 = 0;
    double beta = 0;

    uint xGridNum;
    uint yGridNum;

    uint8_t ActiveGrid = 1;

    std::vector <std::vector <dVectorND <double>>> GridFirst;
    std::vector <std::vector <dVectorND <double>>> GridSecond;

    std::vector <double> f;

    dVectorND <double> X;
    dVectorND <double> Y;
    dVectorND <double> F;
    dVectorND <double> R;

    //----------//

    void InitGrid();

    //----------//

    void CalcX(const dVectorND <double>& SolutionP);
    void CalcY(const dVectorND <double>& SolutionP);
    void CalcF(const dVectorND <double>& SolutionP);
    void CalcR(uint yPosP, const dVectorND <double>& SolutionP);
};
//-----------------------------//
#endif
