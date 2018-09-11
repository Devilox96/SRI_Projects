#include "Calculation.h"
//-----------------------------//
Calculation::Calculation(uint xGridNumP, uint yGridNumP) {
    xGridNum = xGridNumP;
    yGridNum = yGridNumP;

    InitGrid();

    Solver = new dLaxFriedrichs;
}
//-----------------------------//
void Calculation::GetSolution(uint xPosP, uint yPosP, dVectorND <double>& TargetVectorP) {
    TargetVectorP[0] = Grid[ActiveGrid][xPosP][yPosP][0];
    TargetVectorP[1] = Grid[ActiveGrid][xPosP][yPosP][1] / Grid[ActiveGrid][xPosP][yPosP][0];
    TargetVectorP[2] = Grid[ActiveGrid][xPosP][yPosP][2] / Grid[ActiveGrid][xPosP][yPosP][0];
    TargetVectorP[3] = Grid[ActiveGrid][xPosP][yPosP][3] / Grid[ActiveGrid][xPosP][yPosP][0];
    TargetVectorP[4] = Grid[ActiveGrid][xPosP][yPosP][4] / Grid[ActiveGrid][xPosP][yPosP][0];
    TargetVectorP[5] = Grid[ActiveGrid][xPosP][yPosP][5];
}
void Calculation::Solve(uint yPosP) {
    for (uint i = 1; i < xGridNum - 1; i++) {
        for (uint j = 1; j < yGridNum - 1; j++) {
            Grid[!ActiveGrid][i][j] = Solver -> Solve(
                    Grid[ActiveGrid][i][j],
                    CalcX(Grid[ActiveGrid][i + 1][j]),
                    CalcX(Grid[ActiveGrid][i - 1][j]),
                    CalcY(Grid[ActiveGrid][i][j + 1]),
                    CalcY(Grid[ActiveGrid][i][j - 1]),
                    CalcF(Grid[ActiveGrid][i][j]),
                    CalcR(yPosP, Grid[ActiveGrid][i][j]),
                    0.005,
                    0.1,
                    0.1);
        }
    }

    for (uint i = 1; i < xGridNum - 1; i++) {
        Grid[!ActiveGrid][i][0] = Solver -> Solve(
                Grid[ActiveGrid][i][0],
                CalcX(Grid[ActiveGrid][i + 1][0]),
                CalcX(Grid[ActiveGrid][i - 1][0]),
                CalcY(Grid[ActiveGrid][i][1]),
                CalcY(Grid[ActiveGrid][i][yGridNum - 1]),
                CalcF(Grid[ActiveGrid][i][0]),
                CalcR(yPosP, Grid[ActiveGrid][i][0]),
                0.005,
                0.1,
                0.1);
        Grid[!ActiveGrid][i][yGridNum - 1] = Solver -> Solve(
                Grid[ActiveGrid][i][yGridNum - 1],
                CalcX(Grid[ActiveGrid][i + 1][yGridNum - 1]),
                CalcX(Grid[ActiveGrid][i - 1][yGridNum - 1]),
                CalcY(Grid[ActiveGrid][i][0]),
                CalcY(Grid[ActiveGrid][i][yGridNum - 2]),
                CalcF(Grid[ActiveGrid][i][yGridNum - 1]),
                CalcR(yPosP, Grid[ActiveGrid][i][yGridNum - 1]),
                0.005,
                0.1,
                0.1);
    }

    for (uint j = 1; j < yGridNum - 1; j++) {
        Grid[!ActiveGrid][0][j] = Solver -> Solve(
                Grid[ActiveGrid][0][j],
                CalcX(Grid[ActiveGrid][1][j]),
                CalcX(Grid[ActiveGrid][xGridNum - 1][j]),
                CalcY(Grid[ActiveGrid][0][j + 1]),
                CalcY(Grid[ActiveGrid][0][j - 1]),
                CalcF(Grid[ActiveGrid][0][j]),
                CalcR(yPosP, Grid[ActiveGrid][0][j]),
                0.005,
                0.1,
                0.1);
        Grid[!ActiveGrid][xGridNum - 1][j] = Solver -> Solve(
                Grid[ActiveGrid][xGridNum - 1][j],
                CalcX(Grid[ActiveGrid][0][j]),
                CalcX(Grid[ActiveGrid][xGridNum - 2][j]),
                CalcY(Grid[ActiveGrid][xGridNum - 1][j + 1]),
                CalcY(Grid[ActiveGrid][xGridNum - 1][j - 1]),
                CalcF(Grid[ActiveGrid][xGridNum - 1][j]),
                CalcR(yPosP, Grid[ActiveGrid][xGridNum - 1][j]),
                0.005,
                0.1,
                0.1);
    }

    Grid[!ActiveGrid][0][0] = Solver -> Solve(
            Grid[ActiveGrid][0][0],
            CalcX(Grid[ActiveGrid][1][0]),
            CalcX(Grid[ActiveGrid][xGridNum - 1][0]),
            CalcY(Grid[ActiveGrid][0][1]),
            CalcY(Grid[ActiveGrid][0][yGridNum - 1]),
            CalcF(Grid[ActiveGrid][0][0]),
            CalcR(yPosP, Grid[ActiveGrid][0][0]),
            0.005,
            0.1,
            0.1);
    Grid[!ActiveGrid][0][yGridNum - 1] = Solver -> Solve(
            Grid[ActiveGrid][0][yGridNum - 1],
            CalcX(Grid[ActiveGrid][1][yGridNum - 1]),
            CalcX(Grid[ActiveGrid][xGridNum - 1][yGridNum - 1]),
            CalcY(Grid[ActiveGrid][0][0]),
            CalcY(Grid[ActiveGrid][0][yGridNum - 2]),
            CalcF(Grid[ActiveGrid][0][yGridNum - 1]),
            CalcR(yPosP, Grid[ActiveGrid][0][yGridNum - 1]),
            0.005,
            0.1,
            0.1);
    Grid[!ActiveGrid][xGridNum - 1][0] = Solver -> Solve(
            Grid[ActiveGrid][xGridNum - 1][0],
            CalcX(Grid[ActiveGrid][0][0]),
            CalcX(Grid[ActiveGrid][xGridNum - 2][0]),
            CalcY(Grid[ActiveGrid][xGridNum - 1][1]),
            CalcY(Grid[ActiveGrid][xGridNum - 1][yGridNum - 1]),
            CalcF(Grid[ActiveGrid][xGridNum - 1][0]),
            CalcR(yPosP, Grid[ActiveGrid][xGridNum - 1][0]),
            0.005,
            0.1,
            0.1);
    Grid[!ActiveGrid][xGridNum - 1][yGridNum - 1] = Solver -> Solve(
            Grid[ActiveGrid][xGridNum - 1][yGridNum - 1],
            CalcX(Grid[ActiveGrid][0][yGridNum - 1]),
            CalcX(Grid[ActiveGrid][xGridNum - 2][yGridNum - 1]),
            CalcY(Grid[ActiveGrid][xGridNum - 1][0]),
            CalcY(Grid[ActiveGrid][xGridNum - 1][yGridNum - 2]),
            CalcF(Grid[ActiveGrid][xGridNum - 1][yGridNum - 1]),
            CalcR(yPosP, Grid[ActiveGrid][xGridNum - 1][yGridNum - 1]),
            0.005,
            0.1,
            0.1);

    ActiveGrid = !ActiveGrid;
}
//-----------------------------//
void Calculation::InitGrid() {
    Grid[0].resize(yGridNum);
    Grid[1].resize(yGridNum);

    for (uint i = 0; i < yGridNum; i++) {
        Grid[0][i].resize(xGridNum, dVectorND <double>(6));
        Grid[1][i].resize(xGridNum, dVectorND <double>(6));
    }

    for (uint j = 0; j < yGridNum; j++) {
        for (uint i = 0; i < xGridNum; i++) {
            Grid[ActiveGrid][i][j][0] = 1;
            Grid[ActiveGrid][i][j][1] = 0;
            Grid[ActiveGrid][i][j][2] = 0;
            Grid[ActiveGrid][i][j][3] = 0;
            Grid[ActiveGrid][i][j][4] = 0;
            Grid[ActiveGrid][i][j][5] = 0;
        }
    }

    f.resize(yGridNum, 0);
}
//-----------------------------//
dVectorND <double> Calculation::CalcX(const dVectorND <double>& SolutionP) {
    //---SolutionP[0] - h-----//
    //---SolutionP[1] - v_x---//
    //---SolutionP[2] - v_y---//
    //---SolutionP[3] - B_x---//
    //---SolutionP[4] - B_y---//
    //---SolutionP[5] - B_z---//

    dVectorND <double> XL(6);

    XL[0] = SolutionP[0] * SolutionP[1];
    XL[1] = SolutionP[0] * (pow(SolutionP[1], 2.0) - pow(SolutionP[3], 2.0)) + g * pow(SolutionP[0], 2.0) / 2.0;
    XL[2] = SolutionP[0] * (SolutionP[1] * SolutionP[2] - SolutionP[3] * SolutionP[4]);
    XL[3] = 0.0;
    XL[4] = SolutionP[0] * (SolutionP[1] * SolutionP[4] - SolutionP[2] * SolutionP[3]);
    XL[5] = B_0 * SolutionP[1];

    return XL;
}
dVectorND <double> Calculation::CalcY(const dVectorND <double>& SolutionP) {
    //---SolutionP[0] - h-----//
    //---SolutionP[1] - v_x---//
    //---SolutionP[2] - v_y---//
    //---SolutionP[3] - B_x---//
    //---SolutionP[4] - B_y---//
    //---SolutionP[5] - B_z---//

    dVectorND <double> YL(6);

    YL[0] = SolutionP[0] * SolutionP[2];
    YL[1] = SolutionP[0] * (SolutionP[1] * SolutionP[2] - SolutionP[3] * SolutionP[4]);
    YL[2] = SolutionP[0] * (pow(SolutionP[2], 2.0) - pow(SolutionP[4], 2.0)) + g * pow(SolutionP[0], 2.0) / 2.0;
    YL[3] = SolutionP[0] * (SolutionP[2] * SolutionP[3] - SolutionP[1] * SolutionP[4]);
    YL[4] = 0.0;
    YL[5] = B_0 * SolutionP[2];

    return YL;
}
dVectorND <double> Calculation::CalcF(const dVectorND <double>& SolutionP) {
    //---SolutionP[0] - h-----//
    //---SolutionP[1] - v_x---//
    //---SolutionP[2] - v_y---//
    //---SolutionP[3] - B_x---//
    //---SolutionP[4] - B_y---//
    //---SolutionP[5] - B_z---//

    dVectorND <double> FL(6);

    FL[0] = 0.0;
    FL[1] = B_0 * SolutionP[3] - SolutionP[0] * SolutionP[2] * f_0;
    FL[2] = B_0 * SolutionP[4] + SolutionP[0] * SolutionP[1] * f_0;
    FL[3] = -B_0 * SolutionP[1];
    FL[4] = -B_0 * SolutionP[2];
    FL[5] = 0.0;

    return FL;
}
dVectorND <double> Calculation::CalcR(uint yPosP, const dVectorND <double>& SolutionP) {
    //---SolutionP[0] - h-----//
    //---SolutionP[1] - v_x---//
    //---SolutionP[2] - v_y---//
    //---SolutionP[3] - B_x---//
    //---SolutionP[4] - B_y---//
    //---SolutionP[5] - B_z---//

    dVectorND <double> RL(6);

    RL[0] = 0.0;
    RL[1] = -SolutionP[0] * SolutionP[2] * beta * yPosP;
    RL[2] = SolutionP[0] * SolutionP[1] * beta * yPosP;
    RL[3] = 0.0;
    RL[4] = 0.0;
    RL[5] = 0.0;

    return RL;
}