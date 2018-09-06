#include "Calculation.h"
//-----------------------------//
Calculation::Calculation(uint xGridNumP, uint yGridNumP) {
    xGridNum = xGridNumP;
    yGridNum = yGridNumP;

    X.Resize(6);
    Y.Resize(6);
    F.Resize(6);
    R.Resize(6);

    InitGrid();
}
//-----------------------------//
void Calculation::GetSolution(uint xPosP, uint yPosP, dVectorND <double>& TargetVectorP) {
    if (ActiveGrid == 1) {
        TargetVectorP[0] = GridFirst[xPosP][yPosP][0];
        TargetVectorP[1] = GridFirst[xPosP][yPosP][1] / GridFirst[xPosP][yPosP][0];
        TargetVectorP[2] = GridFirst[xPosP][yPosP][2] / GridFirst[xPosP][yPosP][0];
        TargetVectorP[3] = GridFirst[xPosP][yPosP][3] / GridFirst[xPosP][yPosP][0];
        TargetVectorP[4] = GridFirst[xPosP][yPosP][4] / GridFirst[xPosP][yPosP][0];
        TargetVectorP[5] = GridFirst[xPosP][yPosP][5];
    } else {
        TargetVectorP[0] = GridSecond[xPosP][yPosP][0];
        TargetVectorP[1] = GridSecond[xPosP][yPosP][1] / GridSecond[xPosP][yPosP][0];
        TargetVectorP[2] = GridSecond[xPosP][yPosP][2] / GridSecond[xPosP][yPosP][0];
        TargetVectorP[3] = GridSecond[xPosP][yPosP][3] / GridSecond[xPosP][yPosP][0];
        TargetVectorP[4] = GridSecond[xPosP][yPosP][4] / GridSecond[xPosP][yPosP][0];
        TargetVectorP[5] = GridSecond[xPosP][yPosP][5];
    }
}
//-----------------------------//
void Calculation::InitGrid() {
    GridFirst.resize(yGridNum);
    GridSecond.resize(yGridNum);

    for (uint i = 0; i < yGridNum; i++) {
        GridFirst[i].resize(xGridNum, dVectorND <double>(6));
        GridSecond[i].resize(xGridNum, dVectorND <double>(6));
    }

    for (uint j = 0; j < yGridNum; j++) {
        for (uint i = 0; i < xGridNum; i++) {
            GridFirst[i][j][0] = 1;
            GridFirst[i][j][1] = 0;
            GridFirst[i][j][2] = 0;
            GridFirst[i][j][3] = 0;
            GridFirst[i][j][4] = 0;
            GridFirst[i][j][5] = 0;
        }
    }

    f.resize(yGridNum, 0);
}
//-----------------------------//
void Calculation::CalcX(const dVectorND <double>& SolutionP) {
    //---SolutionP[0] - h-----//
    //---SolutionP[1] - v_x---//
    //---SolutionP[2] - v_y---//
    //---SolutionP[3] - B_x---//
    //---SolutionP[4] - B_y---//
    //---SolutionP[5] - B_z---//

    X[0] = SolutionP[0] * SolutionP[1];
    X[1] = SolutionP[0] * (pow(SolutionP[1], 2.0) - pow(SolutionP[3], 2.0)) + g * pow(SolutionP[0], 2.0) / 2.0;
    X[2] = SolutionP[0] * (SolutionP[1] * SolutionP[2] - SolutionP[3] * SolutionP[4]);
    X[3] = 0.0;
    X[4] = SolutionP[0] * (SolutionP[1] * SolutionP[4] - SolutionP[2] * SolutionP[3]);
    X[5] = B_0 * SolutionP[1];
}
void Calculation::CalcY(const dVectorND <double>& SolutionP) {
    //---SolutionP[0] - h-----//
    //---SolutionP[1] - v_x---//
    //---SolutionP[2] - v_y---//
    //---SolutionP[3] - B_x---//
    //---SolutionP[4] - B_y---//
    //---SolutionP[5] - B_z---//

    Y[0] = SolutionP[0] * SolutionP[2];
    Y[1] = SolutionP[0] * (SolutionP[1] * SolutionP[2] - SolutionP[3] * SolutionP[4]);
    Y[2] = SolutionP[0] * (pow(SolutionP[2], 2.0) - pow(SolutionP[4], 2.0)) + g * pow(SolutionP[0], 2.0) / 2.0;
    Y[3] = SolutionP[0] * (SolutionP[2] * SolutionP[3] - SolutionP[1] * SolutionP[4]);
    Y[4] = 0.0;
    Y[5] = B_0 * SolutionP[2];
}
void Calculation::CalcF(const dVectorND <double>& SolutionP) {
    //---SolutionP[0] - h-----//
    //---SolutionP[1] - v_x---//
    //---SolutionP[2] - v_y---//
    //---SolutionP[3] - B_x---//
    //---SolutionP[4] - B_y---//
    //---SolutionP[5] - B_z---//

    F[0] = 0.0;
    F[1] = B_0 * SolutionP[3] - SolutionP[0] * SolutionP[2] * f_0;
    F[2] = B_0 * SolutionP[4] + SolutionP[0] * SolutionP[1] * f_0;
    F[3] = -B_0 * SolutionP[1];
    F[4] = -B_0 * SolutionP[2];
    F[5] = 0.0;
}
void Calculation::CalcR(uint yPosP, const dVectorND <double>& SolutionP) {
    //---SolutionP[0] - h-----//
    //---SolutionP[1] - v_x---//
    //---SolutionP[2] - v_y---//
    //---SolutionP[3] - B_x---//
    //---SolutionP[4] - B_y---//
    //---SolutionP[5] - B_z---//

    R[0] = 0.0;
    R[1] = -SolutionP[0] * SolutionP[2] * beta * yPosP;
    R[2] = SolutionP[0] * SolutionP[1] * beta * yPosP;
    R[3] = 0.0;
    R[4] = 0.0;
    R[5] = 0.0;
}