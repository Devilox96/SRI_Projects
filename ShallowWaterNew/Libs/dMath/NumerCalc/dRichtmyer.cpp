#include "dRichtmyer.h"
//-----------------------------//
void dRichtmyer::SetTimeStep(double TimeStepP) {
    TimeStep = TimeStepP;
}

void dRichtmyer::SetXStep(double xStepP) {
    xStep = xStepP;
}
void dRichtmyer::SetYStep(double yStepP) {
    yStep = yStepP;
}
void dRichtmyer::SetZStep(double zStepP) {
    zStep = zStepP;
}
//-----------------------------//
double dRichtmyer::GetTimeStep() {
    return TimeStep;
}

double dRichtmyer::GetXStep() {
    return xStep;
}
double dRichtmyer::GetYStep() {
    return yStep;
}
double dRichtmyer::GetZStep() {
    return zStep;
}
//-----------------------------//
dVectorND <double> dRichtmyer::Solve1D( dVectorND <double> U,
                                        dVectorND <double> Ux_minus_1,
                                        dVectorND <double> Ux_plus_1) {
    if (TimeStep == 0.0 || xStep == 0.0) {
        std::cout << "Error: set steps!" << std::endl;
        exit(-1);
    }

    dVectorND <double> UxHalf_plusL = UxHalfVector(U, Ux_plus_1);
    dVectorND <double> UxHalf_minusL = UxHalfVector(Ux_minus_1, U);

    return  U -
            (xFunc(UxHalf_plusL) - xFunc(UxHalf_minusL)) * (TimeStep / xStep);
}
dVectorND <double> dRichtmyer::Solve2D( dVectorND <double> U,
                                        dVectorND <double> Ux_minus_1,
                                        dVectorND <double> Ux_plus_1,
                                        dVectorND <double> Uy_minus_1,
                                        dVectorND <double> Uy_plus_1) {
    if (TimeStep == 0.0 || xStep == 0.0 || yStep == 0.0) {
        std::cout << "Error: set steps!" << std::endl;
        exit(-1);
    }

    dVectorND <double> UxHalf_plusL = UxHalfVector(U, Ux_plus_1);
    dVectorND <double> UxHalf_minusL = UxHalfVector(Ux_minus_1, U);

    dVectorND <double> UyHalf_plusL = UyHalfVector(U, Uy_plus_1);
    dVectorND <double> UyHalf_minusL = UyHalfVector(Uy_minus_1, U);

    return  U -
            (xFunc(UxHalf_plusL) - xFunc(UxHalf_minusL)) * (TimeStep / xStep) -
            (yFunc(UyHalf_plusL) - yFunc(UyHalf_minusL)) * (TimeStep / yStep);
}
dVectorND <double> dRichtmyer::Solve3D( dVectorND <double> U,
                                        dVectorND <double> Ux_minus_1,
                                        dVectorND <double> Ux_plus_1,
                                        dVectorND <double> Uy_minus_1,
                                        dVectorND <double> Uy_plus_1,
                                        dVectorND <double> Uz_minus_1,
                                        dVectorND <double> Uz_plus_1) {
    if (TimeStep == 0.0 || xStep == 0.0 || yStep == 0.0 || zStep == 0.0) {
        std::cout << "Error: set steps!" << std::endl;
        exit(-1);
    }

    dVectorND <double> UxHalf_plusL = UxHalfVector(U, Ux_plus_1);
    dVectorND <double> UxHalf_minusL = UxHalfVector(Ux_minus_1, U);

    dVectorND <double> UyHalf_plusL = UyHalfVector(U, Uy_plus_1);
    dVectorND <double> UyHalf_minusL = UyHalfVector(Uy_minus_1, U);

    dVectorND <double> UzHalf_plusL = UzHalfVector(U, Uz_plus_1);
    dVectorND <double> UzHalf_minusL = UzHalfVector(Uz_minus_1, U);

    return  U -
            (xFunc(UxHalf_plusL) - xFunc(UxHalf_minusL)) * (TimeStep / xStep) -
            (yFunc(UyHalf_plusL) - yFunc(UyHalf_minusL)) * (TimeStep / yStep) -
            (zFunc(UzHalf_plusL) - zFunc(UzHalf_minusL)) * (TimeStep / zStep);
}
//-----------------------------//
dVectorND <double> dRichtmyer::UxHalfVector(const dVectorND <double>& Ux, const dVectorND <double>& Ux_plus_1) {
    return (Ux_plus_1 + Ux) / 2.0 - (xFunc(Ux_plus_1) - xFunc(Ux)) * (TimeStep / 2.0 / xStep);
}
dVectorND <double> dRichtmyer::UyHalfVector(const dVectorND <double>& Uy, const dVectorND <double>& Uy_plus_1) {
    return (Uy_plus_1 + Uy) / 2.0 - (yFunc(Uy_plus_1) - yFunc(Uy)) * (TimeStep / 2.0 / yStep);
}
dVectorND <double> dRichtmyer::UzHalfVector(const dVectorND <double>& Uz, const dVectorND <double>& Uz_plus_1) {
    return (Uz_plus_1 + Uz) / 2.0 - (zFunc(Uz_plus_1) - zFunc(Uz)) * (TimeStep / 2.0 / zStep);
}
//-----------------------------//
//-----------------------------//
//-----------------------------//
dRichtmyerSolver::dRichtmyerSolver(double TimeStepP) {
    TimeStep = TimeStepP;
}
dRichtmyerSolver::dRichtmyerSolver(double TimeStepP, double xStepP) {
    TimeStep = TimeStepP;
    xStep = xStepP;
}
dRichtmyerSolver::dRichtmyerSolver(double TimeStepP, double xStepP, double yStepP) {
    TimeStep = TimeStepP;
    xStep = xStepP;
    yStep = yStepP;
}
dRichtmyerSolver::dRichtmyerSolver(double TimeStepP, double xStepP, double yStepP, double zStepP) {
    TimeStep = TimeStepP;
    xStep = xStepP;
    yStep = yStepP;
    zStep = zStepP;
}
//-----------------------------//
dVectorND <double> dRichtmyerSolver::xFunc(const dVectorND <double>& U) {
//    return dVectorND <double>({U[1], pow(U[1], 2.0) / U[0] + 0.5 * g * pow(U[0], 2.0), U[1] * U[2] / U[0]});
    return dVectorND <double> ({U[1],
                                (pow(U[1], 2.0) - pow(U[3], 2.0)) / U[0] + 0.5 * g * pow(U[0], 2.0),
                                (U[1] * U[2] - U[3] * U[4]) / U[0],
                                0,
                                (U[1] * U[4] - U[2] * U[3]) / U[0],
                                B_0 * U[1] / U[0]});
}
dVectorND <double> dRichtmyerSolver::yFunc(const dVectorND <double>& U) {
//    return dVectorND <double>({U[2], U[1] * U[2] / U[0], pow(U[2], 2.0) / U[0] + 0.5 * g * pow(U[0], 2.0)});
    return dVectorND <double> ({U[2],
                                (U[1] * U[2] - U[3] * U[4]) / U[0],
                                (pow(U[2], 2.0) - pow(U[4], 2.0)) / U[0] + 0.5 * g * pow(U[0], 2.0),
                                (U[2] * U[3] - U[1] * U[4]) / U[0],
                                0,
                                B_0 * U[2] / U[0]});
}
dVectorND <double> dRichtmyerSolver::zFunc(const dVectorND <double>& U) {}