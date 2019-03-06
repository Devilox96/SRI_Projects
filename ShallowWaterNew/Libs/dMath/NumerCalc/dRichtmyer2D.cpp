#include "dRichtmyer2D.h"
//-----------------------------//
void dRichtmyer2D::SetTimeStep(double TimeStepP) {
    TimeStep = TimeStepP;
}

void dRichtmyer2D::SetXStep(double xStepP) {
    xStep = xStepP;
}
void dRichtmyer2D::SetYStep(double yStepP) {
    yStep = yStepP;
}
//-----------------------------//
double dRichtmyer2D::GetTimeStep() {
    return TimeStep;
}

double dRichtmyer2D::GetXStep() {
    return xStep;
}
double dRichtmyer2D::GetYStep() {
    return yStep;
}
//-----------------------------//
dVectorND <double> dRichtmyer2D::Solve( const std::vector <std::vector <dVectorND <double>>>& GridP,
                                        long xIndexP,
                                        long yIndexP,
                                        long xSizeP,
                                        long ySizeP) {
    long xIndex_plus_1 = (xIndexP + 1 == xSizeP ? 0 : xIndexP + 1);
    long xIndex_plus_2 = (xIndexP + 2 >= xSizeP ? xSizeP - xIndexP : xIndexP + 2);
    long xIndex_minus_1 = (xIndexP - 1 < 0 ? xSizeP - 1 : xIndexP - 1);
    long xIndex_minus_2 = (xIndexP - 2 < 0 ? xSizeP + xIndexP - 2 : xIndexP - 2);
    long yIndex_plus_1 = (yIndexP + 1 == ySizeP ? 0 : yIndexP + 1);
    long yIndex_plus_2 = (yIndexP + 2 >= ySizeP ? ySizeP - yIndexP : yIndexP + 2);
    long yIndex_minus_1 = (yIndexP - 1 < 0 ? ySizeP - 1 : yIndexP - 1);
    long yIndex_minus_2 = (yIndexP - 2 < 0 ? ySizeP + yIndexP - 2 : yIndexP - 2);

    dVectorND <double> Ux_minus_1L = FirstStepSolve(GridP[xIndex_minus_1][yIndexP],
                                                    GridP[xIndex_minus_2][yIndexP],
                                                    GridP[xIndexP][yIndexP],
                                                    GridP[xIndex_minus_1][yIndex_minus_1],
                                                    GridP[xIndex_minus_1][yIndex_plus_1]);
    dVectorND <double> Ux_plus_1L = FirstStepSolve( GridP[xIndex_plus_1][yIndexP],
                                                    GridP[xIndexP][yIndexP],
                                                    GridP[xIndex_plus_2][yIndexP],
                                                    GridP[xIndex_plus_1][yIndex_minus_1],
                                                    GridP[xIndex_plus_1][yIndex_plus_1]);
    dVectorND <double> Uy_minus_1L = FirstStepSolve(GridP[xIndexP][yIndex_minus_1],
                                                    GridP[xIndex_minus_1][yIndex_minus_1],
                                                    GridP[xIndex_plus_1][yIndex_minus_1],
                                                    GridP[xIndexP][yIndex_minus_2],
                                                    GridP[xIndexP][yIndexP]);
    dVectorND <double> Uy_plus_1L = FirstStepSolve( GridP[xIndexP][yIndex_plus_1],
                                                    GridP[xIndex_minus_1][yIndex_plus_1],
                                                    GridP[xIndex_plus_1][yIndex_plus_1],
                                                    GridP[xIndexP][yIndexP],
                                                    GridP[xIndexP][yIndex_plus_2]);

    if (Homogeneous) {
        return  GridP[xIndexP][yIndexP] -
                TimeStep / xStep * (xFunc(Ux_plus_1L) - xFunc(Ux_minus_1L)) -
                TimeStep / yStep * (yFunc(Uy_plus_1L) - yFunc(Uy_minus_1L));
    } else {
        return  GridP[xIndexP][yIndexP] -
                TimeStep / xStep * (xFunc(Ux_plus_1L) - xFunc(Ux_minus_1L)) -
                TimeStep / yStep * (yFunc(Uy_plus_1L) - yFunc(Uy_minus_1L)) -
                TimeStep * AbsValFunc(GridP[xIndexP][yIndexP]);
    }
}

dVectorND <double> dRichtmyer2D::FirstStepSolve(const dVectorND <double>& U,
                                                const dVectorND <double>& Ux_minus_1,
                                                const dVectorND <double>& Ux_plus_1,
                                                const dVectorND <double>& Uy_minus_1,
                                                const dVectorND <double>& Uy_plus_1) {
    if (Homogeneous) {
        return  0.25 * (Ux_plus_1 + Ux_minus_1 + Uy_plus_1 + Uy_minus_1) -
                TimeStep / (2.0 * xStep) * (xFunc(Ux_plus_1) - xFunc(Ux_minus_1)) -
                TimeStep / (2.0 * yStep) * (yFunc(Uy_plus_1) - yFunc(Uy_minus_1));
    } else {
        return  0.25 * (Ux_plus_1 + Ux_minus_1 + Uy_plus_1 + Uy_minus_1) -
                TimeStep / (2.0 * xStep) * (xFunc(Ux_plus_1) - xFunc(Ux_minus_1)) -
                TimeStep / (2.0 * yStep) * (yFunc(Uy_plus_1) - yFunc(Uy_minus_1)) -
                TimeStep * AbsValFunc(U);
    }
}