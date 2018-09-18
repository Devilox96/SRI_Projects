#include "dRichtmyerMethod.h"
//-----------------------------//
void dRichtmyerMethod::Solve1D( const dVectorND <double>& CurVecP,
                                const dVectorND <double>& xPlusVecP,    const dVectorND <double>& xMinusVecP,
                                dVectorND <double>& TargetVecP) {
    if (TimeStep <= 0.0) {
        std::cout << "Error (dRichtmyerMethod): time step is invalid or was not set!" << std::endl;

        return;
    }

    if (xStep <= 0.0) {
        std::cout << "Error (dRichtmyerMethod): X step is invalid or was not set!" << std::endl;

        return;
    }

    TargetVecP =    CurVecP -
                    (xPlusHalfVec(CurVecP, xPlusVecP) - xMinusHalfVec(CurVecP, xMinusVecP)) * (TimeStep / xStep);
}
void dRichtmyerMethod::Solve2D( const dVectorND <double>& CurVecP,
                                const dVectorND <double>& xPlusVecP,    const dVectorND <double>& xMinusVecP,
                                const dVectorND <double>& yPlusVecP,    const dVectorND <double>& yMinusVecP,
                                dVectorND <double>& TargetVecP) {
    if (TimeStep <= 0.0) {
        std::cout << "Error (dRichtmyerMethod): time step is invalid or was not set!" << std::endl;

        return;
    }

    if (xStep <= 0.0) {
        std::cout << "Error (dRichtmyerMethod): X step is invalid or was not set!" << std::endl;

        return;
    }
    if (yStep <= 0.0) {
        std::cout << "Error (dRichtmyerMethod): Y step is invalid or was not set!" << std::endl;

        return;
    }

    TargetVecP =    CurVecP -
                    (xPlusHalfVec(CurVecP, xPlusVecP) - xMinusHalfVec(CurVecP, xMinusVecP)) * (TimeStep / xStep) -
                    (yPlusHalfVec(CurVecP, yPlusVecP) - yMinusHalfVec(CurVecP, yMinusVecP)) * (TimeStep / yStep);
}
void dRichtmyerMethod::Solve3D( const dVectorND <double>& CurVecP,
                                const dVectorND <double>& xPlusVecP,    const dVectorND <double>& xMinusVecP,
                                const dVectorND <double>& yPlusVecP,    const dVectorND <double>& yMinusVecP,
                                const dVectorND <double>& zPlusVecP,    const dVectorND <double>& zMinusVecP,
                                dVectorND <double>& TargetVecP) {
    if (TimeStep <= 0.0) {
        std::cout << "Error (dRichtmyerMethod): time step is invalid or was not set!" << std::endl;

        return;
    }

    if (xStep <= 0.0) {
        std::cout << "Error (dRichtmyerMethod): X step is invalid or was not set!" << std::endl;

        return;
    }
    if (yStep <= 0.0) {
        std::cout << "Error (dRichtmyerMethod): Y step is invalid or was not set!" << std::endl;

        return;
    }
    if (zStep <= 0.0) {
        std::cout << "Error (dRichtmyerMethod): Z step is invalid or was not set!" << std::endl;

        return;
    }

    TargetVecP =    CurVecP -
                    (xPlusHalfVec(CurVecP, xPlusVecP) - xMinusHalfVec(CurVecP, xMinusVecP)) * (TimeStep / xStep) -
                    (yPlusHalfVec(CurVecP, yPlusVecP) - yMinusHalfVec(CurVecP, yMinusVecP)) * (TimeStep / yStep) -
                    (zPlusHalfVec(CurVecP, zPlusVecP) - zMinusHalfVec(CurVecP, zMinusVecP)) * (TimeStep / zStep);
}

void dRichtmyerMethod::SetTimeStep(double StepP) {
    TimeStep = StepP;
}

void dRichtmyerMethod::SetXStep(double StepP) {
    xStep = StepP;
}
void dRichtmyerMethod::SetYStep(double StepP) {
    yStep = StepP;
}
void dRichtmyerMethod::SetZStep(double StepP) {
    zStep = StepP;
}
//-----------------------------//
dVectorND <double> dRichtmyerMethod::xPlusHalfVec(  const dVectorND <double>& CurVecP,
                                                    const dVectorND <double>& PlusVecP) {
    return (PlusVecP + CurVecP) / 2.0 - (xFunc(PlusVecP) - xFunc(CurVecP)) * (TimeStep / (xStep * 2.0));
}
dVectorND <double> dRichtmyerMethod::xMinusHalfVec( const dVectorND <double>& CurVecP,
                                                    const dVectorND <double>& MinusVecP) {
    return (CurVecP + MinusVecP) / 2.0 - (xFunc(CurVecP) - xFunc(MinusVecP)) * (TimeStep / (xStep * 2.0));
}
dVectorND <double> dRichtmyerMethod::yPlusHalfVec(  const dVectorND <double>& CurVecP,
                                                    const dVectorND <double>& PlusVecP) {
    return (PlusVecP + CurVecP) / 2.0 - (yFunc(PlusVecP) - yFunc(CurVecP)) * (TimeStep / (yStep * 2.0));
}
dVectorND <double> dRichtmyerMethod::yMinusHalfVec( const dVectorND <double>& CurVecP,
                                                    const dVectorND <double>& MinusVecP) {
    return (CurVecP + MinusVecP) / 2.0 - (yFunc(CurVecP) - yFunc(MinusVecP)) * (TimeStep / (yStep * 2.0));
}
dVectorND <double> dRichtmyerMethod::zPlusHalfVec(  const dVectorND <double>& CurVecP,
                                                    const dVectorND <double>& PlusVecP) {
    return (PlusVecP + CurVecP) / 2.0 - (zFunc(PlusVecP) - zFunc(CurVecP)) * (TimeStep / (zStep * 2.0));
}
dVectorND <double> dRichtmyerMethod::zMinusHalfVec( const dVectorND <double>& CurVecP,
                                                    const dVectorND <double>& MinusVecP) {
    return (CurVecP + MinusVecP) / 2.0 - (zFunc(CurVecP) - zFunc(MinusVecP)) * (TimeStep / (zStep * 2.0));
}