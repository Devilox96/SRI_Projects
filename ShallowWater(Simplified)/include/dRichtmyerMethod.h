#ifndef DRICHTMYERMETHOD_H
#define DRICHTMYERMETHOD_H
//-----------------------------//
#include "../Libs/dMath/Core/dVectors.h"
//-----------------------------//
#include <iostream>
//-----------------------------//
class dRichtmyerMethod {
public:
    dRichtmyerMethod() = default;
    ~dRichtmyerMethod() = default;

    void Solve1D(   const dVectorND <double>& CurVecP,
                    const dVectorND <double>& xPlusVecP,    const dVectorND <double>& xMinusVecP,
                    dVectorND <double>& TargetVecP);
    void Solve2D(   const dVectorND <double>& CurVecP,
                    const dVectorND <double>& xPlusVecP,    const dVectorND <double>& xMinusVecP,
                    const dVectorND <double>& yPlusVecP,    const dVectorND <double>& yMinusVecP,
                    dVectorND <double>& TargetVecP);
    void Solve3D(   const dVectorND <double>& CurVecP,
                    const dVectorND <double>& xPlusVecP,    const dVectorND <double>& xMinusVecP,
                    const dVectorND <double>& yPlusVecP,    const dVectorND <double>& yMinusVecP,
                    const dVectorND <double>& zPlusVecP,    const dVectorND <double>& zMinusVecP,
                    dVectorND <double>& TargetVecP);

    void SetTimeStep(double StepP);

    void SetXStep(double StepP);
    void SetYStep(double StepP);
    void SetZStep(double StepP);
private:
    double TimeStep = 0.0;

    double xStep = 0.0;
    double yStep = 0.0;
    double zStep = 0.0;

    virtual dVectorND <double> xFunc(const dVectorND <double>& VecP) = 0;
    virtual dVectorND <double> yFunc(const dVectorND <double>& VecP) = 0;
    virtual dVectorND <double> zFunc(const dVectorND <double>& VecP) = 0;

    dVectorND <double> xPlusHalfVec(    const dVectorND <double>& CurVecP,
                                        const dVectorND <double>& PlusVecP);
    dVectorND <double> xMinusHalfVec(   const dVectorND <double>& CurVecP,
                                        const dVectorND <double>& MinusVecP);
    dVectorND <double> yPlusHalfVec(    const dVectorND <double>& CurVecP,
                                        const dVectorND <double>& PlusVecP);
    dVectorND <double> yMinusHalfVec(   const dVectorND <double>& CurVecP,
                                        const dVectorND <double>& MinusVecP);
    dVectorND <double> zPlusHalfVec(    const dVectorND <double>& CurVecP,
                                        const dVectorND <double>& PlusVecP);
    dVectorND <double> zMinusHalfVec(   const dVectorND <double>& CurVecP,
                                        const dVectorND <double>& MinusVecP);
};
//-----------------------------//
#endif
