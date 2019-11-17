#ifndef DRICHTMYER2D_H
#define DRICHTMYER2D_H
//-----------------------------//
#include <vector>
#include <algorithm>
//-----------------------------//
#include "dVectors.h" //---dMath/Core/dVectors.h---//
//-----------------------------//
template <class T>
class dRichtmyer2D {
public:
    explicit dRichtmyer2D(bool HomegeneousP = true) : Homogeneous(HomegeneousP) {}
    ~dRichtmyer2D() = default;

    //----------//

    void SetTimeStep(double TimeStepP) {
        TimeStep = TimeStepP;
    }

    void SetXStep(double xStepP) {
        xStep = xStepP;
    }
    void SetYStep(double yStepP) {
        yStep = yStepP;
    }

    //----------//

    double GetTimeStep() {
        return TimeStep;
    }

    double GetXStep() {
        return xStep;
    }
    double GetYStep() {
        return yStep;
    }
protected:
    double TimeStep = 0.0;

    double xStep = 0.0;
    double yStep = 0.0;

    bool Homogeneous = true;

    //----------//

    virtual T xFunc(const T& U) = 0;
    virtual T yFunc(const T& U) = 0;
    virtual T AbsValFunc(int xPosP, int yPosP) {
        return dVectorND <double> (1);
    }
    virtual T Viscosity(int xPosP, int yPosP) {
        return dVectorND <double> (1);
    }

    //----------//

    T SecondStepSolve(  const T& tOffsetZero,
                        const T& tX_min_2_Y,        const T& tX_pl_2_Y,
                        const T& tX_Y_min_2,        const T& tX_Y_pl_2,
                        const T& tX_pl_1_Y_pl_1,    const T& tX_min_1_Y_min_1,
                        const T& tX_pl_1_Y_min_1,   const T& tX_min_1_Y_pl_1,
                        int xPos, int yPos) {
        T Ux_minus_1L = FirstStepSolve(tX_min_2_Y, tOffsetZero, tX_min_1_Y_min_1, tX_min_1_Y_pl_1, xPos, yPos);
        T Ux_plus_1L = FirstStepSolve(tOffsetZero, tX_pl_2_Y, tX_pl_1_Y_min_1, tX_pl_1_Y_pl_1, xPos, yPos);
        T Uy_minus_1L = FirstStepSolve(tX_min_1_Y_min_1, tX_pl_1_Y_min_1, tX_Y_min_2, tOffsetZero, xPos, yPos);
        T Uy_plus_1L = FirstStepSolve(tX_min_1_Y_pl_1, tX_pl_1_Y_pl_1, tOffsetZero, tX_Y_pl_2, xPos, yPos);

        if (Homogeneous) {
            return  tOffsetZero -
                    TimeStep / xStep * (xFunc(Ux_plus_1L) - xFunc(Ux_minus_1L)) -
                    TimeStep / yStep * (yFunc(Uy_plus_1L) - yFunc(Uy_minus_1L));
        } else {
            return  tOffsetZero -
                    TimeStep / xStep * (xFunc(Ux_plus_1L) - xFunc(Ux_minus_1L)) -
                    TimeStep / yStep * (yFunc(Uy_plus_1L) - yFunc(Uy_minus_1L)) +
                    TimeStep * AbsValFunc(xPos, yPos) -
                    TimeStep * Viscosity(xPos, yPos) * 0.001;
        }
    }

    T FirstStepSolve(const T& Ux_minus_1, const T& Ux_plus_1, const T& Uy_minus_1, const T& Uy_plus_1, int xStepP, int yStepP) {
        if (Homogeneous) {
            return  0.25 * (Ux_plus_1 + Ux_minus_1 + Uy_plus_1 + Uy_minus_1) -
                    TimeStep / (2.0 * xStep) * (xFunc(Ux_plus_1) - xFunc(Ux_minus_1)) -
                    TimeStep / (2.0 * yStep) * (yFunc(Uy_plus_1) - yFunc(Uy_minus_1));
        } else {
            return  0.25 * (Ux_plus_1 + Ux_minus_1 + Uy_plus_1 + Uy_minus_1) -
                    TimeStep / (2.0 * xStep) * (xFunc(Ux_plus_1) - xFunc(Ux_minus_1)) -
                    TimeStep / (2.0 * yStep) * (yFunc(Uy_plus_1) - yFunc(Uy_minus_1)) +
                    TimeStep * AbsValFunc(xStepP, yStepP) -
                    TimeStep * Viscosity(xStepP, yStepP) * 0.001;
        }
    }
};
//-----------------------------//
#endif
