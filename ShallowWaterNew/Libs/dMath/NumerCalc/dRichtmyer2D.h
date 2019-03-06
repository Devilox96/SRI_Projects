#ifndef DRICHTMYER2D_H
#define DRICHTMYER2D_H
//-----------------------------//
#include <vector>
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

    //----------//

    T Solve(const std::vector <std::vector <T>>& GridP,
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

        T Ux_minus_1L = FirstStepSolve(GridP[xIndex_minus_1][yIndexP],
                                                        GridP[xIndex_minus_2][yIndexP],
                                                        GridP[xIndexP][yIndexP],
                                                        GridP[xIndex_minus_1][yIndex_minus_1],
                                                        GridP[xIndex_minus_1][yIndex_plus_1]);
        T Ux_plus_1L = FirstStepSolve( GridP[xIndex_plus_1][yIndexP],
                                                        GridP[xIndexP][yIndexP],
                                                        GridP[xIndex_plus_2][yIndexP],
                                                        GridP[xIndex_plus_1][yIndex_minus_1],
                                                        GridP[xIndex_plus_1][yIndex_plus_1]);
        T Uy_minus_1L = FirstStepSolve(GridP[xIndexP][yIndex_minus_1],
                                                        GridP[xIndex_minus_1][yIndex_minus_1],
                                                        GridP[xIndex_plus_1][yIndex_minus_1],
                                                        GridP[xIndexP][yIndex_minus_2],
                                                        GridP[xIndexP][yIndexP]);
        T Uy_plus_1L = FirstStepSolve( GridP[xIndexP][yIndex_plus_1],
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
protected:
    double TimeStep = 0.0;

    double xStep = 0.0;
    double yStep = 0.0;

    bool Homogeneous = true;

    //----------//

    virtual T xFunc(const T& U) = 0;
    virtual T yFunc(const T& U) = 0;
    virtual T AbsValFunc(const T& U) {
        return dVectorND <double> (1);
    }

    //----------//

    T FirstStepSolve(   const T& U,
                        const T& Ux_minus_1,
                        const T& Ux_plus_1,
                        const T& Uy_minus_1,
                        const T& Uy_plus_1) {
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
};
//-----------------------------//
#endif
