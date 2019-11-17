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
    explicit dRichtmyer2D(unsigned long xSizeP = 1, unsigned long ySizeP = 1, bool HomegeneousP = true) :
                    xSize(xSizeP), ySize(ySizeP), Homogeneous(HomegeneousP) {}
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

    void Solve() {
        long xIndex_plus_1;
        long xIndex_plus_2;
        long xIndex_minus_1;
        long xIndex_minus_2;
        long yIndex_plus_1;
        long yIndex_plus_2;
        long yIndex_minus_1;
        long yIndex_minus_2;
        
        for (int i = 0; i < xSize; i++) {
            for (int j = 0; j < ySize; j++) {
                xIndex_plus_1 = (i + 1 == xSize ? 0 : i + 1);
                xIndex_plus_2 = (i + 2 >= xSize ? i + 2 - xSize : i + 2);
                xIndex_minus_1 = (i - 1 < 0 ? xSize - 1 : i - 1);
                xIndex_minus_2 = (i - 2 < 0 ? xSize + i - 2 : i - 2);
                yIndex_plus_1 = (j + 1 == ySize ? 0 : j + 1);
                yIndex_plus_2 = (j + 2 >= ySize ? j + 2 - ySize : j + 2);
                yIndex_minus_1 = (j - 1 < 0 ? ySize - 1 : j - 1);
                yIndex_minus_2 = (j - 2 < 0 ? ySize + j - 2 : j - 2);

                T Ux_minus_1L = FirstStepSolve((*CurrentData)[xIndex_minus_1][j],
                                               (*CurrentData)[xIndex_minus_2][j],
                                               (*CurrentData)[i][j],
                                               (*CurrentData)[xIndex_minus_1][yIndex_minus_1],
                                               (*CurrentData)[xIndex_minus_1][yIndex_plus_1], i, j);
                T Ux_plus_1L = FirstStepSolve( (*CurrentData)[xIndex_plus_1][j],
                                               (*CurrentData)[i][j],
                                               (*CurrentData)[xIndex_plus_2][j],
                                               (*CurrentData)[xIndex_plus_1][yIndex_minus_1],
                                               (*CurrentData)[xIndex_plus_1][yIndex_plus_1], i, j);
                T Uy_minus_1L = FirstStepSolve((*CurrentData)[i][yIndex_minus_1],
                                               (*CurrentData)[xIndex_minus_1][yIndex_minus_1],
                                               (*CurrentData)[xIndex_plus_1][yIndex_minus_1],
                                               (*CurrentData)[i][yIndex_minus_2],
                                               (*CurrentData)[i][j], i, j);
                T Uy_plus_1L = FirstStepSolve( (*CurrentData)[i][yIndex_plus_1],
                                               (*CurrentData)[xIndex_minus_1][yIndex_plus_1],
                                               (*CurrentData)[xIndex_plus_1][yIndex_plus_1],
                                               (*CurrentData)[i][j],
                                               (*CurrentData)[i][yIndex_plus_2], i, j);

                if (Homogeneous) {
                    (*TempData)[i][j] = (*CurrentData)[i][j] -
                                        TimeStep / xStep * (xFunc(Ux_plus_1L) - xFunc(Ux_minus_1L)) -
                                        TimeStep / yStep * (yFunc(Uy_plus_1L) - yFunc(Uy_minus_1L));
                } else {
                    (*TempData)[i][j] = (*CurrentData)[i][j] -
                                        TimeStep / xStep * (xFunc(Ux_plus_1L) - xFunc(Ux_minus_1L)) -
                                        TimeStep / yStep * (yFunc(Uy_plus_1L) - yFunc(Uy_minus_1L)) +
                                        TimeStep * AbsValFunc(i, j) -
                                        TimeStep * Viscosity(i, j) * 0.001;
                }
            }
        }

        std::swap(CurrentData, TempData);
    }
protected:
    double TimeStep = 0.0;

    double xStep = 0.0;
    double yStep = 0.0;

    bool Homogeneous = true;
    
    unsigned long xSize;
    unsigned long ySize;
    
    std::vector <std::vector <T>> DataFirst;
    std::vector <std::vector <T>> DataSecond;

    std::vector <std::vector <T>>* CurrentData = &DataFirst;
    std::vector <std::vector <T>>* TempData = &DataSecond;

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

    T FirstStepSolve(const T& U, const T& Ux_minus_1, const T& Ux_plus_1, const T& Uy_minus_1, const T& Uy_plus_1, int xStepP, int yStepP) {
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
