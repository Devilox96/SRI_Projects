#ifndef DRICHTMYER2D_H
#define DRICHTMYER2D_H
//-----------------------------//
#include <vector>
#include "dVectors.h" //---dMath/Core/dVectors.h---//
//-----------------------------//
class dRichtmyer2D {
public:
    explicit dRichtmyer2D(bool HomegeneousP = true) : Homogeneous(HomegeneousP) {}
    ~dRichtmyer2D() = default;

    //----------//

    void SetTimeStep(double TimeStepP);

    void SetXStep(double xStepP);
    void SetYStep(double yStepP);

    //----------//

    double GetTimeStep();

    double GetXStep();
    double GetYStep();

    //----------//

    dVectorND <double> Solve(   const std::vector <std::vector <dVectorND<double>>>& GridP,
                                long xIndexP,
                                long yIndexP,
                                long xSizeP,
                                long ySizeP);
protected:
    double TimeStep = 0.0;

    double xStep = 0.0;
    double yStep = 0.0;

    bool Homogeneous = true;

    //----------//

    virtual dVectorND <double> xFunc(const dVectorND <double>& U) = 0;
    virtual dVectorND <double> yFunc(const dVectorND <double>& U) = 0;
    virtual dVectorND <double> AbsValFunc(const dVectorND <double>& U) {
        return dVectorND <double> (1);
    }

    //----------//

    dVectorND <double> FirstStepSolve(  const dVectorND <double>& U,
                                        const dVectorND <double>& Ux_minus_1,
                                        const dVectorND <double>& Ux_plus_1,
                                        const dVectorND <double>& Uy_minus_1,
                                        const dVectorND <double>& Uy_plus_1);
};
//-----------------------------//
#endif
