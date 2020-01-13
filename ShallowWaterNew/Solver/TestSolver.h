#ifndef TESTSOLVER_H
#define TESTSOLVER_H
//-----------------------------//
#include <cmath>
#include <vector>
#include <fstream>
#include <random>
//-----------------------------//
#include "../Libs/dMath/Core/dVector.h"
#include "../Libs/dMath/NumerCalc/dRichtmyer2D.h"
//-----------------------------//
using dGrid = std::vector <std::vector <dVector <double, 3>>>;
//-----------------------------//
template <class T>
class ClippedRichtmyer {
public:
    ClippedRichtmyer() = default;
    ~ClippedRichtmyer() = default;

    //----------//

    void setTimeStep(double tStep) {
        mStepTime = tStep;

        if (mStepX != 0) {
            mRatioX = mStepTime / mStepX;
        }

        if (mStepY != 0) {
            mRatioY = mStepTime / mStepY;
        }
    }
    void setXStep(double tStep) {
        mStepX = tStep;

        if (mStepX != 0) {
            mRatioX = mStepTime / mStepX;
        }
    }
    void setYStep(double tStep) {
        mStepY = tStep;

        if (mStepY != 0) {
            mRatioY = mStepTime / mStepY;
        }
    }

    //----------//

    T solve(const T& tOffsetZero,
            const T& tX_min_2_Y,        const T& tX_pl_2_Y,
            const T& tX_Y_min_2,        const T& tX_Y_pl_2,
            const T& tX_pl_1_Y_pl_1,    const T& tX_min_1_Y_min_1,
            const T& tX_pl_1_Y_min_1,   const T& tX_min_1_Y_pl_1) {
        if (mStepTime <= 0 || mStepX <= 0 || mStepY <= 0) {
            std::cout << "ERROR! dRichtmyer2D: Incorrect steps!" << std::endl;
        }

        T HalfMinusX = firstStep(tX_min_2_Y, tOffsetZero, tX_min_1_Y_min_1, tX_min_1_Y_pl_1);
        T HalfPlusX = firstStep(tOffsetZero, tX_pl_2_Y, tX_pl_1_Y_min_1, tX_pl_1_Y_pl_1);
        T HalfMinusY = firstStep(tX_min_1_Y_min_1, tX_pl_1_Y_min_1, tX_Y_min_2, tOffsetZero);
        T HalfPlusY = firstStep(tX_min_1_Y_pl_1, tX_pl_1_Y_pl_1, tOffsetZero, tX_Y_pl_2);

        T Solution = tOffsetZero -
                     (funcX(HalfPlusX) - funcX(HalfMinusX)) * mRatioX -
                     (funcY(HalfPlusY) - funcY(HalfMinusY)) * mRatioY;

        return Solution;
    }
protected:
    virtual T funcX(const T& tVal) {
        return tVal;
    }
    virtual T funcY(const T& tVal) {
        return tVal;
    }

    //----------//

    double mStepTime = 0.0;

    double mStepX = 0.0;
    double mStepY = 0.0;

    double mRatioX = 0.0;   //---mStepTime / mStepX---//
    double mRatioY = 0.0;   //---mStepTime / mStepY---//

    //----------//

    T firstStep(const T& tOffsetMinusX, const T& tOffsetPlusX,
                const T& tOffsetMinusY, const T& tOffsetPlusY) {
        return  (tOffsetPlusX + tOffsetMinusX + tOffsetPlusY + tOffsetMinusY) / 4 -
                ((funcX(tOffsetPlusX) - funcX(tOffsetMinusX)) * mRatioX +
                 (funcY(tOffsetPlusY) - funcY(tOffsetMinusY)) * mRatioY) / 2;
    }
};
//-----------------------------//
class TestSolver : public dRichtmyer2D <dVector <double, 3>> {
public:
    TestSolver();

    //----------//

    void setSavePath(const std::string& tPath);
    void openFiles();
    void solve();
private:
    //---Constants----//
    const double mGrav  = 9.81;
    const double mCorParam_0 = 1.0e-04;
    const double mBetaParam = 1.6e-11;
    //---Constants----//

    //----Initials----//
    double mInitEnergy;
    //----Initials----//

    //---Beta-plane---//
    std::vector <double> mCorParam;
    //---Beta-plane---//

    //------Grid------//
    int mGridX = 254;
//    int mGridX = 15;
    int mGridY = 50;
//    int mGridY = 15;
    //------Grid------//

    //-----Arrays-----//
    dGrid mDataFirst;
    dGrid mDataSecond;

    dGrid* CurrentData  = &mDataFirst;
    dGrid* TempData     = &mDataSecond;

    std::vector <std::vector <double>> mGeography;

    dGrid mid_xt;
    dGrid mid_yt;

    std::vector <std::vector <double>> Ux;
    std::vector <std::vector <double>> Uy;

    std::vector <std::vector <double>> Vx;
    std::vector <std::vector <double>> Vy;

    std::vector <std::vector <double>> Ux_mid_xt;
    std::vector <std::vector <double>> Uy_mid_yt;

    std::vector <std::vector <double>> Vx_mid_xt;
    std::vector <std::vector <double>> Vy_mid_yt;
    //-----Arrays-----//

    //-----Saving-----//
    std::ofstream mAmpFile;
    std::ofstream mVelFileX;
    std::ofstream mVelFileY;

    std::string mSavePath = "./";
    //-----Saving-----//

    //----------//

    void initGrid();
    void initGeography();
    void initCoriolis();
    void initConditions();

//    double getFullEnergy();
    void appendData();
    void saveData();
//
//    dVector <double, 3> funcX(const dVector <double, 3>& tVec) override;
//    dVector <double, 3> funcY(const dVector <double, 3>& tVec) override;
    dVector <double, 3> source(int tPosX, int tPosY);
};
//-----------------------------//
#endif
