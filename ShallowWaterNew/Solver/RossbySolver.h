#ifndef ROSSBYSOLVER_H
#define ROSSBYSOLVER_H
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
class RossbySolver : public dRichtmyer2D <dVector <double, 3>> {
public:
    RossbySolver(double tStepTime, double tStepX, double tStepY);

    //----------//

    void initGrid(int tGridX, int tGridY);
    void setSavePath(const std::string& tPath);
    void solve();
private:
    //---Constants----//
    const double mGrav  = 9.81;
    const double mOmega = 7.2921e-05;
    const double mRad   = 6.371e+06;
    //---Constants----//

    //----Initials----//
    double mInitEnergy;
    //----Initials----//

    //---Beta-plane---//
    std::vector <double> mAngle;
    std::vector <double> mCorParam;

    double mAngle_0;
    double mCorParam_0;
    double mBetaParam;
    //---Beta-plane---//

    //------Grid------//
    int mGridX = 1;
    int mGridY = 1;
    //------Grid------//

    //-----Arrays-----//
    dGrid DataFirst;
    dGrid DataSecond;

    dGrid* CurrentData  = &DataFirst;
    dGrid* TempData     = &DataSecond;

    std::vector <dVector <double, 3>> mLeftBoundary;
    std::vector <dVector <double, 3>> mRightBoundary;
    std::vector <dVector <double, 3>> mUpperBoundary;
    std::vector <dVector <double, 3>> mLowerBoundary;

//    std::vector <std::vector <std::vector <double>>> mGradient;
    dGrid mGeography;
    //-----Arrays-----//

    //-----Saving-----//
    std::ofstream mAmpFile;
    std::ofstream mVelFileX;
    std::ofstream mVelFileY;

    std::string mSavePath = "./";
    //-----Saving-----//

    //----------//

    std::vector <std::vector <std::vector <double>>> gradient(const dGrid& tGrid);
    void boundaries();
    double getFullEnergy();
    void appendData();
    void saveData();

    dVector <double, 3> funcX(const dVector <double, 3>& tVec) override;
    dVector <double, 3> funcY(const dVector <double, 3>& tVec) override;
    dVector <double, 3> source(int tPosX, int tPosY);
};
//-----------------------------//
#endif
