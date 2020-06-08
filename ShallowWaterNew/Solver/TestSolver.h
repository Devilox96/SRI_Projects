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
#include "../Libs/dMath/NumerCalc/dLaxFriedrichs2D.h"
//-----------------------------//
using dGrid = std::vector <std::vector <dVector <double, 5>>>;
//-----------------------------//
class TestSolver : public dRichtmyer2D <dVector <double, 5>> {
public:
    TestSolver();

    //----------//

    void setSavePath(const std::string& tPath);
    void openFiles();
    void solveCustom();
private:
    //---Constants----//
    const double mGrav  = 9.81;
    const double mCorParam_0 = 1.0e-04;
    const double mBetaParam = 1.6e-11;
//    const double mBetaParam = 0.0;
    //---Constants----//

    //----Initials----//
    double mInitEnergy;

    //---0 - flat-------//
    //---1 - Earth------//
    //---2 - random-----//
    //---3 - gaussian---//

    uint8_t mGeographyType = 3;

    int mDaysToCalc = 1;
    int mSaveInterval = 1;  //---Hours---//
//    int mSaveInterval = 6;  //---Hours---//

    //----Initials----//

    //---Beta-plane---//
    std::vector <double> mCorParam;
    //---Beta-plane---//

    //-----Field------//
    std::vector <double> mVertField;
    std::vector <double> mHorizFieldY;
    //-----Field------//

    //------Grid------//
    int mGridX = 254;
    int mGridY = 50;
    //------Grid------//

    //-----Arrays-----//
    dGrid mDataFirst;
    dGrid mDataSecond;

    dGrid* CurrentData  = &mDataFirst;
    dGrid* TempData     = &mDataSecond;

    std::vector <std::vector <double>> mGeography;
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
    void initFields();
    void initConditions();

    double getFullEnergy();
    void appendData();
    void saveData();
    double getVelocity();
    void adjustTimeStep();

    void updateBoundaries();

    dVector <double, 5> funcX(const dVector <double, 5>& tVec) override;
    dVector <double, 5> funcY(const dVector <double, 5>& tVec) override;
    dVector <double, 5> source(int tPosX, int tPosY);
    dVector <double, 5> viscosity(int tPosX, int tPosY);
    dVector <double, 5> artVisc(int tPosX, int tPosY, double tFirstParam, double tSecondParam);

    //----------//

    dVector <double, 5> WENO_2(
            const dVector <double, 5>& tVal_minus,
            const dVector <double, 5>& tVal,
            const dVector <double, 5>& tVal_plus,
            bool tNegFlux);
};
//-----------------------------//
#endif
