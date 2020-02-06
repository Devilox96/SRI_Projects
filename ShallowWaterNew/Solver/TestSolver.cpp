#include "TestSolver.h"
//-----------------------------//
TestSolver::TestSolver() {
    mStepTime = 60.0;

    mStepX = 100.0e+03;
    mStepY = 100.0e+03;
    
    mRatioX = mStepTime / mStepX;
    mRatioY = mStepTime / mStepY;

    initGrid();
    initGeography();
    initCoriolis();
    initConditions();
}
void TestSolver::setSavePath(const std::string& tPath) {
    mSavePath = tPath;
}
void TestSolver::openFiles() {
    //---Saving---//

    mAmpFile.open(mSavePath + "Amplitude.dat");
//    mVelFileX.open(mSavePath + "xVelocity.dat");
//    mVelFileY.open(mSavePath + "yVelocity.dat");

    //---Saving---//
}
void TestSolver::solveCustom() {
//    for (int iTime = 0; iTime < 128 * 24 * 60; iTime++) {
    for (int iTime = 0; iTime < 32 * 24 * 60; iTime++) {
        double FullEnergy = getFullEnergy();

//        if (iTime % 960 == 0) {
        if (iTime % 480 == 0) {
            appendData();
            std::cout   << "Step: " << iTime
                                    << " Full energy: " << FullEnergy
                                    << " (limit: " << mInitEnergy * 1.1 << ")"
                                    << std::endl;
        }

        if (FullEnergy > mInitEnergy * 1.1) {
            break;
        }

        //----------//

        for (int i = 1; i < mGridX - 1; i++) {
            for (int j = 1; j < mGridY - 1; j++) {
//                auto Extra = mStepTime * (source(i - 1, j - 1) + viscosity(i, j) * 1000.0);
                
//                dVector <double, 5> SmoothX(0.0, 0.0, 0.0, 0.0, 0.0);
//                dVector <double, 5> SmoothY(0.0, 0.0, 0.0, 0.0, 0.0);
//
//                for (int k = 1; k < 5; k++) {
//                    SmoothX[k] =
//                            fabs((*CurrentData)[i + 1][j][k] / (*CurrentData)[i + 1][j][0] - (*CurrentData)[i][j][k] / (*CurrentData)[i][j][0]) *
//                            ((*CurrentData)[i + 1][j][k] / (*CurrentData)[i + 1][j][0] - (*CurrentData)[i][j][k] / (*CurrentData)[i][j][0]) -
//                            fabs((*CurrentData)[i][j][k] / (*CurrentData)[i][j][0] - (*CurrentData)[i - 1][j][k] / (*CurrentData)[i - 1][j][0]) *
//                            ((*CurrentData)[i][j][k] / (*CurrentData)[i][j][0] - (*CurrentData)[i - 1][j][k] / (*CurrentData)[i - 1][j][0]);
//                    SmoothY[k] =
//                            fabs((*CurrentData)[i][j + 1][k] / (*CurrentData)[i][j + 1][0] - (*CurrentData)[i][j][k] / (*CurrentData)[i][j][0]) *
//                            ((*CurrentData)[i][j + 1][k] / (*CurrentData)[i][j + 1][0] - (*CurrentData)[i][j][k] / (*CurrentData)[i][j][0]) -
//                            fabs((*CurrentData)[i][j][k] / (*CurrentData)[i][j][0] - (*CurrentData)[i][j - 1][k] / (*CurrentData)[i][j - 1][0]) *
//                            ((*CurrentData)[i][j][k] / (*CurrentData)[i][j][0] - (*CurrentData)[i][j - 1][k] / (*CurrentData)[i][j - 1][0]);
//                }

                auto Extra =
                        mStepTime * (source(i - 1, j - 1));

//                auto Extra =
//                        mStepTime * (source(i - 1, j - 1)) +
//                        mStepTime * viscosity(i, j) +
//                        10.0 * mStepTime * (SmoothX / mStepX + SmoothY / mStepY) * (*CurrentData)[i][j][0];
//                auto Extra =
//                        mStepTime * (source(i - 1, j - 1)) +
//                        artVisc(i, j, 0.0, 0.00005);

                (*TempData)[i][j] = solveZwas(
                        (*CurrentData)[i][j],
                        (*CurrentData)[i - 1][j],
                        (*CurrentData)[i + 1][j],
                        (*CurrentData)[i][j - 1],
                        (*CurrentData)[i][j + 1],
                        (*CurrentData)[i + 1][j + 1],
                        (*CurrentData)[i - 1][j - 1],
                        (*CurrentData)[i + 1][j - 1],
                        (*CurrentData)[i - 1][j + 1],
                        Extra);

//                (*TempData)[i][j] = solve(
//                        (*CurrentData)[i - 1][j],
//                        (*CurrentData)[i + 1][j],
//                        (*CurrentData)[i][j - 1],
//                        (*CurrentData)[i][j + 1]) + Extra;
            }
        }

        //----------//

        for (int i = 1; i < mGridY - 1; i++) {
            (*TempData)[0][i][0] = (*TempData)[mGridX - 2][i][0];
            (*TempData)[mGridX - 1][i][0] = (*TempData)[1][i][0];
        }

        for (int i = 1; i < mGridY - 1; i++) {
            (*TempData)[0][i][1] = (*TempData)[mGridX - 2][i][1] / (*TempData)[mGridX - 2][i][0] * (*TempData)[0][i][0];
            (*TempData)[mGridX - 1][i][1] = (*TempData)[1][i][1] / (*TempData)[1][i][0] * (*TempData)[mGridX - 1][i][0];

            (*TempData)[0][i][2] = (*TempData)[mGridX - 2][i][2] / (*TempData)[mGridX - 2][i][0] * (*TempData)[0][i][0];
            (*TempData)[mGridX - 1][i][2] = (*TempData)[1][i][2] / (*TempData)[1][i][0] * (*TempData)[mGridX - 1][i][0];

            (*TempData)[0][i][3] = (*TempData)[mGridX - 2][i][3] / (*TempData)[mGridX - 2][i][0] * (*TempData)[0][i][0];
            (*TempData)[mGridX - 1][i][3] = (*TempData)[3][i][3] / (*TempData)[3][i][0] * (*TempData)[mGridX - 1][i][0];

            (*TempData)[0][i][4] = (*TempData)[mGridX - 2][i][4] / (*TempData)[mGridX - 2][i][0] * (*TempData)[0][i][0];
            (*TempData)[mGridX - 1][i][4] = (*TempData)[1][i][4] / (*TempData)[1][i][0] * (*TempData)[mGridX - 1][i][0];
        }

        for (int i = 0; i < mGridX; i++) {
            (*TempData)[i][0][1] = (*TempData)[i][1][1] / (*TempData)[i][1][0] * (*TempData)[i][0][0];
            (*TempData)[i][mGridY - 1][1] = (*TempData)[i][mGridY - 2][1] / (*TempData)[i][mGridY - 2][0] * (*TempData)[i][mGridY - 1][0];

            (*TempData)[i][0][2] = (*TempData)[i][1][2] / (*TempData)[i][1][0] * (*TempData)[i][0][0];
            (*TempData)[i][mGridY - 1][2] = (*TempData)[i][mGridY - 2][2] / (*TempData)[i][mGridY - 2][0] * (*TempData)[i][mGridY - 1][0];

            (*TempData)[i][0][3] = (*TempData)[i][3][3] / (*TempData)[i][3][0] * (*TempData)[i][0][0];
            (*TempData)[i][mGridY - 1][3] = (*TempData)[i][mGridY - 2][3] / (*TempData)[i][mGridY - 2][0] * (*TempData)[i][mGridY - 1][0];

            (*TempData)[i][0][4] = (*TempData)[i][1][4] / (*TempData)[i][1][0] * (*TempData)[i][0][0];
            (*TempData)[i][mGridY - 1][4] = (*TempData)[i][mGridY - 2][4] / (*TempData)[i][mGridY - 2][0] * (*TempData)[i][mGridY - 1][0];
        }

        for (int i = 0; i < mGridX; i++) {
            (*TempData)[i][0][2] = 0.0;
            (*TempData)[i][mGridY - 1][2] = 0.0;

            (*TempData)[i][0][4] = 0.0;
            (*TempData)[i][mGridY - 1][4] = 0.0;
        }

        //----------//

        std::swap(CurrentData, TempData);
    }

    saveData();
}

void TestSolver::initGrid() {
    mDataFirst.resize(mGridX);
    mDataSecond.resize(mGridX);

    mGeography.resize(mGridX);

    for (unsigned long i = 0; i < mGridX; i++) {
        mDataFirst[i].resize(mGridY);
        mDataSecond[i].resize(mGridY);

        mGeography[i].resize(mGridY);
    }
}
void TestSolver::initGeography() {
    double StdX = 5.0 * mStepX;
    double StdY = 5.0 * mStepY;

    double MeanX = int(mGridX / 2) * mStepX;
    double MeanY = int(mGridY / 2) * mStepY;

    //---Random---//
//    std::random_device rd;
//    std::mt19937 gen(rd());
//    std::uniform_real_distribution <> dis(0.0, 500.0);
    //---Random---//

    //---File---//
//    std::ifstream File("geography.dat", std::ios::in);
    //---File---//

    for (int i = 0; i < mGridX; i++) {
        for (int j = 0; j < mGridY; j++) {
            mGeography[i][j] = 4000 * exp(
                    -0.5 * pow((i * mStepX - MeanX) / StdX, 2.0)
                    -0.5 * pow((j * mStepY - MeanY) / StdY, 2.0));
//
//            mGeography[i][j] = dis(gen);
//
//            mGeography[i][j] = 0.0;

//            File >> mGeography[i][j];
        }
    }

//    File.close();
}
void TestSolver::initCoriolis() {
    double MeanY = int(mGridY / 2) * mStepY;

    mCorParam.resize(mGridY);

    for (int i = 0; i < mGridY; i++) {
        mCorParam[i] = mCorParam_0 + mBetaParam * (i * mStepY - MeanY);
    }
}
void TestSolver::initConditions() {
    double MeanWind = 20.0;
    double MeanY = int(mGridY / 2) * mStepY;

    for (int i = 0; i < mGridX; i++) {
        for (int j = 0; j < mGridY; j++) {
            mDataFirst[i][j] = dVector<double, 5>(10000.0 - (MeanWind * mCorParam_0 / mGrav) * (j * mStepY - MeanY),
                                                  0.0, 0.0);
            mDataSecond[i][j] = dVector<double, 5>(10000.0 - (MeanWind * mCorParam_0 / mGrav) * (j * mStepY - MeanY),
                                                   0.0, 0.0);
        }
    }

    double StdX = 5.0 * mStepX;
    double StdY = 5.0 * mStepY;

    double MeanX = int(mGridX / 2) * mStepX;

    for (int i = 0; i < mGridX; i++) {
        for (int j = 0; j < mGridY; j++) {
            mDataFirst[i][j][0] -= 4000 * exp(
                    -0.5 * pow((i * mStepX - MeanX) / StdX, 2.0)
                    -0.5 * pow((j * mStepY - MeanY) / StdY, 2.0));

            mDataSecond[i][j][0] -= 4000 * exp(
                    -0.5 * pow((i * mStepX - MeanX) / StdX, 2.0)
                    -0.5 * pow((j * mStepY - MeanY) / StdY, 2.0));

////            mGeography[i][j] = dis(gen);
//            mDataFirst[i][j][0] -= mGeography[i][j];
//            mDataSecond[i][j][0] -= mGeography[i][j];
        }
    }

    for (int i = 0; i < mGridX; i++) {
        for (int j = 1; j < mGridY - 1; j++) {
            mDataFirst[i][j][1] = mDataFirst[i][j][0] *
                    (-0.5 * mGrav / (mCorParam[j] * mStepX) * (mDataFirst[i][j + 1][0] - mDataFirst[i][j - 1][0]));
            mDataSecond[i][j][1] = mDataSecond[i][j][0] *
                    (-0.5 * mGrav / (mCorParam[j] * mStepX) * (mDataFirst[i][j + 1][0] - mDataFirst[i][j - 1][0]));
        }
    }

    for (int i = 1; i < mGridX - 1; i++) {
        for (int j = 0; j < mGridY; j++) {
            mDataFirst[i][j][2] = mDataFirst[i][j][0] *
                    (0.5 * mGrav / (mCorParam[j] * mStepX) * (mDataFirst[i + 1][j][0] - mDataFirst[i - 1][j][0]));
            mDataSecond[i][j][2] = mDataSecond[i][j][0] *
                    (0.5 * mGrav / (mCorParam[j] * mStepX) * (mDataFirst[i + 1][j][0] - mDataFirst[i - 1][j][0]));
        }
    }

    for (int i = 0; i < mGridY; i++) {
        mDataFirst[0][i][1] = mDataFirst[1][i][1] / mDataFirst[1][i][0] * mDataFirst[0][i][0];
        mDataFirst[mGridX - 1][i][1] = mDataFirst[mGridX - 2][i][1] / mDataFirst[mGridX - 2][i][0] * mDataFirst[mGridX - 1][i][0];

        mDataSecond[0][i][1] = mDataSecond[1][i][1] / mDataSecond[1][i][0] * mDataSecond[0][i][0];
        mDataSecond[mGridX - 1][i][1] = mDataSecond[mGridX - 2][i][1] / mDataSecond[mGridX - 2][i][0] * mDataSecond[mGridX - 1][i][0];
    }

    for (int i = 0; i < mGridX; i++) {
        mDataFirst[i][0][2] = 0;
        mDataFirst[i][mGridY - 1][2] = 0;

        mDataSecond[i][0][2] = 0;
        mDataSecond[i][mGridY - 1][2] = 0;
    }

    //------//

    mInitEnergy = getFullEnergy();
}

double TestSolver::getFullEnergy() {
    double Energy = 0.0;

    for (const auto& iLine : (*CurrentData)) {
        for (const auto& iVal : iLine) {
            Energy += (
                    mGrav * iVal[0] +
                    pow(iVal[1] / iVal[0], 2.0) +
                    pow(iVal[2] / iVal[0], 2.0) +
                    pow(iVal[3] / iVal[0], 2.0) +
                    pow(iVal[4] / iVal[0], 2.0));
        }
    }

    return Energy;
}
void TestSolver::appendData() {
    for (int j = 0; j < mGridY; j++) {
        for (int i = 0; i < mGridX; i++) {
//            mAmpFile << (*CurrentData)[i][j][0] + mGeography[i][j] << "\t";
            mAmpFile << (*CurrentData)[i][j][0] << "\t";
//            mVelFileX << (*CurrentData)[i][j][1] / (*CurrentData)[i][j][0] << "\t";
//            mVelFileY << (*CurrentData)[i][j][2] / (*CurrentData)[i][j][0] << "\t";
        }

        mAmpFile << std::endl;
//        mVelFileX << std::endl;
//        mVelFileY << std::endl;
    }
}
void TestSolver::saveData() {
    mAmpFile.close();
//    mVelFileX.close();
//    mVelFileY.close();

    system("python3.6 Plotting.py 254 50 0 254 0 50");
}

dVector <double, 5> TestSolver::funcX(const dVector <double, 5>& tVec) {
    return dVector <double, 5> (
            tVec[1],
            (pow(tVec[1], 2.0) - pow(tVec[3], 2.0)) / tVec[0] + 0.5 * mGrav * pow(tVec[0], 2.0),
            (tVec[1] * tVec[2] - tVec[3] * tVec[4]) / tVec[0],
            0.0,
            (tVec[1] * tVec[4] - tVec[2] * tVec[3]) / tVec[0]);
}
dVector <double, 5> TestSolver::funcY(const dVector <double, 5>& tVec) {
    return dVector <double, 5> (
            tVec[2],
            (tVec[1] * tVec[2] - tVec[3] * tVec[4]) / tVec[0],
            (pow(tVec[2], 2.0) - pow(tVec[4], 2.0)) / tVec[0] + 0.5 * mGrav * pow(tVec[0], 2.0),
            (tVec[2] * tVec[3] - tVec[1] * tVec[4]) / tVec[0],
            0.0);
}
dVector <double, 5> TestSolver::source(int tPosX, int tPosY) {
    return dVector <double, 5> (
            0.0,
            -mField_0 * (*CurrentData)[tPosX + 1][tPosY + 1][3] / (*CurrentData)[tPosX + 1][tPosY + 1][0] +
            mCorParam[tPosY + 1] * (*CurrentData)[tPosX + 1][tPosY + 1][2] -
            mGrav / (2.0 * mStepX) * (mGeography[tPosX + 2][tPosY + 1] -
            mGeography[tPosX][tPosY + 1]) * (*CurrentData)[tPosX + 1][tPosY + 1][0],
            -mField_0 * (*CurrentData)[tPosX + 1][tPosY + 1][4] / (*CurrentData)[tPosX + 1][tPosY + 1][0] -
            mCorParam[tPosY + 1] * (*CurrentData)[tPosX + 1][tPosY + 1][1] -
            mGrav / (2.0 * mStepY) * (mGeography[tPosX + 1][tPosY + 2] -
            mGeography[tPosX + 1][tPosY]) * (*CurrentData)[tPosX + 1][tPosY + 1][0],
            mField_0 * (*CurrentData)[tPosX + 1][tPosY + 1][1] / (*CurrentData)[tPosX + 1][tPosY + 1][0],
            mField_0 * (*CurrentData)[tPosX + 1][tPosY + 1][2] / (*CurrentData)[tPosX + 1][tPosY + 1][0]);
}
dVector <double, 5> TestSolver::viscosity(int tPosX, int tPosY) {
    double v_x_xx = (
            (*CurrentData)[tPosX - 1][tPosY][1] / (*CurrentData)[tPosX - 1][tPosY][0] +
            (*CurrentData)[tPosX][tPosY][1] / (*CurrentData)[tPosX][tPosY][0] * 2.0 +
            (*CurrentData)[tPosX + 1][tPosY][1] / (*CurrentData)[tPosX + 1][tPosY][0]) /
            pow(mStepX, 2.0);
    double v_x_yy = (
            (*CurrentData)[tPosX][tPosY - 1][1] / (*CurrentData)[tPosX][tPosY - 1][0] +
            (*CurrentData)[tPosX][tPosY][1] / (*CurrentData)[tPosX][tPosY][0] * 2.0 +
            (*CurrentData)[tPosX][tPosY + 1][1] / (*CurrentData)[tPosX][tPosY + 1][0]) /
            pow(mStepY, 2.0);
    double v_y_xx = (
            (*CurrentData)[tPosX - 1][tPosY][2] / (*CurrentData)[tPosX - 1][tPosY][0] +
            (*CurrentData)[tPosX][tPosY][2] / (*CurrentData)[tPosX][tPosY][0] * 2.0 +
            (*CurrentData)[tPosX + 1][tPosY][2] / (*CurrentData)[tPosX + 1][tPosY][0]) /
            pow(mStepX, 2.0);
    double v_y_yy = (
            (*CurrentData)[tPosX][tPosY - 1][2] / (*CurrentData)[tPosX][tPosY - 1][0] +
            (*CurrentData)[tPosX][tPosY][2] / (*CurrentData)[tPosX][tPosY][0] * 2.0 +
            (*CurrentData)[tPosX][tPosY + 1][2] / (*CurrentData)[tPosX][tPosY + 1][0]) /
            pow(mStepY, 2.0);
//    double hB_x_xx = (
//            (*CurrentData)[tPosX - 1][tPosY][3] +
//            (*CurrentData)[tPosX][tPosY][3] * 2.0 +
//            (*CurrentData)[tPosX + 1][tPosY][3]) /
//            pow(mStepX, 2.0);
//    double hB_x_yy = (
//            (*CurrentData)[tPosX][tPosY - 1][3] +
//            (*CurrentData)[tPosX][tPosY][3] * 2.0 +
//            (*CurrentData)[tPosX][tPosY + 1][3]) /
//            pow(mStepY, 2.0);
//    double hB_y_xx = (
//            (*CurrentData)[tPosX - 1][tPosY][4] +
//            (*CurrentData)[tPosX][tPosY][4] * 2.0 +
//            (*CurrentData)[tPosX + 1][tPosY][4]) /
//            pow(mStepX, 2.0);
//    double hB_y_yy = (
//            (*CurrentData)[tPosX][tPosY - 1][4] +
//            (*CurrentData)[tPosX][tPosY][4] * 2.0 +
//            (*CurrentData)[tPosX][tPosY + 1][4]) /
//            pow(mStepY, 2.0);

    return dVector <double, 5> (
            0.0,
            (*CurrentData)[tPosX][tPosY][0] * (v_x_xx + v_x_yy),
            (*CurrentData)[tPosX][tPosY][0] * (v_y_xx + v_y_yy),
            0.0,
            0.0);
//                hB_x_xx + hB_x_yy,
//                hB_y_xx + hB_y_yy);
}
dVector <double, 5> TestSolver::artVisc(int tPosX, int tPosY, double tFirstParam, double tSecondParam) {
    auto Visc = dVector <double, 5>(0.0, 0.0, 0.0, 0.0, 0.0);

    double dhVx = ((*CurrentData)[tPosX + 1][tPosY][1] - (*CurrentData)[tPosX - 1][tPosY][1]) / (2.0 * mStepX);
    double dhVy = ((*CurrentData)[tPosX][tPosY + 1][2] - (*CurrentData)[tPosX][tPosY - 1][2]) / (2.0 * mStepY);

    if (dhVx > 0.0) {
        Visc[1] = tFirstParam * dhVx * 200.0 + tSecondParam * pow(dhVx, 2.0);
    }

    if (dhVy > 0.0) {
        Visc[2] = tFirstParam * dhVy * 200.0 + tSecondParam * pow(dhVy, 2.0);
    }

    return Visc;
}