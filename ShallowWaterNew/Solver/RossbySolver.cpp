#include "RossbySolver.h"
//-----------------------------//
RossbySolver::RossbySolver(double tStepTime, double tStepX, double tStepY) {
    mStepTime = tStepTime;

    mStepX = tStepX;
    mStepY = tStepY;
}
//-----------------------------//
void RossbySolver::initGrid(int tGridX, int tGridY) {
    mGridX = tGridX;
    mGridY = tGridY;

    DataFirst.resize(mGridX);
    DataSecond.resize(mGridX);

    for (unsigned long i = 0; i < mGridX; i++) {
        DataFirst[i].resize(mGridY);
        DataSecond[i].resize(mGridY);
    }

    //---Prepare angles and parameters---//

    for (int i = 0; i < mGridY; i++) {
        mAngle.emplace_back(i * 88.0 / mGridY * M_PI / 180.0);
    }

    mAngle_0 = mAngle[int(mGridY / 2.0)];

    mBetaParam = mOmega * 2.0 * cos(mAngle_0) / mRad;
    mCorParam_0 = mOmega * 2.0 * sin(mAngle_0);

    double MeanCoordY = int(mGridY / 2.0) * mStepY;

    for (int i = 0; i < mGridY; i++) {
        mCorParam.emplace_back(mCorParam_0 + mBetaParam * (i * mStepY - MeanCoordY));
    }

    //---Initializing depth gradient from 5700m to 5500m---//

    double LatitudeGrad = (5700.0 - 5500.0) / mGridY;

    for (unsigned long i = 0; i < mGridX; i++) {
        for (unsigned long j = 0; j < mGridY; j++) {
            DataFirst[i][j] = dVector <double, 3>(5500 + LatitudeGrad * double(j), 0.0, 0.0);
            DataSecond[i][j] = dVector <double, 3>(5500 + LatitudeGrad * double(j), 0.0, 0.0);
        }
    }

    //---Initializing geostrophic wind---//

    auto Grad = gradient();

    for (unsigned long i = 0; i < mGridX; i++) {
        for (unsigned long j = 0; j < mGridY; j++) {
            DataFirst[i][j][1] = -mGrav / mCorParam[j] * Grad[i][j][1] / mStepX * DataFirst[i][j][0];
            DataSecond[i][j][2] = mGrav / mCorParam[j] * Grad[i][j][0] / mStepY * DataSecond[i][j][0];
        }
    }

    //---Energy---//

    mInitEnergy = getFullEnergy();

    //---Saving---//

    mAmpFile.open(mSavePath + "Amplitude.dat");
    mVelFileX.open(mSavePath + "xVelocity.dat");
    mVelFileY.open(mSavePath + "yVelocity.dat");
}
void RossbySolver::setSavePath(const std::string& tPath) {
    mSavePath = tPath;
}
void RossbySolver::solve() {
    long xIndex_plus_1;
    long xIndex_minus_1;
    long yIndex_plus_1;
    long yIndex_minus_1;

    dVector <double, 3> Extra(0, 0, 0);

    for (int iTime = 0; iTime < 10000; iTime++) {
        double FullEnergy = getFullEnergy();

        if (iTime % 10 == 0) {
            appendData();
            std::cout << "Step: " << iTime << " Full energy: " << FullEnergy << std::endl;
        }

        if (FullEnergy > mInitEnergy * 3.0) {
            saveData();
            break;
        }

        setTimeStep(mStepY / getFullEnergy() * 5000.0);

        for (int i = 0; i < mGridX; i++) {
            for (int j = 0; j < mGridY; j++) {
                xIndex_plus_1 = (i + 1 == mGridX ? 0 : i + 1);
                xIndex_minus_1 = (i - 1 < 0 ? mGridX - 1 : i - 1);
                yIndex_plus_1 = (j + 1 == mGridY ? 0 : j + 1);
                yIndex_minus_1 = (j - 1 < 0 ? mGridY - 1 : j - 1);

                Extra = source(i, j);

                (*TempData)[i][j] = solveZwas(
                        (*CurrentData)[i][j],
                        (*CurrentData)[xIndex_minus_1][j],
                        (*CurrentData)[xIndex_plus_1][j],
                        (*CurrentData)[i][yIndex_minus_1],
                        (*CurrentData)[i][yIndex_plus_1],
                        (*CurrentData)[xIndex_plus_1][yIndex_plus_1],
                        (*CurrentData)[xIndex_minus_1][yIndex_minus_1],
                        (*CurrentData)[xIndex_plus_1][yIndex_minus_1],
                        (*CurrentData)[xIndex_minus_1][yIndex_plus_1],
                        Extra);
            }
        }

        std::swap(CurrentData, TempData);
    }

    saveData();
}
//-----------------------------//
std::vector <std::vector <std::vector <double>>> RossbySolver::gradient() {
    std::vector <std::vector <std::vector <double>>> Grad; //---[0] - x, [1] - y---//

    Grad.resize(mGridX);

    for (unsigned long i = 0; i < mGridX; i++) {
        Grad[i].resize(mGridY);
    }

    //----------//

    for (unsigned long j = 0; j < mGridY; j++) {
        for (unsigned long i = 1; i < mGridX - 1; i++) {
            Grad[i][j].emplace_back((DataFirst[i + 1][j][0] - DataFirst[i - 1][j][0]) / 2);
        }

        Grad[0][j].emplace_back(Grad[1][j][0] * 2.0 - Grad[2][j][0]);
        Grad[mGridX - 1][j].emplace_back(Grad[mGridX - 2][j][0] * 2.0 - Grad[mGridX - 3][j][0]);
    }

    for (unsigned long i = 0; i < mGridX; i++) {
        for (unsigned long j = 1; j < mGridY - 1; j++) {
            Grad[i][j].emplace_back((DataFirst[i][j + 1][0] - DataFirst[i][j - 1][0]) / 2);
        }

        Grad[i][0].emplace_back(Grad[i][1][1] * 2.0 - Grad[i][2][1]);
        Grad[i][mGridY - 1].emplace_back(Grad[i][mGridY - 2][1] * 2.0 - Grad[i][mGridY - 3][1]);
    }

    return Grad;
}
double RossbySolver::getMaxSpeed() {
    double VelMax = 0.0;

//    for (unsigned long i = 0; i < mGridX; i++) {
//        for (unsigned long j = 0; j < mGridY; j++) {
//            double CurVal = 0.0;
//
//            CurVal = sqrt((*CurrentData)[i][j][1] * (*CurrentData)[i][j][1] +
//                              (*CurrentData)[i][j][2] * (*CurrentData)[i][j][2]) / (*CurrentData)[i][j][0];
//
//            if (VelMax < CurVal) {
//                VelMax = CurVal;
//            }
//        }
//    }

    for (unsigned long i = 0; i < mGridX; i++) {
        for (unsigned long j = 0; j < mGridY; j++) {
            double Val =    fabs((*CurrentData)[i][j][2] / (*CurrentData)[i][j][1]) +
                            sqrt((*CurrentData)[i][j][1] * mGrav);

            if (Val > VelMax) {
                VelMax = Val;
            }
        }
    }

    return VelMax;
}
double RossbySolver::getFullEnergy() {
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
void RossbySolver::appendData() {
    for (int j = 0; j < mGridY; j++) {
        for (int i = 0; i < mGridX; i++) {
            mAmpFile << (*CurrentData)[i][j][0] << "\t";
            mVelFileX << (*CurrentData)[i][j][1] << "\t";
            mVelFileY << (*CurrentData)[i][j][2] << "\t";
        }

        mAmpFile << std::endl;
        mVelFileX << std::endl;
        mVelFileY << std::endl;
    }
}
void RossbySolver::saveData() {
    mAmpFile.close();
    mVelFileX.close();
    mVelFileY.close();

    system("python3.6 Plotting.py 180 88 0 180 0 88");
}

dVector <double, 3> RossbySolver::funcX(const dVector <double, 3>& tVec) {
    return dVector <double, 3> (
            tVec[1],
            pow(tVec[1], 2.0) / tVec[0] + 0.5 * mGrav * pow(tVec[0], 2.0),
            tVec[1] * tVec[2] / tVec[0]);
}
dVector <double, 3> RossbySolver::funcY(const dVector <double, 3>& tVec) {
    return dVector <double, 3> (
            tVec[2],
            tVec[1] * tVec[2] / tVec[0],
            pow(tVec[2], 2.0) / tVec[0] + 0.5 * mGrav * pow(tVec[0], 2.0));
}
dVector <double, 3> RossbySolver::source(int tPosX, int tPosY) {
    return dVector <double, 3> (
            0,
            -(*CurrentData)[tPosX][tPosY][2] * mCorParam[tPosY],
            (*CurrentData)[tPosX][tPosY][1] * mCorParam[tPosY]);
}