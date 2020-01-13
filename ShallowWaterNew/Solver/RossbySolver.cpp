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

    //---Geography---//

//    mGeography.resize(mGridX);
//
//    for (unsigned long i = 0; i < mGridX; i++) {
//        mGeography[i].resize(mGridY);
//    }
//
//    for (int i = 0; i < mGridX; i++) {
//        for (int j = 0; j < mGridY; j++) {
//            mGeography[i][j][0] = 1000.0 * exp(-0.5 * (pow((i - mGridX / 2.0) / 5.0, 2.0) + pow((j - mGridY / 2.0) / 5.0, 2.0)));
//        }
//    }

//    mGradient = gradient(mGeography);

    //---Initializing depth gradient from 5700m to 5500m---//

    double LatitudeGrad = (5700.0 - 5500.0) / mGridY;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution <> dis(0.0, 1.0);

    for (unsigned long i = 0; i < mGridX; i++) {
        for (unsigned long j = 0; j < mGridY; j++) {
//            DataFirst[i][j] = dVector <double, 3>(5700 - LatitudeGrad * double(j), 0.0, 0.0);
//            DataSecond[i][j] = dVector <double, 3>(5700 - LatitudeGrad * double(j), 0.0, 0.0);
            DataFirst[i][j] = dVector <double, 3>(5600 - 2 * mCorParam_0 / mGrav * (double(j) * mStepY - MeanCoordY), 0.0, 0.0);
            DataSecond[i][j] = dVector <double, 3>(5600 - 2 * mCorParam_0 / mGrav * (double(j) * mStepY - MeanCoordY), 0.0, 0.0);

            DataFirst[i][j][0] += dis(gen) * 5;
            DataSecond[i][j][0] += dis(gen) * 5;
        }
    }

//    for (unsigned long i = 0; i < mGridX; i++) {
//        for (unsigned long j = 0; j < mGridY; j++) {
//            DataFirst[i][j][0] -= mGeography[i][j][0];
//            DataSecond[i][j][0] -= mGeography[i][j][0];
//        }
//    }

    //---Initializing geostrophic wind---//

    auto Grad = gradient(DataFirst);

    for (unsigned long i = 0; i < mGridX; i++) {
        for (unsigned long j = 0; j < mGridY; j++) {
            DataFirst[i][j][1] = -mGrav / mCorParam[j] * Grad[i][j][1] / mStepX * DataFirst[i][j][0];
            DataSecond[i][j][2] = mGrav / mCorParam[j] * Grad[i][j][0] / mStepY * DataSecond[i][j][0];
        }
    }

    //---Boundaries---//

    mLeftBoundary.resize(mGridY);
    mRightBoundary.resize(mGridY);
    mUpperBoundary.resize(mGridX);
    mLowerBoundary.resize(mGridX);

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
    dVector <double, 3> Extra(0, 0, 0);
    double Timer = 0.0;

    for (int iTime = 0; iTime < 20000; iTime++) {
        double FullEnergy = getFullEnergy();

        if (iTime % 200 == 0) {
            appendData();
            std::cout   << "Step: " << iTime
                        << " Full energy: " << FullEnergy
                        << " (limit: " << mInitEnergy * 1.1
                        << ") Time: " << Timer << "s"
                        << std::endl;
        }

        if (FullEnergy > mInitEnergy * 1.1) {
            saveData();
            break;
        }

//        double Time = mStepY / getFullEnergy() * 1000.0;
        double Time = 14.0;

        Timer += Time;

        setTimeStep(Time);
        boundaries();

        for (int i = 1; i < mGridX - 1; i++) {
            for (int j = 1; j < mGridY - 1; j++) {
                Extra = source(i, j);

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
            }
        }

        for (int i = 1; i < mGridY - 1; i++) {
            Extra = source(0, i);

            (*TempData)[0][i] = solveZwas(
                    (*CurrentData)[0][i],
                    mLeftBoundary[i],
                    (*CurrentData)[1][i],
                    (*CurrentData)[0][i - 1],
                    (*CurrentData)[0][i + 1],
                    (*CurrentData)[1][i + 1],
                    mLeftBoundary[i - 1],
                    (*CurrentData)[1][i - 1],
                    mLeftBoundary[i + 1],
                    Extra);

            Extra = source(mGridX - 1, i);

            (*TempData)[mGridX - 1][i] = solveZwas(
                    (*CurrentData)[mGridX - 1][i],
                    (*CurrentData)[mGridX - 2][i],
                    mRightBoundary[i],
                    (*CurrentData)[mGridX - 1][i - 1],
                    (*CurrentData)[mGridX - 1][i + 1],
                    mRightBoundary[i + 1],
                    (*CurrentData)[mGridX - 2][i - 1],
                    mRightBoundary[i - 1],
                    (*CurrentData)[mGridX - 2][i + 1],
                    Extra);
        }

        for (int i = 1; i < mGridX - 1; i++) {
            Extra = source(i, 0);

            (*TempData)[i][0] = solveZwas(
                    (*CurrentData)[i][0],
                    (*CurrentData)[i - 1][0],
                    (*CurrentData)[i + 1][0],
                    mUpperBoundary[i],
                    (*CurrentData)[i][1],
                    (*CurrentData)[i + 1][1],
                    mUpperBoundary[i - 1],
                    mUpperBoundary[i + 1],
                    (*CurrentData)[i - 1][1],
                    Extra);

            Extra = source(i, mGridY - 1);

            (*TempData)[i][mGridY - 1] = solveZwas(
                    (*CurrentData)[i][mGridY - 1],
                    (*CurrentData)[i - 1][mGridY - 1],
                    (*CurrentData)[i + 1][mGridY - 1],
                    (*CurrentData)[i][mGridY - 2],
                    mLowerBoundary[i],
                    mLowerBoundary[i + 1],
                    (*CurrentData)[i - 1][mGridY - 2],
                    (*CurrentData)[i + 1][mGridY - 2],
                    mLowerBoundary[i - 1],
                    Extra);
        }

        //---Left-Top---//
        Extra = source(0, 0);

        (*TempData)[0][0] = solveZwas(
                (*CurrentData)[0][0],
                mLeftBoundary[0],
                (*CurrentData)[1][0],
                mUpperBoundary[0],
                (*CurrentData)[0][1],
                (*CurrentData)[1][1],
                dVector <double, 3>(5700, mUpperBoundary[1][1], 0),
                mUpperBoundary[1],
                mLeftBoundary[1],
                Extra);
        //---Left-Top---//

        //---Right-Top---//
        Extra = source(mGridX - 1, 0);

        (*TempData)[mGridX - 1][0] = solveZwas(
                (*CurrentData)[mGridX - 1][0],
                (*CurrentData)[mGridX - 2][0],
                mRightBoundary[0],
                mUpperBoundary[mGridX - 1],
                (*CurrentData)[mGridX - 1][1],
                mRightBoundary[1],
                mUpperBoundary[mGridX - 2],
                dVector <double, 3>(5700, mUpperBoundary[1][1], 0),
                (*CurrentData)[mGridX - 2][1],
                Extra);
        //---Right-Top---//

        //---Left-Bottom---//
        Extra = source(0, mGridY - 1);

        (*TempData)[0][mGridY - 1] = solveZwas(
                (*CurrentData)[0][mGridY - 1],
                mLeftBoundary[mGridY - 1],
                (*CurrentData)[1][mGridY - 1],
                (*CurrentData)[0][mGridY - 2],
                mLowerBoundary[0],
                mLowerBoundary[1],
                mLeftBoundary[mGridY - 2],
                (*CurrentData)[1][mGridY - 2],
                dVector <double, 3>(5500, mLowerBoundary[1][1], 0),
                Extra);
        //---Left-Bottom---//

        //---Right-Bottom---//
        Extra = source(mGridX - 1, mGridY - 1);

        (*TempData)[mGridX - 1][mGridY - 1] = solveZwas(
                (*CurrentData)[mGridX - 1][mGridY - 1],
                (*CurrentData)[mGridX - 2][mGridY - 1],
                mRightBoundary[mGridY - 1],
                (*CurrentData)[mGridX - 1][mGridY - 2],
                mLowerBoundary[mGridX - 1],
                dVector <double, 3>(5500, mLowerBoundary[1][1], 0),
                (*CurrentData)[mGridX - 2][mGridY - 2],
                mRightBoundary[mGridY - 2],
                mLowerBoundary[mGridX - 2],
                Extra);
        //---Right-Bottom---//

        std::swap(CurrentData, TempData);
    }

    saveData();
}
//-----------------------------//
std::vector <std::vector <std::vector <double>>> RossbySolver::gradient(const dGrid& tGrid) {
    std::vector <std::vector <std::vector <double>>> Grad; //---[0] - x, [1] - y---//

    Grad.resize(mGridX);

    for (unsigned long i = 0; i < mGridX; i++) {
        Grad[i].resize(mGridY);
    }

    //----------//

    for (unsigned long j = 0; j < mGridY; j++) {
        for (unsigned long i = 1; i < mGridX - 1; i++) {
            Grad[i][j].emplace_back((tGrid[i + 1][j][0] - tGrid[i - 1][j][0]) / 2);
        }

        Grad[0][j].emplace_back(Grad[1][j][0] * 2.0 - Grad[2][j][0]);
        Grad[mGridX - 1][j].emplace_back(Grad[mGridX - 2][j][0] * 2.0 - Grad[mGridX - 3][j][0]);
    }

    for (unsigned long i = 0; i < mGridX; i++) {
        for (unsigned long j = 1; j < mGridY - 1; j++) {
            Grad[i][j].emplace_back((tGrid[i][j + 1][0] - tGrid[i][j - 1][0]) / 2);
        }

        Grad[i][0].emplace_back(Grad[i][1][1] * 2.0 - Grad[i][2][1]);
        Grad[i][mGridY - 1].emplace_back(Grad[i][mGridY - 2][1] * 2.0 - Grad[i][mGridY - 3][1]);
    }

    return Grad;
}
void RossbySolver::boundaries() {
    for (int i = 0; i < mGridY; i++) {
        mLeftBoundary[i] = (*CurrentData)[1][i];
    }

    for (int i = 0; i < mGridX; i++) {
        mUpperBoundary[i][0] = 5700;
        mLowerBoundary[i][0] = 5500;

        mUpperBoundary[i][1] = (*CurrentData)[i][1][1];
        mLowerBoundary[i][1] = (*CurrentData)[i][mGridY - 2][1];
    }

    mRightBoundary = mLeftBoundary;

    for (int i = 0; i < mGridX; i++) {
        mUpperBoundary[i][2] = 0;
        mLowerBoundary[i][2] = 0;
    }
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
            mVelFileX << (*CurrentData)[i][j][1] / (*CurrentData)[i][j][0] << "\t";
            mVelFileY << (*CurrentData)[i][j][2] / (*CurrentData)[i][j][0] << "\t";
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
//            -(*CurrentData)[tPosX][tPosY][2] * mCorParam[tPosY] + mGrav * mGradient[tPosX][tPosY][0] / mStepX,
//            (*CurrentData)[tPosX][tPosY][1] * mCorParam[tPosY]) + mGrav * mGradient[tPosX][tPosY][1] / mStepY;
            -(*CurrentData)[tPosX][tPosY][2] * mCorParam[tPosY],
            (*CurrentData)[tPosX][tPosY][1] * mCorParam[tPosY]);
}