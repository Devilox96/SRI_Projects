#include "TestSolver.h"
//-----------------------------//
TestSolver::TestSolver() {
    mStepTime = 60.0;

    mStepX = 100.0e+03;
    mStepY = 100.0e+03;

    initGrid();
    initGeography();
    initCoriolis();
    initConditions();

    //---Solver---//

    uh.resize(mGridX);
    vh.resize(mGridX);

    for (unsigned long i = 0; i < mGridX; i++) {
        uh[i].resize(mGridY);
        vh[i].resize(mGridY);
    }

    Ux.resize(mGridX);
    Uy.resize(mGridX);

    for (unsigned long i = 0; i < mGridX; i++) {
        Ux[i].resize(mGridY);
        Uy[i].resize(mGridY);
    }

    Vx.resize(mGridX);
    Vy.resize(mGridX);

    for (unsigned long i = 0; i < mGridX; i++) {
        Vx[i].resize(mGridY);
        Vy[i].resize(mGridY);
    }

    //----------//

    mid_xt.resize(mGridX - 1);
    mid_yt.resize(mGridX);

    for (unsigned long i = 0; i < mGridX - 1; i++) {
        mid_xt[i].resize(mGridY);
    }

    for (unsigned long i = 0; i < mGridX; i++) {
        mid_yt[i].resize(mGridY - 1);
    }

    for (unsigned long i = 0; i < mGridX - 1; i++) {
        for (unsigned long j = 0; j < mGridY; j++) {
            mid_xt[i][j] = dVector <double, 3>(0.0, 0.0, 0.0);
        }
    }

    for (unsigned long i = 0; i < mGridX; i++) {
        for (unsigned long j = 0; j < mGridY - 1; j++) {
            mid_yt[i][j] = dVector<double, 3>(0.0, 0.0, 0.0);
        }
    }

    //-----------------------------//

    Ux_mid_xt.resize(mGridX - 1);
    Uy_mid_yt.resize(mGridX);

    for (unsigned long i = 0; i < mGridX - 1; i++) {
        Ux_mid_xt[i].resize(mGridY);
    }

    for (unsigned long i = 0; i < mGridX; i++) {
        Uy_mid_yt[i].resize(mGridY - 1);
    }

    //----------//

    Vx_mid_xt.resize(mGridX - 1);
    Vy_mid_yt.resize(mGridX);

    for (unsigned long i = 0; i < mGridX - 1; i++) {
        Vx_mid_xt[i].resize(mGridY);
    }

    for (unsigned long i = 0; i < mGridX; i++) {
        Vy_mid_yt[i].resize(mGridY - 1);
    }

    //---Solver---//
}
void TestSolver::setSavePath(const std::string& tPath) {
    mSavePath = tPath;
}
void TestSolver::openFiles() {
    //---Saving---//

    mAmpFile.open(mSavePath + "Amplitude.dat");
    mVelFileX.open(mSavePath + "xVelocity.dat");
    mVelFileY.open(mSavePath + "yVelocity.dat");

    //---Saving---//
}
void TestSolver::solve() {
    for (int iTime = 0; iTime < 4 * 24 * 60; iTime++) {
        if (iTime % 60 == 0) {
            appendData();
            std::cout   << "Step: " << iTime << std::endl;
        }

        //----------//

        for (int i = 0; i < mGridX; i++) {
            for (int j = 0; j < mGridY; j++) {
                uh[i][j] = (*CurrentData)[i][j][1] * (*CurrentData)[i][j][0];
                vh[i][j] = (*CurrentData)[i][j][2] * (*CurrentData)[i][j][0];

                Ux[i][j] = pow(uh[i][j], 2.0) / (*CurrentData)[i][j][0] + 0.5 * mGrav * pow((*CurrentData)[i][j][0], 2.0);
                Uy[i][j] = uh[i][j] * vh[i][j] / (*CurrentData)[i][j][0];

                Vx[i][j] = uh[i][j] * vh[i][j] / (*CurrentData)[i][j][0];
                Vy[i][j] = pow(vh[i][j], 2.0) / (*CurrentData)[i][j][0] + 0.5 * mGrav * pow((*CurrentData)[i][j][0], 2.0);
            }
        }

        //----------//

        for (unsigned long i = 0; i < mGridX - 1; i++) {
            for (unsigned long j = 0; j < mGridY; j++) {
                mid_xt[i][j] =
                        0.5 *
                        dVector <double, 3>(
                        (*CurrentData)[i + 1][j][0] + (*CurrentData)[i][j][0],
                        uh[i + 1][j] + uh[i][j],
                        vh[i + 1][j] + vh[i][j]) -
                        0.5 * (mStepTime / mStepX) *
                        dVector <double, 3>(
                        uh[i + 1][j] - uh[i][j],
                        Ux[i + 1][j] - Ux[i][j],
                        Vx[i + 1][j] - Vx[i][j]);

                Ux_mid_xt[i][j] = mid_xt[i][j][1] * mid_xt[i][j][1] / mid_xt[i][j][0] + 0.5 * mGrav * pow(mid_xt[i][j][0], 2.0);
                Vx_mid_xt[i][j] = mid_xt[i][j][1] * mid_xt[i][j][2] / mid_xt[i][j][0];
            }
        }

        for (unsigned long i = 0; i < mGridX; i++) {
            for (unsigned long j = 0; j < mGridY - 1; j++) {
                mid_yt[i][j] =
                        0.5 *
                        dVector <double, 3>(
                                (*CurrentData)[i][j + 1][0] + (*CurrentData)[i][j][0],
                                uh[i][j + 1] + uh[i][j],
                                vh[i][j + 1] + vh[i][j]) -
                        0.5 * (mStepTime / mStepY) *
                        dVector <double, 3>(
                                vh[i][j + 1] - vh[i][j],
                                Uy[i][j + 1] - Uy[i][j],
                                Vy[i][j + 1] - Vy[i][j]);

                Uy_mid_yt[i][j] = mid_yt[i][j][1] * mid_yt[i][j][2] / mid_yt[i][j][0];
                Vy_mid_yt[i][j] = mid_yt[i][j][2] * mid_yt[i][j][2] / mid_yt[i][j][0] + 0.5 * mGrav * pow(mid_yt[i][j][0], 2.0);
            }
        }

        //----------//

        for (int i = 1; i < mGridX - 1; i++) {
            for (int j = 1; j < mGridY - 1; j++) {
                auto Source = source(i - 1, j - 1);

                (*TempData)[i][j][0] = (*CurrentData)[i][j][0] -
                                       (mStepTime / mStepX) * (mid_xt[i][j][1] - mid_xt[i - 1][j][1]) -
                                       (mStepTime / mStepY) * (mid_yt[i][j][2] - mid_yt[i][j - 1][2]);
                (*TempData)[i][j][1] = (uh[i][j] -
                                       (mStepTime / mStepX) * (Ux_mid_xt[i][j] - Ux_mid_xt[i - 1][j]) -
                                       (mStepTime / mStepY) * (Uy_mid_yt[i][j] - Uy_mid_yt[i][j - 1]) +
                                       mStepTime * Source[1] * 0.5 * ((*CurrentData)[i][j][0] + (*TempData)[i][j][0])) /
                                       (*TempData)[i][j][0];
                (*TempData)[i][j][2] = (vh[i][j] -
                                        (mStepTime / mStepX) * (Vx_mid_xt[i][j] - Vx_mid_xt[i - 1][j]) -
                                        (mStepTime / mStepY) * (Vy_mid_yt[i][j] - Vy_mid_yt[i][j - 1]) +
                                        mStepTime * Source[2] * 0.5 * ((*CurrentData)[i][j][0] + (*TempData)[i][j][0])) /
                                       (*TempData)[i][j][0];
            }
        }

        for (int i = 1; i < mGridY - 1; i++) {
            (*TempData)[0][i][0] = (*TempData)[mGridX - 2][i][0];
            (*TempData)[mGridX - 1][i][0] = (*TempData)[1][i][0];
        }

        for (int i = 1; i < mGridY - 1; i++) {
            (*TempData)[0][i][1] = (*TempData)[mGridX - 2][i][1];
            (*TempData)[mGridX - 1][i][1] = (*TempData)[1][i][1];

            (*TempData)[0][i][2] = (*TempData)[mGridX - 2][i][2];
            (*TempData)[mGridX - 1][i][2] = (*TempData)[1][i][2];
        }

        for (int i = 0; i < mGridX; i++) {
            (*TempData)[i][0][1] = (*TempData)[i][1][1];
            (*TempData)[i][mGridY - 1][1] = (*TempData)[i][mGridY - 2][1];

            (*TempData)[i][0][2] = (*TempData)[i][1][2];
            (*TempData)[i][mGridY - 1][2] = (*TempData)[i][mGridY - 2][2];
        }

        for (int i = 0; i < mGridX; i++) {
            (*TempData)[i][0][2] = 0.0;
            (*TempData)[i][mGridY - 1][2] = 0.0;
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

    for (int i = 0; i < mGridX; i++) {
        for (int j = 0; j < mGridY; j++) {
            mGeography[i][j] = 4000 * exp(
                    -0.5 * pow((i * mStepX - MeanX) / StdX, 2.0)
                    -0.5 * pow((j * mStepY - MeanY) / StdY, 2.0));
        }
    }
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
            mDataFirst[i][j] = dVector<double, 3>(10000.0 - (MeanWind * mCorParam_0 / mGrav) * (j * mStepY - MeanY),
                                                  0.0, 0.0);
            mDataSecond[i][j] = dVector<double, 3>(10000.0 - (MeanWind * mCorParam_0 / mGrav) * (j * mStepY - MeanY),
                                                   0.0, 0.0);
        }
    }

    for (int i = 0; i < mGridX; i++) {
        for (int j = 1; j < mGridY - 1; j++) {
            mDataFirst[i][j][1] =
                    -0.5 * mGrav / (mCorParam[j] * mStepX) * (mDataFirst[i][j + 1][0] - mDataFirst[i][j - 1][0]);
            mDataSecond[i][j][1] =
                    -0.5 * mGrav / (mCorParam[j] * mStepX) * (mDataFirst[i][j + 1][0] - mDataFirst[i][j - 1][0]);
        }
    }

    for (int i = 1; i < mGridX - 1; i++) {
        for (int j = 0; j < mGridY; j++) {
            mDataFirst[i][j][2] =
                    0.5 * mGrav / (mCorParam[j] * mStepX) * (mDataFirst[i + 1][j][0] - mDataFirst[i - 1][j][0]);
            mDataSecond[i][j][2] =
                    0.5 * mGrav / (mCorParam[j] * mStepX) * (mDataFirst[i + 1][j][0] - mDataFirst[i - 1][j][0]);
        }
    }

    for (int i = 0; i < mGridY; i++) {
        mDataFirst[0][i][1] = mDataFirst[1][i][1];
        mDataFirst[mGridX - 1][i][1] = mDataFirst[mGridX - 2][i][1];

        mDataSecond[0][i][1] = mDataSecond[1][i][1];
        mDataSecond[mGridX - 1][i][1] = mDataSecond[mGridX - 2][i][1];
    }

    for (int i = 0; i < mGridX; i++) {
        mDataFirst[i][0][2] = 0;
        mDataFirst[i][mGridY - 1][2] = 0;

        mDataSecond[i][0][2] = 0;
        mDataSecond[i][mGridY - 1][2] = 0;
    }
}

void TestSolver::appendData() {
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
void TestSolver::saveData() {
    mAmpFile.close();
    mVelFileX.close();
    mVelFileY.close();

    system("python3.6 Plotting.py 254 50 0 254 0 50");
}

dVector <double, 3> TestSolver::source(int tPosX, int tPosY) {
    return dVector <double, 3> (
                0.0,
                mCorParam[tPosY + 1] * (*CurrentData)[tPosX + 1][tPosY + 1][2] - mGrav / (2.0 * mStepX) * (mGeography[tPosX + 2][tPosY + 1] - mGeography[tPosX][tPosY + 1]),
                -mCorParam[tPosY + 1] * (*CurrentData)[tPosX + 1][tPosY + 1][1] - mGrav / (2.0 * mStepY) * (mGeography[tPosX + 1][tPosY + 2] - mGeography[tPosX + 1][tPosY])
            );
}