#include "mainwindow.h"

void MainWindow :: Calc(unsigned int StepsP) {
    unsigned int PointNum = 100;

    InitGrid(Grid, 10.0, 0.0, 0.0, PointNum);

    Grid[50][50].x = 12.0;

    std :: vector <std :: vector <dVector3D <double>>> TempGrid = Grid;

    for (unsigned int iter = 0; iter < StepsP; iter++) {
        for (unsigned int i = 1; i < PointNum - 1; i++) {
            for (unsigned int j = 1; j < PointNum - 1; j++) {
                dVector3D <double> TempH = NextH(0.001, 0.1, 0.1, Grid, i, j, 9.81);

                TempGrid[i][j].x = TempH.x;
                TempGrid[i][j].y = TempH.y / TempH.x;
                TempGrid[i][j].z = TempH.z / TempH.x;
            }
        }

        for (unsigned int i = 0; i < PointNum; i++) {
            TempGrid[i][0] = TempGrid[i][1];
            TempGrid[i][PointNum - 1] = TempGrid[i][PointNum - 2];
            TempGrid[0][i] = TempGrid[1][i];
            TempGrid[PointNum - 1][i] = TempGrid[PointNum - 2][i];
        }

        Grid = TempGrid;
    }
}
//-----------------------------//
MainWindow::MainWindow(QWidget *parent) : QMainWindow(parent) {
    Calc(500);

    MainWidget = new QWidget;
    setCentralWidget(MainWidget);

    MainLayout = new QGridLayout;
    MainWidget -> setLayout(MainLayout);

    graph = new Q3DSurface;
    QWidget *container = QWidget :: createWindowContainer(graph);
    container -> setMinimumSize(800, 800);
    graph -> setFlags(graph -> flags() ^ Qt :: FramelessWindowHint);

    int sampleCountX = 100;
    int sampleCountZ = 100;

    auto dataArray = new QSurfaceDataArray;
    dataArray->reserve(sampleCountZ);

    for (int i = 0 ; i < sampleCountZ ; i++) {
        QSurfaceDataRow *newRow = new QSurfaceDataRow(sampleCountX);

        float z = i;
        int index = 0;
        for (int j = 0; j < sampleCountX; j++) {
            float x = j;

            auto y = float(Grid[i][j].x);
            (*newRow)[index++].setPosition(QVector3D(x, y, z));
        }
        *dataArray << newRow;
    }

    auto series = new QSurface3DSeries;
    series->setDrawMode(QSurface3DSeries::DrawSurface);
    series->dataProxy()->resetArray(dataArray);
    graph -> addSeries(series);

    MainLayout -> addWidget(container);
}
//-----------------------------//
void InitGrid(std :: vector <std :: vector <dVector3D <double>>>& GridP, double ExcitationP, double VXP, double VYP, int PointsNumP) {
    for (int i = 0; i < PointsNumP; i++) {
        GridP.emplace_back(std :: vector <dVector3D <double>>());

        for (int j = 0; j < PointsNumP; j++) {
            GridP.back().emplace_back(dVector3D <double>(ExcitationP, VXP, VYP));
        }
    }
}

dVector3D <double> GetH(const std :: vector <std :: vector <dVector3D <double>>>& GridP, unsigned int i, unsigned int j) {
    return dVector3D <double>(GridP[i][j].x,
                     GridP[i][j].y * GridP[i][j].x,
                     GridP[i][j].z * GridP[i][j].x);
}
dVector3D <double> GetU(const dVector3D <double>& HP, double gP) {
    return dVector3D <double>(HP.y,
                     pow(HP.y, 2.0) / HP.x + 0.5 * gP * pow(HP.x, 2.0),
                     HP.y * HP.z / HP.x);
}
dVector3D <double> GetV(const dVector3D <double>& HP, double gP) {
    return dVector3D <double>(HP.z,
                     HP.y * HP.z / HP.x,
                     pow(HP.z, 2.0) / HP.x + 0.5 * gP * pow(HP.x, 2.0));
}

dVector3D <double> GetHiHalf(const double DeltaTP, const double DeltaXP, const std :: vector <std :: vector <dVector3D <double>>>& GridP, unsigned int i, unsigned int j, double gP) {
    dVector3D <double> H_iplus1_jL = GetH(GridP, i + 1, j);
    dVector3D <double> H_i_jL = GetH(GridP, i, j);

    dVector3D <double> U_iplus1_jL = GetU(H_iplus1_jL, gP);
    dVector3D <double> U_i_jL = GetU(H_i_jL, gP);

    return (H_iplus1_jL + H_i_jL) / 2.0 - (U_iplus1_jL - U_i_jL) * (DeltaTP / DeltaXP) / 2.0;
}
dVector3D <double> GetHjHalf(const double DeltaTP, const double DeltaYP, const std :: vector <std :: vector <dVector3D <double>>>& GridP, unsigned int i, unsigned int j, double gP) {
    dVector3D <double> H_i_jplus1L = GetH(GridP, i, j + 1);
    dVector3D <double> H_i_jL = GetH(GridP, i, j);

    dVector3D <double> V_i_jplus1L = GetV(H_i_jplus1L, gP);
    dVector3D <double> V_i_jL = GetV(H_i_jL, gP);

    return (H_i_jplus1L + H_i_jL) / 2.0 - (V_i_jplus1L - V_i_jL) * (DeltaTP / DeltaYP) / 2.0;
}

dVector3D <double> NextH(const double DeltaTP, const double DeltaXP, const double DeltaYP, const std :: vector <std :: vector <dVector3D <double>>>& GridP, unsigned int i, unsigned int j, double gP) {
    dVector3D <double> H_iplushalf_jL = GetHiHalf(DeltaTP, DeltaXP, GridP, i, j, gP);
    dVector3D <double> H_iminushalf_jL = GetHiHalf(DeltaTP, DeltaXP, GridP, i - 1, j, gP);

    dVector3D <double> H_i_jplushalfL = GetHjHalf(DeltaTP, DeltaYP, GridP, i, j, gP);
    dVector3D <double> H_i_jminushalfL = GetHjHalf(DeltaTP, DeltaYP, GridP, i, j - 1, gP);

    dVector3D <double> U_iplushalf_jL = GetU(H_iplushalf_jL, gP);
    dVector3D <double> U_iminushalf_jL = GetU(H_iminushalf_jL, gP);

    dVector3D <double> V_i_jplushalfL = GetV(H_i_jplushalfL, gP);
    dVector3D <double> V_i_jminushalfL = GetV(H_i_jminushalfL, gP);

    dVector3D <double> H_i_jL = GetH(GridP, i, j);

    return H_i_jL - (U_iplushalf_jL - U_iminushalf_jL) * (DeltaTP / DeltaXP) - (V_i_jplushalfL - V_i_jminushalfL) * (DeltaTP / DeltaYP);
}
