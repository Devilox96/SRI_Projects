#include "mainwindow.h"

void MainWindow :: Calc(unsigned int StepsP) {
    unsigned int PointNum = 100;

    InitGrid(Grid, 10.0, 0.0, 0.0, PointNum);

    Grid[50][50].h = 12.0;

    Grid[49][50].h = 11.9;
    Grid[50][49].h = 11.9;
    Grid[50][51].h = 11.9;
    Grid[51][50].h = 11.9;

    Grid[49][49].h = 11.7;
    Grid[51][51].h = 11.7;
    Grid[49][51].h = 11.7;
    Grid[51][49].h = 11.7;

    std :: vector <std :: vector <Point>> TempGrid = Grid;

    for (unsigned int iter = 0; iter < StepsP; iter++) {
        for (unsigned int i = 1; i < PointNum - 1; i++) {
            for (unsigned int j = 1; j < PointNum - 1; j++) {
                dVector3D TempH = NextH(0.001, 0.1, 0.1, Grid, i, j, 9.81);

                TempGrid[i][j].h = TempH.x;
                TempGrid[i][j].u = TempH.y / TempH.x;
                TempGrid[i][j].v = TempH.z / TempH.x;
            }
        }

        for (unsigned int i = 0; i < PointNum; i++) {
            TempGrid[i][0].h = TempGrid[i][1].h;
            TempGrid[i][0].u = TempGrid[i][1].u;
            TempGrid[i][0].v = TempGrid[i][1].v;

            TempGrid[i][PointNum - 1].h = TempGrid[i][PointNum - 2].h;
            TempGrid[i][PointNum - 1].u = TempGrid[i][PointNum - 2].u;
            TempGrid[i][PointNum - 1].v = TempGrid[i][PointNum - 2].v;

            TempGrid[0][i].h = TempGrid[1][i].h;
            TempGrid[0][i].u = TempGrid[1][i].u;
            TempGrid[0][i].v = TempGrid[1][i].v;

            TempGrid[PointNum - 1][i].h = TempGrid[PointNum - 2][i].h;
            TempGrid[PointNum - 1][i].u = TempGrid[PointNum - 2][i].u;
            TempGrid[PointNum - 1][i].v = TempGrid[PointNum - 2][i].v;
        }

        Grid = TempGrid;
    }
}
//-----------------------------//
MainWindow::MainWindow(QWidget *parent) : QMainWindow(parent) {
    Calc(2000);

    MainWidget = new QWidget;
    setCentralWidget(MainWidget);

    MainLayout = new QGridLayout;
    MainWidget -> setLayout(MainLayout);

    graph = new Q3DSurface;
    QWidget *container = QWidget :: createWindowContainer(graph);
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

            auto y = float(Grid[i][j].h);
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
void InitGrid(std :: vector <std :: vector <Point>>& GridP, double ExcitationP, double VXP, double VYP, int PointsNumP) {
    for (int i = 0; i < PointsNumP; i++) {
        GridP.emplace_back(std :: vector <Point>());

        for (int j = 0; j < PointsNumP; j++) {
            GridP.back().emplace_back(Point(ExcitationP, VXP, VYP));
        }
    }
}

dVector3D GetH(const std :: vector <std :: vector <Point>>& GridP, unsigned int i, unsigned int j) {
    return dVector3D(GridP[i][j].h,
                     GridP[i][j].u * GridP[i][j].h,
                     GridP[i][j].v * GridP[i][j].h);
}
dVector3D GetU(const dVector3D& HP, double gP) {
    return dVector3D(HP.y,
                     pow(HP.y, 2.0) / HP.x + 0.5 * gP * pow(HP.x, 2.0),
                     HP.y * HP.z / HP.x);
}
dVector3D GetV(const dVector3D& HP, double gP) {
    return dVector3D(HP.z,
                     HP.y * HP.z / HP.x,
                     pow(HP.z, 2.0) / HP.x + 0.5 * gP * pow(HP.x, 2.0));
}

dVector3D GetHiHalf(const double DeltaTP, const double DeltaXP, const std :: vector <std :: vector <Point>>& GridP, unsigned int i, unsigned int j, double gP) {
    dVector3D H_iplus1_jL = GetH(GridP, i + 1, j);
    dVector3D H_i_jL = GetH(GridP, i, j);

    dVector3D U_iplus1_jL = GetU(H_iplus1_jL, gP);
    dVector3D U_i_jL = GetU(H_i_jL, gP);

    return (H_iplus1_jL + H_i_jL) / 2.0 - (U_iplus1_jL - U_i_jL) * (DeltaTP / DeltaXP) / 2.0;
}
dVector3D GetHjHalf(const double DeltaTP, const double DeltaYP, const std :: vector <std :: vector <Point>>& GridP, unsigned int i, unsigned int j, double gP) {
    dVector3D H_i_jplus1L = GetH(GridP, i, j + 1);
    dVector3D H_i_jL = GetH(GridP, i, j);

    dVector3D V_i_jplus1L = GetV(H_i_jplus1L, gP);
    dVector3D V_i_jL = GetV(H_i_jL, gP);

    return (H_i_jplus1L + H_i_jL) / 2.0 - (V_i_jplus1L - V_i_jL) * (DeltaTP / DeltaYP) / 2.0;
}

dVector3D NextH(const double DeltaTP, const double DeltaXP, const double DeltaYP, const std :: vector <std :: vector <Point>>& GridP, unsigned int i, unsigned int j, double gP) {
    dVector3D H_iplushalf_jL = GetHiHalf(DeltaTP, DeltaXP, GridP, i, j, gP);
    dVector3D H_iminushalf_jL = GetHiHalf(DeltaTP, DeltaXP, GridP, i - 1, j, gP);

    dVector3D H_i_jplushalfL = GetHjHalf(DeltaTP, DeltaYP, GridP, i, j, gP);
    dVector3D H_i_jminushalfL = GetHjHalf(DeltaTP, DeltaYP, GridP, i, j - 1, gP);

    dVector3D U_iplushalf_jL = GetU(H_iplushalf_jL, gP);
    dVector3D U_iminushalf_jL = GetU(H_iminushalf_jL, gP);

    dVector3D V_i_jplushalfL = GetV(H_i_jplushalfL, gP);
    dVector3D V_i_jminushalfL = GetV(H_i_jminushalfL, gP);

    dVector3D H_i_jL = GetH(GridP, i, j);

    return H_i_jL - (U_iplushalf_jL - U_iminushalf_jL) * (DeltaTP / DeltaXP) - (V_i_jplushalfL - V_i_jminushalfL) * (DeltaTP / DeltaYP);
}
