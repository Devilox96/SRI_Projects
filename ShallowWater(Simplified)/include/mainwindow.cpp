#include "mainwindow.h"

void MainWindow :: Calc(unsigned int StepsP) {
    unsigned int PointNum = 100;

    InitGrid(Grid, 10.0, 0.0, 0.0, PointNum);

    Grid[50][50][0] = 12.0;

    std :: vector <std :: vector <dVectorND <double>>> TempGrid = Grid;

    for (unsigned int iter = 0; iter < StepsP; iter++) {
        for (unsigned int i = 1; i < PointNum - 1; i++) {
            for (unsigned int j = 1; j < PointNum - 1; j++) {
                dVectorND <double> TempH = NextH(0.001, 0.1, 0.1, Grid, i, j, 9.81);

                TempGrid[i][j][0] = TempH[0];
                TempGrid[i][j][1] = TempH[1] / TempH[0];
                TempGrid[i][j][2] = TempH[2] / TempH[0];
            }
        }

        for (unsigned int i = 0; i < PointNum; i++) {
            TempGrid[i][0] = TempGrid[i][1];
            TempGrid[i][PointNum - 1] = TempGrid[i][PointNum - 2];
            TempGrid[0][i] = TempGrid[1][i];
            TempGrid[PointNum - 1][i] = TempGrid[PointNum - 2][i];
        }

        Grid = TempGrid;

        if (iter % 10 == 0) {
            std::cout << iter << std::endl;
        }
    }
}
//-----------------------------//
MainWindow::MainWindow(QWidget *parent) : QMainWindow(parent) {
    auto start = std::chrono::system_clock::now();

    Calc(50);

    auto end = std::chrono::system_clock::now();

    std::chrono::duration <double> elapsed_seconds = end - start;

    std::cout << elapsed_seconds.count() << std::endl;

    //--------------------------------------------------------------//

    MainWidget = new QWidget;
    setCentralWidget(MainWidget);

    MainLayout = new QGridLayout;
    MainWidget -> setLayout(MainLayout);

    graph = new Q3DSurface;
    QWidget *container = QWidget :: createWindowContainer(graph);
    container -> setMinimumSize(1280, 720);
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

            auto y = float(Grid[i][j][0]);
            (*newRow)[index++].setPosition(QVector3D(x, y, z));
        }
        *dataArray << newRow;
    }

    auto series = new QSurface3DSeries;
    series -> setDrawMode(QSurface3DSeries::DrawSurface);
    series -> dataProxy()->resetArray(dataArray);
    graph -> addSeries(series);

    MainLayout -> addWidget(container);
}
//-----------------------------//
void InitGrid(std :: vector <std :: vector <dVectorND <double>>>& GridP, double ExcitationP, double VXP, double VYP, int PointsNumP) {
    for (int i = 0; i < PointsNumP; i++) {
        GridP.emplace_back(std :: vector <dVectorND <double>>());

        for (int j = 0; j < PointsNumP; j++) {
            GridP.back().emplace_back(dVectorND <double>({ExcitationP, VXP, VYP}));
        }
    }
}

dVectorND <double> GetH(const std :: vector <std :: vector <dVectorND <double>>>& GridP, unsigned int i, unsigned int j) {
    return dVectorND <double>({GridP[i][j][0],
                     GridP[i][j][1] * GridP[i][j][0],
                     GridP[i][j][2] * GridP[i][j][0]});
}
dVectorND <double> GetU(const dVectorND <double>& HP, double gP) {
    return dVectorND <double>({HP[1],
                     pow(HP[1], 2.0) / HP[0] + 0.5 * gP * pow(HP[0], 2.0),
                     HP[1] * HP[2] / HP[0]});
}
dVectorND <double> GetV(const dVectorND <double>& HP, double gP) {
    return dVectorND <double>({HP[2],
                     HP[1] * HP[2] / HP[0],
                     pow(HP[2], 2.0) / HP[0] + 0.5 * gP * pow(HP[0], 2.0)});
}

dVectorND <double> GetHiHalf(const double DeltaTP, const double DeltaXP, const std :: vector <std :: vector <dVectorND <double>>>& GridP, unsigned int i, unsigned int j, double gP) {
    dVectorND <double> H_iplus1_jL = GetH(GridP, i + 1, j);
    dVectorND <double> H_i_jL = GetH(GridP, i, j);

    dVectorND <double> U_iplus1_jL = GetU(H_iplus1_jL, gP);
    dVectorND <double> U_i_jL = GetU(H_i_jL, gP);

    return (H_iplus1_jL + H_i_jL) / 2.0 - (U_iplus1_jL - U_i_jL) * (DeltaTP / DeltaXP) / 2.0;
}
dVectorND <double> GetHjHalf(const double DeltaTP, const double DeltaYP, const std :: vector <std :: vector <dVectorND <double>>>& GridP, unsigned int i, unsigned int j, double gP) {
    dVectorND <double> H_i_jplus1L = GetH(GridP, i, j + 1);
    dVectorND <double> H_i_jL = GetH(GridP, i, j);

    dVectorND <double> V_i_jplus1L = GetV(H_i_jplus1L, gP);
    dVectorND <double> V_i_jL = GetV(H_i_jL, gP);

    return (H_i_jplus1L + H_i_jL) / 2.0 - (V_i_jplus1L - V_i_jL) * (DeltaTP / DeltaYP) / 2.0;
}

dVectorND <double> NextH(const double DeltaTP, const double DeltaXP, const double DeltaYP, const std :: vector <std :: vector <dVectorND <double>>>& GridP, unsigned int i, unsigned int j, double gP) {
    dVectorND <double> H_iplushalf_jL = GetHiHalf(DeltaTP, DeltaXP, GridP, i, j, gP);
    dVectorND <double> H_iminushalf_jL = GetHiHalf(DeltaTP, DeltaXP, GridP, i - 1, j, gP);

    dVectorND <double> H_i_jplushalfL = GetHjHalf(DeltaTP, DeltaYP, GridP, i, j, gP);
    dVectorND <double> H_i_jminushalfL = GetHjHalf(DeltaTP, DeltaYP, GridP, i, j - 1, gP);

    return  GetH(GridP, i, j) -
            (GetU(H_iplushalf_jL, gP) - GetU(H_iminushalf_jL, gP)) * (DeltaTP / DeltaXP) -
            (GetV(H_i_jplushalfL, gP) - GetV(H_i_jminushalfL, gP)) * (DeltaTP / DeltaYP);
}
