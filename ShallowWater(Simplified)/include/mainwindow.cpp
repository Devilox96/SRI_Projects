#include "mainwindow.h"
//-----------------------------//
dRichtmyer::dRichtmyer(double TimeStepP) : TimeStep(TimeStepP) {}
dRichtmyer::dRichtmyer(double TimeStepP, double xStepP) : TimeStep(TimeStepP), xStep(xStepP) {}
dRichtmyer::dRichtmyer(double TimeStepP, double xStepP, double yStepP) :
        TimeStep(TimeStepP), xStep(xStepP), yStep(yStepP) {}
dRichtmyer::dRichtmyer(double TimeStepP, double xStepP, double yStepP, double zStepP) :
        TimeStep(TimeStepP), xStep(xStepP), yStep(yStepP), zStep(zStepP) {}
//-----------------------------//
void dRichtmyer::SetTimeStep(double TimeStepP) {
    TimeStep = TimeStepP;
}

void dRichtmyer::SetXStep(double xStepP) {
    xStep = xStepP;
}
void dRichtmyer::SetYStep(double yStepP) {
    yStep = yStepP;
}
void dRichtmyer::SetZStep(double zStepP) {
    zStep = zStepP;
}
//-----------------------------//
dVector3D <double> dRichtmyer::Solve1D( dVector3D <double> U,
                                        dVector3D <double> Ux_minus_1,
                                        dVector3D <double> Ux_plus_1,
                                        dVector3D <double> X_minus_1,
                                        dVector3D <double> X,
                                        dVector3D <double> X_plus_1) {
    if (TimeStep == 0.0 || xStep == 0.0) {
        std::cout << "Error: set steps!" << std::endl;
        exit(-1);
    }

    dVector3D <double> HalfStepX_plus = HalfStepVector(U, Ux_plus_1, X, X_plus_1, xStep);
    dVector3D <double> HalfStepX_minus = HalfStepVector(Ux_minus_1, U, X_minus_1, X, xStep);

    return  U -
            (HalfStepX_plus - HalfStepX_minus) * (TimeStep / xStep);
}
dVector3D <double> dRichtmyer::Solve2D( dVector3D <double> U,
                                        dVector3D <double> Ux_minus_1,
                                        dVector3D <double> Ux_plus_1,
                                        dVector3D <double> Uy_minus_1,
                                        dVector3D <double> Uy_plus_1,
                                        dVector3D <double> X_minus_1,
                                        dVector3D <double> X,
                                        dVector3D <double> X_plus_1,
                                        dVector3D <double> Y_minus_1,
                                        dVector3D <double> Y,
                                        dVector3D <double> Y_plus_1) {
    if (TimeStep == 0.0 || xStep == 0.0 || yStep == 0.0) {
        std::cout << "Error: set steps!" << std::endl;
        exit(-1);
    }

    dVector3D <double> HalfStepX_plus = HalfStepVector(U, Ux_plus_1, X, X_plus_1, xStep);
    dVector3D <double> HalfStepX_minus = HalfStepVector(Ux_minus_1, U, X_minus_1, X, xStep);

    dVector3D <double> HalfStepY_plus = HalfStepVector(U, Uy_plus_1, Y, Y_plus_1, yStep);
    dVector3D <double> HalfStepY_minus = HalfStepVector(Uy_minus_1, U, Y_minus_1, Y, yStep);

    return  U -
            (HalfStepX_plus - HalfStepX_minus) * (TimeStep / xStep) -
            (HalfStepY_plus - HalfStepY_minus) * (TimeStep / yStep);
}
dVector3D <double> dRichtmyer::Solve3D( dVector3D <double> U,
                                        dVector3D <double> Ux_minus_1,
                                        dVector3D <double> Ux_plus_1,
                                        dVector3D <double> Uy_minus_1,
                                        dVector3D <double> Uy_plus_1,
                                        dVector3D <double> Uz_minus_1,
                                        dVector3D <double> Uz_plus_1,
                                        dVector3D <double> X_minus_1,
                                        dVector3D <double> X,
                                        dVector3D <double> X_plus_1,
                                        dVector3D <double> Y_minus_1,
                                        dVector3D <double> Y,
                                        dVector3D <double> Y_plus_1,
                                        dVector3D <double> Z_minus_1,
                                        dVector3D <double> Z,
                                        dVector3D <double> Z_plus_1) {
    if (TimeStep == 0.0 || xStep == 0.0 || yStep == 0.0 || zStep == 0.0) {
        std::cout << "Error: set steps!" << std::endl;
        exit(-1);
    }

    dVector3D <double> HalfStepX_plus = HalfStepVector(U, Ux_plus_1, X, X_plus_1, xStep);
    dVector3D <double> HalfStepX_minus = HalfStepVector(Ux_minus_1, U, X_minus_1, X, xStep);

    dVector3D <double> HalfStepY_plus = HalfStepVector(U, Uy_plus_1, Y, Y_plus_1, yStep);
    dVector3D <double> HalfStepY_minus = HalfStepVector(Uy_minus_1, U, Y_minus_1, Y, yStep);

    dVector3D <double> HalfStepZ_plus = HalfStepVector(U, Uz_plus_1, Z, Z_plus_1, zStep);
    dVector3D <double> HalfStepZ_minus = HalfStepVector(Uz_minus_1, U, Z_minus_1, Z, zStep);

    return  U -
            (HalfStepX_plus - HalfStepX_minus) * (TimeStep / xStep) -
            (HalfStepY_plus - HalfStepY_minus) * (TimeStep / yStep) -
            (HalfStepZ_plus - HalfStepZ_minus) * (TimeStep / zStep);
}
//-----------------------------//
dVector3D <double> dRichtmyer::HalfStepVector(  dVector3D <double> U,
                                                dVector3D <double> U_plus_1,
                                                dVector3D <double> F,
                                                dVector3D <double> F_plus_1,
                                                double CoordStepP) {
    return (U_plus_1 + U) / 2.0 - (F_plus_1 - F) * (TimeStep / 2.0 / CoordStepP);
}
//-----------------------------//






















Solver::Solver() {
    Test = new dRichtmyer(0.001, 0.1, 0.1);
}
Solver::~Solver() {
    delete Test;
}


void Solver :: Calc(unsigned int StepsP) {
    unsigned int PointNum = 200;

    InitGrid(10.0, 0.0, 0.0, PointNum);

    Grid[100][100].x = 12.0;

    std :: vector <std :: vector <dVector3D <double>>> TempGrid = Grid;

    for (unsigned int iter = 0; iter < StepsP; iter++) {
        double Energy = 0.0;

        for (unsigned int i = 1; i < PointNum - 1; i++) {
            for (unsigned int j = 1; j < PointNum - 1; j++) {
//                dVector3D <double> TempH = NextH(i, j);

                dVector3D <double> TempH = Test -> Solve2D( Grid[i][j],
                                                            Grid[i - 1][j],
                                                            Grid[i + 1][j],
                                                            Grid[i][j - 1],
                                                            Grid[i][j + 1],
                                                            GetU(Grid[i - 1][j]),
                                                            GetU(Grid[i][j]))

                TempGrid[i][j] = TempH;

                Energy += ( pow((Grid[i][j].x - TempGrid[i][j].x) / TimeStep, 2.0) +
                            pow(TempGrid[i][j].y / TempGrid[i][j].x, 2.0) +
                            pow(TempGrid[i][j].z / TempGrid[i][j].x, 2.0));
            }
        }

        for (unsigned int i = 0; i < PointNum; i++) {
            TempGrid[i][0] = TempGrid[i][1];
            TempGrid[i][PointNum - 1] = TempGrid[i][PointNum - 2];
            TempGrid[0][i] = TempGrid[1][i];
            TempGrid[PointNum - 1][i] = TempGrid[PointNum - 2][i];
        }

        for (unsigned int i = 0; i < PointNum; i++) {
            Energy += ( pow((Grid[i][0].x - TempGrid[i][0].x) / TimeStep, 1.0) +
                        pow(TempGrid[i][0].y / TempGrid[i][0].x, 2.0) +
                        pow(TempGrid[i][0].z / TempGrid[i][0].x, 2.0));
        }
        for (unsigned int i = 1; i < PointNum - 1; i++) {
            Energy += ( pow((Grid[0][i].x - TempGrid[0][i].x / TimeStep), 2.0) +
                        pow(TempGrid[0][i].y / TempGrid[0][i].x, 2.0) +
                        pow(TempGrid[0][i].z / TempGrid[0][i].x, 2.0));
        }

        std::cout << Energy << " : " << iter << std::endl;

        Grid = TempGrid;
    }
}
//-----------------------------//
MainWindow::MainWindow(QWidget *parent) : QMainWindow(parent) {
    TestSolver = new Solver;

    TestSolver -> Calc(2000);

    MainWidget = new QWidget;
    setCentralWidget(MainWidget);

    MainLayout = new QGridLayout;
    MainWidget -> setLayout(MainLayout);

    graph = new Q3DSurface;
    QWidget *container = QWidget :: createWindowContainer(graph);
    container -> setMinimumSize(800, 800);
    graph -> setFlags(graph -> flags() ^ Qt :: FramelessWindowHint);

    int sampleCountX = 200;
    int sampleCountZ = 200;

    auto dataArray = new QSurfaceDataArray;
    dataArray->reserve(sampleCountZ);

    for (int i = 0 ; i < sampleCountZ ; i++) {
        QSurfaceDataRow *newRow = new QSurfaceDataRow(sampleCountX);

        float z = i;
        int index = 0;
        for (int j = 0; j < sampleCountX; j++) {
            float x = j;

            auto y = float(TestSolver -> Grid[i][j].x);
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
void Solver::InitGrid(double ExcitationP, double VXP, double VYP, int PointsNumP) {
    for (int i = 0; i < PointsNumP; i++) {
        Grid.emplace_back(std :: vector <dVector3D <double>>());

        for (int j = 0; j < PointsNumP; j++) {
            Grid.back().emplace_back(dVector3D <double>(ExcitationP, VXP * ExcitationP, VYP * ExcitationP));
        }
    }
}

dVector3D <double> Solver::GetU(const dVector3D <double>& HP) {
    return dVector3D <double>(HP.y, pow(HP.y, 2.0) / HP.x + 0.5 * g * pow(HP.x, 2.0), HP.y * HP.z / HP.x);
}
dVector3D <double> Solver::GetV(const dVector3D <double>& HP) {
    return dVector3D <double>(HP.z, HP.y * HP.z / HP.x, pow(HP.z, 2.0) / HP.x + 0.5 * g * pow(HP.x, 2.0));
}

dVector3D <double> Solver::GetHiHalf(unsigned int i, unsigned int j) {
    dVector3D <double> H_iplus1_jL = Grid[i + 1][j];
    dVector3D <double> H_i_jL = Grid[i][j];

    dVector3D <double> U_iplus1_jL = GetU(H_iplus1_jL);
    dVector3D <double> U_i_jL = GetU(H_i_jL);

    return (H_iplus1_jL + H_i_jL) / 2.0 - (U_iplus1_jL - U_i_jL) * (TimeStep / xStep) / 2.0;
}
dVector3D <double> Solver::GetHjHalf(unsigned int i, unsigned int j) {
    dVector3D <double> H_i_jplus1L = Grid[i][j + 1];
    dVector3D <double> H_i_jL = Grid[i][j];

    dVector3D <double> V_i_jplus1L = GetV(H_i_jplus1L);
    dVector3D <double> V_i_jL = GetV(H_i_jL);

    return (H_i_jplus1L + H_i_jL) / 2.0 - (V_i_jplus1L - V_i_jL) * (TimeStep / yStep) / 2.0;
}

dVector3D <double> Solver::NextH(unsigned int i, unsigned int j) {
    dVector3D <double> H_iplushalf_jL = GetHiHalf(i, j);
    dVector3D <double> H_iminushalf_jL = GetHiHalf(i - 1, j);

    dVector3D <double> H_i_jplushalfL = GetHjHalf(i, j);
    dVector3D <double> H_i_jminushalfL = GetHjHalf(i, j - 1);

    dVector3D <double> U_iplushalf_jL = GetU(H_iplushalf_jL);
    dVector3D <double> U_iminushalf_jL = GetU(H_iminushalf_jL);

    dVector3D <double> V_i_jplushalfL = GetV(H_i_jplushalfL);
    dVector3D <double> V_i_jminushalfL = GetV(H_i_jminushalfL);

    dVector3D <double> H_i_jL = Grid[i][j];

    return H_i_jL - (U_iplushalf_jL - U_iminushalf_jL) * (TimeStep / xStep) - (V_i_jplushalfL - V_i_jminushalfL) * (TimeStep / yStep);
}
