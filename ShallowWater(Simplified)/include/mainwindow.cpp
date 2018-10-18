#include "mainwindow.h"
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
double dRichtmyer::GetTimeStep() {
    return TimeStep;
}

double dRichtmyer::GetXStep() {
    return xStep;
}
double dRichtmyer::GetYStep() {
    return yStep;
}
double dRichtmyer::GetZStep() {
    return zStep;
}
//-----------------------------//
dVectorND <double> dRichtmyer::Solve1D( dVectorND <double> U,
                                        dVectorND <double> Ux_minus_1,
                                        dVectorND <double> Ux_plus_1) {
    if (TimeStep == 0.0 || xStep == 0.0) {
        std::cout << "Error: set steps!" << std::endl;
        exit(-1);
    }

    dVectorND <double> UxHalf_plusL = UxHalfVector(U, Ux_plus_1);
    dVectorND <double> UxHalf_minusL = UxHalfVector(Ux_minus_1, U);

    return  U -
            (xFunc(UxHalf_plusL) - xFunc(UxHalf_minusL)) * (TimeStep / xStep);
}
dVectorND <double> dRichtmyer::Solve2D( dVectorND <double> U,
                                        dVectorND <double> Ux_minus_1,
                                        dVectorND <double> Ux_plus_1,
                                        dVectorND <double> Uy_minus_1,
                                        dVectorND <double> Uy_plus_1) {
    if (TimeStep == 0.0 || xStep == 0.0 || yStep == 0.0) {
        std::cout << "Error: set steps!" << std::endl;
        exit(-1);
    }

    dVectorND <double> UxHalf_plusL = UxHalfVector(U, Ux_plus_1);
    dVectorND <double> UxHalf_minusL = UxHalfVector(Ux_minus_1, U);

    dVectorND <double> UyHalf_plusL = UyHalfVector(U, Uy_plus_1);
    dVectorND <double> UyHalf_minusL = UyHalfVector(Uy_minus_1, U);

    return  U -
            (xFunc(UxHalf_plusL) - xFunc(UxHalf_minusL)) * (TimeStep / xStep) -
            (yFunc(UyHalf_plusL) - yFunc(UyHalf_minusL)) * (TimeStep / yStep);
}
dVectorND <double> dRichtmyer::Solve3D( dVectorND <double> U,
                                        dVectorND <double> Ux_minus_1,
                                        dVectorND <double> Ux_plus_1,
                                        dVectorND <double> Uy_minus_1,
                                        dVectorND <double> Uy_plus_1,
                                        dVectorND <double> Uz_minus_1,
                                        dVectorND <double> Uz_plus_1) {
    if (TimeStep == 0.0 || xStep == 0.0 || yStep == 0.0 || zStep == 0.0) {
        std::cout << "Error: set steps!" << std::endl;
        exit(-1);
    }

    dVectorND <double> UxHalf_plusL = UxHalfVector(U, Ux_plus_1);
    dVectorND <double> UxHalf_minusL = UxHalfVector(Ux_minus_1, U);

    dVectorND <double> UyHalf_plusL = UyHalfVector(U, Uy_plus_1);
    dVectorND <double> UyHalf_minusL = UyHalfVector(Uy_minus_1, U);

    dVectorND <double> UzHalf_plusL = UzHalfVector(U, Uz_plus_1);
    dVectorND <double> UzHalf_minusL = UzHalfVector(Uz_minus_1, U);

    return  U -
            (xFunc(UxHalf_plusL) - xFunc(UxHalf_minusL)) * (TimeStep / xStep) -
            (yFunc(UyHalf_plusL) - yFunc(UyHalf_minusL)) * (TimeStep / yStep) -
            (zFunc(UzHalf_plusL) - zFunc(UzHalf_minusL)) * (TimeStep / zStep);
}
//-----------------------------//
dVectorND <double> dRichtmyer::UxHalfVector(const dVectorND <double>& Ux, const dVectorND <double>& Ux_plus_1) {
    return (Ux_plus_1 + Ux) / 2.0 - (xFunc(Ux_plus_1) - xFunc(Ux)) * (TimeStep / 2.0 / xStep);
}
dVectorND <double> dRichtmyer::UyHalfVector(const dVectorND <double>& Uy, const dVectorND <double>& Uy_plus_1) {
    return (Uy_plus_1 + Uy) / 2.0 - (yFunc(Uy_plus_1) - yFunc(Uy)) * (TimeStep / 2.0 / yStep);
}
dVectorND <double> dRichtmyer::UzHalfVector(const dVectorND <double>& Uz, const dVectorND <double>& Uz_plus_1) {
    return (Uz_plus_1 + Uz) / 2.0 - (zFunc(Uz_plus_1) - zFunc(Uz)) * (TimeStep / 2.0 / zStep);
}
//-----------------------------//
//-----------------------------//
//-----------------------------//
dRichtmyerSolver::dRichtmyerSolver(double TimeStepP) {
    TimeStep = TimeStepP;
}
dRichtmyerSolver::dRichtmyerSolver(double TimeStepP, double xStepP) {
    TimeStep = TimeStepP;
    xStep = xStepP;
}
dRichtmyerSolver::dRichtmyerSolver(double TimeStepP, double xStepP, double yStepP) {
    TimeStep = TimeStepP;
    xStep = xStepP;
    yStep = yStepP;
}
dRichtmyerSolver::dRichtmyerSolver(double TimeStepP, double xStepP, double yStepP, double zStepP) {
    TimeStep = TimeStepP;
    xStep = xStepP;
    yStep = yStepP;
    zStep = zStepP;
}
//-----------------------------//
dVectorND <double> dRichtmyerSolver::xFunc(const dVectorND <double>& U) {
//    return dVectorND <double>({U[1], pow(U[1], 2.0) / U[0] + 0.5 * g * pow(U[0], 2.0), U[1] * U[2] / U[0]});
    return dVectorND <double> ({U[1],
                                (pow(U[1], 2.0) - pow(U[3], 2.0)) / U[0] + 0.5 * g * pow(U[0], 2.0),
                                (U[1] * U[2] - U[3] * U[4]) / U[0],
                                0,
                                (U[1] * U[4] - U[2] * U[3]) / U[0],
                                B_0 * U[1] / U[0]});
}
dVectorND <double> dRichtmyerSolver::yFunc(const dVectorND <double>& U) {
//    return dVectorND <double>({U[2], U[1] * U[2] / U[0], pow(U[2], 2.0) / U[0] + 0.5 * g * pow(U[0], 2.0)});
    return dVectorND <double> ({U[2],
                                (U[1] * U[2] - U[3] * U[4]) / U[0],
                                (pow(U[2], 2.0) - pow(U[4], 2.0)) / U[0] + 0.5 * g * pow(U[0], 2.0),
                                (U[2] * U[3] - U[1] * U[4]) / U[0],
                                0,
                                B_0 * U[2] / U[0]});
}
dVectorND <double> dRichtmyerSolver::zFunc(const dVectorND <double>& U) {}



















Solver::Solver() {
    Test = new dRichtmyerSolver(0.001, 0.1, 0.1);
}
Solver::~Solver() {
    delete Test;
}

void Solver :: Calc(unsigned int StepsP) {
    unsigned int PointNum = 200;

    InitGrid(10.0, 0.0, 0.0, PointNum);

    Grid[100][100][0] = 12.0;

    std :: vector <std :: vector <dVectorND <double>>> TempGrid = Grid;

    for (unsigned int iter = 0; iter < StepsP; iter++) {
        double Energy = 0.0;

        for (unsigned int i = 1; i < PointNum - 1; i++) {
            for (unsigned int j = 1; j < PointNum - 1; j++) {
                dVectorND <double> TempH = Test -> Solve2D(Grid[i][j], Grid[i - 1][j], Grid[i + 1][j], Grid[i][j - 1], Grid[i][j + 1]);

                TempGrid[i][j] = TempH;

                Energy += ( pow((Grid[i][j][0] - TempGrid[i][j][0]) / Test -> GetTimeStep(), 2.0) +
                            pow(TempGrid[i][j][1] / TempGrid[i][j][0], 2.0) +
                            pow(TempGrid[i][j][2] / TempGrid[i][j][0], 2.0));
            }
        }

        for (unsigned int i = 0; i < PointNum; i++) {
            TempGrid[i][0] = TempGrid[i][1];
            TempGrid[i][PointNum - 1] = TempGrid[i][PointNum - 2];
            TempGrid[0][i] = TempGrid[1][i];
            TempGrid[PointNum - 1][i] = TempGrid[PointNum - 2][i];
        }

        for (unsigned int i = 0; i < PointNum; i++) {
            Energy += ( pow((Grid[i][0][0] - TempGrid[i][0][0]) / Test -> GetTimeStep(), 1.0) +
                        pow(TempGrid[i][0][1] / TempGrid[i][0][0], 2.0) +
                        pow(TempGrid[i][0][2] / TempGrid[i][0][0], 2.0));
        }
        for (unsigned int i = 1; i < PointNum - 1; i++) {
            Energy += ( pow((Grid[0][i][0] - TempGrid[0][i][0] / Test -> GetTimeStep()), 2.0) +
                        pow(TempGrid[0][i][1] / TempGrid[0][i][0], 2.0) +
                        pow(TempGrid[0][i][2] / TempGrid[0][i][0], 2.0));
        }

        std::cout << Energy << " : " << iter << std::endl;

        Grid = TempGrid;
    }
}
//-----------------------------//
MainWindow::MainWindow(QWidget *parent) : QMainWindow(parent) {
    TestSolver = new Solver;

    TestSolver -> Calc(100);

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

            auto y = float(TestSolver -> Grid[i][j][0]);
            (*newRow)[index++].setPosition(QVector3D(x, y, z));
        }
        *dataArray << newRow;
    }

    auto series = new QSurface3DSeries;
    series->setDrawMode(QSurface3DSeries::DrawSurface);
    series->dataProxy()->resetArray(dataArray);
    graph -> addSeries(series);

    graph -> renderToImage(0, QSize(1000, 1000)).save("Plot.png");

    MainLayout -> addWidget(container);
}
//-----------------------------//
void Solver::InitGrid(double ExcitationP, double VXP, double VYP, int PointsNumP) {
    for (int i = 0; i < PointsNumP; i++) {
        Grid.emplace_back(std :: vector <dVectorND <double>>());

        for (int j = 0; j < PointsNumP; j++) {
            Grid.back().emplace_back(dVectorND <double>({ExcitationP, VXP * ExcitationP, VYP * ExcitationP, 0, 0, 0}));
        }
    }
}