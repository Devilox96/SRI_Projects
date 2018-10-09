#ifndef MAINWINDOW_H
#define MAINWINDOW_H
//-----------------------------//
#include <QMainWindow>
#include <QGridLayout>
#include <QtDataVisualization/Q3DSurface>
#include <QtDataVisualization/QSurfaceDataProxy>
#include <QtDataVisualization/QHeightMapSurfaceDataProxy>
#include <QtDataVisualization/QSurface3DSeries>
#include <QtMath>
//-----------------------------//
#include "Libs/dMath/Core/dVectors.h"
//-----------------------------//
using namespace QtDataVisualization;
//-----------------------------//
class dRichtmyer {
public:
    dRichtmyer() = default;
    explicit dRichtmyer(double TimeStepP);
    dRichtmyer(double TimeStepP, double xStepP);
    dRichtmyer(double TimeStepP, double xStepP, double yStepP);
    dRichtmyer(double TimeStepP, double xStepP, double yStepP, double zStepP);
    ~dRichtmyer() = default;

    //----------//

    void SetTimeStep(double TimeStepP);

    void SetXStep(double xStepP);
    void SetYStep(double yStepP);
    void SetZStep(double zStepP);

    //----------//

    dVector3D <double> Solve1D( dVector3D <double> U,
                                dVector3D <double> Ux_minus_1,
                                dVector3D <double> Ux_plus_1,
                                dVector3D <double> X_minus_1,
                                dVector3D <double> X,
                                dVector3D <double> X_plus_1);
    dVector3D <double> Solve2D( dVector3D <double> U,
                                dVector3D <double> Ux_minus_1,
                                dVector3D <double> Ux_plus_1,
                                dVector3D <double> Uy_minus_1,
                                dVector3D <double> Uy_plus_1,
                                dVector3D <double> X_minus_1,
                                dVector3D <double> X,
                                dVector3D <double> X_plus_1,
                                dVector3D <double> Y_minus_1,
                                dVector3D <double> Y,
                                dVector3D <double> Y_plus_1);
    dVector3D <double> Solve3D( dVector3D <double> U,
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
                                dVector3D <double> Z_plus_1);
private:
    double TimeStep = 0.0;

    double xStep = 0.0;
    double yStep = 0.0;
    double zStep = 0.0;

    //----------//

    dVector3D <double> HalfStepVector(  dVector3D <double> U,
                                        dVector3D <double> U_plus_1,
                                        dVector3D <double> F,
                                        dVector3D <double> F_plus_1,
                                        double CoordStepP);
};
//-----------------------------//
class Solver {
public:
    Solver();
    ~Solver();

    void Calc(unsigned int StepsP);

    std :: vector <std :: vector <dVector3D <double>>> Grid;
    dRichtmyer* Test;
private:
    double TimeStep = 0.001;

    double xStep = 0.1;
    double yStep = 0.1;
    double zStep = 0.0;

    double g = 9.81;




    void InitGrid(double ExcitationP, double VXP, double VYP, int PointsNumP);

    dVector3D <double> GetU(const dVector3D <double>& HP);
    dVector3D <double> GetV(const dVector3D <double>& HP);

    dVector3D <double> GetHiHalf(unsigned int i, unsigned int j);
    dVector3D <double> GetHjHalf(unsigned int i, unsigned int j);

    dVector3D <double> NextH(unsigned int i, unsigned int j);
};
//-----------------------------//
class MainWindow : public QMainWindow {
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = nullptr);
    ~MainWindow() override = default;

    Q3DSurface *graph;
    QWidget* MainWidget;
    QGridLayout* MainLayout;

    Solver* TestSolver;
};
//-----------------------------//
#endif
