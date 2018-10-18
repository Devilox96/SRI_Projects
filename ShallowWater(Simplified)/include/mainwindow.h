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
    ~dRichtmyer() = default;

    //----------//

    void SetTimeStep(double TimeStepP);

    void SetXStep(double xStepP);
    void SetYStep(double yStepP);
    void SetZStep(double zStepP);

    //----------//

    double GetTimeStep();

    double GetXStep();
    double GetYStep();
    double GetZStep();

    //----------//

    dVectorND <double> Solve1D( dVectorND <double> U,
                                dVectorND <double> Ux_minus_1,
                                dVectorND <double> Ux_plus_1);
    dVectorND <double> Solve2D( dVectorND <double> U,
                                dVectorND <double> Ux_minus_1,
                                dVectorND <double> Ux_plus_1,
                                dVectorND <double> Uy_minus_1,
                                dVectorND <double> Uy_plus_1);
    dVectorND <double> Solve3D( dVectorND <double> U,
                                dVectorND <double> Ux_minus_1,
                                dVectorND <double> Ux_plus_1,
                                dVectorND <double> Uy_minus_1,
                                dVectorND <double> Uy_plus_1,
                                dVectorND <double> Uz_minus_1,
                                dVectorND <double> Uz_plus_1);
protected:
    double TimeStep = 0.0;

    double xStep = 0.0;
    double yStep = 0.0;
    double zStep = 0.0;

    //----------//

    virtual dVectorND <double> xFunc(const dVectorND <double>& U) = 0;
    virtual dVectorND <double> yFunc(const dVectorND <double>& U) = 0;
    virtual dVectorND <double> zFunc(const dVectorND <double>& U) = 0;

    //----------//

    dVectorND <double> UxHalfVector(const dVectorND <double>& Ux, const dVectorND <double>& Ux_plus_1);
    dVectorND <double> UyHalfVector(const dVectorND <double>& Uy, const dVectorND <double>& Uy_plus_1);
    dVectorND <double> UzHalfVector(const dVectorND <double>& Uz, const dVectorND <double>& Uz_plus_1);
};
//-----------------------------//
class dRichtmyerSolver : public dRichtmyer {
public:
    dRichtmyerSolver() = default;
    explicit dRichtmyerSolver(double TimeStepP);
    dRichtmyerSolver(double TimeStepP, double xStepP);
    dRichtmyerSolver(double TimeStepP, double xStepP, double yStepP);
    dRichtmyerSolver(double TimeStepP, double xStepP, double yStepP, double zStepP);
    ~dRichtmyerSolver() = default;
private:
    const double g = 9.81;
    const double B_0 = 0;

    //----------//

    dVectorND <double> xFunc(const dVectorND <double>& U) override;
    dVectorND <double> yFunc(const dVectorND <double>& U) override;
    dVectorND <double> zFunc(const dVectorND <double>& U) override;
};
//-----------------------------//
class Solver {
public:
    Solver();
    ~Solver();

    void Calc(unsigned int StepsP);

    std :: vector <std :: vector <dVectorND <double>>> Grid;
    dRichtmyerSolver* Test;
private:
    void InitGrid(double ExcitationP, double VXP, double VYP, int PointsNumP);
};
//-----------------------------//
class MainWindow : public QMainWindow {
    Q_OBJECT
public:
    explicit MainWindow(QWidget *parent = nullptr);
    ~MainWindow() override = default;
private:
    Q3DSurface *graph;
    QWidget* MainWidget;
    QGridLayout* MainLayout;

    Solver* TestSolver;

    //----------//


};
//-----------------------------//
#endif
