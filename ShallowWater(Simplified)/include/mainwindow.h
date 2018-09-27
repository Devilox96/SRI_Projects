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
class Solver {
public:
    Solver() = default;
    ~Solver() = default;

    void Calc(unsigned int StepsP);

    std :: vector <std :: vector <dVector3D <double>>> Grid;
private:
    double TimeStep = 0.001;

    double xStep = 0.1;
    double yStep = 0.1;
    double zStep = 0.0;

    double g = 9.81;




    void InitGrid(double ExcitationP, double VXP, double VYP, int PointsNumP);

    dVector3D <double> GetU(const dVector3D <double>& HP, double gP);
    dVector3D <double> GetV(const dVector3D <double>& HP, double gP);

    dVector3D <double> GetHiHalf(unsigned int i, unsigned int j, double gP);
    dVector3D <double> GetHjHalf(unsigned int i, unsigned int j, double gP);

    dVector3D <double> NextH(unsigned int i, unsigned int j, double gP);
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

    QSurfaceDataProxy* m_sqrtSinProxy;
    QSurface3DSeries* m_sqrtSinSeries;

    Solver* TestSolver;

//    void Calc(unsigned int StepsP);
};
//-----------------------------//
#endif
