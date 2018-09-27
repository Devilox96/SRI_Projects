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
struct Point {
    Point(double hP, double uP, double vP) : h(hP), u(uP), v(vP) {}

    double h, u, v;
};
//-----------------------------//
void InitGrid(std :: vector <std :: vector <dVector3D <double>>>& GridP, double ExcitationP, double VXP, double VYP, int PointsNumP);

dVector3D <double> GetH(const std :: vector <std :: vector <dVector3D <double>>>& GridP, unsigned int i, unsigned int j);
dVector3D <double> GetU(const dVector3D <double>& HP, double gP);
dVector3D <double> GetV(const dVector3D <double>& HP, double gP);

dVector3D <double> GetHiHalf(double DeltaTP, double DeltaXP, const std :: vector <std :: vector <dVector3D <double>>>& GridP, unsigned int i, unsigned int j, double gP);
dVector3D <double> GetHjHalf(double DeltaTP, double DeltaYP, const std :: vector <std :: vector <dVector3D <double>>>& GridP, unsigned int i, unsigned int j, double gP);

dVector3D <double> NextH(double DeltaTP, double DeltaXP, double DeltaYP, const std :: vector <std :: vector <dVector3D <double>>>& GridP, unsigned int i, unsigned int j, double gP);
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

    std :: vector <std :: vector <dVector3D <double>>> Grid;

    void Calc(unsigned int StepsP);
};
//-----------------------------//
#endif
