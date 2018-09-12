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
#include "dvectors.h"
//-----------------------------//
using namespace QtDataVisualization;
//-----------------------------//
struct Point {
    Point(double hP, double uP, double vP) : h(hP), u(uP), v(vP) {}

    double h, u, v;
};
//-----------------------------//
void InitGrid(std :: vector <std :: vector <Point>>& GridP, double ExcitationP, double VXP, double VYP, int PointsNumP);

dVector3D GetH(const std :: vector <std :: vector <Point>>& GridP, unsigned int i, unsigned int j);
dVector3D GetU(const dVector3D& HP, double gP);
dVector3D GetV(const dVector3D& HP, double gP);

dVector3D GetHiHalf(double DeltaTP, double DeltaXP, const std :: vector <std :: vector <Point>>& GridP, unsigned int i, unsigned int j, double gP);
dVector3D GetHjHalf(double DeltaTP, double DeltaYP, const std :: vector <std :: vector <Point>>& GridP, unsigned int i, unsigned int j, double gP);

dVector3D NextH(double DeltaTP, double DeltaXP, double DeltaYP, const std :: vector <std :: vector <Point>>& GridP, unsigned int i, unsigned int j, double gP);
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

    std :: vector <std :: vector <Point>> Grid;

    void Calc(unsigned int StepsP);
};
//-----------------------------//
#endif
