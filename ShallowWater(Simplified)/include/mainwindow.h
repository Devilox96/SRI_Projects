#ifndef MAINWINDOW_H
#define MAINWINDOW_H
//-----------------------------//
#include <iostream>
#include <chrono>
//-----------------------------//
#include <QMainWindow>
#include <QGridLayout>
#include <QtDataVisualization/Q3DSurface>
#include <QtDataVisualization/QSurfaceDataProxy>
#include <QtDataVisualization/QHeightMapSurfaceDataProxy>
#include <QtDataVisualization/QSurface3DSeries>
#include <QtMath>
//-----------------------------//
//#include "dvectors.h"
#include "include/dRichtmyerMethod.h"
#include "Libs/dMath/Core/dVectors.h"
//-----------------------------//
using namespace QtDataVisualization;
//-----------------------------//
class Solver : public dRichtmyerMethod {
public:
    dVectorND <double> xFunc(const dVectorND <double>& VecP) override {
        return dVectorND <double> ({VecP[1],
                                    pow(VecP[1], 2.0) / VecP[0] + 0.5 * 9.81 * pow(VecP[0], 2.0),
                                    VecP[1] * VecP[2] / VecP[0]});
    }
    dVectorND <double> yFunc(const dVectorND <double>& VecP) override {
        return dVectorND <double> ({VecP[2],
                                    VecP[1] * VecP[2] / VecP[0],
                                    pow(VecP[2], 2.0) / VecP[0] + 0.5 * 9.81 * pow(VecP[0], 2.0)});
    }
    dVectorND <double> zFunc(const dVectorND <double>& VecP) override {
        return dVectorND <double>();
    }
};
//-----------------------------//
void InitGrid(std :: vector <std :: vector <dVectorND <double>>>& GridP, double ExcitationP, double VXP, double VYP, int PointsNumP);

dVectorND <double> GetH(const std :: vector <std :: vector <dVectorND <double>>>& GridP, unsigned int i, unsigned int j);
dVectorND <double> GetU(const dVectorND <double>& HP, double gP);
dVectorND <double> GetV(const dVectorND <double>& HP, double gP);

dVectorND <double> GetHiHalf(double DeltaTP, double DeltaXP, const std :: vector <std :: vector <dVectorND <double>>>& GridP, unsigned int i, unsigned int j, double gP);
dVectorND <double> GetHjHalf(double DeltaTP, double DeltaYP, const std :: vector <std :: vector <dVectorND <double>>>& GridP, unsigned int i, unsigned int j, double gP);

dVectorND <double> NextH(double DeltaTP, double DeltaXP, double DeltaYP, const std :: vector <std :: vector <dVectorND <double>>>& GridP, unsigned int i, unsigned int j, double gP);
//-----------------------------//
class MainWindow : public QMainWindow {
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = nullptr);
    ~MainWindow() override = default;

    Q3DSurface *graph;
    QWidget* MainWidget;
    QGridLayout* MainLayout;

    std :: vector <std :: vector <dVectorND <double>>> Grid;

    void Calc(unsigned int StepsP);
};
//-----------------------------//
#endif
