#ifndef MAINWINDOW_H
#define MAINWINDOW_H
//-----------------------------//
#include <iostream>
//-----------------------------//
#include <QMainWindow>
#include <QWidget>
#include <QLineEdit>
#include <QLabel>
#include <QGridLayout>
#include <QRegExpValidator>
#include <QPushButton>
#include <QtDataVisualization>
#include <QSpacerItem>
#include <QFrame>
#include <QObject>
#include <QVector>
//-----------------------------//
class MainWindow : public QMainWindow {
    Q_OBJECT
public:
    explicit MainWindow(QWidget* ParentP = nullptr);
    ~MainWindow() override = default;
private:
    QWidget* MainWidget;
    QGridLayout* MainLayout;

    //----------//

    QLabel* xGridLabel;
    QLabel* yGridLabel;

    QLabel* TimeStepLabel;
    QLabel* xStepLabel;
    QLabel* yStepLabel;

    QLineEdit* xGridLine;
    QLineEdit* yGridLine;

    QLineEdit* TimeStepLine;
    QLineEdit* xStepLine;
    QLineEdit* yStepLine;

    QRegExpValidator* IntValidator;
    QRegExpValidator* DoubleValidator;

    QPushButton* CalculateButton;

    QSpacerItem* ControlSpacer;

    //----------//

    QtDataVisualization::Q3DSurface* Plot;

    QtDataVisualization::QSurfaceDataArray* DisplayData;
    QVector <QtDataVisualization::QSurfaceDataRow*> DisplayDataRows;
    QtDataVisualization::QSurface3DSeries* DataSeries;

    //----------//

    void InitMain();
    void InitValidators();
    void InitControl();
    void InitSurface();
    void InitConnections();
public slots:
    void GetDisplayDataSlot(const QVector <QVector <double>>& DataP, double xGridStepP, double yGridStepP);
private slots:
    void EnableCalculateButtonSlot();
};
//-----------------------------//
#endif
