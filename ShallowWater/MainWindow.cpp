#include "MainWindow.h"
//-----------------------------//
MainWindow::MainWindow(QWidget*) {
    InitMain();
    InitValidators();
    InitControl();
    InitSurface();
    InitConnections();
}
//-----------------------------//
void MainWindow::InitMain() {
    MainWidget = new QWidget;
    MainLayout = new QGridLayout;

    setCentralWidget(MainWidget);
    MainWidget -> setLayout(MainLayout);
}
void MainWindow::InitValidators() {
    IntValidator = new QRegExpValidator(QRegExp("[1-9][0-9]{0,4}"));
    DoubleValidator = new QRegExpValidator(QRegExp("[0-1](\\.[0-9]{1,8})?"));
}
void MainWindow::InitControl() {
    xGridLabel = new QLabel("Grid size X:");
    yGridLabel = new QLabel("Grid size Y:");

    TimeStepLabel = new QLabel("Time step:");
    xStepLabel = new QLabel("X coord step:");
    yStepLabel = new QLabel("Y coord step:");

    xGridLine = new QLineEdit("50");
    yGridLine = new QLineEdit("50");

    TimeStepLine = new QLineEdit("0.005");
    xStepLine = new QLineEdit("0.1");
    yStepLine = new QLineEdit("0.1");

    xGridLine -> setValidator(IntValidator);
    yGridLine -> setValidator(IntValidator);

    TimeStepLine -> setValidator(DoubleValidator);
    xStepLine -> setValidator(DoubleValidator);
    yStepLine -> setValidator(DoubleValidator);

    CalculateButton = new QPushButton("Calculate");
    CalculateButton -> setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);

    ControlSpacer = new QSpacerItem(0, 0, QSizePolicy::Minimum, QSizePolicy::Expanding);

    MainLayout -> addWidget(xGridLabel, 0, 0, 1, 1, Qt::AlignRight);
    MainLayout -> addWidget(xGridLine, 0, 1, 1, 1, Qt::AlignRight);

    MainLayout -> addWidget(yGridLabel, 1, 0, 1, 1, Qt::AlignRight);
    MainLayout -> addWidget(yGridLine, 1, 1, 1, 1, Qt::AlignRight);

    MainLayout -> addWidget(TimeStepLabel, 2, 0, 1, 1, Qt::AlignRight);
    MainLayout -> addWidget(TimeStepLine, 2, 1, 1, 1, Qt::AlignRight);

    MainLayout -> addWidget(xStepLabel, 3, 0, 1, 1, Qt::AlignRight);
    MainLayout -> addWidget(xStepLine, 3, 1, 1, 1, Qt::AlignRight);

    MainLayout -> addWidget(yStepLabel, 4, 0, 1, 1, Qt::AlignRight);
    MainLayout -> addWidget(yStepLine, 4, 1, 1, 1, Qt::AlignRight);

    MainLayout -> addWidget(CalculateButton, 5, 0, 1, 2, Qt::AlignHCenter);

    MainLayout -> addItem(ControlSpacer, 6, 0, 1, 2);
}
void MainWindow::InitSurface() {
    using namespace QtDataVisualization;

    auto surface = new Q3DSurface;
    surface -> setFlags(Qt::FramelessWindowHint);
    auto data = new QSurfaceDataArray;
    QSurfaceDataRow *dataRow1 = new QSurfaceDataRow;
    auto dataRow2 = new QSurfaceDataRow;

    *dataRow1 << QVector3D(0.0f, 0.1f, 0.5f) << QVector3D(1.0f, 0.5f, 0.5f);
    *dataRow2 << QVector3D(0.0f, 1.8f, 1.0f) << QVector3D(1.0f, 1.2f, 1.0f);
    *data << dataRow1 << dataRow2;

    auto series = new QSurface3DSeries;
    series->dataProxy()->resetArray(data);
    surface -> addSeries(series);

    auto SurfaceFrame = new QFrame;
    SurfaceFrame -> setFrameShape(QFrame::Shape::Box);
    SurfaceFrame -> setFrameShadow(QFrame::Shadow::Raised);
    SurfaceFrame -> setLineWidth(0);
    SurfaceFrame -> setMidLineWidth(1);

    auto PlotLayout = new QGridLayout;
    PlotLayout -> setMargin(0);
    SurfaceFrame -> setLayout(PlotLayout);

    auto Container = QWidget::createWindowContainer(surface);
    Container -> setMinimumSize(400, 400);
    Container -> setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
    PlotLayout -> addWidget(Container);

    MainLayout -> addWidget(SurfaceFrame, 0, 2, 7, 1);
}
void MainWindow::InitConnections() {
    connect(xGridLine, &QLineEdit::textChanged, this, &MainWindow::EnableCalculateButtonSlot);
    connect(yGridLine, &QLineEdit::textChanged, this, &MainWindow::EnableCalculateButtonSlot);

    connect(TimeStepLine, &QLineEdit::textChanged, this, &MainWindow::EnableCalculateButtonSlot);
    connect(xStepLine, &QLineEdit::textChanged, this, &MainWindow::EnableCalculateButtonSlot);
    connect(yStepLine, &QLineEdit::textChanged, this, &MainWindow::EnableCalculateButtonSlot);
}
//-----------------------------//
void MainWindow::EnableCalculateButtonSlot() {
    if (TimeStepLine -> text().toDouble()   == 0.0      ||
        xStepLine -> text().toDouble()      == 0.0      ||
        yStepLine -> text().toDouble()      == 0.0      ||
        xGridLine -> text().toInt()         < 10        ||
        yGridLine -> text().toInt()         < 10) {
        CalculateButton -> setDisabled(true);
    } else {
        CalculateButton -> setEnabled(true);
    }
}