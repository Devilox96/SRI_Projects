#include "MainWindow.h"
//-----------------------------//
MainWindow::MainWindow(QWidget* ParentP) : QMainWindow(ParentP) {
    InitMain();

    DivivsionLine = new QLineEdit;
    DivivsionLine -> setFixedSize(100, 25);
    DivisionValidator = new QIntValidator(1, 5);
    DivivsionLine -> setValidator(DivisionValidator);

    DivisionButton = new QPushButton("Divide");
    DivisionButton -> setFixedSize(100, 25);

    DivideSpacer = new QSpacerItem(0, 0, QSizePolicy::Minimum, QSizePolicy::Expanding);

    DivisionDisplay = new DivisionGL;
    DivisionDisplay -> setMinimumSize(800, 600);

    MainLayout -> addWidget(DivivsionLine, 0, 0, 1, 1);
    MainLayout -> addWidget(DivisionButton, 1, 0, 1, 1);
    MainLayout -> addItem(DivideSpacer, 2, 0, 1, 1);

    MainLayout -> addWidget(DivisionDisplay, 0, 1, 3, 1);
}
//-----------------------------//
void MainWindow::InitMain() {
    MainWidget = new QWidget;
    setCentralWidget(MainWidget);

    MainLayout = new QGridLayout;
    MainWidget -> setLayout(MainLayout);
}