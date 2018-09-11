#include <iostream>
#include <vector>
//-----------------------------//
#include <QApplication>
//-----------------------------//
#include "MainWindow.h"
//-----------------------------//
#include "Include/dMath/dMath.h"
#include "Include/Calculation/Calculation.h"
//-----------------------------//
int main(int argc, char** argv) {
    QApplication App(argc, argv);

    MainWindow Window;
    Window.show();

    return QApplication::exec();
}