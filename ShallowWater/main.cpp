#include <iostream>
#include <vector>
//-----------------------------//
#include <QApplication>
#include <QObject>
//-----------------------------//
#include "UI/MainWindow.h"
#include "Calculation/Calculation.h"
//-----------------------------//
int main(int argc, char** argv) {
    QApplication App(argc, argv);

    MainWindow Window;
    Window.show();

    Calculation Calc(50, 50);

    QObject::connect(&Calc, &Calculation::SendDisplayDataSignal, &Window, &MainWindow::GetDisplayDataSlot);

    Calc.Solve(0);
    Calc.DisplayData();

    return QApplication::exec();
}