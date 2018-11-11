#ifndef MAINWINDOW_H
#define MAINWINDOW_H
//-----------------------------//
#include <QMainWindow>
#include <QOpenGLWidget>
#include <QLineEdit>
#include <QIntValidator>
#include <QPushButton>
#include <QWidget>
#include <QGridLayout>
#include <QSpacerItem>
//-----------------------------//
#include "Include/DivisionGL.h"
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

    QLineEdit* DivivsionLine;
    QIntValidator* DivisionValidator;
    QPushButton* DivisionButton;
    QSpacerItem* DivideSpacer;

    DivisionGL* DivisionDisplay;

    //----------//

    void InitMain();
};
//-----------------------------//
#endif
