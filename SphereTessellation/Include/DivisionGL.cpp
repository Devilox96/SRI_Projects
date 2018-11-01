#include "DivisionGL.h"
//-----------------------------//
DivisionGL::DivisionGL(QWidget* ParentP) : QOpenGLWidget(ParentP) {}
//-----------------------------//
void DivisionGL::initializeGL() {
    initializeOpenGLFunctions();
//    connect(this, SIGNAL(frameSwapped()), this, SLOT(update()));

//    glEnable(GL_CULL_FACE);
    glEnable(GL_DEPTH_TEST);
    glClearColor(0.1, 0.1, 0.2, 1.0);

    MakeShader();
    MakeTriangle();

    matrix.ortho(-2.0f, 2.0, -2.0f, 2.0, 2.0, -2.0f);
    matrix.translate(0.0, 0.0, 1.0);
}
void DivisionGL::paintGL() {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

//    QMatrix4x4 matrix;


    ShaderProgram.bind();

    ShaderProgram.setUniformValue(u_modelToWorld, matrix);


    VAO.bind();
    glDrawArrays(GL_TRIANGLES, 0, 3);
    VAO.release();

    ShaderProgram.release();
}
void DivisionGL::resizeGL(int WidthP, int HeightP) {
    int side = qMin(WidthP, HeightP);
    glViewport((WidthP - side) / 2, (HeightP - side) / 2, side, side);
    m_projection.setToIdentity();
    m_projection.perspective(45.0f, WidthP / float(HeightP), 0.0f, 1000.0f);
}
//-----------------------------//
void DivisionGL::MakeShader() {
    if (!ShaderProgram.addShaderFromSourceFile(QOpenGLShader::Vertex, "shaders/vShader.glsl")) {
        qCritical() << QObject::tr("Could not compile vertex shader: ") << ShaderProgram.log();
    }
    if (!ShaderProgram.addShaderFromSourceFile(QOpenGLShader::Fragment, "shaders/fShader.glsl")) {
        qCritical() << QObject::tr("Could not compile fragment shader: ") << ShaderProgram.log();
    }
    if (!ShaderProgram.link()) {
        qCritical() << QObject::tr("Could not link shader program: ") << ShaderProgram.log();
    }

    ShaderProgram.bind();
}
void DivisionGL::MakeTriangle() {
    QVector <float> VerticesL = {
           -1.0,   -1.0,    0.0,
            1.0,   -1.0,    0.0,
            0.0,    1.0,    0.0
    };

    QVector <float> ColorsL = {
            1.0,    0.0,    0.0,
            0.0,    1.0,    0.0,
            0.0,    0.0,    1.0
    };

    u_modelToWorld = ShaderProgram.uniformLocation("modelToWorld");
    u_worldToView = ShaderProgram.uniformLocation("worldToView");

    VAO.create();
    VAO.bind();

    CoordBuffer.create();
    CoordBuffer.bind();
    CoordBuffer.allocate(VerticesL.constData(), VerticesL.count() * sizeof(float));
    ShaderProgram.enableAttributeArray(0);
    ShaderProgram.setAttributeBuffer(0, GL_FLOAT, 0, 3);

    ColorBuffer.create();
    ColorBuffer.bind();
    ColorBuffer.allocate(ColorsL.constData(), ColorsL.count() * sizeof(float));
    ShaderProgram.enableAttributeArray(1);
    ShaderProgram.setAttributeBuffer(1, GL_FLOAT, 0, 3);

    VAO.release();
    CoordBuffer.release();
    ColorBuffer.release();
    ShaderProgram.release();
}