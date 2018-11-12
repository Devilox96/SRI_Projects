#ifndef DIVISIONGL_H
#define DIVISIONGL_H
//-----------------------------//
#include <iostream>
//-----------------------------//
#include <QOpenGLWidget>
#include <QOpenGLShaderProgram>
#include <QOpenGLBuffer>
#include <QOpenGLFunctions_4_5_Core>
#include <QOpenGLVertexArrayObject>
#include <QVector>
#include <QOpenGLContext>
#include <QSurfaceFormat>
//-----------------------------//
class DivisionGL : public QOpenGLWidget, protected QOpenGLFunctions_4_5_Core {
    Q_OBJECT
public:
    explicit DivisionGL(QWidget* ParentP = nullptr);
    ~DivisionGL() override = default;
protected:
    void initializeGL() override;
    void paintGL() override;
    void resizeGL(int WidthP, int HeightP) override;
private:
    QMatrix4x4 ModelMatrix;
    QMatrix4x4 ViewMatrix;
    QMatrix4x4 ProjectionMatrix;

    GLuint ShaderProgram;
    GLuint VertexShader;
    GLuint FragmentShader;

    GLuint VAO;
    GLuint VBO;
    GLuint EBO;

    //----------//

    QVector <GLfloat> vertices;
    QVector <GLuint> indices;

    //----------//

    void SetMatrices();
    void MakeShader();
    void MakeTriangle();
};
//-----------------------------//
#endif
