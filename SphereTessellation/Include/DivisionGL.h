#ifndef DIVISIONGL_H
#define DIVISIONGL_H
//-----------------------------//
#include <QOpenGLWidget>
#include <QOpenGLShaderProgram>
#include <QOpenGLBuffer>
#include <QOpenGLFunctions>
#include <QOpenGLVertexArrayObject>
#include <QVector>
//-----------------------------//
class DivisionGL : public QOpenGLWidget, protected QOpenGLFunctions {
    Q_OBJECT
public:
    explicit DivisionGL(QWidget* ParentP = nullptr);
    ~DivisionGL() override = default;
protected:
    void initializeGL() override;
    void paintGL() override;
    void resizeGL(int WidthP, int HeightP) override;
private:
    QOpenGLShaderProgram ShaderProgram;
    QOpenGLBuffer CoordBuffer;
    QOpenGLBuffer ColorBuffer;
    QOpenGLVertexArrayObject VAO;

    int u_modelToWorld;
    int u_worldToView;
    QMatrix4x4 m_projection;
    QMatrix4x4 matrix;

    //----------//

    void MakeShader();
    void MakeTriangle();
};
//-----------------------------//
#endif
