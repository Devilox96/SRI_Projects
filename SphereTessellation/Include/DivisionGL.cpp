#include "DivisionGL.h"
//-----------------------------//
DivisionGL::DivisionGL(QWidget* ParentP) : QOpenGLWidget(ParentP) {
//    vertices = {
//           -0.5f,  -1.0f,   0.0f,
//            1.0f,   0.0f,   0.0f,
//            0.5f,  -0.5f,   0.0f,
//            0.0f,   1.0f,   0.0f,
//            0.0f,   0.5f,   0.0f,
//            0.0f,   0.0f,   1.0f,
//            0.0f,  -1.5f,   0.0f,
//            0.3f,   0.4f,   0.5f,
//    };
//    indices = {
//            0,  1,  2,
//            0,  1,  3
//    };
    vertices = {
            -1.0f,   1.0f,   -1.0f, //---0---//
            1.0f,   0.0f,    0.0f,
            1.0f,   1.0f,   -1.0f, //---1---//
            0.0f,   1.0f,    0.0f,
            1.0f,   1.0f,    1.0f, //---2---//
            0.0f,   0.0f,    1.0f,
            -1.0f,   1.0f,    1.0f, //---3---//
            1.0f,   1.0f,    0.0f,
            -1.0f,  -1.0f,   -1.0f, //---4---//
            1.0f,   0.0f,    1.0f,
            1.0f,  -1.0f,   -1.0f, //---5---//
            0.0f,   1.0f,    1.0f,
            1.0f,  -1.0f,    1.0f, //---6---//
            1.0f,   1.0f,    1.0f,
            -1.0f,  -1.0f,    1.0f  //---7---//
    };
    indices = {
            0,  1,  2,  3,
            7,  6,  5,  4,
            0,  4,  5,  1,
            3,  2,  6,  7,
            1,  5,  6,  2,
            3,  7,  4,  0
    };

}
//-----------------------------//
void DivisionGL::initializeGL() {
    initializeOpenGLFunctions();
    glEnable(GL_DEPTH_TEST);
    glClearColor(0.1, 0.1, 0.2, 1.0);

    MakeShader();
    MakeTriangle();
}
void DivisionGL::paintGL() {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glEnable(GL_DEPTH_TEST);

    SetMatrices();

    glBindVertexArray(BoxVAO);
    glDrawElements(GL_QUADS, 24, GL_UNSIGNED_INT, 0);
    glBindVertexArray(0);
}
void DivisionGL::resizeGL(int WidthP, int HeightP) {
    int side = qMin(WidthP, HeightP);
    glViewport((WidthP - side) / 2, (HeightP - side) / 2, side, side);
}
//-----------------------------//
void DivisionGL::SetMatrices() {
    ModelMatrix.setToIdentity();
    ViewMatrix.setToIdentity();
    ViewMatrix.translate(0.0, 0.5f, -3.0f);
    ProjectionMatrix.setToIdentity();
    ProjectionMatrix.perspective(45.0f, float(width()) / float(height()), 0.1f, 1000.0f);

    GLint ModelLocationL = glGetUniformLocation(ShaderProgram, "model");
    GLint ViewLocationL = glGetUniformLocation(ShaderProgram, "view");
    GLint ProjectionLocationL = glGetUniformLocation(ShaderProgram, "projection");

    glBindAttribLocation(ShaderProgram, 0, "position");
    glBindAttribLocation(ShaderProgram, 1, "color");

    glLinkProgram(ShaderProgram);

    glUniformMatrix4fv(ModelLocationL, 1, GL_FALSE, ModelMatrix.constData());
    glUniformMatrix4fv(ViewLocationL, 1, GL_FALSE, ViewMatrix.constData());
    glUniformMatrix4fv(ProjectionLocationL, 1, GL_FALSE, ProjectionMatrix.constData());
}
void DivisionGL::MakeShader() {
    std::ifstream vShader("../Shaders/vShader.glsl");
    std::string vShaderStr((std::istreambuf_iterator <char>(vShader)), std::istreambuf_iterator <char>());
    const char* vShaderData = vShaderStr.data();

    std::ifstream fShader("../Shaders/fShader.glsl");
    std::string fShaderStr((std::istreambuf_iterator <char>(fShader)), std::istreambuf_iterator <char>());
    const char* fShaderData = fShaderStr.data();

    int SuccessL;
    char ErrorLogL[512];

    //---Vertex shader---//
    GLuint vShaderInstance = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vShaderInstance, 1, &vShaderData, nullptr);
    glCompileShader(vShaderInstance);

    glGetShaderiv(vShaderInstance, GL_COMPILE_STATUS, &SuccessL);

    if (SuccessL) {
        std::cout << "Vertex shader has been compiled successfully!" << std::endl;
    } else {
        glGetShaderInfoLog(vShaderInstance, 512, nullptr, ErrorLogL);
        std::cout << "Vertex shader compilation error: " << ErrorLogL << std::endl;
    }
    //---Vertex shader---//

    //---Fragment shader---//
    GLuint fShaderInstance = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fShaderInstance, 1, &fShaderData, nullptr);
    glCompileShader(fShaderInstance);

    glGetShaderiv(fShaderInstance, GL_COMPILE_STATUS, &SuccessL);

    if (SuccessL) {
        std::cout << "Fragment shader has been compiled successfully!" << std::endl;
    } else {
        glGetShaderInfoLog(fShaderInstance, 512, nullptr, ErrorLogL);
        std::cout << "Fragment shader compilation error: " << ErrorLogL << std::endl;
    }
    //---Fragment shader---//

    //---Shader program---//
    ShaderProgram = glCreateProgram();

    glAttachShader(ShaderProgram, vShaderInstance);
    glAttachShader(ShaderProgram, fShaderInstance);

    glLinkProgram(ShaderProgram);

    glGetProgramiv(ShaderProgram, GL_LINK_STATUS, &SuccessL);

    if (SuccessL) {
        std::cout << "Shader program has been linked successfully!" << std::endl;
    } else {
        glGetProgramInfoLog(ShaderProgram, 512, nullptr, ErrorLogL);
        std::cout << "Shader program linking error: " << ErrorLogL << std::endl;
    }
    //---Shader program---//

    glDeleteShader(vShaderInstance);
    glDeleteShader(fShaderInstance);
}
void DivisionGL::MakeTriangle() {
    GLuint VBO;
    GLuint EBO;

    glGenVertexArrays(1, &BoxVAO);
    glGenBuffers(1, &VBO);
    glGenBuffers(1, &EBO);

    glBindVertexArray(BoxVAO);

    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, 45 * sizeof(float), vertices.constData(), GL_STATIC_DRAW);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, 24 * sizeof(unsigned int), indices.constData(), GL_STATIC_DRAW);

    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), nullptr);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)(3 * sizeof(float)));
    glEnableVertexAttribArray(1);

    glBindVertexArray(0);

    glUseProgram(ShaderProgram);
}