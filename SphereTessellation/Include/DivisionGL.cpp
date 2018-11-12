#include "DivisionGL.h"
//-----------------------------//
DivisionGL::DivisionGL(QWidget* ParentP) : QOpenGLWidget(ParentP) {
    vertices = {
           -0.5f,  -0.5f,   0.0f,
            1.0f,   0.0f,   0.0f,
            0.5f,  -0.5f,   0.0f,
            0.0f,   1.0f,   0.0f,
            0.0f,   0.5f,   0.0f,
            0.0f,   0.0f,   1.0f,
            0.0f,  -1.5f,   0.0f,
            0.3f,   0.4f,   0.5f,
    };
    indices = {
            0,  1,  2,
            0,  1,  3
    };

    InitMatrices();
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

    glBindVertexArray(VAO);
    glDrawElements(GL_TRIANGLES, indices.size(), GL_UNSIGNED_INT, 0);
    glBindVertexArray(0);
}
void DivisionGL::resizeGL(int WidthP, int HeightP) {
    int side = qMin(WidthP, HeightP);
    glViewport((WidthP - side) / 2, (HeightP - side) / 2, side, side);
    ProjectionMatrix.setToIdentity();
    ProjectionMatrix.perspective(45.0f, WidthP / float(HeightP), 0.0f, 1000.0f);
}
//-----------------------------//
void DivisionGL::InitMatrices() {
    ModelMatrix.setToIdentity();
    ViewMatrix.translate(0.0, 0.0, -3.0);
}
void DivisionGL::MakeShader() {
    const char* VertexShaderSource =
            "#version 450\n"
            "layout(location = 0) in vec3 position;\n"
            "layout(location = 1) in vec3 color;\n"
            "out vec4 vColor;\n"
            "uniform mat4 modelToWorld;\n"
            "uniform mat4 worldToView;\n"
            "void main() {\n"
            "gl_Position = vec4(position, 1.0);\n"
            "vColor = vec4(color, 1.0);\n"
            "}\0";

    const char* FragmentShaderSource =
            "#version 450\n"
            "in vec4 vColor;\n"
            "out vec4 fColor;\n"
            "void main() {\n"
            "fColor = vColor;\n"
            "}\0";

    VertexShader = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(VertexShader, 1, &VertexShaderSource, nullptr);
    glCompileShader(VertexShader);

    int success;
    char infoLog[512];
    glGetShaderiv(VertexShader, GL_COMPILE_STATUS, &success);
    if (!success)
    {
        glGetShaderInfoLog(VertexShader, 512, NULL, infoLog);
        std::cout << "ERROR::SHADER::VERTEX::COMPILATION_FAILED\n" << infoLog << std::endl;
    } else {
        std::cout << "Success!" << std::endl;
    }

    FragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(FragmentShader, 1, &FragmentShaderSource, nullptr);
    glCompileShader(FragmentShader);

    glGetShaderiv(FragmentShader, GL_COMPILE_STATUS, &success);
    if (!success)
    {
        glGetShaderInfoLog(FragmentShader, 512, NULL, infoLog);
        std::cout << "ERROR::SHADER::FRAGMENT::COMPILATION_FAILED\n" << infoLog << std::endl;
    } else {
        std::cout << "Success!" << std::endl;
    }

    ShaderProgram = glCreateProgram();
    glAttachShader(ShaderProgram, VertexShader);
    glAttachShader(ShaderProgram, FragmentShader);
    glLinkProgram(ShaderProgram);

    glGetProgramiv(ShaderProgram, GL_LINK_STATUS, &success);
    if (!success) {
        glGetProgramInfoLog(ShaderProgram, 512, NULL, infoLog);
        std::cout << "ERROR::SHADER::PROGRAM::LINKING_FAILED\n" << infoLog << std::endl;
    } else {
        std::cout << "Success!" << std::endl;
    }

    glDeleteShader(VertexShader);
    glDeleteShader(FragmentShader);
}
void DivisionGL::MakeTriangle() {
    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);
    glGenBuffers(1, &EBO);

    glBindVertexArray(VAO);

    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(float), vertices.constData(), GL_STATIC_DRAW);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(unsigned int), indices.constData(), GL_STATIC_DRAW);

    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), nullptr);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)(3 * sizeof(float)));
    glEnableVertexAttribArray(1);

    glBindVertexArray(0);

    glUseProgram(ShaderProgram);
}