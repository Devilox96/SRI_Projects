#version 450

layout(location = 0) in vec3 Pos;
layout(location = 1) in vec3 Col;

layout(location = 0) out vec3 FragCol;

void main() {
    gl_Position = vec4(Pos, 1.0);
    FragCol = Col;
}