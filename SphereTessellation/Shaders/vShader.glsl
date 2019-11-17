#version 130
//-----------------------------//
in vec3 position;
in vec3 color;
out vec4 vColor;
uniform mat4 model;
uniform mat4 view;
uniform mat4 projection;
//-----------------------------//
void main() {
    gl_Position = projection * view * model * vec4(position, 1.0);
    vColor = vec4(color, 1.0);
}