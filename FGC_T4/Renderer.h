#pragma once
#define NOMINMAX
#include <GL/glew.h>
#include "GLFW/glfw3.h"
#include <string>
#include <fstream>
#include <iostream>
#include <unordered_map>

class Shader {
    public:
        GLuint fShader;
        GLuint vShader;
        GLuint program;
        std::string source;
        std::unordered_map<std::string, GLuint> uniformLocs;
};

class Renderer {
    private:
        Shader rt;
        GLuint VAO;
        GLuint VBO;
        void setupRenderQuad();
        void loadGLSLShader();
        void logShaderCompilerErrors();
        static void keyCallback(GLFWwindow *window, int key, int scancode, int action, int mods);
    public:
        void init();
        void terminate();
        Renderer(void);
        ~Renderer(void);
};

