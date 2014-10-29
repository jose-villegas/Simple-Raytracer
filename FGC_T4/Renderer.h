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
        GLuint frameBuffer;
        GLuint renderedTexture;
        Shader rt;
        GLenum drawBuffer;
        void loadGLSLShader();

        void logShaderCompilerErrors();

        void createRenderTarget();
        void createRenderTexture();
        void configureFrameBuffer();
        static void keyCallback(GLFWwindow *window, int key, int scancode, int action, int mods);
        void setupRenderQuad();
    public:
        void init();
        void terminate();
        Renderer(void);
        ~Renderer(void);
};

