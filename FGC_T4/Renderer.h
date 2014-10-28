#pragma once
#define NOMINMAX
#include <GL/glew.h>
#include "GLFW/glfw3.h"
#include <string>
#include <iostream>

class Renderer {
    private:
        GLuint frameBuffer;
        GLuint renderedTexture;
        GLenum drawBuffer;
        std::string shaderSource;
        void loadGLSLShader();
        void createRenderTarget();
        void createRenderTexture();
        void configureFrameBuffer();
        static void key_callback(GLFWwindow *window, int key, int scancode, int action, int mods);
    public:
        void init();
        void terminate();
        Renderer(void);
        ~Renderer(void);
};

