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
        GLuint VAO;		// Main Vertex Array Object
        GLuint VBO;		// Main Vertex Buffer Object
        GLuint FBO;		// Main Frame Buffer Object
        GLuint CBF;		// Color Buffer
        GLuint RTX;		// Render Texture
        GLuint MFBO;	// Multi-Sample FBO
        GLuint MTX;		// Multi-Sample Texture
        GLuint MCBF;	// Multi-Sample Color Buffer
        GLuint DBF;		// Draw Buffer
        void setupRenderQuad();
        void loadGLSLShader();
        void createRenderTarget(int width, int height);
        void logShaderCompilerErrors();
        static void keyCallback(GLFWwindow *window, int key, int scancode, int action, int mods);
    public:
        void init();
        void terminate();
        void toFile(std::string filename, int renderWidth, int renderHeight);
        Renderer(void);
        ~Renderer(void);
};

