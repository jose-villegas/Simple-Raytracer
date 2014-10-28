#include "Renderer.h"

Renderer::Renderer(void)
{
}


Renderer::~Renderer(void)
{
}

void Renderer::createRenderTarget()
{
    frameBuffer = 0;
    glGenFramebuffers(1, &frameBuffer);
    glBindFramebuffer(GL_FRAMEBUFFER, frameBuffer);
}

void Renderer::createRenderTexture()
{
    glGenTextures(1, &renderedTexture);
    glBindTexture(GL_TEXTURE_2D, renderedTexture);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 1024, 768, 0, GL_RGB, GL_UNSIGNED_BYTE, 0);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
}

void Renderer::configureFrameBuffer()
{
    glFramebufferTexture(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, renderedTexture, 0);
    drawBuffer = GL_COLOR_ATTACHMENT0;
    glDrawBuffers(1, &drawBuffer);
}

void Renderer::init()
{
    if (!glfwInit()) {
        exit(EXIT_FAILURE);
    }

    GLFWwindow *window = glfwCreateWindow(640, 480, "Raytracing - RT", NULL, NULL);

    if (!window) {
        glfwTerminate();
        exit(EXIT_FAILURE);
    }

    glfwMakeContextCurrent(window);
    glewExperimental = GL_FALSE;
    GLenum error = glGetError();

    if (error != GL_NO_ERROR) {
        std::cout << "OpenGL Error: " << error << std::endl;
    }

    GLenum glewinit = glewInit();

    if (glewinit != GLEW_OK) {
        std::cout << "Glew not okay! " << glewinit;
        exit(EXIT_FAILURE);
    }

    glfwSetKeyCallback(window, key_callback);

    while (!glfwWindowShouldClose(window)) {
        glClear(GL_COLOR_BUFFER_BIT);
        // Render
        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    glfwDestroyWindow(window);
    glfwTerminate();
}

void Renderer::terminate()
{
    glfwTerminate();
}

void Renderer::key_callback(GLFWwindow *window, int key, int scancode, int action, int mods)
{
    if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS) {
        glfwSetWindowShouldClose(window, GL_TRUE);
    }
}

void Renderer::loadGLSLShader()
{
}
