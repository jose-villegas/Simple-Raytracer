#include "Renderer.h"
#include "glm.hpp"

Renderer::Renderer(void)
{
}


Renderer::~Renderer(void)
{
}

void Renderer::init()
{
    if (!glfwInit()) {
        exit(EXIT_FAILURE);
    }

    glfwWindowHint(GLFW_SAMPLES, 16); // Anti Alias
    GLFWwindow *window = glfwCreateWindow(1440, 900, "Raytracing - RT", NULL, NULL);

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

    glfwSetKeyCallback(window, keyCallback);
    loadGLSLShader();
    setupRenderQuad();
    // Load Relevant Uniforms
    glUseProgram(rt.program);
    glUniform3fv(rt.uniformLocs["iResolution"], 1, &glm::vec3(1440, 900, 0)[0]);

    while (!glfwWindowShouldClose(window)) {
        glClear(GL_COLOR_BUFFER_BIT);
        glDrawArrays(GL_TRIANGLES, 0, 6);
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

void Renderer::keyCallback(GLFWwindow *window, int key, int scancode, int action, int mods)
{
    if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS) {
        glfwSetWindowShouldClose(window, GL_TRUE);
    }
}

void Renderer::loadGLSLShader()
{
    std::ifstream ifs;
    ifs.open("raytracing_fragment.glsl");
    rt.source.assign((std::istreambuf_iterator<char>(ifs)), (std::istreambuf_iterator<char>()));
    ifs.close();
    const char *c_shaderSource = rt.source.c_str();
    rt.fShader = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(rt.fShader, 1, &c_shaderSource, NULL);
    glCompileShader(rt.fShader);
    rt.program = glCreateProgram();
    glAttachShader(rt.program, rt.fShader);
    glLinkProgram(rt.program);
    rt.uniformLocs["iResolution"] = glGetUniformLocation(rt.program, "iResolution");
    logShaderCompilerErrors();
}

void Renderer::logShaderCompilerErrors()
{
    int bufflen;
    glGetShaderiv(rt.fShader, GL_INFO_LOG_LENGTH, &bufflen);

    if (bufflen > 1) {
        GLchar *log_string = new char[bufflen + 1];
        glGetShaderInfoLog(rt.fShader, bufflen, 0, log_string);
        printf("Log found for '%s.glsl':\n%s", "raytracing_fragment", log_string);
        delete log_string;
    }

    glGetShaderiv(rt.fShader, GL_COMPILE_STATUS, &bufflen);

    if (bufflen != GL_TRUE) {
        printf("Failed to compile vertex shader.");
    }
}

void Renderer::setupRenderQuad()
{
    static const GLfloat quadVertex[] = {
        -1.0f, -1.0f, 0.0f,
        1.0f, -1.0f, 0.0f,
        -1.0f,  1.0f, 0.0f,
        -1.0f,  1.0f, 0.0f,
        1.0f, -1.0f, 0.0f,
        1.0f,  1.0f, 0.0f,
    };
    glGenVertexArrays(1, &VAO);
    glBindVertexArray(VAO);
    glGenBuffers(1, &VBO);
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(quadVertex), quadVertex, GL_STATIC_DRAW);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, (GLvoid *)0);
    glEnableVertexAttribArray(0);
}

