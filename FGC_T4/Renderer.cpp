#include "Renderer.h"
#include "glm.hpp"

Renderer::Renderer(void)
{
}


Renderer::~Renderer(void)
{
}

void Renderer::init(Scene::SceneManager sceneElements, std::string outFile, int rWidth, int rHeight)
{
    this->sceneData = sceneElements;
    const int renderWidth = rWidth;
    const int renderHeight = rHeight;
    const int windowWidth = 1;
    const int windowHeight = 1;

    if (!glfwInit()) {
        exit(EXIT_FAILURE);
    }

    GLFWwindow *window = glfwCreateWindow(windowWidth, windowHeight, "Raytracing - RT", NULL, NULL);

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
    createRenderTarget(renderWidth, renderHeight);
    setupRenderQuad();
    // Load Relevant Uniforms
    glUseProgram(rt.program);
    glUniform3fv(rt.uniformLocs["iResolution"], 1, &glm::vec3(renderWidth, renderHeight, 0)[0]);

    // Draw Single Frame
    for (int i = 0; i < sceneData.cameras.size(); i++) {
        glUniform3fv(rt.uniformLocs["cameraPosition"], 1, &sceneData.cameras[i].position[0]);
        glUniform3fv(rt.uniformLocs["cameraDirection"], 1, &sceneData.cameras[i].direction[0]);
        glUniform1f(rt.uniformLocs["enableAntiAliasing"], sceneData.cameras[i].enableAntiAlias);
        glUniform1f(rt.uniformLocs["enableShadows"], sceneData.cameras[i].enableShadows);
        glBindFramebuffer(GL_FRAMEBUFFER, FBO);
        glViewport(0, 0, renderWidth, renderHeight);
        glDrawArrays(GL_TRIANGLES, 0, 6);
        // Reset Render To Screen
        glBindFramebuffer(GL_FRAMEBUFFER, 0);
        glfwSwapBuffers(window);
        // Save Render
        toFile(outFile + std::to_string(i) + ".tga", renderWidth, renderHeight);
    }

    glViewport(0, 0, windowWidth, windowHeight);
    // Destory OGL Context
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
    size_t GLSL_RECEIVER_BLOCK = rt.source.find("BEGIN:SCENEPARAMS");
    rt.source.insert(GLSL_RECEIVER_BLOCK + std::string("BEGIN:SCENEPARAMS").size(), this->sceneData.toGLSL());
    const char *c_shaderSource = rt.source.c_str();
    rt.fShader = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(rt.fShader, 1, &c_shaderSource, NULL);
    glCompileShader(rt.fShader);
    rt.program = glCreateProgram();
    glAttachShader(rt.program, rt.fShader);
    glLinkProgram(rt.program);
    rt.uniformLocs["iResolution"] = glGetUniformLocation(rt.program, "iResolution");
    rt.uniformLocs["cameraPosition"] = glGetUniformLocation(rt.program, "cameraPosition");
    rt.uniformLocs["cameraDirection"] = glGetUniformLocation(rt.program, "cameraDirection");
    rt.uniformLocs["enableAntiAliasing"] = glGetUniformLocation(rt.program, "enableAntiAliasing");
    rt.uniformLocs["enableShadows"] = glGetUniformLocation(rt.program, "enableShadows");
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

void Renderer::createRenderTarget(int width, int height)
{
    // Create Frame Buffer Object
    glGenFramebuffers(1, &FBO);
    glBindFramebuffer(GL_FRAMEBUFFER, FBO);
    // Create color render buffer
    glGenRenderbuffers(1, &CBF);
    glBindRenderbuffer(GL_RENDERBUFFER, CBF);
    glRenderbufferStorage(GL_RENDERBUFFER, GL_RGBA8, width, height);
    glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_RENDERBUFFER, CBF);
    // Create a Render Texture
    glGenTextures(1, &RTX);
    glBindTexture(GL_TEXTURE_2D, RTX);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, width, height, 0, GL_RGB, GL_UNSIGNED_BYTE, 0);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, RTX, 0);
    // Set the list of draw buffers.
    DBF = GL_COLOR_ATTACHMENT0;
    glDrawBuffers(1, &DBF); // 1 = size of drawbuffer

    if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE) {
        exit(1);
    }
}

void Renderer::toFile(std::string filename, int renderWidth, int renderHeight)
{
    glBindFramebufferEXT(GL_FRAMEBUFFER, FBO);
    glReadBuffer(GL_COLOR_ATTACHMENT0);
    // Reserve Image Size
    long imageSize = renderWidth * renderHeight * 3;
    unsigned char *data = new unsigned char[imageSize];
    glReadPixels(0, 0, renderWidth, renderHeight, GL_RGB, GL_UNSIGNED_BYTE, data);
    // Save Data
    int xa = renderWidth % 256;
    int xb = (renderWidth - xa) / 256; int ya = renderHeight % 256;
    int yb = (renderHeight - ya) / 256; // Assemble the header
    unsigned char header[18] = {0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, (char)xa, (char)xb, (char)ya, (char)yb, 24, 0};
    // Write header and data to file
    std::ofstream ofs(filename, std::ios::out | std::ios::binary);
    ofs.write(reinterpret_cast<char *>(header), sizeof(char) * 18);
    ofs.write(reinterpret_cast<char *>(data), sizeof(char)*imageSize);
    ofs.close();
    delete[] data;
    data = NULL;
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
}

