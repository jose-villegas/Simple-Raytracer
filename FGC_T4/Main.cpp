#include "Scene.h"
#include "Renderer.h"
#include <iostream>

int main(int argc, char *argv[])
{
    Scene::Parser parser(std::ifstream("test.yml"));
    Renderer render;
    Scene::SceneManager sceneElements = Scene::Serializer::serializeParsedScene(parser);
    // Show data readed from YML
    std::cout << sceneElements.toGLSL();
    // Render Scene
    render.init(sceneElements);
    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); std::cin.get();
    return 0;
}