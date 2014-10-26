#include <iostream>
#include "Scene.h"

int main(int argc, char *argv[])
{
    Scene::Parser parser(std::ifstream("test.yml"));
    Scene::SceneManager sceneElements = Scene::Serializer::serializeParsedScene(parser);
    // Wait for User input to close
    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); std::cin.get();
    return 0;
}