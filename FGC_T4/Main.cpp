#include "Scene.h"
#include "Renderer.h"
#include <iostream>

int main(int argc, char *argv[])
{
    if (argc > 1) {
        int renderWidth = 400, renderHeight = 300;
        Scene::Parser *parser;
        Renderer render;
        Scene::SceneManager sceneElements;
        std::string oFile = std::string("salida");

        for (int i = 0; i < argc; i++) {
            std::string stringArg = std::string(argv[i]);

            if (stringArg.substr(0, 2) == "-i") {
                parser = new Scene::Parser(std::ifstream(stringArg.substr(2)));
                sceneElements = Scene::Serializer::serializeParsedScene(*parser);
            }

            if (stringArg.substr(0, 2) == "-o") {
                oFile = stringArg.substr(2);
            }

            if (stringArg.substr(0, 2) == "-s") {
                renderWidth = std::stoi(stringArg.substr(2, stringArg.find("x")));
                renderHeight = std::stoi(stringArg.substr(stringArg.find("x") + 1));
            }
        }

        // Show data readed from YML
        std::cout << sceneElements.toGLSL();
        // Render Scene
        render.init(sceneElements, oFile, renderWidth, renderHeight);
    }

    return 0;
}