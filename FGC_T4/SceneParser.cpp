#include "SceneParser.h"


Scene::Parser::Parser(std::ifstream sFile)
{
    std::string objectType = "";
    std::string token;
    std::string line;
    bool startingBlock = false;

    while (std::getline(sFile, line)) {
        // Object Delimiter
        if (line == "---") { objectType == ""; }

        // Delete Comments
        auto it = line.find_first_of('#');

        if (it != std::string::npos) {
            line[it] = '\0';
        }

        std::stringstream sLine(line);

        while (line != "---" && std::getline(sLine, token, ':')) {
            std::string value;

            if (token == "tipo") {
                sLine >> std::ws;					// Ignore White Spaces
                std::getline(sLine, objectType);
                scene[objectType].push_back(Node());
            }

            if (objectType != "" && token != "tipo") {
                sLine >> std::ws;					// Ignore White Spaces
                std::getline(sLine, value);
                scene[objectType].back()[token] = value;
            }
        }
    }

    return;
}

Scene::Parser::~Parser()
{
}