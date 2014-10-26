#include <iostream>
#include "SceneParser.h"

int main(int argc, char *argv[])
{
    Scene::Parser parser(std::ifstream("test.yml"));
    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); std::cin.get();
    return 0;
}