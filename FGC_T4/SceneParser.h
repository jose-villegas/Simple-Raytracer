#pragma once
#include <iostream>
#include <unordered_map>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>

namespace Scene {

    typedef std::unordered_map<std::string, std::string> Node;

    class Parser {
        public:
            std::unordered_map<std::string, std::vector<Node>> scene;
        public:
            Parser(std::ifstream sFile);
            ~Parser();
    };
}