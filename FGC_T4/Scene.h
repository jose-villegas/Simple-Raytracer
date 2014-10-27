#pragma once
#include <iostream>
#include <unordered_map>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include "glm.hpp"

namespace Scene {
    class Parser {
        protected:
            typedef std::unordered_map<std::string, std::string> Node;
            typedef std::unordered_map<std::string, std::vector<Node>> SceneNodes;
            SceneNodes getScene() const { return scene; }
        private:
            SceneNodes scene;
        public:
            Parser(std::ifstream sFile);
            ~Parser();
            friend class Serializer;
    };

    namespace Figures {
        class Figure {
            public:
                std::string materialName;
        };

        class Sphere : public Figure {
            public:
                glm::vec3 center;
                float radius;
        };

        class Cylinder : public Figure {
            public:
                glm::vec2 center;
                char axis;
                float radius;
                float min;
                float max;
        };

        class Triangle : public Figure {
            public:
                glm::vec3 points[3];
        };

        class Cube : public Figure {
            public:
                glm::vec2 facesX;
                glm::vec2 facesY;
                glm::vec2 facesZ;
        };
    };

    class Camera {
        public:
            glm::vec3 position;
            glm::vec3 direction;
            bool enableShadows;
            bool enableAntiAlias;
    };

    class Material {
        public:
            std::string name;
            glm::vec3 diffuse;
            glm::vec3 specular;
            float reflectionIndex;
            float refractionIndex;
            float alphaValue;
    };

    class Light {
        public:
            glm::vec3 position;
            glm::vec3 color;
    };

    class SceneManager {
        public:
            std::vector<Camera> cameras;
            std::vector<Material> materials;
            std::vector<Light> lights;
            std::vector<Figures::Triangle> triangles;
            std::vector<Figures::Sphere> spheres;
            std::vector<Figures::Cylinder> cylinders;
            std::vector<Figures::Cube> cubes;
    };

    class Serializer : private Parser {
        public:
            static SceneManager serializeParsedScene(Parser parsed);
    };
}