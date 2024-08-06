#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include <vector>
#include <unordered_map>
#include <unordered_set>

typedef unsigned int meshInt_t;
typedef double meshFloat_t;

struct Vertex {
    meshFloat_t x, y, z;
};

struct Edge {
    meshInt_t v1, v2;
};

struct Face {
    meshInt_t v1, v2, v3;
};

class TriangleMesh {
  public:
    meshInt_t addVertex(meshFloat_t x, meshFloat_t y, meshFloat_t z) {
        vertices.push_back({x, y, z});
        return vertices.size() - 1;
    }

    meshInt_t addFace(meshInt_t v1, meshInt_t v2, meshInt_t v3) {
        faces.push_back({v1, v2, v3});
        addEdge(v1, v2);
        addEdge(v2, v3);
        addEdge(v3, v1);
        return faces.size() - 1;
    }

    void addEdge(meshInt_t v1, meshInt_t v2) {
        if (v1 > v2)
            std::swap(v1, v2);
        edgeMap[v1].insert(v2);
    }

    bool loadFromOBJ(const std::string& filename) {
        std::ifstream file(filename);
        if (!file.is_open()) {
            std::cerr << "Failed to open this obj file!" << std::endl;
            return false;
        }

        std::string line;
        while (std::getline(file, line)) {
            std::istringstream iss(line);
            std::string prefix;
            iss >> prefix;

            if (prefix == "v") {
                // Vertex position
                double x, y, z;
                iss >> x >> y >> z;
                addVertex(x, y, z);
            } else if (prefix == "f") {
                int v1, v2, v3;
                char slash;
                std::string vertex1, vertex2, vertex3;
                iss >> vertex1 >> vertex2 >> vertex3;

                // We don't read normals
                v1 = std::stoi(vertex1.substr(0, vertex1.find('/'))) - 1;
                v2 = std::stoi(vertex2.substr(0, vertex2.find('/'))) - 1;
                v3 = std::stoi(vertex3.substr(0, vertex3.find('/'))) - 1;

                addFace(v1, v2, v3);
            } else if (prefix == "o" || prefix == "g" || prefix == "s" || prefix == "mtllib" || prefix == "usemtl") {
                // Skip object names, groups, smoothing groups, materials
                // Parsing only to move to next line, no data storage needed
            } else {
                // Skip other lines (such as comments)
            }
        }

        return true;
    }

    void removeVertex(meshInt_t vertexIndex, bool removeAdjFaces = false) {
        vertices.erase(vertices.begin() + vertexIndex);
    }

    // Remove a face
    void removeFace(meshInt_t faceIndex) {
        if (faceIndex < 0 || faceIndex >= faces.size())
            return;
        auto& face = faces[faceIndex];

        removeEdge(face.v1, face.v2);
        removeEdge(face.v2, face.v3);
        removeEdge(face.v3, face.v1);

        faces.erase(faces.begin() + faceIndex);
    }

  private:
    std::vector<Vertex> vertices;
    std::vector<Face> faces;
    std::unordered_map<meshInt_t, std::unordered_set<meshInt_t>> edgeMap;

    void removeEdge(meshInt_t v1, meshInt_t v2) {
        if (v1 > v2)
            std::swap(v1, v2);
        if (edgeMap.find(v1) != edgeMap.end()) {
            edgeMap[v1].erase(v2);
            if (edgeMap[v1].empty()) {
                edgeMap.erase(v1);
            }
        }
    }
};

int main() {
    TriangleMesh mesh;
    mesh.loadFromOBJ("test.obj");

    return 0;
}
