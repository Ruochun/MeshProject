#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include <cmath>

#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <stack>
#include <tuple>

typedef unsigned int meshInt_t;
typedef double meshFloat_t;

struct myHasher {
    template <class T1, class T2>
    std::size_t operator()(const std::pair<T1, T2>& p) const {
        auto hash1 = std::hash<T1>{}(p.first);
        auto hash2 = std::hash<T2>{}(p.second);
        return hash1 ^ hash2;
    }
};

// Helper function to calculate the angle between two vectors
meshFloat_t calculateAngleBetweenVectors(const std::tuple<meshFloat_t, meshFloat_t, meshFloat_t>& vecA,
                                         const std::tuple<meshFloat_t, meshFloat_t, meshFloat_t>& vecB) {
    // Calculate dot product and magnitudes
    meshFloat_t dotProduct = std::get<0>(vecA) * std::get<0>(vecB) + std::get<1>(vecA) * std::get<1>(vecB) +
                             std::get<2>(vecA) * std::get<2>(vecB);
    meshFloat_t magnitudeA = std::sqrt(std::get<0>(vecA) * std::get<0>(vecA) + std::get<1>(vecA) * std::get<1>(vecA) +
                                       std::get<2>(vecA) * std::get<2>(vecA));
    meshFloat_t magnitudeB = std::sqrt(std::get<0>(vecB) * std::get<0>(vecB) + std::get<1>(vecB) * std::get<1>(vecB) +
                                       std::get<2>(vecB) * std::get<2>(vecB));

    meshFloat_t cosAngle = dotProduct / (magnitudeA * magnitudeB);

    // Clamp cosAngle to [-1, 1] to avoid numerical issues
    cosAngle = std::max(-1.0, std::min(1.0, cosAngle));

    return std::acos(cosAngle);
}

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
    /* Start of 1a */
    // NOTE! It cannot correctly read mesh with quad section and other non-triangle sections.
    // But even so, it is consistent with the requirement of the assignment
    int loadOBJ(const std::string& filename) {
        std::ifstream file(filename);
        if (!file.is_open()) {
            std::cerr << "Failed to open this obj file!" << std::endl;
            return 1;
        }

        std::string line;
        while (std::getline(file, line)) {
            std::istringstream iss(line);
            std::string prefix;
            iss >> prefix;

            if (prefix == "v") {
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
            }

            // if (prefix == "o" || prefix == "g" || prefix == "s" || prefix == "mtllib" || prefix == "usemtl")
            // Don't read them since we don't have to store them anyway....
        }

        file.close();
        return 0;
    }

    int writeOBJ(const std::string& filename) const {
        std::ofstream file(filename);
        if (!file.is_open()) {
            std::cerr << "Failed to open this obj file for output!" << std::endl;
            return 1;
        }

        for (const auto& vertex : vertices) {
            file << "v " << vertex.x << " " << vertex.y << " " << vertex.z << "\n";
        }

        // 1-based indexing
        for (const auto& face : faces) {
            file << "f " << (face.v1 + 1) << " " << (face.v2 + 1) << " " << (face.v3 + 1) << "\n";
        }

        file.close();
        return 0;
    }

    /* Start of 1b */
    // Get adjacent faces for a given vertex ID
    std::vector<meshInt_t> getAdjacentFacesFromVertex(meshInt_t vertexId) const {
        std::vector<meshInt_t> adjacentFaces;
        for (size_t i = 0; i < faces.size(); ++i) {
            const auto& face = faces[i];
            if (face.v1 == vertexId || face.v2 == vertexId || face.v3 == vertexId) {
                adjacentFaces.push_back(i);
            }
        }
        return adjacentFaces;
    }

    // Get adjacent vertices for a given vertex ID
    std::unordered_set<meshInt_t> getAdjacentVertices(meshInt_t vertexId) const {
        std::unordered_set<meshInt_t> adjacentVertices;
        for (const auto& face : faces) {
            if (face.v1 == vertexId) {
                adjacentVertices.insert(face.v2);
                adjacentVertices.insert(face.v3);
            } else if (face.v2 == vertexId) {
                adjacentVertices.insert(face.v1);
                adjacentVertices.insert(face.v3);
            } else if (face.v3 == vertexId) {
                adjacentVertices.insert(face.v1);
                adjacentVertices.insert(face.v2);
            }
        }
        return adjacentVertices;
    }

    // Method to get adjacent faces for a given face ID
    std::vector<meshInt_t> getAdjacentFacesFromFace(meshInt_t faceId) const {
        if (faceId < 0 || faceId >= faces.size())
            return {};

        const auto& face = faces[faceId];
        std::unordered_set<meshInt_t> vertices = {face.v1, face.v2, face.v3};
        std::vector<meshInt_t> adjacentFaces;

        for (size_t i = 0; i < faces.size(); ++i) {
            if (i == faceId)
                continue;
            const auto& otherFace = faces[i];
            std::unordered_set<meshInt_t> otherVertices = {otherFace.v1, otherFace.v2, otherFace.v3};
            std::unordered_set<meshInt_t> both_have;
            for (meshInt_t v : vertices) {
                if (otherVertices.find(v) != otherVertices.end()) {
                    both_have.insert(v);
                }
            }
            // If two faces share at least two vertices, they are adjacent
            if (both_have.size() >= 2) {
                adjacentFaces.push_back(i);
            }
        }
        return adjacentFaces;
    }

    /* Start of 1c */
    std::vector<Vertex>& getVertices() { return vertices; }

    std::vector<Face>& getFaces() { return faces; }

    /* Start of 1d */
    Vertex getVertex(meshInt_t vertexIndex) { return vertices[vertexIndex]; }

    /* Start of 1e */
    // Remove a vertex and optionally all associated faces and edges
    void removeVertex(meshInt_t vertexIndex, bool removeAdjFaces = false) {
        if (vertexIndex < 0 || vertexIndex >= vertices.size())
            return;

        if (removeAdjFaces) {
            std::vector<meshInt_t> adjFaces = std::move(getAdjacentFacesFromVertex(vertexIndex));
            for (meshInt_t i : adjFaces) {
                removeFace(i);
            }
        }

        vertices.erase(vertices.begin() + vertexIndex);

        if (removeAdjFaces)
            vertexRemovalPostProcess(vertexIndex);
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

    /* Start of 1f */
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

    /* Start of 1g */
    // Flip orientation in the sense of RHR
    void flipOrientation(meshInt_t faceIndex) {
        if (faceIndex < 0 || faceIndex >= faces.size()) {
            return;
        }
        std::swap(faces[faceIndex].v1, faces[faceIndex].v2);
    }

    /* Start of P2 */
    bool isOrientationConsistent() const {
        if (faces.empty())
            return true;

        std::unordered_set<meshInt_t> visitedFaces;
        std::queue<meshInt_t> faceQueue;
        std::vector<uint8_t> faceNormals(faces.size(), 0);

        // Start with the first face and use its normal as the reference
        faceQueue.push(0);
        bool flag = true;
        visitedFaces.insert(0);

        while (!faceQueue.empty()) {
            meshInt_t currentFaceIndex = faceQueue.front();
            faceQueue.pop();
            const Face& currentFace = faces[currentFaceIndex];

            // Get adjacent faces and check their orientations
            std::vector<meshInt_t> adjacentFaces = getAdjacentFacesFromFace(currentFaceIndex);
            for (meshInt_t adjFaceIndex : adjacentFaces) {
                // Move on
                if (visitedFaces.find(adjFaceIndex) != visitedFaces.end())
                    continue;

                const Face& adjFace = faces[adjFaceIndex];
                // Find shared edge
                meshInt_t sharedVertices = 0;
                int shared1 = -1, shared2 = -1;

                bool jump = true;
                if (currentFace.v1 == adjFace.v1 || currentFace.v1 == adjFace.v2 || currentFace.v1 == adjFace.v3)
                    sharedVertices++, shared1 = currentFace.v1;
                if (currentFace.v2 == adjFace.v1 || currentFace.v2 == adjFace.v2 || currentFace.v2 == adjFace.v3) {
                    if (sharedVertices == 1) {
                        shared2 = currentFace.v2;
                    } else {
                        shared1 = currentFace.v2;
                    }
                    sharedVertices++;
                }
                if (currentFace.v3 == adjFace.v1 || currentFace.v3 == adjFace.v2 || currentFace.v3 == adjFace.v3) {
                    if (sharedVertices >= 1) {
                        shared2 = currentFace.v3;
                    } else {
                        shared1 = currentFace.v3;
                    }
                }
                if (currentFace.v1 == adjFace.v2 || currentFace.v2 == adjFace.v2 || currentFace.v3 == adjFace.v2)
                    jump = false;

                // std::cout << "currentFace1: " << currentFace.v1 << std::endl;
                // std::cout << "currentFace2: " << currentFace.v2 << std::endl;
                // std::cout << "currentFace3: " << currentFace.v3 << std::endl;
                // std::cout << "Shared1: " << shared1 << std::endl;
                // std::cout << "Shared2: " << shared2 << std::endl;

                // Calculate orientation based on shared edge
                if ((currentFace.v1 == shared1 && currentFace.v2 == shared2) ||
                    (currentFace.v2 == shared1 && currentFace.v3 == shared2) ||
                    (currentFace.v3 == shared1 && currentFace.v1 == shared2)) {
                    if (!jump) {
                        flag = false;
                        break;
                    }
                } else {
                    if (jump) {
                        flag = false;
                        break;
                    }
                }

                visitedFaces.insert(adjFaceIndex);
                faceQueue.push(adjFaceIndex);
            }
        }

        return flag;
    }

    /* Start of P3 */
    meshInt_t getNumLoops() const {
        // Map to track edges and their counts in the mesh
        std::unordered_map<std::pair<meshInt_t, meshInt_t>, meshInt_t, myHasher> edgeCount;

        // Count the number of times each edge appears
        for (const auto& face : faces) {
            // Get edges for each face
            std::pair<meshInt_t, meshInt_t> e1 = std::minmax(face.v1, face.v2);
            std::pair<meshInt_t, meshInt_t> e2 = std::minmax(face.v2, face.v3);
            std::pair<meshInt_t, meshInt_t> e3 = std::minmax(face.v3, face.v1);

            edgeCount[e1]++;
            edgeCount[e2]++;
            edgeCount[e3]++;
        }

        // Find boundary edges (those appearing only once)
        std::unordered_map<meshInt_t, std::unordered_set<meshInt_t>> boundaryEdges;
        for (const auto& [edge, count] : edgeCount) {
            if (count == 1) {
                boundaryEdges[edge.first].insert(edge.second);
                boundaryEdges[edge.second].insert(edge.first);
            }
        }

        // Count loops by performing DFS over boundary edges
        std::unordered_set<meshInt_t> visitedVertices;
        meshInt_t loopCount = 0;

        for (const auto& [vertex, neighbors] : boundaryEdges) {
            if (visitedVertices.find(vertex) == visitedVertices.end()) {
                loopCount++;
                std::stack<meshInt_t> stack;
                stack.push(vertex);

                while (!stack.empty()) {
                    meshInt_t current = stack.top();
                    stack.pop();
                    visitedVertices.insert(current);

                    for (meshInt_t neighbor : boundaryEdges[current]) {
                        if (visitedVertices.find(neighbor) == visitedVertices.end()) {
                            stack.push(neighbor);
                        }
                    }
                }
            }
        }

        return loopCount;
    }

    /* Start of P4 */
    std::vector<meshInt_t> getFacesWithSmallAngles(meshFloat_t angleDeg) const {
        std::vector<meshInt_t> facesWithSmallAngles;

        // Convert threshold from degrees to radians
        meshFloat_t angleRad = angleDeg * M_PI / 180.0;

        for (size_t i = 0; i < faces.size(); ++i) {
            const auto& face = faces[i];
            const Vertex& v1 = vertices[face.v1];
            const Vertex& v2 = vertices[face.v2];
            const Vertex& v3 = vertices[face.v3];

            // Calculate vectors for each edge of the triangle
            auto vec1 = std::make_tuple(v2.x - v1.x, v2.y - v1.y, v2.z - v1.z);
            auto vec2 = std::make_tuple(v3.x - v1.x, v3.y - v1.y, v3.z - v1.z);
            auto vec3 = std::make_tuple(v3.x - v2.x, v3.y - v2.y, v3.z - v2.z);
            auto vec4 = std::make_tuple(v1.x - v2.x, v1.y - v2.y, v1.z - v2.z);

            // Calculate angles between vectors
            meshFloat_t angle1 = calculateAngleBetweenVectors(vec1, vec2);
            meshFloat_t angle2 = calculateAngleBetweenVectors(vec3, vec4);
            meshFloat_t angle3 = M_PI - angle1 - angle2;  // Third angle

            // Check if any angle is below the threshold
            if (angle1 < angleRad || angle2 < angleRad || angle3 < angleRad) {
                facesWithSmallAngles.push_back(i);
            }
        }

        return facesWithSmallAngles;
    }

    /* Start of P5 */
    void collapseEdges(meshFloat_t lengthThreshold) {
        std::unordered_map<meshInt_t, meshInt_t> vertexCollapseMap;

        for (size_t i = 0; i < faces.size(); ++i) {
            const auto& face = faces[i];
            std::vector<std::pair<meshInt_t, meshInt_t>> edges = {
                {face.v1, face.v2}, {face.v2, face.v3}, {face.v3, face.v1}};

            for (const auto& [v1, v2] : edges) {
                if (vertexCollapseMap.count(v1) || vertexCollapseMap.count(v2)) {
                    continue;  // Skip already collapsed vertices
                }

                const Vertex& vert1 = vertices[v1];
                const Vertex& vert2 = vertices[v2];

                meshFloat_t dx = vert1.x - vert2.x;
                meshFloat_t dy = vert1.y - vert2.y;
                meshFloat_t dz = vert1.z - vert2.z;
                meshFloat_t lengthSquared = dx * dx + dy * dy + dz * dz;

                if (lengthSquared < lengthThreshold * lengthThreshold) {
                    // Collapse by averaging
                    vertices[v1].x = (vert1.x + vert2.x) / 2.0;
                    vertices[v1].y = (vert1.y + vert2.y) / 2.0;
                    vertices[v1].z = (vert1.z + vert2.z) / 2.0;

                    vertexCollapseMap[v2] = v1;
                }
            }
        }

        // Now collapse the faces
        for (auto& face : faces) {
            if (vertexCollapseMap.count(face.v1)) {
                face.v1 = vertexCollapseMap[face.v1];
            }
            if (vertexCollapseMap.count(face.v2)) {
                face.v2 = vertexCollapseMap[face.v2];
            }
            if (vertexCollapseMap.count(face.v3)) {
                face.v3 = vertexCollapseMap[face.v3];
            }
        }

        // Degeneracy is not tolerated
        faces.erase(std::remove_if(faces.begin(), faces.end(),
                                   [](const Face& face) {
                                       return face.v1 == face.v2 || face.v2 == face.v3 || face.v3 == face.v1;
                                   }),
                    faces.end());

        // Clean up isolated vertices
        removeIsolatedVertices();
    }

    /* Start of P6 */
    void doDiagonalEdgeSwaps(meshFloat_t angleDeg) {
        // Use what we already have...
        std::vector<meshInt_t> targetFaces = getFacesWithSmallAngles(angleDeg);
        // for (int i=0;i<targetFaces.size();i++) {
        //     std::cout << targetFaces[i] << " ";
        // }

        meshFloat_t angleRad = angleDeg * M_PI / 180.0;

        for (meshInt_t faceIndex : targetFaces) {
            const auto& face = faces[faceIndex];
            const Vertex& v1 = vertices[face.v1];
            const Vertex& v2 = vertices[face.v2];
            const Vertex& v3 = vertices[face.v3];

            // Calculate angles for the current triangle
            auto vec1 = std::make_tuple(v2.x - v1.x, v2.y - v1.y, v2.z - v1.z);
            auto vec2 = std::make_tuple(v3.x - v1.x, v3.y - v1.y, v3.z - v1.z);
            auto vec3 = std::make_tuple(v3.x - v2.x, v3.y - v2.y, v3.z - v2.z);
            auto vec4 = std::make_tuple(v1.x - v2.x, v1.y - v2.y, v1.z - v2.z);

            meshFloat_t angle1 = calculateAngleBetweenVectors(vec1, vec2);
            meshFloat_t angle2 = calculateAngleBetweenVectors(vec3, vec4);
            meshFloat_t angle3 = M_PI - angle1 - angle2;

            if ((angle1 > M_PI / 2 || angle2 > M_PI / 2 || angle3 > M_PI / 2) &&
                (angle1 < angleRad || angle2 < angleRad || angle3 < angleRad)) {
                // Determine which edge to flip
                if (angle1 < angleRad) {
                    tryDiagonalEdgeSwap(face.v1, face.v2, face.v3);
                }
                if (angle2 < angleRad) {
                    tryDiagonalEdgeSwap(face.v2, face.v3, face.v1);
                }
                if (angle3 < angleRad) {
                    tryDiagonalEdgeSwap(face.v3, face.v1, face.v2);
                }
            }
        }
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

    void vertexRemovalPostProcess(meshInt_t vertexIndex) {
        // Update face and edge indices
        for (auto& face : faces) {
            if (face.v1 > vertexIndex)
                --face.v1;
            if (face.v2 > vertexIndex)
                --face.v2;
            if (face.v3 > vertexIndex)
                --face.v3;
        }

        std::unordered_map<meshInt_t, std::unordered_set<meshInt_t>> newEdgeMap;
        for (auto& [v1, neighbors] : edgeMap) {
            if (v1 == vertexIndex)
                continue;
            meshInt_t newV1 = v1 > vertexIndex ? v1 - 1 : v1;
            for (meshInt_t v2 : neighbors) {
                if (v2 == vertexIndex)
                    continue;
                meshInt_t newV2 = v2 > vertexIndex ? v2 - 1 : v2;
                newEdgeMap[newV1].insert(newV2);
            }
        }
        edgeMap = std::move(newEdgeMap);
    }

    void removeIsolatedVertices() {
        std::unordered_set<meshInt_t> usedVertices;
        for (const auto& face : faces) {
            usedVertices.insert(face.v1);
            usedVertices.insert(face.v2);
            usedVertices.insert(face.v3);
        }

        std::vector<Vertex> newVertices;
        std::unordered_map<meshInt_t, meshInt_t> oldToNewIndex;
        meshInt_t newIndex = 0;

        for (size_t i = 0; i < vertices.size(); ++i) {
            if (usedVertices.count(i)) {
                newVertices.push_back(vertices[i]);
                oldToNewIndex[i] = newIndex++;
            }
        }

        vertices = std::move(newVertices);

        for (auto& face : faces) {
            face.v1 = oldToNewIndex[face.v1];
            face.v2 = oldToNewIndex[face.v2];
            face.v3 = oldToNewIndex[face.v3];
        }
    }

    void tryDiagonalEdgeSwap(meshInt_t v1, meshInt_t v2, meshInt_t v3) {
        // Find the adjacent face sharing the edge (v1, v2)
        meshInt_t adjacentFaceIndex;
        bool found = false;
        for (size_t i = 0; i < faces.size(); ++i) {
            const auto& face = faces[i];
            if ((face.v1 == v1 && face.v2 == v2) || (face.v1 == v2 && face.v2 == v1) ||
                (face.v1 == v2 && face.v2 == v3) || (face.v1 == v3 && face.v2 == v2) ||
                (face.v1 == v3 && face.v2 == v1) || (face.v1 == v1 && face.v2 == v3)) {
                adjacentFaceIndex = i;
                found = true;
                break;
            }
        }

        if (!found) {
            return;  // No adjacent face
        }

        const auto& adjFace = faces[adjacentFaceIndex];
        meshInt_t sharedVertex = (adjFace.v1 == v1 || adjFace.v1 == v2)   ? adjFace.v1
                                 : (adjFace.v2 == v1 || adjFace.v2 == v2) ? adjFace.v2
                                                                          : adjFace.v3;
        meshInt_t oppositeVertex = (adjFace.v1 != v1 && adjFace.v1 != v2)   ? adjFace.v1
                                   : (adjFace.v2 != v1 && adjFace.v2 != v2) ? adjFace.v2
                                                                            : adjFace.v3;

        // Perform the edge swap: replace edge (v1, v2) with (v3, oppositeVertex)
        faces[adjacentFaceIndex] = {v1, v3, oppositeVertex};
        faces[faces.size() - 1] = {v2, v3, oppositeVertex};  // Last face was swapped

        // Remove the degenerate triangle if created
        if (faces[adjacentFaceIndex].v1 == faces[adjacentFaceIndex].v2 ||
            faces[adjacentFaceIndex].v2 == faces[adjacentFaceIndex].v3 ||
            faces[adjacentFaceIndex].v3 == faces[adjacentFaceIndex].v1) {
            faces.erase(faces.begin() + adjacentFaceIndex);
        }
        if (faces[faces.size() - 1].v1 == faces[faces.size() - 1].v2 ||
            faces[faces.size() - 1].v2 == faces[faces.size() - 1].v3 ||
            faces[faces.size() - 1].v3 == faces[faces.size() - 1].v1) {
            faces.pop_back();
        }
    }
};

int main() {
    TriangleMesh mesh;
    mesh.loadOBJ("test.obj");

    if (mesh.isOrientationConsistent()) {
        std::cout << "Orientation is consistent." << std::endl;
    } else {
        std::cout << "Orientation is NOT consistent!" << std::endl;
    }

    std::cout << "Num loop: " << mesh.getNumLoops() << std::endl;

    auto smallFaces = mesh.getFacesWithSmallAngles(10);
    std::cout << "Small faces: ";
    for (int i = 0; i < smallFaces.size(); i++) {
        std::cout << smallFaces[i] << " ";
    }
    std::cout << "\n";

    mesh.doDiagonalEdgeSwaps(90);
    mesh.writeOBJ("myOBJ.obj");

    return 0;
}
