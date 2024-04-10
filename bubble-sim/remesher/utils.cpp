#include "utils.h"

std::vector<Vector3> *pointVerticesToHalfedgeVertices(const Eigen::MatrixXd& vertices) {
    assert(vertices.cols() == 3);
    std::vector<Vector3> *vertices_vector = new std::vector<Vector3>;
    for (int row = 0; row < vertices.rows(); row += 1) {
        auto &v = vertices.row(row);
        // std::cout << v.x() << ", " << v.y() << ", " << v.z() << std::endl;
        vertices_vector->push_back(Vector3(v.x(), v.y(), v.z()));
    }
    
    return vertices_vector;
}

std::vector<std::vector<size_t>> *pointFacesToHalfedgeFaces(const Eigen::MatrixXi& faces) {
    assert(faces.cols() == 3);
    std::vector<std::vector<size_t>> *faces_vector = new std::vector<std::vector<size_t>>;

    for (int row = 0; row < faces.rows(); row += 1) {
        auto &v = faces.row(row);
        // std::cout << v.x() << ", " << v.y() << ", " << v.z() << std::endl;
        faces_vector->push_back(std::vector<size_t> {
            (size_t) v.x(), (size_t) v.y(), (size_t) v.z()
        });
    }

    return faces_vector;
}

void halfedgeMeshToPointMesh(IsotropicHalfedgeMesh *mesh, Eigen::MatrixXd &verts, Eigen::MatrixXi &faces, Eigen::MatrixXd &velocities) {
    size_t vertex_count = mesh->vertex_count();
    size_t face_count = mesh->face_count();
    // std::cout << "verts: " << vertex_count << std::endl;
    // std::cout << "faces: " << face_count << std::endl; 
    
    verts.resize(vertex_count, 3);
    faces.resize(face_count, 3);
    velocities.resize(vertex_count, 3);
    
    IsotropicHalfedgeMesh::Face *face = nullptr;

    std::unordered_map<IsotropicHalfedgeMesh::Vertex*, size_t> vertex_map;

    size_t face_index = 0;
    while ((face = mesh->moveToNextFace(face)) != nullptr) {
        IsotropicHalfedgeMesh::Vertex *v0 = face->halfedge->startVertex;
        IsotropicHalfedgeMesh::Vertex *v1 = face->halfedge->nextHalfedge->startVertex;
        IsotropicHalfedgeMesh::Vertex *v2 = face->halfedge->nextHalfedge->nextHalfedge->startVertex;
        
        size_t index0, index1, index2;

        auto find_index_and_insert_if_not_present = [&](IsotropicHalfedgeMesh::Vertex *v, size_t &index) {
            if (vertex_map.find(v) != vertex_map.end()) {
                index = vertex_map[v];
                // std::cout << "saw index again: " << index << " (" << v << ")" << std::endl;
            }
            else {
                index = vertex_map.size();
                // std::cout << "saw index for first time: " << index << " (" << v << ")" << std::endl;
                // std::cout << "position is " << v->position << std::endl;
                vertex_map[v] = index;

                verts(index, 0) = v->position.x();
                verts(index, 1) = v->position.y();
                verts(index, 2) = v->position.z();
                velocities(index, 0) = v->velocity.x();
                velocities(index, 1) = v->velocity.y();
                velocities(index, 2) = v->velocity.z();
            }
        };

        find_index_and_insert_if_not_present(v0, index0);
        find_index_and_insert_if_not_present(v1, index1);
        find_index_and_insert_if_not_present(v2, index2);


        // should preserve winding order
        // std::cout << "inserting face " << face_index << ": " << index0 << ", " << index1 << ", " << index2 << std::endl;
        faces(face_index, 0) = index0;
        faces(face_index, 1) = index1;
        faces(face_index, 2) = index2;
        face_index += 1;
    }
}