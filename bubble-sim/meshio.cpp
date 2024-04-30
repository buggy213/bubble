#include "meshio.h"
#include <Eigen/src/Core/Matrix.h>
#include <igl/per_vertex_normals.h>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>

void save_mesh_as_obj(const std::string &filename, const Eigen::MatrixXd &vertices, const Eigen::MatrixXi &faces) {
    std::ofstream obj_file {filename};
    if (!obj_file.is_open()) {
        std::cerr << "unable to open file " << filename << std::endl;
    }
    // disable scientific notation
    obj_file << std::fixed << std::setprecision(3);

    for (int i = 0; i < vertices.rows(); i += 1) {
        obj_file << "v" << " " << vertices(i, 0) << " " << vertices(i, 1) << " " << vertices(i, 2) << "\n";
    }

    for (int i = 0; i < faces.rows(); i += 1) {
        // vertex indices are 1-indexed in .obj
        obj_file << "f" << " " << faces(i, 0) + 1 << " " << faces(i, 1) + 1 << " " << faces(i, 2) + 1 << "\n";
    }

    // compute vertex normals and save those too
    Eigen::MatrixXd normals;
    igl::per_vertex_normals(
        vertices, 
        faces, 
        igl::PerVertexNormalsWeightingType::PER_VERTEX_NORMALS_WEIGHTING_TYPE_AREA,
        normals
    );

    for (int i = 0; i < normals.rows(); i += 1) {
        obj_file << "vn" << " " << normals(i, 0) << " " << normals(i, 1) << " " << normals(i, 2) << "\n";
    }
}

bool read_obj_as_mesh(const std::string &filename, Eigen::MatrixXd &vertices, Eigen::MatrixXi &faces) {
    std::ifstream obj_file{filename};
    std::string line;

    std::vector<float> vertex_data;
    std::vector<int> face_data;
    int vertex_count = 0;
    int face_count = 0;
    while (std::getline(obj_file, line)) {
        std::stringstream splitter{line};
        std::string s1, s2, s3, s4;
        splitter >> s1 >> s2 >> s3 >> s4;
        if (s1 == "v") {
            float x = std::stof(s2);
            float y = std::stof(s3);
            float z = std::stof(s4);
            vertex_data.push_back(x);
            vertex_data.push_back(y);
            vertex_data.push_back(z);
            vertex_count += 1;
        }
        else if (s1 == "f") {
            int p0 = std::stoi(s2);
            int p1 = std::stoi(s3);
            int p2 = std::stoi(s4);
            face_data.push_back(p0-1); // reverse 1-indexing
            face_data.push_back(p1-1);
            face_data.push_back(p2-1);
            face_count += 1;
        }
        else if (s1 == "#" || s1.find_first_not_of(" \t\n\v\f\r") == std::string::npos || s1 == "s") {
            // pass
        }
        else {
            return false;
        }
        
        vertices.resize(vertex_count, 3);
        int row = 0;
        int col = 0;
        for (float f : vertex_data) {
            vertices(row, col) = f;
            col += 1;
            if (col == 3) {
                row += 1;
                col = 0;
            }
        }

        faces.resize(face_count, 3);
        row = 0;
        col = 0;
        for (int i : face_data) {
            faces(row, col) = i;
            col += 1;
            if (col == 3) {
                row += 1;
                col = 0;
            }
        }
    }

    return true;
}