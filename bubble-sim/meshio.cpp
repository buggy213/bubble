#include "meshio.h"
#include <igl/per_vertex_normals.h>
#include <fstream>

void save_mesh_as_obj(std::string filename, const Eigen::MatrixXd &vertices, const Eigen::MatrixXi &faces) {
    std::ofstream obj_file(filename);

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
        obj_file << "vn" << " " << normals(i, 0) << " " << normals(i, 1) << normals(i, 2) << "\n";
    }
}