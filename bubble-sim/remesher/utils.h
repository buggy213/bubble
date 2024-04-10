#ifndef ISOTROPIC_REMESHER_UTILS_H
#define ISOTROPIC_REMESHER_UTILS_H

#include <cassert>
#include <vector>
#include <unordered_map>
#include "isotropichalfedgemesh.h"
#include "vector3.h"
#include <Eigen/Eigen>

using namespace Remesher;

// convert between point mesh description (vertices + indices) to format that is expected for building up halfedge mesh
// ideally, they would be more tightly coupled so we don't need to convert between on every simulation step
std::vector<Vector3> *pointVerticesToHalfedgeVertices(const Eigen::MatrixXd& vertices);

std::vector<std::vector<size_t>> *pointFacesToHalfedgeFaces(const Eigen::MatrixXi& faces);

void halfedgeMeshToPointMesh(IsotropicHalfedgeMesh *mesh, Eigen::MatrixXd &verts, Eigen::MatrixXi &faces, Eigen::MatrixXd &velocities);

#endif