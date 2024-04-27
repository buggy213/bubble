#pragma once
#include <Eigen/Eigen>
#include <string>

void save_mesh_as_obj(std::string filename, const Eigen::MatrixXd &vertices, const Eigen::MatrixXi &faces);