#pragma once
#include <Eigen/Eigen>
#include <string>

void save_mesh_as_obj(const std::string &filename, const Eigen::MatrixXd &vertices, const Eigen::MatrixXi &faces);
bool read_obj_as_mesh(const std::string &filename, Eigen::MatrixXd &vertices, Eigen::MatrixXi &faces);