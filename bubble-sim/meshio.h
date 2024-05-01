#pragma once
#include <Eigen/Eigen>
#include <string>


void compute_uv_shrinkwrap(const Eigen::MatrixXd &vertices, Eigen::MatrixXd &uv);
void save_mesh_as_obj(const std::string &filename, const Eigen::MatrixXd &vertices, const Eigen::MatrixXi &faces);
bool read_obj_as_mesh(const std::string &filename, Eigen::MatrixXd &vertices, Eigen::MatrixXi &faces);