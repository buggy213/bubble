#pragma once
#include "perlin.h"
#include <Eigen/Eigen>

class SimParameters {
public:
    SimParameters(double beta, double delta_t, double gravity): 
        beta(beta), delta_t(delta_t), gravity(gravity) {};

    double beta; // tunable parameter related to surface tension coefficient
    double delta_t; // timestep

    double gravity;
};

class Simulator {
public:
    // takes ownership of verts, faces
    Simulator(Eigen::MatrixXd&& verts, Eigen::MatrixXi&& faces, SimParameters params);

    void compute_wind_force(const Eigen::MatrixXd& vertices, Eigen::MatrixXd& wind_force);
    void visualize_wind(std::function<void(const Eigen::MatrixXd&, const Eigen::MatrixXd&, const Eigen::MatrixXd&)> edge_cb);
    void step();
    
    void set_params(SimParameters params);
    const Eigen::MatrixXd& get_verts();
    const Eigen::MatrixXi& get_faces();
    void get_curvature(Eigen::MatrixXd &curvaure);
    int get_step();
    void display_stats();

private:
    Eigen::MatrixXd verts;      // n * 3
    Eigen::MatrixXi faces;      // f * 3
    Eigen::MatrixXd velocities; // n * 3
    Eigen::MatrixXd velocities_ext;
    SimParameters params;

    double initial_volume;
    size_t initial_tri_count;

    int current_step;
    double current_time;
    double current_volume;
    double current_ke;

    PerlinPermutation wind_x, wind_y, wind_z;
};