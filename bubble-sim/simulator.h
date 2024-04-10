#include "remesher/isotropicremesher.h"
#include "remesher/utils.h"
#include <Eigen/Eigen>

class SimParameters {
public:
    SimParameters(double beta, double delta_t): beta(beta), delta_t(delta_t) {};

    double beta; // tunable parameter related to surface tension coefficient
    double delta_t; // timestep
};

class Simulator {
public:
    // takes ownership of verts, faces
    Simulator(Eigen::MatrixXd&& verts, Eigen::MatrixXi&& faces, double beta, double delta_t);

    void step();
    void set_params(SimParameters params);
    const Eigen::MatrixXd& get_verts();
    const Eigen::MatrixXi& get_faces();

private:
    Eigen::MatrixXd verts;      // n * 3
    Eigen::MatrixXi faces;      // f * 3
    Eigen::MatrixXd velocities; // n * 3
    SimParameters params;

    double initial_volume;
    size_t initial_tri_count;
};