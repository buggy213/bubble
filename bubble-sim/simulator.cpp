#include "simulator.h"
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <igl/invert_diag.h>
#include <igl/per_vertex_normals.h>

Simulator::Simulator(Eigen::MatrixXd&& verts, Eigen::MatrixXi&& faces, double beta, double delta_t): 
    verts(verts), faces(faces), velocities(), params(beta, delta_t) {
    
    // initialize to zeros
    velocities.setZero(verts.rows(), 3);
};

void Simulator::set_params(SimParameters params) {
    this->params = params;
}

void Simulator::step() {
    // compute mean curvature at each vertex (H is n * 1)
    Eigen::MatrixXd HN;
    Eigen::SparseMatrix<double> L, M, Minv;
    igl::cotmatrix(verts, faces, L);
    igl::massmatrix(verts, faces,igl::MASSMATRIX_TYPE_VORONOI,M);
    igl::invert_diag(M, Minv);
    HN = -Minv*(L*verts);
    

    // d^2{x}/dt^2 = -\beta H(x,t)n(x,t)
    Eigen::MatrixXd accels;
    
    accels = HN;
    accels.array() *= (-params.beta); 

    // compute intermediate velocity
    Eigen::MatrixXd v_next;
    v_next = velocities + (accels * params.delta_t);

    // compute next position
    verts += v_next * params.delta_t;
    velocities = v_next; // this might do a copy
}

const Eigen::MatrixXd& Simulator::get_verts() {
    return verts;
}

const Eigen::MatrixXi& Simulator::get_faces() {
    return faces;
}