#include "simulator.h"
#include "imgui.h"
#include "remesher/isotropicremesher.h"
#include "remesher/utils.h"

#include <igl/doublearea.h>
#include <igl/volume.h>
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <igl/invert_diag.h>
#include <igl/per_vertex_normals.h>

double compute_volume(Eigen::MatrixXd &verts, Eigen::MatrixXi &faces) {
    // igl::volume computes the signed volume of tetrahedrons (by computing 1/6 of a scalar triple product)
    // we can use this by creating a virtual point at origin, creating tetrahedrons from it to every triangle, 
    // and summing their volume. Normal pointing away from origin -> positive volume, normal pointing toward origin -> negative volume
    // even for non-convex shapes, this will work. 
    Eigen::MatrixXd V2(verts.rows() + 1, verts.cols());
    V2.topRows(verts.rows()) = verts;
    V2.bottomRows(1).setZero();
    Eigen::MatrixXi T(faces.rows(), 4);
    T.leftCols(3) = faces;
    T.rightCols(1).setConstant(verts.rows());
    Eigen::VectorXd vol;
    igl::volume(V2, T, vol);
    return std::abs(vol.sum());
}

Simulator::Simulator(Eigen::MatrixXd&& verts, Eigen::MatrixXi&& faces, double beta, double delta_t): 
    verts(verts), faces(faces), velocities(), params(beta, delta_t), current_step(0), current_time(0.0) {
    
    // initialize to zeros
    velocities.setZero(this->verts.rows(), 3);

    // compute initial volume
    initial_volume = compute_volume(this->verts, this->faces);
    initial_tri_count = this->faces.rows();
    
    current_volume = initial_volume; 
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

    // compute next velocity
    velocities += (accels * params.delta_t);

    // compute next position (using next velocity; semi-implicit Euler)
    verts += velocities * params.delta_t;
    
    // volume preservation
    double volume = compute_volume(verts, faces);
    Eigen::VectorXd areas;
    igl::doublearea(verts, faces, areas);
    areas /= 2.0;

    double total_area = areas.sum();
    double correction = (initial_volume - volume) / total_area;

    // std::cout << "volume: " << volume << "\n";
    
    // push each vertex along its normal (N * 3) by `correction`
    Eigen::MatrixXd normals;
    igl::per_vertex_normals(
        verts, 
        faces, 
        igl::PerVertexNormalsWeightingType::PER_VERTEX_NORMALS_WEIGHTING_TYPE_AREA,
        normals
    );

    normals.array() *= correction;
    verts += normals;
    
    // update velocities to account for volume preservation (step 5 of Algorithm 1)
    normals.array() /= params.delta_t; // normals has been scaled by correction already
    velocities += normals;
    
    // remesh
    Remesher::IsotropicRemesher remesher {
        pointVerticesToHalfedgeVertices(verts),
        pointFacesToHalfedgeFaces(faces),
        pointVerticesToHalfedgeVertices(velocities)
    };

    // preserve edges > 15 degrees (can tune this later...)
    remesher.setSharpEdgeIncludedAngle(15.0);
    // try to preserve initial mesh resolution (without this, it seems like mesh just gets simplified too much over time)
    remesher.setTargetTriangleCount(initial_tri_count);
    remesher.remesh(1);

    halfedgeMeshToPointMesh(
        remesher.remeshedHalfedgeMesh(), 
        verts, 
        faces, 
        velocities
    );

    // update internal stats
    current_step += 1;
    current_time += params.delta_t;
    current_volume = volume;
}

const Eigen::MatrixXd& Simulator::get_verts() {
    return verts;
}

const Eigen::MatrixXi& Simulator::get_faces() {
    return faces;
}

void Simulator::display_stats() {
    ImGui::Text("Current step: %d", current_step);
    ImGui::Text("Current time: %.2f", current_time);
    ImGui::Text("Current volume: %.2f", current_volume);
}