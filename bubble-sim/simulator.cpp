#include "simulator.h"
#include "imgui.h"
#include "perlin.h"
#include "remesher/isotropicremesher.h"
#include "remesher/utils.h"

#include <functional>
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

Simulator::Simulator(Eigen::MatrixXd&& verts, Eigen::MatrixXi&& faces, SimParameters params): 
    verts(verts), faces(faces), params(params), current_step(0), current_time(0.0) {
    // initialize to zeros
    velocities.setZero(this->verts.rows(), 3);
    velocities_ext.setZero(velocities.rows(), 3);

    // compute initial volume
    initial_volume = compute_volume(this->verts, this->faces);
    initial_tri_count = this->faces.rows();
    
    current_volume = initial_volume;
    current_ke = 0.0;

    wind_x = generate_permutation(42);
    wind_y = generate_permutation(2);
    wind_z = generate_permutation(77);
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

    // compute external forces
    Eigen::MatrixXd accels_ext;
    accels_ext.setZero(accels.rows(), accels.cols());
    // in opengl coordinate system, -y is down
    accels_ext.col(1).setConstant(params.gravity);

    Eigen::MatrixXd wind_force;
    compute_wind_force(verts, wind_force);
    wind_force.array() *= 0.05;
    accels_ext += wind_force;

    // compute next velocity
    velocities += (accels * params.delta_t) + (accels_ext * params.delta_t);

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

    // compute (something proportional to) total KE for debug purposes
    current_ke = (velocities.array() * velocities.array()).sum();
    
    
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

void Simulator::compute_wind_force(const Eigen::MatrixXd& vertices, Eigen::MatrixXd& wind_force) {
    // ensure it is of the right shape
    wind_force.resizeLike(vertices);
    perlin_curl(
        vertices,
        0.2,
        wind_force, 
        wind_x, 
        wind_y, 
        wind_z
    );
}

// visualizes wind vector field at N*N*N evenly-spaced grid across [-bound, bound] \times [-bound, bound] \times [-bound, bound]
void Simulator::visualize_wind(std::function<void(const Eigen::MatrixXd&, const Eigen::MatrixXd&, const Eigen::MatrixXd&)> edge_cb) {
    int N = 15;
    double bound = 5.0;

    Eigen::MatrixXd points;
    points.resize(N * N * N, 3);

    int i = 0;
    for (int x = 0; x < N; x += 1) {
        for (int y = 0; y < N; y += 1) {
            for (int z = 0; z < N; z += 1) {
                points(i, 0) = -bound + 2.0 * bound * ((double) x / (N - 1));
                points(i, 1) = -bound + 2.0 * bound * ((double) y / (N - 1));
                points(i, 2) = -bound + 2.0 * bound * ((double) z / (N - 1));
                
                i += 1;
            }
        }
    }

    // compute wind force
    Eigen::MatrixXd wind;

    compute_wind_force(points, wind);
    
    Eigen::MatrixXd P2 = points + wind;
    Eigen::MatrixXd C;
    C.resize(1, 3);
    C(0, 0) = 1.0;
    C(0, 1) = 0.0;
    C(0, 2) = 0.0;
    edge_cb(points, P2, C);
}

const Eigen::MatrixXd& Simulator::get_verts() {
    return verts;
}

const Eigen::MatrixXi& Simulator::get_faces() {
    return faces;
}

int Simulator::get_step() {
    return current_step;
}

void Simulator::display_stats() {
    ImGui::Text("Current step: %d", current_step);
    ImGui::Text("Current time: %.2f", current_time);
    ImGui::Text("Current volume: %.2f", current_volume);
    ImGui::Text("Current kinetic energy: %.2f", current_ke);
}