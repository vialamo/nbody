#pragma once
#include <Eigen/Dense>
#include <vector>

constexpr double N_cond_crit = 100.0;

//#define USE_MIDPOINT_QUADRATURE

namespace Reconstruction {

struct FluidStateArrays {
    const double* pos_x;
    const double* pos_y;
    const double* pos_z;
    const double* vel_x;
    const double* vel_y;
    const double* vel_z;
    const double* mass;
    const double* rho;
    const double* pressure;
    const double* h;
    const double* cond_num;
};

struct ParticleState {
    Eigen::Vector3d pos;  // Comoving position
    Eigen::Vector3d vel;  // Comoving velocity
    double mass;
    double rho;       // Comoving density
    double pressure;  // Comoving pressure
    double h;         // Comoving smoothing length
};

// Represents the comoving spatial gradients for a single particle
struct ParticleGradients {
    Eigen::Vector3d grad_rho;
    Eigen::Vector3d grad_p;
    Eigen::Vector3d grad_vx;
    Eigen::Vector3d grad_vy;
    Eigen::Vector3d grad_vz;
    bool ill_conditioned;
};

// The extrapolated states between particle i and j
struct ReconstructedFace {
    double rho_L, rho_R;       // Left (i) and Right (j) densities
    double p_L, p_R;           // Left (i) and Right (j) pressures
    Eigen::Vector3d v_L, v_R;  // Left (i) and Right (j) velocities

    Eigen::Vector3d
        n;  // Unit normal vector pointing from i to j [Dimensionless]
    Eigen::Vector3d area_vec;  // The un-normalized Face Area vector
    double r;       // Comoving distance between particles [Code Length]
    bool is_valid;  // False if particles overlap
    bool used_sph_fallback;  // Diagnostic flag
};

// Computes spatial gradients using least-squares matrix inversion
ParticleGradients compute_single_particle_gradients(
    size_t i,                             // The central particle index
    const int* neighbor_indices,          // Flat array of neighbor indices
    size_t num_neighbors,                 // Number of neighbors
    const FluidStateArrays& system_data,  // Pointers to global SoA
    const Eigen::Matrix3d& B, bool ill_conditioned, double domain_size);

// Extrapolates particle states to the face using spatial gradients
ReconstructedFace compute_face_reconstruction(
    const ParticleState& p_i, const ParticleGradients& grad_i,
    const ParticleState& p_j, const ParticleGradients& grad_j,
    const Eigen::Matrix3d& B_matrix_i, const Eigen::Matrix3d& B_matrix_j,
    double domain_size, double density_floor, double pressure_floor);

}  // namespace Reconstruction