#pragma once
#include <Eigen/Dense>
#include <vector>

namespace Reconstruction {

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
    Eigen::Matrix3d B_matrix;
    Eigen::Vector3d grad_rho;
    Eigen::Vector3d grad_p;
    Eigen::Vector3d grad_vx;
    Eigen::Vector3d grad_vy;
    Eigen::Vector3d grad_vz;
    bool ill_conditioned;

    // Debugging
    double condition_number;
    Eigen::Vector3d raw_sum_p;
};

// The extrapolated states between particle i and j
struct ReconstructedFace {
    double rho_L, rho_R;       // Left (i) and Right (j) densities
    double p_L, p_R;           // Left (i) and Right (j) pressures
    Eigen::Vector3d v_L, v_R;  // Left (i) and Right (j) velocities

    Eigen::Vector3d
        n;          // Unit normal vector pointing from i to j [Dimensionless]
    double r;       // Comoving distance between particles [Code Length]
    bool is_valid;  // False if particles overlap
};

// Computes spatial gradients using least-squares matrix inversion
ParticleGradients compute_single_particle_gradients(
    const ParticleState& p_i, const std::vector<ParticleState>& neighbors,
    double domain_size);

// Extrapolates particle states to the face using spatial gradients
ReconstructedFace compute_face_reconstruction(
    const ParticleState& p_i, const ParticleGradients& grad_i,
    const ParticleState& p_j, const ParticleGradients& grad_j,
    double domain_size, double density_floor, double pressure_floor);

}  // namespace Reconstruction