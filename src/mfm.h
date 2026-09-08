#pragma once
#include <Eigen/Dense>
#include <vector>

#include "config.h"
#include "particles.h"  // Needed for CellList

class Cooling;

class GasParticleSystem {
   public:
    size_t num_particles = 0;

    // Core physical quantities
    std::vector<double> pos_x;  // Comoving coordinates [Code Length]
    std::vector<double> pos_y;
    std::vector<double> pos_z;
    std::vector<double>
        vel_x;  // Peculiar velocity (v = a*x_dot) [Code Velocity]
    std::vector<double> vel_y;
    std::vector<double> vel_z;
    std::vector<double>
        acc_x;  // Peculiar gravitational acceleration [Code Vel / Code Time]
    std::vector<double> acc_y;
    std::vector<double> acc_z;
    std::vector<double> mass;  // Particle mass [Code Mass]

    // MFM-specific quantities
    std::vector<double> hydro_acc_x;  // Peculiar hydrodynamic acceleration
                                      // [Code Vel / Code Time]
    std::vector<double> hydro_acc_y;
    std::vector<double> hydro_acc_z;
    std::vector<double> h;    // Comoving smoothing length [Code Length]
    std::vector<double> rho;  // Comoving density (rho_c = a^3 * rho_phys) [Code
                              // Mass / Code Length^3]
    std::vector<double>
        pressure;  // Comoving pressure (P_c = (gamma-1) * rho_c * u) [Code Mass
                   // / (Code Length * Code Time^2)]
    std::vector<double>
        total_energy;           // Specific total energy [Code Velocity^2]
    std::vector<double> u;      // Specific internal energy [Code Velocity^2]
    std::vector<double> du_dt;  // Rate of change of specific internal energy
                                // [Code Velocity^2 / Code Time]
    std::vector<double> de_dt;  // Rate of change of specific total energy [Code
                                // Velocity^2 / Code Time]
    std::vector<double>
        metal_frac;  // Dimensionless metallicity mass fraction [0.0 to 1.0]
    std::vector<Eigen::Matrix3d>
        B_matrix;  // Geometric inverse matrix [1 / Code Length^2]

    std::vector<double>
        entropy;  // Entropic function (S = P / rho^gamma)
                  // [Code Velocity^2 / (Code Mass / Code Length^3)^(gamma-1)]
    std::vector<double> max_rel_ke;  // Maximum specific relative kinetic energy
                                     // of neighbors [Code Velocity^2]
    std::vector<double>
        delta_E_grav;  // Specific gravitational energy variation across
                       // smoothing length (|a_grav| * h) [Code Velocity^2]

    // Gradients are evaluated with respect to comoving coordinates
    // (d/dx_comoving)
    std::vector<Eigen::Vector3d> grad_rho;  // [Density / Code Length]
    std::vector<Eigen::Vector3d> grad_vx;   // [Velocity / Code Length]
    std::vector<Eigen::Vector3d> grad_vy;
    std::vector<Eigen::Vector3d> grad_vz;
    std::vector<Eigen::Vector3d> grad_p;  // [Pressure / Code Length]

    std::vector<double> zeta;  // Correction term for adaptive gravity softening

    Grid3D gas_rho;  // Gridded comoving gas density for PM gravity/diagnostics
                     // [Code Mass / Code Length^3]

    size_t cooling_failed_cells = 0;
    size_t cooling_total_cycles = 0;
    double accumulated_radiated_energy = 0.0;
    double accumulated_photoheating_energy = 0.0;
    double accumulated_gravitational_work = 0.0;
    double accumulated_expansion_work = 0.0;
    double accumulated_entropy_switch_energy = 0.0;
    double pressure_floor = 0.0;
    size_t ill_conditioned_cases = 0;
    size_t clamped_area_cases = 0;

    // Debugging
    std::vector<double> cond_num;
    std::vector<Eigen::Vector3d> raw_sum_p;

    // Dynamic Spatial Hashing
    CellList sph_cell_list;
    CellList pm_cell_list;
    std::vector<CIC_Data> cic_data;
    double max_h = 0.0;
    int hash_grid_dim = 0;
    double hash_cell_size = 0.0;
    double max_accel_sq = 0.0;

    GasParticleSystem(const Config& config);

    void add_particle(double px, double py, double pz, double vx, double vy,
                      double vz, double m, double initial_u, double initial_h,
                      double z_metal);

    void compute_density_and_h(const Config& config, const ParticleSystem& dm);

    // Gravity methods
    void bin_and_assign_mass(const Config& config);
    void interpolate_cic_forces(const Grid3D& ax_grid, const Grid3D& ay_grid,
                                const Grid3D& az_grid, const Config& config);

    // Short-range PP Gravity
    void compute_and_add_pp_forces(const Config& config, Diagnostics& diag);
    void compute_cross_pp_forces(ParticleSystem& dm, const Config& config,
                                 Diagnostics& diag);

    const Grid3D& get_rho() const { return gas_rho; }

    void apply_cooling(double dt, double a, const Config& config,
                       Cooling& cooling);

    // Time step calculations
    double get_gravity_timestep(const Config& config) const;
    double get_cfl_timestep(const Config& config) const;
    double get_cooling_timestep(double a, const Config& config,
                                Cooling& cooling) const;

    void hydro_step(const Config& config, double a, double dt);

    // Sync tracked internal energy, total energy, and compute pressure
    void update_primitive_variables(const Config& config, double a);

    // Compute the least-squares matrix gradients for all primitive variables
    void compute_gradients(const Config& config);

    // Solve the Riemann problem between neighbors and update hydro_acc and
    // du_dt
    void compute_hydro_forces(const Config& config, double a, double dt);

   private:
    void build_spatial_hash(double domain_size);

    // Cubic Spline Kernel (Monaghan 1992)
    // Support radius is 1h. Returns W (value) and dWdr (derivative).
    inline void kernel_cubic_spline(double r, double h, double& W,
                                    double& dWdh) const;

    void sort_arrays(const std::vector<int>& cell_start,
                     const std::vector<int>& particle_cell_idx);
};

struct ParticleState {
    Eigen::Vector3d pos;  // Comoving position
    Eigen::Vector3d vel;  // Peculiar velocity
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

// Computes spatial gradients using least-squares matrix inversion.
// INPUT: Central particle state, vector of neighbor states, comoving domain
// size.
// OUTPUT: Populated ParticleGradients struct containing the B_matrix and
// local slopes.
ParticleGradients compute_single_particle_gradients(
    const ParticleState& p_i, const std::vector<ParticleState>& neighbors,
    double domain_size);

// The extrapolated states at the exact midpoint between particle i and j
struct ReconstructedFace {
    double rho_L, rho_R;       // Left (i) and Right (j) comoving densities
    double p_L, p_R;           // Left (i) and Right (j) comoving pressures
    Eigen::Vector3d v_L, v_R;  // Left (i) and Right (j) peculiar velocities

    Eigen::Vector3d
        n;          // Unit normal vector pointing from i to j [Dimensionless]
    double r;       // Comoving distance between particles [Code Length]
    bool is_valid;  // False if particles perfectly overlap (r2 < 1e-24)
};

// Limits the reconstructed face value between phi_i and phi_j
double pairwise_slope_limiter(double phi_i, double dphi, double phi_j);

// Extrapolates particle states to the face midpoint using spatial gradients.
// INPUT: Base states and gradients for i and j.
// OUTPUT: Reconstructed left/right states at the face, limited to prevent
// overshoots.
ReconstructedFace compute_face_reconstruction(const ParticleState& p_i,
                                              const ParticleGradients& grad_i,
                                              const ParticleState& p_j,
                                              const ParticleGradients& grad_j,
                                              double domain_size);

// Output of the Riemann Solver
struct MFMFaceFlux {
    Eigen::Vector3d
        flux_mom;  // Momentum flux density vector (P_star * n) [Pressure Units]
    double P_star;  // Resolved face pressure [Pressure Units]
    double S_star;  // Resolved relative contact wave speed [Velocity Units]
};

// Solves the HLLC Riemann problem at the moving face between particles.
// INPUT: Reconstructed face states, relative frame velocity (v_frame), and
// adiabatic index.
// OUTPUT: The resolved pressure (P_star) and relative wave speed (S_star).
MFMFaceFlux solve_mfm_riemann(const ReconstructedFace& face,
                              const Eigen::Vector3d& v_frame, double gamma);