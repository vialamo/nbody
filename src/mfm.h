#pragma once
#include <Eigen/Dense>
#include <vector>

#include "config.h"
#include "lbvh.h"
#include "particles.h"
#include "reconstruction.h"

class Cooling;

struct KickWork {
    double grav_work = 0.0;
    double exp_work = 0.0;
    double hydro_exp_work = 0.0;
};

class GasParticleSystem {
   public:
    size_t num_particles = 0;

    // Core physical quantities
    std::vector<double> pos_x;  // Comoving coordinates [Code Length]
    std::vector<double> pos_y;
    std::vector<double> pos_z;
    std::vector<double> vel_x;  // Comoving velocity [Code Velocity]
    std::vector<double> vel_y;
    std::vector<double> vel_z;
    std::vector<double> acc_x;  // Comoving gravitational acceleration
                                // [Code Vel / Code Time]
    std::vector<double> acc_y;
    std::vector<double> acc_z;
    std::vector<double> mass;  // Particle mass [Code Mass]

    // MFM-specific quantities
    std::vector<double> hydro_acc_x;  // Comoving hydrodynamic acceleration
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

    std::vector<double> v_sig_max; // Cached maximum signal velocity

    size_t cooling_failed_cells = 0;
    size_t cooling_total_cycles = 0;
    double accumulated_radiated_energy = 0.0;
    double accumulated_photoheating_energy = 0.0;
    double accumulated_gravitational_work = 0.0;
    double accumulated_expansion_work = 0.0;
    double accumulated_entropy_switch_energy = 0.0;
    double pressure_floor = 0.0;
    size_t ill_conditioned_cases = 0;
    size_t clamped_h_cases = 0;
    size_t non_converged_h_cases = 0;

    // Debugging
    std::vector<double> cond_num;
    std::vector<Eigen::Vector3d> raw_sum_p;
    std::vector<double> n_enc_final;
    double active_particles_fraction = 0.0;
    size_t active_particles_num_cycles = 0;

    // Spatial Hashing
    std::vector<CIC_Data> cic_data;
    double max_h = 0.0;
    double max_accel_sq = 0.0;
    std::vector<uint64_t> morton_codes;
    std::vector<int> sorted_indices;
    std::vector<BVHNode> bvh_nodes;

    // Time tracking per particle (Hierarchical block time-stepping)
    std::vector<int> time_bin;       // Power-of-two bin level (n)
    std::vector<double> dt_step;     // Actual timestep size (Delta t_i)
    std::vector<double> t_current;   // Time this particle was last drifted to
    std::vector<double> t_end;       // Time this particle's current step ends
    std::vector<uint8_t> is_active;  // Particle is synced with the global clock
    std::vector<uint8_t>
        needs_wakeup;  // Thread-safe flag for waking up sleeping particles
    std::vector<int> active_indices;
    size_t num_active = 0;
    double global_time = 0.0;

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
    void compute_and_add_pp_forces(double a, const Config& config,
                                   Diagnostics& diag);
    void compute_cross_pp_forces(double a, ParticleSystem& dm,
                                 const Config& config, Diagnostics& diag);

    const Grid3D& get_rho() const { return gas_rho; }

    void apply_cooling(double dt, double a, const Config& config,
                       Cooling& cooling);

    // Time step calculations
    double get_gravity_timestep(const Config& config) const;
    double get_cfl_timestep(double a, const Config& config) const;
    double get_cooling_timestep(double a, const Config& config,
                                Cooling& cooling) const;
    // Assigns power-of-two block timesteps to all active particles
    void update_particle_timesteps(double dt_max, double a,
                                   const Config& config, Cooling& cooling);
    // Sync the clock, process wakeups, and flag active particles
    void sync_and_activate(double dt, double a, double H, const Config& config);

    // Reversible single-particle kick physics
    KickWork kick_particle_gravity(size_t i, double dt, double a, double H, const Config& config);
    KickWork kick_particle_hydro(size_t i, double dt, const Config& config);

    void hydro_step(const Config& config, double a, double H, double dt);

    // Sync tracked internal energy, total energy, and compute pressure
    void update_primitive_variables(const Config& config, double a);

    // Compute the least-squares matrix gradients for all primitive variables
    void compute_gradients(const Config& config);

    // Solve the Riemann problem between neighbors and update hydro_acc and
    // du_dt
    void compute_hydro_forces(const Config& config, double a, double dt);

    void build_lbvh(const Config& config);

    double get_active_particles_per_cycle_and_reset(const Config& config);

   private:
    void sort_arrays(const std::vector<int>& sorted_indices);

    // Walks a tree to compute the smoothed particle number density and its
    // derivative
    void evaluate_density_sum(size_t particle_idx, double h_guess,
                              double domain_size, double& out_n,
                              double& out_dn_dh) const;

    // A tree-walker that computes the zeta gravity correction from a target
    // tree
    double compute_zeta_contribution(double p1_x, double p1_y, double p1_z,
                                     double h_i,
                                     const std::vector<double>& target_x,
                                     const std::vector<double>& target_y,
                                     const std::vector<double>& target_z,
                                     const std::vector<double>& target_mass,
                                     const std::vector<BVHNode>& target_bvh,
                                     const Config& config) const;
};

// Output of the Riemann Solver
struct MFMFaceFlux {
    Eigen::Vector3d flux_mom;  // Momentum flux density vector (P_star * n)
                               // [PHYSICAL Pressure Units]
    double P_star;  // Resolved face pressure [PHYSICAL Pressure Units]
    double S_star;  // Resolved relative contact wave speed [PHYSICAL Velocity
                    // Units]
};

// Solves the HLLC Riemann problem at the moving face between particles.
// INPUT: Reconstructed face states, relative frame velocity (v_frame), and
// adiabatic index.
// OUTPUT: The resolved pressure (P_star) and relative wave speed (S_star).
// Note that it should get and return physical units
MFMFaceFlux solve_mfm_riemann(const Reconstruction::ReconstructedFace& face,
                              const Eigen::Vector3d& v_frame, double gamma);