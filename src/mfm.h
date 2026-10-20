#pragma once
#include <Eigen/Dense>
#include <vector>

#include "lbvh.h"
#include "types.h"

constexpr double N_cond_crit = 5.0;
constexpr double density_floor = 1e-12;

// #define USE_ADAPTIVE_SOFTENING
// #define USE_GEOMETRIC_MATRIX

struct Config;
class Cooling;

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

    std::vector<double> cond_num;
    std::vector<double> cond_num_convergence;

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

#ifdef USE_ADAPTIVE_SOFTENING
    std::vector<double> zeta;  // Correction term for adaptive gravity softening
#endif

    Grid3D grid_rho;  // Gridded comoving gas density for PM gravity/diagnostics
                     // [Code Mass / Code Length^3]

    std::vector<double> v_sig_max;  // Cached maximum signal velocity

    size_t cooling_failed_cells = 0;
    size_t cooling_total_cycles = 0;
    double accumulated_radiated_energy = 0.0;
    double accumulated_photoheating_energy = 0.0;
    double accumulated_gravitational_work = 0.0;
    double accumulated_expansion_work = 0.0;
    double accumulated_entropy_switch_energy = 0.0;
    double pressure_floor = 0.0;

    // Debugging
    std::vector<double> n_enc_final;
    double active_particles_fraction = 0.0;
    size_t active_particles_num_cycles = 0;
    size_t ill_conditioned_cases = 0;
    size_t clamped_h_cases = 0;
    size_t non_converged_h_cases = 0;
    size_t not_enough_neighbors_cases = 0;
    size_t riemann_failed_cases = 0;
    size_t neighbor_increased_cases = 0;
    size_t sph_fallback_cases = 0;

    // Spatial Hashing
    std::vector<CIC_Data> cic_data;
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

    std::vector<double> vel_x_pred;
    std::vector<double> vel_y_pred;
    std::vector<double> vel_z_pred;
    std::vector<double> pressure_pred;

    GasParticleSystem(const Config& config);

    void add_particle(double px, double py, double pz, double vx, double vy,
                      double vz, double m, double initial_u, double initial_h,
                      double z_metal);
    
    void build_lbvh(const Config& config);

    // Assigns power-of-two block timesteps to all active particles
    void update_particle_timesteps(double dt_max, double a,
                                   const Config& config, Cooling& cooling);
                                   
    // Sync the clock, process wakeups, and flag active particles
    void sync_and_activate(double dt, double a, double H, const Config& config);
    
    void hydro_step(const Config& config, double a);

    // Time step calculations
    double get_cfl_timestep(double a, const Config& config) const;
    

   private:
    void sort_arrays(const std::vector<int>& sorted_indices);


};
