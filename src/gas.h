#pragma once
#include <Eigen/Dense>
#include <vector>

#include "lbvh.h"
#include "types.h"

struct Config;
class Integrator;

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

    std::vector<double>
        entropy;  // Entropic function (S = P / rho^gamma)
                  // [Code Velocity^2 / (Code Mass / Code Length^3)^(gamma-1)]
    std::vector<double> max_rel_ke;  // Maximum specific relative kinetic energy
                                     // of neighbors [Code Velocity^2]

    // Gradients are evaluated with respect to comoving coordinates
    // (d/dx_comoving)
    std::vector<Eigen::Vector3d> grad_rho;  // [Density / Code Length]
    std::vector<Eigen::Vector3d> grad_vx;   // [Velocity / Code Length]
    std::vector<Eigen::Vector3d> grad_vy;
    std::vector<Eigen::Vector3d> grad_vz;
    std::vector<Eigen::Vector3d> grad_p;  // [Pressure / Code Length]

    std::vector<double> v_sig_max;  // Cached max signal velocity (code units)

    EnergyBudget budget;

    // Spatial Hashing
    std::vector<CIC_Data> cic_data;
    std::vector<uint64_t> morton_codes;
    std::vector<int> sorted_indices;
    std::vector<BVHNode> bvh_nodes;

    // Time tracking per particle (Hierarchical block time-stepping)
    uint64_t ti_global_time = 0;
    std::vector<uint64_t> ti_step;     // Duration of the current step
    std::vector<uint64_t> ti_current;  // Time this particle was last drifted to
    std::vector<uint64_t> ti_end;      // Time this particle's current step ends
    std::vector<uint8_t> changed_timestep;  // Particle has changed timestep

    std::vector<uint8_t> is_active;  // Particle is synced with the global clock
    std::vector<int> active_indices;
    size_t get_num_active() const { return active_indices.size(); }
    std::vector<uint8_t> needs_wakeup;

    std::vector<double> vel_x_pred;
    std::vector<double> vel_y_pred;
    std::vector<double> vel_z_pred;
    std::vector<double> pressure_pred;

    std::vector<double> hydro_dv_x;
    std::vector<double> hydro_dv_y;
    std::vector<double> hydro_dv_z;
    std::vector<double> hydro_du;
    std::vector<double> hydro_de;

    std::vector<double> snap_vel_x;
    std::vector<double> snap_vel_y;
    std::vector<double> snap_vel_z;
    std::vector<double> snap_u;
    std::vector<double> snap_total_energy;
    std::vector<double> snap_pressure;
    std::vector<double> snap_entropy;

    // Debugging
    std::vector<double> n_enc_final;

    struct GasDiagnostics {
        double active_particles_fraction = 0.0;
        size_t active_particles_num_cycles = 0;
        size_t ill_conditioned = 0;
        size_t clamped_h = 0;
        size_t non_converged_h = 0;
        size_t not_enough_neighbors = 0;
        size_t riemann_failed = 0;
        size_t neighbor_increased = 0;
        size_t sph_fallback = 0;
        size_t cooling_failed = 0;
        size_t particles_at_min_dt = 0;
        size_t forced_wakeup = 0;
    };
    GasDiagnostics diag;

    GasParticleSystem(const Config& config);

    void add_particle(double px, double py, double pz, double vx, double vy,
                      double vz, double m, double initial_u, double initial_h,
                      double z_metal);

    void build_lbvh(const Config& config);

    // Assigns power-of-two block timesteps to all active particles
    void update_particle_timesteps(const Config& config,
                                   double softening_length,
                                   double hubble_param);

    // Sync the clock, process wakeups, and flag active particles
    void sync_and_activate(Integrator& integrator, uint64_t ti_dt, double a,
                           double H, const Config& config);

    void reset_diagnostic_counters();
    void reset_accumulators();

   private:
    void enforce_timestep_restrictions(uint64_t ti_dt, const Config& config);
    void sort_arrays(const std::vector<int>& sorted_indices);

    // Reusable scratch buffers for sort_arrays()
    std::vector<double> scratch_double;
    std::vector<uint8_t> scratch_uint8;
    std::vector<uint64_t> scratch_uint64;
    std::vector<Eigen::Vector3d> scratch_vec3;
    std::vector<Eigen::Matrix3d> scratch_mat3;
    std::vector<CIC_Data> scratch_cic;

    template <typename T>
    static void permute_array(std::vector<T>& arr, std::vector<T>& scratch,
                              const std::vector<int>& sorted_indices) {
        const size_t n = arr.size();
        if (scratch.size() < n) {
            scratch.resize(n);
        }

#pragma omp parallel for schedule(static)
        for (size_t i = 0; i < n; ++i) {
            scratch[i] = arr[sorted_indices[i]];
        }
        arr.swap(scratch);
    }
};
