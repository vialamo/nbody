#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

#include "lbvh.h"

struct BVHNode;
struct Config;
struct Diagnostics;
struct SimState;
class GasParticleSystem;
class ParticleSystem;
class StarParticleSystem;

class GravitySolver {
#ifdef USE_GPU
   private:
    size_t current_gpu_capacity = 0;
    size_t current_bvh_capacity = 0;
    double *d_u_px = nullptr, *d_u_py = nullptr, *d_u_pz = nullptr,
           *d_u_m = nullptr;
    double *d_t_px = nullptr, *d_t_py = nullptr, *d_t_pz = nullptr,
           *d_t_m = nullptr;
    double *d_u_ax = nullptr, *d_u_ay = nullptr, *d_u_az = nullptr;
    int *d_u_orig_idx = nullptr, *d_t_orig_idx = nullptr,
        *d_u_sorted_idx = nullptr;
    uint8_t *d_u_type = nullptr, *d_t_type = nullptr, *d_t_is_active = nullptr;
    uint64_t* d_u_morton = nullptr;
    BVHNode* d_u_bvh = nullptr;

    // Adaptive softening
    double *d_u_eps = nullptr, *d_u_gradh_fac = nullptr;
    double *d_t_eps = nullptr, *d_t_gradh_fac = nullptr;
    double* d_u_xi = nullptr;

    // Arrays requiring CPU interaction
    std::vector<uint8_t> u_is_active;
    std::vector<int> u_active_indices;

    LBVH::Workspace lbvh_ws;

    void prepare_arrays(ParticleSystem& dm, GasParticleSystem& gas,
                        StarParticleSystem& star, const Config& config,
                        Diagnostics& diag);

    void resize_gpu_buffers(size_t required_capacity);

    void upload_species_to_gpu(size_t n, size_t offset, uint8_t type,
                               const double* px, const double* py,
                               const double* pz, const double* m,
                               const uint8_t* is_active, const double* h,
                               const double* rho, const double* omega,
                               bool use_adaptive, double eps_min,
                               double eps_max);

    void download_species_accels(size_t n, uint8_t type, double* ax, double* ay,
                                 double* az, const int* d_active_idx,
                                 size_t n_active);
#endif

   public:
    GravitySolver();
    ~GravitySolver();

    GravitySolver(const GravitySolver&) = delete;
    GravitySolver& operator=(const GravitySolver&) = delete;

    void compute_dm_density_and_h(ParticleSystem& dm, const Config& config);
    void update_softening(double a, const Config& config);

    void add_pp_accel(ParticleSystem& dm, GasParticleSystem& gas,
                      StarParticleSystem& star, double a, const Config& config,
                      Diagnostics& diag);

    void map_results(ParticleSystem& dm, GasParticleSystem& gas,
                     StarParticleSystem& star, Diagnostics& diag);

    void reset_accels(SimState& state, const Config& config) const;

    double get_timestep(const ParticleSystem& dm, const GasParticleSystem& gas,
                        const Config& config) const;

    double softening_length;
};
