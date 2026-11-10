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

class GravitySolver {
#ifdef GRAVITY_GPU_MEMORY
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

    // Arrays requiring CPU interaction
    std::vector<uint8_t> u_is_active;
    std::vector<int> u_active_indices;

    LBVH::Workspace lbvh_ws;

    void free_gpu_resources();
    void prepare_arrays(ParticleSystem& dm, GasParticleSystem& gas,
                        const Config& config, Diagnostics& diag);
#endif

   public:
    GravitySolver();
    ~GravitySolver();

    void update_softening(double a, const Config& config);

    void add_pp_accel(ParticleSystem& dm, GasParticleSystem& gas, double a,
                      const Config& config, Diagnostics& diag);

    void map_results(ParticleSystem& dm, GasParticleSystem& gas,
                     Diagnostics& diag);

    void reset_accels(SimState& state, const Config& config) const;

    // Timestep
    double get_timestep(const ParticleSystem& dm, const GasParticleSystem& gas,
                        const Config& config) const;

    double softening_length;
};
