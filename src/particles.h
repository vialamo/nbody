#pragma once
#include <vector>

#include "lbvh.h"
#include "types.h"

struct Config;

class ParticleSystem {
   private:
    friend struct ParticleTestAccess;

   public:
    std::vector<CIC_Data> cic_data;

    std::vector<uint64_t> morton_codes;
    std::vector<int> sorted_indices;
    std::vector<BVHNode> bvh_nodes;

    size_t num_particles = 0;

    std::vector<double> pos_x;
    std::vector<double> pos_y;
    std::vector<double> pos_z;

    std::vector<double> vel_x;
    std::vector<double> vel_y;
    std::vector<double> vel_z;

    std::vector<double> acc_x;  // Comoving gravitational acceleration
                                // [Code Vel / Code Time]
    std::vector<double> acc_y;
    std::vector<double> acc_z;

    std::vector<double> mass;

    EnergyBudget budget;

    uint64_t ti_global_time = 0;
    std::vector<uint64_t> ti_step;     // Duration of the current step
    std::vector<uint64_t> ti_current;  // Time this particle was last drifted to
    std::vector<uint64_t> ti_end;      // Time this particle's current step ends

    std::vector<uint8_t> is_active;
    std::vector<int> active_indices;
    size_t get_num_active() const { return active_indices.size(); }

    ParticleSystem(const Config& config);

    void build_lbvh(const Config& config);

    void add_particle(double px, double py, double pz, double vx, double vy,
                      double vz, double m);

    void update_particle_timesteps(const Config& config,
                                   double softening_length,
                                   double hubble_param);
    void sync_and_activate(uint64_t ti_dt, const Config& config);

   private:
    void sort_arrays(const std::vector<int>& sorted_indices);

    // Reusable scratch buffers for sort_arrays()
    std::vector<double> scratch_double;
    std::vector<uint8_t> scratch_uint8;
    std::vector<uint64_t> scratch_uint64;
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
