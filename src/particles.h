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

    double accumulated_gravitational_work = 0.0;
    double accumulated_expansion_work = 0.0;

    std::vector<uint8_t> is_active;
    std::vector<int> time_bin;
    std::vector<double> dt_step;
    std::vector<double> t_current;
    std::vector<double> t_end;
    std::vector<int> active_indices;
    size_t num_active = 0;
    double global_time = 0.0;

    ParticleSystem(const Config& config);

    void build_lbvh(const Config& config);

    void add_particle(double px, double py, double pz, double vx, double vy,
                      double vz, double m);

    void update_particle_timesteps(double dt_max, const Config& config);
    void sync_and_activate(double dt, const Config& config);

   private:
    void sort_arrays(const std::vector<int>& sorted_indices);
};
