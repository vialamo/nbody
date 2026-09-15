#pragma once
#include <vector>

#include "config.h"
#include "types.h"
#include "lbvh.h"

class Diagnostics;
class GasGrid;

class ParticleSystem {
   private:
    std::vector<CIC_Data> cic_data;

    friend struct ParticleTestAccess;

   public:
    std::vector<uint64_t> morton_codes;
    std::vector<int> sorted_indices;
    std::vector<BVHNode> bvh_nodes;

    Grid3D dm_rho;

    size_t num_particles = 0;

    std::vector<double> pos_x;
    std::vector<double> pos_y;
    std::vector<double> pos_z;

    std::vector<double> vel_x;
    std::vector<double> vel_y;
    std::vector<double> vel_z;

    std::vector<double> acc_x;
    std::vector<double> acc_y;
    std::vector<double> acc_z;

    std::vector<double> mass;

    double max_accel_sq;

    double accumulated_gravitational_work = 0.0;
    double accumulated_expansion_work = 0.0;

    ParticleSystem(const Config& config);

    void build_lbvh(const Config& config);
    void bin_and_assign_mass(const Config& config);

    void interpolate_cic_forces(const Grid3D& ax_grid, const Grid3D& ay_grid,
                                const Grid3D& az_grid, const Config& config);

    void compute_and_add_pp_forces(const Config& config, Diagnostics& diag);

    void compute_gas_dm_pp_forces(const GasGrid& gas, Grid3D& grav_x,
                                  Grid3D& grav_y, Grid3D& grav_z,
                                  const Config& config, Diagnostics& diag);

    double get_gravity_timestep(const Config& config) const;

    const Grid3D& get_rho() const { return dm_rho; }

    void add_particle(double px, double py, double pz, double vx, double vy,
                      double vz, double m);

   private:
    void sort_arrays(const std::vector<int>& sorted_indices);
};