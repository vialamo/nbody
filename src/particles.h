#pragma once
#include <vector>

#include "config.h"
#include "types.h"

class Diagnostics;
class GasGrid;

struct BoundingBox {
    double min_x, min_y, min_z;
    double max_x, max_y, max_z;
};

struct BVHNode {
    int parent;
    int left_child;
    int right_child;
    int particle_idx;  // Refers to the sorted particle array; -1 if internal
                       // node

    BoundingBox bbox;  // Used for intersection tests
    double max_h;      // Maximum smoothing length in this branch

    // Multipole data for Barnes-Hut Gravity
    double mass;
    double com_x, com_y, com_z;
};

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
};

void compute_and_add_generic_pp_forces(
    size_t n_parts, const double* __restrict__ pos_x,
    const double* __restrict__ pos_y, const double* __restrict__ pos_z,
    const double* __restrict__ mass, double* __restrict__ acc_x,
    double* __restrict__ acc_y, double* __restrict__ acc_z,
    const BVHNode* __restrict__ bvh_nodes, const Config& config,
    Diagnostics& diag);