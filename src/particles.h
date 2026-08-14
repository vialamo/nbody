#pragma once
#include <vector>

#include "config.h"
#include "types.h"

class Diagnostics;
class GasGrid;

// Encapsulates a flat-array for particle binning
struct CellList {
    std::vector<int> cell_start;  // Size: num_cells
    std::vector<int> cell_count;  // Size: num_cells

    void resize(int num_cells, int num_parts) {
        cell_start.assign(num_cells, 0);
        cell_count.assign(num_cells, 0);
    }
};

class ParticleSystem {
   private:
    std::vector<CIC_Data> cic_data;

    friend struct ParticleTestAccess;

   public:
    CellList cell_list;
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
    size_t n_parts,
    const double* __restrict__ pos_x,
    const double* __restrict__ pos_y,
    const double* __restrict__ pos_z,
    const double* __restrict__ mass,
    double* __restrict__ acc_x,
    double* __restrict__ acc_y,
    double* __restrict__ acc_z,
    const int* __restrict__ cell_start,
    const int* __restrict__ cell_count,
    const Config& config,
    Diagnostics& diag);