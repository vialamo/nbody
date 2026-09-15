#pragma once
#include <vector>

#include "config.h"
#include "types.h"

namespace CIC {

// Splats particle masses onto a 3D grid and caches the interpolation weights
void bin_and_assign_mass(const Config& config, size_t num_particles,
                         const std::vector<double>& pos_x,
                         const std::vector<double>& pos_y,
                         const std::vector<double>& pos_z,
                         const std::vector<double>& mass,
                         std::vector<CIC_Data>& out_cic_data,
                         Grid3D& out_rho_grid);

// Gathers forces from 3D grids back to the particles using cached weights
void interpolate_forces(const Config& config, size_t num_particles,
                        const std::vector<CIC_Data>& cic_data,
                        const Grid3D& ax_grid, const Grid3D& ay_grid,
                        const Grid3D& az_grid, std::vector<double>& acc_x,
                        std::vector<double>& acc_y, std::vector<double>& acc_z);

}  // namespace CIC