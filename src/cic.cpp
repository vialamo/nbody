#include "cic.h"

#include <omp.h>

#include <cmath>

namespace CIC {

void bin_and_assign_mass(const Config& config, size_t num_particles,
                         const std::vector<double>& pos_x,
                         const std::vector<double>& pos_y,
                         const std::vector<double>& pos_z,
                         const std::vector<double>& mass,
                         std::vector<CIC_Data>& out_cic_data,
                         Grid3D& out_rho_grid) {
    if (num_particles == 0) return;

    out_rho_grid.setZero();
    out_cic_data.assign(num_particles, {});

    int N = config.mesh_size;
    double inv_cell_size = 1.0 / config.cell_size;

    // Calculate cells & densities
    // Note: This loop is intentionally sequential to avoid OpenMP
    // race conditions when multiple particles write to the same grid cell
    for (size_t i = 0; i < num_particles; ++i) {
        double px = pos_x[i], py = pos_y[i], pz = pos_z[i];

        // Cell centered PM grid nodes
        double shifted_x = px - 0.5 * config.cell_size;
        double shifted_y = py - 0.5 * config.cell_size;
        double shifted_z = pz - 0.5 * config.cell_size;

        // Ensure periodic wrap-around bounds[cite: 13, 15]
        shifted_x = fmod(shifted_x + config.domain_size, config.domain_size);
        shifted_y = fmod(shifted_y + config.domain_size, config.domain_size);
        shifted_z = fmod(shifted_z + config.domain_size, config.domain_size);

        int ix = static_cast<int>(shifted_x * inv_cell_size);
        int iy = static_cast<int>(shifted_y * inv_cell_size);
        int iz = static_cast<int>(shifted_z * inv_cell_size);

        double frac_x = (shifted_x * inv_cell_size) - ix;
        double frac_y = (shifted_y * inv_cell_size) - iy;
        double frac_z = (shifted_z * inv_cell_size) - iz;

        // Trilinear interpolation weights[cite: 13, 15]
        double w000 = (1 - frac_x) * (1 - frac_y) * (1 - frac_z);
        double w100 = frac_x * (1 - frac_y) * (1 - frac_z);
        double w010 = (1 - frac_x) * frac_y * (1 - frac_z);
        double w110 = frac_x * frac_y * (1 - frac_z);
        double w001 = (1 - frac_x) * (1 - frac_y) * frac_z;
        double w101 = frac_x * (1 - frac_y) * frac_z;
        double w011 = (1 - frac_x) * frac_y * frac_z;
        double w111 = frac_x * frac_y * frac_z;

        out_cic_data[i] = {ix,   iy,   iz,   w000, w100, w010,
                           w110, w001, w101, w011, w111};

        int ix0 = (ix + N) % N, ix1 = (ix + 1 + N) % N;
        int iy0 = (iy + N) % N, iy1 = (iy + 1 + N) % N;
        int iz0 = (iz + N) % N, iz1 = (iz + 1 + N) % N;

        double m = mass[i];
        out_rho_grid(ix0, iy0, iz0) += m * w000;
        out_rho_grid(ix1, iy0, iz0) += m * w100;
        out_rho_grid(ix0, iy1, iz0) += m * w010;
        out_rho_grid(ix1, iy1, iz0) += m * w110;
        out_rho_grid(ix0, iy0, iz1) += m * w001;
        out_rho_grid(ix1, iy0, iz1) += m * w101;
        out_rho_grid(ix0, iy1, iz1) += m * w011;
        out_rho_grid(ix1, iy1, iz1) += m * w111;
    }

    out_rho_grid.data /= config.cell_volume;
}

void interpolate_forces(const Config& config, size_t num_particles,
                        const std::vector<CIC_Data>& cic_data,
                        const Grid3D& ax_grid, const Grid3D& ay_grid,
                        const Grid3D& az_grid, std::vector<double>& acc_x,
                        std::vector<double>& acc_y,
                        std::vector<double>& acc_z) {
    if (num_particles == 0) return;

    const int N = config.mesh_size;

#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < num_particles; ++i) {
        const auto& cd = cic_data[i];

        int ix0 = (cd.ix + N) % N, ix1 = (cd.ix + 1 + N) % N;
        int iy0 = (cd.iy + N) % N, iy1 = (cd.iy + 1 + N) % N;
        int iz0 = (cd.iz + N) % N, iz1 = (cd.iz + 1 + N) % N;

        // Gather force from the 8 intersecting grid nodes
        auto interp = [&](const Grid3D& grid) {
            return grid(ix0, iy0, iz0) * cd.w000 +
                   grid(ix1, iy0, iz0) * cd.w100 +
                   grid(ix0, iy1, iz0) * cd.w010 +
                   grid(ix1, iy1, iz0) * cd.w110 +
                   grid(ix0, iy0, iz1) * cd.w001 +
                   grid(ix1, iy0, iz1) * cd.w101 +
                   grid(ix0, iy1, iz1) * cd.w011 +
                   grid(ix1, iy1, iz1) * cd.w111;
        };

        acc_x[i] = interp(ax_grid);
        acc_y[i] = interp(ay_grid);
        acc_z[i] = interp(az_grid);
    }
}

}  // namespace CIC