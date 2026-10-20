#include "pm.h"

#include <omp.h>

#include <vector>

#include "config.h"
#include "pocketfft_hdronly.h"
#include "state.h"

namespace PM {

static void bin_and_assign_mass(const Config& config, size_t num_particles,
                                const std::vector<double>& pos_x,
                                const std::vector<double>& pos_y,
                                const std::vector<double>& pos_z,
                                const std::vector<double>& mass,
                                std::vector<CIC_Data>& out_cic_data,
                                Grid3D& out_rho_grid) {
    if (num_particles == 0) return;

    out_rho_grid.setZero();
    out_cic_data.assign(num_particles, {});

    int Nx = config.mesh_size.x;
    int Ny = config.mesh_size.y;
    int Nz = config.mesh_size.z;
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
        shifted_x =
            fmod(shifted_x + config.domain_size.x, config.domain_size.x);
        shifted_y =
            fmod(shifted_y + config.domain_size.y, config.domain_size.y);
        shifted_z =
            fmod(shifted_z + config.domain_size.z, config.domain_size.z);

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

        int ix0 = (ix + Nx) % Nx, ix1 = (ix + 1 + Nx) % Nx;
        int iy0 = (iy + Ny) % Ny, iy1 = (iy + 1 + Ny) % Ny;
        int iz0 = (iz + Nz) % Nz, iz1 = (iz + 1 + Nz) % Nz;

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

void bin_and_assign_mass(SimState& state, const Config& config) {
    ParticleSystem& dm = state.dm;
    bin_and_assign_mass(config, dm.num_particles, dm.pos_x, dm.pos_y, dm.pos_z,
                        dm.mass, dm.cic_data, dm.grid_rho);

    if (config.hydro_method == HydroMethod::MFM) {
        GasParticleSystem& gas = *state.mfm_gas;
        bin_and_assign_mass(config, gas.num_particles, gas.pos_x, gas.pos_y,
                            gas.pos_z, gas.mass, gas.cic_data, gas.grid_rho);
    }

    if (config.hydro_method == HydroMethod::MFM) {
        state.total_rho.data =
            state.dm.grid_rho.data + state.mfm_gas->grid_rho.data;
    } else if (config.hydro_method == HydroMethod::Eulerian) {
        state.total_rho.data =
            state.dm.grid_rho.data + state.gas->get_density().data;
    } else {
        state.total_rho.data = state.dm.grid_rho.data;
    }
}

static void interpolate_forces(const Config& config, size_t num_particles,
                               const std::vector<CIC_Data>& cic_data,
                               const Grid3D& ax_grid, const Grid3D& ay_grid,
                               const Grid3D& az_grid,
                               std::vector<double>& acc_x,
                               std::vector<double>& acc_y,
                               std::vector<double>& acc_z) {
    if (num_particles == 0) return;

    const int Nx = config.mesh_size.x;
    const int Ny = config.mesh_size.y;
    const int Nz = config.mesh_size.z;

#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < num_particles; ++i) {
        const auto& cd = cic_data[i];

        int ix0 = (cd.ix + Nx) % Nx, ix1 = (cd.ix + 1 + Nx) % Nx;
        int iy0 = (cd.iy + Ny) % Ny, iy1 = (cd.iy + 1 + Ny) % Ny;
        int iz0 = (cd.iz + Nz) % Nz, iz1 = (cd.iz + 1 + Nz) % Nz;

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

void interpolate_forces(SimState& state, const Config& config) {
    ParticleSystem& dm = state.dm;
    interpolate_forces(config, dm.num_particles, dm.cic_data,
                       state.pm_gravity_x, state.pm_gravity_y,
                       state.pm_gravity_z, dm.acc_x, dm.acc_y, dm.acc_z);

    if (config.hydro_method == HydroMethod::MFM) {
        GasParticleSystem& gas = *state.mfm_gas;
        interpolate_forces(config, gas.num_particles, gas.cic_data,
                           state.pm_gravity_x, state.pm_gravity_y,
                           state.pm_gravity_z, gas.acc_x, gas.acc_y, gas.acc_z);
    }
}

void compute_acceleration(SimState& state, const Config& config) {
    int Nx = config.mesh_size.x;
    int Ny = config.mesh_size.y;
    int Nz = config.mesh_size.z;

    const GasGrid& gas = *state.gas;
    const Grid3D& dm_rho = state.dm.grid_rho;
    Grid3D& total_rho = state.total_rho;
    Grid3D& acc_x = state.pm_gravity_x;
    Grid3D& acc_y = state.pm_gravity_y;
    Grid3D& acc_z = state.pm_gravity_z;
    Grid3D& phi = state.phi;

    if (config.hydro_method == HydroMethod::MFM) {
        total_rho.data = dm_rho.data + state.mfm_gas->grid_rho.data;
    } else if (config.hydro_method == HydroMethod::Eulerian) {
        total_rho.data = dm_rho.data + gas.get_density().data;
    } else {
        total_rho.data = dm_rho.data;
    }

    // Anisotropic PocketFFT Shape and Strides
    const pocketfft::shape_t shape = {(size_t)Nx, (size_t)Ny, (size_t)Nz};

    const pocketfft::stride_t stride_r = {
        static_cast<ptrdiff_t>((size_t)Ny * Nz * sizeof(double)),
        static_cast<ptrdiff_t>(Nz * sizeof(double)), sizeof(double)};

    // The Z-axis is cut roughly in half for Real-to-Complex transforms
    const pocketfft::stride_t stride_c = {
        static_cast<ptrdiff_t>((size_t)Ny * (Nz / 2 + 1) *
                               sizeof(std::complex<double>)),
        static_cast<ptrdiff_t>((Nz / 2 + 1) * sizeof(std::complex<double>)),
        sizeof(std::complex<double>)};

    size_t complex_size = (size_t)Nx * Ny * (Nz / 2 + 1);
    std::vector<std::complex<double>> rho_k(complex_size);

    pocketfft::r2c(shape, stride_r, stride_c, {0, 1, 2}, true,
                   total_rho.raw_data(), rho_k.data(), 1.0, 0);

    std::vector<std::complex<double>> phi_k(complex_size);

    const double r_s = config.PM_smoothing_cells * config.cell_size;

#pragma omp parallel for collapse(3)
    for (int i = 0; i < Nx; ++i) {
        for (int j = 0; j < Ny; ++j) {
            for (int k = 0; k < Nz / 2 + 1; ++k) {
                if (i == 0 && j == 0 && k == 0) {
                    phi_k[0] = {0.0, 0.0};
                    continue;
                }

                // Frequencies rely on their specific axis sizes
                double kx_freq =
                    static_cast<double>((i < Nx / 2) ? i : (i - Nx));
                double ky_freq =
                    static_cast<double>((j < Ny / 2) ? j : (j - Ny));
                double kz_freq = static_cast<double>(k);

                // K-vectors rely on their specific physical domain sizes
                double kx = kx_freq * 2.0 * M_PI / config.domain_size.x;
                double ky = ky_freq * 2.0 * M_PI / config.domain_size.y;
                double kz = kz_freq * 2.0 * M_PI / config.domain_size.z;
                double k2 = kx * kx + ky * ky + kz * kz;

                size_t idx = static_cast<size_t>(i) * Ny * (Nz / 2 + 1) +
                             static_cast<size_t>(j) * (Nz / 2 + 1) +
                             static_cast<size_t>(k);

                double filter = std::exp(-k2 * r_s * r_s);

                auto sinc = [](double freq_index, int mesh_size) {
                    if (freq_index == 0.0) return 1.0;
                    double arg = M_PI * freq_index / mesh_size;
                    return std::sin(arg) / arg;
                };

                // CIC filter deconvolves using the explicit mesh sizes
                double sinc_x = sinc(kx_freq, Nx);
                double sinc_y = sinc(ky_freq, Ny);
                double sinc_z = sinc(kz_freq, Nz);
                double cic_blur = std::pow(sinc_x * sinc_y * sinc_z, 4);

                phi_k[idx] = rho_k[idx] * (-4.0 * M_PI * config.G / k2) *
                             (filter / cic_blur);
            }
        }
    }

    pocketfft::c2r(shape, stride_c, stride_r, {0, 1, 2}, false, phi_k.data(),
                   const_cast<double*>(phi.raw_data()), 1.0, 0);

    // Normalization and gradients
    double norm = 1.0 / ((double)Nx * Ny * Nz);
    double a_inv3 =
        1.0 / (state.scale_factor * state.scale_factor * state.scale_factor);
    double factor = (-norm / (2.0 * config.cell_size)) * a_inv3;

#pragma omp parallel for collapse(2)
    for (int i = 0; i < Nx; ++i) {
        for (int j = 0; j < Ny; ++j) {
            int ip1 = (i + 1 == Nx) ? 0 : i + 1;
            int im1 = (i == 0) ? Nx - 1 : i - 1;
            int jp1 = (j + 1 == Ny) ? 0 : j + 1;
            int jm1 = (j == 0) ? Ny - 1 : j - 1;

#pragma omp simd
            for (int k = 0; k < Nz; ++k) {
                int kp1 = (k + 1 == Nz) ? 0 : k + 1;
                int km1 = (k == 0) ? Nz - 1 : k - 1;

                acc_x(i, j, k) = (phi(ip1, j, k) - phi(im1, j, k)) * factor;
                acc_y(i, j, k) = (phi(i, jp1, k) - phi(i, jm1, k)) * factor;
                acc_z(i, j, k) = (phi(i, j, kp1) - phi(i, j, km1)) * factor;
            }
        }
    }
}

}  // namespace PM
