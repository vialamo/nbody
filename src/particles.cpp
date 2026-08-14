#include "particles.h"

#include <omp.h>

#include <cmath>
#include <limits>

#include "diagnostics.h"
#include "gas.h"
#include "math_utils.h"

ParticleSystem::ParticleSystem(const Config& config)
    : dm_rho(config.mesh_size),
      cic_data(config.num_dm_particles),
      max_accel_sq(0.0) {
    pos_x.reserve(config.num_dm_particles);
    pos_y.reserve(config.num_dm_particles);
    pos_z.reserve(config.num_dm_particles);

    vel_x.reserve(config.num_dm_particles);
    vel_y.reserve(config.num_dm_particles);
    vel_z.reserve(config.num_dm_particles);

    acc_x.reserve(config.num_dm_particles);
    acc_y.reserve(config.num_dm_particles);
    acc_z.reserve(config.num_dm_particles);

    mass.reserve(config.num_dm_particles);

    int num_cells = config.mesh_size * config.mesh_size * config.mesh_size;
    cell_list.resize(num_cells, config.num_dm_particles);
}

void ParticleSystem::add_particle(double px, double py, double pz, double vx,
                                  double vy, double vz, double m) {
    pos_x.push_back(px);
    pos_y.push_back(py);
    pos_z.push_back(pz);
    vel_x.push_back(vx);
    vel_y.push_back(vy);
    vel_z.push_back(vz);
    acc_x.push_back(0.0);
    acc_y.push_back(0.0);
    acc_z.push_back(0.0);
    mass.push_back(m);
    num_particles++;
}

void ParticleSystem::bin_and_assign_mass(const Config& config) {
    dm_rho.setZero();
    cic_data.assign(num_particles, {});

    int N = config.mesh_size;
    int num_cells = N * N * N;

    std::fill(cell_list.cell_count.begin(), cell_list.cell_count.end(), 0);
    std::fill(cell_list.cell_start.begin(), cell_list.cell_start.end(), 0);

    std::vector<int> particle_cell_idx(num_particles);

    // Calculate cells & densities
    for (size_t i = 0; i < num_particles; ++i) {
        double px = pos_x[i], py = pos_y[i], pz = pos_z[i];

        // Must use physical coordinates to align with the PP Gravity solver
        int hash_ix = static_cast<int>(px / config.cell_size) % N;
        int hash_iy = static_cast<int>(py / config.cell_size) % N;
        int hash_iz = static_cast<int>(pz / config.cell_size) % N;

        // Protect against negative bounds
        hash_ix = (hash_ix + N) % N;
        hash_iy = (hash_iy + N) % N;
        hash_iz = (hash_iz + N) % N;
        int cell_index = hash_iz * N * N + hash_iy * N + hash_ix;

        // Cell centered PM grid nodes
        double shifted_x = px - 0.5 * config.cell_size;
        double shifted_y = py - 0.5 * config.cell_size;
        double shifted_z = pz - 0.5 * config.cell_size;

        // Ensure periodic wrap-around bounds
        shifted_x = fmod(shifted_x + config.domain_size, config.domain_size);
        shifted_y = fmod(shifted_y + config.domain_size, config.domain_size);
        shifted_z = fmod(shifted_z + config.domain_size, config.domain_size);

        int ix = static_cast<int>(shifted_x / config.cell_size);
        int iy = static_cast<int>(shifted_y / config.cell_size);
        int iz = static_cast<int>(shifted_z / config.cell_size);

        double frac_x = (shifted_x / config.cell_size) - ix;
        double frac_y = (shifted_y / config.cell_size) - iy;
        double frac_z = (shifted_z / config.cell_size) - iz;

        double w000 = (1 - frac_x) * (1 - frac_y) * (1 - frac_z);
        double w100 = frac_x * (1 - frac_y) * (1 - frac_z);
        double w010 = (1 - frac_x) * frac_y * (1 - frac_z);
        double w110 = frac_x * frac_y * (1 - frac_z);
        double w001 = (1 - frac_x) * (1 - frac_y) * frac_z;
        double w101 = frac_x * (1 - frac_y) * frac_z;
        double w011 = (1 - frac_x) * frac_y * frac_z;
        double w111 = frac_x * frac_y * frac_z;

        cic_data[i] = {ix,   iy,   iz,   w000, w100, w010,
                       w110, w001, w101, w011, w111};

        int ix0 = (ix + N) % N, ix1 = (ix + 1 + N) % N;
        int iy0 = (iy + N) % N, iy1 = (iy + 1 + N) % N;
        int iz0 = (iz + N) % N, iz1 = (iz + 1 + N) % N;

        double m = mass[i];
        dm_rho(ix0, iy0, iz0) += m * w000;
        dm_rho(ix1, iy0, iz0) += m * w100;
        dm_rho(ix0, iy1, iz0) += m * w010;
        dm_rho(ix1, iy1, iz0) += m * w110;
        dm_rho(ix0, iy0, iz1) += m * w001;
        dm_rho(ix1, iy0, iz1) += m * w101;
        dm_rho(ix0, iy1, iz1) += m * w011;
        dm_rho(ix1, iy1, iz1) += m * w111;

        particle_cell_idx[i] = cell_index;
        cell_list.cell_count[cell_index]++;
    }

    dm_rho.data /= config.cell_volume;

    // Prefix sum: calculate offsets where each cell's block of particles
    // will start in the sorted array
    int current_offset = 0;
    for (int c = 0; c < num_cells; ++c) {
        cell_list.cell_start[c] = current_offset;
        current_offset += cell_list.cell_count[c];
    }

    // Sort the arrays
    std::vector<double> new_px(num_particles), new_py(num_particles),
        new_pz(num_particles);
    std::vector<double> new_vx(num_particles), new_vy(num_particles),
        new_vz(num_particles);
    std::vector<double> new_ax(num_particles), new_ay(num_particles),
        new_az(num_particles);
    std::vector<double> new_m(num_particles);
    std::vector<CIC_Data> new_cic(num_particles);

    std::vector<int> write_offset = cell_list.cell_start;
    for (size_t i = 0; i < num_particles; ++i) {
        int dest = write_offset[particle_cell_idx[i]]++;
        new_px[dest] = pos_x[i];
        new_py[dest] = pos_y[i];
        new_pz[dest] = pos_z[i];
        new_vx[dest] = vel_x[i];
        new_vy[dest] = vel_y[i];
        new_vz[dest] = vel_z[i];
        new_ax[dest] = acc_x[i];
        new_ay[dest] = acc_y[i];
        new_az[dest] = acc_z[i];
        new_m[dest] = mass[i];
        new_cic[dest] = cic_data[i];
    }

    // Move data
    pos_x = std::move(new_px);
    pos_y = std::move(new_py);
    pos_z = std::move(new_pz);
    vel_x = std::move(new_vx);
    vel_y = std::move(new_vy);
    vel_z = std::move(new_vz);
    acc_x = std::move(new_ax);
    acc_y = std::move(new_ay);
    acc_z = std::move(new_az);
    mass = std::move(new_m);
    cic_data = std::move(new_cic);
}

void ParticleSystem::interpolate_cic_forces(const Grid3D& ax_grid,
                                            const Grid3D& ay_grid,
                                            const Grid3D& az_grid,
                                            const Config& config) {
    const int N = config.mesh_size;

#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < num_particles; ++i) {
        const auto& cd = cic_data[i];

        int ix0 = (cd.ix + N) % N, ix1 = (cd.ix + 1 + N) % N;
        int iy0 = (cd.iy + N) % N, iy1 = (cd.iy + 1 + N) % N;
        int iz0 = (cd.iz + N) % N, iz1 = (cd.iz + 1 + N) % N;

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

void ParticleSystem::compute_and_add_pp_forces(const Config& config,
                                               Diagnostics& diag) {
    compute_and_add_generic_pp_forces(
        num_particles, pos_x.data(), pos_y.data(), pos_z.data(), mass.data(),
        acc_x.data(), acc_y.data(), acc_z.data(), cell_list.cell_start.data(),
        cell_list.cell_count.data(), config, diag);
}

void ParticleSystem::compute_gas_dm_pp_forces(const GasGrid& gas,
                                              Grid3D& grav_x, Grid3D& grav_y,
                                              Grid3D& grav_z,
                                              const Config& config,
                                              Diagnostics& diag) {
    if (config.hydro_method != HydroMethod::Eulerian) return;

    // Shared setup variables
    int search_radius =
        config.use_PM
            ? static_cast<int>(ceil(config.cutoff_radius / config.cell_size))
            : config.mesh_size / 2;
    const int N = config.mesh_size;
    const double domain_size = config.domain_size;
    const double G = config.G;
    const double cutoff_sq = config.cutoff_radius_squared;
    const double r_s = config.PM_smoothing_cells * config.cell_size;
    const bool use_pm = config.use_PM;
    const size_t n_parts = num_particles;

    const double cell_vol = config.cell_volume;
    const double cell_size = config.cell_size;

#ifdef USE_GPU
    if (config.enable_GPU) {
        // ========================================================================
        // GPU IMPLEMENTATION
        // ========================================================================
        double* d_px = pos_x.data();
        double* d_py = pos_y.data();
        double* d_pz = pos_z.data();
        double* d_m = mass.data();
        double* d_ax = acc_x.data();
        double* d_ay = acc_y.data();
        double* d_az = acc_z.data();

        const double* d_gas_rho = gas.get_density().raw_data();
        double* d_grav_x = grav_x.data.data();
        double* d_grav_y = grav_y.data.data();
        double* d_grav_z = grav_z.data.data();

#pragma omp target enter data map(to : d_px[0 : n_parts], d_py[0 : n_parts], \
                                      d_pz[0 : n_parts], d_m[0 : n_parts],   \
                                      d_ax[0 : n_parts], d_ay[0 : n_parts],  \
                                      d_az[0 : n_parts])

#pragma omp target teams distribute parallel for
        for (size_t i = 0; i < n_parts; ++i) {
            double p1_x = d_px[i], p1_y = d_py[i], p1_z = d_pz[i];
            int ix = static_cast<int>(p1_x / cell_size);
            int iy = static_cast<int>(p1_y / cell_size);
            int iz = static_cast<int>(p1_z / cell_size);

            double local_acc_x = 0.0, local_acc_y = 0.0, local_acc_z = 0.0;
            double dm_mass = d_m[i];

            for (int dx_cell = -search_radius; dx_cell <= search_radius;
                 ++dx_cell) {
                for (int dy_cell = -search_radius; dy_cell <= search_radius;
                     ++dy_cell) {
                    for (int dz_cell = -search_radius; dz_cell <= search_radius;
                         ++dz_cell) {
                        int neighbor_ix = (((ix + dx_cell) % N) + N) % N;
                        int neighbor_iy = (((iy + dy_cell) % N) + N) % N;
                        int neighbor_iz = (((iz + dz_cell) % N) + N) % N;

                        int cell_idx =
                            neighbor_iz * N * N + neighbor_iy * N + neighbor_ix;

                        double cell_px = (neighbor_ix + 0.5) * cell_size;
                        double cell_py = (neighbor_iy + 0.5) * cell_size;
                        double cell_pz = (neighbor_iz + 0.5) * cell_size;
                        double gas_mass = d_gas_rho[cell_idx] * cell_vol;

                        double dx =
                            periodic_displacement(cell_px - p1_x, domain_size);
                        double dy =
                            periodic_displacement(cell_py - p1_y, domain_size);
                        double dz =
                            periodic_displacement(cell_pz - p1_z, domain_size);

                        double dist_sq = dx * dx + dy * dy + dz * dz;
                        if (use_pm && dist_sq > cutoff_sq) continue;
                        if (dist_sq == 0.0)
                            continue;  // Prevent NaN if particle is dead center

                        double r = std::sqrt(dist_sq);
                        double q = r / cell_size;

                        // Monaghan 1992 Cubic Spline Weighting (Models internal
                        // fluid Shell Theorem)
                        double spline_weight = 1.0;
                        if (q < 1.0) {
                            double q3 = q * q * q;
                            if (q < 0.5) {
                                double q5 = q3 * q * q;
                                double q6 = q5 * q;
                                spline_weight = (32.0 / 3.0) * q3 -
                                                (192.0 / 5.0) * q5 + 32.0 * q6;
                            } else {
                                double q4 = q3 * q;
                                double q5 = q4 * q;
                                double q6 = q5 * q;
                                spline_weight = -1.0 / 15.0 +
                                                (64.0 / 3.0) * q3 - 48.0 * q4 +
                                                (192.0 / 5.0) * q5 -
                                                (32.0 / 3.0) * q6;
                            }
                        }

                        // True Newtonian magnitude, scaled smoothly down to 0
                        // inside the cell
                        double a_pp_mag =
                            (G * gas_mass / dist_sq) * spline_weight;
                        double a_gas_mag =
                            (G * dm_mass / dist_sq) * spline_weight;

                        if (use_pm) {
                            double r_scaled = r / (2.0 * r_s);
                            double pm_filter =
                                (std::erfc(r_scaled) +
                                 (r / (std::sqrt(M_PI) * r_s)) *
                                     std::exp(-r_scaled * r_scaled));
                            a_pp_mag *= pm_filter;
                            a_gas_mag *= pm_filter;
                        }

                        // Accumulate locally for DM particle
                        local_acc_x += a_pp_mag * dx / r;
                        local_acc_y += a_pp_mag * dy / r;
                        local_acc_z += a_pp_mag * dz / r;

                        // Write equal and opposite reaction to Gas grid
                        double gas_ax = -a_gas_mag * dx / r;
                        double gas_ay = -a_gas_mag * dy / r;
                        double gas_az = -a_gas_mag * dz / r;

#pragma omp atomic
                        d_grav_x[cell_idx] += gas_ax;
#pragma omp atomic
                        d_grav_y[cell_idx] += gas_ay;
#pragma omp atomic
                        d_grav_z[cell_idx] += gas_az;
                    }
                }
            }

            d_ax[i] += local_acc_x;
            d_ay[i] += local_acc_y;
            d_az[i] += local_acc_z;
        }
#pragma omp target exit data map(                                         \
        from : d_ax[0 : n_parts], d_ay[0 : n_parts], d_az[0 : n_parts],   \
            d_grav_x[0 : N * N * N], d_grav_y[0 : N * N * N],             \
            d_grav_z[0 : N * N * N])                                      \
    map(delete : d_px[0 : n_parts], d_py[0 : n_parts], d_pz[0 : n_parts], \
            d_m[0 : n_parts])
    } else
#endif
    {
        // ========================================================================
        // CPU IMPLEMENTATION
        // ========================================================================
#pragma omp parallel for schedule(dynamic, 64)
        for (size_t i = 0; i < n_parts; ++i) {
            double p1_x = pos_x[i], p1_y = pos_y[i], p1_z = pos_z[i];

            int ix = static_cast<int>(p1_x / cell_size);
            int iy = static_cast<int>(p1_y / cell_size);
            int iz = static_cast<int>(p1_z / cell_size);

            double local_acc_x = 0.0, local_acc_y = 0.0, local_acc_z = 0.0;
            double dm_mass = mass[i];

            for (int dx_cell = -search_radius; dx_cell <= search_radius;
                 ++dx_cell) {
                for (int dy_cell = -search_radius; dy_cell <= search_radius;
                     ++dy_cell) {
                    for (int dz_cell = -search_radius; dz_cell <= search_radius;
                         ++dz_cell) {
                        int neighbor_ix = (((ix + dx_cell) % N) + N) % N;
                        int neighbor_iy = (((iy + dy_cell) % N) + N) % N;
                        int neighbor_iz = (((iz + dz_cell) % N) + N) % N;

                        int cell_idx =
                            neighbor_iz * N * N + neighbor_iy * N + neighbor_ix;

                        double cell_px = (neighbor_ix + 0.5) * cell_size;
                        double cell_py = (neighbor_iy + 0.5) * cell_size;
                        double cell_pz = (neighbor_iz + 0.5) * cell_size;
                        double gas_mass =
                            gas.get_density().data[cell_idx] * cell_vol;

                        double dx =
                            periodic_displacement(cell_px - p1_x, domain_size);
                        double dy =
                            periodic_displacement(cell_py - p1_y, domain_size);
                        double dz =
                            periodic_displacement(cell_pz - p1_z, domain_size);

                        double dist_sq = dx * dx + dy * dy + dz * dz;
                        if (use_pm && dist_sq > cutoff_sq) continue;
                        if (dist_sq == 0.0)
                            continue;  // Prevent NaN if particle is dead center

                        double r = std::sqrt(dist_sq);
                        double q = r / cell_size;

                        // Monaghan 1992 Cubic Spline Weighting
                        double spline_weight = 1.0;
                        if (q < 1.0) {
                            double q3 = q * q * q;
                            if (q < 0.5) {
                                double q5 = q3 * q * q;
                                double q6 = q5 * q;
                                spline_weight = (32.0 / 3.0) * q3 -
                                                (192.0 / 5.0) * q5 + 32.0 * q6;
                            } else {
                                double q4 = q3 * q;
                                double q5 = q4 * q;
                                double q6 = q5 * q;
                                spline_weight = -1.0 / 15.0 +
                                                (64.0 / 3.0) * q3 - 48.0 * q4 +
                                                (192.0 / 5.0) * q5 -
                                                (32.0 / 3.0) * q6;
                            }
                        }

                        // True Newtonian magnitude, scaled smoothly down to 0
                        // inside the cell
                        double a_pp_mag =
                            (G * gas_mass / dist_sq) * spline_weight;
                        double a_gas_mag =
                            (G * dm_mass / dist_sq) * spline_weight;

                        if (use_pm) {
                            double r_scaled = r / (2.0 * r_s);
                            double pm_filter =
                                (std::erfc(r_scaled) +
                                 (r / (std::sqrt(M_PI) * r_s)) *
                                     std::exp(-r_scaled * r_scaled));
                            a_pp_mag *= pm_filter;
                            a_gas_mag *= pm_filter;
                        }

                        // Accumulate locally for DM particle
                        local_acc_x += a_pp_mag * dx / r;
                        local_acc_y += a_pp_mag * dy / r;
                        local_acc_z += a_pp_mag * dz / r;

                        // Write equal and opposite reaction to Gas grid
                        double gas_ax = -a_gas_mag * dx / r;
                        double gas_ay = -a_gas_mag * dy / r;
                        double gas_az = -a_gas_mag * dz / r;

#pragma omp atomic
                        grav_x.data[cell_idx] += gas_ax;
#pragma omp atomic
                        grav_y.data[cell_idx] += gas_ay;
#pragma omp atomic
                        grav_z.data[cell_idx] += gas_az;
                    }
                }
            }

            acc_x[i] += local_acc_x;
            acc_y[i] += local_acc_y;
            acc_z[i] += local_acc_z;
        }
    }
}

double ParticleSystem::get_gravity_timestep(const Config& config) const {
    if (num_particles == 0) return std::numeric_limits<double>::infinity();

    double epsilon = std::sqrt(config.softening_squared);
    double a_max = std::sqrt(max_accel_sq);
    double dt_grav = std::sqrt(epsilon / a_max);

    return dt_grav * config.gravity_accuracy_eta;
}

void compute_and_add_generic_pp_forces(
    size_t n_parts, const double* __restrict__ pos_x,
    const double* __restrict__ pos_y, const double* __restrict__ pos_z,
    const double* __restrict__ mass, double* __restrict__ acc_x,
    double* __restrict__ acc_y, double* __restrict__ acc_z,
    const int* __restrict__ cell_start, const int* __restrict__ cell_count,
    const Config& config, Diagnostics& diag) {
    const int N = config.mesh_size;
    const double cell_size = config.cell_size;
    const double domain_size = config.domain_size;
    const double G = config.G;
    const double soft_sq = config.softening_squared;
    const double cutoff_sq = config.cutoff_radius_squared;
    const double r_s = config.PM_smoothing_cells * cell_size;
    const bool use_pm = config.use_PM;

    int dx_start =
        use_pm ? -static_cast<int>(ceil(config.cutoff_radius / cell_size)) : 0;
    int dx_end = use_pm
                     ? static_cast<int>(ceil(config.cutoff_radius / cell_size))
                     : N - 1;

#ifdef USE_GPU
    if (config.enable_GPU) {
        // ========================================================================
        // GPU IMPLEMENTATION (Sorted Array)
        // ========================================================================
        const int num_cells = N * N * N;

        // Aliasing pointers
        const double* d_px = pos_x;
        const double* d_py = pos_y;
        const double* d_pz = pos_z;
        const double* d_m = mass;
        double* d_ax = acc_x;
        double* d_ay = acc_y;
        double* d_az = acc_z;
        const int* d_cell_start = cell_start;
        const int* d_cell_count = cell_count;

        auto start_transfer = std::chrono::high_resolution_clock::now();
#pragma omp target enter data map(                                             \
        to : d_px[0 : n_parts], d_py[0 : n_parts], d_pz[0 : n_parts],          \
            d_m[0 : n_parts], d_cell_start[0 : num_cells],                     \
            d_cell_count[0 : num_cells], d_ax[0 : n_parts], d_ay[0 : n_parts], \
            d_az[0 : n_parts])
        auto end_transfer = std::chrono::high_resolution_clock::now();

        auto start_compute = std::chrono::high_resolution_clock::now();
#pragma omp target teams distribute parallel for
        for (size_t i = 0; i < n_parts; ++i) {
            double p1_x = d_px[i], p1_y = d_py[i], p1_z = d_pz[i];
            int ix = static_cast<int>(p1_x / cell_size);
            int iy = static_cast<int>(p1_y / cell_size);
            int iz = static_cast<int>(p1_z / cell_size);

            double local_acc_x = 0.0, local_acc_y = 0.0, local_acc_z = 0.0;

            for (int dx_cell = dx_start; dx_cell <= dx_end; ++dx_cell) {
                for (int dy_cell = dx_start; dy_cell <= dx_end; ++dy_cell) {
                    for (int dz_cell = dx_start; dz_cell <= dx_end; ++dz_cell) {
                        int neighbor_ix =
                            use_pm ? (((ix + dx_cell) % N) + N) % N : dx_cell;
                        int neighbor_iy =
                            use_pm ? (((iy + dy_cell) % N) + N) % N : dy_cell;
                        int neighbor_iz =
                            use_pm ? (((iz + dz_cell) % N) + N) % N : dz_cell;
                        int cell_idx =
                            neighbor_iz * N * N + neighbor_iy * N + neighbor_ix;

                        int start_idx = d_cell_start[cell_idx];
                        int end_idx = start_idx + d_cell_count[cell_idx];

                        for (int j = start_idx; j < end_idx; ++j) {
                            double dx = periodic_displacement(d_px[j] - p1_x,
                                                              domain_size);
                            double dy = periodic_displacement(d_py[j] - p1_y,
                                                              domain_size);
                            double dz = periodic_displacement(d_pz[j] - p1_z,
                                                              domain_size);

                            double dist_sq = dx * dx + dy * dy + dz * dz;

                            if (use_pm && dist_sq > cutoff_sq) continue;

                            double pp_dist_sq = dist_sq + soft_sq;
                            double pp_dist = std::sqrt(pp_dist_sq);
                            double a_pp = G * d_m[j] / pp_dist_sq;

                            if (use_pm) {
                                double r_scaled = pp_dist / (2.0 * r_s);
                                a_pp *= (std::erfc(r_scaled) +
                                         (pp_dist / (std::sqrt(M_PI) * r_s)) *
                                             std::exp(-r_scaled * r_scaled));
                            }

                            // If i == j, dx, dy, and dz are 0, adding 0.0 to
                            // acceleration
                            local_acc_x += a_pp * dx / pp_dist;
                            local_acc_y += a_pp * dy / pp_dist;
                            local_acc_z += a_pp * dz / pp_dist;
                        }
                    }
                }
            }
            d_ax[i] += local_acc_x;
            d_ay[i] += local_acc_y;
            d_az[i] += local_acc_z;
        }
        auto end_compute = std::chrono::high_resolution_clock::now();

        auto start_return = std::chrono::high_resolution_clock::now();
#pragma omp target exit data map(from : d_ax[0 : n_parts], d_ay[0 : n_parts], \
                                     d_az[0 : n_parts])                       \
    map(delete : d_px[0 : n_parts], d_py[0 : n_parts], d_pz[0 : n_parts],     \
            d_m[0 : n_parts], d_cell_start[0 : num_cells],                    \
            d_cell_count[0 : num_cells])
        auto end_return = std::chrono::high_resolution_clock::now();

        double diff_transf =
            std::chrono::duration_cast<std::chrono::microseconds>(
                end_transfer - start_transfer)
                .count();
        double diff_comput =
            std::chrono::duration_cast<std::chrono::microseconds>(end_compute -
                                                                  start_compute)
                .count();
        double diff_return =
            std::chrono::duration_cast<std::chrono::microseconds>(end_return -
                                                                  start_return)
                .count();

        diag.add_prof_time(ProfRegion::Transf, diff_transf);
        diag.add_prof_time(ProfRegion::Compute, diff_comput);
        diag.add_prof_time(ProfRegion::Ret, diff_return);
    } else
#endif
    {
        // ========================================================================
        // CPU IMPLEMENTATION (Sorted Array + Gather Buffer)
        // ========================================================================
// ========================================================================
// CPU IMPLEMENTATION (Sorted Array + Chunked Gather Buffer)
// ========================================================================
#pragma omp parallel for schedule(dynamic, 64)
        for (size_t i = 0; i < n_parts; ++i) {
            double p1_x = pos_x[i], p1_y = pos_y[i], p1_z = pos_z[i];
            int ix = static_cast<int>(p1_x / config.cell_size);
            int iy = static_cast<int>(p1_y / config.cell_size);
            int iz = static_cast<int>(p1_z / config.cell_size);

            constexpr int MAX_NEIGHBORS = 1536;
            double n_x[MAX_NEIGHBORS];
            double n_y[MAX_NEIGHBORS];
            double n_z[MAX_NEIGHBORS];
            double n_m[MAX_NEIGHBORS];
            int num_neighbors = 0;

            double local_acc_x = 0.0, local_acc_y = 0.0, local_acc_z = 0.0;

            // Lambda to execute the SIMD math and flush the buffer
            auto compute_simd_batch = [&]() {
#pragma omp simd reduction(+ : local_acc_x, local_acc_y, local_acc_z)
                for (int k = 0; k < num_neighbors; ++k) {
                    double dx =
                        periodic_displacement(n_x[k] - p1_x, domain_size);
                    double dy =
                        periodic_displacement(n_y[k] - p1_y, domain_size);
                    double dz =
                        periodic_displacement(n_z[k] - p1_z, domain_size);

                    double dist_sq = dx * dx + dy * dy + dz * dz;
                    if (use_pm && dist_sq > cutoff_sq) continue;

                    double pp_dist_sq = dist_sq + soft_sq;
                    double pp_dist = std::sqrt(pp_dist_sq);
                    double a_pp = G * n_m[k] / pp_dist_sq;

                    if (use_pm) {
                        double r_scaled = pp_dist / (2.0 * r_s);
                        a_pp *= (std::erfc(r_scaled) +
                                 (pp_dist / (std::sqrt(M_PI) * r_s)) *
                                     std::exp(-r_scaled * r_scaled));
                    }

                    local_acc_x += a_pp * dx / pp_dist;
                    local_acc_y += a_pp * dy / pp_dist;
                    local_acc_z += a_pp * dz / pp_dist;
                }
                num_neighbors = 0;  // Reset the buffer counter
            };

            // Gather: Linearly pull from sorted memory
            for (int dx_cell = dx_start; dx_cell <= dx_end; ++dx_cell) {
                for (int dy_cell = dx_start; dy_cell <= dx_end; ++dy_cell) {
                    for (int dz_cell = dx_start; dz_cell <= dx_end; ++dz_cell) {
                        int neighbor_ix =
                            use_pm ? (((ix + dx_cell) % N) + N) % N : dx_cell;
                        int neighbor_iy =
                            use_pm ? (((iy + dy_cell) % N) + N) % N : dy_cell;
                        int neighbor_iz =
                            use_pm ? (((iz + dz_cell) % N) + N) % N : dz_cell;
                        int cell_idx =
                            neighbor_iz * N * N + neighbor_iy * N + neighbor_ix;

                        int start_idx = cell_start[cell_idx];
                        int end_idx = start_idx + cell_count[cell_idx];

                        for (int j = start_idx; j < end_idx; ++j) {
                            if (i != j) {
                                n_x[num_neighbors] = pos_x[j];
                                n_y[num_neighbors] = pos_y[j];
                                n_z[num_neighbors] = pos_z[j];
                                n_m[num_neighbors] = mass[j];
                                num_neighbors++;

                                // Flush the buffer if it gets full
                                if (num_neighbors == MAX_NEIGHBORS) {
                                    compute_simd_batch();
                                }
                            }
                        }
                    }
                }
            }

            // Process any remaining particles in the buffer
            if (num_neighbors > 0) {
                compute_simd_batch();
            }

            acc_x[i] += local_acc_x;
            acc_y[i] += local_acc_y;
            acc_z[i] += local_acc_z;
        }
    }
}