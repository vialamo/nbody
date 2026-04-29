#include "particles.h"

#include <omp.h>

#include <cmath>
#include <limits>

#include "math_utils.h"

ParticleSystem::ParticleSystem(const Config& config)
    : dm_rho(config.MESH_SIZE),
      cic_data(config.NUM_DM_PARTICLES),
      max_accel_sq(0.0) {
    pos_x.reserve(config.NUM_DM_PARTICLES);
    pos_y.reserve(config.NUM_DM_PARTICLES);
    pos_z.reserve(config.NUM_DM_PARTICLES);

    vel_x.reserve(config.NUM_DM_PARTICLES);
    vel_y.reserve(config.NUM_DM_PARTICLES);
    vel_z.reserve(config.NUM_DM_PARTICLES);

    acc_x.reserve(config.NUM_DM_PARTICLES);
    acc_y.reserve(config.NUM_DM_PARTICLES);
    acc_z.reserve(config.NUM_DM_PARTICLES);

    mass.reserve(config.NUM_DM_PARTICLES);

    int num_cells = config.MESH_SIZE * config.MESH_SIZE * config.MESH_SIZE;
    cell_list.resize(num_cells, config.NUM_DM_PARTICLES);
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
    cell_list.reset();

    int N = config.MESH_SIZE;

    for (size_t i = 0; i < num_particles; ++i) {
        double px = pos_x[i];
        double py = pos_y[i];
        double pz = pos_z[i];
        double m = mass[i];

        int ix = static_cast<int>(px / config.CELL_SIZE);
        int iy = static_cast<int>(py / config.CELL_SIZE);
        int iz = static_cast<int>(pz / config.CELL_SIZE);

        double frac_x = (px / config.CELL_SIZE) - ix;
        double frac_y = (py / config.CELL_SIZE) - iy;
        double frac_z = (pz / config.CELL_SIZE) - iz;

        // 8 corners of the 3D cell
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

        dm_rho(ix0, iy0, iz0) += m * w000;
        dm_rho(ix1, iy0, iz0) += m * w100;
        dm_rho(ix0, iy1, iz0) += m * w010;
        dm_rho(ix1, iy1, iz0) += m * w110;
        dm_rho(ix0, iy0, iz1) += m * w001;
        dm_rho(ix1, iy0, iz1) += m * w101;
        dm_rho(ix0, iy1, iz1) += m * w011;
        dm_rho(ix1, iy1, iz1) += m * w111;

        int cell_index = iz0 * N * N + iy0 * N + ix0;
        cell_list.next_particle[i] = cell_list.first_in_cell[cell_index];
        cell_list.first_in_cell[cell_index] = static_cast<int>(i);
    }

    dm_rho.data /= config.CELL_VOLUME;
}

void ParticleSystem::interpolate_cic_forces(const Grid3D& ax_grid,
                                            const Grid3D& ay_grid,
                                            const Grid3D& az_grid,
                                            const Config& config) {
    int N = config.MESH_SIZE;

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

void ParticleSystem::compute_pp_forces(const Config& config) {
    int search_radius =
        config.USE_PM
            ? static_cast<int>(ceil(config.CUTOFF_RADIUS / config.CELL_SIZE))
            : config.MESH_SIZE / 2;

    int N = config.MESH_SIZE;
    int num_cells = N * N * N;

    // Pre-calculate constants
    const double domain_size = config.DOMAIN_SIZE;
    const double G = config.G;
    const double soft_sq = config.SOFTENING_SQUARED;
    const double pm_soft_sq =
        (0.5 * config.CELL_SIZE) * (0.5 * config.CELL_SIZE);
    const double r_switch_start = config.R_SWITCH_START;
    const double r_switch_sq = config.R_SWITCH_START_SQ;
    const double cutoff_sq = config.CUTOFF_RADIUS_SQUARED;
    const double trans_width = config.CUTOFF_TRANSITION_WIDTH;
    const bool use_pm = config.USE_PM;
    const size_t n_parts = num_particles;

#ifdef USE_GPU
    // ========================================================================
    // GPU IMPLEMENTATION (SIMT / Target Offload)
    // ========================================================================
    double* d_px = pos_x.data();
    double* d_py = pos_y.data();
    double* d_pz = pos_z.data();
    double* d_m = mass.data();
    double* d_ax = acc_x.data();
    double* d_ay = acc_y.data();
    double* d_az = acc_z.data();
    int* d_head = cell_list.first_in_cell.data();
    int* d_next = cell_list.next_particle.data();

#pragma omp target teams distribute parallel for map(                 \
        to : d_px[0 : n_parts], d_py[0 : n_parts], d_pz[0 : n_parts], \
            d_m[0 : n_parts])                                         \
    map(to : d_head[0 : num_cells], d_next[0 : n_parts])              \
    map(tofrom : d_ax[0 : n_parts], d_ay[0 : n_parts], d_az[0 : n_parts])
    for (size_t i = 0; i < n_parts; ++i) {
        double p1_x = d_px[i];
        double p1_y = d_py[i];
        double p1_z = d_pz[i];

        int ix = static_cast<int>(p1_x / config.CELL_SIZE);
        int iy = static_cast<int>(p1_y / config.CELL_SIZE);
        int iz = static_cast<int>(p1_z / config.CELL_SIZE);

        double local_acc_x = 0.0;
        double local_acc_y = 0.0;
        double local_acc_z = 0.0;

        for (int dx_cell = -search_radius; dx_cell <= search_radius;
             ++dx_cell) {
            for (int dy_cell = -search_radius; dy_cell <= search_radius;
                 ++dy_cell) {
                for (int dz_cell = -search_radius; dz_cell <= search_radius;
                     ++dz_cell) {
                    int neighbor_ix = (ix + dx_cell + N) % N;
                    int neighbor_iy = (iy + dy_cell + N) % N;
                    int neighbor_iz = (iz + dz_cell + N) % N;
                    int neighbor_cell_index =
                        neighbor_iz * N * N + neighbor_iy * N + neighbor_ix;

                    int j = d_head[neighbor_cell_index];

                    while (j != -1) {
                        if (i == j) {
                            j = d_next[j];
                            continue;
                        }

                        double p2_x = d_px[j];
                        double p2_y = d_py[j];
                        double p2_z = d_pz[j];
                        double p2_m = d_m[j];

                        double dx =
                            periodic_displacement(p2_x - p1_x, domain_size);
                        double dy =
                            periodic_displacement(p2_y - p1_y, domain_size);
                        double dz =
                            periodic_displacement(p2_z - p1_z, domain_size);
                        double dist_sq = dx * dx + dy * dy + dz * dz;

                        if (use_pm && dist_sq > cutoff_sq) {
                            j = d_next[j];
                            continue;
                        }

                        double S = 1.0;
                        if (use_pm && dist_sq > r_switch_sq) {
                            double dist = sqrt(dist_sq);
                            double x = (dist - r_switch_start) / trans_width;
                            S = 2.0 * (x * x * x) - 3.0 * (x * x) + 1.0;
                        }

                        double a_pm_short_x = 0.0, a_pm_short_y = 0.0,
                               a_pm_short_z = 0.0;

                        if (use_pm) {
                            double soft_dist_sq = dist_sq + pm_soft_sq;
                            double soft_dist = sqrt(soft_dist_sq);
                            double a_pm_short = G * p2_m / soft_dist_sq;
                            a_pm_short_x = a_pm_short * dx / soft_dist;
                            a_pm_short_y = a_pm_short * dy / soft_dist;
                            a_pm_short_z = a_pm_short * dz / soft_dist;
                        }

                        double pp_dist_sq = dist_sq + soft_sq;
                        double pp_dist = sqrt(pp_dist_sq);
                        double a_pp = G * p2_m / pp_dist_sq;

                        double a_pp_x = a_pp * dx / pp_dist;
                        double a_pp_y = a_pp * dy / pp_dist;
                        double a_pp_z = a_pp * dz / pp_dist;

                        local_acc_x += S * (a_pp_x - a_pm_short_x);
                        local_acc_y += S * (a_pp_y - a_pm_short_y);
                        local_acc_z += S * (a_pp_z - a_pm_short_z);

                        j = d_next[j];
                    }
                }
            }
        }
        d_ax[i] += local_acc_x;
        d_ay[i] += local_acc_y;
        d_az[i] += local_acc_z;
    }

#else
// ========================================================================
// CPU IMPLEMENTATION (Gather Buffer)
// ========================================================================
#pragma omp parallel for schedule(dynamic, 64)
    for (size_t i = 0; i < n_parts; ++i) {
        double p1_x = pos_x[i];
        double p1_y = pos_y[i];
        double p1_z = pos_z[i];

        int ix = static_cast<int>(p1_x / config.CELL_SIZE);
        int iy = static_cast<int>(p1_y / config.CELL_SIZE);
        int iz = static_cast<int>(p1_z / config.CELL_SIZE);

        constexpr int MAX_NEIGHBORS = 1500;
        double n_x[MAX_NEIGHBORS];
        double n_y[MAX_NEIGHBORS];
        double n_z[MAX_NEIGHBORS];
        double n_m[MAX_NEIGHBORS];
        int num_neighbors = 0;

        for (int dx_cell = -search_radius; dx_cell <= search_radius;
             ++dx_cell) {
            for (int dy_cell = -search_radius; dy_cell <= search_radius;
                 ++dy_cell) {
                for (int dz_cell = -search_radius; dz_cell <= search_radius;
                     ++dz_cell) {
                    int neighbor_ix = (ix + dx_cell + N) % N;
                    int neighbor_iy = (iy + dy_cell + N) % N;
                    int neighbor_iz = (iz + dz_cell + N) % N;
                    int neighbor_cell_index =
                        neighbor_iz * N * N + neighbor_iy * N + neighbor_ix;

                    int j = cell_list.first_in_cell[neighbor_cell_index];

                    while (j != -1) {
                        if (i != j && num_neighbors < MAX_NEIGHBORS) {
                            n_x[num_neighbors] = pos_x[j];
                            n_y[num_neighbors] = pos_y[j];
                            n_z[num_neighbors] = pos_z[j];
                            n_m[num_neighbors] = mass[j];
                            num_neighbors++;
                        }
                        j = cell_list.next_particle[j];
                    }
                }
            }
        }

        double local_acc_x = 0.0;
        double local_acc_y = 0.0;
        double local_acc_z = 0.0;

#pragma omp simd reduction(+ : local_acc_x, local_acc_y, local_acc_z)
        for (int k = 0; k < num_neighbors; ++k) {
            double dx = periodic_displacement(n_x[k] - p1_x, domain_size);
            double dy = periodic_displacement(n_y[k] - p1_y, domain_size);
            double dz = periodic_displacement(n_z[k] - p1_z, domain_size);
            double dist_sq = dx * dx + dy * dy + dz * dz;

            if (use_pm && dist_sq > cutoff_sq) {
                continue;
            }

            double S = 1.0;
            if (use_pm && dist_sq > r_switch_sq) {
                double dist = sqrt(dist_sq);
                double x = (dist - r_switch_start) / trans_width;
                S = 2.0 * (x * x * x) - 3.0 * (x * x) + 1.0;
            }

            double a_pm_short_x = 0.0, a_pm_short_y = 0.0, a_pm_short_z = 0.0;

            if (use_pm) {
                double mask = use_pm ? 1 : 0;
                double soft_dist_sq = dist_sq + pm_soft_sq;
                double soft_dist = sqrt(soft_dist_sq);
                double a_pm_short = G * n_m[k] / soft_dist_sq;
                a_pm_short_x = a_pm_short * dx / soft_dist;
                a_pm_short_y = a_pm_short * dy / soft_dist;
                a_pm_short_z = a_pm_short * dz / soft_dist;
            }

            double pp_dist_sq = dist_sq + soft_sq;
            double pp_dist = sqrt(pp_dist_sq);
            double a_pp = G * n_m[k] / pp_dist_sq;

            double a_pp_x = a_pp * dx / pp_dist;
            double a_pp_y = a_pp * dy / pp_dist;
            double a_pp_z = a_pp * dz / pp_dist;

            local_acc_x += S * (a_pp_x - a_pm_short_x);
            local_acc_y += S * (a_pp_y - a_pm_short_y);
            local_acc_z += S * (a_pp_z - a_pm_short_z);
        }

        acc_x[i] += local_acc_x;
        acc_y[i] += local_acc_y;
        acc_z[i] += local_acc_z;
    }
#endif
}

double ParticleSystem::calculate_kinetic_energy(double a) const {
    double kinetic_energy = 0.0;

#pragma omp parallel for reduction(+ : kinetic_energy)
    for (size_t i = 0; i < num_particles; ++i) {
        double vx = vel_x[i];
        double vy = vel_y[i];
        double vz = vel_z[i];

        double proper_vel_sq =
            (a * vx) * (a * vx) + (a * vy) * (a * vy) + (a * vz) * (a * vz);

        kinetic_energy += 0.5 * mass[i] * proper_vel_sq;
    }
    return kinetic_energy;
}

double ParticleSystem::get_gravity_timestep(const Config& config) const {
    if (num_particles == 0) return std::numeric_limits<double>::infinity();

    double dt_grav = sqrt(config.SOFTENING_SQUARED) / sqrt(max_accel_sq);
    return dt_grav * config.GRAVITY_DT_FACTOR;
}