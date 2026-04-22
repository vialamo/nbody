#include "particles.h"

#include <omp.h>

#include <cmath>
#include <limits>

#include "math_utils.h"

ParticleSystem::ParticleSystem(const Config& config)
    : dm_rho(config.MESH_SIZE),
      cic_data(config.NUM_DM_PARTICLES),
      cell_grid(config.MESH_SIZE * config.MESH_SIZE * config.MESH_SIZE) {
    particles.reserve(config.NUM_DM_PARTICLES);
}

void ParticleSystem::bin_and_assign_mass(const Config& config) {
    dm_rho.setZero();
    cic_data.assign(particles.size(), {});

    for (auto& cell : cell_grid) {
        cell.clear();
    }

    int N = config.MESH_SIZE;

    for (size_t i = 0; i < particles.size(); ++i) {
        const auto& p = particles[i];
        int ix = static_cast<int>(p.pos.x / config.CELL_SIZE);
        int iy = static_cast<int>(p.pos.y / config.CELL_SIZE);
        int iz = static_cast<int>(p.pos.z / config.CELL_SIZE);

        double frac_x = (p.pos.x / config.CELL_SIZE) - ix;
        double frac_y = (p.pos.y / config.CELL_SIZE) - iy;
        double frac_z = (p.pos.z / config.CELL_SIZE) - iz;

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

        dm_rho(ix0, iy0, iz0) += p.mass * w000;
        dm_rho(ix1, iy0, iz0) += p.mass * w100;
        dm_rho(ix0, iy1, iz0) += p.mass * w010;
        dm_rho(ix1, iy1, iz0) += p.mass * w110;
        dm_rho(ix0, iy0, iz1) += p.mass * w001;
        dm_rho(ix1, iy0, iz1) += p.mass * w101;
        dm_rho(ix0, iy1, iz1) += p.mass * w011;
        dm_rho(ix1, iy1, iz1) += p.mass * w111;

        int cell_index = iz0 * N * N + iy0 * N + ix0;
        cell_grid[cell_index].push_back(static_cast<int>(i));
    }

    dm_rho.data /= config.CELL_VOLUME;
}

void ParticleSystem::interpolate_cic_forces(const Grid3D& ax_grid,
                                            const Grid3D& ay_grid,
                                            const Grid3D& az_grid,
                                            std::vector<Vec3>& forces,
                                            const Config& config) {
    forces.assign(particles.size(), {0.0, 0.0, 0.0});
    int N = config.MESH_SIZE;

#pragma omp parallel for schedule(dynamic, 64)
    for (size_t i = 0; i < particles.size(); ++i) {
        const auto& p = particles[i];
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

        forces[i].x = interp(ax_grid) * p.mass;
        forces[i].y = interp(ay_grid) * p.mass;
        forces[i].z = interp(az_grid) * p.mass;
    }
}

void ParticleSystem::compute_pp_forces(std::vector<Vec3>& pp_forces,
                                       const Config& config) {
    pp_forces.assign(particles.size(), {0.0, 0.0, 0.0});
    int search_radius;
    if (config.USE_PM) {
        search_radius =
            static_cast<int>(ceil(config.CUTOFF_RADIUS / config.CELL_SIZE));
    } else {
        search_radius = config.MESH_SIZE / 2;
    }
    int N = config.MESH_SIZE;

// Each thread gets its own 'i'.
#pragma omp parallel for schedule(dynamic, 64)
    for (size_t i = 0; i < particles.size(); ++i) {
        auto& p1 = particles[i];
        int ix = static_cast<int>(p1.pos.x / config.CELL_SIZE);
        int iy = static_cast<int>(p1.pos.y / config.CELL_SIZE);
        int iz = static_cast<int>(p1.pos.z / config.CELL_SIZE);

        // Local accumulator to avoid hammering the main array in memory
        Vec3 local_force = {0.0, 0.0, 0.0};

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

                    for (int j : cell_grid[neighbor_cell_index]) {
                        // Because OMP, we must compute force from j on i, even
                        // if we already did i on j
                        if (i == j) continue;

                        auto& p2 = particles[j];

                        double dx = periodic_displacement(p2.pos.x - p1.pos.x,
                                                          config.DOMAIN_SIZE);
                        double dy = periodic_displacement(p2.pos.y - p1.pos.y,
                                                          config.DOMAIN_SIZE);
                        double dz = periodic_displacement(p2.pos.z - p1.pos.z,
                                                          config.DOMAIN_SIZE);
                        double dist_sq = dx * dx + dy * dy + dz * dz;

                        if (config.USE_PM &&
                            dist_sq > config.CUTOFF_RADIUS_SQUARED)
                            continue;

                        double S = 1.0;
                        if (config.USE_PM &&
                            dist_sq > config.R_SWITCH_START_SQ) {
                            double dist = sqrt(dist_sq);
                            double x = (dist - config.R_SWITCH_START) /
                                       config.CUTOFF_TRANSITION_WIDTH;
                            S = 2 * pow(x, 3) - 3 * pow(x, 2) + 1;
                        }

                        double soft_dist_sq =
                            dist_sq + pow(0.5 * config.CELL_SIZE, 2);
                        double f_pm_short =
                            config.G * p1.mass * p2.mass / soft_dist_sq;
                        double soft_dist = sqrt(soft_dist_sq);
                        Vec3 f_pm_short_vec = {f_pm_short * dx / soft_dist,
                                               f_pm_short * dy / soft_dist,
                                               f_pm_short * dz / soft_dist};

                        double pp_dist_sq = dist_sq + config.SOFTENING_SQUARED;
                        double f_pp = config.G * p1.mass * p2.mass / pp_dist_sq;
                        double pp_dist = sqrt(pp_dist_sq);
                        Vec3 f_pp_vec = {f_pp * dx / pp_dist,
                                         f_pp * dy / pp_dist,
                                         f_pp * dz / pp_dist};

                        if (!config.USE_PM) f_pm_short_vec = {0.0, 0.0, 0.0};

                        // Accumulate locally
                        local_force.x += S * (f_pp_vec.x - f_pm_short_vec.x);
                        local_force.y += S * (f_pp_vec.y - f_pm_short_vec.y);
                        local_force.z += S * (f_pp_vec.z - f_pm_short_vec.z);
                    }
                }
            }
        }
        // Write the finalized force to memory ONCE per particle
        pp_forces[i].x += local_force.x;
        pp_forces[i].y += local_force.y;
        pp_forces[i].z += local_force.z;
    }
}

double ParticleSystem::calculate_kinetic_energy(double a) const {
    double kinetic_energy = 0.0;
    
#pragma omp parallel for reduction(+ : kinetic_energy)
    for (size_t i = 0; i < particles.size(); ++i) {
        const auto& p = particles[i];
        
        double proper_vel_sq = (a * p.vel.x) * (a * p.vel.x) +
                               (a * p.vel.y) * (a * p.vel.y) +
                               (a * p.vel.z) * (a * p.vel.z);
                               
        kinetic_energy += 0.5 * p.mass * proper_vel_sq;
    }
    return kinetic_energy;
}

double ParticleSystem::get_gravity_timestep(const Config& config) const {
    if (particles.empty()) return std::numeric_limits<double>::infinity();

    double max_accel_sq = 1e-9;
    for (const auto& p : particles) {
        double accel_sq =
            p.acc.x * p.acc.x + p.acc.y * p.acc.y + p.acc.z * p.acc.z;
        if (accel_sq > max_accel_sq) max_accel_sq = accel_sq;
    }
    double dt_grav = sqrt(config.SOFTENING_SQUARED) / sqrt(max_accel_sq);
    return dt_grav * config.GRAVITY_DT_FACTOR;
}