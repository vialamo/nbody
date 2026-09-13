#include "particles.h"

#include <omp.h>

#include <algorithm>  // For std::sort
#include <cmath>
#include <limits>
#include <numeric>  // For std::iota

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

    size_t num_nodes = 2 * config.num_dm_particles - 1;
    if (config.num_dm_particles > 0) {
        morton_codes.reserve(config.num_dm_particles);
        sorted_indices.reserve(config.num_dm_particles);
        bvh_nodes.reserve(num_nodes);
    }
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

    // Calculate cells & densities for PM grid
    for (size_t i = 0; i < num_particles; ++i) {
        double px = pos_x[i], py = pos_y[i], pz = pos_z[i];

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
    }

    dm_rho.data /= config.cell_volume;
}

void ParticleSystem::build_lbvh(const Config& config) {
    if (num_particles == 0) return;

    morton_codes.resize(num_particles);
    sorted_indices.resize(num_particles);
    bvh_nodes.resize(2 * num_particles - 1);

    double inv_domain = 1.0 / config.domain_size;
    // We use 21 bits per dimension (2^21 = 2097152) to fit in a 64-bit int
    double bound = 2097152.0;

    // Compute Morton codes for all particles
#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < num_particles; ++i) {
        // Normalize coordinates to [0, 1) and scale to the 21-bit integer range
        uint32_t x = static_cast<uint32_t>(
            fmod(pos_x[i] * inv_domain + 1.0, 1.0) * bound);
        uint32_t y = static_cast<uint32_t>(
            fmod(pos_y[i] * inv_domain + 1.0, 1.0) * bound);
        uint32_t z = static_cast<uint32_t>(
            fmod(pos_z[i] * inv_domain + 1.0, 1.0) * bound);

        morton_codes[i] = morton3D(x, y, z);
    }

    // Initialize indices and sort them based on the Morton codes
    std::iota(sorted_indices.begin(), sorted_indices.end(), 0);
    std::sort(sorted_indices.begin(), sorted_indices.end(),
              [&](int a, int b) { return morton_codes[a] < morton_codes[b]; });

    // Rearrange particle arrays to match the new sorted order
    std::vector<double> new_px(num_particles), new_py(num_particles),
        new_pz(num_particles);
    std::vector<double> new_vx(num_particles), new_vy(num_particles),
        new_vz(num_particles);
    std::vector<double> new_ax(num_particles), new_ay(num_particles),
        new_az(num_particles);
    std::vector<double> new_m(num_particles);
    // CIC_Data is only needed if for PM gravity
    // We sort it here just in case
    std::vector<CIC_Data> new_cic(num_particles);
    std::vector<uint64_t> new_morton(num_particles);

    for (size_t i = 0; i < num_particles; ++i) {
        int src = sorted_indices[i];
        new_px[i] = pos_x[src];
        new_py[i] = pos_y[src];
        new_pz[i] = pos_z[src];
        new_vx[i] = vel_x[src];
        new_vy[i] = vel_y[src];
        new_vz[i] = vel_z[src];
        new_ax[i] = acc_x[src];
        new_ay[i] = acc_y[src];
        new_az[i] = acc_z[src];
        new_m[i] = mass[src];

        if (!cic_data.empty()) new_cic[i] = cic_data[src];
        new_morton[i] = morton_codes[src];
    }

    // Move sorted data back
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
    if (!cic_data.empty()) cic_data = std::move(new_cic);
    morton_codes = std::move(new_morton);

    // TREE TOPOLOGY CONSTRUCTION (Karras 2012)

    // Total nodes = 2N - 1.
    // Indices [0, N-2] are internal nodes.
    // Indices [N-1, 2N-2] are leaf nodes.

    // Utility lambda to find the longest common prefix between two Morton codes
    auto delta = [&](int i, int j) -> int {
        if (j < 0 || j >= num_particles) return -1;
        uint64_t code_i = morton_codes[i];
        uint64_t code_j = morton_codes[j];

        if (code_i == code_j) {
            // Tie-breaker for identical coordinates using the original index
            return 64 + __builtin_clzll(static_cast<unsigned long long>(i ^ j));
        }
        return __builtin_clzll(code_i ^ code_j);
    };

// Initialize Leaf Nodes
#pragma omp parallel for schedule(static)
    for (int i = 0; i < num_particles; ++i) {
        int leaf_idx = num_particles - 1 + i;
        bvh_nodes[leaf_idx].particle_idx = i;
        bvh_nodes[leaf_idx].left_child = -1;
        bvh_nodes[leaf_idx].right_child = -1;
        // The parent will be set by the internal node that points to this leaf
    }

// Construct Internal Nodes in parallel
#pragma omp parallel for schedule(static)
    for (int i = 0; i < num_particles - 1; ++i) {
        // Determine direction of the range (+1 or -1)
        int d = (delta(i, i + 1) - delta(i, i - 1)) > 0 ? 1 : -1;

        // Compute upper bound for the length of the range
        int delta_min = delta(i, i - d);
        int l_max = 2;
        while (delta(i, i + l_max * d) > delta_min) {
            l_max *= 2;
        }

        // Find the other end of the range using binary search
        int l = 0;
        for (int t = l_max / 2; t >= 1; t /= 2) {
            if (delta(i, i + (l + t) * d) > delta_min) {
                l += t;
            }
        }
        int j = i + l * d;

        // Find the split position using binary search
        int delta_node = delta(i, j);
        int s = 0;
        int t = l;
        do {
            t = (t + 1) >> 1;  // ceil(t/2)
            if (s + t < l && delta(i, i + (s + t) * d) > delta_node) {
                s += t;
            }
        } while (t > 1);

        int split = i + s * d + std::min(d, 0);
        int min_idx = std::min(i, j);
        int max_idx = std::max(i, j);

        // Assign children
        int left_child, right_child;

        if (min_idx == split) {
            left_child = num_particles - 1 + split;  // Points to leaf
        } else {
            left_child = split;  // Points to internal node
        }

        if (max_idx == split + 1) {
            right_child = num_particles - 1 + split + 1;  // Points to leaf
        } else {
            right_child = split + 1;  // Points to internal node
        }

        bvh_nodes[i].left_child = left_child;
        bvh_nodes[i].right_child = right_child;
        bvh_nodes[i].particle_idx = -1;  // -1 indicates an internal node

        // Assign parent pointers to children
        bvh_nodes[left_child].parent = i;
        bvh_nodes[right_child].parent = i;
    }

    // Set the root node's parent to itself or -1
    bvh_nodes[0].parent = -1;

    // BOTTOM-UP AGGREGATION (Bounding Boxes & Center of Mass)

    // Counter for each internal node to track when both children are processed
    std::vector<int> atomic_flags(num_particles - 1, 0);

// Initialize Leaf Nodes and trigger the walk up
#pragma omp parallel for schedule(static)
    for (int i = 0; i < num_particles; ++i) {
        int leaf_idx = num_particles - 1 + i;
        double px = pos_x[i], py = pos_y[i], pz = pos_z[i];

        double radius = 0.0;

        bvh_nodes[leaf_idx].bbox.min_x = px - radius;
        bvh_nodes[leaf_idx].bbox.max_x = px + radius;
        bvh_nodes[leaf_idx].bbox.min_y = py - radius;
        bvh_nodes[leaf_idx].bbox.max_y = py + radius;
        bvh_nodes[leaf_idx].bbox.min_z = pz - radius;
        bvh_nodes[leaf_idx].bbox.max_z = pz + radius;

        bvh_nodes[leaf_idx].max_h = radius;

        bvh_nodes[leaf_idx].mass = mass[i];
        // Store mass-weighted positions temporarily to make summation easy
        bvh_nodes[leaf_idx].com_x = px * mass[i];
        bvh_nodes[leaf_idx].com_y = py * mass[i];
        bvh_nodes[leaf_idx].com_z = pz * mass[i];

        // Walk up the tree
        int curr = bvh_nodes[leaf_idx].parent;
        while (curr != -1) {
            int old_flag;
#pragma omp atomic capture
            {
                old_flag = atomic_flags[curr];
                atomic_flags[curr]++;
            }

            if (old_flag == 0) {
                // First thread to arrive. The other child isn't ready yet.
                // Terminate.
                break;
            }

            // Second thread to arrive. Both children are ready. Compute parent.
            int left = bvh_nodes[curr].left_child;
            int right = bvh_nodes[curr].right_child;

            // Combine Bounding Boxes
            bvh_nodes[curr].bbox.min_x = std::min(bvh_nodes[left].bbox.min_x,
                                                  bvh_nodes[right].bbox.min_x);
            bvh_nodes[curr].bbox.max_x = std::max(bvh_nodes[left].bbox.max_x,
                                                  bvh_nodes[right].bbox.max_x);
            bvh_nodes[curr].bbox.min_y = std::min(bvh_nodes[left].bbox.min_y,
                                                  bvh_nodes[right].bbox.min_y);
            bvh_nodes[curr].bbox.max_y = std::max(bvh_nodes[left].bbox.max_y,
                                                  bvh_nodes[right].bbox.max_y);
            bvh_nodes[curr].bbox.min_z = std::min(bvh_nodes[left].bbox.min_z,
                                                  bvh_nodes[right].bbox.min_z);
            bvh_nodes[curr].bbox.max_z = std::max(bvh_nodes[left].bbox.max_z,
                                                  bvh_nodes[right].bbox.max_z);

            bvh_nodes[curr].max_h =
                std::max(bvh_nodes[left].max_h, bvh_nodes[right].max_h);

            // Sum Mass and mass-weighted positions
            bvh_nodes[curr].mass = bvh_nodes[left].mass + bvh_nodes[right].mass;
            bvh_nodes[curr].com_x =
                bvh_nodes[left].com_x + bvh_nodes[right].com_x;
            bvh_nodes[curr].com_y =
                bvh_nodes[left].com_y + bvh_nodes[right].com_y;
            bvh_nodes[curr].com_z =
                bvh_nodes[left].com_z + bvh_nodes[right].com_z;

            // Move up to the next parent
            curr = bvh_nodes[curr].parent;
        }
    }

// Normalize Center of Mass
#pragma omp parallel for schedule(static)
    for (int i = 0; i < 2 * num_particles - 1; ++i) {
        if (bvh_nodes[i].mass > 0.0) {
            bvh_nodes[i].com_x /= bvh_nodes[i].mass;
            bvh_nodes[i].com_y /= bvh_nodes[i].mass;
            bvh_nodes[i].com_z /= bvh_nodes[i].mass;
        }
    }
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
    if (num_particles == 0) return;

    compute_and_add_generic_pp_forces(num_particles, pos_x.data(), pos_y.data(),
                                      pos_z.data(), mass.data(), acc_x.data(),
                                      acc_y.data(), acc_z.data(),
                                      bvh_nodes.data(), config, diag);
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
    const BVHNode* __restrict__ bvh_nodes, const Config& config,
    Diagnostics& diag) {
    if (n_parts == 0) return;

    const double domain_size = config.domain_size;
    const double G = config.G;
    const double soft_sq = config.softening_squared;
    const double cutoff_sq = config.cutoff_radius_squared;
    const double r_s = config.PM_smoothing_cells * config.cell_size;
    const bool use_pm = config.use_PM;
    const size_t num_nodes = 2 * n_parts - 1;

    const double search_sq =
        use_pm ? cutoff_sq : std::numeric_limits<double>::infinity();

#ifdef USE_GPU
    if (config.enable_GPU) {
        auto start_transfer = std::chrono::high_resolution_clock::now();

#pragma omp target enter data map(                                           \
        to : pos_x[0 : n_parts], pos_y[0 : n_parts], pos_z[0 : n_parts],     \
            mass[0 : n_parts], bvh_nodes[0 : num_nodes], acc_x[0 : n_parts], \
            acc_y[0 : n_parts], acc_z[0 : n_parts])

        auto end_transfer = std::chrono::high_resolution_clock::now();
        auto start_compute = std::chrono::high_resolution_clock::now();

#pragma omp target teams distribute parallel for
        for (size_t i = 0; i < n_parts; ++i) {
            double p1_x = pos_x[i], p1_y = pos_y[i], p1_z = pos_z[i];
            double local_acc_x = 0.0, local_acc_y = 0.0, local_acc_z = 0.0;

            int stack[128];
            int stack_ptr = 0;
            stack[stack_ptr++] = 0;

            while (stack_ptr > 0) {
                int node_idx = stack[--stack_ptr];
                const BVHNode& node = bvh_nodes[node_idx];

                double aabb_dist_sq =
                    min_periodic_dist_sq(p1_x, node.bbox.min_x, node.bbox.max_x,
                                         domain_size) +
                    min_periodic_dist_sq(p1_y, node.bbox.min_y, node.bbox.max_y,
                                         domain_size) +
                    min_periodic_dist_sq(p1_z, node.bbox.min_z, node.bbox.max_z,
                                         domain_size);

                if (aabb_dist_sq > search_sq) continue;

                if (node.particle_idx != -1) {
                    int j = node.particle_idx;
                    if (static_cast<size_t>(j) == i) continue;

                    double dx = p1_x - pos_x[j];
                    if (dx > 0.5 * domain_size)
                        dx -= domain_size;
                    else if (dx < -0.5 * domain_size)
                        dx += domain_size;

                    double dy = p1_y - pos_y[j];
                    if (dy > 0.5 * domain_size)
                        dy -= domain_size;
                    else if (dy < -0.5 * domain_size)
                        dy += domain_size;

                    double dz = p1_z - pos_z[j];
                    if (dz > 0.5 * domain_size)
                        dz -= domain_size;
                    else if (dz < -0.5 * domain_size)
                        dz += domain_size;

                    dx = -dx;
                    dy = -dy;
                    dz = -dz;

                    double dist_sq = dx * dx + dy * dy + dz * dz;

                    if (use_pm && dist_sq > cutoff_sq) continue;

                    double pp_dist_sq = dist_sq + soft_sq;
                    double pp_dist = std::sqrt(pp_dist_sq);
                    double a_pp = G * node.mass / pp_dist_sq;

                    if (use_pm) {
                        double r_scaled = pp_dist / (2.0 * r_s);
                        a_pp *= (std::erfc(r_scaled) +
                                 (pp_dist / (std::sqrt(M_PI) * r_s)) *
                                     std::exp(-r_scaled * r_scaled));
                    }

                    local_acc_x += a_pp * dx / pp_dist;
                    local_acc_y += a_pp * dy / pp_dist;
                    local_acc_z += a_pp * dz / pp_dist;
                } else {
                    stack[stack_ptr++] = node.left_child;
                    stack[stack_ptr++] = node.right_child;
                }
            }

            acc_x[i] += local_acc_x;
            acc_y[i] += local_acc_y;
            acc_z[i] += local_acc_z;
        }

        auto end_compute = std::chrono::high_resolution_clock::now();
        auto start_return = std::chrono::high_resolution_clock::now();

#pragma omp target exit data map(from : acc_x[0 : n_parts],                  \
                                     acc_y[0 : n_parts], acc_z[0 : n_parts]) \
    map(delete : pos_x[0 : n_parts], pos_y[0 : n_parts], pos_z[0 : n_parts], \
            mass[0 : n_parts], bvh_nodes[0 : num_nodes])

        auto end_return = std::chrono::high_resolution_clock::now();

        diag.add_prof_time(
            ProfRegion::Transf,
            std::chrono::duration_cast<std::chrono::microseconds>(
                end_transfer - start_transfer)
                .count());
        diag.add_prof_time(
            ProfRegion::Compute,
            std::chrono::duration_cast<std::chrono::microseconds>(end_compute -
                                                                  start_compute)
                .count());
        diag.add_prof_time(
            ProfRegion::Ret,
            std::chrono::duration_cast<std::chrono::microseconds>(end_return -
                                                                  start_return)
                .count());

    } else
#endif
    {
        // ========================================================================
        // CPU IMPLEMENTATION
        // ========================================================================
#pragma omp parallel for schedule(dynamic, 64)
        for (size_t i = 0; i < n_parts; ++i) {
            double p1_x = pos_x[i], p1_y = pos_y[i], p1_z = pos_z[i];
            double local_acc_x = 0.0, local_acc_y = 0.0, local_acc_z = 0.0;

            int stack[128];
            int stack_ptr = 0;
            stack[stack_ptr++] = 0;

            while (stack_ptr > 0) {
                int node_idx = stack[--stack_ptr];
                const BVHNode& node = bvh_nodes[node_idx];

                double aabb_dist_sq =
                    min_periodic_dist_sq(p1_x, node.bbox.min_x, node.bbox.max_x,
                                         domain_size) +
                    min_periodic_dist_sq(p1_y, node.bbox.min_y, node.bbox.max_y,
                                         domain_size) +
                    min_periodic_dist_sq(p1_z, node.bbox.min_z, node.bbox.max_z,
                                         domain_size);

                if (aabb_dist_sq > search_sq) continue;

                if (node.particle_idx != -1) {
                    int j = node.particle_idx;
                    if (static_cast<size_t>(j) == i) continue;

                    double dx = p1_x - pos_x[j];
                    if (dx > 0.5 * domain_size)
                        dx -= domain_size;
                    else if (dx < -0.5 * domain_size)
                        dx += domain_size;

                    double dy = p1_y - pos_y[j];
                    if (dy > 0.5 * domain_size)
                        dy -= domain_size;
                    else if (dy < -0.5 * domain_size)
                        dy += domain_size;

                    double dz = p1_z - pos_z[j];
                    if (dz > 0.5 * domain_size)
                        dz -= domain_size;
                    else if (dz < -0.5 * domain_size)
                        dz += domain_size;

                    dx = -dx;
                    dy = -dy;
                    dz = -dz;

                    double dist_sq = dx * dx + dy * dy + dz * dz;

                    if (use_pm && dist_sq > cutoff_sq) continue;

                    double pp_dist_sq = dist_sq + soft_sq;
                    double pp_dist = std::sqrt(pp_dist_sq);
                    double a_pp = G * node.mass / pp_dist_sq;

                    if (use_pm) {
                        double r_scaled = pp_dist / (2.0 * r_s);
                        a_pp *= (std::erfc(r_scaled) +
                                 (pp_dist / (std::sqrt(M_PI) * r_s)) *
                                     std::exp(-r_scaled * r_scaled));
                    }

                    local_acc_x += a_pp * dx / pp_dist;
                    local_acc_y += a_pp * dy / pp_dist;
                    local_acc_z += a_pp * dz / pp_dist;
                } else {
                    stack[stack_ptr++] = node.left_child;
                    stack[stack_ptr++] = node.right_child;
                }
            }

            acc_x[i] += local_acc_x;
            acc_y[i] += local_acc_y;
            acc_z[i] += local_acc_z;
        }
    }
}