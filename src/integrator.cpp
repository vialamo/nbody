#include "integrator.h"

#include <omp.h>

#include <chrono>

#include "particles.h"
#include "pocketfft_hdronly.h"

void update_cosmology(SimState& state, const Config& config) {
    if (!config.EXPANDING_UNIVERSE) {
        state.scale_factor = 1.0;
        state.hubble_param = 0.0;
        return;
    }

    if (config.OMEGA_LAMBDA == 0.0) {
        // Einstein-de Sitter
        double t_start = std::pow(config.START_A, 1.5);
        double current_cosmo_time = t_start + state.total_time;
        state.scale_factor = std::pow(current_cosmo_time, 2.0 / 3.0);
        state.hubble_param = (2.0 / 3.0) / current_cosmo_time;
    } else {
        // Lambda-CDM
        const double Om = config.OMEGA_M;
        const double Ol = config.OMEGA_LAMBDA;
        const double H0 = 2.0 / (3.0 * std::sqrt(Om));
        double factor = std::sqrt(Ol / Om);

        // Reverse-engineer the physical time of the Big Bang from START_A
        double t_start = (2.0 / (3.0 * H0 * std::sqrt(Ol))) *
                         std::asinh(factor * std::pow(config.START_A, 1.5));

        // Add elapsed simulation time
        double current_t = t_start + state.total_time;

        // Calculate current scale factor
        double sinh_term = std::sinh(1.5 * H0 * std::sqrt(Ol) * current_t);
        state.scale_factor =
            std::pow(Om / Ol, 1.0 / 3.0) * std::pow(sinh_term, 2.0 / 3.0);

        // Calculate current Hubble parameter H(a)
        double a3 = std::pow(state.scale_factor, 3.0);
        state.hubble_param = H0 * std::sqrt(Om / a3 + Ol);
    }
}

void compute_gravitational_acceleration(SimState& state, const Config& config) {
    int N = config.MESH_SIZE;
    const GasGrid& gas = state.gas;
    const Grid3D& dm_rho = state.dm.get_rho();
    Grid3D& total_rho = state.total_rho;
    Grid3D& acc_x = state.gravity_x;
    Grid3D& acc_y = state.gravity_y;
    Grid3D& acc_z = state.gravity_z;

    Grid3D& phi = state.phi;
    total_rho.data =
        dm_rho.data + (config.USE_HYDRO ? gas.get_density().data
                                        : Eigen::VectorXd::Zero(N * N * N));

    pocketfft::shape_t shape = {(size_t)N, (size_t)N, (size_t)N};

    pocketfft::stride_t stride_r = {
        static_cast<ptrdiff_t>((size_t)N * N * sizeof(double)),
        static_cast<ptrdiff_t>(N * sizeof(double)), sizeof(double)};
    pocketfft::stride_t stride_c = {
        static_cast<ptrdiff_t>((size_t)N * (N / 2 + 1) *
                               sizeof(std::complex<double>)),
        static_cast<ptrdiff_t>((N / 2 + 1) * sizeof(std::complex<double>)),
        sizeof(std::complex<double>)};

    std::vector<std::complex<double>> rho_k((size_t)N * N * (N / 2 + 1));
    pocketfft::r2c(shape, stride_r, stride_c, {0, 1, 2}, true,
                   total_rho.raw_data(), rho_k.data(), 1.0, 0);

    std::vector<std::complex<double>> phi_k((size_t)N * N * (N / 2 + 1));

#pragma omp parallel for collapse(3)
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            for (int k = 0; k < N / 2 + 1; ++k) {
                if (i == 0 && j == 0 && k == 0) {
                    phi_k[0] = {0.0, 0.0};
                    continue;
                }

                double kx_freq = static_cast<double>((i < N / 2) ? i : (i - N));
                double ky_freq = static_cast<double>((j < N / 2) ? j : (j - N));
                double kz_freq = static_cast<double>(k);

                double kx = kx_freq * 2.0 * M_PI / config.DOMAIN_SIZE;
                double ky = ky_freq * 2.0 * M_PI / config.DOMAIN_SIZE;
                double kz = kz_freq * 2.0 * M_PI / config.DOMAIN_SIZE;
                double k2 = kx * kx + ky * ky + kz * kz;

                size_t idx = static_cast<size_t>(i) * N * (N / 2 + 1) +
                             static_cast<size_t>(j) * (N / 2 + 1) +
                             static_cast<size_t>(k);

                phi_k[idx] = rho_k[idx] * (-4.0 * M_PI * config.G / k2);
            }
        }
    }

    pocketfft::c2r(shape, stride_c, stride_r, {0, 1, 2}, false, phi_k.data(),
                   const_cast<double*>(phi.raw_data()), 1.0, 0);

    double norm = 1.0 / ((double)N * N * N);
    double factor = -norm / (2.0 * config.CELL_SIZE);

#pragma omp parallel for collapse(3)
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            for (int k = 0; k < N; ++k) {
                acc_x(i, j, k) =
                    (phi((i + 1) % N, j, k) - phi((i - 1 + N) % N, j, k)) *
                    factor;
                acc_y(i, j, k) =
                    (phi(i, (j + 1) % N, k) - phi(i, (j - 1 + N) % N, k)) *
                    factor;
                acc_z(i, j, k) =
                    (phi(i, j, (k + 1) % N) - phi(i, j, (k - 1 + N) % N)) *
                    factor;
            }
        }
    }
}

void apply_gas_kick(GasGrid& gas, const Grid3D& grav_x, const Grid3D& grav_y,
                    const Grid3D& grav_z, double dt, double a, double H,
                    const Config& config) {
    if (!config.USE_HYDRO) return;

    gas.update_primitive_variables();
    double a3 = a * a * a;

    Grid3D total_ax_gas(config.MESH_SIZE), total_ay_gas(config.MESH_SIZE),
        total_az_gas(config.MESH_SIZE);
    total_ax_gas.data =
        (grav_x.array() / a3) - (2 * H * gas.get_velocity_x().array());
    total_ay_gas.data =
        (grav_y.array() / a3) - (2 * H * gas.get_velocity_y().array());
    total_az_gas.data =
        (grav_z.array() / a3) - (2 * H * gas.get_velocity_z().array());

    Grid3D g_mom_x_source(config.MESH_SIZE), g_mom_y_source(config.MESH_SIZE),
        g_mom_z_source(config.MESH_SIZE);
    g_mom_x_source.data = gas.get_density().array() * total_ax_gas.array();
    g_mom_y_source.data = gas.get_density().array() * total_ay_gas.array();
    g_mom_z_source.data = gas.get_density().array() * total_az_gas.array();

    Grid3D power_density(config.MESH_SIZE);
    power_density.data = gas.get_velocity_x().array() * g_mom_x_source.array() +
                         gas.get_velocity_y().array() * g_mom_y_source.array() +
                         gas.get_velocity_z().array() * g_mom_z_source.array();

    gas.momentum_x.array() += g_mom_x_source.array() * dt;
    gas.momentum_y.array() += g_mom_y_source.array() * dt;
    gas.momentum_z.array() += g_mom_z_source.array() * dt;
    gas.energy.array() += power_density.array() * dt;
}

static void apply_dm_kick(std::vector<Particle>& particles, double dt, double a,
                          double H) {
    double a3 = a * a * a;
    for (auto& p : particles) {
        double total_ax_p = (p.acc.x / a3) - (2 * H * p.vel.x);
        double total_ay_p = (p.acc.y / a3) - (2 * H * p.vel.y);
        double total_az_p = (p.acc.z / a3) - (2 * H * p.vel.z);
        p.vel.x += total_ax_p * dt;
        p.vel.y += total_ay_p * dt;
        p.vel.z += total_az_p * dt;
    }
}

static void apply_dm_drift(std::vector<Particle>& particles, double dt,
                           double domain_size) {
    for (auto& p : particles) {
        p.pos.x = fmod(p.pos.x + p.vel.x * dt + domain_size, domain_size);
        p.pos.y = fmod(p.pos.y + p.vel.y * dt + domain_size, domain_size);
        p.pos.z = fmod(p.pos.z + p.vel.z * dt + domain_size, domain_size);
    }
}

static void update_particles_accel(ParticleSystem& dm,
                                   const std::vector<Vec3>& pm_forces,
                                   const std::vector<Vec3>& pp_forces) {
    dm.max_accel_sq = 1e-9;
    for (size_t i = 0; i < dm.particles.size(); ++i) {
        auto& p = dm.particles[i];
        p.acc.x = (pp_forces[i].x + pm_forces[i].x) / p.mass;
        p.acc.y = (pp_forces[i].y + pm_forces[i].y) / p.mass;
        p.acc.z = (pp_forces[i].z + pm_forces[i].z) / p.mass;
        double accel_sq =
            p.acc.x * p.acc.x + p.acc.y * p.acc.y + p.acc.z * p.acc.z;
        if (accel_sq > dm.max_accel_sq) {
            dm.max_accel_sq = accel_sq;
        }
    }
}

std::map<std::string, double> KDK_step(SimState& state, double dt,
                                       Config& config) {
    std::map<std::string, double> timings;
    auto start_time = std::chrono::high_resolution_clock::now();
    auto end_time = start_time;

    update_cosmology(state, config);

    // KICK 1 (Half step)
    apply_gas_kick(state.gas, state.gravity_x, state.gravity_y, state.gravity_z,
                   dt / 2.0, state.scale_factor, state.hubble_param, config);
    apply_dm_kick(state.dm.particles, dt / 2.0, state.scale_factor,
                  state.hubble_param);

    // DRIFT
    apply_dm_drift(state.dm.particles, dt, config.DOMAIN_SIZE);

    // HYDRO DYNAMICS
    start_time = std::chrono::high_resolution_clock::now();
    if (config.USE_HYDRO) {
        state.gas.hydro_step(dt);
    }
    end_time = std::chrono::high_resolution_clock::now();
    timings["hydro"] =
        std::chrono::duration_cast<std::chrono::duration<double>>(end_time -
                                                                  start_time)
            .count();

    // UPDATE COSMOLOGY to t + dt
    state.total_time += dt;
    update_cosmology(state, config);

    // COMPUTE FORCES (Update Accelerations)
    start_time = std::chrono::high_resolution_clock::now();
    state.dm.bin_and_assign_mass(config);
    compute_gravitational_acceleration(state, config);

    std::vector<Vec3> pm_forces(state.dm.particles.size(), {0.0, 0.0, 0.0});
    if (config.USE_PM) {
        state.dm.interpolate_cic_forces(state.gravity_x, state.gravity_y,
                                        state.gravity_z, pm_forces, config);
    }
    end_time = std::chrono::high_resolution_clock::now();
    timings["pm"] = std::chrono::duration_cast<std::chrono::duration<double>>(
                        end_time - start_time)
                        .count();

    start_time = std::chrono::high_resolution_clock::now();
    std::vector<Vec3> pp_forces(state.dm.particles.size(), {0.0, 0.0, 0.0});
    if (config.USE_PP) {
        state.dm.compute_pp_forces(pp_forces, config);
    }
    end_time = std::chrono::high_resolution_clock::now();
    timings["pp"] = std::chrono::duration_cast<std::chrono::duration<double>>(
                        end_time - start_time)
                        .count();

    // Assign final accelerations back to particles
    if (!config.STANDING_PARTICLES) {
        update_particles_accel(state.dm, pm_forces, pp_forces);
    }

    // KICK 2 (Half step)
    apply_gas_kick(state.gas, state.gravity_x, state.gravity_y, state.gravity_z,
                   dt / 2.0, state.scale_factor, state.hubble_param, config);
    apply_dm_kick(state.dm.particles, dt / 2.0, state.scale_factor,
                  state.hubble_param);

    return timings;
}