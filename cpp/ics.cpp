#include "ics.h"

#include <random>

#include "integrator.h"
#include "particles.h"
#include "pocketfft_hdronly.h"

struct ZeldovichField {
    std::vector<double> dx;
    std::vector<double> dy;
    std::vector<double> dz;
    double std_x;
    double std_y;
    double std_z;
};

static ZeldovichField compute_zeldovich_field(const Config& config) {
    int M = config.MESH_SIZE;
    size_t M3_real = static_cast<size_t>(M) * M * M;
    size_t M3_complex = static_cast<size_t>(M) * M * (M / 2 + 1);

    std::vector<double> real_space_random_field(M3_real);
    std::default_random_engine generator(config.SEED);
    std::normal_distribution<double> distribution(0.0, 1.0);
    for (auto& val : real_space_random_field) {
        val = distribution(generator);
    }

    pocketfft::shape_t shape_ic = {(size_t)M, (size_t)M, (size_t)M};
    pocketfft::stride_t stride_r_ic = {
        static_cast<ptrdiff_t>(M * M * sizeof(double)),
        static_cast<ptrdiff_t>(M * sizeof(double)), sizeof(double)};
    pocketfft::stride_t stride_c_ic = {
        static_cast<ptrdiff_t>(M * (M / 2 + 1) * sizeof(std::complex<double>)),
        static_cast<ptrdiff_t>((M / 2 + 1) * sizeof(std::complex<double>)),
        sizeof(std::complex<double>)};

    std::vector<std::complex<double>> random_k(M3_complex);
    pocketfft::r2c(shape_ic, stride_r_ic, stride_c_ic, {0, 1, 2}, true,
                   real_space_random_field.data(), random_k.data(), 1.0);

    std::vector<std::complex<double>> disp_x_k(M3_complex);
    std::vector<std::complex<double>> disp_y_k(M3_complex);
    std::vector<std::complex<double>> disp_z_k(M3_complex);

    for (int i = 0; i < M; ++i) {
        for (int j = 0; j < M; ++j) {
            for (int k = 0; k < M / 2 + 1; ++k) {
                if (i == 0 && j == 0 && k == 0) {
                    disp_x_k[0] = {0, 0};
                    disp_y_k[0] = {0, 0};
                    disp_z_k[0] = {0, 0};
                    continue;
                }

                auto kx_freq = static_cast<double>((i < M / 2) ? i : (i - M));
                auto ky_freq = static_cast<double>((j < M / 2) ? j : (j - M));
                auto kz_freq = static_cast<double>(k);

                double kx = kx_freq * 2.0 * M_PI / config.DOMAIN_SIZE;
                double ky = ky_freq * 2.0 * M_PI / config.DOMAIN_SIZE;
                double kz = kz_freq * 2.0 * M_PI / config.DOMAIN_SIZE;
                double k2 = kx * kx + ky * ky + kz * kz;

                size_t idx = static_cast<size_t>(i) * M * (M / 2 + 1) +
                             static_cast<size_t>(j) * (M / 2 + 1) +
                             static_cast<size_t>(k);

                double power_spectrum_sqrt =
                    sqrt(pow(k2, config.INITIAL_POWER_SPECTRUM_INDEX / 2.0));
                std::complex<double> delta_k =
                    random_k[idx] * power_spectrum_sqrt;
                std::complex<double> phi_k = -delta_k / k2;

                disp_x_k[idx] = std::complex<double>(0, -1) * kx * phi_k;
                disp_y_k[idx] = std::complex<double>(0, -1) * ky * phi_k;
                disp_z_k[idx] = std::complex<double>(0, -1) * kz * phi_k;
            }
        }
    }

    ZeldovichField field;
    field.dx.resize(M3_real);
    field.dy.resize(M3_real);
    field.dz.resize(M3_real);

    pocketfft::c2r(shape_ic, stride_c_ic, stride_r_ic, {0, 1, 2}, false,
                   disp_x_k.data(), field.dx.data(), 1.0);
    pocketfft::c2r(shape_ic, stride_c_ic, stride_r_ic, {0, 1, 2}, false,
                   disp_y_k.data(), field.dy.data(), 1.0);
    pocketfft::c2r(shape_ic, stride_c_ic, stride_r_ic, {0, 1, 2}, false,
                   disp_z_k.data(), field.dz.data(), 1.0);

    double norm_ic = 1.0 / static_cast<double>(M3_real);
    field.std_x = 0;
    field.std_y = 0;
    field.std_z = 0;

    for (size_t i = 0; i < M3_real; ++i) {
        field.dx[i] *= norm_ic;
        field.dy[i] *= norm_ic;
        field.dz[i] *= norm_ic;
        field.std_x += field.dx[i] * field.dx[i];
        field.std_y += field.dy[i] * field.dy[i];
        field.std_z += field.dz[i] * field.dz[i];
    }
    field.std_x = sqrt(field.std_x / M3_real);
    field.std_y = sqrt(field.std_y / M3_real);
    field.std_z = sqrt(field.std_z / M3_real);

    return field;
}

void initialize_dm(SimState& state, const Config& config,
                   const ZeldovichField& zf) {
    state.dm.particles.clear();
    int M = config.MESH_SIZE;
    double cell_size = config.DOMAIN_SIZE / M;
    int N_part = config.N_PER_SIDE;
    double spacing = config.DOMAIN_SIZE / N_part;

    for (int i = 0; i < N_part; ++i) {
        for (int j = 0; j < N_part; ++j) {
            for (int k = 0; k < N_part; ++k) {
                double qx = (i + 0.5) * spacing;
                double qy = (j + 0.5) * spacing;
                double qz = (k + 0.5) * spacing;

                int ix = static_cast<int>(qx / cell_size) % M;
                int iy = static_cast<int>(qy / cell_size) % M;
                int iz = static_cast<int>(qz / cell_size) % M;

                size_t idx = static_cast<size_t>(ix) * M * M +
                             static_cast<size_t>(iy) * M +
                             static_cast<size_t>(iz);

                double dx =
                    (zf.dx[idx] / zf.std_x) * state.scale_factor * config.SIGMA_RMS;
                double dy =
                    (zf.dy[idx] / zf.std_y) * state.scale_factor * config.SIGMA_RMS;
                double dz =
                    (zf.dz[idx] / zf.std_z) * state.scale_factor * config.SIGMA_RMS;

                Particle p;
                p.pos.x =
                    fmod(qx + dx + config.DOMAIN_SIZE, config.DOMAIN_SIZE);
                p.pos.y =
                    fmod(qy + dy + config.DOMAIN_SIZE, config.DOMAIN_SIZE);
                p.pos.z =
                    fmod(qz + dz + config.DOMAIN_SIZE, config.DOMAIN_SIZE);

                p.vel.x =
                    config.STANDING_PARTICLES ? 0.0 : state.hubble_param * dx;
                p.vel.y =
                    config.STANDING_PARTICLES ? 0.0 : state.hubble_param * dy;
                p.vel.z =
                    config.STANDING_PARTICLES ? 0.0 : state.hubble_param * dz;
                p.mass = config.DM_PARTICLE_MASS;

                state.dm.particles.push_back(p);
            }
        }
    }
}

void initialize_gas(SimState& state, const Config& config,
                    const ZeldovichField& zf) {
    if (!config.USE_HYDRO) return;

    int M = config.MESH_SIZE;
    double cell_size = config.DOMAIN_SIZE / M;
    size_t M3_real = static_cast<size_t>(M) * M * M;

    double total_dm_mass = state.dm.particles.size() * config.DM_PARTICLE_MASS;
    double mass_ratio = config.GAS_TOTAL_MASS / total_dm_mass;

    auto& gas = state.gas;

    gas.density.data = state.dm.dm_rho.data * mass_ratio;
    gas.density.data = (gas.get_density().array() < 1e-12)
                           .select(1e-12, gas.get_density().data);

    double initial_internal_energy = 1e-6;

    for (size_t i = 0; i < M3_real; ++i) {
        double dx = (zf.dx[i] / zf.std_x) * state.scale_factor * config.SIGMA_RMS;
        double dy = (zf.dy[i] / zf.std_y) * state.scale_factor * config.SIGMA_RMS;
        double dz = (zf.dz[i] / zf.std_z) * state.scale_factor * config.SIGMA_RMS;

        double vx = config.STANDING_PARTICLES ? 0.0 : state.hubble_param * dx;
        double vy = config.STANDING_PARTICLES ? 0.0 : state.hubble_param * dy;
        double vz = config.STANDING_PARTICLES ? 0.0 : state.hubble_param * dz;

        gas.velocity_x.data[i] = vx;
        gas.velocity_y.data[i] = vy;
        gas.velocity_z.data[i] = vz;

        double rho = gas.get_density().data[i];
        gas.momentum_x.data[i] = rho * vx;
        gas.momentum_y.data[i] = rho * vy;
        gas.momentum_z.data[i] = rho * vz;

        double kin_energy = 0.5 * rho * (vx * vx + vy * vy + vz * vz);
        gas.energy.data[i] = (rho * initial_internal_energy) + kin_energy;
    }

    gas.update_primitive_variables();
}

SimState initialize_state(Config& config) {
    SimState state(config);
    state.total_time = 0;
    update_cosmology(state, config);

    // 1. Math Step
    ZeldovichField z_field = compute_zeldovich_field(config);

    // 2. Dark Matter Step
    initialize_dm(state, config, z_field);
    state.dm.bin_and_assign_mass(config);

    // 3. Gas Step
    initialize_gas(state, config, z_field);

    // 4. Forces Setup (UPDATED FOR NEW GRAVITY ARCHITECTURE)
    Grid3D total_rho = compute_gravitational_acceleration(
        state.gravity_x, state.gravity_y, state.gravity_z, state.gas, config,
        state.dm.dm_rho);

    std::vector<Vec3> pp_forces, pm_forces;
    if (config.USE_PM) {
        state.dm.interpolate_cic_forces(state.gravity_x, state.gravity_y,
                                        state.gravity_z, pm_forces, config);
    } else {
        pm_forces.assign(state.dm.particles.size(), {0.0, 0.0, 0.0});
    }

    if (config.USE_PP) {
        state.dm.compute_pp_forces(pp_forces, config);
    } else {
        pp_forces.assign(state.dm.particles.size(), {0.0, 0.0, 0.0});
    }

    // 5. Apply Initial Forces
    for (size_t i = 0; i < state.dm.particles.size(); ++i) {
        state.dm.particles[i].acc.x =
            (pp_forces[i].x + pm_forces[i].x) / state.dm.particles[i].mass;
        state.dm.particles[i].acc.y =
            (pp_forces[i].y + pm_forces[i].y) / state.dm.particles[i].mass;
        state.dm.particles[i].acc.z =
            (pp_forces[i].z + pm_forces[i].z) / state.dm.particles[i].mass;
    }

    return state;
}