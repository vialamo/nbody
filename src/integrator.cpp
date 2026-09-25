#include "integrator.h"

#include <omp.h>

// #include <iostream>

#include "math_utils.h"
#include "particles.h"
#include "pocketfft_hdronly.h"

void update_cosmology(SimState& state, const Config& config) {
    if (!config.expanding_universe) {
        state.scale_factor = 1.0;
        state.hubble_param = 0.0;
        return;
    }

    if (config.omega_lambda == 0.0) {
        // Einstein-de Sitter
        double t_start = std::pow(config.a_start, 1.5);
        double current_t = t_start + state.total_time;
        state.scale_factor = std::pow(current_t, 2.0 / 3.0);
        state.hubble_param = (2.0 / 3.0) / current_t;
    } else {
        // Lambda-CDM
        const double Om = config.omega_m;
        const double Ol = config.omega_lambda;
        const double H0 = 2.0 / (3.0 * std::sqrt(Om));
        double factor = std::sqrt(Ol / Om);

        // Reverse-engineer the initial time from start_a
        double t_start = (2.0 / (3.0 * H0 * std::sqrt(Ol))) *
                         std::asinh(factor * std::pow(config.a_start, 1.5));
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

double get_time_from_scale_factor(double a, const Config& config) {
    if (!config.expanding_universe) {
        // Without expansion, 'a' doesn't change
        return std::numeric_limits<double>::infinity();
    }

    if (config.omega_lambda == 0.0) {
        // Einstein-de Sitter
        return std::pow(a, 1.5);
    } else {
        // Lambda-CDM
        const double Om = config.omega_m;
        const double Ol = config.omega_lambda;
        const double H0 = 2.0 / (3.0 * std::sqrt(Om));
        double factor = std::sqrt(Ol / Om);

        // This is the inverse of the Lambda-CDM solver
        return (2.0 / (3.0 * H0 * std::sqrt(Ol))) *
               std::asinh(factor * std::pow(a, 1.5));
    }
}

void compute_PM_acceleration(SimState& state, const Config& config) {
    int N = config.mesh_size;
    const GasGrid& gas = *state.gas;
    const Grid3D& dm_rho = state.dm.get_rho();
    Grid3D& total_rho = state.total_rho;
    Grid3D& acc_x = state.pm_gravity_x;
    Grid3D& acc_y = state.pm_gravity_y;
    Grid3D& acc_z = state.pm_gravity_z;

    Grid3D& phi = state.phi;
    if (config.hydro_method == HydroMethod::MFM) {
        total_rho.data = dm_rho.data + state.mfm_gas->get_rho().data;
    } else if (config.hydro_method == HydroMethod::Eulerian) {
        total_rho.data = dm_rho.data + gas.get_density().data;
    } else {
        total_rho.data = dm_rho.data;
    }

    const pocketfft::shape_t shape = {(size_t)N, (size_t)N, (size_t)N};

    const pocketfft::stride_t stride_r = {
        static_cast<ptrdiff_t>((size_t)N * N * sizeof(double)),
        static_cast<ptrdiff_t>(N * sizeof(double)), sizeof(double)};
    const pocketfft::stride_t stride_c = {
        static_cast<ptrdiff_t>((size_t)N * (N / 2 + 1) *
                               sizeof(std::complex<double>)),
        static_cast<ptrdiff_t>((N / 2 + 1) * sizeof(std::complex<double>)),
        sizeof(std::complex<double>)};

    std::vector<std::complex<double>> rho_k((size_t)N * N * (N / 2 + 1));
    pocketfft::r2c(shape, stride_r, stride_c, {0, 1, 2}, true,
                   total_rho.raw_data(), rho_k.data(), 1.0, 0);

    std::vector<std::complex<double>> phi_k((size_t)N * N * (N / 2 + 1));

    // Gaussian smoothing scale for Fourier-Split PM
    const double r_s = config.PM_smoothing_cells * config.cell_size;

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

                double kx = kx_freq * 2.0 * M_PI / config.domain_size;
                double ky = ky_freq * 2.0 * M_PI / config.domain_size;
                double kz = kz_freq * 2.0 * M_PI / config.domain_size;
                double k2 = kx * kx + ky * ky + kz * kz;

                size_t idx = static_cast<size_t>(i) * N * (N / 2 + 1) +
                             static_cast<size_t>(j) * (N / 2 + 1) +
                             static_cast<size_t>(k);

                // Calculate the Gaussian filter
                double filter = std::exp(-k2 * r_s * r_s);

                // CIC Deconvolution Filter
                auto sinc = [](double freq_index, int mesh_size) {
                    if (freq_index == 0.0) return 1.0;
                    double arg = M_PI * freq_index / mesh_size;
                    return std::sin(arg) / arg;
                };

                double sinc_x = sinc(kx_freq, N);
                double sinc_y = sinc(ky_freq, N);
                double sinc_z = sinc(kz_freq, N);

                // CIC is applied twice (mass splatting + force gathering).
                // Each application is (sinc_x * sinc_y * sinc_z)^2.
                double cic_blur = std::pow(sinc_x * sinc_y * sinc_z, 4);

                // Apply the Gaussian and remove the CIC blur
                phi_k[idx] = rho_k[idx] * (-4.0 * M_PI * config.G / k2) *
                             (filter / cic_blur);
            }
        }
    }

    pocketfft::c2r(shape, stride_c, stride_r, {0, 1, 2}, false, phi_k.data(),
                   const_cast<double*>(phi.raw_data()), 1.0, 0);

    double norm = 1.0 / ((double)N * N * N);
    double a_inv3 =
        1.0 / (state.scale_factor * state.scale_factor * state.scale_factor);
    // Bake the comoving scaling into the finite difference operator
    double factor = (-norm / (2.0 * config.cell_size)) * a_inv3;

    // Only collapse the outer two loops so we can optimize the index math
#pragma omp parallel for collapse(2)
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            // Calculate i and j boundaries ONCE per column
            int ip1 = (i + 1 == N) ? 0 : i + 1;
            int im1 = (i == 0) ? N - 1 : i - 1;
            int jp1 = (j + 1 == N) ? 0 : j + 1;
            int jm1 = (j == 0) ? N - 1 : j - 1;

#pragma omp simd
            for (int k = 0; k < N; ++k) {
                int kp1 = (k + 1 == N) ? 0 : k + 1;
                int km1 = (k == 0) ? N - 1 : k - 1;

                acc_x(i, j, k) = (phi(ip1, j, k) - phi(im1, j, k)) * factor;
                acc_y(i, j, k) = (phi(i, jp1, k) - phi(i, jm1, k)) * factor;
                acc_z(i, j, k) = (phi(i, j, kp1) - phi(i, j, km1)) * factor;
            }
        }
    }
}

void apply_mesh_gas_gravity_kick(GasGrid& gas, const Grid3D& grav_x,
                                 const Grid3D& grav_y, const Grid3D& grav_z,
                                 double dt, double a, double H,
                                 const Config& config) {
    if (config.hydro_method != HydroMethod::Eulerian) return;

    gas.update_primitive_variables(a);

    Grid3D total_ax_gas(config.mesh_size), total_ay_gas(config.mesh_size),
        total_az_gas(config.mesh_size);
    total_ax_gas.data = grav_x.array() - (2 * H * gas.get_velocity_x().array());
    total_ay_gas.data = grav_y.array() - (2 * H * gas.get_velocity_y().array());
    total_az_gas.data = grav_z.array() - (2 * H * gas.get_velocity_z().array());

    Grid3D g_mom_x_source(config.mesh_size), g_mom_y_source(config.mesh_size),
        g_mom_z_source(config.mesh_size);
    g_mom_x_source.data = gas.get_density().array() * total_ax_gas.array();
    g_mom_y_source.data = gas.get_density().array() * total_ay_gas.array();
    g_mom_z_source.data = gas.get_density().array() * total_az_gas.array();

    // Calculate the change in momentum (dp)
    auto dp_x = g_mom_x_source.array() * dt;
    auto dp_y = g_mom_y_source.array() * dt;
    auto dp_z = g_mom_z_source.array() * dt;

    // Fetch current momentum
    const auto& p_x = gas.get_momentum_x().array();
    const auto& p_y = gas.get_momentum_y().array();
    const auto& p_z = gas.get_momentum_z().array();

    // Kinetic Energy work done by gravity: (P_new^2 - P_old^2) / (2*rho)
    Eigen::ArrayXd gravity_work_density =
        ((p_x + dp_x).square() + (p_y + dp_y).square() + (p_z + dp_z).square() -
         (p_x.square() + p_y.square() + p_z.square())) /
        (2.0 * gas.get_density().array());

    // Apply updates
    gas.momentum_x.array() += dp_x;
    gas.momentum_y.array() += dp_y;
    gas.momentum_z.array() += dp_z;
    gas.energy.array() += gravity_work_density;

    // Adiabatic Expansion Cooling (PdV Work)
    //
    // In physical units, adiabatic expansion causes the thermal energy to decay
    // at a rate of -3(gamma - 1)H.
    //
    // However, we track comoving velocity as v = dx/dt, meaning
    // physical kinetic energy scales as a^2 * v_code^2. To conserve total
    // energy when adding kinetic and internal arrays together, our comoving
    // internal energy must also absorb this scaling: E_phys = a^2 * E_code.
    //
    // Applying the chain rule (d/dt) to this a^2 transformation introduces an
    // extra -2H mathematical decay just to keep the arrays synchronized.
    // The total required numerical decay is -(3*gamma - 1)H.
    double expansion_factor = (3.0 * config.gamma - 1.0) * H * dt;
    Eigen::ArrayXd expansion_cooling =
        expansion_factor * gas.get_internal_energy().array();

    gas.energy.array() -= expansion_cooling;
    gas.internal_energy.array() -= expansion_cooling;

    // Accumulate global work for diagnostics
    // sum() adds up the energy density of all cells. Multiply by volume to get
    // total code energy.
    double step_grav_work = gravity_work_density.sum() * config.cell_volume;
    double step_exp_work = expansion_cooling.sum() * config.cell_volume;

    gas.add_gravitational_work(step_grav_work);
    gas.add_expansion_work(step_exp_work);
}

static void apply_dm_kick(ParticleSystem& dm, double dt, double a, double H,
                          const Config& config, bool is_kick1) {
    double step_grav_work = 0.0;
    double step_exp_work = 0.0;
    const size_t n = dm.num_particles;

    // dt passed to this function is half of the macro step
    double dt_macro = 2.0 * dt;

#pragma omp parallel for reduction(+ : step_grav_work, step_exp_work) \
    schedule(static)
    for (size_t i = 0; i < n; i++) {
        // Gate the kick based on the individual timeline boundaries
        if (config.individual_particle_timesteps) {
            double t_start = dm.t_end[i] - dm.dt_step[i];
            bool do_kick = false;

            if (is_kick1) {
                // Kick 1 only applies if the particle's block step just started
                if (std::abs(dm.global_time - dt_macro - t_start) < 1e-10)
                    do_kick = true;
            } else {
                // Kick 2 only applies if the particle's block step just ended
                if (std::abs(dm.global_time - dm.t_end[i]) < 1e-10)
                    do_kick = true;
            }
            if (!do_kick) continue;
        } else {
            if (!dm.is_active[i]) continue;
        }

        // Use the exact half-step of the individual particle
        double half_dt = dm.dt_step[i] / 2.0;

        // Calculate cosmological damping factor tailored to this specific step
        double drag_factor = 1.0;
        if (H > 0.0 && a > 0.0) {
            double a_next = a + a * H * half_dt;
            drag_factor =
                (a / a_next) * (a / a_next);  // (1 / (1 + H*half_dt))^2
        }

        double m = dm.mass[i];
        double vx_old = dm.vel_x[i];
        double vy_old = dm.vel_y[i];
        double vz_old = dm.vel_z[i];

        // Comoving Gravitational Acceleration
        double gx = dm.acc_x[i];
        double gy = dm.acc_y[i];
        double gz = dm.acc_z[i];

        // Intermediate velocity after pure gravity kick
        double vx_g = vx_old + gx * half_dt;
        double vy_g = vy_old + gy * half_dt;
        double vz_g = vz_old + gz * half_dt;

        // Gravitational Work (dW_grav = dK_grav)
        double ke_old =
            0.5 * m * (vx_old * vx_old + vy_old * vy_old + vz_old * vz_old);
        double ke_g = 0.5 * m * (vx_g * vx_g + vy_g * vy_g + vz_g * vz_g);
        step_grav_work += (ke_g - ke_old);

        // Final velocity after Hubble Drag decay
        double vx_new = vx_g * drag_factor;
        double vy_new = vy_g * drag_factor;
        double vz_new = vz_g * drag_factor;

        // Expansion Loss Work (dW_exp = dK_drag)
        double ke_new =
            0.5 * m * (vx_new * vx_new + vy_new * vy_new + vz_new * vz_new);
        step_exp_work += (ke_g - ke_new);

        dm.vel_x[i] = vx_new;
        dm.vel_y[i] = vy_new;
        dm.vel_z[i] = vz_new;
    }

    dm.accumulated_gravitational_work += step_grav_work;
    dm.accumulated_expansion_work += step_exp_work;
}

static void apply_dm_drift(ParticleSystem& dm, double dt, double domain_size) {
    double inv_domain = 1.0 / domain_size;
    const size_t n = dm.num_particles;

#pragma omp parallel for simd schedule(static)
    for (size_t i = 0; i < n; i++) {
        // Drift strictly from the time the particle was last moved
        double drift_dt = dm.global_time - dm.t_current[i];
        if (drift_dt <= 0.0) continue;

        double nx = dm.pos_x[i] + dm.vel_x[i] * drift_dt;
        double ny = dm.pos_y[i] + dm.vel_y[i] * drift_dt;
        double nz = dm.pos_z[i] + dm.vel_z[i] * drift_dt;

        // SIMD-friendly periodic boundary
        dm.pos_x[i] = nx - domain_size * std::floor(nx * inv_domain);
        dm.pos_y[i] = ny - domain_size * std::floor(ny * inv_domain);
        dm.pos_z[i] = nz - domain_size * std::floor(nz * inv_domain);

        // Mark the particle as drifted up to the current global clock
        dm.t_current[i] = dm.global_time;
    }
}

static void apply_gas_particle_gravity_kick(GasParticleSystem& gas, double dt,
                                            double a, double H,
                                            const Config& config,
                                            bool is_kick1) {
    double step_grav_work = 0.0;
    double step_exp_work = 0.0;
    const size_t n = gas.num_particles;
    double dt_macro = 2.0 * dt;

#pragma omp parallel for reduction(+ : step_grav_work, step_exp_work) \
    schedule(static)
    for (size_t i = 0; i < n; i++) {
        if (config.individual_particle_timesteps) {
            double t_start = gas.t_end[i] - gas.dt_step[i];
            bool do_kick = false;
            if (is_kick1) {
                if (std::abs(gas.global_time - dt_macro - t_start) < 1e-10)
                    do_kick = true;
            } else {
                if (std::abs(gas.global_time - gas.t_end[i]) < 1e-10)
                    do_kick = true;
            }
            if (!do_kick) continue;
        } else {
            if (!gas.is_active[i]) continue;
        }

        double half_dt = gas.dt_step[i] / 2.0;
        KickWork work = gas.kick_particle_gravity(i, half_dt, a, H, config);

        step_grav_work += work.grav_work;
        step_exp_work += work.exp_work;
    }

    gas.accumulated_gravitational_work += step_grav_work;
    gas.accumulated_expansion_work += step_exp_work;
}

static void apply_gas_particle_hydro_kick(GasParticleSystem& gas, double dt,
                                          double a, const Config& config,
                                          bool is_kick1) {
    const size_t n = gas.num_particles;
    double step_hydro_exp_work = 0.0;
    double dt_macro = 2.0 * dt;

#pragma omp parallel for reduction(+ : step_hydro_exp_work) schedule(static)
    for (size_t i = 0; i < n; i++) {
        if (config.individual_particle_timesteps) {
            double t_start = gas.t_end[i] - gas.dt_step[i];
            bool do_kick = false;
            if (is_kick1) {
                if (std::abs(gas.global_time - dt_macro - t_start) < 1e-10)
                    do_kick = true;
            } else {
                if (std::abs(gas.global_time - gas.t_end[i]) < 1e-10)
                    do_kick = true;
            }
            if (!do_kick) continue;
        } else {
            if (!gas.is_active[i]) continue;
        }

        double half_dt = gas.dt_step[i] / 2.0;
        KickWork work = gas.kick_particle_hydro(i, half_dt, config);

        step_hydro_exp_work += work.hydro_exp_work;
    }

    gas.accumulated_expansion_work -= step_hydro_exp_work;
}

static void apply_gravity_kick(SimState& state, double dt, double a, double H,
                               const Config& config, bool is_kick1) {
    // Forward the config and kick state to the updated DM function
    apply_dm_kick(state.dm, dt, a, H, config, is_kick1);

    if (config.hydro_method == HydroMethod::Eulerian) {
        apply_mesh_gas_gravity_kick(*state.gas, state.pm_gravity_x,
                                    state.pm_gravity_y, state.pm_gravity_z, dt,
                                    a, H, config);
    } else if (config.hydro_method == HydroMethod::MFM) {
        apply_gas_particle_gravity_kick(*state.mfm_gas, dt, a, H, config,
                                        is_kick1);
    }
}

static void apply_gas_particle_drift(GasParticleSystem& gas, double dt,
                                     double domain_size) {
    double inv_domain = 1.0 / domain_size;
    const size_t n = gas.num_particles;

    // Drift from t_current to global_time
#pragma omp parallel for simd schedule(static)
    for (size_t i = 0; i < n; i++) {
        double drift_dt = gas.global_time - gas.t_current[i];
        if (drift_dt <= 0.0) continue;  // Skip if already drifted

        double nx = gas.pos_x[i] + gas.vel_x[i] * drift_dt;
        double ny = gas.pos_y[i] + gas.vel_y[i] * drift_dt;
        double nz = gas.pos_z[i] + gas.vel_z[i] * drift_dt;

        // Periodic boundaries
        gas.pos_x[i] = nx - domain_size * std::floor(nx * inv_domain);
        gas.pos_y[i] = ny - domain_size * std::floor(ny * inv_domain);
        gas.pos_z[i] = nz - domain_size * std::floor(nz * inv_domain);

        gas.t_current[i] = gas.global_time;  // Mark as drifted
    }
}

static void calculate_max_acceleration(ParticleSystem& dm) {
    double local_max = 1e-9;
    const size_t n = dm.num_particles;
    const double* ax = dm.acc_x.data();
    const double* ay = dm.acc_y.data();
    const double* az = dm.acc_z.data();

#pragma omp parallel for reduction(max : local_max) schedule(static)
    for (size_t i = 0; i < n; ++i) {
        double accel_sq = ax[i] * ax[i] + ay[i] * ay[i] + az[i] * az[i];
        local_max = (accel_sq > local_max) ? accel_sq : local_max;
    }

    dm.max_accel_sq = local_max;
}

static void calculate_max_acceleration(GasParticleSystem& gas) {
    double local_max = 1e-9;
    const size_t n = gas.num_particles;
    const double* ax = gas.acc_x.data();
    const double* ay = gas.acc_y.data();
    const double* az = gas.acc_z.data();

#pragma omp parallel for reduction(max : local_max) schedule(static)
    for (size_t i = 0; i < n; ++i) {
        double accel_sq = ax[i] * ax[i] + ay[i] * ay[i] + az[i] * az[i];
        local_max = (accel_sq > local_max) ? accel_sq : local_max;
    }
    gas.max_accel_sq = local_max;
}

static void update_softening(SimState& state, Config& config) {
    double a = state.scale_factor;
    double cap_a = config.physical_softening_cap_a;
    double current_comoving_softening = config.base_comoving_softening;

    if (a > cap_a) {
        // Shrink comoving softening so physical softening stays constant
        current_comoving_softening =
            config.base_comoving_softening * (cap_a / a);
    }
    config.softening_squared =
        current_comoving_softening * current_comoving_softening;
}

void compute_forces(SimState& state, Config& config, Diagnostics& diag) {
    update_softening(state, config);

    // Build trees
    {
        ScopedTimer timer(diag, TimerRegion::Tree);
        state.dm.build_lbvh(config);
        if (config.hydro_method == HydroMethod::MFM) {
            state.mfm_gas->build_lbvh(config);
        }
    }

    // PM
    {
        ScopedTimer pm_timer(diag, TimerRegion::PM);
        state.dm.bin_and_assign_mass(config);

        if (config.hydro_method == HydroMethod::MFM) {
            state.mfm_gas->bin_and_assign_mass(config);
            state.total_rho.data =
                state.dm.get_rho().data + state.mfm_gas->get_rho().data;
        } else if (config.hydro_method == HydroMethod::Eulerian) {
            state.total_rho.data =
                state.dm.get_rho().data + state.gas->get_density().data;
        } else {
            state.total_rho.data = state.dm.get_rho().data;
        }

        if (config.use_PM) {
            compute_PM_acceleration(state, config);

            state.dm.interpolate_cic_forces(state.pm_gravity_x,
                                            state.pm_gravity_y,
                                            state.pm_gravity_z, config);
            if (config.hydro_method == HydroMethod::MFM) {
                state.mfm_gas->interpolate_cic_forces(
                    state.pm_gravity_x, state.pm_gravity_y, state.pm_gravity_z,
                    config);
            }
        } else {
            state.pm_gravity_x.setZero();
            state.pm_gravity_y.setZero();
            state.pm_gravity_z.setZero();

            std::fill(state.dm.acc_x.begin(), state.dm.acc_x.end(), 0.0);
            std::fill(state.dm.acc_y.begin(), state.dm.acc_y.end(), 0.0);
            std::fill(state.dm.acc_z.begin(), state.dm.acc_z.end(), 0.0);
            if (config.hydro_method == HydroMethod::MFM) {
                GasParticleSystem& gas = *state.mfm_gas;
                std::fill(gas.acc_x.begin(), gas.acc_x.end(), 0.0);
                std::fill(gas.acc_y.begin(), gas.acc_y.end(), 0.0);
                std::fill(gas.acc_z.begin(), gas.acc_z.end(), 0.0);
            }
        }
    }

    // PP GRAVITY
    if (config.use_PP) {
        ScopedTimer pp_timer(diag, TimerRegion::PP);
        state.dm.compute_and_add_pp_forces(state.scale_factor, config,
                                           diag);  // DM-DM

        if (config.hydro_method == HydroMethod::MFM) {
            state.mfm_gas->compute_and_add_pp_forces(state.scale_factor, config,
                                                     diag);  // Gas-Gas
            state.mfm_gas->compute_cross_pp_forces(state.scale_factor, state.dm,
                                                   config,
                                                   diag);  // Gas-DM and DM-Gas

            // Gas pulling on DM
            state.dm.compute_cross_pp_forces(state.scale_factor, *state.mfm_gas,
                                             config, diag);
        }
    }

    // MFM HYDRODYNAMICS
    if (config.hydro_method == HydroMethod::MFM) {
        ScopedTimer hydro_timer(diag, TimerRegion::Hydro);
        state.mfm_gas->compute_density_and_h(config, state.dm);
    }

    // Finalize
    {
        ScopedTimer pp_timer(diag, TimerRegion::PP);
        calculate_max_acceleration(state.dm);
        if (config.hydro_method == HydroMethod::MFM) {
            calculate_max_acceleration(*state.mfm_gas);
        }
    }
}

void KDK_step(SimState& state, TimestepInfo& ts, Config& config,
              Diagnostics& diag) {
    double dt = ts.dt_macro;

    // Clock sync and wake-ups
    state.dm.sync_and_activate(dt, config);
    if (config.hydro_method == HydroMethod::MFM) {
        state.mfm_gas->sync_and_activate(dt, state.scale_factor,
                                         state.hubble_param, config);
    }

    // KICK 1 (Half step)
    apply_gravity_kick(state, dt / 2.0, state.scale_factor, state.hubble_param,
                       config, true);
    if (config.hydro_method == HydroMethod::MFM) {
        // Note: apply_gas_particle_hydro_kick ignores the (dt / 2.0) passed
        // here and uses gas.dt_step[i] / 2.0 internally
        apply_gas_particle_hydro_kick(*state.mfm_gas, dt / 2.0,
                                      state.scale_factor, config, true);
    }

    // Approximate the scale factor at the half-step (t + dt/2)
    double mid_a = state.scale_factor * (1.0 + 0.5 * state.hubble_param * dt);
    double mid_H = state.hubble_param;

    // DRIFT
    apply_dm_drift(state.dm, dt, config.domain_size);

    if (config.hydro_method == HydroMethod::Eulerian) {
        {
            ScopedTimer hydro_timer(diag, TimerRegion::Hydro);
            state.gas->hydro_step(dt, mid_a);
        }
        if (config.enable_cooling) {
            ScopedTimer cooling_timer(diag, TimerRegion::Cool);
            state.gas->apply_cooling(dt, mid_a, state.cooling);
            diag.add_substeps(SubstepCounter::Cool,
                              state.gas->get_cooling_total_cycles());
        }
    } else if (config.hydro_method == HydroMethod::MFM) {
        {
            ScopedTimer hydro_timer(diag, TimerRegion::Hydro);

            // apply_gas_particle_drift automatically drifts to global_time
            apply_gas_particle_drift(*state.mfm_gas, dt, config.domain_size);
            // These loops skip inactive particles
            state.mfm_gas->compute_density_and_h(config, state.dm);
            state.mfm_gas->hydro_step(config, mid_a, mid_H, dt);
        }

        if (config.enable_cooling) {
            ScopedTimer cooling_timer(diag, TimerRegion::Cool);
            state.mfm_gas->apply_cooling(dt, mid_a, config, state.cooling);
            diag.add_substeps(SubstepCounter::Cool,
                              state.mfm_gas->cooling_total_cycles);
        }
    }

    // UPDATE COSMOLOGY to t + dt
    state.total_time += dt;
    update_cosmology(state, config);
    compute_forces(state, config, diag);

    // KICK 2 (Half step)
    apply_gravity_kick(state, dt / 2.0, state.scale_factor, state.hubble_param,
                       config, false);

    if (config.hydro_method == HydroMethod::MFM) {
        apply_gas_particle_hydro_kick(*state.mfm_gas, dt / 2.0,
                                      state.scale_factor, config, false);
    }

    // Time bin assignment
    if (config.hydro_method == HydroMethod::MFM) {
        state.mfm_gas->update_particle_timesteps(
            config.fixed_dt, state.scale_factor, config, state.cooling);
    }
    state.dm.update_particle_timesteps(config.fixed_dt, config);
}