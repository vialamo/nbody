#include "ics.h"

#include <random>

#include "integrator.h"
#include "math_utils.h"
#include "particles.h"
#include "pocketfft_hdronly.h"

// The pure shape of the BBKS power spectrum (Unnormalized)
static double unnormalized_pk(double k_h, const Config& config) {
    // k_h is the wavenumber in units of h Mpc^-1
    // Cosmological Shape Parameter (Gamma)
    double Gamma = config.OMEGA_M * config.HUBBLE_PARAM;

    // 'q' parameter for BBKS, scaled for Mpc/h units
    double q = k_h / Gamma;

    // BBKS Polynomial
    double T_k = log(1.0 + 2.34 * q) / (2.34 * q) *
                 pow(1.0 + 3.89 * q + pow(16.1 * q, 2) + pow(5.46 * q, 3) +
                         pow(6.71 * q, 4),
                     -0.25);

    // P(k) = k^n_s * T(k)^2
    return pow(k_h, config.SPECTRAL_INDEX) * T_k * T_k;
}

// Fourier transform of a spherical Top-Hat filter
static double window_tophat(double kR) {
    // Protect against division by zero at k=0 using Taylor expansion
    if (kR < 1e-4) return 1.0 - (kR * kR / 10.0);
    return 3.0 * (sin(kR) - kR * cos(kR)) / (kR * kR * kR);
}

ZeldovichField compute_zeldovich_field(double scale_factor,
                                       const Config& config) {
    double R = 8.0;  // The standard 8 Mpc/h scale
    double unnorm_variance = 0.0;
    double dk = 0.001;  // Integration step size

    // Numerically integrate from very large scales to very small scales
    for (double k = 0.001; k < 100.0; k += dk) {
        double pk = unnormalized_pk(k, config);
        double w = window_tophat(k * R);
        unnorm_variance += pk * w * w * k * k * dk;
    }
    unnorm_variance /= (2.0 * M_PI * M_PI);

    // The master normalization constant
    double A = (config.SIGMA_8 * config.SIGMA_8) / unnorm_variance;

    int M = config.MESH_SIZE;
    size_t M3_real = static_cast<size_t>(M) * M * M;
    size_t M3_complex = static_cast<size_t>(M) * M * (M / 2 + 1);

    double box_vol = pow(config.BOX_SIZE_MPC, 3.0);
    double amplitude_scaling =
        std::sqrt(static_cast<double>(M3_real) / box_vol);

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

                size_t idx = static_cast<size_t>(i) * M * (M / 2 + 1) +
                             static_cast<size_t>(j) * (M / 2 + 1) +
                             static_cast<size_t>(k);

                double kx = kx_freq * 2.0 * M_PI / config.BOX_SIZE_MPC;
                double ky = ky_freq * 2.0 * M_PI / config.BOX_SIZE_MPC;
                double kz = kz_freq * 2.0 * M_PI / config.BOX_SIZE_MPC;
                double k_mag = sqrt(kx * kx + ky * ky + kz * kz);
                double k_h = k_mag / config.HUBBLE_PARAM;

                // The BBKS Transfer Function T(k)
                double power_spectrum_sqrt = 0.0;
                if (k_mag > 0.0) {
                    // Multiply the unnormalized shape by A
                    double true_pk = A * unnormalized_pk(k_h, config);
                    power_spectrum_sqrt = sqrt(true_pk) * amplitude_scaling;
                }

                // Apply the scaling to the Fourier amplitudes
                std::complex<double> delta_k =
                    random_k[idx] * power_spectrum_sqrt;

                // Convert density contrast to Zel'dovich displacement
                // potential
                double code_kx = kx_freq * 2.0 * M_PI / config.DOMAIN_SIZE;
                double code_ky = ky_freq * 2.0 * M_PI / config.DOMAIN_SIZE;
                double code_kz = kz_freq * 2.0 * M_PI / config.DOMAIN_SIZE;
                double code_k2 =
                    code_kx * code_kx + code_ky * code_ky + code_kz * code_kz;

                std::complex<double> phi_k = -delta_k / code_k2;

                disp_x_k[idx] = std::complex<double>(0, -1) * code_kx * phi_k;
                disp_y_k[idx] = std::complex<double>(0, -1) * code_ky * phi_k;
                disp_z_k[idx] = std::complex<double>(0, -1) * code_kz * phi_k;
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

    double norm_ic = scale_factor / static_cast<double>(M3_real);

    for (size_t i = 0; i < M3_real; ++i) {
        field.dx[i] *= norm_ic;
        field.dy[i] *= norm_ic;
        field.dz[i] *= norm_ic;
    }

    // Calculate the growth rate 'f' for LambdaCDM velocities
    double a3 = pow(scale_factor, 3.0);
    // Omega_m(a) = the matter density at the current scale factor
    double Om_a = config.OMEGA_M / (config.OMEGA_M + config.OMEGA_LAMBDA * a3);
    field.f = pow(Om_a, 0.55);  // The Peebles approximation

    return field;
}

void initialize_dm(SimState& state, const Config& config,
                   const ZeldovichField& zf) {
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

                double dx = zf.dx[idx];
                double dy = zf.dy[idx];
                double dz = zf.dz[idx];

                double p_x =
                    fmod(qx + dx + config.DOMAIN_SIZE, config.DOMAIN_SIZE);
                double p_y =
                    fmod(qy + dy + config.DOMAIN_SIZE, config.DOMAIN_SIZE);
                double p_z =
                    fmod(qz + dz + config.DOMAIN_SIZE, config.DOMAIN_SIZE);

                double v_x = config.STANDING_PARTICLES
                                 ? 0.0
                                 : state.hubble_param * dx * zf.f;
                double v_y = config.STANDING_PARTICLES
                                 ? 0.0
                                 : state.hubble_param * dy * zf.f;
                double v_z = config.STANDING_PARTICLES
                                 ? 0.0
                                 : state.hubble_param * dz * zf.f;
                double mass = config.DM_PARTICLE_MASS;

                state.dm.add_particle(p_x, p_y, p_z, v_x, v_y, v_z,
                                      config.DM_PARTICLE_MASS);
            }
        }
    }
}

double get_internal_energy_from_temp_k(double T_kelvin, double hubble_param,
                                       double box_size_mpc, double domain_size,
                                       double gamma) {
    // Physical internal energy (Joules/kg)
    const double k_B = 1.380649e-23;
    const double m_p = 1.67262192e-27;
    const double mu = 1.22;
    const double u_physical_m2_s2 =
        (k_B * T_kelvin) / ((gamma - 1.0) * mu * m_p);

    // Convert dimensionless 'h' (e.g., 0.7) to absolute H0 (70.0 km/s/Mpc)
    double H0 = hubble_param * 100.0;

    // Formula: V_unit = 1.5 * H0 * (Physical size of 1 code unit)
    const double code_unit_mpc = box_size_mpc / domain_size;
    const double V_unit_km_s = 1.5 * H0 * code_unit_mpc;
    const double V_unit_m_s = V_unit_km_s * 1000.0;  // Convert km/s to m/s

    // Convert to code units
    return u_physical_m2_s2 / (V_unit_m_s * V_unit_m_s);
}

void initialize_gas(SimState& state, const Config& config,
                    const ZeldovichField& zf) {
    if (!config.USE_HYDRO) return;

    int M = config.MESH_SIZE;
    double cell_size = config.DOMAIN_SIZE / M;
    size_t M3_real = static_cast<size_t>(M) * M * M;

    assert(state.dm.num_particles > 0);
    double total_dm_mass = state.dm.num_particles * config.DM_PARTICLE_MASS;
    double mass_ratio = config.GAS_TOTAL_MASS / total_dm_mass;

    auto& gas = state.gas;

    gas.density.data = state.dm.get_rho().data * mass_ratio;
    gas.density.data = (gas.get_density().array() < 1e-12)
                           .select(1e-12, gas.get_density().data);

    const double initial_internal_energy = get_internal_energy_from_temp_k(
        config.INITIAL_GAS_TEMPERATURE_K, config.HUBBLE_PARAM,
        config.BOX_SIZE_MPC, config.DOMAIN_SIZE, config.GAMMA);

    for (size_t i = 0; i < M3_real; ++i) {
        double dx = zf.dx[i];
        double dy = zf.dy[i];
        double dz = zf.dz[i];

        double vx =
            config.STANDING_PARTICLES ? 0.0 : state.hubble_param * dx * zf.f;
        double vy =
            config.STANDING_PARTICLES ? 0.0 : state.hubble_param * dy * zf.f;
        double vz =
            config.STANDING_PARTICLES ? 0.0 : state.hubble_param * dz * zf.f;

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

    ZeldovichField z_field =
        compute_zeldovich_field(state.scale_factor, config);

    // Dark Matter Step
    initialize_dm(state, config, z_field);
    state.dm.bin_and_assign_mass(config);

    // Gas Step
    initialize_gas(state, config, z_field);

    // Forces
    compute_gravitational_acceleration(state, config);

    if (config.USE_PM) {
        // This directly overwrites state.dm.acc_x/y/z
        state.dm.interpolate_cic_forces(state.gravity_x, state.gravity_y,
                                        state.gravity_z, config);
    } else {
        // If PM is off, we must ensure the acceleration arrays start at zero
        std::fill(state.dm.acc_x.begin(), state.dm.acc_x.end(), 0.0);
        std::fill(state.dm.acc_y.begin(), state.dm.acc_y.end(), 0.0);
        std::fill(state.dm.acc_z.begin(), state.dm.acc_z.end(), 0.0);
    }

    if (config.USE_PP) {
        // This adds PP accelerations directly on top of the PM accelerations
        state.dm.compute_pp_forces(config);
    }

    if (!config.STANDING_PARTICLES) {
        double local_max = 1e-9;

        // Grab raw pointers for vectorization
        const size_t n = state.dm.num_particles;
        const double* ax = state.dm.acc_x.data();
        const double* ay = state.dm.acc_y.data();
        const double* az = state.dm.acc_z.data();

        // Find the maximum acceleration
        for (size_t i = 0; i < n; ++i) {
            double accel_sq = ax[i] * ax[i] + ay[i] * ay[i] + az[i] * az[i];
            local_max = (accel_sq > local_max) ? accel_sq : local_max;
        }

        state.dm.max_accel_sq = local_max;
    } else {
        // If particles are frozen, wipe the accelerations clean
        std::fill(state.dm.acc_x.begin(), state.dm.acc_x.end(), 0.0);
        std::fill(state.dm.acc_y.begin(), state.dm.acc_y.end(), 0.0);
        std::fill(state.dm.acc_z.begin(), state.dm.acc_z.end(), 0.0);
        state.dm.max_accel_sq = 1e-9;
    }

    return state;
}
