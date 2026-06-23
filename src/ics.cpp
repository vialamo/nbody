#include "ics.h"

#include <cassert>
#include <random>

#include "cooling.h"
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

static double compute_normalization_constant(const Config& config) {
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

    // The normalization constant A = (sigma_8^2) / (unnormalized variance)
    return (config.SIGMA_8 * config.SIGMA_8) / unnorm_variance;
}

/* Computes the initial Zel'dovich displacement field and velocity growth rate.
 This function generates a Gaussian random field shaped by the LambdaCDM
 power spectrum (using the BBKS transfer function) and normalized to the
 present-day sigma_8 variance. It solves Poisson's equation in Fourier space
 and takes the negative gradient to produce the displacement vectors.
 Returns a ZeldovichField struct containing:
 * The 3D displacement vectors in internal CODE UNITS (scaled to DOMAIN_SIZE).
 They have already been scaled to the initial growth factor (D)
 and are ready to be directly added to the unperturbed particle lattice.
 * The logarithmic growth rate (d ln D / d ln a) evaluated at the
 initial epoch, required by the caller to calculate peculiar velocities.
 */
ZeldovichField compute_zeldovich_field(double scale_factor,
                                       const Config& config) {
    // The master normalization constant
    double A = compute_normalization_constant(config);

    int M = config.MESH_SIZE;
    size_t M3_real = static_cast<size_t>(M) * M * M;
    size_t M3_complex = static_cast<size_t>(M) * M * (M / 2 + 1);

    // Convert box size to h^-1 Mpc to match the units of true_pk
    double box_h_mpc = config.BOX_SIZE_MPC * config.HUBBLE_PARAM;
    double box_vol_h = pow(box_h_mpc, 3.0);

    double amplitude_scaling =
        std::sqrt(static_cast<double>(M3_real) / box_vol_h);

    // Create Gaussian random field
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

                size_t idx = static_cast<size_t>(i) * M * (M / 2 + 1) +
                             static_cast<size_t>(j) * (M / 2 + 1) +
                             static_cast<size_t>(k);

                auto kx_freq = static_cast<double>((i < M / 2) ? i : (i - M));
                auto ky_freq = static_cast<double>((j < M / 2) ? j : (j - M));
                auto kz_freq = static_cast<double>(k);
                double kx = kx_freq * 2.0 * M_PI / config.BOX_SIZE_MPC;
                double ky = ky_freq * 2.0 * M_PI / config.BOX_SIZE_MPC;
                double kz = kz_freq * 2.0 * M_PI / config.BOX_SIZE_MPC;
                double k_mag = sqrt(kx * kx + ky * ky + kz * kz);
                double k_h = k_mag / config.HUBBLE_PARAM;

                // The BBKS Transfer Function T(k)
                // Multiply the unnormalized shape by A
                double true_pk = A * unnormalized_pk(k_h, config);

                // The power spectrum P(k) represents the variance of the
                // density field. To scale our Gaussian white noise (which has a
                // variance of 1.0), we need the amplitude (standard deviation),
                // hence the square root. We also apply the grid-to-physical
                // volume scaling to bridge the discrete/continuous gap.
                double power_spectrum_sqrt = sqrt(true_pk) * amplitude_scaling;

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

                // Calculate the Zel'dovich displacement by taking the negative
                // gradient of the potential (-nabla Phi). In Fourier space,
                // taking the negative gradient translates to multiplying by -i
                // * k. std::complex<double>(0, -1) is the exact C++
                // representation of -i.
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

    // Calculate the exact linear growth suppression factor g(a) = D(a)/a at z=0
    // (a=1) Using the Carroll, Press, and Turner (1992) fitting formula
    double Om = config.OMEGA_M;
    double Ol = config.OMEGA_LAMBDA;
    double g_1 = (2.5 * Om) / (pow(Om, 4.0 / 7.0) - Ol +
                               (1.0 + Om / 2.0) * (1.0 + Ol / 70.0));

    // Rewind the z=0 density field back to the starting epoch.
    // Because structure growth stalls at late times due to Dark Energy, in a
    // LambdaCDM D(1) is only ~0.77, not 1.0. We must divide by g_1 to recover
    // the true early-universe amplitude.
    double true_growth_factor = scale_factor / g_1;

    // Calculate final normalizer (incorporating the 1/M^3 FFT correction)
    double norm_ic = true_growth_factor / static_cast<double>(M3_real);

    // Finalize the real-space displacement vectors. This loop does two things:
    // * Divides by M^3 to correct the unnormalized output of the inverse FFT.
    // * Multiplies by the initial growth factor (a) to wind the present-day
    //   displacements back to the starting epoch
    for (size_t i = 0; i < M3_real; ++i) {
        field.dx[i] *= norm_ic;
        field.dy[i] *= norm_ic;
        field.dz[i] *= norm_ic;
    }

    // Calculate the growth rate 'f' (d ln D / d ln a)
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

void initialize_gas(SimState& state, const Config& config,
                    const ZeldovichField& zf) {
    if (!config.USE_HYDRO) return;

    int M = config.MESH_SIZE;
    size_t M3_real = static_cast<size_t>(M) * M * M;

    assert(state.dm.num_particles > 0);
    double total_dm_mass = state.dm.num_particles * config.DM_PARTICLE_MASS;
    double mass_ratio = config.GAS_TOTAL_MASS / total_dm_mass;

    auto& gas = state.gas;

    gas.density.data = state.dm.get_rho().data * mass_ratio;
    gas.density.data = (gas.get_density().array() < 1e-12)
                           .select(1e-12, gas.get_density().data);

    const double initial_internal_energy =
        cooling::get_internal_energy_from_temp(config.INITIAL_GAS_TEMPERATURE_K,
                                               state.scale_factor, config);

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
        gas.internal_energy.data[i] = rho * initial_internal_energy;
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

    Diagnostics dummy_diag;
    compute_forces(state, config, dummy_diag);

    return state;
}
