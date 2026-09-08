#include "ics.h"

#include <cassert>
#include <random>

#include "constants.h"
#include "cooling.h"
#include "integrator.h"
#include "math_utils.h"
#include "particles.h"
#include "pocketfft_hdronly.h"

// The pure shape of the BBKS power spectrum (Unnormalized)
static double unnormalized_pk(double k_h, const Config& config) {
    // k_h is the wavenumber in units of h Mpc^-1
    // Cosmological Shape Parameter (Gamma)
    double Gamma = config.omega_m * config.hubble_h;

    // 'q' parameter for BBKS, scaled for Mpc/h units
    double q = k_h / Gamma;

    // BBKS Polynomial
    double T_k = log(1.0 + 2.34 * q) / (2.34 * q) *
                 pow(1.0 + 3.89 * q + pow(16.1 * q, 2) + pow(5.46 * q, 3) +
                         pow(6.71 * q, 4),
                     -0.25);

    // P(k) = k^n_s * T(k)^2
    return pow(k_h, config.spectral_index) * T_k * T_k;
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
    return (config.sigma_8 * config.sigma_8) / unnorm_variance;
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

    int M = config.mesh_size;
    size_t M3_real = static_cast<size_t>(M) * M * M;
    size_t M3_complex = static_cast<size_t>(M) * M * (M / 2 + 1);

    // Convert box size to h^-1 Mpc to match the units of true_pk
    double box_h_mpc = config.box_size_mpc * config.hubble_h;
    double box_vol_h = pow(box_h_mpc, 3.0);

    double amplitude_scaling =
        std::sqrt(static_cast<double>(M3_real) / box_vol_h);

    // Create Gaussian random field
    std::vector<double> real_space_random_field(M3_real);
    std::default_random_engine generator(config.seed);
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
                double kx = kx_freq * 2.0 * M_PI / config.box_size_mpc;
                double ky = ky_freq * 2.0 * M_PI / config.box_size_mpc;
                double kz = kz_freq * 2.0 * M_PI / config.box_size_mpc;
                double k_mag = sqrt(kx * kx + ky * ky + kz * kz);
                double k_h = k_mag / config.hubble_h;

                // The BBKS Transfer Function T(k)
                // Multiply the unnormalized shape by A
                double true_pk = A * unnormalized_pk(k_h, config);

                // The power spectrum P(k) represents the variance of the
                // density field. To scale our Gaussian white noise (which has a
                // variance of 1.0), we need the amplitude (standard deviation),
                // hence the square root. We also apply the grid-to-physical
                // volume scaling to bridge the discrete/continuous gap.
                double power_spectrum_sqrt = sqrt(true_pk) * amplitude_scaling;

                std::complex<double> current_k = random_k[idx];

                // Invert the phase by multiplying the complex number by -1.
                // This turns peaks into voids and voids into peaks.
                if (config.invert_phases) {
                    current_k = -current_k;
                }

                // Apply the scaling to the Fourier amplitudes
                std::complex<double> delta_k;
                if (config.fixed_ics) {
                    // Extract the random phase, discard the random magnitude.
                    // std::abs() on a complex number returns its magnitude.
                    double mag = std::abs(current_k);

                    // Protect against division by zero
                    std::complex<double> phase_only =
                        (mag > 1e-15) ? (current_k / mag)
                                      : std::complex<double>(0, 0);

                    // Restore the expected baseline FFT magnitude (sqrt of
                    // total grid cells)
                    double expected_mag =
                        std::sqrt(static_cast<double>(M3_real));

                    // Multiply the phase (magnitude = 1) by the theoretical
                    // amplitude
                    delta_k = phase_only * expected_mag * power_spectrum_sqrt;
                } else {
                    delta_k = current_k * power_spectrum_sqrt;
                }

                // Convert density contrast to Zel'dovich displacement
                // potential
                double code_kx = kx_freq * 2.0 * M_PI / config.domain_size;
                double code_ky = ky_freq * 2.0 * M_PI / config.domain_size;
                double code_kz = kz_freq * 2.0 * M_PI / config.domain_size;
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
    double Om = config.omega_m;
    double Ol = config.omega_lambda;
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
    double Om_a = config.omega_m / (config.omega_m + config.omega_lambda * a3);
    field.f = pow(Om_a, 0.55);  // The Peebles approximation

    return field;
}

// Simple struct to hold the 3D displacement vectors
struct Vec3 {
    double x, y, z;
};

// Evaluates the continuous Zel'dovich field at arbitrary coordinates
// using Cloud-In-Cell (CIC) interpolation.
static Vec3 sample_zeldovich_cic(double qx, double qy, double qz,
                                 const ZeldovichField& zf, double cell_size,
                                 int M) {
    double x_idx = qx / cell_size;
    double y_idx = qy / cell_size;
    double z_idx = qz / cell_size;

    int ix = static_cast<int>(x_idx);
    int iy = static_cast<int>(y_idx);
    int iz = static_cast<int>(z_idx);

    double frac_x = x_idx - ix;
    double frac_y = y_idx - iy;
    double frac_z = z_idx - iz;

    double w000 = (1.0 - frac_x) * (1.0 - frac_y) * (1.0 - frac_z);
    double w100 = frac_x * (1.0 - frac_y) * (1.0 - frac_z);
    double w010 = (1.0 - frac_x) * frac_y * (1.0 - frac_z);
    double w110 = frac_x * frac_y * (1.0 - frac_z);
    double w001 = (1.0 - frac_x) * (1.0 - frac_y) * frac_z;
    double w101 = frac_x * (1.0 - frac_y) * frac_z;
    double w011 = (1.0 - frac_x) * frac_y * frac_z;
    double w111 = frac_x * frac_y * frac_z;

    int ix0 = (ix + M) % M, ix1 = (ix + 1 + M) % M;
    int iy0 = (iy + M) % M, iy1 = (iy + 1 + M) % M;
    int iz0 = (iz + M) % M, iz1 = (iz + 1 + M) % M;

    auto interp = [&](const std::vector<double>& field) {
        return field[ix0 * M * M + iy0 * M + iz0] * w000 +
               field[ix1 * M * M + iy0 * M + iz0] * w100 +
               field[ix0 * M * M + iy1 * M + iz0] * w010 +
               field[ix1 * M * M + iy1 * M + iz0] * w110 +
               field[ix0 * M * M + iy0 * M + iz1] * w001 +
               field[ix1 * M * M + iy0 * M + iz1] * w101 +
               field[ix0 * M * M + iy1 * M + iz1] * w011 +
               field[ix1 * M * M + iy1 * M + iz1] * w111;
    };

    return {interp(zf.dx), interp(zf.dy), interp(zf.dz)};
}

void initialize_dm(SimState& state, const Config& config,
                   const ZeldovichField& zf) {
    int M = config.mesh_size;
    double cell_size = config.domain_size / M;
    int N_part = config.num_particles_1d;
    double spacing = config.domain_size / N_part;

    for (int i = 0; i < N_part; ++i) {
        for (int j = 0; j < N_part; ++j) {
            for (int k = 0; k < N_part; ++k) {
                double qx = (i + 0.5) * spacing;
                double qy = (j + 0.5) * spacing;
                double qz = (k + 0.5) * spacing;

                Vec3 d = sample_zeldovich_cic(qx, qy, qz, zf, cell_size, M);

                double p_x =
                    fmod(qx + d.x + config.domain_size, config.domain_size);
                double p_y =
                    fmod(qy + d.y + config.domain_size, config.domain_size);
                double p_z =
                    fmod(qz + d.z + config.domain_size, config.domain_size);

                double v_x = config.standing_particles
                                 ? 0.0
                                 : state.hubble_param * d.x * zf.f;
                double v_y = config.standing_particles
                                 ? 0.0
                                 : state.hubble_param * d.y * zf.f;
                double v_z = config.standing_particles
                                 ? 0.0
                                 : state.hubble_param * d.z * zf.f;

                state.dm.add_particle(p_x, p_y, p_z, v_x, v_y, v_z,
                                      config.dm_particle_mass);
            }
        }
    }
}

void initialize_gas(SimState& state, const Config& config,
                    const ZeldovichField& zf) {
    if (config.hydro_method == HydroMethod::None) return;

    int M = config.mesh_size;
    double seed_metallicity =
        config.seed_metallicity_solar * constants::Z_SOLAR;

    const double initial_internal_energy =
        Cooling::get_internal_energy_from_temp(config.initial_gas_temperature_k,
                                               state.scale_factor, config);

    if (config.hydro_method == HydroMethod::Eulerian) {
        auto& gas = *state.gas;
        double mean_gas_rho =
            config.gas_total_mass / std::pow(config.domain_size, 3.0);
        double inv_2dx = 1.0 / (2.0 * config.cell_size);

        for (int i = 0; i < M; ++i) {
            for (int j = 0; j < M; ++j) {
                for (int k = 0; k < M; ++k) {
                    size_t idx = static_cast<size_t>(i) * M * M + j * M + k;

                    // Calculate cell center
                    double qx = (i + 0.5) * config.cell_size;
                    double qy = (j + 0.5) * config.cell_size;
                    double qz = (k + 0.5) * config.cell_size;

                    // Sample local Zel'dovich displacement
                    Vec3 d = sample_zeldovich_cic(qx, qy, qz, zf,
                                                  config.cell_size, M);

                    // Central difference for divergence (sampled at
                    // neighboring cell centers)
                    Vec3 dx_ip1 = sample_zeldovich_cic(
                        qx + config.cell_size, qy, qz, zf, config.cell_size, M);
                    Vec3 dx_im1 = sample_zeldovich_cic(
                        qx - config.cell_size, qy, qz, zf, config.cell_size, M);
                    Vec3 dy_jp1 = sample_zeldovich_cic(
                        qx, qy + config.cell_size, qz, zf, config.cell_size, M);
                    Vec3 dy_jm1 = sample_zeldovich_cic(
                        qx, qy - config.cell_size, qz, zf, config.cell_size, M);
                    Vec3 dz_kp1 = sample_zeldovich_cic(
                        qx, qy, qz + config.cell_size, zf, config.cell_size, M);
                    Vec3 dz_km1 = sample_zeldovich_cic(
                        qx, qy, qz - config.cell_size, zf, config.cell_size, M);

                    double div = (dx_ip1.x - dx_im1.x) * inv_2dx +
                                 (dy_jp1.y - dy_jm1.y) * inv_2dx +
                                 (dz_kp1.z - dz_km1.z) * inv_2dx;

                    // Calculate Density
                    double rho = mean_gas_rho * (1.0 - div);
                    double safe_rho = std::max(rho, 1e-12);
                    gas.density.data[idx] = safe_rho;
                    gas.metal_density.data[idx] = safe_rho * seed_metallicity;

                    // Calculate Velocity and Energy
                    double vx = config.standing_particles
                                    ? 0.0
                                    : state.hubble_param * d.x * zf.f;
                    double vy = config.standing_particles
                                    ? 0.0
                                    : state.hubble_param * d.y * zf.f;
                    double vz = config.standing_particles
                                    ? 0.0
                                    : state.hubble_param * d.z * zf.f;

                    gas.velocity_x.data[idx] = vx;
                    gas.velocity_y.data[idx] = vy;
                    gas.velocity_z.data[idx] = vz;
                    gas.momentum_x.data[idx] = safe_rho * vx;
                    gas.momentum_y.data[idx] = safe_rho * vy;
                    gas.momentum_z.data[idx] = safe_rho * vz;

                    double kin_energy =
                        0.5 * safe_rho * (vx * vx + vy * vy + vz * vz);
                    gas.energy.data[idx] =
                        (safe_rho * initial_internal_energy) + kin_energy;
                    gas.internal_energy.data[idx] =
                        safe_rho * initial_internal_energy;
                }
            }
        }
        gas.update_primitive_variables(state.scale_factor);
    } else if (config.hydro_method == HydroMethod::MFM) {
        int N_part = config.num_gas_particles_1d;
        double spacing = config.domain_size / N_part;
        double gas_particle_mass =
            config.gas_total_mass / (N_part * N_part * N_part);

        // Initial guess for smoothing length
        double initial_h = 1.2 * spacing;

        for (int i = 0; i < N_part; ++i) {
            for (int j = 0; j < N_part; ++j) {
                for (int k = 0; k < N_part; ++k) {
                    // Shift the gas lattice by half a spacing relative to DM to
                    // avoid perfect overlap
                    double qx = i * spacing;
                    double qy = j * spacing;
                    double qz = k * spacing;

                    Vec3 d = sample_zeldovich_cic(qx, qy, qz, zf,
                                                  config.cell_size, M);

                    double p_x =
                        fmod(qx + d.x + config.domain_size, config.domain_size);
                    double p_y =
                        fmod(qy + d.y + config.domain_size, config.domain_size);
                    double p_z =
                        fmod(qz + d.z + config.domain_size, config.domain_size);

                    double v_x = state.hubble_param * d.x * zf.f;
                    double v_y = state.hubble_param * d.y * zf.f;
                    double v_z = state.hubble_param * d.z * zf.f;

                    state.mfm_gas->add_particle(
                        p_x, p_y, p_z, v_x, v_y, v_z, gas_particle_mass,
                        initial_internal_energy, initial_h, seed_metallicity);
                }
            }
        }
    }
}

void initialize_sod_shock_tube(SimState& state, const Config& config) {
    // GADGET/GIZMO Shock Tube parameters (Hernquist & Katz 1989)
    double rho_L = 1.0;
    double rho_R = 0.25;
    double P_L = 1.0;
    double P_R = 0.1795;

    double L = config.domain_size;
    double gamma = config.gamma;
    double u_L = P_L / (rho_L * (gamma - 1.0));
    double u_R = P_R / (rho_R * (gamma - 1.0));
    double v_x = 0.0, v_y = 0.0, v_z = 0.0;
    double seed_metallicity = 0.0;

    if (config.hydro_method == HydroMethod::Eulerian) {
        auto& gas = *state.gas;
        int M = config.mesh_size;

        for (int i = 0; i < M; ++i) {
            for (int j = 0; j < M; ++j) {
                for (int k = 0; k < M; ++k) {
                    size_t idx = static_cast<size_t>(i) * M * M + j * M + k;

                    // Calculate cell center x-coordinate
                    double x = (i + 0.5) * config.cell_size;

                    // Assign Left or Right state based on the midpoint
                    double rho = (x < L / 2.0) ? rho_L : rho_R;
                    double u = (x < L / 2.0) ? u_L : u_R;

                    gas.density.data[idx] = rho;
                    gas.metal_density.data[idx] = rho * seed_metallicity;

                    gas.velocity_x.data[idx] = v_x;
                    gas.velocity_y.data[idx] = v_y;
                    gas.velocity_z.data[idx] = v_z;
                    gas.momentum_x.data[idx] = rho * v_x;
                    gas.momentum_y.data[idx] = rho * v_y;
                    gas.momentum_z.data[idx] = rho * v_z;

                    double kin_energy =
                        0.5 * rho * (v_x * v_x + v_y * v_y + v_z * v_z);
                    gas.energy.data[idx] = (rho * u) + kin_energy;
                    gas.internal_energy.data[idx] = rho * u;
                }
            }
        }
        gas.update_primitive_variables(state.scale_factor);

    } else if (config.hydro_method == HydroMethod::MFM) {
        // Using the 1D particle count to define the left-side resolution
        int N = config.num_gas_particles_1d;

        // Spacing: Right side spacing must scale by the cube root of the
        // density ratio in 3D to maintain equal particle masses across the
        // domain.
        double dx_L = L / N;
        double dx_R =
            dx_L *
            std::cbrt(rho_L / rho_R);  // For a 4:1 ratio, cbrt(4) ≈ 1.5874

        // Since mass is constant, we define it based on the left side
        double particle_mass = rho_L * (dx_L * dx_L * dx_L);

        // Calculate starting offsets to ensure L/2.0 is exactly hit by the grid
        std::vector<double> yz_coords_L;
        double start_L = std::fmod(L / 2.0, dx_L);
        if (start_L < 0.0) start_L += dx_L;  // Prevent negative float quirk
        for (double pos = start_L; pos < L - 1e-5; pos += dx_L) {
            yz_coords_L.push_back(pos);
        }

        std::vector<double> yz_coords_R;
        double start_R = std::fmod(L / 2.0, dx_R);
        if (start_R < 0.0) start_R += dx_R;
        for (double pos = start_R; pos < L - 1e-5; pos += dx_R) {
            yz_coords_R.push_back(pos);
        }

        // LEFT SIDE (High Density, High Pressure)
        // x goes from 0 to 0.5 * L (X stays unchanged, Y/Z use centered arrays)
        for (double x = dx_L / 2.0; x < 0.5 * L; x += dx_L) {
            for (double y : yz_coords_L) {
                for (double z : yz_coords_L) {
                    // Keep the smoothing length slightly larger than the
                    // spacing to ensure overlap
                    double initial_h = 1.2 * dx_L;
                    state.mfm_gas->add_particle(x, y, z, v_x, v_y, v_z,
                                                particle_mass, u_L, initial_h,
                                                seed_metallicity);
                }
            }
        }

        // RIGHT SIDE (Low Density, Low Pressure)
        // x goes from 0.5 * L to L
        for (double x = 0.5 * L + dx_R / 2.0; x < L; x += dx_R) {
            for (double y : yz_coords_R) {
                for (double z : yz_coords_R) {
                    double initial_h = 1.2 * dx_R;
                    state.mfm_gas->add_particle(x, y, z, v_x, v_y, v_z,
                                                particle_mass, u_R, initial_h,
                                                seed_metallicity);
                }
            }
        }
    }
}

void initialize_adiabatic_expansion(SimState& state, const Config& config) {
    const double initial_internal_energy =
        Cooling::get_internal_energy_from_temp(config.initial_gas_temperature_k,
                                               state.scale_factor, config);

    double seed_metallicity = 0.0;
    double v_x = 0.0, v_y = 0.0, v_z = 0.0;

    if (config.hydro_method == HydroMethod::Eulerian) {
        auto& gas = *state.gas;
        int M = config.mesh_size;

        // Calculate the uniform background density
        double mean_gas_rho =
            config.gas_total_mass / std::pow(config.domain_size, 3.0);

        for (int i = 0; i < M; ++i) {
            for (int j = 0; j < M; ++j) {
                for (int k = 0; k < M; ++k) {
                    size_t idx = static_cast<size_t>(i) * M * M + j * M + k;

                    // Uniform density and metallicity
                    gas.density.data[idx] = mean_gas_rho;
                    gas.metal_density.data[idx] =
                        mean_gas_rho * seed_metallicity;

                    // Zero velocity and momentum
                    gas.velocity_x.data[idx] = v_x;
                    gas.velocity_y.data[idx] = v_y;
                    gas.velocity_z.data[idx] = v_z;
                    gas.momentum_x.data[idx] = 0.0;
                    gas.momentum_y.data[idx] = 0.0;
                    gas.momentum_z.data[idx] = 0.0;

                    // Because kinetic energy is zero, total energy = internal
                    // energy
                    gas.energy.data[idx] =
                        mean_gas_rho * initial_internal_energy;
                    gas.internal_energy.data[idx] =
                        mean_gas_rho * initial_internal_energy;
                }
            }
        }
        gas.update_primitive_variables(state.scale_factor);

    } else if (config.hydro_method == HydroMethod::MFM) {
        int N_part = config.num_gas_particles_1d;
        double spacing = config.domain_size / N_part;
        double gas_particle_mass =
            config.gas_total_mass / (N_part * N_part * N_part);

        // Initial guess for smoothing length
        double initial_h = 1.2 * spacing;

        for (int i = 0; i < N_part; ++i) {
            for (int j = 0; j < N_part; ++j) {
                for (int k = 0; k < N_part; ++k) {
                    // Shift the gas lattice by half a spacing relative to DM to
                    // avoid perfect overlap
                    double qx = i * spacing;
                    double qy = j * spacing;
                    double qz = k * spacing;

                    double p_x =
                        fmod(qx + config.domain_size, config.domain_size);
                    double p_y =
                        fmod(qy + config.domain_size, config.domain_size);
                    double p_z =
                        fmod(qz + config.domain_size, config.domain_size);

                    state.mfm_gas->add_particle(
                        p_x, p_y, p_z, v_x, v_y, v_z, gas_particle_mass,
                        initial_internal_energy, initial_h, seed_metallicity);
                }
            }
        }
    }
}

void initialize_sedov_blastwave(SimState& state, const Config& config) {
    // Standard Sedov-Taylor parameters
    double E_total = 1.0;
    double rho_bg = 1.0;
    double P_bg = 1e-5;  // Almost zero pressure for the background to create a
                         // strong shock
    double gamma = config.gamma;
    double u_bg = P_bg / (rho_bg * (gamma - 1.0));
    double seed_metallicity = 0.0;
    double L = config.domain_size;
    double center = L / 2.0;

    // Define the injection radius
    int N_res = (config.hydro_method == HydroMethod::Eulerian)
                    ? config.mesh_size
                    : config.num_gas_particles_1d;
    double dx = L / static_cast<double>(N_res);
    double R_inj = 2.0 * dx;

    if (config.hydro_method == HydroMethod::Eulerian) {
        auto& gas = *state.gas;
        int M = config.mesh_size;

        // Initialize uniform cold background
        for (int i = 0; i < M; ++i) {
            for (int j = 0; j < M; ++j) {
                for (int k = 0; k < M; ++k) {
                    size_t idx = static_cast<size_t>(i) * M * M + j * M + k;
                    gas.density.data[idx] = rho_bg;
                    gas.metal_density.data[idx] = 0.0;
                    gas.velocity_x.data[idx] = 0.0;
                    gas.velocity_y.data[idx] = 0.0;
                    gas.velocity_z.data[idx] = 0.0;
                    gas.momentum_x.data[idx] = 0.0;
                    gas.momentum_y.data[idx] = 0.0;
                    gas.momentum_z.data[idx] = 0.0;
                    gas.internal_energy.data[idx] = rho_bg * u_bg;
                    gas.energy.data[idx] = rho_bg * u_bg;
                }
            }
        }

        // Inject energy into the central region
        double cell_vol = dx * dx * dx;
        double inj_mass = 0.0;
        std::vector<size_t> inj_indices;

        for (int i = 0; i < M; ++i) {
            for (int j = 0; j < M; ++j) {
                for (int k = 0; k < M; ++k) {
                    double x = (i + 0.5) * dx;
                    double y = (j + 0.5) * dx;
                    double z = (k + 0.5) * dx;
                    double r = std::sqrt((x - center) * (x - center) +
                                         (y - center) * (y - center) +
                                         (z - center) * (z - center));

                    if (r <= R_inj) {
                        inj_indices.push_back(static_cast<size_t>(i) * M * M +
                                              j * M + k);
                        inj_mass += rho_bg * cell_vol;
                    }
                }
            }
        }

        // Distribute specific internal energy
        double u_inj = E_total / inj_mass;
        for (size_t idx : inj_indices) {
            gas.internal_energy.data[idx] += rho_bg * u_inj;
            gas.energy.data[idx] += rho_bg * u_inj;
        }
        gas.update_primitive_variables(state.scale_factor);

    } else if (config.hydro_method == HydroMethod::MFM) {
        // Use a Face-Centered Cubic (FCC) lattice.
        // FCC is close-packed (12 equidistant nearest neighbors) yielding good
        // isotropic shock propagation, and tiles inside a cubic box.

        // Calculate how many FCC unit cells we need to approximate the
        // requested resolution
        int N_cell = static_cast<int>(std::round(N_res / std::cbrt(4.0)));
        if (N_cell < 1) N_cell = 1;

        int total_particles = 4 * N_cell * N_cell * N_cell;
        double L_c =
            L / static_cast<double>(N_cell);  // Size of the cubic unit cell

        double gas_particle_mass =
            rho_bg * (L * L * L) / static_cast<double>(total_particles);
        double effective_dx = std::cbrt(gas_particle_mass / rho_bg);
        double initial_h = 1.2 * effective_dx;

        struct PartData {
            double x, y, z, r, weight;
        };
        std::vector<PartData> pdata;
        pdata.reserve(total_particles);

        double w_sum = 0.0;

        // FCC Basis vectors inside a unit cell (scaled by L_c)
        double basis[4][3] = {
            {0.0, 0.0, 0.0}, {0.5, 0.5, 0.0}, {0.5, 0.0, 0.5}, {0.0, 0.5, 0.5}};

        // Lay out FCC particles and calculate energy weights
        for (int i = 0; i < N_cell; ++i) {
            for (int j = 0; j < N_cell; ++j) {
                for (int k = 0; k < N_cell; ++k) {
                    for (int b = 0; b < 4; ++b) {
                        double qx = (i + basis[b][0]) * L_c;
                        double qy = (j + basis[b][1]) * L_c;
                        double qz = (k + basis[b][2]) * L_c;

                        // Microscopic grid-lock breaking noise
                        /*double noise_x = ((rand() / (double)RAND_MAX) - 0.5) *
                                         1e-4 * effective_dx;
                        double noise_y = ((rand() / (double)RAND_MAX) - 0.5) *
                                         1e-4 * effective_dx;
                        double noise_z = ((rand() / (double)RAND_MAX) - 0.5) *
                                         1e-4 * effective_dx;*/

                        double p_x = std::fmod(qx + /*noise_x +*/ L, L);
                        double p_y = std::fmod(qy + /*noise_y +*/ L, L);
                        double p_z = std::fmod(qz + /*noise_z +*/ L, L);

                        double r = std::sqrt((p_x - center) * (p_x - center) +
                                             (p_y - center) * (p_y - center) +
                                             (p_z - center) * (p_z - center));

                        // Gaussian energy smoothing
                        double w = 0.0;
                        if (r < 3.0 * R_inj) {
                            w = std::exp(-(r * r) / (R_inj * R_inj));
                        }
                        w_sum += w;

                        pdata.push_back({p_x, p_y, p_z, r, w});
                    }
                }
            }
        }

        // Inject energy and initialize particles
        constexpr bool smooth_energy_injection = false;

        // Find the central particle for single-particle injection
        double min_r = std::numeric_limits<double>::max();
        size_t central_idx = 0;
        if (!smooth_energy_injection) {
            for (size_t i = 0; i < pdata.size(); ++i) {
                if (pdata[i].r < min_r) {
                    min_r = pdata[i].r;
                    central_idx = i;
                }
            }
        }

        // Inject energy and initialize particles
        for (size_t i = 0; i < pdata.size(); ++i) {
            const auto& pd = pdata[i];
            double particle_u = u_bg;

            if (smooth_energy_injection) {
                // Distribute E_total based on the gaussian weights
                if (pd.weight > 0.0 && w_sum > 0.0) {
                    double injected_E = E_total * (pd.weight / w_sum);
                    particle_u += (injected_E / gas_particle_mass);
                }
            } else {
                if (i == central_idx) {
                    particle_u += (E_total / gas_particle_mass);
                }
            }

            state.mfm_gas->add_particle(pd.x, pd.y, pd.z, 0.0, 0.0, 0.0,
                                        gas_particle_mass, particle_u,
                                        initial_h, seed_metallicity);
        }
    }
}

SimState initialize_state(Config& config) {
    SimState state(config);
    state.total_time = 0;
    bool prev_expanding_universe = config.expanding_universe;
    config.expanding_universe = true;
    update_cosmology(state, config);

    if (config.initial_setup == InitialSetup::Cosmological) {
        ZeldovichField z_field =
            compute_zeldovich_field(state.scale_factor, config);

        // Dark Matter Step
        initialize_dm(state, config, z_field);
        state.dm.bin_and_assign_mass(config);

        // Gas Step
        initialize_gas(state, config, z_field);
    } else if (config.initial_setup == InitialSetup::SodShockTube) {
        initialize_sod_shock_tube(state, config);
    } else if (config.initial_setup == InitialSetup::AdiabaticExpansion) {
        initialize_adiabatic_expansion(state, config);
    } else if (config.initial_setup == InitialSetup::SedovBlastwave) {
        initialize_sedov_blastwave(state, config);
    }

    config.expanding_universe = prev_expanding_universe;
    update_cosmology(state, config);

    // Compute initial MFM properties before forces are calculated
    if (config.hydro_method == HydroMethod::MFM) {
        state.mfm_gas->compute_density_and_h(config, state.dm);
        state.mfm_gas->bin_and_assign_mass(config);
        state.mfm_gas->update_primitive_variables(config, state.scale_factor);
    }

    Diagnostics dummy_diag;
    compute_forces(state, config, dummy_diag);

    if (config.hydro_method == HydroMethod::MFM) {
        // Sync pressure and energies based on the initial density
        state.mfm_gas->update_primitive_variables(config, state.scale_factor);
        state.mfm_gas->compute_gradients(config);
        state.mfm_gas->compute_hydro_forces(config, state.scale_factor, 0.0);
    }

    return state;
}
