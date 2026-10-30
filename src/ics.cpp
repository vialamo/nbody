#include "ics.h"

#include <fftw3.h>

#include <random>

#include "constants.h"
#include "cooling.h"
#include "cosmology.h"
#include "diagnostics.h"
#include "hdf5_reader.h"
#include "integrator.h"
#include "math_utils.h"
#include "particles.h"
#include "thermodynamics.h"

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

    int Mx = config.mesh_size.x;
    int My = config.mesh_size.y;
    int Mz = config.mesh_size.z;

    size_t M3_real = static_cast<size_t>(Mx) * My * Mz;
    size_t M3_complex = static_cast<size_t>(Mx) * My * (Mz / 2 + 1);

    // Anisotropic volumes
    double box_vol_h = (config.box_size_mpc.x * config.hubble_h) *
                       (config.box_size_mpc.y * config.hubble_h) *
                       (config.box_size_mpc.z * config.hubble_h);

    double amplitude_scaling =
        std::sqrt(static_cast<double>(M3_real) / box_vol_h);

    // Create Gaussian random field
    std::vector<double> real_space_random_field(M3_real);
    std::default_random_engine generator(config.seed);
    std::normal_distribution<double> distribution(0.0, 1.0);
    for (auto& val : real_space_random_field) {
        val = distribution(generator);
    }

    std::vector<std::complex<double>> random_k(M3_complex);

    // Initialize FFTW threads
    fftw_init_threads();
    fftw_plan_with_nthreads(omp_get_max_threads());

    fftw_complex* random_k_fftw =
        reinterpret_cast<fftw_complex*>(random_k.data());

    // Create and execute forward plan
    fftw_plan forward_plan =
        fftw_plan_dft_r2c_3d(Mx, My, Mz, real_space_random_field.data(),
                             random_k_fftw, FFTW_ESTIMATE);

    fftw_execute(forward_plan);
    fftw_destroy_plan(forward_plan);

    std::vector<std::complex<double>> disp_x_k(M3_complex);
    std::vector<std::complex<double>> disp_y_k(M3_complex);
    std::vector<std::complex<double>> disp_z_k(M3_complex);

    for (int i = 0; i < Mx; ++i) {
        for (int j = 0; j < My; ++j) {
            for (int k = 0; k < Mz / 2 + 1; ++k) {
                if (i == 0 && j == 0 && k == 0) {
                    disp_x_k[0] = {0, 0};
                    disp_y_k[0] = {0, 0};
                    disp_z_k[0] = {0, 0};
                    continue;
                }

                size_t idx = static_cast<size_t>(i) * My * (Mz / 2 + 1) +
                             static_cast<size_t>(j) * (Mz / 2 + 1) +
                             static_cast<size_t>(k);

                auto kx_freq = static_cast<double>((i < Mx / 2) ? i : (i - Mx));
                auto ky_freq = static_cast<double>((j < My / 2) ? j : (j - My));
                auto kz_freq = static_cast<double>(k);
                double kx = kx_freq * 2.0 * M_PI / config.box_size_mpc.x;
                double ky = ky_freq * 2.0 * M_PI / config.box_size_mpc.y;
                double kz = kz_freq * 2.0 * M_PI / config.box_size_mpc.z;
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
                    // std::abs() on a complex number returns its magnitude
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
                double code_kx = kx_freq * 2.0 * M_PI / config.domain_size.x;
                double code_ky = ky_freq * 2.0 * M_PI / config.domain_size.y;
                double code_kz = kz_freq * 2.0 * M_PI / config.domain_size.z;
                double code_k2 =
                    code_kx * code_kx + code_ky * code_ky + code_kz * code_kz;

                // CIC Deconvolution
                auto sinc = [](double freq, int mesh_size) {
                    if (freq == 0.0) return 1.0;
                    double arg = M_PI * freq / mesh_size;
                    return std::sin(arg) / arg;
                };

                // The CIC interpolation window in Fourier space is sinc^2 per
                // dimension
                double sinc_x = sinc(kx_freq, Mx);
                double sinc_y = sinc(ky_freq, My);
                double sinc_z = sinc(kz_freq, Mz);
                double cic_filter =
                    (sinc_x * sinc_x) * (sinc_y * sinc_y) * (sinc_z * sinc_z);

                std::complex<double> phi_k = -delta_k / code_k2;

                // Divide the potential by the CIC filter
                phi_k /= cic_filter;

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

    fftw_complex* disp_x_k_fftw =
        reinterpret_cast<fftw_complex*>(disp_x_k.data());
    fftw_complex* disp_y_k_fftw =
        reinterpret_cast<fftw_complex*>(disp_y_k.data());
    fftw_complex* disp_z_k_fftw =
        reinterpret_cast<fftw_complex*>(disp_z_k.data());

    fftw_plan backward_plan_x = fftw_plan_dft_c2r_3d(
        Mx, My, Mz, disp_x_k_fftw, field.dx.data(), FFTW_ESTIMATE);
    fftw_plan backward_plan_y = fftw_plan_dft_c2r_3d(
        Mx, My, Mz, disp_y_k_fftw, field.dy.data(), FFTW_ESTIMATE);
    fftw_plan backward_plan_z = fftw_plan_dft_c2r_3d(
        Mx, My, Mz, disp_z_k_fftw, field.dz.data(), FFTW_ESTIMATE);

    fftw_execute(backward_plan_x);
    fftw_execute(backward_plan_y);
    fftw_execute(backward_plan_z);

    fftw_destroy_plan(backward_plan_x);
    fftw_destroy_plan(backward_plan_y);
    fftw_destroy_plan(backward_plan_z);

    fftw_cleanup_threads();

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

// Evaluates the continuous Zel'dovich field at arbitrary coordinates
// using Cloud-In-Cell (CIC) interpolation
static Vec3 sample_zeldovich_cic(double qx, double qy, double qz,
                                 const ZeldovichField& zf, double cell_size,
                                 const Vec3& M) {
    int Mx = static_cast<int>(M.x);
    int My = static_cast<int>(M.y);
    int Mz = static_cast<int>(M.z);

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

    int ix0 = (ix + Mx) % Mx, ix1 = (ix + 1 + Mx) % Mx;
    int iy0 = (iy + My) % My, iy1 = (iy + 1 + My) % My;
    int iz0 = (iz + Mz) % Mz, iz1 = (iz + 1 + Mz) % Mz;

    int stride_x = static_cast<int>(My * Mz);
    int stride_y = static_cast<int>(Mz);

    auto interp = [&](const std::vector<double>& field) {
        return field[ix0 * stride_x + iy0 * stride_y + iz0] * w000 +
               field[ix1 * stride_x + iy0 * stride_y + iz0] * w100 +
               field[ix0 * stride_x + iy1 * stride_y + iz0] * w010 +
               field[ix1 * stride_x + iy1 * stride_y + iz0] * w110 +
               field[ix0 * stride_x + iy0 * stride_y + iz1] * w001 +
               field[ix1 * stride_x + iy0 * stride_y + iz1] * w101 +
               field[ix0 * stride_x + iy1 * stride_y + iz1] * w011 +
               field[ix1 * stride_x + iy1 * stride_y + iz1] * w111;
    };

    return {interp(zf.dx), interp(zf.dy), interp(zf.dz)};
}

void initialize_dm(SimState& state, const Config& config,
                   const ZeldovichField& zf) {
    const double baryon_fraction = config.omega_baryon / config.omega_m;
    const double dm_total_mass = 1.0 - baryon_fraction;
    if (config.use_glass) {
        // Find largest power-of-two divisor to tile the particles
        int N_1d = config.num_dm_particles_1d;
        int N_glass = 1;
        while (N_1d % (N_glass * 2) == 0 && N_glass * 2 <= 128) {
            N_glass *= 2;
        }
        N_glass = std::max(8, N_glass);

        double L_base = std::min(
            {config.domain_size.x, config.domain_size.y, config.domain_size.z});
        int sub_tiles = N_1d / N_glass;
        double L_block = L_base / sub_tiles;

        int n_tiles_x =
            static_cast<int>(std::round(config.domain_size.x / L_block));
        int n_tiles_y =
            static_cast<int>(std::round(config.domain_size.y / L_block));
        int n_tiles_z =
            static_cast<int>(std::round(config.domain_size.z / L_block));

        std::string file_glass = config.glass_folder_path + "/glassCube_" +
                                 std::to_string(N_glass) + ".hdf5";
        GlassData glass;
        try {
            glass = HDF5Reader::read_glass_cube(file_glass);
        } catch (const std::exception& e) {
            exit(EXIT_FAILURE);
        }

        // Derive actual particle mass from the tiled glass
        double total_dm_particles =
            glass.num_particles * n_tiles_x * n_tiles_y * n_tiles_z;
        double dm_particle_mass = dm_total_mass / total_dm_particles;

        double cell_size = config.domain_size.x / config.mesh_size.x;

        for (int tx = 0; tx < n_tiles_x; ++tx) {
            for (int ty = 0; ty < n_tiles_y; ++ty) {
                for (int tz = 0; tz < n_tiles_z; ++tz) {
                    for (size_t p = 0; p < glass.num_particles; ++p) {
                        // Unperturbed glass coordinates for DM
                        double qx = tx * L_block + glass.pos_x[p] * L_block;
                        double qy = ty * L_block + glass.pos_y[p] * L_block;
                        double qz = tz * L_block + glass.pos_z[p] * L_block;

                        Vec3 d = sample_zeldovich_cic(qx, qy, qz, zf, cell_size,
                                                      config.mesh_size);

                        double p_x = std::fmod(qx + d.x + config.domain_size.x,
                                               config.domain_size.x);
                        double p_y = std::fmod(qy + d.y + config.domain_size.y,
                                               config.domain_size.y);
                        double p_z = std::fmod(qz + d.z + config.domain_size.z,
                                               config.domain_size.z);

                        double v_x = state.hubble_param * d.x * zf.f;
                        double v_y = state.hubble_param * d.y * zf.f;
                        double v_z = state.hubble_param * d.z * zf.f;

                        state.dm.add_particle(p_x, p_y, p_z, v_x, v_y, v_z,
                                              dm_particle_mass);
                    }
                }
            }
        }
    } else {
        int M = config.mesh_size_1d;
        double cell_size = config.domain_size_1d / M;
        double spacing = config.domain_size_1d / config.num_dm_particles_1d;
        double dm_particle_mass =
            dm_total_mass / pow(config.num_dm_particles_1d, 3);

        int scale_x =
            static_cast<int>(config.domain_size.x / config.domain_size_1d);
        int scale_y =
            static_cast<int>(config.domain_size.y / config.domain_size_1d);
        int scale_z =
            static_cast<int>(config.domain_size.z / config.domain_size_1d);

        int N_part_x = config.num_dm_particles_1d * scale_x;
        int N_part_y = config.num_dm_particles_1d * scale_y;
        int N_part_z = config.num_dm_particles_1d * scale_z;

        for (int i = 0; i < N_part_x; ++i) {
            for (int j = 0; j < N_part_y; ++j) {
                for (int k = 0; k < N_part_z; ++k) {
                    double qx = (i + 0.5) * spacing;
                    double qy = (j + 0.5) * spacing;
                    double qz = (k + 0.5) * spacing;

                    Vec3 d = sample_zeldovich_cic(qx, qy, qz, zf, cell_size,
                                                  config.mesh_size);

                    double p_x = fmod(qx + d.x + config.domain_size.x,
                                      config.domain_size.x);
                    double p_y = fmod(qy + d.y + config.domain_size.y,
                                      config.domain_size.y);
                    double p_z = fmod(qz + d.z + config.domain_size.z,
                                      config.domain_size.z);

                    double v_x = state.hubble_param * d.x * zf.f;
                    double v_y = state.hubble_param * d.y * zf.f;
                    double v_z = state.hubble_param * d.z * zf.f;

                    state.dm.add_particle(p_x, p_y, p_z, v_x, v_y, v_z,
                                          dm_particle_mass);
                }
            }
        }
    }
}

void initialize_gas(SimState& state, const Config& config,
                    const ZeldovichField& zf) {
    if (config.hydro_method == HydroMethod::None) return;

    int Mx = static_cast<int>(config.mesh_size.x);
    int My = static_cast<int>(config.mesh_size.y);
    int Mz = static_cast<int>(config.mesh_size.z);

    double seed_metallicity =
        config.seed_metallicity_solar * constants::Z_SOLAR;

    const double initial_internal_energy =
        Thermodynamics::get_internal_energy_from_temp(
            config.initial_gas_temperature_k, state.scale_factor, config);

    double gas_total_mass = config.omega_baryon / config.omega_m;

    if (config.hydro_method == HydroMethod::Eulerian) {
        auto& gas = *state.gas;
        double mean_gas_rho = gas_total_mass / pow(config.domain_size_1d, 3);
        double inv_2dx = 1.0 / (2.0 * config.cell_size);

        for (int i = 0; i < Mx; ++i) {
            for (int j = 0; j < My; ++j) {
                for (int k = 0; k < Mz; ++k) {
                    size_t idx = static_cast<size_t>(i) * My * Mz + j * Mz + k;

                    // Calculate cell center
                    double qx = (i + 0.5) * config.cell_size;
                    double qy = (j + 0.5) * config.cell_size;
                    double qz = (k + 0.5) * config.cell_size;

                    // Sample local Zel'dovich displacement
                    Vec3 d = sample_zeldovich_cic(
                        qx, qy, qz, zf, config.cell_size, config.mesh_size);

                    double qx_p =
                        fmod(qx + config.cell_size, config.domain_size.x);
                    double qx_m =
                        fmod(qx - config.cell_size + config.domain_size.x,
                             config.domain_size.x);
                    double qy_p =
                        fmod(qy + config.cell_size, config.domain_size.y);
                    double qy_m =
                        fmod(qy - config.cell_size + config.domain_size.y,
                             config.domain_size.y);
                    double qz_p =
                        fmod(qz + config.cell_size, config.domain_size.z);
                    double qz_m =
                        fmod(qz - config.cell_size + config.domain_size.z,
                             config.domain_size.z);

                    // Central difference for divergence (sampled at
                    // neighboring cell centers)
                    Vec3 dx_ip1 = sample_zeldovich_cic(
                        qx_p, qy, qz, zf, config.cell_size, config.mesh_size);
                    Vec3 dx_im1 = sample_zeldovich_cic(
                        qx_m, qy, qz, zf, config.cell_size, config.mesh_size);
                    Vec3 dy_jp1 = sample_zeldovich_cic(
                        qx, qy_p, qz, zf, config.cell_size, config.mesh_size);
                    Vec3 dy_jm1 = sample_zeldovich_cic(
                        qx, qy_m, qz, zf, config.cell_size, config.mesh_size);
                    Vec3 dz_kp1 = sample_zeldovich_cic(
                        qx, qy, qz_p, zf, config.cell_size, config.mesh_size);
                    Vec3 dz_km1 = sample_zeldovich_cic(
                        qx, qy, qz_m, zf, config.cell_size, config.mesh_size);

                    double div = (dx_ip1.x - dx_im1.x) * inv_2dx +
                                 (dy_jp1.y - dy_jm1.y) * inv_2dx +
                                 (dz_kp1.z - dz_km1.z) * inv_2dx;

                    double rho = mean_gas_rho * (1.0 - div);
                    double safe_rho = std::max(rho, 1e-12);
                    gas.density.data[idx] = safe_rho;
                    gas.metal_density.data[idx] = safe_rho * seed_metallicity;

                    double vx = state.hubble_param * d.x * zf.f;
                    double vy = state.hubble_param * d.y * zf.f;
                    double vz = state.hubble_param * d.z * zf.f;

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
        if (config.use_glass) {
            int N_1d = config.num_gas_particles_1d;
            int N_glass = 1;
            while (N_1d % (N_glass * 2) == 0 && N_glass * 2 <= 128) {
                N_glass *= 2;
            }
            N_glass = std::max(8, N_glass);

            double L_base =
                std::min({config.domain_size.x, config.domain_size.y,
                          config.domain_size.z});
            int sub_tiles = N_1d / N_glass;
            double L_block = L_base / sub_tiles;

            int n_tiles_x =
                static_cast<int>(std::round(config.domain_size.x / L_block));
            int n_tiles_y =
                static_cast<int>(std::round(config.domain_size.y / L_block));
            int n_tiles_z =
                static_cast<int>(std::round(config.domain_size.z / L_block));

            std::string file_glass = config.glass_folder_path + "/glassCube_" +
                                     std::to_string(N_glass) + ".hdf5";
            GlassData glass;
            try {
                glass = HDF5Reader::read_glass_cube(file_glass);
            } catch (const std::exception& e) {
                exit(EXIT_FAILURE);
            }

            double total_gas_particles =
                glass.num_particles * n_tiles_x * n_tiles_y * n_tiles_z;
            double gas_particle_mass = gas_total_mass / total_gas_particles;

            // Shift the gas glass by half a spacing relative to DM to avoid
            // overlap
            double spacing = L_block / N_glass;
            double shift = 0.5 * spacing;

            double initial_h =
                0.5 * spacing * std::cbrt(config.mfm_target_neighbors);
            double cell_size = config.domain_size.x / config.mesh_size.x;

            for (int tx = 0; tx < n_tiles_x; ++tx) {
                for (int ty = 0; ty < n_tiles_y; ++ty) {
                    for (int tz = 0; tz < n_tiles_z; ++tz) {
                        for (size_t p = 0; p < glass.num_particles; ++p) {
                            // Unperturbed glass coordinates for Gas (Shifted
                            // with local block wrapping)
                            double qx =
                                tx * L_block +
                                std::fmod(glass.pos_x[p] * L_block + shift,
                                          L_block);
                            double qy =
                                ty * L_block +
                                std::fmod(glass.pos_y[p] * L_block + shift,
                                          L_block);
                            double qz =
                                tz * L_block +
                                std::fmod(glass.pos_z[p] * L_block + shift,
                                          L_block);

                            Vec3 d = sample_zeldovich_cic(
                                qx, qy, qz, zf, cell_size, config.mesh_size);

                            double p_x =
                                std::fmod(qx + d.x + config.domain_size.x,
                                          config.domain_size.x);
                            double p_y =
                                std::fmod(qy + d.y + config.domain_size.y,
                                          config.domain_size.y);
                            double p_z =
                                std::fmod(qz + d.z + config.domain_size.z,
                                          config.domain_size.z);

                            double v_x = state.hubble_param * d.x * zf.f;
                            double v_y = state.hubble_param * d.y * zf.f;
                            double v_z = state.hubble_param * d.z * zf.f;

                            state.mfm_gas->add_particle(
                                p_x, p_y, p_z, v_x, v_y, v_z, gas_particle_mass,
                                initial_internal_energy, initial_h,
                                seed_metallicity);
                        }
                    }
                }
            }
        } else {
            int scale_x =
                static_cast<int>(config.domain_size.x / config.domain_size_1d);
            int scale_y =
                static_cast<int>(config.domain_size.y / config.domain_size_1d);
            int scale_z =
                static_cast<int>(config.domain_size.z / config.domain_size_1d);

            int N_part_x = config.num_gas_particles_1d * scale_x;
            int N_part_y = config.num_gas_particles_1d * scale_y;
            int N_part_z = config.num_gas_particles_1d * scale_z;

            double spacing =
                config.domain_size_1d / config.num_gas_particles_1d;
            double gas_particle_mass =
                gas_total_mass / pow(config.num_gas_particles_1d, 3);

            // Initial guess for smoothing length
            double initial_h =
                0.5 * spacing * cbrt(config.mfm_target_neighbors);

            for (int i = 0; i < N_part_x; ++i) {
                for (int j = 0; j < N_part_y; ++j) {
                    for (int k = 0; k < N_part_z; ++k) {
                        // Shift the gas lattice by half a spacing relative to
                        // DM to avoid perfect overlap
                        double qx = i * spacing;
                        double qy = j * spacing;
                        double qz = k * spacing;

                        Vec3 d = sample_zeldovich_cic(
                            qx, qy, qz, zf, config.cell_size, config.mesh_size);

                        double p_x = fmod(qx + d.x + config.domain_size.x,
                                          config.domain_size.x);
                        double p_y = fmod(qy + d.y + config.domain_size.y,
                                          config.domain_size.y);
                        double p_z = fmod(qz + d.z + config.domain_size.z,
                                          config.domain_size.z);

                        double v_x = state.hubble_param * d.x * zf.f;
                        double v_y = state.hubble_param * d.y * zf.f;
                        double v_z = state.hubble_param * d.z * zf.f;

                        state.mfm_gas->add_particle(
                            p_x, p_y, p_z, v_x, v_y, v_z, gas_particle_mass,
                            initial_internal_energy, initial_h,
                            seed_metallicity);
                    }
                }
            }
        }
    }
}

void initialize_sod_shock_tube(SimState& state, const Config& config) {
    // Shock Tube parameters (Hernquist & Katz 1989)
    double rho_L = 1.0;
    double rho_R = 0.125;
    double P_L = 1.0;
    double P_R = 0.1;

    double Lx = config.domain_size.x;
    double Ly = config.domain_size.y;
    double Lz = config.domain_size.z;

    // The base 1D reference length for mapping a single glass block
    double L_base = config.domain_size_1d;

    double gamma = config.gamma;
    double u_L = P_L / (rho_L * (gamma - 1.0));
    double u_R = P_R / (rho_R * (gamma - 1.0));
    double v_x = 0.0, v_y = 0.0, v_z = 0.0;
    double seed_metallicity = 0.0;

    if (config.hydro_method == HydroMethod::Eulerian) {
        auto& gas = *state.gas;

        int Mx = static_cast<int>(config.mesh_size.x);
        int My = static_cast<int>(config.mesh_size.y);
        int Mz = static_cast<int>(config.mesh_size.z);

        for (int i = 0; i < Mx; ++i) {
            for (int j = 0; j < My; ++j) {
                for (int k = 0; k < Mz; ++k) {
                    size_t idx = static_cast<size_t>(i) * My * Mz + j * Mz + k;

                    double x = (i + 0.5) * config.cell_size;

                    // Shock interface is at the domain midpoint
                    double rho = (x < Lx / 2.0) ? rho_L : rho_R;
                    double u = (x < Lx / 2.0) ? u_L : u_R;

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
        if (config.use_glass) {
            int target_N = config.num_gas_particles_1d;
            int N_L = 8;
            if (target_N >= 96)
                N_L = 128;
            else if (target_N >= 48)
                N_L = 64;
            else if (target_N >= 24)
                N_L = 32;
            else if (target_N >= 12)
                N_L = 16;

            int N_R = N_L / 2;  // Right side is one tier lower

            std::string file_L = config.glass_folder_path + "/glassCube_" +
                                 std::to_string(N_L) + ".hdf5";
            std::string file_R = config.glass_folder_path + "/glassCube_" +
                                 std::to_string(N_R) + ".hdf5";

            GlassData glass_L, glass_R;
            try {
                glass_L = HDF5Reader::read_glass_cube(file_L);
                glass_R = HDF5Reader::read_glass_cube(file_R);
            } catch (const std::exception& e) {
                exit(EXIT_FAILURE);
            }

            // Particle mass is defined by the density of the base block
            double particle_mass = rho_L * (L_base * L_base * L_base) /
                                   static_cast<double>(glass_L.num_particles);

            double target_neighbors = config.mfm_target_neighbors;
            double radius_factor =
                std::pow((3.0 * target_neighbors) / (4.0 * M_PI), 1.0 / 3.0);

            // True spacing assuming the glass fills L_base
            double dx_L =
                L_base * std::pow(1.0 / glass_L.num_particles, 1.0 / 3.0);
            double dx_R =
                L_base * std::pow(1.0 / glass_R.num_particles, 1.0 / 3.0);

            double h_L = dx_L * radius_factor;
            double h_R = dx_R * radius_factor;

            // Calculate the number of times we must tile the base block
            int n_tiles_x = static_cast<int>(std::ceil(Lx / L_base));
            int n_tiles_y = static_cast<int>(std::ceil(Ly / L_base));
            int n_tiles_z = static_cast<int>(std::ceil(Lz / L_base));

            // LEFT SIDE (High Density, High Pressure)
            // Only tile up to the domain midpoint
            int left_tiles_x = static_cast<int>(std::ceil((Lx / 2.0) / L_base));
            for (int tx = 0; tx < left_tiles_x; ++tx) {
                for (int ty = 0; ty < n_tiles_y; ++ty) {
                    for (int tz = 0; tz < n_tiles_z; ++tz) {
                        for (size_t p = 0; p < glass_L.num_particles; ++p) {
                            double p_x =
                                tx * L_base + glass_L.pos_x[p] * L_base;
                            double p_y =
                                ty * L_base + glass_L.pos_y[p] * L_base;
                            double p_z =
                                tz * L_base + glass_L.pos_z[p] * L_base;

                            // Clip particles crossing the mid-plane or
                            // extending past Y/Z boundaries
                            if (p_x < Lx / 2.0 && p_y < Ly && p_z < Lz) {
                                state.mfm_gas->add_particle(
                                    p_x, p_y, p_z, v_x, v_y, v_z, particle_mass,
                                    u_L, h_L, seed_metallicity);
                            }
                        }
                    }
                }
            }

            // RIGHT SIDE (Low Density, Low Pressure)
            // Start tiling from the block intersecting the domain midpoint
            int start_right_tx =
                static_cast<int>(std::floor((Lx / 2.0) / L_base));
            for (int tx = start_right_tx; tx < n_tiles_x; ++tx) {
                for (int ty = 0; ty < n_tiles_y; ++ty) {
                    for (int tz = 0; tz < n_tiles_z; ++tz) {
                        for (size_t p = 0; p < glass_R.num_particles; ++p) {
                            double p_x =
                                tx * L_base + glass_R.pos_x[p] * L_base;
                            double p_y =
                                ty * L_base + glass_R.pos_y[p] * L_base;
                            double p_z =
                                tz * L_base + glass_R.pos_z[p] * L_base;

                            // Clip particles on the left of the mid-plane or
                            // past the outer boundaries
                            if (p_x >= Lx / 2.0 && p_x < Lx && p_y < Ly &&
                                p_z < Lz) {
                                state.mfm_gas->add_particle(
                                    p_x, p_y, p_z, v_x, v_y, v_z, particle_mass,
                                    u_R, h_R, seed_metallicity);
                            }
                        }
                    }
                }
            }
        } else {
            // LATTICE INITIALIZATION
            double dx_L = L_base / config.num_gas_particles_1d;

            // Stretch ONLY the X-axis by the full 1D density ratio
            double dx_R_x = dx_L * (rho_L / rho_R);

            double particle_mass = rho_L * (dx_L * dx_L * dx_L);

            // Pre-calculate valid Y and Z coordinates to fit exactly within Ly
            // and Lz
            std::vector<double> y_coords;
            for (double y = dx_L / 2.0; y < Ly; y += dx_L)
                y_coords.push_back(y);

            std::vector<double> z_coords;
            for (double z = dx_L / 2.0; z < Lz; z += dx_L)
                z_coords.push_back(z);

            // LEFT SIDE
            for (double x = dx_L / 2.0; x < Lx / 2.0; x += dx_L) {
                for (double y : y_coords) {
                    for (double z : z_coords) {
                        double initial_h = 1.2 * dx_L;
                        state.mfm_gas->add_particle(
                            x, y, z, v_x, v_y, v_z, particle_mass, u_L,
                            initial_h, seed_metallicity);
                    }
                }
            }

            // RIGHT SIDE (X spacing stretches by a factor of 8, Y/Z stay
            // locked)
            for (double x = Lx / 2.0 + dx_R_x / 2.0; x < Lx; x += dx_R_x) {
                for (double y : y_coords) {
                    for (double z : z_coords) {
                        // The smoothing length must encompass the stretched
                        // X-axis
                        double initial_h = 1.2 * std::max(dx_R_x, dx_L);
                        state.mfm_gas->add_particle(
                            x, y, z, v_x, v_y, v_z, particle_mass, u_R,
                            initial_h, seed_metallicity);
                    }
                }
            }
        }
    }
}

void initialize_sedov_blastwave(SimState& state, const Config& config) {
    // Sedov-Taylor parameters
    double E_total = 1.0;
    double rho_bg = 1.0;
    double P_bg = 1e-6;  // Almost zero pressure for the background
    double gamma = config.gamma;
    double u_bg = P_bg / (rho_bg * (gamma - 1.0));
    double seed_metallicity = 0.0;

    // Explicit domain dimensions
    double Lx = config.domain_size.x;
    double Ly = config.domain_size.y;
    double Lz = config.domain_size.z;
    double L_base = config.domain_size_1d;

    // True geometric center of the rectangular domain
    double cx = Lx / 2.0;
    double cy = Ly / 2.0;
    double cz = Lz / 2.0;

    int N_res = (config.hydro_method == HydroMethod::Eulerian)
                    ? config.mesh_size_1d
                    : config.num_gas_particles_1d;
    double dx_base = L_base / static_cast<double>(N_res);
    double R_inj = 2.0 * dx_base;

    if (config.hydro_method == HydroMethod::Eulerian) {
        auto& gas = *state.gas;
        int Mx = static_cast<int>(config.mesh_size.x);
        int My = static_cast<int>(config.mesh_size.y);
        int Mz = static_cast<int>(config.mesh_size.z);

        // Initialize uniform cold background
        for (int i = 0; i < Mx; ++i) {
            for (int j = 0; j < My; ++j) {
                for (int k = 0; k < Mz; ++k) {
                    size_t idx = static_cast<size_t>(i) * My * Mz + j * Mz + k;
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
        double cell_vol = dx_base * dx_base * dx_base;
        double inj_mass = 0.0;
        std::vector<size_t> inj_indices;

        for (int i = 0; i < Mx; ++i) {
            for (int j = 0; j < My; ++j) {
                for (int k = 0; k < Mz; ++k) {
                    double x = (i + 0.5) * dx_base;
                    double y = (j + 0.5) * dx_base;
                    double z = (k + 0.5) * dx_base;

                    // Distance to true rectangular center
                    double r =
                        std::sqrt((x - cx) * (x - cx) + (y - cy) * (y - cy) +
                                  (z - cz) * (z - cz));

                    if (r <= R_inj) {
                        inj_indices.push_back(static_cast<size_t>(i) * My * Mz +
                                              j * Mz + k);
                        inj_mass += rho_bg * cell_vol;
                    }
                }
            }
        }

        double u_inj = E_total / inj_mass;
        for (size_t idx : inj_indices) {
            gas.internal_energy.data[idx] += rho_bg * u_inj;
            gas.energy.data[idx] += rho_bg * u_inj;
        }
        gas.update_primitive_variables(state.scale_factor);

    } else if (config.hydro_method == HydroMethod::MFM) {
        if (!config.use_glass) {
            int scale_x = static_cast<int>(std::round(Lx / L_base));
            int scale_y = static_cast<int>(std::round(Ly / L_base));
            int scale_z = static_cast<int>(std::round(Lz / L_base));

            int N_x = N_res * scale_x;
            int N_y = N_res * scale_y;
            int N_z = N_res * scale_z;

            double gas_particle_mass = rho_bg * (dx_base * dx_base * dx_base);
            double initial_h = 2.0 * dx_base;

            constexpr bool spread_over_8_particles = true;

            // Calculate center indices per axis
            auto get_bounds = [](int N, bool spread) {
                if (spread && N % 2 == 0)
                    return std::make_pair((N / 2) - 1, N / 2);
                return std::make_pair(N / 2, N / 2);
            };

            auto [ix_min, ix_max] = get_bounds(N_x, spread_over_8_particles);
            auto [iy_min, iy_max] = get_bounds(N_y, spread_over_8_particles);
            auto [iz_min, iz_max] = get_bounds(N_z, spread_over_8_particles);

            int num_central = (ix_max - ix_min + 1) * (iy_max - iy_min + 1) *
                              (iz_max - iz_min + 1);

            double offset = spread_over_8_particles ? 0.5 : 0.0;
            double injected_E = E_total / num_central;
            double added_u = injected_E / gas_particle_mass;

            for (int i = 0; i < N_x; ++i) {
                for (int j = 0; j < N_y; ++j) {
                    for (int k = 0; k < N_z; ++k) {
                        double p_x = (i + offset) * dx_base;
                        double p_y = (j + offset) * dx_base;
                        double p_z = (k + offset) * dx_base;

                        double particle_u = u_bg;
                        if (i >= ix_min && i <= ix_max && j >= iy_min &&
                            j <= iy_max && k >= iz_min && k <= iz_max) {
                            particle_u += added_u;
                        }

                        state.mfm_gas->add_particle(
                            p_x, p_y, p_z, 0.0, 0.0, 0.0, gas_particle_mass,
                            particle_u, initial_h, seed_metallicity);
                    }
                }
            }
        } else {
            // Glass Initialization (Two-Pass Approach)
            int N_glass = std::max(
                8, static_cast<int>(std::pow(2, std::round(std::log2(N_res)))));
            N_glass = std::min(N_glass, 128);

            std::string file_glass = config.glass_folder_path + "/glassCube_" +
                                     std::to_string(N_glass) + ".hdf5";
            GlassData glass;
            try {
                glass = HDF5Reader::read_glass_cube(file_glass);
            } catch (const std::exception& e) {
                return;
            }

            double vol_base = L_base * L_base * L_base;
            double gas_particle_mass =
                rho_bg * vol_base / static_cast<double>(glass.num_particles);

            // Tile the glass and collect valid coordinates
            struct ParticleData {
                Vec3 pos;
                double r2;
            };
            std::vector<ParticleData> tiled_particles;

            int n_tiles_x = static_cast<int>(std::ceil(Lx / L_base));
            int n_tiles_y = static_cast<int>(std::ceil(Ly / L_base));
            int n_tiles_z = static_cast<int>(std::ceil(Lz / L_base));

            for (int tx = 0; tx < n_tiles_x; ++tx) {
                for (int ty = 0; ty < n_tiles_y; ++ty) {
                    for (int tz = 0; tz < n_tiles_z; ++tz) {
                        for (size_t p = 0; p < glass.num_particles; ++p) {
                            double px = tx * L_base + glass.pos_x[p] * L_base;
                            double py = ty * L_base + glass.pos_y[p] * L_base;
                            double pz = tz * L_base + glass.pos_z[p] * L_base;

                            if (px < Lx && py < Ly && pz < Lz) {
                                double dx2 = (px - cx) * (px - cx);
                                double dy2 = (py - cy) * (py - cy);
                                double dz2 = (pz - cz) * (pz - cz);
                                tiled_particles.push_back(
                                    {{px, py, pz}, dx2 + dy2 + dz2});
                            }
                        }
                    }
                }
            }

            // Find the central N_inject particles
            int N_inject = 15;
            constexpr bool smooth_injection = true;

            N_inject =
                std::min(N_inject, static_cast<int>(tiled_particles.size()));

            std::partial_sort(tiled_particles.begin(),
                              tiled_particles.begin() + N_inject,
                              tiled_particles.end(),
                              [](const ParticleData& a, const ParticleData& b) {
                                  return a.r2 < b.r2;
                              });

            double centroid_x = 0.0, centroid_y = 0.0, centroid_z = 0.0;
            for (int i = 0; i < N_inject; ++i) {
                centroid_x += tiled_particles[i].pos.x;
                centroid_y += tiled_particles[i].pos.y;
                centroid_z += tiled_particles[i].pos.z;
            }
            centroid_x /= N_inject;
            centroid_y /= N_inject;
            centroid_z /= N_inject;

            // Define the shift required to move the centroid to (cx, cy, cz)
            const double dx_center = centroid_x - cx;
            const double dy_center = centroid_y - cy;
            const double dz_center = centroid_z - cz;

            // Shift ALL particles, apply periodic wrapping, and update r2
            for (size_t i = 0; i < tiled_particles.size(); ++i) {
                tiled_particles[i].pos.x -= dx_center;
                tiled_particles[i].pos.y -= dy_center;
                tiled_particles[i].pos.z -= dz_center;

                tiled_particles[i].pos.x =
                    fmod(tiled_particles[i].pos.x + config.domain_size.x,
                         config.domain_size.x);
                tiled_particles[i].pos.y =
                    fmod(tiled_particles[i].pos.y + config.domain_size.y,
                         config.domain_size.y);
                tiled_particles[i].pos.z =
                    fmod(tiled_particles[i].pos.z + config.domain_size.z,
                         config.domain_size.z);

                // Recalculate distance to the center
                double dx = tiled_particles[i].pos.x - cx;
                double dy = tiled_particles[i].pos.y - cy;
                double dz = tiled_particles[i].pos.z - cz;

                tiled_particles[i].r2 = dx * dx + dy * dy + dz * dz;
            }

            // Re-sort ONLY the injected group based on their new centered r2
            // values to ensure r_max is perfectly accurate for the Gaussian
            // standard deviation
            std::sort(tiled_particles.begin(),
                      tiled_particles.begin() + N_inject,
                      [](const ParticleData& a, const ParticleData& b) {
                          return a.r2 < b.r2;
                      });

            double r_max = std::sqrt(tiled_particles[N_inject - 1].r2);
            double sigma = r_max / 2.5;
            double two_sigma2 = 2.0 * sigma * sigma;

            double total_weight = 0.0;
            std::vector<double> weights(N_inject, 0.0);
            for (int i = 0; i < N_inject; ++i) {
                weights[i] = smooth_injection
                                 ? std::exp(-tiled_particles[i].r2 / two_sigma2)
                                 : 1.0 / static_cast<double>(N_inject);
                total_weight += weights[i];
            }

            // Map injected energy to a flat array aligned with tiled_particles
            std::vector<double> injected_u_map(tiled_particles.size(), 0.0);
            for (int i = 0; i < N_inject; ++i) {
                double E_i = E_total * (weights[i] / total_weight);
                injected_u_map[i] = E_i / gas_particle_mass;
            }

            double dx_grid = L_base / static_cast<double>(N_glass);
            double initial_h = 2.0 * dx_grid;

            // Push to simulation state
            for (size_t i = 0; i < tiled_particles.size(); ++i) {
                double particle_u = u_bg + injected_u_map[i];
                state.mfm_gas->add_particle(
                    tiled_particles[i].pos.x, tiled_particles[i].pos.y,
                    tiled_particles[i].pos.z, 0.0, 0.0, 0.0, gas_particle_mass,
                    particle_u, initial_h, seed_metallicity);
            }
        }
    }
}

void initialize_adiabatic_expansion(SimState& state, const Config& config) {
    const double initial_internal_energy =
        Thermodynamics::get_internal_energy_from_temp(
            config.initial_gas_temperature_k, state.scale_factor, config);

    double seed_metallicity = 0.0;
    double v_x = 0.0, v_y = 0.0, v_z = 0.0;

    double gas_total_mass = config.omega_baryon / config.omega_m;

    if (config.hydro_method == HydroMethod::Eulerian) {
        auto& gas = *state.gas;

        // Extract multi-axis mesh dimensions
        int Mx = static_cast<int>(config.mesh_size.x);
        int My = static_cast<int>(config.mesh_size.y);
        int Mz = static_cast<int>(config.mesh_size.z);

        // Calculate the uniform background density based on the base volume
        double mean_gas_rho =
            gas_total_mass / std::pow(config.domain_size_1d, 3.0);

        for (int i = 0; i < Mx; ++i) {
            for (int j = 0; j < My; ++j) {
                for (int k = 0; k < Mz; ++k) {
                    // Correct 3D flattening for rectangular arrays
                    size_t idx = static_cast<size_t>(i) * My * Mz + j * Mz + k;

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
        // Calculate integer scaling multipliers for each axis
        int scale_x =
            static_cast<int>(config.domain_size.x / config.domain_size_1d);
        int scale_y =
            static_cast<int>(config.domain_size.y / config.domain_size_1d);
        int scale_z =
            static_cast<int>(config.domain_size.z / config.domain_size_1d);

        // Scale the particle loops
        int N_part_1d = config.num_gas_particles_1d;
        int N_part_x = N_part_1d * scale_x;
        int N_part_y = N_part_1d * scale_y;
        int N_part_z = N_part_1d * scale_z;

        double spacing = config.domain_size_1d / N_part_1d;

        // Particle mass relies solely on the base 1D box setup
        double gas_particle_mass =
            gas_total_mass / (N_part_1d * N_part_1d * N_part_1d);

        // Initial guess for smoothing length
        double initial_h = 1.2 * spacing;

        for (int i = 0; i < N_part_x; ++i) {
            for (int j = 0; j < N_part_y; ++j) {
                for (int k = 0; k < N_part_z; ++k) {
                    // Shift the gas lattice by half a spacing relative to DM
                    double qx = i * spacing;
                    double qy = j * spacing;
                    double qz = k * spacing;

                    double p_x =
                        fmod(qx + config.domain_size.x, config.domain_size.x);
                    double p_y =
                        fmod(qy + config.domain_size.y, config.domain_size.y);
                    double p_z =
                        fmod(qz + config.domain_size.z, config.domain_size.z);

                    state.mfm_gas->add_particle(
                        p_x, p_y, p_z, v_x, v_y, v_z, gas_particle_mass,
                        initial_internal_energy, initial_h, seed_metallicity);
                }
            }
        }
    }
}

void initialize_kelvin_helmholtz(SimState& state, const Config& config) {
    if (config.hydro_method != HydroMethod::MFM) return;

    // KH Instability Parameters (McNally et al. 2012 / Hopkins 2015)
    double rho_1 = 1.0;
    double rho_2 = 2.0;
    double delta_rho = 0.5 * (rho_2 - rho_1);

    double v_1 = -0.5;
    double v_2 = 0.5;
    double delta_v = 0.5 * (v_2 - v_1);

    double delta_y = 0.025;
    double P = 2.5;
    double gamma = config.gamma;
    double seed_metallicity = 0.0;

    double Lx = config.domain_size.x;
    double Ly = config.domain_size.y;
    double Lz = config.domain_size.z;
    double L_base = config.domain_size_1d;

    int N_part = config.num_gas_particles_1d;
    double spacing_base = L_base / N_part;

    // Calculate integer scaling multipliers for each axis
    int N_x = N_part * static_cast<int>(std::round(Lx / L_base));
    int N_y = N_part * static_cast<int>(std::round(Ly / L_base));
    int N_z = N_part * static_cast<int>(std::round(Lz / L_base));

    double spacing_x = Lx / N_x;
    double spacing_y = Ly / N_y;
    double spacing_z = Lz / N_z;
    double cell_volume = spacing_x * spacing_y * spacing_z;

    double initial_h = 1.2 * std::max({spacing_x, spacing_y, spacing_z});

    for (int i = 0; i < N_x; ++i) {
        for (int j = 0; j < N_y; ++j) {
            for (int k = 0; k < N_z; ++k) {
                double x = (i + 0.5) * spacing_x;
                double y = (j + 0.5) * spacing_y;
                double z = (k + 0.5) * spacing_z;

                // Normalize coordinates to [0, 1) for the profile definitions
                double x_frac = x / Lx;
                double y_frac = y / Ly;

                double rho, v_x;

                // Evaluate smoothed density and x-velocity profiles
                if (y_frac < 0.25) {
                    double exp_term = std::exp((y_frac - 0.25) / delta_y);
                    rho = rho_1 + delta_rho * exp_term;
                    v_x = v_1 + delta_v * exp_term;
                } else if (y_frac <= 0.5) {
                    double exp_term = std::exp((0.25 - y_frac) / delta_y);
                    rho = rho_2 - delta_rho * exp_term;
                    v_x = v_2 - delta_v * exp_term;
                } else if (y_frac <= 0.75) {
                    double exp_term = std::exp((y_frac - 0.75) / delta_y);
                    rho = rho_2 - delta_rho * exp_term;
                    v_x = v_2 - delta_v * exp_term;
                } else {
                    double exp_term = std::exp((0.75 - y_frac) / delta_y);
                    rho = rho_1 + delta_rho * exp_term;
                    v_x = v_1 + delta_v * exp_term;
                }

                // Apply initial y-velocity seed mode
                double v_y = 0.01 * std::sin(4.0 * M_PI * x_frac);
                double v_z = 0.0;

                double particle_mass = rho * cell_volume;
                double u = P / ((gamma - 1.0) * rho);

                state.mfm_gas->add_particle(x, y, z, v_x, v_y, v_z,
                                            particle_mass, u, initial_h,
                                            seed_metallicity);
            }
        }
    }
}

void initialize_evrard_collapse(SimState& state, const Config& config) {
    bool use_gas_particles =
        (config.num_gas_particles_1d > config.num_dm_particles_1d);

    if (config.num_gas_particles_1d > 0 &&
        config.hydro_method != HydroMethod::MFM)
        return;

    // Evrard Collapse Parameters (Hopkins 2015, Section 4.5.1)
    double M_total =
        use_gas_particles
            ? config.omega_baryon / config.omega_m
            : (config.omega_m - config.omega_baryon) / config.omega_m;
    double R_sphere = 0.35;
    double u_initial = 0.05 * (config.G * M_total) / R_sphere;
    double seed_metallicity = 0.0;

    // Geometric center of the domain
    double cx = config.domain_size.x / 2.0;
    double cy = config.domain_size.y / 2.0;
    double cz = config.domain_size.z / 2.0;

    struct RawParticle {
        double x, y, z, ru;
    };
    std::vector<RawParticle> uniform_points;

    // Define the target 1D resolution for the bounding box
    int N_res =
        std::max(config.num_gas_particles_1d, config.num_dm_particles_1d);

    // Generate a uniform sphere of points (Glass or Lattice)
    if (config.use_glass) {
        int N_glass = std::max(
            8, static_cast<int>(std::pow(2, std::round(std::log2(N_res)))));
        N_glass = std::min(N_glass, 128);

        std::string file_glass = config.glass_folder_path + "/glassCube_" +
                                 std::to_string(N_glass) + ".hdf5";
        GlassData glass;
        try {
            glass = HDF5Reader::read_glass_cube(file_glass);
        } catch (const std::exception& e) {
            return;
        }

        // Scale the length of one glass block so its internal number density
        // matches the target N_res^3 over the (2 * R_sphere)^3 bounding box
        double L_block =
            (2.0 * R_sphere) * (static_cast<double>(N_glass) / N_res);

        // Calculate required tiles to cover the bounding box
        int n_tiles =
            static_cast<int>(std::ceil(static_cast<double>(N_res) / N_glass));

        for (int tx = 0; tx < n_tiles; ++tx) {
            for (int ty = 0; ty < n_tiles; ++ty) {
                for (int tz = 0; tz < n_tiles; ++tz) {
                    for (size_t p = 0; p < glass.num_particles; ++p) {
                        double px = -R_sphere + (tx + glass.pos_x[p]) * L_block;
                        double py = -R_sphere + (ty + glass.pos_y[p]) * L_block;
                        double pz = -R_sphere + (tz + glass.pos_z[p]) * L_block;

                        double ru = std::sqrt(px * px + py * py + pz * pz);
                        // The radius check natively discards the overhang of
                        // the fractional tiles
                        if (ru <= R_sphere && ru > 1e-8) {
                            uniform_points.push_back({px, py, pz, ru});
                        }
                    }
                }
            }
        }
    } else {
        double spacing = (2.0 * R_sphere) / N_res;
        for (int i = 0; i < N_res; ++i) {
            for (int j = 0; j < N_res; ++j) {
                for (int k = 0; k < N_res; ++k) {
                    double px = -R_sphere + (i + 0.5) * spacing;
                    double py = -R_sphere + (j + 0.5) * spacing;
                    double pz = -R_sphere + (k + 0.5) * spacing;

                    double ru = std::sqrt(px * px + py * py + pz * pz);
                    if (ru <= R_sphere && ru > 1e-8) {
                        uniform_points.push_back({px, py, pz, ru});
                    }
                }
            }
        }
    }

    // Map uniform sphere to the 1/r density profile and push to state
    double particle_mass = M_total / static_cast<double>(uniform_points.size());

    for (const auto& p : uniform_points) {
        // Radial transformation to achieve rho(r) ~ 1/r
        double r_target = R_sphere * std::pow(p.ru / R_sphere, 1.5);

        // Normalize direction vectors
        double nx = p.x / p.ru;
        double ny = p.y / p.ru;
        double nz = p.z / p.ru;

        // Shift to geometric center of the domain and apply periodic bounds
        double final_x = std::fmod(cx + nx * r_target + config.domain_size.x,
                                   config.domain_size.x);
        double final_y = std::fmod(cy + ny * r_target + config.domain_size.y,
                                   config.domain_size.y);
        double final_z = std::fmod(cz + nz * r_target + config.domain_size.z,
                                   config.domain_size.z);

        // Theoretical density at the new target radius
        double target_rho =
            M_total / (2.0 * M_PI * R_sphere * R_sphere * r_target);

        // Initial guess for smoothing length
        double initial_h = 0.5 * std::cbrt(config.mfm_target_neighbors *
                                           particle_mass / target_rho);

        if (use_gas_particles) {
            state.mfm_gas->add_particle(final_x, final_y, final_z, 0.0, 0.0,
                                        0.0, particle_mass, u_initial,
                                        initial_h, seed_metallicity);
        } else {
            state.dm.add_particle(final_x, final_y, final_z, 0.0, 0.0, 0.0,
                                  particle_mass);
        }
    }
}

void initialize_soundwave(SimState& state, const Config& config) {
    if (config.hydro_method != HydroMethod::MFM) return;

    // Explicit domain dimensions
    double Lx = config.domain_size.x;
    double Ly = config.domain_size.y;
    double Lz = config.domain_size.z;
    double L_base = config.domain_size_1d;

    // OG3 Soundwave Parameters (Groth et al. 2023)
    double rho_0 = 1.0;
    double delta_rho = 1e-4;
    double c_s = 2.0 / 3.0;
    double gamma = config.gamma;
    double P_0 = (c_s * c_s * rho_0) / gamma;

    // Wavenumber for a single period across the X-axis
    double k_wave = 2.0 * M_PI / Lx;
    double seed_metallicity = 0.0;

    int N_part = config.num_gas_particles_1d;
    double spacing = L_base / N_part;

    // Calculate integer scaling multipliers for each axis
    int N_x = N_part * static_cast<int>(std::round(Lx / L_base));
    int N_y = N_part * static_cast<int>(std::round(Ly / L_base));
    int N_z = N_part * static_cast<int>(std::round(Lz / L_base));

    // Constant particle mass based on the unperturbed base background
    double gas_particle_mass =
        (rho_0 * L_base * L_base * L_base) / (N_part * N_part * N_part);
    double initial_h = 1.2 * spacing;

    for (int i = 0; i < N_x; ++i) {
        for (int j = 0; j < N_y; ++j) {
            for (int k_idx = 0; k_idx < N_z; ++k_idx) {
                // Unperturbed Lagrangian coordinates
                double qx = (i + 0.5) * spacing;
                double qy = (j + 0.5) * spacing;
                double qz = (k_idx + 0.5) * spacing;

                // Apply Lagrangian spatial shift
                double p_x =
                    qx + (delta_rho / (k_wave * rho_0)) * std::cos(k_wave * qx);
                double p_y = qy;
                double p_z = qz;

                // Wrap periodic boundaries along the specific axis
                p_x = std::fmod(p_x + Lx, Lx);
                p_y = std::fmod(p_y + Ly, Ly);
                p_z = std::fmod(p_z + Lz, Lz);

                // Thermodynamic and kinematic perturbations
                double sin_kqx = std::sin(k_wave * qx);

                // Left-moving wave velocity perturbation
                double v_x = -c_s * (delta_rho / rho_0) * sin_kqx;
                double v_y = 0.0;
                double v_z = 0.0;

                // Pressure perturbation
                double P = P_0 + (c_s * c_s) * delta_rho * sin_kqx;
                double rho_local = rho_0 + delta_rho * sin_kqx;

                double u = P / ((gamma - 1.0) * rho_local);

                state.mfm_gas->add_particle(p_x, p_y, p_z, v_x, v_y, v_z,
                                            gas_particle_mass, u, initial_h,
                                            seed_metallicity);
            }
        }
    }
}

SimState initialize_state(Config& config) {
    SimState state(config);
    state.total_time = 0;

    Cosmology::update(state, config);

    if (config.initial_setup == InitialSetup::Cosmological) {
        ZeldovichField z_field =
            compute_zeldovich_field(state.scale_factor, config);

        initialize_dm(state, config, z_field);
        initialize_gas(state, config, z_field);
    } else if (config.initial_setup == InitialSetup::SodShockTube) {
        initialize_sod_shock_tube(state, config);
    } else if (config.initial_setup == InitialSetup::AdiabaticExpansion) {
        initialize_adiabatic_expansion(state, config);
    } else if (config.initial_setup == InitialSetup::SedovBlastwave) {
        initialize_sedov_blastwave(state, config);
    } else if (config.initial_setup == InitialSetup::SoundWave) {
        initialize_soundwave(state, config);
    } else if (config.initial_setup == InitialSetup::KelvinHelmholtz) {
        initialize_kelvin_helmholtz(state, config);
    } else if (config.initial_setup == InitialSetup::EvrardCollapse) {
        initialize_evrard_collapse(state, config);
    }

    Integrator::initialize(state, config);

    return state;
}
