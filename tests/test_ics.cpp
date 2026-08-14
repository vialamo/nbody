#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "config.h"
#include "constants.h"
#include "cooling.h"
#include "ics.h"

// ---------------------------------------------------------------------
// Thermodynamics & Initial Conditions Tests
// ---------------------------------------------------------------------
TEST_CASE("Initial internal energy scales correctly with code units",
          "[hydro][thermodynamics][units]") {
    Config config;
    config.gamma = 5.0 / 3.0;
    config.primordial_mu = 1.22;

    // Let the config initialize its derived units naturally
    config.compute_derived_data();

    double T_kelvin = 50.0;
    double a = 1.0;

    SECTION("Halving the velocity unit quadruples the resulting code energy") {
        // Get the base energy using the default box size (default velocity
        // unit)
        double u_code_base =
            Cooling::get_internal_energy_from_temp(T_kelvin, a, config);

        // To safely halve the velocity unit, we halve the box_size_mpc.
        // As seen in init_derived_units, unit_velocity_kms scales linearly with
        // unit_length_mpc, which scales linearly with box_size_mpc.
        config.box_size_mpc /= 2.0;

        // Force the Config to rebuild the thermodynamics bridge naturally
        config.compute_derived_data();

        // Calculate the new energy
        double u_code_half_v =
            Cooling::get_internal_energy_from_temp(T_kelvin, a, config);

        // Because V is squared in the denominator of the conversion,
        // halving V (via the box size) must exactly quadruple the code energy.
        REQUIRE(u_code_half_v == Catch::Approx(u_code_base * 4.0));
    }

    SECTION("Expected mathematical output is physically accurate") {
        double u_code =
            Cooling::get_internal_energy_from_temp(T_kelvin, a, config);

        // Calculate the exact physical specific internal energy (in erg/g)
        double u_phys_cgs =
            (constants::K_B_CGS * T_kelvin) /
            ((config.gamma - 1.0) * config.primordial_mu * constants::M_P_CGS);

        // Read the EXACT velocity unit the config actually generated
        double v_unit_cgs = config.unit_velocity_kms * 1.0e5;

        // Compute the dynamic expected answer
        double expected_u_code = u_phys_cgs / (v_unit_cgs * v_unit_cgs);

        REQUIRE(u_code == Catch::Approx(expected_u_code).margin(1e-13));
    }
}

TEST_CASE("initialize_state produces physically sound macro-states",
          "[ics][system]") {
    // Manually build the configuration
    Config config;
    config.mesh_size = 32;
    config.box_size_mpc = 100.0;
    config.num_particles_1d = 32;
    config.use_PM = true;
    config.use_PP = false;

    config.compute_derived_data();

    // Initialize the simulation state
    SimState state = initialize_state(config);

    SECTION("Conservation of Mass") {
        double total_dm_mass = state.dm.num_particles * config.dm_particle_mass;

        double total_gas_mass = 0.0;
        int N3 = config.mesh_size * config.mesh_size * config.mesh_size;
        for (int i = 0; i < N3; ++i) {
            // Integrate density over Eulerian volume
            total_gas_mass +=
                state.gas->get_density().data[i] * config.cell_volume;
        }

        // The universe must have exactly 1.0 code mass
        REQUIRE(total_dm_mass + total_gas_mass ==
                Catch::Approx(1.0).margin(1e-5));
        REQUIRE(total_gas_mass ==
                Catch::Approx(config.gas_total_mass).margin(1e-5));
    }

    SECTION("Conservation of Momentum") {
        double p_dm_x = 0.0;
        for (size_t i = 0; i < state.dm.num_particles; ++i) {
            p_dm_x += state.dm.vel_x[i] * state.dm.mass[i];
        }

        double p_gas_x_density = 0.0;
        int N3 = config.mesh_size * config.mesh_size * config.mesh_size;
        for (int i = 0; i < N3; ++i) {
            p_gas_x_density += state.gas->get_momentum_x().data[i];
        }

        // Convert momentum density to actual physical momentum
        double physical_p_gas_x = p_gas_x_density * config.cell_volume;

        // In a periodic box, the sum of all momenta must be zero
        REQUIRE(std::abs(p_dm_x) < 1e-10);
        REQUIRE(std::abs(physical_p_gas_x) < 1e-3);
    }

    SECTION("Velocities are mathematically bound to cosmic expansion") {
        double max_v_dm = 0.0;
        for (size_t i = 0; i < state.dm.num_particles; ++i) {
            max_v_dm = std::max(max_v_dm, std::abs(state.dm.vel_x[i]));
        }

        // At a=0.02, H(a) is ~129 code units.
        // Considering the g_1 Dark Energy suppression fixes,
        // the max velocities are roughly ~0.84 max.
        REQUIRE(max_v_dm > 0.3);
        REQUIRE(max_v_dm < 1.3);
    }
}

TEST_CASE("Zeldovich field amplitudes scale perfectly with sigma_8",
          "[ics][cosmology]") {
    Config config;
    config.mesh_size = 32;
    config.box_size_mpc = 100.0;
    config.omega_m = 0.3;
    config.omega_lambda = 0.7;
    config.hubble_h = 0.7;
    config.spectral_index = 0.96;
    config.seed = 42;

    double start_a = 0.02;

    // Generate field with base sigma_8
    config.sigma_8 = 0.8;
    ZeldovichField zf_base = compute_zeldovich_field(start_a, config);

    // Generate field with EXACTLY double the sigma_8
    config.sigma_8 = 1.6;
    ZeldovichField zf_double = compute_zeldovich_field(start_a, config);

    // Verify mathematical linearity cell-by-cell
    int M3 = config.mesh_size * config.mesh_size * config.mesh_size;
    for (int i = 0; i < M3; ++i) {
        REQUIRE(zf_double.dx[i] ==
                Catch::Approx(2.0 * zf_base.dx[i]).margin(1e-12));
        REQUIRE(zf_double.dy[i] ==
                Catch::Approx(2.0 * zf_base.dy[i]).margin(1e-12));
        REQUIRE(zf_double.dz[i] ==
                Catch::Approx(2.0 * zf_base.dz[i]).margin(1e-12));
    }
}

TEST_CASE("Peculiar velocities exactly obey Zel'dovich kinematics",
          "[ics][kinematics]") {
    Config config;
    config.mesh_size = 16;
    config.num_particles_1d = 16;
    config.box_size_mpc = 100.0;
    config.omega_m = 0.3;
    config.omega_baryon = 0.048;
    config.omega_lambda = 0.7;
    config.hubble_h = 0.7;
    config.a_start = 0.02;
    config.sigma_8 = 0.81;
    config.seed = 123;
    config.compute_derived_data();

    // Generate state and grab the Zeldovich field growth rate
    SimState state = initialize_state(config);
    ZeldovichField zf = compute_zeldovich_field(config.a_start, config);

    double H = state.hubble_param;
    double f = zf.f;
    double expected_ratio = H * f;

    // Verify every single particle obeys the velocity-displacement law
    for (size_t i = 0; i < state.dm.num_particles; ++i) {
        // Skip stationary particles at perfectly zero displacement
        if (std::abs(state.dm.vel_x[i]) < 1e-14) continue;

        // Reverse engineer the displacement from the velocity
        double inferred_dx = state.dm.vel_x[i] / expected_ratio;
        double inferred_dy = state.dm.vel_y[i] / expected_ratio;
        double inferred_dz = state.dm.vel_z[i] / expected_ratio;

        // Check that pos = unperturbed_pos + displacement
        // (Taking into account periodic boundary conditions)
        double unperturbed_x =
            std::fmod(state.dm.pos_x[i] - inferred_dx + config.domain_size,
                      config.domain_size);
        double unperturbed_y =
            std::fmod(state.dm.pos_y[i] - inferred_dy + config.domain_size,
                      config.domain_size);
        double unperturbed_z =
            std::fmod(state.dm.pos_z[i] - inferred_dz + config.domain_size,
                      config.domain_size);

        // The unperturbed positions should perfectly align with the regular
        // lattice spacing
        double spacing = config.domain_size / config.num_particles_1d;

        REQUIRE(std::fmod(unperturbed_x, spacing) ==
                Catch::Approx(0.5 * spacing).margin(1e-7));
        REQUIRE(std::fmod(unperturbed_y, spacing) ==
                Catch::Approx(0.5 * spacing).margin(1e-7));
        REQUIRE(std::fmod(unperturbed_z, spacing) ==
                Catch::Approx(0.5 * spacing).margin(1e-7));
    }
}

TEST_CASE("Zeldovich field has exactly zero mean (k=0 mode is dead)",
          "[ics][math]") {
    Config config;
    config.mesh_size = 32;
    config.box_size_mpc = 100.0;
    config.omega_m = 0.3;
    config.omega_lambda = 0.7;
    config.hubble_h = 0.7;
    config.sigma_8 = 0.81;
    config.spectral_index = 0.96;
    config.seed = 42;

    ZeldovichField zf = compute_zeldovich_field(0.02, config);

    double sum_dx = 0.0, sum_dy = 0.0, sum_dz = 0.0;
    int M3 = config.mesh_size * config.mesh_size * config.mesh_size;

    for (int i = 0; i < M3; ++i) {
        sum_dx += zf.dx[i];
        sum_dy += zf.dy[i];
        sum_dz += zf.dz[i];
    }

    // The average displacement of the entire universe must be zero
    REQUIRE(std::abs(sum_dx / M3) < 1e-14);
    REQUIRE(std::abs(sum_dy / M3) < 1e-14);
    REQUIRE(std::abs(sum_dz / M3) < 1e-14);
}

TEST_CASE("standing_particles flag explicitly zeroes all initial velocities",
          "[ics][config]") {
    Config config;
    config.mesh_size = 16;
    config.num_particles_1d = 16;
    config.hydro_method = HydroMethod::Eulerian;

    // Explicitly set to true
    config.standing_particles = true;
    config.compute_derived_data();

    SimState state = initialize_state(config);

    // Check Dark Matter
    for (size_t i = 0; i < state.dm.num_particles; ++i) {
        REQUIRE(state.dm.vel_x[i] == 0.0);
        REQUIRE(state.dm.vel_y[i] == 0.0);
        REQUIRE(state.dm.vel_z[i] == 0.0);
    }

    // Check Gas Momentum
    int N3 = config.mesh_size * config.mesh_size * config.mesh_size;
    for (int i = 0; i < N3; ++i) {
        REQUIRE(state.gas->get_momentum_x().data[i] == 0.0);
        REQUIRE(state.gas->get_momentum_y().data[i] == 0.0);
        REQUIRE(state.gas->get_momentum_z().data[i] == 0.0);
    }
}