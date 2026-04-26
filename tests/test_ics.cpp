#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "config.h"
#include "ics.h"

TEST_CASE("Initial internal energy scales correctly with code units",
          "[hydro][thermodynamics][units]") {
    // Baseline configuration parameters
    double T_kelvin = 50.0;
    double h_param = 0.7;
    double box_mpc = 100.0;
    double domain = 1.0;
    double gamma = 5.0 / 3.0;

    SECTION("Changing Domain Size dynamically updates the velocity unit") {
        // Baseline: Domain of 1.0
        double u_code_base = get_internal_energy_from_temp_k(
            T_kelvin, h_param, box_mpc, 1.0, gamma);

        // Domain of 2.0 means 1 code unit = 50 Mpc instead of 100 Mpc.
        // This halves V_unit. Since u_code is divided by V_unit^2,
        // halving V_unit should exactly quadruple the final u_code result.
        double u_code_double_domain = get_internal_energy_from_temp_k(
            T_kelvin, h_param, box_mpc, 2.0, gamma);

        REQUIRE(u_code_double_domain == Catch::Approx(u_code_base * 4.0));
    }

    SECTION("Expected mathematical output is physically accurate") {
        double u_code = get_internal_energy_from_temp_k(T_kelvin, h_param,
                                                        box_mpc, domain, gamma);

        // Manual calculation check:
        // u_phys = (1.380649e-23 * 50) / ((2/3) * 1.22 * 1.67262192e-27) ≈
        // 507442.7 m^2/s^2 V_unit = 1.5 * 70.0 * 100.0 = 10500 km/s ->
        // 10,500,000 m/s V_unit^2 = 1.1025e14 u_code = 507442.7 / 1.1025e14
        // ≈ 4.60265e-9

        REQUIRE(u_code == Catch::Approx(4.60265e-9).margin(1e-13));
    }
}

TEST_CASE("initialize_state produces physically sound macro-states",
          "[ics][system]") {
    // 1. Manually build the configuration
    Config config;
    config.DOMAIN_SIZE = 1.0;
    config.MESH_SIZE = 32;
    config.BOX_SIZE_MPC = 100.0;
    config.N_PER_SIDE = 32;
    config.USE_PM = true;
    config.USE_PP = false;

    config.compute_derived_data();

    // 2. Initialize the entire simulation state
    SimState state = initialize_state(config);

    SECTION("Conservation of Mass") {
        double total_dm_mass =
            state.dm.particles.size() * config.DM_PARTICLE_MASS;

        double total_gas_mass = 0.0;
        int N3 = config.MESH_SIZE * config.MESH_SIZE * config.MESH_SIZE;
        for (int i = 0; i < N3; ++i) {
            // Integrate density over Eulerian volume
            total_gas_mass +=
                state.gas.get_density().data[i] * config.CELL_VOLUME;
        }

        // The universe must have exactly 1.0 code mass
        REQUIRE(total_dm_mass + total_gas_mass ==
                Catch::Approx(1.0).margin(1e-6));
        REQUIRE(total_gas_mass ==
                Catch::Approx(config.GAS_TOTAL_MASS).margin(1e-6));
    }

    SECTION("Conservation of Momentum") {
        double p_dm_x = 0.0;
        for (const auto &p : state.dm.particles) {
            p_dm_x += p.vel.x * p.mass;
        }

        double p_gas_x_density = 0.0;
        int N3 = config.MESH_SIZE * config.MESH_SIZE * config.MESH_SIZE;
        for (int i = 0; i < N3; ++i) {
            p_gas_x_density += state.gas.get_momentum_x().data[i];
        }

        // Convert momentum density to actual physical momentum
        double physical_p_gas_x = p_gas_x_density * config.CELL_VOLUME;

        // In a periodic box, the sum of all momenta must be zero
        REQUIRE(std::abs(p_dm_x) < 1e-10);
        REQUIRE(std::abs(physical_p_gas_x) < 1e-3);
    }

    SECTION("Velocities are mathematically bound to cosmic expansion") {
        double max_v_dm = 0.0;
        for (const auto &p : state.dm.particles) {
            max_v_dm = std::max(max_v_dm, std::abs(p.vel.x));
        }

        // At a=0.02, H(a) is ~129 code units.
        // Max scaled displacement is ~0.003 code units.
        // Therefore, absolute max velocity must be near ~0.38
        REQUIRE(max_v_dm > 0.1);
        REQUIRE(max_v_dm < 0.6);
    }
}

TEST_CASE("Zeldovich field amplitudes scale perfectly with sigma_8",
          "[ics][cosmology]") {
    Config config;
    config.DOMAIN_SIZE = 1.0;
    config.MESH_SIZE = 32;
    config.BOX_SIZE_MPC = 100.0;
    config.OMEGA_M = 0.3;
    config.OMEGA_LAMBDA = 0.7;
    config.HUBBLE_PARAM = 0.7;
    config.SPECTRAL_INDEX = 0.96;
    config.SEED = 42;

    double start_a = 0.02;

    // Generate field with base sigma_8
    config.SIGMA_8 = 0.8;
    ZeldovichField zf_base = compute_zeldovich_field(start_a, config);

    // Generate field with EXACTLY double the sigma_8
    config.SIGMA_8 = 1.6;
    ZeldovichField zf_double = compute_zeldovich_field(start_a, config);

    // Verify mathematical linearity cell-by-cell
    int M3 = config.MESH_SIZE * config.MESH_SIZE * config.MESH_SIZE;
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
    config.DOMAIN_SIZE = 1.0;
    config.MESH_SIZE = 16;
    config.N_PER_SIDE = 16;
    config.BOX_SIZE_MPC = 100.0;
    config.OMEGA_M = 0.3;
    config.OMEGA_BARYON = 0.048;
    config.OMEGA_LAMBDA = 0.7;
    config.HUBBLE_PARAM = 0.7;
    config.START_A = 0.02;
    config.SIGMA_8 = 0.81;
    config.SEED = 123;
    config.compute_derived_data();

    // Generate state and grab the Zeldovich field growth rate
    SimState state = initialize_state(config);
    ZeldovichField zf = compute_zeldovich_field(config.START_A, config);

    double H = state.hubble_param;
    double f = zf.f;
    double expected_ratio = H * f;

    // Verify every single particle obeys the velocity-displacement law
    for (const auto &p : state.dm.particles) {
        // Skip stationary particles at perfectly zero displacement
        if (std::abs(p.vel.x) < 1e-14) continue;

        // Reverse engineer the displacement from the velocity
        double inferred_dx = p.vel.x / expected_ratio;
        double inferred_dy = p.vel.y / expected_ratio;
        double inferred_dz = p.vel.z / expected_ratio;

        // Check that pos = unperturbed_pos + displacement
        // (Taking into account periodic boundary conditions)
        double unperturbed_x = std::fmod(
            p.pos.x - inferred_dx + config.DOMAIN_SIZE, config.DOMAIN_SIZE);
        double unperturbed_y = std::fmod(
            p.pos.y - inferred_dy + config.DOMAIN_SIZE, config.DOMAIN_SIZE);
        double unperturbed_z = std::fmod(
            p.pos.z - inferred_dz + config.DOMAIN_SIZE, config.DOMAIN_SIZE);

        // The unperturbed positions should perfectly align with the regular
        // lattice spacing
        double spacing = config.DOMAIN_SIZE / config.N_PER_SIDE;

        REQUIRE(std::fmod(unperturbed_x, spacing) ==
                Catch::Approx(0.5 * spacing).margin(1e-7));
        REQUIRE(std::fmod(unperturbed_y, spacing) ==
                Catch::Approx(0.5 * spacing).margin(1e-7));
        REQUIRE(std::fmod(unperturbed_z, spacing) ==
                Catch::Approx(0.5 * spacing).margin(1e-7));
    }
}

TEST_CASE("Gas density perfectly traces Dark Matter density", "[ics][hydro]") {
    Config config;
    config.DOMAIN_SIZE = 1.0;
    config.MESH_SIZE = 16;
    config.N_PER_SIDE = 16;
    config.BOX_SIZE_MPC = 100.0;
    config.OMEGA_M = 0.3;
    config.OMEGA_BARYON = 0.048;
    config.OMEGA_LAMBDA = 0.7;
    config.HUBBLE_PARAM = 0.7;
    config.START_A = 0.02;
    config.SIGMA_8 = 0.81;
    config.SPECTRAL_INDEX = 0.96;
    config.SEED = 42;
    config.compute_derived_data();

    SimState state = initialize_state(config);

    double total_dm_mass = state.dm.particles.size() * config.DM_PARTICLE_MASS;
    double expected_mass_ratio = config.GAS_TOTAL_MASS / total_dm_mass;

    int N3 = config.MESH_SIZE * config.MESH_SIZE * config.MESH_SIZE;

    for (int i = 0; i < N3; ++i) {
        double dm_rho = state.dm.get_rho().data[i];
        double gas_rho = state.gas.get_density().data[i];

        double expected_gas_rho = dm_rho * expected_mass_ratio;

        if (expected_gas_rho >= 1e-12) {
            REQUIRE(gas_rho == Catch::Approx(expected_gas_rho).margin(1e-10));
        } else {
            // Floor value check
            REQUIRE(gas_rho == Catch::Approx(1e-12));
        }
    }
}

TEST_CASE("Zeldovich field has exactly zero mean (k=0 mode is dead)",
          "[ics][math]") {
    Config config;
    config.DOMAIN_SIZE = 1.0;
    config.MESH_SIZE = 32;
    config.BOX_SIZE_MPC = 100.0;
    config.OMEGA_M = 0.3;
    config.OMEGA_LAMBDA = 0.7;
    config.HUBBLE_PARAM = 0.7;
    config.SIGMA_8 = 0.81;
    config.SPECTRAL_INDEX = 0.96;
    config.SEED = 42;

    ZeldovichField zf = compute_zeldovich_field(0.02, config);

    double sum_dx = 0.0, sum_dy = 0.0, sum_dz = 0.0;
    int M3 = config.MESH_SIZE * config.MESH_SIZE * config.MESH_SIZE;

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

TEST_CASE("STANDING_PARTICLES flag explicitly zeroes all initial velocities",
          "[ics][config]") {
    Config config;
    config.DOMAIN_SIZE = 1.0;
    config.MESH_SIZE = 16;
    config.N_PER_SIDE = 16;
    config.USE_HYDRO = true;

    // Explicitly set to true
    config.STANDING_PARTICLES = true;
    config.compute_derived_data();

    SimState state = initialize_state(config);

    // 1. Check Dark Matter
    for (const auto &p : state.dm.particles) {
        REQUIRE(p.vel.x == 0.0);
        REQUIRE(p.vel.y == 0.0);
        REQUIRE(p.vel.z == 0.0);
    }

    // 2. Check Gas Momentum
    int N3 = config.MESH_SIZE * config.MESH_SIZE * config.MESH_SIZE;
    for (int i = 0; i < N3; ++i) {
        REQUIRE(state.gas.get_momentum_x().data[i] == 0.0);
        REQUIRE(state.gas.get_momentum_y().data[i] == 0.0);
        REQUIRE(state.gas.get_momentum_z().data[i] == 0.0);
    }
}