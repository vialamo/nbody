#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "config.h"
#include "cooling.h"
#include "thermodynamics.h"

TEST_CASE("Thermodynamic conversions are perfectly reversible",
          "[cooling][thermodynamics]") {
    Config config;
    config.compute_derived_data();

    double original_T = 1.5e5;  // 150,000 K
    double a = 1.0;

    // Convert to code units and back
    double u_code =
        Thermodynamics::get_internal_energy_from_temp(original_T, a, config);
    double recovered_T =
        Thermodynamics::get_temp_from_internal_energy(u_code, a, config);

    REQUIRE(recovered_T == Catch::Approx(original_T));
}

TEST_CASE("Cooling rate behaves physically", "[cooling][physics]") {
    Config config;
    config.compute_derived_data();
    double a = 1.0;
    double Z = 0.0;
    double rho = 1.0;
    double dt = 1e10;  // Large timestep to force significant cooling
    int iterations = 0;
    Cooling cooling(config);

    SECTION("compute_du_dt acts as a continuous mathematical curve") {
        double cold_temp = 50.0;
        double u_cold =
            Thermodynamics::get_internal_energy_from_temp(cold_temp, a, config);

        double du_dt = cooling.compute_du_dt(u_cold, rho, Z, a, config);

        // It must NOT be zero. The solver relies on this function
        // returning a valid derivative everywhere.
        REQUIRE(du_dt < 0.0);
    }

    SECTION(
        "Implicit solver pre-check blocks cooling for cold primordial gas") {
        double cold_temp = config.cooling_cutoff_k;
        double u_cold =
            Thermodynamics::get_internal_energy_from_temp(cold_temp, a, config);

        double u_new = cooling.solve_cooling_implicit(
            u_cold, rho, Z, a, dt, u_cold, u_cold, config, iterations);

        // The solver should return immediately without altering the energy
        REQUIRE(u_new == u_cold);
        REQUIRE(iterations == 0);
    }

    SECTION(
        "Implicit solver clamps aggressively cooling gas to the physical "
        "floor") {
        double expected_floor_k = std::max(config.cooling_cutoff_k, config.temp_floor_k);
        double hot_temp = expected_floor_k + 10.0;  // Just above the floor
        double u_hot =
            Thermodynamics::get_internal_energy_from_temp(hot_temp, a, config);

        // Target physical floor defined in the solver
        double target_u_floor =
            Thermodynamics::get_internal_energy_from_temp(expected_floor_k, a, config);

        // With a massive dt, it will attempt to cool well below floor
        double massive_dt = 1e15;
        double u_new = cooling.solve_cooling_implicit(
            u_hot, rho, Z, a, massive_dt, target_u_floor, target_u_floor, config, iterations);

        // The solver must intercept the overshoot and clamp exactly to the
        // threshold
        REQUIRE(u_new == Catch::Approx(target_u_floor).epsilon(1e-5));
    }

    SECTION("Cooling scales with density squared (Bremsstrahlung)") {
        double u_hot = Thermodynamics::get_internal_energy_from_temp(1e6, a, config);

        double du_dt_base = cooling.compute_du_dt(u_hot, rho, Z, a, config);
        double du_dt_double_rho =
            cooling.compute_du_dt(u_hot, 2 * rho, Z, a, config);

        // Since du_dt is specific energy change (du/dt), and volumetric change
        // scales as rho^2, specific change du/dt scales linearly with rho.
        REQUIRE(du_dt_double_rho == Catch::Approx(du_dt_base * 2.0));
    }
}

TEST_CASE("Implicit backward Euler solver is unconditionally stable",
          "[cooling][solver]") {
    Config config;
    config.temp_floor_k = 10.0;
    config.compute_derived_data();
    Cooling cooling(config);

    double u_initial =
        Thermodynamics::get_internal_energy_from_temp(1e7, 1.0, config);  // Very hot
    double rho = 100.0;                                            // Dense
    double a = 1.0;
    double Z = 0.0;
    int iters = 0;

    SECTION("A massive timestep drops to the floor but never goes negative") {
        double dt_massive = 1e20;  // Practically infinite time

        double expected_floor_k = std::max(config.cooling_cutoff_k, config.temp_floor_k);
        double u_floor =
            Thermodynamics::get_internal_energy_from_temp(expected_floor_k, a, config);
        double u_final = cooling.solve_cooling_implicit(
            u_initial, rho, Z, a, dt_massive, u_floor, u_floor, config, iters);
        double T_final =
            Thermodynamics::get_temp_from_internal_energy(u_final, a, config);

        REQUIRE(u_final > 0.0);
        // It should have safely caught itself at the floor
        REQUIRE(T_final == Catch::Approx(expected_floor_k).margin(1e-5));
    }

    SECTION("A normal timestep cools the gas a physical amount") {
        double dt_normal = 0.001;
        double u_final = cooling.solve_cooling_implicit(
            u_initial, rho, Z, a, dt_normal, u_initial * 0.1, u_initial * 0.1, config, iters);

        REQUIRE(u_final < u_initial);
        REQUIRE(u_final > 0.0);
    }
}

