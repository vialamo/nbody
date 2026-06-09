#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "config.h"
#include "cooling.h"
#include "gas.h"

// Helper proxy to inject values into the GasGrid
struct GasGridTestAccess {
    static Grid3D& density(GasGrid& g) { return g.density; }
    static Grid3D& momentum_x(GasGrid& g) { return g.momentum_x; }
    static Grid3D& energy(GasGrid& g) { return g.energy; }
    static Grid3D& internal_energy(GasGrid& g) { return g.internal_energy; }
};

TEST_CASE("Thermodynamic conversions are perfectly reversible",
          "[cooling][thermodynamics]") {
    Config config;
    config.compute_derived_data();

    double original_T = 1.5e5;  // 150,000 K
    double a = 1.0;

    // Convert to code units and back
    double u_code = cooling::get_internal_energy_from_temp(original_T, a, config);
    double recovered_T = cooling::get_temp_from_internal_energy(u_code, a, config);

    REQUIRE(recovered_T == Catch::Approx(original_T));
}

TEST_CASE("Cooling rate behaves physically", "[cooling][physics]") {
    Config config;
    config.compute_derived_data();
    double a = 1.0;
    double rho = 1.0;
    double dt = 1e10; // Large timestep to force significant cooling
    int iterations = 0;

    SECTION("compute_cooling_rate acts as a continuous mathematical curve") {
        double cold_temp = 50.0; 
        double u_cold = cooling::get_internal_energy_from_temp(cold_temp, a, config);
        
        double lambda = cooling::compute_cooling_rate(u_cold, rho, a, config);

        // It must NOT be zero. The solver relies on this function 
        // returning a valid derivative everywhere.
        REQUIRE(lambda > 0.0);
    }

    SECTION("Implicit solver pre-check blocks cooling for cold primordial gas") {
        double cold_temp = 50.0; // Well below the 10,000 K physical floor
        double u_cold = cooling::get_internal_energy_from_temp(cold_temp, a, config);
        
        double u_new = cooling::solve_cooling_implicit(u_cold, rho, a, dt, config, iterations);

        // The solver should return immediately without altering the energy
        REQUIRE(u_new == u_cold);
        REQUIRE(iterations == 0);
    }

    SECTION("Implicit solver clamps aggressively cooling gas to the physical floor") {
        double hot_temp = 10500.0; // Just above the floor
        double u_hot = cooling::get_internal_energy_from_temp(hot_temp, a, config);
        
        // Target physical floor defined in the solver
        double expected_floor_k = std::max(10000.0, config.TEMP_FLOOR_KELVIN);
        double target_u_floor = cooling::get_internal_energy_from_temp(expected_floor_k, a, config);

        // With a massive dt, it will attempt to cool well below 10,000 K
        double massive_dt = 1e15;
        double u_new = cooling::solve_cooling_implicit(u_hot, rho, a, massive_dt, config, iterations);

        // The solver must intercept the overshoot and clamp exactly to the threshold
        REQUIRE(u_new == Catch::Approx(target_u_floor).epsilon(1e-5));
    }

    SECTION("Cooling scales with density squared (Bremsstrahlung)") {
        double u_hot = cooling::get_internal_energy_from_temp(1e6, a, config);

        double lambda_base = cooling::compute_cooling_rate(u_hot, rho, a, config);
        double lambda_double_rho =
            cooling::compute_cooling_rate(u_hot, 2 * rho, a, config);

        // Since lambda is specific energy loss (du/dt), and volumetric loss
        // scales as rho^2, specific loss du/dt scales linearly with rho.
        REQUIRE(lambda_double_rho == Catch::Approx(lambda_base * 2.0));
    }
}

TEST_CASE("Implicit backward Euler solver is unconditionally stable",
          "[cooling][solver]") {
    Config config;
    config.TEMP_FLOOR_KELVIN = 10.0;
    config.compute_derived_data();

    double u_initial =
        cooling::get_internal_energy_from_temp(1e7, 1.0, config);  // Very hot
    double rho = 100.0;                                       // Dense
    double a = 1.0;
    int iters = 0;

    SECTION("A massive timestep drops to the floor but never goes negative") {
        double dt_massive = 1e20;  // Practically infinite time

        double u_final = cooling::solve_cooling_implicit(
            u_initial, rho, a, dt_massive, config, iters);
        double T_final =
            cooling::get_temp_from_internal_energy(u_final, a, config);

        double expected_floor_k = std::max(10000.0, config.TEMP_FLOOR_KELVIN);

        REQUIRE(u_final > 0.0);
        // It should have safely caught itself at the floor
        REQUIRE(T_final == Catch::Approx(expected_floor_k).margin(1e-5));
    }

    SECTION("A normal timestep cools the gas a physical amount") {
        double dt_normal = 0.001;
        double u_final = cooling::solve_cooling_implicit(
            u_initial, rho, a, dt_normal, config, iters);

        REQUIRE(u_final < u_initial);
        REQUIRE(u_final > 0.0);
    }
}

TEST_CASE("GasGrid correctly extracts radiated energy from BOTH arrays",
          "[cooling][gas]") {
    Config config;
    config.MESH_SIZE = 2;
    config.ENABLE_COOLING = true;
    config.compute_derived_data();
    double a = 1.0;

    GasGrid grid(config);

    auto& rho = GasGridTestAccess::density(grid);
    auto& mx = GasGridTestAccess::momentum_x(grid);
    auto& en = GasGridTestAccess::energy(grid);
    auto& ie = GasGridTestAccess::internal_energy(grid);

    // Setup a hot, dense cell
    rho(0, 0, 0) = 10.0;

    // Give it 50.0 Kinetic Energy (KE = 0.5 * rho * v^2)
    // 50.0 = 0.5 * 10.0 * v^2  =>  v = sqrt(10)
    // momentum = rho * v = 10 * sqrt(10)
    mx(0, 0, 0) = 10.0 * std::sqrt(10.0);

    // Give it 100.0 total energy, 50.0 of which is internal, 50.0 is kinetic
    en(0, 0, 0) = 100.0;
    ie(0, 0, 0) = 50.0;

    // Apply cooling
    double dt = 0.05;
    double e_lost = grid.apply_cooling(dt, a);

    // Verify energy was lost
    REQUIRE(e_lost > 0.0);
    REQUIRE(ie(0, 0, 0) < 50.0);
    REQUIRE(en(0, 0, 0) < 100.0);

    // Did total energy and internal energy drop by the EXACT SAME amount?
    // This proves kinetic energy was perfectly preserved.
    double ie_drop = 50.0 - ie(0, 0, 0);
    double en_drop = 100.0 - en(0, 0, 0);

    REQUIRE(ie_drop == Catch::Approx(en_drop));
}