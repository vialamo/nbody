#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "config.h"
#include "diagnostics.h"
#include "gravity.h"
#include "particles.h"
#include "pm.h"
#include "state.h"

// ---------------------------------------------------------------------
// Cloud-in-Cell (CIC) Density Assignment
// ---------------------------------------------------------------------
TEST_CASE("PM Density Assignment maps particles to grid correctly",
          "[pm][density]") {
    Config config;
    config.mesh_size_1d = 4;
    config.hydro_method = HydroMethod::None;
    config.compute_derived_data(); // Sets cell_size to 0.25 and volume to 0.015625

    auto get_total_grid_mass = [&](const SimState& s) {
        double total_mass = 0.0;
        // Evaluate the master density grid used by the Poisson solver
        for (double rho_val : s.pm.rho.data) {
            total_mass += rho_val * config.cell_volume;
        }
        return total_mass;
    };

    SECTION("Particle exactly at a cell center (grid node, 0 fractional part)") {
        SimState state(config);

        // Node is at (i + 0.5) * cell_size
        // x = 1.5 * 0.25 = 0.375 | y = 2.5 * 0.25 = 0.625 | z = 3.5 * 0.25 = 0.875
        state.dm.add_particle(0.375, 0.625, 0.875, 0.0, 0.0, 0.0, 16.0);
        
        // Use the unified PM orchestrator instead of particle-level methods
        state.pm.bin_and_assign_mass(state, config);

        REQUIRE(get_total_grid_mass(state) == Catch::Approx(16.0));
        double expected_density = 16.0 / config.cell_volume;
        REQUIRE(state.pm.rho(1, 2, 3) == Catch::Approx(expected_density));
        REQUIRE(state.pm.rho(2, 2, 3) == Catch::Approx(0.0));
    }

    SECTION("Particle on a cell boundary splits perfectly into 8 grid nodes") {
        SimState state(config);

        // Boundary intersection of 8 cells is at integer multiples of cell_size
        state.dm.add_particle(0.5, 0.5, 0.5, 0.0, 0.0, 0.0, 8.0);
        state.pm.bin_and_assign_mass(state, config);

        REQUIRE(get_total_grid_mass(state) == Catch::Approx(8.0));
        double expected_density = 1.0 / config.cell_volume;

        REQUIRE(state.pm.rho(1, 1, 1) == Catch::Approx(expected_density));
        REQUIRE(state.pm.rho(2, 1, 1) == Catch::Approx(expected_density));
        REQUIRE(state.pm.rho(1, 2, 1) == Catch::Approx(expected_density));
        REQUIRE(state.pm.rho(2, 2, 1) == Catch::Approx(expected_density));
        REQUIRE(state.pm.rho(1, 1, 2) == Catch::Approx(expected_density));
        REQUIRE(state.pm.rho(2, 1, 2) == Catch::Approx(expected_density));
        REQUIRE(state.pm.rho(1, 2, 2) == Catch::Approx(expected_density));
        REQUIRE(state.pm.rho(2, 2, 2) == Catch::Approx(expected_density));
    }

    SECTION("Particle on the periodic boundary wraps to the opposite side") {
        SimState state(config);

        // X=0.0 is the domain boundary. Y and Z are nodes (i=0 -> 0.125)
        state.dm.add_particle(0.0, 0.125, 0.125, 0.0, 0.0, 0.0, 10.0);
        state.pm.bin_and_assign_mass(state, config);

        REQUIRE(get_total_grid_mass(state) == Catch::Approx(10.0));
        double expected_density = 5.0 / config.cell_volume;

        REQUIRE(state.pm.rho(3, 0, 0) == Catch::Approx(expected_density));
        REQUIRE(state.pm.rho(0, 0, 0) == Catch::Approx(expected_density));
    }
}

// ---------------------------------------------------------------------
// Short-Range Gravity (PP / P3M)
// ---------------------------------------------------------------------
TEST_CASE("Short-range gravity calculates Newtonian and P3M forces",
          "[gravity][pp]") {
    Config config;
    config.mesh_size_1d = 4;
    config.hydro_method = HydroMethod::None;
    config.compute_derived_data();
    
    config.G = 1.0;
    config.softening = 0.01;
    config.base_comoving_softening = 0.01;

    Diagnostics dummy_diag;
    SimState state(config);

    // Two particles 0.25 units apart on the X-axis
    state.dm.add_particle(0.25, 0.5, 0.5, 0.0, 0.0, 0.0, 1.0);
    state.dm.add_particle(0.50, 0.5, 0.5, 0.0, 0.0, 0.0, 1.0);

    SECTION("Pure Newtonian PP Force (use_PM = false)") {
        config.use_PM = false;
        config.cutoff_radius_squared = std::numeric_limits<double>::infinity();

        state.dm.build_lbvh(config);
        
        state.gravity.add_pp_accel(state.dm, nullptr, 1.0, config, dummy_diag);

        // Exact Newtonian Gravity: F = G * m1 * m2 / r^2
        // F = 1.0 * 1.0 * 1.0 / (0.25 * 0.25) = 1.0 / 0.0625 = 16.0
        double expected_f = 16.0;

        // Convert to acceleration (a = F/m). Since m=1, a=F.
        REQUIRE(state.dm.acc_x[0] == Catch::Approx(expected_f).margin(1e-1));
        REQUIRE(state.dm.acc_y[0] == Catch::Approx(0.0));
        REQUIRE(state.dm.acc_z[0] == Catch::Approx(0.0));

        // p2 is pulled to the left (-x) by p1
        REQUIRE(state.dm.acc_x[1] == Catch::Approx(-expected_f).margin(1e-1));
    }

    SECTION("P3M Cutoff safely ignores distant particles") {
        config.use_PM = true;
        config.cutoff_radius = 0.15;
        config.cutoff_radius_squared = 0.0225;

        state.dm.build_lbvh(config);
        
        // Reset accelerations to 0 before computing
        std::fill(state.dm.acc_x.begin(), state.dm.acc_x.end(), 0.0);

        state.gravity.add_pp_accel(state.dm, nullptr, 1.0, config, dummy_diag);

        // Because dist (0.25) > cutoff_radius (0.15), the short-range force is 0.0
        REQUIRE(state.dm.acc_x[0] == Catch::Approx(0.0));
        REQUIRE(state.dm.acc_x[1] == Catch::Approx(0.0));
    }
}

// ---------------------------------------------------------------------
// Long-Range Gravity (PM Poisson Solver & Force Interpolation)
// ---------------------------------------------------------------------
TEST_CASE("PM Poisson Solver calculates accurate grid accelerations",
          "[gravity][pm]") {
    Config config;
    config.mesh_size_1d = 4;
    config.hydro_method = HydroMethod::None;
    config.compute_derived_data();
    config.G = 1.0;

    SimState state(config);

    SECTION("A uniform density field produces zero net acceleration") {
        state.pm.rho.fill(5.0);

        state.pm.compute_acceleration(state, config);

        REQUIRE(state.pm.acc_x.maxCoeff() == Catch::Approx(0.0).margin(1e-10));
        REQUIRE(state.pm.acc_x.minCoeff() == Catch::Approx(0.0).margin(1e-10));
        REQUIRE(state.pm.acc_y.maxCoeff() == Catch::Approx(0.0).margin(1e-10));
    }

    SECTION("A single point mass creates symmetric attractive forces") {
        // Place a massive particle on the cell-centered PM node (1, 1, 1)
        state.dm.add_particle(0.375, 0.375, 0.375, 0.0, 0.0, 0.0, 100.0);

        // Maps particles to the rho grid
        state.pm.bin_and_assign_mass(state, config);

        state.pm.compute_acceleration(state, config);

        // Center of mass check (Should be 0 due to central differences)
        REQUIRE(state.pm.acc_x(1, 1, 1) == Catch::Approx(0.0).margin(1e-10));

        // Attraction checks (The cell to the right (+x) should be pulled to the left (-x))
        REQUIRE(state.pm.acc_x(2, 1, 1) < 0.0);

        // Attraction checks (The cell to the left (-x) should be pulled to the right (+x))
        REQUIRE(state.pm.acc_x(0, 1, 1) > 0.0);

        // Magnitude symmetry check
        REQUIRE(state.pm.acc_x(0, 1, 1) == Catch::Approx(-state.pm.acc_x(2, 1, 1)));
    }
}
