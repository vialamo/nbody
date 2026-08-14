#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "config.h"
#include "ics.h"
#include "integrator.h"
#include "particles.h"
#include "state.h"

// ---------------------------------------------------------------------
// Cloud-in-Cell (CIC) Density Assignment
// ---------------------------------------------------------------------
TEST_CASE("CIC Density Assignment maps particles to grid correctly",
          "[cic][particles]") {
    Config config;
    config.mesh_size = 4;
    config.compute_derived_data(); // Automatically sets cell_size to 0.25 and volume to 0.015625
    config.num_dm_particles = 1;

    auto get_total_grid_mass = [&](const ParticleSystem& sys) {
        double total_mass = 0.0;
        for (int i = 0;
             i < config.mesh_size * config.mesh_size * config.mesh_size; ++i) {
            total_mass += sys.get_rho().data[i] * config.cell_volume;
        }
        return total_mass;
    };

    SECTION("Particle exactly at a cell center (grid node, 0 fractional part)") {
        ParticleSystem sys(config);

        // Node is at (i + 0.5) * cell_size
        // x = 1.5 * 0.25 = 0.375 | y = 2.5 * 0.25 = 0.625 | z = 3.5 * 0.25 = 0.875
        sys.add_particle(0.375, 0.625, 0.875, 0.0, 0.0, 0.0, 16.0);
        sys.bin_and_assign_mass(config);

        REQUIRE(get_total_grid_mass(sys) == Catch::Approx(16.0));
        double expected_density = 16.0 / config.cell_volume;
        REQUIRE(sys.get_rho()(1, 2, 3) == Catch::Approx(expected_density));
        REQUIRE(sys.get_rho()(2, 2, 3) == Catch::Approx(0.0));
    }

    SECTION("Particle on a cell boundary splits perfectly into 8 grid nodes") {
        ParticleSystem sys(config);

        // Boundary intersection of 8 cells is at integer multiples of cell_size 
        // e.g., (2 * 0.25) = 0.5
        sys.add_particle(0.5, 0.5, 0.5, 0.0, 0.0, 0.0, 8.0);
        sys.bin_and_assign_mass(config);

        REQUIRE(get_total_grid_mass(sys) == Catch::Approx(8.0));
        double expected_density = 1.0 / config.cell_volume;

        REQUIRE(sys.get_rho()(1, 1, 1) == Catch::Approx(expected_density));
        REQUIRE(sys.get_rho()(2, 1, 1) == Catch::Approx(expected_density));
        REQUIRE(sys.get_rho()(1, 2, 1) == Catch::Approx(expected_density));
        REQUIRE(sys.get_rho()(2, 2, 1) == Catch::Approx(expected_density));
        REQUIRE(sys.get_rho()(1, 1, 2) == Catch::Approx(expected_density));
        REQUIRE(sys.get_rho()(2, 1, 2) == Catch::Approx(expected_density));
        REQUIRE(sys.get_rho()(1, 2, 2) == Catch::Approx(expected_density));
        REQUIRE(sys.get_rho()(2, 2, 2) == Catch::Approx(expected_density));
    }

    SECTION("Particle on the periodic boundary wraps to the opposite side") {
        ParticleSystem sys(config);

        // X=0.0 is the domain boundary. Y and Z are exact nodes (i=0 -> 0.125)
        sys.add_particle(0.0, 0.125, 0.125, 0.0, 0.0, 0.0, 10.0);
        sys.bin_and_assign_mass(config);

        REQUIRE(get_total_grid_mass(sys) == Catch::Approx(10.0));
        double expected_density = 5.0 / config.cell_volume;

        REQUIRE(sys.get_rho()(3, 0, 0) == Catch::Approx(expected_density));
        REQUIRE(sys.get_rho()(0, 0, 0) == Catch::Approx(expected_density));
    }
}

// ---------------------------------------------------------------------
// Short-Range Gravity (PP / P3M)
// ---------------------------------------------------------------------
TEST_CASE("Short-range gravity calculates Newtonian and P3M forces",
          "[gravity][pp]") {
    Config config;
    config.mesh_size = 4;
    config.cell_size = 0.25;
    config.num_dm_particles = 2;

    config.G = 1.0;
    config.softening_squared = 0.0001;

    Diagnostics dummy_diag;

    ParticleSystem sys(config);
    // Two particles exactly 0.25 units apart on the X-axis
    sys.add_particle(0.25, 0.5, 0.5, 0.0, 0.0, 0.0, 1.0);
    sys.add_particle(0.50, 0.5, 0.5, 0.0, 0.0, 0.0, 1.0);

    SECTION("Pure Newtonian PP Force (use_PM = false)") {
        config.use_PM = false;
        config.cutoff_radius_squared =
            (config.domain_size / 2.0) * (config.domain_size / 2.0);

        sys.bin_and_assign_mass(config);
        sys.compute_and_add_pp_forces(config, dummy_diag);

        // Exact Newtonian Gravity: F = G * m1 * m2 / r^2
        // F = 1.0 * 1.0 * 1.0 / (0.25 * 0.25) = 1.0 / 0.0625 = 16.0
        double expected_f = 16.0;

        // Convert to acceleration (a = F/m). Since m=1, a=F.
        // p1 is pulled to the right (+x) by p2
        REQUIRE(sys.acc_x[0] == Catch::Approx(expected_f).margin(1e-1));
        REQUIRE(sys.acc_y[0] == Catch::Approx(0.0));
        REQUIRE(sys.acc_z[0] == Catch::Approx(0.0));

        // p2 is pulled to the left (-x) by p1
        REQUIRE(sys.acc_x[1] == Catch::Approx(-expected_f).margin(1e-1));
    }

    SECTION("P3M Cutoff safely ignores distant particles") {
        config.use_PM = true;
        config.cutoff_radius = 0.15;
        config.cutoff_radius_squared = 0.0225;

        sys.bin_and_assign_mass(config);

        // Reset accelerations to 0 before computing
        std::fill(sys.acc_x.begin(), sys.acc_x.end(), 0.0);

        sys.compute_and_add_pp_forces(config, dummy_diag);

        // Because dist (0.25) > cutoff_radius (0.15), the force should be
        // exactly 0.0
        REQUIRE(sys.acc_x[0] == Catch::Approx(0.0));
        REQUIRE(sys.acc_x[1] == Catch::Approx(0.0));
    }
}

// ---------------------------------------------------------------------
// Long-Range Gravity (PM Poisson Solver & CIC Interpolation)
// ---------------------------------------------------------------------
TEST_CASE("PM Poisson Solver calculates accurate grid accelerations",
          "[gravity][pm]") {
    Config config;
    config.mesh_size = 4;
    config.cell_size = 0.25;
    config.num_dm_particles = 1;
    config.G = 1.0;
    config.hydro_method = HydroMethod::None;

    SimState state(config);

    SECTION("A uniform density field produces zero net acceleration") {
        state.dm.dm_rho.fill(5.0);

        // Compute gravity (it will copy dm_rho into total_rho internally)
        compute_PM_acceleration(state, config);

        REQUIRE(state.pm_gravity_x.maxCoeff() ==
                Catch::Approx(0.0).margin(1e-10));
        REQUIRE(state.pm_gravity_x.minCoeff() ==
                Catch::Approx(0.0).margin(1e-10));
        REQUIRE(state.pm_gravity_y.maxCoeff() ==
                Catch::Approx(0.0).margin(1e-10));
    }

    SECTION("A single point mass creates symmetric attractive forces") {
        // Place a massive particle exactly on the cell-centered PM node (1, 1,
        // 1) x = (1 + 0.5) * 0.25 = 0.375
        state.dm.add_particle(0.375, 0.375, 0.375, 0.0, 0.0, 0.0, 100.0);

        // Run the binning algorithm to populate dm_rho
        state.dm.bin_and_assign_mass(config);

        // Now compute gravity. It will successfully copy the 100.0 spike
        compute_PM_acceleration(state, config);

        // Center of mass check (Should be exactly 0 due to central differences)
        REQUIRE(state.pm_gravity_x(1, 1, 1) ==
                Catch::Approx(0.0).margin(1e-10));

        // Attraction checks (The cell to the right (+x) should be pulled to the
        // left (-x))
        REQUIRE(state.pm_gravity_x(2, 1, 1) < 0.0);

        // Attraction checks (The cell to the left (-x) should be pulled to the
        // right (+x))
        REQUIRE(state.pm_gravity_x(0, 1, 1) > 0.0);

        // Magnitude symmetry check
        REQUIRE(state.pm_gravity_x(0, 1, 1) ==
                Catch::Approx(-state.pm_gravity_x(2, 1, 1)));
    }
}