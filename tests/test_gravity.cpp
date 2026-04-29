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
    // 1. Setup a tiny mock universe
    Config config;
    config.MESH_SIZE = 4;
    config.CELL_SIZE = 2.0;
    config.CELL_VOLUME = 2.0 * 2.0 * 2.0;
    config.NUM_DM_PARTICLES = 1;

    // Helper function to sum total mass in the grid
    auto get_total_grid_mass = [&](const ParticleSystem& sys) {
        double total_mass = 0.0;
        for (int i = 0; i < config.MESH_SIZE * config.MESH_SIZE * config.MESH_SIZE; ++i) {
            total_mass += sys.get_rho().data[i] * config.CELL_VOLUME;
        }
        return total_mass;
    };

    SECTION("Particle exactly on a grid vertex (0 fractional part)") {
        ParticleSystem sys(config);
        
        // Add particle exactly on vertex (ix=1, iy=2, iz=3)
        // px, py, pz, vx, vy, vz, mass
        sys.add_particle(2.0, 4.0, 6.0, 0.0, 0.0, 0.0, 16.0); 
        sys.bin_and_assign_mass(config);

        // 1. Conservation Check
        REQUIRE(get_total_grid_mass(sys) == Catch::Approx(16.0));

        // 2. Spatial Check: 100% of the mass should be in cell (1, 2, 3)
        double expected_density = 16.0 / config.CELL_VOLUME;
        REQUIRE(sys.get_rho()(1, 2, 3) == Catch::Approx(expected_density));

        // A neighboring cell should be strictly zero
        REQUIRE(sys.get_rho()(2, 2, 3) == Catch::Approx(0.0));
    }

    SECTION("Particle in dead center splits perfectly into 8 cells") {
        ParticleSystem sys(config);
        
        // Dead center of cell (ix=1, iy=1, iz=1)
        sys.add_particle(3.0, 3.0, 3.0, 0.0, 0.0, 0.0, 8.0);
        sys.bin_and_assign_mass(config);

        // 1. Conservation Check
        REQUIRE(get_total_grid_mass(sys) == Catch::Approx(8.0));

        // 2. Spatial Check: Mass should be split into exactly 8 chunks of 1.0
        double expected_density = 1.0 / config.CELL_VOLUME;

        // Check the 8 corners of the CIC cube
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
        
        // Straddling the X-axis boundary (x=7 out of 8)
        sys.add_particle(7.0, 0.0, 0.0, 0.0, 0.0, 0.0, 10.0);
        sys.bin_and_assign_mass(config);

        // 1. Conservation Check
        REQUIRE(get_total_grid_mass(sys) == Catch::Approx(10.0));

        // 2. Spatial Wrapping Check
        double expected_density = 5.0 / config.CELL_VOLUME;

        // Half goes to the rightmost cell (ix=3)
        REQUIRE(sys.get_rho()(3, 0, 0) == Catch::Approx(expected_density));

        // Half wraps around to the leftmost cell (ix=0)
        REQUIRE(sys.get_rho()(0, 0, 0) == Catch::Approx(expected_density));
    }
}

// ---------------------------------------------------------------------
// Short-Range Gravity (PP / P3M)
// ---------------------------------------------------------------------
TEST_CASE("Short-range gravity calculates Newtonian and P3M forces",
          "[gravity][pp]") {
    Config config;
    config.MESH_SIZE = 4;
    config.CELL_SIZE = 2.0;
    config.DOMAIN_SIZE = 8.0;
    config.NUM_DM_PARTICLES = 2;

    // Set G=1.0 and softenings to 0.0 to make the mental math exact
    config.G = 1.0;
    config.SOFTENING_SQUARED = 0.0;

    ParticleSystem sys(config);
    // Two particles exactly 3.0 units apart on the X-axis
    sys.add_particle(2.0, 2.0, 2.0, 0.0, 0.0, 0.0, 1.0);
    sys.add_particle(5.0, 2.0, 2.0, 0.0, 0.0, 0.0, 1.0);

    SECTION("Pure Newtonian PP Force (USE_PM = false)") {
        config.USE_PM = false;
        config.CUTOFF_RADIUS_SQUARED = (config.DOMAIN_SIZE / 2.0) * (config.DOMAIN_SIZE / 2.0);

        sys.bin_and_assign_mass(config);
        sys.compute_pp_forces(config);

        // Exact Newtonian Gravity: F = G * m1 * m2 / r^2
        // F = 1.0 * 1.0 * 1.0 / (3.0 * 3.0) = 1.0 / 9.0 = 0.111111...
        double expected_f = 1.0 / 9.0;
        
        // Convert to acceleration (a = F/m). Since m=1, a=F.
        // p1 is pulled to the right (+x) by p2
        REQUIRE(sys.acc_x[0] == Catch::Approx(expected_f));
        REQUIRE(sys.acc_y[0] == Catch::Approx(0.0));
        REQUIRE(sys.acc_z[0] == Catch::Approx(0.0));

        // p2 is pulled to the left (-x) by p1
        REQUIRE(sys.acc_x[1] == Catch::Approx(-expected_f));
    }

    SECTION("P3M Cutoff safely ignores distant particles") {
        config.USE_PM = true;
        config.CUTOFF_RADIUS = 2.0;
        config.CUTOFF_RADIUS_SQUARED = 4.0;
        config.R_SWITCH_START = 1.0;
        config.R_SWITCH_START_SQ = 1.0;
        config.CUTOFF_TRANSITION_WIDTH = 1.0;

        sys.bin_and_assign_mass(config);

        // Reset accelerations to 0 before computing
        std::fill(sys.acc_x.begin(), sys.acc_x.end(), 0.0);
        
        sys.compute_pp_forces(config);

        // Because dist (3.0) > CUTOFF_RADIUS (2.0), the force should be exactly 0.0
        REQUIRE(sys.acc_x[0] == Catch::Approx(0.0));
        REQUIRE(sys.acc_x[1] == Catch::Approx(0.0));
    }
}

/// ---------------------------------------------------------------------
// Long-Range Gravity (PM Poisson Solver & CIC Interpolation)
// ---------------------------------------------------------------------
TEST_CASE("PM Poisson Solver calculates accurate grid accelerations",
          "[gravity][pm]") {
    Config config;
    config.MESH_SIZE = 4;
    config.CELL_SIZE = 2.0;
    config.DOMAIN_SIZE = 8.0;
    config.NUM_DM_PARTICLES = 1;
    config.G = 1.0;
    config.USE_HYDRO = false;

    SimState state(config);

    SECTION("A uniform density field produces zero net acceleration") {
        state.dm.dm_rho.fill(5.0); 
        
        // Compute gravity (it will copy dm_rho into total_rho internally)
        compute_gravitational_acceleration(state, config);

        REQUIRE(state.gravity_x.maxCoeff() == Catch::Approx(0.0).margin(1e-10));
        REQUIRE(state.gravity_x.minCoeff() == Catch::Approx(0.0).margin(1e-10));
        REQUIRE(state.gravity_y.maxCoeff() == Catch::Approx(0.0).margin(1e-10));
    }

    SECTION("A single point mass creates symmetric attractive forces") {
        // Place a massive particle exactly on the grid vertex so 
        // 100% of its mass falls into cell (1, 1, 1) during CIC
        state.dm.add_particle(2.0, 2.0, 2.0, 0.0, 0.0, 0.0, 100.0);
        
        // Run the binning algorithm to populate dm_rho
        state.dm.bin_and_assign_mass(config);

        // Now compute gravity. It will successfully copy the 100.0 spike
        compute_gravitational_acceleration(state, config);

        // Center of mass check (Should be exactly 0 due to central differences)
        REQUIRE(state.gravity_x(1, 1, 1) == Catch::Approx(0.0).margin(1e-10));

        // Attraction checks (The cell to the right (+x) should be pulled to the left (-x))
        REQUIRE(state.gravity_x(2, 1, 1) < 0.0);

        // Attraction checks (The cell to the left (-x) should be pulled to the right (+x))
        REQUIRE(state.gravity_x(0, 1, 1) > 0.0);

        // Magnitude symmetry check
        REQUIRE(state.gravity_x(0, 1, 1) == Catch::Approx(-state.gravity_x(2, 1, 1)));
    }
}

TEST_CASE("CIC Force Interpolation correctly computes F = m * a",
          "[gravity][cic]") {
    Config config;
    config.MESH_SIZE = 4;
    config.CELL_SIZE = 2.0;

    ParticleSystem sys(config);
    sys.add_particle(3.0, 3.0, 3.0, 0.0, 0.0, 0.0, 5.0); // Dead center of cell

    sys.bin_and_assign_mass(config);

    int N = config.MESH_SIZE;
    Grid3D acc_x(N), acc_y(N), acc_z(N);

    // Set a uniform acceleration field pointing diagonally
    acc_x.fill(2.0);
    acc_y.fill(-1.5);
    acc_z.fill(0.0);

    sys.interpolate_cic_forces(acc_x, acc_y, acc_z, config);

    // The interpolated acceleration directly modifies the SoA arrays
    REQUIRE(sys.acc_x[0] == Catch::Approx(2.0));  
    REQUIRE(sys.acc_y[0] == Catch::Approx(-1.5)); 
    REQUIRE(sys.acc_z[0] == Catch::Approx(0.0));   
}

TEST_CASE("Gravitational forces are physically scaled",
          "[gravity][integrator]") {
    Config config;
    config.MESH_SIZE = 4;
    config.N_PER_SIDE = 2;
    config.STANDING_PARTICLES = false;
    config.START_A = 0.5;
    config.compute_derived_data();

    SimState state(config);
    state.total_time = 0;
    state.scale_factor = 0.5;  // Mid-way through the simulation
    update_cosmology(state, config);

    // Two close particles
    state.dm.add_particle(0.5, 0.5, 0.4, 0.0, 0.0, 0.0, 0.5);
    state.dm.add_particle(0.5, 0.5, 0.6, 0.0, 0.0, 0.0, 0.5);

    // Run gravity solvers
    state.dm.bin_and_assign_mass(config);
    state.total_rho = state.dm.get_rho(); // Simulate the copy that integrator does
    
    compute_gravitational_acceleration(state, config);

    // Clear accelerations, then add PM and PP
    std::fill(state.dm.acc_x.begin(), state.dm.acc_x.end(), 0.0);
    std::fill(state.dm.acc_y.begin(), state.dm.acc_y.end(), 0.0);
    std::fill(state.dm.acc_z.begin(), state.dm.acc_z.end(), 0.0);
    
    state.dm.interpolate_cic_forces(state.gravity_x, state.gravity_y, state.gravity_z, config);
    state.dm.compute_pp_forces(config);

    // Check if the total acceleration is physically active
    REQUIRE(std::abs(state.dm.acc_z[0]) > 0.05);
}