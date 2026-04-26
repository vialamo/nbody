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
struct ParticleTestAccess {
    static Grid3D &dm_rho(ParticleSystem &s) { return s.dm_rho; }
    static std::vector<Particle> &particles(ParticleSystem &s) {
        return s.particles;
    }
};
TEST_CASE("CIC Density Assignment maps particles to grid correctly",
          "[cic][particles]") {
    // 1. Setup a tiny mock universe
    Config config;
    config.MESH_SIZE = 4;
    config.CELL_SIZE = 2.0;
    config.CELL_VOLUME = 2.0 * 2.0 * 2.0;
    config.NUM_DM_PARTICLES = 1;

    // Create the system using the config
    ParticleSystem sys(config);

    // Helper function to sum total mass in the grid
    auto get_total_grid_mass = [&]() {
        double total_mass = 0.0;
        for (int i = 0;
             i < config.MESH_SIZE * config.MESH_SIZE * config.MESH_SIZE; ++i) {
            total_mass += sys.get_rho().data[i] * config.CELL_VOLUME;
        }
        return total_mass;
    };

    SECTION("Particle exactly on a grid vertex (0 fractional part)") {
        Particle p;
        p.pos = {2.0, 4.0, 6.0};  // Exactly on vertex (ix=1, iy=2, iz=3)
        p.mass = 16.0;

        ParticleTestAccess::particles(sys).clear();
        ParticleTestAccess::particles(sys).push_back(p);
        sys.bin_and_assign_mass(config);

        // 1. Conservation Check
        REQUIRE(get_total_grid_mass() == Catch::Approx(16.0));

        // 2. Spatial Check: 100% of the mass should be in cell (1, 2, 3)
        double expected_density = 16.0 / config.CELL_VOLUME;
        REQUIRE(sys.get_rho()(1, 2, 3) == Catch::Approx(expected_density));

        // A neighboring cell should be strictly zero
        REQUIRE(sys.get_rho()(2, 2, 3) == Catch::Approx(0.0));
    }

    SECTION("Particle in dead center splits perfectly into 8 cells") {
        Particle p;
        p.pos = {3.0, 3.0, 3.0};  // Dead center of cell (ix=1, iy=1, iz=1)
        p.mass = 8.0;

        ParticleTestAccess::particles(sys).clear();
        ParticleTestAccess::particles(sys).push_back(p);
        sys.bin_and_assign_mass(config);

        // 1. Conservation Check
        REQUIRE(get_total_grid_mass() == Catch::Approx(8.0));

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
        Particle p;
        p.pos = {7.0, 0.0,
                 0.0};  // Straddling the X-axis boundary (x=7 out of 8)
        p.mass = 10.0;

        ParticleTestAccess::particles(sys).clear();
        ParticleTestAccess::particles(sys).push_back(p);
        sys.bin_and_assign_mass(config);

        // 1. Conservation Check
        REQUIRE(get_total_grid_mass() == Catch::Approx(10.0));

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

    // We will place two particles exactly 3.0 units apart on the X-axis
    Particle p1, p2;
    p1.pos = {2.0, 2.0, 2.0};
    p1.mass = 1.0;

    p2.pos = {5.0, 2.0, 2.0};
    p2.mass = 1.0;

    ParticleSystem sys(config);
    ParticleTestAccess::particles(sys) = {p1, p2};

    SECTION("Pure Newtonian PP Force (USE_PM = false)") {
        config.USE_PM = false;

        // Since we aren't using PM, the cutoff should theoretically be half the
        // box size squared
        config.CUTOFF_RADIUS_SQUARED =
            (config.DOMAIN_SIZE / 2.0) * (config.DOMAIN_SIZE / 2.0);

        // Populate the cell_grid so the neighbor search finds the particles
        sys.bin_and_assign_mass(config);

        std::vector<Vec3> forces;
        sys.compute_pp_forces(forces, config);

        // Exact Newtonian Gravity: F = G * m1 * m2 / r^2
        // F = 1.0 * 1.0 * 1.0 / (3.0 * 3.0) = 1.0 / 9.0 = 0.111111...
        double expected_f = 1.0 / 9.0;

        // p1 is pulled to the right (+x) by p2
        REQUIRE(forces[0].x == Catch::Approx(expected_f));
        REQUIRE(forces[0].y == Catch::Approx(0.0));
        REQUIRE(forces[0].z == Catch::Approx(0.0));

        // p2 is pulled to the left (-x) by p1 (Newton's 3rd Law)
        REQUIRE(forces[1].x == Catch::Approx(-expected_f));
    }

    SECTION("P3M Cutoff safely ignores distant particles") {
        config.USE_PM = true;
        // The distance is 3.0, so if we set cutoff to 2.0, the force MUST be
        // zero
        config.CUTOFF_RADIUS = 2.0;
        config.CUTOFF_RADIUS_SQUARED = 4.0;

        // Setup the switch transition variables just in case
        config.R_SWITCH_START = 1.0;
        config.R_SWITCH_START_SQ = 1.0;
        config.CUTOFF_TRANSITION_WIDTH = 1.0;

        sys.bin_and_assign_mass(config);

        std::vector<Vec3> forces;
        sys.compute_pp_forces(forces, config);

        // Because dist (3.0) > CUTOFF_RADIUS (2.0), the force should be exactly
        // 0.0
        REQUIRE(forces[0].x == Catch::Approx(0.0));
        REQUIRE(forces[1].x == Catch::Approx(0.0));
    }
}

// ---------------------------------------------------------------------
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
    config.USE_HYDRO = false;  // We just test the dark matter PM solver

    int N = config.MESH_SIZE;
    SimState state(config);

    SECTION("A uniform density field produces zero net acceleration") {
        // Fill the universe with a uniform density of 5.0
        ParticleTestAccess::dm_rho(state.dm).fill(5.0);

        // Compute the FFT Poisson solve
        compute_gravitational_acceleration(state, config);

        // Because k=0 is set to 0.0, and there are no other waves, acceleration
        // MUST be 0 everywhere
        REQUIRE(state.gravity_x.maxCoeff() == Catch::Approx(0.0).margin(1e-10));
        REQUIRE(state.gravity_x.minCoeff() == Catch::Approx(0.0).margin(1e-10));
        REQUIRE(state.gravity_y.maxCoeff() == Catch::Approx(0.0).margin(1e-10));
        REQUIRE(state.gravity_y.maxCoeff() == Catch::Approx(0.0).margin(1e-10));
    }

    SECTION("A single point mass creates symmetric attractive forces") {
        ParticleTestAccess::dm_rho(state.dm).setZero();

        // Place a massive density spike at (ix=1, iy=1, iz=1)
        ParticleTestAccess::dm_rho(state.dm)(1, 1, 1) = 100.0;

        compute_gravitational_acceleration(state, config);

        // 1. Center of mass check (Should be exactly 0 due to central
        // differences)
        REQUIRE(state.gravity_x(1, 1, 1) == Catch::Approx(0.0).margin(1e-10));

        // 2. Attraction check (The cell to the right (+x) should be pulled to
        // the left (-x))
        REQUIRE(state.gravity_x(2, 1, 1) < 0.0);

        // 3. Attraction check (The cell to the left (-x) should be pulled to
        // the right (+x))
        REQUIRE(state.gravity_x(0, 1, 1) > 0.0);

        // 4. Magnitude symmetry check
        REQUIRE(state.gravity_x(0, 1, 1) ==
                Catch::Approx(-state.gravity_x(2, 1, 1)));
    }
}

TEST_CASE("CIC Force Interpolation correctly computes F = m * a",
          "[gravity][cic]") {
    Config config;
    config.MESH_SIZE = 4;
    config.CELL_SIZE = 2.0;

    ParticleSystem sys(config);

    Particle p;
    p.pos = {3.0, 3.0, 3.0};  // Dead center of cell
    p.mass = 5.0;

    ParticleTestAccess::particles(sys).push_back(p);

    // MUST call this first to populate cic_data inside the ParticleSystem!
    sys.bin_and_assign_mass(config);

    int N = config.MESH_SIZE;
    Grid3D acc_x(N), acc_y(N), acc_z(N);

    // Set a uniform acceleration field pointing diagonally
    acc_x.fill(2.0);
    acc_y.fill(-1.5);
    acc_z.fill(0.0);

    std::vector<Vec3> forces;
    sys.interpolate_cic_forces(acc_x, acc_y, acc_z, forces, config);

    // Force = mass * acceleration
    REQUIRE(forces[0].x == Catch::Approx(10.0));  // 5.0 * 2.0
    REQUIRE(forces[0].y == Catch::Approx(-7.5));  // 5.0 * -1.5
    REQUIRE(forces[0].z == Catch::Approx(0.0));   // 5.0 * 0.0
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

    // Place two particles relatively close to each other
    Particle p1, p2;
    p1.pos = {0.5, 0.5, 0.4};
    p1.mass = 0.5;  // Half the universe mass
    p1.vel = {0, 0, 0};

    p2.pos = {0.5, 0.5, 0.6};
    p2.mass = 0.5;
    p2.vel = {0, 0, 0};

    state.dm.particles.push_back(p1);
    state.dm.particles.push_back(p2);

    // Run your gravity solvers!
    state.dm.bin_and_assign_mass(config);
    compute_gravitational_acceleration(state, config);

    std::vector<Vec3> pm_forces, pp_forces;
    state.dm.interpolate_cic_forces(state.gravity_x, state.gravity_y,
                                    state.gravity_z, pm_forces, config);
    state.dm.compute_pp_forces(pp_forces, config);

    double total_accel_z_p1 = (pm_forces[0].z + pp_forces[0].z) / p1.mass;

    // In Newtonian physics: a = G * m / r^2
    // r = 0.2. r^2 = 0.04.
    // a_newton = 0.0159 * 0.5 / 0.04 ≈ 0.198

    // Check if the force is completely dead (e.g., < 0.001)
    REQUIRE(std::abs(total_accel_z_p1) > 0.05);
}
