#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "config.h"
#include "gas.h"
#include "integrator.h"
#include "math_utils.h"
#include "particles.h"
#include "types.h"

// ---------------------------------------------------------------------
// Particle Spatial Displacement (Periodic Boundaries)
// ---------------------------------------------------------------------
TEST_CASE("Periodic displacement calculates the shortest path",
          "[boundaries][math]") {
    const double box_size = 64.0;

    SECTION("Particles far apart wrap around the boundary") {
        // Particle A is at x=1, Particle B is at x=63
        // If A looks for B, B is to its "left" across the boundary
        double dx = 63.0 - 1.0;

        double dist = periodic_displacement(dx, box_size);

        // The shortest path is -2.0 (wrapping backwards past 0)
        REQUIRE(dist == Catch::Approx(-2.0));

        // If B looks for A (1.0 - 63.0)
        double reverse_dist = periodic_displacement(-dx, box_size);
        REQUIRE(reverse_dist == Catch::Approx(2.0));
    }

    SECTION("Particles close together do not wrap") {
        // Particle A at 32, Particle B at 30
        double dx = 30.0 - 32.0;

        double dist = periodic_displacement(dx, box_size);

        // Normal distance applies
        REQUIRE(dist == Catch::Approx(-2.0));
    }

    SECTION("Exact half-box distance edge case") {
        // Distance is exactly half the box size
        double dx = 32.0;
        double dist = periodic_displacement(dx, box_size);

        // Should resolve to either -32 or +32 depending on fmod behavior
        // We use absolute value to check the magnitude
        REQUIRE(std::abs(dist) == Catch::Approx(32.0));
    }
}

/// ---------------------------------------------------------------------
// Gas Grid Rolling (Periodic Shift)
// ---------------------------------------------------------------------
TEST_CASE("Grid rolling correctly wraps memory boundaries",
          "[boundaries][hydro]") {
    // We use a tiny 4x4x4 grid for testing.
    // It's small enough to do the mental math, but large enough to have a
    // distinct middle.
    const int N = 4;
    Grid3D grid(N);

    // Manually fill the grid with zeros
    for (int i = 0; i < N * N * N; ++i) {
        grid.data[i] = 0.0;
    }

    // Set exactly ONE cell to 1.0.
    // We'll place it on the X-boundary: (x=3, y=0, z=0)
    int x = 3, y = 0, z = 0;
    size_t edge_idx = x * N * N + y * N + z;
    grid.data[edge_idx] = 1.0;

    SECTION("Rolling +1 along X axis wraps to the opposite face") {
        // Roll the grid +1 cell along axis 0 (X-axis)
        Grid3D shifted_grid = grid.roll(1, 0);

        // The value 1.0 should have wrapped around to (x=0, y=0, z=0)
        size_t wrapped_idx = 0 * N * N + 0 * N + 0;

        REQUIRE(shifted_grid.data[wrapped_idx] == Catch::Approx(1.0));

        // The old position at x=3 should now be empty
        REQUIRE(shifted_grid.data[edge_idx] == Catch::Approx(0.0));
    }

    SECTION("Rolling -1 along X axis shifts normally without wrapping") {
        // Roll the grid -1 cell along axis 0 (X-axis)
        Grid3D shifted_grid = grid.roll(-1, 0);

        // The value 1.0 should step backwards to (x=2, y=0, z=0)
        size_t normal_idx = 2 * N * N + 0 * N + 0;

        REQUIRE(shifted_grid.data[normal_idx] == Catch::Approx(1.0));
    }

    SECTION("Rolling along Y axis shifts entire memory strides correctly") {
        // Roll the grid +1 cell along axis 1 (Y-axis)
        Grid3D shifted_grid = grid.roll(1, 1);

        // The value 1.0 should step sideways to (x=3, y=1, z=0)
        size_t stride_idx = 3 * N * N + 1 * N + 0;

        REQUIRE(shifted_grid.data[stride_idx] == Catch::Approx(1.0));
    }
}

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
    config.TOTAL_MASS = 1.0;
    config.NUM_DM_PARTICLES = 1;

    // Create the system using the config
    ParticleSystem sys(config);

    // Helper function to sum total mass in the grid
    auto get_total_grid_mass = [&]() {
        double total_mass = 0.0;
        for (int i = 0;
             i < config.MESH_SIZE * config.MESH_SIZE * config.MESH_SIZE; ++i) {
            total_mass += sys.dm_rho.data[i] * config.CELL_VOLUME;
        }
        return total_mass;
    };

    SECTION("Particle exactly on a grid vertex (0 fractional part)") {
        Particle p;
        p.pos = {2.0, 4.0, 6.0};  // Exactly on vertex (ix=1, iy=2, iz=3)
        p.mass = 16.0;

        sys.particles.clear();
        sys.particles.push_back(p);
        sys.bin_and_assign_mass(config);

        // 1. Conservation Check
        REQUIRE(get_total_grid_mass() == Catch::Approx(16.0));

        // 2. Spatial Check: 100% of the mass should be in cell (1, 2, 3)
        double expected_density = 16.0 / config.CELL_VOLUME;
        REQUIRE(sys.dm_rho(1, 2, 3) == Catch::Approx(expected_density));

        // A neighboring cell should be strictly zero
        REQUIRE(sys.dm_rho(2, 2, 3) == Catch::Approx(0.0));
    }

    SECTION("Particle in dead center splits perfectly into 8 cells") {
        Particle p;
        p.pos = {3.0, 3.0, 3.0};  // Dead center of cell (ix=1, iy=1, iz=1)
        p.mass = 8.0;

        sys.particles.clear();
        sys.particles.push_back(p);
        sys.bin_and_assign_mass(config);

        // 1. Conservation Check
        REQUIRE(get_total_grid_mass() == Catch::Approx(8.0));

        // 2. Spatial Check: Mass should be split into exactly 8 chunks of 1.0
        double expected_density = 1.0 / config.CELL_VOLUME;

        // Check the 8 corners of the CIC cube
        REQUIRE(sys.dm_rho(1, 1, 1) == Catch::Approx(expected_density));
        REQUIRE(sys.dm_rho(2, 1, 1) == Catch::Approx(expected_density));
        REQUIRE(sys.dm_rho(1, 2, 1) == Catch::Approx(expected_density));
        REQUIRE(sys.dm_rho(2, 2, 1) == Catch::Approx(expected_density));
        REQUIRE(sys.dm_rho(1, 1, 2) == Catch::Approx(expected_density));
        REQUIRE(sys.dm_rho(2, 1, 2) == Catch::Approx(expected_density));
        REQUIRE(sys.dm_rho(1, 2, 2) == Catch::Approx(expected_density));
        REQUIRE(sys.dm_rho(2, 2, 2) == Catch::Approx(expected_density));
    }

    SECTION("Particle on the periodic boundary wraps to the opposite side") {
        Particle p;
        p.pos = {7.0, 0.0,
                 0.0};  // Straddling the X-axis boundary (x=7 out of 8)
        p.mass = 10.0;

        sys.particles.clear();
        sys.particles.push_back(p);
        sys.bin_and_assign_mass(config);

        // 1. Conservation Check
        REQUIRE(get_total_grid_mass() == Catch::Approx(10.0));

        // 2. Spatial Wrapping Check
        double expected_density = 5.0 / config.CELL_VOLUME;

        // Half goes to the rightmost cell (ix=3)
        REQUIRE(sys.dm_rho(3, 0, 0) == Catch::Approx(expected_density));

        // Half wraps around to the leftmost cell (ix=0)
        REQUIRE(sys.dm_rho(0, 0, 0) == Catch::Approx(expected_density));
    }
}

// ---------------------------------------------------------------------
// HLL Riemann Solver Core Logic
// ---------------------------------------------------------------------
struct RiemannSolverTestAccess {
    static Grid3D &density(RiemannSolver &s) { return s.density; }
    static Grid3D &S_L(RiemannSolver &s) { return s.S_L; }
    static Grid3D &S_R(RiemannSolver &s) { return s.S_R; }
    static Grid3D &S_R_minus_S_L(RiemannSolver &s) { return s.S_R_minus_S_L; }
};

TEST_CASE("HLL Solver calculates correct numerical fluxes",
          "[hydro][riemann]") {
    int N = 2;
    RiemannSolver solver(N);

    // Use Test Proxy to initialize the private variables
    RiemannSolverTestAccess::density(solver) = Grid3D(N);
    RiemannSolverTestAccess::S_L(solver) = Grid3D(N);
    RiemannSolverTestAccess::S_R(solver) = Grid3D(N);
    RiemannSolverTestAccess::S_R_minus_S_L(solver) = Grid3D(N);

    Grid3D FL(N), FR(N), UL(N), UR(N);

    SECTION("Regime 1: Supersonic flow to the Right (S_L >= 0)") {
        RiemannSolverTestAccess::S_L(solver).fill(1.5);
        RiemannSolverTestAccess::S_R(solver).fill(3.0);

        FL.fill(10.0);  // Expect this to be the answer
        FR.fill(20.0);
        UL.fill(1.0);
        UR.fill(2.0);

        Grid3D flux = solver.solve_hll(FL, FR, UL, UR);
        REQUIRE(flux(0, 0, 0) == Catch::Approx(10.0));
    }

    SECTION("Regime 2: Supersonic flow to the Left (S_R <= 0)") {
        RiemannSolverTestAccess::S_L(solver).fill(-3.0);
        RiemannSolverTestAccess::S_R(solver).fill(-1.5);

        FL.fill(10.0);
        FR.fill(20.0);  // Expect this to be the answer
        UL.fill(1.0);
        UR.fill(2.0);

        Grid3D flux = solver.solve_hll(FL, FR, UL, UR);
        REQUIRE(flux(0, 0, 0) == Catch::Approx(20.0));
    }

    SECTION("Regime 3: Subsonic / Mixed flow (S_L < 0 < S_R)") {
        RiemannSolverTestAccess::S_L(solver).fill(-1.0);
        RiemannSolverTestAccess::S_R(solver).fill(2.0);
        RiemannSolverTestAccess::S_R_minus_S_L(solver).fill(3.0);

        FL.fill(4.0);
        FR.fill(10.0);
        UL.fill(2.0);
        UR.fill(5.0);

        Grid3D flux = solver.solve_hll(FL, FR, UL, UR);
        REQUIRE(flux(0, 0, 0) == Catch::Approx(4.0));
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
    sys.particles = {p1, p2};

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
        state.dm.dm_rho.fill(5.0);

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
        state.dm.dm_rho.setZero();

        // Place a massive density spike at (ix=1, iy=1, iz=1)
        state.dm.dm_rho(1, 1, 1) = 100.0;

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
        REQUIRE(state.gravity_x(0, 1, 1) == Catch::Approx(-state.gravity_x(2, 1, 1)));
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

    sys.particles.push_back(p);

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