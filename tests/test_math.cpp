#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "types.h"
#include "math_utils.h"

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



