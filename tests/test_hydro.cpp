#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "gas.h"

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

    Grid3D FL(N), FR(N), UL(N), UR(N), flux(N);

    SECTION("Regime 1: Supersonic flow to the Right (S_L >= 0)") {
        RiemannSolverTestAccess::S_L(solver).fill(1.5);
        RiemannSolverTestAccess::S_R(solver).fill(3.0);

        FL.fill(10.0);  // Expect this to be the answer
        FR.fill(20.0);
        UL.fill(1.0);
        UR.fill(2.0);

        solver.solve_hll(FL, FR, UL, UR, flux);
        REQUIRE(flux(0, 0, 0) == Catch::Approx(10.0));
    }

    SECTION("Regime 2: Supersonic flow to the Left (S_R <= 0)") {
        RiemannSolverTestAccess::S_L(solver).fill(-3.0);
        RiemannSolverTestAccess::S_R(solver).fill(-1.5);

        FL.fill(10.0);
        FR.fill(20.0);  // Expect this to be the answer
        UL.fill(1.0);
        UR.fill(2.0);

        solver.solve_hll(FL, FR, UL, UR, flux);
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

        solver.solve_hll(FL, FR, UL, UR, flux);
        REQUIRE(flux(0, 0, 0) == Catch::Approx(4.0));
    }
}
