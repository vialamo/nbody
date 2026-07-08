#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "config.h"
#include "gas.h"

// ---------------------------------------------------------------------
// Test Proxies to access private grid data
// ---------------------------------------------------------------------
struct GasGridTestAccess {
    static Grid3D& density(GasGrid& g) { return g.density; }
    static Grid3D& momentum_x(GasGrid& g) { return g.momentum_x; }
    static Grid3D& momentum_y(GasGrid& g) { return g.momentum_y; }
    static Grid3D& momentum_z(GasGrid& g) { return g.momentum_z; }
    static Grid3D& energy(GasGrid& g) { return g.energy; }
    static Grid3D& internal_energy(GasGrid& g) { return g.internal_energy; }
    static Grid3D& pressure(GasGrid& g) { return g.pressure; }
    static void update_primitive_variables(GasGrid& g) {
        g.update_primitive_variables();
    }
};

struct RiemannSolverTestAccess {
    static Grid3D& rho_L(RiemannSolver& s) { return s.rho_L; }
    static Grid3D& rho_R(RiemannSolver& s) { return s.rho_R; }
};

// ---------------------------------------------------------------------
// Dual Energy Formalism Tests
// ---------------------------------------------------------------------
TEST_CASE("Dual Energy Switch protects pressure in hypersonic flows",
          "[hydro][dual_energy]") {
    Config conf;
    conf.mesh_size = 2;
    conf.gamma = 5.0 / 3.0;
    GasGrid grid(conf);

    auto& rho = GasGridTestAccess::density(grid);
    auto& mx = GasGridTestAccess::momentum_x(grid);
    auto& en = GasGridTestAccess::energy(grid);
    auto& ie = GasGridTestAccess::internal_energy(grid);
    auto& pr = GasGridTestAccess::pressure(grid);

    // Setup uniform density and safe internal energy
    rho.fill(1.0);
    ie.fill(1.0);  // P should be (5/3 - 1) * 1.0 = 0.66666...

    SECTION("Thermal Dominated (Subsonic) - Trusts Total Energy") {
        // Kinetic energy is small compared to total energy
        mx.fill(1.0);  // v = 1.0, KE = 0.5
        en.fill(1.5);  // Total energy = IE (1.0) + KE (0.5)

        // Corrupt the tracked internal energy to prove it DOES NOT use it
        ie.fill(999.0);

        GasGridTestAccess::update_primitive_variables(grid);

        // It should have computed ie = total - kinetic
        // P = (gamma-1) * (E_tot - KE) = (2/3) * (1.5 - 0.5) = 0.6666...
        REQUIRE(pr(0, 0, 0) == Catch::Approx(2.0 / 3.0));
        // Tracked IE should be synced BACK to 1.0
        REQUIRE(ie(0, 0, 0) == Catch::Approx(1.0));
    }

    SECTION("Kinetic Dominated (Hypersonic) - Trusts Tracked Internal Energy") {
        // Kinetic energy is > 99.9% of total energy
        mx.fill(1000.0);    // v = 1000, KE = 500000
        ie.fill(1.0);       // Safe internal energy
        en.fill(500001.0);  // True total energy

        // Simulate catastrophic floating point cancellation by slightly
        // corrupting Total Energy
        en.fill(499999.0);  // If E_tot - KE is used, IE is negative!

        GasGridTestAccess::update_primitive_variables(grid);

        // It should have computed totak = ie + kinetic
        // Because KE is dominant, the switch should trigger and use tracked IE
        // instead. P = (gamma-1) * tracked_IE = (2/3) * 1.0 = 0.6666...
        REQUIRE(pr(0, 0, 0) == Catch::Approx(2.0 / 3.0));
        // Total energy should be safely synced to IE + KE
        REQUIRE(en(0, 0, 0) == Catch::Approx(500001.0));
    }
}

// ---------------------------------------------------------------------
// MUSCL Piecewise Linear Method Reconstruction & Minmod Limiter Tests
// ---------------------------------------------------------------------
TEST_CASE("MUSCL PLM applies Minmod limiter to extrema", "[hydro][muscl]") {
    Config conf;
    conf.mesh_size = 4;
    conf.gamma = 5.0 / 3.0;
    GasGrid grid(conf);

    auto& rho = GasGridTestAccess::density(grid);
    GasGridTestAccess::momentum_x(grid).fill(0.0);
    GasGridTestAccess::momentum_y(grid).fill(0.0);
    GasGridTestAccess::momentum_z(grid).fill(0.0);
    GasGridTestAccess::energy(grid).fill(1.0);
    GasGridTestAccess::internal_energy(grid).fill(1.0);

    // Create a density spike (Extremum)
    rho(0, 0, 0) = 1.0;
    rho(1, 0, 0) = 10.0;  // The spike!
    rho(2, 0, 0) = 1.0;
    rho(3, 0, 0) = 1.0;

    GasGridTestAccess::update_primitive_variables(grid);

    RiemannSolver solver(conf.mesh_size);
    solver.compute_fluxes(grid, 0, conf.gamma);  // X-Axis sweep

    auto& rho_L = RiemannSolverTestAccess::rho_L(solver);

    // At i=1, the slope from the left is +9, but from the right is -9.
    // The minmod limiter MUST return 0 for opposite signs.
    // Therefore, the reconstructed left state at the interface should exactly
    // equal the cell center.
    REQUIRE(rho_L(1, 0, 0) == Catch::Approx(10.0));
}

// ---------------------------------------------------------------------
// HLLC Riemann Solver Flux Tests
// ---------------------------------------------------------------------
TEST_CASE("HLLC correctly resolves upwind supersonic fluxes", "[hydro][hllc]") {
    // Upwinding in a Supersonic Flow. Prove that the HLLC solver
    // correctly realizes when a fluid is moving so fast that
    // it outruns its own sound waves,
    // preventing information from physically traveling backward.

    Config conf;
    conf.mesh_size = 4;
    conf.gamma = 5.0 / 3.0;
    GasGrid grid(conf);

    // Setup uniform supersonic flow to the right
    GasGridTestAccess::density(grid).fill(1.0);
    GasGridTestAccess::momentum_x(grid).fill(10.0);  // v = 10.0
    GasGridTestAccess::momentum_y(grid).fill(0.0);
    GasGridTestAccess::momentum_z(grid).fill(0.0);
    GasGridTestAccess::internal_energy(grid).fill(1.0);
    GasGridTestAccess::energy(grid).fill(51.0);  // IE(1.0) + KE(50.0)

    GasGridTestAccess::update_primitive_variables(grid);

    RiemannSolver solver(conf.mesh_size);
    solver.compute_fluxes(grid, 0, conf.gamma);

    // P = (gamma - 1) x ei = (conf.gamma(5/3) - 1) x 1 = 2/3
    // The speed of sound in this gas is:
    // cs = sqrt(gamma * P / rho) = sqrt((5/3)*(2/3)/1) ~ 1.05
    // v = 10.0. Therefore, the leftmost (slowest) wave speed in the
    // Riemann fan is S_L (approx) = v - c_s = 10.0 - 1.05 = 8.95.
    // Since S_L > 0, the entire wave structure moves right.
    // In this regime, HLLC flux must exactly equal F_L.

    // F_dens = rho * v = 1.0 * 10.0 = 10.0
    REQUIRE(solver.get_flux_density()(0, 0, 0) == Catch::Approx(10.0));

    // F_mom_n = rho * v^2 + P = 1.0 * 100.0 + 2/3 = 100.666...
    REQUIRE(solver.get_flux_mom_n()(0, 0, 0) ==
            Catch::Approx(100.0 + 2.0 / 3.0));

    // F_en = (E_tot + P) * v = (51.0 + 2/3) * 10.0 = 516.666...
    REQUIRE(solver.get_flux_energy()(0, 0, 0) ==
            Catch::Approx(516.66666666666));
}