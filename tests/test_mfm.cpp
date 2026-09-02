#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <vector>

#include "mfm.h"

TEST_CASE("MFM Gradient Estimator", "[hydro][mfm_gradients]") {
    // Grid Setup
    double spacing = 0.4;
    double cell_volume = spacing * spacing * spacing;
    double domain_size = 10.0;

    // Setup Central Particle
    ParticleState p_i;
    p_i.pos = Eigen::Vector3d(5.0, 5.0, 5.0);
    p_i.h = 0.7;

    // Evaluate fields at central position for baseline
    // Rho: 10.0 + 2x + 3y + 4z
    // P:   20.0 - 1x + 5y - 2z
    // Vx:  3.0x
    p_i.rho = 10.0 + 2.0 * 5.0 + 3.0 * 5.0 + 4.0 * 5.0;
    p_i.pressure = 20.0 - 1.0 * 5.0 + 5.0 * 5.0 - 2.0 * 5.0;
    p_i.vel = Eigen::Vector3d(3.0 * 5.0, 0.0, 0.0);

    p_i.mass = p_i.rho * cell_volume;

    SECTION("Exact linear field recovery in a uniform 3D distribution") {
        std::vector<ParticleState> neighbors;

        for (int i = -1; i <= 1; ++i) {
            for (int j = -1; j <= 1; ++j) {
                for (int k = -1; k <= 1; ++k) {
                    if (i == 0 && j == 0 && k == 0) continue;  // Skip self

                    ParticleState nj;
                    nj.pos = p_i.pos + Eigen::Vector3d(i * spacing, j * spacing,
                                                       k * spacing);
                    nj.h = 1.0;

                    // Assign properties based on the linear field
                    nj.rho = 10.0 + 2.0 * nj.pos.x() + 3.0 * nj.pos.y() +
                             4.0 * nj.pos.z();
                    nj.pressure = 20.0 - 1.0 * nj.pos.x() + 5.0 * nj.pos.y() -
                                  2.0 * nj.pos.z();
                    nj.vel = Eigen::Vector3d(3.0 * nj.pos.x(), 0.0, 0.0);

                    // Mass must be consistent with density and spacing
                    nj.mass = nj.rho * cell_volume;

                    neighbors.push_back(nj);
                }
            }
        }

        ParticleGradients grads =
            compute_single_particle_gradients(p_i, neighbors, domain_size);

        // E matrix should be successfully inverted
        REQUIRE(grads.B_matrix.norm() > 0.0);

        // Density Gradient should exactly equal (2, 3, 4)
        REQUIRE(grads.grad_rho.x() == Catch::Approx(2.0).margin(1e-12));
        REQUIRE(grads.grad_rho.y() == Catch::Approx(3.0).margin(1e-12));
        REQUIRE(grads.grad_rho.z() == Catch::Approx(4.0).margin(1e-12));

        // Pressure Gradient should exactly equal (-1, 5, -2)
        REQUIRE(grads.grad_p.x() == Catch::Approx(-1.0).margin(1e-12));
        REQUIRE(grads.grad_p.y() == Catch::Approx(5.0).margin(1e-12));
        REQUIRE(grads.grad_p.z() == Catch::Approx(-2.0).margin(1e-12));

        // Vx Gradient should exactly equal (3, 0, 0)
        REQUIRE(grads.grad_vx.x() == Catch::Approx(3.0).margin(1e-12));
        REQUIRE(grads.grad_vx.y() == Catch::Approx(0.0).margin(1e-12));
        REQUIRE(grads.grad_vx.z() == Catch::Approx(0.0).margin(1e-12));
    }

    SECTION("Safely falls back to SPH gradients for pathological geometries") {
        std::vector<ParticleState> neighbors;

        extern double g_pressure_floor;
        g_pressure_floor = 1e-16;

        // Only insert neighbors aligned on the X-axis (1D line)
        for (int i = -1; i <= 1; i += 2) {
            ParticleState nj = p_i;  // Copy base state
            nj.pos.x() += i * spacing;
            nj.rho = 10.0 + 2.0 * nj.pos.x();  // Still linear
            nj.mass = nj.rho * cell_volume;
            neighbors.push_back(nj);
        }

        ParticleGradients grads =
            compute_single_particle_gradients(p_i, neighbors, domain_size);

        // B matrix should evaluate to the Identity matrix during the fallback
        REQUIRE(grads.B_matrix(0, 0) == 1.0);
        REQUIRE(grads.B_matrix(1, 1) == 1.0);
        REQUIRE(grads.B_matrix(2, 2) == 1.0);

        // Off-diagonals should be zero
        REQUIRE(grads.B_matrix(0, 1) == 0.0);

        // The SPH fallback should successfully compute an X-gradient,
        // so we only strictly assert that the Y and Z gradients remain zero.
        REQUIRE(grads.grad_rho.y() == 0.0);
        REQUIRE(grads.grad_rho.z() == 0.0);

        REQUIRE(grads.grad_p.y() == 0.0);
        REQUIRE(grads.grad_p.z() == 0.0);

        REQUIRE(grads.grad_vx.y() == 0.0);
        REQUIRE(grads.grad_vx.z() == 0.0);
    }
}

TEST_CASE("MFM Face Reconstruction and Limiting",
          "[hydro][mfm_reconstruction]") {
    extern double g_pressure_floor;
    g_pressure_floor = 1e-16;

    double domain_size = 10.0;

    // Setup Particle i (Left)
    ParticleState p_i;
    p_i.pos = Eigen::Vector3d(1.0, 0.0, 0.0);
    p_i.rho = 1.0;
    p_i.pressure = 2.0;
    p_i.vel = Eigen::Vector3d::Zero();
    p_i.h = 0.5;

    ParticleGradients grad_i;
    grad_i.grad_rho = Eigen::Vector3d::Zero();
    grad_i.grad_p = Eigen::Vector3d::Zero();
    grad_i.grad_vx = Eigen::Vector3d::Zero();
    grad_i.grad_vy = Eigen::Vector3d::Zero();
    grad_i.grad_vz = Eigen::Vector3d::Zero();

    // Setup Particle j (Right)
    ParticleState p_j;
    p_j.pos = Eigen::Vector3d(3.0, 0.0, 0.0);  // Distance of 2.0
    p_j.rho = 2.0;
    p_j.pressure = 4.0;
    p_j.vel = Eigen::Vector3d::Zero();
    p_j.h = 0.5;

    ParticleGradients grad_j = grad_i;  // Copy zeros

    SECTION("First-order fallback (Zero Gradients)") {
        // With zero gradients, the extrapolated face states should simply equal
        // the base particle states
        ReconstructedFace face =
            compute_face_reconstruction(p_i, grad_i, p_j, grad_j, domain_size);

        REQUIRE(face.is_valid == true);
        REQUIRE(face.r == Catch::Approx(2.0));
        REQUIRE(face.n.x() == Catch::Approx(1.0));

        REQUIRE(face.rho_L == Catch::Approx(1.0));
        REQUIRE(face.rho_R == Catch::Approx(2.0));
    }

    SECTION("Smooth second-order extrapolation") {
        // Set gradients pointing along the correct slope
        // dx_half is 1.0. delta_rho is 1.0.
        grad_i.grad_rho.x() =
            0.5;  // Will extrapolate i's density by +0.5 to 1.5
        grad_j.grad_rho.x() =
            0.5;  // Will extrapolate j's density by -0.5 to 1.5

        ReconstructedFace face =
            compute_face_reconstruction(p_i, grad_i, p_j, grad_j, domain_size);

        // They should meet in the middle
        REQUIRE(face.rho_L == Catch::Approx(1.5));
        REQUIRE(face.rho_R == Catch::Approx(1.5));
    }

    SECTION("Vacuum flooring protects against negative states") {
        // Set negative pressure
        p_i.pressure = -5.0;

        ReconstructedFace face =
            compute_face_reconstruction(p_i, grad_i, p_j, grad_j, domain_size);

        // Should catch the negative pressure and clamp it to 1e-16
        REQUIRE(face.p_L == Catch::Approx(1e-16));
    }
}

TEST_CASE("MFM Riemann Solver", "[hydro][mfm_riemann]") {
    double gamma = 5.0 / 3.0;

    ReconstructedFace face;
    face.n = Eigen::Vector3d(1.0, 0.0, 0.0);
    face.is_valid = true;

    // Default Face Frame Velocity
    Eigen::Vector3d v_frame = Eigen::Vector3d::Zero();

    SECTION("Sod Shock Tube Contact Wave (Subsonic)") {
        // High pressure left, low pressure right, initially at rest
        face.rho_L = 1.0;
        face.p_L = 1.0;
        face.v_L = Eigen::Vector3d::Zero();
        face.rho_R = 0.125;
        face.p_R = 0.1;
        face.v_R = Eigen::Vector3d::Zero();

        MFMFaceFlux flux = solve_mfm_riemann(face, v_frame, gamma);

        // P_star must resolve strictly between P_R (0.1) and P_L (1.0).
        REQUIRE(flux.P_star > 0.1);
        REQUIRE(flux.P_star < 1.0);

        // Momentum flux is P_star * n.
        REQUIRE(flux.flux_mom.x() == Catch::Approx(flux.P_star));
        REQUIRE(flux.flux_mom.y() == 0.0);
        REQUIRE(flux.flux_mom.z() == 0.0);

        // Since the high pressure pushes right, the contact wave (S_star) is
        // positive.
        REQUIRE(flux.S_star > 0.0);
    }

    SECTION("Supersonic flow exclusively trusts Upwind State") {
        // Both states moving to the right (Mach 10+)
        face.rho_L = 1.0;
        face.p_L = 1.0;
        face.v_L = Eigen::Vector3d(100.0, 0.0, 0.0);
        face.rho_R = 1.0;
        face.p_R = 2.0;
        face.v_R = Eigen::Vector3d(100.0, 0.0, 0.0);

        MFMFaceFlux flux = solve_mfm_riemann(face, v_frame, gamma);

        // Because flow is supersonic to the right, S_L > 0.
        // Solver MUST pick exactly the Left state (Upwinding).
        REQUIRE(flux.P_star == Catch::Approx(1.0));  // P_star = P_L = 1.0
        REQUIRE(flux.flux_mom.x() == Catch::Approx(1.0));

        // S_star defaults to the upwind velocity in highly supersonic flow
        REQUIRE(flux.S_star == Catch::Approx(100.0));
    }
}