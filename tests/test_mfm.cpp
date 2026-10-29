#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <vector>

#include "config.h"
#include "density.h"
#include "hydro.h"
#include "kernels.h"
#include "mfm.h"
#include "particles.h"
#include "reconstruction.h"

// ---------------------------------------------------------------------
// MFM Face Reconstruction and Limiting
// ---------------------------------------------------------------------
TEST_CASE("MFM Face Reconstruction and Limiting", "[hydro][mfm_reconstruction]") {
    Config config;
    config.domain_size = Vec3{10.0, 10.0, 10.0};
    
    double density_floor = 1e-12;
    double pressure_floor = 1e-16;

    // Setup Particle i (Left)
    Reconstruction::ParticleState p_i;
    p_i.pos = Eigen::Vector3d(1.0, 0.0, 0.0);
    p_i.rho = 1.0;
    p_i.pressure = 2.0;
    p_i.vel = Eigen::Vector3d::Zero();
    p_i.h = 2.5;
    p_i.mass = 1.0;

    Reconstruction::ParticleGradients grad_i;
    grad_i.grad_rho = Eigen::Vector3d::Zero();
    grad_i.grad_p = Eigen::Vector3d::Zero();
    grad_i.grad_vx = Eigen::Vector3d::Zero();
    grad_i.grad_vy = Eigen::Vector3d::Zero();
    grad_i.grad_vz = Eigen::Vector3d::Zero();
    grad_i.ill_conditioned = false;

    // Setup Particle j (Right)
    Reconstruction::ParticleState p_j;
    p_j.pos = Eigen::Vector3d(3.0, 0.0, 0.0);  // Distance of 2.0
    p_j.rho = 2.0;
    p_j.pressure = 4.0;
    p_j.vel = Eigen::Vector3d::Zero();
    p_j.h = 2.5;
    p_j.mass = 2.0;

    Reconstruction::ParticleGradients grad_j = grad_i;  // Copy zeros

    Eigen::Matrix3d B_identity = Eigen::Matrix3d::Identity();

    SECTION("First-order fallback (Zero Gradients)") {
        // With zero gradients, the extrapolated face states should simply equal
        // the base particle states
        Reconstruction::ReconstructedFace face = Reconstruction::compute_face_reconstruction(
            p_i, grad_i, p_j, grad_j, B_identity, B_identity, config,
            density_floor, pressure_floor);

        REQUIRE(face.is_valid == true);
        REQUIRE(face.r == Catch::Approx(2.0));
        REQUIRE(face.n.x() == Catch::Approx(1.0));

        REQUIRE(face.rho_L == Catch::Approx(1.0));
        REQUIRE(face.rho_R == Catch::Approx(2.0));
    }

    SECTION("Smooth second-order extrapolation") {
        // Set gradients pointing along the correct slope
        // distance dx = 2.0, midpoint is at dx=1.0 from each particle
        grad_i.grad_rho.x() = 0.5;  // Extrapolate i's density by +0.5 to 1.5
        grad_j.grad_rho.x() = 0.5;  // Extrapolate j's density by -0.5 to 1.5

        Reconstruction::ReconstructedFace face = Reconstruction::compute_face_reconstruction(
            p_i, grad_i, p_j, grad_j, B_identity, B_identity, config,
            density_floor, pressure_floor);

        // They should meet perfectly in the middle
        REQUIRE(face.rho_L == Catch::Approx(1.5));
        REQUIRE(face.rho_R == Catch::Approx(1.5));
    }

    SECTION("Vacuum flooring protects against negative states") {
        p_i.pressure = -5.0; // Unphysical negative pressure

        Reconstruction::ReconstructedFace face = Reconstruction::compute_face_reconstruction(
            p_i, grad_i, p_j, grad_j, B_identity, B_identity, config,
            density_floor, pressure_floor);

        // Should catch the negative pressure and clamp it to the floor
        REQUIRE(face.p_L == Catch::Approx(pressure_floor));
    }
}

// ---------------------------------------------------------------------
// MFM Gradient Estimator
// ---------------------------------------------------------------------
TEST_CASE("MFM Gradient Estimator recovers linear fields", "[hydro][mfm_gradients]") {
    Config config;
    config.domain_size = {10.0, 10.0, 10.0};
    config.num_gas_particles = 27;
    config.compute_derived_data();
    
    GasParticleSystem gas(config);
    double spacing = 0.4;
    double cell_volume = spacing * spacing * spacing;
    
    // Create a 3x3x3 grid of particles
    for (int i = -1; i <= 1; ++i) {
        for (int j = -1; j <= 1; ++j) {
            for (int k = -1; k <= 1; ++k) {
                double px = 5.0 + i * spacing;
                double py = 5.0 + j * spacing;
                double pz = 5.0 + k * spacing;
                // Add dummy particle, properties overridden below
                gas.add_particle(px, py, pz, 0, 0, 0, 1.0, 1.0, 0.7, 0.0);
            }
        }
    }
    
    // Build the BVH to establish neighbors
    gas.build_lbvh(config);
    
    // Manually impose a linear field
    int central_idx = -1;
    for (size_t i = 0; i < gas.num_particles; ++i) {
        double px = gas.pos_x[i], py = gas.pos_y[i], pz = gas.pos_z[i];
        
        // Rho: 10.0 + 2x + 3y + 4z
        // P:   20.0 - 1x + 5y - 2z
        // Vx:  3.0x
        gas.rho[i] = 10.0 + 2.0 * px + 3.0 * py + 4.0 * pz;
        gas.pressure[i] = 20.0 - 1.0 * px + 5.0 * py - 2.0 * pz;
        gas.vel_x[i] = 3.0 * px;
        
        gas.pressure_pred[i] = gas.pressure[i];
        gas.vel_x_pred[i] = gas.vel_x[i];
        gas.vel_y_pred[i] = 0.0;
        gas.vel_z_pred[i] = 0.0;
        
        gas.mass[i] = gas.rho[i] * cell_volume;
        gas.cond_num[i] = 1.0;
        
        if (std::abs(px - 5.0) < 1e-5 && std::abs(py - 5.0) < 1e-5 && std::abs(pz - 5.0) < 1e-5) {
            central_idx = i;
        }
    }
    
    REQUIRE(central_idx != -1);

    // Compute the B_matrix for each particle to mock the DensitySolver
    for (size_t i = 0; i < gas.num_particles; ++i) {
        Eigen::Matrix3d E = Eigen::Matrix3d::Zero();
        for (size_t j = 0; j < gas.num_particles; ++j) {
            if (i == j) continue;
            
            double dx = gas.pos_x[j] - gas.pos_x[i];
            double dy = gas.pos_y[j] - gas.pos_y[i];
            double dz = gas.pos_z[j] - gas.pos_z[i];
            
            // Periodic wrap
            if (dx > 5.0) dx -= 10.0; else if (dx < -5.0) dx += 10.0;
            if (dy > 5.0) dy -= 10.0; else if (dy < -5.0) dy += 10.0;
            if (dz > 5.0) dz -= 10.0; else if (dz < -5.0) dz += 10.0;
            
            double r = std::sqrt(dx*dx + dy*dy + dz*dz);
            
            if (r < gas.h[i] && r > 1e-12) {
                double W;
                Kernels::cubic_spline_value(r, gas.h[i], W);
                double V_j = gas.mass[j] / gas.rho[j];
                double weight = W * V_j;
                
                E(0,0) += dx * dx * weight;
                E(0,1) += dx * dy * weight;
                E(0,2) += dx * dz * weight;
                E(1,1) += dy * dy * weight;
                E(1,2) += dy * dz * weight;
                E(2,2) += dz * dz * weight;
            }
        }
        E(1,0) = E(0,1);
        E(2,0) = E(0,2);
        E(2,1) = E(1,2);
        
        gas.B_matrix[i] = E.inverse();
    }
    
    // Compute gradients globally
    HydroSolver::compute_gradients(gas, config);
    
    // Density Gradient should equal (2, 3, 4)
    REQUIRE(gas.grad_rho[central_idx].x() == Catch::Approx(2.0).margin(1e-12));
    REQUIRE(gas.grad_rho[central_idx].y() == Catch::Approx(3.0).margin(1e-12));
    REQUIRE(gas.grad_rho[central_idx].z() == Catch::Approx(4.0).margin(1e-12));

    // Pressure Gradient should equal (-1, 5, -2)
    REQUIRE(gas.grad_p[central_idx].x() == Catch::Approx(-1.0).margin(1e-12));
    REQUIRE(gas.grad_p[central_idx].y() == Catch::Approx(5.0).margin(1e-12));
    REQUIRE(gas.grad_p[central_idx].z() == Catch::Approx(-2.0).margin(1e-12));

    // Vx Gradient should equal (3, 0, 0)
    REQUIRE(gas.grad_vx[central_idx].x() == Catch::Approx(3.0).margin(1e-12));
    REQUIRE(gas.grad_vx[central_idx].y() == Catch::Approx(0.0).margin(1e-12));
    REQUIRE(gas.grad_vx[central_idx].z() == Catch::Approx(0.0).margin(1e-12));
}

// ---------------------------------------------------------------------
// MFM Hydrodynamic Forces and Conservation
// ---------------------------------------------------------------------
TEST_CASE("MFM Hydro Solver conserves momentum and energy", "[hydro][mfm_forces]") {
    Config config;
    config.domain_size = {10.0, 10.0, 10.0};
    config.num_gas_particles = 2;
    config.compute_derived_data();
    
    GasParticleSystem gas(config);
    
    // Particle 0 (Left)
    gas.add_particle(4.8, 5.0, 5.0, 10.0, 0.0, 0.0, 1.0, 5.0, 0.5, 0.0);
    // Particle 1 (Right)
    gas.add_particle(5.2, 5.0, 5.0, -10.0, 0.0, 0.0, 1.0, 5.0, 0.5, 0.0);
    
    gas.build_lbvh(config);
    
    // Impose an arbitrary pressure/density differential
    gas.rho[0] = 1.0; gas.pressure[0] = 2.0;
    gas.rho[1] = 0.5; gas.pressure[1] = 1.0;
    
    // Mock the density solver dependencies
    for (size_t i = 0; i < 2; ++i) {
        gas.B_matrix[i] = Eigen::Matrix3d::Identity();
        gas.cond_num[i] = 1.0;
        gas.vel_x_pred[i] = gas.vel_x[i];
        gas.pressure_pred[i] = gas.pressure[i];
    }
    
    // We must run gradients first to prepare the arrays
    HydroSolver::compute_gradients(gas, config);
    
    // Run the Riemann solver and flux exchange
    double a = 1.0;
    HydroSolver::compute_forces(gas, config, a);
    
    // Momentum Conservation: m0 * a0 = -m1 * a1
    double force_x_0 = gas.mass[0] * gas.hydro_acc_x[0];
    double force_x_1 = gas.mass[1] * gas.hydro_acc_x[1];
    REQUIRE(force_x_0 == Catch::Approx(-force_x_1).margin(1e-14));
    
    // Ensure forces are actually being applied (particles are crashing into each other)
    REQUIRE(std::abs(force_x_0) > 0.0);
    
    // Energy Conservation: Rate of change of total energy must balance
    double power_0 = gas.mass[0] * gas.de_dt[0];
    double power_1 = gas.mass[1] * gas.de_dt[1];
    REQUIRE(power_0 == Catch::Approx(-power_1).margin(1e-14));
}
