#pragma once

#include <Eigen/Dense>

#include "reconstruction.h"

namespace Riemann {
enum class SolverError {
    None,
    BadInputs,
    WaveInversion,
    NegativePressure,
    DerivativeCollapse,
    MaxIterations
};

// Output of the Riemann Solver
struct FaceFlux {
    Eigen::Vector3d flux_mom;  // Momentum flux density vector (P_star * n)
                               // [PHYSICAL Pressure Units]
    double P_star;  // Resolved face pressure [PHYSICAL Pressure Units]
    double S_star;  // Resolved relative contact wave speed [PHYSICAL Velocity
                    // Units]
    bool success = true;
    SolverError error = SolverError::None;
};

// Solves the HLLC Riemann problem at the moving face between particles.
// INPUT: Reconstructed face states, relative frame velocity (v_frame), and
// adiabatic index.
// OUTPUT: The resolved pressure (P_star) and relative wave speed (S_star).
// Note that it should get and return physical units
FaceFlux solve_hllc(const Reconstruction::ReconstructedFace& face,
                    const Eigen::Vector3d& v_frame, double gamma);

FaceFlux solve_exact(const Reconstruction::ReconstructedFace& face,
                     const Eigen::Vector3d& v_frame, double gamma);
}  // namespace Riemann
