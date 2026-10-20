#include "riemann.h"

namespace Riemann {

// Solves the HLLC Riemann problem at the moving face between particles.
// INPUT: Reconstructed face states, relative frame velocity (v_frame), and
// adiabatic index.
// OUTPUT: The resolved pressure (P_star) and relative wave speed (S_star).
// Note that it should get and return physical units
FaceFlux solve_hllc(const Reconstruction::ReconstructedFace& face,
                    const Eigen::Vector3d& v_frame, double gamma) {
    FaceFlux out;
    out.flux_mom = Eigen::Vector3d::Zero();
    out.P_star = 0.0;
    out.S_star = 0.0;
    out.success = true;
    out.error = SolverError::None;

    // Check for unphysical or NaN input states
    if (std::isnan(face.p_L) || std::isnan(face.p_R) || face.p_L < 0.0 ||
        face.p_R < 0.0 || std::isnan(face.rho_L) || std::isnan(face.rho_R) ||
        face.rho_L <= 0.0 || face.rho_R <= 0.0) {
        out.success = false;
        out.error = SolverError::BadInputs;
        return out;
    }

    // Boost and Project Left & Right States to the face frame
    double vn_L = (face.v_L - v_frame).dot(face.n);
    double vn_R = (face.v_R - v_frame).dot(face.n);

    // Compute Sound Speeds
    double cs_L =
        (face.rho_L > 1e-12) ? std::sqrt(gamma * face.p_L / face.rho_L) : 0.0;
    double cs_R =
        (face.rho_R > 1e-12) ? std::sqrt(gamma * face.p_R / face.rho_R) : 0.0;

    // Wave Speed Estimates
    double S_L = std::min(vn_L - cs_L, vn_R - cs_R);
    double S_R = std::max(vn_L + cs_L, vn_R + cs_R);

    // Check for wave speed inversion
    if (S_L >= S_R) {
        out.success = false;
        out.error = SolverError::WaveInversion;
        return out;
    }

    // Supersonic to the right: face sees only the Left state
    if (S_L >= 0.0) {
        out.P_star = face.p_L;
        out.S_star = vn_L;
        out.flux_mom = out.P_star * face.n;
        return out;
    }

    // Supersonic to the left: face sees only the Right state
    if (S_R <= 0.0) {
        out.P_star = face.p_R;
        out.S_star = vn_R;
        out.flux_mom = out.P_star * face.n;
        return out;
    }

    // Subsonic: contact wave lies inside the fan
    // Solve the HLLC Riemann Problem for the Star State
    // MFM always samples the Star State, because the face moves at S_star
    double den_star = face.rho_L * (S_L - vn_L) - face.rho_R * (S_R - vn_R);
    if (std::abs(den_star) < 1e-15) den_star = (den_star > 0 ? 1e-15 : -1e-15);

    double S_star = (face.p_R - face.p_L + face.rho_L * vn_L * (S_L - vn_L) -
                     face.rho_R * vn_R * (S_R - vn_R)) /
                    den_star;

    // Enforce physical bounds on the Riemann fan
    S_star = std::max(S_L, std::min(S_star, S_R));

    double P_star = face.p_L + face.rho_L * (S_L - vn_L) * (S_star - vn_L);
    // Check for strong rarefaction causing negative pressure
    if (P_star < 0.0) {
        out.success = false;
        out.error = SolverError::NegativePressure;
        return out;
    }

    // Assemble MFM Fluxes (Mass flux is zero)
    out.flux_mom = P_star * face.n;
    out.P_star = P_star;
    out.S_star = S_star;

    return out;
}

// MFM Exact Riemann Solver (Toro 2009)
FaceFlux solve_exact(const Reconstruction::ReconstructedFace& face,
                     const Eigen::Vector3d& v_frame, double gamma) {
    FaceFlux out;
    out.flux_mom = Eigen::Vector3d::Zero();
    out.P_star = 0.0;
    out.S_star = 0.0;
    out.success = true;

    // Check for unphysical or NaN input states
    if (std::isnan(face.p_L) || std::isnan(face.p_R) || face.p_L < 0.0 ||
        face.p_R < 0.0 || std::isnan(face.rho_L) || std::isnan(face.rho_R) ||
        face.rho_L <= 0.0 || face.rho_R <= 0.0) {
        out.success = false;
        out.error = SolverError::BadInputs;
        return out;
    }

    // Boost and Project Left & Right States to the face frame
    double u_L = (face.v_L - v_frame).dot(face.n);
    double u_R = (face.v_R - v_frame).dot(face.n);
    double P_L = face.p_L;
    double P_R = face.p_R;
    double rho_L = face.rho_L;
    double rho_R = face.rho_R;

    // Compute Sound Speeds
    double c_L = (rho_L > 1e-12) ? std::sqrt(gamma * P_L / rho_L) : 0.0;
    double c_R = (rho_R > 1e-12) ? std::sqrt(gamma * P_R / rho_R) : 0.0;

    // Vacuum Generation Check
    if (2.0 / (gamma - 1.0) * (c_L + c_R) <= (u_R - u_L)) {
        out.P_star = 0.0;
        out.S_star = 0.5 * (u_L + u_R);
        out.flux_mom = Eigen::Vector3d::Zero();
        return out;
    }

    // Newton-Raphson Iteration for Exact P_star
    double P_star = 0.5 * (P_L + P_R);
    const int MAX_ITER = 100;
    const double TOL = 1e-6;
    double mu2 = (gamma - 1.0) / (gamma + 1.0);

    for (int iter = 0; iter < MAX_ITER; ++iter) {
        double f_L, f_R, df_L, df_R;

        // Left Wave (Shock or Rarefaction)
        if (P_star > P_L) {
            double A_L = 2.0 / ((gamma + 1.0) * rho_L);
            double B_L = mu2 * P_L;
            f_L = (P_star - P_L) * std::sqrt(A_L / (P_star + B_L));
            df_L = std::sqrt(A_L / (P_star + B_L)) *
                   (1.0 - 0.5 * (P_star - P_L) / (P_star + B_L));
        } else {
            double pratio = P_star / P_L;
            double pterm = std::pow(pratio, (gamma - 1.0) / (2.0 * gamma));
            f_L = (2.0 * c_L / (gamma - 1.0)) * (pterm - 1.0);
            df_L = (1.0 / (rho_L * c_L)) *
                   std::pow(pratio, -(gamma + 1.0) / (2.0 * gamma));
        }

        // Right Wave (Shock or Rarefaction)
        if (P_star > P_R) {
            double A_R = 2.0 / ((gamma + 1.0) * rho_R);
            double B_R = mu2 * P_R;
            f_R = (P_star - P_R) * std::sqrt(A_R / (P_star + B_R));
            df_R = std::sqrt(A_R / (P_star + B_R)) *
                   (1.0 - 0.5 * (P_star - P_R) / (P_star + B_R));
        } else {
            double pratio = P_star / P_R;
            double pterm = std::pow(pratio, (gamma - 1.0) / (2.0 * gamma));
            f_R = (2.0 * c_R / (gamma - 1.0)) * (pterm - 1.0);
            df_R = (1.0 / (rho_R * c_R)) *
                   std::pow(pratio, -(gamma + 1.0) / (2.0 * gamma));
        }

        // Root equation f(P) = f_L + f_R + u_R - u_L = 0
        double f = f_L + f_R + u_R - u_L;
        double df = df_L + df_R;

        // Derivative collapse check
        if (std::abs(df) < 1e-15) {
            out.success = false;
            out.error = SolverError::DerivativeCollapse;
            return out;
        }

        double dP = -f / df;
        double P_new = P_star + dP;

        // Instead of snapping to a hard floor, halve the current pressure
        // to smoothly asymptote towards the near-vacuum true root
        if (P_new <= 0.0) {
            P_new = 0.5 * P_star;
        }

        // Use a small epsilon (1e-14) for the relative error denominator
        double rel_err = std::abs(P_new - P_star) / (P_star + 1e-14);
        P_star = P_new;

        if (rel_err < TOL) break;

        // Non-convergence check
        if (iter == MAX_ITER - 1) {
            out.success = false;
            out.error = SolverError::MaxIterations;
            return out;
        }
    }

    // Calculate exact S_star using converged P_star
    double f_L_final;
    if (P_star > P_L) {
        double A_L = 2.0 / ((gamma + 1.0) * rho_L);
        double B_L = mu2 * P_L;
        f_L_final = (P_star - P_L) * std::sqrt(A_L / (P_star + B_L));
    } else {
        double pterm = std::pow(P_star / P_L, (gamma - 1.0) / (2.0 * gamma));
        f_L_final = (2.0 * c_L / (gamma - 1.0)) * (pterm - 1.0);
    }

    // S_star acts as the contact discontinuity velocity
    double S_star = u_L - f_L_final;

    out.P_star = P_star;
    out.S_star = S_star;
    out.flux_mom = P_star * face.n;

    return out;
}

}  // namespace Riemann
