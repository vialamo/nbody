#include "debug_utils.h"

#include <iomanip>
#include <iostream>

#include "math_utils.h"
#include "mfm.h"

namespace DebugUtils {

void print_gradients(const GasParticleSystem& gas, size_t i) {
    // Trigger condition: highly clustered particles where the matrix may be
    // ill-conditioned
    if (gas.n_enc_final[i] > 150.0) {
        static int print_count = 0;

        // Limit to the first 10 occurrences to prevent OpenMP terminal flooding
        if (print_count < 10) {
#pragma omp critical
            {
                if (print_count < 10) {
                    std::cout
                        << "\n=== MFM GRADIENT DIVERGENCE DIAGNOSTIC ===\n"
                        << "Particle Index : " << i << "\n"
                        << "Position       : (" << gas.pos_x[i] << ", "
                        << gas.pos_y[i] << ", " << gas.pos_z[i] << ")\n"
                        << "N_enc          : " << gas.n_enc_final[i] << "\n"
                        << "Condition Num  : " << gas.cond_num[i] << "\n"
                        << std::scientific << std::setprecision(6)
                        << "grad_vx        : (" << gas.grad_vx[i].x() << ", "
                        << gas.grad_vx[i].y() << ", " << gas.grad_vx[i].z()
                        << ")\n"
                        << "grad_vy        : (" << gas.grad_vy[i].x() << ", "
                        << gas.grad_vy[i].y() << ", " << gas.grad_vy[i].z()
                        << ")\n"
                        << "grad_vz        : (" << gas.grad_vz[i].x() << ", "
                        << gas.grad_vz[i].y() << ", " << gas.grad_vz[i].z()
                        << ")\n"
                        << "grad_p         : (" << gas.grad_p[i].x() << ", "
                        << gas.grad_p[i].y() << ", " << gas.grad_p[i].z()
                        << ")\n"
                        << "==========================================\n\n";
                    print_count++;
                }
            }
        }
    }
}

void print_ill_conditioned_diagnostic(bool& debug_printed, size_t i,
                                      const GasParticleSystem& gas,
                                      const std::vector<int>& neighbors,
                                      const Vec3& domain_size) {
    if (!debug_printed) {
        debug_printed = true;
        std::cout << "\n=== ILL-CONDITIONED MATRIX DIAGNOSTIC ===\n";
        std::cout << "Particle Index: " << i << "\n";
        std::cout << "Pos: (" << gas.pos_x[i] << ", " << gas.pos_y[i] << ", "
                  << gas.pos_z[i] << ")\n";
        std::cout << "h_i: " << gas.h[i] << " | rho_i: " << gas.rho[i] << "\n";
        std::cout << "n_enc: " << gas.n_enc_final[i] << "\n";
        std::cout << "Condition Number: " << gas.cond_num[i] << "\n";
        std::cout << "Condition Number (conv): " << gas.cond_num_convergence[i]
                  << "\n";
        std::cout << "Number of neighbors found: " << neighbors.size() << "\n";

        int exact_overlaps = 0;
        int extremely_close = 0;

        std::cout << "Neighbor distances:\n";
        for (size_t n = 0; n < neighbors.size(); ++n) {
            int j = neighbors[n];
            double dx = periodic_displacement(gas.pos_x[j] - gas.pos_x[i],
                                              domain_size.x);
            double dy = periodic_displacement(gas.pos_y[j] - gas.pos_y[i],
                                              domain_size.y);
            double dz = periodic_displacement(gas.pos_z[j] - gas.pos_z[i],
                                              domain_size.z);
            double r = std::sqrt(dx * dx + dy * dy + dz * dz);

            std::cout << "  N[" << n << "] (Idx: " << j << ") r = " << r
                      << " | dx: " << dx << " dy: " << dy << " dz: " << dz
                      << "\n";

            if (r < 1e-12)
                exact_overlaps++;
            else if (r < 0.05 * gas.h[i])
                extremely_close++;
        }
        std::cout << "Diagnostics:\n";
        std::cout << "  Exact overlaps (r < 1e-12): " << exact_overlaps << "\n";
        std::cout << "  Highly clumped (r < 0.05*h): " << extremely_close
                  << "\n";
        std::cout << "=========================================\n";
    }
}

void print_riemann_diagnostic(size_t i, int j,
                              const Reconstruction::ReconstructedFace& face,
                              const Eigen::Vector3d& v_frame, double gamma,
                              Riemann::SolverError error) {
    static int fail_count = 0;

    // Restrict output to avoid OpenMP terminal flooding
    if (fail_count < 20) {
#pragma omp critical
        {
            if (fail_count < 20) {
                std::cout << "\n=== RIEMANN SOLVER FAILURE DIAGNOSTIC ===\n"
                          << "Interacting Pair: Particle " << i
                          << " and Particle " << j << "\n"
                          << "Failure Reason  : ";

                switch (error) {
                    case Riemann::SolverError::BadInputs:
                        std::cout << "Bad Inputs (NaN or Negative)\n";
                        break;
                    case Riemann::SolverError::WaveInversion:
                        std::cout << "Wave Speed Inversion (S_L >= S_R)\n";
                        break;
                    case Riemann::SolverError::NegativePressure:
                        std::cout
                            << "Negative Star Pressure (Strong Rarefaction)\n";
                        break;
                    case Riemann::SolverError::DerivativeCollapse:
                        std::cout << "Newton-Raphson Derivative Collapse\n";
                        break;
                    case Riemann::SolverError::MaxIterations:
                        std::cout << "Newton-Raphson Max Iterations Reached\n";
                        break;
                    default:
                        std::cout << "Unknown\n";
                        break;
                }

                double u_L = (face.v_L - v_frame).dot(face.n);
                double u_R = (face.v_R - v_frame).dot(face.n);
                double c_L = (face.rho_L > 0)
                                 ? std::sqrt(gamma * face.p_L / face.rho_L)
                                 : 0;
                double c_R = (face.rho_R > 0)
                                 ? std::sqrt(gamma * face.p_R / face.rho_R)
                                 : 0;

                std::cout << "--- Physical Face States (Face Frame) ---\n"
                          << std::scientific << std::setprecision(6)
                          << "rho_L : " << face.rho_L
                          << " | rho_R : " << face.rho_R << "\n"
                          << "p_L   : " << face.p_L << " | p_R   : " << face.p_R
                          << "\n"
                          << "u_L   : " << u_L << " | u_R   : " << u_R << "\n"
                          << "c_L   : " << c_L << " | c_R   : " << c_R << "\n"
                          << "Face Normal : (" << face.n.x() << ", "
                          << face.n.y() << ", " << face.n.z() << ")\n"
                          << "=========================================\n\n";
                fail_count++;
            }
        }
    }
}

void print_non_converged_diagnostic(
    const GasParticleSystem& gas, size_t i, double h_guess, double h_low,
    double h_high, double effective_neighbors, double current_target_N,
    double target_N, double tol, int num_neighb, int min_num_neighb,
    int nr_iter, int max_iter, int density_retries, bool h_clamped,
    bool valid_density_found, bool& diagnostic_printed) {
    if (!diagnostic_printed) {
#pragma omp critical
        {
            if (!diagnostic_printed) {
                std::cout
                    << "\n=== PASS 1: NON-CONVERGED PARTICLE DIAGNOSTIC ===\n"
                    << "Global Particle Index : " << i << "\n"
                    << "Position            : (" << gas.pos_x[i] << ", "
                    << gas.pos_y[i] << ", " << gas.pos_z[i] << ")\n"
                    << "Final h_guess       : " << std::setprecision(12)
                    << h_guess << "\n"
                    << "h Bracket [low, high]: [" << h_low << ", " << h_high
                    << "]\n"
                    << "Bracket Width       : " << (h_high - h_low) << "\n"
                    << "Effective Neighbors : " << effective_neighbors << "\n"
                    << "Current Target N    : " << current_target_N
                    << " (Base N: " << target_N << ")\n"
                    << "Difference to Target: "
                    << std::abs(effective_neighbors - current_target_N)
                    << " (Tol: " << tol << ")\n"
                    << "Discrete Neighbors  : " << num_neighb
                    << " (Min required: " << min_num_neighb << ")\n"
                    << "Newton-Raphson Iters: " << nr_iter << " / " << max_iter
                    << "\n"
                    << "Density Retries     : " << density_retries << "\n"
                    << "h_clamped           : "
                    << (h_clamped ? "TRUE" : "FALSE") << "\n"
                    << "valid_density_found : "
                    << (valid_density_found ? "TRUE" : "FALSE") << "\n"
                    << "=================================================\n\n";
                diagnostic_printed = true;
            }
        }
    }
}

}  // namespace DebugUtils
