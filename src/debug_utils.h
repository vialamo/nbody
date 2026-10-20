#pragma once

#include <vector>

#include "riemann.h"
#include "types.h"

class GasParticleSystem;

namespace DebugUtils {
void print_ill_conditioned_diagnostic(bool& debug_printed, size_t i,
                                      const GasParticleSystem& gas,
                                      const std::vector<int>& neighbors,
                                      const Vec3& domain_size);

void print_gradients(const GasParticleSystem& gas, size_t i);

void print_riemann_diagnostic(size_t i, int j,
                              const Reconstruction::ReconstructedFace& face,
                              const Eigen::Vector3d& v_frame, double gamma,
                              Riemann::SolverError error);

void print_non_converged_diagnostic(
    const GasParticleSystem& gas, size_t i, double h_guess, double h_low,
    double h_high, double effective_neighbors, double current_target_N,
    double target_N, double tol, int num_neighb, int min_num_neighb,
    int nr_iter, int max_iter, int density_retries, bool h_clamped,
    bool valid_density_found, bool& diagnostic_printed);

}  // namespace DebugUtils
