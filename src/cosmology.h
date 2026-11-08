#pragma once

struct Config;
struct SimState;

namespace Cosmology {

void update(SimState& state, const Config& config);
double get_time_from_scale_factor(double a, const Config& config);

void compute_step_factors(double a, double H, double dt, double gamma,
                          double& drag_factor, double& exp_factor);

}  // namespace Cosmology
