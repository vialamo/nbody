#pragma once

struct Config;
struct SimState;

namespace Cosmology {

void update(SimState& state, const Config& config);
double get_time_from_scale_factor(double a, const Config& config);

}
