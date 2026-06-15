#pragma once
#include <map>
#include <string>

#include "config.h"
#include "diagnostics.h"
#include "state.h"
#include "types.h"

void update_cosmology(SimState& state, const Config& config);

double get_time_from_scale_factor(double a, const Config& config);

// FFT Gravity Solver
void compute_gravitational_acceleration(SimState& state, const Config& config);

void compute_forces(SimState& state, Config& config, Diagnostics& diag);

void KDK_step(SimState& state, TimestepInfo& ts, Config& config, Diagnostics& diag);
