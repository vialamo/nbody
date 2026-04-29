#pragma once
#include <map>
#include <string>

#include "config.h"
#include "diagnostics.h"
#include "state.h"

void update_cosmology(SimState& state, const Config& config);

// FFT Gravity Solver
void compute_gravitational_acceleration(SimState& state, const Config& config);

void KDK_step(SimState& state, double dt, Config& config, Diagnostics& diag);
