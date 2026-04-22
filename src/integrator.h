#pragma once
#include <map>
#include <string>

#include "config.h"
#include "state.h"

void update_cosmology(SimState& state, const Config& config);

// FFT Gravity Solver
void compute_gravitational_acceleration(SimState& state, const Config& config);

std::map<std::string, double> KDK_step(SimState& state, double dt,
                                       Config& config);
