#pragma once
#include "state.h"
#include "config.h"
#include <map>
#include <string>

void update_cosmology( SimState& state, const Config& config );

// FFT Gravity Solver
Grid3D compute_gravitational_acceleration(
    GasGrid& gas,
    const Config& config,
    const Grid3D& dm_rho );

// The Master Step
std::map<std::string, double> KDK_step(
    SimState& state,
    double dt,
    Config& config );
