#pragma once
#include <map>
#include <string>

#include "config.h"
#include "state.h"

void update_cosmology(SimState& state, const Config& config);

// FFT Gravity Solver
Grid3D compute_gravitational_acceleration(Grid3D& acc_x, Grid3D& acc_y,
                                          Grid3D& acc_z, GasGrid& gas,
                                          const Config& config,
                                          const Grid3D& dm_rho);

std::map<std::string, double> KDK_step(SimState& state, double dt,
                                       Config& config);
