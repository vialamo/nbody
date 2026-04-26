#pragma once
#include "config.h"
#include "state.h"

SimState initialize_state(Config& config);

double get_internal_energy_from_temp_k(double T_kelvin, double hubble_param,
                                       double box_size_mpc, double domain_size,
                                       double gamma);

struct ZeldovichField {
    std::vector<double> dx;
    std::vector<double> dy;
    std::vector<double> dz;
    double f;
};

ZeldovichField compute_zeldovich_field(double scale_factor,
                                       const Config& config);
