#pragma once
#include "config.h"
#include "state.h"

SimState initialize_state(Config& config);

double get_internal_energy_from_temp_k(double T_kelvin, double V_unit_km_s,
                                       double gamma);

struct ZeldovichField {
    std::vector<double> dx;
    std::vector<double> dy;
    std::vector<double> dz;
    double f;
};

ZeldovichField compute_zeldovich_field(double scale_factor,
                                       const Config& config);
