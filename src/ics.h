#pragma once
#include "config.h"
#include "state.h"

SimState initialize_state(Config& config);

struct ZeldovichField {
    std::vector<double> dx;
    std::vector<double> dy;
    std::vector<double> dz;
    double f;
};

ZeldovichField compute_zeldovich_field(double scale_factor,
                                       const Config& config);
