#pragma once
#include <memory>

#include "config.h"
#include "cooling.h"
#include "gas.h"
#include "gravity.h"
#include "particles.h"
#include "pm.h"
#include "star.h"
#include "types.h"

struct SimState {
    uint64_t ti_total_time = 0;
    double total_time = 0.0;
    double scale_factor = 1.0;
    double hubble_param = 0.0;

    ParticleSystem dm;
    GasParticleSystem gas;
    StarParticleSystem star;
    Cooling cooling;
    PM pm;
    GravitySolver gravity;

    SimState(const Config& config);
};
