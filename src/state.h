#pragma once
#include <memory>

#include "config.h"
#include "cooling.h"
#include "gas.h"
#include "gravity.h"
#include "mfm.h"
#include "particles.h"
#include "pm.h"
#include "types.h"

struct SimState {
    double total_time = 0.0;
    double scale_factor = 1.0;
    double hubble_param = 0.0;

    ParticleSystem dm;
    std::unique_ptr<GasGrid> gas;
    std::unique_ptr<GasParticleSystem> mfm_gas;
    Cooling cooling;
    PM pm;
    GravitySolver gravity;

    SimState(const Config& config);
};
