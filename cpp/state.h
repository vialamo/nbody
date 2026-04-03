#pragma once
#include "types.h"
#include "gas.h"
#include "particles.h"

struct SimState {
    ParticleSystem dm;
    GasGrid gas;

    double total_time = 0.0;
    double scale_factor = 1.0;
    double hubble_param = 0.0;

    SimState( const Config& config )
        : dm( config ),
        gas( config )
    {
    }
};