#pragma once
#include "types.h"
#include "gas.h"
#include "particles.h" // Include the new header

struct SimState {
    ParticleSystem dm; // The dark matter class handles its own memory!
    GasGrid gas;

    // Time and cosmology track inside the state now
    double total_time = 0.0;
    double scale_factor = 1.0;
    double hubble_param = 0.0;

    // Constructor initializes the classes
    SimState( const Config& config )
        : dm( config ),
        gas( config )
    {
    }
};