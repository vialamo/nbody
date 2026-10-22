#pragma once
#include <memory>

#include "config.h"
#include "cooling.h"
#include "gas.h"
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

    SimState(const Config& config)
        : dm(config),
          gas(nullptr),
          mfm_gas(nullptr),
          cooling(config),
          pm(config) {
        if (config.hydro_method == HydroMethod::Eulerian) {
            gas = std::make_unique<GasGrid>(config);
        } else if (config.hydro_method == HydroMethod::MFM) {
            mfm_gas = std::make_unique<GasParticleSystem>(config);
        }
    }
};
