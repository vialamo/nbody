#include "state.h"

SimState::SimState(const Config& config)
    : dm(config),
      gas(nullptr),
      mfm_gas(nullptr),
      cooling(config),
      pm(config),
      gravity() {
    if (config.hydro_method == HydroMethod::Eulerian) {
        gas = std::make_unique<GasGrid>(config);
    } else if (config.hydro_method == HydroMethod::MFM) {
        mfm_gas = std::make_unique<GasParticleSystem>(config);
    }
}
