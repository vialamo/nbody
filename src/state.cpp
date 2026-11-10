#include "state.h"

SimState::SimState(const Config& config)
    : dm(config), gas(config), cooling(config), pm(config), gravity() {}
