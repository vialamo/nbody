#pragma once

struct Config;
class GasParticleSystem;

namespace HydroSolver {
constexpr double MIN_SIGNAL_VELOCITY = 1e-9;

void step(GasParticleSystem& gas, const Config& config, double a);

double get_timestep(const GasParticleSystem& gas, const Config& config);

}  // namespace HydroSolver
