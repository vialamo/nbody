#pragma once

struct Config;
class GasParticleSystem;

namespace HydroSolver {

void step(GasParticleSystem& gas, const Config& config, double a);

double get_timestep(GasParticleSystem& gas, const Config& config);

}  // namespace HydroSolver
