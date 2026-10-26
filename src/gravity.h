#pragma once

struct Config;
struct Diagnostics;
struct SimState;
class GasParticleSystem;
class ParticleSystem;

class GravitySolver {
   public:
    GravitySolver();
    ~GravitySolver();

    void update_softening(SimState& state, Config& config);

    void add_pp_accel(ParticleSystem& dm, GasParticleSystem* gas, double a,
                      const Config& config, Diagnostics& diag);

    // Timestep
    double get_timestep(const ParticleSystem& dm, const Config& config) const;
    double get_timestep(const GasParticleSystem& gas,
                        const Config& config) const;
};
