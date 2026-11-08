#pragma once

#include <stddef.h>
#ifdef GPU_CPU_PARALLEL_EXECUTION
#include "worker.h"
#endif

struct Config;
struct Diagnostics;
struct SimState;
class GasParticleSystem;

class Integrator {
   public:
    void initialize(SimState& state, const Config& config);
    void KDK_step(SimState& state, double dt, const Config& config,
                  Diagnostics& diag);

    struct KickWork {
        double grav_work = 0.0;
        double exp_work = 0.0;
        double switch_energy = 0.0;
        double injected_energy = 0.0;
    };

    KickWork kick_single_gas_particle(GasParticleSystem& gas, size_t i,
                                      double dt, double a, double H,
                                      const Config& config);

    Integrator() = default;
    ~Integrator() = default;

   private:
    void apply_kick(SimState& state, double dt, const Config& config,
                    bool is_kick1);
    void apply_gas_kick(GasParticleSystem& gas, double dt, double a, double H,
                        const Config& config, bool is_kick1);
    void compute_forces(SimState& state, const Config& config,
                        Diagnostics& diag);

#ifdef GPU_CPU_PARALLEL_EXECUTION
    Worker pp_worker;
#endif
};
