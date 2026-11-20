#pragma once

#include <stddef.h>

#include <cstdint>
#ifdef GPU_CPU_PARALLEL_EXECUTION
#include "worker.h"
#endif

struct Config;
struct Diagnostics;
struct SimState;
class GasParticleSystem;

class Integrator {
   public:
    struct KickWork {
        double grav_work = 0.0;
        double exp_work = 0.0;
        double def_energy = 0.0;
        double entropy_energy = 0.0;
        double injected_energy = 0.0;
    };

    void initialize(SimState& state, const Config& config);
    void KDK_step(SimState& state, uint64_t ti_dt, const Config& config,
                  Diagnostics& diag);

    void kick_wakeup_catchup(GasParticleSystem& gas, size_t i, double dt_actual,
                             double a, double H, const Config& config);

    Integrator() = default;
    ~Integrator() = default;

   private:
    void apply_kick(SimState& state, uint64_t ti_dt, const Config& config,
                    bool is_kick1);
    void apply_gas_kick(GasParticleSystem& gas, uint64_t ti_dt, double a,
                        double H, const Config& config, bool is_kick1);
    void compute_forces(SimState& state, const Config& config,
                        Diagnostics& diag);

#ifdef GPU_CPU_PARALLEL_EXECUTION
    Worker pp_worker;
#endif
};
